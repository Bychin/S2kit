/*
    FST_semi_memo.c - routines to perform convolutions on the 2-sphere using a combination of
    semi-naive and naive algorithms.

    Assumes that all precomputed data is already in memory, not to be read from disk.

    The primary functions in this package are:
    // TODO rename
    SHTSemiMemo()
    InvSHTSemiMemo()
    ZHTSemiMemo()
    ConvOn2SphereSemiMemo()
    1) FST_semi_memo()         - computes the spherical harmonic transform;
    2) InvFST_semi_memo()      - computes the inverse spherical harmonic transform;
    3) FZT_semi_memo()         - computes the zonal harmonic transform;
    4) Conv2Sphere_semi_memo() - convolves two functins defined on the 2-sphere, using seminaive transform.

    For descriptions on calling these functions, see the documentation preceding each function.
*/

#include "FST_semi_memo.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "legendre_polynomials/cospml.h"
#include "legendre_transform/naive.h"
#include "legendre_transform/seminaive.h"
#include "legendre_transform/weights.h"
#include "util/chebyshev_nodes.h"
#include "util/util.h"

/*
    Performs a spherical harmonic transform using the semi-naive
    and naive algorithms

    bw - bandwidth of problem
    size -> size = 2*bw -> dimension of input array (recall that
           sampling is done at twice the bandwidth)

   The inputs rdata and idata are expected to be pointers to
   size x size arrays. The array rdata contains the real parts
   of the function samples, and idata contains the imaginary
   parts.

   rcoeffs and icoeffs are expected to be pointers to bw x bw arrays,
   and will contain the harmonic coefficients in a "linearized" form.
   The array rcoeffs contains the real parts of the coefficients,
   and icoeffs contains the imaginary parts.

   spharmonic_pml_table should be a (double **) pointer to
   the result of a call to Spharmonic_Pml_Table.  Because this
   table is re-used in the inverse transform, and because for
   timing purposes the computation of the table is not included,
   it is passed in as an argument.  Also, at some point this
   code may be used as par of a series of convolutions, so
   reducing repetitive computation is prioritized.

   spharmonic_pml_table will be an array of (double *) pointers
   the array being of length TableSize(m,bw)

   workspace needs to be a double pointer to an array of size
   (8 * bw^2) + (7 * bw).

   cutoff -> what order to switch from semi-naive to naive
             algorithm.

   dataformat =0 -> samples are complex, =1 -> samples real


   Output Ordering of coeffs f(m,l) is
   f(0,0) f(0,1) f(0,2) ... f(0,bw-1)
          f(1,1) f(1,2) ... f(1,bw-1)
          etc.
                 f(bw-2,bw-2), f(bw-2,bw-1)
                       f(bw-1,bw-1)
                   f(-(bw-1),bw-1)
         f(-(bw-2),bw-2) f(-(bw-2),bw-1)
      etc.
              f(-2,2) ... f(-2,bw-1)
      f(-1,1) f(-1,2) ... f(-1,bw-1)


   This only requires an array of size (bw*bw).  If zero-padding
   is used to make the indexing nice, then you need a an
   (2bw-1) * bw array - but that is not done here.
   Because of the amount of space necessary for doing
   large transforms, it is important not to use any
   more than necessary.

*/
// TODO: dataformat to enum?
void FST_semi_memo(double* rdata, double* idata, double* rcoeffs, double* icoeffs, const int bw,
                   double** seminaive_naive_table, double* workspace, const int dataformat, const int cutoff,
                   fftw_plan* DCT_plan, fftw_plan* FFT_plan, double* weights) {
    int size = 2 * bw;

    // total workspace is (8 * bw^2) + (7 * bw)
    double* rres = workspace;               // needs (size * size) = (4 * bw^2)
    double* ires = rres + (size * size);    // needs (size * size) = (4 * bw^2)
    double* fltres = ires + (size * size);  // needs bw
    double* eval_pts = fltres + bw;         // needs (2 * bw)
    double* scratchpad = eval_pts + (size); // needs (4 * bw)

    // do the FFTs along phi
    fftw_execute_split_dft(*FFT_plan, rdata, idata, rres, ires);

    /*
        normalize

        // TODO rewrite
        note that I'm getting the sqrt(2*pi) in there at
        this point ... to account for the fact that the spherical
        harmonics are of norm 1: I need to account for
        the fact that the associated Legendres are
        of norm 1
    */
    double normed_coeff = sqrt(2. * M_PI) / size;
    for (int i = 0; i < size * size; ++i) {
        rres[i] *= normed_coeff;
        ires[i] *= normed_coeff;
    }

    // point to start of output data buffers
    double* rdataptr = rcoeffs;
    double* idataptr = icoeffs;

    for (int m = 0; m < cutoff; ++m) { // semi-naive part
        // real part
        SemiNaiveReduced(rres + (m * size), bw, m, fltres, scratchpad, seminaive_naive_table[m], weights, DCT_plan);
        // load real part of coefficients into output space
        memcpy(rdataptr, fltres, sizeof(double) * (bw - m));
        rdataptr += bw - m;

        // imaginary part
        SemiNaiveReduced(ires + (m * size), bw, m, fltres, scratchpad, seminaive_naive_table[m], weights, DCT_plan);
        // load imaginary part of coefficients into output space
        memcpy(idataptr, fltres, sizeof(double) * (bw - m));
        idataptr += bw - m;
    }

    for (int m = cutoff; m < bw; ++m) { // naive part
        // real part
        Naive_AnalysisX(rres + (m * size), bw, m, weights, fltres, seminaive_naive_table[m], scratchpad);
        memcpy(rdataptr, fltres, sizeof(double) * (bw - m));
        rdataptr += bw - m;

        // imaginary part
        Naive_AnalysisX(ires + (m * size), bw, m, weights, fltres, seminaive_naive_table[m], scratchpad);
        memcpy(idataptr, fltres, sizeof(double) * (bw - m));
        idataptr += bw - m;
    }

    // upper coefficients
    /*
        If the data is real, we don't have to compute the coeffs whose
        order is less than 0, since the data is real, we know that

        f-hat(l,-m) = (-1)^m * conjugate(f-hat(l,m)),

        so use that to get the rest of the coefficients
    */
    if (dataformat == 1) { // real data
        double coeff = 1.;
        for (int i = 1; i < bw; ++i) {
            coeff *= -1.;
            for (int j = i; j < bw; ++j) {
                int index0 = seanindex(i, j, bw);
                int index1 = seanindex(-i, j, bw);

                rcoeffs[index1] = coeff * rcoeffs[index0];
                icoeffs[index1] = -coeff * icoeffs[index0];
            }
        }

        return;
    }
    
    // complex data
    /* 
        Note that `m` is greater than `bw` here, but this is for
        purposes of indexing the input data arrays. The "true"
        value of `m` as a parameter for Pml is `size-m` 
    */
    // TODO into two cycles?
    for (int m = bw + 1; m < size; ++m) {
        if ((size - m) < cutoff) { // semi-naive part
            // real part
            SemiNaiveReduced(rres + (m * size), bw, size - m, fltres, scratchpad, seminaive_naive_table[size - m],
                             weights, DCT_plan);
            // load real part of coefficients into output space
            if (m % 2) {
                for (int i = 0; i < m - bw; ++i)
                    rdataptr[i] = -fltres[i];
            } else {
                memcpy(rdataptr, fltres, sizeof(double) * (m - bw));
            }
            rdataptr += m - bw;

            // imaginary part
            SemiNaiveReduced(ires + (m * size), bw, size - m, fltres, scratchpad, seminaive_naive_table[size - m],
                             weights, DCT_plan);
            // load imaginary part of coefficients into output space
            if (m % 2) {
                for (int i = 0; i < m - bw; ++i)
                    idataptr[i] = -fltres[i];
            } else {
                memcpy(idataptr, fltres, sizeof(double) * (m - bw));
            }
            idataptr += m - bw;

            continue;
        }
        
        // naive part
        // real part
        Naive_AnalysisX(rres + (m * size), bw, size - m, weights, fltres, seminaive_naive_table[size - m],
                        scratchpad);
        // load real part of coefficients into output space
        if (m % 2) {
            for (int i = 0; i < m - bw; ++i)
                rdataptr[i] = -fltres[i];
        } else {
            memcpy(rdataptr, fltres, sizeof(double) * (m - bw));
        }
        rdataptr += m - bw;

        // imaginary part
        Naive_AnalysisX(ires + (m * size), bw, size - m, weights, fltres, seminaive_naive_table[size - m],
                        scratchpad);
        // load imaginary part of coefficients into output space
        if (m % 2) {
            for (int i = 0; i < m - bw; ++i)
                idataptr[i] = -fltres[i];
        } else {
            memcpy(idataptr, fltres, sizeof(double) * (m - bw));
        }
        idataptr += m - bw;
    }
}

/*
    Inverse spherical harmonic transform.

    bw - bandwidth of problem

    Inputs rcoeffs and icoeffs are harmonic coefficients stored
    in `bw*bw` arrays in the order spec'ed above.

   rdata and idata are (2*bw x 2*bw) arrays with the transformed result.

   transpose_spharmonic_pml_table should be the (double **)
   result of a call to Transpose_Spharmonic_Pml_Table()

   workspace is (8 * bw^2) + (10 * bw)

   dataformat =0 -> samples are complex, =1 -> samples real

*/
// TODO try to use double** -> double*
void InvFST_semi_memo(double* rcoeffs, double* icoeffs, double* rdata, double* idata, const int bw,
                      double** transpose_seminaive_naive_table, double* workspace, const int dataformat,
                      const int cutoff, fftw_plan* inv_DCT_plan, fftw_plan* inv_FFT_plan) {
    int size = 2 * bw;

    // total workspace = (8 * bw^2) + (10 * bw)
    double* rfourdata = workspace;                  // needs (size * size)
    double* ifourdata = rfourdata + (size * size);  // needs (size * size)
    double* rinvfltres = ifourdata + (size * size); // needs (2 * bw)
    double* iminvfltres = rinvfltres + (size);      // needs (2 * bw)
    double* sin_values = iminvfltres + (size);      // needs (2 * bw)
    double* eval_pts = sin_values + (size);         // needs (2 * bw)
    double* scratchpad = eval_pts + (size);         // needs (2 * bw)

    // load up the sin_values array
    AcosOfChebyshevNodes(size, eval_pts);
    for (int i = 0; i < size; ++i)
        sin_values[i] = sin(eval_pts[i]);

    // do all of the inverse Legendre transforms
    double* rdataptr = rcoeffs;
    double* idataptr = icoeffs;

    for (int m = 0; m < cutoff; ++m) { // semi-naive part
        // real part
        InvSemiNaiveReduced(rdataptr, bw, m, rinvfltres, transpose_seminaive_naive_table[m], sin_values,
                            scratchpad, inv_DCT_plan);

        // imaginary part
        InvSemiNaiveReduced(idataptr, bw, m, iminvfltres, transpose_seminaive_naive_table[m], sin_values,
                            scratchpad, inv_DCT_plan);

        // store normal, then tranpose before doing inverse fft
        memcpy(rfourdata + (m * size), rinvfltres, sizeof(double) * size);
        memcpy(ifourdata + (m * size), iminvfltres, sizeof(double) * size);

        rdataptr += bw - m;
        idataptr += bw - m;
    }

    for (int m = cutoff; m < bw; ++m) { // naive part
        // real part
        Naive_SynthesizeX(rdataptr, bw, m, rinvfltres, transpose_seminaive_naive_table[m]);
        // imaginary part
        Naive_SynthesizeX(idataptr, bw, m, iminvfltres, transpose_seminaive_naive_table[m]);

        // store normal, then tranpose before doing inverse fft
        memcpy(rfourdata + (m * size), rinvfltres, sizeof(double) * size);
        memcpy(ifourdata + (m * size), iminvfltres, sizeof(double) * size);

        rdataptr += bw - m;
        idataptr += bw - m;
    }

    // fill in zero values where m = bw (from problem definition)
    memset(rfourdata + (bw * size), 0, sizeof(double) * size);
    memset(ifourdata + (bw * size), 0, sizeof(double) * size);

    /*
        If the data is real, we don't have to compute the coeffs whose order is less than 0,
        i.e. since the data is real, we know that

        invf-hat(l,-m) = conjugate(invf-hat(l,m)),

        so use that to get the rest of the coefficients
    */
    if (dataformat == 0) { // complex data
        // do negative m values
        // TODO into two cycles?
        for (int m = bw + 1; m < size; ++m) {
            if ((size - m) < cutoff) { // semi-naive part
                InvSemiNaiveReduced(rdataptr, bw, size - m, rinvfltres, transpose_seminaive_naive_table[size - m],
                                    sin_values, scratchpad, inv_DCT_plan);
                InvSemiNaiveReduced(idataptr, bw, size - m, iminvfltres, transpose_seminaive_naive_table[size - m],
                                    sin_values, scratchpad, inv_DCT_plan);

                // store normal, then tranpose before doing inverse fft
                if (m % 2)
                    for (int i = 0; i < size; ++i) {
                        rinvfltres[i] *= -1.0;
                        iminvfltres[i] *= -1.0;
                    }

                memcpy(rfourdata + (m * size), rinvfltres, sizeof(double) * size);
                memcpy(ifourdata + (m * size), iminvfltres, sizeof(double) * size);

                rdataptr += bw - (size - m);
                idataptr += bw - (size - m);

                continue;
            }

            // naive part
            Naive_SynthesizeX(rdataptr, bw, size - m, rinvfltres, transpose_seminaive_naive_table[size - m]);
            Naive_SynthesizeX(idataptr, bw, size - m, iminvfltres, transpose_seminaive_naive_table[size - m]);

            // store normal, then tranpose before doing inverse fft
            if (m % 2)
                for (int i = 0; i < size; ++i) {
                    rinvfltres[i] *= -1.0;
                    iminvfltres[i] *= -1.0;
                }

            memcpy(rfourdata + (m * size), rinvfltres, sizeof(double) * size);
            memcpy(ifourdata + (m * size), iminvfltres, sizeof(double) * size);

            rdataptr += bw - (size - m);
            idataptr += bw - (size - m);
        }
    } else { // real data
        for (int m = bw + 1; m < size; ++m) {
            memcpy(rfourdata + (m * size), rfourdata + ((size - m) * size), sizeof(double) * size);
            memcpy(ifourdata + (m * size), ifourdata + ((size - m) * size), sizeof(double) * size);

            for (int i = 0; i < size; ++i)
                ifourdata[(m * size) + i] *= -1.0;
        }
    }

    // normalize
    double normed_coeff = 1. / (sqrt(2. * M_PI));
    for (int i = 0; i < 4 * bw * bw; i++) {
        rfourdata[i] *= normed_coeff;
        ifourdata[i] *= normed_coeff;
    }

    fftw_execute_split_dft(*inv_FFT_plan, ifourdata, rfourdata, idata, rdata);
}

/*
    Zonal Harmonic transform using seminaive algorithm - used in convolutions

  bw -> bandwidth of problem

  size = 2 * bw

  rdata and idata should be pointers to size x size arrays.
  rres and ires should be pointers to double arrays of size bw.

  cos_pml_table contains Legendre coefficients of P(0,l) functions
  and is result of CosPmlTableGen for m = 0;
  FZT_semi only computes spherical harmonics for m=0.

  dataformat =0 -> samples are complex, =1 -> samples real

  workspace needed is (12 * bw)

*/
void FZT_semi_memo(double* rdata, double* idata, double* rres, double* ires, const int bw, double* cos_pml_table,
                   double* workspace, const int dataformat, fftw_plan* DCT_plan, double* weights) {
    int size = 2 * bw;

    double* r0 = workspace;         // needs (2 * bw)
    double* i0 = r0 + size;         // needs (2 * bw)
    double* scratchpad = i0 + size; // needs (4 * bw)

    // total workspace = 13*bw

    double dsize = sqrt(2. * M_PI) / size;
    // compute the m = 0 components
    for (int i = 0; i < size; ++i) {
        double tmpreal = 0.0;
        double tmpimag = 0.0;

        for (int j = 0; j < size; ++j) {
            tmpreal += rdata[(i * size) + j];
            tmpimag += idata[(i * size) + j];
        }

        // normalize
        r0[i] = tmpreal * dsize;
        i0[i] = tmpimag * dsize;
    }

    // real part
    SemiNaiveReduced(r0, bw, 0, rres, scratchpad, cos_pml_table, weights, DCT_plan);

    if (dataformat == 0) // if data is complex, do imaginary part
        SemiNaiveReduced(i0, bw, 0, ires, scratchpad, cos_pml_table, weights, DCT_plan);
    else // otherwise set coefficients to zero
        memset(ires, 0, sizeof(double) * size);
}

/*
    Convolves two functions defined on the 2-sphere.
   Uses seminaive algorithms for spherical harmonic transforms

   size = 2*bw

   Inputs:

   rdata, idata - (size * size) arrays containing real and
                  imaginary parts of sampled function.
   rfilter, ifilter - (size * size) arrays containing real and
                      imaginary parts of sampled filter function.
   rres, ires - (size * size) arrays containing real and
                  imaginary parts of result function.


   Suggestion - if you want to do multiple convolutions,
   don't keep allocating and freeing space with every call,
   or keep recomputing the spharmonic_pml tables.
   Allocate workspace once before you call this function, then
   just set up pointers as first step of this procedure rather
   than mallocing.  And do the same with the FST, FZT, and InvFST functions.

   ASSUMPTIONS:
   1. data is strictly REAL
   2. will do semi-naive algorithm for ALL orders -> change the cutoff
      value if you want it to be different

   Memory requirements for Conv2Sphere

   Need space for spharmonic tables and local workspace and
   scratchpad space for FST_semi

   Let legendre_size = Reduced_Naive_TableSize(bw,cutoff) +
                      Reduced_SpharmonicTableSize(bw,cutoff)

   Then the workspace needs to be this large:

   2 * legendre_size  +
   8 * (bw*bw)  + 10*bw +
   4 * (bw*bw) + 2*bw

   for a total of

   2 * legendre_size  +
   12 * (bw*bw) + 12*bw ;
*/
void Conv2Sphere_semi_memo(double* rdata, double* idata, double* rfilter, double* ifilter, double* rres, double* ires,
                           const int bw, double* workspace) {
    int size = 2 * bw;
    int cutoff = bw; // TODO move to args?
    int legendre_size = Reduced_Naive_TableSize(bw, cutoff) + Reduced_SpharmonicTableSize(bw, cutoff);

    double* spharmonic_result_space = workspace; // needs legendre_size
    double* transpose_spharmonic_result_space = spharmonic_result_space + legendre_size; // needs legendre_size

    double* frres = transpose_spharmonic_result_space + legendre_size; // needs (bw * bw)
    double* fires = frres + (bw * bw);                                 // needs (bw * bw)
    double* trres = fires + (bw * bw);                                 // needs (bw * bw)
    double* tires = trres + (bw * bw);                                 // needs (bw * bw)
    double* filtrres = tires + (bw * bw);                              // needs bw
    double* filtires = filtrres + bw;                                  // needs bw
    double* scratchpad = filtires + bw;                                // needs (8 * bw^2) + (10 * bw)

    double* weights = (double*)malloc(sizeof(double) * 4 * bw);
    makeweights(bw, weights);

    // Make DCT plans. Note that I will be using the GURU interface to execute these plans within the routines
    fftw_plan DCT_plan = fftw_plan_r2r_1d(size, weights, rdata, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_plan inv_DCT_plan = fftw_plan_r2r_1d(size, weights, rdata, FFTW_REDFT01, FFTW_ESTIMATE);

    // fftw "preamble"
    // Note that FFT plan places the output in a transposed array
    int rank = 1;
    fftw_iodim dims[rank];
    dims[0].n = size;
    dims[0].is = 1;
    dims[0].os = size;

    int howmany_rank = 1;
    fftw_iodim howmany_dims[howmany_rank];
    howmany_dims[0].n = size;
    howmany_dims[0].is = size;
    howmany_dims[0].os = 1;

    fftw_plan FFT_plan = fftw_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, rdata, idata, workspace,
                                                  workspace + (4 * bw * bw), FFTW_ESTIMATE);

    // Note that FFT plans assumes that I'm working with a transposed array, e.g. the inputs for a length 2*bw transform
    // are placed every 2*bw apart, the output will be consecutive entries in the array

    rank = 1;
    dims[0].n = size;
    dims[0].is = size;
    dims[0].os = 1;

    howmany_rank = 1;
    howmany_dims[0].n = size;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = size;

    fftw_plan inv_FFT_plan = fftw_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, rdata, idata, workspace,
                                                      workspace + (4 * bw * bw), FFTW_ESTIMATE);

    // precompute the associated Legendre fcts
    double** spharmonic_pml_table = Spharmonic_Pml_Table(bw, spharmonic_result_space, scratchpad);
    double** transpose_spharmonic_pml_table =
        Transpose_Spharmonic_Pml_Table(spharmonic_pml_table, bw, transpose_spharmonic_result_space);

    FST_semi_memo(rdata, idata, frres, fires, bw, spharmonic_pml_table, scratchpad, 1, bw, &DCT_plan, &FFT_plan, weights);
    FZT_semi_memo(rfilter, ifilter, filtrres, filtires, bw, spharmonic_pml_table[0], scratchpad, 1, &DCT_plan, weights);

    TransMult(frres, fires, filtrres, filtires, trres, tires, bw);

    InvFST_semi_memo(trres, tires, rres, ires, bw, transpose_spharmonic_pml_table, scratchpad, 1, bw, &inv_DCT_plan,
                     &inv_FFT_plan);

    free(weights);
    free(spharmonic_pml_table);
    free(transpose_spharmonic_pml_table);

    fftw_destroy_plan(inv_FFT_plan);
    fftw_destroy_plan(FFT_plan);
    fftw_destroy_plan(inv_DCT_plan);
    fftw_destroy_plan(DCT_plan);
}
