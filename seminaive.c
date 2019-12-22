/*
    Source code for computing the Legendre transform where projections are carried out in cosine space,
    i.e., the "seminaive" algorithm.

    For a description, see the related paper or Sean's thesis.
*/

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <fftw3.h>

#include "cospmls.h"

/* 
    Computes the inverse Legendre transform using the transposed seminaive algorithm.
    Note that because the Legendre transform is orthogonal, the inverse can be
    computed by transposing the matrix formulation of the problem.

    The forward transform looks like

    l = PCWf

    where f is the data vector, W is a quadrature matrix, C is a cosine transform matrix,
    P is a matrix full of coefficients of the cosine series representation of each Pml
    function P(m,m) P(m,m+1) ... P(m,bw-1), and l is the (associated) Legendre series 
    representation of f.

    So to do the inverse, you do

    f = trans(C) trans(P) l

    so you need to transpose the matrix P from the forward transform and then do a cosine
    series evaluation. No quadrature matrix is necessary. If order `m` is odd, then there
    is also a sin factor that needs to be accounted for.

    Note that this function was written to be part of a full spherical harmonic transform,
    so a lot of precomputation has been assumed.

    Arguments description:

    coeffs - a double pointer to an array of length `bw-m` containing associated Legendre
             series coefficients. Assumed that first entry contains the P(m,m) coefficient;
    bw - problem bandwidth;
    m - order of the associated Legendre functions;
   result - a double pointer to an array of `2*bw` samples representing the evaluation
            of the Legendre series at `2*bw` Chebyshev nodes;
   trans_cos_pml_table - double pointer to array representing
                         the linearized form of trans(P) above.
                         See cospmls.{h,c} for a description
                         of the function `Transpose_CosPmlTableGen()`
                         which generates this array;
   sin_values - when `m` is odd, need to factor in the sin(x) that
                is factored out of the generation of the values
                in trans(P);
   workspace - a double array of size `2*bw` -> temp space involving
               intermediate array;
   fplan - pointer to `fftw_plan` with input array being fcos
           and output being result; I'll probably use the
           guru interface to execute - that way I can apply the
           same plan to different arrays; the plan should be:
           `fftw_plan_r2r_1d(2*bw, fcos, result, FFTW_REDFT01, FFTW_ESTIMATE)`.
*/
// TODO remove memset?
void InvSemiNaiveReduced(double* coeffs, const int bw, const int m, double* result, double* trans_cos_pml_table,
                         double* sin_values, double* workspace, fftw_plan* fplan) {
    double* fcos = workspace;

    // for paranoia, zero out arrays
    memset(fcos, 0, sizeof(double) * 2 * bw);
    memset(result, 0, sizeof(double) * 2 * bw);

    double* trans_tableptr = trans_cos_pml_table;
    double* p = trans_cos_pml_table;

    /* 
        main loop - compute each value of fcos

        Note that all zeroes have been stripped out of the
        `trans_cos_pml_table`, so indexing is somewhat complicated.
    */

    for (int i = 0; i < bw; ++i) {
        if (i == (bw - 1)) {
            if (m % 2) {
                fcos[bw - 1] = 0.0;
                break;
            }
        }

        double* assoc_offset;
        if (i > m)
            assoc_offset = coeffs + (i - m) + (m % 2);
        else
            assoc_offset = coeffs + (i % 2);

        double fcos0 = 0.0;
        double fcos1 = 0.0;
        double fcos2 = 0.0;
        double fcos3 = 0.0;

        int rowsize = Transpose_RowSize(i, m, bw);
        int j = 0;
        // TODO one item per cycle
        for (; j < rowsize % 4; ++j)
            fcos0 += assoc_offset[2 * j] * trans_tableptr[j];

        for (; j < rowsize; j += 4) {
            fcos0 += assoc_offset[2 * j] * trans_tableptr[j];
            fcos1 += assoc_offset[2 * (j + 1)] * trans_tableptr[j + 1];
            fcos2 += assoc_offset[2 * (j + 2)] * trans_tableptr[j + 2];
            fcos3 += assoc_offset[2 * (j + 3)] * trans_tableptr[j + 3];
        }
        fcos[i] = fcos0 + fcos1 + fcos2 + fcos3;

        trans_tableptr += rowsize;
    }

    /*
        now we have the cosine series for the result,
        so now evaluate the cosine series at `2*bw` Chebyshev nodes
    */

    // scale coefficients prior to taking inverse DCT
    double fudge = 0.5 / sqrt(bw);
    for (int j = 1; j < 2 * bw; ++j)
        fcos[j] *= fudge;
    fcos[0] /= sqrt(2. * bw);

    // take the inverse dct
    // Note that I am using the guru interface
    fftw_execute_r2r(*fplan, fcos, result);

    // if m is odd, then need to multiply by sin(x) at Chebyshev nodes
    if (m % 2) {
        for (int j = 0; j < (2 * bw); j++)
            result[j] *= sin_values[j];
    }

    trans_tableptr = p;
}

/*
    Computes the Legendre transform of data. This function has been designed to be
    a component in a full spherical harmonic transform.

    data - pointer to double array of size `2*bw` containing
           function to be transformed. Assumes sampling at Chebyshev nodes;
   bw   - bandwidth of the problem;
   m    - order of the problem.  0 <= `m` < `bw`;
   result - pointer to double array of length bw for returning computed Legendre coefficients.
            Contains `bw-m` coeffs, with the <f,P(m,m)> coefficient located in `result[0]`;
   cos_pml_table - a pointer to an array containing the cosine series coefficients of
                   the Pmls (or Gmls) for this problem. This table can be compute using the
                   `CosPmlTableGen()` function, and the offset for a particular Pml can be found
                    by calling the function `TableOffset()`. The size of the table is computed using
                    the `TableSize()` function. Note that since the cosine series are always 
                    zero-striped, the zeroes have been removed;
   weights - ptr to double array of size `4*bw` - this array holds the weights for both even
             (starting at `weights[0]`) and odd (`weights[2*bw]`) transforms;
   workspace - tmp space: pointer to double array of size `4*bw`;
   fplan - pointer to `fftw_plan` with input array being weighted_data and output being cos_data;
           I'll probably use the guru interface to execute; the plan should be:
           `fftw_plan_r2r_1d(2*bw, weighted_data, cos_data, FFTW_REDFT10, FFTW_ESTIMATE)`.
*/
void SemiNaiveReduced(double* data, const int bw, const int m, double* result, double* workspace,
                      double* cos_pml_table, double* weights, fftw_plan* fplan) {
    int n = 2 * bw;

    double* weighted_data = workspace;
    double* cos_data = weighted_data + (2 * bw);

    // for paranoia, zero out the result array
    memset(result, 0, sizeof(double) * (bw - m));

    // apply quadrature weights to the data and
    // compute the cosine transform
    if (m % 2)
        for (int i = 0; i < n; ++i)
            weighted_data[i] = data[i] * weights[2 * bw + i];
    else
        for (int i = 0; i < n; ++i)
            weighted_data[i] = data[i] * weights[i];


    // smooth the weighted signal

    fftw_execute_r2r(*fplan, weighted_data, cos_data);

    // normalize
    cos_data[0] *= M_SQRT1_2;
    double fudge = 1. / sqrt(2. * n);
    for (int j = 0; j < n; ++j)
        cos_data[j] *= fudge;

    /*
        Do the projections.
        Note that the `cos_pml_table` has had all the zeroes
        stripped out so the indexing is complicated somewhat
    */

    // TODO bring back the original one
    /******** this is the original loop

    int toggle = 0 ;
    for (int i=m; i<bw; ++i)
    {
    double* pml_ptr = cos_pml_table + TableOffset(m,i);

    if ((m % 2) == 0)
    {
    for (int j=0; j<(i/2)+1; ++j)
    result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
    }
    else
    {
    if (((i-m) % 2) == 0)
    {
    for (int j=0; j<(i/2)+1; ++j)
    result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
    }
    else
    {
    for (int j=0; j<(i/2); ++j)
    result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
    }
    }

    toggle = (toggle+1) % 2;
    }

    *****/

    /******** this is the new loop *********/
    int toggle = 0;
    for (int i = m; i < bw; ++i) {
        double* pml_ptr = cos_pml_table + TableOffset(m, i);

        double result0 = 0.0;
        double result1 = 0.0;
        double result2 = 0.0;
        double result3 = 0.0;

        int j = 0;
        for (; j < ((i / 2) % 4); ++j)
            result0 += cos_data[(2 * j) + toggle] * pml_ptr[j];

        for (; j < (i / 2); j += 4) {
            result0 += cos_data[(2 * j) + toggle] * pml_ptr[j];
            result1 += cos_data[(2 * (j + 1)) + toggle] * pml_ptr[j + 1];
            result2 += cos_data[(2 * (j + 2)) + toggle] * pml_ptr[j + 2];
            result3 += cos_data[(2 * (j + 3)) + toggle] * pml_ptr[j + 3];
        }

        if ((((i - m) % 2) == 0) || ((m % 2) == 0))
            result0 += cos_data[(2 * (i / 2)) + toggle] * pml_ptr[(i / 2)];

        result[i - m] = result0 + result1 + result2 + result3;

        toggle = (toggle + 1) % 2;
    }
}
