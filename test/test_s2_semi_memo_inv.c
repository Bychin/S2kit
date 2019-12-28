/**
 * @file test_s2_semi_memo_inv.c
 * @brief Example of source code to computie @b inverse spherical harmonic transform using the
 * seminaive and naive algorithms.
 *
 * Computed function samples using spherical coefficients of the function, which are given in the
 * named input file, and writes those samples in the named output file.
 *
 * The code will <b>pre-compute</b> associated Legendre functions before doing inverse transform.
 * This @b won't a be part of what's being timed.
 *
 * In its current state, assuming that will seminaive at @b all orders. If you wish to change this,
 * modify the @p cutoff variable in this file, i.e. <tt>cutoff = ...</tt> at what order to switch
 * from seminaive to naive algorithm.
 *
 * Sample call:
 * @code
 * test_s2_semi_memo_inv coeffs_file   output_file     bw
 * test_s2_semi_memo_inv coeff_bw8.dat samples_bw8.dat 8
 * @endcode
 *
 * The format of the input coefficient file is an interleaved real/imaginary parts of the
 * coefficients, where the coefficients are given in "code" order, as defined in
 * test_s2_semi_memo.c:
 * @code
 * f(0,0) f(0,1) f(0,2)       ...       f(0,bw-1)
 *        f(1,1) f(1,2)       ...       f(1,bw-1)
 * etc.
 *                      f(bw-2,bw-2)    f(bw-2,bw-1)
 *                                      f(bw-1,bw-1)
 *                                      f(-(bw-1),bw-1)
 *                      f(-(bw-2),bw-2) f(-(bw-2),bw-1)
 * etc.
 *               f(-2,2)      ...       f(-2,bw-1)
 *       f(-1,1) f(-1,2)      ...       f(-1,bw-1)
 * @endcode
 * To help you out, there is a function IndexOfHarmonicCoeff() which returns the array
 * index of the coefficient <tt>f_{l,m}</tt>.
 *
 * The format of the output sample file is an interleaved real/imaginary parts of the function
 * samples arranged in "latitude-major" format, i.e. the function is sampled in this order:
 * @code
 * (theta_0,      phi_0)
 * (theta_0,      phi_1)
 * (theta_0,      phi_2)
 *           ...
 * (theta_0,      phi_{bw-1})
 * (theta_1,      phi_0)
 * (theta_1,      phi_1)
 *           ...
 * (theta_{bw-1}, phi_{bw-1})
 *
 * where theta_k = pi*(2*j+1)/(4*bw)
 *         phi_j = 2*pi*k/(2*bw)
 * @endcode
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <fftw3.h>

#include "s2kit/FST_semi_memo.h"
#include "s2kit/cospml.h"
#include "s2kit/weights.h"

#include "util/csecond.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stdout, "Usage: test_s2_semi_memo_inv coeffs_file output_file bw\n");
        exit(0);
    }

    int bw = atoi(argv[3]);

    int size = 2 * bw;
    int cutoff = bw; // seminaive all orders
    DataFormat data_format = COMPLEX;

    double* workspace = (double*)malloc(sizeof(double) * ((8 * (bw * bw)) + (10 * bw)));
    double* seminaive_naive_tablespace = (double*)malloc(
        sizeof(double) * (Reduced_Naive_TableSize(bw, cutoff) + Reduced_SpharmonicTableSize(bw, cutoff)));
    double* trans_seminaive_naive_tablespace = (double*)malloc(
        sizeof(double) * (Reduced_Naive_TableSize(bw, cutoff) + Reduced_SpharmonicTableSize(bw, cutoff)));

    // precompute the Legendres (that's what memo suffix for)
    fprintf(stdout, "Generating seminaive_naive tables...\n");
    double** seminaive_naive_table = SemiNaive_Naive_Pml_Table(bw, cutoff, seminaive_naive_tablespace, workspace);

    fprintf(stdout, "Generating trans_seminaive_naive tables...\n");
    double** trans_seminaive_naive_table = Transpose_SemiNaive_Naive_Pml_Table(
        seminaive_naive_table, bw, cutoff, trans_seminaive_naive_tablespace, workspace);

    double* weights = (double*)malloc(sizeof(double) * 4 * bw);
    double* rdata = (double*)malloc(sizeof(double) * (size * size));
    double* idata = (double*)malloc(sizeof(double) * (size * size));

    // Make inverse DCT plan. Note that I will be using the GURU interface to execute these plans within the routines
    fftw_plan inv_DCT_plan = fftw_plan_r2r_1d(2 * bw, weights, rdata, FFTW_REDFT01, FFTW_ESTIMATE);

    // fftw "preamble"
    // Note that FFT plans assumes that I'm working with a transposed array, e.g. the inputs for a length 2*bw transform
    // are placed every 2*bw apart, the output will be consecutive entries in the array

    int rank = 1;
    fftw_iodim dims[rank];
    dims[0].n = 2 * bw;
    dims[0].is = 2 * bw;
    dims[0].os = 1;

    int howmany_rank = 1;
    fftw_iodim howmany_dims[howmany_rank];
    howmany_dims[0].n = 2 * bw;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = 2 * bw;

    fftw_plan inv_FFT_plan = fftw_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, rdata, idata, workspace,
                                                      workspace + (4 * bw * bw), FFTW_ESTIMATE);

    GenerateWeightsForDLT(bw, weights);

    double* rcoeffs = (double*)malloc(sizeof(double) * (bw * bw));
    double* icoeffs = (double*)malloc(sizeof(double) * (bw * bw));

    // read coefficients
    FILE* fp = fopen(argv[1], "r");
    for (int i = 0; i < bw * bw; ++i) {
        fscanf(fp, "%lf", rcoeffs + i); // the real part of the coefficient
        fscanf(fp, "%lf", icoeffs + i); // the imaginary part
    }
    fclose(fp);

    double time_start = csecond();
    // inverse spherical transform
    InvFSTSemiMemo(rcoeffs, icoeffs, rdata, idata, bw, trans_seminaive_naive_table, workspace, data_format, cutoff,
                   &inv_DCT_plan, &inv_FFT_plan);

    fprintf(stderr, "inv time \t = %.4e\n", csecond() - time_start);
    fprintf(stdout, "about to write out samples\n");

    fp = fopen(argv[2], "w");
    for (int i = 0; i < size * size; ++i)
        fprintf(fp, "%.15f\n%.15f\n", rdata[i], idata[i]);
    fclose(fp);

    fprintf(stdout, "finished writing samples\n");

    fftw_destroy_plan(inv_FFT_plan);
    fftw_destroy_plan(inv_DCT_plan);

    free(trans_seminaive_naive_table);
    free(seminaive_naive_table);
    free(trans_seminaive_naive_tablespace);
    free(seminaive_naive_tablespace);
    free(workspace);
    free(weights);
    free(idata);
    free(rdata);
    free(icoeffs);
    free(rcoeffs);

    return 0;
}
