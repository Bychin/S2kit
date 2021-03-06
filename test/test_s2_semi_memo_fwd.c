/**
 * @file test_s2_semi_memo_fwd.c
 * @brief Example of source code to computie @b forward spherical harmonic transform using the
 * seminaive and naive algorithms.
 *
 * Computes the spherical coefficients of the function, whose sample values are given in the named
 * input file, and writes those coefficients in the named output file.
 *
 * The code will <b>pre-compute</b> associated Legendre functions before doing transform. This
 * @b won't a be part of what's being timed.
 *
 * In its current state, assuming that will seminaive at @b all orders. If you wish to change this,
 * modify the @p cutoff variable in this file, i.e. <tt>cutoff = ...</tt> at what order to switch
 * from seminaive to naive algorithm.
 *
 * Sample call:
 * @code
 * test_s2_semi_memo_fwd sample_file output_file  bw [output_format]
 * test_s2_semi_memo_fwd y31_bw8.dat y31_coef.dat 8
 * @endcode
 *
 * The format of the input sample file will be an interleaved real/imaginary parts of the function
 * samples arranged in "latitude-major" format, i.e. the function will be sampled in this order:
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
 *
 * @p output_format is an optional argument:\n
 * @c 0 - coefficients in "code" order (interleaved real/imaginary), i.e. suitable for
 * test_s2_semi_memo_inv.c;\n
 * @c 1 - coefficients in prettier "human" order (verbose), i.e. if <tt>f_{l,m}</tt> is the
 * coefficient of degree @c l, order @c m, then the coefficients will be arranged this way:
 * @code
 * f_{0,0},
 * f_{1,-1}, f_{1,0},  f_{1,1},
 * f_{2,-2}, f_{2,-1}, f_{2,0}, f_{2,1}, f_{2,2},
 * ...
 * E.g. l = 2  m = 1  2.3 + 6 I.
 * @endcode
 * The default output format is "code" order. To help you out, there is a function
 * IndexOfHarmonicCoeff() which returns the array index of the coefficient
 * <tt>f_{l,m}</tt>.
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
#include "s2kit/util.h"

#include "util/csecond.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stdout, "Usage: test_s2_semi_memo_fwd sample_file output_file bw "
                        "[output_format]\n");
        exit(0);
    }

    int bw = atoi(argv[3]);

    int size = 2 * bw;
    int cutoff = bw; // seminaive all orders // TODO?
    DataFormat data_format = COMPLEX;

    double* workspace = (double*)malloc(sizeof(double) * ((8 * (bw * bw)) + (7 * bw)));
    double* seminaive_naive_tablespace = (double*)malloc(
        sizeof(double) * (Reduced_Naive_TableSize(bw, cutoff) + Reduced_SpharmonicTableSize(bw, cutoff)));

    // precompute the Legendres (that's what memo suffix for)
    fprintf(stdout, "Generating seminaive_naive tables...\n");
    double** seminaive_naive_table = SemiNaive_Naive_Pml_Table(bw, cutoff, seminaive_naive_tablespace, workspace);

    double* weights = (double*)malloc(sizeof(double) * 4 * bw);
    double* rdata = (double*)malloc(sizeof(double) * (size * size));
    double* idata = (double*)malloc(sizeof(double) * (size * size));

    // Make DCT plan. Note that I will be using the GURU interface to execute these plans within the routines
    fftw_plan DCT_plan = fftw_plan_r2r_1d(2 * bw, weights, rdata, FFTW_REDFT10, FFTW_ESTIMATE);

    // fftw "preamble"
    // Note that FFT plan places the output in a transposed array
    int rank = 1;
    fftw_iodim dims[rank];
    dims[0].n = 2 * bw;
    dims[0].is = 1;
    dims[0].os = 2 * bw;

    int howmany_rank = 1;
    fftw_iodim howmany_dims[howmany_rank];
    howmany_dims[0].n = 2 * bw;
    howmany_dims[0].is = 2 * bw;
    howmany_dims[0].os = 1;

    fftw_plan FFT_plan = fftw_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, rdata, idata, workspace,
                                                  workspace + (4 * bw * bw), FFTW_ESTIMATE);

    GenerateWeightsForDLT(bw, weights);

    // read samples
    FILE* fp = fopen(argv[1], "r");
    for (int i = 0; i < size * size; ++i) {
        fscanf(fp, "%lf", rdata + i); // the real part of the sample
        fscanf(fp, "%lf", idata + i); // the imaginary part
    }
    fclose(fp);

    double* rcoeffs = (double*)malloc(sizeof(double) * (bw * bw));
    double* icoeffs = (double*)malloc(sizeof(double) * (bw * bw));

    double time_start = csecond();
    // forward spherical transform
    FSTSemiMemo(rdata, idata, rcoeffs, icoeffs, bw, seminaive_naive_table, workspace, data_format, cutoff, &DCT_plan,
                &FFT_plan, weights);

    fprintf(stdout, "forward time \t = %.4e\n", csecond() - time_start);
    fprintf(stdout, "about to write out coefficients\n");

    // choose format for writing coefficients
    int order = 0; // code format
    if (argc == 5)
        order = atoi(argv[4]); // human-readable format

    fp = fopen(argv[2], "w");
    if (order == 0)
        for (int i = 0; i < bw * bw; ++i)
            fprintf(fp, "%.15f\n%.15f\n", rcoeffs[i], icoeffs[i]);
    else
        for (int l = 0; l < bw; ++l)
            for (int m = -l; m < l + 1; ++m) {
                int index = IndexOfHarmonicCoeff(m, l, bw);
                fprintf(fp, "l = %d\t m = %d\t %.15f + %.15f I\n", l, m, rcoeffs[index], icoeffs[index]);
            }
    fclose(fp);

    fprintf(stdout, "finished writing coefficients\n");

    fftw_destroy_plan(FFT_plan);
    fftw_destroy_plan(DCT_plan);

    free(workspace);
    free(seminaive_naive_table);
    free(seminaive_naive_tablespace);
    free(weights);
    free(icoeffs);
    free(rcoeffs);
    free(idata);
    free(rdata);

    return 0;
}
