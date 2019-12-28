/*
    Source code to test full spherical harmonic transform using the seminaive and naive algorithms
    coded up during October, 1995.

    In its current state, assuming that will seminaive at *all* orders.
    If you wish to change this, modify the `cutoff` variable in this file, i.e.
    `cutoff =` at what order to switch from semi-naive to naive algorithm.

    Note: Code will precompute everything in memory before doing transform.

    Idea is to record the execution time for the spectral -> grid -> spectral roundtrip,
    after precomputations have been done.

    The strategy is to generate random coefficients representing a complete spherical harmonic series
    of funcions Y(m,l). The ordering of the these coefficients is assumed to be the same as that of the
    output of the FST_seminaive() routine, namely:

    f(0,0) f(0,1) f(0,2)       ...       f(0,bw-1)
           f(1,1) f(1,2)       ...       f(1,bw-1)
    etc.
                         f(bw-2,bw-2)    f(bw-2,bw-1)
                                         f(bw-1,bw-1)
                                         f(-(bw-1),bw-1)
                         f(-(bw-2),bw-2) f(-(bw-2),bw-1)
    etc.
                  f(-2,2)      ...       f(-2,bw-1)
          f(-1,1) f(-1,2)      ...       f(-1,bw-1)

    This means that there are (bw*bw) coefficients.

    Once the coefficients are generated, the corresponding function is synthesized using InvFSTSemiMemo(),
    then transformed (analyzed) using FSTSemiMemo(). Timing data is printed.

    Sample call:

    test_s2_semi_memo bw loops [error_file]

    Appropriate timing data will be printed out.

    NOTE: In this program, the coefficients generated are such that the grid points (sample values) produced are *real*.
    The routines InvFSTSemiFly() and FSTSemiFly() will take advantage of this. If you wish to change this,
    change the 7th argument of each function further down from "1" to "0".
    This is also documented in the file FST_semi_fly.c.
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
    if (argc < 3) {
        fprintf(stdout, "Usage: test_s2_semi_memo bw loops [error_file]\n");
        exit(0);
    }

    int bw = atoi(argv[1]);
    int loops = atoi(argv[2]);

    int size = 2 * bw;

    // TODO more tests with dataformat 0 and cutoff
    int cutoff = bw / 4; // = bw // seminaive all orders
    DataFormat data_format = COMPLEX;

    double* workspace = (double*)malloc(sizeof(double) * ((8 * (bw * bw)) + (16 * bw)));
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

    time_t seed;
    time(&seed);
    srand48(seed); // init random generator

    double* weights = (double*)malloc(sizeof(double) * 4 * bw);
    double* rdata = (double*)malloc(sizeof(double) * (size * size));
    double* idata = (double*)malloc(sizeof(double) * (size * size));

    // Make DCT plans. Note that I will be using the GURU interface to execute these plans within the routines
    fftw_plan DCT_plan = fftw_plan_r2r_1d(2 * bw, weights, rdata, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_plan inv_DCT_plan = fftw_plan_r2r_1d(2 * bw, weights, rdata, FFTW_REDFT01, FFTW_ESTIMATE);

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

    // Note that FFT plans assumes that I'm working with a transposed array, e.g. the inputs for a length 2*bw transform
    // are placed every 2*bw apart, the output will be consecutive entries in the array

    rank = 1;
    dims[0].n = 2 * bw;
    dims[0].is = 2 * bw;
    dims[0].os = 1;

    howmany_rank = 1;
    howmany_dims[0].n = 2 * bw;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = 2 * bw;

    fftw_plan inv_FFT_plan = fftw_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, rdata, idata, workspace,
                                                      workspace + (4 * bw * bw), FFTW_ESTIMATE);

    GenerateWeightsForDLT(bw, weights);

    double* rcoeffs = (double*)malloc(sizeof(double) * (bw * bw));
    double* icoeffs = (double*)malloc(sizeof(double) * (bw * bw));

    double* rresult = (double*)malloc(sizeof(double) * (bw * bw));
    double* iresult = (double*)malloc(sizeof(double) * (bw * bw));

    double* relerror = (double*)malloc(sizeof(double) * loops);
    double* curmax = (double*)malloc(sizeof(double) * loops);

    double fwd_time = 0.0;
    double inv_time = 0.0;
    double total_error = 0.0;
    double total_relerror = 0.0;

    fprintf(stdout, "about to enter loop\n\n");
    for (int i = 0; i < loops; ++i) {

        // generate spherical harmonic coefficients
        for (int m = 0; m < bw; ++m)
            for (int l = m; l < bw; ++l) {
                double x = 2.0 * (drand48() - 0.5);
                double y = 2.0 * (drand48() - 0.5);

                int index = IndexOfHarmonicCoeff(m, l, bw);
                rcoeffs[index] = x;
                icoeffs[index] = y;

                index = IndexOfHarmonicCoeff(-m, l, bw);
                rcoeffs[index] = pow(-1.0, m) * x;
                icoeffs[index] = pow(-1.0, m + 1) * y;
            }

        // have to zero out the m=0 coefficients, since those are real
        for (int m = 0; m < bw; ++m)
            icoeffs[m] = 0.0;

        double time_start = csecond();
        // inverse spherical transform
        InvFSTSemiMemo(rcoeffs, icoeffs, rdata, idata, bw, trans_seminaive_naive_table, workspace, data_format, cutoff,
                       &inv_DCT_plan, &inv_FFT_plan);

        double duration = csecond() - time_start;
        inv_time += duration;
        fprintf(stdout, "inv time \t = %.4e\n", duration);

        time_start = csecond();
        // forward spherical transform
        FSTSemiMemo(rdata, idata, rresult, iresult, bw, seminaive_naive_table, workspace, data_format, cutoff,
                    &DCT_plan, &FFT_plan, weights);

        duration = csecond() - time_start;
        fwd_time += duration;
        fprintf(stdout, "forward time \t = %.4e\n", duration);

        // compute the error
        relerror[i] = 0.0;
        curmax[i] = 0.0;

        for (int j = 0; j < (bw * bw); ++j) {
            double rtmp = rresult[j] - rcoeffs[j];
            double itmp = iresult[j] - icoeffs[j];
            double origmag = sqrt((rcoeffs[j] * rcoeffs[j]) + (icoeffs[j] * icoeffs[j]));
            double tmpmag = sqrt((rtmp * rtmp) + (itmp * itmp));
            relerror[i] = fmax(relerror[i], tmpmag / (origmag + pow(10.0, -50.0)));
            curmax[i] = fmax(curmax[i], tmpmag);
        }

        fprintf(stdout, "r-o error\t = %.12f\n", curmax[i]);
        fprintf(stdout, "(r-o)/o error\t = %.12f\n\n", relerror[i]);

        total_error += curmax[i];
        total_relerror += relerror[i];
    }

    double avg_error = total_error / loops;
    double avg_relerror = total_relerror / loops;
    double stddev_error = 0.0;
    double stddev_relerror = 0.0;

    for (int i = 0; i < loops; ++i) {
        stddev_error += pow(avg_error - curmax[i], 2.0);
        stddev_relerror += pow(avg_relerror - relerror[i], 2.0);
    }

    if (loops != 1) {
        stddev_error = sqrt(stddev_error / (loops - 1));
        stddev_relerror = sqrt(stddev_relerror / (loops - 1));
    }

    fprintf(stdout, "Program: test_s2_semi_memo\n");
    fprintf(stdout, "Bandwidth = %d\n", bw);

    const char* timing_object = "cpu";
#ifdef WALLCLOCK
    timing_object = "wall";
#endif

    double total_time = inv_time + fwd_time;

    fprintf(stdout, "Total elapsed %s time :\t\t %.4e seconds.\n", timing_object, total_time);
    fprintf(stdout, "Average %s forward per iteration:\t %.4e seconds.\n", timing_object, fwd_time / loops);
    fprintf(stdout, "Average %s inverse per iteration:\t %.4e seconds.\n", timing_object, inv_time / loops);

    fprintf(stdout, "Average r-o error:\t\t %.4e\t", total_error / loops);
    fprintf(stdout, "std dev: %.4e\n", stddev_error);
    fprintf(stdout, "Average (r-o)/o error:\t\t %.4e\t", total_relerror / loops);
    fprintf(stdout, "std dev: %.4e\n\n", stddev_relerror);

    if (argc == 4) {
        FILE* errorsfp = fopen(argv[3], "w");

        for (int m = 0; m < bw; ++m) {
            for (int l = m; l < bw; ++l) {
                int index = IndexOfHarmonicCoeff(m, l, bw);
                fprintf(errorsfp, "index = %d\t m = %d\tl = %d\t%.10f  %.10f\n", index, m, l,
                        fabs(rcoeffs[index] - rresult[index]), fabs(icoeffs[index] - iresult[index]));

                index = IndexOfHarmonicCoeff(-m, l, bw);
                fprintf(errorsfp, "index = %d\t m = %d\tl = %d\t%.10f  %.10f\n", index, -m, l,
                        fabs(rcoeffs[index] - rresult[index]), fabs(icoeffs[index] - iresult[index]));
            }
        }

        fclose(errorsfp);
    }

    fftw_destroy_plan(inv_FFT_plan);
    fftw_destroy_plan(FFT_plan);
    fftw_destroy_plan(inv_DCT_plan);
    fftw_destroy_plan(DCT_plan);

    free(weights);
    free(trans_seminaive_naive_table);
    free(seminaive_naive_table);
    free(curmax);
    free(relerror);
    free(workspace);
    free(trans_seminaive_naive_tablespace);
    free(seminaive_naive_tablespace);
    free(iresult);
    free(rresult);
    free(idata);
    free(rdata);
    free(icoeffs);
    free(rcoeffs);

    return 0;
}
