/*
    test_semi.c - top-level code for computing Legendre transform using the
    seminaive algorithm; will do forward and inverse DLTs, and return
    error and timing results

    m     - order of the problem
    bw    - bandwidth
    loops - number of loops thru timed portion of code. Intended
            to reduce noise due to multiprocessing and
            discretization errors
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <fftw3.h>

#include "cospmls.h"
#include "makeweights.h"
#include "util/primitive.h"
#include "seminaive.h"

#include "util/csecond.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stdout, "Usage: test_semi m bw loops\n");
        exit(0);
    }

    int m = atoi(argv[1]);
    int bw = atoi(argv[2]);
    int loops = atoi(argv[3]);
    int n = 2 * bw;

    double* samples = (double*)malloc(sizeof(double) * n);
    double* coeffs = (double*)malloc(sizeof(double) * (bw - m));
    double* new_coeffs = (double*)malloc(sizeof(double) * (bw - m));

    double* relerror = (double*)malloc(sizeof(double) * loops);
    double* curmax = (double*)malloc(sizeof(double) * loops);

    double* eval_args = (double*)malloc(sizeof(double) * n);
    double* sin_values = (double*)malloc(sizeof(double) * n);

    // generate sin values (need if order of transform is odd)
    AcosOfChebyshevNodes(n, eval_args);
    for (int i = 0; i < n; ++i)
        sin_values[i] = sin(eval_args[i]);

    double* workspace = (double*)malloc(sizeof(double) * 9 * bw);
    double* weights = (double*)malloc(sizeof(double) * 4 * bw);

    // Note the arrays! Since I'll be using the guru fftw interface,
    // I don't care what the arrays are - I just want to make sure they're big enough
    fftw_plan forwardDCT = fftw_plan_r2r_1d(2 * bw, workspace, weights, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_plan inverseDCT = fftw_plan_r2r_1d(2 * bw, workspace, weights, FFTW_REDFT01, FFTW_ESTIMATE);

    makeweights(bw, weights);

    double* cos_pml_table = (double*)malloc(sizeof(double) * TableSize(m, bw));
    double* transpose_cos_pml_table = (double*)malloc(sizeof(double) * TableSize(m, bw));

    // precompute cosine series for Pml (Gml) functions, necessary for the forward transform
    CosPmlTableGen(bw, m, cos_pml_table, workspace);

    // take transpose of the precomputed cosine series, necessary for the inverse transform
    Transpose_CosPmlTableGen(bw, m, cos_pml_table, transpose_cos_pml_table);

    double sum_error = 0.0;
    double sum_relerror = 0.0;

    double total_time_f = 0.0;
    double total_time_i = 0.0;

    long int seed;
    time(&seed);
    srand48(seed); // init random generator

    for (int k = 0; k < loops; ++k) {
        // generate random coefficients
        for (int i = 0; i < (bw - m); ++i)
            coeffs[i] = 2.0 * (drand48() - 0.5);

        double time_start = csecond();
        // inverse semi-naive transform
        InvSemiNaiveReduced(coeffs, bw, m, samples, transpose_cos_pml_table, sin_values, workspace, &inverseDCT);
        total_time_i += (csecond() - time_start);

        time_start = csecond();
        // forward semi-naive transform
        SemiNaiveReduced(samples, bw, m, new_coeffs, workspace, cos_pml_table, weights, &forwardDCT);
        total_time_f += (csecond() - time_start);

        // tally up the error between the original coefficients and the new ones
        relerror[k] = 0.0;
        curmax[k] = 0.0;

        for (int j = 0; j < bw - m; ++j) {
            double tmp_error = fabs(coeffs[j] - new_coeffs[j]);
            double tmp_relerror = tmp_error / (fabs(coeffs[j]) + pow(10.0, -50.0));
            curmax[k] = fmax(curmax[k], tmp_error);
            relerror[k] = fmax(relerror[k], tmp_relerror);
        }

        sum_error += curmax[k];
        sum_relerror += relerror[k];
    }

    double avg_error = sum_error / loops;
    double avg_relerror = sum_relerror / loops;
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

    fprintf(stdout, "bw = %d\tm = %d\n", bw, m);
    fprintf(stdout, "loops = %d\n", loops);

    fprintf(stdout, "Average r-o error:\t\t %.4e\t", sum_error / loops);
    fprintf(stdout, "std dev: %.4e\n", stddev_error);
    fprintf(stdout, "Average (r-o)/o error:\t\t %.4e\t", sum_relerror / loops);
    fprintf(stdout, "std dev: %.4e\n\n", stddev_relerror);

    fprintf(stdout, "average forward time = %.4e\n", total_time_f / loops);
    fprintf(stdout, "average inverse time = %.4e\n", total_time_i / loops);

    fftw_destroy_plan(inverseDCT);
    fftw_destroy_plan(forwardDCT);

    free(weights);
    free(workspace);
    free(curmax);
    free(relerror);
    free(sin_values);
    free(eval_args);
    free(transpose_cos_pml_table);
    free(cos_pml_table);
    free(new_coeffs);
    free(coeffs);
    free(samples);

    return 0;
}
