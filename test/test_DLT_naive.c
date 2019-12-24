/*
    test_naive.c - top-level code for computing forward and inversee
    Legendre transforms using the naive algorithm; will return timing
    and error information

    m     - order of the problem
    bw    - bandwidth
    loops - number of loops thru timed portion of code. Intended
            to reduce noise due to multiprocessing and
            discretization errors

    Sample calls:

    test_naive m bw loops

    test_naive 0 32 10
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "legendre_polynomials/pml.h"
#include "legendre_transform/naive.h"
#include "legendre_transform/weights.h"

#include "util/csecond.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stdout, "Usage: test_DLT_naive m bw loops\n");
        return 0;
    }

    int m = atoi(argv[1]);
    int bw = atoi(argv[2]);
    int loops = atoi(argv[3]);

    double* samples = (double*)malloc(sizeof(double) * 2 * bw);
    double* coeffs = (double*)malloc(sizeof(double) * (bw - m));
    double* new_coeffs = (double*)malloc(sizeof(double) * (bw - m));

    double* relerror = (double*)malloc(sizeof(double) * loops);
    double* curmax = (double*)malloc(sizeof(double) * loops);

    double* plm = (double*)malloc(sizeof(double) * 2 * bw * (bw - m));
    double* workspace = (double*)malloc(sizeof(double) * 18 * bw);
    GeneratePmlTable(bw, m, plm, workspace);

    double* weights = (double*)malloc(sizeof(double) * 4 * bw);
    GenerateWeightsForDLT(bw, weights);

    long int seed;
    time(&seed);
    srand48(seed); // init random generator

    double sum_error = 0.0;
    double sum_relerror = 0.0;

    double total_time_f = 0.0;
    double total_time_i = 0.0;

    for (int k = 0; k < loops; ++k) {
        // generate random coefficients
        for (int i = 0; i < (bw - m); ++i)
            coeffs[i] = 2.0 * (drand48() - 0.5);

        double time_start = csecond();
        // inverse naive transform
        InvDLTNaive(coeffs, bw, m, samples, plm);
        total_time_i += (csecond() - time_start);

        time_start = csecond();
        // forward naive transform
        DLTNaive(samples, bw, m, weights, new_coeffs, plm, workspace);
        total_time_f += (csecond() - time_start);

        // tally up the error between the original coefficients and the new ones
        relerror[k] = 0.0;
        curmax[k] = 0.0;

        for (int i = 0; i < bw - m; ++i) {
            double tmp_error = fabs(coeffs[i] - new_coeffs[i]);
            double tmp_relerror = tmp_error / (fabs(coeffs[i]) + pow(10.0, -50.0));
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

    free(curmax);
    free(relerror);
    free(workspace);
    free(weights);
    free(new_coeffs);
    free(coeffs);
    free(samples);

    return 0;
}
