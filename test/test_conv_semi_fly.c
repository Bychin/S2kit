/*
    Source code to convolve two real-valued functions defined on the 2-sphere. Uses seminaive algorithms.

    Will precompute data as needed, on the fly.

    Reads in a function and filter from files specified at shell level, and dumps output into a specified file.
    Both function and filter must be (size*size) arrays, where size = 2*bandwidth.

    Sample call:

    test_conv_semi_fly signalFile filterFile convolveFile bandwidth

    test_conv_semi_fly s64.dat f64.dat c64.dat 64


    In this example, the signal and filter function samples are stored in the files:

    s64.dat ( signal - for bandwidth = 64 )
    f64.dat ( filter - for bandwidth = 64 )

    s128.dat ( signal - for bandwidth = 128 )
    f128.dat ( filter - for bandwidth = 128 )

    The signal is a "noisey" bump on the sphere, and the filter is a smooth, symmetric bump at the north pole.

    The samples for each are in "latitude-major" format. I.e.
    (theta_0, phi_0)
    (theta_0, phi_1)
    (theta_0, phi_2)
    ...
    (theta_0, phi_{bw-1})
    (theta_1, phi_0)
    (theta_1, phi_1)
    ...
    (theta_{bw-1}, phi_{bw-1})

    where theta_k = pi*(2*j+1)/(4*bw)
          phi_j = 2*pi*k/(2*bw)

    The *location* of the maximum value in the output file tells me where the bump is.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "FST_semi_fly.h"

int main(int argc, char** argv) {

    if (argc < 5) {
        fprintf(stdout, "Usage: test_conv_semi_fly signal_file filter_file "
                        "output_file bw\n");
        exit(0);
    }

    int bw = atoi(argv[4]);
    int size = 2 * bw;

    double* rsignal = (double*)malloc(sizeof(double) * size * size);
    double* rfilter = (double*)malloc(sizeof(double) * size * size);

    // read signal and filter
    fprintf(stdout, "Reading signal file...\n");
    FILE* fp = fopen(argv[1], "r");
    for (int i = 0; i < size * size; ++i)
        fscanf(fp, "%lf", rsignal + i);
    fclose(fp);

    fprintf(stdout, "Reading filter file...\n");
    fp = fopen(argv[2], "r");
    for (int i = 0; i < size * size; ++i)
        fscanf(fp, "%lf", rfilter + i);
    fclose(fp);

    // the imaginary parts are zeros,
    // since the data are strictly real-valued
    double* isignal = (double*)calloc((size_t)(size * size), sizeof(double));
    double* ifilter = (double*)calloc((size_t)(size * size), sizeof(double));

    double* rresult = (double*)malloc(sizeof(double) * size * size);
    double* iresult = (double*)malloc(sizeof(double) * size * size);
    double* workspace = (double*)malloc(sizeof(double) * (14 * bw * bw + 26 * bw));

    fprintf(stdout, "Calling ConvOn2SphereSemiFly()\n");
    ConvOn2SphereSemiFly(rsignal, isignal, rfilter, ifilter, rresult, iresult, bw, workspace);

    // convolving real functions results in real output,
    // so no need to write the imaginary array
    fprintf(stdout, "Writing output file...\n");
    fp = fopen(argv[3], "w");
    for (int i = 0; i < size * size; ++i)
        fprintf(fp, "%.16f\n", rresult[i]);
    fclose(fp);

    free(workspace);
    free(iresult);
    free(rresult);
    free(ifilter);
    free(rfilter);
    free(isignal);
    free(rsignal);

    return 0;
}
