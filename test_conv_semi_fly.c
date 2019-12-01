/*
  Source code to convolve two real-valued functions defined on the 2-sphere.
  Uses seminaive algorithms.

  WILL PRECOMPUTE DATA AS NEEDED, ON THE FLY (that sounds odd,
  "precompute ... the fly")

  Reads in a function and filter from files specified at
  shell level, and dumps output into a
  specified file.
  Both function and filter must be
  (size x size) arrays, where
  size = 2*bandwidth.

  User specifies bandwidth as fourth argument

  Sample call:

  test_conv_semi_fly signalFile filterFile convolveFile bandwidth

  test_conv_semi_fly s64.dat f64.dat c64.dat 64


  In this example, the signal and filter function samples are stored
  in the files

  s64.dat ( signal - for bandwidth = 64 )
  f64.dat ( filter - for bandwidth = 64 )

  s128.dat ( signal - for bandwidth = 128 )
  f128.dat ( filter - for bandwidth = 128 )

  The signal is a "noisey" bump on the sphere, and the filter
  is a smooth, symmetric bump at the north pole

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


  The *location* of the maximum value in the output file tells me
  where the bump is.

*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FST_semi_fly.h"
#include "fftw3.h"

int main(int argc, char** argv) {
    FILE* fp;
    int i, bw, size;
    double *rsignal, *isignal, *rfilter, *ifilter, *rresult, *iresult;
    double* workspace;

    if (argc < 5) {
        fprintf(stdout, "Usage: test_conv_semi_fly signal_file filter_file "
                        "output_file bw\n");
        exit(0);
    }

    bw = atoi(argv[4]);
    size = 2 * bw;

    /*
      allocate memory
    */
    rsignal = (double*)malloc(sizeof(double) * size * size);
    isignal = (double*)malloc(sizeof(double) * size * size);
    rfilter = (double*)malloc(sizeof(double) * size * size);
    ifilter = (double*)malloc(sizeof(double) * size * size);
    rresult = (double*)malloc(sizeof(double) * size * size);
    iresult = (double*)malloc(sizeof(double) * size * size);

    workspace = (double*)malloc(sizeof(double) * (14 * bw * bw + 26 * bw));

    /****
      At this point, check to see if all the memory has been
      allocated. If it has not, there's no point in going further.
      ****/

    if ((rsignal == NULL) || (isignal == NULL) || (rfilter == NULL) ||
        (ifilter == NULL) || (rresult == NULL) || (iresult == NULL) ||
        (workspace == NULL)) {
        perror("Error in allocating memory");
        exit(1);
    }

    /* read in signal and filter */
    fprintf(stdout, "Reading signal file...\n");
    fp = fopen(argv[1], "r");
    for (i = 0; i < size * size; i++)
        fscanf(fp, "%lf", rsignal + i);
    fclose(fp);

    fprintf(stdout, "Reading filter file...\n");
    fp = fopen(argv[2], "r");
    for (i = 0; i < size * size; i++)
        fscanf(fp, "%lf", rfilter + i);
    fclose(fp);

    /*
      since the data are strictly real-valued, I need to zero out
      the imaginary parts
    */
    memset(isignal, 0, sizeof(double) * size * size);
    memset(ifilter, 0, sizeof(double) * size * size);

    /* now convolve */
    fprintf(stdout, "Calling Conv2Sphere_semi_fly()\n");
    Conv2Sphere_semi_fly(rsignal, isignal, rfilter, ifilter, rresult, iresult,
                         bw, workspace);

    /* convolving real functions results in real output,
       so no need to write out the imaginary array */

    fprintf(stdout, "Writing output file...\n");
    fp = fopen(argv[3], "w");
    for (i = 0; i < size * size; i++)
        fprintf(fp, "%.16f\n", rresult[i]);
    fclose(fp);

    /* free memory */
    free(workspace);
    free(iresult);
    free(rresult);
    free(ifilter);
    free(rfilter);
    free(isignal);
    free(rsignal);

    return 0;
}
