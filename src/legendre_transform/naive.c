/*
    Source code to synthesize functions using a naive method based on recurrence.
    This is slow but does not require any precomputed functions, and is stable.
*/

#include "naive.h"

#include <math.h>
#include <string.h>

/*
    Computes the discrete Legendre transform of a function via summing naively.
    I.e. this is the forward discrete Legendre transform.

    bw        - bandwidth;
    m         - order;
    data      - a pointer to double array of size `2*bw` containing the sample points;
    result    - a pointer to double array of size `bw-m` which, at the conclusion of the routine,
                will contains the coefficients;
    plmtable  - a pointer to a double array of size `2*bw*(bw-m)` which contains the precomputed Plms,
                i.e. associated Legendre functions. E.g. Should be generated by a call to `PmlTableGen()`
                (see pml.c). Note that these Legendres are normalized with norm equal to 1;
    workspace - array of size `2*bw`;
*/
void Naive_AnalysisX(double* data, const int bw, const int m, double* weights, double* result,
                     double* plmtable, double* workspace) {
    int size = 2 * bw;

    /*
        Apply quadrature weights.

        We only have to differentiate between even and odd weights when doing something like seminaive,
        something which involves the dct. In this naive case, the parity of the order of the transform
        doesn't matter because we are not dividing by sin(x) when precomputing the Legendres
        (because we are not taking their dct). The plain weights are just fine.
    */
    double* wdata = workspace;
    for (int i = 0; i < size; ++i)
        wdata[i] = data[i] * weights[i];

    for (int i = 0; i < bw - m; ++i) {
        double sum = 0.;
        for (int j = 0; j < size; ++j)
            sum += wdata[j] * plmtable[j];

        result[i] = sum;

        plmtable += size;
    }
}

/*
    Synthesizes a function from a list of coefficients of a Legendre series.
    I.e. this is the inverse discrete Legendre transform.

    Function is synthesized at the `2*bw` Chebyshev nodes. Associated
    Legendre functions are assumed to be precomputed.

    bw       - bandwidth
    m        - order
    coeffs   - a pointer to double array of size `bw-m`. First coefficient is
               coefficient for Pmm
    result   - a pointer to double array of size `2*bw`; at the conclusion
               of the routine, this array will contain the
               synthesized function
    plmtable - a pointer to a double array of size `2*bw*(bw-m)`;
               contains the PRECOMPUTED plms, i.e. associated Legendre
               functions. E.g. Should be generated by a call to `PmlTableGen()`
               (see pml.c).

    Note that these Legendres are normalized with norm equal to 1!
*/
void Naive_SynthesizeX(double* coeffs, const int bw, const int m, double* result, double* plmtable) {
    memset(result, 0, sizeof(double) * 2 * bw); // make sure result is zeroed out

    for (int i = 0; i < bw - m; ++i) {
        double coeff = coeffs[i];

        if (coeff == 0.) {
            plmtable += (2 * bw);
            continue;
        }

        for (int j = 0; j < 2 * bw; ++j)
            result[j] += coeff * plmtable[j];

        plmtable += (2 * bw);
    }
}