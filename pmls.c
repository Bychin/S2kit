/*
    Source code for generating cosine transforms of Pml and Gml functions.
*/

#include "pmls.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "primitive.h"

/*
    Generates all of the Pmls for a specified value of `m`.

    bw        - bandwidth;
    m         - order;
    storeplm  - array of size 2*bw*(bw-m);
    workspace - array of size 16*bw // TODO add size checkers?

    P(m,l,j) respresents the associated Legendre function P_l^m evaluated at the j-th Chebyshev point
    (for the bandwidth `bw`) cos((2 * j + 1) * PI / (2 * bw)).

    The array is placed in storeplm as follows:
    P(m,m,0)      P(m,m,1)    ...  P(m,m,2*bw-1)
    P(m,m+1,0)    P(m,m+1,1)  ...  P(m,m+1,2*bw-1)
    P(m,m+2,0)    P(m,m+2,1)  ...  P(m,m+2,2*bw-1)
    ...
    P(m,bw-1,0)   P(m,bw-1,1) ...  P(m,bw-1,2*bw-1)

    This array will eventually be used by the naive transform algorithm.
    This function will precompute the arrays necessary for the algorithm.
*/
void PmlTableGen(const int bw, const int m, double* storeplm, double* workspace) {
    int size = 2 * bw;

    double* prevprev = workspace;
    double* prev = prevprev + size;
    double* temp1 = prev + size;
    double* temp2 = temp1 + size;
    double* temp3 = temp2 + size;
    double* temp4 = temp3 + size;
    double* x_i = temp4 + size;
    double* eval_args = x_i + size;

    // get the evaluation nodes
    ChebyshevNodes(size, x_i);
    AcosOfChebyshevNodes(size, eval_args);

    // set initial values of first two Pmls
    for (int i = 0; i < size; ++i)
        prevprev[i] = 0.;

    if (!m)
        for (int i = 0; i < size; ++i)
            prev[i] = M_SQRT1_2;
    else
        Pmm_L2(m, eval_args, size, prev);

    memcpy(storeplm, prev, sizeof(double) * size);

    for (int i = 0; i < bw - m - 1; ++i) {
        vec_mul(L2_cn(m, m + i), prevprev, temp1, size);
        vec_dot(prev, x_i, temp2, size);
        vec_mul(L2_an(m, m + i), temp2, temp3, size);
        vec_add(temp3, temp1, temp4, size); // temp4 now contains P(m,m+i+1)

        storeplm += size;
        memcpy(storeplm, temp4, sizeof(double) * size);
        memcpy(prevprev, prev, sizeof(double) * size);
        memcpy(prev, temp4, sizeof(double) * size);
    }
}
