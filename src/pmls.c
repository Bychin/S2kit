/*
    Source code for generating cosine transforms of Pml and Gml functions.
*/

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
    (for the bandwidth bw) cos((2 * j + 1) * PI / (2 * bw)).

    The array is placed in storeplm as follows:
    P(m,m,0)      P(m,m,1)    ...  P(m,m,2*bw-1)
    P(m,m+1,0)    P(m,m+1,1)  ...  P(m,m+1,2*bw-1)
    P(m,m+2,0)    P(m,m+2,1)  ...  P(m,m+2,2*bw-1)
    ...
    P(m,bw-1,0)   P(m,bw-1,1) ...  P(m,bw-1,2*bw-1)

    This array will eventually be used by the naive transform algorithm.
    This function will precompute the arrays necessary for the algorithm.
*/
void PmlTableGen(int bw, int m, double* storeplm, double* workspace) {
    double *prev, *prevprev;
    double *temp1, *temp2, *temp3, *temp4;
    double *x_i, *eval_args;
    int i;

    prevprev = workspace;
    prev = prevprev + (2 * bw);
    temp1 = prev + (2 * bw);
    temp2 = temp1 + (2 * bw);
    temp3 = temp2 + (2 * bw);
    temp4 = temp3 + (2 * bw);
    x_i = temp4 + (2 * bw);
    eval_args = x_i + (2 * bw);

    /* get the evaluation nodes */
    ChebyshevNodes(2 * bw, x_i);
    AcosOfChebyshevNodes(2 * bw, eval_args);

    /* set initial values of first two Pmls */
    for (i = 0; i < 2 * bw; i++)
        prevprev[i] = 0.0;
    if (m == 0)
        for (i = 0; i < 2 * bw; i++)
            prev[i] = M_SQRT1_2;
    else
        Pmm_L2(m, eval_args, 2 * bw, prev);

    memcpy(storeplm, prev, sizeof(double) * 2 * bw);

    for (i = 0; i < bw - m - 1; i++) {
        vec_mul(L2_cn(m, m + i), prevprev, temp1, 2 * bw);
        vec_dot(prev, x_i, temp2, 2 * bw);
        vec_mul(L2_an(m, m + i), temp2, temp3, 2 * bw);
        vec_add(temp3, temp1, temp4, 2 * bw); /* temp4 now contains P(m,m+i+1) */

        storeplm += (2 * bw);
        memcpy(storeplm, temp4, sizeof(double) * 2 * bw);
        memcpy(prevprev, prev, sizeof(double) * 2 * bw);
        memcpy(prev, temp4, sizeof(double) * 2 * bw);
    }
}
