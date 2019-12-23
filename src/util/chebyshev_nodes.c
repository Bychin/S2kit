/*
    Utils for generating Chebyshev nodes
*/

#include "chebyshev_nodes.h"

#include <math.h>

/*
    Returns an array of the angular arguments of `n` Chebyshev nodes into `eval_points`.

    Note: `eval_points` must be an array of length `n`
*/
void AcosOfChebyshevNodes(const int n, double* eval_points) {
    double denominator = 2. * n;

    for (int i = 0; i < n; ++i)
        eval_points[i] = (2. * i + 1.) * M_PI / denominator;
}

/*
    Returns an array of `n` Chebyshev nodes into `eval_points`.

    Note: `eval_points` must be an array of length `n`
*/
void ChebyshevNodes(const int n, double* eval_points) {
    double denominator = 2. * n;

    for (int i = 0; i < n; ++i)
        eval_points[i] = cos((2. * i + 1.) * M_PI / denominator);
}
