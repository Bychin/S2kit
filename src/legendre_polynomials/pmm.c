#include "pmm.h"

#include <math.h>

/*
    Returns L2-normed Pmm.

    The norming constant can be found in Sean's PhD thesis (I didn't find it).
    
    Note: input must be of order `m`, `eval_points` must be an array of length `n` of the angular
    arguments of evaluation points, `result` must be an array of length `n`
*/
void Pmm_L2(const int m, double* eval_points, const int n, double* result) {
    double normed_coeff = sqrt(m + 0.5);

    for (int i = 0; i < m; ++i)
        normed_coeff *= sqrt((m - (i / 2.)) / ((double)m - i));

    if (m)
        normed_coeff *= pow(2., -m / 2.);
    if (m % 2)
        normed_coeff *= -1.;

    for (int i = 0; i < n; ++i)
        result[i] = normed_coeff * pow(sin(eval_points[i]), m);
}