/*
    Some "primitive" functions that are used in cospmls.c
*/

// TODO move to utils?

#include "primitive.h"

#include <math.h>
#include <string.h>

/*
    Recurrence coefficents for L2-normed associated Legendre recurrence.

    When using these coeffs, make sure that inital Pmm function is also L2-normed.

    l - degree;
    m - order.
*/

double L2_an(const int m, const int l) {
    return (
        sqrt(
            ((2. * l + 3.) / (2. * l + 1.)) * ((l - m + 1.) / (l + m + 1.))
        ) * ((2. * l + 1.) / (l - m + 1.))
    );
}

double L2_cn(const int m, const int l) {
    if (!l) {
        return 0.;
    }

    return (
        -1.0 * sqrt(
            ((2. * l + 3.) / (2. * l - 1.)) * ((l - m + 1.) / (l + m + 1.)) * (((double)l - m) / ((double)l + m))
        ) * ((l + m) / (l - m + 1.))
    );
}

// Note: when using the reverse recurrence, instead of `1/L2_cn(m,l)`
double L2_cn_inv(const int m, const int l) {
    return (
        -(1.0 + (1. - 2. * m) / ((double)m + l)) *
            sqrt(
                ((-1. + 2. * l) / (3. + 2. * l)) * (
                    (l + l * l + m + 2. * l * m + m * m) /
                        (l + l * l - m - 2. * l * m + m * m)
                )
            )
    );
}

// Note: when using the reverse recurrence, instead of `-L2_an(m,l)/L2_cn(m,l)`
double L2_ancn(const int m, const int l) {
    return (
        sqrt(4. + ((4. * m * m - 1.) / ((double)l * l - m * m)))
    );
}

/*
    Vector arithmetic operations
*/

/*
    Adds two vectors into a third one.
    `result = v1 + v2`

    Note: `result` and `v{1,2}` must be vectors of length `len`
*/
void vec_add(double* v1, double* v2, double* result, const int len) {
    for (int i = 0; i < len; ++i)
        result[i] = v1[i] + v2[i];
}

/*
    Multiplies the vector `v` by `scalar` and returns in `result`.

    Note: `result` and `v` must be vectors of length `len`
*/
void vec_mul(const double scalar, double* v, double* result, const int len) {
    for (int i = 0; i < len; ++i)
        result[i] = scalar * v[i];
}

/* 
    Returns dot product of `v1` and `v2` in `result`

    Note: `result` and `v{1,2}` must be vectors of length `len`
*/
void vec_dot(double* data1, double* data2, double* result, const int len) {
    for (int i = 0; i < len; ++i)
        result[i] = data1[i] * data2[i];
}

/*
    Utils for Chebyshev nodes
*/

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

    for (int i = 0; i < n; i++)
        eval_points[i] = cos((2. * i + 1.) * M_PI / denominator);
}

// TODO: move somewhere?

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
