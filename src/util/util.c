/*
    Contains some utility functions.
*/

#include "util.h"

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
