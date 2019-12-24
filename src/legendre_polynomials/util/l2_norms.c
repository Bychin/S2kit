/*
    Recurrence coefficents for L2-normed associated Legendre recurrence.

    When using these coeffs, make sure that inital Pmm function is also L2-normed.

    l - degree;
    m - order.
*/

#include "l2_norms.h"

#include <math.h>

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
