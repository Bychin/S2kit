/*
    Some "primitive" functions that are used in cospmls.c
*/

// TODO move to utils?

#include <math.h>
#include <string.h>

/*
    Recurrence coefficents for L2-normed associated Legendre recurrence.

    When using these coeffs, make sure that inital Pmm function is also L2-normed.

    l - degree;
    m - order.
*/

double L2_an(const int m, const int l) {
    double m_d = (double)m;
    double l_d = (double)l;

    return (
        sqrt(
            ((2. * l_d + 3.) / (2. * l_d + 1.)) * ((l_d - m_d + 1.) / (l_d + m_d + 1.))
        ) * ((2. * l_d + 1.) / (l_d - m_d + 1.))
    );
}

// Note: if input `l` is zero, return 0
double L2_cn(const int m, const int l) {
    if (l == 0) {
        return 0.;
    }

    double m_d = (double)m;
    double l_d = (double)l;

    return (
        -1.0 * sqrt(
            ((2. * l_d + 3.) / (2. * l_d - 1.)) * ((l_d - m_d + 1.) / (l_d + m_d + 1.)) *
                ((l_d - m_d) / (l_d + m_d))
        ) * ((l_d + m_d) / (l_d - m_d + 1.))
    );
}

// Note: when using the reverse recurrence, instead of `1/L2_cn(m,l)`
double L2_cn_inv(const int m, const int l) {
    double m_d = (double)m;
    double l_d = (double)l;

    return (
        -(1.0 + (1. - 2. * m_d) / (m_d + l_d)) *
            sqrt(
                ((-1. + 2. * l_d) / (3. + 2. * l_d)) * (
                    (l_d + l_d * l_d + m_d + 2. * l_d * m_d + m_d * m_d) /
                        (l_d + l_d * l_d - m_d - 2. * l_d * m_d + m_d * m_d)
                )
            )
    );
}

// Note: when using the reverse recurrence, instead of `-L2_an(m,l)/L2_cn(m,l)`
// TODO: remove?
double L2_ancn(const int m, const int l) {
    double m_d = (double)m;
    double l_d = (double)l;

    return (
        sqrt(4. + ((4. * m_d * m_d - 1.) / (l_d * l_d - m_d * m_d)))
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
    double m_d = (double)m;
    double norming_const = sqrt(m_d + 0.5);

    for (int i = 0; i < m; ++i)
        norming_const *= sqrt((m_d - (i / 2.)) / (m_d - i));

    if (m != 0)
        norming_const *= pow(2., -m_d / 2.);
    if (m % 2 != 0)
        norming_const *= -1.;

    for (int i = 0; i < n; ++i)
        result[i] = norming_const * pow(sin(eval_points[i]), m_d);
}

/*
    Synthesizes a function which is the weighted sum of L2-normed associated Legendre functions.

    Note: `coeffs` array must contain `bw-m` coefficients ordered from zeroth degree to `bw-1`,
    `eval_points` must be an array of length `2*bw` of the angular (arccos) arguments of the desired
    evaluation points, `result` must be an array of length `2*bw`, `workspace` must be an array of 
    length `16*bw`
*/
// TODO was this code somewhere else? duplicate?
// TODO remove? this func is not used anywhere
// TODO rename EvaluateP
void P_eval(const int m, double* coeffs, double* eval_points, double* result, double* workspace, const int bw) {
    int n = 2 * bw;

    double* prevprev = workspace; // TODO rename
    double* prev = prevprev + n;
    double* temp1 = prev + n;
    double* temp2 = temp1 + n;
    double* temp3 = temp2 + n;
    double* temp4 = temp3 + n;
    double* x_i = temp4 + n;

    // get the evaluation nodes
    ChebyshevNodes(n, x_i);

    for (int i = 0; i < n; ++i)
        prevprev[i] = 0.0;

    if (m == 0) {
        for (int i = 0; i < n; ++i) {
            prev[i] = M_SQRT1_2;

            // mult by first coeff and add to result
            result[i] = coeffs[0] * prev[i];
        }
    } else {
        Pmm_L2(m, eval_points, n, prev);
        double splat = coeffs[0];
        for (int i = 0; i < n; ++i)
            result[i] = splat * prev[i];
    }

    for (int i = 0; i < bw - m - 1; ++i) {
        vec_mul(L2_cn(m, m + i), prevprev, temp1, n);
        vec_dot(prev, x_i, temp2, n);
        vec_mul(L2_an(m, m + i), temp2, temp3, n);
        vec_add(temp3, temp1, temp4, n); // temp4 now contains P(m,m+i+1)

        // add weighted P(m,m+i+1) to the result
        double splat = coeffs[i + 1];
        for (int j = 0; j < n; ++j)
            result[j] += splat * temp4[j];

        memcpy(prevprev, prev, sizeof(double) * n);
        memcpy(prev, temp4, sizeof(double) * n);
    }
}
