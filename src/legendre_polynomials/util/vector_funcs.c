/*
    Vector arithmetic operations
*/

#include "vector_funcs.h"

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
