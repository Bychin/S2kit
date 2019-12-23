// TODO move this file into legendre_polynomials

#ifndef _UTIL_H
#define _UTIL_H

double L2_an(const int, const int);

double L2_cn(const int, const int);

double L2_cn_inv(const int, const int);

double L2_ancn(const int, const int);

void vec_add(double*, double*, double*, const int);

void vec_mul(const double, double*, double*, const int);

void vec_dot(double*, double*, double*, const int);

#endif // _UTIL_H
