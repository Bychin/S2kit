#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H

double L2_an(const int, const int);

double L2_cn(const int, const int);

double L2_cn_inv(const int, const int);

double L2_ancn(const int, const int);

void vec_add(double*, double*, double*, const int);

void vec_mul(const double, double*, double*, const int);

void vec_dot(double*, double*, double*, const int);

void AcosOfChebyshevNodes(const int, double*);

void ChebyshevNodes(const int, double*);

void Pmm_L2(const int, double*, const int, double*);

void P_eval(const int, double*, double*, double*, double*, const int);

#endif // _PRIMITIVE_H
