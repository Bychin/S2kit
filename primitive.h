#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H

double L2_an(int, int);

double L2_cn(int, int);

double L2_cn_inv(int, int);

double L2_ancn(int, int);

void vec_add(double*, double*, double*, int);

void vec_mul(double, double*, double*, int);

void vec_dot(double*, double*, double*, int);

void AcosOfChebyshevNodes(int, double*);

void ChebyshevNodes(int, double*);

void Pmm_L2(int, double*, int, double*);

void P_eval(int, double*, double*, double*, double*, int);
#endif /* _PRIMITIVE_H */
