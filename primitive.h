#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H

extern double L2_an(int, int);

extern double L2_cn(int, int);

extern double L2_cn_inv(int, int);

extern double L2_ancn(int, int);

extern void vec_add(double *, double *, double *, int);

extern void vec_mul(double, double *, double *, int);

extern void vec_pt_mul(double *, double *, double *, int);

extern void ArcCosEvalPts(int, double *);

extern void EvalPts(int, double *);

extern void Pmm_L2(int, double *, int, double *);

extern void P_eval(int, double *, double *, double *, double *, int);
#endif /* _PRIMITIVE_H */
