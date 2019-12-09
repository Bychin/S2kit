#ifndef _FSTSEMI_FLY_H
#define _FSTSEMI_FLY_H

#include <fftw3.h>

#define compmult(a, b, c, d, e, f)                                                                                     \
    (e) = ((a) * (c)) - ((b) * (d));                                                                                   \
    (f) = ((a) * (d)) + ((b) * (c))

extern int seanindex(int, int, int);

extern void TransMult(double*, double*, double*, double*, double*, double*, int);

extern void FST_semi_fly(double*, double*, double*, double*, int, double*, int, int, fftw_plan*, fftw_plan*, double*);

extern void InvFST_semi_fly(double*, double*, double*, double*, int, double*, int, int, fftw_plan*, fftw_plan*);

extern void FZT_semi_fly(double*, double*, double*, double*, int, double*, int, fftw_plan*, double*);

extern void Conv2Sphere_semi_fly(double*, double*, double*, double*, double*, double*, int, double*);

#endif // _FSTSEMI_FLY_H
