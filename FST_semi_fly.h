#ifndef _FST_SEMI_FLY_H
#define _FST_SEMI_FLY_H

#include <fftw3.h>

extern int seanindex(int, int, int);

extern void TransMult(double*, double*, double*, double*, double*, double*, int);

extern void FST_semi_fly(double*, double*, double*, double*, int, double*, int, int, fftw_plan*, fftw_plan*, double*);

extern void InvFST_semi_fly(double*, double*, double*, double*, int, double*, int, int, fftw_plan*, fftw_plan*);

extern void FZT_semi_fly(double*, double*, double*, double*, int, double*, int, fftw_plan*, double*);

extern void Conv2Sphere_semi_fly(double*, double*, double*, double*, double*, double*, int, double*);

#endif // _FST_SEMI_FLY_H
