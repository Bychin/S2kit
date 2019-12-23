#ifndef _FST_SEMI_FLY_H
#define _FST_SEMI_FLY_H

#include <fftw3.h>

int seanindex(const int, const int, const int);

void FST_semi_fly(double*, double*, double*, double*, const int, double*, const int, const int,
                  fftw_plan*, fftw_plan*, double*);

void InvFST_semi_fly(double*, double*, double*, double*, const int, double*, const int, const int,
                     fftw_plan*, fftw_plan*);

void FZT_semi_fly(double*, double*, double*, double*, const int, double*, const int, fftw_plan*, double*);

void TransMult(double*, double*, double*, double*, double*, double*, const int);

void Conv2Sphere_semi_fly(double*, double*, double*, double*, double*, double*, const int, double*);

#endif // _FST_SEMI_FLY_H
