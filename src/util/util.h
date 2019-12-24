#ifndef _UTIL_H
#define _UTIL_H

// DataFormat is used in FST functions to determine samples' data format.
typedef enum {
    COMPLEX = 0,
    REAL
} DataFormat; 

int IndexOfHarmonicCoeff(const int, const int, const int);

void TransMult(double*, double*, double*, double*, double*, double*, const int);

#endif // _UTIL_H
