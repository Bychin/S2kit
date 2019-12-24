/*
    Contains utility functions:

    1) IndexOfHarmonicCoeff(m,l,bw) - gives the position of the coefficient f-hat(m,l) in the one-row array;
    2) TransMult()                  - multiplies harmonic coefficients using Driscoll-Healy result,
                                      dual of convolution in "time" domain.
*/

#include "util.h"

#include <math.h>
#include <stdlib.h>

// Multiplies two complex numbers (x+iy * u+iv) and stores result separately for real and imaginary part.
void inline ComplexMult(const double x, const double y, const double u, const double v,
                        double* real_result, double* imag_result) {
    *real_result = x * u - y * v;
    *imag_result = x * v - y * u;
}

/*
    Returns the position of the coefficient `f-hat(m,l)` in the one-row
    array that Sean stores the spherical coefficients. This is needed
    to help preserve the symmetry that the coefficients have: 
    (`l` = degree, `m` = order, and `abs(m)` <= `l`)

    f-hat(l,-m) = (-1)^m * conjugate( f-hat(l,m) )

    Previous name was `seanindex`
*/
int IndexOfHarmonicCoeff(const int m, const int l, const int bw) {
    int bigL = bw - 1;

    if (m >= 0)
        return (m * (bigL + 1) - ((m * (m - 1)) / 2) + (l - m));

    return (((bigL * (bigL + 3)) / 2) + 1 + ((bigL + m) * (bigL + m + 1) / 2) + (l - abs(m)));
}

/*
    Multiplies harmonic coefficients of a function and a filter.
    See convolution theorem of Driscoll and Healy for details.

    bw - bandwidth of problem

    datacoeffs should be output of an SHT, filtercoeffs the
    output of an ZHT.  There should be (bw * bw) datacoeffs,
    and bw filtercoeffs.
    rres and ires should point to arrays of dimension bw * bw.
*/
void TransMult(double* rdatacoeffs, double* idatacoeffs, double* rfiltercoeffs, double* ifiltercoeffs, double* rres,
               double* ires, const int bw) {
    double* rdptr = rdatacoeffs;
    double* idptr = idatacoeffs;
    double* rrptr = rres;
    double* irptr = ires;

    for (int m = 0; m < bw; ++m) {
        for (int l = m; l < bw; ++l) {
            ComplexMult(rfiltercoeffs[l], ifiltercoeffs[l], rdptr[l - m], idptr[l - m], rrptr + l - m, irptr + l - m);

            rrptr[l - m] *= sqrt(4. * M_PI / (2. * l + 1.));
            irptr[l - m] *= sqrt(4. * M_PI / (2. * l + 1.));
        }
        rdptr += bw - m;
        idptr += bw - m;
        rrptr += bw - m;
        irptr += bw - m;
    }

    int size = 2 * bw;
    for (int m = bw + 1; m < size; ++m) {
        for (int l = size - m; l < bw; ++l) {
            ComplexMult(rfiltercoeffs[l], ifiltercoeffs[l], rdptr[l - size + m], idptr[l - size + m],
                        rrptr + l - size + m, irptr + l - size + m);

            rrptr[l - size + m] *= sqrt(4. * M_PI / (2. * l + 1.));
            irptr[l - size + m] *= sqrt(4. * M_PI / (2. * l + 1.));
        }

        rdptr += m - bw;
        idptr += m - bw;
        rrptr += m - bw;
        irptr += m - bw;
    }
}
