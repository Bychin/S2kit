/*
    Contains the function that generates the weights for a bandwidth bw Legendre transform.
    Basically, it contains the implementation of the formula as defined in the tensor paper,
    and also given in the so(3) paper. It's just mentioned in the s^2 paper!

    This formula is slightly different from the one given in the original DH paper because
    they were sampling at the poles, and now we're not.

    In pseudo-TeX, the formula for the bandwidth B weights is

    w_B(j) = 2/B sin((pi*(2j+1))/(4B)) *
        sum_{k=0}^{B-1} 1/(2k+1)*sin((2j+1)(2k+1)pi/(4B))
    where j = 0, 1, ..., 2 B - 1

    Note that if you want to use these weights for an *odd* order transform, given the way
    the code is set up, you have to multiply the j-th weight by sin(pi*(2j+1)/(4B))
*/

#include "weights.h"

#include <math.h>

/*
    Makes weights for both even *and* odd order Legendre transforms for a given bandwidth bw.

    bw      - bandwidth of transform;
    weights - pointer to double array of length `4*bw`, this array will contain the even and odd weights,
        even weights start at `weights[0]`, and odd weights start at `weights[2*bw]`
*/
// TODO rename GenerateWeightsForDLT()
void makeweights(const int bw, double* weights) {
    double coeff = M_PI / (4. * bw);

    for (int i = 0; i < 2 * bw; ++i) {
        double sum = 0.;
        double k = 2. * i + 1.;

        for (int j = 0; j < bw; ++j)
            sum += 1. / (2. * j + 1.) * sin(k * (2. * j + 1.) * coeff);

        sum *= 2. * sin(k * coeff) / bw;

        weights[i] = sum;
        weights[i + 2 * bw] = sum * sin(k * coeff);
    }
}
