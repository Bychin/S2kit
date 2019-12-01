/*
    source file which contains the function that generates the
    weights for a bandwidth bw legendre transform. Basically,
    it contains the implementation of the formula as defined in
    the tensor paper, and also given in the so(3) paper. It's
    just mentioned in the s^2 paper!

    This formula is slightly different from the one given in the
    original DH paper because they were sampling at the poles in
    that paper, and now we're not.

    In pseudo-TeX, the formula for the bandwidth B weights is

    w_B(j) = 2/B sin((pi*(2j+1))/(4B)) *
        sum_{k=0}^{B-1} 1/(2k+1)*sin((2j+1)(2k+1)pi/(4B))

    where j = 0, 1, ..., 2 B - 1


    Note that if you want to use these weights for an *odd*
    order transform, given the way the code is set up, you
    have to MULTIPLY the j-th weight by sin(pi*(2j+1)/(4B))
*/

#include <math.h>

/*
    makeweights: given a bandwidth bw, make weights for
    both even *and* odd order Legendre transforms.

    bw -> bandwidth of transform
    weights -> pointer to double array of length 4*bw;
        this array will contain the even and odd weights;
        even weights start at weight[0], and odd weights
        start at weights[2*bw]
*/

void makeweights(int bw, double* weights) {
    double fudge = M_PI / ((double)(4 * bw));

    for (int j = 0; j < 2 * bw; ++j) {
        double tmpsum = 0.0;

        for (int k = 0; k < bw; ++k)
            tmpsum += 1. / ((double)(2 * k + 1)) * sin((double)((2 * j + 1) * (2 * k + 1)) * fudge);

        tmpsum *= sin((double)(2 * j + 1) * fudge);
        tmpsum *= 2. / ((double)bw);

        weights[j] = tmpsum;
        weights[j + 2 * bw] = tmpsum * sin((double)(2 * j + 1) * fudge);
    }
}

// TODO: move to test file

// just a hack to test the above function
/*
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    int bw = atoi(argv[1]);

    double* weights = (double*)malloc(sizeof(double) * 4 * bw);

    makeweights(bw, weights);

    for (int i = 0; i < 2 * bw; i++)
        printf("%d\t%f\t%f\n", i, weights[i], weights[i + 2 * bw]);

    free(weights);

    return 0;
}
*/
