/*
    Source code for generating cosine transforms of Pml and Gml functions.
*/

#include "cospml.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "util/chebyshev_nodes.h"

#include "pml.h"
#include "pmm.h"
#include "util/l2_norms.h"
#include "util/vector_funcs.h"

/*
    Utility functions for table management.
*/

/*
    Computes the number of non-zero entries in a table containing
    cosine series coefficients of all the P(m,l) or G(m,l) functions
    necessary for computing the seminaive transform for a given
    bandwidth `bw` and order `m`.

    Assumes the table is generated by `CosPmlTableGen()`.
*/
int TableSize(const int m, const int bw) {
    int fudge;
    int a1, a2;

    int k = bw / 2;

    if (bw % 2) { // if the bandwidth is odd
        fudge = (m + 1) % 2;

        a1 = k * (k + 1);
        a2 = fudge * (k + 1);

    } else { // bandwidth is even
        fudge = m % 2;

        a1 = (k - fudge) * (k - fudge + 1);
        a2 = fudge * k;
    }

    int fudge2 = m / 2;
    int a3 = fudge2 * (fudge2 + 1);

    return (a1 + a2 - a3);
}

/*
    Returns an integer value for the amount of space necessary
    to fill out an entire spharmonic table. Note that in the
    above `TableSize()` formula, you need to sum this formula
    over `m` as `m` ranges from 0 to `(bw-1)`.
    The critical closed form that you need is that:

    \sum_{k=0}^n = \frac{(n(n+1)(2n+1)}{6}

    You also need to account for integer division.
    From this you should derive an upper bound on the
    amount of space.

    Some notes - because of integer division, you need to account for
    a fudge factor - this is the additional " + bw" at the
    end.  This gaurantees that you will always have slightly more
    space than you need, which is clearly better than underestimating!
    Also, if bw > 512, the closed form
    fails because of the bw*bw*bw term (at least on Sun Sparcstations)
    so the loop computation is used instead.

    Also, the transpose is exactly the same size, obviously.
*/
// TODO check bw > 512 (~750 ok?)
int Spharmonic_TableSize(const int bw) {
    if (bw <= 512) {
        return (((4 * bw * bw * bw) + (6 * bw * bw) - (8 * bw)) / 24) + bw;
    }

    int sum = 0;
    for (int m = 0; m < bw; ++m)
        sum += TableSize(m, bw);

    return sum;
}

/* 
    This is a "reduced" version of `Spharmonic_TableSize(m)`.

    Returns an integer value for the amount of space necessary
    to fill out a spharmonic table
    if interesting in using it only for orders up to (but NOT
    including) order `m`.
    This will be used in the hybrid algorithm's call of the
    semi-naive algorithm (which won't need the full table ...
    hopefully this'll cut down on the memory usage).

    Also, the transpose is exactly the same size, obviously.
*/
int Reduced_SpharmonicTableSize(const int bw, const int m) {
    // TODO optimize? or just use Spharmonic_TableSize? (say no to economy)

    int sum = 0;
    for (int i = 0; i < m; ++i)
        sum += TableSize(i, bw);

    return sum;
}

/*
    Computes the location of the first coefficient of Pml for an array
    containing cosine series coefficients of Pml or Gml functions.

    Assumes the table is generated by `CosPmlTableGen()`.
*/
int TableOffset(const int m, const int l) {
    int tm = m;
    int tl = l;

    if (m % 2) {
        tm = m - 1;
        tl = l - 1;
    }

    int offset = ((tl / 2) * ((tl / 2) + 1)) - ((tm / 2) * ((tm / 2) + 1));
    if (tl % 2)
        offset += (tl / 2) + 1;

    return offset;
}

/*
    GenerateS all of the cosine series for L2-normalized Pmls or Gmls for
    a specified value of `m`. Note especially that since series are
    zero-striped, all zeroes have been removed.

    tablespace points to a double array of size TableSize(m,bw);

    Workspace needs to be `9 * bw`

    Let P(m,l,j) represent the j-th coefficient of the
    cosine series representation of Pml. The array
    stuffed into tablespace is organized as follows:

    P(m,m,0)    P(m,m,2)   ... P(m,m,m)
    P(m,m+1,1)  P(m,m+1,3) ... P(m,m+1,m+1)
    P(m,m+2,0)  P(m,m+2,2) ... P(m,m+2,m+2)

    etc.  Appropriate modifications are made for `m` odd (Gml functions).

    NOTE that the Pmls or Gmls are being sampled at bw-many points,
    and not 2*bw-many points. I can get away with this. HOWEVER, I
    need to multiply the coefficients by sqrt(2), because the expected
    input of the seminaive transform of bandwidth bw will be sampled
    at 2-bw many points. So the sqrt(2) is a scaling factor.
*/
void CosPmlTableGen(const int bw, const int m, double* tablespace, double* workspace) {
    double* prevprev = workspace;
    double* prev = prevprev + bw;
    double* temp1 = prev + bw;
    double* temp2 = temp1 + bw;
    double* temp3 = temp2 + bw;
    double* temp4 = temp3 + bw;
    double* x_i = temp4 + bw;
    double* eval_args = x_i + bw;
    double* cosres = eval_args + bw;

    double* tableptr = tablespace;

    fftw_plan p = fftw_plan_r2r_1d(bw, temp4, cosres, FFTW_REDFT10, FFTW_ESTIMATE);

    // set the initial number of evaluation points to appropriate amount

    // get the evaluation nodes
    ChebyshevNodes(bw, x_i);
    AcosOfChebyshevNodes(bw, eval_args);

    // set initial values of first two Pmls
    for (int i = 0; i < bw; ++i)
        prevprev[i] = 0.0;

    if (m == 0)
        for (int i = 0; i < bw; ++i)
            prev[i] = M_SQRT1_2; // sqrt(1/2)
    else
        Pmm_L2(m, eval_args, bw, prev);

    if (m % 2)
        for (int i = 0; i < bw; ++i)
            prev[i] /= sin(eval_args[i]);

    
    int k; // set k to highest degree coefficient
    if ((m % 2) == 0)
        k = m;
    else
        k = m - 1;

    // compute cosine transform
    memcpy(temp4, prev, sizeof(double) * bw);
    fftw_execute(p);
    cosres[0] *= M_SQRT1_2;
    double fudge = 1. / sqrt(bw);
    for (int i = 0; i < bw; ++i)
        cosres[i] *= fudge;

    // store what we've got so far
    for (int i = 0; i <= k; i += 2)
        tableptr[i / 2] = cosres[i];

    // update tableptr
    tableptr += k / 2 + 1;

    // generate remaining Pmls
    for (int i = 0; i < bw - m - 1; ++i) {
        vec_mul(L2_cn(m, m + i), prevprev, temp1, bw);
        vec_dot(prev, x_i, temp2, bw);
        vec_mul(L2_an(m, m + i), temp2, temp3, bw);
        vec_add(temp3, temp1, temp4, bw); // temp4 now contains P(m,m+i+1)

        // compute cosine transform
        fftw_execute(p);
        cosres[0] *= M_SQRT1_2;
        for (int j = 0; j < bw; ++j)
            cosres[j] *= fudge;

        // update degree counter
        ++k;

        for (int j = (i % 2) ? 0 : 1; j <= k; j += 2)
            tableptr[j / 2] = cosres[j];

        // update tableptr
        tableptr += k / 2 + 1;

        // update Pi and P(i+1)
        memcpy(prevprev, prev, sizeof(double) * bw);
        memcpy(prev, temp4, sizeof(double) * bw);
    }

    fftw_destroy_plan(p);
}

/*
    RowSize returns the number of non-zero coefficients in a row of the
    cospmltable if were really in matrix form.  Helpful in transpose
    computations. It is helpful to think of the parameter l as
    the row of the corresponding matrix.
*/
int RowSize(const int m, const int l) {
    if (l < m)
        return 0;

    if (m % 2 == 0)
        return (l / 2) + 1;

    return ((l - 1) / 2) + 1;
}

/* 
    Transposed row size returns the number of non-zero coefficients
    in the transposition of the matrix representing a cospmltable.
    Used for generating arrays for inverse seminaive transform.
    Unlike `RowSize()`, need to know the bandwidth `bw`. Also, in
    the cospml array, the first `m+1` rows are empty, but in
    the transpose, all rows have non-zero entries, and the first
    `m+1` columns are empty. So the input parameters are a bit different
    in the you need to specify the row you want.
*/
int Transpose_RowSize(const int row, const int m, const int bw) {
    if (row >= bw)
        return 0;

    if ((m % 2) == 0) {
        if (row <= m)
            return (bw - m) / 2;

        return  ((bw - row - 1) / 2) + 1;
    }

    if (row == (bw - 1))
        return 0;

    if (row >= m)
        return Transpose_RowSize(row + 1, m - 1, bw);
    else
        return Transpose_RowSize(row + 1, m - 1, bw) - (row % 2);
}

/* 
    Inverse transform is transposition of forward transform.
    Thus, need to provide transposed version of table
    returned by `CosPmlTableGen()`. This function does that
    by taking as input a `cos_pml_table` for a particular value
    of `bw` and `m`, and loads the `result` as a transposed,
    decimated version of it for use by an inverse seminaive
    transform computation.

    `result` needs to be of size `TableSize(m, bw)`
*/

void Transpose_CosPmlTableGen(const int bw, const int m, double* cos_pml_table, double* result) {
    // Recall that `cos_pml_table` has had all the zeroes stripped out,
    // and that if `m` is odd, then it is really a Gml function, which affects indexing a bit.

    // note that the number of non-zero entries is the same as in the non-transposed case
    double* trans_tableptr = result;
    double* tableptr;

    // traverse the `cos_pml_table`, loading appropriate values into the rows of transposed array
    if (m == bw - 1) {
        memcpy(result, cos_pml_table, sizeof(double) * TableSize(m, bw));
        return;
    }

    for (int row = 0; row < bw; ++row) {
        // if `m` odd, no need to do last row - all zeroes
        if (row == (bw - 1) && (m % 2))
            return;

        // get the rowsize for the transposed array
        int rowsize = Transpose_RowSize(row, m, bw);

        // compute the starting point for values in `cos_pml_table`
        if (row <= m) {
            if ((row % 2) == 0)
                tableptr = cos_pml_table + (row / 2);
            else
                tableptr = cos_pml_table + (m / 2) + 1 + (row / 2);
        } else {
            // then the highest degree coefficient of P(m,row) should be the first coefficient loaded
            // into the transposed array, so figure out where this point is
            int offset = 0;
            int end_row = (m % 2) == 0 ? row : row + 1;
            for (int i = m; i <= end_row; ++i)
                offset += RowSize(m, i);
            // we are pointing one element too far, so decrement
            --offset;

            tableptr = cos_pml_table + offset;
        }

        // `stride` is how far we need to jump between values in `cos_pml_table`, i.e.,
        // to traverse the columns of the `cos_pml_table`. Need to set initial value.
        // `stride` always increases by 2 after that.
        int stride;
        if (row <= m)
            stride = m + 2 - (m % 2) + (row % 2);
        else
            stride = row + 2;

        // load up this row of the transposed table
        int costable_offset = 0;
        for (int i = 0; i < rowsize; ++i) {
            trans_tableptr[i] = tableptr[costable_offset];
            costable_offset += stride;
            stride += 2;
        }

        trans_tableptr += rowsize;
    }
}

/*
    Returns all of the (cosine transforms of) Pmls and Gmls necessary
    to do a full spherical harmonic transform, i.e., it calls
    `CosPmlTableGen()` for each value of `m` less than `bw`, returning a
    table of tables (a pointer of type (double**), which points
    to an array of size `m`, each containing a (double*) pointer
    to a set of CosPml or CosGml values, which are the (decimated)
    cosine series representations of Pml (even `m`) or Gml (odd `m`)
    functions. See `CosPmlTableGen()` for further clarification.

    bw - bandwidth of the problem;
    resultspace - need to allocate `Spharmonic_TableSize(bw)` for storing results
    workspace - needs to be `16*bw`

    Note that `resultspace` is necessary and contains the results/values
    so one should be careful about when it is OK to re-use this space.
    workspace, though, does not have any meaning after this function is
    finished executing.
*/
// TODO: 2 fors in 1?
// check by s2_conv_memo
double** Spharmonic_Pml_Table(const int bw, double* resultspace, double* workspace) {
    double** spharmonic_pml_table = (double**)malloc(sizeof(double*) * bw);

    // traverse the array, assigning a location in the `resultspace` to each pointer
    spharmonic_pml_table[0] = resultspace;
    for (int i = 1; i < bw; ++i)
        spharmonic_pml_table[i] = spharmonic_pml_table[i - 1] + TableSize(i - 1, bw);

    // load up the array with CosPml and CosGml values
    for (int i = 0; i < bw; ++i)
        CosPmlTableGen(bw, i, spharmonic_pml_table[i], workspace);

    return spharmonic_pml_table;
}

/* 
    For the inverse semi-naive spharmonic transform, the "transpose"
    of the `spharmonic_pml_table` is needed. Need to be careful because the
    entries in the qspharmonic_pml_tableq have been decimated, i.e.,
    the zeroes have been stripped out.

    spharmonic_pml_table - generated by `Spharmonic_Pml_Table()`;
    bw - bandwidth of the problem;
    resultspace - need to allocate Spharmonic_TableSize(bw) for storing results.

    Allocates memory for the (double**) `resultspace`.
*/
double** Transpose_Spharmonic_Pml_Table(double** spharmonic_pml_table, const int bw, double* resultspace) {
    double** transpose_spharmonic_pml_table = (double**)malloc(sizeof(double*) * bw);

    // load up the `transpose_spharmonic_pml_table` by transposing the tables in `spharmonic_pml_table`
    transpose_spharmonic_pml_table[0] = resultspace;
    for (int i = 0; i < bw; ++i) {
        Transpose_CosPmlTableGen(bw, i, spharmonic_pml_table[i], transpose_spharmonic_pml_table[i]);

        if (i != (bw - 1))
            transpose_spharmonic_pml_table[i + 1] = transpose_spharmonic_pml_table[i] + TableSize(i, bw);
    }

    return transpose_spharmonic_pml_table;
}

/*
    Returns an integer value for the amount of space necessary to fill out
    a reduced naive table of Pmls if interested in using it only for orders `m` through `bw-1`.
*/
int Reduced_Naive_TableSize(const int bw, const int m) {
    int sum = 0;
    for (int i = m; i < bw; ++i)
        sum += 2 * bw * (bw - i);

    return sum;
}

/*
    Just like Spharmonic_Pml_Table(), except generates a table for use
    with the semi-naive and naive algorithms.

    bw - bandwidth of the problem;
    m - the cutoff order, where to switch from semi-naive to naive algorithms;
    resultspace - stores results, must be of size
        `Reduced_Naive_TableSize(bw, m) + Reduced_SpharmonicTableSize(bw, m)`;
*/
double** SemiNaive_Naive_Pml_Table(const int bw, const int m, double* resultspace, double* workspace) {
    double** seminaive_naive_table = (double**)malloc(sizeof(double) * (bw + 1));

    seminaive_naive_table[0] = resultspace;
    for (int i = 1; i < m; ++i)
        seminaive_naive_table[i] = seminaive_naive_table[i - 1] + TableSize(i - 1, bw);

    if (m != 0) {
        seminaive_naive_table[m] = seminaive_naive_table[m - 1] + TableSize(m - 1, bw);
    }

    for (int i = m + 1; i < bw; ++i)
        seminaive_naive_table[i] = seminaive_naive_table[i - 1] + (2 * bw * (bw - (i - 1)));

    // load up the array with CosPml and CosGml values
    for (int i = 0; i < m; ++i)
        CosPmlTableGen(bw, i, seminaive_naive_table[i], workspace);

    // load up Pml values
    for (int i = m; i < bw; ++i)
        PmlTableGen(bw, i, seminaive_naive_table[i], workspace);

    return seminaive_naive_table;
}

/*
    For the inverse seminaive_naive transform, need the "transpose"
    of the seminaive_naive_pml_table. Need to be careful because the
    entries in the seminaive portion have been decimated, i.e.,
    the zeroes have been stripped out.

    seminaive_naive_pml_table - generated by `SemiNaive_Naive_Pml_Table()`;
    bw - bandwidth of the problem;
    m - the cutoff order, where to switch from semi-naive to naive algorithms;
    resultspace - need to allocate
        Reduced_Naive_TableSize(bw, m) + Reduced_SpharmonicTableSize(bw, m) for storing results;
    workspace - size of `16*bw`
*/
double** Transpose_SemiNaive_Naive_Pml_Table(double** seminaive_naive_pml_table, const int bw, const int m,
                                             double* resultspace, double* workspace) {
    double**  trans_seminaive_naive_pml_table = (double**)malloc(sizeof(double*) * (bw + 1));

    // need to load up the `transpose_seminaive_naive_pml_table` by transposing
    // the tables in the seminiave portion of `seminaive_naive_pml_table`

    trans_seminaive_naive_pml_table[0] = resultspace;
    for (int i = 1; i < m; ++i)
        trans_seminaive_naive_pml_table[i] = trans_seminaive_naive_pml_table[i - 1] + TableSize(i - 1, bw);

    if (m != 0) {
        int lastspace = TableSize(m - 1, bw);
        trans_seminaive_naive_pml_table[m] = trans_seminaive_naive_pml_table[m - 1] + lastspace;
    }

    for (int i = m + 1; i < bw; ++i)
        trans_seminaive_naive_pml_table[i] = trans_seminaive_naive_pml_table[i - 1] + (2 * bw * (bw - (i - 1)));

    for (int i = 0; i < m; ++i) {
        Transpose_CosPmlTableGen(bw, i, seminaive_naive_pml_table[i], trans_seminaive_naive_pml_table[i]);

        if (i != (bw - 1)) {
            trans_seminaive_naive_pml_table[i + 1] = trans_seminaive_naive_pml_table[i] + TableSize(i, bw);
        }
    }

    // load up Pml values
    for (int i = m; i < bw; ++i)
        PmlTableGen(bw, i, trans_seminaive_naive_pml_table[i], workspace);

    return trans_seminaive_naive_pml_table;
}
