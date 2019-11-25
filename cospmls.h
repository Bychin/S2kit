#ifndef _COSPMLS_H
#define _COSPMLS_H

extern int TableSize(int, int);

extern int Spharmonic_TableSize(int);

extern int Reduced_SpharmonicTableSize(int, int);

extern int Reduced_Naive_TableSize(int, int);

extern int NewTableOffset(int, int);

extern void CosPmlTableGen(int, int, double *, double *);

extern int RowSize(int, int);

extern int Transpose_RowSize(int, int, int);

extern void Transpose_CosPmlTableGen(int, int, double *, double *);

extern double **Spharmonic_Pml_Table(int, double *, double *);

extern double **Transpose_Spharmonic_Pml_Table(double **, int, double *,
                                               double *);

extern double **SemiNaive_Naive_Pml_Table(int, int, double *, double *);

extern double **Transpose_SemiNaive_Naive_Pml_Table(double **, int, int,
                                                    double *, double *);

#endif
