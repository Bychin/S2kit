#ifndef _COSPMLS_H
#define _COSPMLS_H

// TODO extern?
extern int TableSize(const int, const int);

extern int Spharmonic_TableSize(const int);

extern int Reduced_SpharmonicTableSize(const int, const int);

extern int Reduced_Naive_TableSize(const int, const int);

extern int TableOffset(const int, const int);

extern void CosPmlTableGen(const int, const int, double*, double*);

extern int RowSize(const int, const int);

extern int Transpose_RowSize(const int, const int, const int);

extern void Transpose_CosPmlTableGen(const int, const int, double*, double*);

extern double** Spharmonic_Pml_Table(const int, double*, double*);

extern double** Transpose_Spharmonic_Pml_Table(double**, const int, double*);

extern double** SemiNaive_Naive_Pml_Table(const int, const int, double*, double*);

extern double** Transpose_SemiNaive_Naive_Pml_Table(double**, const int, const int, double*, double*);

#endif // _COSPMLS_H
