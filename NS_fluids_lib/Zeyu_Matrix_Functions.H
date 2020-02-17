#ifndef ZEYU_MATRIX_FUNCTIONS_H
#define ZEYU_MATRIX_FUNCTIONS_H

using namespace std;

// calculate condition number
double CondNum(double **H, const int m, const int n, const int sm, const int sn);

// calculate Householder vector
void House(const double *const x, double *v, double &beta, const int n);

// calculate Givens rotation matrix
//    _     _ T _  _     _ _
//   |  c  s | |  a |   | r |
//                    = 
//   |_-s  c_| |_ b_|   |_0_|
//    
void Givens(const double a, const double b, double &c, double &s);

// Householder Bidiagonalization
void GetBidiag(double **A, double *Bd, double *Bs, const int m, const int n);

// Golub-Kahan SVD step
void GKSVD(double *Bd, double *Bs, const int n);

// SVD algorithm
void SVD(double **A, double *D, const int m, const int n);

// zero 1st row of the bidiagonal matrix when the 1st diagonal entry is 0
void ZeroRow(double *Bd, double *Bs, const int n);

// zero the last column if the last diagonal entry is 0
void ZeroColumn(double *Bd, double *Bs, const int n);

#endif
