#ifndef ZEYU_MATRIX_FUNCTIONS_H
#define ZEYU_MATRIX_FUNCTIONS_H

#define breakdown_free_flag (0)
#define check_step_by_step (0)
#define debug_BF_GMRES (0)
#define debug_MATRIX_ANALYSIS (0)
#define KRYLOV_NORM_CUTOFF (0.0)
#define disable_additional_basis (0)

using namespace std;

// calculate condition number
// if condition number>1/local_tol then a new basis vector is found.
// SUSSMAN add local_tol parameter (local_tol should be a small number
// e.g. 10^{-6})
double CondNum(double **H, const int m, const int n, 
	 const int sm, const int sn,double local_tol);

// calculate Householder vector
void House(const double *const x, double *v, double &beta, const int n);

// calculate Givens rotation matrix
//    _     _ T _  _     _ _
//   |  c  s | |  a |   | r |
//                    = 
//   |_-s  c_| |_ b_|   |_0_|
//    
void Givens(const double a, const double b, double &c, double &s);

// calculate Givens rotation matrix, overload
//    _     _ T _      _  _     _     _    _
//   |  c  s | |  a  b  ||  c  s |   | x  0 |
//                                 = 
//   |_-s  c_| |_ b  c _||_-s  c_|   |_0  y_|
//    
void Givens(const double aii, const double aij, const double ajj, double &c, double &s);

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

// QR factorization for least squares problems
// A: m x n, x: sn, b: sm
// To solve min||Ax-b|| for x:
// 1. A = QR;
// 2. d=QTb;
// 3. solve Rx=d.
void LeastSquaresQR(double **A, double *x, const double *b, const int m, const int n, const int sm, const int sn);

// calculate determinant of a square matrix using PLU decomposition
// det(A) = (-1)^S * det(U)
double GetDeterminant(double **A, const int m);

// PLU decomposition
// A: m x m, is overwritten by L and U, A <= (L-E)+U
// P: m+1, P[m] = S
void PLUDecomposition(double **A, int *P, const int m);

// Jacobi eigenvalue algorithm
// A: m x m, symmetric matrix
// D: m, store eigenvalues
void JacobiEigenvalue(double **A, double *D, const int m);

#endif
