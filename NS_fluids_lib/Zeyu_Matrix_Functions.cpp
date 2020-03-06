#include <iostream>
#include <cmath>

#define STANDALONE 0

#if (STANDALONE==1)
// do nothing
#else
#include <AMReX_Box.H>

#endif

using namespace std;

#include "Zeyu_Matrix_Functions.H"

// calculate condition number
// SUSSMAN add local_tol parameter (local_tol should be a small number
// e.g. 10^{-6})
double CondNum(double **H, const int m, const int n, 
	const int sm, const int sn,double local_tol)
{
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::CondNum(): start" << endl;
    }

    if(sm > m || sn > n || sm<=0 || sn <= 0){
        cout << "Sub-matrix is larger than the orginial matrix!" << endl;
        abort();
    }

    if ((local_tol>0.0)&&(local_tol<1.0)) {
     // do nothing
    } else {
     std::cout << "local_tol invalid\n";
     abort();
    }

    double **M = new double *[sm];
    for(int i = 0; i < sm; ++i)
        M[i] = new double [sn]; // SUSSMAN

    for(int i = 0; i < sm; ++i)
        for(int j = 0; j < sn; ++j)
            M[i][j] = H[i][j];

    double D_QR[sn];
    double D_JAC[sn];

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "calling SVD sm=" << sm << '\n';
     std::cout << "calling SVD sn=" << sn << '\n';
    }

      // modifies M
    SVD(M, D_QR, sm, sn);

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "done with SVD\n";
    }

    double **MTM = new double *[sn];
    for(int i = 0; i < sn; ++i){
     MTM[i] = new double [sn]; //SUSSMAN
    }

    for(int i = 0; i < sn; ++i){
     for(int j = 0; j < sn; ++j){
      MTM[i][j] = 0.0;
      for(int k = 0; k < sm; ++k)
       MTM[i][j] += H[k][i] * H[k][j];
     }
    }

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "calling JacobiEigenvalue\n";
    }

    JacobiEigenvalue(MTM, D_JAC, sn);

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "done with JacobiEigenvalue\n";
    }

    for(int i = 0; i < sn; ++i) {
     if (D_JAC[i]>=0.0) {
      D_JAC[i] = sqrt(D_JAC[i]);
     } else {
      cout << "D_JAC<0" << endl;
      abort();
     }
    }
    for(int i = 0; i < sn; ++i) {
     if (debug_MATRIX_ANALYSIS==1) {
      std::cout << "deleting MTM[i] i=" << i << "\n";
     }
     delete[] MTM[i]; //SUSSMAN
    }
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "deleting MTM\n";
    }
    delete[] MTM;
    
    for(int i = 0; i < sm; ++i)
       delete[] M[i]; // SUSSMAN
    delete[] M;

    double min_JAC = 1.e20;
    double min_QR = 1.e20;
    double max_JAC = 0.0;
    double max_QR = 0.0;
    for(int i = 0; i < sn; ++i){
     if ((D_QR[i]<=0.0)||(D_QR[i]>=0.0)) {
      // do nothing
     } else {
      cout << "D_QR invalid" << endl;
      abort();
     }
     D_QR[i] = std::abs(D_QR[i]);
     if(D_QR[i] > max_QR) max_QR = D_QR[i];
     if(D_QR[i] < min_QR) min_QR = D_QR[i];
     if(D_JAC[i] > max_JAC) max_JAC = D_JAC[i];
     if(D_JAC[i] < min_JAC) min_JAC = D_JAC[i];
    }

/*
    if (1==0) {
     for(int i = 0; i < sn-1; ++i){
        int imax = i;
        for(int j = i+1; j < sn; ++j){
            if(D[j] > D[imax])
                imax = j;
        }
        if(imax != i){
            double temp = D[i];
            D[i] = D[imax];
            D[imax] = temp;
        }
        //cout << "\n" << D[i];
     }
     //cout << "\n" << D[sn-1];
     //cout << endl;

     double Dn = 0.0;
     for(int i = 0; i < sn; ++i)
        Dn += D[i] * D[i];
     Dn = sqrt(Dn);
    
     //sanity check
     double **B = new double *[sn];
     for(int i = 0; i < sn; ++i){
        B[i] = new double [sn]; //SUSSMAN
        for(int j = 0; j < sn; ++j){
            B[i][j] = 0.0;
            for(int k = 0; k < sm; ++k)
                B[i][j] += H[k][i] * H[k][j];
        }
     }
     if(0){
        for(int it = 0; it < sn; ++it){
            for(int i = 0; i < sn; ++i)
                B[i][i] -= D[it] * D[it];

            double detB = GetDeterminant(B, sn);
            for(int i = 0; i < sn; ++i)
                B[i][i] += D[it] * D[it];
            if(std::abs(detB) <= 1.e-16){
              //cout << D[it] << " is singluar value of Matrix H." << endl;
            }
            else{
                cout << "Something is wrong for singluar value "<< D[it] 
                     <<"! det(ATA - sigma^2*I) = " << std::abs(detB) << endl;
            }
        }
     }
     else{
        double Dt[sn] {};
        JacobiEigenvalue(B, Dt, sn);
        for(int i = 0; i < sn-1; ++i){
            int imax = i;
            for(int j = i+1; j < sn; ++j){
                if(Dt[j] > Dt[imax])
                    imax = j;
            }
            if(imax != i){
                double temp = Dt[i];
                Dt[i] = Dt[imax];
                Dt[imax] = temp;
            }
        }
        //cout << endl;
        //for(int i = 0; i < sn; ++i)
        //    cout << sqrt(Dt[i]) << endl;
        for(int i = 0; i < sn; ++i){
            Dt[i] = sqrt(Dt[i]);
            double re = (D[i] - Dt[i])/D[i];
            if(std::abs(re) <= 1.e-12){
                //do nothing!
            }
            else{
                cout << "Something is wrong in singular value" 
                     << "[" << i << "], " 
                     << "relative residual is " << std::abs(re) << endl;
            }
        }
     }

     for(int i = 0; i < sn; ++i)
        delete[] B[i]; //SUSSMAN
     delete[] B;
     //end sanity check
    }
*/
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::CondNum(): end" << endl;
    }
    if (min_QR<0.0) {
     cout << "min_QR<0 invalid." << endl;
     abort();
    } else if (min_QR==0.0) {
     std::cout << "Warning! (QR) Zero singular value." << endl;
    } else if (min_QR>0.0) {
     // do nothing
    } else {
     cout << "min_QR invalid NaN." << endl;
     abort();
    }

    if (min_JAC<0.0) {
     cout << "min_JAC<0 invalid." << endl;
     abort();
    } else if (min_JAC==0.0) {
     cout << "Warning! (JAC) Zero singular value." << endl;
    } else if (min_JAC>0.0) {
     // do nothing
    } else {
     cout << "min_JAC invalid NaN." << endl;
     abort();
    }


    double local_QR_condnum=0.0;
    double local_JAC_condnum=0.0;
    double local_condnum=0.0;
    if ((max_QR>0.0)&&(max_JAC>0.0)) {

     if (min_QR/max_QR<local_tol) {
      local_QR_condnum=1.0+1.0/local_tol;
     } else if (min_QR/max_QR>=local_tol) {
      local_QR_condnum=max_QR/min_QR;
     } else {
      cout << "min_QR, max_QR invalid" << endl;
      abort();
     }
     if (min_JAC/max_JAC<local_tol) {
      local_JAC_condnum=1.0+1.0/local_tol;
     } else if (min_JAC/max_JAC>=local_tol) {
      local_JAC_condnum=max_JAC/min_JAC;
     } else {
      cout << "min_JAC, max_JAC invalid" << endl;
      abort();
     }

     double rel_error=std::abs(2.0*(max_QR-max_JAC)/
		     (local_JAC_condnum+local_QR_condnum));
     double rel_error_min=std::abs(2.0*(min_QR-min_JAC)/
		     (local_JAC_condnum+local_QR_condnum));
     if ((rel_error<=1.0e-5)&&(rel_error_min<=1.0e-5)) {
      local_condnum=0.5*(local_QR_condnum+local_JAC_condnum);
     } else {
      cout << "rel_error too big" << endl;
      std::cout << "local_QR_condnum=" << local_QR_condnum << endl;
      std::cout << "local_JAC_condnum=" << local_JAC_condnum << endl;
      std::cout << "max_QR=" << max_QR << endl;
      std::cout << "max_JAC=" << max_JAC << endl;
      std::cout << "min_QR=" << min_QR << endl;
      std::cout << "min_JAC=" << min_JAC << endl;
      std::cout << "rel_error=" << rel_error << endl;
      std::cout << "rel_error_min=" << rel_error_min << endl;
      std::cout << "m,n,sm,sn = " << m << ' ' << n << ' ' << sm << ' ' <<
       sn << endl;
      for (int i=0;i<sm;i++) {
       for (int j=0;j<sn;j++) {
        std::cout << "i,j,H[i][j] " << i << ' ' << j << ' ' << 
          H[i][j] << endl;
       }
      }
      
#if (STANDALONE==1)
      abort();
#else
      amrex::Error("cannot find condnum");
#endif
     }
    } else {
      cout << "max_QR or max_JAC invalid" << endl;
      abort();
    }
    return local_condnum;
}

// calculate Householder vector
void House(const double *const x, double *v, double &beta, const int n)
{
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::House(): start" << endl;
    }

    double sigma = 0.0;
    double xx[n-1];
    v[0] = 1.0;

    for(int i = 0; i < n-1; ++i){
        xx[i] = x[i+1];
        sigma += xx[i] * xx[i];
        v[i+1] = xx[i];
    }

    if(sigma == 0.0 && x[0] >= 0.0)
        beta = 0.0;
    else if(sigma == 0.0 && x[0] < 0.0)
        beta = -2.0;
    else{
        double mu = sqrt(x[0] * x[0] + sigma);
        if(x[0] <= 0.0)
            v[0] = x[0] - mu;
        else
            v[0] = - sigma / (x[0] + mu);
        beta = 2. * v[0] * v[0] / (sigma + v[0] * v[0]);
        double v0 = v[0];
        for(int i = 0; i < n; ++i)
            v[i] = v[i] / v0;
    }

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::House(): end" << endl;
    }
}

// calculate Givens rotation matrix
//    _     _ T _  _     _ _
//   |  c  s | |  a |   | r |
//                    = 
//   |_-s  c_| |_ b_|   |_0_|
//
void Givens(const double a, const double b, double &c, double &s)
{
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::Givens(): start" << endl;
    }

    if(b == 0.0){
        c = 1.0;
        s = 0.0;
    }
    else{
        if(std::abs(b) > std::abs(a)){
            s = 1.0 / sqrt(1. + a / b * a / b);
            c = s * (- a / b);
        }
        else{
            c = 1.0 / sqrt(1. + b / a * b / a);
            s = c * (- b / a);
        }
    }

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::Givens(): end" << endl;
    }
}

// calculate Givens rotation matrix, overload
//    _     _ T _      _  _     _     _    _
//   |  c  s | |  a  b  ||  c  s |   | x  0 |
//                                 = 
//   |_-s  c_| |_ b  c _||_-s  c_|   |_0  y_|
//    
void Givens(const double aii, const double aij, const double ajj, double &c, double &s)
{
    if(aij == 0.0){
        c = 1.0;
        s = 0.0;
    }
    else{
        double thet = 0.0;
        if(std::abs(ajj-aii) <= 1.e-16){
            if(aij > 0)
                thet = 0.5 * asin(1.0);
            if(aij < 0)
                thet = -0.5 * asin(1.0);
        }
        else{
            thet = 0.5 * atan(2.0 * aij / (ajj - aii));
        }
        c = cos(thet);
        s = sin(thet);
    }
}

// Householder Bidiagonalization
// modifies A
void  GetBidiag(double **A, double *Bd, double *Bs, const int m, const int n)
{
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::GetBidiag(): start" << endl;
    }

    int it;
    for(it = 0; it < n; ++it){
        double* x=new double[m-it];  //SUSSMAN
        for(int i = it; i < m; ++i)
            x[i-it] = A[i][it];
        double* v1=new double[m-it]; //SUSSMAN
        double beta1;
        House(x, v1, beta1, m-it);
        double **u = new double *[m-it];
        for(int i = 0; i < m-it; ++i)
            u[i] = new double [m-it]; //SUSSMAN
        for(int i = 0; i < m-it; ++i){
            for(int j = 0; j < m-it; ++j){
                if(i == j)
                    u[i][j] = 1.0;
                else
                    u[i][j] = 0.0;
                u[i][j] = u[i][j] - beta1 * v1[i] * v1[j];
            }
        }
        delete [] x; //SUSSMAN
        delete [] v1; //SUSSMAN

        for(int j = 0; j < n-it; ++j){
            double* temp=new double[m-it];  //SUSSMAN
            for(int i = 0; i < m-it; ++i)
             temp[i]=0.0;
            for(int i = 0; i < m-it; ++i)
                for(int k = 0; k < m-it; ++k)
                    temp[i] += u[i][k] * A[k+it][j+it];
            for(int i = 0; i < m-it; ++i)
                A[i+it][j+it] = temp[i];
            delete[] temp; //SUSSMAN
        }

        for(int i = 0; i < m-it; ++i)
            delete[] u[i]; //SUSSMAN
        delete[] u;

        Bd[it] = A[it][it];

        if(it < n-2){
            double* v2=new double[n-it-1];
            double beta2;
            double* xr=new double[n-it-1]; //SUSSMAN
            for(int i = 0; i < n-it-1; ++i)
                xr[i] = A[it][i+it+1];
            House(xr, v2, beta2, n-it-1);
            double **v = new double *[n-it-1];
            for(int i = 0; i < n-it-1; ++i)
                v[i] = new double [n-it-1]; //SUSSMAN
            for(int i = 0; i < n-it-1; ++i){
                for(int j = 0; j < n-it-1; ++j){
                    if(i == j)
                        v[i][j] = 1.0;
                    else
                        v[i][j] = 0.0;
                    v[i][j] = v[i][j] - beta2 * v2[i] * v2[j];
                }
            }

            for(int i = 0; i < m-it; ++i){
                double* temp=new double[n-it-1];  //SUSSMAN
                for(int j = 0; j < n-it-1; ++j)
                 temp[j]=0.0;
                for(int j = 0; j < n-it-1; ++j)
                    for(int k = 0; k < n-it-1; ++k)
                        temp[j] += A[i+it][k+it+1] * v[k][j];
                for(int j = 0; j < n-it-1; ++j)
                    A[i+it][j+it+1] = temp[j];
                delete[] temp; //SUSSMAN
            }

            for(int i = 0; i < n-it-1; ++i)
                delete[] v[i];
            delete[] v;
            delete[] v2; //SUSSMAN
            delete[] xr; //SUSSMAN

            Bs[it] = A[it][it+1]; 
        }
        if(it == n-2)
            Bs[it] = A[it][it+1];
    }
 
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::GetBidiag(): end" << endl;
    }
}

// Golub-Kahan SVD step
void GKSVD(double *Bd, double *Bs, const int n)
{
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::GKSVD(): start" << endl;
    }

    double t1;//t(n-1, n-1)
    if(n == 2){
        t1 = Bd[n-2] * Bd[n-2];
    }
    else{
        t1 = Bs[n-3] * Bs[n-3] + Bd[n-2] * Bd[n-2];
    }
    double t2 = Bs[n-2] * Bs[n-2] + Bd[n-1] * Bd[n-1];//t(n, n)
    double t3 = Bs[n-2] * Bd[n-2];//t(n, n-1)
    double d = (t1 - t2) / 2.0;
    double mu = t2 - t3 * t3 / (d + d / std::abs(d) * sqrt(d * d + t3 * t3));

    double y = Bd[0] * Bd[0] - mu;
    double z = Bd[0] * Bs[0];
    double c, s;
    double a1 = 0.0;
    double a2 = 0.0;
    for(int k = 0; k < n - 1; ++k){
        Givens(y, z, c, s);
        double b1[3][2];
        if(k == 0){
            b1[0][0] = 0.0;
            b1[0][1] = 0.0;
        }
        else{
            b1[0][0] = Bs[k-1];
            b1[0][1] = a2;
        }
        b1[1][0] = Bd[k];
        b1[1][1] = Bs[k];
        b1[2][0] = 0.0;
        b1[2][1] = Bd[k+1];

        for(int i = 0; i < 3; ++i){
            double temp1 = b1[i][0] * c + b1[i][1] * (-s);
            double temp2 = b1[i][0] * s + b1[i][1] * c;
            b1[i][0] = temp1;
            b1[i][1] = temp2;
        }

        a1 = b1[2][0];
        if(k != 0) Bs[k-1] = b1[0][0];
        Bd[k] = b1[1][0];
        Bs[k] = b1[1][1];
        Bd[k+1] = b1[2][1];
        y = Bd[k];
        z = a1;
        Givens(y, z, c, s);
        double b2[2][3];
        b2[0][0] = Bd[k];
        b2[0][1] = Bs[k];
        b2[0][2] = 0.0;
        b2[1][0] = a1;
        b2[1][1] = Bd[k+1];
        if(k == n - 2){
            b2[1][2] = 0.0;
        }
        else{
            b2[1][2] = Bs[k+1];
        }

        for(int j = 0; j < 3; ++j){
            double temp1 = c * b2[0][j] + (-s) * b2[1][j];
            double temp2 = s * b2[0][j] + c * b2[1][j];
            b2[0][j] = temp1;
            b2[1][j] = temp2;
        }

        if(k != n - 2) a2 = b2[0][2];
        Bd[k] = b2[0][0];
        Bs[k] = b2[0][1];
        Bd[k+1] = b2[1][1];
        if(k != n - 2) Bs[k+1] = b2[1][2];
        if(k < n - 2){
            y = Bs[k];
            z = a2;
        }
    }

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::GKSVD(): end" << endl;
    }
}

// SVD algorithm
// modifies A
void SVD(double **A, double *D, const int m, const int n)
{
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "in SVD m=" << m << '\n';
     std::cout << "in SVD n=" << n << '\n';
    }
    double **ASAVE = new double *[m];
    for(int i = 0; i < m; ++i)
     ASAVE[i] = new double [n]; // SUSSMAN

    for(int i = 0; i < m; ++i)
     for(int j = 0; j < n; ++j)
      ASAVE[i][j]=A[i][j];

    double* Bd=new double[n];  //SUSSMAN
    for(int i = 0; i < n; ++i)
     Bd[i]=0.0;
    int nn;
    if(n == 1)
        nn = 2;
    else
        nn = n;
    double* Bs=new double[nn-1];  //SUSSMAN
    for(int i = 0; i < nn-1; ++i)
     Bs[i]=0.0;

     // modifies A
    GetBidiag(A, Bd, Bs, m, n);

    if(n == 1) D[0] = Bd[0];

    double Bn = 0.0;
    for(int i = 0; i < n-1; ++i)
        Bn = Bn + Bd[i] * Bd[i] + Bs[i] * Bs[i];
    Bn = Bn + Bd[n-1] * Bd[n-1];
    Bn = sqrt(Bn);
    int p = 0;
    int q = 0;
    int pp = p;
    int qq = q;
    double tol = 1.e-16;
    while(q < n && n > 1){
        if (debug_MATRIX_ANALYSIS==1) {
 	 std::cout << "in q<n loop, q= " << q << " n= " << n << '\n';
 	 std::cout << "in q<n loop, m= " << m << " n= " << n << '\n';
	 for (int ideb=0;ideb<m;ideb++) {
	  for (int jdeb=0;jdeb<n;jdeb++) {
	   std::cout << "i,j,Aij " << ideb << ' ' << jdeb << ' ' <<
	     A[ideb][jdeb] << '\n';
	   std::cout << "i,j,ASAVEij " << ideb << ' ' << jdeb << ' ' <<
	     ASAVE[ideb][jdeb] << '\n';
	  }
	 }
	}

        //set b(i,i+1) to zero, if |b(i,i+1)| <= tol*(|b(i,i)|+|b(i+1,i+1)|)
        for(int i = 0; i < n-p-q-1; ++i){//relative index
            if(std::abs(Bs[p+i]) <= tol * (std::abs(Bd[p+i]) + std::abs(Bd[p+i+1])))
                Bs[p+i] = 0.0;//absolute index
        }
        for(int i = 0; i < n-p-q-1; ++i){
            if(Bs[n-q-2-i] == 0.0){
                ++qq;
                D[n-1-q-i] = Bd[n-q-1-i];
                if(i == n-p-q-2){
                    ++qq;
                    D[n-2-q-i] = Bd[p];
                    qq = qq + pp;
                    pp = 0;
                }
            }
            else
                break;
        }
        q = qq;
        p = pp;
        if(q < n){
            for(int i = 0; i < n-p-q-1; ++i){
                if(Bs[p+i] == 0.0){
                    ++pp;
                    D[p+i] = Bd[p+i];
                }
                else
                    break;
            }
            p = pp;
        }
        else
            break;
        
        int bi = 0;
        for(int i = 0; i < n-p-q; ++i){
            //if any diagonal entry is zero, then zero the superdiagnoal entry in the same row and same column
            if(std::abs(Bd[p+i]) <= tol * Bn){
                Bd[p+i] = 0.0;
                if(i != n-p-q-1 && Bs[p+i] != 0.0)
                    ZeroRow(Bd+p+i, Bs+p+i, n-p-q-i);
                if(i == bi){//if Bd[p+bi] = 0.0
                    bi = i + 1;
                    continue;
                }
                if(Bs[p+i-1] != 0.0)
                    ZeroColumn(Bd+p+bi, Bs+p+bi, i-bi+1);
                if(i-bi == 1){
                    bi = i + 1;
                    continue;
                }
                else{
                    GKSVD(Bd+p+bi, Bs+p+bi, i-bi);
                    bi = i + 1;
                }
            }
        }
        if(bi != n-p-q){
            if(bi == n-p-q-1)
                continue;
            GKSVD(Bd+p+bi, Bs+p+bi, n-p-q-bi);
        }
    }
    delete[] Bd; //SUSSMAN
    delete[] Bs; //SUSSMAN

    for(int i = 0; i < m; ++i)
       delete[] ASAVE[i]; // SUSSMAN
    delete[] ASAVE;

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "after SVD m=" << m << '\n';
     std::cout << "after SVD n=" << n << '\n';
    }
}

// zero 1st row of the bidiagonal matrix when the 1st diagonal entry is 0
void ZeroRow(double *Bd, double *Bs, const int n)
{
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::ZeroRow(): start" << endl;
    }

    double a = Bs[0];
    for(int it = 1; it < n; ++it){
        double y = a;
        double z = Bd[it];
        double c, s;
        Givens(y, z, c, s);
        swap(c, s);
        c = -c;
        double b[2][2];
        b[0][0]=y;
        b[1][0]=z;
        b[0][1]=0.0;
        b[1][1]=0.0; 
        if(it != n-1) b[1][1] = Bs[it];

        for(int j = 0; j < 2; ++j){
            double temp1 = c * b[0][j] + (-s) * b[1][j];
            double temp2 = s * b[0][j] + c * b[1][j];
            b[0][j] = temp1;
            b[1][j] = temp2;
        }

        if(it == 1) Bs[0] = 0.0;
        Bd[it] = b[1][0];
        if(it != n-1) Bs[it] = b[1][1];
        a = b[0][1];
    }

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::ZeroRow(): end" << endl;
    }
}

// zero the last column if the last diagonal entry is 0
void ZeroColumn(double *Bd, double *Bs, const int n)
{
    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::ZeroColumn(): start" << endl;
    }

    double a = Bs[n-2];
    for(int it = n-2; it >= 0; --it){
        double y = Bd[it];
        double z = a;
        double c, s;
        Givens(y, z, c, s);
        double b[2][2];
        b[0][0]=0.0;
        b[1][0]=y;
        b[0][1]=0.0;
        b[1][1]=z; 
        if(it != 0) b[0][0] = Bs[it-1];

        for(int i = 0; i < 2; ++i){
            double temp1 = b[i][0] * c + b[i][1] * (-s);
            double temp2 = b[i][0] * s + b[i][1] * c;
            b[i][0] = temp1;
            b[i][1] = temp2;
        }

        if(it == n-2) Bs[it] = 0.0;
        Bd[it] = b[1][0];
        if(it != 0) Bs[it-1] = b[0][0];
        a = b[0][1];
    }

    if (debug_MATRIX_ANALYSIS==1) {
     std::cout << "Matrix::ZeroColumn(): end" << endl;
    }
}

// QR factorization for least squares problems
// A: m x n, x: sn, b: sm
// To solve min||Ax-b|| for x:
// 1. A = QR;
// 2. d=QTb;
// 3. solve Rx=d.
void LeastSquaresQR(double **A, double *x, const double *b, const int m, const int n, const int sm, const int sn)
{
    //define R(mx(n+1)), d is included in the last column of R
    double **R = new double *[sm];
    for(int i = 0; i < sm; ++i){
        R[i] = new double [sn+1]; //SUSSMAN
    }
    for(int i = 0; i < sm; ++i){
        for(int j = 0; j < sn+1; ++j){
            if(j == sn)
                R[i][j] = b[i];
            else{
                R[i][j] = A[i][j];
            }
        }
    }

    for(int it = 0; it < sn; ++it){
        double y = R[it][it];
        for(int i = it+1; i < sm; ++i){
            double z = R[i][it];
            if(std::abs(z) < 1.0e-16) 
                continue;
            else{
                double c, s;
                Givens(y, z, c, s);

                for(int j = it; j < sn+1; ++j){
                    double temp1 = c * R[it][j] + (-s) * R[i][j];
                    double temp2 = s * R[it][j] + c * R[i][j];
                    R[it][j] = temp1;
                    R[i][j] = temp2;
                }
            }
        }
    }

    if(std::abs(R[sn-1][sn-1]) <= 1.e-16){
        cout << "Warning! R[" << sn-1 << "][" << sn-1 
             << "] is close to zero in QR fatorization!" << endl;
        R[sn-1][sn-1] = 1.e-16;
    }
    x[sn-1] = R[sn-1][sn] / R[sn-1][sn-1];
    for(int i = 0; i < sn-1; ++i){
        x[sn-2-i] = R[sn-2-i][sn];
        for(int j = sn-1-i; j < sn; ++j)
            x[sn-2-i] -= R[sn-2-i][j] * x[j];
        if(std::abs(R[sn-2-i][sn-2-i]) <= 1.e-16){
            cout << "Warn! R[" << sn-2-i << "][" << sn-2-i 
                 << "] is close to zero in QR fatorization!" << endl;
            R[sn-2-i][sn-2-i] = 1.e-16;
        }
        x[sn-2-i] = x[sn-2-i] / R[sn-2-i][sn-2-i];
    }

    for(int i = 0; i < sm; ++i)
        delete[] R[i]; //SUSSMAN
    delete[] R;

    //sanity check
    double **ATA = new double *[sn];
    for(int i = 0; i < sn; ++i){
        ATA[i] = new double [sn]; //SUSSMAN
    }
    for(int i = 0; i < sn; ++i){
        for(int j = 0; j < sn; ++j){
            ATA[i][j] = 0.0;
            for(int k = 0; k < sm; ++k)
                ATA[i][j] += A[k][i] * A[k][j];
        }
    }
    double residual_verify = 0.0;
    double* ATAx=new double[sn];
    double* ATb=new double[sn];
    for(int i = 0; i < sn; ++i){
        ATAx[i] = 0.0;
        ATb[i] = 0.0;
        for(int j = 0; j < sn; ++j){
            ATAx[i] += ATA[i][j] * x[j];
            ATb[i] += A[j][i] * b[j];
        }
        for(int j = sn; j < sm; ++j)
            ATb[i] += A[j][i] * b[j];
        double res_comp = ATAx[i] - ATb[i];;
        residual_verify += res_comp * res_comp;
    }
    delete[] ATAx;
    delete[] ATb;

    residual_verify = sqrt(residual_verify)/sn;
    if(residual_verify <= 1.e-8){
        //cout << "x is the optimization solution." << endl;
    }
    else {
     cout << "Something wrong! ||ATAx - ATb|| = " << residual_verify << endl;
     abort();
    }
    for(int i = 0; i < sn; ++i)
        delete[] ATA[i]; //SUSSMAN
    delete[] ATA;
    //end sanity check
}

// calculate determinant of a square matrixi using PLU decomposition
// det(A) = (-1)^S * det(U)
double GetDeterminant(double **A, const int m)
{
    double detA;

    double **B = new double *[m];
    for(int i = 0; i < m; ++i){
        B[i] = new double [m]; //SUSSMAN
    }
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < m; ++j)
            B[i][j] = A[i][j];
    }

    int P[m+1];

    PLUDecomposition(B, P, m);

    double Bmax = 0.0;
    for(int i = 0; i < m; ++i){
        if(std::abs(B[i][i]) > Bmax)
            Bmax = std::abs(B[i][i]);
    }
    detA = pow(-1, P[m]);
    for(int i = 0; i < m; ++i){
        //cout << B[i][i] << endl;
        if(std::abs(B[m-1-i][m-1-i]) <= 1.e-9 * Bmax){
            //B[m-1-i][m-1-i] = 0.0;
            detA = 0.0;
            break;
        }
        else
            detA *= B[m-1-i][m-1-i];
    }

    //sanity check
    /*double Pt[m][m];
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < m; ++j){
            Pt[j][i] = 0;
            if(j == P[i])
                Pt[j][i] = 1.0;
        }
        //cout << P[i] << " ";
    }
    double Lt[m][m];
    double Ut[m][m];
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < m; ++j){
            if(i == j)
                Lt[i][j] = 1.0;
            else if(i < j)
                Lt[i][j] = 0.0;
            else
                Lt[i][j] = B[i][j];
            if(i > j)
                Ut[i][j] = 0.0;
            else
                Ut[i][j] = B[i][j];
        }
    }
    double At[m][m];
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < m; ++j){
            At[i][j] = 0.0;
            for(int k = 0; k < m; ++k)
                At[i][j] += Pt[i][k] * Lt[k][j];
        }
    }
    for(int i = 0; i < m; ++i){
        double temp[m];
        for(int j = 0; j < m; ++j){
            temp[j] = 0.0;
            for(int k = 0; k < m; ++k)
                temp[j] += At[i][k] * Ut[k][j];
        }
        for(int j = 0; j < m; ++j)
            At[i][j] = temp[j];
    }
    double An = 0.0;
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < m; ++j)
            An += A[i][j] * A[i][j];
    An = sqrt(An);
    double r = 0;
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < m; ++j){
            At[i][j] = At[i][j] - A[i][j];
            if(std::abs(At[i][j]) <= 1.e-14 * An)
                At[i][j] = 0.0;
            cout << Ut[i][j] << " ";
            r += At[i][j] * At[i][j];
        }
        cout << endl;
    }
    r = sqrt(r);
    if(r <= 1.e-16){
        //Do nothing!
    }
    else{
        cout << "Something is wrong in LPU decomposition! ||PLU-A|| = " 
             << r << endl;
    }*/
    //end sanity check

    for(int i = 0; i < m; ++i)
        delete[] B[i]; //SUSSMAN
    delete[] B;

    return detA;
}

// PLU decomposition
// A: m x m, is overwritten by L and U, A <= (L-E)+U
// P: m+1, P[m] = S
void PLUDecomposition(double **A, int *P, const int m)
{
    for(int i = 0; i < m; ++i)
        P[i] = i;
    P[m] = 0;
    
    double An = 0.0;
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < m; ++j)
            An += A[i][j] * A[i][j];
    An = sqrt(An);

    for(int it = 0; it < m; ++it){
        double max = 0.0;
        int imax = it;
        for(int i = it; i < m; ++i){
            if(std::abs(A[i][it]) > max){
                max = std::abs(A[i][it]);
                imax = i;
            }
        }

        if(max <= 1.e-16 * An){
            if(it == m-1){
                //A[it][it] = 0.0;
            }
            else{
                cerr << "Error!PLU decomposition failure!" << endl;
                abort();
            }
        }

        if(imax != it){
            int temp = P[it];
            P[it] = P[imax];
            P[imax] = temp;

            double *pt = A[it];
            A[it] = A[imax];
            A[imax] = pt;

            P[m]++;
        }

        for(int i = it+1; i < m; ++i){
            A[i][it] /= A[it][it];
            for(int j = it+1; j < m; ++j){
                A[i][j] -= A[i][it] * A[it][j];
            }
        }
    }

}

//SUSSMAN
void abort_and_dump(double** A,int m,int n) {

 std::cout << "A is mxn where m,n = " << m << ' ' << n << endl;
 for (int i=0;i<m;i++) {
  for (int j=0;j<n;j++) {
    std::cout << "i,j,A[i][j] " << i << ' ' << j << ' ' <<
	 A[i][j] << endl;
  }
 }
#if (STANDALONE==1)
 abort();
#else
 amrex::Error("abort_and_dump");
#endif
			 
} // abort_and_dump

// Jacobi eigenvalue algorithm
// A: m x m, symmetric matrix
// D: m, store eigenvalues
void JacobiEigenvalue(double **A, double *D, const int m)
{
    double **M = new double *[m];
    for(int i = 0; i < m; ++i){
        M[i] = new double [m]; //SUSSMAN
    }
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < m; ++j)
            M[i][j] = A[i][j];
    }

    int* imax=new int[m]; //SUSSMAM

    int max_i = 0;
    int max_j = 0;
    for(int i = 0; i < m-1; ++i){
        imax[i] = i+1;
        for(int j = i+2; j < m; ++j){
            if(std::abs(M[i][j]) > std::abs(M[i][imax[i]]))
                imax[i] = j;
        }
        if(i == 0){
            max_i = i;
            max_j = imax[i];
        }
        else{
            if(std::abs(M[i][imax[i]]) > std::abs(M[max_i][max_j])){
                max_i = i;
                max_j = imax[i];
            }
        }
    }

    while(1){
        double Mn = 0.0;
        for(int i = 0; i < m; ++i)
            Mn += M[i][i];
        Mn = sqrt(Mn);

        double Dn = 0.0;
        for(int i = 0; i < m; ++i){
            D[i] = M[i][i];
            Dn += D[i] * D[i];
        }
        Dn = sqrt(Dn);

        if(std::abs(M[max_i][max_j]) <= 1.e-16 * Dn || m == 1)
            break;

        double aii = M[max_i][max_i];
        double ajj = M[max_j][max_j];
        double aij = M[max_i][max_j];
        double c, s;
        Givens(aii, aij, ajj, c, s);

        for(int k = 0; k < m; ++k){
            if(k == max_i || k == max_j)
                continue;
            double a = M[max_i][k];
            double b = M[max_j][k];
            M[max_i][k] = c * a - s * b;
            M[k][max_i] = M[max_i][k];
            M[max_j][k] = s * a + c * b;
            M[k][max_j] = M[max_j][k];
        }
        M[max_i][max_i] = c * c * aii - 2. * s * c * aij + s * s * ajj;
        M[max_j][max_j] = s * s * aii + 2. * s * c * aij + c * c * ajj;
        M[max_i][max_j] = (c * c - s * s) * aij + s * c * (aii - ajj);
        M[max_j][max_i] = M[max_i][max_j];

        int mj = max_j;
        if(mj == m-1) mj = mj -1;
        for(int i = 0; i <= mj; ++i){
            imax[i] = i+1;
            for(int j = i+2; j < m; ++j){

	        if ((i>=0)&&(i<m)) {
	         if ((j>=0)&&(j<m)) {
		  if ((imax[i]>=0)&&(imax[i]<m)) {
		   // do nothing    	   
		  } else {
		   std::cout << "imax[i] invalid\n";
		   std::cout << "i,j= " << i << ' ' << j << endl;
		   std::cout << "imax[i]= " << imax[i] << endl;
                   abort_and_dump(A,m,m); //SUSSMAN
		  }
		 } else {
		  std::cout << "j invalid\n";
		  std::cout << "i,j= " << i << ' ' << j << endl;
                  abort_and_dump(A,m,m); //SUSSMAN
		 }
		} else {
		 std::cout << "i invalid\n";
		 std::cout << "i,j= " << i << ' ' << j << endl;
                 abort_and_dump(A,m,m); //SUSSMAN
		}


                if(std::abs(M[i][j]) > std::abs(M[i][imax[i]]))
                    imax[i] = j;
            }
            if(i == 0){
                max_i = i;
                max_j = imax[i];
            }
            else{
	        if ((i>=0)&&(i<m)) {
		 if ((max_i>=0)&&(max_i<m)) {
		  if ((max_j>=0)&&(max_j<m)) {
		   if ((imax[i]>=0)&&(imax[i]<m)) {
		    // do nothing    	   
		   } else {
		    std::cout << "imax[i] invalid\n";
		    std::cout << "i,max_i,max_j= " << i << ' ' <<
		     max_i << ' ' << max_j << endl;
		    std::cout << "imax[i]= " << imax[i] << endl;
                    abort_and_dump(A,m,m); //SUSSMAN
		   }
		  } else {
		   std::cout << "max_j invalid\n";
		   std::cout << "i,max_i,max_j= " << i << ' ' <<
		     max_i << ' ' << max_j << endl;
                   abort_and_dump(A,m,m); //SUSSMAN
		  }
		 } else {
		  std::cout << "max_i invalid\n";
		  std::cout << "i,max_i,max_j= " << i << ' ' <<
		    max_i << ' ' << max_j << endl;
                  abort_and_dump(A,m,m); //SUSSMAN
		 }
		} else {
		 std::cout << "i invalid\n";
		 std::cout << "i,max_i,max_j= " << i << ' ' <<
		    max_i << ' ' << max_j << endl;
                 abort_and_dump(A,m,m); //SUSSMAN
		}

                if(std::abs(M[i][imax[i]]) > std::abs(M[max_i][max_j])){
                    max_i = i;
                    max_j = imax[i];
                }
            }
        }
    }

    delete[] imax; //SUSSMAN

    for(int i = 0; i < m; ++i)
        delete[] M[i];
    delete[] M;
}

#undef STANDALONE
