#include <iostream>
#include <math.h>

using namespace std;

#include "Zeyu_Matrix_Functions.H"

// calculate condition number
double CondNum(double **H, const int m, const int n, const int sm, const int sn)
{
    //if(iprint == 1)
    //    cout << "Matrix::CondNum(): start" << endl;

    if(sm > m || sn > n){
        cout << "Sub-matrix is larger than the orginial matrix!" << endl;
        abort();
    }

    double **M = new double *[sm];
    for(int i = 0; i < sm; ++i)
        *(M+i) = new double [sn];

    for(int i = 0; i < sm; ++i)
        for(int j = 0; j < sn; ++j)
            M[i][j] = H[i][j];

    double D[sn];
    SVD(M, D, sm, sn);
    for(int i = 0; i < sm; ++i)
        delete *(M+i);
    delete[] M;

    double min = 1.e20;
    double max = 0.0;
    for(int i = 0; i < sn; ++i){
        if(fabs(D[i]) > max) max = fabs(D[i]);
        if(fabs(D[i]) < min) min = fabs(D[i]);
    }

    //if(iprint == 1)
    //    cout << "Matrix::CondNum(): end" << endl;
    double Dn = 0.0;
    for(int i = 0; i < sn; ++i)
        Dn = Dn + D[i] * D[i];
    Dn = sqrt(Dn);
    
    if(fabs(min - 0.0) <= 1.e-16 * Dn){
     min=max*(1.0e-100);
    }
    return max / min;
}

// calculate Householder vector
void House(const double *const x, double *v, double &beta, const int n)
{
    //if(iprint == 1)
    //    cout << "Matrix::House(): start" << endl;

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

    //if(iprint == 1)
    //    cout << "Matrix::House(): end" << endl;
}

// calculate Givens rotation matrix
//    _     _ T _  _     _ _
//   |  c  s | |  a |2   | r |
//                    = 
//   |_-s  c_| |_ b_|   |_0_|
//
void Givens(const double a, const double b, double &c, double &s)
{
    //if(iprint == 1)
    //    cout << "Matrix::Givens(): start" << endl;

    if(b == 0.0){
        c = 1.0;
        s = 0.0;
    }
    else{
        if(fabs(b) > fabs(a)){
            s = 1.0 / sqrt(1. + a / b * a / b);
            c = s * (- a / b);
        }
        else{
            c = 1.0 / sqrt(1. + b / a * b / a);
            s = c * (- b / a);
        }
    }

    //if(iprint == 1)
    //    cout << "Matrix::Givens(): end" << endl;
}

// Householder Bidiagonalization
void  GetBidiag(double **A, double *Bd, double *Bs, const int m, const int n)
{
    //if(iprint == 1)
    //    cout << "Matrix::GetBidiag(): start" << endl;

    int it;
    for(it = 0; it < n; ++it){
        double x[m-it];
        for(int i = it; i < m; ++i)
            x[i-it] = A[i][it];
        double v1[m-it];
        double beta1;
        House(x, v1, beta1, m-it);
        double u[m-it][m-it];
        for(int i = 0; i < m-it; ++i){
            for(int j = 0; j < m-it; ++j){
                if(i == j)
                    u[i][j] = 1.0;
                else
                    u[i][j] = 0.0;
                u[i][j] = u[i][j] - beta1 * v1[i] * v1[j];
            }
        }

        double c[m-it][n-it] {};
        for(int i = 0; i < m-it; ++i)
            for(int j = 0; j < n-it; ++j)
                for(int k = 0; k < m-it; ++k)
                    c[i][j] = c[i][j] + u[i][k] * A[k+it][j+it];
        for(int i = 0; i < m-it; ++i)
            for(int j = 0; j < n-it; ++j)
                A[i+it][j+it] = c[i][j];

        Bd[it] = A[it][it];

        if(it < n-2){
            double v2[n-it-1];
            double beta2;
            double xr[n-it-1];
            for(int i = 0; i < n-it-1; ++i)
                xr[i] = A[it][i+it+1];
            House(xr, v2, beta2, n-it-1);
            double v[n-it-1][n-it-1];
            for(int i = 0; i < n-it-1; ++i){
                for(int j = 0; j < n-it-1; ++j){
                    if(i == j)
                        v[i][j] = 1.0;
                    else
                        v[i][j] = 0.0;
                    v[i][j] = v[i][j] - beta2 * v2[i] * v2[j];
                }
            }
            double d[m-it][n-it-1] {};
            for(int i = 0; i < m-it; ++i)
                for(int j = 0; j < n-it-1; ++j)
                    for(int k = 0; k < n-it-1; ++k)
                        d[i][j] = d[i][j] + A[i+it][k+it+1] * v[k][j];
            for(int i = 0; i < m-it; ++i)
                for(int j = 0; j < n-it-1; ++j)
                    A[i+it][j+it+1] = d[i][j];

            Bs[it] = A[it][it+1]; 
        }
        if(it == n-2)
            Bs[it] = A[it][it+1];
    }
 
    //if(iprint == 1)
    //    cout << "Matrix::GetBidiag(): end" << endl;
}

// Golub-Kahan SVD step
void GKSVD(double *Bd, double *Bs, const int n)
{
    //if(iprint == 1)
    //    cout << "Matrix::GKSVD(): start" << endl;

    double t1 = Bs[n-3] * Bs[n-3] + Bd[n-2] * Bd[n-2];//t(n-1, n-1)
    double t2 = Bs[n-2] * Bs[n-2] + Bd[n-1] * Bd[n-1];//t(n, n)
    double t3 = Bs[n-2] * Bd[n-2];//t(n, n-1)
    double d = (t1 - t2) / 2.0;
    double mu = t2 - t3 * t3 / (d + d / fabs(d) * sqrt(d * d + t3 * t3));

    double y = Bd[0] * Bd[0] - mu;
    double z = Bd[0] * Bs[0];
    double c, s;
    double a1 = 0.0;
    double a2 = 0.0;
    double G[2][2];
    for(int k = 0; k < n - 1; ++k){
        Givens(y, z, c, s);
        G[0][0] = c;
        G[0][1] = s;
        G[1][0] = -s;
        G[1][1] = c;
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

        double temp1[3][2] {};
        for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 2; ++j)
                for(int k = 0; k < 2; ++k)
                    temp1[i][j] = temp1[i][j] + b1[i][k] * G[k][j];
        for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 2; ++j)
                b1[i][j] = temp1[i][j];

        a1 = b1[2][0];
        if(k != 0) Bs[k-1] = b1[0][0];
        Bd[k] = b1[1][0];
        Bs[k] = b1[1][1];
        Bd[k+1] = b1[2][1];
        y = Bd[k];
        z = a1;
        Givens(y, z, c, s);
        G[0][0] = c;
        G[0][1] = -s;
        G[1][0] = s;
        G[1][1] = c;
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

        double temp2[2][3] {};
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 3; ++j)
                for(int k = 0; k < 2; ++k)
                    temp2[i][j] = temp2[i][j] + G[i][k] * b2[k][j];
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 3; ++j)
                b2[i][j] = temp2[i][j];

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

    //if(iprint == 1)
    //    cout << "Matrix::GKSVD(): end" << endl;
}

// SVD algorithm
void SVD(double **A, double *D, const int m, const int n)
{
    //if(iprint == 1)
    //    cout << "Matrix::SVD(): start" << endl;

    double Bd[n];
    double Bs[n-1];
    GetBidiag(A, Bd, Bs, m, n);
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
    while(q < n){
        //set b(i,i+1) to zero, if |b(i,i+1)| <= tol*(|b(i,i)|+|b(i+1,i+1)|)
        for(int i = 0; i < n-p-q-1; ++i){//relative index
            if(fabs(Bs[p+i]) <= tol * (fabs(Bd[p+i]) + fabs(Bd[p+i+1])))
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
            if(fabs(Bd[p+i]) <= tol * Bn){
                Bd[p+i] = 0.0;
                if(i != n-p-q-1 && Bs[p+i] != 0.0)
                    ZeroRow(Bd+p+i, Bs+p+i, n-p-q-i);
                if(i == bi){
                    bi = i + 1;
                    continue;
                }
                if(Bs[i-bi-1] != 0.0)
                    ZeroColumn(Bd+p+bi, Bs+p+bi, i-bi+1);
                GKSVD(Bd+p+bi, Bs+p+bi, i-bi);
                //Bs[p+i-1] = 0.0;
                bi = i + 1;
            }
        }
        if(bi != n-p-q){
            if(bi == n-p-q-1)
                continue;
            GKSVD(Bd+p+bi, Bs+p+bi, n-p-q-bi);
        }
    }

    //if(iprint == 1)
    //    cout << "Matrix::SVD(): end" << endl;
}

// zero 1st row of the bidiagonal matrix when the 1st diagonal entry is 0
void ZeroRow(double *Bd, double *Bs, const int n)
{
    //if(iprint == 1)
    //    cout << "Matrix::ZeroRow(): start" << endl;

    double a = Bs[0];
    for(int it = 1; it < n; ++it){
        double y = a;
        double z = Bd[it];
        double c, s;
        Givens(y, z, c, s);
        swap(c, s);
        c = -c;
        double G[2][2] {c, -s, s, c};
        double b[2][2] {y, 0.0, z, 0.0};
        if(it != n-1) b[1][1] = Bs[it];

        double temp[2][2] {};
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                for(int k = 0; k < 2; ++k)
                    temp[i][j] = temp[i][j] + G[i][k] * b[k][j];
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                b[i][j] = temp[i][j];

        if(it == 1) Bs[0] = 0.0;
        Bd[it] = b[1][0];
        if(it != n-1) Bs[it] = b[1][1];
        a = b[0][1];
    }

    //if(iprint == 1)
    //    cout << "Matrix::ZeroRow(): end" << endl;
}

// zero the last column if the last diagonal entry is 0
void ZeroColumn(double *Bd, double *Bs, const int n)
{
    //if(iprint == 1)
    //    cout << "Matrix::ZeroColumn(): start" << endl;

    double a = Bs[n-2];
    for(int it = n-2; it >= 0; --it){
        double y = Bd[it];
        double z = a;
        double c, s;
        Givens(y, z, c, s);
        double G[2][2] {c, s, -s, c};
        double b[2][2] {0.0, 0.0, y, z};
        if(it != 0) b[0][0] = Bs[it-1];

        double temp[2][2] {};
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                for(int k = 0; k < 2; ++k)
                    temp[i][j] = temp[i][j] + b[i][k] * G[k][j];
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                b[i][j] = temp[i][j];

        if(it == n-2) Bs[it] = 0.0;
        Bd[it] = b[1][0];
        if(it != 0) Bs[it-1] = b[0][0];
        a = b[0][1];
    }

    //if(iprint == 1)
    //    cout << "Matrix::ZeroColumn(): end" << endl;
}
