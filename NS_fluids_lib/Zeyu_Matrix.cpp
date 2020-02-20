#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

#include "Zeyu_Matrix.H"

// public member functions
// constructor - construct a new zero matrix
Matrix::Matrix(const int i, const int j)
{
    mm = i;
    nn = j;
    vector<double> row(nn, 0.0);
    for(int iter = 0; iter < mm; ++iter)
        M.push_back(row);
}

Matrix::Matrix(const matrix<double> m)
{
    M = m;
    mm = (int) M.size();
    nn = (int) M[0].size();
}

Matrix::~Matrix()
{
}

// add elements
void Matrix::Add(const int i, const int j, const double a)
{
    if(i > mm){
        vector<double> row(nn, 0.0);
        M.push_back(row);
        ++mm;
    }
    if(j > nn){
        for(int iter = 0; iter < mm; ++iter)
            M[iter].push_back(0.0);
        ++nn;
    }
    M[i-1][j-1] = a;
}

// transpose
void Matrix::Transpose()
{
    M = GetTranspose(M);
}

// calculate condition number
double Matrix::CondNum()
{
    if(iprint == 1)
        cout << "Matrix::CondNum(): start" << endl;

    vector<double> D;
    SVD(M, D);
    int n = (int) D.size();
    double min = 1.e20;
    double max = 0.0;
    for(int i = 0; i < n; ++i){
        if(std::abs(D[i]) > max) max = std::abs(D[i]);
        if(std::abs(D[i]) < min) min = std::abs(D[i]);
        //cout << "diag[" << i + 1 << "] = " << D[i] << endl;
    }

    if(iprint == 1)
        cout << "Matrix::CondNum(): end" << endl;
    
    if(std::abs(min - 0.0) < 1.e-16){
        cerr << "Error in calculating condition number: Zero singular value!" 
             << endl;
        abort();   
    }
    else
        return max / min;
}

// private member functions
// calculate Householder vector
void Matrix::House(const vector<double> x, vector<double> &v, double &beta)
{
    if(iprint == 1)
        cout << "Matrix::House(): start" << endl;

    double sigma;
    vector<double> xx (x);
    v.push_back(1.0);
    xx.erase(xx.begin());
    sigma = xx * xx;
    for(auto it = xx.begin(); it != xx.end(); ++it){
        v.push_back(*it);
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
        v = v / v0;
    }

    if(iprint == 1)
        cout << "Matrix::House(): end" << endl;
}

// calculate Givens rotation matrix
//    _     _ T _  _     _ _
//   |  c  s | |  a |2   | r |
//                    = 
//   |_-s  c_| |_ b_|   |_0_|
//
void Matrix::Givens(const double a, const double b, double &c, double &s)
{
    if(iprint == 1)
        cout << "Matrix::Givens(): start" << endl;

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

    if(iprint == 1)
        cout << "Matrix::Givens(): end" << endl;
}

// Householder Bidiagonalization
void  Matrix::GetBidiag(matrix<double> A, vector<double> &Bd, vector<double> &Bs)
{
    if(iprint == 1)
        cout << "Matrix::GetBidiag(): start" << endl;

    int m = (int) A.size();
    int n = (int) A[0].size();

    int it;
    for(it = 0; it < n; ++it){
         // std::vector< T > is a dynamic structure (SLOW!)
         // std::array<T,N> is a static structure (FAST!)
        vector<double> x;
        for(int i = 0; i < m - it; ++i)
            x.push_back(A[i][0]);
        vector<double> v1;
        double beta1;
        House(x, v1, beta1);
        matrix<double> u;
        u = GetUnity(m - it) - beta1 * (v1 & v1);
        A = u * A;
        Bd.push_back(A[0][0]);
        for(int i = 0; i < m - it; ++i)
            A[i].erase(A[i].begin());
        if(it < n-2){
            vector<double> v2;
            double beta2;
            House(A[0], v2, beta2);
            matrix<double> v;
            v = GetUnity(n - it -1) - beta2 * (v2 & v2);
            A = A * v;
            Bs.push_back(A[0][0]); 
            A.erase(A.begin());
        }
        if(it == n-2){
            Bs.push_back(A[0][0]);
            A.erase(A.begin());
        }
    }
 
    if(iprint == 1)
        cout << "Matrix::GetBidiag(): end" << endl;
}

// Golub-Kahan SVD step
void Matrix::GKSVD(vector<double> &Bd, vector<double> &Bs)
{
    if(iprint == 1)
        cout << "Matrix::GKSVD(): start" << endl;

    int n = (int) Bd.size();
    int ns = (int) Bs.size();
    if(ns != n - 1){
        cerr << "Error! Length of matrix diagonal is " << n
             << ", while length of matrix superdiagonal is " << ns
             << "(should be n - 1)" << endl;
        abort();
    }

    double t1 = Bs[n-3] * Bs[n-3] + Bd[n-2] * Bd[n-2];//t(n-1, n-1)
    double t2 = Bs[n-2] * Bs[n-2] + Bd[n-1] * Bd[n-1];//t(n, n)
    double t3 = Bs[n-2] * Bd[n-2];//t(n, n-1)
    double d = (t1 - t2) / 2.0;
    double mu = t2 - t3 * t3 / (d + d / std::abs(d) * sqrt(d * d + t3 * t3));

    double y = Bd[0] * Bd[0] - mu;
    double z = Bd[0] * Bs[0];
    double c, s;
    double a1 = 0.0;
    double a2 = 0.0;
     // Allocate G as 2x2 up here.
     // Allocate grow1 and grow2 as 2x1 up here for faster speed?
    matrix<double> G;
    vector<double> grow1, grow2;
    for(int k = 0; k < n - 1; ++k){
        Givens(y, z, c, s);
        grow1 = {c, s};
        grow2 = {-s, c};
        G = {grow1, grow2};
        matrix<double> b1;
        if(k == 0){
            vector<double> row1 {Bd[0], Bs[0]};
            vector<double> row2 {0.0, Bd[1]};
            b1 = {row1, row2};
        }
        else{
            vector<double> row1 {Bs[k-1], a2};
            vector<double> row2 {Bd[k], Bs[k]};
            vector<double> row3 {0.0, Bd[k+1]};
            b1 = {row1, row2, row3};
        }
        b1 = b1 * G;
        a1 = b1[b1.size()-1][0];
        if(k != 0) Bs[k-1] = b1[0][0];
        Bd[k] = b1[b1.size()-2][0];
        Bs[k] = b1[b1.size()-2][1];
        Bd[k+1] = b1[b1.size()-1][1];
        y = Bd[k];
        z = a1;
        Givens(y, z, c, s);
        grow1 = {c, -s};
        grow2 = {s, c};
        G = {grow1, grow2};
        matrix<double> b2;
        if(k == n - 2){
            vector<double> row1 {Bd[n-2], Bs[n-2]};
            vector<double> row2 {a1, Bd[n-1]};
            b2 = {row1, row2};
        }
        else{
            vector<double> row1 {Bd[k], Bs[k], 0.0};
            vector<double> row2 {a1, Bd[k+1], Bs[k+1]};
            b2 = {row1, row2};
        }
        b2 = G * b2;
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

    if(iprint == 1)
        cout << "Matrix::GKSVD(): end" << endl;
}

// SVD algorithm
void Matrix::SVD(matrix<double> A, vector<double> &D)
{
    if(iprint == 1)
        cout << "Matrix::SVD(): start" << endl;

    vector<double> Bd;
    vector<double> Bs;
    GetBidiag(A, Bd, Bs);
    D = Bd;
    double Bn = sqrt(Bd * Bd + Bs * Bs);
    int n = (int) Bd.size();
    int p = 0;
    int q = 0;
    int pp = p;
    int qq = q;
    double tol = 1.e-16;
    while(q < n){
        for(int i = 0; i < n-p-q-1; ++i){
            if(std::abs(Bs[i]) <= tol * (std::abs(Bd[i]) + std::abs(Bd[i+1])))
                Bs[i] = 0.0;
        }
        for(int i = 0; i < n-p-q-1; ++i){
            if(Bs[n-p-q-2-i] == 0.0){
                ++qq;
                D[n-1-q-i] = Bd[n-p-q-1-i];
                if(i == n-p-q-2){
                    ++qq;
                    D[n-2-q-i] = Bd[0];
                    qq = qq + pp;
                    pp = 0;
                }
                Bd.pop_back();
                Bs.pop_back();
            }
            else
                break;
        }
        q = qq;
        p = pp;
        if(q < n){
            for(int i = 0; i < n-p-q-1; ++i){
                if(Bs[i] == 0.0){
                    ++pp;
                    D[p+i] = Bd[i];
                    Bd.erase(Bd.begin());
                    Bs.erase(Bs.begin());
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
            if(std::abs(Bd[i]) <= tol * Bn){
                Bd[i] = 0.0;
                if(i != n-p-q-1 && Bs[i] != 0.0)
                    ZeroRow(Bd, Bs, i);
                if(i == bi){
                    bi = i + 1;
                    continue;
                }
                vector<double> Bdt(i-bi+1);
                vector<double> Bst(i-bi);
                for(int ii = bi; ii < i; ++ii){
                    Bdt[ii-bi] = Bd[ii];
                    Bst[ii-bi] = Bs[ii];
                }
                Bdt[i-bi] = Bd[i];
                if(Bst[i-bi-1] != 0.0)
                    ZeroColumn(Bdt, Bst);
                Bdt.pop_back();
                Bst.pop_back();
                GKSVD(Bdt, Bst);
                for(int ii = bi; ii < i-1; ++ii){
                    Bd[ii] = Bdt[ii-bi];
                    Bs[ii] = Bst[ii-bi];
                }
                Bd[i-1] = Bdt[i-bi-1];
                Bs[i-1] = 0.0;
                bi = i + 1;
            }
        }
        if(bi != n-p-q){
            if(bi == n-p-q-1)
                continue;
            vector<double> Bdt(n-p-q-bi);
            vector<double> Bst(n-p-q-1-bi);
            for(int ii = bi; ii < n-p-q-1; ++ii){
                Bdt[ii-bi] = Bd[ii];
                Bst[ii-bi] = Bs[ii];
            }
            Bdt[n-p-q-1-bi] = Bd[n-p-q-1];
            GKSVD(Bdt, Bst);
            for(int ii = bi; ii < n-p-q-1; ++ii){
                Bd[ii] = Bdt[ii-bi];
                Bs[ii] = Bst[ii-bi];
            }
            Bd[n-p-q-1] = Bdt[n-p-q-bi-1];
        }
    }

    if(iprint == 1)
        cout << "Matrix::SVD(): end" << endl;
}

// zero a row of the bidiagonal matrix when the diagonal entry in that row is 0
void Matrix::ZeroRow(vector<double> &Bd, vector<double> &Bs, const int i)
{
    if(iprint == 1)
        cout << "Matrix::ZeroRow(): start" << endl;

    double a = Bs[i];
    int n = (int) Bd.size();
    for(int ii = i+1; ii < n; ++ii){
        double y = a;
        double z = Bd[ii];
        double c, s;
        Givens(y, z, c, s);
        swap(c, s);
        c = -c;
        vector<double> grow1 {c, -s};
        vector<double> grow2 {s, c};
        matrix<double> G {grow1, grow2};
        vector<double> row1 {y, 0.0};
        vector<double> row2;
        if(ii == n-1)
            row2 = {z, 0.0};
        else
            row2 = {z, Bs[ii]};
        matrix<double> b {row1, row2};
        b = G * b;
        if(ii == i+1) Bs[i] = 0.0;
        Bd[ii] = b[1][0];
        if(ii != n-1) Bs[ii] = b[1][1];
        a = b[0][1];
    }

    if(iprint == 1)
        cout << "Matrix::ZeroRow(): end" << endl;
}

// zero the last column if the last diagonal entry is 0
void Matrix::ZeroColumn(vector<double> &Bd, vector<double> &Bs)
{
    if(iprint == 1)
        cout << "Matrix::ZeroColumn(): start" << endl;

    int n = (int) Bd.size();
    double a = Bs[n-2];
    for(int i = n-2; i >= 0; --i){
        double y = Bd[i];
        double z = a;
        double c, s;
        Givens(y, z, c, s);
        vector<double> grow1 {c, s};
        vector<double> grow2 {-s, c};
        matrix<double> G {grow1, grow2};
        vector<double> row1;
        if(i == 0)
            row1 = {0.0, 0.0};
        else
            row1 = {Bs[i-1], 0.0};
        vector<double> row2 {y, z};
        matrix<double> b {row1, row2};
        b = b * G;
        if(i == n-2) Bs[i] = 0.0;
        Bd[i] = b[1][0];
        if(i != 0) Bs[i-1] = b[0][0];
        a = b[0][1];
    }

    if(iprint == 1)
        cout << "Matrix::ZeroColumn(): end" << endl;
}

// get transpose
matrix<double> Matrix::GetTranspose(const matrix<double> A)
{
    if(iprint == 1)
        cout << "Matrix::GetTranspose(): start" << endl;

    int m = (int) A.size();
    int n = (int) A[0].size();
    matrix<double> AT;
    for(int j = 0; j < n; ++j){
        vector<double> temp;
        for(int i = 0; i < m; ++i)
            temp.push_back(A[i][j]);
        AT.push_back(temp);
    }
    if(iprint == 1)
        cout << "Matrix::GetTranspose(): end" << endl;

    return AT;
}

//get unity matrix
matrix<double> Matrix::GetUnity(const int a)
{
    if(iprint == 1)
        cout << "Matrix::GetUnity(): start" << endl;

    matrix<double> unity;
    for(int i = 0; i < a; ++i){
        vector<double> row;
        for(int j = 0; j < a; ++j){
            if(i == j)
                row.push_back(1.0);
            else
                row.push_back(0.0);
        }
        unity.push_back(row);
    }
    if(iprint == 1)
        cout << "Matrix::GetUnity(): end" << endl;

    return unity;
}

/////////////////
// define matrix and vector operator

matrix<double> operator + (const matrix<double> a, const matrix<double> b)
{
    matrix<double> c;
    int ma = (int) a.size();
    int na = (int) a[0].size();
    int mb = (int) b.size();
    int nb = (int) b[0].size();
    
    if(ma != na || mb != nb){
        cerr << "Error in plus operator! Operand types are Matrix(" 
             << ma << " x " << na << ") and Matrix(" 
             << mb << " x " << nb << ")." << endl;
        abort();
    }

    for(int i = 0; i < ma; ++i){
        vector<double> row;
        for(int j = 0; j < na; ++j)
            row.push_back(a[i][j] + b[i][j]);
        c.push_back(row);
    }
    return c;
}

matrix<double> operator - (const matrix<double> a, const matrix<double> b)
{
    matrix<double> c;
    int ma = (int) a.size();
    int na = (int) a[0].size();
    int mb = (int) b.size();
    int nb = (int) b[0].size();
    
    if(ma != na || mb != nb){
        cerr << "Error in minus operator! Operand types are Matrix(" 
             << ma << " x " << na << ") and Matrix(" 
             << mb << " x " << nb << ")." << endl;
        abort();
    }

    for(int i = 0; i < ma; ++i){
        vector<double> row;
        for(int j = 0; j < na; ++j)
            row.push_back(a[i][j] - b[i][j]);
        c.push_back(row);
    }
    return c;
}

matrix<double> operator * (const matrix<double> a, const matrix<double> b)
{
    int ma = (int) a.size();
    int na = (int) a[0].size();
    int mb = (int) b.size();
    int nb = (int) b[0].size();
    if(na != mb){
        cerr << "Error in multiplication! Operand types are Matrix(" 
             << ma << " x " << na << ") and Matrix(" 
             << mb << " x " << nb << ")." << endl;
        abort();
    }
    matrix<double> c;
    for(int i = 0; i < ma; ++i){
        vector<double> row;
        for(int j = 0; j < nb; ++j){
            vector<double> column;
            for(int k = 0; k < mb; ++k)
                column.push_back(b[k][j]);
            row.push_back(a[i] * column);
        }
        c.push_back(row);
    }
    return c;
}

matrix<double> operator * (const matrix<double> a, const double b)
{
    int m = (int) a.size();
    int n = (int) a[0].size();
    matrix<double> c = a;
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < n; ++j)
            c[i][j] = a[i][j] * b;
    return c;
}

matrix<double> operator * (const double b, const matrix<double> a)
{
    return a * b;
}

double operator * (const vector<double> a, const vector<double> b)
{
    if(a.size() != b.size()){
        cerr << "Error in dot product! Operand types are Vector[" 
             << a.size() << "] and Vector[" << b.size() << "]." << endl;
        abort();
    }
    double c = 0.;
    for(int i = 0; i < (int) a.size(); ++i)
        c += a[i] * b[i];
    return c;
}

matrix<double> operator & (const vector<double> a, const vector<double> b)
{
    int m = (int) a.size();
    int n = (int) b.size();
    matrix<double> c;
    for(int i = 0; i < m; ++i){
        vector<double> row;
        for(int j = 0; j < n; ++j)
            row.push_back(a[i] * b[j]);
        c.push_back(row);
    }
    return c;
}

vector<double> operator * (const vector<double> a, const double b)
{
    vector<double> c = a;
    int n = (int) a.size();
    for(int i = 0; i < n; ++i)
        c[i] = a[i] * b;
    return c;
}

vector<double> operator * (const double b, const vector<double> a)
{
    return a * b;
}

vector<double> operator * (const vector<double> a, const matrix<double> b)
{
    int n = (int) a.size();
    int mb = (int) b.size();
    int nb = (int) b.size();
    if(n != mb){
        cerr << "Error in multiplication! Operand types are Vector[" 
             << n << "] and Matrix(" << mb << " x " << nb << ")." << endl;
        abort();
    }
    vector<double> c;
    for(int i = 0; i < nb; ++i){
        double s = 0;
        for(int j = 0; j < mb; ++j)
            s = s + a[j] * b[j][i];
        c.push_back(s);
    }
    return c;
}

vector<double> operator * (const matrix<double> b, const vector<double> a)
{
    int n = (int) a.size();
    int mb = (int) b.size();
    int nb = (int) b.size();
    if(n != nb){
        cerr << "Error in multiplication! Operand types are Matrix(" 
             << mb << " x " << nb << ") and Vector[" << n << "]." << endl;
        abort();
    }
    vector<double> c;
    for(int i = 0; i < mb; ++i){
        double s = 0;
        for(int j = 0; j < nb; ++j)
            s = s + a[j] * b[i][j];
        c.push_back(s);
    }
    return c;
}

matrix<double> operator / (const matrix<double> a, const double b)
{
    if(b == 0){
        cerr << "Error! Denominator is zero!" << endl;
        abort();
    }
    return a * (1 / b);
}

vector<double> operator / (const vector<double> a, const double b)
{
    if(b == 0){
        cerr << "Error! Denominator is zero!" << endl;
        abort();
    }
    return a * (1 / b);
}
