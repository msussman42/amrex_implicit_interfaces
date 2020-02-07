#ifndef ZEYU_MATRIX_H
#define ZEYU_MATRIX_H

using namespace std;

template <class T>
using matrix = vector<vector<T>>;

class Matrix
{
public:

    //////////////
    // public member functions

    // constructor - construct a new zero matrix
    Matrix(const int i, const int j);
    Matrix(const matrix<double> m);
    ~Matrix();

    // add elements
    void Add(const int i, const int j, const double a);

    // return number of rows
    int ni(){return mm;}

    // return number of columns
    int nj(){return nn;}

    // transpose
    void Transpose();

    // calculate condition number
    double CondNum();

    void SetIPrint(const int i){iprint = i;}

private:

    /////////////
    // private member functions

    // calculate Householder vector
    void House(const vector<double> x, vector<double> &v, double &beta);

    // calculate Givens rotation matrix
    //    _     _ T _  _     _ _
    //   |  c  s | |  a |   | r |
    //                    = 
    //   |_-s  c_| |_ b_|   |_0_|
    //    
    void Givens(const double a, const double b, double &c, double &s);

    // Householder Bidiagonalization
    void GetBidiag(matrix<double> A, vector<double> &Bd, vector<double> &Bs);

    // Golub-Kahan SVD step
    void GKSVD(vector<double> &Bd, vector<double> &Bs);

    // SVD algorithm
    void SVD(matrix<double> A, vector<double> &D);

    // zero Row i of the bidiagonal matrix when the diagonal entry in that row is 0
    void ZeroRow(vector<double> &Bd, vector<double> &Bs, const int i);

    // zero the last column if the last diagonal entry is 0
    void ZeroColumn(vector<double> &Bd, vector<double> &Bs);

    // get transpose
    matrix<double> GetTranspose(const matrix<double> A);

    // get unity matrix
    matrix<double> GetUnity(const int a);
 
    ////////////
    // private data members

    int mm; //number of rows
    int nn; //number of column
    matrix<double> M; //the matrix

    int iprint = 0;

};

////////////
// define matrix operator

matrix<double> operator + (const matrix<double> a, const matrix<double> b);
matrix<double> operator - (const matrix<double> a, const matrix<double> b);
matrix<double> operator * (const matrix<double> a, const matrix<double> b);
matrix<double> operator * (const matrix<double> a, const double b);
matrix<double> operator * (const double b, const matrix<double> a);
double operator * (const vector<double> a, const vector<double> b);
matrix<double> operator & (const vector<double> a, const vector<double> b);
vector<double> operator * (const vector<double> a, const double b);
vector<double> operator * (const double b, const vector<double> a);
vector<double> operator * (const vector<double> a, const matrix<double> b);
vector<double> operator * (const matrix<double> b, const vector<double> a);
matrix<double> operator / (const matrix<double> a, const double b);
vector<double> operator / (const vector<double> a, const double b);

#endif
