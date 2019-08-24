#include <fstream>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

void run(int i, double psi);

void barbelldist(double* x, double* y, double* z,
                 double* xcen1, double* ycen1, double* zcen1,
                 double* xcen2, double* ycen2, double* zcen2,
                 double* r1, double* r2, double* r3, double* dist, double* phi);

void spheredist(double* x, double* y, double* z,
                double* xcen, double* ycen, double* zcen,
                double* r, double* phi);

void matmult(double A[][3], double B[][3], double C[][3]);


int main() {
    const static double pi = 3.1415;

    for (int i = 0; i < 100; i++) {
        double psi = (double) i * pi / 180.0;
        run(i, psi);
        cout << "Running step " << i << endl;
    }

    return 0;
}

void run(int i, double psi) {
    ofstream fout;
    ostringstream oss;
    oss << i;
    string file = "data" + oss.str() + ".tec";
    fout.open(file.c_str());

    fout << "TITLE=\"CLSVOF data\"\nVARIABLES=\"X\",\"Y\",\"Z\",\"LS\"" << endl;
    fout << "zone i=65 j=65 k=65 f=point" << endl;

    double dx = 1.0 / 64.0, dy = 1.0 / 64.0, dz = 1.0 / 64.0;
    //double xcen1 = 0.2, ycen1 = 0.5, zcen1 = 0.5;
    //double xcen2 =  0.8, ycen2 = 0.5, zcen2 = 0.5;
    double xcen1 = 0.3, ycen1 = 0.5, zcen1 = 0.5;
    double xcen2 = 0.7, ycen2 = 0.5, zcen2 = 0.5;
    double r1 = 0.1, r2 = 0.1, r3 = 0.03;
    double x, y, z;
    double dist = sqrt(pow(xcen1 - xcen2, 2) + pow(ycen1 - ycen2,
                       2) + pow(zcen1 - zcen2, 2));
    double phi;

    double theta = 0;
    //double psi = 20.0 * pi / 180.0;
    double gamma = 0;

    double xv[3];
    double yv[3];
    double cenMasses[3];

    cenMasses[0] = (xcen1 + xcen2) / 2.0;
    cenMasses[1] = (ycen1 + ycen2) / 2.0;
    cenMasses[2] = (zcen1 + zcen2) / 2.0;

    double A[3][3];
    double B[3][3];
    double C[3][3];
    double D[3][3];
    double E[3][3];

    A[0][0] = cos(theta);
    A[0][1] = -sin(theta);
    A[0][2] = 0;
    A[1][0] = sin(theta);
    A[1][1] = cos(theta);
    A[1][2] = 0;
    A[2][0] = 0;
    A[2][1] = 0;
    A[2][2] = 1.0;

    B[0][0] = cos(psi);
    B[0][1] = 0;
    B[0][2] = -sin(psi);
    B[1][0] = 0;
    B[1][1] = 1.0;
    B[1][2] = 0;
    B[2][0] = sin(psi);
    B[2][1] = 0;
    B[2][2] = cos(psi);

    C[0][0] = 1.0;
    C[0][1] = 0;
    C[0][2] = 0;
    C[1][0] = 0;
    C[1][1] = cos(gamma);
    C[1][2] = -sin(gamma);
    C[2][0] = 0;
    C[2][1] = sin(gamma);
    C[2][2] = cos(gamma);

    for (int i = 1; i <= 65; i++) {
        for (int j = 1; j <= 65; j++) {
            for (int k = 1; k <= 65; k++) {
                x = (i - 1.0) * dx;
                y = (j - 1.0) * dy;
                z = (k - 1.0) * dz;

                matmult(A, B, D);
                matmult(D, C, E);

                xv[0] = x;
                xv[1] = y;
                xv[2] = z;

                for (int l = 0; l < 3; l++) {
                    yv[l] = 0;

                    for (int m = 0; m < 3; m++) {
                        yv[l] += E[l][m] * (xv[m] - cenMasses[m]);
                    }

                    yv[l] += cenMasses[l];
                }

                barbelldist(&yv[0], &yv[1], &yv[2], &xcen1, &ycen1, &zcen1,
                            &xcen2, &ycen2, &zcen2, &r1, &r2, &r3, &dist, &phi);

                fout << x << " " << y << " " << z << " " <<
                     phi << endl;
                //spheredist(x, y, z, xcen1, ycen1, zcen1, r1) << endl;
            }
        }
    }

    fout.close();
}

void barbelldist(double* x, double* y, double* z,
                 double* xcen1, double* ycen1, double* zcen1,
                 double* xcen2, double* ycen2, double* zcen2,
                 double* r1, double* r2, double* r3, double* dist, double* phi) {

    double phi1, phi2, phi3;

    phi1 = sqrt(pow(*x - *xcen1, 2) + pow(*y - *ycen1, 2) + pow(*z - *zcen1,
                2)) - *r1;
    phi2 = sqrt(pow(*x - *xcen2, 2) + pow(*y - *ycen2, 2) + pow(*z - *zcen2,
                2)) - *r2;
    phi3 = sqrt(pow(*y - 0.5, 2) + pow(*z - 0.5, 2)) - *r3;

    if (*x < *xcen1) {
        *phi = phi1;
    }
    else if (*x >= *xcen2) {
        *phi = phi2;
    }
    else if (*xcen1 <= *x && *x <= ((*xcen1 + *xcen2) / 2)) {
        double r = sqrt(pow(*y - *ycen1, 2) + pow(*z - *zcen1, 2));

        if (r > *r3 && *x <= (*xcen1 + *r1)) {
            *phi = phi1;
        }
        else if (phi1 <= 0) {
            *phi = -sqrt(pow(phi1, 2) + pow(phi3, 2));
        }
        else {
            *phi = phi3;
        }
    }
    else if (((*xcen1 + *xcen2) / 2) <= *x && *x <= *xcen2) {
        double r = sqrt(pow(*y - *ycen2, 2) + pow(*z - *zcen2, 2));

        if (r > *r3 && *x >= *xcen2 - *r2) {
            *phi = phi2;
        }
        else if (phi2 < 0) {
            *phi = -sqrt(pow(phi2, 2) + pow(phi3, 2));
        }
        else {
            *phi = phi3;
        }
    }

}

void matmult(double A[][3], double B[][3], double C[][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[i][j] = 0;

            for (int k = 0; k < 3; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}



void spheredist(double* x, double* y, double* z,
                double* xcen, double* ycen, double* zcen,
                double* r, double* phi) {

    *phi = sqrt(pow(*x - *xcen, 2) + pow(*y - *ycen, 2) + pow(*z - *zcen, 2)) - *r;
}