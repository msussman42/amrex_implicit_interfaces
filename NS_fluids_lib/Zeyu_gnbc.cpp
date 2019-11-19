#include <iostream>
#include <string>
#include <math.h>
using namespace std;
#include "Zeyu_gnbc.h"

double d_2d(double r, double l_grid)
{
    double d = 0.0;    

    if(abs(r) <= 2.0 * l_grid)
        d = 1.0 / (4.0 * l_grid) * (1.0 + cos(2.0 * asin(1.0) * r / 2.0 / l_grid));
    else
        d = 0.0;

    return d;
}

double d_3d(double r)
{
    double d = 0.0;

    if(abs(r) <= 1.0)
        d = (3.0 - 2.0 * abs(r) + sqrt(1.0 + 4.0 * abs(r) - 4 * r * r)) / 8.0;
    else if(1.0 < abs(r) && abs(r) <= 2.0)
        d = 0.5 - (3.0 - 2.0 * abs(2.0 - abs(r)) + sqrt(1.0 + 4.0 * abs(2.0 - abs(r)) - 4 * (2.0 - abs(r)) * (2.0 - abs(r)))) / 8.0;
    else
        d = 0.0;

    return d;
}

// x_cl, y_cl, z_cl coordinate of closest point on the contact line (input)
// dim=2,3 input, mu_l=liquid viscosity, sigma=surface tension,
// theta_s=static angle, thet_d_grid grid scale dynamic angle,
// l_grid=characteristic length scale, l_macro (input),
// l_micro (input), chi (input), (x,y,z) coordinate,
// thet_d_macro=thet_d_grid?, thet_d_micro (output),
// Ca (output), u_slip (output - most essential variable)
void gnbc(const int dim, const double mu_l, const double sigma, const double thet_s, const double thet_d_grid, const double l_grid, const double l_macro, const double l_micro, const double chi, const double x_cl, const double y_cl, const double z_cl, const double x, const double y, const double z, double &thet_d_macro, double &thet_d_micro, double &Ca, double &u_slip)
{
    int itermax = 100;
    int iter = 0;
    double erromax = 1.e-4;
    double thet_d_micro_old, Ca_old;
    double a = 9.0 * log(l_macro / l_micro);
    double b = pow(thet_d_grid, 3);
    double c = chi * cos(thet_s);

    thet_d_micro_old = thet_d_grid;
    Ca_old = chi * (cos(thet_s)-cos(thet_d_micro_old));

    while(iter <= itermax){
        double f1 = pow(thet_d_micro_old, 3) + a * Ca_old - b;
        double f2 = Ca_old + chi * cos(thet_d_micro_old) - c;
        double Ja11 = 3.0 * pow(thet_d_micro_old, 2);
        double Ja12 = a;
        double Ja21 = - chi * sin(thet_d_micro_old);
        double Ja22 = 1.0;
        double b1 = Ja11 * thet_d_micro_old + Ja12 * Ca_old - f1;
        double b2 = Ja21 * thet_d_micro_old + Ja22 * Ca_old - f2;
        double u11 = Ja11;
        double u12 = Ja12;
        double l21 = Ja21 / (u11 + 1.e-20);
        double u22 = Ja22 - l21 * u12;
        double y1 = b1;
        double y2 = b2 - l21 * y1;
        Ca = y2 / (u22 + 1.e-20);
        thet_d_micro = (y1 - u12 * Ca) / (u11 + 1.e-20);
        if(abs((thet_d_micro - thet_d_micro_old)/thet_d_micro_old) < erromax && abs((Ca - Ca_old)/Ca_old) < erromax)
            break;
        thet_d_micro_old = thet_d_micro;
        Ca_old = Ca;
        ++iter;
    }

    cout << "number of iteration is: " << iter << endl; 
    thet_d_macro = pow(pow(thet_d_micro, 3) + a * Ca, 1.0/3.0);
    u_slip = Ca * sigma / mu_l;
    if(dim == 2)
        u_slip = u_slip * 2.0 * l_grid * d_2d(y - y_cl, l_grid);
    else if(dim == 3)
        u_slip = u_slip * 4.0 * d_3d((x - x_cl) / l_grid) * d_3d((y - y_cl) / l_grid);
    else
        cerr << "wrong input of dim: " << dim << endl;
}
