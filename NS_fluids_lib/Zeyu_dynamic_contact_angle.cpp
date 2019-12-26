#include <iostream>
#include <string>
#include <math.h>
#include <vector>
using namespace std;
#include "Zeyu_dynamic_contact_angle.h"

double delta(const double r, const int idelta)
{
    double d = 0.0;    

    if(idelta == 1){
        if(abs(r) <= 2.0)
            d = 1.0 / 4.0 * (1.0 + cos(asin(1.0) * r));
        else
            d = 0.0;
    }
    else if(idelta == 2){
        if(abs(r) <= 1.0)
            d = (3.0 - 2.0 * abs(r) + sqrt(1.0 + 4.0 * abs(r) - 4.0 * r * r)) / 8.0;
        else if(1.0 < abs(r) && abs(r) <= 2.0)
            d = (5.0 - 2.0 * abs(r) - sqrt(-7.0 + 12.0 * abs(r) - 4.0 * r * r)) / 8.0;
        else
            d = 0.0;
    }
    else
        cerr << "idelta index is wrong" << endl;

    return d;
}

void dynamic_contact_angle(const int dim, const double mu_l, const double mu_g, const double sigma, const double thet_s, const int imodel, const int ifgnbc, const int idelta, const int ifmicro, vector<double> &param, const vector<double> gridscale, const vector<Xcl> x_cl, const vector<vect3> n_cl, const vector<double> u_cl, const vector<double> thet_d_apparent, const vector<vect3> x_wall, vector<vect3> &u_slip, vector<double> &thet_d)
{
    if(ifgnbc == 0)
        cout << "Implement Dynamic Contact Angle..." << endl;
    else
        cout << "Implement Generalized Navier Boundary Condition..." << endl;

    vector<double> Ca;

    for(auto it = u_cl.begin(); it != u_cl.end(); ++it)
        Ca.push_back(mu_l * (*it) / sigma);

    if(imodel == 1){//Cox1986
        cout << "Implement model 1..." << endl;
        double &lamda = param[0];
        double &l_macro = param[1];
        double &l_micro = param[2];

        //if l_macro = 0, set it as the minimum grid length
        double minGridScale = 1.e3;
        for(auto it = gridscale.begin(); it != gridscale.end(); ++it)
            if(*it < minGridScale) minGridScale = *it;
        if(abs(l_macro - 0.0) < 1.e-9) l_macro = minGridScale;

        if(ifgnbc == 0){ //If don't implement GNBC.
            int index = 0;
            for(auto it = thet_d.begin(); it != thet_d.end(); ++it){
                *it = pow(pow(thet_s,3.0)+9.*Ca[index]*log(l_macro / l_micro),1./3.);
                ++index;
            }
            for(auto it = u_slip.begin(); it != u_slip.end(); ++it)
                *it = {0, 0, 0};
        }
        else{ //Implement GNBC.
            vector<double> thet_d_micro;
            double beta = mu_l / lamda;
            if(ifmicro != 0){ //Implement microscopic dynamic contact angle.
                double chi = (mu_l + mu_g) / (2.0 * beta * minGridScale);
                vector<double> thet_d_micro_old, Ca_old;
                double a = 9.0 * log(l_macro / l_micro);
                vector<double> b;
                double c = chi * cos(thet_s);

                int index = 0;
                for(auto it = thet_d_apparent.begin(); it != thet_d_apparent.end(); ++it){
                    int itermax = 100;
                    int iter = 0;
                    double erromax = 1.e-4;
                    thet_d_micro.push_back(0.0);
                    Ca[index] = 0.0;
                    thet_d_micro_old.push_back(*it);
                    Ca_old.push_back(chi * (cos(thet_s)-cos(*it)));
                    b.push_back(pow(*it, 3));
                
                    while(iter <= itermax){
                        double f1 = pow(thet_d_micro_old[index], 3) + a * Ca_old[index] - b[index];
                        double f2 = Ca_old[index] + chi * cos(thet_d_micro_old[index]) - c;
                        double Ja11 = 3.0 * pow(thet_d_micro_old[index], 2);
                        double Ja12 = a;
                        double Ja21 = - chi * sin(thet_d_micro_old[index]);
                        double Ja22 = 1.0;
                        double b1 = Ja11 * thet_d_micro_old[index] + Ja12 * Ca_old[index] - f1;
                        double b2 = Ja21 * thet_d_micro_old[index] + Ja22 * Ca_old[index] - f2;
                        double u11 = Ja11;
                        double u12 = Ja12;
                        double l21 = Ja21 / (u11 + 1.e-20);
                        double u22 = Ja22 - l21 * u12;
                        double y1 = b1;
                        double y2 = b2 - l21 * y1;
                        Ca[index] = y2 / (u22 + 1.e-20);
                        thet_d_micro[index] = (y1 - u12 * Ca[index]) / (u11 + 1.e-20);
                        if(abs((thet_d_micro[index] - thet_d_micro_old[index])/thet_d_micro_old[index]) < erromax && abs((Ca[index] - Ca_old[index])/Ca_old[index]) < erromax)
                            break;
                        thet_d_micro_old[index] = thet_d_micro[index];
                        Ca_old[index] = Ca[index];
                        ++iter;
                    }
                    cout << "Calculating Ca and thet_d_micro...\nnumber of iteration is: " << iter << endl;
                    ++index;
                }
            }
            int index1 = 0;
            for(auto it1 = u_slip.begin(); it1 != u_slip.end(); ++it1){
                double tempx = 0.;
                double tempy = 0.;
                double tempz = 0.;
                int index2 = 0;
                for(auto it2 = x_cl.begin(); it2 != x_cl.end(); ++it2){
                    double temp = 1. / (beta+1.e-20)  / (minGridScale+1.e-20) * delta(abs(x_wall[index1].x-(*it2).x)/(minGridScale+1.e-20),idelta) * sigma;
                    if(ifmicro == 0)
                        temp = temp * (cos(thet_s) - cos(thet_d_apparent[index2]));
                    else
                        temp = temp * (cos(thet_s) - cos(thet_d_micro[index2]));
                    if(dim == 3)
                        temp = temp * 1. / (minGridScale+1.e-20) * delta(abs(x_wall[index1].y-(*it2).y)/(minGridScale+1.e-20),idelta) * (*it2).l;
                    tempx = tempx + temp * n_cl[index2].x;
                    tempy = tempy + temp * n_cl[index2].y;//If 2D, tempy should be zero
                    tempz = tempz + temp * n_cl[index2].z;//tempz should be zero
                    thet_d[index2] = thet_d_apparent[index2];//In GNBC, dynamic contact angle isn't restrained?
                    ++index2;
                }
                (*it1).x = tempx;// x-component of slip velocity
                (*it1).y = tempy;// y-component of slip velocity, if 2D, should be zero
                (*it1).z = tempz;// z-component of slip velocity, should be zero
                ++index1;
            }
        }
    }
    else if(imodel == 2){//Jiang1970
        cout << "Implement model 2..." << endl;
        int index = 0;
        for(auto it = thet_d.begin(); it != thet_d.end(); ++it){
            *it = acos(cos(thet_s) - (1. + cos(thet_s)) * tanh(4.96 * pow(Ca[index], 0.702)));
            ++index;
        }
        for(auto it = u_slip.begin(); it != u_slip.end(); ++it)
            *it = {0, 0, 0};
    }
    else if(imodel == 3){//Shikmurzaev2008
        cout << "Implement model 3..." << endl;
        double a1, a2, a3, a4;
        vector<double> u, u0, thet_d_old;

        a2 = 0.54;
        a3 = 12.5;
        a4 = 0.07;
        a1 = 1. + (1. - a2) * (cos(thet_s) - a4);
        int index = 0;
        for(auto it = thet_d.begin(); it != thet_d.end(); ++it){
            u.push_back(a3 * Ca[index]);
            thet_d_old.push_back(thet_d_apparent[index]);
            u0.push_back(0.0);
            int iter;
            for(iter = 0; iter < 100; ++iter){
                u0[index] = (sin(thet_d_old[index] - thet_d_old[index] * cos(thet_d_old[index]))) 
                          / (sin(thet_d_old[index]) * cos(thet_d_old[index]) - thet_d_old[index]);
                *it = acos(cos(thet_s) - 2. * u[index] * (a1 + a2 * u0[index]) / (1. - a2) / (sqrt(a1 + u[index] * u[index]) + u[index]));
                if(abs((*it - thet_d_old[index])/thet_d_old[index]) < 1e-4)
                    break;
                else
                    thet_d_old[index] = *it;
            }
            cout << "Calculating thet_d...\nnumber of iteration is: " << iter << endl;
            ++index;
        }
        for(auto it = u_slip.begin(); it != u_slip.end(); ++it)
            *it = {0, 0, 0};
    }
    else if(imodel == 4){//Kalliadasis1994-'abs(tan(the_d))=...', so how to determine the_d<90 or thet_d>90? This one seems have problems...
        cout << "Implement model 4..." << endl;
        double lamda = 1.e-8;
        int index = 0;
        for(auto it = thet_d.begin(); it != thet_d.end(); ++it){
            *it = atan(7.48 * pow(Ca[index], 1./3.) - 3.28 * pow(lamda, 0.04) * pow(Ca[index], 0.293));
            ++index;
        }
        for(auto it = u_slip.begin(); it != u_slip.end(); ++it)
            *it = {0, 0, 0};
    }
    else if(imodel == 5){//Kistler1993
        cout << "Implement model 5..." << endl;
        double temp = pow(atanh((1. - cos(thet_s)) / 2.) / 5.16, 1./0.706);
        double fHI, fHI_old;
        fHI_old = temp / (1. - 1.31 * temp);
        int iter;
        for(iter = 0; iter < 100; ++iter){
            fHI = temp * (1. + 1.31 * pow(fHI_old, 0.99));
            if(abs((fHI - fHI_old)/fHI_old) < 1.e-4)
                break;
            else
                fHI_old = fHI;
        }
        cout << "Calculating fHI...\nnumber of iteration is: " << iter << endl;
        int index = 0;
        for(auto it = thet_d.begin(); it != thet_d.end(); ++it){
            *it = acos(1. - 2. * tanh(5.16 * pow((Ca[index] + fHI)/(1. + 1.31 * pow(Ca[index] + fHI, 0.99)), 0.706)));
            ++index;
        }
        for(auto it = u_slip.begin(); it != u_slip.end(); ++it)
            *it = {0, 0, 0};
    }
    else if(imodel == 6){//Bracke1989
        cout << "Implement model 6..." << endl;
        int index = 0;
        for(auto it = thet_d.begin(); it != thet_d.end(); ++it){
            *it = acos(cos(thet_s) - 2. * (1. + cos(thet_s)) * pow(Ca[index], 0.5));
            ++index;
        }
        for(auto it = u_slip.begin(); it != u_slip.end(); ++it)
            *it = {0, 0, 0};
    }
    else if(imodel == 7){//Blake2006, Popescu2008
        cout << "Implement model 7..." << endl;
        const double sigma_0 = param[0];
        const double v_0 = param[1];
        int index = 0;
        for(auto it = thet_d.begin(); it != thet_d.end(); ++it){
            *it = acos(cos(thet_s) - sigma_0 / sigma * asinh(u_cl[index] / v_0));
            ++index;
        }
        for(auto it = u_slip.begin(); it != u_slip.end(); ++it)
            *it = {0, 0, 0};
    }
    else
        cerr << "imodel index is wrong" << endl;

}
