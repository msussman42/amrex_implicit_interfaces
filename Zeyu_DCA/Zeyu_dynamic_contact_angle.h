#ifndef GNBC_H
#define GNBC_H

struct Xcl{
    //The center point coordinates and length of each contact line segment.
    double x, y, z, l; //In 2D case, l should be set to 1.
};

struct vect3{
    double x, y, z;
};

void dynamic_contact_angle(const int dim, const double mu_l, const double mu_g, const double sigma, const double thet_s, const int imodel, const int ifgnbc, const int idelta, const int ifmicro, vector<double> &param, const vector<double> gridscale, const vector<Xcl> x_cl, const vector<vect3> n_cl, const vector<double> u_cl, const vector<double> thet_d_apparent, const vector<vect3> x_wall, vector<vect3> &u_slip, vector<double> &thet_d);

#endif
