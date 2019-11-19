#include <iostream>
#include <math.h>
using namespace std;

#include "Zeyu_read_input_file.h"
#include "Zeyu_gnbc.h"

/*
    g++ -std=c++11 -O3 -Wall Zeyu_main.cpp Zeyu_read_input_file.cpp  Zeyu_gnbc.cpp -o main_gnbc
*/

int main()
{
    int dim;
    double thet_s, thet_d_grid, thet_d_macro, thet_d_micro, l_grid, l_macro, l_micro, chi, u_slip, Ca;
    double mu_l, sigma;
    double x_cl, y_cl, z_cl;//position of contact line. For 2D case, assume contact line(point) moves along y direction.
                                                      //For 3D case, assmue contact line moves on plane x-z.

    //read input file
    read_input_file(dim, mu_l, sigma, thet_s, thet_d_grid, l_grid, l_macro, l_micro, chi, x_cl, y_cl, z_cl);
    thet_s = thet_s / 180.0 * 2.0 * asin(1.0);
    thet_d_grid = thet_d_grid / 180.0 * 2.0 * asin(1.0);
    if(dim == 2)
        cout << "y_cl = " << y_cl << ", l_grid = " << l_grid << endl;
    else if(dim == 3)
        cout << "x_cl = " << x_cl << ", z_cl = " << z_cl << ", l_grid = " << l_grid << endl;
    else
        cerr << "wrong input of dim: " << dim << endl;

    double x, y, z;//coordinates of boundary wall grid point where you want to calculate slip velocity.
    cout << "Please input coordinates of a point(x y z)\n(If it is a 2D case, input (0 y 0), if it is a 3D case, input (x 0 z)):" << endl;
    cin >> x >> y >> z;

    gnbc(dim, mu_l, sigma, thet_s, thet_d_grid, l_grid, l_macro, l_micro, chi, x_cl, y_cl, z_cl, x, y, z, thet_d_macro, thet_d_micro, Ca, u_slip);

    cout << "\nThe macroscopic dynamic contact angle is: " << thet_d_macro / (2.0 * asin(1.0)) * 180.0 << endl;
    cout << "The microscopic dynamic contact angle is: " << thet_d_micro / (2.0 * asin(1.0)) * 180.0 << endl;
    cout << "The nondimensional contact line velocity Ca is: " << Ca << endl;
    cout << "The magnitude of slip velocity of point (" << x << " " << y << " " << z << ") is: " << u_slip << endl;

    return 0;
}
