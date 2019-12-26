#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

#include "Zeyu_read_input_file_dynamic_contact_angle.h"
#include "Zeyu_dynamic_contact_angle.h"

/*
    g++ -std=c++11 -O3 -Wall Zeyu_main_dynamic_contact_angle.cpp Zeyu_read_input_file_dynamic_contact_angle.cpp  Zeyu_dynamic_contact_angle.cpp -o main_dynamic_contact_angle
*/

int main()
{
    int dim, imodel, ifgnbc, idelta, ifmicro;
    double mu_l, mu_g, sigma, thet_s;
    vector<double> param, u_cl, thet_d_apparent, thet_d;
    vector<Xcl> x_cl;
    vector<vect3> n_cl, x_wall, u_slip;

    //read input file
    read_input_file(dim, mu_l, mu_g, sigma, thet_s, imodel, ifgnbc, idelta, ifmicro, param);

    //suspose values of some variables from simulation, initialize
    vector<double> gridscale(20,3.2e-5);
    if(dim ==2){//2D case
        //u_cl comes from simulation.(Popescu2008: advancing u_cl>0, receding u_cl<0?)
        u_cl.push_back(0.01);
        //If implement GNBC, thet_d_apparent comes from simulation.
        //If implement dynamic contact angle models, we don't use it.
        thet_d_apparent.push_back(50. / 180.0 * 2.0 * asin(1.0));
        //thet_d should be initialized outside of the subroutine.
        //If implement GNBC, we don't need it.
        //IF implement dynamic contact angle, it is the output.
        thet_d.push_back(0);
        //x_cl comes from simulation.
        x_cl.push_back({0.10002, 0, 0, 1});
        //n_cl comes from simulation.
        n_cl.push_back({1, 0, 0});
        for(int num_of_grid_around_interface = 0; num_of_grid_around_interface < 6; ++num_of_grid_around_interface){
            //x_wall comes from simulation.
            x_wall.push_back({0.1+(num_of_grid_around_interface-2)*3.2e-5, 0, 0});
            //u_slip should be initialized outside of the subroutine.
            //If implement GNBC, it is the output.
            //If implement dynamic contact angle, it is zero.
            u_slip.push_back({0, 0, 0});
        }
    }
    else{//3D case
        for(int num_of_interface_segments = 0; num_of_interface_segments < 10; ++num_of_interface_segments){
            u_cl.push_back(0.01);
            thet_d_apparent.push_back(50. / 180.0 * 2.0 * asin(1.0));
            thet_d.push_back(0);
            x_cl.push_back({0.10002, num_of_interface_segments*3.2e-5, 0, 3.2e-5});
            n_cl.push_back({1, 0, 0});
        }
        for(int num_of_grid_around_interface_in_x_direction = 0; num_of_grid_around_interface_in_x_direction < 6; ++num_of_grid_around_interface_in_x_direction){
            for(int num_of_interface_segments = 1; num_of_interface_segments < 9; ++num_of_interface_segments){
                x_wall.push_back({0.1+(num_of_grid_around_interface_in_x_direction-2)*3.2e-5, num_of_interface_segments*3.2e-5, 0});
                u_slip.push_back({0, 0, 0});
            }
        }
    }

    thet_s = thet_s / 180.0 * 2.0 * asin(1.0);

/***************************************/
    //INPUT FROM INPUT FILE:
    //dim: dimension of the simulation domain
    //mu_l: viscocity of liquid
    //mu_g: viscocity of gas
    //sigma: surface tension coefficient
    //thet_s: static contact angle
    //imodel: dynamic angle model index
    //ifgnbc: if implement GNBC
    //idelta: delta function index
    //ifmicro: if implement microscopic dynamic angle
    //param(vector<double>): parameters for different dynamic contact angle
/**************************************/
    //INPUT FROM SIMULATION:
    //gridscale(vector<double>): grid scale(I set lx=ly=lz for each grid, which means there's only one value of gridscale for each grid)
    //x_cl(vector<Xcl>): the coordinates of the center points of each contact line segment, and the length of this segment(x, y, z, l)
    //n_cl(vector<vect3>): unit normal vector of each contact line segment(x, y, z)
    //u_cl(vector<vect3>): velocity of each contact line segment(x, y, z)
    //thet_d_apparent(vector<double>): dynamic contact angle of each contact line segment from simulation
    //x_wall(vector<vect3>): points surrounding the contact line(x, y, z)
/***************************************/
    //OUTPUT:
    //u_slip(vector<vect3>)(for gnbc): slip velocity of each contact line segment(x, y, z)
    //thet_d(vector<vect3>)(for dynamic contact angle models): dynamic contact angle of each contact line segment(x, y, z)
/***************************************/
    //vector containers of x_cl, n_cl, u_cl, thet_d_apparent, and thet_d have the same size
    //vector containers of x_wall and  u_slip have the same size
    //vector container of gridscale is in the size of the quantity of all grids
    dynamic_contact_angle(dim, mu_l, mu_g, sigma, thet_s, imodel, ifgnbc, idelta, ifmicro, param, gridscale, x_cl, n_cl, u_cl, thet_d_apparent, x_wall, u_slip, thet_d);

    for(auto it = thet_d.begin(); it != thet_d.end(); ++it)
        cout << "The dynamic contact angle is: " << (*it) / (2.0 * asin(1.0)) * 180.0 << endl;
    int index = 0;
    for(auto it = u_slip.begin(); it != u_slip.end(); ++it){
        cout << "The slip velocity of point (" << x_wall[index].x << ", " << x_wall[index].y << ", " << x_wall[index].z <<") is: (" << (*it).x << ", " << (*it).y << ", " << (*it).z << ")"<< endl;
        ++index;
    }

    return 0;
}
