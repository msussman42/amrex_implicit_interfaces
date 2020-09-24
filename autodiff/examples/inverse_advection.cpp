#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <adept_source.h>
#include <adept.h>

#define Ncells 64
#define Nnodes (Ncells+1)
#define Nsteps 128

using namespace std;

adept::adouble cost_function(adept::adouble uinit[Nnodes],
		adept::adouble uleft[Nsteps],
		adept::adouble uright[Nsteps]) {

 adept::adouble vn[Nnodes]; // no AD: "double vn[Nnodes]"
 adept::adouble vnp1[Nnodes];

 adept::adouble xlo=0.0;
 adept::adouble xhi=0.0;
 adept::adouble stop_time=1.0;
 adept::adouble h=(xhi-xlo)/Ncells;
 adept::adouble a=1.0;
 adept::adouble k=stop_time/Nsteps;

 adept::adouble y=x[0]*x[0] + x[1]*x[1];
 return y;

}

double algorithm_and_gradient(
  const double x_val[2],
  double dy_dx[2]) {

 adept::Stack stack;
 using adept::adouble;
 adouble x[2] = {x_val[0], x_val[1]};
 stack.new_recording();  // records the tape
 adouble y = cost_function(x);
 y.set_gradient(1.0);
 stack.compute_adjoint();
 dy_dx[0] = x[0].get_gradient();
 dy_dx[1] = x[1].get_gradient();
 return y.value();
}

int main() {

 double x_val[2];
 double dy_dx[2];
 x_val[0]=3.0;
 x_val[1]=4.0; // expecting y=25.0

 double y=algorithm_and_gradient(x_val,dy_dx);
 std::cout << "x_val[0]= " << x_val[0] << '\n';
 std::cout << "x_val[1]= " << x_val[1] << '\n';
 std::cout << "y= " << y << '\n';
  // expecting dy/dx0 = 6
  // expecting dy/dx1 = 8
 std::cout << "dy_dx[0]= " << dy_dx[0] << '\n';
 std::cout << "dy_dx[1]= " << dy_dx[1] << '\n';


}
