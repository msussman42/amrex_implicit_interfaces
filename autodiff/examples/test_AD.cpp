#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <adept_source.h>
#include <adept.h>

using namespace std;

adept::adouble cost_function(adept::adouble x[2]) {

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

adept::adouble dfdx(
  adept::adouble x_val[2]) {

 adept::Stack stack;
 using adept::adouble;
 adouble x[2];
 x[0]=x_val[0];
 x[1]=x_val[1];
 stack.new_recording();  // records the tape
 adouble y = cost_function(x);
 y.set_gradient(1.0);
 stack.compute_adjoint();
 adouble dydx=x[0].get_gradient();
 return dydx;
}


double algorithm_and_dxdx_dxdy(
  const double x_val[2],
  double dy_dx[2]) {

 adept::Stack stack;
 using adept::adouble;
 adouble x[2] = {x_val[0], x_val[1]};
 stack.new_recording();  // records the tape
 adouble y = dfdx(x);
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

 double dydx_dx[2];

 double dydx=algorithm_and_dxdx_dxdy(x_val,dydx_dx);
 std::cout << "x_val[0]= " << x_val[0] << '\n';
 std::cout << "x_val[1]= " << x_val[1] << '\n';
 std::cout << "dydx= " << dydx << '\n';
  // expecting dydx/dx0 = 2
  // expecting dydx/dx1 = 2
 std::cout << "dydx_dx[0]= " << dydx_dx[0] << '\n';
 std::cout << "dydx_dx[1]= " << dydx_dx[1] << '\n';


}
