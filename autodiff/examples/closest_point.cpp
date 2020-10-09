#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <adept_source.h>
#include <adept.h>

#define x0 (2.0)
#define y0 (1.0)

using namespace std;

adept::adouble fx_curve(adept::adouble x) {

	adept::adouble y=3.0*x+5.0;
	return y;

}

adept::adouble cost_function(adept::adouble x) {

 adept::adouble f_val=fx_curve(x);
 adept::adouble y=(x0-x)*(x0-x)+(f_val-y0)*(f_val-y0);
 return y;

}

double algorithm_and_gradient(
  const double x_val,
  double& dy_dx) {

 adept::Stack stack;
 using adept::adouble;
 adouble x = x_val;
 stack.new_recording();  // records the tape
 adouble y = cost_function(x);
 y.set_gradient(1.0);
 stack.compute_adjoint();
 dy_dx = x.get_gradient();
 return y.value();
}

int main() {

 double x_val;
 double dy_dx;
 x_val=3.0;

 double y=algorithm_and_gradient(x_val,dy_dx);
 std::cout << "x_val= " << x_val << '\n';
 std::cout << "y= " << y << '\n';
 std::cout << "dy_dx= " << dy_dx << '\n';


}
