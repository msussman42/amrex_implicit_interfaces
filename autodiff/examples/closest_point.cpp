// instructions:
// 1. get adept-2.0.5.tar.gz
// 2. gunzip adept-2.0.5.tar.gz
// 3. tar -xf adept-2.0.5.tar
// 4. there should be a directory called "adept-2.0.5"
// 5. mkdir examples
// 6. cd examples
// 7. create a program (e.g. closest_point.cpp) and make sure
//    to include <adept_source.h> and <adept.h>
// THIS CODE IS AN EXAMPLE USING AUTOMATIC DIFFERENTIATION FOR MINIMIZING THE FUNCTION
//   1. J(x)=(x-x0)^2 + (f(x)-y0)^2  using the fixed point method (also known as the steepest 
//      descent method)
//   2. x_closest=argmin_x J(x)  is the same problem as:
//      solve J'(x_closest)=0 if J''(x)>0 at the critical point and J is a convex function.
//   3. J'(x)=0  is equivalent to the fixed point problem:
//      x=x - r J'(x)   r= learning rate
//   4. Is g(x)=x-r J'(x) a contraction mapping?
//      g'(x)=1-r J''(x)  we had assumed in (item 2) just above that J''>0 (convex function) which
//      means for "r" sufficiently small, |g'(x)|<1. (prove this rigorously +2)
//
//   5. fixed point method is: xnp1 = xn - r J'(xn) 
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

// notice: "double" is replaced by "adept::adouble"
adept::adouble fx_curve(adept::adouble x) {

	adept::adouble y=3.0*x+5.0;
	return y;

}

// cost function for finding closest point on the curve y=fx_curve(x) from
// the given point (x0,y0)
adept::adouble cost_function(adept::adouble x) {

 adept::adouble f_val=fx_curve(x);
 adept::adouble y=(x0-x)*(x0-x)+(f_val-y0)*(f_val-y0);
 return y;

}

// why is this preferred over the complex step method for finding 
// derivatives? extra credit +2 ???
double algorithm_and_gradient(
  const double x_val,
  double& dy_dx) {

 adept::Stack stack;
 using adept::adouble;
 adouble x = x_val;
 stack.new_recording();  // records the tape
 adouble y = cost_function(x);
 y.set_gradient(1.0);
 stack.compute_adjoint();  // applies the chain rule: dJ/dx=dJ/dxn *
                           // dxn/dxn-1  * dxn-1/dxn-2  * ...  * dx1/dx
 dy_dx = x.get_gradient();
 return y.value();
}

int main() {

 double x_val;
 double dy_dx;
 x_val=3.0;

  // f(3)=3*3+5=14
  // J(3)=(2-3)^2+(14-1)^2=1+169=170
  // J'(x)=2(x0-x)+2(f(x)-y0)f'(x)=-2(2-3)+2(14-1)3=2+13(6)=80
 double y=algorithm_and_gradient(x_val,dy_dx);
 std::cout << "x_val= " << x_val << '\n';
 std::cout << "y= " << y << '\n';
 std::cout << "dy_dx= " << dy_dx << '\n';


}
