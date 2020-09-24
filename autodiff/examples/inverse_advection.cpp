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
using namespace adept;

// This article seemed to have good success with smoothing the
// Heaviside function upwinding advection terms:
// Efficient evaluation of the direct and adjoint 
// linearized dynamics from compressible flow solvers
// MF De Pando, D Sipp, PJ Schmid
// Journal of Computational Physics 231 (23), 7739-7755

adept::adouble H_smooth(adept::adouble x,adept::adouble eps) {

	adept::adouble Hreturn;
	if (x<=-eps) {
		Hreturn=0.0;
	} else if (x>=eps) {
		Hreturn=1.0;
	} else {
		Hreturn=?
	}
	return Hreturn;
}

// this cost function will be:
//  q_t + a q_x =0   0<t<T   0<x<1  (*)
//  control variables:
//  q(0,x)=uinit(x)
//  q(t,0)=uleft(x)
//  q(t,1)=uright(x)
//  J[uinit,uleft,uright]=integral_{x=0 to 1} 
//    w1(x)(q(T,x)-qobservation(T,x))^2 +
//    w2(x)(uinit(x)-uback_init(x))^2 +
//    w3(x)(uleft(x)-uback_left(x))^2 +
//    w4(x)(uright(x)-uback_right(x))^2 
//
//    q(T,x) solves the advection equation (*)
//  
//  OUR APPROACH IN THIS EXAMPLE IS CALLED "discretize then differentiate"
//  referring to the review article by Giles and Pierce.
//  This is in contrast to "differentiate then discretize."
adept::adouble cost_function(adept::adouble uinit[Nnodes],
		adept::adouble uleft[Nsteps],
		adept::adouble uright[Nsteps]) {

 adept::adouble eps=1.0e-10;
 adept::adouble vn[Nnodes]; // no AD: "double vn[Nnodes]"
 adept::adouble vnp1[Nnodes];
 adept::adouble a[Nnodes];

 adept::adouble xlo=0.0;
 adept::adouble xhi=0.0;
 adept::adouble stop_time=1.0;
 adept::adouble h=(xhi-xlo)/Ncells;
 adept::adouble k=stop_time/Nsteps;

 for (int i=0;i<=Ncells;i++) {
	 adept::adouble xi=xlo+i*h;
	 vnp1[i]=unit[i];
	 a[i]=1.0;
 }

 for (int ntime=0;ntime<Nsteps;ntime++) {
	 adept::adouble DPLUS[Nnodes];
	 adept::adouble DMINUS[Nnodes];
	 for (int i=0;i<Ncells;i++) {
		 DPLUS[i]=(vn[i+1]-vn[i])/h;
		 DMINUS[i+1]=(vn[i+1]-vn[i])/h;
	 }
	 DPLUS[Ncells]=0.0;
	 DMINUS[0]=0.0;
	  // upwinding switch
	 adept::adouble H_array[Nnodes];
	 for (int i=0;i<=Ncells;i++) {
          H_array[i]=H_smooth(a[i],eps);  // H_smooth is a differentiable 
	                                  // Heaviside function
	 }
	 for (int i=0;i<=Ncells;i++) {
		 vn[i]=vnp1[i];
                 adept::adouble DM=vn[i]-vn[i-1];
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
