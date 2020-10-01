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

	double local_pi=4.0*atan(1.0;

	adept::adouble Hreturn;
	if (x<=-eps) {
		Hreturn=0.0;
	} else if (x>=eps) {
		Hreturn=1.0;
	} else {
		Hreturn=0.5*(x/eps+sin(local_pi*x/eps)/pi+1.0);
	}
	return Hreturn;
}

// this cost function will be:
//  q_t + (a q)_x =0   0<t<T   0<x<1  (*)
//  control variables:
//  q(0,x)=uinit(x)
//  q(t,0)=uleft(x)  if a(t,0)>0
//  q(t,1)=uright(x) if a(t,1)<0
//  J[uinit,uleft,uright]=integral_{x=0 to 1} 
//    w1(x)(q(T,x)-qobservation(T,x))^2 +
//    w2(x)(uinit(x)-uback_init(x))^2 +
//    w3(x)(uleft(x)-uback_left(x))^2 +
//    w4(x)(uright(x)-uback_right(x))^2 
//
//    q(T,x) solves the linear conservation law equation (*)
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
	 adept::adouble fluxes[Nnodes];
	 adept::adouble H_nodes_array[Nnodes];
	 adept::adouble H_flux_array[Ncells];
	 adept::adouble a_flux_array[Ncells];

	 for (int i=0;i<=Ncells;i++) {
          H_nodes_array[i]=H_smooth(a[i],eps);  // H_smooth is a differentiable 
	                                        // Heaviside function
	 }
	  // nodes: x_i=i*h  i=0...Ncells
	  // flux locations: x_{i+1/2}=(i+1/2)h  i=0..Ncells-1
	 for (int i=0;i<Ncells;i++) {
		 a_flux_array[i]=0.5*(a[i]+a[i+1]);
		 H_flux_array[i]=H_smooth(a_flux_array[i],eps);
		 fluxes[i]=(H_flux_array[i]*vn[i]+
			 (1.0-H_flux_array[i])*vn[i+1])*a_flux_array[i];
	 }
	  // nodes: x_i=i*h  i=0...Ncells
	  // flux locations: x_{i+1/2}=(i+1/2)h  i=0..Ncells-1
	  // i=0 .... i=Ncells  (finite difference domain)
	  // i=1 .... i=Ncells-1 (interior finite difference grid points)
         for (int i=1;i<Ncells;i++) {
		 vnp1[i]=vn[i]-(k/h)*(fluxes[i]-fluxes[i-1]);
	 }

	  // at i=0, if a>0, then vnp1=boundary value
	  //         if a<0, then vnp1=vn-(k/h)(fluxes[0]-a[0]*vn[0])
         adept::adouble left_flux=a[0]*vn[0];
	 adept::adouble numerical_vnp1_left=vn[0]-(k/h)*(fluxes[0]-left_flux);
	 adept::adouble bc_vnp1_left=uleft[ntime+1];
	 vnp1[0]=H_nodes_array[0]*bc_vnp1_left+
		 (1.0-H_nodes_array[0])*numerical_vnp1_left;

         adept::adouble right_flux=a[Ncells]*vn[Ncells];
	 adept::adouble numerical_vnp1_right= 
		 vn[Ncells]-(k/h)*(right_flux-fluxes[Ncells-1]);
	 adept::adouble bc_vnp1_right=uright[ntime+1];
	 vnp1[Ncells]=H_nodes_array[Ncells]*bc_vnp1_right+
		 (1.0-H_nodes_array[Ncells])*numerical_vnp1_right;

         for (int i=0;i<=Ncells;i++) {
                 vn[i]=vnp1[i];
	 }
 }

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
