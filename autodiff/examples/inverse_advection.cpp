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

	double local_pi=4.0*atan(1.0);

	adept::adouble Hreturn;
	if (x<=-eps) {
		Hreturn=0.0;
	} else if (x>=eps) {
		Hreturn=1.0;
	} else {
		Hreturn=0.5*(x/eps+sin(local_pi*x/eps)/local_pi+1.0);
	}
	return Hreturn;
}
// to express the data assimilation problem precisely:
// find min J(q(u),u) under the constraint that
//          N(q(u),u)=0   q(u)=state variables u=control
// The constraint equations are:
//  N(q(u),u):   q_t + (a q)_x =0   0<t<T   0<x<1  (*)
//  (q(t,x)-u(t,x))*wt_background(t,x)=0
//  control variable:
//  u(t,x)
//  J[u,q[u]]=integral_{x=0 to 1} integral_{t=0..T}
//    wt_observation(t,x)(q(t,x)-q_observation(t,x))^2 +
//    wt_background(t,x)(ucontrol(t,x)-uback(t,x))^2
//
//    q(t,x) solves the linear conservation law equation (*)
// 
//
//  OUR APPROACH IN THIS EXAMPLE IS CALLED "discretize then differentiate"
//  referring to the review article by Giles and Pierce.
//  This is in contrast to "differentiate then discretize."
adept::adouble cost_function(
	 vector< vector< adept::adouble > > uinput) {

 adept::adouble local_pi=4.0*atan(1.0);
 adept::adouble eps=1.0e-10;
 vector< vector< adept::adouble > > v;  // approximate solution
 vector< vector< adept::adouble > > a;
 vector< vector< adept::adouble > > ubackground;
 vector< vector< adept::adouble > > wt_background;
 vector< vector< adept::adouble > > q_observation;
 vector< vector< adept::adouble > > wt_observation;

 adept::adouble xlo=0.0;
 adept::adouble xhi=0.0;
 adept::adouble stop_time=1.0;
 adept::adouble h=(xhi-xlo)/Ncells;
 adept::adouble k=stop_time/Nsteps;

  // AMReX has "Array4" designed for GPU, OMP, MPI
  // Can one declare Array4< adept::adouble > u  ???
 ubackground.resize(Nsteps+1);
 wt_background.resize(Nsteps+1);
 q_observation.resize(Nsteps+1);
 wt_observation.resize(Nsteps+1);
 for (int ntime=0;ntime<=Nsteps;ntime++) {
	 adept::adouble local_t=ntime*k;

	 ubackground[ntime].resize(Nnodes);
	 wt_background[ntime].resize(Nnodes);
	 q_observation[ntime].resize(Nnodes);
	 wt_observation[ntime].resize(Nnodes);
         for (int i=0;i<=Ncells;i++) {
	  adept::adouble local_x=xlo+h*i;

          ubackground[ntime][i]=0.0;
          wt_background[ntime][i]=0.0;
	  if (ntime==0) {
           wt_background[ntime][i]=1.0;
	  }
	  if (i==0) {
           wt_background[ntime][i]=1.0;
          }
	  if (i==Ncells) {
           wt_background[ntime][i]=1.0;
          }
          wt_observation[ntime][i]=0.0;
          q_observation[ntime][i]=0.0;
	  if (ntime==Nsteps) {
           q_observation[ntime][i]=sin(2.0*local_pi*local_x);
	  }
	 }
 }

 v.resize(Nsteps+1);
 for (int ntime=0;ntime<=Nsteps;ntime++) {
	 v[ntime].resize(Nnodes);
         for (int i=0;i<=Ncells;i++) {
		 v[ntime][i]=uinput[ntime][i];
 	 }
 }
 a.resize(Nsteps+1);
 for (int ntime=0;ntime<=Nsteps;ntime++) {
	 a[ntime].resize(Nnodes);
         for (int i=0;i<=Ncells;i++) {
		 a[ntime][i]=1.0;
 	 }
 }


 adept::adouble y=0.0;

 for (int ntime=0;ntime<Nsteps;ntime++) {
	 adept::adouble fluxes[Nnodes];
	 adept::adouble H_nodes_array[Nnodes];
	 adept::adouble H_flux_array[Ncells];
	 adept::adouble a_flux_array[Ncells];

	 for (int i=0;i<=Ncells;i++) {
	  // H_smooth is a differentiable Heaviside function
          H_nodes_array[i]=H_smooth(a[ntime][i],eps); 
	 }
	  // nodes: x_i=i*h  i=0...Ncells
	  // flux locations: x_{i+1/2}=(i+1/2)h  i=0..Ncells-1
	 for (int i=0;i<Ncells;i++) {
		 a_flux_array[i]=0.5*(a[ntime][i]+a[ntime][i+1]);
		 H_flux_array[i]=H_smooth(a_flux_array[i],eps);
		 fluxes[i]=(H_flux_array[i]*v[ntime][i]+
  		  (1.0-H_flux_array[i])*v[ntime][i+1])*a_flux_array[i];
	 }
	  // nodes: x_i=i*h  i=0...Ncells
	  // flux locations: x_{i+1/2}=(i+1/2)h  i=0..Ncells-1
	  // i=0 .... i=Ncells  (finite difference domain)
	  // i=1 .... i=Ncells-1 (interior finite difference grid points)
         for (int i=1;i<Ncells;i++) {
		 v[ntime+1][i]=v[ntime][i]-(k/h)*(fluxes[i]-fluxes[i-1]);
	 }

	  // at i=0, if a>0, then vnp1=boundary value
	  //         if a<0, then vnp1=vn-(k/h)(fluxes[0]-a[0]*vn[0])
	  // left_flux is a numerical flux assumed to be located at
	  // x=-h/2
         adept::adouble left_flux=a[ntime][0]*v[ntime][0];
	 adept::adouble numerical_vnp1_left=
		 v[ntime][0]-(k/h)*(fluxes[0]-left_flux);
	 adept::adouble bc_vnp1_left=uinput[ntime+1][0];
	 v[ntime+1][0]=H_nodes_array[0]*bc_vnp1_left+
		 (1.0-H_nodes_array[0])*numerical_vnp1_left;

	  // right_flux is a numerical flux assumed to be located at
	  // x=1+h/2
         adept::adouble right_flux=a[ntime][Ncells]*v[ntime][Ncells];
	 adept::adouble numerical_vnp1_right= 
		 v[ntime][Ncells]-(k/h)*(right_flux-fluxes[Ncells-1]);
	 adept::adouble bc_vnp1_right=uinput[ntime+1][Ncells];
	 v[ntime+1][Ncells]=H_nodes_array[Ncells]*bc_vnp1_right+
		 (1.0-H_nodes_array[Ncells])*numerical_vnp1_right;

         for (int i=0;i<=Ncells;i++) {
		 adept::adouble local_wt=0.0;
		 if (wt_background[ntime+1][i]>0.0)
			 local_wt=1.0;
                 v[ntime+1][i]=local_wt*uinput[ntime+1][i]+ 
			 (1.0-local_wt)*v[ntime+1][i];
     	         y=y+wt_background[ntime+1][i]*h*
                  (uinput[ntime+1][i]-ubackground[ntime+1][i])*
                  (uinput[ntime+1][i]-ubackground[ntime+1][i]);
     	         y=y+wt_observation[ntime+1][i]*h*
	   	  (q_observation[ntime+1][i]-v[ntime+1][i])*
	   	  (q_observation[ntime+1][i]-v[ntime+1][i]);
	 }
 }

 return y;

}

double algorithm_and_gradient(
 vector< vector< double > > x,	// control
 vector< vector< double > > dJdx) { // gradient

 adept::Stack stack;
 using adept::adouble;
 vector< vector< adept::adouble > > adept_x;  

 adept_x.resize(Nsteps+1);
 for (int ntime=0;ntime<=Nsteps;ntime++) {
	 adept_x[ntime].resize(Nnodes);
         for (int i=0;i<=Ncells;i++) {
		 adept_x[ntime][i]=x[ntime][i];
 	 }
 }
 stack.new_recording();  // records the tape
 adouble y = cost_function(adept_x);
 y.set_gradient(1.0);
 stack.compute_adjoint();
 for (int ntime=0;ntime<=Nsteps;ntime++) {
         for (int i=0;i<=Ncells;i++) {
		 dJdx[ntime][i]=adept_x[ntime][i].get_gradient();
	 }
 }
 return y.value();
}

int main() {

 vector< vector< double > > x;  // control
 vector< vector< double > > dJdx;
 dJdx.resize(Nsteps+1);
 x.resize(Nsteps+1);
 for (int ntime=0;ntime<=Nsteps;ntime++) {
	 dJdx[ntime].resize(Nnodes);
	 x[ntime].resize(Nnodes);
         for (int i=0;i<=Ncells;i++) {
		 x[ntime][i]=0.0;
		 dJdx[ntime][i]=0.0;
 	 }
 }

 // future: use Wolfe condition
 // future: use Fletcher Reeve's method (a.k.a. nonlinear CG)
 // If the discretized constraint equation is a linear equation, then
 // the cost function will be convex which implies that an
 // optimal solution exists, AND, so long as the learnining rate
 // is sufficiently small, one can prove, using fixed point theory
 // that the steepest descent method is a ``contraction mapping'' AND
 // given any x, and hypercube [a,b] containing x,
 // then g(x)=x - r * dJdx(x) is also contained within
 // the same hypercube.
 double learning_rate=0.001;  // steepest descent parameter

 int max_iterations=100;
 for (int iter=0;iter<max_iterations;iter++) {
   double y=algorithm_and_gradient(x,dJdx);
   for (int ntime=0;ntime<=Nsteps;ntime++) {
         for (int i=0;i<=Ncells;i++) {
          x[ntime][i]=x[ntime][i]-learning_rate*dJdx[ntime][i];
	 }
   }
 }


}
