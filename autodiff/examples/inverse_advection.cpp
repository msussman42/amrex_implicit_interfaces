#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <adept_source.h>
#include <adept.h>

using namespace std;
using namespace amrex;
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
//
/*
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

*/

Real algorithm_and_gradient(
 vector< MultiFab* > x,	// control
 vector< MultiFab* > dJdx,
 const Real* dx,
 Real dt) { // gradient

 int nsteps=x.size()-1;

 adept::Stack stack;
 using adept::adouble;

 using My_adept_FAB=BaseFab< adept::adouble >;
 using My_adept_MFAB=FabArray< My_adept_FAB >;
 vector< My_adept_MFAB* > adept_x;  

 adept_x.resize(nsteps+1);

 adouble y;

 for (int ntime=0;ntime<=nsteps;ntime++) {
  BoxArray ba_node=x[ntime]->boxArray();
  DistributionMapping dm=x[ntime]->DistributionMap();
  int Ncomp=x[ntime]->nComp();
  int Nghost=x[ntime]->nGrow();

  adept_x[ntime]=new My_adept_MFAB(ba_node,dm,Ncomp,Nghost);

  My_adept_MFAB* cur_frame_adept=adept_x[ntime];
  MultiFab* cur_frame=x[ntime];
  for (MFIter mfi(*cur_frame_adept,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   My_adept_FAB& adept_fab=(*cur_frame_adept)[gridno];
   FArrayBox& Real_fab=(*cur_frame)[gridno];
   Array4< adept::adouble > const& adept_array=adept_fab.array();
   Array4<Real> const& Real_array=Real_fab.array();
   Dim3 lo = lbound(adept_array);
   Dim3 hi = ubound(adept_array);
   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
    adept_array(i,j,k)=Real_array(i,j,k);
   } // i
   } // j
   } // k
  } // mfi
 } // ntime

/*
 stack.new_recording();  // records the tape
 y = cost_function(adept_x);
 y.set_gradient(1.0);
 stack.compute_adjoint();
 for (int ntime=0;ntime<=Nsteps;ntime++) {
         for (int i=0;i<=Ncells;i++) {
		 dJdx[ntime][i]=adept_x[ntime][i].get_gradient();
	 }
 }
*/
 return y.value();
}


int main(int argc,char* argv[]) {

 amrex::Initialize(argc,argv);

 // What time is it now?  We'll use this to compute total run time.
 Real strt_time = amrex::second();

 // AMREX_SPACEDIM: number of dimensions
 int n_cell, max_grid_size, nsteps, plot_int;
 Real xlo=0.0;
 Real xhi=1.0;
 Real stop_time=1.0;

 // DEFAULT: NONE PERIODIC 
 // is_periodic(dir)=0 in all direction by default (dir=0,1,.., sdim-1)
 Vector<int> is_periodic(AMREX_SPACEDIM,0);  
     
 // inputs parameters
 {
  // ParmParse is way of reading inputs from the inputs file
  ParmParse pp;

  // We need to get n_cell from the inputs file - this is the 
  // number of cells on each side of 
  //   a square (or cubic) domain.
  pp.get("n_cell",n_cell);

  // The domain is broken into boxes of size max_grid_size
  pp.get("max_grid_size",max_grid_size);

  // Default plot_int to -1, allow us to set it to 
  // something else in the inputs file
  //  If plot_int < 0 then no plot files will be writtenq
  plot_int = -1;
  pp.query("plot_int",plot_int);

  // Default nsteps to 10, allow us to set it 
  // to something else in the inputs file
  nsteps = 10;
  pp.query("nsteps",nsteps);

  pp.get("xlo",xlo);
  pp.get("xhi",xhi);
  pp.get("stop_time",stop_time);
 }


 // make BoxArray and Geometry
 BoxArray ba;
 Geometry geom;
 {
  IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
  IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
  Box domain(dom_lo, dom_hi);

  // Initialize the boxarray "ba" from the single box "bx"
  ba.define(domain);
  // Break up boxarray "ba" into chunks no larger than 
  // "max_grid_size" along a direction
  ba.maxSize(max_grid_size);

  RealBox real_box({AMREX_D_DECL(xlo,xlo,xlo)},
                   {AMREX_D_DECL(xhi,xhi,xhi)});

  // This defines a Geometry object
  geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
 }
 BoxArray ba_node(ba);
 IndexType node_type=IndexType::NODE;
 ba_node.convert(node_type);

 // Nghost = number of ghost cells for each array 
 int Nghost = 0;
    
 // Ncomp = number of components for each array
 int Ncomp  = 1;
  
 // How Boxes are distrubuted among MPI processes
 DistributionMapping dm(ba);

 const Real* dx = geom.CellSize();
 Real dt=stop_time/nsteps;

/*
    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, 0);
    }
    if (plot_int > 0 && n%plot_int == 0)
    {
      const std::string& pltfile = amrex::Concatenate("plt",n,5);
      WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, n);
    }
*/

 vector< MultiFab* > x;  // control
 vector< MultiFab* > dJdx; 
 x.resize(nsteps+1);
 dJdx.resize(nsteps+1);

// BaseFab< adept::adouble > test_basefab();
// Array4< adept::adouble > const& a=test_basefab.array();


 for (int ntime=0;ntime<=nsteps;ntime++) {
	 x[ntime]=new MultiFab(ba_node,dm,Ncomp,Nghost);
	 dJdx[ntime]=new MultiFab(ba_node,dm,Ncomp,Nghost);

	 x[ntime]->setVal(0.0);
	 dJdx[ntime]->setVal(0.0);
 }

 // future: use Wolfe condition
 // future: use Fletcher Reeve's method (a.k.a. nonlinear CG)
 // future: use Gauss-Newton method
 // If the discretized constraint equation is a linear equation, then
 // the cost function will be convex which implies that an
 // optimal solution exists, AND, so long as the learnining rate
 // is sufficiently small, one can prove, using fixed point theory
 // that the steepest descent method is a ``contraction mapping'' AND
 // given any x, and hypercube [a,b] containing x,
 // then g(x)=x - r * dJdx(x) is also contained within
 // the same hypercube.
 //
 // future: grid stretching parameters and level set advection
 // velocity will be learned, instead
 // of resorting to AMR or standard level set methods.
 Real learning_rate=0.001;  // steepest descent parameter

 int max_iterations=100;
 for (int iter=0;iter<max_iterations;iter++) {
 
   Real y=algorithm_and_gradient(
             x,dJdx,
	     dx,dt);
 
   for (int ntime=0;ntime<=nsteps;ntime++) {
    MultiFab* cur_frame=x[ntime];
    MultiFab* cur_frame_dJdx=dJdx[ntime];
    for (MFIter mfi(*cur_frame,false); mfi.isValid(); ++mfi) {
     const int gridno = mfi.index();
     FArrayBox& x_fab=(*cur_frame)[gridno];
     FArrayBox& dJdx_fab=(*cur_frame_dJdx)[gridno];

     Array4<Real> const& x_array=x_fab.array();
     Array4<Real> const& dJdx_array=dJdx_fab.array();
     Dim3 lo = lbound(x_array);
     Dim3 hi = ubound(x_array);
     for (int k = lo.z; k <= hi.z; ++k) {
     for (int j = lo.y; j <= hi.y; ++j) {
     for (int i = lo.x; i <= hi.x; ++i) {
      x_array(i,j,k)=x_array(i,j,k)-learning_rate*dJdx_array(i,j,k);
     } // i
     } // j
     } // k
    } // mfi
   } // ntime
 } // iter

 Real real_stop_time = amrex::second() - strt_time;
 const int IOProc = ParallelDescriptor::IOProcessorNumber();
 ParallelDescriptor::ReduceRealMax(real_stop_time,IOProc);

 // Tell the I/O Processor to write out the "run time"
 amrex::Print() << "Run time = " << real_stop_time << std::endl;

 amrex::Finalize();
 return 0;

}
