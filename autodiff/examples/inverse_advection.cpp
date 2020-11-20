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
using My_adept_FAB=BaseFab< adept::adouble >;
using My_adept_MFAB=FabArray< My_adept_FAB >;

// Tracking deformable objects with unscented Kalman filtering 
// and geometric active contours
// This article seemed to have good success with smoothing the
// Heaviside function upwinding advection terms:
// Efficient evaluation of the direct and adjoint 
// linearized dynamics from compressible flow solvers
// MF De Pando, D Sipp, PJ Schmid
// Journal of Computational Physics 231 (23), 7739-7755

adept::adouble H_smooth(adept::adouble x,adept::adouble eps) {

//	double local_pi=4.0*atan(1.0);

	adept::adouble Hreturn;
         // "Logistic function"
        Hreturn=1.0/(1.0+exp(-x/eps));
	return Hreturn;
}

// to express the data assimilation problem precisely:
// find min J(q(u),u) under the constraint that
//          N(q(u),u)=0   q(u)=state variables u=control
// The constraint equations are:
//  N(q(u),u):   q_t + (a q)_x =0   0<t<T   0<x<1  (*)
//  (q(t,x)-u(t,x))*wt_back(t,x)=0
//  control variable:
//  u(t,x)
//  J[u,q[u]]=integral_{x=0 to 1} integral_{t=0..T}
//    wt_obs(t,x)(q(t,x)-q_obs(t,x))^2 +
//    wt_back(t,x)(ucontrol(t,x)-uback(t,x))^2
//
//    q(t,x) solves the linear conservation law equation (*)
// 
//
//  OUR APPROACH IN THIS EXAMPLE IS CALLED "discretize then differentiate"
//  referring to the review article by Giles and Pierce.
//  This is in contrast to "differentiate then discretize."
//
// this sample code solves PDE constrained optimization, first pass:
// PDE constraint is q_{t} + a q_{x}=0  (*)
// inverse problem: find q_{0}(x) (initial conditions) and 
// boundary conditions so that
// the solution to (*) at time t=T is q(x,T)=q^{observation}(x)

void plot_adept_data(std::string plt_name_string,
	int ntime,int steepest_descent_iter,
	const Geometry& geom,
	My_adept_MFAB* data_to_plot,double time) {

 const std::string& pltfile1 = amrex::Concatenate(plt_name_string,ntime,5);
 const std::string& pltfile = 
	 amrex::Concatenate(pltfile1,steepest_descent_iter,5);
 BoxArray ba_cell=data_to_plot->boxArray();
 DistributionMapping dm=data_to_plot->DistributionMap();
 int Ncomp=data_to_plot->nComp();
 int Nghost=data_to_plot->nGrow();
 MultiFab* data_plot_mf=new MultiFab(ba_cell,dm,Ncomp,Nghost);

 for (MFIter mfi(*data_plot_mf,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   My_adept_FAB& v_fab=(*data_to_plot)[gridno]; // type: adouble
   FArrayBox& plot_fab=(*data_plot_mf)[gridno]; // type: double

   Array4< adept::adouble > const& v_array=v_fab.array(); // type: adouble

   Array4< Real > const& plot_array=plot_fab.array(); // type: double
   Array4< Real > const& plot_cell_array=plot_cell_fab.array(); //type:double
//   Dim3 lo = lbound(v_array);
//   Dim3 hi = ubound(v_array);
   Dim3 lo=lbound(plot_fab.box());
   Dim3 hi=ubound(plot_fab.box());
   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
    plot_array(i,j,k)=v_array(i,j,k).value();
   }
   }
   }
 }

 std::cout << "plotting " << pltfile << '\n';
 int nghost_cell=1;
 bool local=false;
 std::cout << " max= " << data_plot_mf->max(0,nghost_cell,local) << '\n';
 std::cout << " min= " << data_plot_mf->min(0,nghost_cell,local) << '\n';
 WriteSingleLevelPlotfile(pltfile, *data_plot_mf, 
		 {"data_plot_mf"}, 
		 geom, time, 0);
 delete data_plot_mf;
} // plot_adept_data

adept::adouble cost_function(
	 vector< My_adept_MFAB* > uinput,
	 const Real* dx,
	 Real dt,
	 const Geometry& geom,
	 int plot_int,
	 int steepest_descent_iter) {

 int nsteps=uinput.size()-1;

 adept::adouble local_pi=4.0*atan(1.0);

// recommended eps from MF De Pando, D Sipp, PJ Schmid
// Journal of Computational Physics 231 (23), 7739-7755
 adept::adouble eps=sqrt(1.0e-15);

 vector< My_adept_MFAB* > v;  // approximate solution
 vector< My_adept_MFAB* > a;
 vector< My_adept_MFAB* > uback;
 vector< My_adept_MFAB* > wt_back;
 vector< My_adept_MFAB* > q_obs;
 vector< My_adept_MFAB* > wt_obs;

 v.resize(nsteps+1);
 a.resize(nsteps+1);
 uback.resize(nsteps+1);
 wt_back.resize(nsteps+1);
 q_obs.resize(nsteps+1);
 wt_obs.resize(nsteps+1);

 int ntime=0;
 BoxArray ba=uinput[ntime]->boxArray();  // cell centered box(es)

 BoxArray ba_flux_x(ba);
 IndexType flux_x_type=IndexType::TheUMACType();
 ba_flux_x.convert(flux_x_type);

 DistributionMapping dm=uinput[ntime]->DistributionMap();
 int Ncomp=uinput[ntime]->nComp();
 int Nghost=uinput[ntime]->nGrow();

 adept::adouble xlo=geom.ProbLo(0);
 adept::adouble xhi=geom.ProbHi(1);;

 adept::adouble stop_time=dt*nsteps;

 for (ntime=0;ntime<=nsteps;ntime++) {
  adept::adouble local_t=ntime*dt;

  v[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  a[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  uback[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  wt_back[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  q_obs[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  wt_obs[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);

  My_adept_MFAB* uinput_frame=uinput[ntime];

  My_adept_MFAB* v_frame=v[ntime];
  My_adept_MFAB* a_frame=a[ntime];
  My_adept_MFAB* uback_frame=uback[ntime];
  My_adept_MFAB* wt_back_frame=wt_back[ntime];
  My_adept_MFAB* q_obs_frame=q_obs[ntime];
  My_adept_MFAB* wt_obs_frame=wt_obs[ntime];

  for (MFIter mfi(*v_frame,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   My_adept_FAB& uinput_fab=(*uinput_frame)[gridno];

   My_adept_FAB& v_fab=(*v_frame)[gridno];
   My_adept_FAB& a_fab=(*a_frame)[gridno];
   My_adept_FAB& uback_fab=(*uback_frame)[gridno];
   My_adept_FAB& wt_back_fab=(*wt_back_frame)[gridno];
   My_adept_FAB& q_obs_fab=(*q_obs_frame)[gridno];
   My_adept_FAB& wt_obs_fab=(*wt_obs_frame)[gridno];

   Array4< adept::adouble > const& uinput_array=uinput_fab.array();

   Array4< adept::adouble > const& v_array=v_fab.array();
   Array4< adept::adouble > const& a_array=a_fab.array();
   Array4< adept::adouble > const& uback_array=
	   uback_fab.array();
   Array4< adept::adouble > const& wt_back_array=
	   wt_back_fab.array();
   Array4< adept::adouble > const& q_obs_array=
	   q_obs_fab.array();
   Array4< adept::adouble > const& wt_obs_array=
	   wt_obs_fab.array();

//   Dim3 lo = lbound(v_array);
//   Dim3 hi = ubound(v_array);
   Dim3 lo=lbound(v_fab.box());
   Dim3 hi=ubound(v_fab.box());

   for (int k = lo.z; k <= hi.z[2]; ++k) {
   for (int j = lo.y; j <= hi.y[1]; ++j) {
   for (int i = lo.x; i <= hi.x[0]; ++i) {

    adept::adouble local_x=xlo+dx[0]*(i+0.5);

    uback_array(i,j,k)=0.0;
    wt_back_array(i,j,k)=0.0;
    if (ntime==0) {
     wt_back_array(i,j,k)=1.0;
    }
    if (local_x<=xlo) {
     wt_back_array(i,j,k)=1.0;
    }
    if (local_x>=xhi) {
     wt_back_array(i,j,k)=1.0;
    }
    wt_obs_array(i,j,k)=0.0;
    q_obs_array(i,j,k)=0.0;
    if (ntime==nsteps) {
     q_obs_array(i,j,k)=sin(2.0*local_pi*local_x);
    }
    v_array(i,j,k)=uinput_array(i,j,k);
    a_array(i,j,k)=1.0;
   } // i
   } // j
   } // k
  } // mfi
 } // ntime


 adept::adouble y=0.0;

 double time=0.0;
 ntime = 0;
 My_adept_MFAB* v_frame_plot=v[ntime];
 plot_adept_data("plt",ntime,steepest_descent_iter,geom,v_frame_plot,
		 time);

 for (ntime=0;ntime<nsteps;ntime++) {

  My_adept_MFAB* a_frame=a[ntime];
  My_adept_MFAB* v_frame=v[ntime];
  My_adept_MFAB* vnp1_frame=v[ntime+1];
  My_adept_MFAB* unp1_frame=uinput[ntime+1];

  My_adept_MFAB* wtnp1_back_frame=wt_back[ntime+1];
  My_adept_MFAB* unp1_back_frame=uback[ntime+1];
  My_adept_MFAB* qnp1_obs_frame=q_obs[ntime+1];
  My_adept_MFAB* wtnp1_obs_frame=wt_obs[ntime+1];

   // finite volume method, unknowns located at cell centers.
  My_adept_MFAB* fluxes=new My_adept_MFAB(ba_flux_x,dm,Ncomp,Nghost);
  My_adept_MFAB* H_fluxes=new My_adept_MFAB(ba_flux_x,dm,Ncomp,Nghost);
  My_adept_MFAB* a_fluxes=new My_adept_MFAB(ba_flux_x,dm,Ncomp,Nghost);
  My_adept_MFAB* H_cells=new My_adept_MFAB(ba,dm,Ncomp,Nghost);

  for (MFIter mfi(*H_cells,false); mfi.isValid(); ++mfi) {

   const int gridno = mfi.index();

   My_adept_FAB& a_fab=(*a_frame)[gridno];
   My_adept_FAB& v_fab=(*v_frame)[gridno];
   My_adept_FAB& vnp1_fab=(*vnp1_frame)[gridno];
   My_adept_FAB& unp1_fab=(*unp1_frame)[gridno];

   My_adept_FAB& wtnp1_back_fab=(*wtnp1_back_frame)[gridno];
   My_adept_FAB& unp1_back_fab=(*unp1_back_frame)[gridno];
   My_adept_FAB& qnp1_obs_fab=(*qnp1_obs_frame)[gridno];
   My_adept_FAB& wtnp1_obs_fab=(*wtnp1_obs_frame)[gridno];

   My_adept_FAB& fluxes_fab=(*fluxes)[gridno];
   My_adept_FAB& H_fluxes_fab=(*H_fluxes)[gridno];
   My_adept_FAB& H_nodes_fab=(*H_nodes)[gridno];
   My_adept_FAB& a_fluxes_fab=(*a_fluxes)[gridno];

   Array4< adept::adouble > const& a_array=a_fab.array();
   Array4< adept::adouble > const& v_array=v_fab.array();
   Array4< adept::adouble > const& vnp1_array=vnp1_fab.array();
   Array4< adept::adouble > const& unp1_array=unp1_fab.array();

   Array4< adept::adouble > const& wtnp1_back_array=wtnp1_back_fab.array();
   Array4< adept::adouble > const& unp1_back_array=unp1_back_fab.array();
   Array4< adept::adouble > const& qnp1_obs_array=qnp1_obs_fab.array();
   Array4< adept::adouble > const& wtnp1_obs_array=wtnp1_obs_fab.array();

   Array4< adept::adouble > const& fluxes_array=fluxes_fab.array();
   Array4< adept::adouble > const& H_fluxes_array=H_fluxes_fab.array();
   Array4< adept::adouble > const& H_nodes_array=H_nodes_fab.array();
   Array4< adept::adouble > const& a_fluxes_array=a_fluxes_fab.array();

   Dim3 lo = lbound(v_array);
   Dim3 hi = ubound(v_array);

   for (int k = lo.z-Nghost_array[2]; k <= hi.z+Nghost_array[2]; ++k) {
   for (int j = lo.y-Nghost_array[1]; j <= hi.y+Nghost_array[1]; ++j) {
   for (int i = lo.x-Nghost_array[0]; i <= hi.x+Nghost_array[0]; ++i) {
    // H_smooth is a differentiable Heaviside function
    H_cells_array(i,j,k)=H_smooth(a_array(i,j,k),eps); 
	 
    // nodes: x_i=i*h  i=lo.x ... hi.x
    // flux locations: x_{i+1/2}=(i+1/2)h  i=lo.x ... hi.x-1
    if (i<hi.x) {
     a_fluxes_array(i,j,k)=
        0.5*(a_array(i,j,k)+a_array(i+1,j,k));
     H_fluxes_array(i,j,k)=
	H_smooth(a_fluxes_array(i,j,k),eps);
     fluxes_array(i,j,k)=
        (H_fluxes_array(i,j,k)*v_array(i,j,k)+
 	 (1.0-H_fluxes_array(i,j,k))*v_array(i+1,j,k))*
	a_fluxes_array(i,j,k);
    }
   } //i
   } //j
   } //k
   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
    if ((i>lo.x)&&(i<hi.x)) {
     // nodes: x_i=i*h  i=lo.x .. hi.x
     // flux locations: x_{i+1/2}=(i+1/2)h  i=lo.x .. hi.x-1
     vnp1_array(i,j,k)=v_array(i,j,k)-
	     (dt/dx[0])*(fluxes_array(i,j,k)-fluxes_array(i-1,j,k));
    }
   } //i
   } //j
   } //k

    // at i=0, if a>0, then vnp1=boundary value
    //         if a<0, then vnp1=vn-(k/h)(fluxes[0]-a[0]*vn[0])
    // left_flux is a numerical flux assumed to be located at
    // x=-h/2
   int i=lo.x;
   int j=0;
   int k=0;
   adept::adouble left_flux=a_array(i,j,k)*v_array(i,j,k);
   adept::adouble numerical_vnp1_left=
    v_array(i,j,k)-(dt/dx[0])*(fluxes_array(i,j,k)-left_flux);
   adept::adouble bc_vnp1_left=unp1_array(i,j,k);
   vnp1_array(i,j,k)=
      H_nodes_array(i,j,k)*bc_vnp1_left+
      (1.0-H_nodes_array(i,j,k))*numerical_vnp1_left;

    // right_flux is a numerical flux assumed to be located at
    // x=1+h/2
   i=hi.x;
   adept::adouble right_flux=a_array(i,j,k)*v_array(i,j,k);
   adept::adouble numerical_vnp1_right= 
    v_array(i,j,k)-(dt/dx[0])*(right_flux-fluxes_array(i-1,j,k));
   adept::adouble bc_vnp1_right=unp1_array(i,j,k);
   vnp1_array(i,j,k)= 
     (1.0-H_nodes_array(i,j,k))*bc_vnp1_right+
      H_nodes_array(i,j,k)*numerical_vnp1_right;

   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
    adept::adouble local_wt=0.0;
    if (wtnp1_back_array(i,j,k)>0.0)
     local_wt=1.0;
    vnp1_array(i,j,k)=
	 local_wt*unp1_array(i,j,k)+ 
	 (1.0-local_wt)*vnp1_array(i,j,k);
    y=y+wtnp1_back_array(i,j,k)*dx[0]*
        (unp1_array(i,j,k)-unp1_back_array(i,j,k))*
        (unp1_array(i,j,k)-unp1_back_array(i,j,k));

    y=y+wtnp1_obs_array(i,j,k)*dx[0]*
        (qnp1_obs_array(i,j,k)-vnp1_array(i,j,k))*
        (qnp1_obs_array(i,j,k)-vnp1_array(i,j,k));
   } // i
   } // j
   } // k

   delete fluxes;
   delete H_fluxes;
   delete a_fluxes;
   delete H_nodes;

  } // mfi

  if (ntime==nsteps-1) {
   plot_adept_data("plt_obs",ntime+1,steepest_descent_iter,geom,
		 qnp1_obs_frame,
		 time);
  }

 } //ntime=0..nsteps-1

 for (ntime=0;ntime<=nsteps;ntime++) {
  delete v[ntime];
  delete a[ntime];
  delete uback[ntime];
  delete wt_back[ntime];
  delete q_obs[ntime];
  delete wt_obs[ntime];
 }

 return y;

}


Real algorithm_and_gradient(
 vector< MultiFab* > x,	// control
 vector< MultiFab* > dJdx, // gradient
 const Real* dx,
 Real dt,
 const Geometry& geom,
 int plot_int,
 int steepest_descent_iter) { 

 int nsteps=x.size()-1;

 adept::Stack stack;
 using adept::adouble;

 vector< My_adept_MFAB* > adept_x;  

 adept_x.resize(nsteps+1);

 adouble y;

 for (int ntime=0;ntime<=nsteps;ntime++) {
  BoxArray ba=x[ntime]->boxArray();
  DistributionMapping dm=x[ntime]->DistributionMap();
  int Ncomp=x[ntime]->nComp();
  int Nghost=x[ntime]->nGrow();

  adept_x[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);

  My_adept_MFAB* cur_frame_adept=adept_x[ntime];
  MultiFab* cur_frame=x[ntime];
  for (MFIter mfi(*cur_frame_adept,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   My_adept_FAB& adept_fab=(*cur_frame_adept)[gridno];
   FArrayBox& Real_fab=(*cur_frame)[gridno];
   Array4< adept::adouble > const& adept_array=adept_fab.array();
   Array4<Real> const& Real_array=Real_fab.array();
//   Dim3 lo = lbound(adept_array);
//   Dim3 hi = ubound(adept_array);
   Dim3 lo=lbound(adept_fab.box());
   Dim3 hi=ubound(adept_fab.box());
   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
    adept_array(i,j,k)=Real_array(i,j,k);
   } // i
   } // j
   } // k
  } // mfi
 } // ntime

 stack.new_recording();  // records the tape
 y = cost_function(adept_x,dx,dt,geom,plot_int,steepest_descent_iter);
 y.set_gradient(1.0);
 stack.compute_adjoint();
 for (int ntime=0;ntime<=nsteps;ntime++) {

  My_adept_MFAB* x_frame=adept_x[ntime];
  MultiFab* dJdx_frame=dJdx[ntime];

  for (MFIter mfi(*x_frame,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   My_adept_FAB& x_fab=(*x_frame)[gridno];
   FArrayBox& dJdx_fab=(*dJdx_frame)[gridno];
   Array4< adept::adouble > const& x_array=x_fab.array();
   Array4< double > const& dJdx_array=dJdx_fab.array();

//   Dim3 lo = lbound(x_array);
//   Dim3 hi = ubound(x_array);
   Dim3 lo=lbound(x_fab.box());
   Dim3 hi=ubound(x_fab.box());

   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
    dJdx_array(i,j,k)=x_array(i,j,k).get_gradient();
   } // i
   } // j
   } // k
  } // mfi
 }

 return y.value();
}

// history:
// Fortran was the language of choice for optimized simulation of fluid
// mechanics.
// Fortran is inconvenient for manipulating complex data structures.
// elementary 2d arrays in c++
// double* array2d=new double*[nrows];
// for (int i=0;i<ncols;i++) {
//  array2d[i]=new double[ncols];
//
// array2d[i][j]=...   0<=i<=nrows-1   0<=j<=ncols-1
//
// AMReX:  FarrayBox array2d(box);   box has dimensions of a 2d box.
// Vector  (1d)
// IntVect (Vector of integers)  Vector<int>
// Box: IntVect Boxlo, IntVect Boxhi   
// e.g. Boxlo(0)=0 Boxlo(1)=0
//      Boxhi(0)=nrows-1  Boxhi(1)=ncols-1
// BaseFab<double> = FArrayBox(Box box_in)   (array2d analogy)
// BaseFab<int>
// BaseFab<adept::adouble>
// Fab = Fortran Array Box and stores rows and columns in an order consistent
// with Fortran multidimensional arrays.
// BaseFab<double> my_double_array(box_in,....)
// FabArray<int>,
// FabArray<double>,
// FabArray<adept::adouble>
// FabArray is an array of BaseFab's i.e. Vector<BaseFab<adept::adouble>>
int main(int argc,char* argv[]) {

 amrex::Initialize(argc,argv);

 // What time is it now?  We'll use this to compute total run time.
 Real strt_time = amrex::second();

 // AMREX_SPACEDIM: number of dimensions
 int n_cell[AMREX_SPACEDIM];
 int max_grid_size[AMREX_SPACEDIM];
 int nsteps, plot_int;
 Real xlo[AMREX_SPACEDIM];
 Real xhi[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  xlo[dir]=0.0;
  xhi[dir]=1.0;
 }
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
  pp.get("n_cell",n_cell[0]);
  n_cell[1]=1;

  // The domain is broken into boxes of size max_grid_size
  pp.get("max_grid_size",max_grid_size[0]);

  // Default plot_int to -1, allow us to set it to 
  // something else in the inputs file
  //  If plot_int < 0 then no plot files will be writtenq
  plot_int = -1;
  pp.query("plot_int",plot_int);

  // Default nsteps to 10, allow us to set it 
  // to something else in the inputs file
  nsteps = 10;
  pp.query("nsteps",nsteps);

  pp.get("xlo",xlo[0]);
  pp.get("xhi",xhi[0]);
  xlo[1]=xlo[0];
  xhi[1]=xhi[0];
  pp.get("stop_time",stop_time);
 }


 // make BoxArray and Geometry
 BoxArray ba;
 Geometry geom;
 {
  IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
  IntVect dom_hi(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1));
  Box domain(dom_lo, dom_hi);

  // Initialize the boxarray "ba" from the single box "bx"
  ba.define(domain);
  // Break up boxarray "ba" into chunks no larger than 
  // "max_grid_size" along a direction
  ba.maxSize(max_grid_size[0]);

  RealBox real_box({AMREX_D_DECL(xlo[0],xlo[1],xlo[2])},
                   {AMREX_D_DECL(xhi[0],xhi[1],xhi[2])});

  // This defines a Geometry object
  geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
 }

 // Nghost = number of ghost cells for each array 
 int Nghost = 1;
    
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
	 x[ntime]=new MultiFab(ba,dm,Ncomp,Nghost);
	 dJdx[ntime]=new MultiFab(ba,dm,Ncomp,Nghost);

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

    // compute cost function and gradient of cost function 
   Real y=algorithm_and_gradient(
             x,dJdx,
	     dx,dt,geom,plot_int,iter);
 
   for (int ntime=0;ntime<=nsteps;ntime++) {
    MultiFab* cur_frame=x[ntime];
    MultiFab* cur_frame_dJdx=dJdx[ntime];
    for (MFIter mfi(*cur_frame,false); mfi.isValid(); ++mfi) {
     const int gridno = mfi.index();
     FArrayBox& x_fab=(*cur_frame)[gridno];
     FArrayBox& dJdx_fab=(*cur_frame_dJdx)[gridno];

     Array4<Real> const& x_array=x_fab.array();
     Array4<Real> const& dJdx_array=dJdx_fab.array();
//     Dim3 lo = lbound(x_array);
//     Dim3 hi = ubound(x_array);
     Dim3 lo=lbound(x_fab.box());
     Dim3 hi=ubound(x_fab.box());
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
