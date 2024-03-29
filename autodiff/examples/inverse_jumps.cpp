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
	if (eps==0.0) {
		if (x>=0.0)
			return 1.0;
		else if (x<0.0)
			return 0.0;
		else
			amrex::Error("x invalid");
	} else if (eps>0.0) {
         // "Logistic function"
         Hreturn=1.0/(1.0+exp(-x/eps));
	} else
		amrex::Error("eps invalid");

	return Hreturn;
}

// wave equation:
//   q_tt = c^2 q_xx   c=speed of wave propagation
//   if q(0,x)=f(x) and q_{t}(0,x)=0
//   the exact solution to the wave equation is:
//   q(t,x)=(1/2)(f(x-ct) + f(x+ct))
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

void plot_adept_data(
	std::string plt_name_string,
	int steepest_descent_iter,
	int ntime,
	const Geometry& geom,
	My_adept_MFAB* data_to_plot,double time) {

 const std::string& pltfile1 = 
	 amrex::Concatenate(plt_name_string,steepest_descent_iter,5);
 const std::string& pltfile = 
	 amrex::Concatenate(pltfile1,ntime,5);
 BoxArray ba_cell=data_to_plot->boxArray();
 DistributionMapping dm=data_to_plot->DistributionMap();
 int Ncomp=data_to_plot->nComp();
 IntVect Nghost=data_to_plot->nGrowVect();
 MultiFab* data_plot_mf=new MultiFab(ba_cell,dm,Ncomp,Nghost);

 for (MFIter mfi(*data_plot_mf,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   My_adept_FAB& v_fab=(*data_to_plot)[gridno]; // type: adouble
   FArrayBox& plot_fab=(*data_plot_mf)[gridno]; // type: double

   Array4< adept::adouble > const& v_array=v_fab.array(); // type: adouble

   Array4< Real > const& plot_array=plot_fab.array(); // type: double
   Dim3 lo=lbound(plot_fab.box());
   Dim3 hi=ubound(plot_fab.box());
   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
   for (int n = 0; n<Ncomp; ++n) {
    plot_array(i,j,k,n)=v_array(i,j,k,n).value();
   }
   }
   }
   }
 }

 std::cout << "plotting " << pltfile << '\n';
 int nghost_cell=0; // max and min not upgraded yet to take IntVect nghost
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

 int num_state_vars=1; // the number of state variables for each material
 int num_materials=2;
 int num_interfaces=1;

 adept::adouble local_pi=4.0*atan(1.0);

// recommended eps from MF De Pando, D Sipp, PJ Schmid
// Journal of Computational Physics 231 (23), 7739-7755
 adept::adouble eps=sqrt(1.0e-15);
// eps=0.0;

 vector< My_adept_MFAB* > v;  // approximate solution
 vector< My_adept_MFAB* > uback;
 vector< My_adept_MFAB* > wt_back;
 vector< My_adept_MFAB* > q_obs;
 vector< My_adept_MFAB* > wt_obs;

 v.resize(nsteps+1);
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
 IntVect Nghost=uinput[ntime]->nGrowVect();

 adept::adouble xlo=geom.ProbLo(0);
 adept::adouble xhi=geom.ProbHi(1);;

 adept::adouble stop_time=dt*nsteps;

 for (ntime=0;ntime<=nsteps;ntime++) {

  int Ncomp=uinput[ntime]->nComp();
  adept::adouble local_t=ntime*dt;

  v[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  uback[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  wt_back[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  q_obs[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);
  wt_obs[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);

  My_adept_MFAB* uinput_frame=uinput[ntime];

  My_adept_MFAB* v_frame=v[ntime];
  My_adept_MFAB* uback_frame=uback[ntime];
  My_adept_MFAB* wt_back_frame=wt_back[ntime];
  My_adept_MFAB* q_obs_frame=q_obs[ntime];
  My_adept_MFAB* wt_obs_frame=wt_obs[ntime];

  for (MFIter mfi(*v_frame,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();

   My_adept_FAB& uinput_fab=(*uinput_frame)[gridno];

   My_adept_FAB& v_fab=(*v_frame)[gridno];
   My_adept_FAB& uback_fab=(*uback_frame)[gridno];
   My_adept_FAB& wt_back_fab=(*wt_back_frame)[gridno];
   My_adept_FAB& q_obs_fab=(*q_obs_frame)[gridno];
   My_adept_FAB& wt_obs_fab=(*wt_obs_frame)[gridno];

   Array4< adept::adouble > const& uinput_array=uinput_fab.array();

   Array4< adept::adouble > const& v_array=v_fab.array();
   Array4< adept::adouble > const& uback_array=
	   uback_fab.array();
   Array4< adept::adouble > const& wt_back_array=
	   wt_back_fab.array();
   Array4< adept::adouble > const& q_obs_array=
	   q_obs_fab.array();
   Array4< adept::adouble > const& wt_obs_array=
	   wt_obs_fab.array();

   Dim3 lo=lbound(v_fab.box());
   Dim3 hi=ubound(v_fab.box());

   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
   for (int n = 0; n<Ncomp; ++n) {

    adept::adouble local_x=xlo+dx[0]*(i+0.5);

    adept::adouble wt_back=0.1;
    if (n==Ncomp-1)
     wt_back=0.0;

     // background solution
    uback_array(i,j,k,n)=0.0;
     // background error weight
    wt_back_array(i,j,k,n)=0.0;
    if (ntime==0) {
     wt_back_array(i,j,k,n)=wt_back;
    }
    if (local_x<=xlo) {
     wt_back_array(i,j,k,n)=wt_back;
    }
    if (local_x>=xhi) {
     wt_back_array(i,j,k,n)=wt_back;
    }
     // observation data only considered at t=Tstop
    wt_obs_array(i,j,k,n)=0.0;
    q_obs_array(i,j,k,n)=0.0;
    if (ntime==nsteps) {
     if (n==0) {
      q_obs_array(i,j,k,n)=1.0; //state1
     } else if (n==1) {
      q_obs_array(i,j,k,n)=0.0; //state2
     } else if (n==2) {  
      q_obs_array(i,j,k,n)=0.25-fabs(x-0.5); //phi (levelset function)
     } else
      amrex::Error("n invalid");

     wt_obs_array(i,j,k,n)=1.0;
    }
    v_array(i,j,k,n)=uinput_array(i,j,k,n);
   } // n
   } // i
   } // j
   } // k
  } // mfi
 } // ntime


 adept::adouble y=0.0;

 double time=0.0;

 for (ntime=0;ntime<nsteps;ntime++) {
  int Ncomp=uinput[ntime]->nComp();

  My_adept_MFAB* v_frame=v[ntime];
  My_adept_MFAB* vnp1_frame=v[ntime+1];
  My_adept_MFAB* unp1_frame=uinput[ntime+1];

  My_adept_MFAB* wtnp1_back_frame=wt_back[ntime+1];
  My_adept_MFAB* unp1_back_frame=uback[ntime+1];
  My_adept_MFAB* qnp1_obs_frame=q_obs[ntime+1];
  My_adept_MFAB* wtnp1_obs_frame=wt_obs[ntime+1];

   // finite volume method, unknowns located at cell centers.
  My_adept_MFAB* fluxes=new My_adept_MFAB(ba_flux_x,dm,Ncomp,0);

  if (plot_int>0) {
   int ntime_div=ntime/plot_int;
   if (ntime_div*plot_int==ntime) {
    plot_adept_data("plt",steepest_descent_iter,
		    ntime,geom,v_frame,time);
   }
  }

  for (MFIter mfi(*H_cells,false); mfi.isValid(); ++mfi) {

   const int gridno = mfi.index();

   My_adept_FAB& v_fab=(*v_frame)[gridno];
   My_adept_FAB& vnp1_fab=(*vnp1_frame)[gridno];
   My_adept_FAB& unp1_fab=(*unp1_frame)[gridno];

   My_adept_FAB& wtnp1_back_fab=(*wtnp1_back_frame)[gridno];
   My_adept_FAB& unp1_back_fab=(*unp1_back_frame)[gridno];
   My_adept_FAB& qnp1_obs_fab=(*qnp1_obs_frame)[gridno];
   My_adept_FAB& wtnp1_obs_fab=(*wtnp1_obs_frame)[gridno];

   My_adept_FAB& fluxes_fab=(*fluxes)[gridno];

   Array4< adept::adouble > const& v_array=v_fab.array();
   Array4< adept::adouble > const& vnp1_array=vnp1_fab.array();
   Array4< adept::adouble > const& unp1_array=unp1_fab.array();

   Array4< adept::adouble > const& wtnp1_back_array=wtnp1_back_fab.array();
   Array4< adept::adouble > const& unp1_back_array=unp1_back_fab.array();
   Array4< adept::adouble > const& qnp1_obs_array=qnp1_obs_fab.array();
   Array4< adept::adouble > const& wtnp1_obs_array=wtnp1_obs_fab.array();

   Array4< adept::adouble > const& fluxes_array=fluxes_fab.array();

   Dim3 lo=lbound(v_fab.box());
   Dim3 hi=ubound(v_fab.box());

   Dim3 lo_x=lbound(fluxes_fab.box());
   Dim3 hi_x=ubound(fluxes_fab.box());

    // cells: x_i=(i+1/2)*h  i=0,..,n_cells-1
    // flux locations: x_{i-1/2}=i * h  i=0 .. n_cells
   for (int k = lo_x.z; k <= hi_x.z; ++k) {
   for (int j = lo_x.y; j <= hi_x.y; ++j) {
   for (int i = lo_x.x; i <= hi_x.x; ++i) {
   for (int n = 0; n<Ncomp; ++n) {

    adept::adouble a=1.0;
    adept::adouble H=H_smooth(a,eps);
       // cell(i) lives at xi=(i+1/2)dx
       // face(i) lives at i*dx
    fluxes_array(i,j,k,n)=
        a*( (1.0-H)*v_array(i,j,k,n)+H*v_array(i-1,j,k,n) )
   } //n
   } //i
   } //j
   } //k

   Dim3 lo_int=lbound(ba[gridno]);
   Dim3 hi_int=ubound(ba[gridno]);

   for (int k = lo_int.z; k <= hi_int.z; ++k) {
   for (int j = lo_int.y; j <= hi_int.y; ++j) {
   for (int i = lo_int.x; i <= hi_int.x; ++i) {
   for (int n = 0; n<Ncomp; ++n) {
     vnp1_array(i,j,k)=v_array(i,j,k)-
	     (dt/dx[0])*(fluxes_array(i+1,j,k)-fluxes_array(i,j,k));
   } //n
   } //i
   } //j
   } //k

   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
   for (int n = 0; n<Ncomp; ++n) {
    // ok to have this "if" statement here since the 
    // contents of the conditional rule does not depend on
    // the control "u" or the state "v"
    if ((i<lo_int.x)||(i>hi_int.x)||
        (j<lo_int.y)||(j>hi_int.y)||
        (k<lo_int.z)||(k<hi_int.z)) {
     int iface=i;
     int jface=j;
     int kface=k;
     int i_interior=i;
     int j_interior=j;
     int k_interior=k;
      // face_normal points out of the domain
      // if a * face_normal < 0 => Dirichlet boundary condition
      //   wave moves from outside the domain into the inside.
      // if a * face_normal >0  => numerical boundary condition 
      //   wave moves from inside the domain  to the outside
     double face_normal=0.0;
     if (i<lo_int.x) {
      iface=lo_int.x;
      i_interior=lo_int.x;
      face_normal=-1.0;
     } else if (i>hi_int.x) {
      iface=hi_int.x+1;
      i_interior=hi_int.x;
      face_normal=1.0;
     }
     if (j<lo_int.y) {
      jface=lo_int.y;
      j_interior=lo_int.y;
      face_normal=-1.0;
     } else if (j>hi_int.y) {
      jface=hi_int.y+1;
      j_interior=hi_int.y;
      face_normal=1.0;
     }
     if (k<lo_int.z) {
      kface=lo_int.z;
      k_interior=lo_int.z;
      face_normal=-1.0;
     } else if (k>hi_int.z) {
      kface=hi_int.z+1;
      k_interior=hi_int.z;
      face_normal=1.0;
     }
     adept::adouble a=1.0;
      // a_bias < 0 => Dirichlet boundary
     adept::adouble a_bias=a*face_normal;
      // H_bias = 0 => Dirichlet boundary
     adept::adouble H_bias=H_smooth(a_bias,eps);
     
     vnp1_array(i,j,k,n)=(1.0-H_bias)*unp1_array(i,j,k,n)+
       H_bias*vnp1_array(i_interior,j_interior,k_interior,n);

    }
   }
   }
   }
   }

     // next time, levelset redistancing and state variable extrapolation
     // go here.
     // Note: for nonlinear problems (e.g. Stefan problem, or Burger's
     // equation), the speed for interface propagation must be extrapolated
     // too.
     //
   for (int k = lo_int.z; k <= hi_int.z; ++k) {
   for (int j = lo_int.y; j <= hi_int.y; ++j) {
   for (int i = lo_int.x; i <= hi_int.x; ++i) {

    double xi=(i+0.5)*dx[0];

    adept::adouble local_wt=0.0;
     // conditional does not depend on the control or state.
    if (wtnp1_back_array(i,j,k)>0.0)
     local_wt=1.0;

     // prescribe the control solution where background error weight > 0.0
    if (1==1) {
     vnp1_array(i,j,k)=
	 local_wt*unp1_array(i,j,k)+ 
	 (1.0-local_wt)*vnp1_array(i,j,k);
    }
    
    y=y+wtnp1_back_array(i,j,k)*dx[0]*
        (unp1_array(i,j,k)-unp1_back_array(i,j,k))*
        (unp1_array(i,j,k)-unp1_back_array(i,j,k));

    y=y+wtnp1_obs_array(i,j,k)*dx[0]*
        (qnp1_obs_array(i,j,k)-vnp1_array(i,j,k))*
        (qnp1_obs_array(i,j,k)-vnp1_array(i,j,k));
   } // i
   } // j
   } // k

  } // mfi

  delete fluxes;
  delete H_fluxes;
  delete a_fluxes;
  delete H_cells;

  time=time+dt;

  if (ntime==nsteps-1) {

   if (plot_int>0) {
    plot_adept_data("plt",steepest_descent_iter,
		    ntime,geom,vnp1_frame,time);
   }

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
  IntVect Nghost=x[ntime]->nGrowVect();

  adept_x[ntime]=new My_adept_MFAB(ba,dm,Ncomp,Nghost);

  My_adept_MFAB* cur_frame_adept=adept_x[ntime];
  MultiFab* cur_frame=x[ntime];
  for (MFIter mfi(*cur_frame_adept,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   My_adept_FAB& adept_fab=(*cur_frame_adept)[gridno];
   FArrayBox& Real_fab=(*cur_frame)[gridno];
   Array4< adept::adouble > const& adept_array=adept_fab.array();
   Array4<Real> const& Real_array=Real_fab.array();
   Dim3 lo=lbound(adept_fab.box());
   Dim3 hi=ubound(adept_fab.box());
   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
   for (int n = 0; n<Ncomp; ++n) {
    adept_array(i,j,k,n)=Real_array(i,j,k,n);
   } // n
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
  int Ncomp=x_frame[ntime]->nComp();

  for (MFIter mfi(*x_frame,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   My_adept_FAB& x_fab=(*x_frame)[gridno];
   FArrayBox& dJdx_fab=(*dJdx_frame)[gridno];
   Array4< adept::adouble > const& x_array=x_fab.array();
   Array4< double > const& dJdx_array=dJdx_fab.array();

   Dim3 lo=lbound(x_fab.box());
   Dim3 hi=ubound(x_fab.box());

   for (int k = lo.z; k <= hi.z; ++k) {
   for (int j = lo.y; j <= hi.y; ++j) {
   for (int i = lo.x; i <= hi.x; ++i) {
   for (int n = 0; n<Ncomp; ++n) {
    dJdx_array(i,j,k,n)=x_array(i,j,k,n).get_gradient();
   } // n
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

 int num_state_vars=1; // the number of state variables for each material
 int num_materials=2;
 int num_interfaces=1;

 // DEFAULT: NONE PERIODIC 
 // is_periodic(dir)=0 in all direction by default (dir=0,1,.., sdim-1)
 Vector<int> is_periodic(AMREX_SPACEDIM,0);  
     
 int max_opt_steps=10;
 Real learning_rate=0.01;

 // inputs parameters
 {
  // ParmParse is way of reading inputs from the inputs file
  ParmParse pp;

  // We need to get n_cell from the inputs file - this is the 
  // number of cells on each side of 
  //   a square (or cubic) domain.
  pp.get("n_cell",n_cell[0]);
  n_cell[1]=1;
  n_cell[AMREX_SPACEDIM-1]=1;

  // The domain is broken into boxes of size max_grid_size
  pp.get("max_grid_size",max_grid_size[0]);
  max_grid_size[1]=max_grid_size[0];
  max_grid_size[AMREX_SPACEDIM-1]=max_grid_size[0];

  // Default plot_int to -1, allow us to set it to 
  // something else in the inputs file
  //  If plot_int < 0 then no plot files will be writtenq
  plot_int = -1;
  pp.query("plot_int",plot_int);

  // Default nsteps to 10, allow us to set it 
  // to something else in the inputs file
  nsteps = 10;
  pp.query("nsteps",nsteps);
  pp.query("max_opt_steps",max_opt_steps);
  pp.query("learning_rate",learning_rate);

  pp.get("xlo",xlo[0]);
  pp.get("xhi",xhi[0]);
  xlo[1]=xlo[0];
  xhi[1]=xhi[0];
  xlo[AMREX_SPACEDIM-1]=xlo[0];
  xhi[AMREX_SPACEDIM-1]=xhi[0];
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

 IntVect nGrowVect(AMREX_SPACEDIM);
 nGrowVect[0]=1;
 nGrowVect[1]=0;
 nGrowVect[AMREX_SPACEDIM-1]=0;
    
 // Ncomp = number of components for each array
 // state1  state2   Levelset function
 // if levelset(t,x)<0 => state2
 // if levelset(t,x)>0 => state1
 //
 int Ncomp  = num_materials*num_state_vars+num_interfaces;
  
 // How Boxes are distrubuted among MPI processes
 DistributionMapping dm(ba);

 const Real* dx = geom.CellSize();
 Real dt=stop_time/nsteps;

 vector< MultiFab* > x;  // control
 vector< MultiFab* > dJdx; 
 x.resize(nsteps+1);
 dJdx.resize(nsteps+1);

// BaseFab< adept::adouble > test_basefab();
// Array4< adept::adouble > const& a=test_basefab.array();


 for (int ntime=0;ntime<=nsteps;ntime++) {
    x[ntime]=new MultiFab(ba,dm,Ncomp,nGrowVect);
    dJdx[ntime]=new MultiFab(ba,dm,Ncomp,nGrowVect);

    x[ntime]->setVal(0.0);

    int comp=num_materials*num_state_vars;
    Real smallest_distance=-(xhi[0]-xlo[0]);
    x[ntime]->setVal(smallest_distance,comp,1,nGrowVect);

    dJdx[ntime]->setVal(0.0);

    MultiFab* cur_frame=x[ntime];
    for (MFIter mfi(*cur_frame,false); mfi.isValid(); ++mfi) {
     const int gridno = mfi.index();
     FArrayBox& x_fab=(*cur_frame)[gridno];
     Array4<Real> const& x_array=x_fab.array();
     Dim3 lo=lbound(x_fab.box());
     Dim3 hi=ubound(x_fab.box());
     for (int k = lo.z; k <= hi.z; ++k) {
     for (int j = lo.y; j <= hi.y; ++j) {
     for (int i = lo.x; i <= hi.x; ++i) {
      double xi=xlo[0]+(i+0.5)*dx[0];
       // initially, the control variable state1(t,x)=0 at t=0
       //                                 state1(t,x)=0 at x=0,1
       //                                 state2(t,x)=0 at t=0
       //                                 state2(t,x)=0 at x=0,1
       //                                 cost function:
       //                                 a) vi_t + a vi_x = 0
       //                                   vi(t,x)=ui(t,x) on the boundaries
       //                                   and initial condition
       //                                   vi(T,x)=ui(0,x-aT)=0  if x-aT>0
       //                                         =0  if x-aT<0
       //
       //                                         i=1,2
       //                                    extension PDEs:
       //                                    n points from 2 to 1
       //                                    n = grad phi/|grad phi|
       //       Peng, Merriman, Osher, Kang, Zhou equation (21):
       //                             characteristics point from 1 to 2:
       //                             v1(t,x)_tau - n H(-phi) dot grad v1 =0
       //                             characteristics point from 2 to 1:
       //                             v2(t,x)_tau + n H(phi) dot grad v2 =0
       //                                 b) phi_t + a phi_x = 0
       //                                    redistance:
       //                                    Hamilton-Jacobi equation: 
       //                                    phi_tau + s(phi_tilde) *
       //                                      (1-|grad phi|)=0
       //                                    Sussman, Smereka, Osher 1994
       //                 Peng, Merriman, Osher, Kang, Zhou equation (29)
       //
       //                         c) C=sum_{i=1}^{num_mat} 
       //                            ||vi(T,.)-vi^{obs}(T,.)||_w1i+
       //                            ||ui(t,x)-ui^{back}(t,x)||_w2i +
       //                            sum_{j=1}^{num_interfaces}
       //                            ||phi_j(T,.)-phi_j^{obs}(T,.)||_w3j
       //
       //  v1_observation(T,x)=1  
       //  v2_observation(T,x)=0
       //  phi_observation(T,x)=1/4-|x-1/2|
       //
       //  background solution: u^background = 0
       //  w2 = 0.1 at t=0, x=0, x=1, w2=0.0 otherwise
       //  The background error in the cost function is
       //  "Tikhonov regularization?"

     }
     }
     }
    } // mfi
 } // ntime=0 ..nsteps

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

 for (int iter=0;iter<max_opt_steps;iter++) {

    // compute cost function and gradient of cost function 
   Real y=algorithm_and_gradient(
             x,dJdx,
	     dx,dt,geom,plot_int,iter);

   std::cout << "iter, cost function " << iter << ' ' <<
    y << '\n';
 
   for (int ntime=0;ntime<=nsteps;ntime++) {
    MultiFab* cur_frame=x[ntime];
    MultiFab* cur_frame_dJdx=dJdx[ntime];
    int Ncomp=x[ntime]->nComp();
    for (MFIter mfi(*cur_frame,false); mfi.isValid(); ++mfi) {
     const int gridno = mfi.index();
     FArrayBox& x_fab=(*cur_frame)[gridno];
     FArrayBox& dJdx_fab=(*cur_frame_dJdx)[gridno];

     Array4<Real> const& x_array=x_fab.array();
     Array4<Real> const& dJdx_array=dJdx_fab.array();
     Dim3 lo=lbound(x_fab.box());
     Dim3 hi=ubound(x_fab.box());
     for (int k = lo.z; k <= hi.z; ++k) {
     for (int j = lo.y; j <= hi.y; ++j) {
     for (int i = lo.x; i <= hi.x; ++i) {
     for (int n = 0; n<Ncomp; ++n) {
      x_array(i,j,k,n)=x_array(i,j,k,n)-learning_rate*dJdx_array(i,j,k,n);
     } // n
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
