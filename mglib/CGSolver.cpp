#include <winstd.H>

#include <algorithm>
#include <iomanip>

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <CG_F.H>
#include <CGSolver.H>
#include <MultiGrid.H>

int              CGSolver::initialized            = 0;
int              CGSolver::def_maxiter            = 40;
int              CGSolver::def_verbose            = 0;

void
CGSolver::initialize ()
{
    ParmParse pp("cg");

    pp.query("maxiter", def_maxiter);
    pp.query("v", def_verbose);
    pp.query("verbose", def_verbose);

    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
        std::cout << "CGSolver settings...\n";
	std::cout << "   def_maxiter            = " << def_maxiter << '\n';
    }
    
    initialized = 1;
}

CGSolver::CGSolver (ABecLaplacian& _Lp,
            int    _use_mg_precond,
            int    _lev)
    :
    Lp(_Lp),
    mg_precond(0),
    lev(_lev),
    use_mg_precond(_use_mg_precond)
{
    if (!initialized)
        initialize();
    maxiter = def_maxiter;
    verbose = def_verbose;
    set_mg_precond();
}

void
CGSolver::set_mg_precond ()
{
    if (mg_precond!=0)
     BoxLib::Error("mg_precond set somehow");
    if ((use_mg_precond<0)||(use_mg_precond>1))
     BoxLib::Error("use_mg_precond invalid");

    if (use_mg_precond>0) {
        if (lev!=0)
         BoxLib::Error("mg precond only on top level");
        mg_precond = new MultiGrid(Lp);
        mg_precond->setNumIter(1);
    }
}

CGSolver::~CGSolver ()
{
    if (use_mg_precond>0)
     delete mg_precond;
    else if (mg_precond!=0)
     BoxLib::Error("mg_precond set somehow");
}

// z=K^{-1}r
void
CGSolver::pcg_solve(MultiFab* z,MultiFab* r,
    Real eps_abs,Real bot_atol,
    MultiFab* pbdryhom,
    Array<int> bcpres_array,
    int usecg_at_bottom,
    int smooth_type,int bottom_smooth_type,
    int presmooth,int postsmooth,
    int use_PCG,
    int ncomp,int nghost,int nghostRHS) {

 Lp.project_null_space((*r),lev);

 // z=K^{-1} r
 z->setVal(0.0,0,ncomp,nghost);
 if (use_PCG==0) {
  MultiFab::Copy(*z,*r,0,0,ncomp,nghostRHS);
 } else if (use_PCG==1) {
  if (use_mg_precond==1) {
   mg_precond->
    solve(0,*z,*r,eps_abs,bot_atol,
      usecg_at_bottom,*pbdryhom,bcpres_array,
      smooth_type,bottom_smooth_type,
      presmooth,postsmooth);
  } else if (use_mg_precond==0) {
   for (int j=0;j<presmooth+postsmooth;j++) {
    Lp.smooth((*z),(*r),lev,*pbdryhom,bcpres_array,smooth_type);
   }
  } else
   BoxLib::Error("use_mg_precond invalid");
 } else
  BoxLib::Error("use_PCG invalid");

 Lp.project_null_space((*z),lev);

} // subroutine pcg_solve

void 
CGSolver::check_for_convergence(Real rnorm,Real rnorm0,Real eps_abs,
		Real relative_error,int nit,int& error_close_to_zero) {

	int critical_nit=30;
	Real critical_abs_tol=1.0e-14;
	Real critical_rel_tol=1.0e-14;
	if (critical_abs_tol>eps_abs)
		critical_abs_tol=eps_abs;
	if (critical_rel_tol>relative_error)
		critical_rel_tol=relative_error;
	if (nit>critical_nit) {
         error_close_to_zero=((rnorm<=eps_abs)||
                              (rnorm<=relative_error*rnorm0));
	} else if ((nit>=0)&&(nit<=critical_nit)) {
         error_close_to_zero=((rnorm<=critical_abs_tol)||
                              (rnorm<=critical_rel_tol*rnorm0));
	} else
 	 BoxLib::Error("nit invalid");

} // check_for_convergence

void 
CGSolver::dump_params(Real rnorm,Real rnorm0,
		Real eps_abs,Real relative_error,
                int is_bottom,Real bot_atol,
                int usecg_at_bottom,int smooth_type,
		int bottom_smooth_type,int presmooth,
		int postsmooth,MultiFab& mf1,
		MultiFab& mf2) {

 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "lev= " << lev << '\n';
  std::cout << "rnorm0,rnorm,eps_abs,relative_error " <<
   rnorm0 << ' ' << rnorm << ' ' <<  eps_abs << 
   ' ' << relative_error << '\n';
  std::cout << "is_bottom= " << is_bottom << '\n';
  std::cout << "bot_atol= " << bot_atol << '\n';
  std::cout << "usecg_at_bottom= " << usecg_at_bottom << '\n';
  std::cout << "smooth_type= " << smooth_type << '\n';
  std::cout << "bottom_smooth_type= " << bottom_smooth_type << '\n';
  std::cout << "presmooth= " << presmooth << '\n';
  std::cout << "postsmooth= " << postsmooth << '\n';
  std::cout << "Lp.cfd_level= " << Lp.cfd_level << '\n';
  std::cout << "Lp.cfd_project_option= " << Lp.cfd_project_option << '\n';
  std::cout << "Lp.laplacian_solvability= " << 
          Lp.laplacian_solvability << '\n';
  std::cout << "Lp.check_for_singular= " << 
          Lp.check_for_singular << '\n';
  std::cout << "Lp.get_nsolve()= " << Lp.get_nsolve() << '\n';
  std::cout << "Lp.numGrids()= " << Lp.numGrids() << '\n';
  std::cout << "Lp.numLevels()= " << Lp.numLevels() << '\n';
  std::cout << "mf1.boxArray()= " << mf1.boxArray() << '\n';
  std::cout << "Lp.norm(mf1,lev)= " << Lp.norm(mf1,lev) << '\n';
  std::cout << "mf2.boxArray()= " << mf2.boxArray() << '\n';
  std::cout << "Lp.norm(mf2,lev)= " << Lp.norm(mf2,lev) << '\n';
 }

} // end subroutine dump_params

void
CGSolver::solve(
    int bicgstab_flag,
    int nsverbose,int is_bottom,
    MultiFab& sol,MultiFab& rhs,
    Real eps_abs,Real bot_atol,
    MultiFab& pbdry,
    Array<int> bcpres_array,
    int usecg_at_bottom,
    int& meets_tol,
    int smooth_type,int bottom_smooth_type,
    int presmooth,int postsmooth,
    Real& error0)
{
 //
 // algorithm:
 //
 //   k=0;r=rhs-A*soln_0;
 //   while (||r_k||^2_2 > eps^2*||r_o||^2_2 && k < maxiter {
 //      k++
 //      solve Mz_k-1 = r_k-1 (if preconditioning, else z_k-1 = r_k-1)
 //      rho_k-1 = r_k-1 dot z_k-1
 //      if (k=1) { p_1 = z_0 }
 //      else { beta = rho_k-1/rho_k-2; p = z + beta*p }
 //      Ap = A*p
 //      alpha = rho_k-1/(p dot Ap)
 //      x += alpha p
 //      r = b - A*x
 //   }
 //
 BL_ASSERT(sol.boxArray() == Lp.boxArray(lev));
 BL_ASSERT(rhs.boxArray() == Lp.boxArray(lev));
 BL_ASSERT(pbdry.boxArray() == Lp.boxArray(lev));

 Real relative_error=1.0e-12;

 int nghost = 1; 
 int nghostRHS = 0;

 int ngrow_v_search=0;
 if (bicgstab_flag==0) {
  ngrow_v_search=1;
 } else if (bicgstab_flag==1) {
  ngrow_v_search=0;
 } else
  BoxLib::Error("bicgstab_flag invalid");

 int nsolve=Lp.get_nsolve();
 int ncomp = sol.nComp();
 if (ncomp!=nsolve)
  BoxLib::Error("ncomp invalid");

 MultiFab* delta_sol = new MultiFab(sol.boxArray(), ncomp, nghost, 
   Fab_allocate);
 MultiFab* r  = new MultiFab(sol.boxArray(), ncomp, nghostRHS, 
   Fab_allocate);
 MultiFab* z  = new MultiFab(sol.boxArray(), ncomp, nghost, 
   Fab_allocate);
 MultiFab* Av_search = new MultiFab(sol.boxArray(), ncomp, nghostRHS, 
   Fab_allocate);
 MultiFab* p_search  = new MultiFab(sol.boxArray(), ncomp, nghostRHS, 
   Fab_allocate);
 MultiFab* v_search  = new MultiFab(sol.boxArray(), ncomp, ngrow_v_search, 
   Fab_allocate);
 MultiFab* rhs_resid_cor_form =
   new MultiFab(sol.boxArray(), ncomp, nghostRHS,Fab_allocate);
 MultiFab* pbdryhom  = new MultiFab(sol.boxArray(), ncomp, nghost, 
   Fab_allocate);
 pbdryhom->setVal(0.0,0,ncomp,1);

 Lp.project_null_space(rhs,lev);

 rhs_resid_cor_form->setVal(0.0,0,ncomp,nghostRHS); 
 MultiFab::Copy(*rhs_resid_cor_form,rhs,0,0,ncomp,nghostRHS);

 if ((verbose>0)||(nsverbose>0)||(1==0))
  if (ParallelDescriptor::IOProcessor())
   std::cout << "CGSolver: is_bottom= " << is_bottom << '\n';

  // resid,rhs,soln
 Lp.residual((*r),(*rhs_resid_cor_form),sol,lev,pbdry,bcpres_array);
 Lp.project_null_space((*r),lev);

  // put solution and residual in residual correction form
 delta_sol->setVal(0.0,0,ncomp,1);
 MultiFab::Copy(*rhs_resid_cor_form,*r,0,0,ncomp,nghostRHS);

 Real rnorm = sqrt(Lp.norm(*r, lev));
 Real rnorm0=rnorm;
 error0=rnorm0;

 if ((verbose>0)||(nsverbose>0)) {
  if (ParallelDescriptor::IOProcessor()) {
   if (is_bottom==1)
    std::cout << "CGsolver(BOTTOM):Initial error(error0)="<<rnorm << '\n';
   else
    std::cout << "CGsolver(NOBOT):Initial error(error0)="<<rnorm << '\n';
  }
 }

 int force_restart=0;
 int nit=0;
 int nit_since_restart=0;
 int restart_period;

 int use_PCG=1;
 if (use_PCG==0) 
  restart_period=100;
 else if (use_mg_precond==1)
  restart_period=10;
 else if (use_mg_precond==0)
  restart_period=50; 
 else
  BoxLib::Error("use_mg_precond invalid");

 if (z->nComp()!=nsolve)
  BoxLib::Error("ncomp invalid");

 Real beta=0.0;
 Real rho=1.0;
 Real rho_old=1.0;
 Real omega=1.0;
 Real alpha=1.0;
 p_search->setVal(0.0,0,ncomp,nghostRHS); 
 v_search->setVal(0.0,0,ncomp,ngrow_v_search); 

 Real restart_tol=eps_abs*eps_abs*1.0e-4;

 int error_close_to_zero=0;
 check_for_convergence(rnorm,rnorm0,eps_abs,relative_error,nit,
		 error_close_to_zero);

 if (ParallelDescriptor::IOProcessor()) {
  if (verbose>1) {
   if (is_bottom==1)
    std::cout << "CGSolver(BOT): rnorm0,eps_abs,relative_error " <<
     rnorm0 << ' ' << eps_abs << ' ' << relative_error << '\n';
   else
    std::cout << "CGSolver(NOBOT): rnorm0,eps_abs,relative_error " <<
     rnorm0 << ' ' << eps_abs << ' ' << relative_error << '\n';
  }
 }

 for(nit = 0;((nit < maxiter)&&(error_close_to_zero==0)) ; ++nit) {

  if (nit_since_restart>=restart_period) {
   nit_since_restart=0;
   beta=0.0;
   rho=1.0;
   rho_old=1.0;
   omega=1.0;
   alpha=1.0;
   p_search->setVal(0.0,0,ncomp,nghostRHS); 
   v_search->setVal(0.0,0,ncomp,ngrow_v_search); 
   sol.plus(*delta_sol,0,ncomp,0);
   Lp.project_null_space(sol,lev);
   delta_sol->setVal(0.0,0,ncomp,1);
   MultiFab::Copy(*rhs_resid_cor_form,*r,0,0,ncomp,0);
  }

  rho_old=rho;

  rnorm=sqrt(Lp.norm(*r,lev));
  if (nit==0)
   rnorm0=rnorm;

  check_for_convergence(rnorm,rnorm0,eps_abs,relative_error,nit,
		 error_close_to_zero);

  force_restart=0;

  if (error_close_to_zero==0) {

    // "meets_tol" informs the main solver in NavierStokes3.cpp 
    // whether it needs to continue.
   if ((nit==0)&&(rnorm>eps_abs*10.0)) {
    meets_tol=0;
   } else if ((nit==0)&&(rnorm<=eps_abs*10.0)) {
    // do nothing
   } else if (nit>0) {
    // do nothing
   } else
    BoxLib::Error("nit invalid");

   if (bicgstab_flag==0) {

    // z=K^{-1} r
    pcg_solve(z,r,eps_abs,bot_atol,pbdryhom,bcpres_array,
     usecg_at_bottom,smooth_type,bottom_smooth_type,
     presmooth,postsmooth,use_PCG,ncomp,nghost,nghostRHS);

    // rho=z dot r=z dot Kz>=0 if K SPD
    Lp.LP_dot(*z,*r,lev,rho); 

    if (rho>=0.0) {

     if (rho_old>restart_tol) {
      beta = rho/rho_old;
      // v_search=z+beta v_search
      advance( (*v_search), beta, (*z),(*v_search) );
     } else if (rho_old<=restart_tol) {
      force_restart=1;
     } else {
      std::cout << "rho_old= " << rho_old << '\n';
      dump_params(rnorm,rnorm0,eps_abs,relative_error,
       is_bottom,bot_atol,usecg_at_bottom,smooth_type,
       bottom_smooth_type,presmooth,postsmooth,
       *z,*r);
      BoxLib::Error("rho_old invalid");
     }

    } else {
     std::cout << "rho= " << rho << '\n';
     dump_params(rnorm,rnorm0,eps_abs,relative_error,
      is_bottom,bot_atol,usecg_at_bottom,smooth_type,
      bottom_smooth_type,presmooth,postsmooth,
      *z,*r);
     BoxLib::Error("rho invalid");
    }

    Real vAv=0.0;
    if (force_restart==0) {
     // Av_search=A*v_search 
     Lp.apply(*Av_search,*v_search,lev,*pbdryhom,bcpres_array);
     // vAv=(Av) dot v>=0.0 if A SPD
     Lp.LP_dot(*Av_search,*v_search,lev,vAv);

     if ((vAv>=0.0)&&(vAv<=restart_tol)) {
      force_restart=1;
     } else if (vAv>restart_tol) {
      // do nothing
     } else
      BoxLib::Error("vAv invalid");
     
     if (force_restart==0) {
      alpha = rho/vAv;
       // x=x+alpha v_search
      Lp.LP_update( (*delta_sol), alpha, (*delta_sol),(*v_search),lev );
      Lp.project_null_space((*delta_sol),lev);
      Lp.residual((*r),(*rhs_resid_cor_form),(*delta_sol),
	lev,*pbdryhom,bcpres_array); 
      Lp.project_null_space((*r),lev);
     } else if (force_restart==1) {
      // do nothing
     } else
      BoxLib::Error("force_restart invalid");
    } else if (force_restart==1) {
     // do nothing
    } else
     BoxLib::Error("force_restart invalid");

   } else if (bicgstab_flag==1) {

    Lp.LP_dot(*rhs_resid_cor_form,*r,lev,rho); 
    if ((rho_old>restart_tol)&&(omega>restart_tol)) {
     beta=rho*alpha/(rho_old*omega);
      // p=p - omega v
     advance( (*p_search),-omega,(*p_search),(*v_search) );
      // p=r + beta p
     advance( (*p_search),beta,(*r),(*p_search) );
      // z=K^{-1} p
     pcg_solve(z,p_search,eps_abs,bot_atol,pbdryhom,bcpres_array,
      usecg_at_bottom,smooth_type,bottom_smooth_type,
      presmooth,postsmooth,use_PCG,ncomp,nghost,nghostRHS);
      // v_search=A*z 
     Lp.apply(*v_search,*z,lev,*pbdryhom,bcpres_array);
    } else if ((rho_old<=restart_tol)||(omega<=restart_tol)) {
     force_restart=1;
    } else
     BoxLib::Error("rho_old or omega invalid");

    if (force_restart==0) {
     Lp.LP_dot(*rhs_resid_cor_form,*v_search,lev,alpha);

     if (alpha>restart_tol) {
      alpha=rho/alpha;

       // x=x+alpha z
      Lp.LP_update( (*delta_sol), alpha, (*delta_sol),(*z),lev );
      Lp.project_null_space((*delta_sol),lev);
      Lp.residual((*r),(*rhs_resid_cor_form),(*delta_sol),
	lev,*pbdryhom,bcpres_array); 
      Lp.project_null_space((*r),lev);
      rnorm=sqrt(Lp.norm(*r,lev));

      check_for_convergence(rnorm,rnorm0,eps_abs,relative_error,nit,
		 error_close_to_zero);
      if (error_close_to_zero==0) {
       // z=K^{-1} r
       pcg_solve(z,r,eps_abs,bot_atol,pbdryhom,bcpres_array,
        usecg_at_bottom,smooth_type,bottom_smooth_type,
        presmooth,postsmooth,use_PCG,ncomp,nghost,nghostRHS);
       // Av_search=A*z
       Lp.apply(*Av_search,*z,lev,*pbdryhom,bcpres_array);
       Real rAz=0.0;
       Real zAAz=0.0;
       // rAz=(Az) dot r =z^T A^T r = z^T A^T K z >=0 if A and K SPD.
       Lp.LP_dot(*Av_search,*r,lev,rAz);
       if (rAz>=0.0) {

        Lp.LP_dot(*Av_search,*Av_search,lev,zAAz);
	//Az dot Az >=0 if A SPD
        if (zAAz>restart_tol) {
         omega=rAz/zAAz;

         // x=x+omega z
         Lp.LP_update( (*delta_sol), omega, (*delta_sol),(*z),lev );
         Lp.project_null_space((*delta_sol),lev);
         Lp.residual((*r),(*rhs_resid_cor_form),(*delta_sol),
    	  lev,*pbdryhom,bcpres_array); 
         Lp.project_null_space((*r),lev);
        } else if ((zAAz>=0.0)&&(zAAz<=restart_tol)) {
         force_restart=1;
        } else
 	 BoxLib::Error("zAAz invalid");

       } else
        BoxLib::Error("rAz invalid");

      } else if (error_close_to_zero==1) {
       nit_since_restart++;
      } else
       BoxLib::Error("error_close_to_zero invalid");
     } else if (alpha<=restart_tol) {
      force_restart=1;
     } else
      BoxLib::Error("alpha invalid");
    } else if (force_restart==1) {
     // do nothing
    } else 
     BoxLib::Error("force_restart invalid");

   } else
    BoxLib::Error("bicgstab_flag invalid");

   if (force_restart==1) {

    if (error_close_to_zero==0) {
     nit_since_restart=0;
     beta=0.0;
     rho=1.0;
     rho_old=1.0;
     omega=1.0;
     alpha=1.0;
     p_search->setVal(0.0,0,ncomp,nghostRHS); 
     v_search->setVal(0.0,0,ncomp,ngrow_v_search); 
     sol.plus(*delta_sol,0,ncomp,0);
     Lp.project_null_space(sol,lev);
     delta_sol->setVal(0.0,0,ncomp,1);
     MultiFab::Copy(*rhs_resid_cor_form,*r,0,0,ncomp,0);
    } else if (error_close_to_zero==1) {
     BoxLib::Error("cannot have both force_restart and error_close_to_zero");
    } else
     BoxLib::Error("error_close_to_zero invalid");
   } else if (force_restart==0) {
    if (error_close_to_zero==0) {
     nit_since_restart++;
    } else if (error_close_to_zero==1) {
     // do nothing
    } else
     BoxLib::Error("error_close_to_zero invalid");
   } else 
    BoxLib::Error("force_restart invalid");

  } else if (error_close_to_zero==1) {
   nit_since_restart++;
  } else
   BoxLib::Error("error_close_to_zero invalid");

 }  // end of CGSolver loop

 sol.plus(*delta_sol,0,ncomp,0);
 Lp.project_null_space(sol,lev);

 if ((verbose>0)||(nsverbose>0)) {
  if (ParallelDescriptor::IOProcessor()) {
   if (use_mg_precond==1) {
    if (is_bottom==0)
     std::cout << "mgpcg (mac) nit (NOBOT)" << nit << '\n';
    else
     std::cout << "mgpcg (mac) nit (BOT)" << nit << '\n';
    std::cout << "CGSolver:" << " nit_since_restart " <<
     nit_since_restart << '\n';
   } else if (use_mg_precond==0) {
    if (is_bottom==0) {
     std::cout << "pcg (mac) nit (NOBOT)" << nit << '\n';
     std::cout << "CGSolver:" << " nit_since_restart " <<
      nit_since_restart << '\n';
    }
   } else
    BoxLib::Error("use_mg_precond invalid");
  }
 }

 if (!error_close_to_zero) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "Warning: CGSolver:: failed to converge! \n";
  }
  dump_params(rnorm,rnorm0,eps_abs,relative_error,
    is_bottom,bot_atol,usecg_at_bottom,smooth_type,
    bottom_smooth_type,presmooth,postsmooth,
    sol,rhs);
 }

 if (ncomp!=nsolve)
  BoxLib::Error("ncomp invalid");

 if ((verbose>0)||(nsverbose>0)) {
  Lp.residual((*r),rhs,sol,lev,pbdry,bcpres_array);
  Real testnorm = sqrt(Lp.norm(*r,lev));
  if (ParallelDescriptor::IOProcessor()) {
   if (is_bottom==1)
    std::cout << "residual non-homogeneous bc (BOT) " << testnorm << '\n';  
   else
    std::cout << "residual non-homogeneous bc (NOBOT)"<<testnorm << '\n';  
  }
 }

 delete pbdryhom;
 delete rhs_resid_cor_form;
 delete p_search;
 delete v_search;
 delete Av_search;
 delete z;
 delete r;
 delete delta_sol;
} // subroutine CGSolver::solve

void CGSolver::advance (
       MultiFab& p,
       Real beta, 
       const MultiFab& z,
       MultiFab& y) {
    //
    // Compute p = z  +  beta y
    // only interior cells are updated.
    //
    const BoxArray& gbox = Lp.boxArray(lev);
    int nsolve=Lp.get_nsolve();
    int ncomp = p.nComp();
    if (ncomp!=nsolve)
     BoxLib::Error("p ncomp invalid");
    if (z.nComp()!=nsolve)
     BoxLib::Error("z ncomp invalid");
    if (y.nComp()!=nsolve)
     BoxLib::Error("y ncomp invalid");
    if ((p.nGrow()!=0)&&(p.nGrow()!=1))
     BoxLib::Error("p ngrow invalid");
    if ((z.nGrow()!=0)&&(z.nGrow()!=1))
     BoxLib::Error("z ngrow invalid");
    if ((y.nGrow()!=0)&&(y.nGrow()!=1))
     BoxLib::Error("y ngrow invalid");

    const BoxArray& zbox = z.boxArray();

    int bfact=Lp.get_bfact_array(lev);
    int bfact_top=Lp.get_bfact_array(0);

    bool use_tiling=Lp.cfd_tiling;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(p,use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(zbox[mfi.index()] == gbox[mfi.index()]);
     int gridno=mfi.index();
     const Box& tilegrid=mfi.tilebox();
     const Box& fabgrid=gbox[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();

     for (int veldir=0;veldir<nsolve;veldir++) {
      FORT_CGADVCP(
       p[mfi].dataPtr(veldir),
       ARLIM(p[mfi].loVect()), ARLIM(p[mfi].hiVect()),
       z[mfi].dataPtr(veldir),
       ARLIM(z[mfi].loVect()), ARLIM(z[mfi].hiVect()),
       y[mfi].dataPtr(veldir),
       ARLIM(y[mfi].loVect()), ARLIM(y[mfi].hiVect()),
       &beta,
       tilelo,tilehi,
       fablo,fabhi,&bfact,&bfact_top);
     } // veldir
    } // mfi
} // omp
    ParallelDescriptor::Barrier();
}

