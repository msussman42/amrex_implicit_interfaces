#include <winstd.H>

#include <algorithm>

#if defined(BL_OLD_STL)
#include <stdlib.h>
#else
#include <cstdlib>
#endif

#include <ParmParse.H>
#include <Utility.H>
#include <ParallelDescriptor.H>
#include <CGSolver.H>
#include <MG_F.H>
#include <MultiGrid.H>

bool MultiGrid::initialized     = false;
int MultiGrid::def_nu_0         = 1;
int MultiGrid::def_nu_f         = 8;
int MultiGrid::def_maxiter      = 40;
int MultiGrid::def_maxiter_b    = 400;
int MultiGrid::def_numiter      = -1;
int MultiGrid::def_verbose      = 0;
int MultiGrid::def_nu_b         = 0;
int MultiGrid::def_numLevelsMAX = 1024;
int MultiGrid::def_use_bicgstab_in_mglib = 0;

static
void
Spacer (std::ostream& os, int lev)
{
 for (int k = 0; k < lev; k++) {
  os << "   ";
 }
}

void
MultiGrid::initialize () {

 ParmParse pp("mg");
 ParmParse ppLp("Lp");

 initialized = true;

 pp.query("maxiter", def_maxiter);
 pp.query("maxiter_b", def_maxiter_b);
 pp.query("numiter", def_numiter);
 pp.query("nu_0", def_nu_0);
 pp.query("nu_f", def_nu_f);
 pp.query("v", def_verbose);
 pp.query("verbose", def_verbose);
 pp.query("nu_b", def_nu_b);
 pp.query("numLevelsMAX", def_numLevelsMAX);

 ppLp.query("use_bicgstab_in_mglib",def_use_bicgstab_in_mglib);

 if (ParallelDescriptor::IOProcessor() && def_verbose) {
  std::cout << "MultiGrid settings...\n";
  std::cout << "   def_nu_0 =         " << def_nu_0         << '\n';
  std::cout << "   def_nu_f =         " << def_nu_f         << '\n';
  std::cout << "   def_maxiter =      " << def_maxiter      << '\n';
  std::cout << "   def_maxiter_b =      " << def_maxiter_b  << '\n';
  std::cout << "   def_nu_b =         " << def_nu_b         << '\n';
  std::cout << "   def_numLevelsMAX = " << def_numLevelsMAX << '\n';
  std::cout << "   def_use_bicgstab_in_mglib = " << 
	  def_use_bicgstab_in_mglib << '\n';
 }
}

MultiGrid::MultiGrid (ABecLaplacian &_Lp)
    :
    initialsolution(0),
    Lp(_Lp)
{
 if (!initialized)
     initialize();

 maxiter      = def_maxiter;
 numiter      = def_numiter;
 nu_0         = def_nu_0;
 nu_f         = def_nu_f;
 verbose      = def_verbose;
 maxiter_b   = def_maxiter_b;
 nu_b         = def_nu_b;
 numLevelsMAX = def_numLevelsMAX;
 numlevels    = numLevels();
 if (ParallelDescriptor::IOProcessor() && verbose > 2) {
  std::cout << "MultiGrid: " << numlevels
    << " multigrid levels created for this solve" << '\n';
  std::cout << "Grids: " << '\n';
  BoxArray tmp = Lp.boxArray(0);
  for (int i = 0; i < numlevels; ++i) {
   if (i > 0)
    tmp.coarsen(2);
   std::cout << " Level: " << i << '\n';
   for (int k = 0; k < tmp.size(); k++) {
    const Box& b = tmp[k];
    std::cout << "  [" << k << "]: " << b << "   ";
    for (int j = 0; j < BL_SPACEDIM; j++)
     std::cout << b.length(j) << ' ';
    std::cout << '\n';
   }
  }
 }
}

MultiGrid::~MultiGrid ()
{
 delete initialsolution;

 for (int i = 0; i < cor.size(); ++i) {
  delete res[i];
  delete rhs[i];
  delete cor[i];
 }
}

Real
MultiGrid::errorEstimate(int level,MultiFab& pbdry,Array<int> bcpres_array) {
  
 int nsolve=Lp.get_nsolve();

 Lp.residual(*(res[level]),*(rhs[level]),*(cor[level]), 
     level,pbdry,bcpres_array);
 MultiFab& resid = *(res[level]);
 int ncomp=resid.nComp();
 if (ncomp!=nsolve)
  BoxLib::Error("ncomp invalid");

 Real local_error=sqrt(Lp.norm(resid,level));
     
 return local_error;
}

void
MultiGrid::prepareForLevel (int level)
{
 if (cor.size() > level) return;

 int nghost1=1;
 int nghost0=0;

 res.resize(level+1, (MultiFab*)0);
 rhs.resize(level+1, (MultiFab*)0);
 cor.resize(level+1, (MultiFab*)0);

 Lp.prepareForLevel(level);

 if (cor[level] == 0) {
  int nsolve=Lp.get_nsolve();
  if (nsolve<=0)
   BoxLib::Error("nsolve invalid");
  res[level] = new MultiFab(Lp.boxArray(level),nsolve,nghost0,Fab_allocate);
  res[level]->setVal(0.0,0,nsolve,nghost0);
  rhs[level] = new MultiFab(Lp.boxArray(level),nsolve,nghost0,Fab_allocate);
  rhs[level]->setVal(0.0,0,nsolve,nghost0);
  cor[level] = new MultiFab(Lp.boxArray(level),nsolve,nghost1,Fab_allocate);
  cor[level]->setVal(0.0,0,nsolve,nghost1);
  if (level == 0) {
   initialsolution=new MultiFab(Lp.boxArray(level),nsolve,nghost1,Fab_allocate);
   initialsolution->setVal(0.0,0,nsolve,nghost1);
  }
 }
}

void MultiGrid::residualCorrectionForm (MultiFab& newrhs,
      MultiFab& oldrhs,MultiFab& solnL,
      MultiFab& inisol,MultiFab& pbdry,Array<int> bcpres_array,
      int level) {

 int nsolve=Lp.get_nsolve();
 if (solnL.nComp()!=nsolve)
  BoxLib::Error("ncomp invalid");
 if (solnL.nGrow()!=1)
  BoxLib::Error("solution should have ngrow=1");

 initialsolution->copy(inisol);
 solnL.copy(inisol);
 Lp.residual(newrhs, oldrhs, solnL, level, pbdry,bcpres_array);
 solnL.setVal(0.0,0,nsolve,1);
}

void
MultiGrid::solve (int nsverbose,
  MultiFab& _sol, MultiFab& _rhs,
  Real _eps_abs,Real _atol_b,
  int usecg_at_bottom,MultiFab& pbdry,
  Array<int> bcpres_array,
  int smooth_type,
  int bottom_smooth_type,int presmooth,
  int postsmooth) {

 int nsolve=Lp.get_nsolve();
 int level = 0;
 prepareForLevel(level);
 Lp.project_null_space(_rhs,level);
 residualCorrectionForm(*rhs[level],_rhs,*cor[level],
     _sol,pbdry,bcpres_array,level);
 Lp.project_null_space(*rhs[level],level);
 MultiFab* pbdryhom=new MultiFab(Lp.boxArray(level),nsolve,1,Fab_allocate);
 pbdryhom->setVal(0.0,0,nsolve,1);

 if (!solve_(nsverbose,_sol, _eps_abs, _atol_b, 
     *pbdryhom,bcpres_array,usecg_at_bottom,
     smooth_type,
     bottom_smooth_type,presmooth,postsmooth)) {
     BoxLib::Error("MultiGrid:: failed to converge!");
 }
 delete pbdryhom;
}

// pbdry will always be identically zero since residual correction form.
int
MultiGrid::solve_ (int nsverbose,MultiFab& _sol,
  Real eps_abs,Real atol_b,MultiFab& pbdry,Array<int> bcpres_array,
  int usecg_at_bottom,
  int smooth_type,int bottom_smooth_type,int presmooth,int postsmooth) {
  //
  // Relax system maxiter times, stop if 
  // absolute err <= _abs_eps
  //
  int nsolve=Lp.get_nsolve();
  int ncomp=_sol.nComp();
  if (ncomp!=nsolve)
   BoxLib::Error("ncomp invalid");

  int  level=0;
  int  returnVal = 0;
  const Real error0 = errorEstimate(level,pbdry,bcpres_array);
  Real error = error0;
  if (ParallelDescriptor::IOProcessor() && ((verbose)||(nsverbose>0))) {
   Spacer(std::cout, level);
   std::cout << "MultiGrid: Initial error (error0) = " << error0 << '\n';
  }

  //
  // Initialize correction to zero at this level (auto-filled at levels below)
  //
  (*cor[level]).setVal(0.0);
  int nit = 0;

  Real relative_error=1.0e-12;
  int error_close_to_zero=((error<eps_abs)||
                           (error<=relative_error*error0));

  for(nit = 0;
          	 (nit < maxiter)
	    &&   (nit < numiter || numiter < 0)
            &&   ((!error_close_to_zero)||((numiter>0)&&(numiter<10)) );
	++nit) {

   relax(*cor[level],*rhs[level],level,eps_abs,
     atol_b,usecg_at_bottom,pbdry,bcpres_array,
     smooth_type,bottom_smooth_type,
     presmooth,postsmooth);
   Lp.project_null_space(*cor[level],level);

   error = errorEstimate(level,pbdry,bcpres_array);
   error_close_to_zero=((error<eps_abs)||
                        (error<=relative_error*error0));

   if (ParallelDescriptor::IOProcessor()) {
       if (verbose > 1 ) {
           Spacer(std::cout, level);
           std::cout << "MultiGrid: Iteration "
                     << nit
                     << " error/error0 "
                     << error/error0 
                     << " error " 
                     << error << '\n';
       }
   }
  }
  if (ParallelDescriptor::IOProcessor()) {
   if (((verbose==1)||(nsverbose>0))&&(numiter<0)) {
    std::cout << "mg (linop-no precond) nit " << nit << '\n';
    std::cout << "final Iteration " << nit << " error0 "
      << error0 << " error " << error << '\n';
   }
  }

  if ((nit == numiter)||(error_close_to_zero)) {
   _sol.copy(*cor[level]);
   _sol.plus(*initialsolution,0,_sol.nComp(),0);
   Lp.project_null_space(_sol,level);
   returnVal = 1;
  }
  if (returnVal!=1) {
   std::cout << "nit= " << nit << " numiter= " << numiter << 
    " error= " << error << " relative_error = " << relative_error << 
    " error0=" << error0 << " eps_abs= " << eps_abs << '\n';
   BoxLib::Error("error in multigrid _solve");
  }

  //
  // Otherwise, failed to solve satisfactorily
  //
  return returnVal;
}

int
MultiGrid::numLevels () const
{
 int ng = Lp.numGrids();
 int lv = numLevelsMAX;
 //
 // The routine `falls through' since coarsening and refining
 // a unit box does not yield the initial box.
 //
 const BoxArray& bs = Lp.boxArray(0);

 for (int i = 0; i < ng; ++i) {
  int llv = 0;
  Box tmp = bs[i];
  for (;;) {
      Box ctmp  = tmp;   ctmp.coarsen(2);
      Box rctmp = ctmp; rctmp.refine(2);
      if (tmp != rctmp || ctmp.numPts() == 1)
          break;
      llv++;
      tmp = ctmp;
  }
  //
  // Set number of levels so that every box can be refined to there.
  //
  if (lv >= llv)
      lv = llv;
 }

 return lv+1; // Including coarsest.
}

void
MultiGrid::coarsestSmooth(MultiFab& solL,MultiFab& rhsL,
   int level,Real eps_abs,Real atol_b,int usecg_at_bottom,
   MultiFab& pbdry,Array<int> bcpres_array,
   int smooth_type,int bottom_smooth_type,
   int presmooth,int postsmooth)
{

 prepareForLevel(level);

 int nsolve=Lp.get_nsolve();
 int ncomp=solL.nComp();
 if (ncomp!=nsolve)
  BoxLib::Error("ncomp invalid");

 int use_mg_precond = 0;
 int is_bottom=1;

 if (usecg_at_bottom==0) {
  Real error0;
  if (verbose) {
   error0 = errorEstimate(level,pbdry,bcpres_array);
   if (ParallelDescriptor::IOProcessor())
    std::cout << "   Bottom Smoother: Initial error (error0) = " 
       << error0 << '\n';
  }

  Lp.project_null_space(rhsL,level);

  for (int i = finalSmooth(); i > 0; i--) {
   Lp.smooth(solL,rhsL,level,pbdry,bcpres_array,smooth_type); 

   if (verbose > 1 || (i == 1 && verbose)) {
    Real error = errorEstimate(level,pbdry,bcpres_array);
    if (ParallelDescriptor::IOProcessor())
     std::cout << "   Bottom Smoother: Iteration " << i
       << " error/error0 " << error/error0 << " error " 
       << error << '\n';
   }
  }
 } else {
  CGSolver cg(Lp, use_mg_precond, level);
  cg.setMaxIter(maxiter_b);
  int local_meets_tol=0;
  Real local_error0=0.0;
  int nsverbose=0;

  cg.solve(
    def_use_bicgstab_in_mglib,
    nsverbose,is_bottom,
    solL,rhsL, atol_b, atol_b,
    pbdry,bcpres_array,usecg_at_bottom,
    local_meets_tol,
    bottom_smooth_type,bottom_smooth_type,
    presmooth,postsmooth,local_error0);
 }
}


void
MultiGrid::relax (MultiFab& solL,MultiFab& rhsL,
   int level,Real eps_abs,
   Real atol_b,int usecg_at_bottom,
   MultiFab& pbdry,Array<int> bcpres_array,
   int smooth_type,int bottom_smooth_type,int presmooth,
   int postsmooth) {

 int nsolve=Lp.get_nsolve();
 int ncomp=solL.nComp();
 if (ncomp!=nsolve)
  BoxLib::Error("ncomp invalid");

 if (level < numlevels - 1 ) {
 
  if (presmooth!=postsmooth)
   BoxLib::Error("presmooth must equal postsmooth for mgpcg");

  Lp.project_null_space(rhsL,level);

  for (int i = presmooth ; i > 0 ; i--) {
   Lp.smooth(solL,rhsL,level,pbdry,bcpres_array,smooth_type);
  }
  Lp.residual(*res[level],rhsL,solL,level,pbdry,bcpres_array);
  Lp.project_null_space(*res[level],level);
  prepareForLevel(level+1);
  average(*rhs[level+1], *res[level],level+1,level);
  cor[level+1]->setVal(0.0);

  if (!((usecg_at_bottom==0)||(usecg_at_bottom==1)))
   BoxLib::Error("usecg_at_bottom invalid");

  MultiFab* pbdrycoarser=new MultiFab(Lp.boxArray(level+1),
    nsolve,1,Fab_allocate);
  pbdrycoarser->setVal(0.0,0,nsolve,1); 
  for (int i = cntRelax(); i > 0 ; i--) {
   relax(*(cor[level+1]),*(rhs[level+1]),level+1,
    eps_abs,atol_b,usecg_at_bottom,
    *pbdrycoarser,bcpres_array,
    smooth_type,bottom_smooth_type,presmooth,postsmooth);
  }
  delete pbdrycoarser;

  interpolate(solL, *(cor[level+1]),level+1,level);
  for (int i = postsmooth; i > 0 ; i--) {
   Lp.smooth(solL, rhsL, level,pbdry,bcpres_array,smooth_type);
  }
 } else {
  coarsestSmooth(solL,rhsL,level,eps_abs,atol_b,
   usecg_at_bottom,pbdry,bcpres_array,
   smooth_type,bottom_smooth_type,
   presmooth,postsmooth);
 }
}


void
MultiGrid::average (MultiFab& c,MultiFab& f,
  int clevel,int flevel)
{
 if (clevel!=flevel+1)
  BoxLib::Error("clevel invalid");
 if (flevel<0)
  BoxLib::Error("flevel invalid"); 

 int bfact_coarse=Lp.get_bfact_array(clevel);
 int bfact_fine=Lp.get_bfact_array(flevel);
 int bfact_top=Lp.get_bfact_array(0);

#ifdef _OPENMP
#pragma omp parallel
#endif
 for (MFIter mfi(c); mfi.isValid(); ++mfi) {

  BL_ASSERT(c.boxArray().get(mfi.index()) == mfi.validbox());

  const Box& bx = mfi.validbox();

  int nsolve=Lp.get_nsolve();
  int nc = c.nComp();
  if (nc!=nsolve) {
   std::cout << "nc,nsolve = " << nc << ' ' << nsolve << '\n';
   BoxLib::Error("nc invalid in average");
  }
   // divide by 4 in 2D and 8 in 3D
  int iaverage=1;
  for (int veldir=0;veldir<nsolve;veldir++) {
   FORT_AVERAGE(
    c[mfi].dataPtr(veldir),
    ARLIM(c[mfi].loVect()), ARLIM(c[mfi].hiVect()),
    f[mfi].dataPtr(veldir),
    ARLIM(f[mfi].loVect()), ARLIM(f[mfi].hiVect()),
    bx.loVect(), bx.hiVect(),&iaverage,
    &bfact_coarse,&bfact_fine,&bfact_top);
  }  // veldir
 } // mfi
 ParallelDescriptor::Barrier();
}

void
MultiGrid::interpolate (MultiFab& f,MultiFab& c,
  int clevel,int flevel)
{
 if (clevel!=flevel+1)
  BoxLib::Error("clevel invalid");
 if (flevel<0)
  BoxLib::Error("flevel invalid"); 

 int bfact_coarse=Lp.get_bfact_array(clevel);
 int bfact_fine=Lp.get_bfact_array(flevel);
 int bfact_top=Lp.get_bfact_array(0);


 //
 // Use fortran function to interpolate up (prolong) c to f
 // Note: returns f=f+P(c) , i.e. ADDS interp'd c to f.
 //
#ifdef _OPENMP
#pragma omp parallel
#endif
 for (MFIter mfi(f); mfi.isValid(); ++mfi) {

  const Box& bx = c.boxArray()[mfi.index()];
  int nsolve=Lp.get_nsolve();
  int nc = f.nComp();
  if (nc!=nsolve) {
   std::cout << "nc,nsolve = " << nc << ' ' << nsolve << '\n';
   BoxLib::Error("nc invalid in interpolate");
  }

  for (int veldir=0;veldir<nsolve;veldir++) {
   FORT_INTERP(
     &bfact_coarse,&bfact_fine,&bfact_top,
     f[mfi].dataPtr(veldir),
     ARLIM(f[mfi].loVect()), ARLIM(f[mfi].hiVect()),
     c[mfi].dataPtr(veldir),
     ARLIM(c[mfi].loVect()), ARLIM(c[mfi].hiVect()),
     bx.loVect(), bx.hiVect());
  } // veldir

 } // mfi
 ParallelDescriptor::Barrier();
}
