//#include <winstd.H>
#if defined(BL_OLD_STL)
#include <stdlib.h>
#else
#include <cstdlib>
#endif

#include <algorithm>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <ABecLaplacian.H>
#include <LO_F.H>
#include <ABec_F.H>
#include <CG_F.H>
#include <MG_F.H>

#define SCALAR_WORK_NCOMP 11

namespace amrex{

#define profile_solver 0

int ABecLaplacian::gmres_max_iter = 64;
int ABecLaplacian::gmres_precond_iter_base_mg = 4;
int ABecLaplacian::nghostRHS=0;
int ABecLaplacian::nghostSOLN=1;

Real ABecLaplacian::a_def     = 0.0;
Real ABecLaplacian::b_def     = 1.0;

int ABecLaplacian::CG_def_maxiter = 200;
int ABecLaplacian::CG_def_verbose = 0;

int ABecLaplacian::MG_def_nu_0         = 1;
int ABecLaplacian::MG_def_nu_f         = 8;
int ABecLaplacian::MG_def_verbose      = 0;
int ABecLaplacian::MG_def_nu_b         = 0;

extern void GMRES_MIN_CPP(Real** HH,Real beta, Real* yy,int m,
 int caller_id,int& status);

static
void
Spacer (std::ostream& os, int lev)
{
 for (int k = 0; k < lev; k++) {
  os << "   ";
 }
}

// level==0 is the finest level
void
ABecLaplacian::apply (MultiFab& out,MultiFab& in,
  int level,MultiFab& pbdry,Vector<int> bcpres_array) {

    applyBC(in,level,pbdry,bcpres_array);
    Fapply(out,in,level);
}


// level==0 is the finest level
void
ABecLaplacian::applyBC (MultiFab& inout,int level,
       MultiFab& pbdry,Vector<int> bcpres_array) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::applyBC";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);

 bprof.stop();
#endif

#if (profile_solver==1)
 bprof.start();
#endif

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

 if (inout.nGrow()!=1) {
  std::cout << "inout ngrow= " << inout.nGrow() << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "bfact_top= " << bfact_top << '\n';
  amrex::Error("inout ngrow<>1");
 }
 if (pbdry.nGrow()!=1)
  amrex::Error("pbdry ngrow<>1");
 
 if (inout.nComp()!=nsolve_bicgstab)
  amrex::Error("inout.nComp invalid");
 if (pbdry.nComp()!=nsolve_bicgstab)
  amrex::Error("pbdry.nComp invalid");

 inout.FillBoundary(geomarray[level].periodicity());

  //
  // Fill boundary cells.
  //

 if (bcpres_array.size()!=gbox[0].size()*AMREX_SPACEDIM*2*nsolve_bicgstab)
  amrex::Error("bcpres_array size invalid");

 if (maskvals[level]->nGrow()!=1)
  amrex::Error("maskvals invalid ngrow");

#if (profile_solver==1)
 bprof.stop();
#endif

   // if openmp and no tiling, then tilegrid=validbox
   // and the grids are distributed amongst the threads.
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(inout); mfi.isValid(); ++mfi) {

#if (profile_solver==1)
  std::string subname="APPLYBC";
  std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
  popt_string_stream << cfd_project_option;
  std::string profname=subname+popt_string_stream.str();
  profname=profname+"_";
  std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
  lev_string_stream << level;
  profname=profname+lev_string_stream.str();

  BLProfiler bprof(profname);
#endif

  BL_ASSERT(gbox[level][mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid=mfi.tilebox(); 
  const Box& fabgrid=gbox[level][gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();        

  Vector<int> bcpres;
  bcpres.resize(2*AMREX_SPACEDIM*nsolve_bicgstab);
  int ibase=2*AMREX_SPACEDIM*gridno*nsolve_bicgstab;
  for (int i=0;i<2*AMREX_SPACEDIM*nsolve_bicgstab;i++)
   bcpres[i]=bcpres_array[i+ibase];
  FArrayBox& mfab=(*maskvals[level])[gridno];
  FArrayBox& bfab=pbdry[gridno];
  FORT_APPLYBC( 
   &nsolve_bicgstab,
   inout[gridno].dataPtr(),
   ARLIM(inout[gridno].loVect()),ARLIM(inout[gridno].hiVect()),
   bfab.dataPtr(),ARLIM(bfab.loVect()),ARLIM(bfab.hiVect()),
   mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
   bcpres.dataPtr(),
   tilelo,tilehi,
   fablo,fabhi,&bfact,&bfact_top);

#if (profile_solver==1)
  bprof.stop();
#endif
 }  // mfi
} // omp
 
 // no Barrier needed since FillBoundary called up above.  
}


void
ABecLaplacian::residual (MultiFab& residL,MultiFab& rhsL,
  MultiFab& solnL,int level,
  MultiFab& pbdry,Vector<int> bcpres_array) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::residual";
 std::stringstream popt_string_stream(std::stringstream::in |
   std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
   std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 int mg_coarsest_level=MG_numlevels_var-1;

 bool use_tiling=cfd_tiling;

 if (rhsL.boxArray()==laplacian_ones[level]->boxArray()) {
  // do nothing
 } else
  amrex::Error("rhsL.boxArray()!=laplacian_ones[level]->boxArray()");

 if (laplacian_ones[level]->nComp()==1) {
  // do nothing
 } else
  amrex::Error("laplacian_ones[level]->nComp()!=1");

#if (profile_solver==1)
 bprof.stop();
#endif

 if (residL.nGrow()!=0)
  amrex::Error("residL invalid ngrow");
  
 apply(residL,solnL,level,pbdry,bcpres_array);
 Fdiagsum(*MG_CG_diagsumL[level],level);
 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(solnL,use_tiling); mfi.isValid(); ++mfi) {

#if (profile_solver==1)
   std::string subname="RESIDL";
   std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
   popt_string_stream << cfd_project_option;
   std::string profname=subname+popt_string_stream.str();
   profname=profname+"_";
   std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
   lev_string_stream << level;
   profname=profname+lev_string_stream.str();

   BLProfiler bprof(profname);
#endif

   int nc = residL.nComp();
   if (nc!=nsolve_bicgstab)
    amrex::Error("nc invalid in residual");
   BL_ASSERT(gbox[level][mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid=mfi.tilebox();
   const Box& fabgrid=gbox[level][gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   FArrayBox& ones_FAB=(*laplacian_ones[level])[mfi];

     // in: LO_3D.F90
   FORT_RESIDL(
    &level,
    &mg_coarsest_level,
    &nsolve_bicgstab,
    ones_FAB.dataPtr(), 
    ARLIM(ones_FAB.loVect()),ARLIM(ones_FAB.hiVect()),
    residL[mfi].dataPtr(), 
    ARLIM(residL[mfi].loVect()),ARLIM(residL[mfi].hiVect()),
    rhsL[mfi].dataPtr(), 
    ARLIM(rhsL[mfi].loVect()),ARLIM(rhsL[mfi].hiVect()),
    residL[mfi].dataPtr(), 
    ARLIM(residL[mfi].loVect()),ARLIM(residL[mfi].hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,&bfact_top); 

     // mask off the residual if the off diagonal entries sum to 0.
   FArrayBox& diagfab=(*MG_CG_diagsumL[level])[mfi];
   FArrayBox& residfab=residL[mfi];
   residfab.mult(diagfab,tilegrid,0,0,nsolve_bicgstab); 

#if (profile_solver==1)
   bprof.stop();
#endif
  } // mfi
} // omp

}

void
ABecLaplacian::smooth(MultiFab& solnL,MultiFab& rhsL,
  int level,MultiFab& pbdry,Vector<int> bcpres_array,
  int smooth_type) {

    int nc = solnL.nComp();
    if (nc!=nsolve_bicgstab)
     amrex::Error("nc invalid in smooth");
    int ngrow_soln=solnL.nGrow();
    int ngrow_rhs=rhsL.nGrow();
    if (ngrow_soln<1)
     amrex::Error("ngrow_soln invalid");
    if (ngrow_rhs!=0)
     amrex::Error("ngrow_rhs invalid");

    applyBC(solnL,level,pbdry,bcpres_array);
    Fsmooth(solnL, rhsL, level, smooth_type);

}

// L2 norm
Real
ABecLaplacian::LPnorm(MultiFab &in, int level) const
{


#if (profile_solver==1)
 std::string subname="ABecLaplacian::norm";
 std::stringstream popt_string_stream(std::stringstream::in |
    std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
    std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 int nc = in.nComp();
 if (nc!=nsolve_bicgstab)
  amrex::Error("nc invalid in norm");

 Real mf_norm=0.0;
 for (int n=0;n<nc;n++) {
  Real test_norm=in.norm2(n);
  test_norm*=test_norm;
  mf_norm+=test_norm;
 }

#if (profile_solver==1)
 bprof.stop();
#endif

 return mf_norm;
}

// level==0 is the finest level
// level is the coarse level
// level-1 is the fine level
// avg=0  just sum
// avg=1  take avg
// avg=2  this is the ones_mf variable
void
ABecLaplacian::makeCoefficients (
	MultiFab& crs,
        const MultiFab& fine,
        int             level,
        int             avg)
{

  //
  // Determine index type of incoming MultiFab.
  //
 const IndexType iType(fine.boxArray()[0].ixType());

 const IndexType cType(D_DECL(IndexType::CELL,IndexType::CELL,IndexType::CELL));
 const IndexType xType(D_DECL(IndexType::NODE,IndexType::CELL,IndexType::CELL));
 const IndexType yType(D_DECL(IndexType::CELL,IndexType::NODE,IndexType::CELL));
 const IndexType zType(D_DECL(IndexType::CELL,IndexType::CELL,IndexType::NODE));

 int cdir;
 if (iType == cType) {
  cdir = -1;
 } else if (iType == xType) {
  cdir = 0;
 } else if (iType == yType) {
  cdir = 1;
 } else if ((iType == zType)&&(AMREX_SPACEDIM==3)) {
  cdir = 2;
 } else {
  amrex::Error("ABecLaplacian::makeCoeffients: Bad index type");
 }

#if (profile_solver==1)
 std::string subname="ABecLaplacian::makeCoefficients";
 std::stringstream popt_string_stream(std::stringstream::in |
     std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
     std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();
 profname=profname+"_";
 std::stringstream cdir_string_stream(std::stringstream::in |
     std::stringstream::out);
 cdir_string_stream << cdir;
 profname=profname+cdir_string_stream.str();

 BLProfiler bprof(profname);
#endif

 if (level>0) {
  // do nothing
 } else
  amrex::Error("level invalid");

 int flevel=level-1;
 int clevel=level;
 int bfact_coarse=bfact_array[clevel];
 int bfact_fine=bfact_array[flevel];
 int bfact_top=bfact_array[0];

 int nComp_expect = nsolve_bicgstab;
 if ((avg==0)||(avg==1)) {
  // do nothing
 } else if (avg==2) { // ones_mf variable
  nComp_expect=1;
 } else
  amrex::Error("avg invalid");

 if (nComp_expect==fine.nComp()) {
  // do nothing
 } else
  amrex::Error("nComp_expect!=fine.nComp()");

 BoxArray d(gbox[level]);
 if ((cdir>=0)&&(cdir<AMREX_SPACEDIM)) {
  d.surroundingNodes(cdir);
 } else if (cdir==-1) {
  // do nothing
 } else
  amrex::Error("cdir invalid");

   //
   // Only single-component solves supported (verified) by this class.
   //
 int nGrow=fine.nGrow();

 int ngrow_expect=0;

 if ((avg==0)||(avg==1)) {
  ngrow_expect=0;
 } else if (avg==2) { // ones_mf variable
  ngrow_expect=1;
  if (cdir==-1) {
   // do nothing
  } else
   amrex::Error("cdir invalid");
 } else
  amrex::Error("avg invalid");

 if (nGrow==ngrow_expect) {
  // do nothing
 } else
  amrex::Error("ngrow invalid in makecoeff");

 if (crs.boxArray()==d) {
  // do nothing
 } else
  amrex::Error("crs.boxArray() invalid");

 if (crs.nComp()==nComp_expect) {
  // do nothing
 } else
  amrex::Error("crs.nComp() invalid");

 if (crs.nGrow()==nGrow) {
  // do nothing
 } else
  amrex::Error("crs.nGrow() invalid");

 if ((avg==0)||(avg==1)) {
  crs.setVal(0.0,0,nComp_expect,nGrow); 
 } else if (avg==2) { // ones_mf variable
  crs.setVal(1.0,0,nComp_expect,nGrow); 
 } else
  amrex::Error("avg invalid");

 const BoxArray& grids = gbox[level]; // coarse grids

#if (profile_solver==1)
 bprof.stop();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(crs); mfi.isValid();++mfi) {

#if (profile_solver==1)
  std::string subname="AVERAGE_CC_or_EC";
  std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
  popt_string_stream << cfd_project_option;
  std::string profname=subname+popt_string_stream.str();
  profname=profname+"_";
  std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
  lev_string_stream << level;
  profname=profname+lev_string_stream.str();
  profname=profname+"_";
  std::stringstream cdir_string_stream(std::stringstream::in |
      std::stringstream::out);
  cdir_string_stream << cdir;
  profname=profname+cdir_string_stream.str();

  BLProfiler bprof(profname);
#endif

  const Box& fabgrid=grids[mfi.index()];
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  
  if (cdir==-1) {
   FORT_AVERAGECC(
     &nsolve_bicgstab,
     &nComp_expect,
     crs[mfi].dataPtr(), 
     ARLIM(crs[mfi].loVect()),
     ARLIM(crs[mfi].hiVect()),
     fine[mfi].dataPtr(),
     ARLIM(fine[mfi].loVect()),
     ARLIM(fine[mfi].hiVect()),
     fablo,fabhi,
     &avg,
     &nGrow,
     &bfact_coarse,&bfact_fine,&bfact_top);
  } else if ((cdir>=0)&&(cdir<AMREX_SPACEDIM)) {
   FORT_AVERAGEEC(
     &nComp_expect,
     crs[mfi].dataPtr(), 
     ARLIM(crs[mfi].loVect()),
     ARLIM(crs[mfi].hiVect()),
     fine[mfi].dataPtr(), 
     ARLIM(fine[mfi].loVect()),
     ARLIM(fine[mfi].hiVect()),
     fablo,fabhi,
     &cdir,&avg,
     &bfact_coarse,&bfact_fine,&bfact_top);
  } else
   amrex::Error("ABecLaplacian:: bad coefficient coarsening direction!");

#if (profile_solver==1)
  bprof.stop();
#endif
 }  // mfi
} //omp

 if ((avg==0)||(avg==1)) {
  // do nothing
 } else if (avg==2) { // ones_mf variable
  crs.FillBoundary(geomarray[level].periodicity());
 } else
  amrex::Error("avg invalid");

} // subroutine makeCoefficients


// level==0 is the finest level
void 
ABecLaplacian::buildMatrix() {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::buildMatrix";
 std::stringstream popt_string_stream(std::stringstream::in |
    std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();

 BLProfiler bprof(profname);
 bprof.stop();
#endif

 bool use_tiling=cfd_tiling;
 use_tiling=false;  // two different tile instances might mod the same data.

 int ncwork=AMREX_SPACEDIM*3+SCALAR_WORK_NCOMP;

 for (int level=0;level<MG_numlevels_var;level++) {
#if (profile_solver==1)
  bprof.start();
#endif

  if ((level>=1)&&(level<MG_numlevels_var)) {

   // remark: in the multigrid, average==1 too, so that one has:
   //  rhs[level-1] ~ volume_fine div u/dt
   //  rhs[level] ~ volume_coarse div u/dt  (1/2^d)
   //
   // divide by 4 in 2D and 8 in 3D
   // acoefs[level-1] ~ volume_fine
   // acoefs[level] ~ volume_coarse/ 2^d
   int avg=1;
   makeCoefficients(*acoefs[level],*acoefs[level-1],level,avg);
   makeCoefficients(*a_dual_coefs[level],*a_dual_coefs[level-1],level,avg);

   avg=2;
   makeCoefficients(*laplacian_ones[level],*laplacian_ones[level-1],level,avg);

   for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    // divide by 2 in 2D and 4 in 3D.
    // bcoefs[level-1] ~ area_fine/dx_fine  
    // bcoefs[level] ~ (area_coarse/dx_fine)/2^(d-1) =
    //                 4(area_coarse/dx_coarse)/2^d
    avg=1;
    makeCoefficients(*bcoefs[level][dir],*bcoefs[level-1][dir],level,avg);
    MultiFab& bmf=*bcoefs[level][dir];

    int ncomp=bmf.nComp();
    if (ncomp!=nsolve_bicgstab)
     amrex::Error("ncomp invalid");

    int ngrow=bmf.nGrow();
    if (ngrow!=0)
     amrex::Error("bcoefs should have ngrow=0");

     // after this step: bcoefs[level] ~ (area_coarse/dx_coarse)/2^d
    bmf.mult(0.25);
   } // dir=0..sdim-1

   offdiag_coeff[level] = 0.0;

   Real denom=0.0;
   if (AMREX_SPACEDIM==2) {
    denom=0.25;
   } else if (AMREX_SPACEDIM==3) {
    denom=0.125;
   } else
    amrex::Error("dimension bust");

   offdiag_coeff[level]=denom*offdiag_coeff[level-1];
   if (offdiag_coeff[level]>0.0) {
    // do nothing
   } else
    amrex::Error("offdiag_coeff[level] invalid");

  } else if (level==0) {
   // do nothing
  } else
   amrex::Error("level invalid");

  if (acoefs[level]->nGrow()!=nghostRHS)
   amrex::Error("acoefs[level]->nGrow() invalid");
  if (a_dual_coefs[level]->nGrow()!=nghostRHS)
   amrex::Error("a_dual_coefs[level]->nGrow() invalid");

  if (offdiag_coeff[level]>0.0) {
   // do nothing
  } else
   amrex::Error("offdiag_coeff[level] invalid");

  if (workcoefs[level]->nComp()==ncwork*nsolve_bicgstab) {
   // do nothing
  } else
   amrex::Error("workcoefs[level]->nComp() invalid");

  if (workcoefs[level]->nGrow()==nghostSOLN) {
   // do nothing
  } else 
   amrex::Error("workcoefs[level]->nGrow() invalid");

  workcoefs[level]->setVal(0.0,0,ncwork*nsolve_bicgstab,nghostSOLN);

  int bxleftcomp=0;
  int byleftcomp=bxleftcomp+1;
  int bzleftcomp=byleftcomp+AMREX_SPACEDIM-2;
  int bxrightcomp=bzleftcomp+1;
  int byrightcomp=bxrightcomp+1;
  int bzrightcomp=byrightcomp+AMREX_SPACEDIM-2;
  int icbxcomp=bzrightcomp+1;
  int icbycomp=icbxcomp+1;
  int icbzcomp=icbycomp+AMREX_SPACEDIM-2;

  int diag_non_singcomp=icbzcomp+1;
  int diag_dualcomp=diag_non_singcomp+1;
  int diag_comp=diag_dualcomp+1;
  int maskcomp=diag_comp+1;

  int icdiagcomp=maskcomp+1;
  int icdiagrbcomp=icdiagcomp+1;
  int axcomp=icdiagrbcomp+1;
  int solnsavecomp=axcomp+1;
  int rhssavecomp=solnsavecomp+1;
  int redsolncomp=rhssavecomp+1;
  int blacksolncomp=redsolncomp+1;

  if (workcoefs[level]->nComp()<=blacksolncomp*nsolve_bicgstab)
   amrex::Error("workcoefs[level] invalid");

  int bfact=bfact_array[level];
  int bfact_top=bfact_array[0];

#if (profile_solver==1)
  bprof.stop();
#endif

  for (int isweep=0;isweep<4;isweep++) {
   for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*workcoefs[level],use_tiling); mfi.isValid(); ++mfi) {

#if (profile_solver==1)
     std::string subname="BUILDMAT";
     std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
     popt_string_stream << cfd_project_option;
     std::string profname=subname+popt_string_stream.str();
     profname=profname+"_";
     std::stringstream lev_string_stream(std::stringstream::in |
       std::stringstream::out);
     lev_string_stream << level;
     profname=profname+lev_string_stream.str();
     profname=profname+"_";
     std::stringstream isweep_string_stream(std::stringstream::in |
      std::stringstream::out);
     isweep_string_stream << isweep;
     profname=profname+isweep_string_stream.str();
     profname=profname+"_";
     std::stringstream veldir_string_stream(std::stringstream::in |
      std::stringstream::out);
     veldir_string_stream << veldir;
     profname=profname+veldir_string_stream.str();

     BLProfiler bprof(profname);
#endif

     BL_ASSERT(gbox[level][mfi.index()] == mfi.validbox());
     const int gridno = mfi.index();
     const Box& tilegrid=mfi.tilebox();
     const Box& fabgrid=gbox[level][gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();

     int ofs=veldir*ncwork;

     FArrayBox& workFAB=(*workcoefs[level])[mfi];
     FArrayBox& onesFAB=(*laplacian_ones[level])[mfi];
     FArrayBox& adualFAB=(*a_dual_coefs[level])[mfi];
     FArrayBox& aFAB=(*acoefs[level])[mfi];
     FArrayBox& bxFAB=(*bcoefs[level][0])[mfi];
     FArrayBox& byFAB=(*bcoefs[level][1])[mfi];
     FArrayBox& bzFAB=(*bcoefs[level][AMREX_SPACEDIM-1])[mfi];

      // in: LO_3D.F90
     FORT_BUILDMAT(
      &level, // level==0 is finest
      &veldir,
      &nsolve_bicgstab,
      &isweep,
      &offdiag_coeff[level],
      &check_for_singular,
      &diag_regularization,
      onesFAB.dataPtr(),
      ARLIM(onesFAB.loVect()),ARLIM(onesFAB.hiVect()),
      adualFAB.dataPtr(veldir),
      ARLIM(adualFAB.loVect()),ARLIM(adualFAB.hiVect()),
      aFAB.dataPtr(veldir),
      ARLIM(aFAB.loVect()),ARLIM(aFAB.hiVect()),
      bxFAB.dataPtr(veldir),
      ARLIM(bxFAB.loVect()),ARLIM(bxFAB.hiVect()),
      byFAB.dataPtr(veldir),
      ARLIM(byFAB.loVect()),ARLIM(byFAB.hiVect()),
      bzFAB.dataPtr(veldir),
      ARLIM(bzFAB.loVect()),ARLIM(bzFAB.hiVect()),
      workFAB.dataPtr(diag_non_singcomp+ofs),
      ARLIM(workFAB.loVect()),
      ARLIM(workFAB.hiVect()),
      workFAB.dataPtr(diag_dualcomp+ofs),
      workFAB.dataPtr(diag_comp+ofs),
      workFAB.dataPtr(bxleftcomp+ofs),
      workFAB.dataPtr(bxrightcomp+ofs), 
      workFAB.dataPtr(byleftcomp+ofs),
      workFAB.dataPtr(byrightcomp+ofs), 
      workFAB.dataPtr(bzleftcomp+ofs),
      workFAB.dataPtr(bzrightcomp+ofs), 
      workFAB.dataPtr(icbxcomp+ofs), 
      workFAB.dataPtr(icbycomp+ofs), 
      workFAB.dataPtr(icbzcomp+ofs), 
      workFAB.dataPtr(icdiagcomp+ofs), 
      workFAB.dataPtr(icdiagrbcomp+ofs), 
      workFAB.dataPtr(maskcomp+ofs), 
      ARLIM(workFAB.loVect()),
      ARLIM(workFAB.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,&bfact,&bfact_top);

#if (profile_solver==1)
     bprof.stop();
#endif
    } // mfi
} // omp

   } // veldir=0...nsolve_bicgstab-1
  } // isweep=0..3

  MG_res[level]->setVal(0.0,0,nsolve_bicgstab,nghostRHS);
  MG_rhs[level]->setVal(0.0,0,nsolve_bicgstab,nghostRHS);
  MG_cor[level]->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);
  MG_pbdrycoarser[level]->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);
  if (level == 0) {
   MG_initialsolution->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);
  }

 } // level=0..MG_numlevels_var-1

} // end subroutine buildMatrix

ABecLaplacian::ABecLaplacian (
 const BoxArray& grids,
 const Geometry& geom,
 const DistributionMapping& dmap,
 int bfact,
 int cfd_level_in,
 int cfd_project_option_in,
 int nsolve_in,
 bool ns_tiling_in,
 int _use_mg_precond_at_top) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::ABecLaplacian";
 std::stringstream popt_string_stream(std::stringstream::in |
    std::stringstream::out);
 popt_string_stream << cfd_project_option_in;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
    std::stringstream::out);
 lev_string_stream << cfd_level_in;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 CG_use_mg_precond_at_top=_use_mg_precond_at_top;

 if (CG_use_mg_precond_at_top==0) {
  MG_numlevels_var=1;
  CG_numlevels_var=1;
 } else if (CG_use_mg_precond_at_top==1) {
  MG_numlevels_var = MG_numLevels(grids);
  if (MG_numlevels_var==1) {
   CG_numlevels_var=1;
  } else if (MG_numlevels_var>=2) {
   CG_numlevels_var=2;
  } else
   amrex::Error("MG_numlevels_var invalid");

  if (cfd_level_in==0) {
   // do nothing
  } else
   amrex::Error("cfd_level_in invalid");

 } else
  amrex::Error("CG_use_mg_precond_at_top invalid");

 laplacian_solvability=0;
 check_for_singular=0;
 diag_regularization=0.0;

 cfd_level=cfd_level_in;
 cfd_project_option=cfd_project_option_in;
 cfd_tiling=ns_tiling_in;

 CG_maxiter = CG_def_maxiter;
 CG_verbose = CG_def_verbose;

 CG_error_history.resize(CG_maxiter);
 for (int ehist=0;ehist<CG_error_history.size();ehist++) {
  for (int ih=0;ih<4;ih++)
   CG_error_history[ehist][ih]=0.0;
 }

 nsolve_bicgstab=nsolve_in; 

 gbox.resize(MG_numlevels_var);
 dmap_array.resize(MG_numlevels_var);
 geomarray.resize(MG_numlevels_var);
 bfact_array.resize(MG_numlevels_var);
 maskvals.resize(MG_numlevels_var);
 a_dual_coefs.resize(MG_numlevels_var,(MultiFab*)0);
 acoefs.resize(MG_numlevels_var,(MultiFab*)0);
 bcoefs.resize(MG_numlevels_var);
 for (int lev=0;lev<MG_numlevels_var;lev++) {
  for (int dir=0;dir<AMREX_SPACEDIM;dir++)
   bcoefs[lev][dir]=(MultiFab*)0;
 }
 offdiag_coeff.resize(MG_numlevels_var,0.0);
 workcoefs.resize(MG_numlevels_var,(MultiFab*)0);
 laplacian_ones.resize(MG_numlevels_var,(MultiFab*)0);
 MG_CG_diagsumL.resize(MG_numlevels_var,(MultiFab*)0);
 MG_CG_ones_mf_copy.resize(MG_numlevels_var,(MultiFab*)0);

 MG_res.resize(MG_numlevels_var, (MultiFab*)0);
 MG_rhs.resize(MG_numlevels_var, (MultiFab*)0);
 MG_cor.resize(MG_numlevels_var, (MultiFab*)0);
 MG_pbdrycoarser.resize(MG_numlevels_var, (MultiFab*)0);
 MG_CG_diagsumL.resize(MG_numlevels_var, (MultiFab*)0);
 MG_CG_ones_mf_copy.resize(MG_numlevels_var, (MultiFab*)0);

 GMRES_V_MF.resize(gmres_max_iter*nsolve_bicgstab);
 GMRES_Z_MF.resize(gmres_max_iter*nsolve_bicgstab);

 for (int coarsefine=0;coarsefine<CG_numlevels_var;coarsefine++) {
  for (int m=0;m<gmres_max_iter*nsolve_bicgstab;m++) {
   GMRES_V_MF[m][coarsefine]=(MultiFab*)0;
   GMRES_Z_MF[m][coarsefine]=(MultiFab*)0;
  }
  GMRES_W_MF[coarsefine]=(MultiFab*)0;
  CG_delta_sol[coarsefine]=(MultiFab*)0;
  CG_r[coarsefine]=(MultiFab*)0;
  CG_z[coarsefine]=(MultiFab*)0;
  CG_Av_search[coarsefine]=(MultiFab*)0;
  CG_p_search[coarsefine]=(MultiFab*)0;
  CG_v_search[coarsefine]=(MultiFab*)0;
  CG_rhs_resid_cor_form[coarsefine]=(MultiFab*)0;
  CG_pbdryhom[coarsefine]=(MultiFab*)0;
 } // coarsefine=0,...,CG_numlevels_var-1

 MG_initialsolution=(MultiFab*) 0; 
 
 for (int level=0;level<MG_numlevels_var;level++) {

  offdiag_coeff[level]=0.0;

  if (level==0) {

   gbox[level] = grids;

   dmap_array[level] = dmap;

   geomarray[level] = geom;
   bfact_array[level] = bfact;

   // mask=0 at fine/fine interfaces
   maskvals[level]=new MultiFab(grids,dmap_array[level],1,nghostSOLN,
      MFInfo().SetTag("maskvals level"),FArrayBoxFactory());
   maskvals[level]->setVal(0.0,0,1,nghostSOLN);
   maskvals[level]->setBndry(1.0);
   maskvals[level]->FillBoundary(geom.periodicity());

  } else if (level>0) {

   gbox[level] = gbox[level-1];
   gbox[level].coarsen(2);
   if (gbox[level].size()!=gbox[level-1].size())
    amrex::Error("gbox[level].size()!=gbox[level-1].size()");

   dmap_array[level] = dmap_array[level-1];

   geomarray[level].define(amrex::coarsen(geomarray[level-1].Domain(),2));
   bfact_array[level] = bfact_array[level-1];

   // mask=0 at fine/fine interfaces
   maskvals[level]=new MultiFab(gbox[level],dmap_array[level],1,nghostSOLN,
     MFInfo().SetTag("maskvals level"),FArrayBoxFactory());
   maskvals[level]->setVal(0.0,0,1,nghostSOLN);
   maskvals[level]->setBndry(1.0);
   maskvals[level]->FillBoundary(geomarray[level].periodicity());

  } else
   amrex::Error("level invalid");

  MG_CG_ones_mf_copy[level]=new MultiFab(gbox[level],dmap_array[level],
    1,nghostSOLN,
    MFInfo().SetTag("MG_CG_ones_mf_copy"),FArrayBoxFactory());
  MG_CG_ones_mf_copy[level]->setVal(1.0,0,1,nghostSOLN);

  MG_CG_diagsumL[level]=new MultiFab(gbox[level],dmap_array[level],
	nsolve_bicgstab,nghostRHS,
        MFInfo().SetTag("MG_CG_diagsumL"),FArrayBoxFactory());
  MG_CG_diagsumL[level]->setVal(0.0,0,nsolve_bicgstab,nghostRHS);

  a_dual_coefs[level]=new MultiFab(gbox[level],dmap_array[level],
	nsolve_bicgstab,nghostRHS,
        MFInfo().SetTag("a_dual_coefs"),FArrayBoxFactory());
  a_dual_coefs[level]->setVal(a_def,0,nsolve_bicgstab,nghostRHS);

  acoefs[level]=new MultiFab(gbox[level],dmap_array[level],
	nsolve_bicgstab,nghostRHS,
        MFInfo().SetTag("acoefs"),FArrayBoxFactory());
  acoefs[level]->setVal(a_def,0,nsolve_bicgstab,nghostRHS);

  int ncomp_work=(AMREX_SPACEDIM*3)+SCALAR_WORK_NCOMP;
  workcoefs[level]=new MultiFab(gbox[level],dmap_array[level],
    ncomp_work*nsolve_bicgstab,nghostSOLN,
    MFInfo().SetTag("workcoefs"),FArrayBoxFactory());
  workcoefs[level]->setVal(0.0,0,ncomp_work*nsolve_bicgstab,nghostSOLN);

  laplacian_ones[level]=new MultiFab(gbox[level],dmap_array[level],
	1,nghostSOLN,
        MFInfo().SetTag("laplacian_ones"),FArrayBoxFactory());
  laplacian_ones[level]->setVal(1.0,0,1,nghostSOLN);

  MG_res[level] = new MultiFab(gbox[level],dmap_array[level],
	 nsolve_bicgstab,nghostRHS,
	 MFInfo().SetTag("MG_res"),FArrayBoxFactory());
  MG_res[level]->setVal(0.0,0,nsolve_bicgstab,nghostRHS);

  MG_rhs[level] = new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostRHS,
    MFInfo().SetTag("MG_rhs"),FArrayBoxFactory());
  MG_rhs[level]->setVal(0.0,0,nsolve_bicgstab,nghostRHS);

  MG_cor[level] = new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostSOLN,
    MFInfo().SetTag("MG_cor"),FArrayBoxFactory());
  MG_cor[level]->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);

  MG_pbdrycoarser[level] = 
     new MultiFab(gbox[level],dmap_array[level],
	nsolve_bicgstab,nghostSOLN,
	MFInfo().SetTag("MG_pbdrycoarser"),FArrayBoxFactory());
  MG_pbdrycoarser[level]->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);

  // no ghost cells for edge or node coefficients
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
   BoxArray edge_boxes(gbox[level]);
   edge_boxes.surroundingNodes(dir);
   bcoefs[level][dir]=new MultiFab(edge_boxes,dmap_array[level],
     nsolve_bicgstab,0,
     MFInfo().SetTag("bcoefs"),FArrayBoxFactory());
   bcoefs[level][dir]->setVal(b_def,0,nsolve_bicgstab,0);
  }

 }  // level=0..MG_numlevels_var-1

 int finest_mg_level=0;

 MG_initialsolution = new MultiFab(gbox[finest_mg_level],
   dmap_array[finest_mg_level],
   nsolve_bicgstab,nghostSOLN,
   MFInfo().SetTag("MG_initialsolution"),FArrayBoxFactory());
 MG_initialsolution->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);

 MG_pbdryhom = new MultiFab(gbox[finest_mg_level],
   dmap_array[finest_mg_level],
   nsolve_bicgstab,nghostSOLN,
   MFInfo().SetTag("MG_pbdryhom"),FArrayBoxFactory());
 MG_pbdryhom->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);

 for (int coarsefine=0;coarsefine<CG_numlevels_var;coarsefine++) {

  int level=0;
  if (coarsefine==0) {
   level=0;
  } else if ((coarsefine==1)&&(coarsefine<CG_numlevels_var)) {
   level=MG_numlevels_var-1;
  } else
   amrex::Error("coarsefine invalid");

  for (int j=0;j<gmres_precond_iter_base_mg*nsolve_bicgstab;j++) {
   GMRES_V_MF[j][coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostRHS,
    MFInfo().SetTag("GMRES_V_MF"),FArrayBoxFactory()); 
   GMRES_Z_MF[j][coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostSOLN,
    MFInfo().SetTag("GMRES_Z_MF"),FArrayBoxFactory()); 

   GMRES_V_MF[j][coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostRHS);
   GMRES_Z_MF[j][coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);
  }
  GMRES_W_MF[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
   nsolve_bicgstab,nghostRHS,
   MFInfo().SetTag("GMRES_W_MF"),FArrayBoxFactory());

  CG_delta_sol[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostSOLN,
    MFInfo().SetTag("CG_delta_sol"),FArrayBoxFactory());
  CG_r[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostRHS,
    MFInfo().SetTag("CG_r"),FArrayBoxFactory());
  CG_z[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostSOLN,
    MFInfo().SetTag("CG_z"),FArrayBoxFactory());
  CG_Av_search[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostRHS,
    MFInfo().SetTag("CG_Av_search"),FArrayBoxFactory());
  CG_p_search[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostRHS,
    MFInfo().SetTag("CG_p_search"),FArrayBoxFactory());
  CG_v_search[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostRHS,
    MFInfo().SetTag("CG_v_search"),FArrayBoxFactory());
  CG_rhs_resid_cor_form[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostRHS,
    MFInfo().SetTag("CG_rhs_resid_cor_form"),FArrayBoxFactory());
  CG_pbdryhom[coarsefine]=new MultiFab(gbox[level],dmap_array[level],
    nsolve_bicgstab,nghostSOLN,
    MFInfo().SetTag("CG_pbdryhom"),FArrayBoxFactory());

 } // coarsefine=0..CG_numlevels_var-1


 ParmParse ppcg("cg");

 ppcg.query("maxiter", CG_def_maxiter);
 ppcg.query("v", CG_def_verbose);
 ppcg.query("verbose", CG_def_verbose);

 if (ParallelDescriptor::IOProcessor() && CG_def_verbose) {
  std::cout << "CGSolver settings...\n";
  std::cout << "   CG_def_maxiter   = " << CG_def_maxiter << '\n';
 }
    
 ParmParse ppmg("mg");
 ParmParse ppLp("Lp");

 ppmg.query("nu_0", MG_def_nu_0);
 ppmg.query("nu_f", MG_def_nu_f);
 ppmg.query("v", MG_def_verbose);
 ppmg.query("verbose", MG_def_verbose);
 ppmg.query("nu_b", MG_def_nu_b);

 if ((ParallelDescriptor::IOProcessor())&&(MG_def_verbose)) {
  std::cout << "MultiGrid settings...\n";
  std::cout << "   def_nu_0 =         " << MG_def_nu_0         << '\n';
  std::cout << "   def_nu_f =         " << MG_def_nu_f         << '\n';
  std::cout << "   def_nu_b =         " << MG_def_nu_b         << '\n';
 }

 MG_nu_0         = MG_def_nu_0;
 MG_nu_f         = MG_def_nu_f;
 MG_verbose      = MG_def_verbose;
 MG_nu_b         = MG_def_nu_b;

 if ((ParallelDescriptor::IOProcessor())&&(MG_verbose > 2)) {
  std::cout << "MultiGrid: " << MG_numlevels_var
    << " multigrid levels created for this solve" << '\n';
  std::cout << "Grids: " << '\n';
  BoxArray tmp = LPboxArray(0);
  for (int i = 0; i < MG_numlevels_var; ++i) {
   if (i > 0)
    tmp.coarsen(2);
   std::cout << " Level: " << i << '\n';
   for (int k = 0; k < tmp.size(); k++) {
    const Box& b = tmp[k];
    std::cout << "  [" << k << "]: " << b << "   ";
    for (int j = 0; j < AMREX_SPACEDIM; j++)
     std::cout << b.length(j) << ' ';
    std::cout << '\n';
   }
  }
 }


#if (profile_solver==1)
 bprof.stop();
#endif

}

ABecLaplacian::~ABecLaplacian ()
{
 for (int level=0;level<MG_numlevels_var;level++) {
  delete maskvals[level];
  maskvals[level]=(MultiFab*)0;
  delete MG_CG_ones_mf_copy[level];
  MG_CG_ones_mf_copy[level]=(MultiFab*)0;
  delete MG_CG_diagsumL[level];
  MG_CG_diagsumL[level]=(MultiFab*)0;
  delete a_dual_coefs[level];
  a_dual_coefs[level]=(MultiFab*)0;
  delete acoefs[level];
  acoefs[level]=(MultiFab*)0;
  delete workcoefs[level];
  workcoefs[level]=(MultiFab*)0;
  delete laplacian_ones[level];
  laplacian_ones[level]=(MultiFab*)0;
  delete MG_res[level];
  MG_res[level]=(MultiFab*)0;
  delete MG_rhs[level];
  MG_rhs[level]=(MultiFab*)0;
  delete MG_cor[level];
  MG_cor[level]=(MultiFab*)0;
  delete MG_pbdrycoarser[level];
  MG_pbdrycoarser[level]=(MultiFab*)0;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
   delete bcoefs[level][dir];
   bcoefs[level][dir]=(MultiFab*)0;
  }
 } // level=0..MG_numlevels-1

 delete MG_initialsolution;
 MG_initialsolution=(MultiFab*)0;

 delete MG_pbdryhom;
 MG_pbdryhom=(MultiFab*)0;

 for (int coarsefine=0;coarsefine<CG_numlevels_var;coarsefine++) {
  for (int j=0;j<gmres_precond_iter_base_mg*nsolve_bicgstab;j++) {
   delete GMRES_V_MF[j][coarsefine];
   GMRES_V_MF[j][coarsefine]=(MultiFab*)0;
   delete GMRES_Z_MF[j][coarsefine];
   GMRES_Z_MF[j][coarsefine]=(MultiFab*)0;
  }
  delete GMRES_W_MF[coarsefine];
  GMRES_W_MF[coarsefine]=(MultiFab*)0;

  delete CG_delta_sol[coarsefine];
  CG_delta_sol[coarsefine]=(MultiFab*)0;
  delete CG_r[coarsefine];
  CG_r[coarsefine]=(MultiFab*)0;
  delete CG_z[coarsefine];
  CG_z[coarsefine]=(MultiFab*)0;
  delete CG_Av_search[coarsefine];
  CG_Av_search[coarsefine]=(MultiFab*)0;
  delete CG_p_search[coarsefine];
  CG_p_search[coarsefine]=(MultiFab*)0;
  delete CG_v_search[coarsefine];
  CG_v_search[coarsefine]=(MultiFab*)0;
  delete CG_rhs_resid_cor_form[coarsefine];
  CG_rhs_resid_cor_form[coarsefine]=(MultiFab*)0;
  delete CG_pbdryhom[coarsefine];
  CG_pbdryhom[coarsefine]=(MultiFab*)0;
 } // coarsefine=0..CG_numlevels_var-1

}

void
ABecLaplacian::Fsmooth (MultiFab& solnL,
                        MultiFab& rhsL,
                        int level,
                        int smooth_type) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::Fsmooth";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 int mg_coarsest_level=MG_numlevels_var-1;

 bool use_tiling=cfd_tiling;

 int gsrb_timing=0;
  // double second() 
 double t1=0.0;
 double t2=0.0;

 const BoxArray& bxa = gbox[level];

 Real offdiag_coeff_level=offdiag_coeff[level];
 if (offdiag_coeff_level>0.0) {
  // do nothing
 } else
  amrex::Error("offdiag_coeff_level invalid");

#if (profile_solver==1)
 bprof.stop();
#endif

 const MultiFab & work=*workcoefs[level];

#if (profile_solver==1)
 bprof.start();
#endif

 int ncwork=AMREX_SPACEDIM*3+SCALAR_WORK_NCOMP;

 int nctest = work.nComp();
 if (nctest!=ncwork*nsolve_bicgstab)
  amrex::Error("ncwork invalid");

 int nc = solnL.nComp();
 if (nc!=nsolve_bicgstab)
  amrex::Error("nc bust");
 int ngrow=solnL.nGrow();
 if (ngrow<1)
  amrex::Error("ngrow solnL invalid");
 int ngrow_rhs=rhsL.nGrow();
 if (ngrow_rhs!=0)
  amrex::Error("ngrow rhsL invalid");

#if (profile_solver==1)
 bprof.stop();
#endif

 const MultiFab& ones_mf=*laplacian_ones[level];
 const BoxArray& temp_ba_test=ones_mf.boxArray();
 if ((temp_ba_test==bxa)&&(ones_mf.nComp()==1)) {
  // do nothing
 } else {
  amrex::Error("ones_mf invalid in fsmooth");
 }

 int bxleftcomp=0;
 int byleftcomp=bxleftcomp+1;
 int bzleftcomp=byleftcomp+AMREX_SPACEDIM-2;
 int bxrightcomp=bzleftcomp+1;
 int byrightcomp=bxrightcomp+1;
 int bzrightcomp=byrightcomp+AMREX_SPACEDIM-2;
 int icbxcomp=bzrightcomp+1;
 int icbycomp=icbxcomp+1;
 int icbzcomp=icbycomp+AMREX_SPACEDIM-2;

 int diag_non_singcomp=icbzcomp+1;
 int diag_dualcomp=diag_non_singcomp+1;
 int diag_comp=diag_dualcomp+1;
 int maskcomp=diag_comp+1;

 int icdiagcomp=maskcomp+1;
 int icdiagrbcomp=icdiagcomp+1;
 int axcomp=icdiagrbcomp+1;
 int solnsavecomp=axcomp+1;
 int rhssavecomp=solnsavecomp+1;
 int redsolncomp=rhssavecomp+1;
 int blacksolncomp=redsolncomp+1;

 int number_grids=gbox[level].size();

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

 int num_sweeps=0;
 if (smooth_type==0) // GSRB
  num_sweeps=6;
 else if (smooth_type==3) // Jacobi
  num_sweeps=4;
 else if (smooth_type==2) // ILU
  num_sweeps=7;
 else if (smooth_type==1) // ICRB
  num_sweeps=6;
 else
  amrex::Error("smooth_type invalid");


 for (int isweep=0;isweep<num_sweeps;isweep++) {

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(solnL,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(bxa[mfi.index()] == mfi.validbox());
   int gridno = mfi.index();
   if (gridno>=number_grids)
    amrex::Error("gridno invalid");
   const Box& tilegrid=mfi.tilebox();
   const Box& fabgrid=gbox[level][gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   if (gsrb_timing==1) 
    t1 = ParallelDescriptor::second();

   for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

#if (profile_solver==1)
    std::string subname="GSRB";
    std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
    popt_string_stream << cfd_project_option;
    std::string profname=subname+popt_string_stream.str();
    profname=profname+"_";
    std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
    lev_string_stream << level;
    profname=profname+lev_string_stream.str();
    profname=profname+"_";
    std::stringstream isweep_string_stream(std::stringstream::in |
      std::stringstream::out);
    isweep_string_stream << isweep;
    profname=profname+isweep_string_stream.str();
    profname=profname+"_";
    std::stringstream veldir_string_stream(std::stringstream::in |
      std::stringstream::out);
    veldir_string_stream << veldir;
    profname=profname+veldir_string_stream.str();

    BLProfiler bprof(profname);
#endif

    int ofs=veldir*ncwork;

    FORT_GSRB(
     &level,
     &mg_coarsest_level,
     &isweep,
     &num_sweeps,
     &offdiag_coeff_level,
     &check_for_singular,
     &diag_regularization,
     ones_mf[mfi].dataPtr(), 
     ARLIM(ones_mf[mfi].loVect()),ARLIM(ones_mf[mfi].hiVect()),
     solnL[mfi].dataPtr(veldir), 
     ARLIM(solnL[mfi].loVect()),ARLIM(solnL[mfi].hiVect()),
     rhsL[mfi].dataPtr(veldir), 
     ARLIM(rhsL[mfi].loVect()), ARLIM(rhsL[mfi].hiVect()),

     work[mfi].dataPtr(diag_non_singcomp+ofs),
     ARLIM(work[mfi].loVect()), ARLIM(work[mfi].hiVect()),
     work[mfi].dataPtr(diag_dualcomp+ofs),
     work[mfi].dataPtr(diag_comp+ofs),

     work[mfi].dataPtr(bxleftcomp+ofs),
     work[mfi].dataPtr(bxrightcomp+ofs), 
     work[mfi].dataPtr(byleftcomp+ofs),
     work[mfi].dataPtr(byrightcomp+ofs), 
     work[mfi].dataPtr(bzleftcomp+ofs),
     work[mfi].dataPtr(bzrightcomp+ofs), 
     work[mfi].dataPtr(icbxcomp+ofs), 
     work[mfi].dataPtr(icbycomp+ofs), 
     work[mfi].dataPtr(icbzcomp+ofs), 
     work[mfi].dataPtr(icdiagcomp+ofs), 
     work[mfi].dataPtr(icdiagrbcomp+ofs), 
     work[mfi].dataPtr(maskcomp+ofs), 
     work[mfi].dataPtr(axcomp+ofs), 
     work[mfi].dataPtr(solnsavecomp+ofs), 
     work[mfi].dataPtr(rhssavecomp+ofs), 
     work[mfi].dataPtr(redsolncomp+ofs), 
     work[mfi].dataPtr(blacksolncomp+ofs), 
     tilelo,tilehi,
     fablo,fabhi,&bfact,&bfact_top,
     &smooth_type);

#if (profile_solver==1)
    bprof.stop();
#endif
   } // veldir=0...nsolve_bicgstab-1

   if (gsrb_timing==1) {
    t2 = ParallelDescriptor::second();
    std::cout << "GSRB time, level= " << level << " smooth_type=" <<
     smooth_type << " gridno= " << gridno << " t2-t1=" << t2-t1 << '\n';
   }
  } // mfi
} // omp
 } // isweep

} // subroutine Fsmooth

// y=Ax
void
ABecLaplacian::Fapply (MultiFab& y,
                       MultiFab& x,
                       int level)
{

 int mg_coarsest_level=MG_numlevels_var-1;

 bool use_tiling=cfd_tiling;

 const BoxArray& bxa = gbox[level];

 Real offdiag_coeff_level=offdiag_coeff[level];
 if (offdiag_coeff_level>0.0) {
  // do nothing
 } else
  amrex::Error("offdiag_coeff_level invalid");

 const MultiFab & work=*workcoefs[level];
 int ncwork=AMREX_SPACEDIM*3+SCALAR_WORK_NCOMP;

 int nctest = work.nComp();
 if (nctest!=ncwork*nsolve_bicgstab)
  amrex::Error("ncwork invalid");

 int nc = y.nComp();
 if (nc!=nsolve_bicgstab)
  amrex::Error("nc bust");
 int ngrow_y=y.nGrow();
 if (ngrow_y!=0) {
  std::cout << "ngrow_y= " << ngrow_y << '\n';
  amrex::Error("ngrow y invalid");
 }
 int ngrow_x=x.nGrow();
 if (ngrow_x<1)
  amrex::Error("ngrow x invalid");

 const MultiFab& ones_mf=*laplacian_ones[level];
 const BoxArray& temp_ba_test=ones_mf.boxArray();
 if ((temp_ba_test==bxa)&&(ones_mf.nComp()==1)) {
  // do nothing
 } else {
  amrex::Error("ones_mf invalid in fapply");
 }

 int bxleftcomp=0;
 int byleftcomp=bxleftcomp+1;
 int bzleftcomp=byleftcomp+AMREX_SPACEDIM-2;
 int bxrightcomp=bzleftcomp+1;
 int byrightcomp=bxrightcomp+1;
 int bzrightcomp=byrightcomp+AMREX_SPACEDIM-2;
 int icbxcomp=bzrightcomp+1;
 int icbycomp=icbxcomp+1;
 int icbzcomp=icbycomp+AMREX_SPACEDIM-2;

 int diag_non_singcomp=icbzcomp+1;
 int diag_dualcomp=diag_non_singcomp+1;
 int diag_comp=diag_dualcomp+1;
 int maskcomp=diag_comp+1;

 int icdiagcomp=maskcomp+1;
 int icdiagrbcomp=icdiagcomp+1;
 int axcomp=icdiagrbcomp+1;
 int solnsavecomp=axcomp+1;
 int rhssavecomp=solnsavecomp+1;
 int redsolncomp=rhssavecomp+1;
 int blacksolncomp=redsolncomp+1;

 if (work.nComp()<=blacksolncomp*nsolve_bicgstab)
  amrex::Error("work.nComp() invalid");

 int number_grids=gbox[level].size();
 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(y,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(bxa[mfi.index()] == mfi.validbox());
  int gridno=mfi.index();
  if (gridno>=number_grids)
   amrex::Error("gridno invalid");
  const Box& tilegrid=mfi.tilebox();
  const Box& fabgrid=gbox[level][gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

#if (profile_solver==1)
   std::string subname="ADOTX";
   std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
   popt_string_stream << cfd_project_option;
   std::string profname=subname+popt_string_stream.str();
   profname=profname+"_";
   std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
   lev_string_stream << level;
   profname=profname+lev_string_stream.str();
   profname=profname+"_";
   std::stringstream veldir_string_stream(std::stringstream::in |
      std::stringstream::out);
   veldir_string_stream << veldir;
   profname=profname+veldir_string_stream.str();

   BLProfiler bprof(profname);
#endif

   int ofs=veldir*ncwork;

    // in: ABec_3D.F90
   FORT_ADOTX(
    &level,
    &mg_coarsest_level,
    &offdiag_coeff_level,
    &check_for_singular,
    &diag_regularization,
    ones_mf[mfi].dataPtr(), 
    ARLIM(ones_mf[mfi].loVect()),ARLIM(ones_mf[mfi].hiVect()),
    y[mfi].dataPtr(veldir),
    ARLIM(y[mfi].loVect()), ARLIM(y[mfi].hiVect()),
    x[mfi].dataPtr(veldir),
    ARLIM(x[mfi].loVect()), ARLIM(x[mfi].hiVect()),

    work[mfi].dataPtr(diag_non_singcomp+ofs),
    ARLIM(work[mfi].loVect()),ARLIM(work[mfi].hiVect()),
    work[mfi].dataPtr(diag_dualcomp+ofs),
    work[mfi].dataPtr(diag_comp+ofs),

    work[mfi].dataPtr(bxleftcomp+ofs),
    work[mfi].dataPtr(bxrightcomp+ofs),
    work[mfi].dataPtr(byleftcomp+ofs),
    work[mfi].dataPtr(byrightcomp+ofs),
    work[mfi].dataPtr(bzleftcomp+ofs),
    work[mfi].dataPtr(bzrightcomp+ofs),
    tilelo,tilehi,
    fablo,fabhi,&bfact,&bfact_top);

#if (profile_solver==1)
   bprof.stop();
#endif
  } // veldir=0..nsolve_bicgstab-1
 } // mfi
} // omp

} // end subroutine Fapply



void
ABecLaplacian::LP_update (MultiFab& sol,
   Real alpha,MultiFab& y,
   const MultiFab& p,int level) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::LP_update";
 std::stringstream popt_string_stream(std::stringstream::in |
   std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
   std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

    //
    // compute sol=y+alpha p  
    //
 bool use_tiling=cfd_tiling;

 if (level>=MG_numlevels_var)
  amrex::Error("level invalid in LP_update");

 const BoxArray& gboxlev = gbox[level];
 int ncomp = sol.nComp();
 if (ncomp==nsolve_bicgstab) {
  // do nothing
 } else
  amrex::Error("ncomp invalid");

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

#if (profile_solver==1)
 bprof.stop();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(sol,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(mfi.validbox() == gboxlev[mfi.index()]);
  const int gridno = mfi.index();
  const Box& tilegrid=mfi.tilebox();
  const Box& fabgrid=gbox[level][gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

#if (profile_solver==1)
   std::string subname="CGUPDATE";
   std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
   popt_string_stream << cfd_project_option;
   std::string profname=subname+popt_string_stream.str();
   profname=profname+"_";
   std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
   lev_string_stream << level;
   profname=profname+lev_string_stream.str();
   profname=profname+"_";
   std::stringstream veldir_string_stream(std::stringstream::in |
      std::stringstream::out);
   veldir_string_stream << veldir;
   profname=profname+veldir_string_stream.str();

   BLProfiler bprof(profname);
#endif

   FORT_CGUPDATE(
    sol[mfi].dataPtr(veldir),
    ARLIM(sol[mfi].loVect()), ARLIM(sol[mfi].hiVect()),
    &alpha,
    y[mfi].dataPtr(veldir),
    ARLIM(y[mfi].loVect()), ARLIM(y[mfi].hiVect()),
    p[mfi].dataPtr(veldir),
    ARLIM(p[mfi].loVect()), ARLIM(p[mfi].hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,&bfact_top);

#if (profile_solver==1)
   bprof.stop();
#endif
  } // veldir=0 ... nsolve_bicgstab-1
 }
} // omp

} // end subroutine LP_update


void ABecLaplacian::LP_dot(MultiFab& w,const MultiFab& p,
   int level,Real& result) {

#define profile_dot 0

#if (profile_dot==1)
  std::string subname="ABecLaplacian::LP_dot";
  std::stringstream popt_string_stream(std::stringstream::in |
   std::stringstream::out);
  popt_string_stream << cfd_project_option;
  std::string profname=subname+popt_string_stream.str();
  profname=profname+"_";
  std::stringstream lev_string_stream(std::stringstream::in |
   std::stringstream::out);
  lev_string_stream << level;
  profname=profname+lev_string_stream.str();

  BLProfiler bprof(profname);
#endif
 
 bool use_tiling=cfd_tiling;

 if (level>=MG_numlevels_var)
  amrex::Error("level invalid in LP_dot");

 if (level>=gbox.size()) {
  std::cout << "level= " << level << '\n';
  std::cout << "gboxsize= " << gbox.size() << '\n';
  std::cout << "num levels = " << MG_numlevels_var << '\n';
  amrex::Error("level exceeds gbox size");
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");

#if (profile_dot==1)
  bprof.stop();

  std::string subname2="ABecLaplacian::LP_dot_pw";
  std::string profname2=subname2+popt_string_stream.str();
  profname2=profname2+"_";
  profname2=profname2+lev_string_stream.str();

  BLProfiler bprof2(profname2);
#endif

 Vector<Real> pw;
 pw.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  pw[tid] = 0.0;
 }

#if (profile_dot==1)
  bprof2.stop();

  std::string subname3="ABecLaplacian::LP_dot_gboxlev";
  std::string profname3=subname3+popt_string_stream.str();
  profname3=profname3+"_";
  profname3=profname3+lev_string_stream.str();

  BLProfiler bprof3(profname3);
#endif 

 const BoxArray& gboxlev = gbox[level];
 int ncomp = p.nComp();
 if (ncomp!=nsolve_bicgstab)
  amrex::Error("ncomp invalid p");
 if (w.nComp()!=nsolve_bicgstab)
  amrex::Error("ncomp invalid w");

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

#if (profile_dot==1)
  bprof3.stop();

  std::string subname4="ABecLaplacian::LP_dot_MFIter";
  std::string profname4=subname4+popt_string_stream.str();
  profname4=profname4+"_";
  profname4=profname4+lev_string_stream.str();

  BLProfiler bprof4(profname4);
 
  bprof4.stop();
#endif

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
 for (MFIter mfi(w,use_tiling); mfi.isValid(); ++mfi) {

#if (profile_dot==1)
   std::string subname5="ABecLaplacian::LP_dot_MFIter_tilebox";
   std::string profname5=subname5+popt_string_stream.str();
   profname5=profname5+"_";
   profname5=profname5+lev_string_stream.str();

   BLProfiler bprof5(profname5);
#endif 

  Real tpw=0.0;
  BL_ASSERT(mfi.validbox() == gboxlev[mfi.index()]);
  const int gridno = mfi.index();
  const Box& tilegrid=mfi.tilebox();
  const Box& fabgrid=gbox[level][gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  int tid=0;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#endif
  if ((tid<0)||(tid>=thread_class::nthreads))
   amrex::Error("tid invalid");

#if (profile_dot==1)
   bprof5.stop();

   std::string subname6="CGXDOTY";
   std::string profname6=subname6+popt_string_stream.str();
   profname6=profname6+"_";
   profname6=profname6+lev_string_stream.str();

   BLProfiler bprof6(profname6);
#endif 

  FORT_CGXDOTY(
   &ncomp,
   &tpw, // init to 0.0d0 in CGXDOTY
   p[mfi].dataPtr(),ARLIM(p[mfi].loVect()), ARLIM(p[mfi].hiVect()),
   w[mfi].dataPtr(),ARLIM(w[mfi].loVect()), ARLIM(w[mfi].hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,&bfact_top);
  pw[tid] += tpw;

#if (profile_dot==1)
   bprof6.stop();
#endif 

 } // MFIter
} // omp

#if (profile_dot==1)
  bprof4.start();
#endif

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  pw[0]+=pw[tid];
 }

  // no Barrier needed since all processes must wait in order to receive the
  // reduced value.
 ParallelDescriptor::ReduceRealSum(pw[0]);

 result=pw[0];

#if (profile_dot==1)
  bprof4.stop();
#endif 

#undef profile_dot 

} // subroutine LP_dot

// onesCoefficients:
// =1 if diagonal <> 0
// =0 otherwise
void
ABecLaplacian::project_null_space(MultiFab& rhsL,int level) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::project_null_space";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);

 bprof.stop();
#endif

 if (laplacian_solvability==0) {
  if (check_for_singular==1) {
   if ((rhsL.nComp()==1)&&
       (laplacian_ones[level]->nComp()==1)) {
    MultiFab::Multiply(rhsL,*laplacian_ones[level],0,0,1,0);
   } else {
    std::cout << "laplacian_solvability= " << 
     laplacian_solvability << '\n';
    std::cout << "check_for_singular= " << 
     check_for_singular << '\n';
    amrex::Error("rhsL or ones_mf invalid nComp");
   }
  } else if (check_for_singular==0) {
   // do nothing
  } else
   amrex::Error("check_for_singular invalid");

 } else if (laplacian_solvability==1) {

  if (nsolve_bicgstab!=1)
   amrex::Error("nsolve_bicgstab invalid");

  if (check_for_singular==1) {
   if ((rhsL.nComp()==1)&&
       (laplacian_ones[level]->nComp()==1)) {
    MultiFab::Multiply(rhsL,*laplacian_ones[level],0,0,1,0);
   } else {
    std::cout << "laplacian_solvability= " << 
     laplacian_solvability << '\n';
    std::cout << "check_for_singular= " << 
     check_for_singular << '\n';
    amrex::Error("rhsL or ones_mf invalid nComp");
   }
   if ((laplacian_ones[level]->nComp()==1)&&
       (laplacian_ones[level]->nGrow()==1)) {
    MultiFab::Copy(*MG_CG_ones_mf_copy[level],
       *laplacian_ones[level],0,0,1,1);

    Real result,domainsum;
    LP_dot(rhsL,*laplacian_ones[level],level,result);
    LP_dot(*MG_CG_ones_mf_copy[level],
           *laplacian_ones[level],level,domainsum); 

    double total_cells=gbox[level].d_numPts();
    if (domainsum>total_cells) {
     std::cout << "level= " << level << '\n';
     std::cout << "result= " << result << '\n';
     std::cout << "cfd_level= " << cfd_level << '\n';
     std::cout << "cfd_project_option= " << cfd_project_option << '\n';
     std::cout << "domainsum= " << domainsum << '\n';
     std::cout << "total_cells= " << total_cells << '\n';
     amrex::Error("domainsum too big");
    }
    if (1==0) {
     std::cout << "cfd_level= " << cfd_level << '\n';
     std::cout << "cfd_project_option= " << cfd_project_option << '\n';
     std::cout << "level= " << level << '\n';
     std::cout << "result= " << result << '\n';
     std::cout << "domainsum= " << domainsum << '\n';
     std::cout << "total_cells= " << total_cells << '\n';
    }

    if (domainsum>=1.0) {
     Real coef=-result/domainsum;
      // rhsL=rhsL+coef * ones_mf
     LP_update(rhsL,coef,rhsL,*laplacian_ones[level],level); 
    } else if (domainsum==0.0) {
     // do nothing
    } else
     amrex::Error("domainsum invalid");

    MultiFab::Multiply(rhsL,*laplacian_ones[level],0,0,1,0);

    if (1==0) {
     std::cout << "check rhsL after projection \n";
     LP_dot(rhsL,*laplacian_ones[level],level,result);
     std::cout << "level= " << level << '\n';
     std::cout << "result= " << result << '\n';
    }

   } else
    amrex::Error("laplacian_ones[level]: ncomp or ngrow invalid");

  } else
   amrex::Error("check_for_singular invalid");

 } else
  amrex::Error("laplacian solvability incorrect");

} // subroutine project_null_space


// off diagonal sum flag 
void
ABecLaplacian::Fdiagsum(MultiFab&       y,
                       int             level) {

 bool use_tiling=cfd_tiling;

 const BoxArray& bxa = gbox[level];
FIX ME
 const MultiFab& a   = *acoefs[level];
 const MultiFab& bX  = *bcoefs[level][0];
 const MultiFab& bY  = *bcoefs[level][1];
 const MultiFab& bZ  = *bcoefs[level][AMREX_SPACEDIM-1];
 int nc = y.nComp();
 if (nc!=nsolve_bicgstab)
  amrex::Error("nc bust");

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(y,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(bxa[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid=mfi.tilebox();
  const Box& fabgrid=gbox[level][gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

#if (profile_solver==1)
   std::string subname="DIAGSUM";
   std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
   popt_string_stream << cfd_project_option;
   std::string profname=subname+popt_string_stream.str();
   profname=profname+"_";
   std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
   lev_string_stream << level;
   profname=profname+lev_string_stream.str();
   profname=profname+"_";
   std::stringstream veldir_string_stream(std::stringstream::in |
      std::stringstream::out);
   veldir_string_stream << veldir;
   profname=profname+veldir_string_stream.str();

   BLProfiler bprof(profname);
#endif

   FORT_DIAGSUM(
    y[mfi].dataPtr(veldir),
    ARLIM(y[mfi].loVect()), ARLIM(y[mfi].hiVect()),
    a[mfi].dataPtr(veldir), 
    ARLIM(a[mfi].loVect()), ARLIM(a[mfi].hiVect()),
    bX[mfi].dataPtr(veldir), 
    ARLIM(bX[mfi].loVect()), ARLIM(bX[mfi].hiVect()),
    bY[mfi].dataPtr(veldir), 
    ARLIM(bY[mfi].loVect()), ARLIM(bY[mfi].hiVect()),
    bZ[mfi].dataPtr(veldir), 
    ARLIM(bZ[mfi].loVect()), ARLIM(bZ[mfi].hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,&bfact_top);

#if (profile_solver==1)
   bprof.stop();
#endif
  } // veldir
 }
} // omp

} // end subroutine Fdiagsum



// z=K^{-1}r
void
ABecLaplacian::pcg_solve(
    MultiFab* z_in,
    MultiFab* r_in,
    Real eps_abs,Real bot_atol,
    MultiFab* pbdryhom_in,
    Vector<int> bcpres_array,
    int usecg_at_bottom,
    int smooth_type,int bottom_smooth_type,
    int presmooth,int postsmooth,
    int use_PCG,
    int level) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::pcg_solve";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);

 bprof.stop();
#endif

 project_null_space(*r_in,level);

 // z=K^{-1} r
 z_in->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);
 if (use_PCG==0) {
  MultiFab::Copy(*z_in,
		 *r_in,0,0,nsolve_bicgstab,nghostRHS);
 } else if (use_PCG==1) {
  if ((CG_use_mg_precond_at_top==1)&&
      (level==0)&&
      (MG_numlevels_var-1>0)) {

   if (CG_numlevels_var==2) {
    MG_solve(0,
      *z_in,
      *r_in,
      eps_abs,bot_atol,
      usecg_at_bottom,*pbdryhom_in,bcpres_array,
      smooth_type,bottom_smooth_type,
      presmooth,postsmooth);
   } else
    amrex::Error("CG_numlevels_var invalid");

  } else if ((CG_use_mg_precond_at_top==0)||
	     (level==MG_numlevels_var-1)) {
   for (int j=0;j<presmooth+postsmooth;j++) {
    smooth(*z_in,
	   *r_in,
	   level,
	   *pbdryhom_in,bcpres_array,smooth_type);
   }
  } else
   amrex::Error("use_mg_precond invalid");
 } else
  amrex::Error("use_PCG invalid");

 project_null_space(*z_in,level);

} // subroutine pcg_solve


// z=K^{-1}r
void
ABecLaplacian::pcg_GMRES_solve(
    int gmres_precond_iter,
    MultiFab* z_in,
    MultiFab* r_in,
    Real eps_abs,Real bot_atol,
    MultiFab* pbdryhom_in,
    Vector<int> bcpres_array,
    int usecg_at_bottom,
    int smooth_type,int bottom_smooth_type,
    int presmooth,int postsmooth,
    int use_PCG,
    int level) {

 if (nsolve_bicgstab>=1) {
  // do nothing
 } else
  amrex::Error("nsolve_bicgstab invalid");

 int coarsefine=0;
 if (level==0) {
  coarsefine=0;
 } else if ((level==MG_numlevels_var-1)&&(CG_numlevels_var==2)) {
  coarsefine=1;
 } else
  amrex::Error("level invalid");

 project_null_space((*r_in),level);

 // z=K^{-1} r
 z_in->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);

 if (gmres_precond_iter==0) {
   // Z=M^{-1}R
  pcg_solve(
   z_in,r_in,
   eps_abs,bot_atol,
   pbdryhom_in,
   bcpres_array,
   usecg_at_bottom,
   smooth_type,bottom_smooth_type,
   presmooth,postsmooth,
   use_PCG,
   level);

 } else if ((gmres_precond_iter>0)&&
            (gmres_precond_iter<=gmres_max_iter)) {

  int m=nsolve_bicgstab*gmres_precond_iter;

  Real* yy=new Real[m];

  Real** HH=new Real*[m+1];
  for (int i=0;i<m+1;i++) { 
   HH[i]=new Real[m];
   for (int j=0;j<m;j++) 
    HH[i][j]=0.0;
  }

  Real beta=0.0;
  LP_dot(*r_in,*r_in,level,beta);

  if (beta>0.0) {

   beta=sqrt(beta);

   for (int j=0;j<m;j++) {
    if (j>=gmres_precond_iter_base_mg*nsolve_bicgstab) {
     GMRES_V_MF[j][coarsefine]=new MultiFab(gbox[level],dmap_array[level],
       nsolve_bicgstab,nghostRHS,
       MFInfo().SetTag("GMRES_V_MF"),FArrayBoxFactory()); 
     GMRES_Z_MF[j][coarsefine]=new MultiFab(gbox[level],dmap_array[level],
       nsolve_bicgstab,nghostSOLN,
       MFInfo().SetTag("GMRES_Z_MF"),FArrayBoxFactory()); 
    }

    GMRES_V_MF[j][coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostRHS);
    GMRES_Z_MF[j][coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);
   } 

   Real aa=1.0/beta;
    // V0=V0+aa R
   CG_advance((*GMRES_V_MF[0][coarsefine]),aa,
              (*GMRES_V_MF[0][coarsefine]),(*r_in),level);

   int m_small=m;

   for (int j=0;j<m_small;j++) {
    GMRES_W_MF[coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostRHS);

     // Zj=M^{-1}Vj
    pcg_solve(
     GMRES_Z_MF[j][coarsefine],
     GMRES_V_MF[j][coarsefine],
     eps_abs,bot_atol,
     pbdryhom_in,
     bcpres_array,
     usecg_at_bottom,
     smooth_type,bottom_smooth_type,
     presmooth,postsmooth,
     use_PCG,
     level);

     // w=A Z
    apply(*GMRES_W_MF[coarsefine],
          *GMRES_Z_MF[j][coarsefine],
          level,*pbdryhom_in,bcpres_array);
    for (int i=0;i<=j;i++) {
     LP_dot(*GMRES_W_MF[coarsefine],
            *GMRES_V_MF[i][coarsefine],level,HH[i][j]);

     aa=-HH[i][j];
      // W=W+aa Vi
     CG_advance((*GMRES_W_MF[coarsefine]),aa,
                (*GMRES_W_MF[coarsefine]),
                (*GMRES_V_MF[i][coarsefine]),level);
    } // i=0..j

    LP_dot(*GMRES_W_MF[coarsefine],
           *GMRES_W_MF[coarsefine],level,HH[j+1][j]);

    if (HH[j+1][j]>0.0) {
     HH[j+1][j]=sqrt(HH[j+1][j]);
     if ((j>=0)&&(j<m-1)) {
      aa=1.0/HH[j+1][j];
       // V=V+aa W
      CG_advance((*GMRES_V_MF[j+1][coarsefine]),aa,
                 (*GMRES_V_MF[j+1][coarsefine]),
                 (*GMRES_W_MF[coarsefine]),level);
     } else if (j==m-1) {
      // do nothing
     } else
      amrex::Error("j invalid");
    } else if (HH[j+1][j]==0.0) {
     m_small=j;
    } else {
     std::cout << "HH[j+1][j]= " << HH[j+1][j] << '\n';
     amrex::Error("HH[j+1][j] invalid");
    }

   } // j=0..m-1
 
   int status=1;
   z_in->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);

   if (m_small==m) {
    int caller_id=3;
    GMRES_MIN_CPP(HH,beta,yy,m,caller_id,status);
   } else if ((m_small>=1)&&
              (m_small<m)) {
    Real** HH_small=new Real*[m_small+1];
    for (int i=0;i<m_small+1;i++) { 
     HH_small[i]=new Real[m_small];
     for (int j=0;j<m_small;j++) 
      HH_small[i][j]=HH[i][j];
    }
    int caller_id=4;
    GMRES_MIN_CPP(HH_small,beta,yy,m_small,caller_id,status);

    for (int i=0;i<m_small+1;i++) 
     delete [] HH_small[i];
    delete [] HH_small;
   } else if (m_small==0) {
    // Z=M^{-1}R
    pcg_solve(
     z_in,r_in,
     eps_abs,bot_atol,
     pbdryhom_in,
     bcpres_array,
     usecg_at_bottom,
     smooth_type,bottom_smooth_type,
     presmooth,postsmooth,
     use_PCG,
     level);
   } else {
    std::cout << "m_small= " << m_small << '\n';
    amrex::Error("CGSolver.cpp m_small invalid");
   }
 
   if (status==1) {
    for (int j=0;j<m_small;j++) {
     aa=yy[j];
      // Z=Z+aa Zj
     CG_advance((*z_in),aa,(*z_in),
                (*GMRES_Z_MF[j][coarsefine]),level);
    }
   } else
    amrex::Error("status invalid");

   for (int j=gmres_precond_iter_base_mg*nsolve_bicgstab;j<m;j++) {
    delete GMRES_V_MF[j][coarsefine];
    delete GMRES_Z_MF[j][coarsefine];
    GMRES_V_MF[j][coarsefine]=(MultiFab*)0;
    GMRES_Z_MF[j][coarsefine]=(MultiFab*)0;
   }

  } else if (beta==0.0) {

   // Z=M^{-1}R
   pcg_solve(
    z_in,r_in,
    eps_abs,bot_atol,
    pbdryhom_in,
    bcpres_array,
    usecg_at_bottom,
    smooth_type,bottom_smooth_type,
    presmooth,postsmooth,
    use_PCG,
    level);

  } else
   amrex::Error("beta invalid");

  for (int i=0;i<m+1;i++) 
   delete [] HH[i];
  delete [] HH;

  delete [] yy;
  
 } else {
  std::cout << "gmres_precond_iter= " << gmres_precond_iter << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "MG_numlevels_var= " << MG_numlevels_var << '\n';
  std::cout << "CG_numlevels_var= " << CG_numlevels_var << '\n';

  for (int ehist=0;ehist<CG_error_history.size();ehist++) {
   std::cout << "nit " << ehist << " CG_error_history[nit][0,1] " <<
    CG_error_history[ehist][2*coarsefine+0] << ' ' <<
    CG_error_history[ehist][2*coarsefine+1] << '\n';
  }

  amrex::Error("ABecLaplacian.cpp: gmres_precond_iter invalid");
 }

  // if v in nullspace(A) then we require that the solution satisfies:
  // z dot v = 0 which will desingularize the matrix system.
  // (z=z0 + c v,  z dot v=0 => c=-z0 dot v/(v dot v)
 project_null_space((*z_in),level);

} // subroutine pcg_GMRES_solve

void 
ABecLaplacian::CG_check_for_convergence(
  Real rnorm,Real rnorm_init,Real eps_abs,
  Real relative_error,int nit,int& error_close_to_zero,
  int level) {

 int critical_nit=10;
 Real critical_abs_tol=1.0e-14;
 Real critical_rel_tol=1.0e-14;
 if (critical_abs_tol>eps_abs)
  critical_abs_tol=eps_abs;
 if (critical_rel_tol>relative_error)
  critical_rel_tol=relative_error;

 int base_check=((rnorm<=eps_abs)||
                 (rnorm<=relative_error*rnorm_init));

 if (nit>critical_nit) {
  error_close_to_zero=base_check;
 } else if ((nit>=0)&&(nit<=critical_nit)) {
  error_close_to_zero=((rnorm<=critical_abs_tol)||
                       (rnorm<=critical_rel_tol*rnorm_init));

  if ((error_close_to_zero==0)&&
      (base_check==1))
   error_close_to_zero=2;

 } else
  amrex::Error("nit invalid");

 if (ParallelDescriptor::IOProcessor()) {
  if (CG_verbose>1) {
   std::cout << "in: CG_check_for_convergence nit= " << nit << 
	   " rnorm_init= " <<
	  rnorm_init << " rnorm= " << rnorm << '\n';
  }
 }

} // CG_check_for_convergence

void 
ABecLaplacian::CG_dump_params(Real rnorm,Real rnorm_init,
		Real eps_abs,Real relative_error,
                int is_bottom,Real bot_atol,
                int usecg_at_bottom,int smooth_type,
		int bottom_smooth_type,int presmooth,
		int postsmooth,MultiFab& mf1,
		MultiFab& mf2,int level) {

 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "level= " << level << '\n';
  std::cout << "rnorm_init,rnorm,eps_abs,relative_error " <<
   rnorm_init << ' ' << rnorm << ' ' <<  eps_abs << 
   ' ' << relative_error << '\n';
  std::cout << "is_bottom= " << is_bottom << '\n';
  std::cout << "bot_atol= " << bot_atol << '\n';
  std::cout << "usecg_at_bottom= " << usecg_at_bottom << '\n';
  std::cout << "smooth_type= " << smooth_type << '\n';
  std::cout << "bottom_smooth_type= " << bottom_smooth_type << '\n';
  std::cout << "presmooth= " << presmooth << '\n';
  std::cout << "postsmooth= " << postsmooth << '\n';
  std::cout << "cfd_level= " << cfd_level << '\n';
  std::cout << "cfd_project_option= " << cfd_project_option << '\n';
  std::cout << "laplacian_solvability= " << 
          laplacian_solvability << '\n';
  std::cout << "check_for_singular= " << 
          check_for_singular << '\n';
  std::cout << "nsolve_bicgstab= " << nsolve_bicgstab << '\n';
  std::cout << "gbox[0].size()= " << gbox[0].size() << '\n';
  std::cout << "numLevels()= " << MG_numlevels_var << '\n';
  std::cout << "mf1.boxArray()= " << mf1.boxArray() << '\n';
  std::cout << "LPnorm(mf1,lev)= " << LPnorm(mf1,level) << '\n';
  std::cout << "mf2.boxArray()= " << mf2.boxArray() << '\n';
  std::cout << "LPnorm(mf2,lev)= " << LPnorm(mf2,level) << '\n';
 }

} // end subroutine CG_dump_params

void
ABecLaplacian::CG_solve(
    int& cg_cycles_out,
    int nsverbose,int is_bottom,
    MultiFab& sol,
    MultiFab& rhs,
    Real eps_abs,Real bot_atol,
    MultiFab& pbdry,
    Vector<int> bcpres_array,
    int usecg_at_bottom,
    int& meets_tol,
    int smooth_type,int bottom_smooth_type,
    int presmooth,int postsmooth,
    Real& error_init,int level)
{

#if (profile_solver==1)
 std::string subname="ABecLaplacian::solve";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();
 BLProfiler bprof(profname);
#endif

 int coarsefine=0;
  // finest level of the hierarchy
 if (level==0) {
  coarsefine=0;
  // coarsest level of the hierarchy
 } else if ((level==MG_numlevels_var-1)&&(CG_numlevels_var==2)) {
  coarsefine=1;
 } else
  amrex::Error("level invalid CG_solve");

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
 BL_ASSERT(sol.boxArray() == LPboxArray(level));
 BL_ASSERT(rhs.boxArray() == LPboxArray(level));
 BL_ASSERT(pbdry.boxArray() == LPboxArray(level));

 Real relative_error=1.0e-12;

 int ncomp = sol.nComp();
 if (ncomp!=nsolve_bicgstab)
  amrex::Error("ncomp invalid");

 CG_pbdryhom[coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);

#if (profile_solver==1)
 bprof.stop();
#endif

 project_null_space(rhs,level);
 project_null_space(sol,level);

 CG_rhs_resid_cor_form[coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostRHS); 
 MultiFab::Copy(*CG_rhs_resid_cor_form[coarsefine],
   rhs,0,0,nsolve_bicgstab,nghostRHS);

 if ((CG_verbose>0)||(nsverbose>0)||(1==0))
  if (ParallelDescriptor::IOProcessor())
   std::cout << "CGSolver: is_bottom= " << is_bottom << '\n';

  // resid,rhs,soln
 residual((*CG_r[coarsefine]),(*CG_rhs_resid_cor_form[coarsefine]),
    sol,level,pbdry,bcpres_array);
 project_null_space((*CG_r[coarsefine]),level);

  // put solution and residual in residual correction form
 CG_delta_sol[coarsefine]->setVal(0.0,0,nsolve_bicgstab,1);
 MultiFab::Copy(*CG_rhs_resid_cor_form[coarsefine],
   *CG_r[coarsefine],0,0,nsolve_bicgstab,nghostRHS);

 Real rnorm = LPnorm(*CG_r[coarsefine], level);

#if (profile_solver==1)
 bprof.start();
#endif

 if (rnorm>=0.0) {
  rnorm=sqrt(rnorm);
 } else {
  amrex::Error("rnorm invalid");
 }
	 
 Real rnorm_init=rnorm;
 error_init=rnorm_init;

 if ((CG_verbose>0)||(nsverbose>0)) {
  if (ParallelDescriptor::IOProcessor()) {
   if (is_bottom==1)
    std::cout << "CGsolver(BOTTOM):Initial error(error_init)="<<rnorm << '\n';
   else
    std::cout << "CGsolver(NOBOT):Initial error(error_init)="<<rnorm << '\n';
  }
 }

 int nit=0;

 int use_PCG=1;

 if (use_PCG==0) {
  // do nothing
 } else if (CG_use_mg_precond_at_top==1) {
  // do nothing
 } else if (CG_use_mg_precond_at_top==0) {
  // do nothing
 } else {
  std::cout << "CG_use_mg_precond_at_top=" << 
	  CG_use_mg_precond_at_top << '\n';
  amrex::Error("CG_use_mg_precond_at_top invalid");
 }

 if (CG_z[coarsefine]->nComp()!=nsolve_bicgstab)
  amrex::Error("ncomp invalid");

 Real beta=0.0;
 Real rho=1.0;
 Real rho_old=1.0;
 Real omega=1.0;
 Real alpha=1.0;
 CG_p_search[coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostRHS); 
 CG_v_search[coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostRHS); 

 Real restart_tol=eps_abs*eps_abs*1.0e-4;
 int restart_flag=0;
 int prev_restart_flag=0;
 int gmres_precond_iter=gmres_precond_iter_base_mg;

 int error_close_to_zero=0;
 CG_check_for_convergence(rnorm,rnorm_init,eps_abs,relative_error,nit,
   error_close_to_zero,level);

 if (ParallelDescriptor::IOProcessor()) {
  if (CG_verbose>1) {
   if (is_bottom==1)
    std::cout << "CGSolver(BOT): rnorm_init,eps_abs,relative_error " <<
     rnorm_init << ' ' << eps_abs << ' ' << relative_error << '\n';
   else
    std::cout << "CGSolver(NOBOT): rnorm_init,eps_abs,relative_error " <<
     rnorm_init << ' ' << eps_abs << ' ' << relative_error << '\n';
  }
 }

#if (profile_solver==1)
 bprof.stop();
#endif

 for (int ehist=0;ehist<CG_error_history.size();ehist++) {
  for (int ih=0;ih<2;ih++)
   CG_error_history[ehist][2*coarsefine+ih]=0.0;
 }

 for(nit = 0;((nit < CG_maxiter)&&(error_close_to_zero!=1)); ++nit) {

  restart_flag=0;

  rho_old=rho;

  rnorm=LPnorm(*CG_r[coarsefine],level);
  if (rnorm>=0.0) {
   rnorm=sqrt(rnorm);
  } else {
   amrex::Error("rnorm invalid mglib");
  }
  if (nit==0)
   rnorm_init=rnorm;

  CG_error_history[nit][2*coarsefine]=rnorm;
  CG_error_history[nit][2*coarsefine+1]=eps_abs;

  CG_check_for_convergence(rnorm,rnorm_init,eps_abs,relative_error,nit,
		 error_close_to_zero,level);

  if (error_close_to_zero!=1) {

    // "meets_tol" informs the main solver in NavierStokes3.cpp 
    // whether it needs to continue.
   if ((nit==0)&&(rnorm>eps_abs*10.0)) {
    meets_tol=0;
   } else if ((nit==0)&&(rnorm<=eps_abs*10.0)) {
    // do nothing
   } else if (nit>0) {
    // do nothing
   } else
    amrex::Error("nit invalid");

    // rho=r0 dot r
   LP_dot(*CG_rhs_resid_cor_form[coarsefine],*CG_r[coarsefine],level,rho); 

   if (rho>=0.0) {
     if ((rho_old>restart_tol)&&(omega>restart_tol)) {
      beta=rho*alpha/(rho_old*omega);
       // p=p - omega v
      CG_advance( (*CG_p_search[coarsefine]),-omega,
        (*CG_p_search[coarsefine]),
        (*CG_v_search[coarsefine]),level );
       // p=r + beta p
      CG_advance( (*CG_p_search[coarsefine]),beta,(*CG_r[coarsefine]),
               (*CG_p_search[coarsefine]),level );
      project_null_space((*CG_p_search[coarsefine]),level);
       // z=K^{-1} p
      pcg_GMRES_solve(
       gmres_precond_iter,
       CG_z[coarsefine],
       CG_p_search[coarsefine],
       eps_abs,bot_atol,
       CG_pbdryhom[coarsefine],
       bcpres_array,
       usecg_at_bottom,
       smooth_type,bottom_smooth_type,
       presmooth,postsmooth,
       use_PCG,level);
       // v_search=A*z 
      apply(*CG_v_search[coarsefine],
            *CG_z[coarsefine],level,
            *CG_pbdryhom[coarsefine],bcpres_array);
     } else if ((rho_old<=restart_tol)||(omega<=restart_tol)) {
      restart_flag=1;
     } else
      amrex::Error("rho_old or omega invalid");
   } else if (rho<0.0) {
     restart_flag=1;
   } else
     amrex::Error("rho invalid mglib");

   if (restart_flag==0) {
     LP_dot(*CG_rhs_resid_cor_form[coarsefine],
            *CG_v_search[coarsefine],level,alpha);

     if (alpha>restart_tol) {
      alpha=rho/alpha;

       // x=x+alpha z
      LP_update( (*CG_delta_sol[coarsefine]), alpha, 
                 (*CG_delta_sol[coarsefine]),
		 (*CG_z[coarsefine]),level );
      project_null_space((*CG_delta_sol[coarsefine]),level);
      residual((*CG_r[coarsefine]),(*CG_rhs_resid_cor_form[coarsefine]),
        (*CG_delta_sol[coarsefine]),
	level,
        *CG_pbdryhom[coarsefine],
        bcpres_array); 
      project_null_space((*CG_r[coarsefine]),level);

      rnorm=LPnorm(*CG_r[coarsefine],level);
      if (rnorm>=0.0) {
       rnorm=sqrt(rnorm);
      } else {
       amrex::Error("rnorm invalid mglib");
      }

      CG_check_for_convergence(rnorm,rnorm_init,eps_abs,relative_error,nit,
	 error_close_to_zero,level);

      if (error_close_to_zero!=1) {
       // z=K^{-1} r
       pcg_GMRES_solve(
        gmres_precond_iter,
        CG_z[coarsefine],
	CG_r[coarsefine],
	eps_abs,bot_atol,
	CG_pbdryhom[coarsefine],
	bcpres_array,
        usecg_at_bottom,
	smooth_type,bottom_smooth_type,
        presmooth,postsmooth,
	use_PCG,
	level);
       // Av_search=A*z
       apply(*CG_Av_search[coarsefine],
	     *CG_z[coarsefine],level,
             *CG_pbdryhom[coarsefine],bcpres_array);
       Real rAz=0.0;
       Real zAAz=0.0;
       // rAz=(Az) dot r =z^T A^T r = z^T A^T K z >=0 if A and K SPD.
       LP_dot(*CG_Av_search[coarsefine],*CG_r[coarsefine],level,rAz);
       if (rAz>=0.0) {

        LP_dot(*CG_Av_search[coarsefine],
               *CG_Av_search[coarsefine],level,zAAz);
	//Az dot Az >=0 if A SPD
        if (zAAz>restart_tol) {

         omega=rAz/zAAz;
	 if (omega>=0.0) {
	  // do nothing
	 } else
	  amrex::Error("omega invalid mglib");

         // x=x+omega z
         LP_update( (*CG_delta_sol[coarsefine]), omega, 
                    (*CG_delta_sol[coarsefine]),
		    (*CG_z[coarsefine]),level );
         project_null_space((*CG_delta_sol[coarsefine]),level);
         residual(
	  (*CG_r[coarsefine]),
	  (*CG_rhs_resid_cor_form[coarsefine]),
          (*CG_delta_sol[coarsefine]),
    	  level,
          *CG_pbdryhom[coarsefine],
	  bcpres_array); 
         project_null_space((*CG_r[coarsefine]),level);
        } else if ((zAAz>=0.0)&&(zAAz<=restart_tol)) {
         restart_flag=1;
        } else
 	 amrex::Error("zAAz invalid");

       } else if (rAz<0.0) {
        restart_flag=1;
       } else
        amrex::Error("rAz invalid");

      } else if (error_close_to_zero==1) {
       // do nothing
      } else
       amrex::Error("error_close_to_zero invalid");
     } else if (alpha<=restart_tol) {
      restart_flag=1;
     } else
      amrex::Error("alpha invalid");
   } else if (restart_flag==1) {
     // do nothing
   } else 
     amrex::Error("restart_flag invalid");

   if (restart_flag==1) {

    if ((CG_verbose>0)||(nsverbose>0)) {
     if (ParallelDescriptor::IOProcessor()) {
      std::cout << "WARNING:RESTARTING: nit= " << nit << '\n';
      std::cout << "WARNING:RESTARTING: level= " << level << '\n';
      std::cout << "RESTARTING: gmres_precond_iter= " << 
       gmres_precond_iter << '\n';
      std::cout << "RESTARTING: rnorm= " << 
       rnorm << '\n';
      std::cout << "RESTARTING: nsolve_bicgstab= " << 
        nsolve_bicgstab << '\n';
     }
    } else if ((CG_verbose==0)&&(nsverbose==0)) {
     // do nothing
    } else
     amrex::Error("CG_verbose or nsverbose invalid");


    if (error_close_to_zero!=1) {
     beta=0.0;
     rho=1.0;
     rho_old=1.0;
     omega=1.0;
     alpha=1.0;
     CG_p_search[coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostRHS); 
     CG_v_search[coarsefine]->setVal(0.0,0,nsolve_bicgstab,nghostRHS); 
     sol.plus(*CG_delta_sol[coarsefine],0,nsolve_bicgstab,0);
     project_null_space(sol,level);
     CG_delta_sol[coarsefine]->setVal(0.0,0,nsolve_bicgstab,1);
     MultiFab::Copy(*CG_rhs_resid_cor_form[coarsefine],
       *CG_r[coarsefine],0,0,nsolve_bicgstab,0);
    } else if (error_close_to_zero==1) {
     amrex::Error("cannot have both restart_flag and error_close_to_zero");
    } else
     amrex::Error("error_close_to_zero invalid");

   } else if (restart_flag==0) {
    if (error_close_to_zero!=1) {
     // do nothing
    } else if (error_close_to_zero==1) {
     // do nothing
    } else
     amrex::Error("error_close_to_zero invalid");
   } else 
    amrex::Error("restart_flag invalid");

  } else if (error_close_to_zero==1) {
   // do nothing
  } else
   amrex::Error("error_close_to_zero invalid");

  if ((prev_restart_flag==1)&&(restart_flag==1)) {

   if ((error_close_to_zero==1)||
       (error_close_to_zero==2)) {
    gmres_precond_iter=gmres_precond_iter_base_mg;
   } else if (error_close_to_zero==0) {
    gmres_precond_iter=2*gmres_precond_iter;
   } else
    amrex::Error("error_close_to_zero invalid");

  } else if ((prev_restart_flag==0)&&(restart_flag==1)) {
   gmres_precond_iter=2*gmres_precond_iter_base_mg;
  } else if ((prev_restart_flag==1)&&(restart_flag==0)) {
   gmres_precond_iter=gmres_precond_iter_base_mg;
  } else if ((prev_restart_flag==0)&&(restart_flag==0)) {
   gmres_precond_iter=gmres_precond_iter_base_mg;
  } else
   amrex::Error("prev_restart_flag or restart_flag invalid");
        
  prev_restart_flag=restart_flag;

 }  // end of CGSolver loop

  // src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost
 sol.ParallelAdd(*CG_delta_sol[coarsefine],0,0,nsolve_bicgstab,0,0);
 project_null_space(sol,level);

 cg_cycles_out=nit;

 if ((CG_verbose>0)||(nsverbose>0)) {
  if (ParallelDescriptor::IOProcessor()) {
   if (CG_use_mg_precond_at_top==1) {
    if (is_bottom==0)
     std::cout << "mgpcg (mac) nit (NOBOT)" << nit << '\n';
    else
     std::cout << "mgpcg (mac) nit (BOT)" << nit << '\n';
   } else if (CG_use_mg_precond_at_top==0) {
    if (is_bottom==0) {
     std::cout << "pcg (mac) nit (NOBOT)" << nit << '\n';
    }
   } else
    amrex::Error("CG_use_mg_precond_at_top invalid");
  }
 }

 if (error_close_to_zero!=1) {

  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "Warning: ABecLaplacian:: failed to converge! \n";
   std::cout << "coarsefine= " << coarsefine << '\n';
   for (int ehist=0;ehist<CG_error_history.size();ehist++) {
    std::cout << "nit " << ehist << " CG_error_history[nit][0,1] " <<
     CG_error_history[ehist][2*coarsefine+0] << ' ' <<
     CG_error_history[ehist][2*coarsefine+1] << '\n';
   }
  }

  CG_dump_params(rnorm,rnorm_init,
    eps_abs,relative_error,
    is_bottom,bot_atol,
    usecg_at_bottom,smooth_type,
    bottom_smooth_type,presmooth, 
    postsmooth,
    sol,
    rhs,
    level);
 } else if (error_close_to_zero==1) {
  // do nothing
 } else {
  amrex::Error("error_close_to_zero invalid");
 }

 if (ncomp!=nsolve_bicgstab)
  amrex::Error("ncomp invalid");

 if ((CG_verbose>0)||(nsverbose>0)) {
  residual((*CG_r[coarsefine]),rhs,sol,level,pbdry,bcpres_array);

  Real testnorm=LPnorm(*CG_r[coarsefine],level);
  if (testnorm>=0.0) {
   testnorm=sqrt(testnorm);
  } else {
   amrex::Error("testnorm invalid mglib");
  }

  if (ParallelDescriptor::IOProcessor()) {
   if (is_bottom==1)
    std::cout << "residual non-homogeneous bc (BOT) " << testnorm << '\n';  
   else
    std::cout << "residual non-homogeneous bc (NOBOT)"<<testnorm << '\n';  
  }
 }

} // subroutine ABecLaplacian::CG_solve

// p=z+beta y
void ABecLaplacian::CG_advance (
       MultiFab& p,
       Real beta, 
       const MultiFab& z,
       MultiFab& y, 
       int level) {
    //
    // Compute p = z  +  beta y
    // only interior cells are updated.
    //
    const BoxArray& gbox_local = LPboxArray(level);
    int ncomp = p.nComp();
    if (ncomp!=nsolve_bicgstab)
     amrex::Error("p ncomp invalid");
    if (z.nComp()!=nsolve_bicgstab)
     amrex::Error("z ncomp invalid");
    if (y.nComp()!=nsolve_bicgstab)
     amrex::Error("y ncomp invalid");
    if ((p.nGrow()!=0)&&(p.nGrow()!=1))
     amrex::Error("p ngrow invalid");
    if ((z.nGrow()!=0)&&(z.nGrow()!=1))
     amrex::Error("z ngrow invalid");
    if ((y.nGrow()!=0)&&(y.nGrow()!=1))
     amrex::Error("y ngrow invalid");

    const BoxArray& zbox = z.boxArray();

    int bfact=get_bfact_array(level);
    int bfact_top=get_bfact_array(0);

    bool use_tiling=cfd_tiling;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(p,use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(zbox[mfi.index()] == gbox_local[mfi.index()]);
     int gridno=mfi.index();
     const Box& tilegrid=mfi.tilebox();
     const Box& fabgrid=gbox_local[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();

     for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

#if (profile_solver==1)
      std::string subname="CGADVCP";
      std::stringstream popt_string_stream(std::stringstream::in |
        std::stringstream::out);
      popt_string_stream << cfd_project_option;
      std::string profname=subname+popt_string_stream.str();
      profname=profname+"_";
      std::stringstream lev_string_stream(std::stringstream::in |
        std::stringstream::out);
      lev_string_stream << level;
      profname=profname+lev_string_stream.str();
      profname=profname+"_";
      std::stringstream veldir_string_stream(std::stringstream::in |
        std::stringstream::out);
      veldir_string_stream << veldir;
      profname=profname+veldir_string_stream.str();

      BLProfiler bprof(profname);
#endif

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

#if (profile_solver==1)
      bprof.stop();
#endif
     } // veldir
    } // mfi
} // omp

} // end subroutine advance

Real
ABecLaplacian::MG_errorEstimate(int level,
  MultiFab& pbdry,Vector<int> bcpres_array) {
  
 residual(*(MG_res[level]),*(MG_rhs[level]),*(MG_cor[level]), 
     level,pbdry,bcpres_array);
 MultiFab& resid = *(MG_res[level]);
 int ncomp=resid.nComp();
 if (ncomp!=nsolve_bicgstab)
  amrex::Error("ncomp invalid");

 Real local_error=sqrt(LPnorm(resid,level));
     
 return local_error;
}  // end subroutine MG_errorEstimate

void ABecLaplacian::MG_residualCorrectionForm (MultiFab& newrhs,
      MultiFab& oldrhs,MultiFab& solnL,
      MultiFab& inisol,MultiFab& pbdry,Vector<int> bcpres_array,
      int level) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::MG_residualCorrectionForm";
 std::stringstream popt_string_stream(std::stringstream::in |
   std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
   std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 if (solnL.nComp()!=nsolve_bicgstab)
  amrex::Error("ncomp invalid");
 if (solnL.nGrow()!=1)
  amrex::Error("solution should have ngrow=1");

 MG_initialsolution->ParallelCopy(inisol);
 solnL.ParallelCopy(inisol);

#if (profile_solver==1)
 bprof.stop();
#endif

 residual(newrhs, oldrhs, solnL, level, pbdry,bcpres_array);
 solnL.setVal(0.0,0,nsolve_bicgstab,1);
}

void
ABecLaplacian::MG_solve (int nsverbose,
  MultiFab& _sol, MultiFab& _rhs,
  Real _eps_abs,Real _atol_b,
  int usecg_at_bottom,MultiFab& pbdry,
  Vector<int> bcpres_array,
  int smooth_type,
  int bottom_smooth_type,int presmooth,
  int postsmooth) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::MG_solve0_";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();

 BLProfiler bprof(profname);
#endif

 int level = 0;

#if (profile_solver==1)
 bprof.stop();
#endif

 project_null_space(_rhs,level);
 MG_residualCorrectionForm(*MG_rhs[level],_rhs,*MG_cor[level],
     _sol,pbdry,bcpres_array,level);
 project_null_space(*MG_rhs[level],level);
 MG_pbdryhom->setVal(0.0,0,nsolve_bicgstab,nghostSOLN);

 MG_solve_(nsverbose,_sol, 
   _eps_abs, _atol_b, 
   *MG_pbdryhom,bcpres_array,
   usecg_at_bottom,
   smooth_type,
   bottom_smooth_type,presmooth,postsmooth);

} // subroutine MG_solve

// pbdry will always be identically zero since residual correction form.
void
ABecLaplacian::MG_solve_ (int nsverbose,MultiFab& _sol,
  Real eps_abs,Real atol_b,MultiFab& pbdry,Vector<int> bcpres_array,
  int usecg_at_bottom,
  int smooth_type,int bottom_smooth_type,int presmooth,int postsmooth) {

 int  level=0;

#if (profile_solver==1)
 std::string subname="ABecLaplacian::MG_solve_";
 std::stringstream popt_string_stream(std::stringstream::in |
   std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
   std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

  //
  // Relax system maxiter times, stop if 
  // absolute err <= _abs_eps
  //
 int ncomp=_sol.nComp();
 if (ncomp!=nsolve_bicgstab)
  amrex::Error("ncomp invalid");

 const Real error0 = MG_errorEstimate(level,pbdry,bcpres_array);
 Real error = error0;
 if ((ParallelDescriptor::IOProcessor()) && 
     ((MG_verbose)||(nsverbose>0))) {
   Spacer(std::cout, level);
   std::cout << "MultiGrid: Initial error (error0) = " << error0 << '\n';
 }

  //
  // Initialize correction to zero at this level (auto-filled at levels below)
  //
 (*MG_cor[level]).setVal(0.0);

#if (profile_solver==1)
 bprof.stop();
#endif

 
 MG_relax(*MG_cor[level],*MG_rhs[level],level,eps_abs,
   atol_b,usecg_at_bottom,pbdry,bcpres_array,
   smooth_type,bottom_smooth_type,
   presmooth,postsmooth);
 project_null_space(*MG_cor[level],level);

 error = MG_errorEstimate(level,pbdry,bcpres_array);

 if (ParallelDescriptor::IOProcessor()) {
  if (MG_verbose > 1 ) {
   Spacer(std::cout, level);
   std::cout << "MultiGrid: error/error0 "
             << error/error0 
             << " error " 
             << error << '\n';
  }
 }

 _sol.ParallelCopy(*MG_cor[level]);
   // src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost
 _sol.ParallelAdd(*MG_initialsolution,0,0,_sol.nComp(),0,0);
 project_null_space(_sol,level);

}  // subroutine MG_solve_

int
ABecLaplacian::MG_numLevels (const BoxArray& grids) const
{
 int MG_numLevelsMAX=1024;

 int ng = grids.size();
 int lv = MG_numLevelsMAX;
 //
 // The routine `falls through' since coarsening and refining
 // a unit box does not yield the initial box.
 //

 for (int i = 0; i < ng; ++i) {
  int llv = 0;
  Box tmp = grids[i];
  for (;;) {
      Box ctmp = tmp;   
      ctmp.coarsen(2);
      Box rctmp = ctmp; 
      rctmp.refine(2);
      if ((tmp != rctmp)||(ctmp.numPts() == 1))
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

} // end subroutine MG_numLevels

void
ABecLaplacian::MG_coarsestSmooth(MultiFab& solL,MultiFab& rhsL,
   int level,Real eps_abs,Real atol_b,int usecg_at_bottom,
   MultiFab& pbdry,Vector<int> bcpres_array,
   int smooth_type,int bottom_smooth_type,
   int presmooth,int postsmooth)
{


 int ncomp=solL.nComp();
 if (ncomp!=nsolve_bicgstab)
  amrex::Error("ncomp invalid");

 int is_bottom=1;

 if (usecg_at_bottom==0) {
  Real error0;
  if (MG_verbose) {
   error0 = MG_errorEstimate(level,pbdry,bcpres_array);
   if (ParallelDescriptor::IOProcessor())
    std::cout << "   Bottom Smoother: Initial error (error0) = " 
       << error0 << '\n';
  }

  project_null_space(rhsL,level);

  for (int i = MG_nu_f; i > 0; i--) {
   smooth(solL,rhsL,level,pbdry,bcpres_array,smooth_type); 

   if (MG_verbose > 1 || (i == 1 && MG_verbose)) {
    Real error = MG_errorEstimate(level,pbdry,bcpres_array);
    if (ParallelDescriptor::IOProcessor())
     std::cout << "   Bottom Smoother: Iteration " << i
       << " error/error0 " << error/error0 << " error " 
       << error << '\n';
   }
  }
 } else {
  int local_meets_tol=0;
  Real local_error0=0.0;
  int nsverbose=0;
  int cg_cycles_parm=0;

  CG_solve(
    cg_cycles_parm,
    nsverbose,is_bottom,
    solL,rhsL, atol_b, atol_b,
    pbdry,bcpres_array,usecg_at_bottom,
    local_meets_tol,
    bottom_smooth_type,bottom_smooth_type,
    presmooth,postsmooth,local_error0,level);
 }
}


void
ABecLaplacian::MG_relax (MultiFab& solL,MultiFab& rhsL,
   int level,Real eps_abs,
   Real atol_b,int usecg_at_bottom,
   MultiFab& pbdry,Vector<int> bcpres_array,
   int smooth_type,int bottom_smooth_type,int presmooth,
   int postsmooth) {

#if (profile_solver==1)
 std::string subname="ABecLaplacian::MG_relax";
 std::stringstream popt_string_stream(std::stringstream::in |
   std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
   std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 int ncomp=solL.nComp();
 if (ncomp!=nsolve_bicgstab)
  amrex::Error("ncomp invalid");

#if (profile_solver==1)
 bprof.stop();
#endif

 if (level < MG_numlevels_var - 1 ) {
 
  if (presmooth!=postsmooth)
   amrex::Error("presmooth must equal postsmooth for mgpcg");

  project_null_space(rhsL,level);

  for (int i = presmooth ; i > 0 ; i--) {
   smooth(solL,rhsL,level,pbdry,bcpres_array,smooth_type);
  }
  residual(*MG_res[level],rhsL,solL,level,pbdry,bcpres_array);
  project_null_space(*MG_res[level],level);
  MG_average(*MG_rhs[level+1], *MG_res[level],level+1,level);
  MG_cor[level+1]->setVal(0.0);

  if (!((usecg_at_bottom==0)||(usecg_at_bottom==1)))
   amrex::Error("usecg_at_bottom invalid");

  MG_pbdrycoarser[level+1]->setVal(0.0,0,nsolve_bicgstab,nghostSOLN); 
  for (int i = MG_def_nu_0; i > 0 ; i--) {
   MG_relax(*(MG_cor[level+1]),*(MG_rhs[level+1]),level+1,
    eps_abs,atol_b,usecg_at_bottom,
    *(MG_pbdrycoarser[level+1]),bcpres_array,
    smooth_type,bottom_smooth_type,presmooth,postsmooth);
  }

  MG_interpolate(solL, *(MG_cor[level+1]),level+1,level);
  for (int i = postsmooth; i > 0 ; i--) {
   smooth(solL, rhsL, level,pbdry,bcpres_array,smooth_type);
  }
 } else {
  MG_coarsestSmooth(solL,rhsL,level,eps_abs,atol_b,
   usecg_at_bottom,pbdry,bcpres_array,
   smooth_type,bottom_smooth_type,
   presmooth,postsmooth);
 }
} // subroutine MG_relax

void
ABecLaplacian::MG_average (MultiFab& c,MultiFab& f,
  int clevel,int flevel)
{

#if (profile_solver==1)
 std::string subname="ABecLaplacian::MG_average";
 std::stringstream popt_string_stream(std::stringstream::in |
   std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
   std::stringstream::out);
 lev_string_stream << clevel;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 if (clevel!=flevel+1)
  amrex::Error("clevel invalid");
 if (flevel<0)
  amrex::Error("flevel invalid"); 

 int bfact_coarse=get_bfact_array(clevel);
 int bfact_fine=get_bfact_array(flevel);
 int bfact_top=get_bfact_array(0);

#ifdef _OPENMP
#pragma omp parallel
#endif
 for (MFIter mfi(c); mfi.isValid(); ++mfi) {

  BL_ASSERT(c.boxArray().get(mfi.index()) == mfi.validbox());

  const Box& bx = mfi.validbox();

  int nc = c.nComp();
  if (nc!=nsolve_bicgstab) {
   std::cout << "nc,nsolve_bicgstab = " << nc << ' ' << 
    nsolve_bicgstab << '\n';
   amrex::Error("nc invalid in average");
  }
   // divide by 4 in 2D and 8 in 3D
  int iaverage=1;
  for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {
   FORT_AVERAGE(
    c[mfi].dataPtr(veldir),
    ARLIM(c[mfi].loVect()), ARLIM(c[mfi].hiVect()),
    f[mfi].dataPtr(veldir),
    ARLIM(f[mfi].loVect()), ARLIM(f[mfi].hiVect()),
    bx.loVect(), bx.hiVect(),&iaverage,
    &bfact_coarse,&bfact_fine,&bfact_top);
  }  // veldir
 } // mfi
#if (profile_solver==1)
 bprof.stop();
#endif
}  // end subroutine average

void
ABecLaplacian::MG_interpolate (MultiFab& f,MultiFab& c,
  int clevel,int flevel)
{

#if (profile_solver==1)
 std::string subname="ABecLaplacian::MG_interpolate";
 std::stringstream popt_string_stream(std::stringstream::in |
   std::stringstream::out);
 popt_string_stream << cfd_project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
   std::stringstream::out);
 lev_string_stream << clevel;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 if (clevel!=flevel+1)
  amrex::Error("clevel invalid");
 if (flevel<0)
  amrex::Error("flevel invalid"); 

 int bfact_coarse=get_bfact_array(clevel);
 int bfact_fine=get_bfact_array(flevel);
 int bfact_top=get_bfact_array(0);


 //
 // Use fortran function to interpolate up (prolong) c to f
 // Note: returns f=f+P(c) , i.e. ADDS interp'd c to f.
 //
#ifdef _OPENMP
#pragma omp parallel
#endif
 for (MFIter mfi(f); mfi.isValid(); ++mfi) {

  const Box& bx = c.boxArray()[mfi.index()];
  int nc = f.nComp();
  if (nc!=nsolve_bicgstab) {
   std::cout << "nc,nsolve_bicgstab = " << nc << ' ' << 
    nsolve_bicgstab << '\n';
   amrex::Error("nc invalid in interpolate");
  }

  for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {
   FORT_INTERP(
     &bfact_coarse,&bfact_fine,&bfact_top,
     f[mfi].dataPtr(veldir),
     ARLIM(f[mfi].loVect()), ARLIM(f[mfi].hiVect()),
     c[mfi].dataPtr(veldir),
     ARLIM(c[mfi].loVect()), ARLIM(c[mfi].hiVect()),
     bx.loVect(), bx.hiVect());
  } // veldir

 } // mfi

#if (profile_solver==1)
 bprof.stop();
#endif
}  // end subroutine interpolate


#undef profile_solver

}/* namespace amrex */

