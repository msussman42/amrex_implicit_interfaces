#include <winstd.H>
#if defined(BL_OLD_STL)
#include <stdlib.h>
#else
#include <cstdlib>
#endif

#include <algorithm>

#include <ParmParse.H>
#include <ABecLaplacian.H>
#include <LO_F.H>
#include <ABec_F.H>
#include <CG_F.H>
#include <ParallelDescriptor.H>

Real ABecLaplacian::a_def     = 0.0;
Real ABecLaplacian::b_def     = 1.0;

// level==0 is the finest level
ABecLaplacian::ABecLaplacian (const ABecLaplacian& _lp,int level) {

    nsolve_bicgstab=_lp.get_nsolve();

    cfd_level=_lp.cfd_level;
    cfd_project_option=_lp.cfd_project_option;
    cfd_tiling=_lp.cfd_tiling;

    gbox.resize(1);
    gbox[0] = _lp.boxArray(level);

    bfact_array.resize(1);
    bfact_array[0]=_lp.get_bfact_array(level);


    geomarray.resize(1);
    geomarray[0] = _lp.getGeom(level);
    laplacian_solvability=_lp.laplacian_solvability;
    check_for_singular=_lp.check_for_singular;
    diag_regularization=_lp.diag_regularization;

    if (_lp.numLevels()<=level)
     BoxLib::Error("numlevels too small");
}

// level==0 is the finest level
void
ABecLaplacian::apply (MultiFab& out,MultiFab& in,
  int level,MultiFab& pbdry,Array<int> bcpres_array) {

    applyBC(in,level,pbdry,bcpres_array);
    Fapply(out,in,level);
}


// level==0 is the finest level
void
ABecLaplacian::applyBC (MultiFab& inout,int level,
       MultiFab& pbdry,Array<int> bcpres_array) {

 prepareForLevel(level);

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

 if (inout.nGrow()!=1) {
  std::cout << "inout ngrow= " << inout.nGrow() << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "bfact_top= " << bfact_top << '\n';
  BoxLib::Error("inout ngrow<>1");
 }
 if (pbdry.nGrow()!=1)
  BoxLib::Error("pbdry ngrow<>1");
 
 if (inout.nComp()!=nsolve_bicgstab)
  BoxLib::Error("inout.nComp invalid");
 if (pbdry.nComp()!=nsolve_bicgstab)
  BoxLib::Error("pbdry.nComp invalid");

 inout.FillBoundary(geomarray[level].periodicity());

  //
  // Fill boundary cells.
  //

 if (bcpres_array.size()!=numGrids()*BL_SPACEDIM*2*nsolve_bicgstab)
  BoxLib::Error("bcpres_array size invalid");

 if (maskvals[level]->nGrow()!=1)
  BoxLib::Error("maskvals invalid ngrow");

   // if openmp and no tiling, then tilegrid=validbox
   // and the grids are distributed amongst the threads.
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(inout); mfi.isValid(); ++mfi) {
  BL_ASSERT(gbox[level][mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid=mfi.tilebox(); 
  const Box& fabgrid=gbox[level][gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();        

  Array<int> bcpres;
  bcpres.resize(2*BL_SPACEDIM*nsolve_bicgstab);
  int ibase=2*BL_SPACEDIM*gridno*nsolve_bicgstab;
  for (int i=0;i<2*BL_SPACEDIM*nsolve_bicgstab;i++)
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
 }  // mfi
} // omp
 ParallelDescriptor::Barrier();
}


void
ABecLaplacian::residual (MultiFab& residL,MultiFab& rhsL,
  MultiFab& solnL,int level,
  MultiFab& pbdry,Array<int> bcpres_array) {

 bool use_tiling=cfd_tiling;

 const BoxArray& temp_ba=rhsL.boxArray();
 MultiFab* diagsumL=new MultiFab(temp_ba,nsolve_bicgstab,0,Fab_allocate); 
  
 if (residL.nGrow()!=0)
  BoxLib::Error("residL invalid ngrow");
  
 apply(residL,solnL,level,pbdry,bcpres_array);
 Fdiagsum(*diagsumL,level);
 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(solnL,use_tiling); mfi.isValid(); ++mfi) {
   int nc = residL.nComp();
   if (nc!=nsolve_bicgstab)
    BoxLib::Error("nc invalid in residual");
   BL_ASSERT(gbox[level][mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid=mfi.tilebox();
   const Box& fabgrid=gbox[level][gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   FORT_RESIDL(
    &nsolve_bicgstab,
    residL[mfi].dataPtr(), 
    ARLIM(residL[mfi].loVect()),ARLIM(residL[mfi].hiVect()),
    rhsL[mfi].dataPtr(), 
    ARLIM(rhsL[mfi].loVect()),ARLIM(rhsL[mfi].hiVect()),
    residL[mfi].dataPtr(), 
    ARLIM(residL[mfi].loVect()),ARLIM(residL[mfi].hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,&bfact_top); 

     // mask off the residual if the off diagonal entries sum to 0.
   FArrayBox& diagfab=(*diagsumL)[mfi];
   FArrayBox& residfab=residL[mfi];
   residfab.mult(diagfab,tilegrid,0,0,nsolve_bicgstab); 
  } // mfi
} // omp
  ParallelDescriptor::Barrier();

  delete diagsumL;
}

void
ABecLaplacian::smooth(MultiFab& solnL,MultiFab& rhsL,
  int level,MultiFab& pbdry,Array<int> bcpres_array,
  int smooth_type) {

    int nc = solnL.nComp();
    if (nc!=nsolve_bicgstab)
     BoxLib::Error("nc invalid in smooth");
    int ngrow_soln=solnL.nGrow();
    int ngrow_rhs=rhsL.nGrow();
    if (ngrow_soln<1)
     BoxLib::Error("ngrow_soln invalid");
    if (ngrow_rhs!=0)
     BoxLib::Error("ngrow_rhs invalid");

    applyBC(solnL,level,pbdry,bcpres_array);
    Fsmooth(solnL, rhsL, level, smooth_type);

}

// L2 norm
Real
ABecLaplacian::norm(MultiFab &in, int level) const
{

 int nc = in.nComp();
 if (nc!=nsolve_bicgstab)
  BoxLib::Error("nc invalid in norm");

 Real mf_norm=0.0;
 for (int n=0;n<nc;n++) {
  Real test_norm=in.norm2(n);
  test_norm*=test_norm;
  mf_norm+=test_norm;
 }
 return mf_norm;
}

// level==0 is the finest level
// level is the coarse level
// level-1 is the fine level
// avg=0  just sum
// avg=1  take avg
// avg=2  this is the ones_mf variable
void
ABecLaplacian::makeCoefficients (MultiFab& cs,
                         const MultiFab& fn,
                         int             level,
                         int             avg)
{


 if (level<=0)
  BoxLib::Error("level invalid");

 int flevel=level-1;
 int clevel=level;
 int bfact_coarse=bfact_array[clevel];
 int bfact_fine=bfact_array[flevel];
 int bfact_top=bfact_array[0];

 int nComp_expect = nsolve_bicgstab;
 if ((avg==0)||(avg==1)) {
  // do nothing
 } else if (avg==2) {
  nComp_expect=1;
 } else
  BoxLib::Error("avg invalid");

 if (nComp_expect!=fn.nComp())
  BoxLib::Error("nComp_expect!=fn.nComp()");

  //
  // Determine index type of incoming MultiFab.
  //
 const IndexType iType(fn.boxArray()[0].ixType());

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
 } else if ((iType == zType)&&(BL_SPACEDIM==3)) {
  cdir = 2;
 } else {
  BoxLib::Error("ABecLaplacian::makeCoeffients: Bad index type");
 }

 BoxArray d(gbox[level]);
 if (cdir >= 0)
  d.surroundingNodes(cdir);
   //
   // Only single-component solves supported (verified) by this class.
   //
 int nGrow=fn.nGrow();

 int ngrow_expect=0;

 if ((avg==0)||(avg==1)) {
  ngrow_expect=0;
 } else if (avg==2) {
  ngrow_expect=1;
  if (cdir!=-1)
   BoxLib::Error("cdir invalid");
 } else
  BoxLib::Error("avg invalid");

 if (nGrow!=ngrow_expect)
  BoxLib::Error("ngrow invalid in makecoeff");

 cs.define(d, nComp_expect, nGrow, Fab_allocate);
 if ((avg==0)||(avg==1)) {
  // do nothing
 } else if (avg==2) {
  cs.setVal(1.0,0,nComp_expect,nGrow); 
 } else
  BoxLib::Error("avg invalid");

 const BoxArray& grids = gbox[level]; // coarse grids

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(cs); mfi.isValid();++mfi) {

  const Box& fabgrid=grids[mfi.index()];
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  
  if (cdir==-1) {
   FORT_AVERAGECC(
     &nsolve_bicgstab,
     &nComp_expect,
     cs[mfi].dataPtr(), 
     ARLIM(cs[mfi].loVect()),
     ARLIM(cs[mfi].hiVect()),
     fn[mfi].dataPtr(),
     ARLIM(fn[mfi].loVect()),
     ARLIM(fn[mfi].hiVect()),
     fablo,fabhi,
     &avg,
     &nGrow,
     &bfact_coarse,&bfact_fine,&bfact_top);
  } else if ((cdir>=0)&&(cdir<BL_SPACEDIM)) {
   FORT_AVERAGEEC(
     &nComp_expect,
     cs[mfi].dataPtr(),
     ARLIM(cs[mfi].loVect()),
     ARLIM(cs[mfi].hiVect()),
     fn[mfi].dataPtr(), 
     ARLIM(fn[mfi].loVect()),
     ARLIM(fn[mfi].hiVect()),
     fablo,fabhi,
     &cdir,&avg,
     &bfact_coarse,&bfact_fine,&bfact_top);
  } else
   BoxLib::Error("ABecLaplacian:: bad coefficient coarsening direction!");
 }  // mfi
} //omp
 ParallelDescriptor::Barrier();

 if ((avg==0)||(avg==1)) {
  // do nothing
 } else if (avg==2) {
  cs.FillBoundary(geomarray[level].periodicity());
 } else
  BoxLib::Error("avg invalid");

} // subroutine makeCoefficients


// level==0 is the finest level
void 
ABecLaplacian::buildMatrix(MultiFab& ones_mf,MultiFab& work,
     Real offdiag_coeff_level,
     MultiFab& a,MultiFab& bx,MultiFab& by,MultiFab& bz,
     int define_flag,int level) {

 bool use_tiling=cfd_tiling;
 use_tiling=false;  // two different tile instances might mod the same data.

 BoxArray d(gbox[level]);

 int nGrow_work=1;
 int nGrow_a=0;
 if (a.nGrow()!=nGrow_a)
  BoxLib::Error("a.ngrow invalid");

 if (offdiag_coeff_level>0.0) {
  // do nothing
 } else
  BoxLib::Error("offdiag_coeff_level invalid");

 int ncwork=BL_SPACEDIM*3+10;

 if (define_flag==1) {
    // bxs,nvar,ngrow,mem_mode
  work.define(d,ncwork*nsolve_bicgstab,nGrow_work,Fab_allocate);
 }
 work.setVal(0.0,0,ncwork*nsolve_bicgstab,nGrow_work);

 int bxleftcomp=0;
 int byleftcomp=bxleftcomp+1;
 int bzleftcomp=byleftcomp+BL_SPACEDIM-2;
 int bxrightcomp=bzleftcomp+1;
 int byrightcomp=bxrightcomp+1;
 int bzrightcomp=byrightcomp+BL_SPACEDIM-2;
 int icbxcomp=bzrightcomp+1;
 int icbycomp=icbxcomp+1;
 int icbzcomp=icbycomp+BL_SPACEDIM-2;

 int diag_non_singcomp=icbzcomp+1;
 int diag_singcomp=diag_non_singcomp+1;
 int maskcomp=diag_singcomp+1;

 int icdiagcomp=maskcomp+1;
 int icdiagrbcomp=icdiagcomp+1;
 int axcomp=icdiagrbcomp+1;
 int solnsavecomp=axcomp+1;
 int rhssavecomp=solnsavecomp+1;
 int redsolncomp=rhssavecomp+1;
 int blacksolncomp=redsolncomp+1;

 if (work.nComp()<=blacksolncomp*nsolve_bicgstab)
  BoxLib::Error("work.nComp() invalid");

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

 for (int isweep=0;isweep<4;isweep++) {
  for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(work,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(gbox[level][mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid=mfi.tilebox();
    const Box& fabgrid=gbox[level][gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();

    int ofs=veldir*ncwork;

    FORT_BUILDMAT(
     &level, // level==0 is finest
     &veldir,
     &nsolve_bicgstab,
     &isweep,
     &offdiag_coeff_level,
     &check_for_singular,
     &diag_regularization,
     ones_mf[mfi].dataPtr(),
     ARLIM(ones_mf[mfi].loVect()),ARLIM(ones_mf[mfi].hiVect()),
     a[mfi].dataPtr(veldir),
     ARLIM(a[mfi].loVect()),ARLIM(a[mfi].hiVect()),
     bx[mfi].dataPtr(veldir),
     ARLIM(bx[mfi].loVect()),ARLIM(bx[mfi].hiVect()),
     by[mfi].dataPtr(veldir),
     ARLIM(by[mfi].loVect()),ARLIM(by[mfi].hiVect()),
     bz[mfi].dataPtr(veldir),
     ARLIM(bz[mfi].loVect()),ARLIM(bz[mfi].hiVect()),

     work[mfi].dataPtr(diag_non_singcomp+ofs),
     ARLIM(work[mfi].loVect()),ARLIM(work[mfi].hiVect()),
     work[mfi].dataPtr(diag_singcomp+ofs),

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
     ARLIM(work[mfi].loVect()),ARLIM(work[mfi].hiVect()),
     tilelo,tilehi,
     fablo,fabhi,&bfact,&bfact_top);
   } // mfi
} // omp
   ParallelDescriptor::Barrier();

  } // veldir=0...nsolve_bicgstab-1
 } // isweep=0..3

} // end subroutine buildMatrix

ABecLaplacian::ABecLaplacian (const BoxArray& grids,const Geometry& geom,
 int bfact,
 int cfd_level_in,
 int cfd_project_option_in,
 int nsolve,
 bool ns_tiling_in) {

    laplacian_solvability=0;
    check_for_singular=0;
    diag_regularization=0.0;

    cfd_level=cfd_level_in;
    cfd_project_option=cfd_project_option_in;
    cfd_tiling=ns_tiling_in;

    nsolve_bicgstab=nsolve; 

    int level = 0;
    gbox.resize(1);
    gbox[level] = grids;
    geomarray.resize(1);
    geomarray[level] = geom;
    bfact_array.resize(1);
    bfact_array[level] = bfact;

    int ngrow_mask=1;

      // mask=0 at fine/fine interfaces
    maskvals.resize(1);
    maskvals[level]=new MultiFab(grids,1,ngrow_mask,Fab_allocate);
    maskvals[level]->setVal(0.0,0,1,ngrow_mask);
    maskvals[level]->setBndry(1.0);
    maskvals[level]->FillBoundary(geom.periodicity());
    initCoefficients(grids);
}

ABecLaplacian::~ABecLaplacian ()
{
    for (int i = 0; i < maskvals.size(); ++i) {
     delete maskvals[i];
    }
    clearToLevel(-1);
}


void
ABecLaplacian::clearToLevel (int level)
{
    BL_ASSERT(level >= -1);

    for (int i = level+1; i < numLevels(); ++i)
    {
        delete workcoefs[i]; // must follow cleangpu_data
        delete laplacian_ones[i];
        delete acoefs[i];

        a_valid[i] = false;

        for (int j = 0; j < BL_SPACEDIM; ++j) {
         delete bcoefs[i][j];
        }
        b_valid[i] = false;
        non_sing_valid[i] = false;
    }
}

// level==0 is the finest level

void
ABecLaplacian::prepareForLevel (int level)
{
 if (level<0)
  BoxLib::Error("level invalid");

 if (level == 0 )
  return;

 prepareForLevel(level-1);

 if (gbox.size() > level) return;

 if (gbox.size()!=level)
  BoxLib::Error("gbox size problem");

 geomarray.resize(level+1);
 geomarray[level].define(BoxLib::coarsen(geomarray[level-1].Domain(),2));

 gbox.resize(level+1);
 gbox[level] = gbox[level-1];
 gbox[level].coarsen(2);

 bfact_array.resize(level+1);
 bfact_array[level] = bfact_array[level-1];
 
 if (maskvals.size()!=level)
  BoxLib::Error("maskvals size problem");
 maskvals.resize(level+1);
 maskvals[level]=new MultiFab(gbox[level],1,1,Fab_allocate);
 maskvals[level]->setVal(0.0,0,1,1);
 maskvals[level]->setBndry(1.0);
 maskvals[level]->FillBoundary(geomarray[level].periodicity());

 int need_to_build=0;

 if (level >= a_valid.size() || a_valid[level] == false) {
  need_to_build=1;

  if (acoefs.size() < level+1) {
   acoefs.resize(level+1);

   workcoefs.resize(level+1);
   laplacian_ones.resize(level+1);
  } else {
   delete acoefs[level];

   delete workcoefs[level]; // must follow cleangpu_data
   delete laplacian_ones[level];
  }
  acoefs[level] = new MultiFab;
  workcoefs[level] = new MultiFab;
  laplacian_ones[level] = new MultiFab;

   // remark: in the multigrid, average==1 too, so that one has:
   //  rhs[level-1] ~ volume_fine div u/dt
   //  rhs[level] ~ volume_coarse div u/dt  (1/2^d)
   //
   // divide by 4 in 2D and 8 in 3D
   // acoefs[level-1] ~ volume_fine
   // acoefs[level] ~ volume_coarse/ 2^d
  int avg=1;
  makeCoefficients(*acoefs[level],*acoefs[level-1],level,avg);
  a_valid.resize(level+1);
  a_valid[level] = true;

  avg=2;
  makeCoefficients(*laplacian_ones[level],*laplacian_ones[level-1],level,avg);
 }

 if (level >= b_valid.size() || b_valid[level] == false) {

  if (bcoefs.size() < level+1) {
   bcoefs.resize(level+1);
   for(int i = 0; i < BL_SPACEDIM; ++i)
    bcoefs[level][i] = new MultiFab;
  } else {
   for(int i = 0; i < BL_SPACEDIM; ++i) {
    delete bcoefs[level][i];
    bcoefs[level][i] = new MultiFab;
   }
  }
  for (int i = 0; i < BL_SPACEDIM; ++i) {
    // divide by 2 in 2D and 4 in 3D.
    // bcoefs[level-1] ~ area_fine/dx_fine  
    // bcoefs[level] ~ (area_coarse/dx_fine)/2^(d-1) =
    //                 4(area_coarse/dx_coarse)/2^d
   int avg=1;
   makeCoefficients(*bcoefs[level][i],*bcoefs[level-1][i],level,avg);
   MultiFab& bmf=*bcoefs[level][i];

   int ncomp=bmf.nComp();
   if (ncomp!=nsolve_bicgstab)
    BoxLib::Error("ncomp invalid");

   int ngrow=bmf.nGrow();
   if (ngrow!=0)
    BoxLib::Error("bcoefs should have ngrow=0");

    // after this step: bcoefs[level] ~ (area_coarse/dx_coarse)/2^d
   bmf.mult(0.25);
  }
  b_valid.resize(level+1);
  b_valid[level] = true;
 }

 if ((level >= non_sing_valid.size())|| 
     (non_sing_valid[level] == false)) {

  if (offdiag_coeff.size() < level+1) {
   offdiag_coeff.resize(level+1);
   offdiag_coeff[level] = 0.0;
  } else {
   offdiag_coeff[level] = 0.0;
  }
  Real denom=0.0;
  if (BL_SPACEDIM==2) {
   denom=0.25;
  } else if (BL_SPACEDIM==3) {
   denom=0.125;
  } else
   BoxLib::Error("dimension bust");

  offdiag_coeff[level]=denom*offdiag_coeff[level-1];
  if (offdiag_coeff[level]>0.0) {
   // do nothing
  } else
   BoxLib::Error("offdiag_coeff[level] invalid");

  non_sing_valid.resize(level+1);
  non_sing_valid[level] = true;
 }

  // acoefs[level] and bcoefs[level] might be modified depending on
  // the contents of laplacian_ones[level]
 if (need_to_build==1) {
  int define_flag=1;
  buildMatrix(*laplacian_ones[level],*workcoefs[level],
   offdiag_coeff[level],
   *acoefs[level],*bcoefs[level][0],*bcoefs[level][1],
   *bcoefs[level][BL_SPACEDIM-1],define_flag,level);
 }
} // subroutine prepareForLevel

void
ABecLaplacian::initCoefficients (const BoxArray& _ba)
{
 const int nComp=nsolve_bicgstab;
 const int nGrow=1;
 const int nGrow_acoef=0;

 acoefs.resize(1);
 bcoefs.resize(1);
 offdiag_coeff.resize(1);

 offdiag_coeff[0]=0.0;

 workcoefs.resize(1);
 laplacian_ones.resize(1);

 acoefs[0] = new MultiFab(_ba, nComp, nGrow_acoef,Fab_allocate);
 acoefs[0]->setVal(a_def,0,nComp,nGrow_acoef);

 int ncomp_work=(BL_SPACEDIM*3)+10;
 workcoefs[0] = new MultiFab(_ba, ncomp_work*nComp, nGrow,Fab_allocate);
 workcoefs[0]->setVal(0.0,0,ncomp_work*nComp,nGrow);
 laplacian_ones[0] = new MultiFab(_ba, 1, nGrow,Fab_allocate);
 laplacian_ones[0]->setVal(1.0,0,1,nGrow);

 a_valid.resize(1);
 a_valid[0] = true;

  // no ghost cells for edge or node coefficients
 for (int i = 0; i < BL_SPACEDIM; ++i) {
  BoxArray edge_boxes(_ba);
  edge_boxes.surroundingNodes(i);
  bcoefs[0][i] = new MultiFab(edge_boxes,nComp,0,Fab_allocate);
  bcoefs[0][i]->setVal(b_def,0,nComp,0);
 }

 b_valid.resize(1);
 b_valid[0] = true;

 non_sing_valid.resize(1);
 non_sing_valid[0]=true;
}

void
ABecLaplacian::invalidate_a_to_level (int lev)
{
 lev = (lev >= 0 ? lev : 0);
 for (int i = lev; i < numLevels(); i++)
  a_valid[i] = false;
}

void
ABecLaplacian::invalidate_b_to_level (int lev)
{
 lev = (lev >= 0 ? lev : 0);
 for (int i = lev; i < numLevels(); i++)
  b_valid[i] = false;
}


void
ABecLaplacian::invalidate_non_sing_to_level (int lev)
{
 lev = (lev >= 0 ? lev : 0);
 for (int i = lev; i < numLevels(); i++)
  non_sing_valid[i] = false;
}

//
// Must be defined for MultiGrid/CGSolver to work.
//

void
ABecLaplacian::Fsmooth (MultiFab& solnL,
                        MultiFab& rhsL,
                        int level,
                        int smooth_type) {

 bool use_tiling=cfd_tiling;

 int gsrb_timing=0;
 Real t1=0.0;
 Real t2=0.0;

 const BoxArray& bxa = gbox[level];

 Real offdiag_coeff_level=offdiag_coeff[level];
 if (offdiag_coeff_level>0.0) {
  // do nothing
 } else
  BoxLib::Error("offdiag_coeff_level invalid");

 const MultiFab & work=matCoefficients(level);
 int ncwork=BL_SPACEDIM*3+10;

 int nctest = work.nComp();
 if (nctest!=ncwork*nsolve_bicgstab)
  BoxLib::Error("ncwork invalid");

 int nc = solnL.nComp();
 if (nc!=nsolve_bicgstab)
  BoxLib::Error("nc bust");
 int ngrow=solnL.nGrow();
 if (ngrow<1)
  BoxLib::Error("ngrow solnL invalid");
 int ngrow_rhs=rhsL.nGrow();
 if (ngrow_rhs!=0)
  BoxLib::Error("ngrow rhsL invalid");

 int bxleftcomp=0;
 int byleftcomp=bxleftcomp+1;
 int bzleftcomp=byleftcomp+BL_SPACEDIM-2;
 int bxrightcomp=bzleftcomp+1;
 int byrightcomp=bxrightcomp+1;
 int bzrightcomp=byrightcomp+BL_SPACEDIM-2;
 int icbxcomp=bzrightcomp+1;
 int icbycomp=icbxcomp+1;
 int icbzcomp=icbycomp+BL_SPACEDIM-2;

 int diag_non_singcomp=icbzcomp+1;
 int diag_singcomp=diag_non_singcomp+1;
 int maskcomp=diag_singcomp+1;

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
  BoxLib::Error("smooth_type invalid");


 for (int isweep=0;isweep<num_sweeps;isweep++) {

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(solnL,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(bxa[mfi.index()] == mfi.validbox());
   int gridno = mfi.index();
   if (gridno>=number_grids)
    BoxLib::Error("gridno invalid");
   const Box& tilegrid=mfi.tilebox();
   const Box& fabgrid=gbox[level][gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   if (gsrb_timing==1) 
    t1 = ParallelDescriptor::second();

   for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

    int ofs=veldir*ncwork;

    FORT_GSRB(
     &isweep,
     &num_sweeps,
     &offdiag_coeff_level,
     &check_for_singular,
     &diag_regularization,
     solnL[mfi].dataPtr(veldir), 
     ARLIM(solnL[mfi].loVect()),ARLIM(solnL[mfi].hiVect()),
     rhsL[mfi].dataPtr(veldir), 
     ARLIM(rhsL[mfi].loVect()), ARLIM(rhsL[mfi].hiVect()),

     work[mfi].dataPtr(diag_non_singcomp+ofs),
     ARLIM(work[mfi].loVect()), ARLIM(work[mfi].hiVect()),
     work[mfi].dataPtr(diag_singcomp+ofs),

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

   } // veldir

   if (gsrb_timing==1) {
    t2 = ParallelDescriptor::second();
    std::cout << "GSRB time, level= " << level << " smooth_type=" <<
     smooth_type << " gridno= " << gridno << " t2-t1=" << t2-t1 << '\n';
   }
  } // mfi
} // omp
  ParallelDescriptor::Barrier();
 } // isweep
}

// y=Ax
void
ABecLaplacian::Fapply (MultiFab& y,
                       MultiFab& x,
                       int level)
{

 bool use_tiling=cfd_tiling;

 const BoxArray& bxa = gbox[level];

 Real offdiag_coeff_level=offdiag_coeff[level];
 if (offdiag_coeff_level>0.0) {
  // do nothing
 } else
  BoxLib::Error("offdiag_coeff_level invalid");

 const MultiFab & work=matCoefficients(level);
 int ncwork=BL_SPACEDIM*3+10;

 int nctest = work.nComp();
 if (nctest!=ncwork*nsolve_bicgstab)
  BoxLib::Error("ncwork invalid");

 int nc = y.nComp();
 if (nc!=nsolve_bicgstab)
  BoxLib::Error("nc bust");
 int ngrow_y=y.nGrow();
 if (ngrow_y!=0) {
  std::cout << "ngrow_y= " << ngrow_y << '\n';
  BoxLib::Error("ngrow y invalid");
 }
 int ngrow_x=x.nGrow();
 if (ngrow_x<1)
  BoxLib::Error("ngrow x invalid");

 int bxleftcomp=0;
 int byleftcomp=bxleftcomp+1;
 int bzleftcomp=byleftcomp+BL_SPACEDIM-2;
 int bxrightcomp=bzleftcomp+1;
 int byrightcomp=bxrightcomp+1;
 int bzrightcomp=byrightcomp+BL_SPACEDIM-2;
 int icbxcomp=bzrightcomp+1;
 int icbycomp=icbxcomp+1;
 int icbzcomp=icbycomp+BL_SPACEDIM-2;

 int diag_non_singcomp=icbzcomp+1;
 int diag_singcomp=diag_non_singcomp+1;
 int maskcomp=diag_singcomp+1;

 int icdiagcomp=maskcomp+1;
 int icdiagrbcomp=icdiagcomp+1;
 int axcomp=icdiagrbcomp+1;
 int solnsavecomp=axcomp+1;
 int rhssavecomp=solnsavecomp+1;
 int redsolncomp=rhssavecomp+1;
 int blacksolncomp=redsolncomp+1;

 if (work.nComp()<=blacksolncomp*nsolve_bicgstab)
  BoxLib::Error("work.nComp() invalid");

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
   BoxLib::Error("gridno invalid");
  const Box& tilegrid=mfi.tilebox();
  const Box& fabgrid=gbox[level][gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  for (int veldir=0;veldir<nsolve_bicgstab;veldir++) {

   int ofs=veldir*ncwork;

   FORT_ADOTX(
    &offdiag_coeff_level,
    &check_for_singular,
    &diag_regularization,
    y[mfi].dataPtr(veldir),
    ARLIM(y[mfi].loVect()), ARLIM(y[mfi].hiVect()),
    x[mfi].dataPtr(veldir),
    ARLIM(x[mfi].loVect()), ARLIM(x[mfi].hiVect()),

    work[mfi].dataPtr(diag_non_singcomp+ofs),
    ARLIM(work[mfi].loVect()),ARLIM(work[mfi].hiVect()),
    work[mfi].dataPtr(diag_singcomp+ofs),

    work[mfi].dataPtr(bxleftcomp+ofs),
    work[mfi].dataPtr(bxrightcomp+ofs),
    work[mfi].dataPtr(byleftcomp+ofs),
    work[mfi].dataPtr(byrightcomp+ofs),
    work[mfi].dataPtr(bzleftcomp+ofs),
    work[mfi].dataPtr(bzrightcomp+ofs),
    tilelo,tilehi,
    fablo,fabhi,&bfact,&bfact_top);

  } // veldir
 } // mfi
} // omp
 ParallelDescriptor::Barrier();
}



void
ABecLaplacian::LP_update (MultiFab& sol,
   Real alpha,MultiFab& y,
   const MultiFab& p,int level) {
    //
    // compute sol=y+alpha p  
    //
 bool use_tiling=cfd_tiling;

 if (level>=numLevels())
  BoxLib::Error("level invalid in LP_update");

 const BoxArray& gboxlev = gbox[level];
 int ncomp = sol.nComp();
 if (ncomp!=nsolve_bicgstab)
  BoxLib::Error("ncomp invalid");

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

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
  } // veldir
 }
} // omp
 ParallelDescriptor::Barrier();
}


void ABecLaplacian::LP_dot(MultiFab& w,const MultiFab& p,
   int level,Real& result) {
 
 bool use_tiling=cfd_tiling;

 if (level>=numLevels())
  BoxLib::Error("level invalid in LP_dot");
 if (level>=gbox.size()) {
  std::cout << "level= " << level << '\n';
  std::cout << "gboxsize= " << gbox.size() << '\n';
  std::cout << "num levels = " << numLevels() << '\n';
  BoxLib::Error("level exceeds gbox size");
 }

 if (thread_class::nthreads<1)
  BoxLib::Error("thread_class::nthreads invalid");

 Array<Real> pw;
 pw.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  pw[tid] = 0.0;
 }
 const BoxArray& gboxlev = gbox[level];
 int ncomp = p.nComp();
 if (ncomp!=nsolve_bicgstab)
  BoxLib::Error("ncomp invalid p");
 if (w.nComp()!=nsolve_bicgstab)
  BoxLib::Error("ncomp invalid w");

 int bfact=bfact_array[level];
 int bfact_top=bfact_array[0];

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
 for (MFIter mfi(w,use_tiling); mfi.isValid(); ++mfi) {
  Real tpw;
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
   BoxLib::Error("tid invalid");

  FORT_CGXDOTY(
   &ncomp,
   &tpw,
   p[mfi].dataPtr(),ARLIM(p[mfi].loVect()), ARLIM(p[mfi].hiVect()),
   w[mfi].dataPtr(),ARLIM(w[mfi].loVect()), ARLIM(w[mfi].hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,&bfact_top);
  pw[tid] += tpw;
 }
} // omp
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  pw[0]+=pw[tid];
 }
 ParallelDescriptor::Barrier();
 ParallelDescriptor::ReduceRealSum(pw[0]);

 result=pw[0];
}


void
ABecLaplacian::project_null_space(MultiFab& rhsL,int level) {

 if (laplacian_solvability==0) {
  if (check_for_singular==1) {
   const MultiFab& ones_mf=onesCoefficients(level);
   MultiFab::Multiply(rhsL,ones_mf,0,0,1,0);
  } else if (check_for_singular==0) {
   // do nothing
  } else
   BoxLib::Error("check_for_singular invalid");

 } else if (laplacian_solvability==1) {

  if (nsolve_bicgstab!=1)
   BoxLib::Error("nsolve_bicgstab invalid");

  if (check_for_singular==1) {
   const MultiFab& ones_mf=onesCoefficients(level);
   MultiFab::Multiply(rhsL,ones_mf,0,0,1,0);

   const BoxArray& temp_ba=ones_mf.boxArray();
   int ncomp=ones_mf.nComp();
   int ngrow=ones_mf.nGrow();
   MultiFab* ones_mf_copy=new MultiFab(temp_ba,ncomp,ngrow,Fab_allocate);
   MultiFab::Copy(*ones_mf_copy,ones_mf,0,0,ncomp,ngrow);

   Real result,domainsum;
   LP_dot(rhsL,ones_mf,level,result);
   LP_dot(*ones_mf_copy,ones_mf,level,domainsum); 
   delete ones_mf_copy;

   double total_cells=temp_ba.d_numPts();
   if (domainsum>total_cells) {
    std::cout << "domainsum= " << domainsum << '\n';
    std::cout << "total_cells= " << total_cells << '\n';
    BoxLib::Error("domainsum too big");
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
    LP_update(rhsL,coef,rhsL,ones_mf,level); 
   } else if (domainsum==0.0) {
    // do nothing
   } else
    BoxLib::Error("domainsum invalid");

   MultiFab::Multiply(rhsL,ones_mf,0,0,1,0);

   if (1==0) {
    std::cout << "check rhsL after projection \n";
    LP_dot(rhsL,ones_mf,level,result);
    std::cout << "level= " << level << '\n';
    std::cout << "result= " << result << '\n';
   }

  } else
   BoxLib::Error("check_for_singular invalid");

 } else
  BoxLib::Error("laplacian solvability incorrect");

} // subroutine project_null_space


// off diagonal sum flag 
void
ABecLaplacian::Fdiagsum(MultiFab&       y,
                       int             level) {

 bool use_tiling=cfd_tiling;

 const BoxArray& bxa = gbox[level];
 const MultiFab& a   = aCoefficients(level);
 const MultiFab& bX  = bCoefficients(0,level);
 const MultiFab& bY  = bCoefficients(1,level);
 const MultiFab& bZ  = bCoefficients(BL_SPACEDIM-1,level);
 int nc = y.nComp();
 if (nc!=nsolve_bicgstab)
  BoxLib::Error("nc bust");

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
  } // veldir
 }
} // omp
    ParallelDescriptor::Barrier();
}
