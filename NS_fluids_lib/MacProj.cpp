#include <winstd.H>

#include <CoordSys.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <NavierStokes.H>
#include <LEVEL_F.H>
#include <GODUNOV_F.H>
#include <NAVIERSTOKES_F.H>
#include <DIFFUSION_F.H>
#include <MACOPERATOR_F.H>
#include <MG_F.H>
#include <CGSolver.H>
#include <MultiGrid.H>

// residual correction form:
// alpha(p-p*)-div beta grad p = -div u*
// dp=p-p0
// alpha(dp+p0-p*)-div beta grad (dp+p0) = -div u*
// alpha dp - div beta grad dp = -div V + alpha(p*-p0)
// V=u*-beta grad p0
//
// prescribed solid algorithm:
// k=0
// pk given
// q=p*-pk
// p=pk
// V=H(u*-beta grad pk)+(1-H)Usolid^k
// repeat
//  alpha dp - div beta grad dp = -div V + alpha q
//  k=k+1
//  q=q-dp
//  V=H(V-beta grad dp)+(1-H)Usolid^k 
//  p=p+dp
// until converge
// 
// at convergence:
//  q=p*-p
//  V=u* - beta grad p
//
// alpha=(1/c)^2 / (rho dt^2)
// alpha=avgdown(alpha)*volume
//
// bx=(1/rho)(1-H_solid)
// bx=avgdown(bx)
// bx*=(area/dx)
//
// RHS= -(1/dt)(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)+
//  alpha_{i}poldhold_{i}
//
// alpha_{i}p_{i}-
//   (bx_{i+1/2} (p_{i+1}-p_{i})-bx_{i-1/2} (p_{i}-p_{i-1})) = RHS
//
// called from: NavierStokes::update_SEM_forcesALL
//              NavierStokes::multiphase_project
void
NavierStokes::allocate_maccoef(int project_option,int nsolve) {

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==13)||  //FSI_material_exists 1st project
     (project_option==11)) { //FSI_material_exists 2nd project

  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
  if (singular_possible!=1)
   BoxLib::Error("singular_possible invalid");

 } else if (project_option==12) {  // extension project

  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
  if (singular_possible!=1)
   BoxLib::Error("singular_possible invalid");

 } else if (project_option==3) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
  if (singular_possible!=0)
   BoxLib::Error("singular_possible invalid");
 } else if ((project_option==2)||     // thermal diffusion
            ((project_option>=100)&&  // species
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
  if (singular_possible!=0)
   BoxLib::Error("singular_possible invalid");
 } else
  BoxLib::Error("project_option invalid60");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;

 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  BoxLib::Error("num_materials_face invalid");
 
 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  BoxLib::Error("num_state_base invalid");

 int bfact=parent->Space_blockingFactor(level);

 Array<int> scomp;
 Array<int> ncomp;
 int state_index;
 int ncomp_check;

 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (ncomp_check!=nsolveMM)
  BoxLib::Error("ncomp_check invalid");

 Real dt_diffuse=dt_slab;

 if (project_option==2) { // temperature diffusion
  if (lower_slab_time==0.0) {
   if (SDC_outer_sweeps!=0)
    BoxLib::Error("SDC_outer_sweeps invalid");
   if (slab_step==0) {
    if (initial_temperature_diffuse_duration>dt_diffuse)
     dt_diffuse=initial_temperature_diffuse_duration; 
   } else if ((slab_step>0)&&(slab_step<ns_time_order)) {
    // do nothing
   } else
    BoxLib::Error("slab_step invalid");
  } else if (lower_slab_time>0.0) {
   // do nothing
  } else
   BoxLib::Error("lower_slab_time invalid");
 } // project_option==2 (temperature diffusion)

 const Real* dx = geom.CellSize();
 const BoxArray& gridparm=grids;
 const Geometry& geomparm=geom;

 mac_op=new ABecLaplacian(gridparm,geomparm,
  bfact,
  level,
  project_option,
  nsolveMM,
  ns_tiling);

 mac_op->laplacian_solvability=solvability_level_flag;
 mac_op->check_for_singular=singular_possible;
 mac_op->diag_regularization=NS_diag_regularization;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "allocate_maccoef level= " << level << 
     " solvability_level_flag= " << solvability_level_flag << 
     " singular_possible= " << singular_possible << '\n';
  }
 }

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(BL_SPACEDIM+1))
  BoxLib::Error("LS_new.nComp()!=nmat*(BL_SPACEDIM+1)");
 
 new_localMF(ALPHACOEF_MF,nsolveMM,0,-1);
 new_localMF(ALPHANOVOLUME_MF,nsolveMM,0,-1);
 localMF[ALPHANOVOLUME_MF]->setVal(0.0,0,nsolveMM,0);

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,200);
 debug_ngrow(FACE_VAR_MF,0,201);

 debug_ngrow(ONES_MF,0,202);
 if (localMF[ONES_MF]->nComp()!=num_materials_face)
  BoxLib::Error("localMF[ONES_MF]->nComp() invalid");

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,202);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,113);
 debug_ngrow(CELL_SOUND_MF,0,113);
 debug_ngrow(CELL_DEN_MF,1,113);
 debug_ngrow(CELL_VISC_MF,1,113);
 debug_ngrow(CELL_DEDT_MF,1,113);
 debug_ngrow(OFF_DIAG_CHECK_MF,0,113);

 if (localMF[OFF_DIAG_CHECK_MF]->nComp()!=nsolveMM)
  BoxLib::Error("localMF[OFF_DIAG_CHECK_MF]->nComp() invalid");
 if (localMF[CELL_DEN_MF]->nComp()!=nmat+1)
  BoxLib::Error("localMF[CELL_DEN_MF]->nComp() invalid");
 if (localMF[CELL_VISC_MF]->nComp()!=nmat+1)
  BoxLib::Error("localMF[CELL_VISC_MF]->nComp() invalid");
 if (localMF[CELL_DEDT_MF]->nComp()!=nmat+1)
  BoxLib::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 if (localMF[CELL_SOUND_MF]->nComp()!=2)
  BoxLib::Error("localMF[CELL_SOUND_MF]->nComp() invalid");

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[ALPHANOVOLUME_MF],use_tiling);mfi.isValid();++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  const Real* xlo = grid_loc[gridno].lo();

  FArrayBox& offdiagcheck=(*localMF[OFF_DIAG_CHECK_MF])[mfi];

   // 1/(c^2 dt^2)  (first component)
   // p^advect      (2nd component)
  FArrayBox& c2fab=(*localMF[CELL_SOUND_MF])[mfi];
  FArrayBox& denfab=(*localMF[CELL_DEN_MF])[mfi];  // inverse of density
   // 1/(rho cv)   (DeDT=cv)
  FArrayBox& DeDTfab=(*localMF[CELL_DEDT_MF])[mfi];  
  FArrayBox& cterm = (*localMF[ALPHANOVOLUME_MF])[mfi];
  FArrayBox& lsfab = LS_new[mfi];
  FArrayBox& mufab=(*localMF[CELL_VISC_MF])[mfi];
  FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];

  int rzflag=0;
  if (CoordSys::IsRZ())
   rzflag=1;
  else if (CoordSys::IsCartesian())
   rzflag=0;
  else if (CoordSys::IsCYLINDRICAL())
   rzflag=3;
  else
   BoxLib::Error("CoordSys bust 51");

   // defined in MACOPERATOR_3D.F90
  FORT_SCALARCOEFF(
    &num_materials_face,
    &nsolve,
    &nsolveMM,
    &nmat,
    xlo,dx,
    offdiagcheck.dataPtr(),
    ARLIM(offdiagcheck.loVect()),ARLIM(offdiagcheck.hiVect()),
    cterm.dataPtr(),ARLIM(cterm.loVect()),ARLIM(cterm.hiVect()),
    c2fab.dataPtr(),ARLIM(c2fab.loVect()),ARLIM(c2fab.hiVect()),
    DeDTfab.dataPtr(),ARLIM(DeDTfab.loVect()),ARLIM(DeDTfab.hiVect()),
    reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
    mufab.dataPtr(),ARLIM(mufab.loVect()),ARLIM(mufab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    &finest_level,
    &visc_coef,&angular_velocity,
    &dt_diffuse,
    &project_option,&rzflag, 
    &solidheat_flag);

 }  // mfi
} // omp
 ParallelDescriptor::Barrier();

  // alpha T - div beta grad T = f
 if (project_option==2) {
  if (is_phasechange==1) {
    // alphanovolume=(rho cv/(dt*fact))+(1/vol) sum_face Aface k_m/(theta dx)
   int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
   int adjust_temperature=0; 
   int GFM_flag=0;
   for (int im=0;im<2*nten;im++) {
    if (latent_heat[im]!=0.0)
     if ((freezing_model[im]==0)||
         (freezing_model[im]==5))
      GFM_flag=1;
   }
   if (GFM_flag==1)
    stefan_solver_init(localMF[ALPHANOVOLUME_MF],adjust_temperature);
  } else if (is_phasechange==0) {
   // do nothing
  } else
   BoxLib::Error("is_phasechange invalid");
 }

 avgDown_localMF(ALPHANOVOLUME_MF,0,nsolveMM,0);
 Copy_localMF(ALPHACOEF_MF,ALPHANOVOLUME_MF,0,0,nsolveMM,0);
  // dest,source,scomp,dcomp,ncomp,ngrow
 for (int veldir=0;veldir<nsolveMM;veldir++)
  Mult_localMF(ALPHACOEF_MF,VOLUME_MF,0,veldir,1,0);

 for (int dir=0;dir<BL_SPACEDIM;dir++) {

   if (localMF[AREA_MF+dir]->boxArray()!=
       localMF[FACE_WEIGHT_MF+dir]->boxArray())
    BoxLib::Error("face_weight_stable boxarrays do not match");

   if (localMF[AREA_MF+dir]->boxArray()!=
       localMF[FACE_VAR_MF+dir]->boxArray())
    BoxLib::Error("face_var boxarrays do not match");

   new_localMF(BXCOEFNOAREA_MF+dir,nsolveMM,0,dir);
   new_localMF(BXCOEF_MF+dir,nsolveMM,0,dir);
   localMF[BXCOEFNOAREA_MF+dir]->setVal(1.0,0,nsolveMM,0);

   if (localMF[FACE_WEIGHT_MF+dir]->nComp()!=nsolveMM_FACE) 
    BoxLib::Error("localMF[FACE_WEIGHT_MF+dir]->nComp() invalid");
   if (localMF[BXCOEFNOAREA_MF+dir]->nComp()!=nsolveMM) 
    BoxLib::Error("localMF[BXCOEFNOAREA_MF+dir]->nComp() invalid");

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[ALPHACOEF_MF],use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox & bxfab=(*localMF[BXCOEFNOAREA_MF+dir])[mfi];
    FArrayBox & facefab=(*localMF[FACE_WEIGHT_MF+dir])[mfi];

     // BXCOEFNOAREA *= (facewtL + facewtR)/2
    FORT_MULT_FACEWT(
     &num_materials_face,
     &nsolve,
     &nsolveMM,
     &nsolveMM_FACE,
     bxfab.dataPtr(),ARLIM(bxfab.loVect()),ARLIM(bxfab.hiVect()),
     facefab.dataPtr(),ARLIM(facefab.loVect()),ARLIM(facefab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     xlo,dx,&dir);
   } // mfi
} // omp
   ParallelDescriptor::Barrier();

   int ncomp_edge=-1;
   int scomp=0;
   int ncomp_mf=1;
   avgDownEdge_localMF(BXCOEFNOAREA_MF,scomp,ncomp_edge,dir,ncomp_mf,0,17);
   Copy_localMF(BXCOEF_MF+dir,BXCOEFNOAREA_MF+dir,0,0,nsolveMM,0);
    // dest,source,scomp,dcomp,ncomp,ngrow
   for (int veldir=0;veldir<nsolveMM;veldir++)
    Mult_localMF(BXCOEF_MF+dir,AREA_MF+dir,0,veldir,1,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[ALPHACOEF_MF],use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox & bxfab=(*localMF[BXCOEF_MF+dir])[mfi];
    FORT_DIVIDEDX(
     &num_materials_face,
     &nsolve,
     &nsolveMM,
     bxfab.dataPtr(),ARLIM(bxfab.loVect()),ARLIM(bxfab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     xlo,dx,&dir);
   } // mfi
} // omp
   ParallelDescriptor::Barrier();

   mac_op->bCoefficients(*localMF[BXCOEF_MF+dir],dir);
 }  // dir=0...sdim-1

 mac_op->aCoefficients(*localMF[ALPHACOEF_MF]);

 OFFDIAG_NONSING_LEVEL=max_face_wt[0][1];
 if (OFFDIAG_NONSING_LEVEL>0.0) {
  // do nothing
 } else
  BoxLib::Error("OFFDIAG_NONSING_LEVEL invalid");

 if (BL_SPACEDIM==2) {
  // do nothing
 } else if (BL_SPACEDIM==3) {
  OFFDIAG_NONSING_LEVEL*=dx[0];  // ~ area/dx
 } else
  BoxLib::Error("dimension bust");

 mac_op->non_sing_coefficients(OFFDIAG_NONSING_LEVEL);

  // generateCoefficients calls buildMatrix on the finest mglib level (lev=0)
  // buildMatrix calls FORT_BUILDMAT
 mac_op->generateCoefficients();

 new_localMF(DIAG_NON_SING_MF,nsolveMM,0,-1);
 new_localMF(DIAG_SING_MF,nsolveMM,0,-1);
 new_localMF(MASK_DIV_RESIDUAL_MF,1,0,-1);
 new_localMF(MASK_RESIDUAL_MF,1,0,-1);
 localMF[DIAG_NON_SING_MF]->setVal(0.0,0,nsolveMM,0);
 localMF[DIAG_SING_MF]->setVal(0.0,0,nsolveMM,0);
 localMF[MASK_DIV_RESIDUAL_MF]->setVal(0.0,0,1,0);
 localMF[MASK_RESIDUAL_MF]->setVal(0.0,0,1,0);

 int mm_areafrac_index=FACE_VAR_MF;
 int mm_cell_areafrac_index=SLOPE_RECON_MF;
 if (num_materials_face==nmat) {
  mm_areafrac_index=FACEFRAC_SOLVE_MM_MF;
  mm_cell_areafrac_index=CELLFRAC_MM_MF;
 } else if (num_materials_face==1) {
  // do nothing
 } else
  BoxLib::Error("num_materials_face invalid");

 // (ml,mr,2) frac_pair(ml,mr), dist_pair(ml,mr)  
 int nfacefrac=nmat*nmat*2; 
 // im_inside,im_outside,3+sdim -->
 //   area, dist_to_line, dist, line normal.
 int ncellfrac=nmat*nmat*(3+BL_SPACEDIM);

 for (int dir=0;dir<BL_SPACEDIM;dir++)
  debug_ngrow(mm_areafrac_index+dir,0,111);
 debug_ngrow(mm_cell_areafrac_index,0,113);


#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[DIAG_SING_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  FArrayBox& ones_fab=(*localMF[ONES_MF])[mfi];

   // mask=tag if not covered by level+1 or outside the domain.
  FArrayBox& maskcov = (*localMF[MASKCOEF_MF])[mfi];
  FArrayBox& alphafab = (*localMF[ALPHACOEF_MF])[mfi];

  FArrayBox& offdiagcheck=(*localMF[OFF_DIAG_CHECK_MF])[mfi];
  FArrayBox& diagnonsingfab = (*localMF[DIAG_NON_SING_MF])[mfi];
  FArrayBox& diagsingfab = (*localMF[DIAG_SING_MF])[mfi];

  FArrayBox& maskdivresidfab = (*localMF[MASK_DIV_RESIDUAL_MF])[mfi];
  FArrayBox& maskresidfab = (*localMF[MASK_RESIDUAL_MF])[mfi];
  FArrayBox& bxfab = (*localMF[BXCOEF_MF])[mfi];
  FArrayBox& byfab = (*localMF[BXCOEF_MF+1])[mfi];
  FArrayBox& bzfab = (*localMF[BXCOEF_MF+BL_SPACEDIM-1])[mfi];
  FArrayBox& fwtxfab = (*localMF[FACE_WEIGHT_MF])[mfi];
  FArrayBox& fwtyfab = (*localMF[FACE_WEIGHT_MF+1])[mfi];
  FArrayBox& fwtzfab = (*localMF[FACE_WEIGHT_MF+BL_SPACEDIM-1])[mfi];

  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];  
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];  
  FArrayBox& zface=(*localMF[FACE_VAR_MF+BL_SPACEDIM-1])[mfi];  

  FArrayBox& xfacemm=(*localMF[mm_areafrac_index])[mfi];  
  FArrayBox& yfacemm=(*localMF[mm_areafrac_index+1])[mfi];  
  FArrayBox& zfacemm=(*localMF[mm_areafrac_index+BL_SPACEDIM-1])[mfi];  
  FArrayBox& cellfracmm=(*localMF[mm_cell_areafrac_index])[mfi];  

  Array<int> bc;
  getBCArray_list(bc,state_index,gridno,scomp,ncomp);
  if (bc.size()!=nsolveMM*BL_SPACEDIM*2)
   BoxLib::Error("bc.size() invalid");

   // FORT_NSGENERATE is in MACOPERATOR_3D.F90
   // initializes DIAG_SING
  FORT_NSGENERATE(
    &num_materials_face,
    &level,
    &finest_level,
    &nsolve,
    &nsolveMM,
    &nsolveMM_FACE,
    &nfacefrac,
    &ncellfrac,
    &nmat,
    &project_option,
    &ncphys,
    cellfracmm.dataPtr(),
    ARLIM(cellfracmm.loVect()),ARLIM(cellfracmm.hiVect()),
    xfacemm.dataPtr(),ARLIM(xfacemm.loVect()),ARLIM(xfacemm.hiVect()),
    yfacemm.dataPtr(),ARLIM(yfacemm.loVect()),ARLIM(yfacemm.hiVect()),
    zfacemm.dataPtr(),ARLIM(zfacemm.loVect()),ARLIM(zfacemm.hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
    yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
    zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
    ones_fab.dataPtr(),ARLIM(ones_fab.loVect()),ARLIM(ones_fab.hiVect()),
    maskcov.dataPtr(),ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    alphafab.dataPtr(),ARLIM(alphafab.loVect()),ARLIM(alphafab.hiVect()),
    offdiagcheck.dataPtr(),
    ARLIM(offdiagcheck.loVect()),ARLIM(offdiagcheck.hiVect()),

    diagnonsingfab.dataPtr(),
    ARLIM(diagnonsingfab.loVect()),ARLIM(diagnonsingfab.hiVect()),
    diagsingfab.dataPtr(),
    ARLIM(diagsingfab.loVect()),ARLIM(diagsingfab.hiVect()),

    maskdivresidfab.dataPtr(),
    ARLIM(maskdivresidfab.loVect()),ARLIM(maskdivresidfab.hiVect()),
    maskresidfab.dataPtr(),
    ARLIM(maskresidfab.loVect()),ARLIM(maskresidfab.hiVect()),
    bxfab.dataPtr(),ARLIM(bxfab.loVect()),ARLIM(bxfab.hiVect()),
    byfab.dataPtr(),ARLIM(byfab.loVect()),ARLIM(byfab.hiVect()),
    bzfab.dataPtr(),ARLIM(bzfab.loVect()),ARLIM(bzfab.hiVect()),
    fwtxfab.dataPtr(),ARLIM(fwtxfab.loVect()),ARLIM(fwtxfab.hiVect()),
    fwtyfab.dataPtr(),ARLIM(fwtyfab.loVect()),ARLIM(fwtyfab.hiVect()),
    fwtzfab.dataPtr(),ARLIM(fwtzfab.loVect()),ARLIM(fwtzfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    bc.dataPtr());
 }// mfi
} // omp
 ParallelDescriptor::Barrier();

}  // subroutine allocate_maccoef

void 
NavierStokes::restore_active_pressure(int save_mf) {

 int finest_level=parent->finestLevel();

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 int nsolveMM=1;

 bool use_tiling=ns_tiling;
 int bfact=parent->Space_blockingFactor(level);

 debug_ngrow(save_mf,0,113);
 if (localMF[save_mf]->nComp()!=1)
  BoxLib::Error("localMF[save_mf]->nComp() invalid");
 debug_ngrow(OFF_DIAG_CHECK_MF,0,113);
 if (localMF[OFF_DIAG_CHECK_MF]->nComp()!=nsolveMM)
  BoxLib::Error("localMF[OFF_DIAG_CHECK_MF]->nComp() invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int pcomp=num_materials_vel*BL_SPACEDIM;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[OFF_DIAG_CHECK_MF],use_tiling);mfi.isValid();++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  FArrayBox& offdiagcheck=(*localMF[OFF_DIAG_CHECK_MF])[mfi];
  FArrayBox& savepres=(*localMF[save_mf])[mfi];
  FArrayBox& newpres=S_new[mfi];

   // defined in MACOPERATOR_3D.F90
  FORT_RESTORE_PRES(
   offdiagcheck.dataPtr(),
   ARLIM(offdiagcheck.loVect()),ARLIM(offdiagcheck.hiVect()),
   savepres.dataPtr(),
   ARLIM(savepres.loVect()),ARLIM(savepres.hiVect()),
   newpres.dataPtr(pcomp),
   ARLIM(newpres.loVect()),ARLIM(newpres.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &finest_level);

 }  // mfi
} // omp
 ParallelDescriptor::Barrier();

} // end subroutine restore_active_pressure

void
NavierStokes::deallocate_maccoef(int project_option) {

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists 2nd project
     (project_option==13)|| //FSI_material_exists 1st project
     (project_option==12)||
     (project_option==2)||
     (project_option==3)||
     ((project_option>=100)&&
      (project_option<100+num_species_var))) {
  delete mac_op;
 } else
  BoxLib::Error("project_option invalid deallocate_maccoef");

 delete_localMF(DIAG_NON_SING_MF,1);
 delete_localMF(DIAG_SING_MF,1);
 delete_localMF(MASK_DIV_RESIDUAL_MF,1);
 delete_localMF(MASK_RESIDUAL_MF,1);
 delete_localMF(ALPHANOVOLUME_MF,1);
 delete_localMF(ALPHACOEF_MF,1);
 delete_localMF(BXCOEFNOAREA_MF,BL_SPACEDIM);
 delete_localMF(BXCOEF_MF,BL_SPACEDIM);

} // end subroutine deallocate_maccoef


// interpolates coarse data and ADDS it to the fine data.
// interpolates where cdiagsing>0
void
NavierStokes::AllinterpScalarMAC(
  MultiFab* coarsedata,MultiFab* finedata,
  BoxArray& cgridscen,BoxArray& fgridscen,
  MultiFab* cdiagsing,MultiFab* fdiagsing,
  int nsolve,int project_option) {

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==13)|| //FSI_material_exists 1st project
     (project_option==11)|| //FSI_material_exists 2nd project
     (project_option==12)) {

  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");

 } else if (project_option==3) {  // viscosity

  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");

 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid61");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face; 

 int finest_level=parent->finestLevel();
 if ((level>finest_level)||(level<1))
  BoxLib::Error("level invalid AllinterpScalarMAC");

 if ((coarsedata->nComp()!=nsolveMM)||
     (finedata->nComp()!=nsolveMM)||
     (cdiagsing->nComp()!=1)||
     (fdiagsing->nComp()!=1))
  BoxLib::Error("invalid ncomp");

 BoxArray crse_cen_fine_BA(fgridscen.size());
 for (int i = 0; i < fgridscen.size(); ++i) {
  crse_cen_fine_BA.set(i,BoxLib::coarsen(fgridscen[i],2));
 }

 MultiFab& S_fine = *finedata;
 MultiFab& pcoarse = *coarsedata;
 const BoxArray& fgrids=S_fine.boxArray();

 BoxArray crse_S_fine_BA(fgrids.size());

 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,BoxLib::coarsen(fgrids[i],2));
 }

 MultiFab crse_S_fine(crse_S_fine_BA,nsolveMM,0);
 crse_S_fine.copy(pcoarse,0,0,nsolveMM);

 MultiFab crse_diagsing_fine(crse_S_fine_BA,1,0);
 crse_diagsing_fine.copy(*cdiagsing,0,0,1);

 int bfact_f=parent->Space_blockingFactor(level);
 int bfact=parent->Space_blockingFactor(level-1);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgridscen[mfi.index()] == mfi.validbox());
  const int i = mfi.index();

  FArrayBox& crse_fab = crse_S_fine[mfi];
  const Box& cbox = crse_cen_fine_BA[i];

  FArrayBox& fine_fab = S_fine[mfi];

  FArrayBox& cdiagfab=crse_diagsing_fine[mfi];
  FArrayBox& fdiagfab=(*fdiagsing)[mfi];
  for (int veldir=0;veldir<nsolveMM;veldir++) {
   FORT_INTERPMAC(
    &bfact,&bfact_f,
    fine_fab.dataPtr(veldir),
    ARLIM(fine_fab.loVect()),ARLIM(fine_fab.hiVect()),
    crse_fab.dataPtr(veldir),
    ARLIM(crse_fab.loVect()),ARLIM(crse_fab.hiVect()),
    cbox.loVect(),cbox.hiVect(),
    cdiagfab.dataPtr(0),
    ARLIM(cdiagfab.loVect()),ARLIM(cdiagfab.hiVect()),
    fdiagfab.dataPtr(0),
    ARLIM(fdiagfab.loVect()),ARLIM(fdiagfab.hiVect()));  
  } // veldir
 }   // mfi
} //omp
 ParallelDescriptor::Barrier();

}  // subroutine AllinterpScalarMAC

void
NavierStokes::interpScalarMAC(MultiFab* coarsedata,MultiFab* finedata,
 int nsolve,int project_option) {

  int finest_level=parent->finestLevel();
  if ((level>finest_level)||(level<1))
   BoxLib::Error("level invalid interpScalarMAC");

  NavierStokes& coarse_lev = getLevel(level-1);
  BoxArray& fgridscen=grids;
  BoxArray& cgridscen=coarse_lev.grids;
  AllinterpScalarMAC(
    coarsedata,finedata,cgridscen,fgridscen,
    coarse_lev.localMF[MASK_RESIDUAL_MF],
    localMF[MASK_RESIDUAL_MF],
    nsolve,project_option);

}  // subroutine interpScalarMAC

void
NavierStokes::Allaverage(MultiFab* coarsedata,MultiFab* finedata,
  BoxArray& cgridscen,BoxArray& fgridscen,int iaverage,
  int scomp,int dcomp) {
 

 int finest_level=parent->finestLevel();
 if (level>=finest_level)
  BoxLib::Error("level invalid Allaverage");
 int bfact_fine=parent->Space_blockingFactor(level+1);
 int bfact_coarse=parent->Space_blockingFactor(level);

 BoxArray crse_cen_fine_BA(fgridscen.size());
 for (int i = 0; i < fgridscen.size(); ++i) {
  crse_cen_fine_BA.set(i,BoxLib::coarsen(fgridscen[i],2));
 }

 MultiFab& S_crse = *coarsedata;
 MultiFab& S_fine = *finedata;
 const BoxArray& fgrids=S_fine.boxArray();

 BoxArray crse_S_fine_BA(fgrids.size());

 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,BoxLib::coarsen(fgrids[i],2));
 }

 MultiFab crse_S_fine(crse_S_fine_BA,1,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const int i = mfi.index();

  FArrayBox& crse_fab = crse_S_fine[mfi];
  const Box& cbox = crse_cen_fine_BA[i];

  FArrayBox& fine_fab = S_fine[mfi];

  FORT_AVERAGE(
    crse_fab.dataPtr(0),
    ARLIM(crse_fab.loVect()),ARLIM(crse_fab.hiVect()),
    fine_fab.dataPtr(scomp),
    ARLIM(fine_fab.loVect()),ARLIM(fine_fab.hiVect()),
    cbox.loVect(),cbox.hiVect(),
    &iaverage, 
    &bfact_coarse,
    &bfact_fine,&bfact_fine);
 }   // mfi
} //omp
 ParallelDescriptor::Barrier();

// src,scomp,dcomp,ncomp
 S_crse.copy(crse_S_fine,0,dcomp,1);

} // Allaverage

void
NavierStokes::averageRhs(int idx_MF,int nsolve,int project_option) {

  int nmat=num_materials;

  int num_materials_face=num_materials_vel;
  if ((project_option==0)||
      (project_option==1)||
      (project_option==10)||
      (project_option==11)|| //FSI_material_exists 2nd project
      (project_option==13)|| //FSI_material_exists 1st project
      (project_option==12)||
      (project_option==3)) {  // viscosity
   if (num_materials_face!=1)
    BoxLib::Error("num_materials_face invalid");
  } else if ((project_option==2)||  // thermal diffusion
             ((project_option>=100)&&
              (project_option<100+num_species_var))) {
   num_materials_face=num_materials_scalar_solve;
  } else
   BoxLib::Error("project_option invalid62");

  if (num_materials_vel!=1)
   BoxLib::Error("num_materials_vel invalid");

  if ((num_materials_face!=1)&&
      (num_materials_face!=nmat))
   BoxLib::Error("num_materials_face invalid");

  int nsolveMM=nsolve*num_materials_face;

  if ((nsolve!=1)&&(nsolve!=BL_SPACEDIM))
   BoxLib::Error("nsolve invalid25");
  if (localMF[idx_MF]->nComp()!=nsolveMM)
   BoxLib::Error("nsolve invalid26");

  int finest_level=parent->finestLevel();
  int iavg=0;
  if (level>=finest_level)
   BoxLib::Error("level invalid averageRhs");

  NavierStokes& fine_lev = getLevel(level+1);
  BoxArray& fgridscen=fine_lev.grids;
  BoxArray& cgridscen=grids;

  for (int veldir=0;veldir<nsolveMM;veldir++) {
   Allaverage(localMF[idx_MF],fine_lev.localMF[idx_MF],
     cgridscen,fgridscen,iavg,veldir,veldir); 
  }

}  // averageRhs

// x^{new}=x^{old}+D^{-1}(b-Ax^{old})
// on finest level,
// x^{new}=x^{old}+ILU^{-1}(b-Ax^{old})
void NavierStokes::JacobiALL(
 int idx_resid,int idx_rhs,int idx_xnew,
 int project_option,int nsolve) {

 int finest_level=parent->finestLevel();
 int ilev=finest_level;

 NavierStokes& ns_level = getLevel(ilev);
 ns_level.applyBC_LEVEL(project_option,idx_xnew,nsolve);
 int apply_lev=0;

  // the smoother uses A_LOW: e.g.
  // z^{k+1}=z^{k}+D_LOW^{-1}(r-A_LOW z^{k})
 
 ns_level.mac_op->Fsmooth(
  *ns_level.localMF[idx_xnew],
  *ns_level.localMF[idx_rhs],
  apply_lev,smooth_type);

 for (ilev=finest_level-1;ilev>=0;ilev--) {
  NavierStokes& ns_level = getLevel(ilev);
  ns_level.DiagInverse(
   ns_level.localMF[idx_resid],
   ns_level.localMF[idx_xnew],nsolve,project_option);
 } 

}  // subroutine JacobiALL

void NavierStokes::DiagInverse(
  MultiFab* resid,MultiFab* xnew,int nsolve,int project_option) {
 
 bool use_tiling=ns_tiling;

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,845);

 debug_ngrow(FACE_VAR_MF,0,205);

 int finest_level=parent->finestLevel();
 if (level > finest_level)
  BoxLib::Error("level too big");

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists 2nd project
     (project_option==13)|| //FSI_material_exists 1st project
     (project_option==12)||
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid63");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;
 if ((nsolve!=1)&&(nsolve!=BL_SPACEDIM))
  BoxLib::Error("nsolve invalid250");
 if (resid->nComp()!=nsolveMM)
  BoxLib::Error("resid->nComp() invalid");
 if (xnew->nComp()!=nsolveMM)
  BoxLib::Error("xnew->nComp() invalid");

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*resid,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& residfab = (*resid)[mfi];
  FArrayBox& xnewfab  = (*xnew)[mfi];
  FArrayBox& mfab = (*localMF[MASKCOEF_MF])[mfi];
  FArrayBox& xoldfab = (*xnew)[mfi];

  FArrayBox& diagnonsingfab = (*localMF[DIAG_NON_SING_MF])[mfi];
  if (diagnonsingfab.nComp()!=nsolveMM)
   BoxLib::Error("diagnonsingfab.nComp() invalid");
  FArrayBox& diagsingfab = (*localMF[DIAG_SING_MF])[mfi];
  if (diagsingfab.nComp()!=nsolveMM)
   BoxLib::Error("diagsingfab.nComp() invalid");

  for (int veldir=0;veldir<nsolveMM;veldir++) {
    // FORT_DIAGINV is in NAVIERSTOKES_3D.F90
   FORT_DIAGINV(
    &singular_possible,
    &OFFDIAG_NONSING_LEVEL,
    &NS_diag_regularization,
    diagnonsingfab.dataPtr(veldir),
    ARLIM(diagnonsingfab.loVect()),ARLIM(diagnonsingfab.hiVect()),
    diagsingfab.dataPtr(veldir),
    ARLIM(diagsingfab.loVect()),ARLIM(diagsingfab.hiVect()),
    residfab.dataPtr(veldir),ARLIM(residfab.loVect()),ARLIM(residfab.hiVect()),
    xnewfab.dataPtr(veldir),ARLIM(xnewfab.loVect()),ARLIM(xnewfab.hiVect()),
    xoldfab.dataPtr(veldir),ARLIM(xoldfab.loVect()),ARLIM(xoldfab.hiVect()),
    mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact);
  } // veldir
 } // mfi
} // omp
 ParallelDescriptor::Barrier();

} // subroutine DiagInverse


// RHS=p^n*alpha-(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+... )/dt+diffusionRHS
// solving p*alpha-(bx_{i+1/2}(p_{i+1}-p_{i})+...)=RHS
// resid=rhs-Aphi
void NavierStokes::residALL(
  int project_option,
  int idx_rhs,int idx_resid,int idx_phi,int nsolve) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  BoxLib::Error("level invalid residALL");

 applyALL(project_option,idx_phi,idx_resid,nsolve);

  // resid_array=rhs_array-resid_array
 Real beta=-1.0;
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.levelCombine(
   project_option,
   ns_level.localMF[idx_rhs],
   ns_level.localMF[idx_resid],
   ns_level.localMF[idx_resid],
   beta,nsolve);
  if (ilev<finest_level)
   ns_level.averageRhs(idx_resid,nsolve,project_option);
 }
  // in: residALL
 project_right_hand_side(idx_resid,project_option);
} // subroutine residALL

void NavierStokes::applyALL(
  int project_option,
  int idx_phi,int idx_Aphi,int nsolve) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  BoxLib::Error("level invalid applyALL");
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  ns_level.avgDown_localMF(idx_phi,0,nsolve,1);
  ns_level.applyBC_LEVEL(project_option,idx_phi,nsolve);

   // gradpedge=-dt W grad p
  int homflag=1;
  int energyflag=0;
  int simple_AMR_BC_flag=0;
  ns_level.apply_pressure_grad(
   simple_AMR_BC_flag,
   homflag,
   energyflag,
   GRADPEDGE_MF,
   idx_phi,
   project_option,nsolve);

  if (ilev<finest_level) {
   int ncomp_edge=-1;
   int scomp=0;
   ns_level.avgDownEdge_localMF(GRADPEDGE_MF,scomp,ncomp_edge,0,
		   BL_SPACEDIM,1,18);
  }

  MultiFab* mdot_local=ns_level.localMF[DIFFUSIONRHS_MF];

  if (1==0) {
   Real nrm=mdot_local->norm2(0);
   std::cout << "ilev= " << ilev << "mdot norm2 = " << nrm << '\n';
  }

// Aphi=phi*alpha+(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)/dt
  homflag=1;
  ns_level.apply_div(
   project_option,homflag,
   idx_phi,
   ns_level.localMF[idx_Aphi],
   mdot_local,
   GRADPEDGE_MF,
   nsolve);
 } // ilev=finest_level ... level
} // subroutine applyALL

// called from JacobiALL, applyALL, applyGradALL
// applyALL is called from residALL
// JacobiALL is called from jacobi_cycles
void NavierStokes::applyBC_LEVEL(int project_option,int idx_phi,int nsolve) {

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists 2nd project
     (project_option==13)|| //FSI_material_exists 1st project
     (project_option==12)||
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid64");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;

 if (override_bc_to_homogeneous!=1)
  BoxLib::Error("override_bc_to_homogeneous invalid");

 if (localMF[idx_phi]->nComp()!=nsolveMM)
  BoxLib::Error("invalid ncomp");
 if (localMF[idx_phi]->nGrow()!=1)
  BoxLib::Error("invalid ngrow");

 localMF[idx_phi]->setBndry(0.0);

 int bfact=parent->Space_blockingFactor(level);
 if ((bfact<1)||(bfact>64))
  BoxLib::Error("bfact out of range");

 int state_index;
 Array<int> scomp;
 Array<int> ncomp;
 int ncomp_check;
 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp, 
  ncomp_check);
 if (ncomp_check!=nsolveMM)
  BoxLib::Error("nsolveMM invalid 898");

 Array<int> scompBC_map;
 scompBC_map.resize(nsolveMM);

 int dcomp=0; 
 for (int ilist=0;ilist<scomp.size();ilist++) {
  for (int nc=0;nc<ncomp[ilist];nc++) {
   scompBC_map[dcomp]=scomp[ilist]+nc;
   dcomp++;
  }
 }
 if (dcomp!=nsolveMM)
  BoxLib::Error("dcomp invalid"); 

 MultiFab* cmf;

 if (level>0) {
  NavierStokes& ns_coarse=getLevel(level-1);
  ns_coarse.localMF[idx_phi]->setBndry(0.0);
  ns_coarse.localMF[idx_phi]->FillBoundary(ns_coarse.geom.periodicity());
  cmf=ns_coarse.localMF[idx_phi];
 } else if (level==0) {
  cmf=localMF[idx_phi];
 } else
  BoxLib::Error("level invalid applyBC_LEVEL");

 InterpBorders(
   *cmf,
   *localMF[idx_phi],
   cur_time_slab,
   state_index,
   0,  // scomp=0
   scompBC_map,
   nsolveMM);

}  // subroutine applyBC_LEVEL

//homflag==1 => going down the V-cycle
//homflag==0 => going up the V-cycle
void NavierStokes::applyBC_MGLEVEL(int idx_phi,
 MultiFab* pbdry,int homflag,int nsolve,int project_option) {

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists 2nd project
     (project_option==13)|| //FSI_material_exists 1st project
     (project_option==12)||
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid65");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;

 if (override_bc_to_homogeneous!=1)
  BoxLib::Error("override_bc_to_homogeneous invalid");

 int bfact=parent->Space_blockingFactor(level);
 if ((bfact<1)||(bfact>64))
  BoxLib::Error("bfact out of range");

 if ((homflag!=1)&&(homflag!=0))
  BoxLib::Error("homflag invalid");

 if (localMF[idx_phi]->nComp()!=nsolveMM)
  BoxLib::Error("invalid ncomp");
 if (localMF[idx_phi]->nGrow()!=1)
  BoxLib::Error("invalid ngrow");

  // down the V-cycle
 if (homflag==1) 
  setVal_localMF(idx_phi,0.0,0,nsolveMM,1);

 localMF[idx_phi]->setBndry(0.0);

 int state_index;
 Array<int> scomp;
 Array<int> ncomp;
 int ncomp_check;
 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 if (ncomp_check!=nsolveMM)
  BoxLib::Error("nsolveMM invalid 976");

  // up the V-cycle
 if ((homflag==0)&&(level>0)) {
  NavierStokes& ns_coarse=getLevel(level-1);
  ns_coarse.localMF[idx_phi]->setBndry(0.0);
  ns_coarse.localMF[idx_phi]->FillBoundary(ns_coarse.geom.periodicity());

  Array<int> scompBC_map;
  scompBC_map.resize(nsolveMM);
  int dcomp=0;
  for (int ilist=0;ilist<scomp.size();ilist++) {
   for (int nc=0;nc<ncomp[ilist];nc++) {
    scompBC_map[dcomp]=scomp[ilist]+nc;
    dcomp++;
   }
  }
  if (dcomp!=nsolveMM)
   BoxLib::Error("dcomp invalid");

  InterpBorders(*ns_coarse.localMF[idx_phi],
    *localMF[idx_phi],
    cur_time_slab,
    state_index,
    0,  // scomp=0
    scompBC_map,
    nsolveMM);
  interpScalarMAC(ns_coarse.localMF[idx_phi],
                  localMF[idx_phi],nsolve,project_option);
 }
 
 MultiFab::Copy(*pbdry,*localMF[idx_phi],0,0,nsolveMM,1);

} // applyBC_MGLEVEL


   // gradpedge=-dt W grad p
void NavierStokes::applyGradALL(
  int project_option,int idx_phi,int nsolve) {

 int homflag=1;
 int energyflag=0;
 int finest_level=parent->finestLevel();
 if (level!=0)
  BoxLib::Error("level invalid applyGradALL");
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.applyBC_LEVEL(project_option,idx_phi,nsolve);

   // gradpedge=-dt W grad p
  int simple_AMR_BC_flag=0;
  ns_level.apply_pressure_grad(
   simple_AMR_BC_flag,
   homflag,
   energyflag,
   GRADPEDGE_MF,
   idx_phi,
   project_option,nsolve);
  if (ilev<finest_level) {
   int ncomp_edge=-1;
   int scomp=0;
   ns_level.avgDownEdge_localMF(GRADPEDGE_MF,scomp,ncomp_edge,0,
		   BL_SPACEDIM,1,19);
  }
 } // ilev
} // subroutine applyGradALL

// homflag=0 =>
// RHS=p*alpha-(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2})/dt+diffusionRHS
// homflag=1 =>
// RHS=p*alpha+(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)/dt
// homflag=2 =>
// RHS=-p*alpha-(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2})/dt+diffusionRHS
// aka: residmf=-idx_phi*alpha-(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2})/dt+
//              rhsmf 
// homflag=3 (same as homflag=0) =>
// rhsmf=p*alpha-(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2})/dt+idx_rhs
//
// homflag=4: rhsmf=div u
//


void NavierStokes::apply_div(
  int project_option,int homflag,
  int idx_phi,
  MultiFab* rhsmf, 
  MultiFab* diffusionRHScell,
  int idx_gphi,
  int nsolve) {

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists 2nd project
     (project_option==13)|| //FSI_material_exists 1st project
     (project_option==12)||
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid66");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int finest_level=parent->finestLevel();

 int nsolveMM=nsolve*num_materials_face;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  BoxLib::Error("num_materials_face invalid");
 
 bool use_tiling=ns_tiling;

 Array<int> scomp;
 Array<int> ncomp;
 int state_index;
 int ncomp_check;

 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (ncomp_check!=nsolveMM)
  BoxLib::Error("ncomp_check invalid");

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>=nmat))
  BoxLib::Error("nparts invalid");
 Array<int> im_solid_map_null;
 im_solid_map_null.resize(1);

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=nmat-1)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  BoxLib::Error("nparts invalid");

 resize_FSI_GHOST_MF(1);
 if (localMF[FSI_GHOST_MF]->nGrow()!=1)
  BoxLib::Error("localMF[FSI_GHOST_MF]->nGrow()!=1");
 if (localMF[FSI_GHOST_MF]->nComp()!=nparts_def*BL_SPACEDIM)
  BoxLib::Error("localMF[FSI_GHOST_MF]->nComp()!=nparts_def*BL_SPACEDIM");

 const Real* dx = geom.CellSize();

 if ((homflag!=0)&&
     (homflag!=1)&&
     (homflag!=2)&&
     (homflag!=3)&&
     (homflag!=4))
  BoxLib::Error("homflag invalid");

 int fluxvel_index=0;
 int fluxden_index=BL_SPACEDIM;

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,250);
 debug_ngrow(FACE_VAR_MF,0,251);
 debug_ngrow(FSI_GHOST_MF,1,252);
 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 debug_ngrow(MASKCOEF_MF,1,253); // maskcoef=1 if not covered by finer lev.
 debug_ngrow(MASK_NBR_MF,1,253); // mask_nbr=1 at fine-fine bc.

 debug_ngrow(DOTMASK_MF,0,253); 

 if (diffusionRHScell->nGrow()<0)
  BoxLib::Error("diffusionRHScell invalid");
 if ((localMF[idx_phi]->nComp()!=nsolveMM)||
     (rhsmf->nComp()!=nsolveMM)||
     (localMF[DOTMASK_MF]->nComp()!=num_materials_face)||
     (diffusionRHScell->nComp()!=nsolveMM)||
     (localMF[idx_gphi]->nComp()!=nsolveMM_FACE)||
     (localMF[idx_gphi+1]->nComp()!=nsolveMM_FACE)||
     (localMF[idx_gphi+BL_SPACEDIM-1]->nComp()!=nsolveMM_FACE)) 
  BoxLib::Error("invalid nComp");

 VOF_Recon_resize(1,SLOPE_RECON_MF);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[FSI_GHOST_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& ax = (*localMF[AREA_MF])[mfi];
  FArrayBox& ay = (*localMF[AREA_MF+1])[mfi];
  FArrayBox& az = (*localMF[AREA_MF+BL_SPACEDIM-1])[mfi];

  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
  FArrayBox& zface=(*localMF[FACE_VAR_MF+BL_SPACEDIM-1])[mfi];

  FArrayBox& vol = (*localMF[VOLUME_MF])[mfi];

  FArrayBox& ux = (*localMF[idx_gphi])[mfi];
  FArrayBox& uy = (*localMF[idx_gphi+1])[mfi];
  FArrayBox& uz = (*localMF[idx_gphi+BL_SPACEDIM-1])[mfi];

  FArrayBox& rhs = (*rhsmf)[mfi];

   // mask=1.0 at interior fine bc ghost cells
  FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   //1=not cov  0=cov
  FArrayBox& maskcoef = (*localMF[MASKCOEF_MF])[mfi];

  FArrayBox& poldfab = (*localMF[idx_phi])[mfi];
  FArrayBox& diffusionRHSfab = (*diffusionRHScell)[mfi];

  FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& solfab = (*localMF[FSI_GHOST_MF])[mfi];
  FArrayBox& cterm = (*localMF[ALPHACOEF_MF])[mfi];
  FArrayBox& maskdivresfab = (*localMF[MASK_DIV_RESIDUAL_MF])[mfi];
  FArrayBox& maskresfab = (*localMF[MASK_RESIDUAL_MF])[mfi];
  FArrayBox& maskSEMfab = (*localMF[MASKSEM_MF])[mfi];

  const Real* xlo = grid_loc[gridno].lo();

  Array<int> presbc;
  getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
  if (presbc.size()!=nsolveMM*BL_SPACEDIM*2)
   BoxLib::Error("presbc.size() invalid");

  Array<int> velbc=getBCArray(State_Type,gridno,0,
   num_materials_vel*BL_SPACEDIM);

// for heat equation:
// vol*(T-T^n)*(rho cv)/dt-vol*grad dot k grad T = -1/dt vol*div u+
//   diffusionRHS
//
// residual correction form: let T=dT+T0,
// vol*(dT+T0-T^n)*(rho cv)/dt-vol*div k grad dT = -1/dt vol*div u+
//   diffusionRHS+vol div k grad T0
//   u=u-dt k grad T0
// 
// vol*(p-p^n)/(rho c^2 dt*dt)-vol*grad dot grad p/rho = -1/dt vol*div u+
//   diffusionRHS
// RHS=p*alpha-(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)/dt+diffusionRHS
//
// we start with:
// alpha(p-p^n)- DGp=(-1/dt) D u+diffusionRHS
// suppose we have an initial guess pstar
// alpha(p-pstar-p^n)-DG(p-pstar)=(-1/dt) D u - alpha pstar + DG pstar+
//  diffusionRHS
// define dp=p-pstar
// alpha(dp)-DG(dp)=(-1/dt) Du -alpha pstar+alpha p^n + DG pstar + diffusionRHS
// alpha(dp)-DG(dp)=(-1/dt)D(u-dt G pstar)-alpha(pstar-p^n)+diffusionRHS=
//   (-1/dt)D(u-dt G pstar)+alpha(p^n-pstar)+diffusionRHS
// solution becomes: p=pstar+dp
// u gets u-dt G dp
// poldhold gets poldhold-dp
//
// poldhold holds p^n-pstar 
//

  int operation_flag=0;
  int energyflag=0; // not used when operation_flag==0
  int local_enable_spectral=enable_spectral;
  int use_VOF_weight=0;

   // in: NavierStokes::apply_div
  FORT_MAC_TO_CELL(
   &nsolveMM_FACE, 
   &num_materials_face,
   &ns_time_order, 
   &divu_outer_sweeps, 
   &num_divu_outer_sweeps, 
   &operation_flag,  // operation_flag=0
   &energyflag,
   temperature_primitive_variable.dataPtr(),
   &nmat,
   &nparts,
   &nparts_def,
   im_solid_map_ptr,
   prescribed_solid_scale.dataPtr(),
   added_weight.dataPtr(),
   &nten,
   &level, 
   &finest_level,
   &face_flag,
   &make_interface_incomp,
   &local_solvability_projection,
   &project_option,
   &local_enable_spectral,
   &fluxvel_index,
   &fluxden_index,
   &facevel_index,
   &facecut_index,
   &icefacecut_index,
   &curv_index,
   &conservative_tension_force,
   &pforce_index,
   &faceden_index,
   &icemask_index,
   &massface_index,
   &vofface_index,
   &ncphys,
   velbc.dataPtr(),
   presbc.dataPtr(), 
   &cur_time_slab, 
   &slab_step,
   &dt_slab,
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   ux.dataPtr(),ARLIM(ux.loVect()),ARLIM(ux.hiVect()),//xp
   uy.dataPtr(),ARLIM(uy.loVect()),ARLIM(uy.hiVect()),//yp
   uz.dataPtr(),ARLIM(uz.loVect()),ARLIM(uz.hiVect()),//zp
   ux.dataPtr(),ARLIM(ux.loVect()),ARLIM(ux.hiVect()),
   uy.dataPtr(),ARLIM(uy.loVect()),ARLIM(uy.hiVect()),
   uz.dataPtr(),ARLIM(uz.loVect()),ARLIM(uz.hiVect()),
   xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
   yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
   zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
   ax.dataPtr(),ARLIM(ax.loVect()),ARLIM(ax.hiVect()),
   ay.dataPtr(),ARLIM(ay.loVect()),ARLIM(ay.hiVect()),
   az.dataPtr(),ARLIM(az.loVect()),ARLIM(az.hiVect()),
   vol.dataPtr(),ARLIM(vol.loVect()),ARLIM(vol.hiVect()),
   rhs.dataPtr(),ARLIM(rhs.loVect()),ARLIM(rhs.hiVect()),
   rhs.dataPtr(),ARLIM(rhs.loVect()),ARLIM(rhs.hiVect()), // veldest
   rhs.dataPtr(),ARLIM(rhs.loVect()),ARLIM(rhs.hiVect()), // dendest
   maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
   ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   maskcoef.dataPtr(), // 1=not covered  0=covered
   ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
   maskSEMfab.dataPtr(), 
   ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
   reconfab.dataPtr(), // levelPC
   ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
   solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
   cterm.dataPtr(),ARLIM(cterm.loVect()),ARLIM(cterm.hiVect()),
   poldfab.dataPtr(),ARLIM(poldfab.loVect()),ARLIM(poldfab.hiVect()),
   poldfab.dataPtr(),ARLIM(poldfab.loVect()),ARLIM(poldfab.hiVect()),//denold
   poldfab.dataPtr(),ARLIM(poldfab.loVect()),ARLIM(poldfab.hiVect()),//ustar
   reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
   diffusionRHSfab.dataPtr(),
   ARLIM(diffusionRHSfab.loVect()),ARLIM(diffusionRHSfab.hiVect()),
   maskdivresfab.dataPtr(),
   ARLIM(maskdivresfab.loVect()),ARLIM(maskdivresfab.hiVect()),
   maskresfab.dataPtr(),
   ARLIM(maskresfab.loVect()),ARLIM(maskresfab.hiVect()),
   &SDC_outer_sweeps,
   &homflag,
   &use_VOF_weight,
   &nsolve,
   &SEM_advection_algorithm);
 } // mfi
} // omp
 ParallelDescriptor::Barrier();

 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {

  if (nsolve!=1)
   BoxLib::Error("nsolve invalid");

  for (int nc=0;nc<nsolveMM;nc++) {
   int nc_mask=nc;
 
   if ((nc_mask<0)||(nc_mask>=num_materials_face))
    BoxLib::Error("nc_mask invalid");
 
   MultiFab::Multiply(*rhsmf,*localMF[DOTMASK_MF],nc_mask,nc,1,0);
  } // nc 
    
 } else
  BoxLib::Error("num_materials_face invalid");

} // subroutine apply_div


// if temperature: -div(k grad T)-THERMAL_FORCE_MF
// if viscosity  : -div(2 mu D)-HOOP_FORCE_MARK_MF
// if mom force  : NEG_MOM_FORCE_MF  (project_option==4)
// called from: NavierStokes::do_the_advance
//              NavierStokes::veldiffuseALL
// called if project_option==0,2,3,4
void NavierStokes::update_SEM_forcesALL(int project_option,
 int idx_source,int update_spectral,int update_stable) {

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  BoxLib::Error("SDC_outer_sweeps invalid update_SEM_forcesALL");

 Real save_dt=dt_slab;
 dt_slab=1.0;

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==4)||   // NEG_MOM_FORCE
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
 } else if (project_option==2) {  // thermal diffusion
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid67");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nsolve=1;
 if (project_option==0) { // grad p, div(u p)
  nsolve=1;
 } else if (project_option==2) { // -div(k grad T)-THERMAL_FORCE_MF
  nsolve=1;
 } else if (project_option==3) { // -div(2 mu D)-HOOP_FORCE_MARK_MF
  nsolve=BL_SPACEDIM;
 } else if (project_option==4) { // NEG_MOM_FORCE_MF
  nsolve=BL_SPACEDIM;
 } else
  BoxLib::Error("project_option invalid68"); 

 if ((project_option==0)||   // grad p, div(u p)
     (project_option==2)||   // -div(k grad T)-THERMAL_FORCE_MF
     (project_option==3)) {  // -div(2 mu D)-HOOP_FORCE_MARK_MF

   // allocate and initialize to 0.0
  allocate_MAC_velocityALL(nsolve,UMAC_MF);
  allocate_MAC_velocityALL(nsolve,GP_DEST_FACE_MF);

  min_face_wt.resize(thread_class::nthreads);
  max_face_wt.resize(thread_class::nthreads);
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   min_face_wt[tid].resize(4);
   max_face_wt[tid].resize(4);
   for (int iwt=0;iwt<4;iwt++) {
    min_face_wt[tid][iwt]=1.0e+20;
    max_face_wt[tid][iwt]=-1.0e+20;
   }
  } // tid

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   
   ns_level.new_localMF(ONES_MF,num_materials_face,0,-1);
   ns_level.setVal_localMF(ONES_MF,1.0,0,num_materials_face,0);
   ns_level.ones_sum_global=0.0;

   ns_level.makeDotMask(nsolve,project_option);
   ns_level.allocate_FACE_WEIGHT(nsolve,project_option);
   ns_level.allocate_pressure_work_vars(nsolve,project_option);
  }

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.allocate_maccoef(project_option,nsolve);
  }

   // automatically initializes GP_DEST_CELL=0.0
  allocate_array(0,BL_SPACEDIM*num_materials_face,-1,GP_DEST_CELL_MF);

 } else if (project_option==4) { // NEG_MOM_FORCE_MF
  // do nothing
 } else
  BoxLib::Error("project_option invalid69");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.update_SEM_forces(project_option,
    idx_source,update_spectral,update_stable);
 }

 if ((project_option==0)||  // grad p, div(u p)
     (project_option==2)||  // -div(k grad T)-THERMAL_FORCE_MF
     (project_option==3)) { // -div(2 mu D)-HOOP_FORCE_MARK_MF

  delete_array(GP_DEST_CELL_MF);

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.deallocate_maccoef(project_option);
  }

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.remove_pressure_work_vars();
   ns_level.delete_localMF(FACE_WEIGHT_MF,BL_SPACEDIM);
   ns_level.delete_localMF(OFF_DIAG_CHECK_MF,1);
   ns_level.delete_localMF(DOTMASK_MF,1);
   ns_level.delete_localMF(ONES_MF,1);
  } // ilev=finest_level ... level

  remove_MAC_velocityALL(UMAC_MF);
  remove_MAC_velocityALL(GP_DEST_FACE_MF);

 } else if (project_option==4) { // NEG_MOM_FORCE_MF
  // do nothing
 } else
  BoxLib::Error("project_option invalid70");

 dt_slab=save_dt;

} // subroutine update_SEM_forcesALL

// called if project_option==0,2,3,4
void NavierStokes::update_SEM_forces(int project_option,
 int idx_source,int update_spectral,int update_stable) {

 if (dt_slab!=1.0)
  BoxLib::Error("dt_slab invalid5");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  BoxLib::Error("SDC_outer_sweeps invalid update_SEM_forces");

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==4)||   // -momforce
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
 } else if (project_option==2) { // thermal diffusion
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid71");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int local_idx_gp=GP_DEST_CELL_MF;
 int local_idx_gpmac=GP_DEST_FACE_MF;
 int local_idx_div=MACDIV_MF;
 
 int nsolve=1;
 if (project_option==0) { // grad p, div(u p)
  nsolve=1;
 } else if (project_option==2) { // -div(k grad T)-THERMAL_FORCE_MF
  nsolve=1;
 } else if (project_option==3) { // -div(2 mu D)-HOOP_FORCE_MARK_MF
  nsolve=BL_SPACEDIM;
 } else if (project_option==4) { // -mom force
  nsolve=BL_SPACEDIM;
  local_idx_gp=NEG_MOM_FORCE_MF;
  local_idx_gpmac=NEG_MOM_FORCE_MF;
  local_idx_div=NEG_MOM_FORCE_MF;
 } else
  BoxLib::Error("project_option invalid72"); 

 int nsolveMM=nsolve*num_materials_face;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  BoxLib::Error("num_materials_face invalid");

 if ((project_option==0)||   // grad p, div(u p)
     (project_option==2)||   // -div(k grad T)-THERMAL_FORCE_MF
     (project_option==3)) {  // -div(2 mu D)-HOOP_FORCE_MARK_MF

  if (localMF_grow[MACDIV_MF]>=0)
   BoxLib::Error("local div data not previously deleted");
  new_localMF(MACDIV_MF,nsolveMM,0,-1);

  int homflag=0;

  if ((project_option==2)||  // thermal diffusion
      (project_option==3)) { // viscosity

   // UMAC=-k grad T  (project_option=2)
   // UMAC=-2 mu D  (project_option=3)
   int energyflag=0;
   int simple_AMR_BC_flag=0;
   apply_pressure_grad(
    simple_AMR_BC_flag,
    homflag,energyflag,UMAC_MF,
    idx_source,
    project_option,nsolve);

   if (localMF[UMAC_MF]->nComp()!=nsolveMM_FACE)
    BoxLib::Error("localMF[UMAC_MF]->nComp() invalid");

   int ncomp_edge=-1;
   int scomp=0;
   avgDownEdge_localMF(UMAC_MF,scomp,ncomp_edge,0,BL_SPACEDIM,1,20);

   MultiFab* sourcemf=localMF[idx_source];
   MultiFab* rhs=localMF[MACDIV_MF];
   homflag=4;
   // rhs=div u
   apply_div(
    project_option,homflag,
    idx_source,
    rhs,
    sourcemf,
    UMAC_MF,
    nsolve);

  } else if (project_option==0) {

   // energyflag=2: 
   //   get grad p instead of \pm dt grad p/rho
   //   get div(up) instead of -dt div(up)/rho
   int energyflag=2;

   // GP_DEST_FACE=grad p instead of -dt grad p/rho
   int simple_AMR_BC_flag=0;
   apply_pressure_grad(
    simple_AMR_BC_flag,
    homflag,energyflag,GP_DEST_FACE_MF,
    idx_source,
    project_option,nsolve);

   for (int dir=0;dir<BL_SPACEDIM;dir++) {
    MultiFab* macvel=
     getStateMAC(0,dir,0,nsolveMM_FACE,cur_time_slab); 
    MultiFab::Copy(*localMF[UMAC_MF+dir],*macvel,0,0,nsolveMM_FACE,0);
    delete macvel;
   } 
   // grad p: GP_DEST_CELL_MF
   // div(up) : MACDIV_MF 
   apply_cell_pressure_gradient(
    project_option,
    energyflag,
    idx_source,
    UMAC_MF,
    GP_DEST_CELL_MF,
    MACDIV_MF);
  } else
   BoxLib::Error("project_option invalid73");

 } else if (project_option==4) { // -momforce
  // do nothing
 } else
  BoxLib::Error("project_option invalid74");

  // f=-div 2 mu D - HOOP_FORCE_MARK_MF  (project_option==3) or
  // f=NEG_MOM_FORCE_MF                  (project_option==4) or
  // f=-div k grad T - THERMAL_FORCE_MF  (project_option==2) or
  // f=grad p                            (project_option==0) or
  // f=div (up)
  // NavierStokes::update_SEM_delta_force (NavierStokes.cpp)
  // calls: FORT_UPDATESEMFORCE
  // does not look at enable_spectral
 update_SEM_delta_force(project_option,
   local_idx_gp,
   local_idx_gpmac,
   local_idx_div,
   update_spectral,update_stable,nsolve);

 if ((project_option==0)||  // grad p, div(u p)
     (project_option==2)||  // -div(k grad T)-THERMAL_FORCE_MF
     (project_option==3)) { // -div(2 mu D)-HOOP_FORCE_MARK_MF

  delete_localMF(MACDIV_MF,1);

 } else if (project_option==4) { // -momforce
  // do nothing
 } else
  BoxLib::Error("project_option invalid75");

} // subroutine update_SEM_forces

// if compressible: DIV_new=-dt(pnew-padv)/(rho c^2 dt^2)+MDOT_MF dt
// if incompressible: DIV_new=MDOT_MF dt
// called from NavierStokes::multiphase_project if project_option==11
// and called from NavierStokes::do_the_advance if advance_status==1.
void NavierStokes::ADVECT_DIV_ALL() {

 if (level!=0)
  BoxLib::Error("level invalid ADVECT_DIV_ALL");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel!=1");

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.ADVECT_DIV();
   // spectral_override==0 => always low order.
  int spectral_override=1;
  ns_level.avgDown(DIV_Type,0,num_materials_vel,spectral_override);
 }

} // subroutine ADVECT_DIV_ALL


// if compressible: DIV_new=-dt(pnew-padv)/(rho c^2 dt^2)+MDOT_MF dt
// if incompressible: DIV_new=MDOT_MF dt
void NavierStokes::ADVECT_DIV() {
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel!=1");

 debug_ngrow(FACE_VAR_MF,0,660);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,0,661);

 debug_ngrow(CELL_SOUND_MF,0,113);

 if (localMF[CELL_SOUND_MF]->nComp()!=2)
  BoxLib::Error("localMF[CELL_SOUND_MF]->nComp() invalid");

 if (num_state_base!=2)
  BoxLib::Error("num_state_base invalid");

 const Real* dx = geom.CellSize();

 MultiFab& DIV_new=get_new_data(DIV_Type,slab_step+1);
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int pcomp=num_materials_vel*BL_SPACEDIM;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(DIV_new,use_tiling); mfi.isValid(); ++mfi) {

  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  const Real* xlo = grid_loc[gridno].lo();

  // coeff_avg,padvect_avg,coeff1,padvect1, ... ,coeffn,padvectn 
  FArrayBox& csoundfab=(*localMF[CELL_SOUND_MF])[mfi];
  FArrayBox& mdotfab=(*localMF[MDOT_MF])[mfi];
  FArrayBox& snewfab=S_new[mfi];
  FArrayBox& divfab=DIV_new[mfi];

  FORT_UPDATE_DIV(
   xlo,dx,
   &dt_slab,
   csoundfab.dataPtr(),
   ARLIM(csoundfab.loVect()),ARLIM(csoundfab.hiVect()),
   mdotfab.dataPtr(),ARLIM(mdotfab.loVect()),ARLIM(mdotfab.hiVect()),
   snewfab.dataPtr(pcomp),
   ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
   divfab.dataPtr(),
   ARLIM(divfab.loVect()),ARLIM(divfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &nmat);

 } // mfi
} // omp
 ParallelDescriptor::Barrier();

} // subroutine ADVECT_DIV

void NavierStokes::getStateDIV_ALL(int idx,int ngrow) {

 if (level!=0)
  BoxLib::Error("level invalid getStateDIV_ALL");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  BoxLib::Error("SDC_outer_sweeps invalid");

 int finest_level=parent->finestLevel();

 int save_enable_spectral=enable_spectral;
 override_enable_spectral(projection_enable_spectral);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateDIV(idx,ngrow);
  int scomp=0;
  int ncomp=ns_level.localMF[idx]->nComp();
  ns_level.avgDown_localMF(idx,scomp,ncomp,0);
 }

 override_enable_spectral(save_enable_spectral);

} // subroutine getStateDIV_ALL 

void NavierStokes::getStateDIV(int idx,int ngrow) {

 bool use_tiling=ns_tiling;

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  BoxLib::Error("SDC_outer_sweeps invalid");

 if (localMF_grow[idx]>=0)
  BoxLib::Error("local div data not previously deleted");
 new_localMF(idx,num_materials_vel,ngrow,-1);

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int nsolve=1;
 int nsolveMM_FACE=nsolve*num_materials_vel;

 int fluxvel_index=0;
 int fluxden_index=BL_SPACEDIM;

 MultiFab* velmac[BL_SPACEDIM];
 for (int dir=0;dir<BL_SPACEDIM;dir++)
  velmac[dir]=getStateMAC(0,dir,0,nsolveMM_FACE,cur_time_slab);

 const Real* dx = geom.CellSize();

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,250);
 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,80);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,610);

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>=nmat))
  BoxLib::Error("nparts invalid");
 Array<int> im_solid_map_null;
 im_solid_map_null.resize(1);

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=nmat-1)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  BoxLib::Error("nparts invalid");

 resize_FSI_GHOST_MF(1);
 if (localMF[FSI_GHOST_MF]->nGrow()!=1)
  BoxLib::Error("localMF[FSI_GHOST_MF]->nGrow()!=1");
 if (localMF[FSI_GHOST_MF]->nComp()!=nparts_def*BL_SPACEDIM)
  BoxLib::Error("localMF[FSI_GHOST_MF]->nComp()!=nparts_def*BL_SPACEDIM");

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[idx],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   FArrayBox& ax = (*localMF[AREA_MF])[mfi];
   FArrayBox& ay = (*localMF[AREA_MF+1])[mfi];
   FArrayBox& az = (*localMF[AREA_MF+BL_SPACEDIM-1])[mfi];

   FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
   FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
   FArrayBox& zface=(*localMF[FACE_VAR_MF+BL_SPACEDIM-1])[mfi];

   FArrayBox& maskcoef = (*localMF[MASKCOEF_MF])[mfi];//1=not cov  0=cov
   FArrayBox& maskSEMfab = (*localMF[MASKSEM_MF])[mfi];

   FArrayBox& vol = (*localMF[VOLUME_MF])[mfi];
   FArrayBox& ux = (*velmac[0])[mfi];
   FArrayBox& uy = (*velmac[1])[mfi];
   FArrayBox& uz = (*velmac[BL_SPACEDIM-1])[mfi];
   FArrayBox& rhs = (*localMF[idx])[mfi];
   FArrayBox& solfab = (*localMF[FSI_GHOST_MF])[mfi];
   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
   const Real* xlo = grid_loc[gridno].lo();
   Array<int> velbc=getBCArray(State_Type,gridno,0,
    num_materials_vel*BL_SPACEDIM);

// RHS=(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)/vol_ij
//


   int operation_flag=1;
   int energyflag=0;
   int project_option=0;
   int homflag=0; // default
   int local_enable_spectral=enable_spectral;
   int num_materials_face=num_materials_vel;
   int use_VOF_weight=0;

    // in: NavierStokes::getStateDIV
   FORT_MAC_TO_CELL(
    &nsolveMM_FACE,
    &num_materials_face,
    &ns_time_order,
    &divu_outer_sweeps,
    &num_divu_outer_sweeps,
    &operation_flag, // operation_flag==1
    &energyflag,
    temperature_primitive_variable.dataPtr(),
    &nmat,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    prescribed_solid_scale.dataPtr(),
    added_weight.dataPtr(),
    &nten,
    &level, 
    &finest_level,
    &face_flag,
    &make_interface_incomp,
    &solvability_projection,
    &project_option,
    &local_enable_spectral,
    &fluxvel_index,
    &fluxden_index,
    &facevel_index,
    &facecut_index,
    &icefacecut_index,
    &curv_index,
    &conservative_tension_force,
    &pforce_index,
    &faceden_index,
    &icemask_index,
    &massface_index,
    &vofface_index,
    &ncphys,
    velbc.dataPtr(),
    velbc.dataPtr(), // presbc
    &cur_time_slab, 
    &slab_step,
    &dt_slab, 
    xlo,dx,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    ux.dataPtr(),ARLIM(ux.loVect()),ARLIM(ux.hiVect()),//xp
    uy.dataPtr(),ARLIM(uy.loVect()),ARLIM(uy.hiVect()),//yp
    uz.dataPtr(),ARLIM(uz.loVect()),ARLIM(uz.hiVect()),//zp
    ux.dataPtr(),ARLIM(ux.loVect()),ARLIM(ux.hiVect()),
    uy.dataPtr(),ARLIM(uy.loVect()),ARLIM(uy.hiVect()),
    uz.dataPtr(),ARLIM(uz.loVect()),ARLIM(uz.hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
    yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
    zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
    ax.dataPtr(),ARLIM(ax.loVect()),ARLIM(ax.hiVect()),
    ay.dataPtr(),ARLIM(ay.loVect()),ARLIM(ay.hiVect()),
    az.dataPtr(),ARLIM(az.loVect()),ARLIM(az.hiVect()),
    vol.dataPtr(),ARLIM(vol.loVect()),ARLIM(vol.hiVect()),
    rhs.dataPtr(),ARLIM(rhs.loVect()),ARLIM(rhs.hiVect()),
    rhs.dataPtr(),ARLIM(rhs.loVect()),ARLIM(rhs.hiVect()), // veldest
    rhs.dataPtr(),ARLIM(rhs.loVect()),ARLIM(rhs.hiVect()), // dendest
    maskcoef.dataPtr(),
    ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
    maskcoef.dataPtr(), // 1=not covered  0=covered
    ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
    maskSEMfab.dataPtr(), 
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    solfab.dataPtr(), //levelPC
    ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
    solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
    solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),//cterm
    solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),//pold
    solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),//denold
    solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),//ustar
    reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),// mdot
    solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),//maskdivres
    solfab.dataPtr(),ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),//maskres
    &SDC_outer_sweeps,
    &homflag,
    &use_VOF_weight,
    &nsolve,
    &SEM_advection_algorithm);
 } // mfi
} // omp
 ParallelDescriptor::Barrier();

 for (int dir=0;dir<BL_SPACEDIM;dir++)
  delete velmac[dir];

} // subroutine getStateDIV



// initializes righthand side and mac_phi for doing jacobi method.
void NavierStokes::mac_project_rhs(int project_option,
 int idx_mac_phi_crse,int idx_mac_rhs_crse,int nsolve) {

 int finest_level=parent->finestLevel();

 if (level>finest_level)
  BoxLib::Error("level invalid mac_project_rhs");

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists 2nd project
     (project_option==13)|| //FSI_material_exists 1st project
     (project_option==12)||
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid76");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;
 if (localMF[idx_mac_phi_crse]->nComp()!=nsolveMM)
  BoxLib::Error("localMF[idx_mac_phi_crse]->nComp()  invalid");
 if (localMF[idx_mac_rhs_crse]->nComp()!=nsolveMM)
  BoxLib::Error("localMF[idx_mac_rhs_crse]->nComp()  invalid");

 localMF[idx_mac_phi_crse]->setVal(0.0,0,nsolveMM,1);
 localMF[idx_mac_phi_crse]->setBndry(0.0);

   // residual correction form:
   // ap - div grad p = f
   // a(dp) - div grad dp= -a p* + div grad p* + f
   // (POLDHOLD=-p*)
   // rhs=POLDHOLD*alpha- vol div u/dt + diffusionRHS
 int homflag=0;
 apply_div(
   project_option,homflag,
   POLDHOLD_MF,
   localMF[idx_mac_rhs_crse], 
   localMF[DIFFUSIONRHS_MF],
   UMAC_MF,
   nsolve);

}  // mac_project_rhs

// gradpedge=-dt W grad p
void NavierStokes::mac_update(MultiFab* mac_phi_crse,int project_option,
  int nsolve) {

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists 2nd project
     (project_option==13)|| //FSI_material_exists 1st project
     (project_option==12)||
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   BoxLib::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  BoxLib::Error("project_option invalid77");

 if (num_materials_vel!=1)
  BoxLib::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  BoxLib::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;

 if (mac_phi_crse->nGrow()!=1)
  BoxLib::Error("mac_phi_crse has invalid ngrow");
 if (mac_phi_crse->nComp()!=nsolveMM)
  BoxLib::Error("mac_phi_crse has invalid ncomp");

 if (num_state_base!=2)
  BoxLib::Error("num_state_base invalid");

 Array<int> scomp;
 Array<int> ncomp;
 int state_index;
 int ncomp_check;

 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 
 if (ncomp_check!=nsolveMM)
  BoxLib::Error("ncomp_check invalid");

 MultiFab::Subtract(*localMF[POLDHOLD_MF],*mac_phi_crse,0,0,nsolveMM,0);

   // UMAC=UMAC+GRADPEDGE
 correct_velocity(project_option,
   UMAC_MF,UMAC_MF,GRADPEDGE_MF,nsolve);

 MultiFab& P_new=get_new_data(state_index,slab_step+1);

 int scomp_ilist=0;
 for (int ilist=0;ilist<scomp.size();ilist++) {
  MultiFab::Add(P_new,*mac_phi_crse,scomp_ilist,
   scomp[ilist],ncomp[ilist],1);
  scomp_ilist+=ncomp[ilist];
 }
 if (scomp_ilist!=nsolveMM)
  BoxLib::Error("scomp_ilist invalid");
 
} // subroutine mac_update

// adjust tolerance if too stringent.
void NavierStokes::adjust_tolerance(Real& error0,Real& error0_max,
  Real& parm_eps_abs,Real& parm_bot_atol,Real& parm_min_rel_error,
  int project_option) {

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| // FSI_material_exists 2nd project
     (project_option==13)|| // FSI_material_exists 1st project
     (project_option==12)|| // pressure extension project
     (project_option==3)||  // viscosity
     (project_option==2)||  // thermal diffusion
     ((project_option>=100)&&  // species
      (project_option<100+num_species_var))) {

  if (error0>error0_max)
   error0_max=error0;

  Real bot_rel_error=parm_min_rel_error*parm_bot_atol/parm_eps_abs;

  if (parm_eps_abs<parm_min_rel_error*error0_max) {
    parm_eps_abs=1.01*parm_min_rel_error*error0_max;
    parm_bot_atol=1.01*bot_rel_error*error0_max;
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "adjusting the tolerance project_option=" << 
       project_option << '\n';
     std::cout << "parm_min_rel_error, error0_max " << 
       parm_min_rel_error << ' ' <<
       error0_max << '\n';
     std::cout << "new tolerance: " << parm_eps_abs << '\n';
     std::cout << "new bottom tolerance: " << parm_bot_atol << '\n';
    } // ioproc?
  } // parm_eps_abs<parm_min_rel_error*error0_max

 } else
  BoxLib::Error("project_option invalid");

}  // adjust_tolerance


