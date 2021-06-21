//#include <winstd.H>
#include <AMReX_CoordSys.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BoxDomain.H>

#include <NavierStokes.H>
#include <LEVEL_F.H>
#include <GODUNOV_F.H>
#include <NAVIERSTOKES_F.H>
#include <DIFFUSION_F.H>
#include <MACOPERATOR_F.H>
#include <MG_F.H>

// residual correction form:
// alpha(p-p*)-div beta grad p = -div u*
// dp=p-p0
// alpha(dp+p0-p*)-div beta grad (dp+p0) = -div u*
// alpha dp - div beta grad dp = -div V + alpha(p*-p0)
// V=u*-beta grad p0
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
// called from: NavierStokes::update_SEM_forcesALL (MacProj.cpp)
//              NavierStokes::multiphase_project
namespace amrex{

void
NavierStokes::allocate_maccoefALL(int project_option,int nsolve,
		int create_hierarchy) {

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.allocate_maccoef(project_option,nsolve,create_hierarchy);
 }

  // FORT_EXTRAPFILL, pc_interp
 Vector<int> scompBC_map;
 scompBC_map.resize(1);
 scompBC_map[0]=0;

  // ngrow=1 scomp=0 ncomp=1 
 PCINTERP_fill_bordersALL(ONES_GROW_MF,1,0,1,State_Type,scompBC_map);

} // end subroutine allocate_maccoefALL

void
NavierStokes::allocate_maccoef(int project_option,int nsolve,
		int create_hierarchy) {

 int nmat=num_materials;

 if (project_option_projection(project_option)==1) {

  // do nothing
  
 } else if (project_option==12) {  // extension project

  // do nothing

 } else if (project_option==3) {  // viscosity

  // do nothing
  
 } else if ((project_option==2)||     // thermal diffusion
            ((project_option>=100)&&  // species
             (project_option<100+num_species_var))||
            (project_option==200)) {  // smooth temperature

  // do nothing
  
 } else
  amrex::Error("project_option invalid60");

 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int bfact=parent->Space_blockingFactor(level);

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;
 int ncomp_check;

  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (ncomp_check!=nsolve)
  amrex::Error("ncomp_check invalid");

 Real dt_diffuse=dt_slab;

 if (project_option==2) { // temperature diffusion
  if (lower_slab_time==0.0) {
   if (SDC_outer_sweeps!=0)
    amrex::Error("SDC_outer_sweeps invalid");
   if (slab_step==0) {
    if (initial_temperature_diffuse_duration>dt_diffuse)
     dt_diffuse=initial_temperature_diffuse_duration; 
   } else if ((slab_step>0)&&(slab_step<ns_time_order)) {
    // do nothing
   } else
    amrex::Error("slab_step invalid");
  } else if (lower_slab_time>0.0) {
   // do nothing
  } else
   amrex::Error("lower_slab_time invalid");
 } // project_option==2 (temperature diffusion)

 const Real* dx = geom.CellSize();
 const BoxArray& gridparm=grids;
 const Geometry& geomparm=geom;
 const DistributionMapping dmapparm=dmap;

 int local_use_mg_precond=0;

 if ((create_hierarchy==1)&&
     (use_mg_precond_in_mglib==1)&&
     (level==0)) {

  int bfact_mg_cutoff=1024;

  if ((bfact>=1)&&(bfact<=bfact_mg_cutoff)) {
   local_use_mg_precond=1;
  } else if (bfact>bfact_mg_cutoff) {
   local_use_mg_precond=0;
  } else
   amrex::Error("bfact invalid in allocate_maccoef");

 } else if ((create_hierarchy==0)||
   	    (use_mg_precond_in_mglib==0)||
	    ((level>=1)&&(level<=finest_level))) {
  local_use_mg_precond=0;
 } else
  amrex::Error("create_hierarchy, use_mg_precond, or level invalid");

 mac_op=new ABecLaplacian(
  ns_max_grid_size,
  gridparm,
  geomparm,
  dmapparm,
  bfact,
  level,
  project_option,
  nsolve,
  ns_tiling,
  local_use_mg_precond);

 mac_op->laplacian_solvability=0; // nonsingular matrix

 if (create_hierarchy==0) {
  // do nothing
 } else if (create_hierarchy==1) {

  if ((coarsest_ONES_level>=0)&&
      (coarsest_ONES_level<=finest_level)) {
   // do nothing
  } else
   amrex::Error("coarsest_ONES_level invalid");

  int all_singular_patches=1;

  if (level>coarsest_ONES_level) {
   mac_op->laplacian_solvability=0; // nonsingular matrix
   all_singular_patches=0;
  } else if ((level>=0)&&(level<=coarsest_ONES_level)) {
   if (color_ONES_count>0) {
    for (int icolor=0;icolor<color_ONES_count;icolor++) {
     if (singular_patch_flag[icolor]==0) {
      // do nothing
     } else if (singular_patch_flag[icolor]==1) {
      // do nothing
     
      //compressible cell or Dirichlet cell.
     } else if (singular_patch_flag[icolor]==2) {
      all_singular_patches=0;
     } else 
      amrex::Error("invalid singular_patch_flag[icolor]");
    } // icolor=0;icolor<color_ONES_count;icolor++
    if (all_singular_patches==0) {
     mac_op->laplacian_solvability=0; //nonsingular matrix
    } else if (all_singular_patches==1) {
     if (project_option_singular_possible(project_option)==1) {
      mac_op->laplacian_solvability=1; //singular matrix
     } else
      amrex::Error("all_singular_patches invalid");
    } else
     amrex::Error("all_singular_patches invalid");
   } else {
    std::cout << "level=" << level << '\n';
    std::cout << "coarsest_ONES_level=" << coarsest_ONES_level << '\n';
    std::cout << "color_ONES_count=" << color_ONES_count << '\n';
    amrex::Error("color_ONES_count invalid 1");
   }
  } else
   amrex::Error("level invalid");

  if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "allocate_maccoef level= " << level << 
     " all_singular_patches= " << all_singular_patches << '\n';
   }
  } //verbose>0

 } else
  amrex::Error("create_hierarchy invalid");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)");
 
 new_localMF(ALPHACOEF_MF,nsolve,0,-1);
 new_localMF(ALPHANOVOLUME_MF,nsolve,0,-1);
 localMF[ALPHANOVOLUME_MF]->setVal(0.0,0,nsolve,0);

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,200);
 debug_ngrow(FACE_VAR_MF,0,201);

  // ONES_MF=1 if diag>0  ONES_MF=0 if diag==0.
 debug_ngrow(ONES_MF,0,202);
 if (localMF[ONES_MF]->nComp()!=1)
  amrex::Error("localMF[ONES_MF]->nComp() invalid");
 debug_ngrow(ONES_GROW_MF,1,202);
 if (localMF[ONES_GROW_MF]->nComp()!=1)
  amrex::Error("localMF[ONES_GROW_MF]->nComp() invalid");

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,202);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,134);
 debug_ngrow(CELL_SOUND_MF,0,135);
 debug_ngrow(CELL_DEN_MF,1,136);
 debug_ngrow(CELL_VISC_MF,1,137);
 debug_ngrow(CELL_DEDT_MF,1,138);
 debug_ngrow(OFF_DIAG_CHECK_MF,0,139);

 if (localMF[OFF_DIAG_CHECK_MF]->nComp()!=nsolve)
  amrex::Error("localMF[OFF_DIAG_CHECK_MF]->nComp() invalid");
 if (localMF[CELL_DEN_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");
 if (localMF[CELL_VISC_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_VISC_MF]->nComp() invalid");
 if (localMF[CELL_DEDT_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 if (localMF[CELL_SOUND_MF]->nComp()!=2)
  amrex::Error("localMF[CELL_SOUND_MF]->nComp() invalid");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[ALPHANOVOLUME_MF]->boxArray().d_numPts());

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
  if (geom.IsRZ())
   rzflag=1;
  else if (geom.IsCartesian())
   rzflag=0;
  else if (geom.IsCYLINDRICAL())
   rzflag=3;
  else
   amrex::Error("CoordSys bust 51");

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // defined in MACOPERATOR_3D.F90
  FORT_SCALARCOEFF(
    &nsolve,
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
    &visc_coef,
    &angular_velocity,
    &dt_diffuse,
    &cur_time_slab,
    &project_option,
    &rzflag, 
    &solidheat_flag);

 }  // mfi
} // omp
 ns_reconcile_d_num(32);

 int GFM_flag=0;
 int adjust_temperature=0; 
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

  // alpha T - div beta grad T = f
 if (project_option==2) {

  if (is_phasechange==1) {
    // alphanovolume=(rho cv/(dt*fact))+(1/vol) sum_face Aface k_m/(theta dx)
   for (int im=0;im<2*nten;im++) {
    if (latent_heat[im]!=0.0) {
     if ((freezing_model[im]==0)|| // fully saturated
         (freezing_model[im]==5)||
         (freezing_model[im]==6))  // Palmore and Desjardins
      GFM_flag=1;
    } else if (latent_heat[im]==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] invalid");
   } // im=0..2 nten-1
  } else if (is_phasechange==0) {
   // do nothing
  } else
   amrex::Error("is_phasechange invalid");
 }

 if ((project_option>=100)&&
     (project_option<100+num_species_var)) {

  if (is_phasechange==1) {

   for (int im=0;im<2*nten;im++) {
    if (latent_heat[im]!=0.0) {
     if (freezing_model[im]==6) {  // Palmore and Desjardins
      int ispec=mass_fraction_id[im];
      if ((ispec>=1)&&(ispec<=num_species_var)) {
       if (ispec==project_option-100+1)
        GFM_flag=1;
      } else
       amrex::Error("ispec invalid");
     }
    } else if (latent_heat[im]==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] invalid");
   } // im=0.. 2 nten -1

  } else if (is_phasechange==0) {
   // do nothing
  } else
   amrex::Error("is_phasechange invalid");
 }

 if (GFM_flag==1) {
  stefan_solver_init(localMF[ALPHANOVOLUME_MF],
		  adjust_temperature,
		  project_option);
 } else if (GFM_flag==0) {
  // do nothing
 } else
  amrex::Error("GFM_flag invalid");

  // average down from level+1 to level.
 avgDown_localMF(ALPHANOVOLUME_MF,0,nsolve,0);
 Copy_localMF(ALPHACOEF_MF,ALPHANOVOLUME_MF,0,0,nsolve,0);

 for (int veldir=0;veldir<nsolve;veldir++) {

   // dest,source,scomp,dcomp,ncomp,ngrow
  Mult_localMF(ALPHACOEF_MF,VOLUME_MF,0,veldir,1,0);
 } // veldir=0,..,nsolve-1

 mac_op->aCoefficients(*localMF[ALPHACOEF_MF]);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[FACE_WEIGHT_MF+dir]->boxArray())
   amrex::Error("face_weight_stable boxarrays do not match");

  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[FACE_VAR_MF+dir]->boxArray())
   amrex::Error("face_var boxarrays do not match");

  new_localMF(BXCOEFNOAREA_MF+dir,nsolve,0,dir);
  new_localMF(BXCOEF_MF+dir,nsolve,0,dir);
  localMF[BXCOEFNOAREA_MF+dir]->setVal(1.0,0,nsolve,0);

  if (localMF[FACE_WEIGHT_MF+dir]->nComp()!=nsolve) 
   amrex::Error("localMF[FACE_WEIGHT_MF+dir]->nComp() invalid");
  if (localMF[BXCOEFNOAREA_MF+dir]->nComp()!=nsolve) 
   amrex::Error("localMF[BXCOEFNOAREA_MF+dir]->nComp() invalid");

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[ALPHACOEF_MF]->boxArray().d_numPts());

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

   int tid_current=ns_thread();
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // FACE_WEIGHT_MF initialized in BUILDFACEWT (LEVELSET_3D.F90)
    // BXCOEFNOAREA *= (facewtL + facewtR)/2
    // MULT_FACEWT is in MACOPERATOR_3D.F90
   FORT_MULT_FACEWT(
    &nsolve,
    bxfab.dataPtr(),ARLIM(bxfab.loVect()),ARLIM(bxfab.hiVect()),
    facefab.dataPtr(),ARLIM(facefab.loVect()),ARLIM(facefab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    xlo,dx,&dir);
  } // mfi
} // omp
  ns_reconcile_d_num(33);
 }  // dir=0..sdim-1

 new_localMF(DIAG_REGULARIZE_MF,nsolve,0,-1);
 new_localMF(MASK_DIV_RESIDUAL_MF,1,0,-1);
 new_localMF(MASK_RESIDUAL_MF,1,0,-1);
 localMF[DIAG_REGULARIZE_MF]->setVal(0.0,0,nsolve,0);
 localMF[MASK_DIV_RESIDUAL_MF]->setVal(0.0,0,1,0);
 localMF[MASK_RESIDUAL_MF]->setVal(0.0,0,1,0);

 debug_ngrow(DIFFUSIONRHS_MF,0,141);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(
   localMF[DIAG_REGULARIZE_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[DIAG_REGULARIZE_MF],use_tiling); 
		 mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

   // ONES_MF=1 if diag>0  ONES_MF=0 if diag==0.
  FArrayBox& ones_fab=(*localMF[ONES_MF])[mfi];

  // mask=tag if not covered by level+1 or outside the domain.
  FArrayBox& maskcov = (*localMF[MASKCOEF_MF])[mfi];
  FArrayBox& alpha = (*localMF[ALPHACOEF_MF])[mfi];

  FArrayBox& offdiagcheck=(*localMF[OFF_DIAG_CHECK_MF])[mfi];

  FArrayBox& maskdivresidfab = (*localMF[MASK_DIV_RESIDUAL_MF])[mfi];
  FArrayBox& maskresidfab = (*localMF[MASK_RESIDUAL_MF])[mfi];
  FArrayBox& mdotfab=(*localMF[DIFFUSIONRHS_MF])[mfi];

  FArrayBox& bxfab = (*localMF[BXCOEFNOAREA_MF])[mfi];
  FArrayBox& byfab = (*localMF[BXCOEFNOAREA_MF+1])[mfi];
  FArrayBox& bzfab = (*localMF[BXCOEFNOAREA_MF+AMREX_SPACEDIM-1])[mfi];
  FArrayBox& fwtxfab = (*localMF[FACE_WEIGHT_MF])[mfi];
  FArrayBox& fwtyfab = (*localMF[FACE_WEIGHT_MF+1])[mfi];
  FArrayBox& fwtzfab = (*localMF[FACE_WEIGHT_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];  
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];  
  FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];  

  Vector<int> bc;
  getBCArray_list(bc,state_index,gridno,scomp,ncomp);
  if (bc.size()!=nsolve*AMREX_SPACEDIM*2)
   amrex::Error("bc.size() invalid");

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
 
   // FORT_INIT_MASK_SING is in MACOPERATOR_3D.F90
   // initializes MASK_DIV_RESIDUAL and MASK_RESIDUAL
  FORT_INIT_MASK_SING(
    &level,
    &finest_level,
    &nsolve,
    &nmat,
    &project_option,
    &ncphys,
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
    yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
    zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
     // ONES_MF=1 if diag>0  ONES_MF=0 if diag==0.
    ones_fab.dataPtr(),ARLIM(ones_fab.loVect()),ARLIM(ones_fab.hiVect()),
    maskcov.dataPtr(),ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    alpha.dataPtr(),
    ARLIM(alpha.loVect()),ARLIM(alpha.hiVect()),
    offdiagcheck.dataPtr(),
    ARLIM(offdiagcheck.loVect()),ARLIM(offdiagcheck.hiVect()),
    maskdivresidfab.dataPtr(),
    ARLIM(maskdivresidfab.loVect()),ARLIM(maskdivresidfab.hiVect()),
    maskresidfab.dataPtr(),
    ARLIM(maskresidfab.loVect()),ARLIM(maskresidfab.hiVect()),
    mdotfab.dataPtr(),
    ARLIM(mdotfab.loVect()),ARLIM(mdotfab.hiVect()),
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
 ns_reconcile_d_num(35);

  // dest,soucre,scomp,dcomp,ncomp,ngrow
 Copy_localMF(ONES_GROW_MF,ONES_MF,0,0,1,0);

 Real min_interior_coeff=0.0;
 if (denconst_max>0.0) {
  if (mglib_min_coeff_factor>=1.0) {
   min_interior_coeff=1.0/(denconst_max*mglib_min_coeff_factor);
  } else
   amrex::Error("mglib_min_coeff_factor invalid");
 } else
  amrex::Error("denconst_max invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[BXCOEFNOAREA_MF+dir]->boxArray())
   amrex::Error("BXCOEFNOAREA boxarrays do not match");

  if (localMF[BXCOEFNOAREA_MF+dir]->nComp()!=nsolve) 
   amrex::Error("localMF[BXCOEFNOAREA_MF+dir]->nComp() invalid");

   // if project_option==0,
   //    project_option==1,
   //    project_option==11,  FSI_material_exists last project
   //    project_option==12,  extension project
  if (project_option_singular_possible(project_option)==1) {

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[ALPHACOEF_MF]->boxArray().d_numPts());

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

    int tid_current=ns_thread();
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // BXCOEFNOAREA = min_interior_coeff if not on the
    // edge of the domain and BXCOEFNOAREA previously = 0.0
    // REGULARIZE_BX is in MACOPERATOR_3D.F90
    FORT_REGULARIZE_BX(
     &nsolve,
     bxfab.dataPtr(),ARLIM(bxfab.loVect()),ARLIM(bxfab.hiVect()),
     &min_interior_coeff,
     domlo,domhi,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     xlo,dx,&dir);
   } // mfi
 } // omp
   ns_reconcile_d_num(33);

  } else if (project_option_singular_possible(project_option)==0) {
   // do nothing 
  } else
   amrex::Error("project_option_singular_possible invalid");

  int ncomp_edge=-1;
  int scomp_bx=0;
  int ncomp_mf=1;
  avgDownEdge_localMF(BXCOEFNOAREA_MF,scomp_bx,ncomp_edge,dir,ncomp_mf,0,17);
  Copy_localMF(BXCOEF_MF+dir,BXCOEFNOAREA_MF+dir,0,0,nsolve,0);
   // dest,source,scomp,dcomp,ncomp,ngrow
  for (int veldir=0;veldir<nsolve;veldir++)
   Mult_localMF(BXCOEF_MF+dir,AREA_MF+dir,0,veldir,1,0);

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[ALPHACOEF_MF]->boxArray().d_numPts());

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

    int tid_current=ns_thread();
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FORT_DIVIDEDX(
     &nsolve,
     bxfab.dataPtr(),ARLIM(bxfab.loVect()),ARLIM(bxfab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     xlo,dx,&dir);
  } // mfi
} // omp
  ns_reconcile_d_num(34);

  mac_op->bCoefficients(*localMF[BXCOEF_MF+dir],dir);
 }  // dir=0...sdim-1

  // generateCoefficients calls buildMatrix ranging from 
  // the finest mglib level (lev=0) down to the coarsest 
  // (lev=MG_numlevels_var-1)
  //
  // buildMatrix calls FORT_BUILDMAT
 mac_op->generateCoefficients();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(
   localMF[DIAG_REGULARIZE_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[DIAG_REGULARIZE_MF],use_tiling); 
		 mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  FArrayBox& alpha = (*localMF[ALPHACOEF_MF])[mfi];
  FArrayBox& diagfab = (*localMF[DIAG_REGULARIZE_MF])[mfi];
  FArrayBox& bxfab = (*localMF[BXCOEF_MF])[mfi];
  FArrayBox& byfab = (*localMF[BXCOEF_MF+1])[mfi];
  FArrayBox& bzfab = (*localMF[BXCOEF_MF+AMREX_SPACEDIM-1])[mfi];

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
 
   // FORT_NSGENERATE is in MACOPERATOR_3D.F90
   // initializes DIAG_REGULARIZE
  FORT_NSGENERATE(
    &level,
    &finest_level,
    &nsolve,
    &nmat,
    &project_option,
    &ncphys,
    alpha.dataPtr(),
    ARLIM(alpha.loVect()),ARLIM(alpha.hiVect()),
    diagfab.dataPtr(),
    ARLIM(diagfab.loVect()),ARLIM(diagfab.hiVect()),
    bxfab.dataPtr(),ARLIM(bxfab.loVect()),ARLIM(bxfab.hiVect()),
    byfab.dataPtr(),ARLIM(byfab.loVect()),ARLIM(byfab.hiVect()),
    bzfab.dataPtr(),ARLIM(bzfab.loVect()),ARLIM(bzfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact);
 }// mfi
} // omp
 ns_reconcile_d_num(35);

}  // subroutine allocate_maccoef

// called at end of pressure extrapolation.
void 
NavierStokes::restore_active_pressure(int save_mf) {

 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;
 int bfact=parent->Space_blockingFactor(level);

 debug_ngrow(save_mf,0,142);
 if (localMF[save_mf]->nComp()!=1)
  amrex::Error("localMF[save_mf]->nComp() invalid");
 debug_ngrow(OFF_DIAG_CHECK_MF,0,143);
 if (localMF[OFF_DIAG_CHECK_MF]->nComp()!=1)
  amrex::Error("localMF[OFF_DIAG_CHECK_MF]->nComp() invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int pcomp=AMREX_SPACEDIM;

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[OFF_DIAG_CHECK_MF]->boxArray().d_numPts());

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

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

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
 ns_reconcile_d_num(36);

} // end subroutine restore_active_pressure

void
NavierStokes::deallocate_maccoefALL(int project_option) {

 int finest_level=parent->finestLevel();

 if (level==0) {
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.deallocate_maccoef(project_option);
  }
 } else
  amrex::Error("level must be 0 in deallocate_maccoefALL");

} // end subroutine deallocate_maccoefALL

void
NavierStokes::deallocate_maccoef(int project_option) {

 if (project_option_is_valid(project_option)==1) {
  delete mac_op;
 } else
  amrex::Error("project_option invalid deallocate_maccoef");

 delete_localMF(DIAG_REGULARIZE_MF,1);
 delete_localMF(MASK_DIV_RESIDUAL_MF,1);
 delete_localMF(MASK_RESIDUAL_MF,1);
 delete_localMF(ALPHANOVOLUME_MF,1);
 delete_localMF(ALPHACOEF_MF,1);
 delete_localMF(BXCOEFNOAREA_MF,AMREX_SPACEDIM);
 delete_localMF(BXCOEF_MF,AMREX_SPACEDIM);

} // end subroutine deallocate_maccoef


// interpolates coarse data and ADDS it to the fine data.
// interpolates where cdiagsing>0
void
NavierStokes::AllinterpScalarMAC(
  MultiFab* coarsedata,MultiFab* finedata,
  BoxArray& cgridscen,BoxArray& fgridscen,
  DistributionMapping& fdmap,
  MultiFab* cdiagsing,MultiFab* fdiagsing,
  int nsolve,int project_option) {

 if (project_option_momeqn(project_option)==1) {

  // do nothing
  
 } else if (project_option_momeqn(project_option)==0) {

  // do nothing

 } else
  amrex::Error("project_option invalid61");

 int finest_level=parent->finestLevel();
 if ((level>finest_level)||(level<1))
  amrex::Error("level invalid AllinterpScalarMAC");

 if ((coarsedata->nComp()!=nsolve)||
     (finedata->nComp()!=nsolve)||
     (cdiagsing->nComp()!=1)||
     (fdiagsing->nComp()!=1))
  amrex::Error("invalid ncomp");

 BoxArray crse_cen_fine_BA(fgridscen.size());
 for (int i = 0; i < fgridscen.size(); ++i) {
  crse_cen_fine_BA.set(i,amrex::coarsen(fgridscen[i],2));
 }

 MultiFab& S_fine = *finedata;
 MultiFab& pcoarse = *coarsedata;
 const BoxArray& fgrids=S_fine.boxArray();

 BoxArray crse_S_fine_BA(fgrids.size());

 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }

 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,nsolve,0,
   MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());
 crse_S_fine.copy(pcoarse,0,0,nsolve);

 MultiFab crse_diagsing_fine(crse_S_fine_BA,crse_dmap,1,0,
   MFInfo().SetTag("crse_diagsing_fine"),FArrayBoxFactory());
 crse_diagsing_fine.copy(*cdiagsing,0,0,1);

 int bfact_f=parent->Space_blockingFactor(level);
 int bfact=parent->Space_blockingFactor(level-1);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgridscen[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  const int i = mfi.index();

  FArrayBox& crse_fab = crse_S_fine[mfi];
  const Box& cbox = crse_cen_fine_BA[i];

  FArrayBox& fine_fab = S_fine[mfi];

  FArrayBox& cdiagfab=crse_diagsing_fine[mfi];
  FArrayBox& fdiagfab=(*fdiagsing)[mfi];

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  for (int veldir=0;veldir<nsolve;veldir++) {
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
  } // veldir=0..nsolve-1
 }   // mfi
} //omp
 ns_reconcile_d_num(37);

}  // subroutine AllinterpScalarMAC

void
NavierStokes::interpScalarMAC(MultiFab* coarsedata,MultiFab* finedata,
 int nsolve,int project_option) {

  int finest_level=parent->finestLevel();
  if ((level>finest_level)||(level<1))
   amrex::Error("level invalid interpScalarMAC");

  NavierStokes& coarse_lev = getLevel(level-1);
  BoxArray& fgridscen=grids;
  DistributionMapping& fdmap=dmap;
  BoxArray& cgridscen=coarse_lev.grids;
  AllinterpScalarMAC(
    coarsedata,finedata,
    cgridscen,fgridscen,
    fdmap,
    coarse_lev.localMF[MASK_RESIDUAL_MF],
    localMF[MASK_RESIDUAL_MF],
    nsolve,project_option);

}  // subroutine interpScalarMAC

void
NavierStokes::Allaverage(
  MultiFab* coarsedata,MultiFab* finedata,
  BoxArray& cgridscen,BoxArray& fgridscen,
  DistributionMapping& fdmap,
  int iaverage,
  int scomp,int dcomp) {
 

 int finest_level=parent->finestLevel();
 if (level>=finest_level)
  amrex::Error("level invalid Allaverage");
 int bfact_fine=parent->Space_blockingFactor(level+1);
 int bfact_coarse=parent->Space_blockingFactor(level);

 BoxArray crse_cen_fine_BA(fgridscen.size());
 for (int i = 0; i < fgridscen.size(); ++i) {
  crse_cen_fine_BA.set(i,amrex::coarsen(fgridscen[i],2));
 }

 MultiFab& S_crse = *coarsedata;
 MultiFab& S_fine = *finedata;
 const BoxArray& fgrids=S_fine.boxArray();

 BoxArray crse_S_fine_BA(fgrids.size());

 if (fgrids.size()==fgridscen.size()) {
  // do nothing
 } else {
  amrex::Error("expecting: fgrids.size()==fgridscen.size()");
 }

 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }

 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,1,0,
   MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  const int i = mfi.index();

  FArrayBox& crse_fab = crse_S_fine[mfi];
  const Box& cbox = crse_cen_fine_BA[i];

  FArrayBox& fine_fab = S_fine[mfi];

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: MG_3D.F90; a low order average down.
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
 ns_reconcile_d_num(38);

// src,scomp,dcomp,ncomp
 S_crse.copy(crse_S_fine,0,dcomp,1);

} // Allaverage

void
NavierStokes::averageRhs(int idx_MF,int nsolve,int project_option) {

  if (project_option_momeqn(project_option)==1) {
   //do nothing
  } else if (project_option_momeqn(project_option)==0) {
   //do nothing
  } else
   amrex::Error("project_option_momeqn invalid62");

  if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
   amrex::Error("nsolve invalid25");
  if (localMF[idx_MF]->nComp()!=nsolve)
   amrex::Error("nsolve invalid26");

  int finest_level=parent->finestLevel();
  int iavg=0;
  if (level>=finest_level)
   amrex::Error("level invalid averageRhs");

  NavierStokes& fine_lev = getLevel(level+1);
  BoxArray& fgridscen=fine_lev.grids;
  DistributionMapping& fdmap=fine_lev.dmap;
  BoxArray& cgridscen=grids;

  for (int veldir=0;veldir<nsolve;veldir++) {
   Allaverage(
     localMF[idx_MF],fine_lev.localMF[idx_MF],
     cgridscen,fgridscen,
     fdmap,
     iavg,
     veldir,veldir); 
  }

}  // averageRhs

// x^{new}=x^{old}+D^{-1}(b-Ax^{old})
// on finest level,
// x^{new}=x^{old}+ILU^{-1}(b-Ax^{old})
void NavierStokes::JacobiALL(
 int idx_resid,int idx_rhs,int idx_xnew,
 int project_option,int nsolve) {

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("JacobiALL should only be called from level==0");

 int change_flag=0;
 project_right_hand_side(idx_xnew,project_option,change_flag);

 NavierStokes& ns_finest = getLevel(finest_level);
 ns_finest.applyBC_LEVEL(project_option,idx_xnew,nsolve);
 int apply_lev=0;

  // the smoother uses A_LOW: e.g.
  // z^{k+1}=z^{k}+D_LOW^{-1}(r-A_LOW z^{k})
 
 ns_finest.mac_op->Fsmooth(
  *ns_finest.localMF[idx_xnew],
  *ns_finest.localMF[idx_rhs],
  apply_lev,smooth_type);

 for (int ilev=finest_level-1;ilev>=0;ilev--) {
  NavierStokes& ns_level = getLevel(ilev);
  ns_level.DiagInverse(
   ns_level.localMF[idx_resid],
   ns_level.localMF[idx_xnew],nsolve,project_option);
 } 
 project_right_hand_side(idx_xnew,project_option,change_flag);

}  // subroutine JacobiALL

void NavierStokes::DiagInverse(
  MultiFab* resid,MultiFab* xnew,int nsolve,int project_option) {
 
 bool use_tiling=ns_tiling;

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,845);

 debug_ngrow(FACE_VAR_MF,0,205);

 int finest_level=parent->finestLevel();
 if (level > finest_level)
  amrex::Error("level too big");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn invalid63");

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid250");
 if (resid->nComp()!=nsolve)
  amrex::Error("resid->nComp() invalid");
 if (xnew->nComp()!=nsolve)
  amrex::Error("xnew->nComp() invalid");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(resid->boxArray().d_numPts());

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
  // mask=tag if not covered by level+1 or outside the domain.
  FArrayBox& mfab = (*localMF[MASKCOEF_MF])[mfi];
  FArrayBox& xoldfab = (*xnew)[mfi];

  FArrayBox& diagfab = (*localMF[DIAG_REGULARIZE_MF])[mfi];
  if (diagfab.nComp()!=nsolve)
   amrex::Error("diagfab.nComp() invalid");

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  for (int veldir=0;veldir<nsolve;veldir++) {
    // FORT_DIAGINV is in NAVIERSTOKES_3D.F90
   FORT_DIAGINV(
    diagfab.dataPtr(veldir),
    ARLIM(diagfab.loVect()),ARLIM(diagfab.hiVect()),
    residfab.dataPtr(veldir),ARLIM(residfab.loVect()),ARLIM(residfab.hiVect()),
    xnewfab.dataPtr(veldir),ARLIM(xnewfab.loVect()),ARLIM(xnewfab.hiVect()),
    xoldfab.dataPtr(veldir),ARLIM(xoldfab.loVect()),ARLIM(xoldfab.hiVect()),
    mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact);
  } // veldir
 } // mfi
} // omp
 ns_reconcile_d_num(39);

} // subroutine DiagInverse


// RHS=p^n*alpha-(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+... )/dt+diffusionRHS
// solving p*alpha-(bx_{i+1/2}(p_{i+1}-p_{i})+...)=RHS
// resid=rhs-Aphi
void NavierStokes::residALL(
  int project_option,
  int idx_rhs,int idx_resid,int idx_phi,int nsolve) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid residALL");

  // 1. (start) calls project_right_hand_side(idx_phi) 
  // 2. (end)   calls project_right_hand_side(idx_resid)
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
  // project_right_hand_side declared in NavierStokes.cpp.
 int change_flag=0;
 project_right_hand_side(idx_resid,project_option,change_flag);
} // subroutine residALL

void NavierStokes::applyALL(
  int project_option,
  int idx_phi,int idx_Aphi,int nsolve) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid applyALL");

   // in: applyALL
 int change_flag=0;
 project_right_hand_side(idx_phi,project_option,change_flag);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  ns_level.avgDown_localMF(idx_phi,0,nsolve,1);
  ns_level.applyBC_LEVEL(project_option,idx_phi,nsolve);

   // gradpedge=-dt W grad p
  int homflag=1;
  int energyflag=0;
  int simple_AMR_BC_flag=0;
  int simple_AMR_BC_flag_viscosity=0;
  ns_level.apply_pressure_grad(
   simple_AMR_BC_flag,
   simple_AMR_BC_flag_viscosity,
   homflag,
   energyflag,
   GRADPEDGE_MF,
   idx_phi,
   project_option,nsolve);

  if (ilev<finest_level) {
   int ncomp_edge=-1;
   int scomp=0;
   ns_level.avgDownEdge_localMF(GRADPEDGE_MF,scomp,ncomp_edge,0,
		   AMREX_SPACEDIM,1,18);
  }

  MultiFab* mdot_local=ns_level.localMF[DIFFUSIONRHS_MF];

  if (1==0) {
   Real nrm=mdot_local->norm2(0);
   std::cout << "ilev= " << ilev << "mdot norm2 = " << nrm << '\n';
  }

// Aphi=phi*alpha_dual+(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)/dt
  homflag=1;
  ns_level.apply_div(
   project_option,
   homflag,
   idx_phi,
   ns_level.localMF[idx_Aphi],
   mdot_local,
   GRADPEDGE_MF,
   nsolve);
 } // ilev=finest_level ... level

 project_right_hand_side(idx_Aphi,project_option,change_flag);

} // subroutine applyALL

// called from JacobiALL, applyALL, applyGradALL
// applyALL is called from residALL
// JacobiALL is called from jacobi_cycles
void NavierStokes::applyBC_LEVEL(int project_option,int idx_phi,int nsolve) {

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option invalid64");

 if (override_bc_to_homogeneous!=1) {
  std::cout << "override_bc_to_homogeneous= " <<
	  override_bc_to_homogeneous << '\n';
  amrex::Error("override_bc_to_homogeneous invalid1");
 }

 if (localMF[idx_phi]->nComp()!=nsolve)
  amrex::Error("invalid ncomp");
 if (localMF[idx_phi]->nGrow()!=1)
  amrex::Error("invalid ngrow");

 localMF[idx_phi]->setBndry(0.0);

 int bfact=parent->Space_blockingFactor(level);
 if ((bfact<1)||(bfact>64))
  amrex::Error("bfact out of range");

 int state_index;
 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp, 
  ncomp_check);
 if (ncomp_check!=nsolve)
  amrex::Error("nsolve invalid 898");

 Vector<int> scompBC_map;
 scompBC_map.resize(nsolve);

 int dcomp=0; 
 for (int ilist=0;ilist<scomp.size();ilist++) {
  for (int nc=0;nc<ncomp[ilist];nc++) {
   scompBC_map[dcomp]=scomp[ilist]+nc;
   dcomp++;
  }
 }
 if (dcomp!=nsolve)
  amrex::Error("dcomp invalid"); 

 MultiFab* cmf=nullptr;

 if (level>0) {
  NavierStokes& ns_coarse=getLevel(level-1);
  ns_coarse.localMF[idx_phi]->setBndry(0.0);
  ns_coarse.localMF[idx_phi]->FillBoundary(ns_coarse.geom.periodicity());
  cmf=ns_coarse.localMF[idx_phi];
 } else if (level==0) {
  cmf=localMF[idx_phi];
 } else
  amrex::Error("level invalid applyBC_LEVEL");

 InterpBorders(
   *cmf,
   *localMF[idx_phi],
   cur_time_slab,
   state_index,
   0,  // scomp=0
   scompBC_map,
   nsolve,
   debug_fillpatch);

}  // subroutine applyBC_LEVEL

//homflag==1 => going down the V-cycle
//homflag==0 => going up the V-cycle
void NavierStokes::applyBC_MGLEVEL(int idx_phi,
 MultiFab* pbdry,int homflag,int nsolve,int project_option) {

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn invalid65");

 if (override_bc_to_homogeneous!=1) {
  std::cout << "override_bc_to_homogeneous= " <<
	  override_bc_to_homogeneous << '\n';
  amrex::Error("override_bc_to_homogeneous invalid2");
 }

 int bfact=parent->Space_blockingFactor(level);
 if ((bfact<1)||(bfact>64))
  amrex::Error("bfact out of range");

 if ((homflag!=1)&&(homflag!=0))
  amrex::Error("homflag invalid");

 if (localMF[idx_phi]->nComp()!=nsolve)
  amrex::Error("invalid ncomp");
 if (localMF[idx_phi]->nGrow()!=1)
  amrex::Error("invalid ngrow");

  // down the V-cycle
 if (homflag==1) 
  setVal_localMF(idx_phi,0.0,0,nsolve,1);

 localMF[idx_phi]->setBndry(0.0);

 int state_index;
 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 if (ncomp_check!=nsolve)
  amrex::Error("nsolve invalid 976");

  // up the V-cycle
 if ((homflag==0)&&(level>0)) {
  NavierStokes& ns_coarse=getLevel(level-1);
  ns_coarse.localMF[idx_phi]->setBndry(0.0);
  ns_coarse.localMF[idx_phi]->FillBoundary(ns_coarse.geom.periodicity());

  Vector<int> scompBC_map;
  scompBC_map.resize(nsolve);
  int dcomp=0;
  for (int ilist=0;ilist<scomp.size();ilist++) {
   for (int nc=0;nc<ncomp[ilist];nc++) {
    scompBC_map[dcomp]=scomp[ilist]+nc;
    dcomp++;
   }
  }
  if (dcomp!=nsolve)
   amrex::Error("dcomp invalid");

  InterpBorders(*ns_coarse.localMF[idx_phi],
    *localMF[idx_phi],
    cur_time_slab,
    state_index,
    0,  // scomp=0
    scompBC_map,
    nsolve,
    debug_fillpatch);
  interpScalarMAC(ns_coarse.localMF[idx_phi],
                  localMF[idx_phi],nsolve,project_option);
 }
 
 MultiFab::Copy(*pbdry,*localMF[idx_phi],0,0,nsolve,1);

} // applyBC_MGLEVEL


   // gradpedge=-dt W grad p
void NavierStokes::applyGradALL(
  int project_option,int idx_phi,int nsolve) {

 int homflag=1;
 int energyflag=0;
 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid applyGradALL");
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.applyBC_LEVEL(project_option,idx_phi,nsolve);

   // gradpedge=-dt W grad p
  int simple_AMR_BC_flag=0;
  int simple_AMR_BC_flag_viscosity=0;

  ns_level.apply_pressure_grad(
   simple_AMR_BC_flag,
   simple_AMR_BC_flag_viscosity,
   homflag,
   energyflag,
   GRADPEDGE_MF,
   idx_phi,
   project_option,nsolve);
  if (ilev<finest_level) {
   int ncomp_edge=-1;
   int scomp_edge=0;
   int start_dir=0;
   int spectral_override=1; // order determined from enable_spectral
   ns_level.avgDownEdge_localMF(
    GRADPEDGE_MF,
    scomp_edge,ncomp_edge,
    start_dir,AMREX_SPACEDIM,spectral_override,19);
  }
 } // ilev
} // subroutine applyGradALL

// NOTE: POLDHOLD=p^advect-p^init
// homflag=0 =>
// called from mac_project_rhs =>
// RHS=POLDHOLD*alpha-
//     (a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2})/dt+diffusionRHS
//
// homflag=1 =>
// called from applyALL =>
// RHS=phi*alpha+(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)/dt
//
// homflag=2 =>
// called from relaxLEVEL which is called from mg_cycleALL =>
// RHS=-phi*alpha-
// (a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2})/dt+diffusionRHS
// aka: residmf=-idx_phi*alpha-
//              (a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2})/dt+
//              rhsmf 
//
// homflag=3 =>
// called from relaxLEVEL which is called from mg_cycleALL
// (phi=phi_dual=0, u=0 on finest level) =>
// rhsmf=-(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2})/dt+idx_rhs
//
// homflag=4 =>
// called from update_SEM_forces =>
// rhsmf=div u
// (phi=phi_dual)
//
void NavierStokes::apply_div(
  int project_option,int homflag,
  int idx_phi,
  MultiFab* rhsmf, 
  MultiFab* diffusionRHScell,
  int idx_gphi,
  int nsolve) {

 int nmat=num_materials;

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn invalid66");

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;
 int ncomp_check;

  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (ncomp_check!=nsolve)
  amrex::Error("ncomp_check invalid");

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>nmat))
  amrex::Error("nparts invalid");
 Vector<int> im_solid_map_null;
 im_solid_map_null.resize(1);

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=nmat)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  amrex::Error("nparts invalid");

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

 const Real* dx = geom.CellSize();

 if (homflag==0) {
  if (idx_phi==POLDHOLD_MF) {
   // do nothing
  } else
   amrex::Error("expecting POLDHOLD");
 } else if (homflag==1) {
  // do nothing
 } else if (homflag==2) {
  // do nothing
 } else if (homflag==3) {
  // do nothing
 } else if (homflag==4) {
  // do nothing
 } else
  amrex::Error("homflag invalid");

 int fluxvel_index=0;
 int fluxden_index=AMREX_SPACEDIM;

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,250);
 debug_ngrow(FACE_VAR_MF,0,251);
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }
 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 debug_ngrow(MASKCOEF_MF,1,253); // maskcoef=1 if not covered by finer lev.
 debug_ngrow(MASK_NBR_MF,1,253); // mask_nbr=1 at fine-fine bc.

 debug_ngrow(DIAG_REGULARIZE_MF,0,253); 
 if (localMF[DIAG_REGULARIZE_MF]->nComp()!=nsolve)
  amrex::Error("localMF[DIAG_REGULARIZE_MF]->nComp()!=nsolve");

 debug_ngrow(ALPHACOEF_MF,0,253); 
 if (localMF[ALPHACOEF_MF]->nComp()!=nsolve)
  amrex::Error("localMF[ALPHACOEF_MF]->nComp()!=nsolve");

 if (diffusionRHScell->nGrow()<0)
  amrex::Error("diffusionRHScell invalid");
 if ((localMF[idx_phi]->nComp()!=nsolve)||
     (rhsmf->nComp()!=nsolve)||
     (diffusionRHScell->nComp()!=nsolve)||
     (localMF[idx_gphi]->nComp()!=nsolve)||
     (localMF[idx_gphi+1]->nComp()!=nsolve)||
     (localMF[idx_gphi+AMREX_SPACEDIM-1]->nComp()!=nsolve)) 
  amrex::Error("invalid nComp");

 VOF_Recon_resize(1,SLOPE_RECON_MF);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[SLOPE_RECON_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[SLOPE_RECON_MF],use_tiling); mfi.isValid(); ++mfi) {
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
  FArrayBox& az = (*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
  FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& vol = (*localMF[VOLUME_MF])[mfi];

  FArrayBox& ux = (*localMF[idx_gphi])[mfi];
  FArrayBox& uy = (*localMF[idx_gphi+1])[mfi];
  FArrayBox& uz = (*localMF[idx_gphi+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& rhs = (*rhsmf)[mfi];

   // mask=1.0 at interior fine bc ghost cells
  FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   //1=not cov  0=cov
  FArrayBox& maskcoef = (*localMF[MASKCOEF_MF])[mfi];

   // DIAG_REGULARIZE_MF used for sanity check only.
  FArrayBox& diagfab = (*localMF[DIAG_REGULARIZE_MF])[mfi];

  FArrayBox& poldfab = (*localMF[idx_phi])[mfi];

  FArrayBox& diffusionRHSfab = (*diffusionRHScell)[mfi];

  FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];

  FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
  FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
  FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& cterm = (*localMF[ALPHACOEF_MF])[mfi];

  FArrayBox& maskdivresfab = (*localMF[MASK_DIV_RESIDUAL_MF])[mfi];
  FArrayBox& maskresfab = (*localMF[MASK_RESIDUAL_MF])[mfi];
  FArrayBox& maskSEMfab = (*localMF[MASKSEM_MF])[mfi];

  const Real* xlo = grid_loc[gridno].lo();

  Vector<int> presbc;
  getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
  if (presbc.size()!=nsolve*AMREX_SPACEDIM*2)
   amrex::Error("presbc.size() invalid");

  Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

// for heat equation:
// vol*(T-T^n)*(rho cv)/dt-vol*grad dot k grad T = -1/dt vol*div u+
//   diffusionRHS
//
// residual correction form: let T=dT+T0,
// vol*(dT+T0-T^n)*(rho cv)/dt-vol*div k grad dT = -1/dt vol*div u+
//   diffusionRHS+vol div k grad T0
//   u=u-dt k grad T0
// 
// a=vol/(rho c^2 dt^2)
// a p-div beta grad p=- div u/dt+a p^adv+diffusionRHS
// p=dp+pguess
// pguess=-POLDHOLD+p^adv
// p^last=-POLDHOLD+p^adv
// V=-dt beta grad pguess
// a (dp+pguess) - div beta grad(dp+pguess)= - div u/dt+a p^adv+diffusionRHS
// a dp - div beta grad dp = 
//   -div u/dt + a p^adv -a pguess + diffusionRHS + div beta grad pguess
// =-div (u+V)/dt + a(p^adv-pguess)+diffusionRHS
// =-div (u+V)/dt + a POLDHOLD +diffusionRHS
// RHS=POLDHOLD*a-
//     (a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)/dt+diffusionRHS
//

  int operation_flag=0;
  int energyflag=0; // not used when operation_flag==0
  int local_enable_spectral=enable_spectral;
  int use_VOF_weight=0;

  int ncomp_denold=nsolve;
  int ncomp_veldest=cterm.nComp();
  int ncomp_dendest=poldfab.nComp();

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: NavierStokes::apply_div
   // fort_mac_to_cell declared in: LEVELSET_3D.F90
  fort_mac_to_cell(
   &ns_time_order, 
   &divu_outer_sweeps, 
   &num_divu_outer_sweeps, 
   &operation_flag,  // operation_flag=0
   &energyflag,
   temperature_primitive_variable.dataPtr(),
   constant_density_all_time.dataPtr(),
   &nmat,
   &nparts,
   &nparts_def,
   im_solid_map_ptr,
   added_weight.dataPtr(),
   &nten,
   &level, 
   &finest_level,
   &face_flag,
   &project_option,
   &local_enable_spectral,
   &fluxvel_index,
   &fluxden_index,
   &facevel_index,
   &facecut_index,
   &icefacecut_index,
   &curv_index,
   &conservative_tension_force,
   &conservative_div_uu,
   filter_velocity.dataPtr(),
   &ignore_div_up,
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
   cterm.dataPtr(),
   ARLIM(cterm.loVect()),ARLIM(cterm.hiVect()), // veldest
   poldfab.dataPtr(),
   ARLIM(poldfab.loVect()),ARLIM(poldfab.hiVect()), // dendest
   maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
   ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   maskcoef.dataPtr(), // 1=not covered  0=covered
   ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
   maskSEMfab.dataPtr(), 
   ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
   reconfab.dataPtr(), // levelPC
   ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
   solxfab.dataPtr(),
   ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
   solyfab.dataPtr(),
   ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
   solzfab.dataPtr(),
   ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
   cterm.dataPtr(),ARLIM(cterm.loVect()),ARLIM(cterm.hiVect()),
   poldfab.dataPtr(),ARLIM(poldfab.loVect()),ARLIM(poldfab.hiVect()),
   diagfab.dataPtr(),
   ARLIM(diagfab.loVect()),ARLIM(diagfab.hiVect()),//denold
   poldfab.dataPtr(),ARLIM(poldfab.loVect()),ARLIM(poldfab.hiVect()),//ustar
   reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
   diffusionRHSfab.dataPtr(), //mdotcell
   ARLIM(diffusionRHSfab.loVect()),ARLIM(diffusionRHSfab.hiVect()),
   maskdivresfab.dataPtr(),
   ARLIM(maskdivresfab.loVect()),ARLIM(maskdivresfab.hiVect()),
   maskresfab.dataPtr(),
   ARLIM(maskresfab.loVect()),ARLIM(maskresfab.hiVect()),
   &SDC_outer_sweeps,
   &homflag,
   &use_VOF_weight,
   &nsolve,
   &ncomp_denold,
   &ncomp_veldest,
   &ncomp_dendest,
   &SEM_advection_algorithm);
 } // mfi
} // omp
 ns_reconcile_d_num(40);

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
  amrex::Error("SDC_outer_sweeps invalid update_SEM_forcesALL");

 Real save_dt=dt_slab;
 dt_slab=1.0;

 int finest_level=parent->finestLevel();

 if ((project_option==0)||
     (project_option==4)||   // NEG_MOM_FORCE
     (project_option==3)) {  // viscosity
  //do nothing
 } else if (project_option==2) {  // thermal diffusion
  //do nothing
 } else
  amrex::Error("project_option invalid67");

 int nsolve=1;
 if (project_option==0) { // grad p, div(u p)
  nsolve=1;
 } else if (project_option==2) { // -div(k grad T)-THERMAL_FORCE_MF
  nsolve=1;
 } else if (project_option==3) { // -div(2 mu D)-HOOP_FORCE_MARK_MF
  nsolve=AMREX_SPACEDIM;
 } else if (project_option==4) { // NEG_MOM_FORCE_MF
  nsolve=AMREX_SPACEDIM;
 } else
  amrex::Error("project_option invalid68"); 

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
   
   ns_level.new_localMF(DIFFUSIONRHS_MF,nsolve,0,-1);
   ns_level.setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolve,0);

    // in: NavierStokes::update_SEM_forcesALL
    // ONES_MF=1 if diag>0  ONES_MF=0 if diag==0.
   ns_level.new_localMF(ONES_MF,1,0,-1);
   ns_level.new_localMF(ONES_GROW_MF,1,1,-1);
   ns_level.setVal_localMF(ONES_MF,1.0,0,1,0);
   ns_level.setVal_localMF(ONES_GROW_MF,1.0,0,1,1);

   ns_level.allocate_FACE_WEIGHT(nsolve,project_option);
   ns_level.allocate_pressure_work_vars(nsolve,project_option);
  } // ilev=finest_level ... level

  int create_hierarchy=0;
  allocate_maccoefALL(project_option,nsolve,create_hierarchy);

   // automatically initializes GP_DEST_CELL=0.0
  allocate_array(0,AMREX_SPACEDIM,-1,GP_DEST_CELL_MF);

 } else if (project_option==4) { // NEG_MOM_FORCE_MF
  // do nothing
 } else
  amrex::Error("project_option invalid69");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.update_SEM_forces(project_option,
    idx_source,update_spectral,update_stable);
 }

 if ((project_option==0)||  // grad p, div(u p)
     (project_option==2)||  // -div(k grad T)-THERMAL_FORCE_MF
     (project_option==3)) { // -div(2 mu D)-HOOP_FORCE_MARK_MF

  delete_array(GP_DEST_CELL_MF);

  deallocate_maccoefALL(project_option);

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.remove_pressure_work_vars();
   ns_level.delete_localMF(FACE_WEIGHT_MF,AMREX_SPACEDIM);
   ns_level.delete_localMF(OFF_DIAG_CHECK_MF,1);
   ns_level.delete_localMF(ONES_MF,1);
   ns_level.delete_localMF(ONES_GROW_MF,1);
   ns_level.delete_localMF(TYPE_ONES_MF,1);
   ns_level.delete_localMF(COLOR_ONES_MF,1);
   ns_level.delete_localMF(DIFFUSIONRHS_MF,1);
  } // ilev=finest_level ... level

  remove_MAC_velocityALL(UMAC_MF);
  remove_MAC_velocityALL(GP_DEST_FACE_MF);

 } else if (project_option==4) { // NEG_MOM_FORCE_MF
  // do nothing
 } else
  amrex::Error("project_option invalid70");

 dt_slab=save_dt;

} // subroutine update_SEM_forcesALL

// called if project_option==0,2,3,4
void NavierStokes::update_SEM_forces(int project_option,
 int idx_source,int update_spectral,int update_stable) {

 if (dt_slab==1.0) {
  // do nothing
 } else
  amrex::Error("dt_slab invalid in update_SEM_forces (5)");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid update_SEM_forces");

 if ((project_option==0)||
     (project_option==4)||   // -momforce
     (project_option==3)) {  // viscosity
  //do nothing
 } else if (project_option==2) { // thermal diffusion
  //do nothing
 } else
  amrex::Error("project_option invalid71");

 int local_idx_gp=GP_DEST_CELL_MF;
 int local_idx_gpmac=GP_DEST_FACE_MF;
 int local_idx_div=MACDIV_MF;
 
 int nsolve=1;
 if (project_option==0) { // grad p, div(u p)
  nsolve=1;
 } else if (project_option==2) { // -div(k grad T)-THERMAL_FORCE_MF
  nsolve=1;
 } else if (project_option==3) { // -div(2 mu D)-HOOP_FORCE_MARK_MF
  nsolve=AMREX_SPACEDIM;
 } else if (project_option==4) { // -mom force
  nsolve=AMREX_SPACEDIM;
  local_idx_gp=NEG_MOM_FORCE_MF;
  local_idx_gpmac=NEG_MOM_FORCE_MF;
  local_idx_div=NEG_MOM_FORCE_MF;
 } else
  amrex::Error("project_option invalid72"); 

 if ((project_option==0)||   // grad p, div(u p)
     (project_option==2)||   // -div(k grad T)-THERMAL_FORCE_MF
     (project_option==3)) {  // -div(2 mu D)-HOOP_FORCE_MARK_MF

  if (localMF_grow[MACDIV_MF]==-1) {
   new_localMF(MACDIV_MF,nsolve,0,-1);
  } else
   amrex::Error("localMF_grow[MACDIV_MF] invalid");

  int homflag=0;

  if ((project_option==2)||  // thermal diffusion
      (project_option==3)) { // viscosity

   // note: dt_slab=1 in update_SEM_forcesALL
   // UMAC=-dt_slab k grad T  (project_option=2)
   // UMAC=-dt_slab 2 mu D  (project_option=3)
   int energyflag=0;
   int simple_AMR_BC_flag=0;
   int simple_AMR_BC_flag_viscosity=0;
   apply_pressure_grad(
    simple_AMR_BC_flag,
    simple_AMR_BC_flag_viscosity,
    homflag,energyflag,UMAC_MF,
    idx_source,
    project_option,nsolve);

   if (localMF[UMAC_MF]->nComp()!=nsolve)
    amrex::Error("localMF[UMAC_MF]->nComp() invalid");

   int ncomp_edge=-1;
   int scomp=0;
   avgDownEdge_localMF(UMAC_MF,scomp,ncomp_edge,0,AMREX_SPACEDIM,1,20);

   MultiFab* sourcemf=localMF[idx_source];
   MultiFab* rhs=localMF[MACDIV_MF];
   homflag=4;
   // rhs=div u
   apply_div(
    project_option,
    homflag,
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
   int simple_AMR_BC_flag_viscosity=0;
   apply_pressure_grad(
    simple_AMR_BC_flag,
    simple_AMR_BC_flag_viscosity,
    homflag,energyflag,GP_DEST_FACE_MF,
    idx_source,
    project_option,nsolve);

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    MultiFab* macvel=
     getStateMAC(Umac_Type,0,dir,0,nsolve,cur_time_slab); 
    MultiFab::Copy(*localMF[UMAC_MF+dir],*macvel,0,0,nsolve,0);
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
   amrex::Error("project_option invalid73");

 } else if (project_option==4) { // -momforce
  // do nothing
 } else
  amrex::Error("project_option invalid74");

  // f=-div 2 mu D - HOOP_FORCE_MARK_MF  (project_option==3) or
  // f=NEG_MOM_FORCE_MF                  (project_option==4) or
  // f=-div k grad T - THERMAL_FORCE_MF  (project_option==2) or
  // f=grad p                            (project_option==0) or
  // f=div (up)
  // NavierStokes::update_SEM_delta_force (NavierStokes.cpp)
  // calls: FORT_UPDATESEMFORCE
  // does not look at enable_spectral
 if ((update_spectral+update_stable>=1)&&
     (update_spectral+update_stable<=2)) {
  update_SEM_delta_force(project_option,
   local_idx_gp,
   local_idx_gpmac,
   local_idx_div,
   update_spectral,update_stable,nsolve);
 } else
  amrex::Error("update_spectral+update_stable invalid");

 if ((project_option==0)||  // grad p, div(u p)
     (project_option==2)||  // -div(k grad T)-THERMAL_FORCE_MF
     (project_option==3)) { // -div(2 mu D)-HOOP_FORCE_MARK_MF

  delete_localMF(MACDIV_MF,1);

 } else if (project_option==4) { // -momforce
  // do nothing
 } else
  amrex::Error("project_option invalid75");

} // subroutine update_SEM_forces

// if compressible: DIV_new=-dt(pnew-padv)/(rho c^2 dt^2)+MDOT_MF dt/vol
// if incompressible: DIV_new=MDOT_MF dt/vol
// called from NavierStokes::multiphase_project if project_option==11
// and called from NavierStokes::do_the_advance if advance_status==1.
void NavierStokes::ADVECT_DIV_ALL() {

 if (level!=0)
  amrex::Error("level invalid ADVECT_DIV_ALL");

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.ADVECT_DIV();
   // spectral_override==0 => always low order.
  int spectral_override=1;
  ns_level.avgDown(DIV_Type,0,1,spectral_override);
 }

} // subroutine ADVECT_DIV_ALL


// if compressible: DIV_new=-dt(pnew-padv)/(rho c^2 dt^2)+MDOT_MF dt/vol
// if incompressible: DIV_new=MDOT_MF dt/vol
void NavierStokes::ADVECT_DIV() {
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,700);

 debug_ngrow(FACE_VAR_MF,0,660);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,0,661);

 debug_ngrow(CELL_SOUND_MF,0,144);

 if (localMF[CELL_SOUND_MF]->nComp()!=2)
  amrex::Error("localMF[CELL_SOUND_MF]->nComp() invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 const Real* dx = geom.CellSize();

 MultiFab& DIV_new=get_new_data(DIV_Type,slab_step+1);
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int pcomp=AMREX_SPACEDIM;

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(DIV_new.boxArray().d_numPts());

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

  FArrayBox& volumefab=(*localMF[VOLUME_MF])[mfi];

  // coeff_avg,padvect_avg,coeff1,padvect1, ... ,coeffn,padvectn 
  FArrayBox& csoundfab=(*localMF[CELL_SOUND_MF])[mfi];
  FArrayBox& mdotfab=(*localMF[MDOT_MF])[mfi];
  FArrayBox& snewfab=S_new[mfi];
  FArrayBox& divfab=DIV_new[mfi];

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: NAVIERSTOKES_3D.F90
  FORT_UPDATE_DIV(
   xlo,dx,
   &dt_slab,
   volumefab.dataPtr(),
   ARLIM(volumefab.loVect()),ARLIM(volumefab.hiVect()),
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
 ns_reconcile_d_num(41);

} // subroutine ADVECT_DIV

void NavierStokes::getStateDIV_ALL(int idx,int ngrow) {

 if (level!=0)
  amrex::Error("level invalid getStateDIV_ALL");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

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

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 if (localMF_grow[idx]==-1) {
  new_localMF(idx,1,ngrow,-1);
 } else
  amrex::Error("local div data not previously deleted");

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int nsolve=1;

 int fluxvel_index=0;
 int fluxden_index=AMREX_SPACEDIM;

 MultiFab* velmac[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  velmac[dir]=getStateMAC(Umac_Type,0,dir,0,nsolve,cur_time_slab);

 const Real* dx = geom.CellSize();

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,250);
 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,80);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,610);

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>nmat))
  amrex::Error("nparts invalid");
 Vector<int> im_solid_map_null;
 im_solid_map_null.resize(1);

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=nmat)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  amrex::Error("nparts invalid");

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[idx]->boxArray().d_numPts());

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
   FArrayBox& az = (*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

   FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
   FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
   FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

   FArrayBox& maskcoef = (*localMF[MASKCOEF_MF])[mfi];//1=not cov  0=cov
   FArrayBox& maskSEMfab = (*localMF[MASKSEM_MF])[mfi];

   FArrayBox& vol = (*localMF[VOLUME_MF])[mfi];
   FArrayBox& ux = (*velmac[0])[mfi];
   FArrayBox& uy = (*velmac[1])[mfi];
   FArrayBox& uz = (*velmac[AMREX_SPACEDIM-1])[mfi];
   FArrayBox& rhs = (*localMF[idx])[mfi];
   FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
   FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
   FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
   const Real* xlo = grid_loc[gridno].lo();
   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

// RHS=(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)/vol_ij
//


   int operation_flag=1;
   int energyflag=0;
   int project_option=0;
   int homflag=0; // default
   int local_enable_spectral=enable_spectral;
   int use_VOF_weight=0;

   int ncomp_denold=vol.nComp();
   int ncomp_veldest=rhs.nComp();
   int ncomp_dendest=rhs.nComp();

   int tid_current=ns_thread();
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: NavierStokes::getStateDIV
   fort_mac_to_cell(
    &ns_time_order,
    &divu_outer_sweeps,
    &num_divu_outer_sweeps,
    &operation_flag, // operation_flag==1
    &energyflag,
    temperature_primitive_variable.dataPtr(),
    constant_density_all_time.dataPtr(),
    &nmat,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    added_weight.dataPtr(),
    &nten,
    &level, 
    &finest_level,
    &face_flag,
    &project_option,
    &local_enable_spectral,
    &fluxvel_index,
    &fluxden_index,
    &facevel_index,
    &facecut_index,
    &icefacecut_index,
    &curv_index,
    &conservative_tension_force,
    &conservative_div_uu,
    filter_velocity.dataPtr(),
    &ignore_div_up,
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
    reconfab.dataPtr(), //levelPC
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    solxfab.dataPtr(),
    ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
    solyfab.dataPtr(),
    ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
    solzfab.dataPtr(),
    ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
    reconfab.dataPtr(), //cterm
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    reconfab.dataPtr(), //pold
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    vol.dataPtr(), // denold
    ARLIM(vol.loVect()),ARLIM(vol.hiVect()),
    reconfab.dataPtr(), // ustar
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    reconfab.dataPtr(),
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    reconfab.dataPtr(), //mdotcell
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    reconfab.dataPtr(), //maskdivres
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    reconfab.dataPtr(), //maskres
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    &SDC_outer_sweeps,
    &homflag,
    &use_VOF_weight,
    &nsolve,
    &ncomp_denold,
    &ncomp_veldest,
    &ncomp_dendest,
    &SEM_advection_algorithm);
 } // mfi
} // omp
 ns_reconcile_d_num(42);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  delete velmac[dir];

} // subroutine getStateDIV



// initializes righthand side and mac_phi for doing jacobi method.
void NavierStokes::mac_project_rhs(int project_option,
 int idx_mac_phi_crse,int idx_mac_rhs_crse,int nsolve) {

 int finest_level=parent->finestLevel();

 if (level>finest_level)
  amrex::Error("level invalid mac_project_rhs");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn invalid76");

 if (localMF[idx_mac_phi_crse]->nComp()!=nsolve)
  amrex::Error("localMF[idx_mac_phi_crse]->nComp()  invalid");
 if (localMF[idx_mac_rhs_crse]->nComp()!=nsolve)
  amrex::Error("localMF[idx_mac_rhs_crse]->nComp()  invalid");

 localMF[idx_mac_phi_crse]->setVal(0.0,0,nsolve,1);
 localMF[idx_mac_phi_crse]->setBndry(0.0);

   // residual correction form,
   // p=dp+p^init, POLDHOLD=p^adv-p^{init}
   //         
   // a (p-p^adv) - div grad p = f
   // a (dp+p^init-p^adv) - div grad (dp+p^init) = f
   // a dp - div grad dp = div grad p^init + f + a(p^adv-p^init)
   //
   // rhs=POLDHOLD*a  - vol div u/dt + diffusionRHS
 int homflag=0;
 apply_div(
   project_option,homflag,
   POLDHOLD_MF,
   localMF[idx_mac_rhs_crse], 
   localMF[DIFFUSIONRHS_MF],
   UMAC_MF,
   nsolve);

}  // mac_project_rhs

// GRADPEDGE=-dt W grad p
// UMAC=UMAC+GRADPEDGE
// pnew+=mac_phi_crse
// POLDHOLD-=mac_phi_crse
void NavierStokes::mac_update(MultiFab* mac_phi_crse,int project_option,
  int nsolve) {

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option invalid77");

 if (mac_phi_crse->nGrow()!=1)
  amrex::Error("mac_phi_crse has invalid ngrow");
 if (mac_phi_crse->nComp()!=nsolve)
  amrex::Error("mac_phi_crse has invalid ncomp");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;
 int ncomp_check;

  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 
 if (ncomp_check!=nsolve)
  amrex::Error("ncomp_check invalid");

 MultiFab::Subtract(*localMF[POLDHOLD_MF],*mac_phi_crse,0,0,nsolve,0);

  // in: NavierStokes3.cpp
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
 if (scomp_ilist!=nsolve)
  amrex::Error("scomp_ilist invalid");
 
} // subroutine mac_update

// adjust tolerance if too stringent.
void NavierStokes::adjust_tolerance(Real& error0,Real& error0_max,
  int project_option) {

 if (project_option_is_valid(project_option)==1) {

  if (error0>error0_max)
   error0_max=error0;

  Real bot_rel_error=0.0;
  if (save_mac_abs_tol>0.0) {
   bot_rel_error=save_min_rel_error*save_atol_b/save_mac_abs_tol;
  } else
   amrex::Error("save_mac_abs_tol invalid");

  if (save_mac_abs_tol<save_min_rel_error*error0_max) {
    save_mac_abs_tol=1.01*save_min_rel_error*error0_max;
    save_atol_b=1.01*bot_rel_error*error0_max;
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "adjusting the tolerance project_option=" << 
       project_option << '\n';
     std::cout << "save_min_rel_error, error0_max " << 
       save_min_rel_error << ' ' <<
       error0_max << '\n';
     std::cout << "new tolerance: " << save_mac_abs_tol << '\n';
     std::cout << "new bottom tolerance: " << save_atol_b << '\n';
    } // ioproc?
  } // save_mac_abs_tol<save_min_rel_error*error0_max

 } else
  amrex::Error("project_option_is_valid() invalid");

}  // adjust_tolerance

}/* namespace amrex */
