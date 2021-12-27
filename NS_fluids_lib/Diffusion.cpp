// #include <winstd.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_CoordSys.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParmParse.H>

#include <NavierStokes.H>

#include <DIFFUSION_F.H>
#include <MACOPERATOR_F.H>
#include <LEVEL_F.H>
#include <GODUNOV_F.H>
#include <PROB_F.H>

#include <algorithm>

#include <cfloat>
#include <cmath>

namespace amrex{
	
void NavierStokes::diffuse_hoopALL(int idx_vel,int idx_thermal,
 int idx_force,int update_state) {

 int finest_level=parent->finestLevel();

//ux,vx,wx,uy,vy,wy,uz,vz,wz
// grad u in cylindrical coordinates:
//
// S= (grad u + grad u^T)/2 
//
// grad u=| u_r  u_t/r-v/r  u_z  |
//        | v_r  v_t/r+u/r  v_z  |
//        | w_r  w_t/r      w_z  |
//
// S=
//   |u_r     (u_t/r+v_r-v/r)/2   (u_z+w_r)/2   |
//   | .      v_t/r+u/r           (v_z+w_t/r)/2 |
//   | .      .                   w_z           |
//
// note: v_r-v/r=r(v/r)_r
//
// 2S=
//
//   |2u_r     (u_t/r+v_r-v/r)   (u_z+w_r)   |
//   | .       2v_t/r+2u/r       (v_z+w_t/r) |
//   | .      .                    2w_z      |
//
// 
// div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
//         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
//         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |
// 
// ur =     costheta u + sintheta v
// utheta = -sintheta u + costheta v
// 
// u = costheta ur - sintheta utheta
// v = sintheta ur + costheta utheta
//
// e.g. theta=pi/2  ur=0 
//   u=-utheta  v=0
// if constant viscosity:
// div u = (ru)_r/r + v_t/r + w_z= u_r + u/r +v_t/r + w_z=0
// (div u)_r=u_rr+u_r/r-u/r^2+v_tr/r-v_t/r^2+w_zr=0
// (div u)_t/r=u_rt/r+u_t/r^2+v_tt/r^2+w_zt/r=0
// (div u)_z=u_rz+u_z/r+v_tz/r+w_zz=0
//
// div(2 S)=
// |2u_rr+2u_r/r+u_tt/r^2+v_tr/r-v_t/r^2-2v_t/r^2-2u/r^2+u_zz+w_rz |
// |u_tr/r+v_rr+v_r/r-v_r/r+2v_tt/r^2+2u_t/r^2+u_t/r^2+v_r/r-v/r^2+v_zz+w_tz/r|
// |u_zr+u_z/r+w_rr+w_r/r+v_zt/r+w_tt/r^2 + 2w_zz |=
//
// |u_rr+u_r/r-u/r^2+u_tt/r^2-2v_t/r^2+u_zz |    
// |v_rr+v_r/r+v_tt/r^2+2u_t/r^2-v/r^2+v_zz |
// |w_rr+w_r/r+w_tt/r^2+w_zz                |
//
// compromise: 
//
// GU=| u_r       u_t/r  u_z  |
//    | v_r       v_t/r  v_z  |
//    | w_r       w_t/r  w_z  |
//
// hoop term 1st component:  -3 v_t/r^2 - 2 u/r^2
// hoop term 2nd component:   3 u_t/r^2 - v/r^2
// 
// If uncoupled_viscosity==true:
// hoop term 1st component:  -2 v_t/r^2 - u/r^2
// hoop term 2nd component:   2 u_t/r^2 - v/r^2
// No coupling terms.
// Diagonal terms not multiplied by 2.


 int simple_AMR_BC_flag_viscosity=1; 
 int do_alloc=0;
 init_gradu_tensorALL(
   idx_vel,
   do_alloc,
   CELLTENSOR_MF,
   FACETENSOR_MF,
   simple_AMR_BC_flag_viscosity);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.diffuse_hoop(idx_vel,idx_thermal,idx_force,update_state);
 } // ilev

 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

} // subroutine diffuse_hoopALL

// 1. compute explicit terms: u^* = u^n + dt F1(u^n)
// 2. compute implicit terms: u^** = u^* + dt F2(u^**)
void NavierStokes::diffuse_hoop(int idx_vel,int idx_thermal,
 int idx_force,int update_state) {
 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 60");

 Real gravity_normalized=0.0;
 gravity_normalized=std::abs(gravity);
 if (invert_gravity==1)
  gravity_normalized=-gravity_normalized;
 else if (invert_gravity==0) {
  // do nothing
 } else
  amrex::Error("invert_gravity invalid");


 const Real* dx = geom.CellSize();

 int nmat=num_materials;

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
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }

 int nsolve=AMREX_SPACEDIM;

 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;

 debug_ngrow(FACE_VAR_MF,0,810);

 debug_ngrow(idx_vel,1,812);
 if (localMF[idx_vel]->nComp()!=nsolve)
  amrex::Error("localMF[idx_vel]->nComp() invalid");

 debug_ngrow(idx_thermal,1,812);
 if (localMF[idx_thermal]->nComp()!=1)
  amrex::Error("localMF[idx_thermal]->nComp() invalid");

 debug_ngrow(idx_force,1,812);
 if (localMF[idx_force]->nComp()!=nsolve)
  amrex::Error("localMF[idx_force]->nComp() invalid");

 debug_ngrow(CELLTENSOR_MF,1,6);
 if (localMF[CELLTENSOR_MF]->nComp()!=ntensor)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 debug_ngrow(CELL_DEN_MF,1,811);
 debug_ngrow(CELL_VISC_MF,1,811);

 if (localMF[CELL_DEN_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");
 if (localMF[CELL_VISC_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_VISC_MF]->nComp() invalid");

 MultiFab* Un=localMF[idx_vel];

 MultiFab& U_new=get_new_data(State_Type,slab_step+1);
 int nstate=STATE_NCOMP;
 if (U_new.nComp()!=nstate) 
  amrex::Error("U_new.nComp()!=nstate");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)");

 Real local_visc_coef=visc_coef;

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(U_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(U_new,use_tiling); mfi.isValid(); ++mfi) {
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

  FArrayBox& uoldfab=(*Un)[mfi];
  FArrayBox& unewfab=U_new[mfi];
  FArrayBox& lsfab=LS_new[mfi];
  FArrayBox& denfab=(*localMF[CELL_DEN_MF])[mfi]; 
  FArrayBox& mufab=(*localMF[CELL_VISC_MF])[mfi];
  FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
  FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
  FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];
  FArrayBox& thermalfab=(*localMF[idx_thermal])[mfi];
  FArrayBox& forcefab=(*localMF[idx_force])[mfi];
  FArrayBox& tensorfab=(*localMF[CELLTENSOR_MF])[mfi];

  int tid_current=0;
#ifdef _OPENMP
  tid_current = omp_get_thread_num();
#endif
  if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
   // do nothing
  } else
   amrex::Error("tid_current invalid");

  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: DIFFUSION_3D.F90
  fort_hoopimplicit(
   override_density.dataPtr(), 
   constant_density_all_time.dataPtr(),
   &nstate,
   &gravity_normalized,
   &gravity_dir,
   forcefab.dataPtr(),
   ARLIM(forcefab.loVect()),ARLIM(forcefab.hiVect()),
   tensorfab.dataPtr(),
   ARLIM(tensorfab.loVect()),ARLIM(tensorfab.hiVect()),
   thermalfab.dataPtr(),
   ARLIM(thermalfab.loVect()),ARLIM(thermalfab.hiVect()),
   reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
   solxfab.dataPtr(),
   ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
   solyfab.dataPtr(),
   ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
   solzfab.dataPtr(),
   ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
   xlo,dx,
   uoldfab.dataPtr(),ARLIM(uoldfab.loVect()),ARLIM(uoldfab.hiVect()),
   unewfab.dataPtr(),ARLIM(unewfab.loVect()),ARLIM(unewfab.hiVect()),
   lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   mufab.dataPtr(),ARLIM(mufab.loVect()),ARLIM(mufab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &finest_level,
   &local_visc_coef,
   &angular_velocity,
   &uncoupled_viscosity,
   &update_state,
   &dt_slab,
   &rzflag,
   &nmat,
   &nparts,
   &nparts_def,
   im_solid_map_ptr,
   &ntensor,
   &nsolve);
 } // mfi
} // omp
 ns_reconcile_d_num(25);

}  // diffuse_hoop


// fluxes expected to approximate -dt k grad S
void NavierStokes::viscous_boundary_fluxes(
 int project_option,
 MultiFab* xflux,MultiFab* yflux,MultiFab* zflux,
 int nsolve) {
  
  bool use_tiling=ns_tiling;

  int nmat=num_materials;
  int nten=num_interfaces;

  if (project_option==SOLVETYPE_VISC) { 
   if (nsolve!=AMREX_SPACEDIM)
    amrex::Error("nsolve invalid");
  } else if (project_option==SOLVETYPE_HEAT) { 
   if (nsolve!=1)
    amrex::Error("nsolve invalid");
  } else
   amrex::Error("project_option not supported");

  const Real* dx = geom.CellSize();
     
  const Box& domain = geom.Domain();
  const int* domlo = domain.loVect();
  const int* domhi = domain.hiVect();

  VOF_Recon_resize(1,SLOPE_RECON_MF);
  debug_ngrow(SLOPE_RECON_MF,1,820);

  resize_levelsetLO(2,LEVELPC_MF);
  debug_ngrow(LEVELPC_MF,2,120);
  if (localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("levelpc mf has incorrect ncomp");

  resize_metrics(1);
 
  Vector<int> temp_dombc(2*AMREX_SPACEDIM);
  const BCRec& descbc = get_desc_lst()[State_Type].getBC(STATECOMP_STATES+1);
  const int* b_rec=descbc.vect();
  for (int m=0;m<2*AMREX_SPACEDIM;m++)
   temp_dombc[m]=b_rec[m];

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   MultiFab* fluxdir=nullptr;
   if (dir==0)
    fluxdir=xflux;
   else if (dir==1)
    fluxdir=yflux;
   else if ((dir==2)&&(AMREX_SPACEDIM==3))
    fluxdir=zflux;
   else
    amrex::Error("dir invalid viscous boundary fluxes");

   if (fluxdir->nComp()!=nsolve)
    amrex::Error("conserve_fluxes invalid ncomp");

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[LEVELPC_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[LEVELPC_MF],use_tiling);mfi.isValid(); ++mfi) {
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

    FArrayBox& xfluxfab=(*fluxdir)[mfi];
    FArrayBox& areafab=(*localMF[AREA_MF+dir])[mfi];
    FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];

    Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);
    Vector<int> tempbc=getBCArray(State_Type,gridno,STATECOMP_STATES+1,1);

    int tid_current=0;
#ifdef _OPENMP
    tid_current = omp_get_thread_num();
#endif
    if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
     // do nothing
    } else
     amrex::Error("tid_current invalid");

    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // if solidheat_flag=2, then inhomogeneous Neumann BC are
      // prescribed.
      // declared in: PROB.F90
    fort_viscfluxfill(
     macrolayer_size.dataPtr(),
     microlayer_substrate.dataPtr(),
     microlayer_temperature_substrate.dataPtr(),
     latent_heat.dataPtr(),
     freezing_model.dataPtr(),
     saturation_temp.dataPtr(),
     &nsolve,
     &dir,
     xlo,dx,
     velbc.dataPtr(), 
     tempbc.dataPtr(), 
     temp_dombc.dataPtr(),
     lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     areafab.dataPtr(),ARLIM(areafab.loVect()),ARLIM(areafab.hiVect()),
     xfluxfab.dataPtr(),
     ARLIM(xfluxfab.loVect()),ARLIM(xfluxfab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     domlo,domhi,
     &dt_slab,
     &nmat,
     &nten,
     &solidheat_flag,
     &project_option,
     &cur_time_slab);
   } // mfi
} //omp

   ns_reconcile_d_num(28);

  } // dir=0..sdim-1

} // subroutine viscous_boundary_fluxes

// project_option = SOLVETYPE_PRES (flux var=mac velocity)
// project_option = SOLVETYPE_HEAT (temp)
// project_option = SOLVETYPE_VISC (cell centered velocity) 
// project_option = SOLVETYPE_VISC,... (species)
// combine_flag==0 (FVM -> GFM) (T[im]=T_interp im=1..nmat)
// combine_flag==1 (GFM -> FVM)
// combine_flag==2 (combine if vfrac<VOFTOL)
// interface_cond_avail==1 if interface temperature and mass fraction are
//   available.
void NavierStokes::combine_state_variable(
 int project_option,
 int combine_idx,
 int combine_flag,
 int hflag,
 int update_flux,
 int interface_cond_avail) {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 int nmat=num_materials;

 if ((project_option==SOLVETYPE_PRES)||  // mac velocity
     (project_option==SOLVETYPE_VISC)) { // cell velocity
  if (interface_cond_avail==0) {
   // do nothing
  } else
   amrex::Error("interface_cond_avail invalid");
 } else if ((project_option==SOLVETYPE_HEAT)||   
            ((project_option>=SOLVETYPE_SPEC)&&
             (project_option<SOLVETYPE_SPEC+num_species_var))) {//mass fraction
  if ((interface_cond_avail==0)||
      (interface_cond_avail==1)) {
   // do nothing
  } else
   amrex::Error("interface_cond_avail invalid");
 } else
  amrex::Error("project_option invalid in combine_state_variable");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)");
 int nstate=STATE_NCOMP;
 if (S_new.nComp()!=nstate)
  amrex::Error("(S_new.nComp()!=nstate)");

 int nten=num_interfaces;
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (hydrate_flag==1) {
  if (num_species_var!=1)
   amrex::Error("num_species_var invalid");
 }

 int num_materials_combine=1;

 int nsolve=1;
 if ((project_option==SOLVETYPE_INITPROJ)||   
     (project_option==SOLVETYPE_PRESCOR)) { 
  amrex::Error("project_option invalid in combine_state_variable");
 } else if (project_option==SOLVETYPE_PRES) {    // regular projection
  nsolve=1;
  if (combine_flag!=2)
   amrex::Error("combine_flag invalid");
  num_materials_combine=1;
  if (update_flux!=1)
   amrex::Error("update_flux invalid");
 } else if (project_option==SOLVETYPE_HEAT) { 
  nsolve=1;
  num_materials_combine=nmat;
 } else if (project_option==SOLVETYPE_VISC) { 
  nsolve=AMREX_SPACEDIM;
  if (combine_flag!=2)  //combine if vfrac<VOFTOL
   amrex::Error("combine_flag invalid");
  num_materials_combine=1;
 } else if ((project_option>=SOLVETYPE_SPEC)&&
            (project_option<SOLVETYPE_SPEC+num_species_var)) { 
  nsolve=1;
  num_materials_combine=nmat;
 } else
  amrex::Error("project_option invalid81");

 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;
 int ncomp_check;
 get_mm_scomp_solver(
  num_materials_combine,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (state_index!=State_Type)
  amrex::Error("state_index invalid");

 if (ncomp_check!=nsolve*num_materials_combine) {
  std::cout << "nsolve= " << nsolve << '\n';
  std::cout << "num_materials_combine= " << num_materials_combine << '\n';
  std::cout << "ncomp_check= " << ncomp_check << '\n';
  print_project_option(project_option);
  std::cout << "scomp.size " << scomp.size() << '\n';
  std::cout << "ncomp.size " << ncomp.size() << '\n';
  std::cout << "scomp[0] " << scomp[0] << '\n';
  std::cout << "ncomp[0] " << ncomp[0] << '\n';
  amrex::Error("nsolve invalid 313");
 }

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
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }

 int ntsat=EXTRAP_NCOMP_TSAT;

 MultiFab* LEVEL_COMBINE=nullptr;

 MultiFab* STATE_INTERFACE=nullptr;

 if ((combine_flag==0)||  // FVM -> GFM
     (combine_flag==1)) { // GFM -> FVM

  resize_levelsetLO(2,LEVELPC_MF);
  LEVEL_COMBINE=localMF[LEVELPC_MF];

  if (interface_cond_avail==1) {
   if (is_phasechange==1) {
    STATE_INTERFACE=localMF[SATURATION_TEMP_MF];
    debug_ngrow(SATURATION_TEMP_MF,ngrow_make_distance,830);
    if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
     amrex::Error("localMF[SATURATION_TEMP_MF]->nComp()!=ntsat");
   } else if (is_phasechange==0) {
    STATE_INTERFACE=&LS_new;  // placeholder
   } else
    amrex::Error("is_phasechange invalid");
  } else if (interface_cond_avail==0) {
   STATE_INTERFACE=&LS_new;  // placeholder
  } else
   amrex::Error("interface_cond_avail invalid");


  debug_ngrow(LEVELPC_MF,2,830);

 } else if (combine_flag==2) { // combine if vfrac<VOFTOL

  STATE_INTERFACE=&LS_new; // placeholder
  if (update_flux==0) {
   LEVEL_COMBINE=&LS_new;
  } else if (update_flux==1) {
   resize_levelsetLO(2,LEVELPC_MF);
   LEVEL_COMBINE=localMF[LEVELPC_MF];
   debug_ngrow(LEVELPC_MF,2,830);
  } else
   amrex::Error("update_flux invalid");

 } else
  amrex::Error("combine_flag invalid");
  
 if (LEVEL_COMBINE->nGrow()<1) {
  std::cout << "LEVEL_COMBINE->nGrow()= " << LEVEL_COMBINE->nGrow() << '\n';
  std::cout << "combine_flag= " << combine_flag << '\n';
  std::cout << "update_flux= " << update_flux << '\n';
  amrex::Error("LEVEL_COMBINE->nGrow()<1");
 }
 if (LEVEL_COMBINE->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LEVEL_COMBINE->nComp()!=nmat*(AMREX_SPACEDIM+1)");

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,832);
 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,830);
 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,6001);

 MultiFab* signed_distance=nullptr;
 if ((combine_flag==0)||  // FVM -> GFM
     (combine_flag==1)) { // GFM -> FVM

  signed_distance=LEVEL_COMBINE;
  if (signed_distance->nComp()!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("signed_distance->nComp() invalid");

 } else if (combine_flag==2) { // combine if F==0

  signed_distance=localMF[SLOPE_RECON_MF];
  if (signed_distance->nComp()!=nmat*ngeom_recon)
   amrex::Error("signed_distance->nComp() invalid");

 } else
  amrex::Error("combine_flag invalid");

 if ((signed_distance->nGrow()==1)||
     (signed_distance->nGrow()==2)) {
  // do nothing
 } else {
  std::cout << "signed_distance->nGrow()= " <<
   signed_distance->nGrow() << '\n';
  amrex::Error("signed_distance->nGrow() invalid");
 }
  
 if ((update_flux==0)&&(combine_idx==-1)&&
     ((combine_flag==0)||(combine_flag==1))) {
  init_boundary_list(scomp,ncomp);
 } else if ((update_flux==1)||(combine_idx>=0)||(combine_flag==2)) {
  // do nothing
 } else
  amrex::Error("update_flux or combine_idx invalid");

 const Real* dx = geom.CellSize();

 if (update_flux==1) {

  if (combine_flag==2) {
   // do nothing
  } else
   amrex::Error("combine_flag invalid");

  if (project_option==SOLVETYPE_PRES) { // mac velocity

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

    debug_ngrow(FACE_VAR_MF+dir,0,832);

    MultiFab* face_mf=nullptr;

    if (combine_idx==-1) {

     if (project_option==SOLVETYPE_PRES) {
      face_mf=&get_new_data(Umac_Type+dir,slab_step+1);
     } else
      amrex::Error("project_option invalid82");

    } else if (combine_idx>=0) {

     face_mf=localMF[combine_idx];

    } else
     amrex::Error("combine_idx invalid");

    if (face_mf->nComp()!=1)
     amrex::Error("face_mf->nComp() invalid");

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
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

     FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];  

     FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
     FArrayBox& macfab=(*face_mf)[mfi];
     FArrayBox& lsfab=(*LEVEL_COMBINE)[mfi];
     FArrayBox& solfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];
     Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

     int tid_current=ns_thread();
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
 
      // declared in: GODUNOV_3D.F90
     fort_combinevelface(
      &tid_current,
      &hflag,
      &nmat,
      &nparts,
      &nparts_def,
      im_solid_map_ptr,
      &nten,
      &combine_idx,
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      &level,
      &finest_level,
      velbc.dataPtr(),
      voffab.dataPtr(),
      ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
      macfab.dataPtr(),
      ARLIM(macfab.loVect()),ARLIM(macfab.hiVect()),
      xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
      lsfab.dataPtr(),
      ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
      solfab.dataPtr(),
      ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
      xlo,dx,
      &dir,
      &cur_time_slab);

    }  // mfi
} // omp
    ns_reconcile_d_num(29);
   } // dir=0..sdim-1

  } else
   amrex::Error("project_option invalid");

 } else if (update_flux==0) {

  MultiFab* cell_mf=nullptr;
  if (combine_idx==-1) {
   cell_mf=&S_new;
   if (cell_mf->nComp()!=nstate)
    amrex::Error("cell_mf->nComp() invalid");
  } else if (combine_idx>=0) {
   cell_mf=localMF[combine_idx];
   if (cell_mf->nComp()!=nsolve*num_materials_combine)
    amrex::Error("cell_mf->nComp() invalid");
  } else
   amrex::Error("combine_idx invalid");

  int ncomp_cell=cell_mf->nComp();

  if (combine_flag==1) {  // center -> centroid

   if (combine_idx!=-1)
    amrex::Error("combine_idx invalid");

    // GFM solver only solves for first component.
   for (int im=1;im<nmat;im++) 
    MultiFab::Copy(*cell_mf,*cell_mf,scomp[0],scomp[im],1,0);

  } else if (combine_flag==0) { // centroid -> center
   // do nothing
  } else if (combine_flag==2) { // init zero VFRAC areas.
   // do nothing
  } else
   amrex::Error("combine_flag invalid");

  MultiFab* new_combined;
  if ((combine_flag==0)||  // FVM -> GFM
      (combine_flag==1)) { // GFM -> FVM
   if (nsolve!=1)
    amrex::Error("nsolve invalid");
   new_combined=new MultiFab(grids,dmap,nsolve*nmat,0,
     MFInfo().SetTag("new_combined"),FArrayBoxFactory());
   if (combine_idx==-1) {
    for (int im=0;im<nmat;im++) 
     MultiFab::Copy(*new_combined,*cell_mf,scomp[im],im,1,0);
   } else if (combine_idx>=0) {
    MultiFab::Copy(*new_combined,*cell_mf,0,0,nmat,0);
   } else
    amrex::Error("combine_idx invalid");
  } else if (combine_flag==2) {
   new_combined=cell_mf;
  } else
   amrex::Error("combine_flag invalid");

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
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

    // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& levelpcfab=(*signed_distance)[mfi];
   FArrayBox& lsnewfab=LS_new[mfi];
   FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& cellfab=(*cell_mf)[mfi];
   FArrayBox& newcell=(*new_combined)[mfi];
   FArrayBox& statefab=S_new[mfi];

   // ntsat=nten*(ncomp_per_tsat+1)
   // e.g. for interface 12,
   //  component 1=0 if T_gamma,Y_gamma not defined
   //             =1 if T_gamma,Y_gamma defined in RATEMASSCHANGE
   //             =2 if T_gamma,Y_gamma defined after extrapolation
   //             =-1 or -2 for condensation case.
   //  component 2=T_gamma
   //  component 3=Y_gamma
   //  repeats ....
   FArrayBox& Tsatfab=(*STATE_INTERFACE)[mfi];

   if ((combine_flag==0)||
       (combine_flag==1)) {
    if ((is_phasechange==1)&&(interface_cond_avail==1)) {
     if (Tsatfab.nComp()==ntsat) {
      //do nothing
     } else {
      std::cout << "Tsatfab.nComp()=" << Tsatfab.nComp() << 
        " ntsat=" << ntsat << '\n';
      amrex::Error("Tsatfab.nComp()!=ntsat 1");
     }
    } else if ((is_phasechange==0)||
               (interface_cond_avail==0)) {
     // do nothing
    } else
     amrex::Error("is_phasechange or interface_cond_avail invalid");
   } else if (combine_flag==2) {
    // do nothing
   } else
    amrex::Error("combine_flag invalid");

   FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
   FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
   FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);
   Vector<int> listbc;
   getBCArray_list(listbc,state_index,gridno,scomp,ncomp);
   if (listbc.size()!=nsolve*num_materials_combine*AMREX_SPACEDIM*2)
    amrex::Error("listbc.size() invalid");

   int scomp_size=scomp.size();

   int tid_current=ns_thread();
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: GODUNOV_3D.F90
   fort_combinevel(
    &tid_current,
    &hflag,
    &num_materials_combine,
    mass_fraction_id.dataPtr(),
    latent_heat.dataPtr(),
    freezing_model.dataPtr(),
    Tanasawa_or_Schrage_or_Kassemi.dataPtr(),
    distribute_from_target.dataPtr(),
    saturation_temp.dataPtr(),
    &hydrate_flag,
    &nmat,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    &nten,
    &nsolve,
    &project_option,
    &combine_idx,
    &combine_flag,
    &interface_cond_avail,
    &nstate,
    &ncomp_cell,
    scomp.dataPtr(),
    ncomp.dataPtr(),
    &scomp_size,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    &finest_level,
    &ntsat,
    Tsatfab.dataPtr(),
    ARLIM(Tsatfab.loVect()),ARLIM(Tsatfab.hiVect()),
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    solxfab.dataPtr(),
    ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
    solyfab.dataPtr(),
    ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
    solzfab.dataPtr(),
    ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
    lsnewfab.dataPtr(),
    ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
    voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    cellfab.dataPtr(),ARLIM(cellfab.loVect()),ARLIM(cellfab.hiVect()),
    newcell.dataPtr(),ARLIM(newcell.loVect()),ARLIM(newcell.hiVect()),
    statefab.dataPtr(),
    ARLIM(statefab.loVect()),ARLIM(statefab.hiVect()),
    velbc.dataPtr(),
    listbc.dataPtr(),
    xlo,dx,
    &cur_time_slab);

  }  // mfi
} // omp
  ns_reconcile_d_num(30);

  if ((combine_flag==0)||  // FVM->GFM
      (combine_flag==1)) { // GFM->FVM
   if (combine_idx==-1) {
    for (int im=0;im<nmat;im++)
     MultiFab::Copy(*cell_mf,*new_combined,im,scomp[im],1,0);
   } else if (combine_idx>=0) {
    MultiFab::Copy(*cell_mf,*new_combined,0,0,nmat,0);
   } else
    amrex::Error("combine_idx invalid");

   delete new_combined;

  } else if (combine_flag==2) {
   // do nothing
  } else
   amrex::Error("combine_flag invalid");

  if ((combine_idx==-1)&&
      (hflag==0)) {

    // if level<finest_level, then average down from level+1.
   int spectral_override=1;
   avgDown_list(state_index,scomp,ncomp,spectral_override); 

  } else if ((combine_idx==0)||
             (hflag==1)) {

   // do nothing

  } else
   amrex::Error("combine_idx or hflag invalid");

 } else
  amrex::Error("update_flux invalid"); 

}  // combine_state_variable


// viscous heating (viscoelastic heating not done here)
void NavierStokes::diffusion_heating(int source_idx,int idx_heat) {
 
 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;
 int nden=nmat*num_state_material;
 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;

 int nsolve=AMREX_SPACEDIM;

 debug_ngrow(FACE_VAR_MF,0,2);
 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,3);
 debug_ngrow(idx_heat,0,4);
 debug_ngrow(source_idx,1,842);
 debug_ngrow(CELLTENSOR_MF,1,6);
 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,5);
 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,845);

 if (localMF[source_idx]->nComp()!=nsolve)
  amrex::Error("localMF[source_idx]->nComp() invalid");

 if (localMF[idx_heat]->nComp()!=1)
  amrex::Error("localMF[idx_heat]->nComp()!=1");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(CONSERVE_FLUXES_MF+dir,0,7);
  if (localMF[CONSERVE_FLUXES_MF+dir]->nComp()!=nsolve)
   amrex::Error("localMF[CONSERVE_FLUXES_MF+dir]->nComp() invalid");
  setVal_localMF(CONSERVE_FLUXES_MF+dir,0.0,0,nsolve,0);
 }
 int homflag=0;
 int energyflag=0;
 int project_option=SOLVETYPE_VISC; 
 int simple_AMR_BC_flag=1; 
 int simple_AMR_BC_flag_viscosity=1; 
 apply_pressure_grad(
  simple_AMR_BC_flag, 
  simple_AMR_BC_flag_viscosity, 
  homflag,energyflag,CONSERVE_FLUXES_MF,
  source_idx,
  project_option,nsolve);

 int ncomp_edge=-1;
 avgDownEdge_localMF(
  CONSERVE_FLUXES_MF,
  0,ncomp_edge,0,AMREX_SPACEDIM,1,21);

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 2");

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[CELLTENSOR_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[CELLTENSOR_MF],use_tiling); mfi.isValid(); ++mfi) {
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

  FArrayBox& xstress=(*localMF[CONSERVE_FLUXES_MF])[mfi];
  FArrayBox& ystress=(*localMF[CONSERVE_FLUXES_MF+1])[mfi];
  FArrayBox& zstress=(*localMF[CONSERVE_FLUXES_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& DeDTinversefab=(*localMF[CELL_DEDT_MF])[mfi]; // 1/(rho cv)
  if (DeDTinversefab.nComp()!=1)
   amrex::Error("DeDTinversefab.nComp() invalid");

  FArrayBox& gradufab=(*localMF[CELLTENSOR_MF])[mfi];
  if (gradufab.nComp()!=ntensor)
   amrex::Error("gradufab.nComp() invalid");

  FArrayBox& heatfab=(*localMF[idx_heat])[mfi];
  FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];

  Real local_dt_slab=dt_slab;
  if (hold_dt_factors[0]==1.0) {
   // do nothing
  } else if ((hold_dt_factors[0]>0.0)&&
             (hold_dt_factors[0]<1.0)) {
   local_dt_slab*=hold_dt_factors[0];
  } else
   amrex::Error("hold_dt_factors[0] invalid");

  int tid_current=0;
#ifdef _OPENMP
  tid_current = omp_get_thread_num();
#endif
  if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
   // do nothing
  } else
   amrex::Error("tid_current invalid");

  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   //declared in: GODUNOV_3D.F90
  fort_visctensorheat(
   &ntensor,
   &nsolve,
   &nstate,
   xlo,dx,
   lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
   DeDTinversefab.dataPtr(),
   ARLIM(DeDTinversefab.loVect()),ARLIM(DeDTinversefab.hiVect()),
   heatfab.dataPtr(),
   ARLIM(heatfab.loVect()),ARLIM(heatfab.hiVect()),
   xstress.dataPtr(),
   ARLIM(xstress.loVect()),ARLIM(xstress.hiVect()),
   ystress.dataPtr(),
   ARLIM(ystress.loVect()),ARLIM(ystress.hiVect()),
   zstress.dataPtr(),
   ARLIM(zstress.loVect()),ARLIM(zstress.hiVect()),
   gradufab.dataPtr(),
   ARLIM(gradufab.loVect()),ARLIM(gradufab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,&level,
   &local_dt_slab,&rzflag,&nmat,&nden);
 }  // mfi  
} // omp

 ns_reconcile_d_num(31);

}   // subroutine diffusion_heating

}/* namespace amrex */
