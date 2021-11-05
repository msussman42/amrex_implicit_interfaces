// get rid of autoindent   :setl noai nocin nosi inde=
#include <AMReX_FArrayBox.H>
#include <AMReX_CoordSys.H>
#include <AMReX_ParmParse.H>

#include <NavierStokes.H>
#include <PROB_F.H>
#include <DERIVE_F.H>
#include <NAVIERSTOKES_F.H>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

namespace amrex{

//
// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
//


// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int norm_vel_bc[] =
{ INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, EXT_DIR, EXT_DIR };

static int norm_vel_extrap_bc[] =
{ INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_ODD, FOEXTRAP, FOEXTRAP };

static int tang_vel_bc[] =
{ INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, EXT_DIR };

static int tang_vel_extrap_bc[] =
{ INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP };

// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int grad_dot_n_norm_vel_bc[] =
{ INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP };
static int grad_dot_t_norm_vel_bc[] =
{ INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_ODD, FOEXTRAP, FOEXTRAP };

static int grad_dot_n_tang_vel_bc[] =
{ INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_ODD, FOEXTRAP, FOEXTRAP };
static int grad_dot_t_tang_vel_bc[] =
{ INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP };


// Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int scalar_bc[] =
{ INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP };

// Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int temperature_bc[] =
{ INT_DIR, EXT_DIR, EXT_DIR, REFLECT_EVEN, FOEXTRAP, FOEXTRAP };

// Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int custom_temperature_bc[] =
{ INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, EXT_DIR, EXT_DIR };

// Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int custom2_temperature_bc[] =
{ INT_DIR, EXT_DIR, EXT_DIR, REFLECT_EVEN, EXT_DIR, EXT_DIR };

// Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int scalar_vof_bc[] =
{ INT_DIR, EXT_DIR, EXT_DIR, REFLECT_EVEN, EXT_DIR, EXT_DIR };

// Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int press_bc[] =
{ INT_DIR, FOEXTRAP, EXT_DIR, REFLECT_EVEN, FOEXTRAP, FOEXTRAP };

// Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
static int extrap_bc[] =
{ INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP };

static
void
set_tensor_bc (BCRec&       bc,
               const BCRec& phys_bc,int dir1,int dir2)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
   
     // if s=REFLECT EVEN => grad s dot n=REFLECT ODD
     //                      grad s dot t=REFLECT EVEN
     // if s=FOEXTRAP     => grad s dot n=FOEXTRAP
     //                      grad s dot t=FOEXTRAP
     // if s=REFLECT ODD  => grad s dot n=REFLECT EVEN
     //                      grad s dot t=REFLECT_ODD
     // if s=EXT_DIR => grad s dot n=FOEXTRAP
     //                 grad s dot t=FOEXTRAP
     // inflow:     u dot n=EXT_DIR,     u dot t=EXT_DIR
     // outflow:    u dot n=FOEXTRAP,    u dot t=FOEXTRAP
     // symmetry:   u dot n=REFLECT_ODD, u dot t=REFLECT_EVEN
     // SlipWall:   u dot n=EXT_DIR,     u dot t=FOEXTRAP
     // NoSlipWall: u dot n=EXT_DIR,     u dot t=EXT_DIR
     //
     // grad dot n:
     // inflow:     u dot n=EXT_DIR      -> FOEXTRAP
     // outflow:    u dot n=FOEXTRAP,    -> FOEXTRAP
     // symmetry:   u dot n=REFLECT_ODD, -> REFLECT_EVEN
     // SlipWall:   u dot n=EXT_DIR,     -> FOEXTRAP
     // NoSlipWall: u dot n=EXT_DIR,     -> FOEXTRAP
     // grad dot t:
     // inflow:     u dot n=EXT_DIR      -> FOEXTRAP
     // outflow:    u dot n=FOEXTRAP,    -> FOEXTRAP
     // symmetry:   u dot n=REFLECT_ODD, -> REFLECT_ODD
     // SlipWall:   u dot n=EXT_DIR,     -> FOEXTRAP
     // NoSlipWall: u dot n=EXT_DIR,     -> FOEXTRAP

     //
     // grad dot n:
     // inflow:     u dot t=EXT_DIR      -> FOEXTRAP
     // outflow:    u dot t=FOEXTRAP,    -> FOEXTRAP
     // symmetry:   u dot t=REFLECT_EVEN -> REFLECT_ODD
     // SlipWall:   u dot t=FOEXTRAP,    -> FOEXTRAP
     // NoSlipWall: u dot t=EXT_DIR,     -> FOEXTRAP
     // grad dot t:
     // inflow:     u dot t=EXT_DIR      -> FOEXTRAP
     // outflow:    u dot t=FOEXTRAP,    -> FOEXTRAP
     // symmetry:   u dot t=REFLECT_EVEN -> REFLECT_EVEN
     // SlipWall:   u dot t=FOEXTRAP,    -> FOEXTRAP
     // NoSlipWall: u dot t=EXT_DIR,     -> FOEXTRAP

    if ((dir1>=0)&&(dir1<AMREX_SPACEDIM)) {
     // do nothing
    } else
     amrex::Error("dir1 invalid");
    if ((dir2>=0)&&(dir2<AMREX_SPACEDIM)) {
     // do nothing
    } else
     amrex::Error("dir2 invalid");

    for (int dir3=0;dir3<AMREX_SPACEDIM;dir3++) {

     if ((dir1==dir3)&&(dir2!=dir3)) {
      bc.setLo(dir3,grad_dot_t_norm_vel_bc[lo_bc[dir3]]);
      bc.setHi(dir3,grad_dot_t_norm_vel_bc[hi_bc[dir3]]);
     } else if ((dir1!=dir3)&&(dir2==dir3)) {
      bc.setLo(dir3,grad_dot_n_tang_vel_bc[lo_bc[dir3]]);
      bc.setHi(dir3,grad_dot_n_tang_vel_bc[hi_bc[dir3]]);
     } else if ((dir1!=dir3)&&(dir2!=dir3)) {
      bc.setLo(dir3,grad_dot_t_tang_vel_bc[lo_bc[dir3]]);
      bc.setHi(dir3,grad_dot_t_tang_vel_bc[hi_bc[dir3]]);
     } else if ((dir1==dir3)&&(dir2==dir3)) {
      bc.setLo(dir3,grad_dot_n_norm_vel_bc[lo_bc[dir3]]);
      bc.setHi(dir3,grad_dot_n_norm_vel_bc[hi_bc[dir3]]);
     } else
      amrex::Error("dir1,dir2,dir3 bust");

    } // dir3=0..sdim-1

} // subroutine set_tensor_bc


// in RZ: the hoop term is 2 u/r
//  as r->0, u/r ->u_r u_r has REFLECT_EVEN BC
static
void
set_hoop_bc (BCRec& bc,const BCRec& phys_bc)
{
 if (AMREX_SPACEDIM==2) {

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
    //dir1="theta"  dir2="theta"  dir3="r" 
  bc.setLo(0,grad_dot_t_tang_vel_bc[lo_bc[0]]); //REFLECT_EVEN if symmetric
  bc.setHi(0,grad_dot_t_tang_vel_bc[hi_bc[0]]);
    //dir1="theta"  dir2="theta"  dir3="z" 
  bc.setLo(1,grad_dot_t_tang_vel_bc[lo_bc[1]]); //REFLECT_EVEN if symmetric
  bc.setHi(1,grad_dot_t_tang_vel_bc[hi_bc[1]]);

 } else
  amrex::Error("bl_spacedim invalid");

} // subroutine set_hoop_bc


static
void
set_x_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static void
set_y_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}


static
void
set_x_vel_extrap_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_extrap_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_extrap_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_extrap_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_extrap_bc[hi_bc[1]]);
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,tang_vel_extrap_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_extrap_bc[hi_bc[2]]);
#endif
}

static void
set_y_vel_extrap_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_extrap_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_extrap_bc[hi_bc[0]]);
    bc.setLo(1,norm_vel_extrap_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_extrap_bc[hi_bc[1]]);
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,tang_vel_extrap_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_extrap_bc[hi_bc[2]]);
#endif
}


static
void
set_z_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
#endif
}


static
void
set_z_vel_extrap_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_extrap_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_extrap_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_extrap_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_extrap_bc[hi_bc[1]]);
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,norm_vel_extrap_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_extrap_bc[hi_bc[2]]);
#endif
}

static
void
set_scalar_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_bc[lo_bc[i]]);
        bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}


static
void
set_custom_temperature_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,custom_temperature_bc[lo_bc[i]]);
        bc.setHi(i,custom_temperature_bc[hi_bc[i]]);
    }
}


static
void
set_custom2_temperature_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,custom2_temperature_bc[lo_bc[i]]);
        bc.setHi(i,custom2_temperature_bc[hi_bc[i]]);
    }
}



static
void
set_temperature_bc (BCRec&       bc,
           const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,temperature_bc[lo_bc[i]]);
        bc.setHi(i,temperature_bc[hi_bc[i]]);
    }
}

static
void
set_extrap_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo(); 
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,extrap_bc[lo_bc[i]]);
        bc.setHi(i,extrap_bc[hi_bc[i]]);
    }
}


static
void
set_scalar_vof_bc (BCRec&       bc,
                   const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_vof_bc[lo_bc[i]]);
        bc.setHi(i,scalar_vof_bc[hi_bc[i]]);
    }
}



static
void
set_pressure_bc (BCRec&       bc,
                 const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i,press_bc[lo_bc[i]]);
        bc.setHi(i,press_bc[hi_bc[i]]);
    }
}

void
NavierStokes::override_enable_spectral(int enable_spectral_in) {

 int nmat=num_materials;
 int scomp_states=AMREX_SPACEDIM+1;
 int ncomp_vel_pres=AMREX_SPACEDIM+1;

 enable_spectral=enable_spectral_in;

 if ((enable_spectral_in==1)||(enable_spectral_in==2)) {

  sem_interp_HIGH_PARM.interp_enable_spectral=enable_spectral_in;
  for (int im=0;im<nmat;im++) {
   int ibase=im*num_state_material;
   desc_lst.resetMapper(State_Type,scomp_states+ibase+DenVar,
     &sem_interp_HIGH_PARM);
   desc_lst.resetMapper(State_Type,scomp_states+ibase+TemperatureVar,
     &sem_interp_HIGH_PARM);
  } // im=0..nmat-1

  for (int imvel=0;imvel<ncomp_vel_pres;imvel++) {
   desc_lst.resetMapper(State_Type,imvel,&sem_interp_HIGH_PARM);
  }

 } else if ((enable_spectral_in==0)||(enable_spectral_in==3)) {

  sem_interp_LOW_PARM.interp_enable_spectral=enable_spectral_in;
  for (int im=0;im<nmat;im++) {
   int ibase=im*num_state_material;
   desc_lst.resetMapper(State_Type,scomp_states+ibase+DenVar,
     &sem_interp_LOW_PARM);
   desc_lst.resetMapper(State_Type,scomp_states+ibase+TemperatureVar,
     &sem_interp_LOW_PARM);
  } // im=0..nmat-1

  for (int imvel=0;imvel<ncomp_vel_pres;imvel++) {
   desc_lst.resetMapper(State_Type,imvel,&sem_interp_LOW_PARM);
  }
 } else
  amrex::Error("enable_spectral_in invalid");

}  // subroutine override_enable_spectral

void
NavierStokes::override_enable_spectralGHOST(
  int scomp,int ncomp,int enable_spectral_in) {

 if ((scomp>=0)&&(ncomp>=1)) {

  for (int nc=0;nc<ncomp;nc++) {
   int target_comp=scomp+nc;
   if ((enable_spectral_in==0)||(enable_spectral_in==3)) {
    desc_lstGHOST.resetMapper(State_Type,target_comp,&pc_interp);
   } else if ((enable_spectral_in==1)||(enable_spectral_in==2)) {
    sem_interp_HIGH_PARM.interp_enable_spectral=enable_spectral_in;
    desc_lstGHOST.resetMapper(State_Type,target_comp,&sem_interp_HIGH_PARM);
   } else
    amrex::Error("enable_spectral_in invalid");
  } // nc=0..ncomp-1

 } else
  amrex::Error("scomp or ncomp invalid");

}  // subroutine override_enable_spectralGHOST

// called from: NavierStokes::allocate_levelsetLO
void
NavierStokes::override_LS_HO(int Interp_LO) { // 0=use normals 1=PC

 if ((Interp_LO!=0)&&(Interp_LO!=1))
  amrex::Error("Interp_LO invalid");

 int nmat=num_materials;
 int ncomp_ls=(AMREX_SPACEDIM+1)*nmat;

 if (Interp_LO==0) {
  ls_ho_interp_HIGH_PARM.LSHOInterp_nmat=nmat;
  ls_ho_interp_HIGH_PARM.LSHOInterp_LO=Interp_LO;
  for (int im=0;im<ncomp_ls;im++) {
   desc_lst.resetMapper(LS_Type,im,&ls_ho_interp_HIGH_PARM);
  } // im
 } else if (Interp_LO==1) {
  ls_ho_interp_LOW_PARM.LSHOInterp_nmat=nmat;
  ls_ho_interp_LOW_PARM.LSHOInterp_LO=Interp_LO;
  for (int im=0;im<ncomp_ls;im++) {
   desc_lst.resetMapper(LS_Type,im,&ls_ho_interp_LOW_PARM);
  } // im
 } else
  amrex::Error("Interp_LO invalid");

}  // subroutine override_LS_HO

void
NavierStokes::set_tensor_extrap_components(
  int coord,std::string postfix,int indx,int ibase_tensor) {

 BCRec bc;

 if (ibase_tensor>=0) {
  // do nothing
 } else
  amrex::Error("ibase_tensor invalid");

 int ibase_tensor_local=ibase_tensor;

 if (NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
  // do nothing
 } else
  amrex::Error("NUM_TENSOR_TYPE invalid");

  // no EXT_DIR BCs
 set_tensor_bc(bc,phys_bc,0,0);
 std::string T11_strE="T11extrap"+postfix; 
  // low order extrapolation (if EXT_DIR BCs were present)
 desc_lstGHOST.setComponent(indx,ibase_tensor_local,
   T11_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

 ibase_tensor_local++;
     
  // no EXT_DIR BCs
 set_tensor_bc(bc,phys_bc,0,1);
 std::string T12_strE="T12extrap"+postfix; 
 desc_lstGHOST.setComponent(indx,ibase_tensor_local,
   T12_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

 ibase_tensor_local++;
     
  // no EXT_DIR BCs
 set_tensor_bc(bc,phys_bc,1,1);
 std::string T22_strE="T22extrap"+postfix; 
 desc_lstGHOST.setComponent(indx,ibase_tensor_local,
   T22_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

 ibase_tensor_local++;
    
 if (AMREX_SPACEDIM==2) {
  if ((CoordSys::CoordType) coord == CoordSys::RZ) {
    // no EXT_DIR BCs
   set_hoop_bc(bc,phys_bc);
  } else if ((CoordSys::CoordType) coord == CoordSys::cartesian) {
   // placeholder: Q33 should always be 0
   // no EXT_DIR BCs
   set_hoop_bc(bc,phys_bc);
  } else if ((CoordSys::CoordType) coord == CoordSys::CYLINDRICAL) {
   // placeholder: Q33 should always be 0
   // no EXT_DIR BCs
   set_hoop_bc(bc,phys_bc);
  } else
   amrex::Error("(CoordSys::CoordType) coord invalid");
 } else if (AMREX_SPACEDIM==3) {
  if ((CoordSys::CoordType) coord == CoordSys::cartesian) {
   // no EXT_DIR BCs
   set_tensor_bc(bc,phys_bc,2,2);
  } else if ((CoordSys::CoordType) coord == CoordSys::CYLINDRICAL) {
   // no EXT_DIR BCs
   set_tensor_bc(bc,phys_bc,2,2);
  } else
   amrex::Error("(CoordSys::CoordType) coord invalid");
  } else
   amrex::Error("sdim invalid");
 
  std::string T33_strE="T33extrap"+postfix; 
  desc_lstGHOST.setComponent(indx,ibase_tensor_local,
    T33_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

#if (AMREX_SPACEDIM == 3)
  ibase_tensor_local++;

   // no EXT_DIR BCs
  set_tensor_bc(bc,phys_bc,0,2);
  std::string T13_strE="T13extrap"+postfix; 
  desc_lstGHOST.setComponent(indx,ibase_tensor_local,
    T13_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);
     
  ibase_tensor_local++;
     
   // no EXT_DIR BCs
  set_tensor_bc(bc,phys_bc,1,2);
  std::string T23_strE="T23extrap"+postfix; 
  desc_lstGHOST.setComponent(indx,ibase_tensor_local,
    T23_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);
#endif

  if (ibase_tensor_local==ibase_tensor+NUM_TENSOR_TYPE-1) {
   // do nothing
  } else {
   std::cout << "ibase_tensor_local=" << ibase_tensor_local << '\n';
   amrex::Error("ibase_tensor_local invalid");
  }

} // end subroutine set_tensor_extrap_components


void
NavierStokes::set_tensor_extrap_components_main(
  int coord,std::string postfix,int indx) {

 BCRec bc;

 int ibase_tensor_local=0;

 if (NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
  // do nothing
 } else
  amrex::Error("NUM_TENSOR_TYPE invalid");

  // no EXT_DIR BCs
 set_tensor_bc(bc,phys_bc,0,0);
 std::string T11_strE="T11main"+postfix; 
  // low order extrapolation (if EXT_DIR BCs were present)
 desc_lst.setComponent(indx,ibase_tensor_local,
   T11_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

 ibase_tensor_local++;
     
  // no EXT_DIR BCs
 set_tensor_bc(bc,phys_bc,0,1);
 std::string T12_strE="T12main"+postfix; 
 desc_lst.setComponent(indx,ibase_tensor_local,
   T12_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

 ibase_tensor_local++;
     
  // no EXT_DIR BCs
 set_tensor_bc(bc,phys_bc,1,1);
 std::string T22_strE="T22main"+postfix; 
 desc_lst.setComponent(indx,ibase_tensor_local,
   T22_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

 ibase_tensor_local++;
    
 if (AMREX_SPACEDIM==2) {
  if ((CoordSys::CoordType) coord == CoordSys::RZ) {
    // no EXT_DIR BCs
   set_hoop_bc(bc,phys_bc);
  } else if ((CoordSys::CoordType) coord == CoordSys::cartesian) {
   // placeholder: Q33 should always be 0
   // no EXT_DIR BCs
   set_hoop_bc(bc,phys_bc);
  } else if ((CoordSys::CoordType) coord == CoordSys::CYLINDRICAL) {
   // placeholder: Q33 should always be 0
   // no EXT_DIR BCs
   set_hoop_bc(bc,phys_bc);
  } else
   amrex::Error("(CoordSys::CoordType) coord invalid");
 } else if (AMREX_SPACEDIM==3) {
  if ((CoordSys::CoordType) coord == CoordSys::cartesian) {
   // no EXT_DIR BCs
   set_tensor_bc(bc,phys_bc,2,2);
  } else if ((CoordSys::CoordType) coord == CoordSys::CYLINDRICAL) {
   // no EXT_DIR BCs
   set_tensor_bc(bc,phys_bc,2,2);
  } else
   amrex::Error("(CoordSys::CoordType) coord invalid");
  } else
   amrex::Error("sdim invalid");
 
  std::string T33_strE="T33main"+postfix; 
  desc_lst.setComponent(indx,ibase_tensor_local,
    T33_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

#if (AMREX_SPACEDIM == 3)
  ibase_tensor_local++;

   // no EXT_DIR BCs
  set_tensor_bc(bc,phys_bc,0,2);
  std::string T13_strE="T13main"+postfix; 
  desc_lst.setComponent(indx,ibase_tensor_local,
    T13_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);
     
  ibase_tensor_local++;
     
   // no EXT_DIR BCs
  set_tensor_bc(bc,phys_bc,1,2);
  std::string T23_strE="T23main"+postfix; 
  desc_lst.setComponent(indx,ibase_tensor_local,
    T23_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);
#endif

  if (ibase_tensor_local==NUM_TENSOR_TYPE-1) {
   // do nothing
  } else {
   std::cout << "ibase_tensor_local=" << ibase_tensor_local << '\n';
   amrex::Error("ibase_tensor_local invalid");
  }

} // end subroutine set_tensor_extrap_components_main



void
NavierStokes::variableSetUp ()
{

     // AmrLevel.H, protected:
     // static DescriptorList desc_lst
     // static DescriptorList desc_lstGHOST
    BL_ASSERT(desc_lst.size() == 0);
    BL_ASSERT(desc_lstGHOST.size() == 0);

     // static variable
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        phys_bc.setLo(dir,NoSlipWall);
        phys_bc.setHi(dir,NoSlipWall);
    }

    std::cout << "entering variableSetUp before read_params\n";

    read_params();

    std::cout << "entering variableSetUp after read_params\n";
    std::cout << "num_materials= " << num_materials << '\n';
    std::cout << "num_state_material= " << num_state_material << '\n';
    std::cout << "num_state_base= " << num_state_base << '\n';
    std::cout << "num_species_var= " << num_species_var << '\n';
    std::cout << "num_materials_viscoelastic= " << 
     num_materials_viscoelastic << '\n';
    std::cout << "prescribe_temperature_outflow= " << 
     prescribe_temperature_outflow << '\n';
    
    int nmat=num_materials;
    int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

    int null_state_holds_data=0;
    int state_holds_data=1;

    if ((nmat<1)||(nmat>999)) {
     std::cout << "nmat= " << nmat << '\n';
     amrex::Error("nmat invalid in ns setup variable setup");
    }

    int nc=(AMREX_SPACEDIM+1)+
     nmat*(ngeom_raw+num_state_material)+1;

    int ncomp_per_burning=AMREX_SPACEDIM;
    int ncomp_per_tsat=2; // interface temperature and mass fraction
      // first nten components represent a status.
    int nburning=nten*(ncomp_per_burning+1);
    int ncomp_tsat=nten*(ncomp_per_tsat+1);

     // 11,12,22,33,13,23,XD,YD,ZD
    ncghost_elastic=NUM_TENSOR_TYPE+AMREX_SPACEDIM;
    ncghost_state=1+AMREX_SPACEDIM+nmat*ngeom_recon+1+nburning+ncomp_tsat+
         ncghost_elastic;

    std::string CC_postfix_str="CC";

    BCRec bc;

    ParmParse ppgeom("geometry");
    int coord;
    ppgeom.get("coord_sys",coord);
    
    if ((CoordSys::CoordType) coord == CoordSys::cartesian) {
     //rz_flag=0 
    } else if ((CoordSys::CoordType) coord == CoordSys::RZ) {
     //rz_flag=1
     if (AMREX_SPACEDIM!=2)
      amrex::Error("RZ only in 2D");
    } else if ((CoordSys::CoordType) coord == CoordSys::CYLINDRICAL) {
     //rz_flag=3
    } else
     amrex::Error("coord_sys invalid");

    BCRec phys_bc_pres;
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     phys_bc_pres.setLo(dir,phys_bc.lo(dir));
     phys_bc_pres.setHi(dir,phys_bc.hi(dir));
     if (phys_bc.lo(dir)==SlipWall)
      amrex::Error("SlipWall not allowed; use NoSlipWall instead");
     if (phys_bc.hi(dir)==SlipWall)
      amrex::Error("SlipWall not allowed; use NoSlipWall instead");
    }

// Umac_Type  -------------------------------------------

    umac_interp.interp_enable_spectral=enable_spectral;

     // AmrLevel.H, protected: static DescriptorList desc_lst
     // ngrow=0
    desc_lst.addDescriptor(Umac_Type,IndexType::TheUMACType(),
       0,1,&umac_interp,state_holds_data);
    desc_lstGHOST.addDescriptor(Umac_Type,IndexType::TheUMACType(),
       0,1,&umac_interp,null_state_holds_data);
    set_x_vel_bc(bc,phys_bc);

     // if new parameters added to the FILL routine, then
     // the following files must be modified:
     //  1. PROB_F.H
     //  2. PROB_3D.F90
     //  3. StateDescriptor.H  (BndryFuncDefaultSUSSMAN, 
     //      virtual void operator())
     //  4. StateDescriptor.cpp (operator(), m_func(...))
     //  5. StateData.cpp (calls to bndryFill)
    std::string u_mac_str="umac"; 
    desc_lst.setComponent(Umac_Type,0,u_mac_str,bc,FORT_UMACFILL,
      &umac_interp);

// Vmac_Type  -------------------------------------------

     // ngrow=0
    desc_lst.addDescriptor(Vmac_Type,IndexType::TheVMACType(),
      0,1,&umac_interp,state_holds_data);
    desc_lstGHOST.addDescriptor(Vmac_Type,IndexType::TheVMACType(),
      0,1,&umac_interp,null_state_holds_data);
    set_y_vel_bc(bc,phys_bc);

    std::string v_mac_str="vmac"; 
    desc_lst.setComponent(Vmac_Type,0,v_mac_str,bc,FORT_UMACFILL,
      &umac_interp);

// Wmac_Type  -------------------------------------------

    set_z_vel_bc(bc,phys_bc); // prevent warnings.

#if (AMREX_SPACEDIM == 3)

      // ngrow=0
    desc_lst.addDescriptor(Wmac_Type,IndexType::TheWMACType(),
      0,1,&umac_interp,state_holds_data);
    desc_lstGHOST.addDescriptor(Wmac_Type,IndexType::TheWMACType(),
      0,1,&umac_interp,null_state_holds_data);
    set_z_vel_bc(bc,phys_bc);

    std::string w_mac_str="wmac";
    desc_lst.setComponent(Wmac_Type,0,w_mac_str,bc,FORT_UMACFILL,
       &umac_interp);
#endif

    sem_interp_DEFAULT.interp_enable_spectral=enable_spectral;

    if ((enable_spectral==1)||(enable_spectral==2)) {
     sem_interp_HIGH_PARM.interp_enable_spectral=enable_spectral;
     sem_interp_LOW_PARM.interp_enable_spectral=0;
    } else if ((enable_spectral==0)||(enable_spectral==3)) {
     sem_interp_LOW_PARM.interp_enable_spectral=enable_spectral;
     sem_interp_HIGH_PARM.interp_enable_spectral=1;
    } else
     amrex::Error("enable_spectral invalid");

// DIV -------------------------------------------

    desc_lst.addDescriptor(DIV_Type,IndexType::TheCellType(),
     1,1,&sem_interp_DEFAULT,state_holds_data);

    desc_lstGHOST.addDescriptor(DIV_Type,IndexType::TheCellType(),
     1,1,&sem_interp_DEFAULT,null_state_holds_data);

    set_extrap_bc(bc,phys_bc);
    std::string divghost_str="divghost"; 
    desc_lstGHOST.setComponent(DIV_Type,0,
      divghost_str,bc,FORT_EXTRAPFILL,&sem_interp_DEFAULT);

    set_pressure_bc(bc,phys_bc_pres);
    std::string div_pres_str="div_pressure";
    desc_lst.setComponent(DIV_Type,0,
      div_pres_str,bc,FORT_PRESSUREFILL,&sem_interp_DEFAULT);


// Solid_State_Type  -------------------------------------------

    int nparts=im_solid_map.size();

    if ((nparts>=1)&&(nparts<nmat)) {
 
     desc_lst.addDescriptor(Solid_State_Type,IndexType::TheCellType(),
      1,nparts*AMREX_SPACEDIM,&pc_interp,state_holds_data);

     int ncghost_solid=1+AMREX_SPACEDIM;

     desc_lstGHOST.addDescriptor(Solid_State_Type,IndexType::TheCellType(),
      1,ncghost_solid,&pc_interp,null_state_holds_data);

     int dcomp=0;
     set_extrap_bc(bc,phys_bc);
     std::string extrap_str_solid="extrap_solid"; 
     desc_lstGHOST.setComponent(Solid_State_Type,dcomp,
      extrap_str_solid,bc,FORT_EXTRAPFILL,&pc_interp);

     dcomp++;
     std::string u_extrap_str_solid="u_extrap_solid";
     set_x_vel_extrap_bc(bc,phys_bc);
     desc_lstGHOST.setComponent(Solid_State_Type,dcomp,
      u_extrap_str_solid,bc,FORT_EXTRAPFILL,&pc_interp);

     dcomp++;
     std::string v_extrap_str_solid="v_extrap_solid";
     set_y_vel_extrap_bc(bc,phys_bc);
     desc_lstGHOST.setComponent(Solid_State_Type,dcomp,
      v_extrap_str_solid,bc,FORT_EXTRAPFILL,&pc_interp);

     if (AMREX_SPACEDIM==3) {
      dcomp++;
      std::string w_extrap_str_solid="w_extrap_solid";
      set_z_vel_extrap_bc(bc,phys_bc);
      desc_lstGHOST.setComponent(Solid_State_Type,dcomp,
       w_extrap_str_solid,bc,FORT_EXTRAPFILL,&pc_interp);
     }
     if (dcomp!=ncghost_solid-1)
      amrex::Error("dcomp invalid");

     for (int partid=0;partid<nparts;partid++) {

      int im_part=im_solid_map[partid];
      if ((im_part<0)||(im_part>=nmat))
       amrex::Error("im_part invalid");

      std::stringstream im_string_stream(std::stringstream::in |
       std::stringstream::out);
      im_string_stream << im_part+1;
      std::string im_string=im_string_stream.str();

      Vector<std::string> MOFvelocity_names_solid;
      MOFvelocity_names_solid.resize(AMREX_SPACEDIM);

      Vector<BCRec> MOFvelocity_bcs_solid;
      MOFvelocity_bcs_solid.resize(AMREX_SPACEDIM);

      int ibase_solid=0;

      std::string xvel_str_solid="x_velocity_solid"; 
      xvel_str_solid+=im_string;
      MOFvelocity_names_solid[ibase_solid]=xvel_str_solid;
      set_x_vel_bc(MOFvelocity_bcs_solid[ibase_solid],phys_bc);

      ibase_solid++;
     
      std::string yvel_str_solid="y_velocity_solid"; 
      yvel_str_solid+=im_string;
      MOFvelocity_names_solid[ibase_solid]=yvel_str_solid;
      set_y_vel_bc(MOFvelocity_bcs_solid[ibase_solid],phys_bc);
     
#if (AMREX_SPACEDIM == 3)
      ibase_solid++;

      std::string zvel_str_solid="z_velocity_solid"; 
      zvel_str_solid+=im_string;
      MOFvelocity_names_solid[ibase_solid]=zvel_str_solid;
      set_z_vel_bc(MOFvelocity_bcs_solid[ibase_solid],phys_bc);
#endif

      StateDescriptor::BndryFunc MOFvelocity_fill_class_solid(FORT_SOLVFILL,
       FORT_GROUP_SOLVFILL);

      desc_lst.setComponent(Solid_State_Type,
       partid*AMREX_SPACEDIM,
       MOFvelocity_names_solid,
       MOFvelocity_bcs_solid,
       MOFvelocity_fill_class_solid,
       &sem_interp_DEFAULT);
     } // partid=0..nparts-1

    } else if (nparts==0) {
     // do nothing
    } else
     amrex::Error("nparts invalid");

// XDmac_Type ... ZDmac_Type

    int displacement_enable_spectral=0;

    xd_mac_interp.interp_enable_spectral=displacement_enable_spectral;
    xd_mac_lo_interp.interp_enable_spectral=0;

     // ngrow=0
    desc_lst.addDescriptor(XDmac_Type,IndexType::TheUMACType(),
     0,1,&xd_mac_interp,state_holds_data);
    // ngrow=0,ncomp=1
    desc_lstGHOST.addDescriptor(XDmac_Type,IndexType::TheUMACType(),
     0,1,&xd_mac_lo_interp,null_state_holds_data);

     // same as x_vel_bc except that EXT_DIR => FOEXTRAP
    set_x_vel_extrap_bc(bc,phys_bc);
    std::string extrap_str="xd_extrap"; 
      //dcomp=0
    desc_lstGHOST.setComponent(XDmac_Type,0,
      extrap_str,bc,FORT_X_EXTRAPFILL,&xd_mac_lo_interp);

     // ngrow=0
    desc_lst.addDescriptor(YDmac_Type,IndexType::TheVMACType(),
     0,1,&xd_mac_interp,state_holds_data);
     // ngrow=0,ncomp=1
    desc_lstGHOST.addDescriptor(YDmac_Type,IndexType::TheVMACType(),
     0,1,&xd_mac_lo_interp,null_state_holds_data);

     // same as y_vel_bc except that EXT_DIR => FOEXTRAP
    set_y_vel_extrap_bc(bc,phys_bc);
    extrap_str="yd_extrap"; 
      //dcomp=0
    desc_lstGHOST.setComponent(YDmac_Type,0,
      extrap_str,bc,FORT_X_EXTRAPFILL,&xd_mac_lo_interp);

#if (AMREX_SPACEDIM == 3)
     // ngrow=0
    desc_lst.addDescriptor(ZDmac_Type,IndexType::TheWMACType(),
     0,1,&xd_mac_interp,state_holds_data);
     // ngrow=0,ncomp=1
    desc_lstGHOST.addDescriptor(ZDmac_Type,IndexType::TheWMACType(),
     0,1,&xd_mac_lo_interp,null_state_holds_data);

     // same as z_vel_bc except that EXT_DIR => FOEXTRAP
    set_z_vel_extrap_bc(bc,phys_bc);
    extrap_str="zd_extrap"; 
      //dcomp=0
    desc_lstGHOST.setComponent(ZDmac_Type,0,
      extrap_str,bc,FORT_X_EXTRAPFILL,&xd_mac_lo_interp);
#endif

    std::string xd_mac_name="XDMAC";
    BCRec xd_mac_bcs;
     // same as x_vel_bc except that EXT_DIR => FOEXTRAP
    set_x_vel_extrap_bc(xd_mac_bcs,phys_bc);
    desc_lst.setComponent(XDmac_Type,0,xd_mac_name,xd_mac_bcs,
       FORT_XDMACFILL,&xd_mac_interp);

    std::string yd_mac_name="YDMAC";
    BCRec yd_mac_bcs;
     // same as y_vel_bc except that EXT_DIR => FOEXTRAP
    set_y_vel_extrap_bc(yd_mac_bcs,phys_bc);
    desc_lst.setComponent(YDmac_Type,0,yd_mac_name,yd_mac_bcs,
       FORT_XDMACFILL,&xd_mac_interp);

#if (AMREX_SPACEDIM == 3)
    std::string zd_mac_name="ZDMAC";
    BCRec zd_mac_bcs;
     // same as z_vel_bc except that EXT_DIR => FOEXTRAP
    set_z_vel_extrap_bc(zd_mac_bcs,phys_bc);
    desc_lst.setComponent(ZDmac_Type,0,zd_mac_name,zd_mac_bcs,
       FORT_XDMACFILL,&xd_mac_interp);
#endif

// Tensor_Type  -------------------------------------------

    if (num_materials_viscoelastic!=im_elastic_map.size())
     amrex::Error("num_materials_viscoelastic!=im_elastic_map.size()");

    if (NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
     // do nothing
    } else
     amrex::Error("NUM_TENSOR_TYPE invalid");

    if (NUM_CELL_ELASTIC==num_materials_viscoelastic*NUM_TENSOR_TYPE) {
     // do nothing
    } else
     amrex::Error("NUM_CELL_ELASTIC invalid");

    if ((num_materials_viscoelastic>=1)&&
        (num_materials_viscoelastic<=nmat)) {

     desc_lst.addDescriptor(Tensor_Type,IndexType::TheCellType(),
      1,NUM_CELL_ELASTIC,
      &pc_interp,
      state_holds_data);

      // ngrow=1
     desc_lstGHOST.addDescriptor(Tensor_Type,IndexType::TheCellType(),
      1,ncghost_elastic,&pc_interp,null_state_holds_data);

      // setComponent: 0..NUM_TENSOR_TYPE-1
     set_tensor_extrap_components(coord,CC_postfix_str,Tensor_Type,0);

     int ibase_tensor=NUM_TENSOR_TYPE;

      // same as x_vel_bc except that EXT_DIR => FOEXTRAP
     set_x_vel_extrap_bc(bc,phys_bc);
     std::string xdisplace_strE="XDISPLACEextrap"; 
     desc_lstGHOST.setComponent(Tensor_Type,ibase_tensor,
      xdisplace_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

     ibase_tensor++;

      // same as y_vel_bc except that EXT_DIR => FOEXTRAP
     set_y_vel_extrap_bc(bc,phys_bc);
     std::string ydisplace_strE="YDISPLACEextrap"; 
     desc_lstGHOST.setComponent(Tensor_Type,ibase_tensor,
      ydisplace_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

#if (AMREX_SPACEDIM == 3)
     if (AMREX_SPACEDIM==3) {
      ibase_tensor++;

      // same as z_vel_bc except that EXT_DIR => FOEXTRAP
      set_z_vel_extrap_bc(bc,phys_bc);
      std::string zdisplace_strE="ZDISPLACEextrap"; 
      desc_lstGHOST.setComponent(Tensor_Type,ibase_tensor,
       zdisplace_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);
     }
#endif

     if (ibase_tensor!=ncghost_elastic-1)
      amrex::Error("ibase_tensor invalid");

     for (int partid=0;partid<num_materials_viscoelastic;partid++) {

      int im_part=im_elastic_map[partid];
      if ((im_part<0)||(im_part>=nmat))
       amrex::Error("im_part invalid");

      std::stringstream im_string_stream(std::stringstream::in |
       std::stringstream::out);
      im_string_stream << im_part+1;
      std::string im_string=im_string_stream.str();

      Vector<std::string> MOFvelocity_names_tensor;
      MOFvelocity_names_tensor.resize(NUM_TENSOR_TYPE);

      Vector<BCRec> MOFvelocity_bcs_tensor;
      MOFvelocity_bcs_tensor.resize(NUM_TENSOR_TYPE);

      ibase_tensor=0;

       // analogous to u_x
       // reflect at symmetric BC x-faces
       // negative reflect at symmetric BC y or z-faces
       // reflect at outflow
       // reflect at walls.
      std::string T11_str="T11"; 
      T11_str+=im_string;
      MOFvelocity_names_tensor[ibase_tensor]=T11_str;
      set_tensor_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc,0,0);

      ibase_tensor++;
     
      std::string T12_str="T12"; 
      T12_str+=im_string; 
      MOFvelocity_names_tensor[ibase_tensor]=T12_str;
      set_tensor_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc,0,1);
    
      ibase_tensor++;
     
      std::string T22_str="T22"; 
      T22_str+=im_string; 
      MOFvelocity_names_tensor[ibase_tensor]=T22_str;
      set_tensor_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc,1,1);

      ibase_tensor++;
     
      std::string T33_str="T33"; 
      T33_str+=im_string; 
      MOFvelocity_names_tensor[ibase_tensor]=T33_str;

      if (AMREX_SPACEDIM==2) {
       if ((CoordSys::CoordType) coord == CoordSys::RZ) {
        set_hoop_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc);
       } else if ((CoordSys::CoordType) coord == CoordSys::cartesian) {
	   // placeholder: Q33 should always be 0
        set_hoop_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc);
       } else if ((CoordSys::CoordType) coord == CoordSys::CYLINDRICAL) {
	   // placeholder: Q33 should always be 0
        set_hoop_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc);
       } else
        amrex::Error("(CoordSys::CoordType) coord invalid");
      } else if (AMREX_SPACEDIM==3) {
       if ((CoordSys::CoordType) coord == CoordSys::cartesian) {
        set_tensor_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc,2,2);
       } else if ((CoordSys::CoordType) coord == CoordSys::CYLINDRICAL) {
        set_tensor_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc,2,2);
       } else
        amrex::Error("(CoordSys::CoordType) coord invalid");
      } else
       amrex::Error("sdim invalid");

#if (AMREX_SPACEDIM == 3)
      ibase_tensor++;
     
      std::string T13_str="T13"; 
      T13_str+=im_string; 
      MOFvelocity_names_tensor[ibase_tensor]=T13_str;
      set_tensor_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc,0,2);

      ibase_tensor++;
     
      std::string T23_str="T23"; 
      T23_str+=im_string; 
      MOFvelocity_names_tensor[ibase_tensor]=T23_str;
      set_tensor_bc(MOFvelocity_bcs_tensor[ibase_tensor],phys_bc,1,2);

#endif

      StateDescriptor::BndryFunc MOFvelocity_fill_class_tensor(
       FORT_TENSORFILL,
       FORT_GROUP_TENSORFILL);

      desc_lst.setComponent(Tensor_Type,
       partid*NUM_TENSOR_TYPE,
       MOFvelocity_names_tensor,
       MOFvelocity_bcs_tensor,
       MOFvelocity_fill_class_tensor,
       &pc_interp);

     } // partid=0..nparts-1

     //ngrow=1
     desc_lst.addDescriptor(TensorXU_Type,IndexType::TheCellType(),
      1,NUM_TENSOR_TYPE,&tensor_pc_interp,
      null_state_holds_data);
     desc_lstGHOST.addDescriptor(TensorXU_Type,IndexType::TheCellType(),
      1,NUM_TENSOR_TYPE,&tensor_pc_interp,
      null_state_holds_data);

     std::string MAC_postfix_str="XU";
     set_tensor_extrap_components(coord,MAC_postfix_str,TensorXU_Type,0);
     set_tensor_extrap_components_main(coord,MAC_postfix_str,TensorXU_Type);

     //ngrow=1
     desc_lst.addDescriptor(TensorYU_Type,IndexType::TheYUMACType(),
      1,NUM_TENSOR_TYPE,&tensor_pc_interp,
      null_state_holds_data);
     desc_lstGHOST.addDescriptor(TensorYU_Type,IndexType::TheYUMACType(),
      1,NUM_TENSOR_TYPE,&tensor_pc_interp,
      null_state_holds_data);

     MAC_postfix_str="YU";
     set_tensor_extrap_components(coord,MAC_postfix_str,TensorYU_Type,0);
     set_tensor_extrap_components_main(coord,MAC_postfix_str,TensorYU_Type);

      //ngrow=1
     desc_lst.addDescriptor(TensorZU_Type,IndexType::TheZUMACType(),
      1,NUM_TENSOR_TYPE,&tensor_pc_interp,
      null_state_holds_data);
     desc_lstGHOST.addDescriptor(TensorZU_Type,IndexType::TheZUMACType(),
      1,NUM_TENSOR_TYPE,&tensor_pc_interp,
      null_state_holds_data);

     MAC_postfix_str="ZU";
     set_tensor_extrap_components(coord,MAC_postfix_str,TensorZU_Type,0);
     set_tensor_extrap_components_main(coord,MAC_postfix_str,TensorZU_Type);

     //ngrow=1
     desc_lst.addDescriptor(TensorZV_Type,IndexType::TheZVMACType(),
      1,NUM_TENSOR_TYPE,&tensor_pc_interp,
      null_state_holds_data);
     desc_lstGHOST.addDescriptor(TensorZV_Type,IndexType::TheZVMACType(),
      1,NUM_TENSOR_TYPE,&tensor_pc_interp,
      null_state_holds_data);

     MAC_postfix_str="ZV";
     set_tensor_extrap_components(coord,MAC_postfix_str,TensorZV_Type,0);
     set_tensor_extrap_components_main(coord,MAC_postfix_str,TensorZV_Type);

    } else if (num_materials_viscoelastic==0) {
     // do nothing
    } else
     amrex::Error("num_materials_viscoelastic invalid");


// LEVELSET ------------------------------------------------- 

    int ncomp_ls=(AMREX_SPACEDIM+1)*nmat;

    desc_lst.addDescriptor(LS_Type,IndexType::TheCellType(),
     1,ncomp_ls,&pc_interp,state_holds_data);

     // components 0..nmat * AMREX_SPACEDIM-1 are for interface normal vectors.
     // components nmat * AMREX_SPACEDIM .. nmat * AMREX_SPACEDIM + 
     //   nmat * (AMREX_SPACEDIM+1) are the same (except for the string name)
     //   as for dest_lst.
    int ncomp_LS_ghost=(2*AMREX_SPACEDIM+1)*nmat;

    desc_lstGHOST.addDescriptor(LS_Type,IndexType::TheCellType(),
     1,ncomp_LS_ghost,&pc_interp,null_state_holds_data);

    int dcomp=0;
    for (int imls=0;imls<nmat;imls++) { 

     std::stringstream im_string_stream(std::stringstream::in |
      std::stringstream::out);
     im_string_stream << imls+1;
     std::string im_string=im_string_stream.str();

     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

      std::string nrm_extrap_str=" ";
      if (dir==0) {
       set_x_vel_extrap_bc(bc,phys_bc);
       nrm_extrap_str="x_norm"; 
       nrm_extrap_str+=im_string;
      } else if (dir==1) {
       set_y_vel_extrap_bc(bc,phys_bc);
       nrm_extrap_str="y_norm"; 
       nrm_extrap_str+=im_string;
      } else if ((dir==2)&&(AMREX_SPACEDIM==3)) {
       set_z_vel_extrap_bc(bc,phys_bc);
       nrm_extrap_str="z_norm"; 
       nrm_extrap_str+=im_string;
      } else 
       amrex::Error("dir invalid ns_setup");

      desc_lstGHOST.setComponent(LS_Type,dcomp,
        nrm_extrap_str,bc,FORT_EXTRAPFILL,&pc_interp);
 
      dcomp++;
     } // dir=0..AMREX_SPACEDIM-1

    } // imls=0..nmat-1

    if (dcomp!=AMREX_SPACEDIM*nmat)
     amrex::Error("dcomp invalid");

    Vector<std::string> LS_HO_names;
    LS_HO_names.resize(ncomp_ls);
    Vector<BCRec> LS_HO_bcs;
    LS_HO_bcs.resize(ncomp_ls);

    dcomp=0;
  
    for (int imls=0;imls<nmat;imls++) {

     std::stringstream im_string_stream(std::stringstream::in |
      std::stringstream::out);
     im_string_stream << imls+1;
     std::string im_string=im_string_stream.str();

     std::string LS_str="LS_HO"; 
     LS_str+=im_string;
     LS_HO_names[imls]=LS_str;
     set_scalar_vof_bc(LS_HO_bcs[imls],phys_bc);

     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

      std::string nrm_extrap_str=" ";
      if (dir==0) {
       set_x_vel_extrap_bc(LS_HO_bcs[nmat+dcomp],phys_bc);
       nrm_extrap_str="x_norm_HO"; 
       nrm_extrap_str+=im_string;
      } else if (dir==1) {
       set_y_vel_extrap_bc(LS_HO_bcs[nmat+dcomp],phys_bc);
       nrm_extrap_str="y_norm_HO"; 
       nrm_extrap_str+=im_string;
      } else if ((dir==2)&&(AMREX_SPACEDIM==3)) {
       set_z_vel_extrap_bc(LS_HO_bcs[nmat+dcomp],phys_bc);
       nrm_extrap_str="z_norm_HO"; 
       nrm_extrap_str+=im_string;
      } else 
       amrex::Error("dir invalid ns_setup");

      LS_HO_names[nmat+dcomp]=nrm_extrap_str;

      dcomp++;

     } // dir=0..AMREX_SPACEDIM-1

    }  // imls=0...nmat-1

    if (dcomp!=AMREX_SPACEDIM*nmat)
     amrex::Error("dcomp invalid");
    if (dcomp+nmat!=ncomp_ls)
     amrex::Error("dcomp invalid");

     // GROUP_LS_HO_FILL: grouplsBC for components 1..nmat
     //                   extrapBC for components nmat+1..nmat * (sdim+1)
    StateDescriptor::BndryFunc LS_HO_fill_class(FORT_LS_HO_FILL,
       FORT_GROUP_LS_HO_FILL);

    ls_ho_interp_HIGH_PARM.LSHOInterp_nmat=nmat;
    ls_ho_interp_HIGH_PARM.LSHOInterp_LO=0; // 0=use normals 1=piecewise const 

    ls_ho_interp_LOW_PARM.LSHOInterp_nmat=nmat;
    ls_ho_interp_LOW_PARM.LSHOInterp_LO=1; // 0=use normals 1=piecewise const 

    desc_lstGHOST.setComponent(LS_Type,
      AMREX_SPACEDIM*nmat,LS_HO_names,
      LS_HO_bcs,LS_HO_fill_class,&ls_ho_interp_HIGH_PARM);

    Vector<std::string> LS_main_names;
    LS_main_names.resize(ncomp_ls);
    Vector<BCRec> LS_main_bcs;
    LS_main_bcs.resize(ncomp_ls);

    dcomp=0;

    for (int imls=0;imls<nmat;imls++) {

     std::stringstream im_string_stream(std::stringstream::in |
      std::stringstream::out);
     im_string_stream << imls+1;
     std::string im_string=im_string_stream.str();

     std::string LS_str="LS"; 
     LS_str+=im_string;
     LS_main_names[imls]=LS_str;
     set_scalar_vof_bc(LS_main_bcs[imls],phys_bc);

     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

      std::string nrm_extrap_str=" ";
      if (dir==0) {
       set_x_vel_extrap_bc(LS_main_bcs[nmat+dcomp],phys_bc);
       nrm_extrap_str="x_norm_main"; 
       nrm_extrap_str+=im_string;
      } else if (dir==1) {
       set_y_vel_extrap_bc(LS_main_bcs[nmat+dcomp],phys_bc);
       nrm_extrap_str="y_norm_main"; 
       nrm_extrap_str+=im_string;
      } else if ((dir==2)&&(AMREX_SPACEDIM==3)) {
       set_z_vel_extrap_bc(LS_main_bcs[nmat+dcomp],phys_bc);
       nrm_extrap_str="z_norm_main"; 
       nrm_extrap_str+=im_string;
      } else 
       amrex::Error("dir invalid ns_setup");

      LS_main_names[nmat+dcomp]=nrm_extrap_str;

      dcomp++;
     } // dir=0..sdim-1

    }  // imls=0...nmat-1

    if (dcomp!=AMREX_SPACEDIM*nmat)
     amrex::Error("dcomp invalid");
    if (dcomp+nmat!=ncomp_ls)
     amrex::Error("dcomp invalid");

     // GROUP_LS_HO_FILL: grouplsBC for components 1..nmat
     //                   extrapBC for components nmat+1..nmat * (sdim+1)
    StateDescriptor::BndryFunc LS_main_fill_class(FORT_LS_HO_FILL,
       FORT_GROUP_LS_HO_FILL);

    ls_ho_interp_HIGH_PARM.LSHOInterp_nmat=nmat;
    ls_ho_interp_HIGH_PARM.LSHOInterp_LO=0; // 0=use normals 1=piecewise const 

    ls_ho_interp_LOW_PARM.LSHOInterp_nmat=nmat;
    ls_ho_interp_LOW_PARM.LSHOInterp_LO=1; // 0=use normals 1=piecewise const 

      //default is low order (piecewise constant)
    desc_lst.setComponent(LS_Type,0,LS_main_names,
      LS_main_bcs,LS_main_fill_class,&ls_ho_interp_LOW_PARM);


// State_Type  ------------------------------------------------- 
// newdata FABS have ncomp=desc->nComp() components.
//
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
     1,nc,&pc_interp,state_holds_data);

    desc_lstGHOST.addDescriptor(State_Type,IndexType::TheCellType(),
     1,ncghost_state,&pc_interp,null_state_holds_data);

    dcomp=0;
    set_extrap_bc(bc,phys_bc);
    std::string state_extrap_str="extrap"; 
    desc_lstGHOST.setComponent(State_Type,0,
      state_extrap_str,bc,FORT_EXTRAPFILL,&pc_interp);

    dcomp++;
    std::string u_extrap_str="u_extrap";
    set_x_vel_extrap_bc(bc,phys_bc);
    desc_lstGHOST.setComponent(State_Type,dcomp,
     u_extrap_str,bc,FORT_EXTRAPFILL,&pc_interp);

    dcomp++;
    std::string v_extrap_str="v_extrap";
    set_y_vel_extrap_bc(bc,phys_bc);
    desc_lstGHOST.setComponent(State_Type,dcomp,
     v_extrap_str,bc,FORT_EXTRAPFILL,&pc_interp);

    if (AMREX_SPACEDIM==3) {
     dcomp++;
     std::string w_extrap_str="w_extrap";
     set_z_vel_extrap_bc(bc,phys_bc);
     desc_lstGHOST.setComponent(State_Type,dcomp,
      w_extrap_str,bc,FORT_EXTRAPFILL,&pc_interp);
    }


    // in NavierStokes::VOF_Recon
    // 1. get MOF data with 1 ghost cell (so that CMOF can be chosen)
    // 2. reconstruct interior cells only.
    // 3. do extended filpatch; MOF used for coarse/fine and ext_dir cells.
    Vector<std::string> EXTMOF_names;
    EXTMOF_names.resize(nmat*ngeom_recon);
    Vector<BCRec> EXTMOF_bcs;
    EXTMOF_bcs.resize(nmat*ngeom_recon);

    for (int im=0;im<nmat;im++) {

     int ibase_extmof=im*ngeom_recon;

     std::stringstream im_string_stream(std::stringstream::in |
        std::stringstream::out);

     im_string_stream << im+1;
     std::string im_string=im_string_stream.str();

     std::string vof_str="vofE"; 
     vof_str+=im_string; 
     EXTMOF_names[ibase_extmof]=vof_str;
     set_scalar_vof_bc(EXTMOF_bcs[ibase_extmof],phys_bc);

     ibase_extmof++;
     std::string cenx_str="cenxE"; 
     cenx_str+=im_string; 
     EXTMOF_names[ibase_extmof]=cenx_str;
     set_x_vel_bc(EXTMOF_bcs[ibase_extmof],phys_bc);

     ibase_extmof++;
     std::string ceny_str="cenyE"; 
     ceny_str+=im_string; 
     EXTMOF_names[ibase_extmof]=ceny_str;
     set_y_vel_bc(EXTMOF_bcs[ibase_extmof],phys_bc);

#if (AMREX_SPACEDIM==3)
     ibase_extmof++;
     std::string cenz_str="cenzE"; 
     cenz_str+=im_string; 
     EXTMOF_names[ibase_extmof]=cenz_str;
     set_z_vel_bc(EXTMOF_bcs[ibase_extmof],phys_bc);
#endif    

     ibase_extmof++;
     std::string order_str="orderE"; 
     order_str+=im_string; 
     EXTMOF_names[ibase_extmof]=order_str;
     set_scalar_vof_bc(EXTMOF_bcs[ibase_extmof],phys_bc);

     ibase_extmof++;
     std::string nrmx_str="nrmxE"; 
     nrmx_str+=im_string; 
     EXTMOF_names[ibase_extmof]=nrmx_str;
     set_x_vel_bc(EXTMOF_bcs[ibase_extmof],phys_bc);

     ibase_extmof++;
     std::string nrmy_str="nrmyE"; 
     nrmy_str+=im_string; 
     EXTMOF_names[ibase_extmof]=nrmy_str;
     set_y_vel_bc(EXTMOF_bcs[ibase_extmof],phys_bc);

#if (AMREX_SPACEDIM==3)
     ibase_extmof++;
     std::string nrmz_str="nrmzE"; 
     nrmz_str+=im_string; 
     EXTMOF_names[ibase_extmof]=nrmz_str;
     set_z_vel_bc(EXTMOF_bcs[ibase_extmof],phys_bc);
#endif    

     ibase_extmof++;
     std::string intercept_str="interceptE"; 
     intercept_str+=im_string; 
     EXTMOF_names[ibase_extmof]=intercept_str;
     set_scalar_vof_bc(EXTMOF_bcs[ibase_extmof],phys_bc);

     if (ibase_extmof!=(im+1)*ngeom_recon-1)
      amrex::Error("ibase_extmof invalid");

    }  // im=0..nmat-1  (vfrac, cen, order, slope,int)

    StateDescriptor::BndryFunc EXTMOF_fill_class(FORT_EXTMOFFILL,
       FORT_GROUP_EXTMOFFILL);

    int extrecon_start_pos=1+AMREX_SPACEDIM; //state extrap, vel extrap

    multi_extmof_interp.multiMOFInterp_nmat=nmat;
    multi_extmof_interp.multiMOFInterp_ngeom_raw=ngeom_raw;
    multi_extmof_interp.multiMOFInterp_ngeom_recon=ngeom_recon;

    desc_lstGHOST.setComponent(State_Type,extrecon_start_pos,EXTMOF_names,
     EXTMOF_bcs,EXTMOF_fill_class,&multi_extmof_interp);

    set_extrap_bc(bc,phys_bc);
    std::string maskextrap_str="maskSEMextrap"; 

    int mask_scomp=extrecon_start_pos+nmat*ngeom_recon;
    desc_lstGHOST.setComponent(State_Type,mask_scomp,
      maskextrap_str,bc,FORT_EXTRAPFILL,&mask_sem_interp);

    Vector<std::string> BURNVEL_names;
    BURNVEL_names.resize(nburning);
    Vector<BCRec> BURNVEL_bcs;
    BURNVEL_bcs.resize(nburning);

    for (int im=0;im<nten;im++) {

     std::stringstream im_string_stream(std::stringstream::in |
        std::stringstream::out);

     im_string_stream << im+1;
     std::string im_string=im_string_stream.str();

     std::string status_str="burnstat"; 
     status_str+=im_string; 
     BURNVEL_names[im]=status_str;
     set_extrap_bc(BURNVEL_bcs[im],phys_bc);

    }  // im=0..nten-1  (status for burning velocity)

    for (int im=0;im<nten;im++) {

     int ibase_burnvel=nten+im*ncomp_per_burning;

     std::stringstream im_string_stream(std::stringstream::in |
        std::stringstream::out);

     im_string_stream << im+1;
     std::string im_string=im_string_stream.str();

     std::string burnxvel_str="burnxvel"; 
     burnxvel_str+=im_string; 
     BURNVEL_names[ibase_burnvel]=burnxvel_str;
     set_x_vel_extrap_bc(BURNVEL_bcs[ibase_burnvel],phys_bc);

     ibase_burnvel++;
     std::string burnyvel_str="burnyvel"; 
     burnyvel_str+=im_string; 
     BURNVEL_names[ibase_burnvel]=burnyvel_str;
     set_y_vel_extrap_bc(BURNVEL_bcs[ibase_burnvel],phys_bc);

#if (AMREX_SPACEDIM==3)
     ibase_burnvel++;
     std::string burnzvel_str="burnzvel"; 
     burnzvel_str+=im_string; 
     BURNVEL_names[ibase_burnvel]=burnzvel_str;
     set_z_vel_extrap_bc(BURNVEL_bcs[ibase_burnvel],phys_bc);
#endif    

     if (ibase_burnvel!=nten+(im+1)*ncomp_per_burning-1)
      amrex::Error("ibase_burnvel invalid");

    }  // im=0..nten-1  (burning velocity)

    StateDescriptor::BndryFunc BURNVEL_fill_class(FORT_EXTRAPFILL,
       FORT_GROUP_EXTRAPFILL);

    int burnvel_start_pos=mask_scomp+1;

    burnvel_interp.burnvel_nmat=nmat;
    burnvel_interp.burnvel_nten=nten;
    burnvel_interp.burnvel_ncomp_per=ncomp_per_burning;

    desc_lstGHOST.setComponent(State_Type,burnvel_start_pos,BURNVEL_names,
     BURNVEL_bcs,BURNVEL_fill_class,&burnvel_interp);

    Vector<std::string> TSAT_names;
    TSAT_names.resize(ncomp_tsat);
    Vector<BCRec> TSAT_bcs;
    TSAT_bcs.resize(ncomp_tsat);

    for (int im=0;im<nten;im++) {

     std::stringstream im_string_stream(std::stringstream::in |
        std::stringstream::out);

     im_string_stream << im+1;
     std::string im_string=im_string_stream.str();

     std::string status_str="tsatstat"; 
     status_str+=im_string; 
     TSAT_names[im]=status_str;
     set_extrap_bc(TSAT_bcs[im],phys_bc);

    }  // im=0..nten-1  (status for TSAT)

    for (int im=0;im<nten;im++) {

     int ibase_tsat=nten+im*ncomp_per_tsat;

     std::stringstream im_string_stream(std::stringstream::in |
        std::stringstream::out);

     im_string_stream << im+1;
     std::string im_string=im_string_stream.str();

     std::string tsat_str="tsat"; 
     tsat_str+=im_string; 
     TSAT_names[ibase_tsat]=tsat_str;
     set_extrap_bc(TSAT_bcs[ibase_tsat],phys_bc);

     ibase_tsat++;

     std::string massfrac_str="massfracI"; 
     massfrac_str+=im_string; 
     TSAT_names[ibase_tsat]=massfrac_str;
     set_extrap_bc(TSAT_bcs[ibase_tsat],phys_bc);

     if (ibase_tsat!=nten+(im+1)*ncomp_per_tsat-1)
      amrex::Error("ibase_tsat invalid");

    }  // im=0..nten-1  (TSAT)

    StateDescriptor::BndryFunc TSAT_fill_class(FORT_EXTRAPFILL,
       FORT_GROUP_EXTRAPFILL);

    int tsat_start_pos=burnvel_start_pos+nburning;

    tsat_interp.burnvel_nmat=nmat;
    tsat_interp.burnvel_nten=nten;
    //interface temperature and mass fraction
    tsat_interp.burnvel_ncomp_per=ncomp_per_tsat; 

    desc_lstGHOST.setComponent(State_Type,tsat_start_pos,TSAT_names,
     TSAT_bcs,TSAT_fill_class,&tsat_interp);

    if (tsat_start_pos==ncghost_state-ncghost_elastic-ncomp_tsat) {
     // do nothing
    } else
     amrex::Error("tsat_start_pos invalid");

     // setComponent: 0..NUM_TENSOR_TYPE-1
    set_tensor_extrap_components(coord,CC_postfix_str,State_Type,
		    ncghost_state-ncghost_elastic);

    int ibase_state_tensor=ncghost_state-ncghost_elastic+NUM_TENSOR_TYPE;

     // same as x_vel_bc except that EXT_DIR => FOEXTRAP
    set_x_vel_extrap_bc(bc,phys_bc);
    std::string xdisplace_state_strE="XDISPLACEextrap"; 
    desc_lstGHOST.setComponent(State_Type,ibase_state_tensor,
      xdisplace_state_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

    ibase_state_tensor++;

     // same as y_vel_bc except that EXT_DIR => FOEXTRAP
    set_y_vel_extrap_bc(bc,phys_bc);
    std::string ydisplace_state_strE="YDISPLACEextrap"; 
    desc_lstGHOST.setComponent(State_Type,ibase_state_tensor,
      ydisplace_state_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);

#if (AMREX_SPACEDIM == 3)
    if (AMREX_SPACEDIM==3) {
     ibase_state_tensor++;

     // same as z_vel_bc except that EXT_DIR => FOEXTRAP
     set_z_vel_extrap_bc(bc,phys_bc);
     std::string zdisplace_state_strE="ZDISPLACEextrap"; 
     desc_lstGHOST.setComponent(State_Type,ibase_state_tensor,
       zdisplace_state_strE,bc,FORT_EXTRAPFILL,&tensor_pc_interp);
    }
#endif

    if (ibase_state_tensor!=ncghost_state-1)
     amrex::Error("ibase_state_tensor invalid");


    // boundary routines are of type BndryFuncDefaultSUSSMAN
    // setComponent expects a parameter of type 
    // StateDescriptor::BndryFunc&
    // Then it gives the following line:
    // bc_func.set(comp,func.clone());
    // where clone() is a member of BndryFunc:
    // StateDescriptor::BndryFunc*
    // StateDescriptor::BndryFunc::clone () const
    // { return new BndryFunc(*this); }
    // one of the constructors of BndryFunc is:
    // BndryFunc (BndryFuncDefaultSUSSMAN inFunc);
    // another constructor:
    // BndryFunc (BndryFuncDefaultSUSSMAN inFunc,BndryFuncDefaultSUSSMAN gFunc);


    Vector<std::string> MOFvelocity_names;
    MOFvelocity_names.resize(AMREX_SPACEDIM);

    Vector<BCRec> MOFvelocity_bcs;
    MOFvelocity_bcs.resize(AMREX_SPACEDIM);

    int ibase_state=0;

    std::string xvel_str="x_velocity"; 
    MOFvelocity_names[ibase_state]=xvel_str;
    set_x_vel_bc(MOFvelocity_bcs[ibase_state],phys_bc);

    ibase_state++;
     
    std::string yvel_str="y_velocity"; 
    MOFvelocity_names[ibase_state]=yvel_str;
    set_y_vel_bc(MOFvelocity_bcs[ibase_state],phys_bc);
     
#if (AMREX_SPACEDIM == 3)
    ibase_state++;

    std::string zvel_str="z_velocity"; 
    MOFvelocity_names[ibase_state]=zvel_str;
    set_z_vel_bc(MOFvelocity_bcs[ibase_state],phys_bc);
#endif

    StateDescriptor::BndryFunc MOFvelocity_fill_class(FORT_VELFILL,
       FORT_GROUP_VELFILL);

    desc_lst.setComponent(State_Type,0,
      MOFvelocity_names,
      MOFvelocity_bcs,MOFvelocity_fill_class,&sem_interp_DEFAULT);

     // pressure
    set_pressure_bc(bc,phys_bc_pres);
    std::string pres_str="pressure"; 
    desc_lst.setComponent(State_Type,AMREX_SPACEDIM,
      pres_str,bc,FORT_PRESSUREFILL,&sem_interp_DEFAULT);

    Vector<std::string> MOFstate_names;
    MOFstate_names.resize(num_state_material*nmat);

    Vector<BCRec> MOFstate_bcs;
    MOFstate_bcs.resize(num_state_material*nmat);

    for (int im=0;im<nmat;im++) {

     int ibase_transport=im*num_state_material;

     std::stringstream im_string_stream(std::stringstream::in |
      std::stringstream::out);

     im_string_stream << im+1;
     std::string im_string=im_string_stream.str();

     std::string density_str="density"; 
     density_str+=im_string; 
     MOFstate_names[ibase_transport]=density_str;
     set_scalar_bc(MOFstate_bcs[ibase_transport],phys_bc);

     ibase_transport++;

     std::string temperature_str="temperature"; 
     temperature_str+=im_string; 
     MOFstate_names[ibase_transport]=temperature_str;
     
      // 0=Dirichlet at inflow, insulating at walls and outflow (default)
      // 1=Dirichlet at inflow and outflow, insulating at walls.
      // 2=Dirichlet at inflow and at walls, insulating at outflow. 
      // 3=Dirichlet at inflow, outflow, and walls.
     if (prescribe_temperature_outflow==2)  
      set_custom_temperature_bc(MOFstate_bcs[ibase_transport],
		      temperature_phys_bc);
     else if (prescribe_temperature_outflow==1)
      set_temperature_bc(MOFstate_bcs[ibase_transport],
		      temperature_phys_bc);
     else if (prescribe_temperature_outflow==0)
      set_scalar_bc(MOFstate_bcs[ibase_transport],
		      temperature_phys_bc);
     else if (prescribe_temperature_outflow==3)
      set_custom2_temperature_bc(MOFstate_bcs[ibase_transport],
		      temperature_phys_bc);
     else 
      amrex::Error("prescribe_temperature_outflow invalid");

     for (int spec_comp=0;spec_comp<num_species_var;spec_comp++) {
      
      std::stringstream spec_string_stream(std::stringstream::in |
        std::stringstream::out);

      spec_string_stream << spec_comp+1;
      std::string spec_string_num=spec_string_stream.str();

      ibase_transport++;

      std::string species_str="species";
      species_str+=spec_string_num;
      species_str+=im_string;

      MOFstate_names[ibase_transport]=species_str;
      set_scalar_bc(MOFstate_bcs[ibase_transport],species_phys_bc); 

     }  // spec_comp

     if (ibase_transport!=(im+1)*num_state_material-1)
      amrex::Error("ibase_transport bust");

    } // im (scalar state variables)

    StateDescriptor::BndryFunc MOFstate_fill_class(FORT_STATEFILL,
       FORT_GROUP_STATEFILL);

    int scomp_states=(AMREX_SPACEDIM+1);

    desc_lst.setComponent(State_Type,
     scomp_states,
     MOFstate_names,
     MOFstate_bcs,
     MOFstate_fill_class,&pc_interp);

     // reset the interpolation properties of 
     // velocity, pressure, density, and temperature. 
    override_enable_spectral(enable_spectral);

    Vector<std::string> MOF_names;
    MOF_names.resize(nmat*ngeom_raw);
    Vector<BCRec> MOF_bcs;
    MOF_bcs.resize(nmat*ngeom_raw);

    for (int im=0;im<nmat;im++) {

     int ibase_mof=im*ngeom_raw;

     std::stringstream im_string_stream(std::stringstream::in |
        std::stringstream::out);

     im_string_stream << im+1;
     std::string im_string=im_string_stream.str();

     std::string vof_str="vof"; 
     vof_str+=im_string; 
     MOF_names[ibase_mof]=vof_str;
     set_scalar_vof_bc(MOF_bcs[ibase_mof],phys_bc);

     ibase_mof++;
     std::string cenx_str="cenx"; 
     cenx_str+=im_string; 
     MOF_names[ibase_mof]=cenx_str;
     set_x_vel_bc(MOF_bcs[ibase_mof],phys_bc);

     ibase_mof++;
     std::string ceny_str="ceny"; 
     ceny_str+=im_string; 
     MOF_names[ibase_mof]=ceny_str;
     set_y_vel_bc(MOF_bcs[ibase_mof],phys_bc);

#if (AMREX_SPACEDIM==3)
     ibase_mof++;
     std::string cenz_str="cenz"; 
     cenz_str+=im_string; 
     MOF_names[ibase_mof]=cenz_str;
     set_z_vel_bc(MOF_bcs[ibase_mof],phys_bc);
#endif    

     if (ngeom_raw==NUM_MOF_VAR) {
      amrex::Error("cannot have ngeom_raw=ngeom_recon");
     } else if (ngeom_raw==AMREX_SPACEDIM+1) {
      // do nothing
     } else
      amrex::Error("ngeom_raw invalid");

     if (ibase_mof!=(im+1)*ngeom_raw-1)
      amrex::Error("ibase_mof invalid");

    }  // im  (volume fractions and centroids)

     // GROUP_MOFFILL uses probtype to specify ext_dir.
     // FORT_MOFFILL should never be called.
    StateDescriptor::BndryFunc MOF_fill_class(FORT_MOFFILL,
       FORT_GROUP_MOFFILL);

    int recon_start_pos=scomp_states+nmat*num_state_material;

    multi_mof_interp.multiMOFInterp_nmat=nmat;
    multi_mof_interp.multiMOFInterp_ngeom_raw=ngeom_raw;
    multi_mof_interp.multiMOFInterp_ngeom_recon=ngeom_recon;

    desc_lst.setComponent(State_Type,recon_start_pos,MOF_names,
     MOF_bcs,MOF_fill_class,&multi_mof_interp);

    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,nc-1,"errorind",bc,
      FORT_SCALARFILL,&pc_interp_null);

}  // end subroutine variableSetUp

// post_init_flag==0 if called from post_timestep
// post_init_flag==0 if called from writePlotFile
// post_init_flag==1 if called from post_init_state
// post_init_flag==2 if called from post_restart
void
NavierStokes::sum_integrated_quantities (int post_init_flag) {

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1;

 SDC_outer_sweeps=0;
 SDC_setup_step();

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 if ((adv_dir<1)||(adv_dir>2*AMREX_SPACEDIM+1))
  amrex::Error("adv_dir invalid");
 if (upper_slab_time<0.0)
  amrex::Error("times should be positive");

 std::fflush(NULL);
  // call FLUSH(6), unit=6 screen
 fort_flush_fortran();

 ParallelDescriptor::Barrier();
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();

 if (ParallelDescriptor::IOProcessor()) {
   std::cout << "Starting: sum_integrated_quantities\n";
   std::cout << "This routine takes a while.\n";
   std::cout << "ns.sum_int determines the frequency this\n";
   std::cout << "routine is called\n";
   std::cout << "post_init_flag= " << post_init_flag << '\n';
   std::cout << "upper_slab_time= " << upper_slab_time << '\n';
   std::cout << "adapt_quad_depth= " << adapt_quad_depth << '\n';
 }
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();
 ParallelDescriptor::Barrier();
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();

 if (level!=0)
  amrex::Error("level invalid in sum_integrated_quantities");

 int finest_level = parent->finestLevel();

 NavierStokes& ns_fine = getLevel(finest_level);

 const Real* fine_dx = ns_fine.geom.CellSize();

 int rz_flag=0;
 if (geom.IsRZ())
  rz_flag=1;
 else if (geom.IsCartesian())
  rz_flag=0;
 else if (geom.IsCYLINDRICAL())
  rz_flag=3;
 else
  amrex::Error("geom bust 1");

 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];
 Real problen[AMREX_SPACEDIM];
 Real prob_volume=1.0;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
  problen[dir]=probhi[dir]-problo[dir];
  prob_volume*=problen[dir];
 }

 if (rz_flag==0) {
  // do nothing
 } else if (rz_flag==1) {
  prob_volume*=(2.0*NS_PI*problen[0]/2.0);
 } else if (rz_flag==3) {
  prob_volume*=(2.0*NS_PI);
 } else
  amrex::Error("rz_Flag invalid");


 metrics_dataALL(1);
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);

   // mask=tag if not covered by level+1
   // mask=tag outside the domain
   // mask=tag at coarse-fine ghost cell
   // mask=tag at uncovered fine-fine ghost cell.
   // mask=1-tag if covered by level+1 
  Real tag=1.0;
  int clearbdry=0; 
  ns_level.maskfiner_localMF(MASKCOEF_MF,1,tag,clearbdry);
  ns_level.prepare_mask_nbr(1);

 } // ilev=level..finest_level

 init_FSI_GHOST_MAC_MF_ALL(5);

 build_masksemALL();

 Vector<Real> dt_min;
 dt_min.resize(n_scales+1);
 for (int iscale=0;iscale<dt_min.size();iscale++) {
  dt_min[iscale]=1.0E+30;
 }

 Real vel_max_estdt[AMREX_SPACEDIM+1];
 Real vel_max_cap_wave=0.0;
 int caller_id=3;
 MaxAdvectSpeedALL(dt_min,vel_max_estdt,vel_max_cap_wave,caller_id);

  // 0 empty
  // F,E  2 x nmat
  //
  // drag (3 x nmat comp)
  //
  // min interface location 3 x nmat  (x1,y1,z1   x2,y2,z2  ...)
  // max interface location 3 x nmat  (x1,y1,z1   x2,y2,z2  ...)
  //
  // pressure drag (3 x nmat comp)
  // viscoelastic drag (3 x nmat comp)
  //
  // min den,temp 2 x nmat
  // max den,temp 2 x nmat
  // x=0 amplitude
  // centroid 3 x nmat (x1,y1,z1  x2,y2,z2  ... )
  // min dist from centroid  nmat
  // max dist from centroid  nmat
  // mass      nmat
  // momentum  3 x nmat
  // energy    nmat
  // left pressure, right pressure, left weight, right weight
  // kinetic energy derived  nmat
  // LS F  nmat
  // LS centroid 3 x nmat (x1,y1,z1  x2,y2,z2 ... )
  //
  // torque (3 x nmat comp)
  // pressure torque (3 x nmat comp)
  // viscoelastic torque (3 x nmat comp)
  // perimeter (rasterized) (nmat comp)
  //
  // min interface extent on slice (nmat comp)
  // max interface extent on slice (nmat comp)
  // integral of vorticity (3 comp)
  // vort_error (1 comp)
  // vel_error (1 comp)
  // energy_moment (1 comp)
  // enstrophy (nmat comp)
  // user defined (ncomp_sum_int_user1+ncomp_sum_int_user2 comp)
  // species mass (num_species_var * nmat comp)

 int filler_comp=0;
 int FE_sum_comp=filler_comp+1;

 int drag_sum_comp=FE_sum_comp+2*nmat;

 int minint_sum_comp=drag_sum_comp+3*nmat;
 int maxint_sum_comp=minint_sum_comp+3*nmat;

 int pdrag_sum_comp=maxint_sum_comp+3*nmat;
 int visco_sum_comp=pdrag_sum_comp+3*nmat;

 int minden_sum_comp=visco_sum_comp+3*nmat;
 int maxden_sum_comp=minden_sum_comp+2*nmat;
 int xnot_amp_sum_comp=maxden_sum_comp+2*nmat;
 int cen_sum_comp=xnot_amp_sum_comp+1;
 int mincen_sum_comp=cen_sum_comp+3*nmat;
 int maxcen_sum_comp=mincen_sum_comp+nmat;
 int mass_sum_comp=maxcen_sum_comp+nmat;
 int mom_sum_comp=mass_sum_comp+nmat;
 int energy_sum_comp=mom_sum_comp+3*nmat;
 int left_pressure_sum=energy_sum_comp+nmat;
 int kinetic_energy_sum_comp=left_pressure_sum+4;
 int LS_F_sum_comp=kinetic_energy_sum_comp+nmat;
 int LS_cen_sum_comp=LS_F_sum_comp+nmat;

 int torque_sum_comp=LS_cen_sum_comp+3*nmat;
 int ptorque_sum_comp=torque_sum_comp+3*nmat;
 int viscotorque_sum_comp=torque_sum_comp+3*nmat;
 int step_perim_sum_comp=viscotorque_sum_comp+3*nmat;

 int minint_slice=step_perim_sum_comp+nmat;
 int maxint_slice=minint_slice+nmat;
 int vort_sum_comp=maxint_slice+nmat;
 int vort_error=vort_sum_comp+3;
 int vel_error=vort_error+1;
 int energy_moment=vel_error+1; 
 int enstrophy=energy_moment+1; // integral of w dot w
 int user_comp=enstrophy+nmat;
 int species_mass_comp=user_comp+ncomp_sum_int_user1+ncomp_sum_int_user2;
 int total_comp=species_mass_comp+num_species_var*nmat; 

 Vector<Real> sumdata;
 sumdata.resize(total_comp);
 Vector<int> sumdata_type;
 sumdata_type.resize(total_comp);
 Vector<int> sumdata_sweep;
 sumdata_sweep.resize(total_comp);

 Vector<Real> F_MAT;
 Vector<Real> MASS_MAT;
 F_MAT.resize(nmat);
 MASS_MAT.resize(nmat);

 for (int isum=0;isum<total_comp;isum++) {
  sumdata[isum]=0.0;
  sumdata_type[isum]=1;  // reduce real sum
  sumdata_sweep[isum]=0;  // update first sweep
 }
 for (int im=0;im<ncomp_sum_int_user1;im++) {
  sumdata_sweep[user_comp+im]=0; //update 1st sweep
 }
 for (int im=0;im<ncomp_sum_int_user2;im++) {
  sumdata_sweep[user_comp+ncomp_sum_int_user1+im]=1; //update 2nd sweep
 }

 sumdata_type[vort_error]=3;  // reduce real max (-1.0E+6)
 sumdata_type[vel_error]=3;   // reduce real max (-1.0E+6)

 for (int idir=0;idir<3;idir++) {
  for (int im=0;im<nmat;im++) {
   sumdata_type[idir+minint_sum_comp+3*im]=2; // reduce real min (1.0E+6)
   sumdata_type[idir+maxint_sum_comp+3*im]=3; // reduce real max (-1.0E+6)
  }
 }
 for (int im=0;im<nmat;im++) {
  sumdata_type[minint_slice+im]=2; // reduce real min (1.0E+6)
  sumdata_type[maxint_slice+im]=3; // reduce real max (-1.0E+6)
 }
 for (int idir=0;idir<2*nmat;idir++) {
  sumdata_type[idir+minden_sum_comp]=2;  // reduce real min
  sumdata_type[idir+maxden_sum_comp]=3;  // reduce real max
 }

 sumdata_type[xnot_amp_sum_comp]=3;  // x=0 amplitude  material 1

 for (int idir=0;idir<nmat;idir++) {
  sumdata_type[idir+mincen_sum_comp]=2;  // min dist from centroid
  sumdata_type[idir+maxcen_sum_comp]=3;  // max dist from centroid
  sumdata_sweep[idir+mincen_sum_comp]=1;  
  sumdata_sweep[idir+maxcen_sum_comp]=1; 
 }

 for (int isum=0;isum<total_comp;isum++) {
  sumdata[isum]=0.0;
  if (sumdata_type[isum]==2) // min
   sumdata[isum]=1.0E+6;
  else if (sumdata_type[isum]==3)  // max
   sumdata[isum]=-1.0E+6;
  else if (sumdata_type[isum]==1)
   sumdata[isum]=0.0;
  else
   amrex::Error("sumdata_type invalid");
 } // isum

 int dirx=AMREX_SPACEDIM-1;
 int diry=0;
 int cut_flag=0;
 
 if ((AMREX_SPACEDIM==2)&&
     (probtype==41)&&
     (axis_dir==4)) {
  dirx=AMREX_SPACEDIM-1;
  diry=0;
  cut_flag=1;
 }
 if ((AMREX_SPACEDIM==3)&&(probtype==53)&&(axis_dir==0)) {
  dirx=0;
  diry=AMREX_SPACEDIM-1;
  cut_flag=1;
 }

 Vector<Real> ZZ;
 Vector<Real> FF;
 const Box& fdomain = ns_fine.geom.Domain();
 const int* fdomlo = fdomain.loVect();
 const int* fdomhi = fdomain.hiVect();
 int NN=fdomhi[dirx]-fdomlo[dirx]+1;
 ZZ.resize(NN+1);
 FF.resize(NN+1);
 for (int iz=0;iz<=NN;iz++) {
  ZZ[iz]=0.0;
  FF[iz]=0.0;
 }

 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();
 ParallelDescriptor::Barrier();
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();

 for (int isweep=0;isweep<2;isweep++) {
   // VOF_Recon_ALL 
   // make_physics_varsALL
   // FORT_SUMMASS -> stackerror -> get_symmetric_error -> uses mofdata_tess
  volWgtSumALL(
    post_init_flag,
    sumdata,sumdata_type,sumdata_sweep,
    ZZ,FF,dirx,diry,cut_flag,isweep);

  if (isweep==0) {
   for (int im=0;im<nmat;im++) {
    Real volmat=sumdata[FE_sum_comp+2*im];
    Real LSvolmat=sumdata[LS_F_sum_comp+im];
    if (volmat>0.0) {
     for (int dir=0;dir<AMREX_SPACEDIM;dir++)
      sumdata[3*im+cen_sum_comp+dir]=sumdata[3*im+cen_sum_comp+dir]/volmat;
    }
    if (LSvolmat>0.0) {
     for (int dir=0;dir<AMREX_SPACEDIM;dir++)
      sumdata[3*im+LS_cen_sum_comp+dir]=
       sumdata[3*im+LS_cen_sum_comp+dir]/LSvolmat;
    }
   } // im
  }  // isweep=0
 }

 int f_js=0;
 int f_je=NN-1;
 if (ParallelDescriptor::IOProcessor()) {
  FORT_COFLOW(&upper_slab_time,&f_js,&f_je,&NN,ZZ.dataPtr(),FF.dataPtr(),
    &dirx,&diry,&cut_flag);
 }

 ParallelDescriptor::Barrier();

 Real minpres,maxpres;
 Real maxvel;
 Real maxvel_collide;
 MaxPressureVelocityALL(minpres,maxpres,maxvel,maxvel_collide);

 ParallelDescriptor::Barrier();
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();
 ParallelDescriptor::Barrier();
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();

 Vector<blobclass> blobdata;
 Vector< Vector<Real> > mdot_data;
 Vector< Vector<Real> > mdot_comp_data;
 Vector< Vector<Real> > mdot_data_redistribute;
 Vector< Vector<Real> > mdot_comp_data_redistribute;
 Vector<int> type_flag;

 if (output_drop_distribution==1) {
  int color_count=0;
  int coarsest_level=0;
  int tessellate=1;
  int idx_mdot=-1; //idx_mdot==-1 => do not collect auxiliary data.
  int operation_flag=0;

   // declared in NavierStokes3.cpp
  ColorSumALL(
    operation_flag, // =0
    tessellate, // =1
    coarsest_level,
    color_count,
    TYPE_MF,
    COLOR_MF,
    idx_mdot,
    idx_mdot,
    type_flag,
    blobdata,
    mdot_data,
    mdot_comp_data,
    mdot_data_redistribute,
    mdot_comp_data_redistribute
    );

  if (color_count!=blobdata.size())
   amrex::Error("color_count!=blobdata.size()");
  delete_array(TYPE_MF);
  delete_array(COLOR_MF);
 } else if (output_drop_distribution==0) {
  // do nothing
 } else
  amrex::Error("output_drop_distribution invalid");

 ParallelDescriptor::Barrier();
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();
 ParallelDescriptor::Barrier();

 if (ParallelDescriptor::IOProcessor()) {

  Real smallest_dx=fine_dx[0];
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   std::cout << "TIME= "<<upper_slab_time<<" dir= " << 
     dir << " finest dx=" << fine_dx[dir] << '\n';
   if (smallest_dx>fine_dx[dir])
    smallest_dx=fine_dx[dir];
  }

  if (output_drop_distribution==1) {

   for (int im=0;im<nmat;im++) {
    std::cout << "-----------------------------------------------\n";
    Vector<int> sorted_blob_list;
    int material_blob_count=0;
    for (int iblob=0;iblob<blobdata.size();iblob++) 
     if (blobdata[iblob].im==im+1)
      material_blob_count++;
    if (material_blob_count>0) {
     std::cout << "material im= " << im+1 << "---------------------\n";
     std::cout << "material_blob_count= " << material_blob_count << 
        "---------------------\n";
     sorted_blob_list.resize(material_blob_count);
     int local_count=0;
     for (int iblob=0;iblob<blobdata.size();iblob++) {
      if (blobdata[iblob].im==im+1) {
       sorted_blob_list[local_count]=iblob;
       local_count++;
      }
     }
     if (local_count!=material_blob_count)
      amrex::Error("local_count!=material_blob_count");
     for (int isort1=0;isort1<material_blob_count;isort1++) {
      for (int isort2=0;isort2<material_blob_count-1;isort2++) {
       if (blobdata[sorted_blob_list[isort2]].blob_volume<
           blobdata[sorted_blob_list[isort2+1]].blob_volume) {
        int swapvar=sorted_blob_list[isort2+1];
        sorted_blob_list[isort2+1]=sorted_blob_list[isort2];
        sorted_blob_list[isort2]=swapvar;
       }
      }
     }
     for (int isort1=0;isort1<material_blob_count;isort1++) {
      int iblob=sorted_blob_list[isort1]; 
      Real gvol=blobdata[iblob].blob_volume;
      Real gvol_modify=gvol;
      if (probtype==5700) {
       if (phys_bc.lo(AMREX_SPACEDIM-1)==Symmetry)
        gvol_modify=2.0*gvol_modify;
      }

       // 1..nmat
      int imbase=blobdata[iblob].im;
      if ((imbase<1)||(imbase>nmat))
       amrex::Error("imbase invalid");

      // r=(3V/(4pi))^(1/3)
      // r^3=3V/(4pi)
      // V=4 pi r^3/3
      Real gdiam=2.0*exp(log(3.0*gvol_modify/(4.0*NS_PI))/3.0);
      std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 << 
       " im= " << imbase <<
       " volume = " << gvol_modify << '\n';
      std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 << 
        " im= " << imbase <<     
        " diameter= " << gdiam << " perim= " <<
        blobdata[iblob].blob_perim << '\n';

      for (int imnbr=0;imnbr<nmat;imnbr++) {
       std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
       " im= " << imbase << " imnbr= " << imnbr+1 << 
       " perimnbr= " << 
       blobdata[iblob].blob_perim_mat[imnbr] << '\n';
      } // imnbr=0..nmat-1
      for (int im1=0;im1<nmat;im1++) {
       for (int im2=0;im2<im1;im2++) {
        if ((im1+1!=imbase)&&(im2+1!=imbase)) {
         Real triple_perim=blobdata[iblob].blob_triple_perim[im1][im2];
         if (triple_perim!=0.0) {
          std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
          " im= " << imbase << " im1= " << im1+1 << " im2= " << im2+1 <<
          " triple_perim= " << triple_perim << '\n';
         }  // triple_perim!=0
        }  // im1 and im2 <> imbase
       }  // im2
      }  // im1

      if (gvol>0.0) {
       for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
        std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 << 
        " im= " << imbase << 
        " center dir = " << dir << " coordinate=" << 
        blobdata[iblob].blob_center_integral[dir]/gvol << '\n';
       }
      }  // gvol>0.0

      for (int veltype=0;veltype<3;veltype++) {
       std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
        " im= " << imbase <<
        " veltype= " << veltype << " blob_mass_for_velocity= " <<
        blobdata[iblob].blob_mass_for_velocity[veltype] << '\n';
       for (int dir=0;dir<2*AMREX_SPACEDIM;dir++) {
        std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
        " im= " << imbase <<
        " veltype= " << veltype << " dir= " << dir << " vel= " <<
        blobdata[iblob].blob_velocity[2*AMREX_SPACEDIM*veltype+dir] << '\n';
       } // dir=0..2 sdim-1
      } // veltype=0..2
      for (int dir=0;dir<2*AMREX_SPACEDIM;dir++) {

       std::string momstring="momentum";
       if (dir==0) 
        momstring="momx";
       else if (dir==1)
        momstring="momy";
       else if ((dir==2)&&(AMREX_SPACEDIM==3))
        momstring="momz";
       else if (dir==AMREX_SPACEDIM)
        momstring="momxy";
       else if (dir==AMREX_SPACEDIM+1)
        momstring="momxz";
       else if (dir==AMREX_SPACEDIM+2)
        momstring="momyz";
       else
        amrex::Error("dir invalid");

       std::string velstring="vel";
       if (dir==0) 
        velstring="velx";
       else if (dir==1)
        velstring="vely";
       else if ((dir==2)&&(AMREX_SPACEDIM==3))
        velstring="velz";
       else if (dir==AMREX_SPACEDIM)
        velstring="velxy";
       else if (dir==AMREX_SPACEDIM+1)
        velstring="velxz";
       else if (dir==AMREX_SPACEDIM+2)
        velstring="velyz";
       else
        amrex::Error("dir invalid");


       std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
        " im= " << imbase <<
        " dir= " << dir << " " << momstring << "= " <<
        blobdata[iblob].blob_integral_momentum[dir] 
	<< '\n';
       Real numerator=blobdata[iblob].blob_integral_momentum[dir];
       Real denom=blobdata[iblob].blob_integral_momentum[2*AMREX_SPACEDIM+dir];
       Real avg_vel=numerator;
       if (denom>0.0) 
	avg_vel/=denom;
       std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
        " im= " << imbase <<
        " dir= " << dir << " average " << velstring << "= " << 
        avg_vel << '\n';
       std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
        " im= " << imbase <<
        " dir= " << dir << " " << momstring << " divisor= " << denom << '\n';
      } // dir=0..2 sdim-1
      std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
       " im= " << imbase <<
       " energy= " <<
        blobdata[iblob].blob_energy << '\n';

      Real cellvol=blobdata[iblob].blob_cellvol_count;
      std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
       " im= " << imbase <<
       " blob_cellvol= " <<
        cellvol << '\n';
      if (cellvol>0.0) {
       std::cout << "TIME= " << upper_slab_time << " isort= " << isort1 <<
        " im= " << imbase <<
        " blob_pressure= " << blobdata[iblob].blob_pressure/cellvol
         << '\n';
      } else if (cellvol==0.0) {
       // do nothing
      } else
       amrex::Error("cellvol invalid");

     }  // isort1=0..material_blob_count-1
     std::cout << "-----------------------------------------------\n";
    } // material_blob_count>0
   } // im=0..nmat-1

  } else if (output_drop_distribution==0) {
   // do nothing
  } else
   amrex::Error("output_drop_distribution invalid");

  Real UMACH=0.0;
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    //vel_max_estdt takes into account extra terms from the 
    //cell velocity when spectral element method is turned on.
   std::cout << "TIME= "<<upper_slab_time<<" dir= " << dir << 
     " vel_max_estdt=" << vel_max_estdt[dir] << '\n';
   if (std::abs(vel_max_estdt[dir])>UMACH)
    UMACH=std::abs(vel_max_estdt[dir]);
   std::cout << "TIME= "<<upper_slab_time<<" dir= " << dir <<
	   " AMR_max_phase_change_rate=" <<
	   parent->AMR_max_phase_change_rate[dir] << '\n';
   std::cout << "TIME= "<<upper_slab_time<<" dir= " << dir <<
	   " AMR_min_phase_change_rate=" <<
	   parent->AMR_min_phase_change_rate[dir] << '\n';
  } // dir=0..sdim-1

  std::cout << "TIME= "<<upper_slab_time<< 
     " vel_max_cap_wave=" << vel_max_cap_wave << '\n';

  for (int iten=0;iten<nten;iten++) {
   std::cout << "TIME= "<<upper_slab_time<<" iten= " << iten <<
     " cap_wave_speed=" << cap_wave_speed[iten] << '\n';
  }
  for (int im=0;im<nmat;im++) {
   if (denconst[im]>0.0) {
    Real elastic_wave_speed=elastic_viscosity[im]/denconst[im];
    elastic_wave_speed=sqrt(elastic_wave_speed);
    std::cout << "TIME= "<<upper_slab_time<<" im= " << im <<
     " elastic_wave_speed=" << elastic_wave_speed << '\n';
   } else
    amrex::Error("denconst[im] invalid");
  }

  Real ccsqr=vel_max_estdt[AMREX_SPACEDIM];
  if (ccsqr<0.0)
   amrex::Error("cannot have negative c^2");
  Real USOUND=sqrt(ccsqr);
  if (USOUND>0.0)
   UMACH=UMACH/USOUND;
  else
   UMACH=0.0;

  std::cout << "TIME= "<<upper_slab_time<< " sound_max=" << USOUND << '\n';
  std::cout << "TIME= "<<upper_slab_time<< " max|U|/max|C|=" << UMACH << '\n';

  for (int im=0;im<nmat;im++) {
   F_MAT[im]=sumdata[2*im+FE_sum_comp];

   std::cout <<"TIME= "<< upper_slab_time << " MAT="<<im<<" F=" << 
             F_MAT[im] << '\n';
   std::cout <<"TIME= "<< upper_slab_time << " MAT="<<im<<" LS F=" <<
      sumdata[im+LS_F_sum_comp] << '\n';
   std::cout <<"TIME= "<< upper_slab_time << " MAT="<<im<<" E=" <<
      sumdata[2*im+FE_sum_comp+1] << '\n';
  }
  if (parent->AMR_volume_history_recorded==0) {
   parent->AMR_volume_history.resize(nmat);
   for (int im=0;im<nmat;im++) {
    parent->AMR_volume_history[im]=F_MAT[im];
   } 
   parent->AMR_volume_history_recorded=1;
  } else if (parent->AMR_volume_history_recorded==1) {
   // do nothing
  } else
   amrex::Error("AMR_volume_history_recorded invalid");

  for (int im=0;im<nmat;im++) {
   Real F_ratio=0.0;
   if (parent->AMR_volume_history[im]==0.0) {
    // do nothing
   } else if (parent->AMR_volume_history[im]>0.0) {
    F_ratio=F_MAT[im]/parent->AMR_volume_history[im];
   } else
    amrex::Error("parent->AMR_volume_history[im] invalid");

   std::cout <<"TIME= "<< upper_slab_time << " MAT="<<im<<" FRATIO=" << 
             F_ratio << '\n';

   Real A_ratio=F_ratio;
   if (F_ratio>0.0) {
    if (AMREX_SPACEDIM==2) {
     if ((rz_flag==1)||(rz_flag==3)) {
      A_ratio=exp(log(A_ratio)*2.0/3.0);
     } else if (rz_flag==0) {
      A_ratio=exp(log(A_ratio)/2.0);
     } else
      amrex::Error("rz_flag invalid"); 
    } else if (AMREX_SPACEDIM==3) {
     if ((rz_flag==0)||(rz_flag==3)) {
      A_ratio=exp(log(A_ratio)*2.0/3.0);
     } else
      amrex::Error("rz_flag invalid"); 
    } else
     amrex::Error("dimension bust");
   } else if (F_ratio==0.0) {
    // do nothing
   } else
    amrex::Error("F_ratio invalid");

   std::cout <<"TIME= "<< upper_slab_time << " MAT="<<im<<" ARATIO=" << 
             A_ratio << '\n';
  } // im=0..nmat-1

  for (int im=1;im<=nmat;im++) {
   for (int im_opp=im+1;im_opp<=nmat;im_opp++) {
    for (int ireverse=0;ireverse<=1;ireverse++) {
     if ((im>nmat)||(im_opp>nmat))
      amrex::Error("im or im_opp bust 200cpp");
     int iten,im_source,im_dest;
     get_iten_cpp(im,im_opp,iten,nmat);
     if (iten<1)
      amrex::Error("iten invalid");
     Real LL=latent_heat[iten+ireverse*nten-1];
     if ((ns_is_rigid(im-1)==1)||
         (ns_is_rigid(im_opp-1)==1)) {
      // do nothing
     } else if (LL!=0.0) {
      im_source=im;im_dest=im_opp;
      if (ireverse==0) {
       im_source=im;im_dest=im_opp;
      } else if (ireverse==1) {
       im_source=im_opp;im_dest=im;
      } else
       amrex::Error("ireverse invalid");
      Real ratio=0.0;
      Real denom=parent->AMR_volume_history[im_source-1]-F_MAT[im_source-1];
      if (denom!=0.0) 
       ratio=(F_MAT[im_dest-1]-parent->AMR_volume_history[im_dest-1])/denom;
      std::cout <<"TIME= "<< upper_slab_time << " MAT="<<im_dest-1<<" RATIO="<<
        ratio << '\n';
     } // LL!=0
    } // ireverse
   } // im_opp
  } // im=1..nmat

  for (int im=0;im<nmat;im++) {
   std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" min den=" <<
      sumdata[2*im+minden_sum_comp] << '\n';
   std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" min temp=" <<
      sumdata[2*im+minden_sum_comp+1] << '\n';
   std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" max den=" <<
      sumdata[2*im+maxden_sum_comp] << '\n';
   std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" max temp=" <<
      sumdata[2*im+maxden_sum_comp+1] << '\n';

   MASS_MAT[im]=sumdata[im+mass_sum_comp];

   std::cout<<"TIME= "<<upper_slab_time<<" MAT="<<im<<
	   " mass="<<MASS_MAT[im]<< '\n';

   for (int ispec=0;ispec<num_species_var;ispec++) {
    Real mass_spec=sumdata[im+species_mass_comp+ispec*nmat];
    std::cout<<"TIME= "<<upper_slab_time<<" MAT="<<im<<
     " ispec="<<ispec<< " species mass="<<mass_spec<< '\n';
   }

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" dir= "<<
      dir << " mom=" << sumdata[3*im+dir+mom_sum_comp] << '\n';
   }
   std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" energy=" <<
      sumdata[im+energy_sum_comp] << '\n';
  }
  for (int im=0;im<nmat;im++) { 
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" cendir=" << dir << 
     " centroid=" << sumdata[cen_sum_comp+3*im+dir] << '\n';
   }
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" cendir=" << dir << 
     " LS centroid=" << sumdata[LS_cen_sum_comp+3*im+dir] << '\n';
   }
   std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" mindistcen=" <<
     sumdata[mincen_sum_comp+im] << '\n';
   std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<" maxdistcen=" <<
     sumdata[maxcen_sum_comp+im] << '\n';

   std::cout << "TIME= " << upper_slab_time << " MAT="<<im<<
     " KINETIC ENERGY=" <<
     sumdata[kinetic_energy_sum_comp+im] << '\n';
  } // im=0..nmat-1

  for (int dir=0;dir<3;dir++) {
   std::cout << "TIME= "<<upper_slab_time<<" DIR= " << dir << " VORT SUM " << 
     sumdata[vort_sum_comp+dir] << '\n';
  }
  for (int im=0;im<nmat;im++) {
   std::cout << "TIME= "<<upper_slab_time<<
    "material id (1..nmat) " << im+1 <<
    " ENSTROPHY " << sumdata[enstrophy+im] << '\n';
  }
  for (int im=0;im<ncomp_sum_int_user1;im++) {
   std::cout << "TIME= "<<upper_slab_time<<
    "user_comp1 (1..ncomp_sum_int_user1) " << im+1 <<
    " sum_int_user1 " << sumdata[user_comp+im] << '\n';
  }
  for (int im=0;im<ncomp_sum_int_user2;im++) {
   std::cout << "TIME= "<<upper_slab_time<<
    "user_comp2 (1..ncomp_sum_int_user2) " << im+1 <<
    " sum_int_user2 " << sumdata[user_comp+ncomp_sum_int_user1+im] << '\n';
  }

  std::cout << "TIME= "<<upper_slab_time<<" VORT ERR= " << 
    sumdata[vort_error] << '\n';
  std::cout << "TIME= "<<upper_slab_time<<" VEL ERR= " << 
    sumdata[vel_error] << '\n';

  Real r_moment=0.0;
  Real energy_first_mat=sumdata[kinetic_energy_sum_comp]; 
  if (energy_first_mat>0.0) {
   r_moment=sumdata[energy_moment]/energy_first_mat;
  } else if (energy_first_mat==0.0) {
   // do nothing
  } else
   amrex::Error("energy_first_mat invalid");
  std::cout << "TIME= "<<upper_slab_time<<" ENERGY MOMENT= " << 
    r_moment << '\n';

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   std::cout << "TIME= "<<upper_slab_time<<" DIR= " << dir << " DRAG " << 
     sumdata[drag_sum_comp+dir] << '\n';
  }

  if ((probtype==55)&&(axis_dir==5)&&(AMREX_SPACEDIM==2)&&(nmat==4)) {
   std::cout << "TIME= "<<upper_slab_time<<" F1+F3= " << 
     F_MAT[0]+F_MAT[2] << '\n';
   std::cout << "TIME= "<<upper_slab_time<<" M1+M3= " << 
    MASS_MAT[0]+MASS_MAT[2] << '\n';
  }

   // melting of block of ice.
  if ((probtype==59)&&(AMREX_SPACEDIM==2)&&(nmat==4)) {
   std::cout << "TIME= "<<upper_slab_time<<" F1+F3= " << 
     F_MAT[0]+F_MAT[2] << '\n';
   std::cout << "TIME= "<<upper_slab_time<<" M1+M3= " << 
    MASS_MAT[0]+MASS_MAT[2] << '\n';
  }

  if ((probtype==563)&&(axis_dir==2)&&(AMREX_SPACEDIM==3)) {
   std::cout << "thickness of gear is 3 cm \n";
   Real thick=probhi[1]-problo[1];
   std::cout << "comp thick of gear is " << thick << '\n';
   Real power=sumdata[drag_sum_comp+1]*(3.0/thick)/1.0E7;
   std::cout << "TIME= "<<upper_slab_time<<" predicted power loss " << power << '\n';
   std::cout << "expected power loss 500 watts \n";
  }
  if (probtype==32) {
   Real ff=0.0;
   Real UU=std::abs(adv_vel);
   if (std::abs(advbot)>UU)
    UU=std::abs(advbot);
   if (xblob4>0.0) {
    ff=1.0/xblob4;
    if (std::abs(ff)>UU)
     UU=std::abs(ff);
   }
   if (radblob4>0.0)
    UU=radblob4;

   if ((adv_dir<1)||(adv_dir>AMREX_SPACEDIM))
    amrex::Error("adv_dir invalid");

   Real dcoef=denconst[0]*UU*UU*radblob;
     // pressure and viscous drag
   Real dragcoeff=sumdata[drag_sum_comp+adv_dir-1];  
     // pressure drag
   Real pdragcoeff=sumdata[pdrag_sum_comp+adv_dir-1];  

   if (dcoef!=0.0) {
    dragcoeff/=dcoef;
    pdragcoeff/=dcoef;
    int symmetry_flag=0;

    if (AMREX_SPACEDIM==3) {
     if (adv_dir==1) {
      if (zblob==0.0) {
       if (phys_bc.lo(AMREX_SPACEDIM-1)==Symmetry)
        symmetry_flag=1;
      } else if (yblob==0.0) {
       if (phys_bc.lo(1)==Symmetry)
        symmetry_flag=1;
      } else
       amrex::Error("always run with symmetric bc");
     } else if (adv_dir==2) {
      if (zblob==0.0) {
       if (phys_bc.lo(AMREX_SPACEDIM-1)==Symmetry)
        symmetry_flag=1;
      } else if (xblob==0.0) {
       if (phys_bc.lo(0)==Symmetry)
        symmetry_flag=1;
      } else
       amrex::Error("always run with symmetric bc");
     } else if (adv_dir==3) {
      if (xblob==0.0) {
       if (phys_bc.lo(0)==Symmetry)
        symmetry_flag=1;
      } else if (yblob==0.0) {
       if (phys_bc.lo(1)==Symmetry)
        symmetry_flag=1;
      } else
       amrex::Error("always run with symmetric bc");
     } else
      amrex::Error("adv_dir invalid");
    } else if (AMREX_SPACEDIM==2) {
     if (adv_dir==1) {
      if (yblob==0.0) {
       if (phys_bc.lo(AMREX_SPACEDIM-1)==Symmetry)
        symmetry_flag=1;
      } else
       amrex::Error("always run with symmetric bc");
     } else if (adv_dir==2) {
      if (xblob==0.0) {
       if (phys_bc.lo(0)==Symmetry)
        symmetry_flag=1;
      } else
       amrex::Error("always run with symmetric bc");
     } else
      amrex::Error("adv_dir invalid");
    } else
     amrex::Error("dimension bust"); 

    if (symmetry_flag==1) {
     dragcoeff*=2.0;
     pdragcoeff*=2.0;
    } else if (symmetry_flag!=0)
     amrex::Error("symmetry_flag invalid");

    std::cout << "TIME= " << upper_slab_time << " Cd " << dragcoeff << '\n';
    std::cout << "TIME= " << upper_slab_time << " PCd " << pdragcoeff << '\n';
    std::cout << "Cd computed as F/(rho U^2 diam/2 ) \n";
    std::cout << "(rho)denconst0=" << denconst[0] << '\n';
    std::cout << "U=" << UU << '\n';
    std::cout << "(A/2)radblob=" << radblob << '\n';
   }  // dcoef<>0
  } // probtype=32

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   std::cout << "TIME= "<<upper_slab_time<<" DIR= " << dir << " PDRAG " << 
     sumdata[pdrag_sum_comp+dir] << '\n';
  }
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   std::cout << "TIME= "<<upper_slab_time<<" DIR= " << dir << " VDRAG " << 
     sumdata[drag_sum_comp+dir]-sumdata[pdrag_sum_comp+dir] << '\n';
  }
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   std::cout << "TIME= "<<upper_slab_time<<" DIR= " << dir << " PTORQUE " <<
     sumdata[ptorque_sum_comp+dir] << '\n';
  }
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   std::cout << "TIME= "<<upper_slab_time<<" DIR= " << dir << " VTORQUE " <<
     sumdata[torque_sum_comp+dir]-sumdata[ptorque_sum_comp+dir] << '\n';
  }
  std::cout << "TIME= "<<upper_slab_time<< " HEIGHT_PERIM " <<
    sumdata[step_perim_sum_comp] << '\n';


  for (int im=0;im<nmat;im++) {
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

    std::cout << "TIME=" << upper_slab_time << " MAT="<<im<<
     " dir= " << dir << " GLOBAL MIN INT=" << 
     sumdata[dir+3*im+minint_sum_comp] << '\n';

    std::cout << "TIME=" << upper_slab_time << " MAT="<<im<<
     " dir= " << dir << " GLOBAL MAX INT=" << 
     sumdata[dir+3*im+maxint_sum_comp] << '\n';

    std::cout << "TIME=" << upper_slab_time << " MAT="<<im<<
     " dir= " << dir << " SEPARATION=" << 
     sumdata[dir+3*im+maxint_sum_comp]-
     sumdata[dir+3*im+minint_sum_comp] << '\n';

   }  // dir


   std::cout << "TIME=" << upper_slab_time << " MAT="<<im<<
    " SLICE MIN INT=" << 
    sumdata[im+minint_slice] << '\n';

   std::cout << "TIME=" << upper_slab_time << " MAT="<<im<<
    " SLICE MAX INT=" << 
    sumdata[im+maxint_slice] << '\n';

   std::cout << "TIME=" << upper_slab_time << " MAT="<<im<<
    " SLICE SEPARATION=" << 
    sumdata[im+maxint_slice]-
    sumdata[im+minint_slice] << '\n';

  }  // im
  Real offset=yblob;
  if (AMREX_SPACEDIM==3)
   offset=zblob;

  std::cout << "TIME=" << upper_slab_time << " FREE AMPLITUDE=" <<
      sumdata[xnot_amp_sum_comp]-offset << '\n';

  std::cout << "TIME= " << upper_slab_time << " MINPRES=  " << minpres << '\n';
  std::cout << "TIME= " << upper_slab_time << " MAXPRES=  " << maxpres << '\n';
  std::cout << "TIME= " << upper_slab_time << " MAXVEL=  " << maxvel << '\n';
  std::cout << "TIME= " << upper_slab_time << " MAXVEL COLLIDE=  " 
   << maxvel_collide << '\n';

  for (int iscale=0;iscale<dt_min.size();iscale++) {
   std::cout << "TIME= " << upper_slab_time << " iscale=  " << iscale << ' '
	   << " dt_min= " << dt_min[iscale] << '\n';
  }

  Real leftwt=sumdata[left_pressure_sum+2];
  Real rightwt=sumdata[left_pressure_sum+3];
  if ((leftwt<=0.0)||(rightwt<=0.0))
   amrex::Error("leftwt or rightwt are invalid");
  Real leftpres=sumdata[left_pressure_sum]/leftwt;
  Real rightpres=sumdata[left_pressure_sum+1]/rightwt;
  std::cout << "TIME= " << upper_slab_time << 
   " LEFTPRES=  " << leftpres << '\n';
  std::cout << "TIME= " << upper_slab_time << 
   " RIGHTPRES=  " << rightpres << '\n';

   // bubble jetting problem
  if (probtype==42) {
   Real bubble_volume=sumdata[FE_sum_comp+2];
   Real radbubble=exp(log(3.0*bubble_volume/(4.0*NS_PI))/3.0);
   std::cout << "TIME= " << upper_slab_time << " JETTINGVOL=  " << 
    radbubble << '\n';
  }

   // inputs.circular_freeze
   // pi r^2/4=F
   // r=sqrt(4F/pi)
  if ((probtype==801)&&(axis_dir==3)) {
   Real bubble_volume=sumdata[FE_sum_comp+2];
   Real radbubble=0.0;
   if (rz_flag==0) {
    radbubble=sqrt(4.0*bubble_volume/NS_PI);

   // 4/3 pi r^3 = 2V
   // r=(3V/(2 pi))^{1/3}
   } else if (rz_flag==1) {
    radbubble=exp(log(3.0*bubble_volume/(2.0*NS_PI))/3.0);
   } else
    amrex::Error("rz_flag invalid");

   std::cout << "TIME= " << upper_slab_time << " EFFECTIVE RAD=  " << 
    radbubble << '\n';
  }

   // Sato and Niceno or  Tryggvason and Lu test cases.
  if ((probtype==55)&&
      ((axis_dir==6)||    // incompressible boiling
       (axis_dir==7))) {  // compressible boiling

   Real bubble_volume=sumdata[FE_sum_comp+2];
   Real radbubble=0.0;

   if (AMREX_SPACEDIM==2) {
    // 4/3 pi r^3 = V
    // r=(3V/(4 pi))^{1/3}
    if (rz_flag==1) {
     radbubble=exp(log(3.0*bubble_volume/(4.0*NS_PI))/3.0);
    } else if (rz_flag==0) {
    // pi r^2 = 2V
     radbubble=sqrt(2.0*bubble_volume/NS_PI);
    } else
     amrex::Error("rz_flag invalid");
   } else if (AMREX_SPACEDIM==3) {
    radbubble=exp(log(3.0*bubble_volume/(4.0*NS_PI))/3.0);
   } else
    amrex::Error("sdim bust");

   std::cout << "TIME= " << upper_slab_time << " EFFECTIVE RAD=  " << 
    radbubble << '\n';
  } // probtype==55 and axis_dir==6 or 7

  std::cout << "TIME= " << upper_slab_time << " number_mfiter_loops=  " << 
    thread_class::number_mfiter_loops << '\n';

 }  // io processor
 ParallelDescriptor::Barrier();
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();
 ParallelDescriptor::Barrier();

 delete_array(MASKCOEF_MF);

 ParallelDescriptor::Barrier();
 std::fflush(NULL);
  // call FLUSH(6)
 fort_flush_fortran();
 ParallelDescriptor::Barrier();

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "done sum_integrated_quantities\n";
  }
} // end subroutine sum_integrated_quantities

}/* namespace amrex */
