#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#include "AMReX_ArrayLim.H"
#include "EXTRAP_COMP.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

! probtype==534 (see run2d/inputs.sinking_particle)
module sinking_particle_module
use amrex_fort_module, only : amrex_real

implicit none                   

contains

  ! fluids tessellate the domain, solids are immersed. 
  ! The FSI rigid material (FSI_flag==FSI_RIGID_NOTPRESCRIBED)         
  ! is treated as a fluid.
 subroutine sinking_FSI_LS(x,t,LS)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 integer im
 real(amrex_real) LS(num_materials)
 real(amrex_real) LS_particle

 if ((num_materials.eq.2).and.(probtype.eq.534)) then

  if (SDIM.eq.2) then
   LS_particle=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
  else if (SDIM.eq.3) then
   LS_particle=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
  else
   print *,"dimension bust"
   stop
  endif

  do im=1,num_materials
   if (im.eq.1) then 
    LS(im)=-LS_particle
   else if (im.eq.2) then ! rigid object
    LS(im)=LS_particle
   else
    print *,"im invalid"
    stop
   endif
  enddo ! im=1..num_materials
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine sinking_FSI_LS

 subroutine sinking_FSI_VEL(x,t,LS,VEL,velsolid_flag)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) VEL(SDIM)
 integer dir
 integer velsolid_flag

 if ((velsolid_flag.eq.0).or. &
     (velsolid_flag.eq.1)) then
  ! do nothing
 else
  print *,"velsolid_flag invalid"
  stop
 endif

 do dir=1,SDIM
  VEL(dir)=zero
 enddo

 return 
 end subroutine sinking_FSI_VEL

 subroutine sinking_FSI_PRES(x,t,LS,PRES)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) PRES

 PRES=outflow_pressure

 return 
 end subroutine sinking_FSI_PRES


 subroutine sinking_FSI_STATE(x,t,LS,STATE)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) STATE(num_materials*num_state_material)
 integer im,ibase,n

 if ((num_materials.eq.2).and. &
     (num_state_material.ge.2).and. &
     (probtype.eq.534)) then
  do im=1,num_materials
   ibase=(im-1)*num_state_material
   STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im)
   if (t.eq.zero) then
    STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im)
   else if (t.gt.zero) then
    STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
   else
    print *,"t invalid"
    stop
   endif
   do n=1,num_species_var
    STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
   enddo
  enddo ! im=1..num_materials
 else
  print *,"num_materials,num_state_material, or probtype invalid"
  stop
 endif
  
 return
 end subroutine sinking_FSI_STATE

  ! dir=1..sdim  side=1..2
 subroutine sinking_FSI_LS_BC(xwall,xghost,t,LS, &
    LS_in,dir,side,dx)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) xwall
 real(amrex_real) xghost(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) LS_in(num_materials)
 integer dir,side
 real(amrex_real) dx(SDIM)

 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2)) then
  call sinking_FSI_LS(xghost,t,LS)
 else
  print *,"dir or side invalid"
  stop
 endif
 
 return
 end subroutine sinking_FSI_LS_BC


  ! dir=1..sdim  side=1..2 veldir=1..sdim
 subroutine sinking_FSI_VEL_BC(xwall,xghost,t,LS, &
    VEL,VEL_in,veldir,dir,side,dx)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) xwall
 real(amrex_real) xghost(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) VEL
 real(amrex_real) VEL_in
 integer veldir,dir,side
 real(amrex_real) dx(SDIM)
 real(amrex_real) local_VEL(SDIM)
 integer velsolid_flag

 velsolid_flag=0
 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2).and. &
     (veldir.ge.1).and.(veldir.le.SDIM)) then

  call sinking_FSI_VEL(xghost,t,LS,local_VEL,velsolid_flag)
  VEL=local_VEL(veldir)

 else
  print *,"dir,side, or veldir invalid"
  stop
 endif

 return
 end subroutine sinking_FSI_VEL_BC


  ! dir=1..sdim  side=1..2
 subroutine sinking_FSI_PRES_BC(xwall,xghost,t,LS, &
    PRES,PRES_in,dir,side,dx)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) xwall
 real(amrex_real) xghost(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) PRES
 real(amrex_real) PRES_in
 integer dir,side
 real(amrex_real) dx(SDIM)

 if ((dir.ge.1).and.(dir.le.SDIM).and. &
     (side.ge.1).and.(side.le.2)) then

  call sinking_FSI_PRES(xghost,t,LS,PRES)

 else
  print *,"dir or side invalid"
  stop
 endif

 return
 end subroutine sinking_FSI_PRES_BC

  ! dir=1..sdim  side=1..2
 subroutine sinking_FSI_STATE_BC(xwall,xghost,t,LS, &
    STATE,STATE_merge,STATE_in,im,istate,dir,side,dx)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE

 real(amrex_real) xwall
 real(amrex_real) xghost(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) local_STATE(num_materials*num_state_material)
 real(amrex_real) STATE
 real(amrex_real) STATE_merge
 real(amrex_real) STATE_in
 integer dir,side
 real(amrex_real) dx(SDIM)
 integer istate,im
 integer ibase,im_crit

 if ((istate.ge.1).and. &
     (istate.le.num_state_material).and. &
     (im.ge.1).and. &
     (im.le.num_materials)) then
  call sinking_FSI_STATE(xghost,t,LS,local_STATE)
  ibase=(im-1)*num_state_material
  STATE=local_STATE(ibase+istate)
  call get_primary_material(LS,im_crit)
  ibase=(im_crit-1)*num_state_material
  STATE_merge=local_STATE(ibase+istate)
 else
  print *,"istate invalid"
  stop
 endif

 return
 end subroutine sinking_FSI_STATE_BC

end module sinking_particle_module
