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

! probtype==311 (see run2d/inputs.Dhir_user_defined)
module USERDEF_module
use amrex_fort_module, only : amrex_real

implicit none                   

contains

   ! do any initial preparation needed
 subroutine INIT_USERDEF_MODULE()
 IMPLICIT NONE

 return
 end subroutine INIT_USERDEF_MODULE

  ! fluids tessellate the domain, solids are immersed. 
 subroutine USERDEF_LS(x,t,LS)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 integer im
 real(amrex_real) LS(num_materials)

 if ((num_materials.eq.4).and.(probtype.eq.311)) then
  do im=1,num_materials
   if (im.eq.1) then !liquid
    LS(im)=sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)-radblob
   else if (im.eq.2) then ! vapor
    LS(im)=-(sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)-radblob)
   else if (im.eq.3) then ! substrate
    LS(im)=yblob2-x(2)
   else if (im.eq.4) then ! torus
    LS(im)=radblob5-sqrt((x(1)-xblob5)**2+(x(2)-yblob5)**2)
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
 end subroutine USERDEF_LS

 subroutine USERDEF_VEL(x,t,LS,VEL,velsolid_flag)
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
 end subroutine USERDEF_VEL

 subroutine USERDEF_PRES(x,t,LS,PRES)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) PRES
 real(amrex_real) :: gravity_dz
 integer :: gravity_dir

 call fort_derive_gravity_dir(gravity_vector,gravity_dir)

 if (gravity_dir.eq.1) then
  gravity_dz=x(gravity_dir)-probhix
 else if (gravity_dir.eq.2) then
  gravity_dz=x(gravity_dir)-probhiy
 else if ((gravity_dir.eq.SDIM).and.(SDIM.eq.3)) then
  gravity_dz=x(gravity_dir)-probhiz
 else
  print *,"dimension bust"
  stop
 endif

 PRES=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*gravity_dz

 return 
 end subroutine USERDEF_PRES


 subroutine USERDEF_STATE(x,t,LS,STATE)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) STATE(num_materials*num_state_material)
 integer im,ibase

 if ((num_materials.eq.4).and. &
     (num_state_material.eq.2).and. &
     (probtype.eq.311)) then
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
  enddo ! im=1..num_materials
 else
  print *,"num_materials,num_state_material, or probtype invalid"
  stop
 endif
  
 return
 end subroutine USERDEF_STATE

  ! dir=1..sdim  side=1..2
 subroutine USERDEF_LS_BC(xwall,xghost,t,LS, &
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
  call USERDEF_LS(xghost,t,LS)
 else
  print *,"dir or side invalid"
  stop
 endif
 
 return
 end subroutine USERDEF_LS_BC


  ! dir=1..sdim  side=1..2 veldir=1..sdim
 subroutine USERDEF_VEL_BC(xwall,xghost,t,LS, &
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

  call USERDEF_VEL(xghost,t,LS,local_VEL,velsolid_flag)
  VEL=local_VEL(veldir)

 else
  print *,"dir,side, or veldir invalid"
  stop
 endif

 return
 end subroutine USERDEF_VEL_BC


  ! dir=1..sdim  side=1..2
 subroutine USERDEF_PRES_BC(xwall,xghost,t,LS, &
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

  call USERDEF_PRES(xghost,t,LS,PRES)

 else
  print *,"dir or side invalid"
  stop
 endif

 return
 end subroutine USERDEF_PRES_BC

  ! dir=1..sdim  side=1..2
 subroutine USERDEF_STATE_BC(xwall,xghost,t,LS, &
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
  call USERDEF_STATE(xghost,t,LS,local_STATE)
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
 end subroutine USERDEF_STATE_BC

 subroutine USERDEF_HEATSOURCE(im,VFRAC,time,x,temp, &
      heat_source,den,CV,dt)
 use probcommon_module
 IMPLICIT NONE

 integer im
 real(amrex_real) VFRAC(num_materials)
 real(amrex_real) time
 real(amrex_real) x(SDIM)
 real(amrex_real) temp(num_materials)
 real(amrex_real) den(num_materials)
 real(amrex_real) CV(num_materials)
 real(amrex_real) dt
 real(amrex_real) heat_source
 real(amrex_real) userdef_temperature
 integer userdef_im

 if ((num_materials.eq.4).and.(probtype.eq.311)) then
  userdef_im=4
  userdef_temperature=fort_tempconst(userdef_im)
  if ((VFRAC(userdef_im).gt.VOFTOL).and. &
      (VFRAC(userdef_im).le.one+EPS1)) then
   if (temp(im).lt.userdef_temperature) then
    heat_source=(userdef_temperature-temp(im))* &
             den(im)*CV(im)/dt 
   else if (temp(im).ge.userdef_temperature) then
    heat_source=zero
   else
    print *,"temp(im) or userdef_temp invalid"
    stop
   endif 
  else if (abs(VFRAC(userdef_im)).le.VOFTOL) then
   heat_source=zero
  else
   print *,"VFRAC invalid"
   stop
  endif
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine USERDEF_HEATSOURCE

end module USERDEF_module
