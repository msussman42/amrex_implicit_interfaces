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

! probtype==401 (see run3d/inputs.HELIX)  3 materials
module HELIX_module
use amrex_fort_module, only : amrex_real

implicit none                   

contains

   ! do any initial preparation needed
 subroutine INIT_HELIX_MODULE()
 IMPLICIT NONE

 return
 end subroutine INIT_HELIX_MODULE


  ! fluids tessellate the domain, solids are immersed. 
 subroutine HELIX_LS(x,t,LS)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 integer im
 real(amrex_real) LS(num_materials)

  ! fluid materials tessellate the domain
 if ((num_materials.eq.3).and.(probtype.eq.401)) then
  do im=1,num_materials
   if (im.eq.1) then !liquid
    LS(im)=sqrt( (x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)-radblob
   else if (im.eq.2) then !gas
    LS(im)=-sqrt( (x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)+radblob
   else if (im.eq.3) then ! geometry (placeholder)
    LS(im)=-99999.0
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
 end subroutine HELIX_LS

 subroutine HELIX_VEL(x,t,LS,VEL,velsolid_flag)
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

 if ((adv_dir.ge.1).and.(adv_dir.le.SDIM)) then
  do dir=1,SDIM
   VEL(dir)=zero
  enddo
  VEL(adv_dir)=adv_vel
 else
  print *,"adv_dir invalid in HELIX_VEL"
  stop
 endif

 return 
 end subroutine HELIX_VEL

 subroutine HELIX_PRES(x,t,LS,PRES)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) PRES

 PRES=outflow_pressure

 return 
 end subroutine HELIX_PRES


 subroutine HELIX_STATE(x,t,LS,STATE)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real), INTENT(in) :: x(SDIM)
 real(amrex_real), INTENT(in) :: t
 real(amrex_real), INTENT(in) :: LS(num_materials)
 real(amrex_real), INTENT(out) :: STATE(num_materials*num_state_material)
 integer im,ibase,n

 if ((num_materials.eq.3).and. &
     (num_state_material.ge.2).and. &
     (probtype.eq.401)) then
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
 end subroutine HELIX_STATE

  ! dir=1..sdim  side=1..2
 subroutine HELIX_LS_BC(xwall,xghost,t,LS, &
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
  call HELIX_LS(xghost,t,LS)
 else
  print *,"dir or side invalid"
  stop
 endif
 
 return
 end subroutine HELIX_LS_BC


  ! dir=1..sdim  side=1..2 veldir=1..sdim
 subroutine HELIX_VEL_BC(xwall,xghost,t,LS, &
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

  call HELIX_VEL(xghost,t,LS,local_VEL,velsolid_flag)
  VEL=local_VEL(veldir)

 else
  print *,"dir,side, or veldir invalid"
  stop
 endif

 return
 end subroutine HELIX_VEL_BC


  ! dir=1..sdim  side=1..2
 subroutine HELIX_PRES_BC(xwall,xghost,t,LS, &
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

  call HELIX_PRES(xghost,t,LS,PRES)

 else
  print *,"dir or side invalid"
  stop
 endif

 return
 end subroutine HELIX_PRES_BC

 function is_HELIX_overlay(nmat,im)
 use probcommon_module
 IMPLICIT NONE

 integer is_HELIX_overlay
 integer nmat,im

 if (nmat.eq.num_materials) then
  if (num_materials.eq.3) then
   if ((im.ge.1).and.(im.le.nmat)) then
    if (im.eq.3) then 
     is_HELIX_overlay=1
    else
     is_HELIX_overlay=0
    endif
   else
    print *,"im invalid in is_HELIX_overlay"
    stop
   endif
  else
   print *,"num_materials invalid in is_HELIX_overlay"
   stop
  endif
 else
  print *,"nmat invalid in is_HELIX_overlay"
  stop
 endif
  
 return
 end function is_HELIX_overlay

  ! dir=1..sdim  side=1..2
 subroutine HELIX_STATE_BC(xwall,xghost,t,LS, &
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
 integer ibase,im_crit,im_loop

 if ((istate.ge.1).and. &
     (istate.le.num_state_material).and. &
     (im.ge.1).and. &
     (im.le.num_materials)) then
  call HELIX_STATE(xghost,t,LS,local_STATE)
  ibase=(im-1)*num_state_material
  STATE=local_STATE(ibase+istate)
  call get_primary_material(LS,im_crit)

  do im_loop=1,num_materials
   if (is_HELIX_overlay(num_materials,im_loop).eq.1) then
    if (LS(im_loop).ge.zero) then
     im_crit=im_loop
    else if (LS(im_loop).le.zero) then
     ! do nothing
    else
     print *,"LS(im_loop) invalid"
     stop
    endif
   else if (is_HELIX_overlay(num_materials,im_loop).eq.0) then
    ! do nothing
   else
    print *,"is_HELIX_overlay(num_materials,im_loop) invalid"
    stop
   endif
  enddo ! im_loop=1,num_materials

  ibase=(im_crit-1)*num_state_material
  STATE_merge=local_STATE(ibase+istate)
 else
  print *,"istate invalid"
  stop
 endif

 return
 end subroutine HELIX_STATE_BC

 subroutine HELIX_HEATSOURCE(im,VFRAC,time,x,temp, &
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

 if ((num_materials.eq.3).and.(probtype.eq.401)) then
  heat_source=zero
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine HELIX_HEATSOURCE

end module HELIX_module
