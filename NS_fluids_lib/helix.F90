#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#include "AMReX_ArrayLim.H"


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

 REAL_T x(SDIM)
 REAL_T t
 INTEGER_T im
 REAL_T LS(num_materials)

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

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T VEL(SDIM)
 INTEGER_T dir
 INTEGER_T velsolid_flag

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

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T PRES

 PRES=outflow_pressure

 return 
 end subroutine HELIX_PRES


 subroutine HELIX_STATE(x,t,LS,STATE)
 use probcommon_module
 IMPLICIT NONE

 REAL_T, intent(in) :: x(SDIM)
 REAL_T, intent(in) :: t
 REAL_T, intent(in) :: LS(num_materials)
 REAL_T, intent(out) :: STATE(num_materials*num_state_material)
 INTEGER_T im,ibase,n

 if ((num_materials.eq.3).and. &
     (num_state_material.ge.2).and. &
     (probtype.eq.401)) then
  do im=1,num_materials

   ibase=(im-1)*num_state_material
   STATE(ibase+1)=fort_denconst(im)

   if (t.eq.zero) then
    STATE(ibase+2)=fort_initial_temperature(im)
   else if (t.gt.zero) then
    STATE(ibase+2)=fort_tempconst(im)
   else
    print *,"t invalid"
    stop
   endif

   do n=1,num_species_var
    STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
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

 REAL_T xwall
 REAL_T xghost(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T LS_in(num_materials)
 INTEGER_T dir,side
 REAL_T dx(SDIM)

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

 REAL_T xwall
 REAL_T xghost(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T VEL
 REAL_T VEL_in
 INTEGER_T veldir,dir,side
 REAL_T dx(SDIM)
 REAL_T local_VEL(SDIM)
 INTEGER_T velsolid_flag

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

 REAL_T xwall
 REAL_T xghost(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T PRES
 REAL_T PRES_in
 INTEGER_T dir,side
 REAL_T dx(SDIM)

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

 INTEGER_T is_HELIX_overlay
 INTEGER_T nmat,im

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
 IMPLICIT NONE

 REAL_T xwall
 REAL_T xghost(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T local_STATE(num_materials*num_state_material)
 REAL_T STATE
 REAL_T STATE_merge
 REAL_T STATE_in
 INTEGER_T dir,side
 REAL_T dx(SDIM)
 INTEGER_T istate,im
 INTEGER_T ibase,im_crit,im_loop

 if ((istate.ge.1).and. &
     (istate.le.num_state_material).and. &
     (im.ge.1).and. &
     (im.le.num_materials)) then
  call HELIX_STATE(xghost,t,LS,local_STATE)
  ibase=(im-1)*num_state_material
  STATE=local_STATE(ibase+istate)
  im_crit=1
  do im_loop=2,num_materials
   if (LS(im_loop).gt.LS(im_crit)) then
    im_crit=im_loop
   endif
  enddo

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

 INTEGER_T im
 REAL_T VFRAC(num_materials)
 REAL_T time
 REAL_T x(SDIM)
 REAL_T temp(num_materials)
 REAL_T den(num_materials)
 REAL_T CV(num_materials)
 REAL_T dt
 REAL_T heat_source

 if ((num_materials.eq.3).and.(probtype.eq.401)) then
  heat_source=zero
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine HELIX_HEATSOURCE

end module HELIX_module
