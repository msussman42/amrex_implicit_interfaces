#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"


#if (BL_SPACEDIM==3)
#define SDIM 3
#elif (BL_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

! probtype==533 (see run2d/inputs.splashing_rigid_objectRZ or
!                    run3d/inputs.splashing_rigid_objectXYZ)
module rigid_FSI_module

implicit none                   

contains

   ! do any initial preparation needed
 subroutine INIT_rigid_FSI_MODULE()
 IMPLICIT NONE

 return
 end subroutine INIT_rigid_FSI_MODULE

  ! fluids tessellate the domain, solids are immersed. 
  ! The FSI rigid material (FSI_flag==5) is treated as a fluid.
 subroutine rigid_FSI_LS(x,t,LS)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 INTEGER_T im
 REAL_T LS(num_materials)
 REAL_T LS_particle

 if ((num_materials.eq.3).and.(probtype.eq.533)) then

  if (SDIM.eq.2) then
   LS_particle=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
  else if (SDIM.eq.3) then
   LS_particle=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
  else
   print *,"dimension bust"
   stop
  endif

  do im=1,num_materials
   if (im.eq.1) then !liquid
    LS(im)=radblob2-x(SDIM)
   else if (im.eq.2) then ! gas
    LS(im)=x(SDIM)-radblob2
    if (abs(LS(im)).ge.abs(LS_particle)) then
     LS(im)=-LS_particle
    endif
   else if (im.eq.3) then ! rigid object
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
 end subroutine rigid_FSI_LS

 subroutine rigid_FSI_VEL(x,t,LS,VEL,velsolid_flag)
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
  if ((LS(3).ge.-radblob/100.0D0).or. &
      (velsolid_flag.eq.1)) then
   VEL(SDIM)=-abs(vinletgas)
  else if ((LS(3).le.-radblob/100.0D0).and. &
           (velsolid_flag.eq.0)) then
   ! do nothing
  else
   print *,"LS(3) invalid"
   stop
  endif
 else
  print *,"adv_dir invalid in rigid_FSI_VEL"
  stop
 endif

 return 
 end subroutine rigid_FSI_VEL

 subroutine rigid_FSI_PRES(x,t,LS,PRES)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T PRES

 PRES=outflow_pressure

 return 
 end subroutine rigid_FSI_PRES


 subroutine rigid_FSI_STATE(x,t,LS,STATE)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T STATE(num_materials*num_state_material)
 INTEGER_T im,ibase,n

 if ((num_materials.eq.3).and. &
     (num_state_material.ge.2).and. &
     (probtype.eq.533)) then
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
 end subroutine rigid_FSI_STATE

  ! dir=1..sdim  side=1..2
 subroutine rigid_FSI_LS_BC(xwall,xghost,t,LS, &
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
  call rigid_FSI_LS(xghost,t,LS)
 else
  print *,"dir or side invalid"
  stop
 endif
 
 return
 end subroutine rigid_FSI_LS_BC


  ! dir=1..sdim  side=1..2 veldir=1..sdim
 subroutine rigid_FSI_VEL_BC(xwall,xghost,t,LS, &
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

  call rigid_FSI_VEL(xghost,t,LS,local_VEL,velsolid_flag)
  VEL=local_VEL(veldir)

 else
  print *,"dir,side, or veldir invalid"
  stop
 endif

 return
 end subroutine rigid_FSI_VEL_BC


  ! dir=1..sdim  side=1..2
 subroutine rigid_FSI_PRES_BC(xwall,xghost,t,LS, &
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

  call rigid_FSI_PRES(xghost,t,LS,PRES)

 else
  print *,"dir or side invalid"
  stop
 endif

 return
 end subroutine rigid_FSI_PRES_BC

  ! dir=1..sdim  side=1..2
 subroutine rigid_FSI_STATE_BC(xwall,xghost,t,LS, &
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
  call rigid_FSI_STATE(xghost,t,LS,local_STATE)
  ibase=(im-1)*num_state_material
  STATE=local_STATE(ibase+istate)
  im_crit=1
  do im_loop=2,num_materials
   if (LS(im_loop).gt.LS(im_crit)) then
    im_crit=im_loop
   endif
  enddo
  ibase=(im_crit-1)*num_state_material
  STATE_merge=local_STATE(ibase+istate)
 else
  print *,"istate invalid"
  stop
 endif

 return
 end subroutine rigid_FSI_STATE_BC

 subroutine rigid_FSI_HEATSOURCE(im,VFRAC,time,x,temp, &
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

 if ((num_materials.eq.3).and.(probtype.eq.533)) then
  heat_source=zero
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine rigid_FSI_HEATSOURCE

end module rigid_FSI_module
