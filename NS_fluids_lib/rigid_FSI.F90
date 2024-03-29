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

! probtype==533 (see run2d/inputs.splashing_rigid_objectRZ or
!                    run3d/inputs.splashing_rigid_objectXYZ)
module rigid_FSI_module
use amrex_fort_module, only : amrex_real

implicit none                   

contains

   ! do any initial preparation needed
 subroutine INIT_rigid_FSI_MODULE()
 IMPLICIT NONE

 return
 end subroutine INIT_rigid_FSI_MODULE

  ! fluids tessellate the domain, solids are immersed. 
  ! The FSI rigid material (FSI_flag==FSI_RIGID_NOTPRESCRIBED) 
  ! is treated as a fluid.
 subroutine rigid_FSI_LS(x,t,LS)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 integer im
 real(amrex_real) LS(num_materials)
 real(amrex_real) LS_particle

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

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) PRES

 PRES=outflow_pressure

 return 
 end subroutine rigid_FSI_PRES


 subroutine rigid_FSI_STATE(x,t,LS,STATE)
 use probcommon_module
 IMPLICIT NONE

 real(amrex_real) x(SDIM)
 real(amrex_real) t
 real(amrex_real) LS(num_materials)
 real(amrex_real) STATE(num_materials*num_state_material)
 integer im,ibase,n

 if ((num_materials.eq.3).and. &
     (num_state_material.ge.2).and. &
     (probtype.eq.533)) then
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
 end subroutine rigid_FSI_STATE

  ! dir=1..sdim  side=1..2
 subroutine rigid_FSI_LS_BC(xwall,xghost,t,LS, &
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
  call rigid_FSI_STATE(xghost,t,LS,local_STATE)
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
 end subroutine rigid_FSI_STATE_BC

 subroutine rigid_FSI_HEATSOURCE(im,VFRAC,time,x,temp, &
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

 if ((num_materials.eq.3).and.(probtype.eq.533)) then
  heat_source=zero
 else
  print *,"num_materials or probtype invalid"
  stop
 endif

 return
 end subroutine rigid_FSI_HEATSOURCE

end module rigid_FSI_module
