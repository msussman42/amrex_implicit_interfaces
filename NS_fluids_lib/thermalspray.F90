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

! probtype==402 (see run3d/inputs.thermalspray)  
! material 1 = melt
! material 2 = air
! material 3 = solidified melt
! material 4 = substrate
module TSPRAY_module

implicit none                   


        REAL_T, allocatable, dimension(:,:) :: drop_data
        REAL_T, dimension(3) :: TDOMAIN
        INTEGER_T :: T_num_drops

contains

   ! do any initial preparation needed
 subroutine INIT_TSPRAY_MODULE()
 use probcommon_module
 IMPLICIT NONE

 INTEGER_T nd

 if ((num_materials.eq.4).and. &
     (probtype.eq.402).and. &
     (SDIM.eq.3)) then
  print *,"opening: inputs_data_file_Zeyu"
   ! this file has to be able to be opened by multiple processes.
  open(unit=2,file='inputs_data_file_Zeyu')
  print *,"reading expected domain size TDOMAIN"
  read(2,*) TDOMAIN(1),TDOMAIN(2),TDOMAIN(3)
  print *,"TDOMAIN= ",TDOMAIN(1),TDOMAIN(2),TDOMAIN(3)
  if ((abs(TDOMAIN(1)-probhix).le.1.0E-10).and. &
      (abs(TDOMAIN(2)-probhiy).le.1.0E-10).and. &
      (abs(TDOMAIN(3)-probhiz).le.1.0E-10)) then
   ! do nothing
  else
   print *,"TDOMAIN invalid"
   stop
  endif

  print *,"reading the number of drops"
  read(2,*) T_num_drops
  print *,"number of drops, T_num_drops=",T_num_drops
  if (T_num_drops>0) then
   allocate(drop_data(T_num_drops,4))
   do nd=1,T_num_drops
    print *,"reading drop number nd=",nd
    read(2,*) drop_data(nd,1),drop_data(nd,2), &
          drop_data(nd,3),drop_data(nd,4)
    print *,"drop nd,x,y,z,r= ",nd,drop_data(nd,1), &
          drop_data(nd,2),drop_data(nd,3), &
          drop_data(nd,4)
   enddo ! nd=1,T_num_drops
  else
   print *,"T_num_drops invalid"
   stop
  endif

  close(2)

 else
  print *,"num_materials, probtype, or sdim invalid"
  stop
 endif

 return
 end subroutine INIT_TSPRAY_MODULE


  ! fluids tessellate the domain, solids are immersed. 
 subroutine TSPRAY_LS(x,t,LS)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 INTEGER_T nd
 REAL_T LS(num_materials)
 INTEGER_T in_droplet
 REAL_T minLS,tempLS

 if ((num_materials.eq.4).and. &
     (probtype.eq.402).and.  &
     (SDIM.eq.3)) then
  in_droplet=0
  minLS=1.0D+20
  tempLS=0.0d0
  do nd=1,T_num_drops
   tempLS=drop_data(nd,4)- &
       sqrt((x(1)-drop_data(nd,1))**2+ &
            (x(2)-drop_data(nd,2))**2+ &
            (x(SDIM)-drop_data(nd,SDIM))**2) 
   if (in_droplet.eq.0) then
    if (tempLS.ge.zero) then
     LS(1)=tempLS
     in_droplet=1
    else if (tempLS.le.zero) then
     if (abs(tempLS).lt.minLS) then
      minLS=abs(tempLS)
     endif
    else
     print *,"tempLS invalid"
     stop
    endif
   else if (in_droplet.eq.1) then
    ! do nothing
   else
    print *,"in_droplet invalid"
    stop
   endif
  enddo !nd=1,T_num_drops

  if (in_droplet.eq.0) then
   LS(1)=-minLS
  else if (in_droplet.eq.1) then
   ! do nothing
  else
   print *,"in_droplet invalid"
   stop
  endif
  LS(4)=zblob-x(SDIM)  ! substrate
  LS(3)=radblob2-abs(x(SDIM)-(zblob+radblob2)) ! ice layer
   ! gas
  tempLS=x(SDIM)-(zblob+two*radblob2)
  if (tempLS.le.zero) then
   LS(2)=tempLS
  else if (tempLS.ge.zero) then
   if (LS(1).ge.zero) then
    LS(2)=-LS(1)
   else if (LS(1).le.zero) then
    if (abs(LS(1)).le.abs(tempLS)) then
     LS(2)=-LS(1)
    else if (abs(LS(1)).ge.abs(tempLS)) then
     LS(2)=tempLS
    else
     print *,"LS(1) bust"
     stop
    endif
   else
    print *,"LS(1) bust"
    stop
   endif
  else
   print *,"tempLS bust"
   stop
  endif
 else
  print *,"num_materials, probtype, or sdim invalid"
  stop
 endif

 return
 end subroutine TSPRAY_LS

 subroutine TSPRAY_VEL(x,t,LS,VEL,velsolid_flag)
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
  print *,"adv_dir invalid in TSPRAY_VEL"
  stop
 endif

 return 
 end subroutine TSPRAY_VEL

 subroutine TSPRAY_PRES(x,t,LS,PRES)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T PRES

 PRES=outflow_pressure

 return 
 end subroutine TSPRAY_PRES


 subroutine TSPRAY_STATE(x,t,LS,STATE)
 use probcommon_module
 IMPLICIT NONE

 REAL_T x(SDIM)
 REAL_T t
 REAL_T LS(num_materials)
 REAL_T STATE(num_materials*num_state_material)
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
 end subroutine TSPRAY_STATE

  ! dir=1..sdim  side=1..2
 subroutine TSPRAY_LS_BC(xwall,xghost,t,LS, &
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
  call TSPRAY_LS(xghost,t,LS)
 else
  print *,"dir or side invalid"
  stop
 endif
 
 return
 end subroutine TSPRAY_LS_BC


  ! dir=1..sdim  side=1..2 veldir=1..sdim
 subroutine TSPRAY_VEL_BC(xwall,xghost,t,LS, &
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

  call TSPRAY_VEL(xghost,t,LS,local_VEL,velsolid_flag)
  VEL=local_VEL(veldir)

 else
  print *,"dir,side, or veldir invalid"
  stop
 endif

 return
 end subroutine TSPRAY_VEL_BC


  ! dir=1..sdim  side=1..2
 subroutine TSPRAY_PRES_BC(xwall,xghost,t,LS, &
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

  call TSPRAY_PRES(xghost,t,LS,PRES)

 else
  print *,"dir or side invalid"
  stop
 endif

 return
 end subroutine TSPRAY_PRES_BC

  ! dir=1..sdim  side=1..2
 subroutine TSPRAY_STATE_BC(xwall,xghost,t,LS, &
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
  call TSPRAY_STATE(xghost,t,LS,local_STATE)
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
 end subroutine TSPRAY_STATE_BC

 subroutine TSPRAY_HEATSOURCE(im,VFRAC,time,x,temp, &
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
 end subroutine TSPRAY_HEATSOURCE

end module TSPRAY_module
