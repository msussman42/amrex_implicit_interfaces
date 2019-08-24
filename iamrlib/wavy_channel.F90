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

  ! probtype==915
 module WAVY_Channel_module

  implicit none                   

 contains

  subroutine INIT_WAVY_MODULE() 
   IMPLICIT NONE
   return
  end subroutine INIT_WAVY_MODULE

  !****************************************************
  ! level set initial value for vapor material
  subroutine WAVY_INIT_LS(x,t,LS)
   use probcommon_module
   IMPLICIT NONE

   REAL_T x(SDIM)
   REAL_T t
   INTEGER_T im
   REAL_T LS(num_materials)

    ! 12 x 12 domain
   if (SDIM.eq.2) then
    if ((num_materials.eq.3).and.(probtype.eq.915)) then
     do im=1,num_materials 
      if (im.eq.1) then !vapor
       LS(im)=x(SDIM)-7.0
      else if (im.eq.2) then ! liquid
       LS(im)=7.0 - x(SDIM)
      else if (im.eq.3) then ! wavy substrate
         ! does not have to be an exact distance.  Main thing
         ! is that the function is smooth and the zero LS is correct.
         ! also, since RZ, use cos(x(1) pi/3) since this is symmetric
         ! at r=0.
       LS(im)=three+sin(two*Pi*t)+cos(x(1)*Pi/three)-x(SDIM)
      else
       print *,"im invalid"
       stop
      endif
     enddo ! im=1..num_materials
    else
     print *,"num_materials or probtype invalid"
     stop
    endif
   else
    print *,"not setup for this dimension yet dim=",SDIM
    stop
   endif

   return
  end subroutine WAVY_INIT_LS


  subroutine WAVY_INIT_STATE(x,t,LS,STATE)
   use probcommon_module
   IMPLICIT NONE

   REAL_T x(SDIM)
   REAL_T t
   REAL_T LS(num_materials)
   REAL_T STATE(num_materials*num_state_material)
   INTEGER_T im,ibase

   if ((num_materials.eq.3).and. &
       (num_state_material.eq.2).and. &
       (num_species_var.eq.0).and. &
       (probtype.eq.915)) then
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
    enddo ! im=1..num_materials
   else
    print *,"num_materials,num_state_material, num_species_var,"
    print *,"or probtype invalid"
    stop
   endif

   return
  end subroutine WAVY_INIT_STATE

  !********************************************************   

  subroutine WAVY_INIT_VEL(x,t,LS,VEL,velsolid_flag)
   use probcommon_module
   IMPLICIT NONE

   REAL_T x(SDIM)
   REAL_T t
   REAL_T LS(num_materials)
   REAL_T VEL(SDIM)
   INTEGER_T dir,im_solid_wavy
   INTEGER_T velsolid_flag

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
   im_solid_wavy=3
   if ((LS(im_solid_wavy).ge.zero).or. &
       (velsolid_flag.eq.1)) then
    VEL(SDIM)=two*Pi*cos(two*Pi*t)
   else if ((LS(im_solid_wavy).le.zero).and. &
            (velsolid_flag.eq.0)) then
    ! do nothing
   else
    print *,"LS(im_solid_wavy) bust"
    stop
   endif

   return 
  end subroutine WAVY_INIT_VEL

  !********************************************************  

  subroutine WAVY_INIT_PRES(x,t,LS,PRES)
   use probcommon_module
   IMPLICIT NONE

   REAL_T x(SDIM)
   REAL_T t
   REAL_T LS(num_materials)
   REAL_T PRES
   REAL_T gravity_dz

   if (SDIM.eq.2) then
    gravity_dz=x(SDIM)-probhiy
   else if (SDIM.eq.3) then
    gravity_dz=x(SDIM)-probhiz
   else
    print *,"dimension bust"
    stop
   endif

   PRES=-fort_denconst(1)*abs(gravity)*gravity_dz

   return 
  end subroutine WAVY_INIT_PRES

  !******************************************************* 

  ! Boundary condition
  ! Left: Inflow
  !---------------------
  ! Density         : Dirichlet  (*)
  ! Temperature     : Dirichlet  (*)
  ! Velocity        : Dirichlet  (*)
  ! Level set       : Dirichlet  (*) (WAVY => switched to Extrap)
  ! Volume fraction : Dirichlet  (*) (WAVY => switched to Extrap)
  ! Pressure        : N  (*)

  ! Right: Outflow
  !---------------------
  ! Density         : Extrapolation (*)
  ! Temperature     : Extrapolation (*)
  ! Velocity        : Extrapolation (*)
  ! Level set       : Extrapolation (*)
  ! Volume fraction : Extrapolation (*)
  ! Pressure        : D (*) ?

  ! Top and bottom: Inflow
  !---------------------
  ! Density         : Dirichlet  (*)
  ! Temperature     : Dirichlet  (*) 
  ! Velocity        : Dirichlet  (*)
  ! Level set       : Dirichlet  (*)
  ! Volume fraction : Dirichlet  (*)
  ! Pressure        : N (*) Dp/dy=0

  ! (*) : user function needed
  ! Extrapolation are low order (?)

  !******************************************************* 

  subroutine WAVY_LS_BC(xwall,xghost,t,LS, &
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
    call WAVY_INIT_LS(xghost,t,LS)
   else
    print *,"dir or side invalid"
    stop
   endif

   return
  end subroutine WAVY_LS_BC

  !******************************************************* 

  subroutine WAVY_VEL_BC(xwall,xghost,t,LS, &
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

    call WAVY_INIT_VEL(xghost,t,LS,local_VEL,velsolid_flag)
    VEL=local_VEL(veldir)

   else
    print *,"dir,side, or veldir invalid"
    stop
   endif

   return
  end subroutine WAVY_VEL_BC

  !******************************************************* 

  subroutine WAVY_PRES_BC(xwall,xghost,t,LS, &
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

    call WAVY_INIT_PRES(xghost,t,LS,PRES)

   else
    print *,"dir or side invalid"
    stop
   endif

   return
  end subroutine WAVY_PRES_BC

  !******************************************************* 

  subroutine WAVY_STATE_BC(xwall,xghost,t,LS, &
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
    call WAVY_INIT_STATE(xghost,t,LS,local_STATE)
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
  end subroutine WAVY_STATE_BC

 end module WAVY_Channel_module

