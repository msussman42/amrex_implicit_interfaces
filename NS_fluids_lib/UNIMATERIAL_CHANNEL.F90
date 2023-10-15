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


 module unimaterialChannel_module
 use amrex_fort_module, only : amrex_real

  implicit none                   

  real(amrex_real) :: SIZE_L        ! Length of the channel
  real(amrex_real) :: SIZE_H        ! Height of the channel
  real(amrex_real) :: DENS_MAT      ! density
  real(amrex_real) :: TEMP_MAT      ! temperature
  real(amrex_real) :: PRES_MAT      ! pressure
  real(amrex_real) :: VELO_AVG      ! average inflow velocity  
  real(amrex_real) :: RAMP_TIME

 contains

   !called from PROB_CPP_PARMS.F90 (subroutine fort_override)
  subroutine UNIMAT_INIT_MODULE( &
   DENS_IN, &
   TEMP_IN, &
   PRES_IN, &
   VELO_IN, &
   RAMP_IN, &
   SIZE_L_IN, &
   SIZE_H_IN)
   IMPLICIT NONE

   real(amrex_real), INTENT(in) :: DENS_IN
   real(amrex_real), INTENT(in) :: TEMP_IN
   real(amrex_real), INTENT(in) :: PRES_IN
   real(amrex_real), INTENT(in) :: VELO_IN
   real(amrex_real), INTENT(in) :: RAMP_IN
   real(amrex_real), INTENT(in) :: SIZE_L_IN,SIZE_H_IN

   DENS_MAT=DENS_IN
   TEMP_MAT=TEMP_IN
   PRES_MAT=PRES_IN
   VELO_AVG=VELO_IN
   RAMP_TIME=RAMP_IN
   SIZE_L=SIZE_L_IN
   SIZE_H=SIZE_H_IN

   return
  end subroutine UNIMAT_INIT_MODULE

  !****************************************************
  ! level set initial value for material
  subroutine UNIMAT_INIT_LS_MAT(dx,LS)
   IMPLICIT NONE
   real(amrex_real) dx,LS

   LS = 10.0 * dx
  end subroutine UNIMAT_INIT_LS_MAT

  !****************************************************
  ! level set initial value for ghost material
  subroutine UNIMAT_INIT_LS_GST(dx,LS)
   IMPLICIT NONE
   real(amrex_real) dx,LS
   LS = -10.0 * dx
  end subroutine UNIMAT_INIT_LS_GST

  !****************************************************
  subroutine UNIMAT_INIT_VEL(x,y,z,vel)
   real(amrex_real), INTENT(in) :: x,y,z
   real(amrex_real), INTENT(out) :: vel(SDIM)

    if (RAMP_TIME.eq.zero) then
     vel(1)= VELO_AVG 
    else if (RAMP_TIME.gt.zero) then
     vel(1)=zero 
    else
     print *,"RAMP_TIME invalid"
     stop
    endif

    vel(2)= zero

  end subroutine UNIMAT_INIT_VEL

  ! Boundary condition
  ! Left: Inflow
  !---------------------
  ! Density         : Dirichlet  (*)
  ! Temperature     : Dirichlet  (*)
  ! Velocity        : Dirichlet  (*)
  ! Level set       : Dirichlet  (*) (unimaterial => switched to Extrap)
  ! Volume fraction : Dirichlet  (*) (unimaterial => switched to Extrap)
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

  !**************************************************** 
  subroutine UNIMAT_DENS_BC(&
   time, &                  ! time
   dir, &                   ! direction
   side, &                  ! side 
   q_out, &                 ! value outisde 
   x_wall, &                ! wall position
   q_in, &                  ! value inside
   x,y,z, &                 ! boundary point position
   dx, &                    ! dx
   im)                      ! material indicator

   integer dir,side,istate,im
   real(amrex_real) time,x_wall,x,y,z
   real(amrex_real) q_out,q_in
   real(amrex_real) dx(SDIM)

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(UNIMATERIAL_CHANNEL->DENS_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     q_out = DENS_MAT
    else if (side.eq.2) then
     q_out = q_in
    else 
     print *,"invalid side in x dirction &
      &(UNIMATERIAL_CHANNEL->DENS_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if ((side.eq.1).or.(side.eq.2)) then
     q_out = DENS_MAT
    else 
     print *,"invalid side in x dirction &
      &(UNIMATERIAL_CHANNEL->DENS_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (UNIMATERIAL_CHANNEL->DENS_BC)"
    stop
   endif
  end subroutine UNIMAT_DENS_BC

  !**************************************************** 
  subroutine UNIMAT_TEMP_BC(&
   time, &                  ! time
   dir, &                   ! direction
   side, &                  ! side 
   q_out, &                 ! value outisde 
   x_wall, &                ! wall position
   q_in, &                  ! value inside
   x,y,z, &                 ! boundary point position
   dx, &                    ! dx
   im)                      ! material indicator

   integer dir,side,im
   real(amrex_real) time,x_wall,x,y,z
   real(amrex_real) q_out,q_in
   real(amrex_real) dx(SDIM)

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(UNIMATERIAL_CHANNEL->TEMP_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     q_out = TEMP_MAT
    else if (side.eq.2) then
     q_out = q_in
    else 
     print *,"invalid side in x dirction &
      &(UNIMATERIAL_CHANNEL->TEMP_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if ((side.eq.1).or.(side.eq.2)) then
     q_out = TEMP_MAT
    else 
     print *,"invalid side in x dirction &
      &(UNIMATERIAL_CHANNEL->TEMP_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (UNIMATERIAL_CHANNEL->TEMP_BC)"
    stop
   endif
  end subroutine UNIMAT_TEMP_BC

  !**************************************************** 
  subroutine UNIMAT_VELO_BC(&
   time, &                  ! time
   dir, &                   ! direction
   side, &                  ! side 
   veldir, &
   q_out, &                 ! value outisde 
   q_in, &                  ! value inside
   x,y,z, &                 ! boundary point position
   dx)                      ! dx

   integer, INTENT(in) :: dir,side,veldir
   real(amrex_real), INTENT(in) :: time,x,y,z
   real(amrex_real), INTENT(out) :: q_out
   real(amrex_real), INTENT(in) :: q_in
   real(amrex_real), INTENT(in) :: dx(SDIM)
   real(amrex_real) rad,R


   if (SDIM.ne.2) then
    print *,"invalid system dimension (UNIMATERIAL_CHANNEL->VELO_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then ! xlo
     if (veldir.eq.1) then

      if (RAMP_TIME.eq.zero) then
       q_out = VELO_AVG 
      else if (RAMP_TIME.gt.zero) then
       if (time.ge.RAMP_TIME) then
        q_out = VELO_AVG 
       else if ((time.ge.zero).and.(time.le.RAMP_TIME)) then
        q_out=VELO_AVG*time/RAMP_TIME
       else
        print *,"time invalid"
        stop
       endif
      else
       print *,"RAMP_TIME invalid"
       stop
      endif

     else if (veldir.eq.2) then 
      q_out = zero  
    else if (side.eq.2) then ! xhi
      q_out = q_in  
     else 
      print *,"veldir invalid"
      stop
     endif
    else 
     print *,"invalid side in x dirction (UNIMATERIAL_CHANNEL->VELO_BC)"
     stop
    endif
   else if (dir.eq.2) then  ! SLIP WALL
    if (side.eq.1) then
     if (veldir.eq.1) then

      if (RAMP_TIME.eq.zero) then
       q_out = VELO_AVG 
      else if (RAMP_TIME.gt.zero) then
       if (time.ge.RAMP_TIME) then
        q_out = VELO_AVG 
       else if ((time.ge.zero).and.(time.le.RAMP_TIME)) then
        q_out=VELO_AVG*time/RAMP_TIME
       else
        print *,"time invalid"
        stop
       endif
      else
       print *,"RAMP_TIME invalid"
       stop
      endif

     else if (veldir.eq.2) then
      q_out = zero
     else 
      print *,"veldir invalid"
      stop
     endif
    else if (side.eq.2) then
     if (veldir.eq.1) then

      if (RAMP_TIME.eq.zero) then
       q_out = VELO_AVG 
      else if (RAMP_TIME.gt.zero) then
       if (time.ge.RAMP_TIME) then
        q_out = VELO_AVG 
       else if ((time.ge.zero).and.(time.le.RAMP_TIME)) then
        q_out=VELO_AVG*time/RAMP_TIME
       else
        print *,"time invalid"
        stop
       endif
      else
       print *,"RAMP_TIME invalid"
       stop
      endif

     else if (veldir.eq.2) then
      q_out = zero
     else 
      print *,"veldir invalid"
      stop
     endif
    else 
     print *,"invalid side in x dirction (UNIMATERIAL_CHANNEL->VELO_BC)"
     stop
    endif
   else
    print *,"invalid direction (UNIMATERIAL_CHANNEL->VELO_BC)"
    stop
   endif
  end subroutine UNIMAT_VELO_BC
  !**************************************************** 
  subroutine UNIMAT_LVLS_BC(&
   time, &                  ! time
   dir, &                   ! direction
   side, &                  ! side 
   q_out, &                 ! value outisde 
   x_wall, &                ! wall position
   q_in, &                  ! value inside
   x,y,z, &                 ! boundary point position
   dx, &                    ! dx
   im)                      ! material indicator

   integer dir,side,istate,im
   real(amrex_real) time,x_wall,x,y,z
   real(amrex_real) q_out,q_in
   real(amrex_real) dx(SDIM)

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(UNIMATERIAL_CHANNEL->LVLS_BC)"
    stop
   endif
   
   if(((dir.eq.1).or.(dir.eq.2)).and.((side.eq.1).or.(side.eq.2))) then 
    q_out = q_in
   else 
    print *,"invalid dir/side in &
     &(UNIMATERIAL_CHANNEL->LVLS_BC)"
    stop
   endif
  end subroutine UNIMAT_LVLS_BC

  !**************************************************** 
  subroutine UNIMAT_VOLF_BC(&
   time, &                  ! time
   dir, &                   ! direction
   side, &                  ! side 
   q_out, &                 ! value outisde 
   x_wall, &                ! wall position
   q_in, &                  ! value inside
   x,y,z, &                 ! boundary point position
   dx, &                    ! dx
   im)                      ! material indicator

   integer dir,side,istate,im
   real(amrex_real) time,x_wall,x,y,z
   real(amrex_real) q_out,q_in
   real(amrex_real) dx(SDIM)

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(UNIMATERIAL_CHANNEL->VOLF_BC)"
    stop
   endif

   if(((dir.eq.1).or.(dir.eq.2)).and.((side.eq.1).or.(side.eq.2))) then
    q_out = q_in
   else 
    print *,"invalid dir/side&
     &(UNIMATERIAL_CHANNEL->VOLF_BC)"
    stop
   endif
  end subroutine UNIMAT_VOLF_BC

  !**************************************************** 
  subroutine UNIMAT_MOMF_BC(&
   time, &                  ! time
   dir, &                   ! direction
   side, &                  ! side 
   q_out, &                 ! value outisde 
   x_wall, &                ! wall position
   q_in, &                  ! value inside
   x,y,z, &                 ! boundary point position
   dx, &                    ! dx
   im)                      ! material indicator

   integer dir,side,istate,im
   real(amrex_real) time,x_wall,x,y,z
   real(amrex_real) q_out(SDIM), q_in(SDIM)
   real(amrex_real) dx(SDIM)

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(UNIMATERIAL_CHANNEL->MOMF_BC)"
    stop
   endif

   if(dir.eq.1) then
    if ((side.eq.1).or.(side.eq.2)) then
     q_out(1)=-q_in(1)
     q_out(2)= q_in(2)
    else 
     print *,"invalid side in x dirction &
      &(UNIMATERIAL_CHANNEL->MOMF_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if ((side.eq.1).or.(side.eq.2)) then
     q_out(1)= q_in(1)
     q_out(2)=-q_in(2)
    else 
     print *,"invalid side in x dirction &
      &(UNIMATERIAL_CHANNEL->MOMF_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (UNIMATERIAL_CHANNEL->MOMF_BC)"
    stop
   endif
  end subroutine UNIMAT_MOMF_BC

  !**************************************************** 
  subroutine UNIMAT_PRES_BC(&
   time, &                  ! time
   dir, &                   ! direction
   side, &                  ! side 
   q_out, &                 ! value outisde 
   x_wall, &                ! wall position
   q_in, &                  ! value inside
   x,y,z, &                 ! boundary point position
   dx)                      ! dx

   integer dir,side,istate
   real(amrex_real) time,x_wall,x,y,z
   real(amrex_real) q_out, q_in
   real(amrex_real) dx(SDIM)

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(UNIMATERIAL_CHANNEL->PRES_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     !! do nothing
     !! Nuemann BC this should not be called
     print *,"This should not be called &
      &(UNIMATERIAL_CHANNEL->PRES_BC)"
     stop
    else if (side.eq.2) then
     !! It should not matter at this version of the code
     q_out=0
    else 
     print *,"invalid side in x dirction &
      &(UNIMATERIAL_CHANNEL->PRES_BC)"
     stop
    endif
   else if (dir.eq.2) then
    !! do nothing
    !! Nuemann BC this should not be called
     print *,"This should not be called &
      &(UNIMATERIAL_CHANNEL->PRES_BC)"
     stop
!!$    if (side.eq.1) then
!!$     q_out= q_in
!!$    else if (side.eq.2) then
!!$     q_out=q_in
!!$    else 
!!$     print *,"invalid side in x dirction &
!!$      &(UNIMATERIAL_CHANNEL->PRES_BC)"
!!$     stop
!!$    endif
    else
     print *,"invalid direction&
      & (UNIMATERIAL_CHANNEL->PRES_BC)"
     stop
    endif
  end subroutine UNIMAT_PRES_BC

  !**************************************************** 
 end module unimaterialChannel_module
