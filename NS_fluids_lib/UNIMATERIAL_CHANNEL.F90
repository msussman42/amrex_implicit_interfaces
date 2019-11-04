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

  implicit none                   

  REAL_T :: SIZE_L        ! Length of the channel
  REAL_T :: SIZE_H        ! Height of the channel
  REAL_T :: DENS_MAT      ! density
  REAL_T :: TEMP_MAT      ! temperature
  REAL_T :: PRES_MAT      ! pressure
  REAL_T :: VELO_AVG      ! average inflow velocity  

 contains

  subroutine UNIMAT_INIT_MODULE( &
   DENS_IN, &
   TEMP_IN, &
   PRES_IN, &
   VELO_IN, &
   SIZE_L_IN, &
   SIZE_H_IN)
   IMPLICIT NONE

   REAL_T DENS_IN
   REAL_T TEMP_IN
   REAL_T PRES_IN
   REAL_T VELO_IN
   REAL_T SIZE_L_IN,SIZE_H_IN

   DENS_MAT=DENS_IN
   TEMP_MAT=TEMP_IN
   PRES_MAT=PRES_IN
   VELO_AVG=VELO_IN
   SIZE_L=SIZE_L_IN
   SIZE_H=SIZE_H_IN

   return
  end subroutine UNIMAT_INIT_MODULE

  !****************************************************
  ! level set initial value for material
  subroutine UNIMAT_INIT_LS_MAT(dx,LS)
   IMPLICIT NONE
   REAL_T dx,LS

   LS = 10.0 * dx
  end subroutine UNIMAT_INIT_LS_MAT

  !****************************************************
  ! level set initial value for ghost material
  subroutine UNIMAT_INIT_LS_GST(dx,LS)
   IMPLICIT NONE
   REAL_T dx,LS
   LS = -10.0 * dx
  end subroutine UNIMAT_INIT_LS_GST

  !****************************************************
  subroutine UNIMAT_INIT_VEL(x,y,z,vel)
   REAL_T x,y,z,v_avg
   REAL_T vel(SDIM)
   
    vel(1)= VELO_AVG 
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

   INTEGER_T dir,side,istate,im
   REAL_T time,x_wall,x,y,z
   REAL_T q_out,q_in
   REAL_T dx(SDIM)

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

   INTEGER_T dir,side,im
   REAL_T time,x_wall,x,y,z
   REAL_T q_out,q_in
   REAL_T dx(SDIM)

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

   INTEGER_T dir,side,veldir
   REAL_T time,x,y,z
   REAL_T q_out,q_in
   REAL_T dx(SDIM)
   REAL_T rad,R


   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(UNIMATERIAL_CHANNEL->VELO_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then ! xlo
     if (veldir.eq.1) then
      q_out = VELO_AVG 
     else if (veldir.eq.2) then 
      q_out = zero  
    else if (side.eq.2) then ! xhi
      q_out = q_in  
     else 
      print *,"veldir invalid"
      stop
     endif
    else 
     print *,"invalid side in x dirction &
      &(UNIMATERIAL_CHANNEL->VELO_BC)"
     stop
    endif
   else if (dir.eq.2) then  ! SLIP WALL
    if (side.eq.1) then
     if (veldir.eq.1) then
      q_out = VELO_AVG ! q_in
     else if (veldir.eq.2) then
      q_out = zero
     else 
      print *,"veldir invalid"
      stop
     endif
    else if (side.eq.2) then
     if (veldir.eq.1) then
      q_out = VELO_AVG ! q_in
     else if (veldir.eq.2) then
      q_out = zero
     else 
      print *,"veldir invalid"
      stop
     endif
    else 
     print *,"invalid side in x dirction &
      &(HYDRATE_REACTOR->VELO_BC)"
     stop
    endif
   else
    print *,"invalid direction &
     & (UNIMATERIAL_CHANNEL->VELO_BC)"
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

   INTEGER_T dir,side,istate,im
   REAL_T time,x_wall,x,y,z
   REAL_T q_out,q_in
   REAL_T dx(SDIM)

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

   INTEGER_T dir,side,istate,im
   REAL_T time,x_wall,x,y,z
   REAL_T q_out,q_in
   REAL_T dx(SDIM)

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

   INTEGER_T dir,side,istate,im
   REAL_T time,x_wall,x,y,z
   REAL_T q_out(SDIM), q_in(SDIM)
   REAL_T dx(SDIM)

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

   INTEGER_T dir,side,istate
   REAL_T time,x_wall,x,y,z
   REAL_T q_out, q_in
   REAL_T dx(SDIM)

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
