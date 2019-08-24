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


 !    __________
 !   |          |            
 !   |   Air    |            
 !   |_____     |            
 !   |     \    |            
 !   |__Ice |   |           ^
 !   |Ai\   |   |           |
 !   |__/   |   |           Z_ICE  ^ 
 !   |_____/    |           |      |
 !   |          |           |      Z_AIR
 !   |__________|           |      |
 !   |   Water  | ^         |      |
 !   |          | Z_WATER   |      |
 !   |__________| v         v      v
 !
 !   <-R_ICE->
 !   <R_A>
 !    

 module snow_module

  implicit none                   

  REAL_T :: Z_WATER      = 4.00   ! height of water (yblob)
  REAL_T :: Z_ICE        = 6.00   ! center of snow sphere (yblob2)
  REAL_T :: Z_AIR        = 6.00   ! center of air bubble (yblob3)

  REAL_T :: R_ICE        = 1.00   ! radius of ice sphere (radblob)
  REAL_T :: R_AIR        = 0.50   ! radius of air bubble (radblob2)

  REAL_T :: TEMP_WATER   = 277    ! initial water temperature 
  REAL_T :: TEMP_ICE     = 273    ! initial ice temperature 
  REAL_T :: TEMP_AIR     = 277    ! initial air temperature 

  REAL_T :: VEL_SNOW     = 1.000  ! snow (ice and air inside) particle velocity (radblob3)
  REAL_T :: PRESSURE_AIR = 101.325! air pressure 
 contains

  subroutine INIT_SNOW_MODULE( &
   Z_WATER_IN, &
   Z_ICE_IN, &
   Z_AIR_IN, &
   R_ICE_IN, &
   R_AIR_IN, &
   TEMP_WATER_IN, &
   TEMP_ICE_IN, &
   TEMP_AIR_IN, &
   VEL_SNOW_IN, &
   PRESSURE_AIR_IN)

   IMPLICIT NONE
   REAL_T Z_WATER_IN
   REAL_T Z_ICE_IN
   REAL_T Z_AIR_IN
   REAL_T R_ICE_IN
   REAL_T R_AIR_IN
   REAL_T TEMP_WATER_IN
   REAL_T TEMP_ICE_IN
   REAL_T TEMP_AIR_IN
   REAL_T VEL_SNOW_IN
   REAL_T PRESSURE_AIR_IN

   Z_WATER=Z_WATER_IN
   Z_ICE=Z_ICE_IN
   Z_AIR=Z_AIR_IN
   R_ICE=R_ICE_IN
   R_AIR=R_AIR_IN
   TEMP_WATER=TEMP_WATER_IN
   TEMP_ICE=TEMP_ICE_IN
   TEMP_AIR=TEMP_AIR_IN
   VEL_SNOW=VEL_SNOW_IN
   PRESSURE_AIR=PRESSURE_AIR_IN

   return
  end subroutine INIT_SNOW_MODULE

  !****************************************************
  ! level set initial value for snow porblme
  subroutine INIT_LS_SNOW(x,y,z,t,LS_W,LS_I,LS_A)
   IMPLICIT NONE
   REAL_T x,y,z,t
   REAL_T LS_W,LS_I,LS_A
   REAL_T temp_ls_ice, temp_ls_air
   if(z.le.Z_WATER) then
    ! inside water region
    LS_W = Z_WATER - z
    LS_I = CIRCLE_LS(Z_ICE,R_ICE,x,z) 
    LS_A = z - Z_WATER
   else
    temp_ls_ice = CIRCLE_LS(Z_ICE,R_ICE,x,z)
    if(temp_ls_ice.le.zero) then
     ! in outise air
     LS_W = Z_WATER - z
     LS_I = temp_ls_ice
     LS_A = min(-LS_W,-LS_I)
    else
     ! in snow region
     temp_ls_air = CIRCLE_LS(Z_AIR,R_AIR,x,z)
     if(temp_ls_ice.le.zero) then
      ! in ice
      LS_W = Z_WATER - z
      LS_I = min(temp_ls_ice, -temp_ls_air)
      LS_A = -LS_I
     else
      ! in air bubble in ice
      LS_W = Z_WATER - z
      LS_I = -temp_ls_air
      LS_A = temp_ls_air
     end if ! in ice
    end if ! in outise air
   end if ! in water 
  end subroutine INIT_LS_SNOW

  !****************************************************
  ! velocity initial value for snow porblme
  subroutine SNOW_INIT_VEL(x,y,z,velcell)
   REAL_T x,y,z
   REAL_T velcell(SDIM)
   REAL_T ls_snow

   velcell(1)= zero
   velcell(2)= zero
   velcell(SDIM)= zero

   ls_snow = CIRCLE_LS(Z_ICE,R_ICE,x,z)

   if(ls_snow.ge.zero) then
    velcell(2)= VEL_SNOW
   end if

  end subroutine SNOW_INIT_VEL

  !****************************************************
  ! Evalute the levelset function for a given circle cordinate and
  ! radius and a given point in RZ axisymmetric coordinates. Center
  ! of the circle is located on r=0
  ! Input : circle_z  : Circle center z component
  !       : circle_rad: Circel radius
  !       : r         : point r componet
  !       : z         : point z componet
  ! Output: circle_ls : levelset respect to circle circumfrance
  !                     (positive inside, negative outise)

  REAL_T function CIRCLE_LS(circle_z,circle_rad,r,z) 
   REAL_T circle_z,circle_rad,r,z
   
   circle_ls = circle_rad - sqrt((circle_z-z)**2+r**2)

   return
  end function CIRCLE_LS
  !**************************************************** 

  ! Boundary condition
  ! Left : symmetry
  !---------------------
  ! Density         : Reflect even
  ! Temperature     : Reflect even
  ! Velocity        : U_n  rfelcet odd, U_t Reflect even
  ! Level set       : Reflect even
  ! Volume fraction : Reflect even
  ! Pressure        : Reflect even

  ! Right : No slip wall + Adibatic temperature
  !---------------------
  ! Density         : Extrapolation (*)
  ! Temperature     : Reflect even
  ! Velocity        : Dirichlet U_n=U_t=0 (*)
  ! Level set       : Extrapolation (*)
  ! Volume fraction : Extrapolation (*)
  ! Pressure        : Reflect even

  ! Bottom : Noslip wall + fixed temperature
  !---------------------
  ! Density         : Extrapolation (*)
  ! Temperature     : Dirichlet  (*)
  ! Velocity        : Dirichlet U_n=U_t=0 (*)
  ! Level set       : Extrapolation (*)
  ! Volume fraction : Extrapolation (*)
  ! Pressure        : Reflect even

  ! Top : outflow
  !---------------------
  ! Density         : Extrapolation (*)
  ! Temperature     : Extrapolation (*)
  ! Velocity        : Extrapolation (*)
  ! Level set       : Extrapolation (*)
  ! Volume fraction : Extrapolation (*)
  ! Pressure        : Dirichlet (*)

  ! (*) : user function needed
  ! Extrapolation are low order (?)

  !**************************************************** 
  subroutine SNOW_DENS_BC(&
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
     &(HYDRATE_REACTOR->DENS_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     print *,"This should not be called, &
      &lo side in direction 1 is symmetry BC &
      &(SNOW->DENS_BC)"
     stop
    else if (side.eq.2) then
     q_out = q_in
    else 
     print *,"invalid side in dirction 1 &
      &(SNOW->DENS_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if ((side.eq.1).or.(side.eq.2)) then
     q_out = q_in
    else 
     print *,"invalid side in dirction 2 &
      &(SNOW->DENS_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (SNOW->DENS_BC)"
    stop
   endif
  end subroutine SNOW_DENS_BC

  !**************************************************** 
  subroutine SNOW_TEMP_BC(&
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
     &(SNOW->TEMP_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     print *,"This should not be called, &
      &lo side in direction 1 is symmetry BC &
      &(SNOW->TEMP_BC)"
     stop
    else if (side.eq.2) then
     q_out = q_in
    else 
     print *,"invalid side in dirction 1 &
      &(SNOW->TEMP_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if (side.eq.1) then
     q_out = TEMP_WATER
    else if (side.eq.2) then
     q_out = TEMP_AIR
    else
     print *,"invalid side in dirction 2 &
      &(SNOW->TEMP_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (SNOW->TEMP_BC)"
    stop
   endif
  end subroutine SNOW_TEMP_BC

  !**************************************************** 
  subroutine SNOW_VELO_BC(&
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

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(SNOW->VELO_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     print *,"This should not be called, &
      &lo side in direction 1 is symmetry BC &
      &(SNOW->VELO_BC)"
     stop
    else if (side.eq.2) then ! xhi
     if (veldir.eq.1) then
      q_out = zero
     else if (veldir.eq.2) then ! xhi, y velocity
      q_out = zero  
     else 
      print *,"veldir invalid"
      stop
     endif
    else 
     print *,"invalid side in dirction 1 &
      &(SNOW->VELO_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if (side.eq.1) then
     if (veldir.eq.1) then
      q_out = zero
     else if (veldir.eq.2) then
      q_out = zero
     else 
      print *,"veldir invalid"
      stop
     endif
    else if (side.eq.2) then
     if (veldir.eq.1) then
      q_out = q_in ! q_in
     else if (veldir.eq.2) then
      q_out = q_in ! q_in
     else 
      print *,"veldir invalid"
      stop
     endif
    else 
     print *,"invalid side in dirction 2 &
      &(SNOW->VELO_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (SNOW->VELO_BC)"
    stop
   endif
  end subroutine SNOW_VELO_BC
  !**************************************************** 
  subroutine SNOW_LVLS_BC(&
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
     &(SNOW->LVLS_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     print *,"This should not be called, &
      &lo side in direction 1 is symmetry BC &
      &(SNOW->LVLS_BC)"
     stop
    else if (side.eq.2) then
     q_out = q_in
    else 
     print *,"invalid side in dirction 1 &
      &(SNOW->LVLS_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if ((side.eq.1).or.(side.eq.2)) then
     q_out = q_in
    else 
     print *,"invalid side in dirction 2 &
      &(SNOW->LVLS_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (SNOW->LVLS_BC)"
    stop
   endif
  end subroutine SNOW_LVLS_BC
  
  !**************************************************** 
  subroutine SNOW_VOLF_BC(&
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
     &(SNOW->VOLF_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     print *,"This should not be called, &
      &lo side in direction 1 is symmetry BC &
      &(SNOW->VOLF_BC)"
     stop
    else if (side.eq.2) then
     q_out = q_in
    else 
     print *,"invalid side in dirction 1 &
      &(SNOW->VOLF_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if ((side.eq.1).or.(side.eq.2)) then
     q_out = q_in
    else 
     print *,"invalid side in dirction 2 &
      &(SNOW->VOLF_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (SNOW->VOLF_BC)"
    stop
   endif
  end subroutine SNOW_VOLF_BC

  !**************************************************** 
  subroutine SNOW_MOMF_BC(&
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
   REAL_T q_out(SDIM), q_in(SDIM)
   REAL_T dx(SDIM)

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(SNOW->MOMF_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     print *,"This should not be called, &
      &lo side in direction 1 is symmetry BC &
      &(SNOW->MOMF_BC)"
     stop
    else if (side.eq.2) then
     q_out(1)=-q_in(1)
     q_out(2)= q_in(2)
    else 
     print *,"invalid side in dirction 1 &
      &(SNOW->MOMF_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if ((side.eq.1).or.(side.eq.2)) then
     q_out(1)= q_in(1)
     q_out(2)=-q_in(2)
    else 
     print *,"invalid side in dirction 2 &
      &(SNOW->MOMF_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (SNOW->MOMF_BC)"
    stop
   endif
  end subroutine SNOW_MOMF_BC

  !**************************************************** 
  subroutine SNOW_PRES_BC(&
   time, &                  ! time
   dir, &                   ! direction
   side, &                  ! side 
   q_out, &                 ! value outisde 
   x_wall, &                ! wall position
   q_in, &                  ! value inside
   x,y,z, &                 ! boundary point position
   dx)                      ! dx

   INTEGER_T dir,side
   REAL_T time,x_wall,x,y,z
   REAL_T q_out, q_in
   REAL_T dx(SDIM)

   if (SDIM.ne.2) then
    print *,"invalid system dimension &
     &(SNOW->PRES_BC)"
    stop
   endif

   if(dir.eq.1) then
    if (side.eq.1) then
     print *,"This should not be called, &
      &lo side in direction 1 is symmetry BC &
      &(SNOW->PRES_BC)"
     stop
    else if (side.eq.2) then
     q_out=q_in
    else 
     print *,"invalid side in x dirction &
      &(SNOW->PRES_BC)"
     stop
    endif
   else if (dir.eq.2) then
    if (side.eq.1) then
     q_out= q_in
    else if (side.eq.2) then
     q_out=PRESSURE_AIR
    else 
     print *,"invalid side in x dirction &
      &(SNOW->PRES_BC)"
     stop
    endif
   else
    print *,"invalid direction&
     & (SNOW->PRES_BC)"
    stop
   endif
  end subroutine SNOW_PRES_BC

  !**************************************************** 
 end module snow_module
