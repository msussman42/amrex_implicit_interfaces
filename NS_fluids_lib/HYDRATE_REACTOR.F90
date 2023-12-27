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

! From the experiment done by (JammaludinETAL1991 - Hydrate plugging
! problems in undersea natural gas pipelines under shutdown
! conditions). The reactor is a cylinder with internal volume 500
! cm^3.  (Note: The schematic figure in the paper results in a larger
! volume (629.3307 cm^3), however the instruments and equipments may
! have reduced the effective volume). Based on these assumptions we
! have:

!    __________
!   |          |           ^
!   |   Gas    |           |
!   |          |           |
!   |__________|           |
!   |          |           |
!   | Hydrate  |           H
!   |__________|           |
!   |          | ^         |
!   |          | |         |
!   |   Water  | H_water   |
!   |          | |         |
!   |__________| v         v
!
!   <--- R ---->
!
!    R = 3.81  =>   A = 45.6037
!                   H = 10.9640
!                   H_water = 6.5784
!

! If doing a full domain simulation, to make the domain size ratio an
! integer number, we pick
!
!      R = 3.81 , H = 3*R = 11.43

module hydrateReactor_module
use amrex_fort_module, only : amrex_real

implicit none                   

real(amrex_real) :: SIZE_H       = 11.43  ! Height of the reactor 
real(amrex_real) :: SIZE_R       = 3.81   ! Radius of the reactor


real(amrex_real) :: H_WATER      = 6.5780 ! initial water level [cm]
real(amrex_real) :: H_HYDRATE    = 7      ! initial hydrate level [cm]
                                ! (TEMPORARY!  THIS SHOULD BE PICKED
                                ! BASED ON THE EXPERIMENT INITIAL
                                ! RESULT)

real(amrex_real) :: TEMP_WATER   = 274    ! initial water temperature [K]
real(amrex_real) :: TEMP_HYDRATE = 274    ! initial hydrate temperature [K]
real(amrex_real) :: TEMP_GAS     = 274    ! initial methane temperature [K]

real(amrex_real) :: DENS_WATER   = 1      ! initial water density [g.cm^-3]
real(amrex_real) :: DENS_HYDRATE = 0.91   ! initial hydrate density [g.cm^-3]
                                ! (Source: SloanKoh2008 p. 269)
real(amrex_real) :: DENS_GAS     = 0.0352 ! initial gas density [g.cm^-3]

! Evaluated by the ideal gas law based on the prescribed pressure and
! volume of the gas phase in the reactor: 
! H_HYDRATE = 7 cm, P = 4 MPa => rho_gas = 0.0282 g.cm^-3
! H_HYDRATE = 7 cm, P = 5 MPa => rho_gas = 0.0352 g.cm^-3
! H_HYDRATE = 7 cm, P = 7 MPa => rho_gas = 0.0493 g.cm^-3

! MATLAB CODE:
! H = 11.43; % cm
! H_h = 7; % cm
! R_r = 3.81; % cm
! P = 7; % MPa
! T = 274; % K
! R = 8.314459848 ; %cm^3 MPa K^−1 mol^−1
! M = 16.0425 ; % g mol^-1
!
! V = (H-H_h)*pi*R_r^2; % cm^3
! n = (P*V)/(R*T) ; % mol
! m = n * M; % g
! rho = m/V; % g cm^-3


real(amrex_real) :: PRESSURE_GAS = 5.0e7           ! gas pressure during the
                                         ! experiment [Ba] (1Pa=10Ba)

real(amrex_real) :: MOLAR_MASS_GAS = 16.0425        ! gas molar mass [g.mol^-1]
real(amrex_real) :: MOLAR_MASS_HYDRATE = 119.630475 ! gas molar mass [g.mol^-1]

contains

 subroutine INIT_HYDRATE_MODULE(TEMP_WATER_IN,TEMP_HYDRATE_IN, &
  TEMP_GAS_IN,DENS_WATER_IN,DENS_HYDRATE_IN,DENS_GAS_IN, &
  PRESSURE_GAS_IN,SIZE_H_IN,SIZE_R_IN, &
  H_WATER_IN,H_HYDRATE_IN)
 IMPLICIT NONE

 real(amrex_real) TEMP_WATER_IN,TEMP_HYDRATE_IN,TEMP_GAS_IN
 real(amrex_real) DENS_WATER_IN,DENS_HYDRATE_IN,DENS_GAS_IN
 real(amrex_real) PRESSURE_GAS_IN
 real(amrex_real) SIZE_H_IN,SIZE_R_IN
 real(amrex_real) H_WATER_IN,H_HYDRATE_IN

 TEMP_WATER=TEMP_WATER_IN
 TEMP_HYDRATE=TEMP_HYDRATE_IN
 TEMP_GAS=TEMP_GAS_IN
 DENS_WATER=DENS_WATER_IN
 DENS_HYDRATE=DENS_HYDRATE_IN
 DENS_GAS=DENS_GAS_IN

 PRESSURE_GAS=PRESSURE_GAS_IN
 SIZE_H=SIZE_H_IN
 SIZE_R=SIZE_R_IN
 H_WATER=H_WATER_IN
 H_HYDRATE=H_HYDRATE_IN

 return
 end subroutine INIT_HYDRATE_MODULE
  
!****************************************************
! level set initial value for water
  subroutine INIT_LS_WATER(x,y,z,t,LS)
    IMPLICIT NONE
    real(amrex_real) x,y,z,t,LS
    
    LS = H_WATER - z
  end subroutine INIT_LS_WATER

!****************************************************
! level set initial value for hydrate
  subroutine INIT_LS_HYDRATE(x,y,z,t,LS)
    IMPLICIT NONE
    real(amrex_real) x,y,z,t,LS
    real(amrex_real) mid

    mid = (H_WATER + H_HYDRATE)/two
    if(z.le.mid) then
       LS = z - H_WATER
    else
       LS = H_HYDRATE - z
    endif
  end subroutine INIT_LS_HYDRATE

!****************************************************
! level set initial value for gas
  subroutine INIT_LS_GAS(x,y,z,t,LS)
    IMPLICIT NONE
    real(amrex_real) x,y,z,t,LS
    
    LS = z - H_HYDRATE
  end subroutine INIT_LS_GAS

!****************************************************
! initial values for state variables for water phase
  subroutine INIT_STATE_WATER(x,y,z,t,vel,temp,dens,ccnt)
    real(amrex_real) x,y,z,t,temp,dens,ccnt
    real(amrex_real) vel(SDIM)
    integer dir

    do dir=1,SDIM
     vel(dir)=zero
    enddo 

    temp = TEMP_WATER
    dens = DENS_WATER
    ! We assume the gas concentration is at the saturation level
    ! because of the initial agitation by magnetic stir bar
    ! It should be evaluated by the gas law equation 
    ccnt =  PRESSURE_GAS * HCP_GAS(TEMP_GAS) * MOLAR_MASS_GAS ! [g.cm^-3]
  end subroutine INIT_STATE_WATER

!****************************************************
! initial values for state variables for hydrate phase
  subroutine INIT_STATE_HYDRATE(x,y,z,t,vel,temp,dens,ccnt)
    real(amrex_real) x,y,z,t,temp,dens,ccnt
    real(amrex_real) vel(SDIM)
    integer dir

    do dir=1,SDIM
     vel(dir)=zero
    enddo 
    temp = TEMP_HYDRATE
    dens = DENS_HYDRATE
    ! The concentration in hydrate layer is set to the equlibrium
    ! value in water
    ! ccnt =  PRESSURE_GAS * HCP_GAS(TEMP_GAS) * MOLAR_MASS_GAS ! [g.cm^-3]
    ccnt = DENS_GAS    

  end subroutine INIT_STATE_HYDRATE

!****************************************************
! initial values for state variables for gas phase
  subroutine INIT_STATE_GAS(x,y,z,t,vel,temp,dens,ccnt)
    real(amrex_real) x,y,z,t,temp,dens,ccnt
    real(amrex_real) vel(SDIM)
    integer dir

    do dir=1,SDIM
     vel(dir)=zero
    enddo 
    temp = TEMP_GAS
    dens = DENS_GAS
    ccnt = DENS_GAS    
  end subroutine INIT_STATE_GAS

!****************************************************
! Henry solubility defined via concentration (HCP) 
! Source: Sander2015 - Compilation of Henry’s law constants for water
! as solvent 
! Input : Temperature [K]
! Output: HCP         [mol.cm^-3.Ba^-1]
  real(amrex_real) function HCP_GAS(T) 
    real(amrex_real) T ! Temperature in K
    
! For methane gas
    real(amrex_real) :: HCP_0 = 1.4e-5   ! mol.m^-3.Pa^-1
    real(amrex_real) :: dd = 1600         ! K
    
    HCP_GAS = HCP_0 * exp(dd *(1/T - 1/298.15)) !  mol.m^-3.Pa^-1
    HCP_GAS = HCP_GAS * 1.0e-7                  !  mol.cm^-3.Ba^-1

    return
  end function HCP_GAS
!****************************************************
! Henry solubility defined via molality (HBP) 
! Source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=10#Solubility
! Input : Temperature [K]
! Output: HBP         [mol.g^-1.Ba^-1]
  real(amrex_real) function HBP_GAS(T) 
    real(amrex_real) T ! Temperature [K]
    
! For methane gas
    real(amrex_real) :: HBP_0 = 1.4e-3    ! [mol.kg^-1.bar^-1]
    real(amrex_real) :: dd = 1600         ! [K]
    
    HBP_GAS = HBP_0 * exp(dd *(1/T - 1/298.15)) !  mol.kg^-1.bar^-1
    HBP_GAS = HBP_GAS * 1e-9                    !  mol.g^-1.Ba^-1

    return
  end function HBP_GAS
!**************************************************** 
! Converting concentration to fugacity 
! Source: JammaludinETAL1991 - Hydrate plugging
! problems in undersea natural gas pipelines under shutdown
! conditions 
! Input : Concentration of gas in hydrate shell [g.cm^-3]
! Input : Initial concentration of water [g.cm^-3]
! Input : Temperature [K]
! Output: Fugacity    [Ba]
  real(amrex_real) function CCNT_TO_FUG(C_i,C_w0,T) 
    real(amrex_real) C_i       ! concentration of gas in hydrate shell [g.cm^-3]
    real(amrex_real) C_w0      ! initial concentration of water [g.cm^-3]
    real(amrex_real) T         ! Temperature [K]

    real(amrex_real) H
    
    H = 1 / (MOLAR_MASS_GAS * HBP_GAS(T)) ! [Ba]

    CCNT_TO_FUG = H * C_i / C_w0           ! [Ba]

    return
  end function CCNT_TO_FUG
!**************************************************** 
! Hydrate formation rate at the hydrate/water interface
!
! Input: Concentration of gas in hydarte shell at the hydrate/water
!        interface [g.cm^-3]
! Input: Initial concentration of water [g.cm^-3]
! Input: Temperature at the interface [K]
! Input: Pressure at the interface [Ba]
! Output: Hydrate-water front speed [cm/s]  
! K_f = 6.5e-8  Initrinsic rate constant [mol.cm^-2.Ba^-1.s^-1]
  subroutine HYDRATE_FORMATION_RATE(time,C_i,C_w0,T,p,V,HYD_SAT_TEMP, &
     K_f,verb)
   real(amrex_real) C_i,C_w0,T,p,time
   real(amrex_real) V
   real(amrex_real) HYD_SAT_TEMP

   real(amrex_real) H, p_eq, C_eq
   real(amrex_real) K_f

   integer verb

   if (K_f.le.zero) then
    print *,"K_f invalid"
    stop
   endif
   if (HYD_SAT_TEMP.lt.zero) then
    print *,"HYD_SAT_TEMP invalid"
    stop
   else if (HYD_SAT_TEMP.gt.zero) then
    p_eq = P_EQU(T)

! Checking the hydrate formation condition
    if ( ((p.ge.p_eq).and.(time.gt.zero)).or. &
     ((p.ge.1.0D+19).and.(time.eq.zero))) then
     H = 1 / (MOLAR_MASS_GAS * HBP_GAS(T)) ! [Ba]
     C_eq = p_eq * HCP_GAS(T)
!       dndt = K_f * A * H * (C_i - C_eq) / C_w0 
     V = K_f * H * (C_i - C_eq) / (C_w0 * DENS_HYDRATE)
     if(verb.ge.1) then
      print*,"H    = ", H
      print*,"C_eq = ", C_eq
      print*,"C_i  = ", C_i
      print*,"C_W0 = ", C_w0
     endif

    else
!       dndt = zero
     V = zero
    end if
   else if (HYD_SAT_TEMP.eq.zero) then
    V=zero
   else
    print *,"HYD_SAT_TEMP failure"
    stop
   endif

  end subroutine HYDRATE_FORMATION_RATE
!**************************************************** 
! assume cell_volume=1
  subroutine Methane_usage(dF,dt,D,amount_used)
    real(amrex_real) df,dt,D,amount_used

    real(amrex_real) gen_hydrate,cell_volume

    cell_volume=one

    ! amount of hydrate generated [g]
    gen_hydrate = dF * cell_volume * DENS_HYDRATE

    ! 1 mole of methane hydarte is composed of 
    ! 1 mole of methane and 5.75 mole of water
    amount_used  = gen_hydrate / MOLAR_MASS_HYDRATE * MOLAR_MASS_GAS 

    ! Mass Balance 
    ! C^old * vf^old * V - amount_used [g] = C^new * vf*new * V
    ! From GODUNOV_3D.F90
    ! C^old * vf^old  - amount_used [g.cm^-3] = C^new * vf*new 
    ! So we devide by V to make the units and mass balance consistent 
    amount_used  = amount_used / cell_volume 

  end subroutine Methane_usage
!**************************************************** 
! cell_volume not needed.
  subroutine Hydrate_energy_source_term( &
     dF,dt, &
     ksource_physical, & ! fort_heatviscconst(im_source)
     ksource_model, & ! conductstate * fort_heatflux_factor(im_source)
     energy_source,LL)
    real(amrex_real), INTENT(in) :: dF,dt
    real(amrex_real), INTENT(in) :: ksource_physical
    real(amrex_real), INTENT(in) :: ksource_model
    real(amrex_real), INTENT(in) :: LL
    real(amrex_real), INTENT(out) :: energy_source

    energy_source  = abs(LL) * dF ! [erg.g^-1]
  end subroutine Hydrate_energy_source_term

!**************************************************** 
! Temperature and pressure for Liquid_Water-Hydrate-Vapor and
! Ice-Hydrate-Vapor three-phase equilibrium.  
!
! Source: Moridis2003 - Numerical studies of gas production from
! methane hydrates. (The formula is brought in "USER’S MANUAL FOR THE
! HYDRATE v1.5 OPTION OF TOUGH+ v1.5: A CODE FOR THE SIMULATION OF
! SYSTEM BEHAVIOR IN HYDRATE-BEARING GEOLOGIC MEDIA"
!
! Input: Temperature [K]
! Output: Pressure [Ba]

  real(amrex_real) function P_EQU(T)
    real(amrex_real) T ! Temperature [K]
    real(amrex_real) y ! ln(P_e) where P_e [MPa]
    
    if(T.le.273.2) then
       y =-1.94138504464560e+5 +3.31018213397926e+3 * T  &
       -2.25540264493806e+1 * T**2 +7.67559117787059e-2 * T**3  &
       -1.30465829788791e-4 * T**4 +8.86065316687571e-8 * T**5
    else 
       y = -4.38921173434628e+1 +7.76302133739303e-1 * T  &
       -7.27291427030502e-3 * T**2 +3.85413985900724e-5 * T**3  &
       -1.03669656828834e-7 * T**4 +1.09882180475307e-10* T**5
    end if
    P_EQU = exp(Y) * 1e7 ! transform from ln(MPa) to Ba
    return
  end function P_EQU
!**************************************************** 

! Boundary condition
! Left   : symmetry
!---------------------
! Density         : Reflect even
! Temperature     : Reflect even
! Velocity        : U_n  rfelcet odd, U_t Reflect even
! Concentration   : Reflect even
! Level set       : Reflect even
! Volume fraction : Reflect even
! Pressure        : Reflect even

! Right, and bottom : slip wall + fixed temperature
!---------------------
! Density         : Extrapolation (*)
! Temperature     : Dirichlet  (*)
! Velocity        : Dirichlet U_n=0 (*), U_t Extrapolation (*)
! Concentration   : Extrapolation (*)
! Level set       : Extrapolation (*)
! Volume fraction : Extrapolation (*)
! Pressure        : Reflect even

! Top : outflow
!---------------------
! Density         : Extrapolation (*)
! Temperature     : Extrapolation (*)
! Velocity        : Extrapolation (*)
! Concentration   : Extrapolation (*)
! Level set       : Extrapolation (*)
! Volume fraction : Extrapolation (*)
! Pressure        : Dirichlet (*)

! (*) : user function needed
! Extrapolation are low order (?)

!**************************************************** 
  subroutine HYD_DENS_BC(&
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
            &(HYDRATE_REACTOR->DENS_BC)"
       stop
    endif

    if(dir.eq.1) then
       if (side.eq.1) then
          print *,"This should not be called, &
               &lo side in direction 1 is symmetry BC &
               &(HYDRATE_REACTOR->DENS_BC)"
          stop
       else if (side.eq.2) then
          q_out = q_in
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->DENS_BC)"
          stop
       endif
    else if (dir.eq.2) then
       if ((side.eq.1).or.(side.eq.2)) then
          q_out = q_in
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->DENS_BC)"
          stop
       endif
    else
       print *,"invalid direction&
            & (HYDRATE_REACTOR->DENS_BC)"
       stop
    endif
  end subroutine HYD_DENS_BC

!**************************************************** 
  subroutine HYD_TEMP_BC(&
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
            &(HYDRATE_REACTOR->TEMP_BC)"
       stop
    endif

    if(dir.eq.1) then
       if (side.eq.1) then
          print *,"This should not be called, &
               &lo side in direction 1 is symmetry BC &
               &(HYDRATE_REACTOR->TEMP_BC)"
          stop
       else if (side.eq.2) then
          q_out = TEMP_WATER
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->TEMP_BC)"
          stop
       endif
    else if (dir.eq.2) then
       if ((side.eq.1).or.(side.eq.2)) then
          q_out = TEMP_WATER
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->TEMP_BC)"
          stop
       endif
    else
       print *,"invalid direction&
            & (HYDRATE_REACTOR->TEMP_BC)"
       stop
    endif
  end subroutine HYD_TEMP_BC

!**************************************************** 
  subroutine HYD_VELO_BC(&
       time, &                  ! time
       dir, &                   ! direction
       side, &                  ! side 
       veldir, &
       q_out, &                 ! value outisde 
       q_in, &                  ! value inside
       x,y,z, &                 ! boundary point position
       dx)                      ! dx

    integer dir,side,veldir
    real(amrex_real) time,x,y,z
    real(amrex_real) q_out,q_in
    real(amrex_real) dx(SDIM)

    if (SDIM.ne.2) then
       print *,"invalid system dimension &
            &(HYDRATE_REACTOR->VELO_BC)"
       stop
    endif

    if(dir.eq.1) then
       if (side.eq.1) then
          print *,"This should not be called, &
               &lo side in direction 1 is symmetry BC &
               &(HYDRATE_REACTOR->VELO_BC)"
          stop
       else if (side.eq.2) then ! xhi
          if (veldir.eq.1) then
           q_out = zero
          else if (veldir.eq.2) then ! xhi, y velocity
           q_out = q_in  ! q_in 
          else 
           print *,"veldir invalid"
           stop
          endif
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->VELO_BC)"
          stop
       endif
    else if (dir.eq.2) then
       if (side.eq.1) then
          if (veldir.eq.1) then
           q_out = q_in ! q_inn
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
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->VELO_BC)"
          stop
       endif
    else
       print *,"invalid direction&
            & (HYDRATE_REACTOR->VELO_BC)"
       stop
    endif
  end subroutine HYD_VELO_BC
!**************************************************** 
  subroutine HYD_CCNT_BC(&
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
            &(HYDRATE_REACTOR->CCNT_BC)"
       stop
    endif

    if(dir.eq.1) then
       if (side.eq.1) then
          print *,"This should not be called, &
               &lo side in direction 1 is symmetry BC &
               &(HYDRATE_REACTOR->CCNT_BC)"
          stop
       else if (side.eq.2) then
          q_out = q_in
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->CCNT_BC)"
          stop
       endif
    else if (dir.eq.2) then
       if (side.eq.1) then
          q_out = q_in
       else if (side.eq.2) then
          q_out = DENS_GAS
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->CCNT_BC)"
          stop
       endif
    else
       print *,"invalid direction&
            & (HYDRATE_REACTOR->CCNT_BC)"
       stop
    endif
  end subroutine HYD_CCNT_BC

!**************************************************** 
  subroutine HYD_LVLS_BC(&
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
            &(HYDRATE_REACTOR->LVLS_BC)"
       stop
    endif

    if(dir.eq.1) then
       if (side.eq.1) then
          print *,"This should not be called, &
               &lo side in direction 1 is symmetry BC &
               &(HYDRATE_REACTOR->LVLS_BC)"
          stop
       else if (side.eq.2) then
          q_out = q_in
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->LVLS_BC)"
          stop
       endif
    else if (dir.eq.2) then
       if ((side.eq.1).or.(side.eq.2)) then
          q_out = q_in
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->LVLS_BC)"
          stop
       endif
    else
       print *,"invalid direction&
            & (HYDRATE_REACTOR->LVLS_BC)"
       stop
    endif
  end subroutine HYD_LVLS_BC


!**************************************************** 
  subroutine HYD_VOLF_BC(&
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
            &(HYDRATE_REACTOR->VOLF_BC)"
       stop
    endif

    if(dir.eq.1) then
       if (side.eq.1) then
          print *,"This should not be called, &
               &lo side in direction 1 is symmetry BC &
               &(HYDRATE_REACTOR->VOLF_BC)"
          stop
       else if (side.eq.2) then
          q_out = q_in
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->VOLF_BC)"
          stop
       endif
    else if (dir.eq.2) then
       if ((side.eq.1).or.(side.eq.2)) then
          q_out = q_in
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->VOLF_BC)"
          stop
       endif
    else
       print *,"invalid direction&
            & (HYDRATE_REACTOR->VOLF_BC)"
       stop
    endif
  end subroutine HYD_VOLF_BC

!**************************************************** 
  subroutine HYD_MOMF_BC(&
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
    real(amrex_real) q_out(SDIM), q_in(SDIM)
    real(amrex_real) dx(SDIM)

    if (SDIM.ne.2) then
       print *,"invalid system dimension &
            &(HYDRATE_REACTOR->MOMF_BC)"
       stop
    endif

    if(dir.eq.1) then
       if (side.eq.1) then
          print *,"This should not be called, &
               &lo side in direction 1 is symmetry BC &
               &(HYDRATE_REACTOR->MOMF_BC)"
          stop
       else if (side.eq.2) then
          q_out(1)=-q_in(1)
          q_out(2)= q_in(2)
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->MOMF_BC)"
          stop
       endif
    else if (dir.eq.2) then
       if ((side.eq.1).or.(side.eq.2)) then
          q_out(1)= q_in(1)
          q_out(2)=-q_in(2)
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->MOMF_BC)"
          stop
       endif
    else
       print *,"invalid direction&
            & (HYDRATE_REACTOR->MOMF_BC)"
       stop
    endif
  end subroutine HYD_MOMF_BC

!**************************************************** 
  subroutine HYD_PRES_BC(&
       time, &                  ! time
       dir, &                   ! direction
       side, &                  ! side 
       q_out, &                 ! value outisde 
       x_wall, &                ! wall position
       q_in, &                  ! value inside
       x,y,z, &                 ! boundary point position
       dx)                      ! dx

    integer dir,side
    real(amrex_real) time,x_wall,x,y,z
    real(amrex_real) q_out, q_in
    real(amrex_real) dx(SDIM)

    if (SDIM.ne.2) then
       print *,"invalid system dimension &
            &(HYDRATE_REACTOR->PRES_BC)"
       stop
    endif

    if(dir.eq.1) then
       if (side.eq.1) then
          print *,"This should not be called, &
               &lo side in direction 1 is symmetry BC &
               &(HYDRATE_REACTOR->PRES_BC)"
          stop
       else if (side.eq.2) then
          q_out=q_in
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->PRES_BC)"
          stop
       endif
    else if (dir.eq.2) then
       if (side.eq.1) then
          q_out= q_in
       else if (side.eq.2) then
          q_out=PRESSURE_GAS
       else 
          print *,"invalid side in x dirction &
               &(HYDRATE_REACTOR->PRES_BC)"
          stop
       endif
    else
       print *,"invalid direction&
            & (HYDRATE_REACTOR->PRES_BC)"
       stop
    endif
  end subroutine HYD_PRES_BC

!**************************************************** 
end module hydrateReactor_module
