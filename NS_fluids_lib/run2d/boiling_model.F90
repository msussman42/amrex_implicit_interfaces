      program main
      IMPLICIT NONE

!      real*8,parameter      :: macrolayer = 2.8D-4
!      real*8,parameter      :: microlayer = 1.0D-9
      real*8,parameter      :: macrolayer = 3.91D-5  ! 0.732
!      real*8,parameter      :: macrolayer = 7.82D-5  ! 0.390
!      real*8,parameter      :: macrolayer = 1.95D-5  ! 1.3
!      real*8,parameter      :: microlayer = 1.0D-9 ! 0.732
      real*8,parameter      :: microlayer = 1.0D-8 ! 0.573
      real*8,parameter      :: TSAT = 373.0D0
      real*8,parameter      :: LL = 2.257D+6
!      real*8,parameter      :: T_HOT = 381.5D0  !0.732
       ! QUESTION FOR MITSUHIRO: is there heat resistance between
       ! the heater and the liquid? in otherwords,
       ! does the water touch directly the heater?
       ! if V=velocity_vapor + mdot/density_vapor, then
       ! better agreement between analytical results and simulation.
!      real*8,parameter      :: T_HOT = 383.0D0  ! 0.815
      real*8,parameter      :: T_HOT = 381.5D0  !0.732
!      real*8,parameter      :: T_HOT = 376.5D0  ! 0.285
      real*8,parameter      :: kwater = 0.68D0
      real*8,parameter      :: micro_angle=0.663 ! 0.732 velsrc
!      real*8,parameter      :: micro_angle=0.563  ! 0.721 velsrc
!      real*8,parameter      :: micro_angle=0.01 ! 0.692 velsrc

      real*8 factor,micro_slope,psi_upper,psi_lower,velsrc

      micro_slope=tan(0.5d0*micro_angle)
      psi_upper=0.5d0*macrolayer/micro_slope
      psi_lower=0.5d0*microlayer/micro_slope
      velsrc=( &
        0.5d0*(TSAT+T_HOT)-Tsat)* &
              log(psi_upper/psi_lower)* &
              sqrt(1.0d0+micro_slope**2)/ &
              (micro_slope*(psi_upper-psi_lower))
      velsrc=abs(kwater*velsrc/LL)

      print *,"angle=",micro_angle
      print *,"macro=",macrolayer
      print *,"micro=",microlayer
      print *,"TSAT=",TSAT
      print *,"T_HOT=",T_HOT
      print *,"kwater=",kwater
      print *,"LL=",LL
      print *,"velsrc=",velsrc 
      return
      end

