      subroutine get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      IMPLICIT NONE

      real(8), intent(out) :: A,B,GAMMA,R1,R2,RHOI

      if (1.eq.0) then  ! spherical explosion
       A=5.484D+12
       B=0.09375D+12
       R1=4.94D0
       R2=1.21D0
       GAMMA=1.28D0
       RHOI=1.63D0
      else if (1.eq.1) then !hydrobulge
!      A=6.17D+12 !cgs
       A=6.1327D+12 !cgs
!      B=1.69D+11 !cgs 
       B=1.5069D+11 !cgs 
       R1=4.4D0
       R2=1.2D0
       GAMMA=1.25D0
       RHOI=1.765D0 !cgs
      else if (1.eq.0) then  ! bubble jetting
       A=3.712D+12
       B=0.03231D+12
       R1=4.15D0
       R2=0.95D0
       GAMMA=1.3D0
       RHOI=1.63D0
      else if (1.eq.0) then ! cavitation
       A=3.712D+12
       B=0.03231D+12
       R1=4.15D0
       R2=0.95D0
       GAMMA=1.3D0
       RHOI=1.63D0
      else
       print *,"get_jwl_constants failure"
       stop
      endif

      return
      end subroutine get_jwl_constants


      subroutine INTERNAL_jwl(rho,temperature,internal_energy)
      IMPLICIT NONE

      real(8), intent(in) :: rho,temperature
      real(8), intent(out) :: internal_energy
      real(8) :: cv

      if (rho.gt.0.0d0) then
       !do nothing
      else
       print *,"density negative INTERNAL_jwl:",rho
       stop
      endif
      if (temperature.gt.0.0d0) then
       !do nothing
      else
       print *,"temperature <=0 INTERNAL_jwl: ",temperature
       stop
      endif

      cv=4.1855D+7
      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_jwl

      subroutine TEMPERATURE_jwl(rho,temperature,internal_energy)
      IMPLICIT NONE

      real(8), intent(in) :: rho,internal_energy
      real(8), intent(out) :: temperature
      real(8) :: cv

      if (rho.gt.0.0d0) then
       !do nothing
      else
       print *,"density negative TEMPERATURE_jwl:",rho
       stop
      endif
      if (internal_energy.gt.0.0d0) then
       !do nothing
      else
       print *,"internal energy <=0 TEMPERATURE_jwl:",internal_energy
       stop
      endif

      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_jwl

! e=(E/rho) - (1/2) (u^2 + v^2)
      subroutine EOS_NAjwl(rho,internal_energy,pressure)
      IMPLICIT NONE

      real(8), intent(in) :: rho,internal_energy
      real(8), intent(out) :: pressure
      real(8) A,B,R1,R2,GAMMA,RHOI,OMEGA

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      OMEGA=GAMMA-1.0d0

      if (rho.gt.0.0d0) then
       !do nothing
      else
       print *,"rho invalid EOS_NAjwl: ",rho
       stop
      endif
      if (internal_energy.gt.0.0d0) then
       !do nothing
      else
       print *,"e invalid EOS_NAjwl: ",internal_energy
       stop
      endif
      pressure= &
        A*(1.0d0-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(1.0d0-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)+ &
        OMEGA*rho*internal_energy

      if (pressure.gt.0.0d0) then
       !do nothing
      else
       print *,"vacuum error in NA JWL: ",pressure
       stop
      endif

      return
      end subroutine EOS_NAjwl

! initial sound speed is:
! C=7.8039D+10-5.484D+12 e^(-4.94)-0.09375D+12 e^(-1.21)=
      subroutine SOUNDSQR_NAjwl(rho,internal_energy,soundsqr)
      IMPLICIT NONE

      real(8), intent(in) :: rho,internal_energy
      real(8), intent(out) :: soundsqr
      real(8) A,B,R1,R2,GAMMA,RHOI,OMEGA
      real(8) pressure,dp_de,dp_drho

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      OMEGA=GAMMA-1.0d0

      if (rho.gt.0.0d0) then
       !do nothing
      else
       print *,"rho invalid SOUNDSQR_NAjwl: ",rho
       stop
      endif
      if (internal_energy.gt.0.0d0) then
       !do nothing
      else
       print *,"e invalid SOUNDSQR_NAjwl: ",internal_energy
       stop
      endif

      call EOS_NAjwl(rho,internal_energy,pressure)
      dp_de=OMEGA*rho
      dp_drho= &
        A*(1.0d0-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)* &
        R1*RHOI/(rho**2)- &
        (A*OMEGA/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(1.0d0-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)* &
        R2*RHOI/(rho**2)- &
        B*(OMEGA/(R2*RHOI))*exp(-R2*RHOI/rho)+ &
        OMEGA*internal_energy
    
      soundsqr=(pressure*dp_de)/(rho**2)+dp_drho
 
      if (soundsqr.gt.0.0d0) then
       !do nothing
      else
       print *,"soundsqr invalid in SOUNDSQR_NAjwl:",soundsqr
       stop
      endif

      return
      end subroutine SOUNDSQR_NAjwl

      program main
      IMPLICIT NONE
      real*8 E0,E0_per_mass,RHOI,T0,P0,C2,C

      RHOI=1.765d0 ! cgs
      E0=10.1D+10 ! cgs ergs/cm^3=g cm^2/s^2/cm^3=g/(s^2 cm)
      E0_per_mass=E0/RHOI
      call TEMPERATURE_jwl(RHOI,T0,E0_per_mass)
      call EOS_NAjwl(RHOI,E0_per_mass,P0)
      call SOUNDSQR_NAjwl(RHOI,E0_per_mass,C2)
      C=sqrt(C2)

      print *,"E0=",E0
      print *,"E0_per_mass=",E0_per_mass
      print *,"T0=",T0
      print *,"RHOI=",RHOI
      print *,"P0=",P0
      print *,"C=",C

      return
      end

