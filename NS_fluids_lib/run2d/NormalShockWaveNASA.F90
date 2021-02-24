      program main
      IMPLICIT NONE
       ! shock wave is stationary
      real*8 :: gamma
      real*8 :: rho0,u0,c0,p0,T0,M0  ! upstream supersonic
      real*8 :: rho1,u1,c1,p1,T1,M1  ! downstream subsonic
      real*8 :: cv
      real*8 :: jump
      integer :: MKS_or_CGS

      gamma=1.4d0
      M0=2.0d0  ! Mach 2 upstream
      T0=273.0d0  ! Kelvin

      MKS_or_CGS=1

      if (MKS_or_CGS.eq.0) then
       rho0=0.001d0
       cv=7.18d+6 ! erg/g
      else if (MKS_or_CGS.eq.1) then
       rho0=1.0d0
       cv=7.18d+2 ! J/kg
      else
       print *,"MKS_or_CGS invalid"
       stop
      endif

      p0=(gamma-1.0d0)*rho0*cv*T0
      c0=sqrt(gamma*p0/rho0)
      u0=M0*c0

      p1=p0*(2.0d0*gamma*M0*M0-(gamma-1.0d0))/(gamma+1.0d0)
      rho1=rho0*(gamma+1.0d0)*M0*M0/((gamma-1.0d0)*M0*M0+2.0d0)
      T1=T0*(2.0d0*gamma*M0*M0-(gamma-1.0d0))* &
            ((gamma-1.0d0)*M0*M0+2.0d0)/((M0*(gamma+1.0d0))**2)
      M1=sqrt(((gamma-1.0d0)*M0*M0+2.0d0)/ &
              (2.0d0*gamma*M0*M0-(gamma-1.0d0)))
      c1=sqrt(gamma*p1/rho1)
      u1=M1*c1

      jump=rho0*u0-rho1*u1
      print *,"rho0 u0 - rho1 u1 = ",jump
      jump=rho0*u0*u0+p0-(rho1*u1*u1+p1)
      print *,"rho0 u0^2+p0 - (rho1 u1^2+p1) = ",jump
      jump=(2.0d0/(gamma-1.0d0))*c0*c0+u0*u0- &
           ( (2.0d0/(gamma-1.0d0))*c1*c1+u1*u1 )
      print *,"[2 c^2 /(gamma-1)+u^2] = ",jump

      print *,"UPSTREAM rho,u,p,c,T,M ",rho0,u0,p0,c0,T0,M0
      print *,"DOWNSTREAM rho,u,p,c,M ",rho1,u1,p1,c1,T1,M1

      return
      end

