      subroutine ptillotson(p,rho)
      IMPLICIT NONE

      real*8 rho,p
      real*8 A_TILLOTSON,B_TILLOTSON,C_TILLOTSON
      real*8 RHOBAR_TILLOTSON,OMEGA_TILLOTSON,P0_TILLOTSON,mu,rho0
      real*8 internal_energy,e0
      real*8 cv

      cv=4.1855D+7
      e0=293.0d0*cv
      internal_energy=293.0d0*cv

      P0_TILLOTSON=1.0D+6
      A_TILLOTSON=2.2D+10
      B_TILLOTSON=9.94D+10
      C_TILLOTSON=1.457D+11
      RHOBAR_TILLOTSON=1.0D0
      OMEGA_TILLOTSON=0.28D0

      rho0=1.0d0
      mu=rho/rho0-1.0d0
      print *,"mu (ptillotson): ",mu
      print *,"e0 (ptillotson): ",e0
      p=P0_tillotson+omega_tillotson*rho* &
              (internal_energy-e0)+ &
              A_tillotson*mu+ &
              B_tillotson*(mu**2)+ &
              C_tillotson*(mu**3)

      return
      end

      subroutine soundtillotson(p,rho,sound)
      IMPLICIT NONE

      real*8 rho,p
      real*8 A_TILLOTSON,B_TILLOTSON,C_TILLOTSON
      real*8 RHOBAR_TILLOTSON,OMEGA_TILLOTSON,P0_TILLOTSON,mu,rho0
      real*8 internal_energy,e0,sound,dmu,dpdrho,dpde
      real*8 cv

      cv=4.1855D+7
      e0=293.0d0*cv
      internal_energy=293.0d0*cv
      P0_TILLOTSON=1.0D+6
      A_TILLOTSON=2.2D+10
      B_TILLOTSON=9.94D+10
      C_TILLOTSON=1.457D+11
      RHOBAR_TILLOTSON=1.0D0
      OMEGA_TILLOTSON=0.28D0

      rho0=1.0d0
      mu=rho/rho0-1.0d0
      print *,"mu (soundtillotson): ",mu
      print *,"e0 (soundtillotson): ",e0

      dmu=1.0d0/rho0
      dpdrho=omega_tillotson*(internal_energy-e0)+ &
              A_TILLOTSON*dmu+ &
              2.0d0*B_TILLOTSON*mu*dmu+ &
              3.0d0*C_TILLOTSON*mu*mu*dmu
      dpde=omega_tillotson*rho
      print *,"dpdrho=",dpdrho
      print *,"dpde=",dpde
      call ptillotson(p,rho)
      sound=dpdrho+p*dpde/(rho**2)
      sound=sqrt(sound) 

      return
      end
 
      program main
      IMPLICIT NONE
      real*8 minimum_pressure,minimum_density
      real*8 minimum_sound

      minimum_pressure=5.0d+4
      minimum_density=0.995d0
      minimum_density=1.0d0
      minimum_density=0.99995775d0
      minimum_density=0.9999569d0
      minimum_density=0.99995685d0
      minimum_density=0.99995681d0
      call ptillotson(minimum_pressure,minimum_density)
      print *,"minimum_pressure= ",minimum_pressure
      print *,"minimum_density= ",minimum_density
 
      call soundtillotson(minimum_pressure,minimum_density,minimum_sound)
      print *,"minimum_sound= ",minimum_sound

      return
      end

