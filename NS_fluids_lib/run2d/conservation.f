      program main
      IMPLICIT NONE

      integer select
      real*8 A,T,gamma
      real*8 mass_flux,mom_flux,energy_flux
      real*8 rhoL,uL,pL,eL,energyL
      real*8 rhoR,uR,pR,eR,energyR
      real*8 init_mass,init_mom,init_energy

      select=1
      if (select.eq.1) then
       A=0.4d0
       T=1.8d0
       gamma=1.4d0
       rhoL=3.857148d0
       uL=2.629369d0
       pL=10.3333333d0
       rhoR=1.0d0+0.2d0*sin(45.0d0)
       uR=0.0d0
       pR=1.0d0
       eL=pL/(rhoL*(gamma-1.0d0))
       eR=pR/(rhoR*(gamma-1.0d0))
       energyL=0.5d0*uL*uL+eL
       energyR=0.5d0*uR*uR+eR
       mass_flux=A*T*(rhoL*uL-rhoR*uR)
       mom_flux=A*T*(rhoL*uL*uL+pL-rhoR*uR*uR-pR)
       energy_flux=A*T*(rhoL*uL*energyL+uL*pL-
     &    (rhoR*uR*energyR+uR*pR))
       init_mass=5.1504589949959d0
       init_mom=4.056746151845d0
       init_energy=24.66667428626d0
      endif
      print *,"mass_flux ",mass_flux
      print *,"mom_flux ",mom_flux
      print *,"energy_flux ",energy_flux
      print *,"mass ",mass_flux+init_mass
      print *,"mom ",mom_flux+init_mom
      print *,"energy ",energy_flux+init_energy
      end
