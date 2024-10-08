      program main
      IMPLICIT NONE

      real*8 mypi,dx,sigma,rho,mu_cutoff

      mypi=4.0d0*atan(1.0d0)
      dx=0.2d0/48.0
      sigma=76.0D-5
      rho=1.0d0
c dt_tension=dx^(3/2) sqrt(rho_bar)/sqrt(sigma pi)
c dt_viscosity=dx^2 rho/(6 mu)
      mu_cutoff=sqrt(dx)/(6.0d0*sqrt(1.0d0/(mypi*sigma*rho)))

      print *,"dx=",dx
      print *,"sigma=",sigma
      print *,"rho=",rho
      print *,"mu_cutoff=",mu_cutoff
  
      end

