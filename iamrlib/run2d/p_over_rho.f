      program main
      IMPLICIT NONE

      real*8 K0,n,rho,rho0,P0,rhocav
      integer i
      real*8 P,Pbyrho

      K0=2.15E+9
      n=7.15
      rho0=1000.0
      P0=1.01325E+5

      rhocav=999.9577

      do i=0,1000
       rho=rhocav+(rho0-rhocav)*i/1000.0
       P=((rho/rho0)**n-1.0)*K0/n+P0
       Pbyrho=P/rho
       print *,"rho,Pbyrho ",rho,Pbyrho
      enddo
 
      end

