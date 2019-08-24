      program main
      IMPLICIT NONE

      real*8 Z,L,theta,psi,X,R,PI,term,costerm,R0,V0

      Z=1.02
      L=2.0
      X=(L*L-Z*Z)/(2.0*Z)
      R=X+Z
      theta=asin(X/R) 
      PI=4.0*atan(1.0)
      psi=0.5*PI-theta
      costerm=cos(psi)
      term=(2.0+costerm**3)/3.0-costerm
      term=(2.0/(3.0*term))
      term=term**(1.0/3.0)
      R0=R/term 
      V0=(2.0/3.0)*PI*(R0**3)
      print *,"R= ",R
      print *,"R0= ",R0
      print *,"angle (radians) =",psi
      print *,"angle (degrees) =",180.0*psi/PI
      print *,"V0= ",V0
  
      end

