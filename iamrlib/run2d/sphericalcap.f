      program main
      IMPLICIT NONE

      real*8 mypi
      real*8 r,h,vol,effective_r

      mypi=4.d0*atan(1.d0)
      r=2.0E-3
      h=(2.5E-4)-(1.536E-3)+r
      vol=mypi*h*h*(3.0*r-h)/3.0
      vol=4.0*mypi*(r**3.0)/3.0 - vol
       ! V=(4/3) pi r^3
       ! 3V/(4 pi) = r^3
      effective_r=(3.0*vol/(4.0*mypi))**0.33333

      print *,"r= ",r
      print *,"h=",h
      print *,"vol= ",vol
      print *,"effective_r= ",effective_r
  
      end

