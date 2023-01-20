      program main
      IMPLICIT NONE

      real*8 dxmin,k,TDIFF,LL,den
      real*8 wavespeed
      real*8 mypi

      mypi=4.0d0*atan(1.0d0)
      dxmin=1.0d0/128.0d0
c erg/(cm s K)=(g cm^2/s^2)/(cm s K)=
c g cm/(s^3 K)

      k=64.0d+5 
c erg/g 
      LL=18.05d+9
c degrees Kelvin
      TDIFF=20.0d0
      den=2.3d0
c erg/(cm s K) * (K/cm) /(erg/g) * (cm^3/g) =cm/s
      wavespeed=k*(TDIFF/dxmin)/(LL*den)

      print *,"dx=",dxmin
      print *,"wavespeed=",wavespeed
  
      end

