      program main
      IMPLICIT NONE

      real*8 rho,cv,k,thick
      real*8 q,TX,THEAT,alpha

      rho=4.0E+3
      cv=7.5E+2
      k=3.5E+1
      thick=0.25E-3
      q=30.0e+3
 
      rho=4.0
      cv=7.5E+6
      k=3.5E+6
      thick=0.25E-1
      q=30.0e+6

      rho=1.0
      cv=4.22E+7
      k=6.79E+4
      thick=1.0
      q=30.0e+6

      alpha=k/(rho*cv)
      TX=q/k
      THEAT=TX*thick

      print *,"THEAT=",THEAT
      end

