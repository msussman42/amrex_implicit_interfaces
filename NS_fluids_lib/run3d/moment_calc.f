      program main
      IMPLICIT NONE

      integer ncell
      integer i,j,k
      real*8 xlo,xhi,h,radblob
      real*8 x,y,z,vol,r2_mom

      xlo=-1.0d0
      xhi=1.0d0
      ncell=128
      h=(xhi-xlo)/ncell
      radblob=1.0d0

      vol=0.0d0
      r2_mom=0.0d0

      do k=0,ncell-1
      do j=0,ncell-1
      do i=0,ncell-1
       x=xlo+(i+0.5d0)*h
       y=xlo+(j+0.5d0)*h
       z=xlo+(k+0.5d0)*h
       if (x**2+y**2+z**2.le.radblob**2) then
        vol=vol+h**3
        r2_mom=r2_mom+(x**2+y**2)*(h**3)
       endif
      enddo
      enddo
      enddo
      print *,"vol=",vol
      print *,"r2_mom=",r2_mom
  
      end

