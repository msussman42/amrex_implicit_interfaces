      program main
      IMPLICIT NONE

      real*8 dxmin,tension,den1,den2,omega,k,tension_target
      real*8 dt_target
      real*8 dt
      real*8 wavespeed
      real*8 mypi

      mypi=4.0d0*atan(1.0d0)
      dxmin=1.0d0/128.0d0
      den1=1.0d+6
      den2=1.0d+6
      k=2.0d0*mypi/dxmin
      tension=831.0d0
      omega=(k**1.5d0)*sqrt(tension/(den1+den2))
      wavespeed=omega/k
      dt=dxmin/wavespeed

c     dt=dxmin/(omega/k)=k dxmin/omega=k dxmin/(omega(1)*sqrt(tension)) 
c     sqrt(tension)=k dxmin/(omega(1)*dt)

      dt_target=1.0D+9
      tension_target=k*dxmin/(omega*dt_target)
      tension_target=tension_target*tension_target
      print *,"dx=",dxmin
      print *,"wavespeed=",wavespeed
      print *,"tension_target=",tension_target
  
      end

