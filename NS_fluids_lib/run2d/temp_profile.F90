      program main
      IMPLICIT NONE

      real*8,parameter      :: p1 = -0.94381
      real*8,parameter      :: p2 = 10.921
      real*8,parameter      :: p3 = -43.593
      real*8,parameter      :: p4 = 51.946
      real*8,parameter      :: p5 = 73.773
      real*8,parameter      :: p6 = -216.18
      real*8,parameter      :: p7 = 465.15
      integer N,i
      real*8 h,y,T

      N=100
      h=3.0/N
      do i=0,N
       y=i*h
       T=p1*(y**6.0d0) + p2*(y**5.0d0) + p3*(y**4.0d0) &
        +p4*(y**3.0d0) + p5*(y**2.0d0) + p6*y + p7
       print *,"y,T ",y,T
      enddo
 
      return
      end

