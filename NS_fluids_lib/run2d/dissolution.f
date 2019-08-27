     
      subroutine tridiag(a,b,c,r,u,n)
      IMPLICIT NONE

      integer n
      integer nmax
 
      real*8 a(n)
      real*8 b(n)
      real*8 c(n)
      real*8 r(n)
      real*8 u(n)

      parameter (nmax = 8000)

      integer j
      real*8 bet
      real*8 gam(nmax)

      if (n .gt. nmax ) then
       print *,"tridiag: size exceeded"
       stop
      endif
      if (b(1) .eq. 0) then
       print *,"tridiag: CANT HAVE B(1) = ZERO"
       stop
      endif

      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) then
         print *,"tridiag: TRIDIAG FAILED"
         stop
        endif
        u(j) = (r(j)-a(j)*u(j-1))/bet
      end do

      do j = n-1,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      end do

      return
      end
 
      program main
      IMPLICIT NONE
      real*8 dt,k,dx,L,tstop,layer,t
      real*8 ulo,uhi,umiddle,alpha
      integer i,j,ncell,ntri
      real*8 u(0:8000),x(0:8000),unew(0:8000)
      real*8 a(8000),b(8000),c(8000),r(8000),utri(8000)
      real*8 rstep,vel,lower,upper

      print *,"this program solves u_t+vel u_x=ku_xx   0<=x<=L"
      
c dissolution coefficient
      k=4.0
      vel=0.0
      L=512.0
      layer=0.0
      ulo=2.0
      uhi=1.0
      umiddle=1.0
      ncell=4000
      dx=L/ncell
      tstop=52.0
      dt=tstop/ncell
     
      print *,"ulo,uhi: ",ulo,uhi
      print *,"ncell,dx,dt ",ncell,dx,dt
      print *,"k,L,tstop,layer ",k,L,tstop,layer 

      do i=0,ncell
       x(i)=i*dx
       u(i)=umiddle
       if (x(i).le.layer) then
        u(i)=ulo
       endif
      enddo
      u(0)=ulo
      u(ncell)=uhi 

      t=0.0
      rstep=0.0
10    continue 
      do i=0,ncell
        unew(i)=u(i)
      enddo
      alpha=k*dt/(dx*dx)
      do i=1,ncell
       lower=-alpha
       upper=-alpha
       if (vel.gt.0.0) then
        lower=lower-dt*vel/dx
       else 
        upper=upper+dt*vel/dx
       endif
       a(i)=lower
       c(i)=upper
       b(i)=1.0-lower-upper
c right hand side
       r(i)=u(i)
       if (i.eq.1) then
        r(1)=u(1)-lower*u(0)
       endif
       if (i.eq.ncell-1) then
        r(ncell-1)=u(ncell-1)-upper*u(ncell)
       endif
      enddo
      ntri=ncell-1
      call tridiag(a,b,c,r,utri,ntri)
      do i=1,ncell-1
       unew(i)=utri(i)
      enddo
      do i=0,ncell
       u(i)=unew(i)
      enddo
      t=t+dt
      rstep=rstep+1.0
      if (t.lt.tstop) then
       goto 10
      endif

      print *,"rstep,t ",rstep,t
      do i=0,ncell
       print *,x(i),u(i)
      enddo

      end

