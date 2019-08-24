      subroutine getweight(a,b,mu,x,wt)
      IMPLICIT NONE

      real*8 a,b,mu,x,wt

      wt=(abs(x-a)**(1.0-mu)-abs(x-b)**(1.0-mu))/(1.0-mu)
      wt=abs(wt)

      return
      end


      subroutine fractalsmooth(ux,ncell,dx,mu,ulo,uhi,x)
      IMPLICIT NONE
 
      integer ncell,i,j
      real*8 ulo,uhi,dx,mu,a,b,wt,totalweight,xface
      real*8 ux(-1:8000),x(-1:8000),uxsmooth(-1:8000)

      do i=0,ncell
       if (mu.ge.1.0) then
        uxsmooth(i)=ux(i)
       else
        xface=i*dx
        totalweight=0.0
        uxsmooth(i)=0.0
        do j=0,ncell
         if (j.gt.0) then
          a=(j-0.5)*dx
          b=j*dx
          call getweight(a,b,mu,xface,wt) 
          uxsmooth(i)=uxsmooth(i)+wt*ux(j)
          totalweight=totalweight+wt
         endif
         if (j.lt.ncell) then
          a=(j+0.5)*dx
          b=j*dx
          call getweight(a,b,mu,xface,wt) 
          uxsmooth(i)=uxsmooth(i)+wt*ux(j)
          totalweight=totalweight+wt
         endif
        enddo
        uxsmooth(i)=uxsmooth(i)/totalweight
       endif
      enddo

      do i=0,ncell
       ux(i)=uxsmooth(i)
      enddo

      return
      end

      program main
      IMPLICIT NONE
      real*8 dt,k,dx,L,tstop,layer,t
      real*8 ulo,uhi,umiddle,alpha
      integer i,j,ncell,ntri
      real*8 u(-1:8000),x(-1:8000),unew(-1:8000)
      real*8 ux(-1:8000)
      real*8 rstep,vel,lower,upper,mu
      integer step

      print *,"this program solves u_t=D^{1+mu}u   0<=x<=L"
     
      L=1.0
      ncell=100
      dx=L/ncell
      tstop=0.01
      dt=0.25*dx*dx
      ulo=0.0
      uhi=1.0 
      mu=0.1 
     
      print *,"ncell,dx,dt ",ncell,dx,dt
      print *,"L,tstop ",L,tstop
      print *,"ulo,uhi ",ulo,uhi

      do i=-1,ncell
       x(i)=(i+0.5)*dx
       if (x(i).le.0.5) then
        u(i)=ulo
       else
        u(i)=uhi
       endif
      enddo

      do i=-1,ncell
        unew(i)=u(i)
      enddo

      t=0.0 
      step=0
      do while (t.lt.tstop)
       do i=0,ncell
        ux(i)=(u(i)-u(i-1))/dx 
       enddo
       call fractalsmooth(ux,ncell,dx,mu,ulo,uhi,x)
       do i=0,ncell-1
        unew(i)=u(i)+dt*(ux(i+1)-ux(i))/dx
       enddo
       unew(-1)=ulo
       unew(ncell)=uhi
       do i=-1,ncell
        u(i)=unew(i)
       enddo
       t=t+dt
       step=step+1
       print *,"step,t,tstop ",step,t,tstop
      enddo

      do i=0,ncell-1
        print *,x(i),u(i)
      enddo

      end

