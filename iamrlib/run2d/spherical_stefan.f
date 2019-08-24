c  solves Ax=b, status=1 if ok, A is an ndim x ndim matrix 
c  in fortran, array dimensions go from 1..n; thus
c  valid regions of AA are 1..ndim x 1..ndim
c  NOTE: AA and bb are OVERWRITTEN after call to matrix_solve.

        subroutine matrix_solve(AA,xx,bb,ndim,status)
        integer ndim
        real*8 AA(ndim,ndim),xx(ndim),bb(ndim)
        real*8 alpha,holdvalue
        integer i,j,k,holdj,status

        status=1
        do i=1,ndim-1
         holdj=i
         holdvalue=AA(i,i)
         do j=i+1,ndim 
          if (abs(AA(j,i)).gt.abs(holdvalue)) then
           holdj=j
           holdvalue=abs(AA(j,i))
          endif
         enddo
         if (holdj.ne.i) then
          do j=i,ndim
           holdvalue=AA(i,j)
           AA(i,j)=AA(holdj,j)
           AA(holdj,j)=holdvalue
          enddo
         endif
         holdvalue=bb(i)
         bb(i)=bb(holdj)
         bb(holdj)=holdvalue
         if (abs(AA(i,i)).lt.1.0E-13) then
          status=0
         else
          do j=i+1,ndim
           alpha=AA(j,i)/AA(i,i)
           do k=i,ndim
            AA(j,k)=AA(j,k)-alpha*AA(i,k)
           enddo
           bb(j)=bb(j)-alpha*bb(i)
          enddo
         endif
        enddo

        do i=ndim,1,-1
         if (status.ne.0) then
          holdvalue=bb(i)
          do j=i+1,ndim
           holdvalue=holdvalue-AA(i,j)*xx(j)
          enddo
          if (abs(AA(i,i)).lt.1.0E-13) then
           status=0
          else
           xx(i)=holdvalue/AA(i,i)
          endif
         endif
        enddo

        return
        end

      subroutine tridiag_solve(l,u,d,n,f,soln)
      IMPLICIT NONE

      integer n,i
      real*8 l(n),u(n),d(n),f(n),soln(n)
      real*8 ll(n),uu(n),dd(n),z(n)

      dd(1)=d(1)
      uu(1)=u(1)/dd(1)
      z(1)=f(1)/dd(1)
      do i=2,n-1
       ll(i)=l(i)
       dd(i)=d(i)-ll(i)*uu(i-1)
       uu(i)=u(i)/dd(i)
       z(i)=(f(i)-ll(i)*z(i-1))/dd(i)
      enddo
      ll(n)=l(n)
      dd(n)=d(n)-ll(n)*uu(n-1)
      z(n)=(f(n)-ll(n)*z(n-1))/dd(n)
      soln(n)=z(n)
      do i=n-1,1,-1
       soln(i)=z(i)-uu(i)*soln(i+1)
      enddo

      return
      end

      real*8 function erfmine(x)
      IMPLICIT NONE

      real*8 x,pi,dx,t,xhead,xtail
      real*8 sumhead,sumtail
      integer i,nhead,ntail

      pi=4.d0*atan(1.d0)
      if (x.le.3.0) then
       xhead=x
       nhead=1200
       xtail=x
       ntail=0
      else 
       xhead=3.0
       nhead=1000
       xtail=x
       ntail=200
      endif

      dx=xhead/nhead
      sumhead=0.0
      do i=0,nhead-1
       t=(i+0.5)*dx
       sumhead=sumhead+exp(-t*t)*dx
      enddo
      sumtail=0.0
      if ((ntail.gt.0).and.(xtail-xhead.gt.0.0)) then
       dx=(xtail-xhead)/ntail
       do i=0,ntail-1
        t=xhead+(i+0.5)*dx
        sumtail=sumtail+exp(-t*t)*dx
       enddo
      endif
      erfmine=(sumtail+sumhead)*2.0/sqrt(pi)

      return
      end
       
      real*8 function transf(x,cp)
      IMPLICIT NONE

      real*8 erfmine
      real*8 x,pi,cp

      pi=4.d0*atan(1.d0)
      transf=x*exp(x*x)*erfmine(x)-cp/sqrt(pi)

      return
      end

      subroutine integrand_function(x,eps,beta,ff)
      IMPLICIT NONE

      real*8 x,eps,beta,ff

      ff=2.0*(beta**3)*exp((beta-x)*(beta+x-2.0*eps*beta*beta/x))/(x*x)

      return
      end

      subroutine exp_integral(lobound,hibound,beta,eps,N,ff)
      IMPLICIT NONE

      integer N,i
      real*8 lobound,hibound,beta,eps,h,ff,midpt,x

      h=(hibound-lobound)/N
      ff=0.0 
      do i=1,N
       x=lobound+(i-0.5)*h
       call integrand_function(x,eps,beta,midpt)
       ff=ff+midpt*h
      enddo

      return
      end

      subroutine bisection_function(beta,eps,JA,ff_bisect)
      IMPLICIT NONE

      real*8 beta,eps,JA,hibound,ff_bisect,ff_int
      integer N

      if (beta.le.0.0) then
       print *,"beta cannot be <= 0"
       stop
      endif

      hibound=10.0+beta
      N=1000
      
      call exp_integral(beta,hibound,beta,eps,N,ff_int)
      ff_bisect=JA-ff_int

      return
      end
 
      program main
      IMPLICIT NONE

      character*12 namestr

      real*8 k_liq,k_vap,cp_liq,cp_vap,rho_liq,rho_vap,sigma
      real*8 mu_liq,mu_vap,L_vap,T_sat,Jacob_number,eps,alpha
      real*8 T_infinity
      real*8 beta
      real*8 ff_a,ff_b,ff_c
      real*8 a,b,c,r_start,r_end,t_start,t_end
      integer niter,N,i
      real*8 ti,dt

      rho_liq=958.0
      rho_vap=0.59
      sigma=0.059  ! not used
      mu_liq=2.82E-4
      mu_vap=1.23E-6 ! not used
      k_liq=0.6
      k_vap=0.026  ! not used
      cp_liq=4216.0
      cp_vap=2034.0  ! not used
      L_vap=2.257E+6
      T_sat=373.0
      eps=1.0-rho_vap/rho_liq 
      alpha=k_liq/(rho_liq*cp_liq)
      Jacob_number=3.0
      T_infinity=(Jacob_number*rho_vap*L_vap)/(rho_liq*cp_liq)+T_sat

      print *,"Jacob_number= ",Jacob_number
      print *,"T_infinity= ",T_infinity
      print *,"T_sat= ",T_sat
      print *,"eps=1-rho_vap/rho_liq= ",eps
      beta=1.0
      call bisection_function(beta,eps,Jacob_number,ff_a)
      niter=0
      do while (ff_a.lt.0.0)
       beta=beta/2.0
       call bisection_function(beta,eps,Jacob_number,ff_a)
       niter=niter+1
      enddo
      a=beta
      print *,"niter,a,ff_a ",niter,a,ff_a

      b=beta
      call bisection_function(beta,eps,Jacob_number,ff_b)
      niter=0
      print *,"niter,b ",niter,b,ff_b
      do while (ff_b.gt.0.0) 
       beta=2.0*beta
       call bisection_function(beta,eps,Jacob_number,ff_b)
       b=beta
       niter=niter+1
       print *,"niter,b ",niter,b,ff_b
      enddo
      b=beta
      print *,"niter,b ",niter,b,ff_b
 
      do niter=1,64
       c=(a+b)/2.0
       beta=c
       call bisection_function(beta,eps,Jacob_number,ff_c)
       if (ff_a*ff_c.le.0.0) then
        b=c
        ff_b=ff_c
       else if (ff_c*ff_b.le.0.0) then
        a=c
        ff_a=ff_c
       else
        print *,"bisection failure"
        print *,"niter=",niter
        print *,"a,b,c ",a,b,c
        print *,"fa,fb,fc ",ff_a,ff_b,ff_c
        stop
       endif
      enddo
      print *,"niter,c,ff_c ",niter,c,ff_c
      print *,"R(t)=2 beta sqrt(alpha * t)"
      print *,"R'(t)=alpha * beta/sqrt(alpha * t)"
      print *,"alpha= ",alpha
      print *,"beta= ",beta
      r_start=1.0E-3
      r_end=2.0*r_start
      t_start=(r_start/(2.0*beta))**2/alpha
      t_end=(r_end/(2.0*beta))**2/alpha
      print *,"r_start,r_end ",r_start,r_end
      print *,"t_start,t_end ",t_start,t_end
      print *,"t_end-t_start= ",t_end-t_start
      print *,"velocity at t=t_start: ",alpha*beta/sqrt(alpha*t_start)
      namestr='radius_front'
      open(unit=11,file=namestr)
      N=200
      dt=(t_end-t_start)/N
      do i=0,N-1
       ti=t_start+(i+0.5)*dt
       write(11,*) ti-t_start,2.0*beta*sqrt(alpha*ti)
      enddo
      close(11)
      print *,"interface file is radius_front"

      end

