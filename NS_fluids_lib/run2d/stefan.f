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
     
c 1d stefan problem
c Welch and Wilson, JCP 2000
c lambda exp(lambda^2)erf(lambda)=
c   cp (theta_wall - theta_sat)/(h_{lg}Pi^.5)
      real*8 function transf(x,cp)
      IMPLICIT NONE

      real*8 erfmine
      real*8 x,pi,cp

      pi=4.d0*atan(1.d0)
      transf=x*exp(x*x)*erfmine(x)-cp/sqrt(pi)

      return
      end
 
      program main
      IMPLICIT NONE
      real*8 deltat,cp,hlg,kvapor,denvapor,alpha
      real*8 kliquid,denliquid,psat,tsat
      real*8 transf
      real*8 erfmine
      real*8 a,b,fa,fb,c,fc
      real*8 fact,xlo,xhi,xstop,tstop,dx,dx2,xi,theta,dt
      real*8 CC,beta,xint,current_slope
      real*8 time,normdiff,ax
      integer i,testcase,N,N2
      integer itime,ispace
      real*8 l(5000),u(5000),d(5000),f(5000),soln(5000)
      real*8 old_theta(5000)
      character*12 namestr
      integer is_freezing
      real*8 xint_output,time_output
      integer RZFLAG,simple_parm
      real*8 rplus,rminus,rcen
      real*8 xcfd_offset
      real*8 xcalibrate
      real*8 tcalibrate

      xcfd_offset=0.0d0
      xcalibrate=0.0d0
      tcalibrate=0.0d0

c theta_t=(1/r)(kr theta_r)_r
c r'=r-int_0^t V(s) ds=r-X(t)
c theta(r'(r,t),t)_t=theta_r' r'_t + theta_t
c theta_t=V(t)theta_r'+(1/(r'+X(t)))(k (r'+X(t))theta_r')_r'
      RZFLAG=1
      simple_parm=1 
c testcase=0 1D stefan (analytical solution known, k=0 in source material)
c Welch and Wilson, JCP 2000
c testcase=1 sucking (1d numerical method, k=0 in destination material)
c sucking problem has option of cylindrical coordinate systems.
      testcase=0
c is_freezing=1 => tstop=700
c is_freezing=2 => tstop=0.002 
c is_freezing=3 => tstop=1.75 ? (SEE H. Hu and Z. Jin)
      is_freezing=3

c stefan problem
      if (testcase.eq.0) then
c When running this test, we set kliquid=0
       cp=4.1855E+7
       if (is_freezing.gt.0) then
        hlg=-3.34E+9
        if (is_freezing.eq.3) then
         hlg=-65.0E+9
        endif
       else
        hlg=2.26E+10
       endif
       denliquid=1.0

c 1W=10^7 erg/s
c 1W/(m K)=10^5 cgs
       kliquid=58.0E+3
       if (is_freezing.gt.0) then
        kvapor=218.0E+3
       else
        kvapor=2.4E+3
       endif
       if (is_freezing.eq.3) then
        kliquid=54.6E+3
c ice
        kvapor=220.0E+3
       endif

       if (is_freezing.gt.0) then
        denvapor=0.934
        tsat=273.0
        deltat=-10.0
        if (is_freezing.eq.3) then
         denvapor=0.917d0
         tsat=273.0
c See H. Hu and Z. Jin
c Tfreeze=32.4
c        deltat=-2.0
c Tfreeze=21.5
c        deltat=-3.0
c Tfreeze=16.2
c        deltat=-4.0
c Tfreeze=13.0
         deltat=-5.0
c ice
         cp=2.03E+7
         xcalibrate=1.5625e-3
         xcfd_offset=3.125e-3
        endif
       else
        psat=101.0
        if (psat.eq.101.0) then
         denvapor=denliquid/1605.2
        else if (psat.eq.571.0) then
         denvapor=denliquid/301.4
        else if (psat.eq.14044.0) then
         denvapor=denliquid/7.08
        else
         print *,"psat invalid"
         stop
        endif
        tsat=373.0
        deltat=25.0
       endif

       fact=cp*deltat/hlg
       alpha=kvapor/(denvapor*cp)
 
       a=0.0
       b=2.0
       fb=transf(b,fact)
       do while (fb.le.0.0) 
        b=2.0*b
        fb=transf(b,fact)
       enddo 
       fa=transf(a,fact)
       fb=transf(b,fact)
       print *,"a,b,fa,fb ",a,b,fa,fb
       if (fa*fb.gt.0.0) then
        print *,"a,b invalid"
       endif
       do i=1,100
        c=(a+b)/2.0  
        fc=transf(c,fact)
        if (fa*fc.gt.0.0) then
         a=c
         fa=fc
        else
         b=c
         fb=fc
        endif
        print *,"c,fc ",c,fc
       enddo

       print *,"pos=2 lambda sqrt(alpha t)"
       print *,"lambda= ",c
       print *,"2 lambda= ",2.0*c
       print *,"alpha= ",alpha 

       tcalibrate=(xcalibrate/(2.0*c))**2/alpha
       print *,"xcalibrate ",xcalibrate
       print *,"tcalibrate ",tcalibrate

       if (is_freezing.eq.3) then
        xlo=xcalibrate
c See H. Hu and Z. Jin
        tstop=33.0d0  
        xhi=2.0d0*c*sqrt(alpha*(tcalibrate+tstop))
        xstop=xhi
       else
        xlo=0.0
        xhi=2.0
        xstop=1.0
        tstop=(xstop/(2.0*c))**2/alpha
       endif
       print *,"xlo+xcfd_offset,xhi+xcfd_offset ",  
     &    xlo+xcfd_offset,xhi+xcfd_offset
       print *,"xstop+xcfd_offset,tstop ",xstop+xcfd_offset,tstop
       namestr='temp_profile'
       open(unit=11,file=namestr)
       N=200
       dx=(xhi-xlo)/N
       do i=0,N-1
        xi=xlo+(i+0.5)*dx
        if (xi.ge.xstop) then
         theta=tsat
        else
         theta=tsat+deltat-
     &    (deltat/erfmine(c))* 
     &    erfmine(xi/(2.0*sqrt(alpha*tstop)))
        endif
        write(11,*) xi+xcfd_offset,theta
       enddo
       
       close(11)
       print *,"temperature file is temp_profile"

       print *,"xcfd_offset ",xcfd_offset
       print *,"xcalibrate ",xcalibrate
       print *,"tcalibrate ",tcalibrate
       time=0.0
       open(unit=12,file='interface')
       xint=0.0
       xint=2.0*c*sqrt(alpha*(time+tcalibrate))
       write(12,*) time,xint+xcfd_offset
       itime=0
       do while (time.lt.tstop)
        dt=tstop/256.0
        time=time+dt
        xint=2.0*c*sqrt(alpha*(time+tcalibrate))
        write(12,*) time,xint+xcfd_offset
       enddo
       close(12)  
       print *,"interface file is: interface"

c sucking problem
      else if (testcase.eq.1) then

c When running this test, we ignore kvapor (the destination material)
c 1 WATT=10^7 cm^2 g/s^3= 1 J/s  1 J=1 N M
c 1 Joule=1 N M = kg m/s^2  m =kg m^2/s^2
       cp=4.1855E+7
c 2260 J/g=2260x10^7 erg/g=2.26x10^10 erg/g
       if (is_freezing.gt.0) then
        hlg=-3.34E+9
       else
        hlg=2.26E+10
       endif
       denliquid=1.0
c 0.58 W/(M K)
       kliquid=58.0E+3
c 0.024 W/(M K)
       if (is_freezing.gt.0) then
        kvapor=218.0E+3
       else
        kvapor=2.4E+3
       endif
       if (is_freezing.gt.0) then
        denvapor=0.934
        tsat=273.0
        deltat=-20.0
       else
        psat=101.0
        if (psat.eq.101.0) then
         denvapor=denliquid/1605.2
        else if (psat.eq.571.0) then
         denvapor=denliquid/301.4
        else if (psat.eq.14044.0) then
         denvapor=denliquid/7.08
        else
         print *,"psat invalid"
         stop
        endif
        tsat=373.0
        deltat=5.0
       endif

       if (is_freezing.eq.1) then
        tstop=700.00
       else if (is_freezing.eq.2) then
        tstop=0.002
       else
        tstop=1.0
       endif
       xhi=4.0
       if (is_freezing.eq.2) then
        xhi=0.01
       endif
       xint=0.0
       xlo=0.0

       if (simple_parm.eq.1) then
        cp=1.0
        denliquid=1.0
        denvapor=1.0
        kliquid=1.0
        kvapor=1.0
        hlg=-8.0
        tsat=273.0
        deltat=-1.0 
        tstop=0.1
       endif
 
       alpha=kliquid/(denliquid*cp)
       CC=kliquid/(denvapor*hlg)
       beta=denvapor/denliquid

       print *,"alpha,CC,beta ",alpha,CC,beta

       N=2000
       dx=(xhi-xlo)/N

       time=0.0
       open(unit=11,file='interface')
       write(11,*) time,xint

       do ispace=1,N-1
        old_theta(ispace)=tsat+deltat
        soln(ispace)=old_theta(ispace)
       enddo
       itime=0
       do while (time.lt.tstop)
        current_slope=(soln(1)-tsat)*CC/dx
        if (abs(current_slope).lt.1.0e-10) then
         dt=0.25*dx
        else
         dt=0.25*dx/abs(current_slope)
        endif
        if (time+dt.gt.tstop) then
         dt=tstop-time+1.0E-10
        endif
        do ispace=1,N-1 

         if (RZFLAG.eq.0) then
          rplus=1.0
          rminus=1.0
          rcen=1.0
         else if (RZFLAG.eq.1) then
          rcen=xlo+ispace*dx+xint
          rplus=rcen+0.5*dx
          rminus=rcen-0.5*dx
         else
          print *,"RZFLAG invalid"
          stop
         endif 
         l(ispace)=rminus*alpha/(dx**2)
         u(ispace)=rplus*alpha/(dx**2)+rcen*beta*current_slope/dx
         d(ispace)=-l(ispace)-u(ispace)-rcen/dt
         f(ispace)=-old_theta(ispace)*rcen/dt
         if (ispace.eq.1) then
          f(ispace)=f(ispace)-l(ispace)*tsat
         endif
         if (ispace.eq.N-1) then
          f(ispace)=f(ispace)-u(ispace)*(tsat+deltat)
         endif
        enddo
        call tridiag_solve(l,u,d,N-1,f,soln)
        do ispace=1,N-1 
         old_theta(ispace)=soln(ispace)
        enddo
        normdiff=0.0
        do ispace=1,N-1
          ax=0.0
          if (ispace.gt.1) then
           ax=ax+l(ispace)*soln(ispace-1)
          endif
          if (ispace.lt.N-1) then
           ax=ax+u(ispace)*soln(ispace+1)
          endif
          ax=ax+d(ispace)*soln(ispace)
          normdiff=normdiff+(ax-f(ispace))**2
        enddo
        normdiff=sqrt(normdiff)
        if (1.eq.0) then
         print *,"itime,normdiff ",itime,normdiff
        endif
        current_slope=(soln(1)-tsat)*CC/dx
        print *,"time,velocity ",time+dt,current_slope
        xint=xint+dt*current_slope 
        xint_output=xint
        time_output=time+dt
        if ((is_freezing.eq.1).and.(testcase.eq.1).and. 
     &      (simple_parm.eq.0)) then
         xint_output=xint_output+0.5
        endif
        if ((is_freezing.eq.1).and.(testcase.eq.1).and. 
     &      (simple_parm.eq.1)) then
         time_output=time_output
        endif
        write(11,*) time_output,xint_output

        time=time+dt
       enddo
       close(11)

       print *,"x(t) file is: interface"

       namestr='temp_profile'
       open(unit=11,file=namestr)
       N2=200 
       dx2=(xint-xlo)/N2
       do i=0,N2-1
        xi=xlo+(i+0.5)*dx2
        theta=tsat
        write(11,*) xi,theta
       enddo
       do ispace=0,N
        xi=xint+ispace*dx
        if (ispace.eq.0) then
         theta=tsat
        else if (ispace.eq.N) then
         theta=tsat+deltat
        else
         theta=soln(ispace)
        endif
        write(11,*) xi,theta
       enddo 

       close(11)
       print *,"theta(x) file is: temp_profile"

      else
       print *,"testcase invalid"
       stop
      endif

      end

