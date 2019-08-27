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
       
      subroutine uexact(x,eps,u)
      IMPLICIT NONE
 
      real*8 x,u,eps

      u=1.0-exp(-x/eps) 

      return
      end


      subroutine uexact1prime(x,eps,u)
      IMPLICIT NONE
 
      real*8 x,u,eps

      u=exp(-x/eps)/eps 

      return
      end

      subroutine uexact2prime(x,eps,u)
      IMPLICIT NONE
 
      real*8 x,u,eps

      u=-exp(-x/eps)/(eps*eps) 

      return
      end

      program main
      IMPLICIT NONE
      real*8 erfmine
      real*8 xlo,xhi,dx,pi,eps,dxelem,xloj,xhij
      real*8 dxvol,dxleft,dxright,maxerr,l1err,exactsoln,localerr
      integer i,j,N,M,jbase
      real*8 cells(0:9000),nodes(0:9000)
      real*8 l(9000),u(9000),d(9000),f(9000),soln(9000)
      character*12 namestr

      pi=4.d0*atan(1.d0)
c number of subcells (gauss lobatto points) per element
      N=160  
c number of elements
      M=2

      xlo=0.0
      xhi=1.0
      eps=0.01
      dxelem=(xhi-xlo)/M
      nodes(0)=xlo
      do j=0,M-1
       xloj=xlo+j*dxelem
       xhij=xloj+dxelem
       jbase=j*N
       do i=0,N-1
        if (1.eq.1) then
         cells(jbase+N-1-i)=
     &    0.5*(xhij-xloj)*cos((2.0*i+1.0)*pi/(2*N))+0.5*(xloj+xhij)
         nodes(jbase+N-i)= 
     &    0.5*(xhij-xloj)*cos(i*pi/N)+0.5*(xloj+xhij)
        else
         cells(jbase+i)=xloj+(i+0.5)*dxelem/N
         nodes(jbase+i)=xloj+(i+1.0)*dxelem/N
        endif
       enddo
      enddo

      do j=0,M-1
       jbase=j*N
       do i=0,N-1
        dxvol=nodes(jbase+i+1)-nodes(jbase+i)
        if ((i.eq.0).and.(j.eq.0)) then
         dxleft=cells(jbase+i)-xlo
        else
         dxleft=cells(jbase+i)-cells(jbase+i-1)
        endif
        if ((i.eq.N-1).and.(j.eq.M-1)) then
         dxright=xhi-cells(jbase+i)
        else
         dxright=cells(jbase+i+1)-cells(jbase+i)
        endif
        l(jbase+i+1)=1.0/(dxleft*dxvol)
        u(jbase+i+1)=1.0/(dxright*dxvol)
        d(jbase+i+1)=-l(jbase+i+1)-u(jbase+i+1)
        call uexact2prime(cells(jbase+i),eps,f(jbase+i+1))
        if ((i.eq.0).and.(j.eq.0)) then
         f(jbase+i+1)=f(jbase+i+1)-l(jbase+i+1)*0.0
        endif
        if ((i.eq.N-1).and.(j.eq.M-1)) then
         f(jbase+i+1)=f(jbase+i+1)-u(jbase+i+1)*1.0
        endif
       enddo
      enddo
      call tridiag_solve(l,u,d,M*N,f,soln)

      namestr='temp_profile'
      open(unit=11,file=namestr)
      maxerr=0.0
      l1err=0.0
      do j=0,M-1
       jbase=j*N
       do i=0,N-1
        write(11,*) cells(jbase+i),soln(jbase+i+1)
        call uexact(cells(jbase+i),eps,exactsoln)
        localerr=abs(exactsoln-soln(jbase+i+1))
        if (localerr.gt.maxerr) then
         maxerr=localerr
        endif
        dxvol=nodes(jbase+i+1)-nodes(jbase+i)
        l1err=l1err+localerr*dxvol
       enddo
      enddo

      close(11)
      print *,"soln(x) file is: temp_profile"
      print *,"number of elements M= ",M
      print *,"number of sub elements N= ",N
      print *,"max error= ",maxerr
      print *,"L1 error= ",l1err

      end

