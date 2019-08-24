      subroutine tridiag_solve(l,u,d,n,f,soln)
      IMPLICIT NONE

      integer n,i
      real(kind=8) :: l(n),u(n),d(n),f(n),soln(n)
      real(kind=8) ::  ll(n),uu(n),dd(n),z(n)

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
      end subroutine



      subroutine get_k(x,y,k,probtype)
      IMPLICIT NONE

      REAL*8 x,y
      COMPLEX*16 k
      integer probtype

      if (probtype.eq.0) then
       k=(40.0,0.0)
      else 
       k=(0.0,0.0)
      endif

      return
      end     


      subroutine get_alpha(x,y,alpha,probtype)
      IMPLICIT NONE

      COMPLEX*16 alpha
      REAL*8 x,y,dx,dy,dist,buffersize
      integer probtype

      buffersize=0.0
      if (probtype.eq.0) then
       alpha=(0.0,0.0)
       if (x.gt.0.5) then
        dx=1.0-x
       else
        dx=x
       endif
       if (y.gt.0.5) then
        dy=1.0-y
       else
        dy=y
       endif
       dist=min(dx,dy) 
       if (dist.lt.buffersize) then
        alpha=cmplx(0.25*(buffersize-dist)/buffersize,0.0)
       endif
      else 
       alpha=(0.0,0.0)
      endif

      return
      end     


      subroutine get_G(x,y,h,G,probtype)
      IMPLICIT NONE
       
      REAL*8 x,y,h,dist
      COMPLEX*16 G
      integer probtype
       
      dist=sqrt( (x-0.5)**2 + (y-0.5)**2 )
      if (dist.le.h) then
       if (probtype.eq.0) then
        G=(1.0,0.0)
       else
        G=(1.0,2.0)
       endif
      else
       G=(0.0,0.0)
      endif 

      return
      end     

    
      subroutine WALLBC(x,y,UWALL,UIN,coeff1,coeff2, &
        h,hflag,probtype,gridtype,dir,side)
      IMPLICIT NONE

      integer dir,side
      REAL*8 x,y,h,eps,slope_sign
      integer hflag,probtype,gridtype
      COMPLEX*16 UWALL,UIN,coeff1,coeff2
      COMPLEX*16 k
      COMPLEX*16 A,B,C,jhat
  
       ! in spherical coordinates:
       ! if p=u exp(-i w t), then
       !  p=exp(i(kr-wt))/r + exp(i(-kr-wt))/r
       !  disregard the latter term:
       !  p_r=ikp-p/r  p_r - ikp = -p/r -> 0  p_r - ikp=0.
       ! if p=u exp(i w t), then
       !  p=exp(i(kr+wt))/r + exp(i(-kr+wt))/r
       !  disregard the former term:
       !  p_r=-ikp-p/r p_r + ikp = -p/r -> 0  p_r + ikp=0.

       ! this would be -1.0 if we had p=u exp(i w t) instead of
       ! p=u exp(-i w t)
      slope_sign=1.0
       ! gridtype=0:
       ! (OUT-IN)/h - jhat k (OUT+IN)/2 = -u/r
       ! A * OUT + B * IN = C
       ! gridtype=1:
       ! (OUT-IN)/h - jhat k (eps OUT+IN)/(1+eps) = -u/r
       ! A * OUT + B * IN = C
      jhat=(0.0,1.0)
      call get_k(x,y,k,probtype)

      if ((dir.eq.1).or.(dir.eq.2)) then
       if ((side.eq.1).or.(side.eq.2)) then
        if (gridtype.eq.0) then
         eps=1.0
        else if (gridtype.eq.1) then
         eps=1.0E-10
        else
         print *,"gridtype invalid"
         stop
        endif

        if (probtype.eq.0) then
         A=cmplx(slope_sign/h,0.0)-jhat*k*eps/(1.0+eps)
         B=cmplx(-slope_sign/h,0.0)-jhat*k/(1.0+eps)
         C=(0.0,0.0)
        else
         A=eps/(1.0+eps)
         B=1.0/(1.0+eps)
         C=(0.0,0.0)
        endif
        coeff1=C/A
        coeff2=-B/A
        UWALL=coeff1+coeff2*UIN
       else
        print *,"side invalid"
        stop
       endif
      else
       print *,"dir invalid"
       stop
      endif

      return
      end
 
      subroutine set_boundary(nx,ny,lox,loy,hix,hiy, &
        U,h,hflag,probtype,gridtype)
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,dir,side,gridtype
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      REAL*8 h,x,y
      COMPLEX*16 coeff1,coeff2
      integer hflag,probtype
      integer ilo,ihi,jlo,jhi,ii,jj

      do dir=1,2
      do side=1,2
        ! xlo
       if ((dir.eq.1).and.(side.eq.1)) then
        ilo=lox-1
        ihi=lox-1
        jlo=loy
        jhi=hiy
        ii=1
        jj=0
        ! xhi
       else if ((dir.eq.1).and.(side.eq.2)) then
        ilo=hix+1
        ihi=hix+1
        jlo=loy
        jhi=hiy
        ii=-1
        jj=0
        ! jlo
       else if ((dir.eq.2).and.(side.eq.1)) then
        jlo=loy-1
        jhi=loy-1
        ilo=lox
        ihi=hix
        ii=0
        jj=1
        ! jhi
       else if ((dir.eq.2).and.(side.eq.2)) then
        jlo=hiy+1
        jhi=hiy+1
        ilo=lox
        ihi=hix
        ii=0
        jj=-1
       else
        print *,"dir or side invalid"
        stop
       endif

       do i=ilo,ihi
       do j=jlo,jhi
        if (gridtype.eq.0) then
         x=(i+0.5)*h+0.5*ii*h
         y=(j+0.5)*h+0.5*jj*h 
        else if (gridtype.eq.1) then
         x=i*h
         y=j*h 
        endif
        call WALLBC(x,y,U(i,j),U(i+ii,j+jj),coeff1,coeff2, &
         h,hflag,probtype,gridtype,dir,side)
       enddo 
       enddo 
      enddo 
      enddo 

      return
      end

      subroutine DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,U,V,dsum)
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,gridtype
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 V(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 dsum

      dsum=(0.0,0.0)
      do i=lox,hix
      do j=loy,hiy
       dsum=dsum+dconjg(U(i,j))*V(i,j)
      enddo
      enddo 
  
      return
      end

      subroutine NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,U,dnorm)
      IMPLICIT NONE

      integer nx,ny,lox,loy,hix,hiy,gridtype
      REAL*8 dnorm
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 dsum

      call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,U,U,dsum) 
      dnorm=sqrt(REAL(dsum)**2+IMAG(dsum)**2) 
      dnorm=sqrt(dnorm)

      return
      end

       ! W=alpha*U + beta*V
      subroutine LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,U,V,W,alpha,beta) 
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,gridtype
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 V(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 W(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 alpha,beta
      do i=lox,hix
      do j=loy,hiy
       W(i,j)=alpha*U(i,j)+beta*V(i,j)
      enddo
      enddo

      return
      end

        ! V=U
      subroutine COPYVEC(nx,ny,lox,loy,hix,hiy,U,V) 
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 V(lox-1:hix+1,loy-1:hiy+1)
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       V(i,j)=U(i,j)
      enddo
      enddo

      return
      end

      subroutine ZAPVEC(nx,ny,lox,loy,hix,hiy,U) 
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       U(i,j)=(0.0,0.0)
      enddo
      enddo

      return
      end

        ! shift: -laplacian-(beta1-jhat beta2)k^2
        ! p_tt + eps p_t - c^2 lap p = 0
        ! if p=u e^{-jhat w t},
        ! -w^2 u - eps jhat w u - c^2 lap u =0
        ! -k^2 u - eps1 k^2 jhat u - lap u = 0
        ! (-lap - (1+jhat eps1)k^{2})u=0  beta2=-eps1  beta1=1
        ! if p=u e^{jhat w t},
        ! -w^2 u + eps jhat w u - c^2 lap u =0
        ! -k^2 u + eps1 k^2 jhat u - lap u = 0
        ! (-lap - (1-jhat eps1)k^{2})u=0  beta2=eps1  beta1=1
      subroutine JACPRECOND(nx,ny,lox,loy,hix,hiy,Z,R,DMINUS,KSQR, &
       h,hflag,probtype,gridtype,do_shift,beta1,beta2,nsmooth)
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,nsmooth,iter,gridtype
      COMPLEX*16 R(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 Z(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 DMINUS(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 KSQR(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 DIAGCOEFF
      REAL*8 h,x,y,weight
      integer hflag,probtype
      COMPLEX*16 OFFDIAG
      integer dir,side,do_shift
      COMPLEX*16 coeff1,coeff2,beta1,beta2,DMINUSTERM,jhat

      COMPLEX*16, dimension(:,:), allocatable :: ZN
      COMPLEX*16, dimension(:,:), allocatable :: ZNP1
      COMPLEX*16, dimension(:,:), allocatable :: ZSTAR

      jhat=(0.0,1.0)

      weight=0.5

      allocate(ZN(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(ZNP1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(ZSTAR(lox-1:hix+1,loy-1:hiy+1)) 

       ! ZN=Z
      call COPYVEC(nx,ny,lox,loy,hix,hiy,Z,ZN)

      do iter=1,nsmooth
       call set_boundary(nx,ny,lox,loy,hix,hiy, &
        ZN,h,hflag,probtype,gridtype)
       do i=lox,hix
       do j=loy,hiy
        if (gridtype.eq.0) then
         x=(i+0.5)*h
         y=(j+0.5)*h
        else
         x=i*h
         y=j*h
        endif
        DIAGCOEFF=dcmplx(4.0/(h*h),0.0)
        OFFDIAG=(0.0,0.0)
        if (i.gt.lox) then
         OFFDIAG=OFFDIAG+ZN(i-1,j)/(h*h)
        else
         dir=1
         side=1
         x=0.0
         call WALLBC(x,y,ZN(i-1,j),ZN(i,j),coeff1,coeff2, &
          h,hflag,probtype,gridtype,dir,side)
         DIAGCOEFF=DIAGCOEFF-coeff2/(h*h)
        endif
        if (i.lt.hix) then
         OFFDIAG=OFFDIAG+ZN(i+1,j)/(h*h)
        else
         dir=1
         side=2
         x=nx*h
         call WALLBC(x,y,ZN(i+1,j),ZN(i,j),coeff1,coeff2, &
          h,hflag,probtype,gridtype,dir,side)
         DIAGCOEFF=DIAGCOEFF-coeff2/(h*h)
        endif
        if (j.gt.loy) then
         OFFDIAG=OFFDIAG+ZN(i,j-1)/(h*h)
        else
         dir=2
         side=1
         y=0.0
         call WALLBC(x,y,ZN(i,j-1),ZN(i,j),coeff1,coeff2, &
          h,hflag,probtype,gridtype,dir,side)
         DIAGCOEFF=DIAGCOEFF-coeff2/(h*h)
        endif
        if (j.lt.hiy) then
         OFFDIAG=OFFDIAG+ZN(i,j+1)/(h*h)
        else
         dir=2
         side=2
         y=ny*h
         call WALLBC(x,y,ZN(i,j+1),ZN(i,j),coeff1,coeff2, &
          h,hflag,probtype,gridtype,dir,side)
         DIAGCOEFF=DIAGCOEFF-coeff2/(h*h)
        endif
        if (do_shift.eq.0) then
         DMINUSTERM=DMINUS(i,j)
        else if (do_shift.eq.1) then
         DMINUSTERM=(beta1-jhat*beta2)*KSQR(i,j)
        else
         print *,"do_shift invalid"
         stop
        endif
        DIAGCOEFF=DIAGCOEFF-DMINUSTERM

        ZSTAR(i,j)=(1.0/DIAGCOEFF)*(R(i,j)+OFFDIAG)
        ZNP1(i,j)=weight*ZSTAR(i,j)+(1.0-weight)*ZN(i,j)
       enddo
       enddo
        ! ZN=ZNP1
       call COPYVEC(nx,ny,lox,loy,hix,hiy,ZNP1,ZN)
      enddo  ! iter
        ! Z=ZN
      call COPYVEC(nx,ny,lox,loy,hix,hiy,ZN,Z)

      deallocate(ZN) 
      deallocate(ZNP1) 
      deallocate(ZSTAR) 

      return
      end


      subroutine ATIMESU(nx,ny,lox,loy,hix,hiy, &
       AU,U,DMINUS,KSQR,h,hflag,probtype, &
       gridtype,do_shift,beta1,beta2)
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,gridtype
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 AU(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 DMINUS(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 KSQR(lox-1:hix+1,loy-1:hiy+1)
      REAL*8 h
      integer hflag,probtype,do_shift
      COMPLEX*16 beta1,beta2,jhat,DMINUSTERM

      jhat=(0.0,1.0)
      call set_boundary(nx,ny,lox,loy,hix,hiy, &
       U,h,hflag,probtype,gridtype)
      do i=lox,hix
      do j=loy,hiy
       if (do_shift.eq.0) then
        DMINUSTERM=DMINUS(i,j)
       else if (do_shift.eq.1) then
        DMINUSTERM=(beta1-jhat*beta2)*KSQR(i,j)
       else
        print *,"do_shift invalid"
        stop
       endif
 
       AU(i,j)=-(U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1)-4.0*U(i,j))/(h*h)- &
         DMINUSTERM*U(i,j)
      enddo
      enddo
        
      return
      end

        ! R=G-A U
      subroutine RESID(nx,ny,lox,loy,hix,hiy, &
        R,U,G,DMINUS,KSQR,h,hflag,probtype, &
        gridtype,do_shift,beta1,beta2)
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,gridtype
      COMPLEX*16 R(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 G(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 DMINUS(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 KSQR(lox-1:hix+1,loy-1:hiy+1)
      REAL*8 h,x,y
      integer hflag,probtype,do_shift
      COMPLEX*16, dimension(:,:), allocatable :: AU
      COMPLEX*16 :: alpha,beta,beta1,beta2

      allocate(AU(lox-1:hix+1,loy-1:hiy+1)) 

      call ATIMESU(nx,ny,lox,loy,hix,hiy, &
       AU,U,DMINUS,KSQR,h,hflag,probtype, &
       gridtype,do_shift,beta1,beta2)
      alpha=(1.0,0.0)
      beta=(-1.0,0.0)
      call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,G,AU,R,alpha,beta)

      deallocate(AU)

      return
      end

      subroutine checkzero(val,idx)
      IMPLICIT NONE

      COMPLEX*16 val
      integer idx

      if (sqrt(real(val)**2+imag(val)**2).lt.1.0E-14) then
       print *,"value is almost 0"
       print *,"val = ",val
       print *,"idx = ",idx
       stop
      endif

      return
      end

      Real*8 function mag_complex(x)
      IMPLICIT NONE

      COMPLEX*16 x
   
      mag_complex=sqrt(REAL(x)**2+IMAG(x)**2)

      return
      end

      subroutine RESTRICT(nx,ny,lox,loy,hix,hiy,nxC,nyC, &
        loxC,loyC,hixC,hiyC,gridtype,RF,RC)
      IMPLICIT NONE

      integer loxC,loyC,hixC,hiyC
      integer nx,ny,lox,loy,hix,hiy,nxC,nyC,gridtype,i,j,if,jf
      COMPLEX*16, dimension(lox-1:hix+1,loy-1:hiy+1) :: RF
      COMPLEX*16, dimension(loxC-1:hixC+1,loyC-1:hiyC+1) :: RC

      do i=loxC,hixC
      do j=loyC,hiyC
       if=2*i
       jf=2*j
       if (gridtype.eq.0) then
        RC(i,j)=(RF(if,jf)+RF(if+1,jf)+ &
                 RF(if,jf+1)+RF(if+1,jf+1))/4.0
       else if (gridtype.eq.1) then
        if ((i.eq.loxC).or.(i.eq.hixC).or. &
            (j.eq.loyC).or.(j.eq.hiyC)) then
         RC(i,j)=RF(if,jf)
        else
         RC(i,j)=( &
          4.0*RF(if,jf)+ &
          2.0*RF(if+1,jf)+ &
          2.0*RF(if-1,jf)+ &
          2.0*RF(if,jf+1)+ &
          2.0*RF(if,jf-1)+ &
          1.0*RF(if+1,jf+1)+ &
          1.0*RF(if-1,jf+1)+ &
          1.0*RF(if+1,jf-1)+ &
          1.0*RF(if-1,jf-1))/16.0
        endif
       else
        print *,"gridtype invalid"
        stop
       endif

      enddo
      enddo

      return
      end

      recursive subroutine RELAX(nx,ny,lox,loy,hix,hiy,Z,R,DMINUS,KSQR, &
       h,hflag,probtype,gridtype,do_shift,beta1,beta2,nsmooth)
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,irelax,nsmooth,gridtype
      COMPLEX*16 R(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 Z(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 DMINUS(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 KSQR(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 DIAGCOEFF
      integer hflag,probtype
      integer do_shift
      COMPLEX*16 beta1,beta2
      integer if,jf,nxC,nyC,loxC,loyC,hixC,hiyC
      REAL*8 h,hC
      integer count_relax,coarsest_nx

      COMPLEX*16, dimension(:,:), allocatable :: ZN
      COMPLEX*16, dimension(:,:), allocatable :: RF
      COMPLEX*16, dimension(:,:), allocatable :: RC
      COMPLEX*16, dimension(:,:), allocatable :: ZC
      COMPLEX*16, dimension(:,:), allocatable :: DMINUSC
      COMPLEX*16, dimension(:,:), allocatable :: KSQRC


      count_relax=2
      coarsest_nx=1

      allocate(ZN(lox-1:hix+1,loy-1:hiy+1)) 

      if (hflag.ne.1) then
       print *,"hflag invalid"
       stop
      endif

      call COPYVEC(nx,ny,lox,loy,hix,hiy,Z,ZN)

      call JACPRECOND(nx,ny,lox,loy,hix,hiy, &
        ZN,R,DMINUS,KSQR,h,hflag,probtype, &
        gridtype,do_shift,beta1,beta2,nsmooth)

      nxC=nx/2
      nyC=ny/2
      loxC=lox/2
      loyC=loy/2
      if (gridtype.eq.0) then
       hixC=nxC-1
       hiyC=nyC-1
      else if (gridtype.eq.1) then
       hixC=nxC
       hiyC=nyC
      else
       print *,"gridtype invalid"
       stop
      endif 
      if ((2*nxC.eq.nx).and.(2*nyC.eq.ny).and.(nxC.ge.coarsest_nx)) then
       hC=2.0*h      
       allocate(RF(lox-1:hix+1,loy-1:hiy+1)) 
        ! RF=R-A ZN
       call RESID(nx,ny,lox,loy,hix,hiy, &
        RF,ZN,R,DMINUS,KSQR,h,hflag,probtype, &
        gridtype,do_shift,beta1,beta2)

       allocate(ZC(loxC-1:hixC+1,loyC-1:hiyC+1)) 
       allocate(RC(loxC-1:hixC+1,loyC-1:hiyC+1)) 
       allocate(DMINUSC(loxC-1:hixC+1,loyC-1:hiyC+1)) 
       allocate(KSQRC(loxC-1:hixC+1,loyC-1:hiyC+1)) 
       call ZAPVEC(nxC,nyC,loxC,loyC,hixC,hiyC,ZC)
       call RESTRICT(nx,ny,lox,loy,hix,hiy,nxC,nyC, &
         loxC,loyC,hixC,hiyC,gridtype,RF,RC)
       call RESTRICT(nx,ny,lox,loy,hix,hiy,nxC,nyC, &
         loxC,loyC,hixC,hiyC,gridtype,DMINUS,DMINUSC)
       call RESTRICT(nx,ny,lox,loy,hix,hiy,nxC,nyC, &
         loxC,loyC,hixC,hiyC,gridtype,KSQR,KSQRC)

       do irelax=1,count_relax
        call RELAX(nxC,nyC,loxC,loyC,hixC,hiyC,ZC,RC,DMINUSC,KSQRC, &
         hC,hflag,probtype,gridtype,do_shift,beta1,beta2,nsmooth)
       enddo
       call set_boundary(nxC,nyC,loxC,loyC,hixC,hiyC, &
        ZC,hC,hflag,probtype,gridtype)
       if (gridtype.eq.0) then
        do i=loxC,hixC
        do j=loyC,hiyC
         if=2*i
         jf=2*j
         ZN(if,jf)=ZN(if,jf)+ZC(i,j)
         ZN(if+1,jf)=ZN(if+1,jf)+ZC(i,j)
         ZN(if+1,jf+1)=ZN(if+1,jf+1)+ZC(i,j)
         ZN(if,jf+1)=ZN(if,jf+1)+ZC(i,j)
        enddo
        enddo
       else if (gridtype.eq.1) then
        do i=loxC+1,hixC
        do j=loyC+1,hiyC
         if=2*i
         jf=2*j
         ZN(if,jf)=ZN(if,jf)+ZC(i,j)
         ZN(if-1,jf)=ZN(if-1,jf)+ &
          0.5*(ZC(i,j)+ZC(i-1,j))
         ZN(if,jf-1)=ZN(if,jf-1)+ &
          0.5*(ZC(i,j)+ZC(i,j-1))
         ZN(if-1,jf-1)=ZN(if-1,jf-1)+ &
          0.25*(ZC(i,j)+ZC(i-1,j)+ZC(i,j-1)+ZC(i-1,j-1))
        enddo
        enddo
       else
        print *,"gridtype invalid"
        stop
       endif

       deallocate(ZC)
       deallocate(RC)
       deallocate(RF)
       deallocate(DMINUSC)
       deallocate(KSQRC)
      endif
      call JACPRECOND(nx,ny,lox,loy,hix,hiy, &
        ZN,R,DMINUS,KSQR,h,hflag,probtype, &
        gridtype,do_shift,beta1,beta2,nsmooth)

      call COPYVEC(nx,ny,lox,loy,hix,hiy,ZN,Z)

      deallocate(ZN) 

      return
      end

      subroutine preconditioner( &
       nx,ny,lox,loy,hix,hiy,Z,R, &
       DMINUS,KSQR,h,hflag,probtype, &
       gridtype,do_shift,beta1,beta2,nsmooth, &
       precond_type)
      IMPLICIT NONE

      integer nx,ny,lox,loy,hix,hiy,nsmooth,gridtype,precond_type
      COMPLEX*16 R(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 Z(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 DMINUS(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 KSQR(lox-1:hix+1,loy-1:hiy+1)
      integer hflag,probtype
      integer do_shift
      COMPLEX*16 beta1,beta2
      REAL*8 h

      call ZAPVEC(nx,ny,lox,loy,hix,hiy,Z)
      if (precond_type.eq.0) then
       call COPYVEC(nx,ny,lox,loy,hix,hiy,R,Z)
      else if (precond_type.eq.1) then
       call JACPRECOND(nx,ny,lox,loy,hix,hiy, &
        Z,R,DMINUS,KSQR,h,hflag,probtype, &
        gridtype,do_shift,beta1,beta2,nsmooth)
      else if (precond_type.eq.2) then
       call RELAX(nx,ny,lox,loy,hix,hiy, &
        Z,R,DMINUS,KSQR,h,hflag,probtype, &
        gridtype,do_shift,beta1,beta2,nsmooth)
      else
       print *,"precond_type invalid" 
       stop
      endif

      return
      end

        ! DMINUS contains (1-alpha)k^2
        ! precond_type=0 M=I, =1 Jacobi, =2 MG
      subroutine bicgstab(nx,ny,lox,loy,hix,hiy, &
         U,G,DMINUS,KSQR,h,hflag, &
         precond_type,probtype,gridtype,tol)
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,hflag,precond_type,probtype,iter
      integer maxiter,gridtype
      integer hflagcg,restart_flag
      COMPLEX*16 U(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 G(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 DMINUS(lox-1:hix+1,loy-1:hiy+1)
      COMPLEX*16 KSQR(lox-1:hix+1,loy-1:hiy+1)
      REAL*8 h,tol,dnorm,dnorm0,mag_complex

      COMPLEX*16, dimension(:,:), allocatable :: U0
      COMPLEX*16, dimension(:,:), allocatable :: V0
      COMPLEX*16, dimension(:,:), allocatable :: P0
      COMPLEX*16, dimension(:,:), allocatable :: R0
      COMPLEX*16, dimension(:,:), allocatable :: U1
      COMPLEX*16, dimension(:,:), allocatable :: V1
      COMPLEX*16, dimension(:,:), allocatable :: P1
      COMPLEX*16, dimension(:,:), allocatable :: R1
      COMPLEX*16, dimension(:,:), allocatable :: R0hat
      COMPLEX*16, dimension(:,:), allocatable :: UINIT
      COMPLEX*16, dimension(:,:), allocatable :: RHS
      COMPLEX*16, dimension(:,:), allocatable :: Y
      COMPLEX*16, dimension(:,:), allocatable :: Hvec
      COMPLEX*16, dimension(:,:), allocatable :: S
      COMPLEX*16, dimension(:,:), allocatable :: T
      COMPLEX*16, dimension(:,:), allocatable :: Z

      COMPLEX*16 rho0,w0,rho1,alpha,w1,beta,a1,a2
      COMPLEX*16 beta1,beta2
      integer do_shift,nsmooth

      if (precond_type.eq.1) then
       nsmooth=4
      else if (precond_type.eq.2) then
       nsmooth=2
      else
       print *,"precond_type invalid"
       stop
      endif

      if (probtype.eq.0) then
       do_shift=1
       beta1=(1.0,0.0)
       beta2=(-0.5,0.0)
      else
       do_shift=0
       beta1=(0.0,0.0)
       beta2=(0.0,0.0)
      endif

      allocate(U0(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(V0(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(P0(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(R0(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(U1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(V1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(P1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(R1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(R0hat(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(UINIT(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(RHS(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(Y(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(Hvec(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(S(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(T(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(Z(lox-1:hix+1,loy-1:hiy+1)) 

        ! U0=V0=P0=0 
      alpha=(0.0,0.0) 
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       U0(i,j)=alpha
       V0(i,j)=alpha
       P0(i,j)=alpha
      enddo
      enddo

       ! R0=G-A U0
      call RESID(nx,ny,lox,loy,hix,hiy, &
        R0,U0,G,DMINUS,KSQR,h,hflag,probtype, &
        gridtype,0,beta1,beta2) 
       ! R0hat=R0
      call COPYVEC(nx,ny,lox,loy,hix,hiy,R0,R0hat)
       ! UINIT=U0
      call COPYVEC(nx,ny,lox,loy,hix,hiy,U0,UINIT)
      call ZAPVEC(nx,ny,lox,loy,hix,hiy,U0)
       ! RHS=R0
      call COPYVEC(nx,ny,lox,loy,hix,hiy,R0,RHS)

       ! rho0=alpha=w0=1
      rho0=(1.0,0.0)
      alpha=(1.0,0.0)
      w0=(1.0,0.0)

      call NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,R0,dnorm0)
      dnorm=1.0
      maxiter=2000
      iter=0

      hflagcg=1
      do while ((dnorm.gt.tol).and.(iter.lt.maxiter))
       print *,"iter,dnorm ",iter,dnorm

         ! rho1= R0hat^H R0
       call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,R0hat,R0,rho1)
       
       restart_flag=0
       if ((sqrt(mag_complex(rho0)).lt.tol*0.01).or. &
           (sqrt(mag_complex(w0)).lt.tol*0.01)) then
        restart_flag=1
       endif 

       if (restart_flag.eq.0) then
          ! (R0hat^H R0)/(R0hat^H dot R0_before)  *   (alpha/w0)
        beta=(rho1/rho0)*(alpha/w0)
        a1=(1.0,0.0)
        a2=-w0
 
         ! P1=P0-w0 V0
        call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,P0,V0,P1,a1,a2)
         ! P1=R0+beta P1
        call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,R0,P1,P1,a1,beta)
        ! Y=M^{-1}P1
        call preconditioner( &
         nx,ny,lox,loy,hix,hiy,Y,P1, &
         DMINUS,KSQR,h,hflagcg,probtype, &
         gridtype,do_shift,beta1,beta2,nsmooth, &
         precond_type)
         
         ! V1=A Y
        call ATIMESU(nx,ny,lox,loy,hix,hiy, &
         V1,Y,DMINUS,KSQR,h,hflagcg,probtype, &
         gridtype,0,beta1,beta2)

         ! alpha=rho1/R0hat dot V1
        call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,R0hat,V1,alpha)

        if (sqrt(mag_complex(alpha)).lt.tol*0.01) then
         restart_flag=1
        endif

        if (restart_flag.eq.0) then
         alpha=rho1/alpha

         ! Hvec=U0+alpha Y
         a1=(1.0,0.0)
         a2=alpha
         call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,U0,Y,Hvec,a1,a2) 
         ! U1=Hvec
         call COPYVEC(nx,ny,lox,loy,hix,hiy,Hvec,U1)
         ! R1=RHS-A U1
         call RESID(nx,ny,lox,loy,hix,hiy, &
          R1,U1,RHS,DMINUS,KSQR,h,hflagcg,probtype, &
          gridtype,0,beta1,beta2) 
         call NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,R1,dnorm)
         dnorm=dnorm/dnorm0

         if (dnorm.gt.tol) then
          ! S=R0-alpha V1
          a1=(1.0,0.0)
          a2=-alpha
          call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,R0,V1,S,a1,a2) 

           ! Z=M^{-1}S
          call preconditioner( &
           nx,ny,lox,loy,hix,hiy,Z,S, &
           DMINUS,KSQR,h,hflagcg,probtype, &
           gridtype,do_shift,beta1,beta2,nsmooth, &
           precond_type)

           ! T=A Z
          call ATIMESU(nx,ny,lox,loy,hix,hiy, &
            T,Z,DMINUS,KSQR,h,hflagcg,probtype, &
            gridtype,0,beta1,beta2)

           ! simple case is: (T,S)/(T,T)=(AZ,S)/(AZ,AZ)   (MZ=S)
          call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,T,S,a1)
          call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,T,T,a2)
          if (sqrt(mag_complex(a2)).lt.tol*0.01) then
           restart_flag=1
          endif

          if (restart_flag.eq.0) then
           w1=a1/a2
           ! U1=Hvec+w1 Z
           a1=(1.0,0.0)
           a2=w1
           call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,Hvec,Z,U1,a1,a2) 
          endif
         endif ! dnorm>tol
          ! R1=RHS-A U1
         call RESID(nx,ny,lox,loy,hix,hiy, &
          R1,U1,RHS,DMINUS,KSQR,h,hflagcg,probtype, &
          gridtype,0,beta1,beta2) 
         call NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,R1,dnorm)
         dnorm=dnorm/dnorm0
         rho0=rho1
         w0=w1
          ! R0=R1
         call COPYVEC(nx,ny,lox,loy,hix,hiy,R1,R0) 
         call COPYVEC(nx,ny,lox,loy,hix,hiy,P1,P0) 
         call COPYVEC(nx,ny,lox,loy,hix,hiy,V1,V0) 
         call COPYVEC(nx,ny,lox,loy,hix,hiy,U1,U0) 
        endif  ! restart_flag=0
       endif ! restart_flag=0

       if (restart_flag.eq.0) then
        ! do nothing
       else if (restart_flag.eq.1) then
        call RESID(nx,ny,lox,loy,hix,hiy, &
          R0,U0,RHS,DMINUS,KSQR,h,hflagcg,probtype, &
          gridtype,0,beta1,beta2) 
         ! R0hat=R0
        call COPYVEC(nx,ny,lox,loy,hix,hiy,R0,R0hat)
        call COPYVEC(nx,ny,lox,loy,hix,hiy,U0,U1) 
         ! rho0=alpha=w0=1
        rho0=(1.0,0.0)
        alpha=(1.0,0.0)
        w0=(1.0,0.0)
        call NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,R0,dnorm)
        dnorm=dnorm/dnorm0
        call ZAPVEC(nx,ny,lox,loy,hix,hiy,V0)
        call ZAPVEC(nx,ny,lox,loy,hix,hiy,P0)
       else
        print *,"restart_flag invalid"
        stop
       endif

       iter=iter+1
      enddo
      print *,"at the end: iter,dnorm ",iter,dnorm
       ! U=UINIT+U1
      a1=(1.0,0.0)
      a2=a1
      call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,UINIT,U1,U,a1,a2)
 
      deallocate(U0) 
      deallocate(V0) 
      deallocate(P0) 
      deallocate(R0) 
      deallocate(U1) 
      deallocate(V1) 
      deallocate(P1) 
      deallocate(R1) 
      deallocate(R0hat) 
      deallocate(UINIT) 
      deallocate(RHS) 
      deallocate(Y) 
      deallocate(Hvec) 
      deallocate(S) 
      deallocate(T) 
      deallocate(Z) 

      return
      end


      program main
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,hflag
      REAL*8 h,x,y
      COMPLEX*16 alpha,k
      COMPLEX*16 Gpoint
      integer precond_type,probtype,gridtype
      REAL*8 tol
      COMPLEX*16, dimension(:,:), allocatable :: U
      COMPLEX*16, dimension(:,:), allocatable :: G
      COMPLEX*16, dimension(:,:), allocatable :: DMINUS
      COMPLEX*16, dimension(:,:), allocatable :: KSQR
      character*8 wavedatafile

      probtype=0
      gridtype=0  ! 0=CELL  1=NODE
      nx=64
      ny=64
      h=1.0/nx
      tol=1.0E-7
      precond_type=2 ! 0 M=I  1=Jacobi precond.  2=multigrid precond.
      hflag=0
 
      if (gridtype.eq.0) then
       lox=0
       loy=0
       hix=nx-1
       hiy=ny-1
      else if (gridtype.eq.1) then
       lox=0
       loy=0
       hix=nx
       hiy=ny
      else
       print *,"gridtype invalid"
       stop
      endif

      allocate(U(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(G(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(DMINUS(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(KSQR(lox-1:hix+1,loy-1:hiy+1)) 
     
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       if (gridtype.eq.0) then
        x=(i+0.5)*h
        y=(j+0.5)*h
       else
        x=i*h
        y=j*h
       endif
       call get_k(x,y,k,probtype)
       call get_alpha(x,y,alpha,probtype)
       DMINUS(i,j)=(cmplx(1.0,0.0)-alpha)*k*k
       KSQR(i,j)=k*k
       U(i,j)=(0.0,0.0) 
       call get_G(x,y,h,Gpoint,probtype)
       G(i,j)=Gpoint
      enddo
      enddo
      call bicgstab(nx,ny,lox,loy,hix,hiy,U,G,DMINUS,KSQR,h,hflag, &
        precond_type,probtype,gridtype,tol)
      call set_boundary(nx,ny,lox,loy,hix,hiy,U,h,0,probtype,gridtype)

      write(wavedatafile,'(A8)') 'wavedata'
      print *,"wavedatafile ",wavedatafile
      open(unit=11,file=wavedatafile)
      if (1.eq.0) then
       write(11,*) '# X,Y,REALU,IMAGU'
      else 
       write(11,*) 'VARIABLES="X","Y","REALU","IMAGU"'
       write(11,*) 'zone i=',hix-lox+3,' j=',hiy-loy+3,' f=point'
      endif
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       if (gridtype.eq.0) then
        x=(i+0.5)*h
        y=(j+0.5)*h
       else if (gridtype.eq.1) then
        x=i*h
        y=j*h
       endif
       write(11,*) x,y,REAL(U(i,j)),IMAG(U(i,j))
      enddo
      enddo

      close(11)
 
      deallocate(U) 
      deallocate(G) 
      deallocate(DMINUS) 
      deallocate(KSQR) 
 
      return
      end









      ! DMINUS contains (1-alpha)k^2
        ! precond_type=0 M=I, =1 Jacobi, =2 MG
      subroutine bicgstab(nx,ny,lox,loy,hix,hiy, &
         U,G,DMINUS,KSQR,h,hflag, &
         precond_type,probtype,gridtype,tol)
      IMPLICIT NONE

      integer i,j,nx,ny,lox,loy,hix,hiy,hflag,precond_type,probtype,iter
      integer maxiter,gridtype
      integer hflagcg,restart_flag
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1)
      REAL*8 G(lox-1:hix+1,loy-1:hiy+1)
      REAL*8 DMINUS(lox-1:hix+1,loy-1:hiy+1)
      REAL*8 KSQR(lox-1:hix+1,loy-1:hiy+1)
      REAL*8 h,tol,dnorm,dnorm0,mag_complex

      REAL*8, dimension(:,:), allocatable :: U0
      REAL*8, dimension(:,:), allocatable :: V0
      REAL*8, dimension(:,:), allocatable :: P0
      REAL*8, dimension(:,:), allocatable :: R0
      REAL*8, dimension(:,:), allocatable :: U1
      REAL*8, dimension(:,:), allocatable :: V1
      REAL*8, dimension(:,:), allocatable :: P1
      REAL*8, dimension(:,:), allocatable :: R1
      REAL*8, dimension(:,:), allocatable :: R0hat
      REAL*8, dimension(:,:), allocatable :: UINIT
      REAL*8, dimension(:,:), allocatable :: RHS
      REAL*8, dimension(:,:), allocatable :: Y
      REAL*8, dimension(:,:), allocatable :: Hvec
      REAL*8, dimension(:,:), allocatable :: S
      REAL*8, dimension(:,:), allocatable :: T
      REAL*8, dimension(:,:), allocatable :: Z

      REAL*8 rho0,w0,rho1,alpha,w1,beta,a1,a2
      REAL*8 beta1,beta2
      integer do_shift,nsmooth

      if (precond_type.eq.1) then
       nsmooth=8
      else if (precond_type.eq.2) then
       nsmooth=2
      else
       print *,"precond_type invalid"
       stop
      endif

      if (probtype.eq.0) then
       do_shift=1
       beta1=(1.0,0.0)
       beta2=(-0.5,0.0)
      else
       do_shift=0
       beta1=(0.0,0.0)
       beta2=(0.0,0.0)
      endif

      allocate(U0(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(V0(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(P0(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(R0(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(U1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(V1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(P1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(R1(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(R0hat(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(UINIT(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(RHS(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(Y(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(Hvec(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(S(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(T(lox-1:hix+1,loy-1:hiy+1)) 
      allocate(Z(lox-1:hix+1,loy-1:hiy+1)) 

        ! U0=V0=P0=0 
      alpha=(0.0,0.0) 
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       U0(i,j)=alpha
       V0(i,j)=alpha
       P0(i,j)=alpha
      enddo
      enddo

       ! R0=G-A U0
      call RESID(nx,ny,lox,loy,hix,hiy, &
        R0,U0,G,DMINUS,KSQR,h,hflag,probtype, &
        gridtype,0,beta1,beta2) 
       ! R0hat=R0
      call COPYVEC(nx,ny,lox,loy,hix,hiy,R0,R0hat)
       ! UINIT=U0
      call COPYVEC(nx,ny,lox,loy,hix,hiy,U0,UINIT)
      call ZAPVEC(nx,ny,lox,loy,hix,hiy,U0)
       ! RHS=R0
      call COPYVEC(nx,ny,lox,loy,hix,hiy,R0,RHS)

       ! rho0=alpha=w0=1
      rho0=(1.0,0.0)
      alpha=(1.0,0.0)
      w0=(1.0,0.0)

      call NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,R0,dnorm0)
      dnorm=1.0
      maxiter=2000
      iter=0

      hflagcg=1
      do while ((dnorm.gt.tol).and.(iter.lt.maxiter))
       print *,"iter,dnorm ",iter,dnorm

         ! rho1= R0hat^H R0
       call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,R0hat,R0,rho1)
       
       restart_flag=0
       if ((sqrt(mag_complex(rho0)).lt.tol*0.01).or. &
           (sqrt(mag_complex(w0)).lt.tol*0.01)) then
        restart_flag=1
       endif 

       if (restart_flag.eq.0) then
          ! (R0hat^H R0)/(R0hat^H dot R0_before)  *   (alpha/w0)
        beta=(rho1/rho0)*(alpha/w0)
        a1=(1.0,0.0)
        a2=-w0
 
         ! P1=P0-w0 V0
        call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,P0,V0,P1,a1,a2)
         ! P1=R0+beta P1
        call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,R0,P1,P1,a1,beta)
        ! Y=M^{-1}P1
        call preconditioner( &
         nx,ny,lox,loy,hix,hiy,Y,P1, &
         DMINUS,KSQR,h,hflagcg,probtype, &
         gridtype,do_shift,beta1,beta2,nsmooth, &
         precond_type)
         
         ! V1=A Y
        call ATIMESU(nx,ny,lox,loy,hix,hiy, &
         V1,Y,DMINUS,KSQR,h,hflagcg,probtype, &
         gridtype,0,beta1,beta2)

         ! alpha=rho1/R0hat dot V1
        call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,R0hat,V1,alpha)

        if (sqrt(mag_complex(alpha)).lt.tol*0.01) then
         restart_flag=1
        endif

        if (restart_flag.eq.0) then
         alpha=rho1/alpha

         ! Hvec=U0+alpha Y
         a1=(1.0,0.0)
         a2=alpha
         call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,U0,Y,Hvec,a1,a2) 
         ! U1=Hvec
         call COPYVEC(nx,ny,lox,loy,hix,hiy,Hvec,U1)
         ! R1=RHS-A U1
         call RESID(nx,ny,lox,loy,hix,hiy, &
          R1,U1,RHS,DMINUS,KSQR,h,hflagcg,probtype, &
          gridtype,0,beta1,beta2) 
         call NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,R1,dnorm)
         dnorm=dnorm/dnorm0

         if (dnorm.gt.tol) then
          ! S=R0-alpha V1
          a1=(1.0,0.0)
          a2=-alpha
          call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,R0,V1,S,a1,a2) 

           ! Z=M^{-1}S
          call preconditioner( &
           nx,ny,lox,loy,hix,hiy,Z,S, &
           DMINUS,KSQR,h,hflagcg,probtype, &
           gridtype,do_shift,beta1,beta2,nsmooth, &
           precond_type)

           ! T=A Z
          call ATIMESU(nx,ny,lox,loy,hix,hiy, &
            T,Z,DMINUS,KSQR,h,hflagcg,probtype, &
            gridtype,0,beta1,beta2)

           ! simple case is: (T,S)/(T,T)=(AZ,S)/(AZ,AZ)   (MZ=S)
          call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,T,S,a1)
          call DOTPROD(nx,ny,lox,loy,hix,hiy,gridtype,T,T,a2)
          if (sqrt(mag_complex(a2)).lt.tol*0.01) then
           restart_flag=1
          endif

          if (restart_flag.eq.0) then
           w1=a1/a2
           ! U1=Hvec+w1 Z
           a1=(1.0,0.0)
           a2=w1
           call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,Hvec,Z,U1,a1,a2) 
          endif
         endif ! dnorm>tol
          ! R1=RHS-A U1
         call RESID(nx,ny,lox,loy,hix,hiy, &
          R1,U1,RHS,DMINUS,KSQR,h,hflagcg,probtype, &
          gridtype,0,beta1,beta2) 
         call NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,R1,dnorm)
         dnorm=dnorm/dnorm0
         rho0=rho1
         w0=w1
          ! R0=R1
         call COPYVEC(nx,ny,lox,loy,hix,hiy,R1,R0) 
         call COPYVEC(nx,ny,lox,loy,hix,hiy,P1,P0) 
         call COPYVEC(nx,ny,lox,loy,hix,hiy,V1,V0) 
         call COPYVEC(nx,ny,lox,loy,hix,hiy,U1,U0) 
        endif  ! restart_flag=0
       endif ! restart_flag=0

       if (restart_flag.eq.0) then
        ! do nothing
       else if (restart_flag.eq.1) then
        call RESID(nx,ny,lox,loy,hix,hiy, &
          R0,U0,RHS,DMINUS,KSQR,h,hflagcg,probtype, &
          gridtype,0,beta1,beta2) 
         ! R0hat=R0
        call COPYVEC(nx,ny,lox,loy,hix,hiy,R0,R0hat)
        call COPYVEC(nx,ny,lox,loy,hix,hiy,U0,U1) 
         ! rho0=alpha=w0=1
        rho0=(1.0,0.0)
        alpha=(1.0,0.0)
        w0=(1.0,0.0)
        call NORMPROD(nx,ny,lox,loy,hix,hiy,gridtype,R0,dnorm)
        dnorm=dnorm/dnorm0
        call ZAPVEC(nx,ny,lox,loy,hix,hiy,V0)
        call ZAPVEC(nx,ny,lox,loy,hix,hiy,P0)
       else
        print *,"restart_flag invalid"
        stop
       endif

       iter=iter+1
      enddo
      print *,"at the end: iter,dnorm ",iter,dnorm
       ! U=UINIT+U1
      a1=(1.0,0.0)
      a2=a1
      call LINCOMB(nx,ny,lox,loy,hix,hiy,gridtype,UINIT,U1,U,a1,a2)
 
      deallocate(U0) 
      deallocate(V0) 
      deallocate(P0) 
      deallocate(R0) 
      deallocate(U1) 
      deallocate(V1) 
      deallocate(P1) 
      deallocate(R1) 
      deallocate(R0hat) 
      deallocate(UINIT) 
      deallocate(RHS) 
      deallocate(Y) 
      deallocate(Hvec) 
      deallocate(S) 
      deallocate(T) 
      deallocate(Z) 

      return
      end



















