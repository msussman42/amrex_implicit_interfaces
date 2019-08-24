c multiply coefficients by 1/(1+eps)^(lam * alpha) 
      program main
      IMPLICIT NONE
      integer nterms,nplot,i,j
      real*8 eps,lam,alpha,L,dx,pi,x,original
      real*8 rapiddecay,slowdecay,base

      nterms=100
      nplot=1000
      eps=0.01
      alpha=0.001
      open(unit=12,file='decay')
      print *,"file is decay"
      L=1.0
      dx=L/nplot
      pi=4.0*atan(1.0)
      do i=0,nplot-1
       x=(i+0.5)*dx
       original=0.0
       rapiddecay=0.0
       slowdecay=0.0
       do j=0,nterms
        lam=(2.0*j+1)*pi/L
        base=-2.0*(L/((2.0*j+1.0)*pi))*sin(lam*x)
        rapiddecay=rapiddecay+base/exp(lam*lam*eps)
        original=original+base
        slowdecay=slowdecay+base/( (1.0+eps)**(lam*lam*alpha) )
       enddo
       write(12,*) x,original,rapiddecay,slowdecay
      enddo
      close(12)

      end

