!============================
! Spectral deferred correction Gauss Quadrature weights
!===========================

 subroutine SDC_GQweights(qbfact,tbfact,GQws)
 IMPLICIT NONE
 INTEGER_T qbfact,tbfact
 INTEGER_T i1,j1,k
 REAL_T GQws(0:tbfact,0:qbfact-1,1:tbfact) 
 REAL_T y(0:qbfact-1)
 REAL_T yGL(0:tbfact)
 REAL_T, dimension(:),allocatable :: bwGL,l
 REAL_T xpt

 do i1=0,qbfact-1
  y(i1)=cache_gauss(qbfact,i1)
 enddo
 do i1=0,tbfact
  yGL(i1)=cache_gauss_lobatto(tbfact,i1)
 enddo

 allocate(bwGL(0:tbfact),l(0:tbfact))
 
 call BarycentricWeights(tbfact,yGL,bwGL)
 
 do k=0,tbfact-1
   do i1=0,qbfact-1
    xpt=yGL(k)+(yGL(k+1)-yGL(k))*0.5D0*(y(i1)+1.0D0)
    call LagrangeInterpolatingPolynomial(tbfact,xpt,yGL,bwGL,l)
    do j1=0,tbfact
     GQws(j1,i1,k+1)=l(j1)
    enddo
   enddo !i1
 enddo !k
 
 deallocate(bwGL,l)

 return
 end subroutine SDC_GQweights
