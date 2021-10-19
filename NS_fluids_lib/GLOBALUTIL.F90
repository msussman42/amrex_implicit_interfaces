#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 0

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#include "AMReX_ArrayLim.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif



      module LagrangeInterpolationPolynomial
      IMPLICIT NONE

      CONTAINS
      SUBROUTINE AlmostEqual(a,b,var1)
      IMPLICIT NONE
      REAL_T, INTENT(IN) :: a,b
      LOGICAL, INTENT(OUT)      :: var1    !!!!!!!!!! logical: true and false
      REAL_T             :: tol

      tol = EPSILON(tol)                        !!!!!!!!!!! epsilon
      IF ((a == 0.0D0).OR.(b == 0.0D0)) THEN
      IF (ABS(a-b) <= 2.0D0*tol) THEN
         var1 = .TRUE.
      ELSE
         var1 = .FALSE.
      END IF
      ELSE
         IF ((ABS(a-b) <= ABS(a)*tol).and.(ABS(a-b) <= ABS(b)*tol)) THEN
            var1 = .TRUE.
         ELSE
            var1 = .FALSE.
         END IF
      END IF
      RETURN
      END SUBROUTINE AlmostEqual 


      SUBROUTINE MxVDerivative(N,D,Phi,Phi_prime)
      IMPLICIT NONE
      INTEGER, INTENT(IN)                           :: N
      REAL_T, DIMENSION(0:N), INTENT(IN)     :: Phi       !values at the nodes
      REAL_T, DIMENSION(0:N,0:N), INTENT(IN) :: D         !polynomial derivative matrix that has been pre-computed
      REAL_T, DIMENSION(0:N), INTENT(OUT)    :: Phi_prime !derivative of the interpolent at the nodes      
      INTEGER                                       :: i,j
      REAL_T                                 :: t

      DO i = 0,N
         t = 0.0D0
         DO j = 0,N
            t = t + D(i,j)*Phi(j)
         END DO
         Phi_prime(i) = t
      END DO
      RETURN
      END SUBROUTINE MxVDerivative


      SUBROUTINE BarycentricWeights(N,x,w)
      IMPLICIT NONE
      INTEGER, INTENT(IN)                        :: N
      REAL_T, DIMENSION(0:N), INTENT(IN)  :: x
      REAL_T, DIMENSION(0:N), INTENT(OUT) :: w
      INTEGER                                    :: j,k
   
      w = 1.0D0
      DO j = 1,N
         DO k = 0,j-1
            w(k) = w(k)*(x(k)-x(j))
            w(j) = w(j)*(x(j)-x(k))
         END DO
      END DO
      w = 1.0D0/w
      RETURN
      END SUBROUTINE BarycentricWeights


      SUBROUTINE LagrangeInterpolation(M,x,w,y,f,p) ! w is barycentric weights 
      IMPLICIT NONE
      INTEGER, INTENT(IN)                        :: M
      REAL_T, DIMENSION(0:M), INTENT(IN)  :: x,w,f
      REAL_T, INTENT(IN)                  :: y
      REAL_T, INTENT(OUT)                 :: p

      INTEGER                                    :: j
      REAL_T                              :: numerator,denominator,t
      LOGICAL                                    :: k

      numerator=0.0D0
      denominator=0.0D0
      DO j=0,M
         CALL AlmostEqual(y,x(j),k)
         IF (k) THEN
            p=f(j)
            RETURN
         ELSE
            t=w(j)/(y-x(j))
            numerator=numerator+t*f(j)
            denominator=denominator+t
         END IF
      END DO
      p=numerator/denominator
      RETURN
      END SUBROUTINE LagrangeInterpolation


!lagrange interpolation l(x)
! Algorithm 34
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!w is barycentric weights
      SUBROUTINE LagrangeInterpolatingPolynomial(N,xpt,x,w,l) 
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: N
      REAL_T, INTENT(IN) :: xpt
      REAL_T, DIMENSION(0:N), INTENT(IN)  :: x,w
      REAL_T, DIMENSION(0:N), INTENT(OUT) :: l

      INTEGER :: j,jmatch
      REAL_T  :: s,t
      LOGICAL :: var1, xMatchesNode

      xMatchesNode = .FALSE.
      jmatch=-1

      DO j = 0,N
         CALL AlmostEqual(xpt,x(j),var1)
         IF (var1) THEN
            l(j)         = one
            xMatchesNode = .TRUE.
            jmatch=j
         END IF
      END DO
      
      IF (.NOT.xMatchesNode) THEN
         s = zero
         DO j = 0,N
            t    = w(j)/(xpt-x(j))
            l(j) = t
            s    = s + t
         END DO
         DO j = 0,N
            l(j) = l(j)/s
         END DO
      else
         if ((jmatch.lt.0).or.(jmatch.gt.N)) then
          print *,"jmatch invalid"
          stop
         endif
         DO j=0,N
          if (j.ne.jmatch) then
           l(j)=zero
          endif
         enddo ! j
      END IF
      RETURN
      END SUBROUTINE LagrangeInterpolatingPolynomial


! Algorithm 36
!lagrange interpolation derivative;
!direct computation of the polynomial derivative in barycentric form
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!!!!--- w is barycentric weights
      SUBROUTINE LagInterpolantDerivative(M,x,y,f,w,p1) 
      IMPLICIT NONE
      INTEGER,INTENT(IN)                         :: M
      REAL_T, DIMENSION(0:M), INTENT(IN)  :: x,w,f
      REAL_T, INTENT(IN)                  :: y      
      REAL_T, INTENT(OUT)                 :: p1
      INTEGER                                    :: i,j
      REAL_T                              :: numerator,denominator,t,p
      LOGICAL                                    :: var1, atNode

      atNode = .FALSE.
      numerator=0.0D0
      
      DO j=0,M
         CALL AlmostEqual(y,x(j),var1)
         IF (var1) THEN
            atNode=.TRUE.
            p=f(j)
            denominator=-w(j)
            i=j
          END IF
      END DO

      IF (atnode) THEN
         DO j=0,M
            IF (j.NE.i) THEN
                numerator=numerator+w(j)*(p-f(j))/(y-x(j)) 
            END IF
         END DO
      ELSE
         denominator=0.0D0
         CALL LagrangeInterpolation(M,x,w,y,f,p) !!! ALGORITHM 31
         DO j=0,M
            t=w(j)/(y-x(j))
            numerator=numerator+t*(p-f(j))/(y-x(j))
            denominator=denominator+t
         END DO
      END IF
      
      p1=numerator/denominator
      
      RETURN
      END SUBROUTINE LagInterpolantDerivative 



!!!!!!!--- w is barycentric weights
      SUBROUTINE PolyDerivativeMatrix(N,x,w,dl)  
      IMPLICIT NONE
      INTEGER,INTENT(IN)                         :: N
      REAL_T, DIMENSION(0:N), INTENT(IN)  :: x,w
      REAL_T, DIMENSION(0:N,0:N), INTENT(OUT) :: DL
      INTEGER                                    :: i,j

      dl=0.0D0
      DO i=0,N
       ! dl(i,i)=0.0D0
         DO j=0,N
            IF (j.NE.i) THEN
               dl(i,j)=w(j)/(w(i)*(x(i)-x(j)))
               dl(i,i)=dl(i,i)-dl(i,j)
            END IF
         END DO
      END DO
      
      RETURN
      END SUBROUTINE  PolyDerivativeMatrix!LIDerivativeMatrixNodes


      end module LagrangeInterpolationPolynomial

       module LegendreNodes
       IMPLICIT NONE
   
!___________________________________________________________________


      REAL_T, dimension(:,:,:), allocatable :: cache_gauss 
      REAL_T, dimension(:,:,:), allocatable :: cache_gauss_lobatto 
      REAL_T, dimension(:,:,:), allocatable :: cache_gauss_w 
      REAL_T, dimension(:,:,:), allocatable :: cache_gauss_lobatto_w 

      REAL_T , dimension(:,:,:,:), allocatable :: cache_wMATGL
      REAL_T , dimension(:,:,:,:), allocatable :: cache_wMAT
      REAL_T , dimension(:,:,:,:), allocatable :: cache_w_right
      REAL_T , dimension(:,:,:,:), allocatable :: cache_w_left

      REAL_T , dimension(:,:,:,:), allocatable :: cache_wMAT_extend
      REAL_T , dimension(:,:,:,:), allocatable :: cache_wRT
      REAL_T , dimension(:,:,:,:), allocatable :: cache_wLT
      REAL_T , dimension(:,:,:,:), allocatable :: cache_wLRT

      REAL_T , dimension(:,:,:,:), allocatable :: cache_wRT_EXT
      REAL_T , dimension(:,:,:,:), allocatable :: cache_wLT_EXT

      INTEGER_T, PARAMETER :: SPTYPE = 0  ! Legendre
      INTEGER_T, PARAMETER :: TMTYPE = 0  ! Legendre
      INTEGER_T, PARAMETER :: GQTYPE = 0  ! Legendre

CONTAINS

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 1
! legendre polynomial & derivative
! Algorithm 22
      SUBROUTINE LegendrePolyAndDeri(k,x,P,dl)
     !k=the degree of Legendre polynomial
     !P=Legendre polynomial; 
     !dl=its derivative;
      IMPLICIT NONE 
      INTEGER,INTENT(IN)          :: k
      REAL_T, INTENT(IN)   :: x
      REAL_T, INTENT(OUT)  :: P,dl
      INTEGER              :: j
      REAL_T               :: P1,P2,dl1,dl2

      if (k.eq.0) then ! constant 
          P=one
          dl=zero
      else if (k.eq.1) then
              P=x
              dl=one
      else
             P1=one
             P2=x
             dl1=zero
             dl2=one
             do j=2,k
                P=(2*j-1)*x*P2/j-(j-1)*P1/j
                dl=(2*j-1)*P2+dl1
                P1=P2
                P2=P
                dl1=dl2
                dl2=dl
             end do
      end if
      RETURN
      END SUBROUTINE LegendrePolyAndDeri

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 2
! Algorithm 23
! Legendre Gauss
      SUBROUTINE LegendreGaussNodesAndWeights(n,y,weight)
      !y=Legendre Gauss Nodes
      !weight is Legendre Gauss weights
      !return y(0:n),weight(0:n);
      INTEGER,INTENT(IN)          :: n
      REAL_T, DIMENSION(0:n), INTENT(OUT)  :: y,weight
      INTEGER                     :: j,k
      REAL_T               :: P1,delta
      REAL_T               :: dP1,localPI

        localPI=four*atan(one)

        if (n.eq.0) then
            y(0)=zero
            weight(0)=two
        else if (n.eq.1) then
                y(0)=-sqrt(one/three)
                weight(0)=one
                y(1)=-y(0)
                weight(1)=weight(0)
        else
             do j=0,int((n+1)/two)-1
                 y(j)=-cos(((2*j+1)*localPI)/(2*n+2))
                 do k=0,4
                    call LegendrePolyAndDeri(n+1,y(j),P1,dP1)
                    delta=-(P1/dP1)
                    y(j)=y(j)+delta
                    if (abs(delta).le.epsilon(four)*abs(y(j))) then
                                 exit
                    end if
                 end do
                 call LegendrePolyAndDeri(n+1,y(j),P1,dP1)
                 y(n-j)=-y(j)
                 weight(j)=two/((1-y(j)**2)*(dP1**2))
                 weight(n-j)=weight(j)
             end do
        end if

        if (mod(n,2).eq.0) then
            call LegendrePolyAndDeri(n+1,zero,P1,dP1)
            y(n/2)=zero
            weight(n/2)=two/(dP1**2)
        end if
      RETURN
      END SUBROUTINE LegendreGaussNodesAndWeights


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 3
! legendre polynomial & derivative; Q=L_N+1-L_N-1
! Algorithm 24
      SUBROUTINE qAndLEvaluation(n,x,q,dq,ln) 
      !k=the degree of Legendre polynomial
      !P=Legendre polynomial; 
      !dl=its derivative;
      !n>=2
      INTEGER,INTENT(IN)          :: n
      REAL_T, INTENT(IN)   :: x
      REAL_T, INTENT(OUT)  :: q,dq,ln
      INTEGER                     :: k
      REAL_T               :: l1,l2,l3,dl1,dl2,dl,dl3
      !l1=l_{n-2}
      !l2=l_{n-1}
      !l3=l_{n+1}

         k=2
         l1=one
         l2=x
         dl1=zero
         dl2=one

         do k=2,n
            ln=((2*k-1)*x*l2-(k-1)*l1)/k
            dl=dl1+(2*k-1)*l2
            l1=l2
            l2=ln
            dl1=dl2
            dl2=dl
         end do

         k=n+1
         l3=((2*k-1)*x*ln-(k-1)*l1)/k
         dl3=dl1+(2*k-1)*l2
         q=l3-l1
         dq=dl3-dl1
      RETURN
      END SUBROUTINE qAndLEvaluation


       ! y1 and y2 are GL points
      subroutine intersect_table(n1,n2,y1,y2,w)
      IMPLICIT NONE

      INTEGER_T n1,n2,i,j
      REAL_T y1(0:n1),y2(0:n2)
      REAL_T w(0:n1-1,0:n2-1),wt,low,high

      do i=0,n1-1
       do j=0,n2-1
        if (y1(i).ge.y2(j+1)) then
         wt=zero
        else if (y1(i+1).le.y2(j)) then
         wt=zero
        else
         low=max(y1(i),y2(j))
         high=min(y1(i+1),y2(j+1))
         if (high-low.gt.zero) then
          wt=high-low
         else
          wt=zero
         endif
        endif
        w(i,j)=wt
       enddo
      enddo

      return
      end subroutine intersect_table

       ! y1 and y2 are G points
      subroutine intersect_tableGL(n1,n2,y1,y2,w)
      IMPLICIT NONE

      INTEGER_T n1,n2,i,j
      REAL_T y1(0:n1),y2(0:n2)
      REAL_T y1ext(0:n1+2),y2ext(0:n2+2)
      REAL_T w(0:n1+1,0:n2+1),wt,low,high

      do i=0,n1
       y1ext(i+1)=y1(i)
      enddo
      y1ext(0)=-one
      y1ext(n1+2)=one

      do i=0,n2
       y2ext(i+1)=y2(i)
      enddo
      y2ext(0)=-one
      y2ext(n2+2)=one

      do i=0,n1+1
       do j=0,n2+1
        if (y1ext(i).ge.y2ext(j+1)) then
         wt=zero
        else if (y1ext(i+1).le.y2ext(j)) then
         wt=zero
        else
         low=max(y1ext(i),y2ext(j))
         high=min(y1ext(i+1),y2ext(j+1))
         if (high-low.gt.zero) then
          wt=high-low
         else
          wt=zero
         endif
        endif
        w(i,j)=wt
       enddo
      enddo

      return
      end subroutine intersect_tableGL


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 4
! Algorithm 25
! Legendre Gauss-Lobatto
      SUBROUTINE LegendreGaussLobattoNodesAndWeights(n,y,weight)
      !y=Legendre Gauss Lobatto Nodes
      !weight is Legendre Gauss Lobatto weights
      !return y(0:n),weight(0:n);
      INTEGER,INTENT(IN)          :: n
      REAL_T, DIMENSION(0:n), INTENT(OUT)  :: y,weight
      INTEGER                     :: j,k
      REAL_T               :: q,dq,ln,delta,localPI

        localPI=four*atan(one)
        if (n.eq.1) then
           y(0)=-one
           weight(0)=one
           y(1)=one
           weight(1)=weight(0)
        else if (n.gt.1) then
           y(0)=-one
           weight(0)=two/(n*(n+1))
           y(n)=one
           weight(n)=weight(0)
           do j=1,int((n+1)/two)-1
              y(j)=-cos((j+fourth)*localPI/n-three/(8*n*localPI*(j+fourth)))
              do k=0,4
                 call qAndLEvaluation(n,y(j),q,dq,ln) 
                 delta=-(q/dq)
                 y(j)=y(j)+delta
                 if (abs(delta).le.epsilon(four)*abs(y(j))) then
                             exit
                 end if
              end do
              call qAndLEvaluation(n,y(j),q,dq,ln) 
              y(n-j)=-y(j)
              weight(j)=two/(n*(n+1)*ln**2)
              weight(n-j)=weight(j)
           end do
        else
         print *,"n invalid"
         stop
        end if

        if (mod(n,2).eq.zero) then
           call qAndLEvaluation(n,zero,q,dq,ln) 
           k=int(n/two)
           y(k)=zero
           weight(k)=two/(n*(n+1)*ln**2)
        end if
      RETURN
      END SUBROUTINE LegendreGaussLobattoNodesAndWeights 
    
      subroutine do_polyinterp(order_r,w,data,sum)
      IMPLICIT NONE

      INTEGER_T i,order_r
      REAL_T w(0:order_r)
      REAL_T data(0:order_r)
      REAL_T sum

      sum=zero
      do i=0,order_r
       sum=sum+data(i)*w(i)
      enddo

      return
      end subroutine do_polyinterp
       

        ! r is the order of the polynomial interpolant
      subroutine polyinterp_weights(order_r,x,w,xtarget)
      IMPLICIT NONE

      INTEGER_T order_r,i,j
      REAL_T xtarget,lag
      REAL_T x(0:order_r)
      REAL_T w(0:order_r)

      do i=0,order_r
       lag=one
       do j=0,order_r
        if (j.ne.i) then
         if (abs(x(i)-x(j)).eq.zero) then
          print *,"x corruption"
          stop
         endif
         lag=lag*(xtarget-x(j))/(x(i)-x(j))
        endif
       enddo
       w(i)=lag
      enddo

      return
      end subroutine polyinterp_weights


      subroutine ClenshawGaussLobattoNodesAndWeights(n,y,weight)
      INTEGER,INTENT(IN)          :: n
      REAL_T, DIMENSION(0:n), INTENT(OUT)  :: y,weight
      INTEGER                     :: j,k
      REAL_T               :: localPI
      REAL_T w(0:n)
      REAL_T yL(0:n)
      REAL_T wL(0:n)

      localPI=four*atan(one)

      if (n.gt.33) then
       print *,"n too big"
       stop
      endif

      do j=0,n
       y(j)=-cos(j*localPI/n)
       weight(j)=zero
      enddo

      call LegendreGaussNodesAndWeights(n,yL,wL)
      do k=0,n
       call polyinterp_weights(n,y,w,yL(k))
       do j=0,n
        weight(j)=weight(j)+wL(k)*w(j)
       enddo
      enddo

      return
      end subroutine ClenshawGaussLobattoNodesAndWeights


      subroutine ClenshawGaussNodesAndWeights(n,y,weight)
      INTEGER,INTENT(IN)          :: n
      REAL_T, DIMENSION(0:n), INTENT(OUT)  :: y,weight
      INTEGER                     :: j,k
      REAL_T               :: localPI
      REAL_T w(0:n)
      REAL_T yL(0:n)
      REAL_T wL(0:n)

      localPI=four*atan(one)

      if (n.gt.32) then
       print *,"n too big"
       stop
      endif

      do j=0,n
       y(j)=-cos((two*j+1)*localPI/(2*n+2))
       weight(j)=zero
      enddo

      call LegendreGaussNodesAndWeights(n,yL,wL)
      do k=0,n
       call polyinterp_weights(n,y,w,yL(k))
       do j=0,n
        weight(j)=weight(j)+wL(k)*w(j)
       enddo
      enddo

      return
      end subroutine ClenshawGaussNodesAndWeights



        ! (order_r+1)/2 is the order of each polynomial interpolant
      subroutine polyinterp_weights_split(order_r,x,w,xtarget)
      IMPLICIT NONE

      INTEGER_T order_r,i,j,bfact,bfactC
      REAL_T xtarget,lag
      REAL_T x(0:order_r)
      REAL_T w(0:order_r)

      bfact=order_r+1
      bfactC=bfact/2
      if (bfactC*2.ne.bfact) then
       print *,"bfact invalid12"
       stop
      endif
      
      do i=0,bfactC-1
       lag=one
       do j=0,bfactC-1
        if (j.ne.i) then
         if (abs(x(i)-x(j)).eq.zero) then
          print *,"x corruption"
          stop
         endif
         lag=lag*(xtarget-x(j))/(x(i)-x(j))
        endif
       enddo
       if (xtarget.gt.zero) then
        w(i)=zero
        w(i+bfactC)=lag
       else
        w(i)=lag
        w(i+bfactC)=zero
       endif
      enddo

      return
      end subroutine polyinterp_weights_split



      subroutine polyinterp_weights_even(order_r,x,w,xtarget)
      IMPLICIT NONE

      INTEGER_T order_r,i,bfact
      REAL_T xtarget
      REAL_T x(0:order_r)
      REAL_T w(0:order_r)
      REAL_T dx,totalw


      bfact=order_r+1
      dx=two/bfact

      if (bfact.eq.1) then
       w(0)=one
      else if (bfact.gt.1) then
       totalw=zero
       do i=0,order_r
        if ((xtarget.ge.x(i)-half*dx).and. &
            (xtarget.le.x(i)+half*dx)) then
         w(i)=one
        else
         w(i)=zero
        endif
        totalw=totalw+w(i)
       enddo
       if (totalw.lt.one) then
        print *,"totalw invalid"
        stop
       else
        do i=0,order_r
         w(i)=w(i)/totalw
        enddo
       endif
      else
       print *,"bfact invalid13"
       stop
      endif
     
      return
      end subroutine polyinterp_weights_even


       ! w_ij=L_i' (x_j)
       ! L_i=PI_m<>i (x-xm)/(xi-xm)
       ! L_i'=sum_l<>i (1/(xi-xl)) PI_m<>l,i (x-xm)/(xi-xm)
       ! L_i'(x_j)=sum_l<>i (1/(xi-xl)) PI_m<>l,i (xj-xm)/(xi-xm)
       ! L_i'(x_i)=sum_l<>i (1/(xi-xl)) 
       ! if j<>i,
       ! L_i'(x_j)=sum_l<>i (1/(xi-xl)) PI_m<>l,i (xj-xm)/(xi-xm)
       !          =(1/(xi-xj))PI_m<>i,j (xj-xm)/(xi-xm)
       ! if f(x)=sum_i f_i L_i(x)  then
       ! f'(x)=sum_i f_i L_i'(x)
       ! f'(x_j)=sum_i f_i L_i'(x_j)
       ! f'(x)=sum_i f'(x_i)L_i(x)=
       !       sum_i (sum_k f_k L_k'(x_i))L_i(x)
       ! suppose f_i=constant then
       ! sum_i L_i'(x_j)=0
       ! wij'=wij - lambda
       ! 0=sum_i wij - (order_r+1)lambda
      subroutine polyinterp_Dmatrix(order_r,x,w)
      IMPLICIT NONE

      INTEGER_T order_r,i,j,l,m
      REAL_T x(0:order_r)
      REAL_T w(0:order_r,0:order_r)
      REAL_T sum

      if (order_r.eq.0) then
       w(0,0)=zero
      else if (order_r.eq.1) then
       w(0,0)=one/(x(0)-x(1))
       w(0,1)=w(0,0)
       w(1,0)=-w(0,0)
       w(1,1)=w(1,0)
      else if (order_r.gt.1) then
       do i=0,order_r
       do j=0,order_r
        if (i.eq.j) then
         sum=zero
         do l=0,order_r
          if (l.ne.i) then
           sum=sum+one/(x(i)-x(l))
          endif
         enddo
         w(i,j)=sum
        else if (i.ne.j) then 
         sum=one
         do m=0,order_r
          if ((m.ne.i).and.(m.ne.j)) then
           sum=sum*(x(j)-x(m))/(x(i)-x(m))
          endif
         enddo
         sum=sum/(x(i)-x(j))
        else
         print *,"i,j bust"
         stop
        endif
        w(i,j)=sum
       enddo
       enddo
       do j=0,order_r
        sum=zero
        do i=0,order_r
         sum=sum+w(i,j) 
        enddo
        sum=sum/(order_r+one)
        do i=0,order_r
         w(i,j)=w(i,j)-sum
        enddo 
       enddo
      else
       print *,"order_r invalid"
       stop
      endif

      return
      end subroutine polyinterp_Dmatrix

      subroutine poly_change_basis(order_r1,order_r2,data1,data2,y1,y2)
      IMPLICIT NONE

      INTEGER_T order_r1,order_r2,i
      REAL_T y1(0:order_r1)
      REAL_T y2(0:order_r2)
      REAL_T data1(0:order_r1)
      REAL_T w(0:order_r1)
      REAL_T data2(0:order_r2)

      do i=0,order_r2
       call polyinterp_weights(order_r1,y1,w,y2(i))
       call do_polyinterp(order_r1,w,data1,data2(i))
      enddo

      return
      end subroutine poly_change_basis

      subroutine deriv_change_basis(order_r1,order_r2,data,datader, &
        wMAT,y1,y2,dx_element)
      IMPLICIT NONE

      INTEGER_T order_r1,order_r2
      REAL_T y1(0:order_r1)
      REAL_T y2(0:order_r2)
      REAL_T data(0:order_r1)
      REAL_T datader1(0:order_r1)
      REAL_T datader(0:order_r2)
      REAL_T wMAT(0:order_r1,0:order_r1)
      INTEGER_T i,j
      REAL_T sum,dx_element
   
      do i=0,order_r1
       sum=zero
       do j=0,order_r1
        sum=sum+data(j)*wMAT(j,i)
       enddo
       datader1(i)=sum
      enddo
      call poly_change_basis(order_r1,order_r2,datader1,datader,y1,y2)
      do i=0,order_r2
       datader(i)=datader(i)*two/dx_element
      enddo

      return
      end subroutine deriv_change_basis

      subroutine delete_cache()
      IMPLICIT NONE
 
      deallocate(cache_gauss)
      deallocate(cache_gauss_lobatto)
      deallocate(cache_gauss_w)
      deallocate(cache_gauss_lobatto_w)

      deallocate(cache_wMATGL)
      deallocate(cache_wMAT)
      deallocate(cache_w_right)
      deallocate(cache_w_left)

      deallocate(cache_wMAT_extend)
      deallocate(cache_wRT)
      deallocate(cache_wLT)
      deallocate(cache_wLRT)

      deallocate(cache_wRT_EXT)
      deallocate(cache_wLT_EXT)

      return
      end subroutine delete_cache

       ! typ=0  Legendre
       ! typ=1  Clenshaw Curtis 
       ! r is the largest "order" aka number of points that
       ! one can prescribe.
       ! 
      subroutine init_cache(order_r)
      IMPLICIT NONE

      INTEGER_T order_r
      REAL_T yleft(0:order_r)
      REAL_T yright(0:order_r)
      REAL_T y(0:order_r)
      REAL_T yGL(0:order_r)
      REAL_T y_extend(0:order_r+1)
      REAL_T yGL_extend(0:order_r+1)
      REAL_T yLT(0:order_r+1)
      REAL_T yRT(0:order_r+1)
      REAL_T yLRT(0:order_r+1)

      REAL_T yLRTextrap(0:order_r-1)
      REAL_T yLTextrap(0:order_r)
      REAL_T yRTextrap(0:order_r)

      REAL_T yLTextrapEXT(0:order_r)
      REAL_T yRTextrapEXT(0:order_r)

      REAL_T, dimension(:,:), allocatable :: deriv_matrix
      REAL_T, dimension(:), allocatable :: tempx,tempw

      INTEGER_T i,j,typ,i1,j1

      if (order_r.lt.1) then
       print *,"order_r invalid"
       stop
      endif

      allocate(cache_gauss(1:order_r+1,0:order_r+1,0:1))
      allocate(cache_gauss_lobatto(1:order_r+1,0:order_r+1,0:1))
      allocate(cache_gauss_w(1:order_r+1,0:order_r+1,0:1))
      allocate(cache_gauss_lobatto_w(1:order_r+1,0:order_r+1,0:1))

      allocate(cache_wMATGL(1:order_r,0:order_r,0:order_r,0:1))
      allocate(cache_wMAT(1:order_r,0:order_r,0:order_r,0:1))
      allocate(cache_w_left(1:order_r,0:order_r,0:order_r,0:1))
      allocate(cache_w_right(1:order_r,0:order_r,0:order_r,0:1))

      allocate(cache_wMAT_extend(1:order_r,0:order_r+1,0:order_r+1,0:1))
      allocate(cache_wRT(1:order_r,0:order_r+1,0:order_r+1,0:1))
      allocate(cache_wLT(1:order_r,0:order_r+1,0:order_r+1,0:1))
      allocate(cache_wLRT(1:order_r,0:order_r+1,0:order_r+1,0:1))

      allocate(cache_wRT_EXT(1:order_r,0:order_r+1,0:order_r+1,0:1))
      allocate(cache_wLT_EXT(1:order_r,0:order_r+1,0:order_r+1,0:1))

      do typ=0,1

       do i=1,order_r+1

        allocate(tempx(0:i))
        allocate(tempw(0:i))

         ! 0..i
        if (typ.eq.0) then
         call LegendreGaussLobattoNodesAndWeights(i,tempx,tempw) 
        else if (typ.eq.1) then
         call ClenshawGaussLobattoNodesAndWeights(i,tempx,tempw) 
        else
         print *,"typ invalid"
         stop
        endif

        do j=0,i
         cache_gauss_lobatto(i,j,typ)=tempx(j)
         cache_gauss_lobatto_w(i,j,typ)=tempw(j)
        enddo

         ! 0..i-1
        if (typ.eq.0) then
         call LegendreGaussNodesAndWeights(i-1,tempx,tempw) 
        else if (typ.eq.1) then
         call ClenshawGaussNodesAndWeights(i-1,tempx,tempw) 
        else
         print *,"typ invalid"
         stop
        endif

        do j=0,i-1
         cache_gauss(i,j,typ)=tempx(j)
         cache_gauss_w(i,j,typ)=tempw(j)
        enddo

        deallocate(tempx,tempw)

       enddo ! i=1.. order_r+1

       do i=1,order_r

        do i1=0,i-1
         y(i1)=cache_gauss(i,i1,typ)
        enddo
        do i1=0,i
         yGL(i1)=cache_gauss_lobatto(i,i1,typ)
        enddo
        yright(0)=-one
        yleft(i)=one
        do i1=1,i
         yright(i1)=y(i1-1)
         yleft(i1-1)=y(i1-1)
        enddo
 
        y_extend(0)=-one-abs(y(0)+one)
        do i1=1,i
         y_extend(i1)=y(i1-1)
        enddo
        y_extend(i+1)=one+abs(one-y(i-1))
        do i1=0,i+1
         yGL_extend(i1)=cache_gauss_lobatto(i+1,i1,typ)
        enddo
        yLT(0)=-one
        yLT(i+1)=y_extend(i+1)
        yRT(0)=y_extend(0)
        yRT(i+1)=one
        yLRT(0)=-one
        yLRT(i+1)=one
        do i1=1,i
         yLT(i1)=y(i1-1)
         yRT(i1)=y(i1-1)
         yLRT(i1)=y(i1-1)
        enddo ! i1

        yLTextrap(i)=y_extend(i+1)
        yLTextrapEXT(i)=one
        do i1=0,i-1
         yLTextrap(i1)=y(i1)
         yLTextrapEXT(i1)=y(i1)
        enddo

        yRTextrap(0)=y_extend(0)
        yRTextrapEXT(0)=-one
        do i1=0,i-1
         yRTextrap(i1+1)=y(i1)
         yRTextrapEXT(i1+1)=y(i1)
        enddo

        do i1=0,i-1
         yLRTextrap(i1)=y(i1)
        enddo

         ! N=Neumann D=Dirichlet E=extrap I=Interior
        allocate(deriv_matrix(0:i,0:i))
        call polyinterp_Dmatrix(i,yRTextrapEXT,deriv_matrix)
        do i1=0,i
        do j1=0,i
         cache_wLT_EXT(i,i1,j1,typ)=deriv_matrix(i1,j1) !N=left E=right
        enddo
        enddo
        deallocate(deriv_matrix)

        allocate(deriv_matrix(0:i,0:i))
        call polyinterp_Dmatrix(i,yLTextrapEXT,deriv_matrix)
        do i1=0,i
        do j1=0,i
         cache_wRT_EXT(i,i1,j1,typ)=deriv_matrix(i1,j1) !E=left N=right
        enddo
        enddo
        deallocate(deriv_matrix)

        allocate(deriv_matrix(0:i+1,0:i+1))
        call polyinterp_Dmatrix(i+1,yGL_extend,deriv_matrix)
        do i1=0,i+1
        do j1=0,i+1
         cache_wMAT_extend(i,i1,j1,typ)=deriv_matrix(i1,j1)
        enddo
        enddo
        deallocate(deriv_matrix)


        allocate(deriv_matrix(0:i+1,0:i+1))
        call polyinterp_Dmatrix(i+1,yRT,deriv_matrix)
        do i1=0,i+1
        do j1=0,i+1
         cache_wRT(i,i1,j1,typ)=deriv_matrix(i1,j1)
        enddo
        enddo
        deallocate(deriv_matrix)

        allocate(deriv_matrix(0:i+1,0:i+1))
        call polyinterp_Dmatrix(i+1,yLT,deriv_matrix)
        do i1=0,i+1
        do j1=0,i+1
         cache_wLT(i,i1,j1,typ)=deriv_matrix(i1,j1)
        enddo
        enddo
        deallocate(deriv_matrix)

        allocate(deriv_matrix(0:i+1,0:i+1))
        call polyinterp_Dmatrix(i+1,yLRT,deriv_matrix)
        do i1=0,i+1
        do j1=0,i+1
         cache_wLRT(i,i1,j1,typ)=deriv_matrix(i1,j1)
        enddo
        enddo
        deallocate(deriv_matrix)


        allocate(deriv_matrix(0:i,0:i))
        call polyinterp_Dmatrix(i,yGL,deriv_matrix)
        do i1=0,i
        do j1=0,i
         cache_wMATGL(i,i1,j1,typ)=deriv_matrix(i1,j1)
        enddo
        enddo
        deallocate(deriv_matrix)

        allocate(deriv_matrix(0:i-1,0:i-1))
        call polyinterp_Dmatrix(i-1,y,deriv_matrix)
        do i1=0,i-1
        do j1=0,i-1
         cache_wMAT(i,i1,j1,typ)=deriv_matrix(i1,j1)
        enddo
        enddo
        deallocate(deriv_matrix)
 
        allocate(deriv_matrix(0:i,0:i))
        call polyinterp_Dmatrix(i,yright,deriv_matrix)
        do i1=0,i
        do j1=0,i
         cache_w_right(i,i1,j1,typ)=deriv_matrix(i1,j1)
        enddo
        enddo
        deallocate(deriv_matrix)

        allocate(deriv_matrix(0:i,0:i))
        call polyinterp_Dmatrix(i,yleft,deriv_matrix)
        do i1=0,i
        do j1=0,i
         cache_w_left(i,i1,j1,typ)=deriv_matrix(i1,j1)
        enddo
        enddo
        deallocate(deriv_matrix)

       enddo ! i

      enddo ! typ

      return
      end subroutine init_cache
 
      subroutine sanity_check(rend)
      IMPLICIT NONE
 
      INTEGER_T typ,order_r,i,p,rstart,rend
      REAL_T, dimension(:), allocatable :: y,weight,data
      REAL_T, dimension(:), allocatable :: w,datader
      REAL_T, dimension(:,:), allocatable :: wMAT
      REAL_T exact,sum,xtarget,dx_element

        ! sanity checks will fail at r=20
      do typ=0,3
       if ((typ.eq.0).or.(typ.eq.2)) then
        rstart=0
       else
        rstart=1
       endif
       do order_r=rstart,rend
        allocate(y(0:order_r))
        allocate(weight(0:order_r))
        allocate(w(0:order_r))
        allocate(wMAT(0:order_r,0:order_r))
        allocate(data(0:order_r))
        allocate(datader(0:order_r))
        if (typ.eq.0) then
         call LegendreGaussNodesAndWeights(order_r,y,weight)
        else if (typ.eq.1) then
         call LegendreGaussLobattoNodesAndWeights(order_r,y,weight)
        else if (typ.eq.2) then
         call ClenshawGaussNodesAndWeights(order_r,y,weight)
        else if (typ.eq.3) then
         call ClenshawGaussLobattoNodesAndWeights(order_r,y,weight)
        endif
        do i=0,order_r-1
         if (y(i).ge.y(i+1)) then
          print *,"nodes out of order"
          stop
         endif
        enddo
        do p=0,order_r
         do i=0,order_r
          data(i)=y(i)**p
         enddo
         sum=zero
         do i=0,order_r
          sum=sum+weight(i)*data(i)
         enddo
         if ((p/2)*2.eq.p) then
          exact=two/(p+one)
         else
          exact=zero
         endif
         if (abs(sum-exact).gt.1.0D-13) then
          print *,"sanity check failed integral"
          print *,"typ,order_r,p,exact,sum ",typ,order_r,p,exact,sum
          stop
         endif
         do i=0,order_r
          if (order_r.eq.0) then
           xtarget=zero
          else if (order_r.gt.0) then
           xtarget=-one+two*i/order_r
          else
           print *,"order_r invalid"
           stop
          endif
          call polyinterp_weights(order_r,y,w,xtarget)
          call do_polyinterp(order_r,w,data,sum)
          exact=xtarget**p
          if (abs(sum-exact).gt.1.0D-13) then
           print *,"sanity check failed polyinterp"
           print *,"r,p ",order_r,p
           stop
          endif
          call polyinterp_Dmatrix(order_r,y,wMAT)
          dx_element=two
          call deriv_change_basis(order_r,order_r,data,datader, &
           wMAT,y,y,dx_element)

          call polyinterp_weights(order_r,y,w,xtarget)
          call do_polyinterp(order_r,w,datader,sum)
          if (p.eq.0) then
           exact=zero
          else
           exact=p*(xtarget**(p-1))
          endif
          if (abs(sum-exact).gt.1.0D-10) then
           print *,"sanity check failed polyinterp_Dmatrix"
           print *,"r,p,typ ",order_r,p,typ
           print *,"sum-exact= ",sum-exact
           stop
          endif
         enddo ! i
        enddo ! p
        deallocate(y)
        deallocate(weight)
        deallocate(w)
        deallocate(wMAT)
        deallocate(data)
        deallocate(datader)
       enddo ! r
      enddo ! typ

      end subroutine sanity_check 
 
      end module LegendreNodes


module global_utility_module

implicit none

      type nucleation_parm_type_inout
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LSnew
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: Snew
      end type nucleation_parm_type_inout

REAL_T :: MOF_PI=zero

contains


! nfree points towards the liquid
! nsolid points away from the solid 
! (nsolid is gradient of psi where psi<0 in solid)
! angle=pi => dry  contact=0 => wetting
! nghost points towards the liquid
! nperp points into the solid tangent to the contours of the extrapolated
!  level set functions.
      subroutine ghostnormal(nfree,nsolid,cos_angle,nghost,nperp)
      IMPLICIT NONE

      REAL_T nfree(SDIM)
      REAL_T nsolid(SDIM)
      REAL_T nghost(SDIM)
      REAL_T nperp(SDIM)
      REAL_T e2(3),e3(3),ntemp(3)
      INTEGER_T i
      REAL_T ss,cos_angle,sin_angle,dist
      REAL_T nsolid_new(SDIM)
      INTEGER_T dir
      REAL_T mag


        ! points away from solid
      do dir=1,SDIM
        nsolid_new(dir)=nsolid(dir)
      enddo

      mag=zero
      do dir=1,SDIM
       mag=mag+nsolid_new(dir)**2
      enddo
      mag=sqrt(mag)
      if (mag.le.zero) then
       print *,"solid normal mag=0 error"
       stop
      else
       do dir=1,SDIM
        nsolid_new(dir)=nsolid_new(dir)/mag
       enddo
      endif

      if ((cos_angle.lt.-one).or.(cos_angle.gt.one)) then
       print *,"cosine out of range"
       stop
      endif
      sin_angle=sqrt(one-cos_angle**2)

! e2 is tangent to the contact line.
! in 2d, e2 comes into/out of the paper.
#if (AMREX_SPACEDIM==3)
      call crossprod(nsolid_new,nfree,e2)
#elif (AMREX_SPACEDIM==2)
      call crossprod2d(nsolid_new,nfree,e2)
#else
      print *,"dimension bust"
      stop
#endif

      dist=sqrt( e2(1)**2+e2(2)**2+e2(3)**2+1.0D-14 )
      do i=1,3
       e2(i)=e2(i)/dist
      enddo

! nsolid_new is the normalized version of nsolid.
! e2 is tangent to the contact line.
! e3 is parallel to the solid interface in a plane that is orthogonal to the 
! contact line.

      do i=1,SDIM
       ntemp(i)=nsolid_new(i)
      enddo
      if (SDIM.eq.2) then
       ntemp(3)=zero
      else if (SDIM.eq.3) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif
      call crossprod(ntemp,e2,e3)

      dist=sqrt( e3(1)**2+e3(2)**2+e3(3)**2+1.0D-14 )
      do i=1,3
       e3(i)=e3(i)/dist
      enddo
      if (SDIM.eq.3) then
       ss=nfree(1)*e3(1)+nfree(2)*e3(2)+nfree(SDIM)*e3(SDIM)
      else if (SDIM.eq.2) then
       ss=nfree(1)*e3(1)+nfree(2)*e3(2)
      else
       print *,"dimension bust"
       stop
      endif
      if (ss.ge.zero) then
       ss=one
      else
       ss=-one
      endif

! in 2d with phi>0 in the droplet and phi<0 outside the droplet,
! and on a solid substrate,
! for the right side of the droplet, 
! nfree=(-nx,ny,0) where nx is positive and ny is pos or neg.
! nsolid=(0 1 0)
! e2=nsolid x nfree=| i   j k |=(0 0 nx)  (nx is positive on right side)
!                   | 0   1 0 |
!                   |-nx ny 0 |
! e2 normalized=(0 0 1)  (right side) 
! e3=nsolid x e2 = | i j k | = (1 0 0)
!                  | 0 1 0 |
!                  | 0 0 1 | 
! ss=nfree dot e3 = (normalized) -1
! nghost=(-sin(theta),-cos(theta)) 
! nperp=(cos(theta),-sin(theta))
!
! for the left side of the droplet,
! nfree=(nx,ny,0) where nx is positive and ny is pos or neg.
! nsolid=(0 1 0)
! e2=nsolid x nfree=| i   j k |=(0 0 -nx)  (nx is positive on left side)
!                   | 0   1 0 |
!                   | nx ny 0 |
! e2 normalized=(0 0 -1)  (left side)
! e3=nsolid x e2 = | i j k  | = (-1 0 0)
!                  | 0 1 0  |
!                  | 0 0 -1 | 
! ss=nfree dot e3 = (normalized) -1
! nghost=(sin(theta),-cos(theta)) 
! nperp=(-cos(theta),-sin(theta)) 
!
! NOTE IF theta=pi/2:
!  nghost=ss * e3=sign(nfree dot e3) * e3=normalize( (nfree dot e3) x e3 )
!  e3=nsolid x (nsolid x nfree)
!  nperp=ss * nsolid

       ! e2=ns cross nfree
       ! e3 is orthogonal to contact line projected to solid surface.
       ! e3=ns cross e2
       ! ss=sign(nfree dot e3)
      do i=1,SDIM
       nghost(i)=ss*sin_angle*e3(i)-cos_angle*nsolid_new(i)
       nperp(i)=ss*sin_angle*nsolid_new(i)+cos_angle*e3(i)
      enddo 

      if (SDIM.eq.3) then
       dist=sqrt( nghost(1)**2+nghost(2)**2+nghost(SDIM)**2+1.0D-14 )
      else if (SDIM.eq.2) then
       dist=sqrt( nghost(1)**2+nghost(2)**2+1.0D-14 )
      else
       print *,"dimension bust"
       stop
      endif
      do i=1,SDIM
       nghost(i)=nghost(i)/dist
      enddo

      if (SDIM.eq.3) then
       dist=sqrt( nperp(1)**2+nperp(2)**2+nperp(SDIM)**2+1.0D-14 )
      else if (SDIM.eq.2) then
       dist=sqrt( nperp(1)**2+nperp(2)**2+1.0D-14 )
      else
       print *,"dimension bust"
       stop
      endif
      do i=1,SDIM
       nperp(i)=nperp(i)/dist
      enddo

      return
      end subroutine ghostnormal

      subroutine safe_data_single(i,j,k,datafab,data_out)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: i,j,k
      REAL_T, intent(in), pointer :: datafab(D_DECL(:,:,:))
      REAL_T, intent(out) :: data_out
      INTEGER_T datalo,datahi
      INTEGER_T dir
      INTEGER_T idata(3)

      idata(1)=i
      idata(2)=j
      idata(3)=k
      do dir=1,SDIM
       datalo=LBOUND(datafab,dir)
       datahi=UBOUND(datafab,dir)
       if (idata(dir).lt.datalo) then
        idata(dir)=datalo
       endif
       if (idata(dir).gt.datahi) then
        idata(dir)=datahi
       endif
      enddo ! dir=1..sdim
      data_out=datafab(D_DECL(idata(1),idata(2),idata(3)))

      return
      end subroutine safe_data_single

      subroutine safe_data(i,j,k,n,datafab,data_out)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: i,j,k,n
      REAL_T, intent(in), pointer :: datafab(D_DECL(:,:,:),:)
      REAL_T, intent(out) :: data_out
      INTEGER_T datalo,datahi
      INTEGER_T dir
      INTEGER_T idata(3)

      idata(1)=i
      idata(2)=j
      idata(3)=k
      do dir=1,SDIM
       datalo=LBOUND(datafab,dir)
       datahi=UBOUND(datafab,dir)
       if (idata(dir).lt.datalo) then
        idata(dir)=datalo
       endif
       if (idata(dir).gt.datahi) then
        idata(dir)=datahi
       endif
      enddo ! dir=1..sdim
      data_out=datafab(D_DECL(idata(1),idata(2),idata(3)),n)

      return
      end subroutine safe_data


      subroutine safe_data_index(i,j,k,i_safe,j_safe,k_safe,datafab)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: i,j,k
      INTEGER_T, intent(out) :: i_safe,j_safe,k_safe
      REAL_T, intent(in), pointer :: datafab(D_DECL(:,:,:),:)
      INTEGER_T datalo,datahi
      INTEGER_T dir
      INTEGER_T idata(3)

      idata(1)=i
      idata(2)=j
      idata(3)=k
      do dir=1,SDIM
       datalo=LBOUND(datafab,dir)
       datahi=UBOUND(datafab,dir)
       if (idata(dir).lt.datalo) then
        idata(dir)=datalo
       endif
       if (idata(dir).gt.datahi) then
        idata(dir)=datahi
       endif
      enddo ! dir=1..sdim
      i_safe=idata(1)
      j_safe=idata(2)
      k_safe=idata(3)

      return
      end subroutine safe_data_index


      ! Added by Guibo 11-12-2012
      subroutine dumpstring(instring)
      implicit none
      character*80 instring
      INTEGER_T len, ii, I

      len = LEN_TRIM(instring)
      do ii = 1,len
       I = ICHAR(instring(ii:ii))
       write(11) I
      end do
      write(11) 0

      return
      end subroutine dumpstring

      INTEGER_T function is_compressible_mat(im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: im
      INTEGER_T nmat

      nmat=num_materials
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im out of range"
       stop
      endif

      if ((fort_material_type(im).ge.1).and. &
          (fort_material_type(im).le.MAX_NUM_EOS)) then
       is_compressible_mat=1
      else if ((fort_material_type(im).eq.0).or. &
               (fort_material_type(im).eq.999)) then
       is_compressible_mat=0
      else
       print *,"fort_material_type invalid"
       stop
       is_compressible_mat=0
      endif

      return
      end function is_compressible_mat

      subroutine debug_EOS(im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T im
      INTEGER_T verbose_EOS
      INTEGER_T mat_type
      INTEGER_T nden
      INTEGER_T i,iden
      REAL_T temperature,denlo,denhi,den
      REAL_T internal_energy
      REAL_T pressure
      REAL_T soundsqr
      REAL_T massfrac_parm(num_species_var+1)
      character*2 im_str
      character*4 filename4
      character*5 filename5


      verbose_EOS=0
      mat_type=fort_material_type(im)
      if ((mat_type.gt.0).and.(mat_type.le.MAX_NUM_EOS)) then
       if (verbose_EOS.eq.1) then
        temperature=fort_tempconst(im)
        print *,"EOS and C2 files when temperature=",temperature
        denlo=fort_density_floor(im)
        denhi=fort_density_ceiling(im)
        nden=1000
        print *,"EOS and C2 files when temperature=",temperature
        print *,"im=",im
        print *,"denlo=",denlo
        print *,"denhi=",denhi
        print *,"nden=",nden

        write(im_str,'(I2)') im
        do i=1,2
         if (im_str(i:i).eq.' ') then
          im_str(i:i)='0'
         endif
        enddo
        write(filename5,'(A3,A2)') 'EOS',im_str
        open(unit=11,file=filename5)
        write(filename4,'(A2,A2)') 'C2',im_str
        open(unit=12,file=filename4)

        do iden=0,nden
         den=denlo+iden*(denhi-denlo)/nden
         call init_massfrac_parm(den,massfrac_parm,im)
         call INTERNAL_material(den,massfrac_parm, &
          temperature,internal_energy, &
          mat_type,im)
         call EOS_material(den,massfrac_parm, &
          internal_energy,pressure,mat_type,im)
         call SOUNDSQR_material(den,massfrac_parm, &
          internal_energy,soundsqr,mat_type,im)
         write (11,*) den,pressure
         write (12,*) den,soundsqr
        enddo ! iden=0..nden

        close(11)
        close(12)
       else if (verbose_EOS.eq.0) then
        ! do nothing
       else
        print *,"verbose_EOS invalid"
        stop
       endif
      else
       print *,"mat_type invalid"
       stop
      endif 
       
      end subroutine debug_EOS

      REAL_T function get_user_temperature(time,bcflag,im)
      use probcommon_module
      IMPLICIT NONE

      REAL_T, intent(in) :: time
      REAL_T startup_time
      INTEGER_T, intent(in) :: bcflag,im
      INTEGER_T nmat

      nmat=num_materials
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im out of range"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in get_user_temperature"
      else if (time.lt.zero) then
       print *,"time invalid in get_user_temperature"
       stop
      else
       print *,"time bust in get_user_temperature"
       stop
      endif

      if ((fort_initial_temperature(im).le.zero).or. &
          (fort_tempconst(im).le.zero)) then
       print *,"fort_initial_temperature(im) invalid or"
       print *,"fort_tempconst(im) invalid"
       stop
      endif

       ! in: get_user_temperature
      if ((probtype.eq.55).and.(axis_dir.eq.5)) then
       startup_time=zero
       if ((startup_time.eq.zero).or.(im.ne.4)) then
        if (bcflag.eq.0) then
         get_user_temperature=fort_initial_temperature(im)
        else if (bcflag.eq.1) then
         get_user_temperature=fort_tempconst(im)
        else
         print *,"bcflag invalid"
         stop
        endif
       else if ((startup_time.gt.zero).and.(im.eq.4)) then
        if (time.ge.startup_time) then
         get_user_temperature=fort_tempconst(im)
        else if ((time.ge.zero).and.(time.le.startup_time)) then
         get_user_temperature=fort_initial_temperature(im)+ &
          (fort_tempconst(im)-fort_initial_temperature(im))* &
          time/startup_time
        else
         print *,"time invalid in get_user_temperature"
         stop
        endif
       else
        print *,"startup_time or im invalid"
        stop
       endif
      else

       if (bcflag.eq.0) then
        get_user_temperature=fort_initial_temperature(im)
       else if (bcflag.eq.1) then
        get_user_temperature=fort_tempconst(im)
       else
        print *,"bcflag invalid"
        stop
       endif

      endif

      return
      end function get_user_temperature


        ! temperature initialized with a default value
        ! called from denBC and process_initdata
        ! probtype.eq.59, probtype.eq.55, or probtype.eq.710
        ! bcflag==0 => called from initdata
        ! bcflag==1 => called from denbc
      subroutine outside_temperature(time,x,y,z,temperature,im,bcflag)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: im,bcflag
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: x,y,z
      REAL_T, intent(out) :: temperature
      REAL_T :: temp_slope
      REAL_T substrate_height
      REAL_T ice_vertical
      INTEGER_T im_solid_temperature
      REAL_T zprime

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in outside_temperature"
      else if (time.lt.zero) then
       print *,"time invalid in outside_temperature"
       stop
      else
       print *,"time bust in outside_temperature"
       stop
      endif

      im_solid_temperature=im_solid_primary()

      if (SDIM.eq.2) then
       if (abs(z-y).gt.VOFTOL) then
        print *,"expecting z=y in 2d routine outside_temperature"
        stop
       endif
      endif
      zprime=z
 
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid62"
       stop
      endif
      if ((bcflag.ne.0).and. & ! called from initdata
          (bcflag.ne.1)) then  ! called from denBC (density, temp, spec bc)
       print *,"bcflag invalid"
       stop
      endif

       ! block of ice melting on substrate
      if ((probtype.eq.59).or. & ! called from PROB.F90
          (probtype.eq.414).or. &
          (probtype.eq.2001)) then! called from MITSUHIRO_MELTING.F90 

       if (SDIM.eq.2) then
        substrate_height=yblob2
        ice_vertical=yblob
       else if (SDIM.eq.3) then 
        substrate_height=zblob2
        ice_vertical=zblob
       else
        print *,"dimension bust"
        stop
       endif

       if (substrate_height.gt.zero) then 

        if (abs(substrate_height-(ice_vertical-half*radblob)).gt.VOFTOL) then
         print *,"init bot of water+ice block should coincide with substrate"
         stop
        endif
        if (radblob3.ge.radblob) then
         print *,"already melted portion exceeds original ice dimensions"
         stop
        endif

         ! called from initdata
        if (bcflag.eq.0) then
         if ((im.eq.1).or.(im.eq.2)) then ! water or air
          temperature=get_user_temperature(time,bcflag,im)
         else if (im.eq.3) then ! ice 
          temperature=get_user_temperature(time,bcflag,im)
         else if (im.eq.4) then ! substrate
          if (im_solid_temperature.ne.4) then
           print *,"im_solid_temperature invalid"
           stop
          endif
          if (num_materials.lt.4) then
           print *,"num_materials invalid"
           stop
          endif
           
          if (zprime.ge.yblob2) then ! above the substrate
           temperature=get_user_temperature(time,bcflag,3) ! ice
          else if ((zprime.ge.zero).and. &
                   (zprime.le.yblob2)) then
           temperature= &
            get_user_temperature(time,bcflag,im_solid_temperature)+ &
            (get_user_temperature(time,bcflag,3)- &
             get_user_temperature(time,bcflag,im_solid_temperature))* &
            zprime/yblob2 
          else if (zprime.le.zero) then
           temperature= &
            get_user_temperature(time,bcflag,im_solid_temperature) !substrate
          else
           print *,"zprime failure"
           stop
          endif
         else
          print *,"im invalid63"
          stop
         endif
        else if (bcflag.eq.1) then ! called from denBC
         if (im_solid_temperature.ne.4) then
          print *,"im_solid_temperature invalid"
          stop
         endif
          ! radblob3=thickness of underside of block already melted.
         if (zprime.ge.yblob2+radblob3) then
          temperature=get_user_temperature(time,bcflag,3) ! ice
         else if (zprime.ge.yblob2) then
          temperature=get_user_temperature(time,bcflag,3) ! water
         else if ((zprime.ge.zero).and. &
                  (zprime.le.yblob2)) then
          temperature= &
           get_user_temperature(time,bcflag,im_solid_temperature)+ &
           (get_user_temperature(time,bcflag,3)- &
            get_user_temperature(time,bcflag,im_solid_temperature))* &
            zprime/yblob2
         else if (zprime.le.zero) then
          ! substrate
          temperature=get_user_temperature(time,bcflag,im_solid_temperature) 
         else
          print *,"zprime failure"
          stop
         endif
        else
         print *,"bcflag invalid"
         stop
        endif
       endif ! yblob2>0

       ! in: outside_temperature
      else if (probtype.eq.55) then

       if (axis_dir.eq.5) then ! freezing: solid, ice, water, air
         ! substrate: 0<y<yblob2
        if (yblob2.gt.zero) then  

         if (zblob4.eq.yblob2) then !transition thermal layer
          ! do nothing
         else if (zblob4.eq.zero) then ! no transition, T=T_substrate
          ! do nothing
         else
          print *,"zblob4 invalid"
          stop
         endif

          ! called from initdata
         if (bcflag.eq.0) then
          if ((im.eq.1).or.(im.eq.2)) then ! water or air
           temperature=get_user_temperature(time,bcflag,im)
          else if (im.eq.3) then ! ice 
           temperature=get_user_temperature(time,bcflag,im)
          else if (im.eq.4) then ! substrate
           if (im_solid_temperature.ne.4) then
            print *,"im_solid_temperature invalid"
            stop
           endif
           if (num_materials.lt.4) then
            print *,"num_materials invalid"
            stop
           endif
           if (zprime.ge.yblob2) then
            temperature=get_user_temperature(time,bcflag,3) ! ice
           else if ((zprime.ge.zero).and. &
                    (zprime.le.yblob2)) then

            if (zblob4.eq.zero) then
             temperature=get_user_temperature(time,bcflag,im_solid_temperature)
            else if (zblob4.eq.yblob2) then
             temperature= &
              get_user_temperature(time,bcflag,im_solid_temperature)+ &
              (get_user_temperature(time,bcflag,3)- &
               get_user_temperature(time,bcflag,im_solid_temperature))* &
              zprime/yblob2 
            else
             print *,"zblob4 invalid"
             stop
            endif

           else if (zprime.le.zero) then
            temperature= &
             get_user_temperature(time,bcflag,im_solid_temperature) !substrate
           else
            print *,"zprime failure"
            stop
           endif
          else
           print *,"im invalid64"
           stop
          endif
         else if (bcflag.eq.1) then ! called from denBC

          if (im_solid_temperature.ne.4) then
           print *,"im_solid_temperature invalid"
           stop
          endif
          if (zprime.ge.yblob2+radblob3) then
           temperature=get_user_temperature(time,bcflag,1) ! water
          else if (zprime.ge.yblob2) then
           temperature=get_user_temperature(time,bcflag,3) ! ice
          else if ((zprime.ge.zero).and. &
                   (zprime.le.yblob2)) then

           if (zblob4.eq.zero) then
            temperature=get_user_temperature(time,bcflag,im_solid_temperature)
           else if (zblob4.eq.yblob2) then
            temperature= &
             get_user_temperature(time,bcflag,im_solid_temperature)+ &
             (get_user_temperature(time,bcflag,3)- &
              get_user_temperature(time,bcflag,im_solid_temperature))* &
             zprime/yblob2
           else
            print *,"zblob4 invalid"
            stop
           endif

          else if (zprime.le.zero) then
           ! substrate
           temperature=get_user_temperature(time,bcflag,im_solid_temperature) 
          else
           print *,"zprime failure"
           stop
          endif
         else
          print *,"bcflag invalid"
          stop
         endif
        else if (yblob2.eq.zero) then
         ! do nothing
        else
         print *,"not expecting yblob2<0"
         print *,"probtype,axis_dir,bcflag ",probtype,axis_dir,bcflag
         print *,"im=",im
         print *,"time=",time
         print *,"zprime=",zprime
         print *,"yblob2=",yblob2
         stop
        endif ! yblob2>0

        ! boiling sites problem
        ! For Sato and Niceno problem:
        ! solidheat_flag=0 (diffuse in solid)
        ! zblob3<0.0 => heat source in solid cells that adjoin fluid cells.

       else if ((axis_dir.eq.6).or. &
                (axis_dir.eq.7)) then

         ! liquid, vapor, substrate,  or,
         ! liquid, vapor, gas, substrate
        if (im_solid_temperature.ne.num_materials) then 
         print *,"im_solid_temperature invalid"
         stop
        endif
          ! thermal layer thickness
        if (yblob3.gt.zero) then
         ! do nothing
        else
         print *,"yblob3 invalid"
         stop
        endif
        if (SDIM.eq.2) then
         substrate_height=yblob2
        else if (SDIM.eq.3) then 
         substrate_height=zblob2
        else
         print *,"dimension bust"
         stop
        endif

        if (radblob2.lt.zero) then
         print *,"radblob2 invalid"
        else if (radblob2.eq.zero) then
         ! do nothing
        else if (radblob2.gt.zero) then ! angle of incline
         if (levelrz.eq.1) then
          if ((SDIM.ne.2).or.(xblob2.ne.zero)) then
           print *,"parameters not supported"
           stop
          endif
          zprime=yblob2+(z-yblob2)*cos(radblob2)-x*sin(radblob2)
         else if (levelrz.eq.0) then
          if (SDIM.eq.2) then
           zprime=yblob2+(z-yblob2)*cos(radblob2)-(x-xblob2)*sin(radblob2)
          else if (SDIM.eq.3) then
           if (radblob3.ne.zero) then
            print *,"radblob3.ne.zero is not supported"
            stop
           endif
           zprime=zblob2+(z-zblob2)*cos(radblob2)-(x-xblob2)*sin(radblob2)
          else
           print *,"dimension bust"
           stop
          endif
         else
          print *,"levelrz not supported"
          stop
         endif
        else
         print *,"radblob2 invalid"
         stop
        endif

         ! (xblob2,yblob2,zblob2) is the "center" of the heated plate.
        if (zprime.ge.substrate_height+yblob3) then
         temperature=get_user_temperature(time,bcflag,1)
        else if (zprime.le.substrate_height) then
         temperature=get_user_temperature(time,bcflag,im_solid_temperature)
        else if ((substrate_height.le.zprime).and. &
                 (zprime.le.substrate_height+yblob3)) then
         temp_slope=(get_user_temperature(time,bcflag,im_solid_temperature)- &
                     get_user_temperature(time,bcflag,1))/yblob3
         temperature= &
          get_user_temperature(time,bcflag,im_solid_temperature)- &
          temp_slope*(zprime-substrate_height)
        else
         print *,"zprime or substrate_height invalid"
         stop
        endif
       endif

        ! in: outside_temperature
      else if (probtype.eq.710) then

       if (im_solid_temperature.ne.3) then
        print *,"im_solid_temperature invalid"
        stop
       endif
       temperature=get_user_temperature(time,bcflag,1)
          ! thermal layer thickness
       if (yblob3.le.zero) then
        print *,"yblob3 invalid"
        stop
       endif
       if (zprime.le.yblob3) then
        temperature=get_user_temperature(time,bcflag,im_solid_temperature)
       endif

      else
       print *,"expecting probtype=55, 59, 414, 2001, or 710 "
       print *,"in outside_temperature"
       stop
      endif
 
      return 
      end subroutine outside_temperature



      subroutine get_rigid_velocity( &
        FSI_prescribed_flag, &
        color,dir,vel,xrigid, &
        blob_array,blob_array_size,num_colors)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(out) :: FSI_prescribed_flag
      INTEGER_T, intent(in) :: color,dir,blob_array_size,num_colors
      REAL_T, intent(in) :: xrigid(SDIM)
      REAL_T xrigid3D(3)
      REAL_T, intent(out) :: vel
      REAL_T, intent(in) :: blob_array(blob_array_size)
      INTEGER_T nmat
      INTEGER_T cen_comp,vel_comp,vol_comp
      INTEGER_T ibase
      INTEGER_T veldir,irow
      REAL_T blob_mass_for_velocity(3)
      REAL_T blob_center(3)
      REAL_T phi_N(3)
      REAL_T vel_local(SDIM) 

      nmat=num_materials

      if (blob_array_size.ne.num_colors*num_elements_blobclass) then
       print *,"blob_array_size invalid"
       stop
      endif

      FSI_prescribed_flag=0

      cen_comp=BLB_CEN_ACT
      vel_comp=BLB_VEL 
      vol_comp=BLB_MASS_VEL
      
      if ((color.ge.1).and.(color.le.num_colors)) then

       if ((dir.ge.1).and.(dir.le.SDIM)) then

        ibase=(color-1)*num_elements_blobclass
        do veldir=1,3
         blob_mass_for_velocity(veldir)=blob_array(ibase+vol_comp+veldir)
        enddo
        do veldir=1,3
         xrigid3D(veldir)=zero
         blob_center(veldir)=zero
        enddo
        do veldir=1,SDIM
         blob_center(veldir)=blob_array(ibase+cen_comp+veldir)
         xrigid3D(veldir)=xrigid(veldir)
        enddo
        do veldir=1,SDIM
         vel_local(veldir)=zero
        enddo 
        do irow=1,2*SDIM
         call init_basis(blob_center,xrigid3D,irow,phi_N)
         do veldir=1,SDIM
          if (blob_mass_for_velocity(2).gt.zero) then ! ice touches a solid.
           vel_local(veldir)=vel_local(veldir)+ &
            blob_array(ibase+vel_comp+2*SDIM+irow)*phi_N(veldir)
           FSI_prescribed_flag=1
          else if (blob_mass_for_velocity(2).eq.zero) then
           if (blob_mass_for_velocity(1).gt.zero) then
            vel_local(veldir)=vel_local(veldir)+ &
             blob_array(ibase+vel_comp+irow)*phi_N(veldir)
           else if (blob_mass_for_velocity(1).eq.zero) then
            if (blob_mass_for_velocity(3).gt.zero) then
             vel_local(veldir)=vel_local(veldir)+ &
              blob_array(ibase+vel_comp+2*(2*SDIM)+irow)*phi_N(veldir)
            else if (blob_mass_for_velocity(3).eq.zero) then
             ! do nothing
            else
             print *,"blob_mass_for_velocity(3) invalid"
             stop
            endif
           else
            print *,"blob_mass_for_velocity(1) invalid"
            stop
           endif
          else
           print *,"blob_mass_for_velocity(2) invalid"
           stop
          endif
         enddo ! veldir=1..sdim
        enddo ! irow=1,2*SDIM
        if ((dir.ge.1).and.(dir.le.SDIM)) then
         vel=vel_local(dir)
        else
         print *,"dir invalid"
         stop
        endif
       else
        print *,"dir invalid"
        stop
       endif
      else
       print *,"color invalid"
       stop
      endif

      return
      end subroutine get_rigid_velocity


       ! nbasis=1..2 * sdim 
      subroutine init_basis(rigid_centroid,x0,nbasis,phi_N)
      IMPLICIT NONE

      INTEGER_T nbasis
      REAL_T rigid_centroid(3)
      REAL_T x0(3)
      REAL_T phi_N(3)
      INTEGER_T dir

      if ((nbasis.lt.1).or.(nbasis.gt.2*SDIM)) then
       print *,"nbasis invalid"
       stop
      endif

      do dir=1,3
       phi_N(dir)=zero
      enddo
      if ((nbasis.ge.1).and.(nbasis.le.SDIM)) then
        phi_N(nbasis)=one
      else if (nbasis.eq.SDIM+1) then
        phi_N(1)=x0(2)-rigid_centroid(2)
        phi_N(2)=-(x0(1)-rigid_centroid(1)) 
      else if ((nbasis.eq.SDIM+2).and.(SDIM.eq.3)) then
        phi_N(1)=x0(3)-rigid_centroid(3)
        phi_N(3)=-(x0(1)-rigid_centroid(1)) 
      else if ((nbasis.eq.SDIM+3).and.(SDIM.eq.3)) then
        phi_N(2)=x0(3)-rigid_centroid(3)
        phi_N(3)=-(x0(2)-rigid_centroid(2)) 
      else if ((nbasis.eq.2*SDIM).and.(SDIM.eq.2)) then
       ! do nothing
      else
        print *,"nbasis invalid"
        print *,"nbasis= ",nbasis
        stop
      endif

      return
      end subroutine init_basis


      subroutine check_stencil(i,j,k,loface,hiface,velbc, &
        sten_flag)
      IMPLICIT NONE

      INTEGER_T dir,i,j,k
      INTEGER_T loface(SDIM),hiface(SDIM)
      INTEGER_T velbc(SDIM,2)
      INTEGER_T sten_flag,idx,side

      do dir=1,SDIM
       if (dir.eq.1) then
        idx=i
       else if (dir.eq.2) then
        idx=j
       else if ((dir.eq.3).and.(SDIM.eq.3)) then
        idx=k
       else
        print *,"dir invalid check stencil"
        stop
       endif
     
       side=0
       if (idx.lt.loface(dir)) then
        side=1
       else if (idx.gt.hiface(dir)) then
        side=2
       endif

       if ((side.eq.1).or.(side.eq.2)) then 
        if (velbc(dir,side).ne.INT_DIR) then
         sten_flag=0
        else if (velbc(dir,side).eq.INT_DIR) then
         ! do nothing
        endif
       else if (side.eq.0) then
        ! do nothing
       else
        print *,"side invalid"
        stop
       endif
      enddo ! dir

      return
      end subroutine check_stencil
   
      subroutine get_crossing(xcross,xcen,xside,LScen,LSside)
      IMPLICIT NONE

      REAL_T xcross,xcen,xside,LScen,LSside

      if ((LScen.eq.zero).and.(LSside.eq.zero)) then
       xcross=half*(xcen+xcross)
      else if (LScen.eq.zero) then
       xcross=xcen
      else if (LSside.eq.zero) then
       xcross=xside
      else if (abs(LSside)+abs(LScen).gt.zero) then
       xcross= &
        (abs(LSside)*xcen+abs(LScen)*xside)/ &
        (abs(LSside)+abs(LScen))
      else
       print *,"LSside or LScen bust"
       stop
      endif 

      end subroutine get_crossing

      subroutine get_physical_dist(xtarget,phi,gradphi,newphi)
      use probcommon_module
   
      IMPLICIT NONE

      REAL_T xtarget(SDIM)
      REAL_T xfoot(SDIM)
      REAL_T xtargetRT(SDIM)
      REAL_T xfootRT(SDIM)
      REAL_T gradphi(SDIM) 
      REAL_T phi,newphi
      REAL_T mag
      INTEGER_T dir

      if (levelrz.eq.0) then
       newphi=phi
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       newphi=phi
      else if (levelrz.eq.3) then
       mag=zero
       do dir=1,SDIM
        mag=mag+gradphi(dir)**2
       enddo
       mag=sqrt(mag)

       if (mag.gt.zero) then
        do dir=1,SDIM
         xfoot(dir)=xtarget(dir)-phi*gradphi(dir)/mag
        enddo
        call RT_transform(xfoot,xfootRT)
        call RT_transform(xtarget,xtargetRT)
        newphi=zero
        do dir=1,SDIM
         newphi=newphi+(xtargetRT(dir)-xfootRT(dir))**2
        enddo
        newphi=sign_funct(phi)*sqrt(newphi)
       else if (mag.eq.zero) then
        newphi=phi
       else
        print *,"mag invalid in get_physical_dist"
        stop
       endif
      else
       print *,"levelrz invalid"
       stop
      endif

      return
      end subroutine get_physical_dist

       ! mag=|grad phi|
      subroutine prepare_normal(gradphi,rval,mag)
      use probcommon_module

      IMPLICIT NONE

      REAL_T gradphi(SDIM)
      REAL_T rval
      REAL_T mag
      INTEGER_T dir

      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.3) then
       if (rval.gt.zero) then
        gradphi(2)=gradphi(2)/rval
       else
        print *,"rval invalid in prepare_normal; rval=",rval
        stop
       endif
      else
       print *,"levelrz invalid"
       stop
      endif

      mag=zero
      do dir=1,SDIM
       mag=mag+gradphi(dir)**2
      enddo
      mag=sqrt(mag)
      if (mag.gt.zero) then
       do dir=1,SDIM
        gradphi(dir)=gradphi(dir)/mag
       enddo
      endif

      end subroutine prepare_normal
 
      subroutine normalize_vector(vin)
      IMPLICIT NONE

      REAL_T, intent(inout) :: vin(SDIM)
      INTEGER_T dir
      REAL_T mag

      mag=zero
      do dir=1,SDIM
       mag=mag+vin(dir)**2
      enddo
      mag=sqrt(mag)
      if (mag.gt.zero) then
       do dir=1,SDIM
        vin(dir)=vin(dir)/mag
       enddo
      else if (mag.eq.zero) then
       do dir=1,SDIM
        vin(dir)=zero
       enddo
      else
       print *,"mag bust"
       stop
      endif

      return 
      end subroutine normalize_vector

      subroutine normalize_LS_normals(nmat,LS)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(inout) :: LS(nmat*(1+SDIM))
      INTEGER_T :: im,dir
      REAL_T :: local_normal(SDIM)

      do im=1,nmat
       do dir=1,SDIM
        local_normal(dir)=LS(nmat+(im-1)*SDIM+dir)
       enddo
       call normalize_vector(local_normal)
       do dir=1,SDIM
        LS(nmat+(im-1)*SDIM+dir)=local_normal(dir)
       enddo
      enddo ! im=1..nmat

      return 
      end subroutine normalize_LS_normals

      subroutine crossprod2d(a,b,c)
      IMPLICIT NONE
      
      REAL_T a(SDIM),b(SDIM),c(3)

      c(1)=zero
      c(2)=zero
      c(3)=a(1)*b(2)-b(1)*a(2)

      return
      end subroutine crossprod2d

      subroutine crossprod(a,b,c)
      IMPLICIT NONE

      REAL_T, intent(in) :: a(3),b(3)
      REAL_T, intent(out) :: c(3)

      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=-(a(1)*b(3)-b(1)*a(3))
      c(3)=a(1)*b(2)-b(1)*a(2)

      return
      end subroutine crossprod



      subroutine get_crse_index(i,j,k,ic,jc,kc,dir)
      IMPLICIT NONE

      INTEGER_T i,j,k,ic,jc,kc,dir

      ic=i
      jc=j
      kc=k
      if (dir.eq.SDIM) then
       if (i.ge.0) then
        ic=i/2
       else
        ic=(i-1)/2
       endif
       if (j.ge.0) then
        jc=j/2
       else
        jc=(j-1)/2
       endif
       if (k.ge.0) then
        kc=k/2
       else
        kc=(k-1)/2
       endif
      else if (dir.eq.0) then
       if (i.ge.0) then
        ic=i/2
       else
        ic=(i-1)/2
       endif
      else if (dir.eq.1) then
       if (j.ge.0) then
        jc=j/2
       else
        jc=(j-1)/2
       endif
      else if (dir.eq.SDIM-1) then
       if (k.ge.0) then
        kc=k/2
       else
        kc=(k-1)/2
       endif
      else
       print *,"dir invalid get crse index"
       stop
      endif

      return
      end subroutine get_crse_index

      subroutine matrix_solve(AA,xx,bb,matstatus,numelem)
      IMPLICIT NONE
      INTEGER_T numelem
      REAL_T AA(numelem,numelem)
      REAL_T xx(numelem)
      REAL_T bb(numelem)
      REAL_T alpha,holdvalue
      INTEGER_T i,j,k,holdj,matstatus
      REAL_T rowsum,maxrowsum

      matstatus=1

        ! first we normalize the matrix system by the infinity norm
        ! of AA
        ! this means the maximum eigenvalue of the normalized
        ! system has magnitude <2.
      maxrowsum=zero
      do i=1,numelem
       rowsum=zero
       do j=1,numelem
        rowsum=rowsum+abs(AA(i,j))
       enddo
       if (rowsum.gt.maxrowsum) then
        maxrowsum=rowsum
       endif
       if (rowsum.le.zero) then
        matstatus=0
       endif
      enddo
      if (maxrowsum.le.zero) then
       matstatus=0
      endif
      
      if (matstatus.eq.1) then
 
        do i=1,numelem
         do j=1,numelem
          AA(i,j)=AA(i,j)/maxrowsum
         enddo
         bb(i)=bb(i)/maxrowsum
        enddo

        do i=1,numelem-1
         holdj=i
         holdvalue=abs(AA(i,i))
         do j=i+1,numelem 
          if (abs(AA(j,i)).gt.holdvalue) then
           holdj=j
           holdvalue=abs(AA(j,i))
          endif
         enddo
         if (holdj.ne.i) then
          do j=i,numelem
           holdvalue=AA(i,j)
           AA(i,j)=AA(holdj,j)
           AA(holdj,j)=holdvalue
          enddo
         endif
         holdvalue=bb(i)
         bb(i)=bb(holdj)
         bb(holdj)=holdvalue
         if (abs(AA(i,i)).lt.1.0D-10) then
          matstatus=0
         else
          do j=i+1,numelem
           alpha=AA(j,i)/AA(i,i)
           AA(j,i)=zero
           do k=i+1,numelem
            AA(j,k)=AA(j,k)-alpha*AA(i,k)
           enddo
           bb(j)=bb(j)-alpha*bb(i)
          enddo
         endif
        enddo

        do i=numelem,1,-1
         if (matstatus.eq.1) then
          holdvalue=bb(i)
          do j=i+1,numelem
           holdvalue=holdvalue-AA(i,j)*xx(j)
          enddo
          if (abs(AA(i,i)).lt.1.0D-10) then
           matstatus=0
          else
           xx(i)=holdvalue/AA(i,i)
          endif
         else if (matstatus.eq.0) then
          ! do nothing
         else
          print *,"matstatus invalid"
          stop
         endif
        enddo  ! back solve
     
      else if (matstatus.eq.0) then
       ! do nothing
      else
       print *,"matstatus invalid"
       stop
      endif

      return
      end subroutine matrix_solve

      subroutine print_matrix(AA,numelem)
      IMPLICIT NONE
      INTEGER_T numelem
      REAL_T AA(numelem,numelem)

      INTEGER_T i

      do i=1,numelem
       print *,AA(i,:)
       print *,"endrow"
      enddo

      return
      end subroutine print_matrix

      subroutine matrix_inverse(AA,xx,matstatus,numelem)
      IMPLICIT NONE
      INTEGER_T numelem
      REAL_T AA(numelem,numelem)
      REAL_T AAhold(numelem,numelem)
      REAL_T xx(numelem,numelem)
      REAL_T bb(numelem,numelem)
      REAL_T alpha,holdvalue
      INTEGER_T i,j,k,holdj,matstatus

      do i=1,numelem
       do j=1,numelem
        AAhold(i,j)=AA(i,j)
        bb(i,j)=zero
       enddo
       bb(i,i)=one
      enddo
 
      matstatus=1
      do i=1,numelem-1
         holdj=i
         holdvalue=abs(AA(i,i))
         do j=i+1,numelem 
          if (abs(AA(j,i)).gt.holdvalue) then
           holdj=j
           holdvalue=abs(AA(j,i))
          endif
         enddo
         if (holdj.ne.i) then
          do j=i,numelem
           holdvalue=AA(i,j)
           AA(i,j)=AA(holdj,j)
           AA(holdj,j)=holdvalue
          enddo
         endif
         do k=1,numelem
          holdvalue=bb(i,k)
          bb(i,k)=bb(holdj,k)
          bb(holdj,k)=holdvalue
         enddo
         if (abs(AA(i,i)).lt.1.0D-32) then
          matstatus=0
         else
          do j=i+1,numelem
           alpha=AA(j,i)/AA(i,i)
           do k=i,numelem
            AA(j,k)=AA(j,k)-alpha*AA(i,k)
           enddo
           do k=1,numelem
            bb(j,k)=bb(j,k)-alpha*bb(i,k)
           enddo
          enddo
         endif
      enddo

      do i=numelem,1,-1
         if (matstatus.ne.0) then
          do k=1,numelem
           holdvalue=bb(i,k)
           do j=i+1,numelem
            holdvalue=holdvalue-AA(i,j)*xx(j,k)
           enddo
           if (abs(AA(i,i)).lt.1.0D-32) then
            matstatus=0
           else
            xx(i,k)=holdvalue/AA(i,i)
           endif
          enddo
         endif
      enddo

      if (1.eq.1) then
       do i=1,numelem
        do j=1,numelem
         holdvalue=zero
         do k=1,numelem
          holdvalue=holdvalue+AAhold(i,k)*xx(k,j)
         enddo 
         if (i.ne.j) then
          if (abs(holdvalue).gt.1.0D-12) then
           print *,"inverse failed1"
           print *,"AAhold="
           call print_matrix(AAhold,numelem)
           print *,"xx="
           call print_matrix(xx,numelem)
           stop
          endif
         else if (i.eq.j) then
          if (abs(holdvalue-one).gt.1.0D-12) then
           print *,"inverse failed2"
           print *,"AAhold="
           call print_matrix(AAhold,numelem)
           print *,"xx="
           call print_matrix(xx,numelem)
           stop
          endif
         endif
        enddo
       enddo
      endif

      return
      end subroutine matrix_inverse



      subroutine least_squares_interp(npoints,x0,xpos,wts,vals, &
       order,linearflag,coeffs,ncoeffs,sdim_parm)
      IMPLICIT NONE

      INTEGER_T npoints,order,ncoeffs,sdim_parm
      INTEGER_T linearflag
      REAL_T x0(sdim_parm)
      REAL_T xpos(npoints,sdim_parm)
      REAL_T wts(npoints)
      REAL_T vals(npoints)
      REAL_T coeffs(ncoeffs)
      REAL_T AA(ncoeffs,ncoeffs)
      REAL_T BB(ncoeffs)
      INTEGER_T ncoeffs_test
      INTEGER_T orderx,ordery,orderz
      INTEGER_T i,i2,ipoints,idx,dir
      INTEGER_T ipower,jpower,kpower
      REAL_T basisfn(ncoeffs)
      REAL_T basisdir(0:order,sdim_parm)
      INTEGER_T matstatus

      if ((sdim_parm.ne.2).and.(sdim_parm.ne.3)) then
       print *,"sdim_parm invalid in least_squares_interp"
       stop
      endif
      if (npoints.le.0) then
       print *,"npoints invalid"
       stop
      endif
      if (linearflag.eq.1) then
       if (order.ne.1) then
        print *,"order invalid"
        stop
       endif
       if (ncoeffs.ne.sdim_parm+1) then
        print *,"ncoeffs invalid"
        stop
       endif
      else if (linearflag.eq.0) then
       if (order.lt.0) then
        print *,"order invalid"
        stop
       endif
       if (sdim_parm.eq.2) then
        ncoeffs_test=(order+1)*(order+1)
       else if (sdim_parm.eq.3) then
        ncoeffs_test=(order+1)*(order+1)*(order+1)
       else
        print *,"sdim_parm invalid"
        stop
       endif
       if (ncoeffs.ne.ncoeffs_test) then
        print *,"ncoeffs invalid"
        stop
       endif
      else
       print *,"linearflag invalid"
       stop
      endif

       ! order>=0
      orderx=order
      ordery=order
      if (sdim_parm.eq.2) then
       orderz=0
      else if (sdim_parm.eq.3) then
       orderz=order
      else
       print *,"sdim_parm invalid"
       stop
      endif
 
      do i=1,ncoeffs

       BB(i)=zero

       do i2=1,ncoeffs
 
        AA(i,i2)=zero

       enddo

      enddo  
      
! E=sum_i w_i ( sum_j a_j phi_j(x_i) - val_i )^2
! dE/da_k=0 => sum_i w_i ( sum_j a_j phi_j(x_i)phi_k(x_i) )=
!              sum_i w_i ( sum_j phi_k(x_i) val_i )
! A_kj=sum_i w_i phi_j(x_i)phi_k(x_i)
! B_k =sum_i w_i val_i phi_k(x_i)
 
      do ipoints=1,npoints
       
       if (linearflag.eq.1) then
        basisfn(1)=one
        do dir=1,sdim_parm
         basisfn(dir+1)=xpos(ipoints,dir)-x0(dir) 
        enddo
        if (ncoeffs.ne.sdim_parm+1) then
         print *,"ncoeffs invalid"
         stop
        endif
       else if (linearflag.eq.0) then
        do dir=1,sdim_parm
         basisdir(0,dir)=one
         do ipower=1,order
          basisdir(ipower,dir)=basisdir(ipower-1,dir)* &
           (xpos(ipoints,dir)-x0(dir))
         enddo
        enddo
        idx=1
        do ipower=0,orderx
        do jpower=0,ordery
        do kpower=0,orderz
         if (idx.gt.ncoeffs) then
          print *,"idx invalid"
          stop
         endif
         basisfn(idx)= &
          basisdir(ipower,1)*basisdir(jpower,2)*basisdir(kpower,sdim_parm)
         idx=idx+1
        enddo
        enddo
        enddo
       else
        print *,"linearflag invalid"
        stop
       endif

       do i=1,ncoeffs
        BB(i)=BB(i)+wts(ipoints)*vals(ipoints)*basisfn(i)
        do i2=1,ncoeffs
         AA(i,i2)=AA(i,i2)+wts(ipoints)*basisfn(i2)*basisfn(i)
        enddo
       enddo
 
      enddo  ! ipoints

      call matrix_solve(AA,coeffs,BB,matstatus,ncoeffs)

      if (matstatus.eq.1) then
       ! do nothing
      else 
       print *,"matrix solve failed in least squares routine"
       stop
      endif

      return
      end subroutine least_squares_interp

      subroutine derive_dist( &
       xsten,nhalf, &
       dist, &
       i,j,k,im,LS)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: i,j,k,im
      INTEGER_T, intent(in) :: nhalf
      INTEGER_T :: dir
      REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
      REAL_T, intent(inout) :: LS
      REAL_T, intent(in), pointer :: dist(D_DECL(:,:,:),:)
      REAL_T n(SDIM)
      REAL_T nsave(SDIM)
      REAL_T RR,mag

      if (nhalf.ne.3) then
       print *,"nhalf invalid"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid15"
       stop
      endif

      LS=dist(D_DECL(i,j,k),im)

      if ((levelrz.eq.0).or.(levelrz.eq.1)) then
       ! do nothing
      else if (levelrz.eq.3) then
 
       dir=1
       n(dir)=(dist(D_DECL(i+1,j,k),im)-dist(D_DECL(i-1,j,k),im))/ &
              (xsten(2,dir)-xsten(-2,dir))
       dir=2
       n(dir)=(dist(D_DECL(i,j+1,k),im)-dist(D_DECL(i,j-1,k),im))/ &
              (xsten(2,dir)-xsten(-2,dir))
       if (SDIM.eq.3) then
        dir=SDIM
        n(dir)=(dist(D_DECL(i,j,k+1),im)-dist(D_DECL(i,j,k-1),im))/ &
               (xsten(2,dir)-xsten(-2,dir))
       endif

       RR=one
       call prepare_normal(n,RR,mag)
       if (mag.gt.zero) then
        RR=xsten(0,1)
        do dir=1,SDIM
         nsave(dir)=n(dir)
        enddo
        call prepare_normal(nsave,RR,mag)
        if (mag.gt.zero) then
         mag=zero
         do dir=1,SDIM
          mag=mag+n(dir)*nsave(dir)
         enddo
         if (abs(mag).gt.one+VOFTOL) then
          print *,"mag invalid in derive_dist"
          stop
         endif
         if (abs(mag).gt.zero) then
          LS=LS*abs(mag)
         endif
        endif
       endif
      else
       print *,"levelrz invalid"
       stop
      endif
      
      return
      end subroutine derive_dist

       ! grid_type=-1..5
      subroutine checkbound(lo,hi, &
       DIMS(data), &
       ngrow,grid_type,id)
      IMPLICIT NONE

      INTEGER_T, intent(in) ::  lo(SDIM), hi(SDIM)
      INTEGER_T, intent(in) ::  DIMDEC(data)
      INTEGER_T, intent(in) ::  ngrow
      INTEGER_T, intent(in) ::  grid_type
      INTEGER_T, intent(in) ::  id

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      INTEGER_T box_type(SDIM)

      INTEGER_T    hidata(SDIM)
      INTEGER_T    lodata(SDIM)
      INTEGER_T    dir2

      hidata(1)=ARG_H1(data)
      hidata(2)=ARG_H2(data)
      lodata(1)=ARG_L1(data)
      lodata(2)=ARG_L2(data)
#if (AMREX_SPACEDIM==3)
      hidata(SDIM)=ARG_H3(data)
      lodata(SDIM)=ARG_L3(data)
#endif
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound id=",id
        print *,"grid_type=",grid_type
        print *,"dir2=",dir2
        stop
       endif
       box_type(dir2)=0
      enddo
       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       if (lo(dir2).lt.0) then
        print *,"lo invalid in checkbound id=",id
        print *,"dir2,dataxlo ",dir2,lodata(dir2)
        print *,"dir2,dataxhi ",dir2,hidata(dir2)
        print *,"dir2,lo,ngrow ",dir2,lo(dir2),ngrow
        print *,"dir2,hi,ngrow ",dir2,hi(dir2),ngrow
        print *,"grid_type=",grid_type
        stop
       endif
      enddo
      if (ngrow.lt.0) then
       print *,"ngrow invalid in checkbound"
       stop
      endif
      if (id.lt.0) then
       print *,"id invalid in checkbound"
       stop
      endif

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound:lo mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"box_type(dir2) ",box_type(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif

      enddo ! dir2=1..SDIM

      return
      end subroutine checkbound

       ! grid_type=-1..5
      subroutine checkbound_array(lo,hi, &
       data_array, &
       ngrow,grid_type,id)
      IMPLICIT NONE

      INTEGER_T, intent(in) ::  lo(SDIM), hi(SDIM)
        ! intent(in) means the pointer cannot be reassigned.
        ! The data itself inherits the intent attribute from the
        ! target.
      REAL_T, intent(in), pointer :: data_array(D_DECL(:,:,:),:)
      INTEGER_T, intent(in) ::  ngrow
      INTEGER_T, intent(in) ::  grid_type
      INTEGER_T, intent(in) ::  id

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      INTEGER_T box_type(SDIM)

      INTEGER_T    hidata(SDIM+1)
      INTEGER_T    lodata(SDIM+1)
      INTEGER_T    dir2

      hidata=UBOUND(data_array)
      lodata=LBOUND(data_array)
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound_array id=",id
        print *,"grid_type=",grid_type
        print *,"dir2=",dir2
        stop
       endif
       box_type(dir2)=0
      enddo
       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       if (lo(dir2).lt.0) then
        print *,"lo invalid in checkbound_array id=",id
        print *,"dir2,dataxlo ",dir2,lodata(dir2)
        print *,"dir2,dataxhi ",dir2,hidata(dir2)
        print *,"dir2,lo,ngrow ",dir2,lo(dir2),ngrow
        print *,"dir2,hi,ngrow ",dir2,hi(dir2),ngrow
        print *,"grid_type=",grid_type
        stop
       endif
      enddo
      if (ngrow.lt.0) then
       print *,"ngrow invalid in checkbound_array"
       stop
      endif
      if (id.lt.0) then
       print *,"id invalid in checkbound_array"
       stop
      endif

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound_array:lo mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"box_type(dir2) ",box_type(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif

      enddo ! dir2=1..SDIM

      return
      end subroutine checkbound_array


       ! grid_type=-1..5
      subroutine checkbound_array_INTEGER(lo,hi, &
       data_array, &
       ngrow,grid_type,id)
      IMPLICIT NONE

      INTEGER_T, intent(in) ::  lo(SDIM), hi(SDIM)
        ! intent(in) means the pointer cannot be reassigned.
        ! The data itself inherits the intent attribute from the
        ! target.
      INTEGER_T, intent(in), pointer :: data_array(D_DECL(:,:,:),:)
      INTEGER_T, intent(in) ::  ngrow
      INTEGER_T, intent(in) ::  grid_type
      INTEGER_T, intent(in) ::  id

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      INTEGER_T box_type(SDIM)

      INTEGER_T    hidata(SDIM+1)
      INTEGER_T    lodata(SDIM+1)
      INTEGER_T    dir2

      hidata=UBOUND(data_array)
      lodata=LBOUND(data_array)
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound_array id=",id
        print *,"grid_type=",grid_type
        print *,"dir2=",dir2
        stop
       endif
       box_type(dir2)=0
      enddo
       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       if (lo(dir2).lt.0) then
        print *,"lo invalid in checkbound_array id=",id
        print *,"dir2,dataxlo ",dir2,lodata(dir2)
        print *,"dir2,dataxhi ",dir2,hidata(dir2)
        print *,"dir2,lo,ngrow ",dir2,lo(dir2),ngrow
        print *,"dir2,hi,ngrow ",dir2,hi(dir2),ngrow
        print *,"grid_type=",grid_type
        stop
       endif
      enddo
      if (ngrow.lt.0) then
       print *,"ngrow invalid in checkbound_array"
       stop
      endif
      if (id.lt.0) then
       print *,"id invalid in checkbound_array"
       stop
      endif

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound_array:lo mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"box_type(dir2) ",box_type(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif

      enddo ! dir2=1..SDIM

      return
      end subroutine checkbound_array_INTEGER



       ! grid_type=-1..5
      subroutine checkbound_array1(lo,hi, &
       data_array1, &
       ngrow,grid_type,id)
      IMPLICIT NONE

      INTEGER_T, intent(in) ::  lo(SDIM), hi(SDIM)
        ! intent(in) means the pointer cannot be reassigned.
        ! The data itself inherits the intent attribute from the
        ! target.
      REAL_T, intent(in), pointer :: data_array1(D_DECL(:,:,:))
      INTEGER_T, intent(in) ::  ngrow
      INTEGER_T, intent(in) ::  grid_type
      INTEGER_T, intent(in) ::  id

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      INTEGER_T box_type(SDIM)

      INTEGER_T    hidata(SDIM)
      INTEGER_T    lodata(SDIM)
      INTEGER_T    dir2

      hidata=UBOUND(data_array1)
      lodata=LBOUND(data_array1)
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound_array1 id=",id
        print *,"grid_type=",grid_type
        print *,"dir2=",dir2
        stop
       endif
       box_type(dir2)=0
      enddo
       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       if (lo(dir2).lt.0) then
        print *,"lo invalid in checkbound_array1 id=",id
        print *,"dir2,dataxlo ",dir2,lodata(dir2)
        print *,"dir2,dataxhi ",dir2,hidata(dir2)
        print *,"dir2,lo,ngrow ",dir2,lo(dir2),ngrow
        print *,"dir2,hi,ngrow ",dir2,hi(dir2),ngrow
        print *,"grid_type=",grid_type
        stop
       endif
      enddo
      if (ngrow.lt.0) then
       print *,"ngrow invalid in checkbound_array1"
       stop
      endif
      if (id.lt.0) then
       print *,"id invalid in checkbound_array1"
       stop
      endif

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound_array:lo mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"box_type(dir2) ",box_type(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif

      enddo ! dir2=1..SDIM

      return
      end subroutine checkbound_array1


       ! grid_type=-1..5
      subroutine checkbound_int_array1(lo,hi, &
       data_array1, &
       ngrow,grid_type,id)
      IMPLICIT NONE

      INTEGER_T, intent(in) ::  lo(SDIM), hi(SDIM)
        ! intent(in) means the pointer cannot be reassigned.
        ! The data itself inherits the intent attribute from the
        ! target.
      INTEGER_T, intent(in), pointer :: data_array1(D_DECL(:,:,:))
      INTEGER_T, intent(in) ::  ngrow
      INTEGER_T, intent(in) ::  grid_type
      INTEGER_T, intent(in) ::  id

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      INTEGER_T box_type(SDIM)

      INTEGER_T    hidata(SDIM)
      INTEGER_T    lodata(SDIM)
      INTEGER_T    dir2

      hidata=UBOUND(data_array1)
      lodata=LBOUND(data_array1)
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound_array1 id=",id
        print *,"grid_type=",grid_type
        print *,"dir2=",dir2
        stop
       endif
       box_type(dir2)=0
      enddo
       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       if (lo(dir2).lt.0) then
        print *,"lo invalid in checkbound_array1 id=",id
        print *,"dir2,dataxlo ",dir2,lodata(dir2)
        print *,"dir2,dataxhi ",dir2,hidata(dir2)
        print *,"dir2,lo,ngrow ",dir2,lo(dir2),ngrow
        print *,"dir2,hi,ngrow ",dir2,hi(dir2),ngrow
        print *,"grid_type=",grid_type
        stop
       endif
      enddo
      if (ngrow.lt.0) then
       print *,"ngrow invalid in checkbound_array1"
       stop
      endif
      if (id.lt.0) then
       print *,"id invalid in checkbound_array1"
       stop
      endif

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound_array:lo mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"box_type(dir2) ",box_type(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif

      enddo ! dir2=1..SDIM

      return
      end subroutine checkbound_int_array1

! find reconstruction from cell averaged data.
! point value lives at the center of the cell in 3D.
      subroutine minmod3D(data,xpos,slopes,intercept)
      IMPLICIT NONE

      REAL_T intercept
      REAL_T data(3,3,3)
      REAL_T xpos(3,3,3,3)
      REAL_T slopes(3)
      INTEGER_T ii,jj,kk,dir
      REAL_T slope_minus,slope_plus
      REAL_T xm,xc,xp,dm,dc,dp

      do dir=1,3
       ii=0
       jj=0
       kk=0
       if (dir.eq.1) then
        ii=1
       else if (dir.eq.2) then
        jj=1
       else if (dir.eq.3) then
        kk=1
       else 
        print *,"dir out of range"
        stop
       endif

       dp=data(2+ii,2+jj,2+kk)
       dc=data(2,2,2)
       dm=data(2-ii,2-jj,2-kk)
       xp=xpos(2+ii,2+jj,2+kk,dir)
       xc=xpos(2,2,2,dir)
       xm=xpos(2-ii,2-jj,2-kk,dir)
       if ((xc-xm.le.zero).or.(xp-xc.le.zero)) then
        print *,"dx became zero in minmod_stencil"
        stop
       endif
       slope_minus=(dc-dm)/(xc-xm) 
       slope_plus=(dp-dc)/(xp-xc)
       
       call minmod(slope_minus,slope_plus,slopes(dir))
      enddo
      intercept=dc

      return 
      end subroutine minmod3D



! find reconstruction from cell averaged data.
! point value lives at the center of the cell in 3D.
      subroutine minmod_stencil(data,xpos,slopes,intercept)
      IMPLICIT NONE

      REAL_T intercept
      REAL_T data(D_DECL(3,3,3))
      REAL_T xpos(D_DECL(3,3,3),SDIM)
      REAL_T slopes(SDIM)
      INTEGER_T ii,jj,kk,dir
      REAL_T slope_minus,slope_plus
      REAL_T xm,xc,xp,dm,dc,dp

      do dir=1,SDIM
       ii=0
       jj=0
       kk=0
       if (dir.eq.1) then
        ii=1
       else if (dir.eq.2) then
        jj=1
       else if ((dir.eq.3).and.(SDIM.eq.3)) then
        kk=1
       else 
        print *,"dir out of range"
        stop
       endif

       dp=data(D_DECL(2+ii,2+jj,2+kk))
       dc=data(D_DECL(2,2,2))
       dm=data(D_DECL(2-ii,2-jj,2-kk))
       xp=xpos(D_DECL(2+ii,2+jj,2+kk),dir)
       xc=xpos(D_DECL(2,2,2),dir)
       xm=xpos(D_DECL(2-ii,2-jj,2-kk),dir)
       if ((xc-xm.le.zero).or.(xp-xc.le.zero)) then
        print *,"dx became zero in minmod_stencil"
        stop
       endif
       slope_minus=(dc-dm)/(xc-xm) 
       slope_plus=(dp-dc)/(xp-xc)
       
       call minmod(slope_minus,slope_plus,slopes(dir))
      enddo
      intercept=dc

      return 
      end subroutine minmod_stencil


      subroutine distfunc(bfact,dx,xsten0,nhalf0, &
        phi0,nn,x,dist,sdim_parm)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: sdim_parm
      INTEGER_T :: dir
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: nhalf0
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: phi0
      REAL_T, intent(out) :: dist
      REAL_T, intent(in) :: nn(sdim_parm)
      REAL_T, intent(in) :: xsten0(-nhalf0:nhalf0,sdim_parm)
      REAL_T, intent(in) :: x(sdim_parm)
 
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid in distfunc"
       stop
      endif
      if ((sdim_parm.ne.3).and.(sdim_parm.ne.2)) then
       print *,"sdim_parm bust distfunc"
       stop
      endif
      dist=phi0
      do dir=1,sdim_parm
       dist=dist+nn(dir)*(x(dir)-xsten0(0,dir))
      enddo

      return
      end subroutine distfunc

        ! -pi/2 < angle < pi/2
      REAL_T function atan_verify(x)
      use probcommon_module
      IMPLICIT NONE
      REAL_T x

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      atan_verify=atan(x)
      if ((atan_verify.le.-half*MOF_PI-VOFTOL).or. &
          (atan_verify.ge.half*MOF_PI+VOFTOL)) then
       print *,"atan out of range"
       stop
      endif
 
      return
      end function atan_verify

! -pi < angle < pi 
! if (x,y) in first quadrant, then 0<theta<pi/2
! if (x,y) in fourth quadrant, then -pi/2<theta<0
! if (x,y) in second quadrant, then pi/2<theta<pi
! if (x,y) in third quadrant, then -pi<theta<-pi/2
      subroutine arctan2(y,x,angle)
      IMPLICIT NONE

      REAL_T y,x,angle

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      if ((y.gt.zero).and.(x.eq.zero)) then
       angle=half*MOF_PI
      else if ((y.lt.zero).and.(x.eq.zero)) then
       angle=-half*MOF_PI
        ! pi/4 <angle<pi/2 (first quadrant)
      else if ((y.gt.zero).and.(x.gt.zero).and.(y.ge.x)) then
       angle=atan_verify(y/x)
        ! pi/2 < angle < 3pi/4 (second quadrant)
      else if ((y.gt.zero).and.(x.lt.zero).and.(y.ge.abs(x))) then
       angle=atan_verify(y/x)+MOF_PI
        ! -pi/2<angle<-pi/4 (fourth quadrant)
      else if ((y.lt.zero).and.(x.gt.zero).and.(abs(y).ge.x)) then
       angle=atan_verify(y/x)
        ! -3pi/4 < angle < -pi/2  (3rd quadrant)
      else if ((y.lt.zero).and.(x.lt.zero).and.(abs(y).ge.abs(x))) then
       angle=atan_verify(y/x)-MOF_PI
      else if ((y.eq.zero).and.(x.gt.zero)) then
       angle=zero
      else if ((y.eq.zero).and.(x.lt.zero)) then
       angle=MOF_PI
       ! 0<angle<pi/4 (1st quadrant)
      else if ((x.gt.zero).and.(y.gt.zero).and.(y.le.x)) then
       angle=atan_verify(y/x)
       ! -pi/4<angle<0 (4th quadrant)
      else if ((x.gt.zero).and.(y.lt.zero).and.(abs(y).le.x)) then
       angle=atan_verify(y/x)
       ! 3pi/4<angle<pi (second quadrant)
      else if ((x.lt.zero).and.(y.gt.zero).and.(y.le.abs(y))) then
       angle=atan_verify(y/x)+MOF_PI
       ! -pi<angle<-3pi/4 (3rd quadrant)
      else if ((x.lt.zero).and.(y.lt.zero).and.(abs(y).le.abs(x))) then
       angle=atan_verify(y/x)-MOF_PI
      else
       angle=zero
      endif

      return
      end subroutine arctan2


      subroutine minmod(x,y,z)
      IMPLICIT NONE
      
      REAL_T x,y,z

      if (x*y.le.zero) then
       z=zero
      else if (abs(x).lt.abs(y)) then
       z=x
      else
       z=y
      endif
 
      return
      end subroutine minmod

      subroutine intersectbox(xlo1,xhi1,xlo2,xhi2,xloint,xhiint,vol)
      use probcommon_module

      IMPLICIT NONE

      REAL_T vol
      REAL_T xlo1(SDIM),xhi1(SDIM)
      REAL_T xlo2(SDIM),xhi2(SDIM)
      REAL_T xloint(SDIM),xhiint(SDIM)
      REAL_T vol1,vol2
      INTEGER_T dir

      vol=one
      vol1=one
      vol2=one
      do dir=1,SDIM
       if ((xlo1(dir).ge.xhi1(dir)).or. &
           (xlo2(dir).ge.xhi2(dir))) then
        print *,"xlo>xhi intersectbox"
        print *,"xlo1 ",xlo1(1),xlo1(2),xlo1(SDIM)
        print *,"xhi1 ",xhi1(1),xhi1(2),xhi1(SDIM)
        print *,"xlo2 ",xlo2(1),xlo2(2),xlo2(SDIM)
        print *,"xhi2 ",xhi2(1),xhi2(2),xhi2(SDIM)
        stop
       endif
       if (xlo1(dir).lt.xlo2(dir)) then
        xloint(dir)=xlo2(dir)
       else
        xloint(dir)=xlo1(dir)
       endif
       if (xhi1(dir).gt.xhi2(dir)) then
        xhiint(dir)=xhi2(dir)
       else
        xhiint(dir)=xhi1(dir)
       endif
       if (xloint(dir).ge.xhiint(dir)) then
        vol=zero
       else
        vol=vol*(xhiint(dir)-xloint(dir))
       endif
       vol1=vol1*(xhi1(dir)-xlo1(dir))
       vol2=vol2*(xhi2(dir)-xlo2(dir))
      enddo

      if ((vol1.le.zero).or.(vol2.le.zero)) then
       print *,"vol1,vol2 invalid"
       stop
      else
       if (vol/vol1.le.VOFTOL) then
        vol=zero
       endif
      endif

      return
      end subroutine intersectbox


      subroutine RT_transform(x,xT)
      use probcommon_module

      IMPLICIT NONE

      REAL_T x(SDIM)
      REAL_T xT(SDIM)
      INTEGER_T dir
     
      do dir=1,SDIM
       xT(dir)=x(dir)
      enddo
      if ((levelrz.eq.0).or.(levelrz.eq.1)) then
       ! do nothing
      else if (levelrz.eq.3) then
       xT(1)=x(1)*cos(x(2))
       xT(2)=x(1)*sin(x(2))
      else
       print *,"levelrz invalid RT transform"
       stop
      endif

      return
      end subroutine RT_transform



      subroutine RT_transform_offset(x,ofs,xT)
      use probcommon_module

      IMPLICIT NONE

      REAL_T, intent(in) :: x(SDIM)
      REAL_T x_ofs(SDIM)
      REAL_T, intent(in) :: ofs(SDIM)
      REAL_T, intent(out) :: xT(SDIM)
      INTEGER_T dir
     
      do dir=1,SDIM
       x_ofs(dir)=x(dir)+ofs(dir)
       xT(dir)=x_ofs(dir)
      enddo
      if ((levelrz.eq.0).or.(levelrz.eq.1)) then
       ! do nothing
      else if (levelrz.eq.3) then
       xT(1)=x_ofs(1)*cos(x_ofs(2))
       xT(2)=x_ofs(1)*sin(x_ofs(2))
      else
       print *,"levelrz invalid RT transform offset"
       stop
      endif

      return
      end subroutine RT_transform_offset


      subroutine get_dxmin(dx,bfact,dxmin)
      use probcommon_module
      use LegendreNodes

      REAL_T,intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T dir
      REAL_T, intent(out) :: dxmin
      REAL_T delta
      REAL_T delta_gauss
      REAL_T RR

      if (bfact.lt.1) then
       print *,"bfact invalid14"
       stop
      endif

      dxmin=zero
      do dir=1,SDIM
       if (bfact.eq.1) then
        delta=dx(dir)
       else if (bfact.gt.1) then
        delta=bfact*(cache_gauss_lobatto(bfact,1,SPTYPE)+one)*dx(dir)/two
        delta_gauss=bfact*(cache_gauss(bfact,0,SPTYPE)+one)*dx(dir)
        if (delta.gt.delta_gauss) then
         delta=delta_gauss
        endif
       else
        print *,"bfact invalid15"
        stop
       endif
       RR=one
       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        ! do nothing
       else if (levelrz.eq.3) then
        if (dir.eq.2) then
         RR=problox
        endif
       else
        print *,"levelrz invalid"
        stop
       endif
       delta=delta*RR

       if (delta.le.zero) then
        print *,"delta invalid in get_dxmin"
        stop
       endif

       if ((delta.lt.dxmin).or.(dir.eq.1)) then
        dxmin=delta
       endif
      enddo ! dir

      end subroutine get_dxmin


      subroutine get_dxmax(dx,bfact,dxmax)
      use probcommon_module
      use LegendreNodes

      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T dir
      REAL_T, intent(out) :: dxmax
      REAL_T delta,RR,xL,xR

      if (bfact.lt.1) then
       print *,"bfact invalid16"
       stop
      endif

      dxmax=zero
      do dir=1,SDIM
       if (bfact.eq.1) then
        delta=dx(dir)
       else if (bfact.gt.1) then
        xL=cache_gauss(bfact,bfact/2-1,SPTYPE)
        xR=cache_gauss(bfact,bfact/2,SPTYPE)
        delta=bfact*(xR-xL)*dx(dir)/two
       else
        print *,"bfact invalid17"
        stop
       endif
       RR=one
       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        ! do nothing
       else if (levelrz.eq.3) then
        if (dir.eq.2) then
         RR=probhix
        endif
       else
        print *,"levelrz invalid"
        stop
       endif
       delta=delta*RR

       if (delta.le.zero) then
        print *,"delta invalid in get_dxmax"
        stop
       endif

       if ((delta.gt.dxmax).or.(dir.eq.1)) then
        dxmax=delta
       endif
      enddo ! dir

      end subroutine get_dxmax


      subroutine get_dxmaxLS(dx,bfact,dxmax)
      use probcommon_module
      use LegendreNodes

      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T dir
      REAL_T, intent(out) :: dxmax
      REAL_T delta,xL,xR

      if (bfact.lt.1) then
       print *,"bfact invalid18"
       stop
      endif

      dxmax=zero
      do dir=1,SDIM
       if (bfact.eq.1) then
        delta=dx(dir)
       else if (bfact.gt.1) then
        xL=cache_gauss(bfact,bfact/2-1,SPTYPE)
        xR=cache_gauss(bfact,bfact/2,SPTYPE)
        delta=bfact*(xR-xL)*dx(dir)/two
       else
        print *,"bfact invalid19"
        stop
       endif

       if (delta.le.zero) then
        print *,"delta invalid get_dxmaxLS"
        stop
       endif

       if ((delta.gt.dxmax).or.(dir.eq.1)) then
        dxmax=delta
       endif
      enddo ! dir

      end subroutine get_dxmaxLS

      subroutine set_dimdec(DIMS(fabdim), &
                      fablo,fabhi,ngrow)
      IMPLICIT NONE

      INTEGER_T, intent(out) :: DIMDEC(fabdim)
      INTEGER_T, intent(in) :: fablo(SDIM)
      INTEGER_T, intent(in) :: fabhi(SDIM)
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T dir

      dir=1
      ARG_L1(fabdim)=fablo(dir)-ngrow
      ARG_H1(fabdim)=fabhi(dir)+ngrow
      dir=2
      ARG_L2(fabdim)=fablo(dir)-ngrow
      ARG_H2(fabdim)=fabhi(dir)+ngrow
      dir=SDIM
#if (AMREX_SPACEDIM==3)
      if (SDIM.eq.3) then
       ARG_L3(fabdim)=fablo(dir)-ngrow
       ARG_H3(fabdim)=fabhi(dir)+ngrow
      else
       print *,"dimension bust"
       stop
      endif
#elif (AMREX_SPACEDIM==2)
      if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif
#else
      print *,"dimension bust"
      stop
#endif
      return
      end subroutine set_dimdec

      subroutine RT_transformVEL(x,vel,velT)
      use probcommon_module

      IMPLICIT NONE

      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: vel(SDIM)
      REAL_T, intent(out) :: velT(SDIM)
      INTEGER_T dir
     
      do dir=1,SDIM
       velT(dir)=vel(dir)
      enddo
      if ((levelrz.eq.0).or.(levelrz.eq.1)) then
       ! do nothing
      else if (levelrz.eq.3) then
       velT(1)=vel(1)*cos(x(2))+vel(2)*sin(x(2))
       velT(2)=-vel(1)*sin(x(2))+vel(2)*cos(x(2))
      else
       print *,"levelrz invalid rt transform vel"
       stop
      endif

      return
      end subroutine RT_transformVEL

 
      subroutine MinModGridInterp(data,xpos,x,H)
      IMPLICIT NONE

      REAL_T data(D_DECL(3,3,3))
      REAL_T xpos(D_DECL(3,3,3),SDIM)
      REAL_T x(SDIM)
      REAL_T x0(SDIM)
      REAL_T H
      REAL_T intercept
      REAL_T slopes(SDIM)
      INTEGER_T dir,i1,j1,k1,k1lo,k1hi
      REAL_T mindata,maxdata

      if (SDIM.eq.3) then
       k1lo=1
       k1hi=3
      else if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else
       print *,"dimension bust"
       stop
      endif

      do dir=1,SDIM
       x0(dir)=xpos(D_DECL(2,2,2),dir)
      enddo
      call minmod_stencil(data,xpos,slopes,intercept)
      H=intercept
      do dir=1,SDIM
       H=H+slopes(dir)*(x(dir)-x0(dir))
      enddo

      mindata=data(D_DECL(2,2,2))
      maxdata=data(D_DECL(2,2,2))
      do i1=1,3
      do j1=1,3
      do k1=k1lo,k1hi
       if (data(D_DECL(i1,j1,k1)).lt.mindata) then
        mindata=data(D_DECL(i1,j1,k1))
       endif
       if (data(D_DECL(i1,j1,k1)).gt.maxdata) then
        maxdata=data(D_DECL(i1,j1,k1))
       endif
      enddo 
      enddo 
      enddo 
      if (mindata.gt.H) then
       H=mindata
      endif
      if (maxdata.lt.H) then
       H=maxdata
      endif
 
      return 
      end subroutine MinModGridInterp


      subroutine cramers_rule(AA,XX,BB,matstatus)
      IMPLICIT NONE

      REAL_T AA(2,2)
      REAL_T XX(2)
      REAL_T BB(2)
      REAL_T det
      INTEGER_T matstatus

      matstatus=0
      XX(1)=zero
      XX(2)=zero

      det=AA(1,1)*AA(2,2)-AA(1,2)*AA(2,1)
      if (det.ne.zero) then
       matstatus=1
       XX(1)=(AA(2,2)*BB(1)-BB(2)*AA(1,2))/det
       XX(2)=(AA(1,1)*BB(2)-BB(1)*AA(2,1))/det
      endif

      return
      end subroutine cramers_rule


      subroutine get_iten(im1,im2,iten,nmat)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: im1,im2,nmat
      INTEGER_T, intent(out) :: iten
      INTEGER_T im,im_opp
      INTEGER_T im_iter
      INTEGER_T previous_count

      if ((im1.lt.1).or.(im1.gt.nmat).or. & 
          (im2.lt.1).or.(im2.gt.nmat).or. &
          (im1.eq.im2)) then
       print *,"im1,im2 mismatch im1,im2=",im1,im2
       print *,"num_materials=",nmat
       stop
      endif

      if (im1.lt.im2) then
       im=im1
       im_opp=im2
      else if (im2.lt.im1) then
       im=im2
       im_opp=im1
      else
       print *,"im1,im2 mismatch2 im1,im2=",im1,im2
       print *,"num_materials=",nmat
       stop
      endif
      if (im_opp.gt.im) then
       previous_count=0
       do im_iter=1,im-1
        previous_count=previous_count+nmat-im_iter
       enddo
       iten=previous_count+im_opp-im
      else
       print *,"im or im_opp invalid"
       stop
      endif

      return  
      end subroutine get_iten    


      subroutine get_inverse_iten(im1,im2,iten,nmat)
      IMPLICIT NONE

      INTEGER_T im1,im2,iten,nmat
      INTEGER_T im,im_opp,iten_test

      im1=0
      im2=0
      do im=1,nmat
      do im_opp=im+1,nmat
       call get_iten(im,im_opp,iten_test,nmat)
       if (iten.eq.iten_test) then
        im1=im
        im2=im_opp
       endif
      enddo
      enddo

      if ((im1.lt.1).or.(im1.gt.nmat).or. & 
          (im2.lt.1).or.(im2.gt.nmat).or. &
          (im1.eq.im2)) then
       print *,"im1,im2 mismatch im1,im2=",im1,im2
       print *,"num_materials=",nmat
       stop
      endif

      return  
      end subroutine get_inverse_iten    

        subroutine check_inbox(xx,xsten,nhalf,inboxflag)
        use probcommon_module

        IMPLICIT NONE

        INTEGER_T nhalf
        REAL_T dx
        REAL_T xsten(-nhalf:nhalf,SDIM)
        REAL_T xx(SDIM)
        INTEGER_T inboxflag,dir

        if (nhalf.lt.1) then
         print *,"nhalf invalid check inbox"
         stop
        endif
        inboxflag=1
        do dir=1,SDIM 
         dx=xsten(1,dir)-xsten(-1,dir)
         if ((xx(dir).lt.xsten(-1,dir)-VOFTOL*dx).or. &
             (xx(dir).gt.xsten(1,dir)+VOFTOL*dx)) then
          inboxflag=0
         endif
        enddo
 
        return
        end subroutine check_inbox


      subroutine planearea(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,perim)
      IMPLICIT NONE

      REAL_T x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL_T perim,perim1,perim2

      perim1=((y2-y1)*(z3-z1)-(y3-y1)*(z2-z1))**2
      perim1=perim1+((x2-x1)*(z3-z1)-(z2-z1)*(x3-x1))**2
      perim1=perim1+((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))**2

      perim2=((y2-y4)*(z3-z4)-(y3-y4)*(z2-z4))**2
      perim2=perim2+((x2-x4)*(z3-z4)-(z2-z4)*(x3-x4))**2
      perim2=perim2+((x2-x4)*(y3-y4)-(x3-x4)*(y2-y4))**2

      perim=half*(sqrt(perim1)+sqrt(perim2))

      return
      end subroutine planearea

 

      REAL_T function getcutvol(r,L)
      IMPLICIT NONE

      REAL_T r,L
 
      getcutvol=Pi*r*r*(r+sqrt(r*r-L*L))- &
        (Pi/three)*(r**3+(r*r-L*L)**(3.0/2.0))

      return
      end function getcutvol

      subroutine getsphere(L,V,angle)
      IMPLICIT NONE

      REAL_T L,V,angle
      REAL_T r1,r2,v1,v2,rmid,vmid
      INTEGER_T i

      r1=L
      r2=ten
      v1=getcutvol(r1,L) 
      v2=getcutvol(r2,L) 
      if (v1.gt.V) then
       rmid=L
       vmid=V
       angle=90.0
      else
       do i=1,100
        rmid=(r1+r2)*half
        vmid=getcutvol(rmid,L)
        if ((V.ge.v1).and.(V.le.vmid)) then
         r2=rmid
         v2=vmid
        else if ((V.le.v2).and.(V.ge.vmid)) then
         r1=rmid
         v1=vmid
        else
         print *,"cannot find angle r1,r2,v1,v2,rmid,vmid,L,V ", &
          r1,r2,v1,v2,rmid,vmid,L,V
         stop
        endif
       enddo   
       angle=90.0+atan(sqrt(rmid**2-L**2)/L)*180.0/Pi
      endif

      return
      end subroutine getsphere

      REAL_T function hsprime(phi,cutoff)
      IMPLICIT NONE
      REAL_T phi,cutoff

      if ((phi.ge.cutoff).or.(phi.le.-cutoff)) then
       hsprime=0
      else
       hsprime=(1.0/(2.0*cutoff))*(1.0+cos(Pi*phi/cutoff))
      endif
      return
      end function hsprime

      INTEGER_T function sign_funct(LS)
      IMPLICIT NONE

      REAL_T LS

      if (LS.lt.zero) then
       sign_funct=-1
      else if (LS.ge.zero) then
       sign_funct=1
      else
       print *,"LS bust LS=",LS
       stop
      endif

      end function sign_funct

      REAL_T function hssign(phi,cutoff)
      IMPLICIT NONE
      REAL_T phi,cutoff

      hssign=two*hs(phi,cutoff)-one
      return
      end function hssign

      REAL_T function bellcurve(phi,cutoff)
      IMPLICIT NONE
      REAL_T phi,cutoff

      if (cutoff.le.zero) then
       print *,"cutoff must be positive"
       stop
      endif
      bellcurve=exp(-(phi/cutoff)**2)

      return
      end function bellcurve

      REAL_T function slopewt(phi,cutoff)
      IMPLICIT NONE
      REAL_T phi,cutoff

      if ((phi.ge.cutoff).or.(phi.le.-cutoff)) then
       slopewt=zero
      else
       slopewt=one+cos(Pi*phi/cutoff)
      endif

      return
      end function slopewt

      REAL_T function hs(phi,cutoff)
      IMPLICIT NONE
      REAL_T phi,cutoff,EPS

      EPS=1.0D-6
      if (phi.ge.cutoff) then
        hs=one
      else if (phi.le.-cutoff) then
        hs=zero
      else
        hs=phi/(two*cutoff)+sin(Pi*phi/cutoff)/(two*Pi)+half
        if ((hs.lt.-EPS).or.(hs.gt.one+EPS)) then
         print *,"hs sanity check failed"
         stop
        else if (hs.lt.zero) then
         hs=zero
        else if (hs.gt.one) then
         hs=one
        endif
      endif

      return
      end function hs

      REAL_T function hs_scale(phi,cutoff)
      IMPLICIT NONE
      REAL_T phi,cutoff,EPS

      EPS=1.0D-6
      if (phi.ge.cutoff) then
        hs_scale=one
      else if (phi.le.-cutoff) then
        hs_scale=zero
      else
        hs_scale=half*(half+(phi/cutoff)+half*((phi/cutoff)**2)- &
         (cos(two*Pi*phi/cutoff)-one)/(four*Pi*Pi)+ &
         (phi+cutoff)*sin(Pi*phi/cutoff)/(cutoff*Pi))

        if ((hs_scale.lt.-EPS).or.(hs_scale.gt.one+EPS)) then
         print *,"hs_scale sanity check failed"
         stop
        else if (hs_scale.lt.zero) then
         hs_scale=zero
        else if (hs_scale.gt.one) then
         hs_scale=one
        endif
      endif

      return
      end function hs_scale

        REAL_T function approximate(start_x,end_x,start_y,end_y,my_x)
        IMPLICIT NONE

        REAL_T start_x,end_x,start_y,end_y,my_x

        if (abs(start_x-end_x).lt.1.0D-15) then
         approximate=(end_y+start_y)/2.0
        else
         approximate=start_y+(end_y-start_y)*(my_x-start_x)/ &
           (end_x-start_x)
        endif
        return
        end function approximate


      subroutine consistent_materials(vfrac,cen,nmat)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T nmat
      REAL_T vfrac(nmat)
      REAL_T cen(SDIM,nmat)

      INTEGER_T imaterial,dir
      REAL_T voftotal


      if (nmat.ne.num_materials) then
       print *,"nmat bust"
       stop
      endif

      voftotal=zero

      do imaterial=1,nmat
       if (vfrac(imaterial).le.VOFTOL) then
        vfrac(imaterial)=zero
        do dir=1,SDIM
         cen(dir,imaterial)=zero
        enddo
       endif
       if (vfrac(imaterial).ge.one-VOFTOL) then
        vfrac(imaterial)=one
        do dir=1,SDIM
         cen(dir,imaterial)=zero
        enddo
       endif
       if (is_rigid(nmat,imaterial).eq.0) then
        voftotal=voftotal+vfrac(imaterial)
       else if (is_rigid(nmat,imaterial).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif
      enddo ! imaterial

      if (voftotal.le.zero) then
       print *,"vacuum bust in consistent materials"
       print *,"voftotal= ",voftotal
       stop
      endif

      do imaterial=1,nmat
       vfrac(imaterial)=vfrac(imaterial)/voftotal
      enddo 

      return
      end subroutine consistent_materials

      subroutine boxrefine(clo,chi,flo,fhi,dir)
      IMPLICIT NONE

      INTEGER_T clo(SDIM),chi(SDIM)
      INTEGER_T flo(SDIM),fhi(SDIM)
      INTEGER_T boxdir,dir

      if ((dir.lt.-1).or.(dir.gt.SDIM)) then
       print *,"dir invalid boxrefine"
       stop
      endif
      do boxdir=1,SDIM
       if ((dir+1.eq.boxdir).or.(dir.eq.SDIM)) then
        flo(boxdir)=2*clo(boxdir)
        fhi(boxdir)=2*chi(boxdir)+1
       else 
        flo(boxdir)=clo(boxdir)
        fhi(boxdir)=chi(boxdir)
       endif
      enddo
      do boxdir=1,SDIM
       if (flo(boxdir).gt.fhi(boxdir)) then
        print *,"flo,fhi sanity check failed"
        stop
       endif
      enddo

      return
      end subroutine boxrefine

      subroutine dim_to_box( &
        DIMS(data), &
        datalo,datahi)
      IMPLICIT NONE

      INTEGER_T DIMDEC(data)
      INTEGER_T datalo(SDIM),datahi(SDIM)
      INTEGER_T dir

      datalo(1)=ARG_L1(data)
      datahi(1)=ARG_H1(data)
      datalo(2)=ARG_L2(data)
      datahi(2)=ARG_H2(data)
#if (AMREX_SPACEDIM==3)
      datalo(SDIM)=ARG_L3(data)
      datahi(SDIM)=ARG_H3(data)
#endif
      do dir=1,SDIM
       if (datalo(dir).gt.datahi(dir)) then
        print *,"datalo, datahi sanity check failed"
        stop
       endif
      enddo

      return
      end subroutine dim_to_box


      subroutine box_to_dim( &
        DIMS(data), &
        datalo,datahi)
      IMPLICIT NONE

      INTEGER_T DIMDEC(data)
      INTEGER_T datalo(SDIM),datahi(SDIM)
      INTEGER_T dir

      ARG_L1(data)=datalo(1)
      ARG_H1(data)=datahi(1)
      ARG_L2(data)=datalo(2)
      ARG_H2(data)=datahi(2)
#if (AMREX_SPACEDIM==3)
      ARG_L3(data)=datalo(SDIM)
      ARG_H3(data)=datahi(SDIM)
#endif
      do dir=1,SDIM
       if (datalo(dir).gt.datahi(dir)) then
        print *,"datalo, datahi sanity check failed"
        stop
       endif
      enddo

      return
      end subroutine box_to_dim

      subroutine gridvol(xsten,nhalf,rzflag,vol)
      IMPLICIT NONE

      REAL_T vol
      INTEGER_T nhalf
      REAL_T xsten(-nhalf:nhalf,SDIM)
      INTEGER_T rzflag
      REAL_T RCENTER

      if (nhalf.lt.1) then
       print *,"nhalf invalid gridvol"
       stop
      endif

      if (rzflag.eq.0) then
       if (SDIM.eq.3) then
        vol=(xsten(1,1)-xsten(-1,1))*(xsten(1,2)-xsten(-1,2))* &
         (xsten(1,SDIM)-xsten(-1,SDIM))
       else if (SDIM.eq.2) then
        vol=(xsten(1,1)-xsten(-1,1))*(xsten(1,2)-xsten(-1,2))
       else
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       RCENTER=half*abs(xsten(-1,1)+xsten(1,1))
       vol=two*Pi*RCENTER* &
         (xsten(1,1)-xsten(-1,1))*(xsten(1,2)-xsten(-1,2))
      else if (rzflag.eq.3) then
       RCENTER=half*abs(xsten(-1,1)+xsten(1,1))
       if (SDIM.eq.3) then
        vol=RCENTER* &
         (xsten(1,1)-xsten(-1,1))*(xsten(1,2)-xsten(-1,2))* &
         (xsten(1,SDIM)-xsten(-1,SDIM))
       else if (SDIM.eq.2) then
        vol=RCENTER* &
         (xsten(1,1)-xsten(-1,1))*(xsten(1,2)-xsten(-1,2))
       else
        print *,"dimension bust"
        stop
       endif
      else 
       print *,"rzflag invalid"
       stop
      endif

      return
      end subroutine gridvol



      subroutine gridarea(xsten,nhalf,rzflag,dir,side,area,areacen)
      IMPLICIT NONE

      REAL_T, intent(out) :: area
      REAL_T, intent(out) :: areacen(SDIM)
      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
      INTEGER_T, intent(in) :: rzflag
      INTEGER_T, intent(in) :: dir,side
      INTEGER_T iside,dir2
      REAL_T RR,R1,R2,RC

      if (side.eq.0) then
       iside=-1
      else if (side.eq.1) then
       iside=1
      else
       print *,"side invalid"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid gridarea"
       stop
      endif

      areacen(dir+1)=xsten(iside,dir+1)

      area=one
      do dir2=1,SDIM
       if (dir2.ne.dir+1) then
        area=area*(xsten(1,dir2)-xsten(-1,dir2))
        areacen(dir2)=half*(xsten(1,dir2)+xsten(-1,dir2))
        if (area.le.zero) then
         print *,"area bust"
         stop
        endif
       endif
      enddo  ! dir2

      if (rzflag.eq.0) then
       ! do nothing
      else if (rzflag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       if (dir.eq.0) then  ! R face
        RR=abs(xsten(iside,1))
       else if (dir.eq.1) then  ! Z face
        dir2=1
        RR=abs(xsten(0,1))
        R1=xsten(-1,1)
        R2=xsten(1,1)
        if (R2+R1.eq.zero) then
         print *,"R2+R1 invalid"
         stop
        endif
        RC=two*(R2*R2+R1*R1+R2*R1)/(three*(R2+R1))
        areacen(dir2)=RC
       else
        print *,"dir invalid gridarea 2"
        stop
       endif
       area=area*two*Pi*RR

      else if (rzflag.eq.3) then

       dir2=1
       R1=xsten(-1,1)
       R2=xsten(1,1)
       if (R2+R1.eq.zero) then
        print *,"R2+R1 invalid"
        stop
       endif
       RC=two*(R2*R2+R1*R1+R2*R1)/(three*(R2+R1))

       if (dir.eq.0) then ! R face
        RR=abs(xsten(iside,1))
       else if (dir.eq.1) then ! Theta face
        RR=one
        areacen(dir2)=RC
       else if ((dir.eq.2).and.(SDIM.eq.3)) then ! Z face
        RR=abs(xsten(0,1))
        areacen(dir2)=RC
       else
        print *,"dir invalid gridarea 3"
        stop
       endif
       area=area*RR

      else 
       print *,"rzflag invalid"
       stop
      endif

      return
      end subroutine gridarea



      subroutine element_nodes1D(xnodes,xlo,ielem,lo_e,dx,bfact)
      use LegendreNodes

      IMPLICIT NONE

      INTEGER_T bfact
      REAL_T xnodes(bfact)
      REAL_T xlo
      INTEGER_T ielem
      INTEGER_T lo_e,inodes
      REAL_T dx
  
      if (bfact.eq.1) then
       do inodes=1,bfact
        xnodes(inodes)=xlo+(ielem-lo_e)*bfact*dx+(inodes-half)*dx 
       enddo
      else if (bfact.gt.1) then
       do inodes=1,bfact
        xnodes(inodes)=xlo+(ielem-lo_e)*bfact*dx+ &
         bfact*dx*(cache_gauss(bfact,inodes-1,SPTYPE)+one)/two
       enddo
      else
       print *,"bfact invalid20"
       stop
      endif

      end subroutine element_nodes1D


       ! returns the Gauss Lobatto points in element "ielem"
      subroutine element_GLnodes1D(xnodes,xlo,ielem,lo_e,dx,bfact)
      use LegendreNodes

      IMPLICIT NONE

      INTEGER_T bfact
      REAL_T xnodes(bfact+1)
      REAL_T xlo
      INTEGER_T ielem
      INTEGER_T lo_e,inodes
      REAL_T dx
  
      if (bfact.eq.1) then  ! evenly spaced points
       do inodes=1,bfact+1
        xnodes(inodes)=xlo+(ielem-lo_e)*bfact*dx+(inodes-one)*dx 
       enddo
      else if (bfact.gt.1) then  ! gauss lobatto points
       do inodes=1,bfact+1
        xnodes(inodes)=xlo+(ielem-lo_e)*bfact*dx+ &
         bfact*dx*(cache_gauss_lobatto(bfact,inodes-1,SPTYPE)+one)/two
       enddo
      else
       print *,"bfact invalid21"
       stop
      endif

      end subroutine element_GLnodes1D


      subroutine gridsten(x,xlo,i,j,k,fablo,bfact,dx,nhalf)
      IMPLICIT NONE 

      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(out) :: x(-nhalf:nhalf,SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM)
      INTEGER_T, intent(in) :: i,j,k,bfact
      INTEGER_T dir,icrit
      REAL_T, dimension(:), allocatable :: xsub

      if (bfact.lt.1) then
       print *,"bfact invalid22"
       stop
      endif

      do dir=1,SDIM
       allocate(xsub(-nhalf:nhalf))
       if (dir.eq.1) then
        icrit=i
       else if (dir.eq.2) then
        icrit=j
       else if ((dir.eq.3).and.(SDIM.eq.3)) then
        icrit=k
       else
        print *,"dir invalid gridsten"
        stop
       endif
       call gridsten1D(xsub,xlo,icrit,fablo,bfact,dx,dir,nhalf)
       do icrit=-nhalf,nhalf
        x(icrit,dir)=xsub(icrit)
       enddo

       deallocate(xsub)
      enddo ! dir

      end subroutine gridsten

        ! -1<=grid_type<=5
      subroutine gridstenMAC(x,xlo,i,j,k,fablo,bfact,dx,nhalf,grid_type, &
                      caller_id)
      IMPLICIT NONE 

      INTEGER_T, intent(in) :: caller_id
      INTEGER_T, intent(in) :: nhalf
      INTEGER_T, intent(in) :: grid_type
      REAL_T, intent(out) :: x(-nhalf:nhalf,SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM)
      INTEGER_T, intent(in) :: i,j,k,bfact
      INTEGER_T dir,icrit
      REAL_T, dimension(:), allocatable :: xsub
      INTEGER_T :: box_type(SDIM)

      if ((grid_type.ge.-1).and.(grid_type.le.5)) then
       ! do nothing
      else
       print *,"grid_type invalid gridstenMAC: grid_type,caller_id", &
          grid_type,caller_id
       stop
      endif

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type(grid_type,box_type)

      if (bfact.lt.1) then
       print *,"bfact invalid23"
       stop
      endif

      do dir=1,SDIM
       allocate(xsub(-nhalf:nhalf))
       if (dir.eq.1) then
        icrit=i
       else if (dir.eq.2) then
        icrit=j
       else if ((dir.eq.3).and.(SDIM.eq.3)) then
        icrit=k
       else
        print *,"dir invalid gridsten mac"
        stop
       endif 
        ! NODE
       if (box_type(dir).eq.1) then
        call gridsten1DMAC(xsub,xlo,icrit,fablo,bfact,dx,dir,nhalf)
       else if (box_type(dir).eq.0) then ! CELL
        call gridsten1D(xsub,xlo,icrit,fablo,bfact,dx,dir,nhalf)
       else
        print *,"box_type(dir) invalid"
        stop
       endif
       do icrit=-nhalf,nhalf
        x(icrit,dir)=xsub(icrit)
       enddo

       deallocate(xsub)
      enddo ! dir=1..sdim

      end subroutine gridstenMAC



      subroutine gridstenND(x,xlo,i,j,k,fablo,bfact,dx,nhalf)
      IMPLICIT NONE 

      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(out) :: x(-nhalf:nhalf,SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM)
      INTEGER_T, intent(in) :: i,j,k,bfact
      INTEGER_T dir,icrit
      REAL_T, dimension(:), allocatable :: xsub

      if (bfact.lt.1) then
       print *,"bfact invalid24"
       stop
      endif

      do dir=1,SDIM
       allocate(xsub(-nhalf:nhalf))
       if (dir.eq.1) then
        icrit=i
       else if (dir.eq.2) then
        icrit=j
       else if ((dir.eq.3).and.(SDIM.eq.3)) then
        icrit=k
       else
        print *,"dir invalid gridsten nd"
        stop
       endif
       call gridsten1DMAC(xsub,xlo,icrit,fablo,bfact,dx,dir,nhalf)
       do icrit=-nhalf,nhalf
        x(icrit,dir)=xsub(icrit)
       enddo

       deallocate(xsub)
      enddo ! dir

      end subroutine gridstenND

      subroutine get_element_index(i,ielem,isub,bfact)
      IMPLICIT NONE

      INTEGER_T i,ielem,isub,bfact,ipos

      if (bfact.eq.1) then
       ielem=i
       isub=0
      else if (bfact.gt.1) then
       if (i.ge.0) then
        ielem=i/bfact
        isub=i-bfact*ielem
       else if (i.lt.0) then
          ! e.g. bfact=2 i=-1 ipos=0 ielem=-1 isub=1
          ! e.g. bfact=2 i=-2 ipos=1 ielem=-1 isub=0
          ! e.g. bfact=2 i=-3 ipos=2 ielem=-2 isub=1
          ! e.g. bfact=2 i=-4 ipos=3 ielem=-2 isub=0 ....
        ipos=-(i+1)
        ielem=ipos/bfact
        isub=ipos-bfact*ielem
        ielem=-ielem-1
        isub=bfact-isub-1
       else
        print *,"i invalid"
        stop
       endif
      else
       print *,"bfact invalid25"
       stop
      endif

      end subroutine get_element_index

       ! cache_gauss(i,j)    j=0..i-1
       ! cache_gauss_lobatto(i,j)  j=0..i
       ! it is assumed that i=0 corresponds to the first cell in the domain.
       ! x(0) is a gauss node
       ! x(-1) and x(1) are gauss-Lobatto nodes....
       ! 1<=dir<=sdim
      subroutine gridsten1D(x,xlo,i,fablo,bfact,dx,dir,nhalf)
      use LegendreNodes

      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(out) :: x(-nhalf:nhalf)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM)
      REAL_T xdomlo,xofs,xabs,xnode
      INTEGER_T, intent(in) :: i,bfact,dir
      INTEGER_T side
      INTEGER_T icell,imac,isten,signflag,ishift,i_e,i_n

      if (bfact.lt.1) then
       print *,"bfact invalid26"
       stop
      endif
      if (nhalf.lt.0) then
       print *,"nhalf invalid gridsten1d"
       stop
      endif
      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid gridsten 1d"
       stop
      endif

      xdomlo=xlo(dir)-fablo(dir)*dx(dir)

      do side=-1,1,2

       do icell=0,nhalf/2

        isten=i+side*icell
        if (isten.lt.0) then
         signflag=-1
         ishift=-1-isten
        else
         signflag=1
         ishift=isten
        endif

        if (bfact.eq.1) then
         xofs=(ishift+half)*dx(dir)
         xabs=xdomlo+signflag*xofs
        else if (bfact.gt.1) then
         i_e=ishift/bfact
         i_n=ishift-i_e*bfact
         if ((i_n.ge.0).and.(i_n.lt.bfact)) then
          xnode=cache_gauss(bfact,i_n,SPTYPE)
          xofs=i_e*bfact*dx(dir)+(xnode+one)*bfact*dx(dir)/two
          xabs=xdomlo+signflag*xofs
         else
          print *,"i_n invalid"
          stop
         endif
        else
         print *,"bfact invalid27"
         stop
        endif
        x(side*2*icell)=xabs

       enddo ! icell

       do imac=1,(nhalf+1)/2

        if (side.eq.1) then
         isten=i+side*imac
        else if (side.eq.-1) then
         isten=i+side*imac+1
        else
         print *,"side invalid"
         stop
        endif

        if (isten.lt.0) then
         signflag=-1
         ishift=-isten
        else 
         signflag=1
         ishift=isten
        endif

        if (bfact.eq.1) then
         xofs=ishift*dx(dir)
         xabs=xdomlo+signflag*xofs
        else if (bfact.gt.1) then
         i_e=ishift/bfact
         i_n=ishift-i_e*bfact
         if ((i_n.ge.0).and.(i_n.lt.bfact)) then
          xnode=cache_gauss_lobatto(bfact,i_n,SPTYPE)
          xofs=i_e*bfact*dx(dir)+(xnode+one)*bfact*dx(dir)/2.0
          xabs=xdomlo+signflag*xofs
         else
          print *,"i_n invalid"
          stop
         endif
        else
         print *,"bfact invalid28"
         stop
        endif
    
        x(side*(2*imac-1))=xabs
       enddo ! imac
      enddo ! side

      return
      end subroutine gridsten1D

       ! x(0) is a gauss-Lobatto node
       ! x(-1) and x(1) are gauss nodes....
      subroutine gridsten1DMAC(x,xlo,i,fablo,bfact,dx,dir,nhalf)
      use LegendreNodes

      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(out) :: x(-nhalf:nhalf)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM)
      REAL_T xdomlo,xofs,xabs,xnode
      INTEGER_T, intent(in) :: i,bfact,dir
      INTEGER_T side
      INTEGER_T icell,imac,isten,signflag,ishift,i_e,i_n

      if (bfact.lt.1) then
       print *,"bfact invalid29"
       stop
      endif
      if (nhalf.lt.0) then
       print *,"nhalf invalid gridsten 1d mac"
       stop
      endif
      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid gridsten1d mac"
       stop
      endif

      xdomlo=xlo(dir)-fablo(dir)*dx(dir)

      do side=-1,1,2

       do icell=1,(nhalf+1)/2

         ! i is a MAC index; isten is a cell index
        if (side.eq.1) then
         isten=i+side*icell-1
        else if (side.eq.-1) then
         isten=i+side*icell
        else
         print *,"side invalid"
         stop
        endif

        if (isten.lt.0) then
         signflag=-1
         ishift=-1-isten
        else
         signflag=1
         ishift=isten
        endif

        if (bfact.eq.1) then
         xofs=(ishift+half)*dx(dir)
         xabs=xdomlo+signflag*xofs
        else if (bfact.gt.1) then
         i_e=ishift/bfact
         i_n=ishift-i_e*bfact
         if ((i_n.ge.0).and.(i_n.lt.bfact)) then
          xnode=cache_gauss(bfact,i_n,SPTYPE)
          xofs=i_e*bfact*dx(dir)+(xnode+one)*bfact*dx(dir)/2.0
          xabs=xdomlo+signflag*xofs
         else
          print *,"i_n invalid"
          stop
         endif
        else
         print *,"bfact invalid30"
         stop
        endif
        x(side*(2*icell-1))=xabs

       enddo ! icell

       do imac=0,nhalf/2

         ! isten and i are bot MAC indexes.
        isten=i+side*imac

        if (isten.lt.0) then
         signflag=-1
         ishift=-isten
        else 
         signflag=1
         ishift=isten
        endif

        if (bfact.eq.1) then
         xofs=ishift*dx(dir)
         xabs=xdomlo+signflag*xofs
        else if (bfact.gt.1) then
         i_e=ishift/bfact
         i_n=ishift-i_e*bfact
         if ((i_n.ge.0).and.(i_n.lt.bfact)) then
          xnode=cache_gauss_lobatto(bfact,i_n,SPTYPE)
          xofs=i_e*bfact*dx(dir)+(xnode+one)*bfact*dx(dir)/2.0
          xabs=xdomlo+signflag*xofs
         else
          print *,"i_n invalid"
          stop
         endif
        else
         print *,"bfact invalid31"
         stop
        endif
    
        x(side*2*imac)=xabs
       enddo ! imac
      enddo ! side

      return
      end subroutine gridsten1DMAC





      subroutine get_istar(icomp,istar)
      IMPLICIT NONE

      INTEGER_T icomp,nstar,icomplocal,dir
      INTEGER_T istar(3)

      nstar=9
      if (SDIM.eq.3) then
       nstar=nstar*3
      endif
      if ((icomp.lt.1).or.(icomp.gt.nstar)) then
       print *,"icomp invalid"
       stop
      endif

       ! 3,3,3
       ! icomp=(istar-1)+(jstar-1)*3+(kstar-1)*9+1 
      istar(3)=0
      icomplocal=icomp-1
      if (SDIM.eq.3) then
       istar(3)=icomplocal/9
       icomplocal=icomplocal-9*istar(3)
      endif
      istar(2)=icomplocal/3
      icomplocal=icomplocal-3*istar(2)
      istar(1)=icomplocal
      icomplocal=icomplocal-istar(1)
      do dir=1,SDIM
       istar(dir)=istar(dir)-1
       if ((istar(dir).lt.-1).or.(istar(dir).gt.1)) then
        print *,"istar invalid"
        stop
       endif
      enddo

      return
      end subroutine get_istar 


      subroutine put_istar(icomp,istar)
      IMPLICIT NONE

      INTEGER_T icomp,nstar,dir
      INTEGER_T istar(3)

      nstar=9
      if (SDIM.eq.3) then
       nstar=nstar*3
      endif
      do dir=1,3
       if ((istar(dir).lt.-1).or.(istar(dir).gt.1)) then
        print *,"istar invalid"
        stop
       endif
      enddo

       ! 3,3,3
       ! icomp=(istar-1)+(jstar-1)*3+(kstar-1)*9+1 
      icomp=0
      if (SDIM.eq.3) then
       icomp=icomp+9*(istar(3)+1)
      endif
      icomp=icomp+3*(istar(2)+1)
      icomp=icomp+(istar(1)+1)
      icomp=icomp+1
      if ((icomp.lt.1).or.(icomp.gt.nstar)) then
       print *,"icomp invalid"
       stop
      endif

      return
      end subroutine put_istar 

      subroutine growntilebox( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,ng)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(out) :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: ng
      INTEGER_T dir2

      growlo(3)=0
      growhi(3)=0

      do dir2=1,SDIM
       growlo(dir2)=tilelo(dir2)
       if (tilelo(dir2).eq.fablo(dir2)) then
        growlo(dir2)=tilelo(dir2)-ng
       else if (tilelo(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilelo(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
       growhi(dir2)=tilehi(dir2)
       if (tilehi(dir2).eq.fabhi(dir2)) then
        growhi(dir2)=tilehi(dir2)+ng
       else if (tilehi(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilehi(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
      enddo ! dir2
      
      return
      end subroutine growntilebox


      subroutine growntilebox_TILE( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,ng)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(out) :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: ng
      INTEGER_T dir2

      growlo(3)=0
      growhi(3)=0

      do dir2=1,SDIM
       growlo(dir2)=tilelo(dir2)-ng
       if (tilelo(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilelo(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
       growhi(dir2)=tilehi(dir2)+ng
       if (tilehi(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilehi(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
      enddo ! dir2
      
      return
      end subroutine growntilebox_TILE


      subroutine intersect_weight_avg(ic,i,bfact_c,bfact_f,wt)
      use probcommon_module

      INTEGER_T ic,i,bfact_c,bfact_f
      REAL_T wt
      INTEGER_T nhalf,dir_index
      INTEGER_T fablo(SDIM)
      REAL_T intlo,inthi
      REAL_T dxc(SDIM)
      REAL_T dxf(SDIM)
      REAL_T xlo(SDIM)
      REAL_T xc(-1:1)
      REAL_T xf(-1:1)

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      nhalf=1

      dxc(1)=one
      dxf(1)=half

      xlo(1)=zero
      dir_index=1
      fablo(1)=0
      call gridsten1D(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
      call gridsten1D(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

      if (bfact_f.eq.1) then
       if (abs(xf(1)-xf(-1)-dxf(1)).ge.VOFTOL) then
        print *,"xf invalid 1"
        stop
       endif
       if (abs(xf(1)+xf(-1)-two*xf(0)).ge.VOFTOL) then
        print *,"xf invalid 2"
        stop
       endif
      else if (bfact_f.gt.1) then
       if ((xf(1)-xf(0).le.zero).or. &
           (xf(0)-xf(-1).le.zero)) then
        print *,"xf invalid 3"
        stop
       endif
      else
       print *,"bfact_f invalid"
       stop
      endif

      if (bfact_c.eq.1) then
       if (abs(xc(1)-xc(-1)-dxc(1)).ge.VOFTOL) then
        print *,"xc invalid"
        stop
       endif
       if (abs(xc(1)+xc(-1)-two*xc(0)).ge.VOFTOL) then
        print *,"xc invalid"
        stop
       endif
      else if (bfact_c.gt.1) then
       if ((xc(1)-xc(0).le.zero).or. &
           (xc(0)-xc(-1).le.zero)) then
        print *,"xc invalid 3"
        stop
       endif
      else
       print *,"bfact_c invalid"
       stop
      endif

      intlo=max(xf(-1),xc(-1))
      inthi=min(xf(1),xc(1))
      if (intlo.ge.inthi) then
       wt=zero
      else
       wt=(inthi-intlo)/(xc(1)-xc(-1)) ! average down
      endif
      if (abs(wt).le.VOFTOL) then
       wt=zero
      else if (abs(wt-one).le.VOFTOL) then
       wt=one
      else if ((wt.gt.zero).and.(wt.lt.one)) then
       ! do nothing
      else
       print *,"wt invalid 1"
       stop
      endif

      return
      end subroutine intersect_weight_avg


      subroutine intersect_weight_avg_COPY(ic,i,bfact_c,bfact_f,wt, &
         dir,dir_current)
      use probcommon_module

      INTEGER_T ic,i,bfact_c,bfact_f,dir,dir_current
      REAL_T wt

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      if ((dir.ge.0).and.(dir.lt.SDIM).and. &
          (dir_current.ge.0).and.(dir_current.lt.SDIM)) then
       if (dir.eq.dir_current) then
        wt=one
       else if (dir.ne.dir_current) then
        call intersect_weight_avg(ic,i,bfact_c,bfact_f,wt)
       else
        print *,"dir or dir_current bust"
        stop
       endif
      else
       print *,"dir or dir_current bust"
       stop
      endif

      return
      end subroutine intersect_weight_avg_COPY


      subroutine intersect_weight_interp(ic,i,bfact_c,bfact_f,wt)
      use probcommon_module

      INTEGER_T ic,i,bfact_c,bfact_f
      REAL_T wt
      INTEGER_T nhalf,dir_index
      INTEGER_T fablo(SDIM)
      REAL_T intlo,inthi
      REAL_T dxc(SDIM)
      REAL_T dxf(SDIM)
      REAL_T xlo(SDIM)
      REAL_T xc(-1:1)
      REAL_T xf(-1:1)

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      nhalf=1

      dxc(1)=one
      dxf(1)=half

      xlo(1)=zero
      dir_index=1
      fablo(1)=0
      call gridsten1D(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
      call gridsten1D(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

      if (bfact_f.eq.1) then
       if (abs(xf(1)-xf(-1)-dxf(1)).ge.VOFTOL) then
        print *,"xf invalid 1"
        stop
       endif
       if (abs(xf(1)+xf(-1)-two*xf(0)).ge.VOFTOL) then
        print *,"xf invalid 2"
        stop
       endif
      else if (bfact_f.gt.1) then
       if ((xf(1)-xf(0).le.zero).or. &
           (xf(0)-xf(-1).le.zero)) then
        print *,"xf invalid 3"
        stop
       endif
      else
       print *,"bfact_f invalid"
       stop
      endif

      if (bfact_c.eq.1) then
       if (abs(xc(1)-xc(-1)-dxc(1)).ge.VOFTOL) then
        print *,"xc invalid"
        stop
       endif
       if (abs(xc(1)+xc(-1)-two*xc(0)).ge.VOFTOL) then
        print *,"xc invalid"
        stop
       endif
      else if (bfact_c.gt.1) then
       if ((xc(1)-xc(0).le.zero).or. &
           (xc(0)-xc(-1).le.zero)) then
        print *,"xc invalid 3"
        stop
       endif
      else
       print *,"bfact_c invalid"
       stop
      endif

      intlo=max(xf(-1),xc(-1))
      inthi=min(xf(1),xc(1))
      if (intlo.ge.inthi) then
       wt=zero
      else
       wt=(inthi-intlo)/(xf(1)-xf(-1)) ! interpolate
      endif
      if (abs(wt).le.VOFTOL) then
       wt=zero
      else if (abs(wt-one).le.VOFTOL) then
       wt=one
      else if ((wt.gt.zero).and.(wt.lt.one)) then
       if ((bfact_c.eq.1).and.(bfact_f.eq.1)) then
        print *,"coarse cell should completely cover fine cell"
        stop
       endif 
       if ((abs(intlo-xf(-1)).lt.VOFTOL*dxf(1)).or. &
           (abs(inthi-xf(1)).lt.VOFTOL*dxf(1))) then
        ! do nothing
       else
        print *,"coarse cell cannot be contained within fine cell"
        stop
       endif  
      else
       print *,"wt invalid 2"
       stop
      endif

      return
      end subroutine intersect_weight_interp


      subroutine intersect_weight_interp_COPY(ic,i,bfact_c,bfact_f,wt, &
         dir,dir_current)
      use probcommon_module

      INTEGER_T ic,i,bfact_c,bfact_f,dir,dir_current
      REAL_T wt

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      if ((dir.ge.0).and.(dir.lt.SDIM).and. &
          (dir_current.ge.0).and.(dir_current.lt.SDIM)) then
       if (dir.eq.dir_current) then
        wt=one
       else if (dir.ne.dir_current) then
        call intersect_weight_interp(ic,i,bfact_c,bfact_f,wt)
       else
        print *,"dir or dir_current bust"
        stop
       endif
      else
       print *,"dir or dir_current bust"
       stop
      endif

      return
      end subroutine intersect_weight_interp_COPY

      subroutine intersect_weightMAC_interp(ic,i,bfact_c,bfact_f,wt)
      use probcommon_module

      INTEGER_T ic,i,bfact_c,bfact_f
      REAL_T wt
      INTEGER_T nhalf,dir_index
      INTEGER_T fablo(SDIM)
      REAL_T intlo,inthi
      REAL_T dxc(SDIM)
      REAL_T dxf(SDIM)
      REAL_T xlo(SDIM)
      REAL_T xc(-1:1)
      REAL_T xf(-1:1)
      INTEGER_T denom

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      if (bfact_c.eq.1) then  ! intersect_weight_interpMAC
       if ((i/2)*2.eq.i) then
        if (2*ic.eq.i) then
         wt=one
        else
         wt=zero
        endif
       else if (((i/2)*2.eq.i-1).or.((i/2)*2.eq.i+1)) then
        if (2*ic.eq.i-1) then
         wt=half
        else if (2*ic.eq.i+1) then
         wt=half
        else if (2*ic.eq.i) then
         print *,"ic or i bust"
         stop
        else
         wt=zero
        endif
       else
        print *,"i/2 failed"
        stop
       endif
      else if (bfact_c.gt.1) then
       denom=2*bfact_c
       if (denom*(i/denom).eq.i) then
        if (2*ic.eq.i) then
         wt=one
        else
         wt=zero
        endif
       else

        nhalf=1

        dxc(1)=one
        dxf(1)=half

        xlo(1)=zero
        dir_index=1
        fablo(1)=0
        call gridsten1DMAC(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
        call gridsten1DMAC(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

        if (bfact_f.eq.1) then
         if (abs(xf(1)-xf(-1)-dxf(1)).ge.VOFTOL) then
          print *,"xf invalid 1"
          stop
         endif
         if (abs(xf(1)+xf(-1)-two*xf(0)).ge.VOFTOL) then
          print *,"xf invalid 2"
          stop
         endif
        else if (bfact_f.gt.1) then
         if ((xf(1)-xf(0).le.zero).or. &
             (xf(0)-xf(-1).le.zero)) then
          print *,"xf invalid 3"
          stop
         endif
        else
         print *,"bfact_f invalid"
         stop
        endif

        if (bfact_c.eq.1) then
         if (abs(xc(1)-xc(-1)-dxc(1)).ge.VOFTOL) then
          print *,"xc invalid"
          stop
         endif
         if (abs(xc(1)+xc(-1)-two*xc(0)).ge.VOFTOL) then
          print *,"xc invalid"
          stop
         endif
        else if (bfact_c.gt.1) then
         if ((xc(1)-xc(0).le.zero).or. &
             (xc(0)-xc(-1).le.zero)) then
          print *,"xc invalid 3"
          stop
         endif
        else
         print *,"bfact_c invalid"
         stop
        endif

        intlo=max(xf(-1),xc(-1))
        inthi=min(xf(1),xc(1))
        if (intlo.ge.inthi) then
         wt=zero
        else
         wt=(inthi-intlo)/(xf(1)-xf(-1))  ! interpolate
        endif
        if (abs(wt).le.VOFTOL) then
         wt=zero
        else if (abs(wt-one).le.VOFTOL) then
         wt=one
        else if ((wt.gt.zero).and.(wt.lt.one)) then
         if ((abs(intlo-xf(-1)).lt.VOFTOL*dxf(1)).or. &
             (abs(inthi-xf(1)).lt.VOFTOL*dxf(1))) then
          ! do nothing
         else
          print *,"coarse cell cannot be contained within fine cell"
          stop
         endif

        else
         print *,"wt invalid 3"
         stop
        endif

       endif

      else
       print *,"bfact_c invalid"
       stop
      endif

      return
      end subroutine intersect_weightMAC_interp

      subroutine gridsten_level(x,i,j,k,level,nhalf)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(out) :: x(-nhalf:nhalf,SDIM)
      INTEGER_T, intent(in) :: i,j,k,level
      INTEGER_T isten,dir
      INTEGER_T dummy_input
 
      if (nhalf.lt.0) then
       print *,"nhalf invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.cache_max_level)) then
       print *,"level invalid gridsten_level"
       stop
      endif
      if ((2*i-nhalf.lt.cache_index_low).or. &
          (2*i+nhalf.gt.cache_index_high)) then
       print *,"i out of range gridsten level"
       print *,"i,nhalf,level ",i,nhalf,level
       print *,"cache_index_low ",cache_index_low
       print *,"cache_index_high ",cache_index_high

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:5613"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read (*,*) dummy_input
       stop
      endif
      if ((2*j-nhalf.lt.cache_index_low).or. &
          (2*j+nhalf.gt.cache_index_high)) then
       print *,"j out of range"
       stop
      endif
      if (SDIM.eq.3) then
       if ((2*k-nhalf.lt.cache_index_low).or. &
           (2*k+nhalf.gt.cache_index_high)) then
        print *,"k out of range"
        stop
       endif
      endif
      do isten=-nhalf,nhalf
       dir=1
       x(isten,dir)=grid_cache(level,2*i+isten,dir)
       dir=2
       x(isten,dir)=grid_cache(level,2*j+isten,dir)
       if (SDIM.eq.3) then
        dir=SDIM
        x(isten,dir)=grid_cache(level,2*k+isten,dir)
       endif
      enddo
       
      return
      end subroutine gridsten_level

       ! 1<=dir<=sdim
      subroutine gridsten1D_level(x,i,level,dir,nhalf)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhalf,dir
      REAL_T, intent(out) :: x(-nhalf:nhalf)
      INTEGER_T, intent(in) :: i,level
      INTEGER_T isten
 
      if (nhalf.lt.0) then
       print *,"nhalf invalid"
       stop
      endif
      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid gridsten 1d level"
       stop
      endif
      if ((level.lt.0).or.(level.gt.cache_max_level)) then
       print *,"level invalid gridsten1D_level"
       stop
      endif
      if ((2*i-nhalf.lt.cache_index_low).or. &
          (2*i+nhalf.gt.cache_index_high)) then
       print *,"i out of range gridsten1D_level"
       print *,"i,nhalf,level ",i,nhalf,level
       stop
      endif
      do isten=-nhalf,nhalf
       x(isten)=grid_cache(level,2*i+isten,dir)
      enddo
       
      return
      end subroutine gridsten1D_level



       ! -1<=grid_type<=5
      subroutine gridstenMAC_level(x,i,j,k,level,nhalf,grid_type,caller_id)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: caller_id
      INTEGER_T, intent(in) :: nhalf
      INTEGER_T, intent(in) :: grid_type
      REAL_T, intent(out) :: x(-nhalf:nhalf,SDIM)
      INTEGER_T, intent(in) :: i,j,k,level
      INTEGER_T isten,dir,ii,jj,kk
      INTEGER_T box_type(SDIM)
      INTEGER_T dummy_input

      if ((grid_type.ge.-1).and.(grid_type.le.5)) then
       ! do nothing
      else
       print *,"grid_type invalid gridstenMAC_level: grid_type,caller_id", &
          grid_type,caller_id
       stop
      endif

      call grid_type_to_box_type(grid_type,box_type)
      ii=0
      jj=0
      kk=0
      ii=box_type(1)
      jj=box_type(2)
      if (SDIM.eq.3) then
       kk=box_type(SDIM)
      else if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif

      if (nhalf.lt.0) then
       print *,"nhalf invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.cache_max_level)) then
       print *,"level invalid gridstenMAC_level"
       print *,"level,cache_max_level ",level,cache_max_level
       print *,"nhalf,grid_type,caller_id ", &
               nhalf,grid_type,caller_id
       print *,"x(0,1..sdim),i,j,k ",x(0,1),x(0,2),x(0,SDIM),i,j,k
       print *,"box_type ",box_type(1),box_type(2),box_type(SDIM)

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:5735"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read (*,*) dummy_input
       stop
      endif
      if (1.eq.0) then
       print *,"before: i,j,k,ii,jj,kk,nhalf ",i,j,k,ii,jj,kk,nhalf
      endif
      if ((2*i-ii-nhalf.lt.cache_index_low).or. &
          (2*i-ii+nhalf.gt.cache_index_high)) then
       print *,"i out of range gridstenMAC_level"
       print *,"i,nhalf,level,grid_type ",i,nhalf,level,grid_type
       stop
      endif
      if ((2*j-jj-nhalf.lt.cache_index_low).or. &
          (2*j-jj+nhalf.gt.cache_index_high)) then
       print *,"j out of range"
       stop
      endif
      if (SDIM.eq.3) then
       if ((2*k-kk-nhalf.lt.cache_index_low).or. &
           (2*k-kk+nhalf.gt.cache_index_high)) then
        print *,"k out of range"
        stop
       endif
      endif
      if (1.eq.0) then
       print *,"after: i,j,k,ii,jj,kk,nhalf ",i,j,k,ii,jj,kk,nhalf
      endif
      do isten=-nhalf,nhalf
       dir=1
       x(isten,dir)=grid_cache(level,2*i-ii+isten,dir)
       dir=2
       x(isten,dir)=grid_cache(level,2*j-jj+isten,dir)
       if (SDIM.eq.3) then
        dir=SDIM
        x(isten,dir)=grid_cache(level,2*k-kk+isten,dir)
       endif
      enddo ! isten=-nhalf,nhalf
       
      return
      end subroutine gridstenMAC_level


      subroutine gridstenND_level(x,i,j,k,level,nhalf)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(out) :: x(-nhalf:nhalf,SDIM)
      INTEGER_T, intent(in) :: i,j,k,level
      INTEGER_T isten,dir

      if (nhalf.lt.0) then
       print *,"nhalf invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.cache_max_level)) then
       print *,"level invalid gridstenND_level"
       stop
      endif
      if ((2*i-1-nhalf.lt.cache_index_low).or. &
          (2*i-1+nhalf.gt.cache_index_high)) then
       print *,"i out of range gridstenND_level"
       print *,"i,nhalf,level ",i,nhalf,level
       stop
      endif
      if ((2*j-1-nhalf.lt.cache_index_low).or. &
          (2*j-1+nhalf.gt.cache_index_high)) then
       print *,"j out of range"
       stop
      endif
      if (SDIM.eq.3) then
       if ((2*k-1-nhalf.lt.cache_index_low).or. &
           (2*k-1+nhalf.gt.cache_index_high)) then
        print *,"k out of range"
        stop
       endif
      endif
      do isten=-nhalf,nhalf
       dir=1
       x(isten,dir)=grid_cache(level,2*i-1+isten,dir)
       dir=2
       x(isten,dir)=grid_cache(level,2*j-1+isten,dir)
       if (SDIM.eq.3) then
        dir=SDIM
        x(isten,dir)=grid_cache(level,2*k-1+isten,dir)
       endif
      enddo  ! isten
       
      return
      end subroutine gridstenND_level

       ! 1<=normdir<=sdim
      subroutine gridsten1DMAC_level(x,i,level,normdir,nhalf)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhalf,normdir
      REAL_T, intent(out) :: x(-nhalf:nhalf)
      INTEGER_T, intent(in) :: i,level
      INTEGER_T isten,ii

      ii=1
      if ((normdir.lt.1).or.(normdir.gt.SDIM)) then
       print *,"normdir invalid"
       stop
      endif 
      if (nhalf.lt.0) then
       print *,"nhalf invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.cache_max_level)) then
       print *,"level invalid gridsten1DMAC_level"
       stop
      endif
      if ((2*i-ii-nhalf.lt.cache_index_low).or. &
          (2*i-ii+nhalf.gt.cache_index_high)) then
       print *,"i out of range gridsten1DMAC_level"
       print *,"i,nhalf,level,normdir ",i,nhalf,level,normdir
       stop
      endif
      do isten=-nhalf,nhalf
       x(isten)=grid_cache(level,2*i-ii+isten,normdir)
      enddo
       
      return
      end subroutine gridsten1DMAC_level


         ! 1<=dir<=sdim
      subroutine intersect_weightMAC_avg(ic,i,bfact_c,bfact_f,wt)
      use probcommon_module

      INTEGER_T ic,i,bfact_c,bfact_f
      REAL_T wt
      INTEGER_T nhalf,dir_index
      INTEGER_T fablo(SDIM)
      REAL_T intlo,inthi
      REAL_T dxc(SDIM)
      REAL_T dxf(SDIM)
      REAL_T xlo(SDIM)
      REAL_T xc(-1:1)
      REAL_T xf(-1:1)

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      if (bfact_c.eq.1) then ! intersect_weightMAC_avg
       if (2*ic.eq.i) then
        wt=one
       else
        wt=zero
       endif
      else if (bfact_c.gt.1) then  
       if (bfact_c*(ic/bfact_c).eq.ic) then
        if (2*ic.eq.i) then
         wt=one
        else 
         wt=zero
        endif
       else

        nhalf=1

        dxc(1)=one
        dxf(1)=half

        xlo(1)=zero
        dir_index=1
        fablo(1)=0
        call gridsten1DMAC(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
        call gridsten1DMAC(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

        if (bfact_f.eq.1) then
         if (abs(xf(1)-xf(-1)-dxf(1)).ge.VOFTOL) then
          print *,"xf invalid 1"
          stop
         endif
         if (abs(xf(1)+xf(-1)-two*xf(0)).ge.VOFTOL) then
          print *,"xf invalid 2"
          stop
         endif
        else if (bfact_f.gt.1) then
         if ((xf(1)-xf(0).le.zero).or. &
             (xf(0)-xf(-1).le.zero)) then
          print *,"xf invalid 3"
          stop
         endif
        else
         print *,"bfact_f invalid"
         stop
        endif

        if (bfact_c.eq.1) then
         if (abs(xc(1)-xc(-1)-dxc(1)).ge.VOFTOL) then
          print *,"xc invalid"
          stop
         endif
         if (abs(xc(1)+xc(-1)-two*xc(0)).ge.VOFTOL) then
          print *,"xc invalid"
          stop
         endif
        else if (bfact_c.gt.1) then
         if ((xc(1)-xc(0).le.zero).or. &
             (xc(0)-xc(-1).le.zero)) then
          print *,"xc invalid 3"
          stop
         endif
        else
         print *,"bfact_c invalid"
         stop
        endif

        intlo=max(xf(-1),xc(-1))
        inthi=min(xf(1),xc(1))
        if (intlo.ge.inthi) then
         wt=zero
        else
         wt=(inthi-intlo)/(xc(1)-xc(-1))  ! average down
        endif
        if (abs(wt).le.VOFTOL) then
         wt=zero
        else if (abs(wt-one).le.VOFTOL) then
         wt=one
        else if ((wt.gt.zero).and.(wt.lt.one)) then
         ! do nothing
        else
         print *,"wt invalid 4"
         stop
        endif

       endif
      else
       print *,"bfact_c invalid"
       stop
      endif

      return
      end subroutine intersect_weightMAC_avg

       ! either a/b or (a+1)/b-1
      INTEGER_T function DIV_FLOOR(a,b)
      IMPLICIT NONE

      INTEGER_T a,b,apos,div1

      if (b.le.0) then
       print *,"b invalid"
       stop
      endif
      if (a.ge.0) then
       DIV_FLOOR=a/b
      else
       apos=-(a+1)
       div1=apos/b
       DIV_FLOOR=-div1-1
      endif

      return
      end function DIV_FLOOR
 
      subroutine fine_subelement_stencil( &
       ic,jc,kc,stenlo,stenhi,bfact_c,bfact_f)
      IMPLICIT NONE

      INTEGER_T ic,jc,kc
      INTEGER_T coarse_index(3)
      INTEGER_T stenlo(3),stenhi(3)
      INTEGER_T bfact_c,bfact_f
      INTEGER_T dir2

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      stenlo(3)=0
      stenhi(3)=0

      coarse_index(1)=ic
      coarse_index(2)=jc
      coarse_index(3)=kc

      do dir2=1,SDIM
       if (bfact_c.eq.1) then
        stenlo(dir2)=2*coarse_index(dir2)
        stenhi(dir2)=stenlo(dir2)+1
       else
        stenlo(dir2)=DIV_FLOOR(coarse_index(dir2),bfact_c)
        stenlo(dir2)=stenlo(dir2)*bfact_c*2
        stenhi(dir2)=stenlo(dir2)+bfact_c*2-1
       endif
      enddo ! dir2
      
      return
      end subroutine fine_subelement_stencil


      subroutine coarse_subelement_stencil( &
       ifine,jfine,kfine,stenlo,stenhi,bfact_c,bfact_f)
      IMPLICIT NONE

      INTEGER_T ifine,jfine,kfine
      INTEGER_T fine_index(3)
      INTEGER_T stenlo(3),stenhi(3)
      INTEGER_T bfact_c,bfact_f
      INTEGER_T dir2

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      stenlo(3)=0
      stenhi(3)=0

      fine_index(1)=ifine
      fine_index(2)=jfine
      fine_index(3)=kfine

      do dir2=1,SDIM
       if (bfact_c.eq.1) then
        stenlo(dir2)=DIV_FLOOR(fine_index(dir2),2)
        stenhi(dir2)=stenlo(dir2)
       else
        stenlo(dir2)=DIV_FLOOR(fine_index(dir2),2*bfact_c)
        stenlo(dir2)=stenlo(dir2)*bfact_c
        stenhi(dir2)=stenlo(dir2)+bfact_c-1
       endif
      enddo ! dir2
      
      return
      end subroutine coarse_subelement_stencil

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      subroutine grid_type_to_box_type(grid_type,box_type)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: grid_type
      INTEGER_T, intent(out) :: box_type(SDIM)
      INTEGER_T dir

      do dir=1,SDIM
       box_type(dir)=0  ! default to CELL
      enddo
      if (grid_type.eq.-1) then
       ! do nothing
      else if ((grid_type.ge.0).and. &
               (grid_type.lt.SDIM)) then
       box_type(grid_type+1)=1  ! NODE
      else if (grid_type.eq.3) then
       box_type(1)=1 ! NODE
       box_type(2)=1 ! NODE
      else if ((grid_type.eq.4).and.(SDIM.eq.3)) then
       box_type(1)=1 ! NODE
       box_type(SDIM)=1 ! NODE
      else if ((grid_type.eq.5).and.(SDIM.eq.3)) then
       box_type(2)=1 ! NODE
       box_type(SDIM)=1 ! NODE
      else
       print *,"grid_type invalid"
       stop
      endif
     
      return 
      end subroutine grid_type_to_box_type

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      subroutine box_type_to_grid_type(grid_type,box_type)
      IMPLICIT NONE

      INTEGER_T, intent(out) :: grid_type
      INTEGER_T, intent(in) :: box_type(SDIM)

      if ((box_type(1).eq.0).and. &
          (box_type(2).eq.0).and. &
          (box_type(SDIM).eq.0)) then
          grid_type=-1  ! cell centered
      else if ((box_type(1).eq.1).and. &
               (box_type(2).eq.0).and. &
               (box_type(SDIM).eq.0)) then
       grid_type=0  !UMAC
      else if ((box_type(1).eq.0).and. &
               (box_type(2).eq.1).and. &
               (SDIM.eq.2)) then
       grid_type=1  !VMAC
      else if ((box_type(1).eq.0).and. &
               (box_type(2).eq.1).and. &
               (box_type(SDIM).eq.0).and. &
               (SDIM.eq.3)) then
       grid_type=1  !VMAC
      else if ((box_type(1).eq.0).and. &
               (box_type(2).eq.0).and. &
               (box_type(SDIM).eq.1).and. &
               (SDIM.eq.3)) then
       grid_type=2  !WMAC
      else if ((box_type(1).eq.1).and. &
               (box_type(2).eq.1).and. &
               (box_type(SDIM).eq.0).and. &
               (SDIM.eq.3)) then
       grid_type=3  !XY
      else if ((box_type(1).eq.1).and. &
               (box_type(2).eq.1).and. &
               (SDIM.eq.2)) then
       grid_type=3  !XY
      else if ((box_type(1).eq.1).and. &
               (box_type(2).eq.0).and. &
               (box_type(SDIM).eq.1).and. &
               (SDIM.eq.3)) then
       grid_type=4  !XZ
      else if ((box_type(1).eq.0).and. &
               (box_type(2).eq.1).and. &
               (box_type(SDIM).eq.1).and. &
               (SDIM.eq.3)) then
       grid_type=5  !YZ
      else
       print *,"box_type not recognizable"
       stop
      endif

      return 
      end subroutine box_type_to_grid_type


      subroutine coarse_subelement_stencilMAC( &
       ifine,jfine,kfine,stenlo,stenhi,bfact_c,bfact_f, &
       grid_type)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: grid_type
      INTEGER_T, intent(in) :: ifine,jfine,kfine
      INTEGER_T fine_index(3)
      INTEGER_T, intent(out) :: stenlo(3),stenhi(3)
      INTEGER_T, intent(in) :: bfact_c,bfact_f
      INTEGER_T dir2,denom
      INTEGER_T :: box_type(SDIM)

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      call grid_type_to_box_type(grid_type,box_type)

      stenlo(3)=0
      stenhi(3)=0

      fine_index(1)=ifine
      fine_index(2)=jfine
      fine_index(3)=kfine

      do dir2=1,SDIM
       if (bfact_c.eq.1) then
        stenlo(dir2)=DIV_FLOOR(fine_index(dir2),2)
        stenhi(dir2)=stenlo(dir2)
        if (box_type(dir2).eq.1) then ! NODE
         if (2*(fine_index(dir2)/2).ne.fine_index(dir2)) then
          stenhi(dir2)=stenhi(dir2)+1
         endif
        else if (box_type(dir2).eq.0) then ! CELL
         ! do nothing
        else
         print *,"box_type invalid"
         stop
        endif
       else if (bfact_c.gt.1) then
        stenlo(dir2)=DIV_FLOOR(fine_index(dir2),2*bfact_c)
        stenlo(dir2)=stenlo(dir2)*bfact_c
        stenhi(dir2)=stenlo(dir2)+bfact_c-1
        if (box_type(dir2).eq.1) then ! NODE
         denom=2*bfact_c
         if (denom*(fine_index(dir2)/denom).eq.fine_index(dir2)) then
          stenlo(dir2)=fine_index(dir2)/2
          stenhi(dir2)=stenlo(dir2)
         else
          stenhi(dir2)=stenhi(dir2)+1
         endif
        else if (box_type(dir2).eq.0) then ! CELL
         ! do nothing
        else
         print *,"box_type invalid"
         stop
        endif
       else
        print *,"bfact_c invalid"
        stop
       endif
      enddo ! dir2
      
      return
      end subroutine coarse_subelement_stencilMAC


        ! 1<=dir<=sdim
      subroutine fine_subelement_stencilMAC( &
       ic,jc,kc,stenlo,stenhi,bfact_c,bfact_f,grid_type)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: ic,jc,kc
      INTEGER_T, intent(in) :: grid_type
      INTEGER_T coarse_index(3)
      INTEGER_T, intent(out) :: stenlo(3),stenhi(3)
      INTEGER_T, intent(in) :: bfact_c,bfact_f
      INTEGER_T dir2
      INTEGER_T :: box_type(SDIM)

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      call grid_type_to_box_type(grid_type,box_type)

      stenlo(3)=0
      stenhi(3)=0

      coarse_index(1)=ic
      coarse_index(2)=jc
      coarse_index(3)=kc

      do dir2=1,SDIM
       if (bfact_c.eq.1) then
        stenlo(dir2)=2*coarse_index(dir2)
        stenhi(dir2)=stenlo(dir2)+1
        if (box_type(dir2).eq.1) then ! NODE
         stenhi(dir2)=stenlo(dir2)
        else if (box_type(dir2).eq.0) then ! CELL
         ! do nothing
        else
         print *,"box_type invalid"
         stop
        endif
       else if (bfact_c.gt.1) then
        stenlo(dir2)=DIV_FLOOR(coarse_index(dir2),bfact_c)
        stenlo(dir2)=stenlo(dir2)*bfact_c*2
        stenhi(dir2)=stenlo(dir2)+bfact_c*2-1
        if (box_type(dir2).eq.1) then ! NODE
         if (bfact_c*(coarse_index(dir2)/bfact_c).eq.coarse_index(dir2)) then
          stenlo(dir2)=2*coarse_index(dir2)
          stenhi(dir2)=stenlo(dir2)
         else
          stenlo(dir2)=DIV_FLOOR(coarse_index(dir2),bfact_c)
          stenlo(dir2)=stenlo(dir2)*bfact_c*2
          stenhi(dir2)=stenlo(dir2)+bfact_c*2
         endif
        else if (box_type(dir2).eq.0) then ! CELL
         ! do nothing
        else
         print *,"box_type invalid"
         stop
        endif
       else
        print *,"bfact_c invalid"
        stop
       endif
      enddo ! dir2
      
      return
      end subroutine fine_subelement_stencilMAC


       ! -1<=grid_type<=5
      subroutine growntileboxMAC( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,ng,grid_type,caller_id)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: caller_id
      INTEGER_T, intent(in) :: grid_type
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(out) :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: ng
      INTEGER_T dir2
      INTEGER_T :: box_type(SDIM)

      growlo(3)=0
      growhi(3)=0

      if ((grid_type.ge.-1).and.(grid_type.le.5)) then
       ! do nothing
      else
       print *,"grid_type invalid growntilebox mac"
       print *,"ng=",ng
       print *,"grid_type=",grid_type
       print *,"caller_id=",caller_id
       stop
      endif 

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       growlo(dir2)=tilelo(dir2)
       if (tilelo(dir2).eq.fablo(dir2)) then
        growlo(dir2)=tilelo(dir2)-ng
       else if (tilelo(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilelo(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
       growhi(dir2)=tilehi(dir2)
       if (tilehi(dir2).eq.fabhi(dir2)) then
        growhi(dir2)=tilehi(dir2)+ng
       else if (tilehi(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilehi(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
      enddo ! dir2=1..sdim

      do dir2=1,SDIM
        ! NODE
       if (box_type(dir2).eq.1) then
        if (tilehi(dir2).eq.fabhi(dir2)) then
         growhi(dir2)=growhi(dir2)+1
        endif
       else if (box_type(dir2).eq.0) then ! CELL
        ! do nothing
       else
        print *,"box_type(dir2) invalid"
        stop
       endif
      enddo !dir2=1..sdim

      return
      end subroutine growntileboxMAC


      subroutine growntileboxNODE( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,ng)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(out) :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: ng
      INTEGER_T dir2

      growlo(3)=0
      growhi(3)=0

      do dir2=1,SDIM
       growlo(dir2)=tilelo(dir2)
       if (tilelo(dir2).eq.fablo(dir2)) then
        growlo(dir2)=tilelo(dir2)-ng
       else if (tilelo(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilelo(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
       growhi(dir2)=tilehi(dir2)
       if (tilehi(dir2).eq.fabhi(dir2)) then
        growhi(dir2)=tilehi(dir2)+ng+1
       else if (tilehi(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilehi(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
      enddo ! dir2
      
      return
      end subroutine growntileboxNODE


      subroutine tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)
      IMPLICIT NONE

      INTEGER_T ux,uy,uz
      INTEGER_T vx,vy,vz
      INTEGER_T wx,wy,wz

      ux=1
      vx=ux+1

      if (SDIM.eq.3) then
       wx=vx+1
      else if (SDIM.eq.2) then
       wx=vx
      else
       print *,"dimension bust"
       stop
      endif

      uy=wx+1
      vy=uy+1

      if (SDIM.eq.3) then
       wy=vy+1
      else if (SDIM.eq.2) then
       wy=vy
      else
       print *,"dimension bust"
       stop
      endif

      uz=wy+1
      vz=uz+1
      wz=vz+1

      return
      end subroutine tensorcomp_matrix

      subroutine growntileboxTENSOR( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,dir)
      IMPLICIT NONE

      INTEGER_T dir
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T dir2

      growlo(3)=0
      growhi(3)=0

      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir invalid growntilebox tensor"
       stop
      endif 
      do dir2=1,SDIM
       growlo(dir2)=tilelo(dir2)
       if (tilelo(dir2).eq.fablo(dir2)) then
        growlo(dir2)=tilelo(dir2)
       else if (tilelo(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilelo(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
       growhi(dir2)=tilehi(dir2)
       if (tilehi(dir2).eq.fabhi(dir2)) then
        growhi(dir2)=tilehi(dir2)
       else if (tilehi(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilehi(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
      enddo ! dir2

      if (tilehi(dir+1).eq.fabhi(dir+1)) then
       growhi(dir+1)=growhi(dir+1)+1
      endif

      do dir2=1,SDIM
       if (dir2.ne.dir+1) then
        if (tilelo(dir2).eq.fablo(dir2)) then
         growlo(dir2)=growlo(dir2)-1
        endif
        if (tilehi(dir2).eq.fabhi(dir2)) then
         growhi(dir2)=growhi(dir2)+1
        endif
       endif
      enddo ! dir2
      
      return
      end subroutine growntileboxTENSOR


      subroutine growntileboxTENSOR_SEM( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,dir)
      IMPLICIT NONE

      INTEGER_T dir
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T dir2

      growlo(3)=0
      growhi(3)=0

      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir invalid growntilebox tensor sem"
       stop
      endif 
      do dir2=1,SDIM
       growlo(dir2)=tilelo(dir2)
       if (tilelo(dir2).eq.fablo(dir2)) then
        growlo(dir2)=tilelo(dir2)
       else if (tilelo(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilelo(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
       growhi(dir2)=tilehi(dir2)
       if (tilehi(dir2).eq.fabhi(dir2)) then
        growhi(dir2)=tilehi(dir2)
       else if (tilehi(dir2).gt.fabhi(dir2)) then
        print *,"tile box incorrect"
        stop
       else if (tilehi(dir2).lt.fablo(dir2)) then
        print *,"tile box incorrect"
        stop
       endif
      enddo ! dir2
      do dir2=1,SDIM
       if (dir2.ne.dir+1) then
        if (tilelo(dir2).eq.fablo(dir2)) then
         growlo(dir2)=growlo(dir2)-1
        endif
        if (tilehi(dir2).eq.fabhi(dir2)) then
         growhi(dir2)=growhi(dir2)+1
        endif
       endif
      enddo ! dir2
      
      return
      end subroutine growntileboxTENSOR_SEM

      subroutine strip_status(i,j,k,bfact,stripstat)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: i,j,k,bfact
      INTEGER_T, intent(out) :: stripstat

      if (bfact.lt.1) then
       print *,"bfact invalid32"
       stop
      endif

      stripstat=1

      if ((i/bfact)*bfact.ne.i) then
       stripstat=0
      endif
      if ((j/bfact)*bfact.ne.j) then
       stripstat=0
      endif
      if (SDIM.eq.3) then
       if ((k/bfact)*bfact.ne.k) then
        stripstat=0
       endif
      endif

      return
      end subroutine strip_status

      subroutine elementbox(i,j,k,bfact,dir,elemlo,elemhi)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: i,j,k,bfact,dir
      INTEGER_T dir2
      INTEGER_T, intent(out) :: elemlo(3)
      INTEGER_T, intent(out) :: elemhi(3)

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid in elementbox"
       stop
      endif

      elemlo(1)=i
      elemlo(2)=j
      elemlo(3)=0
      elemhi(3)=0

      if (SDIM.eq.3) then
       elemlo(SDIM)=k
      endif
      do dir2=1,SDIM
       elemhi(dir2)=elemlo(dir2)+bfact-1
      enddo
      elemhi(dir+1)=elemlo(dir+1)

      return
      end subroutine elementbox

      subroutine strip_status_dir(i,j,k,bfact,dir,stripstat)
      IMPLICIT NONE

      INTEGER_T, intent(in)  :: i,j,k,bfact,dir
      INTEGER_T, intent(out) :: stripstat

      if (bfact.lt.2) then
       print *,"bfact invalid33"
       stop
      endif

      stripstat=1

      if (dir.eq.0) then
       if ((i/bfact)*bfact.ne.i) then
        stripstat=0
       endif
      else if (dir.eq.1) then
       if ((j/bfact)*bfact.ne.j) then
        stripstat=0
       endif
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       if ((k/bfact)*bfact.ne.k) then
        stripstat=0
       endif
      else
       print *,"dir invalid strip status dir"
       stop
      endif

      return
      end subroutine strip_status_dir



!============================
! Spectral deferred correction Gauss Quadrature weights
!===========================

! for interval k=1..tbfact : 
!  f(xpt(i_sub))=sum_{jsten=0}^{tbfact} GQws(jsten,isub,k)f_jsten
! integral is:
!  sum_{iquad=0}^{tbfact-1} w_{iquad} 
!   sum_{jsten=0}^{tbfact} GQws(jsten,iquad,k)f_jsten=
!  sum_{jsten=0}^{tbfact} sum_{iquad=0}^{tbfact-1} 
!    w_{iquad}GQws(jsten,iquad,k)f_{jsten}

      subroutine SDC_GQweights(qbfact,tbfact,GQws)
      use LagrangeInterpolationPolynomial
      use LegendreNodes

      IMPLICIT NONE
      INTEGER_T qbfact,tbfact
      INTEGER_T i1,j1,k
      REAL_T GQws(0:tbfact,0:qbfact-1,1:tbfact) 
      REAL_T y(0:qbfact-1)
      REAL_T yGL(0:tbfact)
      REAL_T, dimension(:),allocatable :: bwGL,l
      REAL_T xpt

      do i1=0,qbfact-1
       y(i1)=cache_gauss(qbfact,i1,TMTYPE)
      enddo
      do i1=0,tbfact
       yGL(i1)=cache_gauss_lobatto(tbfact,i1,TMTYPE)
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

       ! xfine is relative to the lower left hand corner of the
       ! grid element.
      subroutine SEM_INTERP_ELEMENT( &
       nvar,bfact,grid_type, &
       stenhi,dx,xfine,fcoarse,fxfine,caller_id)
      use LagrangeInterpolationPolynomial
      use LegendreNodes

      IMPLICIT NONE
      INTEGER_T, intent(in) :: nvar
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: grid_type ! -1:ggg; 0:lgg; 1:glg; 2:ggl
      INTEGER_T, intent(in) :: stenhi(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xfine(SDIM)
      REAL_T, intent(in) ::  &
        fcoarse(D_DECL(0:stenhi(1),0:stenhi(2),0:stenhi(3)),nvar)
      REAL_T, intent(out) :: fxfine(nvar)
 
      INTEGER_T i1
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      INTEGER_T khi
      INTEGER_T dir
      INTEGER_T n
      REAL_T y(0:bfact-1)
      REAL_T yGL(0:bfact)
      REAL_T wt
      REAL_T, dimension(:),allocatable :: temp,tempGL
      REAL_T, dimension(:),allocatable :: ypoints,bwG
      REAL_T, dimension(:),allocatable :: ypointsGL,bwGL
      REAL_T, dimension(:,:),allocatable :: lg
      REAL_T INTERP_TOL
      REAL_T wtsum
      REAL_T local_data
      INTEGER_T caller_id
      INTEGER_T :: box_type(SDIM)

      INTERP_TOL=1.0D-10

      call grid_type_to_box_type(grid_type,box_type)

      ii=0
      jj=0
      kk=0
      ii=box_type(1)
      jj=box_type(2)
      if (SDIM.eq.3) then
       kk=box_type(SDIM)
      endif

      if (nvar.le.0) then
       print *,"nvar invalid"
       stop
      endif
 
      if (bfact.lt.2) then
       print *,"Not valid order for high order interpolation"
       stop
      endif

      if ((SDIM.ne.2).and.(SDIM.ne.3)) then
       print *,"dimension bust"
       stop
      endif 

      do dir=1,SDIM
       if (box_type(dir).eq.1) then
        if (stenhi(dir).ne.bfact) then
         print *,"stenhi invalid"
         stop
        endif
       else if (box_type(dir).eq.0) then
        if (stenhi(dir).ne.bfact-1) then
         print *,"stenhi invalid"
         stop
        endif
       else
        print *,"box_type invalid"
        stop
       endif
       if (abs(dx(dir)).le.1.0E+20) then
        ! do nothing
       else
        print *,"abs(dx(dir)) overflow"
        stop
       endif
      enddo ! dir=1..sdim


      do i1=0,bfact-1
       y(i1)=cache_gauss(bfact,i1,SPTYPE)
       if (abs(y(i1)).le.1.01) then
        ! do nothing
       else
        print *,"abs(y(i1)) overflow"
        stop
       endif
      enddo
      do i1=0,bfact
       yGL(i1)=cache_gauss_lobatto(bfact,i1,SPTYPE)
      enddo
 
      allocate(lg(0:bfact,SDIM))

      allocate(ypoints(0:bfact-1))
      allocate(bwG(0:bfact-1))
  
      do dir=1,SDIM
       do i=0,bfact-1
        ypoints(i)=(y(i)+one)*bfact*dx(dir)/two
        if (abs(ypoints(i)).le.(1.01)*bfact*dx(dir)) then
         ! do nothing
        else
         print *,"abs(ypoints(i)) overflow"
         stop
        endif
       enddo
       allocate(temp(0:bfact-1))
       call BarycentricWeights(bfact-1,ypoints,bwG)
       if (abs(xfine(dir)).le.1.0E+20) then
        call LagrangeInterpolatingPolynomial(bfact-1, &
         xfine(dir),ypoints,bwG,temp)
        do i=0,bfact-1
         if (abs(temp(i)).le.1.0E+20) then
          lg(i,dir)=temp(i)
         else
          print *,"abs(temp(i)) overflow"
          print *,"i,dir,temp ",i,dir,temp(i)
          stop
         endif
        enddo ! i=0..bfact-1
       else
        print *,"abs(xfine(dir)) overflow"
        stop
       endif
       deallocate(temp)
      enddo ! dir=1..sdim

      deallocate(bwG)
      deallocate(ypoints)

      do dir=1,SDIM
       if (box_type(dir).eq.0) then
        ! do nothing
       else if (box_type(dir).eq.1) then
 
        allocate(ypointsGL(0:bfact))
        allocate(bwGL(0:bfact))
        allocate(tempGL(0:bfact))

        do i=0,bfact
         ypointsGL(i)=(yGL(i)+one)*bfact*dx(dir)/two
        enddo
        call BarycentricWeights(bfact,ypointsGL,bwGL)
        call LagrangeInterpolatingPolynomial(bfact, &
         xfine(dir),ypointsGL,bwGL,tempGL)
        do i=0,bfact
         lg(i,dir)=tempGL(i)
        enddo

        deallocate(tempGL)
        deallocate(bwGL)
        deallocate(ypointsGL)

       else
        print *,"box_type(dir) invalid"
        stop
       endif

      enddo ! dir=1..sdim

      if (SDIM.eq.2) then
       khi=0
      else if (SDIM.eq.3) then
       khi=stenhi(SDIM)
      else
       print *,"dimension bust"
       stop
      endif 

      do n=1,nvar 
       fxfine(n)=zero
      enddo
      wtsum=zero
      do k=0,khi
      do j=0,stenhi(2)
      do i=0,stenhi(1)
       wt=lg(i,1)*lg(j,2)
       if (SDIM.eq.3) then
        wt=wt*lg(k,SDIM)
       endif
       if (abs(wt).le.1.0E+20) then
        do n=1,nvar 
         local_data=fcoarse(D_DECL(i,j,k),n)
         if (abs(local_data).le.1.0E+20) then
          fxfine(n)=fxfine(n)+local_data*wt
         else
          print *,"abs(local_data) overflow SEM_INTERP_ELEMENT"
          print *,"n,nvar,local_data ",n,nvar,local_data
          print *,"caller_id=",caller_id
          stop
         endif
        enddo ! n=1..nvar
       else
        print *,"abs(wt) overflow"
        print *,"i,j,k,stenhi ",i,j,k,stenhi(1),stenhi(2),khi
        stop
       endif
       wtsum=wtsum+wt
      enddo
      enddo
      enddo
      if (abs(wtsum-one).gt.INTERP_TOL) then
       print *,"wtsum invalid"
       print *,"wtsum=",wtsum
       print *,"nvar=",nvar
       print *,"bfact=",bfact
       print *,"grid_type=",grid_type
       do dir=1,SDIM
        print *,"dir,xfine(dir),dx(dir)*bfact ",dir, &
         xfine(dir),dx(dir)*bfact
       enddo
       do i1=0,stenhi(1)
        print *,"i1,lg(i1,1) ",i1,lg(i1,1)
       enddo
       do i1=0,stenhi(2)
        print *,"i1,lg(i1,2) ",i1,lg(i1,2)
       enddo
       do i1=0,stenhi(SDIM)
        print *,"i1,lg(i1,sdim) ",i1,lg(i1,SDIM)
       enddo
       stop
      endif
      do n=1,nvar 
       fxfine(n)=fxfine(n)/wtsum
      enddo

      deallocate(lg)

      return
      end subroutine SEM_INTERP_ELEMENT


      subroutine lineGRAD( &
       levelrz_in, &
       dir, &
       nc, &
       RRface, &
       bctype, &
       bcvalue, &
       vel, &  ! MAC variable
       source, &
       source_side, &
       dest_grad, &
       dest_interp, &
       bfact, &
       dx, &
       x_sep, &
       operation_flag)
      use LegendreNodes
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: levelrz_in
      INTEGER_T, intent(in) :: dir,nc
      REAL_T, intent(in) :: RRface(0:bfact)
    
      INTEGER_T, intent(in) :: operation_flag
      REAL_T, intent(in) :: dx
      INTEGER_T, intent(in) :: bctype(2)
      INTEGER_T local_bctype(2)
      REAL_T, intent(in) :: x_sep(2)
      REAL_T, intent(in) :: bcvalue(2)
      REAL_T, intent(in) :: vel(0:bfact)
      REAL_T, intent(in) :: source_side(2)
      REAL_T local_source_side(2)
      REAL_T, intent(in) :: source(0:bfact-1) 
      REAL_T, intent(out) :: dest_grad(0:bfact)
      REAL_T, intent(out) ::  dest_interp(0:bfact)

      INTEGER_T i1,j1,isten
      REAL_T y(0:bfact-1)
      REAL_T y_extend(0:bfact+1)
      REAL_T yGL(0:bfact)  ! Gauss-Lobatto Legendre points
      REAL_T yGL_extend(0:bfact+1)
      REAL_T yLT(0:bfact+1)
      REAL_T yRT(0:bfact+1)
      REAL_T yLRT(0:bfact+1)
      REAL_T yLRTextrap(0:bfact-1)
      REAL_T yLTextrap(0:bfact)
      REAL_T yRTextrap(0:bfact)

      REAL_T yLTextrapEXT(0:bfact)
      REAL_T yRTextrapEXT(0:bfact)

      REAL_T wMAT(0:bfact-1,0:bfact-1)
      REAL_T wMATGL(0:bfact,0:bfact)
      REAL_T wMAT_extend(0:bfact+1,0:bfact+1)
      REAL_T wLT(0:bfact+1,0:bfact+1)
      REAL_T wRT(0:bfact+1,0:bfact+1)

      REAL_T wLT_EXT(0:bfact,0:bfact)
      REAL_T wRT_EXT(0:bfact,0:bfact)

      REAL_T wLRT(0:bfact+1,0:bfact+1)
      REAL_T dx_element
      REAL_T dx_ends
      REAL_T sum1,sum2,A,B,C,D,det
      REAL_T PLINE(0:bfact+1)
      REAL_T PLINE2(0:bfact+1)
      INTEGER_T AMR_boundary_flag

      if ((operation_flag.lt.0).or.(operation_flag.gt.11)) then
       print *,"operation_flag invalid1"
       stop
      endif
      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid lineGRAD"
       stop
      endif
       ! in: lineGRAD
      if (bfact.lt.2) then
       print *,"bfact invalid34"
       stop
      endif
      if (levelrz_in.ne.levelrz) then
       print *,"levelrz_in invalid"
       stop
      endif

      dx_element=dx*bfact

      if (levelrz_in.eq.0) then
       ! do nothing
      else if (levelrz_in.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz_in.eq.3) then

       if (operation_flag.eq.6) then ! tensor derivatives

        if ((nc.lt.1).or. &
            (nc.gt.SDIM)) then
         print *,"nc invalid"
         stop
        endif

       else if (operation_flag.eq.0) then ! grad p_MAC
        ! do nothing
       else if (operation_flag.eq.2) then ! grad ppot_MAC/rho_pot_MAC
        ! do nothing
       else if (operation_flag.eq.3) then ! u^{Cell->Mac}
        ! do nothing
       else if (operation_flag.eq.5) then ! u^MAC=u^MAC+(DU)^{cell->mac}
        ! do nothing
       else if (operation_flag.eq.11) then!u^MAC=u^MAC+(DU)^{cell->mac}
        ! do nothing
       else if (operation_flag.eq.7) then ! advection
        ! do nothing
       else if (operation_flag.eq.8) then ! coupling
        ! do nothing
       else
        print *,"operation_flag invalid2"
        stop
       endif
      else
       print *,"levelrz_in invalid"
       stop
      endif

      do i1=0,bfact-1
       y(i1)=cache_gauss(bfact,i1,SPTYPE)
      enddo
      do i1=0,bfact
       yGL(i1)=cache_gauss_lobatto(bfact,i1,SPTYPE)
      enddo

      dx_ends=dx_element*(y(0)+one)/two
      if (dx_ends.ge.half*dx) then
       print *,"(dx_ends.gt.half*dx)"
       stop
      endif

      do j1=0,bfact-1
       do i1=0,bfact-1
        wMAT(i1,j1)=cache_wMAT(bfact,i1,j1,SPTYPE)
       enddo
      enddo
      do j1=0,bfact
       do i1=0,bfact
        wMATGL(i1,j1)=cache_wMATGL(bfact,i1,j1,SPTYPE)
       enddo
      enddo

      if ((one+y(0).gt.zero).and. &
          (one-y(bfact-1).gt.zero)) then
       y_extend(0)=-one-abs(y(0)+one)
       do i1=1,bfact
        y_extend(i1)=y(i1-1)
       enddo
       y_extend(bfact+1)=one+abs(one-y(bfact-1))
      else
       print *,"y(0) or y(bfact-1) invalid"
       stop
      endif

      do i1=0,bfact+1
       yGL_extend(i1)=cache_gauss_lobatto(bfact+1,i1,SPTYPE)
      enddo
      yLT(0)=-one
      yLT(bfact+1)=y_extend(bfact+1)
      yRT(0)=y_extend(0)
      yRT(bfact+1)=one
      yLRT(0)=-one
      yLRT(bfact+1)=one
      do i1=1,bfact
       yLT(i1)=y(i1-1)
       yRT(i1)=y(i1-1)
       yLRT(i1)=y(i1-1)
      enddo ! i1

      yLTextrap(bfact)=y_extend(bfact+1)
      yLTextrapEXT(bfact)=one
      do i1=0,bfact-1
       yLTextrap(i1)=y(i1)
       yLTextrapEXT(i1)=y(i1)
      enddo

      yRTextrap(0)=y_extend(0)
      yRTextrapEXT(0)=-one
      do i1=0,bfact-1
       yRTextrap(i1+1)=y(i1)
       yRTextrapEXT(i1+1)=y(i1)
      enddo

      do i1=0,bfact-1
       yLRTextrap(i1)=y(i1)
      enddo


       ! wij=Li'(xj)
!      call polyinterp_Dmatrix(bfact+1,yGL_extend,wMAT_extend)
!      call polyinterp_Dmatrix(bfact+1,yRT,wRT)
!      call polyinterp_Dmatrix(bfact+1,yLT,wLT)
!      call polyinterp_Dmatrix(bfact+1,yLRT,wLRT)

      do j1=0,bfact+1
       do i1=0,bfact+1
        wMAT_extend(i1,j1)=cache_wMAT_extend(bfact,i1,j1,SPTYPE)
        wRT(i1,j1)=cache_wRT(bfact,i1,j1,SPTYPE)
        wLT(i1,j1)=cache_wLT(bfact,i1,j1,SPTYPE)
        wLRT(i1,j1)=cache_wLRT(bfact,i1,j1,SPTYPE)
       enddo
      enddo

      do j1=0,bfact
       do i1=0,bfact
        wRT_EXT(i1,j1)=cache_wRT_EXT(bfact,i1,j1,SPTYPE) !E=left N=right
        wLT_EXT(i1,j1)=cache_wLT_EXT(bfact,i1,j1,SPTYPE) !N=left E=right
       enddo
      enddo


       ! at this point, PLINE(0) and PLINE(bfact+1) need to be filled.
      do i1=0,bfact-1
       PLINE(i1+1)=source(i1)
      enddo

       ! -1=extrap,0=interior,1=dirichlet,2=neumann,
       !  3=reflect_even,4=reflect_odd
       ! -2=coarse next to fine
       ! -3=fine next to coarse
      local_source_side(1)=source_side(1)
      local_source_side(2)=source_side(2)
      local_bctype(1)=bctype(1)
      local_bctype(2)=bctype(2)
      if (bctype(1).eq.3) then ! reflect even
       local_bctype(1)=0
       local_source_side(1)=PLINE(1)
      endif
      if (bctype(2).eq.3) then ! reflect even
       local_bctype(2)=0
       local_source_side(2)=PLINE(bfact)
      endif
      if (bctype(1).eq.4) then ! reflect odd
       local_bctype(1)=0
       local_source_side(1)=-PLINE(1)
      endif
      if (bctype(2).eq.4) then ! reflect odd
       local_bctype(2)=0
       local_source_side(2)=-PLINE(bfact)
      endif

       ! SEM_IMAGE_BC_ALG is defined in PROBCOMMON.F90
      if (SEM_IMAGE_BC_ALG.eq.1) then

       if (bctype(1).eq.1) then ! dirichlet -> reflect odd
        local_bctype(1)=0
        local_source_side(1)=two*bcvalue(1)-PLINE(1)
       endif
       if (bctype(2).eq.1) then ! dirichlet -> reflect odd
        local_bctype(2)=0
        local_source_side(2)=two*bcvalue(2)-PLINE(bfact)
       endif
       if (bctype(1).eq.2) then ! neumann -> reflect even
        local_bctype(1)=0
        local_source_side(1)=PLINE(1)-two*dx_ends*bcvalue(1)
       endif
       if (bctype(2).eq.2) then ! neumann -> reflect even
        local_bctype(2)=0
        local_source_side(2)=PLINE(bfact)+two*dx_ends*bcvalue(2)
       endif

      else if (SEM_IMAGE_BC_ALG.eq.0) then
       ! do nothing
      else
       print *,"SEM_IMAGE_BC_ALG invalid"
       stop
      endif

      AMR_boundary_flag=0

      if ((local_bctype(1).eq.-2).or.(local_bctype(1).eq.-3)) then
       AMR_boundary_flag=1
       local_bctype(1)=0
       if ((x_sep(1).gt.zero).and. &
           (x_sep(1).lt.two)) then
        y_extend(0)=-one-x_sep(1)
        yRT(0)=y_extend(0)
        yRTextrap(0)=y_extend(0)
       else
        print *,"x_sep invalid"
        stop
       endif
      endif

      if ((local_bctype(2).eq.-2).or.(local_bctype(2).eq.-3)) then
       AMR_boundary_flag=1
       local_bctype(2)=0
       if ((x_sep(2).gt.zero).and. &
           (x_sep(2).lt.two)) then
        y_extend(bfact+1)=one+x_sep(2)
        yLT(bfact+1)=y_extend(bfact+1)
        yLTextrap(bfact)=y_extend(bfact+1)
       else
        print *,"x_sep invalid"
        stop
       endif
      endif

      if ((local_bctype(1).eq.0).and. &
          (local_bctype(2).eq.0)) then ! I=left I=right
       PLINE(0)=local_source_side(1)
       PLINE(bfact+1)=local_source_side(2)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        y_extend,yGL_extend)
      else if ((local_bctype(1).eq.0).and. &
               (local_bctype(2).eq.-1)) then !I=left E=right
       PLINE(0)=local_source_side(1)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yRTextrap,yGL_extend)
      else if ((local_bctype(1).eq.1).and. &
               (local_bctype(2).eq.-1)) then !D=left E=right
       PLINE(0)=bcvalue(1)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yRTextrapEXT,yGL_extend)
      else if ((local_bctype(1).eq.2).and. &
               (local_bctype(2).eq.-1)) then !N=left E=right
       sum1=bcvalue(1)*dx_element/two
       do i1=1,bfact
        sum1=sum1-PLINE(i1)*wLT_EXT(i1,0) 
       enddo
       if (wLT_EXT(0,0).eq.zero) then
        print *,"wLT_EXT bust"
        stop
       endif
       PLINE(0)=sum1/wLT_EXT(0,0)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yRTextrapEXT,yGL_extend)
      else if ((local_bctype(1).eq.-1).and. &
               (local_bctype(2).eq.0)) then ! E=left I=right
       do i1=0,bfact-1
        PLINE(i1)=source(i1)
       enddo
       PLINE(bfact)=local_source_side(2)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yLTextrap,yGL_extend)
      else if ((local_bctype(1).eq.-1).and. &
               (local_bctype(2).eq.1)) then ! E=left D=right
       do i1=0,bfact-1
        PLINE(i1)=source(i1)
       enddo
       PLINE(bfact)=bcvalue(2)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yLTextrapEXT,yGL_extend)
      else if ((local_bctype(1).eq.-1).and. &
               (local_bctype(2).eq.2)) then ! E=left N=right
       do i1=0,bfact-1
        PLINE(i1)=source(i1)
       enddo
       sum2=bcvalue(2)*dx_element/two
       do i1=0,bfact-1
        sum2=sum2-PLINE(i1)*wRT_EXT(i1,bfact)
       enddo
       if (wRT_EXT(bfact,bfact).eq.zero) then
        print *,"wRT_EXT bust"
        stop
       endif
       PLINE(bfact)=sum2/wRT_EXT(bfact,bfact)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yLTextrapEXT,yGL_extend)
      else if ((local_bctype(1).eq.-1).and. &
               (local_bctype(2).eq.-1)) then ! E=left E=right
       do i1=0,bfact-1
        PLINE(i1)=source(i1)
       enddo
       call poly_change_basis(bfact-1,bfact+1,PLINE,PLINE2, &
        yLRTextrap,yGL_extend)
      else if ((local_bctype(1).eq.2).and. &
               (local_bctype(2).eq.2)) then ! N=left N=right
       sum1=bcvalue(1)*dx_element/two
       sum2=bcvalue(2)*dx_element/two
       do i1=1,bfact
        sum1=sum1-PLINE(i1)*wLRT(i1,0) 
        sum2=sum2-PLINE(i1)*wLRT(i1,bfact+1)
       enddo
       A=wLRT(0,0)
       B=wLRT(bfact+1,0)
       C=wLRT(0,bfact+1)
       D=wLRT(bfact+1,bfact+1)
       det=A*D-B*C
       if (det.eq.zero) then
        print *,"determinent bust"
        stop
       endif
       PLINE(0)=(D*sum1-B*sum2)/det
       PLINE(bfact+1)=(A*sum2-C*sum1)/det 
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLRT,yGL_extend)
      else if ((local_bctype(1).eq.2).and. &
               (local_bctype(2).eq.0)) then ! N=left I=right
       if (AMR_boundary_flag.ne.0) then
        print *,"this code must be modified for AMR"
        stop
       endif
       PLINE(bfact+1)=local_source_side(2)
       sum1=bcvalue(1)*dx_element/two
       do i1=1,bfact+1
        sum1=sum1-PLINE(i1)*wLT(i1,0) 
       enddo
       if (wLT(0,0).eq.zero) then
        print *,"wLT bust"
        stop
       endif
       PLINE(0)=sum1/wLT(0,0)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLT,yGL_extend)
      else if ((local_bctype(1).eq.2).and. &
               (local_bctype(2).eq.1)) then ! N=left D=right
       PLINE(bfact+1)=bcvalue(2)
       sum1=bcvalue(1)*dx_element/two
       do i1=1,bfact+1
        sum1=sum1-PLINE(i1)*wLRT(i1,0)
       enddo
       if (wLRT(0,0).eq.zero) then
        print *,"wLT bust"
        stop
       endif
       PLINE(0)=sum1/wLRT(0,0)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLRT,yGL_extend)
      else if ((local_bctype(1).eq.1).and. &
               (local_bctype(2).eq.2)) then ! D=left N=right
       PLINE(0)=bcvalue(1)
       sum2=bcvalue(2)*dx_element/two
       do i1=0,bfact
        sum2=sum2-PLINE(i1)*wLRT(i1,bfact+1)
       enddo
       if (wLRT(bfact+1,bfact+1).eq.zero) then
        print *,"wLT bust"
        stop
       endif
       PLINE(bfact+1)=sum2/wLRT(bfact+1,bfact+1)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLRT,yGL_extend)
      else if ((local_bctype(1).eq.0).and. &
               (local_bctype(2).eq.2)) then ! I=left N=right
       if (AMR_boundary_flag.ne.0) then
        print *,"this code must be modified for AMR"
        stop
       endif
       PLINE(0)=local_source_side(1)
       sum2=bcvalue(2)*dx_element/two
       do i1=0,bfact
        sum2=sum2-PLINE(i1)*wRT(i1,bfact+1)
       enddo
       if (wRT(bfact+1,bfact+1).eq.zero) then
        print *,"wRT bust"
        stop
       endif
       PLINE(bfact+1)=sum2/wRT(bfact+1,bfact+1)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yRT,yGL_extend)
      else if ((local_bctype(1).eq.1).and. &
               (local_bctype(2).eq.1)) then ! D=left D=right
       PLINE(0)=bcvalue(1)
       PLINE(bfact+1)=bcvalue(2)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLRT,yGL_extend)
      else if ((local_bctype(1).eq.0).and. &
               (local_bctype(2).eq.1)) then ! I=left D=right
       PLINE(0)=local_source_side(1)
       PLINE(bfact+1)=bcvalue(2)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yRT,yGL_extend)
      else if ((local_bctype(1).eq.1).and. &
               (local_bctype(2).eq.0)) then ! D=left  I=right
       PLINE(0)=bcvalue(1)
       PLINE(bfact+1)=local_source_side(2)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLT,yGL_extend)
      else
       print *,"bctypes invalid"
       stop
      endif

       ! yGL: 0..bfact Gauss Lobatto Legendre (SPTYPE=0)
      call poly_change_basis(bfact+1,bfact,PLINE2,dest_interp,yGL_extend,yGL)
      call deriv_change_basis(bfact+1,bfact,PLINE2,dest_grad, &
        wMAT_extend,yGL_extend,yGL,dx_element)

      if (operation_flag.eq.7) then ! advection

       do i1=0,bfact
        if ((nc.ge.1).and.(nc.le.SDIM)) then
         dest_interp(i1)=dest_interp(i1)*vel(i1)
        else if (nc.eq.SDIM+1) then ! temperature
         dest_interp(i1)=dest_interp(i1)*vel(i1)
        else
         print *,"nc invalid in lineGRAD"
         stop
        endif
        dest_grad(i1)=zero
       enddo ! i1=0..bfact

      else if (operation_flag.eq.0) then ! grad p_MAC
        ! do nothing
      else if (operation_flag.eq.2) then ! grad ppot_MAC/rho_pot_MAC
        ! do nothing
      else if (operation_flag.eq.3) then ! u^{Cell->Mac}
        ! do nothing
      else if (operation_flag.eq.5) then ! u^MAC=u^MAC+(DU)^{cell->mac}
        ! do nothing
      else if (operation_flag.ge.6) then ! rate of strain

       ! do nothing

      else if (operation_flag.eq.11) then

       ! do nothing

      else if (operation_flag.eq.8) then ! coupling terms

       ! do nothing

      else
       print *,"operation_flag invalid3"
       stop
      endif

      do isten=0,bfact

       if (levelrz_in.eq.0) then
        ! do nothing
       else if (levelrz_in.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
       else if (levelrz_in.eq.3) then

        if (operation_flag.eq.6) then

         if ((nc.lt.1).or.(nc.gt.SDIM)) then
          print *,"nc invalid"
          stop
         endif
         if (dir.eq.1) then ! r direction cylindrical coordinates
          if (RRface(isten).le.zero) then
           print *,"RRface invalid"
           stop
          endif
         endif 

        else if (operation_flag.eq.0) then ! grad p_MAC
         ! do nothing
        else if (operation_flag.eq.2) then ! grad ppot_MAC/rho_pot_MAC
         ! do nothing
        else if (operation_flag.eq.3) then ! u^{Cell->Mac}
         ! do nothing
        else if (operation_flag.eq.5) then ! u^MAC=u^MAC+(DU)^{cell->mac}
         ! do nothing
        else if (operation_flag.eq.11) then
         ! do nothing
        else if (operation_flag.eq.7) then
         ! do nothing (advection)
        else if (operation_flag.eq.8) then
         ! do nothing (coupling)
        else
         print *,"operation_flag invalid4"
         stop
        endif

        if (dir.eq.2) then ! theta direction cylindrical coordinates.
         if (RRface(isten).le.zero) then
          print *,"RRface invalid"
          stop
         endif
         dest_grad(isten)=dest_grad(isten)/RRface(isten)
        endif ! dir=2
       else
        print *,"levelrz_in invalid"
        stop
       endif

      enddo ! isten=0..bfact

      return
      end subroutine lineGRAD


       ! mask=0 piecewise FV  mask>0 SEM
      subroutine line_MAC_TO_CELL(source,dest,destdiv,bfact,maskSEM,dx)
      use LegendreNodes

      IMPLICIT NONE

      REAL_T, intent(in) :: dx
      INTEGER_T, intent(in) :: maskSEM
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: source(0:bfact) 
      REAL_T, intent(out) :: dest(0:bfact-1)
      REAL_T, intent(out) :: destdiv(0:bfact-1)

      INTEGER_T i1,j1
      REAL_T y(0:bfact-1)
      REAL_T yGL(0:bfact)
      REAL_T wMATGL(0:bfact,0:bfact)
      REAL_T wtsum,wt
      REAL_T intlo,inthi
      REAL_T dx_element

       ! in: line_MAC_TO_CELL
      if (bfact.lt.1) then
       print *,"bfact invalid35"
       stop
      endif
      dx_element=dx*bfact

      do i1=0,bfact-1
       y(i1)=cache_gauss(bfact,i1,SPTYPE)
      enddo
      do i1=0,bfact
       yGL(i1)=cache_gauss_lobatto(bfact,i1,SPTYPE)
      enddo
      do j1=0,bfact
       do i1=0,bfact
        wMATGL(i1,j1)=cache_wMATGL(bfact,i1,j1,SPTYPE)
       enddo
      enddo
      if ((maskSEM.gt.0).and.(bfact.gt.1)) then
       call poly_change_basis(bfact,bfact-1,source,dest,yGL,y)
       call deriv_change_basis(bfact,bfact-1,source,destdiv, &
        wMATGL,yGL,y,dx_element)
      else if ((maskSEM.eq.0).or.(bfact.eq.1)) then
       do i1=0,bfact-1

        destdiv(i1)=(source(i1+1)-source(i1))/ &
                    (yGL(i1+1)-yGL(i1))
        destdiv(i1)=destdiv(i1)*two/dx_element

        wtsum=zero
        dest(i1)=zero
        do j1=0,bfact
         if (j1.eq.0) then
          intlo=max(zero,yGL(i1))
         else
          intlo=max(y(j1-1),yGL(i1))
         endif
         if (j1.eq.bfact) then
          inthi=min(one,yGL(i1+1))
         else
          inthi=min(y(j1),yGL(i1+1))
         endif
         if (inthi.gt.intlo) then
          wt=inthi-intlo
         else
          wt=zero
         endif
         wtsum=wtsum+wt
         dest(i1)=dest(i1)+wt*source(j1)
        enddo ! j1
        if (wtsum.le.zero) then
         print *,"wtsum invalid"
         stop
        endif
        dest(i1)=dest(i1)/wtsum
       enddo ! i1   
      else
       print *,"maskSEM or bfact invalid"
       stop
      endif

      return
      end subroutine line_MAC_TO_CELL


      subroutine stencilbox( &
       i,j,k,fablo,fabhi,stenlo,stenhi,nsten)
      IMPLICIT NONE

      INTEGER_T i,j,k
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T stenlo(3),stenhi(3)
      INTEGER_T nsten
      INTEGER_T dir2
      INTEGER_T ii(3)

      ii(1)=i
      ii(2)=j
      ii(3)=k
      stenlo(3)=ii(3)
      stenhi(3)=ii(3)

      do dir2=1,SDIM
       stenlo(dir2)=ii(dir2)-nsten
       stenhi(dir2)=ii(dir2)+nsten
      enddo ! dir2
      
      return
      end subroutine stencilbox


      subroutine dimrefine( &
        DIMS(data), &
        DIMS(dataf),dir)
      IMPLICIT NONE

      INTEGER_T dir
      INTEGER_T DIMDEC(data)
      INTEGER_T DIMDEC(dataf)
      INTEGER_T datalo(SDIM),datahi(SDIM)
      INTEGER_T datalof(SDIM),datahif(SDIM)

      call dim_to_box(DIMS(data),datalo,datahi)
      call boxrefine(datalo,datahi,datalof,datahif,dir)
      call box_to_dim(DIMS(dataf),datalof,datahif)

      return
      end subroutine dimrefine

      subroutine get_longdir(lo,hi,longdir,longlo,longhi)
      IMPLICIT NONE

      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T longdir
      INTEGER_T longlo,longhi
      INTEGER_T maxlen,tempmax,dir

      maxlen=0
      longdir=0
      do dir=1,SDIM
       tempmax=hi(dir)-lo(dir)+1
       if (longdir.eq.0) then
        longdir=dir
        maxlen=tempmax
       else if (tempmax.gt.maxlen) then
        longdir=dir 
        maxlen=tempmax
       endif
      enddo
      longlo=lo(longdir)
      longhi=hi(longdir)

      return
      end subroutine get_longdir

       ! nhalf>=3
       ! xsten_center(0,dir)  dir=1..sdim is the coordinate of the center
       !   of the containing cell.
       ! xsten_center(i,dir)  dir=1 is the x coordinate for the 1/2 cell
       !  distances, if uniform, then xsten_center(i,dir)=xcenter+i*dx/2
       ! stencil_offset(dir)= -1,0, or 1.
       ! 
      subroutine bilinear_interp_WT(xsten_center,nhalf,stencil_offset, &
        xtarget,WT)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(in) :: xsten_center(-nhalf:nhalf,SDIM)
      INTEGER_T, intent(in) :: stencil_offset(SDIM)
      REAL_T, intent(in) :: xtarget(SDIM)
      REAL_T, intent(out) :: WT
      INTEGER_T :: WT_ok
      INTEGER_T :: dir
      REAL_T :: denom
      REAL_T theta(SDIM)

      if (nhalf.ge.3) then
       ! do nothing
      else
       print *,"nhalf invalid"
       stop
      endif

      WT_ok=1
      do dir=1,SDIM
       theta(dir)=zero
       if (xtarget(dir).le.xsten_center(0,dir)) then
        if ((stencil_offset(dir).eq.0).or. &
            (stencil_offset(dir).eq.-1)) then
         denom=xsten_center(0,dir)-xsten_center(-2,dir)
         if (denom.gt.zero) then
          theta(dir)=(xsten_center(0,dir)-xtarget(dir))/denom
         else
          print *,"denom invalid"
          stop
         endif
        else if (stencil_offset(dir).eq.1) then
         WT_ok=0
        else
         print *,"stencil_offset(dir) invalid"
         stop
        endif
       else if (xtarget(dir).ge.xsten_center(0,dir)) then
        if ((stencil_offset(dir).eq.0).or. &
            (stencil_offset(dir).eq.1)) then
         denom=xsten_center(0,dir)-xsten_center(2,dir)
         if (denom.lt.zero) then
          theta(dir)=(xsten_center(0,dir)-xtarget(dir))/denom
         else
          print *,"denom invalid"
          stop
         endif
        else if (stencil_offset(dir).eq.-1) then
         WT_ok=0
        else
         print *,"stencil_offset(dir) invalid"
         stop
        endif
       else
        print *,"xtarget or xsten_center invalid"
        stop
       endif
       if (theta(dir).gt.one) then
        theta(dir)=one
       else if ((theta(dir).ge.zero).and. &
                (theta(dir).le.one)) then
        ! do nothing
       else
        print *,"theta(dir) invalid"
        stop
       endif
      enddo ! dir=1..sdim
      if (WT_ok.eq.0) then
       WT=zero
      else if (WT_ok.eq.1) then
       WT=one
       do dir=1,SDIM
        if (stencil_offset(dir).eq.0) then
         WT=WT*(one-theta(dir))
        else if ((stencil_offset(dir).eq.1).or. &
                 (stencil_offset(dir).eq.-1)) then
         WT=WT*theta(dir)
        else
         print *,"stencil_offset(dir) invalid"
         stop
        endif
       enddo !dir=1..sdim
      else
       print *,"WT_ok invalid"
       stop
      endif
     
      return
      end subroutine bilinear_interp_WT 

#define SANITY_CHECK 1

       ! data_in%dir_deriv=1..sdim (derivative)
       ! data_in%dir_deriv=-1 (interp)
      subroutine deriv_from_grid_util(data_in,data_out)
      use probcommon_module
      IMPLICIT NONE

      type(deriv_from_grid_parm_type), intent(in) :: data_in 
      type(interp_from_grid_out_parm_type), intent(out) :: data_out

      INTEGER_T caller_id
      INTEGER_T dir_local
      INTEGER_T ilo(SDIM)
      INTEGER_T ihi(SDIM)
      INTEGER_T istep(SDIM)
      INTEGER_T klosten,khisten,kstep
      INTEGER_T nc
      INTEGER_T data_comp
      INTEGER_T isten,jsten,ksten
      INTEGER_T ii(3)
      REAL_T SGN_FACT
      REAL_T wt_top,wt_bot
      REAL_T xflux_sten(-3:3,SDIM)
      REAL_T xhi_sten(-3:3,SDIM)
      REAL_T xlo_sten(-3:3,SDIM)
      REAL_T dx_sten(SDIM)
      REAL_T dx_top
      REAL_T, target :: xtarget(SDIM)
      INTEGER_T nhalf
      REAL_T, pointer :: local_data_fab(D_DECL(:,:,:),:)
      REAL_T :: local_data_out
      REAL_T, allocatable, dimension(:) :: local_data_max
      REAL_T, allocatable, dimension(:) :: local_data_min
      REAL_T :: scaling
      INTEGER_T :: dummy_input

#ifdef SANITY_CHECK
      type(interp_from_grid_parm_type) :: data_in2 
      type(interp_from_grid_out_parm_type) :: data_out2
#endif

#define dir_FD data_in%dir_deriv
#define ilocal data_in%index_flux

      if ((dir_FD.ge.1).and.(dir_FD.le.SDIM)) then
       ! do nothing (derivative operation)
      else if (dir_FD.eq.-1) then
       ! do nothing (interpolation operation)
      else
       print *,"dir_FD invalid"
       stop
      endif

      local_data_fab=>data_in%disp_data
 
      nhalf=3 
      caller_id=10
      call gridstenMAC_level(xflux_sten,ilocal(1),ilocal(2),ilocal(SDIM), &
        data_in%level,nhalf,data_in%grid_type_flux,caller_id)

      do dir_local=1,SDIM 
       xtarget(dir_local)=xflux_sten(0,dir_local)

       if (data_in%box_type_flux(dir_local).eq. &
           data_in%box_type_data(dir_local)) then
        if (dir_local.eq.dir_FD) then
         ilo(dir_local)=ilocal(dir_local)-1
         ihi(dir_local)=ilocal(dir_local)+1
        else if (dir_local.ne.dir_FD) then
         ilo(dir_local)=ilocal(dir_local)
         ihi(dir_local)=ilocal(dir_local)
        else
         print *,"dir_local bust"
         stop
        endif
       else if ((data_in%box_type_flux(dir_local).eq.0).and. & !CELL
                (data_in%box_type_data(dir_local).eq.1)) then  !NODE
        ilo(dir_local)=ilocal(dir_local)
        ihi(dir_local)=ilocal(dir_local)+1
       else if ((data_in%box_type_flux(dir_local).eq.1).and. & !NODE
                (data_in%box_type_data(dir_local).eq.0)) then  !CELL
        ilo(dir_local)=ilocal(dir_local)-1
        ihi(dir_local)=ilocal(dir_local)
       else
        print *,"box_type corruption"
        stop
       endif
       istep(dir_local)=ihi(dir_local)-ilo(dir_local)
       if (istep(dir_local).eq.0) then
        istep(dir_local)=1
       else if (istep(dir_local).eq.2) then
        ! do nothing
       else if (istep(dir_local).eq.1) then
        ! do nothing
       else
        print *,"istep invalid"
        stop
       endif
      enddo ! dir_local=1..sdim

      call gridstenMAC_level(xhi_sten,ihi(1),ihi(2),ihi(SDIM), &
        data_in%level,nhalf,data_in%grid_type_data,caller_id)
      call gridstenMAC_level(xlo_sten,ilo(1),ilo(2),ilo(SDIM), &
        data_in%level,nhalf,data_in%grid_type_data,caller_id)

      do dir_local=1,SDIM 
       if (ilo(dir_local).eq.ihi(dir_local)) then
        dx_sten(dir_local)=one
       else if (ilo(dir_local).lt.ihi(dir_local)) then
        dx_sten(dir_local)=xhi_sten(0,dir_local)-xlo_sten(0,dir_local)
       else
        print *,"ilo or ihi invalid"
        stop
       endif
      enddo ! dir_local=1..sdim

      if ((data_in%level.ge.0).and. &
          (data_in%level.le.data_in%finest_level).and. &
          (data_in%bfact.ge.1)) then
       ! do nothing
      else
       print *,"level, finest_level or bfact invalid"
       stop
      endif
      do dir_local=1,SDIM
       if (data_in%dx(dir_local).gt.zero) then
        ! do nothing
       else
        print *,"data_in%dx invalid"
        stop
       endif
       if (data_in%fablo(dir_local).ge.0) then
        ! do nothing
       else
        print *,"data_in%fablo invalid"
        stop
       endif
      enddo  ! dir_local=1..sdim
      if (data_in%ncomp.ge.1) then
       ! do nothing
      else
       print *,"ncomp invalid"
       stop
      endif
      if (data_in%scomp.ge.1) then
       ! do nothing
      else
       print *,"scomp invalid"
       stop
      endif

      if (SDIM.eq.2) then
       klosten=0
       khisten=0
       kstep=1
      else if (SDIM.eq.3) then
       klosten=ilo(SDIM)
       khisten=ihi(SDIM)
       kstep=istep(SDIM)
      else
       print *,"sdim invalid"
       stop
      endif

      allocate(local_data_max(data_in%ncomp))
      allocate(local_data_min(data_in%ncomp))

      do nc=1,data_in%ncomp

       local_data_max(nc)=-1.0D-20
       local_data_min(nc)=1.0D-20

       data_out%data_interp(nc)=zero

       do isten=ilo(1),ihi(1),istep(1)
       do jsten=ilo(2),ihi(2),istep(2)
       do ksten=klosten,khisten,kstep
        ii(1)=isten
        ii(2)=jsten
        ii(3)=ksten

        if (dir_FD.eq.-1) then
         SGN_FACT=one
        else if ((dir_FD.ge.1).and.(dir_FD.le.SDIM)) then
         if (ii(dir_FD).eq.ilo(dir_FD)) then
          SGN_FACT=-one
         else if (ii(dir_FD).eq.ihi(dir_FD)) then
          SGN_FACT=one
         else
          print *,"ii(dir_FD) invalid"
          stop
         endif
        else 
         print *,"dir_FD invalid"
         stop
        endif

        wt_top=one
        wt_bot=one
        
        do dir_local=1,SDIM
         if (dir_local.eq.dir_FD) then
          wt_top=wt_top*SGN_FACT
          wt_bot=wt_bot*dx_sten(dir_local)
          if (ihi(dir_local).gt.ilo(dir_local)) then
           ! do nothing
          else
           print *,"ihi or ilo bust"
           stop
          endif
         else if (ihi(dir_local).eq.ilo(dir_local)) then
          ! do nothing
         else if (ihi(dir_local).gt.ilo(dir_local)) then
          if (ii(dir_local).eq.ihi(dir_local)) then
           dx_top=(xtarget(dir_local)-xlo_sten(0,dir_local))
           wt_bot=wt_bot*dx_sten(dir_local)
          else if (ii(dir_local).eq.ilo(dir_local)) then
           dx_top=(xhi_sten(0,dir_local)-xtarget(dir_local))
           wt_bot=wt_bot*dx_sten(dir_local)
          else
           print *,"ii(dir_local) invalid"
           stop
          endif
          if (dx_top.ge.zero) then
           wt_top=wt_top*dx_top
          else
           print *,"dx_top bust"
           stop
          endif
         else
          print *,"ihi or ilo invalid"
          stop
         endif
        enddo ! dir_local=1..sdim

        if ((wt_bot.gt.zero).and.(abs(wt_top).ge.zero)) then 
         data_comp=data_in%scomp+nc-1
         call safe_data(isten,jsten,ksten,data_comp, &
           local_data_fab,local_data_out)

         if (local_data_out.gt.local_data_max(nc)) then
          local_data_max(nc)=local_data_out
         endif
         if (local_data_out.lt.local_data_min(nc)) then
          local_data_min(nc)=local_data_out
         endif

         data_out%data_interp(nc)=data_out%data_interp(nc)+ &
            wt_top*local_data_out/wt_bot
        else
         print *,"wt_bot or wt_top invalid (deriv_from_grid_util):", &
                 wt_bot,wt_top
         print *,"isten,jsten,ksten ",isten,jsten,ksten
         print *,"dir_FD ",dir_FD
         print *,"ilo,ihi ",ilo(1),ilo(2),ilo(SDIM),ihi(1),ihi(2),ihi(SDIM)
         stop
        endif
       enddo !ksten
       enddo !jsten
       enddo !isten
      enddo ! nc=1 .. data_in%ncomp

#ifdef SANITY_CHECK
      if (dir_FD.eq.-1) then
       if (data_in%grid_type_data.eq.-1) then
        data_in2%scomp=data_in%scomp
        data_in2%ncomp=data_in%ncomp
        data_in2%level=data_in%level
        data_in2%finest_level=data_in%finest_level
        data_in2%bfact=data_in%bfact
        data_in2%nmat=num_materials
        data_in2%interp_foot_flag=0
        data_in2%xtarget=>xtarget
        if ((data_in%dx(1).gt.zero).and. &
            (data_in%dx(2).gt.zero).and. &
            (data_in%dx(SDIM).gt.zero)) then
         ! do nothing
        else
         print *,"data_in%dx bust"
         stop
        endif
        data_in2%dx=>data_in%dx
        data_in2%xlo=>data_in%xlo
        data_in2%fablo=>data_in%fablo
        data_in2%fabhi=>data_in%fabhi
        data_in2%state=>data_in%disp_data
        allocate(data_out2%data_interp(data_in2%ncomp))
        call interp_from_grid_util(data_in2,data_out2)
        do nc=1,data_in%ncomp

         scaling=abs(local_data_max(nc))
         if (scaling.lt.abs(local_data_min(nc))) then
          scaling=abs(local_data_min(nc))
         endif        

         if ((scaling.ge.zero).and.(scaling.le.one)) then
          scaling=one
         else if (scaling.ge.one) then
          ! do nothing
         else
          print *,"scaling is NaN"
          stop
         endif

         if (abs(data_out%data_interp(nc)- &
                 data_out2%data_interp(nc)).le.1.0D-12*scaling) then
          ! do nothing
         else
          print *,"data_out%data_interp(nc) invalid"
          print *,"nc=",nc
          print *,"data_in%ncomp=",data_in%ncomp
          print *,"data_in%grid_type_flux=",data_in%grid_type_flux
          print *,"data_in%box_type_flux=",data_in%box_type_flux(1), &
           data_in%box_type_flux(2),data_in%box_type_flux(SDIM)
          print *,"data_in%index_flux=",data_in%index_flux(1), &
           data_in%index_flux(2),data_in%index_flux(SDIM)
          print *,"xtarget=",xtarget(1),xtarget(2),xtarget(SDIM)
          print *,"data_out%data_interp(nc) ",data_out%data_interp(nc)
          print *,"data_out2%data_interp(nc) ",data_out2%data_interp(nc)

          print *,"local_data_max(nc)=",local_data_max(nc)
          print *,"local_data_min(nc)=",local_data_min(nc)
          print *,"scaling=",scaling

          print *,"(breakpoint) break point and gdb: "
          print *,"(1) compile with the -g option"
          print *,"(2) break GLOBALUTIL.F90:8164"
          print *,"By pressing <CTRL C> during this read statement, the"
          print *,"gdb debugger will produce a stacktrace."
          print *,"type 0 then <enter> to exit the program"
          read (*,*) dummy_input
          stop
         endif
        enddo ! nc=1..data_in%ncomp
        deallocate(data_out2%data_interp)
       else if ((data_in%grid_type_data.ge.0).and. &
                (data_in%grid_type_data.le.5)) then
        ! do nothing
       else
        print *,"data_in%grid_type_data invalid"
        stop
       endif
      else if ((dir_FD.ge.1).and.(dir_FD.le.SDIM)) then
       ! do nothing
      else
       print *,"dir_FD invalid"
       stop
      endif
#endif

      deallocate(local_data_max)
      deallocate(local_data_min)

#undef ilocal
#undef dir_FD

      return
      end subroutine deriv_from_grid_util


       ! data_in%dir_deriv=1..sdim (derivative)
       ! data_in%dir_deriv=-1 (interp)
      subroutine single_deriv_from_grid_util(data_in,data_out)
      use probcommon_module
      IMPLICIT NONE

      type(single_deriv_from_grid_parm_type), intent(in) :: data_in 
      type(interp_from_grid_out_parm_type), intent(out) :: data_out

      INTEGER_T caller_id
      INTEGER_T dir_local
      INTEGER_T ilo(SDIM)
      INTEGER_T ihi(SDIM)
      INTEGER_T istep(SDIM)
      INTEGER_T klosten,khisten,kstep
      INTEGER_T isten,jsten,ksten
      INTEGER_T ii(3)
      REAL_T SGN_FACT
      REAL_T wt_top,wt_bot
      REAL_T xflux_sten(-3:3,SDIM)
      REAL_T xhi_sten(-3:3,SDIM)
      REAL_T xlo_sten(-3:3,SDIM)
      REAL_T dx_sten(SDIM)
      REAL_T dx_top
      REAL_T, target :: xtarget(SDIM)
      INTEGER_T nhalf
      REAL_T, pointer :: local_data_fab(D_DECL(:,:,:))
      REAL_T :: local_data_out
      REAL_T :: local_data_max
      REAL_T :: local_data_min
      REAL_T :: scaling
      INTEGER_T :: dummy_input

#ifdef SANITY_CHECK
      type(single_interp_from_grid_parm_type) :: data_in2 
      type(interp_from_grid_out_parm_type) :: data_out2
#endif
  
#define dir_FD data_in%dir_deriv
#define ilocal data_in%index_flux

      if ((dir_FD.ge.1).and.(dir_FD.le.SDIM)) then
       ! do nothing (derivative operation)
      else if (dir_FD.eq.-1) then
       ! do nothing (interpolation operation)
      else
       print *,"dir_FD invalid"
       stop
      endif

      local_data_fab=>data_in%disp_data

      nhalf=3 
      caller_id=10
      call gridstenMAC_level(xflux_sten,ilocal(1),ilocal(2),ilocal(SDIM), &
        data_in%level,nhalf,data_in%grid_type_flux,caller_id)

      do dir_local=1,SDIM 
       xtarget(dir_local)=xflux_sten(0,dir_local)

       if (data_in%box_type_flux(dir_local).eq. &
           data_in%box_type_data(dir_local)) then
        if (dir_local.eq.dir_FD) then
         ilo(dir_local)=ilocal(dir_local)-1
         ihi(dir_local)=ilocal(dir_local)+1
        else if (dir_local.ne.dir_FD) then
         ilo(dir_local)=ilocal(dir_local)
         ihi(dir_local)=ilocal(dir_local)
        else
         print *,"dir_local bust"
         stop
        endif
       else if ((data_in%box_type_flux(dir_local).eq.0).and. & !CELL
                (data_in%box_type_data(dir_local).eq.1)) then  !NODE
        ilo(dir_local)=ilocal(dir_local)
        ihi(dir_local)=ilocal(dir_local)+1
       else if ((data_in%box_type_flux(dir_local).eq.1).and. & !NODE
                (data_in%box_type_data(dir_local).eq.0)) then  !CELL
        ilo(dir_local)=ilocal(dir_local)-1
        ihi(dir_local)=ilocal(dir_local)
       else
        print *,"box_type corruption"
        stop
       endif
       istep(dir_local)=ihi(dir_local)-ilo(dir_local)
       if (istep(dir_local).eq.0) then
        istep(dir_local)=1
       else if (istep(dir_local).eq.2) then
        ! do nothing
       else if (istep(dir_local).eq.1) then
        ! do nothing
       else
        print *,"istep invalid"
        stop
       endif
      enddo ! dir_local=1..sdim

      call gridstenMAC_level(xhi_sten,ihi(1),ihi(2),ihi(SDIM), &
        data_in%level,nhalf,data_in%grid_type_data,caller_id)
      call gridstenMAC_level(xlo_sten,ilo(1),ilo(2),ilo(SDIM), &
        data_in%level,nhalf,data_in%grid_type_data,caller_id)

      do dir_local=1,SDIM 
       if (ilo(dir_local).eq.ihi(dir_local)) then
        dx_sten(dir_local)=one
       else if (ilo(dir_local).lt.ihi(dir_local)) then
        dx_sten(dir_local)=xhi_sten(0,dir_local)-xlo_sten(0,dir_local)
       else
        print *,"ilo or ihi invalid"
        stop
       endif
      enddo ! dir_local=1..sdim

      if ((data_in%level.ge.0).and. &
          (data_in%level.le.data_in%finest_level).and. &
          (data_in%bfact.ge.1)) then
       ! do nothing
      else
       print *,"level, finest_level or bfact invalid"
       stop
      endif
      do dir_local=1,SDIM
       if (data_in%dx(dir_local).gt.zero) then
        ! do nothing
       else
        print *,"data_in%dx invalid"
        stop
       endif
       if (data_in%fablo(dir_local).ge.0) then
        ! do nothing
       else
        print *,"data_in%fablo invalid"
        stop
       endif
      enddo  ! dir_local=1..sdim

      if (SDIM.eq.2) then
       klosten=0
       khisten=0
       kstep=1
      else if (SDIM.eq.3) then
       klosten=ilo(SDIM)
       khisten=ihi(SDIM)
       kstep=istep(SDIM)
      else
       print *,"sdim invalid"
       stop
      endif
      
      local_data_max=-1.0D-20
      local_data_min=1.0D-20

      data_out%data_interp(1)=zero

      do isten=ilo(1),ihi(1),istep(1)
      do jsten=ilo(2),ihi(2),istep(2)
      do ksten=klosten,khisten,kstep
        ii(1)=isten
        ii(2)=jsten
        ii(3)=ksten

        if (dir_FD.eq.-1) then
         SGN_FACT=one
        else if ((dir_FD.ge.1).and.(dir_FD.le.SDIM)) then
         if (ii(dir_FD).eq.ilo(dir_FD)) then
          SGN_FACT=-one
         else if (ii(dir_FD).eq.ihi(dir_FD)) then
          SGN_FACT=one
         else
          print *,"ii(dir_FD) invalid"
          stop
         endif
        else 
         print *,"dir_FD invalid"
         stop
        endif

        wt_top=one
        wt_bot=one
        
        do dir_local=1,SDIM
         if (dir_local.eq.dir_FD) then
          wt_top=wt_top*SGN_FACT
          wt_bot=wt_bot*dx_sten(dir_local)
          if (ihi(dir_local).gt.ilo(dir_local)) then
           ! do nothing
          else
           print *,"ihi or ilo bust"
           stop
          endif
         else if (ihi(dir_local).eq.ilo(dir_local)) then
          ! do nothing
         else if (ihi(dir_local).gt.ilo(dir_local)) then
          if (ii(dir_local).eq.ihi(dir_local)) then
           dx_top=(xtarget(dir_local)-xlo_sten(0,dir_local))
           wt_bot=wt_bot*dx_sten(dir_local)
          else if (ii(dir_local).eq.ilo(dir_local)) then
           dx_top=(xhi_sten(0,dir_local)-xtarget(dir_local))
           wt_bot=wt_bot*dx_sten(dir_local)
          else
           print *,"ii(dir_local) invalid"
           stop
          endif
          if (dx_top.ge.zero) then
           wt_top=wt_top*dx_top
          else
           print *,"dx_top bust"
           stop
          endif
         else
          print *,"ihi or ilo invalid"
          stop
         endif
        enddo ! dir_local=1..sdim

        if ((wt_bot.gt.zero).and.(abs(wt_top).ge.zero)) then 
         call safe_data_single(isten,jsten,ksten,local_data_fab,local_data_out)

         if (local_data_out.gt.local_data_max) then
          local_data_max=local_data_out
         endif
         if (local_data_out.lt.local_data_min) then
          local_data_min=local_data_out
         endif

         if (1.eq.0) then
          print *,"dir_FD,wt_top,wt_bot,isten,jsten,ksten,local_data_out ", &
            dir_FD,wt_top,wt_bot,isten,jsten,ksten,local_data_out
         endif

         data_out%data_interp(1)=data_out%data_interp(1)+ &
             wt_top*local_data_out/wt_bot
        else
         print *,"wt_bot or wt_top invalid(single_deriv_from_grid_util):", &
                 wt_bot,wt_top
         print *,"isten,jsten,ksten ",isten,jsten,ksten
         print *,"dir_FD ",dir_FD
         print *,"ilo,ihi ",ilo(1),ilo(2),ilo(SDIM),ihi(1),ihi(2),ihi(SDIM)
         stop
        endif
      enddo !ksten
      enddo !jsten
      enddo !isten

#ifdef SANITY_CHECK
      if (dir_FD.eq.-1) then
        ! do sanity check if data is at cell centers
       if (data_in%grid_type_data.eq.-1) then 
        data_in2%level=data_in%level
        data_in2%finest_level=data_in%finest_level
        data_in2%bfact=data_in%bfact
        data_in2%interp_foot_flag=0
        data_in2%interp_dir=0 ! not used if interp_foot_flag==0
        data_in2%xtarget=>xtarget
        data_in2%dx=>data_in%dx
        data_in2%xlo=>data_in%xlo
        data_in2%fablo=>data_in%fablo
        data_in2%fabhi=>data_in%fabhi
        data_in2%state=>data_in%disp_data
        allocate(data_out2%data_interp(1))
        call single_interp_from_grid_util(data_in2,data_out2)

        scaling=abs(local_data_max)
        if (scaling.lt.abs(local_data_min)) then
         scaling=abs(local_data_min)
        endif 

        if ((scaling.ge.zero).and.(scaling.le.one)) then
         scaling=one
        else if (scaling.ge.one) then
         ! do nothing
        else
         print *,"scaling is NaN"
         stop
        endif

        if (abs(data_out%data_interp(1)- &
                data_out2%data_interp(1)).le.1.0D-12*scaling) then
         ! do nothing
        else
         print *,"data_in%grid_type_data=",data_in%grid_type_data
         print *,"data_in%grid_type_flux=",data_in%grid_type_flux
         print *,"data_in%box_type_flux=",data_in%box_type_flux(1), &
          data_in%box_type_flux(2),data_in%box_type_flux(SDIM)
         print *,"data_in%index_flux=",data_in%index_flux(1), &
          data_in%index_flux(2),data_in%index_flux(SDIM)
         print *,"data_out%data_interp(1) ",data_out%data_interp(1)
         print *,"data_out2%data_interp(1) ",data_out2%data_interp(1)
         print *,"data_out%data_interp(1) invalid(single_deriv_from_grid_util)"

         print *,"local_data_max=",local_data_max
         print *,"local_data_min=",local_data_min
         print *,"scaling=",scaling

         print *,"(breakpoint) break point and gdb: "
         print *,"(1) compile with the -g option"
         print *,"(2) break GLOBALUTIL.F90:8494"
         print *,"By pressing <CTRL C> during this read statement, the"
         print *,"gdb debugger will produce a stacktrace."
         print *,"type 0 then <enter> to exit the program"
         read (*,*) dummy_input
         stop
        endif
        deallocate(data_out2%data_interp)
       else if ((data_in%grid_type_data.ge.0).and. &
                (data_in%grid_type_data.le.SDIM-1)) then
        allocate(data_out2%data_interp(1))
        data_in2%interp_foot_flag=0

        call interpfab_XDISP( &
          data_in%grid_type_data, & ! start_dir=0..sdim-1
          data_in%grid_type_data, & ! end_dir=0..sdim-1
          data_in2%interp_foot_flag, &
          data_in%bfact, &
          data_in%level, &
          data_in%finest_level, &
          data_in%dx, &
          data_in%xlo, &
          xtarget, &
          data_in%fablo, &
          data_in%fabhi, &
          data_in%disp_data, &
          data_in%disp_data, &
          data_in%disp_data, &
          data_out2%data_interp)

        scaling=abs(local_data_max)
        if (scaling.lt.abs(local_data_min)) then
         scaling=abs(local_data_min)
        endif        

        if ((scaling.ge.zero).and.(scaling.le.one)) then
         scaling=one
        else if (scaling.ge.one) then
         ! do nothing
        else
         print *,"scaling is NaN"
         stop
        endif
        if (abs(data_out%data_interp(1)- &
                data_out2%data_interp(1)).le.1.0D-12*scaling) then
         ! do nothing
        else
         print *,"dir_FD=",dir_FD
         print *,"data_in%grid_type_flux=",data_in%grid_type_flux
         print *,"data_in%grid_type_data=",data_in%grid_type_data
         print *,"data_in%box_type_flux=", &
          data_in%box_type_flux(1), &
          data_in%box_type_flux(2), &
          data_in%box_type_flux(SDIM)
         print *,"data_in%index_flux=",data_in%index_flux(1), &
          data_in%index_flux(2),data_in%index_flux(SDIM)
         print *,"data_in%box_type_data=", &
          data_in%box_type_data(1), &
          data_in%box_type_data(2), &
          data_in%box_type_data(SDIM)
         print *,"data_out%data_interp(1) ",data_out%data_interp(1)
         print *,"data_out2%data_interp(1) ",data_out2%data_interp(1)
         print *,"data_out%data_interp(1) invalid(single_deriv_from_grid_util)"
         print *,"xtarget ",xtarget(1),xtarget(2),xtarget(SDIM)
         print *,"ilocal ",ilocal(1),ilocal(2),ilocal(SDIM)
         print *,"ilo: ",ilo(1),ilo(2),ilo(SDIM)
         print *,"ihi: ",ihi(1),ihi(2),ihi(SDIM)
         print *,"istep: ",istep(1),istep(2),istep(SDIM)
         print *,"klosten,khisten,kstep: ",klosten,khisten,kstep
         print *,"dx_sten ",dx_sten(1),dx_sten(2),dx_sten(SDIM)

         print *,"local_data_max=",local_data_max
         print *,"local_data_min=",local_data_min
         print *,"scaling=",scaling

         print *,"(breakpoint) break point and gdb: "
         print *,"(1) compile with the -g option"
         print *,"(2) break GLOBALUTIL.F90:8571"
         print *,"By pressing <CTRL C> during this read statement, the"
         print *,"gdb debugger will produce a stacktrace."
         print *,"type 0 then <enter> to exit the program"
         read (*,*) dummy_input
         stop
        endif
        deallocate(data_out2%data_interp)
       else if ((data_in%grid_type_data.ge.3).and. &
                (data_in%grid_type_data.le.5)) then
        ! do nothing
       else
        print *,"data_in%grid_type_data invalid"
        stop
       endif
      else if ((dir_FD.ge.1).and.(dir_FD.le.SDIM)) then
       ! do nothing
      else
       print *,"dir_FD invalid"
       stop
      endif
#endif

#undef ilocal
#undef dir_FD

      return
      end subroutine single_deriv_from_grid_util

#undef SANITY_CHECK


       !  This routine calculates the bilinear interpolant at a given point
       !  "x" 
      subroutine interpfab_XDISP( &
       start_dir, &  ! 0<=start_dir<=sdim-1
       end_dir, &    ! 0<=start_dir<=end_dir<=sdim-1
       interp_foot_flag, &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       lo,hi, &
       xdata, &
       ydata, &
       zdata, &
       dest) ! 1..SDIM
      IMPLICIT NONE

      INTEGER_T, intent(in) :: start_dir
      INTEGER_T, intent(in) :: end_dir
      INTEGER_T, intent(in) :: interp_foot_flag
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xtarget(SDIM)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
       ! pointers are always intent(in).
       ! the intent attribute of the data itself is inherited from the
       ! target.
       ! datalox:datahix,dataloy:datahiy,dataloz:datahiz
      REAL_T, intent(in), pointer :: xdata(D_DECL(:,:,:))
      REAL_T, intent(in), pointer :: ydata(D_DECL(:,:,:))
      REAL_T, intent(in), pointer :: zdata(D_DECL(:,:,:))
      REAL_T, intent(out) :: dest(1:end_dir-start_dir+1)

      INTEGER_T dir_disp_comp  ! 0..sdim-1
      INTEGER_T dir_local
      INTEGER_T mac_cell_index(SDIM)
      INTEGER_T istenlo(3),istenhi(3)
      INTEGER_T stencil_offset(SDIM)
      REAL_T WT,total_WT
      INTEGER_T isten,jsten,ksten
      REAL_T local_data

      REAL_T, pointer :: local_data_fab(D_DECL(:,:,:))

      INTEGER_T nhalf
      REAL_T xsten_center(-3:3,SDIM)
      REAL_T xsten_offset(-3:3,SDIM)

      nhalf=3

      if (start_dir.le.end_dir) then
       ! do nothing
      else
       print *,"start_dir or end_dir invalid"
       stop
      endif
      if ((start_dir.ge.0).and.(start_dir.le.SDIM-1)) then
       ! do nothing
      else
       print *,"start_dir invalid"
       stop
      endif
      if ((end_dir.ge.0).and.(end_dir.le.SDIM-1)) then
       ! do nothing
      else
       print *,"end_dir invalid"
       stop
      endif

       ! dir_disp_comp==0 => xdata interpolation
       ! dir_disp_comp==1 => ydata interpolation
       ! dir_disp_comp==2 => zdata interpolation
      do dir_disp_comp=start_dir,end_dir

        ! strategy:
        !   1. determine 3x3x3 MAC grid stencil about x
        !   2. determine the bilinear interpolation weights
        !      (weights are only nonzero in the appropriate 2x2x2 MAC
        !       grid stencil) 
        ! containing_MACcell declared in GLOBALUTIL.F90
       call containing_MACcell(bfact,dx,xlo,lo,xtarget, &
         dir_disp_comp,mac_cell_index)

       istenlo(3)=0
       istenhi(3)=0
       do dir_local=1,SDIM
        istenlo(dir_local)=mac_cell_index(dir_local)-1
        istenhi(dir_local)=mac_cell_index(dir_local)+1
       enddo ! dir_local=1..sdim

       total_WT=zero
       dest(dir_disp_comp-start_dir+1)=zero

       isten=mac_cell_index(1)
       jsten=mac_cell_index(2)
       ksten=mac_cell_index(SDIM)

       call gridstenMAC_level(xsten_center,isten,jsten,ksten,level,nhalf, &
              dir_disp_comp,81)

       if (1.eq.0) then
        print *,"dir_disp_comp,xsten_center,xtarget ", &
         dir_disp_comp,xsten_center(0,1),xsten_center(0,2), &
         xsten_center(0,SDIM),xtarget(1),xtarget(2), &
         xtarget(SDIM)
       endif

       do isten=istenlo(1),istenhi(1)
       do jsten=istenlo(2),istenhi(2)
       do ksten=istenlo(3),istenhi(3)
        call gridstenMAC_level(xsten_offset,isten,jsten,ksten,level,nhalf, &
              dir_disp_comp,81)

        stencil_offset(1)=isten-mac_cell_index(1)
        stencil_offset(2)=jsten-mac_cell_index(2)
        if (SDIM.eq.3) then
         stencil_offset(SDIM)=ksten-mac_cell_index(SDIM)
        endif
        call bilinear_interp_WT(xsten_center,nhalf,stencil_offset,xtarget,WT)
        if ((WT.ge.zero).and.(WT.le.one)) then
         ! do nothing
        else
         print *,"WT invalid"
         stop
        endif

        if (dir_disp_comp.eq.0) then
         local_data_fab=>xdata
        else if (dir_disp_comp.eq.1) then
         local_data_fab=>ydata
        else if ((dir_disp_comp.eq.2).and.(SDIM.eq.3)) then
         local_data_fab=>zdata
        else
         print *,"dir_disp_comp invalid"
         stop
        endif
        call safe_data_single(isten,jsten,ksten,local_data_fab,local_data)

        if ((local_data.ge.-1.0D+30).and. &
            (local_data.le.1.0D+30)) then

         if (interp_foot_flag.eq.0) then
          ! do nothing
         else if (interp_foot_flag.eq.1) then
           ! xdisplace=x-xfoot    xfoot=x-xdisplace
          local_data=xsten_offset(0,dir_disp_comp+1)-local_data
         else
          print *,"interp_foot_flag invalid"
          stop
         endif

         if (1.eq.0) then
          print *,"dir_disp_comp,WT,local_data ", &
            dir_disp_comp,WT,local_data
          print *,"isten,jsten,ksten ",isten,jsten,ksten
         endif

         dest(dir_disp_comp-start_dir+1)= &
            dest(dir_disp_comp-start_dir+1)+WT*local_data
        else
         print *,"local_data overflow"
         print *,"local_data ",local_data
         stop
        endif

        total_WT=total_WT+WT

       enddo ! ksten
       enddo ! jsten
       enddo ! isten

       if (total_WT.gt.zero) then

        dest(dir_disp_comp-start_dir+1)= &
           dest(dir_disp_comp-start_dir+1)/total_WT

       else
        print *,"total_WT invalid"
        stop
       endif

      enddo ! dir_disp_comp=start_dir ... end_dir (0<=dir_disp_comp<=sdim-1)

      return 
      end subroutine interpfab_XDISP


      subroutine interp_from_grid_util(data_in,data_out)
      use probcommon_module
      IMPLICIT NONE
 
      type(interp_from_grid_parm_type), intent(in) :: data_in 
      type(interp_from_grid_out_parm_type), intent(out) :: data_out
      REAL_T :: xsten(-3:3,SDIM)
      REAL_T :: xsten_center(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T dir
      INTEGER_T cell_index(SDIM)
      INTEGER_T stencil_offset(SDIM)
      INTEGER_T istenlo(3)
      INTEGER_T istenhi(3)
      REAL_T WT,total_WT
      INTEGER_T isten,jsten,ksten
      INTEGER_T im
      REAL_T, allocatable, dimension(:) :: local_data
      REAL_T, pointer :: local_data_fab(D_DECL(:,:,:),:)
     
      if ((data_in%level.ge.0).and. &
          (data_in%level.le.data_in%finest_level).and. &
          (data_in%bfact.ge.1)) then
       ! do nothing
      else
       print *,"level, finest_level or bfact invalid"
       stop
      endif
      do dir=1,SDIM
       if (data_in%dx(dir).gt.zero) then
        ! do nothing
       else
        print *,"data_in%dx invalid"
        stop
       endif
       if (data_in%fablo(dir).ge.0) then
        ! do nothing
       else
        print *,"data_in%fablo invalid"
        stop
       endif
      enddo 
      if (data_in%ncomp.ge.1) then
       ! do nothing
      else
       print *,"ncomp invalid"
       stop
      endif
      if (data_in%scomp.ge.1) then
       ! do nothing
      else
       print *,"scomp invalid"
       stop
      endif
      if (data_in%nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"data_in%nmat invalid"
       stop
      endif
      if (data_in%interp_foot_flag.eq.0) then
       ! do nothing
      else if (data_in%interp_foot_flag.eq.1) then
       if (data_in%ncomp.eq.SDIM) then
        ! do nothing
       else
        print *,"ncomp or interp_foot_flag invalid"
        stop
       endif
      else
       print *,"interp_foot_flag invalid"
       stop
      endif
       
      if (data_in%nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif

      allocate(local_data(data_in%ncomp))

      local_data_fab=>data_in%state

      nhalf=3

      call containing_cell( &
        data_in%bfact, &
        data_in%dx, &
        data_in%xlo, &
        data_in%fablo, &
        data_in%xtarget, &
        cell_index)

      istenlo(3)=0
      istenhi(3)=0
      do dir=1,SDIM
       istenlo(dir)=cell_index(dir)-1
       istenhi(dir)=cell_index(dir)+1
      enddo ! dir=1..sdim

      total_WT=zero
      do im=1,data_in%ncomp
       data_out%data_interp(im)=zero
      enddo 

      isten=cell_index(1)
      jsten=cell_index(2)
      ksten=cell_index(SDIM)

      call gridsten_level(xsten_center,isten,jsten,ksten,data_in%level,nhalf)

      do isten=istenlo(1),istenhi(1)
      do jsten=istenlo(2),istenhi(2)
      do ksten=istenlo(3),istenhi(3)

       call gridsten_level(xsten,isten,jsten,ksten,data_in%level,nhalf)
       stencil_offset(1)=isten-cell_index(1)
       stencil_offset(2)=jsten-cell_index(2)
       if (SDIM.eq.3) then
        stencil_offset(SDIM)=ksten-cell_index(SDIM)
       endif
       call bilinear_interp_WT(xsten_center,nhalf,stencil_offset, &
        data_in%xtarget,WT)
       if ((WT.ge.zero).and.(WT.le.one)) then
        ! do nothing
       else
        print *,"WT invalid"
        stop
       endif

       do im=1,data_in%ncomp
        call safe_data(isten,jsten,ksten,data_in%scomp+im-1, &
                local_data_fab,local_data(im))
        if ((local_data(im).ge.-1.0D+30).and. &
            (local_data(im).le.1.0D+30)) then

         if (data_in%interp_foot_flag.eq.0) then
          ! do nothing
         else if (data_in%interp_foot_flag.eq.1) then
           ! xdisplace=x-xfoot    xfoot=x-xdisplace
          local_data(im)=xsten(0,im)-local_data(im)
         else
          print *,"interp_foot_flag invalid"
          stop
         endif

         data_out%data_interp(im)=data_out%data_interp(im)+WT*local_data(im)

        else
         print *,"local_data(im) overflow"
         print *,"im,local_data(im) ",im,local_data(im)
         stop
        endif

       enddo ! im=1..data_in%ncomp

       total_WT=total_WT+WT

      enddo ! ksten
      enddo ! jsten
      enddo ! isten

      if (total_WT.gt.zero) then

       do im=1,data_in%ncomp
        data_out%data_interp(im)=data_out%data_interp(im)/total_WT
       enddo

      else
       print *,"total_WT invalid"
       stop
      endif

      deallocate(local_data)

      return
      end subroutine interp_from_grid_util


      subroutine single_interp_from_grid_util(data_in,data_out)
      use probcommon_module
      IMPLICIT NONE
 
      type(single_interp_from_grid_parm_type), intent(in) :: data_in 
      type(interp_from_grid_out_parm_type), intent(out) :: data_out
      REAL_T :: xsten(-3:3,SDIM)
      REAL_T :: xsten_center(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T dir
      INTEGER_T cell_index(SDIM)
      INTEGER_T stencil_offset(SDIM)
      INTEGER_T istenlo(3)
      INTEGER_T istenhi(3)
      REAL_T WT,total_WT
      INTEGER_T isten,jsten,ksten
      REAL_T :: local_data
      REAL_T, pointer :: local_data_fab(D_DECL(:,:,:))
     
      if ((data_in%level.ge.0).and. &
          (data_in%level.le.data_in%finest_level).and. &
          (data_in%bfact.ge.1)) then
       ! do nothing
      else
       print *,"level, finest_level or bfact invalid"
       stop
      endif
      do dir=1,SDIM
       if (data_in%dx(dir).gt.zero) then
        ! do nothing
       else
        print *,"data_in%dx invalid"
        stop
       endif
       if (data_in%fablo(dir).ge.0) then
        ! do nothing
       else
        print *,"data_in%fablo invalid"
        stop
       endif
      enddo 
      if ((data_in%interp_dir.ge.0).and. &
          (data_in%interp_dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"data_in%interp_dir invalid"
       stop
      endif
      if (data_in%interp_foot_flag.eq.0) then
       ! do nothing
      else if (data_in%interp_foot_flag.eq.1) then
       ! do nothing
      else
       print *,"interp_foot_flag invalid"
       stop
      endif
       
      local_data_fab=>data_in%state

      nhalf=3

      call containing_cell( &
        data_in%bfact, &
        data_in%dx, &
        data_in%xlo, &
        data_in%fablo, &
        data_in%xtarget, &
        cell_index)

      istenlo(3)=0
      istenhi(3)=0
      do dir=1,SDIM
       istenlo(dir)=cell_index(dir)-1
       istenhi(dir)=cell_index(dir)+1
      enddo ! dir=1..sdim

      total_WT=zero
      data_out%data_interp(1)=zero

      isten=cell_index(1)
      jsten=cell_index(2)
      ksten=cell_index(SDIM)

      call gridsten_level(xsten_center,isten,jsten,ksten,data_in%level,nhalf)

      do isten=istenlo(1),istenhi(1)
      do jsten=istenlo(2),istenhi(2)
      do ksten=istenlo(3),istenhi(3)

       call gridsten_level(xsten,isten,jsten,ksten,data_in%level,nhalf)
       stencil_offset(1)=isten-cell_index(1)
       stencil_offset(2)=jsten-cell_index(2)
       if (SDIM.eq.3) then
        stencil_offset(SDIM)=ksten-cell_index(SDIM)
       endif
       call bilinear_interp_WT(xsten_center,nhalf,stencil_offset, &
        data_in%xtarget,WT)
       if ((WT.ge.zero).and.(WT.le.one)) then
        ! do nothing
       else
        print *,"WT invalid"
        stop
       endif

       call safe_data_single(isten,jsten,ksten,local_data_fab,local_data)
       if ((local_data.ge.-1.0D+30).and. &
           (local_data.le.1.0D+30)) then

         if (data_in%interp_foot_flag.eq.0) then
          ! do nothing
         else if (data_in%interp_foot_flag.eq.1) then
           ! xdisplace=x-xfoot    xfoot=x-xdisplace
          local_data=xsten(0,data_in%interp_dir+1)-local_data
         else
          print *,"interp_foot_flag invalid"
          stop
         endif

         data_out%data_interp(1)=data_out%data_interp(1)+WT*local_data

       else
         print *,"local_data overflow"
         print *,"local_data ",local_data
         stop
       endif

       total_WT=total_WT+WT

      enddo ! ksten
      enddo ! jsten
      enddo ! isten

      if (total_WT.gt.zero) then

       data_out%data_interp(1)=data_out%data_interp(1)/total_WT

      else
       print *,"total_WT invalid"
       stop
      endif

      return
      end subroutine single_interp_from_grid_util


      subroutine bilinear_interp_stencil(data_stencil,wt_dist, &
                      ncomp,data_interp,caller_id)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: caller_id
      INTEGER_T, intent(in) :: ncomp
      REAL_T, dimension(D_DECL(2,2,2),ncomp), intent(in) :: data_stencil
      REAL_T, intent(in) :: wt_dist(SDIM)
      REAL_T, intent(out) :: data_interp(ncomp)
      INTEGER_T :: dir
      REAL_T :: c00,c01,c10,c11,c0,c1

      do dir=1,SDIM
       if ((wt_dist(dir).ge.-VOFTOL).and. &
           (wt_dist(dir).le.one+VOFTOL)) then
        ! do nothing
       else
        print *,"wt_dist out of range"
        print *,"dir,wt_dist ",dir,wt_dist(dir)
        print *,"caller_id= ",caller_id
        print *,"ncomp= ",ncomp
        stop
       endif
      enddo ! dir=1..sdim

      do dir=1,ncomp
       c00 = data_stencil(D_DECL(1,1,1),dir)*(one-wt_dist(1)) + &
             data_stencil(D_DECL(2,1,1),dir)*wt_dist(1)
       c01 = data_stencil(D_DECL(1,1,2),dir)*(one-wt_dist(1)) + &
             data_stencil(D_DECL(2,1,2),dir)*wt_dist(1)
       c10 = data_stencil(D_DECL(1,2,1),dir)*(one-wt_dist(1)) + &
             data_stencil(D_DECL(2,2,1),dir)*wt_dist(1)
       c11 = data_stencil(D_DECL(1,2,2),dir)*(one-wt_dist(1)) + &
             data_stencil(D_DECL(2,2,2),dir)*wt_dist(1)

       c0 = c00*(one-wt_dist(2))+c10*wt_dist(2)
       c1 = c01*(one-wt_dist(2))+c11*wt_dist(2)
    
       if (SDIM.eq.3) then
        data_interp(dir) = c0*(one-wt_dist(SDIM))+c1*wt_dist(SDIM)
       else if (SDIM.eq.2) then
        data_interp(dir) = c0
       else
        print *,"macro defined dimension invalid"
        stop
       endif
      enddo ! dir=1..ncomp

      end subroutine bilinear_interp_stencil

      subroutine curverr(curv,LS,xsten,nhalf)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T nhalf
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T LS(D_DECL(-1:1,-1:1,-1:1))
      REAL_T curv
      INTEGER_T inode,jnode,knode
      INTEGER_T ii,jj,kk
      INTEGER_T dir
      INTEGER_T knlo,knhi
      REAL_T mag,delx,rp,rm,rc
      REAL_T nrm(SDIM)
      INTEGER_T interface_found
      REAL_T nn(D_DECL(0:1,0:1,0:1),SDIM)

      if (nhalf.lt.2) then
       print *,"nhalf invalid"
       stop
      endif 
      if (SDIM.eq.3) then
       if ((levelrz.eq.0).or. &
           (levelrz.eq.3)) then
        ! do nothing
       else
        print *,"levelrz invalid"
        stop
       endif
       knlo=0
       knhi=1
      else if (SDIM.eq.2) then
       if ((levelrz.eq.0).or. &
           (levelrz.eq.1).or. &
           (levelrz.eq.3)) then
        ! do nothing
       else
        print *,"levelrz invalid"
        stop
       endif
       knlo=0
       knhi=0
      else
       print *,"dimension bust"
       stop
      endif

      interface_found=0
      do dir=1,SDIM
       ii=0
       jj=0
       kk=0
       if (dir.eq.1) then
        ii=1
       else if (dir.eq.2) then
        jj=1
       else if ((dir.eq.3).and.(SDIM.eq.3)) then
        kk=1
       else
        print *,"dir invalid curverr"
        stop
       endif
       if (LS(D_DECL(0,0,0))*LS(D_DECL(ii,jj,kk)).le.zero) then
        interface_found=1
       endif
       if (LS(D_DECL(0,0,0))*LS(D_DECL(-ii,-jj,-kk)).le.zero) then
        interface_found=1
       endif
      enddo ! dir

      if (1.eq.0) then
       print *,"interface_found=",interface_found
      endif

      if (interface_found.eq.0) then
       curv=zero
      else if (interface_found.eq.1) then
       do inode=0,1
       do jnode=0,1
       do knode=knlo,knhi
        do dir=1,SDIM
         if (dir.eq.1) then
          nrm(dir)= &
           LS(D_DECL(inode,jnode,knode))+ &
           LS(D_DECL(inode,jnode-1,knode))+ &
           LS(D_DECL(inode,jnode,knode-1))+ &
           LS(D_DECL(inode,jnode-1,knode-1))- &
           LS(D_DECL(inode-1,jnode,knode))- &
           LS(D_DECL(inode-1,jnode-1,knode))- &
           LS(D_DECL(inode-1,jnode,knode-1))- &
           LS(D_DECL(inode-1,jnode-1,knode-1))
          delx=xsten(2*inode,dir)-xsten(2*(inode-1),dir)
          if (delx.gt.zero) then
           nrm(dir)=nrm(dir)/delx
          else
           print *,"delx invalid"
           stop
          endif
         else if (dir.eq.2) then
          nrm(dir)= &
           LS(D_DECL(inode,jnode,knode))- &
           LS(D_DECL(inode,jnode-1,knode))+ &
           LS(D_DECL(inode,jnode,knode-1))- &
           LS(D_DECL(inode,jnode-1,knode-1))+ &
           LS(D_DECL(inode-1,jnode,knode))- &
           LS(D_DECL(inode-1,jnode-1,knode))+ &
           LS(D_DECL(inode-1,jnode,knode-1))- &
           LS(D_DECL(inode-1,jnode-1,knode-1))
          delx=xsten(2*jnode,dir)-xsten(2*(jnode-1),dir)
          if (delx.gt.zero) then
           nrm(dir)=nrm(dir)/delx
          else
           print *,"delx invalid"
           stop
          endif
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          nrm(dir)= &
           LS(D_DECL(inode,jnode,knode))+ &
           LS(D_DECL(inode,jnode-1,knode))- &
           LS(D_DECL(inode,jnode,knode-1))- &
           LS(D_DECL(inode,jnode-1,knode-1))+ &
           LS(D_DECL(inode-1,jnode,knode))+ &
           LS(D_DECL(inode-1,jnode-1,knode))- &
           LS(D_DECL(inode-1,jnode,knode-1))- &
           LS(D_DECL(inode-1,jnode-1,knode-1))
          delx=xsten(2*knode,dir)-xsten(2*(knode-1),dir)
          if (delx.gt.zero) then
           nrm(dir)=nrm(dir)/delx
          else
           print *,"delx invalid"
           stop
          endif
         else
          print *,"dir invalid curverr"
          stop
         endif
        enddo ! dir

        if (levelrz.eq.0) then
         rc=one
        else if (levelrz.eq.1) then
         rc=one
        else if (levelrz.eq.3) then
         rc=xsten(0,1)
         if (rc.le.zero) then
          print *,"rc invalid"
          stop
         endif 
        else
         print *,"levelrz invalid"
         stop
        endif 
        call prepare_normal(nrm,rc,mag)
        do dir=1,SDIM
         nn(D_DECL(inode,jnode,knode),dir)=nrm(dir)
        enddo
       enddo
       enddo
       enddo ! inode,jnode,knode

       curv=zero
       do dir=1,SDIM

         ! delx=(x(2)+x(0))/2 - (x(0)+x(-2))/2
        delx=half*(xsten(2,dir)-xsten(-2,dir))
        if (delx.le.zero) then
         print *,"delx invalid"
         stop
        endif
        if (dir.eq.1) then
         if (levelrz.eq.0) then
          rp=one
          rm=one
          rc=one
         else if ((levelrz.eq.1).or. &
                  (levelrz.eq.3)) then
          rp=half*(xsten(2,1)+xsten(0,1))
          rm=half*(xsten(-2,1)+xsten(0,1))
          rc=half*(rp+rm)
          if (levelrz.eq.3) then
           if (rc.le.zero) then
            print *,"rc invalid"
            stop
           endif
          endif
         else
          print *,"levelrz invalid"
          stop
         endif
         curv=curv+( &
          rp*nn(D_DECL(1,1,1),dir)+ &
          rp*nn(D_DECL(1,0,1),dir)+ & 
          rp*nn(D_DECL(1,1,0),dir)+ & 
          rp*nn(D_DECL(1,0,0),dir)- & 
          rm*nn(D_DECL(0,1,1),dir)- &
          rm*nn(D_DECL(0,0,1),dir)- & 
          rm*nn(D_DECL(0,1,0),dir)- & 
          rm*nn(D_DECL(0,0,0),dir)  & 
          )/(rc*delx*four)
        else if (dir.eq.2) then
         if ((levelrz.eq.0).or. &
             (levelrz.eq.1)) then
          rp=one
          rm=one
          rc=one
         else if (levelrz.eq.3) then
          rp=half*(xsten(2,1)+xsten(0,1))
          rm=half*(xsten(-2,1)+xsten(0,1))
          rc=half*(rp+rm)
          if (rc.le.zero) then
           print *,"rc invalid"
           stop
          endif
         else
          print *,"levelrz invalid"
          stop
         endif
         curv=curv+( &
          nn(D_DECL(1,1,1),dir)- &
          nn(D_DECL(1,0,1),dir)+ &
          nn(D_DECL(1,1,0),dir)- &
          nn(D_DECL(1,0,0),dir)+ &
          nn(D_DECL(0,1,1),dir)- &
          nn(D_DECL(0,0,1),dir)+ &
          nn(D_DECL(0,1,0),dir)- &
          nn(D_DECL(0,0,0),dir)  &
          )/(rc*delx*four)
        else if ((dir.eq.3).and.(SDIM.eq.3)) then
         curv=curv+( &
          nn(D_DECL(1,1,1),dir)+ &
          nn(D_DECL(1,0,1),dir)- &
          nn(D_DECL(1,1,0),dir)- &
          nn(D_DECL(1,0,0),dir)+ &
          nn(D_DECL(0,1,1),dir)+ &
          nn(D_DECL(0,0,1),dir)- &
          nn(D_DECL(0,1,0),dir)- &
          nn(D_DECL(0,0,0),dir)  &
          )/(delx*four)
        else
         print *,"dir invalid curverr 2"
         stop
        endif

       enddo ! dir
 
      else
       print *,"interface_found invalid"
       stop
      endif
      if (1.eq.0) then
       print *,"interface_found,curv=",interface_found,curv
      endif

      return
      end subroutine curverr

       ! datatype=0 scalar or vector
       ! datatype=1 tensor face
       ! datatype=2 tensor cell
      subroutine aggressive_worker( &
       datatype, &
       warning_cutoff, &
       tilelo,tilehi, &
       fablo,fabhi, &
       growlo,growhi, &
       bfact, &
       dx, &
       scomp, &
       ncomp, &
       ndefined, &
       ngrow, &
       dir, &
       id, &
       verbose, &
       force_check, &
       gridno,ngrid,level,finest_level, &
       mf, &
       critical_cutoff_low, & ! e.g. -1.0D+99
       critical_cutoff_high)  ! e.g. 1.0D+99

      IMPLICIT NONE

      INTEGER_T, intent(in) :: datatype
      REAL_T, intent(in) :: warning_cutoff
      REAL_T, intent(in) :: critical_cutoff_low
      REAL_T, intent(in) :: critical_cutoff_high
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: growlo(SDIM),growhi(SDIM)
      INTEGER_T growlotest(3),growhitest(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: scomp,ncomp,ndefined
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: dir
      INTEGER_T, intent(in) :: id
      INTEGER_T, intent(in) :: verbose
      INTEGER_T, intent(in) :: force_check
      INTEGER_T, intent(in) :: gridno,ngrid,level,finest_level
      REAL_T, intent(in), pointer :: mf(D_DECL(:,:,:),:)
      INTEGER_T i,j,k,dir2
      INTEGER_T box_type(SDIM)
      INTEGER_T n
      INTEGER_T n_singlelayer,n_interior,n_side,n_corner
      REAL_T sum_interior,sum_side,sum_corner,sum_singlelayer
      REAL_T max_interior,max_side,max_corner,max_singlelayer
      REAL_T val
      INTEGER_T noutside,noutside_single


      if (bfact.lt.1) then
       print *,"bfact invalid36"
       stop
      endif

      if (critical_cutoff_low.lt.critical_cutoff_high) then
       ! do nothing
      else
       print *,"critical_cutoff_low or critical_cutoff_high invalid"
       stop
      endif

      if (((verbose.eq.0).or.(verbose.eq.1)).and. &
          (force_check.eq.0)) then
       ! do nothing
      else if ((verbose.eq.2).or. &
               (force_check.eq.1)) then
       call FLUSH(6) ! unit=6 (screen)

       call grid_type_to_box_type(dir,box_type)

       if (datatype.eq.0) then
        ! do nothing
       else if (datatype.eq.1) then
        ! do nothing
       else if (datatype.eq.2) then
        ! do nothing
       else
        print *,"datatype invalid"
        stop
       endif

       if (scomp+ncomp.gt.ndefined) then
        print *,"scomp invalid aggressive worker"
        print *,"scomp=",scomp
        print *,"ncomp=",ncomp
        print *,"ndefined=",ndefined
        stop
       endif
       if ((gridno.lt.0).or.(gridno.ge.ngrid)) then
        print *,"gridno invalid"
        stop
       endif
       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid aggressive_worker"
        stop
       endif
       if (ngrow.lt.0) then
        print *,"ngrow invalid"
        stop
       endif
       if (warning_cutoff.le.zero) then
        print *,"warning_cutoff invalid"
        stop
       endif

       if ((datatype.eq.0).or.(datatype.eq.1)) then
        call checkbound_array(fablo,fabhi, &
         mf, &
         ngrow,dir,id)
       else if (datatype.eq.2) then
        call checkbound_array(fablo,fabhi, &
         mf, &
         1,-1,id)
       else
        print *,"datatype invalid"
        stop
       endif

       if (verbose.eq.2) then
        print *,"AGGRESSIVE WORKER id= ",id
        print *,"AGGRESSIVE WORKER gridno,ngrid,level,finest_level= ", &
         gridno,ngrid,level,finest_level
        print *,"AGGRESSIVE WORKER scomp,ncomp,ndefined,ngrow=", &
         scomp,ncomp,ndefined,ngrow
        do dir2=1,SDIM
         print *,"dir2,lo,hi,dx ",dir2,fablo(dir2),fabhi(dir2), &
          dx(dir2)
        enddo
       endif

       if (datatype.eq.0) then

        if (dir.eq.-1) then
         call growntilebox(tilelo,tilehi,fablo,fabhi, &
          growlotest,growhitest,ngrow)
        else if ((dir.ge.0).and.(dir.le.5)) then
         call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
          growlotest,growhitest,ngrow,dir,1)
        else
         print *,"dir invalid aggressive worker 2"
         stop
        endif

       else if (datatype.eq.1) then
        if ((ngrow.ne.0).or.(ncomp.ne.1)) then
         print *,"ngrow or ncomp invalid"
         stop
        endif
        call growntileboxTENSOR( &
          tilelo,tilehi,fablo,fabhi, &
          growlotest,growhitest,dir)
       else if (datatype.eq.2) then
        if ((ngrow.ne.0).or.(ncomp.ne.1)) then
         print *,"ngrow or ncomp invalid"
         stop
        endif
        call growntileboxTENSOR_SEM( &
          tilelo,tilehi,fablo,fabhi, &
          growlotest,growhitest,dir)
       else
        print *,"datatype invalid"
        stop
       endif

       do dir2=1,SDIM
        if (growlotest(dir2).ne.growlo(dir2)) then
         print *,"growlo mismatch"
         stop
        endif
        if (growhitest(dir2).ne.growhi(dir2)) then
         print *,"growhi mismatch"
         stop
        endif
       enddo ! dir2


       do n=1,ncomp

        n_singlelayer=0
        n_interior=0
        n_side=0
        n_corner=0
        sum_interior=zero
        sum_side=zero
        sum_corner=zero
        sum_singlelayer=zero
        max_interior=zero
        max_side=zero
        max_corner=zero
        max_singlelayer=zero

        do i=growlotest(1),growhitest(1)
        do j=growlotest(2),growhitest(2)
        do k=growlotest(3),growhitest(3)
         noutside=0
         noutside_single=0
         if ((i.lt.fablo(1)).or.(i.gt.fabhi(1)+box_type(1))) then
          noutside=noutside+1
         endif
         if ((i.lt.fablo(1)-1).or.(i.gt.fabhi(1)+box_type(1)+1)) then
          noutside_single=noutside_single+1
         endif
         if ((j.lt.fablo(2)).or.(j.gt.fabhi(2)+box_type(2))) then
          noutside=noutside+1
         endif
         if ((j.lt.fablo(2)-1).or.(j.gt.fabhi(2)+box_type(2)+1)) then
          noutside_single=noutside_single+1
         endif
         if (SDIM.eq.3) then
          if ((k.lt.fablo(SDIM)).or.(k.gt.fabhi(SDIM)+box_type(SDIM))) then
           noutside=noutside+1
          endif
          if ((k.lt.fablo(SDIM)-1).or.(k.gt.fabhi(SDIM)+box_type(SDIM)+1)) then
           noutside_single=noutside_single+1
          endif
         endif
         val=mf(D_DECL(i,j,k),n+scomp)

         if (val.ge.critical_cutoff_high) then
          print *,"val.ge.critical_cutoff_high ",val,critical_cutoff_high
          print *,"val overflow val,dir,i,j,k,n,scomp,id ", &
           val,dir,i,j,k,n,scomp,id
          print *,"bfact,level,finest_level ",bfact,level,finest_level
          do dir2=1,SDIM
           print *,"dir2,fablo,fabhi ",dir2,fablo(dir2),fabhi(dir2)
           print *,"dir2,tilelo,tilehi ",dir2,tilelo(dir2),tilehi(dir2)
          enddo
          stop
         else if ((val.lt.critical_cutoff_high).and. &
                  (val.gt.critical_cutoff_low)) then
          ! do nothing
         else if (val.le.critical_cutoff_low) then
          print *,"val.le.critical_cutoff_low ",val,critical_cutoff_low
          print *,"val out of bounds val,dir,i,j,k,n,scomp,id ", &
           val,dir,i,j,k,n,scomp,id
          print *,"bfact,level,finest_level ",bfact,level,finest_level
          do dir2=1,SDIM
           print *,"dir2,fablo,fabhi ",dir2,fablo(dir2),fabhi(dir2)
           print *,"dir2,tilelo,tilehi ",dir2,tilelo(dir2),tilehi(dir2)
          enddo
          stop
         else
          print *,"val undefined val,dir,i,j,k,n,scomp,id ", &
           val,dir,i,j,k,n,scomp,id
          print *,"bfact,level,finest_level ",bfact,level,finest_level
          do dir2=1,SDIM
           print *,"dir2,fablo,fabhi ",dir2,fablo(dir2),fabhi(dir2)
           print *,"dir2,tilelo,tilehi ",dir2,tilelo(dir2),tilehi(dir2)
          enddo
          stop
         endif 

         if (noutside.eq.0) then
          sum_interior=sum_interior+val
          n_interior=n_interior+1
          if (abs(val).gt.abs(max_interior)) then
           max_interior=val
          endif 
         endif
         if (noutside_single.eq.0) then
          sum_singlelayer=sum_singlelayer+val
          n_singlelayer=n_singlelayer+1
          if (abs(val).gt.abs(max_singlelayer)) then
           max_singlelayer=val
          endif 
         endif
         if (noutside.eq.1) then
          sum_side=sum_side+val
          n_side=n_side+1
          if (abs(val).gt.abs(max_side)) then
           max_side=val
          endif 
         endif
         if (noutside.gt.1) then
          sum_corner=sum_corner+val
          n_corner=n_corner+1
          if (abs(val).gt.abs(max_corner)) then
           max_corner=val
          endif 
         endif
        enddo
        enddo
        enddo ! i,j,k

        if (n_interior.le.0) then
         print *,"n_interior invalid"
         stop
        endif
        if (n_singlelayer.le.0) then
         print *,"n_singlelayer invalid"
         stop
        endif
        sum_interior=sum_interior/n_interior
        sum_singlelayer=sum_singlelayer/n_singlelayer
        if (n_corner.gt.0) then
         sum_corner=sum_corner/n_corner
        endif
        if (n_side.gt.0) then
         sum_side=sum_side/n_side
        endif

        if (verbose.eq.2) then
         print *,"n+scomp,max_interior,sum_interior,n_interior ", &
          n+scomp,max_interior,sum_interior,n_interior
         print *,"n+scomp,max_singlelayer,sum_singlelayer,n_singlelayer ", &
          n+scomp,max_singlelayer,sum_singlelayer,n_singlelayer
         print *,"n+scomp,max_corner,sum_corner,n_corner ", &
          n+scomp,max_corner,sum_corner,n_corner
         print *,"n+scomp,max_side,sum_side,n_side ", &
          n+scomp,max_side,sum_side,n_side
         if (abs(max_interior).gt.warning_cutoff) then
          print *,"WARNING WARNING MAX_INTERIOR"
         endif
         if (abs(max_singlelayer).gt.warning_cutoff) then
          print *,"WARNING WARNING MAX_SINGLELAYER"
         endif
         if (abs(max_corner).gt.warning_cutoff) then
          print *,"WARNING WARNING MAX_CORNER"
         endif
         if (abs(max_side).gt.warning_cutoff) then
          print *,"WARNING WARNING MAX_SIDE"
         endif
        endif ! verbose==2

       enddo ! n
       call FLUSH(6) ! unit=6
      else
       print *,"verbose invalid"
       stop
      endif

      return
      end subroutine aggressive_worker

      subroutine abs_array_index4(i,j,k,l,Ni,Nj,Nk,Nl,abs_index)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: i,j,k,l,Ni,Nj,Nk,Nl
      INTEGER_T, intent(out) :: abs_index

      if ((i.lt.1).or.(i.gt.Ni).or. &
          (j.lt.1).or.(j.gt.Nj).or. &
          (k.lt.1).or.(k.gt.Nk).or. &
          (l.lt.1).or.(l.gt.Nl)) then
       print *,"index bust abs_array_index4"
       stop
      endif
      abs_index=(i-1)*Nj*Nk*Nl+(j-1)*Nk*Nl+(k-1)*Nl+l

      return
      end subroutine abs_array_index4

       ! finds the cell that contains "x" 
      subroutine containing_cell( &
       bfact,dx,xlo,lo,x,cell_index)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: lo(SDIM)
      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      INTEGER_T, intent(out) :: cell_index(SDIM)

      INTEGER_T lo_e,i1,e_index,dir
      REAL_T xnodes(bfact+1)

       ! NINT=nearest int
      do dir=1,SDIM
       if (bfact.eq.1) then  ! evenly spaced points
        ! x=(i-lo+1/2)dx+xlo  i=(x-xlo)/dx+lo-1/2
        cell_index(dir)=NINT( (x(dir)-xlo(dir))/dx(dir)-half )+lo(dir)
       else if (bfact.gt.1) then
        lo_e=lo(dir)/bfact
        if (lo_e*bfact.ne.lo(dir)) then
         print *,"bfact invalid37"
         stop
        endif
         ! find element whose center is closest to x (i.e. find the
         ! element that contains x)
         ! dx_e=bfact*dx
         ! x=(e_index-lo_e+1/2)dx_e+xlo
        e_index=NINT( (x(dir)-xlo(dir))/(bfact*dx(dir))-half )+lo_e
         ! returns the Gauss Lobatto points in element e_index
        call element_GLnodes1D(xnodes,xlo(dir),e_index,lo_e, &
         dx(dir),bfact)
        if (x(dir).le.xnodes(1)) then
         cell_index(dir)=e_index*bfact
        else if (x(dir).ge.xnodes(bfact+1)) then 
         cell_index(dir)=e_index*bfact+bfact-1
        else
         do i1=1,bfact
          if ((x(dir).ge.xnodes(i1)).and.(x(dir).le.xnodes(i1+1))) then
           cell_index(dir)=e_index*bfact+i1-1
          endif
         enddo
        endif
       else
        print *,"bfact invalid38"
        stop
       endif
      enddo ! dir=1..sdim

      return
      end subroutine containing_cell
   
      ! future work:
      ! turn all the fortran code into an optimized library of routines
      ! that can be called from c++ or python.
      ! NOTE: Antoine Lemoine, JCP, Notus, has developed the fastest 
      !  subroutines for MOF reconstruction.
      !This routine finds the MAC cell that contains "x" 
      ! dir_mac=0,..,sdim-1
      subroutine containing_MACcell( &
       bfact,dx,xlo,lo,x,dir_mac,mac_cell_index)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: dir_mac
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: lo(SDIM)
      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      INTEGER_T, intent(out) :: mac_cell_index(SDIM)

      INTEGER_T lo_e,e_index
      INTEGER_T i1
      INTEGER_T i1crit
      INTEGER_T dir_local
      REAL_T xnodes(bfact+1)

      if ((dir_mac.ge.0).and.(dir_mac.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir_mac invalid"
       stop
      endif

       ! NINT=nearest int
      do dir_local=1,SDIM
       if (bfact.eq.1) then  ! evenly spaced points
        ! dir_local!=dir_mac+1: x=(i-lo+1/2)dx+xlo  i=(x-xlo)/dx+lo-1/2
        ! dir_local==dir_mac+1: x=(i-lo)dx+xlo  i=(x-xlo)/dx+lo
        if (dir_local.ne.dir_mac+1) then
         mac_cell_index(dir_local)= &
          NINT( (x(dir_local)-xlo(dir_local))/dx(dir_local)-half )+ &
          lo(dir_local)
        else if (dir_local.eq.dir_mac+1) then
         mac_cell_index(dir_local)= &
          NINT( (x(dir_local)-xlo(dir_local))/dx(dir_local) )+ &
          lo(dir_local)
        else
         print *,"dir_local or dir_mac bust"
         stop
        endif
       else if (bfact.gt.1) then
        lo_e=lo(dir_local)/bfact
        if (lo_e*bfact.ne.lo(dir_local)) then
         print *,"bfact invalid37"
         stop
        endif
         ! find element whose center is closest to x (i.e. find the
         ! element that contains x)
         ! dx_e=bfact*dx
         ! x=(e_index-lo_e+1/2)dx_e+xlo
        e_index= &
          NINT( (x(dir_local)-xlo(dir_local))/(bfact*dx(dir_local))-half )+lo_e
         ! returns the Gauss Lobatto points in element e_index
        call element_GLnodes1D(xnodes,xlo(dir_local),e_index,lo_e, &
         dx(dir_local),bfact)
       
        if (dir_local.ne.dir_mac+1) then
         if (x(dir_local).le.xnodes(1)) then
          mac_cell_index(dir_local)=e_index*bfact
         else if (x(dir_local).ge.xnodes(bfact+1)) then 
          mac_cell_index(dir_local)=e_index*bfact+bfact-1
         else
          do i1=1,bfact
           if ((x(dir_local).ge.xnodes(i1)).and. &
               (x(dir_local).le.xnodes(i1+1))) then
            mac_cell_index(dir_local)=e_index*bfact+i1-1
           endif
          enddo
         endif
        else if (dir_local.eq.dir_mac+1) then

         i1crit=1
         do i1=2,bfact+1
          if (xnodes(i1).gt.xnodes(i1crit)) then
           if (abs(x(dir_local)-xnodes(i1)).le. &
               abs(x(dir_local)-xnodes(i1crit))) then
            i1crit=i1
           endif
          else
           print *,"expecting xnodes(i1).gt.xnodes(i1crit)"
           stop
          endif
         enddo ! i1=2..bfact+1
         mac_cell_index(dir_local)=e_index*bfact+i1crit-1

        else
         print *,"dir_local or dir_mac bust"
         stop
        endif
       else
        print *,"bfact invalid38"
        stop
       endif
      enddo ! dir_local=1..sdim

      return
      end subroutine containing_MACcell


      subroutine containing_node( &
       bfact,dx,xlo,lo,x,node_index)
      IMPLICIT NONE

      INTEGER_T bfact
      INTEGER_T lo(SDIM)
      REAL_T x(SDIM)
      REAL_T dx(SDIM)
      REAL_T xlo(SDIM)
      INTEGER_T node_index(SDIM)
      INTEGER_T lo_e
      INTEGER_T i1
      INTEGER_T i1crit
      INTEGER_T e_index,dir
      REAL_T xnodes(bfact+1)

       ! NINT=nearest int
      do dir=1,SDIM
       if (bfact.eq.1) then ! evenly spaced points
        ! x=(i-lo)dx+xlo  i=(x-xlo)/dx+lo
        node_index(dir)=NINT( (x(dir)-xlo(dir))/dx(dir) )+lo(dir)
       else if (bfact.gt.1) then
        lo_e=lo(dir)/bfact
        if (lo_e*bfact.ne.lo(dir)) then
         print *,"bfact invalid39"
         stop
        endif
         ! find element whose center is closest to x (i.e. find the
         ! element that contains x)
         ! dx_e=bfact*dx
         ! x=(e_index-lo_e+1/2)dx_e+xlo
        e_index=NINT( (x(dir)-xlo(dir))/(bfact*dx(dir))-half )+lo_e
         ! returns the Gauss Lobatto points in element e_index
         !  e.g. bfact=4
         !  element e_index=0
         !  | .   .    .  .  |
         !  | . x .  x . x . |  (gauss lobatto includes the boundary)
        call element_GLnodes1D(xnodes,xlo(dir),e_index,lo_e, &
         dx(dir),bfact)
        i1crit=1
        do i1=2,bfact+1
         if (xnodes(i1).gt.xnodes(i1crit)) then
          if (abs(x(dir)-xnodes(i1)).le.abs(x(dir)-xnodes(i1crit))) then
           i1crit=i1
          endif
         else
          print *,"expecting xnodes(i1).gt.xnodes(i1crit)"
          stop
         endif
        enddo
        node_index(dir)=e_index*bfact+i1crit-1
       else
        print *,"bfact invalid40"
        stop
       endif
      enddo ! dir=1..sdim

      return
      end subroutine containing_node

      function CTML_FSI_flagF(nmat)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T CTML_FSI_flagF
      INTEGER_T, intent(in) :: nmat
      INTEGER_T im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid CTML_FSI_flagF"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      CTML_FSI_flagF=0
      do im=1,nmat
       if ((FSI_flag(im).eq.4).or. &
           (FSI_flag(im).eq.8)) then
#ifdef MVAHABFSI
        CTML_FSI_flagF=1
#else
        print *,"CTML(F): define MEHDI_VAHAB_FSI in GNUmakefile"
        stop
#endif
       else if ((FSI_flag(im).eq.0).or. &
                (FSI_flag(im).eq.7).or. &
                (FSI_flag(im).eq.1).or. &
                (FSI_flag(im).eq.2).or. &
                (FSI_flag(im).eq.3).or. &
                (FSI_flag(im).eq.6).or. &
                (FSI_flag(im).eq.5)) then
        ! do nothing
       else
        print *,"FSI_flag invalid in CTML_FSI_flagF"
        stop
       endif
      enddo ! im=1..nmat

      return
      end function CTML_FSI_flagF


      function CTML_FSI_mat(nmat,im)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T CTML_FSI_mat
      INTEGER_T, intent(in) :: nmat,im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid CTML_FSI_MAT"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid16"
       stop
      endif
      CTML_FSI_mat=0
      if ((FSI_flag(im).eq.4).or. &
          (FSI_flag(im).eq.8)) then
#ifdef MVAHABFSI
       CTML_FSI_mat=1
#else
       print *,"CTML(F): define MEHDI_VAHAB_FSI in GNUmakefile"
       stop
#endif
      else if ((FSI_flag(im).eq.0).or. &
               (FSI_flag(im).eq.7).or. &
               (FSI_flag(im).eq.1).or. &
               (FSI_flag(im).eq.2).or. &
               (FSI_flag(im).eq.3).or. &
               (FSI_flag(im).eq.6).or. &
               (FSI_flag(im).eq.5)) then
       CTML_FSI_mat=0
      else
       print *,"FSI_flag invalid in CTML_FSI_mat"
       stop
      endif

      return
      end function CTML_FSI_mat

      function fort_FSI_flag_valid(nmat,im)
      use probcommon_module
      IMPLICIT NONE
      INTEGER_T fort_FSI_flag_valid
      INTEGER_T, intent(in) :: nmat,im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid fort_FSI_flag_valid"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid16"
       stop
      endif
      if ((FSI_flag(im).ge.0).and.(FSI_flag(im).le.8)) then
       fort_FSI_flag_valid=1
      else
       print *,"FSI_flag invalid"
       stop
       fort_FSI_flag_valid=0
      endif

      return
      end function fort_FSI_flag_valid

      function is_ice(nmat,im)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T is_ice
      INTEGER_T, intent(in) :: nmat,im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid is_ice"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid16"
       stop
      endif
      is_ice=0
      if ((FSI_flag(im).eq.3).or. &
          (FSI_flag(im).eq.6)) then
       is_ice=1
      else if ((FSI_flag(im).eq.0).or. &
               (FSI_flag(im).eq.7).or. &
               (FSI_flag(im).eq.1).or. &
               (FSI_flag(im).eq.2).or. &
               (FSI_flag(im).eq.4).or. &
               (FSI_flag(im).eq.8).or. &
               (FSI_flag(im).eq.5)) then
       is_ice=0
      else
       print *,"FSI_flag invalid in is_ice"
       stop
      endif

      return
      end function is_ice

      function is_FSI_rigid(nmat,im)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T :: is_FSI_rigid
      INTEGER_T, intent(in) :: nmat,im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid is_FSI_rigid"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid16"
       stop
      endif
! FSI_flag:
! 0 fluid, tessellating (default)
! 1 prescribed rigid solid, non-tessellating (PROB.F90)
! 2 prescribed rigid solid, non-tessellating (sci_clsvof.F90)
! 3 FSI ice,tessellating
! 4 FSI, non-tessellating link w/Kourosh Shoele
! 5 FSI rigid solid, tessellating (PROB.F90)
! 6 FSI ice, tessellating (initial geometry: sci_clsvof.F90)
! 7 fluid, tessellating (initial geometry: sci_clsvof.F90)
! 8 FSI, non-tessellating link w/Kourosh Shoele, pres-vel coupling
      is_FSI_rigid=0
      if (FSI_flag(im).eq.5) then
       is_FSI_rigid=1
      else if ((FSI_flag(im).eq.0).or. &
               (FSI_flag(im).eq.7).or. &
               (FSI_flag(im).eq.1).or. &
               (FSI_flag(im).eq.2).or. &
               (FSI_flag(im).eq.4).or. &
               (FSI_flag(im).eq.8).or. &
               (FSI_flag(im).eq.3).or. &
               (FSI_flag(im).eq.6)) then
       is_FSI_rigid=0
      else
       print *,"FSI_flag invalid in is_FSI_rigid"
       stop
      endif

      return
      end function is_FSI_rigid


      function is_lag_part(nmat,im)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T is_lag_part
      INTEGER_T nmat,im

      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid17 in is_lag_part: im=",im
       print *,"nmat=",nmat
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid is_lag_part"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       error stop
      endif

      if ((FSI_flag(im).eq.1).or. & ! prescribed rigid solid (PROB.F90)
          (FSI_flag(im).eq.2).or. & ! prescribed rigid solid (CAD)
          (FSI_flag(im).eq.4).or. & ! FSI link w/Kourosh Shoele
          (FSI_flag(im).eq.8).or. & ! FSI link w/Kourosh Shoele, pres-vel
          (FSI_flag(im).eq.6).or. & ! lag ice (CAD)
          (FSI_flag(im).eq.7)) then ! lag fluid (CAD)
       is_lag_part=1
      else if ((FSI_flag(im).eq.0).or. &
               (FSI_flag(im).eq.3).or. & ! ice
               (FSI_flag(im).eq.5)) then ! FSI rigid
       is_lag_part=0
      else
       print *,"FSI_flag invalid in is_lag_part"
       stop
      endif

      return
      end function is_lag_part
 
      function is_rigid(nmat,im)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T is_rigid
      INTEGER_T, intent(in) :: nmat,im
      INTEGER_T dummy_input

      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid17 in is_rigid: im=",im
       print *,"nmat=",nmat

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:10214"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read(*,*) dummy_input
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid is_rigid"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      if ((FSI_flag(im).eq.1).or. & ! prescribed rigid solid (PROB.F90)
          (FSI_flag(im).eq.2).or. & ! prescribed rigid solid (sci_clsvof.F90)
          (FSI_flag(im).eq.8).or. & ! FSI pres-vel, Kourosh Shoele
          (FSI_flag(im).eq.4)) then ! FSI link w/Kourosh Shoele
       is_rigid=1  ! non-tessellating material
      else if ((FSI_flag(im).eq.0).or. &
               (FSI_flag(im).eq.7).or. & ! fluid
               (FSI_flag(im).eq.3).or. & ! ice
               (FSI_flag(im).eq.6).or. & ! ice
               (FSI_flag(im).eq.5)) then ! FSI rigid
       is_rigid=0  ! tessellating material
      else
       print *,"FSI_flag invalid in is_rigid"
       stop
      endif

      return
      end function is_rigid

      function fort_is_eulerian_elastic_model(elastic_visc_in, &
        viscoelastic_model_in) 
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T fort_is_eulerian_elastic_model
      REAL_T, intent(in) :: elastic_visc_in
      INTEGER_T, intent(in) :: viscoelastic_model_in

      if (elastic_visc_in.gt.zero) then
       if ((viscoelastic_model_in.eq.0).or. &
           (viscoelastic_model_in.eq.1).or. &
           (viscoelastic_model_in.eq.2).or. &
           (viscoelastic_model_in.eq.3)) then
        fort_is_eulerian_elastic_model=1
       else if (viscoelastic_model_in.eq.4) then
        fort_is_eulerian_elastic_model=0
       else
        print *,"viscoelastic_model_in invalid"
        stop
        fort_is_eulerian_elastic_model=0
       endif
      else if (elastic_visc_in.eq.zero) then
       fort_is_eulerian_elastic_model=0
      else
       print *,"elastic_visc_in invalid"
       stop
       fort_is_eulerian_elastic_model=0
      endif

      return
      end function fort_is_eulerian_elastic_model

      function is_prescribed(nmat,im)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T is_prescribed
      INTEGER_T nmat,im

      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid18"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid is_prescribed"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      if (is_rigid(nmat,im).eq.0) then
       is_prescribed=0
      else if (is_rigid(nmat,im).eq.1) then
       if (CTML_FSI_mat(nmat,im).eq.0) then
        is_prescribed=1
       else if (CTML_FSI_mat(nmat,im).eq.1) then
        if (FSI_flag(im).eq.4) then ! Goldstein et al
         is_prescribed=0
        else if (FSI_flag(im).eq.8) then ! pres-vel coupling
         is_prescribed=1
        else
         print *,"FSI_flag invalid"
         stop
        endif
       else 
        print *,"CTML_FSI_mat(nmat,im) invalid"
        stop
       endif
      else
       print *,"is_rigid invalid"
       stop
      endif

      return
      end function is_prescribed


      function rigid_exists(nmat)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T rigid_exists
      INTEGER_T nmat,im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid rigid_exists"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      rigid_exists=0

      do im=1,nmat
       if (is_rigid(nmat,im).eq.1) then
        rigid_exists=1
       else if (is_rigid(nmat,im).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid in rigid_exists"
        stop
       endif
      enddo ! im

      return
      end function rigid_exists

      function prescribed_exists(nmat)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T prescribed_exists
      INTEGER_T nmat,im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid prescribed_exists"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      prescribed_exists=0

      do im=1,nmat
       if (is_prescribed(nmat,im).eq.1) then
        if (is_rigid(nmat,im).eq.1) then
         prescribed_exists=1
        else
         print *,"is_rigid(nmat,im) invalid"
         stop
        endif
       else if (is_prescribed(nmat,im).eq.0) then
        ! do nothing
       else
        print *,"is_prescribed invalid in prescribed_exists"
        stop
       endif
      enddo ! im=1..nmat

      return
      end function prescribed_exists
 
       ! solid_dist>0 in the solid
      subroutine combine_solid_LS(LS,nmat,solid_dist,im_solid_primary)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T nmat
      REAL_T solid_dist
      REAL_T LS(nmat)
      INTEGER_T im
      INTEGER_T im_solid_primary

      if (nmat.ne.num_materials) then
       print *,"nmat invalid combine_solid_LS"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      solid_dist=-99999.0
      im_solid_primary=0

      do im=1,nmat
       if (is_rigid(nmat,im).eq.0) then
        ! do nothing
       else if (is_rigid(nmat,im).eq.1) then
        if (im_solid_primary.eq.0) then
         solid_dist=LS(im)
         im_solid_primary=im
        else if ((im_solid_primary.ge.1).and.(im_solid_primary.le.nmat)) then
         if (LS(im).gt.solid_dist) then
          solid_dist=LS(im)
          im_solid_primary=im
         else if (LS(im).le.solid_dist) then
          ! do nothing
         else
          print *,"LS(im) bust"
          stop
         endif
        else
         print *,"im_solid_primary invalid"
         stop
        endif 
       else
        print *,"is_rigid invalid"
        stop
       endif
      enddo ! im=1..nmat

      end subroutine combine_solid_LS


      subroutine combine_solid_VOF(VOF,nmat,solid_vof)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T nmat
      REAL_T solid_vof
      REAL_T VOF(nmat)
      INTEGER_T im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid combine_solid_VOF"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      solid_vof=zero

      do im=1,nmat
       if (is_rigid(nmat,im).eq.0) then
        ! do nothing
       else if (is_rigid(nmat,im).eq.1) then
        solid_vof=solid_vof+VOF(im)
       else
        print *,"is_rigid invalid"
        stop
       endif
      enddo ! im=1..nmat

      end subroutine combine_solid_VOF

      subroutine combine_prescribed_VOF(VOF,nmat,solid_vof, &
                      im_prescribed_primary)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T nmat
      REAL_T solid_vof
      REAL_T VOF(nmat)
      INTEGER_T im_prescribed_primary
      INTEGER_T im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid combine_prescribed_VOF"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      solid_vof=zero
      im_prescribed_primary=0

      do im=1,nmat
       if (is_prescribed(nmat,im).eq.0) then
        ! do nothing
       else if (is_prescribed(nmat,im).eq.1) then
        if (is_rigid(nmat,im).eq.1) then
         solid_vof=solid_vof+VOF(im)
         if (im_prescribed_primary.eq.0) then
          im_prescribed_primary=im
         else if ((im_prescribed_primary.ge.1).and. &
                  (im_prescribed_primary.le.nmat)) then
          if (VOF(im).gt.VOF(im_prescribed_primary)) then
           im_prescribed_primary=im
          endif
         else
          print *,"im_prescribed_primary invalid"
          stop
         endif
        else
         print *,"is_rigid(nmat,im) invalid"
         stop
        endif
       else
        print *,"is_prescribed invalid"
        stop
       endif
      enddo ! im=1..nmat

      end subroutine combine_prescribed_VOF

      function rigid_count()
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T rigid_count,im

      rigid_count=0
      do im=1,num_materials
       if (is_rigid(num_materials,im).eq.1) then
        rigid_count=rigid_count+1
       else if (is_rigid(num_materials,im).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif
      enddo ! im

      return
      end function rigid_count

      function is_in_probtype_list()
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T is_in_probtype_list
      INTEGER_T iprob

      is_in_probtype_list=0
      iprob=probtype_list_size
      do while ((is_in_probtype_list.eq.0).and.(iprob.ge.1))
       if (used_probtypes(iprob).eq.probtype) then
        is_in_probtype_list=1
       endif
       iprob=iprob-1
      enddo

      return
      end function is_in_probtype_list

      function im_solid_primary()
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T im_solid_primary
      INTEGER_T im

      im_solid_primary=0
      do im=1,num_materials
       if ((is_rigid(num_materials,im).eq.1).and. &
           (im_solid_primary.eq.0)) then 
        im_solid_primary=im
       else if ((is_rigid(num_materials,im).eq.0).or. &
                (im_solid_primary.gt.0)) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif
      enddo ! im

      return
      end function im_solid_primary

      subroutine debug_im_solid(im)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T im
   
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid19"
       stop
      endif 
      if (is_rigid(num_materials,im).eq.1) then
       ! do nothing
      else if (is_rigid(num_materials,im).eq.0) then
       print *,"expecting solid im=",im
       stop
      else
       print *,"is_rigid invalid"
       stop
      endif

      return 
      end subroutine debug_im_solid

      subroutine get_LS_extend(LS,nmat,iten,LS_extend)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: nmat,iten
      REAL_T, intent(in) :: LS(nmat)
      REAL_T, intent(out) :: LS_extend
      INTEGER_T im,im_opp

      if (nmat.ne.num_materials) then
       print *,"nmat invalid get_LS_extend"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif
      call get_inverse_iten(im,im_opp,iten,nmat)
      if (im.ge.im_opp) then
       print *,"im or im_opp invalid"
       stop
      endif
      if (LS(im).gt.LS(im_opp)) then
       LS_extend=-LS(im_opp)
      else if (LS(im_opp).gt.LS(im)) then
       LS_extend=LS(im)
      else if (LS(im).eq.LS(im_opp)) then
       LS_extend=zero
      else
       print *,"LS bust"
       stop
      endif
      
      return
      end subroutine get_LS_extend

      subroutine get_LSNRM_extend(LS,NRM,nmat,iten,NRM_extend)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nmat,iten
      REAL_T, intent(in) :: LS(nmat)
      REAL_T, intent(in) :: NRM(nmat*SDIM)
      REAL_T, intent(out) :: NRM_extend(SDIM)
      INTEGER_T im,im_opp,dir

      if (nmat.ne.num_materials) then
       print *,"nmat invalid get_LSNRM_extend"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif
      call get_inverse_iten(im,im_opp,iten,nmat)
      if (im.ge.im_opp) then
       print *,"im or im_opp invalid"
       stop
      endif
      if (LS(im).gt.LS(im_opp)) then
       do dir=1,SDIM
        NRM_extend(dir)=-NRM(SDIM*(im_opp-1)+dir)
       enddo
      else if (LS(im_opp).gt.LS(im)) then
       do dir=1,SDIM
        NRM_extend(dir)=NRM(SDIM*(im-1)+dir)
       enddo
      else if (LS(im).eq.LS(im_opp)) then
       do dir=1,SDIM
        NRM_extend(dir)= &
         half*(NRM(SDIM*(im-1)+dir)- &
               NRM(SDIM*(im_opp-1)+dir))
       enddo ! dir=1..sdim
      else
       print *,"LS bust"
       stop
      endif
      
      return
      end subroutine get_LSNRM_extend


      subroutine get_VOF_extend(VOF,nmat,iten,VOF_extend)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: nmat,iten
      REAL_T, intent(in) :: VOF(nmat)
      REAL_T, intent(out) :: VOF_extend
      INTEGER_T im,im_opp

      if (nmat.ne.num_materials) then
       print *,"nmat invalid get_VOF_extend"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif
      call get_inverse_iten(im,im_opp,iten,nmat)
      if (im.ge.im_opp) then
       print *,"im or im_opp invalid"
       stop
      endif
      if ((VOF(im).ge.-VOFTOL).and. &
          (VOF(im).le.one+VOFTOL).and. &
          (VOF(im_opp).ge.-VOFTOL).and. &
          (VOF(im_opp).le.one+VOFTOL)) then 

       if (VOF(im).gt.VOF(im_opp)) then
        VOF_extend=one-VOF(im_opp)
       else if (VOF(im_opp).gt.VOF(im)) then
        VOF_extend=VOF(im)
       else if (abs(VOF(im)-VOF(im_opp)).le.VOFTOL)  then
        VOF_extend=half
       else
        print *,"VOF bust"
        stop
       endif

      else
       print *,"VOF bust"
       stop
      endif
      
      return
      end subroutine get_VOF_extend



        ! "minus" refers to the cell containing the interface with the
        ! largest surface area.
      subroutine gradient_at_dirichlet_boundary( &
         gradT, &
         facearea_minus, &
         facearea_plus, &
         xminus,xplus, &
         Tminus,Tplus, &
         nminus,xminusI, &
         nplus,xplusI, &  ! nplus,xplusI=zero vector if no interface
         x2_I_test, &
         x3_I_test, &
         Tdata) 
      use probcommon_module
      IMPLICIT NONE

      REAL_T gradT(3)
      REAL_T facearea_minus
      REAL_T facearea_plus
      REAL_T facearea_total
      REAL_T xminus(3)
      REAL_T xplus(3)
      REAL_T Tminus,Tplus
      REAL_T nminus(3)
      REAL_T xminusI(3)
      REAL_T nplus(3)
      REAL_T xplusI(3)

      REAL_T mag_I,mag1_I,mag2_I
      REAL_T mag_ntilde,mag1_tilde,mag2_tilde
      REAL_T mag_nplus
      REAL_T mag_test
      INTEGER_T dir,dircrit_tilde,dircrit_I
      INTEGER_T i,j,k
      REAL_T s_I(3)
      REAL_T s_tilde(3)
      REAL_T ntilde(3)
      REAL_T t1_tilde(3)
      REAL_T t2_tilde(3)
      REAL_T t1_I(3)
      REAL_T t2_I(3)
      REAL_T, intent(in) :: x2_I_test(3)
      REAL_T, intent(in) :: x3_I_test(3)
      REAL_T x2_I(3)
      REAL_T x3_I(3)
      REAL_T y(4,3) ! y(j,dir)  j=1..4 dir=1..3
      REAL_T basis_vec(3,3) ! basis_vec(k,dir) k=1..3 dir=1..3
       ! T2I,T3I,TminusI,TplusI
      REAL_T, intent(in) :: Tdata(4) ! Tdata(j) j=1..4
       ! T_{-}+((T_{+}-T_{-})/mag_ntilde)((y_{j}-x_{-}) dot ntilde)
      REAL_T Tstar(4) ! Tstar(j) j=1..4 
      REAL_T wt_epsilon
      REAL_T wt(4)
      REAL_T basis2D(4,4)  ! basis2D(i,j)=g_{i}(y_{j})
      REAL_T AA(2,2) ! AA(k,i)=sum_j w_j g_{i}(y_{j})g_{k}(y_{j})
      REAL_T BB(2)   ! BB(k)=sum_j w_j g_{k}(y_{j})(T_j - T^{*}(y_{j}))
      REAL_T coeffs(2)
      REAL_T fixed_coeff
      INTEGER_T matstatus

      mag_I=zero
      mag_ntilde=zero

      dircrit_I=0
      dircrit_tilde=0
      do dir=1,3
       ntilde(dir)=xplus(dir)-xminus(dir)
       mag_I=mag_I+nminus(dir)**2
       mag_ntilde=mag_ntilde+ntilde(dir)**2

       if (dircrit_I.eq.0) then
        dircrit_I=dir
       else if (abs(nminus(dir)).ge.abs(nminus(dircrit_I))) then
        dircrit_I=dir
       else if (abs(nminus(dir)).le.abs(nminus(dircrit_I))) then
        ! do nothing
       else
        print *,"nminus bust"
        stop
       endif

       if (dircrit_tilde.eq.0) then
        dircrit_tilde=dir
       else if (abs(ntilde(dir)).ge.abs(ntilde(dircrit_tilde))) then
        dircrit_tilde=dir
       else if (abs(ntilde(dir)).le.abs(ntilde(dircrit_tilde))) then
        ! do nothing
       else
        print *,"ntilde bust"
        stop
       endif

      enddo ! dir=1..3

      mag_I=sqrt(mag_I)
      mag_ntilde=sqrt(mag_ntilde)

      if ((abs(mag_I-one).le.VOFTOL).and. &
          (mag_ntilde.gt.zero)) then

       if ((dircrit_I.ge.1).and.(dircrit_I.le.3).and. &
           (dircrit_tilde.ge.1).and.(dircrit_tilde.le.3)) then
        do dir=1,3
         s_I(dir)=one
         s_tilde(dir)=one
         ntilde(dir)=ntilde(dir)/mag_ntilde
        enddo
        s_I(dircrit_I)=zero
        s_tilde(dircrit_tilde)=zero
        call crossprod(nminus,s_I,t1_I)
        call crossprod(ntilde,s_tilde,t1_tilde)
        mag1_I=zero
        mag1_tilde=zero
        do dir=1,3
         mag1_I=mag1_I+t1_I(dir)**2
         mag1_tilde=mag1_tilde+t1_tilde(dir)**2
        enddo
        mag1_I=sqrt(mag1_I)
        mag1_tilde=sqrt(mag1_tilde)
        if ((mag1_I.gt.zero).and.(mag1_tilde.gt.zero)) then

         mag_test=zero
         do dir=1,3
          t1_I(dir)=t1_I(dir)/mag1_I
          t1_tilde(dir)=t1_tilde(dir)/mag1_tilde
          x2_I(dir)=xminusI(dir)+TANGENT_EPS*mag_ntilde*t1_I(dir)
          mag_test=mag_test+(x2_I(dir)-x2_I_test(dir))**2
         enddo
         mag_test=sqrt(mag_test)
         if (mag_test.le.VOFTOL*mag_ntilde) then
          ! do nothing
         else
          print *,"mag_test too big: x2_I"
          stop
         endif

         call crossprod(nminus,t1_I,t2_I)
         call crossprod(ntilde,t1_tilde,t2_tilde)
         mag2_I=zero
         mag2_tilde=zero
         do dir=1,3
          mag2_I=mag2_I+t2_I(dir)**2
          mag2_tilde=mag2_tilde+t2_tilde(dir)**2
         enddo
         mag2_I=sqrt(mag2_I)
         mag2_tilde=sqrt(mag2_tilde)

          ! the magnitude of the cross product of two orthogonal unit
          ! vectors should be one.
         if ((abs(mag2_I-one).le.VOFTOL).and. &
             (abs(mag2_tilde-one).le.VOFTOL)) then

          mag_test=zero
          do dir=1,3
           t2_I(dir)=t2_I(dir)/mag2_I
           t2_tilde(dir)=t2_tilde(dir)/mag2_tilde
           x3_I(dir)=xminusI(dir)+TANGENT_EPS*mag_ntilde*t2_I(dir)
           mag_test=mag_test+(x3_I(dir)-x3_I_test(dir))**2
          enddo
          mag_test=sqrt(mag_test)
          if (mag_test.le.VOFTOL*mag_ntilde) then
           ! do nothing
          else
           print *,"mag_test too big: x3_I"
           stop
          endif

          mag_nplus=zero
          do dir=1,3
           mag_nplus=mag_nplus+nplus(dir)**2
           y(1,dir)=x2_I(dir)
           y(2,dir)=x3_I(dir)
           y(3,dir)=xminusI(dir)
           y(4,dir)=xplusI(dir)
           basis_vec(1,dir)=ntilde(dir)
           basis_vec(2,dir)=t1_tilde(dir)
           basis_vec(3,dir)=t2_tilde(dir)
          enddo
          mag_nplus=sqrt(mag_nplus)
          wt_epsilon=0.01
          wt(1)=wt_epsilon
          wt(2)=wt_epsilon
          facearea_total=facearea_minus+facearea_plus
          if (facearea_total.gt.zero) then
           if (facearea_minus.ge.facearea_plus) then
            wt(3)=facearea_minus/facearea_total
            if (mag_nplus.eq.zero) then
             wt(4)=zero
            else if (mag_nplus.gt.zero) then
             wt(4)=facearea_plus/facearea_total
            else
             print *,"mag_nplus invalid"
             stop
            endif
           else
            print *,"facearea_minus or facearea_plus bad"
            stop
           endif
          else
           print *,"facearea_total cannot be zero"
           stop
          endif

          if (mag_ntilde.gt.zero) then
           fixed_coeff=(Tplus-Tminus)/mag_ntilde
          else
           print *,"mag_ntilde invalid"
           stop
          endif

          do j=1,4
           Tstar(j)=Tminus
           do dir=1,3
            Tstar(j)=Tstar(j)+ &
              fixed_coeff*(y(j,dir)-xminus(dir))*ntilde(dir)
           enddo
          enddo

          do i=1,4 ! 4 basis functions
           do j=1,4 ! 5 data points
            if (wt(j).eq.zero) then
             basis2D(i,j)=zero  !  g_{i}(y_{j})
            else if (wt(j).gt.zero) then
             if (i.eq.1) then
              basis2D(i,j)=one  !  g_{i}(y_{j})
             else if ((i.ge.2).and.(i.le.4)) then
              basis2D(i,j)=zero
              do dir=1,3
               basis2D(i,j)=basis2D(i,j)+ &
                (y(j,dir)-xminus(dir))*basis_vec(i-1,dir)
              enddo
             else
              print *,"i invalid"
              stop
             endif
            else
             print *,"wt invalid"
             stop
            endif
           enddo ! j=1..4
          enddo ! i=1..4

          do k=1,2  ! equation for a_{k}
           do i=1,2  ! basis index
            AA(k,i)=zero
            do j=1,4
             AA(k,i)=AA(k,i)+wt(j)*basis2D(k+2,j)*basis2D(i+2,j)
            enddo
           enddo
           BB(k)=zero
           do j=1,4
            BB(k)=BB(k)+wt(j)*basis2D(k+2,j)*(Tdata(j)-Tstar(j))
           enddo
          enddo !k=1..2

          call matrix_solve(AA,coeffs,BB,matstatus,2)
          if (matstatus.eq.1) then

           do dir=1,3
            gradT(dir)=fixed_coeff*basis_vec(1,dir)
            do k=2,3
             gradT(dir)=gradT(dir)+coeffs(k-1)*basis_vec(k,dir)
            enddo
           enddo

          else
           print *,"matstatus invalid: ",matstatus
           stop
          endif
         else
          print *,"mag2_I or mag2_tilde invalid"
          stop
         endif

        else
         print *,"mag1_I or mag1_tilde invalid"
         stop
        endif
       else
        print *,"dircrit_I or dircrit_tilde invalid"
        stop
       endif
      else
       print *,"mag_I or mag_ntilde invalid"
       stop
      endif

      end subroutine gradient_at_dirichlet_boundary

      subroutine tridiag_solve(l,u,d,n,f,soln)
      IMPLICIT NONE

      INTEGER_T n,i
      REAL_T l(n),u(n),d(n),f(n),soln(n)
      REAL_T ll(n),uu(n),dd(n),z(n)

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
      end subroutine tridiag_solve

      subroutine patterned_substrates(x,y,z,dist,time,im)
      use probcommon_module
      IMPLICIT NONE

      REAL_T, intent(in) :: x,y,z,time
      INTEGER_T, intent(in) :: im
      REAL_T, intent(out) :: dist
      REAL_T :: xprime
      REAL_T :: yprime
      REAL_T :: zprime
      REAL_T :: xvec(SDIM)
      REAL_T :: local_pi
      REAL_T :: pitch,ptb_f,ptb_disbtx,ptb_disbty,ptb_dist,rPillar

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid10"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in patterned dist"
      else if (time.lt.zero) then
       print *,"time invalid in patterned dist"
       stop
      else
       print *,"time bust in patterned dist"
       stop
      endif

      if (is_rigid(num_materials,im).ne.1) then
       print *,"is_rigid invalid"
       stop
      endif

      if ((adv_dir.lt.1).or.(adv_dir.gt.2*SDIM+1)) then
       print *,"adv_dir invalid patterned dist (1)"
       stop
      endif

      xprime=x
      yprime=y
      zprime=z
      xvec(1)=x
      xvec(2)=y
      if (SDIM.eq.3) then
       xvec(SDIM)=z
      endif

      ! ysg patterned surface 
      ! ptb_f=66 !3333.3333
      ! pi=3.14159265
      ! ptb_dist=(sin(pi/2+pi*ptb_f*x)-1.0)**2+(sin(pi/2+pi*ptb_f*y)-1.0)**2
      !  if (ptb_dist.ge.3.95) then
      !          ptb_dist=0.14
      !  else 
      !          ptb_dist=0.0
      !  endif
      !------------------------------------------------------------------------
      ! Ahmed -- this is the first set of run with modified air density to see
      ! what happens with high penetration (more SE to KE conversion)
      ! Do not work on this anymore, copy and paste it to second one to change
      ! and study the pitch , call this study 1 - Ahmed, 11/8/2020
      !------------------------------------------------------------------------
      ! Study 1 - Pillar Penetration with pitch of 0.027
      
      !pi=3.14159265
      !pitch=0.027
      !ptb_f=2/pitch
      ! ptb_disbtx=MOD(x,pitch)
      ! ptb_disbty=MOD(y,pitch)
      !if((ptb_disbtx.le.0.016.or.ptb_disbtx.ge.0.044).and. &
      ! (ptb_disbty.le.0.016.or.ptb_disbty.ge.0.044))then
      ! if(ptb_disbtx.ge.0.044)then
      !  ptb_disbtx=pitch-ptb_disbtx
      ! endif
      ! if(ptb_disbty.ge.0.044)then
      !  ptb_disbty=pitch-ptb_disbty
      ! endif
      ! rPillar=SQRT(ptb_disbtx**2+ptb_disbty**2)
      ! if(rPillar.le.0.016)then
      !  !ptb_dist=SQRT(0.00016**2-rPillar**2)-SQRT(0.00016**2-0.00015**2)+0.002
      !  ptb_dist=0.1
      ! else
      !  ptb_dist=0
      ! endif
      !else
      ! ptb_dist=0
      ! endif
      
      !------------------------------------------------------------------------
      ! Study 2
      !ptb_dist=0.01*sin(20000*pi*x)+(1*10**(-8))*y
      ! this works like a sine function - Ahmed 10-31-2020
      !ptb_dist = 0.03+0.5*0.04*sin((2*pi/0.030)*x) <---- this works
      !ptb_dist = 0.03+0.5*0.04*sin((2*pi/0.031)*x) <---- did not work
      !------------------------------------------------------------------------
      ! Study 3
      ! Implement Paper Guo et al. Droplet Impact on Anisotropic
      ! Superhydrophobic Surfaces, Langmuir
      !ptb_dist=0.05+(8*(2**0.5)/(pi**2))*((0.04*sin(200*x) &
      !         +(0.04*sin(200*3*x)/9)-(0.04*sin(200*5*x)/25) &
      !         -(0.04*sin(200*7*x)/49)+(0.04*sin(200*9*x)/81) &
      !         +(0.04*sin(200*11*x)/121)))    
      !---------------------------------------------------------------
      ! Study 4 - Pillar Penetration with pitch of 0.02 and 0.035 
      !           Pillar radius is held at 0.0125

      local_pi=four*atan(one)

       ! pitch value large yields larger radii posts
      pitch=0.02d0 !0.0350 !0.080 !0.040 !0.030 !0.06 
                   !< this controls the pitch
      ptb_f=two/pitch
      ptb_disbtx=MOD(x,pitch)
      ptb_disbty=MOD(y,pitch)
       ! ! if pitch value cross the limit below, 
       ! ! then the radius of the pillar becomes
       ! ! smaller (less than the .le.) or
       ! ! bigger (if more than .ge.)
      if((ptb_disbtx.le.0.020.or.ptb_disbtx.ge.0.040).and. &
         (ptb_disbty.le.0.020.or.ptb_disbty.ge.0.040))then

        if(ptb_disbtx.ge.0.040)then
         ptb_disbtx=pitch-ptb_disbtx
        endif
        if(ptb_disbty.ge.0.040)then
         ptb_disbty=pitch-ptb_disbty
        endif
        rPillar=SQRT(ptb_disbtx**2+ptb_disbty**2)
        if(rPillar.le.0.01250)then !0.01250 !0.025 !0.025 !0.0125 
                                   !< this controls the radius
         
         !ptb_dist=SQRT(0.00016**2-rPillar**2)-SQRT(0.00016**2-0.00015**2)+ &
         ! 0.002
         ptb_dist=0.05 !< Was 0.0005 to make a dimple surface, now turned it
                         !to 0.001 to make all flat surface
        else
         ptb_dist=0.1
        endif
      else
       ptb_dist=0.1
      endif
       
       !---------------------------------------------------------------------
       !---------------------------------------------------------------------
       ! Study 5 - Pancake Bouncing on Superhydrophobic Surfaces 
       !   Liu,Moevius,Xu,Qian 
       !           Pancake Bouncing: Simulations and Theory and 
       !           Experimental Verification
       !			by Moevius, Liu, Wang, Yeomans
       !pi=3.14159265
       !! pitch value large yields larger radii posts
       !pitch=0.030  !0.080 !0.040 !0.030 !0.06 !< this controls the pitch
       !!ptb_f=2/pitch
       ! ptb_disbtx=MOD(x,pitch)
       ! ptb_disbty=MOD(y,pitch)
       ! ! if pitch value cross the limit below, 
       ! ! then the radius of the pillar becomes
       ! ! smaller (less than the .le.) or
       ! ! bigger (if more than .ge.)
       !if((ptb_disbtx.le.0.010.or.ptb_disbtx.ge.0.0440).and. &
       ! (ptb_disbty.le.0.010.or.ptb_disbty.ge.0.0440))then
       ! if(ptb_disbtx.ge.0.040)then
       !  ptb_disbtx=pitch-ptb_disbtx
       ! endif
       ! if(ptb_disbty.ge.0.040)then
       !  ptb_disbty=pitch-ptb_disbty
       ! endif
       ! rPillar=SQRT(ptb_disbtx**2+ptb_disbty**2)
       ! if(rPillar.le.0.010)then !0.01250 !0.025 !0.025 !0.0125 
                                  !< this controls the radius
         
       !  !ptb_dist=SQRT(0.00016**2-rPillar**2)- &
       !  SQRT(0.00016**2-0.00015**2)+0.002
       !  ptb_dist=0.14
       ! else
       !  ptb_dist=0
       ! endif
       !else
       ! ptb_dist=0
       ! endif
       
       !------------------------------------------------------------------------
      dist=z-ptb_dist
      
      return
      end subroutine patterned_substrates

       ! negative on the inside of the square
      subroutine squaredist(x,y,xlo,xhi,ylo,yhi,dist)
      IMPLICIT NONE

      REAL_T, intent(in) :: x,y,xlo,xhi,ylo,yhi
      REAL_T, intent(out) :: dist
      REAL_T dist1
      REAL_T xmid,ymid
 
      if ((xlo.ge.xhi-1.0D-10).or.(ylo.ge.yhi-1.0D-10)) then 
       print *,"invalid parameters squaredist",xlo,xhi,ylo,yhi
       stop
      endif
      if ((x.le.xlo).and.(y.ge.ylo).and.(y.le.yhi)) then
       dist=xlo-x
      else if ((x.le.xlo).and.(y.ge.yhi)) then
       dist=sqrt( (x-xlo)**2 + (y-yhi)**2 )
      else if ((x.le.xlo).and.(y.le.ylo)) then
       dist=sqrt( (x-xlo)**2 + (y-ylo)**2 )
      else if ((x.ge.xhi).and.(y.ge.ylo).and.(y.le.yhi)) then
       dist=x-xhi
      else if ((x.ge.xhi).and.(y.ge.yhi)) then
       dist=sqrt( (x-xhi)**2 + (y-yhi)**2 )
      else if ((x.ge.xhi).and.(y.le.ylo)) then
       dist=sqrt( (x-xhi)**2 + (y-ylo)**2 )
      else if (y.ge.yhi) then
       dist=y-yhi
      else if (y.le.ylo) then
       dist=ylo-y
      else 
       xmid=half*(xlo+xhi)
       ymid=half*(ylo+yhi)

       if ((x.ge.xlo).and.(x.le.xmid)) then
        dist=x-xlo
       else if ((x.ge.xmid).and.(x.le.xhi)) then
        dist=xhi-x
       else
        print *,"dist invalid in squaredist"
        print *,"x,y,xlo,xhi,ylo,yhi,xmid,ymid ", &
              x,y,xlo,xhi,ylo,yhi,xmid,ymid
         ! gdb: break GLOBALUTIL.F90:9188
        stop
       endif

       if ((y.ge.ylo).and.(y.le.ymid)) then
        dist1=y-ylo
       else if ((y.ge.ymid).and.(y.le.yhi)) then
        dist1=yhi-y
       else
        print *,"dist1 invalid in squaredist"
        stop
       endif
       if (dist.lt.dist1) then
        dist=-dist
       else
        dist=-dist1
       endif
      endif

      return
      end subroutine squaredist

! negative on the inside
      subroutine cubedist(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,dist)
      IMPLICIT NONE

      REAL_T, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
      REAL_T, intent(in) :: x,y,z
      REAL_T, intent(out) :: dist
      REAL_T xcen,ycen,zcen,xrad,yrad,zrad
      REAL_T xdist,ydist,zdist

      if ((xmax.gt.xmin).and. &
          (ymax.gt.ymin).and. &
          (zmax.gt.zmin).and. &
          (abs(x)+abs(y)+abs(z).le.1.0D+20)) then
       ! do nothing
      else
       print *,"xmin,xmax ",xmin,xmax
       print *,"ymin,ymax ",ymin,ymax
       print *,"zmin,zmax ",zmin,zmax
       print *,"cubedist failed"
       stop
      endif

      xcen=half*(xmin+xmax)
      ycen=half*(ymin+ymax)
      zcen=half*(zmin+zmax)
      xrad=xmax-xcen
      yrad=ymax-ycen
      zrad=zmax-zcen

      xdist=abs(x-xcen)-xrad
      ydist=abs(y-ycen)-yrad
      zdist=abs(z-zcen)-zrad

      if ((xdist.le.zero).and.(ydist.le.zero).and.(zdist.le.zero)) then
       dist=xdist
       if (dist.lt.ydist) then
        dist=ydist
       endif
       if (dist.lt.zdist) then
        dist=zdist
       endif
      else
       if (xdist.lt.zero) then
        xdist=zero
       endif
       if (ydist.lt.zero) then
        ydist=zero
       endif
       if (zdist.lt.zero) then
        zdist=zero
       endif
       dist=sqrt(xdist**2+ydist**2+zdist**2)
      endif

      return
      end subroutine cubedist

      subroutine dumpstring_headers_sanity(plot_sdim,ncomp)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: plot_sdim
      INTEGER_T, intent(in) :: ncomp
      character*80 Varname
      character*3 matstr
      INTEGER_T ih,im,i

      if ((plot_sdim.ne.2).and.(plot_sdim.ne.3)) then
       print *,"plot_sdim invalid"
       stop
      endif
      if (ncomp.ge.1) then
       ! do nothing
      else
       print *,"ncomp invalid"
       stop
      endif

      Varname='x_pos'
      call dumpstring(Varname)
      Varname='y_pos'
      call dumpstring(Varname)

      if (plot_sdim.eq.3) then
       Varname='z_pos'
       call dumpstring(Varname)
      endif

      do im=1,ncomp

       write(matstr,'(I3)') im
       do i=1,3
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       Varname='sanity_var'  ! 17..26=>26-17+1=10
       ih=11
       do i=1,3
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)

      enddo ! im=1..ncomp

      return
      end subroutine dumpstring_headers_sanity

      subroutine zones_revolve_sanity( &
       root_char_array, &
       n_root, &
       data_dir, &
       plot_sdim, &
       total_number_grids, &
       grids_per_level_array, &
       levels_array, &
       bfact_array, &
       gridno_array, &
       gridlo_array, &
       gridhi_array, &
       finest_level, &
       SDC_outer_sweeps, &
       slab_step, &
       data_id, &
       nsteps, &
       num_levels, &
       time, &
       visual_revolve, &
       ncomp)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: n_root
      character, dimension(n_root), intent(in) :: root_char_array
      INTEGER_T, intent(in) :: data_dir
      INTEGER_T, intent(in) :: plot_sdim
      INTEGER_T klo_plot,khi_plot

      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: total_number_grids
      INTEGER_T, intent(in) :: num_levels
      INTEGER_T, intent(in) :: grids_per_level_array(num_levels)
      INTEGER_T, intent(in) :: levels_array(total_number_grids)
      INTEGER_T, intent(in) :: bfact_array(total_number_grids)
      INTEGER_T, intent(in) :: gridno_array(total_number_grids)
      INTEGER_T, intent(in) :: gridlo_array(total_number_grids*SDIM)
      INTEGER_T, intent(in) :: gridhi_array(total_number_grids*SDIM)
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: SDC_outer_sweeps
      INTEGER_T, intent(in) :: slab_step
      INTEGER_T, intent(in) :: data_id
      INTEGER_T, intent(in) :: nsteps
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: visual_revolve

      INTEGER_T strandid

      INTEGER_T nwrite3d,nwrite2d,index3d,index2d

      character*3 levstr
      character*5 gridstr
      character*32 filename32
      character*80 rmcommand

      character*6 stepstr
      character*3 outerstr
      character*3 slabstr
      character*3 idstr

      character(len=n_root) :: root_char_str
      character(len=n_root+36) :: newfilename40
      character(len=36) :: fname_extend
      character(len=4) :: step_chars
      character(len=2) :: dir_chars
      character(len=5) :: outer_chars
      character(len=4) :: slab_chars
      character(len=2) :: id_chars
      character(len=4) :: plt_chars

      INTEGER_T i,j,k,dir
      INTEGER_T ilev,igrid
      INTEGER_T lo(plot_sdim),hi(plot_sdim)
      INTEGER_T sysret
      INTEGER_T hi_index_shift(3)

! Guibo

      character*80 Title,Zonename
      REAL*4 ZONEMARKER,EOHMARKER
      integer*4 :: nzones_gb,iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb
      INTEGER_T bfact,testlev,testgridno
      INTEGER_T testlo(SDIM),testhi(SDIM)

      ! define zone structure
      type zone3d_t
         real*8, pointer :: var(:,:,:,:)
      end type zone3d_t
      type(zone3d_t), dimension(:), allocatable :: zone3d_gb

      type zone2d_t
         real*8, pointer :: var(:,:,:)
      end type zone2d_t
      type(zone2d_t), dimension(:), allocatable :: zone2d_gb

      REAL_T theta,rr,zz,xx,yy

      if (plot_sdim.ne.3) then
       print *,"plot_sdim invalid"
       stop
      endif

      if (SDIM.ne.2) then
       print *,"must be called only in 2D"
       stop
      endif

      if (levelrz.ne.1) then
       print *,"levelrz invalid zones revolve"
       stop
      endif
      if (visual_revolve.lt.1) then
       print *,"visual_revolve: ",visual_revolve
       print *,"visual_revolve invalid zones_revolve"
       stop
      endif

      if (ncomp.ge.1) then
       ! do nothing
      else
       print *,"ncomp invalid"
       stop
      endif

      nwrite2d=SDIM+ncomp
      nwrite3d=plot_sdim+ncomp

      if (num_levels.ne.finest_level+1) then
       print *,"num_levels invalid"
       stop
      endif

      write(stepstr,'(I6)') nsteps
      do i=1,6
       if (stepstr(i:i).eq.' ') then
        stepstr(i:i)='0'
       endif
      enddo
      write(outerstr,'(I3)') SDC_outer_sweeps
      write(slabstr,'(I3)') slab_step
      write(idstr,'(I3)') data_id
      do i=1,3
       if (outerstr(i:i).eq.' ') then
        outerstr(i:i)='0'
       endif
       if (slabstr(i:i).eq.' ') then
        slabstr(i:i)='0'
       endif
       if (idstr(i:i).eq.' ') then
        idstr(i:i)='0'
       endif
      enddo

      do i=1,n_root
        root_char_str(i:i)=root_char_array(i)
      enddo

      if (data_dir.eq.-1) then
        dir_chars='CC'
      else if (data_dir.eq.0) then
        dir_chars='XC'
      else if (data_dir.eq.1) then
        dir_chars='YC'
      else if ((data_dir.eq.SDIM-1).and.(SDIM.eq.3)) then
        dir_chars='ZC'
      else if (data_dir.eq.3) then
        dir_chars='XY'
      else if ((data_dir.eq.4).and.(SDIM.eq.3)) then
        dir_chars='XZ'
      else if ((data_dir.eq.5).and.(SDIM.eq.3)) then
        dir_chars='YZ'
      else
        print *,"data_dir invalid"
        stop
      endif
      step_chars='step'
      outer_chars='outer'
      slab_chars='slab'
      id_chars='id'
      plt_chars='.plt'

      newfilename40(1:n_root)=root_char_str
      write(fname_extend,'(A2,A2,A3,A4,A6,A5,A3,A4,A3,A4)') &
               dir_chars,id_chars,idstr, &
               step_chars,stepstr,outer_chars,outerstr, &
               slab_chars,slabstr,plt_chars

      newfilename40(n_root+1:n_root+36)=fname_extend

      print *,"newfilename40 ",newfilename40

      !--------------------------------------------------
      ! Determine nzones_gb and allocate zone_gb, lo_gb, hi_gb
      nzones_gb=0
      do ilev=0,finest_level
      do igrid=0,grids_per_level_array(ilev+1)-1
        nzones_gb = nzones_gb+1
      enddo
      enddo

      if (nzones_gb.ne.total_number_grids) then
       print *,"nzones_gb.ne.total_number_grids: aborting"
       stop
      endif

      print *,"allocating grid structure for tecplot binary"
      print *,"finest_level= ",finest_level
      print *,"number of grids on finest level ", &
         grids_per_level_array(finest_level+1) 
      print *,"total number of grid blocks: ",nzones_gb

      allocate(zone2d_gb(nzones_gb))
      allocate(zone3d_gb(nzones_gb))
      allocate(lo_gb(nzones_gb,plot_sdim))
      allocate(hi_gb(nzones_gb,plot_sdim))

       ! Determine lo_gb, hi_gb
      iz_gb=0
      do ilev=0,finest_level
      do igrid=0,grids_per_level_array(ilev+1)-1
       write(levstr,'(I3)') ilev
       write(gridstr,'(I5)') igrid
       do i=1,3
        if (levstr(i:i).eq.' ') then
         levstr(i:i)='0'
        endif
       enddo
       do i=1,5
        if (gridstr(i:i).eq.' ') then
         gridstr(i:i)='0'
        endif
       enddo
#if (STANDALONE==0)
       write(filename32,'(A14,A10,A3,A5)') &
              './temptecplot/','tempnddata',levstr,gridstr
#elif (STANDALONE==1)
       write(filename32,'(A14,A10,A3,A5)') &
              './temptecplot_','tempnddata',levstr,gridstr
#else
      print *,"STANDALONE invalid"
      stop
#endif
       open(unit=4,file=filename32)

       do dir=1,plot_sdim
        if ((dir.eq.1).or.(dir.eq.2)) then
         read(4,*) lo(dir),hi(dir)
        else if (dir.eq.3) then
         lo(dir)=0
         hi(dir)=visual_revolve-1
        else
         print *,"dir invalid zones revolve"
         stop
        endif
        lo_gb(iz_gb+1,dir)=lo(dir)
        hi_gb(iz_gb+1,dir)=hi(dir)
       enddo ! dir

       bfact=bfact_array(iz_gb+1)
       testlev=levels_array(iz_gb+1)
       testgridno=gridno_array(iz_gb+1)

       if (bfact.lt.1) then
        print *,"bfact invalid150"
        stop
       endif
       if (testlev.ne.ilev) then
        print *,"testlev invalid"
        stop
       endif
       if (testgridno.ne.igrid) then
        print *,"testgridno invalid"
        stop
       endif

       do dir=1,SDIM
        testlo(dir)=gridlo_array(SDIM*iz_gb+dir)
        testhi(dir)=gridhi_array(SDIM*iz_gb+dir)
        if (testlo(dir).ne.lo(dir)) then
         print *,"testlo invalid"
         stop
        endif
        if (testhi(dir).ne.hi(dir)) then
         print *,"testhi invalid"
         stop
        endif
        if ((lo(dir)/bfact)*bfact.ne.lo(dir)) then
         print *,"lo not divisible by bfact"
         stop
        endif
        if (((hi(dir)+1)/bfact)*bfact.ne.hi(dir)+1) then
         print *,"hi+1 not divisible by bfact"
         stop
        endif
       enddo ! dir
 
       close(4)

       iz_gb=iz_gb+1
      enddo ! igrid
      enddo ! ilev

      if (iz_gb.ne.total_number_grids) then
       print *,"iz_gb or total_number_grids invalid"
       stop
      endif

      !-------------------------------------------------------------
      ZONEMARKER = 299.0
      EOHMARKER  = 357.0
      open(unit=11,file=newfilename40,form="unformatted",access="stream")

       ! +++++++ HEADER SECTION ++++++

       ! i.  Magic number, Version number
      write(11) "#!TDV112"

       ! ii. Integer value of 1.
      write(11) 1

       ! iii. Title and variable names.
       ! File type 0 = FULL,1 = GRID,2 = SOLUTION
      write(11) 0
       ! The TITLE
      Title = "DEBUGC data"
      call dumpstring(Title)
       ! Number of variables 
      write(11) nwrite3d

      ! Variable names: zones_revolve_sanity
      call dumpstring_headers_sanity(plot_sdim,ncomp)

       ! Zones
      do iz_gb=1,nzones_gb
       ! Zone marker. Value = 299.0
       write(11) ZONEMARKER
       ! Zone name
       Zonename = "ZONE"
       call dumpstring(Zonename)

       strandid=1    

       write(11) -1   ! Parent Zone
       write(11) 0    ! StrandID
       write(11) time ! Solution time
       write(11) -1   ! Not used. Set to -1
       write(11) 0    ! Zone Type
       write(11) 0    ! Specify Var Location. 0 = Don't specify, 
                      ! all data is located at the nodes.
       write(11) 0    ! Are raw local 1-to-1 face neighbors supplied?
       write(11) 0    ! Number of miscellaneous user-defined  
                      !face neighbor connections

       hi_index_shift(1)=3*(hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1)+1
       hi_index_shift(2)=3*(hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1)+1
       if (plot_sdim.eq.3) then
        hi_index_shift(3)=(hi_gb(iz_gb,plot_sdim)-lo_gb(iz_gb,plot_sdim)+1)+1
       else
        print *,"plot_sdim invalid in zones_revolve_sanity"
        stop
       endif

        ! ----- IMax,JMax,KMax
       write(11) hi_index_shift(1)
       write(11) hi_index_shift(2)
       if (plot_sdim.eq.3) then
        write(11) hi_index_shift(3)
       else
        print *,"plot_sdim invalid in zones_revolve_sanity"
        stop
       endif
 
       write(11) 0
      enddo  ! iz_gb

      write(11) EOHMARKER

      print *,"finished writing header section, now data section"
      print *,"nzones_gb= ",nzones_gb
      print *,"nwrite3d= ",nwrite3d

       ! +++++++ DATA SECTION ++++++

      ilev=0
      igrid=0

      do iz_gb=1,nzones_gb

       if (ilev.gt.finest_level) then
        print *,"ilev invalid"
        print *,"ilev=",ilev
        print *,"finest_level=",finest_level
        stop
       endif
      
       do while (grids_per_level_array(ilev+1).eq.0)
        ilev=ilev+1

        if (igrid.ne.0) then
         print *,"igrid should be 0"
         print *,"igrid=",igrid
         stop
        endif
        if (ilev.gt.finest_level) then
         print *,"ilev invalid in grids_per_level loop"
         print *,"ilev=",ilev
         print *,"finest_level=",finest_level
         stop
        endif
       enddo  ! while grids_per_level_array==0
 
       if (igrid.gt.grids_per_level_array(ilev+1)-1) then
        print *,"igrid invalid"
        print *,"igrid= ",igrid
        print *,"ilev= ",ilev
        print *,"finest_level=",finest_level
        print *,"grids_per_level_array= ", &
         grids_per_level_array(ilev+1)
        print *,"iz_gb = ",iz_gb
        print *,"nzones_gb= ",nzones_gb
        print *,"num_levels= ",num_levels
        print *,"nwrite3d= ",nwrite3d
        stop
       endif

       if (igrid.ne.gridno_array(iz_gb)) then
        print *,"igrid invalid"
        stop
       endif
       if (ilev.ne.levels_array(iz_gb)) then
        print *,"ilev invalid"
        stop
       endif

       do dir=1,plot_sdim
        lo(dir)=lo_gb(iz_gb,dir)
        hi(dir)=hi_gb(iz_gb,dir)
       enddo

       if (plot_sdim.eq.3) then
        klo_plot=lo(plot_sdim)
        khi_plot=hi(plot_sdim)+1
       else if (plot_sdim.eq.2) then
        klo_plot=0
        khi_plot=0
       else
        print *,"plot_sdim invalid"
        stop
       endif

       hi_index_shift(1)=3*(hi(1)-lo(1)+1)
       hi_index_shift(2)=3*(hi(2)-lo(2)+1)
       hi_index_shift(3)=khi_plot-klo_plot

       allocate(zone3d_gb(iz_gb)% &
        var(nwrite3d, &
            0:hi_index_shift(1), &
            0:hi_index_shift(2), &
            0:hi_index_shift(3)))

       allocate(zone2d_gb(iz_gb)% &
        var(nwrite2d, &
            0:hi_index_shift(1), &
            0:hi_index_shift(SDIM)))

       write(levstr,'(I3)') ilev
       write(gridstr,'(I5)') igrid
       do i=1,3
        if (levstr(i:i).eq.' ') then
         levstr(i:i)='0'
        endif
       enddo
       do i=1,5
        if (gridstr(i:i).eq.' ') then
         gridstr(i:i)='0'
        endif
       enddo

#if (STANDALONE==0)
       write(filename32,'(A14,A10,A3,A5)') &
              './temptecplot/','tempnddata',levstr,gridstr
#elif (STANDALONE==1)
       write(filename32,'(A14,A10,A3,A5)') &
              './temptecplot_','tempnddata',levstr,gridstr
#else
       print *,"STANDALONE invalid"
       stop
#endif
       open(unit=4,file=filename32)
       print *,"filename32 ",filename32

       do dir=1,SDIM
        read(4,*) lo(dir),hi(dir)

        if (lo(dir).ne.lo_gb(iz_gb,dir)) then
         print *,"lo and lo_gb different"
         print *,"iz_gb,dir,lo,lo_gb ",iz_gb,dir,lo(dir), &
          lo_gb(iz_gb,dir)
         stop
        endif
        if (hi(dir).ne.hi_gb(iz_gb,dir)) then
         print *,"hi and hi_gb different"
         print *,"iz_gb,dir,hi,hi_gb ",iz_gb,dir,hi(dir), &
          hi_gb(iz_gb,dir)
         stop
        endif

       enddo  ! dir = 1..sdim

       dir=plot_sdim
       lo(dir)=0
       hi(dir)=visual_revolve-1

       if (lo(dir).ne.lo_gb(iz_gb,dir)) then
        print *,"lo and lo_gb different"
        print *,"iz_gb,dir,lo,lo_gb ",iz_gb,dir,lo(dir), &
         lo_gb(iz_gb,dir)
        stop
       endif
       if (hi(dir).ne.hi_gb(iz_gb,dir)) then
        print *,"hi and hi_gb different"
        print *,"iz_gb,dir,hi,hi_gb ",iz_gb,dir,hi(dir), &
         hi_gb(iz_gb,dir)
        stop
       endif

        ! order is IMPORTANT.
       do j=0,hi_index_shift(2)
       do i=0,hi_index_shift(1)
        read(4,*) (zone2d_gb(iz_gb)%var(ivar_gb,i,j),ivar_gb=1,nwrite2d)

        do k=0,hi_index_shift(3)

         index3d=0
         index2d=0

         theta=two*Pi*k/visual_revolve
         rr=zone2d_gb(iz_gb)%var(1,i,j)
         zz=zone2d_gb(iz_gb)%var(2,i,j)
         xx=rr*cos(theta)
         yy=rr*sin(theta)
         zone3d_gb(iz_gb)%var(1,i,j,k)=xx
         zone3d_gb(iz_gb)%var(2,i,j,k)=yy
         zone3d_gb(iz_gb)%var(3,i,j,k)=zz

         index3d=index3d+plot_sdim
         index2d=index2d+SDIM

         do ivar_gb=1,ncomp
          index3d=index3d+1
          index2d=index2d+1
          zone3d_gb(iz_gb)%var(index3d,i,j,k)= &
            zone2d_gb(iz_gb)%var(index2d,i,j)
         enddo

         if ((index3d.ne.nwrite3d).or. &
             (index2d.ne.nwrite2d)) then
          print *,"index mismatch in zone_revolve_sanity"
          stop
         endif

        enddo ! k=0,hi_index_shift(3)

       enddo ! i=0,hi_index_shift(1)
       enddo ! j=0,hi_index_shift(2)

       close(4)

       ! Zone marker  Value = 299.0
       write(11) ZONEMARKER
       ! Data format
       do i=1,nwrite3d
        write(11) 2
       enddo
       write(11) 0  ! Has passive variables: 0 = no, 1 = yes.
       write(11) 0  ! Has variable sharing 0 = no, 1 = yes.
       write(11) -1 ! Share connectivity list (-1 = no sharing). 
   
       do ivar_gb=1,nwrite3d
        write(11) minval(zone3d_gb(iz_gb)%var(ivar_gb,:,:,:))
        write(11) maxval(zone3d_gb(iz_gb)%var(ivar_gb,:,:,:))
       enddo

       ! order is IMPORTANT.
       do ivar_gb=1,nwrite3d
        hi_index_shift(1)=3*(hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1)
        hi_index_shift(2)=3*(hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1)
        hi_index_shift(3)=khi_plot-klo_plot
        do k=0,hi_index_shift(3)
        do j=0,hi_index_shift(2)
        do i=0,hi_index_shift(1)
         write(11) zone3d_gb(iz_gb)%var(ivar_gb,i,j,k)
        enddo
        enddo
        enddo
       enddo ! do ivar_gb=1,nwrite3d

       deallocate(zone3d_gb(iz_gb)%var)
       deallocate(zone2d_gb(iz_gb)%var)

       igrid=igrid+1
       if (igrid.gt.grids_per_level_array(ilev+1)-1) then
        ilev=ilev+1
        igrid=0
       endif

      enddo  ! iz_gb=1,nzones_gp

      deallocate(zone3d_gb)
      deallocate(zone2d_gb)
      deallocate(lo_gb)
      deallocate(hi_gb)

      close(11)
     
      rmcommand='rm ./temptecplot_tempnddata*'
      sysret=0

#if (STANDALONE==0)
      ! do nothing
#elif (STANDALONE==1)
      print *,"issuing command ",rmcommand
      call execute_command_line(rmcommand,exitstat=sysret)

      if (sysret.ne.0) then
       print *,"execute_command_line has sysret=",sysret
       stop
      endif
#else
      print *,"STANDALONE invalid"
      stop
#endif

      return
      end subroutine zones_revolve_sanity


       ! normal is contact line normal pointing towards "im" material
       ! (material 0 liquid)
       ! vel_n is velocity in normal direction.
       ! cos_thetae is the cosine of the static angle inbetween
       ! material 0 and material 2 (the solid material)
       ! cos_thetad is the output cosine of the dynamic angle.
       ! vis is viscosity of material 0 (liquid)
       ! imodel=0 static imodel=1 Jiang   imodel=2 Kistler
      subroutine DCA_select_model(normal,vel_n,cos_thetae,vis, &
        user_tension_scalar, &
        cos_thetad,imodel)
      implicit none

      INTEGER_T imodel
      REAL_T normal(SDIM)
      REAL_T vel_n
      REAL_T cos_thetae  
      REAL_T cos_thetad  
      REAL_T vis 
      REAL_T user_tension_scalar
      REAL_T capillary,f_Hoff_inver,temp,temp1 

      complex(kind=8) :: temp2

      if (user_tension_scalar.le.zero) then
       print *,"user_tension_scalar should be positive"
       stop
      endif

      f_Hoff_inver = 0.0276
      capillary = vel_n*vis/user_tension_scalar

      if (imodel.eq.0) then
       cos_thetad=cos_thetae
      else if (imodel.eq.1) then ! Jiang's model
       temp2 = cmplx(capillary,0.0)
       temp2 = temp2**0.702
       temp = real(temp2)
       temp = tanh(4.96*temp)
       temp = temp*(1.0+cos_thetae)
       temp = cos_thetae-temp
       if (temp.gt.1.0) temp = 1.0  
       if (temp.lt.-1.0) temp = -1.0  
       cos_thetad=temp
      else if (imodel.eq.2) then ! Kistler's model
       temp2 = cmplx(capillary+f_Hoff_inver,0.0)
       temp2 = temp2**0.99
       temp = real(temp2)
       temp1 = temp
       temp = 1.0+1.31*temp
       temp = temp1/temp
       temp2 = cmplx(temp,0.0)
       temp2 = temp2*0.706
       temp = real(temp2)
       temp = 1.0-2.0*tanh(5.16*temp)
       if (temp.gt.1.0) temp = 1.0  
       if (temp.lt.-1.0) temp = -1.0  
       cos_thetad=temp
      else
       print*, 'dynamic contact angle model type not valid'
       stop
      end if

      return
      end subroutine DCA_select_model



       ! Dynamic Contact Angle
      subroutine get_use_DCA(use_DCA)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T use_DCA
      INTEGER_T im_solid_DCA

      im_solid_DCA=im_solid_primary()

      if ((probtype.eq.5501).and.(SDIM.eq.3)) then
       if (im_solid_DCA.ne.3) then
        print *,"im_solid_DCA invalid"
        stop
       endif
       use_DCA=NINT(xblob10)
      else
       use_DCA=-1
       if (fort_ZEYU_DCA_SELECT.eq.-1) then
        ! do nothing
       else if ((fort_ZEYU_DCA_SELECT.ge.1).and. &
                (fort_ZEYU_DCA_SELECT.le.8)) then
        use_DCA=fort_ZEYU_DCA_SELECT+100
       else
        print *,"fort_ZEYU_DCA_SELECT invalid"
        stop
       endif 
      endif

      return
      end subroutine get_use_DCA

      subroutine VISC_dodecane(rho_in,T_in,visccoef,eta)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho_in,T_in
      ! From Caudwell et al., Int. J. Thermophysics (25) 5, 2004
      REAL_T rho,T,eta,visccoef
      REAL_T rmin,rmax
      REAL_T T2,V0,V,VV,VV2,eta_star

      rho=rho_in
      T=T_in

      if (T.le.zero) then
       print *,"VISC_dodecane: temperature zero or negative ",T
       stop
      endif
      T = max(T,298.15)
      T = min(T,473.15)

      T2 = T*T
       ! molar core volume
      V0 = 191.54e-6-0.441338e-6*T+8.98744e-10*T2-6.7792e-13*T2*T 

      if (rho.le.zero) then
       print *,"VISC_dodecane: density zero or negative ",rho
       stop
      endif
      rmin = 0.6089 ! at 0.1 MPa, T = 473.15
      rmax = 0.8091 ! at 161.33 MPa, T = 298.15
      rho = min(rho,rmax)
      rho = max(rho,rmin)

      V  = 0.17034e-3/rho ! in m3/mol
      VV = V/V0
      VV2 = VV*VV
      eta_star = 1/(0.321621-0.4803715*VV+0.222206*VV2-2.964626e-2*VV2*VV)
      eta = 1.9720e-8*V**(-0.666667)*T**0.5*eta_star  ! in g/cm-s

      if (eta.le.zero) then
       print *,"VISC_dodecane: negative or zero viscosity ", eta,&
        "  from ",rho,T
       stop
      endif

      eta = min(eta,0.03677)
      eta = max(eta,0.00218)

      if (OLD_DODECANE.eq.1) then
       eta=visccoef
      endif
 
      return
      end subroutine VISC_dodecane


      subroutine geom_avg(local_plus,local_minus,wt_plus,wt_minus, &
                      coeff)
      use probcommon_module
      IMPLICIT NONE

      REAL_T, intent(in) :: local_plus,local_minus
      REAL_T, intent(in) :: wt_plus,wt_minus
      REAL_T, intent(out) :: coeff

      if ((wt_plus.ge.zero).and.(wt_minus.ge.zero).and. &
          (abs(wt_plus+wt_minus-one).le.VOFTOL)) then
       ! do nothing
      else
       print *,"wt_plus+wt_minus should be 1"
       stop
      endif

      if (1.eq.0) then
       coeff=wt_plus*local_plus+wt_minus*local_minus
      else
       if ((local_plus.eq.zero).and.(local_minus.eq.zero)) then
        coeff=zero
       else if ((local_plus.eq.zero).and.(local_minus.gt.zero)) then
        coeff=zero
       else if ((local_plus.gt.zero).and.(local_minus.eq.zero)) then
        coeff=zero
       else if ((local_plus.gt.zero).and.(local_minus.gt.zero)) then
        coeff=local_plus*local_minus/(wt_plus*local_minus+wt_minus*local_plus)
       else
        print *,"local_plus or local_minus invalid"
        stop
       endif
      endif

      return
      end subroutine geom_avg


        ! CONTAINER ROUTINE FOR MEHDI VAHAB, MITSUHIRO OHTA, and
        ! MARCO ARIENTI
      REAL_T function get_user_viscconst(im,density,temperature)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T im,nmat,ibase,stage
      REAL_T density,temperature,mu

      nmat=num_materials
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im out of range"
       stop
      endif

      if (temperature.le.zero) then
       print *,"temperature invalid"
       stop
      endif
      if (density.le.zero) then
       print *,"density invalid in get_user_viscconst"
       stop
      endif

      if ((fort_viscconst(im).lt.zero).or. &
          (fort_prerecalesce_viscconst(im).lt.zero)) then
       print *,"fortran viscconst invalid"
       stop
      endif

      get_user_viscconst=fort_viscconst(im)

      if (recalesce_material(im).eq.0) then
       get_user_viscconst=fort_viscconst(im)
      else if ((recalesce_material(im).eq.1).or. &
               (recalesce_material(im).eq.2)) then
       ibase=(im-1)*recalesce_num_state
       stage=NINT(recalesce_state_old(ibase+1))
         ! stage=-1 (init)
         ! stage=0 (cooling)
         ! stage=1 (nucleation)
         ! stage=2 (recalesce in progress)
         ! stage=3 (recalesce finished)
         ! stage=4 (frost)
         ! stage=5 (regular freezing starts)
       if (stage.lt.3) then
        get_user_viscconst=fort_prerecalesce_viscconst(im)
       else if ((stage.ge.3).and.(stage.le.5)) then
        get_user_viscconst=fort_viscconst(im)
       else
        print *,"stage invalid"
        stop
       endif
      else
       print *,"recalesce_material invalid"
       stop
      endif

      if (fort_viscosity_state_model(im).eq.0) then
       ! do nothing
      else if (fort_viscosity_state_model(im).eq.1) then
       call VISC_dodecane(density,temperature,fort_viscconst(im),mu)
       get_user_viscconst=mu
      else if (fort_viscosity_state_model(im).eq.2) then
       print *,"viscosity_state_model(im)=2 is for Mitsuhiro."
       stop
      else
       print *,"fort_viscosity_state_model invalid"
       stop
      endif

      return
      end function get_user_viscconst

        ! CONTAINER ROUTINE FOR MEHDI VAHAB
      REAL_T function get_user_heatviscconst(im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T im,nmat,ibase,stage

      nmat=num_materials
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im out of range"
       stop
      endif
      if ((fort_heatviscconst(im).lt.zero).or. &
          (fort_prerecalesce_heatviscconst(im).lt.zero)) then
       print *,"fortran heatviscconst invalid"
       stop
      endif
      get_user_heatviscconst=fort_heatviscconst(im)
      if (recalesce_material(im).eq.0) then
       get_user_heatviscconst=fort_heatviscconst(im)
      else if ((recalesce_material(im).eq.1).or. &
               (recalesce_material(im).eq.2)) then
       ibase=(im-1)*recalesce_num_state
       stage=NINT(recalesce_state_old(ibase+1))
         ! stage=-1 (init)
         ! stage=0 (cooling)
         ! stage=1 (nucleation)
         ! stage=2 (recalesce in progress)
         ! stage=3 (recalesce finished)
         ! stage=4 (frost)
         ! stage=5 (regular freezing starts)
       if (stage.lt.3) then
        get_user_heatviscconst=fort_prerecalesce_heatviscconst(im)
       else if ((stage.ge.3).and.(stage.le.5)) then
        get_user_heatviscconst=fort_heatviscconst(im)
       else
        print *,"stage invalid"
        stop
       endif
      else
       print *,"recalesce_material invalid"
       stop
      endif

      return
      end function get_user_heatviscconst


        ! CONTAINER ROUTINE FOR MEHDI VAHAB
      REAL_T function get_user_stiffCP(im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: im
      INTEGER_T :: nmat,ibase,stage

      nmat=num_materials
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im out of range"
       stop
      endif
      if ((fort_stiffCP(im).lt.zero).or. &
          (fort_prerecalesce_stiffCP(im).lt.zero)) then
       print *,"fortran stiffCP invalid"
       stop
      endif
      get_user_stiffCP=fort_stiffCP(im)
      if (recalesce_material(im).eq.0) then
       get_user_stiffCP=fort_stiffCP(im)
      else if ((recalesce_material(im).eq.1).or. &
               (recalesce_material(im).eq.2)) then
       ibase=(im-1)*recalesce_num_state
       stage=NINT(recalesce_state_old(ibase+1))
         ! stage=-1 (init)
         ! stage=0 (cooling)
         ! stage=1 (nucleation)
         ! stage=2 (recalesce in progress)
         ! stage=3 (recalesce finished)
         ! stage=4 (frost)
         ! stage=5 (regular freezing starts)
       if (stage.lt.3) then
        get_user_stiffCP=fort_prerecalesce_stiffCP(im)
       else if ((stage.ge.3).and.(stage.le.5)) then
        get_user_stiffCP=fort_stiffCP(im)
       else
        print *,"stage invalid"
        stop
       endif
      else
       print *,"recalesce_material invalid"
       stop
      endif

      return
      end function get_user_stiffCP

      subroutine get_user_tension(xpos,time, &
        tension,new_tension, &
        temperature, &
        nmat,nten,caller_id)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: caller_id 
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: xpos(SDIM)
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: nten
      INTEGER_T nten_test
      REAL_T, intent(in) :: temperature(nmat)
      REAL_T, intent(in) :: tension(nten)
      REAL_T, intent(out) :: new_tension(nten)
      REAL_T avgtemp
      INTEGER_T iten,im,im_opp,ibase,stage

      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid get_user_tension nten nten test", &
         nten,nten_test
       print *,"nmat=",nmat
       print *,"caller_id=",caller_id
       stop
      endif

       ! fort_tension
       ! fort_prefreeze_tension
      do iten=1,nten
       new_tension(iten)=tension(iten)
       if (fort_tension_min(iten).lt.zero) then
        print *,"fort_tension_min invalid"
        stop
       endif
       if (fort_tension_slope(iten).ne.zero) then
        if (fort_tension_T0(iten).le.zero) then
         print *,"T0 invalid"
         stop
        endif
         ! im<im_opp
        call get_inverse_iten(im,im_opp,iten,nmat)
        if ((temperature(im).le.zero).or. &
            (temperature(im_opp).le.zero)) then
         print *,"temperature must be positive"
         stop
        endif
        avgtemp=half*(temperature(im)+temperature(im_opp))
        new_tension(iten)=new_tension(iten)+fort_tension_slope(iten)* &
         (avgtemp-fort_tension_T0(iten))
        if (new_tension(iten).lt.fort_tension_min(iten)) then
         new_tension(iten)=fort_tension_min(iten)
        endif
        if (new_tension(iten).lt.zero) then
         print *,"new_tension invalid"
         stop
        endif
       endif
      enddo ! iten
      do im=1,nmat
       if (recalesce_material(im).eq.0) then
        ! do nothing
       else if ((recalesce_material(im).eq.1).or. &
                (recalesce_material(im).eq.2)) then
        ibase=(im-1)*recalesce_num_state
        stage=NINT(recalesce_state_old(ibase+1))
         ! stage=-1 (init)
         ! stage=0 (cooling)
         ! stage=1 (nucleation)
         ! stage=2 (recalesce in progress)
         ! stage=3 (recalesce finished)
         ! stage=4 (frost)
         ! stage=5 (regular freezing starts)
        if (stage.lt.5) then
         do iten=1,nten
          new_tension(iten)=fort_prefreeze_tension(iten)
         enddo
        else if ((stage.eq.5).or.(stage.eq.6)) then
         ! do nothing
        else
         print *,"stage invalid"
         stop
        endif
       else
        print *,"recalesce_material invalid"
        stop
       endif
      enddo ! im=1..nmat

      return
      end subroutine get_user_tension


      subroutine TEMPERATURE_default(rho,temperature,internal_energy, &
        imattype,im)

      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: imattype
      INTEGER_T, intent(in) :: im
      INTEGER_T :: nmat
      REAL_T, intent(in) :: rho
      REAL_T, intent(out) :: temperature
      REAL_T, intent(in) :: internal_energy
      REAL_T cv

      nmat=num_materials
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid70"
       stop
      endif
!      cv=4.1855D+7
      cv=get_user_stiffCP(im)
      if (cv.le.zero) then
       print *,"cv invalid in temperature default"
       stop
      endif
      if (imattype.eq.999) then
       if (is_rigid(nmat,im).eq.0) then
        print *,"is_rigid invalid"
        stop
       endif
!       temperature=fort_tempconst(im)
       temperature=internal_energy/cv
      else if (imattype.eq.0) then
       temperature=internal_energy/cv
      else
       print *,"imattype invalid TEMPERATURE_default"
       stop
      endif

      return
      end subroutine TEMPERATURE_default


      subroutine INTERNAL_default(rho,temperature,internal_energy, &
        imattype,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho
      REAL_T, intent(in) :: temperature
      REAL_T, intent(out) :: internal_energy
      REAL_T cv


      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid in internal_default"
       stop
      endif
      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (temperature.le.zero) then
       print *,"T invalid"
       stop
      endif

!      cv=4.1855D+7
      cv=get_user_stiffCP(im)
      if (cv.le.zero) then
       print *,"cv invalid in internal default"
       stop
      endif
      if (imattype.eq.999) then
       internal_energy=cv*temperature
      else if (imattype.eq.0) then
       internal_energy=cv*temperature
      else
       print *,"imattype invalid in internal default"
       stop
      endif

      return
      end subroutine INTERNAL_default


        ! ice behaves like rigid solid where dist>0
      subroutine ice_substrate_distance(x,y,z,dist)
      use probcommon_module
      IMPLICIT NONE
      
      REAL_T, intent(in) :: x,y,z
      REAL_T, intent(out) :: dist
      INTEGER_T nmat
      REAL_T aspect,yprime,zprime,aspect2

      if (SDIM.eq.2) then
       if (abs(y-z).gt.VOFTOL) then
        print *,"y=z in 2d expected"
        stop
       endif
      endif

      nmat=num_materials

      aspect=tan(radblob2)
      if (SDIM.eq.2) then
        yprime=aspect*(x-xblob2)+yblob2
        dist=y-yprime  ! vertical distance
      else if (SDIM.eq.3) then
        aspect2=tan(radblob3)
        zprime=aspect*(x-xblob2)+aspect2*(y-yblob2)+zblob2
        dist=z-zprime  ! vertical distance
      else
        print *,"dimension bust"
        stop
      endif
      dist=-dist

      end subroutine ice_substrate_distance

        ! Sato and Niceno or Tryggvason
        ! (probtype.eq.55) and ((axis_dir.eq.6).or.(axis_dir.eq.7))
        ! negative if x_point in a bubble.
      subroutine nucleation_sites(x_point,dist,nucleate_pos)
      use probcommon_module

      IMPLICIT NONE

      REAL_T, intent(in) :: x_point(SDIM)
      REAL_T, intent(out) :: dist
      INTEGER_T icomp
      REAL_T distarr(n_sites)
      REAL_T, intent(in) :: nucleate_pos(4*n_sites)
      REAL_T hugedist
      REAL_T xx(SDIM)
      REAL_T rr
      INTEGER_T dir
 
      if (n_sites.lt.1) then
       print *,"n_sites invalid (GLOBALUTIL.F90), n_sites=",n_sites
       stop
      endif
 
      hugedist=99999.0

      dist=hugedist

      do icomp=1,n_sites
       distarr(icomp)=hugedist
       rr=nucleate_pos(4*(icomp-1)+4)
       if (rr.gt.zero) then
        do dir=1,SDIM
         xx(dir)=x_point(dir)-nucleate_pos(4*(icomp-1)+dir)
        enddo
        distarr(icomp)=zero
        do dir=1,SDIM
         distarr(icomp)=distarr(icomp)+xx(dir)**2
        enddo
        distarr(icomp)=sqrt(distarr(icomp))-rr
        
        dist=min(dist,distarr(icomp))
       else
        print *,"rr invalid"
        stop
       endif 
      enddo ! icomp=1..n_sites
    
      return
      end subroutine nucleation_sites


      subroutine get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      use probcommon_module
      IMPLICIT NONE

      REAL_T A,B,GAMMA,R1,R2,RHOI


      if ((probtype.eq.36).and.(axis_dir.eq.2)) then  ! spherical explosion
       A=5.484D+12
       B=0.09375D+12
       R1=4.94D0
       R2=1.21D0
       GAMMA=1.28D0
       RHOI=1.63D0
      else if (probtype.eq.42) then  ! bubble jetting
       A=3.712D+12
       B=0.03231D+12
       R1=4.15D0
       R2=0.95D0
       GAMMA=1.3D0
       RHOI=1.63D0
      else if (probtype.eq.46) then ! cavitation
       A=3.712D+12
       B=0.03231D+12
       R1=4.15D0
       R2=0.95D0
       GAMMA=1.3D0
       RHOI=1.63D0
      else
       print *,"probtype not supported for the jwl material"
       stop
      endif

      return
      end subroutine get_jwl_constants



      subroutine INTERNAL_jwl(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif

      cv=4.1855D+7
      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_jwl

      subroutine TEMPERATURE_jwl(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif

      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine

      subroutine ENTROPY_jwl(rho,internal_energy,entropy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure,entropy,press_adiabat

      call EOS_NAjwl(rho,internal_energy,pressure)
      call EOS_jwlADIABAT(rho,internal_energy,press_adiabat)
      if (press_adiabat.le.zero) then
       print *,"press_adiabat invalid"
       stop
      endif
      entropy=pressure/press_adiabat

      return
      end subroutine ENTROPY_jwl

      subroutine INTERNAL_ENTROPY_jwl(rho,entropy,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,entropy,internal_energy
      REAL_T A,B,R1,R2,GAMMA,RHOI,OMEGA,pressure_part,press_adiabat

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      OMEGA=GAMMA-one
      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (entropy.le.zero) then
       print *,"entropy invalid"
       stop
      endif
      pressure_part= &
        A*(one-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(one-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)
      call EOS_jwlADIABAT(rho,internal_energy,press_adiabat)
      internal_energy=(press_adiabat*entropy-pressure_part)/ &
        (OMEGA*rho)
      if (internal_energy.le.zero) then
       print *,"internal_energy invalid"
       stop
      endif

      return
      end subroutine INTERNAL_ENTROPY_jwl

 
! e=(E/rho) - (1/2) (u^2 + v^2)
      subroutine EOS_NAjwl(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,R1,R2,GAMMA,RHOI,OMEGA

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      OMEGA=GAMMA-one

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif
      pressure= &
        A*(one-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(one-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)+ &
        OMEGA*rho*internal_energy

      if (pressure.le.zero) then
       print *,"vacuum error in NA JWL"
       stop
      endif

      return
      end subroutine

! initial sound speed is:
! C=7.8039D+10-5.484D+12 e^(-4.94)-0.09375D+12 e^(-1.21)=
      subroutine SOUNDSQR_NAjwl(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,R1,R2,GAMMA,RHOI,OMEGA
      REAL_T pressure,dp_de,dp_drho

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      OMEGA=GAMMA-one

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      call EOS_NAjwl(rho,internal_energy,pressure)
      dp_de=OMEGA*rho
      dp_drho= &
        A*(one-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)* &
        R1*RHOI/(rho**2)- &
        (A*OMEGA/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(one-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)* &
        R2*RHOI/(rho**2)- &
        B*(OMEGA/(R2*RHOI))*exp(-R2*RHOI/rho)+ &
        OMEGA*internal_energy
    
      soundsqr=(pressure*dp_de)/(rho**2)+dp_drho
 
      if (soundsqr.le.zero) then
       print *,"cannot have 0 sound speed"
       stop
      endif

      return
      end subroutine


! e=(E/rho) - (1/2) (u^2 + v^2)
      subroutine EOS_jwlADIABAT(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,R1,R2,GAMMA,RI,PI,RHOI,C,OMEGA

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      if (RHOI.le.zero) then
       print *,"RHOI invalid"
       stop
      endif

      OMEGA=GAMMA-one
      RI=16.0D0        ! cm
      PI=7.8039D+10  ! dyne/cm^2
      C=PI-OMEGA*(A*exp(-R1)/R1+B*exp(-R2)/R2)
      if (C.le.zero) then
       print *,"c invalid"
       stop
      endif
      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif
      pressure=(OMEGA*rho/RHOI)*( &
       A*exp(-R1*RHOI/rho)/R1+ &
       B*exp(-R2*RHOI/rho)/R2   )+ &
       C*( (RHOI/rho)**(-GAMMA) )
      if (pressure.le.zero) then
       print *,"vacuum error"
       stop
      endif

      return
      end subroutine

      subroutine SOUNDSQR_jwlADIABAT(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,R1,R2,GAMMA,RI,PI,RHOI,C,OMEGA

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      if (RHOI.le.zero) then
       print *,"RHOI invalid"
       stop
      endif

      OMEGA=GAMMA-one
      RI=16.0D0        ! cm
      PI=7.8039D+10  ! dyne/cm^2
      C=PI-OMEGA*(A*exp(-R1)/R1+B*exp(-R2)/R2)
      if (C.le.zero) then
       print *,"c invalid"
       stop
      endif
      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif
      soundsqr=(OMEGA/RHOI)*(  &
       A*exp(-R1*RHOI/rho)/R1+ &
       B*exp(-R2*RHOI/rho)/R2)
      soundsqr=soundsqr+ &
       (OMEGA*rho/RHOI)*( &
       A*(RHOI/(rho**2))*exp(-R1*RHOI/rho)+ &
       B*(RHOI/(rho**2))*exp(-R2*RHOI/rho)   )+ &
       (GAMMA*C/RHOI)*( (rho/RHOI)**(GAMMA-one) )

      if (soundsqr.le.zero) then
       print *,"cannot have 0 sound speed"
       stop
      endif

      return
      end subroutine SOUNDSQR_jwlADIABAT

! PR without shift:
! https://en.wikipedia.org/wiki/Equation_of_state
! PR with shift:
! Monnery et al, Ind. Eng. Chem. Res. 1998, volume 37, 1663-1672
! table: c_1, c_2, c_3, m, n
! rho  g/cm^3=rho/1000 kg/cm^3=rho 10^6/10^3  kg/m^3 = 1000 rho kg/m^3
      subroutine EOS_peng_robinson(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T temperature
      REAL_T R_pr,water_molar_mass,rhoMKS,rho_molar
      REAL_T Vm,water_critical_temperature,water_critical_pressure
      REAL_T water_critical_molar_volume,a,b,Tr,kappa,alpha
      REAL_T water_acentric_factor
      REAL_T water_c1,water_c2,water_c3,water_m,water_n
      REAL_T Monnery_c4
      REAL_T Vm_shift,adjusted_Vm

      if ((rho.gt.zero).and. &
          (internal_energy.gt.zero)) then

       call TEMPERATURE_peng_robinson(rho,temperature,internal_energy)
       R_pr=8.3144621 ! J/(mol K)  (ideal gas constant)
       water_c1=-0.19963 ! m^3/kmol
       water_c2=2.03032
       water_c3=0.38268
       Monnery_c4=1.0  ! m^3/kmol
       water_m=0.91132
       water_n=-2.285D-1
       water_molar_mass=18.01528 ! g/mol
       rhoMKS=1000.0*rho  ! kg/m^3
       rho_molar=rhoMKS/(0.001*water_molar_mass)  ! mol/m^3
        ! if rho=1g/cm^3 then rhoMKS=1000 kg/m^3 , rho_molar=1D+6/18 mol/m^3
        !  Vm=18D-6 m^3/mole 
       Vm=one/rho_molar     ! m^3/mole

       water_critical_temperature=647.0  ! degrees Kelvin
       water_critical_pressure=22.064D+6 ! Pascals=N/m^2
        ! 55.9 cm^3/mol=55.9D-6 m^3/mole=5.59D-5 m^3/mole
       water_critical_molar_volume=5.59D-5 ! m^3/mole
        ! units: J^2/(mol^2 N/m^2)=J^2 m^2/(mol^2 N)
        !  =(kg^2 m^4/s^4)m^2/(mol^2 kg m/s^2)=
        !   kg m^5/(mol^2 s^2)
        ! note: units of a/Vm^2=kg m^5/(mol^2 s^2)/(m^6/mole^2)=
        !  kg/(m s^2)=(kg m/s^2)(1/m^2)=N/m^2 
       a=0.45724*(R_pr**2)*(water_critical_temperature**2)/ &
             water_critical_pressure
       b=0.07780*R_pr*water_critical_temperature/water_critical_pressure
       Tr=temperature/water_critical_temperature
       water_acentric_factor=0.344
       kappa=0.37464+1.54226*water_acentric_factor- &
             0.26992*(water_acentric_factor**2)
!      alpha=(one+kappa*(one-sqrt(Tr)))**2 ! original Peng Robinson model.
        ! this is from Monnery et al
       alpha=(one+water_m*(one-sqrt(Tr))-water_n*((one-sqrt(Tr))**2))**2 

       Vm_shift=water_c1+(Monnery_c4/(sqrt(two*Pi)*water_c2))* &
            exp(-half*( ((Tr-water_c3)/water_c2)**2 )) !m^3/kmol
       Vm_shift=Vm_shift/1000.0  ! m^3/mol
       adjusted_Vm=Vm-Vm_shift

       if (adjusted_Vm.gt.b) then
        ! N/m^2
        pressure=R_pr*temperature/(adjusted_Vm-b)- &
          a*alpha/(adjusted_Vm**2+two*b*adjusted_Vm-(b**2))
        if (pressure.le.zero) then
         print *,"pressure must be positive"
         print *,"rho,temperature ",rho,temperature
         print *,"Vm= ",Vm
         print *,"Vm offset = ",Vm_shift
         print *,"adjusted_Vm=",adjusted_Vm
         print *,"a=",a
         print *,"b=",b
         print *,"alpha=",alpha
         print *,"part1=",R_pr*temperature/(adjusted_Vm-b)
         print *,"part2=",a*alpha/(adjusted_Vm**2+two*b*adjusted_Vm-(b**2))
         stop
        endif
       else 
        print *,"adjusted_Vm-b must be positive"
        print *,"rho,temperature ",rho,temperature
        stop
       endif

        ! N/m^2 = kg m/s^2 / m^2 = kg/(m s^2)=1000 g/(100 cm s^2)=
        ! 10 g/(cm s^2)=10 dyne/cm^2
       pressure=pressure*ten

      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine EOS_peng_robinson

! c^2 = dp/drho + p dp/de / rho^2
      subroutine SOUNDSQR_peng_robinson(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr,pressure
      REAL_T eps,drho,de,rho_plus,rho_minus,p_plus,p_minus
      REAL_T e_plus,e_minus,dp_drho,dp_de

      if ((rho.gt.zero).and. &
          (internal_energy.gt.zero)) then

       call EOS_peng_robinson(rho,internal_energy,pressure)
       eps=1.0D-6
       drho=eps*rho
       de=eps*internal_energy
       rho_plus=rho+half*drho
       rho_minus=rho-half*drho
       call EOS_peng_robinson(rho_plus,internal_energy,p_plus)
       call EOS_peng_robinson(rho_minus,internal_energy,p_minus)
       dp_drho=(p_plus-p_minus)/drho
       e_plus=internal_energy+half*de
       e_minus=internal_energy-half*de
       call EOS_peng_robinson(rho,e_plus,p_plus)
       call EOS_peng_robinson(rho,e_minus,p_minus)
       dp_de=(p_plus-p_minus)/de
       soundsqr=dp_drho+pressure*dp_de/(rho**2)
      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine SOUNDSQR_peng_robinson

      subroutine INTERNAL_peng_robinson(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_peng_robinson

      subroutine TEMPERATURE_peng_robinson(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_peng_robinson

! A. Brundage, ``Implementation of Tillotson Equation of State for 
! Hypervelocity Impact of Metals, Geologic Materials, and Liquids.''
! Procedia Engineering 58 (2013) 461-470
! material_type=22
! NOTE: the internal energy that is input to this routine is e=cv T
! but the internal energy as a function of temperature for the Tillotson
! model has to be defined such that P=P0 when rho=rho0.
! a) 
! one possible modification to the Tillotson EOS (but the sound speed is
! way too small):
! e'=e+alpha
! P0=(a +  b/(e0'/E0 + 1))rho0 e0' + A mu + B mu^2
! e0'=cv T0 + alpha
! assume mu=0
! P0=(a+b/(e0'/E0 + 1))rho0 e0'
! P0(e0'/E0 + 1) = a rho0 e0' (e0'/E0 + 1) + b rho0 e0'
! let x=e0'
! (a rho0/E0) x^2 + (b rho0 + a rho0 - P0/E0)x - P0 = 0
! b)
! another modification: P=P(rho,e)-P(rho0,e0)+P0
      subroutine EOS_tillotson(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure,strain
      REAL_T rho0,eta,E_ratio,denom,a_term,b_term,beta_term,alpha_term
      REAL_T pressure1,pressure3
      REAL_T e0,E_tillotson,T0_tillotson
      REAL_T E_offset,P_offset
      REAL_T E0_ratio,denom0,a0_term,b0_term

      rho0=fort_denconst(1)
      T0_tillotson=fort_tempconst(1)
      strain=zero

      if ((rho.gt.zero).and. &
          (internal_energy.gt.zero).and. &
          (rho0.gt.zero).and. &
          (rho0.gt.rho_IV_tillotson).and. &
          (T0_tillotson.gt.zero).and. &
          (E0_tillotson.gt.zero).and. &
          (E_CV_tillotson.gt.E_IV_tillotson)) then

       call INTERNAL_tillotson(rho0,T0_tillotson,e0)

       eta=rho/rho0
       E_tillotson=internal_energy
       E_offset=E_tillotson
       E_ratio=E_tillotson/E0_tillotson
       denom=one+E_ratio/(eta**2)

       E0_ratio=e0/E0_tillotson
       denom0=one+E0_ratio
       a0_term=a_hydro_tillotson*rho0*e0
       b0_term=b_hydro_tillotson*rho0*e0

       P_offset=P0_tillotson-(a0_term+b0_term/denom0)

       a_term=a_hydro_tillotson*rho*E_offset
       b_term=b_hydro_tillotson*rho*E_offset
       beta_term=exp(beta_tillotson*(one-one/eta))
       alpha_term=exp(-alpha_tillotson*((one-one/eta)**2))

       pressure1=a_term+b_term/denom+A_strain_tillotson*strain+ &
         B_strain_tillotson*(strain**2)+P_offset

       pressure3=a_term+(b_term/denom+ &
        A_strain_tillotson*strain*beta_term)*alpha_term+P_offset

       if ((rho.ge.rho0).or. &
           ((rho.ge.rho_IV_tillotson).and. &
            (E_tillotson.le.E_IV_tillotson))) then
        pressure=pressure1
       else if ((rho.le.rho0).and. &
                (E_tillotson.ge.E_CV_tillotson)) then
        pressure=pressure3
       else if ((rho.ge.rho_IV_tillotson).and. &
                (rho.le.rho0).and. &
                (E_tillotson.ge.E_IV_tillotson).and. &
                (E_tillotson.le.E_CV_tillotson)) then
        pressure=((E_tillotson-E_IV_tillotson)*pressure3+ &
                  (E_CV_tillotson-E_tillotson)*pressure1)/ &
                 (E_CV_tillotson-E_IV_tillotson)
       else
        pressure=a_term+b_term/denom+A_strain_tillotson*strain+P_offset
       endif
      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine EOS_tillotson

! c^2 = dp/drho + p dp/de / rho^2
      subroutine SOUNDSQR_tillotson(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr,pressure
      REAL_T rho0,eps,drho,de,rho_plus,rho_minus,p_plus,p_minus
      REAL_T e_plus,e_minus,dp_drho,dp_de

      rho0=fort_denconst(1)

      if ((rho.gt.zero).and. &
          (internal_energy.gt.zero).and. &
          (rho0.gt.zero).and. &
          (rho0.gt.rho_IV_tillotson).and. &
          (E0_tillotson.gt.zero).and. &
          (E_CV_tillotson.gt.E_IV_tillotson)) then

       call EOS_tillotson(rho,internal_energy,pressure)
       eps=1.0D-6
       drho=eps*rho
       de=eps*internal_energy
       rho_plus=rho+half*drho
       rho_minus=rho-half*drho
       call EOS_tillotson(rho_plus,internal_energy,p_plus)
       call EOS_tillotson(rho_minus,internal_energy,p_minus)
       dp_drho=(p_plus-p_minus)/drho
       e_plus=internal_energy+half*de
       e_minus=internal_energy-half*de
       call EOS_tillotson(rho,e_plus,p_plus)
       call EOS_tillotson(rho,e_minus,p_minus)
       dp_de=(p_plus-p_minus)/de
       soundsqr=dp_drho+pressure*dp_de/(rho**2)
      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine SOUNDSQR_tillotson


      subroutine INTERNAL_tillotson(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if ((rho.gt.zero).and.(temperature.gt.zero)) then
       cv=4.1855D+7
       internal_energy=cv*temperature
      else
       print *,"rho or temperature invalid"
       stop
      endif

      return
      end subroutine INTERNAL_tillotson


      subroutine TEMPERATURE_tillotson(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if ((rho.gt.zero).and.(internal_energy.gt.zero)) then
       cv=4.1855D+7
       temperature=internal_energy/cv
      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine TEMPERATURE_tillotson


! A fully compressible, two-dimensional model of small, high-speed, cavitating
! nozzles, Schmidt et al, Atomization and Sprays, vol 9, 255-276 (1999)
! units must be CGS
      subroutine EOS_cav_nozzle(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T rho_g,rho_l,soundsqr_g,soundsqr_l,alpha,pgl,numer,denom

      rho_g=fort_cavdenconst(1)
      rho_l=fort_denconst(1)

      if ((rho.gt.zero).and.(internal_energy.gt.zero)) then
       if ((rho_l.gt.zero).and.(rho_g.gt.zero)) then
        call SOUNDSQR_tait(rho_l,internal_energy,soundsqr_l)
        call SOUNDSQR_air(rho_g,internal_energy,soundsqr_g)

        if ((rho_l.gt.rho_g).and. &
            (soundsqr_l.gt.soundsqr_g).and. &
            (soundsqr_l.gt.zero).and. &
            (soundsqr_g.gt.zero)) then

          ! if rho=rho_l then alpha=0, numer=rho_g * rho_l*soundsqr_g
          !  denom=rho_l * rho_g *soundsqr_g
          !  pressure=PCAV_TAIT
          ! if rho=rho_g then alpha=1, numer=rho_g*rho_g*soundsqr_g
          !  denom=rho_l*rho_l*soundsqr_l
          !  pressure=PCAV_TAIT+pgl*log(numer/denom) 
         if (rho.ge.rho_l) then
          alpha=zero
         else if (rho.le.rho_g) then
          alpha=one
         else
          alpha=(rho-rho_l)/(rho_g-rho_l)
         endif
         pgl=rho_g*soundsqr_g*rho_l*soundsqr_l*(rho_g-rho_l)/ &
                 (soundsqr_g*(rho_g**2)-soundsqr_l*(rho_l**2))
         if (pgl.gt.zero) then
          numer=rho_g*soundsqr_g*(alpha*rho_g+(one-alpha)*rho_l)
          denom=rho_l*(alpha*rho_l*soundsqr_l+(one-alpha)* &
                  rho_g*soundsqr_g)
          if ((denom.gt.zero).and.(numer.gt.zero)) then
           pressure=PCAV_TAIT+pgl*log(numer/denom)
           if (rho.ge.rho_l) then
            pressure=pressure+soundsqr_l*(rho-rho_l)
           else if (rho.le.rho_g) then
            pressure=pressure+soundsqr_g*(rho-rho_g)
           else if ((rho.ge.rho_g).and.(rho.le.rho_l)) then
            ! do nothing (pressure already set)
           else
            print *,"rho invalid"
            stop
           endif  
          else
           print *,"denom or numer invalid"
           stop
          endif
         else
          print *,"pgl invalid"
          stop
         endif
        else
         print *,"rho or soundsqr invalid"
         stop
        endif
       else
        print *,"rho invalid"
        stop
       endif
      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine EOS_cav_nozzle


      subroutine SOUNDSQR_cav_nozzle(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T rho_g,rho_l,soundsqr_g,soundsqr_l,alpha,rhomix,inv_soundsqr_mix

      rho_g=fort_cavdenconst(1)
      rho_l=fort_denconst(1)

      if ((rho.gt.zero).and.(internal_energy.gt.zero)) then
       if ((rho_l.gt.zero).and.(rho_g.gt.zero)) then
        call SOUNDSQR_tait(rho_l,internal_energy,soundsqr_l)
        call SOUNDSQR_air(rho_g,internal_energy,soundsqr_g)

        if ((rho_l.gt.rho_g).and. &
            (soundsqr_l.gt.soundsqr_g).and. &
            (soundsqr_l.gt.zero).and. &
            (soundsqr_g.gt.zero)) then
         if (rho.ge.rho_l) then
          alpha=zero
         else if (rho.le.rho_g) then
          alpha=one
         else
          alpha=(rho-rho_l)/(rho_g-rho_l)
         endif
         rhomix=alpha*rho_g+(one-alpha)*rho_l
         inv_soundsqr_mix=alpha/(rho_g*soundsqr_g)+ &
            (one-alpha)/(rho_l*soundsqr_l)
         if ((rhomix.gt.zero).and.(inv_soundsqr_mix.gt.zero)) then 
          soundsqr=one/(rhomix*inv_soundsqr_mix)
         else
          print *,"rhomix or inv_soundsqr invalid"
          stop
         endif
        else
         print *,"rho or soundsqr invalid"
         stop
        endif
       else
        print *,"rho invalid"
        stop
       endif
      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine SOUNDSQR_cav_nozzle

      subroutine INTERNAL_cav_nozzle(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if ((rho.gt.zero).and.(temperature.gt.zero)) then
       cv=4.1855D+7
       internal_energy=cv*temperature
      else
       print *,"rho or temperature invalid"
       stop
      endif

      return
      end subroutine INTERNAL_cav_nozzle


      subroutine TEMPERATURE_cav_nozzle(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if ((rho.gt.zero).and.(internal_energy.gt.zero)) then
       cv=4.1855D+7
       temperature=internal_energy/cv
      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine TEMPERATURE_cav_nozzle


      subroutine INTERNAL_tait(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_tait


      subroutine INTERNAL_tait_vacuum(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_tait_vacuum


      subroutine TEMPERATURE_tait(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_tait


      subroutine TEMPERATURE_tait_vacuum(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_tait_vacuum


      subroutine EOS_tait(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,rhobar,GAMMA,pcav

      A=A_TAIT
      B=B_TAIT
      rhobar=RHOBAR_TAIT ! see PROBCOMMON.F90 for definition
      GAMMA=GAMMA_TAIT
      pcav=PCAV_TAIT

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine EOS_tait

      subroutine EOS_tait_vacuum(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,pressure
      REAL_T A,B,rhobar,GAMMA,pcav

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=RHOBAR_TAIT ! g/cm^3
      GAMMA=GAMMA_TAIT
      pcav=PCAV_TAIT_VACUUM 

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine EOS_tait_vacuum




      subroutine SOUNDSQR_tait(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,rhobar,pcav,rhocav,pressure
      REAL_T rho_sound

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=RHOBAR_TAIT ! g/cm^3
      pcav=PCAV_TAIT
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (rhocav.le.zero) then
       print *,"rhocav invalid"
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_TAIT*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_TAIT-one) )

      return
      end subroutine SOUNDSQR_tait



      subroutine SOUNDSQR_tait_vacuum(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,rhobar,pcav,rhocav,pressure
      REAL_T rho_sound

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=RHOBAR_TAIT ! g/cm^3
      pcav=PCAV_TAIT_VACUUM
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (rhocav.le.zero) then
       print *,"rhocav invalid SOUNDSQR_tait_vacuum"
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_TAIT*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_TAIT-one) )

      return
      end subroutine SOUNDSQR_tait_vacuum


! begin dodecane routines (material type=15)


      subroutine DeDT_dodecane(rho_in,T_in,DeDT)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,rho_in,T,T_in,DeDT
      REAL_T cv,Tc
      REAL_T A,B,rho0,rmax,rmin,P0,pressure
      REAL_T cd,cc,cb,ca,T2,T3

      rho=rho_in
      T=T_in

      P0 = 1D+6    ! dyne/cm^2 reference pressure (0.1 MPa)
      A = 0.08998
      rmin = 0.500 ! at 0.1 MPa, T = 433
      rmax = 0.850 ! at 400 MPa, T = 298
      Tc = 658.2  ! critical temperature
      cv = 2.116D+7 ! for dodecane at 106 MPa and 363 K

      if (rho.le.1e-6) then
       print *,"DeDT_dodecane: density negative or near zero ",rho
       !stop
       rho = rmin
      endif
      if (T.le.zero) then
       print *,"DeDT_dodecane: temperature <= 0"
       stop
      else if (T.lt.200.) then
!      print *,"DeDT_dodecane: temperature < 250 K: ",T,rho
       T = 200.
      else if (T.gt.Tc) then
!      print *,"DeDT_dodecane: temperature above critical: ",T,rho
       T = Tc
      endif
      rho = min(rho,rmax)
      rho = max(rho,rmin)

      T2 = T*T
      T3 = T2*T
      rho0 = 0.9291654-0.5174730e-3*T-3.338672e-7*T2
      B=(345.1-1.1458*T+0.9837e-3*T2)*1e7
      pressure=-B+(B+P0)*exp((1.0-rho0/rho)/A)
      pressure=max(pressure,0.01*P0)

      cd=2.273845+7.701613e-20*pressure*(pressure-1e6);
      cc=-2.279889e-3-3.654273e-13*(pressure-1e6);
      cb=6.106366e-06;
      ca=-3.266302e-09;

      DeDT=(4.*ca*T3+3.*cb*T2+2.*cc*T+cd)*1e7 ! in erg/cm3

      if(DeDT.le.zero) then
       print *,"DeDT_dodecane: negative derivative: ",DeDT
       DeDT = cv
      endif

      return
      end subroutine DeDT_dodecane


      subroutine INTERNAL_dodecane(rho_in,T_in,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho_in,T_in
      REAL_T rho,T,internal_energy,cv,Tc
      REAL_T A,B,rho0,rmax,rmin,P0,pressure
      REAL_T ce,cd,cc,cb,ca,T2,T3,T4

      rho=rho_in
      T=T_in

      P0 = 1D+6    ! dyne/cm^2 reference pressure (0.1 MPa)
      A = 0.08998
      rmin = 0.500 ! at 0.1 MPa, T = 433
      rmax = 0.850 ! at 400 MPa, T = 298
      Tc = 658.2  ! critical temperature
      cv = 2.116D+7 ! for dodecane at 106 MPa and 363 K

      if (rho.le.1e-6) then
       print *,"INTERNAL_dodecane: density negative or near zero ",rho
       !stop
       rho = rmin
      endif
      if (T.le.zero) then
       print *,"INTERNAL_dodecane: temperature <= 0",T
       T = 250.
!      stop
      endif
      rho = min(rho,rmax)
      rho = max(rho,rmin)

      T2 = T*T
      T3 = T2*T
      T4 = T3*T
      rho0 = 0.9291654-0.5174730e-3*T-3.338672e-7*T2
      B=(345.1-1.1458*T+0.9837e-3*T2)*1e7
      pressure=-B+(B+P0)*exp((1.0-rho0/rho)/A)
      pressure=max(pressure,0.0001*P0)

      ce=19.94245;
      cd=2.273845+7.701613e-20*pressure*(pressure-1e6);
      cc=-2.279889e-3-3.654273e-13*(pressure-1e6);
      cb=6.106366e-06;
      ca=-3.266302e-09;

      internal_energy=(ca*T4+cb*T3+cc*T2+cd*T+ce)*1e7 ! in erg/cm3

      if (OLD_DODECANE.eq.1) then
       cv=2.116D+7 ! for dodecane at 106 MPa and 363 K
       internal_energy=T*cv
      endif

      return
      end subroutine INTERNAL_dodecane


      subroutine TEMPERATURE_dodecane(rho_in,T,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho_in
      REAL_T rho,T,internal_energy,cv
      REAL_T T0,T1,Tc,ie1,ie0,dT,dedT,rmin,rmax
      INTEGER iter

      rho=rho_in

      Tc = 658.2  ! critical temperature
      cv = 2.116D+7 ! at 106 MPa and 363 K
      rmin = 0.500 ! at 0.1 MPa, T = 433
      rmax = 0.850 ! at 400 MPa, T = 298

      if (rho.le.1e-6) then
       print *,"TEMPERATURE_dodecane: density negative or near zero ",rho
       !stop
       rho = rmin
      endif
      rho = min(rho,rmax)
      rho = max(rho,rmin)

      ! internal energy can actually be < 0, but this value was
      ! originally translated
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif

      iter = 1
      T0 = min(internal_energy/cv,Tc-50.) ! first guess
      T0 = max(T0,250.)
      call INTERNAL_dodecane(rho,T0,ie0)

      do while (abs(ie0/internal_energy-1.).gt.1e-4) 

       T1 = T0+2.
       call INTERNAL_dodecane(rho,T1,ie1)
       if (abs(ie1-ie0).lt.1e-6) then
        print *,"TEMPERATURE_dodecane: internal energies too close"
        stop
       endif

       dedT = 0.5*(ie1-ie0) ! T1-T0 = 2K
       dT = (ie0-internal_energy)/dedT
       dT = sign(min(abs(dT),10.),dT)
       T0 = T0-dT
       call INTERNAL_dodecane(rho,T0,ie0)

       if (iter.gt.50) then
        print *,"TEMPERATURE_dodecane: temperature not converged",&
         rho,T0,internal_energy,dT
        goto 10
       endif
       iter=iter+1
      enddo

  10  T = T0

      if (T.gt.Tc) then
!      print *,"TEMPERATURE_dodecane: temperature above critical: ",T,&
!       rho,internal_energy/cv
       T = Tc
      else if (T.lt.200.) then
       print *,"TEMPERATURE_dodecane: temperature < 200 K: ",T,&
        rho,internal_energy/cv
       T = 200.
      endif

!     call INTERNAL_dodecane(rho,T,internal_energy)

      if (OLD_DODECANE.eq.1) then
       cv=2.116D+7 ! for dodecane at 106 MPa and 363 K
       T=internal_energy/cv
      endif

      return
      end subroutine TEMPERATURE_dodecane


      subroutine EOS_dodecane(rho_in,internal_energy,T,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho_in
      REAL_T rho,internal_energy,T,pressure
      REAL_T A,B,rho0,rmax,rmin,P0,pcav,cv,Tc,Tred

      rho=rho_in

      P0 = 1D+6    ! dyne/cm^2 reference pressure (0.1 MPa)
      A = 0.08998
      rmin = 0.500 ! at 0.1 MPa, T = 433
      rmax = 0.850 ! at 400 MPa, T = 298

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      rho = min(rho,rmax)
      rho = max(rho,rmin)

      call TEMPERATURE_dodecane(rho,T,internal_energy)

      rho0 = 0.9291654-0.5174730e-3*T-3.338672e-7*T*T
      B=(345.1-1.1458*T+0.9837e-3*T*T)*1e7
      pressure=-B+(B+P0)*exp((1.0-rho0/rho)/A)
      pressure=max(pressure,0.00001*P0)
      pressure=min(pressure,1.8e9)  ! VIP limiter

      if (OLD_DODECANE.eq.1) then
       pcav=PCAV_TAIT
       cv=2.116D+7 ! for dodecane at 106 MPa and 363 K
       P0 = 1D+6    ! dyne/cm^2 reference pressure
       rho0= 0.6973 ! dodecane at 0.1 MPa and 363 K
       A = 0.0881
       B = 5.657D+8 ! dyne/cm^2 at 363 K
       Tc = 658.6
       T = internal_energy/cv
       rho0 = 0.9291654-0.5174730*1e-3*T-3.338672*1D-7*T*T
       Tred = T/Tc
       pressure=-B+(B+P0)*exp((1.0-rho0/rho)/A)
       if (pressure.lt.pcav) then
        pressure=pcav
       endif
      endif

      return
      end subroutine EOS_dodecane

      subroutine EOS_dodecane_ADIABATIC(rho,pressure)
      use probcommon_module
      IMPLICIT NONE
      
      REAL_T rho,pressure
      REAL_T A,B,rhobar,pcav

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=RHOBAR_TAIT ! g/cm^3
      pcav=PCAV_TAIT

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      
      return
      end subroutine EOS_dodecane_ADIABATIC


      subroutine SOUNDSQR_dodecane(rho_in,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho_in
      REAL_T rmin,rmax
      REAL_T rho,internal_energy,soundsqr
      REAL_T c,c0,D,E,pressure,T,T0p5,T1p5,T3p0
      REAL_T A,B,rhobar,pcav,rhocav,rho_sound

      rho=rho_in

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      rmin = 0.500 ! at 0.1 MPa, T = 433
      rmax = 0.850 ! at 400 MPa, T = 298
      rho = min(rho,rmax)
      rho = max(rho,rmin)

      call EOS_dodecane(rho,internal_energy,T,pressure)

      T0p5 = T**0.5
      T1p5 = T0p5*T
      T3p0 = T1p5*T1p5
      ! sound speed at atmospheric conditions - cm/s
      c0 = (4094-183.21*T0p5+0.07974*T1p5-2.348e-6*T3p0)*100.
      if (c0.le.zero) then
       print *,"invalid sound speed at temperature ",T,": c0 =",c0,rho,T,pressure
       c0 = 2.25e10
!      stop
      endif
      ! sound speed at pressure 
     !D = 0.005208-5.1495e-4*T-5.55e-14*T*pressure; ! pressure needed in MPa
      D = 0.1652083+2.5e-3*T-5.85e-13*T*pressure; ! modified correlation
      E = (-56.91+7.3674e-5*T*T+0.02260*T+463.5*exp(-0.001687*T))*1e7 !  convert from dyne/cm^2 to MPa
      c=c0/(1-D*log((E+pressure)/(E+1D6))); ! already in cm/s
      soundsqr=c*c

      if (OLD_DODECANE.eq.1) then
       A=A_TAIT   ! dyne/cm^2
       B=B_TAIT  ! dyne/cm^2
       rhobar=RHOBAR_TAIT ! g/cm^3
       pcav=PCAV_TAIT
       rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
       pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

       if (rhocav.le.zero) then
        print *,"rhocav invalid"
        stop
       endif

       if (pressure.lt.pcav) then
        rho_sound=rhocav
       else
        rho_sound=rho
       endif
       soundsqr=(GAMMA_TAIT*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_TAIT-one) )
      endif

      return
      end subroutine SOUNDSQR_dodecane


! end dodecane routines (material type=15)



      subroutine INTERNAL_tait_rho(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=temperature*cv

      return
      end subroutine INTERNAL_tait_rho

      subroutine TEMPERATURE_tait_rho(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_tait_rho

      subroutine EOS_tait_rho(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,rhobar,pcav


      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos tait rho"
       stop
      endif

      pcav=PCAV_TAIT

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif

      return
      end subroutine EOS_tait_rho


      subroutine EOS_tait_ADIABATIC_rho(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,pressure
      REAL_T A,B,rhobar,pcav


      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos tait adiabatic rho"
       stop
      endif

      pcav=PCAV_TAIT

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine EOS_tait_ADIABATIC_rho



      subroutine SOUNDSQR_tait_rho(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,rhobar,pcav,rhocav,pressure
      REAL_T rho_sound


      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in soundsqr tait rho"
       stop
      endif

      pcav=PCAV_TAIT
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (rhocav.le.zero) then
       print *,"rhocav invalid"
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_TAIT*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_TAIT-one) )

      return
      end subroutine SOUNDSQR_tait_rho



      subroutine INTERNAL_tait_rhohydro(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=temperature*cv

      return
      end subroutine INTERNAL_tait_rhohydro

      subroutine TEMPERATURE_tait_rhohydro(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_tait_rhohydro

      subroutine EOS_tait_rhohydro(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,rhobar,pcav


      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos tait rho"
       stop
      endif

      pcav=PCAV_TAIT

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine EOS_tait_rhohydro


      subroutine SOUNDSQR_tait_rhohydro(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,rhobar,pcav,rhocav,pressure
      REAL_T rho_sound


      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in soundsqr tait rho"
       stop
      endif

      pcav=PCAV_TAIT
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (rhocav.le.zero) then
       print *,"rhocav invalid"
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_TAIT*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_TAIT-one) )

      return
      end subroutine SOUNDSQR_tait_rhohydro




      subroutine INTERNAL_tait_rho3(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=temperature*cv

      return
      end subroutine INTERNAL_tait_rho3

      subroutine TEMPERATURE_tait_rho3(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_tait_rho3

      subroutine EOS_tait_rho3(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,rhobar,pcav


      if (num_materials.lt.3) then
       print *,"num materials must be at least 3"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(3) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos tait rho"
       stop
      endif

      pcav=PCAV_TAIT

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine EOS_tait_rho3


      subroutine EOS_tait_ADIABATIC_rho3(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,pressure
      REAL_T A,B,rhobar,pcav


      if (num_materials.lt.3) then
       print *,"num materials must be at least 3"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(3) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos tait adiabatic rho"
       stop
      endif

      pcav=PCAV_TAIT

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine EOS_tait_ADIABATIC_rho3



      subroutine SOUNDSQR_tait_rho3(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,rhobar,pcav,rhocav,pressure
      REAL_T rho_sound


      if (num_materials.lt.3) then
       print *,"num materials must be at least 3"
       stop
      endif

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(3) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in soundsqr tait rho"
       stop
      endif

      pcav=PCAV_TAIT  ! dyne/cm^2
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (rhocav.le.zero) then
       print *,"rhocav invalid"
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_TAIT*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_TAIT-one) )

      return
      end subroutine SOUNDSQR_tait_rho3


! --------------- NEW EOS STUFF 10,11,12 --------------------


      subroutine INTERNAL_tait_rho2(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=temperature*cv

      return
      end subroutine INTERNAL_tait_rho2

      subroutine TEMPERATURE_tait_rho2(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_tait_rho2

      subroutine EOS_tait_rho2(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,rhobar,pcav


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(2) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos tait rho"
       stop
      endif

      pcav=PCAV_TAIT  ! dyne/cm^2

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine EOS_tait_rho2


      subroutine EOS_tait_ADIABATIC_rho2(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,pressure
      REAL_T A,B,rhobar,pcav


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(2) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos tait adiabatic rho"
       stop
      endif

      pcav=PCAV_TAIT  ! dyne/cm^2

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine EOS_tait_ADIABATIC_rho2



      subroutine SOUNDSQR_tait_rho2(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,rhobar,pcav,rhocav,pressure
      REAL_T rho_sound


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(2) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in soundsqr tait rho"
       stop
      endif

      pcav=PCAV_TAIT  ! dyne/cm^2
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (rhocav.le.zero) then
       print *,"rhocav invalid"
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_TAIT*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_TAIT-one) )

      return
      end subroutine SOUNDSQR_tait_rho2



      subroutine INTERNAL_koren_rho2(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=temperature*cv

      return
      end subroutine

      subroutine TEMPERATURE_koren_rho2(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine

      subroutine EOS_koren_rho2(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,rhobar,pcav,GAMMA_KOREN


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif
! Koren:
! p=(rho/rhoref)^gamma (1+B)pref -B pref
!  = (1+B)pref ( (rho/rhoref)^gamma - 1 ) + pref
! So: A -> pref
!     B -> (1+B)pref

      A=1.0
      B=1.0 
      rhobar=fort_denconst(2) ! g/cm^3

      if (rhobar.lt.0.0001) then
       print *,"rhobar invalid in eos koren rho2"
       stop
      endif

      GAMMA_KOREN=1.4
      pcav=0.001  ! dyne/cm^2

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

! Koren:
! p=(rho/rhoref)^gamma (1+B)pref -B pref
!  = (1+B)pref ( (rho/rhoref)^gamma - 1 ) + pref
! So: A -> pref
!     B -> (1+B)pref

      pressure=B*( (rho/rhobar)**GAMMA_KOREN - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine


      subroutine EOS_koren_ADIABATIC_rho2(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,pressure
      REAL_T A,B,rhobar,GAMMA_KOREN,pcav


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif

      A=1.0
      B=1.0
      rhobar=fort_denconst(2) ! g/cm^3

      if (rhobar.lt.0.0001) then
       print *,"rhobar invalid in eos koren adiabatic rho"
       stop
      endif

      GAMMA_KOREN=1.4
      pcav=0.001

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_KOREN - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine



      subroutine SOUNDSQR_koren_rho2(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,rhobar,GAMMA_KOREN,pcav,rhocav,pressure
      REAL_T rho_sound


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      A=1.0
      B=1.0
      rhobar=fort_denconst(2) ! g/cm^3

      if (rhobar.lt.0.0001) then
       print *,"rhobar invalid in soundsqr koren rho"
       stop
      endif

      GAMMA_KOREN=1.4
      pcav=0.001
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_KOREN) )
      pressure=B*( (rho/rhobar)**GAMMA_KOREN - one ) + A
      
      if (rhocav.le.zero) then
       print *,"rhocav invalid"
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_KOREN*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_KOREN-one) )

      return
      end subroutine


      subroutine INTERNAL_koren_rho1(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature <=0"
       stop
      endif
      cv=4.1855D+7
      internal_energy=temperature*cv

      return
      end subroutine

      subroutine TEMPERATURE_koren_rho1(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,cv


      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy <=0"
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine

      subroutine EOS_koren_rho1(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure
      REAL_T A,B,rhobar,GAMMA_KOREN,pcav


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif
! Koren:
! p=(rho/rhoref)^gamma (1+B)pref -B pref
!  = (1+B)pref ( (rho/rhoref)^gamma - 1 ) + pref
! So: A -> pref
!     B -> (1+B)pref

      A=1.0
      B=3001.0
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos koren rho"
       stop
      endif

      GAMMA_KOREN=7.0
      pcav=0.001  ! dyne/cm^2

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

! Koren:
! p=(rho/rhoref)^gamma (1+B)pref -B pref
!  = (1+B)pref ( (rho/rhoref)^gamma - 1 ) + pref
! So: A -> pref
!     B -> (1+B)pref

      pressure=B*( (rho/rhobar)**GAMMA_KOREN - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine


      subroutine EOS_koren_ADIABATIC_rho1(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,pressure
      REAL_T A,B,rhobar,GAMMA_KOREN,pcav


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif

      A=1.0
      B=3001.0
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos koren adiabatic rho"
       stop
      endif

      GAMMA_KOREN=7.0
      pcav=0.001

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_KOREN - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif
      

      return
      end subroutine



      subroutine SOUNDSQR_koren_rho1(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,soundsqr
      REAL_T A,B,rhobar,GAMMA_KOREN,pcav,rhocav,pressure
      REAL_T rho_sound


      if (num_materials.lt.2) then
       print *,"num materials must be at least 2"
       stop
      endif

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      A=1.0
      B=3001.0
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in soundsqr koren rho"
       stop
      endif

      GAMMA_KOREN=7.0
      pcav=0.001
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_KOREN) )
      pressure=B*( (rho/rhobar)**GAMMA_KOREN - one ) + A
      
      if (rhocav.le.zero) then
       print *,"rhocav invalid"
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_KOREN*B/rhobar)*( (rho_sound/rhobar)**(GAMMA_KOREN-one) )

      return
      end subroutine SOUNDSQR_koren_rho1



! -------------------- END OF EOS 10,11,12 ----------------


       ! BENCHMARK TESTS ONLY REPORT ENERGY, PRESSURE, DENSITY, and VELOCITY,
       ! THESE tests are inviscid and have zero thermal diffusion,
       ! therefore not important what cp or cv are so just set them
       ! to 1 for (a) inputs.shockturbulence, (b) Gary Sod test problem
       ! inputs.sod, (c) inputs.strong, (d) inputs.mach4
      subroutine simple_air_parms(R,cp,cv,gamma_constant,omega)
      use probcommon_module
      IMPLICIT NONE
      REAL_T R,cp,cv,gamma_constant,omega

      cv=one
      gamma_constant=GAMMA_SIMPLE_PARMS
      cp=gamma_constant*cv
      R=cp-cv
      omega=gamma_constant-one

      return
      end subroutine simple_air_parms


      subroutine air_parms(R,cp,cv,gamma_constant,omega)
      use probcommon_module
      IMPLICIT NONE
      REAL_T, intent(out) :: R,cp,cv,gamma_constant,omega

       ! PARAMETERS declared in PROBCOMMON.F90
      R=R_AIR_PARMS  ! ergs/(Kelvin g)
      cv=CV_AIR_PARMS ! ergs/(Kelvin g)
      cp=cv+R  ! ergs/(Kelvin g)
      gamma_constant=cp/cv
      omega=gamma_constant-one

      return
      end subroutine air_parms

       ! density_at_depth previously initialized by:
       ! init_density_at_depth() 
       ! called from:  
       ! boundary_hydrostatic, EOS_air_rho2, EOS_air_rho2_ADIABAT,
       ! SOUNDSQR_air_rho2, EOS_error_ind, presBDRYCOND, FORT_INITDATA 
      subroutine tait_hydrostatic_pressure_density( &
        xpos,rho,pres,from_boundary_hydrostatic)
      use probcommon_module
      IMPLICIT NONE

      REAL_T, intent(in) :: xpos(SDIM)
      REAL_T, intent(inout) :: rho
      REAL_T, intent(inout) :: pres
      INTEGER_T, intent(in) :: from_boundary_hydrostatic
      REAL_T denfree,zfree
      REAL_T z_at_depth

      ! in tait_hydrostatic_pressure_density
      if ((probtype.eq.36).and.(axis_dir.eq.2)) then  ! spherical explosion
       rho=one
       call EOS_tait_ADIABATIC(rho,pres)
      ! JICF nozzle+pressure bc
      else if ((probtype.eq.53).and.(axis_dir.eq.2)) then
       rho=one
       call EOS_tait_ADIABATIC(rho,pres)
       ! JICF
      else if ((probtype.eq.53).and.(fort_material_type(1).eq.7)) then
       rho=fort_denconst(1)
       call EOS_tait_ADIABATIC_rho(rho,pres)
      ! impinging jets
      else if ((probtype.eq.530).and.(axis_dir.eq.1).and. &
               (fort_material_type(1).eq.7).and.(SDIM.eq.3)) then
       rho=fort_denconst(1)
       call EOS_tait_ADIABATIC_rho(rho,pres)

       ! in: tait_hydrostatic_pressure_density
      else if ((probtype.eq.42).and.(SDIM.eq.2)) then  ! bubble jetting

       if (probloy.ne.zero) then
        print *,"probloy must be 0 for bubble jetting problem"
        stop
       endif
       ! yblob is distance from domain bottom of charge
       ! zblob is depth of charge
       denfree=one
       zfree=zblob+yblob  ! relative to computational grid
       z_at_depth=yblob
       if (xpos(SDIM).gt.zfree) then
        rho=denfree
       else
        rho= &
          ((density_at_depth-denfree)/ &
           (z_at_depth-zfree))*(xpos(SDIM)-zfree)+denfree
       endif
       call EOS_tait_ADIABATIC(rho,pres)

      else if ((probtype.eq.46).and.(SDIM.eq.2)) then  ! cavitation

       if (probloy.ne.zero) then
        print *,"probloy must be 0 for cavitation problem"
        stop
       endif
       ! yblob is distance from domain bottom of charge/sphere
       ! zblob is depth of charge (for jwl problem)
       denfree=one
       if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
        zfree=zblob+yblob  ! relative to computational grid
        z_at_depth=yblob
       else if (axis_dir.eq.10) then
        zfree=zblob
        z_at_depth=zero
       else if (axis_dir.eq.20) then
        print *,"there is no gravity for the CODY ESTEBE created test problem"
        stop
       else
        print *,"axis_dir out of range"
        stop
       endif
       if (xpos(SDIM).gt.zfree) then
        rho=denfree
       else
        rho= &
          ((density_at_depth-denfree)/ &
           (z_at_depth-zfree))*(xpos(SDIM)-zfree)+denfree
       endif
       call EOS_tait_ADIABATIC(rho,pres)

      else if (fort_material_type(1).eq.13) then

       denfree=fort_denconst(1)
       if (SDIM.eq.2) then
        zfree=probhiy
        z_at_depth=probloy
       else if (SDIM.eq.3) then
        zfree=probhiz
        z_at_depth=probloz
       else
        print *,"dimension bust"
        stop
       endif

        ! density_at_depth is found so that
        ! (p(density_at_depth)-p(rho_0))/(rho_0 (z_at_depth-zfree))=g
        !
       if (xpos(SDIM).gt.zfree) then
        rho=denfree
       else
        rho= &
          ((density_at_depth-denfree)/ &
           (z_at_depth-zfree))*(xpos(SDIM)-zfree)+denfree
       endif
       call EOS_tait_ADIABATIC_rhohydro(rho,pres)
      else
       print *,"probtype invalid tait_hydrostatic_pressure_density"
       stop
      endif

      return
      end subroutine tait_hydrostatic_pressure_density


      subroutine EOS_air_rho2(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure,omega
      REAL_T cp,cv,R,pressure_adjust,preshydro,rhohydro
      REAL_T xpos(SDIM)
      INTEGER_T from_boundary_hydrostatic

      from_boundary_hydrostatic=0

      call air_parms(R,cp,cv,gamma_constant,omega) 
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
       ! (gamma-1)rho*cv T=(cp/cv -1)*rho*cv T=
       ! (cp-cv)*rho*T=(R_universal/MolarMass)*rho*T
       ! R_universal ergs/(mol K)
       ! MolarMass g/mol
       ! erg=g cm^2/s^2
       ! R_universal/MolarMass=g cm^2 / (s^2 g K)  
       ! (g cm^2 / (s^2 g K)) *( g/cm^3 )  * K = (g/cm)/s^2
       ! pressure=dyne/m^2=g (cm/s^2 )/cm^2 = (g/cm)/s^2
      pressure_adjust=omega*fort_denconst(2)* &
         cv*fort_tempconst(2)

       ! first material uses "EOS_tait_rhohydro" ?
      if (fort_material_type(1).eq.13) then
       if (SDIM.eq.2) then
        xpos(SDIM)=yblob
       else if (SDIM.eq.3) then
        xpos(SDIM)=zblob
       else
        print *,"dimension bust"
        stop
       endif
       call tait_hydrostatic_pressure_density(xpos,rhohydro,preshydro, &
               from_boundary_hydrostatic)
      else
       preshydro=1.0D+6
      endif
      pressure_adjust=preshydro/pressure_adjust

      pressure=omega*rho*internal_energy*pressure_adjust

      return
      end subroutine EOS_air_rho2



      subroutine EOS_air_rho2_ADIABAT(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,gamma_constant,pressure
      REAL_T cp,cv,R,rhohydro,omega
      REAL_T RHOI,PI
      REAL_T xpos(SDIM)
      INTEGER_T from_boundary_hydrostatic

      from_boundary_hydrostatic=0

      RHOI=fort_denconst(2)
      call general_hydrostatic_pressure(PI)
    
      call air_parms(R,cp,cv,gamma_constant,omega)
      if (RHOI.le.zero) then
       print *,"RHOI invalid"
       stop
      endif
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      gamma_constant=cp/cv
      if (fort_material_type(1).eq.13) then
       if (SDIM.eq.2) then
        xpos(SDIM)=yblob
       else if (SDIM.eq.3) then
        xpos(SDIM)=zblob
       else
        print *,"dimension bust"
        stop
       endif
       call tait_hydrostatic_pressure_density(xpos,rhohydro,PI, &
               from_boundary_hydrostatic)
      endif

      pressure=PI*((rho/RHOI)**gamma_constant)

      return
      end subroutine EOS_air_rho2_ADIABAT


      subroutine ENTROPY_air_rho2(rho,internal_energy,entropy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure,entropy,press_adiabat

      call EOS_air_rho2(rho,internal_energy,pressure)
      call EOS_air_rho2_ADIABAT(rho,press_adiabat)
      if (press_adiabat.le.zero) then
       print *,"press_adiabat invalid"
       stop
      endif
      entropy=pressure/press_adiabat

      return
      end subroutine ENTROPY_air_rho2

      subroutine INTERNAL_ENTROPY_air_rho2(rho,entropy,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,entropy,internal_energy
      REAL_T press_adiabat,unit_internal_energy,unit_press

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (entropy.le.zero) then
       print *,"entropy invalid"
       stop
      endif
      call EOS_air_rho2_ADIABAT(rho,press_adiabat)
      unit_internal_energy=one
      call EOS_air_rho2(rho,unit_internal_energy,unit_press)
      internal_energy=press_adiabat*entropy/unit_press
      if (internal_energy.le.zero) then
       print *,"internal_energy invalid"
       stop
      endif

      return
      end subroutine INTERNAL_ENTROPY_air_rho2




      subroutine SOUNDSQR_air_rho2(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure,omega
      REAL_T soundsqr
      REAL_T cp,cv,R,pressure_adjust,preshydro,rhohydro
      REAL_T xpos(SDIM)
      INTEGER_T from_boundary_hydrostatic

      from_boundary_hydrostatic=0

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure_adjust=omega*fort_denconst(2)* &
         cv*fort_tempconst(2)
      if (fort_material_type(1).eq.13) then
       if (SDIM.eq.2) then
        xpos(SDIM)=yblob
       else if (SDIM.eq.3) then
        xpos(SDIM)=zblob
       else
        print *,"dimension bust"
        stop
       endif
       call tait_hydrostatic_pressure_density(xpos,rhohydro,preshydro, &
               from_boundary_hydrostatic)
      else
       preshydro=1.0D+6
      endif
      pressure_adjust=preshydro/pressure_adjust
      pressure=omega*rho*internal_energy*pressure_adjust
      soundsqr=gamma_constant*pressure/rho

      return
      end subroutine SOUNDSQR_air_rho2


      subroutine INTERNAL_air_rho2(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0 in internal_air_rho2"
       print *,temperature, rho, R, cp, cv
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      internal_energy=cv*temperature

      return
      end subroutine



      subroutine TEMPERATURE_air_rho2(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_air_rho2



! e=c_v T
! gamma=cp/cv
! p=(cp/cv-1)rho e=(cp-cv)rho e/cv=(cp-cv)rho T  R=cp-cv
! 1 erg=1 dyne cm=1 g cm^2/s^2 = 10^{-7} Joules
! 1 J/(g K)=10^7 erg/(g K)
! p has units dyne/cm^2=g (cm/s^2)/cm^2=g/(s^2 cm)
! T has units of Kelvin
! cp,cv have units erg/(g K)= (cm^2/s^2)/K
! (cp-cv) rho T units=( (cm^2/s^2)/K ) (g/cm^3) K=g/(cm s^2)
!
! for air, R=287.0 m^2/(s^2 K)=0.287D+7 cm^2/(s^2 K)
! cv=0.72D+7 cm^2/(s^2 K)
!
! for SF6, R=56.92 m^2/(s^2 K)=0.05692D+7 cm^2/(s^2 K)


      subroutine EOS_air(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
        ! e=cv T
        ! rho=g/(cm^3)
        ! internal_energy has units ergs/g
        ! pressure has units of g/(cm^3)  ergs/g =ergs/cm^3
        ! 1erg=1 g cm^2/s^2
        ! erg/cm^3=(g cm^2/s^2)/cm^3 = g/(s^2 cm)
        ! pressure also described as dyne/cm^2 = (g cm/s^2)/cm^2=g/(cm s^2) 
      pressure=omega*rho*internal_energy  ! omega=gamma-1

      return
      end subroutine EOS_air


      subroutine EOS_simple_air(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T cp,cv,R,omega

      call simple_air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=omega*rho*internal_energy

      return
      end subroutine EOS_simple_air


      subroutine EOS_air_ADIABAT(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,gamma_constant,pressure
      REAL_T cp,cv,R,RHOI,PI,omega

      RHOI=fort_denconst(2)
        ! PI=1.0D+6
      call general_hydrostatic_pressure(PI)
        ! gamma=1.399
      call air_parms(R,cp,cv,gamma_constant,omega)

      if (RHOI.le.zero) then
       print *,"RHOI invalid"
       stop
      endif
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=PI*(rho/RHOI)**gamma_constant

      return
      end subroutine EOS_air_ADIABAT

      subroutine get_rhocav(rho)
      use probcommon_module
      IMPLICIT NONE
       
      REAL_T rho
      REAL_T term1,term2

        ! p=B((rho/rho0)^gamma - 1 ) + A
        ! rho=rho0(((p-A)/B+1)^(1/gamma))
      term1=(PCAV_TAIT-A_TAIT)/B_TAIT+one
      term2=term1**(one/GAMMA_TAIT)
      rho=RHOBAR_TAIT*term2

      return
      end subroutine get_rhocav

      subroutine EOS_vacuum(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,gamma_constant,pressure
      REAL_T cp,cv,R,RHOI,PI,omega

      call get_rhocav(RHOI)
      PI=PCAV_TAIT
      call air_parms(R,cp,cv,gamma_constant,omega)

      if (RHOI.le.zero) then
       print *,"RHOI invalid"
       stop
      endif
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=PI*(rho/RHOI)**gamma_constant

      return
      end subroutine EOS_vacuum



      subroutine ENTROPY_air(rho,internal_energy,entropy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure,entropy,press_adiabat

      call EOS_air(rho,internal_energy,pressure)
      call EOS_air_ADIABAT(rho,press_adiabat)
      if (press_adiabat.le.zero) then
       print *,"press_adiabat invalid"
       stop
      endif
      entropy=pressure/press_adiabat

      return
      end subroutine ENTROPY_air

      subroutine INTERNAL_ENTROPY_air(rho,entropy,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,entropy,internal_energy
      REAL_T press_adiabat,unit_internal_energy,unit_press

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (entropy.le.zero) then
       print *,"entropy invalid"
       stop
      endif
      call EOS_air_ADIABAT(rho,press_adiabat)
      unit_internal_energy=one
      call EOS_air(rho,unit_internal_energy,unit_press)
      internal_energy=press_adiabat*entropy/unit_press
      if (internal_energy.le.zero) then
       print *,"internal_energy invalid"
       stop
      endif

      return
      end subroutine INTERNAL_ENTROPY_air




         ! He+28 percent air

      subroutine EOS_Marquina(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T cp,cv,R

    
      R=1.578D+7  
      cv=2.44D+7  

      cp=cv+R 
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      gamma_constant=cp/cv
      pressure=(gamma_constant-one)*rho*internal_energy

      return
      end subroutine

! e=(E/rho) - (1/2) (u^2 + v^2)
! c^2=dp/drho + (p/rho^2) dp/de
! p=(gamma-1)rho e
! c^2=(gamma-1)e+((gamma-1)rho e/rho^2)(gamma-1)rho=
!   (gamma-1)e+(gamma-1) e (gamma-1)=(p/rho)(1+(gamma-1))=gamma p/rho
!


      subroutine SOUNDSQR_air(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T soundsqr
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=omega*rho*internal_energy
      soundsqr=gamma_constant*pressure/rho

      return
      end subroutine


      subroutine SOUNDSQR_simple_air(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T soundsqr
      REAL_T cp,cv,R,omega

      call simple_air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=omega*rho*internal_energy
      soundsqr=gamma_constant*pressure/rho

      return
      end subroutine


 
         ! He+28 percent air

      subroutine SOUNDSQR_Marquina(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T soundsqr
      REAL_T cp,cv,R

    
      R=1.578D+7  
      cv=2.44D+7  

      cp=cv+R 
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      gamma_constant=cp/cv
      pressure=(gamma_constant-one)*rho*internal_energy
      soundsqr=gamma_constant*pressure/rho

      return
      end subroutine


      subroutine INTERNAL_air(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

        ! cv has units of ergs/(Kelvin g)
        ! internal energy has units of ergs/g
      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_air


      subroutine INTERNAL_simple_air(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call simple_air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      internal_energy=cv*temperature

      return
      end subroutine

       ! He+28 percent air

      subroutine INTERNAL_Marquina(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R

    
      R=1.578D+7  
      cv=2.44D+7  

      cp=cv+R 
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      gamma_constant=cp/cv

      internal_energy=cv*temperature

      return
      end subroutine


      subroutine TEMPERATURE_air(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

    
      call air_parms(R,cp,cv,gamma_constant,omega)
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine


      subroutine TEMPERATURE_simple_air(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

    
      call simple_air_parms(R,cp,cv,gamma_constant,omega)
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine


        ! He+28 percent air

      subroutine TEMPERATURE_Marquina(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R

    
      R=1.578D+7  
      cv=2.44D+7  

      cp=cv+R 
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      gamma_constant=cp/cv

      temperature=internal_energy/cv

      return
      end subroutine


! e=(E/rho) - (1/2) (u^2 + v^2)
! c^2=dp/drho + (p/rho^2) dp/de
! p=(gamma-1)rho e
! c^2=(gamma-1)e+((gamma-1)rho e/rho^2)(gamma-1)rho=
!   (gamma-1)e+(gamma-1) e (gamma-1)=(p/rho)(1+(gamma-1))=gamma p/rho
!

      subroutine SOUNDSQR_airADIABAT(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T soundsqr
      REAL_T cp,cv,R,phyd,rho0,omega


      phyd=1.0D+6 
      rho0=fort_denconst(2)

      if (rho0.le.zero) then
       print *,"rho0 negative"
       stop
      endif
      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=phyd*((rho/rho0)**gamma_constant)
      soundsqr=gamma_constant*pressure/rho

      return
      end subroutine SOUNDSQR_airADIABAT


      subroutine SOUNDSQR_vacuum(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T soundsqr
      REAL_T cp,cv,R,phyd,rho0,omega


      phyd=PCAV_TAIT
      call get_rhocav(rho0)

      if (rho0.le.zero) then
       print *,"rho0 negative"
       stop
      endif
      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=phyd*((rho/rho0)**gamma_constant)
      soundsqr=gamma_constant*pressure/rho

      return
      end subroutine SOUNDSQR_vacuum


 
      subroutine INTERNAL_airADIABAT(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_airADIABAT


      subroutine INTERNAL_vacuum(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_vacuum

      subroutine TEMPERATURE_airADIABAT(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_airADIABAT


      subroutine TEMPERATURE_vacuum(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_vacuum

! e=c_v T
! gamma=cp/cv
! p=(cp/cv-1)rho e=(cp-cv)rho e/cv=(cp-cv)rho T  R=cp-cv
! 1 erg=1 dyne cm=1 g cm^2/s^2 = 10^{-7} Joules
! 1 J/(g K)=10^7 erg/(g K)
! p has units dyne/cm^2=g (cm/s^2)/cm^2=g/(s^2 cm)
! T has units of Kelvin
! cp,cv have units erg/(g K)= (cm^2/s^2)/K
! (cp-cv) rho T units=( (cm^2/s^2)/K ) (g/cm^3) K=g/(cm s^2)
!
! for air, R=287.0 m^2/(s^2 K)=0.287D+7 cm^2/(s^2 K)
! cv=0.72D+7 cm^2/(s^2 K)
!
! for SF6, R=56.92 m^2/(s^2 K)=0.05692D+7 cm^2/(s^2 K)
! cp=0.097 kJ/(mol K)=0.097 1000/146.05 J/(g K)=
! 0.097 1000/146.05 10^7  erg/(g K)=0.664D+7 cm^2/(s^2 K)
! R=cp-cv
! cv=cp-R=(0.664-0.05692)=0.607D+7 cm^2/(s^2 K)
! gamma=cp/cv=0.664/0.607=1.09

      subroutine SF6_parms(R,cp,cv,gamma_constant,omega)
      use probcommon_module
      IMPLICIT NONE
      REAL_T R,cp,cv,gamma_constant,omega

      R=0.05692D+7
      cp=0.664D+7
      cv=0.607D+7
      gamma_constant=cp/cv
      omega=gamma_constant-one

      return
      end subroutine SF6_parms

      subroutine EOS_SF6(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T cp,cv,R,omega

      call SF6_parms(R,cp,cv,gamma_constant,omega)
       
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=omega*rho*internal_energy

      return
      end subroutine EOS_SF6


      subroutine EOS_SF6_ADIABAT(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,gamma_constant,pressure
      REAL_T cp,cv,R,RHOI,PI,omega

      RHOI=fort_denconst(2)
      call general_hydrostatic_pressure(PI)
      call SF6_parms(R,cp,cv,gamma_constant,omega)
       
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=PI*(rho/RHOI)**gamma_constant

      return
      end subroutine EOS_SF6_ADIABAT

      subroutine ENTROPY_SF6(rho,internal_energy,entropy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,pressure,entropy,press_adiabat

      call EOS_SF6(rho,internal_energy,pressure)
      call EOS_SF6_ADIABAT(rho,press_adiabat)
      if (press_adiabat.le.zero) then
       print *,"press_adiabat invalid"
       stop
      endif
      entropy=pressure/press_adiabat

      return
      end subroutine ENTROPY_SF6

      subroutine INTERNAL_ENTROPY_SF6(rho,entropy,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,entropy,internal_energy
      REAL_T press_adiabat,unit_internal_energy,unit_press

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (entropy.le.zero) then
       print *,"entropy invalid"
       stop
      endif
      call EOS_SF6_ADIABAT(rho,press_adiabat)
      unit_internal_energy=one
      call EOS_SF6(rho,unit_internal_energy,unit_press)
      internal_energy=press_adiabat*entropy/unit_press
      if (internal_energy.le.zero) then
       print *,"internal_energy invalid"
       stop
      endif

      return
      end subroutine INTERNAL_ENTROPY_SF6

      subroutine TEMPERATURE_ENTROPY_SF6(rho,entropy,temperature)
      use probcommon_module
      IMPLICIT NONE
      
      REAL_T R,cp,cv,gamma_constant,omega 
      REAL_T rho,entropy,temperature,RHOI

      call SF6_parms(R,cp,cv,gamma_constant,omega)
      RHOI=fort_denconst(2)
      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (entropy.le.zero) then
       print *,"entropy invalid"
       stop     
      endif    
      temperature=entropy*((rho/RHOI)**omega)
      if (temperature.le.zero) then
       print *,"temperature invalid"
       stop
      endif
               
      return   
      end subroutine TEMPERATURE_ENTROPY_SF6

      subroutine ENTROPY_TEMPERATURE_SF6(rho,temperature,entropy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T R,cp,cv,gamma_constant,omega 
      REAL_T rho,temperature,entropy,RHOI

      call SF6_parms(R,cp,cv,gamma_constant,omega)
      RHOI=fort_denconst(2)
      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature invalid"
       stop     
      endif    
      entropy=temperature/((rho/RHOI)**omega)
      if (entropy.le.zero) then
       print *,"entropy invalid"
       stop
      endif

      return
      end subroutine ENTROPY_TEMPERATURE_SF6


! e=(E/rho) - (1/2) (u^2 + v^2)
! c^2=dp/drho + (p/rho^2) dp/de
! p=(gamma-1)rho e
! c^2=(gamma-1)e+((gamma-1)rho e/rho^2)(gamma-1)rho=
!   (gamma-1)e+(gamma-1) e (gamma-1)=(p/rho)(1+(gamma-1))=gamma p/rho
! 
      subroutine SOUNDSQR_SF6(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T soundsqr
      REAL_T cp,cv,R,omega

      call SF6_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=omega*rho*internal_energy
      soundsqr=gamma_constant*pressure/rho

      return
      end subroutine


      subroutine INTERNAL_SF6(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call SF6_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      internal_energy=cv*temperature

      return
      end subroutine


      subroutine TEMPERATURE_SF6(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call SF6_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine


      subroutine EOS_stiffened(rho,internal_energy,pressure,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T im
      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T cp,cv,PP

      if (im.lt.1) then
       print *,"im invalid65"
       stop
      endif
      cp=get_user_stiffCP(im)
      if (cp.le.zero) then
       print *,"cp invalid"
       stop
      endif
      gamma_constant=fort_stiffGAMMA(im)
      if (gamma_constant.le.zero) then
       print *,"gamma_constant invalid"
       stop
      endif
      cv=cp/gamma_constant

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      PP=fort_stiffPINF(im)
      if (PP.lt.zero) then
       print *,"fort_stiff PINF invalid"
       stop
      endif
      pressure=(gamma_constant-one)*rho*internal_energy- &
       gamma_constant*PP
      if (pressure.lt.VOFTOL) then
       pressure=VOFTOL
      endif

      return
      end subroutine

      subroutine SOUNDSQR_stiffened(rho,internal_energy,soundsqr,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T im
      REAL_T rho,internal_energy,gamma_constant,pressure
      REAL_T soundsqr
      REAL_T cp,cv,PP

      if (im.lt.1) then
       print *,"im invalid66"
       stop
      endif
      cp=get_user_stiffCP(im)
      if (cp.le.zero) then
       print *,"cp invalid"
       stop
      endif
      gamma_constant=fort_stiffGAMMA(im)
      if (gamma_constant.le.zero) then
       print *,"gamma_constant invalid"
       stop
      endif
      cv=cp/gamma_constant

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      PP=fort_stiffPINF(im)
      if (PP.lt.zero) then
       print *,"fort_stiff PINF invalid"
       stop
      endif

      pressure=(gamma_constant-one)*rho*internal_energy- &
       gamma_constant*PP
      if (pressure.lt.VOFTOL) then
       pressure=VOFTOL
      endif

      soundsqr=gamma_constant*(pressure+PP)/rho

      return
      end subroutine


      subroutine INTERNAL_stiffened(rho,temperature,internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T im
      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,PP

      if (im.lt.1) then
       print *,"im invalid67"
       stop
      endif
      cp=get_user_stiffCP(im)
      if (cp.le.zero) then
       print *,"cp invalid"
       stop
      endif
      gamma_constant=fort_stiffGAMMA(im)
      if (gamma_constant.le.zero) then
       print *,"gamma_constant invalid"
       stop
      endif
      cv=cp/gamma_constant

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      PP=fort_stiffPINF(im)
      if (PP.lt.zero) then
       print *,"fort_stiff PINF invalid"
       stop
      endif

      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif

      internal_energy=cv*temperature+PP/rho

      return
      end subroutine INTERNAL_stiffened


      subroutine TEMPERATURE_stiffened(rho,temperature, &
        internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T im
      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,PP

      if (im.lt.1) then
       print *,"im invalid68"
       stop
      endif
      cp=get_user_stiffCP(im)
      if (cp.le.zero) then
       print *,"cp invalid"
       stop
      endif
      gamma_constant=fort_stiffGAMMA(im)
      if (gamma_constant.le.zero) then
       print *,"gamma_constant invalid"
       stop
      endif
      cv=cp/gamma_constant

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      PP=fort_stiffPINF(im)
      if (PP.lt.zero) then
       print *,"fort_stiff PINF invalid"
       stop
      endif

      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif

      temperature=(internal_energy-PP/rho)/cv
      if (temperature.lt.VOFTOL) then
       temperature=VOFTOL
      endif

      return
      end subroutine TEMPERATURE_stiffened



! e=(E/rho) - (1/2) (u^2 + v^2)
! c^2=dp/drho + (p/rho^2) dp/de
! p=(gamma-1)rho e
! c^2=(gamma-1)e+((gamma-1)rho e/rho^2)(gamma-1)rho=
!   (gamma-1)e+(gamma-1) e (gamma-1)=(p/rho)(1+(gamma-1))=gamma p/rho
! 
      subroutine SOUNDSQR_SF6ADIABAT(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,internal_energy,gamma_constant,pressure,omega
      REAL_T soundsqr
      REAL_T cp,cv,R


      call SF6_parms(R,cp,cv,gamma_constant,omega)

      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      call EOS_SF6_ADIABAT(rho,pressure)
      soundsqr=gamma_constant*pressure/rho

      return
      end subroutine SOUNDSQR_SF6ADIABAT


      subroutine INTERNAL_SF6ADIABAT(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call SF6_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      internal_energy=cv*temperature

      return
      end subroutine


      subroutine TEMPERATURE_SF6ADIABAT(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      REAL_T rho,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call SF6_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_SF6ADIABAT

      subroutine init_massfrac_parm(den,massfrac_parm,im)
      use probcommon_module
      IMPLICIT NONE

      REAL_T, intent(in) :: den
      REAL_T, intent(out) :: massfrac_parm(num_species_var+1)
      INTEGER_T, intent(in) :: im
      INTEGER_T :: ispec,var_comp

      if (num_species_var.eq.0) then
       massfrac_parm(1)=den
      else if (num_species_var.ge.1) then
       do ispec=1,num_species_var
        var_comp=(ispec-1)*num_materials+im
        massfrac_parm(ispec)=fort_speciesconst(var_comp)
       enddo
      else
       print *,"num_species_var invalid"
       stop
      endif

      return
      end subroutine init_massfrac_parm


      subroutine EOS_material_CORE(rho,massfrac_var, &
        internal_energy,pressure, &
        imattype,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho,internal_energy
      REAL_T, intent(in) :: massfrac_var(num_species_var+1)
      REAL_T, intent(inout) :: pressure
      REAL_T :: T

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      if (imattype.eq.1) then
       call EOS_tait(rho,internal_energy,pressure)
      else if (imattype.eq.2) then
       call EOS_jwlADIABAT(rho,internal_energy,pressure)
      else if (imattype.eq.3) then
       call EOS_NAjwl(rho,internal_energy,pressure)
      else if (imattype.eq.4) then
       call EOS_SF6(rho,internal_energy,pressure)
      else if (imattype.eq.5) then
       call EOS_air(rho,internal_energy,pressure)
      else if (imattype.eq.14) then
       call EOS_air_rho2(rho,internal_energy,pressure)
      else if (imattype.eq.6) then
       call EOS_Marquina(rho,internal_energy,pressure)
      else if (imattype.eq.7) then
       call EOS_tait_rho(rho,internal_energy,pressure)
      else if (imattype.eq.13) then
       call EOS_tait_rhohydro(rho,internal_energy,pressure)
      else if (imattype.eq.8) then
       call EOS_air_ADIABAT(rho,pressure)
      else if (imattype.eq.9) then
       call EOS_tait_rho3(rho,internal_energy,pressure)
      else if (imattype.eq.10) then
       call EOS_tait_rho2(rho,internal_energy,pressure)
      else if (imattype.eq.11) then
       call EOS_koren_rho1(rho,internal_energy,pressure)
      else if (imattype.eq.12) then
       call EOS_koren_rho2(rho,internal_energy,pressure)
      else if (imattype.eq.15) then
       call EOS_dodecane(rho,internal_energy,T,pressure)
      else if (imattype.eq.16) then
       call EOS_SF6_ADIABAT(rho,pressure)
      else if (imattype.eq.17) then
       call EOS_stiffened(rho,internal_energy,pressure,im)
      else if (imattype.eq.18) then
       call EOS_simple_air(rho,internal_energy,pressure)
      else if (imattype.eq.19) then
       call EOS_vacuum(rho,pressure)
      else if (imattype.eq.20) then
       call EOS_tait_vacuum(rho,pressure)
      else if (imattype.eq.21) then
       call EOS_cav_nozzle(rho,internal_energy,pressure)
      else if (imattype.eq.22) then
       call EOS_tillotson(rho,internal_energy,pressure)
      else if (imattype.eq.23) then
       call EOS_peng_robinson(rho,internal_energy,pressure)
      else
       print *,"imattype invalid EOS_material_CORE"
       stop
      endif

      return
      end subroutine EOS_material_CORE

      subroutine dVdT_material_CORE(dVdT,massfrac_var, &
        pressure,temperature, &
        imattype,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: pressure,temperature
      REAL_T, intent(in) :: massfrac_var(num_species_var+1)
      REAL_T, intent(out) :: dVdT


      if (pressure.gt.zero) then
       ! do nothing
      else
       print *,"pressure invalid"
       stop
      endif
      if (temperature.gt.zero) then
       ! do nothing
      else
       print *,"temperature invalid"
       stop
      endif
      if ((im.ge.1).and.(im.le.num_materials)) then
       ! do nothing
      else
       print *,"im invalid 16827:",im
       stop
      endif
      if (fort_stiffGAMMA(im).ge.one) then
       ! do nothing
      else
       print *,"fort_stiffGAMMA(im) invalid"
       stop
      endif
      if (fort_stiffCV(im).gt.zero) then
       ! do nothing
      else
       print *,"fort_stiffCV(im) invalid"
       stop
      endif

      if ((imattype.ge.1).and. &
          (imattype.le.23)) then
       dVdT=(fort_stiffGAMMA(im)-one) * fort_stiffCV(im)/pressure
      else if (imattype.eq.0) then
       dVdT=zero
      else if (imattype.eq.999) then
       dVdT=zero
      else
       print *,"imattype invalid in dVdT_material_CORE"
       stop
      endif

      return
      end subroutine dVdT_material_CORE


        ! in general for gas: e=cv T
        !                     p=(gamma-1)rho e=(gamma-1)rho cv T
        !                      =(cp-cv) rho T=rho R T
        ! total energy per unit mass? = (1/2)u dot u  + e
        ! returns c^2(e*scale)/scale
        ! sound squared=c^2(density=rho,internal_energy)
      subroutine SOUNDSQR_material_CORE(rho,massfrac_var, &
        internal_energy,soundsqr, &
        imattype,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho,internal_energy
      REAL_T, intent(in) :: massfrac_var(num_species_var+1)
      REAL_T, intent(out) :: soundsqr


      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"e invalid"
       stop
      endif

      if (imattype.eq.1) then
       call SOUNDSQR_tait(rho,internal_energy,soundsqr)
      else if (imattype.eq.2) then
       call SOUNDSQR_jwlADIABAT(rho,internal_energy,soundsqr)
      else if (imattype.eq.3) then
       call SOUNDSQR_NAjwl(rho,internal_energy,soundsqr)
      else if (imattype.eq.4) then
       call SOUNDSQR_SF6(rho,internal_energy,soundsqr)
      else if (imattype.eq.5) then
       call SOUNDSQR_air(rho,internal_energy,soundsqr)
      else if (imattype.eq.14) then
       call SOUNDSQR_air_rho2(rho,internal_energy,soundsqr)
      else if (imattype.eq.6) then
       call SOUNDSQR_Marquina(rho,internal_energy,soundsqr)
      else if (imattype.eq.7) then
       call SOUNDSQR_tait_rho(rho,internal_energy,soundsqr)
      else if (imattype.eq.13) then
       call SOUNDSQR_tait_rhohydro(rho,internal_energy,soundsqr)
      else if (imattype.eq.8) then
       call SOUNDSQR_airADIABAT(rho,internal_energy,soundsqr)
      else if (imattype.eq.9) then
       call SOUNDSQR_tait_rho3(rho,internal_energy,soundsqr)
      else if (imattype.eq.10) then
       call SOUNDSQR_tait_rho2(rho,internal_energy,soundsqr)
      else if (imattype.eq.11) then
       call SOUNDSQR_koren_rho1(rho,internal_energy,soundsqr)
      else if (imattype.eq.12) then
       call SOUNDSQR_koren_rho2(rho,internal_energy,soundsqr)
      else if (imattype.eq.15) then
       call SOUNDSQR_dodecane(rho,internal_energy,soundsqr)
      else if (imattype.eq.16) then
       call SOUNDSQR_SF6ADIABAT(rho,internal_energy,soundsqr)
      else if (imattype.eq.17) then
       call SOUNDSQR_stiffened(rho,internal_energy,soundsqr,im)
      else if (imattype.eq.18) then
       call SOUNDSQR_simple_air(rho,internal_energy,soundsqr)
      else if (imattype.eq.19) then
       call SOUNDSQR_vacuum(rho,internal_energy,soundsqr)
      else if (imattype.eq.20) then
       call SOUNDSQR_tait_vacuum(rho,internal_energy,soundsqr)
      else if (imattype.eq.21) then
       call SOUNDSQR_cav_nozzle(rho,internal_energy,soundsqr)
      else if (imattype.eq.22) then
       call SOUNDSQR_tillotson(rho,internal_energy,soundsqr)
      else if (imattype.eq.23) then
       call SOUNDSQR_peng_robinson(rho,internal_energy,soundsqr)
      else
       print *,"imattype invalid SOUNDSQR_material_CORE"
       stop
      endif

      return
      end subroutine SOUNDSQR_material_CORE


        ! returns e/scale
        ! internal energy = e(temperature,density=rho)
      subroutine INTERNAL_material_CORE(rho,massfrac_var, &
        temperature,internal_energy, &
        imattype,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho,temperature
      REAL_T, intent(in) :: massfrac_var(num_species_var+1)
      REAL_T, intent(out) :: internal_energy
      REAL_T local_internal_energy

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (temperature.le.zero) then
       print *,"T invalid"
       stop
      endif

      if (imattype.eq.999) then
       call INTERNAL_default(rho,temperature,local_internal_energy,imattype,im)
      else if (imattype.eq.0) then
       call INTERNAL_default(rho,temperature,local_internal_energy,imattype,im)
      else if (imattype.eq.1) then
       call INTERNAL_tait(rho,temperature,local_internal_energy)
      else if (imattype.eq.2) then
       call INTERNAL_jwl(rho,temperature,local_internal_energy)
      else if (imattype.eq.3) then
       call INTERNAL_jwl(rho,temperature,local_internal_energy)
      else if (imattype.eq.4) then
       call INTERNAL_SF6(rho,temperature,local_internal_energy)
      else if (imattype.eq.5) then
       call INTERNAL_air(rho,temperature,local_internal_energy)
      else if (imattype.eq.14) then
       call INTERNAL_air_rho2(rho,temperature,local_internal_energy)
      else if (imattype.eq.6) then
       call INTERNAL_Marquina(rho,temperature,local_internal_energy)
      else if (imattype.eq.7) then
       call INTERNAL_tait_rho(rho,temperature,local_internal_energy)
      else if (imattype.eq.13) then
       call INTERNAL_tait_rhohydro(rho,temperature,local_internal_energy)
      else if (imattype.eq.8) then
       call INTERNAL_airADIABAT(rho,temperature,local_internal_energy)
      else if (imattype.eq.9) then
       call INTERNAL_tait_rho3(rho,temperature,local_internal_energy)
      else if (imattype.eq.10) then
       call INTERNAL_tait_rho2(rho,temperature,local_internal_energy)
      else if (imattype.eq.11) then
       call INTERNAL_koren_rho1(rho,temperature,local_internal_energy)
      else if (imattype.eq.12) then
       call INTERNAL_koren_rho2(rho,temperature,local_internal_energy)
      else if (imattype.eq.15) then
       call INTERNAL_dodecane(rho,temperature,local_internal_energy)
      else if (imattype.eq.16) then
       call INTERNAL_SF6ADIABAT(rho,temperature,local_internal_energy)
      else if (imattype.eq.17) then
       call INTERNAL_stiffened(rho,temperature,local_internal_energy,im)
      else if (imattype.eq.18) then
       call INTERNAL_simple_air(rho,temperature,local_internal_energy)
      else if (imattype.eq.19) then
       call INTERNAL_vacuum(rho,temperature,local_internal_energy)
      else if (imattype.eq.20) then
       call INTERNAL_tait_vacuum(rho,temperature,local_internal_energy)
      else if (imattype.eq.21) then
       call INTERNAL_cav_nozzle(rho,temperature,local_internal_energy)
      else if (imattype.eq.22) then
       call INTERNAL_tillotson(rho,temperature,local_internal_energy)
      else if (imattype.eq.23) then
       call INTERNAL_peng_robinson(rho,temperature,local_internal_energy)
      else
       print *,"imattype invalid INTERNAL_material_CORE"
       stop
      endif

      internal_energy=local_internal_energy

      return
      end subroutine INTERNAL_material_CORE


      subroutine TEMPERATURE_material_CORE(rho,massfrac_var, &
        temperature,internal_energy, &
        imattype,im)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho,internal_energy
      REAL_T, intent(in) :: massfrac_var(num_species_var+1)
      REAL_T, intent(out) :: temperature

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid71"
       stop
      endif
      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy invalid in temperature material_CORE"
       print *,"rho,energy,imat ",rho,internal_energy,imattype 
       stop
      endif

      if (imattype.eq.999) then 
       call TEMPERATURE_default(rho,temperature,internal_energy, &
         imattype,im)
      else if (imattype.eq.0) then
       call TEMPERATURE_default(rho,temperature,internal_energy, &
         imattype,im)
      else if (imattype.eq.1) then
       call TEMPERATURE_tait(rho,temperature,internal_energy)
      else if (imattype.eq.2) then
       call TEMPERATURE_jwl(rho,temperature,internal_energy)
      else if (imattype.eq.3) then
       call TEMPERATURE_jwl(rho,temperature,internal_energy)
      else if (imattype.eq.4) then
       call TEMPERATURE_SF6(rho,temperature,internal_energy)
      else if (imattype.eq.5) then
       call TEMPERATURE_air(rho,temperature,internal_energy)
      else if (imattype.eq.14) then
       call TEMPERATURE_air_rho2(rho,temperature,internal_energy)
      else if (imattype.eq.6) then
       call TEMPERATURE_Marquina(rho,temperature,internal_energy)
      else if (imattype.eq.7) then
       call TEMPERATURE_tait_rho(rho,temperature,internal_energy)
      else if (imattype.eq.13) then
       call TEMPERATURE_tait_rhohydro(rho,temperature,internal_energy)
      else if (imattype.eq.8) then
       call TEMPERATURE_airADIABAT(rho,temperature,internal_energy)
      else if (imattype.eq.9) then
       call TEMPERATURE_tait_rho3(rho,temperature,internal_energy)
      else if (imattype.eq.10) then
       call TEMPERATURE_tait_rho2(rho,temperature,internal_energy)
      else if (imattype.eq.11) then
       call TEMPERATURE_koren_rho1(rho,temperature,internal_energy)
      else if (imattype.eq.12) then
       call TEMPERATURE_koren_rho2(rho,temperature,internal_energy)
      else if (imattype.eq.15) then
       call TEMPERATURE_dodecane(rho,temperature,internal_energy)
      else if (imattype.eq.16) then
       call TEMPERATURE_SF6ADIABAT(rho,temperature,internal_energy)
      else if (imattype.eq.17) then
       call TEMPERATURE_stiffened(rho,temperature,internal_energy,im)
      else if (imattype.eq.18) then
       call TEMPERATURE_simple_air(rho,temperature,internal_energy)
      else if (imattype.eq.19) then
       call TEMPERATURE_vacuum(rho,temperature,internal_energy)
      else if (imattype.eq.20) then
       call TEMPERATURE_tait_vacuum(rho,temperature,internal_energy)
      else if (imattype.eq.21) then
       call TEMPERATURE_cav_nozzle(rho,temperature,internal_energy)
      else if (imattype.eq.22) then
       call TEMPERATURE_tillotson(rho,temperature,internal_energy)
      else if (imattype.eq.23) then
       call TEMPERATURE_peng_robinson(rho,temperature,internal_energy)
      else
       print *,"imattype invalid TEMPERATURE_material_CORE"
       print *,"imattype= ",imattype
       stop
      endif

      return
      end subroutine TEMPERATURE_material_CORE


      subroutine general_hydrostatic_pressure(pres)
      use probcommon_module
      IMPLICIT NONE

      REAL_T pres


      pres=1.0D+6

      return
      end subroutine general_hydrostatic_pressure


      subroutine rampvel(time,vel)
      use probcommon_module
      IMPLICIT NONE

      REAL_T, intent(in) :: time
      REAL_T, intent(out) :: vel
      REAL_T :: tcutoff

      if ((SDIM.eq.2).and. &
          ((probtype.eq.63).or. &
           (probtype.eq.64))) then
       tcutoff=1.0
       if (time.gt.tcutoff) then
        vel=adv_vel
       else
        vel=time*adv_vel/tcutoff
       endif
      else if ((SDIM.eq.3).and.(probtype.eq.9)) then
       tcutoff=zero
       if ((time.gt.tcutoff).or.(tcutoff.eq.zero)) then
        vel=adv_vel
       else if (tcutoff.gt.zero) then
        vel=time*adv_vel/tcutoff
       else
        print *,"tcutoff invalid"
        stop
       endif
      else
       vel=adv_vel
      endif

      return
      end subroutine rampvel

      subroutine default_rampvel(time,xx_vel,yy_vel,zz_vel)
      use probcommon_module
      IMPLICIT NONE

      REAL_T, intent(in) :: time
      REAL_T, intent(out) :: xx_vel,yy_vel,zz_vel

      xx_vel=zero 
      yy_vel=zero 
      zz_vel=zero 

      if ((probtype.eq.26).and.(axis_dir.eq.11)) then
       ! do nothing 
      else if (adv_dir.eq.1) then  ! init velocity
        if (SDIM.eq.2) then
         call rampvel(time,xx_vel)
        else if (SDIM.eq.3) then
         xx_vel = adv_vel
          ! initboat
         if (probtype.eq.9) then
          call rampvel(time,xx_vel)
         endif
        else
         print *,"dimension bust"
         stop
        endif
        yy_vel=zero
        zz_vel=zero
      else if (adv_dir .eq. 2) then
        xx_vel = zero
        yy_vel = adv_vel
        zz_vel=zero
      else if ((adv_dir.eq.3).and.(SDIM.eq.3)) then
        xx_vel = zero
        yy_vel = zero
        zz_vel = adv_vel
      else if ((adv_dir.eq.3).and.(SDIM.eq.2)) then
        xx_vel = adv_vel
        yy_vel = adv_vel
        zz_vel = zero
      else if ((adv_dir.eq.4).and.(SDIM.eq.3)) then
        xx_vel = adv_vel
        yy_vel = adv_vel
        zz_vel = adv_vel
      else if (probtype.eq.40) then
       ! do nothing
      else if (probtype.eq.530) then
       ! do nothing
      else
        print *,"adv_dir invalid: ",adv_dir
        stop
      endif

      return
      end subroutine default_rampvel

      ! sigma_{i,j}cos(theta_{i,k})=sigma_{j,k}-sigma_{i,k}
      ! theta_{ik}=0 => material i wets material k.
      ! im is material "i"  ("fluid" material)
      ! im_opp is material "j"
      ! im_3 is material "k"
      ! iten_13 corresponds to "ik"
      ! iten_23 corresponds to "jk"
      ! if nmat=4, 12 13 14 23 24 34

      subroutine get_CL_iten(im,im_opp,im_3,iten_13,iten_23, &
       user_tension,nten,cos_angle,sin_angle)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T im,im_opp,im_3,iten_13,iten_23,nten,nten_test
      INTEGER_T iten
      INTEGER_T nmat
      REAL_T user_tension(nten)
      REAL_T cos_angle,sin_angle


      nmat=num_materials
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten: get_CL_iten nten nten_test ",nten,nten_test
       stop
      endif

      if ((im.lt.1).or.(im.gt.num_materials).or. &
          (im_opp.lt.1).or.(im_opp.gt.num_materials).or. &
          (im_3.lt.1).or.(im_3.gt.num_materials).or. &
          (im.eq.im_opp).or.(im.eq.im_3).or. &
          (im_opp.eq.im_3)) then
       print *,"im mismatch"
       print *,"im=",im
       print *,"im_opp=",im_opp
       print *,"im_3=",im_3
       stop
      endif

      if (nmat.le.2) then
       print *,"nmat too small for CL treatment"
       stop
       ! 12 13 23
      else if (nmat.eq.3) then
       if ((im.eq.1).and.(im_opp.eq.2).and.(im_3.eq.3)) then
        iten_13=2
        iten_23=3
       else if ((im.eq.2).and.(im_opp.eq.3).and.(im_3.eq.1)) then
        iten_13=1
        iten_23=2
       else if ((im.eq.1).and.(im_opp.eq.3).and.(im_3.eq.2)) then
        iten_13=1
        iten_23=3
       else
        print *,"combination of im,im_opp,im_3 invalid nmat=",nmat
        print *,"im=",im
        print *,"im_opp=",im_opp
        print *,"im_3=",im_3
        stop
       endif
       ! 12 13 14 23 24 34
      else if (nmat.eq.4) then
       if ((im.eq.1).and.(im_opp.eq.2).and.(im_3.eq.3)) then
        iten_13=2
        iten_23=4
       else if ((im.eq.2).and.(im_opp.eq.3).and.(im_3.eq.1)) then
        iten_13=1
        iten_23=2
       else if ((im.eq.1).and.(im_opp.eq.3).and.(im_3.eq.2)) then
        iten_13=1
        iten_23=4
       else if ((im.eq.1).and.(im_opp.eq.2).and.(im_3.eq.4)) then
        iten_13=3
        iten_23=5
       else if ((im.eq.2).and.(im_opp.eq.4).and.(im_3.eq.1)) then
        iten_13=1
        iten_23=3
       else if ((im.eq.1).and.(im_opp.eq.4).and.(im_3.eq.2)) then
        iten_13=1
        iten_23=5
       else if ((im.eq.2).and.(im_opp.eq.3).and.(im_3.eq.4)) then
        iten_13=5
        iten_23=6
       else if ((im.eq.2).and.(im_opp.eq.4).and.(im_3.eq.3)) then
        iten_13=4
        iten_23=6
       else if ((im.eq.3).and.(im_opp.eq.4).and.(im_3.eq.2)) then
        iten_13=4
        iten_23=5
       else if ((im.eq.3).and.(im_opp.eq.4).and.(im_3.eq.1)) then
        iten_13=2
        iten_23=3
        ! 12 13 14 23 24 34
       else if ((im.eq.1).and.(im_opp.eq.3).and.(im_3.eq.4)) then
        iten_13=3
        iten_23=6
       else if ((im.eq.1).and.(im_opp.eq.4).and.(im_3.eq.3)) then
        iten_13=2
        iten_23=6
       else
        print *,"combination of im,im_opp,im_3 invalid nmat=",nmat
        print *,"im=",im
        print *,"im_opp=",im_opp
        print *,"im_3=",im_3
        stop
       endif
      else
       print *,"combination not supported"
       stop
      endif

      call get_iten(im,im_opp,iten,nmat)

      if (user_tension(iten).eq.zero) then  ! default extrapolation
       cos_angle=zero
      else if (user_tension(iten).gt.zero) then
       cos_angle=(user_tension(iten_23)-user_tension(iten_13))/ &
                 user_tension(iten)
      else
       print *,"user_tension(iten) invalid"
       stop
      endif

      if (cos_angle.gt.one) then 
       cos_angle=one
      else if (cos_angle.lt.-one) then
       cos_angle=-one
      endif
      sin_angle=sqrt(one-cos_angle**2)

      return
      end subroutine get_CL_iten

  ! in 2D, "y" is ignored
  ! distance>0 in the drop
  ! if maxtall<radnew-vert then the distance
  ! to the ice part of the droplet is returned in dist
  ! and the distance to the remaining liquid part *above* the ice
  ! is returned in dist_truncate.
  ! this routine is called if probtype=55, axis_dir=0,1, or 5
subroutine drop_slope_dist(x,y,z,time,nmat, &
   maxtall,dist,dist_truncate)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x,y,z,time
REAL_T, intent(out) :: dist,dist_truncate
REAL_T, intent(in) :: maxtall
INTEGER_T im,im_opp,im_3,iten_13,iten_23,imloop
INTEGER_T iten
REAL_T cos_angle,sin_angle
REAL_T term1,Vtarget,radnew,vert,test_angle
REAL_T xprime,yprime,zprime,rprime,rtop,rbot
REAL_T xcheck,ycheck,zcheck
INTEGER_T nten
REAL_T xvec(SDIM)
REAL_T marangoni_temp(nmat)
INTEGER_T im_solid_substrate
REAL_T, allocatable, dimension(:) :: user_tension

if (probtype.eq.55) then

 im_solid_substrate=im_solid_primary()

 xvec(1)=x
 xvec(2)=y
 if (SDIM.eq.3) then
  xvec(SDIM)=z
 endif

 if (nmat.ne.num_materials) then
  print *,"nmat invalid"
  stop
 endif
 nten=( (nmat-1)*(nmat-1)+nmat-1 )/2
 allocate(user_tension(nten))

 if (SDIM.eq.2) then
  if (abs(y-z).gt.VOFTOL) then
   print *,"y=z in 2d expected: drop_slope_dist"
   print *,"x,y,z= ",x,y,z
   stop
  endif
 endif
 if (maxtall.le.zero) then
  print *,"maxtall invalid"
  stop
 endif

  ! in: drop_slope_dist 
 if ((axis_dir.eq.0).or. &
     (axis_dir.eq.5)) then

  xcheck=xblob-xblob2
  ycheck=yblob-yblob2
  zcheck=zero
  if (SDIM.eq.3) then
   zcheck=zblob-zblob2
  endif

  if ((num_materials.ge.3).and. &
      (im_solid_substrate.ge.3).and. &
      (abs(xcheck).lt.1.0D-7).and. &
      (abs(ycheck).lt.1.0D-7).and. &
      (abs(zcheck).lt.1.0D-7)) then
   im=1
   im_opp=2
   im_3=im_solid_substrate
   call get_iten(im,im_opp,iten,num_materials)
   do imloop=1,nmat
    marangoni_temp(imloop)=293.0
   enddo
   call get_user_tension(xvec,time, &
     fort_tension,user_tension, &
     marangoni_temp, &
     nmat,nten,1)
     ! find angle between materials "im" and "im_3"
   call get_CL_iten(im,im_opp,im_3,iten_13,iten_23, &
    user_tension,nten,cos_angle,sin_angle)

    ! angles other than 0 or pi are supported:
    ! 0 < angle < pi
   if (abs(cos_angle).lt.one-1.0D-2) then 

    if (((SDIM.eq.3).and.(levelrz.eq.0)).or. &
        ((SDIM.eq.2).and.(levelrz.eq.1))) then
     term1=two/three-cos_angle+(cos_angle**3)/three
     if (term1.le.zero) then
      print *,"term1 invalid"
      stop
     endif
         
     Vtarget=half*(four/three)*Pi*(radblob**3)
     radnew=(Vtarget/(Pi*term1))**(one/three)
     vert=-radnew*cos_angle
    else if ((SDIM.eq.2).and.(levelrz.eq.0)) then
     test_angle=acos(abs(cos_angle))  ! 0<test_angle<=pi/2
     if (cos_angle.ge.zero) then
      term1=test_angle-half*sin(two*test_angle)
     else
      term1=Pi-test_angle+half*sin(two*test_angle)
     endif
     if (term1.le.zero) then
      print *,"term1 invalid"
      stop
     endif
     Vtarget=half*Pi*(radblob**2)
     radnew=sqrt(Vtarget/term1)
     vert=-radnew*cos_angle
    else
     print *,"dimension bust"
     stop
    endif
      ! rotate clockwise
      ! and shift "center" of inclined plane to origin
      ! need to modify if rotate about y-z plane.
    if (SDIM.eq.2) then
     xprime=(x-xblob2)*cos(radblob2)+(z-yblob2)*sin(radblob2)
     zprime=-(x-xblob2)*sin(radblob2)+(z-yblob2)*cos(radblob2)
     rprime=abs(xprime)
    else if (SDIM.eq.3) then
     xprime=(x-xblob2)*cos(radblob2)+(z-zblob2)*sin(radblob2)
     yprime=y-yblob2
     zprime=-(x-xblob2)*sin(radblob2)+(z-zblob2)*cos(radblob2)
     rprime=sqrt(xprime**2+yprime**2)
    else
     print *,"dimension bust"
     stop
    endif

     ! dist>0 in the liquid drop
    dist=radnew-sqrt(rprime**2+(zprime-vert)**2)
    dist_truncate=dist

     ! find distance to ice part of this droplet; also
     ! find distance to remaining liquid part above the ice.
    if (maxtall-vert.lt.radnew) then
     rtop=sqrt(radnew**2-(maxtall-vert)**2)
     rbot=sqrt(radnew**2-vert**2)

      ! outside drop, and above the ice.
     if ((dist.le.zero).and.(zprime.ge.maxtall)) then
      dist_truncate=dist

      ! inside the original drop.
     else if ((zprime.le.maxtall).and.(rprime.le.rtop)) then
      dist_truncate=zprime-maxtall

      ! outside the original drop, off to side of ice.
     else if ((zprime.le.maxtall).and.(rprime.ge.rtop)) then
      dist_truncate=-sqrt((rprime-rtop)**2+(zprime-maxtall)**2)
     else if (dist.ge.zero) then
      if (dist.lt.zprime-maxtall) then
       dist_truncate=dist
      else
       dist_truncate=zprime-maxtall
      endif 
     else
      print *,"dist invalid drop_slope_dist"
      stop
     endif

     if ((dist.lt.zero).and.(zprime.gt.vert+radnew)) then
      dist=maxtall-zprime
     else if ((dist.ge.zero).and.(zprime.ge.maxtall)) then
      dist=maxtall-zprime
     else if ((dist.ge.zero).and.(zprime.le.maxtall).and. &
              (zprime.ge.half*maxtall)) then
      if (dist.lt.maxtall-zprime) then
       ! do nothing
      else
       dist=maxtall-zprime
      endif
     else if ((dist.ge.zero).and.(zprime.le.half*maxtall).and. &
              (zprime.ge.zero)) then
      if (dist.lt.zprime) then
       ! do nothing
      else
       dist=zprime
      endif
     else if ((dist.ge.zero).and.(zprime.le.zero)) then
      dist=zprime
     else if ((dist.lt.zero).and.(zprime.lt.vert-radnew)) then
      dist=zprime
     else if ((dist.lt.zero).and.(zprime.ge.maxtall)) then
      dist=-sqrt((rprime-rtop)**2+(zprime-maxtall)**2)
     else if ((dist.lt.zero).and.(zprime.ge.zero)) then
      ! do nothing
     else if ((dist.lt.zero).and.(zprime.le.zero)) then
      dist=-sqrt((rprime-rbot)**2+zprime**2)
     else
      print *,"dist or zprime invalid"
      stop
     endif 
    endif !  maxtall-vert<radnew
     
   else
    print *,"contact angle too close to 0 or pi for drop on slope"
    print *,"probtype=",probtype
    print *,"radblob=",radblob
    print *,"radblob2=",radblob2
    print *,"radblob4=",radblob4
    print *,"radblob5=",radblob5
    print *,"radblob6=",radblob6
    print *,"radblob7=",radblob7
    stop
   endif

  else
   print *,"parameter conflict for probtype=55"
   stop
  endif

 else
  print *,"axis_dir incorrect   probtype,axis_dir=",probtype,axis_dir
  stop
 endif

 deallocate(user_tension)

else
 print *,"probtype invalid in drop_slope_dist"
 stop
endif

return
end subroutine drop_slope_dist

subroutine volfrac_from_massfrac(X,Y,WA,WV)
IMPLICIT NONE

REAL_T, intent(in) :: Y,WA,WV
REAL_T, intent(out) :: X

if ((Y.ge.zero).and.(Y.le.one)) then
 if ((WA.gt.zero).and.(WV.gt.zero)) then
  X=WA*Y/(WV*(one-Y)+WA*Y)
 else
  print *,"WA or WV invalid"
  stop
 endif
else
 print *,"Y invalid in volfrac_from_massfrac"
 stop
endif

return
end subroutine volfrac_from_massfrac

! Kassemi, Kartuzova, Hylton
! Cryogenics 89(2018) 1-15, equation (6)
subroutine MDOT_Kassemi(sigma,MolarMassFluid,R,Pgamma,Pvapor_probe, &
  Tgamma,Tvapor_probe,MDOT)
IMPLICIT NONE

REAL_T, intent(in) :: sigma
REAL_T, intent(in) :: MolarMassFluid
REAL_T, intent(in) :: R
REAL_T, intent(in) :: Pgamma
REAL_T, intent(in) :: Pvapor_probe
REAL_T, intent(in) :: Tgamma
REAL_T, intent(in) :: Tvapor_probe
REAL_T, intent(out) :: MDOT

if ((sigma.ge.zero).and.(sigma.lt.two).and. &
    (MolarMassFluid.gt.zero).and. &
    (R.gt.zero).and. &
    (Pgamma.gt.zero).and. &
    (Pvapor_probe.gt.zero).and. &
    (Tgamma.gt.zero).and. &
    (Tvapor_probe.gt.zero)) then
 MDOT=(2.0d0*sigma/(2.0d0-sigma))* &
   sqrt(MolarMassFluid/(2.0d0*Pi*R))* &
   (1.0d0/sqrt(Tgamma))* &
   (Pgamma-Pvapor_probe)

else
 print *,"parameter problems in MDOT_Kassemi"
 stop
endif

return
end subroutine MDOT_Kassemi
 
! TODO: 1. inverse CC equation
!       2. MDOTY equation (either PD or Kassemi)
!       3. X from PSAT,Pgamma,Tgamma,TSAT,WV,WA,Tprobe,Yprobe?
!       4. For now, X=Pgamma/PSAT and focus on validating the
!          results found in Kassemi et al.
! When comparing with Dodd and Ferrante and others:
! TSAT corresponds to TBOIL
! PSAT corresponds to PBOIL
! Tgamma corresponds to Tsat
! Pgamma corresponds to Psat
subroutine Pgamma_Clausius_Clapyron(Pgamma,PSAT,Tgamma,TSAT,L,R,WV)
IMPLICIT NONE

REAL_T, intent(in) :: PSAT,Tgamma,TSAT,L,R,WV
REAL_T, intent(out) :: Pgamma

if ((Tgamma.gt.zero).and.(TSAT.gt.zero)) then
 if (PSAT.gt.zero) then
  if (R.gt.zero) then
   if (L.ne.zero) then
    if (WV.gt.zero) then
     Pgamma=PSAT*exp(-(L*WV/R)*(one/Tgamma-one/TSAT))
    else
     print *,"WV invalid in Pgamma_Clausius_Clapyron"
     stop
    endif
   else
    print *,"L invalid in Pgamma_Clausius_Clapyron"
    stop
   endif
  else
   print *,"R invalid in Pgamma_Clausius_Clapyron"
   stop
  endif
 else
  print *,"PSAT invalid in Pgamma_Clausius_Clapyron"
  stop
 endif
else
 print *,"Tgamma or TSAT invalid in Pgamma_Clausius_Clapyron"
 print *,"Tgamma= ",Tgamma
 print *,"TSAT= ",TSAT
 stop
endif

return
end subroutine Pgamma_Clausius_Clapyron

!TODO:
!XV=(PSAT_REF/PMIX)e^(-(L WV/R)(1/T_GAMMA-1/T_SAT_REF)
!PMIX=(gamma_MIX(YPROBE)-1)rho_MIX CV_MIX(YPROBE) TPROBE
! new parameters: YPROBE,TPROBE,RHO_PROBE in the gas.
subroutine X_from_Tgamma(X,Tgamma,TSAT, &
  L,R,WV)
IMPLICIT NONE

REAL_T, intent(in) :: Tgamma,TSAT,L,R,WV
REAL_T, intent(out) :: X

if ((Tgamma.gt.zero).and.(TSAT.gt.zero)) then
 if (R.gt.zero) then
  if (L.ne.zero) then
   if (WV.gt.zero) then
    if (L.gt.zero) then  ! evaporation
     if (Tgamma.le.TSAT) then
      ! do nothing
     else
      print *,"Tgamma exceeds TSAT"
      print *,"Tgamma,TSAT,R,L,WV ",Tgamma,TSAT,R,L,WV
      stop
     endif
    else if (L.lt.zero) then ! condensation
     if (Tgamma.ge.TSAT) then
      ! do nothing
     else
      print *,"Tgamma is below TSAT"
      stop
     endif
    else
     print *,"L invalid"
     stop
    endif
    X=exp(-(L*WV/R)*(one/Tgamma-one/TSAT))
   else
    print *,"WV invalid in X_from_Tgamma"
    stop
   endif
  else
   print *,"L invalid in X_from_Tgamma"
   stop
  endif
 else
  print *,"R invalid in X_from_Tgamma"
  stop
 endif
else
 print *,"Tgamma or TSAT invalid in X_from_Tgamma"
 print *,"Tgamma= ",Tgamma
 print *,"TSAT= ",TSAT
 stop
endif

return
end subroutine X_from_Tgamma

subroutine XMIN_from_TSAT(X,TSAT, &
  L,R,WV)
IMPLICIT NONE

REAL_T, intent(in) :: TSAT,L,R,WV
REAL_T, intent(out) :: X

if (TSAT.gt.zero) then
 if (R.gt.zero) then
  if (L.ne.zero) then
   if (WV.gt.zero) then
    if (L.gt.zero) then
     X=zero
    else if (L.lt.zero) then
     X=exp((L*WV/R)/TSAT)
    else
     print *,"L invalid in XMIN_from_TSAT1"
     stop
    endif
   else
    print *,"WV invalid in XMIN_from_TSAT"
    stop
   endif
  else
   print *,"L invalid in XMIN_from_TSAT2"
   stop
  endif
 else
  print *,"R invalid in XMIN_from_TSAT"
  stop
 endif
else
 print *,"TSAT invalid in XMIN_from_TSAT"
 print *,"TSAT= ",TSAT
 stop
endif

return
end subroutine XMIN_from_TSAT

! Dodd and Ferrante:
! (30) psat/pboil=exp((-L WV/R)(1/Tsat - 1/Tboil))
! Dodd and Ferrante after nondimensionalization by pboil:
! (33) Ysat=psat WV/(psat WV + (1-psat)WA)
! in otherwords,
! Ysat=X WV/(X WV + (1-X)WA) where X=psat/pboil
!TODO:
!XV=(PSAT_REF/PMIX)e^(-(L WV/R)(1/T_GAMMA-1/T_SAT_REF)
!PMIX=(gamma_MIX(YPROBE)-1)rho_MIX CV_MIX(YPROBE) TPROBE
! new parameters: YPROBE,TPROBE,RHO_PROBE in the gas.
subroutine Tgamma_from_TSAT_and_X(Tgamma,TSAT, &
  X,L,R,WV,Tgamma_min,Tgamma_max)
IMPLICIT NONE

REAL_T, intent(in) :: TSAT,X,L,R,WV,Tgamma_min,Tgamma_max
REAL_T, intent(out) :: Tgamma

if ((X.ge.zero).and.(X.le.one)) then
 if (TSAT.gt.zero) then

   ! X=exp(-(L WV/R)(1/Tgamma - 1/Tsat) )
  if (X.eq.zero) then
   if (L.gt.zero) then
    Tgamma=Tgamma_min
   else if (L.lt.zero) then
    Tgamma=Tgamma_max
   else
    print *,"L invalid in Tgamma_from_TSAT_and_X"
    stop
   endif
  else if (X.eq.one) then
   Tgamma=TSAT
  else if ((X.gt.zero).and.(X.lt.one)) then
   Tgamma=-log(X)*R/(L*WV)
   Tgamma=Tgamma+one/TSAT
   if (Tgamma.gt.zero) then
    Tgamma=one/Tgamma
   else if (Tgamma.le.zero) then
    Tgamma=Tgamma_min
   else
    print *,"Tgamma invalid in Tgamma_from_TSAT_and_X"
    stop
   endif
   if (L.gt.zero) then
    if (Tgamma.le.TSAT) then
     ! do nothing
    else
     print *,"expecting Tgamma<=TSAT"
     stop
    endif
   else if (L.lt.zero) then
    if (Tgamma.ge.TSAT) then
     ! do nothing
    else
     print *,"expecting Tgamma>=TSAT"
     stop
    endif
   else
    print *,"L invalid"
    stop
   endif

  else
   print *,"X invalid in Tgamma_from_TSAT_and_X"
   stop
  endif
  if (Tgamma.lt.Tgamma_min) then
   Tgamma=Tgamma_min
  endif
  if (Tgamma.gt.Tgamma_max) then
   Tgamma=Tgamma_max
  endif

 else
  print *,"TSAT invalid in Tgamma_from_TSAT_and_X"
  stop
 endif
else
 print *,"X invalid in Tgamma_from_TSAT_and_X"
 stop
endif
 
return
end subroutine Tgamma_from_TSAT_and_X


subroutine massfrac_from_volfrac(X,Y,WA,WV)
IMPLICIT NONE

REAL_T, intent(in) :: X,WA,WV
REAL_T, intent(out) :: Y

if ((X.ge.zero).and.(X.le.one)) then
 if ((WA.gt.zero).and.(WV.gt.zero)) then
  Y=WV*X/(WA*(one-X)+WV*X)
 else
  print *,"WA or WV invalid"
  stop
 endif
else
 print *,"X invalid in massfrac_from_volfrac"
 print *,"X,WA,WV ",X,WA,WV
 stop
endif

return
end subroutine massfrac_from_volfrac

  ! returns p(e*scale)/scale
  ! pressure=p(density=rho,internal_energy)
subroutine EOS_material(rho,massfrac_parm, &
  internal_energy_in,pressure, &
  imattype,im)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: imattype,im
REAL_T, intent(in) :: rho
REAL_T, intent(in) :: massfrac_parm(num_species_var+1)
REAL_T :: internal_energy
REAL_T, intent(out) :: pressure
REAL_T, intent(in) :: internal_energy_in
INTEGER_T :: ispec


internal_energy=internal_energy_in*global_pressure_scale

do ispec=1,num_species_var
 if (massfrac_parm(ispec).ge.zero) then
  ! do nothing
 else
  print *,"massfrac_parm invalid"
  stop
 endif
enddo ! ispec
if (rho.gt.zero) then
 ! do nothing
else
 print *,"rho invalid"
 stop
endif
if (internal_energy.gt.zero) then
 ! do nothing
else
 print *,"e invalid"
 stop
endif

if (is_in_probtype_list().eq.1) then
 call SUB_EOS(rho,massfrac_parm, &
   internal_energy,pressure, &
   imattype,im,num_species_var)
else 
 call EOS_material_CORE(rho,massfrac_parm, &
         internal_energy,pressure,imattype,im)
endif

pressure=pressure/global_pressure_scale

return
end subroutine EOS_material


subroutine dVdT_material(dVdT,massfrac_parm, &
  pressure_in,temperature, &
  imattype,im)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: imattype,im
REAL_T, intent(in) :: pressure_in,temperature
REAL_T, intent(in) :: massfrac_parm(num_species_var+1)
REAL_T, intent(out) :: dVdT
REAL_T :: pressure
INTEGER_T :: ispec

if (pressure_in.gt.zero) then
 ! do nothing
else
 print *,"pressure_in invalid"
 stop
endif
if (temperature.gt.zero) then
 ! do nothing
else
 print *,"temperature invalid"
 stop
endif
if ((im.ge.1).and.(im.le.num_materials)) then
 ! do nothing
else
 print *,"im invalid 17977:",im
 stop
endif
if (fort_stiffGAMMA(im).ge.one) then
 ! do nothing
else
 print *,"fort_stiffGAMMA(im) invalid"
 stop
endif
if (fort_stiffCV(im).gt.zero) then
 ! do nothing
else
 print *,"fort_stiffCV(im) invalid"
 stop
endif

do ispec=1,num_species_var
 if (massfrac_parm(ispec).ge.zero) then
  ! do nothing
 else
  print *,"massfrac_parm invalid"
  stop
 endif
enddo ! ispec

pressure=pressure_in*global_pressure_scale

if (is_in_probtype_list().eq.1) then
 call SUB_dVdT(dVdT,massfrac_parm, &
   pressure,temperature, &
   imattype,im,num_species_var)
else 
 call dVdT_material_CORE(dVdT,massfrac_parm, &
         pressure,temperature, &
         imattype,im)
endif

return
end subroutine dVdT_material


  ! returns De/DT / scale 
  ! T=(e/scale)*scale/cv
subroutine DeDT_material(rho,massfrac_parm, &
  temperature,DeDT,imattype,im)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: imattype,im
REAL_T, intent(in) :: rho,temperature
REAL_T, intent(in) :: massfrac_parm(num_species_var+1)
REAL_T, intent(out) :: DeDT
REAL_T :: DT,T2,e1,e2


if ((rho.gt.zero).and.(temperature.gt.zero)) then
 ! do nothing
else
 print *,"rho or temperature invalid in DeDT_material"
 print *,"rho=",rho
 print *,"temperature=",temperature
 print *,"im=",im
 print *,"imattype=",imattype
 stop
endif

if (imattype.eq.18) then ! simple_air
 DeDT=one 
else if ((imattype.eq.15).and.(probtype.eq.541).and. &
    (OLD_DODECANE.eq.0)) then
 call DeDT_dodecane(rho,temperature,DeDT)
else if ((imattype.eq.999).or. &
         (imattype.eq.0).or. &
         ((imattype.ge.1).and.(imattype.le.MAX_NUM_EOS))) then
 call INTERNAL_material(rho,massfrac_parm, &
   temperature,e1,imattype,im)
 DT=temperature*1.0D-6
 T2=temperature+DT
 call INTERNAL_material(rho,massfrac_parm, &
   T2,e2,imattype,im)
 DeDT=(e2-e1)/DT
 if (DeDT.le.zero) then
  print *,"e must be increasing function of T"
  stop
 endif
else
 print *,"imattype invalid DE_dt__material"
 stop
endif

return
end subroutine DeDT_material

  ! in general for gas: e=cv T
  !                     p=(gamma-1)rho e=(gamma-1)rho cv T
  !                      =(cp-cv) rho T=rho R T
  ! total energy per unit mass? = (1/2)u dot u  + e
  ! returns c^2(e*scale)/scale
  ! sound squared=c^2(density=rho,internal_energy)
subroutine SOUNDSQR_material(rho,massfrac_parm, &
  internal_energy_in,soundsqr, &
  imattype,im)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: imattype,im
REAL_T, intent(in) :: rho
REAL_T, intent(in) :: massfrac_parm(num_species_var+1)
REAL_T :: internal_energy
REAL_T, intent(out) :: soundsqr
REAL_T, intent(in) :: internal_energy_in
INTEGER_T :: ispec


internal_energy=internal_energy_in*global_pressure_scale

if (rho.gt.zero) then
 ! do nothing
else
 print *,"rho invalid"
 stop
endif
if (internal_energy.gt.zero) then
 ! do nothing
else
 print *,"e invalid"
 stop
endif
do ispec=1,num_species_var
 if (massfrac_parm(ispec).ge.zero) then
  ! do nothing
 else
  print *,"massfrac_parm invalid"
  stop
 endif
enddo

if (is_in_probtype_list().eq.1) then
 call SUB_SOUNDSQR(rho,massfrac_parm, &
   internal_energy,soundsqr, &
   imattype,im,num_species_var)
else 
 call SOUNDSQR_material_CORE(rho,massfrac_parm, &
   internal_energy,soundsqr, &
   imattype,im)
endif

soundsqr=soundsqr/global_pressure_scale

return
end subroutine SOUNDSQR_material


  ! returns e/scale
subroutine INTERNAL_material_cutoff(rho,internal_energy, &
  imattype,im)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: imattype,im
REAL_T, intent(in) :: rho
REAL_T, intent(out) :: internal_energy
REAL_T :: massfrac_parm(num_species_var+1)


if (rho.gt.zero) then
 ! do nothing
else
 print *,"rho invalid"
 stop
endif
if ((im.lt.1).or.(im.gt.num_materials)) then
 print *,"im invalid69"
 stop
endif
call init_massfrac_parm(rho,massfrac_parm,im)

  ! returns e/scale
call INTERNAL_material(rho,massfrac_parm, &
 fort_tempcutoff(im), &
 internal_energy,imattype,im)

return
end subroutine INTERNAL_material_cutoff

  ! returns e/scale
  ! internal energy = e(temperature,density=rho)
subroutine INTERNAL_material(rho,massfrac_parm, &
  temperature,internal_energy, &
  imattype,im)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: imattype,im
REAL_T, intent(in) :: rho,temperature
REAL_T, intent(in) :: massfrac_parm(num_species_var+1)
REAL_T, intent(out) :: internal_energy
REAL_T local_internal_energy  ! this is an output

if (rho.gt.zero) then
 ! do nothing
else
 print *,"rho invalid"
 stop
endif
if (temperature.gt.zero) then
 ! do nothing
else
 print *,"T invalid"
 stop
endif

if (is_in_probtype_list().eq.1) then
 call SUB_INTERNAL(rho,massfrac_parm, &
   temperature,local_internal_energy, &
   imattype,im,num_species_var)
else 
 call INTERNAL_material_CORE(rho,massfrac_parm, &
  temperature,local_internal_energy, &
  imattype,im)
endif

internal_energy=local_internal_energy/global_pressure_scale

return
end subroutine INTERNAL_material


 ! called from fort_override 
subroutine init_density_at_depth()
use probcommon_module
IMPLICIT NONE

REAL_T depth,pgrad,a,b,c,tol
REAL_T surface_den,depth_den
REAL_T surface_pressure,depth_pressure


density_at_depth=one

 ! water density where the charge is initially located.
if ((probtype.eq.42).and.(SDIM.eq.2)) then  ! bubble jetting
 density_at_depth=1.00039080D0
  ! cavitation
else if ((probtype.eq.46).and.(SDIM.eq.2)) then 
 if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
  density_at_depth=1.00008343D0
 else if (axis_dir.eq.10) then
  density_at_depth=1.0000423520369408D0 ! 10.22 meters
 else if (axis_dir.eq.20) then
  ! do nothing
 else
  print *,"axis_dir out of range"
  stop
 endif
else if (fort_material_type(1).eq.13) then

 if (abs(gravity).eq.zero) then
  print *,"gravity invalid"
  stop
 endif

 if (1.eq.0) then
  print *,"abs(gravity)= ",abs(gravity)
  print *,"denconst(1)= ",fort_denconst(1)
 endif

 surface_den=fort_denconst(1)  ! density at top of domain
 depth_den=surface_den

 if (SDIM.eq.2) then
  depth=probhiy-probloy
 else if (SDIM.eq.3) then
  depth=probhiz-probloz
 else
  print *,"dimension bust"
  stop
 endif
 
!       print *,"depth= ",depth

 if (depth.le.zero) then
  print *,"depth invalid"
  stop
 endif
 if (surface_den.le.zero) then
  print *,"surface den invalid"
  stop
 endif
 call EOS_tait_ADIABATIC_rhohydro(surface_den,surface_pressure)
 call EOS_tait_ADIABATIC_rhohydro(depth_den,depth_pressure)
 pgrad=abs(depth_pressure-surface_pressure)/(depth*surface_den)
 do while (pgrad.le.abs(gravity))
  depth_den=two*depth_den
  call EOS_tait_ADIABATIC_rhohydro(depth_den,depth_pressure)
  pgrad=abs(depth_pressure-surface_pressure)/(depth*surface_den)
  if (1.eq.0) then
   print *,"depth_den,pgrad ",depth_den,pgrad
  endif
 enddo ! while (pgrad.le.abs(gravity))
 a=surface_den
 b=depth_den
 c=half*(surface_den+depth_den)
 call EOS_tait_ADIABATIC_rhohydro(c,depth_pressure)
 pgrad=abs(depth_pressure-surface_pressure)/(depth*surface_den)
 tol=abs(gravity)*1.0D-3
 do while (abs(pgrad-abs(gravity)).gt.tol)
  if (1.eq.0) then
   print *,"a,b,c ",a,b,c
   print *,"pgrad ",pgrad
  endif
  if (pgrad.gt.abs(gravity)) then
   b=c
   c=half*(a+b)
  else if (pgrad.lt.abs(gravity)) then
   a=c
   c=half*(a+b)
  else
   print *,"pgrad bust"
   stop
  endif
  call EOS_tait_ADIABATIC_rhohydro(c,depth_pressure)
  pgrad=abs(depth_pressure-surface_pressure)/(depth*surface_den)
 enddo
 density_at_depth=c
endif  ! material_type(1)=13

return
end subroutine init_density_at_depth


subroutine get_left_velocityIOWA(t,u_left)
use probcommon_module
IMPLICIT NONE
INTEGER_T k1,k2
REAL_T t,u_left

if (probtype.ne.110) then
 print *,"probtype invalid get left velocity"
 stop
endif
if (inflow_count.le.0) then
 print *,"inflow_count invalid"
 stop
endif

do while ((inflow_time(last_inflow_index).gt.t).and. &
          (last_inflow_index.gt.1)) 
 last_inflow_index=last_inflow_index-1
enddo
do while ((inflow_time(last_inflow_index).lt.t).and. &
          (last_inflow_index.lt.inflow_count))
 last_inflow_index=last_inflow_index+1
enddo
if (t.le.inflow_time(1)) then
 u_left=inflow_velocity(1)
else if (t.ge.inflow_time(inflow_count)) then
 u_left=inflow_velocity(inflow_count)
else 

 if (t.le.inflow_time(last_inflow_index)) then
  k1=last_inflow_index-1
 else
  k1=last_inflow_index
 endif
 k2=k1+1
 u_left=inflow_velocity(k1)+ &
  ( (inflow_velocity(k2)-inflow_velocity(k1))/ &
    (inflow_time(k2)-inflow_time(k1)) )*(t-inflow_time(k1))
endif
u_left=u_left*100.0  ! convert to cm/s

return
end subroutine get_left_velocityIOWA


subroutine get_right_velocityIOWA(t,u_right)
use probcommon_module
IMPLICIT NONE
INTEGER_T k1,k2
REAL_T t,u_right

if (probtype.ne.110) then
 print *,"probtype invalid get_right_velocity"
 stop
endif
if (outflow_count.le.0) then
 print *,"outflow_count invalid"
 stop
endif

do while ((outflow_time(last_outflow_index).gt.t).and. &
          (last_outflow_index.gt.1)) 
 last_outflow_index=last_outflow_index-1
enddo
do while ((outflow_time(last_outflow_index).lt.t).and. &
          (last_outflow_index.lt.outflow_count))
 last_outflow_index=last_outflow_index+1
enddo
if (t.le.outflow_time(1)) then
 u_right=outflow_velocity(1)
else if (t.ge.outflow_time(outflow_count)) then
 u_right=outflow_velocity(outflow_count)
else 

 if (t.le.outflow_time(last_outflow_index)) then
  k1=last_outflow_index-1
 else
  k1=last_outflow_index
 endif
 k2=k1+1
 u_right=outflow_velocity(k1)+ &
  ( (outflow_velocity(k2)-outflow_velocity(k1))/ &
    (outflow_time(k2)-outflow_time(k1)) )*(t-outflow_time(k1))
endif
u_right=u_right*100.0  ! convert to cm/s

return
end subroutine get_right_velocityIOWA


subroutine get_left_elevationIOWA(t,elevation_left)
use probcommon_module
IMPLICIT NONE
INTEGER_T k1,k2
REAL_T t,elevation_left

if (probtype.ne.110) then
 print *,"probtype invalid get_left_elevation"
 stop
endif
if (inflow_count.le.0) then
 print *,"inflow_count invalid"
 stop
endif

do while ((inflow_time(last_inflow_index).gt.t).and. &
          (last_inflow_index.gt.1)) 
 last_inflow_index=last_inflow_index-1
enddo
do while ((inflow_time(last_inflow_index).lt.t).and. &
          (last_inflow_index.lt.inflow_count))
 last_inflow_index=last_inflow_index+1
enddo
if (t.le.inflow_time(1)) then
 elevation_left=inflow_elevation(1)
else if (t.ge.inflow_time(inflow_count)) then
 elevation_left=inflow_elevation(inflow_count)
else 

 if (t.le.inflow_time(last_inflow_index)) then
  k1=last_inflow_index-1
 else
  k1=last_inflow_index
 endif
 k2=k1+1
 elevation_left=inflow_elevation(k1)+ &
  ( (inflow_elevation(k2)-inflow_elevation(k1))/ &
    (inflow_time(k2)-inflow_time(k1)) )*(t-inflow_time(k1))
endif
elevation_left=elevation_left*100.0  ! convert to cm

return
end subroutine get_left_elevationIOWA


subroutine get_right_elevationIOWA(t,elevation_right)
use probcommon_module
IMPLICIT NONE
INTEGER_T k1,k2
REAL_T t,elevation_right

if (probtype.ne.110) then
 print *,"probtype invalid get_right_elevation"
 stop
endif
if (outflow_count.le.0) then
 print *,"outflow_count invalid"
 stop
endif

do while ((outflow_time(last_outflow_index).gt.t).and. &
          (last_outflow_index.gt.1)) 
 last_outflow_index=last_outflow_index-1
enddo
do while ((outflow_time(last_outflow_index).lt.t).and. &
          (last_outflow_index.lt.outflow_count))
 last_outflow_index=last_outflow_index+1
enddo
if (t.le.outflow_time(1)) then
 elevation_right=outflow_elevation(1)
else if (t.ge.outflow_time(outflow_count)) then
 elevation_right=outflow_elevation(outflow_count)
else 

 if (t.le.outflow_time(last_outflow_index)) then
  k1=last_outflow_index-1
 else
  k1=last_outflow_index
 endif
 k2=k1+1
 elevation_right=outflow_elevation(k1)+ &
  ( (outflow_elevation(k2)-outflow_elevation(k1))/ &
    (outflow_time(k2)-outflow_time(k1)) )*(t-outflow_time(k1))
endif
elevation_right=elevation_right*100.0  ! convert to cm

return
end subroutine get_right_elevationIOWA


subroutine get_bottom_elevation(x,elevation)
use probcommon_module
IMPLICIT NONE

real*8 x,elevation
real*8 h,l

 ! inlet x/H=-52 outlet x/H=44
 ! inlet x=-52H=594.36 cm  outlet x=44H=502.92 cm
 ! zhi=5H=57.15 cm
 ! (x/l)=(x/(2.5 h))=2/5 (x/h)
 ! at inlet, (x/l)=-52 * 2/5 = -20.8
 ! at outlet, (x/l)=44 * 2/5 = 17.6
 ! 0<z/h<5
 ! dx_min=0.3mm  dz_min=0.2mm ?
 ! -1<x/l<1
h=11.43  ! h=H  maximum bump height in centimeters
l=2.5*h
if (abs(x/l).le.1.0) then
 elevation=h*(1.0-2.0*((x/l)**2)+(x/l)**4)
else if (x/l.lt.0.0) then
 elevation=0.0
else if (x/l.gt.0.0) then
 elevation=0.0
else
 print *,"bust"
 stop
endif

return
end subroutine get_bottom_elevation

subroutine local_shallow_water_elevation(time,x,dist)
use probcommon_module
IMPLICIT NONE

REAL_T time,x,dist
REAL_T xright,xleft,shallow_tstop
REAL_T thetax,thetat,tgrid,xgrid
REAL_T delta_t_grid,delta_x_grid
REAL_T new_time
INTEGER_T igrid,jgrid

xleft=-594.36
xright=502.92
shallow_tstop=15.0
new_time=time+SHALLOW_TIME

delta_t_grid=shallow_tstop/SHALLOW_M
delta_x_grid=(xright-xleft)/SHALLOW_N

  !  NINT->round to nearest whole number

if (new_time.ge.shallow_tstop) then
 igrid=SHALLOW_M-1
else if (new_time.le.zero) then
 igrid=0
else
 igrid=NINT(new_time/delta_t_grid)-2
 if (igrid.lt.0) then
  igrid=0
 endif
 tgrid=igrid*delta_t_grid
 do while ((tgrid.lt.new_time).and.(igrid.lt.SHALLOW_M))
  igrid=igrid+1
  tgrid=igrid*delta_t_grid
 enddo
 igrid=igrid-1
endif

if (x.le.xleft) then
 jgrid=0
else if (x.ge.xright) then
 jgrid=SHALLOW_N-1
else
 jgrid=NINT((x-xleft)/delta_x_grid)-2
 if (jgrid.lt.0) then
  jgrid=0
 endif
 xgrid=xleft+jgrid*delta_x_grid
 do while ((xgrid.lt.x).and.(jgrid.lt.SHALLOW_N))
  jgrid=jgrid+1
  xgrid=xleft+jgrid*delta_x_grid
 enddo
 jgrid=jgrid-1
endif

xgrid=xleft+jgrid*delta_x_grid
tgrid=igrid*delta_t_grid
if (x.le.xgrid) then
 thetax=zero
else if (x.ge.xgrid+delta_x_grid) then
 thetax=one
else
 thetax=(x-xgrid)/delta_x_grid
endif

if (new_time.le.tgrid) then
 thetat=zero
else if (new_time.ge.tgrid+delta_t_grid) then
 thetat=one
else
 thetat=(new_time-tgrid)/delta_t_grid
endif

dist= &
 (one-thetax)*(one-thetat)*shallow_water_data(igrid,jgrid,1)+ &
 (thetax)*(one-thetat)*shallow_water_data(igrid,jgrid+1,1)+ &
 (thetax)*(thetat)*shallow_water_data(igrid+1,jgrid+1,1)+ &
 (one-thetax)*(thetat)*shallow_water_data(igrid+1,jgrid,1)

return
end subroutine local_shallow_water_elevation

subroutine local_shallow_water_velocity(time,x,vel)
use probcommon_module
IMPLICIT NONE

REAL_T time,x,vel
REAL_T xright,xleft,shallow_tstop
REAL_T thetax,thetat,tgrid,xgrid
REAL_T delta_t_grid,delta_x_grid
REAL_T new_time
INTEGER_T igrid,jgrid

xleft=-594.36
xright=502.92
shallow_tstop=15.0
new_time=time+SHALLOW_TIME

delta_t_grid=shallow_tstop/SHALLOW_M
delta_x_grid=(xright-xleft)/SHALLOW_N

  !  NINT->round to nearest whole number

if (new_time.ge.shallow_tstop) then
 igrid=SHALLOW_M-1
else if (new_time.le.zero) then
 igrid=0
else
 igrid=NINT(new_time/delta_t_grid)-2
 if (igrid.lt.0) then
  igrid=0
 endif
 tgrid=igrid*delta_t_grid
 do while ((tgrid.lt.new_time).and.(igrid.lt.SHALLOW_M))
  igrid=igrid+1
  tgrid=igrid*delta_t_grid
 enddo
 igrid=igrid-1
endif

if (x.le.xleft) then
 jgrid=0
else if (x.ge.xright) then
 jgrid=SHALLOW_N-1
else
 jgrid=NINT((x-xleft)/delta_x_grid)-2
 if (jgrid.lt.0) then
  jgrid=0
 endif
 xgrid=xleft+jgrid*delta_x_grid
 do while ((xgrid.lt.x).and.(jgrid.lt.SHALLOW_N))
  jgrid=jgrid+1
  xgrid=xleft+jgrid*delta_x_grid
 enddo
 jgrid=jgrid-1
endif

xgrid=xleft+jgrid*delta_x_grid
tgrid=igrid*delta_t_grid
if (x.le.xgrid) then
 thetax=zero
else if (x.ge.xgrid+delta_x_grid) then
 thetax=one
else
 thetax=(x-xgrid)/delta_x_grid
endif

if (new_time.le.tgrid) then
 thetat=zero
else if (new_time.ge.tgrid+delta_t_grid) then
 thetat=one
else
 thetat=(new_time-tgrid)/delta_t_grid
endif

vel= &
 (one-thetax)*(one-thetat)*shallow_water_data(igrid,jgrid,2)+ &
 (thetax)*(one-thetat)*shallow_water_data(igrid,jgrid+1,2)+ &
 (thetax)*(thetat)*shallow_water_data(igrid+1,jgrid+1,2)+ &
 (one-thetax)*(thetat)*shallow_water_data(igrid+1,jgrid,2)

return
end subroutine local_shallow_water_velocity




subroutine doit(problo,probhi,ncell,dx,tstop)
use probcommon_module
IMPLICIT NONE

real*8 problo,probhi,dx,tstop
integer ncell
real*8 SSold(-1:ncell)
real*8 SS(-1:ncell)
real*8 QQold(-1:ncell)
real*8 QQ(-1:ncell)
real*8 DD(-1:ncell)  ! bottom topography
real*8 xx(-1:ncell)
real*8 SSflux(0:ncell)
real*8 QQflux(0:ncell)
real*8 start_elevation
integer skip,nstep,i
real*8 local_gravity,time
real*8 elevation_right,elevation_left,u_right,u_left
real*8 maxu,maxc,den,mom,uu,cc,dt,lambda
real*8 denleft,denright,momleft,momright,pleft,pright
real*8 minz
integer icrit

integer igrid,jgrid  ! t index x index
REAL_T delta_t_grid,delta_x_grid,t_grid,x_grid
REAL_T thetax,thetat
integer hitgrid,last_index

delta_t_grid=tstop/SHALLOW_M
delta_x_grid=(probhi-problo)/SHALLOW_N

skip=2000

local_gravity=980.0
start_elevation=22.862

time=0.0
nstep=0

do i=-1,ncell
 xx(i)=problo+(i+0.5)*dx
 call get_bottom_elevation(xx(i),DD(i))
 QQ(i)=0.0
 SS(i)=start_elevation-DD(i)
enddo

do while (time.le.tstop-1.0D-10) 

 call get_right_elevationIOWA(time,elevation_right)
 call get_left_elevationIOWA(time,elevation_left)
 call get_right_velocityIOWA(time,u_right)
 call get_left_velocityIOWA(time,u_left)
 SS(-1)=elevation_left
 SS(ncell)=elevation_right
 QQ(-1)=elevation_left*u_left
 QQ(ncell)=elevation_right*u_right

 do i=-1,ncell
  SSold(i)=SS(i)
  QQold(i)=QQ(i)
 enddo

 maxu=0.0
 maxc=0.0 
 do i=0,ncell-1
  den=SS(i)
  mom=QQ(i)
  if (den.le.0.0) then
   print *,"density must be positive (subroutine doit shallow_water)"
   print *,"i,den ",i,den
   stop
  endif
  cc=sqrt(local_gravity*den)
  if (cc.gt.maxc) then
   maxc=cc
  endif
  uu=abs(mom/den)
  if (uu.gt.maxu) then
   maxu=abs(uu)
  endif
 enddo
 if (maxu+maxc.le.0.0) then
  print *,"must have acoustic waves"
  stop
 endif
 dt=0.8*dx/(maxu+maxc)
 if (time+dt.ge.tstop) then
  dt=tstop-time
 endif
 lambda=dt/dx

 do i=0,ncell
  denleft=SS(i-1)
  denright=SS(i)
  momleft=QQ(i-1)
  momright=QQ(i)
  pleft=0.5*local_gravity*(denleft**2)
  pright=0.5*local_gravity*(denright**2)
  SSflux(i)=0.5*(momleft+momright)-0.5*(denright-denleft)/lambda
  if (denright.le.0.0) then
   print *,"hydraulic section must be positive"
   stop
  endif
  if (denleft.le.0.0) then
   print *,"hydraulic section must be positive"
   stop
  endif
  QQflux(i)=0.5*(momleft**2/denleft+momright**2/denright+ &
    pleft+pright)-0.5*(momright-momleft)/lambda
 enddo
 do i=0,ncell-1
  SS(i)=SS(i)-lambda*(SSflux(i+1)-SSflux(i))
  QQ(i)=QQ(i)-lambda*(QQflux(i+1)-QQflux(i))- &
    lambda*SS(i)*local_gravity*0.5*(DD(i+1)-DD(i-1))
 enddo

 hitgrid=0
 do igrid=0,SHALLOW_M
  t_grid=igrid*delta_t_grid
  if (hitgrid.eq.0) then
   if ((time.le.t_grid).and.(time+dt.ge.t_grid)) then
    hitgrid=1
    last_index=-1
    do jgrid=0,SHALLOW_N
     x_grid=problo+jgrid*delta_x_grid
     do while (xx(last_index).le.x_grid)
      last_index=last_index+1
     enddo
     last_index=last_index-1
     if (t_grid.le.time) then
      thetat=zero
     else if (t_grid.ge.time+dt) then
      thetat=one
     else 
      thetat=(t_grid-time)/dt
     endif
     if (x_grid.le.xx(last_index)) then
      thetax=zero
     else if (x_grid.ge.xx(last_index+1)) then
      thetax=one
     else
      thetax=(x_grid-xx(last_index))/dx
     endif
     shallow_water_data(igrid,jgrid,1)= &
       (one-thetax)*(one-thetat)* &
          (DD(last_index)+SSold(last_index))+ &
       (one-thetax)*thetat* &
          (DD(last_index)+SS(last_index))+ &
       (one-thetat)*thetax* &
          (DD(last_index+1)+SSold(last_index+1))+ &
       thetat*thetax* &
         (SS(last_index+1)+DD(last_index+1))
     shallow_water_data(igrid,jgrid,2)= &
       (one-thetax)*(one-thetat)* &
       QQold(last_index)/SSold(last_index)+ &
       (one-thetax)*thetat* &
       QQ(last_index)/SS(last_index)+ &
       (one-thetat)*thetax* &
       QQold(last_index+1)/SSold(last_index+1)+ &
       thetat*thetax* &
       QQ(last_index+1)/SS(last_index+1)
    enddo  ! jgrid
   endif ! time<t_grid<time+dt
  endif ! hitgrid=0
 enddo ! igrid

 time=time+dt
 nstep=nstep+1

 if (nstep-skip*(nstep/skip).eq.0) then
  print *,"time=",time
  print *,"dt=",dt
  icrit=0
  minz=SS(0)+DD(0)
  do i=0,ncell-1
   if (SS(i)+DD(i).lt.minz) then
    minz=SS(i)+DD(i)
    icrit=i
   endif
  enddo
  print *,"icrit,x,minz ",icrit,xx(icrit),minz
 endif
enddo ! while time<tstop

return
end subroutine doit

subroutine shallow_water_solve()
use probcommon_module
IMPLICIT NONE

real*8 problo,probhi,tstop,dx
integer ncell

problo=-594.36
probhi=502.92

 ! 0.1 mm  (0.2mm in Fred Stern's paper)
 ! 100000
ncell=5000
dx=(probhi-problo)/ncell
tstop=15.0

call doit(problo,probhi,ncell,dx,tstop)

return
end subroutine shallow_water_solve

! AUTHOR: Dr. Yang Liu November 2020
subroutine smooth_init(center, r, Tsat, Tinf, x_in, Tout,delta)
! center, r:   center and radius of the circle interface
! Tsat:    saturation temp      Tinf: ambient temp
! x_in: coordinate in
! Tout: temperature out
implicit none

REAL_T,intent(in)  :: center(SDIM),x_in(SDIM)
REAL_T,intent(in)  :: r
REAL_T,intent(in)  :: Tinf,Tsat
REAL_T             :: phi
REAL_T,intent(in)  :: delta   ! size of smooth transition region
REAL_T,intent(out) :: Tout
INTEGER_T          :: i
REAL_T  :: H_local
REAL_T  :: phi_shift

phi=zero
do i=1,SDIM
 phi=phi+(x_in(i)-center(i))**2.0d0
enddo
phi=sqrt(phi)-r

! phi=0,  T=Tsat
! phi=delta,  T=Tinf

phi_shift=phi-half*delta
H_local=hs(phi_shift,half*delta)

Tout=Tinf*H_local+Tsat*(one-H_local)

end subroutine smooth_init

subroutine stress_from_strain( &
 im_elastic, & ! =1..nmat
 x_stress, & ! 1..sdim (x,y,z)
 dx, &  ! representative dx values 1..sdim
 gradu, &  ! dir_x (displace),dir_space
 xdisplace, &
 ydisplace, &
 DISP_TEN, &  ! dir_x (displace),dir_space
 hoop_22) ! if RZ, no "theta,theta" component, no place to put hoop_22
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: im_elastic
REAL_T, intent(in) :: x_stress(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: gradu(SDIM,SDIM)
REAL_T, intent(in) :: xdisplace
REAL_T, intent(in) :: ydisplace
REAL_T, intent(out) :: DISP_TEN(SDIM,SDIM)
REAL_T, intent(out) :: hoop_22
REAL_T :: gradu_local(SDIM,SDIM)
INTEGER_T :: dir_x,dir_space,dir_inner

REAL_T :: hoop_12
REAL_T strain_displacement(SDIM,SDIM)
REAL_T F(SDIM,SDIM)
REAL_T C(SDIM,SDIM)
REAL_T B(SDIM,SDIM)
REAL_T E(SDIM,SDIM)
REAL_T scale_factor
REAL_T Identity_comp,trace_E,trace_SD,bulk_modulus,lame_coefficient
INTEGER_T linear_elastic_model

if ((im_elastic.ge.1).and.(im_elastic.le.num_materials)) then
 ! do nothing
else
 print *,"im_elastic invalid"
 stop
endif

 ! gradient of the Displacement: GRAD X
do dir_x=1,SDIM  ! dir_x (displace)
do dir_space=1,SDIM
 gradu_local(dir_x,dir_space)=gradu(dir_x,dir_space)
enddo
enddo

do dir_space=1,SDIM
 if (dx(dir_space).gt.zero) then
  ! do nothing
 else
  print *,"dx invalid"
  stop
 endif
enddo

hoop_12=0.0d0
hoop_22=0.0d0
if (SDIM.eq.2) then
 if (levelrz.eq.0) then
  ! do nothing
 else if (levelrz.eq.1) then
  if (x_stress(1).gt.VOFTOL*dx(1)) then
   hoop_22=xdisplace/x_stress(1)  ! xdisplace/r
  else if (abs(x_stress(1)).le.VOFTOL*dx(1)) then
   hoop_22=zero
  else 
   print *,"xsten(0,1) invalid"
   stop
  endif
 else if (levelrz.eq.3) then
  if (x_stress(1).gt.VOFTOL*dx(1)) then
   hoop_12=-ydisplace/x_stress(1)  ! -ydisplace/r
   hoop_22=xdisplace/x_stress(1)  ! xdisplace/r
   do dir_x=1,SDIM
    gradu_local(dir_x,2)=gradu_local(dir_x,2)/x_stress(1)
   enddo
   gradu_local(1,2)=gradu_local(1,2)+hoop_12
   gradu_local(2,2)=gradu_local(2,2)+hoop_22
  else if (abs(x_stress(1)).le.VOFTOL*dx(1)) then
   hoop_12=zero
   hoop_22=zero
   do dir_x=1,SDIM
    gradu_local(dir_x,2)=zero
   enddo
  else 
   print *,"x_stress(1) invalid"
   stop
  endif
 else
  print *,"levelrz invalid"
  stop
 endif
else if (SDIM.eq.3) then
 if (levelrz.eq.0) then
  ! do nothing
 else if (levelrz.eq.3) then
  if (x_stress(1).gt.VOFTOL*dx(1)) then
   hoop_12=-ydisplace/x_stress(1)  ! -ydisplace/r
   hoop_22=xdisplace/x_stress(1)  ! xdisplace/r
   do dir_x=1,SDIM
    gradu_local(dir_x,2)=gradu_local(dir_x,2)/x_stress(1)
   enddo
   gradu_local(1,2)=gradu_local(1,2)+hoop_12
   gradu_local(2,2)=gradu_local(2,2)+hoop_22
  else if (abs(x_stress(1)).le.VOFTOL*dx(1)) then
   hoop_12=zero
   hoop_22=zero
   do dir_x=1,SDIM
    gradu_local(dir_x,2)=zero
   enddo
  else 
   print *,"x_stress(1) invalid"
   stop
  endif
 else
  print *,"levelrz invalid"
  stop
 endif
else
 print *,"dimension bust"
 stop
endif


scale_factor=zero

 ! gradu(i,j)=partial XD_{i}/partial x_j
 ! F=grad X + I
 ! strain_displacement=(grad X + grad X^T)/2
do dir_x=1,SDIM 
do dir_space=1,SDIM
 if (dir_x.eq.dir_space) then
  Identity_comp=one
 else
  Identity_comp=zero
 endif

 F(dir_x,dir_space)=gradu_local(dir_x,dir_space)+Identity_comp

 if (scale_factor.le.abs(F(dir_x,dir_space))) then
  scale_factor=abs(F(dir_x,dir_space))
 endif

 C(dir_x,dir_space)=zero
 B(dir_x,dir_space)=zero
  ! look for ``linear elasticity'' on wikipedia (eij)
 strain_displacement(dir_x,dir_space)=half* &
    (gradu_local(dir_x,dir_space)+ &
     gradu_local(dir_space,dir_x))
  ! isotropic
  ! Cijkl=K dij dkl + mu(dik djl + dil djk -(2/3)dij dkl)
  ! Cijkl ekl=K dij ekk + mu(dik djl ekl + dil djk ekl -
  !  (2/3)dij dkl ekl )=
  ! K dij ekk + mu(eij+eji-(2/3)dij ekk)=
  ! K dij ekk + 2mu(eij-(1/3)dij ekk)

enddo
enddo

if (scale_factor.lt.one) then
 scale_factor=one
endif
scale_factor=scale_factor*scale_factor

 ! F=grad X + I
 ! C=F^T F = right cauchy green tensor
 ! E=(1/2)*(C-I)  Green Lagrange strain tensor
do dir_x=1,SDIM 
do dir_space=1,SDIM
 do dir_inner=1,SDIM
   ! C=F^T F
  C(dir_x,dir_space)=C(dir_x,dir_space)+ &
          F(dir_inner,dir_x)*F(dir_inner,dir_space)
   ! B=F F^T
  B(dir_x,dir_space)=B(dir_x,dir_space)+ &
          F(dir_x,dir_inner)*F(dir_space,dir_inner)
 enddo
enddo
enddo
do dir_x=1,SDIM 
do dir_space=1,SDIM
 if (abs(C(dir_x,dir_space)-C(dir_space,dir_x)).le. &
     1.0D-5*scale_factor) then
  ! do nothing
 else
  print *,"scale_factor = ",scale_factor
  print *,"x=",x_stress(1),x_stress(2),x_stress(SDIM)
  print *,"dir_x,dir_space ",dir_x,dir_space
  print *,"C(dir_x,dir_space)=",C(dir_x,dir_space)
  print *,"C(dir_space,dir_x)=",C(dir_space,dir_x)
  print *,"expecting C^T=C"
  stop
 endif 
 if (abs(B(dir_x,dir_space)-B(dir_space,dir_x)).le. &
     1.0D-5*scale_factor) then
  ! do nothing
 else
  print *,"scale_factor = ",scale_factor
  print *,"x=",x_stress(1),x_stress(2),x_stress(SDIM)
  print *,"dir_x,dir_space ",dir_x,dir_space
  print *,"B(dir_x,dir_space)=",B(dir_x,dir_space)
  print *,"B(dir_space,dir_x)=",B(dir_space,dir_x)
  print *,"expecting B^T=B"
  stop
 endif

enddo
enddo
trace_E=zero
trace_SD=zero ! trace of the strain displacement
do dir_x=1,SDIM 
do dir_space=1,SDIM
 if (dir_x.eq.dir_space) then
  Identity_comp=one
 else
  Identity_comp=zero
 endif
  ! grad X is synonomous with grad u (gradient of displacement vector)
  ! (note: grad V=gradient of velocity vector)
  ! F=grad X + I
  ! C=F^T F=right cauchy green tensor
  ! strain_displacement=(1/2)(grad u + (grad u)^T) = eps_ij
  ! E=(C-I)/2=( (grad u + I)^T(grad u +I) - I)/2=eps_ij +
  !  grad u^T gradu/2
 E(dir_x,dir_space)=half*(C(dir_x,dir_space)-Identity_comp)
 trace_E=trace_E+Identity_comp*E(dir_x,dir_space)
 trace_SD=trace_SD+Identity_comp*strain_displacement(dir_x,dir_space)
enddo
enddo 
 ! E=(C-I)/2=(F^T F -I)/2=( (grad X^T+I)(grad X+I)-I )/2=
 ! (1/2)(grad X +grad X^T + grad X^T * grad X)
 ! Sigma=2 mu_s E + lambda Tr(E) I
 ! structure force is div Sigma=div mu_s (Sigma/mu_s)
 ! Richter, JCP, 2013
 ! MKS: mu_s=1E+4   lambda=4E+4  density_struct=1E+3
 ! bulk modulus units:
 ! steel 160 giga Pa=160 * 1E+9  N/m^2
 ! N/m^2=kg m/s^2 / m^2 = kg /(m s^2) = 1000/100 g/(cm s^2)
do dir_x=1,SDIM 
do dir_space=1,SDIM
 if (dir_x.eq.dir_space) then
  Identity_comp=one
 else
  Identity_comp=zero
 endif

 bulk_modulus=fort_elastic_viscosity(im_elastic)
 lame_coefficient=fort_lame_coefficient(im_elastic)
 linear_elastic_model=fort_linear_elastic_model(im_elastic)

 if (bulk_modulus.gt.zero) then

  if (linear_elastic_model.eq.1) then
    ! only valid for small deformations:
    !  strain displacement=SD= (1/2)(grad XD + grad XD^T)=( (F + F^T)/2 - I )
    !  2 mu_S SD + lambda TRACE(SD)
   DISP_TEN(dir_x,dir_space)=( &
       two*bulk_modulus*strain_displacement(dir_x,dir_space)+ &
          lame_coefficient*trace_SD*Identity_comp)/bulk_modulus
  else if (linear_elastic_model.eq.0) then
    ! E=(C-I)/2  C=F^T F = right cauchy green tensor
   DISP_TEN(dir_x,dir_space)=(two*bulk_modulus*E(dir_x,dir_space)+ &
          lame_coefficient*trace_E*Identity_comp)/bulk_modulus
  else
   print *,"linear_elastic_model invalid"
   stop
  endif
 else 
  print *,"bulk_modulus invalid ",bulk_modulus
  stop
 endif
enddo
enddo

end subroutine stress_from_strain

subroutine get_primary_material(LS,nmat,im_primary)
use probcommon_module

IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: LS(nmat)
INTEGER_T, intent(out) :: im_primary
INTEGER_T im,imtest
INTEGER_T tessellate
INTEGER_T is_rigid_local(nmat)

tessellate=0

do im=1,nmat
 is_rigid_local(im)=is_rigid(nmat,im)
 if (tessellate.eq.2) then
  is_rigid_local(im)=0
  print *,"expecting tessellate==0"
  stop
 else if (tessellate.eq.0) then
  ! do nothing
 else if (tessellate.eq.1) then
  print *,"expecting tessellate==0"
  stop
 else if (tessellate.eq.3) then
  print *,"expecting tessellate==0"
  stop
 else
  print *,"tessellate invalid38"
  stop
 endif
enddo ! im=1..nmat

if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
 print *,"nmat invalid get_primary_material"
 print *,"nmat= ",nmat
 stop
endif

im_primary=0
do im=1,nmat
 if (is_rigid_local(im).eq.1) then
  if (LS(im).ge.zero) then
   if (im_primary.ne.0) then
    print *,"cannot have two rigid materials in same place"
    do imtest=1,nmat
     print *,"imtest,LS(imtest) ",imtest,LS(imtest)
    enddo
    stop
   endif
   im_primary=im
  else if (LS(im).le.zero) then
   ! do nothing
  else
   print *,"LS bust"
   stop
  endif
 else if (is_rigid_local(im).eq.0) then
  ! do nothing
 else
  print *,"is_rigid invalid"
  stop
 endif
enddo !im=1..nmat

if (im_primary.eq.0) then

 do im=1,nmat
   if (im_primary.eq.0) then
    im_primary=im
   else if ((im_primary.ge.1).and.(im_primary.lt.im)) then
    if (LS(im).gt.LS(im_primary)) then
     im_primary=im
    else if (LS(im).le.LS(im_primary)) then
     ! do nothing
    else
     print *,"LS bust"
     stop
    endif
   else
    print *,"im_primary invalid"
    stop
   endif
 enddo !im=1..nmat

else if (is_rigid_local(im_primary).eq.1) then
 ! do nothing
else
 print *,"is_rigid or im_primary invalid"
 stop
endif

end subroutine get_primary_material


subroutine tensor_Heaviside( &
    dxmin, &
    im_parm, & ! 0..nmat-1
    mask1,mask2, &
    LS1,LS2, &
    HVAL)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: dxmin
INTEGER_T, intent(in) :: im_parm ! 0..nmat-1
REAL_T, intent(out) :: HVAL
INTEGER_T, intent(in) :: mask1,mask2
REAL_T, intent(in) :: LS1(num_materials)
REAL_T, intent(in) :: LS2(num_materials)
INTEGER_T im1,im2
REAL_T LS_avg

if ((im_parm.ge.0).and.(im_parm.lt.num_materials)) then
 ! do nothing
else
 print *,"im_parm invalid"
 stop
endif

if ((mask1.eq.0).or.(mask2.eq.0)) then
 HVAL=zero
else if ((mask1.eq.1).and.(mask2.eq.1)) then
 call get_primary_material(LS1,num_materials,im1)
 call get_primary_material(LS2,num_materials,im2)
 if ((im1.eq.im_parm+1).and.(im2.eq.im_parm+1)) then
  LS_avg=half*(LS1(im_parm+1)+LS2(im_parm+1))-dxmin
  HVAL=hs(LS_avg,dxmin)
 else
  print *,"im1 or im2 invalid"
  stop
 endif
else
 print *,"mask1 or mask2 invalid"
 stop
endif

end subroutine tensor_Heaviside


INTEGER_T function project_option_is_validF(project_option) 
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: project_option

if (project_option_momeqnF(project_option).eq.1) then
 project_option_is_validF=1
else if (project_option_momeqnF(project_option).eq.0) then
 project_option_is_validF=1
else
 print *,"project_option not valid"
 stop
 project_option_is_validF=0
endif

end function project_option_is_validF

INTEGER_T function project_option_momeqnF(project_option) 
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: project_option

 if ((project_option.eq.0).or. & ! regular project
     (project_option.eq.1).or. & ! initial project
     (project_option.eq.11).or.& ! FSI_material_exists (last project)
     (project_option.eq.12).or.& ! pressure extrapolation
     (project_option.eq.3)) then ! viscosity
  project_option_momeqnF=1
 else if ((project_option.eq.2).or. & ! thermal diffusion
          ((project_option.ge.100).and. & ! species
           (project_option.lt.100+num_species_var)).or. &
          (project_option.eq.200)) then ! smooth temperature
  project_option_momeqnF=0
 else
  print *,"project_option invalid"
  stop
  project_option_momeqnF=0
 endif

end function project_option_momeqnF


INTEGER_T function project_option_singular_possibleF(project_option) 
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: project_option

 if ((project_option.eq.0).or. & ! regular project
     (project_option.eq.1).or. & ! initial project
     (project_option.eq.11).or. & !FSI_material_exists (last project)
     (project_option.eq.12)) then ! pressure extension
  project_option_singular_possibleF=1
 else if ((project_option.eq.2).or. & ! thermal diffusion
          (project_option.eq.3).or. & ! viscosity
          ((project_option.ge.100).and. &
           (project_option.lt.100+num_species_var)).or. & !species
          (project_option.eq.200)) then !smoothing
  project_option_singular_possibleF=0
 else
  print *,"project_option invalid"
  stop
  project_option_singular_possibleF=0
 endif

end function project_option_singular_possibleF

INTEGER_T function project_option_olddata_neededF(project_option) 
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: project_option

 if ((project_option.eq.0).or. & ! regular project
     (project_option.eq.1).or. & ! initial project
     (project_option.eq.11).or. & ! FSI_material_exists (last project)
     (project_option.eq.12)) then ! pressure extension
  project_option_olddata_neededF=0
 else if ((project_option.eq.2).or. & ! thermal diffusion
          (project_option.eq.3).or. & ! viscosity
          ((project_option.ge.100).and. &
           (project_option.lt.100+num_species_var)).or. & !species
          (project_option.eq.200)) then !smoothing
  project_option_olddata_neededF=1
 else 
  print *,"project_option invalid"
  stop
  project_option_olddata_neededF=0
 endif

end function project_option_olddata_neededF

INTEGER_T function project_option_pressureF(project_option)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: project_option

 if ((project_option.eq.0).or. &
     (project_option.eq.1).or. &
     (project_option.eq.12)) then ! pressure extrapolation
  project_option_pressureF=1
 else if ((project_option.eq.11).or. & ! FSI_material_exists (last project)
          (project_option.eq.2).or. &  ! temperature
          (project_option.eq.3).or. &  ! viscosity
          ((project_option.ge.100).and. &
           (project_option.lt.100+num_species_var)).or. & ! species
          (project_option.eq.200)) then ! smoothing of temperature
  project_option_pressureF=0
 else
  print *,"project_option invalid"
  stop
  project_option_pressureF=0
 endif 

end function project_option_pressureF


INTEGER_T function project_option_needs_scalingF(project_option)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: project_option

 if ((project_option.eq.0).or. & ! regular project
     (project_option.eq.11).or. & !FSI_material_exists last project
     (project_option.eq.12)) then ! pressure extrapolation
  project_option_needs_scalingF=1
 else if ((project_option.eq.1).or. & ! initial project
          (project_option.eq.2).or. &  ! temperature
          (project_option.eq.3).or. &  ! viscosity
          ((project_option.ge.100).and. &
           (project_option.lt.100+num_species_var)).or. & ! species
          (project_option.eq.200)) then ! smoothing of temperature
  project_option_needs_scalingF=0
 else
  print *,"project_option invalid"
  stop
  project_option_needs_scalingF=0
 endif 

end function project_option_needs_scalingF


INTEGER_T function project_option_projectionF(project_option)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: project_option

 if ((project_option.eq.0).or. & ! regular project
     (project_option.eq.11).or. & !FSI_material_exists last project
     (project_option.eq.1)) then ! initial_project
  project_option_projectionF=1
 else if ((project_option.eq.12).or. & ! pressure extrapolation
          (project_option.eq.2).or. &  ! temperature
          (project_option.eq.3).or. &  ! viscosity
          ((project_option.ge.100).and. &
           (project_option.lt.100+num_species_var)).or. & ! species
          (project_option.eq.200)) then ! smoothing of temperature
  project_option_projectionF=0
 else
  print *,"project_option invalid"
  stop
  project_option_projectionF=0
 endif 

end function project_option_projectionF

INTEGER_T function is_GFM_freezing_modelF(freezing_model) 
IMPLICIT NONE

INTEGER_T, intent(in) :: freezing_model

 if ((freezing_model.eq.0).or. &   !fully saturated
     (freezing_model.eq.5).or. &   !stefan model evap or condensation
     (freezing_model.eq.6)) then   !Palmore and Desjardins
  is_GFM_freezing_modelF=1
 else if (is_valid_freezing_modelF(freezing_model).eq.1) then
  is_GFM_freezing_modelF=0
 else
  print *,"freezing_model bust(F)"
  stop
  is_GFM_freezing_modelF=0
 endif

end function is_GFM_freezing_modelF 

INTEGER_T function is_hydrate_freezing_modelF(freezing_model) 
IMPLICIT NONE

INTEGER_T, intent(in) :: freezing_model

 if (freezing_model.eq.2) then
  is_hydrate_freezing_modelF=1
 else if (is_valid_freezing_modelF(freezing_model).eq.1) then
  is_hydrate_freezing_modelF=0
 else
  print *,"freezing_model invalid (F)"
  stop
  is_hydrate_freezing_modelF=0
 endif
end function is_hydrate_freezing_modelF

INTEGER_T function is_valid_freezing_modelF(freezing_model) 
IMPLICIT NONE

INTEGER_T, intent(in) :: freezing_model

 if ((freezing_model.eq.4).or. & !Tannasawa or Schrage 
     (freezing_model.eq.5).or. & !Stefan model evaporation or condensation
     (freezing_model.eq.6).or. & !Palmore and Desjardins
     (freezing_model.eq.7)) then !cavitation
  is_valid_freezing_modelF=1
 else if ((freezing_model.eq.0).or. & !Energy jump model
          (freezing_model.eq.1).or. & !source term
          (freezing_model.eq.2).or. & !hydrate
          (freezing_model.eq.3)) then !wildfire
  is_valid_freezing_modelF=1
 else 
  print *,"freezing_model invalid (F)"
  stop
  is_valid_freezing_modelF=0
 endif

end function is_valid_freezing_modelF

INTEGER_T function is_multi_component_evapF(freezing_model, &
   evap_flag,latent_heat) 
IMPLICIT NONE

INTEGER_T, intent(in) :: freezing_model
INTEGER_T, intent(in) :: evap_flag
REAL_T, intent(in) :: latent_heat

 if (latent_heat.eq.zero) then
  is_multi_component_evapF=0
 else if (latent_heat.ne.zero) then

  if ((freezing_model.eq.4).or. & !Tannasawa or Schrage 
      (freezing_model.eq.5).or. & !Stefan model evaporation or condensation
      (freezing_model.eq.6).or. & !Palmore and Desjardins
      (freezing_model.eq.7)) then !cavitation

   if (evap_flag.eq.0) then !Palmore and Desjardins
    is_multi_component_evapF=1
   else if ((evap_flag.eq.1).or. & !Tanasawa
            (evap_flag.eq.2).or. & !Schrage
            (evap_flag.eq.3).or. & !Kassemi
            (evap_flag.eq.4)) then !Tanguy recommendation.
    is_multi_component_evapF=0
   else
    print *,"evap_flag invalid (F) "
    stop
    is_multi_component_evapF=0
   endif

  else if (freezing_model.eq.2) then !hydrate

   is_multi_component_evapF=1

  else if ((freezing_model.eq.0).or. & !Energy jump model
           (freezing_model.eq.1).or. & !source term
           (freezing_model.eq.3)) then !wildfire
   is_multi_component_evapF=0
  else
   print *,"freezing_model invalid (F) "
   stop
   is_multi_component_evapF=0
  endif
 else
  print *,"latent_heat invalid (F)"
  stop
  is_multi_component_evapF=0
 endif

end function is_multi_component_evapF


end module global_utility_module

#undef STANDALONE
