#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif


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

#include "EXTRAP_COMP.H"
#include "TOWER_MF_CODES.H"


      module LagrangeInterpolationPolynomial
      use amrex_fort_module, only : amrex_real
      IMPLICIT NONE

      CONTAINS
      SUBROUTINE AlmostEqual(a,b,var1)
      IMPLICIT NONE
      real(amrex_real), INTENT(IN) :: a,b
      LOGICAL, INTENT(OUT)      :: var1    !!!!!!!!!! logical: true and false
      real(amrex_real)             :: tol

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
      integer, INTENT(IN)                           :: N
      real(amrex_real), DIMENSION(0:N), INTENT(IN)     :: Phi       !values at the nodes
      real(amrex_real), DIMENSION(0:N,0:N), INTENT(IN) :: D         !polynomial derivative matrix that has been pre-computed
      real(amrex_real), DIMENSION(0:N), INTENT(OUT)    :: Phi_prime !derivative of the interpolent at the nodes      
      integer                                       :: i,j
      real(amrex_real)                                 :: t

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
      integer, INTENT(IN)                        :: N
      real(amrex_real), DIMENSION(0:N), INTENT(IN)  :: x
      real(amrex_real), DIMENSION(0:N), INTENT(OUT) :: w
      integer                                    :: j,k
   
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
      integer, INTENT(IN)                        :: M
      real(amrex_real), DIMENSION(0:M), INTENT(IN)  :: x,w,f
      real(amrex_real), INTENT(IN)                  :: y
      real(amrex_real), INTENT(OUT)                 :: p

      integer                                    :: j
      real(amrex_real)                              :: numerator,denominator,t
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
      integer,INTENT(IN) :: N
      real(amrex_real), INTENT(IN) :: xpt
      real(amrex_real), DIMENSION(0:N), INTENT(IN)  :: x,w
      real(amrex_real), DIMENSION(0:N), INTENT(OUT) :: l

      integer :: j,jmatch
      real(amrex_real)  :: s,t
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
      integer,INTENT(IN)                         :: M
      real(amrex_real), DIMENSION(0:M), INTENT(IN)  :: x,w,f
      real(amrex_real), INTENT(IN)                  :: y      
      real(amrex_real), INTENT(OUT)                 :: p1
      integer                                    :: i,j
      real(amrex_real)                              :: numerator,denominator,t,p
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
      integer,INTENT(IN)                         :: N
      real(amrex_real), DIMENSION(0:N), INTENT(IN)  :: x,w
      real(amrex_real), DIMENSION(0:N,0:N), INTENT(OUT) :: DL
      integer                                    :: i,j

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
       use amrex_fort_module, only : amrex_real
       IMPLICIT NONE
   
!___________________________________________________________________


      real(amrex_real), dimension(:,:,:), allocatable :: cache_gauss 
      real(amrex_real), dimension(:,:,:), allocatable :: cache_gauss_lobatto 
      real(amrex_real), dimension(:,:,:), allocatable :: cache_gauss_w 
      real(amrex_real), dimension(:,:,:), allocatable :: cache_gauss_lobatto_w 

      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_wMATGL
      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_wMAT
      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_w_right
      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_w_left

      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_wMAT_extend
      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_wRT
      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_wLT
      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_wLRT

      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_wRT_EXT
      real(amrex_real) , dimension(:,:,:,:), allocatable :: cache_wLT_EXT

      integer, PARAMETER :: SPTYPE = 0  ! Legendre
      integer, PARAMETER :: TMTYPE = 0  ! Legendre
      integer, PARAMETER :: GQTYPE = 0  ! Legendre

CONTAINS

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 1
! legendre polynomial & derivative
! Algorithm 22
      SUBROUTINE LegendrePolyAndDeri(k,x,P,dl)
     !k=the degree of Legendre polynomial
     !P=Legendre polynomial 
     !dl=its derivative
      IMPLICIT NONE 
      integer,INTENT(IN)          :: k
      real(amrex_real), INTENT(IN)   :: x
      real(amrex_real), INTENT(OUT)  :: P,dl
      integer              :: j
      real(amrex_real)               :: P1,P2,dl1,dl2

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
      integer,INTENT(IN)          :: n
      real(amrex_real), DIMENSION(0:n), INTENT(OUT)  :: y,weight
      integer                     :: j,k
      real(amrex_real)               :: P1,delta
      real(amrex_real)               :: dP1,localPI
      real(amrex_real), PARAMETER :: stub_zero=zero

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
            call LegendrePolyAndDeri(n+1,stub_zero,P1,dP1)
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
      integer,INTENT(IN)          :: n
      real(amrex_real), INTENT(IN)   :: x
      real(amrex_real), INTENT(OUT)  :: q,dq,ln
      integer                     :: k
      real(amrex_real)               :: l1,l2,l3,dl1,dl2,dl,dl3
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

      integer n1,n2,i,j
      real(amrex_real) y1(0:n1),y2(0:n2)
      real(amrex_real) w(0:n1-1,0:n2-1),wt,low,high

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

      integer n1,n2,i,j
      real(amrex_real) y1(0:n1),y2(0:n2)
      real(amrex_real) y1ext(0:n1+2),y2ext(0:n2+2)
      real(amrex_real) w(0:n1+1,0:n2+1),wt,low,high

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
      integer,INTENT(IN)          :: n
      real(amrex_real), DIMENSION(0:n), INTENT(OUT)  :: y,weight
      integer                     :: j,k
      real(amrex_real)               :: q,dq,ln,delta,localPI
      real(amrex_real), PARAMETER :: stub_zero=zero

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
           call qAndLEvaluation(n,stub_zero,q,dq,ln) 
           k=int(n/two)
           y(k)=zero
           weight(k)=two/(n*(n+1)*ln**2)
        end if
      RETURN
      END SUBROUTINE LegendreGaussLobattoNodesAndWeights 
    
      subroutine do_polyinterp(order_r,w,data,sum)
      IMPLICIT NONE

      integer i,order_r
      real(amrex_real) w(0:order_r)
      real(amrex_real) data(0:order_r)
      real(amrex_real) sum

      sum=zero
      do i=0,order_r
       sum=sum+data(i)*w(i)
      enddo

      return
      end subroutine do_polyinterp
       

        ! r is the order of the polynomial interpolant
      subroutine polyinterp_weights(order_r,x,w,xtarget)
      IMPLICIT NONE

      integer order_r,i,j
      real(amrex_real) xtarget,lag
      real(amrex_real) x(0:order_r)
      real(amrex_real) w(0:order_r)

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
      integer,INTENT(IN)          :: n
      real(amrex_real), DIMENSION(0:n), INTENT(OUT)  :: y,weight
      integer                     :: j,k
      real(amrex_real)               :: localPI
      real(amrex_real) w(0:n)
      real(amrex_real) yL(0:n)
      real(amrex_real) wL(0:n)

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
      integer,INTENT(IN)          :: n
      real(amrex_real), DIMENSION(0:n), INTENT(OUT)  :: y,weight
      integer                     :: j,k
      real(amrex_real)               :: localPI
      real(amrex_real) w(0:n)
      real(amrex_real) yL(0:n)
      real(amrex_real) wL(0:n)

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

      integer order_r,i,j,bfact,bfactC
      real(amrex_real) xtarget,lag
      real(amrex_real) x(0:order_r)
      real(amrex_real) w(0:order_r)

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

      integer order_r,i,bfact
      real(amrex_real) xtarget
      real(amrex_real) x(0:order_r)
      real(amrex_real) w(0:order_r)
      real(amrex_real) dx,totalw


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

      integer order_r,i,j,l,m
      real(amrex_real) x(0:order_r)
      real(amrex_real) w(0:order_r,0:order_r)
      real(amrex_real) sum

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

      integer order_r1,order_r2,i
      real(amrex_real) y1(0:order_r1)
      real(amrex_real) y2(0:order_r2)
      real(amrex_real) data1(0:order_r1)
      real(amrex_real) w(0:order_r1)
      real(amrex_real) data2(0:order_r2)

      do i=0,order_r2
       call polyinterp_weights(order_r1,y1,w,y2(i))
       call do_polyinterp(order_r1,w,data1,data2(i))
      enddo

      return
      end subroutine poly_change_basis

      subroutine deriv_change_basis(order_r1,order_r2,data,datader, &
        wMAT,y1,y2,dx_element)
      IMPLICIT NONE

      integer order_r1,order_r2
      real(amrex_real) y1(0:order_r1)
      real(amrex_real) y2(0:order_r2)
      real(amrex_real) data(0:order_r1)
      real(amrex_real) datader1(0:order_r1)
      real(amrex_real) datader(0:order_r2)
      real(amrex_real) wMAT(0:order_r1,0:order_r1)
      integer i,j
      real(amrex_real) sum,dx_element
   
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
       ! called from PROB_CPP_PARMS.F90 (order_r=18): fort_override
       ! fort_override is called from 
       ! NavierStokes.cpp: void fortran_parameters
       ! fortran_parameters is called from main.cpp prior to:
       ! 1. AmrCore* amrptr=new AmrCore()
       ! 2. amrptr->init(strt_time,stop_time)
      subroutine init_cache(order_r)
      IMPLICIT NONE

      integer order_r
      real(amrex_real) yleft(0:order_r)
      real(amrex_real) yright(0:order_r)
      real(amrex_real) y(0:order_r)
      real(amrex_real) yGL(0:order_r)
      real(amrex_real) y_extend(0:order_r+1)
      real(amrex_real) yGL_extend(0:order_r+1)
      real(amrex_real) yLT(0:order_r+1)
      real(amrex_real) yRT(0:order_r+1)
      real(amrex_real) yLRT(0:order_r+1)

      real(amrex_real) yLRTextrap(0:order_r-1)
      real(amrex_real) yLTextrap(0:order_r)
      real(amrex_real) yRTextrap(0:order_r)

      real(amrex_real) yLTextrapEXT(0:order_r)
      real(amrex_real) yRTextrapEXT(0:order_r)

      real(amrex_real), dimension(:,:), allocatable :: deriv_matrix
      real(amrex_real), dimension(:), allocatable :: tempx,tempw

      integer i,j,typ,i1,j1

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
         ! exterior BC left, extend right
        yLT(0)=-one
        yLT(i+1)=y_extend(i+1)
         ! exterior BC right, extend left
        yRT(0)=y_extend(0)
        yRT(i+1)=one
         ! exterior BC left and right
        yLRT(0)=-one
        yLRT(i+1)=one
        do i1=1,i
         yLT(i1)=y(i1-1)
         yRT(i1)=y(i1-1)
         yLRT(i1)=y(i1-1)
        enddo ! i1

         ! extrap left, extend right
        yLTextrap(i)=y_extend(i+1)
         ! extrap left, exterior BC right
        yLTextrapEXT(i)=one
        do i1=0,i-1
         yLTextrap(i1)=y(i1)
         yLTextrapEXT(i1)=y(i1)
        enddo

         ! extrap right, extend left
        yRTextrap(0)=y_extend(0)
         ! extrap right, exterior BC left
        yRTextrapEXT(0)=-one
        do i1=0,i-1
         yRTextrap(i1+1)=y(i1)
         yRTextrapEXT(i1+1)=y(i1)
        enddo

         ! extrap left and right
        do i1=0,i-1
         yLRTextrap(i1)=y(i1)
        enddo

        allocate(deriv_matrix(0:i,0:i))
        call polyinterp_Dmatrix(i,yRTextrapEXT,deriv_matrix)
        do i1=0,i
        do j1=0,i
          ! extrap right, exterior BC left
         cache_wLT_EXT(i,i1,j1,typ)=deriv_matrix(i1,j1) 
        enddo
        enddo
        deallocate(deriv_matrix)

        allocate(deriv_matrix(0:i,0:i))
        call polyinterp_Dmatrix(i,yLTextrapEXT,deriv_matrix)
        do i1=0,i
        do j1=0,i
          ! Extrap left, exterior BC right
         cache_wRT_EXT(i,i1,j1,typ)=deriv_matrix(i1,j1) 
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
      use probcommon_module
      IMPLICIT NONE
 
      integer typ,order_r,i,p,rstart,rend
      real(amrex_real), dimension(:), allocatable :: y,weight,data
      real(amrex_real), dimension(:), allocatable :: w,datader
      real(amrex_real), dimension(:,:), allocatable :: wMAT
      real(amrex_real) exact,sum,xtarget,dx_element

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
         if (abs(sum-exact).le.EPS_12_4) then
          !do nothing
         else
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
          if (abs(sum-exact).le.EPS_12_4) then
           !do nothing
          else
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
          if (abs(sum-exact).le.EPS_10_3) then
           !do nothing
          else
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
use amrex_fort_module, only : amrex_real

implicit none

      type nucleation_parm_type_inout
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: LSnew
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: Snew
      end type nucleation_parm_type_inout

real(amrex_real) :: MOF_PI=zero

contains

      subroutine inverse_stress_index(one_dim_index,i,j)
      implicit none

      integer, INTENT(out) :: one_dim_index
      integer, INTENT(in) :: i
      integer, INTENT(in) :: j

      if (((i.eq.2).and.(j.eq.3)).or. &
          ((i.eq.3).and.(j.eq.2))) then
       one_dim_index=6
      else if (((i.eq.1).and.(j.eq.3)).or. &
               ((i.eq.3).and.(j.eq.1))) then
       one_dim_index=5
      else if (((i.eq.1).and.(j.eq.2)).or. &
               ((i.eq.2).and.(j.eq.1))) then
       one_dim_index=2
      else if ((i.eq.3).and.(j.eq.3)) then
       one_dim_index=4
      else if ((i.eq.2).and.(j.eq.2)) then
       one_dim_index=3
      else if ((i.eq.1).and.(j.eq.1)) then
       one_dim_index=1
      else
       print *,"i or j invalid: ",i,j
       stop
      endif

      end subroutine inverse_stress_index

      subroutine stress_index(one_dim_index,i,j)
      implicit none

      integer, INTENT(in) :: one_dim_index
      integer, INTENT(out) :: i
      integer, INTENT(out) :: j

      if (one_dim_index.eq.1) then
       i=1
       j=1
      else if (one_dim_index.eq.2) then
       i=1
       j=2
      else if (one_dim_index.eq.3) then
       i=2
       j=2
      else if (one_dim_index.eq.4) then
       i=3
       j=3
      else if (one_dim_index.eq.5) then
       i=1
       j=3
      else if (one_dim_index.eq.6) then
       i=2
       j=3
      else
       print *,"one_dim_index invalid: ",one_dim_index
       stop
      endif

      end subroutine stress_index

      function wallfunc(u_tau,u,y,K,B,rho_w,mu_w) result(f)
      implicit none
        !incompressible wall function, to be called in 
        !Newton's method, solving for friction 
        !velocity u_tau at the image point
        !u_tau: friction velocity, u: velocity parallel to wall, 
        !y: coord direction normal to wall
       real(amrex_real), INTENT(in) :: u_tau, u, y 
        !rho_w: wall density, mu_w: wall molecular viscosity 
       real(amrex_real), INTENT(in) :: K, B, rho_w, mu_w 
       real(amrex_real) :: f
       real(amrex_real) :: u_plus, y_plus

       if (u_tau.eq.zero) then
        print *,"u_tau.eq.zero"
        stop
       else if (u_tau.ne.zero) then
        ! do nothing
       else
        print *,"u_tau is NaN"
        stop
       endif
       if (mu_w.eq.zero) then
        print *,"mu_w.eq.zero"
        stop
       else if (mu_w.ne.zero) then
        ! do nothing
       else
        print *,"mu_w is NaN"
        stop
       endif
       if (rho_w.gt.zero) then
        ! do nothing
       else
        print *,"rho_w invalid"
        stop
       endif
       if ((K.gt.zero).and.(B.gt.zero)) then
        ! do nothing
       else
        print *,"K or B invalid"
        stop
       endif
       if (y.gt.zero) then
        ! do nothing
       else
        print *,"y invalid"
        stop
       endif
       if ((u.le.zero).or.(u.ge.zero)) then
        ! do nothing
       else
        print *,"u is NaN"
        stop
       endif

       u_plus = u/u_tau
       y_plus = rho_w*u_tau*y/mu_w

       f = -y_plus + u_plus + exp(-K*B)* &
        ( exp(K*u_plus)-one-K*u_plus-half*(K*u_plus)**2-(K*u_plus)**3/six )
      end function wallfunc
  
      function wallfuncderiv(u_tau,u,y,K,B,rho_w,mu_w) result(fprime)
      implicit none
        !derivative of incompressible wall function w.r.t. u_tau, to be 
        !called in Newton's method
       real(amrex_real), INTENT(in) :: u_tau, u, y 
       real(amrex_real), INTENT(in) :: K, B, rho_w, mu_w
       real(amrex_real) :: fprime
       real(amrex_real) :: u_plus
   
       if (u_tau.eq.zero) then
        print *,"u_tau cannot be 0"
        stop
       endif 
       u_plus = u/u_tau
    
       fprime = -rho_w*y/mu_w - u_plus/u_tau + &
         exp(-K*B)*( -(K*u_plus/u_tau)*exp(K*u_plus)+ &
                     K*u_plus/u_tau+(K*u_plus)**2/u_tau+ &
                     (K*u_plus)**3/(two*u_tau) )
      end function wallfuncderiv
      
      subroutine wallfunc_newtonsmethod( &
        dir, & ! =1,2,3
        data_dir, & ! =0,1,2
        dxmin, &
        x_projection_raster, &
        dx, &
        n_raster, & ! points to solid
        u, & !INTENT(in) uimage_raster_solid_frame(dir)
        uimage_tngt_mag, & !INTENT(in) 
        wall_model_velocity, & ! INTENT(in)
        dist_probe, & ! INTENT(in)
        dist_fluid, & ! INTENT(in)
        temperature_image, & !INTENT(in) 
        temperature_wall, & ! INTENT(in)      
        temperature_wall_max, & ! INTENT(in)
        viscosity_molecular, & ! INTENT(in)      
        viscosity_eddy_wall, & ! INTENT(in)      
        y, & !INTENT(in) distance from image to wall
        ughost_tngt, & ! INTENT(out)
        im_fluid, &  ! INTENT(in)
        critical_length) ! INTENT(in) used for sanity check
      use probcommon_module
      implicit none
      integer, INTENT(in) :: dir ! 1,2,3
      integer, INTENT(in) :: data_dir ! 0,1,2
      real(amrex_real), INTENT(in) :: dxmin
      real(amrex_real), INTENT(in), pointer :: x_projection_raster(:)
      real(amrex_real), INTENT(in), pointer :: dx(:)
      real(amrex_real), INTENT(in), pointer :: n_raster(:) ! points to solid
      integer, INTENT(in) :: im_fluid
      real(amrex_real), INTENT(in) :: u !uimage_raster_solid_frame(dir)
      real(amrex_real), INTENT(in) :: uimage_tngt_mag
      real(amrex_real), INTENT(in) :: wall_model_velocity
      real(amrex_real), INTENT(in) :: dist_probe
      real(amrex_real), INTENT(in) :: dist_fluid
      real(amrex_real), INTENT(in) :: temperature_image
      real(amrex_real), INTENT(in) :: temperature_wall
      real(amrex_real), INTENT(in) :: temperature_wall_max
      real(amrex_real), INTENT(in) :: viscosity_molecular
      real(amrex_real), INTENT(in) :: viscosity_eddy_wall
      real(amrex_real), INTENT(in) :: y !delta_r
      real(amrex_real), INTENT(in) :: critical_length
      real(amrex_real), INTENT(out) :: ughost_tngt  ! dir direction

      real(amrex_real) :: u_tau,tau_w,x_n,x_np1 !x_n, x_(n+1) --> u_tau
      real(amrex_real) :: iter_diff
      real(amrex_real) :: f, fprime
      integer, parameter :: iter_max=1000
      integer :: iter
       !initialize wall function parameters here
      real(amrex_real), parameter :: K=0.41d0, B=5.5d0
      real(amrex_real) :: rho_w !wall density
      real(amrex_real) :: mu_w  !mu_w: wall molecular viscosity
      real(amrex_real) :: predict_deriv_utan
      real(amrex_real) :: max_deriv_utan
      real(amrex_real) :: ughost_tngt_mag
      real(amrex_real) :: thermal_conductivity
      real(amrex_real) :: thermal_diffusivity
      real(amrex_real) :: Cp
      real(amrex_real) :: u_abs

      if ((im_fluid.lt.1).or.(im_fluid.gt.num_materials)) then
       print *,"im_fluid invalid in wallfunc_newtonsmethod"
       stop
      endif

      mu_w=fort_viscconst(im_fluid) 
      rho_w=fort_denconst(im_fluid)

      if (rho_w.gt.zero) then
       ! do nothing
      else
       print *,"rho_w invalid"
       stop
      endif

      if (y.gt.zero) then
       ! do nothing
      else
       print *,"y should be positive"
       stop
      endif
      if (dxmin.gt.zero) then
       ! do nothing
      else
       print *,"dxmin should be positive"
       stop
      endif

      thermal_conductivity=fort_heatviscconst(im_fluid)
      Cp=fort_stiffCP(im_fluid)

      if (Cp.gt.zero) then
       ! do nothing
      else
       print *,"Cp invalid"
       stop
      endif

      thermal_diffusivity=thermal_conductivity/(rho_w*Cp)

      if (mu_w.eq.viscosity_molecular) then
       ! do nothing
      else
       print *,"mu_w.eq.viscosity_molecular == false"
       stop
      endif
      if (viscosity_eddy_wall.eq.fort_viscconst_eddy_wall(im_fluid)) then
       ! do nothing
      else
       print *,"viscosity_eddy_wall.eq.fort_viscconst_eddy_wall(im_fluid)"
       print *,"evaluates to false."
       stop
      endif

      if (critical_length.lt.dxmin) then

       u_abs=abs(u)
       x_n=u_abs !initial guess for u_tau
       x_np1=u_abs
       iter_diff=one
       iter=0
       do while ((iter_diff.gt.VOFTOL).and.(iter.lt.iter_max))
        f = wallfunc(x_n,u_abs,y,K,B,rho_w,mu_w)
        fprime = wallfuncderiv(x_n,u_abs,y,K,B,rho_w,mu_w)
   
        if (abs(fprime).le.EPS15) then 
         print *, "divide by zero error: wallfunc_newtonsmethod" !avoid /0
         stop
        else
         x_np1 = x_n-f/fprime
        endif
        iter_diff=abs(x_np1-x_n)
        x_n=x_np1

        iter = iter+1
        if(iter.ge.iter_max)then
         print *, "wallfunc_newtonsmethod: no convergence"
         print *, "u (uimage_raster_solid_frame(dir)) = ",u
         print *, "y (delta_r) = ",y
         print *, "im_fluid= ",im_fluid
         print *, "mu_w= ",mu_w
         print *, "rho_w= ",rho_w
         print *, "f= ",f
         print *, "fprime= ",fprime
         print *, "x_np1= ",x_np1
         print *, "x_n= ",x_n
         print *, "iter=",iter
         print *, "iter_diff = ",iter_diff
         stop
        endif
       enddo ! while (iter_diff>VOFTOL .and. iter<iter_max)
    
       u_tau = x_np1
       tau_w = rho_w*(u_tau**2)

        ! MKS units of viscosity: Pascal * seconds= 
        !  (kg m/s^2)*(1/m^2)*s=kg /(m s)
        ! MKS units of Pascal: N/m^2
        ! MKS units of tau: pascal
        ! This will not be the velocity at the ghost point, it will be
        ! a velocity at the projection (wall) point.
       ughost_tngt_mag=u_abs-tau_w*y/ &
          (viscosity_molecular+viscosity_eddy_wall)

       predict_deriv_utan=abs(ughost_tngt_mag-u_abs)/y
       max_deriv_utan=four*uimage_tngt_mag/critical_length
       if (predict_deriv_utan.lt.max_deriv_utan) then
        ! do nothing
       else
        print *,"predict_deriv_utan or max_deriv_utan invalid"
        print *,"predict_deriv_utan= ",predict_deriv_utan
        print *,"max_deriv_utan= ",max_deriv_utan
        print *,"ughost_tngt_mag=",ughost_tngt_mag
        print *,"u=",u
        print *,"uimage_tngt_mag=",uimage_tngt_mag
        print *,"critical_length= ",critical_length
        print *,"y= ",y
        stop
       endif

      else if (critical_length.ge.dxmin) then

       ughost_tngt_mag=zero
       u_tau=zero
       tau_w=u_abs*(viscosity_molecular+viscosity_eddy_wall)/y

      else
       print *,"critical_length is NaN"
       stop
      endif

      ughost_tngt=ughost_tngt_mag
      if (u.lt.zero) then
       ughost_tngt=-ughost_tngt
      else if (u.ge.zero) then
       ! do nothing
      else
       print *,"u is NaN"
       stop
      endif

      if (1.eq.0) then
       print *, "u (uimage_raster_solid_frame(dir)) = ",u
       print *, "y (delta_r) = ",y
       print *, "im_fluid= ",im_fluid
       print *, "mu_w= ",mu_w
       print *, "rho_w= ",rho_w
       print *, "u_tau= ",u_tau
       print *, "ughost_tngt= ",ughost_tngt
       print *, "tau_w= ",tau_w
       print *, "iter=",iter
       print *, "iter_diff = ",iter_diff
      endif

      end subroutine wallfunc_newtonsmethod



      subroutine wallfunc_general( &
        dir, & ! =1,2,3
        data_dir, & ! =0,1,2
        dxmin, &
        x_projection_raster, &
        dx, &
        n_raster, & ! points to solid
        u, & !INTENT(in) uimage_raster_solid_frame(dir)
        uimage_tngt_mag, & !INTENT(in) 
        wall_model_velocity, & ! INTENT(in)
        dist_probe, & ! INTENT(in)
        dist_fluid, & ! INTENT(in)
        temperature_image, & !INTENT(in) 
        temperature_wall, & ! INTENT(in)      
        temperature_wall_max, & ! INTENT(in)      
        viscosity_molecular, & ! INTENT(in)      
        viscosity_eddy_wall, & ! INTENT(in)      
        y, & !INTENT(in) distance from image to wall
        ughost_tngt, & ! INTENT(out)
        im_fluid, &  ! INTENT(in)
        critical_length) ! INTENT(in) used for sanity check
      use probcommon_module
      implicit none
      integer, INTENT(in) :: dir ! 1,2,3
      integer, INTENT(in) :: data_dir ! 0,1,2
      real(amrex_real), INTENT(in) :: dxmin
      real(amrex_real), INTENT(in), pointer :: x_projection_raster(:)
      real(amrex_real), INTENT(in), pointer :: dx(:)
      real(amrex_real), INTENT(in), pointer :: n_raster(:) ! points to solid
      integer, INTENT(in) :: im_fluid
      real(amrex_real), INTENT(in) :: u !uimage_raster_solid_frame(dir)
      real(amrex_real), INTENT(in) :: uimage_tngt_mag
      real(amrex_real), INTENT(in) :: wall_model_velocity
      real(amrex_real), INTENT(in) :: dist_probe
      real(amrex_real), INTENT(in) :: dist_fluid
      real(amrex_real), INTENT(in) :: temperature_image
      real(amrex_real), INTENT(in) :: temperature_wall
      real(amrex_real), INTENT(in) :: temperature_wall_max
      real(amrex_real), INTENT(in) :: viscosity_molecular
      real(amrex_real), INTENT(in) :: viscosity_eddy_wall
      real(amrex_real), INTENT(in) :: y !delta_r
      real(amrex_real), INTENT(in) :: critical_length
      real(amrex_real), INTENT(out) :: ughost_tngt  ! dir direction

      if (is_in_probtype_list().eq.1) then
       call SUB_wallfunc( &
        dir, & ! =1,2,3
        data_dir, & ! =0,1,2
        dxmin, &
        x_projection_raster, &
        dx, &
        n_raster, & ! points to solid
        u, & !INTENT(in) uimage_raster_solid_frame(dir)
        uimage_tngt_mag, & !INTENT(in)
        wall_model_velocity, & ! INTENT(in)
        dist_probe, & ! INTENT(in)
        dist_fluid, & ! INTENT(in)
        temperature_image, & !INTENT(in) 
        temperature_wall, & ! INTENT(in)      
        temperature_wall_max, & ! INTENT(in)      
        viscosity_molecular, & ! INTENT(in)      
        viscosity_eddy_wall, & ! INTENT(in)      
        y, & !INTENT(in) distance from image to wall
        ughost_tngt, & ! INTENT(out)
        im_fluid, &  ! INTENT(in)
        critical_length) ! INTENT(in) used for sanity check
      else
       call wallfunc_newtonsmethod( &
        dir, & ! =1,2,3
        data_dir, & ! =0,1,2
        dxmin, &
        x_projection_raster, &
        dx, &
        n_raster, & ! points to solid
        u, & !INTENT(in) uimage_raster_solid_frame(dir)
        uimage_tngt_mag, & !INTENT(in)
        wall_model_velocity, & ! INTENT(in)
        dist_probe, & ! INTENT(in)
        dist_fluid, & ! INTENT(in)
        temperature_image, & !INTENT(in) 
        temperature_wall, & ! INTENT(in)      
        temperature_wall_max, & ! INTENT(in)      
        viscosity_molecular, & ! INTENT(in)      
        viscosity_eddy_wall, & ! INTENT(in)      
        y, & !INTENT(in) distance from image to wall
        ughost_tngt, & ! INTENT(out)
        im_fluid, &  ! INTENT(in)
        critical_length) ! INTENT(in) used for sanity check
      endif

      end subroutine wallfunc_general

      subroutine interp_from_fluid( &
       LOW, &
       x_fluid, &
       im_secondary_image, &
       thermal_interp, &
       im_fluid, &
       im_solid, &
       LS_interp)
      use probcommon_module
      implicit none
       
      type(law_of_wall_parm_type), INTENT(in) :: LOW
      integer, INTENT(in) :: im_fluid
      integer, INTENT(in) :: im_solid
      real(amrex_real), INTENT(in) :: x_fluid(SDIM)
      integer, INTENT(inout) :: im_secondary_image
      real(amrex_real), INTENT(out) :: thermal_interp(num_materials)
      real(amrex_real), INTENT(out) :: LS_interp(num_materials*(1+SDIM))

      real(amrex_real) :: xsten(-3:3,SDIM)
      real(amrex_real) :: xsten_center(-3:3,SDIM)
      integer nhalf
      integer im
      integer dir
      integer cell_index(SDIM)
      integer stencil_offset(SDIM)
      integer istenlo(3)
      integer istenhi(3)
      real(amrex_real) WT,total_WT
      integer isten,jsten,ksten
      real(amrex_real) LS_sten(num_materials*(SDIM+1))
      integer im_primary_sten
      real(amrex_real) local_temperature
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: local_data_fab_LS(D_DECL(:,:,:),:)

      local_data_fab=>LOW%state
      local_data_fab_LS=>LOW%LSCP

      nhalf=3

      call containing_cell( &
        LOW%bfact, &
        LOW%dx, &
        LOW%xlo, &
        LOW%fablo, &
        x_fluid, &
        cell_index)

      istenlo(3)=0
      istenhi(3)=0
      do dir=1,SDIM
       istenlo(dir)=cell_index(dir)-1
       istenhi(dir)=cell_index(dir)+1
      enddo ! dir=1..sdim

      total_WT=zero
      do im=1,num_materials*(1+SDIM)
       LS_interp(im)=zero
      enddo 
      do im=1,num_materials
       thermal_interp(im)=zero
      enddo

      isten=cell_index(1)
      jsten=cell_index(2)
      ksten=cell_index(SDIM)

      call gridsten_level(xsten_center,isten,jsten,ksten,LOW%level,nhalf)

      do ksten=istenlo(3),istenhi(3)
      do jsten=istenlo(2),istenhi(2)
      do isten=istenlo(1),istenhi(1)

       call gridsten_level(xsten,isten,jsten,ksten,LOW%level,nhalf)
       stencil_offset(1)=isten-cell_index(1)
       stencil_offset(2)=jsten-cell_index(2)
       if (SDIM.eq.3) then
        stencil_offset(SDIM)=ksten-cell_index(SDIM)
       endif
       call bilinear_interp_WT(xsten_center,nhalf,stencil_offset, &
        x_fluid,WT)
       if ((WT.ge.zero).and.(WT.le.one)) then
        ! do nothing
       else
        print *,"WT invalid"
        stop
       endif
 
       do im=1,num_materials*(1+SDIM)
        call safe_data(isten,jsten,ksten,im,local_data_fab_LS,LS_sten(im))
       enddo

       call get_primary_material(LS_sten,im_primary_sten)
       if (im_primary_sten.eq.im_fluid) then
        ! do nothing
       else if (im_primary_sten.eq.im_solid) then
        ! do nothing
       else if ((im_primary_sten.ge.1).and. &
                (im_primary_sten.le.num_materials)) then
        if (is_rigid(im_primary_sten).eq.1) then
         ! do nothing
        else if (is_rigid(im_primary_sten).eq.0) then
         if (abs(LS_sten(im_fluid)).le.LOW%dxmin*GNBC_RADIUS) then
          if (im_secondary_image.eq.0) then
           im_secondary_image=im_primary_sten
          endif
         else if (abs(LS_sten(im_fluid)).ge.LOW%dxmin*GNBC_RADIUS) then
          ! do nothing
         else
          print *,"LS_sten became corrupt"
          stop
         endif
        else
         print *,"is_rigid(im_primary_sten) invalid"
         stop
        endif
       else
        print *,"im_primary_sten invalid"
        stop
       endif
       do im=1,num_materials*(1+SDIM)
        LS_interp(im)=LS_interp(im)+WT*LS_sten(im)
       enddo
       do im=1,num_materials
        call safe_data(isten,jsten,ksten, &
          (im-1)*num_state_material+2, &
          local_data_fab,local_temperature)
        if (local_temperature.ge.zero) then
         thermal_interp(im)=thermal_interp(im)+WT*local_temperature
        else
         print *,"local_temperature invalid"
         stop
        endif
       enddo ! im=1..num_materials
       total_WT=total_WT+WT

      enddo ! ksten
      enddo ! jsten
      enddo ! isten

      if (total_WT.gt.zero) then

       do im=1,num_materials*(1+SDIM)
        LS_interp(im)=LS_interp(im)/total_WT
       enddo
       do im=1,num_materials
        thermal_interp(im)=thermal_interp(im)/total_WT
       enddo

      else
       print *,"total_WT invalid"
       stop
      endif

      return
      end subroutine interp_from_fluid

      function ZEYU_delta(s)
      use probcommon_module
      implicit none

      real(amrex_real) ZEYU_delta
      real(amrex_real), INTENT(in) :: s
      real(amrex_real) :: r
      
      ! integral_r=0 to r=2 d(r)dr =1 
      ! s=alpha (r/2)
      ! integral_s=0 to s= alpha d(2s/alpha) 2 ds/alpha = 1
      
      if (GNBC_RADIUS.ge.1.0d0) then
       ! do nothing
      else
       print *,"GNBC_RADIUS invalid"
       stop
      endif
      
      r = 2.0d0*s/GNBC_RADIUS
      
      if (abs(r) <= 1.0d0) then
          ZEYU_delta = (3.0d0 - 2.0d0 * abs(r) + &
            sqrt(1.0d0 + 4.0d0 * abs(r) - 4.0d0 * r * r)) / 8.0d0
      else if (1.0d0 < abs(r) .AND. abs(r) <= 2.0d0) then
          ZEYU_delta = (5.0d0 - 2.0d0 * abs(r) -  &
            sqrt(-7.0d0 + 12.0d0 * abs(r) - 4.0d0 * r * r)) / 8.0d0
      else
          ZEYU_delta = 0.0d0
      end if
      
      ZEYU_delta=ZEYU_delta * 2.0d0/GNBC_RADIUS
      
      end function ZEYU_delta

! ZEYU HUANG
!mu_l: dynamic viscocity of liquid
!mu_g: dynamic viscocity of gas
!sigma: surface tension coeffient
!thet_s: static contact angle (liquid region)
!imodel: model index, can be 1 ~ 7
!ifgnbc: if use gnbc, only works in model 1
!lambda: slip length, equal to 8.e-7 in Yamamoto2013, depends on 
!specific problems
!l_macro: parameter in gnbc, can be set as grid length
!l_micro: parameter in gnbc, can be set as 1.e-9
!dgrid: grid length = dxmin when called from getGhostVel
!d_closest: closest distance to the contact line
!thet_d_apparent: dynamic contact angle from simulation (input in gnbc)
!(liquid region)
!u_cl: velocity of contact line (input in dynamic contact angle models)
! (unused for GNBC)
! u_cl is positive if the contact line is advancing into the gas.
!u_slip: slip velocity of wall (output in gnbc)
! u_slip>0 if CL advancing into the vapor region.
!thet_d: dynamic contact angle (output in dynamic contact angle models)
!(liquid region) (unused, GNBC)
!For the test results of different dynamic contact angle models, model 4 and 
!model 7 have large difference between other models.
!
! called from GLOBALUTIL.F90 and LEVELSET_3D.F90
subroutine dynamic_contact_angle(mu_l, mu_g, sigma, &
   thet_s, &
   imodel, ifgnbc, lambda, &
   l_macro, l_micro, &
   dgrid, d_closest, thet_d_apparent, &
   u_cl, &
   u_slip, &
   thet_d)
use probcommon_module
implicit none

integer imodel, ifgnbc
real(amrex_real) mu_l, mu_g, sigma, thet_s, dgrid, d_closest
real(amrex_real) thet_d_apparent, u_cl, u_slip, thet_d
real(amrex_real) lambda, l_macro, l_micro !parameter of gnbc
integer iter
real(amrex_real) Ca
real(amrex_real) sign_Ca
real(amrex_real) thet_d_micro, beta, chi !parameters of model1
real(amrex_real) thet_d_micro_old, Ca_old, a, b, c
real(amrex_real) f1, f2, Ja11, Ja12, Ja21, Ja22
real(amrex_real) b1, b2, u11, u12, u22, l21, y1, y2
real(amrex_real) a1, a2, a3, a4, u, u0, thet_d_old !parameters of model3
real(amrex_real) temp,temp2
real(amrex_real) fHI, fHI_old !parameters of model5
real(amrex_real) sigma_0, v_0 !parameters of model7
integer diag_output

diag_output=1

if (ifgnbc.eq.0) then
        if (diag_output.eq.1) then
         print *, "Implement Dynamic Contact Angle ..."
        endif
else if (ifgnbc.eq.1) then
        if (diag_output.eq.1) then
         print *, "Implement Generalized Navier Boundary Condition ..."
        endif
else
        print *,"ifgnbc invalid"
        stop
end if

 ! sanity check
if (sigma.gt.0.0d0) then
 Ca = mu_l * u_cl / sigma
 if (Ca.ge.0.0d0) then !advancing side of the liquid drop.
  sign_Ca=1.0d0
 else if (Ca.lt.0.0d0) then !receding side of the liquid drop.
  sign_Ca=-1.0d0
 else
  print *,"Ca corrupt"
  stop
 endif
else
 print *,"sigma invalid: ",sigma
 stop
endif

select case (imodel)
    case (1) !GNBC
        if (diag_output.eq.1) then
         print *, "Implement model 1 ..."
        endif
       if (abs(l_macro) < EPS9) then
           l_macro = dgrid
       end if

        if (ifgnbc.eq.0) then !If don't implement GNBC.
            print *,"imodel and ifgnbc mismatch"
            stop
        else if (ifgnbc.eq.1) then !Implement GNBC
            beta = mu_l / lambda
            chi = (mu_l + mu_g) / (2.0d0 * beta * dgrid)
            a = 9.0d0 * log(l_macro / l_micro)
            b = thet_d_apparent**3.0d0
            c = chi * cos(thet_s)
            
            iter = 0
            thet_d_micro = 0.0d0
            Ca = 0.0d0
            thet_d_micro_old = thet_d_apparent
            Ca_old = chi * (cos(thet_s) - cos(thet_d_apparent))
            do iter = 0, 1000
                f1 = thet_d_micro_old**3.0d0 + a * Ca_old - b
                f2 = Ca_old + chi * cos(thet_d_micro_old) - c
                Ja11 = 3.0d0 * thet_d_micro_old**2.0d0
                Ja12 = a
                Ja21 = - chi * sin(thet_d_micro_old)
                Ja22 = 1.0d0
                b1 = Ja11 * thet_d_micro_old + Ja12 * Ca_old - f1
                b2 = Ja21 * thet_d_micro_old + Ja22 * Ca_old - f2
                u11 = Ja11
                u12 = Ja12
                l21 = Ja21 / (u11 + EPS20)
                u22 = Ja22 - l21 * u12
                y1 = b1
                y2 = b2 - l21 * y1
                Ca = y2 / (u22 + EPS20)
                thet_d_micro = (y1 - u12 * Ca) / (u11 + EPS20)
                if (abs((thet_d_micro - thet_d_micro_old)/ &
                        (thet_d_micro_old+EPS20)) < EPS4 &
                   .AND. abs((Ca - Ca_old)/(Ca_old+EPS20)) < EPS4) then
                    exit
                end if
                thet_d_micro_old = thet_d_micro
                Ca_old = Ca
            end do
            print *, "Calculating Ca and thet_d_micro..."
            print *, "number of iteration is: ", iter
             ! beta = mu_l / lambda  e.g. lambda=8.0D-7
             ! integral_{-2 dx}^{2 dx} u_slip = u_CL_JIANG * 4 dx ?
            u_slip = 1.0d0 / (beta * dgrid + EPS20) *  &
                  ZEYU_delta(d_closest/dgrid) * sigma *  &
                  (cos(thet_s) - cos(thet_d_micro))
            thet_d = thet_d_apparent
        else
            print *,"ifgnbc invalid"
            stop
        end if

    case (8) !Cox1986
        if (diag_output.eq.1) then
         print *, "Implement model 1 ..."
        endif
       if (abs(l_macro) < 1.e-9) then
           l_macro = dgrid
       end if

        if (ifgnbc.eq.0) then !If don't implement GNBC.
             !Ca>0 if advancing => thet_d>thet_s
            thet_d = (thet_s**3 + 9. * Ca * log(l_macro / l_micro))**(1./3.)
            u_slip = 0.0d0
        else if (ifgnbc.eq.1) then !Implement GNBC
            print *,"imodel and ifgnbc mismatch"
            stop
        else
            print *,"ifgnbc invalid"
            stop
        end if

    case (2) !Jiang1970
        if (diag_output.eq.1) then
         print *, "Implement model 2, Ca= ...",Ca
        endif
        if (Ca.eq.0.0d0) then
         thet_d=thet_s
        else if (Ca.ne.0.0d0) then
         temp2=abs(Ca)**0.702d0
         temp=tanh(4.96*temp2)
         temp=temp*(1.0d0+cos(thet_s)) 
           ! if advancing side of the liquid drop, then Ca>0 =>
           ! cos(theta_dynamic)<cos(theta_static) =>
           ! theta_dynamic > theta_static    0<=theta<=pi
         temp=cos(thet_s)-(Ca/abs(Ca))*temp
         if (temp.gt.1.0d0) temp = 1.0d0 
         if (temp.lt.-1.0d0) temp = -1.0d0 
         thet_d=acos(temp)
        else
         print *,"Ca invalid"
         stop
        endif
        u_slip = 0.0d0
        if (diag_output.eq.1) then
         print *, "End Implement model 2, thet_d= ...",thet_d
        endif 

    case (3) !Shikmurzaev2008
        if (diag_output.eq.1) then
         print *, "Implement model 3 ..."
        endif
        a2 = 0.54
        a3 = 12.5
        a4 = 0.07
        a1 = 1. + (1. - a2) * (cos(thet_s) - a4)
        u = a3 * Ca
        thet_d_old = thet_d_apparent
        u0 = 0.
        do iter = 0, 1000
            u0 = (sin(thet_d_old - thet_d_old * cos(thet_d_old))) / (sin(thet_d_old) * cos(thet_d_old) - thet_d_old)
            thet_d = acos(cos(thet_s) - 2. * u * (a1 + a2 * u0) / (1. - a2) / (sqrt(a1 + u * u) + u))
            if (abs((thet_d - thet_d_old)/(thet_d_old + 1.e-20)) < 1.e-4) then
                exit
            else
               thet_d_old = thet_d
            end if
        end do
        if (diag_output.eq.1) then
         print *, "Calculating thet_d..."
         print *, "number of iteration is: ", iter
        endif
        u_slip = 0.0d0

    case (4) !Kalliadasis1994-'abs(tan(the_d))=...', so how to determine the_d<90 or thet_d>90? It seems have problems...
        print *,"Kalliadasis 1994 fails if Ca<=0, also no theta_s input!"
        stop

        if (diag_output.eq.1) then
         print *, "Implement model 4 ..."
        endif
        thet_d = atan(7.48 * Ca**(1./3.) - 3.28 * 1.e-8**0.04 * Ca**0.293)
        u_slip = 0.0d0
        
    case (5) !Kistler1993
        if (diag_output.eq.1) then
         print *, "Implement model 5 ..."
        endif 
        temp = (tanh((1. - cos(thet_s)) / 2.) / 5.16)**(1./0.706)
        fHI_old = temp / (1. - 1.31 * temp)
        do iter = 0, 1000
            fHI = temp * (1. + 1.31 * fHI_old**0.99)
            if (abs((fHI - fHI_old)/(fHI_old + 1.e-20)) < 1.e-4) then
                exit
            else
                fHI_old = fHI
            end if
        end do
        if (diag_output.eq.1) then
         print *, "Calculating fHI..."
         print *, "number of iteration is: ", iter
        endif
        thet_d = acos(1. - 2. * tanh(5.16 * ((Ca + fHI)/(1. + 1.31 * (Ca + fHI)**0.99))**0.706))
        u_slip = 0.0d0

    case (6) !modified Bracke1989

           ! Ca<0 then the model is not validated.
        if (diag_output.eq.1) then
         print *, "Implement model 6 ..."
        endif 
        thet_d = acos(cos(thet_s) -  &
           sign_Ca * 2.0d0 * (1.0d0 + cos(thet_s)) * abs(Ca)**0.5)
        u_slip = 0.0d0

    case (7) !Blake2006, Popescu2008
        print *,"Popescu 2008 model has dimensional parameters; disallowed!"
        stop

        if (diag_output.eq.1) then
         print *, "Implement model 7 ..."
        endif
         ! sigma_0=1.7e-2=0.017  n/m  for water
         ! sigma_0=1.0e-2=0.01   n/m  for oil
         ! v_0=5.0e-3 m/s for water
         ! v_0=1.35e-5 m/s for oil
        sigma_0 = 1.7e-2 !sigma_0 = 1.7e-2 for water and 1.e-2 for oil in Popescu2008
        v_0 = 5.e-3 !v_0 = 5e-3 for water and 1.35e-5 for oil in Popescu2008
        thet_d = acos(cos(thet_s) - sigma_0 / sigma * asinh(u_cl / v_0))
        u_slip = 0.0d0

    case default
        print *, "unknown model index"
        stop
end select


end subroutine dynamic_contact_angle



       ! This routine transforms all velocities to a frame of reference with
       ! respect to the solid.
       ! getGhostVel is called from "fort_wallfunction" (GODUNOV_3D.F90)
      subroutine getGhostVel( &
       LOW, &
       law_of_the_wall, &
       iSOLID,jSOLID,kSOLID, & ! index for the solid cell
       iFLUID,jFLUID,kFLUID, & ! index for the fluid cell
       i_probe,j_probe,k_probe, & ! index for the fluid probe cell
       side_solid, &  ! =0 if solid on the left
       side_image, &  ! =0 if image on the left
       data_dir, &  ! normal dir=0..sdim-1
       uimage_raster, & ! in (at the probe)
       wall_model_velocity, & ! INTENT(in)
       dist_probe, & ! INTENT(in)
       dist_fluid, & ! INTENT(in)
       temperature_image, & ! in (at the probe)
       temperature_wall, & ! in
       temperature_wall_max, & ! in
       usolid_law_of_wall, & ! out
       angle_ACT, & ! aka angle_ACT_cell "out"
       im_fluid, &
       im_solid)
       use probcommon_module
       implicit none
       
       type(law_of_wall_parm_type), INTENT(in) :: LOW
       integer, INTENT(in) :: law_of_the_wall
       integer, INTENT(in) :: data_dir ! normal dir=0..sdim-1
       integer, INTENT(in) :: im_fluid
       integer, INTENT(in) :: im_solid
       integer, INTENT(in) :: side_solid
       integer, INTENT(in) :: side_image
       integer, INTENT(in) :: iSOLID,jSOLID,kSOLID
       integer, INTENT(in) :: iFLUID,jFLUID,kFLUID
       integer, INTENT(in) :: i_probe,j_probe,k_probe
       real(amrex_real), dimension(SDIM), INTENT(in) :: uimage_raster
       real(amrex_real), INTENT(in) :: wall_model_velocity
       real(amrex_real), INTENT(in) :: dist_probe
       real(amrex_real), INTENT(in) :: dist_fluid
       real(amrex_real), INTENT(in) :: temperature_image
       real(amrex_real), INTENT(in) :: temperature_wall
       real(amrex_real), INTENT(in) :: temperature_wall_max
       real(amrex_real), dimension(SDIM), INTENT(out) :: usolid_law_of_wall
       real(amrex_real), INTENT(out) :: angle_ACT
      
       !D_DECL is defined in SPACE.H in the BoxLib/Src/C_BaseLib

       real(amrex_real), dimension(SDIM) :: uimage_raster_solid_frame

        ! uimage_tngt_mag used for estimating the boundary layer thickness
       real(amrex_real) :: uimage_tngt_mag 

       real(amrex_real) :: ughost_tngt(SDIM)
       real(amrex_real) :: viscosity_molecular, viscosity_eddy_wall
       real(amrex_real) :: density_fluid
       integer :: dir
       real(amrex_real) :: nrm_sanity
       real(amrex_real) :: nrm_sanity_crossing
       real(amrex_real) :: critical_length
       integer :: im
       integer :: im_fluid_crossing
       integer :: im_primary_image
       integer :: im_secondary_image
       integer :: im_fluid1,im_fluid2
       integer :: iten
       integer :: iten_13,iten_23
       real(amrex_real) :: cos_angle,sin_angle
       integer :: near_contact_line
       real(amrex_real) :: mag
       real(amrex_real) :: sinthetaACT
       real(amrex_real) :: costhetaACT
       real(amrex_real) :: dist_to_CL
       real(amrex_real) :: nf_dot_ns
       real(amrex_real) :: nf_crossing_dot_ns
       real(amrex_real) :: nf_dot_nCL_perp
       real(amrex_real) :: nf_crossing_dot_nCL_perp
       real(amrex_real), dimension(3) :: nCL
       real(amrex_real), dimension(3) :: nCL_crossing
       real(amrex_real), dimension(3) :: nCL_raster
       real(amrex_real), dimension(3) :: nCL_perp
       real(amrex_real), dimension(3) :: nCL_perp_crossing
       real(amrex_real), dimension(3) :: nCL_perp2
       real(amrex_real), dimension(3) :: nCL_perp2_crossing
       real(amrex_real), dimension(3) :: nf_prj
       real(amrex_real), dimension(3) :: nf_prj_crossing
       real(amrex_real) :: ZEYU_mu_l, ZEYU_mu_g, ZEYU_sigma
       real(amrex_real) :: ZEYU_thet_s,ZEYU_lambda,ZEYU_l_macro, ZEYU_l_micro
       real(amrex_real) :: ZEYU_dgrid, ZEYU_d_closest, ZEYU_thet_d_apparent
       real(amrex_real) :: ZEYU_u_cl, ZEYU_u_slip, ZEYU_thet_d
       real(amrex_real) :: angle_im1
       integer :: ZEYU_imodel
       integer :: ZEYU_ifgnbc
       integer :: im_vapor,im_liquid
       real(amrex_real) :: delta_r_raster
       real(amrex_real) :: nCL_dot_n_raster
       integer :: nrad
       integer nhalf
       real(amrex_real) :: xsten_probe(-3:3,SDIM)
       real(amrex_real) :: xstenFLUID(-3:3,SDIM)
       real(amrex_real) :: xstenSOLID(-3:3,SDIM)
       real(amrex_real) :: thermal_interp(num_materials)
       real(amrex_real) :: LSPLUS_interp(num_materials*(1+SDIM))
       real(amrex_real) :: LSMINUS_interp(num_materials*(1+SDIM))
       real(amrex_real) :: LSTRIPLE_interp(num_materials*(1+SDIM))
       real(amrex_real) :: LS_probe(num_materials*(1+SDIM))
       real(amrex_real) :: LS_fluid(num_materials*(1+SDIM))
       real(amrex_real) :: LS_solid(num_materials*(1+SDIM))
       real(amrex_real) :: LS_crossing(num_materials*(1+SDIM))
       real(amrex_real) :: LS_triple(num_materials*(1+SDIM))
       real(amrex_real) :: nrm_solid(3)
       real(amrex_real) :: nrm_fluid(3)
       real(amrex_real) :: nrm_fluid_crossing(3)
       real(amrex_real), allocatable, dimension(:) :: user_tension
       real(amrex_real) :: cross_denom
       real(amrex_real) :: cross_factor
       real(amrex_real) :: cross_factorMINUS
       real(amrex_real) :: cross_factorPLUS
       integer :: cross_factor_flag
       real(amrex_real) :: xcrossing(SDIM)
       real(amrex_real) :: xprobeMINUS_crossing(SDIM)
       real(amrex_real) :: xprobePLUS_crossing(SDIM)
       real(amrex_real) :: xprobe_triple(SDIM)
       real(amrex_real) :: xtriple(SDIM)
       integer :: debug_slip_velocity_enforcement

       if (1.eq.0) then
        print *,"in getGhostVel"
        print *,"ijkSOLID ",iSOLID,jSOLID,kSOLID
        print *,"ijkFLUID ",iFLUID,jFLUID,kFLUID
       endif

       debug_slip_velocity_enforcement=0
    
       nhalf=3 
       allocate(user_tension(num_interfaces))
        
       if ((data_dir.ge.0).and.(data_dir.lt.SDIM)) then
        ! do nothing
       else
        print *,"data_dir invalid"
        stop
       endif
       if ((im_fluid.lt.1).or.(im_fluid.gt.num_materials)) then
        print *,"im_fluid invalid in getGhostVel"
        stop
       endif
       if ((im_solid.lt.1).or.(im_solid.gt.num_materials)) then
        print *,"im_solid invalid in getGhostVel"
        stop
       endif
       if (is_rigid(im_solid).eq.1) then
        ! do nothing
       else
        print *,"is_rigid(im_solid) invalid"
        stop
       endif
       if (LOW%dt.gt.zero) then
        ! do nothing
       else
        print *,"dt invalid"
        stop
       endif 
       if (LOW%time.ge.zero) then
        ! do nothing
       else
        print *,"time invalid"
        stop
       endif 
       if (LOW%visc_coef.ge.zero) then
        ! do nothing
       else
        print *,"visc_coef invalid: ",LOW%visc_coef
        stop
       endif 
       if ((law_of_the_wall.eq.1).or. &
           (law_of_the_wall.eq.2)) then
        ! do nothing
       else
        print *,"law_of_the_wall invalid"
        stop
       endif
       if (abs(wall_model_velocity).le.1.0D+20) then
        ! do nothing
       else
        print *,"wall_model_velocity became corrupt"
        stop
       endif
       if (LOW%dxmin.gt.zero) then
        ! do nothing
       else
        print *,"LOW%dxmin invalid"
        stop
       endif

       nrad=3

        !real(amrex_real), dimension(SDIM), INTENT(in) :: uimage_raster
        !type(law_of_wall_parm_type), INTENT(in) :: LOW
        !LOW%usolid_raster 
       do dir=1,SDIM
        usolid_law_of_wall(dir)=zero
        ughost_tngt(dir)=zero
       enddo

       angle_ACT=zero

       call gridsten_level(xstenFLUID,iFLUID,jFLUID,kFLUID,LOW%level,nhalf)
       call gridsten_level(xsten_probe,i_probe,j_probe,k_probe,LOW%level,nhalf)
       do im=1,num_materials*(1+SDIM)
        LS_probe(im)=LOW%LSCP(D_DECL(i_probe,j_probe,k_probe),im)
        LS_fluid(im)=LOW%LSCP(D_DECL(iFLUID,jFLUID,kFLUID),im)
       enddo
       call gridsten_level(xstenSOLID,iSOLID,jSOLID,kSOLID,LOW%level,nhalf)
       do im=1,num_materials*(1+SDIM)
        LS_solid(im)=LOW%LSCP(D_DECL(iSOLID,jSOLID,kSOLID),im)
       enddo

       near_contact_line=0
       im_primary_image=im_fluid

       if (LS_solid(im_solid).ge.zero) then
        if (LS_fluid(im_fluid).ge.zero) then ! im_fluid dominates the fluids
         if (LS_fluid(im_solid).le.zero) then
          cross_denom=LS_solid(im_solid)-LS_fluid(im_solid)
          if (cross_denom.gt.zero) then
           cross_factor=LS_solid(im_solid)/cross_denom
          else if (cross_denom.eq.zero) then
           cross_factor=half
          else
           print *,"cross_factor invalid"
           stop
          endif
          if ((cross_factor.ge.zero).and. &
              (cross_factor.le.one)) then
              ! xcrossing is where the solid interface passes inbetween
              ! the ghost (solid) point and the fluid (image) point.
           do dir=1,SDIM
            xcrossing(dir)=cross_factor*xstenFLUID(0,dir)+ &
                    (one-cross_factor)*xstenSOLID(0,dir)
           enddo

           do im=1,num_materials*(1+SDIM)
            LS_crossing(im)=cross_factor*LS_fluid(im)+ &
                     (one-cross_factor)*LS_solid(im)
           enddo ! im=1..num_materials*(1+SDIM)

           im_fluid_crossing=-1
           do im=1,num_materials
            if (is_rigid(im).eq.1) then
             ! do nothing
            else if (is_rigid(im).eq.0) then
             if (im_fluid_crossing.eq.-1) then
              im_fluid_crossing=im
             else if ((im_fluid_crossing.ge.1).and. &
                      (im_fluid_crossing.le.num_materials)) then
              if (LS_crossing(im_fluid_crossing).le.LS_crossing(im)) then
               im_fluid_crossing=im
              endif
             else
              print *,"im_fluid_crossing invalid"
              stop
             endif
            else 
             print *,"is_rigid invalid GLOBALUTIL.F90"
             stop
            endif
           enddo ! im=1..num_materials

           call normalize_LS_normals(LS_crossing)

           if ((im_fluid_crossing.ge.1).and. &
               (im_fluid_crossing.le.num_materials)) then
            ! do nothing
           else
            print *,"im_fluid_crossing invalid"
            stop
           endif

           nf_dot_ns=zero
           nf_crossing_dot_ns=zero
           nrm_fluid(3)=zero
           nrm_fluid_crossing(3)=zero
           nrm_solid(3)=zero
           do dir=1,SDIM
             ! points into im_fluid material
            nrm_fluid(dir)= &
                LS_crossing(num_materials+(im_fluid-1)*SDIM+dir)
            nrm_fluid_crossing(dir)= &
                LS_crossing(num_materials+(im_fluid_crossing-1)*SDIM+dir)
             ! points into im_solid material
            nrm_solid(dir)=LS_crossing(num_materials+(im_solid-1)*SDIM+dir)
            nf_dot_ns=nf_dot_ns+ &
                  nrm_fluid(dir)*nrm_solid(dir)
            nf_crossing_dot_ns=nf_crossing_dot_ns+ &
                  nrm_fluid_crossing(dir)*nrm_solid(dir)
           enddo ! dir=1..sdim

              !    \
              !nCL  \
              !  <-- \
              ! ------O----
              ! 
              ! nCL is tangential to the solid/fluid interface
              ! since nrm_solid has been projected away from nrm_fluid
              ! nCL_perp points out of the screen.
           nrm_sanity=zero
           nrm_sanity_crossing=zero
           do dir=1,3
            nCL(dir)=nrm_fluid(dir)- &
                    nf_dot_ns*nrm_solid(dir)
            nrm_sanity=nrm_sanity+nCL(dir)**2
            nCL_crossing(dir)=nrm_fluid_crossing(dir)- &
                    nf_crossing_dot_ns*nrm_solid(dir)
            nrm_sanity_crossing=nrm_sanity_crossing+nCL_crossing(dir)**2
           enddo 
           nrm_sanity=sqrt(nrm_sanity)
           nrm_sanity_crossing=sqrt(nrm_sanity_crossing)

           if (nrm_sanity.gt.zero) then
            do dir=1,3
             nCL(dir)=nCL(dir)/nrm_sanity
            enddo
           else if (nrm_sanity.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity invalid"
            stop
           endif

           if (nrm_sanity_crossing.gt.zero) then
            do dir=1,3
             nCL_crossing(dir)=nCL_crossing(dir)/nrm_sanity_crossing
            enddo
           else if (nrm_sanity_crossing.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity_crossing invalid"
            stop
           endif

           ! nCL_perp is tangent to the contact line in the substrate plane.
           ! nCL is normal to the contact line in the substrate plane
           ! nrm_solid is normal to the substrate
           call crossprod(nCL,nrm_solid,nCL_perp)
           call crossprod(nCL_crossing,nrm_solid,nCL_perp_crossing)

           nrm_sanity=zero
           nrm_sanity_crossing=zero
           do dir=1,3
            nrm_sanity=nrm_sanity+nCL_perp(dir)**2
            nrm_sanity_crossing=nrm_sanity_crossing+nCL_perp_crossing(dir)**2
           enddo 
           nrm_sanity=sqrt(nrm_sanity)
           nrm_sanity_crossing=sqrt(nrm_sanity_crossing)

           if (nrm_sanity.gt.zero) then
            do dir=1,3
             nCL_perp(dir)=nCL_perp(dir)/nrm_sanity
            enddo
           else if (nrm_sanity.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity invalid"
            stop
           endif

           if (nrm_sanity_crossing.gt.zero) then
            do dir=1,3
             nCL_perp_crossing(dir)=nCL_perp_crossing(dir)/nrm_sanity_crossing
            enddo
           else if (nrm_sanity_crossing.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity_crossing invalid"
            stop
           endif

           nf_dot_nCL_perp=zero
           nf_crossing_dot_nCL_perp=zero
           do dir=1,3

            nf_dot_nCL_perp= &
               nf_dot_nCL_perp+ &
               nrm_fluid(dir)*nCL_perp(dir)

            nf_crossing_dot_nCL_perp= &
               nf_crossing_dot_nCL_perp+ &
               nrm_fluid_crossing(dir)*nCL_perp_crossing(dir)

           enddo
           ! nCL_perp is tangent to the contact line in the substrate plane.
           ! nrm_fluid points into im_fluid material.
           ! nf_prj is the fluid normal with the tangent contact line
           ! vector projected away; nf_prj, in a plane perpendicular
           ! to the substrate and perpendicular to the contact line,
           ! is the normal point into im_primary
           ! "im_primary_image"
           nrm_sanity=zero
           nrm_sanity_crossing=zero
           do dir=1,3
            nf_prj(dir)=nrm_fluid(dir)- &
                 nf_dot_nCL_perp*nCL_perp(dir)
            nrm_sanity=nrm_sanity+nf_prj(dir)**2

            nf_prj_crossing(dir)=nrm_fluid_crossing(dir)- &
                 nf_crossing_dot_nCL_perp*nCL_perp_crossing(dir)
            nrm_sanity_crossing=nrm_sanity_crossing+nf_prj_crossing(dir)**2
           enddo

           nrm_sanity=sqrt(nrm_sanity)
           nrm_sanity_crossing=sqrt(nrm_sanity_crossing)

           if (nrm_sanity.gt.zero) then
            do dir=1,3
             nf_prj(dir)=nf_prj(dir)/nrm_sanity
            enddo
           else if (nrm_sanity.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity invalid"
            stop
           endif

           if (nrm_sanity_crossing.gt.zero) then
            do dir=1,3
             nf_prj_crossing(dir)=nf_prj_crossing(dir)/nrm_sanity_crossing
            enddo
           else if (nrm_sanity_crossing.eq.zero) then
            ! do nothing
           else
            print *,"nrm_sanity_crossing invalid"
            stop
           endif


           ! nCL_perp2 will also be tangent to the contact line.
           call crossprod(nrm_solid,nf_prj,nCL_perp2)
           call crossprod(nrm_solid,nf_prj_crossing,nCL_perp2_crossing)

            ! now we update nrm_fluid by finding the triple point
            ! and measuring nrm_fluid just above the solid.
            ! LS_crossing is the levelset function at the point in 
            ! between the ghost point and the image point at which
            ! LS_crossing(im_solid)=0
            ! NOTE: LS_IMAGE(im_fluid) > 0
            !      x  image point in the fluid , LS(im_fluid)>0 here
            !   --------   substrate interface
            !      x  ghost point in the substrate, LS(im_fluid) expected to
            !      be positive here too since the level set function is
            !      extrapolated into the substrate normal to the substrate.
            ! This routine is only called if LS_IMAGE(im_solid) <0 and
            ! LS_GHOST(im_solid)>0.
            ! xprobe_crossing is another point on the substrate:
            !
            !                       n_im_fluid_project
            !    xprobe_crossing     <----         x image  im_fluid
            !   ----x------------------------------x---------- crossing
            !                                      x ghost  im_fluid
            !  if LS_{im_fluid}(xprobe_crossing) *
            !     LS_{im_fluid}(xcrossing) < 0 =>
            !  a triple points exists in between.
            !  if LS_{im_fluid}(xcrossing) * LS_{im_fluid}(xghost) < 0 =>
            !  the triple point is already in the immediate vicinity, no
            !  need to search for it.
           if ((GNBC_RADIUS.ge.one).and. &
               (GNBC_RADIUS.le.three)) then
            do dir=1,SDIM
             xprobeMINUS_crossing(dir)=xcrossing(dir)- &
               GNBC_RADIUS*LOW%dxmin*nCL_crossing(dir)
             xprobePLUS_crossing(dir)=xcrossing(dir)+ &
               GNBC_RADIUS*LOW%dxmin*nCL_crossing(dir)
            enddo ! dir=1..sdim

            im_secondary_image=0

            call interp_from_fluid( &
             LOW, &
             xprobeMINUS_crossing, &
             im_secondary_image, &
             thermal_interp, &
             im_fluid, &
             im_solid, &
             LSMINUS_interp)

            call normalize_LS_normals(LSMINUS_interp)

            call interp_from_fluid( &
             LOW, &
             xprobePLUS_crossing, &
             im_secondary_image, &
             thermal_interp, &
             im_fluid, &
             im_solid, &
             LSPLUS_interp)

            call normalize_LS_normals(LSPLUS_interp)

            cross_factorMINUS=-one
            cross_factorPLUS=-one

            if (LSMINUS_interp(im_fluid_crossing)* &
                LS_crossing(im_fluid_crossing).le.zero) then
             cross_denom=LS_crossing(im_fluid_crossing)- &
                         LSMINUS_interp(im_fluid_crossing)
             if (cross_denom.ne.zero) then
              cross_factorMINUS=LS_crossing(im_fluid_crossing)/cross_denom
             else if (cross_denom.eq.zero) then
              cross_factorMINUS=half
             else
              print *,"cross_denom invalid"
              stop
             endif
            endif

            if (LSPLUS_interp(im_fluid_crossing)* &
                LS_crossing(im_fluid_crossing).le.zero) then
             cross_denom=LS_crossing(im_fluid_crossing)- &
                         LSPLUS_interp(im_fluid_crossing)
             if (cross_denom.ne.zero) then
              cross_factorPLUS=LS_crossing(im_fluid_crossing)/cross_denom
             else if (cross_denom.eq.zero) then
              cross_factorPLUS=half
             else
              print *,"cross_denom invalid"
              stop
             endif
            endif

            if ((cross_factorPLUS.eq.-one).and. &
                (cross_factorMINUS.eq.-one)) then
             cross_factor_flag=0
            else if ((cross_factorPLUS.eq.-one).and. &
                     (cross_factorMINUS.ge.zero)) then
             cross_factor_flag=-1
            else if ((cross_factorMINUS.eq.-one).and. &
                     (cross_factorPLUS.ge.zero)) then
             cross_factor_flag=1
            else if ((cross_factorMINUS.ge.zero).and. &
                     (cross_factorPLUS.ge.zero)) then
             if (cross_factorPLUS.le.cross_factorMINUS) then
              cross_factor_flag=1
             else if (cross_factorPLUS.ge.cross_factorMINUS) then
              cross_factor_flag=-1
             else
              print *,"cross_factor bust"
              stop
             endif
            else
             print *,"cross_factor bust"
             stop
            endif
           
            if (cross_factor_flag.eq.1) then 

             if ((cross_factorPLUS.ge.zero).and. &
                 (cross_factorPLUS.le.one)) then
              do dir=1,SDIM
               xtriple(dir)=cross_factorPLUS*xprobePLUS_crossing(dir)+ &
                  (one-cross_factorPLUS)*xcrossing(dir)
              enddo
              do im=1,num_materials*(1+SDIM)
               LS_triple(im)=cross_factorPLUS*LSPLUS_interp(im)+ &
                   (one-cross_factorPLUS)*LS_crossing(im)
              enddo
             else
              print *,"cross_factorPLUS bust"
              stop
             endif
             
            else if (cross_factor_flag.eq.-1) then 

             if ((cross_factorMINUS.ge.zero).and. &
                 (cross_factorMINUS.le.one)) then
              do dir=1,SDIM
               xtriple(dir)=cross_factorMINUS*xprobeMINUS_crossing(dir)+ &
                  (one-cross_factorMINUS)*xcrossing(dir)
              enddo
              do im=1,num_materials*(1+SDIM)
               LS_triple(im)=cross_factorMINUS*LSMINUS_interp(im)+ &
                   (one-cross_factorMINUS)*LS_crossing(im)
              enddo
             else
              print *,"cross_factorMINUS bust"
              stop
             endif

            else if (cross_factor_flag.eq.0) then
             ! do nothing
            else
             print *,"cross_factor_flag invalid"
             stop
            endif

            if ((cross_factor_flag.eq.1).or. &
                (cross_factor_flag.eq.-1)) then
              call normalize_LS_normals(LS_triple)

              nrm_solid(3)=zero
              do dir=1,SDIM
               nrm_solid(dir)=LS_triple(num_materials+(im_solid-1)*SDIM+dir)
               xprobe_triple(dir)=xtriple(dir)-LOW%dxmin*nrm_solid(dir)
              enddo
              
              call interp_from_fluid( &
               LOW, &
               xprobe_triple, &
               im_secondary_image, &
               thermal_interp, &
               im_fluid, &
               im_solid, &
               LSTRIPLE_interp)

              call normalize_LS_normals(LSTRIPLE_interp)
              nrm_fluid(3)=zero
              do dir=1,SDIM
               nrm_fluid(dir)=LSTRIPLE_interp(num_materials+(im_fluid-1)*SDIM+dir)
              enddo

              nf_dot_ns=zero
              do dir=1,SDIM
               nf_dot_ns=nf_dot_ns+nrm_fluid(dir)*nrm_solid(dir)
              enddo ! dir=1..sdim
 
              nrm_sanity=zero
              do dir=1,3
               nCL(dir)=nrm_fluid(dir)-nf_dot_ns*nrm_solid(dir)
               nrm_sanity=nrm_sanity+nCL(dir)**2
              enddo 
              nrm_sanity=sqrt(nrm_sanity)
              if (nrm_sanity.gt.zero) then
               do dir=1,3
                nCL(dir)=nCL(dir)/nrm_sanity
               enddo
              else if (nrm_sanity.eq.zero) then
               ! do nothing
              else
               print *,"nrm_sanity invalid"
               stop
              endif

              call crossprod(nCL,nrm_solid,nCL_perp)

              nrm_sanity=zero
              do dir=1,3
               nrm_sanity=nrm_sanity+nCL_perp(dir)**2
              enddo 
              nrm_sanity=sqrt(nrm_sanity)
              if (nrm_sanity.gt.zero) then
               do dir=1,3
                nCL_perp(dir)=nCL_perp(dir)/nrm_sanity
               enddo
              else if (nrm_sanity.eq.zero) then
               ! do nothing
              else
               print *,"nrm_sanity invalid"
               stop
              endif

              nf_dot_nCL_perp=zero
              do dir=1,3
               nf_dot_nCL_perp=nf_dot_nCL_perp+nrm_fluid(dir)*nCL_perp(dir)
              enddo
              nrm_sanity=zero
              do dir=1,3
               nf_prj(dir)=nrm_fluid(dir)-nf_dot_nCL_perp*nCL_perp(dir)
               nrm_sanity=nrm_sanity+nf_prj(dir)**2
              enddo
              nrm_sanity=sqrt(nrm_sanity)
              if (nrm_sanity.gt.zero) then
               do dir=1,3
                nf_prj(dir)=nf_prj(dir)/nrm_sanity
               enddo
              else if (nrm_sanity.eq.zero) then
               ! do nothing
              else
               print *,"nrm_sanity invalid"
               stop
              endif
              call crossprod(nrm_solid,nf_prj,nCL_perp2)

              near_contact_line=1
              if (im_secondary_image.eq.0) then
               near_contact_line=0
              else if ((im_secondary_image.ge.1).and. &
                       (im_secondary_image.le.num_materials)) then
               if (is_rigid(im_secondary_image).eq.1) then
                near_contact_line=0
               else if (is_rigid(im_secondary_image).eq.0) then
                ! do nothing
               else
                print *,"is_rigid(im_secondary_image) invalid"
                stop
               endif
              else
               print *,"im_secondary_image invalid"
               stop
              endif

              if (near_contact_line.eq.1) then
               if (im_primary_image.lt.im_secondary_image) then
                im_fluid1=im_primary_image
                im_fluid2=im_secondary_image
               else if (im_primary_image.gt.im_secondary_image) then
                im_fluid2=im_primary_image
                im_fluid1=im_secondary_image
               else
                print *,"im_primary_image or im_secondary_image invalid"
                stop
               endif
               call get_iten(im_fluid1,im_fluid2,iten)
                ! in: subroutine getGhostVel
               call get_user_tension(xprobe_triple,LOW%time, &
                 fort_tension,user_tension,thermal_interp)
               ! cos_angle and sin_angle correspond to the angle in im_fluid1
               call get_CL_iten(im_fluid1,im_fluid2,im_solid, &
                 iten_13,iten_23, &
                 user_tension,cos_angle,sin_angle)
               sinthetaACT=zero
               costhetaACT=zero
                ! because nrm_solid and nf_prj have unit magnitude,
                ! sinthetaACT is the sine of the angle between nrm_solid
                ! and nf_prj; fluid normal in plane perpendicular to 
                ! contact line and solid.
               do dir=1,3
                sinthetaACT=sinthetaACT+nCL_perp2(dir)**2
                costhetaACT=costhetaACT+nrm_solid(dir)*nf_prj(dir)
               enddo
               sinthetaACT=sqrt(sinthetaACT)

               dist_to_CL=zero
               do dir=1,SDIM
                dist_to_CL=dist_to_CL+(xtriple(dir)-xcrossing(dir))**2
               enddo
               dist_to_CL=sqrt(dist_to_CL)
                
               if ((sinthetaACT.ge.zero).and.(costhetaACT.ge.zero)) then
                angle_ACT=asin(sinthetaACT)
               else if ((sinthetaACT.ge.zero).and.(costhetaACT.le.zero)) then
                angle_ACT=Pi-asin(sinthetaACT)
               else
                print *,"sinthetaACT or costhetaACT invalid"
                stop
               endif
               if (DEBUG_DYNAMIC_CONTACT_ANGLE.eq.1) then
                print *,"xcrossing ",xcrossing(1),xcrossing(2),xcrossing(SDIM)
                print *,"xtriple ",xtriple(1),xtriple(2),xtriple(SDIM)
                print *,"nrm_solid ",nrm_solid(1),nrm_solid(2),nrm_solid(SDIM)
                print *,"nrm_fluid ",nrm_fluid(1),nrm_fluid(2),nrm_fluid(SDIM)
                print *,"im_fluid,angle_ACT(rad,deg) ",im_fluid,angle_ACT, &
                        angle_ACT*180.0d0/Pi
                print *,"dx(1),dist_to_CL ",LOW%dx(1),dist_to_CL
                print *,"im_primary_image,im_secondary_image ", &
                        im_primary_image,im_secondary_image
               endif
              else if (near_contact_line.eq.0) then
               ! do nothing
              else
               print *,"near_contact_line invalid"
               stop
              endif
            else if (cross_factor_flag.eq.0) then
              ! do nothing
            else
             print *,"cross_factor_flag invalid"
             stop
            endif
           else
            print *,"GNBC_RADIUS invalid"
            stop
           endif
          else
           print *,"cross_factor invalid"
           stop
          endif
         else
          print *,"expecting LS_fluid(im_solid)<=0"
          stop
         endif
        else
         print *,"expecting LS_fluid(im_fluid)>=0"
         stop
        endif
       else
        print *,"expecting LS_solid(im_solid)>=0"
        stop
       endif

       viscosity_molecular=fort_viscconst(im_fluid)
       viscosity_eddy_wall=fort_viscconst_eddy_wall(im_fluid)
       density_fluid=fort_denconst(im_fluid)

       if (density_fluid.gt.zero) then
        ! do nothing
       else
        print *,"density_fluid invalid"
        stop
       endif

       !x_projection is closest point on the fluid/solid interface. 
       delta_r_raster=abs(LOW%x_probe_raster(data_dir+1)- &
                          LOW%x_projection_raster(data_dir+1))
       if (delta_r_raster.gt.zero) then
        ! do nothing
       else
        print *,"delta_r_raster invalid"
        stop
       endif
      
       ! convert to solid velocity frame of reference.
       uimage_tngt_mag=zero
       do dir=1,SDIM
        uimage_raster_solid_frame(dir)= &
           uimage_raster(dir)-LOW%usolid_raster(dir)
        if (dir.ne.data_dir+1) then
         if (LOW%n_raster(dir).eq.zero) then
          ! do nothing
         else
          print *,"LOW%n_raster(dir) invalid"
          stop
         endif
         uimage_tngt_mag=uimage_tngt_mag+uimage_raster_solid_frame(dir)**2
        else if (dir.eq.data_dir+1) then
         if ((side_solid.eq.1).and.(side_image.eq.0)) then
          if (LOW%n_raster(dir).eq.one) then
           ! do nothing
          else
           print *,"LOW%n_raster(dir) invalid"
           stop
          endif
         else if ((side_solid.eq.0).and.(side_image.eq.1)) then
          if (LOW%n_raster(dir).eq.-one) then
           ! do nothing
          else
           print *,"LOW%n_raster(dir) invalid"
           stop
          endif
         else
          print *,"side_solid or side_image invalid"
          stop
         endif
        else
         print *,"dir or data_dir corrupt"
         stop
        endif
       enddo
       uimage_tngt_mag=sqrt(uimage_tngt_mag)

       ! laminar boundary layer thickness:
       ! delta=4.91 (mu x / (U rho))^(1/2)
       ! dx=4.91 (mu dx/ (U rho))^(1/2)
       ! dx=4.91^2 mu/(U rho)
       !   =(g/(cm s))/(cm/s  g/cm^3)=(g/(cm s))/(g/(cm^2 s))=cm
       if (uimage_tngt_mag.gt.zero) then
        critical_length=((4.91D0)**2)*viscosity_molecular/ &
              (density_fluid*uimage_tngt_mag)
       else if (uimage_tngt_mag.eq.zero) then
        critical_length=LOW%dxmin*1.0D+10
       else
        print *,"uimage_tngt_mag invalid"
        stop
       endif

       if ((viscosity_molecular.gt.zero).and. &
           (viscosity_eddy_wall.ge.zero)) then
        ! do nothing
       else
        print *,"viscosity_molecular.le.zero or viscosity_eddy_wall.lt.zero"
        stop
       endif

       ! From Spaldings' paper, a representative size for the linear 
       ! region is:
       ! y+ = 10.0
       ! y+ = y sqrt(tau rho)/mu_molecular
       !  or y+ is the intersection point of the linear (viscous) 
       !  and log layer profiles.

       if (law_of_the_wall.eq.1) then ! turbulence modeling here.

        do dir=1,SDIM
         if (dir.ne.data_dir+1) then

          if (viscosity_eddy_wall.gt.zero) then

           !obtain ughost_tngt(dir)
           !delta_r_raster is distance from image point to the wall.
           if (delta_r_raster.gt.zero) then
            call wallfunc_general( &
             dir, & ! = 1,2,3
             data_dir, & ! = 0,1,2
             LOW%dxmin, &
             LOW%x_projection_raster, & ! coordinate on the wall
             LOW%dx, &
             LOW%n_raster, & ! points to solid
             uimage_raster_solid_frame(dir), &
             uimage_tngt_mag, &
             wall_model_velocity, & ! INTENT(in)
             dist_probe, & ! INTENT(in)
             dist_fluid, & ! INTENT(in)
             temperature_image, &
             temperature_wall, &
             temperature_wall_max, &
             viscosity_molecular, &
             viscosity_eddy_wall, &
             delta_r_raster, &
             ughost_tngt(dir), &
             im_fluid, &
             critical_length) 

            if (1.eq.0) then
             print *,"after wallfunc_general"
             print *,"dir,ughost_tngt (projection point vel)=", &
                     dir,ughost_tngt(dir)
             print *,"delta_r_raster=",delta_r_raster
             print *,"viscosity_molecular ",viscosity_molecular
             print *,"viscosity_eddy_wall ",viscosity_eddy_wall
            endif

           else
            print *,"delta_r_raster invalid"
            stop
           endif

          else if (viscosity_eddy_wall.eq.zero) then
           ! ghost velocity lives *on* the rasterized interface.
           ughost_tngt(dir) = zero
          else
           print *,"viscosity_eddy_wall invalid"
           stop
          endif
         else if (dir.eq.data_dir+1) then
          ! do nothing
         else
          print *,"dir or data_dir bust"
          stop
         endif
        enddo !dir=1..sdim

       else if (law_of_the_wall.eq.2) then ! GNBC model

        if (near_contact_line.eq.1) then

         if ((fort_denconst(im_fluid1).gt.zero).and. &
             (fort_denconst(im_fluid2).gt.zero)) then
 
          if ((sin_angle.ge.zero).and.(cos_angle.ge.zero)) then
           angle_im1=asin(sin_angle)
           ! sin_angle=sin(a)  cos_angle=cos(a)
           ! a=pi-asin(sin_angle)
          else if ((sin_angle.ge.zero).and.(cos_angle.le.zero)) then
           angle_im1=Pi-asin(sin_angle)
          else
           print *,"sin_angle or cos_angle invalid"
           stop
          endif

          ZEYU_imodel=1 ! GNBC
          ZEYU_ifgnbc=1 ! GNBC
          ZEYU_lambda=8.0D-7  ! slip length
          ZEYU_lambda=LOW%dxmin  ! slip length
          ZEYU_l_macro=LOW%dxmin
          ZEYU_l_micro=EPS9
          ZEYU_dgrid=LOW%dxmin 
          ZEYU_d_closest=abs(dist_to_CL)

          if (fort_denconst(im_fluid1).ge. &
              fort_denconst(im_fluid2)) then
           im_liquid=im_fluid1
           im_vapor=im_fluid2
           ZEYU_thet_s=angle_im1  ! thet_s in the liquid.
          else if (fort_denconst(im_fluid2).ge. &
                   fort_denconst(im_fluid1)) then
           im_liquid=im_fluid2
           im_vapor=im_fluid1
           ZEYU_thet_s=Pi-angle_im1
          else
           print *,"fort_denconst bust"
           stop
          endif
          ZEYU_mu_l=fort_viscconst(im_liquid)
          ZEYU_mu_g=fort_viscconst(im_vapor)
          ZEYU_sigma=user_tension(iten)
          if (im_primary_image.eq.im_liquid) then
           ZEYU_thet_d_apparent=angle_ACT
          else if (im_primary_image.eq.im_vapor) then
           ZEYU_thet_d_apparent=Pi-angle_ACT
          else
           print *,"im_primary_image or im_vapor invalid"
           stop
          endif
          ZEYU_u_cl=zero
          ZEYU_u_slip=zero
          ZEYU_thet_d=ZEYU_thet_d_apparent

          call dynamic_contact_angle(ZEYU_mu_l, ZEYU_mu_g, ZEYU_sigma, &
           ZEYU_thet_s, &
           ZEYU_imodel, ZEYU_ifgnbc, ZEYU_lambda, &
           ZEYU_l_macro, ZEYU_l_micro, &
           ZEYU_dgrid, ZEYU_d_closest, ZEYU_thet_d_apparent, &
           ZEYU_u_cl, &
           ZEYU_u_slip, &
           ZEYU_thet_d)

          ! nCL is normal to the contact line in the substrate plane.
          ! nCL points to the im_primary_image material.
          ! ZEYU_u_slip is positive if the contact line is advancing into
          ! the gas.
          ! NOTE: if the velocity, ZEYU_u_slip=0.0, the interface might
          ! still move slowly since the curvature will not be numerically
          ! a constant on the interface.  This is what people call 
          ! "parasitic currents" when the interface moves due to surface
          ! tension, even though the curvature = constant.

          nCL_dot_n_raster=zero
          do dir=1,SDIM
           nCL_dot_n_raster=nCL_dot_n_raster+nCL(dir)*LOW%n_raster(dir)
          enddo
          mag=zero
          do dir=1,SDIM
           nCL_raster(dir)=nCL(dir)-nCL_dot_n_raster*LOW%n_raster(dir)
           mag=mag+nCL_raster(dir)**2
          enddo
          mag=sqrt(mag)
          if (mag.gt.zero) then
           do dir=1,SDIM
            nCL_raster(dir)=nCL_raster(dir)/mag
           enddo
          else if (mag.eq.zero) then
           ! do nothing
          else
           print *,"mag is NaN (getGhostVel): ",mag
           stop
          endif
          do dir=1,SDIM
           if (dir.ne.data_dir+1) then
            if (im_primary_image.eq.im_liquid) then
             ughost_tngt(dir)=-ZEYU_u_slip*nCL_raster(dir)
            else if (im_primary_image.eq.im_vapor) then
             ughost_tngt(dir)=ZEYU_u_slip*nCL_raster(dir)
            else
             print *,"im_primary_image or im_vapor invalid"
             stop
            endif
            if (LOW%n_raster(dir).eq.zero) then
             ! do nothing
            else
             print *,"LOW%n_raster(dir) invalid"
             stop
            endif
           else if (dir.eq.data_dir+1) then
            if (abs(nCL_raster(dir)).le.VOFTOL) then
             ! do nothing
            else
             print *,"nCL_raster(dir) invalid"
             stop
            endif
            if (abs(LOW%n_raster(dir)).eq.one) then
             ! do nothing
            else
             print *,"LOW%n_raster(dir) invalid"
             stop
            endif
           else
            print *,"dir or data_dir corrupt"
            stop
           endif
          enddo ! dir=1..sdim
         else
          print *,"fort_denconst invalid"
          stop
         endif
         if (DEBUG_DYNAMIC_CONTACT_ANGLE.eq.1) then
          print *,"xcrossing ",xcrossing(1),xcrossing(2),xcrossing(SDIM)
          print *,"xtriple ",xtriple(1),xtriple(2),xtriple(SDIM)
          print *,"nrm_solid ",nrm_solid(1),nrm_solid(2),nrm_solid(SDIM)
          print *,"nrm_fluid ",nrm_fluid(1),nrm_fluid(2),nrm_fluid(SDIM)
          print *,"im_fluid,angle_ACT(rad,deg) ",im_fluid,angle_ACT, &
                    angle_ACT*180.0d0/Pi
          print *,"im_liquid,ZEYU_thet_d_apparent(rad,deg) ",im_liquid, &
               ZEYU_thet_d_apparent,ZEYU_thet_d_apparent*180.0d0/Pi
          print *,"dx(1),dist_to_CL ",LOW%dx(1),dist_to_CL
          print *,"im_primary_image,im_secondary_image ", &
               im_primary_image,im_secondary_image
          print *," ZEYU_thet_s(rad,deg) ", &
           ZEYU_thet_s,ZEYU_thet_s*180.0d0/Pi
          print *,"nCL_raster=",nCL_raster(1),nCL_raster(2),nCL_raster(SDIM)
          print *,"ughost_tngt=",ughost_tngt(1),ughost_tngt(2), &
                  ughost_tngt(SDIM)
         endif

        else if (near_contact_line.eq.0) then
          ! ghost velocity lives *on* the rasterized interface.
         do dir=1,SDIM
          ughost_tngt(dir)=zero
         enddo
        else
         print *,"near_contact_line invalid"
         stop
        endif

        if (debug_slip_velocity_enforcement.eq.1) then
         ughost_tngt(1)=ten
         ughost_tngt(2)=zero
         ughost_tngt(SDIM)=zero
        endif

       else
        print *,"law_of_the_wall invalid"
        stop
       endif

       do dir=1,SDIM
        if (dir.eq.data_dir+1) then
         usolid_law_of_wall(dir)=LOW%usolid_raster(dir)
        else
         usolid_law_of_wall(dir)=ughost_tngt(dir)+LOW%usolid_raster(dir)
        endif
       enddo
      
       deallocate(user_tension)
 
      end subroutine getGhostVel

! nfree points towards the liquid (material im)
! nsolid points away from the solid 
! (nsolid is gradient of psi where psi<0 in solid)
! angle=pi => dry  contact=0 => wetting
! nghost points towards the liquid
! nperp points into the solid tangent to the contours of the extrapolated
!  level set functions.
      subroutine ghostnormal(nfree,nsolid,cos_angle,nghost,nperp)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: nfree(SDIM)
      real(amrex_real), INTENT(in) :: nsolid(SDIM)
      real(amrex_real), INTENT(out) :: nghost(SDIM)
      real(amrex_real), INTENT(out) :: nperp(SDIM)
      real(amrex_real), INTENT(in)  :: cos_angle

      real(amrex_real) e2(3),e3(3),ntemp(3)
      integer i
      real(amrex_real) ss,sin_angle,dist
      real(amrex_real) nsolid_new(SDIM)
      integer dir
      real(amrex_real) mag


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

      dist=sqrt( e2(1)**2+e2(2)**2+e2(3)**2+EPS14 )
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

      dist=sqrt( e3(1)**2+e3(2)**2+e3(3)**2+EPS14 )
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
       dist=sqrt( nghost(1)**2+nghost(2)**2+nghost(SDIM)**2+EPS14 )
      else if (SDIM.eq.2) then
       dist=sqrt( nghost(1)**2+nghost(2)**2+EPS14 )
      else
       print *,"dimension bust"
       stop
      endif
      do i=1,SDIM
       nghost(i)=nghost(i)/dist
      enddo

      if (SDIM.eq.3) then
       dist=sqrt( nperp(1)**2+nperp(2)**2+nperp(SDIM)**2+EPS14 )
      else if (SDIM.eq.2) then
       dist=sqrt( nperp(1)**2+nperp(2)**2+EPS14 )
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

      integer, INTENT(in) :: i,j,k
      real(amrex_real), INTENT(in), pointer :: datafab(D_DECL(:,:,:))
      real(amrex_real), INTENT(out) :: data_out
      integer datalo,datahi
      integer dir
      integer idata(3)

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

      integer, INTENT(in) :: i,j,k,n
      real(amrex_real), INTENT(in), pointer :: datafab(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out) :: data_out
      integer datalo,datahi
      integer dir
      integer idata(3)

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

      subroutine safe_data3D(i,j,k,n,datafab,data_out)
      IMPLICIT NONE

      integer, INTENT(in) :: i,j,k,n
      real(amrex_real), INTENT(in), pointer :: datafab(:,:,:,:)
      real(amrex_real), INTENT(out) :: data_out
      integer datalo,datahi
      integer dir
      integer idata(3)

      idata(1)=i
      idata(2)=j
      idata(3)=k
      do dir=1,3
       datalo=LBOUND(datafab,dir)
       datahi=UBOUND(datafab,dir)
       if (idata(dir).lt.datalo) then
        idata(dir)=datalo
       endif
       if (idata(dir).gt.datahi) then
        idata(dir)=datahi
       endif
      enddo ! dir=1..3
      data_out=datafab(idata(1),idata(2),idata(3),n)

      return
      end subroutine safe_data3D

      subroutine safe_data_index(i,j,k,i_safe,j_safe,k_safe,datafab)
      IMPLICIT NONE

      integer, INTENT(in) :: i,j,k
      integer, INTENT(out) :: i_safe,j_safe,k_safe
      real(amrex_real), INTENT(in), pointer :: datafab(D_DECL(:,:,:),:)
      integer datalo,datahi
      integer dir
      integer idata(3)

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

      integer function fort_pattern_test(source,source_len, &
          pattern,pattern_len)
      IMPLICIT NONE

      integer, INTENT(in) :: source_len,pattern_len
      CHARACTER(len=source_len), INTENT(in) :: source
      CHARACTER(len=pattern_len), INTENT(in) :: pattern
      integer i,j,local_test,source_char,pattern_char

      fort_pattern_test=0
      do i=1,source_len-pattern_len+1

       local_test=1
       do j=1,pattern_len
        source_char=ICHAR(source(i+j-1:i+j-1))
        pattern_char=ICHAR(pattern(j:j))
        if (source_char.eq.pattern_char) then
         ! do nothing
        else if (source_char.ne.pattern_char) then
         local_test=0
        else
         print *,"source_char or pattern_char bust"
         stop
        endif
       enddo !j=1,pattern_len
       if (local_test.eq.1) then
        fort_pattern_test=1
       endif
      enddo !enddo i=1,source_len-pattern_len+1

      return
      end function fort_pattern_test

      ! Added by Guibo 11-12-2012
      subroutine dumpstring(instring)
      implicit none
      character*80 instring
      integer len, ii, I

      len = LEN_TRIM(instring)
      do ii = 1,len
       I = ICHAR(instring(ii:ii))
       write(11) I
      end do
      write(11) 0

      return
      end subroutine dumpstring

      integer function is_compressible_mat(im)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im out of range"
       stop
      endif

       ! MAX_NUM_EOS=998
      if ((fort_material_type(im).ge.1).and. &
          (fort_material_type(im).le.MAX_NUM_EOS)) then
       is_compressible_mat=1
      else if ((fort_material_type(im).eq.0).or. &
               (fort_material_type(im).eq.999)) then
       is_compressible_mat=0
      else
       print *,"fort_material_type invalid: ",im,fort_material_type(im)
       stop
       is_compressible_mat=0
      endif

      return
      end function is_compressible_mat

      subroutine debug_EOS(im)
      use probcommon_module
      IMPLICIT NONE

      integer im
      integer verbose_EOS
      integer mat_type
      integer nden
      integer i,iden
      real(amrex_real) temperature,denlo,denhi,den
      real(amrex_real) internal_energy
      real(amrex_real) pressure
      real(amrex_real) soundsqr
      real(amrex_real) massfrac_parm(num_species_var+1)
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

      real(amrex_real) function get_user_temperature(time,bcflag,im)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      real(amrex_real) startup_time
      integer, INTENT(in) :: bcflag,im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im out of range"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in get_user_temperature: ",time
      else if (time.lt.zero) then
       print *,"time invalid in get_user_temperature: ",time
       stop
      else
       print *,"time bust in get_user_temperature: ",time
       stop
      endif

      if ((fort_initial_temperature(im).gt.zero).and. &
          (fort_tempconst(im).gt.zero)) then
       ! do nothing
      else
       print *,"fort_initial_temperature(im) invalid or"
       print *,"fort_tempconst(im) invalid"
       print *,"im,fort_initial_temperature(im),fort_tempconst(im) ", &
        im,fort_initial_temperature(im),fort_tempconst(im)
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
         print *,"time invalid in get_user_temperature: ",time
         stop
        endif
       else
        print *,"startup_time or im invalid: ",startup_time,im
        stop
       endif
      else

       if (bcflag.eq.0) then
        get_user_temperature=fort_initial_temperature(im)
       else if (bcflag.eq.1) then
        get_user_temperature=fort_tempconst(im)
       else
        print *,"bcflag invalid: ",bcflag
        stop
       endif

      endif

      return
      end function get_user_temperature


        ! temperature initialized with a default value
        ! called from denBC, process_initdata,
        ! GENERAL_PHASE_CHANGE_STATE,
        ! probtype.eq.59, probtype.eq.55, or probtype.eq.710
        ! bcflag==0 => called from initdata or GENERAL_PHASE_CHANGE_STATE
        ! bcflag==1 => called from denbc
      subroutine outside_temperature(time,x,y,z,temperature,im,bcflag)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: im,bcflag
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: x,y,z
      real(amrex_real), INTENT(out) :: temperature
      real(amrex_real) :: temp_slope
      real(amrex_real) substrate_height
      real(amrex_real) ice_vertical
      integer im_solid_temperature
      real(amrex_real) z_shift
      real(amrex_real) seed_thickness

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in outside_temperature: ",time
      else if (time.lt.zero) then
       print *,"time invalid in outside_temperature: ",time
       stop
      else
       print *,"time bust in outside_temperature: ",time
       stop
      endif

      im_solid_temperature=im_solid_primary()

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"expecting z=y in 2d routine outside_temperature: ",y,z
        stop
       endif
      endif
      z_shift=z
 
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid62"
       stop
      endif
      if ((bcflag.ne.0).and. & ! called from initdata
          (bcflag.ne.1)) then  ! called from denBC (density, temp, spec bc)
       print *,"bcflag invalid: ",bcflag
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

        if (abs(substrate_height-(ice_vertical-half*radblob)).le.EPS2) then
         !do nothing
        else
         print *,"init bot of water+ice block should coincide with substrate"
         stop
        endif
        if ((radblob3.lt.radblob).and. &
            (radblob3.ge.zero)) then
         !do nothing
        else
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
           
           ! block of ice melting on substrate
          if (z_shift.ge.substrate_height) then ! above the substrate
           temperature=get_user_temperature(time,bcflag,3) ! ice
          else if ((z_shift.ge.zero).and. &
                   (z_shift.le.substrate_height)) then
           temperature= &
            get_user_temperature(time,bcflag,im_solid_temperature)+ &
            (get_user_temperature(time,bcflag,3)- &
             get_user_temperature(time,bcflag,im_solid_temperature))* &
            z_shift/substrate_height 
          else if (z_shift.le.zero) then
           temperature= &
            get_user_temperature(time,bcflag,im_solid_temperature) !substrate
          else
           print *,"z_shift failure"
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
         if (z_shift.ge.substrate_height+radblob3) then
          temperature=get_user_temperature(time,bcflag,3) ! ice region
         else if (z_shift.ge.substrate_height) then
          temperature=get_user_temperature(time,bcflag,3) ! water region
         else if ((z_shift.ge.zero).and. &
                  (z_shift.le.substrate_height)) then
          temperature= &
           get_user_temperature(time,bcflag,im_solid_temperature)+ &
           (get_user_temperature(time,bcflag,3)- &
            get_user_temperature(time,bcflag,im_solid_temperature))* &
            z_shift/substrate_height
         else if (z_shift.le.zero) then
          ! substrate
          temperature=get_user_temperature(time,bcflag,im_solid_temperature) 
         else
          print *,"z_shift failure"
          stop
         endif
        else
         print *,"bcflag invalid: ",bcflag
         stop
        endif
       endif ! substrate_height>0

       ! in: outside_temperature
      else if (probtype.eq.55) then

       if (SDIM.eq.2) then
        substrate_height=yblob2
       else if (SDIM.eq.3) then 
        substrate_height=zblob2
       else
        print *,"dimension bust: ",SDIM
        stop
       endif

       if (axis_dir.eq.5) then ! freezing: solid, ice, water, air

        seed_thickness=radblob3

         ! substrate: 0<y<substrate_height (yblob2 or zblob2)
        if (substrate_height.gt.zero) then  

         if (zblob4.eq.substrate_height) then !transition thermal layer
          ! do nothing
         else if (zblob4.eq.zero) then ! no transition, T=T_substrate
          ! do nothing
         else
          print *,"zblob4 invalid: ",zblob4
          print *,"substrate_height: ",substrate_height
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
            print *,"im_solid_temperature invalid: ",im_solid_temperature
            stop
           endif
           if (num_materials.lt.4) then
            print *,"num_materials invalid: ",num_materials
            stop
           endif

           if (z_shift.ge.substrate_height) then

            if (zblob4.eq.zero) then
             temperature=get_user_temperature(time,bcflag,im_solid_temperature)
            else if (zblob4.eq.substrate_height) then
             temperature=get_user_temperature(time,bcflag,3) ! ice
            else
             print *,"zblob4 invalid: ",zblob4
             stop
            endif

           else if ((z_shift.ge.zero).and. &
                    (z_shift.le.substrate_height)) then

            if (zblob4.eq.zero) then
             temperature=get_user_temperature(time,bcflag,im_solid_temperature)
            else if (zblob4.eq.substrate_height) then
             temperature= &
              get_user_temperature(time,bcflag,im_solid_temperature)+ &
              (get_user_temperature(time,bcflag,3)- &
               get_user_temperature(time,bcflag,im_solid_temperature))* &
              z_shift/substrate_height 
            else
             print *,"zblob4 invalid: ",zblob4
             stop
            endif

           else if (z_shift.le.zero) then

            temperature= &
             get_user_temperature(time,bcflag,im_solid_temperature) !substrate

           else
            print *,"z_shift failure: ",z_shift
            stop
           endif
          else
           print *,"im invalid64: ",im
           stop
          endif

         else if (bcflag.eq.1) then ! called from denBC

          if (im_solid_temperature.ne.4) then
           print *,"im_solid_temperature invalid"
           stop
          endif
           ! radblob3=thickness of underside of drop that is already
           ! frozen.
          if (z_shift.ge.substrate_height+seed_thickness) then
           temperature=get_user_temperature(time,bcflag,1) ! water
          else if (z_shift.ge.substrate_height) then
           temperature=get_user_temperature(time,bcflag,3) ! ice
          else if ((z_shift.ge.zero).and. &
                   (z_shift.le.substrate_height)) then

           if (zblob4.eq.zero) then
            temperature=get_user_temperature(time,bcflag,im_solid_temperature)
           else if (zblob4.eq.substrate_height) then
            temperature= &
             get_user_temperature(time,bcflag,im_solid_temperature)+ &
             (get_user_temperature(time,bcflag,3)- &
              get_user_temperature(time,bcflag,im_solid_temperature))* &
             z_shift/substrate_height
           else
            print *,"zblob4 invalid"
            stop
           endif

          else if (z_shift.le.zero) then
           ! substrate
           temperature=get_user_temperature(time,bcflag,im_solid_temperature) 
          else
           print *,"z_shift failure"
           stop
          endif
         else
          print *,"bcflag invalid: ",bcflag
          stop
         endif
         ! substrate: 0<y<substrate_height
        else if (substrate_height.eq.zero) then
         ! do nothing
        else
         print *,"not expecting substrate_height<0"
         print *,"probtype,axis_dir,bcflag ",probtype,axis_dir,bcflag
         print *,"im=",im
         print *,"time=",time
         print *,"z_shift=",z_shift
         print *,"substrate_height=",substrate_height
         stop
        endif ! substrate_height>0

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
         print *,"radblob2 invalid: ",radblob2
        else if (radblob2.eq.zero) then
         ! do nothing
        else if (radblob2.gt.zero) then ! angle of incline
         if (levelrz.eq.COORDSYS_RZ) then
          if ((SDIM.ne.2).or.(xblob2.ne.zero)) then
           print *,"parameters not supported"
           stop
          endif
          z_shift=substrate_height+ &
                  (z-substrate_height)*cos(radblob2)-x*sin(radblob2)
         else if (levelrz.eq.COORDSYS_CARTESIAN) then
          z_shift=substrate_height+ &
                  (z-substrate_height)*cos(radblob2)-(x-xblob2)*sin(radblob2)
         else
          print *,"levelrz not supported"
          stop
         endif
        else
         print *,"radblob2 invalid"
         stop
        endif

         ! (xblob2,yblob2,zblob2) is the "center" of the heated plate.
        if (z_shift.ge.substrate_height+yblob3) then
         temperature=get_user_temperature(time,bcflag,1)
        else if (z_shift.le.substrate_height) then
         temperature=get_user_temperature(time,bcflag,im_solid_temperature)
        else if ((substrate_height.le.z_shift).and. &
                 (z_shift.le.substrate_height+yblob3)) then
         temp_slope=(get_user_temperature(time,bcflag,im_solid_temperature)- &
                     get_user_temperature(time,bcflag,1))/yblob3
         temperature= &
          get_user_temperature(time,bcflag,im_solid_temperature)- &
          temp_slope*(z_shift-substrate_height)
        else
         print *,"z_shift or substrate_height invalid"
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
       if (yblob3.gt.zero) then
        ! do nothing
       else
        print *,"yblob3 invalid"
        stop
       endif
       if (z_shift.le.yblob3) then
        temperature=get_user_temperature(time,bcflag,im_solid_temperature)
       endif

      else
       print *,"expecting probtype=55, 59, 414, 2001, or 710 "
       print *,"in outside_temperature"
       print *,"probtype: ",probtype
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

      integer, INTENT(out) :: FSI_prescribed_flag
      integer, INTENT(in) :: color,dir,blob_array_size,num_colors
      real(amrex_real), INTENT(in) :: xrigid(SDIM)
      real(amrex_real) xrigid3D(3)
      real(amrex_real), INTENT(out) :: vel
      real(amrex_real), INTENT(in) :: blob_array(blob_array_size)
      integer ibase
      integer veldir,irow
      real(amrex_real) blob_mass_for_velocity(3)
      real(amrex_real) blob_center(3)
      real(amrex_real) phi_N(3)
      real(amrex_real) vel_local(SDIM) 

      if (blob_array_size.ne.num_colors*num_elements_blobclass) then
       print *,"blob_array_size invalid"
       stop
      endif

      FSI_prescribed_flag=0

      if ((color.ge.1).and.(color.le.num_colors)) then

       if ((dir.ge.1).and.(dir.le.SDIM)) then

        ibase=(color-1)*num_elements_blobclass
        do veldir=1,3
         blob_mass_for_velocity(veldir)=blob_array(ibase+BLB_MASS_VEL+veldir)
        enddo
        do veldir=1,3
         xrigid3D(veldir)=zero
         blob_center(veldir)=zero
        enddo
        do veldir=1,SDIM
         blob_center(veldir)=blob_array(ibase+BLB_CEN_ACT+veldir)
         xrigid3D(veldir)=xrigid(veldir)
        enddo
        do veldir=1,SDIM
         vel_local(veldir)=zero
        enddo 
        do irow=1,2*SDIM
         call init_basis(blob_center,xrigid3D,irow,phi_N)
         do veldir=1,SDIM
           !ice touches a solid
          if (blob_mass_for_velocity(2).gt.zero) then 
           vel_local(veldir)=vel_local(veldir)+ &
            blob_array(ibase+BLB_VEL+2*SDIM+irow)*phi_N(veldir)
           FSI_prescribed_flag=1
          else if (blob_mass_for_velocity(2).eq.zero) then
           if (blob_mass_for_velocity(1).gt.zero) then
            vel_local(veldir)=vel_local(veldir)+ &
             blob_array(ibase+BLB_VEL+irow)*phi_N(veldir)
           else if (blob_mass_for_velocity(1).eq.zero) then
            if (blob_mass_for_velocity(3).gt.zero) then
             vel_local(veldir)=vel_local(veldir)+ &
              blob_array(ibase+BLB_VEL+2*(2*SDIM)+irow)*phi_N(veldir)
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

      integer nbasis
      real(amrex_real) rigid_centroid(3)
      real(amrex_real) x0(3)
      real(amrex_real) phi_N(3)
      integer dir

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

      integer dir,i,j,k
      integer loface(SDIM),hiface(SDIM)
      integer velbc(SDIM,2)
      integer sten_flag,idx,side

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

      real(amrex_real) xcross,xcen,xside,LScen,LSside

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

      real(amrex_real) xtarget(SDIM)
      real(amrex_real) xfoot(SDIM)
      real(amrex_real) xtargetRT(SDIM)
      real(amrex_real) xfootRT(SDIM)
      real(amrex_real) gradphi(SDIM) 
      real(amrex_real) phi,newphi
      real(amrex_real) mag
      integer dir

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       newphi=phi
      else if (levelrz.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       newphi=phi
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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
      subroutine prepare_normal(gradphi,rval,mag,sdim_in)
      use probcommon_module

      IMPLICIT NONE

      integer, intent(in) :: sdim_in
      real(amrex_real), intent(inout) :: gradphi(sdim_in)
      real(amrex_real), intent(in)    :: rval
      real(amrex_real), intent(out)   :: mag
      integer dir

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (sdim_in.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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
      do dir=1,sdim_in
       mag=mag+gradphi(dir)**2
      enddo
      mag=sqrt(mag)
      if (mag.gt.zero) then
       do dir=1,sdim_in
        gradphi(dir)=gradphi(dir)/mag
       enddo
      else if (mag.eq.zero) then
       ! do nothing
      else
       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:4476"
       print *,"mag is NaN (prepare_normal): ",mag
       stop
      endif

      end subroutine prepare_normal
 
      subroutine normalize_vector(vin)
      IMPLICIT NONE

      real(amrex_real), INTENT(inout) :: vin(SDIM)
      integer dir
      real(amrex_real) mag

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

      subroutine normalize_LS_normals(LS)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), INTENT(inout) :: LS(num_materials*(1+SDIM))
      integer :: im,dir
      real(amrex_real) :: local_normal(SDIM)

      do im=1,num_materials
       do dir=1,SDIM
        local_normal(dir)=LS(num_materials+(im-1)*SDIM+dir)
       enddo
       call normalize_vector(local_normal)
       do dir=1,SDIM
        LS(num_materials+(im-1)*SDIM+dir)=local_normal(dir)
       enddo
      enddo ! im=1..num_materials

      return 
      end subroutine normalize_LS_normals

      subroutine crossprod2d(a,b,c)
      IMPLICIT NONE
      
      real(amrex_real), INTENT(in) :: a(SDIM),b(SDIM)
      real(amrex_real), INTENT(out) :: c(3)

      c(1)=zero
      c(2)=zero
      c(3)=a(1)*b(2)-b(1)*a(2)

      return
      end subroutine crossprod2d

      subroutine crossprod(a,b,c)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: a(3),b(3)
      real(amrex_real), INTENT(out) :: c(3)

      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=-(a(1)*b(3)-b(1)*a(3))
      c(3)=a(1)*b(2)-b(1)*a(2)

      return
      end subroutine crossprod

subroutine global_xdist(x1,x2,dist)
IMPLICIT NONE
 
real(amrex_real), dimension(3),INTENT(in) :: x1,x2
real(amrex_real), INTENT(out) :: dist

 dist=sqrt( (x1(1)-x2(1))**2+ &
            (x1(2)-x2(2))**2+ &
            (x1(3)-x2(3))**2 )

return
end subroutine global_xdist

subroutine global_checkinline(nnode,xnode,tol,xc, &
             inplane,unsigned_mindist,xclosest,normal_closest)
IMPLICIT NONE

integer, INTENT(inout) :: inplane
real(amrex_real), INTENT(inout) :: unsigned_mindist
real(amrex_real), dimension(3), INTENT(inout) :: xclosest
real(amrex_real), dimension(3), INTENT(inout) :: normal_closest
real(amrex_real), INTENT(in) :: tol
real(amrex_real), dimension(2,3), INTENT(in) :: xnode,nnode
real(amrex_real), dimension(3), INTENT(in) :: xc
real(amrex_real) :: dottop,dotbot,t,curdist
real(amrex_real), dimension(3) :: xnot,normal
integer :: dir

   ! x(t)=x2+t*(x1-x2)
   ! C(t)=||x(t)-xc||^2
   ! C'(t)=0
   ! (x(t)-xc) dot (x1-x2)=0
   ! (-(xc-x2)+t(x1-x2)) dot (x1-x2)=0
   ! t=(xc-x2) dot (x1-x2)/||x1-x2||^2
   ! since x(t)=x2+t (x1-x2), if t=0 => x=x2
   ! if t=1 => x=x1
  dottop=zero
  dotbot=zero
  do dir=1,3
   dottop=dottop+(xnode(2,dir)-xnode(1,dir))*(xnode(2,dir)-xc(dir))
   dotbot=dotbot+(xnode(2,dir)-xnode(1,dir))**2
  enddo
  if (dotbot.eq.zero) then
   ! do nothing
  else if (dotbot.gt.zero) then
   t=dottop/dotbot
   if ((t.ge.-tol).and.(t.le.one+tol)) then

    t=min(t,one-tol)
    t=max(t,tol)

    do dir=1,3
     xnot(dir)=t*xnode(1,dir)+(one-t)*xnode(2,dir)
     normal(dir)=t*nnode(1,dir)+(one-t)*nnode(2,dir)
    enddo

    call global_xdist(xnot,xc,curdist)
    if (curdist.ge.zero) then
     ! do nothing
    else
     print *,"curdist invalid"
     stop
    endif

    if ((curdist.lt.unsigned_mindist).or. &
        (inplane.eq.0)) then
     inplane=1
     unsigned_mindist=curdist
     do dir=1,3
      xclosest(dir)=xnot(dir)
      normal_closest(dir)=normal(dir)
     enddo
    else if ((curdist.ge.unsigned_mindist).and. &
             (inplane.eq.1)) then
     ! do nothing
    else
     print *,"curdist or inplane invalid"
     stop
    endif
   else if ((t.lt.-tol).or. &
            (t.gt.one+tol)) then
    ! do nothing
   else 
    print *,"t is NaN"
    stop
   endif 
  else
   print *,"dotbot is NaN"
   stop
  endif

return
end subroutine global_checkinline

subroutine global_checkinplane(xnode,xclosest,tol, &
              xclosest_project,inplane)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: tol
real(amrex_real), dimension(3), INTENT(in) :: xclosest
real(amrex_real), dimension(3), INTENT(out) :: xclosest_project
integer, INTENT(out) :: inplane
real(amrex_real), dimension(3,3), INTENT(in) :: xnode ! (ipoint,dir)
real(amrex_real), dimension(3,3) :: AA,AINVERSE
real(amrex_real) :: det
real(amrex_real) :: tx_sum,tx_sum_new
real(amrex_real), dimension(3) :: tx
real(amrex_real), dimension(3) :: tx_project
real(amrex_real), dimension(3) :: v1,v2,v1xv2
integer :: dir,i,j,k
real(amrex_real) :: sanity_tol

 ! xnode(1)-xnode(1) is mapped to (0,0,0)
 ! xnode(2)-xnode(1)=v1 is mapped to (1,0,0)
 ! xnode(3)-xnode(1)=v2 is mapped to (0,1,0) 
 ! v1 x v2 is mapped to (0,0,1)
 ! Let A map from unit space to real space
 ! A (0,0,0) = (0,0,0)
 ! A (1,0,0) = v1 => first column of A is v1
 ! A (0,1,0) = v2 => second column of A is v2
 ! A (0,0,1) = v1 x v2 => third column of A is v1 x v2
 ! A^{-1} maps from real space back to unit space

 if ((tol.ge.zero).and.(tol.lt.one)) then
  ! do nothing
 else
  print *,"tol out of range: ",tol
  stop
 endif

 inplane=1

 do dir=1,3
  v1(dir)=xnode(2,dir)-xnode(1,dir)
  v2(dir)=xnode(3,dir)-xnode(1,dir)
 enddo  
 v1xv2(1)=v1(2)*v2(3)-v1(3)*v2(2)  
 v1xv2(2)=v1(3)*v2(1)-v1(1)*v2(3)  
 v1xv2(3)=v1(1)*v2(2)-v1(2)*v2(1)  
 do dir=1,3
  AA(dir,1)=v1(dir)
  AA(dir,2)=v2(dir)
  AA(dir,3)=v1xv2(dir)
 enddo
  
 det=AA(1,1)*(AA(2,2)*AA(3,3)-AA(2,3)*AA(3,2))- &
     AA(1,2)*(AA(2,1)*AA(3,3)-AA(2,3)*AA(3,1))+ &
     AA(1,3)*(AA(2,1)*AA(3,2)-AA(2,2)*AA(3,1))

 if (abs(det).eq.zero) then
  inplane=0
 else if (abs(det).gt.zero) then
  AINVERSE(1,1)=+(AA(2,2)*AA(3,3)-AA(2,3)*AA(3,2))
  AINVERSE(2,1)=-(AA(2,1)*AA(3,3)-AA(2,3)*AA(3,1))
  AINVERSE(3,1)=+(AA(2,1)*AA(3,2)-AA(2,2)*AA(3,1))
  AINVERSE(1,2)=-(AA(1,2)*AA(3,3)-AA(3,2)*AA(1,3))
  AINVERSE(2,2)=+(AA(1,1)*AA(3,3)-AA(1,3)*AA(3,1))
  AINVERSE(3,2)=-(AA(1,1)*AA(3,2)-AA(3,1)*AA(1,2))
  AINVERSE(1,3)=+(AA(1,2)*AA(2,3)-AA(2,2)*AA(1,3))
  AINVERSE(2,3)=-(AA(1,1)*AA(2,3)-AA(2,1)*AA(1,3))
  AINVERSE(3,3)=+(AA(1,1)*AA(2,2)-AA(1,2)*AA(2,1))
  do i=1,3
  do j=1,3
   AINVERSE(i,j)=AINVERSE(i,j)/det
  enddo
  enddo
  
  do i=1,3
   tx(i)=0.0
   do k=1,3
    tx(i)=tx(i)+AINVERSE(i,k)*(xclosest(k)-xnode(1,k))
   enddo
  enddo
 
  if ((tx(1).lt.-tol).or. &
      (tx(1).gt.one+tol).or. &
      (tx(2).lt.-tol).or. &
      (tx(2).gt.one+tol).or. &
      (tx(1)+tx(2).gt.one+tol)) then
   inplane=0
  else if ((tx(1).ge.-tol).and. &
           (tx(1).le.one+tol).and. &
           (tx(2).ge.-tol).and. &
           (tx(2).le.one+tol).and. &
           (tx(1)+tx(2).le.one+tol)) then

   sanity_tol=VOFTOL*max(1.0d0,1.0d0/det)

   if (abs(tx(3)).le.sanity_tol) then
    ! do nothing
   else if (abs(tx(3)).ge.sanity_tol) then
    print *,"something wrong with transformation"
    print *,"tx(1)= ",tx(1)
    print *,"tx(2)= ",tx(2)
    print *,"tx(3)= ",tx(3)
    print *,"tol= ",tol
    print *,"VOFTOL= ",VOFTOL
    print *,"det= ",det
    print *,"sanity_tol= ",sanity_tol
    stop
   else
    print *,"tx(3) is NaN"
    stop
   endif

   tx(1)=min(tx(1),one-tol)
   tx(1)=max(tx(1),tol)
   tx(2)=min(tx(2),one-tol)
   tx(2)=max(tx(2),tol)
   tx_sum=tx(1)+tx(2)
   if (tx_sum.gt.zero) then
    tx_sum_new=tx_sum
    tx_sum_new=min(tx_sum_new,one-tol)
    tx(1)=tx(1)*tx_sum_new/tx_sum
    tx(2)=tx(2)*tx_sum_new/tx_sum
   else if (tx_sum.eq.zero) then
    tx_sum_new=zero
    tx(1)=zero
    tx(2)=zero
   else
    print *,"tx_sum is NaN"
    stop
   endif
   tx(3)=zero

   do i=1,3
    xclosest_project(i)=xnode(1,i)
    do k=1,3
     xclosest_project(i)=xclosest_project(i)+AA(i,k)*tx(k)
    enddo
   enddo

   do i=1,3
    tx_project(i)=0.0
    do k=1,3
     tx_project(i)=tx_project(i)+AINVERSE(i,k)*(xclosest_project(k)-xnode(1,k))
    enddo
   enddo

   do i=1,3
    if (abs(tx_project(i)-tx(i)).le.sanity_tol) then
     ! do nothing
    else if (abs(tx_project(i)-tx(i)).ge.sanity_tol) then
     print *,"sanity check failed tx_project and tx differ too much"
     print *,"i=",i
     print *,"tx(i)= ",tx(i)
     print *,"tx_project(i)= ",tx_project(i)
     print *,"tol= ",tol
     print *,"VOFTOL= ",VOFTOL
     print *,"det= ",det
     print *,"sanity_tol= ",sanity_tol
     stop
    else
     print *,"tx or tx_project are NaN"
     stop
    endif
   enddo !i=1,3

  else
   print *,"checking in triangle bust"
   stop
  endif
 else
  print *,"det is NaN"
  stop
 endif 

return
end subroutine global_checkinplane

subroutine check_outside_box(xcell,BB,LS,MASK)
IMPLICIT NONE

real(amrex_real), INTENT(in) :: xcell(3)
real(amrex_real), INTENT(in) :: BB(3,2)
integer, INTENT(out) :: MASK
real(amrex_real), INTENT(out) :: LS
real(amrex_real) :: lenscale
real(amrex_real) :: test_lenscale
integer outside_flag
integer dir

lenscale=zero
outside_flag=0

do dir=1,3
 if (BB(dir,2).ge.BB(dir,1)) then
  test_lenscale=(BB(dir,2)-BB(dir,1))/four
  if (test_lenscale.gt.lenscale) then
   lenscale=test_lenscale
  endif
  if (xcell(dir).le.BB(dir,1)) then
   outside_flag=1
  endif
  if (xcell(dir).ge.BB(dir,2)) then
   outside_flag=1
  endif
 else
  print *,"check outside, BB invalid dir,1,2: ",BB(dir,1),BB(dir,2)
  stop
 endif
enddo !dir=1..3

if (outside_flag.eq.1) then
 LS=-lenscale
 MASK=FSI_FINE_SIGN_VEL_VALID
else if (outside_flag.eq.0) then
 ! do nothing
else
 print *,"outside_flag invalid"
 stop
endif

return
end subroutine check_outside_box

subroutine check_inside_box(xcell,BB,LS,MASK)
IMPLICIT NONE

real(amrex_real), INTENT(in) :: xcell(3)
real(amrex_real), INTENT(in) :: BB(3,2)
integer, INTENT(out) :: MASK
real(amrex_real), INTENT(out) :: LS
real(amrex_real) :: lenscale
real(amrex_real) :: test_lenscale
integer inside_flag
integer dir

lenscale=zero
inside_flag=1

do dir=1,3
 if (BB(dir,2).ge.BB(dir,1)) then
  test_lenscale=(BB(dir,2)-BB(dir,1))/four
  if (test_lenscale.gt.lenscale) then
   lenscale=test_lenscale
  endif
  if (xcell(dir).lt.BB(dir,1)) then
   inside_flag=0
  endif
  if (xcell(dir).gt.BB(dir,2)) then
   inside_flag=0
  endif
 else
  print *,"check inside, BB invalid dir,1,2: ",BB(dir,1),BB(dir,2)
  stop
 endif
enddo !dir=1..3

if (inside_flag.eq.1) then
 LS=-lenscale
 MASK=FSI_FINE_SIGN_VEL_VALID
else if (inside_flag.eq.0) then
 ! do nothing
else
 print *,"inside_flag invalid"
 stop
endif

return
end subroutine check_inside_box


subroutine get_crse_index(i,j,k,ic,jc,kc,dir)
IMPLICIT NONE

integer, INTENT(in) :: i,j,k
integer, INTENT(out) :: ic,jc,kc
integer, INTENT(in) :: dir

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

subroutine maxind(k,S,n,m_out)
IMPLICIT NONE
integer, INTENT(in) :: n
integer, INTENT(in) :: k
real(amrex_real), INTENT(in) :: S(n,n)
integer, INTENT(out) :: m_out
integer :: i

if ((n.ge.2).and.(n.le.3)) then
 ! do nothing
else
 print *,"n out of range"
 stop
endif
if ((k.ge.1).and.(k.le.n)) then
 ! do nothing
else
 print *,"k out of range"
 stop
endif
m_out=k+1
do i=k+2,n
 if (abs(S(k,i))>abs(S(k,m_out))) then
  m_out=i
 endif
enddo

return
end subroutine maxind

subroutine EVAL_update(k,t,y,changed,evals,state,n)
IMPLICIT NONE

integer, INTENT(in) :: n
integer, INTENT(in) :: k
real(amrex_real), INTENT(inout) :: y
integer, INTENT(inout) :: state
integer, INTENT(inout) :: changed(n)
real(amrex_real), INTENT(inout) :: evals(n)
real(amrex_real), INTENT(in) :: t

if ((n.ge.2).and.(n.le.3)) then
 ! do nothing
else
 print *,"n out of range"
 stop
endif
if ((k.ge.1).and.(k.le.n)) then
 ! do nothing
else
 print *,"k invalid"
 stop
endif
y=evals(k) 
evals(k)=y+t
if ((changed(k).eq.1).and.(y.eq.evals(k))) then
 changed(k)=0
 state=state-1
else if ((changed(k).eq.0).and.(y.ne.evals(k))) then
 changed(k)=1
 state=state+1
else if ((changed(k).eq.1).and.(y.ne.evals(k))) then
 ! do nothing
else if ((changed(k).eq.0).and.(y.eq.evals(k))) then
 ! do nothing
else 
 print *,"changed or evals invalid"
 stop
endif

return
end subroutine EVAL_update

subroutine EVAL_rotate(k,l,i,j,S,n,sinrot,cosrot)
IMPLICIT NONE

integer, INTENT(in) :: n
real(amrex_real), INTENT(inout) :: S(n,n)
real(amrex_real), INTENT(in) :: sinrot,cosrot
integer, INTENT(in) :: k,l,i,j
real(amrex_real) skl,sij 

if ((n.ge.2).and.(n.le.3)) then
 ! do nothing
else
 print *,"expecting n=2,3"
 stop
endif

skl=S(k,l)
sij=S(i,j)
S(k,l)=cosrot*skl-sinrot*sij
S(i,j)=sinrot*skl+cosrot*sij

return
end subroutine EVAL_rotate

!columns of "evecs" are the eigenvectors.
!evecs(i,j) flattened index is s=n*(j-1)+i
subroutine fort_jacobi_eigenvalue(S,evals,evecs,n) &
bind(c,name='fort_jacobi_eigenvalue')
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: n
real(amrex_real), INTENT(in) :: S(n,n)
real(amrex_real), INTENT(out) :: evals(n)
real(amrex_real), INTENT(out) :: evecs(n,n)
real(amrex_real) :: S_MOD(n,n)
integer :: i,j,k,k_in,l,m,state
real(amrex_real) :: sinrot,cosrot,t,p,y,d,r
integer :: ind(n)
integer :: changed(n)
real(amrex_real) :: eik,eil
real(amrex_real) :: XL(n,n)
real(amrex_real) :: XLXT(n,n)
real(amrex_real) :: max_S
real(amrex_real) :: sanity_err
real(amrex_real) :: swap_hold
real(amrex_real) :: sum_off_diag

if (1.eq.0) then
 print *,"begin of jacobi eval"
endif

if ((n.ge.2).and.(n.le.3)) then
 ! do nothing
else
 print *,"expecting n=2,3"
 stop
endif

sum_off_diag=zero

do i=1,n
do j=1,n
 if (i.eq.j) then
  evecs(i,j)=one
 else if (i.ne.j) then
  evecs(i,j)=zero
  sum_off_diag=sum_off_diag+abs(S(i,j))
 else
  print *,"i or j corrupt"
  stop
 endif
 S_MOD(i,j)=S(i,j)
enddo
enddo

max_S=zero
do i=1,n
do j=1,n  
 if (abs(S_MOD(i,j)).gt.max_S) then
  max_S=abs(S_MOD(i,j))
 endif
enddo
enddo
max_S=max(max_S,one)

state=n

do k=1,n
 call maxind(k,S,n,ind(k))
 evals(k)=S(k,k)
 changed(k)=1
enddo

if ((sum_off_diag.le.EPS_12_4*max_S).and. &
    (sum_off_diag.ge.zero)) then
 state=0
 sanity_err=zero
else if (sum_off_diag.gt.zero) then
 ! do nothing
else
 print *,"sum_off_diag invalid: ",sum_off_diag
 print *,"max_S=",max_S
 stop
endif

do while (state.ne.0)
 m=1
 do k=2,n-1
  if ((ind(k).ge.1).and.(ind(m).ge.1).and. &
      (ind(k).le.n).and.(ind(m).le.n).and. &
      (k.ge.1).and.(k.le.n).and.(m.ge.1).and.(m.le.n)) then
   if (abs(S_MOD(k,ind(k))).gt.abs(S_MOD(m,ind(m)))) then
    m=k
   endif
  else
   print *,"k,m,ind(k),ind(m) bad ",k,m,ind(k),ind(m)
   stop
  endif
 enddo
 k=m
 l=ind(m)
 if ((l.ge.1).and.(l.le.n)) then
  ! do nothing
 else
  print *,"l invalid ",l
  print *,"m=",m
  stop
 endif

 p=S_MOD(k,l)
 y=(evals(l)-evals(k))/two
 d=abs(y)+sqrt(p**2+y**2)
 r=sqrt(p**2+d**2)
 if (r.gt.zero) then
  cosrot=d/r
  sinrot=p/r
 else
  print *,"expecting r>0 in fort_jacobi_eigenvalue"
  print *,"k,l,p,y,d,r ",k,l,p,y,d,r
  print *,"max_S=",max_S
  do i=1,n
  do j=1,n
   print *,"i,j,S(i,j) ",i,j,S(i,j)
  enddo
  enddo
  stop
 endif
 if (d.gt.zero) then
  t=p**2/d
 else
  print *,"expecting d>0 in fort_jacobi_eigenvalue"
  print *,"max_S=",max_S
  stop
 endif
 if (y.lt.zero) then
  sinrot=-sinrot
  t=-t
 else if (y.ge.zero) then
  ! do nothing
 else
  print *,"y invalid"
  stop
 endif
 S_MOD(k,l)=zero
 call EVAL_update(k,-t,y,changed,evals,state,n)       
 call EVAL_update(l,t,y,changed,evals,state,n)       
 do i=1,k-1
  call EVAL_rotate(i,k,i,l,S_MOD,n,sinrot,cosrot)
 enddo
 do i=k+1,l-1
  call EVAL_rotate(k,i,i,l,S_MOD,n,sinrot,cosrot)
 enddo
 do i=l+1,n
  call EVAL_rotate(k,i,l,i,S_MOD,n,sinrot,cosrot)
 enddo
 do i=1,n
  eik=evecs(i,k)
  eil=evecs(i,l)
  evecs(i,k)=cosrot*eik-sinrot*eil
  evecs(i,l)=sinrot*eik+cosrot*eil
 enddo
 call maxind(k,S_MOD,n,ind(k))
 call maxind(l,S_MOD,n,ind(l))

  ! AX=X Lambda
  ! A=X Lambda X^{-1}=X Lambda X^T
 do i=1,n
 do j=1,n  
  XL(i,j)=evecs(i,j)*evals(j)
 enddo
 enddo

 do i=1,n
 do j=1,n  
  XLXT(i,j)=zero
  do k_in=1,n
   XLXT(i,j)=XLXT(i,j)+XL(i,k_in)*evecs(j,k_in)
  enddo
 enddo
 enddo

 sanity_err=zero
 do i=1,n
 do j=1,n  
  if (abs(XLXT(i,j)-S(i,j)).gt.sanity_err) then
   sanity_err=abs(XLXT(i,j)-S(i,j))
  endif
 enddo
 enddo
 if (sanity_err.ge.EPS_10_4*max_S) then
  ! do nothing
 else if (sanity_err.le.EPS_10_4*max_S) then
  state=0
 else
  print *,"sanity_err became corrupt: ",sanity_err
  stop
 endif

enddo !do while (state.ne.0)

if (sanity_err.ge.EPS_8_3*max_S) then
 print *,"sanity_err too large(1): ",sanity_err,max_S
 stop
else if (sanity_err.le.EPS_8_3*max_S) then
 ! do nothing
else
 print *,"sanity_err is NaN: ",sanity_err,max_S
 stop
endif

 ! restore S_MOD
do k=1,n-1
 do l=k+1,n
  S_MOD(k,l)=S_MOD(l,k)
 enddo
enddo
do i=1,n
do j=1,n
 if (abs(S_MOD(i,j)-S(i,j)).le.EPS_10_4*max_S) then
  ! do nothing
 else
  print *,"S_MOD not properly restored"
  print *,"i,j,n,S_MOD(i,j),S(i,j),abs(S_MOD-S): ", &
    i,j,n,S_MOD(i,j),S(i,j),abs(S_MOD(i,j)-S(i,j))
  stop
 endif
enddo
enddo

do k=1,n-1
 m=k 
 do l=k+1,n
  if (abs(evals(l)).gt.abs(evals(m))) then
   m=l
  else if (abs(evals(l)).le.abs(evals(m))) then
   ! do nothing
  else
   print *,"evals NaN error"
   stop
  endif
 enddo
 if (k.ne.m) then
  swap_hold=evals(m)
  evals(m)=evals(k)
  evals(k)=swap_hold
  do i=1,n
   swap_hold=evecs(i,m)
   evecs(i,m)=evecs(i,k)
   evecs(i,k)=swap_hold
  enddo
 endif
enddo

 ! AX=X Lambda
 ! A=X Lambda X^{-1}=X Lambda X^T
do i=1,n
do j=1,n  
 XL(i,j)=evecs(i,j)*evals(j)
enddo
enddo

do i=1,n
do j=1,n  
 XLXT(i,j)=zero
 do k_in=1,n
  XLXT(i,j)=XLXT(i,j)+XL(i,k_in)*evecs(j,k_in)
 enddo
enddo
enddo

sanity_err=zero
do i=1,n
do j=1,n  
 if (abs(XLXT(i,j)-S(i,j)).gt.sanity_err) then
  sanity_err=abs(XLXT(i,j)-S(i,j))
 endif
enddo
enddo

if (sanity_err.ge.EPS_8_3*max_S) then
 print *,"sanity_err too large(2): ",sanity_err,max_S
 stop
else if (sanity_err.le.EPS_8_3*max_S) then
 ! do nothing
else
 print *,"sanity_err became corrupt: ",sanity_err,max_S
 stop
endif

if (1.eq.0) then
 print *,"end of jacobi eval"
endif

return
end subroutine fort_jacobi_eigenvalue

subroutine abs_value_determinant(S,n,determinant_out)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: n
real(amrex_real), INTENT(in) :: S(n,n)
real(amrex_real), INTENT(out) :: determinant_out

real(amrex_real) :: S_local(n,n)
real(amrex_real) :: STS(n,n)
real(amrex_real) :: evals_S(n)
real(amrex_real) :: evecs_S(n,n)
real(amrex_real) :: evals_STS(n)
real(amrex_real) :: evecs_STS(n,n)
integer :: i,j,k
real(amrex_real) :: max_eval_sqr

if (n.ge.2) then
 ! do nothing
else
 print *,"expecting n>=2(abs_value_determinant)"
 stop
endif

do i=1,n
do j=1,n
 S_local(i,j)=S(i,j)
 STS(i,j)=zero
 do k=1,n
  STS(i,j)=STS(i,j)+S(k,i)*S(k,j)
 enddo
enddo
enddo

call fort_jacobi_eigenvalue(S_local,evals_S,evecs_S,n)
call fort_jacobi_eigenvalue(STS,evals_STS,evecs_STS,n)

max_eval_sqr=-1.0D+20

do i=1,n

 if (evals_STS(i).lt.-EPS_10_5) then
  print *,"evals_STS(i) cannot be negative(abs_value_determinant): ", &
     i,n,evals_STS(i)
 else if (evals_STS(i).lt.zero) then
  evals_STS(i)=zero
 else if (evals_STS(i).ge.zero) then
  !do nothing
 else
  print *,"evals_STS(i) bust: ",i,n,evals_STS(i)
  stop
 endif

 if (evals_STS(i).gt.max_eval_sqr) then
  max_eval_sqr=evals_STS(i)
 else if (evals_STS(i).le.max_eval_sqr) then
  ! do nothing
 else
  print *,"evals_STS or max_eval_sqr invalid(abs_value_determinant): ", &
          i,n,evals_STS(i), &
          max_eval_sqr
  stop
 endif

 if (evals_STS(i).lt.zero) then
  print *,"evals_STS(i) cannot be negative(abs_value_determinant): ", &
     i,n,evals_STS(i)
  stop
 else if (evals_STS(i).ge.zero) then
  ! do nothing
 else
  print *,"evals_STS(i) is NaN(abs_value_determinant): ",i,n,evals_STS(i)
  stop
 endif

enddo ! i=1,n

if (max_eval_sqr.lt.zero) then
 print *,"max_eval_sqr cannot be negative: ",max_eval_sqr
 stop
else if (max_eval_sqr.ge.zero) then
 ! do nothing
else
 print *,"max_eval_sqr is NaN: ",max_eval_sqr
 stop
endif

max_eval_sqr=max(max_eval_sqr,one)
do i=1,n
 if ((abs(evals_S(i)**2-evals_STS(i)).le.EPS_8_4*max_eval_sqr).or. &
     (1.eq.1)) then
  ! do nothing
 else
  print *,"evals_S and evals_STS inconsistent"
  print *,"max_eval_sqr= ",max_eval_sqr
  print *,"i,n = ",i,n
  print *,"evals_S(i)= ",evals_S(i)
  print *,"evals_STS(i)= ",evals_STS(i)
  stop
 endif
enddo

determinant_out=one
do i=1,n
 determinant_out=determinant_out*evals_S(i)
enddo
determinant_out=abs(determinant_out)

end subroutine abs_value_determinant

subroutine project_to_traceless(S,n)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: n
real(amrex_real), INTENT(inout) :: S(n,n)

real(amrex_real) :: S_local(n,n)
real(amrex_real) :: STS(n,n)
real(amrex_real) :: evals_project(n)
real(amrex_real) :: XL(n,n)
real(amrex_real) :: evals_S(n)
real(amrex_real) :: evecs_S(n,n)
real(amrex_real) :: evals_STS(n)
real(amrex_real) :: evecs_STS(n,n)
integer :: i,j,k
real(amrex_real) :: max_eval_sqr
real(amrex_real) :: trace_S

if (n.ge.2) then
 ! do nothing
else
 print *,"expecting n>=2 (project_to_traceless)"
 stop
endif

do i=1,n
do j=1,n
 S_local(i,j)=S(i,j)
 STS(i,j)=zero
 do k=1,n
  STS(i,j)=STS(i,j)+S(k,i)*S(k,j)
 enddo
enddo
enddo

call fort_jacobi_eigenvalue(S_local,evals_S,evecs_S,n)
call fort_jacobi_eigenvalue(STS,evals_STS,evecs_STS,n)

max_eval_sqr=-1.0D+20
do i=1,n

 if (evals_STS(i).lt.-EPS_10_5) then
  print *,"evals_STS(i) cannot be negative(project_to_traceless): ", &
     i,n,evals_STS(i)
 else if (evals_STS(i).lt.zero) then
  evals_STS(i)=zero
 else if (evals_STS(i).ge.zero) then
  !do nothing
 else
  print *,"evals_STS(i) bust: ",i,n,evals_STS(i)
  stop
 endif

 if (evals_STS(i).gt.max_eval_sqr) then
  max_eval_sqr=evals_STS(i)
 else if (evals_STS(i).le.max_eval_sqr) then
  ! do nothing
 else
  print *,"evals_STS or max_eval_sqr invalid(project_to_traceless): ", &
          i,n,evals_STS(i), &
          max_eval_sqr
  stop
 endif

 if (evals_STS(i).lt.zero) then
  print *,"evals_STS(i) cannot be negative(project_to_Traceless): ", &
     i,n,evals_STS(i)
  stop
 else if (evals_STS(i).ge.zero) then
  ! do nothing
 else
  print *,"evals_STS(i) is NaN(project_to_traceless): ",i,n,evals_STS(i)
  stop
 endif

enddo !i=1,n

if (max_eval_sqr.lt.zero) then
 print *,"max_eval_sqr cannot be negative: ",max_eval_sqr
 stop
else if (max_eval_sqr.ge.zero) then
 ! do nothing
else
 print *,"max_eval_sqr is NaN: ",max_eval_sqr
 stop
endif

max_eval_sqr=max(max_eval_sqr,one)
do i=1,n
 if ((abs(evals_S(i)**2-evals_STS(i)).le.EPS_8_4*max_eval_sqr).or. &
     (1.eq.1)) then
  ! do nothing
 else
  print *,"evals_S and evals_STS inconsistent"
  print *,"max_eval_sqr= ",max_eval_sqr
  print *,"i,n = ",i,n
  print *,"evals_S(i)= ",evals_S(i)
  print *,"evals_STS(i)= ",evals_STS(i)
  stop
 endif
enddo

trace_S=zero
do i=1,n
 trace_S=trace_S+evals_S(i)
enddo

do i=1,n
 evals_project(i)=evals_S(i)-trace_S/n
enddo
 ! AX=X Lambda
 ! A=X Lambda X^T
do i=1,n
do j=1,n
 XL(i,j)=evecs_STS(i,j)*evals_project(j)
enddo
enddo
do i=1,n
do j=1,n
 S(i,j)=zero
 do k=1,n
  S(i,j)=S(i,j)+XL(i,k)*evecs_STS(j,k)
 enddo
enddo
enddo

end subroutine project_to_traceless



subroutine project_to_positive_definite(S,n,max_condition_number,unity_det)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in)    :: max_condition_number
integer, INTENT(in) :: unity_det
integer, INTENT(in) :: n
real(amrex_real), INTENT(inout) :: S(n,n)

real(amrex_real) :: S_local(n,n)
real(amrex_real) :: STS(n,n)
real(amrex_real) :: evals_project(n)
real(amrex_real) :: XL(n,n)
real(amrex_real) :: evals_S(n)
real(amrex_real) :: evecs_S(n,n)
real(amrex_real) :: evals_STS(n)
real(amrex_real) :: evecs_STS(n,n)
integer :: i,j,k
real(amrex_real) :: min_eval
real(amrex_real) :: max_eval_sqr
real(amrex_real) :: product_evals

if (max_condition_number.gt.one) then
 ! do nothing
else
 print *,"in project_to_positive_definite"
 print *,"max_condition_number invalid: ",max_condition_number
 stop
endif

if ((n.eq.2).or.(n.eq.3)) then
 ! do nothing
else
 print *,"expecting n==2 or 3: ",n
 stop
endif

do i=1,n
do j=1,n
 S_local(i,j)=S(i,j)
 STS(i,j)=zero
 do k=1,n
  STS(i,j)=STS(i,j)+S(k,i)*S(k,j)
 enddo
enddo
enddo

call fort_jacobi_eigenvalue(S_local,evals_S,evecs_S,n)
call fort_jacobi_eigenvalue(STS,evals_STS,evecs_STS,n)

max_eval_sqr=-1.0D+20
do i=1,n

 if (evals_STS(i).lt.-EPS_10_5) then
  print *,"evals_STS(i) cannot be negative(project_to_positive_def): ", &
     i,n,evals_STS(i)
 else if (evals_STS(i).lt.zero) then
  evals_STS(i)=zero
 else if (evals_STS(i).ge.zero) then
  !do nothing
 else
  print *,"evals_STS(i) bust: ",i,n,evals_STS(i)
  stop
 endif

 if (evals_STS(i).gt.max_eval_sqr) then
  max_eval_sqr=evals_STS(i)
 else if (evals_STS(i).le.max_eval_sqr) then
  ! do nothing
 else
  print *,"evals_STS or max_eval_sqr invalid(project_to_positive_definite):", &
    i,n,evals_STS(i),max_eval_sqr
  stop
 endif

 if (evals_STS(i).lt.zero) then
  print *,"evals_STS(i) cannot be negative(project_to_positive_definite): ", &
     i,n,evals_STS(i)
  stop
 else if (evals_STS(i).ge.zero) then
  ! do nothing
 else
  print *,"evals_STS(i) is NaN(project_to_pos_definite): ",i,n,evals_STS(i)
  stop
 endif

enddo !i=1,n

if (max_eval_sqr.lt.zero) then
 print *,"max_eval_sqr cannot be negative: ",max_eval_sqr
 stop
else if (max_eval_sqr.ge.zero) then
 ! do nothing
else
 print *,"max_eval_sqr is NaN: ",max_eval_sqr
 stop
endif

do i=1,n
 if (abs(evals_S(i)**2-evals_STS(i)).le.max_eval_sqr) then
  ! do nothing
 else
  print *,"evals_S and evals_STS inconsistent"
  print *,"max_eval_sqr= ",max_eval_sqr
  print *,"i,n = ",i,n
  print *,"evals_S(i)= ",evals_S(i)
  print *,"evals_STS(i)= ",evals_STS(i)
  stop
 endif
enddo !i=1,n

product_evals=one
min_eval=sqrt(max_eval_sqr)/max_condition_number
do i=1,n
 evals_project(i)=max(min_eval,evals_S(i))
 product_evals=product_evals*evals_project(i)
enddo

if (product_evals.gt.zero) then

 if (unity_det.eq.1) then
  product_evals=product_evals**(one/n)
  if (product_evals.gt.zero) then
   do i=1,n
    evals_project(i)=evals_project(i)/product_evals
   enddo
  else
   print *,"product_evals invalid: ",product_evals
   stop
  endif
 else if (unity_det.eq.0) then
  ! do nothing
 else
  print *,"unity_det invalid: ",unity_det
  stop
 endif

else
 print *,"product_evals invalid: ",product_evals
 stop
endif

 ! AX=X Lambda
 ! A=X Lambda X^T
do i=1,n
do j=1,n
 XL(i,j)=evecs_STS(i,j)*evals_project(j)
enddo
enddo
do i=1,n
do j=1,n
 S(i,j)=zero
 do k=1,n
  S(i,j)=S(i,j)+XL(i,k)*evecs_STS(j,k)
 enddo
enddo
enddo

end subroutine project_to_positive_definite

 ! called from point_updatetensor which is declared in GLOBALUTIL.F90
 ! point_updatetensor called from fort_updatetensor 
 ! which is declared in GODUNOV_3D.F90.
 ! A=Q+I must be symmetric and positive definite.
 ! polymer_factor = 1/L
subroutine project_A_to_positive_definite_or_traceless(A, &
   viscoelastic_model,polymer_factor,unity_det)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(inout) :: A(3,3)
integer, INTENT(in) :: viscoelastic_model
real(amrex_real), INTENT(in) :: polymer_factor
integer, INTENT(in) :: unity_det
integer A_dim
integer i,j
real(amrex_real), parameter :: max_condition_number=1000.0
real(amrex_real) local_diag
real(amrex_real), dimension(:,:), allocatable :: A_local

if ((polymer_factor.ge.zero).and. &
    (polymer_factor.lt.one)) then
 !do nothing
else
 print *,"polymer_factor invalid: ",polymer_factor
 stop
endif

if (SDIM.eq.2) then
 if (levelrz.eq.COORDSYS_RZ) then
  A_dim=2
 else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  A_dim=3
 else if (levelrz.eq.COORDSYS_CARTESIAN) then
  A_dim=2
 else
  print *,"levelrz invalid: ",levelrz
  stop
 endif
else if (SDIM.eq.3) then
 if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  A_dim=3
 else if (levelrz.eq.COORDSYS_CARTESIAN) then
  A_dim=3
 else
  print *,"levelrz invalid: ",levelrz
  stop
 endif
else
 print *,"dimension bust"
 stop
endif

allocate(A_local(A_dim,A_dim))

if (A_dim.eq.3) then
 ! do nothing
else if (A_dim.eq.2) then
 do j=1,3
 do i=1,3
  if ((i.eq.3).or.(j.eq.3)) then
   if (i.eq.j) then
    if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental
     !do nothing
    else
     if (A(i,j).gt.zero) then
      ! do nothing
     else
      print *,"A(i,j) failed sanity check"
      print *,"i,j,A(i,j) ",i,j,A(i,j)
      print *,"viscoelastic_model=",viscoelastic_model
      stop
     endif
    endif
   else if (i.ne.j) then
    if (A(i,j).eq.zero) then
     ! do nothing
    else
     print *,"A(i,j) failed sanity check"
     print *,"i,j,A(i,j) ",i,j,A(i,j)
     print *,"viscoelastic_model=",viscoelastic_model
     stop
    endif
   else
    print *,"i,j invalid: ",i,j
    stop
   endif
  else if ((i.ne.3).and.(j.ne.3)) then
   ! do nothing
  else
   print *,"i,j bust"
   stop
  endif
 enddo !i=1,3
 enddo !j=1,3

else
 print *,"A_dim invalid: ",A_dim
 stop
endif

do j=1,A_dim
do i=1,A_dim
 A_local(i,j)=A(i,j)
enddo
enddo
 
if ((viscoelastic_model.eq.NN_FENE_CR).or. & !FENE-CR
    (viscoelastic_model.eq.NN_OLDROYD_B).or. & !OLDROYD-B
    (viscoelastic_model.eq.NN_FENE_P).or. & !FENE-P
    (viscoelastic_model.eq.NN_LINEAR_PTT)) then !linear PTT

 call project_to_positive_definite(A_local,A_dim,max_condition_number,unity_det)

else if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental
 ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
 ! (grad U)_{ij}^{T}=U_{i,j}
 ! W_{ij}=(U_{j,i}-U_{i,j})/2
 ! Tran and Udaykumar: Omega_{ij}=(U_{i,j}-U_{j,i})/2=-W_{ij}
 ! DQ/Dt=2(D0-D^P)+QW-WQ   Q=S/mu Q=zero matrix at t=0 W=(grad U-grad U^T)/2
 ! Q^n+1 = (I+dt W)Q^{*}(I+dt W)^T + dt * 2(D0-D^P)  trace(Q)=0
 ! trace(Q)=sum lambda(Q)  lambda(Q)=eigenvalues of Q
 ! Q is traceless if trace(Q)=0 at t=0.
 if (1.eq.0) then
  call project_to_traceless(A_local,A_dim)
 endif
else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental Neo-Hookean
 ! Xia, Lu, Tryggvason 2018 Rapid Prototyping Journal
 ! Seungwon Shin, Jalel Chergui, Damir Juric
 ! f=dX/dx=F^{-1}  (dx/dX)_ij = (x_{i})_{j}
 ! F=dx/dX
 ! D X/Dt = 0
 ! Df/Dt + f grad U=0  
 ! Df^T/Dt + grad U^T f^T=0  
 ! Left Cauchy Green tensor B=F F^T=(f^T f)^{-1}
 ! D(f^T f)/Dt=
 !   f^T Df/Dt + Df^T/Dt f =
 !   f^T(-f grad U)+(-grad U^T f^T)f  
 ! let Binv=f^T f
 ! D Binv/Dt + Binv grad U + grad U^T Binv = 0
 ! D (Binv B)/Dt=
 !   D Binv/Dt B + Binv DB/Dt=
 !   (-Binv grad U - grad U^T Binv)B + Binv DB/Dt = 0
 ! -(grad U)B-B grad U^T + DB/Dt = 0
 ! DB/Dt = (grad U)B + B(grad U)^T
 ! equilibrium is B=I
 ! discretely, B should maintain as positive definite:
 ! B^n+1 = (I+dt grad U)Bstar(I+dt grad U)^T
 !
 call project_to_positive_definite(A_local,A_dim,max_condition_number,unity_det)
else
 print *,"viscoelastic_model invalid: ",viscoelastic_model
 stop
endif

do j=1,A_dim
do i=1,A_dim
 A(i,j)=A_local(i,j)
enddo
enddo

if (A_dim.eq.2) then

 local_diag=A(3,3)

 if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental
  !do nothing
 else
  if (local_diag.gt.zero) then
   !do nothing
  else
   print *,"local_diag must be positive: ",local_diag
   print *,"viscoelastic_model=",viscoelastic_model
   stop
  endif
 endif

 if (unity_det.eq.1) then
  if (local_diag.gt.zero) then
   if (levelrz.eq.COORDSYS_RZ) then
    do j=1,A_dim
    do i=1,A_dim
     A(i,j)=A(i,j)/sqrt(local_diag)
    enddo
    enddo
   else if (levelrz.eq.COORDSYS_CARTESIAN) then
    if (local_diag.eq.one) then
     ! do nothing
    else
     print *,"local_diag invalid: ",local_diag
     stop
    endif
   else
    print *,"levelrz invalid: ",levelrz
    stop
   endif
  else
   print *,"local_diag invalid: ",local_diag
   stop
  endif
 else if (unity_det.eq.0) then
  ! do nothing
 else
  print *,"unity_det invalid: ",unity_det
  stop
 endif

else if (A_dim.eq.3) then
 !do nothing
else
 print *,"A_dim invalid: ",A_dim
 stop
endif

deallocate(A_local)

return
end subroutine project_A_to_positive_definite_or_traceless

subroutine matrix_solve(AA,xx,bb,matstatus,numelem)
use probcommon_module
IMPLICIT NONE
integer numelem
real(amrex_real) AA(numelem,numelem)
real(amrex_real) xx(numelem)
real(amrex_real) bb(numelem)
real(amrex_real) alpha,holdvalue
integer i,j,k,holdj,matstatus
real(amrex_real) rowsum,maxrowsum

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
   if (abs(AA(i,i)).lt.EPS_10_5) then
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
    if (abs(AA(i,i)).lt.EPS_10_5) then
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
integer numelem
real(amrex_real) AA(numelem,numelem)

integer i

do i=1,numelem
 print *,AA(i,:)
 print *,"endrow"
enddo

return
end subroutine print_matrix

subroutine matrix_inverse(AA,xx,matstatus,numelem)
use probcommon_module
IMPLICIT NONE
integer, INTENT(in) :: numelem
real(amrex_real), INTENT(inout) :: AA(numelem,numelem)
real(amrex_real) AAhold(numelem,numelem)
real(amrex_real), INTENT(out) :: xx(numelem,numelem)
real(amrex_real) bb(numelem,numelem)
real(amrex_real) alpha,holdvalue
integer i,j,k,holdj
integer, INTENT(out) :: matstatus

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
   if (abs(AA(i,i)).lt.EPS30) then
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
     if (abs(AA(i,i)).lt.EPS30) then
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
    if (abs(holdvalue).gt.EPS12) then
     print *,"inverse failed1"
     print *,"AAhold="
     call print_matrix(AAhold,numelem)
     print *,"xx="
     call print_matrix(xx,numelem)
     stop
    endif
   else if (i.eq.j) then
    if (abs(holdvalue-one).gt.EPS12) then
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

subroutine project_to_tet(sdim,xtarget,xtet)
IMPLICIT NONE

integer, INTENT(in) :: sdim
real(amrex_real), INTENT(inout) :: xtarget(sdim)
real(amrex_real), INTENT(in) :: xtet(sdim+1,sdim)
real(amrex_real) :: mapmat(sdim,sdim)
real(amrex_real) :: mapmat_inv(sdim,sdim)
real(amrex_real) :: mapmat_scratch(sdim,sdim)
real(amrex_real) :: xcomp(sdim)
real(amrex_real) :: xcomp_sum

integer i,j
integer matstatus

! xphys = A xcomp + x0
! xcomp=Ainv(xphys-x0)
! A=x1-x0 y1-y0 z1-z0
!   x2-x0 y2-y0 z2-z0
!   x3-x0 y3-y0 z3-z0
do i=1,sdim
do j=1,sdim
 mapmat(i,j)=xtet(i+1,j)-xtet(1,j)
 mapmat_inv(i,j)=mapmat(i,j)
 mapmat_scratch(i,j)=mapmat(i,j)
enddo
enddo

call matrix_inverse(mapmat_scratch,mapmat_inv,matstatus,sdim)

if (matstatus.eq.1) then
 do i=1,sdim
  xcomp(i)=zero
  do j=1,sdim
   xcomp(i)=xcomp(i)+mapmat_inv(i,j)*(xtarget(j)-xtet(1,j))
  enddo
 enddo

 do j=1,sdim
  if (xcomp(j).gt.one) then
   xcomp(j)=one
  endif
  if (xcomp(j).lt.zero) then
   xcomp(j)=zero
  endif
 enddo

 xcomp_sum=zero
 do j=1,sdim
  xcomp_sum=xcomp_sum+xcomp(j)
 enddo
 if (xcomp_sum.gt.one) then
  do j=1,sdim
   xcomp(j)=xcomp(j)+(one-xcomp_sum)/dble(sdim)
  enddo
 else if ((xcomp_sum.ge.zero).and.(xcomp_sum.le.one)) then
  ! do nothing
 else
  print *,"xcomp_sum invalid project_to_tet: ",xcomp_sum
  stop
 endif

 do j=1,sdim
  if (xcomp(j).gt.one) then
   xcomp(j)=one
  endif
  if (xcomp(j).lt.zero) then
   xcomp(j)=zero
  endif
 enddo

 do i=1,sdim
  xtarget(i)=zero
  do j=1,sdim
   xtarget(i)=xtarget(i)+mapmat(i,j)*xcomp(j)
  enddo
  xtarget(i)=xtarget(i)+xtet(1,i)
 enddo
else
 print *,"matstatus invalid project_to_tet: ",matstatus
 stop
endif

end subroutine project_to_tet

subroutine print_visual_descriptor(im,n_fortran)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: im,n_fortran
integer :: n_cpp

if ((im.ge.1).and.(im.le.num_materials)) then
 ! do nothing
else
 print *,"im invalid"
 stop
endif

print *,"im=1..num_materials; im,num_materials:",im,num_materials

n_cpp=n_fortran-1

if (n_cpp.eq.VISUALCOMP_X) then
 print *,"VISUAL X"
else if (n_cpp.eq.VISUALCOMP_Y) then
 print *,"VISUAL Y"
else if ((n_cpp.eq.VISUALCOMP_Z).and.(SDIM.eq.3)) then
 print *,"VISUAL Z"
else if (n_cpp.eq.VISUALCOMP_U) then
 print *,"VISUAL U"
else if (n_cpp.eq.VISUALCOMP_V) then
 print *,"VISUAL V"
else if ((n_cpp.eq.VISUALCOMP_W).and.(SDIM.eq.3)) then
 print *,"VISUAL W"
else if (n_cpp.eq.VISUALCOMP_PMG) then
 print *,"VISUAL PMG"
else if (n_cpp.eq.VISUALCOMP_DEN) then
 print *,"VISUAL DEN"
else if (n_cpp.eq.VISUALCOMP_TEMP) then
 print *,"VISUAL TEMP"
else if ((n_cpp.ge.VISUALCOMP_Y1).and. &
         (n_cpp.le.VISUALCOMP_Y1+num_species_var-1)) then
 print *,"VISUAL Y species component(1..num_species_var)=", &
         n_cpp-VISUALCOMP_Y1+1
else if (n_cpp.eq.VISUALCOMP_VORTMAG) then
 print *,"VISUAL VORTMAG"
else if ((n_cpp.ge.VISUALCOMP_LS).and. &
         (n_cpp.le.VISUALCOMP_LS+num_materials-1)) then
 print *,"VISUAL LS material component(1..num_materials)=",n_cpp-VISUALCOMP_LS+1
else if (n_cpp.eq.VISUALCOMP_MAGVEL) then
 print *,"VISUAL MAGVEL"
else
 print *,"n_cpp out of range"
 stop
endif

return
end subroutine print_visual_descriptor

      subroutine least_squares_interp(npoints,x0,xpos,wts,vals, &
       order,linearflag,coeffs,ncoeffs,sdim_parm)
      IMPLICIT NONE

      integer, INTENT(in) :: ncoeffs
      integer, INTENT(in) :: order
      integer, INTENT(in) :: npoints
      integer, INTENT(in) :: sdim_parm
      integer, INTENT(in) :: linearflag
      real(amrex_real), INTENT(in) :: x0(sdim_parm)
      real(amrex_real), INTENT(in) :: xpos(npoints,sdim_parm)
      real(amrex_real), INTENT(in) :: wts(npoints)
      real(amrex_real), INTENT(in) :: vals(npoints)
      real(amrex_real), INTENT(out) :: coeffs(ncoeffs)
      real(amrex_real) AA(ncoeffs,ncoeffs)
      real(amrex_real) BB(ncoeffs)
      integer ncoeffs_test
      integer orderx,ordery,orderz
      integer i,i2,ipoints,idx,dir
      integer ipower,jpower,kpower
      real(amrex_real) basisfn(ncoeffs)
      real(amrex_real) basisdir(0:order,sdim_parm)
      integer matstatus

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

      integer, INTENT(in) :: i,j,k,im
      integer, INTENT(in) :: nhalf
      integer :: dir
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(inout) :: LS
      real(amrex_real), INTENT(in), pointer :: dist(D_DECL(:,:,:),:)
      real(amrex_real) n(SDIM)
      real(amrex_real) nsave(SDIM)
      real(amrex_real) RR,mag

      if (nhalf.ne.3) then
       print *,"nhalf invalid"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid15"
       stop
      endif

      LS=dist(D_DECL(i,j,k),im)

      if ((levelrz.eq.COORDSYS_CARTESIAN).or.(levelrz.eq.COORDSYS_RZ)) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
 
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
       call prepare_normal(n,RR,mag,SDIM)
       if (mag.gt.zero) then
        RR=xsten(0,1)
        do dir=1,SDIM
         nsave(dir)=n(dir)
        enddo
        call prepare_normal(nsave,RR,mag,SDIM)
        if (mag.gt.zero) then
         mag=zero
         do dir=1,SDIM
          mag=mag+n(dir)*nsave(dir)
         enddo
         if (abs(mag).le.one+EPS2) then
          !do nothing
         else
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


       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      subroutine grid_type_to_box_type3D(grid_type,box_type)
      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(out) :: box_type(3)
      integer dir

      do dir=1,3
       box_type(dir)=0  ! default to CELL
      enddo
      if (grid_type.eq.-1) then
       ! do nothing
      else if ((grid_type.ge.0).and. &
               (grid_type.lt.3)) then
       box_type(grid_type+1)=1  ! NODE
      else if (grid_type.eq.3) then
       box_type(1)=1 ! NODE
       box_type(2)=1 ! NODE
      else if (grid_type.eq.4) then
       box_type(1)=1 ! NODE
       box_type(3)=1 ! NODE
      else if (grid_type.eq.5) then
       box_type(2)=1 ! NODE
       box_type(3)=1 ! NODE
      else
       print *,"grid_type invalid"
       stop
      endif
     
      return 
      end subroutine grid_type_to_box_type3D

      subroutine checkbound3D_array(lo,hi, &
      data_array, &
      ngrow,grid_type)
      IMPLICIT NONE

      integer, INTENT(in) :: lo(3), hi(3)
       ! INTENT(in) means the pointer cannot be reassigned.
       ! The data itself inherits the INTENT attribute from the
       ! target.
      real(amrex_real), INTENT(in), pointer :: data_array(:,:,:,:)
      integer, INTENT(in) :: ngrow,grid_type

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      integer box_type(3)

      integer hidata(4)
      integer lodata(4)
      integer dir2

      hidata=UBOUND(data_array)
      lodata=LBOUND(data_array)

      do dir2=1,3
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound 3d put breakpoint here"
        print *,"dir2=",dir2
        stop
       endif
       box_type(dir2)=0
      enddo
       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type3D(grid_type,box_type)

      do dir2=1,3
       if (lo(dir2).lt.0) then
        print *,"lo invalid in checkbound3D_array put breakpoint here"
        print *,"dir2,dataxlo ",dir2,lodata(dir2)
        print *,"dir2,dataxhi ",dir2,hidata(dir2)
        print *,"dir2,lo,ngrow ",dir2,lo(dir2),ngrow
        print *,"dir2,hi,ngrow ",dir2,hi(dir2),ngrow
        print *,"grid_type=",grid_type
        stop
       endif
      enddo

      if (ngrow.lt.0) then
       print *,"ngrow invalid in checkbound3D_array"
       stop
      endif

      do dir2=1,3

       if (lodata(dir2).gt.lo(dir2)-ngrow) then
        print *,"lo mismatch put breakpoint here"
        print *,"dir2=",dir2
        stop
       endif
       if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
        print *,"hi mismatch put breakpoint here1"
        print *,"dir2=",dir2
        stop
       endif

      enddo ! dir2=1..3

      return
      end subroutine checkbound3D_array

       ! grid_type=-1..5
      subroutine checkbound_array(lo,hi, &
       data_array, &
       ngrow,grid_type)
      IMPLICIT NONE

      integer, INTENT(in) ::  lo(SDIM), hi(SDIM)
        ! INTENT(in) means the pointer cannot be reassigned.
        ! The data itself inherits the INTENT attribute from the
        ! target.
      real(amrex_real), INTENT(in), pointer :: data_array(D_DECL(:,:,:),:)
      integer, INTENT(in) ::  ngrow
      integer, INTENT(in) ::  grid_type

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      integer box_type(SDIM)

      integer    hidata(SDIM+1)
      integer    lodata(SDIM+1)
      integer    dir2

      hidata=UBOUND(data_array)
      lodata=LBOUND(data_array)
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound_array put breakpoint here"
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
        print *,"lo invalid in checkbound_array put breakpoint here"
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

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound_array:lo mismatch put breakpoint here1"
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch put breakpoint here2 GLOBALUTIL.F90:6795"
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
      subroutine checkbound_array_integer(lo,hi, &
       data_array, &
       ngrow,grid_type)
      IMPLICIT NONE

      integer, INTENT(in) ::  lo(SDIM), hi(SDIM)
        ! INTENT(in) means the pointer cannot be reassigned.
        ! The data itself inherits the INTENT attribute from the
        ! target.
      integer, INTENT(in), pointer :: data_array(D_DECL(:,:,:),:)
      integer, INTENT(in) ::  ngrow
      integer, INTENT(in) ::  grid_type

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      integer box_type(SDIM)

      integer    hidata(SDIM+1)
      integer    lodata(SDIM+1)
      integer    dir2

      hidata=UBOUND(data_array)
      lodata=LBOUND(data_array)
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound_array put breakpoint here"
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
        print *,"lo invalid in checkbound_array put breakpoint here"
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

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound_array:lo mismatch put breakpoint here2"
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch put breakpoint here3"
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"box_type(dir2) ",box_type(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif

      enddo ! dir2=1..SDIM

      return
      end subroutine checkbound_array_integer



       ! grid_type=-1..5
      subroutine checkbound_array1(lo,hi, &
       data_array1, &
       ngrow,grid_type)
      IMPLICIT NONE

      integer, INTENT(in) ::  lo(SDIM), hi(SDIM)
        ! INTENT(in) means the pointer cannot be reassigned.
        ! The data itself inherits the INTENT attribute from the
        ! target.
      real(amrex_real), INTENT(in), pointer :: data_array1(D_DECL(:,:,:))
      integer, INTENT(in) ::  ngrow
      integer, INTENT(in) ::  grid_type

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      integer box_type(SDIM)

      integer    hidata(SDIM)
      integer    lodata(SDIM)
      integer    dir2

      hidata=UBOUND(data_array1)
      lodata=LBOUND(data_array1)
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound_array1 put breakpoint here"
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
        print *,"lo invalid in checkbound_array1 put breakpoint here"
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

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound_array:lo mismatch put breakpoint here3"
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch put breakpoint here4 6851"
         print *,"break GLOBALUTIL.F90:6851"
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
       ngrow,grid_type)
      IMPLICIT NONE

      integer, INTENT(in) ::  lo(SDIM), hi(SDIM)
        ! INTENT(in) means the pointer cannot be reassigned.
        ! The data itself inherits the INTENT attribute from the
        ! target.
      integer, INTENT(in), pointer :: data_array1(D_DECL(:,:,:))
      integer, INTENT(in) ::  ngrow
      integer, INTENT(in) ::  grid_type

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      integer box_type(SDIM)

      integer    hidata(SDIM)
      integer    lodata(SDIM)
      integer    dir2

      hidata=UBOUND(data_array1)
      lodata=LBOUND(data_array1)
 
      do dir2=1,SDIM
       if (lodata(dir2).gt.hidata(dir2)) then
        print *,"swapped bounds in checkbound_array1 put breakpoint here"
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
        print *,"lo invalid in checkbound_array1 put breakpoint here"
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

      do dir2=1,SDIM

        if (lodata(dir2).gt.lo(dir2)-ngrow) then
         print *,"checkbound_array:lo mismatch put breakpoint here4"
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dataxlo ",lodata(dir2)
         print *,"dataxhi ",hidata(dir2)
         print *,"grid_type=",grid_type
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+box_type(dir2)) then
         print *,"hi mismatch put breakpoint here5"
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

      subroutine h_coeffXY(r1,r2,ri,coeffs)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: r1,r2,ri
      real(amrex_real), INTENT(out) :: coeffs(3)

      if (r2.gt.r1) then
   
        ! ((r2-ri)^3/3 - (r1-ri)^3/3)/(r2-r1)=
        ! ((r2-ri)^3/3 - (r1-ri)^3/3)/((r2-ri)-(r1-ri)) 
        ! A^3-B^3=(A-B)*(A^2 + B^2 + AB)=A^3 + AB^2+A^2 B-BA^2-AB^2-B^3
        ! (r1-ri)^2+(r2-ri)^2+(r1-ri)(r2-ri)
       coeffs(1)=(r1**2)/three+r1*r2/three-r1*ri+(r2**2)/three- &
                 r2*ri+(ri**2)

        !A^2-B^2=(A-B)(A+B)
       coeffs(2)=half*(r1+r2)-ri

       coeffs(3)=one

      else
       print *,"r1 or r2 invalid"
       print *,"r1,r2,ri ",r1,r2,ri
       stop
      endif

      return
      end subroutine h_coeffXY

      subroutine h_coeffXYZ(x1,x2,y1,y2,xi,yi,coeffs)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x1,x2,y1,y2,xi,yi
      real(amrex_real), INTENT(out) :: coeffs(9)
      integer :: i,j,row
      real(amrex_real) :: coeffs_x(3)
      real(amrex_real) :: coeffs_y(3)

      if ((x2.gt.x1).and.(y2.gt.y1)) then
   
       ! h(x,y)=sum_ij aij (x-xi)^i (y-yi)^j 

       call h_coeffXY(x1,x2,xi,coeffs_x)
       call h_coeffXY(y1,y2,yi,coeffs_y)

       do i=0,2
       do j=0,2
        row=3*i+j+1
        coeffs(row)=coeffs_x(3-i)*coeffs_y(3-j)
       enddo
       enddo

      else
       print *,"x1,x2,y1, or y2 invalid"
       print *,"x1,x2,y1,y2,xi,yi ",x1,x2,y1,y2,xi,yi
       stop
      endif

      return
      end subroutine h_coeffXYZ



      subroutine hvertical_coeffRZ(r1,r2,ri,coeffs)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: r1,r2,ri
      real(amrex_real), INTENT(out) :: coeffs(3)
      integer dir

      if ((r1.ge.zero).and.(r2.ge.zero).and.(ri.gt.zero).and. &
          (r2.gt.r1)) then
     
       coeffs(1)=three*(r1**3)+three*(r1**2)*r2-eight*(r1**2)*ri+ &
        three*r1*(r2**2)-eight*r1*r2*ri+six*r1*(ri**2)+ &
        three*(r2**3)-eight*(r2**2)*ri+six*r2*(ri**2)

       coeffs(2)=four*(r1**2)+four*r1*r2-six*r1*ri+four*(r2**2)- &
         six*r2*ri

       coeffs(3)=six*(r1+r2)

       do dir=1,3
        coeffs(dir)=coeffs(dir)/(six*(r1+r2))
       enddo

      else
       print *,"r1,r2, or ri invalid"
       print *,"r1,r2,ri ",r1,r2,ri
       stop
      endif

      return
      end subroutine hvertical_coeffRZ


      subroutine hhorizontal_coeffRZ(z1,z2,zi,coeffs)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: z1,z2,zi
      real(amrex_real), INTENT(out) :: coeffs(3)

      if (z2.gt.z1) then
     
       coeffs(1)=(z1**2)/three+z1*z2/three-z1*zi+(z2**2)/three- &
           z2*zi+(zi**2)

       coeffs(2)=z1/two+z2/two-zi

       coeffs(3)=one

      else
       print *,"z1 or z2,invalid"
       stop
      endif

      return
      end subroutine hhorizontal_coeffRZ


      subroutine analyze_heights( &
        htfunc_LS, & !intent(in)
        htfunc_VOF, & !intent(in)
        xsten, &
        nhalf, &
        itan,jtan, &
        curvHT_LS, & !intent(out)
        curvHT_VOF, & !intent(out)
        curvHT_choice, & !intent(out)
        normal_dir, &
        xcenter, &
        n1d, &
        overall_crossing_status, & !intent(in) (1=success 0=not found)
        vof_height_function)
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: overall_crossing_status !1=success 0=not found
      integer, INTENT(in) :: vof_height_function
      real(amrex_real), INTENT(in) :: htfunc_LS(-1:1,-1:1)
      real(amrex_real), INTENT(in) :: htfunc_VOF(-1:1,-1:1)
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: itan,jtan
      integer, INTENT(in) :: normal_dir
      real(amrex_real), INTENT(out) :: curvHT_LS
      real(amrex_real), INTENT(out) :: curvHT_VOF
      real(amrex_real), INTENT(out) :: curvHT_choice
      real(amrex_real), INTENT(in) :: xcenter(SDIM)
      real(amrex_real), INTENT(in) :: n1d

      real(amrex_real) hxR,hxL,hx,hxx
      real(amrex_real) hyR,hyL,hy,hyy
      real(amrex_real) hxy
      real(amrex_real) arclen,arclenx,arcleny,arclenr,g,gx,gy,gr
      real(amrex_real) RR

      integer rowx,rowy,rowxy
      integer interval_x,interval_y
      real(amrex_real) x1,x2,x1_3D,x2_3D,xnot,ynot

      integer num_coeffs

      real(amrex_real) RHS_2D(3)
      real(amrex_real) AA_2D(3,3)
      real(amrex_real) coeffs_2D(3)
      real(amrex_real) xx_2D(3)

      real(amrex_real) RHS_3D(9)
      real(amrex_real) AA_3D(9,9)
      real(amrex_real) coeffs_3D(9)
      real(amrex_real) xx_3D(9)

      integer matstatus
      real(amrex_real) dr
      real(amrex_real) h_of_z,hprime_of_z,hdprime_of_z
      real(amrex_real) hprime_of_r,hdprime_of_r
      integer dir2

      if (nhalf.eq.2*ngrow_distance+1) then
       ! do nothing
      else
       print *,"nhalf invalid"
       stop
      endif
      if ((n1d.eq.-one).or.(n1d.eq.one)) then
       ! do nothing
      else
       print *,"n1d invalid"
       stop
      endif
      if ((normal_dir.ge.1).and.(normal_dir.le.SDIM)) then
       ! do nothing
      else
       print *,"normal_dir invalid"
       stop
      endif

       ! dark material on the bottom
       ! phi=h(x,y)-z  z=h(x,y)   (normal_dir=2)
       ! n=grad phi/|grad phi| = (h_x,h_y,-1)/sqrt(h_x^2+h_y^2+1)
       ! div n = (n_x)_x + (n_y)_y

      hxR=(htfunc_LS(1,0)-htfunc_LS(0,0))/(xsten(2,itan)-xsten(0,itan))
      hxL=(htfunc_LS(0,0)-htfunc_LS(-1,0))/(xsten(0,itan)-xsten(-2,itan))
      hx=(htfunc_LS(1,0)-htfunc_LS(-1,0))/(xsten(2,itan)-xsten(-2,itan))
      hxx=(hxR-hxL)/(xsten(1,itan)-xsten(-1,itan))
      hyR=zero
      hyL=zero
      hy=zero
      hxy=zero
      hyy=zero

      if (SDIM.eq.3) then
       hyR=(htfunc_LS(0,1)-htfunc_LS(0,0))/(xsten(2,jtan)-xsten(0,jtan))
       hyL=(htfunc_LS(0,0)-htfunc_LS(0,-1))/(xsten(0,jtan)-xsten(-2,jtan))
       hy=(htfunc_LS(0,1)-htfunc_LS(0,-1))/(xsten(2,jtan)-xsten(-2,jtan))
       hxy=(htfunc_LS(1,1)-htfunc_LS(-1,1)- &
            htfunc_LS(1,-1)+htfunc_LS(-1,-1))/ &
           ( (xsten(2,itan)-xsten(-2,itan))* &
             (xsten(2,jtan)-xsten(-2,jtan)) )
       hyy=(hyR-hyL)/(xsten(1,jtan)-xsten(-1,jtan))
      endif
      arclen=one+hx**2+hy**2
      arclenx=two*hx*hxx+two*hy*hxy
      arcleny=two*hx*hxy+two*hy*hyy
      g=one/sqrt(arclen)
      gx=-half*(arclen**(-1.5d0))*arclenx
      gy=-half*(arclen**(-1.5d0))*arcleny

      if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! phi=h(x,y)-z  z=h(x,y)   (normal_dir=2)
        !arclen=1+hx^2+hy^2
        !g(x,y)=arclen^{-1/2}  
        !n=grad phi/|grad phi|=(hx g,hy g,-g)
        !div n=hxx g + hx gx +hyy g + hy gy - gz
        !gx=(-1/2)arclen^{-3/2}arclenx
       curvHT_LS=hxx*g+hx*gx+hyy*g+hy*gy
      else if (levelrz.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       if (normal_dir.eq.1) then
         ! phi=h(z)-r  
         ! arclen=1+hz^2 
         ! g(z)=arclen^(-1/2)
         ! n=(-g,hz g) 
         ! div n =-(r g)_r/r + (hz g)_z=-g/r+hzz g + gz hz
        RR=htfunc_LS(0,0)
        if (RR.gt.zero) then
         ! do nothing
        else
         print *,"RR invalid"
         stop
        endif
        curvHT_LS=-g/RR+hxx*g+gx*hx
       else if (normal_dir.eq.2) then
         ! phi=h(r)-z
         ! arclen=1+hr^2
         ! g(r)=arclen^(-1/2)
         ! n=(hr g,-g)
         ! div n = (r hr g)_r/r = hr g/r+hrr g + hr gr
        RR=xcenter(1) 
        curvHT_LS=hx*g/RR+hxx*g+hx*gx
       else
        print *,"normal_dir invalid"
        stop
       endif
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
       if (normal_dir.eq.1) then
         ! phi=h(y,z)-r
         ! grad phi=(-1,hy/r,hz)
         ! arclen=1+(hy/r)^2+hz^2
         ! g(r,y,z)=arclen^(-1/2)
         ! n=(-g,hy g/r,hz g)
         ! div n=(-g r)_r/r+(hy g/r)_y/r+(hz g)_z=
         ! -g_r-g/r+(hyy g + gy hy)/r^2+hzz g + hz gz
         ! g_r=-1/2 arclen^(-3/2) arclen_r
         ! arclen_r=hy^2 (-2/r^3)
         ! g_y=-1/2 arclen^(-3/2) arclen_y
         ! arclen_y=2hy hyy/r^2 +2hz hzy
         ! g_z=-1/2 arclen^(-3/2) arclen_z
         ! arclen_z=2hy hyz/r^2 +2hz hzz
        RR=htfunc_LS(0,0)
        arclen=one+(hx/RR)**2+hy**2
        g=one/sqrt(arclen)
        arclenr=-two*(hx**2)/(RR**3)
        gr=-half*(arclen**(-1.5))*arclenr
        arclenx=two*hx*hxx/(RR**2)+two*hy*hxy
        arcleny=two*hx*hxy/(RR**2)+two*hy*hyy
        gx=-half*(arclen**(-1.5))*arclenx
        gy=-half*(arclen**(-1.5))*arcleny
        curvHT_LS=-gr-g/RR+(hxx*g+hx*gx)/(RR**2)+hyy*g+hy*gy
       else if (normal_dir.eq.2) then 
         ! phi=h(r,z)-y
         ! grad phi=(hr,-1/r,hz)
         ! arclen=hr^2+(1/r)^2+hz^2
         ! g(r,z)=arclen^(-1/2)
         ! n=(g hr,-g/r,hz g)
         ! div n=(g hr r)_r/r-(g/r)_y/r+(hz g)_z=
         !  g_r hr + g hrr + g hr/r +hzz g + hz gz=
         ! g_r=-1/2 arclen^(-3/2) arclen_r
         ! arclen_r=2 hr hrr -2/r^3 + 2hz hzr
         ! g_z=-1/2 arclen^(-3/2) arclen_z
         ! arclen_z=2 hr hrz +2hz hzz
        RR=xcenter(1) 
        arclen=hx**2+hy**2+(one/RR)**2
        g=one/sqrt(arclen)
        arclenx=two*hx*hxx-two/(RR**3)+two*hy*hxy
        arcleny=two*hx*hxy+two*hy*hyy
        gx=-half*(arclen**(-1.5))*arclenx
        gy=-half*(arclen**(-1.5))*arcleny
        curvHT_LS=gx*hx+g*hxx+g*hx/RR+hyy*g+hy*gy
       else if ((normal_dir.eq.3).and.(SDIM.eq.3)) then
         ! phi=h(r,y)-z
         ! grad phi=(hr,hy/r,-1)
         ! arclen=hr^2+(hy/r)^2+1
         ! g(r,y)=arclen^(-1/2)
         ! n=(g hr,g hy/r,-g)
         ! div n=(g hr r)_r/r+(g hy/r)_y/r=
         !  g_r h_r + g h_rr + g hr/r +(gy hy+g hyy)/r^2
         ! g_r=-1/2 arclen^(-3/2) arclen_r
         ! arclen_r=2 hr hrr + 2(hy/r)(hry r-hy)/r^2 =
         !          2 hr hrr + 2 hy hry/r^2 - 2 hy^2/r^3
         ! g_y=-1/2 arclen^(-3/2) arclen_y
         ! arclen_y=2 hr hry +2hy hyy/r^2
        RR=xcenter(1) 
        arclen=hx**2+(hy/RR)**2+one
        g=one/sqrt(arclen)
        arclenx=two*hx*hxx+two*hy*hxy/(RR**2)- &
                two*(hy**2)/(RR**3)
        arcleny=two*hx*hxy+two*hy*hyy/(RR**2)
        gx=-half*(arclen**(-1.5))*arclenx
        gy=-half*(arclen**(-1.5))*arcleny
        curvHT_LS=gx*hx+g*hxx+g*hx/RR+(hyy*g+hy*gy)/(RR**2)
       else
        print *,"normal_dir invalid"
        stop
       endif
      else
       print *,"levelrz invalid init height ls 5"
       stop
      endif

      curvHT_VOF=curvHT_LS
      curvHT_choice=curvHT_LS

      if (overall_crossing_status.eq.0) then
       ! do nothing
      else if (overall_crossing_status.eq.1) then

       if (vof_height_function.eq.1) then

        if (levelrz.eq.COORDSYS_RZ) then

         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif

         num_coeffs=3

         interval_x=-3

         do rowx=1,num_coeffs

          x1=xsten(interval_x,itan)
          x2=xsten(interval_x+2,itan)
          xnot=xsten(0,itan)
          RHS_2D(rowx)=htfunc_VOF(rowx-2,0)

          if (normal_dir.eq.1) then  ! horizontal column in the r direction
           call hhorizontal_coeffRZ(x1,x2,xnot,coeffs_2D)
           RHS_2D(rowx)=htfunc_VOF(rowx-2,0)**2
          else if (normal_dir.eq.2) then  ! vertical column in the z direction
           if (xnot.gt.zero) then
            if (x1.ge.zero) then
             call hvertical_coeffRZ(x1,x2,xnot,coeffs_2D)
            else if (x1.lt.zero) then
             RHS_2D(rowx)=zero
             coeffs_2D(1)=-two*xnot
             coeffs_2D(2)=one
             coeffs_2D(3)=zero
            else
             print *,"x1 invalid"
             stop
            endif
           else
            print *,"something wrong, xnot<=0 ",xnot
            stop
           endif
          else
           print *,"normal_dir invalid"
           stop
          endif
  
          do dir2=1,num_coeffs
           AA_2D(rowx,dir2)=coeffs_2D(dir2)
          enddo
          interval_x=interval_x+2
         enddo ! rowx=1..num_coeffs

           ! xx(1)*(x-x0)^2+xx(2)*(x-x0)+xx(3)
           !matstatus=1 => ok.
         call matrix_solve(AA_2D,xx_2D,RHS_2D,matstatus,num_coeffs) 

         if (matstatus.eq.1) then
          if (normal_dir.eq.1) then  ! horizontal column in the r direction
           dr=xsten(1,1)-xsten(-1,1)
           if (dr.gt.zero) then
            if (xx_2D(3).gt.zero) then
             h_of_z=sqrt(xx_2D(3))
             if (h_of_z.ge.half*dr) then
              hprime_of_z=xx_2D(2)/(two*h_of_z)
              hdprime_of_z=(xx_2D(1)-hprime_of_z**2)/h_of_z
              g=sqrt(one+hprime_of_z**2)
              curvHT_VOF=-(one/(h_of_z*g))+hdprime_of_z/(g**3)
             else if (h_of_z.gt.zero) then
              curvHT_VOF=curvHT_LS
             else
              print *,"h_of_z too small, h_of_z=",h_of_z
              stop
             endif
            else if (xx_2D(3).le.zero) then
             curvHT_VOF=curvHT_LS
            else
             print *,"xx_2D(3) bust: ",xx_2D(3)
             stop
            endif
           else
            print *,"dr invalid, dr=",dr
            stop
           endif
          else if (normal_dir.eq.2) then  ! vertial column in the z direction
           hprime_of_r=xx_2D(2)
           hdprime_of_r=two*xx_2D(1)
           g=sqrt(one+hprime_of_r**2)
           curvHT_VOF=(hprime_of_r/(xnot*g))+hdprime_of_r/(g**3)
          else
           print *,"normal_dir invalid"
           stop
          endif
         else if (matstatus.eq.0) then
          curvHT_VOF=curvHT_LS
         else
          print *,"matstatus corrupt for RZ curvature coeff"
          stop
         endif

        else if (levelrz.eq.COORDSYS_CYLINDRICAL) then

         print *,"VFRAC height function invalid:levelrz=COORDSYS_CYLINDRICAL"

        else if (levelrz.eq.COORDSYS_CARTESIAN) then

         if (SDIM.eq.2) then

          num_coeffs=3

          interval_x=-3

          do rowx=1,num_coeffs

           x1=xsten(interval_x,itan)
           x2=xsten(interval_x+2,itan)
           xnot=xsten(0,itan)
           RHS_2D(rowx)=htfunc_VOF(rowx-2,0)

           call h_coeffXY(x1,x2,xnot,coeffs_2D)
  
           do dir2=1,num_coeffs
            AA_2D(rowx,dir2)=coeffs_2D(dir2)
           enddo
           interval_x=interval_x+2
          enddo ! rowx=1..num_coeffs

           ! xx(1)*(x-x0)^2+xx(2)*(x-x0)+xx(3)
           !matstatus=1 => ok.
          call matrix_solve(AA_2D,xx_2D,RHS_2D,matstatus,num_coeffs) 

          if (matstatus.eq.1) then
           hprime_of_r=xx_2D(2)
           hdprime_of_r=two*xx_2D(1)
           g=sqrt(one+hprime_of_r**2)
           curvHT_VOF=hdprime_of_r/(g**3)
          else if (matstatus.eq.0) then
           curvHT_VOF=curvHT_LS
          else
           print *,"matstatus corrupt for XY curvature coeff"
           stop
          endif

         else if (SDIM.eq.3) then

           !h(x,y)=sum_{i,j=0..2}  aij(x-x0)^i(y-y0)^j 
           !in 3D, low order derivatives are first.
           ! flattening: row=3*i+j+1
          num_coeffs=9

          interval_x=-3

          rowxy=1
          do rowx=1,3

           interval_y=-3

           do rowy=1,3

            x1=xsten(interval_x,itan)
            x2=xsten(interval_x+2,itan)

            x1_3D=xsten(interval_y,jtan)
            x2_3D=xsten(interval_y+2,jtan)

            xnot=xsten(0,itan)
            ynot=xsten(0,jtan)
            if (rowxy.eq.(3*(rowx-1)+rowy)) then
             ! do nothing
            else
             print *,"rowxy fails flattening sanity check"
             stop
            endif
            RHS_3D(rowxy)=htfunc_VOF(rowx-2,rowy-2)

            call h_coeffXYZ(x1,x2,x1_3D,x2_3D,xnot,ynot,coeffs_3D)
  
            do dir2=1,num_coeffs
             AA_3D(rowxy,dir2)=coeffs_3D(dir2)
            enddo
            interval_y=interval_y+2
            rowxy=rowxy+1
           enddo ! rowy=1..3

           interval_x=interval_x+2

          enddo ! rowx=1..3

          if (rowxy.eq.num_coeffs+1) then
           ! do nothing
          else
           print *,"rowxy invalid"
           stop
          endif

           ! 2D: xx(1)*(x-x0)^2+xx(2)*(x-x0)+xx(3)
           ! 3D: h(x,y)=sum_{i,j=0..2}  aij(x-x0)^i(y-y0)^j 
           ! flattening: row=3*i+j+1
           !in 3D, low order derivatives are first.
           !matstatus=1 => ok.
          call matrix_solve(AA_3D,xx_3D,RHS_3D,matstatus,num_coeffs) 

          if (matstatus.eq.1) then
           rowx=1
           rowy=0
           rowxy=3*rowy+rowx+1 
           hx=xx_3D(rowxy)
           rowx=0
           rowy=1
           rowxy=3*rowy+rowx+1 
           hy=xx_3D(rowxy)
           rowx=2
           rowy=0
           rowxy=3*rowy+rowx+1 
           hxx=two*xx_3D(rowxy)
           rowx=0
           rowy=2
           rowxy=3*rowy+rowx+1 
           hyy=two*xx_3D(rowxy)
           rowx=1
           rowy=1
           rowxy=3*rowy+rowx+1 
           hxy=xx_3D(rowxy)

           arclen=one+hx**2+hy**2
           arclenx=two*hx*hxx+two*hy*hxy
           arcleny=two*hx*hxy+two*hy*hyy
           g=one/sqrt(arclen)
           gx=-half*(arclen**(-1.5d0))*arclenx
           gy=-half*(arclen**(-1.5d0))*arcleny

           ! phi=h(x,y)-z  z=h(x,y)   (normal_dir=2)
           !arclen=1+hx^2+hy^2
           !g(x,y)=arclen^{-1/2}  
           !n=grad phi/|grad phi|=(hx g,hy g,-g)
           !div n=hxx g + hx gx +hyy g + hy gy - gz
           !gx=(-1/2)arclen^{-3/2}arclenx
           curvHT_VOF=hxx*g+hx*gx+hyy*g+hy*gy
          else if (matstatus.eq.0) then
           curvHT_VOF=curvHT_LS
          else
           print *,"matstatus corrupt for XYZ curvature coeff"
           stop
          endif

         else
          print *,"dimension bust"
          stop
         endif

        else
         print *,"levelrz invalid"
         stop
        endif

        curvHT_choice=curvHT_VOF

       else if (vof_height_function.eq.0) then
        ! do nothing
       else
        print *,"vof_height_function invalid"
        stop
       endif
   
      else 
       print *,"overall_crossing_status invalid"
       stop
      endif

      if (n1d.eq.one) then ! im material on top
       curvHT_choice=-curvHT_choice
      else if (n1d.eq.-one) then ! im material on bottom
       ! do nothing
      else 
       print *,"n1d invalid"
       stop
      endif

      end subroutine analyze_heights
       


! find reconstruction from cell averaged data.
! point value lives at the center of the cell in 3D.
      subroutine minmod3D(data,xpos,slopes,intercept)
      IMPLICIT NONE

      real(amrex_real) intercept
      real(amrex_real) data(3,3,3)
      real(amrex_real) xpos(3,3,3,3)
      real(amrex_real) slopes(3)
      integer ii,jj,kk,dir
      real(amrex_real) slope_minus,slope_plus
      real(amrex_real) xm,xc,xp,dm,dc,dp

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

      real(amrex_real) intercept
      real(amrex_real) data(D_DECL(3,3,3))
      real(amrex_real) xpos(D_DECL(3,3,3),SDIM)
      real(amrex_real) slopes(SDIM)
      integer ii,jj,kk,dir
      real(amrex_real) slope_minus,slope_plus
      real(amrex_real) xm,xc,xp,dm,dc,dp

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

      integer, INTENT(in) :: sdim_parm
      integer :: dir
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: nhalf0
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: phi0
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real), INTENT(in) :: nn(sdim_parm)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,sdim_parm)
      real(amrex_real), INTENT(in) :: x(sdim_parm)
 
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
      real(amrex_real) function atan_verify(x)
      use probcommon_module
      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: x

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      atan_verify=atan(x)
      if ((atan_verify.ge.-half*MOF_PI-EPS2).and. &
          (atan_verify.le.half*MOF_PI+EPS2)) then
       !do nothing
      else
       print *,"atan out of range"
       stop
      endif
 
      return
      end function atan_verify

      subroutine put_angle_in_range(angle_init,angle_init_range,sdim)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: sdim
      real(amrex_real), INTENT(in) :: angle_init(sdim-1)
      real(amrex_real), INTENT(out) :: angle_init_range(sdim-1)
      integer :: dir

      do dir=1,sdim-1
       angle_init_range(dir)=angle_init(dir)
      enddo

      if (sdim.eq.2) then
       if ((angle_init(1).ge.-Pi-EPS2).and. &
           (angle_init(1).le.Pi+EPS2)) then
        ! do nothing
       else
        print *,"angle_init(1) out of range"
        stop
       endif
      else if (sdim.eq.3) then
       if ((angle_init(1).ge.-Pi-EPS2).and. &
           (angle_init(1).le.Pi+EPS2)) then
        ! do nothing
       else
        print *,"angle_init(1) out of range"
        stop
       endif
       if ((angle_init(sdim-1).ge.-Pi-EPS2).and. &
           (angle_init(sdim-1).le.Pi+EPS2)) then
        if (angle_init(sdim-1).ge.zero) then
         ! do nothing
        else if (angle_init(sdim-1).le.zero) then
         angle_init_range(sdim-1)=-angle_init_range(sdim-1)
         if (angle_init(1).le.zero) then
          angle_init_range(1)=angle_init(1)+Pi
         else if (angle_init(1).ge.zero) then
          angle_init_range(1)=angle_init(1)-Pi
         else
          print *,"angle_init(1) is NaN"
          stop
         endif
        else
         print *,"angle_init(sdim-1) is NaN"
         stop
        endif
       else
        print *,"angle_init(sdim-1) out of range"
        stop
       endif
      else
       print *,"sdim invalid"
       stop
      endif
                      
      end subroutine put_angle_in_range

! -pi < angle < pi 
! if (x,y) in first quadrant, then 0<theta<pi/2
! if (x,y) in fourth quadrant, then -pi/2<theta<0
! if (x,y) in second quadrant, then pi/2<theta<pi
! if (x,y) in third quadrant, then -pi<theta<-pi/2
      subroutine arctan2(y,x,angle)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: y,x
      real(amrex_real), INTENT(out) :: angle

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
      
      real(amrex_real), INTENT(in) :: x,y
      real(amrex_real), INTENT(out) :: z

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

      real(amrex_real), INTENT(out) :: vol
      real(amrex_real), INTENT(in) :: xlo1(SDIM),xhi1(SDIM)
      real(amrex_real), INTENT(in) :: xlo2(SDIM),xhi2(SDIM)
      real(amrex_real), INTENT(out) :: xloint(SDIM),xhiint(SDIM)
      real(amrex_real) vol1,vol2
      integer dir

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

      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(out) :: xT(SDIM)
      integer dir
     
      do dir=1,SDIM
       xT(dir)=x(dir)
      enddo
      if ((levelrz.eq.COORDSYS_CARTESIAN).or.(levelrz.eq.COORDSYS_RZ)) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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

      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real) x_ofs(SDIM)
      real(amrex_real), INTENT(in) :: ofs(SDIM)
      real(amrex_real), INTENT(out) :: xT(SDIM)
      integer dir
     
      do dir=1,SDIM
       x_ofs(dir)=x(dir)+ofs(dir)
       xT(dir)=x_ofs(dir)
      enddo
      if ((levelrz.eq.COORDSYS_CARTESIAN).or.(levelrz.eq.COORDSYS_RZ)) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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

      real(amrex_real),INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: bfact
      integer dir
      real(amrex_real), INTENT(out) :: dxmin
      real(amrex_real) delta
      real(amrex_real) delta_gauss
      real(amrex_real) RR

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
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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

      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: bfact
      integer dir
      real(amrex_real), INTENT(out) :: dxmax
      real(amrex_real) delta,RR,xL,xR

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
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        if (dir.eq.2) then
         RR=probhix
        endif
       else
        print *,"levelrz invalid"
        stop
       endif
       delta=delta*RR

       if (delta.gt.zero) then
        ! do nothing
       else
        print *,"delta invalid in get_dxmax"
        stop
       endif

       if ((delta.gt.dxmax).or.(dir.eq.1)) then
        dxmax=delta
       endif
      enddo ! dir=1..sdim

      end subroutine get_dxmax


      subroutine get_dxmaxLS(dx,bfact,dxmax)
      use probcommon_module
      use LegendreNodes

      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: bfact
      integer dir
      real(amrex_real), INTENT(out) :: dxmax
      real(amrex_real) delta,xL,xR

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

       if (delta.gt.zero) then
        ! do nothing
       else
        print *,"delta invalid get_dxmaxLS"
        stop
       endif

       if ((delta.gt.dxmax).or.(dir.eq.1)) then
        dxmax=delta
       endif
      enddo ! dir=1..sdim

      end subroutine get_dxmaxLS

      subroutine set_dimdec(DIMS(fabdim), &
                      fablo,fabhi,ngrow)
      IMPLICIT NONE

      integer, INTENT(out) :: DIMDEC(fabdim)
      integer, INTENT(in) :: fablo(SDIM)
      integer, INTENT(in) :: fabhi(SDIM)
      integer, INTENT(in) :: ngrow
      integer dir

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

      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: vel(SDIM)
      real(amrex_real), INTENT(out) :: velT(SDIM)
      integer dir
     
      do dir=1,SDIM
       velT(dir)=vel(dir)
      enddo
      if ((levelrz.eq.COORDSYS_CARTESIAN).or.(levelrz.eq.COORDSYS_RZ)) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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

      real(amrex_real) data(D_DECL(3,3,3))
      real(amrex_real) xpos(D_DECL(3,3,3),SDIM)
      real(amrex_real) x(SDIM)
      real(amrex_real) x0(SDIM)
      real(amrex_real) H
      real(amrex_real) intercept
      real(amrex_real) slopes(SDIM)
      integer dir,i1,j1,k1,k1lo,k1hi
      real(amrex_real) mindata,maxdata

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
      do k1=k1lo,k1hi
      do j1=1,3
      do i1=1,3
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

      real(amrex_real) AA(2,2)
      real(amrex_real) XX(2)
      real(amrex_real) BB(2)
      real(amrex_real) det
      integer matstatus

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
  
       ! 1<=gravity_dir<=sdim
      subroutine fort_derive_gravity_dir(gravity_vector_in,gravity_dir) &
      bind(c,name='fort_derive_gravity_dir')
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: gravity_vector_in(SDIM)
      integer, INTENT(out)  :: gravity_dir
      integer :: dir

      gravity_dir=1
      do dir=2,SDIM
       if (abs(gravity_vector_in(dir)).gt. &
           abs(gravity_vector_in(gravity_dir))) then
        gravity_dir=dir
       else if (abs(gravity_vector_in(dir)).le. &
                abs(gravity_vector_in(gravity_dir))) then
        !do nothing
       else
        print *,"gravity_vector_in invalid"
        stop
       endif
      enddo !dir==2..sdim

      return
      end subroutine fort_derive_gravity_dir

      subroutine fort_derive_gravity(gravity_vector_in,gravity) 
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: gravity_vector_in(SDIM)
      real(amrex_real), INTENT(out)  :: gravity
      integer :: gravity_dir

      call fort_derive_gravity_dir(gravity_vector_in,gravity_dir)
      gravity=gravity_vector_in(gravity_dir)

      return
      end subroutine fort_derive_gravity

      subroutine fort_check_operation_flag_MAC(operation_flag) &
      bind(c,name='fort_check_operation_flag_MAC')

      IMPLICIT NONE
      integer, INTENT(in) :: operation_flag

      if ((operation_flag.eq.OP_PRESGRAD_MAC).or. &
          (operation_flag.eq.OP_PRES_CELL_TO_MAC).or. &
          (operation_flag.eq.OP_POTGRAD_TO_MAC).or. &
          (operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. &
          (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC).or. &
          (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. &
          (operation_flag.eq.OP_UGRAD_MAC).or. &
          (operation_flag.eq.OP_ISCHEME_MAC).or. &
          (operation_flag.eq.OP_UGRAD_COUPLING_MAC).or. &
          (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then
       ! do nothing
      else
       print *,"operation_flag invalid"
       stop
      endif

      return
      end subroutine fort_check_operation_flag_MAC

      subroutine get_iten(im1,im2,iten)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: im1,im2
      integer, INTENT(out) :: iten
      integer im,im_opp
      integer im_iter
      integer previous_count

      if ((im1.lt.1).or.(im1.gt.num_materials).or. & 
          (im2.lt.1).or.(im2.gt.num_materials).or. &
          (im1.eq.im2)) then
       print *,"im1,im2 mismatch(A) im1,im2=",im1,im2
       print *,"num_materials=",num_materials
       stop
      endif

      if (im1.lt.im2) then
       im=im1
       im_opp=im2
      else if (im2.lt.im1) then
       im=im2
       im_opp=im1
      else
       print *,"im1,im2 mismatch(A2) im1,im2=",im1,im2
       print *,"num_materials=",num_materials
       stop
      endif
      if (im_opp.gt.im) then
       previous_count=0
       do im_iter=1,im-1
        previous_count=previous_count+num_materials-im_iter
       enddo
       iten=previous_count+im_opp-im
      else
       print *,"im or im_opp invalid: ",im,im_opp
       stop
      endif

      return  
      end subroutine get_iten    


      subroutine get_inverse_iten(im1,im2,iten)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(out) :: im1,im2
      integer, INTENT(in) :: iten
      integer im,im_opp,iten_test

      im1=0
      im2=0
      do im=1,num_materials
      do im_opp=im+1,num_materials
       call get_iten(im,im_opp,iten_test)
       if (iten.eq.iten_test) then
        im1=im
        im2=im_opp
       endif
      enddo
      enddo

      if ((im1.lt.1).or.(im1.gt.num_materials).or. & 
          (im2.lt.1).or.(im2.gt.num_materials).or. &
          (im1.eq.im2)) then
       print *,"im1,im2 mismatch(B) im1,im2=",im1,im2
       print *,"num_materials=",num_materials
       stop
      endif

      return  
      end subroutine get_inverse_iten    

        subroutine check_inbox(xx,xsten,nhalf,inboxflag)
        use probcommon_module

        IMPLICIT NONE

        integer nhalf
        real(amrex_real) dx
        real(amrex_real) xsten(-nhalf:nhalf,SDIM)
        real(amrex_real) xx(SDIM)
        integer inboxflag,dir

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

      real(amrex_real) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      real(amrex_real) perim,perim1,perim2

      perim1=((y2-y1)*(z3-z1)-(y3-y1)*(z2-z1))**2
      perim1=perim1+((x2-x1)*(z3-z1)-(z2-z1)*(x3-x1))**2
      perim1=perim1+((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))**2

      perim2=((y2-y4)*(z3-z4)-(y3-y4)*(z2-z4))**2
      perim2=perim2+((x2-x4)*(z3-z4)-(z2-z4)*(x3-x4))**2
      perim2=perim2+((x2-x4)*(y3-y4)-(x3-x4)*(y2-y4))**2

      perim=half*(sqrt(perim1)+sqrt(perim2))

      return
      end subroutine planearea

      real(amrex_real) function angle_err(a1,a2) 
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: a1,a2

      angle_err=min(abs(a1-a2),abs(a1-a2+two*Pi))
      angle_err=min(angle_err,abs(a1-a2-two*Pi))

      return
      end function angle_err

      real(amrex_real) function getcutvol(r,L)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: r,L
 
      getcutvol=Pi*r*r*(r+sqrt(r*r-L*L))- &
        (Pi/three)*(r**3+(r*r-L*L)**(3.0/2.0))

      return
      end function getcutvol

      subroutine getsphere(L,V,angle)
      IMPLICIT NONE

      real(amrex_real) L,V,angle
      real(amrex_real) r1,r2,v1,v2,rmid,vmid
      integer i

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

      real(amrex_real) function hsprime(phi,cutoff)
      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: phi,cutoff

      if ((phi.ge.cutoff).or.(phi.le.-cutoff)) then
       hsprime=0
      else
       hsprime=(1.0/(2.0*cutoff))*(1.0+cos(Pi*phi/cutoff))
      endif
      return
      end function hsprime

      integer function sign_funct(LS)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: LS

      if (LS.lt.zero) then
       sign_funct=-1
      else if (LS.ge.zero) then
       sign_funct=1
      else
       print *,"LS bust in sign_funct  LS=",LS
       stop
      endif

      end function sign_funct

      real(amrex_real) function hssign(phi,cutoff)
      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: phi,cutoff

      hssign=two*hs(phi,cutoff)-one
      return
      end function hssign

      real(amrex_real) function bellcurve(phi,cutoff)
      IMPLICIT NONE
      real(amrex_real) phi,cutoff

      if (cutoff.le.zero) then
       print *,"cutoff must be positive"
       stop
      endif
      bellcurve=exp(-(phi/cutoff)**2)

      return
      end function bellcurve

      real(amrex_real) function slopewt(phi,cutoff)
      IMPLICIT NONE
      real(amrex_real) phi,cutoff

      if ((phi.ge.cutoff).or.(phi.le.-cutoff)) then
       slopewt=zero
      else
       slopewt=one+cos(Pi*phi/cutoff)
      endif

      return
      end function slopewt

      real(amrex_real) function hs(phi,cutoff)
      use probcommon_module
      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: phi
      real(amrex_real), INTENT(in) :: cutoff

      if (phi.ge.cutoff) then
        hs=one
      else if (phi.le.-cutoff) then
        hs=zero
      else
        hs=phi/(two*cutoff)+sin(Pi*phi/cutoff)/(two*Pi)+half
        if ((hs.lt.-EPS_3_2).or.(hs.gt.one+EPS_3_2)) then
         print *,"hs sanity check failed: ",hs
         stop
        else if (hs.lt.zero) then
         hs=zero
        else if (hs.gt.one) then
         hs=one
        else if ((hs.ge.zero).and.(hs.le.one)) then
         ! do nothing
        else
         print *,"hs() is NaN: ",hs
         stop
        endif
      endif

      return
      end function hs

      real(amrex_real) function hs_smooth(phi,cutoff)
      use probcommon_module
      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: phi
      real(amrex_real) :: phi_scale
      real(amrex_real), INTENT(in) :: cutoff

      if (phi.ge.cutoff) then
        hs_smooth=one
      else if (phi.le.-cutoff) then
        hs_smooth=zero
      else if (abs(phi).le.cutoff) then
       phi_scale=half*(phi/cutoff+one)
       hs_smooth=(phi_scale**4)*(-20.0d0*(phi_scale**3)+ &
           70.0d0*(phi_scale**2)-84.0d0*phi_scale+35.0d0)
       if ((hs_smooth.lt.-EPS_3_2).or.(hs_smooth.gt.one+EPS_3_2)) then
        print *,"hs_smooth sanity check failed: ",hs_smooth
        stop
       else if (hs_smooth.lt.zero) then
        hs_smooth=zero
       else if (hs_smooth.gt.one) then
        hs_smooth=one
       else if ((hs_smooth.ge.zero).and.(hs_smooth.le.one)) then
        ! do nothing
       else
        print *,"hs_smooth() is NaN: ",hs_smooth
        stop
       endif
      else
       print *,"phi corrupt hs_smooth: ",phi
       stop
      endif

      return
      end function hs_smooth

      real(amrex_real) function hsprime_smooth(phi,cutoff)
      use probcommon_module
      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: phi
      real(amrex_real) :: phi_scale
      real(amrex_real), INTENT(in) :: cutoff

      if (phi.ge.cutoff) then
        hsprime_smooth=zero
      else if (phi.le.-cutoff) then
        hsprime_smooth=zero
      else if (abs(phi).le.cutoff) then
       phi_scale=half*(phi/cutoff+one)
       hsprime_smooth=(phi_scale**3)*(-140.0d0*(phi_scale**3)+ &
           420.0d0*(phi_scale**2)-420.0d0*phi_scale+140.0d0)
       if (hsprime_smooth.lt.-EPS_3_2) then
        print *,"hsprime_smooth sanity check failed: ",hsprime_smooth
        stop
       else if (hsprime_smooth.lt.zero) then
        hsprime_smooth=zero
       else if (hsprime_smooth.ge.zero) then
        ! do nothing
       else
        print *,"hsprime_smooth() is NaN: ",hsprime_smooth
        stop
       endif
      else
       print *,"phi corrupt hsprime_smooth: ",phi
       stop
      endif

      return
      end function hsprime_smooth


      real(amrex_real) function hs_scale(phi,cutoff)
      use probcommon_module
      IMPLICIT NONE
      real(amrex_real) phi,cutoff

      if (phi.ge.cutoff) then
        hs_scale=one
      else if (phi.le.-cutoff) then
        hs_scale=zero
      else
        hs_scale=half*(half+(phi/cutoff)+half*((phi/cutoff)**2)- &
         (cos(two*Pi*phi/cutoff)-one)/(four*Pi*Pi)+ &
         (phi+cutoff)*sin(Pi*phi/cutoff)/(cutoff*Pi))

        if ((hs_scale.lt.-EPS_3_2).or.(hs_scale.gt.one+EPS_3_2)) then
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

      real(amrex_real) function approximate(start_x,end_x,start_y,end_y,my_x)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) start_x,end_x,start_y,end_y,my_x

      if (abs(start_x-end_x).lt.EPS15) then
       approximate=(end_y+start_y)/2.0
      else
       approximate=start_y+(end_y-start_y)*(my_x-start_x)/ &
         (end_x-start_x)
      endif
      return
      end function approximate


      subroutine consistent_materials(vfrac,cen)
      use probcommon_module

      IMPLICIT NONE

      real(amrex_real), INTENT(inout) :: vfrac(num_materials)
      real(amrex_real), INTENT(inout) :: cen(SDIM,num_materials)

      integer imaterial,dir
      real(amrex_real) voftotal


      voftotal=zero

      do imaterial=1,num_materials
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
       if (is_rigid(imaterial).eq.0) then
        voftotal=voftotal+vfrac(imaterial)
       else if (is_rigid(imaterial).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid GLOBALUTIL.F90"
        stop
       endif
      enddo ! imaterial

      if (voftotal.gt.zero) then
       !do nothing
      else
       print *,"vacuum bust in consistent materials"
       print *,"voftotal= ",voftotal
       stop
      endif

      do imaterial=1,num_materials
       vfrac(imaterial)=vfrac(imaterial)/voftotal
      enddo 

      return
      end subroutine consistent_materials

      subroutine boxrefine(clo,chi,flo,fhi,dir)
      IMPLICIT NONE

      integer clo(SDIM),chi(SDIM)
      integer flo(SDIM),fhi(SDIM)
      integer boxdir,dir

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

      integer DIMDEC(data)
      integer datalo(SDIM),datahi(SDIM)
      integer dir

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

      integer DIMDEC(data)
      integer datalo(SDIM),datahi(SDIM)
      integer dir

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

      real(amrex_real), INTENT(out) :: vol
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: rzflag
      real(amrex_real) RCENTER

      if (nhalf.lt.1) then
       print *,"nhalf invalid gridvol"
       stop
      endif

      if (rzflag.eq.COORDSYS_CARTESIAN) then
       if (SDIM.eq.3) then
        vol=(xsten(1,1)-xsten(-1,1))*(xsten(1,2)-xsten(-1,2))* &
         (xsten(1,SDIM)-xsten(-1,SDIM))
       else if (SDIM.eq.2) then
        vol=(xsten(1,1)-xsten(-1,1))*(xsten(1,2)-xsten(-1,2))
       else
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       RCENTER=half*abs(xsten(-1,1)+xsten(1,1))
       vol=two*Pi*RCENTER* &
         (xsten(1,1)-xsten(-1,1))*(xsten(1,2)-xsten(-1,2))
      else if (rzflag.eq.COORDSYS_CYLINDRICAL) then
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



      subroutine gridarea(xsten,nhalf,rzflag,dir,side,area)
      IMPLICIT NONE

      real(amrex_real), INTENT(out) :: area
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: rzflag
      integer, INTENT(in) :: dir,side
      integer iside,dir2
      real(amrex_real) RR,R1,R2

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

      area=one
      do dir2=1,SDIM
       if (dir2.ne.dir+1) then
        area=area*(xsten(1,dir2)-xsten(-1,dir2))
        if (area.gt.zero) then
         ! do nothing
        else
         print *,"area bust"
         stop
        endif
       endif
      enddo  ! dir2

      if (rzflag.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (rzflag.eq.COORDSYS_RZ) then
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
       else
        print *,"dir invalid gridarea 2"
        stop
       endif
       area=area*two*Pi*RR

      else if (rzflag.eq.COORDSYS_CYLINDRICAL) then

       dir2=1
       R1=xsten(-1,1)
       R2=xsten(1,1)
       if (R2+R1.eq.zero) then
        print *,"R2+R1 invalid"
        stop
       endif

       if (dir.eq.0) then ! R face
        RR=abs(xsten(iside,1))
       else if (dir.eq.1) then ! Theta face
        RR=one
       else if ((dir.eq.2).and.(SDIM.eq.3)) then ! Z face
        RR=abs(xsten(0,1))
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

      integer bfact
      real(amrex_real) xnodes(bfact)
      real(amrex_real) xlo
      integer ielem
      integer lo_e,inodes
      real(amrex_real) dx
  
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
       print *,"bfact invalid20 (element_nodes1D): ",bfact
       stop
      endif

      end subroutine element_nodes1D


       ! returns the Gauss Lobatto points in element "ielem"
      subroutine element_GLnodes1D(xnodes,xlo,ielem,lo_e,dx,bfact)
      use LegendreNodes

      IMPLICIT NONE

      integer bfact
      real(amrex_real) xnodes(bfact+1)
      real(amrex_real) xlo
      integer ielem
      integer lo_e,inodes
      real(amrex_real) dx
  
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

      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: fablo(SDIM)
      integer, INTENT(in) :: i,j,k,bfact
      integer dir,icrit
      real(amrex_real), dimension(:), allocatable :: xsub

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
      subroutine gridstenMAC(x,xlo,i,j,k,fablo,bfact,dx,nhalf,grid_type)
      IMPLICIT NONE 

      integer, INTENT(in) :: nhalf
      integer, INTENT(in) :: grid_type
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: fablo(SDIM)
      integer, INTENT(in) :: i,j,k,bfact
      integer dir,icrit
      real(amrex_real), dimension(:), allocatable :: xsub
      integer :: box_type(SDIM)

      if ((grid_type.ge.-1).and.(grid_type.le.5)) then
       ! do nothing
      else
       print *,"put breakpoint here to see the caller"
       print *,"grid_type invalid gridstenMAC: grid_type ", &
          grid_type
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

      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: fablo(SDIM)
      integer, INTENT(in) :: i,j,k,bfact
      integer dir,icrit
      real(amrex_real), dimension(:), allocatable :: xsub

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

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: i

      integer, INTENT(out) :: isub,ielem

      integer :: ipos

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

      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: fablo(SDIM)
      real(amrex_real) xdomlo,xofs,xabs,xnode
      integer, INTENT(in) :: i,bfact,dir
      integer side
      integer icell,imac,isten,signflag,ishift,i_e,i_n

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
        x(side*2*icell)=xphys_of_xcomp(dir-1,xabs)

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
    
        x(side*(2*imac-1))=xphys_of_xcomp(dir-1,xabs)
       enddo ! imac
      enddo ! side

      return
      end subroutine gridsten1D

       ! x(0) is a gauss-Lobatto node
       ! x(-1) and x(1) are gauss nodes....
       ! 1<=dir<=sdim
      subroutine gridsten1DMAC(x,xlo,i,fablo,bfact,dx,dir,nhalf)
      use LegendreNodes

      IMPLICIT NONE

      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: fablo(SDIM)
      real(amrex_real) xdomlo,xofs,xabs,xnode
      integer, INTENT(in) :: i,bfact,dir
      integer side
      integer icell,imac,isten,signflag,ishift,i_e,i_n

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
        x(side*(2*icell-1))=xphys_of_xcomp(dir-1,xabs)

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
    
        x(side*2*imac)=xphys_of_xcomp(dir-1,xabs)
       enddo ! imac
      enddo ! side

      return
      end subroutine gridsten1DMAC



      subroutine get_istar(icomp,istar)
      IMPLICIT NONE

      integer icomp,icomplocal,dir
      integer istar(3)

      if (NCOMP_STENCIL/2.eq.27) then
       !do nothing
      else
       print *,"ncomp_stencil invalid get_istar"
       stop
      endif

      if ((icomp.lt.1).or.(icomp.gt.NCOMP_STENCIL/2)) then
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

      integer icomp,dir
      integer istar(3)

      if (NCOMP_STENCIL/2.eq.27) then
       !do nothing
      else
       print *,"ncomp_stencil invalid put_istar"
       stop
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
      if ((icomp.lt.1).or.(icomp.gt.NCOMP_STENCIL/2)) then
       print *,"icomp invalid"
       stop
      endif

      return
      end subroutine put_istar 

      subroutine growntilebox( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,ng)
      IMPLICIT NONE

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(out) :: growlo(3),growhi(3)
      integer, INTENT(in) :: ng
      integer dir2

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
      enddo ! dir2=1..sdim
      
      return
      end subroutine growntilebox


      subroutine growntilebox_TILE( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,ng)
      IMPLICIT NONE

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(out) :: growlo(3),growhi(3)
      integer, INTENT(in) :: ng
      integer dir2

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

      integer, intent(in) :: ic,i,bfact_c,bfact_f
      real(amrex_real), intent(out) :: wt
      integer, parameter :: nhalf=1
      integer dir_index
      integer fablo(SDIM)
      real(amrex_real) intlo,inthi
      real(amrex_real) dxc(SDIM)
      real(amrex_real) dxf(SDIM)
      real(amrex_real) xlo(SDIM)
      real(amrex_real) xc(-nhalf:nhalf)
      real(amrex_real) xf(-nhalf:nhalf)

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid: ",bfact_c
       stop
      endif

      dxc(1)=one
      dxf(1)=half

      xlo(1)=zero
      dir_index=1
      fablo(1)=0
      call gridsten1D(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
      call gridsten1D(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

      if (bfact_f.eq.1) then
       if (abs(xf(1)-xf(-1)-dxf(1)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xf invalid 1: ",abs(xf(1)-xf(-1)-dxf(1))
        stop
       endif
       if (abs(xf(1)+xf(-1)-two*xf(0)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xf invalid 2: ",abs(xf(1)+xf(-1)-two*xf(0))
        stop
       endif
      else if (bfact_f.gt.1) then
       if ((xf(1)-xf(0).gt.zero).and. &
           (xf(0)-xf(-1).gt.zero)) then
        !do nothing
       else
        print *,"xf invalid 3"
        stop
       endif
      else
       print *,"bfact_f invalid"
       stop
      endif

      if (bfact_c.eq.1) then
       if (abs(xc(1)-xc(-1)-dxc(1)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xc invalid1: ",abs(xc(1)-xc(-1)-dxc(1))
        stop
       endif
       if (abs(xc(1)+xc(-1)-two*xc(0)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xc invalid2: ",abs(xc(1)+xc(-1)-two*xc(0))
        stop
       endif
      else if (bfact_c.gt.1) then
       if ((xc(1)-xc(0).gt.zero).and. &
           (xc(0)-xc(-1).gt.zero)) then
        !do nothing
       else
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



      subroutine intersect_weight_avg_refine(ic,ic2,i,i2,bfact_c,bfact_f,wt)
      use probcommon_module

      integer, intent(in) :: ic,ic2,i,i2,bfact_c,bfact_f
      real(amrex_real), intent(out) :: wt
      integer, parameter :: nhalf=1
      integer dir_index
      integer fablo(SDIM)
      real(amrex_real) intlo,inthi
      real(amrex_real) dxc(SDIM)
      real(amrex_real) dxf(SDIM)
      real(amrex_real) xlo(SDIM)
      real(amrex_real) xc(-nhalf:nhalf)
      real(amrex_real) xf(-nhalf:nhalf)
      real(amrex_real) xclo,xchi,xflo,xfhi

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid: ",bfact_c
       stop
      endif

      dxc(1)=one
      dxf(1)=half

      xlo(1)=zero
      dir_index=1
      fablo(1)=0
      call gridsten1D(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
      call gridsten1D(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

      if (bfact_f.eq.1) then
       if (abs(xf(1)-xf(-1)-dxf(1)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xf invalid 1: ",abs(xf(1)-xf(-1)-dxf(1))
        stop
       endif
       if (abs(xf(1)+xf(-1)-two*xf(0)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xf invalid 2: ",abs(xf(1)+xf(-1)-two*xf(0))
        stop
       endif
      else if (bfact_f.gt.1) then
       if ((xf(1)-xf(0).gt.zero).and. &
           (xf(0)-xf(-1).gt.zero)) then
        !do nothing
       else
        print *,"xf invalid 3"
        stop
       endif
      else
       print *,"bfact_f invalid"
       stop
      endif

      if (bfact_c.eq.1) then
       if (abs(xc(1)-xc(-1)-dxc(1)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xc invalid1: ",abs(xc(1)-xc(-1)-dxc(1))
        stop
       endif
       if (abs(xc(1)+xc(-1)-two*xc(0)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xc invalid2: ",abs(xc(1)+xc(-1)-two*xc(0))
        stop
       endif
      else if (bfact_c.gt.1) then
       if ((xc(1)-xc(0).gt.zero).and. &
           (xc(0)-xc(-1).gt.zero)) then
        !do nothing
       else
        print *,"xc invalid 3"
        stop
       endif
      else
       print *,"bfact_c invalid"
       stop
      endif

      if (ic2.eq.0) then
       xclo=xc(-1)
       xchi=xc(0)
      else if (ic2.eq.1) then
       xclo=xc(0)
       xchi=xc(1)
      else
       print *,"ic2 invalid: ",ic2
       stop
      endif

      if (i2.eq.0) then
       xflo=xf(-1)
       xfhi=xf(0)
      else if (i2.eq.1) then
       xflo=xf(0)
       xfhi=xf(1)
      else
       print *,"i2 invalid: ",i2
       stop
      endif

      intlo=max(xflo,xclo)
      inthi=min(xfhi,xchi)
      if (intlo.ge.inthi) then
       wt=zero
      else
       wt=(inthi-intlo)/(xchi-xclo) ! average down
      endif
      if (abs(wt).le.VOFTOL) then
       wt=zero
      else if (abs(wt-one).le.VOFTOL) then
       wt=one
      else if ((wt.gt.zero).and.(wt.lt.one)) then
       !do nothing
      else
       print *,"wt invalid 2"
       stop
      endif

      return
      end subroutine intersect_weight_avg_refine


      subroutine intersect_weight_avg_COPY(ic,i,bfact_c,bfact_f,wt, &
         dir,dir_current)
      use probcommon_module

      integer, intent(in) :: ic,i,bfact_c,bfact_f,dir,dir_current
      real(amrex_real), intent(out) ::  wt

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

      integer, intent(in) :: ic,i,bfact_c,bfact_f
      real(amrex_real), intent(out) :: wt
      integer, parameter :: nhalf=1
      integer dir_index
      integer fablo(SDIM)
      real(amrex_real) intlo,inthi
      real(amrex_real) dxc(SDIM)
      real(amrex_real) dxf(SDIM)
      real(amrex_real) xlo(SDIM)
      real(amrex_real) xc(-nhalf:nhalf)
      real(amrex_real) xf(-nhalf:nhalf)

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid: ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      dxc(1)=one
      dxf(1)=half

      xlo(1)=zero
      dir_index=1
      fablo(1)=0
      call gridsten1D(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
      call gridsten1D(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

      if (bfact_f.eq.1) then
       if (abs(xf(1)-xf(-1)-dxf(1)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xf invalid 1: ",abs(xf(1)-xf(-1)-dxf(1))
        stop
       endif
       if (abs(xf(1)+xf(-1)-two*xf(0)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xf invalid 2: ",abs(xf(1)+xf(-1)-two*xf(0))
        stop
       endif
      else if (bfact_f.gt.1) then
       if ((xf(1)-xf(0).gt.zero).and. &
           (xf(0)-xf(-1).gt.zero)) then
        !do nothing
       else
        print *,"xf invalid 3"
        stop
       endif
      else
       print *,"bfact_f invalid"
       stop
      endif

      if (bfact_c.eq.1) then
       if (abs(xc(1)-xc(-1)-dxc(1)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xc invalid1: ",abs(xc(1)-xc(-1)-dxc(1))
        stop
       endif
       if (abs(xc(1)+xc(-1)-two*xc(0)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xc invalid2: ",abs(xc(1)+xc(-1)-two*xc(0))
        stop
       endif
      else if (bfact_c.gt.1) then
       if ((xc(1)-xc(0).gt.zero).and. &
           (xc(0)-xc(-1).gt.zero)) then
        !do nothing
       else
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
       if ((abs(intlo-xf(-1)).lt.EPS_8_4*dxf(1)).or. &
           (abs(inthi-xf(1)).lt.EPS_8_4*dxf(1))) then
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


      subroutine intersect_weight_interp_refine(ic,ic2,i,i2,bfact_c,bfact_f,wt)
      use probcommon_module

      integer, intent(in) :: ic,ic2,i,i2,bfact_c,bfact_f
      real(amrex_real), intent(out) :: wt
      integer, parameter :: nhalf=1
      integer dir_index
      integer fablo(SDIM)
      real(amrex_real) intlo,inthi
      real(amrex_real) dxc(SDIM)
      real(amrex_real) dxf(SDIM)
      real(amrex_real) xlo(SDIM)
      real(amrex_real) xc(-nhalf:nhalf)
      real(amrex_real) xf(-nhalf:nhalf)
      real(amrex_real) xclo,xchi,xflo,xfhi

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid: ",bfact_c
       stop
      endif

      dxc(1)=one
      dxf(1)=half

      xlo(1)=zero
      dir_index=1
      fablo(1)=0
      call gridsten1D(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
      call gridsten1D(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

      if (bfact_f.eq.1) then
       if (abs(xf(1)-xf(-1)-dxf(1)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xf invalid 1: ",abs(xf(1)-xf(-1)-dxf(1))
        stop
       endif
       if (abs(xf(1)+xf(-1)-two*xf(0)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xf invalid 2: ",abs(xf(1)+xf(-1)-two*xf(0))
        stop
       endif
      else if (bfact_f.gt.1) then
       if ((xf(1)-xf(0).gt.zero).and. &
           (xf(0)-xf(-1).gt.zero)) then
        !do nothing
       else
        print *,"xf invalid 3"
        stop
       endif
      else
       print *,"bfact_f invalid"
       stop
      endif

      if (bfact_c.eq.1) then
       if (abs(xc(1)-xc(-1)-dxc(1)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xc invalid1: ",abs(xc(1)-xc(-1)-dxc(1))
        stop
       endif
       if (abs(xc(1)+xc(-1)-two*xc(0)).lt.EPS_8_2) then
        !do nothing
       else
        print *,"xc invalid2: ",abs(xc(1)+xc(-1)-two*xc(0))
        stop
       endif
      else if (bfact_c.gt.1) then
       if ((xc(1)-xc(0).gt.zero).and. &
           (xc(0)-xc(-1).gt.zero)) then
        !do nothing
       else
        print *,"xc invalid 3"
        stop
       endif
      else
       print *,"bfact_c invalid"
       stop
      endif

      if (ic2.eq.0) then
       xclo=xc(-1)
       xchi=xc(0)
      else if (ic2.eq.1) then
       xclo=xc(0)
       xchi=xc(1)
      else
       print *,"ic2 invalid: ",ic2
       stop
      endif

      if (i2.eq.0) then
       xflo=xf(-1)
       xfhi=xf(0)
      else if (i2.eq.1) then
       xflo=xf(0)
       xfhi=xf(1)
      else
       print *,"i2 invalid: ",i2
       stop
      endif

      intlo=max(xflo,xclo)
      inthi=min(xfhi,xchi)
      if (intlo.ge.inthi) then
       wt=zero
      else
       wt=(inthi-intlo)/(xfhi-xflo) ! interpolate
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
       if ((abs(intlo-xflo).lt.EPS_8_4*dxf(1)).or. &
           (abs(inthi-xfhi).lt.EPS_8_4*dxf(1))) then
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
      end subroutine intersect_weight_interp_refine



      subroutine intersect_weight_interp_COPY(ic,i,bfact_c,bfact_f,wt, &
         dir,dir_current)
      use probcommon_module

      integer, intent(in) :: ic,i,bfact_c,bfact_f,dir,dir_current
      real(amrex_real), intent(out) :: wt

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

      integer, intent(in) :: ic,i,bfact_c,bfact_f
      real(amrex_real), intent(out) :: wt
      integer, parameter :: nhalf=1
      integer dir_index
      integer fablo(SDIM)
      real(amrex_real) intlo,inthi
      real(amrex_real) dxc(SDIM)
      real(amrex_real) dxf(SDIM)
      real(amrex_real) xlo(SDIM)
      real(amrex_real) xc(-nhalf:nhalf)
      real(amrex_real) xf(-nhalf:nhalf)
      integer denom

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


        dxc(1)=one
        dxf(1)=half

        xlo(1)=zero
        dir_index=1
        fablo(1)=0
        call gridsten1DMAC(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
        call gridsten1DMAC(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

        if (bfact_f.eq.1) then
         if (abs(xf(1)-xf(-1)-dxf(1)).ge.EPS_8_2) then
          print *,"xf invalid 1: ",abs(xf(1)-xf(-1)-dxf(1))
          stop
         endif
         if (abs(xf(1)+xf(-1)-two*xf(0)).ge.EPS_8_2) then
          print *,"xf invalid 2: ",abs(xf(1)+xf(-1)-two*xf(0))
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
         if (abs(xc(1)-xc(-1)-dxc(1)).ge.EPS_8_2) then
          print *,"xc invalid1: ",abs(xc(1)-xc(-1)-dxc(1))
          stop
         endif
         if (abs(xc(1)+xc(-1)-two*xc(0)).ge.EPS_8_2) then
          print *,"xc invalid2: ",abs(xc(1)+xc(-1)-two*xc(0))
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
         if ((abs(intlo-xf(-1)).lt.EPS_8_4*dxf(1)).or. &
             (abs(inthi-xf(1)).lt.EPS_8_4*dxf(1))) then
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

      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: i,j,k,level
      integer isten,dir
      integer dummy_input
 
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

      integer, INTENT(in) :: nhalf,dir
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf)
      integer, INTENT(in) :: i,level
      integer isten
 
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
      subroutine gridstenMAC_level(x,i,j,k,level,nhalf,grid_type)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: nhalf
      integer, INTENT(in) :: grid_type
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: i,j,k,level
      integer isten,dir,ii,jj,kk
      integer box_type(SDIM)
      integer dummy_input

      if ((grid_type.ge.-1).and.(grid_type.le.5)) then
       ! do nothing
      else
       print *,"put breakpoint here to see the caller"
       print *,"grid_type invalid gridstenMAC_level: grid_type ", &
          grid_type
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
       print *,"put breakpoint here to see the caller"
       print *,"level invalid gridstenMAC_level"
       print *,"level,cache_max_level ",level,cache_max_level
       print *,"nhalf,grid_type ",nhalf,grid_type
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

      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: i,j,k,level
      integer isten,dir

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

      integer, INTENT(in) :: nhalf,normdir
      real(amrex_real), INTENT(out) :: x(-nhalf:nhalf)
      integer, INTENT(in) :: i,level
      integer isten,ii

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

      integer, intent(in) :: ic,i,bfact_c,bfact_f
      real(amrex_real), intent(out) :: wt
      integer, parameter :: nhalf=1
      integer dir_index
      integer fablo(SDIM)
      real(amrex_real) intlo,inthi
      real(amrex_real) dxc(SDIM)
      real(amrex_real) dxf(SDIM)
      real(amrex_real) xlo(SDIM)
      real(amrex_real) xc(-nhalf:nhalf)
      real(amrex_real) xf(-nhalf:nhalf)

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

        dxc(1)=one
        dxf(1)=half

        xlo(1)=zero
        dir_index=1
        fablo(1)=0
        call gridsten1DMAC(xc,xlo,ic,fablo,bfact_c,dxc,dir_index,nhalf)
        call gridsten1DMAC(xf,xlo,i,fablo,bfact_f,dxf,dir_index,nhalf)

        if (bfact_f.eq.1) then
         if (abs(xf(1)-xf(-1)-dxf(1)).ge.EPS_8_2) then
          print *,"xf invalid 1: ",abs(xf(1)-xf(-1)-dxf(1))
          stop
         endif
         if (abs(xf(1)+xf(-1)-two*xf(0)).ge.EPS_8_2) then
          print *,"xf invalid 2: ",abs(xf(1)+xf(-1)-two*xf(0))
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
         if (abs(xc(1)-xc(-1)-dxc(1)).ge.EPS_8_2) then
          print *,"xc invalid1: ",abs(xc(1)-xc(-1)-dxc(1))
          stop
         endif
         if (abs(xc(1)+xc(-1)-two*xc(0)).ge.EPS_8_2) then
          print *,"xc invalid2: ",abs(xc(1)+xc(-1)-two*xc(0))
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
      integer function DIV_FLOOR(a,b)
      IMPLICIT NONE

      integer a,b,apos,div1

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

      integer ic,jc,kc
      integer coarse_index(3)
      integer stenlo(3),stenhi(3)
      integer bfact_c,bfact_f
      integer dir2

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

      integer ifine,jfine,kfine
      integer fine_index(3)
      integer stenlo(3),stenhi(3)
      integer bfact_c,bfact_f
      integer dir2

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

      integer, INTENT(in) :: grid_type
      integer, INTENT(out) :: box_type(SDIM)
      integer dir

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

      integer, INTENT(out) :: grid_type
      integer, INTENT(in) :: box_type(SDIM)

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

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: ifine,jfine,kfine
      integer fine_index(3)
      integer, INTENT(out) :: stenlo(3),stenhi(3)
      integer, INTENT(in) :: bfact_c,bfact_f
      integer dir2,denom
      integer :: box_type(SDIM)

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

      integer, INTENT(in) :: ic,jc,kc
      integer, INTENT(in) :: grid_type
      integer coarse_index(3)
      integer, INTENT(out) :: stenlo(3),stenhi(3)
      integer, INTENT(in) :: bfact_c,bfact_f
      integer dir2
      integer :: box_type(SDIM)

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
       tilelo,tilehi,fablo,fabhi,growlo,growhi,ng,grid_type)
      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(out) :: growlo(3),growhi(3)
      integer, INTENT(in) :: ng
      integer dir2
      integer :: box_type(SDIM)

      growlo(3)=0
      growhi(3)=0

      if ((grid_type.ge.-1).and.(grid_type.le.5)) then
       ! do nothing
      else
       print *,"put breakpoint here to see the caller"
       print *,"grid_type invalid growntilebox mac"
       print *,"ng=",ng
       print *,"grid_type=",grid_type
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

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(out) :: growlo(3),growhi(3)
      integer, INTENT(in) :: ng
      integer dir2

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

      subroutine growntileboxTENSOR( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,dir)
      IMPLICIT NONE

      integer dir
      integer tilelo(SDIM),tilehi(SDIM)
      integer fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer dir2

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

      integer dir
      integer tilelo(SDIM),tilehi(SDIM)
      integer fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer dir2

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

      integer, INTENT(in) :: i,j,k,bfact
      integer, INTENT(out) :: stripstat

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

      integer, INTENT(in) :: i,j,k,bfact,dir
      integer dir2
      integer, INTENT(out) :: elemlo(3)
      integer, INTENT(out) :: elemhi(3)

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

      integer, INTENT(in)  :: i,j,k,bfact,dir
      integer, INTENT(out) :: stripstat

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
      integer qbfact,tbfact
      integer i1,j1,k
      real(amrex_real) GQws(0:tbfact,0:qbfact-1,1:tbfact) 
      real(amrex_real) y(0:qbfact-1)
      real(amrex_real) yGL(0:tbfact)
      real(amrex_real), dimension(:),allocatable :: bwGL,l
      real(amrex_real) xpt

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
       stenhi,dx,xfine,fcoarse,fxfine)
      use probcommon_module
      use LagrangeInterpolationPolynomial
      use LegendreNodes

      IMPLICIT NONE
      integer, INTENT(in) :: nvar
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: grid_type ! -1:ggg; 0:lgg; 1:glg; 2:ggl
      integer, INTENT(in) :: stenhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xfine(SDIM)
      real(amrex_real), INTENT(in) ::  &
        fcoarse(D_DECL(0:stenhi(1),0:stenhi(2),0:stenhi(3)),nvar)
      real(amrex_real), INTENT(out) :: fxfine(nvar)
 
      integer i1
      integer i,j,k
      integer ii,jj,kk
      integer khi
      integer dir
      integer n
      real(amrex_real) y(0:bfact-1)
      real(amrex_real) yGL(0:bfact)
      real(amrex_real) wt
      real(amrex_real), dimension(:),allocatable :: temp,tempGL
      real(amrex_real), dimension(:),allocatable :: ypoints,bwG
      real(amrex_real), dimension(:),allocatable :: ypointsGL,bwGL
      real(amrex_real), dimension(:,:),allocatable :: lg
      real(amrex_real) wtsum
      real(amrex_real) local_data
      integer :: box_type(SDIM)

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
          print *,"put breakpoint here to see the caller"
          print *,"abs(local_data) overflow SEM_INTERP_ELEMENT"
          print *,"n,nvar,local_data ",n,nvar,local_data
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
      if (abs(wtsum-one).le.EPS_10_4) then
       !do nothing
      else
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

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: levelrz_in
      integer, INTENT(in) :: dir,nc
      real(amrex_real), INTENT(in) :: RRface(0:bfact)
    
      integer, INTENT(in) :: operation_flag
      real(amrex_real), INTENT(in) :: dx
      integer, INTENT(in) :: bctype(2)
      integer local_bctype(2)
      real(amrex_real), INTENT(in) :: x_sep(2)
      real(amrex_real), INTENT(in) :: bcvalue(2)
      real(amrex_real), INTENT(in) :: vel(0:bfact)
      real(amrex_real), INTENT(in) :: source_side(2)
      real(amrex_real) local_source_side(2)
      real(amrex_real), INTENT(in) :: source(0:bfact-1) 
      real(amrex_real), INTENT(out) :: dest_grad(0:bfact)
      real(amrex_real), INTENT(out) ::  dest_interp(0:bfact)

      integer i1,j1,isten
      real(amrex_real) y(0:bfact-1)
      real(amrex_real) y_extend(0:bfact+1)
      real(amrex_real) yGL(0:bfact)  ! Gauss-Lobatto Legendre points
      real(amrex_real) yGL_extend(0:bfact+1)
      real(amrex_real) yLT(0:bfact+1)
      real(amrex_real) yRT(0:bfact+1)
      real(amrex_real) yLRT(0:bfact+1)
      real(amrex_real) yLRTextrap(0:bfact-1)
      real(amrex_real) yLTextrap(0:bfact)
      real(amrex_real) yRTextrap(0:bfact)

      real(amrex_real) yLTextrapEXT(0:bfact)
      real(amrex_real) yRTextrapEXT(0:bfact)

      real(amrex_real) wMAT(0:bfact-1,0:bfact-1)
      real(amrex_real) wMATGL(0:bfact,0:bfact)
      real(amrex_real) wMAT_extend(0:bfact+1,0:bfact+1)
      real(amrex_real) wLT(0:bfact+1,0:bfact+1)
      real(amrex_real) wRT(0:bfact+1,0:bfact+1)

      real(amrex_real) wLT_EXT(0:bfact,0:bfact)
      real(amrex_real) wRT_EXT(0:bfact,0:bfact)

      real(amrex_real) wLRT(0:bfact+1,0:bfact+1)
      real(amrex_real) dx_element
      real(amrex_real) dx_ends
      real(amrex_real) sum1,sum2,A,B,C,D,det
      real(amrex_real) PLINE(0:bfact+1)
      real(amrex_real) PLINE2(0:bfact+1)
      integer AMR_boundary_flag

      call fort_check_operation_flag_MAC(operation_flag)

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

      if (levelrz_in.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz_in.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz_in.eq.COORDSYS_CYLINDRICAL) then

       if (operation_flag.eq.OP_UGRAD_MAC) then ! tensor derivatives

        if ((nc.lt.1).or. &
            (nc.gt.SDIM)) then
         print *,"nc invalid"
         stop
        endif

       else if (operation_flag.eq.OP_PRESGRAD_MAC) then ! grad p_MAC
        ! do nothing
       else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 
        ! do nothing
       else if (operation_flag.eq.OP_UNEW_CELL_TO_MAC) then ! u^{Cell->Mac}
        ! do nothing
       else if (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC) then 
        ! do nothing
       else if (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC) then
        ! do nothing
       else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection
        ! do nothing
       else if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then ! coupling
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

       ! dx_element=dx * bfact
       ! since bfact>=2 => -1<y(0)<-1+1/bfact => 
       ! 0<y(0)+1<1/bfact => 0<(y(0)+1)/2<1/(2 bfact) =>
       ! 0<dx_ends<dx * bfact /(2 bfact) =>
       ! 0<dx_ends<dx/2
      dx_ends=dx_element*(y(0)+one)/two
      if (dx_ends.lt.half*dx) then
       ! do nothing
      else
       print *,"(dx_ends.ge.half*dx)"
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
       ! exterior BC left, extend right
      yLT(0)=-one
      yLT(bfact+1)=y_extend(bfact+1)
       ! exterior BC right, extend left
      yRT(0)=y_extend(0)
      yRT(bfact+1)=one
       ! exterior BC left and right
      yLRT(0)=-one
      yLRT(bfact+1)=one
      do i1=1,bfact
       yLT(i1)=y(i1-1)
       yRT(i1)=y(i1-1)
       yLRT(i1)=y(i1-1)
      enddo ! i1

      yLTextrap(bfact)=y_extend(bfact+1) ! extrap left, extend right
      yLTextrapEXT(bfact)=one ! extrap left, exterior BC right
      do i1=0,bfact-1
       yLTextrap(i1)=y(i1)
       yLTextrapEXT(i1)=y(i1)
      enddo

      yRTextrap(0)=y_extend(0) ! extrap right, extend left
      yRTextrapEXT(0)=-one ! extrap right, exterior BC left
      do i1=0,bfact-1
       yRTextrap(i1+1)=y(i1)
       yRTextrapEXT(i1+1)=y(i1)
      enddo

       ! extrap left and right
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
         ! extrap left, exterior BC right
        wRT_EXT(i1,j1)=cache_wRT_EXT(bfact,i1,j1,SPTYPE) 
         ! extrap right, exterior BC left
        wLT_EXT(i1,j1)=cache_wLT_EXT(bfact,i1,j1,SPTYPE) 
       enddo
      enddo


       ! at this point, PLINE(0) and PLINE(bfact+1) need to be filled.
      do i1=0,bfact-1
       PLINE(i1+1)=source(i1)
      enddo

      local_source_side(1)=source_side(1)
      local_source_side(2)=source_side(2)
      local_bctype(1)=bctype(1)
      local_bctype(2)=bctype(2)
      if (bctype(1).eq.SEM_REFLECT_EVEN) then ! reflect even
       local_bctype(1)=SEM_INTERIOR
       local_source_side(1)=PLINE(1)
      endif
      if (bctype(2).eq.SEM_REFLECT_EVEN) then ! reflect even
       local_bctype(2)=SEM_INTERIOR
       local_source_side(2)=PLINE(bfact)
      endif
      if (bctype(1).eq.SEM_REFLECT_ODD) then ! reflect odd
       local_bctype(1)=SEM_INTERIOR
       local_source_side(1)=-PLINE(1)
      endif
      if (bctype(2).eq.SEM_REFLECT_ODD) then ! reflect odd
       local_bctype(2)=SEM_INTERIOR
       local_source_side(2)=-PLINE(bfact)
      endif

       ! SEM_IMAGE_BC_ALG is defined in PROBCOMMON.F90
      if (SEM_IMAGE_BC_ALG.eq.1) then

        ! this approximation is exact for linear functions:
        ! y=ax+b
        ! PLINE(1)=a(xlo+h/2)+b
        ! local_source_side(1)=2(a(xlo)+b)-[a(xlo+h/2)+b]=
        ! a(xlo)+b-a(h/2)=a(xlo-h/2)+b
       if (bctype(1).eq.SEM_DIRICHLET) then ! dirichlet -> reflect odd
        local_bctype(1)=SEM_INTERIOR
        local_source_side(1)=two*bcvalue(1)-PLINE(1)
       endif
       if (bctype(2).eq.SEM_DIRICHLET) then ! dirichlet -> reflect odd
        local_bctype(2)=SEM_INTERIOR
        local_source_side(2)=two*bcvalue(2)-PLINE(bfact)
       endif
       if (bctype(1).eq.SEM_NEUMANN) then ! neumann -> reflect even
        local_bctype(1)=SEM_INTERIOR
        local_source_side(1)=PLINE(1)-two*dx_ends*bcvalue(1)
       endif
       if (bctype(2).eq.SEM_NEUMANN) then ! neumann -> reflect even
        local_bctype(2)=SEM_INTERIOR
        local_source_side(2)=PLINE(bfact)+two*dx_ends*bcvalue(2)
       endif

      else if (SEM_IMAGE_BC_ALG.eq.0) then
       ! do nothing
       ! i.e. do not replace SEM_NEUMANN and SEM_DIRICHLET with
       ! SEM_INTERIOR with reflection ghost value.
      else
       print *,"SEM_IMAGE_BC_ALG invalid: ",SEM_IMAGE_BC_ALG
       stop
      endif

      AMR_boundary_flag=0

      if ((local_bctype(1).eq.SEM_COARSE_NEXT_TO_FINE).or. &
          (local_bctype(1).eq.SEM_FINE_NEXT_TO_COARSE)) then
       AMR_boundary_flag=1
       local_bctype(1)=SEM_INTERIOR
       if ((x_sep(1).gt.zero).and. &
           (x_sep(1).lt.two)) then
        y_extend(0)=-one-x_sep(1)
        yRT(0)=y_extend(0)
        yRTextrap(0)=y_extend(0)
       else
        print *,"x_sep invalid"
        stop
       endif
        ! cannot use the cache; cache is for AMR_boundary_flag=0
       call polyinterp_Dmatrix(bfact+1,yRT,wRT)
      endif

      if ((local_bctype(2).eq.SEM_COARSE_NEXT_TO_FINE).or. &
          (local_bctype(2).eq.SEM_FINE_NEXT_TO_COARSE)) then
       AMR_boundary_flag=1
       local_bctype(2)=SEM_INTERIOR
       if ((x_sep(2).gt.zero).and. &
           (x_sep(2).lt.two)) then
        y_extend(bfact+1)=one+x_sep(2)
        yLT(bfact+1)=y_extend(bfact+1)
        yLTextrap(bfact)=y_extend(bfact+1)
       else
        print *,"x_sep invalid"
        stop
       endif
        ! cannot use the cache; cache is for AMR_boundary_flag=0
       call polyinterp_Dmatrix(bfact+1,yLT,wLT)
      endif

      if ((local_bctype(1).eq.SEM_INTERIOR).and. &
          (local_bctype(2).eq.SEM_INTERIOR)) then 
       PLINE(0)=local_source_side(1)
       PLINE(bfact+1)=local_source_side(2)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        y_extend,yGL_extend)
      else if ((local_bctype(1).eq.SEM_INTERIOR).and. &
               (local_bctype(2).eq.SEM_EXTRAP)) then 
       PLINE(0)=local_source_side(1)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yRTextrap,yGL_extend)
      else if ((local_bctype(1).eq.SEM_DIRICHLET).and. &
               (local_bctype(2).eq.SEM_EXTRAP)) then 
       PLINE(0)=bcvalue(1)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yRTextrapEXT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_NEUMANN).and. &
               (local_bctype(2).eq.SEM_EXTRAP)) then 
        ! P'(xj)=sum wij * P(xi) * 2/dx_element
        ! bcvalue(1)=sum wi0 * P(i) * 2/dx_element
        ! P(0)=(bcvalue(1)-sum_{i\ne 0} wi0 * P(i) * 2/dx_element)*
        !      (dx_element/2) * (1/w00)
       sum1=bcvalue(1)*dx_element/two
       do i1=1,bfact
        sum1=sum1-PLINE(i1)*wLT_EXT(i1,0) 
       enddo
       if (wLT_EXT(0,0).ne.zero) then
        ! do nothing
       else
        print *,"wLT_EXT bust"
        stop
       endif
       PLINE(0)=sum1/wLT_EXT(0,0)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yRTextrapEXT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_EXTRAP).and. &
               (local_bctype(2).eq.SEM_INTERIOR)) then 
       do i1=0,bfact-1
        PLINE(i1)=source(i1)
       enddo
       PLINE(bfact)=local_source_side(2)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yLTextrap,yGL_extend)
      else if ((local_bctype(1).eq.SEM_EXTRAP).and. &
               (local_bctype(2).eq.SEM_DIRICHLET)) then 
       do i1=0,bfact-1
        PLINE(i1)=source(i1)
       enddo
       PLINE(bfact)=bcvalue(2)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yLTextrapEXT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_EXTRAP).and. &
               (local_bctype(2).eq.SEM_NEUMANN)) then 
        ! P'(xj)=sum wij * P(xi) * 2/dx_element
        ! bcvalue(2)=sum wir * P(i) * 2/dx_element
        ! P(r)=(bcvalue(2)-sum_{i\ne r} wir * P(i) * 2/dx_element)*
        !      (dx_element/2) * (1/wrr)
       do i1=0,bfact-1
        PLINE(i1)=source(i1)
       enddo
       sum2=bcvalue(2)*dx_element/two
       do i1=0,bfact-1
        sum2=sum2-PLINE(i1)*wRT_EXT(i1,bfact)
       enddo
       if (wRT_EXT(bfact,bfact).ne.zero) then
        ! do nothing
       else
        print *,"wRT_EXT bust"
        stop
       endif
       PLINE(bfact)=sum2/wRT_EXT(bfact,bfact)
       call poly_change_basis(bfact,bfact+1,PLINE,PLINE2, &
        yLTextrapEXT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_EXTRAP).and. &
               (local_bctype(2).eq.SEM_EXTRAP)) then 
       do i1=0,bfact-1
        PLINE(i1)=source(i1)
       enddo
       call poly_change_basis(bfact-1,bfact+1,PLINE,PLINE2, &
        yLRTextrap,yGL_extend)
      else if ((local_bctype(1).eq.SEM_NEUMANN).and. &
               (local_bctype(2).eq.SEM_NEUMANN)) then 
        ! P'(xj)=sum wij * P(xi) * 2/dx_element
        ! bcvalue(1)=sum wi0 * P(i) * 2/dx_element
        ! bcvalue(2)=sum wir * P(i) * 2/dx_element
        ! bcvalue(1)*dx_element/2=sum wi0 * P(i) 
        ! bcvalue(2)*dx_element/2=sum wir * P(i) 
        ! sum1=w00 * P(0) + w_{r+1,0}P(r+1)
        ! sum2=w0,r+1 * P(r+1) + w_{r+1,r+1}P(r+1)
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
       if (det.ne.zero) then
        ! do nothing
       else
        print *,"determinent bust"
        stop
       endif
        ! Cramer's rule: Ax=B
        ! x_i = det(A_i)/det(A)
        ! A_i is the matrix formed by replacing the ith column of A
        ! by the column vector B.
       PLINE(0)=(D*sum1-B*sum2)/det
       PLINE(bfact+1)=(A*sum2-C*sum1)/det 
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLRT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_NEUMANN).and. &
               (local_bctype(2).eq.SEM_INTERIOR)) then 
       PLINE(bfact+1)=local_source_side(2)
       sum1=bcvalue(1)*dx_element/two
       do i1=1,bfact+1
        sum1=sum1-PLINE(i1)*wLT(i1,0) 
       enddo
       if (wLT(0,0).ne.zero) then
        ! do nothing
       else
        print *,"wLT bust"
        stop
       endif
       PLINE(0)=sum1/wLT(0,0)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_NEUMANN).and. &
               (local_bctype(2).eq.SEM_DIRICHLET)) then 
       PLINE(bfact+1)=bcvalue(2)
       sum1=bcvalue(1)*dx_element/two
       do i1=1,bfact+1
        sum1=sum1-PLINE(i1)*wLRT(i1,0)
       enddo
       if (wLRT(0,0).ne.zero) then
        ! do nothing
       else
        print *,"wLRT bust"
        stop
       endif
       PLINE(0)=sum1/wLRT(0,0)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLRT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_DIRICHLET).and. &
               (local_bctype(2).eq.SEM_NEUMANN)) then 
       PLINE(0)=bcvalue(1)
       sum2=bcvalue(2)*dx_element/two
       do i1=0,bfact
        sum2=sum2-PLINE(i1)*wLRT(i1,bfact+1)
       enddo
       if (wLRT(bfact+1,bfact+1).ne.zero) then
        ! do nothing
       else
        print *,"wLRT bust"
        stop
       endif
       PLINE(bfact+1)=sum2/wLRT(bfact+1,bfact+1)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLRT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_INTERIOR).and. &
               (local_bctype(2).eq.SEM_NEUMANN)) then
       PLINE(0)=local_source_side(1)
       sum2=bcvalue(2)*dx_element/two
       do i1=0,bfact
        sum2=sum2-PLINE(i1)*wRT(i1,bfact+1)
       enddo
       if (wRT(bfact+1,bfact+1).ne.zero) then
        ! do nothing
       else
        print *,"wRT bust"
        stop
       endif
       PLINE(bfact+1)=sum2/wRT(bfact+1,bfact+1)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yRT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_DIRICHLET).and. &
               (local_bctype(2).eq.SEM_DIRICHLET)) then 
       PLINE(0)=bcvalue(1)
       PLINE(bfact+1)=bcvalue(2)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yLRT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_INTERIOR).and. &
               (local_bctype(2).eq.SEM_DIRICHLET)) then 
       PLINE(0)=local_source_side(1)
       PLINE(bfact+1)=bcvalue(2)
       call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
        yRT,yGL_extend)
      else if ((local_bctype(1).eq.SEM_DIRICHLET).and. &
               (local_bctype(2).eq.SEM_INTERIOR)) then 
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

      if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

       do i1=0,bfact
        if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then
         dest_interp(i1)=dest_interp(i1)*vel(i1)
        else if (nc.eq.SEM_T+1) then ! temperature
         dest_interp(i1)=dest_interp(i1)*vel(i1)
        else
         print *,"nc invalid in lineGRAD"
         stop
        endif
        dest_grad(i1)=zero
       enddo ! i1=0..bfact

      else if (operation_flag.eq.OP_PRESGRAD_MAC) then ! grad p_MAC
        ! do nothing
      else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 
        ! do nothing
      else if (operation_flag.eq.OP_UNEW_CELL_TO_MAC) then ! u^{Cell->Mac}
        ! do nothing
      else if (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC) then 
        ! do nothing
      else if (operation_flag.eq.OP_UGRAD_MAC) then ! rate of strain
       ! do nothing
      else if (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC) then
       ! do nothing
      else if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then ! coupling terms
       ! do nothing
      else
       print *,"operation_flag invalid3"
       stop
      endif

      do isten=0,bfact

       if (levelrz_in.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz_in.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
       else if (levelrz_in.eq.COORDSYS_CYLINDRICAL) then

        if (operation_flag.eq.OP_UGRAD_MAC) then

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

        else if (operation_flag.eq.OP_PRESGRAD_MAC) then ! grad p_MAC
         ! do nothing
        else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 
         ! do nothing
        else if (operation_flag.eq.OP_UNEW_CELL_TO_MAC) then ! u^{Cell->Mac}
         ! do nothing
        else if (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC) then 
         ! do nothing
        else if (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC) then
         ! do nothing
        else if (operation_flag.eq.OP_ISCHEME_MAC) then
         ! do nothing (advection)
        else if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then
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

      real(amrex_real), INTENT(in) :: dx
      integer, INTENT(in) :: maskSEM
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: source(0:bfact) 
      real(amrex_real), INTENT(out) :: dest(0:bfact-1)
      real(amrex_real), INTENT(out) :: destdiv(0:bfact-1)

      integer i1,j1
      real(amrex_real) y(0:bfact-1)
      real(amrex_real) yGL(0:bfact)
      real(amrex_real) wMATGL(0:bfact,0:bfact)
      real(amrex_real) wtsum,wt
      real(amrex_real) intlo,inthi
      real(amrex_real) dx_element

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

      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(out) :: stenlo(3),stenhi(3)
      integer, INTENT(in) :: nsten
      integer dir2
      integer ii(3)

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

      integer dir
      integer DIMDEC(data)
      integer DIMDEC(dataf)
      integer datalo(SDIM),datahi(SDIM)
      integer datalof(SDIM),datahif(SDIM)

      call dim_to_box(DIMS(data),datalo,datahi)
      call boxrefine(datalo,datahi,datalof,datahif,dir)
      call box_to_dim(DIMS(dataf),datalof,datahif)

      return
      end subroutine dimrefine

      subroutine get_longdir(lo,hi,longdir,longlo,longhi)
      IMPLICIT NONE

      integer lo(SDIM),hi(SDIM)
      integer longdir
      integer longlo,longhi
      integer maxlen,tempmax,dir

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

      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten_center(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: stencil_offset(SDIM)
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      real(amrex_real), INTENT(out) :: WT
      integer :: WT_ok
      integer :: dir
      real(amrex_real) :: denom
      real(amrex_real) theta(SDIM)

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
      subroutine deriv_from_grid_util(data_in,disp_dataptr,data_out)
      use probcommon_module
      IMPLICIT NONE

      type(deriv_from_grid_parm_type), INTENT(in) :: data_in 
      real(amrex_real), INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: disp_dataptr
      type(interp_from_grid_out_parm_type), INTENT(out) :: data_out

      integer dir_local
      integer ilo(SDIM)
      integer ihi(SDIM)
      integer istep(SDIM)
      integer klosten,khisten,kstep
      integer nc
      integer data_comp
      integer isten,jsten,ksten
      integer ii(3)
      real(amrex_real) SGN_FACT
      real(amrex_real) wt_top,wt_bot
      real(amrex_real) xflux_sten(-3:3,SDIM)
      real(amrex_real) xhi_sten(-3:3,SDIM)
      real(amrex_real) xlo_sten(-3:3,SDIM)
      real(amrex_real) dx_sten(SDIM)
      real(amrex_real) dx_top
      real(amrex_real), target :: xtarget(SDIM)
      integer nhalf
      real(amrex_real) :: local_data_out
      real(amrex_real), allocatable, dimension(:) :: local_data_max
      real(amrex_real), allocatable, dimension(:) :: local_data_min
      real(amrex_real) :: scaling
      integer :: dummy_input

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

      nhalf=3 
      call gridstenMAC_level(xflux_sten,ilocal(1),ilocal(2),ilocal(SDIM), &
        data_in%level,nhalf,data_in%grid_type_flux)

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
        data_in%level,nhalf,data_in%grid_type_data)
      call gridstenMAC_level(xlo_sten,ilo(1),ilo(2),ilo(SDIM), &
        data_in%level,nhalf,data_in%grid_type_data)

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
       print *,"ncomp invalid in deriv_from_grid_util"
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

       local_data_max(nc)=-1.0D+20
       local_data_min(nc)=1.0D+20

       data_out%data_interp(nc)=zero

       do ksten=klosten,khisten,kstep
       do jsten=ilo(2),ihi(2),istep(2)
       do isten=ilo(1),ihi(1),istep(1)
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
           disp_dataptr,local_data_out)

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
        data_in2%xtarget=xtarget
        if ((data_in%dx(1).gt.zero).and. &
            (data_in%dx(2).gt.zero).and. &
            (data_in%dx(SDIM).gt.zero)) then
         ! do nothing
        else
         print *,"data_in%dx bust"
         stop
        endif
        data_in2%dx=data_in%dx
        data_in2%xlo=data_in%xlo
        data_in2%fablo=data_in%fablo
        data_in2%fabhi=data_in%fabhi
        allocate(data_out2%data_interp(data_in2%ncomp))
        call interp_from_grid_util(data_in2,disp_dataptr,data_out2)
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
                 data_out2%data_interp(nc)).le.EPS_10_4*scaling) then
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
      subroutine single_deriv_from_grid_util(data_in,disp_dataptr,data_out)
      use probcommon_module
      IMPLICIT NONE

      type(single_deriv_from_grid_parm_type), INTENT(in) :: data_in 
      real(amrex_real), INTENT(in), pointer, dimension(D_DECL(:,:,:)) :: &
        disp_dataptr
      type(interp_from_grid_out_parm_type), INTENT(out) :: data_out

      integer dir_local
      integer ilo(SDIM)
      integer ihi(SDIM)
      integer istep(SDIM)
      integer klosten,khisten,kstep
      integer isten,jsten,ksten
      integer ii(3)
      real(amrex_real) SGN_FACT
      real(amrex_real) wt_top,wt_bot
      real(amrex_real) xflux_sten(-3:3,SDIM)
      real(amrex_real) xhi_sten(-3:3,SDIM)
      real(amrex_real) xlo_sten(-3:3,SDIM)
      real(amrex_real) dx_sten(SDIM)
      real(amrex_real) dx_top
      real(amrex_real), target :: xtarget(SDIM)
      integer nhalf
      real(amrex_real) :: local_data_out
      real(amrex_real) :: local_data_max
      real(amrex_real) :: local_data_min
      real(amrex_real) :: scaling
      integer :: dummy_input

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

      nhalf=3 
      call gridstenMAC_level(xflux_sten,ilocal(1),ilocal(2),ilocal(SDIM), &
        data_in%level,nhalf,data_in%grid_type_flux)

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
        data_in%level,nhalf,data_in%grid_type_data)
      call gridstenMAC_level(xlo_sten,ilo(1),ilo(2),ilo(SDIM), &
        data_in%level,nhalf,data_in%grid_type_data)

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
      
      local_data_max=-1.0D+20
      local_data_min=1.0D+20

      data_out%data_interp(1)=zero

      do ksten=klosten,khisten,kstep
      do jsten=ilo(2),ihi(2),istep(2)
      do isten=ilo(1),ihi(1),istep(1)
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
         call safe_data_single(isten,jsten,ksten,disp_dataptr,local_data_out)

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
        data_in2%xtarget=xtarget
        data_in2%dx=data_in%dx
        data_in2%xlo=data_in%xlo
        data_in2%fablo=data_in%fablo
        data_in2%fabhi=data_in%fabhi
        allocate(data_out2%data_interp(1))
        call single_interp_from_grid_util(data_in2,disp_dataptr,data_out2)

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
                data_out2%data_interp(1)).le.EPS_10_4*scaling) then
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
        ! do nothing
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

      subroutine interp_from_grid_util(data_in,stateptr,data_out)
      use probcommon_module
      IMPLICIT NONE
 
      type(interp_from_grid_parm_type), INTENT(in) :: data_in 
      real(amrex_real), INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: &
          stateptr 
      type(interp_from_grid_out_parm_type), INTENT(out) :: data_out
      real(amrex_real) :: xsten(-3:3,SDIM)
      real(amrex_real) :: xsten_center(-3:3,SDIM)
      integer nhalf
      integer dir
      integer cell_index(SDIM)
      integer stencil_offset(SDIM)
      integer istenlo(3)
      integer istenhi(3)
      real(amrex_real) WT,total_WT
      integer isten,jsten,ksten
      integer im
      real(amrex_real), allocatable, dimension(:) :: local_data
     
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
       print *,"ncomp invalid in interp_from_grid_util"
       stop
      endif
      if (data_in%scomp.ge.1) then
       ! do nothing
      else
       print *,"scomp invalid"
       stop
      endif

      allocate(local_data(data_in%ncomp))

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

      do ksten=istenlo(3),istenhi(3)
      do jsten=istenlo(2),istenhi(2)
      do isten=istenlo(1),istenhi(1)

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
                stateptr,local_data(im))
        if ((local_data(im).ge.-1.0D+30).and. &
            (local_data(im).le.1.0D+30)) then

         data_out%data_interp(im)=data_out%data_interp(im)+WT*local_data(im)

        else
         print *,"local_data(im) overflow"
         print *,"isten,jsten,ksten ",isten,jsten,ksten
         print *,"stateptr ",stateptr
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


      subroutine single_interp_from_grid_util(data_in,stateptr,data_out)
      use probcommon_module
      IMPLICIT NONE
 
      type(single_interp_from_grid_parm_type), INTENT(in) :: data_in 
      real(amrex_real), INTENT(in), pointer, dimension(D_DECL(:,:,:)) :: stateptr 
      type(interp_from_grid_out_parm_type), INTENT(out) :: data_out
      real(amrex_real) :: xsten(-3:3,SDIM)
      real(amrex_real) :: xsten_center(-3:3,SDIM)
      integer nhalf
      integer dir
      integer cell_index(SDIM)
      integer stencil_offset(SDIM)
      integer istenlo(3)
      integer istenhi(3)
      real(amrex_real) WT,total_WT
      integer isten,jsten,ksten
      real(amrex_real) :: local_data
     
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

      do ksten=istenlo(3),istenhi(3)
      do jsten=istenlo(2),istenhi(2)
      do isten=istenlo(1),istenhi(1)

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

       call safe_data_single(isten,jsten,ksten,stateptr,local_data)
       if ((local_data.ge.-1.0D+30).and. &
           (local_data.le.1.0D+30)) then

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


       ! finds the cell that contains "x" 
      subroutine containing_cell_aux( &
         auxcomp,x,cell_index)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: auxcomp
      real(amrex_real), INTENT(in) :: x(3)
      integer, INTENT(out) :: cell_index(3)

      integer dir

      if ((auxcomp.ge.1).and.(auxcomp.le.fort_num_local_aux_grids)) then
       ! NINT=nearest int
       do dir=1,3
        ! x=(i-lo+1/2)dx+xlo  i=(x-xlo)/dx+lo-1/2
        cell_index(dir)=NINT( (x(dir)-contain_aux(auxcomp)%xlo3D(dir))/ &
          contain_aux(auxcomp)%dx3D-half )+contain_aux(auxcomp)%lo3D(dir)
       enddo ! dir=1..sdim
      else
       print *,"auxcomp invalid in interp_from_aux_grid"
       stop
      endif

      return
      end subroutine containing_cell_aux


      subroutine interp_from_aux_grid(auxcomp,x,LS)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: auxcomp
      real(amrex_real), INTENT(in) :: x(3)
      real(amrex_real), INTENT(out) :: LS
      integer, dimension(3) :: cell_index
      integer, dimension(3) :: cell_lo
      real(amrex_real) :: local_dx
      real(amrex_real) :: xgrid
      real(amrex_real), dimension(2,2,2,1) :: data_stencil
      real(amrex_real) :: wt_dist(3)
      integer :: nc
      real(amrex_real) :: LS_local(1)
      integer :: dir
      integer :: i,j,k
      integer :: i1,j1,k1

      nc=1

      if ((auxcomp.ge.1).and.(auxcomp.le.fort_num_local_aux_grids)) then
       call containing_cell_aux(auxcomp,x,cell_index)
       local_dx=contain_aux(auxcomp)%dx3D
       if (local_dx.gt.zero) then
        do dir=1,3
         xgrid=contain_aux(auxcomp)%xlo3D(dir)+ &
          (cell_index(dir)-contain_aux(auxcomp)%lo3D(dir)+half)*local_dx
         if (x(dir).le.xgrid) then
          cell_lo(dir)=cell_index(dir)-1
          wt_dist(dir)=(x(dir)-(xgrid-local_dx))/local_dx
         else if (x(dir).ge.xgrid) then
          cell_lo(dir)=cell_index(dir)
          wt_dist(dir)=(x(dir)-xgrid)/local_dx
         else
          print *,"auxcomp=",auxcomp
          print *,"fort_num_local_aux_grids=", &
            fort_num_local_aux_grids
          print *,"x or xgrid is NaN"
          print *,"dir=",dir
          print *,"local_dx=",local_dx
          print *,"x(1),x(2),x(3) ",x(1),x(2),x(3)
          print *,"xgrid = ",xgrid
          stop
         endif
        enddo ! dir=1..3
        do k1=0,1
        do j1=0,1
        do i1=0,1
         dir=1
         i=cell_lo(dir)+i1
         dir=2
         j=cell_lo(dir)+j1
         dir=3
         k=cell_lo(dir)+k1
         call safe_data3D(i,j,k,nc, &
           contain_aux(auxcomp)%LS3D, &
           data_stencil(i1+1,j1+1,k1+1,nc))
        enddo  ! i1
        enddo  ! j1
        enddo  ! k1
        call bilinear_interp_stencil3D(data_stencil,wt_dist,nc, &
          LS_local(1))
        LS=LS_local(1)
       else
        print *,"local_dx invalid"
        stop
       endif

      else
       print *,"auxcomp invalid in interp_from_aux_grid"
       stop
      endif

      end subroutine interp_from_aux_grid

      subroutine bilinear_interp_stencil3D(data_stencil,wt_dist_in, &
                      ncomp,data_interp)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: ncomp
      real(amrex_real), dimension(2,2,2,ncomp), INTENT(in) :: data_stencil
      real(amrex_real), INTENT(in) :: wt_dist_in(3)
      real(amrex_real) :: wt_dist(3)
      real(amrex_real), INTENT(out) :: data_interp(ncomp)
      integer :: dir
      real(amrex_real) :: c00,c01,c10,c11,c0,c1

      do dir=1,3
       wt_dist(dir)=wt_dist_in(dir)
       if ((wt_dist(dir).ge.zero).and. &
           (wt_dist(dir).le.one)) then
        !do nothing
       else if ((wt_dist(dir).ge.-EPS2).and. &
                (wt_dist(dir).le.zero)) then
        wt_dist(dir)=zero
       else if ((wt_dist(dir).le.one+EPS2).and. &
                (wt_dist(dir).ge.one)) then
        wt_dist(dir)=one
       else
        print *,"put breakpoint here to see the caller"
        print *,"wt_dist out of range bilinear_interp_stencil3D"
        print *,"dir,wt_dist ",dir,wt_dist(dir)
        print *,"ncomp= ",ncomp
        stop
       endif
      enddo ! dir=1..3

      do dir=1,ncomp
       c00 = data_stencil(1,1,1,dir)*(one-wt_dist(1)) + &
             data_stencil(2,1,1,dir)*wt_dist(1)
       c01 = data_stencil(1,1,2,dir)*(one-wt_dist(1)) + &
             data_stencil(2,1,2,dir)*wt_dist(1)
       c10 = data_stencil(1,2,1,dir)*(one-wt_dist(1)) + &
             data_stencil(2,2,1,dir)*wt_dist(1)
       c11 = data_stencil(1,2,2,dir)*(one-wt_dist(1)) + &
             data_stencil(2,2,2,dir)*wt_dist(1)

       c0 = c00*(one-wt_dist(2))+c10*wt_dist(2)
       c1 = c01*(one-wt_dist(2))+c11*wt_dist(2)
    
       data_interp(dir) = c0*(one-wt_dist(3))+c1*wt_dist(3)
      enddo ! dir=1..ncomp

      end subroutine bilinear_interp_stencil3D

      subroutine bilinear_interp_stencil(data_stencil,wt_dist_in, &
                      ncomp,data_interp)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: ncomp
      real(amrex_real), dimension(D_DECL(2,2,2),ncomp), INTENT(in) :: &
              data_stencil
      real(amrex_real), INTENT(in) :: wt_dist_in(SDIM)
      real(amrex_real) :: wt_dist(SDIM)
      real(amrex_real), INTENT(out) :: data_interp(ncomp)
      integer :: dir
      real(amrex_real) :: c00,c01,c10,c11,c0,c1

      do dir=1,SDIM

       wt_dist(dir)=wt_dist_in(dir)
       if ((wt_dist(dir).ge.zero).and. &
           (wt_dist(dir).le.one)) then
        !do nothing
       else if ((wt_dist(dir).ge.-EPS2).and. &
                (wt_dist(dir).le.zero)) then
        wt_dist(dir)=zero
       else if ((wt_dist(dir).le.one+EPS2).and. &
                (wt_dist(dir).ge.one)) then
        wt_dist(dir)=one
       else
        print *,"put breakpoint here to see the caller"
        print *,"wt_dist out of range (bilinear_interp_stencil)"
        print *,"dir,wt_dist ",dir,wt_dist(dir)
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

      integer nhalf
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) LS(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) curv
      integer inode,jnode,knode
      integer ii,jj,kk
      integer dir
      integer knlo,knhi
      real(amrex_real) mag,delx,rp,rm,rc
      real(amrex_real) nrm(SDIM)
      integer interface_found
      real(amrex_real) nn(D_DECL(0:1,0:1,0:1),SDIM)

      if (nhalf.lt.2) then
       print *,"nhalf invalid"
       stop
      endif 
      if (SDIM.eq.3) then
       if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
           (levelrz.eq.COORDSYS_CYLINDRICAL)) then
        ! do nothing
       else
        print *,"levelrz invalid"
        stop
       endif
       knlo=0
       knhi=1
      else if (SDIM.eq.2) then
       if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
           (levelrz.eq.COORDSYS_RZ).or. &
           (levelrz.eq.COORDSYS_CYLINDRICAL)) then
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
       do knode=knlo,knhi
       do jnode=0,1
       do inode=0,1
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

        if (levelrz.eq.COORDSYS_CARTESIAN) then
         rc=one
        else if (levelrz.eq.COORDSYS_RZ) then
         rc=one
        else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
         rc=xsten(0,1)
         if (rc.gt.zero) then
          !do nothing
         else
          print *,"rc invalid: ",rc
          stop
         endif 
        else
         print *,"levelrz invalid"
         stop
        endif 
        call prepare_normal(nrm,rc,mag,SDIM)
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
        if (delx.gt.zero) then
         !do nothing
        else
         print *,"delx invalid"
         stop
        endif
        if (dir.eq.1) then
         if (levelrz.eq.COORDSYS_CARTESIAN) then
          rp=one
          rm=one
          rc=one
         else if ((levelrz.eq.COORDSYS_RZ).or. &
                  (levelrz.eq.COORDSYS_CYLINDRICAL)) then
          rp=half*(xsten(2,1)+xsten(0,1))
          rm=half*(xsten(-2,1)+xsten(0,1))
          rc=half*(rp+rm)
          if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           if (rc.gt.zero) then
            !do nothing
           else
            print *,"rc invalid: ",rc
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
         if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
             (levelrz.eq.COORDSYS_RZ)) then
          rp=one
          rm=one
          rc=one
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
          rp=half*(xsten(2,1)+xsten(0,1))
          rm=half*(xsten(-2,1)+xsten(0,1))
          rc=half*(rp+rm)
          if (rc.gt.zero) then
           !do nothing
          else
           print *,"rc invalid: ",rc
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
       fort_caller_string, &
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
       verbose, &
       force_check, &
       gridno,ngrid,level,finest_level, &
       mf, &
       critical_cutoff_low, & ! e.g. -1.0D+30
       critical_cutoff_high)  ! e.g. 1.0D+30

      IMPLICIT NONE
  
      CHARACTER(:), ALLOCATABLE, intent(in) :: fort_caller_string

      integer, INTENT(in) :: datatype
      real(amrex_real), INTENT(in) :: warning_cutoff
      real(amrex_real), INTENT(in) :: critical_cutoff_low
      real(amrex_real), INTENT(in) :: critical_cutoff_high
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: growlo(SDIM),growhi(SDIM)
      integer growlotest(3),growhitest(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: scomp,ncomp,ndefined
      integer, INTENT(in) :: ngrow
      integer, INTENT(in) :: dir
      integer, INTENT(in) :: verbose
      integer, INTENT(in) :: force_check
      integer, INTENT(in) :: gridno,ngrid,level,finest_level
      real(amrex_real), INTENT(in), pointer :: mf(D_DECL(:,:,:),:)
      integer i,j,k,dir2
      integer box_type(SDIM)
      integer n
      integer n_singlelayer,n_interior,n_side,n_corner
      real(amrex_real) sum_interior,sum_side,sum_corner,sum_singlelayer
      real(amrex_real) max_interior,max_side,max_corner,max_singlelayer
      real(amrex_real) val
      integer noutside,noutside_single


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
         ngrow,dir)
       else if (datatype.eq.2) then
        call checkbound_array(fablo,fabhi, &
         mf, &
         1,-1)
       else
        print *,"datatype invalid"
        stop
       endif

       if (verbose.eq.2) then
        print *,"fort_caller_string=",fort_caller_string
        print *,"using gdb, put a break statement here to see the caller"
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
          growlotest,growhitest,ngrow,dir)
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

        do k=growlotest(3),growhitest(3)
        do j=growlotest(2),growhitest(2)
        do i=growlotest(1),growhitest(1)
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
          print *,"fort_caller_string=",fort_caller_string
          print *,"val.ge.critical_cutoff_high ",val,critical_cutoff_high
          print *,"val overflow val,dir,i,j,k,n,scomp ", &
           val,dir,i,j,k,n,scomp
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
          print *,"fort_caller_string=",fort_caller_string
          print *,"val.le.critical_cutoff_low ",val,critical_cutoff_low
          print *,"val out of bounds val,dir,i,j,k,n,scomp ", &
           val,dir,i,j,k,n,scomp
          print *,"bfact,level,finest_level ",bfact,level,finest_level
          do dir2=1,SDIM
           print *,"dir2,fablo,fabhi ",dir2,fablo(dir2),fabhi(dir2)
           print *,"dir2,tilelo,tilehi ",dir2,tilelo(dir2),tilehi(dir2)
          enddo
          stop
         else
          print *,"fort_caller_string=",fort_caller_string
          print *,"val undefined val,dir,i,j,k,n,scomp ", &
           val,dir,i,j,k,n,scomp
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
         print *,"fort_caller_string=",fort_caller_string
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

      subroutine abs_array_index3(i,j,k,Ni,Nj,Nk,abs_index)
      IMPLICIT NONE

      integer, INTENT(in) :: i,j,k,Ni,Nj,Nk
      integer, INTENT(out) :: abs_index

      if ((i.lt.1).or.(i.gt.Ni).or. &
          (j.lt.1).or.(j.gt.Nj).or. &
          (k.lt.1).or.(k.gt.Nk)) then
       print *,"index bust abs_array_index3"
       stop
      endif
      abs_index=(i-1)*Nj*Nk+(j-1)*Nk+k

      return
      end subroutine abs_array_index3

       ! finds the cell that contains "x" 
      subroutine containing_cell( &
       bfact,dx,xlo,lo,x,cell_index)
      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: lo(SDIM)
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      integer, INTENT(out) :: cell_index(SDIM)

      integer lo_e,i1,e_index,dir
      real(amrex_real) xnodes(bfact+1)
      real(amrex_real) xmap(SDIM)

       ! NINT=nearest int
      do dir=1,SDIM

       xmap(dir)=xcomp_of_xphys(dir-1,x(dir))

       if (bfact.eq.1) then  ! evenly spaced points
        ! x=(i-lo+1/2)dx+xlo  i=(x-xlo)/dx+lo-1/2
        cell_index(dir)=NINT( (xmap(dir)-xlo(dir))/dx(dir)-half )+lo(dir)
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
        e_index=NINT( (xmap(dir)-xlo(dir))/(bfact*dx(dir))-half )+lo_e
         ! returns the Gauss Lobatto points in element e_index
        call element_GLnodes1D(xnodes,xlo(dir),e_index,lo_e, &
         dx(dir),bfact)
        if (xmap(dir).le.xnodes(1)) then
         cell_index(dir)=e_index*bfact
        else if (xmap(dir).ge.xnodes(bfact+1)) then 
         cell_index(dir)=e_index*bfact+bfact-1
        else
         do i1=1,bfact
          if ((xmap(dir).ge.xnodes(i1)).and. &
              (xmap(dir).le.xnodes(i1+1))) then
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

      integer, INTENT(in) :: dir_mac
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: lo(SDIM)
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      integer, INTENT(out) :: mac_cell_index(SDIM)

      integer lo_e,e_index
      integer i1
      integer i1crit
      integer dir_local
      real(amrex_real) xnodes(bfact+1)
      real(amrex_real) xmap(SDIM)

      if ((dir_mac.ge.0).and.(dir_mac.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir_mac invalid"
       stop
      endif

       ! NINT=nearest int
      do dir_local=1,SDIM

       xmap(dir_local)=xcomp_of_xphys(dir_local-1,x(dir_local))

       if (bfact.eq.1) then  ! evenly spaced points
        ! dir_local!=dir_mac+1: x=(i-lo+1/2)dx+xlo  i=(x-xlo)/dx+lo-1/2
        ! dir_local==dir_mac+1: x=(i-lo)dx+xlo  i=(x-xlo)/dx+lo
        if (dir_local.ne.dir_mac+1) then
         mac_cell_index(dir_local)= &
          NINT( (xmap(dir_local)-xlo(dir_local))/dx(dir_local)-half )+ &
          lo(dir_local)
        else if (dir_local.eq.dir_mac+1) then
         mac_cell_index(dir_local)= &
          NINT( (xmap(dir_local)-xlo(dir_local))/dx(dir_local) )+ &
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
          NINT( (xmap(dir_local)-xlo(dir_local))/ &
          (bfact*dx(dir_local))-half )+lo_e
         ! returns the Gauss Lobatto points in element e_index
        call element_GLnodes1D(xnodes,xlo(dir_local),e_index,lo_e, &
         dx(dir_local),bfact)
       
        if (dir_local.ne.dir_mac+1) then
         if (xmap(dir_local).le.xnodes(1)) then
          mac_cell_index(dir_local)=e_index*bfact
         else if (xmap(dir_local).ge.xnodes(bfact+1)) then 
          mac_cell_index(dir_local)=e_index*bfact+bfact-1
         else
          do i1=1,bfact
           if ((xmap(dir_local).ge.xnodes(i1)).and. &
               (xmap(dir_local).le.xnodes(i1+1))) then
            mac_cell_index(dir_local)=e_index*bfact+i1-1
           endif
          enddo
         endif
        else if (dir_local.eq.dir_mac+1) then

         i1crit=1
         do i1=2,bfact+1
          if (xnodes(i1).gt.xnodes(i1crit)) then
           if (abs(xmap(dir_local)-xnodes(i1)).le. &
               abs(xmap(dir_local)-xnodes(i1crit))) then
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

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: lo(SDIM)
      real(amrex_real), INTENT(in) ::  x(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      integer, INTENT(out) :: node_index(SDIM)
      integer lo_e
      integer i1
      integer i1crit
      integer e_index,dir
      real(amrex_real) xnodes(bfact+1)
      real(amrex_real) xmap(SDIM)

       ! NINT=nearest int
      do dir=1,SDIM

       xmap(dir)=xcomp_of_xphys(dir-1,x(dir))

       if (bfact.eq.1) then ! evenly spaced points
        ! x=(i-lo)dx+xlo  i=(x-xlo)/dx+lo
        node_index(dir)=NINT( (xmap(dir)-xlo(dir))/dx(dir) )+lo(dir)
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
        e_index=NINT( (xmap(dir)-xlo(dir))/(bfact*dx(dir))-half )+lo_e
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
          if (abs(xmap(dir)-xnodes(i1)).le. &
              abs(xmap(dir)-xnodes(i1crit))) then
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

      function fort_CTML_FSI_mat_base(FSI_flag_local,im) &
      bind(c,name='fort_CTML_FSI_mat_base')
      use probcommon_module

      IMPLICIT NONE

      integer fort_CTML_FSI_mat_base
      integer, INTENT(in) :: FSI_flag_local
      integer, INTENT(in) :: im ! 1<=im<=num_materials

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16 fort_CTML_FSI_mat_base im:",im
       stop
      endif
      fort_CTML_FSI_mat_base=0

      if (FSI_flag_local.eq.FSI_SHOELE_CTML) then
#ifdef MVAHABFSI
       fort_CTML_FSI_mat_base=1
#else
       print *,"CTML(F): define MEHDI_VAHAB_FSI in GNUmakefile"
       stop
#endif
      else if ((FSI_flag_local.eq.FSI_FLUID).or. &
               (FSI_flag_local.eq.FSI_FLUID_NODES_INIT).or. &
               (FSI_flag_local.eq.FSI_PRESCRIBED_PROBF90).or. &
               (FSI_flag_local.eq.FSI_PRESCRIBED_NODES).or. &
               (FSI_flag_local.eq.FSI_ICE_PROBF90).or. &
               (FSI_flag_local.eq.FSI_ICE_STATIC).or. &
               (FSI_flag_local.eq.FSI_ICE_NODES_INIT).or. &
               (FSI_flag_local.eq.FSI_EULERIAN_ELASTIC).or. &
               (FSI_flag_local.eq.FSI_ICE_EULERIAN_ELASTIC).or. &
               (FSI_flag_local.eq.FSI_RIGID_NOTPRESCRIBED)) then
       fort_CTML_FSI_mat_base=0
      else
       print *,"FSI_flag_local invalid in fort_CTML_FSI_mat_base"
       print *,"im,FSI_flag_local ",im,FSI_flag_local
       stop
      endif

      return
      end function fort_CTML_FSI_mat_base

      function CTML_FSI_flagF()
      use probcommon_module

      IMPLICIT NONE

      integer CTML_FSI_flagF
      integer im

      CTML_FSI_flagF=0
      do im=1,num_materials
       if (fort_CTML_FSI_mat_base(FSI_flag(im),im).eq.1) then
#ifdef MVAHABFSI
        CTML_FSI_flagF=1
#else
        print *,"CTML(F): define MEHDI_VAHAB_FSI in GNUmakefile"
        stop
#endif
       else if (fort_CTML_FSI_mat_base(FSI_flag(im),im).eq.0) then
        ! do nothing
       else
        print *,"FSI_flag invalid in CTML_FSI_flagF"
        stop
       endif
      enddo ! im=1..num_materials

      return
      end function CTML_FSI_flagF


      function CTML_FSI_mat(im)
      use probcommon_module

      IMPLICIT NONE

      integer CTML_FSI_mat
      integer, INTENT(in) :: im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16 CTML_FSI_mat"
       stop
      endif
      CTML_FSI_mat=fort_CTML_FSI_mat_base(FSI_flag(im),im)

      return
      end function CTML_FSI_mat

      function fort_FSI_flag_valid_base(FSI_flag_local,im) &
      bind(c,name='fort_FSI_flag_valid_base')
      use probcommon_module

      IMPLICIT NONE

      integer fort_FSI_flag_valid_base
      integer, INTENT(in) :: FSI_flag_local
      integer, INTENT(in) :: im ! 1<=im<=num_materials

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16 fort_FSI_flag_valid_base: ",im
       stop
      endif

      fort_FSI_flag_valid_base=0

      if ((FSI_flag_local.eq.FSI_FLUID).or. &
          (FSI_flag_local.eq.FSI_FLUID_NODES_INIT).or. &
          (FSI_flag_local.eq.FSI_PRESCRIBED_PROBF90).or. &
          (FSI_flag_local.eq.FSI_PRESCRIBED_NODES).or. &
          (FSI_flag_local.eq.FSI_ICE_PROBF90).or. &
          (FSI_flag_local.eq.FSI_ICE_STATIC).or. &
          (FSI_flag_local.eq.FSI_ICE_NODES_INIT).or. &
          (FSI_flag_local.eq.FSI_RIGID_NOTPRESCRIBED).or. &
          (FSI_flag_local.eq.FSI_EULERIAN_ELASTIC).or. &
          (FSI_flag_local.eq.FSI_ICE_EULERIAN_ELASTIC).or. &
          (FSI_flag_local.eq.FSI_SHOELE_CTML)) then
       fort_FSI_flag_valid_base=1
      else
       print *,"FSI_flag_local invalid in fort_FSI_flag_valid_base"
       print *,"im,FSI_flag_local=",im,FSI_flag_local
       stop
       fort_FSI_flag_valid_base=0
      endif

      return
      end function fort_FSI_flag_valid_base


      function fort_FSI_flag_valid(im)
      use probcommon_module
      IMPLICIT NONE
      integer fort_FSI_flag_valid
      integer, INTENT(in) :: im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16"
       stop
      endif

      fort_FSI_flag_valid=fort_FSI_flag_valid_base(FSI_flag(im),im)

      return
      end function fort_FSI_flag_valid


      function fort_is_ice_base(FSI_flag_local,im) &
      bind(c,name='fort_is_ice_base')
      use probcommon_module

      IMPLICIT NONE

      integer fort_is_ice_base
      integer, INTENT(in) :: FSI_flag_local
      integer, INTENT(in) :: im ! 1<=im<=num_materials
      integer dummy_input

      if ((im.lt.1).or.(im.gt.num_materials)) then

       print *,"im invalid16_a in fort_is_ice_base: im=",im
       print *,"num_materials=",num_materials

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:14083"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read(*,*) dummy_input
       stop

      endif

      fort_is_ice_base=0
      if ((FSI_flag_local.eq.FSI_ICE_PROBF90).or. &
          (FSI_flag_local.eq.FSI_ICE_STATIC).or. &
          (FSI_flag_local.eq.FSI_ICE_EULERIAN_ELASTIC).or. &
          (FSI_flag_local.eq.FSI_ICE_NODES_INIT)) then
       fort_is_ice_base=1
      else if ((FSI_flag_local.eq.FSI_FLUID).or. &
               (FSI_flag_local.eq.FSI_FLUID_NODES_INIT).or. &
               (FSI_flag_local.eq.FSI_PRESCRIBED_PROBF90).or. &
               (FSI_flag_local.eq.FSI_PRESCRIBED_NODES).or. &
               (FSI_flag_local.eq.FSI_RIGID_NOTPRESCRIBED).or. &
               (FSI_flag_local.eq.FSI_EULERIAN_ELASTIC).or. &
               (FSI_flag_local.eq.FSI_SHOELE_CTML)) then
       fort_is_ice_base=0
      else
       print *,"FSI_flag_local invalid in fort_is_ice_base"
       print *,"im,FSI_flag_local ",im,FSI_flag_local
       stop
      endif

      return
      end function fort_is_ice_base

      function is_ice(im)
      use probcommon_module

      IMPLICIT NONE

      integer is_ice
      integer, INTENT(in) :: im
      integer dummy_input

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16 is_ice: im=",im

       print *,"num_materials=",num_materials

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:14132"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read(*,*) dummy_input
       stop
      endif

      is_ice=fort_is_ice_base(FSI_flag(im),im)

      return
      end function is_ice

      function fort_is_FSI_rigid_base(FSI_flag_local,im) &
      bind(c,name='fort_is_FSI_rigid_base')
      use probcommon_module

      IMPLICIT NONE

      integer fort_is_FSI_rigid_base
      integer, INTENT(in) :: FSI_flag_local
      integer, INTENT(in) :: im ! 1<=im<=num_materials

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16 fort_is_FSI_rigid_base: im:",im
       stop
      endif
      fort_is_FSI_rigid_base=0
      if (FSI_flag_local.eq.FSI_RIGID_NOTPRESCRIBED) then
       fort_is_FSI_rigid_base=1
      else if ((FSI_flag_local.eq.FSI_FLUID).or. &
               (FSI_flag_local.eq.FSI_FLUID_NODES_INIT).or. &
               (FSI_flag_local.eq.FSI_PRESCRIBED_PROBF90).or. &
               (FSI_flag_local.eq.FSI_PRESCRIBED_NODES).or. &
               (FSI_flag_local.eq.FSI_ICE_PROBF90).or. &
               (FSI_flag_local.eq.FSI_ICE_STATIC).or. &
               (FSI_flag_local.eq.FSI_ICE_NODES_INIT).or. &
               (FSI_flag_local.eq.FSI_EULERIAN_ELASTIC).or. &
               (FSI_flag_local.eq.FSI_ICE_EULERIAN_ELASTIC).or. &
               (FSI_flag_local.eq.FSI_SHOELE_CTML)) then
       fort_is_FSI_rigid_base=0
      else
       print *,"FSI_flag_local invalid in fort_is_FSI_rigid_base"
       print *,"im,FSI_flag_local ",im,FSI_flag_local
       stop
      endif

      return
      end function fort_is_FSI_rigid_base


      function fort_is_FSI_elastic_base(FSI_flag_local,im) &
      bind(c,name='fort_is_FSI_elastic_base')
      use probcommon_module

      IMPLICIT NONE

      integer fort_is_FSI_elastic_base
      integer, INTENT(in) :: FSI_flag_local
      integer, INTENT(in) :: im ! 1<=im<=num_materials

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16 fort_is_FSI_elastic_base: im:",im
       stop
      endif
      fort_is_FSI_elastic_base=0
      if ((FSI_flag_local.eq.FSI_EULERIAN_ELASTIC).or. &
          (FSI_flag_local.eq.FSI_ICE_EULERIAN_ELASTIC)) then
       fort_is_FSI_elastic_base=1
      else if ((FSI_flag_local.eq.FSI_FLUID).or. &
               (FSI_flag_local.eq.FSI_FLUID_NODES_INIT).or. &
               (FSI_flag_local.eq.FSI_PRESCRIBED_PROBF90).or. &
               (FSI_flag_local.eq.FSI_PRESCRIBED_NODES).or. &
               (FSI_flag_local.eq.FSI_ICE_PROBF90).or. &
               (FSI_flag_local.eq.FSI_ICE_STATIC).or. &
               (FSI_flag_local.eq.FSI_ICE_NODES_INIT).or. &
               (FSI_flag_local.eq.FSI_RIGID_NOTPRESCRIBED).or. &
               (FSI_flag_local.eq.FSI_SHOELE_CTML)) then
       fort_is_FSI_elastic_base=0
      else
       print *,"FSI_flag_local invalid in fort_is_FSI_elastic_base"
       print *,"im,FSI_flag_local ",im,FSI_flag_local
       stop
      endif

      return
      end function fort_is_FSI_elastic_base


      function is_FSI_rigid(im)
      use probcommon_module

      IMPLICIT NONE

      integer :: is_FSI_rigid
      integer, INTENT(in) :: im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16 is_FSI_rigid: im=",im
       stop
      endif
      is_FSI_rigid=fort_is_FSI_rigid_base(FSI_flag(im),im)

      return
      end function is_FSI_rigid

      function is_FSI_elastic(im)
      use probcommon_module

      IMPLICIT NONE

      integer :: is_FSI_elastic
      integer, INTENT(in) :: im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid16 is_FSI_elastic: im=",im
       stop
      endif
      is_FSI_elastic=fort_is_FSI_elastic_base(FSI_flag(im),im)

      return
      end function is_FSI_elastic


      function swap1_0(in_flag,control_flag)
      IMPLICIT NONE
      integer :: swap1_0
      integer, INTENT(in) :: in_flag,control_flag

      if (control_flag.eq.0) then
       swap1_0=in_flag
      else if (control_flag.eq.1) then
       swap1_0=1-in_flag
      else
       print *,"control_flag invalid"
       stop
      endif

      if ((swap1_0.eq.0).or.(swap1_0.eq.1)) then
       ! do nothing
      else
       print *,"swap1_0 invalid"
       stop
      endif

      return
      end function swap1_0


      function is_ice_or_FSI_rigid_material(im) &
      bind(c,name='is_ice_or_FSI_rigid_material')

      use probcommon_module

      IMPLICIT NONE

      integer :: is_ice_or_FSI_rigid_material
      integer, INTENT(in) :: im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid (is_ice_or_FSI_rigid_material): ",im
       stop
      endif

      is_ice_or_FSI_rigid_material=0
      if ((is_ice(im).eq.1).or. &
          (is_FSI_rigid(im).eq.1)) then
       is_ice_or_FSI_rigid_material=1
      else if ((is_ice(im).eq.0).and. &
               (is_FSI_rigid(im).eq.0)) then
       is_ice_or_FSI_rigid_material=0
      else
       print *,"is_ice or is_FSI_rigid bad"
       stop
      endif

      return
      end function is_ice_or_FSI_rigid_material

      function is_ice_or_FSI_rigid_material_project(im) &
      bind(c,name='is_ice_or_FSI_rigid_material_project')

      use probcommon_module

      IMPLICIT NONE

      integer :: is_ice_or_FSI_rigid_material_project
      integer, INTENT(in) :: im
      integer :: elastic_flag

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid (is_ice_or_FSI_rigid_material_project): ",im
       stop
      endif

      is_ice_or_FSI_rigid_material_project=0

      elastic_flag=is_FSI_elastic(im)

      if (((is_ice(im).eq.1).and. &
           (elastic_flag.eq.0)).or. &
          (is_FSI_rigid(im).eq.1)) then
       is_ice_or_FSI_rigid_material_project=1
      else if (((is_ice(im).eq.0).or. &
                (elastic_flag.eq.1)).and. &
               (is_FSI_rigid(im).eq.0)) then
       is_ice_or_FSI_rigid_material_project=0
      else
       print *,"is_ice or is_FSI_rigid or elastic_flag bad"
       print *,"im,is_ice(im),elastic_flag,is_FSI_rigid(im) ", &
         im,is_ice(im),elastic_flag,is_FSI_rigid(im)
       stop
      endif

      return
      end function is_ice_or_FSI_rigid_material_project


      function fort_is_lag_part_base(FSI_flag_local,im) &
      bind(c,name='fort_is_lag_part_base')
      use probcommon_module

      IMPLICIT NONE

      integer fort_is_lag_part_base
      integer, INTENT(in) :: FSI_flag_local
      integer, INTENT(in) :: im ! 1<=im<=num_materials
      integer dummy_input

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid17 in fort_is_lag_part_base: im=",im
       print *,"num_materials=",num_materials

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:14649"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read(*,*) dummy_input
       stop
      endif

      if ((FSI_flag_local.eq.FSI_PRESCRIBED_PROBF90).or. & 
          (FSI_flag_local.eq.FSI_PRESCRIBED_NODES).or. & 
          (FSI_flag_local.eq.FSI_SHOELE_CTML).or. & 
          (FSI_flag_local.eq.FSI_ICE_NODES_INIT).or. & 
          (FSI_flag_local.eq.FSI_FLUID_NODES_INIT)) then 
       fort_is_lag_part_base=1
      else if ((FSI_flag_local.eq.FSI_FLUID).or. &
               (FSI_flag_local.eq.FSI_ICE_PROBF90).or. & 
               (FSI_flag_local.eq.FSI_ICE_STATIC).or. & 
               (FSI_flag_local.eq.FSI_EULERIAN_ELASTIC).or. & 
               (FSI_flag_local.eq.FSI_ICE_EULERIAN_ELASTIC).or. & 
               (FSI_flag_local.eq.FSI_RIGID_NOTPRESCRIBED)) then 
       fort_is_lag_part_base=0
      else
       print *,"FSI_flag_local invalid in fort_is_lag_part_base"
       print *,"im,FSI_flag_local ",im,FSI_flag_local
       stop
      endif

      return
      end function fort_is_lag_part_base


      function is_lag_part(im)
      use probcommon_module

      IMPLICIT NONE

      integer is_lag_part
      integer, INTENT(in) :: im
      integer dummy_input

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid17 in is_lag_part: im=",im
       print *,"num_materials=",num_materials

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:14247"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read(*,*) dummy_input
       stop
      endif

      is_lag_part=fort_is_lag_part_base(FSI_flag(im),im)

      return
      end function is_lag_part
     
      function fort_read_from_CAD(fsi_flag_local) &
      bind(c,name='fort_read_from_CAD')
      IMPLICIT NONE
      integer, INTENT(in) :: fsi_flag_local
      integer fort_read_from_CAD

      fort_read_from_CAD=0
      if ((fsi_flag_local.eq.FSI_PRESCRIBED_NODES).or. & 
          (fsi_flag_local.eq.FSI_SHOELE_CTML).or. & 
          (fsi_flag_local.eq.FSI_ICE_NODES_INIT).or. & 
          (fsi_flag_local.eq.FSI_FLUID_NODES_INIT)) then 
       fort_read_from_CAD=1
      else if ((fsi_flag_local.eq.FSI_FLUID).or. & 
               (fsi_flag_local.eq.FSI_PRESCRIBED_PROBF90).or. & 
               (fsi_flag_local.eq.FSI_ICE_PROBF90).or. & 
               (fsi_flag_local.eq.FSI_ICE_STATIC).or. & 
               (fsi_flag_local.eq.FSI_EULERIAN_ELASTIC).or. & 
               (fsi_flag_local.eq.FSI_ICE_EULERIAN_ELASTIC).or. & 
               (fsi_flag_local.eq.FSI_RIGID_NOTPRESCRIBED)) then 
       ! do nothing
      else
       print *,"fsi_flag_local invalid in fort_read_from_CAD"
       print *,"fsi_flag_local=",fsi_flag_local
       stop
      endif

      return
      end function fort_read_from_CAD

       ! drag_comp>=0 and drag_comp<N_DRAG_IQ
       ! drag_im>=0 and drag_im<num_materials 
       ! drag_type>=0 and drag_type<DRAG_TYPE_IQ_NEXT
      function fort_drag_IQ_type(drag_comp,drag_im) &
      bind(c,name='fort_drag_IQ_type')

      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: drag_comp
      integer, INTENT(out) :: drag_im
      integer fort_drag_IQ_type
      integer drag_mod

      drag_im=-1
      fort_drag_IQ_type=-1

      if (num_materials.ge.2) then
       ! do nothing
      else
       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:12838"
       print *,"expecting num_materials>=2 num_materials=",num_materials
       stop
      endif

      if ((drag_comp.ge.DRAGCOMP_IQ_BODYFORCE).and. &
          (drag_comp.lt.DRAGCOMP_IQ_FORCE)) then
       drag_im=drag_comp/3
       drag_mod=MOD(drag_comp,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_FORCE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_PFORCE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_FORCE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_PFORCE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_VISCOUSFORCE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_PFORCE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_VISCOUSFORCE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_VISCOUS0FORCE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_VISCOUSFORCE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_VISCOUS0FORCE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_VISCOFORCE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_VISCOUS0FORCE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_VISCOFORCE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_BODYTORQUE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_VISCOFORCE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_BODYTORQUE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_TORQUE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_BODYTORQUE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_TORQUE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_PTORQUE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_TORQUE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_PTORQUE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_VISCOUSTORQUE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_PTORQUE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_VISCOUSTORQUE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_VISCOUS0TORQUE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_VISCOUSTORQUE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_VISCOUS0TORQUE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_VISCOTORQUE)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_VISCOUS0TORQUE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_VISCOTORQUE).and. &
               (drag_comp.lt.DRAGCOMP_IQ_COM)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_VISCOTORQUE
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_COM).and. &
               (drag_comp.lt.DRAGCOMP_IQ_MOMINERTIA)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_COM
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_MOMINERTIA).and. &
               (drag_comp.lt.DRAGCOMP_IQ_MASS)) then
       drag_mod=drag_comp-DRAGCOMP_IQ_MOMINERTIA
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_UFORCE
       else if (drag_mod.eq.1) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_VFORCE
       else if (drag_mod.eq.2) then
        fort_drag_IQ_type=DRAG_TYPE_IQ_WFORCE
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_IQ_MASS).and. &
               (drag_comp.lt.DRAGCOMP_IQ_PERIM)) then
       drag_im=drag_comp-DRAGCOMP_IQ_MASS
       fort_drag_IQ_type=DRAG_TYPE_IQ_SCALAR
      else if ((drag_comp.ge.DRAGCOMP_IQ_PERIM).and. &
               (drag_comp.lt.N_DRAG_IQ)) then
       drag_im=drag_comp-DRAGCOMP_IQ_PERIM
       fort_drag_IQ_type=DRAG_TYPE_IQ_SCALAR
      else
       print *,"drag_comp invalid GLOBALUTIL, drag_comp=",drag_comp
       stop
      endif

      if ((drag_im.ge.0).and.(drag_im.lt.num_materials)) then
       ! do nothing
      else
       print *,"drag_im invalid"
       stop
      endif
      if ((fort_drag_IQ_type.ge.0).and. &
          (fort_drag_IQ_type.lt.DRAG_TYPE_IQ_NEXT)) then
       ! do nothing
      else
       print *,"fort_drag_IQ_type invalid"
       stop
      endif

      return
      end function fort_drag_IQ_type
      

       ! drag_comp>=0 and drag_comp<N_DRAG
       ! drag_im>=0 and drag_im<num_materials 
       ! drag_type>=0 and drag_type<DRAG_TYPE_NEXT
      function fort_drag_type(drag_comp,drag_im) &
      bind(c,name='fort_drag_type')

      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: drag_comp
      integer, INTENT(out) :: drag_im
      integer fort_drag_type
      integer drag_mod

      drag_im=-1
      fort_drag_type=-1

      if (num_materials.ge.2) then
       ! do nothing
      else
       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:12838"
       print *,"expecting num_materials>=2 num_materials=",num_materials
       stop
      endif

      if ((drag_comp.ge.DRAGCOMP_STRESS).and. &
          (drag_comp.lt.DRAGCOMP_PSTRESS)) then
       drag_im=drag_comp/6
       drag_mod=MOD(drag_comp,6)
       if (drag_mod.eq.0) then
        fort_drag_type=DRAG_TYPE_T11
       else if (drag_mod.eq.1) then
        fort_drag_type=DRAG_TYPE_T12
       else if (drag_mod.eq.2) then
        fort_drag_type=DRAG_TYPE_T22
       else if (drag_mod.eq.3) then
        fort_drag_type=DRAG_TYPE_T33
       else if (drag_mod.eq.4) then
        fort_drag_type=DRAG_TYPE_T13
       else if (drag_mod.eq.5) then
        fort_drag_type=DRAG_TYPE_T23
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_PSTRESS).and. &
               (drag_comp.lt.DRAGCOMP_VISCOUSSTRESS)) then
       drag_mod=drag_comp-DRAGCOMP_PSTRESS
       drag_im=drag_mod/6
       drag_mod=MOD(drag_mod,6)
       if (drag_mod.eq.0) then
        fort_drag_type=DRAG_TYPE_T11
       else if (drag_mod.eq.1) then
        fort_drag_type=DRAG_TYPE_T12
       else if (drag_mod.eq.2) then
        fort_drag_type=DRAG_TYPE_T22
       else if (drag_mod.eq.3) then
        fort_drag_type=DRAG_TYPE_T33
       else if (drag_mod.eq.4) then
        fort_drag_type=DRAG_TYPE_T13
       else if (drag_mod.eq.5) then
        fort_drag_type=DRAG_TYPE_T23
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_VISCOUSSTRESS).and. &
               (drag_comp.lt.DRAGCOMP_VISCOUS0STRESS)) then
       drag_mod=drag_comp-DRAGCOMP_VISCOUSSTRESS
       drag_im=drag_mod/6
       drag_mod=MOD(drag_mod,6)
       if (drag_mod.eq.0) then
        fort_drag_type=DRAG_TYPE_T11
       else if (drag_mod.eq.1) then
        fort_drag_type=DRAG_TYPE_T12
       else if (drag_mod.eq.2) then
        fort_drag_type=DRAG_TYPE_T22
       else if (drag_mod.eq.3) then
        fort_drag_type=DRAG_TYPE_T33
       else if (drag_mod.eq.4) then
        fort_drag_type=DRAG_TYPE_T13
       else if (drag_mod.eq.5) then
        fort_drag_type=DRAG_TYPE_T23
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_VISCOUS0STRESS).and. &
               (drag_comp.lt.DRAGCOMP_VISCOSTRESS)) then
       drag_mod=drag_comp-DRAGCOMP_VISCOUS0STRESS
       drag_im=drag_mod/6
       drag_mod=MOD(drag_mod,6)

       if (drag_mod.eq.0) then
        fort_drag_type=DRAG_TYPE_T11
       else if (drag_mod.eq.1) then
        fort_drag_type=DRAG_TYPE_T12
       else if (drag_mod.eq.2) then
        fort_drag_type=DRAG_TYPE_T22
       else if (drag_mod.eq.3) then
        fort_drag_type=DRAG_TYPE_T33
       else if (drag_mod.eq.4) then
        fort_drag_type=DRAG_TYPE_T13
       else if (drag_mod.eq.5) then
        fort_drag_type=DRAG_TYPE_T23
       else
        print *,"drag_mod invalid"
        stop
       endif

      else if ((drag_comp.ge.DRAGCOMP_VISCOSTRESS).and. &
               (drag_comp.lt.DRAGCOMP_TORQUE_ARM)) then
       drag_mod=drag_comp-DRAGCOMP_VISCOSTRESS
       drag_im=drag_mod/6
       drag_mod=MOD(drag_mod,6)

       if (drag_mod.eq.0) then
        fort_drag_type=DRAG_TYPE_T11
       else if (drag_mod.eq.1) then
        fort_drag_type=DRAG_TYPE_T12
       else if (drag_mod.eq.2) then
        fort_drag_type=DRAG_TYPE_T22
       else if (drag_mod.eq.3) then
        fort_drag_type=DRAG_TYPE_T33
       else if (drag_mod.eq.4) then
        fort_drag_type=DRAG_TYPE_T13
       else if (drag_mod.eq.5) then
        fort_drag_type=DRAG_TYPE_T23
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_TORQUE_ARM).and. &
               (drag_comp.lt.DRAGCOMP_FLAG)) then
       drag_mod=drag_comp-DRAGCOMP_TORQUE_ARM
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_type=DRAG_TYPE_UVEC
       else if (drag_mod.eq.1) then
        fort_drag_type=DRAG_TYPE_VVEC
       else if (drag_mod.eq.2) then
        fort_drag_type=DRAG_TYPE_WVEC
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_VISCOUS0STRESSMAG).and. &
               (drag_comp.lt.DRAGCOMP_VISCOUS0STRESSTAN)) then
       drag_mod=drag_comp-DRAGCOMP_VISCOUS0STRESSMAG
       drag_im=drag_mod
       fort_drag_type=DRAG_TYPE_SCALAR
      else if ((drag_comp.ge.DRAGCOMP_VISCOUS0STRESSTAN).and. &
               (drag_comp.lt.N_DRAG)) then
       drag_mod=drag_comp-DRAGCOMP_VISCOUS0STRESSTAN
       drag_im=drag_mod/3
       drag_mod=MOD(drag_mod,3)
       if (drag_mod.eq.0) then
        fort_drag_type=DRAG_TYPE_UVEC
       else if (drag_mod.eq.1) then
        fort_drag_type=DRAG_TYPE_VVEC
       else if (drag_mod.eq.2) then
        fort_drag_type=DRAG_TYPE_WVEC
       else
        print *,"drag_mod invalid"
        stop
       endif
      else if ((drag_comp.ge.DRAGCOMP_FLAG).and. &
               (drag_comp.lt.DRAGCOMP_VISCOUS0STRESSMAG)) then
       drag_mod=drag_comp-DRAGCOMP_FLAG
       drag_im=drag_mod
       fort_drag_type=DRAG_TYPE_FLAG
      else
       print *,"drag_comp invalid GLOBALUTIL, drag_comp=",drag_comp
       stop
      endif

      if ((drag_im.ge.0).and.(drag_im.lt.num_materials)) then
       ! do nothing
      else
       print *,"drag_im invalid: ",drag_im
       stop
      endif
      if ((fort_drag_type.ge.0).and. &
          (fort_drag_type.lt.DRAG_TYPE_NEXT)) then
       ! do nothing
      else
       print *,"fort_drag_type invalid: ",fort_drag_type
       stop
      endif

      return
      end function fort_drag_type

      function fort_is_passive_advect_test() &
      bind(c,name='fort_is_passive_advect_test')
      use probcommon_module

      IMPLICIT NONE
      integer fort_is_passive_advect_test
      real(amrex_real) dx_placeholder(SDIM)
      real(amrex_real) x_placeholder(SDIM)
      real(amrex_real) time_placeholder
      real(amrex_real) LS
      real(amrex_real) vel(SDIM)
      real(amrex_real) temperature
      integer prescribed_flag
      integer dir

      do dir=1,SDIM
       dx_placeholder(dir)=one
       x_placeholder(dir)=zero
      enddo
      time_placeholder=zero
      call SUB_clamped_LS_no_scale(x_placeholder,time_placeholder,LS,vel, &
         temperature,prescribed_flag,dx_placeholder)
      if (LS.eq.CLAMPED_EVERYWHERE_LS) then
       fort_is_passive_advect_test=1
      else if (LS.ne.CLAMPED_EVERYWHERE_LS) then
       fort_is_passive_advect_test=0
      else
       print *,"LS is NaN"
       stop
       fort_is_passive_advect_test=0
      endif

      return
      end function fort_is_passive_advect_test

      function fort_is_rigid_base(FSI_flag_local,im) &
      bind(c,name='fort_is_rigid_base')
      use probcommon_module

      IMPLICIT NONE

      integer fort_is_rigid_base
      integer, INTENT(in) :: FSI_flag_local
      integer, INTENT(in) :: im
      integer dummy_input

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid17 in fort_is_rigid_base: im=",im
       print *,"num_materials=",num_materials

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:10822"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read(*,*) dummy_input
       stop
      endif

      if ((FSI_flag_local.eq.FSI_PRESCRIBED_PROBF90).or. & 
          (FSI_flag_local.eq.FSI_PRESCRIBED_NODES).or. & 
          (FSI_flag_local.eq.FSI_SHOELE_CTML)) then 
       fort_is_rigid_base=1  ! non-tessellating material
      else if ((FSI_flag_local.eq.FSI_FLUID).or. &
               (FSI_flag_local.eq.FSI_FLUID_NODES_INIT).or. & 
               (FSI_flag_local.eq.FSI_ICE_PROBF90).or. & 
               (FSI_flag_local.eq.FSI_ICE_STATIC).or. & 
               (FSI_flag_local.eq.FSI_ICE_NODES_INIT).or. & 
               (FSI_flag_local.eq.FSI_EULERIAN_ELASTIC).or. & 
               (FSI_flag_local.eq.FSI_ICE_EULERIAN_ELASTIC).or. & 
               (FSI_flag_local.eq.FSI_RIGID_NOTPRESCRIBED)) then 
       fort_is_rigid_base=0  ! tessellating material
      else
       print *,"FSI_flag_local invalid in fort_is_rigid_base"
       print *,"FSI_flag_local=",FSI_flag_local
       print *,"im=",im
       stop
       fort_is_rigid_base=0  ! prevent compiler warnings
      endif

      return
      end function fort_is_rigid_base

      function is_rigid(im) 
      use probcommon_module

      IMPLICIT NONE

      integer is_rigid
      integer, INTENT(in) :: im
      integer dummy_input

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid17 in is_rigid: im=",im
       print *,"num_materials=",num_materials

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:10864"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read(*,*) dummy_input
       stop
      endif

      is_rigid=fort_is_rigid_base(FSI_flag(im),im)

      return
      end function is_rigid

      function is_rigid_CL(im)
      use probcommon_module
      IMPLICIT NONE

      integer is_rigid_CL
      integer, INTENT(in) :: im
      integer dummy_input

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid in is_rigid_CL: im=",im
       print *,"num_materials=",num_materials

       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:15446"
       print *,"By pressing <CTRL C> during this read statement, the"
       print *,"gdb debugger will produce a stacktrace."
       print *,"type 0 then <enter> to exit the program"

       read(*,*) dummy_input
       stop
      endif

      if ((is_rigid(im).eq.1).or. &
          (is_FSI_elastic(im).eq.1).or. &
          (is_ice_or_FSI_rigid_material(im).eq.1)) then
       is_rigid_CL=1
      else if ((is_rigid(im).eq.0).and. &
               (is_FSI_elastic(im).eq.0).and. &
               (is_ice_or_FSI_rigid_material(im).eq.0)) then
       is_rigid_CL=0
      else
       print *,"is_rigid, FSI_elastic,or is_ice_or_FSI_rigid_material invalid"
       print *,"im,is_rigid,is_FSI_elastic,is_ice_or_FSI_rigid_material ", &
        im,is_rigid(im),is_FSI_elastic(im),is_ice_or_FSI_rigid_material(im)
       stop
      endif

      return
      end function is_rigid_CL

      function fort_built_in_elastic_model(elastic_visc_in, &
        viscoelastic_model_in) &
      bind(c,name='fort_built_in_elastic_model')
      IMPLICIT NONE

      integer fort_built_in_elastic_model
      real(amrex_real), INTENT(in) :: elastic_visc_in
      integer, INTENT(in) :: viscoelastic_model_in

      if (elastic_visc_in.gt.zero) then
       if ((viscoelastic_model_in.eq.NN_FENE_CR).or. & ! FENE-CR
           (viscoelastic_model_in.eq.NN_OLDROYD_B).or. & ! Oldroyd B
           (viscoelastic_model_in.eq.NN_MAIRE_ABGRALL_ETAL).or. & !incremental
           (viscoelastic_model_in.eq.NN_NEO_HOOKEAN).or. & ! incremental 
           (viscoelastic_model_in.eq.NN_FENE_P).or. & ! FENE-P
           (viscoelastic_model_in.eq.NN_LINEAR_PTT)) then ! linear PTT
        fort_built_in_elastic_model=1
       else
        print *,"viscoelastic_model_in invalid: ",viscoelastic_model_in
        stop
        fort_built_in_elastic_model=0
       endif
      else if (elastic_visc_in.eq.zero) then
       fort_built_in_elastic_model=0
      else
       print *,"elastic_visc_in invalid: ",elastic_visc_in
       stop
       fort_built_in_elastic_model=0
      endif

      return
      end function fort_built_in_elastic_model

      function rigid_exists()
      use probcommon_module

      IMPLICIT NONE

      integer rigid_exists
      integer im

      rigid_exists=0

      do im=1,num_materials
       if (is_rigid(im).eq.1) then
        rigid_exists=1
       else if (is_rigid(im).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid GLOBALUTIL.F90 in rigid_exists"
        stop
       endif
      enddo ! im

      return
      end function rigid_exists

       ! solid_dist>0 in the solid
      subroutine combine_solid_LS(LS,solid_dist,local_im_solid_primary)
      use probcommon_module

      IMPLICIT NONE

      real(amrex_real), INTENT(out) :: solid_dist
      real(amrex_real), INTENT(in) :: LS(num_materials)
      integer im
      integer, INTENT(out) :: local_im_solid_primary

      solid_dist=-99999.0
      local_im_solid_primary=0

      do im=1,num_materials
       if (is_rigid(im).eq.0) then
        ! do nothing
       else if (is_rigid(im).eq.1) then
        if (local_im_solid_primary.eq.0) then
         solid_dist=LS(im)
         local_im_solid_primary=im
        else if ((local_im_solid_primary.ge.1).and. &
                 (local_im_solid_primary.le.num_materials)) then
         if (LS(im).gt.solid_dist) then
          solid_dist=LS(im)
          local_im_solid_primary=im
         else if (LS(im).le.solid_dist) then
          ! do nothing
         else
          print *,"LS(im) bust"
          stop
         endif
        else
         print *,"local_im_solid_primary invalid: ",local_im_solid_primary
         stop
        endif 
       else
        print *,"is_rigid invalid GLOBALUTIL.F90"
        stop
       endif
      enddo ! im=1..num_materials

      end subroutine combine_solid_LS


      subroutine combine_solid_VOF(VOF,solid_vof,local_im_solid_primary)
      use probcommon_module

      IMPLICIT NONE

      real(amrex_real), INTENT(out) :: solid_vof
      real(amrex_real), INTENT(in) :: VOF(num_materials)
      integer, INTENT(out) :: local_im_solid_primary
      integer im

      solid_vof=zero
      local_im_solid_primary=0

      do im=1,num_materials
       if (is_rigid(im).eq.0) then
        ! do nothing
       else if (is_rigid(im).eq.1) then
        solid_vof=solid_vof+VOF(im)

        if (local_im_solid_primary.eq.0) then
         local_im_solid_primary=im
        else if ((local_im_solid_primary.ge.1).and. &
                 (local_im_solid_primary.le.num_materials)) then
         if (VOF(im).gt.VOF(local_im_solid_primary)) then
          local_im_solid_primary=im
         endif
        else
         print *,"local_im_solid_primary invalid: ",local_im_solid_primary
         stop
        endif

       else
        print *,"is_rigid invalid GLOBALUTIL.F90"
        stop
       endif
      enddo ! im=1..num_materials

      end subroutine combine_solid_VOF

      function rigid_count()
      use probcommon_module

      IMPLICIT NONE

      integer rigid_count,im

      rigid_count=0
      do im=1,num_materials
       if (is_rigid(im).eq.1) then
        rigid_count=rigid_count+1
       else if (is_rigid(im).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid GLOBALUTIL.F90"
        stop
       endif
      enddo ! im

      return
      end function rigid_count

      function is_in_probtype_list()
      use probcommon_module
      IMPLICIT NONE

      integer is_in_probtype_list
      integer iprob

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

      integer im_solid_primary
      integer im

      im_solid_primary=0
      do im=1,num_materials
       if ((is_rigid(im).eq.1).and. &
           (im_solid_primary.eq.0)) then 
        im_solid_primary=im
       else if ((is_rigid(im).eq.0).or. &
                (im_solid_primary.gt.0)) then
        ! do nothing
       else
        print *,"is_rigid invalid GLOBALUTIL.F90"
        stop
       endif
      enddo ! im

      return
      end function im_solid_primary

      subroutine debug_im_solid(im)
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: im
   
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid19"
       stop
      endif 
      if (is_rigid(im).eq.1) then
       ! do nothing
      else if (is_rigid(im).eq.0) then
       print *,"expecting solid im=",im
       stop
      else
       print *,"is_rigid invalid GLOBALUTIL.F90"
       stop
      endif

      return 
      end subroutine debug_im_solid

      subroutine get_LS_extend(LS,iten,LS_extend)
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: iten
      real(amrex_real), INTENT(in) :: LS(num_materials)
      real(amrex_real), INTENT(out) :: LS_extend
      integer im,im_opp

      call get_inverse_iten(im,im_opp,iten)
      if (im.ge.im_opp) then
       print *,"im or im_opp invalid"
       stop
      endif
      if ((LS(im).gt.zero).and. &
          (LS(im_opp).gt.zero)) then
       print *,"cannot have LS(im)>0 and LS(im_opp)>0"
       print *,"im,im_opp,LS(im),LS(im_opp) ", &
         im,im_opp,LS(im),LS(im_opp)
       print *,"is_rigid(im) ",is_rigid(im)
       print *,"is_rigid(im_opp) ",is_rigid(im_opp)
       stop
      else if (LS(im).gt.LS(im_opp)) then
       LS_extend=-LS(im_opp)
      else if (LS(im_opp).gt.LS(im)) then
       LS_extend=LS(im)
      else if (LS(im).eq.LS(im_opp)) then
       LS_extend=zero
      else
       print *,"LS bust in get_LS_extend"
       print *,"im,im_opp,LS(im),LS(im_opp) ", &
         im,im_opp,LS(im),LS(im_opp)
       print *,"is_rigid(im) ",is_rigid(im)
       print *,"is_rigid(im_opp) ",is_rigid(im_opp)
       stop
      endif
      
      return
      end subroutine get_LS_extend

      subroutine get_LSNRM_extend(LS,NRM,iten,NRM_extend)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: iten
      real(amrex_real), INTENT(in) :: LS(num_materials)
      real(amrex_real), INTENT(in) :: NRM(num_materials*SDIM)
      real(amrex_real), INTENT(out) :: NRM_extend(SDIM)
      integer im,im_opp,dir

      call get_inverse_iten(im,im_opp,iten)
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
       print *,"LS bust in get_LSNRM_extend: ",im,im_opp,LS(im),LS(im_opp)
       stop
      endif
      
      return
      end subroutine get_LSNRM_extend

      subroutine get_VOF_extend(VOF,iten,VOF_extend)
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: iten
      real(amrex_real), INTENT(in) :: VOF(num_materials)
      real(amrex_real), INTENT(out) :: VOF_extend
      integer im,im_opp

      call get_inverse_iten(im,im_opp,iten)
      if (im.ge.im_opp) then
       print *,"im or im_opp invalid"
       stop
      endif
      if ((VOF(im).ge.-EPS1).and. &
          (VOF(im).le.one+EPS1).and. &
          (VOF(im_opp).ge.-EPS1).and. &
          (VOF(im_opp).le.one+EPS1)) then 

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
       print *,"get_VOF_extend"
       print *,"VOF bust im,im_opp: ",im,im_opp
       print *,"VOF(im),VOF(im_opp) ",VOF(im),VOF(im_opp)
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

      real(amrex_real) gradT(3)
      real(amrex_real) facearea_minus
      real(amrex_real) facearea_plus
      real(amrex_real) facearea_total
      real(amrex_real) xminus(3)
      real(amrex_real) xplus(3)
      real(amrex_real) Tminus,Tplus
      real(amrex_real) nminus(3)
      real(amrex_real) xminusI(3)
      real(amrex_real) nplus(3)
      real(amrex_real) xplusI(3)

      real(amrex_real) mag_I,mag1_I,mag2_I
      real(amrex_real) mag_ntilde,mag1_tilde,mag2_tilde
      real(amrex_real) mag_nplus
      real(amrex_real) mag_test
      integer dir,dircrit_tilde,dircrit_I
      integer i,j,k
      real(amrex_real) s_I(3)
      real(amrex_real) s_tilde(3)
      real(amrex_real) ntilde(3)
      real(amrex_real) t1_tilde(3)
      real(amrex_real) t2_tilde(3)
      real(amrex_real) t1_I(3)
      real(amrex_real) t2_I(3)
      real(amrex_real), INTENT(in) :: x2_I_test(3)
      real(amrex_real), INTENT(in) :: x3_I_test(3)
      real(amrex_real) x2_I(3)
      real(amrex_real) x3_I(3)
      real(amrex_real) y(4,3) ! y(j,dir)  j=1..4 dir=1..3
      real(amrex_real) basis_vec(3,3) ! basis_vec(k,dir) k=1..3 dir=1..3
       ! T2I,T3I,TminusI,TplusI
      real(amrex_real), INTENT(in) :: Tdata(4) ! Tdata(j) j=1..4
       ! T_{-}+((T_{+}-T_{-})/mag_ntilde)((y_{j}-x_{-}) dot ntilde)
      real(amrex_real) Tstar(4) ! Tstar(j) j=1..4 
      real(amrex_real) wt_epsilon
      real(amrex_real) wt(4)
      real(amrex_real) basis2D(4,4)  ! basis2D(i,j)=g_{i}(y_{j})
      real(amrex_real) AA(2,2) ! AA(k,i)=sum_j w_j g_{i}(y_{j})g_{k}(y_{j})
      real(amrex_real) BB(2)   ! BB(k)=sum_j w_j g_{k}(y_{j})(T_j - T^{*}(y_{j}))
      real(amrex_real) coeffs(2)
      real(amrex_real) fixed_coeff
      integer matstatus

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

      integer, INTENT(in) :: n
      integer :: i
      real(amrex_real), INTENT(in) :: l(n),u(n),d(n),f(n)
      real(amrex_real), INTENT(out) :: soln(n)
      real(amrex_real) :: ll(n),uu(n),dd(n),z(n)

      dd(1)=d(1)

      if (dd(1).ne.zero) then
       ! do nothing
      else
       print *,"dd(1) invalid ",dd(1)
       stop
      endif

      uu(1)=u(1)/dd(1)
      z(1)=f(1)/dd(1)

      do i=2,n-1

       ll(i)=l(i)
       dd(i)=d(i)-ll(i)*uu(i-1)

       if (dd(i).ne.zero) then
        ! do nothing
       else
        print *,"dd(i) invalid ",dd(i)
        stop
       endif

       uu(i)=u(i)/dd(i)
       z(i)=(f(i)-ll(i)*z(i-1))/dd(i)
      enddo
      ll(n)=l(n)
      dd(n)=d(n)-ll(n)*uu(n-1)

      if (dd(n).ne.zero) then
       ! do nothing
      else
       print *,"dd(n) invalid ",dd(n)
       stop
      endif

      z(n)=(f(n)-ll(n)*z(n-1))/dd(n)
      soln(n)=z(n)
      do i=n-1,1,-1
       soln(i)=z(i)-uu(i)*soln(i+1)
      enddo

      return
      end subroutine tridiag_solve
  
       ! dir=0..sdim-1
      real(amrex_real) function xphys_of_xcomp(dir,xcomp)
      use probcommon_module
      IMPLICIT NONE
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: xcomp
      integer :: nelement
      integer :: icrit
      real(amrex_real) :: xlo,xhi,comp_dx
      real(amrex_real) :: xcell_lo,xcell_hi
      real(amrex_real) :: xphys_lo,xphys_hi
 
      nelement=mapping_n_cell(dir)
      if (nelement.ge.1) then
       ! do nothing
      else
       print *,"nelement invalid xphys_of_xcomp nelement=",nelement
       stop
      endif
      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir invalid: ",dir
       stop
      endif
      xlo=problo_array(dir+1)
      xhi=probhi_array(dir+1)
      comp_dx=(xhi-xlo)/nelement
      if (comp_dx.gt.zero) then
       ! do nothing
      else
       print *,"comp_dx invalid: ",comp_dx
       stop
      endif

      if (use_identity_mapping.eq.1) then
       xphys_of_xcomp=xcomp
      else if (use_identity_mapping.eq.0) then

       if (xcomp.le.xlo) then
        xphys_of_xcomp=xcomp
       else if (xcomp.ge.xhi) then
        xphys_of_xcomp=xcomp
       else if ((xcomp.gt.xlo).and. &
                (xcomp.lt.xhi)) then
         ! xcomp=xlo+(i+1/2)*dx
         ! i=NINT((xcomp-xlo)/dx-1/2)
        icrit=NINT((xcomp-xlo)/comp_dx-half)
        if ((icrit.ge.0).and.(icrit.le.nelement-1)) then
         xcell_lo=xlo+icrit*comp_dx
         xcell_hi=xcell_lo+comp_dx
         xphys_lo=mapping_comp_to_phys(icrit,dir)
         xphys_hi=mapping_comp_to_phys(icrit+1,dir)

         if (xphys_hi.gt.xphys_lo) then
          ! do nothing
         else
          print *,"xphys_hi-xphys_lo invalid: ",xphys_lo,xphys_hi
          stop
         endif

         if (abs(xcomp-xcell_lo).le.comp_dx*EPS_12_4) then
          xphys_of_xcomp=xphys_lo
         else if (abs(xcomp-xcell_hi).le.comp_dx*EPS_12_4) then
          xphys_of_xcomp=xphys_hi
         else if ((xcomp.gt.xcell_lo-comp_dx*EPS_12_2).and. &
                  (xcomp.le.xcell_lo)) then
          xphys_of_xcomp=xphys_lo
         else if ((xcomp.lt.xcell_hi+comp_dx*EPS_12_2).and. &
                  (xcomp.ge.xcell_hi)) then
          xphys_of_xcomp=xphys_hi
         else if ((xcomp.ge.xcell_lo).and. &
                  (xcomp.le.xcell_hi)) then
          xphys_of_xcomp=xphys_lo+ &
             (xphys_hi-xphys_lo)* &
             (xcomp-xcell_lo)/(xcell_hi-xcell_lo)
         else
          print *,"xcomp invalid: ",xcomp,xcell_lo,xcell_hi,comp_dx
          stop
         endif
        else if (icrit.lt.0) then
         xphys_of_xcomp=xcomp
        else if (icrit.ge.nelement) then
         xphys_of_xcomp=xcomp
        else
         print *,"icrit invalid"
         stop
        endif
       else
        print *,"xcomp is NaN"
        stop
       endif

      else 
       print *,"use_identity_mapping invalid:",use_identity_mapping
       stop
      endif

      end function xphys_of_xcomp

       !dir=0..sdim-1
      real(amrex_real) function xcomp_of_xphys(dir,xphys)
      use probcommon_module
      IMPLICIT NONE
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: xphys
      integer :: nelement
      integer :: icrit
      real(amrex_real) :: xlo,xhi,comp_dx
      real(amrex_real) :: xcell_lo,xcell_hi
      real(amrex_real) :: xcomp_lo,xcomp_hi
 
      nelement=mapping_n_cell(dir)
      if (nelement.ge.1) then
       ! do nothing
      else
       print *,"nelement invalid xcomp_of_xphys nelement=",nelement
       stop
      endif
      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir invalid: ",dir
       stop
      endif
      xlo=problo_array(dir+1)
      xhi=probhi_array(dir+1)
      comp_dx=(xhi-xlo)/nelement
      if (comp_dx.gt.zero) then
       ! do nothing
      else
       print *,"comp_dx invalid: ",comp_dx
       stop
      endif

      if (use_identity_mapping.eq.1) then
       xcomp_of_xphys=xphys
      else if (use_identity_mapping.eq.0) then

       if (xphys.le.xlo) then
        xcomp_of_xphys=xphys
       else if (xphys.ge.xhi) then
        xcomp_of_xphys=xphys
       else if ((xphys.gt.xlo).and. &
                (xphys.lt.xhi)) then
         ! xphys=xlo+(i+1/2)*dx
         ! i=NINT((xphys-xlo)/dx-1/2)
        icrit=NINT((xphys-xlo)/comp_dx-half)
        if ((icrit.ge.0).and.(icrit.le.nelement-1)) then
         xcell_lo=xlo+icrit*comp_dx
         xcell_hi=xcell_lo+comp_dx
         xcomp_lo=mapping_phys_to_comp(icrit,dir)
         xcomp_hi=mapping_phys_to_comp(icrit+1,dir)

         if (xcomp_hi.gt.xcomp_lo) then
          ! do nothing
         else
          print *,"xcomp_hi-xcomp_lo invalid: ",xcomp_lo,xcomp_hi
          stop
         endif

         if (abs(xphys-xcell_lo).le.EPS_12_4*comp_dx) then
          xcomp_of_xphys=xcomp_lo
         else if (abs(xphys-xcell_hi).le.EPS_12_4*comp_dx) then
          xcomp_of_xphys=xcomp_hi
         else if ((xphys.gt.xcell_lo-comp_dx*EPS_12_2).and. &
                  (xphys.le.xcell_lo)) then
          xcomp_of_xphys=xcomp_lo
         else if ((xphys.lt.xcell_hi+comp_dx*EPS_12_2).and. &
                  (xphys.ge.xcell_hi)) then
          xcomp_of_xphys=xcomp_hi
         else if ((xphys.ge.xcell_lo).and. &
                  (xphys.le.xcell_hi)) then
          xcomp_of_xphys=xcomp_lo+ &
             (xcomp_hi-xcomp_lo)* &
             (xphys-xcell_lo)/(xcell_hi-xcell_lo)
         else
          print *,"xphys invalid: ",xphys
          print *,"xcell_lo=",xcell_lo
          print *,"xcell_hi=",xcell_hi
          print *,"xcomp_lo=",xcomp_lo
          print *,"xcomp_hi=",xcomp_hi
          print *,"comp_dx=",comp_dx
          stop
         endif
        else if (icrit.lt.0) then
         xcomp_of_xphys=xphys
        else if (icrit.ge.nelement) then
         xcomp_of_xphys=xphys
        else
         print *,"icrit invalid"
         stop
        endif
       else
        print *,"xphys is NaN"
        stop
       endif

      else 
       print *,"use_identity_mapping invalid:",use_identity_mapping
       stop
      endif
      
      end function xcomp_of_xphys

       ! dir=0..sdim-1
       ! solve x(X)-xstar=0 for Xstar using Newton's method.
      subroutine inverse_mapping(phys_coord,comp_coord,dir)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) :: xlo,xhi
      real(amrex_real), INTENT(in) :: phys_coord 
      real(amrex_real), INTENT(inout) :: comp_coord 
      integer :: nelement
      integer, INTENT(in) :: dir
      real(amrex_real) :: conv_err
      integer :: conv_iter
      integer, PARAMETER :: conv_iter_max=100
      real(amrex_real) :: comp_dx
      real(amrex_real) :: comp_lo,comp_hi
      real(amrex_real) :: fprime,ff
      real(amrex_real) :: comp_coord_new

      nelement=mapping_n_cell(dir)
      xlo=problo_array(dir+1)
      xhi=probhi_array(dir+1)
      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir invalid"
       stop
      endif
      comp_dx=(xhi-xlo)/nelement
      if (comp_dx.gt.zero) then
       ! do nothing
      else
       print *,"comp_dx invalid: ",comp_dx
       stop
      endif
      if (nelement.ge.1) then
       ! do nothing
      else
       print *,"nelement invalid inverse_mapping nelement=",nelement
       stop
      endif
      conv_iter=0
       ! solve x(X)-xstar=0 for Xstar using Newton's method.
      do while (conv_iter.lt.conv_iter_max-1)
       comp_hi=comp_coord+comp_dx*half
       comp_lo=comp_coord-comp_dx*half
       fprime=(xphys_of_xcomp(dir,comp_hi)- &
               xphys_of_xcomp(dir,comp_lo))/comp_dx
       if (fprime.gt.zero) then
        ! do nothing
       else
        print *,"fprime must be positive"
        stop
       endif
 
       ff=xphys_of_xcomp(dir,comp_coord)-phys_coord
       comp_coord_new=comp_coord-ff/fprime 
 
       conv_err=abs(comp_coord-comp_coord_new)
       comp_coord=comp_coord_new
       conv_iter=conv_iter+1
       if ((conv_iter.ge.1).and. &
           (conv_iter.le.conv_iter_max)) then
        ! do nothing
       else
        print *,"conv_iter out of range"
        stop
       endif

      enddo ! while (conv_iter<conv_iter_max-1)

      end subroutine inverse_mapping

       ! dir=0,1,2
      subroutine single_dimension_grid_mapping(dir)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: dir
      real(amrex_real) :: xlo,xhi
      integer :: nelement

      real(amrex_real) :: comp_coord
      real(amrex_real) :: phys_coord

      real(amrex_real), dimension(:), allocatable :: phys_coord_new
      real(amrex_real), dimension(:), allocatable :: wt_coord
      real(amrex_real), dimension(:), allocatable :: tri_l
      real(amrex_real), dimension(:), allocatable :: tri_u
      real(amrex_real), dimension(:), allocatable :: tri_d
      real(amrex_real), dimension(:), allocatable :: tri_f
      real(amrex_real), dimension(:), allocatable :: tri_soln

      real(amrex_real) :: conv_err
      integer :: conv_iter
      integer, PARAMETER :: conv_iter_max=100
      real(amrex_real) :: comp_dx
      real(amrex_real) :: local_phys
      real(amrex_real) :: local_wt
      integer :: nsolve
      integer :: i

      nelement=mapping_n_cell(dir)
      xlo=problo_array(dir+1)
      xhi=probhi_array(dir+1)

      allocate(phys_coord_new(0:nelement))
      allocate(wt_coord(0:nelement-1))

      nsolve=nelement-1

      allocate(tri_l(1:nsolve))
      allocate(tri_u(1:nsolve))
      allocate(tri_d(1:nsolve))
      allocate(tri_f(1:nsolve))
      allocate(tri_soln(1:nsolve))

      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir invalid"
       stop
      endif
      if (nelement.ge.1) then
       ! do nothing
      else
       print *,"nelement invalid single_dimension_grid_mapping nelement=", &
               nelement
       stop
      endif

      comp_dx=(xhi-xlo)/nelement
      if (comp_dx.gt.zero) then
       ! do nothing
      else
       print *,"comp_dx invalid: ",comp_dx
       stop
      endif

      do i=0,nelement
       comp_coord=xlo+i*comp_dx
       mapping_comp_to_phys(i,dir)=comp_coord
      enddo
       ! get rid of floating point round-off err
      mapping_comp_to_phys(nelement,dir)=xhi 

      conv_iter=0

       ! x(X) maps computational grid to physical grid.
       ! solve: div_X (1/w(x)) grad_X x=0  x(Xlo)=Xlo   x(Xhi)=Xhi
       !
      do while (conv_iter.lt.conv_iter_max-1)

       do i=0,nelement-1 
        local_phys=half*(mapping_comp_to_phys(i,dir)+ &
                         mapping_comp_to_phys(i+1,dir))
         ! returns (1/w) where w>>1 in "trouble" regions
        call SUB_MAPPING_WEIGHT_COEFF(dir,local_wt,local_phys)
        if (local_wt.gt.zero) then
         wt_coord(i)=local_wt/comp_dx
        else 
         print *,"local_wt invalid"
         stop
        endif
       enddo !do i=0,nelement-1 

       nsolve=nelement-1

       do i=1,nsolve
        tri_l(i)=-wt_coord(i-1)
        tri_u(i)=-wt_coord(i)
        tri_d(i)=-(tri_l(i)+tri_u(i))
        tri_f(i)=zero
       enddo 
       tri_f(1)=-tri_l(1)*xlo
       tri_f(nsolve)=-tri_u(nsolve)*xhi

       call tridiag_solve(tri_l,tri_u,tri_d,nsolve,tri_f,tri_soln)

       phys_coord_new(0)=xlo
       phys_coord_new(nelement)=xhi

       conv_err=zero

       do i=1,nsolve
        phys_coord_new(i)=tri_soln(i)
        conv_err=conv_err+ &
          (mapping_comp_to_phys(i,dir)- &
           phys_coord_new(i))**2
        mapping_comp_to_phys(i,dir)=phys_coord_new(i)
       enddo ! do i=1,nsolve

       do i=0,nelement-1
        if (phys_coord_new(i+1)-phys_coord_new(i).gt.zero) then
         ! do nothing
        else
         print *,"phys_coord_new not monotonic"
         stop
        endif
       enddo !i=0,nelement-1

       conv_err=sqrt(conv_err/nelement)
       conv_iter=conv_iter+1
       if ((conv_iter.ge.1).and. &
           (conv_iter.le.conv_iter_max)) then
        ! do nothing
       else
        print *,"conv_iter out of range"
        stop
       endif
       
      enddo ! do while conv_iter<conv_iter_max-1

      deallocate(tri_l)
      deallocate(tri_u)
      deallocate(tri_d)
      deallocate(tri_f)
      deallocate(tri_soln)

      deallocate(wt_coord)
      deallocate(phys_coord_new)

      do i=0,nelement
       phys_coord=xlo+i*comp_dx
       mapping_phys_to_comp(i,dir)=phys_coord
      enddo
       ! get rid of floating point round-off err
      mapping_phys_to_comp(nelement,dir)=xhi 

      do i=1,nelement-1
       phys_coord=xlo+i*comp_dx
        ! initial guess.
       comp_coord=mapping_phys_to_comp(i-1,dir)
        ! solve x(X)=xstar
       call inverse_mapping(phys_coord,comp_coord,dir)
       mapping_phys_to_comp(i,dir)=comp_coord
      enddo ! do i=1,nelement-1
 
      do i=0,nelement-1
       if (mapping_phys_to_comp(i+1,dir)- &
           mapping_phys_to_comp(i,dir).gt.zero) then
        ! do nothing
       else
        print *,"mapping_phys_to_comp not monotonic"
        stop
       endif
      enddo !i=0,nelement-1

      return
      end subroutine single_dimension_grid_mapping

      subroutine patterned_substrates(x,y,z,dist,time,im_substrate, &
                      ptb_dist_low,ptb_dist_high)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x,y,z,time
      real(amrex_real), INTENT(in) :: ptb_dist_low,ptb_dist_high
      integer, INTENT(in) :: im_substrate
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real) :: xprime
      real(amrex_real) :: yprime
      real(amrex_real) :: zprime
      real(amrex_real) :: xvec(SDIM)
      real(amrex_real) :: local_pi
      real(amrex_real) :: pitch,ptb_f,ptb_disbtx,ptb_disbty,ptb_dist,rPillar

      if ((im_substrate.lt.1).or. &
          (im_substrate.gt.num_materials)) then
       print *,"im_substrate invalid10"
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
       print *,"time bust in patterned dist: ",time
       stop
      endif

      if (is_rigid(im_substrate).ne.1) then
       print *,"is_rigid invalid GLOBALUTIL.F90"
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
         !< Was 0.0005 to make a dimple surface, now turned it
         !to 0.001 to make all flat surface
         ptb_dist=ptb_dist_low 
        else
         ptb_dist=ptb_dist_high
        endif
      else
       ptb_dist=ptb_dist_high
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

! dist>0 outside of cylinder
      subroutine cylinderdist(x,y,z,xcen,ycen,rad,zmin,zmax,dist)
      use probcommon_module
      IMPLICIT NONE
    
      real(amrex_real), INTENT(in) :: x,y,z,xcen,ycen,rad
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real) :: dist2
      real(amrex_real), INTENT(in) :: zmin,zmax

      if (zmin.ge.zmax-EPS10) then 
       print *,"invalid parameters ",zmin,zmax
       stop
      endif
      dist=sqrt((x-xcen)**2+(y-ycen)**2)-rad
      if (z.ge.zmax) then
       if (dist.le.zero) then
        dist=z-zmax
       else
        dist=sqrt(dist**2+(z-zmax)**2)
       endif
      else if (z.le.zmin) then
       if (dist.le.zero) then
        dist=zmin-z
       else
        dist=sqrt(dist**2+(zmin-z)**2)
       endif
      else if ((z.ge.zmin).and.(z.le.zmax)) then
       if (dist.ge.zero) then
        !do nothing
       else if (dist.le.zero) then
        dist2=max(zmin-z,z-zmax)
        dist=max(dist2,dist)
       else
        print *,"dist invalid ",dist
        stop
       endif
      else
       print *,"z,zmin or zmax invalid: ",z,zmin,zmax
       stop
      endif

      return 
      end subroutine cylinderdist


! dist>0 outside of annulus
      subroutine annulusdist(x,y,z,xcen,ycen,radmean,radthick,zmin,zmax,dist)
      use probcommon_module
      IMPLICIT NONE
    
      real(amrex_real), INTENT(in) :: x,y,z,xcen,ycen,radmean,radthick
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real) :: dist2
      real(amrex_real), INTENT(in) :: zmin,zmax

      if (zmin.ge.zmax-EPS10) then 
       print *,"invalid parameters ",zmin,zmax
       stop
      endif
      dist=abs(sqrt((x-xcen)**2+(y-ycen)**2)-radmean)-radthick
      if (z.ge.zmax) then
       if (dist.le.zero) then
        dist=z-zmax
       else
        dist=sqrt(dist**2+(z-zmax)**2)
       endif
      else if (z.le.zmin) then
       if (dist.le.zero) then
        dist=zmin-z
       else
        dist=sqrt(dist**2+(zmin-z)**2)
       endif
      else if ((z.ge.zmin).and.(z.le.zmax)) then
       if (dist.ge.zero) then
        !do nothing
       else if (dist.le.zero) then
        dist2=max(zmin-z,z-zmax)
        dist=max(dist2,dist)
       else
        print *,"dist invalid ",dist
        stop
       endif
      else
       print *,"z,zmin or zmax invalid: ",z,zmin,zmax
       stop
      endif

      return 
      end subroutine annulusdist


       ! negative on the inside of the square
      subroutine squaredist(x,y,xlo,xhi,ylo,yhi,dist)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x,y,xlo,xhi,ylo,yhi
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real) dist1
      real(amrex_real) xmid,ymid
 
      if ((xlo.ge.xhi-EPS10).or.(ylo.ge.yhi-EPS10)) then 
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

      real(amrex_real), INTENT(in) :: xmin,xmax,ymin,ymax,zmin,zmax
      real(amrex_real), INTENT(in) :: x,y,z
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real) xcen,ycen,zcen,xrad,yrad,zrad
      real(amrex_real) xdist,ydist,zdist

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

      subroutine dumpstring_headers(plot_sdim,add_sub_cells)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: plot_sdim
      integer, INTENT(in) :: add_sub_cells
      character*80 rootname
      character*80 Varname
      character*2 specstr
      character*2 refinestr
      character*2 matstr
      character*2 matstropp
      integer igrad
      integer im,imls,im_opp
      integer ih,ih_root
      integer ispec
      integer irefine
      integer dir,i
      integer nparts,nparts_def,partid
      integer plot_sdim_macro
      integer test_nwrite

      if ((plot_sdim.ne.2).and.(plot_sdim.ne.3)) then
       print *,"plot_sdim invalid"
       stop
      endif
      if (plot_sdim.ge.SDIM) then
       ! do nothing
      else
       print *,"expecting plot_sdim>=sdim in dumpstring_headers"
       stop
      endif
      plot_sdim_macro=plot_sdim

      nparts=0
      do im=1,num_materials
       if (is_lag_part(im).eq.1) then
        nparts=nparts+1
       else if (is_lag_part(im).eq.0) then
        ! do nothing
       else
        print *,"is_lag_part(im) invalid dumpstring_headers"
        stop
       endif
      enddo !im=1..num_materials

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid dumpstring_headers"
       stop
      endif
      nparts_def=nparts
      if (nparts_def.eq.0) then
       nparts_def=1
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid dumpstring_headers2"
       stop
      endif

      if ((add_sub_cells.eq.2).or. &
          (add_sub_cells.eq.3)) then
       Varname='x_sub'
       call dumpstring(Varname)
       Varname='y_sub'
       call dumpstring(Varname)
       if (add_sub_cells.eq.3) then
        Varname='z_sub'
        call dumpstring(Varname)
       endif
      else if (add_sub_cells.eq.0) then
       ! do nothing
      else
       print *,"add_sub_cells invalid"
       stop
      endif

      test_nwrite=0

      Varname='X'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

      Varname='Y'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

      if (plot_sdim.eq.3) then
       Varname='Z'
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1
      endif

      Varname='x_velocity' 
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1
      Varname='y_velocity'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

      if (plot_sdim.eq.3) then
       Varname='z_velocity'
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1
      endif

      if (test_nwrite.eq.PLOTCOMP_PRES) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_PRES)"
       stop
      endif

       ! multigrid pressure  "PRES_MG"
      Varname='PRES_MG'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

       ! EOS pressure: "PRES_EOS"
      Varname='PRES_EOS'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

       ! Divergence derived from the velocity: "DIV"
       ! see MacProj.cpp: NavierStokes::getStateDIV_ALL
      Varname='DIV_DERIVED'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

       ! expected divergence: "DIV_EXPECT"
       ! see NavierStokes.cpp: NavierStokes::getStateDIV_DATA
       ! "DIV_Type"
      Varname='DIV_EXPECT'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

      Varname='MACH'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

      if (test_nwrite.eq.PLOTCOMP_VFRAC) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_VFRAC)"
       stop
      endif

       !VFRACS 
      do im=1,num_materials

       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       ih=1
       Varname='F'
       ih=ih+1
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1
      enddo  ! im (volume fractions)

      if (test_nwrite.eq.PLOTCOMP_LS) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_LS)"
       stop
      endif

        ! levelset
      do imls=1,num_materials
       im=imls
       im_opp=imls

       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo
       write(matstropp,'(I2)') im_opp
       do i=1,2
        if (matstropp(i:i).eq.' ') then
         matstropp(i:i)='0'
        endif
       enddo

       ih=1
       Varname='L'
       ih=ih+1
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       do i=1,2
        Varname(ih:ih)=matstropp(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1
      enddo  ! imls (levelset variables)

       ! levelset normals
      do imls=1,num_materials
       im=imls
       im_opp=imls

       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo
       write(matstropp,'(I2)') im_opp
       do i=1,2
        if (matstropp(i:i).eq.' ') then
         matstropp(i:i)='0'
        endif
       enddo

       do dir=1,plot_sdim
        ih=1
        if (dir.eq.1) then
         Varname='x_normal'
        else if (dir.eq.2) then
         Varname='y_normal'
        else if ((dir.eq.3).and.(plot_sdim.eq.3)) then
         Varname='z_normal'
        else
         print *,"dir invalid dumpstring_headers"
         stop
        endif
        ih=9

        do i=1,2
         Varname(ih:ih)=matstr(i:i)
         ih=ih+1
        enddo
        do i=1,2
         Varname(ih:ih)=matstropp(i:i)
         ih=ih+1
        enddo
        call dumpstring(Varname)
        test_nwrite=test_nwrite+1
       enddo  ! dir=1..plot_sdim
      enddo  ! imls (levelset normal variables)

      if (test_nwrite.eq.PLOTCOMP_SCALARS) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_SCALARS)"
       stop
      endif

       ! density, temperature, mass fractions
      do im=1,num_materials

       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

        ! density
       ih=1
       Varname='D'
       ih=ih+1
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

        ! temperature
       ih=1
       Varname='T'
       ih=ih+1
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

       do ispec=1,num_species_var

        write(specstr,'(I2)') ispec
        do i=1,2
         if (specstr(i:i).eq.' ') then
          specstr(i:i)='0'
         endif
        enddo

        ih=1
        Varname='S'
        ih=ih+1
        do i=1,2
         Varname(ih:ih)=specstr(i:i)
         ih=ih+1
        enddo
        Varname(ih:ih)='-'
        ih=ih+1
        do i=1,2
         Varname(ih:ih)=matstr(i:i)
         ih=ih+1
        enddo
        call dumpstring(Varname)
        test_nwrite=test_nwrite+1
   
       enddo  ! ispec=1..num_species_var

      enddo  ! im (state variables); do im=1,num_materials

      if (test_nwrite.eq.PLOTCOMP_SCALARS_MERGE) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_SCALARS_MERGE)"
       stop
      endif

      rootname='MERGE'

       ! density
      ih=1
      Varname='D'
      ih=ih+1
      do ih_root=1,5
       Varname(ih:ih)=rootname(ih_root:ih_root)
       ih=ih+1
      enddo
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

       ! temperature
      ih=1
      Varname='T'
      ih=ih+1
      do ih_root=1,5
       Varname(ih:ih)=rootname(ih_root:ih_root)
       ih=ih+1
      enddo
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

      do ispec=1,num_species_var

       write(specstr,'(I2)') ispec
       do i=1,2
        if (specstr(i:i).eq.' ') then
         specstr(i:i)='0'
        endif
       enddo

       ih=1
       Varname='S'
       ih=ih+1

       do ih_root=1,5
        Varname(ih:ih)=rootname(ih_root:ih_root)
        ih=ih+1
       enddo

       Varname(ih:ih)='-'
       ih=ih+1
       do i=1,2
        Varname(ih:ih)=specstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1
   
      enddo  ! ispec=1..num_species_var

      if (test_nwrite.eq.PLOTCOMP_MOMDEN) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_MOMDEN)"
       stop
      endif

       ! mom_density
      do im=1,num_materials

       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

        ! mom_density
       ih=1
       Varname='MOMDEN'
       ih=ih+6
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

      enddo  ! im=1..num_materials mom_den

      do partid=1,num_materials_viscoelastic
       im=fort_im_viscoelastic_map(partid)+1
       if ((im.ge.1).and.(im.le.num_materials)) then
        write(matstr,'(I2)') im
        do i=1,2
         if (matstr(i:i).eq.' ') then
          matstr(i:i)='0'
         endif
        enddo

        do ispec=1,ENUM_NUM_TENSOR_TYPE

         write(specstr,'(I2)') ispec
         do i=1,2
          if (specstr(i:i).eq.' ') then
           specstr(i:i)='0'
          endif
         enddo

         do irefine=1,ENUM_NUM_REFINE_DENSITY_TYPE

          write(refinestr,'(I2)') irefine
          do i=1,2
           if (refinestr(i:i).eq.' ') then
            refinestr(i:i)='0'
           endif
          enddo

          ih=1
          Varname='CT-ijk-'
          ih=ih+7

          do i=1,2
           Varname(ih:ih)=refinestr(i:i)
           ih=ih+1
          enddo
          Varname(ih:ih)='-'
          ih=ih+1

          do i=1,2
           Varname(ih:ih)=specstr(i:i)
           ih=ih+1
          enddo
          Varname(ih:ih)='-'
          ih=ih+1
          Varname(ih:ih)='i'
          ih=ih+1
          Varname(ih:ih)='m'
          ih=ih+1
          Varname(ih:ih)='-'
          ih=ih+1

          do i=1,2
           Varname(ih:ih)=matstr(i:i)
           ih=ih+1
          enddo
          call dumpstring(Varname)
          test_nwrite=test_nwrite+1

         enddo  ! irefine=1..ENUM_NUM_REFINE_DENSITY_TYPE

        enddo  ! ispec=1..num_tensor_type
       else
        print *,"im invalid60"
        stop
       endif
      enddo ! partid=1..num_materials_viscoelastic 

      do partid=1,num_materials_compressible

       im=fort_im_refine_density_map(partid)+1
       if ((im.ge.1).and.(im.le.num_materials)) then
        write(matstr,'(I2)') im
        do i=1,2
         if (matstr(i:i).eq.' ') then
          matstr(i:i)='0'
         endif
        enddo

        do ispec=1,ENUM_NUM_REFINE_DENSITY_TYPE

         write(specstr,'(I2)') ispec
         do i=1,2
          if (specstr(i:i).eq.' ') then
           specstr(i:i)='0'
          endif
         enddo

         ih=1
         Varname='REFINEDEN'
         ih=ih+9
         do i=1,2
          Varname(ih:ih)=specstr(i:i)
          ih=ih+1
         enddo
         Varname(ih:ih)='-'
         ih=ih+1
         do i=1,2
          Varname(ih:ih)=matstr(i:i)
          ih=ih+1
         enddo
         call dumpstring(Varname)
         test_nwrite=test_nwrite+1

        enddo  ! ispec=1..num_refine_density_type
       else
        print *,"im invalid (REFINEDEN varname)"
        stop
       endif
      enddo ! partid=1..num_materials_compressible

      if (test_nwrite.eq.PLOTCOMP_VISC) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_VISC)"
       stop
      endif

       ! viscosity
      do im=1,num_materials
       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       ih=1
       Varname='MU'
       ih=ih+2
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1
      enddo  ! im (viscosity variables)

       ! thermal conductivity
      do im=1,num_materials
       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       ih=1
       Varname='K_THERMAL'
       ih=ih+9
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1
      enddo  ! im (thermal conductivity variables)

      if (test_nwrite.eq.PLOTCOMP_TRACE_A_VORT) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_TRACE_A_VORT)"
       stop
      endif

       ! gamma_dot, TR(A), TR(A)*shear thinning factor, TR(A)*thin*f(A),
       ! vorticity
      do im=1,num_materials

       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       ih=1
       Varname='DT'  ! sqrt(2 * D:D)
       ih=ih+2
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

       ih=1
       Varname='TR' ! Tr(A)
       ih=ih+2
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

       ih=1
       Varname='TRT' ! Tr(A) * (mu_L - eta_S)/eta_P if FENE-CR+Carreau
       ih=ih+3
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

       ih=1
       Varname='TRTF'  ! TRT * f(A)
       ih=ih+4
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

       ih=1
       Varname='VORT'
       ih=ih+4
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

      enddo  ! im (trace variables)

      if (test_nwrite.eq.PLOTCOMP_F_ELASTIC_X) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_F_ELASTIC_X)"
       stop
      endif

      Varname='x_ELSTCFORCE'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1
      Varname='y_ELSTCFORCE'
      call dumpstring(Varname)
      test_nwrite=test_nwrite+1

      if (plot_sdim.eq.3) then
       Varname='z_ELSTCFORCE'
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1
      endif

      if (test_nwrite.eq.PLOTCOMP_GRAD_VELOCITY) then
       ! do nothing
      else
       print *,"(test_nwrite.ne.PLOTCOMP_GRAD_VELOCITY)"
       stop
      endif

      do igrad=1,AMREX_SPACEDIM_SQR

       write(matstr,'(I2)') igrad
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       !ux,vx,wx,uy,vy,wy,uz,vz,wz
       ih=1
       Varname='GRADVEL'
       ih=ih+7
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
       test_nwrite=test_nwrite+1

      enddo  ! igrad=1,AMREX_SPACEDIM_SQR

      if (test_nwrite.eq.PLOTCOMP_NCOMP) then
       ! do nothing
      else
       print *,"test_nwrite != PLOTCOMP_NCOMP"
       stop
      endif

      return
      end subroutine dumpstring_headers


      subroutine dumpstring_headers_sanity(plot_sdim,ncomp)
      IMPLICIT NONE

      integer, INTENT(in) :: plot_sdim
      integer, INTENT(in) :: ncomp
      character*80 Varname
      character*3 matstr
      integer ih,im,i

      if ((plot_sdim.ne.2).and.(plot_sdim.ne.3)) then
       print *,"plot_sdim invalid"
       stop
      endif
      if (ncomp.ge.1) then
       ! do nothing
      else
       print *,"ncomp invalid in dumpstring_headers_sanity"
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
       tower_mf_id, &
       nsteps, &
       num_levels, &
       time, &
       visual_revolve, &
       ncomp)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: n_root
      character, dimension(n_root), INTENT(in) :: root_char_array
      integer, INTENT(in) :: data_dir
      integer, INTENT(in) :: plot_sdim
      integer :: plot_sdim_macro
      integer klo_plot,khi_plot

      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: total_number_grids
      integer, INTENT(in) :: num_levels
      integer, INTENT(in) :: grids_per_level_array(num_levels)
      integer, INTENT(in) :: levels_array(total_number_grids)
      integer, INTENT(in) :: bfact_array(total_number_grids)
      integer, INTENT(in) :: gridno_array(total_number_grids)
      integer, INTENT(in) :: gridlo_array(total_number_grids*SDIM)
      integer, INTENT(in) :: gridhi_array(total_number_grids*SDIM)
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: SDC_outer_sweeps
      integer, INTENT(in) :: slab_step
      integer, INTENT(in) :: tower_mf_id
      integer, INTENT(in) :: nsteps
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: visual_revolve

      integer strandid

      integer nwrite3d,nwrite2d,index3d,index2d

      character*3 levstr
      character*5 gridstr
      character*32 filename32
      character*80 rmcommand

      character(len=plotfile_digits) :: stepstr
      character*3 outerstr
      character*3 slabstr
      character*3 idstr

      character(len=n_root) :: root_char_str
      character(len=n_root+30+plotfile_digits) :: newfilename40
      character(len=30+plotfile_digits) :: fname_extend
      character(len=4) :: step_chars
      character(len=2) :: dir_chars
      character(len=5) :: outer_chars
      character(len=4) :: slab_chars
      character(len=2) :: id_chars
      character(len=4) :: plt_chars

      integer i,j,k,dir
      integer ilev,igrid
      integer lo(plot_sdim),hi(plot_sdim)
      integer sysret
      integer hi_index_shift(3)

! Guibo

      character*80 Title,Zonename
      real(tecplot_real_short) :: ZONEMARKER,EOHMARKER
      integer*4 :: nzones_gb,iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb
      integer bfact,testlev,testgridno
      integer testlo(SDIM),testhi(SDIM)
      integer add_sub_cells

      ! define zone structure
      type zone3d_t
         real(tecplot_real), pointer :: var(:,:,:,:)
      end type zone3d_t
      type(zone3d_t), dimension(:), allocatable :: zone3d_gb

      type zone2d_t
         real(tecplot_real), pointer :: var(:,:,:)
      end type zone2d_t
      type(zone2d_t), dimension(:), allocatable :: zone2d_gb

      real(amrex_real) theta,rr,zz,xx,yy

      plot_sdim_macro=SDIM


      if (tower_mf_id.ge.0) then
       !do nothing
      else if (tower_mf_id-GET_NEW_DATA_OFFSET.ge.0) then
       !do nothing
      else
       print *,"tower_mf_id out of range"
       stop
      endif

      if (plot_sdim.ne.3) then
       print *,"plot_sdim invalid"
       stop
      endif

      if (SDIM.ne.2) then
       print *,"must be called only in 2D"
       stop
      endif

      if (levelrz.ne.COORDSYS_RZ) then
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
       print *,"ncomp invalid zones_revolve_sanity"
       stop
      endif

      nwrite2d=SDIM+ncomp
      nwrite3d=plot_sdim+ncomp

      if (tower_mf_id.eq.MULTIFAB_TOWER_PLT_MF) then
       if (ncomp.eq.PLOTCOMP_NCOMP) then
        ! do nothing
       else
        print *,"ncomp.ne.PLOTCOMP_NCOMP zones_revolve_sanity"
        stop
       endif
      endif

      if (num_levels.ne.finest_level+1) then
       print *,"num_levels invalid"
       stop
      endif

      write(stepstr,18200) nsteps
18200 format(I plotfile_digits)

      do i=1,plotfile_digits
       if (stepstr(i:i).eq.' ') then
        stepstr(i:i)='0'
       endif
      enddo
      write(outerstr,'(I3)') SDC_outer_sweeps
      write(slabstr,'(I3)') slab_step

      if (tower_mf_id.ge.0) then
       write(idstr,'(I3)') tower_mf_id
      else if (tower_mf_id-GET_NEW_DATA_OFFSET.ge.0) then
       write(idstr,'(I3)') tower_mf_id-GET_NEW_DATA_OFFSET
      else
       print *,"tower_mf_id invalid"
       stop
      endif

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
      write(fname_extend,18224) &
               dir_chars,id_chars,idstr, &
               step_chars,stepstr,outer_chars,outerstr, &
               slab_chars,slabstr,plt_chars
18224 format(A2,A2,A3,A4,A plotfile_digits ,A5,A3,A4,A3,A4)

      newfilename40(n_root+1:n_root+30+plotfile_digits)=fname_extend

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
       write(filename32,'(A14,A10,A3,A5)') &
              './temptecplot/','tempnddata',levstr,gridstr
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

      if (tower_mf_id.eq.MULTIFAB_TOWER_PLT_MF) then
       add_sub_cells=plot_sdim
       call dumpstring_headers(plot_sdim_macro,add_sub_cells)
      else
       ! Variable names: zones_revolve_sanity
       call dumpstring_headers_sanity(plot_sdim,ncomp)
      endif

       ! Zones
      do iz_gb=1,nzones_gb
       ! Zone marker. Value = 299.0
       write(11) ZONEMARKER
       ! Zone name
       Zonename = "ZONE"
       call dumpstring(Zonename)

       strandid=1

       write(11) -1   ! Parent Zone
       write(11) strandid-1   ! StrandID (this does not work)
       write(11) round_time(time) ! Solution time
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

       write(filename32,'(A14,A10,A3,A5)') &
              './temptecplot/','tempnddata',levstr,gridstr
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

      return
      end subroutine zones_revolve_sanity


       ! normal is contact line normal pointing towards "im" material
       ! (material 0 liquid)
       ! vel_n is the advancing velocity if vel_n>0 and
       ! the receding velocity if vel_n<0.
       ! cos_thetae is the cosine of the static angle inbetween
       ! material 0 and material 2 (the solid material)
       ! cos_thetad is the output cosine of the dynamic angle.
       ! vis is viscosity of material 0 (liquid)
       ! imodel=0 static 
       ! imodel=1 Jiang   
       ! imodel=2 Kistler (obsolete: use Zeyu's code, DCA_select=105,
       ! i.e. ZEYU_DCA_SELECT=5)
       ! DCA_select_model was written by Yongsheng Lian and
       ! students of his.
      subroutine DCA_select_model(normal,vel_n,cos_thetae,vis, &
        user_tension_scalar, &
        cos_thetad,imodel)
      implicit none

      integer, intent(in) :: imodel !0=static  1=Jiang  2=Kistler
      real(amrex_real), intent(in) :: normal(SDIM)
      real(amrex_real), intent(in) :: vel_n
      real(amrex_real), intent(in) :: cos_thetae  
      real(amrex_real), intent(out) :: cos_thetad  
      real(amrex_real), intent(in) :: vis 
      real(amrex_real), intent(in) :: user_tension_scalar
      real(amrex_real) capillary,f_Hoff_inver,temp

      if (user_tension_scalar.gt.zero) then
       ! do nothing
      else
       print *,"user_tension_scalar invalid: ",user_tension_scalar
       stop
      endif
      if (vis.ge.zero) then
       ! do nothing
      else
       print *,"vis invalid: ",vis
       stop
      endif

      f_Hoff_inver = 0.0276
       ! vel_n is the advancing velocity if vel_n>0 and
       ! the receding velocity if vel_n<0.
      capillary = abs(vel_n)*vis/user_tension_scalar

      if (capillary.eq.zero) then
       cos_thetad=cos_thetae !static
      else if (capillary.ne.zero) then

       if (imodel.eq.0) then ! static
        cos_thetad=cos_thetae
       else if (imodel.eq.1) then ! Jiang's model
        temp = capillary**0.702d0
        temp = tanh(4.96d0*temp)
        temp = temp*(1.0d0+cos_thetae)
        if (vel_n.gt.zero) then
          !advancing => cos_thetad<cos_thetae => thetad>thetae
         temp = cos_thetae-temp
        else if (vel_n.lt.zero) then
          !receding => cos_thetad>cos_thetae => thetad<thetae
         temp = cos_thetae+temp
        else
         print *,"expecting vel_n<>0: ",vel_n
         stop
        endif
        if (temp.gt.1.0) temp = 1.0  
        if (temp.lt.-1.0) temp = -1.0  
        cos_thetad=temp
       else if (imodel.eq.2) then ! Kistler's model
        print *,"use Zeyu's code instead for Kistler!"
        print *,"make use_DCA=105 for Kistler (i.e. ZEYU_DCA_SELECT=5)"
        stop
       else
        print*,"dynamic contact angle model type not valid: ",imodel
        stop
       end if

      else
       print *,"capillary invalid: ",capillary
       stop
      endif

      return
      end subroutine DCA_select_model



       ! Dynamic Contact Angle
      subroutine get_use_DCA(use_DCA) &
      bind(c,name='get_use_DCA')

      use probcommon_module

      IMPLICIT NONE

      integer, intent(out) :: use_DCA
      integer im_solid_DCA

      im_solid_DCA=im_solid_primary()

      if ((probtype.eq.5501).and.(SDIM.eq.3)) then
       if (im_solid_DCA.ne.3) then
        print *,"im_solid_DCA invalid"
        stop
       endif
        ! xblob10=0 static
        ! xblob10=1 Jiang
        ! xblob10=2 Kistler
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

      real(amrex_real) rho_in,T_in
      ! From Caudwell et al., Int. J. Thermophysics (25) 5, 2004
      real(amrex_real) rho,T,eta,visccoef
      real(amrex_real) rmin,rmax
      real(amrex_real) T2,V0,V,VV,VV2,eta_star

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

      real(amrex_real), INTENT(in) :: local_plus,local_minus
      real(amrex_real), INTENT(in) :: wt_plus,wt_minus
      real(amrex_real), INTENT(out) :: coeff

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
      real(amrex_real) function get_user_viscconst(im,density,temperature)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: im
      real(amrex_real), INTENT(in) :: density,temperature
      real(amrex_real) mu

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im out of range"
       stop
      endif

      if (temperature.gt.zero) then
       ! do nothing
      else
       print *,"temperature invalid"
       stop
      endif
      if (density.gt.zero) then
       ! do nothing
      else
       print *,"density invalid in get_user_viscconst"
       stop
      endif

      if ((fort_viscconst(im).ge.zero).and. &
          (fort_prerecalesce_viscconst(im).ge.zero)) then
       ! do nothing
      else
       print *,"fortran viscconst invalid"
       stop
      endif

      get_user_viscconst=fort_viscconst(im)

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
      real(amrex_real) function get_user_heatviscconst(im)
      use probcommon_module
      IMPLICIT NONE

      integer im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im out of range"
       stop
      endif

      if ((fort_heatviscconst(im).ge.zero).and. &
          (fort_prerecalesce_heatviscconst(im).ge.zero)) then
       ! do nothing
      else
       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break GLOBALUTIL.F90:18205"
       print *,"fortran heatviscconst invalid"
       stop
      endif

      get_user_heatviscconst=fort_heatviscconst(im)

      return
      end function get_user_heatviscconst


        ! CONTAINER ROUTINE FOR MEHDI VAHAB
      real(amrex_real) function get_user_stiffCP(im)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: im

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im out of range"
       stop
      endif
      if ((fort_stiffCP(im).ge.zero).and. &
          (fort_prerecalesce_stiffCP(im).ge.zero)) then
       ! do nothing
      else
       print *,"fortran stiffCP invalid"
       stop
      endif
      get_user_stiffCP=fort_stiffCP(im)

      return
      end function get_user_stiffCP

       ! 1<=iten<=2 * num_interfaces
       ! default_flag=1 => only the sign is needed
       ! default_flag=0 => the value is important too.
      real(amrex_real) function get_user_latent_heat(iten,temperature, &
                default_flag) &
      bind(c,name='get_user_latent_heat')

      use probcommon_module
      IMPLICIT NONE

      integer, value, INTENT(in) :: iten
      real(amrex_real), value, INTENT(in) :: temperature
      integer, value, INTENT(in) :: default_flag
      real(amrex_real) new_latent_heat
      real(amrex_real) sign_latent_heat

      if ((iten.ge.1).and.(iten.le.2*num_interfaces)) then
       ! do nothing
      else
       print *,"iten invalid"
       stop
      endif

      if ((default_flag.eq.0).or.(default_flag.eq.1)) then
       ! do nothing
      else
       print *,"default_flag invalid"
       stop
      endif

      if (fort_latent_heat(iten).ge.zero) then
       sign_latent_heat=one
      else if (fort_latent_heat(iten).lt.zero) then
       sign_latent_heat=-one
      else
       print *,"fort_latent_heat invalid"
       stop
      endif

      new_latent_heat=abs(fort_latent_heat(iten))
      if (fort_latent_heat_slope(iten).eq.zero) then
       ! do nothing
      else if (fort_latent_heat_slope(iten).lt.zero) then
       if (default_flag.eq.1) then
        ! do nothing
       else if (default_flag.eq.0) then
        new_latent_heat=abs(fort_latent_heat(iten))+ &
           fort_latent_heat_slope(iten)* &
           (temperature-fort_latent_heat_T0(iten))

        if (is_in_probtype_list().eq.1) then
         call SUB_VARIABLE_LATENT_HEAT(iten,temperature,new_latent_heat)
        endif

        if (new_latent_heat.lt.abs(fort_latent_heat_min(iten))) then
         new_latent_heat=abs(fort_latent_heat_min(iten))
        else if (new_latent_heat.ge. &
                 abs(fort_latent_heat_min(iten))) then
         ! do nothing
        else
         print *,"new_latent_heat or fort_latent_heat_min NaN"
         stop
        endif
        if (new_latent_heat.ge.zero) then
         ! do nothing
        else
         print *,"new_latent_heat invalid"
         stop
        endif
       else
        print *,"default_flag invalid"
        stop
       endif
      else
       print *,"fort_latent_heat_slope(iten) should be non-positive"
       stop
      endif

      get_user_latent_heat=sign_latent_heat*new_latent_heat

      return
      end function get_user_latent_heat


      subroutine get_user_tension(xpos,time, &
        tension,new_tension,temperature)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: xpos(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: temperature(num_materials)
      real(amrex_real), INTENT(in) :: tension(num_interfaces)
      real(amrex_real), INTENT(out) :: new_tension(num_interfaces)
      real(amrex_real) avgtemp
      integer iten,im,im_opp

       ! fort_tension
       ! fort_prefreeze_tension
      do iten=1,num_interfaces
       new_tension(iten)=tension(iten)
       if (fort_tension_min(iten).ge.zero) then
        ! do nothing
       else
        print *,"fort_tension_min invalid"
        stop
       endif

        ! im<im_opp
       call get_inverse_iten(im,im_opp,iten)

       if (fort_tension_slope(iten).eq.zero) then

        if ((temperature(im).gt.zero).and. &
            (temperature(im_opp).gt.zero)) then
         ! do nothing
        else
         print *,"temperature must be positive"
         print *,"im,im_opp ",im,im_opp
         print *,"temperature(im) ",temperature(im)
         print *,"temperature(im_opp) ",temperature(im_opp)
         stop
        endif
        avgtemp=half*(temperature(im)+temperature(im_opp))

        if (is_in_probtype_list().eq.1) then
         call SUB_VARIABLE_SURFACE_TENSION(xpos,time,iten, &
           avgtemp,new_tension(iten))
        endif

        if (new_tension(iten).ge.zero) then
         ! do nothing
        else
         print *,"new_tension invalid: ",iten,new_tension(iten)
         stop
        endif

       else if (fort_tension_slope(iten).lt.zero) then

        if (fort_tension_T0(iten).gt.zero) then
         ! do nothing
        else
         print *,"T0 invalid"
         stop
        endif
        if ((temperature(im).gt.zero).and. &
            (temperature(im_opp).gt.zero)) then
         ! do nothing
        else
         print *,"temperature must be positive"
         stop
        endif
        avgtemp=half*(temperature(im)+temperature(im_opp))
        new_tension(iten)=new_tension(iten)+fort_tension_slope(iten)* &
         (avgtemp-fort_tension_T0(iten))

        if (is_in_probtype_list().eq.1) then
         call SUB_VARIABLE_SURFACE_TENSION(xpos,time,iten, &
           avgtemp,new_tension(iten))
        endif

        if (new_tension(iten).lt.fort_tension_min(iten)) then
         new_tension(iten)=fort_tension_min(iten)
        else if (new_tension(iten).ge.fort_tension_min(iten)) then
         ! do nothing
        else
         print *,"new_tension or fort_tension_min NaN"
         stop
        endif
        if (new_tension(iten).ge.zero) then
         ! do nothing
        else
         print *,"new_tension invalid: ",iten,new_tension(iten)
         stop
        endif
       else
        print *,"fort_tension_slope(iten) should be non-positive"
        stop
       endif
      enddo ! iten

      return
      end subroutine get_user_tension


      subroutine TEMPERATURE_default(rho,temperature,internal_energy, &
        imattype,im)

      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: imattype
      integer, INTENT(in) :: im
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(out) :: temperature
      real(amrex_real), INTENT(in) :: internal_energy
      real(amrex_real) cv

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid70"
       stop
      endif
!      cv=4.1855D+7
      cv=get_user_stiffCP(im)
      if (cv.gt.zero) then
       ! do nothing
      else
       print *,"cv invalid in temperature default (must be positive) cv=",cv
       stop
      endif
      if (imattype.eq.999) then
       if (is_rigid(im).eq.0) then
        print *,"is_rigid invalid GLOBALUTIL.F90"
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

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(in) :: temperature
      real(amrex_real), INTENT(out) :: internal_energy
      real(amrex_real) cv


      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid in internal_default"
       stop
      endif
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

!      cv=4.1855D+7
      cv=get_user_stiffCP(im)
      if (cv.gt.zero) then
       ! do nothing
      else
       print *,"cv invalid in internal default (must be positive)"
       print *,"cv=get_user_stiffCP(im)=",cv
       print *,"im=",im
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
      subroutine ice_substrate_distance(x,y,z,angle_x,angle_y,dist)
      use probcommon_module
      IMPLICIT NONE
      
      real(amrex_real), INTENT(in) :: x,y,z,angle_x,angle_y
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real) aspect,yprime,zprime,aspect2

      if (SDIM.eq.2) then
       if (abs(y-z).le.EPS2) then
        ! do nothing
       else
        print *,"y=z in 2d expected"
        stop
       endif
      endif

      aspect=tan(angle_x)
      if (SDIM.eq.2) then
        yprime=aspect*(x-xblob2)+yblob2
        dist=y-yprime  ! vertical distance
      else if (SDIM.eq.3) then
        aspect2=tan(angle_y)
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

      real(amrex_real), INTENT(in) :: x_point(SDIM)
      real(amrex_real), INTENT(out) :: dist
      integer icomp
      real(amrex_real) distarr(n_sites)
      real(amrex_real), INTENT(in) :: nucleate_pos(4*n_sites)
      real(amrex_real) hugedist
      real(amrex_real) xx(SDIM)
      real(amrex_real) rr
      integer dir
 
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

      real(amrex_real), intent(out) :: A,B,GAMMA,R1,R2,RHOI


      if ((probtype.eq.36).and.(axis_dir.eq.2)) then  ! spherical explosion
       A=5.484D+12
       B=0.09375D+12
       R1=4.94D0
       R2=1.21D0
       GAMMA=1.28D0
       RHOI=1.63D0
      else if ((probtype.eq.36).and.(axis_dir.eq.310)) then !hydrobulge
       A=6.17D+12 !cgs
!      A=6.1327D+12 !cgs
       B=1.69D+11 !cgs 
!      B=1.5069D+11 !cgs 
       R1=4.4D0
       R2=1.2D0
       GAMMA=1.25D0
       RHOI=1.765D0 !cgs
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
       print *,"probtype not supported for the jwl material: ",probtype
       print *,"axis_dir=",axis_dir
       stop
      endif

      return
      end subroutine get_jwl_constants


      subroutine INTERNAL_jwl(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,temperature
      real(amrex_real), intent(out) :: internal_energy
      real(amrex_real) :: cv

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative INTERNAL_jwl:",rho
       stop
      endif
      if (temperature.gt.zero) then
       !do nothing
      else
       print *,"temperature <=0 INTERNAL_jwl: ",temperature
       stop
      endif

      cv=4.1855D+7
      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_jwl

      subroutine TEMPERATURE_jwl(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: temperature
      real(amrex_real) :: cv

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative TEMPERATURE_jwl:",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy <=0 TEMPERATURE_jwl:",internal_energy
       stop
      endif

      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_jwl

      subroutine ENTROPY_jwl(rho,internal_energy,entropy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: entropy
      real(amrex_real) :: pressure,press_adiabat

      call EOS_NAjwl(rho,internal_energy,pressure)
      call EOS_jwlADIABAT(rho,internal_energy,press_adiabat)
      if (press_adiabat.gt.zero) then
       !do nothing
      else
       print *,"press_adiabat invalid"
       stop
      endif
      entropy=pressure/press_adiabat

      return
      end subroutine ENTROPY_jwl

      subroutine INTERNAL_ENTROPY_jwl(rho,entropy,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,entropy,internal_energy
      real(amrex_real) A,B,R1,R2,GAMMA,RHOI,OMEGA,pressure_part,press_adiabat

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      OMEGA=GAMMA-one
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid"
       stop
      endif
      if (entropy.gt.zero) then
       !do nothing
      else
       print *,"entropy invalid"
       stop
      endif
      pressure_part= &
        A*(one-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(one-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)
      call EOS_jwlADIABAT(rho,internal_energy,press_adiabat)
      internal_energy=(press_adiabat*entropy-pressure_part)/ &
        (OMEGA*rho)
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal_energy invalid"
       stop
      endif

      return
      end subroutine INTERNAL_ENTROPY_jwl

 
! e=(E/rho) - (1/2) (u^2 + v^2)
      subroutine EOS_NAjwl(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: pressure
      real(amrex_real) A,B,R1,R2,GAMMA,RHOI,OMEGA

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      OMEGA=GAMMA-one

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid EOS_NAjwl: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid EOS_NAjwl: ",internal_energy
       stop
      endif
      pressure= &
        A*(one-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(one-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)+ &
        OMEGA*rho*internal_energy

      if (pressure.gt.zero) then
       !do nothing
      else
       print *,"vacuum error in NA JWL: ",pressure
       stop
      endif

      return
      end subroutine EOS_NAjwl

! initial sound speed is:
! C=7.8039D+10-5.484D+12 e^(-4.94)-0.09375D+12 e^(-1.21)=
      subroutine SOUNDSQR_NAjwl(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: soundsqr
      real(amrex_real) A,B,R1,R2,GAMMA,RHOI,OMEGA
      real(amrex_real) pressure,dp_de,dp_drho

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      OMEGA=GAMMA-one

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid SOUNDSQR_NAjwl: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid SOUNDSQR_NAjwl: ",internal_energy
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
 
      if (soundsqr.gt.zero) then
       !do nothing
      else
       print *,"soundsqr invalid in SOUNDSQR_NAjwl:",soundsqr
       stop
      endif

      return
      end subroutine SOUNDSQR_NAjwl


! e=(E/rho) - (1/2) (u^2 + v^2)
      subroutine EOS_Sandusky_jwl(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: pressure
      real(amrex_real) A,B,R1,R2,GAMMA,RHOI,OMEGA
      real(amrex_real) :: E0

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
       ! 1 J = 10^7 ergs
       ! 1 J/m^3=10 erg/cm^3
      E0=10.1D+10 !ergs/cm^3
      OMEGA=GAMMA-one

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid EOS_Sandusky_jwl: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid EOS_Sandusky_jwl: ",internal_energy
       stop
      endif
      pressure= &
        A*(one-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(one-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)+ &
        OMEGA*E0*rho/rhoI

      if (pressure.gt.zero) then
       !do nothing
      else
       print *,"vacuum error in Sandusky JWL: ",pressure
       stop
      endif

      return
      end subroutine EOS_Sandusky_jwl

! initial sound speed is:
! C=7.8039D+10-5.484D+12 e^(-4.94)-0.09375D+12 e^(-1.21)=
      subroutine SOUNDSQR_Sandusky_jwl(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: soundsqr
      real(amrex_real) A,B,R1,R2,GAMMA,RHOI,OMEGA
      real(amrex_real) pressure,dp_de,dp_drho
      real(amrex_real) :: E0

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
       ! 1 J = 10^7 ergs
       ! 1 J/m^3=10 erg/cm^3
      E0=10.1D+10 !ergs/cm^3
      OMEGA=GAMMA-one

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid SOUNDSQR_Sandusky_jwl: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid SOUNDSQR_Sandusky_jwl: ",internal_energy
       stop
      endif

      call EOS_Sandusky_jwl(rho,internal_energy,pressure)
      dp_de=zero
      dp_drho= &
        A*(one-OMEGA*rho/(R1*RHOI))*exp(-R1*RHOI/rho)* &
        R1*RHOI/(rho**2)- &
        (A*OMEGA/(R1*RHOI))*exp(-R1*RHOI/rho)+ &
        B*(one-OMEGA*rho/(R2*RHOI))*exp(-R2*RHOI/rho)* &
        R2*RHOI/(rho**2)- &
        B*(OMEGA/(R2*RHOI))*exp(-R2*RHOI/rho)+ &
        OMEGA*E0/rhoI
    
      soundsqr=(pressure*dp_de)/(rho**2)+dp_drho
 
      if (soundsqr.gt.zero) then
       !do nothing
      else
       print *,"soundsqr invalid in SOUNDSQR_Sandusky_jwl:",soundsqr
       stop
      endif

      return
      end subroutine SOUNDSQR_Sandusky_jwl

! e=(E/rho) - (1/2) (u^2 + v^2)
      subroutine EOS_jwlADIABAT(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: pressure
      real(amrex_real) A,B,R1,R2,GAMMA,RI,PI,RHOI,C,OMEGA

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      if (RHOI.gt.zero) then
       !do nothing
      else
       print *,"RHOI invalid"
       stop
      endif

      OMEGA=GAMMA-one
      RI=16.0D0        ! cm
      PI=7.8039D+10  ! dyne/cm^2
      C=PI-OMEGA*(A*exp(-R1)/R1+B*exp(-R2)/R2)
      if (C.gt.zero) then
       !do nothing
      else
       print *,"c invalid"
       stop
      endif
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid"
       stop
      endif
      pressure=(OMEGA*rho/RHOI)*( &
       A*exp(-R1*RHOI/rho)/R1+ &
       B*exp(-R2*RHOI/rho)/R2   )+ &
       C*( (RHOI/rho)**(-GAMMA) )
      if (pressure.gt.zero) then
       !do nothing
      else
       print *,"vacuum error"
       stop
      endif

      return
      end subroutine EOS_jwlADIABAT

      subroutine SOUNDSQR_jwlADIABAT(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: soundsqr
      real(amrex_real) A,B,R1,R2,GAMMA,RI,PI,RHOI,C,OMEGA

      call get_jwl_constants(A,B,GAMMA,R1,R2,RHOI)
      if (RHOI.gt.zero) then 
       !do nothing
      else
       print *,"RHOI invalid"
       stop
      endif

      OMEGA=GAMMA-one
      RI=16.0D0        ! cm
      PI=7.8039D+10  ! dyne/cm^2
      C=PI-OMEGA*(A*exp(-R1)/R1+B*exp(-R2)/R2)
      if (C.gt.zero) then
       !do nothing
      else
       print *,"c invalid"
       stop
      endif
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
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

      if (soundsqr.gt.zero) then
       !do nothing
      else
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

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) temperature
      real(amrex_real) R_pr,water_molar_mass,rhoMKS,rho_molar
      real(amrex_real) Vm,water_critical_temperature,water_critical_pressure
      real(amrex_real) water_critical_molar_volume,a,b,Tr,kappa,alpha
      real(amrex_real) water_acentric_factor
      real(amrex_real) water_c1,water_c2,water_c3,water_m,water_n
      real(amrex_real) Monnery_c4
      real(amrex_real) Vm_shift,adjusted_Vm

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

      real(amrex_real) rho,internal_energy,soundsqr,pressure
      real(amrex_real) drho,de,rho_plus,rho_minus,p_plus,p_minus
      real(amrex_real) e_plus,e_minus,dp_drho,dp_de

      if ((rho.gt.zero).and. &
          (internal_energy.gt.zero)) then

       call EOS_peng_robinson(rho,internal_energy,pressure)
       drho=EPS6*rho
       de=EPS6*internal_energy
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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,internal_energy,pressure,strain
      real(amrex_real) rho0,eta,E_ratio,denom,a_term,b_term,beta_term,alpha_term
      real(amrex_real) pressure1,pressure3
      real(amrex_real) e0,E_tillotson,T0_tillotson
      real(amrex_real) E_offset,P_offset
      real(amrex_real) E0_ratio,denom0,a0_term,b0_term

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

      real(amrex_real) rho,internal_energy,soundsqr,pressure
      real(amrex_real) rho0,drho,de,rho_plus,rho_minus,p_plus,p_minus
      real(amrex_real) e_plus,e_minus,dp_drho,dp_de

      rho0=fort_denconst(1)

      if ((rho.gt.zero).and. &
          (internal_energy.gt.zero).and. &
          (rho0.gt.zero).and. &
          (rho0.gt.rho_IV_tillotson).and. &
          (E0_tillotson.gt.zero).and. &
          (E_CV_tillotson.gt.E_IV_tillotson)) then

       call EOS_tillotson(rho,internal_energy,pressure)
       drho=EPS6*rho
       de=EPS6*internal_energy
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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,temperature,internal_energy,cv

      if ((rho.gt.zero).and.(internal_energy.gt.zero)) then
       cv=4.1855D+7
       temperature=internal_energy/cv
      else
       print *,"rho or internal_energy invalid"
       stop
      endif

      return
      end subroutine TEMPERATURE_tillotson


      subroutine EOS_wardlaw_tillotson(rho,internal_energy,pressure,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: pressure
      real(amrex_real) :: rho0,T0,e0,mu

      rho0=fort_denconst(im)
      T0=fort_tempconst(im)

      call INTERNAL_wardlaw_tillotson(rho0,T0,e0,im)

      if ((rho0.gt.zero).and. &
          (rho.gt.zero).and. &
          (T0.ge.zero).and. &
          (internal_energy.ge.zero).and. &
          (e0.ge.zero).and. &
          (im.ge.1).and. &
          (im.le.num_materials)) then
       !do nothing
      else
       print *,"EOS_wardlaw_tillotson invalid parms"
       print *,"rho0=",rho0
       print *,"rho=",rho
       print *,"T0=",T0
       print *,"internal_energy=",internal_energy
       print *,"e0=",e0
       print *,"im=",im
       stop
      endif

       !see PROBCOMMON.F90 for "wardlaw_tillotson" parameters.
      mu=rho/rho0-one
      pressure=P0_tillotson+omega_wardlaw_tillotson*rho* &
              (internal_energy-e0)+ &
              A_wardlaw_tillotson*mu+ &
              B_wardlaw_tillotson*(mu**2)+ &
              C_wardlaw_tillotson*(mu**3)
      if (pressure.ge.P_cav_tillotson) then
       !do nothing
      else if (pressure.lt.P_cav_tillotson) then
       pressure=P_cav_tillotson
      else
       print *,"pressure or pcav invalid (eos_wardlaw_tillotson): ", &
          pressure,P_cav_tillotson
       stop
      endif

      return
      end subroutine EOS_wardlaw_tillotson

! c^2 = dp/drho + p dp/de / rho^2
      subroutine SOUNDSQR_wardlaw_tillotson(rho,internal_energy,soundsqr,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: soundsqr
      real(amrex_real) :: rho0,T0,e0,mu,pressure,local_rho
      real(amrex_real) :: dmu,dpdrho,dpde,T0_sanity,e0_sanity
      real(amrex_real) :: pressure_sanity,sound_sanity

      rho0=fort_denconst(im)
      T0=fort_tempconst(im)

      call INTERNAL_wardlaw_tillotson(rho0,T0,e0,im)

      T0_sanity=293.0d0
      call INTERNAL_wardlaw_tillotson(rho0,T0_sanity,e0_sanity,im)

      if ((rho0.gt.zero).and. &
          (rho.gt.zero).and. &
          (T0.ge.zero).and. &
          (T0_sanity.gt.zero).and. &
          (internal_energy.ge.zero).and. &
          (e0.ge.zero).and. &
          (e0_sanity.gt.zero).and. &
          (im.ge.1).and. &
          (im.le.num_materials)) then
       !do nothing
      else
       print *,"soundsqr_wardlaw_tillotson invalid parms"
       print *,"rho0=",rho0
       print *,"rho=",rho
       print *,"T0=",T0
       print *,"T0_sanity=",T0_sanity
       print *,"internal_energy=",internal_energy
       print *,"e0=",e0
       print *,"e0_sanity=",e0_sanity
       print *,"im=",im
       stop
      endif

      call EOS_wardlaw_tillotson(rho,internal_energy,pressure,im)

      call EOS_wardlaw_tillotson(rho_cav_wardlaw_tillotson, &
        e0_sanity,pressure_sanity,im)

      if (abs(pressure_sanity-P_cav_tillotson).le. &
          0.01d0*P_cav_tillotson) then
       !do nothing
      else
       print *,"P_cav_tillotson test failed"
       print *,"pressure_sanity: ",pressure_sanity
       print *,"P_cav_tillotson: ",P_cav_tillotson
       stop
      endif

      if (pressure.le.P_cav_tillotson) then
       local_rho=rho_cav_wardlaw_tillotson
      else if (pressure.ge.P_cav_tillotson) then
       local_rho=rho
      else
       print *,"pressure out of range: ",pressure
       stop
      endif

      if (local_rho.lt.rho_cav_wardlaw_tillotson) then
       local_rho=rho_cav_wardlaw_tillotson
      else if (local_rho.ge.rho_cav_wardlaw_tillotson) then
       !do nothing
      else
       print *,"local_rho invalid: ",local_rho
       stop
      endif

      mu=local_rho/rho0-one
      dmu=one/rho0
      dpdrho=omega_wardlaw_tillotson*(internal_energy-e0)+ &
              A_wardlaw_tillotson*dmu+ &
              two*B_wardlaw_tillotson*mu*dmu+ &
              three*C_wardlaw_tillotson*mu*mu*dmu
      dpde=omega_wardlaw_tillotson*local_rho
      soundsqr=dpdrho+pressure*dpde/(local_rho**2)

      if (soundsqr.gt.zero) then
       !do nothing
      else
       print *,"soundsqr invalid: ",soundsqr
       stop
      endif

      mu=rho_cav_wardlaw_tillotson/rho0-one
      dpdrho=omega_wardlaw_tillotson*(e0_sanity-e0)+ &
              A_wardlaw_tillotson*dmu+ &
              two*B_wardlaw_tillotson*mu*dmu+ &
              three*C_wardlaw_tillotson*mu*mu*dmu
      dpde=omega_wardlaw_tillotson*rho_cav_wardlaw_tillotson
       ! c^2 = p_rho + p p_e/rho^2
      sound_sanity=dpdrho+pressure_sanity*dpde/(rho_cav_wardlaw_tillotson**2)
      sound_sanity=sqrt(sound_sanity)
      if (abs(sound_sanity-sound_cav_wardlaw_tillotson).le. &
          0.01d0*sound_cav_wardlaw_tillotson) then
       ! do nothing
      else
       print *,"sound_sanity invalid"
       print *,"dpdrho=",dpdrho
       print *,"dpde=",dpde
       print *,"mu=",mu
       print *,"dmu=",dmu
       print *,"omega_wardlaw_tillotson ",omega_wardlaw_tillotson
       print *,"A_wardlaw_tillotson ",A_wardlaw_tillotson
       print *,"B_wardlaw_tillotson ",B_wardlaw_tillotson
       print *,"C_wardlaw_tillotson ",C_wardlaw_tillotson
       print *,"rho0=",rho0
       print *,"rho_cav_wardlaw_tillotson=", &
         rho_cav_wardlaw_tillotson
       print *,"T0_sanity ",T0_sanity
       print *,"e0_sanity ",e0_sanity
       print *,"e0 ",e0
       print *,"sound_sanity=",sound_sanity
       print *,"sound_cav_wardlaw_tillotson=",sound_cav_wardlaw_tillotson
       stop
      endif 

      return
      end subroutine SOUNDSQR_wardlaw_tillotson


      subroutine INTERNAL_wardlaw_tillotson(rho,temperature, &
                      internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,temperature
      real(amrex_real), intent(out) :: internal_energy
      real(amrex_real) :: cv

      if ((rho.gt.zero).and.(temperature.ge.zero)) then
       cv=4.1855D+7
       internal_energy=cv*temperature
      else
       print *,"rho or temperature invalid:",rho,temperature
       stop
      endif

      return
      end subroutine INTERNAL_wardlaw_tillotson


      subroutine TEMPERATURE_wardlaw_tillotson(rho,temperature, &
                      internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: temperature
      real(amrex_real) :: cv

      if ((rho.gt.zero).and.(internal_energy.ge.zero)) then
       cv=4.1855D+7
       temperature=internal_energy/cv
      else
       print *,"rho or internal_energy invalid: ",rho,internal_energy
       stop
      endif

      return
      end subroutine TEMPERATURE_wardlaw_tillotson

       !material type=36
      subroutine EOS_Mie_Gruneison(rho,internal_energy,pressure,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: pressure
      real(amrex_real) :: rho0,T0,e0,mu,p0_factor
      real(amrex_real) :: V0,V,Gamma0,c0,s

      rho0=fort_denconst(im)
      T0=fort_tempconst(im)

      call INTERNAL_Mie_Gruneison(rho0,T0,e0,im)

      if ((rho0.gt.zero).and. &
          (rho.gt.zero).and. &
          (T0.ge.zero).and. &
          (internal_energy.ge.zero).and. &
          (e0.ge.zero).and. &
          (im.ge.1).and. &
          (im.le.num_materials)) then
       !do nothing
      else
       print *,"EOS_Mie_Gruneison invalid parms"
       print *,"rho0=",rho0
       print *,"rho=",rho
       print *,"T0=",T0
       print *,"internal_energy=",internal_energy
       print *,"e0=",e0
       print *,"im=",im
       stop
      endif

      V0=one/rho0
      V=one/rho
      Gamma0=fort_stiffGAMMA(im)
      c0=fort_stiff_sound_speed(im)
      if ((Gamma0.ge.one).and.(Gamma0.le.1.5d0).and.(c0.gt.zero)) then
       !do nothing
      else
       print *,"Gamma0 or c0 invalid:",Gamma0,c0
       stop
      endif
      s=1.49d0
      mu=one-rho0/rho
       !V0-V=1/rho0-1/rho=mu/rho0
       !V0-s(V0-V)=(1/rho0)(1-s mu)
       !p=c^2 rho0(mu/(1-s mu)^2)+
       !  gamma0 rho0 *
       !  (e-(1/2)c^2(mu/(1-s mu))^2)=
       !  rho0 (c^2 mu(1-gamma0 mu/2)/(1-s mu)^2)+gamma0 rho0 e

      p0_factor=1.0D+6/(Gamma0*rho0*e0)

      if (mu.ge.zero) then !rho>=rho0

       if ((mu.ge.zero).and.(mu.lt.one)) then
        !do nothing
       else
        print *,"mu invalid: ",mu
        stop
       endif

       pressure=rho0*(c0**2)*mu*(one-half*Gamma0*mu)/((one-s*mu)**2)+ &
        Gamma0*rho0*internal_energy*p0_factor

       if (pressure.gt.P_cav_mie_gruneisen) then
        !do nothing
       else if (pressure.le.P_cav_mie_gruneisen) then

        pressure=P_cav_mie_gruneisen

        if (pressure.gt.zero) then
         !do nothing
        else if (pressure.le.zero) then

         print *,"pressure invalid(a): ",pressure
         print *,"c0=",c0
         print *,"e0=",e0
         print *,"rho=",rho
         print *,"rho0=",rho0
         print *,"Gamma0=",Gamma0
         print *,"internal_energy=",internal_energy
         print *,"p0_factor=",p0_factor
         stop
        else
         print *,"pressure corrupt: ",pressure
         stop
        endif

       else
        print *,"pressure corrupt: ",pressure
        stop
       endif

      else if (mu.lt.zero) then !rho<rho0
       pressure=(c0**2)*(rho-rho0)+Gamma0*rho0*internal_energy*p0_factor
       if (pressure.gt.P_cav_mie_gruneisen) then
        !do nothing
       else if (pressure.le.P_cav_mie_gruneisen) then

        pressure=P_cav_mie_gruneisen

        if (pressure.gt.zero) then
         !do nothing
        else if (pressure.le.zero) then

         print *,"pressure invalid(b): ",pressure
         print *,"c0=",c0
         print *,"e0=",e0
         print *,"rho=",rho
         print *,"rho0=",rho0
         print *,"Gamma0=",Gamma0
         print *,"internal_energy=",internal_energy
         print *,"p0_factor=",p0_factor
         stop
        else
         print *,"pressure corrupt: ",pressure
         stop
        endif
       else
        print *,"pressure corrupt: ",pressure
        stop
       endif
      else
       print *,"mu invalid: ",mu
       stop
      endif 
              
      return
      end subroutine EOS_Mie_Gruneison

! c^2 = dp/drho + p dp/de / rho^2
      subroutine SOUNDSQR_Mie_Gruneison(rho,internal_energy,soundsqr,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: soundsqr
      real(amrex_real) :: rho0,T0,e0,mu,pressure,p0_factor
      real(amrex_real) :: dmu,dpdrho,dpde
      real(amrex_real) :: w,dtop_bottom,dbottom_top,dpdmu
      real(amrex_real) :: V0,V,Gamma0,c0,s

      rho0=fort_denconst(im)
      T0=fort_tempconst(im)

      call INTERNAL_Mie_Gruneison(rho0,T0,e0,im)

      if ((rho0.gt.zero).and. &
          (rho.gt.zero).and. &
          (T0.ge.zero).and. &
          (internal_energy.ge.zero).and. &
          (e0.ge.zero).and. &
          (im.ge.1).and. &
          (im.le.num_materials)) then
       !do nothing
      else
       print *,"soundsqr_Mie_Gruneison invalid parms"
       print *,"rho0=",rho0
       print *,"rho=",rho
       print *,"T0=",T0
       print *,"internal_energy=",internal_energy
       print *,"e0=",e0
       print *,"im=",im
       stop
      endif

      V0=one/rho0
      V=one/rho
      Gamma0=fort_stiffGAMMA(im)
      c0=fort_stiff_sound_speed(im)

      if ((Gamma0.ge.one).and.(Gamma0.le.1.5d0).and.(c0.gt.zero)) then
       !do nothing
      else
       print *,"Gamma0 or c0 invalid:",Gamma0,c0
       stop
      endif
      s=1.49d0
      mu=one-rho0/rho

      call EOS_Mie_Gruneison(rho,internal_energy,pressure,im)
      p0_factor=1.0D+6/(Gamma0*rho0*e0)

      dmu=rho0/(rho**2)
      w=one-s*mu
      dtop_bottom=rho0*(c0**2)*(w**2)*(one-Gamma0*mu)
      dbottom_top=rho0*(c0**2)*mu*(one-Gamma0*half*mu)*two*(-s)*w
      dpdmu=(dtop_bottom-dbottom_top)/(w**4)
      dpdrho=dpdmu*dmu
      dpde=Gamma0*rho0*p0_factor
      if (mu.ge.zero) then
       !do nothing
      else if (mu.lt.zero) then
       dpdrho=c0**2
      else
       print *,"mu invalid: ",mu
       stop
      endif 

      soundsqr=dpdrho+pressure*dpde/(rho**2)

      if (soundsqr.lt.c0**2) then
       soundsqr=c0**2
      else if (soundsqr.ge.c0**2) then
       !do nothing
      else
       print *,"soundsqr invalid: ",soundsqr
       stop
      endif

      if (soundsqr.gt.zero) then
       !do nothing
      else
       print *,"soundsqr invalid: ",soundsqr
       stop
      endif

      return
      end subroutine SOUNDSQR_Mie_Gruneison


      subroutine INTERNAL_Mie_Gruneison(rho,temperature, &
                      internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,temperature
      real(amrex_real), intent(out) :: internal_energy
      real(amrex_real) :: cv

      if ((rho.gt.zero).and.(temperature.ge.zero)) then
!      cv=4.1855D+7
       cv=fort_stiffCV(im)
       internal_energy=cv*temperature
      else
       print *,"rho or temperature invalid:",rho,temperature
       stop
      endif

      return
      end subroutine INTERNAL_Mie_Gruneison


      subroutine TEMPERATURE_Mie_Gruneison(rho,temperature, &
                      internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: temperature
      real(amrex_real) :: cv

      if ((rho.gt.zero).and.(internal_energy.ge.zero)) then
!      cv=4.1855D+7
       cv=fort_stiffCV(im)
       temperature=internal_energy/cv
      else
       print *,"rho or internal_energy invalid: ",rho,internal_energy
       stop
      endif

      return
      end subroutine TEMPERATURE_Mie_Gruneison



! A fully compressible, two-dimensional model of small, high-speed, cavitating
! nozzles, Schmidt et al, Atomization and Sprays, vol 9, 255-276 (1999)
! units must be CGS
      subroutine EOS_cav_nozzle(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) rho_g,rho_l,soundsqr_g,soundsqr_l,alpha,pgl,numer,denom

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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) rho_g,rho_l,soundsqr_g,soundsqr_l,alpha,rhomix,inv_soundsqr_mix

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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,temperature,internal_energy,cv

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

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) A,B,rhobar,GAMMA,pcav

      A=A_TAIT
      B=B_TAIT
      rhobar=RHOBAR_TAIT ! see PROBCOMMON.F90 for definition
      GAMMA=GAMMA_TAIT
      pcav=PCAV_TAIT

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid EOS_tait: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid EOS_tait ",internal_energy
       stop
      endif
      if (A.gt.zero) then
       !do nothing
      else
       print *,"A invalid EOS_tait ",A
       stop
      endif
      if (B.gt.zero) then
       !do nothing
      else
       print *,"B invalid EOS_tait ",B
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

      real(amrex_real) rho,pressure
      real(amrex_real) A,B,rhobar,GAMMA,pcav

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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,pcav,rhocav,pressure
      real(amrex_real) rho_sound

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid SOUNDSQR_tait ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid SOUNDSQR_tait ",internal_energy
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=RHOBAR_TAIT ! g/cm^3
      pcav=PCAV_TAIT
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (rhocav.gt.zero) then
       !do nothing
      else
       print *,"rhocav invalid SOUNDSQR_TAIT ",rhocav
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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,pcav,rhocav,pressure
      real(amrex_real) rho_sound

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

      real(amrex_real) rho,rho_in,T,T_in,DeDT
      real(amrex_real) cv,Tc
      real(amrex_real) A,B,rho0,rmax,rmin,P0,pressure
      real(amrex_real) cd,cc,cb,ca,T2,T3

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

      cd=2.273845+7.701613e-20*pressure*(pressure-1e6)
      cc=-2.279889e-3-3.654273e-13*(pressure-1e6)
      cb=6.106366e-06
      ca=-3.266302e-09

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

      real(amrex_real) rho_in,T_in
      real(amrex_real) rho,T,internal_energy,cv,Tc
      real(amrex_real) A,B,rho0,rmax,rmin,P0,pressure
      real(amrex_real) ce,cd,cc,cb,ca,T2,T3,T4

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

      ce=19.94245
      cd=2.273845+7.701613e-20*pressure*(pressure-1e6)
      cc=-2.279889e-3-3.654273e-13*(pressure-1e6)
      cb=6.106366e-06
      ca=-3.266302e-09

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

      real(amrex_real) rho_in
      real(amrex_real) rho,T,internal_energy,cv
      real(amrex_real) T0,T1,Tc,ie1,ie0,dT,dedT,rmin,rmax
      real(amrex_real), PARAMETER :: stub_ten=ten
      integer iter

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
       dT = sign(min(abs(dT),stub_ten),dT)
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

      real(amrex_real) rho_in
      real(amrex_real) rho,internal_energy,T,pressure
      real(amrex_real) A,B,rho0,rmax,rmin,P0,pcav,cv,Tc,Tred

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
      pressure=max(pressure,1.0D-5*P0)
      pressure=min(pressure,1.8D+9)  ! VIP limiter

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
      
      real(amrex_real) rho,pressure
      real(amrex_real) A,B,rhobar,pcav

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

      real(amrex_real) rho_in
      real(amrex_real) rmin,rmax
      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) c,c0,D,E,pressure,T,T0p5,T1p5,T3p0
      real(amrex_real) A,B,rhobar,pcav,rhocav,rho_sound

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
     !D = 0.005208-5.1495e-4*T-5.55e-14*T*pressure ! pressure needed in MPa
      D = 0.1652083+2.5e-3*T-5.85e-13*T*pressure ! modified correlation
      E = (-56.91+7.3674e-5*T*T+0.02260*T+463.5*exp(-0.001687*T))*1e7 !  convert from dyne/cm^2 to MPa
      c=c0/(1-D*log((E+pressure)/(E+1D6))) ! already in cm/s
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

      real(amrex_real) rho,temperature,internal_energy,cv


      if (rho.gt.zero) then 
       !do nothing
      else
       print *,"density negative"
       stop
      endif
      if (temperature.gt.zero) then
       !do nothing
      else
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

      real(amrex_real) rho,temperature,internal_energy,cv


      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative"
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
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

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) A,B,rhobar,pcav


      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.ge.0.001) then
       !do nothing
      else
       print *,"rhobar invalid in eos tait rho: ",rhobar
       stop
      endif

      pcav=PCAV_TAIT

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid (internal_energy) ",internal_energy
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

      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(out) :: pressure
      real(amrex_real) A,B,rhobar,pcav


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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,pcav,rhocav,pressure
      real(amrex_real) rho_sound


      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid"
       stop
      endif

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.ge.0.001) then
       !do nothing
      else
       print *,"rhobar invalid in soundsqr tait rho"
       stop
      endif

      pcav=PCAV_TAIT
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_TAIT) )
      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A
      
      if (rhocav.gt.zero) then
       !do nothing
      else
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


      subroutine INTERNAL_galinstan_rho(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,temperature,internal_energy,cv


      if (rho.gt.zero) then 
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (temperature.gt.zero) then
       !do nothing
      else
       print *,"temperature <=0: ",temperature
       stop
      endif
      cv=4.1855D+7
      internal_energy=temperature*cv

      return
      end subroutine INTERNAL_galinstan_rho

      subroutine TEMPERATURE_galinstan_rho(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,temperature,internal_energy,cv


      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy <=0: ",internal_energy
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_galinstan_rho

      subroutine EOS_galinstan_rho(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) A,B,rhobar,pcav


      A=A_GALINSTAN   ! dyne/cm^2
      B=B_GALINSTAN  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.ge.0.001) then
       !do nothing
      else
       print *,"rhobar invalid in eos galinstan rho: ",rhobar
       stop
      endif

      pcav=PCAV_TAIT

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid: ",internal_energy
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_GALINSTAN - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif

      return
      end subroutine EOS_galinstan_rho


      subroutine SOUNDSQR_galinstan_rho(rho,internal_energy,soundsqr)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,pcav,rhocav,pressure
      real(amrex_real) rho_sound


      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid: ", internal_energy
       stop
      endif

      A=A_GALINSTAN   ! dyne/cm^2
      B=B_GALINSTAN  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.ge.0.001) then
       !do nothing
      else
       print *,"rhobar invalid in soundsqr galinstan rho: ",rhobar
       stop
      endif

      pcav=PCAV_TAIT
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_GALINSTAN) )
      pressure=B*( (rho/rhobar)**GAMMA_GALINSTAN - one ) + A
      
      if (rhocav.gt.zero) then
       !do nothing
      else
       print *,"rhocav invalid: ",rhocav
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_GALINSTAN*B/rhobar)* &
        ( (rho_sound/rhobar)**(GAMMA_GALINSTAN-one) )

      return
      end subroutine SOUNDSQR_galinstan_rho


      subroutine INTERNAL_elastic_rho(rho,temperature,internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,temperature
      real(amrex_real), intent(out) :: internal_energy
      real(amrex_real) :: cv

      if ((im.ge.1).and.(im.le.num_materials)) then
       !do nothing
      else
       print *,"im invalid: ",im
       stop
      endif

      if (rho.gt.zero) then 
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (temperature.gt.zero) then
       !do nothing
      else
       print *,"temperature <=0: ",temperature
       stop
      endif
      cv=4.1855D+7
      internal_energy=temperature*cv

      return
      end subroutine INTERNAL_elastic_rho

      subroutine TEMPERATURE_elastic_rho(rho,temperature,internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: temperature
      real(amrex_real) :: cv


      if ((im.ge.1).and.(im.le.num_materials)) then
       !do nothing
      else
       print *,"im invalid: ",im
       stop
      endif

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy <=0: ",internal_energy
       stop
      endif
      cv=4.1855D+7
      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_elastic_rho

      subroutine EOS_elastic_rho(rho,internal_energy,pressure,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real),intent(in):: rho,internal_energy
      real(amrex_real),intent(out):: pressure
      real(amrex_real) A,B,rhobar,pcav

      if ((im.ge.1).and.(im.le.num_materials)) then
       !do nothing
      else
       print *,"im invalid: ",im
       stop
      endif

      A=A_ELASTIC   ! dyne/cm^2
      B=B_ELASTIC  ! dyne/cm^2
      rhobar=fort_denconst(im) ! g/cm^3

      if (rhobar.ge.0.001) then
       !do nothing
      else
       print *,"rhobar invalid in eos elastic rho: ",rhobar
       stop
      endif

      pcav=PCAV_ELASTIC

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid: ",internal_energy
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_ELASTIC - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif

      return
      end subroutine EOS_elastic_rho


      subroutine SOUNDSQR_elastic_rho(rho,internal_energy,soundsqr,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real),intent(in) :: rho,internal_energy
      real(amrex_real),intent(out) :: soundsqr
      real(amrex_real) A,B,rhobar,pcav,rhocav,pressure
      real(amrex_real) rho_sound


      if ((im.ge.1).and.(im.le.num_materials)) then
       !do nothing
      else
       print *,"im invalid: ",im
       stop
      endif
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid: ",internal_energy
       stop
      endif

      A=A_ELASTIC ! dyne/cm^2
      B=B_ELASTIC ! dyne/cm^2
      rhobar=fort_denconst(im) ! g/cm^3

      if (rhobar.ge.0.001) then
       !do nothing
      else
       print *,"rhobar invalid in soundsqr elastic rho: ",rhobar
       stop
      endif

      pcav=PCAV_ELASTIC
      rhocav=rhobar*( ((pcav-A)/B+one)**(one/GAMMA_ELASTIC) )
      pressure=B*( (rho/rhobar)**GAMMA_ELASTIC - one ) + A
      
      if (rhocav.gt.zero) then
       !do nothing
      else
       print *,"rhocav invalid: ",rhocav
       stop
      endif

      if (pressure.lt.pcav) then
       rho_sound=rhocav
      else
       rho_sound=rho
      endif
      soundsqr=(GAMMA_ELASTIC*B/rhobar)* &
        ( (rho_sound/rhobar)**(GAMMA_ELASTIC-one) )

      return
      end subroutine SOUNDSQR_ELASTIC_rho

      subroutine INTERNAL_tait_rhohydro(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) A,B,rhobar,pcav


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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,pcav,rhocav,pressure
      real(amrex_real) rho_sound


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) A,B,rhobar,pcav


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

      real(amrex_real) rho,pressure
      real(amrex_real) A,B,rhobar,pcav


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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,pcav,rhocav,pressure
      real(amrex_real) rho_sound


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) A,B,rhobar,pcav


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

      real(amrex_real) rho,pressure
      real(amrex_real) A,B,rhobar,pcav


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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,pcav,rhocav,pressure
      real(amrex_real) rho_sound


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) A,B,rhobar,pcav,GAMMA_KOREN


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

      real(amrex_real) rho,pressure
      real(amrex_real) A,B,rhobar,GAMMA_KOREN,pcav


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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,GAMMA_KOREN,pcav,rhocav,pressure
      real(amrex_real) rho_sound


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,temperature,internal_energy,cv


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

      real(amrex_real) rho,internal_energy,pressure
      real(amrex_real) A,B,rhobar,GAMMA_KOREN,pcav


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

      real(amrex_real) rho,pressure
      real(amrex_real) A,B,rhobar,GAMMA_KOREN,pcav


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

      real(amrex_real) rho,internal_energy,soundsqr
      real(amrex_real) A,B,rhobar,GAMMA_KOREN,pcav,rhocav,pressure
      real(amrex_real) rho_sound


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
      real(amrex_real) R,cp,cv,gamma_constant,omega

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
      real(amrex_real), INTENT(out) :: R,cp,cv,gamma_constant,omega

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

      real(amrex_real), INTENT(in) :: xpos(SDIM)
      real(amrex_real), INTENT(inout) :: rho
      real(amrex_real), INTENT(inout) :: pres
      integer, INTENT(in) :: from_boundary_hydrostatic
      real(amrex_real) denfree,zfree
      real(amrex_real) z_at_depth

      ! in tait_hydrostatic_pressure_density
      if ((probtype.eq.36).and.(axis_dir.eq.2)) then  ! spherical explosion
       rho=one
       call EOS_tait_ADIABATIC(rho,pres)
      else if ((probtype.eq.36).and.(axis_dir.eq.310)) then !hydrobulge
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

       if (probloy.eq.zero) then
        !do nothing
       else
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
       else if (xpos(SDIM).le.zfree) then
        if (z_at_depth.ne.zfree) then
         if (density_at_depth.gt.zero) then
          rho= &
           ((density_at_depth-denfree)/ &
            (z_at_depth-zfree))*(xpos(SDIM)-zfree)+denfree
         else 
          print *,"density_at_depth invalid"
          stop
         endif
        else 
         print *,"z_at_depth invalid"
         stop
        endif
       else
        print *,"xpos invalid"
        stop
       endif
       call EOS_tait_ADIABATIC(rho,pres)

      else if ((probtype.eq.46).and.(SDIM.eq.2)) then  ! cavitation

       if (probloy.eq.zero) then
        !do nothing
       else
        print *,"probloy must be 0 for cavitation problem:",probloy
        stop
       endif

       ! yblob is distance from domain bottom of charge/sphere
       ! zblob is depth of charge (for jwl problem)
       denfree=one

       if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
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

       else if (axis_dir.eq.10) then
        rho=denfree
        pres=zero
       else if (axis_dir.eq.11) then
        rho=denfree
        pres=zero
       else if (axis_dir.eq.20) then
        print *,"there is no gravity for the CODY ESTEBE created test problem"
        stop
       else
        print *,"axis_dir out of range: ",axis_dir
        stop
       endif

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

       !material_type=14
      subroutine EOS_air_rho2(rho,internal_energy,pressure,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real),intent(in) :: rho,internal_energy
      real(amrex_real),intent(out) :: pressure
      real(amrex_real) omega
      real(amrex_real) gamma_constant
      real(amrex_real) cp,cv,R,pressure_adjust,preshydro,rhohydro
      real(amrex_real) xpos(SDIM)
      integer, PARAMETER :: from_boundary_hydrostatic=0

      call air_parms(R,cp,cv,gamma_constant,omega) 
      if (rho.gt.zero) then
       ! do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be negative: ",internal_energy
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error: ",cv
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error: ",cp
       stop
      endif
      if (omega.gt.zero) then
       !do nothing
      else
       print *,"omega error: ",omega
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
      pressure_adjust=omega*fort_denconst(im)* &
         cv*fort_tempconst(im)

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


      subroutine EOS_air_rho2_ADIABAT(rho,pressure,im)
      use probcommon_module
      IMPLICIT NONE

      integer,intent(in) :: im
      real(amrex_real),intent(in) :: rho
      real(amrex_real),intent(out) :: pressure
      real(amrex_real) gamma_constant
      real(amrex_real) cp,cv,R,rhohydro,omega
      real(amrex_real) RHOI,PI
      real(amrex_real) xpos(SDIM)
      integer, PARAMETER :: from_boundary_hydrostatic=0

      RHOI=fort_denconst(im)
      call general_hydrostatic_pressure(PI)
    
      call air_parms(R,cp,cv,gamma_constant,omega)
      if (RHOI.gt.zero) then
       !do nothing
      else
       print *,"RHOI invalid"
       stop
      endif
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative"
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error"
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
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


      subroutine ENTROPY_air_rho2(rho,internal_energy,entropy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real),intent(in) :: rho,internal_energy
      real(amrex_real),intent(out) :: entropy
      real(amrex_real) pressure,press_adiabat

      call EOS_air_rho2(rho,internal_energy,pressure,im)
      call EOS_air_rho2_ADIABAT(rho,press_adiabat,im)
      if (press_adiabat.gt.zero) then
       !do nothing
      else
       print *,"press_adiabat invalid"
       stop
      endif
      entropy=pressure/press_adiabat

      return
      end subroutine ENTROPY_air_rho2

      subroutine INTERNAL_ENTROPY_air_rho2(rho,entropy,internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real),intent(in) :: rho,entropy
      real(amrex_real),intent(out) :: internal_energy
      real(amrex_real) press_adiabat,unit_internal_energy,unit_press

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid"
       stop
      endif
      if (entropy.gt.zero) then
       !do nothing
      else
       print *,"entropy invalid"
       stop
      endif
      call EOS_air_rho2_ADIABAT(rho,press_adiabat,im)
      unit_internal_energy=one
      call EOS_air_rho2(rho,unit_internal_energy,unit_press,im)
      internal_energy=press_adiabat*entropy/unit_press
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal_energy invalid"
       stop
      endif

      return
      end subroutine INTERNAL_ENTROPY_air_rho2


      subroutine SOUNDSQR_air_rho2(rho,internal_energy,soundsqr,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real),intent(in) :: rho,internal_energy
      real(amrex_real) gamma_constant,pressure,omega
      real(amrex_real), intent(out) :: soundsqr
      real(amrex_real) cp,cv,R,pressure_adjust,preshydro,rhohydro
      real(amrex_real) xpos(SDIM)
      integer, PARAMETER :: from_boundary_hydrostatic=0

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative"
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error"
       stop
      endif
      if (cp.gt.zero) then 
       !do nothing
      else
       print *,"cp error"
       stop
      endif
      pressure_adjust=omega*fort_denconst(im)* &
         cv*fort_tempconst(im)
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


      subroutine INTERNAL_air_rho2(rho,temperature,internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real),intent(in) :: rho,temperature
      real(amrex_real),intent(out) :: internal_energy
      real(amrex_real) gamma_constant
      real(amrex_real) cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (temperature.gt.zero) then
       !do nothing
      else
       print *,"temperature cannot be <=0 in internal_air_rho2"
       print *,temperature, rho, R, cp, cv
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error: ",cv
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error: ",cp
       stop
      endif

      internal_energy=cv*temperature

      return
      end subroutine INTERNAL_air_rho2


      subroutine TEMPERATURE_air_rho2(rho,temperature,internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer, intent(in) :: im
      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) ::  temperature
      real(amrex_real) gamma_constant
      real(amrex_real) cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be <=0 ",internal_energy
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error: ",cv
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error: ",cp
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

        !see ``A level set approach to Eulerian-Lagrangian Coupling''
        !Arienti, Hung, Morano, Shepherd
        !see Appendix C ``perfect gas analytical solutions''
        !Joeseph E. Shepherd
        !The dynamics and thermodynamics of Compressible Flow
        !Ascher H. Shapiro
        !material_type=5
      subroutine postshock_air(Mn,rho1,T1,rho2,T2,P2)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: Mn,rho1,T1
      real(amrex_real), intent(out) :: rho2,T2,P2
      real(amrex_real) :: cp,cv,gamma_constant,gap1,R,P1,Msq,omega

      ! check compatibility with gas EOS
      if (Mn.ge.1d0) then
       ! do nothing
      else
       print *,"ERROR - postshock_air: Mn",Mn
       stop
      endif

      call air_parms(R,cp,cv,gamma_constant,omega)

      gap1 = one+gamma_constant
      Msq=Mn*Mn
      P1=R*rho1*T1
      rho2=rho1/((gamma_constant-one)/gap1+two/gap1/Msq)
      P2=P1*(one+two*gamma_constant/gap1*(Msq-one))
      T2=P2/rho2/R

      return
      end subroutine postshock_air

      subroutine vnpostshock_air(T1,Mn,vn)
      IMPLICIT NONE

      real(amrex_real), intent(in) :: T1,Mn
      real(amrex_real), intent(out) :: vn
      real(amrex_real) :: cn
      REAL_T cv,gamma_constant,R,omega,cp

      call air_parms(R,cp,cv,gamma_constant,omega)

      cn=sqrt(gamma_constant*R*T1)
      vn=cn*(Mn-two/(gamma_constant+one)*(Mn-one/Mn))

      return
      end subroutine vnpostshock_air

        !material_type=5
      subroutine EOS_air(rho,internal_energy,pressure)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: pressure
      real(amrex_real) gamma_constant
      real(amrex_real) cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be <=0: ",internal_energy
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error: ",cv
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error: ",cp
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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) cp,cv,R,omega

      call simple_air_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative"
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error"
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error"
       stop
      endif
      pressure=omega*rho*internal_energy

      return
      end subroutine EOS_simple_air


      subroutine EOS_air_ADIABAT(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,gamma_constant,pressure
      real(amrex_real) cp,cv,R,RHOI,PI,omega

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
       
      real(amrex_real) rho
      real(amrex_real) term1,term2

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

      real(amrex_real) rho,gamma_constant,pressure
      real(amrex_real) cp,cv,R,RHOI,PI,omega

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

      real(amrex_real) rho,internal_energy,pressure,entropy,press_adiabat

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

      real(amrex_real) rho,entropy,internal_energy
      real(amrex_real) press_adiabat,unit_internal_energy,unit_press

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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) cp,cv,R

    
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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) soundsqr
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) soundsqr
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) soundsqr
      real(amrex_real) cp,cv,R

    
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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R

    
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

      real(amrex_real), intent(in) :: rho,internal_energy
      real(amrex_real), intent(out) :: temperature
      real(amrex_real) gamma_constant
      real(amrex_real) cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
      if (rho.gt.zero) then
       ! do nothing
      else
       print *,"rho invalid: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       ! do nothing
      else
       print *,"internal energy invalid: ",internal_energy
       stop
      endif
      if (cv.gt.zero) then
       ! do nothing
      else
       print *,"cv invalid: ",cv
       stop
      endif
      if (cp.gt.zero) then
       ! do nothing
      else
       print *,"cp invalid: ",cp
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine TEMPERATURE_air


      subroutine TEMPERATURE_simple_air(rho,temperature,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

    
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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R

    
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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) soundsqr
      real(amrex_real) cp,cv,R,phyd,rho0,omega


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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) soundsqr
      real(amrex_real) cp,cv,R,phyd,rho0,omega


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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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
      real(amrex_real) R,cp,cv,gamma_constant,omega

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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) cp,cv,R,omega

      call SF6_parms(R,cp,cv,gamma_constant,omega)
       
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative"
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error"
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error"
       stop
      endif
      pressure=omega*rho*internal_energy

      return
      end subroutine EOS_SF6


      subroutine EOS_SF6_ADIABAT(rho,pressure)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,gamma_constant,pressure
      real(amrex_real) cp,cv,R,RHOI,PI,omega

      RHOI=fort_denconst(2)
      call general_hydrostatic_pressure(PI)
      call SF6_parms(R,cp,cv,gamma_constant,omega)
       
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative"
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error"
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error"
       stop
      endif
      pressure=PI*(rho/RHOI)**gamma_constant

      return
      end subroutine EOS_SF6_ADIABAT

      subroutine ENTROPY_SF6(rho,internal_energy,entropy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,internal_energy,pressure,entropy,press_adiabat

      call EOS_SF6(rho,internal_energy,pressure)
      call EOS_SF6_ADIABAT(rho,press_adiabat)
      if (press_adiabat.gt.zero) then
       !do nothing
      else
       print *,"press_adiabat invalid"
       stop
      endif
      entropy=pressure/press_adiabat

      return
      end subroutine ENTROPY_SF6

      subroutine INTERNAL_ENTROPY_SF6(rho,entropy,internal_energy)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real) rho,entropy,internal_energy
      real(amrex_real) press_adiabat,unit_internal_energy,unit_press

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid"
       stop
      endif
      if (entropy.gt.zero) then
       !do nothing
      else
       print *,"entropy invalid"
       stop
      endif
      call EOS_SF6_ADIABAT(rho,press_adiabat)
      unit_internal_energy=one
      call EOS_SF6(rho,unit_internal_energy,unit_press)
      internal_energy=press_adiabat*entropy/unit_press
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal_energy invalid"
       stop
      endif

      return
      end subroutine INTERNAL_ENTROPY_SF6

      subroutine TEMPERATURE_ENTROPY_SF6(rho,entropy,temperature)
      use probcommon_module
      IMPLICIT NONE
      
      real(amrex_real) R,cp,cv,gamma_constant,omega 
      real(amrex_real) rho,entropy,temperature,RHOI

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

      real(amrex_real) R,cp,cv,gamma_constant,omega 
      real(amrex_real) rho,temperature,entropy,RHOI

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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) soundsqr
      real(amrex_real) cp,cv,R,omega

      call SF6_parms(R,cp,cv,gamma_constant,omega)
    
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative"
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be <=0"
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error"
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      integer im
      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) cp,cv,PP

      if (im.lt.1) then
       print *,"im invalid65"
       stop
      endif
      cp=get_user_stiffCP(im)
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp invalid: ",cp
       stop
      endif
      gamma_constant=fort_stiffGAMMA(im)
      if (gamma_constant.gt.zero) then
       !do nothing
      else
       print *,"gamma_constant invalid: ",gamma_constant
       stop
      endif
      cv=cp/gamma_constant

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be <=0: ",internal_energy
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error: ",cv
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error: ",cp
       stop
      endif
      PP=fort_stiffPINF(im)
      if (PP.ge.zero) then
       !do nothing
      else
       print *,"fort_stiff PINF invalid: ",PP
       stop
      endif
      pressure=(gamma_constant-one)*rho*internal_energy- &
       gamma_constant*PP
      if (pressure.lt.VOFTOL) then
       pressure=VOFTOL
      endif

      return
      end subroutine EOS_stiffened

      subroutine SOUNDSQR_stiffened(rho,internal_energy,soundsqr,im)
      use probcommon_module
      IMPLICIT NONE

      integer im
      real(amrex_real) rho,internal_energy,gamma_constant,pressure
      real(amrex_real) soundsqr
      real(amrex_real) cp,cv,PP

      if (im.lt.1) then
       print *,"im invalid66"
       stop
      endif
      cp=get_user_stiffCP(im)
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp invalid: ",cp
       stop
      endif
      gamma_constant=fort_stiffGAMMA(im)
      if (gamma_constant.gt.zero) then
       !do nothing
      else
       print *,"gamma_constant invalid: ",gamma_constant
       stop
      endif
      cv=cp/gamma_constant

      if (rho.gt.zero) then
       !do nothing
      else
       print *,"density negative: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"internal energy cannot be <=0: ",internal_energy
       stop
      endif
      if (cv.gt.zero) then
       !do nothing
      else
       print *,"cv error: ",cv
       stop
      endif
      if (cp.gt.zero) then
       !do nothing
      else
       print *,"cp error: ",cp
       stop
      endif
      PP=fort_stiffPINF(im)
      if (PP.ge.zero) then
       !do nothing
      else
       print *,"fort_stiff PINF invalid: ",PP
       stop
      endif

      pressure=(gamma_constant-one)*rho*internal_energy- &
       gamma_constant*PP
      if (pressure.lt.VOFTOL) then
       pressure=VOFTOL
      endif

      soundsqr=gamma_constant*(pressure+PP)/rho

      return
      end subroutine SOUNDSQR_stiffened


      subroutine INTERNAL_stiffened(rho,temperature,internal_energy,im)
      use probcommon_module
      IMPLICIT NONE

      integer im
      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,PP

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

      integer im
      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,PP

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

      real(amrex_real) rho,internal_energy,gamma_constant,pressure,omega
      real(amrex_real) soundsqr
      real(amrex_real) cp,cv,R


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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real) rho,temperature,internal_energy,gamma_constant
      real(amrex_real) cp,cv,R,omega

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

      real(amrex_real), INTENT(in) :: den
      real(amrex_real), INTENT(out) :: massfrac_parm(num_species_var+1)
      integer, INTENT(in) :: im
      integer :: ispec,var_comp

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

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho,internal_energy
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var+1)
      real(amrex_real), INTENT(inout) :: pressure
      real(amrex_real) :: T

      if ((im.ge.1).and.(im.le.num_materials)) then
       !do nothing
      else
       print *,"im invalid: ",im
       stop
      endif
      if (rho.gt.zero) then
       !do nothing
      else
       print *,"rho invalid: ",rho
       stop
      endif
      if (internal_energy.gt.zero) then
       !do nothing
      else
       print *,"e invalid: ",internal_energy
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
       call EOS_air_rho2(rho,internal_energy,pressure,im)
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
      else if (imattype.eq.25) then !uses denconst(num_materials)
       call EOS_elastic_rho(rho,internal_energy,pressure,im)
      else if (imattype.eq.34) then
       call EOS_galinstan_rho(rho,internal_energy,pressure)
      else if (imattype.eq.35) then
       call EOS_wardlaw_tillotson(rho,internal_energy,pressure,im)
      else if (imattype.eq.36) then
       call EOS_Mie_Gruneison(rho,internal_energy,pressure,im)
      else if (imattype.eq.37) then
       call EOS_Sandusky_jwl(rho,internal_energy,pressure)
      else
       print *,"imattype invalid EOS_material_CORE: ",imattype
       stop
      endif

      return
      end subroutine EOS_material_CORE

      subroutine dVdT_material_CORE(dVdT,massfrac_var, &
        pressure,temperature, &
        imattype,im)
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: pressure,temperature
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var+1)
      real(amrex_real), INTENT(out) :: dVdT


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
          (imattype.le.35)) then
       dVdT=(fort_stiffGAMMA(im)-one) * fort_stiffCV(im)/pressure
      else if (imattype.eq.0) then
       dVdT=zero
      else if (imattype.eq.999) then
       dVdT=zero
      else
       print *,"imattype invalid in dVdT_material_CORE: ",imattype
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

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho,internal_energy
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var+1)
      real(amrex_real), INTENT(out) :: soundsqr

      if ((im.ge.1).and.(im.le.num_materials)) then
       !do nothing
      else
       print *,"im invalid: ",im
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
       call SOUNDSQR_air_rho2(rho,internal_energy,soundsqr,im)
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
      else if (imattype.eq.25) then
       call SOUNDSQR_elastic_rho(rho,internal_energy,soundsqr,im)
      else if (imattype.eq.34) then
       call SOUNDSQR_galinstan_rho(rho,internal_energy,soundsqr)
      else if (imattype.eq.35) then
       call SOUNDSQR_wardlaw_tillotson(rho,internal_energy,soundsqr,im)
      else if (imattype.eq.36) then
       call SOUNDSQR_Mie_Gruneison(rho,internal_energy,soundsqr,im)
      else if (imattype.eq.37) then
       call SOUNDSQR_Sandusky_jwl(rho,internal_energy,soundsqr)
      else
       print *,"imattype invalid SOUNDSQR_material_CORE: ",imattype
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

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho,temperature
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var+1)
      real(amrex_real), INTENT(out) :: internal_energy
      real(amrex_real) local_internal_energy

      if ((im.ge.1).and.(im.le.num_materials)) then
       !do nothing
      else
       print *,"im invalid: ",im
       stop
      endif
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
       call INTERNAL_air_rho2(rho,temperature,local_internal_energy,im)
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
      else if (imattype.eq.25) then
       call INTERNAL_elastic_rho(rho,temperature,local_internal_energy,im)
      else if (imattype.eq.34) then
       call INTERNAL_galinstan_rho(rho,temperature,local_internal_energy)
      else if (imattype.eq.35) then
       call INTERNAL_wardlaw_tillotson(rho,temperature,local_internal_energy,im)
      else if (imattype.eq.36) then
       call INTERNAL_Mie_Gruneison(rho,temperature,local_internal_energy,im)
      else if (imattype.eq.37) then
       call INTERNAL_jwl(rho,temperature,local_internal_energy)
      else
       print *,"imattype invalid INTERNAL_material_CORE: ",imattype
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

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho,internal_energy
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var+1)
      real(amrex_real), INTENT(out) :: temperature

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
       call TEMPERATURE_air_rho2(rho,temperature,internal_energy,im)
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
      else if (imattype.eq.25) then
       call TEMPERATURE_elastic_rho(rho,temperature,internal_energy,im)
      else if (imattype.eq.34) then
       call TEMPERATURE_galinstan_rho(rho,temperature,internal_energy)
      else if (imattype.eq.35) then
       call TEMPERATURE_wardlaw_tillotson(rho,temperature,internal_energy,im)
      else if (imattype.eq.36) then
       call TEMPERATURE_Mie_Gruneison(rho,temperature,internal_energy,im)
      else if (imattype.eq.37) then
       call TEMPERATURE_jwl(rho,temperature,internal_energy)
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

      real(amrex_real) pres


      pres=1.0D+6

      return
      end subroutine general_hydrostatic_pressure


      subroutine rampvel(time,vel)
      use probcommon_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(out) :: vel
      real(amrex_real) :: tcutoff

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

      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(out) :: xx_vel,yy_vel,zz_vel

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
      ! if num_materials=4, 12 13 14 23 24 34

      subroutine get_CL_iten(im,im_opp,im_3,iten_13,iten_23, &
       user_tension,cos_angle,sin_angle)
      use probcommon_module
      IMPLICIT NONE

      integer,intent(in) :: im,im_opp,im_3
      integer,intent(out) :: iten_13,iten_23
      integer iten
      real(amrex_real), intent(in) :: user_tension(num_interfaces)
      real(amrex_real), intent(out) :: cos_angle,sin_angle

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

      if (num_materials.le.2) then
       print *,"num_materials too small for CL treatment"
       stop
       ! 12 13 23
      else if (num_materials.eq.3) then
       if ((im.eq.1).and.(im_opp.eq.2).and.(im_3.eq.3)) then
        iten_13=2
        iten_23=3
       ! 12 13 23
       else if ((im.eq.2).and.(im_opp.eq.3).and.(im_3.eq.1)) then
        iten_13=1
        iten_23=2
       else if ((im.eq.1).and.(im_opp.eq.3).and.(im_3.eq.2)) then
        iten_13=1
        iten_23=3
       else
        print *,"combination of im,im_opp,im_3 invalid num_materials=", &
                num_materials
        print *,"im=",im
        print *,"im_opp=",im_opp
        print *,"im_3=",im_3
        stop
       endif
       ! 12 13 14 23 24 34
      else if (num_materials.eq.4) then
       if ((im.eq.1).and.(im_opp.eq.2).and.(im_3.eq.3)) then
        iten_13=2
        iten_23=4
       ! 12 13 14 23 24 34
       else if ((im.eq.2).and.(im_opp.eq.3).and.(im_3.eq.1)) then
        iten_13=1
        iten_23=2
       else if ((im.eq.1).and.(im_opp.eq.3).and.(im_3.eq.2)) then
        iten_13=1
        iten_23=4
       ! 12 13 14 23 24 34
       else if ((im.eq.1).and.(im_opp.eq.2).and.(im_3.eq.4)) then
        iten_13=3
        iten_23=5
       else if ((im.eq.2).and.(im_opp.eq.4).and.(im_3.eq.1)) then
        iten_13=1
        iten_23=3
       ! 12 13 14 23 24 34
       else if ((im.eq.1).and.(im_opp.eq.4).and.(im_3.eq.2)) then
        iten_13=1
        iten_23=5
       else if ((im.eq.2).and.(im_opp.eq.3).and.(im_3.eq.4)) then
        iten_13=5
        iten_23=6
       ! 12 13 14 23 24 34
       else if ((im.eq.2).and.(im_opp.eq.4).and.(im_3.eq.3)) then
        iten_13=4
        iten_23=6
       else if ((im.eq.3).and.(im_opp.eq.4).and.(im_3.eq.2)) then
        iten_13=4
        iten_23=5
       ! 12 13 14 23 24 34
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
        print *,"combination of im,im_opp,im_3 invalid num_materials=", &
                num_materials
        print *,"im=",im
        print *,"im_opp=",im_opp
        print *,"im_3=",im_3
        stop
       endif
      else
       print *,"combination not supported"
       stop
      endif

      call get_iten(im,im_opp,iten)

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
subroutine drop_slope_dist(x,y,z,time, &
   maxtall,dist,dist_truncate)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x,y,z,time
real(amrex_real), INTENT(out) :: dist,dist_truncate
real(amrex_real), INTENT(in) :: maxtall
integer im,im_opp,im_3,iten_13,iten_23,imloop
integer iten
real(amrex_real) cos_angle,sin_angle
real(amrex_real) term1,Vtarget,radnew,vert,test_angle
real(amrex_real) xprime,yprime,zprime,rprime,rtop,rbot
real(amrex_real) xcheck,ycheck,zcheck
real(amrex_real) xvec(SDIM)
real(amrex_real) marangoni_temp(num_materials)
integer im_solid_substrate
real(amrex_real), allocatable, dimension(:) :: user_tension

if (probtype.eq.55) then

 im_solid_substrate=im_solid_primary()

 if (is_rigid(im_solid_substrate).eq.1) then
  !do nothing
 else
  print *,"expecting is_rigid(im_solid_substrate).eq.1"
  stop
 endif

 xvec(1)=x
 xvec(2)=y
 if (SDIM.eq.3) then
  xvec(SDIM)=z
 endif

 allocate(user_tension(num_interfaces))

 if (SDIM.eq.2) then
  if (abs(y-z).le.EPS2) then
   !do nothing
  else
   print *,"y=z in 2d expected: drop_slope_dist"
   print *,"x,y,z= ",x,y,z
   stop
  endif
 endif
 if (maxtall.gt.zero) then
  ! do nothing
 else
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
      (abs(xcheck).lt.EPS7).and. &
      (abs(ycheck).lt.EPS7).and. &
      (abs(zcheck).lt.EPS7)) then
   im=1 ! liquid
   im_opp=2 ! gas
   im_3=im_solid_substrate
   call get_iten(im,im_opp,iten)
   do imloop=1,num_materials
    marangoni_temp(imloop)=room_temperature ! 293.0d0 if double precision
   enddo
   call get_user_tension(xvec,time, &
     fort_tension_init,user_tension,marangoni_temp)
     ! find the angle between the "im,im_3" interface and the
     ! "im,im_opp" interface.
     ! i.e. between the liquid/substrate and liquid/gas interfaces.
   call get_CL_iten(im,im_opp,im_3,iten_13,iten_23, &
    user_tension,cos_angle,sin_angle)

    ! angles other than 0 or pi are supported:
    ! 0 < angle < pi
   if (abs(cos_angle).lt.one-EPS2) then 

    if (((SDIM.eq.3).and.(levelrz.eq.COORDSYS_CARTESIAN)).or. &
        ((SDIM.eq.2).and.(levelrz.eq.COORDSYS_RZ))) then
     term1=two/three-cos_angle+(cos_angle**3)/three
     if (term1.gt.zero) then
      ! do nothing
     else
      print *,"term1 invalid"
      stop
     endif
         
     Vtarget=half*(four/three)*Pi*(radblob**3)
     radnew=(Vtarget/(Pi*term1))**(one/three)
     vert=-radnew*cos_angle
    else if ((SDIM.eq.2).and.(levelrz.eq.COORDSYS_CARTESIAN)) then
     test_angle=acos(abs(cos_angle))  ! 0<test_angle<=pi/2
     if (cos_angle.ge.zero) then
      term1=test_angle-half*sin(two*test_angle)
     else
      term1=Pi-test_angle+half*sin(two*test_angle)
     endif
     if (term1.gt.zero) then
      ! do nothing
     else
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
     ! if radblob2=0, then in RZ:
     ! zprime=z-yblob2
     ! yblob2_new=yblob2+vert
     ! if cos_angle>0 then vert<0 (i.e. center of drop shifted down)
     ! if cos_angle<0 then vert>0 (i.e. center of drop shifted up)
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
    else if (maxtall-vert.ge.radnew) then
     ! do nothing
    else
     print *,"maxtall, vert, or radnew invalid"
     stop
    endif 
     
   else
    print *,"contact angle too close to 0 or pi for drop on slope"
    print *,"probtype=",probtype
    print *,"radblob=",radblob
    print *,"radblob2=",radblob2
    print *,"radblob5=",radblob5
    print *,"radblob6=",radblob6
    print *,"radblob7=",radblob7
    stop
   endif

  else
   print *,"parameter conflict for probtype=55"
   print *,"axis_dir = ",axis_dir
   print *,"xblob,yblob,zblob = ",xblob,yblob,zblob
   print *,"xblob2,yblob2,zblob2 = ",xblob2,yblob2,zblob2
   print *,"radblob ",radblob
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

real(amrex_real), INTENT(in) :: Y,WA,WV
real(amrex_real), INTENT(out) :: X

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
subroutine MDOT_Kassemi(sigma,MolarMassFluid,R,Pgamma, &
  Pvapor_probe, &
  Tgamma, &
  Tvapor_probe,MDOT)
IMPLICIT NONE

real(amrex_real), INTENT(in) :: sigma
real(amrex_real), INTENT(in) :: MolarMassFluid
real(amrex_real), INTENT(in) :: R
real(amrex_real), INTENT(in) :: Pgamma
real(amrex_real), INTENT(in) :: Pvapor_probe
real(amrex_real), INTENT(in) :: Tgamma
real(amrex_real), INTENT(in) :: Tvapor_probe
real(amrex_real), INTENT(out) :: MDOT

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

real(amrex_real), INTENT(in) :: PSAT,Tgamma,TSAT,L,R,WV
real(amrex_real), INTENT(out) :: Pgamma

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
  print *,"PSAT invalid in Pgamma_Clausius_Clapyron: ",PSAT
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

real(amrex_real), INTENT(in) :: Tgamma,TSAT,L,R,WV
real(amrex_real), INTENT(out) :: X

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

real(amrex_real), INTENT(in) :: TSAT,L,R,WV
real(amrex_real), INTENT(out) :: X

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

real(amrex_real), INTENT(in) :: TSAT,X,L,R,WV,Tgamma_min,Tgamma_max
real(amrex_real), INTENT(out) :: Tgamma

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

real(amrex_real), INTENT(in) :: X,WA,WV
real(amrex_real), INTENT(out) :: Y

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

integer, INTENT(in) :: imattype,im
real(amrex_real), INTENT(in) :: rho
real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)
real(amrex_real) :: internal_energy
real(amrex_real), INTENT(out) :: pressure
real(amrex_real), INTENT(in) :: internal_energy_in
integer :: ispec


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
 print *,"rho invalid (EOS_material) ",rho
 stop
endif
if (internal_energy.gt.zero) then
 ! do nothing
else
 print *,"e invalid (EOS_material) ",internal_energy
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

! This routine is used by the Tryggvason low Mach number model,
! Xia, Liu, Tryggvason,
! "Fully resolved numerical simulations of fused deposition modeling."
! Part II -
! solidification, residual stresses and modeling of the nozzle". 
! In: Rapid Prototyping Journal 24.6 (2018), pp. 973-987.
! Set ns.material_type_lowmach=0 to disable the Tryggvason low mach
! term.
subroutine dVdT_material(dVdT,massfrac_parm, &
  pressure_in,temperature, &
  imattype,im)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: imattype,im
real(amrex_real), INTENT(in) :: pressure_in,temperature
real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)
real(amrex_real), INTENT(out) :: dVdT
real(amrex_real) :: pressure
integer :: ispec

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

integer, INTENT(in) :: imattype,im
real(amrex_real), INTENT(in) :: rho,temperature
real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)
real(amrex_real), INTENT(out) :: DeDT
real(amrex_real) :: DT,T2,e1,e2


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
 DT=temperature*EPS6
 T2=temperature+DT
 call INTERNAL_material(rho,massfrac_parm, &
   T2,e2,imattype,im)
 DeDT=(e2-e1)/DT
 if (DeDT.gt.zero) then
  ! do nothing
 else
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

integer, INTENT(in) :: imattype,im
real(amrex_real), INTENT(in) :: rho
real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)
real(amrex_real) :: internal_energy
real(amrex_real), INTENT(out) :: soundsqr
real(amrex_real), INTENT(in) :: internal_energy_in
integer :: ispec


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

integer, INTENT(in) :: imattype,im
real(amrex_real), INTENT(in) :: rho
real(amrex_real), INTENT(out) :: internal_energy
real(amrex_real) :: massfrac_parm(num_species_var+1)


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

integer, INTENT(in) :: imattype,im
real(amrex_real), INTENT(in) :: rho,temperature
real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)
real(amrex_real), INTENT(out) :: internal_energy
real(amrex_real) local_internal_energy  ! this is an output

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

real(amrex_real) depth,pgrad,a,b,c,tol
real(amrex_real) surface_den,depth_den
real(amrex_real) surface_pressure,depth_pressure
real(amrex_real) gravity

call fort_derive_gravity(gravity_vector,gravity)

density_at_depth=one

 ! water density where the charge is initially located.
if ((probtype.eq.42).and.(SDIM.eq.2)) then  ! bubble jetting
 density_at_depth=1.00039080D0
  ! cavitation
else if ((probtype.eq.46).and.(SDIM.eq.2)) then 
 if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
  density_at_depth=1.00008343D0
 else if (axis_dir.eq.10) then
  density_at_depth=1.000003 !52x52x120 tank that is pressureized
 else if (axis_dir.eq.11) then
  density_at_depth=one
 else if (axis_dir.eq.20) then
  ! do nothing
 else
  print *,"axis_dir out of range: ",axis_dir
  stop
 endif
else if (fort_material_type(1).eq.13) then

 if (abs(gravity).gt.zero) then
  !do nothing
 else
  print *,"gravity invalid: ",gravity
  stop
 endif

 if (1.eq.0) then
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
 
! print *,"depth= ",depth

 if (depth.gt.zero) then
  !do nothing
 else
  print *,"depth invalid: ",depth
  stop
 endif
 if (surface_den.gt.zero) then
  !do nothing
 else
  print *,"surface den invalid: ",surface_den
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
 tol=abs(gravity)*EPS3
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
integer k1,k2
real(amrex_real) t,u_left

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
integer k1,k2
real(amrex_real) t,u_right

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
integer k1,k2
real(amrex_real) t,elevation_left

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
integer k1,k2
real(amrex_real) t,elevation_right

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

real(amrex_real) x,elevation
real(amrex_real) h,l

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

real(amrex_real) time,x,dist
real(amrex_real) xright,xleft,shallow_tstop
real(amrex_real) thetax,thetat,tgrid,xgrid
real(amrex_real) delta_t_grid,delta_x_grid
real(amrex_real) new_time
integer igrid,jgrid

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

real(amrex_real) time,x,vel
real(amrex_real) xright,xleft,shallow_tstop
real(amrex_real) thetax,thetat,tgrid,xgrid
real(amrex_real) delta_t_grid,delta_x_grid
real(amrex_real) new_time
integer igrid,jgrid

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

! negative on the inside of the triangle!
subroutine triangledist(x,y,xlo,xhi,ylo,yhi,dist)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x,y,xlo,xhi,ylo,yhi
real(amrex_real), INTENT(out) :: dist
real(amrex_real) dist1,dist2,dist3
real(amrex_real) m,b

if ((xlo.ge.xhi-EPS10).or.(ylo.ge.yhi-EPS10)) then 
 print *,"invalid parameters triangle dist",xlo,xhi,ylo,yhi
 stop
endif
dist1=xlo-x
dist2=ylo-y
m=(yhi-ylo)/(xlo-xhi)
b=yhi-m*xlo
dist3=y-(m*x+b)
dist=dist1
if (dist2.gt.dist) then
 dist=dist2
endif
if (dist3.gt.dist) then
 dist=dist3
endif

return
end subroutine triangledist

! negative on the inside of the polygon!
subroutine polygondist(x,y,xlo,xhi,ylo,yhi,xwid,ywid,dist)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x,y,xlo,xhi,ylo,yhi,xwid,ywid
real(amrex_real), INTENT(out) :: dist
real(amrex_real) :: dist1,dist2,dist3
real(amrex_real) dist4,dist5
real(amrex_real) m,b

if ((xlo.ge.xhi-EPS10).or.(ylo.ge.yhi-EPS10)) then 
 print *,"invalid parameters triangle dist",xlo,xhi,ylo,yhi
 stop
endif
dist1=xlo-xwid-x
dist2=ylo-ywid-y
m=(yhi-ylo)/(xlo-xhi)
b=yhi-m*xlo
dist3=y-(m*x+b)
dist4=y-yhi
dist5=x-xhi
dist=dist1
if (dist2.gt.dist) then
 dist=dist2
endif
if (dist3.gt.dist) then
 dist=dist3
endif
if (dist4.gt.dist) then
 dist=dist4
endif
if (dist5.gt.dist) then
 dist=dist5
endif

return
end subroutine polygondist

subroutine zalesakdist(dist,xx,yy)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: xx,yy
real(amrex_real), INTENT(out) :: dist
real(amrex_real) x,y
real(amrex_real) dist1,dist2

x=xx
y=yy
if (axis_dir.eq.0) then
 if (levelrz.eq.COORDSYS_CARTESIAN) then
  ! do nothing
 else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  x=xx*cos(yy)+50.0d0
  y=xx*sin(yy)+50.0d0
 else
  print *,"levelrz invalid in zalesakdist"
  stop
 endif
 dist=sqrt((x-50.0d0)**2+(y-75.0d0)**2)-15.0d0
 if ((x.ge.47.5d0).and.(x.le.52.5d0)) then
  if (y.le.60.0d0) then
   if (x.lt.50.0d0) then
    dist=sqrt( (y-60.0d0)**2+(x-47.5d0)**2 )
   else
    dist=sqrt( (y-60.0d0)**2+(x-52.5d0)**2 )
   endif
  else if (y.le.85.0d0) then
   if (x.lt.50.0d0) then
    dist1=x-47.5d0
   else
    dist1=52.5d0-x
   endif
   dist2=85.0d0-y
   dist=min(dist1,dist2)
  else if ((y.le.90.0d0).and.(dist.le.zero)) then
   dist=max(dist,85.0d0-y)
  endif
 else if ((dist.lt.zero).and.(x.lt.47.5d0)) then
  if (y.le.85.0d0) then
   dist=max(dist,(x-47.5d0))
  else
   dist=max(dist,-sqrt( (x-47.5d0)**2+(y-85.0d0)**2 ) )
  endif
 else if ((dist.lt.zero).and.(x.gt.52.5d0)) then
  if (y.le.85.0d0) then
   dist=max(dist,(52.5d0-x))
  else
   dist=max(dist,-sqrt( (x-52.5d0)**2+(y-85.0d0)**2 ) )
  endif
 endif

else if (axis_dir.eq.1) then
 if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  dist=x-radblob
 else if (levelrz.eq.COORDSYS_CARTESIAN) then
  dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
 else
  print *,"levelrz invalid zalesak dist"
  stop
 endif
else if (axis_dir.eq.2) then
 dist=sqrt((x-50.0d0)**2+(y-75.0d0)**2)-15.0d0
else
 print *,"axis_dir invalid zalesakdist"
 stop
endif
 
return
end subroutine zalesakdist

SUBROUTINE Adist(xx, yy, dist)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT (IN) :: xx
real(amrex_real), INTENT (IN) :: yy
real(amrex_real), INTENT (INOUT) :: dist

real(amrex_real), DIMENSION(4) :: xvec
real(amrex_real), DIMENSION(4) :: yvec
real(amrex_real), DIMENSION(4) :: nx
real(amrex_real), DIMENSION(4) :: ny
real(amrex_real), DIMENSION(4) :: m

real(amrex_real) :: px, py
real(amrex_real) :: vx, vy
real(amrex_real) :: phi_i, maxval
real(amrex_real) :: mhat
real(amrex_real) :: dist1, dist2, dist3
integer :: i

dist = -999.9e9
dist1 = -999.9e9
dist2 = -999.9e9
dist3 = -999.9e9

! Big Triangle !
xvec(1) = 0.0
xvec(2) = 1.0
xvec(3) = -1.0

yvec(1) = 3.0
yvec(2) = 0.0
yvec(3) = 0.0

nx(1) = 3.0/sqrt(10.0)
nx(2) = 0.0
nx(3) = -3.0/sqrt(10.0)

ny(1) = 1.0/sqrt(10.0)
ny(2) = -1.0
ny(3) = 1.0/sqrt(10.0)

m(1) = -3.0
m(2) = 0.0
m(3) = 3.0

maxval = -999.9
do i = 1, 3
	if (i /= 2) then
		!write(*,*) 'i = ', i
		mhat = -1.0/m(i)
		px = (m(i)*xvec(i) - mhat*xx + yy - yvec(i))/(m(i) - mhat)
		py = m(i)*(px - xvec(i)) + yvec(i)
	else	! i == 2, horizontal edge
		!write(*,*) 'i = ', i
		px = xx
		py = yvec(i)
	endif
	
	vx = xx - px
	vy = yy - py
	phi_i = vx*nx(i) + vy*ny(i)
	
	!perr = abs(m(i) - (py - yvec(i))/(px - xvec(i)))
	!write(*,*) 'Projection error = ', perr
	
	if(phi_i > maxval) then
		maxval = phi_i
		!write(*,*) 'i = ', i, 'new phi_i = ', phi_i
	endif
	
	!write(*,*) 'i = ', i, 'nx = ', nx(i), 'ny = ', ny(i), 'vx = ', vx, 'vy = ', vy
enddo

dist1 = maxval
! End Big Triangle !

! Trapezoid !
xvec(1) = 4.0/9.0
xvec(2) = (2.0 + EPS2)/3.0
xvec(3) = -(2.0 + EPS2)/3.0
xvec(4) = -4.0/9.0

yvec(1) = 2.0/3.0
yvec(2) = -EPS2
yvec(3) = -EPS2
yvec(4) = 2.0/3.0

nx(1) = -3.0/sqrt(10.0)
nx(2) = 0.0
nx(3) = 3.0/sqrt(10.0)
nx(4) = 0.0

ny(1) = -1.0/sqrt(10.0)
ny(2) = 1.0
ny(3) = -1.0/sqrt(10.0)
ny(4) = -1.0

m(1) = -3.0
m(2) = 0.0
m(3) = 3.0
m(4) = 0.0

maxval = 999.9
do i = 1, 4
	if ((i == 1) .OR. (i == 3)) then
		!write(*,*) 'i = ', i
		mhat = -1.0/m(i)
		px = (m(i)*xvec(i) - mhat*xx + yy - yvec(i))/(m(i) - mhat)
		py = m(i)*(px - xvec(i)) + yvec(i)
	else	! i == 2 || 4, horizontal edge
		!write(*,*) 'i = ', i
		px = xx
		py = yvec(i)
	endif
	
	vx = xx - px
	vy = yy - py
	phi_i = vx*nx(i) + vy*ny(i)
	
	!perr = abs(m(i) - (py - yvec(i))/(px - xvec(i)))
	!write(*,*) 'Projection error = ', perr
	
	if(phi_i < maxval) then
		maxval = phi_i
		!write(*,*) 'i = ', i, 'new phi_i = ', phi_i
	endif
	
	!write(*,*) 'i = ', i, 'nx = ', nx(i), 'ny = ', ny(i), 'vx = ', vx, 'vy = ', vy
enddo

dist2 = maxval
! End Trapezoid !

! Little Triangle !
xvec(1) = 0.0
xvec(2) = 1.0/3.0
xvec(3) = -1.0/3.0

yvec(1) = 2.0
yvec(2) = 1.0
yvec(3) = 1.0

nx(1) = -3.0/sqrt(10.0)
nx(2) = 0.0
nx(3) = 3.0/sqrt(10.0)

ny(1) = -1.0/sqrt(10.0)
ny(2) = 1.0
ny(3) = -1.0/sqrt(10.0)

m(1) = -3.0
m(2) = 0.0
m(3) = 3.0

maxval = 999.9
do i = 1, 3
	if ((i == 1) .OR. (i == 3)) then
		!write(*,*) 'i = ', i
		mhat = -1.0/m(i)
		px = (m(i)*xvec(i) - mhat*xx + yy - yvec(i))/(m(i) - mhat)
		py = m(i)*(px - xvec(i)) + yvec(i)
	else	! i == 2, horizontal edge
		!write(*,*) 'i = ', i
		px = xx
		py = yvec(i)
	endif
	
	vx = xx - px
	vy = yy - py
	phi_i = vx*nx(i) + vy*ny(i)
	
	!perr = abs(m(i) - (py - yvec(i))/(px - xvec(i)))
	!write(*,*) 'Projection error = ', perr
	
	if(phi_i < maxval) then
		maxval = phi_i
		!write(*,*) 'i = ', i, 'new phi_i = ', phi_i
	endif
	
	!write(*,*) 'i = ', i, 'nx = ', nx(i), 'ny = ', ny(i), 'vx = ', vx, 'vy = ', vy
enddo

dist3 = maxval
! End Little Triangle !

dist = max(dist1, dist2, dist3)

END SUBROUTINE Adist

! Cervone et al 2009, page 416
subroutine deformdist(dist,x,y)
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y
real(amrex_real), INTENT(out) :: dist

dist=sqrt( (x-half)**2 + (y-0.75d0)**2 )-0.15d0

return
end subroutine deformdist

subroutine streamdeform3D(s,x,y,z)
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x,y,z
real(amrex_real), INTENT(out) :: s
real(amrex_real) :: sx,sy,sz

sx=sin(Pi*x)
sy=sin(Pi*y)
sz=sin(Pi*z)

s=sx*sx*sy*sy*sz*sz

return
end subroutine streamdeform3D

subroutine deform3duu(u,x,y,z,t,dx)
IMPLICIT NONE
real(amrex_real), INTENT(in) ::  x,y,z,t
real(amrex_real), INTENT(out) ::  u
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) aa


aa=cos(Pi*t/three)
u=two*(sin(Pi*x)**2)*sin(two*Pi*y)*sin(two*Pi*z)*aa

return
end subroutine deform3duu

subroutine deform3dvv(u,x,y,z,t,dx)
IMPLICIT NONE
real(amrex_real), INTENT(in) ::  x,y,z,t
real(amrex_real), INTENT(out) ::  u
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) aa


aa=cos(Pi*t/three)
u=-(sin(Pi*y)**2)*sin(two*Pi*x)*sin(two*Pi*z)*aa

return
end subroutine deform3dvv

subroutine deform3dww(u,x,y,z,t,dx)
IMPLICIT NONE
real(amrex_real), INTENT(in) ::  x,y,z,t
real(amrex_real), INTENT(out) ::  u
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) aa


aa=cos(Pi*t/three)
u=-(sin(Pi*z)**2)*sin(two*Pi*x)*sin(two*Pi*y)*aa

return 
end subroutine deform3dww

! u = -s_y  
! v = s_x
! u_i+1/2 = -(s_i+1,j+1 + s_i,j+1 - s_i+1,j-1 - s_i,j-1)/(4 dy)
! u_i-1/2 = -(s_i,j+1 + s_i-1,j+1 - s_i,j-1 - s_i-1,j-1)/(4 dy)
! v_j+1/2 = (s_i+1,j+1 + s_i+1,j - s_i-1,j+1 - s_i-1,j)/(4 dx)
! v_j-1/2 = (s_i+1,j + s_i+1,j-1 - s_i-1,j - s_i-1,j-1)/(4 dx)
! div u= -(s_i+1,j+1 + s_i,j+1 - s_i+1,j-1 - s_i,j-1)/(4 dy dx)+
!         (s_i,j+1 + s_i-1,j+1 - s_i,j-1 - s_i-1,j-1)/(4 dy dx)+
!         (s_i+1,j+1 + s_i+1,j - s_i-1,j+1 - s_i-1,j)/(4 dx dy)+
!        -(s_i+1,j + s_i+1,j-1 - s_i-1,j - s_i-1,j-1)/(4 dx dy)=0
! 
subroutine streamdeform(s,x,y)
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x,y
real(amrex_real), INTENT(out) :: s
real(amrex_real) sx,sy

sx=sin(Pi*x)
sy=sin(Pi*y)
s=sx*sx*sy*sy

return
end subroutine streamdeform

subroutine deformuu(u,x,y,t,dx)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y,t
real(amrex_real), INTENT(out) :: u
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) aa,s1,s2,s3,s4
real(amrex_real) x1,x2,y1,y2


if (probtype.ne.29) then
 print *,"probtype should be 29"
 stop
endif

if (fort_stop_time.gt.zero) then
 ! do nothing
else
 print *,"fort_stop_time invalid"
 stop
endif
if (period_time.gt.zero) then
 ! do nothing
else
 print *,"period_time invalid"
 stop
endif
aa=cos(Pi*t/period_time)
if (axis_dir.eq.0) then
 u=sin(four*Pi*(x+half))*sin(four*Pi*(y+half))*aa
else if ((axis_dir.eq.1).or.(axis_dir.eq.3).or. &
         (axis_dir.eq.4)) then
! when t=T, object is back to circular
! Cervone et al 2009, page 413
! psi=(1/pi)sin^2(pi x)sin^2(pi y)
! u=-psi_y=-sin^2(pi x)sin(2 pi y)=-2 sin^2(pi x)sin(pi y)cos(pi y)

 if ((axis_dir.eq.3).or.(axis_dir.eq.1)) then
  s1=sin(Pi*y)
  s2=cos(Pi*y)
  s3=sin(Pi*x)
  u=-two*s1*s2*s3*s3
 else if (axis_dir.eq.4) then
  x1=x+half*dx(1)
  x2=x-half*dx(1)
  y1=y+dx(2)
  y2=y-dx(2)
  call streamdeform(s1,x1,y1)
  call streamdeform(s2,x1,y2)
  call streamdeform(s3,x2,y1)
  call streamdeform(s4,x2,y2)
  u=-(s1-s2+s3-s4)/(four*dx(2)*Pi)
 else
  print *,"bust"
  stop
 endif

 if ((axis_dir.eq.3).or.(axis_dir.eq.4)) then
  u=u*aa
 endif
else if (axis_dir.eq.2) then
 s1=sin(Pi*x)
 s2=sin(two*Pi*y)
 u=s1*s1*s2
 if (t.ge.one) then
  u=-u
 endif
else
 print *,"axis_dir invalid deformuu"
 stop
endif

return
end subroutine deformuu

subroutine deformvv(v,x,y,t,dx)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y,t
real(amrex_real), INTENT(out) :: v
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) aa,s1,s2,s3,s4
real(amrex_real) x1,x2,y1,y2


if (probtype.ne.29) then
 print *,"probtype should be 29"
 stop
endif

if (fort_stop_time.gt.zero) then
 ! do nothing
else
 print *,"fort_stop_time invalid"
 stop
endif
if (period_time.gt.zero) then
 ! do nothing
else
 print *,"period_time invalid"
 stop
endif
aa=cos(Pi*t/period_time)
if (axis_dir.eq.0) then
 v=cos(four*Pi*(x+half))*cos(four*Pi*(y+half))*aa
else if ((axis_dir.eq.1).or.(axis_dir.eq.3).or. &
         (axis_dir.eq.4)) then
! when t=T, object is back to circular
! Cervone et al 2009, page 413
! psi=(1/pi)sin^2(pi x)sin^2(pi y)
! v=psi_x=2 sin(pi x)cos(pi x)sin^2 (pi y)

 if ((axis_dir.eq.3).or.(axis_dir.eq.1)) then
  s1=sin(Pi*x)
  s2=cos(Pi*x)
  s3=sin(Pi*y)
  v=two*s1*s2*s3*s3
 else if (axis_dir.eq.4) then
  x1=x+dx(1)
  x2=x-dx(1)
  y1=y+half*dx(2)
  y2=y-half*dx(2)
  call streamdeform(s1,x1,y1)
  call streamdeform(s2,x2,y1)
  call streamdeform(s3,x1,y2)
  call streamdeform(s4,x2,y2)
  v=(s1-s2+s3-s4)/(four*dx(1)*Pi)
 else
  print *,"bust"
  stop
 endif

 if ((axis_dir.eq.3).or.(axis_dir.eq.4)) then
  v=v*aa
 endif
else if (axis_dir.eq.2) then
 s1=sin(Pi*y)
 s2=sin(two*Pi*x)
 v=-s1*s1*s2
 if (t.ge.one) then
  v=-v
 endif
else
 print *,"axis_dir invalid deformvv"
 stop
endif

return
end subroutine deformvv

subroutine zalesakuu(u,x,y,z,time,dx)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y,z,time
real(amrex_real), INTENT(out) :: u
real(amrex_real), INTENT(in) :: dx(SDIM)

if (probtype.ne.28) then
 print *,"probtype invalid"
 stop
endif
if (SDIM.eq.2) then
 if (abs(z-y).le.EPS2) then
  !do nothing
 else
  print *,"abs(z-y) bust"
  stop
 endif 
endif 

if (adv_vel.eq.zero) then

 if (levelrz.eq.COORDSYS_CARTESIAN) then
  u=-(Pi/314.0)*(y-50.0)
 else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  u=zero
 else
  print *,"zalesakuu: levelrz invalid"
  stop
 endif

else if (adv_vel.ne.zero) then

 if ((adv_dir.eq.1).or.(adv_dir.eq.SDIM+1)) then
  u=adv_vel
 else if (adv_dir.eq.2) then
  u=zero
 else if (adv_dir.eq.SDIM) then
  u=zero
 else
  print *,"adv_dir invalid zalesakuu (7)"
  stop
 endif

else
 print *,"adv_vel invalid zalesakuu: ",adv_vel
 stop
endif

return
end subroutine zalesakuu

subroutine zalesakvv(v,x,y,z,time,dx)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y,z,time
real(amrex_real), INTENT(out) :: v
real(amrex_real), INTENT(in) :: dx(SDIM)

if (probtype.ne.28) then
 print *,"probtype invalid"
 stop
endif
if (SDIM.eq.2) then
 if (abs(z-y).le.EPS2) then
  !do nothing
 else
  print *,"abs(z-y) bust"
  stop
 endif 
endif 

if (adv_vel.eq.zero) then

 if (levelrz.eq.COORDSYS_CARTESIAN) then
  v=(Pi/314.0)*(x-50.0)
 else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  v=(Pi/314.0)*x
 else
  print *,"zalesakvv: levelrz invalid"
  stop
 endif

else if (adv_vel.ne.zero) then

 if ((adv_dir.eq.2).or.(adv_dir.eq.SDIM+1)) then
  v=adv_vel
 else if (adv_dir.eq.1) then
  v=zero
 else if ((adv_dir.eq.SDIM).and.(SDIM.eq.3)) then
  v=zero
 else
  print *,"adv_dir invalid zalesakvv (8)"
  stop
 endif

else
 print *,"adv_vel invalid zalesakvv: ",adv_vel
 stop
endif

return
end subroutine zalesakvv


subroutine zalesakww(w,x,y,z,time,dx)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y,z,time
real(amrex_real), INTENT(out) :: w
real(amrex_real), INTENT(in) :: dx(SDIM)

if (probtype.ne.28) then
 print *,"probtype invalid"
 stop
endif
if (levelrz.ne.COORDSYS_CARTESIAN) then
 print *,"levelrz invalid zalesakvv"
 stop
endif
if (SDIM.eq.2) then
 if (abs(z-y).le.EPS2) then
  !do nothing
 else
  print *,"abs(z-y) bust"
  stop
 endif 
endif 

if (adv_vel.eq.zero) then

 w=zero

else if (adv_vel.ne.zero) then

 if ((adv_dir.eq.SDIM).or.(adv_dir.eq.SDIM+1)) then
  w=adv_vel
 else if (adv_dir.eq.1) then
  w=zero
 else if (adv_dir.eq.2) then
  w=zero
 else
  print *,"adv_dir invalid zalesakww (9)"
  stop
 endif

else
 print *,"adv_vel invalid zalesakww: ",adv_vel
 stop
endif

return
end subroutine zalesakww


subroutine circleuu(u,x,y,z)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y,z
real(amrex_real), INTENT(out) :: u

u=zero
if (adv_dir .eq. 1) then
   u = adv_vel
else if (adv_dir .eq. 2) then
   u = zero
else if (adv_dir.eq.SDIM) then
   u = zero
else if (adv_dir.eq.SDIM+1) then
   u = adv_vel
endif

return
end subroutine circleuu

subroutine circlevv(v,x,y,z)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y,z
real(amrex_real), INTENT(out) :: v

v=zero
if (adv_dir .eq. 2) then
   v = adv_vel
else if (adv_dir .eq. 1) then
   v = zero
else if ((adv_dir.eq.SDIM).and.(SDIM.eq.3)) then
   v = zero
else if (adv_dir.eq.SDIM+1) then
   v = adv_vel
endif

return
end subroutine circlevv

subroutine circleww(w,x,y,z)
use probcommon_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: x,y,z
real(amrex_real), INTENT(out) :: w

if (SDIM.ne.3) then
 print *,"dimension bust circleww"
 stop
endif

w=zero
if ((adv_dir.eq.3).or.(adv_dir.eq.4)) then
   w = adv_vel
else
   w = zero
endif

return
end subroutine circleww

subroutine JohnsonCookSoftening( &
 xsten,nhalf, &
 i,j,k, &
 irefine,jrefine,krefine, &
 level, &
 finest_level, &
 im_critical, &  ! 0<=im_critical<=num_materials-1
 base_yield_stress, &
 T, &
 TM, & !yield_temperature
 T0, & !tempconst
 alpha, &  !fort_yield_alpha
 ref_eps_p, & !ref_plastic_strain
 ref_dot_eps_p, & !ref_plastic_strain_dot
 Johnson_Cook_C, & !fort_Johnson_Cook_C
 hardening_coeff, & ! hardening_coefficient
 eps_p, & !plastic_strain_old
 dot_eps_p, & !dot_plastic_strain
 n, & !yield_n 
 yield_stress)
use probcommon_module
IMPLICIT NONE

integer, intent(in) :: nhalf
real(amrex_real), intent(in) :: xsten(-nhalf:nhalf,SDIM)
integer, intent(in) :: i,j,k
integer, intent(in) :: irefine,jrefine,krefine
integer, intent(in) :: level,finest_level
integer, intent(in) :: im_critical
real(amrex_real), INTENT(in) :: base_yield_stress
real(amrex_real), INTENT(in) :: T,TM,T0,alpha,ref_eps_p,ref_dot_eps_p
real(amrex_real), INTENT(in) :: Johnson_Cook_C
real(amrex_real), INTENT(in) :: hardening_coeff
real(amrex_real), INTENT(in) :: eps_p,dot_eps_p,n
real(amrex_real), INTENT(out) :: yield_stress

if (ref_eps_p.ge.zero) then
 !do nothing
else
 print *,"ref_eps_p invalid ",ref_eps_p
 stop
endif
if (eps_p.ge.zero) then
 !do nothing
else
 print *,"eps_p invalid ",eps_p
 stop
endif

if (ref_dot_eps_p.ge.zero) then
 !do nothing
else
 print *,"ref_dot_eps_p invalid ",ref_dot_eps_p
 stop
endif
if (dot_eps_p.ge.zero) then
 !do nothing
else
 print *,"dot_eps_p invalid ",dot_eps_p
 stop
endif
if (hardening_coeff.ge.zero) then
 !do nothing
else
 print *,"hardening_coeff invalid ",hardening_coeff
 stop
endif
if (Johnson_Cook_C.ge.zero) then
 !do nothing
else
 print *,"Johnson_Cook_C ",Johnson_Cook_C
 stop
endif


if (base_yield_stress.gt.zero) then
 !do nothing
else
 print *,"base_yield_stress invalid: ",base_yield_stress
 stop
endif
if (T.le.TM) then
 !do nothing
else if (T.ge.TM) then
 print *,"T.ge.TM !"
 print *,"T=",T
 print *,"TM=",TM
 print *,"ref_eps_p=",ref_eps_p
 print *,"eps_p=",eps_p
 print *,"i,j,k= ",i,j,k
 print *,"irefine,jrefine,krefine= ",irefine,jrefine,krefine
 print *,"level,finest_level= ",level,finest_level
 print *,"im_critical(0...nmat-1)=",im_critical
 print *,"nhalf ",nhalf
 print *,"xsten(0,1) xsten(0,2) xsten(0,sdim) ",xsten(0,1),xsten(0,2), &
         xsten(0,SDIM)
 stop
else
 print *,"T and/or TM corruption: ",T,TM
 stop
endif
if ((T.gt.zero).and.(TM.gt.zero).and.(T0.gt.zero).and. &
    (alpha.gt.zero).and.(ref_eps_p.gt.zero).and. &
    (ref_dot_eps_p.gt.zero).and. &
    (eps_p.ge.zero).and.(n.gt.zero)) then
 !do nothing
else
 print *,"Johnson Cook Softening parameters corrupt"
 print *,"eps_p ",eps_p
 print *,"ref_eps_p ",ref_eps_p
 print *,"ref_dot_eps_p ",ref_dot_eps_p
 print *,"T ",T
 print *,"TM ",TM
 print *,"T0 ",T0
 print *,"n ",n
 print *,"alpha ",alpha
 stop
endif

yield_stress=base_yield_stress+hardening_coeff*(eps_p**n)
if (dot_eps_p.lt.ref_dot_eps_p) then
 !do nothing
else if (dot_eps_p.ge.ref_dot_eps_p) then
 yield_stress=yield_stress*(one+Johnson_Cook_C*log(dot_eps_p/ref_dot_eps_p))
else
 print *,"dot_eps_p invalid: ",dot_eps_p
 stop
endif

if (T.lt.T0) then
 !do nothing
else if (T.ge.T0) then
 ! note: (13) from Tran and Udaykumar, have 
 ! sigma_y=(A+B(eps^p)^n)(1+C ln(epsdot^p/epsdot^0))(1-theta^m)
 ! see also
 ! Camacho and Ortiz (1997) (31), (33), and just below (33)
 ! Also for heating due to fracture: see (16) for Camacho and Ortiz:
 ! "beta (dot W)^p"
 ! see also:
 ! "Plasticity induced heating in the fracture
 !  and cutting of metals"
 ! From Zehnder, Potdar, and Bhalla:
 ! rho c T_t = div(k grad T) - alpha(3 lambda + 2 mu)T0 eps_dot_kk +
 !   beta sigma_ij eps_dot^p_ij
 ! 3 lambda + 2mu = 3(bulk modulus)
 ! reference [5] from that paper:
 ! Maugin, G.A., "The Thermomechanics of plasticity and fracture"
 ! Cambridge University Press, 1992.
 ! if T=TM then yield_stress=0
 ! note: alpha ~ 1.2
 if (T.ge.TM) then
  yield_stress=zero
 else if ((T.le.TM).and.(T.ge.zero)) then
  yield_stress=yield_stress*(one-((T-T0)/(TM-T0))**alpha) 
 else
  print *,"T or TM invalid: ",T,TM
  stop
 endif
else
 print *,"T invalid: ",T
 stop
endif

return
end subroutine JohnsonCookSoftening

  ! called from fort_updatetensor() and fort_update_particle_tensor()
  ! in GODUNOV_3D.F90
  ! vel is the advective velocity
  ! Plastic algorithm resource: 
  ! An extension of the radial return algorithm to account for rate-dependent
  ! effects in frictional contact and visco-plasticity.
  ! Journal of Materials Processing Technology 
  ! author: Jean-Philippe Ponthot 1998
subroutine point_updatetensor( &
 xsten_in,nhalf_in, &
 i,j,k, &
 irefine,jrefine,krefine, &
 level, &
 finest_level, &
 im_critical, &  ! 0<=im_critical<=num_materials-1
 ncomp_visc, & 
 visc, &
 one_over_den, &
 tendata, & !tendata:fort_getshear,only_scalar=0
 dx,xlo, &
 vel, &
 xmac, &
 ymac, &
 zmac, &
 tnew, &
 told, &
 cell_temperature, &
 tilelo, tilehi,  &
 fablo, fabhi, &
 bfact,  &
 dt, &
 elastic_time, &
 viscoelastic_model, &
 polymer_factor, &
 elastic_viscosity, &
 yield_stress, &
 hardening_coefficient, &
 plastic_work, &
 irz, &
 bc) 

use probcommon_module
IMPLICIT NONE

integer, intent(in) :: nhalf_in
real(amrex_real), intent(in) :: xsten_in(-nhalf_in:nhalf_in,SDIM)
integer, INTENT(in) :: i,j,k
integer, INTENT(in) :: irefine,jrefine,krefine
integer, INTENT(in) :: level
integer, INTENT(in) :: finest_level
integer, INTENT(in) :: im_critical
integer, INTENT(in) :: ncomp_visc
integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
integer, INTENT(in) :: bfact
real(amrex_real), INTENT(in) :: dx(SDIM),xlo(SDIM)

real(amrex_real), INTENT(in), pointer :: visc(D_DECL(:,:,:),:)
real(amrex_real), INTENT(in), pointer :: one_over_den(D_DECL(:,:,:))

! D=(1/2)(gradU + gradU^Transpose)
! DERIVE_TENSOR_MAG+1: sqrt(2 * D : D)
! DERIVE_TENSOR_RATE_DEFORM+1: D11,D12,D13,D21,D22,D23,D31,D32,D33
! DERIVE_TENSOR_GRAD_VEL+1: ux,uy,uz,vx,vy,vz,wx,wy,wz
! grad u = | ux vx wx |  (grad u)_ij=u_{j,i}
!          | uy vy wy |
!          | uz vz wz |
! (\partial u)/(\partial x)_ij=(grad u)^T=u_{i,j}
!
real(amrex_real), INTENT(in), pointer :: tendata(D_DECL(:,:,:),:)
real(amrex_real), INTENT(in), pointer :: vel(D_DECL(:,:,:),:)
real(amrex_real), INTENT(in), pointer :: xmac(D_DECL(:,:,:))
real(amrex_real), INTENT(in), pointer :: ymac(D_DECL(:,:,:))
real(amrex_real), INTENT(in), pointer :: zmac(D_DECL(:,:,:))

real(amrex_real), INTENT(out) :: tnew(ENUM_NUM_TENSOR_TYPE)
real(amrex_real), INTENT(in) :: told(ENUM_NUM_TENSOR_TYPE)
real(amrex_real), INTENT(in) :: cell_temperature

integer :: n
real(amrex_real), INTENT(in) :: dt,elastic_time
integer, INTENT(in) :: viscoelastic_model
real(amrex_real), INTENT(in) :: polymer_factor
real(amrex_real), INTENT(in) :: elastic_viscosity
real(amrex_real), INTENT(in) :: yield_stress
real(amrex_real), INTENT(in) :: hardening_coefficient
integer, INTENT(in) :: bc(SDIM,2,SDIM)
integer, INTENT(in) :: irz
integer ii,jj,kk
integer iofs,jofs,kofs
integer iofs2,jofs2,kofs2
real(amrex_real) visctensor(3,3) ! derived from CC velocity
real(amrex_real) visctensorMAC(3,3) !derived from MAC velocity
  !deviatoric strain-rate
real(amrex_real) visctensorMAC_traceless(3,3) !derived from MAC velocity
real(amrex_real) trace_MAC
integer trace_dim
real(amrex_real) gradV_MAC
real(amrex_real) gradV_transpose_FENECR(3,3)
 !(grad V)_ji=(partial V/partial x)_ij=v_i,j
real(amrex_real) gradV_transpose(3,3) 
real(amrex_real) Q(3,3)
! W=(1/2)(grad V - (grad V)^T)  derived from the MAC velocity
real(amrex_real) W_Jaumann(3,3)  
real(amrex_real) Aadvect(3,3)
real(amrex_real) Smult_left(3,3)
real(amrex_real) Smult_right(3,3)
real(amrex_real) SA(3,3)
real(amrex_real) SAS(3,3)
real(amrex_real) NP(3,3)
real(amrex_real) shear
real(amrex_real) modtime,trace_A
real(amrex_real) equilibrium_diagonal
real(amrex_real) inverse_tol
real(amrex_real) dui(3,3)
real(amrex_real) dxj(3,3)
real(amrex_real) signcoeff
real(amrex_real) weightcoeff
real(amrex_real) total_weight
real(amrex_real) weight_prev

integer dir_local

integer DErelaxation_model !Deborah number defines rate of relaxation
real(amrex_real) magA,NP_dotdot_D,Y_plastic_parm_scaled,f_plastic
real(amrex_real) gamma_not
real(amrex_real) force_coef
real(amrex_real) one_over_den_local
real(amrex_real) r_hoop

real(amrex_real) plastic_strain_old,plastic_strain_dot
real(amrex_real), intent(out) :: plastic_work

integer Johnson_iter
integer force_unity_determinant
integer unity_det

integer, parameter :: nhalf=3
real(amrex_real) xsten(-nhalf:nhalf,SDIM)

if (irz.ne.levelrz) then
 print *,"irz invalid"
 stop
endif
if (bfact.lt.1) then
 print *,"bfact invalid60"
 stop
endif
if ((level.lt.0).or.(level.gt.finest_level)) then
 print *,"level invalid 34"
 stop
endif

if (polymer_factor.ge.zero) then !1/L
 ! do nothing
else
 print *,"polymer_factor out of range"
 stop
endif
if (elastic_viscosity.gt.zero) then 
 ! do nothing
else
 print *,"elastic_viscosity out of range"
 stop
endif
if (yield_stress.gt.zero) then 
 ! do nothing
else
 print *,"yield_stress out of range: ",yield_stress
 stop
endif
if (hardening_coefficient.ge.zero) then 
 ! do nothing
else
 print *,"hardening_coefficient out of range: ",hardening_coefficient
 stop
endif
if (cell_temperature.gt.zero) then 
 ! do nothing
else
 print *,"cell_temperature out of range: ",cell_temperature
 stop
endif

if (fort_material_type(im_critical+1).eq.0) then
 force_unity_determinant=1
else if (fort_material_type(im_critical+1).eq.999) then
 print *,"999 is unexpected material_type"
 print *,"im_critical,fort_material_type: ",im_critical, &
    fort_material_type(im_critical+1)
 stop
else if ((fort_material_type(im_critical+1).gt.0).and. &
         (fort_material_type(im_critical+1).le.MAX_NUM_EOS)) then
 force_unity_determinant=0
else
 print *,"unexpected material_type"
 print *,"im_critical,fort_material_type: ",im_critical, &
    fort_material_type(im_critical+1)
 stop
endif

if (viscoelastic_model.eq.NN_FENE_CR) then ! FENE-CR
 ! coeff=(visc-etaS)/(modtime+dt)
 ! modtime=max(0.0,elastic_time*(1-Tr(A)/L^2))
 force_unity_determinant=0
else if (viscoelastic_model.eq.NN_OLDROYD_B) then ! Oldroyd-B
 ! coeff=(visc-etaS)/(modtime+dt)
 ! modtime=elastic_time
 force_unity_determinant=0
else if (viscoelastic_model.eq.NN_FENE_P) then ! FENE-P
 ! coeff=(visc-etaS)/(modtime+dt)
 ! modtime=max(0.0,elastic_time*(1-Tr(A)/L^2))
 force_unity_determinant=0
else if (viscoelastic_model.eq.NN_LINEAR_PTT) then ! linear PTT
 ! coeff=(visc-etaS)/(modtime+dt)
 ! modtime=elastic_time
 force_unity_determinant=0
else if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then !incremental model
 ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
 ! coeff=elastic_viscosity
 force_unity_determinant=0
else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental 
 ! Xia, Lu, Tryggvason 2018
 ! coeff=elastic_viscosity
 force_unity_determinant=0
else
 print *,"viscoelastic_model invalid: ",viscoelastic_model
 stop
endif

if (levelrz.eq.COORDSYS_CARTESIAN) then
 !do nothing
else if (levelrz.eq.COORDSYS_RZ) then
 !do nothing
else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
 !do nothing
else
 print *,"levelrz invalid: ",levelrz
 stop
endif

if ((im_critical.lt.0).or. &
    (im_critical.ge.num_materials)) then
 print *,"im_critical invalid27: ",im_critical
 stop
endif
if (ncomp_visc.ne.3*num_materials) then
 print *,"ncomp visc invalid: ",num_materials,ncomp_visc
 stop
endif
if (dt.gt.zero) then
 ! do nothing
else
 print *,"dt invalid: ",dt
 stop
endif

call checkbound_array(fablo,fabhi,visc,0,-1)
call checkbound_array1(fablo,fabhi,one_over_den,0,-1)
call checkbound_array(fablo,fabhi,tendata,0,-1)
call checkbound_array(fablo,fabhi,vel,1,-1)
call checkbound_array1(fablo,fabhi,xmac,1,0)
call checkbound_array1(fablo,fabhi,ymac,1,1)
call checkbound_array1(fablo,fabhi,zmac,1,SDIM-1)

if ((viscoelastic_model.eq.NN_FENE_CR).or. & !FENE-CR
    (viscoelastic_model.eq.NN_OLDROYD_B).or. & !OLDROYD-B
    (viscoelastic_model.eq.NN_FENE_P).or. & !FENE-P
    (viscoelastic_model.eq.NN_LINEAR_PTT)) then !linear PTT
 DErelaxation_model=1
else if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental model
 ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
 DErelaxation_model=0
else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental 
 ! Xia, Lu, Tryggvason 2018
 ! Seungwon Shin, Jalel Chergui, Damir Juric
 DErelaxation_model=0
else
 print *,"viscoelastic_model invalid: ",viscoelastic_model
 stop
endif

call gridsten_level(xsten,i,j,k,level,nhalf)
 ! (grad u)_{ji}=
 ! (grad u)^T)_{ij}=
 !               | u_r  u_t/r-v/r  u_z  |=(\partial u)/(\partial x)=u_{i,j}
 !               | v_r  v_t/r+u/r  v_z  |
 !               | w_r  w_t/r      w_z  |
 ! if levelrz==COORDSYS_RZ,
 !  gradU(3,3)=u/|r|
 ! if levelrz==COORDSYS_CYLINDRICAL,
 !  gradU(2,2)+=u/|r|
 !  (gradU)^T(1,2)-=v/|r|
 ! tendata is intialized in: fort_getshear
 ! D=(1/2)(gradU + gradU^Transpose)
 ! tendata has: |D|, D, grad U
 ! DERIVE_TENSOR_MAG+1: sqrt(2 * D : D)   
 ! DERIVE_TENSOR_RATE_DEFORM+1: D11,D12,D13,D21,D22,D23,D31,D32,D33
 ! DERIVE_TENSOR_GRAD_VEL+1: ux,uy,uz,vx,vy,vz,wx,wy,wz
 ! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
 !         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
 !         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |
shear=tendata(D_DECL(i,j,k),DERIVE_TENSOR_MAG+1) ! sqrt(2 D:D)
n=DERIVE_TENSOR_RATE_DEFORM+1
do ii=1,3
do jj=1,3
 ! (1/2) (grad U + (grad U)^T)
 visctensor(ii,jj)=tendata(D_DECL(i,j,k),n)
 n=n+1
enddo 
enddo 

if (n.eq.DERIVE_TENSOR_GRAD_VEL+1) then
 ! do nothing
else
 print *,"n.eq.DERIVE_TENSOR_GRAD_VEL+1 failed"
 stop
endif

do ii=1,3 !veldir
do jj=1,3 !deriv dir

 dui(ii,jj)=zero
 dxj(ii,jj)=zero

 total_weight=zero

 kofs=0
#if (AMREX_SPACEDIM==3)
 do kofs=0,1
#endif
 do jofs=0,1
 do iofs=0,1
   
  iofs2=iofs
  jofs2=jofs
  kofs2=kofs

  if (ii.eq.1) then !veldir (ii) xmac (umac)

   if ((jj.ge.1).and.(jj.le.SDIM)) then

    if (jrefine.eq.0) then
     jofs2=jofs2-1
    else if (jrefine.eq.1) then
     !do nothing
    else
     print *,"jrefine invalid: ",jrefine
     stop
    endif
#if (AMREX_SPACEDIM==3)
    if (krefine.eq.0) then
     kofs2=kofs2-1
    else if (krefine.eq.1) then
     !do nothing
    else
     print *,"krefine invalid: ",krefine
     stop
    endif
#endif

    signcoeff=one
    weightcoeff=one

    if (jj.eq.1) then !du/dx

     if (jofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(jofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"jofs2 invalid"
      stop
     endif
#if (AMREX_SPACEDIM==3)
     if (kofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(kofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"kofs2 invalid"
      stop
     endif
#endif
     if (iofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (iofs.eq.1) then
      signcoeff=one
     else
      print *,"iofs invalid ",iofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*iofs2-1,jj)

    else if (jj.eq.2) then !du/dy

     if (iofs2.eq.irefine) then
      weightcoeff=weightcoeff*0.75d0
     else if (iofs2.eq.1-irefine) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"iofs2 invalid"
      stop
     endif
#if (AMREX_SPACEDIM==3)
     if (kofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(kofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"kofs2 invalid"
      stop
     endif
#endif
     if (jofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (jofs.eq.1) then
      signcoeff=one
     else
      print *,"jofs invalid ",jofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*jofs2,jj)

    else if ((jj.eq.3).and.(SDIM.eq.3)) then !du/dz

     if (iofs2.eq.irefine) then
      weightcoeff=weightcoeff*0.75d0
     else if (iofs2.eq.1-irefine) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"iofs2 invalid"
      stop
     endif
     if (jofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(jofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"jofs2 invalid"
      stop
     endif
     if (kofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (kofs.eq.1) then
      signcoeff=one
     else
      print *,"kofs invalid ",kofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*kofs2,jj)

    else
     print *,"jj invalid"
     stop
    endif

    dui(ii,jj)=dui(ii,jj)+ &
     weightcoeff*signcoeff*xmac(D_DECL(i+iofs2,j+jofs2,k+kofs2))

   else if ((jj.eq.3).and.(SDIM.eq.2)) then
    !do nothing
   else
    print *,"jj invalid(b): ",jj
    stop
   endif 

  else if (ii.eq.2) then !veldir="v" (ii) ymac (vmac)

   if ((jj.ge.1).and.(jj.le.SDIM)) then

    if (irefine.eq.0) then
     iofs2=iofs2-1
    else if (irefine.eq.1) then
     !do nothing
    else
     print *,"irefine invalid: ",irefine
     stop
    endif
#if (AMREX_SPACEDIM==3)
    if (krefine.eq.0) then
     kofs2=kofs2-1
    else if (krefine.eq.1) then
     !do nothing
    else
     print *,"krefine invalid: ",krefine
     stop
    endif
#endif

    signcoeff=one
    weightcoeff=one

    if (jj.eq.1) then !dv/dx

     if (jofs2.eq.jrefine) then
      weightcoeff=weightcoeff*0.75d0
     else if (jofs2.eq.1-jrefine) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"jofs2 invalid"
      stop
     endif
#if (AMREX_SPACEDIM==3)
     if (kofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(kofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"kofs2 invalid"
      stop
     endif
#endif
     if (iofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (iofs.eq.1) then
      signcoeff=one
     else
      print *,"iofs invalid ",iofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*iofs2,jj)

    else if (jj.eq.2) then !dv/dy

     if (iofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(iofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"iofs2 invalid"
      stop
     endif
#if (AMREX_SPACEDIM==3)
     if (kofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(kofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"kofs2 invalid"
      stop
     endif
#endif
     if (jofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (jofs.eq.1) then
      signcoeff=one
     else
      print *,"jofs invalid ",jofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*jofs2-1,jj)

    else if ((jj.eq.3).and.(SDIM.eq.3)) then !dv/dz

     if (jofs2.eq.jrefine) then
      weightcoeff=weightcoeff*0.75d0
     else if (jofs2.eq.1-jrefine) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"jofs2 invalid"
      stop
     endif
     if (iofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(iofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"iofs2 invalid"
      stop
     endif
     if (kofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (kofs.eq.1) then
      signcoeff=one
     else
      print *,"kofs invalid ",kofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*kofs2,jj)

    else
     print *,"jj invalid"
     stop
    endif

    dui(ii,jj)=dui(ii,jj)+ &
     weightcoeff*signcoeff*ymac(D_DECL(i+iofs2,j+jofs2,k+kofs2))

   else if ((jj.eq.3).and.(SDIM.eq.2)) then
    !do nothing
   else
    print *,"jj invalid(b): ",jj
    stop
   endif 

  else if ((ii.eq.SDIM).and.(SDIM.eq.3)) then ! zmac (wmac) (ii)

   if ((jj.ge.1).and.(jj.le.SDIM)) then

    if (irefine.eq.0) then
     iofs2=iofs2-1
    else if (irefine.eq.1) then
     !do nothing
    else
     print *,"irefine invalid: ",irefine
     stop
    endif
    if (jrefine.eq.0) then
     jofs2=jofs2-1
    else if (jrefine.eq.1) then
     !do nothing
    else
     print *,"jrefine invalid: ",jrefine
     stop
    endif

    signcoeff=one
    weightcoeff=one

    if (jj.eq.1) then !dw/dx

     if (kofs2.eq.krefine) then
      weightcoeff=weightcoeff*0.75d0
     else if (kofs2.eq.1-krefine) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"kofs2 invalid"
      stop
     endif
     if (jofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(jofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"jofs2 invalid"
      stop
     endif
     if (iofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (iofs.eq.1) then
      signcoeff=one
     else
      print *,"iofs invalid ",iofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*iofs2,jj)

    else if (jj.eq.2) then !dw/dy

     if (iofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(iofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"iofs2 invalid"
      stop
     endif
     if (kofs2.eq.krefine) then
      weightcoeff=weightcoeff*0.75d0
     else if (kofs2.eq.1-krefine) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"kofs2 invalid"
      stop
     endif
     if (jofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (jofs.eq.1) then
      signcoeff=one
     else
      print *,"jofs invalid ",jofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*jofs2,jj)

    else if (jj.eq.3) then !dw/dz

     if (jofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(jofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"jofs2 invalid"
      stop
     endif
     if (iofs2.eq.0) then
      weightcoeff=weightcoeff*0.75d0
     else if (abs(iofs2).eq.1) then
      weightcoeff=weightcoeff*0.25d0
     else
      print *,"iofs2 invalid"
      stop
     endif
     if (kofs.eq.0) then
      signcoeff=-one
      total_weight=total_weight+weightcoeff
     else if (kofs.eq.1) then
      signcoeff=one
     else
      print *,"kofs invalid ",kofs
      stop
     endif
     dxj(ii,jj)=dxj(ii,jj)+weightcoeff*signcoeff*xsten(2*kofs2-1,jj)

    else
     print *,"jj invalid"
     stop
    endif

    dui(ii,jj)=dui(ii,jj)+ &
     weightcoeff*signcoeff*zmac(D_DECL(i+iofs2,j+jofs2,k+kofs2))

   else if ((jj.eq.3).and.(SDIM.eq.2)) then
    !do nothing
   else
    print *,"jj invalid(b): ",jj
    stop
   endif 

  else if ((ii.eq.3).and.(SDIM.eq.2)) then
   !do nothing
  else
   print *,"ii invalid: ",ii
   stop
  endif

 enddo !iofs=0,1
 enddo !jofs=0,1
#if (AMREX_SPACEDIM==3)
 enddo !kofs=0,1
#endif
 
 gradV_MAC=zero
 if ((ii.ge.1).and.(ii.le.SDIM).and. &
     (jj.ge.1).and.(jj.le.SDIM)) then

  if (abs(total_weight-one).le.EPS10) then
   !do nothing
  else
   print *,"total_weight invalid: ",total_weight
   stop
  endif

  if (dxj(ii,jj).gt.zero) then
   gradV_MAC=dui(ii,jj)/dxj(ii,jj)  !v_i,j
  else
   print *,"dxj(ii,jj) invalid: ",ii,jj,dxj(ii,jj)
   stop
  endif

 else if ((ii.eq.3).or.(jj.eq.3)) then
  ! do nothing
 else
  print *,"ii or jj invalid: ",ii,jj
  stop
 endif

 r_hoop=xsten(0,1)
 if (levelrz.eq.COORDSYS_RZ) then
  if (SDIM.eq.2) then
   !do nothing
  else
   print *,"dimension bust"
   stop
  endif
  if (r_hoop.gt.zero) then
   !do nothing
  else
   print *,"r_hoop invalid: ",r_hoop
   stop
  endif
  if ((ii.eq.3).and.(jj.eq.3)) then
   gradV_MAC=gradV_MAC+half*(xmac(D_DECL(i,j,k))+ &
     xmac(D_DECL(i+1,j,k)))/r_hoop
  endif
 else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  if (r_hoop.gt.zero) then
   !do nothing
  else
   print *,"r_hoop invalid: ",r_hoop
   stop
  endif
  if (jj.eq.2) then
   if (ii.eq.1) then
    gradV_MAC=gradV_MAC-half*(ymac(D_DECL(i,j,k))+ &
     ymac(D_DECL(i,j+1,k)))
   else if (ii.eq.2) then
    gradV_MAC=gradV_MAC+half*(xmac(D_DECL(i,j,k))+ &
     xmac(D_DECL(i+1,j,k)))
   else if (ii.eq.3) then
    !do nothing
   else
    print *,"ii invalid: ",ii
    stop
   endif
   gradV_MAC=gradV_MAC/r_hoop
  endif
 else if (levelrz.eq.COORDSYS_CARTESIAN) then
  !do nothing
 else
  print *,"levelrz invalid: ",levelrz
  stop
 endif
     
!gradV_transpose_{i,j}=(\partial V_{i})/(\partial x_{j})
!gradV_transpose(ii,jj)=tendata(D_DECL(i,j,k),n) !(vel dir,deriv dir)

 !((grad V)^{T})_{ij}=
 !(grad V)_ji=(partial V/partial x)_ij=v_i,j
 gradV_transpose(ii,jj)=gradV_MAC !(vel dir,deriv dir)  v_i,j

 if ((viscoelastic_model.eq.NN_FENE_CR).or. & !FENE-CR
     (viscoelastic_model.eq.NN_OLDROYD_B).or. & !OLDROYD-B
     (viscoelastic_model.eq.NN_FENE_P).or. & !FENE-P
     (viscoelastic_model.eq.NN_LINEAR_PTT)) then !linear PTT
  !gradV_transpose
  gradV_transpose_FENECR(ii,jj)=gradV_transpose(ii,jj)
 else if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then !incremental model
  !stub
  gradV_transpose_FENECR(ii,jj)=gradV_transpose(ii,jj)
 else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental 
  gradV_transpose_FENECR(ii,jj)=gradV_transpose(ii,jj)
 else
  print *,"viscoelastic_model invalid: ",viscoelastic_model
  stop
 endif

 n=n+1
enddo !jj=1,3
enddo !ii=1,3

if (n.eq.DERIVE_TENSOR_NCOMP+1) then
 ! do nothing
else
 print *,"n.eq.DERIVE_TENSOR_NCOMP+1 failed: ",n
 stop
endif

trace_MAC=zero

do ii=1,3
do jj=1,3
  !W=(grad V - grad V^T)/2
  !W_{ij}=(V_{j,i}-V_{i,j})/2=-Omega_ij  (see Tran and Udaykumar)
 W_Jaumann(ii,jj)=half*(gradV_transpose(jj,ii)-gradV_transpose(ii,jj))
 visctensorMAC(ii,jj)=half*(gradV_transpose(jj,ii)+gradV_transpose(ii,jj))
 visctensorMAC_traceless(ii,jj)=visctensorMAC(ii,jj)
 if (ii.eq.jj) then
  trace_MAC=trace_MAC+visctensorMAC(ii,jj)
 endif
enddo 
enddo 

!hydrostatic strain rate
if ((SDIM.eq.2).and.(levelrz.eq.COORDSYS_RZ)) then
 trace_MAC=trace_MAC/3.0d0
 trace_dim=3
else if ((SDIM.eq.2).and.(levelrz.eq.COORDSYS_CYLINDRICAL)) then
 trace_MAC=trace_MAC/3.0d0
 trace_dim=3
else if ((SDIM.eq.2).and.(levelrz.eq.COORDSYS_CARTESIAN)) then
 trace_MAC=trace_MAC/2.0d0
 trace_dim=2
else if (SDIM.eq.3) then
 trace_MAC=trace_MAC/3.0d0
 trace_dim=3
else
 print *,"sdim or levelrz invalid: ",SDIM,levelrz
 stop
endif

do ii=1,trace_dim
 visctensorMAC_traceless(ii,ii)=visctensorMAC_traceless(ii,ii)-trace_MAC
enddo
 
 
do ii=1,3
do jj=1,3
 Q(ii,jj)=zero
enddo
enddo
do dir_local=1,ENUM_NUM_TENSOR_TYPE_BASE
 call stress_index(dir_local,ii,jj)
 Q(ii,jj)=told(dir_local)
enddo

plastic_strain_old=told(ENUM_NUM_TENSOR_TYPE_BASE+1)
plastic_strain_dot=zero
plastic_work=zero

Q(2,1)=Q(1,2)
Q(3,1)=Q(1,3)
Q(3,2)=Q(2,3)

 ! modtime=lambda/f(A)
modtime=visc(D_DECL(i,j,k),2*num_materials+im_critical+1)
if (modtime.ge.zero) then
 ! do nothing
else
 print *,"modtime invalid: ",modtime
 stop
endif

! e.g. visc_coef * elastic_viscosity/(modtime+dt)
force_coef=visc(D_DECL(i,j,k),num_materials+im_critical+1)

if (force_coef.gt.zero) then
 ! do nothing
else
 print *,"expecting force_coef>0: ",force_coef
 print *,"im_critical=",im_critical
 print *,"fort_visc_coef=",fort_visc_coef
 stop
endif

one_over_den_local=one_over_den(D_DECL(i,j,k))

if (one_over_den_local.gt.zero) then
 ! do nothing
else
 print *,"one_over_den_local invalid: ",one_over_den_local
 stop
endif

trace_A=zero
do ii=1,3
 trace_A=trace_A+Q(ii,ii)+one
 !NN_FENE_CR,NN_OLDROYD_B,NN_FENE_P,NN_LINEAR_PTT
 if (DErelaxation_model.eq.1) then
  if (Q(ii,ii)+one.gt.zero) then
   ! do nothing
  else
   print *,"A=Q+I should be positive definite"
   print *,"ii,Q(ii,ii) ",ii,Q(ii,ii)
   stop
  endif
  if (Q(3,3).gt.-one) then !hoop term
   ! do nothing
  else
   print *,"Q(3,3) invalid (hoop): ",Q(3,3)
   stop
  endif
 else if (DErelaxation_model.eq.0) then ! e.g. incremental model
  ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
  !   or
  ! Xia, Lu, Tryggvason 2018

  if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental model
   ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
   ! do nothing
  else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental 
   ! Xia, Lu, Tryggvason 2018
   if (Q(ii,ii)+one.gt.zero) then
    ! do nothing
   else
    print *,"A=Q+I should be positive definite"
    print *,"ii,Q(ii,ii) ",ii,Q(ii,ii)
    stop
   endif
   if (Q(3,3).gt.-one) then !hoop term
    ! do nothing
   else
    print *,"Q(3,3) invalid (hoop): ",Q(3,3)
    stop
   endif
  else
   print *,"viscoelastic_model invalid: ",viscoelastic_model
   stop
  endif

 else
  print *,"DErelaxation_model invalid: ",DErelaxation_model
  stop
 endif
enddo ! ii=1,3

 ! f(A)=1/(1-trace(A)/L^2)
 ! (1/lambda)*(f(A)A-I)=
 ! (f(A)/lambda)*(A-I/f(A))=
 ! (f(A)/lambda)*(Q+I-I/f(A))
 ! 1-1/f(A)=1-(1-trac(A)/L^2)=trac(A)/L^2
 ! (1/lambda)*(f(A)A-I)=(f(A)/lambda)*(Q+trac(A)I/L^{2})
if (viscoelastic_model.eq.NN_FENE_P) then !FENE-P          
 if (trace_A.gt.zero) then
  if (polymer_factor.gt.zero) then !1/L
   equilibrium_diagonal=min(trace_A*(polymer_factor**2),one)
  else
   print *,"polymer_factor out of range for FENE-P: ",polymer_factor
   stop
  endif
 else
  print *,"trace_A should be positive for FENE-P"
  stop
 endif
else if (viscoelastic_model.eq.NN_FENE_CR) then !FENE-CR
 equilibrium_diagonal=zero
else if (viscoelastic_model.eq.NN_OLDROYD_B) then !OLDROYD-B
 equilibrium_diagonal=zero
else if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then !incremental model
 ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
 equilibrium_diagonal=zero
else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental 
 ! Xia, Lu, Tryggvason 2018
 equilibrium_diagonal=zero
else if (viscoelastic_model.eq.NN_LINEAR_PTT) then !linearPTT
 equilibrium_diagonal=zero
 if (trace_A.gt.zero) then
  if (polymer_factor.ge.zero) then 
   if (SDIM*(polymer_factor**2).lt.one) then
    modtime=modtime/(one+(polymer_factor**2)*(trace_A-three))
   else
    print *,"eps=polymer_factor**2"
    print *,"need eps * sdim < 1"
    stop
   endif
  else
   print *,"polymer_factor invalid: ",polymer_factor
   stop
  endif
 else
  print *,"trace_A invalid: ",trace_A
  stop
 endif
else
 print *,"viscoelastic_model invalid: ",viscoelastic_model
 stop
endif

if ((viscoelastic_model.eq.NN_FENE_CR).or. & !FENE-CR
    (viscoelastic_model.eq.NN_OLDROYD_B).or. & !OLDROYD_B
    (viscoelastic_model.eq.NN_FENE_P).or. & !FENE-P
    (viscoelastic_model.eq.NN_LINEAR_PTT).or. & !linear PTT
    (viscoelastic_model.eq.NN_NEO_HOOKEAN).or. & !incremental Neo-Hookean
    (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL)) then !incremental

 do ii=1,3 
 do jj=1,3 

  Aadvect(ii,jj)=Q(ii,jj)

   !cfl cond: |u|dt<dx and dt|gradu|<1
  if (DErelaxation_model.eq.1) then
   Smult_left(ii,jj)=dt*gradV_transpose_FENECR(ii,jj) 
   Smult_right(ii,jj)=Smult_left(ii,jj)
  else if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then !incremental
   if (DErelaxation_model.eq.0) then
     !W_{ij}=(V_{j,i}-V_{i,j})/2=-OMEGA (see Tran and Udaykumar)
     !(I-dt W)Q(I-dt W^T)=Q-dt WQ-dt QW^T=Q-dt WQ+dt QW
    Smult_left(ii,jj)=-dt*W_Jaumann(ii,jj) ! dt * OMEGA
    Smult_right(ii,jj)=Smult_left(ii,jj)   ! dt * OMEGA
   else
    print *,"DErelaxation_model invalid: ",DErelaxation_model
    stop
   endif
  else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then !incremental Neo-Hookean
   ! Xia, Lu, Tryggvason 2018
   ! Seungwon Shin, Jalel Chergui, Damir Juric
   ! (grad U)^T_{ki}=(partial u_k)/(partial x_i)
   ! Df/Dt + f (grad U)^T=0  Left Cauchy Green tensor B=F F^T=(f^T f)^{-1}
   ! D(f^T f)/Dt=f^T Df/Dt + Df^T/Dt f =
   !             f^T(-f grad U^T)+(-grad U f^T)f  
   ! let Binv=f^T f
   ! D Binv/Dt + Binv grad U^T + grad U Binv = 0
   ! D (Binv B)/Dt=D Binv/Dt B + Binv DB/Dt=
   ! (-Binv grad U^T - grad U Binv)B + Binv DB/Dt = 0
   ! -(grad U^T)B-B grad U + DB/Dt = 0
   ! DB/Dt = (grad U)^T B + B(grad U)
   ! equilibrium is B=I
   ! discretely, B should maintain as positive definite:
   ! B^n+1 = (I+dt grad U^T)Bstar(I+dt grad U)
   if (DErelaxation_model.eq.0) then
    Smult_left(ii,jj)=dt*gradV_transpose_FENECR(ii,jj) ! dt v_{i,j}
    Smult_right(ii,jj)=Smult_left(ii,jj)
   else
    print *,"DErelaxation_model invalid: ",DErelaxation_model
    stop
   endif
  else
   print *,"viscoelastic_model invalid: ",viscoelastic_model
   stop
  endif

  inverse_tol=0.05d0

   !if the CFL condition is satified, then we expect
   !dt |gradu| <=1
  if (Smult_left(ii,jj).le.-one+inverse_tol) then
   Smult_left(ii,jj)=-one+inverse_tol
  else if (Smult_left(ii,jj).ge.one-inverse_tol) then
   Smult_left(ii,jj)=one-inverse_tol
  else if (abs(Smult_left(ii,jj)).lt.one) then
   ! do nothing
  else
   print *,"Smult_left(ii,jj) became corrupt: ",ii,jj, &
           Smult_left(ii,jj)
   stop
  endif

  if (Smult_right(ii,jj).le.-one+inverse_tol) then
   Smult_right(ii,jj)=-one+inverse_tol
  else if (Smult_right(ii,jj).ge.one-inverse_tol) then
   Smult_right(ii,jj)=one-inverse_tol
  else if (abs(Smult_right(ii,jj)).lt.one) then
   ! do nothing
  else
   print *,"Smult_right(ii,jj) became corrupt: ",ii,jj, &
           Smult_right(ii,jj)
   stop
  endif

 enddo ! jj=1,3
 enddo ! ii=1,3

 do ii=1,3

  Smult_left(ii,ii)=Smult_left(ii,ii)+one
  Smult_right(ii,ii)=Smult_right(ii,ii)+one

   ! Aadvect <-- Q+I
  if (DErelaxation_model.eq.1) then
   Aadvect(ii,ii)=Aadvect(ii,ii)+one
  else if (DErelaxation_model.eq.0) then ! e.g. incremental model
   if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental model
    ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
    ! do nothing
   else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental 
    ! Xia, Lu, Tryggvason 2018
    ! Seungwon Shin, Jalel Chergui, Damir Juric
    Aadvect(ii,ii)=Aadvect(ii,ii)+one
   else
    print *,"viscoelastic_model invalid: ",viscoelastic_model
    stop
   endif
  else
   print *,"DErelaxation_model invalid: ",DErelaxation_model
   stop
  endif
 enddo  ! ii=1,3

 unity_det=0
 call project_A_to_positive_definite_or_traceless(Aadvect, &
         viscoelastic_model,polymer_factor,unity_det)

  ! "Deborah number" relaxation 
  ! e.g. D^{triangle}/Dt A = -(1/De) (A-I)
 if (DErelaxation_model.eq.1) then
  ! do nothing
 else if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then !plastic model
   !Aadvect corresponds to "s" in Tran and Udaykumar
  if (DErelaxation_model.eq.0) then

    ! now Aadvect = s_{trial} (Tran and Udaykumar)
   do ii=1,3
   do jj=1,3
    Aadvect(ii,jj)=Aadvect(ii,jj)+dt*two*visctensorMAC_traceless(ii,jj) 
   enddo
   enddo

   magA=zero
   do ii=1,3
   do jj=1,3
    magA=magA+Aadvect(ii,jj)**2
   enddo
   enddo
   magA=sqrt(magA) !sqrt(A : A)
    !NP=A/sqrt(A:A)  (NP : NP=1)
   NP_dotdot_D=zero
   do ii=1,3
   do jj=1,3
    if (magA.gt.zero) then
     NP(ii,jj)=Aadvect(ii,jj)/magA
    else if (magA.eq.zero) then
     NP(ii,jj)=zero
    else
     print *,"magA invalid (point_updatetensor): ",magA
     stop
    endif
    NP_dotdot_D=NP_dotdot_D+NP(ii,jj)*visctensorMAC(ii,jj)
   enddo
   enddo


   do Johnson_iter=0,1

    if (Johnson_iter.eq.0) then
     plastic_strain_dot=zero
    else if (Johnson_iter.eq.1) then
     !do nothing
    else
     print *,"Johnson_iter invalid"
     stop
    endif

    ! "S" from Maire et al corresponds to "Aadvect" times the 
    ! shear modulus.
    ! see: Udaykumar, Tran, Belk, Vanden JCP 2003
    ! note: (13) from Tran and Udaykumar, have 
    ! sigma_y=(A+B(eps^p)^n)(1+C ln(epsdot^p/epsdot^0))(1-theta^m)
    call JohnsonCookSoftening( &
     xsten_in,nhalf_in, &
     i,j,k, &
     irefine,jrefine,krefine, &
     level, &
     finest_level, &
     im_critical, &  ! 0<=im_critical<=num_materials-1
     yield_stress, &
     cell_temperature, &
     fort_yield_temperature(im_critical+1), &
     fort_tempconst(im_critical+1), &
     fort_yield_alpha(im_critical+1), &
     fort_ref_plastic_strain(im_critical+1), &
     fort_ref_plastic_strain_dot(im_critical+1), &
     fort_Johnson_Cook_C(im_critical+1), &
     hardening_coefficient, &
     plastic_strain_old, &
     plastic_strain_dot, &
     fort_yield_n(im_critical+1), &
     gamma_not) !"yield_stress" intent(out)
  
    Y_plastic_parm_scaled=(gamma_not/elastic_viscosity)*sqrt(2.0d0/3.0d0)
     !magA=sqrt(A:A)
     !assume G=1 Tr(A)=0
     !=> f_plastic=sqrt(2/3)f   (see equation (2) from Ponthot)
    f_plastic=magA-Y_plastic_parm_scaled

    do ii=1,3
    do jj=1,3

     if ((f_plastic.lt.zero).or. &
         ((f_plastic.ge.zero).and.(NP_dotdot_D.le.zero))) then

      plastic_strain_dot=zero
      plastic_work=zero

      !NP=A/sqrt(A:A)  (NP : NP=1)
     else if ((f_plastic.ge.zero).and.(NP_dotdot_D.gt.zero)) then

       !Anew-Aold=sqrt(2/3) sigma_v A/sqrt(A:A) - A=
       !NP sqrt(A:A)(sqrt(2/3) sigma_v/sqrt(A:A)-1)=
       !NP(sqrt(2/3) sigma_v-NP sqrt(A:A))=-2 Gamma NP
       !Gamma=(1/2)(sqrt(A:A)-sqrt(2/3) sigma_v)
       !Note: Ponthot (21) has sqrt(sigma_v) by mistake.

       ! hardening_coefficient = h (Tran and Udaykumar)
      weight_prev=hardening_coefficient/(three*elastic_viscosity)

       ! Tran and Udaykumar:
       ! psi=(magA-sqrt(2/3) sigma_{v}^{0})/(2 G (1+h/(3G)))
      Aadvect(ii,jj)= &
        (weight_prev*magA+Y_plastic_parm_scaled)*NP(ii,jj)/ &
        (weight_prev+one)

        !https://www.brown.edu/Departments/Engineering/Courses/En1750/Notes/Plasticity/Plasticity.htm
        !search for "hardening"
        !Y(\bar{eps}^{p})=Y_{0} + h \bar{eps}^{p}
        !see "K" in the table under the wikipedia cite 
        !"Strain hardening exponent"
        !stainless steel=1275 MPa
        !copper=325 MPa
        !Table 3 Tran and Udaykumar:
        !B=177MPa n=0.12 C=0.016 m=1.0 Tungsten 
        !B=569MPa n=0.22 C=0.003 m=1.17 Steel
        !f_plastic=magA-sqrt(2/3) sigma_{v}^{0}
        !sigma_{v}^{1}=sigma_{v}^{0}+dt \sqrt{2/3} h\Gamma
      plastic_strain_dot=sqrt(f_plastic/(three*(one+weight_prev)))/dt

       ! (16) from Camacho and Ortiz.
       ! just below (11) in Tran and Udaykumar.
       ! !magA=sqrt(A:A)
      plastic_work=fort_mechanical_to_thermal(im_critical+1)* &
         plastic_strain_dot*magA*sqrt(3.0d0/2.0d0)*elastic_viscosity

      if (plastic_strain_dot.ge.zero) then
       !do nothing
      else
       print *,"plastic_strain_dot out of range ",plastic_strain_dot
       stop
      endif
      if (plastic_work.ge.zero) then
       !do nothing
      else
       print *,"plastic_work out of range ",plastic_work
       stop
      endif

     else
      print *,"f_plastic or NP_dotdot_D invalid"
      print *,"f_plastic=",f_plastic
      print *,"NP_dotdot_D=",NP_dotdot_D
      stop
     endif

    enddo !jj=1,3
    enddo !ii=1,3

   enddo ! Johnson_iter=0,1

  else
   print *,"DErelaxation_model invalid: ",DErelaxation_model
   stop
  endif
 else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental
  if (DErelaxation_model.eq.0) then
   ! do nothing
  else
   print *,"DErelaxation_model invalid: ",DErelaxation_model
   stop
  endif
 else
  print *,"viscoelastic_model invalid: ",viscoelastic_model
  stop
 endif
 
 !if NN_MAIRE_ABGRALL:
 !Maire and Abgrall paper:
 !  W=(grad V - (grad V)^T)/2  (grad V)_{ij} = v_{j,i} = jacobian_{ji}
 !  W_ij=(v_{j,i}-v_{i,j})/2
 !Tran and Udaykumar:
 !  Omega=-W    Omega_ij=(v_i,j - v_j,i)/2
 !Maire and Abgrall
 ! S_t=2 G(D0-Dp)-SW+WS
 !Tran and Udaykumar:
 ! S_t=2 G(Dbar-Dp)-S Omega + OMEGA S
 !
 !then Smult_left=Smult_right=I-dt W=I+dt OMEGA
 !DQ/Dt=2(D0-D^P)+QW-WQ=
 !      2(D0-D^P)+OMEGA Q - Q OMEGA
 !   Q=S/mu Q=zero matrix at t=0 W=(grad U-grad U^T)/2
 !grad U_{ij} = U_{j,i}
 !W=(grad V - grad V^T)/2
 !W_{ij}=(V_{j,i}-V_{i,j})/2=-OMEGA (see Tran and Udaykumar)
 !Q^{n+1}=(I-dt WLEFT)Q^{n}(I-dt WRIGHT^T)=
 !        Q^{n}-dt WLEFT Q - dt Q WRIGHT^T + 
 !        dt^2 WLEFT Q WRIGHT^T \approx
 !        Q^{n}-dt WLEFT Q + dt Q WRIGHT
 do ii=1,3
 do jj=1,3
   !SA=S_left * A
  SA(ii,jj)=zero
  do kk=1,3
   SA(ii,jj)=SA(ii,jj)+Smult_left(ii,kk)*Aadvect(kk,jj)
  enddo
 enddo
 enddo
 
 do ii=1,3
 do jj=1,3
   !SAS=S_left * A * S_right^T
  SAS(ii,jj)=zero

  do kk=1,3
   SAS(ii,jj)=SAS(ii,jj)+SA(ii,kk)*Smult_right(jj,kk)
  enddo
  Q(ii,jj)=SAS(ii,jj)

 enddo  ! jj=1..3
 enddo  ! ii=1..3

  !NOTE: gradU already has hoop terms built in.
 if (SDIM.eq.3) then
  ! do nothing
 else if (SDIM.eq.2) then
  ! do nothing
 else
  print *,"dimension bust"
  stop
 endif

  ! Q=S A S^T at this stage
 unity_det=0
 call project_A_to_positive_definite_or_traceless(Q, &
    viscoelastic_model,polymer_factor,unity_det)

 if (force_unity_determinant.eq.1) then
  unity_det=1
  call project_A_to_positive_definite_or_traceless(Q, &
    viscoelastic_model,polymer_factor,unity_det)
 else if (force_unity_determinant.eq.0) then
  ! do nothing
 else
  print *,"force_unity_determinant invalid: ",force_unity_determinant
  stop
 endif

 do ii=1,3

  if (DErelaxation_model.eq.1) then
   Q(ii,ii)=Q(ii,ii)-one  ! Q <--  A-I
  else if (DErelaxation_model.eq.0) then ! e.g. incremental model
   if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental model
    ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
    ! do nothing
   else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental 
    ! Xia, Lu, Tryggvason (2018)
    ! Seungwon Shin, Jalel Chergui, Damir Juric
    Q(ii,ii)=Q(ii,ii)-one  ! Q <--  A-I
   else
    print *,"viscoelastic_model invalid: ",viscoelastic_model
    stop
   endif
  else
   print *,"DErelaxation_model invalid: ",DErelaxation_model
   stop
  endif

 enddo !ii=1,3

 ! note: for viscoelastic_model==
 !   NN_MAIRE_ABGRALL_ETAL,
 !   NN_NEO_HOOKEAN,
 !  modtime=lambda_tilde=elastic_time >> 1
 !
 ! lambda_tilde=f(A)/lambda=(1/(lambda(1-tr(A)/L^2)))
 ! ofs \equiv equilibrium_diagonal
 ! lambda_tilde (Q^n+1-Q^n)=-dt (Q^n+1 + ofs I)
 ! (Q^n+1-Q^n)/dt = -(Q^n+1+ofs)/lambda_tilde
 ! Q^n+1 * (1/dt+1/lambda_tilde) = Q_n/dt-ofs/lambda_tilde
 ! Q^n+1 * (1+dt/lambda_tilde) = Q_n-ofs*dt/lambda_tilde
 ! Q^{n+1}=(lambda_tilde/(lambda_tilde+dt))*(Q_n-ofs*dt/lambda_tilde)
 ! Q^n+1 = (Q_n/dt-ofs/lambda_tilde)*dt*lambda_tilde/
 !         (lambda_tilde+dt)
 ! Q^n+1=(Q^n * lambda_tilde - ofs * dt)/(lambda+dt) 
 ! Note for determinants:
 ! Anp1=alpha An + (1-alpha)I  alpha=1-dt
 ! new eigenvalues: alpha lambda + (1-alpha)
 ! in 2D the product is:
 ! (alpha l1 + 1-alpha)(alpha l2+1-alpha)= 
 ! (alpha l1 + 1-alpha)(alpha/l1+1-alpha)=
 ! alpha^2 + (1-alpha)(alpha)(l1+1/l1)+(1-alpha)^2
 do ii=1,3
 do jj=1,3
  if (ii.eq.jj) then
   Q(ii,jj)=(modtime*Q(ii,jj)-equilibrium_diagonal*dt)/(modtime+dt)
  else if (ii.ne.jj) then
   Q(ii,jj)=modtime*Q(ii,jj)/(modtime+dt)
  else
   print *,"ii or jj invalid"
   stop
  endif
 enddo !jj=1,3
 enddo !ii=1,3

 do ii=1,3
 do jj=1,3
  Aadvect(ii,jj)=Q(ii,jj)
 enddo
 enddo

 do ii=1,3
  if (DErelaxation_model.eq.1) then
   Aadvect(ii,ii)=Aadvect(ii,ii)+one
  else if (DErelaxation_model.eq.0) then ! incremental model
   if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental model
    ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
    ! do nothing
   else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental
    ! Xia, Lu, Tryggvason 2018
    ! Seungwon Shin, Jalel Chergui, Damir Juric
    Aadvect(ii,ii)=Aadvect(ii,ii)+one
   else
    print *,"viscoelastic_model invalid: ",viscoelastic_model
    stop
   endif
  else
   print *,"DErelaxation_model invalid: ",DErelaxation_model
   stop
  endif
 enddo ! do ii=1,3

 unity_det=0
 call project_A_to_positive_definite_or_traceless(Aadvect, &
   viscoelastic_model,polymer_factor,unity_det)

 do ii=1,3
 do jj=1,3
  Q(ii,jj)=Aadvect(ii,jj)
 enddo
 enddo

 do ii=1,3
  if (DErelaxation_model.eq.1) then
   Q(ii,ii)=Q(ii,ii)-one
  else if (DErelaxation_model.eq.0) then ! e.g. incremental model
   if (viscoelastic_model.eq.NN_MAIRE_ABGRALL_ETAL) then ! incremental model
    ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
    ! do nothing
   else if (viscoelastic_model.eq.NN_NEO_HOOKEAN) then ! incremental
    ! Xia, Lu, Tryggvason (2018)
    ! Seungwon Shin, Jalel Chergui, Damir Juric
    Q(ii,ii)=Q(ii,ii)-one  ! Q <--  A-I
   else
    print *,"viscoelastic_model invalid: ",viscoelastic_model
    stop
   endif
  else
   print *,"DErelaxation_model invalid: ",DErelaxation_model
   stop
  endif
 enddo

else 
 print *,"viscoelastic_model invalid: ",viscoelastic_model
 stop
endif

do dir_local=1,ENUM_NUM_TENSOR_TYPE_BASE
 call stress_index(dir_local,ii,jj)
 tnew(dir_local)=Q(ii,jj)
enddo
tnew(ENUM_NUM_TENSOR_TYPE_BASE+1)=plastic_strain_old+dt*plastic_strain_dot

return
end subroutine point_updatetensor




subroutine doit(problo,probhi,ncell,dx,tstop)
use probcommon_module
IMPLICIT NONE

real(amrex_real) problo,probhi,dx,tstop
integer ncell
real(amrex_real) SSold(-1:ncell)
real(amrex_real) SS(-1:ncell)
real(amrex_real) QQold(-1:ncell)
real(amrex_real) QQ(-1:ncell)
real(amrex_real) DD(-1:ncell)  ! bottom topography
real(amrex_real) xx(-1:ncell)
real(amrex_real) SSflux(0:ncell)
real(amrex_real) QQflux(0:ncell)
real(amrex_real) start_elevation
integer skip,nstep,i
real(amrex_real) local_gravity
real(amrex_real) time
real(amrex_real) elevation_right,elevation_left,u_right,u_left
real(amrex_real) maxu,maxc,den,mom,uu,cc,dt,lambda
real(amrex_real) denleft,denright,momleft,momright,pleft,pright
real(amrex_real) minz
integer icrit

integer igrid,jgrid  ! t index x index
real(amrex_real) delta_t_grid,delta_x_grid,t_grid,x_grid
real(amrex_real) thetax,thetat
integer hitgrid,last_index

delta_t_grid=tstop/SHALLOW_M
delta_x_grid=(probhi-problo)/SHALLOW_N

skip=2000

local_gravity=980.0d0
start_elevation=22.862d0

time=0.0
nstep=0

do i=-1,ncell
 xx(i)=problo+(i+0.5)*dx
 call get_bottom_elevation(xx(i),DD(i))
 QQ(i)=0.0
 SS(i)=start_elevation-DD(i)
enddo

do while (time.le.tstop-EPS10) 

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

real(amrex_real) problo,probhi,tstop,dx
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

real(amrex_real),INTENT(in)  :: center(SDIM),x_in(SDIM)
real(amrex_real),INTENT(in)  :: r
real(amrex_real),INTENT(in)  :: Tinf,Tsat
real(amrex_real)             :: phi
real(amrex_real),INTENT(in)  :: delta   ! size of smooth transition region
real(amrex_real),INTENT(out) :: Tout
integer          :: i
real(amrex_real)  :: H_local
real(amrex_real)  :: phi_shift
real(amrex_real)  :: half_delta

phi=zero
do i=1,SDIM
 phi=phi+(x_in(i)-center(i))**2.0d0
enddo
phi=sqrt(phi)-r

! phi=0,  T=Tsat
! phi=delta,  T=Tinf

half_delta=half*delta

phi_shift=phi-half_delta

H_local=hs(phi_shift,half_delta)

Tout=Tinf*H_local+Tsat*(one-H_local)

end subroutine smooth_init

subroutine get_primary_material(LS,im_primary)
use probcommon_module

IMPLICIT NONE

real(amrex_real), INTENT(in) :: LS(num_materials)
integer, INTENT(out) :: im_primary
integer im,imtest
integer, parameter :: tessellate=0
integer is_rigid_local(num_materials)

do im=1,num_materials
 is_rigid_local(im)=is_rigid(im)
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
  print *,"tessellate invalid38: ",tessellate
  stop
 endif
enddo ! im=1..num_materials

im_primary=0
do im=1,num_materials
 if (is_rigid_local(im).eq.1) then
  if (LS(im).ge.zero) then
   if (im_primary.ne.0) then
    print *,"cannot have two rigid materials in same place"
    do imtest=1,num_materials
     print *,"imtest,LS(imtest) ",imtest,LS(imtest)
    enddo
    stop
   endif
   im_primary=im
  else if (LS(im).le.zero) then
   ! do nothing
  else
   print *,"LS bust in get_primary_material: ",im,LS(im)
   stop
  endif
 else if (is_rigid_local(im).eq.0) then
  ! do nothing
 else
  print *,"is_rigid_local invalid GLOBALUTIL.F90: ",im,is_rigid_local(im)
  stop
 endif
enddo !im=1..num_materials

if (im_primary.eq.0) then

 do im=1,num_materials
   if (im_primary.eq.0) then
    im_primary=im
   else if ((im_primary.ge.1).and.(im_primary.lt.im)) then
    if (LS(im).gt.LS(im_primary)) then
     im_primary=im
    else if (LS(im).le.LS(im_primary)) then
     ! do nothing
    else
     print *,"LS bust in get_primary_material(2): ", &
             im,im_primary,LS(im),LS(im_primary)
     stop
    endif
   else
    print *,"im_primary invalid: ",im_primary
    stop
   endif
 enddo !im=1..num_materials

else if ((im_primary.ge.1).and. &
         (im_primary.le.num_materials).and. &
         (is_rigid_local(im_primary).eq.1)) then
 ! do nothing
else
 print *,"is_rigid or im_primary invalid"
 print *,"im_primary=",im_primary
 print *,"is_rigid_local(im_primary)=",is_rigid_local(im_primary)
 stop
endif

end subroutine get_primary_material


integer function project_option_is_validF(project_option) 
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: project_option

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

integer function project_option_momeqnF(project_option) 
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: project_option

 if ((project_option.eq.SOLVETYPE_PRES).or. & ! regular project
     (project_option.eq.SOLVETYPE_INITPROJ).or. & ! initial project
     (project_option.eq.SOLVETYPE_PRESEXTRAP).or.& ! pressure extrapolation
     (project_option.eq.SOLVETYPE_VISC)) then      ! viscosity
  project_option_momeqnF=1
 else if ((project_option.eq.SOLVETYPE_HEAT).or. & ! thermal diffusion
          ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
           (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then
  project_option_momeqnF=0
 else
  print *,"project_option invalid"
  stop
  project_option_momeqnF=0
 endif

end function project_option_momeqnF

integer function project_option_singular_possibleF(project_option) &
bind(c,name='project_option_singular_possibleF')
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: project_option

 if ((project_option.eq.SOLVETYPE_PRES).or. & ! regular project
     (project_option.eq.SOLVETYPE_INITPROJ).or. & ! initial project
     (project_option.eq.SOLVETYPE_PRESEXTRAP)) then ! pressure extension
  project_option_singular_possibleF=1
 else if ((project_option.eq.SOLVETYPE_HEAT).or. & ! thermal diffusion
          (project_option.eq.SOLVETYPE_VISC).or. & ! viscosity
          ((project_option.ge.SOLVETYPE_SPEC).and. &
           (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then !species
  project_option_singular_possibleF=0
 else
  print *,"project_option invalid"
  stop
  project_option_singular_possibleF=0
 endif

end function project_option_singular_possibleF

integer function project_option_olddata_neededF(project_option) 
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: project_option

 if ((project_option.eq.SOLVETYPE_PRES).or. & ! regular project
     (project_option.eq.SOLVETYPE_INITPROJ).or. & ! initial project
     (project_option.eq.SOLVETYPE_PRESEXTRAP)) then ! pressure extension
  project_option_olddata_neededF=0
 else if ((project_option.eq.SOLVETYPE_HEAT).or. & ! thermal diffusion
          (project_option.eq.SOLVETYPE_VISC).or. & ! viscosity
          ((project_option.ge.SOLVETYPE_SPEC).and. &
           (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then !species
  project_option_olddata_neededF=1
 else 
  print *,"project_option invalid"
  stop
  project_option_olddata_neededF=0
 endif

end function project_option_olddata_neededF

integer function project_option_pressureF(project_option)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: project_option

 if ((project_option.eq.SOLVETYPE_PRES).or. &
     (project_option.eq.SOLVETYPE_INITPROJ).or. &
     (project_option.eq.SOLVETYPE_PRESEXTRAP)) then  !pressure extrap
  project_option_pressureF=1
 else if ((project_option.eq.SOLVETYPE_HEAT).or. &  ! temperature
          (project_option.eq.SOLVETYPE_VISC).or. &  ! viscosity
          ((project_option.ge.SOLVETYPE_SPEC).and. &
           (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then!species
  project_option_pressureF=0
 else
  print *,"project_option invalid"
  stop
  project_option_pressureF=0
 endif 

end function project_option_pressureF


function project_option_needs_scalingF(project_option) &
bind(c,name='project_option_needs_scalingF')
use probcommon_module
IMPLICIT NONE

integer :: project_option_needs_scalingF
integer, INTENT(in) :: project_option

 if ((project_option.eq.SOLVETYPE_PRES).or. & 
     (project_option.eq.SOLVETYPE_PRESEXTRAP)) then 
  project_option_needs_scalingF=1
 else if ((project_option.eq.SOLVETYPE_INITPROJ).or. & 
          (project_option.eq.SOLVETYPE_HEAT).or. & 
          (project_option.eq.SOLVETYPE_VISC).or. &  
          ((project_option.ge.SOLVETYPE_SPEC).and. &
           (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then 
  project_option_needs_scalingF=0
 else
  print *,"project_option invalid"
  stop
  project_option_needs_scalingF=0
 endif 

end function project_option_needs_scalingF

! initialize blobdata if project_option_FSI_rigid==1
integer function project_option_FSI_rigid(project_option) &
bind(c,name='project_option_FSI_rigid')

use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: project_option

 if ((project_option.eq.SOLVETYPE_PRES).or. & ! regular project
     (project_option.eq.SOLVETYPE_INITPROJ)) then ! initial project
  project_option_FSI_rigid=1
 else if ((project_option.eq.SOLVETYPE_PRESEXTRAP).or. &
          (project_option.eq.SOLVETYPE_VISC).or. & ! viscosity
          (project_option.eq.SOLVETYPE_HEAT).or. & ! thermal diffusion
          ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
           (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then
  project_option_FSI_rigid=0
 else
  print *,"project_option invalid"
  stop
  project_option_FSI_rigid=0
 endif

end function project_option_FSI_rigid



function project_option_projectionF(project_option) &
bind(c,name='project_option_projectionF')
use probcommon_module
IMPLICIT NONE

integer :: project_option_projectionF
integer, INTENT(in) :: project_option

 if ((project_option.eq.SOLVETYPE_PRES).or. & 
     (project_option.eq.SOLVETYPE_INITPROJ)) then 
  project_option_projectionF=1
 else if ((project_option.eq.SOLVETYPE_PRESEXTRAP).or. & 
          (project_option.eq.SOLVETYPE_HEAT).or. &  ! temperature
          (project_option.eq.SOLVETYPE_VISC).or. &  ! viscosity
          ((project_option.ge.SOLVETYPE_SPEC).and. &
           (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then ! species
  project_option_projectionF=0
 else
  print *,"project_option invalid"
  stop
  project_option_projectionF=0
 endif 

end function project_option_projectionF

integer function is_GFM_freezing_modelF(freezing_model) 
IMPLICIT NONE

integer, INTENT(in) :: freezing_model

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

integer function is_hydrate_freezing_modelF(freezing_model) 
IMPLICIT NONE

integer, INTENT(in) :: freezing_model

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

integer function is_valid_freezing_modelF(freezing_model) 
IMPLICIT NONE

integer, INTENT(in) :: freezing_model

 if ((freezing_model.eq.5).or. & !Stefan model evaporation or condensation
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

! this has to do base 10 rounding in such a way 
! so that the ascii output matches the floating point
! output.
real(tecplot_real) function round_time(time)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: time
integer :: int_time
integer :: power
integer :: i
real(tecplot_real) :: local_time

round_time=time
if (time.lt.zero) then
 print *,"time cannot be negative in round_time ",time
 stop
else if (time.eq.zero) then
 round_time=time
else if (time.gt.zero) then
 power=0
 local_time=time
 do while (local_time.lt.one)
  local_time=local_time*ten
  power=power-1
 enddo
 do while (local_time.gt.one)
  local_time=local_time/ten
  power=power+1
 enddo
 do i=1,7
  local_time=local_time*ten
  power=power-1
 enddo 
 int_time=NINT(local_time)
 local_time=int_time
 do while (power.lt.0)
  local_time=local_time/ten
  power=power+1
 enddo
 do while (power.gt.0)
  local_time=local_time*ten
  power=power-1
 enddo
 round_time=local_time
else
 print *,"time is NaN"
 stop
endif

end function round_time

integer function is_multi_component_evapF(freezing_model, &
   evap_flag,latent_heat) 
IMPLICIT NONE

integer, INTENT(in) :: freezing_model
integer, INTENT(in) :: evap_flag
real(amrex_real), INTENT(in) :: latent_heat

 if (latent_heat.eq.zero) then
  is_multi_component_evapF=0
 else if (latent_heat.ne.zero) then

  if ((freezing_model.eq.5).or. & !Stefan model evaporation or condensation
      (freezing_model.eq.6).or. & !Palmore and Desjardins
      (freezing_model.eq.7)) then !cavitation

   if (evap_flag.eq.0) then !Palmore and Desjardins
    is_multi_component_evapF=1
   else if ((evap_flag.eq.1).or. & !Tanasawa
            (evap_flag.eq.2).or. & !Schrage
            (evap_flag.eq.3)) then !Kassemi
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

subroutine TopDownMergeSortReal(data_to_sort,A,B,n)
IMPLICIT NONE
integer, INTENT(in) :: n
real(amrex_real), allocatable, INTENT(in) :: data_to_sort(:)
integer, allocatable, INTENT(inout) :: A(:)
integer, allocatable, INTENT(inout) :: B(:)

 call CopyArrayReal(A,0,n,B)
 call TopDownSplitMergeReal(data_to_sort,B,0,n,A)

end subroutine TopDownMergeSortReal


recursive subroutine TopDownSplitMergeReal(data_to_sort,B,iBegin,iEnd,A)
IMPLICIT NONE
real(amrex_real), allocatable, INTENT(in) :: data_to_sort(:)
integer, INTENT(in) :: iBegin
integer, INTENT(in) :: iEnd
integer, allocatable, INTENT(inout) :: A(:)
integer, allocatable, INTENT(inout) :: B(:)
integer :: iMiddle

 if (iEnd-iBegin.le.1) then
  ! do nothing
 else
  iMiddle=(iEnd+iBegin)/2
  call TopDownSplitMergeReal(data_to_sort,A,iBegin,iMiddle,B)
  call TopDownSplitMergeReal(data_to_sort,A,iMiddle,iEnd,B)
  call TopDownMergeReal(data_to_sort,B,iBegin,iMiddle,iEnd,A)
 endif

end subroutine TopDownSplitMergeReal

subroutine TopDownMergeReal(data_to_sort,A,iBegin,iMiddle,iEnd,B)
IMPLICIT NONE
real(amrex_real), allocatable, INTENT(in) :: data_to_sort(:)
integer, INTENT(in) :: iBegin
integer, INTENT(in) :: iMiddle
integer, INTENT(in) :: iEnd
integer, allocatable, INTENT(inout) :: A(:)
integer, allocatable, INTENT(inout) :: B(:)
integer :: i,j,k,compare_flag

 i=iBegin
 j=iMiddle
 k=iBegin

 do while (k.lt.iEnd)
   !Ai<Aj compare_flag=-1
   !Ai>Aj compare_flag=1

  compare_flag=0

  if ((i.lt.iMiddle).and.(j.lt.iEnd)) then

   if (data_to_sort(A(i+1)).lt.data_to_sort(A(j+1))) then
    compare_flag=-1
   else if (data_to_sort(A(i+1)).gt.data_to_sort(A(j+1))) then
    compare_flag=1
   else if (data_to_sort(A(i+1)).eq.data_to_sort(A(j+1))) then
    compare_flag=0
   else
    print *,"data_to_sort NaN"
    stop
   endif

  else if ((i.ge.iMiddle).or.(j.ge.iEnd)) then
   ! do nothing
  else
   print *,"i,j bust"
   stop
  endif

  if ((i.lt.iMiddle).and. &
      ((j.ge.iEnd).or.(compare_flag.le.0))) then
   B(k+1)=A(i+1)
   i=i+1
  else
   B(k+1)=A(j+1)
   j=j+1
  endif
  k=k+1 
 enddo ! do while (k.lt.iEnd)

end subroutine TopDownMergeReal


subroutine CopyArrayReal(A,iBegin,iEnd,B)
IMPLICIT NONE
integer, INTENT(in) :: iBegin
integer, INTENT(in) :: iEnd
integer, allocatable, INTENT(inout) :: A(:)
integer, allocatable, INTENT(inout) :: B(:)
integer :: k

 do k=iBegin,iEnd-1
  B(k+1)=A(k+1)
 enddo

end subroutine CopyArrayReal

subroutine sort_branch_data(source_branch,splittingrule, &
     median_index,median_value)
use probcommon_module

Type(branch_type), INTENT(inout) :: source_branch
integer, INTENT(in) :: splittingrule
integer, INTENT(out) :: median_index
real(amrex_real), INTENT(out) :: median_value
integer, allocatable :: A_list(:)
integer, allocatable :: B_list(:)
real(amrex_real), allocatable :: data_to_sort(:)
real(amrex_real), allocatable :: save_data_decisions(:,:)
real(amrex_real), allocatable :: save_data_classify(:,:)
integer :: datalo(2)
integer :: datahi(2)
integer :: datalo_classify(2)
integer :: datahi_classify(2)
integer :: idata,dir
real(amrex_real)    :: data1,data2

 datalo=LBOUND(source_branch%data_decisions)
 datahi=UBOUND(source_branch%data_decisions)

 allocate(save_data_decisions(datahi(1),datahi(2)))
 save_data_decisions=source_branch%data_decisions

 datalo_classify=LBOUND(source_branch%data_classify)
 datahi_classify=UBOUND(source_branch%data_classify)

 allocate(save_data_classify(datahi_classify(1),datahi_classify(2)))
 save_data_classify=source_branch%data_classify

 if ((datalo(1).eq.1).and. &
     (datalo(2).eq.1)) then
  if ((datahi(1).eq.source_branch%ndata).and. &
      (datahi(1).ge.2).and.  &
      (datahi(2).ge.splittingrule)) then
   allocate(A_list(datahi(1)))
   allocate(B_list(datahi(1)))
   allocate(data_to_sort(datahi(1)))
   do idata=1,datahi(1)
    A_list(idata)=idata
    B_list(idata)=idata
    data_to_sort(idata)=source_branch%data_decisions(idata,splittingrule)
   enddo
   call TopDownMergeSortReal(data_to_sort,A_list,B_list,datahi(1))

    !sanity check
   do idata=1,datahi(1)-1
    if (data_to_sort(A_list(idata)).le.data_to_sort(A_list(idata+1))) then
     ! do nothing
    else
     print *,"data_to_sort not sorted properly"
     stop
    endif
   enddo

   do idata=1,datahi(1)
    do dir=1,datahi(2)
     source_branch%data_decisions(idata,dir)= &
            save_data_decisions(A_list(idata),dir)
    enddo
   enddo

   if ((datalo_classify(1).eq.1).and. &
       (datalo_classify(2).eq.1)) then
    if ((datahi_classify(1).eq.source_branch%ndata).and. &
        (datahi_classify(1).ge.2)) then

     do idata=1,datahi_classify(1)
      do dir=1,datahi_classify(2)
       source_branch%data_classify(idata,dir)= &
            save_data_classify(A_list(idata),dir)
      enddo
     enddo

    else
     print *,"datahi invalid"
     stop
    endif
   else
    print *,"datalo invalid"
    stop
   endif

   !sanity check
   do idata=1,datahi(1)-1
    data1=source_branch%data_decisions(idata,splittingrule)
    data2=source_branch%data_decisions(idata+1,splittingrule)
    if (data1.le.data2) then
     ! do nothing
    else
     print *,"source_branch not sorted properly"
     stop
    endif
   enddo

   median_index=datahi(1)/2
   source_branch%median_index=median_index

   if ((median_index.ge.1).and. &
       (median_index.lt.datahi(1))) then
    median_value= &
     half*(source_branch%data_decisions(median_index,splittingrule)+ &
           source_branch%data_decisions(median_index+1,splittingrule))
    source_branch%median_value=median_value
   else
    print *,"median_index invalid"
    stop
   endif

   deallocate(A_list)
   deallocate(B_list)
   deallocate(data_to_sort)

   deallocate(save_data_decisions)
   deallocate(save_data_classify)

  else
   print *,"datahi invalid"
   stop
  endif
 else
  print *,"datalo invalid"
  stop
 endif

end subroutine sort_branch_data

subroutine statistics_branch_data( &
   ndim_classify, &
   source_branch, &
   splittingrule, &
   variance_reduction)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: ndim_classify
integer, INTENT(in) :: splittingrule
Type(branch_type), INTENT(in) :: source_branch
real(amrex_real), INTENT(out) :: variance_reduction
integer :: datalo(2)
integer :: datahi(2)
integer :: datalo_classify(2)
integer :: datahi_classify(2)
integer :: idata,dir
real(amrex_real) :: mean(ndim_classify)
real(amrex_real) :: mean_child1(ndim_classify)
real(amrex_real) :: mean_child2(ndim_classify)
real(amrex_real) :: variance
real(amrex_real) :: variance_scale
real(amrex_real) :: variance_child1
real(amrex_real) :: variance_child2

 datalo=LBOUND(source_branch%data_decisions)
 datahi=UBOUND(source_branch%data_decisions)

 datalo_classify=LBOUND(source_branch%data_classify)
 datahi_classify=UBOUND(source_branch%data_classify)

 if ((datalo(1).eq.1).and. &
     (datalo(2).eq.1).and. &
     (datalo_classify(1).eq.1).and. &
     (datalo_classify(2).eq.1)) then
  if ((datahi(1).eq.source_branch%ndata).and. &
      (datahi(1).ge.2).and. &
      (datahi(2).ge.splittingrule).and. &
      (datahi_classify(1).eq.source_branch%ndata).and. &
      (datahi_classify(2).eq.ndim_classify)) then

   do dir=1,ndim_classify
    mean(dir)=zero
    mean_child1(dir)=zero
    mean_child2(dir)=zero
   enddo

   do idata=1,datahi_classify(1)
    do dir=1,ndim_classify
     mean(dir)=mean(dir)+source_branch%data_classify(idata,dir)
     if (idata.le.source_branch%median_index) then
      mean_child1(dir)=mean_child1(dir)+ &
          source_branch%data_classify(idata,dir)
     else if ((idata.gt.source_branch%median_index).and. &
              (idata.le.source_branch%ndata)) then
      mean_child2(dir)=mean_child2(dir)+ &
          source_branch%data_classify(idata,dir)
     else
      print *,"idata invalid"
      stop
     endif
    enddo
   enddo
   do dir=1,ndim_classify
    mean(dir)=mean(dir)/source_branch%ndata

    if ((source_branch%median_index.ge.1).and. &
        (source_branch%median_index.lt.source_branch%ndata)) then
     mean_child1(dir)=mean_child1(dir)/source_branch%median_index
     mean_child2(dir)=mean_child2(dir)/ &
        (source_branch%ndata-source_branch%median_index)
    else
     print *,"median_index invalid"
     stop
    endif
   enddo
   variance=zero
   variance_child1=zero
   variance_child2=zero
    
   do idata=1,datahi_classify(1)
    do dir=1,ndim_classify
     variance=variance+ &
        (source_branch%data_classify(idata,dir)-mean(dir))**2
     if (idata.le.source_branch%median_index) then
      variance_child1=variance_child1+ &
        (source_branch%data_classify(idata,dir)-mean_child1(dir))**2
     else if ((idata.gt.source_branch%median_index).and. &
              (idata.le.source_branch%ndata)) then
      variance_child2=variance_child2+ &
        (source_branch%data_classify(idata,dir)-mean_child2(dir))**2
     else
      print *,"idata invalid"
      stop
     endif
    enddo
   enddo

  else
   print *,"datahi invalid"
   stop
  endif
 else
  print *,"datalo invalid"
  stop
 endif

 variance_reduction=variance-(variance_child1+variance_child2)
 variance_scale=max(variance,one)

 if (variance_reduction.ge.zero) then
  !do nothing
 else if (variance_reduction.ge.-EPS_8_4*variance_scale) then
  variance_reduction=zero
 else
  print *,"variance_reduction invalid: ",variance_reduction
  print *,"variance: ",variance
  print *,"variance_scale: ",variance_scale
  print *,"variance_child1: ",variance_child1
  print *,"variance_child2: ",variance_child2
  stop
 endif

end subroutine statistics_branch_data


subroutine children_branch_data( &
   tree_var, &
   ndim_decisions, &
   ndim_classify, &
   current_level, &
   current_branch_id, &
   source_branch, &
   splittingrule)
use probcommon_module
IMPLICIT NONE

Type(tree_type), INTENT(inout) :: tree_var
integer, INTENT(in) :: ndim_decisions
integer, INTENT(in) :: ndim_classify
integer, INTENT(in) :: splittingrule
Type(branch_type), INTENT(inout) :: source_branch
integer, INTENT(in) :: current_level
integer, INTENT(inout) :: current_branch_id
integer :: datalo(2)
integer :: datahi(2)
integer :: datalo_classify(2)
integer :: datahi_classify(2)
integer :: idata,dir
real(amrex_real), allocatable :: save_data_decisions(:,:)
real(amrex_real), allocatable :: save_data_classify(:,:)

integer :: local_ndata
integer :: local_parent_id
integer :: local_parent_level
integer :: local_current_id
integer :: local_current_level
integer :: local_splittingrule
integer :: local_median_index
real(amrex_real) :: local_median_value
integer :: local_child1_id
integer :: local_child_level
integer :: local_child2_id

 datalo=LBOUND(source_branch%data_decisions)
 datahi=UBOUND(source_branch%data_decisions)

 datalo_classify=LBOUND(source_branch%data_classify)
 datahi_classify=UBOUND(source_branch%data_classify)

 if ((datalo(1).eq.1).and. &
     (datalo(2).eq.1).and. &
     (datalo_classify(1).eq.1).and. &
     (datalo_classify(2).eq.1)) then
  if ((datahi(1).eq.source_branch%ndata).and. &
      (datahi(1).ge.2).and. &
      (datahi(2).ge.splittingrule).and. &
      (datahi(2).eq.ndim_decisions).and. &
      (datahi_classify(1).eq.source_branch%ndata).and. &
      (datahi_classify(2).eq.ndim_classify)) then

   source_branch%splittingrule=splittingrule

   current_branch_id=current_branch_id+2
   tree_var%branch_list_level(current_level+1)%nbranches=current_branch_id

   source_branch%child1_id=current_branch_id-1
   source_branch%child2_id=current_branch_id
   source_branch%child_level=current_level+1

   if (source_branch%current_level.eq.current_level) then
    ! do nothing
   else
    print *,"current_level invalid(1)"
    print *,"current_level=",current_level
    print *,"source_branch%current_level=",source_branch%current_level
    stop
   endif

   local_ndata=source_branch%median_index
   local_parent_id=source_branch%current_id
   local_parent_level=current_level
   local_current_id=current_branch_id-1
   local_current_level=current_level+1
   local_splittingrule=-1
   local_median_index=-1
   local_median_value=zero
   local_child1_id=-1
   local_child2_id=-1
   local_child_level=-1

   allocate(save_data_decisions(local_ndata,ndim_decisions))
   allocate(save_data_classify(local_ndata,ndim_classify))

   do idata=1,local_ndata
    do dir=1,ndim_decisions
     save_data_decisions(idata,dir)=source_branch%data_decisions(idata,dir)
    enddo
   enddo
   do idata=1,local_ndata
    do dir=1,ndim_classify
     save_data_classify(idata,dir)=source_branch%data_classify(idata,dir)
    enddo
   enddo

   call init_branch( &
    tree_var%branch_list_level(local_current_level)% &
       branch_list(local_current_id), &
    ndim_decisions, &
    ndim_classify, &
    save_data_decisions, &
    save_data_classify, &
    local_ndata, &
    local_parent_id, &
    local_parent_level, &
    local_current_id, &
    local_current_level, &
    local_splittingrule, &
    local_median_index, &
    local_median_value, &
    local_child1_id, &
    local_child2_id, &
    local_child_level)

   deallocate(save_data_decisions)
   deallocate(save_data_classify)

   local_ndata=source_branch%ndata-source_branch%median_index
   local_parent_id=source_branch%current_id
   local_parent_level=current_level
   local_current_id=current_branch_id
   local_current_level=current_level+1
   local_splittingrule=-1
   local_median_index=-1
   local_median_value=zero
   local_child1_id=-1
   local_child2_id=-1
   local_child_level=-1

   allocate(save_data_decisions(local_ndata,ndim_decisions))
   allocate(save_data_classify(local_ndata,ndim_classify))

   do idata=1,local_ndata
    do dir=1,ndim_decisions
     save_data_decisions(idata,dir)= &
          source_branch%data_decisions(idata+source_branch%median_index,dir)
    enddo
   enddo
   do idata=1,local_ndata
    do dir=1,ndim_classify
     save_data_classify(idata,dir)= &
          source_branch%data_classify(idata+source_branch%median_index,dir)
    enddo
   enddo

   call init_branch( &
    tree_var%branch_list_level(local_current_level)% &
       branch_list(local_current_id), &
    ndim_decisions, &
    ndim_classify, &
    save_data_decisions, &
    save_data_classify, &
    local_ndata, &
    local_parent_id, &
    local_parent_level, &
    local_current_id, &
    local_current_level, &
    local_splittingrule, &
    local_median_index, &
    local_median_value, &
    local_child1_id, &
    local_child2_id, &
    local_child_level)

   deallocate(save_data_decisions)
   deallocate(save_data_classify)

   deallocate(source_branch%data_decisions)
   deallocate(source_branch%data_classify)

  else
   print *,"datahi invalid"
   stop
  endif
 else
  print *,"datalo invalid"
  stop
 endif

end subroutine children_branch_data

subroutine test_variance_results( &
   nbranches, &
   ndim_classify, &
   splittingrule, &
   tree_var, &
   current_level, &
   variance_reduction)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nbranches
integer, INTENT(in) :: ndim_classify
integer, INTENT(in) :: splittingrule
Type(tree_type), INTENT(inout) :: tree_var
integer, INTENT(in) :: current_level
real(amrex_real), INTENT(out) :: variance_reduction(nbranches)
integer :: ibranch
integer :: local_ndata
integer :: local_median_index
real(amrex_real) :: local_median_value
real(amrex_real) :: local_variance_reduction

 if ((nbranches.eq.tree_var%nbranches_level(current_level)).and. &
     (nbranches.eq.tree_var%branch_list_level(current_level)%nbranches)) then
  ! do nothing
 else
  print *,"nbranches corrupt(1)"
  print *,"nbranches=",nbranches
  print *,"tree_var%nbranches_level(current_level)= ", &
    tree_var%nbranches_level(current_level)
  stop
 endif

 do ibranch=1,nbranches
  variance_reduction(ibranch)=zero
 enddo

 do ibranch=1,nbranches
  local_ndata=tree_var%branch_list_level(current_level)% &
          branch_list(ibranch)%ndata
  if (local_ndata.eq.1) then
   ! do nothing
  else if (local_ndata.ge.2) then
   call sort_branch_data( &
    tree_var%branch_list_level(current_level)%branch_list(ibranch), &
    splittingrule, & 
    local_median_index, &
    local_median_value)
   call statistics_branch_data( &
    ndim_classify, &
    tree_var%branch_list_level(current_level)%branch_list(ibranch), &
    splittingrule, & 
    local_variance_reduction)
   variance_reduction(ibranch)=variance_reduction(ibranch)+ &
       local_variance_reduction
  else
   print *,"local_ndata invalid"
   stop
  endif
 enddo ! do ibranch=1,nbranches

end subroutine test_variance_results

subroutine init_new_tree_level( &
   nbranches, &
   ndim_decisions, &
   ndim_classify, &
   splittingrule, &
   tree_var, &
   current_level)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nbranches
integer, INTENT(in) :: ndim_decisions
integer, INTENT(in) :: ndim_classify
integer, INTENT(in) :: splittingrule(nbranches)
Type(tree_type), INTENT(inout) :: tree_var
integer, INTENT(in) :: current_level
integer :: current_branch_id
integer :: ibranch
integer :: local_ndata
integer :: local_median_index
real(amrex_real) :: local_median_value

 if ((nbranches.eq.tree_var%nbranches_level(current_level)).and. &
     (nbranches.eq.tree_var%branch_list_level(current_level)%nbranches)) then
  ! do nothing
 else
  print *,"nbranches corrupt(2)"
  print *,"nbranches=",nbranches
  print *,"tree_var%nbranches_level(current_level)= ", &
    tree_var%nbranches_level(current_level)
  stop
 endif

 current_branch_id=0
 do ibranch=1,nbranches
  local_ndata=tree_var%branch_list_level(current_level)% &
          branch_list(ibranch)%ndata
  if (local_ndata.eq.1) then
   ! do nothing
  else if (local_ndata.ge.2) then
   call sort_branch_data( &
    tree_var%branch_list_level(current_level)%branch_list(ibranch), &
    splittingrule(ibranch), & 
    local_median_index, &
    local_median_value)
   call children_branch_data( &
    tree_var, &
    ndim_decisions, &
    ndim_classify, &
    current_level, &
    current_branch_id, &
    tree_var%branch_list_level(current_level)%branch_list(ibranch), &
    splittingrule(ibranch))
  else
   print *,"local_ndata invalid"
   stop
  endif
 enddo ! ibranch=1,nbranches

end subroutine init_new_tree_level

subroutine init_branch( &
   dest_branch, &
   ndim_decisions, &
   ndim_classify, &
   data_decisions, &
   data_classify, &
   ndata, &
   parent_id, &
   parent_level, &
   current_id, &
   current_level, &
   splittingrule, &
   median_index, &
   median_value, &
   child1_id, &
   child2_id, &
   child_level)
use probcommon_module
IMPLICIT NONE

Type(branch_type), INTENT(out) :: dest_branch
integer, INTENT(in) :: ndim_decisions
integer, INTENT(in) :: ndim_classify
integer, INTENT(in) :: ndata
real(amrex_real), INTENT(in) :: data_decisions(ndata,ndim_decisions)
real(amrex_real), INTENT(in) :: data_classify(ndata,ndim_classify)
integer, INTENT(in) :: parent_id
integer, INTENT(in) :: parent_level
integer, INTENT(in) :: current_id
integer, INTENT(in) :: current_level
integer, INTENT(in) :: splittingrule
integer, INTENT(in) :: median_index
real(amrex_real), INTENT(in) :: median_value
integer, INTENT(in) :: child1_id
integer, INTENT(in) :: child2_id
integer, INTENT(in) :: child_level

 if (ndata.ge.1) then
  ! do nothing
 else
  print *,"ndata invalid"
  stop
 endif

 dest_branch%ndata=ndata
 dest_branch%parent_id=parent_id
 dest_branch%parent_level=parent_level
 dest_branch%current_id=current_id
 dest_branch%current_level=current_level
 dest_branch%splittingrule=splittingrule
 dest_branch%median_index=median_index
 dest_branch%median_value=median_value
 dest_branch%child1_id=child1_id
 dest_branch%child2_id=child2_id
 dest_branch%child_level=child_level

 allocate(dest_branch%data_decisions(ndata,ndim_decisions))
 dest_branch%data_decisions=data_decisions
 allocate(dest_branch%data_classify(ndata,ndim_classify))
 dest_branch%data_classify=data_classify
    
end subroutine init_branch

subroutine initialize_decision_tree(data_decisions,data_classify, &
        nsamples,ndim_decisions,ndim_classify,tree_var)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nsamples
integer, INTENT(in) :: ndim_decisions
integer, INTENT(in) :: ndim_classify
real(amrex_real), INTENT(in) :: data_decisions(nsamples,ndim_decisions)
real(amrex_real), INTENT(in) :: data_classify(nsamples,ndim_classify)
Type(tree_type), INTENT(out) :: tree_var

integer :: nsamples_copy
integer :: datahi_decisions(2)
integer :: datahi_classify(2)

integer :: local_ndata
integer :: local_parent_id
integer :: local_parent_level
integer :: local_current_id
integer :: local_current_level
integer :: local_splittingrule
integer :: local_median_index
real(amrex_real) :: local_median_value
integer :: local_child1_id
integer :: local_child_level
integer :: local_child2_id

real(amrex_real), allocatable :: local_variance_reduction(:)
real(amrex_real), allocatable :: max_variance_reduction(:)
integer, allocatable :: max_splittingrule(:)
real(amrex_real) :: total_variance_reduction

integer :: local_nbranches
integer :: local_nbranches_next_level
integer :: ibranch

 tree_var%max_number_tree_levels=1
 tree_var%number_tree_levels=1
 
 nsamples_copy=nsamples
 do while (nsamples_copy.gt.0) 
  nsamples_copy=nsamples_copy/2
  tree_var%max_number_tree_levels=tree_var%max_number_tree_levels+1
 enddo

 print *,"nsamples,max_number_tree_levels ",nsamples, &
         tree_var%max_number_tree_levels
 print *,"ndim_decisions,ndim_classify ",ndim_decisions,ndim_classify

 allocate(tree_var%nbranches_level(tree_var%max_number_tree_levels))
 allocate(tree_var%branch_list_level(tree_var%max_number_tree_levels))

 datahi_decisions=UBOUND(data_decisions)
 datahi_classify=UBOUND(data_classify)
 if ((datahi_decisions(2).eq.ndim_decisions).and. &
     (datahi_decisions(1).eq.nsamples)) then
  ! do nothing
 else
  print *,"datahi_decisions failed"
  stop
 endif
 if ((datahi_classify(2).eq.ndim_classify).and. &
     (datahi_classify(1).eq.nsamples)) then
  ! do nothing
 else
  print *,"datahi_classify failed"
  stop
 endif

 tree_var%nbranches_level(1)=1
 tree_var%branch_list_level(1)%nbranches=1
 allocate(tree_var%branch_list_level(1)%branch_list(1))

 local_ndata=nsamples
 local_parent_id=-1
 local_parent_level=-1
 local_current_id=1
 local_current_level=1
 local_splittingrule=-1
 local_median_index=-1
 local_median_value=zero
 local_child1_id=-1
 local_child2_id=-1
 local_child_level=2

 call init_branch( &
  tree_var%branch_list_level(local_current_level)%branch_list(1), &
  ndim_decisions, &
  ndim_classify, &
  data_decisions, &
  data_classify, &
  local_ndata, &
  local_parent_id, &
  local_parent_level, &
  local_current_id, &
  local_current_level, &
  local_splittingrule, &
  local_median_index, &
  local_median_value, &
  local_child1_id, &
  local_child2_id, &
  local_child_level)

 do while (local_current_level.lt.tree_var%max_number_tree_levels)

   local_nbranches=tree_var%branch_list_level(local_current_level)%nbranches
   allocate(tree_var%branch_list_level(local_current_level+1)% &
        branch_list(2*local_nbranches))

   local_nbranches_next_level=0

   do ibranch=1,local_nbranches
    local_ndata=tree_var%branch_list_level(local_current_level)% &
            branch_list(ibranch)%ndata
    if (local_ndata.ge.2) then
     local_nbranches_next_level=local_nbranches_next_level+2
    else if (local_ndata.eq.1) then
     ! do nothing
    else
     print *,"local_ndata invalid"
     stop
    endif
   enddo !ibranch=1,local_nbranches

   tree_var%nbranches_level(local_current_level+1)= &
       local_nbranches_next_level
   tree_var%branch_list_level(local_current_level+1)%nbranches= &
       local_nbranches_next_level

   if ((local_nbranches_next_level.ge.0).and. &
       (local_nbranches_next_level.le.2*local_nbranches)) then

    if (local_nbranches_next_level.gt.0) then
     tree_var%number_tree_levels=local_current_level+1
    endif

   else
    print *,"local_nbranches_next_level invalid"
    stop
   endif

   allocate(max_variance_reduction(local_nbranches))
   allocate(local_variance_reduction(local_nbranches))
   allocate(max_splittingrule(local_nbranches))

   do ibranch=1,local_nbranches
    max_variance_reduction(ibranch)=zero
    local_variance_reduction(ibranch)=zero
    max_splittingrule(ibranch)=0
   enddo
   total_variance_reduction=zero

   do local_splittingrule=1,ndim_decisions

    call test_variance_results( &
     local_nbranches, &
     ndim_classify, &
     local_splittingrule, &
     tree_var, &
     local_current_level, &
     local_variance_reduction)

    total_variance_reduction=zero
    do ibranch=1,local_nbranches
     if (local_variance_reduction(ibranch).ge.zero) then
      total_variance_reduction=total_variance_reduction+ &
       local_variance_reduction(ibranch)
     else
      print *,"local_variance_reduction is NaN"
      print *,"ibranch,local_nbranches: ",ibranch,local_nbranches
      print *,"local_variance_reduction(ibranch): ", &
         local_variance_reduction(ibranch)
      stop
     endif
    enddo

    if (total_variance_reduction.ge.zero) then
     ! do nothing
    else
     print *,"total_variance_reduction is NaN: ",total_variance_reduction
     stop
    endif

    if (1.eq.1) then
     do ibranch=1,local_nbranches
      local_variance_reduction(ibranch)=total_variance_reduction
     enddo
    endif

    if (local_splittingrule.eq.1) then
     do ibranch=1,local_nbranches
      max_variance_reduction(ibranch)=local_variance_reduction(ibranch)
      max_splittingrule(ibranch)=local_splittingrule
     enddo
    else if ((local_splittingrule.gt.1).and. &
             (local_splittingrule.le.ndim_decisions)) then

     do ibranch=1,local_nbranches

      if (local_variance_reduction(ibranch).gt. &
          max_variance_reduction(ibranch)) then
       max_variance_reduction(ibranch)=local_variance_reduction(ibranch)
       max_splittingrule(ibranch)=local_splittingrule
      else if (local_variance_reduction(ibranch).le. &
               max_variance_reduction(ibranch)) then
       ! do nothing
      else
       print *,"NaN encountered"
       stop
      endif

     enddo ! do ibranch=1,local_nbranches

    else
     print *,"local_splittingrule invalid"
     stop
    endif

   enddo !do local_splittingrule=1,ndim_decisions

   call init_new_tree_level( &
     local_nbranches, &
     ndim_decisions, &
     ndim_classify, &
     max_splittingrule, &
     tree_var, &
     local_current_level)

   deallocate(max_variance_reduction)
   deallocate(local_variance_reduction)
   deallocate(max_splittingrule)

   local_current_level=local_current_level+1
 enddo !while (local_current_level.lt.tree_var%max_number_tree_levels)

end subroutine initialize_decision_tree

subroutine decision_tree_predict(data_decision,data_classified, &
  ndim_decisions,ndim_classify,tree_var)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: ndim_decisions
integer, INTENT(in) :: ndim_classify
real(amrex_real), INTENT(in) :: data_decision(ndim_decisions)
real(amrex_real), INTENT(out) :: data_classified(ndim_classify)
Type(tree_type), INTENT(in) :: tree_var
integer :: prev_id
integer :: prev_level
integer :: current_level
integer :: current_id
integer :: current_ndata
integer :: previous_ndata
integer :: splittingrule
integer :: median_index
integer :: dir
real(amrex_real) :: median_value
integer :: datalo(2)
integer :: datahi(2)

 if ((tree_var%number_tree_levels.ge.1).and. &
     (tree_var%number_tree_levels.le. &
      tree_var%max_number_tree_levels)) then
  ! do nothing
 else
  print *,"number_tree_levels invalid"
  stop
 endif

 current_level=1
 current_id=1
 current_ndata=tree_var%branch_list_level(current_level)% &
    branch_list(current_id)%ndata

 do while (current_ndata.ge.2)

  if (current_level.le.tree_var%number_tree_levels) then
   ! do nothing
  else
   print *,"current_level invalid(2)"
   print *,"current_level=",current_level
   print *,"tree_var%number_tree_levels=", &
     tree_var%number_tree_levels
   stop
  endif

  if (current_id.le.tree_var%branch_list_level(current_level)%nbranches) then
   ! do nothing
  else
   print *,"current_id invalid"
   stop
  endif

  datalo=LBOUND(tree_var%branch_list_level(current_level)% &
          branch_list(current_id)%data_decisions)
  datahi=UBOUND(tree_var%branch_list_level(current_level)% &
          branch_list(current_id)%data_decisions)
  if ((datalo(1).eq.1).and. &
      (datalo(2).eq.1).and. &
      (datahi(1).eq.current_ndata).and. &
      (datahi(2).eq.ndim_decisions)) then
   ! do nothing
  else
   print *,"datalo or datahi invalid"
   stop
  endif

  datalo=LBOUND(tree_var%branch_list_level(current_level)% &
          branch_list(current_id)%data_classify)
  datahi=UBOUND(tree_var%branch_list_level(current_level)% &
          branch_list(current_id)%data_classify)

  if ((datalo(1).eq.1).and. &
      (datalo(2).eq.1).and. &
      (datahi(1).eq.current_ndata).and. &
      (datahi(2).eq.ndim_classify)) then
   ! do nothing
  else
   print *,"datalo or datahi invalid"
   stop
  endif

  splittingrule=tree_var%branch_list_level(current_level)% &
          branch_list(current_id)%splittingrule
  if ((splittingrule.ge.1).and. &
      (splittingrule.le.ndim_decisions)) then
   median_index=tree_var%branch_list_level(current_level)% &
          branch_list(current_id)%median_index
   median_value=tree_var%branch_list_level(current_level)% &
          branch_list(current_id)%median_value
   if ((median_index.ge.1).and. &
       (median_index.lt.current_ndata)) then

    prev_id=current_id
    prev_level=current_level

    if (data_decision(splittingrule).le.median_value) then
     current_id=tree_var%branch_list_level(prev_level)% &
          branch_list(prev_id)%child1_id
     current_level=tree_var%branch_list_level(prev_level)% &
          branch_list(prev_id)%child_level
    else if (data_decision(splittingrule).gt.median_value) then
     current_id=tree_var%branch_list_level(prev_level)% &
          branch_list(prev_id)%child2_id
     current_level=tree_var%branch_list_level(prev_level)% &
          branch_list(prev_id)%child_level
    else
     print *,"data_decision bust"
     stop
    endif

    if (current_level.eq.prev_level+1) then
     ! do nothing
    else
     print *,"current_level invalid(3)"
     print *,"current_level=",current_level
     print *,"prev_level=",prev_level
     stop
    endif

    previous_ndata=current_ndata
    current_ndata=tree_var%branch_list_level(current_level)% &
         branch_list(current_id)%ndata
    if ((current_ndata.le.previous_ndata/2+1).and. &
        (current_ndata.ge.1).and. &
        (current_ndata.ge.previous_ndata/2)) then
     !do nothing
    else
     print *,"current_ndata too large"
     stop
    endif
   else
    print *,"median_index invalid"
    stop
   endif
  else
   print *,"splittingrule invalid"
   stop
  endif
 enddo !do while (current_ndata.ge.2)

 if (current_ndata.eq.1) then
  ! do nothing
 else
  print *,"current_ndata invalid"
  stop
 endif

 do dir=1,ndim_classify
  data_classified(dir)=tree_var%branch_list_level(current_level)% &
     branch_list(current_id)%data_classify(1,dir)
 enddo

end subroutine decision_tree_predict

end module global_utility_module

