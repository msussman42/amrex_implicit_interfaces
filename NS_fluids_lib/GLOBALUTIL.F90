#undef BL_LANG_CC
#define BL_LANG_FORT

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



      MODULE LagrangeInterpolationPolynomial
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


      END MODULE LagrangeInterpolationPolynomial

       MODULE LegendreNodes
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
         if (abs(sum-exact).gt.1.0E-13) then
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
          if (abs(sum-exact).gt.1.0E-13) then
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
          if (abs(sum-exact).gt.1.0E-10) then
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
 
      END MODULE LegendreNodes


module global_utility_module

implicit none

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
          (probtype.eq.414)) then! called from MITSUHIRO_MELTING.F90 

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
        if (yblob3.le.zero) then
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
       print *,"expecting probtype=55, 59, or 710 in outside_temperature"
       stop
      endif
 
      return 
      end subroutine outside_temperature



      subroutine get_rigid_velocity( &
        FSI_prescribed_flag, &
        color,dir,vel,xrigid, &
        blob_array,blob_array_size,num_colors,num_elements_blobclass)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(out) :: FSI_prescribed_flag
      INTEGER_T, intent(in) :: color,dir,blob_array_size,num_colors
      INTEGER_T, intent(in) :: num_elements_blobclass
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

      !blob_matrix,blob_RHS,blob_velocity,
      !blob_integral_momentum,blob_energy,
      !blob_mass_for_velocity (3 comp)
      !blob_volume, 
      !blob_center_integral,blob_center_actual
      !blob_perim, blob_perim_mat, blob_triple_perim, 
      if (num_elements_blobclass.ne. &
          3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
          2*(2*SDIM)+1+ &
          3+1+2*SDIM+1+nmat+nmat*nmat) then
       print *,"num_elements_blobclass invalid"
       stop
      endif

      FSI_prescribed_flag=0

      cen_comp=3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
        2*(2*SDIM)+1+ &
        3+1+SDIM
      vel_comp=3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)
      vol_comp=vel_comp+3*(2*SDIM)+2*(2*SDIM)+1
      
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
        print *,"mag invalid"
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

      REAL_T vin(SDIM)
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
      else
       do dir=1,SDIM
        vin(dir)=zero
       enddo
      endif

      return 
      end subroutine normalize_vector


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

      REAL_T a(3),b(3),c(3)

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
          if (abs(holdvalue).gt.1.0E-12) then
           print *,"inverse failed1"
           print *,"AAhold="
           call print_matrix(AAhold,numelem)
           print *,"xx="
           call print_matrix(xx,numelem)
           stop
          endif
         else if (i.eq.j) then
          if (abs(holdvalue-one).gt.1.0E-12) then
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
       DIMS(dist), &
       i,j,k,im,LS)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T i,j,k,im,nhalf,dir
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T LS
      INTEGER_T DIMDEC(dist)
      REAL_T dist(DIMV(dist),num_materials)
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
          print *,"mag invalid"
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

      subroutine checkbound(lo,hi, &
       DIMS(data), &
       ngrow,dir,id)
      IMPLICIT NONE

      INTEGER_T, intent(in) ::  lo(SDIM), hi(SDIM)
      INTEGER_T, intent(in) ::  DIMDEC(data)
      INTEGER_T, intent(in) ::  ngrow,dir,id

      INTEGER_T ii(SDIM)

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
        print *,"dir=",dir
        print *,"dir2=",dir2
        stop
       endif
       ii(dir2)=0
      enddo
      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ii(dir+1)=1
      else if (dir.eq.-1) then
       ! do nothing
      else
       print *,"dir invalid checkbound"
       stop
      endif

      do dir2=1,SDIM
       if (lo(dir2).lt.0) then
        print *,"lo invalid in checkbound id=",id
        print *,"dir2,dataxlo ",dir2,lodata(dir2)
        print *,"dir2,dataxhi ",dir2,hidata(dir2)
        print *,"dir2,lo,ngrow ",dir2,lo(dir2),ngrow
        print *,"dir2,hi,ngrow ",dir2,hi(dir2),ngrow
        print *,"dir=",dir
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
         print *,"dir=",dir
         print *,"dir2=",dir2
         stop
        endif
        if (hidata(dir2).lt.hi(dir2)+ngrow+ii(dir2)) then
         print *,"hi mismatch id=",id
         print *,"datalo,datahi ",lodata(dir2),hidata(dir2)
         print *,"ii(dir2) ",ii(dir2)
         print *,"lo,hi,ngrow ",lo(dir2),hi(dir2),ngrow
         print *,"dir=",dir
         print *,"dir2=",dir2
         stop
        endif

      enddo ! dir2

      return
      end subroutine checkbound


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

      INTEGER_T sdim_parm,dir,bfact,nhalf0
      REAL_T dx(SDIM)
      REAL_T phi0,dist
      REAL_T nn(sdim_parm)
      REAL_T xsten0(-nhalf0:nhalf0,sdim_parm)
      REAL_T x(sdim_parm)
 
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

      EPS=1.0E-6
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

      EPS=1.0E-6
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
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
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

      REAL_T area
      REAL_T areacen(SDIM)
      INTEGER_T nhalf
      REAL_T xsten(-nhalf:nhalf,SDIM)
      INTEGER_T rzflag
      INTEGER_T dir,side,iside,itan,jtan,dir2
      REAL_T RR,R1,R2,RC

      if (dir.eq.0) then
       itan=2
       jtan=SDIM
      else if (dir.eq.1) then
       itan=1
       jtan=SDIM
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       itan=1
       jtan=2
      else
       print *,"dir invalid gridarea"
       stop
      endif
 
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

        ! 1<=normdir<=sdim
      subroutine gridstenMAC(x,xlo,i,j,k,fablo,bfact,dx,nhalf,normdir)
      IMPLICIT NONE 

      INTEGER_T, intent(in) :: nhalf,normdir
      REAL_T, intent(out) :: x(-nhalf:nhalf,SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM)
      INTEGER_T, intent(in) :: i,j,k,bfact
      INTEGER_T dir,icrit
      REAL_T, dimension(:), allocatable :: xsub

      if (bfact.lt.1) then
       print *,"bfact invalid23"
       stop
      endif
      if ((normdir.lt.1).or.(normdir.gt.SDIM)) then
       print *,"normdir invalid"
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
       if (dir.eq.normdir) then
        call gridsten1DMAC(xsub,xlo,icrit,fablo,bfact,dx,dir,nhalf)
       else
        call gridsten1D(xsub,xlo,icrit,fablo,bfact,dx,dir,nhalf)
       endif
       do icrit=-nhalf,nhalf
        x(icrit,dir)=xsub(icrit)
       enddo

       deallocate(xsub)
      enddo ! dir

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



       ! 1<=normdir<=sdim
      subroutine gridstenMAC_level(x,i,j,k,level,nhalf,normdir)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhalf,normdir
      REAL_T, intent(out) :: x(-nhalf:nhalf,SDIM)
      INTEGER_T, intent(in) :: i,j,k,level
      INTEGER_T isten,dir,ii,jj,kk

      ii=0
      jj=0
      kk=0
      if (normdir.eq.1) then
       ii=1
      else if (normdir.eq.2) then
       jj=1
      else if ((normdir.eq.3).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif
 
      if (nhalf.lt.0) then
       print *,"nhalf invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.cache_max_level)) then
       print *,"level invalid gridstenMAC_level"
       stop
      endif
      if ((2*i-ii-nhalf.lt.cache_index_low).or. &
          (2*i-ii+nhalf.gt.cache_index_high)) then
       print *,"i out of range gridstenMAC_level"
       print *,"i,nhalf,level,normdir ",i,nhalf,level,normdir
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
      do isten=-nhalf,nhalf
       dir=1
       x(isten,dir)=grid_cache(level,2*i-ii+isten,dir)
       dir=2
       x(isten,dir)=grid_cache(level,2*j-jj+isten,dir)
       if (SDIM.eq.3) then
        dir=SDIM
        x(isten,dir)=grid_cache(level,2*k-kk+isten,dir)
       endif
      enddo
       
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


      subroutine coarse_subelement_stencilMAC( &
       ifine,jfine,kfine,stenlo,stenhi,bfact_c,bfact_f,dir)
      IMPLICIT NONE

      INTEGER_T dir
      INTEGER_T ifine,jfine,kfine
      INTEGER_T fine_index(3)
      INTEGER_T stenlo(3),stenhi(3)
      INTEGER_T bfact_c,bfact_f
      INTEGER_T dir2,denom

      if (bfact_c.lt.1) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid coarse subelement stencil mac"
       print *,"dir=",dir
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
        if (dir2.eq.dir+1) then
         if (2*(fine_index(dir2)/2).ne.fine_index(dir2)) then
          stenhi(dir2)=stenhi(dir2)+1
         endif
        endif
       else if (bfact_c.gt.1) then
        stenlo(dir2)=DIV_FLOOR(fine_index(dir2),2*bfact_c)
        stenlo(dir2)=stenlo(dir2)*bfact_c
        stenhi(dir2)=stenlo(dir2)+bfact_c-1
        if (dir2.eq.dir+1) then
         denom=2*bfact_c
         if (denom*(fine_index(dir2)/denom).eq.fine_index(dir2)) then
          stenlo(dir2)=fine_index(dir2)/2
          stenhi(dir2)=stenlo(dir2)
         else
          stenhi(dir2)=stenhi(dir2)+1
         endif
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
       ic,jc,kc,stenlo,stenhi,bfact_c,bfact_f,dir)
      IMPLICIT NONE

      INTEGER_T ic,jc,kc,dir
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
      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid fine subelement stencil mac"
       print *,"dir=",dir
       print *,"ic,jc,kc=",ic,jc,kc
       print *,"bfactc, bfactf=",bfact_c,bfact_f
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
        if (dir2.eq.dir) then
         stenhi(dir2)=stenlo(dir2)
        endif
       else if (bfact_c.gt.1) then
        stenlo(dir2)=DIV_FLOOR(coarse_index(dir2),bfact_c)
        stenlo(dir2)=stenlo(dir2)*bfact_c*2
        stenhi(dir2)=stenlo(dir2)+bfact_c*2-1
        if (dir2.eq.dir) then
         if (bfact_c*(coarse_index(dir2)/bfact_c).eq.coarse_index(dir2)) then
          stenlo(dir2)=2*coarse_index(dir2)
          stenhi(dir2)=stenlo(dir2)
         else
          stenlo(dir2)=DIV_FLOOR(coarse_index(dir2),bfact_c)
          stenlo(dir2)=stenlo(dir2)*bfact_c*2
          stenhi(dir2)=stenlo(dir2)+bfact_c*2
         endif
        endif
       else
        print *,"bfact_c invalid"
        stop
       endif
      enddo ! dir2
      
      return
      end subroutine fine_subelement_stencilMAC


       ! 0<=dir<sdim
      subroutine growntileboxMAC( &
       tilelo,tilehi,fablo,fabhi,growlo,growhi,ng,dir)
      IMPLICIT NONE

      INTEGER_T, intent(in) :: dir
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(out) :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: ng
      INTEGER_T dir2

      growlo(3)=0
      growhi(3)=0

      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir invalid growntilebox mac"
       stop
      endif 
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
      if (tilehi(dir+1).eq.fabhi(dir+1)) then
       growhi(dir+1)=growhi(dir+1)+1
      endif
      
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
       nvar,bfact,gridtype, &
       stenhi,dx,xfine,fcoarse,fxfine,caller_id)
      use LagrangeInterpolationPolynomial
      use LegendreNodes

      IMPLICIT NONE
      INTEGER_T nvar
      INTEGER_T bfact
      INTEGER_T gridtype ! 0:ggg; 1:lgg; 2:glg; 3:ggl
      INTEGER_T stenhi(SDIM)
      REAL_T dx(SDIM)
      REAL_T xfine(SDIM)
      REAL_T fcoarse(D_DECL(0:stenhi(1),0:stenhi(2),0:stenhi(3)),nvar)
      REAL_T fxfine(nvar)
 
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

      INTERP_TOL=1.0E-10

      ii=0
      jj=0
      kk=0
      if (gridtype.eq.0) then
       ! do nothing
      else if (gridtype.eq.1) then
       ii=1
      else if (gridtype.eq.2) then
       jj=1
      else if ((gridtype.eq.3).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"gridtype invalid"
       stop
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

      if ((gridtype.lt.0).or.(gridtype.gt.SDIM)) then
       print *,"gridtype invalid"
       stop
      endif 

      do dir=1,SDIM
       if (dir.eq.gridtype) then
        if (stenhi(dir).ne.bfact) then
         print *,"stenhi invalid"
         stop
        endif
       else
        if (stenhi(dir).ne.bfact-1) then
         print *,"stenhi invalid"
         stop
        endif
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

      if (gridtype.eq.0) then
       ! do nothing
      else if ((gridtype.ge.1).and.(gridtype.le.SDIM)) then
 
       allocate(ypointsGL(0:bfact))
       allocate(bwGL(0:bfact))
       allocate(tempGL(0:bfact))

       do i=0,bfact
        ypointsGL(i)=(yGL(i)+one)*bfact*dx(gridtype)/two
       enddo
       call BarycentricWeights(bfact,ypointsGL,bwGL)
       call LagrangeInterpolatingPolynomial(bfact, &
         xfine(gridtype),ypointsGL,bwGL,tempGL)
       do i=0,bfact
        lg(i,gridtype)=tempGL(i)
       enddo

       deallocate(tempGL)
       deallocate(bwGL)
       deallocate(ypointsGL)

      else
       print *,"gridtype invalid"
       stop
      endif

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
       print *,"gridtype=",gridtype
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
       conservative_div_uu, &
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

      INTEGER_T, intent(in) :: conservative_div_uu
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

       else if ((operation_flag.ge.0).and. &
                (operation_flag.le.5)) then
        ! do nothing
       else if ((operation_flag.eq.10).or. &
                (operation_flag.eq.11)) then
        ! do nothing
       else if (operation_flag.eq.7) then ! advection
        ! do nothing
       else if (operation_flag.eq.8) then ! coupling
        ! do nothing
       else if (operation_flag.eq.9) then ! den: cell->mac
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
         if ((conservative_div_uu.eq.1).or. &
             (conservative_div_uu.eq.2)) then
          dest_interp(i1)=dest_interp(i1)*vel(i1)
         else if (conservative_div_uu.eq.0) then
          ! do nothing
         else
          print *,"conservative_div_uu invalid"
          stop
         endif
        else if (nc.eq.SDIM+1) then ! density: NONCONSERVATIVE
         ! do nothing
        else if (nc.eq.SDIM+2) then ! temperature: NONCONSERVATIVE
         ! do nothing
        else
         print *,"nc invalid in lineGRAD"
         stop
        endif
        dest_grad(i1)=zero
       enddo ! i1=0..bfact

      else if ((operation_flag.ge.0).and. &
               (operation_flag.le.6)) then

       ! do nothing

      else if ((operation_flag.eq.10).or. &
               (operation_flag.eq.11)) then

       ! do nothing

      else if (operation_flag.eq.8) then ! coupling terms

       ! do nothing

      else if (operation_flag.eq.9) then ! den: cell->mac

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

        else if ((operation_flag.ge.0).and. &
                 (operation_flag.le.5)) then
         ! do nothing
        else if ((operation_flag.eq.10).or. &
                 (operation_flag.eq.11)) then
         ! do nothing
        else if (operation_flag.eq.7) then
         ! do nothing (advection)
        else if (operation_flag.eq.8) then
         ! do nothing (coupling)
        else if (operation_flag.eq.9) then
         ! do nothing (den: cell-to-mac)
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

      enddo ! isten

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

      subroutine bilinear_interp_stencil(data_stencil,wt_dist, &
                      ncomp,data_interp)
      use probcommon_module
      IMPLICIT NONE

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
       ngrow,dir,id, &
       verbose, &
       force_check, &
       gridno,ngrid,level,finest_level, &
       mf,DIMS(mf))

      IMPLICIT NONE

      INTEGER_T datatype
      REAL_T warning_cutoff
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(SDIM),growhi(SDIM)
      INTEGER_T growlotest(3),growhitest(3)
      INTEGER_T bfact
      REAL_T dx(SDIM)
      INTEGER_T scomp,ncomp,ndefined
      INTEGER_T ngrow
      INTEGER_T dir
      INTEGER_T id
      INTEGER_T verbose
      INTEGER_T force_check
      INTEGER_T gridno,ngrid,level,finest_level
      INTEGER_T DIMDEC(mf)
      REAL_T mf(DIMV(mf),ndefined)
      INTEGER_T i,j,k,ii,jj,kk,dir2
      INTEGER_T n
      INTEGER_T n_singlelayer,n_interior,n_side,n_corner
      REAL_T sum_interior,sum_side,sum_corner,sum_singlelayer
      REAL_T max_interior,max_side,max_corner,max_singlelayer
      REAL_T val
      INTEGER_T noutside,noutside_single
      REAL_T critical_cutoff

      critical_cutoff=1.0D+99

      if (bfact.lt.1) then
       print *,"bfact invalid36"
       stop
      endif

      if (((verbose.eq.0).or.(verbose.eq.1)).and.(force_check.eq.0)) then
       ! do nothing
      else if ((verbose.eq.2).or.(force_check.eq.1)) then
       call FLUSH(6) ! unit=6 (screen)

       ii=0
       jj=0
       kk=0

       if (datatype.eq.0) then
        if (dir.eq.-1) then
         ! do nothing
        else if (dir.eq.0) then
         ii=1
        else if (dir.eq.1) then
         jj=1
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         kk=1
        else
         print *,"dir invalid aggressive worker"
         stop
        endif
       else if (datatype.eq.1) then

        if (dir.eq.0) then
         ii=1
        else if (dir.eq.1) then
         jj=1
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         kk=1
        else
         print *,"dir invalid aggressive worker"
         stop
        endif

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
        call checkbound(fablo,fabhi, &
         DIMS(mf), &
         ngrow,dir,id)
       else if (datatype.eq.2) then
        call checkbound(fablo,fabhi, &
         DIMS(mf), &
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
        else if ((dir.ge.0).and.(dir.lt.SDIM)) then
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

        do i=growlotest(1),growhitest(1)
        do j=growlotest(2),growhitest(2)
        do k=growlotest(3),growhitest(3)
         noutside=0
         noutside_single=0
         if ((i.lt.fablo(1)).or.(i.gt.fabhi(1)+ii)) then
          noutside=noutside+1
         endif
         if ((i.lt.fablo(1)-1).or.(i.gt.fabhi(1)+ii+1)) then
          noutside_single=noutside_single+1
         endif
         if ((j.lt.fablo(2)).or.(j.gt.fabhi(2)+jj)) then
          noutside=noutside+1
         endif
         if ((j.lt.fablo(2)-1).or.(j.gt.fabhi(2)+jj+1)) then
          noutside_single=noutside_single+1
         endif
         if (SDIM.eq.3) then
          if ((k.lt.fablo(SDIM)).or.(k.gt.fabhi(SDIM)+kk)) then
           noutside=noutside+1
          endif
          if ((k.lt.fablo(SDIM)-1).or.(k.gt.fabhi(SDIM)+kk+1)) then
           noutside_single=noutside_single+1
          endif
         endif
         val=mf(D_DECL(i,j,k),n+scomp)

         if (val.ge.critical_cutoff) then
          print *,"val overflow val,dir,i,j,k,n,scomp,id ", &
           val,dir,i,j,k,n,scomp,id
          print *,"bfact,level,finest_level ",bfact,level,finest_level
          do dir2=1,SDIM
           print *,"dir2,fablo,fabhi ",dir2,fablo(dir2),fabhi(dir2)
           print *,"dir2,tilelo,tilehi ",dir2,tilelo(dir2),tilehi(dir2)
          enddo
          stop
         else if ((val.lt.critical_cutoff).and. &
                  (val.gt.-critical_cutoff)) then
          ! do nothing
         else if (val.le.-critical_cutoff) then
          print *,"-val overflow val,dir,i,j,k,n,scomp,id ", &
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
      enddo ! dir

      return
      end subroutine containing_cell
    

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
      INTEGER_T nmat,im

      if (nmat.ne.num_materials) then
       print *,"nmat invalid CTML_FSI_flagF"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      CTML_FSI_flagF=0
      do im=1,nmat
       if (FSI_flag(im).eq.4) then
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
      INTEGER_T nmat,im

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
      if (FSI_flag(im).eq.4) then
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


      function is_ice(nmat,im)
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T is_ice
      INTEGER_T nmat,im

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

      INTEGER_T is_FSI_rigid
      INTEGER_T nmat,im

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
      is_FSI_rigid=0
      if (FSI_flag(im).eq.5) then
       is_FSI_rigid=1
      else if ((FSI_flag(im).eq.0).or. &
               (FSI_flag(im).eq.7).or. &
               (FSI_flag(im).eq.1).or. &
               (FSI_flag(im).eq.2).or. &
               (FSI_flag(im).eq.4).or. &
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
       print *,"im invalid17: im=",im
       print *,"nmat=",nmat
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid is_lag_part"
       print *,"nmat=",nmat
       print *,"num_materials=",num_materials
       stop
      endif

      if ((FSI_flag(im).eq.1).or. & ! prescribed rigid solid (PROB.F90)
          (FSI_flag(im).eq.2).or. & ! prescribed rigid solid (CAD)
          (FSI_flag(im).eq.4).or. & ! FSI link w/Kourosh Shoele
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
      INTEGER_T nmat,im

      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid17: im=",im
       print *,"nmat=",nmat
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

      if ((is_rigid(nmat,im).eq.1).and. &
          (CTML_FSI_mat(nmat,im).eq.0)) then
       is_prescribed=1
      else if ((is_rigid(nmat,im).eq.0).or. &
               (CTML_FSI_mat(nmat,im).eq.1)) then
       is_prescribed=0
      else
       print *,"is_rigid or CTML_FSI_mat invalid"
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

      INTEGER_T nmat,iten
      REAL_T LS(nmat)
      REAL_T LS_extend
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

      INTEGER_T nmat,iten
      REAL_T LS(nmat)
      REAL_T NRM(nmat*SDIM)
      REAL_T NRM_extend(SDIM)
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

      subroutine make_mixture_density(massfrac, &
        den_ambient,denvapor)
      IMPLICIT NONE

      REAL_T, intent(in) :: massfrac
      REAL_T, intent(in) :: denvapor
      REAL_T, intent(inout) :: den_ambient

      if ((massfrac.ge.zero).and.(massfrac.le.one)) then
       if (denvapor.gt.zero) then
        if (den_ambient.gt.zero) then
         den_ambient=den_ambient*denvapor/ &
           (denvapor*(one-massfrac)+den_ambient*massfrac)
        else
         print *,"den_ambient invalid"
         stop
        endif
       else
        print *,"denvapor invalid"
        stop
       endif
      else
       print *,"massfrac invalid"
       stop
      endif

      return
      end subroutine make_mixture_density

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

      Varname='X'
      call dumpstring(Varname)
      Varname='Y'
      call dumpstring(Varname)

      if (plot_sdim.eq.3) then
       Varname='Z'
       call dumpstring(Varname)
      endif

      do im=1,ncomp

       write(matstr,'(I3)') im
       do i=1,3
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       ih=1
       Varname='U'
       ih=ih+1
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
       visual_option, &
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
      INTEGER_T, intent(in) :: visual_option
      INTEGER_T, intent(in) :: visual_revolve

      INTEGER_T strandid

      INTEGER_T nwrite3d,nwrite2d,index3d,index2d

      character*3 levstr
      character*5 gridstr
      character*18 filename18
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
      else if (data_dir.eq.SDIM) then
        dir_chars='ND'
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
       write(filename18,'(A10,A3,A5)') 'tempnddata',levstr,gridstr
       open(unit=4,file=filename18)

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

       write(filename18,'(A10,A3,A5)') 'tempnddata',levstr,gridstr
       open(unit=4,file=filename18)
       print *,"filename18 ",filename18

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

       enddo  ! dir

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

        enddo ! k

       enddo
       enddo
!      enddo

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
       enddo

       deallocate(zone3d_gb(iz_gb)%var)
       deallocate(zone2d_gb(iz_gb)%var)

       igrid=igrid+1
       if (igrid.gt.grids_per_level_array(ilev+1)-1) then
        ilev=ilev+1
        igrid=0
       endif

      enddo  ! iz_gb

      deallocate(zone3d_gb)
      deallocate(zone2d_gb)
      deallocate(lo_gb)
      deallocate(hi_gb)

      close(11)
     
      rmcommand='rm tempnddata*'

      print *,"issuing command ",rmcommand

      sysret=0

#ifdef PGIFORTRAN
      call system(rmcommand)
#else
      call execute_command_line(rmcommand,exitstat=sysret)
#endif
      if (sysret.ne.0) then
       print *,"execute_command_line has sysret=",sysret
       stop
      endif

      return
      end subroutine zones_revolve_sanity



end module global_utility_module

