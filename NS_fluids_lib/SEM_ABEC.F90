#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"
#include "BC_TYPES.H"
#include "SEM_F.H"
#include "ArrayLim.H"

#if (BL_SPACEDIM==3)
#define SDIM 3
#elif (BL_SPACEDIM==2)
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


      SUBROUTINE LagrangeInterpolation(M,x,w,y,f,p)  !!!!--- w is barycentric weights 
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
      SUBROUTINE LagrangeInterpolatingPolynomial(N,xpt,x,w,l) !!!!!!--- w is barycentric weights
      IMPLICIT NONE
      INTEGER,INTENT(IN)                         :: N
      REAL_T, INTENT(IN)                  :: xpt
      REAL_T, DIMENSION(0:N), INTENT(IN)  :: x,w
      REAL_T, DIMENSION(0:N), INTENT(OUT) :: l

      INTEGER                                    :: j
      REAL_T                              :: s,t
      LOGICAL                                    :: var1, xMatchesNode

      xMatchesNode = .FALSE.
      DO j = 0,N
         CALL AlmostEqual(xpt,x(j),var1)
         IF (var1) THEN
            l(j)         = 1.0D0
            xMatchesNode = .TRUE.
         END IF
      END DO
      
      IF (.NOT.xMatchesNode) THEN
         s = 0.0D0
         DO j = 0,N
            t    = w(j)/(xpt-x(j))
            l(j) = t
            s    = s + t
         END DO
         DO j = 0,N
            l(j) = l(j)/s
         END DO
      END IF
      RETURN
      END SUBROUTINE LagrangeInterpolatingPolynomial


! Algorithm 36
!lagrange interpolation derivative;
!direct computation of the polynomial derivative in barycentric form
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
      SUBROUTINE LagInterpolantDerivative(M,x,y,f,w,p1) !!!!--- w is barycentric weights
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



      SUBROUTINE PolyDerivativeMatrix(N,x,w,dl)  !!!!!!!--- w is barycentric weights
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
   
      REAL_T, dimension(:,:), allocatable :: cache_gauss 
      REAL_T, dimension(:,:), allocatable :: cache_gauss_lobatto 
      REAL_T, dimension(:,:), allocatable :: cache_gauss_w 
      REAL_T, dimension(:,:), allocatable :: cache_gauss_lobatto_w 

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
      INTEGER                     :: i,j
      REAL_T               :: P1,P2,dl1,dl2

      if (k.eq.0) then ! constant 
          P=one
          dl=zero
      elseif (k.eq.1) then
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

      !return y(0:n),weight(0:n);
      !two sets of y(0:(n+1)/2-1)
      subroutine SplitNodes(n,y)
      INTEGER,INTENT(IN)          :: n
      REAL_T, DIMENSION(0:n), INTENT(OUT)  :: y
      REAL_T ysplit(0:n)
      integer i,i1,bfact,bfactC

      bfact=n+1
      bfactC=(n+1)/2
      if (2*bfactC.ne.bfact) then
       print *,"n invalid"
       stop
      endif
      do i1=0,bfactC-1
       ysplit(i1)=cache_gauss(bfactC,i1)
      enddo
      do i=0,bfactC-1
       y(i)=(ysplit(i)-one)/two
       y(i+bfactC)=y(i)+one
      enddo

      return 
      end subroutine SplitNodes


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
        elseif (n.eq.1) then
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


      !return y(0:n),weight(0:n);
      subroutine SplitLobattoNodes(n,y)
      INTEGER,INTENT(IN)          :: n
      REAL_T, DIMENSION(0:n), INTENT(OUT)  :: y
      REAL_T ysplit(0:n)
      integer i,i1,bfact,bfactC

      bfact=n
      bfactC=n/2
      if (2*bfactC.ne.bfact) then
       print *,"n invalid"
       stop
      endif
      do i1=0,bfactC
       ysplit(i1)=cache_gauss_lobatto(bfactC,i1)
      enddo
      do i=0,bfactC
       y(i)=(ysplit(i)-one)/two
       y(i+bfactC)=y(i)+one
      enddo

      return 
      end subroutine SplitLobattoNodes




       ! y1 and y2 are GL points
      subroutine intersect_table(n1,n2,y1,y2,w)
      IMPLICIT NONE

      integer n1,n2,i,j
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

      integer n1,n2,i,j
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
    
      subroutine do_polyinterp(r,w,data,sum)
      IMPLICIT NONE

      integer i,r
      REAL_T w(0:r)
      REAL_T data(0:r)
      REAL_T sum

      sum=zero
      do i=0,r
       sum=sum+data(i)*w(i)
      enddo

      return
      end subroutine do_polyinterp
       

        ! r is the order of the polynomial interpolant
      subroutine polyinterp_weights(r,x,w,xtarget)
      IMPLICIT NONE

      integer r,i,j
      REAL_T xtarget,lag
      REAL_T x(0:r)
      REAL_T w(0:r)

      do i=0,r
       lag=one
       do j=0,r
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



        ! (r+1)/2 is the order of each polynomial interpolant
      subroutine polyinterp_weights_split(r,x,w,xtarget)
      IMPLICIT NONE

      integer r,i,j,bfact,bfactC
      REAL_T xtarget,lag
      REAL_T x(0:r)
      REAL_T w(0:r)

      bfact=r+1
      bfactC=bfact/2
      if (bfactC*2.ne.bfact) then
       print *,"bfact invalid"
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



      subroutine polyinterp_weights_even(r,x,w,xtarget)
      IMPLICIT NONE

      integer r,i,bfact
      REAL_T xtarget,lag
      REAL_T x(0:r)
      REAL_T w(0:r)
      REAL_T dx,totalw


      bfact=r+1
      dx=two/bfact

      if (bfact.eq.1) then
       w(0)=one
      else if (bfact.gt.1) then
       totalw=zero
       do i=0,r
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
        do i=0,r
         w(i)=w(i)/totalw
        enddo
       endif
      else
       print *,"bfact invalid"
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
       ! 0=sum_i wij - (r+1)lambda
      subroutine polyinterp_Dmatrix(r,x,w)
      IMPLICIT NONE

      integer r,i,j,l,m
      REAL_T x(0:r)
      REAL_T w(0:r,0:r)
      REAL_T sum

      if (r.eq.0) then
       w(0,0)=zero
      else if (r.eq.1) then
       w(0,0)=one/(x(0)-x(1))
       w(0,1)=w(0,0)
       w(1,0)=-w(0,0)
       w(1,1)=w(1,0)
      else if (r.gt.1) then
       do i=0,r
       do j=0,r
        if (i.eq.j) then
         sum=zero
         do l=0,r
          if (l.ne.i) then
           sum=sum+one/(x(i)-x(l))
          endif
         enddo
         w(i,j)=sum
        else if (i.ne.j) then 
         sum=one
         do m=0,r
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
       do j=0,r
        sum=zero
        do i=0,r
         sum=sum+w(i,j) 
        enddo
        sum=sum/(r+one)
        do i=0,r
         w(i,j)=w(i,j)-sum
        enddo 
       enddo
      else
       print *,"r invalid"
       stop
      endif

      return
      end subroutine polyinterp_Dmatrix

      subroutine poly_change_basis(r1,r2,data1,data2,y1,y2)
      IMPLICIT NONE

      integer r1,r2,i
      REAL_T y1(0:r1)
      REAL_T y2(0:r2)
      REAL_T data1(0:r1)
      REAL_T w(0:r1)
      REAL_T data2(0:r2)

      do i=0,r2
       call polyinterp_weights(r1,y1,w,y2(i))
       call do_polyinterp(r1,w,data1,data2(i))
      enddo

      return
      end subroutine poly_change_basis

      subroutine deriv_change_basis(r1,r2,data,datader, &
        wMAT,y1,y2,dx_element)
      IMPLICIT NONE

      integer r1,r2
      REAL_T y1(0:r1)
      REAL_T y2(0:r2)
      REAL_T data(0:r1)
      REAL_T datader1(0:r1)
      REAL_T datader(0:r2)
      REAL_T wMAT(0:r1,0:r1)
      integer i,j
      REAL_T sum,dx_element
      REAL_T w(0:r1)
   
      do i=0,r1
       sum=zero
       do j=0,r1
        sum=sum+data(j)*wMAT(j,i)
       enddo
       datader1(i)=sum
      enddo
      call poly_change_basis(r1,r2,datader1,datader,y1,y2)
      do i=0,r2
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

      return
      end subroutine delete_cache

       ! typ=0  Legendre
       ! typ=1  Clenshaw Curtis 
      subroutine init_cache(r,typ)
      IMPLICIT NONE

      REAL_T tempx(0:r)
      REAL_T tempw(0:r)

      integer r,i,j,typ

      if (r.lt.1) then
       print *,"r invalid"
       stop
      endif

      allocate(cache_gauss(1:r,0:r))
      allocate(cache_gauss_lobatto(1:r,0:r))
      allocate(cache_gauss_w(1:r,0:r))
      allocate(cache_gauss_lobatto_w(1:r,0:r))

      do i=1,r
       if (typ.eq.0) then
        call LegendreGaussLobattoNodesAndWeights(i,tempx,tempw) 
       else if (typ.eq.1) then
        call ClenshawGaussLobattoNodesAndWeights(i,tempx,tempw) 
       else
        print *,"typ invalid"
        stop
       endif

       do j=0,i
        cache_gauss_lobatto(i,j)=tempx(j)
        cache_gauss_lobatto_w(i,j)=tempw(j)
       enddo

       if (typ.eq.0) then
        call LegendreGaussNodesAndWeights(i-1,tempx,tempw) 
       else if (typ.eq.1) then
        call ClenshawGaussNodesAndWeights(i-1,tempx,tempw) 
       else
        print *,"typ invalid"
        stop
       endif

       do j=0,i-1
        cache_gauss(i,j)=tempx(j)
        cache_gauss_w(i,j)=tempw(j)
       enddo
      enddo

      return
      end subroutine init_cache
 
      subroutine sanity_check(rend)
      IMPLICIT NONE
 
      integer typ,r,i,j,k,p,rstart,rend
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
       do r=rstart,rend
        allocate(y(0:r))
        allocate(weight(0:r))
        allocate(w(0:r))
        allocate(wMAT(0:r,0:r))
        allocate(data(0:r))
        allocate(datader(0:r))
        if (typ.eq.0) then
         call LegendreGaussNodesAndWeights(r,y,weight)
        else if (typ.eq.1) then
         call LegendreGaussLobattoNodesAndWeights(r,y,weight)
        else if (typ.eq.2) then
         call ClenshawGaussNodesAndWeights(r,y,weight)
        else if (typ.eq.3) then
         call ClenshawGaussLobattoNodesAndWeights(r,y,weight)
        endif
        do i=0,r-1
         if (y(i).ge.y(i+1)) then
          print *,"nodes out of order"
          stop
         endif
        enddo
        do p=0,r
         do i=0,r
          data(i)=y(i)**p
         enddo
         sum=zero
         do i=0,r
          sum=sum+weight(i)*data(i)
         enddo
         if ((p/2)*2.eq.p) then
          exact=two/(p+one)
         else
          exact=zero
         endif
         if (abs(sum-exact).gt.1.0E-13) then
          print *,"sanity check failed integral"
          print *,"typ,r,p,exact,sum ",typ,r,p,exact,sum
          stop
         endif
         do i=0,r
          xtarget=-one+two*i/r
          call polyinterp_weights(r,y,w,xtarget)
          call do_polyinterp(r,w,data,sum)
          exact=xtarget**p
          if (abs(sum-exact).gt.1.0E-13) then
           print *,"sanity check failed polyinterp"
           print *,"r,p ",r,p
           stop
          endif
          call polyinterp_Dmatrix(r,y,wMAT)
          dx_element=two
          call deriv_change_basis(r,r,data,datader,wMAT,y,y,dx_element)

          call polyinterp_weights(r,y,w,xtarget)
          call do_polyinterp(r,w,datader,sum)
          if (p.eq.0) then
           exact=zero
          else
           exact=p*(xtarget**(p-1))
          endif
          if (abs(sum-exact).gt.1.0E-10) then
           print *,"sanity check failed polyinterp_Dmatrix"
           print *,"r,p,typ ",r,p,typ
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


      MODULE Nodal2DClass
       ! one element
      USE LegendreNodes
      USE LagrangeInterpolationPolynomial
      IMPLICIT NONE

      TYPE DGprecomputed
      INTEGER                                     :: N
      REAL_T, ALLOCATABLE, DIMENSION(:)    :: x,Qweights 
      REAL_T, ALLOCATABLE, DIMENSION(:,:)  :: Dlx 
      REAL_T, ALLOCATABLE, DIMENSION(:)    :: L1,Lminus1
      END TYPE DGprecomputed

      CONTAINS

      SUBROUTINE ConstructNodal2D(N,this,typenodes) 
      IMPLICIT NONE
      INTEGER,INTENT(IN)                             :: N
      TYPE(DGprecomputed), INTENT(OUT)               :: this
      REAL_T, ALLOCATABLE, DIMENSION(:)       :: BaryWeights 
      REAL_T, ALLOCATABLE, DIMENSION(:,:)     :: D
      CHARACTER(LEN=2), INTENT(IN)                   :: typenodes  
      INTEGER                                        :: i,j

      this%N = N

      ALLOCATE(this%x(0:N),this%Qweights(0:N),this%Dlx(0:N,0:N), &
         & this%Lminus1(0:N),this%L1(0:N))

      this%x=0.0D0
      this%Qweights=0.0D0
      this%Dlx=0.0D0
      this%Lminus1 = 0.0D0
      this%L1 = 0.0D0
     
      ALLOCATE(BaryWeights(0:N))
      ALLOCATE(D(0:N,0:N))

      IF (typenodes == 'GL') THEN
          CALL LegendreGaussLobattoNodesAndWeights(N, &
                                        & this%x,this%Qweights)
      ELSEIF (typenodes == 'Ga') THEN    
          CALL LegendreGaussNodesAndWeights(N,this%x,this%Qweights)
      END IF   
   
      CALL BarycentricWeights(N,this%x,BaryWeights)  

      CALL LagrangeInterpolatingPolynomial(N, &
                 & -1.0D0,this%x,BaryWeights,this%Lminus1)
      CALL LagrangeInterpolatingPolynomial(N, &
                 & 1.0D0,this%x,BaryWeights,this%L1)      
 
       
      CALL PolyDerivativeMatrix(N,this%x,BaryWeights,D)                
      DO j = 0, N
        DO i = 0, N
        !!!! the derivative is taken on basis function
           this%Dlx(i,j) = D(j,i)*this%Qweights(j)   
        END DO
      END DO
           
      DEALLOCATE(BaryWeights)
      DEALLOCATE(D)
     
      RETURN
      END SUBROUTINE ConstructNodal2D



      SUBROUTINE DestructNodal2D(this)
      IMPLICIT NONE
      TYPE(DGprecomputed) :: this

      this%N = -1
      DEALLOCATE(this%x,this%Qweights,this%Dlx, &
           & this%Lminus1,this%L1)
       
      END SUBROUTINE DestructNodal2D
      
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
      SUBROUTINE InterpolateToBoundary(N,u,l,ubc)
      IMPLICIT NONE
      INTEGER,INTENT(IN)                                :: N
      REAL_T, DIMENSION(0:N), INTENT(IN)         :: u,l
      REAL_T,INTENT(OUT)                         :: ubc
      INTEGER                                           :: j

      ubc=0.0D0
      do j=0,N
        ubc=ubc+l(j)*u(j)
      end do

      RETURN
      END SUBROUTINE InterpolateToBoundary   
   

!!!!!!
      SUBROUTINE GradientWeakFormCompute(spA,u,IuDw,ubc)
      IMPLICIT NONE
      TYPE(DGprecomputed), INTENT(IN)               :: spA
      REAL_T, DIMENSION(0:spA%N), INTENT(IN) :: u
      REAL_T, DIMENSION(0:spA%N), INTENT(OUT):: IuDw
      REAL_T, DIMENSION(0:1), INTENT(OUT)    :: ubc


      !!!! compute - \int u*(w_i)_x for each w_i      
      IuDw = 0.0D0
      call MxVDerivative(spA%N,spA%Dlx,u,IuDw)      
      IuDw = -IuDw 
          
      
      !!!!! at cell interface
      !!! w_i(-1): precomputed and stored in spA%Lminus1
      !!! w_i(1):  precomputed and stored in spA%L1
      
      !!! u(-1): ubc(0)
      !!! u(1):  ubc(1)
      CALL InterpolateToBoundary(spA%N,u,spA%Lminus1,ubc(0))
      CALL InterpolateToBoundary(spA%N,u,spA%L1,ubc(1))
      
      

      

      RETURN
      END SUBROUTINE GradientWeakFormCompute  


      SUBROUTINE AffineMap(psi,x_k,x_k_minus_one,x)
      IMPLICIT NONE
      REAL_T, INTENT(IN)  :: psi,x_k,x_k_minus_one
      REAL_T, INTENT(OUT) :: x

      x = x_k_minus_one + ((psi+1.0D0)/2.0D0)*(x_k-x_k_minus_one)

      RETURN
      END SUBROUTINE AffineMap   
  

      SUBROUTINE InternalPenalty_GradientFlux(spA,u,xlo,xhi,uxbc)
      IMPLICIT NONE
      TYPE(DGprecomputed), INTENT(IN)               :: spA
      REAL_T, DIMENSION(0:spA%N), INTENT(IN) :: u
      REAL_T, INTENT(IN)                     :: xlo,xhi
      REAL_T, DIMENSION(0:1), INTENT(OUT)    :: uxbc
      REAL_T, ALLOCATABLE, DIMENSION(:)         :: BaryWeights


      ALLOCATE(BaryWeights(0:spA%N))

      CALL BarycentricWeights(spA%N,spA%x,BaryWeights)
      
      call LagInterpolantDerivative(spA%N,spA%x,-1.0D0, &
                  & u,BaryWeights,uxbc(0))      
      call LagInterpolantDerivative(spA%N,spA%x,1.0D0, &
                  & u,BaryWeights,uxbc(1)) 
                  
      uxbc(0)=uxbc(0)*(2.0D0/(xhi-xlo))                       
      uxbc(1)=uxbc(1)*(2.0D0/(xhi-xlo))
      
      DEALLOCATE(BaryWeights)
   
   
      RETURN
      END SUBROUTINE InternalPenalty_GradientFlux     
 
      END MODULE Nodal2DClass



        ! R=G-A U
      subroutine RESID(ncell,lo,hi,problo,probhi, &
        R,U,G,betax,betay,betacell,dx,bfact,hflag,probtype, &
        piecewise,weak_flag)
      IMPLICIT NONE

      integer ncell(SDIM),lo(SDIM),hi(SDIM),i,j,piecewise,weak_flag
      REAL_T problo(SDIM),probhi(SDIM),dx(SDIM)
      integer bfact,hflag,probtype
      REAL_T betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T betacell(lo(1):hi(1),lo(SDIM):hi(SDIM))
      REAL_T R(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T U(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T G(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)

      REAL_T, dimension(:,:), allocatable :: AU
      REAL_T alpha,beta

      allocate(AU(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      if ((weak_flag.lt.0).or.(weak_flag.gt.4)) then
       print *,"weak_flag invalid in RESID"
       stop
      endif

      call ATIMESU(ncell,lo,hi,problo,probhi, &
       dx,bfact,U,betax,betay,betacell,hflag,probtype,AU, &
       piecewise,weak_flag)
      alpha=1.0
      beta=-1.0
      call LINCOMB(lo,hi,G,AU,R,alpha,beta)

      deallocate(AU)

      return
      end subroutine RESID

      subroutine checkzero(val,idx)
      IMPLICIT NONE

      REAL_T val
      integer idx

      if (abs(val).lt.1.0E-14) then
       print *,"value is almost 0"
       print *,"val = ",val
       print *,"idx = ",idx
       stop
      endif

      return
      end



      recursive subroutine RELAX(lo,hi,Z,R,betax,betay, &
       dx,hflag,problo,probhi,probtype,nsmooth,tol,bfact)
      IMPLICIT NONE

      integer bfact,bfactC
      integer i,j,lo(SDIM),hi(SDIM),irelax,nsmooth,hflag
      REAL_T tol
      REAL_T betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T R(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T Z(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T dx(SDIM)
      REAL_T problo(SDIM),probhi(SDIM)
      integer probtype,piecewise

      integer loC(SDIM),hiC(SDIM)
      REAL_T dxC(SDIM)
      integer count_relax,coarsest_nx,ncellC(SDIM)
      integer ncell(SDIM)
      integer dir,local_precond_type,bottomflag,weak_flag
      REAL_T bottom_tol

      REAL_T, dimension(:,:), allocatable :: ZN
      REAL_T, dimension(:,:), allocatable :: RF
      REAL_T, dimension(:,:), allocatable :: RC
      REAL_T, dimension(:,:), allocatable :: ZC
      REAL_T, dimension(:,:), allocatable :: betaxC
      REAL_T, dimension(:,:), allocatable :: betayC

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

      bottom_tol=0.01*tol

      weak_flag=0
      piecewise=1
      count_relax=1
      coarsest_nx=1

      allocate(ZN(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 

      if (hflag.ne.1) then
       print *,"hflag invalid"
       stop
      endif

      call COPYVEC(lo,hi,Z,ZN)

      call JACPRECOND(lo,hi, &
        ZN,R,betax,betay,dx,hflag,probtype,nsmooth,problo,probhi,bfact)

      do dir=1,SDIM 
       ncell(dir)=hi(dir)-lo(dir)+1
       ncellC(dir)=ncell(dir)/2
       loC(dir)=lo(dir)/2
       hiC(dir)=loC(dir)+ncellC(dir)-1
      enddo
      if ((2*ncellC(1).eq.ncell(1)).and. &
          (2*ncellC(2).eq.ncell(2)).and. &
          (ncellC(1).ge.coarsest_nx).and. &
          (ncellC(2).ge.coarsest_nx)) then
       bfactC=bfact
       do dir=1,SDIM
        if ((ncellC(dir)/bfactC)*bfactC.ne.ncellC(dir)) then
         bfactC=bfact/2
         if (bfactC*2.ne.bfact) then
          print *,"bfact must be powers of 2"
          stop
         endif
        endif
       enddo

       do dir=1,SDIM
        dxC(dir)=2.0*dx(dir)
       enddo
       allocate(RF(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
        ! RF=R-A ZN
       call RESID(ncell,lo,hi,problo,probhi, &
        RF,ZN,R,betax,betay,betay,dx,bfact,hflag,probtype, &
        piecewise,weak_flag)

       allocate(ZC(loC(1)-1:hiC(1)+1,loC(SDIM)-1:hiC(SDIM)+1)) 
       allocate(RC(loC(1)-1:hiC(1)+1,loC(SDIM)-1:hiC(SDIM)+1)) 
       allocate(betaxC(loC(1):hiC(1)+1,loC(SDIM):hiC(SDIM))) 
       allocate(betayC(loC(1):hiC(1),loC(SDIM):hiC(SDIM)+1)) 
       call ZAPVEC(loC,hiC,ZC)
       call RESTRICT(RF,RC,betax,betay,betaxC,betayC, &
        lo,hi,loC,hiC,bfact,bfactC)

       do irelax=1,count_relax
        call RELAX(loC,hiC,ZC,RC,betaxC,betayC, &
         dxC,hflag,problo,probhi,probtype,nsmooth,tol,bfactC)
       enddo
       if ((weak_flag.ne.0).or.(piecewise.ne.1)) then
        print *,"weak_flag or piecewise became corrupted"
        stop
       endif
       call set_boundary(loC,hiC,dxC,hflag,probtype,ZC,bfactC, &
        problo,probhi,piecewise,weak_flag)
       call PROLONGATE(ZN,ZC,lo,hi,loC,hiC,bfact,bfactC)

       deallocate(ZC)
       deallocate(RC)
       deallocate(RF)
       deallocate(betaxC)
       deallocate(betayC)

      else
       bottomflag=1
       local_precond_type=1  ! Jacobi
       call pcg(ncell,lo,hi,bfact, &
        ZN,R,betax,betay,betay,dx,hflag,local_precond_type,  &
        problo,probhi,probtype,bottom_tol,piecewise,weak_flag, &
        bottomflag)
      endif
      call JACPRECOND(lo,hi, &
        ZN,R,betax,betay,dx,hflag,probtype,nsmooth,problo,probhi,bfact)

      call COPYVEC(lo,hi,ZN,Z)

      deallocate(ZN) 

      return
      end subroutine RELAX

      subroutine abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
      IMPLICIT NONE

      integer i,j,i1,j1,i2,j2,bfact,dir
      integer lo_e(SDIM),hi_e(SDIM)
      integer lo(SDIM),hi(SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

      i2=(i-lo_e(1))*bfact+i1+lo(1)
      j2=(j-lo_e(2))*bfact+j1+lo(2)

      if (bfact.eq.1) then
       if (((i2.ne.i).or.(j2.ne.j)).and.(dir.eq.-1)) then
        print *,"abs_index bust"
        print *,"lo,hi ",lo(1),hi(1),lo(2),hi(2)
        print *,"bfact= ",bfact
        print *,"lo_e,hi_e ",lo_e(1),hi_e(1),lo_e(2),hi_e(2)
        print *,"i1,j1 ",i1,j1
        stop
       endif
      endif

      return
      end subroutine abs_index



      subroutine RESTRICT(RF,RC,betax,betay,betaxC,betayC, &
        lo,hi,loC,hiC,bfact,bfactC)
      use LegendreNodes

      IMPLICIT NONE

      integer ncell(SDIM)
      integer ncell_e(SDIM)
      integer ncellC(SDIM)
      integer ncellC_e(SDIM)
      integer lo(SDIM),hi(SDIM)
      integer loC(SDIM),hiC(SDIM)
      integer bfact,bfactC,i,j,dir
      REAL_T RF(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T RC(loC(1)-1:hiC(1)+1,loC(SDIM)-1:hiC(SDIM)+1)
      REAL_T RSOURCE,BETASOURCE
      REAL_T yG(0:2*bfact)
      REAL_T yE(0:2*bfact)
      REAL_T yGL(0:2*bfact)
      REAL_T yEL(0:2*bfact)
      integer i1C,j1C,i1,j1,i2,j2,if,jf,i1f,j1f

      integer lo_e(SDIM),hi_e(SDIM)
      integer loC_e(SDIM),hiC_e(SDIM)
      REAL_T betaxC(loC(1):hiC(1)+1,loC(SDIM):hiC(SDIM))
      REAL_T betayC(loC(1):hiC(1),loC(SDIM):hiC(SDIM)+1)
      REAL_T betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T, dimension(:,:), allocatable :: wtx
      REAL_T, dimension(:,:), allocatable :: wty
      REAL_T, dimension(:,:), allocatable :: wtcen
      REAL_T, dimension(:,:), allocatable :: intvolG
      REAL_T, dimension(:,:), allocatable :: intvolGL
      REAL_T localvol

      if ((bfact.lt.1).or.(bfactC.lt.1)) then
       print *,"bfact or bfactC invalid"
       stop
      endif

      do dir=1,SDIM
       lo_e(dir)=lo(dir)/bfact
       hi_e(dir)=((hi(dir)+1)/bfact)-1
       ncell(dir)=hi(dir)-lo(dir)+1
       ncell_e(dir)=hi_e(dir)-lo_e(dir)+1
      enddo
      do dir=1,SDIM
       loC_e(dir)=loC(dir)/bfactC
       hiC_e(dir)=((hiC(dir)+1)/bfactC)-1
       ncellC(dir)=hiC(dir)-loC(dir)+1
       ncellC_e(dir)=hiC_e(dir)-loC_e(dir)+1
       if (2*ncellC(dir).ne.ncell(dir)) then
        print *,"ncell bust"
        stop
       endif
       if (bfact.eq.bfactC) then
        if (2*ncellC_e(dir).ne.ncell_e(dir)) then
         print *,"ncell_e bust"
         stop
        endif
       else if (bfact.eq.2*bfactC) then
        if (ncellC_e(dir).ne.ncell_e(dir)) then
         print *,"ncell_e bust"
         stop
        endif
       else
        print *,"bfact bust"
        stop
       endif
      enddo ! dir

      allocate(wtx(loC(1):hiC(1)+1,loC(SDIM):hiC(SDIM))) 
      allocate(wty(loC(1):hiC(1),loC(SDIM):hiC(SDIM)+1)) 
      allocate(wtcen(loC(1):hiC(1),loC(SDIM):hiC(SDIM))) 

      do i1=0,bfactC-1
       yG(i1)=cache_gauss(bfactC,i1)
      enddo
      do i1=0,bfactC
       yGL(i1)=cache_gauss_lobatto(bfactC,i1)
      enddo
 
      if (bfact.eq.bfactC) then
       call SplitNodes(2*bfact-1,yE)
       call SplitLobattoNodes(2*bfact,yEL)
       allocate(intvolG(0:bfactC-1,0:2*bfact-1))
       allocate(intvolGL(0:bfactC,0:2*bfact))
       call intersect_table(bfactC,2*bfact,yGL,yEL,intvolG)
       call intersect_tableGL(bfactC-1,2*bfact-1,yG,yE,intvolGL)
      else if (bfact.eq.2*bfactC) then
       do i1=0,bfact-1
        yE(i1)=cache_gauss(bfact,i1)
       enddo
       do i1=0,bfact
        yEL(i1)=cache_gauss_lobatto(bfact,i1)
       enddo
       allocate(intvolG(0:bfactC-1,0:bfact-1))
       allocate(intvolGL(0:bfactC,0:bfact))
       call intersect_table(bfactC,bfact,yGL,yEL,intvolG)
       call intersect_tableGL(bfactC-1,bfact-1,yG,yE,intvolGL)
      else
       print *,"bfact invalid"
       stop
      endif


      do i=loC(1),hiC(1)
      do j=loC(2),hiC(2)
       wtcen(i,j)=zero
       RC(i,j)=zero
      enddo
      enddo

      do i=loC_e(1),hiC_e(1)
      do j=loC_e(SDIM),hiC_e(SDIM)
       do j1C=0,bfactC-1
       do i1C=0,bfactC-1
        do j1=0,2*bfactC-1
        do i1=0,2*bfactC-1
         dir=-1
         if (bfact.eq.bfactC) then
          if=2*i
          jf=2*j
          i1f=i1
          j1f=j1
          if (i1.ge.bfactC) then
           if=if+1
           i1f=i1-bfactC
          endif
          if (j1.ge.bfactC) then
           jf=jf+1
           j1f=j1-bfactC
          endif
         else if (bfact.eq.2*bfactC) then
          if=i
          jf=j
          i1f=i1
          j1f=j1
         else
          print *,"bfact invalid"
          stop
         endif 
         call abs_index(if,jf,i1f,j1f,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
         localvol=intvolG(i1C,i1)*intvolG(j1C,j1)
         RSOURCE=localvol*RF(i2,j2)
         call abs_index(i,j,i1C,j1C,i2,j2,loC_e,hiC_e, &
            bfactC,loC,hiC,dir)
         RC(i2,j2)=RC(i2,j2)+RSOURCE
         wtcen(i2,j2)=wtcen(i2,j2)+localvol
        enddo 
        enddo 
       enddo 
       enddo 
      enddo
      enddo
      do i=loC(1),hiC(1)
      do j=loC(2),hiC(2)
       localvol=wtcen(i,j)
       if (localvol.le.zero) then
        print *,"localvol invalid"
        stop
       endif
       RC(i,j)=RC(i,j)/localvol
      enddo
      enddo
       
      do i=loC(1),hiC(1)+1
      do j=loC(2),hiC(2)
       wtx(i,j)=zero
       betaxC(i,j)=zero
      enddo
      enddo

      do i=loC_e(1),hiC_e(1)
      do j=loC_e(SDIM),hiC_e(SDIM)
       do j1C=0,bfactC-1
       do i1C=0,bfactC
        do j1=0,2*bfactC-1
        do i1=0,2*bfactC
         dir=0
         if (bfact.eq.bfactC) then
          if=2*i
          jf=2*j
          i1f=i1
          j1f=j1
          if (i1.ge.bfactC) then
           if=if+1
           i1f=i1-bfactC
          endif
          if (j1.ge.bfactC) then
           jf=jf+1
           j1f=j1-bfactC
          endif
         else if (bfact.eq.2*bfactC) then
          if=i
          jf=j
          i1f=i1
          j1f=j1
         else
          print *,"bfact invalid"
          stop
         endif 
         call abs_index(if,jf,i1f,j1f,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
         localvol=intvolGL(i1C,i1)*intvolG(j1C,j1)
         BETASOURCE=localvol*betax(i2,j2)
         call abs_index(i,j,i1C,j1C,i2,j2,loC_e,hiC_e, &
             bfactC,loC,hiC,dir)
         betaxC(i2,j2)=betaxC(i2,j2)+BETASOURCE
         wtx(i2,j2)=wtx(i2,j2)+localvol
        enddo 
        enddo 
       enddo 
       enddo 
      enddo
      enddo
      do i=loC(1),hiC(1)+1
      do j=loC(2),hiC(2)
       localvol=wtx(i,j)
       if (localvol.le.zero) then
        print *,"localvol invalid"
        stop
       endif
       betaxC(i,j)=betaxC(i,j)/localvol
      enddo
      enddo

      do i=loC(1),hiC(1)
      do j=loC(2),hiC(2)+1
       wty(i,j)=zero
       betayC(i,j)=zero
      enddo
      enddo

      do i=loC_e(1),hiC_e(1)
      do j=loC_e(SDIM),hiC_e(SDIM)
       do j1C=0,bfactC
       do i1C=0,bfactC-1
        do j1=0,2*bfactC
        do i1=0,2*bfactC-1
         dir=1
         if (bfact.eq.bfactC) then
          if=2*i
          jf=2*j
          i1f=i1
          j1f=j1
          if (i1.ge.bfactC) then
           if=if+1
           i1f=i1-bfactC
          endif
          if (j1.ge.bfactC) then
           jf=jf+1
           j1f=j1-bfactC
          endif
         else if (bfact.eq.2*bfactC) then
          if=i
          jf=j
          i1f=i1
          j1f=j1
         else
          print *,"bfact invalid"
          stop
         endif 
         call abs_index(if,jf,i1f,j1f,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
         localvol=intvolG(i1C,i1)*intvolGL(j1C,j1)
         BETASOURCE=localvol*betay(i2,j2)
         call abs_index(i,j,i1C,j1C,i2,j2,loC_e,hiC_e, &
             bfactC,loC,hiC,dir)
         betayC(i2,j2)=betayC(i2,j2)+BETASOURCE
         wty(i2,j2)=wty(i2,j2)+localvol
        enddo 
        enddo 
       enddo 
       enddo 
      enddo
      enddo
      do i=loC(1),hiC(1)
      do j=loC(2),hiC(2)+1
       localvol=wty(i,j)
       if (localvol.le.zero) then
        print *,"localvol invalid"
        stop
       endif
       betayC(i,j)=betayC(i,j)/localvol
      enddo
      enddo

      deallocate(wtx)
      deallocate(wty)
      deallocate(wtcen)
      deallocate(intvolG)
      deallocate(intvolGL)

      return
      end subroutine RESTRICT



      subroutine PROLONGATE(ZF,ZC,lo,hi,loC,hiC,bfact,bfactC)
      use LegendreNodes
      IMPLICIT NONE

      integer ncell(SDIM)
      integer ncell_e(SDIM)
      integer ncellC(SDIM)
      integer ncellC_e(SDIM)
      integer lo(SDIM),hi(SDIM)
      integer loC(SDIM),hiC(SDIM)
      integer bfact,bfactC,i,j,dir
      REAL_T ZF(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T ZC(loC(1)-1:hiC(1)+1,loC(SDIM)-1:hiC(SDIM)+1)
      REAL_T ZSOURCE
      REAL_T yG(0:2*bfact)
      REAL_T yE(0:2*bfact)
      REAL_T yGL(0:2*bfact)
      REAL_T yEL(0:2*bfact)
      integer i1C,j1C,i1,j1,i2,j2,if,jf,i1f,j1f

      integer lo_e(SDIM),hi_e(SDIM)
      integer loC_e(SDIM),hiC_e(SDIM)
      REAL_T, dimension(:,:), allocatable :: ZCOR
      REAL_T, dimension(:,:), allocatable :: wtcen
      REAL_T, dimension(:,:), allocatable :: intvolG
      REAL_T localvol



      if ((bfact.lt.1).or.(bfactC.lt.1)) then
       print *,"bfact or bfactC invalid"
       stop
      endif

      do dir=1,SDIM
       lo_e(dir)=lo(dir)/bfact
       hi_e(dir)=((hi(dir)+1)/bfact)-1
       ncell(dir)=hi(dir)-lo(dir)+1
       ncell_e(dir)=hi_e(dir)-lo_e(dir)+1
      enddo
      do dir=1,SDIM
       loC_e(dir)=loC(dir)/bfactC
       hiC_e(dir)=((hiC(dir)+1)/bfactC)-1
       ncellC(dir)=hiC(dir)-loC(dir)+1
       ncellC_e(dir)=hiC_e(dir)-loC_e(dir)+1
       if (2*ncellC(dir).ne.ncell(dir)) then
        print *,"ncell bust"
        stop
       endif
       if (bfact.eq.bfactC) then
        if (2*ncellC_e(dir).ne.ncell_e(dir)) then
         print *,"ncell_e bust"
         stop
        endif
       else if (bfact.eq.2*bfactC) then
        if (ncellC_e(dir).ne.ncell_e(dir)) then
         print *,"ncell_e bust"
         stop
        endif
       else
        print *,"bfact bust"
        stop
       endif
      enddo ! dir

      allocate(ZCOR(lo(1):hi(1),lo(SDIM):hi(SDIM))) 
      allocate(wtcen(lo(1):hi(1),lo(SDIM):hi(SDIM))) 

      do i1=0,bfactC-1
       yG(i1)=cache_gauss(bfactC,i1)
      enddo
      do i1=0,bfactC
       yGL(i1)=cache_gauss_lobatto(bfactC,i1)
      enddo
      if (bfact.eq.bfactC) then
       call SplitNodes(2*bfact-1,yE)
       call SplitLobattoNodes(2*bfact,yEL)
       allocate(intvolG(0:bfactC-1,0:2*bfact-1))
       call intersect_table(bfactC,2*bfact,yGL,yEL,intvolG)
      else if (bfact.eq.2*bfactC) then
       do i1=0,bfact-1
        yE(i1)=cache_gauss(bfact,i1)
       enddo
       do i1=0,bfact
        yEL(i1)=cache_gauss_lobatto(bfact,i1)
       enddo
       allocate(intvolG(0:bfactC-1,0:bfact-1))
       call intersect_table(bfactC,bfact,yGL,yEL,intvolG)
      else
       print *,"bfact invalid"
       stop
      endif

      do i=lo(1),hi(1)
      do j=lo(2),hi(2)
       wtcen(i,j)=zero
       ZCOR(i,j)=zero
      enddo
      enddo

      do i=loC_e(1),hiC_e(1)
      do j=loC_e(SDIM),hiC_e(SDIM)
       do j1C=0,bfactC-1
       do i1C=0,bfactC-1
        do j1=0,2*bfactC-1
        do i1=0,2*bfactC-1
         dir=-1
         if (bfact.eq.bfactC) then
          if=2*i
          jf=2*j
          i1f=i1
          j1f=j1
          if (i1.ge.bfactC) then
           if=if+1
           i1f=i1-bfactC
          endif
          if (j1.ge.bfactC) then
           jf=jf+1
           j1f=j1-bfactC
          endif
         else if (bfact.eq.2*bfactC) then
          if=i
          jf=j
          i1f=i1
          j1f=j1
         else
          print *,"bfact invalid"
          stop
         endif 

         call abs_index(i,j,i1C,j1C,i2,j2,loC_e,hiC_e, &
            bfactC,loC,hiC,dir)
         localvol=intvolG(i1C,i1)*intvolG(j1C,j1)
         ZSOURCE=localvol*ZC(i2,j2)

         call abs_index(if,jf,i1f,j1f,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
         ZCOR(i2,j2)=ZCOR(i2,j2)+ZSOURCE
         wtcen(i2,j2)=wtcen(i2,j2)+localvol
        enddo 
        enddo 
       enddo 
       enddo 
      enddo
      enddo
      do i=lo(1),hi(1)
      do j=lo(2),hi(2)
       localvol=wtcen(i,j)
       if (localvol.le.zero) then
        print *,"localvol invalid"
        stop
       endif
       ZCOR(i,j)=ZCOR(i,j)/localvol
       ZF(i,j)=ZF(i,j)+ZCOR(i,j)
      enddo
      enddo

      deallocate(ZCOR)
      deallocate(wtcen)
      deallocate(intvolG)

      return
      end subroutine PROLONGATE



      recursive subroutine preconditioner( &
       ncell, &
       lo,hi,bfact,Z,R, &
       betax,betay,dx,hflag,probtype, &
       problo,probhi, &
       nsmooth,precond_type,tol)
      IMPLICIT NONE

      REAL_T tol
      integer lo(SDIM),hi(SDIM),bfact,nsmooth,precond_type
      integer ncell(SDIM)
      REAL_T R(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T Z(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T dx(SDIM),problo(SDIM),probhi(SDIM)
      integer hflag,probtype
      integer bottomflag

      if (hflag.ne.1) then
       print *,"hflag invalid"
       stop
      endif

      call ZAPVEC(lo,hi,Z)
      if (precond_type.eq.0) then
       call COPYVEC(lo,hi,R,Z)
      else if (precond_type.eq.1) then
       call JACPRECOND(lo,hi,Z,R,betax,betay, &
        dx,hflag,probtype,nsmooth,problo,probhi,bfact)
      else if (precond_type.eq.2) then
       call RELAX(lo,hi,Z,R,betax,betay, &
         dx,hflag,problo,probhi,probtype,nsmooth,tol,bfact)
      else
       print *,"precond_type invalid" 
       print *,"precond_type= ",precond_type
       stop
      endif

      return
      end subroutine preconditioner

        ! precond_type=0 M=I, =1 Jacobi, =2 MG
      recursive subroutine bicgstab(ncell,lo,hi,bfact, &
         U,G,betax,betay,betacell,dx,hflag, &
         precond_type,problo,probhi,probtype,tol,piecewise, &
         weak_flag)
      IMPLICIT NONE

      integer bfact,piecewise,weak_flag
      REAL_T problo(SDIM),probhi(SDIM)
      integer i,j,ncell(SDIM),lo(SDIM),hi(SDIM),hflag
      integer precond_type,probtype,iter
      integer maxiter_bicgstab
      integer hflagcg,restart_flag
      REAL_T U(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T G(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T betacell(lo(1):hi(1),lo(SDIM):hi(SDIM))
      REAL_T dx(SDIM),tol,dnorm,dnorm0

      REAL_T, dimension(:,:), allocatable :: U0
      REAL_T, dimension(:,:), allocatable :: V0
      REAL_T, dimension(:,:), allocatable :: P0
      REAL_T, dimension(:,:), allocatable :: R0
      REAL_T, dimension(:,:), allocatable :: U1
      REAL_T, dimension(:,:), allocatable :: V1
      REAL_T, dimension(:,:), allocatable :: P1
      REAL_T, dimension(:,:), allocatable :: R1
      REAL_T, dimension(:,:), allocatable :: R0hat
      REAL_T, dimension(:,:), allocatable :: UINIT
      REAL_T, dimension(:,:), allocatable :: RHS
      REAL_T, dimension(:,:), allocatable :: Y
      REAL_T, dimension(:,:), allocatable :: Hvec
      REAL_T, dimension(:,:), allocatable :: S
      REAL_T, dimension(:,:), allocatable :: T
      REAL_T, dimension(:,:), allocatable :: Z

      REAL_T rho0,w0,rho1,alpha,w1,beta,a1,a2
      integer nsmooth

      if ((weak_flag.lt.0).or.(weak_flag.gt.4)) then
       print *,"weak_flag invalid in bicgstab"
       stop
      endif
      if (precond_type.eq.1) then
       nsmooth=4
      else if (precond_type.eq.2) then
       nsmooth=2
      else if (precond_type.eq.0) then
       nsmooth=0
      else
       print *,"precond_type invalid in bicgstab"
       stop
      endif

      allocate(U0(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(V0(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(P0(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(R0(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(U1(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(V1(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(P1(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(R1(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(R0hat(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(UINIT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(RHS(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(Y(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(Hvec(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(S(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(T(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(Z(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 

        ! U1=U0=V0=P0=0 
      alpha=0.0 
      do i=lo(1)-1,hi(1)+1
      do j=lo(SDIM)-1,hi(SDIM)+1
       U0(i,j)=U(i,j)
       V0(i,j)=alpha
       P0(i,j)=alpha
       U1(i,j)=alpha
      enddo
      enddo

       ! R0=G-A U0
      call RESID(ncell,lo,hi,problo,probhi, &
        R0,U0,G,betax,betay,betacell,dx,bfact,hflag,probtype, &
        piecewise,weak_flag)
       ! R0hat=R0
      call COPYVEC(lo,hi,R0,R0hat)
       ! UINIT=U0
      call COPYVEC(lo,hi,U0,UINIT)
      call ZAPVEC(lo,hi,U0)
       ! RHS=R0
      call COPYVEC(lo,hi,R0,RHS)

       ! rho0=alpha=w0=1
      rho0=1.0
      alpha=1.0
      w0=1.0

      call NORMPROD(lo,hi,R0,dnorm0)
      dnorm=1.0
      if (dnorm0.lt.tol) then
       dnorm=dnorm0
      endif
      maxiter_bicgstab=400
!     maxiter_bicgstab=2
      iter=0

      hflagcg=1
      do while ((dnorm.gt.tol).and.(iter.lt.maxiter_bicgstab))
       print *,"bicgstab bfact,iter,dnorm ",bfact,iter,dnorm

         ! rho1= R0hat^H R0
       call DOTPROD(lo,hi,R0hat,R0,rho1)
       
       restart_flag=0
       if ((sqrt(abs(rho0)).lt.tol*0.01).or. &
           (sqrt(abs(w0)).lt.tol*0.01)) then
        restart_flag=1
       endif 

       if (restart_flag.eq.0) then
          ! (R0hat^H R0)/(R0hat^H dot R0_before)  *   (alpha/w0)
        beta=(rho1/rho0)*(alpha/w0)
        a1=1.0
        a2=-w0
 
         ! P1=P0-w0 V0
        call LINCOMB(lo,hi,P0,V0,P1,a1,a2)
         ! P1=R0+beta P1
        call LINCOMB(lo,hi,R0,P1,P1,a1,beta)
        ! Y=M^{-1}P1
        call preconditioner( &
         ncell, &
         lo,hi,bfact,Y,P1, &
         betax,betay,dx,hflagcg,probtype, &
         problo,probhi, &
         nsmooth,precond_type,tol)
         
         ! V1=A Y
        call ATIMESU(ncell,lo,hi,problo,probhi, &
         dx,bfact,Y,betax,betay,betacell,hflagcg,probtype,V1, &
         piecewise,weak_flag)

         ! alpha=rho1/R0hat dot V1
        call DOTPROD(lo,hi,R0hat,V1,alpha)

        if (sqrt(abs(alpha)).lt.tol*0.01) then
         restart_flag=1
        endif

        if (restart_flag.eq.0) then
         alpha=rho1/alpha

         ! Hvec=U0+alpha Y
         a1=1.0
         a2=alpha
         call LINCOMB(lo,hi,U0,Y,Hvec,a1,a2) 
         ! U1=Hvec
         call COPYVEC(lo,hi,Hvec,U1)
         ! R1=RHS-A U1
         call RESID(ncell,lo,hi,problo,probhi, &
          R1,U1,RHS,betax,betay,betacell,dx,bfact,hflagcg,probtype, &
          piecewise,weak_flag)
         call NORMPROD(lo,hi,R1,dnorm)
         dnorm=dnorm/dnorm0

         if (dnorm.gt.tol) then
          ! S=R0-alpha V1
          a1=1.0
          a2=-alpha
          call LINCOMB(lo,hi,R0,V1,S,a1,a2) 

           ! Z=M^{-1}S
          call preconditioner( &
           ncell, &
           lo,hi,bfact,Z,S, &
           betax,betay,dx,hflagcg,probtype, &
           problo,probhi,nsmooth,precond_type,tol)

           ! T=A Z
          call ATIMESU(ncell,lo,hi,problo,probhi, &
            dx,bfact,Z,betax,betay,betacell,hflagcg,probtype,T, &
            piecewise,weak_flag)

           ! simple case is: (T,S)/(T,T)=(AZ,S)/(AZ,AZ)   (MZ=S)
          call DOTPROD(lo,hi,T,S,a1)
          call DOTPROD(lo,hi,T,T,a2)
          if (sqrt(abs(a2)).lt.tol*0.01) then
           restart_flag=1
          endif

          if (restart_flag.eq.0) then
           w1=a1/a2
           ! U1=Hvec+w1 Z
           a1=1.0
           a2=w1
           call LINCOMB(lo,hi,Hvec,Z,U1,a1,a2) 
          endif
         endif ! dnorm>tol
          ! R1=RHS-A U1
         call RESID(ncell,lo,hi,problo,probhi, &
          R1,U1,RHS,betax,betay,betacell,dx,bfact,hflagcg,probtype, &
          piecewise,weak_flag)
         call NORMPROD(lo,hi,R1,dnorm)
         dnorm=dnorm/dnorm0
         rho0=rho1
         w0=w1
          ! R0=R1
         call COPYVEC(lo,hi,R1,R0) 
         call COPYVEC(lo,hi,P1,P0) 
         call COPYVEC(lo,hi,V1,V0) 
         call COPYVEC(lo,hi,U1,U0) 
        endif  ! restart_flag=0
       endif ! restart_flag=0

       if (restart_flag.eq.0) then
        ! do nothing
       else if (restart_flag.eq.1) then
        print *,"restarting bicgstab"
        call RESID(ncell,lo,hi,problo,probhi, &
          R0,U0,RHS,betax,betay,betacell,dx,bfact,hflagcg,probtype, &
          piecewise,weak_flag)
         ! R0hat=R0
        call COPYVEC(lo,hi,R0,R0hat)
        call COPYVEC(lo,hi,U0,U1) 
         ! rho0=alpha=w0=1
        rho0=1.0
        alpha=1.0
        w0=1.0
        call NORMPROD(lo,hi,R0,dnorm)
        dnorm=dnorm/dnorm0
        call ZAPVEC(lo,hi,V0)
        call ZAPVEC(lo,hi,P0)
       else
        print *,"restart_flag invalid"
        stop
       endif

       iter=iter+1
      enddo
      print *,"bicgstab at the end: bfact,iter,dnorm ",bfact,iter,dnorm
       ! U=UINIT+U1
      a1=1.0
      a2=a1
      call LINCOMB(lo,hi,UINIT,U1,U,a1,a2)
 
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
      end subroutine bicgstab


        ! precond_type=0 M=I, =1 Jacobi, =2 MG
      recursive subroutine pcg(ncell,lo,hi,bfact, &
         U,G,betax,betay,betacell,dx,hflag, &
         precond_type,problo,probhi,probtype,tol, &
         piecewise,weak_flag,bottomflag)
      IMPLICIT NONE

      integer bfact,bottomflag,piecewise,weak_flag
      REAL_T problo(SDIM),probhi(SDIM)
      integer i,j,ncell(SDIM),lo(SDIM),hi(SDIM),hflag
      integer precond_type,probtype,iter
      integer maxiter_pcg
      integer hflagcg,restart_flag
      REAL_T U(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T G(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T betacell(lo(1):hi(1),lo(SDIM):hi(SDIM))
      REAL_T dx(SDIM),tol,dnorm,dnorm0

      REAL_T, dimension(:,:), allocatable :: U0
      REAL_T, dimension(:,:), allocatable :: P0
      REAL_T, dimension(:,:), allocatable :: R0
      REAL_T, dimension(:,:), allocatable :: U1
      REAL_T, dimension(:,:), allocatable :: R1
      REAL_T, dimension(:,:), allocatable :: UINIT
      REAL_T, dimension(:,:), allocatable :: RHS
      REAL_T, dimension(:,:), allocatable :: S
      REAL_T, dimension(:,:), allocatable :: T

      REAL_T rho0,rho1,alpha,w1,beta,a1,a2
      integer nsmooth

      if (bottomflag.eq.0) then
       if ((weak_flag.ne.4).and.(piecewise.ne.1)) then
        print *,"weak_flag or piecewise invalid"
        stop
       endif
      else if (bottomflag.eq.1) then
       if ((weak_flag.ne.0).or.(piecewise.ne.1)) then
        print *,"weak_flag or piecewise invalid"
        stop
       endif
      else
       print *,"bottomflag invalid"
       stop
      endif

      if (precond_type.eq.1) then
       nsmooth=4
      else if (precond_type.eq.2) then
       nsmooth=2
      else if (precond_type.eq.0) then
       nsmooth=0
      else
       print *,"precond_type invalid in pcg"
       stop
      endif

      allocate(U0(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(P0(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(R0(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(U1(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(R1(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(UINIT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(RHS(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(S(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(T(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 

      alpha=0.0 
      do i=lo(1)-1,hi(1)+1
      do j=lo(SDIM)-1,hi(SDIM)+1
       U0(i,j)=U(i,j)
       P0(i,j)=alpha
       U1(i,j)=alpha
      enddo
      enddo

       ! R0=G-A U0
      call RESID(ncell,lo,hi,problo,probhi, &
        R0,U0,G,betax,betay,betacell,dx,bfact,hflag,probtype, &
        piecewise,weak_flag)
       ! UINIT=U0
      call COPYVEC(lo,hi,U0,UINIT)
      call ZAPVEC(lo,hi,U0)
       ! RHS=R0
      call COPYVEC(lo,hi,R0,RHS)
      call COPYVEC(lo,hi,R0,R1)

      hflagcg=1

        ! P0=M^{-1}R0
      call preconditioner( &
         ncell, &
         lo,hi,bfact,P0,R0, &
         betax,betay,dx,hflagcg,probtype, &
         problo,probhi, &
         nsmooth,precond_type,tol)

       ! rho1= R0^H P0
      call DOTPROD(lo,hi,R0,P0,rho1)
      rho0=rho1

      alpha=1.0

      call NORMPROD(lo,hi,R0,dnorm0)
      dnorm=1.0
      if (dnorm0.lt.tol) then
       dnorm=dnorm0
      endif
      maxiter_pcg=400
!     maxiter_pcg=2
      iter=0

      do while ((dnorm.gt.tol).and.(iter.lt.maxiter_pcg))
       if (bottomflag.eq.0) then
        print *,"pcg iter,dnorm ",iter,dnorm
       else
        print *,"pcg (bottom) iter,dnorm ",iter,dnorm
       endif

         ! T=A P0
       call ATIMESU(ncell,lo,hi,problo,probhi, &
         dx,bfact,P0,betax,betay,betacell,hflagcg,probtype,T, &
         piecewise,weak_flag)

         ! alpha= P0^H T
       call DOTPROD(lo,hi,P0,T,alpha)
       
       restart_flag=0
       if (sqrt(abs(alpha)).lt.tol*0.000001) then
        restart_flag=1
       endif 

       if (restart_flag.eq.0) then
        alpha=rho1/alpha

         ! U1=U0+alpha P0
        a1=1.0
        a2=alpha
        call LINCOMB(lo,hi,U0,P0,U1,a1,a2) 

         ! R1=RHS-A U1
        call RESID(ncell,lo,hi,problo,probhi, &
         R1,U1,RHS,betax,betay,betacell,dx,bfact,hflagcg,probtype, &
         piecewise,weak_flag)

         ! S=M^{-1}R1
        call preconditioner( &
         ncell, &
         lo,hi,bfact,S,R1, &
         betax,betay,dx,hflagcg,probtype, &
         problo,probhi, &
         nsmooth,precond_type,tol)

        rho0=rho1

         ! rho1=R1 dot S
        call DOTPROD(lo,hi,R1,S,rho1)

        beta=rho1/rho0 

        a1=1.0
        a2=beta
 
         ! P0=S+beta P0
        call LINCOMB(lo,hi,S,P0,P0,a1,a2)
       endif ! restart_flag=0

       if (restart_flag.eq.0) then
        call COPYVEC(lo,hi,R1,R0) 
        call COPYVEC(lo,hi,U1,U0) 
        call NORMPROD(lo,hi,R0,dnorm)
       else if (restart_flag.eq.1) then
        print *,"restarting pcg"
        call RESID(ncell,lo,hi,problo,probhi, &
          R0,U0,RHS,betax,betay,betacell, &
          dx,bfact,hflagcg,probtype, &
          piecewise,weak_flag)
        call COPYVEC(lo,hi,U0,U1) 

         ! P0=M^{-1}R0
        call preconditioner( &
         ncell, &
         lo,hi,bfact,P0,R0, &
         betax,betay,dx,hflagcg,probtype, &
         problo,probhi, &
         nsmooth,precond_type,tol)

         ! rho1= R0^H P0
        call DOTPROD(lo,hi,R0,P0,rho1)
        rho0=rho1

        alpha=1.0

        call NORMPROD(lo,hi,R0,dnorm)
       else
        print *,"restart_flag invalid"
        stop
       endif

       iter=iter+1
      enddo
      if (bottomflag.eq.0) then
       print *,"pcg at the end: iter,dnorm ",iter,dnorm
      else
       print *,"pcg (bottom) at the end: iter,dnorm ",iter,dnorm
      endif
       ! U=UINIT+U1
      a1=1.0
      a2=a1
      call LINCOMB(lo,hi,UINIT,U1,U,a1,a2)
 
      deallocate(U0) 
      deallocate(P0) 
      deallocate(R0) 
      deallocate(U1) 
      deallocate(R1) 
      deallocate(UINIT) 
      deallocate(RHS) 
      deallocate(S) 
      deallocate(T) 

      return
      end subroutine pcg




      subroutine SEM_DIV(ncell,lo,hi,problo,probhi,dx,bfact, &
        P,U,V,betax,betay,ubcLEFT,ubcRIGHT, &
        G,hflag,probtype,for_source,piecewise,weak_flag)
      use LegendreNodes
      use Nodal2DClass
      IMPLICIT NONE

      integer for_source,piecewise,weak_flag
      integer lo_e(SDIM),hi_e(SDIM)
      integer ncell(SDIM),lo(SDIM),hi(SDIM),i,j
      REAL_T problo(SDIM),probhi(SDIM),dx(SDIM)
      integer bfact,hflag,probtype,dir
      REAL_T U(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T V(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T ubcLEFT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1,SDIM)
      REAL_T ubcRIGHT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1,SDIM)
      REAL_T G(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T P(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T PLINE(0:bfact)
      REAL_T ULINE(0:bfact)
      REAL_T VLINE(0:bfact)
      REAL_T DU(0:bfact-1)
      REAL_T y(0:bfact)
      REAL_T yGL(0:bfact)
      REAL_T wMAT(0:bfact,0:bfact)
      REAL_T dx_element
      integer i1,j1,i2,j2
      REAL_T AX,AY,vol
      REAL_T dxL(SDIM)
      REAL_T dxR(SDIM)
      REAL_T, dimension(:,:), allocatable :: IuDwCELL
      REAL_T, dimension(:,:,:), allocatable :: qbcLEFT
      REAL_T, dimension(:,:,:), allocatable :: qbcRIGHT

      REAL_T qbc(0:1)
      REAL_T IuDw(0:bfact)
      REAL_T xlo_e(SDIM),xhi_e(SDIM)
      integer dirbc,bctype,side
      REAL_T bc_sign
      REAL_T QLEFT,QRIGHT,FLUXLEFT,FLUXRIGHT
      REAL_T UJUMPLEFT,UJUMPRIGHT,tau

      
      TYPE (DGprecomputed) spA
      CHARACTER(LEN=2) :: typenodes  

      typenodes="Ga"
 
      if ((weak_flag.lt.0).or.(weak_flag.gt.4)) then
       print *,"weak_flag invalid in SEM_DIV"
       stop
      endif 
      do dir=1,SDIM
       lo_e(dir)=lo(dir)/bfact
       hi_e(dir)=(hi(dir)+1)/bfact-1
      enddo

      if ((weak_flag.eq.0).or.(weak_flag.eq.4)) then
 
       if (probtype.eq.0) then
        ! left BC
        i=lo(1)
        do j=lo(SDIM),hi(SDIM)
         if (for_source.eq.1) then
          U(i,j)=1.0
         else
          U(i,j)=0.0
         endif
        enddo
        ! top and bottom BC
        do i=lo(1),hi(1)
         j=lo(SDIM)
         V(i,j)=0.0
         j=hi(SDIM)+1
         V(i,j)=0.0
        enddo
       else if (probtype.eq.1) then
        ! do nothing
       else if (probtype.eq.2) then
        ! left BC
        i=lo(1)
        do j=lo(SDIM),hi(SDIM)
         if (for_source.eq.1) then
          U(i,j)=1.0
         else
          U(i,j)=0.0
         endif
        enddo
        ! top and bottom BC
        do i=lo(1),hi(1)
         j=lo(SDIM)
         V(i,j)=0.0
         j=hi(SDIM)+1
         V(i,j)=0.0
        enddo
       else if (probtype.eq.3) then
        ! do nothing
       else
        print *,"probtype invalid"
        stop
       endif

      else if ((weak_flag.ge.1).and.(weak_flag.le.3)) then
       ! do nothing
      else
       print *,"weak_flag invalid in SEM_DIV weak_flag=",weak_flag
       stop
      endif

      do i1=0,bfact-1
       y(i1)=cache_gauss(bfact,i1)
      enddo
      do i1=0,bfact
       yGL(i1)=cache_gauss_lobatto(bfact,i1)
      enddo
      call polyinterp_Dmatrix(bfact,yGL,wMAT)

      if (piecewise.eq.0) then

       if ((weak_flag.eq.0).or.(weak_flag.eq.4)) then

        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

          ! G holds the divergence and is at Gauss points.
         do i1=0,bfact-1
         do j1=0,bfact-1
          dir=-1
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
          G(i2,j2)=0.0
         enddo
         enddo
         do j1=0,bfact-1
          do i1=0,bfact
           dir=0
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           ULINE(i1)=U(i2,j2)
          enddo
           ! ULINE: Legendre gauss lobatto
           ! DU   : Legendre gauss
          dx_element=bfact*dx(1)
          call deriv_change_basis(bfact,bfact-1,ULINE,DU, &
           wMAT,yGL,y,dx_element)
 
          do i1=0,bfact-1
           dir=-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           G(i2,j2)=G(i2,j2)+DU(i1)
          enddo
         enddo ! j1
 
         do i1=0,bfact-1
          do j1=0,bfact
           dir=1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           ULINE(j1)=V(i2,j2)
          enddo
           ! ULINE: Legendre gauss lobatto
           ! DU   : Legendre gauss
          dx_element=bfact*dx(2)
          call deriv_change_basis(bfact,bfact-1,ULINE,DU, &
             wMAT,yGL,y,dx_element)

          do j1=0,bfact-1
           dir=-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           G(i2,j2)=G(i2,j2)+DU(j1)
          enddo
         enddo ! i1
        
        enddo
        enddo

       else if ((weak_flag.ge.1).and.(weak_flag.le.3)) then


        call ConstructNodal2D(bfact-1,spA,typenodes)
        allocate(qbcLEFT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1,SDIM)) 
        allocate(qbcRIGHT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1,SDIM)) 
        allocate(IuDwCELL(lo(1):hi(1),lo(SDIM):hi(SDIM))) 



         ! ux 
        dir=1
        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

         call element_size(problo,i,j,lo_e,hi_e,bfact, &
           lo,hi,dx,xlo_e,xhi_e)
         do j1=0,bfact-1
          do i1=0,bfact-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           ULINE(i1)=U(i2,j2)
           PLINE(i1)=P(i2,j2)
          enddo

          call GradientWeakFormCompute(spA,ULINE,IuDw,qbc)

          if ((weak_flag.eq.3).and.(for_source.eq.0)) then
           call InternalPenalty_GradientFlux(spA,PLINE, &
             xlo_e(dir),xhi_e(dir),qbc) 
           i1=0
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           qbc(0)=qbc(0)*betax(i2,j2)
           i1=bfact
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           qbc(1)=qbc(1)*betax(i2,j2)
          endif

          do i1=0,bfact-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           IuDwCELL(i2,j2)=IuDw(i1)
          enddo

          i1=0
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
          qbcLEFT(i,j2,dir)=qbc(0)
          qbcRIGHT(i,j2,dir)=qbc(1)
          if (i.eq.lo_e(dir)) then
           dirbc=dir
           side=1
           call get_bctype(dirbc,side,bctype,bc_sign,probtype)
           if (bctype.eq.0) then ! dirichlet
            qbcRIGHT(i-1,j2,dir)=qbc(0)
            if (for_source.eq.1) then
             if ((probtype.eq.0).or.(probtype.eq.2)) then
              qbcRIGHT(i-1,j2,dir)=one
             endif
            endif
           else
            print *,"bctype invalid"
            stop
           endif
          endif
          if (i.eq.hi_e(dir)) then
           dirbc=dir
           side=2
           call get_bctype(dirbc,side,bctype,bc_sign,probtype)
           if (bctype.eq.0) then ! dirichlet
            qbcLEFT(i+1,j2,dir)=qbc(1)
           else if (bctype.eq.1) then  ! neumann
            qbcLEFT(i+1,j2,dir)=zero
           else
            print *,"bctype invalid"
            stop
           endif
          endif

         enddo !j1

        enddo
        enddo

        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

         call element_size(problo,i,j,lo_e,hi_e,bfact, &
           lo,hi,dx,xlo_e,xhi_e)
         do j1=0,bfact-1
          i1=0
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)

          if (weak_flag.eq.1) then
           QLEFT=qbcLEFT(i,j2,dir)
           QRIGHT=qbcLEFT(i+1,j2,dir)
           if (i.eq.lo_e(dir)) then
            QLEFT=qbcRIGHT(i-1,j2,dir)
           endif
          else if ((weak_flag.eq.2).or.(weak_flag.eq.3)) then
           if (i.eq.lo_e(dir)) then
            QLEFT=qbcRIGHT(i-1,j2,dir)
           else
            QLEFT=half*(qbcLEFT(i,j2,dir)+qbcRIGHT(i-1,j2,dir))
           endif

           if (i.eq.hi_e(dir)) then
            QRIGHT=qbcLEFT(i+1,j2,dir)
           else
            QRIGHT=half*(qbcLEFT(i+1,j2,dir)+qbcRIGHT(i,j2,dir))
           endif
          else
           print *,"weak_flag invalid"
           stop
          endif

          if (for_source.eq.1) then
           UJUMPLEFT=zero
           UJUMPRIGHT=zero
          else if (for_source.eq.0) then
           UJUMPLEFT=ubcRIGHT(i-1,j2,dir)-ubcLEFT(i,j2,dir)
           UJUMPRIGHT=ubcRIGHT(i,j2,dir)-ubcLEFT(i+1,j2,dir)
           if ((weak_flag.eq.1).or.(weak_flag.eq.2)) then
            tau=bfact/(xhi_e(dir)-xlo_e(dir))
           else if (weak_flag.eq.3) then
            tau=bfact*bfact/(xhi_e(dir)-xlo_e(dir))
           else
            print *,"weak_flag invalid"
            stop
           endif
           UJUMPLEFT=tau*UJUMPLEFT
           UJUMPRIGHT=tau*UJUMPRIGHT
          else
           print *,"for_source invalid"
           stop
          endif

          do i1=0,bfact-1
           FLUXLEFT=(QLEFT-UJUMPLEFT)*spA%Lminus1(i1)
           FLUXRIGHT=(QRIGHT-UJUMPRIGHT)*spA%L1(i1)
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           G(i2,j2)=(IuDwCELL(i2,j2)+(FLUXRIGHT-FLUXLEFT))/ &
             (half*(xhi_e(dir)-xlo_e(dir))*spA%Qweights(i1))
          enddo

         enddo !j1
        enddo
        enddo
          

         ! vy
        dir=2
        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

         call element_size(problo,i,j,lo_e,hi_e,bfact, &
           lo,hi,dx,xlo_e,xhi_e)
         do i1=0,bfact-1
          do j1=0,bfact-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           VLINE(j1)=V(i2,j2)
           PLINE(j1)=P(i2,j2)
          enddo

          call GradientWeakFormCompute(spA,VLINE,IuDw,qbc)

          if ((weak_flag.eq.3).and.(for_source.eq.0)) then
           call InternalPenalty_GradientFlux(spA,PLINE, &
             xlo_e(dir),xhi_e(dir),qbc) 
           j1=0
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           qbc(0)=qbc(0)*betay(i2,j2)
           j1=bfact
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           qbc(1)=qbc(1)*betay(i2,j2)
          endif

          do j1=0,bfact-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           IuDwCELL(i2,j2)=IuDw(j1)
          enddo

          j1=0
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
          qbcLEFT(i2,j,dir)=qbc(0)
          qbcRIGHT(i2,j,dir)=qbc(1)
          if (j.eq.lo_e(dir)) then
           dirbc=dir
           side=1
           call get_bctype(dirbc,side,bctype,bc_sign,probtype)
           if (bctype.eq.0) then ! dirichlet
            qbcRIGHT(i2,j-1,dir)=qbc(0)
           else if (bctype.eq.1) then ! neumann
            qbcRIGHT(i2,j-1,dir)=zero
           else
            print *,"bctype invalid"
            stop
           endif
          endif
          if (j.eq.hi_e(dir)) then
           dirbc=dir
           side=2
           call get_bctype(dirbc,side,bctype,bc_sign,probtype)
           if (bctype.eq.0) then ! dirichlet
            qbcLEFT(i2,j+1,dir)=qbc(1)
           else if (bctype.eq.1) then  ! neumann
            qbcLEFT(i2,j+1,dir)=zero
           else
            print *,"bctype invalid"
            stop
           endif
          endif

         enddo !i1

        enddo
        enddo

        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

         call element_size(problo,i,j,lo_e,hi_e,bfact, &
           lo,hi,dx,xlo_e,xhi_e)
         do i1=0,bfact-1
          j1=0
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)

          if (weak_flag.eq.1) then

           QLEFT=qbcLEFT(i2,j,dir)
           QRIGHT=qbcLEFT(i2,j+1,dir)
           if (j.eq.lo_e(dir)) then
            QLEFT=qbcRIGHT(i2,j-1,dir)
           endif

          else if ((weak_flag.eq.2).or.(weak_flag.eq.3)) then
           if (i.eq.lo_e(dir)) then
            QLEFT=qbcRIGHT(i2,j-1,dir)
           else
            QLEFT=half*(qbcLEFT(i2,j,dir)+qbcRIGHT(i2,j-1,dir))
           endif

           if (i.eq.hi_e(dir)) then
            QRIGHT=qbcLEFT(i2,j+1,dir)
           else
            QRIGHT=half*(qbcLEFT(i2,j+1,dir)+qbcRIGHT(i2,j,dir))
           endif
          else
           print *,"weak_flag invalid"
           stop
          endif



          if (for_source.eq.1) then
           UJUMPLEFT=zero
           UJUMPRIGHT=zero
          else if (for_source.eq.0) then
           UJUMPLEFT=ubcRIGHT(i2,j-1,dir)-ubcLEFT(i2,j,dir)
           UJUMPRIGHT=ubcRIGHT(i2,j,dir)-ubcLEFT(i2,j+1,dir)
           if ((weak_flag.eq.1).or.(weak_flag.eq.2)) then
            tau=bfact/(xhi_e(dir)-xlo_e(dir))
           else if (weak_flag.eq.3) then
            tau=bfact*bfact/(xhi_e(dir)-xlo_e(dir))
           else
            print *,"weak_flag invalid"
            stop
           endif
           UJUMPLEFT=tau*UJUMPLEFT
           UJUMPRIGHT=tau*UJUMPRIGHT
          else
           print *,"for_source invalid"
           stop
          endif

          do j1=0,bfact-1
           FLUXLEFT=(QLEFT-UJUMPLEFT)*spA%Lminus1(j1)
           FLUXRIGHT=(QRIGHT-UJUMPRIGHT)*spA%L1(j1)
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           G(i2,j2)=G(i2,j2)+(IuDwCELL(i2,j2)+(FLUXRIGHT-FLUXLEFT))/ &
             (half*(xhi_e(dir)-xlo_e(dir))*spA%Qweights(j1))
          enddo

         enddo !i1
        enddo
        enddo


        deallocate(qbcLEFT)
        deallocate(qbcRIGHT)
        deallocate(IuDwCELL)
        call DestructNodal2D(spA)
       else
        print *,"weak_flag invalid"
        stop
       endif

      else if (piecewise.eq.1) then
       do i=lo(1),hi(1)
       do j=lo(SDIM),hi(SDIM)
        call FVgeom(problo,i,j,lo,bfact,dx,AX,AY,dxL,dxR,vol)
        if (bfact.eq.1) then
         if (abs(vol/AX-dx(1)).gt.1.0E-13) then
          print *,"FVgeom failed"
          stop
         endif
         if (abs(vol/AY-dx(2)).gt.1.0E-13) then
          print *,"FVgeom failed"
          stop
         endif
        endif
        G(i,j)=(AX*(U(i+1,j)-U(i,j))+AY*(V(i,j+1)-V(i,j)))/vol
       enddo
       enddo
      else
       print *,"piecewise invalid"
       stop
      endif

      return
      end subroutine SEM_DIV


      subroutine TFRCELL(lo,hi,problo,probhi,dx,bfact, &
        U,UCELL,V,VCELL,piecewise,weak_flag)
      use LegendreNodes
      IMPLICIT NONE

      integer weak_flag
      integer lo_e(SDIM),hi_e(SDIM)
      integer lo(SDIM),hi(SDIM),i,j,piecewise
      REAL_T problo(SDIM),probhi(SDIM),dx(SDIM)
      integer bfact,dir
      REAL_T U(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T V(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T UCELL(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T VCELL(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T ULINE(0:bfact)
      REAL_T DU(0:bfact-1)
      REAL_T y(0:bfact)
      REAL_T yGL(0:bfact)
      integer i1,j1,i2,j2
      REAL_T theta
      REAL_T x(SDIM)
      REAL_T xL(SDIM)
      REAL_T xR(SDIM)
 
      do dir=1,SDIM
       lo_e(dir)=lo(dir)/bfact
       hi_e(dir)=(hi(dir)+1)/bfact-1
      enddo
 
      do i1=0,bfact-1
       y(i1)=cache_gauss(bfact,i1)
      enddo
      do i1=0,bfact
       yGL(i1)=cache_gauss_lobatto(bfact,i1)
      enddo

      if (piecewise.eq.0) then

       if ((weak_flag.eq.0).or.(weak_flag.eq.4)) then

        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)
         do j1=0,bfact-1
          do i1=0,bfact
           dir=0
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           ULINE(i1)=U(i2,j2)
          enddo
          call poly_change_basis(bfact,bfact-1,ULINE,DU,yGL,y)
          do i1=0,bfact-1
           dir=-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           UCELL(i2,j2)=DU(i1)
          enddo
         enddo
 
         do i1=0,bfact-1
          do j1=0,bfact
           dir=1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           ULINE(j1)=V(i2,j2)
          enddo
          call poly_change_basis(bfact,bfact-1,ULINE,DU,yGL,y)
          do j1=0,bfact-1
           dir=-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           VCELL(i2,j2)=DU(j1)
          enddo
         enddo
        
        enddo
        enddo

       else if ((weak_flag.ge.1).and.(weak_flag.le.3)) then

        do i=lo(1),hi(1)
        do j=lo(SDIM),hi(SDIM)
         UCELL(i,j)=U(i,j) 
         VCELL(i,j)=V(i,j) 
        enddo
        enddo

       else 
        print *,"weak_flag invalid"
        stop
       endif

      else if (piecewise.eq.1) then
       do i=lo(1),hi(1)
       do j=lo(SDIM),hi(SDIM)
        call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1)
        dir=1
        theta=(x(dir)-xL(dir))/(xR(dir)-xL(dir))
        UCELL(i,j)=(one-theta)*U(i,j)+theta*U(i+1,j)
        dir=2
        theta=(x(dir)-xL(dir))/(xR(dir)-xL(dir))
        VCELL(i,j)=(one-theta)*V(i,j)+theta*V(i,j+1)
       enddo
       enddo

      else
       print *,"piecewise invalid"
       stop
      endif

      return
      end subroutine TFRCELL

      subroutine SEM_Neumann(PIN,BCVAL,POUT,y,bfact,piecewise, &
         weak_flag)
      use LegendreNodes
      IMPLICIT NONE

      integer bfact,i,piecewise,weak_flag
      REAL_T PIN(0:bfact-1)
      REAL_T BCVAL
      REAL_T POUT,sum
      REAL_T y(0:bfact-1)
      REAL_T yextend(0:bfact)
      REAL_T pextend(0:bfact)
      REAL_T wextend(0:bfact)
      REAL_T w(0:bfact,0:bfact)
      REAL_T xtarget

      if (BCVAL.ne.zero) then
       print *,"inhomogeneous Neumann BC not ready yet (need dx)"
       stop
      endif

      if (piecewise.eq.0) then

       if ((weak_flag.ge.0).and.(weak_flag.le.4)) then
        yextend(0)=-one
        do i=1,bfact
         yextend(i)=y(i-1)
         pextend(i)=PIN(i-1)
        enddo
        call polyinterp_Dmatrix(bfact,yextend,w)
        sum=zero
        do i=1,bfact
         sum=sum+PIN(i-1)*w(i,0)
        enddo
        if (w(0,0).eq.zero) then
         print *,"Dmatrix corruption"
         stop
        endif
        pextend(0)=(BCVAL-sum)/w(0,0)
        if ((weak_flag.eq.0).or. &
            ((weak_flag.ge.1).and.(weak_flag.le.3))) then
         xtarget=-one-abs(yextend(1)+one)
        else if (weak_flag.eq.4) then
         xtarget=-one
        else
         print *,"weak_flag invalid"
         stop
        endif
        call polyinterp_weights(bfact,yextend,wextend,xtarget)
        POUT=zero
        do i=0,bfact
         POUT=POUT+pextend(i)*wextend(i)
        enddo
 
        if (bfact.eq.1) then
         if ((weak_flag.eq.0).or. &
             ((weak_flag.ge.1).and.(weak_flag.le.3))) then
          if (abs(half*(PIN(0)-POUT)-BCVAL).gt.1.0E-13) then
           print *,"SEM_Neumann corruption"
           stop
          endif
         else if (weak_flag.eq.4) then
          if (abs((PIN(0)-POUT)-BCVAL).gt.1.0E-13) then
           print *,"SEM_Neumann corruption"
           stop
          endif
         else
          print *,"weak_flag invalid"
          stop
         endif
        endif
       else
        print *,"weak_flag invalid SEM_Neumann",weak_flag
        stop
       endif

      else if (piecewise.eq.1) then
       POUT=PIN(0)-two*BCVAL
      else
       print *,"piecewise invalid"
       stop
      endif

      return
      end subroutine SEM_Neumann


      subroutine SEM_Dirichlet(PIN,BCVAL,POUT,y,bfact,piecewise, &
        weak_flag)
      use LegendreNodes
      IMPLICIT NONE

      integer bfact,i,piecewise,weak_flag
      REAL_T PIN(0:bfact-1)
      REAL_T BCVAL
      REAL_T POUT
      REAL_T y(0:bfact-1)
      REAL_T yextend(0:bfact)
      REAL_T pextend(0:bfact)
      REAL_T wextend(0:bfact)
      REAL_T xtarget

      if (piecewise.eq.0) then

       if ((weak_flag.eq.0).or.(weak_flag.eq.4)) then
        yextend(0)=-one
        do i=1,bfact
         yextend(i)=y(i-1)
         pextend(i)=PIN(i-1)
        enddo
        pextend(0)=BCVAL
        if (weak_flag.eq.0) then
         xtarget=-one-abs(yextend(1)+one)
        else if (weak_flag.eq.4) then
         xtarget=-one
        else
         print *,"weak_flag invalid"
         stop
        endif
        call polyinterp_weights(bfact,yextend,wextend,xtarget)
        POUT=zero
        do i=0,bfact
         POUT=POUT+pextend(i)*wextend(i)
        enddo

        if (bfact.eq.1) then
         if (weak_flag.eq.0) then
          if (abs(half*(PIN(0)+POUT)-BCVAL).gt.1.0E-13) then
           print *,"SEM_Dirichlet corruption"
           stop
          endif
         else if (weak_flag.eq.4) then
          if (abs(POUT-BCVAL).gt.1.0E-13) then
           print *,"SEM_Dirichlet corruption"
           stop
          endif
         else
          print *,"weak_flag invalid"
          stop
         endif
        endif
       else if ((weak_flag.ge.1).and.(weak_flag.le.3)) then
        POUT=BCVAL
       else
        print *,"weak_flag invalid SEM_Dirichlet"
        stop
       endif

      else if (piecewise.eq.1) then
       POUT=two*BCVAL-PIN(0)
      else
       print *,"piecewise invalid"
       stop
      endif

      return
      end subroutine SEM_Dirichlet


      subroutine set_boundary(lo,hi,dx,hflag,probtype,P,bfact, &
        problo,probhi,piecewise,weak_flag)
      use LegendreNodes
      IMPLICIT NONE

      integer piecewise,weak_flag
      REAL_T problo(SDIM),probhi(SDIM)
      integer lo(SDIM),hi(SDIM),i,j,hflag,probtype,bfact
      REAL_T dx(SDIM)
      REAL_T x(SDIM)
      REAL_T xL(SDIM)
      REAL_T xR(SDIM)
      REAL_T P(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      integer i1,j1
      REAL_T PIN(0:bfact-1)
      REAL_T y(0:bfact-1)
      REAL_T BCVAL

      if ((weak_flag.lt.0).or.(weak_flag.gt.4)) then
       print *,"weak_flag invalid set_boundary",weak_flag
       stop
      endif
      do i1=0,bfact-1
       y(i1)=cache_gauss(bfact,i1)
      enddo

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      BCVAL=zero
       ! left and right BC
      do j=lo(SDIM),hi(SDIM)
       i=lo(1)
       do i1=0,bfact-1
        PIN(i1)=P(i+i1,j)
       enddo
       if (probtype.eq.0) then
        BCVAL=zero
        call SEM_Neumann(PIN,BCVAL,P(i-1,j),y,bfact,piecewise,weak_flag)
       else if (probtype.eq.1) then
        BCVAL=zero
        if (hflag.eq.0) then
         call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1) 
         xL(2)=x(2)
         call get_PEXACT(xL,dx,BCVAL,probtype,bfact,weak_flag)
        endif
        call SEM_Dirichlet(PIN,BCVAL,P(i-1,j),y, &
         bfact,piecewise,weak_flag)
       else if (probtype.eq.2) then
        BCVAL=zero
        call SEM_Neumann(PIN,BCVAL,P(i-1,j),y,bfact,piecewise,weak_flag)
       else if (probtype.eq.3) then
        BCVAL=zero
        call SEM_Dirichlet(PIN,BCVAL,P(i-1,j),y, &
         bfact,piecewise,weak_flag)
       else
        print *,"probtype invalid"
        stop
       endif
       i=hi(1)
       do i1=0,bfact-1
        PIN(i1)=P(i-i1,j)
       enddo
       if (probtype.eq.0) then
        BCVAL=zero
        call SEM_Dirichlet(PIN,BCVAL,P(i+1,j),y, &
         bfact,piecewise,weak_flag)
       else if (probtype.eq.1) then
        BCVAL=zero
        if (hflag.eq.0) then
         call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1)
         xR(2)=x(2)
         call get_PEXACT(xR,dx,BCVAL,probtype,bfact,weak_flag)
        endif
        call SEM_Dirichlet(PIN,BCVAL,P(i+1,j),y, &
         bfact,piecewise,weak_flag)
       else if (probtype.eq.2) then
        BCVAL=zero
        call SEM_Dirichlet(PIN,BCVAL,P(i+1,j),y, &
         bfact,piecewise,weak_flag)
       else if (probtype.eq.3) then
        BCVAL=zero
        call SEM_Dirichlet(PIN,BCVAL,P(i+1,j),y, &
         bfact,piecewise,weak_flag)
       else
        print *,"probtype invalid"
        stop
       endif
      enddo ! j

       ! bottom and top BC
      do i=lo(1),hi(1)
       j=lo(SDIM)
       do j1=0,bfact-1
        PIN(j1)=P(i,j+j1)
       enddo
       if (probtype.eq.0) then
        BCVAL=zero
        call SEM_Neumann(PIN,BCVAL,P(i,j-1),y,bfact,piecewise,weak_flag)
       else if (probtype.eq.1) then
        BCVAL=zero
        if (hflag.eq.0) then
         call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1)
         xL(1)=x(1)
         call get_PEXACT(xL,dx,BCVAL,probtype,bfact,weak_flag)
        endif
        call SEM_Dirichlet(PIN,BCVAL,P(i,j-1),y, &
         bfact,piecewise,weak_flag)
       else if (probtype.eq.2) then
        BCVAL=zero
        call SEM_Neumann(PIN,BCVAL,P(i,j-1),y, &
          bfact,piecewise,weak_flag)
       else if (probtype.eq.3) then
        BCVAL=zero
        call SEM_Dirichlet(PIN,BCVAL,P(i,j-1),y, &
          bfact,piecewise,weak_flag)
       else
        print *,"probtype invalid"
        stop
       endif
       j=hi(SDIM) 
       do j1=0,bfact-1
        PIN(j1)=P(i,j-j1)
       enddo
       if (probtype.eq.0) then
        BCVAL=zero
        call SEM_Neumann(PIN,BCVAL,P(i,j+1),y, &
          bfact,piecewise,weak_flag)
       else if (probtype.eq.1) then
        BCVAL=zero
        if (hflag.eq.0) then
         call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1)
         xR(1)=x(1)
         call get_PEXACT(xR,dx,BCVAL,probtype,bfact,weak_flag)
        endif
        call SEM_Dirichlet(PIN,BCVAL,P(i,j+1), &
          y,bfact,piecewise,weak_flag)
       else if (probtype.eq.2) then
        BCVAL=zero
        call SEM_Neumann(PIN,BCVAL,P(i,j+1), &
          y,bfact,piecewise,weak_flag)
       else if (probtype.eq.3) then
        BCVAL=zero
        call SEM_Dirichlet(PIN,BCVAL,P(i,j+1), &
         y,bfact,piecewise,weak_flag)
       else
        print *,"probtype invalid"
        stop
       endif
      enddo ! i

       ! corners
      P(lo(1)-1,lo(2)-1)=P(lo(1),lo(2))
      P(lo(1)-1,hi(2)+1)=P(lo(1),hi(2))
      P(hi(1)+1,lo(2)-1)=P(hi(1),lo(2))
      P(hi(1)+1,hi(2)+1)=P(hi(1),hi(2))
      return
      end subroutine set_boundary


      subroutine element_size(problo,i,j,lo_e,hi_e,bfact, &
         lo,hi,dx,xlo_e,xhi_e)
      IMPLICIT NONE

      integer i,j,bfact,dir,i1,j1,index(SDIM)
      integer lo_e(SDIM),hi_e(SDIM)
      integer lo(SDIM),hi(SDIM)
      REAL_T xlo_e(SDIM),xhi_e(SDIM),dx(SDIM),problo(SDIM)

      do dir=0,SDIM-1
       i1=0
       j1=0
       call abs_index(i,j,i1,j1,index(1),index(2), &
         lo_e,hi_e,bfact,lo,hi,dir)
       call gridloc1D(xlo_e(dir+1),problo(dir+1), &
         index(dir+1),lo(dir+1),bfact,dx(dir+1),1)
       i1=bfact
       j1=bfact
       call abs_index(i,j,i1,j1,index(1),index(2), &
         lo_e,hi_e,bfact,lo,hi,dir)
       call gridloc1D(xhi_e(dir+1),problo(dir+1), &
         index(dir+1),lo(dir+1),bfact,dx(dir+1),1)
      enddo  ! dir

      return
      end subroutine element_size

 
      subroutine SEM_GRAD(ncell,lo,hi,problo,probhi,dx,bfact, &
        P,gx,gy,ubcLEFT,ubcRIGHT, &
        betax,betay,betacell,hflag,probtype,piecewise,weak_flag)
      use LegendreNodes
      use Nodal2DClass

      IMPLICIT NONE

      integer piecewise,weak_flag
      integer lo_e(SDIM),hi_e(SDIM)
      integer ncell(SDIM),lo(SDIM),hi(SDIM)
      REAL_T problo(SDIM),probhi(SDIM),dx(SDIM)
      integer bfact,hflag,probtype,dir
      REAL_T gx(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T gy(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))
      REAL_T betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)
      REAL_T betacell(lo(1):hi(1),lo(SDIM):hi(SDIM))
      REAL_T P(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)
      REAL_T ubcLEFT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1,SDIM)
      REAL_T ubcRIGHT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1,SDIM)
      REAL_T, dimension(:,:), allocatable :: gxR
      REAL_T, dimension(:,:), allocatable :: gyR
      REAL_T, dimension(:,:), allocatable :: IuDwCELL
      REAL_T y(0:bfact)
      REAL_T yGL(0:bfact)
      REAL_T y_extend(0:bfact+1)
      REAL_T yGL_extend(0:bfact+1)
      REAL_T wMAT_extend(0:bfact+1,0:bfact+1)
      REAL_T wMATGL(0:bfact,0:bfact)
      REAL_T wMAT(0:bfact-1,0:bfact-1)
      REAL_T PLINE(0:bfact+1)
      REAL_T PLINE2(0:bfact+1)
      REAL_T PLINE3(0:bfact-1)
      REAL_T PLINE3GL(0:bfact)
      REAL_T PLINE_R(0:bfact)
      REAL_T PLINE_L(0:bfact)
      REAL_T yleft(0:bfact)
      REAL_T yright(0:bfact)
      REAL_T pleft(0:bfact)
      REAL_T pright(0:bfact)
      REAL_T w_right(0:bfact,0:bfact)
      REAL_T w_left(0:bfact,0:bfact)
      REAL_T sum,denom
      REAL_T DP(0:bfact)
      REAL_T dx_element,xL,xR,dxSEM
      integer i,j,i1,j1,i2,j2
      REAL_T ubc(0:1)
      REAL_T IuDw(0:bfact)
      REAL_T xlo_e(SDIM),xhi_e(SDIM)
      integer dirbc,bctype,side
      REAL_T bc_sign
      REAL_T ULEFT,URIGHT,FLUXLEFT,FLUXRIGHT

      TYPE (DGprecomputed) spA
      CHARACTER(LEN=2) :: typenodes  

      typenodes="Ga"

      if ((weak_flag.lt.0).or.(weak_flag.gt.4)) then
       print *,"weak_flag invalid in SEM_GRAD"
       stop
      endif 
      if ((hflag.ne.0).and.(hflag.ne.1)) then
       print *,"hflag invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
 
      allocate(gxR(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))) 
      allocate(gyR(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)) 

      do dir=1,SDIM
       lo_e(dir)=lo(dir)/bfact
       hi_e(dir)=((hi(dir)+1)/bfact)-1
      enddo
 
      call set_boundary(lo,hi,dx,hflag,probtype,P,bfact, &
        problo,probhi,piecewise,weak_flag)

      do i1=0,bfact-1
       y(i1)=cache_gauss(bfact,i1)
      enddo
      do i1=0,bfact
       yGL(i1)=cache_gauss_lobatto(bfact,i1)
      enddo
      y_extend(0)=-one-abs(y(0)+one)
      do i=1,bfact
       y_extend(i)=y(i-1)
      enddo
      y_extend(bfact+1)=one+abs(one-y(bfact-1))
      do i1=0,bfact+1
       yGL_extend(i1)=cache_gauss_lobatto(bfact+1,i1)
      enddo

      yright(0)=-one
      yleft(bfact)=one
      do i1=1,bfact
       yright(i1)=y(i1-1)
       yleft(i1-1)=y(i1-1)
      enddo


       ! polyinterp_Dmatrix(r,x,w)
       ! x(0:r), w(0:r,0:r), w_ij=L_i'(x_j)
      call polyinterp_Dmatrix(bfact+1,yGL_extend,wMAT_extend)
      call polyinterp_Dmatrix(bfact,yGL,wMATGL)
      call polyinterp_Dmatrix(bfact-1,y,wMAT)
      call polyinterp_Dmatrix(bfact,yright,w_right)
      call polyinterp_Dmatrix(bfact,yleft,w_left)

      if (piecewise.eq.0) then

       if ((weak_flag.eq.0).or.(weak_flag.eq.4)) then

        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

          ! DP/DX LOOP
         do j1=0,bfact-1

          if (weak_flag.eq.0) then

           do i1=-1,bfact
            dir=-2 
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            PLINE(i1+1)=P(i2,j2)
           enddo

           ! PLINE: extended Legendre gauss 
           ! DP   : Legendre gauss lobatto
           call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
            y_extend,yGL_extend)
           dx_element=dx(1)*bfact
           call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
            wMAT_extend,yGL_extend,yGL,dx_element)  

           do i1=0,bfact-1
            dir=0
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            gxR(i2,j2)=DP(i1)
           enddo

           do i1=1,bfact
            dir=0
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            gx(i2,j2)=DP(i1)
           enddo

            ! DP/DX LOOP
          else if (weak_flag.eq.4) then

           do i1=0,bfact-1
            dir=-1 
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            PLINE3(i1)=P(i2,j2)
           enddo

           ! PLINE3: Legendre gauss 
           ! PLINE3GL: Legendre gauss Lobatto
           ! DP   : Legendre gauss lobatto
           call poly_change_basis(bfact-1,bfact,PLINE3,PLINE3GL,y,yGL)
           dx_element=dx(1)*bfact
           call deriv_change_basis(bfact,bfact,PLINE3GL,DP, &
            wMATGL,yGL,yGL,dx_element)  

           do i1=0,bfact
            dir=0
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            gxR(i2,j2)=DP(i1)
            gx(i2,j2)=DP(i1)
           enddo

            ! LEFT FACE

           do i1=0,bfact-1
             PLINE_R(i1+1)=PLINE3(i1)
           enddo

           if (i.eq.lo_e(1)) then
             i1=-1
             dir=-2
             call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
             PLINE_R(0)=P(i2,j2)
           else 
             do i1=1,bfact
              pright(i1)=PLINE_R(i1) 
              dir=-1 
              call abs_index(i-1,j,i1-1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
              pleft(i1-1)=P(i2,j2)
             enddo
               ! w_ij=L_i'(x_j)
             sum=zero
             do i1=1,bfact
              sum=sum+pright(i1)*w_right(i1,0)- &
                      pleft(i1-1)*w_left(i1-1,bfact)
             enddo
             denom=w_right(0,0)-w_left(bfact,bfact)
             if (denom.eq.zero) then
              print *,"denom invalid"
              stop
             endif
             PLINE_R(0)=-sum/denom 
           endif

           sum=zero
           do i1=0,bfact
             sum=sum+PLINE_R(i1)*w_right(i1,0)
           enddo
           sum=sum*two/dx_element

           i1=0
           dir=0
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           gxR(i2,j2)=sum
           gx(i2,j2)=sum

             ! RIGHT FACE

           do i1=0,bfact-1
             PLINE_L(i1)=PLINE3(i1)
           enddo

           if (i.eq.hi_e(1)) then
             i1=bfact
             dir=-2
             call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
             PLINE_L(bfact)=P(i2,j2)
           else 
             do i1=1,bfact
              pleft(i1-1)=PLINE_L(i1-1) 
              dir=-1 
              call abs_index(i+1,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
              pright(i1)=P(i2,j2)
             enddo
             sum=zero
             do i1=1,bfact
              sum=sum+pright(i1)*w_right(i1,0)- &
                      pleft(i1-1)*w_left(i1-1,bfact)
             enddo
             denom=w_right(0,0)-w_left(bfact,bfact)
             if (denom.eq.zero) then
              print *,"denom invalid"
              stop
             endif
             PLINE_L(bfact)=-sum/denom 
           endif

           sum=zero
           do i1=0,bfact
             sum=sum+PLINE_L(i1)*w_left(i1,bfact)
           enddo
           sum=sum*two/dx_element

           i1=bfact
           dir=0
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           gxR(i2,j2)=sum
           gx(i2,j2)=sum

          else
           print *,"weak_flag invalid"
           stop
          endif

         enddo ! j1

          ! DP/DY LOOP
         do i1=0,bfact-1

          if (weak_flag.eq.0) then

           do j1=-1,bfact
            dir=-3
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            PLINE(j1+1)=P(i2,j2)
           enddo

           ! PLINE: extended Legendre gauss
           ! DU   : Legendre gauss lobatto
           call poly_change_basis(bfact+1,bfact+1,PLINE,PLINE2, &
            y_extend,yGL_extend)
           dx_element=dx(SDIM)*bfact
           call deriv_change_basis(bfact+1,bfact,PLINE2,DP, &
            wMAT_extend,yGL_extend,yGL,dx_element)

           do j1=0,bfact-1
            dir=1
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            gyr(i2,j2)=DP(j1)
           enddo

           do j1=1,bfact
            dir=1
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            gy(i2,j2)=DP(j1)
           enddo

           ! DP/DY LOOP
          else if (weak_flag.eq.4) then

           do j1=0,bfact-1
            dir=-1 
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            PLINE3(j1)=P(i2,j2)
           enddo

           ! PLINE3: Legendre gauss 
           ! PLINE3GL: Legendre gauss Lobatto
           ! DP   : Legendre gauss lobatto
           call poly_change_basis(bfact-1,bfact,PLINE3,PLINE3GL,y,yGL)
           dx_element=dx(SDIM)*bfact
           call deriv_change_basis(bfact,bfact,PLINE3GL,DP, &
            wMATGL,yGL,yGL,dx_element)  

           do j1=0,bfact
            dir=SDIM-1
            call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
            gyR(i2,j2)=DP(j1)
            gy(i2,j2)=DP(j1)
           enddo

            ! LEFT FACE

           do j1=0,bfact-1
             PLINE_R(j1+1)=PLINE3(j1)
           enddo

           if (j.eq.lo_e(SDIM)) then
             j1=-1
             dir=-3
             call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
             PLINE_R(0)=P(i2,j2)
           else 
             do j1=1,bfact
              pright(j1)=PLINE_R(j1) 
              dir=-1 
              call abs_index(i,j-1,i1,j1-1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
              pleft(j1-1)=P(i2,j2)
             enddo
             sum=zero
             do j1=1,bfact
              sum=sum+pright(j1)*w_right(j1,0)- &
                      pleft(j1-1)*w_left(j1-1,bfact)
             enddo
             denom=w_right(0,0)-w_left(bfact,bfact)
             if (denom.eq.zero) then
              print *,"denom invalid"
              stop
             endif
             PLINE_R(0)=-sum/denom 
           endif

           sum=zero
           do j1=0,bfact
             sum=sum+PLINE_R(j1)*w_right(j1,0)
           enddo
           sum=sum*two/dx_element

           j1=0
           dir=SDIM-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           gyR(i2,j2)=sum
           gy(i2,j2)=sum

             ! RIGHT FACE

           do j1=0,bfact-1
             PLINE_L(j1)=PLINE3(j1)
           enddo

           if (j.eq.hi_e(SDIM)) then
             j1=bfact
             dir=-3
             call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
             PLINE_L(bfact)=P(i2,j2)
           else 
             do j1=1,bfact
              pleft(j1-1)=PLINE_L(j1-1) 
              dir=-1 
              call abs_index(i,j+1,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
              pright(j1)=P(i2,j2)
             enddo
             sum=zero
             do j1=1,bfact
              sum=sum+pright(j1)*w_right(j1,0)- &
                      pleft(j1-1)*w_left(j1-1,bfact)
             enddo
             denom=w_right(0,0)-w_left(bfact,bfact)
             if (denom.eq.zero) then
              print *,"denom invalid"
              stop
             endif
             PLINE_L(bfact)=-sum/denom 
           endif

           sum=zero
           do j1=0,bfact
             sum=sum+PLINE_L(j1)*w_left(j1,bfact)
           enddo
           sum=sum*two/dx_element

           j1=bfact
           dir=SDIM
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           gyR(i2,j2)=sum
           gy(i2,j2)=sum

          else
           print *,"weak_flag invalid"
           stop
          endif

         enddo ! i1
        
        enddo
        enddo

        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)
         do j1=0,bfact-1
          do i1=0,bfact
           dir=0
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           if ((i1.eq.0).and.(i.eq.lo_e(1))) then
            gx(i2,j2)=gxR(i2,j2)
           else if ((i1.eq.bfact).and.(i.eq.hi_e(1))) then
            ! do nothing
           else if ((i1.gt.0).and.(i1.lt.bfact)) then
            ! do nothing
           else
            gx(i2,j2)=0.5*(gx(i2,j2)+gxR(i2,j2))
           endif
          enddo
         enddo
         do i1=0,bfact-1
          do j1=0,bfact
           dir=1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,dir)
           if ((j1.eq.0).and.(j.eq.lo_e(SDIM))) then
            gy(i2,j2)=gyR(i2,j2)
           else if ((j1.eq.bfact).and.(j.eq.hi_e(SDIM))) then
            ! do nothing
           else if ((j1.gt.0).and.(j1.lt.bfact)) then
            ! do nothing
           else
            gy(i2,j2)=0.5*(gy(i2,j2)+gyR(i2,j2))
           endif
          enddo
         enddo
        enddo
        enddo

       else if ((weak_flag.ge.1).and.(weak_flag.le.3)) then

        call ConstructNodal2D(bfact-1,spA,typenodes)
 
        allocate(IuDwCELL(lo(1):hi(1),lo(SDIM):hi(SDIM))) 

         ! gx 
        dir=1
        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

         call element_size(problo,i,j,lo_e,hi_e,bfact, &
           lo,hi,dx,xlo_e,xhi_e)
         do j1=0,bfact-1
          do i1=0,bfact-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           PLINE(i1)=P(i2,j2)
          enddo
          call GradientWeakFormCompute(spA,PLINE,IuDw,ubc)
          do i1=0,bfact-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           IuDwCELL(i2,j2)=IuDw(i1)
          enddo

          i1=0
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
          ubcLEFT(i,j2,dir)=ubc(0)
          ubcRIGHT(i,j2,dir)=ubc(1)

          if (i.eq.lo_e(dir)) then
           dirbc=dir
           side=1
           call get_bctype(dirbc,side,bctype,bc_sign,probtype)
           if (bctype.eq.0) then ! dirichlet
            ubcRIGHT(i-1,j2,dir)=P(lo(dir)-1,j2)
           else if (bctype.eq.1) then ! neumann
            ubcRIGHT(i-1,j2,dir)=ubc(0)
           else
            print *,"bctype invalid"
            stop
           endif
          endif
          if (i.eq.hi_e(dir)) then
           dirbc=dir
           side=2
           call get_bctype(dirbc,side,bctype,bc_sign,probtype)
           if (bctype.eq.0) then ! dirichlet
            ubcLEFT(i+1,j2,dir)=P(hi(dir)+1,j2)
           else if (bctype.eq.1) then  ! neumann
            ubcLEFT(i+1,j2,dir)=ubc(1)
           else
            print *,"bctype invalid"
            stop
           endif
          endif

         enddo !j1

        enddo
        enddo

        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

         call element_size(problo,i,j,lo_e,hi_e,bfact, &
           lo,hi,dx,xlo_e,xhi_e)
         do j1=0,bfact-1
          i1=0
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)

          if (weak_flag.eq.1) then
           ULEFT=ubcRIGHT(i-1,j2,dir)
           URIGHT=ubcRIGHT(i,j2,dir)
           if (i.eq.hi_e(dir)) then
            URIGHT=ubcLEFT(i+1,j2,dir)
           endif
          else if ((weak_flag.eq.2).or.(weak_flag.eq.3)) then
           if (i.eq.lo_e(dir)) then
            ULEFT=ubcRIGHT(i-1,j2,dir)
           else
            ULEFT=half*(ubcRIGHT(i-1,j2,dir)+ubcLEFT(i,j2,dir))
           endif
           if (i.eq.hi_e(dir)) then
            URIGHT=ubcLEFT(i+1,j2,dir)
           else
            URIGHT=half*(ubcRIGHT(i,j2,dir)+ubcLEFT(i+1,j2,dir))
           endif
          else
           print *,"weak_flag invalid"
           stop
          endif 


          do i1=0,bfact-1
           FLUXLEFT=ULEFT*spA%Lminus1(i1)
           FLUXRIGHT=URIGHT*spA%L1(i1)
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           gx(i2,j2)=(IuDwCELL(i2,j2)+(FLUXRIGHT-FLUXLEFT))/ &
             (half*(xhi_e(dir)-xlo_e(dir))*spA%Qweights(i1))
          enddo

         enddo !j1
        enddo
        enddo
          
         ! gy 
        dir=2
        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

         call element_size(problo,i,j,lo_e,hi_e,bfact, &
           lo,hi,dx,xlo_e,xhi_e)
         do i1=0,bfact-1
          do j1=0,bfact-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           PLINE(j1)=P(i2,j2)
          enddo
          call GradientWeakFormCompute(spA,PLINE,IuDw,ubc)
          do j1=0,bfact-1
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           IuDwCELL(i2,j2)=IuDw(j1)
          enddo

          j1=0
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
          ubcLEFT(i2,j,dir)=ubc(0)
          ubcRIGHT(i2,j,dir)=ubc(1)

          if (j.eq.lo_e(dir)) then
           dirbc=dir
           side=1
           call get_bctype(dirbc,side,bctype,bc_sign,probtype)
           if (bctype.eq.0) then ! dirichlet
            ubcRIGHT(i2,j-1,dir)=P(i2,lo(dir)-1)
           else if (bctype.eq.1) then  ! neumann
            ubcRIGHT(i2,j-1,dir)=ubc(0)
           else
            print *,"bctype invalid"
            stop
           endif
          endif
          if (j.eq.hi_e(dir)) then
           dirbc=dir
           side=2
           call get_bctype(dirbc,side,bctype,bc_sign,probtype)
           if (bctype.eq.0) then  ! dirichlet
            ubcLEFT(i2,j+1,dir)=P(i2,hi(dir)+1)
           else if (bctype.eq.1) then ! neumann
            ubcLEFT(i2,j+1,dir)=ubc(1)
           else
            print *,"bctype invalid"
            stop
           endif
          endif

         enddo !i1

        enddo
        enddo

        do i=lo_e(1),hi_e(1)
        do j=lo_e(SDIM),hi_e(SDIM)

         call element_size(problo,i,j,lo_e,hi_e,bfact, &
            lo,hi,dx,xlo_e,xhi_e)
         do i1=0,bfact-1
          j1=0
          call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)

          if (weak_flag.eq.1) then
           ULEFT=ubcRIGHT(i2,j-1,dir)
           URIGHT=ubcRIGHT(i2,j,dir)
           if (j.eq.hi_e(dir)) then
            URIGHT=ubcLEFT(i2,j+1,dir)
           endif
          else if ((weak_flag.eq.2).or.(weak_flag.eq.3)) then
           if (j.eq.lo_e(dir)) then
            ULEFT=ubcRIGHT(i2,j-1,dir)
           else
            ULEFT=half*(ubcRIGHT(i2,j-1,dir)+ubcLEFT(i2,j,dir))
           endif
           if (j.eq.hi_e(dir)) then
            URIGHT=ubcLEFT(i2,j+1,dir)
           else
            URIGHT=half*(ubcRIGHT(i2,j,dir)+ubcLEFT(i2,j+1,dir))
           endif
          else
           print *,"weak_flag invalid"
           stop
          endif 

          do j1=0,bfact-1
           FLUXLEFT=ULEFT*spA%Lminus1(j1)
           FLUXRIGHT=URIGHT*spA%L1(j1)
           call abs_index(i,j,i1,j1,i2,j2,lo_e,hi_e,bfact,lo,hi,-1)
           gy(i2,j2)=(IuDwCELL(i2,j2)+(FLUXRIGHT-FLUXLEFT))/ &
             (half*(xhi_e(dir)-xlo_e(dir))*spA%Qweights(j1))
          enddo
         enddo !i1
        enddo
        enddo

        deallocate(IuDwCELL)
        call DestructNodal2D(spA)

       else 
        print *,"weak_flag invalid SEM_grad"
        stop
       endif

      else if (piecewise.eq.1) then

       dir=1
       do i=lo(1),hi(1)+1
       do j=lo(SDIM),hi(SDIM)
        call gridloc1D(xL,problo(dir),i-1,lo(dir),bfact,dx(dir),0)
        call gridloc1D(xR,problo(dir),i,lo(dir),bfact,dx(dir),0)
        dxSEM=xR-xL
        if (bfact.eq.1) then
         if (abs(dxSEM-dx(dir)).gt.1.0E-13) then
          print *,"dxSEM invalid"
          print *,"dir,dxSEM,dx(dir) ",dir,dxSEM,dx(dir)
          print *,"dir,lo,hi ",dir,lo(dir),hi(dir)
          print *,"i,j,xR,xL ",i,j,xR,xL
          stop
         endif
        endif
        gx(i,j)=(P(i,j)-P(i-1,j))/dxSEM
       enddo
       enddo
       dir=2
       do i=lo(1),hi(1)
       do j=lo(SDIM),hi(SDIM)+1
        call gridloc1D(xL,problo(dir),j-1,lo(dir),bfact,dx(dir),0)
        call gridloc1D(xR,problo(dir),j,lo(dir),bfact,dx(dir),0)
        dxSEM=xR-xL
        if (bfact.eq.1) then
         if (abs(dxSEM-dx(dir)).gt.1.0E-13) then
          print *,"dxSEM invalid"
          print *,"dir,dxSEM,dx(dir) ",dir,dxSEM,dx(dir)
          print *,"dir,lo,hi ",dir,lo(dir),hi(dir)
          print *,"i,j,xR,xL ",i,j,xR,xL
          stop
         endif
        endif
        gy(i,j)=(P(i,j)-P(i,j-1))/dxSEM
       enddo
       enddo
      else
       print *,"piecewise invalid"
       stop
      endif

      if ((weak_flag.eq.0).or.(weak_flag.eq.4).or.(piecewise.eq.1)) then

       do i=lo(1),hi(1)+1
       do j=lo(SDIM),hi(SDIM)
        gx(i,j)=gx(i,j)*betax(i,j)
       enddo
       enddo
       do i=lo(1),hi(1)
       do j=lo(SDIM),hi(SDIM)+1
        gy(i,j)=gy(i,j)*betay(i,j)
       enddo
       enddo

      else if (((weak_flag.ge.1).and.(weak_flag.le.3)).and. &
               (piecewise.eq.0)) then

       do i=lo(1),hi(1)
       do j=lo(SDIM),hi(SDIM)
        gx(i,j)=gx(i,j)*betacell(i,j)
       enddo
       enddo
       do i=lo(1),hi(1)
       do j=lo(SDIM),hi(SDIM)
        gy(i,j)=gy(i,j)*betacell(i,j)
       enddo
       enddo

      else
       print *,"invalid weak_flag or piecewise"
       stop
      endif

      deallocate(gxR) 
      deallocate(gyR)
 
      return
      end subroutine SEM_GRAD



      REAL_T function hs(phi,cutoff)
      use LegendreNodes
      IMPLICIT NONE
      REAL_T phi,cutoff,EPS,local_pi

      local_pi=four*atan(one)
      EPS=1.0E-6
      if (phi.ge.cutoff) then
        hs=one
      elseif (phi.le.-cutoff) then
        hs=zero
      else
        hs=phi/(two*cutoff)+ &
         sin(local_pi*phi/cutoff)/(two*local_pi)+half
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
   
          
      subroutine get_U(x,dx,U,probtype,bfact,weak_flag)
      use LegendreNodes
      IMPLICIT NONE

      REAL_T x(SDIM),dx(SDIM),U,phi,H,cutoff,hs,localPI,cutoffbeta
      integer probtype,bfact,weak_flag

      localPI=four*atan(one)
      cutoff=dx(1)*THICKVEL
      cutoffbeta=dx(1)*THICKBETA
      if (probtype.eq.0) then
       U=1.0
       phi=sqrt( (x(1)-0.5)**2+(x(2)-0.5)**2 )-0.25
       H=hs(phi-cutoffbeta+cutoff,cutoff)
       U=H
      else if (probtype.eq.1) then
       if (weak_flag.eq.0) then
        U=bfact*(x(1)**(bfact-one))
       else if ((weak_flag.ge.1).and.(weak_flag.le.4)) then
        if (bfact.lt.2) then
         U=zero
        else
         U=(bfact-one)*(x(1)**(bfact-two))
        endif
       else
        print *,"weak_flag invalid"
        stop
       endif
      else if (probtype.eq.2) then
       U=1.0
       phi=sqrt( (x(1)-0.5)**2+(x(2)-0.25)**2 )-0.15
       H=hs(phi,cutoff)
       U=H
      else if (probtype.eq.3) then
       U=localPI*cos(localPI*x(1))*sin(localPI*x(2))
      else
       print *,"probtype invalid" 
       stop
      endif

      return
      end subroutine get_U


      subroutine get_PEXACT(x,dx,P,probtype,bfact,weak_flag)
      use LegendreNodes
      IMPLICIT NONE

      REAL_T x(SDIM),dx(SDIM),P,localPI
      integer probtype,bfact,weak_flag

      localPI=four*atan(one)

      if (probtype.eq.0) then
       P=0.0
      else if (probtype.eq.1) then
       if (weak_flag.eq.0) then
        P=x(1)**bfact+x(2)**bfact
       else if ((weak_flag.ge.1).and.(weak_flag.le.4)) then
        P=x(1)**(bfact-one)+x(2)**(bfact-one)
       else
        print *,"weak_flag invalid"
        stop
       endif
      else if (probtype.eq.2) then
       P=0.0
      else if (probtype.eq.3) then
       P=sin(localPI*x(1))*sin(localPI*x(2))
       P=0.0
      else
       print *,"probtype invalid"
       stop
      endif
      
      return
      end subroutine get_PEXACT

      subroutine get_V(x,dx,V,probtype,bfact,weak_flag)
      use LegendreNodes
      IMPLICIT NONE

      REAL_T x(SDIM),dx(SDIM),V,localPI
      integer probtype,bfact,weak_flag

      localPI=four*atan(one)
      if (probtype.eq.0) then
       V=0.0
      else if (probtype.eq.1) then
       if (weak_flag.eq.0) then
        V=bfact*(x(2)**(bfact-one))
       else if ((weak_flag.ge.1).and.(weak_flag.le.4)) then
        if (bfact.lt.2) then
         V=zero
        else
         V=(bfact-one)*(x(2)**(bfact-two))
        endif
       else
        print *,"weak_flag invalid"
        stop
       endif
      else if (probtype.eq.2) then
       V=0.0
      else if (probtype.eq.3) then
       V=localPI*cos(localPI*x(2))*sin(localPI*x(1))
      else
       print *,"probtype invalid"
       stop
      endif
      
      return
      end subroutine get_V



      subroutine get_beta(x,dx,beta,probtype)
      use LegendreNodes
      IMPLICIT NONE

      REAL_T x(SDIM),dx(SDIM),beta,cutoff,phi,H,hs
      integer probtype

      cutoff=dx(1)*THICKBETA

      if (probtype.eq.0) then
       beta=one
       phi=sqrt( (x(1)-0.5)**2+(x(2)-0.5)**2 )-0.25
       H=hs(phi,cutoff)
       beta=one/(H+DEN_OBSTACLE*(one-H))
      else if (probtype.eq.1) then
       beta=one
      else if (probtype.eq.2) then
       beta=1.0
       phi=sqrt( (x(1)-0.5)**2+(x(2)-0.25)**2 )-0.15
       H=hs(phi,cutoff)
       beta=one/(H+DEN_OBSTACLE*(one-H))
      else if (probtype.eq.3) then
       beta=one
      else
       print *,"probtype invalid"
       stop
      endif
      
      return
      end subroutine get_beta

       ! ii=0 Gauss, ii=1 Gauss Lobatto
      subroutine gridloc1D(x,problo,i,lo,bfact,dx,ii)
      use LegendreNodes

      IMPLICIT NONE

      REAL_T x,problo,dx,xnode
      integer i,lo,bfact,ii,i1,i_e,ireflect

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

      if (i.ge.lo) then
       i_e=(i-lo)/bfact
       i1=i-lo-i_e*bfact

       if (ii.eq.0) then
        xnode=cache_gauss(bfact,i1)
       else if (ii.eq.1) then
        xnode=cache_gauss_lobatto(bfact,i1)
       else
        print *,"ii invalid"
        stop
       endif

       x=problo+i_e*bfact*dx+(xnode+1.0)*bfact*dx/2.0
      else if (i.lt.lo) then
       ireflect=lo+(lo-i-1+ii)
       i_e=(ireflect-lo)/bfact
       i1=ireflect-lo-i_e*bfact

       if (ii.eq.0) then
        xnode=cache_gauss(bfact,i1)
       else if (ii.eq.1) then
        xnode=cache_gauss_lobatto(bfact,i1)
       else
        print *,"ii invalid"
        stop
       endif

       x=problo-i_e*bfact*dx-(xnode+1.0)*bfact*dx/2.0
      endif

      return
      end subroutine gridloc1D

      subroutine gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,dir)
      IMPLICIT NONE

      REAL_T xL(SDIM),xR(SDIM),x(SDIM),problo(SDIM)
      integer i,j,lo(SDIM),bfact,dir,dir2,index
      integer ii(SDIM)
      REAL_T dx(SDIM)

      do dir2=1,SDIM
       if (bfact*(lo(dir2)/bfact).ne.lo(dir2)) then
        print *,"lo invalid"
        stop
       endif
      enddo

      if ((dir.lt.-1).or.(dir.ge.SDIM)) then
       print *,"dir invalid"
       stop
      endif
      ii(1)=0
      ii(2)=0
      if (dir.ge.0) then
       ii(dir+1)=1
      endif

      do dir2=1,SDIM
       if (dir2.eq.1) then
        index=i
       else
        index=j
       endif

       call gridloc1D(x(dir2),problo(dir2),index,lo(dir2),bfact, &
        dx(dir2),ii(dir2))
       if (ii(dir2).eq.1) then
        call gridloc1D(xL(dir2),problo(dir2),index-1,lo(dir2),bfact, &
         dx(dir2),1-ii(dir2))
        call gridloc1D(xR(dir2),problo(dir2),index,lo(dir2),bfact, &
         dx(dir2),1-ii(dir2))
       else
        call gridloc1D(xL(dir2),problo(dir2),index,lo(dir2),bfact, &
         dx(dir2),1-ii(dir2))
        call gridloc1D(xR(dir2),problo(dir2),index+1,lo(dir2),bfact, &
         dx(dir2),1-ii(dir2))
       endif
      enddo

      return
      end subroutine gridloc

      subroutine FVgeom(problo,i,j,lo,bfact,dx, &
        AX,AY,dxL,dxR,vol)
      use LegendreNodes
      IMPLICIT NONE

      REAL_T problo(SDIM)
      integer i,j,lo(SDIM)
      integer bfact,dir
      REAL_T dx(SDIM)
      REAL_T dxL(SDIM)
      REAL_T dxR(SDIM)
      REAL_T AX,AY,vol,xnbr
      REAL_T x(SDIM)
      REAL_T xL(SDIM)
      REAL_T xR(SDIM)
      
      dir=-1
      call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,dir)
      AX=xR(SDIM)-xL(SDIM)
      AY=xR(1)-xL(1)
      vol=AX*AY
      dir=1
      call gridloc1D(xnbr,problo(dir),i-1,lo(dir),bfact,dx(dir),0)
      dxL(dir)=x(dir)-xnbr
      call gridloc1D(xnbr,problo(dir),i+1,lo(dir),bfact,dx(dir),0)
      dxR(dir)=xnbr-x(dir)
      dir=2
      call gridloc1D(xnbr,problo(dir),j-1,lo(dir),bfact,dx(dir),0)
      dxL(dir)=x(dir)-xnbr
      call gridloc1D(xnbr,problo(dir),j+1,lo(dir),bfact,dx(dir),0)
      dxR(dir)=xnbr-x(dir)

      return
      end subroutine FVgeom 

      subroutine den_cell_to_mac( &
           
      subroutine ABEC_solver( &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       presbc, &
       xlo,dx, &
       vol,DIMS(vol), &
  
      program main
      use LegendreNodes

      IMPLICIT NONE
    
      integer i,j,hflag,precond_type,dir,for_source
      integer lo(SDIM),hi(SDIM),bfact,ncell(SDIM)
      REAL_T  dx(SDIM)
      REAL_T  x(SDIM)
      REAL_T  xL(SDIM)
      REAL_T  xR(SDIM)
      REAL_T  problo(SDIM)
      REAL_T  probhi(SDIM)
      integer probtype
      REAL_T tol
      REAL_T, dimension(:,:), allocatable :: P
      REAL_T, dimension(:,:), allocatable :: U
      REAL_T, dimension(:,:), allocatable :: V
      REAL_T, dimension(:,:), allocatable :: G
      REAL_T, dimension(:,:), allocatable :: gx
      REAL_T, dimension(:,:), allocatable :: gy
      REAL_T, dimension(:,:), allocatable :: betax
      REAL_T, dimension(:,:), allocatable :: betay
      REAL_T, dimension(:,:), allocatable :: betacell
      REAL_T, dimension(:,:), allocatable :: UCELL
      REAL_T, dimension(:,:), allocatable :: VCELL
      REAL_T, dimension(:,:), allocatable :: BXCELL
      REAL_T, dimension(:,:), allocatable :: BYCELL
      REAL_T, dimension(:,:), allocatable :: GXCELL
      REAL_T, dimension(:,:), allocatable :: GYCELL
      REAL_T, dimension(:,:,:), allocatable :: ubcLEFT
      REAL_T, dimension(:,:,:), allocatable :: ubcRIGHT
      integer bottomflag,piecewise_finest,weak_flag,typ_nodes,nplot
      REAL_T hplot,xplot,yplot,cutoff,hs
      REAL_T y(0:64)
      REAL_T w(0:64)
      REAL_T data(0:64)
     
      character*9 gibbsdatafile
      character*8 wavedatafile
      character*9 debugdatafile
      character*8 celldatafile
      character*9 cell2datafile


       ! things to try:
       ! 1. Helmholtz equation
       ! 2. weighted Jacobi smoother
       ! 3. GSRB smoother
       ! 4. ILURB smoother
       ! 5. stokes solver

       ! probtype=0 flow past cylinder
       ! probtype=1 sanity check
       ! probtype=3 p=sin(pi x) sin(pi y)
      probtype=0
      problo(1)=0.0
      problo(2)=0.0
      probhi(1)=1.0
      probhi(2)=1.0
      ncell(1)=128
      ncell(2)=128
      bfact=16
      call sanity_check(bfact)
      typ_nodes=0  ! 0=Legendre  1=Clenshaw Curtis
      piecewise_finest=0
       ! weak_flag=0 structured grid pseudo spectral element method.
       ! weak_flag=1,2,3 Discontinuous Galerkin spectral element method.
       ! weak_flag=4 structured or unstructured grid Pseudo spectral element
       !             method.
      weak_flag=4

      do dir=1,SDIM
       dx(dir)=(probhi(dir)-problo(dir))/ncell(dir)
       if (bfact*(ncell(dir)/bfact).ne.ncell(dir)) then
        print *,"ncell must be a multiple of bfact"
        stop
       endif
      enddo
      tol=1.0E-7

      do dir=1,SDIM
       lo(dir)=0
       hi(dir)=ncell(dir)-1
      enddo

      precond_type=2 ! 0 M=I  1=Jacobi precond.  2=multigrid precond.
      hflag=0

      print *,"Running amazing code with the following parameters:"
      print *,"weak_flag= ",weak_flag
      print *,"bfact= ",bfact
      print *,"typ_nodes= ",typ_nodes
      print *,"ncell= ",ncell(1),ncell(2)
      print *,"piecewise=",piecewise_finest
      print *,"precond_type= ",precond_type
      print *,"probtype= ",probtype
      print *,"density ratio (if probtype=0) ",DEN_OBSTACLE
      print *,"beta thickness parameter (if probtype=0) ",THICKBETA
      print *,"vel thickness parameter (if probtype=0) ",THICKVEL
 
      allocate(P(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(G(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(U(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))) 
      allocate(V(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)) 
      allocate(gx(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))) 
      allocate(gy(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)) 
      allocate(betax(lo(1):hi(1)+1,lo(SDIM):hi(SDIM))) 
      allocate(betay(lo(1):hi(1),lo(SDIM):hi(SDIM)+1)) 
      allocate(betacell(lo(1):hi(1),lo(SDIM):hi(SDIM))) 
      allocate(UCELL(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(VCELL(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(BXCELL(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(BYCELL(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(GXCELL(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(GYCELL(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1)) 
      allocate(ubcLEFT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1,SDIM)) 
      allocate(ubcRIGHT(lo(1)-1:hi(1)+1,lo(SDIM)-1:hi(SDIM)+1,SDIM)) 

      call init_cache(bfact+1,typ_nodes)

      write(gibbsdatafile,'(A9)') 'gibbsdata'
      print *,"gibbsdatafile ",gibbsdatafile
      open(unit=25,file=gibbsdatafile)

      cutoff=two/bfact
      cutoff=zero*cutoff
      do i=0,bfact-1
       y(i)=cache_gauss(bfact,i)
       data(i)=hs(y(i),cutoff)
      enddo
      nplot=1000
      hplot=two/nplot
      do i=0,nplot
       xplot=-one+i*hplot 
       call polyinterp_weights(bfact-1,y,w,xplot)
       call do_polyinterp(bfact-1,w,data,yplot)
       write(25,*) xplot,yplot
      enddo 
      close(25)

      dir=-1
      do i=lo(1),hi(1)
      do j=lo(2),hi(2)
       call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,dir)
       call get_beta(x,dx,betacell(i,j),probtype)
      enddo
      enddo

      dir=0
      do i=lo(1),hi(1)+1
      do j=lo(2),hi(2)
       call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,dir) 
       call get_beta(x,dx,betax(i,j),probtype)
      enddo
      enddo
      dir=1
      do i=lo(1),hi(1)
      do j=lo(2),hi(2)+1
       call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,dir) 
       call get_beta(x,dx,betay(i,j),probtype)
      enddo
      enddo

      if (((weak_flag.ge.1).and.(weak_flag.le.3)).and. &
          (piecewise_finest.eq.0)) then
       dir=-1
       do i=lo(1),hi(1)
       do j=lo(2),hi(2)
        call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,dir)
        call get_U(x,dx,U(i,j),probtype,bfact,weak_flag)
        call get_V(x,dx,V(i,j),probtype,bfact,weak_flag)
       enddo
       enddo
      else if ((weak_flag.eq.0).or.(weak_flag.eq.4).or. &
               (piecewise_finest.eq.1)) then

       dir=0
       do i=lo(1),hi(1)+1
       do j=lo(2),hi(2)
        call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,dir) 
        call get_U(x,dx,U(i,j),probtype,bfact,weak_flag)
       enddo
       enddo
       dir=1
       do i=lo(1),hi(1)
       do j=lo(2),hi(2)+1
        call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,dir) 
        call get_V(x,dx,V(i,j),probtype,bfact,weak_flag)
       enddo
       enddo

      else
       print *,"weak_flag or piecewise_finest invalid"
       stop
      endif

      do i=lo(1),hi(1)
      do j=lo(2),hi(2)
       call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1) 
       call get_PEXACT(x,dx,P(i,j),probtype,bfact,weak_flag)
      enddo 
      enddo 

      hflag=0
      call SEM_GRAD(ncell,lo,hi,problo,probhi,dx,bfact, &
       P,gx,gy,ubcLEFT,ubcRIGHT,betax,betay, &
       betacell,hflag,probtype,piecewise_finest,weak_flag)
      call TFRCELL(lo,hi,problo,probhi,dx,bfact,gx,GXCELL, &
        gy,GYCELL,piecewise_finest,weak_flag)

      hflag=0
      for_source=1
      call SEM_DIV(ncell,lo,hi,problo,probhi,dx,bfact, &
       P,U,V,betax,betay,ubcLEFT,ubcRIGHT, &
       G,hflag,probtype,for_source,piecewise_finest, &
       weak_flag)
      call NEGVEC(lo,hi,G)

      if ((weak_flag.eq.0).or.(weak_flag.eq.4)) then
       call TFRCELL(lo,hi,problo,probhi,dx,bfact,betax,BXCELL, &
        betay,BYCELL,piecewise_finest,weak_flag)
      else
       do i=lo(1),hi(1)
       do j=lo(2),hi(2)
        BXCELL(i,j)=betacell(i,j)
        BYCELL(i,j)=betacell(i,j)
       enddo
       enddo
      endif

      write(debugdatafile,'(A9)') 'debugdata'
      print *,"debugdatafile ",debugdatafile
      open(unit=12,file=debugdatafile)
      if (1.eq.0) then
       write(12,*) '# X,Y,P,G,BX,BY,gx,gy'
      else 
       write(12,*) 'VARIABLES="X","Y","P","G","BX","BY","gx","gy"'
       write(12,*) 'zone i=',hi(1)-lo(1)+1, &
               ' j=',hi(2)-lo(2)+1,' f=point'
      endif

      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
       call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1) 
       write(12,*) x(1),x(2),P(i,j),G(i,j),BXCELL(i,j),BYCELL(i,j), &
         GXCELL(i,j),GYCELL(i,j)
      enddo
      enddo

      close(12)

      write(celldatafile,'(A8)') 'celldata'
      print *,"celldatafile ",celldatafile
      open(unit=23,file=celldatafile)
      if (1.eq.0) then
       write(23,*) '# X,Y,P'
      else
       write(23,*) 'VARIABLES="X","Y","P"'
       write(23,*) 'zone i=',hi(1)-lo(1)+3, &
               ' j=',hi(2)-lo(2)+3,' f=point'
      endif

      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
       call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1)
       write(23,*) x(1),x(2),P(i,j)
      enddo
      enddo

      close(23)

      if (1.eq.1) then
       call bicgstab(ncell,lo,hi,bfact,P,G,betax,betay, &
        betacell,dx,hflag,precond_type,problo,probhi,probtype,tol, &
        piecewise_finest,weak_flag)
      else
       bottomflag=0
       call pcg(ncell,lo,hi,bfact,P,G,betax,betay, &
        betacell,dx,hflag,precond_type,problo,probhi,probtype,tol, &
        piecewise_finest,weak_flag,bottomflag)
      endif

      hflag=0
      call SEM_GRAD(ncell,lo,hi,problo,probhi,dx,bfact, &
       P,gx,gy,ubcLEFT,ubcRIGHT,betax,betay, &
       betacell,hflag,probtype,piecewise_finest,weak_flag)

      if (((weak_flag.ge.1).and.(weak_flag.le.3)).and. &
          (piecewise_finest.eq.0)) then
       do i=lo(1),hi(1)
       do j=lo(2),hi(2)
        U(i,j)=U(i,j)-gx(i,j)
        V(i,j)=V(i,j)-gy(i,j)
       enddo
       enddo
      else if ((weak_flag.eq.0).or.(weak_flag.eq.4).or. &
               (piecewise_finest.eq.1)) then
       do i=lo(1),hi(1)+1
       do j=lo(2),hi(2)
        U(i,j)=U(i,j)-gx(i,j)
       enddo
       enddo
       do i=lo(1),hi(1)
       do j=lo(2),hi(2)+1
        V(i,j)=V(i,j)-gy(i,j)
       enddo
       enddo
      else
       print *,"weak_flag or piecewise_finest invalid"
       stop
      endif
         ! this is done so that the velocity BC are set and check
         ! that div u=0
      hflag=0
      for_source=1
      call SEM_DIV(ncell,lo,hi,problo,probhi,dx,bfact, &
        P,U,V,betax,betay,ubcLEFT,ubcRIGHT, &
        G,hflag,probtype,for_source,piecewise_finest, &
        weak_flag)
      write(wavedatafile,'(A8)') 'wavedata'
      print *,"wavedatafile ",wavedatafile
      open(unit=11,file=wavedatafile)
      if (1.eq.0) then
       write(11,*) '# X,Y,P,U,V,DIV'
      else 
       write(11,*) 'VARIABLES="X","Y","P","U","V","DIV"'
       write(11,*) 'zone i=',hi(1)-lo(1)+1, &
               ' j=',hi(2)-lo(2)+1,' f=point'
      endif
        ! transfer from MAC grid to CELL grid.
      call TFRCELL(lo,hi,problo,probhi,dx,bfact,U,UCELL, &
        V,VCELL,piecewise_finest,weak_flag)

      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
       call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1) 
       write(11,*) x(1),x(2),P(i,j),UCELL(i,j),VCELL(i,j),G(i,j)
      enddo
      enddo

      close(11)


      write(cell2datafile,'(A9)') 'cell2data'
      print *,"cell2datafile ",cell2datafile
      open(unit=24,file=cell2datafile)
      if (1.eq.0) then
       write(24,*) '# X,Y,P'
      else
       write(24,*) 'VARIABLES="X","Y","P"'
       write(24,*) 'zone i=',hi(1)-lo(1)+3, &
               ' j=',hi(2)-lo(2)+3,' f=point'
      endif

      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
       call gridloc(xL,xR,x,problo,i,j,lo,bfact,dx,-1)
       write(24,*) x(1),x(2),P(i,j)
      enddo
      enddo

      close(24)


 
      deallocate(P) 
      deallocate(G) 
      deallocate(U) 
      deallocate(V) 
      deallocate(UCELL) 
      deallocate(VCELL) 
      deallocate(BXCELL) 
      deallocate(BYCELL) 
      deallocate(GXCELL) 
      deallocate(GYCELL) 
      deallocate(gx) 
      deallocate(gy) 
      deallocate(betax) 
      deallocate(betay) 
      deallocate(betacell) 
      deallocate(ubcLEFT) 
      deallocate(ubcRIGHT) 

      call delete_cache()
 
      return
      end
