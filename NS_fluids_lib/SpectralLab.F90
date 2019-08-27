#define one 1.0d0
#define two 2.0d0
#define three 3.0d0
#define four 4.0d0
#define fourth 0.25d0
#define zero 0.0d0

      MODULE LagrangeInterpolationPolynomial
      IMPLICIT NONE

      CONTAINS
      SUBROUTINE AlmostEqual(a,b,var1)
      IMPLICIT NONE
      real*8, INTENT(IN) :: a,b
      LOGICAL, INTENT(OUT)      :: var1    !!!!!!!!!! logical: true and false
      real*8             :: tol

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
      real*8, DIMENSION(0:N), INTENT(IN)     :: Phi       !values at the nodes
      real*8, DIMENSION(0:N,0:N), INTENT(IN) :: D         !polynomial derivative matrix that has been pre-computed
      real*8, DIMENSION(0:N), INTENT(OUT)    :: Phi_prime !derivative of the interpolent at the nodes      
      INTEGER                                       :: i,j
      real*8                                 :: t

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
      real*8, DIMENSION(0:N), INTENT(IN)  :: x
      real*8, DIMENSION(0:N), INTENT(OUT) :: w
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
      real*8, DIMENSION(0:M), INTENT(IN)  :: x,w,f
      real*8, INTENT(IN)                  :: y
      real*8, INTENT(OUT)                 :: p

      INTEGER                                    :: j
      real*8                              :: numerator,denominator,t
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
      real*8, INTENT(IN) :: xpt
      real*8, DIMENSION(0:N), INTENT(IN)  :: x,w
      real*8, DIMENSION(0:N), INTENT(OUT) :: l

      INTEGER :: j,jmatch
      real*8  :: s,t
      LOGICAL :: var1, xMatchesNode

      xMatchesNode = .FALSE.
      jmatch=-1

      DO j = 0,N
         CALL AlmostEqual(xpt,x(j),var1)
         IF (var1) THEN
            l(j)         = 1.0D0
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
           l(j)=0.0D0
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
      real*8, DIMENSION(0:M), INTENT(IN)  :: x,w,f
      real*8, INTENT(IN)                  :: y      
      real*8, INTENT(OUT)                 :: p1
      INTEGER                                    :: i,j
      real*8                              :: numerator,denominator,t,p
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
      real*8, DIMENSION(0:N), INTENT(IN)  :: x,w
      real*8, DIMENSION(0:N,0:N), INTENT(OUT) :: DL
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


      real*8, dimension(:,:), allocatable :: cache_gauss 
      real*8, dimension(:,:), allocatable :: cache_gauss_lobatto 
      real*8, dimension(:,:), allocatable :: cache_gauss_w 
      real*8, dimension(:,:), allocatable :: cache_gauss_lobatto_w 

      real*8 , dimension(:,:,:), allocatable :: cache_wMATGL
      real*8 , dimension(:,:,:), allocatable :: cache_wMAT
      real*8 , dimension(:,:,:), allocatable :: cache_w_right
      real*8 , dimension(:,:,:), allocatable :: cache_w_left

      real*8 , dimension(:,:,:), allocatable :: cache_wMAT_extend
      real*8 , dimension(:,:,:), allocatable :: cache_wRT
      real*8 , dimension(:,:,:), allocatable :: cache_wLT
      real*8 , dimension(:,:,:), allocatable :: cache_wLRT

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
      real*8, INTENT(IN)   :: x
      real*8, INTENT(OUT)  :: P,dl
      INTEGER                     :: i,j
      real*8               :: P1,P2,dl1,dl2

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
      real*8, DIMENSION(0:n), INTENT(OUT)  :: y,weight
      INTEGER                     :: j,k
      real*8               :: P1,delta
      real*8               :: dP1,localPI

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
      real*8, INTENT(IN)   :: x
      real*8, INTENT(OUT)  :: q,dq,ln
      INTEGER                     :: k
      real*8               :: l1,l2,l3,dl1,dl2,dl,dl3
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




!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 4
! Algorithm 25
! Legendre Gauss-Lobatto
      SUBROUTINE LegendreGaussLobattoNodesAndWeights(n,y,weight)
      !y=Legendre Gauss Lobatto Nodes
      !weight is Legendre Gauss Lobatto weights
      !return y(0:n),weight(0:n);
      INTEGER,INTENT(IN)          :: n
      real*8, DIMENSION(0:n), INTENT(OUT)  :: y,weight
      INTEGER                     :: j,k
      real*8               :: q,dq,ln,delta,localPI

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

      integer i,order_r
      real*8 w(0:order_r)
      real*8 data(0:order_r)
      real*8 sum

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
      real*8 xtarget,lag
      real*8 x(0:order_r)
      real*8 w(0:order_r)

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
      real*8, DIMENSION(0:n), INTENT(OUT)  :: y,weight
      INTEGER                     :: j,k
      real*8               :: localPI
      real*8 w(0:n)
      real*8 yL(0:n)
      real*8 wL(0:n)

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
      real*8, DIMENSION(0:n), INTENT(OUT)  :: y,weight
      INTEGER                     :: j,k
      real*8               :: localPI
      real*8 w(0:n)
      real*8 yL(0:n)
      real*8 wL(0:n)

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
      real*8 x(0:order_r)
      real*8 w(0:order_r,0:order_r)
      real*8 sum

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
      real*8 y1(0:order_r1)
      real*8 y2(0:order_r2)
      real*8 data1(0:order_r1)
      real*8 w(0:order_r1)
      real*8 data2(0:order_r2)

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
      real*8 y1(0:order_r1)
      real*8 y2(0:order_r2)
      real*8 data(0:order_r1)
      real*8 datader1(0:order_r1)
      real*8 datader(0:order_r2)
      real*8 wMAT(0:order_r1,0:order_r1)
      integer i,j
      real*8 sum,dx_element
      real*8 w(0:order_r1)
   
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

      return
      end subroutine delete_cache

       ! typ=0  Legendre
       ! typ=1  Clenshaw Curtis 
       ! r is the largest "order" aka number of points that
       ! one can prescribe.
      subroutine init_cache(order_r,typ)
      IMPLICIT NONE

      integer order_r
      real*8 yleft(0:order_r)
      real*8 yright(0:order_r)
      real*8 y(0:order_r)
      real*8 yGL(0:order_r)
      real*8 y_extend(0:order_r+1)
      real*8 yGL_extend(0:order_r+1)
      real*8 yLT(0:order_r+1)
      real*8 yRT(0:order_r+1)
      real*8 yLRT(0:order_r+1)

      real*8, dimension(:,:), allocatable :: deriv_matrix
      real*8, dimension(:), allocatable :: tempx,tempw

      integer i,j,typ,i1,j1

      if (order_r.lt.1) then
       print *,"order_r invalid"
       stop
      endif

      allocate(cache_gauss(1:order_r+1,0:order_r+1))
      allocate(cache_gauss_lobatto(1:order_r+1,0:order_r+1))
      allocate(cache_gauss_w(1:order_r+1,0:order_r+1))
      allocate(cache_gauss_lobatto_w(1:order_r+1,0:order_r+1))

      allocate(cache_wMATGL(1:order_r,0:order_r,0:order_r))
      allocate(cache_wMAT(1:order_r,0:order_r,0:order_r))
      allocate(cache_w_left(1:order_r,0:order_r,0:order_r))
      allocate(cache_w_right(1:order_r,0:order_r,0:order_r))

      allocate(cache_wMAT_extend(1:order_r,0:order_r+1,0:order_r+1))
      allocate(cache_wRT(1:order_r,0:order_r+1,0:order_r+1))
      allocate(cache_wLT(1:order_r,0:order_r+1,0:order_r+1))
      allocate(cache_wLRT(1:order_r,0:order_r+1,0:order_r+1))

      do i=1,order_r+1

       allocate(tempx(0:i))
       allocate(tempw(0:i))

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

       deallocate(tempx,tempw)

      enddo ! i=1.. order_r+1

      do i=1,order_r

       do i1=0,i-1
        y(i1)=cache_gauss(i,i1)
       enddo
       do i1=0,i
        yGL(i1)=cache_gauss_lobatto(i,i1)
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
        yGL_extend(i1)=cache_gauss_lobatto(i+1,i1)
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

       allocate(deriv_matrix(0:i+1,0:i+1))
       call polyinterp_Dmatrix(i+1,yGL_extend,deriv_matrix)
       do i1=0,i+1
       do j1=0,i+1
        cache_wMAT_extend(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i+1,0:i+1))
       call polyinterp_Dmatrix(i+1,yRT,deriv_matrix)
       do i1=0,i+1
       do j1=0,i+1
        cache_wRT(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i+1,0:i+1))
       call polyinterp_Dmatrix(i+1,yLT,deriv_matrix)
       do i1=0,i+1
       do j1=0,i+1
        cache_wLT(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i+1,0:i+1))
       call polyinterp_Dmatrix(i+1,yLRT,deriv_matrix)
       do i1=0,i+1
       do j1=0,i+1
        cache_wLRT(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)


       allocate(deriv_matrix(0:i,0:i))
       call polyinterp_Dmatrix(i,yGL,deriv_matrix)
       do i1=0,i
       do j1=0,i
        cache_wMATGL(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i-1,0:i-1))
       call polyinterp_Dmatrix(i-1,y,deriv_matrix)
       do i1=0,i-1
       do j1=0,i-1
        cache_wMAT(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i,0:i))
       call polyinterp_Dmatrix(i,yright,deriv_matrix)
       do i1=0,i
       do j1=0,i
        cache_w_right(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

       allocate(deriv_matrix(0:i,0:i))
       call polyinterp_Dmatrix(i,yleft,deriv_matrix)
       do i1=0,i
       do j1=0,i
        cache_w_left(i,i1,j1)=deriv_matrix(i1,j1)
       enddo
       enddo
       deallocate(deriv_matrix)

      enddo ! i

      return
      end subroutine init_cache
 
      subroutine sanity_check(rend)
      IMPLICIT NONE
 
      integer typ,order_r,i,j,k,p,rstart,rend
      real*8, dimension(:), allocatable :: y,weight,data
      real*8, dimension(:), allocatable :: w,datader
      real*8, dimension(:,:), allocatable :: wMAT
      real*8 exact,sum,xtarget,dx_element

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

      PROGRAM spectral_lab
      use LegendreNodes
      IMPLICIT NONE

      integer order_r,typ,rend,i,j,nplot
      real*8, dimension(:), allocatable :: y,weight
      real*8 yextend(0:1000)
      real*8 hplot,prod,xplot

      order_r=18
      typ=0
      call init_cache(order_r,typ)
      rend=18
      call sanity_check(rend)

      typ=2
      order_r=20
      allocate(y(0:order_r))
      allocate(weight(0:order_r))
      if (typ.eq.0) then
       call LegendreGaussNodesAndWeights(order_r,y,weight)
      else if (typ.eq.2) then
       call ClenshawGaussNodesAndWeights(order_r,y,weight)
      endif

      do i=0,order_r
       print *,"i,y,weight ",i,y(i),weight(i)
      enddo

      yextend(0)=-1.0d0-abs(-1.0d0-y(0))
      yextend(order_r+2)=1.0d0+abs(1.0d0-y(order_r))
      do i=0,order_r
       yextend(i+1)=y(i)
      enddo

      do i=0,order_r+2
       print *,"i,yextend ",i,yextend(i)
      enddo
   
      open(unit=4,file="interp_error1")
      open(unit=5,file="interp_error2")
      nplot=1000
      hplot=(yextend(order_r+2)-yextend(0))/nplot 
      do j=0,nplot
       xplot=yextend(0)+j*hplot

       prod=1.0d0
       do i=0,order_r
        prod=prod*(xplot-y(i))
       enddo
       write(4,*) xplot,prod

       prod=1.0d0
       do i=0,order_r+2
        prod=prod*(xplot-yextend(i))
       enddo
       write(5,*) xplot,prod
      enddo
      close(4)
      close(5)

      print *,"file1 `interp_error1'"
      print *,"file2 `interp_error2'"

      deallocate(y)
      deallocate(weight)

      END PROGRAM
