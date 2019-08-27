MODULE GMRes_Andrew
  !v2	21 Nov	implementing pointers
  !v3	24 Nov	removing pointers, changing subroutine minimize
  !v4	29 Nov	changing to Jerry's newest, removing 'minimize'
  !v5	06 Dec	v4 working! Changing to single row allocatable arrays
  !				  to ease stack handling

implicit none

public	:: gmresIter
contains

!-------------------------------------------------------
! Title:        Iterative matrix solver (GMRES) 
! Purpose:      Solve Ax=B
! By:           Jerry Emhoff
! Rewritten by: Anton VanderWyst
!
!  A=a[num][num]
!  b=b[num]
!  x=x[num]
!
! INPUT: influence matrix A, RHS b, initial matrix guess x (can be zero)
!
! OUTPUT: new x
!
! CALLING PROGRAM: MAT_Greens

subroutine gmresIter(a, b, x, num) 

! incoming variables
  integer, intent(in)					:: num    ! size of array
  real*8, dimension(num), intent(in)	:: b
  real*8, dimension(num*num),intent(in) :: a

! outgoing variable
  real*8, dimension(num), intent(inout)	:: x

! program variables
  integer							:: go ! formerly logical
  integer							:: i, j, n, n1, numplus
    
  real*8							:: temp, errlim, bnorm

  real*8, dimension(:), allocatable	:: v, y, r, b1, qfact, givens, q, h, rfact

!  real*8, dimension(num,num)		:: q	! Q satisfying the equation AQ(n)=Q(n+1)H(n)
!  real*8, dimension(num+1,num)     :: h	! Hessenberg matrix formed by orthogonal transformation of A
!  real*8, dimension(num)			:: v	! Vector used for computation of Aq(n) and other quantities
!  real*8, dimension(num)			:: y	! Vector satisfying the least squares problem ||(Hy-||b||e1)||
!  real*8, dimension(num)			:: r	! Residual vector, which is solved for given an initial guess
!  real*8, dimension(num)			:: b1	! Vector equaling b-Ax(0), so that Ar=b1 can be solved

!  real*8, dimension(num+1)			:: qfact  ! The first row of the Q matrix of H=QR
!  real*8, dimension(num,num)		:: rfact  ! The R matrix of H's QR Decomposition
!  real*8, dimension(2*num)			:: givens ! The Givens Rotation coefficients for each iteration

continue

  allocate(q(num*num)); allocate(h((num+1)*num)); allocate(v(num)); allocate(y(num))
  allocate(r(num)); allocate(b1(num)); allocate(qfact(num+1)); allocate(rfact(num*num))
  allocate(givens(2*num))

  go=1; n=1; n1=2
  temp=0.; errlim=1e-7; bnorm=0.
  q(:)=0
  qfact(:)=0; rfact(:)=0; givens(:)=0

  ! Set the initial values of qfact
  qfact(1)=1.0;
  qfact(2)=1.0;

  ! Set the b1-vector using the initial guess of x, b1=b-Ax(0)
  do i=1,num
    b1(i)=b(i)
	do j=1,num
	  b1(i) = b1(i) - a((i-1)*num+j)*x(j)
	enddo
  enddo

  ! Compute the norm of the b1 vector
  bnorm=0.0;
  do i=1,num
    bnorm = bnorm + b1(i)**2
  enddo
  bnorm = sqrt(bnorm)

  ! Compute q(0)=b1/||b1||
  do i=1,num
    q(0*num+i)=b1(i)/bnorm
  enddo

  GoWhile: do while (1==go) 
    ! Compute v = A*q(n)
	do i=1,num
	  v(i)=0.
	  do j=1,num
	    v(i) = v(i) + a((i-1)*num+j)*q((n-1)*num+j)
	  enddo
	enddo

    do j=1, n1
	  ! Compute h(j,n) = q(j)*v
      h((j-1)*num+n) = 0.
	  do i=1, num
	    h((j-1)*num+n) = h((j-1)*num+n) + q((j-1)*num+i)*v(i)
      enddo

      ! Compute v = v-h(j,n)*q(j)
      do i=1,num
	    v(i) = v(i) - h((j-1)*num+n)*q((j-1)*num+i)
      enddo
    enddo

    ! Set h(n+1,n) = ||v||
    temp=0.0;
	do i=1,num
	  temp = temp + v(i)**2
    enddo
	h((n1-1)*num+n)=sqrt(temp);

    if (n1 < num) then
      ! Compute q(n+1) = v/h(n+1,n)
      do i=1,num
        q((n1-1)*num+i)=v(i)/h((n1-1)*num+n)
      enddo
    endif

    ! Do the QR Factorization using a Givens Rotation
	rfact(n)=h(n)

    ! Apply the previous rotations to the new column of rfact
	do i=1,n-1
	! do i=1,n
	  temp=rfact((i-1)*num+n)
	  rfact((i-1)*num+n) = givens(2*i)*temp-givens(2*i+1)*h(i*num+n)
	  rfact(i*num+n) = givens(2*i+1)*temp+givens(2*i)*h(i*num+n)
    enddo

    ! Compute the Givens Rotation coefficients
    temp=sqrt( rfact((n-1)*num+n)**2 + h((n1-1)*num+n)**2 )
	givens(2*n) = rfact((n-1)*num+n)/temp
    givens(2*n+1) = -h((n1-1)*num+n)/temp

    ! Compute the new diagonal element of rfact
    temp=rfact((n-1)*num+n)
	rfact((n-1)*num+n) = givens(2*n)*temp - givens(2*n+1)*h((n1-1)*num+n)

    ! Compute the top row of Q
	qfact(n1) = qfact(n)*givens(2*i+1)
	qfact(n) = qfact(n)*givens(2*i)

    ! Use back substitution to calculate y from the system Ry=Q*b
    do i=n, 1, -1
	  y(i) = bnorm*qfact(i)
	  do j=i+1, n1
	    y(i) = y(i) - rfact((i-1)*num+j)*y(j)
	  enddo
	  y(i) = y(i)/rfact((i-1)*num+i)
	enddo
   
    ! Calculate x=x+Q*y when converged
	if (n1==num .OR. abs(y(n)) < errlim .AND. abs(y(n)) >1e-20 ) then
      r(:)=0.
	  
	  ! Set x(i) = q(i)*y
	  do j=1,n1
	    do i=1,num
		  r(i) = r(i) + q((j-1)*num+i)*y(j)
		enddo
	  enddo

      do i=1,num
	    x(i) = x(i) + r(i)
	  enddo

      go=0
	  !write(*, fmt=('(a,i4)')) 'Total iterations n= ', n
    endif

    n=n+1; n1=n+1
  enddo GoWhile

  write(*, fmt=10), '   Done after ', n-1, ' GM_Res steps.'

  deallocate(q); deallocate(h); deallocate(v); deallocate(y)
  deallocate(r); deallocate(b1); deallocate(qfact); deallocate(rfact)
  deallocate(givens)

goto 999		! actual program end, just lists formats and errors below
  ! ***********
! FORMAT listings
10  FORMAT(a, i3, a)

999 end SUBROUTINE gmresIter

end MODULE GMRes_Andrew
