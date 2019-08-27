! ***************
! Program CLSVOF[x]
! ***************
! DESCRIPTION: main calling program to test call my boundary integral subroutines
!
! INPUT: <na>
!
! OUTPUT: Potential field at a given set of points
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   17 Sep 2004 begun
!	07 Dec	v2	Passing sigma values back and forth to reduce GMRes times.
!					Also, determining xarray, yarray directly through global vars
!	08 Dec	v3	Seems not to make a timing difference; reverting to before.
!					Clearing out coding and print statements

PROGRAM CLSVOF3
  use MAT_Greens16,	only : MAT_Greens ! subroutines
  use Wrapper6,		only : Intersection_Points, xarray, yarray

  ! PARAMETER
  real, parameter			:: pi=3.14159

  integer					:: gridlo1, gridlo2, gridhi1, gridhi2, & 
								llo1, llo2, lhi1, lhi2, sizeLevel
  real*8					:: dx, dy, domlo1, domlo2, domhi1, domhi2
  real*8, dimension(:,:), allocatable	:: level

! outgoing variable
  real*8, dimension(:), allocatable		:: qE

! routine local variables
  integer								:: i, j
  real									:: randomNum
  real, dimension(:), allocatable		:: x, y
  real*8, dimension(:), allocatable		:: xarray2, yarray2


continue

  domlo1=0.0
  domhi1=1.1
  domlo2=0.0
  domhi2=1.0
  dx=1./128.; dy=1./128.
  gridlo1=0
  gridlo2=0
  gridhi1=127
  gridhi2=127
  llo1=-1
  llo2=-1   
  lhi1=126   
  lhi2=126   

  allocate(x(llo1:lhi1)); allocate(y(llo2:lhi2))
  allocate(level(llo1:lhi1, llo2:lhi2))

  do i=llo1,lhi1
    x(i)=(i+0.5)*dx + domlo1
    do j=llo2,lhi2
	  y(j)=(j+0.5)*dy + domlo2
      level(i,j)= 0.5 + 0.25*cos ( 2*pi*(x(i)-0.5) )-y(j)
    end do
  end do
deallocate(x); deallocate(y)

! determines the number of intersections and their locations. Builds xarray,
!   yarray and sizeLevel
call Intersection_Points( llo1, lhi1, llo2, lhi2, level, dx, dy, &
  sizeLevel, domlo1, domlo2 )

! using boundary integeral method, calculates electric field at given locations
allocate(qE(sizeLevel))

call MAT_Greens( level, llo1, llo2, lhi1, lhi2, gridlo1, gridlo2, &
  gridhi1, gridhi2, dx, dy, domlo1, domlo2, domhi1, domhi2, xarray, &
  yarray, sizeLevel, qE ) 

goto 999
! FORMAT listings
10 FORMAT(3f12.3)
11 FORMAT(2f10.3, e12.3)

! ERROR listings
990 print *, 'Error open sigmas'
991 print *, 'Error read sigmas'
992 print *, 'Error write sigmas'
993 print *, 'Error close sigmas'

999 deallocate(level); deallocate(qE)

end PROGRAM CLSVOF3
