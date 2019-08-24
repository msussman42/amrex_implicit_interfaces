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
!	05 Jan 2005 v4. Changing to CGS and adding 'eheight' as passable variable.
!	06 Jan	v5	Removing extra variables that are not used.
!	09 Jan  v6  Calls Mark's sample xarray, yarray instead of generating it
!	16 Jan  v7  Updates sizeLevel as a function of step number

PROGRAM CLSVOF
  use MAT_Greens,	only : BoundaryForceCalc ! subroutines
  use SortNodePts, only : SortPts

  ! PARAMETER
  integer, parameter	:: XYArrayunit=10

  real*8				:: domlo1,domhi1,eheight

! outgoing variable

  real*8, dimension(:), allocatable		:: qE, xarray, yarray, xarrayNew, yarrayNew

! routine local variables
  integer								:: i, j, sizeLevel
  real*8								:: t1, t2, tempx, tempy, timeCount

continue
  call cpu_time(t1)

  domlo1=0.0	! minX, in cm
  domhi1=1.5	! maxX, in cm
  eheight=1.5	! maxY, in cm

! determines the size of the file
sizeLevel=1	
open(unit=XYArrayunit, file='crashFracDist-4.dat', status='old', action='read', err=990)
do
  read(unit=XYArrayunit, fmt='(2e1.2)', err=994) tempx, tempy
  sizeLevel=sizeLevel+1 ! number of boundary points
enddo
994 close (unit=XYArrayunit, status='keep', err=993)
sizeLevel=sizeLevel-1

! reads input xarray, yarray values
allocate(xarray(sizeLevel)); allocate(yarray(sizeLevel))
allocate(xarrayNew(sizeLevel)); allocate(yarrayNew(sizeLevel))

open(unit=XYArrayunit, file='crashFracDist-4.dat', status='old', action='read', err=990)
do i =1,sizeLevel
  read(unit=XYArrayunit, fmt='(e8.2,e11.2)', err=991) xarray(i), yarray(i)
enddo
close (unit=XYArrayunit, status='keep', err=993)

! sorts xarray, yarray into an nearest-neighbor listing from lower left corner
call SortPts(xarrayNew, yarrayNew, sizeLevel, xarray, yarray )

! using boundary integeral method, calculates electric field at given locations
allocate(qE(sizeLevel))

!yleft,yright, domlo1, domhi1, eheight, and the LS function
! MAT_Greens17: domlo1,domhi1,eheight,sizeLevel -  0.0, 1.5, 1.5, 74
call BoundaryForceCalc( domlo1, domhi1, eheight, xarrayNew, yarrayNew, sizeLevel, qE ) 
write(*, fmt=11) 'Done Step 1. Qe(1)= ', qE(1)

deallocate(xarray); deallocate(yarray)
deallocate(xarrayNew); deallocate(yarrayNew)

goto 999

!------------
! FORMAT listings
10 FORMAT(a,f6.2,a)
11 FORMAT(a,e15.6)

! ERROR listings
990 print *, 'Error open xyarray'
991 print *, 'Error read xyarray'
992 print *, 'Error write xyarray'
993 print *, 'Error close xyarray'

999   call cpu_time(t2); timeCount=t2-t1
write(*, fmt=10) 'UMich program call took ', timeCount, ' seconds.'
deallocate(qE)

end PROGRAM CLSVOF
