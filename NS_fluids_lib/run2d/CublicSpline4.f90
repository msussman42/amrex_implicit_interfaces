MODULE CubicSpline4

  use SortNodePts3, only : GaussQuad_data

implicit none
! used to restrict variables to prevent accidental calling outside
private      ! everything defaults to 'private' except those expressly labeled otherwise.
public								:: CurveFit, NewCoords
real*8, parameter					:: epsil=1e-6
real*8, dimension(:,:), allocatable :: NewCoords

contains 

! ***************
! Subroutine CurveFit
! ***************
! DESCRIPTION: given N points and S(1):S(N) and S'(1)=S'(N)=0, determines the linear or 
!   cubic spline fit between those points so the N equations match on x, x' and x'' 
!
! INPUT: knots, panel lengths
!
! OUTPUT: C[x] coordinate 
!
! CALLING PROGRAM: MAT_Greens
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   03 Dec 2004 begun
!	07 Dec	v4	changing to calculate the total number of 

SUBROUTINE CurveFit ( SplinePts, knotsx, knotsy, hX, hY, sizeknot, &
  PanelVals, ConeSize, CublicSplineFit )

! incoming variables
integer, intent(in)							:: sizeknot, ConeSize, CublicSplineFit
real*8, dimension(sizeknot), intent(in)		:: knotsx, knotsy, hX, hY
real*8, dimension(sizeknot,2), intent(in)	:: PanelVals

! outsourced variables
integer, intent(out)						:: SplinePts

! subroutine entirely variables
logical										:: XYFit
integer										:: i, j
real*8, dimension(sizeknot)					:: totDist
real*8, dimension(ConeSize,4)				:: Sx, Sy

continue ! body of program
!*******
! the final equation will look like Si[x]=ai+bi(x-xi)+ci(x-xi)^2+di(x-xi)^3
! setting up matricies and constants for 
! CubicSplineFit = 1 for cubic interpolation, = 0 for linear interpolation

! finds the spline coefficient fits (linear or cubic) for the knotsx/knotsy points 
!   using x vs. y. Use only Sx
!   NOTE: 'totDist' is meaningless in the first iteration
call LinearCubic ( Sx, Sy, knotsx, knotsy, hX, hY, sizeknot, ConeSize, totDist, CublicSplineFit, 1 )

! using Sx, calculates the total distance around the shape and the corresponding
!   number of total panels for a balance between accuracy and speed. Use totDist, SplinePts
call totDistCalc ( totDist, SplinePts, Sx, knotsx, knotsy, ConeSize, sizeknot, 10 )

! using totDist, finds the spline coefficient fits (linear or cubic) for the 
!   knotsx/knotsy points using S vs. x and S vs. y. Use new Sx, Sy
call LinearCubic ( Sx, Sy, knotsx, knotsy, hX, hY, sizeknot, ConeSize, totDist, CublicSplineFit, 0 )

! using Sx, Sy for the new curve, returns the evenly space points SplineCoords
!   along the surface, Finally(!) get 'SplineCoords' we desire
allocate(NewCoords(SplinePts, 4))
call SplineLengthFit (NewCoords, Sx, Sy, knotsx, knotsy, SplinePts, PanelVals, sizeknot, &
  totDist, ConeSize, CublicSplineFit)

END SUBROUTINE CurveFit

! ***************
! Subroutine LinearCubic
! ***************
! DESCRIPTION: given N points and S(1):S(N) and S'(1)=S'(N)=0, determines the linear or
!   cubic spline fit between those points so the N equations match on x, x' and x'' 
!
! INPUT: knots, end derivatives, panel lengths
!
! OUTPUT: C[x] coordinate 
!
! CALLING PROGRAM: MAT_Greens
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   14 Nov 2004 begun
!   27 Nov		changing distance calculation to go along spline, not linear
!	02 Dec		doing an x-y spline fit, then calculating distance afterwards
!   03 Dec		pulling linear and cubic into 1 subroutine linearcubic

SUBROUTINE LinearCubic ( Sx, Sy, knotsx, knotsy, hX, hY, sizeknot, ConeSize, totDist, &
  CubicSplineFit, XYFit )

! incoming variables
integer, intent(in)							:: sizeknot, ConeSize, CubicSplineFit, XYFit
real*8, dimension(sizeknot), intent(in)		:: knotsx, knotsy, hX, hY, totDist

! outsourced variables
real*8, dimension(ConeSize,4), intent(out)  :: Sx, Sy

! subroutine entirely variables
integer										:: i, j, CubicXYFit
real*8, dimension(ConeSize)					:: a, b, c, d, RHSCubic, aDiag, bDiag, &
												cDiag
real*8, dimension(sizeknot)					:: hDist, distChange

continue ! body of program
!*******
! the final equation will look like Si[x]=ai+bi(x-xi)+ci(x-xi)^2+di(x-xi)^3
! setting up matricies and constants for 
CubicXYFit = 10*XYFit + CubicSplineFit

select case (CubicXYFit)
case (11)				! First loop, x vs y, cubic interpolation
!-------------------------------------------------------------------
  a(1)=knotsy(1); a(ConeSize)=knotsy(ConeSize)
  aDiag(1)=0; bDiag(1) = 1; cDiag(1)=0
  aDiag(ConeSize) = 0; bDiag(ConeSize)=1; cDiag(ConeSize)=0

  do i=2,ConeSize-1
    a(i)=knotsy(i)
    aDiag(i)=hX(i-1)/3.
    bDiag(i)=2.*hX(i-1)/3. + 2.*hX(i)/3.
    cDiag(i)=hX(i)/3.
  enddo

  do i=2,ConeSize-1
    if ( abs(hX(i)) > epsil .AND. a(i+1)-a(i) > epsil/10 ) then
      RHSCubic(i)=1.*(a(i+1)-a(i))/hX(i) - (a(i)-a(i-1))/hX(i-1)
    elseif ( abs(hX(i)) <= epsil .OR. a(i+1)-a(i) <= epsil/10) then
      RHSCubic(i)=0
    else
      stop 'wrong RHSCubic'
    endif
  enddo

  !RHSCubic(1)=(a(2)-a(1))/hX(1)
  RHSCubic(1)=0; RHSCubic(ConeSize)=0

  call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

  do i=1,ConeSize-1
    b(i) = ( a(i+1)-a(i) )/hX(i) - (2.*hX(i)*c(i)/3. + hX(i)*c(i+1)/3.)
    d(i) = ( c(i+1)-c(i) )/(3.*hX(i))
  enddo

  Sx(:,1)=a(:); Sx(:,2)=b(:); Sx(:,3)=c(:); Sx(:,4)=d(:)

case (10)				! First loop, x vs y, linear interpolation
!-------------------------------------------------------------------
  Sx(:,3)=0; Sx(:,4)=0

  do i=1,ConeSize
    Sx(i,1)=knotsy(i)
  enddo
  do i=1,ConeSize-1
    Sx(i,2) = ( Sx(i+1,1)-Sx(i,1) )/hX(i)
  enddo

case (01)				! Second loop, x,y vs s, cubic interpolation
!-------------------------------------------------------------------
  distChange(1) = totDist(1)
  a(1)=knotsx(1); a(ConeSize)=knotsx(ConeSize)
  aDiag(1)=0; bDiag(1) = 1; cDiag(1)=0
  aDiag(ConeSize) = 0; bDiag(ConeSize)=1; cDiag(ConeSize)=0

  do i=2,ConeSize-1
    distChange(i) = totDist(i)-totDist(i-1)
    a(i)=knotsx(i)
    aDiag(i)=distChange(i-1)/3.
    bDiag(i)=2.*distChange(i-1)/3. + 2.*distChange(i)/3.
    cDiag(i)=distChange(i)/3.
  enddo

  do i=2,ConeSize-1
    if ( abs(distChange(i)) > epsil .AND. a(i+1)-a(i) > epsil/10 ) then
      RHSCubic(i)=1.*(a(i+1)-a(i))/distChange(i) - (a(i)-a(i-1))/distChange(i-1)
    elseif ( abs(distChange(i)) <= epsil .OR. a(i+1)-a(i) <= epsil/10) then
      RHSCubic(i)=0
    else
      stop 'wrong RHSCubic'
    endif
  enddo

  !RHSCubic(1)=(a(2)-a(1))/hX(1)
  RHSCubic(1)=0; RHSCubic(ConeSize)=0

  call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

  do i=1,ConeSize-1
    b(i) = ( a(i+1)-a(i) )/distChange(i) - (2.*distChange(i)*c(i)/3. + distChange(i)*c(i+1)/3.)
    d(i) = ( c(i+1)-c(i) )/(3.*distChange(i))
  enddo

  Sx(:,1)=a(:); Sx(:,2)=b(:); Sx(:,3)=c(:); Sx(:,4)=d(:)
  !******
  ! Same section, now for y vs s.
  do i=1,ConeSize
    a(i)=knotsy(i)
  enddo

  do i=2,ConeSize-1
    if ( abs(distChange(i)) > epsil .AND. a(i+1)-a(i) > epsil/10 ) then
      RHSCubic(i)=1.*(a(i+1)-a(i))/distChange(i) - (a(i)-a(i-1))/distChange(i-1)
    elseif ( abs(distChange(i)) <= epsil .OR. a(i+1)-a(i) <= epsil/10) then
      RHSCubic(i)=0
    else
      stop 'wrong RHSCubic'
    endif
  enddo

  !RHSCubic(1)=(a(2)-a(1))/hY(1)
  RHSCubic(1)=0; RHSCubic(ConeSize)=0

  call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

  do i=1,ConeSize-1
    b(i) = ( a(i+1)-a(i) )/distChange(i) - (2.*distChange(i)*c(i)/3. + distChange(i)*c(i+1)/3.)
    d(i) = ( c(i+1)-c(i) )/(3.*distChange(i))
  enddo

  Sy(:,1)=a(:); Sy(:,2)=b(:); Sy(:,3)=c(:); Sy(:,4)=d(:)

case (00)				! Second loop, x,y vs s, linear interpolation
!-------------------------------------------------------------------
  Sx(:,3)=0; Sx(:,4)=0
  Sy(:,3)=0; Sy(:,4)=0

  do i=1,ConeSize
    if (i /= 1) then
  	  distChange(i) = totDist(i)-totDist(i-1)
	else
  	  distChange(i) = totDist(i)
	endif
	Sx(i,1)=knotsx(i)
	Sy(i,1)=knotsy(i)
  enddo

  do i=1,ConeSize-1
	Sx(i,2) = ( Sx(i+1,1)-Sx(i,1) )/distChange(i)
	Sy(i,2) = ( Sy(i+1,1)-Sy(i,1) )/distChange(i)
  enddo
case default
  stop 'Error on CubicXYFit'
end select ! switch CubicXYFit

END SUBROUTINE LinearCubic

!----------
! ***************
! Subroutine SplineLengthFit
! ***************
! DESCRIPTION: takes linear and cubic line fits and places points along them evenly
!
! INPUT: Sx, Sy
!
! OUTPUT: SplineCoords
!
! CALLING PROGRAM: Linear, Cubic
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   21 Nov 2004 begun

SUBROUTINE SplineLengthFit (SplineCoords, Sx, Sy, knotsx, knotsy, SplinePts, &
  PanelVals, sizeknot, totDist, ConeSize, CubicSplineFit )

! incoming variables
integer, intent(in)							:: SplinePts, sizeknot, ConeSize, CubicSplineFit
! CubicSplineFit = 1 for cubic interpolation, = 0 for linear interpolation
real*8, dimension(ConeSize,4), intent(in)	:: Sx, Sy
real*8, dimension(sizeknot,2), intent(in)	:: PanelVals
real*8, dimension(sizeknot), intent(in)		:: totDist, knotsx, knotsy

! outsourced variables
real*8, dimension(SplinePts,4), intent(out) :: SplineCoords

! subroutine entirely variables
integer										:: Splineunit=10, tempJ, i, j
real*8										:: temptotDist, tempfracDist

real*8, dimension(SplinePts)				:: SplineLength

continue

tempJ=1; SplineLength(:)=0
do i=1,SplinePts
  SplineLength(i) = 1.*(i-1)*totDist(sizeknot)/SplinePts
  if (tempJ /= 1) then
    j = tempJ-1
  else
    j = 1
  endif
  
  ChoosePanel: do while ( totDist(j)+epsil <= SplineLength(i) .AND. j<=sizeknot ) 
    j=j+1
  enddo ChoosePanel
  tempJ = j

  if (1==j) then
    temptotDist = 0
  else
    temptotDist = totDist(j-1)
  endif

  if (j<ConeSize) then
    tempfracDist = SplineLength(i)-temptotDist
    SplineCoords(i,1) = Sx(j,1) + Sx(j,2)*(tempfracDist) + &
	  Sx(j,3)*tempfracDist**2 + Sx(j,4)*tempfracDist**3
    SplineCoords(i,2) = Sy(j,1) + Sy(j,2)*(tempfracDist) + &
      Sy(j,3)*(tempfracDist)**2 + Sy(j,4)*(tempfracDist)**3
  else
    tempfracDist = ( SplineLength(i)-totDist(j-1)) / (totDist(j)-totDist(j-1))
	if (j /=sizeknot) then
  	  SplineCoords(i,1) = tempfracDist*(knotsx(j+1)-knotsx(j))+knotsx(j)
      SplineCoords(i,2) = tempfracDist*(knotsy(j+1)-knotsy(j))+knotsy(j)
	else
  	  SplineCoords(i,1) = tempfracDist*(knotsx(1)-knotsx(j))+knotsx(j)
      SplineCoords(i,2) = tempfracDist*(knotsy(1)-knotsy(j))+knotsy(j)
	endif
  endif
  SplineCoords(i,3) = PanelVals(j,1)
  SplineCoords(i,4) = PanelVals(j,2)
enddo

  open(unit=Splineunit, file='SplinePts.dat', status='replace')
  if (CubicSplineFit==1) then
    write(unit=Splineunit, fmt='(a,i5,a)') 'TITLE="Boundary shape using cubic interpolation and ', SplinePts, ' panels"'
  elseif (CubicSplineFit==0) then
    write(unit=Splineunit, fmt='(a,i5,a)') 'TITLE="Boundary shape using linear interpolation and ', SplinePts, ' panels"'
  else
    stop 'Spline unit label error'
  endif
  write(unit=Splineunit, fmt='(a)') 'VARIABLES="X"'
  write(unit=Splineunit, fmt='(a)') '"Y"'
  write(unit=Splineunit, fmt='(a)') '"Spline Num"'
  write(unit=Splineunit, fmt='(a)') 'ZONE T="main"'
  do i=1,SplinePts
    write(unit=Splineunit, fmt='(2e20.6, i8)') SplineCoords(i,1), SplineCoords(i,2), i
  enddo
  close (unit=Splineunit, status='keep')

END SUBROUTINE SplineLengthFit

!----------
! ***************
! Subroutine Tridiag
! ***************
! DESCRIPTION: given a tridiagonal matrix, rapidly solves for Ax=b
!
! INPUT: three bands, rhs
!
! OUTPUT: x
!
! CALLING PROGRAM: Cubic
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   14 Nov 2004 begun

SUBROUTINE Tridiag (u, leftrow, centerrow, rightrow, rhsTri, N)

! incoming variables
integer, intent(in)					:: N
real*8, dimension(N), intent(in)	:: leftrow, centerrow, rightrow, rhsTri

! outsourced variables
real*8, dimension(N), intent(out)	:: u

! subroutine entirely variables
integer								:: j
real*8								:: bet
real*8, dimension(N)				:: gam

continue ! body of program
!*******
if (0.0==centerrow(1)) then
  stop 'Error 1 in tridiag'
endif

bet=centerrow(1)
u(1)=rhsTri(1)/bet

do j=2,N		! decomposition and forward substitution
  gam(j) = rightrow(j-1)/bet
  bet=centerrow(j)-leftrow(j)*gam(j)
  if (0.0==bet) then
    stop 'Error 2 in tridiag'
  endif
  u(j)=( rhsTri(j)-leftrow(j)*u(j-1) )/bet
enddo

do j=N-1,1,-1	! back substitution
  u(j) = u(j) - gam(j+1)*u(j+1)
enddo

END SUBROUTINE Tridiag

!----------
! ***************
! Subroutine totDistCalc
! ***************
! DESCRIPTION: given a spline fitting coefficients, numerically integrates
!   along surface to get an accurate total distance around shape
!
! INPUT: knot locations, cubic spline coeffs, # of knots along test shape and 
!	overall
!
! OUTPUT: total distance around entire border, calced along spline arclength
!
! CALLING PROGRAM: Cubic
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   27 Nov 2004 begun

SUBROUTINE totDistCalc ( totDist1, SplinePts1, Sx1, knotsx1, knotsy1, ConeSize1, sizeknot1, M1 )

! incoming variables
integer, intent(in)							:: ConeSize1, sizeknot1, M1
real*8, dimension(ConeSize1,4), intent(in)	:: Sx1
real*8, dimension(sizeknot1), intent(in)	:: knotsx1, knotsy1

! outsourced variables
integer, intent(out)						:: SplinePts1
real*8, dimension(sizeknot1), intent(out)	:: totDist1

! subroutine entirely variables
integer										:: i,j, panel, totPan, ConePan
real*8										:: tempXDist, tempArcDist, tempXChange, &
												tempYChange
real*8, dimension(sizeknot1)				:: xPanLen, yPanLen
real*8, dimension(sizeknot1*M1)				:: xPan, yPan, dydx, PanGaussWeight

continue
 
totPan = sizeknot1*M1; 
ConePan = ConeSize1*M1; tempArcDist=0.

! dePanel = floor(1.*(i-1)/M)+1
do i = 1, sizeknot1
  if (i /= sizeknot1) then
    xPanLen(i) = knotsx1(i+1)-knotsx1(i)
	yPanLen(i) = knotsy1(i+1)-knotsy1(i)
  elseif (i == sizeknot1) then
    xPanLen(i) = knotsx1(1)-knotsx1(i)
	yPanLen(i) = knotsy1(1)-knotsy1(i)
  else
    stop 'xPanLen error'
  endif

  do j=1, M1
    panel= (i-1)*M1+j
    if (panel<=ConePan) then
      call GaussQuad_data ( xPan(panel), yPan(panel), PanGaussWeight(panel), M1, &
   	    knotsx1(i), xPanLen(i), knotsy1(i), tempArcDist, j )

	  tempXDist = xPan(panel)-knotsx1(i)

      yPan(panel) = Sx1(i,1) + Sx1(i,2)*(tempXDist) + Sx1(i,3)*(tempXDist)**2 + &
  	    Sx1(i,4)*(tempXDist)**3   
    else
      call GaussQuad_data ( xPan(panel), yPan(panel), PanGaussWeight(panel), M1, &
   	    knotsx1(i), xPanLen(i), knotsy1(i), yPanLen(i), j )
    endif
  enddo
enddo

do i = 1, sizeknot1
  tempArcDist = 0.
  do j=1, M1
    panel= (i-1)*M1+j

    if (panel /= totPan) then
      tempXChange = xPan(panel+1)-xPan(panel); tempYChange = yPan(panel+1)-yPan(panel)
    else
      tempXChange = xPan(1)-xPan(panel); tempYChange = yPan(1)-yPan(panel)  
    endif
    if ( abs(tempXChange) < epsil ) then
      dydx(panel) = 0
    else
      dydx(panel) = tempYChange / tempXChange
    endif

    if ( xPanLen(i) /= 0 ) then
  	  tempArcDist = tempArcDist + abs( xPanLen(i) )*0.5*PanGaussWeight(panel)*( (1+dydx(panel)**2) )**0.5
	else
	  tempArcDist = tempArcDist + abs( yPanLen(i) )*0.5*PanGaussWeight(panel)
	endif
  enddo

  if (i /= 1) then
    totDist1(i) = totDist1(i-1) + tempArcDist
  else
    totDist1(i) = tempArcDist
  endif
enddo

! now calculate the optimal number of total panels for the calculation. For a balance between
!   speed and accuracy, it should be large enough so there are between 0.5 and 2x the number 
!   of spline points in totDist(ConeSize1) as ConeSize1 ( i.e. 0.5*ConeSize1 < SplinePts1 < 2*ConeSize1 )
SplinePts1 = 0.4*ConeSize1*(totDist1(sizeknot1)/totDist1(ConeSize1))

END SUBROUTINE totDistCalc

END MODULE CubicSpline4
