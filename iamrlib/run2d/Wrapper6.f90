MODULE Wrapper6
! 03 Jun v1
! 10 Jun v2
! 14 Jun v3
! 19 Sep v4
! 11 Oct v5
! 08 Dec v6 Making 'xarray','yarray' global

  use SortNodePts3, only : TrapInit

implicit none

! used to restrict variables to prevent accidental calling outside
private      ! everything defaults to 'private' except those expressly labeled otherwise.
public					:: Intersection_Points, ForceZero
real, parameter, public :: pi=3.14159
real*8, dimension(:), allocatable, public :: xarray, yarray

contains 

! ***************
! SUBROUTINE Intersection_Points
! ***************
! DESCRIPTION: takes in rectangular grid and levelset values at each and outputs
!   electric fields at the points around wherelevelset changes sign
!
! INPUT:  loX,hiX, loY, hiY  (dimensions of computational domain)
!         levelset (2d array of values for the levelset function)
!         levello,levelhi (dimensions of levelset array)
!         dx,dy (grid size)
!
! OUTPUT: ElectricfieldS (2d array of values for the electric field)
!         elo,ehi (dimensions of electricfield array)
!
! CALLING PROGRAMS: CLSVOF
! 
! ALGORITHMN: 
!   1) finds the c0 (surface) from the levelset through interpolation
!   2) using Green's function, computes the electric field on that surface
!   3) extrapolates the field back onto the original node points
!
! Primary author: Anton VanderWyst, University of Michigan 
!
! Version history:
!   03 Jun 2004 begun
!   09 Jun got rectangular grids and linear c0 interpolation working. 283 lines
!   10 Jun Green's functions for all the points in the electrode file
!   12 Jun plotting point electric fields back to orig. grid. Involves tracking
!            each surface point's original grid location. 
!   14 Jun new approach of constructing the surface as you go. Instead of sweeping
!            through n^2, do the 8n along the normal
!   19 Sep skipping v3, customing for only marching triangles and CLSVOF

SUBROUTINE Intersection_Points( loXS, hiXS, loYS, hiYS, levelsetS, &
	dxS, dyS, c0Count, minDomainX, minDomainY )

! incoming variables
integer, intent(in)                 :: loXS, loYS, hiXS, hiYS
real*8, dimension(loXS:hiXS,loYS:hiYS), intent(in) :: levelsetS
real*8, intent (in)					:: dxS, dyS, minDomainX, minDomainY

! outsourced variables
integer, intent(out) 				:: c0Count

! subroutine entirely variables
integer								:: i,j,k, sizeDomain, GridPos
integer								:: levelL, levelR, levelU, levelD
real*8, dimension( loXS:hiXS )		:: xS
real*8, dimension( loYS:hiYS )		:: yS
real*8, dimension(5000)       		:: c0X, c0Y

real*8								:: xMidPhi, yMidPhi

continue ! body of program
!*******

do i=loXS, hiXS
  xS(i)=minDomainX+(i-loXS)*dxS
end do
do j=loYS, hiYS
  yS(j)=minDomainY+(j-loYS)*dyS
end do

! ALGORITHMN(1) finding c0
Rowc0: do j=loYS,hiYS
  Colc0: do i=loXS,hiXS
    sizeDomain=hiYS*(j-1)+i

	if (1==sizeDomain) then						! LL-1
	  GridPos=1
	else if (hiXS==i .and. 1==j) then			! LR-2
	  GridPos=2
    else if (1==j) then							! bottom-3
      GridPos=3
	else if (hiXS==i .and. j<hiYS) then			! RHS-4
	  GridPos=4
	else if (1==i .and. j<hiYS) then			! LHS-5
	  GridPos=5
	else if (1==i .and. hiYS==j) then			! UL-6
	  GridPos=6 
	else if (hiXS==i .and. hiYS==j) then		! UR-7
	  GridPos=7
    else if (hiYS==j) then 						! top-8
	  GridPos=8
	else										! center blocks-9
      GridPos=9
	end if

	LevelNode: select case (GridPos)
	  case (1)
	    levelD=j; levelU=j+1; levelR=i+1; levelL=i
	  case (2)
	    levelD=j; levelU=j+1; levelR=i; levelL=i-1
	  case (3)
	  	levelD=j; levelU=j+1; levelR=i+1; levelL=i-1
	  case (4)
	  	levelD=j-1; levelU=j+1; levelR=i; levelL=i-1
	  case (5)								
	  	levelD=j-1; levelU=j+1; levelR=i+1; levelL=i
	  case (6)
	  	levelD=j-1; levelU=j; levelR=i+1; levelL=i
	  case (7)
	  	levelD=j-1; levelU=j; levelR=i; levelL=i-1
	  case (8)
	  	levelD=j-1; levelU=j; levelR=i+1; levelL=i-1
	  case (9)
	  	levelD=j-1; levelU=j+1; levelR=i+1; levelL=i-1
	  case default
	    print *, 'Error, no LevelNode case selected.'
		exit
    end select LevelNode

	UPLevel: if (levelsetS(i,j)>0) then	  
	  if (levelsetS(i,levelU)<0 .and. j<hiYS) then			! neg. up
	    c0Count=c0Count+1
	    c0X(c0Count)= xS(i)
		c0Y(c0Count)= yS(j) - levelsetS(i,j)*(yS(levelU)-yS(j))/ &
		  (levelsetS(i,levelU)-levelsetS(i,j))
	  end if

	  if (levelsetS(levelR,levelU)<0 .and. i<hiXS .and. j<hiYS) then ! neg. upper right
	    call TwoD_Phi( xS(levelR),yS(levelU),levelsetS(levelR,levelU), xS(i), &
		  yS(j), levelsetS(i,j), xMidPhi, yMidPhi )

	    c0Count=c0Count+1
	    c0X(c0Count)= xMidPhi
	    c0Y(c0Count)= yMidPhi
	  end if

	  if (levelsetS(levelR,j)<0 .and. i<hiXS) then			! neg. right
	    c0Count=c0Count+1
		c0X(c0Count)= xS(i) - levelsetS(i,j)*(xS(levelR)-xS(i))/ & 
		  (levelsetS(levelR,j)-levelsetS(i,j))
	    c0Y(c0Count)= yS(j)
	  end if

	  if (levelsetS(levelR,levelD)<0 .and. i>1 .and. j>1) then ! neg. lower right
	    call TwoD_Phi( xS(levelR),yS(levelD),levelsetS(levelR,levelD), xS(i), &
		  yS(j), levelsetS(i,j), xMidPhi, yMidPhi )

	    c0Count=c0Count+1
	    c0X(c0Count)= xMidPhi
	    c0Y(c0Count)= yMidPhi
	  end if

	  if (levelsetS(i,levelD)<0 .and. j>1) then				! neg. down
	    c0Count=c0Count+1
	    c0X(c0Count)= xS(i)
		c0Y(c0Count)= yS(j) - levelsetS(i,j)*(yS(j)-yS(levelD))/ &
		  (levelsetS(i,j)-levelsetS(i,levelD))
	  end if

	  if (levelsetS(levelL,levelD)<0 .and. i>1 .and. j>1) then ! neg. lower left
	    call TwoD_Phi( xS(i),yS(j),levelsetS(i,j), xS(levelL), &
		  yS(levelD), levelsetS(levelL,levelD), xMidPhi, yMidPhi )

	    c0Count=c0Count+1
	    c0X(c0Count)= xMidPhi
	    c0Y(c0Count)= yMidPhi
	  end if

	  if (levelsetS(levelL,j)<0 .and. i>1) then				! neg. left
	    c0Count=c0Count+1
		c0X(c0Count)= xS(i) - levelsetS(i,j)*(xS(i)-xS(levelL))/ &
		  (levelsetS(i,j)-levelsetS(levelL,j))
	    c0Y(c0Count)= yS(j)
	  end if

	  if (levelsetS(levelL,levelU)<0 .and. i>1 .and. j<hiYS) then ! neg. upper left
	    call TwoD_Phi( xS(i),yS(j),levelsetS(i,j), xS(levelL), &
		  yS(levelU), levelsetS(levelL,levelU), xMidPhi, yMidPhi )

	    c0Count=c0Count+1
	    c0X(c0Count)= xMidPhi
		c0Y(c0Count)= yMidPhi	
	  end if
    end if UPLevel
  end do Colc0
end do Rowc0

allocate(xarray(c0Count)); allocate(yarray(c0Count))

do i=1,c0Count
  xarray(i)=c0X(i)
  yarray(i)=c0Y(i)
end do

goto 999		! actual program end, just lists formats and errors below
! ***********

! FORMAT listings
10 FORMAT(3f12.3)
11 FORMAT(2f10.3, e12.3)
12 FORMAT(10f10.3)
13 FORMAT(f12.3)

999 print *, 'Done with subroutine <EField>'
end SUBROUTINE Intersection_Points

! ***************
! SUBROUTINE TwoD_Phi
! ***************
! DESCRIPTION: takes two grid points and opposite sign value of phi at each
!   and calculates the (x,y) point where phi=0
!
! INPUT:  (x,y,phi)1, (x,y,phi)2
!
! OUTPUT: (x*,y*) where phi=0
! 
! ALGORITHMN: 
!   1) line between (x,y)1 and (x,y)2
!   2) slope between line and delta(phi)
!   3) move from phi1 over line with decreasing slope until phi*=0. 
!        distance d* along line
!   4) solve new distance equation =d* (a) and line 1)
!   5) gives x, use 1) to get y
!
! Primary author: Anton VanderWyst, University of Michigan 

! Version history:
!   09 Jun 2004 begun
!	10 Jun working for all 4 diagonal directions

SUBROUTINE TwoD_Phi (x1, y1, phi1, x2, y2, phi2, xstar, ystar)

! incoming variables
real*8, intent(in)		:: x1, y1, phi1, x2, y2, phi2

! outsourced variables
real*8, intent(out)		:: xstar, ystar

! local subroutine variables
real					:: mP, bP, dP, dstar, A, B, C, tempMid

continue ! body of program
!*******
mP = (y2-y1) / (x2-x1)
bP = (x2*y1-x1*y2) / (x2-x1)				!(1)
dP = ( (x2-x1)**2 + (y2-y1)**2 )**0.5	
dstar = dP*phi1 / (phi1-phi2)				!(2)(3)

A = 1+mP**2									!(4)
B = -2*x1 + 2*mP*bP - 2*y1*mP
C = bP**2 - 2*y1*bP + y1**2 + x1**2 - dstar**2

tempMid = B**2 - 4*A*C 

!if (tempMid>0) then
!  xstar = (-B - ( tempMid )**0.5 )/(2*A) !(5)
!  ystar = mP*xstar + bP
!else
xstar = -phi1/(phi2-phi1) * (x2-x1) + x1
ystar = -phi1/(phi2-phi1) * (y2-y1) + y1
!  print *, 'tempMid <0. xy1, xy2, xy*, phi1, phi2:'
!  write(*,fmt='(a, 6f7.2, 2e10.2)'), '  ', x1, y1, x2, y2, &
!    xstar, ystar, phi1, phi2
!end if

end SUBROUTINE TwoD_Phi

! ***************
! SUBROUTINE VarSpreadOldNodes
! ***************
! DESCRIPTION: takes any point (x*,y*) with variable value and spreads it out 
!   to grid points (x1,y1) ... (x4,y4)
!
! INPUT:  (x*,y*,var), (x,y)1 -> (x,y)4, phi1->phi4, IgnoreVal, xVar, yVar, size
!
! OUTPUT: var1 -> var4
!
! ALGORITHMN: 
!   1) calc. distance from * pt to each grid point
!   2) divy up charge with more going to closer grid points
!
! Primary author: Anton VanderWyst, University of Michigan 
!
! Version history:
!   12 Jun 2004 begun

SUBROUTINE VarSpreadOldNodes (xstar, ystar, var, Axx, Ayy, &
	Bxx, Byy, Cxx, Cyy, Dxx, Dyy, Avar, Bvar, Cvar, Dvar, IgnoreVal, &
	xVar, yVar, sizexVar)

! incoming variables
integer, intent(in)		:: sizexvar, Axx, Ayy, Bxx, Byy, Cxx, Cyy, &
							Dxx, Dyy, IgnoreVal
real*8, intent(in)		:: xstar, ystar, var
real*8, dimension(sizexvar), intent(in) :: xVar, yVar

! outsourced variables
real*8, intent(out)		:: Avar, Bvar, Cvar, Dvar
! local subroutine variables
real					:: Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, &
							dA, dB, dC, dD, sumD
continue ! body of program
!*******

Ax=xVar(Axx); Ay=yVar(Ayy)
Bx=xVar(Bxx); By=yVar(Byy)
    
dA=min( 1./( (Ax-xstar)**2 + (Ay-ystar)**2 )**0.5, 1.e20 )
dB=min( 1./( (Bx-xstar)**2 + (By-ystar)**2 )**0.5, 1.e20 )

if (IgnoreVal==Cxx) then 
  dC= 0
else
  Cx=xVar(Cxx); Cy=yVar(Cyy)
  dC=min( 1./( (Cx-xstar)**2 + (Cy-ystar)**2 )**0.5, 1.e20 )
end if

if (IgnoreVal==Dxx) then 
  dD= 0
else
  Dx=xVar(Dxx); Dy=yVar(Dyy)
  dD=min( 1./( (Dx-xstar)**2 + (Dy-ystar)**2 )**0.5, 1.e20 )
end if

sumD = dA+dB+dC+dD

Avar = var*dA/sumD		! dA/( dA+dB+dC+dD)
Bvar = var*dB/sumD
Cvar = var*dC/sumD
Dvar = var*dD/sumD

end SUBROUTINE VarSpreadOldNodes

! ***************
! SUBROUTINE CalcE_Field
! ***************
! DESCRIPTION: given a location and panel strengths, calcs the EField there
! INPUT:  
!
! OUTPUT: 
!
! CALLING PROGRAMS: MAT_Greens
!
! Primary author: Anton VanderWyst, University of Michigan 
!
! Version history:
!   24 Sept 2004 begun
!	11 Oct	v5	changed to direct computation of field without interpolation

SUBROUTINE CalcE_Field( i, k, M, xOff, yOff, xPanEF, yPanEF, &
  normX, normY, PanWeight, PanLength, PanFT, sigEF, EfieldEF )

! incoming variables
integer, intent(in)		:: i, k, M, PanFT
real*8, intent(in)		:: xPanEF, yPanEF, normX, normY, &
							PanWeight, sigEF, PanLength, xOff, yOff
! outgoing variables
real*8, dimension(2), intent(inout)	:: EfieldEF

! routine variables
real*8					:: delX, delY, GreenPart, fullNorm, noNorm

continue

fullNorm=1.; noNorm=0.
!do k = 1,M
!  panel=(i-1)*M + k

delX=xOff-xPanEF
delY=yOff-yPanEF

! calculate Ex since normPanX(j) = 1

call TrapInit ( GreenPart, delX, delY, normX, normY, &
  fullNorm, noNorm, PanWeight, PanLength, PanFT, 0, k, M )

EfieldEF(1)= EfieldEF(1) + GreenPart*sigEF

! calculate Ey since normPanY(j) = 1
call TrapInit ( GreenPart, delX, delY, normX, normY, &
  noNorm, fullNorm, PanWeight, PanLength, PanFT, 0, k, M )

EfieldEF(2)= EfieldEF(2) + GreenPart*sigEF

end SUBROUTINE CalcE_Field

! ***************
! SUBROUTINE ForceZero
! ***************
! DESCRIPTION: forces an variable to zero if they are too close
! INPUT: variable, how close
!
! OUTPUT: variable 
!
! CALLING PROGRAMS: MAT_Greens
!
! Primary author: Anton VanderWyst, University of Michigan 
!
! Version history:
!   07 Oct 2004 begun

SUBROUTINE ForceZero (ForcedVariable, HowClose)

! incoming variables
real*8, intent(inout)		:: ForcedVariable
real*8, intent(in)			:: HowClose

continue

if ( abs(ForcedVariable)-HowClose < 0 ) then
	ForcedVariable =0
! else no change
endif

end SUBROUTINE ForceZero


! ***************
! SUBROUTINE ExtrapEdgePoints
! ***************
! DESCRIPTION: smoothes out edges between Neumann and Dirichlet boundary panels
! 
! INPUT: panel strengths
!
! OUTPUT: extrapolated panel strengths, steady value panel number
!
! CALLING PROGRAMS: MAT_Greens
!
! Primary author: Anton VanderWyst, University of Michigan 
!
! Version history:
!   06 Nov 2004 begun

SUBROUTINE ExtrapEdgePoints (PanelNum, PrevSig, Sig, NextSig, StartSteady, ExtrapSigma, CloseSig)

! incoming variables
integer, intent(in)			:: PanelNum
real*8, intent(in)			:: NextSig, PrevSig, Sig, CloseSig

! outgoing variables
integer, intent(out)		:: StartSteady
real*8, intent(out)			:: ExtrapSigma

continue

if (abs(PrevSig-Sig) < CloseSig*Sig .AND. abs(Sig-NextSig) < CloseSig*Sig) then	
  ! You're close enough to call it steady state
  StartSteady = PanelNum

endif

ExtrapSigma = Sig
end SUBROUTINE ExtrapEdgePoints

end module Wrapper6