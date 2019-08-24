MODULE SortNodePts3
implicit none

! used to restrict variables to prevent accidental calling outside
private      ! everything defaults to 'private' except those expressly labeled otherwise.
public			:: SortPts, CountCoords, TrapInit, GaussQuad_data, &
	InSideYesNo, ExactGreenSolve, pi	! subroutines
real, parameter :: pi=3.14159265358979323846264338327950288419716939937510582097494459, epsil = 1e-6

contains 

! ***************
! SUBROUTINE SortPts[x]
! ***************
! DESCRIPTION: takes unsorted list of c0 surface points and produces
!   a list where one point links to the following one
!
! INPUT: LevelPtsX, LevelPtsY 
!
! OUTPUT: SortedPts
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
!   18 Aug 2004 begun
!   28 Aug      added neumann panel locations in for ring electrode
!	31 Aug		MatLab code frozen, moving to Fortran
!	19 Sept v2  changing to fit Sussman code. Using passed in grid value

SUBROUTINE SortPts( Coords, sizeLevel, xarraySP, yarraySP, minX, &
	maxX, minY, maxY )

! incoming variables
integer, intent(in)				    :: sizeLevel
real*8, dimension(sizeLevel), intent(in) :: xarraySP, yarraySP
real*8, intent(in)					:: minX, minY, maxX, maxY

! outgoing variables
real*8, dimension(sizeLevel+7,4), intent(out) :: Coords

! subroutine entirely variables
integer, dimension(sizeLevel+7)		:: indexIn
integer								:: nearPtNow, nearPtNext, i, j, xPosition

real*8								:: outrangeVar, dist, distNeighbor, &
									   tempXarray

continue
!*************
! main program commands begin

outrangeVar=1e6
do i=2,sizeLevel
  indexIn(i)=1	! forces the variable in the loop to be all points
end do

! description of coordinate columns
! [1-2] coordinates to form a trapezoid around to test slanted panels
! [3] potential of boundary panel (if 1==[4]) or slope of panel (if 0==[4])
! [4] 1=dirch; 0=neumann boundary conditions

tempXarray=1/epsil
do i=1,sizeLevel
  if (xarraySP(i) < tempXarray) then
    tempXarray=xarraySP(i)
	xPosition=i
  endif
enddo
indexIn(xPosition)=0; nearPtNow=xPosition
Coords(1,1)=xarraySP(xPosition)
Coords(1,2)=yarraySP(xPosition)

! loops through and organizes 
outer: do i=1,sizeLevel-1
    dist=outrangeVar; nearPtNext=0
	do j=1,sizeLevel
        if ( 1==indexIn(j) ) then
			distNeighbor = (xarraySP(nearPtNow)-xarraySP(j))**2 + &
				(yarraySP(nearPtNow)-yarraySP(j))**2

            if ( distNeighbor < dist ) then
                dist = distNeighbor
                nearPtNext=j
            endif
        endif
    end do
    
    Coords(i+1,1)=xarraySP(nearPtNext)
    Coords(i+1,2)=yarraySP(nearPtNext)
    indexIn(nearPtNext)=0; nearPtNow=nearPtNext
end do outer 

  ! adding wrap-around points
  Coords(i+1,1)= maxX;		Coords(i+1,2)= Coords(i,2)
  Coords(i+2,1)= maxX;		Coords(i+2,2)= maxY
  Coords(i+3,1)= 0.8*maxX;	Coords(i+3,2)= maxY
  Coords(i+4,1)= 0.7*maxX;	Coords(i+4,2)= maxY
  Coords(i+5,1)= 0.3*maxX;	Coords(i+5,2)= maxY
  Coords(i+6,1)= 0.2*maxX;	Coords(i+6,2)= maxY
  Coords(i+7,1)= minX;		Coords(i+7,2)= maxY
  !Coords(i+8,1)= minX;		Coords(i+8,2)= Coords(i,2)
  ! define the fluid panels as potential 2.0 and dirichlet(==1)
  Coords(1:sizeLevel,3) = -2.0
  Coords(1:sizeLevel,4) = 1

  ! now define the added points as N/D and their fluxes...
  Coords(i+1,3)= 0;       Coords(i+1,4)= 0
  Coords(i+2,3)= 0;       Coords(i+2,4)= 0
  Coords(i+3,3)= -1000;   Coords(i+3,4)= 1
  Coords(i+4,3)= 0;       Coords(i+4,4)= 0
  Coords(i+5,3)= -1000;   Coords(i+5,4)= 1
  Coords(i+6,3)= 0;       Coords(i+6,4)= 0
  Coords(i+7,3)= 0;       Coords(i+7,4)= 0
  
goto 999		! actual program end, just lists formats and errors below
! ***********

! FORMAT listings
10 FORMAT(e20.8)

! ERROR listings

999 print *, 'Done with subroutine <SortPts>'

end SUBROUTINE SortPts

! ***************
! SUBROUTINE CountCoords
! ***************
! DESCRIPTION: counts the number of coordinates in the interpolated free surface
!
! INPUT: < >
!
! OUTPUT: NumCoords
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
!   31 Aug 2004 begun

SUBROUTINE CountCoords( NumCoords )

integer, intent(out)	:: NumCoords
integer, parameter		:: XDataPts=10
real					:: a

continue
NumCoords=0
! find size of levelset surface coordinate list
open(unit=XDataPts, file='LevelPtsX.dat', status='old', action='read', err=986)
do	
  NumCoords=NumCoords+1
  read(unit=XDataPts, fmt=11, err=987, end=100) a
end do
100 NumCoords=NumCoords-1
close (unit=XDataPts, err=989, status='keep')

goto 999		! actual program end, just lists formats and errors below
! ***********

! FORMAT listings
11 FORMAT(f10.3)

! ERROR listings
986 print *, 'Error open Count Coords'
987 print *, 'Error read Count Coords'
988 print *, 'Error write Count Coords'
989 print *, 'Error close Count Coords'

999 print *, 'Done with subroutine <CountCoords>'
end SUBROUTINE CountCoords

! ***************
! SUBROUTINE GaussQuad_data
! ***************
! DESCRIPTION: acts as a lookup to provide the panel points M locations and weights
!
! INPUT: # points / panel, panel beginning and end points
!
! OUTPUT: location of Gaussian quadrature points along panel
!
! CALLING PROGRAM: MAT_Greens
!
! Primary author: Anton VanderWyst, University of Michigan 
!
! Version history:
!   01 Jul 2004 begun  
!	07 Oct   v3	removing the trapezoid component. Only Gaussian used

!SUBROUTINE [PtlocationX, PtlocationY, Ptweight] = GaussQuad_data (PtsPerPanel, ...
!    PanBeginX, PanLenX, PanstepSizeX, PanBeginY, PanLenY, PanstepSizeY, j, IntMethod)

SUBROUTINE GaussQuad_data ( PtlocationX, PtlocationY, Ptweigh, PtsPerPanel, &
    PanBeginX, PanLenX, PanBeginY, PanLenY, j )

! incoming variable
integer, intent(in)		:: PtsPerPanel, j
real*8, intent(in)		:: PanBeginX, PanLenX, PanBeginY, PanLenY

! outgoing variables
real*8, intent(out)		:: PtlocationX, PtlocationY, Ptweigh

! subroutine entirely variables
real					:: PanEndX, PanEndY
real, dimension(PtsPerPanel) :: Mlocation, Mweight
character				:: IntMethod

continue
IntMethod = 'g'

Method: select case( IntMethod )
    case ('g')
        PanEndX=PanBeginX+PanLenX
        PanEndY=PanBeginY+PanLenY
		PtsPan: select case( PtsPerPanel )
            case (2)
                Mlocation(2)= 0.5773502692; Mweight(2)=1.
                Mlocation(1)= -0.5773502692; Mweight(1)=1.
            case (3)
                Mlocation(3)= 0.7745966692; Mweight(3)= 5/9
                Mlocation(2)= 0;            Mweight(2)= 8/9  
                Mlocation(1)= -0.7745966692; Mweight(1)= 5/9
            case (4)
                Mlocation(4)= 0.8611363116; Mweight(4)= 0.3478548451
                Mlocation(3)= 0.3399810436; Mweight(3)= 0.6521451549  
                Mlocation(2)= -0.3399810436; Mweight(2)= 0.6521451549
                Mlocation(1)= -0.8611363116; Mweight(1)= 0.3478548451
            case (5)
                Mlocation(5)= 0.90617985; Mweight(5)= 0.23692689
                Mlocation(4)= 0.53846931; Mweight(4)= 0.47862867
                Mlocation(3)= 0.0;         Mweight(3)= 0.568888889  
                Mlocation(2)= -0.53846931; Mweight(2)= 0.47862867
                Mlocation(1)= -0.90617985; Mweight(1)= 0.23692689
            case (6)
                Mlocation(6)= 0.93246951; Mweight(6)= 0.17132449
                Mlocation(5)= 0.66120939; Mweight(5)= 0.36076157  
                Mlocation(4)= 0.23861918; Mweight(4)= 0.46791393
                Mlocation(3)= -0.23861918; Mweight(3)= 0.46791393  
                Mlocation(2)= -0.66120939; Mweight(2)= 0.36076157
                Mlocation(1)= -0.93246951; Mweight(1)= 0.17132449        
             case (8)  
                Mlocation(8)= 0.96028986; Mweight(8)= 0.10122854
                Mlocation(7)= 0.79666648; Mweight(7)= 0.22238103  
                Mlocation(6)= 0.52553241; Mweight(6)= 0.31370665
                Mlocation(5)= 0.18343464; Mweight(5)= 0.36268378  
                Mlocation(4)= -0.18343464; Mweight(4)= 0.36268378
                Mlocation(3)= -0.52553241; Mweight(3)= 0.31370665  
                Mlocation(2)= -0.79666648; Mweight(2)= 0.22238103
                Mlocation(1)= -0.96028986; Mweight(1)= 0.10122854    
            case (10)
                Mlocation(10)= 0.97390653; Mweight(10)= 0.06667134
                Mlocation(9)= 0.86506337; Mweight(9)= 0.14945135  
                Mlocation(8)= 0.67940957; Mweight(8)= 0.21908636
                Mlocation(7)= 0.43339539; Mweight(7)= 0.26926672  
                Mlocation(6)= 0.14887434; Mweight(6)= 0.29552422
                Mlocation(5)= -0.14887434; Mweight(5)= 0.29552422  
                Mlocation(4)= -0.43339539; Mweight(4)= 0.26926672
                Mlocation(3)= -0.67940957; Mweight(3)= 0.21908636  
                Mlocation(2)= -0.86506337; Mweight(2)= 0.14945135
                Mlocation(1)= -0.97390653; Mweight(1)= 0.06667134            
            case default
				stop 'Data not entered for this number of points per panel.'
        end select PtsPan

        PtlocationX=(PanEndX+PanBeginX)/2 + Mlocation(j)*(PanEndX-PanBeginX)/2
        PtlocationY=(PanEndY+PanBeginY)/2 + Mlocation(j)*(PanEndY-PanBeginY)/2
        Ptweigh=Mweight(j)
           
	case default
		stop 'GaussQuad integration method not chosen'
end select Method

end SUBROUTINE GaussQuad_data

! ***************
! SUBROUTINE TrapInit
! ***************
! DESCRIPTION: integrates green's function for a 2D cartesian field along panels 
!   using the trapezoid rule and Gaussian quadrature
!
! INPUT: check coord, comparison coord, loop values, potential points
!
! OUTPUT: ( grad(G).n ) ds, ( grad(phi).n ) ds = aN(ij)
!
! CALLING PROGRAM: MAT_Greens
!
! Primary author: Anton VanderWyst, 
!	University of Michigan 
!
! Version history:
!   21 Jun 2004 begun
!           v2. removed neumann components
!   29 Jun  v3. checking trapezoid approx vs. v2's exact sum  
!   01 Jul  v4. (v3) is correct, adding Gaussian quadrature   
!   06 Jul  v5. (v4) is correct, adding mixed BC
!   08 Jul  v6. v5 is _not_ working, jumping to arbitrary panel
!                   orientation case where z,a means dissolve
!   13 Aug  v7. (v6) is correct, adding flag to do trapezoid rule 
!                   for 2nd part of potential calc, after c0 strengths
!                   determined. Changes if pt gets too close to panel

SUBROUTINE TrapInit ( GreenPart, delX, delY, normXj, normYj, normXi, normYi, &
    PanWeight, PanLen, PanFromTo, flagTrap, k, Mpts )

! incoming variables
integer, intent(in)		:: PanFromTo, flagTrap, k, Mpts
real*8, intent(in)		:: PanLen, delX, delY, normXj, normYj, normXi, normYi, PanWeight

! outgoing variables
real*8, intent(out)		:: GreenPart

! subroutine entirely variables
real					:: HalfLen							

continue
GreenPart = 0.

if ( 0==flagTrap .OR. 1==k .OR. Mpts==k ) then
	HalfLen = 0.5
else
	HalfLen = 1.0 
end if

if ( abs(delX)<epsil .AND. abs(delY)<epsil ) then          
    ! do nothing, on the border                
else
    select case (PanFromTo)    ! 1=dirichlet, 0=neumann type BC
    case (11)
        GreenPart = PanLen * HalfLen * PanWeight * GradGreens(delX, delY, normXj, normYj)
    case (10) 
        GreenPart = PanLen * HalfLen * PanWeight * Greens(delX, delY)
	case (00)
        GreenPart = PanLen * HalfLen * PanWeight * GradGreens(delX, delY, normXi, normYi)
    case (01)
        GreenPart = PanLen * HalfLen * PanWeight * GradGradGreens(delX, delY, &
            normXj, normYj, normXi, normYi)
    case default
		stop 'Error on PanFromTo'
    end select ! switch PanFromTo
end if    ! on border

end SUBROUTINE TrapInit

	!-------------------------
	FUNCTION Greens( zzz, aaa )
	! Greens subfunction
	
	real*8, intent(in)	:: zzz, aaa
	real*8				:: Greens
		
	continue

	Greens = -log( (zzz)**2. + (aaa)**2. ) / (4.*pi)
	end FUNCTION Greens

	!-------------------------
	FUNCTION GradGreens( delX, delY, normX, normY )
	! grad(G) subfunction

	real*8				:: GradGreens
	real*8, intent(in)	:: delX, delY, normX, normY
 
	continue
	!   d/dn {ln (r_ij) } = grad_i {ln (r_ij) } . n_i
	! = (nx*dG/dx+ny*dG/dy)*dS

	GradGreens = ( delX*normX + delY*normY ) / ( ( delX**2 + delY**2 )*(2.*pi) )
	end FUNCTION GradGreens

	!-------------------------
	FUNCTION GradGradGreens(delX, delY, normXj, normYj, normXi, normYi)
	!neumann-dirichlet grad.grad subfunction
	! d/dn_i (d/dn_j {ln (r_ij) }) = grad_i (grad_j {ln (r_ij) } . n_j) . n_i

	real*8, intent(in)	:: delX, delY, normXj, normYj, normXi, normYi
	real*8				:: GradGradGreens, R, gamma

	continue
	R = ( delX**2 + delY**2 )    ! = r^2

	if ( 0 /= delX .AND. 0 /= delY) then 
		gamma = -R/(2*delX*delY)

		GradGradGreens = (-1/(2*gamma*pi*R))*( ( (gamma+delX/delY)*normXj + normYj)*normXi + &
        ( (gamma+delY/delX)*normYj + normXj)*normYi )
	elseif ( 0==delX ) then
		GradGradGreens = (-1/(2*pi))*( ( (R**-1)*normXj )*normXi + &
        ( (R**-1 - 2*delY**2/R**2)*normYj)*normYi )
    elseif ( 0==delY ) then
		GradGradGreens = (-1/(2*pi))*( ( (R**-1 - 2*delX**2/R**2)*normXj )*normXi + &
			( (R**-1)*normYj )*normYi )
	else
		stop 'Error on gradgradG'
	end if

	end FUNCTION GradGradGreens

! ***************
! FUNCTION ExactGreenSolve
! ***************
! DESCRIPTION: integrates green's function for a 2D cartesian field exactly. 
!   Similar to 'TrapInit', but uses calculus instead of numerics
!
! INPUT: check coord, panel center
!
! OUTPUT: int(G)
!
! CALLING PROGRAM: MAT_Greens
!
! Primary author: Anton VanderWyst, 
!	University of Michigan 
!
! Version history:
!   18 Oct 2004 begun

FUNCTION ExactGreenSolve ( x, a )
! incoming variables
real*8, intent(in)		:: x, a
real*8					:: ExactGreenSolve

continue
! make sure that the 'X' term is the one that you want to integrate with respect to.
! int(ln(x^2+a^2) dx = x ln(x^2+a^2) -2x +2a tan^-1(x/a)

ExactGreenSolve = (x*log(x**2 + a**2) - 2.*x + 2.*a*atan2(a,x))/(4.*pi)

end FUNCTION ExactGreenSolve

! ***************
! SUBROUTINE InSideYesNo
! ***************
! DESCRIPTION: Determines the horizontal and vertical intersections for pt
!   P onto the line between pts A,B
!
! CALLING PROGRAM: MAT_Greens[x]
!
! INPUT: X,Y sweep direction, examine pt, segment end coordinates, # 
!   previous crossings
!
! OUTPUT: # and location of intersections
!
! Primary author: Anton VanderWyst, University of Michigan 
!
! Version history:
!   11 Jul 2004 begun
!   18 Jul trying to use # of intersections as in/out determiner

!SUBROUTINE [CrossCount, CrossYesNo, IntersecX, IntersecY] = InSideYesNo(CrossDir, ...
!    Ax, Ay, Bx, By, Px, Py, CrossCount)

SUBROUTINE InSideYesNo( CrossCount, CrossYesNo, IntersecX, IntersecY, CrossDir, &
    Ax, Ay, Bx, By, Px, Py )

! incoming variables
real, intent(in)		:: Ax, Ay, Bx, By, Px, Py
character, intent(in)	:: CrossDir

! outgoing variables
logical, intent(out)	:: CrossYesNo
integer, intent(inout)	:: CrossCount
real, intent(out)		:: IntersecX, IntersecY

! subroutine entirely variables
real					:: delBAx, delBAy, t

continue

CrossYesNo=.FALSE.
delBAx = Bx-Ax; delBAy = By-Ay

!write(*, fmt='(a, f7.3, a, f7.3)') 'delBAx= ', delBAx, ', delBAy= ', delBAy  

select case( CrossDir )
case ('Y')
    if ( 0 /= delBAy ) then 
        t = (Py - Ay)/delBAy
        IntersecX = Ax + (Py-Ay)*delBAx/delBAy
        IntersecY = t*delBAy + Ay
    else 
        t = (Px - Ax)/(delBAx)                
        IntersecX = Px; IntersecY = Ay
    endif
case ('X')
    if ( 0 /= delBAx ) then                
        t = (Px - Ax)/delBAx  
        IntersecX = t*delBAx + Ax
        IntersecY = Ay + (Px-Ax)*delBAy/delBAx
    else 
        t = (Py - Ay)/(delBAy)   
        IntersecX = Ax; IntersecY = Py
    endif
case default
	stop 'Incorrect CrossDir'
end select

if ( t>=0-epsil .AND. t<=1+epsil ) then
    CrossCount = CrossCount+1; CrossYesNo=.TRUE.
endif

end SUBROUTINE InSideYesNo

end MODULE SortNodePts3