MODULE MAT_Greens16

  use SortNodePts3,  only : SortPts, TrapInit, GaussQuad_data, &
       InSideYesNo, ExactGreenSolve, pi
  use Inverters,	 only : LEGS
  use GMRes_Andrew5, only : gmresIter
  use CubicSpline5,  only : CurveFit, NewCoords ! subroutine, variable

implicit none

! used to restrict variables to prevent accidental calling outside
private      ! everything defaults to 'private' except those expressly labeled otherwise.
public		:: MAT_Greens, epsil	! subroutines
real*8, parameter :: epsil = 1e-6, toolarge=1e8, epsil_0 = 8.85e-12

contains 

! ***************
! Program MAT_Greens[x]
! ***************
! DESCRIPTION: test program that calculates the boundary integral
!   weights necessary to give accurate potentials
!
! INPUT: Coords, SelectedGrid
!
! OUTPUT: sigmas, potentials at SelectedGrid
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   15 Jun 2004 begun
!   18 Jun  v2. Swapping the RHS calculation to a 4x1 sum. Calcing sides 
!                   separately
!   21 Jun  v5. introduced trapezoid integration and ds normalizing
!           v6. tried direct integration for boundaries
!   30 Jun  v7. scrapped, trying Jerry's method to get sum(weights)=0.
!   01 Jul  v8. (v7) working, implementing optional Gaussian quadrature
!   03 Jul  v9. (v8) working, implementing mixed Neumann/Dirichlet BC
!   07 Jul  v10. skipping over v9, and changing so there is a relatively 
!                   constant length of panel, _not_ constant #/side
!   08 Jul  v11. (v10) working, implementing arbitrary angled panels
!                   with just dirichlet boundary conditions. Removing trapezoid.
!                   Allows more than 4 sides total. Normals facing inward
!   14 Jul  v12. difficulties with the inside/outside part for an arbitrary shape.
!                   Trying to implement case of drawing subfigures in/out
!   13 Aug  v14. skipped 13, computing close points to boundaries (used when 
!                   diagonal panels are present) with increasing panel colocation
!                   density. Switching to trapezoid rule for arb. potential calcs.
!   18 Aug  v15. (v14) almost totally working, but 'PanLen' is too arbitrary whether
!                   it is +/-. Trying to recast so aN will be correct with an
!                   always + 'PanLen'. xPanLen, yPanLen unchanged. 
!   19 Aug  v16. (v15) working, now that the diagonal is always +0.5(!). PanLen
!                   is constant +, so aN is also all +. Current problem is huge 
!                   number of panels (14k+) on test problem, so attempting strategies
!                   like NPan = 1 for dx < minimum delta_x 
!   21 Aug  v17. (v16) working, testing alternative strategy of removing close points. 
!   25 Aug  v18. (v17) tentatively working. Adding neumann BC panels
!   29 Aug  v19. (v18) almost working. Debugging the automatic points remover, NewCoords
!					this is being ported over from Matlab v6.0
!	07 Sep	v03. (v02) working! Uses Gaussian elimination with partial pivoting to solve
!					aN*sigma = RHS. However, as that is O(N^3), changing to iterative
!					solver GMRES here. Uses NAG library calls for F11BDF and F11BEF
!	19 Sep  v05. (v04) skipped. Moving from v02 to now compute potentials only along
!					supplied array. Changed to subroutine to use with CLSVOF from FL.
!					Completely re-doing potential calculation to determine Ex, Ey fields
!	22 Sep	v06. (v05) partially working. Changing E-field calc to offset distance, back
!					to Gaussian colocation points. Removing MEffect
!	04 Oct	v07. (v06) partially working. Implementing 4 term source so can use direct
!					A,D solving _without_ interpolation
!	07 Oct	v08. (v07) not working. Changing panel centers to a Chebyshev distribution
!					between fixed NewCoords
!	19 Oct	v09. (v08) not working. Continuing shift to 4 term sources. C. dist correct
!	26 Oct  v10. (v09) WORKING!!!! Stripping out extra calculations and computations not 
!						relevant to returning only the electric field. Calculated where 
!						the end of 'sizeLevel' is after removing points
!	29 Oct  v11. (v10) working. Adding in extrapolating for the end points to smooth out 
!						the ripples near N/D panel boundaries. Also, calc'ing F=qE instead
!						of just E_n
!   12 Nov	v12. (v11) working. Now forming cublic splines to smooth out original points.
!	21 Nov  v13. (v12) working. At large numbers of panels, the solution oscillates. Also, trying
!						to implement pointers for faster array iterative inversions.
!	07 Dec  v14. (v13) working. Moving some 2D arrays to NxN 1D to speed up run times.
!						Now have working cubic splines, GMRes routines
!	08 Dec	v15. (v14) working. However, it isn't any faster, so reverting to near-v13. Removing 
!						print statements 
!	21 Dec  v16. (w15) working. Changing variable definition to 'parameter' or defined incode to
!						allow repetitive calling. Removing more write statements

SUBROUTINE MAT_Greens( level, llo1, llo2, lhi1, lhi2, gridlo1, gridlo2, &
  gridhi1, gridhi2, dx, dy, domlo1, domlo2, domhi1, domhi2, xarray, &
  yarray, sizeLevel, Efield ) 

! incoming variables
  integer, intent(in)	:: gridlo1, gridlo2, gridhi1, gridhi2, &
    llo1, llo2, lhi1, lhi2, sizeLevel
  real*8, intent(in)	:: dx, dy, domlo1, domlo2, domhi1, domhi2
  real*8, dimension(llo1:lhi1, llo2:lhi2), intent(in) :: level
  real*8, dimension(sizeLevel), intent(in)	:: xarray, yarray

! outgoing variable
  real*8, dimension(sizeLevel), intent(out) :: Efield

! program variables
  ! fixed variables
  integer, parameter	:: M=10, N=1, & ! # pts / panel; # panels / shortest side
	CubicSplineFit=0  			! if '0', linear spline fit. if '1', cubic spline fit.

  ! regular variables
  logical		:: LoopTotPan

  integer		:: NumCoords, i, j, k, aLinePan, panel, PanFromTo, CloseCoord, NextPan, totPanels
  integer, dimension(:), allocatable		:: indxLegs
 	  
  real*8		:: t1, t2, timeCount, delX, delY, xCoordDist, yCoordDist, &
	GreenPart, GreenPart2, ai, bi, ci, di, SizeLevelDistPrev, SizeLevelDist, tempEpsil, tempCoeff
  real*8, dimension(:), allocatable			:: xPanLen, yPanLen, PanLen, xPan, yPan, &
	xPanCen, yPanCen, xPanBegin, yPanBegin, delXCoord, delYCoord, sigma, distNewCoord, &
	PanGaussWeight, normPanX, normPanY, RHS_gamma_phi, RHS_gamma_dphi, &
	RHS_alpha_phi, RHS_alpha_dphi, RHStemp, aLine
  real*8, dimension(:,:), allocatable		:: Coords, RHS

  continue
  !**********
  ! main program commands begin
  !**********
  call cpu_time(t1)

  allocate(Coords(sizeLevel+7,4))
  call SortPts ( Coords, sizeLevel, xarray, yarray, domlo1, domhi1, domlo2, domhi2 )
  NumCoords = sizeLevel+7

  ! basic error checking to use Gaussian quadrature as a method of numerically
  !   integrating the panels later
  if ( mod(M,2) /= 0 ) then
    stop '"M" needs to be even to avoid being the center of the panel'
  endif
  if ( M>10 )  then
    stop 'FYI: Gaussian quadrature not coded up to that level of precision'
  endif

  allocate(delXCoord(NumCoords)); allocate(delYCoord(NumCoords))

  do i=1,NumCoords
    if ( i /= NumCoords) then
      delXCoord(i)=Coords(i+1,1)-Coords(i,1) 
      delYCoord(i)=Coords(i+1,2)-Coords(i,2)
    else
      delXCoord(i)=Coords(1,1)-Coords(i,1)         
      delYCoord(i)=Coords(1,2)-Coords(i,2)
    endif
  enddo

  ! create spaced points along the spline. Returns totPanels, NewCoords(global)
  call CurveFit (totPanels, Coords(:,1), Coords(:,2), delXCoord, delYCoord, &
    NumCoords, Coords(:,3:4), sizeLevel, CubicSplineFit)

  ! generate panels locations, normals and colocation pt. weights
  allocate(xPan(totPanels*M)); allocate(yPan(totPanels*M))
  allocate(PanGaussWeight(totPanels*M)); allocate(distNewCoord(totPanels))
  allocate(sigma(totPanels))

  allocate(xPanBegin(totPanels)); allocate(yPanBegin(totPanels))
  allocate(xPanLen(totPanels)); allocate(yPanLen(totPanels))
  allocate(PanLen(totPanels))
  allocate(xPanCen(totPanels)); allocate(yPanCen(totPanels))

  allocate(normPanX(totPanels)); allocate(normPanY(totPanels))
  allocate(RHS(totPanels,2)); allocate(RHStemp(totPanels))
  allocate(aLine(totPanels*totPanels)); allocate(indxLegs(totPanels))
  allocate(RHS_gamma_phi(totPanels)); allocate(RHS_gamma_dphi(totPanels)) 
  allocate(RHS_alpha_phi(totPanels)); allocate(RHS_alpha_dphi(totPanels))

  ! define parts of panels
  PanDefineFrom: do i = 1,totPanels
	if (i /= totPanels) then
	  delX = NewCoords(i+1,1) - NewCoords(i,1)
	  delY = NewCoords(i+1,2) - NewCoords(i,2)
	else
	  delX = NewCoords(1,1) - NewCoords(i,1)
	  delY = NewCoords(1,2) - NewCoords(i,2)
	endif
    distNewCoord(i) = (delX**2 + delY**2)**0.5
    normPanX(i)=-delY/distNewCoord(i); normPanY(i)=delX/distNewCoord(i)

	xPanBegin(i)=NewCoords(i,1)
	yPanBegin(i)=NewCoords(i,2)

    ! interpolating potentials on the RHS if NPan > 1 between coordinates      
    RHS(i,2) = NewCoords(i,4)
    RHS(i,1)= NewCoords(i,3)   

    ! finding the centers of the panels   
    if ( i>1 ) then 
      xPanCen(i-1)=(xPanBegin(i)-xPanBegin(i-1))/2 + xPanBegin(i-1)
      yPanCen(i-1)=(yPanBegin(i)-yPanBegin(i-1))/2 + yPanBegin(i-1)
    endif
  end do PanDefineFrom  ! for i

  ! Spread the collocation points IN EACH PANEL according to Gaussian quadrature
  PanDefineFrom2: do i = 1,totPanels
	if (i<totPanels) then
	  xPanLen(i) = xPanBegin(i+1)-xPanBegin(i)
	  yPanLen(i) = yPanBegin(i+1)-yPanBegin(i)
	else
	  xPanLen(i) = xPanBegin(1)-xPanBegin(i)
	  yPanLen(i) = yPanBegin(1)-yPanBegin(i)
    endif
	
	PanLen(i)= (xPanLen(i)**2 + yPanLen(i)**2)**0.5
	
	if (0==PanLen(i)) then
	  stop 'PanLen error'
	endif

    PanDefineTo2: do j = 1,M
      panel= (i-1)*M+j

      call GaussQuad_data (xPan(panel), yPan(panel), PanGaussWeight(panel), M, &
		xPanBegin(i), xPanLen(i), yPanBegin(i), yPanLen(i), j )
	enddo PanDefineTo2
  enddo PanDefineFrom2

  deallocate(distNewCoord); deallocate(delXCoord); deallocate(delYCoord)

  xPanCen(totPanels)=(xPanBegin(1)-xPanBegin(totPanels))/2 + xPanBegin(totPanels)
  yPanCen(totPanels)=(yPanBegin(1)-yPanBegin(totPanels))/2 + yPanBegin(totPanels)  

  ! green's function from (i) to (j) through (k) pts per panel. LHS calc
  GreenJ: do j = 1,totPanels		
    RHS_gamma_phi(j)=0; RHS_gamma_dphi(j)=0
    RHS_alpha_phi(j)=0; RHS_alpha_dphi(j)=0
	GreenI: do i = 1,totPanels
	  PanFromTo = RHS(j,2)*10 + RHS(i,2)
	  aLinePan = (j-1)*totPanels+i
	  aLine(aLinePan)=0.
	  ai=0; bi=0; ci=0; di=0
      do k = 1,M
	    panel=(i-1)*M + k
		delX=xPanCen(j)-xPan(panel)
		delY=yPanCen(j)-yPan(panel)

		select case (PanFromTo)		! 1=dirichlet, 0=neumann type BC
		case (11)
		  call TrapInit ( GreenPart, delX, delY, normPanX(i), &
   		    normPanY(i), normPanX(j), normPanY(j), PanGaussWeight(panel), &
		    PanLen(i), 10, 0, k, M )

		  ai = ai + GreenPart

		  call TrapInit ( GreenPart, delX, delY, normPanX(i), &
		    normPanY(i), normPanX(j), normPanY(j), PanGaussWeight(panel), &
		    PanLen(i), 11, 0, k, M )

		  RHS_gamma_phi(j) = RHS_gamma_phi(j) + RHS(i,1)*GreenPart

		case (10)
		  call TrapInit ( GreenPart, delX, delY, normPanX(i), &
		    normPanY(i), normPanX(j), normPanY(j), PanGaussWeight(panel), &
			PanLen(i), 10, 0, k, M )

		  RHS_gamma_dphi(j) = RHS_gamma_dphi(j) + RHS(i,1)*GreenPart

		  call TrapInit ( GreenPart, delX, delY, normPanX(i), &
		    normPanY(i), normPanX(j), normPanY(j), PanGaussWeight(panel), &
			PanLen(i), 11, 0, k, M )

		  bi = bi - GreenPart 		  

		case (01)
		  call TrapInit ( GreenPart, delX, delY, normPanX(i), &
		    normPanY(i), normPanX(j), normPanY(j), PanGaussWeight(panel), &
			PanLen(i), 11, 0, k, M )

		  RHS_alpha_phi(j) = RHS_alpha_phi(j) + RHS(i,1)*GreenPart

		  call TrapInit ( GreenPart, delX, delY, normPanX(i), &
		    normPanY(i), normPanX(j), normPanY(j), PanGaussWeight(panel), &
		    PanLen(i), 10, 0, k, M )

		  ci = ci + GreenPart 

     	case (00)
		  call TrapInit ( GreenPart, delX, delY, normPanX(i), &
		    normPanY(i), normPanX(j), normPanY(j), PanGaussWeight(panel), &
			PanLen(i), 10, 0, k, M )

	      RHS_alpha_dphi(j) = RHS_alpha_dphi(j) + RHS(i,1)*GreenPart

		  if (i /=j) then
		  	call TrapInit ( GreenPart, delX, delY, normPanX(i), &
		      normPanY(i), normPanX(j), normPanY(j), PanGaussWeight(panel), &
		      PanLen(i), 00, 0, k, M )

			di = di + GreenPart
		  elseif ( M==k ) then
			di = 0.5
		  endif 
		case default
		  stop 'Error on PanFromTo aLine'
		end select ! switch PanFromTo
	  end do		! for k
	  aLine(aLinePan) = ai + bi + ci + di
    end do GreenI   ! for i   
  end do GreenJ		! for j
  !RHS(:,1) = RHS(:,1) + RHS_gamma_phi(:) - RHS_gamma_dphi(:) + RHS_alpha_phi(:) - RHS_alpha_dphi(:)
  do j = 1,totPanels
    if (1==RHS(j,2)) then
	  RHStemp(j) = RHS_gamma_phi(j) - RHS_gamma_dphi(j) -0.5*RHS(j,1)
	elseif (0==RHS(j,2)) then
	  RHStemp(j) = 1.*( RHS_alpha_phi(j) - RHS_alpha_dphi(j) )
	else 
	  stop 'RHS does not exist'
	endif
  enddo

  ! solves for {X} in Ax=b through Gaussian inversion {LEGS} or interatively {gmresIter}
  !call LEGS( aN, totPanels, RHStemp, sigma, indxLegs )
  !----
  call gmresIter( aLine, RHStemp, sigma, totPanels)
  !----

  ! Smoothes transition panels out so the neumann/dirichlet corners and extrapolated
  !   back from steady state
  do i=1,totPanels
    if (i<totPanels) then
	  NextPan = i+1
	else
	  NextPan = 1
	endif

    if (RHS(i,2) /= RHS(NextPan,2)) then
	  if (RHS(i,2) == 1) then	
	    tempCoeff = 0.5
	  else 
	    tempCoeff = 0.5
	  endif

	  j=i; LoopTotPan = .TRUE.
  	  PanExtrapLeft: do while ( LoopTotPan == .TRUE. .AND. j>2 )
		j=j-1; tempEpsil = tempCoeff*abs(sigma(j))
	    if (abs(sigma(j-1)-sigma(j)) < tempEpsil .AND. abs(sigma(j)-sigma(j+1)) &
		  < tempEpsil) then	! You're close enough to call it steady state
		  LoopTotPan = .FALSE.
        endif
	  enddo PanExtrapLeft
      j = j+1
	  do while (j<=i)
	    sigma(j)= sigma(j-1)
        j=j+1
	  end do
      
	  if (RHS(NextPan,2) == 1) then	
	    tempCoeff = 0.01
	  else 
	    tempCoeff = 0.05
	  endif

	  j=i-1; LoopTotPan = .TRUE.
	  PanExtrapRight: do while ( LoopTotPan == .TRUE. .AND. j<totPanels-1 )
	    j=j+1; tempEpsil = tempCoeff*abs(sigma(j))
	    if (abs(sigma(j-1)-sigma(j)) < tempEpsil .AND. abs(sigma(j)-sigma(j+1)) &
		  < tempEpsil) then	! You're close enough to call it steady state
		  LoopTotPan = .FALSE.
        endif
  	  enddo PanExtrapRight
	  j=j-1
	  do while (j>=i+1)
	    sigma(j)= sigma(j+1)
	    j=j-1
	  end do
    endif
  enddo
  
  ! converts from the NewCoords back to the original 'xarray, yarray'
  do i=1,sizeLevel
    SizeLevelDistPrev=toolarge
    do j=1,totPanels
	  delY = xarray(i) - xPanCen(j)
	  delX = yarray(i) - yPanCen(j)
      SizeLevelDist = delX**2 + delY**2

	  if (SizeLevelDist<SizeLevelDistPrev) then ! you are at a closer point
        SizeLevelDistPrev=SizeLevelDist
		CloseCoord = j
	  endif
    enddo
	Efield(i) = -epsil_0*PanLen(CloseCoord)* &
	  (NormPanY(CloseCoord) + NormPanX(CloseCoord))* &
	  sigma(CloseCoord)**2
  enddo

  deallocate(NewCoords); deallocate(xPanBegin); deallocate(yPanBegin); 
  deallocate(xPanLen); deallocate(yPanLen); deallocate(aLine)
  deallocate(PanLen); deallocate(xPan); deallocate(yPan); deallocate(PanGaussWeight)
  deallocate(normPanX); deallocate(normPanY); deallocate(RHS); deallocate(RHStemp);
  deallocate(indxLegs); deallocate(xPanCen); deallocate(yPanCen)
  deallocate(RHS_gamma_phi); deallocate(RHS_gamma_dphi)
  deallocate(RHS_alpha_phi); deallocate(RHS_alpha_dphi)
  deallocate(Coords)

  goto 999		! actual program end, just lists formats and errors below
  ! ***********

! FORMAT listings
10  FORMAT(a, f8.2)
11  FORMAT(a, i1, a, e8.2, a, i5)
12  FORMAT(a, e8.3)
13  FORMAT(3e14.3)
14  FORMAT(f8.0, a, f8.2, a)
15  FORMAT(4e12.3)
16  FORMAT(2i6, e14.3)
17  FORMAT(a, f8.3, a, i4, a, i4)

! ERROR listings

999 end SUBROUTINE MAT_Greens

end MODULE MAT_Greens16
