MODULE PtFitting

! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   1 Aug 2005 begun

	USE ShadingGlobal
	IMPLICIT NONE

! used to restrict variables to prevent accidental calling outside
	PRIVATE      ! everything defaults to 'private' except those expressly labeled otherwise.
	PUBLIC		:: CurveFit, gaussvals, GAUSS, TrapInit, LEGS, Nexti, CalcPDF, DoPtsCross, &
		polygon_area_2d_2, PerturbWave

CONTAINS 

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
! CALLING PROGRAM: BoundaryForceCalc
	SUBROUTINE CurveFit ( NumPanels, SplinePts, knotsx, knotsy, hX, Numknots, &
			PanelVals, whichblob, flatyes, NumBlob, MaxNumBlob)

		! incoming variables
		INTEGER, INTENT(in)							:: Numknots, flatyes, NumBlob, MaxNumBlob
		  ! how many known points, droplet or surface/boundary, blob #, # of non-boundary blobs
		REAL*8, DIMENSION(Numknots), INTENT(in)		:: knotsx, knotsy
		REAL*8, DIMENSION(Numknots,3), INTENT(in) 	:: whichblob
		REAL*8, DIMENSION(Numknots,2), INTENT(in)	:: PanelVals

		! outsourced variables
		INTEGER, INTENT(inout)						:: NumPanels, SplinePts
		REAL*8, DIMENSION(Numknots), INTENT(inout)	:: hX

		! subroutine entirely variables
		REAL*8, DIMENSION(Numknots)					:: totDist
		REAL*8, DIMENSION(Numknots,4)				:: Sx, Sy
		REAL*8, DIMENSION(:,:), ALLOCATABLE			:: temp7

		CONTINUE ! body of program
   		!*******
		! the final equation will look like Si[x]=ai+bi(x-xi)+ci(x-xi)^2+di(x-xi)^3

		! finds the spline coefficient fits (linear or cubic) for the knotsx/knotsy points 
		!   using x vs. y. Use only Sx
		!   NOTE: 'totDist' is meaningless in the first iteration
		CALL LinearCubic ( Sx, Sy, knotsx, knotsy, hX, Numknots, totDist, 1, flatyes )

		! using Sx, calculates the total distance around the shape and the corresponding
		!   number of total panels for a balance between accuracy and speed. Use totDist, SplinePts
		totDist(:)=0; SplinePts=0
		CALL totDistCalc ( totDist, SplinePts, knotsx, knotsy, Numknots, flatyes, NumBlob, MaxNumBlob)

		! using totDist, finds the spline coefficient fits (linear or cubic) for the 
		!   knotsx/knotsy points using S vs. x and S vs. y. Use new Sx, Sy
		CALL LinearCubic ( Sx, Sy, knotsx, knotsy, hX, Numknots, totDist, 0, flatyes )

		! using Sx, Sy for the new curve, returns the evenly space points SplineCoords
		!   along the surface, Finally(!) get 'NewCoords' we desire
		ALLOCATE(temp7(SplinePts,8)); temp7(:,:)=-1
		CALL SplineLengthFit (temp7, Sx, Sy, SplinePts, PanelVals, totDist, Numknots, &
			whichblob, flatyes)
		NewCoords(NumPanels+1:NumPanels+SplinePts,:)=temp7
		DEALLOCATE(temp7)

		! update the total number of panels
		NumPanels=NumPanels+SplinePts
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
! CALLING PROGRAM: CurveFit
	SUBROUTINE LinearCubic ( Sx, Sy, knotsx, knotsy, hX, Conesize, totDist, &
			XYFit, flatyes )

		! incoming variables
		INTEGER, INTENT(in)							:: Conesize, XYFit, flatyes
		  ! # known pts, cubic spline fit?, which iteration, droplet or needle/boundary
		REAL*8, DIMENSION(Conesize), INTENT(in) 	:: knotsx, knotsy, totDist

		! outsourced variables
		REAL*8, DIMENSION(Conesize), INTENT(inout) 	:: hX
		REAL*8, DIMENSION(ConeSize,4),INTENT(inout) :: Sx, Sy

		! subroutine entirely variables
		INTEGER										:: i, CubicXYFit, iPlus
		REAL*8, DIMENSION(ConeSize)	:: a, b, c, d, RHSCubic, aDiag, bDiag, &
			cDiag, distChange

		CONTINUE ! body of program
		!*******
		! the final equation will look like Si[x]=ai+bi(x-xi)+ci(x-xi)^2+di(x-xi)^3
		! setting up matricies and constants for 
		CubicXYFit = 10*XYFit + CubicSplineFit

		DO i=1,Conesize
			iPlus=Nexti(i,Conesize,flatyes)
			IF (iPlus==i) THEN
				hX(i)=knotsx(i)-knotsx(1)
			ELSE
				! don't change hX
			ENDIF
		ENDDO

		SELECT CASE (CubicXYFit)
		CASE (11)				! First loop, x vs y, cubic interpolation
			!-------------------------------------------------------------------
			RHSCubic(1:ConeSize)=0

			a(1)=knotsy(1)
			aDiag(1)=0; bDiag(1) = 1.; cDiag(1)=0

			DO i=2,Conesize
				a(i)=1.*knotsy(i)
				aDiag(i)=1.*hX(i-1)/3.
				bDiag(i)=2.*hX(i-1)/3. + 2.*hX(i)/3.
				cDiag(i)=1.*hX(i)/3.
			ENDDO

			RHSCubic(1:ConeSize)=0
			do i=2,Conesize
				iPlus=Nexti(i,Conesize,flatyes)
				if ( abs(hX(i)) > smallepsil .AND. abs(hX(i-1))> smallepsil .AND. a(iPlus)-a(i) > epsil/10. ) then
					RHSCubic(i)=1.*(a(iPlus)-a(i))/hX(i) - 1.*(a(i)-a(i-1))/hX(i-1)
				endif
			enddo

			c(1:ConeSize)=0
			call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

			b(1:ConeSize)=0; d(1:ConeSize)=0
			do i=1,Conesize
				if ( abs(hX(i)) > smallepsil ) then
					iPlus=Nexti(i,Conesize,flatyes)
					b(i) = 1.*( a(iPlus)-a(i) )/hX(i) - (2.*hX(i)*c(i)/3. + hX(i)*c(iPlus)/3.)
					d(i) = 1.*( c(iPlus)-c(i) )/(3.*hX(i))
				endif
			enddo

			Sx(1:ConeSize,1)=a(1:ConeSize); Sx(1:ConeSize,2)=b(1:ConeSize)
			Sx(1:ConeSize,3)=c(1:ConeSize); Sx(1:ConeSize,4)=d(1:ConeSize)

		case (10)				! First loop, x vs y, linear interpolation
			!-------------------------------------------------------------------
			Sx(1:ConeSize,3)=0; Sx(1:ConeSize,4)=0

			do i=1,ConeSize
				Sx(i,1)=1.*knotsy(i)
			enddo
			do i=1,Conesize
				iPlus=Nexti(i,Conesize,flatyes)

				if ( abs(hX(i)) > smallepsil ) then
					Sx(i,2) = 1.*( Sx(iPlus,1)-Sx(i,1) )/hX(i)
				else
					Sx(i,2) = 0
				endif
			enddo

		case (01)				! Second loop, x,y vs s, cubic interpolation
			!-------------------------------------------------------------------
			distChange(1) = 1.*totDist(1)
			a(1)=1.*knotsx(1); 
			aDiag(1)=0; bDiag(1) = 1.; cDiag(1)=0

			do i=2,Conesize
				distChange(i) = 1.*(totDist(i)-totDist(i-1))

				a(i)=1.*knotsx(i)
				aDiag(i)=distChange(i-1)/3.
				bDiag(i)=2.*distChange(i-1)/3. + 2.*distChange(i)/3.
				cDiag(i)=distChange(i)/3.
			enddo

			RHSCubic(1:ConeSize)=0
			do i=2,Conesize
				iPlus=Nexti(i,Conesize,flatyes)
				if ( abs(distChange(i)) > smallepsil .AND. a(iPlus)-a(i) > epsil/10. ) then
					RHSCubic(i)=1.*(a(iPlus)-a(i))/distChange(i) - 1.*(a(i)-a(i-1))/distChange(i-1)
				endif
			enddo

			c(1:ConeSize)=0
			call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

			b(1:ConeSize)=0; d(1:ConeSize)=0
			do i=1,Conesize
				if ( abs(distChange(i)) > smallepsil ) then
					iPlus=Nexti(i,Conesize,flatyes)
					b(i) = 1.*( a(iPlus)-a(i) )/distChange(i) - (2.*distChange(i)*c(i)/3. + distChange(i)*c(iPlus)/3.)
					d(i) = ( c(iPlus)-c(i) )/(3.*distChange(i))
				endif
			enddo

			Sx(1:ConeSize,1)=a(1:ConeSize); Sx(1:ConeSize,2)=b(1:ConeSize)
			Sx(1:ConeSize,3)=c(1:ConeSize); Sx(1:ConeSize,4)=d(1:ConeSize)
 			!******
  			! Same section, now for y vs s.
			do i=1,ConeSize
				a(i)=1.*knotsy(i)
			enddo

			RHSCubic(1:ConeSize)=0.
			do i=2,Conesize
				if ( abs(distChange(i)) > epsil .AND. ABS(a(i+1)-a(i)) > epsil/10. ) then
					RHSCubic(i)=1.*(a(i+1)-a(i))/distChange(i) - 1.*(a(i)-a(i-1))/distChange(i-1)
				endif
			enddo

			c(1:ConeSize)=0
			call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

			b(1:ConeSize)=0; d(1:ConeSize)=0
			do i=1,Conesize
				if ( abs(distChange(i)) > epsil ) then
					iPlus=Nexti(i,Conesize,flatyes)
					b(i) = 1.*( a(iPlus)-a(i) )/distChange(i) - (2.*distChange(i)*c(i)/3. + distChange(i)*c(iPlus)/3.)
					d(i) = 1.*( c(iPlus)-c(i) )/(3.*distChange(i))
				endif
			enddo

			Sy(1:ConeSize,1)=a(1:ConeSize); Sy(1:ConeSize,2)=b(1:ConeSize)
			Sy(1:ConeSize,3)=c(1:ConeSize); Sy(1:ConeSize,4)=d(1:ConeSize)

		case (00)				! Second loop, x,y vs s, linear interpolation
			!-------------------------------------------------------------------
			Sx(:,3)=0; Sx(:,4)=0
			Sy(:,3)=0; Sy(:,4)=0

			do i=1,ConeSize
				if (i /= 1) then
					distChange(i) = 1.*( totDist(i)-totDist(i-1) )
				else
					distChange(i) = 1.*totDist(i)
				endif

				Sx(i,1)=1.*knotsx(i)
				Sy(i,1)=1.*knotsy(i)
			enddo

			do i=1,Conesize
				if ( abs(distChange(i)) > smallepsil ) then
					iPlus=Nexti(i,Conesize,flatyes)
					Sx(i,2) = 1.*( Sx(iPlus,1)-Sx(i,1) )/distChange(i)
					Sy(i,2) = 1.*( Sy(iPlus,1)-Sy(i,1) )/distChange(i)
				else 
					Sx(i,2)=0; Sy(i,2)=0
				endif
			enddo
		case default
			stop 'CubicSpline: Error on CubicXYFit'
		end select ! switch CubicXYFit

	END SUBROUTINE LinearCubic

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

	SUBROUTINE SplineLengthFit (SplineCoords, Sx, Sy, SplinePts, PanelVals, totDist, ConeSize, &
			whichblob1, flatyes2)

		! incoming variables
		integer, intent(in)							:: SplinePts, ConeSize, flatyes2
		  ! CubicSplineFit = 1 for cubic interpolation, = 0 for linear interpolation
		  ! 'SplinePts'=# of new coords. 'ConeSize'=# of old coords
		real*8, dimension(ConeSize,5), intent(in)	:: Sx, Sy
		real*8, dimension(ConeSize,3),intent(in)	:: whichblob1
		real*8, dimension(ConeSize,2), intent(in)	:: PanelVals
		real*8, dimension(ConeSize), intent(in)		:: totDist

		! outsourced variables
		real*8,dimension(SplinePts,8),intent(inout) :: SplineCoords

		! subroutine entirely variables
		integer										:: tempJ, i, j
		real*8										:: temptotDist, tempfracDist

		real*8, dimension(SplinePts)				:: SplineLength

		continue
		tempJ=1; SplineLength(:)=0
		do i=1,SplinePts
			SplineLength(i) = 1.*(i-1)*totDist(ConeSize)/(SplinePts-1)
			IF (tempJ /=1) THEN
				j = tempJ-1
			ELSE
				j = 1
			ENDIF

			ChoosePanel: DO WHILE ( totDist(j) < SplineLength(i) ) 
				j=j+1
				IF (j>ConeSize) THEN
					j=ConeSize
					EXIT
				ENDIF
			ENDDO ChoosePanel
			tempJ = j

			IF (1==j) THEN
				temptotDist = 0
			ELSE
				temptotDist = 1.*totDist(j-1)
			ENDIF

			tempfracDist = SplineLength(i)-temptotDist
			SplineCoords(i,1) = Sx(j,1) + Sx(j,2)*(tempfracDist) + &
				Sx(j,3)*(tempfracDist**2) + Sx(j,4)*(tempfracDist**3)
			SplineCoords(i,2) = Sy(j,1) + Sy(j,2)*(tempfracDist) + &
				Sy(j,3)*(tempfracDist**2) + Sy(j,4)*(tempfracDist**3)

			SplineCoords(i,3) = PanelVals(j,1) ! panel value
			SplineCoords(i,4) = PanelVals(j,2) ! panel type
			SplineCoords(i,5) = whichblob1(j,1)	! knot {j} is in blob {i}
			SplineCoords(i,6) = whichblob1(j,2)	! knot {j} has blob_xc {i}
			SplineCoords(i,7) = whichblob1(j,3)	! knot {j} has blob_yc {i}
			SplineCoords(i,8) = flatyes2		! surface is open(flat)=1 or droplet=0
		enddo
	END SUBROUTINE SplineLengthFit

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
		real*8, dimension(N), intent(inout) :: u

		! subroutine entirely variables
		integer								:: j
		real*8								:: bet
		real*8, dimension(N)				:: gam

		continue ! body of program
		!*******
		if (0.0==centerrow(1)) then
			stop 'CubicSpline: Error 1 in tridiag'
		endif

		bet=centerrow(1)
		u(1)=rhsTri(1)/bet

		do j=2,N		! decomposition and forward substitution
			gam(j) = rightrow(j-1)/bet
			bet=centerrow(j)-leftrow(j)*gam(j)
			if (0==bet) then
				if (centerrow(j)==0 .AND. leftrow(j)==0 .AND. gam(j)==0 .AND. &
					rhsTri(j)==0) then
					u(j)=0
				else
					print *, 'CubicSpline: Error 2 in tridiag. j=', j
					print *, '  centerrow(j), leftrow, gam=', centerrow(j), leftrow(j), gam(j)
					print *, '  rhs(j)=', rhsTri(j)
					stop 
				endif
			else
				u(j)=( rhsTri(j)-leftrow(j)*u(j-1) )/bet
			endif
		enddo

		do j=N-1,1,-1	! back substitution
			u(j) = u(j) - gam(j+1)*u(j+1)
		enddo

	END SUBROUTINE Tridiag

! ***************
! Subroutine totDistCalc
! ***************
! DESCRIPTION: given a spline fitting coefficients, numerically integrates
!   along surface to get an accurate total distance around shape
!
! INPUT: knot locations, cubic spline coeffs, # of knots along test shape and 
!	overall
!
! OUTPUT: total distance around entire border, # of points needed from global distance param
!
! CALLING PROGRAM: CurveFit

	SUBROUTINE totDistCalc ( totDist1, SplinePts1, knotsx1, knotsy1, ConeSize1, flatyes3, &
			NumBlob3, MaxNumBlob3)
		! incoming variables
		INTEGER, INTENT(in)							:: ConeSize1, flatyes3, NumBlob3, MaxNumBlob3
		REAL*8, DIMENSION(ConeSize1), INTENT(in)	:: knotsx1, knotsy1

		! outsourced variables
		INTEGER, INTENT(out)						:: SplinePts1
		REAL*8, DIMENSION(ConeSize1), INTENT(inout)	:: totDist1

		! subroutine entirely variables
		INTEGER										:: i, iPlus
		REAL*8										:: tempArcDist, one
		REAL*8, DIMENSION(ConeSize1)				:: xPanLen, yPanLen, dydx

		CONTINUE		
		tempArcDist=0.0; one=1.0

		! dePanel = floor(1.*(i-1)/M)+1
		DO i = 1, ConeSize1
			iPlus=Nexti(i,Conesize1,flatyes3)

			xPanLen(i) = knotsx1(iPlus)-knotsx1(i)
			yPanLen(i) = knotsy1(iPlus)-knotsy1(i)

			IF ( ABS(xPanLen(i)) > smallepsil ) THEN
				dydx(i) = yPanLen(i) / xPanLen(i)
				tempArcDist = ABS(xPanLen(i))*SQRT( 1+dydx(i)*dydx(i) )
			ELSE
				tempArcDist = ABS(yPanLen(i))
			ENDIF

			IF (i /= 1) THEN
				totDist1(i) = totDist1(i-1) + tempArcDist
			ELSE
				totDist1(i) = tempArcDist
			ENDIF
		ENDDO

		! given the total distance around, calculates the number of points needed
		!   'MinNumPts_blob' is global
		IF ( NumBlob3 > MaxNumBlob3 ) THEN	! have fewer points around the boundary
			SplinePts1 = CEILING(MAX( FLOOR(totDist1(ConeSize1)/maxDistStep)+one, one*MinNumPtsBlob ))
		ELSE
			SplinePts1 = CEILING(MAX( StepMult*FLOOR(totDist1(ConeSize1)/maxDistStep)+one, one*MinNumPtsBlob ))
		ENDIF
	END SUBROUTINE totDistCalc

! ***************
! SUBROUTINE CalcPDF
! ***************
! DESCRIPTION: Once the first blob has formed, calculate snapoff PDF
!
! INPUT: <>
!
! OUTPUT: Droplet charge around surface
!
! CALLING PROGRAM: Shading:BFC

	SUBROUTINE CalcPDF ( totCharge, DropLen, DropEn, sizeDrop )
		! incoming variables
		INTEGER, INTENT(in)			:: sizeDrop
		REAL*8, DIMENSION(sizeDrop), INTENT(in)	:: DropLen, DropEn

		! outgoing variables
		REAL*8, INTENT(out)	:: totCharge	! in [C]

		! subroutine entirely variables
		INTEGER						:: i
		REAL*8				:: temptotCharge

		CONTINUE
		totCharge=0;

		DO i=1,sizeDrop
			temptotCharge=DropEn(i)*(DropLen(i)/100.0)
			totCharge=totCharge+temptotCharge
		ENDDO
		totCharge=ABS(epsil0*totCharge)
	END SUBROUTINE CalcPDF

! ***************
! SUBROUTINE PerturbWave
! ***************
! DESCRIPTION: Takes a set original plot and adds a perpendicular low-amplitude sine wave to it
!
! INPUT: smooth 'phi'
!
! OUTPUT: perturbed phi
!	
! Author: Anton VanderWyst, with heavy input from Matt McNenly
!
! CALLING PROGRAM: Shading:BFC

	SUBROUTINE PerturbWave ( NewPhi, x, y, MaxX, MaxY, dx, dy )	
		! incoming
		INTEGER, INTENT(in) 	:: MaxX, MaxY
		REAL*8, INTENT(in) 		:: dx, dy

		! outgoing
		REAL*8,DIMENSION(MaxX, MaxY),INTENT(inout) 	:: NewPhi
		REAL*8, DIMENSION(MaxX), INTENT(out) 		:: x
		REAL*8, DIMENSION(MaxY), INTENT(out) 		:: y

		! internal parameters
		REAL*8, PARAMETER		:: PerturbAmp=0.02, PerturbFreq=17.0, xblob=0.5, &
			yblob=0.5, radblob=0.25

		! internal variables
		INTEGER					:: i, j
		REAL*8					:: CosConst, NormF, NormX, dddx, DispH, DispX, DispD
		REAL*8, DIMENSION(MaxY) :: S
		REAL*8, DIMENSION(MaxX, MaxY) 	:: d

		CONTINUE 

		CosConst = 2.0*Pi/xblob
		S(:)=0; NewPhi(:,:)=0
		DO i=1,MaxX
			x(i) = (i+0.5)*dx
			DO j=1,MaxY
				IF (1==i) THEN
					y(j) = (j+0.5)*dy
				ENDIF
				d(i,j)=yblob+radblob*cos(x(i)*CosConst)-y(j)

				IF (1==j) THEN
					dddx=CosConst*sin( CosConst*x(i) )
					NormX=dddx/SQRT( 1+dddx*dddx )
					NormF=-1.0/SQRT( 1+dddx*dddx )
				ENDIF
				IF (j>1) THEN
					S(j)=S(j-1)+SQRT( (xblob/(MaxX-1))*(xblob/(MaxX-1)) + &
						(d(i,j)-d(i,j-1))*(d(i,j)-d(i,j-1)) )
				ENDIF

				DispH = PerturbAmp*cos(PerturbFreq*S(j))
				DispX = x(i) + DispH*NormX
				DispD = y(j) + DispH*NormF

				NewPhi(i,j)= yblob+radblob*cos(DispX*CosConst)-DispD;
			ENDDO
		ENDDO
	END SUBROUTINE PerturbWave

! ***************
! SUBROUTINE POLYGON_AREA_2D_2
! ***************
! DESCRIPTION: Computes the area of a polygon in 2D. The area is the sum of the 
! areas of the triangles formed by node N with consecutive pairs of nodes.
!
! INPUT: integer: N, the number of vertices of the polygon. real:V(2,N), the vertices.
!
! OUTPUT: real: AREA, the absolute area of the polygon.
!	
! Author: John Burkardt
! http://www.csit.fsu.edu/~burkardt/f_src/geometry/geometry.html
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
! CALLING PROGRAM: Shading:BFC

	SUBROUTINE polygon_area_2d_2 ( areaPoly, v, n )
		! incoming
		INTEGER, INTENT(in) 	:: n
		REAL*8, DIMENSION(n,2), INTENT(in) :: v

		! outgoing
		REAL*8, INTENT(out) 	:: areaPoly

		! internal
		INTEGER					:: i
		REAL*8					:: areat
		REAL*8, DIMENSION(3,2) 	:: t

		CONTINUE 
		areaPoly = 0.0

		DO i = 1, n - 2
			!t = RESHAPE ( (/ v(1:2,i), v(1:2,i+1), v(1:2,n) /), (/ 2, 3 /) )
			t(1,:) = v(i,:)
			t(2,:) = v(i+1,:)
			t(3,:) = v(n,:)

			areat = 0.5 * ABS ( t(1,1) * ( t(2,2) - t(3,2) ) &
				+ t(2,1) * ( t(3,2) - t(1,2) ) + t(3,1) * ( t(1,2) - t(2,2) ) )

			areaPoly = areaPoly + areat
		ENDDO

		RETURN
	END SUBROUTINE polygon_area_2d_2

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

	SUBROUTINE TrapInit ( GreenPart, delX, delY, normXj, normYj, normXi, normYi, &
			PanWeight, PanLen, PanFromTo, flagTrap, k, Mpts )

		! incoming variables
		INTEGER, INTENT(in)		:: PanFromTo, flagTrap, k, Mpts
		REAL*8, INTENT(in)		:: PanLen, delX, delY, normXj, normYj, normXi, normYi, PanWeight

		! outgoing variables
		REAL*8, INTENT(inout)	:: GreenPart

		! subroutine entirely variables
		REAL*8					:: HalfLen							

		CONTINUE
		GreenPart = 0.

		IF ( 0==flagTrap .OR. 1==k .OR. Mpts==k ) then
			HalfLen = 0.5
		ELSE
			HalfLen = 1.0 
		ENDIF

		IF ( abs(delX)<epsil .AND. abs(delY)<epsil ) then          
    		! do nothing, on the border                
		ELSE
			SELECT CASE (PanFromTo)    ! 1=dirichlet, 0=neumann type BC
			CASE (11)
				GreenPart = PanLen * HalfLen * PanWeight * GradGreens(delX, delY, normXj, normYj)
			CASE (10) 
				GreenPart = PanLen * HalfLen * PanWeight * Greens(delX, delY)
			CASE (0)
				GreenPart = PanLen * HalfLen * PanWeight * GradGreens(delX, delY, normXi, normYi)
			CASE (01)
				GreenPart = PanLen * HalfLen * PanWeight * GradGradGreens(delX, delY, &
					normXj, normYj, normXi, normYi)
			CASE DEFAULT
				PRINT *, 'PtFitting:TrapInit err; PanFromTo=', PanFromTo
				STOP 
			END SELECT ! switch PanFromTo
		ENDIF    ! on border
	END SUBROUTINE TrapInit

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
		real*8				:: GradGradGreens, R, gamma, tempXY

		continue
		R = ( delX**2 + delY**2 )    ! = r^2

		if ( 0 /= delX .AND. 0 /= delY) then 
			gamma = -R/(2*delX*delY)

			GradGradGreens = (-1/(2*gamma*pi*R))*( ( (gamma+delX/delY)*normXj + normYj)*normXi + &
				( (gamma+delY/delX)*normYj + normXj)*normYi )
		elseif ( 0==delX ) then
			tempXY=(1.0/R-2*delY**2/R**2)*normYj
			GradGradGreens = (-1/(2*pi))*((normXj/R)*normXi + tempXY*normYi )
		elseif ( 0==delY ) then
			tempXY=(1.0/R-2*delX**2/R**2)*normXj
			GradGradGreens = (-1/(2*pi))*(tempXY*normXi + (normYj/R)*normYi )
		else
			stop 'Error on gradgradG'
		end if

	end FUNCTION GradGradGreens

! ***************
! SUBROUTINE DoPtsCross
! ***************
! DESCRIPTION: takes 4 data points where the connectivity is known and determines if they 
!   cross. Online=cross.
!
! INPUT: x1y1 -> x2y2, x3y3 -> x4y4
!
! OUTPUT: do they cross (yes/no)
!
! CALLING PROGRAM: Shading:SortNN_NoCross
!
! ALGORITHM: compare sign of ( (ad x ab) v. (ac x ab) ). if same, no cross
	SUBROUTINE DoPtsCross(  DoCross, x1, y1, x2, y2, x3, y3, x4, y4 )
		! incoming variables
		REAL*8, INTENT(in)	:: x1, y1, x2, y2, x3, y3, x4, y4

		! outgoing variables
		LOGICAL,INTENT(OUT) :: DoCross

		! start regular variables
		REAL*8				:: one, ax4, ay4, ax3, ay3, bx, by, k41, k31, sgn4, sgn3

		CONTINUE 
		one=1.0
		ax4=x4-x1; ay4=y4-y1
		ax3=x3-x1; ay3=y3-y1
		bx=x2-x1; by=y2-y1

		k41=ax4*by-bx*ay4
		sgn4=SIGN(one,k41)

		k31=ax3*by-bx*ay3
		sgn3=SIGN(one,k31)

		IF (sgn3==sgn4) THEN	! no cross
			DoCross=.FALSE.
		ELSEIF ( 2*one ==ABS(sgn3-sgn4) ) THEN
			DoCross=.TRUE.
		ELSE
			PRINT *, 'PtFitting:DoPtsCross sign err. sgn3,sgn4=',sgn3,sgn4
			STOP
		ENDIF
	END SUBROUTINE DoPtsCross

! ***************
! SUBROUTINE gaussvals
! ***************
! DESCRIPTION: Calculates the Gaussian quadrature points and weights
!
! INPUT: # points / panel, panel beginning and end points
!
! OUTPUT: location of Gaussian quadrature points along panel
!
! CALLING PROGRAM: 
!
! 	20 Jun 2005 using Jerry Emhoff's automatic number generation code

	SUBROUTINE gaussvals(points, weights)
		! incoming variables - only global number of g. points
		! outgoing variables
		REAL*8, DIMENSION(M), INTENT(OUT) :: points, weights

		! start regular variables
		INTEGER					:: i, j, ind, counter
		REAL*8					:: error, oldval
		REAL*8, DIMENSION(M) 	:: beta
		REAL*8, DIMENSION(M*M) 	:: z, q, r

		CONTINUE 
		error=1.; oldval=2.; z=0.

		DO i=2,M
			beta(i-1)=0.5/sqrt(1.0-1.0/(4.0*(i-1)*(i-1)))
		ENDDO

		DO i=1,M-1
			z((i-1)*M+i+1)=beta(i)
			z(i*M+i)=beta(i)
		ENDDO

		counter=0
		DoErr: DO WHILE (error > 1e-14)  
			counter=counter+1
			CALL mgs(z, q, r)			! Compute the QR factorization of z

			DO i=1,M				! Multiply a and q to get z
				DO j=1,M
					z((i-1)*M+j)=0.0

					IF (i==1) THEN
						z((i-1)*M+j) = z((i-1)*M+j)+ beta(i)*q(i*M+j)
					ELSEIF (i==M) THEN
						z((i-1)*M+j) = z((i-1)*M+j) + beta(i-1)*q((i-2)*M+j)
					ELSE
						z((i-1)*M+j) = z((i-1)*M+j) + beta(i-1)*q((i-2)*M+j) + &
							beta(i)*q(i*M+j)
					ENDIF
				ENDDO
			ENDDO

			! Compute change in first eigenvalue
			error=ABS((r(1)-oldval)/oldval)
			oldval=r(1)
		ENDDO DoErr

		! Assign points from the diagonal of r (the eigenvalues)
		ind=1
		DO i=1,M,2 
			points(ind) = -r((i-1)*M+i)
			points(M-ind+1) = r((i-1)*M+i)
			ind = ind+1
		ENDDO		

		! Assign weights from the square of the non-zero leading values in the 
		!   columns of q (the eigenvectors)
		ind=1
		DO i=1,M
			IF (q(i) /= 0.0) THEN
				weights(ind)=q(i)*q(i)
				weights(M-ind+1)=q(i)*q(i)
				ind = ind+1
			ENDIF
		ENDDO

		! Take care of the odd numbered case
		!sum=0.0
		IF (MOD(M,2) == 1) THEN
			PRINT *, 'Err: AxiSym:gaussvals odd num'
			STOP
		!points(M/2)=0.0
		!
		!DO i=1,M/2
		!sum = sum + weights(i)
		!ENDDO
		!weights(M/2)=2.0-2.0*sum
		ENDIF
	END SUBROUTINE gaussvals

! ***************
! SUBROUTINE mgs
! ***************
! DESCRIPTION: Computes QR factorization
! CALLING PROGRAM: gaussvals
! Primary author: J. Emhoff, University of Michigan 
!
! Version history:
!   21 Jun 2005 begun 
	SUBROUTINE mgs(a, q, r)
		! incoming variables - only global number of g. points
		REAL*8, DIMENSION(M*M), INTENT(in) 	:: a

		! outgoing variables
		REAL*8, DIMENSION(M*M), INTENT(OUT) :: q, r

		! start regular variables
		INTEGER								:: i, j, k, n
		REAL*8, DIMENSION(M*M) 				:: v

		CONTINUE 
		n=M; v=a; r=0;

		ido: DO i=1,n
			r((i-1)*n+i)=0
			DO k=1,n
				r((i-1)*n+i) = r((i-1)*n+i) + v((k-1)*n+i)*v((k-1)*n+i)
			ENDDO
			r((i-1)*n+i)=sqrt(r((i-1)*n+i))		! Calculate rjj

			DO k=1,n					! Calculate q
				q((k-1)*n+i)=v((k-1)*n+i)/r((i-1)*n+i)
			ENDDO

			DO j=i+1,n
				DO k=1,n				! Do the dot product...
					r((i-1)*n+j) = r((i-1)*n+j) + q((k-1)*n+i)*v((k-1)*n+j)
				ENDDO

				DO k=1,n				! Calculate v
					v((k-1)*n+j) = v((k-1)*n+j) - r((i-1)*n+j)*q((k-1)*n+i)
				ENDDO
			ENDDO
		ENDDO iDo
	END SUBROUTINE mgs

! Updated 10/24/2001. Please Note:                                      !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997. 
!
! http://www.physics.unlv.edu/~pang/comp3/code43.f90
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! An example of solving linear equation set A(N,N)*X(N) = B(N)
! with the partial-pivoting Gaussian elimination scheme. 
! Copyright (c) Tao Pang 2001.

	SUBROUTINE LEGS (A,N,B,X,INDX)

		INTEGER, INTENT (IN) :: N
		INTEGER :: I,J
		INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
		REAL*8, INTENT (INOUT), DIMENSION (N,N) :: A
		REAL*8, INTENT (INOUT), DIMENSION (N) :: B
		REAL*8, INTENT (INOUT), DIMENSION (N) :: X

		CONTINUE 

		CALL ELGS (A,N,INDX)

		DO I = 1, N-1
			DO J = I+1, N
				B(INDX(J)) = B(INDX(J))-A(INDX(J),I)*B(INDX(I))
			END DO
		END DO
		X(N) = B(INDX(N))/A(INDX(N),N)

		IF (ABS(X(N))>1e15) THEN
			PRINT *, 'PtFitting:LEGS X(N) NaN/Inf err. X(N),N=', X(N), N
			STOP
		ENDIF

		DO I = N-1, 1, -1
			X(I) = B(INDX(I))
			DO J = I+1, N
				X(I) = X(I)-A(INDX(I),J)*X(J)
			END DO
			X(I) =  X(I)/A(INDX(I),I)
		END DO
	END SUBROUTINE LEGS

	SUBROUTINE ELGS (A,N,INDX)

! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.

		INTEGER, INTENT (IN) 	:: N
		INTEGER 				:: I, J, K, L, ITMP
		INTEGER, INTENT (OUT), DIMENSION (N) 	:: INDX
		REAL*8 					:: C1,PI,PI1,PJ
		REAL*8, INTENT (INOUT), DIMENSION (N,N) :: A
		REAL*8, DIMENSION (N) 	:: C

		CONTINUE 

		! Initialize the index
		DO I = 1, N
			INDX(I) = I
		END DO

		! Find the rescaling factors, one from each row
		DO I = 1, N
			C1= 0.0
			DO J = 1, N
				C1 = DMAX1(C1,ABS(A(I,J)))
			END DO
			C(I) = C1
			IF (0==C1) THEN
				PRINT *, 'PtFitting:ELGS C(I)=0/NaN err. i=', I
				STOP
			ENDIF
		END DO

		! Search the pivoting (largest) element from each column
		DO J = 1, N-1
			PI1 = 0.0; K=0
			DO I = J, N
				PI = ABS(A(INDX(I),J))/C(INDX(I))

				IF (PI > PI1) THEN
					PI1 = PI
					K   = I
				ENDIF
			END DO
			IF ( 0==K ) THEN
				PRINT *, 'PtFitting:ELGS no swap. j,n=',j,n, A(INDX(N-1),849),A(INDX(N),849)
				PRINT *, '  INDX(I)=',INDX(849),INDX(850)
				STOP
			ENDIF

			! Interchange the rows via INDX(N) to record pivoting order
			ITMP    = INDX(J)
			INDX(J) = INDX(K)
			INDX(K) = ITMP

			DO I = J+1, N
				PJ  = A(INDX(I),J)/A(INDX(J),J)

				! Record pivoting ratios below the diagonal
				A(INDX(I),J) = PJ

				! Modify other elements accordingly
				DO L = J+1, N
					A(INDX(I),L) = A(INDX(I),L)-PJ*A(INDX(J),L)

					IF ( 0==A(INDX(I),L)) THEN
						PRINT *, 'zero A. i,j,l,pj,indx(i)=',I,J,L,PJ,INDX(I)
					ENDIF
				END DO
			END DO
		END DO
	END SUBROUTINE ELGS

! ***************
! Function Nexti
! ***************
! DESCRIPTION: gives the next coordinate for the surface if it is open(flat) as
!   the needle/electrodes are or rounded/closed as the droplets. 
!
! INPUT: i, # panels, if open/closed
!
! OUTPUT: next i

	FUNCTION Nexti( i1, numknots, flatyes1, currentpt2 )
		INTEGER, INTENT(in)	:: i1, numknots, flatyes1
		INTEGER, OPTIONAL, INTENT(in)	:: currentpt2

		! function
		INTEGER				:: Nexti

		! internal
		INTEGER				:: flagcurrentNull, tempcurrent
		CONTINUE

		! changes to summing from the beginning instead of just within the blob
		IF ( PRESENT(currentpt2) ) THEN
			flagcurrentNull=0	! sum all blobs
			tempcurrent=currentpt2
		ELSE
			flagcurrentNull=1	! sum blob only
			tempcurrent=0
		ENDIF

		IF (1==flatyes1) THEN
			! needle surface or boundary. Leave last pt free
			IF (numknots==i1) THEN
				Nexti=i1*flagcurrentNull+(1-flagcurrentNull)*tempcurrent
			ELSE
				Nexti=i1*flagcurrentNull+(1-flagcurrentNull)*tempcurrent+1
			ENDIF
		ELSEIF (0==flatyes1) THEN	
			! droplet, with wraparound needed from last-> first
			IF (numknots==i1) THEN
				Nexti=flagcurrentNull+(1-flagcurrentNull)*(tempcurrent-numknots+1)
			ELSE
				Nexti=i1*flagcurrentNull+(1-flagcurrentNull)*tempcurrent+1
			ENDIF
		ELSE
			PRINT *, 'PtFitting:Nexti err; i1, numknots, flatyes1, currentpt2=', i1, numknots, flatyes1, currentpt2
			STOP
		ENDIF
	end FUNCTION Nexti

!**GAUSS****************************************************************
! Subroutine to find solution of a linear system of N equations in N   *
! unknowns using Gaussian elimination, provided a unique solution      *
! exists.  The coefficients and constants of the linear system are     *
! stored in the matrix LIN, which has N rows and N+1 columns.		   *
! If the system is singular, SINGUL is returned as true, and the       *
! solution X is undefined.  Local identifiers used are:                *                        
!     I,J,K  : subscripts                                              *
!     MULT   : multiplier used to eliminate an unknown                 *
!     ABSPIV : absolute value of pivot element                         *
!     PIVROW : row containing pivot element                            *
!     EPSIL  : a small positive real value ("almost zero")             *
!     TEMP   : used to interchange rows of matrix                      *
!																	   *
! Accepts: Two-dimensional array A, integer LIMROW				   	   *
! Returns: One-dimensional array X           						   *
!***********************************************************************

	SUBROUTINE GAUSS(A, LIMROW, X, B)

		! incoming variables
		INTEGER, INTENT(in)						:: LIMROW
		REAL*8, DIMENSION(LIMROW, LIMROW), INTENT(in) :: A
		REAL*8, DIMENSION(LIMROW), INTENT(in)	:: B

		! outsourced variables
		REAL*8, DIMENSION(LIMROW),INTENT(inout) :: X

		INTEGER									:: PIVROW, I, J, K
		REAL*8, PARAMETER						:: EPSIL=1e-15
		REAL*8, DIMENSION(LIMROW, LIMROW+1)		:: LIN
		REAL*8									:: TEMP, MULT, ABSPIV

		DO I = 1, LIMROW
			DO J = 1, LIMROW
				LIN(I,J) = A(I,J)
			END DO
			LIN(I,LIMROW+1) = B(I)
		END DO

		ILoop: DO I = 1, LIMROW ! Locate pivot element
			ABSPIV = abs(LIN(I,I))
			PIVROW = I
			DO K = I + 1, LIMROW
				IF (ABS(LIN(K,I)) > ABSPIV) THEN
					ABSPIV = ABS(LIN(K,I))
					PIVROW = K
				ENDIF
			END DO

			! Check if matrix is (nearly) singular
			IF (ABSPIV < EPSIL) THEN
				PRINT *, '** PtFitting:Gauss err; nearly singular matrix. i, lin(i,i)=', I, LIN(I,I)
				STOP
			ENDIF

			! It isn't, so interchange rows PIVROW and I if necessary
			IF (PIVROW /= I) THEN
				DO J = 1, LIMROW + 1
					TEMP = LIN(I,J)
					LIN(I,J) = LIN(PIVROW,J)
					LIN(PIVROW,J) = TEMP
				END DO
			ENDIF

			! Eliminate Ith unknown from equations I + 1, ..., N
			JLoop: DO J = I + 1, LIMROW
				MULT = -LIN(J,I) / LIN(I,I)
				KLoop: DO K = I, LIMROW + 1
					LIN(J,K) = LIN(J,K) +  MULT * LIN(I,K)
				END DO KLoop
			END DO JLoop
		END DO ILoop

		! Find the solutions by back substitution
		X(LIMROW) = LIN(LIMROW, LIMROW + 1) / LIN(LIMROW,LIMROW)
		DO J = LIMROW - 1, 1, -1
			X(J) = LIN(J, LIMROW + 1)
			DO K = J + 1, LIMROW
				X(J) = X(J) - LIN(J,K) * X(K)
			END DO
			X(J) = X(J) / LIN(J,J)
		END DO
	END SUBROUTINE GAUSS
END MODULE PtFitting
