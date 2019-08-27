MODULE Panels
	USE global, ONLY	: CaseShape, NewCoords, epsil0, pi, SmallInitAry, npoints, &
		PotMeshSize, MaxNumInpotIter, cmtom, AformDouble, toolarge, eps, CalcPotBack, fluidPot, &
		ShowInLoopCoords, elecPotential, recastDist, ShowGMResConverge, smallepsil, &
		maxDistStep, StepMult, MinNumPtsBlob, CubicSplineFit, PANEL

	IMPLICIT NONE

	INTERFACE
		SUBROUTINE calcpotdnwrap(a, b, c, d, e, f)
			USE global, ONLY	: npoints, PANEL
			INTEGER, INTENT(in)	:: d
			REAL*8, INTENT(in) 	:: b, c
			REAL*8, DIMENSION(npoints), INTENT(in) 	:: e
			TYPE(PANEL), DIMENSION(d), INTENT(in) 	:: f

			REAL*8, INTENT(out)	:: a
		END SUBROUTINE calcpotdnwrap

		SUBROUTINE calcpotwrap(a, b, c, d, e, f)
			USE global, ONLY	: npoints, PANEL
			INTEGER, INTENT(in)	:: d
			REAL*8, INTENT(in) 	:: b, c
			REAL*8, DIMENSION(npoints), INTENT(in) 	:: e
			TYPE(PANEL), DIMENSION(d), INTENT(in) 	:: f

			REAL*8, INTENT(out)	:: a
		END SUBROUTINE calcpotwrap

	END INTERFACE

	! used to restrict variables to prevent accidental calling outside
	PRIVATE      ! everything defaults to 'private' except those expressly labeled otherwise.
	PUBLIC		:: AddBCandAround, CurveFit, Nexti, CalcGridPot, &
		polygon_area_2d_2, gmresIter, CalcDropletPot, GetPotGrid, & 
	ForceToMeters, ZeroorCopyPanel, DoPtsCross, CalcPDF

CONTAINS 

! ***************
! SUBROUTINE AddBCandAround[x]
! ***************
! DESCRIPTION: adds 3rd column (boundary type) and 4th column (boundary value) to
!	pre-sorted data
!
! INPUT: grouped blob x,y surface data 
!
! OUTPUT: additional exterior surface
!
! CALLING PROGRAM: BoundaryForceCalc

	SUBROUTINE AddBCandAround( Coords, sizeLevel1, NumBlobCoords, blobPlus, blobNum, &
			elecHeight, elecWidth, elecEdge, maxR, minR, maxZ, minZ, numBorderCoord, &
			NeedleNum)

		! incoming variables
		INTEGER, INTENT(in)		:: NumBlobCoords, blobNum, numBorderCoord, NeedleNum
		REAL*8, INTENT(in)		:: elecHeight, elecWidth, maxR

		! outgoing variables
		INTEGER, INTENT(out)	:: blobPlus
		INTEGER, DIMENSION(SmallInitAry), INTENT(out) :: sizeLevel1
		REAL*8, INTENT(out)		:: elecEdge, minR, maxZ, minZ
		REAL*8, DIMENSION(NumBlobCoords+numBorderCoord,8), INTENT(inout) :: Coords
		  ! [x, y, pt val, pt type, blob#, blob_rc, blob_zc, flatpan]	

		! subroutine entirely variables
		INTEGER					:: i
		REAL*8					:: elecThick, templow, temphigh, templowY, needlez

		CONTINUE
		!*************
		! main program commands begin
		templow=1e10; templowY=1e10; temphigh=0
		DO i=1,NumBlobCoords
			templow  = MIN(templow, Coords(i,1)) 
			temphigh = MAX(temphigh, Coords(i,1))
			templowY = MIN(templowY, Coords(i,2))

			IF (INT(Coords(i,5))==NeedleNum) THEN
				needlez=Coords(i,2)
			ENDIF
		ENDDO
		minR = templow; 	!maxR = 0.3
		minZ=needlez 								!!minZ =  Coords(NumBlobCoords,2);	
		maxZ = elecHeight
		i=NumBlobCoords
		elecThick = 0.1; elecEdge=0.001
		! 'elecpotential' is a global variable
		! "CaseShape" is also global.

		SELECT CASE(CaseShape)
		CASE(1)				! test case around single box
			Coords(i+1,1)= temphigh+0.005;			Coords(i+1,2)= minZ
			Coords(i+2,1)= maxR;					Coords(i+2,2)= minZ
			Coords(i+3,1)= maxR;					Coords(i+3,2)= 0.7*elecHeight
			Coords(i+4,1)= maxR;					Coords(i+4,2)= 0.8*elecHeight
			Coords(i+5,1)= maxR;					Coords(i+5,2)= 0.9*elecHeight
			Coords(i+6,1)= maxR;					Coords(i+6,2)= elecHeight
			Coords(i+7,1)= 0.8*maxR;				Coords(i+7,2)= elecHeight
			Coords(i+8,1)= 0.7*maxR;				Coords(i+8,2)= elecHeight
			Coords(i+9,1)= 0.6*maxR;				Coords(i+9,2)= elecHeight
			Coords(i+10,1)= 0;						Coords(i+10,2)= elecHeight

			! debugging, test shape, force to all Dirichlet to check code bit by bit
			Coords(i+1:i+numBorderCoord,3)=elecPotential; Coords(i+1:i+numBorderCoord,4)= 0
			! debugging, test shape, force to all Neumann to check code bit by bit
			!Coords(i+1:i+numBorderCoord,3)=0; Coords(i+1:i+numBorderCoord,4)= 1

			! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatpan]
			Coords(i+1:i+numBorderCoord,5)= blobNum+1;
			blobPlus=blobNum+1; sizeLevel1(blobNum+1)=numBorderCoord

			Coords(i+1:i+numBorderCoord,6)=0.5*(maxR+minR); Coords(i+1:i+numBorderCoord,7)=0.5*(maxZ+minZ);
			Coords(i+1:i+numBorderCoord,8)=1

		CASE(2) 			!=2 single, real pot
			Coords(i+1,1)= temphigh;				Coords(i+1,2)= minZ
			Coords(i+2,1)= maxR;					Coords(i+2,2)= minZ
			Coords(i+3,1)= maxR;					Coords(i+3,2)= elecHeight
			Coords(i+4,1)= elecEdge+elecWidth;		Coords(i+4,2)= elecHeight !*
			Coords(i+5,1)= elecEdge+elecWidth;		Coords(i+5,2)= elecHeight 
			Coords(i+6,1)= elecEdge+elecWidth;		Coords(i+6,2)= elecHeight-elecThick
			Coords(i+7,1)= elecEdge;				Coords(i+7,2)= elecHeight-elecThick
			Coords(i+8,1)= elecEdge;				Coords(i+8,2)= elecHeight+0.15 !*
			Coords(i+9,1)= elecEdge;				Coords(i+9,2)= elecHeight+0.15
			Coords(i+10,1)= 0;						Coords(i+10,2)= elecHeight+0.15

			! sigma value; panel type. 1=Neumann
			Coords(i+1,3)= 0;						Coords(i+1,4)= 1
			Coords(i+2,3)= 0;						Coords(i+2,4)= 1
			Coords(i+3,3)= 0;						Coords(i+3,4)= 1
			Coords(i+4,3)= 0;						Coords(i+4,4)= 1
			Coords(i+5,3)= elecPotential;			Coords(i+5,4)= 0
			Coords(i+6,3)= elecPotential;			Coords(i+6,4)= 0
			Coords(i+7,3)= elecPotential;			Coords(i+7,4)= 0
			Coords(i+8,3)= elecPotential;			Coords(i+8,4)= 0
			Coords(i+9,3)= 0;						Coords(i+9,4)= 1
			Coords(i+10,3)= 0;						Coords(i+10,4)= 1

			! which blob does it belong to?
			Coords(i+1,5)= blobNum+1
			Coords(i+2,5)= blobNum+1
			Coords(i+3,5)= blobNum+1
			Coords(i+4,5)= blobNum+1; sizeLevel1(blobNum+1)=4
			Coords(i+5,5)= blobNum+2
			Coords(i+6,5)= blobNum+2
			Coords(i+7,5)= blobNum+2
			Coords(i+8,5)= blobNum+2; sizeLevel1(blobNum+2)=4 
			Coords(i+9,5)= blobNum+3
			Coords(i+10,5)=blobNum+3; sizeLevel1(blobNum+3)=2; blobPlus=blobNum+3; 

			! blob (r,z) center 
			Coords(i+1:i+4,6)=0.6 
			Coords(i+1:i+4,7)=0.5*(Coords(i+1,2)+Coords(i+4,2)); 
			Coords(i+5:i+8,6)=0.5*(Coords(i+5,1)+Coords(i+8,1)); 
			Coords(i+5:i+8,7)=0.5*(Coords(i+5,2)+Coords(i+8,2)); 
			Coords(i+9:i+10,6)=0.5*(Coords(i+9,1)+Coords(i+10,1)); 
			Coords(i+9:i+10,7)=0.5*(Coords(i+9,2)+Coords(i+10,2)); 

			! open or closed panel types
			Coords(i+1:i+numBorderCoord,8)=1

		CASE(3)	! Suvorov test case with point electrode on axis
			Coords(i+1,1)= temphigh;			Coords(i+1,2)= minZ
			Coords(i+2,1)= maxR;				Coords(i+2,2)= minZ
			Coords(i+3,1)= maxR;				Coords(i+3,2)= 0.2*(elecHeight-minZ)+minZ
			Coords(i+4,1)= maxR;				Coords(i+4,2)= 0.4*(elecHeight-minZ)+minZ
			Coords(i+5,1)= maxR;				Coords(i+5,2)= 0.6*(elecHeight-minZ)+minZ 
			Coords(i+6,1)= maxR;				Coords(i+6,2)= 0.8*(elecHeight-minZ)+minZ
			Coords(i+7,1)= maxR;				Coords(i+7,2)= elecHeight
			Coords(i+8,1)= elecWidth;			Coords(i+8,2)= elecHeight !*
			Coords(i+9,1)= elecWidth;			Coords(i+9,2)= elecHeight
			Coords(i+10,1)= 0;					Coords(i+10,2)= elecHeight

			! sigma value; panel type, 1=Neumann
			Coords(i+1,3)= 0;					Coords(i+1,4)= 1
			Coords(i+2,3)= 0;					Coords(i+2,4)= 1
			Coords(i+3,3)= 0;					Coords(i+3,4)= 1
			Coords(i+4,3)= 0;					Coords(i+4,4)= 1
			Coords(i+5,3)= 0;					Coords(i+5,4)= 1
			Coords(i+6,3)= 0;					Coords(i+6,4)= 1
			Coords(i+7,3)= 0;					Coords(i+7,4)= 1
			Coords(i+8,3)= 0;					Coords(i+8,4)= 1
			Coords(i+9,3)= elecPotential;		Coords(i+9,4)= 0
			Coords(i+10,3)= elecPotential;		Coords(i+10,4)= 0

			! which blob does it belong to?
			Coords(i+1,5)= blobNum+1
			Coords(i+2,5)= blobNum+1
			Coords(i+3,5)= blobNum+1
			Coords(i+4,5)= blobNum+1
			Coords(i+5,5)= blobNum+1
			Coords(i+6,5)= blobNum+1
			Coords(i+7,5)= blobNum+1
			Coords(i+8,5)= blobNum+1; sizeLevel1(blobNum+1)=8 
			Coords(i+9,5)= blobNum+2
			Coords(i+10,5)=blobNum+2; sizeLevel1(blobNum+2)=2; blobPlus=blobNum+2; 

			! blob (r,z) center 
			Coords(i+1:i+8,6)=0.5*(temphigh+temphigh); 
			Coords(i+1:i+8,7)=0.5*(minZ+elecHeight); 
			Coords(i+9:i+10,6)=0.5*(0+elecWidth); 
			Coords(i+9:i+10,7)=0.5*(elecHeight+elecHeight);

			! open or closed panel types
			Coords(i+1:i+numBorderCoord,8)=1

		CASE(4)	! real case with curved electrodes and higher partial "T"
			PRINT *, 'The case is currently incorrectly coded. Please remedy.'
			STOP

			Coords(i+1,1)= maxR;					Coords(i+1,2)= minZ
			Coords(i+2,1)= maxR;					Coords(i+2,2)= elecHeight
			Coords(i+3,1)= elecEdge+elecWidth;		Coords(i+3,2)= elecHeight
			Coords(i+4,1)= elecEdge+elecWidth;		Coords(i+4,2)= elecHeight-elecThick
			Coords(i+5,1)= elecEdge;				Coords(i+5,2)= elecHeight-elecThick
			Coords(i+6,1)= elecEdge;				Coords(i+6,2)= elecHeight+elecThick
			Coords(i+7,1)= elecEdge+0.5*elecWidth;	Coords(i+7,2)= elecHeight+elecThick
			Coords(i+8,1)= elecEdge+0.5*elecWidth;	Coords(i+8,2)= maxZ

			Coords(i,3)  = 0;						Coords(i,4)  = 1
			Coords(i+1,3)= 0;						Coords(i+1,4)= 1
			Coords(i+2,3)= 0;						Coords(i+2,4)= 1
			Coords(i+3,3)= elecPotential;			Coords(i+3,4)= 0
			Coords(i+4,3)= elecPotential;			Coords(i+4,4)= 0
			Coords(i+5,3)= elecPotential;			Coords(i+5,4)= 0
			Coords(i+6,3)= elecPotential;			Coords(i+6,4)= 0
			Coords(i+7,3)= 0;						Coords(i+7,4)= 1
			Coords(i+8,3)= 0;						Coords(i+8,4)= 1
		case default
			stop 'Err CaseShape'
		end select

		PRINT *, 'Done with subroutine <AddBCandAround>'
	end SUBROUTINE AddBCandAround

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
			PanelVals, whichblob, flatyes, NumBlob, MaxNumBlob, maxDistStep2)

		! incoming variables
		INTEGER, INTENT(in)							:: Numknots, flatyes, NumBlob, MaxNumBlob
		  ! how many known points, droplet or surface/boundary, blob #, # of non-boundary blobs, 
		  !   # of LS surface points. Used if EqualLsBemNumCoord=1
		REAL*8, INTENT(in)							:: maxDistStep2
		REAL*8, DIMENSION(Numknots), INTENT(in)		:: knotsx, knotsy
		REAL*8, DIMENSION(Numknots,3), INTENT(in) 	:: whichblob
		REAL*8, DIMENSION(Numknots,2), INTENT(in)	:: PanelVals

		! outsourced variables
		INTEGER, INTENT(inout)						:: NumPanels, SplinePts
		REAL*8, DIMENSION(Numknots), INTENT(inout)	:: hX

		! subroutine entirely variables
		INTEGER										:: i, j
		REAL*8, DIMENSION(Numknots)					:: totDist
		REAL*8, DIMENSION(Numknots,4)				:: Sx, Sy
		REAL*8, DIMENSION(:,:), ALLOCATABLE			:: temp8

		CONTINUE ! body of program
   		!*******
		! the final equation will look like Si[x]=ai+bi(x-xi)+ci(x-xi)^2+di(x-xi)^3

		! finds the spline coefficient fits (linear or cubic) for the knotsx/knotsy points 
		!   using x vs. y. Use only Sx
		!   NOTE: 'totDist' is meaningless in the first iteration
		CALL LinearCubic ( Sx, Sy, knotsx, knotsy, hX, Numknots, totDist, 1, flatyes, &
			NumBlob, MaxNumBlob )

		! using Sx, calculates the total distance around the shape and the corresponding
		!   number of total panels for a balance between accuracy and speed. Use totDist, SplinePts
		DO i=1,NumKnots
			totDist(i)=0; 
		ENDDO
		SplinePts=0
		CALL totDistCalc ( totDist, SplinePts, knotsx, knotsy, Numknots, flatyes, NumBlob, &
			MaxNumBlob, maxDistStep2 )

		! using totDist, finds the spline coefficient fits (linear or cubic) for the 
		!   knotsx/knotsy points using S vs. x and S vs. y. Use new Sx, Sy
		CALL LinearCubic ( Sx, Sy, knotsx, knotsy, hX, Numknots, totDist, 0, flatyes, &
			NumBlob, MaxNumBlob )

		! using Sx, Sy for the new curve, returns the evenly space points SplineCoords
		!   along the surface, Finally(!) get 'NewCoords' we desire
		IF (ALLOCATED(temp8)) DEALLOCATE(temp8)
		ALLOCATE(temp8(SplinePts,8)); 
		DO i=1,SplinePts
			temp8(i,1)=-1.; temp8(i,2)=-1.; temp8(i,3)=-1.; temp8(i,4)=-1.
			temp8(i,5)=-1.; temp8(i,6)=-1.; temp8(i,7)=-1.; temp8(i,8)=-1.
		ENDDO

		CALL SplineLengthFit (temp8, Sx, Sy, SplinePts, PanelVals, totDist, Numknots, &
			whichblob, flatyes)

		j=0
		DO i=NumPanels+1,NumPanels+SplinePts
			j=j+1
			NewCoords(i,1)=temp8(j,1); NewCoords(i,2)=temp8(j,2); NewCoords(i,3)=temp8(j,3) 
			NewCoords(i,4)=temp8(j,4); NewCoords(i,5)=temp8(j,5); NewCoords(i,6)=temp8(j,6) 
			NewCoords(i,7)=temp8(j,7); NewCoords(i,8)=temp8(j,8)
		ENDDO
		DEALLOCATE(temp8)

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
			XYFit, flatyes, BlobNum, BlobmaxNum )

		! incoming variables
		INTEGER, INTENT(in)							:: Conesize, XYFit, flatyes, &
			BlobmaxNum, BlobNum
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

		IF (BlobNum>BlobmaxNum) THEN	! force linear fit with 'CSF==0'
			CubicXYFit = 10*XYFit
		ELSE
			CubicXYFit = 10*XYFit + CubicSplineFit
		ENDIF

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
			DO i=1,ConeSize
				RHSCubic(i)=0.
				b(i)=0; c(i)=0.; d(i)=0.
			ENDDO

			a(1)=knotsy(1)
			aDiag(1)=0; bDiag(1) = 1.; cDiag(1)=0

			DO i=2,Conesize
				a(i)=1.*knotsy(i)
				aDiag(i)=1.*hX(i-1)/3.
				bDiag(i)=2.*hX(i-1)/3. + 2.*hX(i)/3.
				cDiag(i)=1.*hX(i)/3.
			ENDDO

			do i=2,Conesize
				iPlus=Nexti(i,Conesize,flatyes)
				if ( abs(hX(i)) > smallepsil .AND. abs(hX(i-1))> smallepsil .AND. a(iPlus)-a(i) > eps/10. ) then
					RHSCubic(i)=1.*(a(iPlus)-a(i))/hX(i) - 1.*(a(i)-a(i-1))/hX(i-1)
				endif
			enddo

			call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

			do i=1,Conesize
				if ( abs(hX(i)) > smallepsil ) then
					iPlus=Nexti(i,Conesize,flatyes)
					b(i) = 1.*( a(iPlus)-a(i) )/hX(i) - (2.*hX(i)*c(i)/3. + hX(i)*c(iPlus)/3.)
					d(i) = 1.*( c(iPlus)-c(i) )/(3.*hX(i))
				endif

				Sx(i,1)=a(i); Sx(i,2)=b(i)
				Sx(i,3)=c(i); Sx(i,4)=d(i)
			enddo
		case (10)				! First loop, x vs y, linear interpolation
			!-------------------------------------------------------------------
			do i=1,ConeSize
				Sx(i,1)=1.*knotsy(i)
				Sx(i,3)=0.; Sx(i,4)=0.
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

			b(1)=0.; c(1)=0.; d(1)=0.
			do i=2,Conesize
				b(i)=0.; c(i)=0.; d(i)=0.
				iPlus=Nexti(i,Conesize,flatyes)
				if ( abs(distChange(i)) > smallepsil .AND. a(iPlus)-a(i) > eps/10. ) then
					RHSCubic(i)=1.*(a(iPlus)-a(i))/distChange(i) - 1.*(a(i)-a(i-1))/distChange(i-1)
				endif
			enddo

			call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

			do i=1,Conesize
				if ( abs(distChange(i)) > smallepsil ) then
					iPlus=Nexti(i,Conesize,flatyes)
					b(i) = 1.*( a(iPlus)-a(i) )/distChange(i) - (2.*distChange(i)*c(i)/3. + distChange(i)*c(iPlus)/3.)
					d(i) = ( c(iPlus)-c(i) )/(3.*distChange(i))
				endif

				Sx(i,1)=a(i); Sx(i,2)=b(i)
				Sx(i,3)=c(i); Sx(i,4)=d(i)
			enddo

 			!******
  			! Same section, now for y vs s.
			do i=1,ConeSize
				a(i)=1.*knotsy(i)
				b(i)=0.; c(i)=0.; d(i)=0.
				RHSCubic(i)=0.
			enddo

			do i=2,Conesize
				iPlus=Nexti(i,Conesize,flatyes)
				if ( ABS(distChange(i)) > eps .AND. ABS(a(iPlus)-a(i)) > eps/10. ) then
					RHSCubic(i)=1.*(a(iPlus)-a(i))/distChange(i) - 1.*(a(i)-a(i-1))/distChange(i-1)
				endif
			enddo

			call Tridiag( c, aDiag, bDiag, cDiag, RHSCubic, ConeSize )

			do i=1,Conesize
				if ( abs(distChange(i)) > eps ) then
					iPlus=Nexti(i,Conesize,flatyes)
					b(i) = 1.*( a(iPlus)-a(i) )/distChange(i) - (2.*distChange(i)*c(i)/3. + distChange(i)*c(iPlus)/3.)
					d(i) = 1.*( c(iPlus)-c(i) )/(3.*distChange(i))
				endif

				Sy(i,1)=a(i); Sy(i,2)=b(i)
				Sy(i,3)=c(i); Sy(i,4)=d(i)
			enddo
		case (00)				! Second loop, x,y vs s, linear interpolation
			!-------------------------------------------------------------------
			DO i=1,ConeSize
				Sx(i,3)=0; Sx(i,4)=0
				Sy(i,3)=0; Sy(i,4)=0

				if (i /= 1) then
					distChange(i) = 1.*( totDist(i)-totDist(i-1) )
				else
					distChange(i) = 1.*totDist(i)
				endif

				Sx(i,1)=1.*knotsx(i)
				Sy(i,1)=1.*knotsy(i)
			ENDDO

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
		INTEGER, INTENT(in)							:: SplinePts, ConeSize, flatyes2
		  ! 'SplinePts'=# of new coords. 'ConeSize'=# of old coords. 'flatyes2' open or closed
		REAL*8, DIMENSION(ConeSize,5), INTENT(in)	:: Sx, Sy
		REAL*8, DIMENSION(ConeSize,3),INTENT(in)	:: whichblob1
		REAL*8, DIMENSION(ConeSize,2), INTENT(in)	:: PanelVals
		REAL*8, DIMENSION(ConeSize), INTENT(in)		:: totDist

		! outsourced variables
		REAL*8,DIMENSION(SplinePts,8),INTENT(inout) :: SplineCoords

		! subroutine entirely variables
		INTEGER										:: tempJ, i, j
		REAL*8										:: temptotDist, tempfracDist

		REAL*8, DIMENSION(SplinePts)				:: SplineLength

		CONTINUE
		tempJ=1
		DO i=1,SplinePts
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
		ENDDO
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
			NumBlob3, MaxNumBlob3, maxDistStep2)
		! incoming variables
		INTEGER, INTENT(in)							:: ConeSize1, flatyes3, NumBlob3, MaxNumBlob3
		REAL*8, INTENT(in)							:: maxDistStep2
		REAL*8, DIMENSION(ConeSize1), INTENT(in)	:: knotsx1, knotsy1

		! outsourced variables
		INTEGER, INTENT(out)						:: SplinePts1
		REAL*8, DIMENSION(ConeSize1), INTENT(inout)	:: totDist1

		! subroutine entirely variables
		INTEGER										:: i, iPlus
		REAL*8										:: tempArcDist
		REAL*8, DIMENSION(ConeSize1)				:: xPanLen, yPanLen, dydx

		CONTINUE		
		tempArcDist=0.0

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
		!PRINT *, 'totDist1(ConeSize1)=',totDist1(ConeSize1)
		IF ( NumBlob3 > MaxNumBlob3 ) THEN	! have fewer points around the boundary
			SplinePts1 = INT(MAX( StepMult*FLOOR(totDist1(ConeSize1)/maxDistStep2)+1.0, 1.0*MinNumPtsBlob ))
		ELSE
			SplinePts1 = INT(MAX( FLOOR(totDist1(ConeSize1)/maxDistStep2)+1.0, 1.0*MinNumPtsBlob ))
		ENDIF
	END SUBROUTINE totDistCalc

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
! DESCRIPTION: iteratively solves the matrix for x. 
!
! INPUT: influence matrix A, RHS b, initial matrix guess x (can be zero)
!
! OUTPUT: new x
!
! CALLING PROGRAM: MAT_Greens

	subroutine gmresIter(PanStr, a, b, num)

		! incoming variables
		integer, intent(in)				:: num
		real*8, dimension(num*num), intent(in) 	:: a
		real*8, dimension(num), intent(in) 		:: b

		! outgoing variable array
		real*8, dimension(num), intent(inout)	:: PanStr

		! program variables
		integer							:: go, x0_zero ! formerly logical
		integer							:: i, j, n, n1, RestartLimit, countRestart, &
			countRestartReset, finalsigunit, intersigunit
		real*8							:: temp, errlim, bnorm, tempYn, ynorm, ynormOld, normerr
		real*8, dimension(:), allocatable	:: v, y, r, b1, qfact, givens, q, h, rfact, &
			finalans, tempPanStr, tempR

		!  real*8, dimension(num,num)		:: q	! Q satisfying the equation AQ(n)=Q(n+1)H(n)
		!  real*8, dimension(num+1,num)     :: h	! Hessenberg matrix formed by orthogonal transformation of A
		!  real*8, dimension(num)			:: v	! Vector used for computation of Aq(n) and other quantities
		!  real*8, dimension(num)			:: y	! Vector satisfying the least squares problem ||(Hy-||b||e1)||
		!  real*8, dimension(num)			:: r	! Residual vector, which is solved for given an initial guess
		!  real*8, dimension(num)			:: b1	! Vector equaling b-Ax(0), so that Ar=b1 can be solved

		!  real*8, dimension(num+1)			:: qfact  ! The first row of the Q matrix of H=QR
		!  real*8, dimension(num,num)		:: rfact  ! The R matrix of H's QR Decomposition
		!  real*8, dimension(2*num)			:: givens ! The Givens Rotation coefficients for each iteration

		CONTINUE
		finalsigunit=10; intersigunit=11

		allocate(q(num*num))
		allocate(h((num+1)*num))
		allocate(v(num))
		allocate(r(num))
		allocate(b1(num))
		allocate(qfact(num+1))
		allocate(rfact(num*num))
		allocate(givens(2*num))
		ALLOCATE(finalans(num))

		j=num*num
		DO i=1,j
			q(i)=0.; h(i)=0.; rfact(i)=0.

			IF (i<=num) THEN
				v(i)=0.; r(i)=0.; b1(i)=0.
				qfact(i)=0. 
				givens(i)=0.; givens(num+i)=0.
				finalans(i)=-1.
			ENDIF
		ENDDO
		qfact(num+1)=0.
		n=(num+1)*num
		DO i=j+1,n
			h(i)=0.
		ENDDO

		IF (1==ShowGMResConverge) THEN
			OPEN(unit=finalsigunit, file='Sigma.dat', action='read')
			DO i=1,num
				READ(unit=finalsigunit, fmt='(e15.6)', END=798) finalans(i)
			ENDDO
			798 CLOSE (unit=finalsigunit, status='keep')

			OPEN(unit=intersigunit, file='GMconverge2.dat', status='replace')
		ENDIF

  		!RestartLimit=min(floor(num/10.)+1,35); countRestart=0; countRestartReset=0
		RestartLimit=num+5; countRestart=0; countRestartReset=0
		allocate(y(RestartLimit+1))
		DO i=1,RestartLimit+1
			y(i)=0.
		ENDDO

  		! Set the b1-vector using the initial guess of x, b1=b-Ax(0)
		x0_zero=1				! assume initial vector all zeros
		do i=1,num
			if (PanStr(i) /= 0) then
				x0_zero=0
				exit
			endif
		enddo

		errlim=1e-7; go=1; tempYn=1e6; ynorm=0.
		RestartLoop: do while (.TRUE.)
			n=1; n1=2

			q(:)=0; y(:)=0; v(:)=0
			qfact(:)=0; rfact(:)=0; givens(:)=0

  			! Set the initial values of qfact
			qfact(1)=1.0; qfact(2)=1.0

			if (0==x0_zero) then	! vector isn't all zeros
				do i=1,num
					b1(i)=b(i)
					do j=1,num
						b1(i) = b1(i) - a((i-1)*num+j)*PanStr(j)
					enddo
				enddo
			elseif (1==x0_zero) then
				b1=b
			else 
				print *, 'GMRes, treeyes error'
				stop
			endif

  			! Compute the norm of the b1 vector
			bnorm=0.0;
			do i=1,num
				bnorm = bnorm + b1(i)*b1(i)
			enddo
			bnorm = sqrt(bnorm)

  			! Compute q(i,1)=b1/||b1||
			do i=1,num
				q(0*num+i)=b1(i)/bnorm
			enddo

			GoWhile: do while (1==go .AND. n<RestartLimit)  
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
					temp = temp + v(i)*v(i)
				enddo
				h((n1-1)*num+n)=sqrt(temp);

				if (n1 < num) then
      				! Compute q(n+1) = v/h(n+1,n) = v/||v||
					do i=1,num
						q(n*num+i)=v(i)/h(n*num+n)
					enddo
				endif

    			! Do the QR Factorization using a Givens Rotation
				rfact(n)=h(n)

    			! Apply the previous rotations to the new column of rfact
				do i=1,n-1
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

				ynormOld=ynorm
				ynorm=abs(y(n)/bnorm)

				IF (1==ShowGMResConverge) THEN
					ALLOCATE(tempR(num)); ALLOCATE(tempPanStr(num))

					! Set x(i) = q(i)*y
					DO i=1,num
						tempR(i)=0.
						tempPanStr(i)=0
					ENDDO

					do j=1,n1
						do i=1,num
							tempR(i) = tempR(i) + q((j-1)*num+i)*y(j)
						enddo
					enddo

					normerr=0
					do i=1,num
						tempPanStr(i) = tempPanStr(i) + tempR(i)
						normerr=normerr + (tempPanStr(i)-finalans(i))/finalans(i)
					enddo
					normerr=ABS(normerr/num)

					WRITE(unit=intersigunit, fmt='(i8, e13.4)') n, normerr
					DEALLOCATE(tempR); DEALLOCATE(tempPanStr)
				ELSE

				ENDIF

    			! Calculate x=x+Q*y when converged
				if (n1==num .OR. (abs(ynorm) < errlim .AND. &
					ynorm < ynormOld .AND. abs(y(n)) >1e-20 .AND. n>num/20 )) then
					! Set x(i) = q(i)*y
					DO i=1,num
						r(i)=0.; PanStr(i)=0.
					ENDDO
					do j=1,n1
						do i=1,num
							r(i) = r(i) + q((j-1)*num+i)*y(j)
						enddo
					enddo

					normerr=0
					do i=1,num
						PanStr(i) = PanStr(i) + r(i)
					enddo

					go=0
				endif
				n=n+1; n1=n+1
			enddo GoWhile		! while go

			IF (1==ShowGMResConverge) THEN
				CLOSE (unit=intersigunit, status='keep')
			ENDIF

  			! exited loop. Decide to exit or restart
			if (0==go) then
				exit
			else
				countRestart=countRestart+1
				print *, countRestart
				x0_zero=0
				DO i=1,num
					r(i)=0.
				ENDDO

				if (abs(y(n-1))>tempYn .AND. RestartLimit<=num) then
					countRestartReset=countRestartReset+countRestart
					countRestart=0
					tempYn=1e6
					deallocate(y)
					RestartLimit=RestartLimit*2
					allocate(y(RestartLimit+1)); 
					DO i=1,RestartLimit+1
						y(i)=0.
					ENDDO
				else 
					tempYn=abs(y(n-1))
				endif

    			! Set x(i) = q(i)*y
				do j=1,n1
					do i=1,num
						r(i) = r(i) + q((j-1)*num+i)*y(j)
					enddo
				enddo

				do i=1,num
					PanStr(i) = PanStr(i)+r(i)
				enddo
			endif
		enddo RestartLoop	! while Restart

		write(*, fmt=10), '   Done after ', n-1, ' GM_Res steps.'
		!  write(*, fmt=10), '     with ', countRestartReset, ' restarts.'

		deallocate(q); deallocate(h); deallocate(v); deallocate(y)
		deallocate(r); deallocate(b1); deallocate(qfact); deallocate(rfact)
		deallocate(givens); DEALLOCATE(finalans)

		goto 999		! actual program end, just lists formats and errors below

  		! ***********
		! FORMAT listings
		10  FORMAT(a, i4, a)

	999 end SUBROUTINE gmresIter

! ***************
! SUBROUTINE CalcPDF
! ***************
! DESCRIPTION: Once the first blob has formed, calculate snapoff PDF
!
! INPUT: <>
!
! OUTPUT: Droplet charge around surface
!
! CALLING PROGRAM: MAT_Greens:BFC

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
			temptotCharge=DropEn(i)*DropLen(i)
			totCharge=totCharge+temptotCharge
		ENDDO
		totCharge=ABS(epsil0*totCharge)
	END SUBROUTINE CalcPDF

! ***************
! SUBROUTINE POLYGON_AREA_2D_2
! ***************
! DESCRIPTION: Computes the area of a polygon in 2D. The area is the sum of the 
! areas of the triangles formed by node N with consecutive pairs of nodes.
!
! INPUT: integer: N, the number of vertices of the polygon. real:V(N,2), the vertices.
!
! OUTPUT: areaPoly, the absolute area of the polygon.
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
! CALLING PROGRAM: MAT_Greens:BFC

	SUBROUTINE polygon_area_2d_2 ( areaPoly, v, n )
		! incoming
		INTEGER, INTENT(in) 	:: n
		REAL*8, DIMENSION(n,2), INTENT(in) :: v

		! outgoing
		REAL*8, INTENT(inout) 	:: areaPoly

		! internal
		INTEGER					:: i, j
		REAL*8					:: areat
		REAL*8, DIMENSION(3,2) 	:: t

		CONTINUE 
		areaPoly = 0.0

		DO i = 1, n - 2
			!t = RESHAPE ( (/ v(1:2,i), v(1:2,i+1), v(1:2,n) /), (/ 2, 3 /) )
			DO j=1,2
				t(1,j) = v(i,j)
				t(2,j) = v(i+1,j)
				t(3,j) = v(n,j)
			ENDDO

			areat = 0.5 * ABS ( t(1,1) * ( t(2,2) - t(3,2) ) &
				+ t(2,1) * ( t(3,2) - t(1,2) ) + t(3,1) * ( t(1,2) - t(2,2) ) )

			areaPoly = areaPoly + areat
		ENDDO

		RETURN
	END SUBROUTINE polygon_area_2d_2

! ***************
! SUBROUTINE CalcDropletPot
! ***************
! DESCRIPTION: Computes the correct potential boundary condition for detached droplets
!
! INPUT: Coords, droplet charges and centers
!
! OUTPUT: updated dirichlet potentials
!
! CALLING PROGRAM: MAT_Greens:BFC
!
! ALGORITHM: 1) linearly superimpose droplet charges as a potential on other droplets
!	2) load and add background potential to droplet value

	SUBROUTINE CalcDropletPot (tempCoord, blobq, PotErr, LastPot, &
			numBlobErr, blob_rc, blob_zc, sizeblob, numneedle, numblob, &
			numblobplusBc, numpan, newdropNum, wts, RecalcPot, redoneGP, &
			numnewcoords, sizenewcoords, cumsizeGP, whichBlob, Areab, &
			DropPotSet, z2)

		! incoming
		LOGICAL, INTENT(in)		:: RecalcPot
		INTEGER, INTENT(in) 	:: numneedle, numblob, numblobplusBc, &
			numpan, newdropNum, numnewcoords, numBlobErr, z2
		INTEGER, DIMENSION(numblob), INTENT(in)			:: whichBlob !which previous blob #
		INTEGER, DIMENSION(numblobplusBc), INTENT(in) 	:: sizeblob, sizenewcoords, cumsizeGP
			! #Coords/blob, #NewCoords/blob, cumulative starting panel # for each blob
		REAL*8, DIMENSION(numblob), INTENT(in)  		:: Areab
		REAL*8, DIMENSION(numblobplusBc), INTENT(in) 	:: blob_rc, blob_zc, blobq
		REAL*8, DIMENSION(npoints), INTENT(in) 			:: wts

		! outgoing
		REAL*8, INTENT(inout)	:: PotErr, LastPot
		REAL*8, DIMENSION(numpan*3), INTENT(inout)		:: tempCoord
			! Coords(:,3:5), "Pt val", "Pt type", "Blob#"
		REAL*8, DIMENSION(numblob, MaxNumInpotIter), INTENT(inout) 	:: DropPotSet
		TYPE(PANEL), DIMENSION(numnewcoords), INTENT(inout)			:: redoneGP

		! internal
		CHARACTER				:: aC
		CHARACTER (LEN=12)		:: PotCase
		INTEGER, PARAMETER		:: Potunit=33, Coordsunit2=34, tempPanunit=35
		REAL*8, PARAMETER		:: k=1.0/(4.0*pi*epsil0)	! coulomb's const

		INTEGER					:: i, j, kk, m, n, filestat, nearM, nearM2, &
			nearN, nearN2, sizeRedo, tempBlobNum, startblob, pan
		INTEGER, DIMENSION(:), ALLOCATABLE 	:: InvwhichBlob
		REAL*8					:: tempPot, temp, ri, dist, xp, yp, distm, &
			lowz, highz, minPot, maxPot
		REAL*8, DIMENSION(4)	:: xdiv, ydiv, phidiv
		REAL*8, DIMENSION(:), ALLOCATABLE	:: meshR, meshZ
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: GridPot
		TYPE(PANEL),DIMENSION(numnewcoords) :: tempGP

		CONTINUE 
		temp=0.; minPot=MIN(fluidPot,elecPotential)
		maxPot=MAX(fluidPot,elecPotential)

		ALLOCATE(meshR(PotMeshSize))
		ALLOCATE(meshZ(PotMeshSize))
		ALLOCATE(GridPot(PotMeshSize,PotMeshSize))
		ALLOCATE(InvwhichBlob(numblob))

		DO i = 1, PotMeshSize
			meshR(i)=-1.; meshZ(i)=-1.
			DO j = 1, PotMeshSize
				GridPot(i,j)=-1.
			ENDDO
		ENDDO
		DO i=1,numblob
			InvwhichBlob(whichBlob(i))=i
		ENDDO

		ShowLoopDetails: IF (1==ShowInLoopCoords) THEN
			PRINT *, '  -Writing "loopCoords.dat" '
			OPEN(unit=Coordsunit2, file='loopCoords.dat', status='replace')
			WRITE(unit=Coordsunit2, fmt='(a)') 'VARIABLES="Pt val"'
			WRITE(unit=Coordsunit2, fmt='(a)') '"Pt type"'
			WRITE(unit=Coordsunit2, fmt='(a)') '"Blob#"'
			WRITE(unit=Coordsunit2, fmt='(a)') '"Pan #"'
			DO i=1,numpan
				WRITE(unit=Coordsunit2, fmt='(e12.4, 3i4)') tempCoord(i), &
					INT(tempCoord(i+numpan)), INT(tempCoord(i+2*numpan)), i
			ENDDO
			CLOSE (unit=Coordsunit2, status='keep')

			IF (numnewcoords>0) THEN
				PRINT *, '  -Writing "loopPanDetail.dat" '
				OPEN(unit=tempPanunit, file='loopPanDetail.dat', status='replace')
				WRITE(unit=tempPanunit, fmt='(a)') 'VARIABLES="R [cm]"'
				WRITE(unit=tempPanunit, fmt='(a)') '"Z [cm]"'
				WRITE(unit=tempPanunit, fmt='(a)') '"Pt val"'
				WRITE(unit=tempPanunit, fmt='(a)') '"Pt type"'
				WRITE(unit=tempPanunit, fmt='(a)') '"Blob Num"'
				WRITE(unit=tempPanunit, fmt='(a)') '"Rnorm"'
				WRITE(unit=tempPanunit, fmt='(a)') '"Znorm"'
				WRITE(unit=tempPanunit, fmt='(a)') '"Pan #"'
				WRITE(unit=tempPanunit, fmt='(a)') 'ZONE T="main"'

				DO m = 1,numblobplusBc
					kk=sizenewcoords(m)
					tempBlobNum=m

					DO j=1,kk
						pan=cumsizeGP(tempBlobNum)+j-1
						WRITE(unit=tempPanunit, fmt='(3e12.4, 2i3, 2e12.4, i5)') &
							redoneGP(pan)%midr, redoneGP(pan)%midz, &
							redoneGP(pan)%midPot, INT(redoneGP(pan)%type1), &
							tempBlobNum, redoneGP(pan)%rnrm, redoneGP(pan)%znrm, pan
						!DO n=1,npoints
						!WRITE(unit=tempPanunit, fmt='(3e12.4, 2i2, 2e12.4, i5)') &
						!	redoneGP(pan)%rpoints(n), redoneGP(pan)%zpoints(n), &
						!	redoneGP(pan)%str, INT(redoneGP(pan)%type1), &
						!	tempBlobNum, redoneGP(pan)%rnrm, redoneGP(pan)%znrm, pan
						!ENDDO
					ENDDO
				ENDDO
				CLOSE (unit=tempPanunit, status='keep')
			ENDIF
		ENDIF ShowLoopDetails
		DO i=1,numnewcoords
			!.....
			CALL ZeroorCopyPanel(tempGP(i), 0, redoneGP(i))
			!.....
		ENDDO

		IF ( .NOT. RecalcPot) THEN
			! load background electric field, computed beforehand
			!.....
			CALL GetPotGrid(temp, PotCase, temp, temp, temp, 1)
			!.....
			OPEN(UNIT=Potunit, file=PotCase, status='old', IOSTAT=filestat)
			IF (0==filestat) THEN	! base background file exists
				DO i=1,6			! skip header
					READ(UNIT=Potunit, FMT='(a)', END=1257) aC
				ENDDO

				1257 DO i=1,PotMeshSize
					DO j=1,PotMeshSize
						READ(UNIT=Potunit, FMT='(3e14.6)', END=1262) meshR(i), meshZ(j), GridPot(i,j)
					ENDDO
				ENDDO
			ELSE	
				PRINT *, '  Panels:CalcDropletPot missing background electric field.'
			ENDIF
			1262 CLOSE(unit=Potunit, status='keep')

			PRINT *, '    Computing efield from base setup on drop center'	! used later, but right spot
		ENDIF

		! *CHANGE Coord boundary types and potentials!
		! #1,2) use these + detached blob charges to calc potential on detached blobs
		tempPot=0.
		ChangeBlob: DO i=1,numblob	! base blob comparing from
			ChangeBlobIf: IF (i /= numneedle) THEN
				Recalcif: IF (RecalcPot) THEN	! directly calculate the potential on the surface
					! only using the NewCoords. Requires it to be iteratively called and already
					! defined, with both position and 'sigma'. Otherwise, this section has no meaning.

					IF (1==RecastDist ) THEN ! changes matrix to a meters-based system instead of cm.
						CALL ForceToMeters(redoneGP, numnewcoords, 1);
					ENDIF

					! set the %str of this blob to zero to effectively remove it
					kk=sizenewcoords(i)
					DO m=0,kk-1
						tempGP(cumsizeGP(i)+m)%str=0.
					ENDDO
					sizeRedo=kk

					!.....	calcing pot 'temp' at blob center, minus the blob
					temp=0.
					CALL calcpotwrap( temp, blob_zc(i), blob_rc(i), numnewCoords, &
						wts, tempGP )

					! return to previous %str
					DO m=0,kk-1
						tempGP(cumsizeGP(i)+m)%str=redoneGP(cumsizeGP(i)+m)%str
					ENDDO

					IF (1==RecastDist ) THEN ! returns matrix to a cm-based system instead of m.
						CALL ForceToMeters(redoneGP, numnewcoords, 0);
					ENDIF

					tempPot=temp
				ELSE	! (.NOT. Recalcif): calc potential from background efield. 
						!    NewCoords not yet created
					IF (1==RecastDist ) THEN ! changes matrix to a meters-based system instead of cm.
						distm=100.
					ELSE
						distm=1.
					ENDIF
					tempPot=0.
					DO j=1,numblob	! all other detached blobs' affects
						IF (j /= i .AND. j /= numneedle) THEN
							ri=SQRT( (blob_rc(i)-blob_rc(j))*(blob_rc(i)-blob_rc(j)) + &
								(blob_zc(i)-blob_zc(j))*(blob_zc(i)-blob_zc(j)) )
							tempPot = tempPot + blobq(j)/(distm*ri)
						ENDIF
					ENDDO
					tempPot = tempPot*k

					PotFileIf: IF (0==filestat) THEN	! base background file exists
						!------------------
						! add in background electric charges
						!------------------
						! find nearest 4 predefined points
						temp=1e6
						DO m=1,PotMeshSize	! closest R
							dist=(blob_rc(i)-meshR(m))*(blob_rc(i)-meshR(m))
							IF (dist<=temp) THEN
								temp=dist
								nearM=m

								IF (blob_rc(i)-meshR(m)>0) THEN
									nearM2=m+1
								ELSE
									nearM2=m-1
								ENDIF
							ENDIF
						ENDDO
						IF (1==nearM) THEN
							nearM2=nearM+1
						ELSEIF(PotMeshSize==nearM) THEN
							nearM2=nearM-1
						ENDIF
						temp=1e6
						DO n=1,PotMeshSize	! closest Z
							dist=(blob_zc(i)-meshZ(n))*(blob_zc(i)-meshZ(n))
							IF (dist<=temp) THEN
								temp=dist
								nearN=n

								IF (blob_zc(i)-meshR(n)>0) THEN
									nearN2=n-1
								ELSE
									nearN2=n+1
								ENDIF
							ENDIF
						ENDDO
						IF (1==nearN) THEN
							nearN2=nearN+1
						ELSEIF(PotMeshSize==nearN) THEN
							nearN2=nearN-1
						ENDIF

						! add in linearly interpolated background potential
						!.....
						xdiv(1)=meshR(nearM);	ydiv(1)=meshZ(nearN);	phidiv(1)=GridPot(nearM, nearN)
						xdiv(2)=meshR(nearM2);	ydiv(2)=meshZ(nearN);	phidiv(2)=GridPot(nearM2, nearN)
						xdiv(3)=meshR(nearM2);	ydiv(3)=meshZ(nearN2);	phidiv(3)=GridPot(nearM2, nearN2)
						xdiv(4)=meshR(nearM);	ydiv(4)=meshZ(nearN2);	phidiv(4)=GridPot(nearM, nearN2)
						xp = blob_rc(i);		yp = blob_zc(i)

						!.....
						temp=0.
						CALL LinInterpolate(xdiv, ydiv, phidiv, xp, yp, temp, 4)
						!.....
						IF (temp<minPot .OR. temp>maxPot) THEN 
							temp= (maxPot-minPot)*ABS((blob_zc(i)-lowz)/ &
								(highz-lowz)) + minPot
						ENDIF

						tempPot=tempPot+temp
						IF ( tempPot<MIN(fluidPot,elecPotential)-1 .OR. &
							tempPot>MAX(fluidPot,elecPotential)+1 ) THEN	! uh-oh.
							PRINT *, '** Panels:CalcDropletPot, pot background calc err. fluidPot,elecPotential, newPot=', &
								fluidPot,elecPotential,tempPot
							STOP
						ENDIF
					ELSE	! no pot. file exists. Make rudimentary first guess	
						lowz=1e6; highz=0
						DO kk=1,numblobplusBc
							IF (blob_zc(kk)<lowz) 	lowz=blob_zc(kk)
							IF (blob_zc(kk)>highz) 	highz=blob_zc(kk)
						ENDDO
						temp= (maxPot-minPot)*ABS((blob_zc(i)-lowz)/ &
							(highz-lowz)) + minPot
						!tempPot=tempPot + 0.58177*(elecPotential-fluidPot)*EXP(LOG(blob_zc(i)-lowz+eps)/ &
						!	LOG(highz-lowz)) + fluidPot
					ENDIF PotFileIf
					sizeRedo=sizeblob(InvwhichBlob(i))
				ENDIF Recalcif

				DropPotSet(i,z2)=tempPot
				WRITE(*, fmt='(a, i2, a, e9.3, a, e10.4, a)') '   BC droplet(', i, '): Charge=', &
					blobq(i), ' [C] & Potential=', tempPot, ' [V]'
				WRITE(*,fmt='(a,3e11.3,i7)') '       rc, zc, area, numCoord=',blob_rc(i), &
					blob_zc(i), Areab(i), sizeRedo

				IF (Areab(i)>1e3) THEN
					PRINT *, 'ISNAN break'
					STOP
				ENDIF

				IF ( tempPot<MIN(fluidPot,elecPotential)-1 .OR. &
					tempPot>MAX(fluidPot,elecPotential)+1 ) THEN	! uh-oh. Force to position interpolation
					!tempPot=0.58177*(elecPotential-fluidPot)*EXP(LOG(blob_zc(i)-lowz+eps)/ &
					!	LOG(highz-lowz)) + fluidPot
					PRINT *, '** Panels:CalcDropletPot, pot calc err. fluidPot,elecPotential, newPot=', &
						fluidPot,elecPotential,tempPot
					STOP
				ENDIF

				IF (i==numBlobErr) THEN
					IF (0==lastPot) THEN
						PotErr=100.
					ELSE
						PotErr=ABS((tempPot-lastPot)/lastPot)
					ENDIF

					WRITE(*, fmt='(a, e8.2, a)') '     Potential iterate err = ', 100.*PotErr, '%.'
					LastPot=tempPot
				ENDIF

				! change from Neumann to Dirichlet boundary conditions; update fixed pot.
				n=sizeblob(i); startblob=1
				! get to starting number
				DO WHILE(INT(tempCoord(startblob+2*numpan))/=i)
					IF(INT(tempCoord(startblob+2*numpan))==i) EXIT
					startblob=startblob+1
				ENDDO

				ChangeCoords: DO m=startblob,n+startblob-1
					tempCoord(m) = tempPot
					tempCoord(m+numpan) = 0.
				ENDDO ChangeCoords
			ENDIF ChangeBlobIf
			pan=pan+sizenewcoords(i)
		ENDDO ChangeBlob

		PotCase=aC; i=newdropNum		! doesn't do anything, but suppresses warning
		DEALLOCATE(meshR); DEALLOCATE(meshZ); DEALLOCATE(GridPot)
		DEALLOCATE(InvwhichBlob)
	END SUBROUTINE CalcDropletPot

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
	END FUNCTION Nexti

! ***************
! FUNCTION SolveCubicEq
! ***************
! DESCRIPTION: solves for the most positive root of (x^3) + a1(x^2) + a2(x) + a3=0
!
! INPUT: a+b(delX)+c(delX^2)+d(delX^3)=y
!
! OUTPUT: x
!
! CALLING PROGRAM: SplineLengthFit
!
! Primary author: Anton VanderWyst, 
!	University of Michigan 
!
! Version history:
!   18 Jan 2005 begun

	FUNCTION SolveCubicEq ( a,b,c,d,y )
		! incoming variables
		REAL*8, intent(in)	:: a,b,c,d,y

		! outgoing function variable
		REAL*8				:: SolveCubicEq

		! internal variables	
		REAL*8				:: a1,a2,a3, Q, R, S, T

		continue

		! change input into form (x^3) + a1(x^2) + a2(x) + a3=0
		a1=1.*c/d; a2=1.*b/d; a3=1.*(a-y)/d
		Q = (3.*a2-a1**2)/9.
		R = (9.*a1*a2-27.*a3-2.*a1**3)/54.
		S = ( R+ ( Q**3 + R**2 )**0.5   )**(1./3)
		T = ( R- ( Q**3 + R**2 )**0.5   )**(1./3)

		SolveCubicEq = S + T -a1/3.
	END FUNCTION SolveCubicEq

! ***************
! SUBROUTINE GetPotGrid
! ***************
! DESCRIPTION: returns mesh points for a potential grid, depending on which shape is chosen in 
!   global:CaseShape and defined in Panels:AddBCandAround
!
! INPUT: shape type, max & min values. Global values of 'CalcPotBack', 'CaseShape' & 'PotMeshSize'
!
! OUTPUT: grid point value
	SUBROUTINE GetPotGrid(gridPt, potcase2, offset, minVal, maxVal, index)
		! incoming
		INTEGER, INTENT(in)		:: index
		REAL*8, INTENT(in)		:: offset, minVal, maxVal

		! outgoing		
		CHARACTER(len=12), INTENT(out) :: potcase2
		REAL*8, INTENT(inout)	:: gridPt

		CONTINUE

		IF (CalcPotBack/=1) THEN
			gridPt = (offset*(maxVal-minVal)+minVal) + ((2*offset*(minVal-maxVal)+maxVal-minVal)* &
				(index-1)/(PotMeshSize-1))	
		ELSE
			gridPt = minVal + (maxVal-minVal)*(index-1)/(PotMeshSize-1)
		ENDIF

		SELECT CASE(CaseShape)
		CASE(1)			! single box, simple pot
			potcase2='PotCase1.dat'
		CASE(2)			! ripple box, real pot
			potcase2='PotCase2.dat'				
		CASE(3)			! Suvorov test case with electrode on axis
			potcase2='PotCase3.dat'
		CASE(4)			! 'T' top shape				
			potcase2='PotCase4.dat'
		CASE DEFAULT
			PRINT *, 'Err in MAT_Greens:CaseShape'
		END SELECT

	END SUBROUTINE GetPotGrid

! ***************
! SUBROUTINE LinInterpolate
! ***************
! DESCRIPTION: for any straight-sided shape, returns linear interpolation of potential
!   for the interior or edge point
!
! INPUT: (x,y,phi) for n vertices + (xp,yp) of desired point. If using 'WeightMethod'== 1,
!   points need to be entered consecutively
!
! OUTPUT: phi at desired point(p)
!	
! AUTHOR: Anton VanderWyst
! BASE ALGORITHM CONCEPTS: John Yim(1), Yongjun Choi(2) and Leo Scalabrin(4)

	SUBROUTINE LinInterpolate(x, y, phi, xp, yp, phip, n)
		! incoming
		INTEGER, INTENT(in)	:: n
		REAL*8, INTENT(in)	:: xp, yp
		REAL*8, DIMENSION(n), INTENT(in) :: x, y, phi

		! outgoing		
		REAL*8,INTENT(inout):: phip

		! internal
		INTEGER, PARAMETER	:: WeightMethod=1
			! =1 for bisected inverse area-weighted, =2 for distance-weighted, 
			!= 3 for dist^2 weighted, =4 for inverse distance-weighted
			! recommended to use only (=1) or (=4) as they are more accurate
		REAL*8, PARAMETER	:: dmin=1e-10

		INTEGER				:: i, j, jj
		REAL*8				:: dtot, ddtot, da, db
		REAL*8,DIMENSION(4) :: d, w, xbi, ybi
		REAL*8,DIMENSION(3,2):: t, t2

		CONTINUE

		SELECT CASE(WeightMethod)
		CASE(1)
			! find bisected {x,y} points
			DO i=0,n-1
				xbi(i+1)=0.5*( x(MOD(i,n)+1) + x(MOD(i+1,n)+1) )
				ybi(i+1)=0.5*( y(MOD(i,n)+1) + y(MOD(i+1,n)+1) )
			ENDDO

			! find area of each triangle part.
			dtot=0.
			DO i = 1, n
				j=MOD(n+i-2,n)+1; jj=MOD(n+i-1,n)+1  ! note: jj==i

				! first triangle from bisect(i-1) to vertex(i)
				t(1,1) = xbi(j);	t(1,2) = ybi(j)
				t(2,1) = x(jj);		t(2,2) = y(jj)
				t(3,1) = xp;		t(3,2) = yp

				! second triangle from vertex(i) to bisect(i)
				t2(1,1) = x(jj);	t2(1,2) = y(jj)
				t2(2,1) = xbi(jj);	t2(2,2) = ybi(jj)
				t2(3,1) = xp;		t2(3,2) = yp

				da = 0.5 * ABS( t(1,1) * ( t(2,2) - t(3,2) ) &
					+ t(2,1) * ( t(3,2) - t(1,2) ) + t(3,1) * ( t(1,2) - t(2,2) ) )

				db = 0.5 * ABS( t2(1,1) * ( t2(2,2) - t2(3,2) ) &
					+ t2(2,1) * ( t2(3,2) - t2(1,2) ) + t2(3,1) * ( t2(1,2) - t2(2,2) ) )

				d(i)=da+db

				IF (d(i)<dmin) THEN
					d(i)=dmin
				ENDIF
				dtot=dtot+d(i)
			ENDDO

			! weights are relative fraction area away from each edge
			ddtot=0.
			DO i=1,n
				w(i)=1./d(i)
				ddtot=ddtot+w(i)
			ENDDO

			DO i=1,n
				w(i)=w(i)/ddtot
			ENDDO

			! new potential is area-weighted sum of edge potentials
			phip=0.
			DO i=1,n
				phip=phip+w(i)*phi(i)
			ENDDO
		CASE(2)
			! get distance from {xp,yp} to each edge
			dtot=0
			DO i = 1, n
				d(i) = SQRT((xp-x(i))*(xp-x(i)) + (yp-y(i))*(yp-y(i)) )
				dtot=dtot+d(i)
			ENDDO

			! weights are relative fraction away from each edge
			DO i=1,n
				w(i)=dtot-d(i)
			ENDDO

			! new potential is distance-weighted sum of edge potentials
			phip=0.
			DO i=1,n
				phip=phip+w(i)*phi(i)
			ENDDO
			phip=phip/((n-1)*dtot)
			!phip = (w(1)*phi1 + w(2)*phi2 + w(3)*phi3 + w(4)*phi4)/(w(1)+w(2)+w(3)+w(4))
		CASE(3)
			! get distance from {xp,yp} to each edge
			dtot=0
			DO i = 1, n
				d(i) = SQRT((xp-x(i))*(xp-x(i)) + (yp-y(i))*(yp-y(i)) )
				dtot=dtot+d(i)
			ENDDO

			! weights are relative squared fraction away from each edge
			ddtot=0
			DO i=1,n
				w(i)=1.*(dtot-d(i))/dtot
				w(i)=w(i)*w(i)
				ddtot=ddtot+w(i)
			ENDDO

			! renorm the weights from the squared setup
			DO i=1,n
				w(i)=w(i)/ddtot
			ENDDO

			! new potential is distance-weighted sum of edge potentials
			phip=0.
			DO i=1,n
				phip=phip+w(i)*phi(i)
			ENDDO
		CASE(4)
			! get distance from {xp,yp} to each edge
			dtot=0
			DO i = 1, n
				d(i) = SQRT((xp-x(i))*(xp-x(i)) + (yp-y(i))*(yp-y(i)) )
				IF (d(i)<=dmin) THEN
					d(i)=dmin
				ENDIF
				dtot = dtot + 1./d(i)
			ENDDO

			! weights are relative inverse fraction away from each edge
			DO i=1,n
				w(i)=1./(d(i)*dtot)
			ENDDO

			! new potential is inverse distance-weighted sum of edge potentials
			phip=0.
			DO i=1,n
				phip=phip+w(i)*phi(i)
			ENDDO
		CASE DEFAULT
			PRINT *, 'Incorrect LinInterpolate:WeightMethod'
			STOP
		END SELECT
	END SUBROUTINE LinInterpolate

! ***************
! SUBROUTINE InsidePot
! ***************
!	
! DESCRIPTION: calculates whether something is inside or outside a shape

	SUBROUTINE InsidePot(InsideYes, NumNeed, meshi, meshj, NumGrid, NumPan, &
			PanX, PanY, PanBlob, NumBlobPlus)
		! incoming
		INTEGER, INTENT(in)		:: NumGrid, NumPan, NumNeed, NumBlobPlus
		INTEGER, DIMENSION(NumBlobPlus), INTENT(inout)	:: PanBlob
		REAL*8, DIMENSION(NumPan), INTENT(in) 			:: PanX, PanY
		REAL*8, DIMENSION(NumGrid), INTENT(in) 			:: meshi, meshj

		! outgoing
		INTEGER,DIMENSION(NumGrid,NumGrid,2),INTENT(inout):: InsideYes 
			!{r,z,(in/out,WhichBlob)}

		! internal
		INTEGER					:: i, j, z, n, CrossCount, CL, CrossBeforePx, CrossAfterPx, &
			CrossBeforePy, CrossAfterPy, countInt, countIntersecOut, sizeCross, &
			NoDouble, xx, yy, zz, CrossYes, countIntersecIn, countDeadXY, countPan, &
			WhichBlob, tempBlob
		INTEGER, DIMENSION(NumPan,4) 	:: BlobCrossDir
		REAL*8							:: Ax, Bx, Ay, By, Px, Py, TriangleArea, &
			IntersecX, IntersecY
		REAL*8, DIMENSION(NumPan,4) 	:: CrossLocation, NewCrossLocation

		CONTINUE

		XXLoop: do xx = 1,NumGrid
			IF ( 0==MOD(xx,10) ) THEN
				PRINT *, '    % done with "InsidePot"= ', 100.*(xx-1)/NumGrid
			ENDIF
			YYLoop: DO yy = 1,NumGrid
				!IF ( 0==MOD(xx,1) ) THEN
				!PRINT *, '    % done with "InsidePot"= ', 100.*((xx-1)*NumGrid+yy)/ &
				!	(NumGrid*NumGrid)
				!ENDIF
				DO j = 1,4
					DO i = 1,countInt
						NewCrossLocation(i,j)=0.
						BlobCrossDir(i,j)=0
					ENDDO
					DO i = 1,CL
						CrossLocation(i,j)=0. 
					ENDDO
				ENDDO

				CrossCount=0; CL=1; 
				Px = meshi(xx); Py = meshj(yy)
				countPan=0; WhichBlob=1

				Zdo: DO z = 1,NumPan 
					countPan=countPan+1

					IF ( NumPan==z ) THEN 
						Ax=PanX(z); Ay=PanY(z)
						Bx=PanX(1); By=PanY(z)            
					ELSE ! note the shift - in axi coords don't drop to beginning point
						Ax=PanX(z); Ay=PanY(z)   
						Bx=PanX(z+1); By=PanY(z+1)   
					ENDIF

					!TriangeAreaDet(1,:)=[1. 1. 1.];  
					!TriangeAreaDet(2,:)=[Ax Bx Px];  
					!TriangeAreaDet(3,:)=[Ay By Py];
					TriangleArea = 0.5*ABS(Py*Bx - Px*By - Py*Ax + Ax*By + Ay*Px - Ay*Bx)
					TriIf: IF (ABS(TriangleArea)-2.*eps<0) THEN	! if yes, pt is within segment or parallel to it
						CALL InSideYesNo( CrossCount, CrossYes, IntersecX, IntersecY, 'Y', Ax, Ay, &
							Bx, By, Px, Py )
						Cross1: if ( 1==CrossYes ) THEN			! within segment 
							DO zz = 1,NumPan-1
								IF ( ((PanX(zz)-IntersecX <eps .AND. PanX(zz+1)-IntersecX >-eps) .OR. &
									(PanX(zz)-IntersecX >-eps .AND. PanX(zz+1)-IntersecX <eps)) .AND. &
									((PanY(zz)-IntersecY <eps .AND. PanY(zz+1)-IntersecY >-eps) .OR. &
									(PanY(zz)-IntersecY >-eps .AND. PanY(zz+1)-IntersecY <eps)) ) THEN
									InsideYes(xx,yy,1)=0
									InsideYes(xx,yy,2)=WhichBlob
									EXIT
								ENDIF
							END DO
						ENDIF Cross1
					ELSE						! general case   
						CALL InSideYesNo( CrossCount, CrossYes, IntersecX, IntersecY, 'Y', Ax, Ay, &
							Bx, By, Px, Py ) 
						Cross2: IF ( 1==CrossYes ) THEN	! inside shape, crossing left-right
							CrossLocation(CL, 1)=IntersecX
							CrossLocation(CL, 2)=IntersecY
							CrossLocation(CL, 3)=z
							CrossLocation(CL, 4)=WhichBlob
							CL=CL+1
						ENDIF Cross2
						CALL InSideYesNo( CrossCount, CrossYes, IntersecX, IntersecY, 'X', Ax, Ay, &
							Bx, By, Px, Py )
						Cross3: if ( 1==CrossYes ) THEN	 	! inside shape, crossing up-down
							CrossLocation(CL, 1)=IntersecX
							CrossLocation(CL, 2)=IntersecY
							CrossLocation(CL, 3)=z  
							CrossLocation(CL, 4)=WhichBlob
							CL=CL+1
						ENDIF Cross3
					ENDIF TriIf

					IF (countPan==PanBlob(WhichBlob)) THEN
						countPan=0
						WhichBlob=WhichBlob+1
					ENDIF
				ENDDO Zdo					! for z=            
				sizeCross = CL-1

				CheckCrossing: IF ( sizeCross>1 .AND. 1 == InsideYes(xx,yy,1) ) THEN
					CrossBeforePx=0; CrossAfterPx=0
					CrossBeforePy=0; CrossAfterPy=0

					countInt=0; 
					NewCrossLocation(1,1)=0.; NewCrossLocation(1,2)=0.; NewCrossLocation(1,3)=0.
					cIO: DO countIntersecOut = 1,sizeCross		! remove duplicate entries
						NoDouble=1
						IF (sizeCross /= countIntersecOut) THEN
							DO	countIntersecIn = countIntersecOut+1,sizeCross
								IF ( abs(CrossLocation(countIntersecIn,1) - CrossLocation(countIntersecOut,1)) &
									<eps .AND. abs(CrossLocation(countIntersecIn,2) - &
									CrossLocation(countIntersecOut,2))<eps ) THEN
									NoDouble=0
								ENDIF
							END DO
						ENDIF
						CrossTest: IF	( 1==NoDouble ) THEN 
							countInt=countInt+1
							DO j = 1,4
								NewCrossLocation(countInt,j)=CrossLocation(countIntersecOut,j)
							ENDDO
						ENDIF CrossTest
					ENDDO cIO
					DO	countDeadXY = 1,countInt			! vertical crossing
						tempBlob=INT(NewCrossLocation(countDeadXY,4))
						IF ( ABS(Px -NewCrossLocation(countDeadXY,1))<eps ) THEN 
							IF ( Py < NewCrossLocation(countDeadXY,2) ) THEN
								CrossAfterPx=CrossAfterPx+1
								BlobCrossDir(tempBlob,1)=BlobCrossDir(tempBlob,1)+1
							ELSEIF ( Py > NewCrossLocation(countDeadXY,2) ) THEN
								CrossBeforePx=CrossBeforePx+1
								BlobCrossDir(tempBlob,2)=BlobCrossDir(tempBlob,2)+1
							ELSE 
								PRINT *,"Cross dead on X"
								EXIT
							ENDIF
						ELSE                 				! horizontal crossing
							IF( Px < NewCrossLocation(countDeadXY,1) ) THEN
								CrossAfterPy=CrossAfterPy+1
								BlobCrossDir(tempBlob,3)=BlobCrossDir(tempBlob,3)+1
							ELSEIF ( Px > NewCrossLocation(countDeadXY,1) ) THEN 
								CrossBeforePy=CrossBeforePy+1
								BlobCrossDir(tempBlob,4)=BlobCrossDir(tempBlob,4)+1
							ELSE 
								PRINT *,"Cross dead on Y"
								EXIT
							ENDIF
						ENDIF
					END DO

					! outside because even number of crossings
					!IF ( 0==CrossAfterPy .OR. 0==CrossAfterPx .OR. 0==CrossBeforePx .OR. 0==CrossBeforePy ) THEN               
					IF ( 0==CrossAfterPx .OR. 0==CrossBeforePx .OR. MOD(CrossAfterPy,2)/=1) THEN  
						InsideYes(xx,yy,1)=0						

						! find which blob the crossing is in 
						DO n=1,countInt
							IF ( NumNeed==NewCrossLocation(n,4) .AND. &
								1==MOD(BlobCrossDir(n,1)-BlobCrossDir(n,2),2) ) THEN		! needle, crosses up 
								InsideYes(xx,yy,2)=INT(NewCrossLocation(n,4))
							ELSEIF (1==BlobCrossDir(n,1) .AND. 1==BlobCrossDir(n,2)) THEN 	! within, left/right	
								InsideYes(xx,yy,2)=INT(NewCrossLocation(n,4))
							ENDIF
						ENDDO
					ENDIF
				ENDIF CheckCrossing	! sizeCross
			ENDDO YYLoop			! for yy
		ENDDO XXLoop           		! for xx

		xx=CrossAfterPy; xx=CrossBeforePy		! useless, but supresses warnings
	END SUBROUTINE InsidePot

! ***************
! SUBROUTINE InSideYesNo
! ***************
! DESCRIPTION: Determines the horizontal and vertical intersections for pt
!   P onto the line between pts A,B
!
! CALLING PROGRAM: BoundaryForceCalc
!
! INPUT: X,Y sweep direction, examine pt, segment end coordinates, # 
!   previous crossings
!
! OUTPUT: # and location of intersections

	SUBROUTINE InSideYesNo( CrossCount, CrossYes, IntersecX, IntersecY, CrossDir, &
			Ax, Ay, Bx, By, Px, Py )
		! incoming variables
		REAL*8, INTENT(in)		:: Ax, Ay, Bx, By, Px, Py
		CHARACTER, INTENT(in)	:: CrossDir

 		! outgoing variables
		INTEGER, INTENT(out)	:: CrossYes
		INTEGER, INTENT(inout)	:: CrossCount
		REAL*8, INTENT(out)		:: IntersecX, IntersecY

		! subroutine entirely variables
		REAL*8					:: delBAx, delBAy, t

		CONTINUE

		CrossYes=0; t = toolarge
		delBAx = Bx-Ax; delBAy = By-Ay

		SELECT CASE( CrossDir )
		CASE ('Y')
			IF ( ABS(delBAy)-eps>0 ) THEN 
				t = (Py - Ay)/delBAy
				IntersecX = Ax + t*delBAx
				IntersecY = Ay + t*delBAy
			ENDIF
		CASE ('X')
			IF ( ABS(delBAx)-eps>0 ) THEN                
				t = (Px - Ax)/delBAx  
				IntersecX = Ax + t*delBAx
				IntersecY = Ay + t*delBAy
			ENDIF
		CASE DEFAULT
			STOP 'Incorrect CrossDir'
		END SELECT

		IF ( t>=0-eps .AND. t<=1+eps ) THEN
			CrossCount = CrossCount+1; CrossYes=1
		ENDIF
	END SUBROUTINE InSideYesNo

! ***************
! SUBROUTINE CalcGridPot
! ***************
! DESCRIPTION: Calculates the grid points location and potential values
!
! CALLING PROGRAM: BoundaryForceCalc, CalcDropletPot
!
! INPUT: 
!
! OUTPUT: printed potential file
	SUBROUTINE CalcGridPot(meshR, meshZ, numblob, numneed, minR, maxR, minZ, &
			maxZ, numGaussPan, wts, GaussPan, runArbPot, numPtsPan, numblobplus, &
			DistMult)
		! incoming variables
		INTEGER, INTENT(in)		:: numblob, numGaussPan, runArbPot, numneed, numblobplus
		INTEGER, DIMENSION(numblobplus), INTENT(in) 	:: numPtsPan
		REAL*8, INTENT(in)		:: minR, maxR, minZ, maxZ, DistMult
		REAL*8, DIMENSION(npoints), INTENT(in) 			:: wts
		TYPE(PANEL),DIMENSION(numGaussPan),INTENT(in) 	:: GaussPan

 		! outgoing variables - none
		REAL*8, DIMENSION(PotMeshSize), INTENT(inout) 	:: meshR, meshZ

		! subroutine entirely variables
		CHARACTER(LEN=12) 		:: PotCase
		INTEGER, PARAMETER		:: Potunit=12
		INTEGER					:: i, j, k, m
		INTEGER, DIMENSION(:), ALLOCATABLE	:: tempWhichBlob, cumGPpan
		INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: InsideYes
		REAL*8								:: PotOffset
		REAL*8, DIMENSION(:), ALLOCATABLE	:: tempr, tempz
		REAL*8, DIMENSION(:,:), ALLOCATABLE	:: ArbPot

		CONTINUE
		ALLOCATE(ArbPot(PotMeshSize,PotMeshSize))
		ALLOCATE(InsideYes(PotMeshSize,PotMeshSize,2))
		DO i=1,PotMeshSize
			DO j=1,PotMeshSize
				ArbPot(i,j)=0.
				InsideYes(i,j,1)=1; InsideYes(i,j,2)=1
			ENDDO
		ENDDO

		PotOffset=0.01	! fraction of min-max range for calc'ed potentials
		PotDo: DO i=1,PotMeshSize
			!.....	
			CALL GetPotGrid(meshR(i), PotCase, PotOffset, minR, maxR, i)
			!.....	
			Jdo: DO j=1,PotMeshSize
				IF (1==i) THEN
					CALL GetPotGrid(meshZ(j), PotCase, PotOffset, minZ, maxZ, j)
				ENDIF

				IF (1==runArbPot .AND. 0==AformDouble) THEN	
					!.....	// only do if writing entire 'Arbpot.dat'. Otherwise, force interior=Needle
					CALL calcpotwrap(ArbPot(i,j), meshZ(j), meshR(i), numGaussPan, wts, GaussPan)
					!.....	
				ELSEIF (1==runArbPot .AND. 1==AformDouble) THEN
					!.....	// same as above, except different 'GaussPan' formulation
					CALL calcpotdnwrap(ArbPot(i,j), meshZ(j), meshR(i), numGaussPan, wts, GaussPan)
					!.....	
				ENDIF
			ENDDO Jdo
		ENDDO PotDo

		BackIf: IF (0==runArbPot) THEN ! Calcs base background E field, with corrected outside values
			! ** NOTE, will produce incorrect values if droplets present
			IF (numblob>1) THEN
				PRINT *, '   WARNING: Check electric field in ', PotCase, ' for incorrect potentials.'
				PRINT *, '     Error coming from multiple droplets in MatGreens:OptionalCalcPot.'
			ENDIF
			! calculates whether something is inside or outside a shape
			ALLOCATE(tempr(numGaussPan))
			ALLOCATE(tempz(numGaussPan))
			ALLOCATE(tempWhichBlob(numblobplus))
			ALLOCATE(cumGPpan(numblobplus))
			DO i=1,numGaussPan
				tempr(i)=0.; tempz(i)=0.
			ENDDO
			DO i=1,numGaussPan
				tempWhichBlob(i)=0.; cumGPpan(i)=0.
			ENDDO

			m=0
			DO i=1,numblobplus
				j=numPtsPan(i)
				IF (i>1) THEN
					cumGPpan(i)=cumGPpan(i-1)+j
				ELSE
					cumGPpan(i)=j
				ENDIF

				DO k=1,j
					m=m+1
					tempr(m)=GaussPan(m)%r0
					tempz(m)=GaussPan(m)%z0
					tempWhichBlob(i)=j
				ENDDO
			ENDDO

			!.....
			PRINT *, '  NOTE: "InsidePot" changed for axisymmetric version'
			CALL InsidePot(InsideYes, numNeed, meshR, meshZ, PotMeshSize, numGaussPan, &
				tempr, tempz, tempWhichBlob, numblobplus)
			!.....
			DEALLOCATE(tempr); DEALLOCATE(tempz); DEALLOCATE(tempWhichBlob)
			PRINT *, '   Found which spots inside shape'

			PotDo2: DO i=1,PotMeshSize
				Jdo2: DO j=1,PotMeshSize
					IF (1==InsideYes(i,j,1)) THEN
						IF (0==AformDouble) THEN	
							CALL calcpotwrap(ArbPot(i,j), meshZ(j), meshR(i), numGaussPan, wts, GaussPan)
						ELSEIF (1==AformDouble) THEN
							CALL calcpotdnwrap(ArbPot(i,j), meshZ(j), meshR(i), numGaussPan, wts, GaussPan)
						ENDIF
					ELSE	! determine which blob it belongs to and assing that pot
						!PRINT *, InsideYes(i,j,2), cumGPpan(INT(InsideYes(i,j,2))), &
						ArbPot(i,j) = GaussPan(cumGPpan(INT(InsideYes(i,j,2))))%midPot
					ENDIF
				ENDDO Jdo2
			ENDDO PotDo2

			DEALLOCATE(cumGPpan)
		ENDIF BackIf

		IF (1==runArbPot) THEN	
			PRINT *, '  -Writing ArbPot.dat'
			OPEN(unit=Potunit, file='ArbPot.dat', status='replace')
		ELSE
			PRINT *, '  -Writing ', PotCase
			OPEN(unit=Potunit, file=PotCase, status='replace')
		ENDIF
		WRITE(unit=Potunit, fmt='(a,i5,a)') 'TITLE="Potentials using boundary matrix with ', numGaussPan, ' panels"'
		WRITE(unit=Potunit, fmt='(a)') 'VARIABLES="R [cm]"'
		WRITE(unit=Potunit, fmt='(a)') '"Z [cm]"'
		WRITE(unit=Potunit, fmt='(a)') '"Potential [V]"'
		WRITE(unit=Potunit, fmt='(a)') 'ZONE T="main"'
		WRITE(unit=Potunit, fmt='(a,i3,a,i3,a)') 'I=', PotMeshSize, ', J=', PotMeshSize, ', K=1, F=POINT'

		DO i=1,PotMeshSize
			DO j=1,PotMeshSize
				WRITE(unit=Potunit, fmt=13) meshR(i)*DistMult, meshZ(j)*DistMult, ArbPot(i,j)
			ENDDO
		ENDDO
		CLOSE (unit=Potunit, status='keep')

		GOTO 231

		!-----------------------------
		! FORMAT listings
		13 FORMAT(3e14.6)

		231 DEALLOCATE(ArbPot); DEALLOCATE(InsideYes)
	END SUBROUTINE CalcGridPot

! ***************
! SUBROUTINE ForceToMeters
! ***************
! DESCRIPTION: recalculates the total charge and forces on each panel to keep electron
!   *density* (via area) constant between timesteps
!
! INPUT: New panel numbering matrix, array sizes, normal electric fields 
!
! OUTPUT: new panel charges and forces 
	SUBROUTINE ForceToMeters(panarray2, ttlpanels, tomyes)
		! incoming/outgoing
		INTEGER, INTENT(in)							  	:: ttlpanels, tomyes
		TYPE(PANEL), DIMENSION(ttlpanels),INTENT(inout) :: panarray2

		! internal
		INTEGER				:: i, j

		CONTINUE

		SELECT CASE (tomyes)
		CASE (1)	! change cm to m
			DO i=1,ttlpanels
				panarray2(i)%z0=panarray2(i)%z0*cmtom;
				panarray2(i)%z1=panarray2(i)%z1*cmtom;
				panarray2(i)%r0=panarray2(i)%r0*cmtom;
				panarray2(i)%r1=panarray2(i)%r1*cmtom;
				panarray2(i)%midz=panarray2(i)%midz*cmtom;
				panarray2(i)%midr=panarray2(i)%midr*cmtom;
				panarray2(i)%length=panarray2(i)%length*cmtom;
				!panarray2(i)%str=panarray2(i)%str*cmtom;

				DO j=1,npoints
					panarray2(i)%zpoints(j)=panarray2(i)%zpoints(j)*cmtom
					panarray2(i)%rpoints(j)=panarray2(i)%rpoints(j)*cmtom

				ENDDO
			ENDDO
		CASE (0)	! change m to cm
			DO i=1,ttlpanels
				panarray2(i)%z0=panarray2(i)%z0/cmtom;
				panarray2(i)%z1=panarray2(i)%z1/cmtom;
				panarray2(i)%r0=panarray2(i)%r0/cmtom;
				panarray2(i)%r1=panarray2(i)%r1/cmtom;
				panarray2(i)%midz=panarray2(i)%midz/cmtom;
				panarray2(i)%midr=panarray2(i)%midr/cmtom;
				panarray2(i)%length=panarray2(i)%length/cmtom;
				!panarray2(i)%str=panarray2(i)%str/cmtom;

				DO j=1,npoints
					panarray2(i)%zpoints(j)=panarray2(i)%zpoints(j)/cmtom
					panarray2(i)%rpoints(j)=panarray2(i)%rpoints(j)/cmtom
				ENDDO
			ENDDO
		CASE DEFAULT
			PRINT *, 'Err, MAT_Greens:ForceToMeters:tomyes'
			STOP
		END SELECT
	END SUBROUTINE ForceToMeters

! ***************
! SUBROUTINE ZeroorCopyPanel
! ***************
! DESCRIPTION: for a derived type of structure 'panel', initializes it to zeros
!   or copies from one instance to another
!
! INPUT: old and new Panel, nullify or copy
!
! OUTPUT: zeroed or copies instance
	SUBROUTINE ZeroorCopyPanel(NewPanel, Nullyes, OldPanel, reverseyes)
		INTEGER, INTENT(in)			:: Nullyes
		INTEGER, INTENT(in), OPTIONAL		:: reverseyes
		TYPE(PANEL), INTENT(in), OPTIONAL 	:: OldPanel

		TYPE(PANEL), INTENT(inout) 	:: NewPanel

		INTEGER						:: j
		REAL						:: flippedPan

		CONTINUE

		flippedPan=1.
		IF (PRESENT(reverseyes) .AND. 1==reverseyes) flippedPan=-1.

		SELECT CASE(Nullyes)
		CASE(0)	! copy the panel component, single instance
			NewPanel%type1	=OldPanel%type1
			NewPanel%z0		=OldPanel%z0
			NewPanel%z1		=OldPanel%z1
			NewPanel%r0		=OldPanel%r0
			NewPanel%r1		=OldPanel%r1
			NewPanel%str	=OldPanel%str
			NewPanel%znrm	=flippedPan*OldPanel%znrm
			NewPanel%rnrm	=flippedPan*OldPanel%rnrm
			NewPanel%midPot	=OldPanel%midPot
			NewPanel%midz	=OldPanel%midz
			NewPanel%midr	=OldPanel%midr
			NewPanel%length	=OldPanel%length

			IF (PRESENT(reverseyes) .AND. 1==reverseyes) THEN
				DO j=1,npoints
					NewPanel%zpoints(j)=OldPanel%zpoints(npoints-j+1)
					NewPanel%rpoints(j)=OldPanel%rpoints(npoints-j+1)
					NewPanel%rscaled(j)=OldPanel%rscaled(npoints-j+1)
				ENDDO
			ELSE
				DO j=1,npoints
					NewPanel%zpoints(j)=OldPanel%zpoints(j)
					NewPanel%rpoints(j)=OldPanel%rpoints(j)
					NewPanel%rscaled(j)=OldPanel%rscaled(j)
				ENDDO
			ENDIF

		CASE(1)	! nullify the panel
			IF (PRESENT(reverseyes) .OR. PRESENT(OldPanel)) THEN
				PRINT *, 'Err Panels:ZeroorCopyPanel case choice.'
				STOP
			ENDIF
			NewPanel%type1=0.
			NewPanel%z0=0.
			NewPanel%z1=0.
			NewPanel%r0=0.
			NewPanel%r1=0.
			NewPanel%str=0.
			NewPanel%znrm=0.
			NewPanel%rnrm=0.
			NewPanel%midPot=0.
			NewPanel%midz=0.
			NewPanel%midr=0.
			NewPanel%length=0.

			DO j=1,npoints
				NewPanel%zpoints(j)=0.
				NewPanel%rpoints(j)=0.
				NewPanel%rscaled(j)=0.
			ENDDO
		CASE DEFAULT
			PRINT *, 'Panels:ZeroorCopyPanel err'
			STOP
		END SELECT
	END SUBROUTINE ZeroorCopyPanel

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
! CALLING PROGRAM: MAT_Greens:SortNN_NoCross
!
! ALGORITHM: compare sign of ( (ad x ab) v. (ac x ab) ). if same, no cross

	SUBROUTINE DoPtsCross(  DoCross, x1, y1, x2, y2, x3, y3, x4, y4 )
		! incoming variables
		REAL*8, INTENT(in)	:: x1, y1, x2, y2, x3, y3, x4, y4

		! outgoing variables
		LOGICAL,INTENT(OUT) :: DoCross

		! start regular variables
		REAL*8				:: s1lx, s1rx, s1ly, s1ry, lsign, rsign, &
			s2lx, s2rx, s2ly, s2ry

		CONTINUE 

		IF ( (x1==x3 .AND. y1==y3) .OR. (x1==x4 .AND. y1==y4) .OR. &
			(x2==x3 .AND. y2==y3) .OR. (x2==x4 .AND. y2==y4) ) THEN
			DoCross=.FALSE.				! no cross, since only 3 points instead of 4
		ELSE
			! http://softsurfer.com/Archive/algorithm_0108/algorithm_0108.htm
			IF (x1<=x2) THEN
				s1lx=x1; s1rx=x2;
				s1ly=y1; s1ry=y2;
			ELSE
				s1lx=x2; s1rx=x1;
				s1ly=y2; s1ry=y1;
			ENDIF

			IF (x3<=x4) THEN
				s2lx=x3; s2rx=x4;
				s2ly=y3; s2ry=y4;
			ELSE
				s2lx=x4; s2rx=x3;
				s2ly=y4; s2ry=y3;
			ENDIF

			! test for existence of an intersect point
			CALL isLeft(lsign, s1lx, s1ly, s1rx, s1ry, s2lx, s2ly)
			CALL isLeft(rsign, s1lx, s1ly, s1rx, s1ry, s2rx, s2ry)

			IF (lsign * rsign >= 0) THEN !s2 endpoints have same sign relative to s1
				DoCross=.FALSE.			! => on same side => no intersect is possible
			ELSE
				CALL isLeft(lsign, s2lx, s2ly, s2rx, s2ry, s1lx, s1ly)
				CALL isLeft(rsign, s2lx, s2ly, s2rx, s2ry, s1rx, s1ry)

				IF (lsign * rsign >= 0) THEN ! s1 endpoints have same sign relative to s2
					DoCross=.FALSE.     ! => on same side => no intersect is possible
				ELSE					! the segments s1 and s2 straddle each other
					DoCross=.TRUE.   	! => an intersect exists
				ENDIF
			ENDIF
		ENDIF
CONTAINS
		! tests if point P2 is Left|On|Right of the line P0 to P1.
		!   returns: >0 for left, 0 for on, and <0 for right of the line.
		SUBROUTINE isLeft (LR, P0x, P0y, P1x, P1y, P2x, P2y)
			REAL*8, INTENT(in)	:: P0x, P0y, P1x, P1y, P2x, P2y
			REAL*8, INTENT(out) :: LR
			CONTINUE

			LR= (P1x - P0x)*(P2y - P0y) - (P2x - P0x)*(P1y - P0y)
		END SUBROUTINE isLeft
	END SUBROUTINE DoPtsCross
END MODULE Panels
