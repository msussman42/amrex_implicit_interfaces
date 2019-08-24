MODULE Shading

! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   24 Jul 2005 begun

	USE ShadingGlobal
	USE PtFitting, only 	: CurveFit, GAUSS, gaussvals, TrapInit, LEGS, Nexti, CalcPDF, &
		DoPtsCross, polygon_area_2d_2, PerturbWave
	USE GMRes_Andrew, only 	: gmresIter

	IMPLICIT NONE

! used to restrict variables to prevent accidental calling outside
	PRIVATE      ! everything defaults to 'private' except those expressly labeled otherwise.
	PUBLIC		:: BoundaryForceCalc
CONTAINS 

! ***************
! Subroutine BoundaryForceCalc
! ***************
! DESCRIPTION: calculates the force on panels using the boundary integral method.
!
! INPUT: Level set grid, grid sizing and edges
!
! OUTPUT: roughly equally spaced points along the surface (zero LS) with the applied force
	SUBROUTINE BoundaryForceCalc( BlobNum, blob_xc, blob_yc, NewNumPtsPan, MxPt_blob, radBlob, &
			thetaBlob, forceBlob, tempphi, LSMaxX, LSMaxY, is_symmetric, xlo, ylo, &
			dx, dy, elecHeight,ptime)	! formal full version

		! incoming variables
		INTEGER, INTENT(in)					:: LSMaxX, LSMaxY, is_symmetric
		REAL*8, INTENT(in)					:: xlo, ylo, dx, dy, elecHeight, ptime
		REAL*8, DIMENSION(LSMaxX, LSMaxY), INTENT(in) 	:: tempphi

		! outgoing variables
		INTEGER, INTENT(out)				:: BlobNum, MxPt_blob
		INTEGER, DIMENSION(SmallInitAry), INTENT(out)	:: NewNumPtsPan
		REAL*8, DIMENSION(SmallInitAry), INTENT(out) 	:: blob_xc, blob_yc
		REAL*8, DIMENSION(SmallInitAry, InitArySize), INTENT(out)	:: radBlob, thetaBlob, forceBlob

		! program variables
  		! fixed variables
		INTEGER			:: Potunit, Forceunit, Coordsunit, NewCoordsunit, NumBlobunit, &
			tempPanunit, NewDropunit, NewDropunit2, tempLSunit, Phiunit, Phititleunit, &
			PDFsingleunit, PDFhighunit, Resunit
		REAL*8	:: electrodeWidth	! 0.3cm = 0.003m

  		! regular variables
		CHARACTER(len=15) 					:: LSstep
		LOGICAL								:: NotFirst
		INTEGER			:: NumCoords, i, j, k, z, panel, PanFromTo, NeedleNum, &
			BlobPlus, NumBlobCoords, iPlus, i2, ClosePanCount, SkipPanTo, tempRange, &
			totPanels, NewdropPts, xLeftSmoothEdge
		INTEGER, DIMENSION(SmallInitAry) 	:: sizeLevel
		INTEGER, DIMENSION(:), ALLOCATABLE	:: indxLegs, NumPanels
		INTEGER, DIMENSION(:,:),ALLOCATABLE :: InsideYes

		REAL*8							:: delX, delY, GreenPart, GreenPart2, ai, bi, ci, di, t2, &
			t1, t3, timeCount, elecLedge, elecRedge, maxX, minX, maxY, minY, dropVol, chargeMassArea, &
			chargeMassVol, chi, MTCR, droprad, aa, ref, c, ThreePtAngle, x, y, minX2, totalTime, zero
		REAL*8, DIMENSION(2)			:: PtRef, PtC, PtA
		REAL*8, DIMENSION(M)			:: gpts, gwgt
		REAL*8, DIMENSION(:), ALLOCATABLE		:: xPanLen, yPanLen, PanLen, xPan, yPan, &
			xPanCen, yPanCen, xPanBegin, yPanBegin, delXCoord, normPanX, normPanY, x1, y1, &
			RHS_gamma_phi, RHS_gamma_dphi, RHS_alpha_phi, RHS_alpha_dphi, RHStemp, sigmaGauss, &
			meshX, meshY, PanCharge, PanArea, Along, temp1, temp2, temp3
		REAL*8, DIMENSION(:,:), ALLOCATABLE		:: Coords, RHS, A, ArbPot, PtsForce, &
			temp4, temp5
			! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatpan]

		CONTINUE
  		!**********
  		! main program commands begin
  		!**********
  		! defining parameters in main program 
		PRINT *, '********************'
		CALL cpu_time(t3)
		PRINT *, 'LSMaxX, LSMaxY, is_symmetric =',LSMaxX, LSMaxY, is_symmetric
                PRINT *, 'ptime=',ptime
		WRITE(*, fmt='(a,3e11.4)') '  dx, dy, elecHeight =', dx, dy, elecHeight

		! handle definitions
		Potunit=11; Forceunit=12; Coordsunit=13; NewCoordsunit=14; NumBlobunit=15
		tempPanunit=16; NewDropunit=17; NewDropunit2=18; tempLSunit=19; Phiunit=20
		Phititleunit=21; PDFsingleunit=22; PDFhighunit=23; Resunit=24

		totalTime=0
		electrodeWidth=0.2	! 0.3cm = 0.003m

		counter=counter+1		! makes it possible to show how many times routine called
		ALLOCATE(phi(LSMaxX, LSMaxY))

		!-----------------------------
		IF (1==LoadPastStep) THEN
			! debugging step, loading latest number step
			! writes a data file with the step number and electric field
			counter=500		! makes it possible to show how many times routine called

			! convert data number to string file name
			CALL Int2Str(LSstep, 'visual00000', counter,'',0)
			!LSstep='Phi_pass.dat'

			OPEN(unit=tempLSunit, file=LSstep, status='old')
			panel=0
			DO i =1,LSMaxX
				DO j=1,LSMaxY
					panel=panel+1
					READ(UNIT=tempLSunit, FMT=12) k, k, ai, phi(i,j), ai
				ENDDO
			ENDDO	
			CLOSE (UNIT=tempLSunit, STATUS='keep')
			PRINT *, '      ', LSstep, ' read.'
		ELSE
			phi=tempphi 	! use unaltered original data
		ENDIF
		ShowPhi:IF (1==ShowPassedPhi) THEN
				! write inputted 'phi to a file'
			PRINT *, '  -Writing "Phi_pass.dat" '
			OPEN(unit=Phiunit, file='Phi_pass.dat', status='replace')
			OPEN(unit=Phititleunit, file='Phi_pass2.dat', status='replace')	

			WRITE(unit=Phititleunit, fmt='(a)') 'VARIABLES="X [cm]"'
			WRITE(unit=Phititleunit, fmt='(a)') '"Y [cm]"'
			WRITE(unit=Phititleunit, fmt='(a)') '"phi"'
			WRITE(unit=Phititleunit, fmt='(a)') 'ZONE T="main2"'
			WRITE(unit=Phititleunit, fmt='(a,i4,a,i4,a)') 'I=', LSMaxY, ', J=', LSMaxX, ', K=1, F=POINT'

			DO i=1,LSMaxX
				x = (i+0.5)*dx + xlo
				DO j=1,LSMaxY
					y = (j+0.5)*dy + ylo
					WRITE(unit=Phiunit, FMT=12) i, j, x, phi(i,j), y 
					WRITE(unit=Phititleunit, fmt='(3e14.6)') x, y, phi(i,j)
				ENDDO
			ENDDO
			CLOSE (unit=Phiunit, status='keep')
			CLOSE (unit=Phititleunit, status='keep')
		ENDIF ShowPhi

		!-----------------------------
		! add random fluctuation to the surface
		PerturbIf: IF (UsePerturb .AND. 1==counter) THEN
			! confirm that it is the first file and not a restart by confirming there is no visual00001
			OPEN(unit=Phititleunit, file='visual00001', status='old', err=801)
			PRINT *, '    Skipping perturbing phi because appears to not be the simulation start.'
			GOTO 802	! clunky, but effectively skips this section if file present

			801 ALLOCATE(x1(LSMaxX)); ALLOCATE(y1(LSMaxY)); 
			CALL PerturbWave ( phi, x1, y1, LSMaxX, LSMaxY, dx, dy )

			PRINT *, '  -Writing "Phi_bump.dat" '
			OPEN(unit=Phiunit, file='Phi_bump.dat', status='replace')
			WRITE(unit=Phiunit, fmt='(a)') 'VARIABLES="X [cm]"'
			WRITE(unit=Phiunit, fmt='(a)') '"Y [cm]"'
			WRITE(unit=Phiunit, fmt='(a)') '"phi"'
			WRITE(unit=Phiunit, fmt='(a)') 'ZONE T="main2"'
			WRITE(unit=Phiunit, fmt='(a,i4,a,i4,a)') 'I=', LSMaxY, ', J=', LSMaxX, ', K=1, F=POINT'

			DO i=LSMaxX,1,-1
				DO j=1,LSMaxY
					WRITE(unit=Phiunit, fmt='(3e14.6)') x1(i), y1(j), phi(i,j)
				ENDDO
			ENDDO	

			DEALLOCATE(x1);	DEALLOCATE(y1)
			CLOSE (unit=Phiunit, status='keep')
		ENDIF PerturbIf
		802 ALLOCATE(Coords(InitArySize,8)); Coords(:,:)=-1; sizeLevel(:)=-1
		ALLOCATE(alreadyChecked(LSMaxX,LSMaxY)); alreadyChecked(:,:)=-1 ! nothing is initially checked
		ALLOCATE(insideMarker(LSMaxX,LSMaxY)); insideMarker(:,:)=-1 ! nothing is initially marked inside
		ChokeCoord(:,:)=0; minX2=1e6

		! determine all the surface points from passed levelset
		!.....
		CALL ShapeGrow( Coords, BlobNum, NeedleNum, MxPt_blob, sizeLevel, NumBlobCoords, &
			LSMaxX, LSMaxY, xlo, ylo, dx, dy, is_symmetric, minX2 )
		!.....
		DEALLOCATE(insideMarker); DEALLOCATE(alreadyChecked)
		NumCoords= NumBlobCoords+20

		!-----------------------------	
		! add side N. BC and electrodes to all blobs
		BlobPlus=-1; elecLedge=-1.0; elecRedge=-1.0
		ALLOCATE(temp5(NumCoords,8)); temp5(:,:)=Coords(1:NumCoords,:)
	    !.....	
		CALL AddBCandAround ( temp5, sizeLevel, NumBlobCoords, BlobPlus, BlobNum, &
			elecHeight, electrodeWidth, elecLedge, elecRedge, maxX, minX, maxY, minY )
		!.....
		Coords(1:NumCoords,:)=temp5
		DEALLOCATE(temp5)

		!-----------------------------	
		! switch so the points go COUNTERclockwise around. 
		! Also, define the fluid coordinates with xarray, yarray and panels 
  		!   with fixed potential (-2.0) and dirichlet boundary conditions(==1) and
		!   fixed flux (0) and neumann bc (==0) for detached droplets
		k=0; ALLOCATE(temp1(8))
		DO i=1,BlobNum
			i2=k
			DO j=1,sizeLevel(i)
				k=k+1
				IF (NeedleNum==i) THEN
					Coords(k,3) = fluidPot
					Coords(k,4) = 1
				ELSE
					Coords(k,3) = 0
					Coords(k,4) = 0

					IF (j<=INT( FLOOR(sizeLevel(i)/2.0) )) THEN
						temp1(:)=Coords(k,:)
						z=i2+sizeLevel(i)-j+1

						Coords(k,:)=Coords(z,:)
						Coords(z,:)=temp1(:)
					ENDIF
				ENDIF
			ENDDO
		ENDDO
		DEALLOCATE(temp1)
		WRITE(*, fmt=17) '   The left, right electrode middle edges are at x=', elecRedge, elecLedge

		IF (1==ShowCoords) THEN
			! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, open shape]
			PRINT *, '  -Writing "Coords.dat" '
			OPEN(unit=Coordsunit, file='Coords.dat', status='replace')
			WRITE(unit=Coordsunit, fmt='(a)') 'VARIABLES="X [cm]"'
			WRITE(unit=Coordsunit, fmt='(a)') '"Y [cm]"'
			WRITE(unit=Coordsunit, fmt='(a)') '"Pt val"'
			WRITE(unit=Coordsunit, fmt='(a)') '"Pt type"'
			WRITE(unit=Coordsunit, fmt='(a)') '"Blob#"'
			WRITE(unit=Coordsunit, fmt='(a)') '"Blob_xc"'
			WRITE(unit=Coordsunit, fmt='(a)') '"Blob_yc"'
			WRITE(unit=Coordsunit, fmt='(a)') '"OpenShape"'
			WRITE(unit=Coordsunit, fmt='(a)') '"Pan #"'
			DO i=1,NumCoords
				WRITE(unit=Coordsunit, fmt='(3e14.6, 2i4, 2e14.6, 2i5)') Coords(i,1), Coords(i,2), &
					Coords(i,3), INT(Coords(i,4)), INT(Coords(i,5)), Coords(i,6), Coords(i,7), &
					INT(Coords(i,8)), i
			ENDDO
			CLOSE (unit=Coordsunit, status='keep')
		ENDIF		

		ALLOCATE(delXCoord(NumCoords))
		DO i=1,NumCoords
			iPlus=Nexti( i, sizeLevel(INT(Coords(i,5))), INT(Coords(i,8)) )
			delXCoord(i)=Coords(iPlus,1)-Coords(i,1) 
		ENDDO

		!-----------------------------	
		! determine the gaussian quadrature weights and locations along each {-1:1} panel
		CALL gaussvals(gpts, gwgt)

		!-----------------------------	
  		! create spaced points along the spline. Returns NewCoords
		ALLOCATE(NumPanels(BlobPlus)); NumPanels(:)=-1
		NewNumPtsPan(:)=-1; NewCoords(:,:)=-1
		i2=0
		DO i=1,BlobPlus
			IF (i>1) THEN
				NumPanels(i)=NumPanels(i-1)
				i2=i2+sizeLevel(i)
			ELSE
				NumPanels(i)=0
				i2=sizeLevel(i)
			ENDIF

			j=1
			DO WHILE (INT(Coords(j,5))/=i)
				j=j+1
				IF (j>=NumCoords) THEN
					PRINT *, 'Shading:BoundaryForceCalc Call curvefit err. i=',i
					STOP
				ENDIF
			ENDDO

			! vary the last parameter by droplet (closed surf) or needle/boundary (open) type
			IF (i>BlobNum .OR. NeedleNum==i) THEN	! open(flat) type
				Coords(j,8)=1
			ELSE		
				Coords(j,8)=0						! closed type
			ENDIF

			ALLOCATE(temp1(sizeLevel(i))); ALLOCATE(temp2(sizeLevel(i)))
			ALLOCATE(temp3(sizeLevel(i))); ALLOCATE(temp4(sizeLevel(i),2))
			ALLOCATE(temp5(sizeLevel(i),3))

			tempRange=i2
			temp1=Coords(j:tempRange,1); temp2=Coords(j:tempRange,2)
			temp3=delXCoord(j:tempRange); temp4=Coords(j:tempRange,3:4)
			temp5=Coords(j:tempRange,5:7)

			PRINT *, '    blob#,pt#, sizeblob=',i,j,sizeLevel(i)
			!-------------
			! creates 'newcoords', # of points/panel and linear/cubic fit
			CALL CurveFit ( NumPanels(i), NewNumPtsPan(i), temp1, &
				temp2, temp3, sizeLevel(i), temp4, temp5, &
				INT(Coords(j,8)), i, BlobNum )

			DEALLOCATE(temp1); DEALLOCATE(temp2); DEALLOCATE(temp3)
			DEALLOCATE(temp4); DEALLOCATE(temp5)
		ENDDO
		totPanels=NumPanels(BlobPlus)
		PRINT *, 'Done after CurveFit!'

		SNC: IF (1==ShowNewCoords) THEN
	  		! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatpan]
			PRINT *, '  -Writing "NewCoords.dat" '
			OPEN(unit=NewCoordsunit, file='NewCoords.dat', status='replace')
			WRITE(unit=NewCoordsunit, fmt='(a)') 'VARIABLES="X [cm]"'
			WRITE(unit=NewCoordsunit, fmt='(a)') '"Y [cm]"'
			WRITE(unit=NewCoordsunit, fmt='(a)') '"Pt val"'
			WRITE(unit=NewCoordsunit, fmt='(a)') '"Pt type"'
			WRITE(unit=NewCoordsunit, fmt='(a)') '"Blob#"'
			WRITE(unit=NewCoordsunit, fmt='(a)') '"Blob_xc"'
			WRITE(unit=NewCoordsunit, fmt='(a)') '"Blob_yc"'
			WRITE(unit=NewCoordsunit, fmt='(a)') '"Open"'
			WRITE(unit=NewCoordsunit, fmt='(a)') '"Panel #"'
			panel=0
			DO i=1,BlobPlus
				DO j=1,NewNumPtsPan(i)
					panel=panel+1
					WRITE(unit=NewCoordsunit, fmt='(3e14.6, 2i4, 2e14.6, 2i5)') NewCoords(panel,1), &
						NewCoords(panel,2), NewCoords(panel,3), INT(NewCoords(panel,4)), &
						INT(NewCoords(panel,5)), NewCoords(panel,6), NewCoords(panel,7), &
						INT(NewCoords(panel,8)), panel
				ENDDO
			ENDDO
			CLOSE (unit=NewCoordsunit, status='keep')
		ENDIF SNC
		DEALLOCATE(delXCoord)

		!-----------------------------	
  		! define parts of panels, locations, normals and colocation pt. weights		
		ALLOCATE(xPan(totPanels*M)); ALLOCATE(yPan(totPanels*M))
		ALLOCATE(xPanBegin(totPanels)); xPanBegin(:)=0.
		ALLOCATE(yPanBegin(totPanels)); yPanBegin(:)=0.
		ALLOCATE(xPanLen(totPanels)); ALLOCATE(yPanLen(totPanels))
		ALLOCATE(PanLen(totPanels)); PanLen(:)=0.
		ALLOCATE(xPanCen(totPanels)); xPanCen(:)=0.
		ALLOCATE(yPanCen(totPanels)); yPanCen(:)=0.
		ALLOCATE(normPanX(totPanels)); normPanX(:)=0.
		ALLOCATE(normPanY(totPanels)); normPanY(:)=0.
		ALLOCATE(RHS(totPanels,2))

		ClosePanCount=0; NotFirst=.FALSE.; xLeftSmoothEdge=0
		PanDefineI: DO i = 1,totPanels
			ClosePanCount=ClosePanCount+1
			itemp: IF ( i<totPanels ) THEN
				i2=i+1
			ELSE
				i2=1
			ENDIF itemp
			IF ( INT(NewCoords(i,5)) == INT(NewCoords(i2,5)) ) THEN		! middle of blob
				iPlus=i+1
			ELSE
				iPlus=i

				NewNumPtsPan(INT(NewCoords(i,5)))=NewNumPtsPan(INT(NewCoords(i,5)))-1
				xPanCen(ClosePanCount-1)= 0.5*(NewCoords(i,1)+NewCoords(i-1,1))
				yPanCen(ClosePanCount-1)= 0.5*(NewCoords(i,2)+NewCoords(i-1,2))
				ClosePanCount=ClosePanCount-1
			ENDIF
			IF ( NeedleNum==INT(NewCoords(i,5)) .AND. NewCoords(i,1)<-0.23) THEN
				! finds where the left boundary of arbitrary forced smoothing in the needle is located. Later run
				xLeftSmoothEdge=xLeftSmoothEdge+1
			ENDIF

			OpenEdge: IF ( iPlus /= i ) THEN
				delX = NewCoords(iPlus,1) - NewCoords(i,1)
				delY = NewCoords(iPlus,2) - NewCoords(i,2)

				xPanLen(ClosePanCount)  = delX; yPanLen(ClosePanCount) = delY
				PanLen(ClosePanCount)   = SQRT( xPanLen(ClosePanCount)*xPanLen(ClosePanCount) + &
					yPanLen(ClosePanCount)*yPanLen(ClosePanCount) )
				normPanX(ClosePanCount) = -delY/PanLen(ClosePanCount) 
				normPanY(ClosePanCount) = delX/PanLen(ClosePanCount)

				xPanBegin(ClosePanCount)= NewCoords(i,1)
				yPanBegin(ClosePanCount)= NewCoords(i,2)

				IF (0==PanLen(ClosePanCount)) THEN
					PRINT *, 'Shading:BFC PanLen error. i,iPlus,k,ClosePanCount=', i,iPlus,k,ClosePanCount
					STOP 
				ENDIF

				RHS(ClosePanCount,1) = NewCoords(i,3) ! pan value				
				RHS(ClosePanCount,2) = NewCoords(i,4) ! pan type

				! finding the centers of the panels
				IF (i>1) THEN
					IF (INT(NewCoords(i,5)) == INT(NewCoords(i-1,5))) THEN
						NotFirst=.TRUE.
					ENDIF
				ENDIF
				IF (NotFirst) THEN
					xPanCen(ClosePanCount-1)=(xPanBegin(ClosePanCount)+xPanBegin(ClosePanCount-1))/2
					yPanCen(ClosePanCount-1)=(yPanBegin(ClosePanCount)+yPanBegin(ClosePanCount-1))/2 

					IF ( INT(NewCoords(i,5)) /= INT(NewCoords(i2,5)) ) THEN	! reached end of blob; wrap around
						IF ( 0==INT(NewCoords(i,5)) ) THEN
							xPanCen(ClosePanCount)= (xPanBegin(ClosePanCount) + &
								xPanBegin(ClosePanCount-NewNumPtsPan( INT(NewCoords(i,5)) )+2))/2
							yPanCen(ClosePanCount)= (yPanBegin(ClosePanCount) + &
								yPanBegin(ClosePanCount-NewNumPtsPan( INT(NewCoords(i,5)) )+2))/2
						ELSE
							xPanCen(ClosePanCount)= xPanBegin(ClosePanCount) + &
								xPanLen(ClosePanCount)/2.
							yPanCen(ClosePanCount)= yPanBegin(ClosePanCount) + &
								yPanLen(ClosePanCount)/2.
						ENDIF
					ENDIF
					NotFirst=.FALSE.
				ENDIF

				! Spread the collocation points IN EACH PANEL according to Gaussian quadrature
				PanDefineTo2: DO j = 1,M
					panel= (ClosePanCount-1)*M+j

					xPan(panel)=(xPanLen(ClosePanCount)+2.*xPanBegin(ClosePanCount))/2. + &
						gpts(j)*xPanLen(ClosePanCount)/2.
					yPan(panel)=(yPanLen(ClosePanCount)+2.*yPanBegin(ClosePanCount))/2. + &
						gpts(j)*yPanLen(ClosePanCount)/2.
				ENDDO PanDefineTo2
			ENDIF OpenEdge
		ENDDO PanDefineI  ! for i
		PRINT *, '  Done after full panel definition with ClosePanCount=',ClosePanCount

		ShowPanDetails: IF (1==ShowSubPan) THEN
			PRINT *, '  -Writing "PanDetail.dat" '
			OPEN(unit=tempPanunit, file='PanDetail.dat', status='replace')
			WRITE(unit=tempPanunit, fmt='(a)') 'VARIABLES="X [cm]"' 
			WRITE(unit=tempPanunit, fmt='(a)') '"Y [cm]"'	
			WRITE(unit=tempPanunit, fmt='(a)') '"PanValue"'
			WRITE(unit=tempPanunit, fmt='(a)') '"PanType"'
			WRITE(unit=tempPanunit, fmt='(a)') '"Panel #"'		
			WRITE(unit=tempPanunit, fmt='(a)') 'ZONE T="main"'
			panel=0
			DO i = 1,ClosePanCount
				!WRITE(unit=tempPanunit, fmt='(3e14.6, 2i6)') xPanCen(i), yPanCen(i), &
				!	RHS(i,1), INT(RHS(i,2)), i
				DO k=1,M
					panel=panel+1
					WRITE(unit=tempPanunit, fmt='(3e14.6, 2i6)') xPan(panel), yPan(panel), &
						RHS(i,1), INT(RHS(i,2)), panel
				ENDDO
			ENDDO
			CLOSE (unit=tempPanunit, status='keep')
		ENDIF ShowPanDetails

		ALLOCATE( A(ClosePanCount,ClosePanCount) ); ALLOCATE(RHStemp(ClosePanCount))
		ALLOCATE(RHS_gamma_phi(ClosePanCount)); ALLOCATE(RHS_gamma_dphi(ClosePanCount)) 
		ALLOCATE(RHS_alpha_phi(ClosePanCount)); ALLOCATE(RHS_alpha_dphi(ClosePanCount))

		RHS_gamma_phi(:)=0; RHS_gamma_dphi(:)=0
		RHS_alpha_phi(:)=0; RHS_alpha_dphi(:)=0
		GreenJ: DO j = 1,ClosePanCount
			GreenI: DO i = 1,ClosePanCount
				PanFromTo = int( RHS(j,2)*10 + RHS(i,2) )
				A(j,i)=0.
				ai=0; bi=0; ci=0; di=0

				DO k = 1,M
					panel=(i-1)*M + k
					delX=xPanCen(j)-xPan(panel)
					delY=yPanCen(j)-yPan(panel)

					SELECT CASE (PanFromTo)		! 1=dirichlet, 0=neumann type BC
					CASE (11)
						GreenPart=0
						CALL TrapInit ( GreenPart, delX, delY, normPanX(i), &
							normPanY(i), normPanX(j), normPanY(j), gwgt(k), &
							PanLen(i), 10, 0, k, M )

						ai = ai + GreenPart

						CALL TrapInit ( GreenPart, delX, delY, normPanX(i), &
							normPanY(i), normPanX(j), normPanY(j), gwgt(k), &
							PanLen(i), 11, 0, k, M )

						RHS_gamma_phi(j) = RHS_gamma_phi(j) + RHS(i,1)*GreenPart

					CASE (10)
						GreenPart=0
						CALL TrapInit ( GreenPart, delX, delY, normPanX(i), &
							normPanY(i), normPanX(j), normPanY(j), gwgt(k), &
							PanLen(i), 10, 0, k, M )

						RHS_gamma_dphi(j) = RHS_gamma_dphi(j) + RHS(i,1)*GreenPart

						CALL TrapInit ( GreenPart, delX, delY, normPanX(i), &
							normPanY(i), normPanX(j), normPanY(j), gwgt(k), &
							PanLen(i), 11, 0, k, M )

						bi = bi - GreenPart 		  

					CASE (1)
						GreenPart=0
						CALL TrapInit ( GreenPart, delX, delY, normPanX(i), &
							normPanY(i), normPanX(j), normPanY(j), gwgt(k), &
							PanLen(i), 11, 0, k, M )

						RHS_alpha_phi(j) = RHS_alpha_phi(j) + RHS(i,1)*GreenPart

						CALL TrapInit ( GreenPart, delX, delY, normPanX(i), &
							normPanY(i), normPanX(j), normPanY(j), gwgt(k), &
							PanLen(i), 10, 0, k, M )

						ci = ci + GreenPart 

					CASE (0)
						GreenPart=0
						CALL TrapInit ( GreenPart, delX, delY, normPanX(i), &
							normPanY(i), normPanX(j), normPanY(j), gwgt(k), &
							PanLen(i), 10, 0, k, M )

						RHS_alpha_dphi(j) = RHS_alpha_dphi(j) + RHS(i,1)*GreenPart

						IF (i /=j) THEN
							CALL TrapInit ( GreenPart, delX, delY, normPanX(i), &
								normPanY(i), normPanX(j), normPanY(j), gwgt(k), &
								PanLen(i), 00, 0, k, M )

							di = di + GreenPart
						ELSEIF ( M==k ) THEN
							di = 0.5
						ENDIF 
					CASE DEFAULT
						PRINT *, 'Shading:BFC PanFromTo aLine err. j,i,k,PanFromTo=', j,i,k,PanFromTo
						STOP
					END SELECT ! switch PanFromTo
				ENDDO		! for k

				A(j,i) = ai + bi + ci + di
			ENDDO GreenI   ! for i   
		ENDDO GreenJ		! for j

		!RHS(:,1) = RHS(:,1) + RHS_gamma_phi(:) - RHS_gamma_dphi(:) + RHS_alpha_phi(:) - RHS_alpha_dphi(:)
		DO j = 1,ClosePanCount
			IF (1==INT(RHS(j,2))) THEN
				RHStemp(j) = RHS_gamma_phi(j) - RHS_gamma_dphi(j) -0.5*RHS(j,1)
			ELSEIF (0==INT(RHS(j,2))) THEN
				RHStemp(j) = 1.*( RHS_alpha_phi(j) - RHS_alpha_dphi(j) )
			ELSE 
				PRINT *, 'Shading: BFC RHS empty err. j,ClosePanCount,RHS=',j,ClosePanCount,RHS(j,2)
				STOP
			ENDIF
		ENDDO

		!-----------------------------
  		! solves for {X} in Ax=b
		WRITE(*, fmt=18), '   Done up to matrix Ax=b with ', ClosePanCount, ' panels.'

		CALL cpu_time(t1)
		SELECT CASE( MatSolveMethod )
		CASE (1)
    		!---- gmres iterative matrix solver
			!----------------------------------
			ALLOCATE( sigma(ClosePanCount) ) 

			! if you want to call past calculated sigmas to make a good guess about future values
			IF (1==ShowSigmas) THEN
				PRINT *, '  -Loading "PastSigmas.dat" '
				OPEN(unit=Resunit, file='PastSigmas.dat')
				DO i=1,ClosePanCount
					READ(unit=Resunit, fmt=21) sigma(i)
				ENDDO
				CLOSE (unit=Resunit, status='keep')
			ELSE
				sigma(:)=0
			ENDIF

			ALLOCATE(Along(ClosePanCount*ClosePanCount))
			DO i=1, ClosePanCount
				DO j=1,ClosePanCount
					panel=(i-1)*ClosePanCount+j
					Along(panel)=A(i,j)
				ENDDO
			ENDDO
			CALL  gmresIter(RHStemp, sigma, Along, ClosePanCount)
			DEALLOCATE(Along)

		CASE (2)  
    			!---- Gaussian inversion with partial pivoting
			ALLOCATE( indxLegs(ClosePanCount) ); indxLegs(:)=-1
			ALLOCATE( sigmaGauss(ClosePanCount) ); sigmaGauss(:)=0.
			ALLOCATE( sigma(ClosePanCount)); sigma(:)=-1.
			CALL LEGS (A, ClosePanCount, RHStemp, sigmaGauss, indxLegs)
			PRINT *, '    Done with inversion'

			DO j=1,ClosePanCount	! correct the shuffled sigma locations
				sigma(indxLegs(j))=sigmaGauss(j)
			ENDDO
			DEALLOCATE(indxLegs); DEALLOCATE(sigmaGauss)
			PRINT *, '  Done after full Gaussian inversion'

		CASE (3)
				!---- Gaussian inversion WITHOUT partial pivoting
			ALLOCATE( sigma(ClosePanCount)); sigma(:)=-1.
			CALL GAUSS(A, ClosePanCount, sigma, RHStemp)
			PRINT *, '  Done after regular Gaussian inversion'
		CASE DEFAULT
			STOP 'Wrong method for MatSolveMethod'
		END SELECT
		CALL cpu_time(t2); timeCount=t2-t1

		WRITE(*, fmt=16) '      sigma(1)=', sigma(1)
		WRITE(*, fmt=19) '   Solving for sigma took ', timeCount, ' seconds.'

				! write out sigma values to make future calls to GMRES faster
		IF (1==ShowSigmas) THEN
			PRINT *, '  -Writing "PastSigmas.dat" '
			OPEN(unit=Resunit, file='PastSigmas.dat', status='replace')
			DO i=1,ClosePanCount
				WRITE(unit=Resunit, fmt=21) sigma(i)+1e5
			ENDDO
			CLOSE (unit=Resunit, status='keep')
		ENDIF

		!-----------------------------
		! assigns all above values to output vector 'Pts_and_force' for the non-boundary pts
		  ! [xpt, ypt, xnormal, ynormal, radius to blob center, (0-2pi radians) theta to 
		  ! blob center, blob number, blob_x, blob_y, panel force] for each panel
		panel=1; SkipPanTo=1; MxPt_blob=0; ALLOCATE(PtsForce( ClosePanCount, 5) )

		iBlobDo: DO i=1,BlobPlus
			MxPt_blob=MAX(MxPt_blob,NewNumPtsPan(i))
			blob_xc(i)=NewCoords(SkipPanTo,6); blob_yc(i)=NewCoords(SkipPanTo,7)

			! Create temp horizontal reference line
			aa=0.1; PtC(1)=NewCoords(SkipPanTo,6); PtC(2)=NewCoords(SkipPanTo,7) 
			PtRef(1)=PtC(1)+0.1; PtRef(2)=PtC(2)

			jsizeBlobDo: DO j=1,NewNumPtsPan(i)
				IF (panel>ClosePanCount) THEN
					PRINT *, 'Shading:BFC Pts_force panel err. panel,i,j=', panel, i,j
					STOP
				ENDIF

				PtsForce(panel,1)=NewCoords(SkipPanTo,1); PtsForce(panel,2)=NewCoords(SkipPanTo,2)
				PtsForce(panel,3)=normPanX(panel);	PtsForce(panel,4)=normPanY(panel)
				PtsForce(panel,5)=i

				PtA(1)=NewCoords(SkipPanTo,1); PtA(2)=NewCoords(SkipPanTo,2)

				ref = SQRT( (PtA(1)-PtC(1))*(PtA(1)-PtC(1)) + (PtA(2)-PtC(2))*(PtA(2)-PtC(2)) )
				c = SQRT( (PtA(1)-PtRef(1))*(PtA(1)-PtRef(1)) + (PtA(2)-PtRef(2))*(PtA(2)-PtRef(2)) )

				! 0 < acos(x) < pi. Move to correct quadrant
				IF (0==ref) THEN
					ThreePtAngle=0
				ELSE
					ThreePtAngle=( (aa*aa) + (ref*ref) - (c*c) )/(2*aa*ref)
				ENDIF
				radBlob(i,j)=ref

				Flatcheck: IF ( ABS(ThreePtAngle+1)<epsil ) THEN	! too close to -1
					thetaBlob(i,j)=pi
				ELSEIF ( ABS(ThreePtAngle-1)<epsil ) THEN			! too close to 1
					thetaBlob(i,j)=0
				ELSE
					IF ( PtA(2)>=PtC(2) ) THEN
						thetaBlob(i,j)=acos( ThreePtAngle )
					ELSE
						thetaBlob(i,j)=-acos( ThreePtAngle )
					ENDIF

					IF ( thetaBlob(i,j) < 0 ) THEN
						thetaBlob(i,j)=thetaBlob(i,j)+2.0*pi
					ENDIF

				ENDIF Flatcheck

			 	! Force/m^2 = -e_0 * En^2
				forceBlob(i,j)= -epsil0*sigma(panel)*sigma(panel)

				panel=panel+1; SkipPanTo=SkipPanTo+1
			ENDDO jsizeBlobDo
			SkipPanTo=SkipPanTo+1
		ENDDO iBlobDo

		!-----------------------------
		! Smooth the edges out for the needle shape
		PRINT *, '    Forced smoothing for ',xLeftSmoothEdge,' panels.'
		DO j=1,xLeftSmoothEdge
			forceBlob(NeedleNum,j)=0
			forceBlob(NeedleNum,NewNumPtsPan(NeedleNum)-j+1)=0
		ENDDO

		! writing the data file of how many blobs existed. Used for seeing when another one breaks off
		PRINT *, '  -Writing "NumBlob.dat" '
		OPEN(unit=NumBlobunit, file='NumBlob.dat', status='old')
		READ(unit=NumBlobunit, fmt='(i4)') z
		REWIND(unit=NumBlobunit)
		WRITE(unit=NumBlobunit, fmt='(i4)') BlobNum 				
		CLOSE(unit=NumBlobunit, status='keep')

		!-----------------------------
		! Write number of blobs and size/charge in scenario
		highBlobNumIf: IF ( BlobNum>1 ) THEN	
			! Create snapoff PDF
			PRINT *, 'There are snapped off blobs! Number detached=', BlobNum-1
			ALLOCATE(PanCharge(BlobNum)); ALLOCATE(PanArea(BlobNum))
			tempRange=0; PanCharge(:)=0; PanArea(:)=0
			DO i=1,BlobNum
				IF (i>1) THEN
					tempRange = tempRange+NewNumPtsPan(i-1)
				ENDIF

				RunCharge: IF (i /=NeedleNum) THEN
					NewdropPts=NewNumPtsPan(i)
					ALLOCATE(temp1(NewdropPts)); ALLOCATE(temp4(NewdropPts,2))

					temp1=sigma(tempRange+1:tempRange+NewdropPts)
					temp4=PtsForce(tempRange+1:tempRange+NewdropPts,1:2)

					!..... charge in [C]
					CALL CalcPDF( PanCharge(i), PanLen, temp1, NewdropPts )
					!..... area in [cm^2]
					CALL polygon_area_2d_2( PanArea(i), temp4, NewdropPts )
					!.....	

					DEALLOCATE(temp1); DEALLOCATE(temp4)
				ENDIF RunCharge
			ENDDO

			PRINT *, '  -Writing PDFBlob.dat'
			OPEN(unit=NewDropunit, file='PDFBlob.dat', status='old', access='sequential', position='append')
			OPEN(unit=NewDropunit2, file='PDFBlob_comma.dat', status='old', access='sequential', position='append')
			DO i=1,BlobNum
				PrintCharge: IF (i /=NeedleNum) THEN
					! changing pi_r^2 into 4/3*pi_r^3. volume in [m^3]
					droprad = 0.01*SQRT( PanArea(i)/pi )			! radius [m]
					!dropVol = (4/3.0)*pi*r2*droprad
					dropVol = pi*droprad*droprad*1					! [m^3] now, using cylinder of charge

					chargeMassArea = PanCharge(i)/(PanArea(i)*q)	! [#charges/cm^2]
					chargeMassVol = PanCharge(i)/dropVol			! [C/m^3]
					!chi = PanCharge(i)*PanCharge(i)/(64*pi*pi*epsil0*surften_in*aa) (sqrt?)
					chi = 58780*PanCharge(i)*PanCharge(i)/SQRT(droprad)	! 3.455e9 otherwise

					!MTCR = 871080*SQRT( rho_in*rho_in*aa/(36*epsil0*surften_in) )	! 
					MTCR = dropVol*rho_in*5.068e24*q/PanCharge(i)	! [#in particle/1 short e-]
					WRITE(unit=NewDropunit, fmt=11) PanCharge(i), PanArea(i), dropVol, &
						chargeMassArea, chargeMassVol, blob_xc(i), blob_yc(i), chi, MTCR, i, &
						counter-1, totalTime
					WRITE(unit=NewDropunit2, fmt=20) PanCharge(i),',', PanArea(i),',', dropVol,',', &
						chargeMassArea,',', chargeMassVol,',', blob_xc(i),',', blob_yc(i),',', &
						chi,',', MTCR,',', i,',', counter-1,',', totalTime

					! mark each blob # once as it passes about a global height
					IF (highPDFpt<blob_yc(i) .AND. i>highblobNum) THEN
						PRINT *, '  -Writing "PDFhigh.dat"'
						highblobNum=highblobNum+1
						OPEN(unit=PDFhighunit, file='PDFhigh.dat', status='old', access='sequential', position='append')
						WRITE(unit=PDFhighunit, fmt=20) PanCharge(i),',', PanArea(i),',', dropVol,',', &
							chargeMassArea,',', chargeMassVol,',', blob_xc(i),',', blob_yc(i),',', &
							chi,',', MTCR,',', i,',', counter-1,',', totalTime
						CLOSE (unit=PDFhighunit, status='keep')
					ENDIF
				ENDIF PrintCharge
			ENDDO
			CLOSE (unit=NewDropunit, status='keep')
			CLOSE (unit=NewDropunit2, status='keep')

			IF (z<BlobNum) THEN
				PRINT *, '  -Writing "PDFsingle.dat"'
				DO i=z+1,BlobNum
					IF (i/=NeedleNum) THEN
						OPEN(unit=PDFsingleunit, file='PDFsingle.dat', status='old', access='sequential', position='append')
						WRITE(unit=PDFsingleunit, fmt=20) PanCharge(i),',', PanArea(i),',', dropVol,',', &
							chargeMassArea,',', chargeMassVol,',', blob_xc(i),',', blob_yc(i),',', &
							chi,',', MTCR,',', i,',', counter-1,',', totalTime
						CLOSE (unit=PDFsingleunit, status='keep')
					ENDIF
				ENDDO
			ENDIF
			DEALLOCATE(PanCharge); DEALLOCATE(PanArea)
		ENDIF highBlobNumIf

		!-----------------------------		
		! run optional arbitrary point potential analysis to see plots
		OptionalCalcPot: IF (1==ShowPotentials) THEN
			CALL cpu_time(t1)
			ALLOCATE(ArbPot(PotMeshSize,PotMeshSize))
			ALLOCATE(InsideYes(PotMeshSize,PotMeshSize)); InsideYes(:,:)=-1 
			ALLOCATE(meshX(PotMeshSize)); meshX(:)=-1.
			ALLOCATE(meshY(PotMeshSize)); meshY(:)=-1.

			! calculates whether something is inside or outside a shape
			CALL CalcPot(InsideYes, meshX, meshY, PotMeshSize, ClosePanCount, xPanBegin, yPanBegin)

			! -- debugging line, to see how potentials inside droplets behave
			InsideYes(:,:)=1; zero=0.0
			!-------

			PotDo: DO i=1,PotMeshSize
				Jdo: DO j=1,PotMeshSize
					ArbPot(i,j)=0

					IF (1==InsideYes(i,j)) THEN
						Kdo: DO k=1,ClosePanCount
							Zdo: DO z=1,M
								panel=(k-1)*M + z
								delX=meshX(i)-xPan(panel)
								delY=meshY(j)-yPan(panel)

								CALL TrapInit ( GreenPart, delX, delY, normPanX(k), normPanY(k), &
									zero, zero, gwgt(z), PanLen(k), 11, 0, k, M)

								CALL TrapInit ( GreenPart2, delX, delY, normPanX(k), normPanY(k), &
									zero, zero, gwgt(z), PanLen(k), 10, 0, k, M)

								IF (1==RHS(k,2)) THEN
									ArbPot(i,j)=ArbPot(i,j) + RHS(k,1)*GreenPart - sigma(k)*GreenPart2
								ELSEIF (0==RHS(k,2)) THEN
									ArbPot(i,j)=ArbPot(i,j) + sigma(k)*GreenPart - RHS(k,1)*GreenPart2
								ELSE
									PRINT *, 'Shading:BFC Arbpot err, k, RHS(k)=', k, RHS(k,2)
									STOP 
								ENDIF
							ENDDO Zdo
						ENDDO Kdo
					ELSEIF (0==InsideYes(i,j)) THEN
						ArbPot(i,j)=-2;
					ELSE
						PRINT *, 'MAT_Greens: InsideYes err'
						STOP
					ENDIF
				ENDDO Jdo
			ENDDO PotDo

			PRINT *, '  -Writing "ArbPot.dat" '
			OPEN(unit=Potunit, file='ArbPot.dat', status='replace')
			WRITE(unit=Potunit, fmt='(a,i5,a)') 'TITLE="Potentials using boundary matrix with ', totPanels, ' panels"'
			WRITE(unit=Potunit, fmt='(a)') 'VARIABLES="X [cm]"'
			WRITE(unit=Potunit, fmt='(a)') '"Y [cm]"'
			WRITE(unit=Potunit, fmt='(a)') '"Potential [V]"'
			WRITE(unit=Potunit, fmt='(a)') 'ZONE T="main"'
			WRITE(unit=Potunit, fmt='(a,i3,a,i3,a)') 'I=', PotMeshSize, ', J=', PotMeshSize, ', K=1, F=POINT'

			DO i=1,PotMeshSize
				DO j=1,PotMeshSize
					WRITE(unit=Potunit, fmt=13) meshX(i), meshY(j), ArbPot(i,j)
				ENDDO
			ENDDO
			CLOSE (unit=Potunit, status='keep')

			CALL cpu_time(t2); timeCount=t2-t1
			WRITE(*, fmt=19) '     Solving for potentials took ', timeCount, ' seconds.'

			DEALLOCATE(ArbPot); DEALLOCATE(InsideYes)
			DEALLOCATE(meshX); DEALLOCATE(meshY); 
		ENDIF OptionalCalcPot

		!-----------------------------
		! write out last data file
		IF (1==ShowPt_Force) THEN
			! convert data number to string file name
			CALL Int2Str(LSstep, 'PForce', counter-1,'.dat',4)
			PRINT *, '  -Writing ', LSstep
			OPEN(unit=Forceunit, file=LSstep, status='replace')
			WRITE(unit=Forceunit, fmt='(a)') 'VARIABLES="X [cm]"' !#1
			WRITE(unit=Forceunit, fmt='(a)') '"Y [cm]"'		!#2
			WRITE(unit=Forceunit, fmt='(a)') '"Xnorm"'		!#3
			WRITE(unit=Forceunit, fmt='(a)') '"Ynorm"'		!#4
			WRITE(unit=Forceunit, fmt='(a)') '"Radius"'		!#5
			WRITE(unit=Forceunit, fmt='(a)') '"Theta"'		!#6
			WRITE(unit=Forceunit, fmt='(a)') '"Blob #"'		!#7
			WRITE(unit=Forceunit, fmt='(a)') '"Blob xc"' 	!#8
			WRITE(unit=Forceunit, fmt='(a)') '"Blob yc"'	!#9	
			WRITE(unit=Forceunit, fmt='(a)') '"Force"'		!#10
			WRITE(unit=Forceunit, fmt='(a)') '"Pan #"'		
			WRITE(unit=Forceunit, fmt='(a)') 'ZONE T="main"'

			k=0
			DO i = 1,BlobPlus
				DO j=1,NewNumPtsPan(i)
					k=k+1
					IF (PtsForce(k,2) /= -1) THEN
						WRITE(unit=Forceunit, fmt=14) PtsForce(k,1), PtsForce(k,2), &
							PtsForce(k,3), PtsForce(k,4), radBlob(i,j), &
							thetaBlob(i,j), INT(PtsForce(k,5)), blob_xc(i), &
							blob_yc(i), forceBlob(i,j), k
					ENDIF
				ENDDO
			ENDDO
			CLOSE (unit=Forceunit, status='keep')

			PRINT *,'    BlobNum, MxPt_blob=',BlobNum, MxPt_blob
		ENDIF

		!-----------------------------		
		! deallocates everything - you're done
		DEALLOCATE(A)
		DEALLOCATE(xPanLen); DEALLOCATE(yPanLen)
		DEALLOCATE(xPanBegin); DEALLOCATE(yPanBegin)
		DEALLOCATE(PanLen); DEALLOCATE(xPan)
		DEALLOCATE(yPan); DEALLOCATE(normPanX)
		DEALLOCATE(normPanY); DEALLOCATE(RHS)
		DEALLOCATE(RHStemp); DEALLOCATE(xPanCen)
		DEALLOCATE(yPanCen); DEALLOCATE(RHS_gamma_phi)
		DEALLOCATE(RHS_gamma_dphi); DEALLOCATE(RHS_alpha_phi)
		DEALLOCATE(RHS_alpha_dphi); DEALLOCATE(Coords)
		DEALLOCATE(NumPanels); DEALLOCATE(sigma) 
		DEALLOCATE(PtsForce); DEALLOCATE(phi)

		PRINT *, 'Done with main routine <BoundaryForceCalc>!'

		GOTO 999		! actual program end, just lists formats and errors below

		!-----------------------------
		! FORMAT listings
		11 FORMAT(5e12.4, 3e13.5, e12.4, i4, i5, e12.4)
		12 FORMAT(2i5, 3f20.10)
		13 FORMAT(3e14.6)
		14 FORMAT(6e13.5, i3, 3e13.5, i5)
		16 FORMAT(a,e11.4)
		17 FORMAT(a,2e11.4)
		18 FORMAT(a,i4,a)
		19 FORMAT(a,f8.2,a)
		20 FORMAT(e12.4, a, e12.4, a, e12.4, a, e12.4, a, e12.4, a, e13.5, a, e13.5, &
			a, e12.4, a, e12.4, a, i4, a, i5, a, e12.4)
		21 FORMAT(e14.6)

		999 CALL cpu_time(t2); timeCount=t2-t3
		WRITE(*, fmt=19) '   Entire UMich routine took ', timeCount, ' seconds.'

		IF (counter>MaxStep) THEN
			PRINT *, 'Ending Routine at step = ',MaxStep
			STOP
		ELSE
			PRINT *, '********************'
		ENDIF
	END SUBROUTINE BoundaryForceCalc

! ***************
! Subroutine ShapeGrow
! ***************
! DESCRIPTION: for the Sussman level set visual file, finds and sorts the number and vertex members 
!   of each blob
!
! INPUT: visual data
!
! OUTPUT: members of each shape, number of shapes, ordered surface around
	SUBROUTINE ShapeGrow( SortedBlobPts, Countblob, NeedleBlob, MxPt_blob, sizeNeedleBlob, NumPt, &
			LSMaxX, LSMaxY, xLeft, yBottom, delX, delY, is_symmetric1, minX ) 
		! incoming variables
		INTEGER, INTENT(in)			:: LSMaxX, LSMaxY, is_symmetric1
		REAL*8, INTENT(in)			:: xLeft, yBottom, delX, delY

		! in/outgoing variables
		INTEGER, INTENT(inout)		:: NumPt, Countblob, NeedleBlob, MxPt_blob
		INTEGER, DIMENSION(SmallInitAry),INTENT(inout) 			:: sizeNeedleBlob
		REAL*8, INTENT(inout)		:: minX
		REAL*8, DIMENSION(InitArySize,8), INTENT(inout)	:: SortedBlobPts
	      ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatyesno]

		! program variables
		INTEGER						:: i, j, k, aa, bb, cc, tempSize, CountUp, &
			Crossunit, Sortunit, startCase, prevNumPt, tempMaxSize, tempPrevStart, &
			Countblobbump
		INTEGER, DIMENSION(:), ALLOCATABLE	:: InPts, temp3
		REAL*8, DIMENSION(InitArySize,4)  	:: tempInPts, swapInPts
		  ! [xpt, ypt, blob#, particle#]
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: temp1, temp2

		CONTINUE

		! variable setup............
		Crossunit=12; Sortunit=13

		! Countblob=number of discrete blobs. NumPt=total number of points in all blobs
		Countblob=0; NumPt=0; tempInPts(:,:)=0
		Countblobbump=0

		! begins computation............
		! groups each point into a blob family
		PRINT *, '  Entering surface determination'
		tempInPts(:,:)=-1; minX=1e6
		DO i=StartNum,LSMaxX+StopNum
			DO j=LSMaxY+StopNum,StartNum,-1
				insidePhiIf: IF ( phi(i,j) >= 0 .AND. -1==insideMarker(i,j) ) THEN
					Countblob=Countblob+1; prevNumPt=NumPt
					startCase=8

					!.....
					CALL BlobSeparate( LSMaxX, LSMaxY, Countblob, i, j )
					!.....
					CALL SurfaceTrack (tempInPts, NumPt, Countblob, i, j, &
						startCase, delX, delY, xLeft, yBottom, LSMaxX, LSMaxY)
					!.....	

					! removes repeat Coordinate points
					RmCoordDo: DO aa=prevNumPt+1,NumPt-1
						Do bb=aa+1,NumPt
							IF ( tempInPts(aa,1)==tempInPts(bb,1) .AND. &
								tempInPts(aa,2)==tempInPts(bb,2) ) THEN
								IF (bb<NumPt) THEN 
									DO cc=bb+1,NumPt
										tempInPts(cc-1,:)=tempInPts(cc,:)
									ENDDO
								ELSE
									tempInPts(NumPt,:)=tempInPts(NumPt,:)
								ENDIF
								NumPt=NumPt-1
								PRINT *, '    Remove Coord pt; aa,bb=',aa,bb
							ENDIF
						ENDDO
					ENDDO RmCoordDo

					! create mirror image of all blobs if symmetric
					symDo: IF ( 1 == is_symmetric1 ) THEN
						tempSize=NumPt-prevNumPt

						EnoughPtsIf: IF (tempSize>1) THEN
							IF (Countblob>1) THEN
								PRINT *, '  Multiple blobs. Countblob=', Countblob
							ENDIF

							DO aa=prevNumPt+1,NumPt	! find x minimum
								IF ( tempInPts(aa,1)<minX ) THEN
									minX=tempInPts(aa,1)
								ENDIF
							ENDDO
							swapInPts(prevNumPt+1:NumPt,:)=tempInPts(prevNumPt+1:NumPt,:); k=0
							tempMaxSize=tempSize*2+prevNumPt; tempPrevStart=prevNumPt+1
							DO aa=tempPrevStart,tempMaxSize	! now add mirrored points, {x,y, blob#, particle#}
								IF (aa<=tempSize+prevNumPt) THEN
									k=k+1
									tempInPts(aa,1)=2*minX-swapInPts(NumPt-k+1,1)-1e-5
									tempInPts(aa,2)=swapInPts(NumPt-k+1,2)
									tempInPts(aa,3)=swapInPts(NumPt-k+1,3)+Countblobbump
								ELSE
									IF (i/=1 .AND. tempSize+prevNumPt+1==aa) THEN
										Countblobbump=Countblobbump+1
									ENDIF
									tempInPts(aa,1)=swapInPts(aa-NumPt+prevNumPt,1)
									tempInPts(aa,2)=swapInPts(aa-NumPt+prevNumPt,2)
									tempInPts(aa,3)=swapInPts(aa-NumPt+prevNumPt,3)+Countblobbump
								ENDIF
								tempInPts(aa,4)=aa
							ENDDO
							NumPt=NumPt+tempSize
							PRINT *, '  Done flipping symmetric shape. NumPt=',NumPt
						ELSE
							Countblob=Countblob-1
							NumPt=prevNumPt
						ENDIF EnoughPtsIf
					ELSE
						PRINT *, '  NumPt=',NumPt
					ENDIF symDo
				ENDIF insidePhiIf 
			ENDDO
		ENDDO
		PRINT *, '  After SurfaceTrack, there are ',NumPt,' coordinates.'
		Countblob=Countblob+Countblobbump

		! =begin blob1, begin blob2
		ALLOCATE(InPts(Countblob))
		InPts(:)=-1
		CountUp=1
		DO i=1,NumPt
			IF ( CountUp==INT(tempInPts(i,3)) ) THEN
				InPts(CountUp)=i
				CountUp=CountUp+1
			ENDIF
		ENDDO
		CountUp=CountUp-1
		PRINT *, '  Done grouping points into blobs. NumBlob=', CountUp

		! ***************
		IF (1==ShowNewCross) THEN
			PRINT *, '  -Writing "BlobCross.dat" '
			OPEN(unit=Crossunit, file='BlobCross.dat', status='replace', err=990)
			WRITE(unit=Crossunit, fmt='(a)', err=992) 'VARIABLES="X [cm]"'
			WRITE(unit=Crossunit, fmt='(a)', err=992) '"Y [cm]"'
			WRITE(unit=Crossunit, fmt='(a)', err=992) '"BlobNum"'
			WRITE(unit=Crossunit, fmt='(a)', err=992) '"ParticleNum"'
			DO i=1,NumPt
				WRITE(unit=Crossunit, fmt='(2e18.7, 2i8)') tempInPts(i,1), tempInPts(i,2), &
					INT(tempInPts(i,3)), INT(tempInPts(i,4))
			ENDDO
			CLOSE (unit=Crossunit, status='keep', err=993)
		ENDIF

		ALLOCATE(temp1(NumPt,8)); ALLOCATE(temp2(NumPt,4)); ALLOCATE(temp3(CountBlob))

		temp1=SortedBlobPts(1:NumPt,:)
		temp2=tempInPts(1:NumPt,:)
		temp3=sizeNeedleBlob(1:CountBlob)
		!.....
		CALL UpdateBlobList( temp1, NeedleBlob, MxPt_blob, temp3, InPts, temp2, Countblob, NumPt )
		!.....

		SortedBlobPts(1:NumPt,:)=temp1
		sizeNeedleBlob(1:CountBlob)=temp3;

		DEALLOCATE(temp1); DEALLOCATE(temp2); DEALLOCATE(temp3)

		! ***************
		IF (1==ShowSort) THEN
			PRINT *, '  -Writing "BlobSort.dat" '
			OPEN(unit=Sortunit, file='BlobSort.dat', status='replace', err=990)
			WRITE(unit=Sortunit, fmt='(a)', err=992) 'VARIABLES="X [cm]"'
			WRITE(unit=Sortunit, fmt='(a)', err=992) '"Y [cm]"'
			WRITE(unit=Sortunit, fmt='(a)', err=992) '"BlobNum"'
			WRITE(unit=Sortunit, fmt='(a)', err=992) '"Blob_xc"'
			WRITE(unit=Sortunit, fmt='(a)', err=992) '"Blob_yc"'
			DO i=1,NumPt
				WRITE(unit=Sortunit, fmt='(2e14.6, i8, 2e14.6)') SortedBlobPts(i,1), SortedBlobPts(i,2), &
					INT(SortedBlobPts(i,5)), SortedBlobPts(i,6), SortedBlobPts(i,7)
			ENDDO
			CLOSE (unit=Sortunit, status='keep', err=993)
		ENDIF

		GOTO 998

	    ! ERROR listings
		990 print *, 'Error open blobarray'
		992 print *, 'Error write blobarray'
		993 print *, 'Error close blobarray'

		998 DEALLOCATE(InPts)

		PRINT *, 'Done with "ShapeGrow" in <Shading>'

	END SUBROUTINE ShapeGrow

! ***************
! SUBROUTINE SurfaceTrack
! ***************
!
! DESCRIPTION: tracks the blob from the left edge around the surface. Promises to simplify later
!   sorting tasks
!
! INPUT: ls, ls blob number. first point on surface
!
! OUTPUT: sorted members in each blob
!	
! APPROACH: move clockwise around from straight up. If phi(new)=negative, mark interpolated 
!   crossing point and continue looking around clockwise. If =positive, don't mark, but re-look
!   around
!
! CALLING ROUTINE: Shading:ShapeGrow

	RECURSIVE SUBROUTINE SurfaceTrack ( tempIn, NumPt, NumBlob, WhichI, &
			WhichJ, prevCase, delX, delY, xLeft, yBottom, LSMaxX, LSMaxY )
		! incoming variables
		INTEGER, INTENT(in)			:: WhichI, WhichJ, NumBlob, LSMaxX, LSMaxY
		REAL*8, INTENT(in)			:: delX, delY, xLeft, yBottom

		! outgoing variable
		INTEGER, INTENT(inout)		:: NumPt, prevCase
		REAL*8, DIMENSION(InitArySize,4), INTENT(inout) 	:: tempIn
		  ! {x,y, blob#, particle#}

		! internal variables
		LOGICAL						:: onEdge, negPhi, diagDir, LoopChoke
		INTEGER						:: i, iNext, jNext, caseplus2i, caseplus2j, &
			countLoop, plusPt, minusPt
		CONTINUE

		onEdge=.FALSE.; alreadyChecked(WhichI, WhichJ)=1
		! check to make sure at least one of the 8points is phi<0 (i.e. you are on the boundary)
		DO i=1,8
			CALL CaseLoop ( iNext, jNext, caseplus2i, caseplus2j, i, &
				WhichI, WhichJ, LSMaxX, LSMaxY )

			IF (phi(iNext, jNext)<=0) THEN
				onEdge=.TRUE.
				EXIT
			ENDIF
		ENDDO
		IF (.NOT.(onEdge)) THEN	! not on edge
			PRINT *, 'Shading:SurfaceTrack end of shape at i,j,prevCase,NumPt=', &
				WhichI, WhichJ,prevCase,NumPt
		ENDIF

		negPhi=.TRUE.; countLoop=0
		negLoop: DO WHILE (negPhi .AND. onEdge)
			countLoop=countLoop+1
			jMax: IF (iNext<LSMaxX .AND. jNext<LSMaxY) THEN
				IF (1==mod(prevCase,2)) THEN
					diagDir=.TRUE.
				ELSE
					diagDir=.FALSE.
				ENDIF
				CALL CaseLoop ( iNext, jNext, caseplus2i, caseplus2j, prevCase+1, &
					WhichI, WhichJ, LSMaxX, LSMaxY )
				MoveMark: IF ( phi(iNext, jNext)<0 .AND. .NOT. (diagDir) .AND. &
					phi(WhichI,WhichJ)/=0 ) THEN	! * MARK *
					NumPt=NumPt+1

					IF ( 9==prevCase+1 .OR. 5==prevCase+1 ) THEN		! if up/down
						tempIn(NumPt,1)=(WhichI-1)*delX+xLeft
						tempIn(NumPt,2)=LinFrac(phi(WhichI,WhichJ),phi(iNext, jNext), &
							(WhichJ-1)*delY+yBottom, (jNext-1)*delY+yBottom)

						! check to see if you are on a horizontal choke point
						IF (WhichJ<LSMaxY) THEN
							plusPt=WhichJ+1
						ELSE
							plusPt=WhichJ-1
						ENDIF
						IF (WhichJ/=1) THEN
							minusPt=WhichJ-1
						ELSE
							minusPt=WhichJ+1
						ENDIF
						IF (phi(WhichI, plusPt)<0 .AND. phi(WhichI, minusPt)<0) THEN
							LoopChoke=.TRUE.; i=0
							DO WHILE (LoopChoke)
								i=i+1
								IF (WhichI==ChokeCoord(i,1) .AND. WhichJ==ChokeCoord(i,2)) THEN
									ChokeCoord(i,3)=ChokeCoord(i,3)+1
									LoopChoke=.FALSE.
								ELSEIF (0==ChokeCoord(i,1)) THEN
									ChokeCoord(i,1)=WhichI; ChokeCoord(i,2)=WhichJ
									ChokeCoord(i,3)=ChokeCoord(i,3)+1
									alreadyChecked(WhichI,WhichJ)=0	! allow one come-back node landing
									LoopChoke=.FALSE.
								ENDIF
							ENDDO
						ENDIF
					ELSEIF ( 3==prevCase+1 .OR. 7==prevCase+1 ) THEN	! if left/right
						tempIn(NumPt,1)=LinFrac(phi(WhichI,WhichJ),phi(iNext, jNext), &
							(WhichI-1)*delX+xLeft, (iNext-1)*delX+xLeft)
						tempIn(NumPt,2)=(WhichJ-1)*delY+yBottom

						! check to see if you are on a vertical choke point
						IF (WhichI<LSMaxX) THEN
							plusPt=WhichI+1
						ELSE
							plusPt=WhichI-1
						ENDIF
						IF (WhichI/=1) THEN
							minusPt=WhichI-1
						ELSE
							minusPt=WhichI+1
						ENDIF
						IF (phi(plusPt, WhichJ)<0 .AND. phi(minusPt, WhichJ)<0) THEN
							LoopChoke=.TRUE.; i=0
							DO WHILE (LoopChoke)
								i=i+1
								IF (WhichI==ChokeCoord(i,1) .AND. WhichJ==ChokeCoord(i,2)) THEN
									ChokeCoord(i,3)=ChokeCoord(i,3)+1
									LoopChoke=.FALSE.
								ELSEIF (0==ChokeCoord(i,1)) THEN
									ChokeCoord(i,1)=WhichI; ChokeCoord(i,2)=WhichJ
									ChokeCoord(i,3)=ChokeCoord(i,3)+1
									alreadyChecked(WhichI,WhichJ)=0	! allow one come-back node landing
									LoopChoke=.FALSE.
								ENDIF
							ENDDO
						ENDIF
					ELSE
						PRINT *, 'Shading:SurfaceTrack mark location err.'
						STOP
					ENDIF
					tempIn(NumPt,3)=1.0*NumBlob
					tempIn(NumPt,4)=1.0*NumPt

					prevCase = CaseNumShift(prevCase,1)
				ELSEIF (0==phi(iNext,jNext) .AND. alreadyChecked(iNext,jNext) /=1) THEN		! catching unusual case where no marking would occur
					NumPt=NumPt+1

					tempIn(NumPt,1)=(iNext-1)*delX+xLeft
					tempIn(NumPt,2)=(jNext-1)*delY+yBottom
					tempIn(NumPt,3)=1.0*NumBlob
					tempIn(NumPt,4)=1.0*NumPt

					prevCase = CaseNumShift(prevCase,-3)
					CALL SurfaceTrack ( tempIn, NumPt, NumBlob, iNext, jNext, prevCase, &
						delX, delY, xLeft, yBottom, LSMaxX, LSMaxY )
				ELSEIF ( (phi(iNext, jNext)>=0 .AND. alreadyChecked(iNext,jNext) /=1) .AND. &
					((phi(caseplus2i, caseplus2j)>=0 .AND. (diagDir)) .OR. .NOT. (diagDir)) ) THEN	! * MOVE *
					negPhi=.FALSE.
					prevCase = CaseNumShift(prevCase,-3)
					CALL SurfaceTrack ( tempIn, NumPt, NumBlob, iNext, jNext, prevCase, &
						delX, delY, xLeft, yBottom, LSMaxX, LSMaxY )
				ELSE
					prevCase = CaseNumShift(prevCase,1)
				ENDIF MoveMark
			ELSE
				onEdge=.FALSE.
			ENDIF jMax
			IF (countLoop>8) THEN	! you're at the end of the shape
				EXIT
			ENDIF
		ENDDO negLoop

		! Error checking
		IF (NumPt-4>InitArySize) THEN
			PRINT *, 'Shading:BlobSeparate InitArySize overflow'
			STOP
		ENDIF
	END SUBROUTINE SurfaceTrack

! ***************
! SUBROUTINE CaseLoop
! ***************
!
! DESCRIPTION: Defines before and after points clockwise around the grid
!
! INPUT: case point
!
! OUTPUT: next case point
!
! CALLING PROGRAM: SurfaceTrack
	SUBROUTINE CaseLoop ( iNext, jNext, caseplus2i, caseplus2j, caseFrom, &
			iPt, jPt, xMax, yMax )
		! incoming variables
		INTEGER, INTENT(in)		:: iPt, jPt, caseFrom, xMax, yMax

		! outgoing variables
		INTEGER, INTENT(out)	:: iNext, jNext, caseplus2i, caseplus2j

		! program variables
		CONTINUE

		SELECT CASE (caseFrom)	
		CASE (1)	! up
			iNext=iPt; jNext=jPt+1
			caseplus2i=iPt+1; caseplus2j=jPt+1
		CASE (2)	! up-right
			iNext=iPt+1; jNext=jPt+1
			caseplus2i=iPt+1; caseplus2j=jPt
		CASE (3)	! right
			iNext=iPt+1; jNext=jPt
			caseplus2i=iPt+1; caseplus2j=jPt-1
		CASE (4)	! down-right
			iNext=iPt+1; jNext=jPt-1
			caseplus2i=iPt; caseplus2j=jPt-1
		CASE (5)	! down
			iNext=iPt; jNext=jPt-1
			caseplus2i=iPt-1; caseplus2j=jPt-1
		CASE (6)	! down-left
			iNext=iPt-1; jNext=jPt-1
			caseplus2i=iPt-1; caseplus2j=jPt
		CASE (7)	! left
			iNext=iPt-1; jNext=jPt
			caseplus2i=iPt; caseplus2j=jPt+1
		CASE (8)	! left-up
			iNext=iPt-1; jNext=jPt+1
			caseplus2i=iPt; caseplus2j=jPt+1
		CASE (9)	! wrap-around, repeated up
			iNext=iPt; jNext=jPt+1
			caseplus2i=iPt+1; caseplus2j=jPt+1
		CASE (10)	! wrap-around, repeated up-right
			iNext=iPt+1; jNext=jPt+1
			caseplus2i=iPt+1; caseplus2j=jPt
		CASE DEFAULT 
			PRINT *, 'Shading:CaseLoop _CaseFrom_ err. iPt,jPt,iNext,jNext=', &
				iPt, jPt, iNext, jNext
			PRINT *, 'xMax, yMax, caseFrom=',xMax, yMax, caseFrom
			STOP
		END SELECT

		! checking too small boundary values
		IF (caseplus2i<1) THEN
			caseplus2i=1
		ENDIF
		IF (caseplus2j<1) THEN
			caseplus2j=1
		ENDIF
		IF (caseplus2i>xMax) THEN
			caseplus2i=xMax
		ENDIF
		IF (caseplus2j>yMax) THEN
			caseplus2j=yMax
		ENDIF
		IF (iNext<1) THEN
			iNext=1
		ENDIF
		IF (jNext<1) THEN
			jNext=1
		ENDIF
		IF (iNext>xMax) THEN
			iNext=xMax
		ENDIF
		IF (jNext>yMax) THEN
			jNext=yMax
		ENDIF
	END SUBROUTINE CaseLoop

! ***************
! FUNCTION CaseNumShift
! ***************
!
! DESCRIPTION: Moves the case number forward or backward [x] points. 1=up, goes clockwise ->8
!
! INPUT: case number, number forward/back, if opposite
!
! OUTPUT: modified case number 

	FUNCTION CaseNumShift (CaseNum,ShiftAmt)
		! incoming variables
		INTEGER, INTENT(in) :: CaseNum,ShiftAmt

		! outgoing function
		INTEGER				:: CaseNumShift

		! internal variables
		INTEGER				:: tempNum
		CONTINUE

		tempNum = CaseNum+ShiftAmt

		IF (tempNum<1) THEN
			tempNum=8+tempNum
		ELSEIF (tempNum>8) THEN
			tempNum=MOD(tempNum,8)
		ENDIF

		CaseNumShift=tempNum

		IF (CaseNumShift >8 .OR. CaseNumShift<1) THEN
			PRINT *, 'Shading:CaseNumShift size err'
			STOP
		ENDIF
	END FUNCTION CaseNumShift

! ***************
! SUBROUTINE UpdateBlobList
! ***************
!
! DESCRIPTION: Sorts points in each blob clockwise. Finds circle center, sorts on decreasing
!   angles from an arbitrary horizontal point
!
! INPUT: array of n {x,y} pts in clump N, example list of clump pts z{3:3+n}
!
! OUTPUT: sorted clumped array

	RECURSIVE SUBROUTINE UpdateBlobList ( SortedBlobPts, NeedleBlob1, MxPt_blob, sizeBlob, InPts, &
			tempInPts, Countblob, NumPts )
		! incoming variables
		INTEGER, INTENT(in)								:: Countblob, NumPts
		INTEGER, DIMENSION(Countblob), INTENT(in)		:: InPts
		REAL*8, DIMENSION(NumPts,4), INTENT(in) 	:: tempInPts

		! outgoing variable
		INTEGER, INTENT(inout)							:: NeedleBlob1, MxPt_blob
		INTEGER, DIMENSION(Countblob), INTENT(inout)	:: sizeBlob
		REAL*8,DIMENSION(NumPts,8),INTENT(inout) :: SortedBlobPts
  	      ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, type panel]

		! program variables
		INTEGER											:: i, k, CurrentPt, ptSum, iPlus
		REAL*8, DIMENSION(2)						:: PtC
		REAL*8, DIMENSION(:,:), ALLOCATABLE		:: tempArray2, tempSorted 

		CONTINUE

		ptSum=0
		iDo: DO i=1,CountBlob
			IF (i<Countblob) THEN
				sizeBlob(i)=InPts(i+1)-InPts(i)
			ELSE
				sizeBlob(i)=NumPts-InPts(i)+1
			ENDIF
			IF (i>1) THEN
				ptSum=ptSum+sizeBlob(i-1)
			ENDIF
			PtC(1)=SUM( tempInPts(InPts(i):InPts(i)+sizeBlob(i)-1, 1) )
			PtC(1)=PtC(1)/sizeBlob(i)

			PtC(2)=SUM( tempInPts(InPts(i):InPts(i)+sizeBlob(i)-1, 2) )
			PtC(2)=PtC(2)/sizeBlob(i)

			DO k=1,sizeBlob(i)
			! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, open blob]
				CurrentPt=InPts(i)+k-1
				SortedBlobPts(CurrentPt,1)=tempInPts(k+ptSum,1)
				SortedBlobPts(CurrentPt,2)=tempInPts(k+ptSum,2)
				SortedBlobPts(CurrentPt,3)=0
				SortedBlobPts(CurrentPt,4)=0
				SortedBlobPts(CurrentPt,5)=tempInPts(k+ptSum,3)
				SortedBlobPts(CurrentPt,6)=PtC(1)
				SortedBlobPts(CurrentPt,7)=PtC(2)

				IF ( sizeBlob(i)<0.3*NumPts ) THEN
					SortedBlobPts(CurrentPt,8)=0
				ELSE
					SortedBlobPts(CurrentPt,8)=1
					NeedleBlob1=i
				ENDIF
			ENDDO
		ENDDO iDo

		BigblobFind: DO i=1,Countblob
			iPlus=InPts(i)+sizeBlob(i)-1

			ALLOCATE(tempArray2(sizeBlob(i),8)); 
			ALLOCATE(tempSorted(sizeBlob(i),8)); tempSorted(:,:)=0
			tempArray2=SortedBlobPts( InPts(i):iPlus, : )

			! changes tempSorted if scattered points
			IF (6==flagNN) THEN
				! do nothing; previously sorted
				tempSorted=tempArray2
			ELSEIF (5==flagNN) THEN
				CALL SortWeighted_AngDist(tempSorted, sizeBlob(i), tempArray2)
			ELSEIF (4==flagNN) THEN
				CALL SortNN_NoLong(tempSorted, sizeBlob(i), tempArray2)
			ELSEIF (3==flagNN) THEN
				CALL SortNN_NoCross(tempSorted, sizeBlob(i), tempArray2)
			ELSEIF (2==flagNN) THEN
				CALL SortNN_Norms(tempSorted, sizeBlob(i), tempArray2)
			ELSEIF (1==flagNN) THEN
				CALL SortNN(tempSorted, sizeBlob(i), tempArray2)
			ELSEIF (0==flagNN) THEN
				PRINT *, 'Shading:UBL flagNN err. Currently cannot sort by theta'
				STOP
				!ALLOCATE(tempArray3(sizeBlob(i),3))
				!tempArray3=AngleConnect( InPts(i):iPlus, : )
				!CALL CheckCloseAngle(tempSorted, tempArray2, tempArray3, sizeBlob(i) )
				!DEALLOCATE(tempArray3)
			ELSE
				PRINT *, 'Shading:UpdateBlobList flagNN err'
			ENDIF

			SortedBlobPts( InPts(i):iPlus, : ) = tempSorted
			DEALLOCATE(tempArray2); DEALLOCATE(tempSorted)
		ENDDO BigblobFind

		IF (6==flagNN) THEN
			PRINT *, '  No additional sorting routine called.'
		ELSEIF (5==flagNN) THEN
			PRINT *, '  Done with sorting subroutines <SortWeighted_AngDist>.'
		ELSEIF (4==flagNN) THEN
			PRINT *, '  Done with sorting subroutines <SortNN_NoLong>.'
		ELSEIF (3==flagNN) THEN
			PRINT *, '  Done with sorting subroutines <SortNN_NoCross>.'
		ELSEIF (2==flagNN) THEN
			PRINT *, '  Done with sorting subroutines <SortNN_Norms>.'
		ELSEIF (1==flagNN) THEN
			PRINT *, '  Done with sorting subroutines <SortNN>.'
		ELSE
			PRINT *, '  Shading:UpdateBlobList flagNN err.'
			STOP
		ENDIF
	END SUBROUTINE UpdateBlobList

! ***************
! SUBROUTINE CheckCloseAngle
! ***************
!
! DESCRIPTION: checks to make sure that sorted angle points don't bounce back and forth 
!   too much due to shape curving and overlap
!
! INPUT: sorted points, length from center and theta from arbitray horizontal line
!
! OUTPUT: possibly updated sorted list
	SUBROUTINE CheckCloseAngle (NewSort, OldSort, ASortData, sBlob)
		! incoming variables
		INTEGER, INTENT(in)			:: sBlob
		REAL*8, DIMENSION(sBlob,3), INTENT(in) :: ASortData
		  ! [r, theta, blob#]
		REAL*8, DIMENSION(sBlob,4), INTENT(in) :: OldSort
		  ! [xpt, ypt, blob#, particle#]

		! outgoing variable
		REAL*8, DIMENSION(sBlob,4), INTENT(out) :: NewSort
		  ! [xpt, ypt, blob#, particle#]

		! program variables
		INTEGER				:: i, kup, kdown, iPlus, tempUp, tempDown, flagNoSwitch
		REAL*8				:: closeAngle, loopclose, delAngle, &
			delR, relcloseR

		CONTINUE
		closeAngle=0.05; loopclose=2.0*closeAngle
		relcloseR=0.05	! sorted points ok if within [relcloseR*100]% change in radius

		iiDo: DO i=1,sBlob
			IF (i<sBlob) THEN
				iPlus=i+1
			ELSE
				iPlus=1
			ENDIF
			IF ( (OldSort(i,2) > 1.172) .AND. (OldSort(i,2) < 1.173) ) THEN
				PRINT *, OldSort(i,1), OldSort(i,2)
			ENDIF

			delAngle=ABS( ASortData(i,2)-ASortData(iPlus,2) )

			CheckClose: IF (delAngle < closeAngle) THEN		! need check around pt, up
				! go forward 3x and check all r values to find closest 
				kup=i; kdown=i
				DO WHILE ( delAngle<loopclose .AND. ASortData(i,3)==ASortData(kup,3) )
					kup=kup+1
					delAngle=ABS( ASortData(i,2)-ASortData(kup,2) )
				ENDDO
				kup=kup-1

				delAngle=ABS( ASortData(i,2)-ASortData(kdown,2) )
				DO WHILE ( delAngle<loopclose .AND. ASortData(i,3)==ASortData(kdown,3) )
					kdown=kdown-1
					delAngle=ABS( ASortData(i,2)-ASortData(kdown,2) ) 
				ENDDO
				kdown=kdown+1

				! now that you have the bounds set via angle, check distance for nearby
				delR= ABS(( ASortData(i,1)-ASortData(iPlus,1) ) / ASortData(i,1))
				tempUp=iPlus; flagNoSwitch=1
				DO WHILE ( tempUp<=kup .AND. 1==flagNoSwitch)
					IF ( delR>relcloseR ) THEN	! points are far apart - scan for nearer
						delR= ABS(( ASortData(i,1)-ASortData(tempUp,1) ) / ASortData(i,1))
						tempUp=tempUp+1
					ELSE
						flagNoSwitch=0
						NewSort(i,:)=OldSort(i,:)
						NewSort(iPlus,:)=OldSort(tempUp,:)
						EXIT ! fine, points are close together both in theta and r
					ENDIF
				ENDDO

				CheckBegin: IF (i>1) THEN
					delAngle=ABS( ASortData(i,2)-ASortData(kdown,2) )
					delR= ABS(( ASortData(i,1)-ASortData(i-1,1) ) / ASortData(i,1))
					tempDown=i-1; flagNoSwitch=1
					DO WHILE ( tempDown>=kdown .AND. 1==flagNoSwitch)
						IF ( delR>relcloseR ) THEN	! points are far apart - scan for nearer
							delR= ABS(( ASortData(i,1)-ASortData(tempDown,1) ) / ASortData(i,1))
							tempDown=tempDown-1
						ELSE
							flagNoSwitch=0
							NewSort(i,:)=OldSort(i,:)
							NewSort(i-1,:)=OldSort(tempDown,:)
							EXIT ! fine, points are close together both in theta and r
						ENDIF
					ENDDO			
				ENDIF CheckBegin
			ENDIF CheckClose
		ENDDO iiDo
	END SUBROUTINE CheckCloseAngle

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

	SUBROUTINE AddBCandAround( Coords, sizeLevel1, NumBlobCoords, blobPlus, blobNum, elecHeight, &
			elecWidth, elecLedge, elecRedge, maxX, minX, maxY, minY)

		! incoming variables
		INTEGER, INTENT(in)			:: NumBlobCoords, blobNum
		REAL*8, INTENT(in)			:: elecHeight, elecWidth

		! outgoing variables
		INTEGER, INTENT(out)		:: blobPlus
		INTEGER, DIMENSION(SmallInitAry), INTENT(inout) 	:: sizeLevel1
		REAL*8, INTENT(out)			:: elecLedge, elecRedge, maxX, minX, maxY, minY 
		REAL*8, DIMENSION(NumBlobCoords+20,8), INTENT(inout) :: Coords
		  ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatpan]	

		! subroutine entirely variables
		INTEGER						:: i, CaseShape
		REAL*8						:: elecThick, JustAbove, templow, temphigh, templowY

		CONTINUE
		!*************
		! main program commands begin
		templow=1e10; templowY=1e10; temphigh=0
		DO i=1,NumBlobCoords
			templow  = MIN(templow, Coords(i,1)) 
			temphigh = MAX(temphigh, Coords(i,1))
			templowY = MIN(templowY, Coords(i,2))
		ENDDO
		minX = templow - 0.4; maxX = temphigh+0.4
		minY = templowY-0.1; maxY = elecHeight+0.8; 
		i=NumBlobCoords

		elecThick = 0.1; JustAbove=0.01;
		elecLedge=0.1; elecRedge=-0.1 
		! 'elecpotential' is a global variable

		CaseShape=4			!=1, single box, simple pot; =2 ripple box, real pot; 
							!=3, full 'T', real pot; =4 is a short 'T' top shape
		SELECT CASE(CaseShape)
		CASE(1)				! test case around single box
			Coords(i+1,1)= maxX;					Coords(i+1,2)= minY
			Coords(i+2,1)= maxX;					Coords(i+2,2)= 0.5*elecHeight
			Coords(i+3,1)= maxX;					Coords(i+3,2)= 0.6*elecHeight
			Coords(i+4,1)= maxX;					Coords(i+4,2)= 0.7*elecHeight
			Coords(i+5,1)= maxX;					Coords(i+5,2)= 0.8*elecHeight
			Coords(i+6,1)= maxX;					Coords(i+6,2)= 0.9*elecHeight
			Coords(i+7,1)= maxX;					Coords(i+7,2)= elecHeight
			Coords(i+8,1)= elecLedge+elecWidth;		Coords(i+8,2)= elecHeight
			Coords(i+9,1)= elecLedge;				Coords(i+9,2)= elecHeight
			Coords(i+10,1)= elecRedge;				Coords(i+10,2)= elecHeight
			Coords(i+11,1)= elecRedge-elecWidth;	Coords(i+11,2)= elecHeight
			Coords(i+12,1)= minX;					Coords(i+12,2)= elecHeight
			Coords(i+13,1)= minX;					Coords(i+13,2)= 0.9*elecHeight
			Coords(i+14,1)= minX;					Coords(i+14,2)= 0.8*elecHeight
			Coords(i+15,1)= minX;					Coords(i+15,2)= 0.7*elecHeight
			Coords(i+16,1)= minX;					Coords(i+16,2)= 0.6*elecHeight
			Coords(i+17,1)= minX;					Coords(i+17,2)= 0.5*elecHeight
			Coords(i+18,1)= minX;					Coords(i+18,2)= 0.4*elecHeight
			Coords(i+19,1)= minX;					Coords(i+19,2)= 0.3*elecHeight
			Coords(i+20,1)= minX;					Coords(i+20,2)= minY

			! debugging, test shape, force to all Dirichlet to check code bit by bit
			Coords(i+1,3)= elecPotential;			Coords(i+1,4)= 1
			Coords(i+2,3)= elecPotential;			Coords(i+2,4)= 1
			Coords(i+3,3)= elecPotential;			Coords(i+3,4)= 1
			Coords(i+4,3)= elecPotential;			Coords(i+4,4)= 1
			Coords(i+5,3)= elecPotential;			Coords(i+5,4)= 1
			Coords(i+6,3)= elecPotential;			Coords(i+6,4)= 1
			Coords(i+7,3)= elecPotential;			Coords(i+7,4)= 1
			Coords(i+8,3)= elecPotential;			Coords(i+8,4)= 1
			Coords(i+9,3)= elecPotential;			Coords(i+9,4)= 1
			Coords(i+10,3)= elecPotential;			Coords(i+10,4)= 1
			Coords(i+11,3)= elecPotential;			Coords(i+11,4)= 1
			Coords(i+12,3)= elecPotential;			Coords(i+12,4)= 1
			Coords(i+13,3)= elecPotential;			Coords(i+13,4)= 1
			Coords(i+14,3)= elecPotential;			Coords(i+14,4)= 1
			Coords(i+15,3)= elecPotential;			Coords(i+15,4)= 1
			Coords(i+16,3)= elecPotential;			Coords(i+16,4)= 1
			Coords(i+17,3)= elecPotential;			Coords(i+17,4)= 1
			Coords(i+18,3)= elecPotential;			Coords(i+18,4)= 1
			Coords(i+19,3)= elecPotential;			Coords(i+19,4)= 1
			Coords(i+20,3)= elecPotential;			Coords(i+20,4)= 1

			Coords(i+1,5)= blobNum+1
			Coords(i+2,5)= blobNum+1
			Coords(i+3,5)= blobNum+1
			Coords(i+4,5)= blobNum+1
			Coords(i+5,5)= blobNum+1
			Coords(i+6,5)= blobNum+1
			Coords(i+7,5)= blobNum+1
			Coords(i+8,5)= blobNum+1
			Coords(i+9,5)= blobNum+1
			Coords(i+10,5)=blobNum+1
			Coords(i+11,5)=blobNum+1
			Coords(i+12,5)=blobNum+1
			Coords(i+13,5)=blobNum+1
			Coords(i+14,5)=blobNum+1
			Coords(i+15,5)=blobNum+1
			Coords(i+16,5)=blobNum+1
			Coords(i+17,5)=blobNum+1
			Coords(i+18,5)=blobNum+1
			Coords(i+19,5)=blobNum+1
			Coords(i+20,5)=blobNum+1; blobPlus=blobNum+1; sizeLevel1(blobNum+1)=20

			Coords(i+1:i+20,6)=0.5*(maxX+minX); Coords(i+1:i+20,7)=0.5*(maxY+minY);
			Coords(i+1:i+20,8)=1

		CASE(2) 			!=2 single, real pot. Duplicating points with 2 types of BC
			Coords(i+1,1)= temphigh;				Coords(i+1,2)= minY			
			Coords(i+2,1)= maxX;					Coords(i+2,2)= minY			!*
			Coords(i+3,1)= maxX;					Coords(i+3,2)= minY
			Coords(i+4,1)= maxX;					Coords(i+4,2)= elecHeight
			Coords(i+5,1)= elecLedge+elecWidth;		Coords(i+5,2)= elecHeight 	!*
			Coords(i+6,1)= elecLedge+elecWidth;		Coords(i+6,2)= elecHeight
			Coords(i+7,1)= elecLedge+elecWidth;		Coords(i+7,2)= elecHeight-elecThick
			Coords(i+8,1)= elecLedge;				Coords(i+8,2)= elecHeight-elecThick
			Coords(i+9,1)= elecLedge;				Coords(i+9,2)= elecHeight 	!* 
			Coords(i+10,1)= elecLedge;				Coords(i+10,2)= elecHeight
			Coords(i+11,1)= elecRedge;				Coords(i+11,2)= elecHeight 	!*
			Coords(i+12,1)= elecRedge;				Coords(i+12,2)= elecHeight
			Coords(i+13,1)= elecRedge;				Coords(i+13,2)= elecHeight-elecThick 
			Coords(i+14,1)= elecRedge-elecWidth;	Coords(i+14,2)= elecHeight-elecThick
			Coords(i+15,1)= elecRedge-elecWidth;	Coords(i+15,2)= elecHeight	!*
			Coords(i+16,1)= elecRedge-elecWidth;	Coords(i+16,2)= elecHeight
			Coords(i+17,1)= minX;					Coords(i+17,2)= elecHeight
			Coords(i+18,1)= minX;					Coords(i+18,2)= minY		!*
			Coords(i+19,1)= minX;					Coords(i+19,2)= minY
			Coords(i+20,1)= templow;				Coords(i+20,2)= minY

			Coords(i+1,3)= fluidPot;				Coords(i+1,4)= 1
			Coords(i+2,3)= fluidPot;				Coords(i+2,4)= 1
			Coords(i+3,3)= 0;						Coords(i+3,4)= 0			
			Coords(i+4,3)= 0;						Coords(i+4,4)= 0
			Coords(i+5,3)= 0;						Coords(i+5,4)= 0
			Coords(i+6,3)= elecPotential;			Coords(i+6,4)= 1
			Coords(i+7,3)= elecPotential;			Coords(i+7,4)= 1
			Coords(i+8,3)= elecPotential;			Coords(i+8,4)= 1
			Coords(i+9,3)= elecPotential;			Coords(i+9,4)= 1
			Coords(i+10,3)= 0;						Coords(i+10,4)= 0
			Coords(i+11,3)= 0;						Coords(i+11,4)= 0
			Coords(i+12,3)= elecPotential;			Coords(i+12,4)= 1
			Coords(i+13,3)= elecPotential;			Coords(i+13,4)= 1
			Coords(i+14,3)= elecPotential;			Coords(i+14,4)= 1
			Coords(i+15,3)= elecPotential;			Coords(i+15,4)= 1
			Coords(i+16,3)= 0;						Coords(i+16,4)= 0
			Coords(i+17,3)= 0;						Coords(i+17,4)= 0
			Coords(i+18,3)= 0;						Coords(i+18,4)= 0
			Coords(i+19,3)= fluidPot;				Coords(i+19,4)= 1
			Coords(i+20,3)= fluidPot;				Coords(i+20,4)= 1

			Coords(i+1,5)= blobNum+1
			Coords(i+2,5)= blobNum+1; sizeLevel1(blobNum+1)=2
			Coords(i+3,5)= blobNum+2
			Coords(i+4,5)= blobNum+2
			Coords(i+5,5)= blobNum+2; sizeLevel1(blobNum+2)=3 
			Coords(i+6,5)= blobNum+3
			Coords(i+7,5)= blobNum+3 
			Coords(i+8,5)= blobNum+3
			Coords(i+9,5)= blobNum+3; sizeLevel1(blobNum+3)=4
			Coords(i+10,5)=blobNum+4
			Coords(i+11,5)=blobNum+4; sizeLevel1(blobNum+4)=2 
			Coords(i+12,5)=blobNum+5
			Coords(i+13,5)=blobNum+5
			Coords(i+14,5)=blobNum+5
			Coords(i+15,5)=blobNum+5; sizeLevel1(blobNum+5)=4
			Coords(i+16,5)=blobNum+6
			Coords(i+17,5)=blobNum+6
			Coords(i+18,5)=blobNum+6; sizeLevel1(blobNum+6)=3
			Coords(i+19,5)=blobNum+7
			Coords(i+20,5)=blobNum+7; blobPlus=blobNum+7; sizeLevel1(blobNum+7)=2

			Coords(i+1:i+2,6)=0.5*(temphigh+maxX); Coords(i+1:i+2,7)=minY
			Coords(i+3:i+5,6)=0.5*(maxX+elecLedge+elecWidth); Coords(i+3:i+5,7)=0.5*(elecHeight+minY)
			Coords(i+6:i+9,6)=0.5*(2.*elecLedge+elecWidth); Coords(i+6:i+9,7)=0.5*(2.*elecHeight-elecThick)
			Coords(i+10:i+11,6)=0.25*(3.*elecLedge+elecRedge); Coords(i+10:i+11,7)=0.5*(2.*elecHeight-elecThick)
			Coords(i+12:i+15,6)=0.5*(2.*elecRedge-elecWidth); Coords(i+12:i+15,7)=0.5*(2.*elecHeight-elecThick)
			Coords(i+16:i+18,6)=minX+0.1*(maxX-minX); Coords(i+16:i+18,7)=0.5*(minY+maxY)
			Coords(i+19:i+20,6)=0.5*(minX+templow); Coords(i+19:i+20,7)=minY

			Coords(i+1:i+20,8)=1
		CASE(3)				! real case with curved electrodes and higher full "T"
			Coords(i+1,1)= maxX;					Coords(i+1,2)= minY
			Coords(i+2,1)= maxX;					Coords(i+2,2)= elecHeight
			Coords(i+3,1)= elecLedge+elecWidth;		Coords(i+3,2)= elecHeight
			Coords(i+4,1)= elecLedge+elecWidth;		Coords(i+4,2)= elecHeight-elecThick
			Coords(i+5,1)= elecLedge;				Coords(i+5,2)= elecHeight-elecThick
			Coords(i+6,1)= elecLedge;				Coords(i+6,2)= elecHeight+JustAbove
			Coords(i+7,1)= elecLedge+elecWidth;		Coords(i+7,2)= elecHeight+JustAbove
			Coords(i+8,1)= maxX;					Coords(i+8,2)= elecHeight+JustAbove
			Coords(i+9,1)= maxX;					Coords(i+9,2)= elecHeight+0.2
			Coords(i+10,1)= minX;					Coords(i+10,2)= elecHeight+0.2
			Coords(i+11,1)= minX;					Coords(i+11,2)= elecHeight+JustAbove
			Coords(i+12,1)= elecRedge-elecWidth;	Coords(i+12,2)= elecHeight+JustAbove
			Coords(i+13,1)= elecRedge;				Coords(i+13,2)= elecHeight+JustAbove
			Coords(i+14,1)= elecRedge;				Coords(i+14,2)= elecHeight-elecThick
			Coords(i+15,1)= elecRedge-elecWidth;	Coords(i+15,2)= elecHeight-elecThick
			Coords(i+16,1)= elecRedge-elecWidth;	Coords(i+16,2)= elecHeight
			Coords(i+17,1)= minX;					Coords(i+17,2)= elecHeight
			Coords(i+18,1)= minX;					Coords(i+18,2)= 0.9*elecHeight
			Coords(i+19,1)= minX;					Coords(i+19,2)= 0.8*elecHeight
			Coords(i+20,1)= minX;					Coords(i+20,2)= minY

			Coords(i+1,3)= 0;						Coords(i+1,4)= 0
			Coords(i+2,3)= 0;						Coords(i+2,4)= 0
			Coords(i+3,3)= elecPotential;			Coords(i+3,4)= 1
			Coords(i+4,3)= elecPotential;			Coords(i+4,4)= 1
			Coords(i+5,3)= elecPotential;			Coords(i+5,4)= 1
			Coords(i+6,3)= elecPotential;			Coords(i+6,4)= 1
			Coords(i+7,3)= 0;						Coords(i+7,4)= 0
			Coords(i+8,3)= 0;						Coords(i+8,4)= 0
			Coords(i+9,3)= 0;						Coords(i+9,4)= 0
			Coords(i+10,3)= 0;						Coords(i+10,4)= 0
			Coords(i+11,3)= 0;						Coords(i+11,4)= 0
			Coords(i+12,3)= elecPotential;			Coords(i+12,4)= 1
			Coords(i+13,3)= elecPotential;			Coords(i+13,4)= 1
			Coords(i+14,3)= elecPotential;			Coords(i+14,4)= 1
			Coords(i+15,3)= elecPotential;			Coords(i+15,4)= 1
			Coords(i+16,3)= 0;						Coords(i+16,4)= 0
			Coords(i+17,3)= 0;						Coords(i+17,4)= 0
			Coords(i+18,3)= 0;						Coords(i+18,4)= 0
			Coords(i+19,3)= 0;						Coords(i+19,4)= 0
			Coords(i+20,3)= 0;						Coords(i+20,4)= 0

		CASE(4)	! real case with curved electrodes and higher partial "T"
			Coords(i+1,1)= maxX;					Coords(i+1,2)= minY
			Coords(i+2,1)= maxX;					Coords(i+2,2)= elecHeight
			Coords(i+3,1)= elecLedge+elecWidth;		Coords(i+3,2)= elecHeight	!*
			Coords(i+4,1)= elecLedge+elecWidth;		Coords(i+4,2)= elecHeight
			Coords(i+5,1)= elecLedge+elecWidth;		Coords(i+5,2)= elecHeight-elecThick
			Coords(i+6,1)= elecLedge;				Coords(i+6,2)= elecHeight-elecThick
			Coords(i+7,1)= elecLedge;				Coords(i+7,2)= elecHeight+elecThick
			Coords(i+8,1)= elecLedge+elecWidth;		Coords(i+8,2)= elecHeight+elecThick !*	
			Coords(i+9,1)= elecLedge+elecWidth;		Coords(i+9,2)= elecHeight+elecThick
			Coords(i+10,1)= elecLedge+elecWidth;	Coords(i+10,2)= elecHeight+0.2
			Coords(i+11,1)= elecRedge-elecWidth;	Coords(i+11,2)= elecHeight+0.2
			Coords(i+12,1)= elecRedge-elecWidth;	Coords(i+12,2)= elecHeight+elecThick !*
			Coords(i+13,1)= elecRedge-elecWidth;	Coords(i+13,2)= elecHeight+elecThick 
			Coords(i+14,1)= elecRedge;				Coords(i+14,2)= elecHeight+elecThick
			Coords(i+15,1)= elecRedge;				Coords(i+15,2)= elecHeight-elecThick
			Coords(i+16,1)= elecRedge-elecWidth;	Coords(i+16,2)= elecHeight-elecThick
			Coords(i+17,1)= elecRedge-elecWidth;	Coords(i+17,2)= elecHeight !*
			Coords(i+18,1)= elecRedge-elecWidth;	Coords(i+18,2)= elecHeight 
			Coords(i+19,1)= minX;					Coords(i+19,2)= elecHeight
			Coords(i+20,1)= minX;					Coords(i+20,2)= minY

			Coords(i+1,3)= 0;						Coords(i+1,4)= 0
			Coords(i+2,3)= 0;						Coords(i+2,4)= 0
			Coords(i+3,3)= 0;						Coords(i+3,4)= 0
			Coords(i+4,3)= elecPotential;			Coords(i+4,4)= 1
			Coords(i+5,3)= elecPotential;			Coords(i+5,4)= 1
			Coords(i+6,3)= elecPotential;			Coords(i+6,4)= 1
			Coords(i+7,3)= elecPotential;			Coords(i+7,4)= 1
			Coords(i+8,3)= elecPotential;			Coords(i+8,4)= 1
			Coords(i+9,3)= 0;						Coords(i+9,4)= 0
			Coords(i+10,3)= 0;						Coords(i+10,4)= 0
			Coords(i+11,3)= 0;						Coords(i+11,4)= 0
			Coords(i+12,3)= 0;						Coords(i+12,4)= 0
			Coords(i+13,3)= elecPotential;			Coords(i+13,4)= 1
			Coords(i+14,3)= elecPotential;			Coords(i+14,4)= 1
			Coords(i+15,3)= elecPotential;			Coords(i+15,4)= 1
			Coords(i+16,3)= elecPotential;			Coords(i+16,4)= 1
			Coords(i+17,3)= elecPotential;			Coords(i+17,4)= 1
			Coords(i+18,3)= 0;						Coords(i+18,4)= 0
			Coords(i+19,3)= 0;						Coords(i+19,4)= 0
			Coords(i+20,3)= 0;						Coords(i+20,4)= 0

			Coords(i+1,5)= blobNum+1
			Coords(i+2,5)= blobNum+1
			Coords(i+3,5)= blobNum+1; sizeLevel1(blobNum+1)=3
			Coords(i+4,5)= blobNum+2 
			Coords(i+5,5)= blobNum+2 
			Coords(i+6,5)= blobNum+2 
			Coords(i+7,5)= blobNum+2
			Coords(i+8,5)= blobNum+2; sizeLevel1(blobNum+2)=5
			Coords(i+9,5)= blobNum+3
			Coords(i+10,5)=blobNum+3
			Coords(i+11,5)=blobNum+3
			Coords(i+12,5)=blobNum+3; sizeLevel1(blobNum+3)=4
			Coords(i+13,5)=blobNum+4
			Coords(i+14,5)=blobNum+4
			Coords(i+15,5)=blobNum+4
			Coords(i+16,5)=blobNum+4
			Coords(i+17,5)=blobNum+4; sizeLevel1(blobNum+4)=5
			Coords(i+18,5)=blobNum+5
			Coords(i+19,5)=blobNum+5
			Coords(i+20,5)=blobNum+5; blobPlus=blobNum+5; sizeLevel1(blobNum+5)=3

			Coords(i+1:i+3,6)=0.5*(maxX+elecLedge+elecWidth); Coords(i+1:i+3,7)=0.5*(elecHeight+minY)
			Coords(i+4:i+8,6)=0.5*(2.*elecLedge+elecWidth); Coords(i+4:i+8,7)=0.5*(2.*elecHeight-elecThick)
			Coords(i+9:i+12,6)=0.5*(2.*elecLedge); Coords(i+9:i+12,7)=0.5*(2.*elecHeight-elecThick)
			Coords(i+13:i+17,6)=0.5*(2.*elecRedge-elecWidth); Coords(i+13:i+17,7)=0.5*(2.*elecHeight-elecThick)
			Coords(i+18:i+20,6)=minX; Coords(i+18:i+20,7)=0.5*(elecHeight+minY)

			Coords(i+1:i+20,8)=1
		CASE DEFAULT
			STOP 'Err CaseShape'
		END SELECT

		PRINT *, '  Done with subroutine <AddBCandAround>'
	END SUBROUTINE AddBCandAround

	!---------------------------
	! calculates whether something is inside or outside a shape
	SUBROUTINE CalcPot(InsideYes, meshi, meshj, NumGrid, NumPan, PanX, PanY)
		! incoming
		INTEGER, INTENT(in)		:: NumGrid, NumPan
		REAL*8, DIMENSION(NumPan),INTENT(in) 			:: PanX, PanY

		! outgoing
		INTEGER, DIMENSION(NumGrid,NumGrid),INTENT(out) :: InsideYes 
		REAL*8, DIMENSION(NumGrid), INTENT(out) 		:: meshi, meshj

		! internal
		INTEGER					:: z, CrossCount, CL, CrossBeforePx, CrossAfterPx, CrossBeforePy, &
			CrossAfterPy, countInt, countIntersecOut, sizeCross, NoDouble, xx, yy, zz, CrossYes, &
			countIntersecIn, countDeadXY
		REAL*8					:: Ax, Bx, Ay, By, Px, Py, TriangleArea, IntersecX, IntersecY
		REAL*8, DIMENSION(4) 	:: XYPos ! xmin, xmax, ymin, ymax
		REAL*8, DIMENSION(NumPan,3) 	:: CrossLocation, NewCrossLocation


		CONTINUE
		!XYPos(1)=minX; XYPos(2)=maxX; XYPos(3)=minY-0.1; XYPos(4)=maxY-0.5
		XYPos(1)=-0.001; XYPos(2)=0.001; XYPos(3)=0.77; XYPos(4)=0.81

		XXLoop: do xx = 1,NumGrid
			meshi(xx) = XYPos(1)+(xx-1)*(XYPos(2)-XYPos(1))/(1.*NumGrid-1.)
			YYLoop: do yy = 1,NumGrid
				IF (1==xx) THEN
					meshj(yy) = XYPos(3)+(yy-1)*(XYPos(4)-XYPos(3))/(1.*NumGrid-1.)
				ENDIF
				CrossCount=0; CL=1; CrossLocation=0
				InsideYes(xx,yy)=1
				Px = meshi(xx); Py = meshj(yy)
				Zdo: DO z = 1,NumPan 
					if ( NumPan==z ) then 
						Ax=PanX(z); Ay=PanY(z)
						Bx=PanX(1); By=PanY(1)            
					else 
						Ax=PanX(z); Ay=PanY(z)   
						Bx=PanX(z+1); By=PanY(z+1)   
					endif

					!TriangeAreaDet(1,:)=[1. 1. 1.];  
					!TriangeAreaDet(2,:)=[Ax Bx Px];  
					!TriangeAreaDet(3,:)=[Ay By Py];
					TriangleArea = 0.5*ABS(Py*Bx - Px*By - Py*Ax + Ax*By + Ay*Px - Ay*Bx)
					if (ABS(TriangleArea)-2.*epsil<0) then	! if yes, pt is within segment or parallel to it
						call InSideYesNo( CrossCount, CrossYes, IntersecX, IntersecY, 'Y', Ax, Ay, &
							Bx, By, Px, Py )
						Cross1: if ( 1==CrossYes ) then			! within segment 
							do zz = 1,NumPan-1
								if ( ((PanX(zz)-IntersecX <epsil .AND. PanX(zz+1)-IntersecX >-epsil) .OR. &
									(PanX(zz)-IntersecX >-epsil .AND. PanX(zz+1)-IntersecX <epsil)) .AND. &
									((PanY(zz)-IntersecY <epsil .AND. PanY(zz+1)-IntersecY >-epsil) .OR. &
									(PanY(zz)-IntersecY >-epsil .AND. PanY(zz+1)-IntersecY <epsil)) ) then
									InsideYes(xx,yy)=0
									EXIT
								endif
							end do
						endif Cross1
					else						! general case   
						call InSideYesNo( CrossCount, CrossYes, IntersecX, IntersecY, 'Y', Ax, Ay, &
							Bx, By, Px, Py ) 
						Cross2: if ( 1==CrossYes ) then			! inside shape 
							CrossLocation(CL, 1)=IntersecX
							CrossLocation(CL, 2)=IntersecY
							CrossLocation(CL, 3)=z
							CL=CL+1
						endif Cross2
						call InSideYesNo( CrossCount, CrossYes, IntersecX, IntersecY, 'X', Ax, Ay, &
							Bx, By, Px, Py )
						Cross3: if ( 1==CrossYes ) then	        ! inside shape 
							CrossLocation(CL, 1)=IntersecX
							CrossLocation(CL, 2)=IntersecY
							CrossLocation(CL, 3)=z  
							CL=CL+1
						endif Cross3
					endif			        ! if abs(TriangleArea)
				ENDDO Zdo					! for z=            
				sizeCross = CL-1

				CheckCrossing: if ( sizeCross>1 .AND. 1 == InsideYes(xx,yy) ) then
					CrossBeforePx=0; CrossAfterPx=0
					CrossBeforePy=0; CrossAfterPy=0

					countInt=0; NewCrossLocation(1,1:3)=0
					do countIntersecOut = 1,sizeCross		! remove duplicate entries
						NoDouble=1
						if (sizeCross == countIntersecOut) then
                   			! do nothing
						else
							do countIntersecIn = countIntersecOut+1,sizeCross
								if ( abs(CrossLocation(countIntersecIn,1) - CrossLocation(countIntersecOut,1)) &
									<epsil .AND. abs(CrossLocation(countIntersecIn,2) - &
									CrossLocation(countIntersecOut,2))<epsil ) then
									NoDouble=0
								endif
							end do
						endif
						if ( 1==NoDouble ) then 
							countInt=countInt+1
							NewCrossLocation(countInt,:)=CrossLocation(countIntersecOut,:)
						endif
					end do
					do countDeadXY = 1,countInt
						if ( abs(Px -NewCrossLocation(countDeadXY,1))<epsil ) then ! vertical crossing
							if ( Py < NewCrossLocation(countDeadXY,2) ) then
								CrossAfterPx=CrossAfterPx+1
							elseif ( Py > NewCrossLocation(countDeadXY,2) ) then
								CrossBeforePx=CrossBeforePx+1
							else 
								print *, 'Cross dead on X'
								exit
							endif
						else                                        ! horizontal crossing
							if ( Px < NewCrossLocation(countDeadXY,1) ) then
								CrossAfterPy=CrossAfterPy+1
							elseif ( Px > NewCrossLocation(countDeadXY,1) ) then 
								CrossBeforePy=CrossBeforePy+1
							else 
								print *, 'Cross dead on Y'
								exit
							endif
						endif
					end do

					if ( 0==CrossAfterPy .OR. 0==CrossAfterPx .OR. 0==CrossBeforePx .OR. 0==CrossBeforePy ) then               
                		! outside because even number of crossings
						InsideYes(xx,yy)=0
					endif
				endif CheckCrossing	!sizeCross
			ENDDO YYLoop		! for yy
			IF ( 0==MOD(xx,10) ) THEN
				!PRINT *, 100.*xx/NumGrid, '% done'
			ENDIF
		ENDDO XXLoop           ! for xx
	END SUBROUTINE CalcPot

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
		real*8, intent(in)		:: Ax, Ay, Bx, By, Px, Py
		character, intent(in)	:: CrossDir

 	! outgoing variables
		INTEGER, intent(out)	:: CrossYes
		integer, intent(inout)	:: CrossCount
		real*8, intent(out)		:: IntersecX, IntersecY

	! subroutine entirely variables
		real*8					:: delBAx, delBAy, t

		continue

		CrossYes=0
		delBAx = Bx-Ax; delBAy = By-Ay

		select case( CrossDir )
		case ('Y')
			if ( 0 /= delBAy ) then 
				t = (Py - Ay)/delBAy
				IntersecX = Ax + t*delBAx
				IntersecY = Ay + t*delBAy
			else 
				t = (Px - Ax)/(delBAx)                
				IntersecX = Px; IntersecY = Ay
			endif
		case ('X')
			if ( 0 /= delBAx ) then                
				t = (Px - Ax)/delBAx  
				IntersecX = Ax + t*delBAx
				IntersecY = Ay + t*delBAy
			else 
				t = (Py - Ay)/(delBAy)   
				IntersecX = Ax; IntersecY = Py
			endif
		case default
			stop 'Incorrect CrossDir'
		end select

		if ( t>=0-epsil .AND. t<=1+epsil ) then
			CrossCount = CrossCount+1; CrossYes=1
		endif
	end SUBROUTINE InSideYesNo

! ***************
! SUBROUTINE LinFrac
! ***************
!
! DESCRIPTION: Gives the linear interpolation between two points where the levelset =0
!
! INPUT: levelset value at 2 points, position of 2 points
!
! OUTPUT: new intermediate crossing location
	FUNCTION LinFrac(LS1, LS2, Pos1, Pos2)
		! incoming variables
		REAL*8, INTENT(in)	:: LS1, LS2, Pos1, Pos2

		! outgoing variable
		REAL*8 				:: LinFrac

		! program variables
		REAL*8				:: m, b
		CONTINUE
		m=(LS2-LS1)/(Pos2-Pos1)
		b=LS1-m*Pos1

		LinFrac = -b/m

		IF (abs(LinFrac)>toolarge) THEN
			PRINT *, 'Error Shading:LinFrac'
			STOP
		ENDIF
	END FUNCTION LinFrac

! ***************
! SUBROUTINE BubbleSort
! ***************
!
! DESCRIPTION: Bubble sorts an array and makes the same interchanges in
!   an auxiliary array.  The array is sorted in decreasing order.
! 	
! http://www.personal.psu.edu/faculty/j/h/jhm/f90/examples/sort/sort2.f
! INPUT: array, size of array
!
! OUTPUT: sorted array
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      IY - array to be carried with X (all swaps of X elements are
!          matched in IY .  After the sort IY(J) contains the original
!          postition of the value X(J) in the unsorted X array.
!      N - number of values in array X to be sorted
!
!   REVISION HISTORY  (YYMMDD)
!   950310  DATE WRITTEN
!   John Mahaffy
!	050725 Anton VanderWyst - Added r, CurrentPos to integrate into my code

	SUBROUTINE BubbleSort (X, IY, N, CurrentPos)
		! incoming variables
		INTEGER, INTENT(in)	:: N, CurrentPos

		! outgoing variable
		INTEGER, DIMENSION(N), INTENT(out) :: IY
		REAL*8, DIMENSION(N), INTENT(inout) :: X

		! program variables
		INTEGER				:: I, J, JMAX, ITEMP
		REAL*8				:: TEMP

		CONTINUE

		DO I=1,N
			IY(I)=CurrentPos+I-1
		ENDDO

		JMAX=N-1
		DO I=1,N-1
			TEMP=1e20
			DO J=1,JMAX
				IF(X(J) > X(J+1)) THEN
					! do nothing
				ELSE
					TEMP=X(J)
					X(J)=X(J+1)
					X(J+1)=TEMP
					ITEMP=IY(J)
					IY(J)=IY(J+1)
					IY(J+1)=ITEMP
				ENDIF
			ENDDO
			IF(1e20==TEMP) THEN
				EXIT
			ENDIF
			JMAX=JMAX-1
		ENDDO
	END SUBROUTINE BubbleSort

! ***************
! SUBROUTINE SortNN[x]
! ***************
! DESCRIPTION: takes unsorted list of points and produces a sorted list, based on 
! nearest neighbor
!
! INPUT: xarraySP, yarraySP
!
! OUTPUT: sorted xarraySP, yarraySP
!
! CALLING PROGRAM: CheckCloseAngle
!
! Version history:
!   29 Jul 2005 copied over from SortNodePts8.f90

	SUBROUTINE SortNN( ArrayNew, sizeLevel5, ArrayOld)		
		! incoming variables
		INTEGER, INTENT(in)				    :: sizeLevel5
		REAL*8, DIMENSION(sizeLevel5,8), INTENT(in) :: ArrayOld

		! outgoing variables
		REAL*8, DIMENSION(sizeLevel5,8), INTENT(out) :: ArrayNew

		! subroutine entirely variables
		INTEGER, DIMENSION(sizeLevel5)		:: indexIn
		INTEGER								:: nearPtNow, nearPtNext, i, j, tempLeft

		REAL*8							:: outrangeVar, dist, distNeighbor
		REAL*8,DIMENSION(sizeLevel5)	:: xarrayOld, yarrayOld

		CONTINUE
		!*************
		! main program commands begin
		xarrayOld=ArrayOld(:,1); yarrayOld=ArrayOld(:,2)

		! left most point gotten as furthest left x		
		outrangeVar=1e6; tempLeft=1; indexIn(:)=1
		DO i=1,sizeLevel5
			IF (xarrayOld(i)< xarrayOld(tempLeft) ) THEN
				tempLeft=i
			ENDIF
		ENDDO
		indexIn(tempLeft)=0; nearPtNow=tempLeft
		ArrayNew(1,:)=ArrayOld(tempLeft,:)

		! loops through and organizes 
		outer: DO i=1,sizeLevel5-1
			dist=outrangeVar; nearPtNext=0
			DO j=1,sizeLevel5
				IF ( 1==indexIn(j) ) then
					distNeighbor = (xarrayOld(nearPtNow)-xarrayOld(j))*(xarrayOld(nearPtNow)-xarrayOld(j)) + &
						(yarrayOld(nearPtNow)-yarrayOld(j))*(yarrayOld(nearPtNow)-yarrayOld(j))

					IF ( distNeighbor < dist ) THEN
						dist = distNeighbor
						nearPtNext=j
					ENDIF
				ENDIF
			ENDDO

			indexIn(nearPtNext)=0; nearPtNow=nearPtNext
			ArrayNew(i+1,:)=ArrayOld(nearPtNext,:)
		ENDDO outer 
	END SUBROUTINE SortNN

! ***************
! SUBROUTINE SortNN_Norms
! ***************
! DESCRIPTION: takes theta sorted list of points and produces a new sorted list, based on 
! nearest neighbor AND change in panel normals
!
! INPUT: Coords(x8)
!
! OUTPUT: sorted Coords(x8)
!
! CALLING PROGRAM: UpdateBlobList
	SUBROUTINE SortNN_Norms( ArrayNew, sizeBubble, ArrayOld)		
		! incoming variables
		INTEGER, INTENT(in)				:: sizeBubble
		REAL*8, DIMENSION(sizeBubble,8), INTENT(in) :: ArrayOld
		  ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, type panel]

		! outgoing variables
		REAL*8, DIMENSION(sizeBubble,8), INTENT(out) :: ArrayNew

		! subroutine entirely variables
		INTEGER, DIMENSION(sizeBubble)	:: indexIn
		INTEGER							:: nearPtNow, nearPtNext, i, j, tempLeft
		REAL*8							:: outrangeVar, dist, distNeighbor, &
			delX, delY, PanLen2, deltaNorm, maxNorm, minNorm, xnormNow, ynormNow, &
			xnormNext, ynormNext, xNormtemp, yNormtemp
		REAL*8, DIMENSION(sizeBubble)	:: xOld, yOld

		CONTINUE
		!*************
		! main program commands begin
		xOld=ArrayOld(:,1); yOld=ArrayOld(:,2)

		! left-most point gotten as furthest left x		
		outrangeVar=1e6; tempLeft=1; indexIn(:)=1
		DO i=1,sizeBubble
			IF (xOld(i)< xOld(tempLeft) ) THEN
				tempLeft=i
			ENDIF
		ENDDO
		indexIn(tempLeft)=0; nearPtNow=tempLeft
		ArrayNew(1,:)=ArrayOld(tempLeft,:)
		xnormNow=-1.0; ynormNow=0.0	! b/c starts directly to right, assume a vertical panel going down

		! loops through and organizes 
		outer: DO i=1,sizeBubble-1
			dist=outrangeVar; nearPtNext=0; maxNorm=0; minNorm=1000
			DO j=1,sizeBubble
				IF ( 1==indexIn(j) ) then
					distNeighbor = (xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
						(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j))

					! calculate norm to pt
					delX = xOld(j) - xOld(nearPtNow)
					delY = yOld(j) - yOld(nearPtNow)

					PanLen2  = SQRT( delX*delX + delY*delY )
					xNormtemp= -delY/PanLen2
					yNormtemp= delX/PanLen2

					deltaNorm=MAX(ABS(xNormtemp-xnormNow), ABS(yNormtemp-ynormNow))

					IF ( maxNorm<deltaNorm ) THEN
						maxNorm=deltaNorm
					ENDIF
					IF ( minNorm>deltaNorm ) THEN
						minNorm=deltaNorm
					ENDIF
					IF ( distNeighbor < dist .AND. deltaNorm<toofarNorm ) THEN
						dist = distNeighbor
						nearPtNext=j
						xNormNext=xNormtemp; yNormNext=yNormtemp
					ENDIF
				ENDIF
			ENDDO

			IF ( 0==nearPtNext ) THEN	! no connection, just do NN
				DO j=1,sizeBubble
					IF ( 1==indexIn(j) ) then
						distNeighbor = (xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
							(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j))

						IF ( distNeighbor < dist ) THEN
							dist = distNeighbor
							nearPtNext=j
							xNormNext=xNormtemp; yNormNext=yNormtemp
						ENDIF
					ENDIF
				ENDDO
			ENDIF
			indexIn(nearPtNext)=0; 
			nearPtNow=nearPtNext; xnormNow=xnormNext; ynormNow=ynormNext
			ArrayNew(i+1,:)=ArrayOld(nearPtNext,:)
		ENDDO outer 
	END SUBROUTINE SortNN_Norms

! ***************
! SUBROUTINE SortNN_NoCross
! ***************
! DESCRIPTION: takes theta sorted list of points and produces a new sorted list, based on 
! nearest neighbor, change in panel normals AND not crossing over lines
!
! INPUT: Coords(x8)
!
! OUTPUT: sorted Coords(x8)
!
! CALLING PROGRAM: UpdateBlobList

	SUBROUTINE SortNN_NoCross( ArrayNew, sizeBubble, ArrayOld )		
		! incoming variables
		INTEGER, INTENT(in)				:: sizeBubble
		REAL*8, DIMENSION(sizeBubble,8), INTENT(in) 	:: ArrayOld
		  ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, type panel]

		! outgoing variables
		REAL*8, DIMENSION(sizeBubble,8), INTENT(out) :: ArrayNew

		! subroutine entirely variables
		LOGICAL							:: Crossed_Unknown, LoopCrossed
		INTEGER, DIMENSION(sizeBubble)	:: indexIn
		INTEGER							:: nearPtNow, nearPtNext, i, j, tempLeft
		REAL*8							:: outrangeVar, dist, distNeighbor
		REAL*8,DIMENSION(sizeBubble)	:: xOld, yOld

		CONTINUE
		!*************
		! main program commands begin
		! ** not debugged!!!

		xOld=ArrayOld(:,1); yOld=ArrayOld(:,2)

		! left-most point gotten as furthest left x		
		outrangeVar=1e6; tempLeft=1; indexIn(:)=1; Crossed_Unknown=.TRUE.
		DO i=1,sizeBubble
			IF (xOld(i)< xOld(tempLeft) ) THEN
				tempLeft=i
			ENDIF
		ENDDO

		NC_do: DO WHILE (Crossed_Unknown)
			indexIn(:)=1
			indexIn(tempLeft)=0; nearPtNow=tempLeft
			ArrayNew(1,:)=ArrayOld(tempLeft,:)

			! loops through and organizes 
			outer: DO i=1,sizeBubble-1
				dist=outrangeVar; nearPtNext=0
				DO j=1,sizeBubble
					IF ( 1==indexIn(j) ) then
						distNeighbor = ( xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
							(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j) )

						IF ( distNeighbor < dist ) THEN
							dist = distNeighbor
							nearPtNext=j
						ENDIF
					ENDIF
				ENDDO

				indexIn(nearPtNext)=0; nearPtNow=nearPtNext
				ArrayNew(i+1,:)=ArrayOld(nearPtNext,:)
			ENDDO outer 

			!check to see if there is a crossing with this setup
			DO i=1,sizeBubble-1
				DO j=1,i
					LoopCrossed=.FALSE.
					CALL DoPtsCross( LoopCrossed, ArrayNew(i,1),ArrayNew(i,2),ArrayOld(nearPtNext,1), &
						ArrayOld(nearPtNext,2), ArrayNew(j,1),ArrayNew(j,2), ArrayNew(j+1,1), &
						ArrayNew(j+1,2) )

					IF (LoopCrossed) THEN		
						! crossed, loop again. take (incorrect) final pt and force it to the nearest 
						!   open neighbor 
						PRINT *, '  Crossed link - try again'

						EXIT
					ELSE
						Crossed_Unknown=.FALSE.
					ENDIF

					indexIn(nearPtNext)=0; nearPtNow=nearPtNext
					ArrayNew(i+1,:)=ArrayOld(nearPtNext,:)
				ENDDO
			ENDDO
		ENDDO NC_do
	END SUBROUTINE SortNN_NoCross

! ***************
! SUBROUTINE SortNN_NoLong
! ***************
! DESCRIPTION: takes theta sorted list of points and produces a new sorted list, based on 
! nearest neighbor, norm change minimization AND limiting the point-point distance to a max
!
! INPUT: Coords(x8)
!
! OUTPUT: sorted Coords(x8)
!
! CALLING PROGRAM: UpdateBlobList

	SUBROUTINE SortNN_NoLong( ArrayNew, sizeBubble, ArrayOld )		
		! incoming variables
		INTEGER, INTENT(in)				    :: sizeBubble
		REAL*8, DIMENSION(sizeBubble,8), INTENT(in) 	:: ArrayOld
		  ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, type panel]

		! outgoing variables
		REAL*8, DIMENSION(sizeBubble,8), INTENT(out) :: ArrayNew

		! subroutine entirely variables
		LOGICAL								:: ReRun
		INTEGER, DIMENSION(sizeBubble)		:: indexIn, forceLink
		INTEGER								:: Crossunit, nearPtNow, nearPtNext, nearPtLast, &
			i, j, k, tempLeft, countRun, countRepeat
		REAL*8, PARAMETER				:: MultDist=4.0
		REAL*8							:: outrangeVar, dist, distNeighbor, refLength, &
			delX, delY, PanLen2, deltaNorm, maxNorm, minNorm, xnormNow, ynormNow, &
			xnormNext, ynormNext, xnormTemp, ynormTemp
		REAL*8, DIMENSION(sizeBubble)	:: xOld, yOld

		CONTINUE
		!*************
		! main program commands begin
		! ***** works only most of the time *****
		xOld=ArrayOld(:,1); yOld=ArrayOld(:,2); Crossunit=19

		! left-most point gotten as furthest left x		
		outrangeVar=1e6; tempLeft=1; forceLink(:)=-1
		DO i=1,sizeBubble
			IF (xOld(i)< xOld(tempLeft) ) THEN
				tempLeft=i
			ENDIF
		ENDDO
		nearPtNow=tempLeft; nearPtLast=nearPtNow

		ReRun=.TRUE.; countRun=0; countRepeat=0
		ReRunDo: DO WHILE (ReRun)
			102 ReRun=.FALSE.; countRun=countRun+1	! hate doing this, but need to exit out 3 loops at once
			refLength=ABS(xOld(tempLeft)-xOld(tempLeft+1))

			indexIn(:)=1; indexIn(tempLeft)=0; nearPtNow=tempLeft
			ArrayNew(1,:)=ArrayOld(tempLeft,:)

			! loops through and organizes 
			outer: DO i=1,sizeBubble-1
				forcelinkIF: IF ( forceLink(nearPtNow) /= -1 ) THEN	! you are on an iterated run and need to forcelink
					! to a particular point to avoid overall error
					nearPtNext= forceLink(nearPtNow)

					! calculate norm to pt
					delX = xOld(nearPtNext) - xOld(nearPtNow)
					delY = yOld(nearPtNext) - yOld(nearPtNow)

					PanLen2  = SQRT( delX*delX + delY*delY )
					xnormNext= -delY/PanLen2
					ynormNext= delX/PanLen2
				ELSE	
					dist=outrangeVar; nearPtNext=0; maxNorm=0; minNorm=1000
					NNNormDo: DO j=1,sizeBubble
						indexIf: IF ( 1==indexIn(j) ) THEN
							distNeighbor = (xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
								(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j))

							! calculate norm to pt
							delX = xOld(j) - xOld(nearPtNow)
							delY = yOld(j) - yOld(nearPtNow)

							PanLen2  = SQRT( delX*delX + delY*delY )
							xnormTemp= -delY/PanLen2
							ynormTemp= delX/PanLen2

							deltaNorm=MAX(ABS(xnormTemp-xnormNow), ABS(ynormTemp-ynormNow))

							IF ( maxNorm<deltaNorm ) THEN
								maxNorm=deltaNorm
							ENDIF
							IF ( minNorm>deltaNorm ) THEN
								minNorm=deltaNorm
							ENDIF
							IF ( distNeighbor < dist .AND. deltaNorm<toofarNorm .AND. &
								distNeighbor<MultDist*refLength) THEN
								dist = distNeighbor
								nearPtNext=j
								xnormNext=xnormTemp; ynormNext=ynormTemp
							ENDIF
						ENDIF indexIf
					ENDDO NNNormDo

					! *********
					! Now checking nested cases to see if the above shape is 'good'
					IF (dist<outrangeVar) THEN
						refLength=((i-1)*refLength+dist)/i	! compute a running average on the reference length
					ELSE
						! do nothing - refLength is tentatively correct as-is
					ENDIF
					NN_noNormIF: IF ( 0==nearPtNext ) THEN	! if no connection, just do nearby NN
						DO j=1,sizeBubble
							IF ( 1==indexIn(j) ) then
								distNeighbor = (xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
									(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j))

								IF ( distNeighbor < dist .AND. distNeighbor<MultDist*refLength) THEN
									dist = distNeighbor
									nearPtNext=j

									! calculate norm to pt
									delX = xOld(j) - xOld(nearPtNow)
									delY = yOld(j) - yOld(nearPtNow)

									PanLen2  = SQRT( delX*delX + delY*delY )
									xnormNext= -delY/PanLen2
									ynormNext= delX/PanLen2
								ENDIF
							ENDIF
						ENDDO
					ENDIF NN_noNormIF

					NN_noDistIF: IF ( 0==nearPtNext ) THEN	
						! if still no connection, *major* miss. Next point is wrong, because too far away; 
						!   1) find it. 2) it's NN already _taken_. 3) rerun
						dist=outrangeVar; countRepeat=countRepeat+1
						DO j=1,sizeBubble	! doing (1)
							IF ( 1==indexIn(j) ) THEN
								distNeighbor = (xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
									(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j))

								IF ( distNeighbor < dist ) THEN
									dist = distNeighbor
									nearPtNext=j
								ENDIF
							ENDIF
						ENDDO

						dist=outrangeVar
						DO j=1,sizeBubble	! doing (2)
							IF ( indexIn(j)==0 .AND. j/=nearPtNext ) THEN
								distNeighbor = (xOld(nearPtNext)-xOld(j))*(xOld(nearPtNext)-xOld(j)) + &
									(yOld(nearPtNext)-yOld(j))*(yOld(nearPtNext)-yOld(j))

								IF ( distNeighbor < dist .AND. distNeighbor<MultDist*refLength) THEN
									dist = distNeighbor
									k=j
								ENDIF
							ENDIF
						ENDDO
						IF ( outrangeVar==dist ) THEN
							PRINT *, 'Shading:SortNN_NoLong (2) err. i,j=',i,j
							STOP
						ENDIF

						ReRun=.TRUE.		! doing (3)
						forceLink(k)=nearPtNext

						PRINT *, 'Running NN_noDist. xFromOld, yFromOld, xFromNew, yFromNew, xTo, yTo=', &
							xOld(nearPtNow), yOld(nearPtNow), xOld(k), yOld(k), &
							xOld(nearPtNext), yOld(nearPtNext)

						PRINT *, '  -Writing "BlobTempCross.dat" '
						OPEN(unit=Crossunit, file='BlobTempCross.dat', status='replace')
						WRITE(unit=Crossunit, fmt='(a)') 'VARIABLES="X [cm]"'
						WRITE(unit=Crossunit, fmt='(a)') '"Y [cm]"'
						DO k=1,i
							WRITE(unit=Crossunit, fmt='(2e18.6)') ArrayNew(k,1), ArrayNew(k,2)
						ENDDO
						CLOSE (unit=Crossunit, status='keep')

						IF (countRepeat>5) THEN	! you are stuck in an interative loop. Find
							! nearest non-prior point and force connect
							dist=outrangeVar
							DO j=1,sizeBubble	
								IF ( j /= nearPtLast) THEN
									distNeighbor = (xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
										(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j))

									IF ( distNeighbor < dist ) THEN
										dist = distNeighbor
										k=j
									ENDIF
								ENDIF
							ENDDO

							! force connect point to current hanging point
							forceLink(k)=nearPtNow
						ELSE
							GOTO 102
						ENDIF
					ENDIF NN_noDistIF
				ENDIF forcelinkIF	

				indexIn(nearPtNext)=0; 
				nearPtLast=nearPtNow; nearPtNow=nearPtNext; 
				xnormNow=xnormNext; ynormNow=ynormNext
				ArrayNew(i+1,:)=ArrayOld(nearPtNext,:)
			ENDDO outer
		ENDDO ReRunDo
	END SUBROUTINE SortNN_NoLong

! ***************
! SUBROUTINE SortWeighted_AngDist
! ***************
! DESCRIPTION: takes unsorted list of points and produces a sorted list, based on 
! weighting distance and angle. Looks over only roughly nearby points
!
! INPUT: Coords(x8)
!
! OUTPUT: sorted Coords(x8)
!
! CALLING PROGRAM: UpdateBlobList

	SUBROUTINE SortWeighted_AngDist( ArrayNew, sizeBubble, ArrayOld )		
		! incoming variables
		INTEGER, INTENT(in)				    :: sizeBubble
		REAL*8, DIMENSION(sizeBubble,8), INTENT(in) 	:: ArrayOld
		  ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, type panel]

		! outgoing variables
		REAL*8, DIMENSION(sizeBubble,8), INTENT(out) :: ArrayNew

		! subroutine entirely variables
		INTEGER								:: i, j, tempLeft, nearPtNow, nearPtNext
		INTEGER, DIMENSION(sizeBubble)		:: indexIn
		REAL*8							:: WeightPt, LowestWeightPt, outrangeVar, xnormNow, &
			ynormNow, distNeighbor, delX, delY, PanLen2, xnormTemp, ynormTemp, xnormNext, &
			ynormNext, angNext, angNow
		REAL*8, PARAMETER				:: MultDist=1.0
		REAL*8, DIMENSION(sizeBubble)	:: xOld, yOld
		REAL*8, DIMENSION(sizeBubble,8) :: LocalNNCand, LocalAngleCand
			! 4 NN points, with ,1 being the closest and ,4 being the furthest
			!   and 5=j of ,1 with 8=j of ,4. 'Angle(,1)'=weight dist. 'Angle(,5)'=weight angle norm

		CONTINUE
		! main program commands begin
		! ***** does NOT work yet *****
		xOld=ArrayOld(:,1); yOld=ArrayOld(:,2); indexIn(:)=1

		! left-most point gotten as furthest left x		
		outrangeVar=1e6; tempLeft=1
		DO i=1,sizeBubble
			IF (xOld(i)< xOld(tempLeft) ) THEN
				tempLeft=i
			ENDIF
		ENDDO
		indexIn(tempLeft)=0; nearPtNow=tempLeft
		ArrayNew(1,:)=ArrayOld(tempLeft,:)
		xnormNow=-1.0; ynormNow=0.0	
		  ! b/c starts directly to right, assume a vertical panel going down

		LocalNNCand(:,:)=outrangeVar
		! find 4 nearest neighbors at each point
		distangDo: DO i=1,sizeBubble-1
			DO j=1,sizeBubble
				IF ( 1==indexIn(j) ) THEN
					distNeighbor = (xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
						(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j))

					IF ( distNeighbor < LocalNNCand(nearPtNow,1) ) THEN
						! cascade old values down
						IF ( LocalNNCand(nearPtNow,1)<LocalNNCand(nearPtNow,2) ) THEN
							IF ( LocalNNCand(nearPtNow,2)<LocalNNCand(nearPtNow,3) ) THEN
								IF ( LocalNNCand(nearPtNow,3)<LocalNNCand(nearPtNow,4) ) THEN
									LocalNNCand(nearPtNow,4) = LocalNNCand(nearPtNow,3)
									LocalNNCand(nearPtNow,8) = LocalNNCand(nearPtNow,7)
								ENDIF
								LocalNNCand(nearPtNow,3) = LocalNNCand(nearPtNow,2)
								LocalNNCand(nearPtNow,7) = LocalNNCand(nearPtNow,6)
							ENDIF
							LocalNNCand(nearPtNow,2) = LocalNNCand(nearPtNow,1)
							LocalNNCand(nearPtNow,6) = LocalNNCand(nearPtNow,5)
						ENDIF

						! update with new, closer distances
						LocalNNCand(nearPtNow,1) = distNeighbor
						LocalNNCand(nearPtNow,5) = j
					ENDIF
					IF ( distNeighbor < LocalNNCand(nearPtNow,2) .AND. &
						distNeighbor > LocalNNCand(nearPtNow,1) .AND. i<sizeBubble-1) THEN
						! cascade old values down

						IF ( LocalNNCand(nearPtNow,2)<LocalNNCand(nearPtNow,3) ) THEN
							IF ( LocalNNCand(nearPtNow,3)<LocalNNCand(nearPtNow,4) ) THEN
								LocalNNCand(nearPtNow,4) = LocalNNCand(nearPtNow,3)
								LocalNNCand(nearPtNow,8) = LocalNNCand(nearPtNow,7)
							ENDIF
							LocalNNCand(nearPtNow,3) = LocalNNCand(nearPtNow,2)
							LocalNNCand(nearPtNow,7) = LocalNNCand(nearPtNow,6)
						ENDIF

						! update with new, closer distances
						LocalNNCand(nearPtNow,2) = distNeighbor
						LocalNNCand(nearPtNow,6) = j
					ELSEIF ( i>=sizeBubble-1 ) THEN
						LocalNNCand(nearPtNow,2) = outrangeVar
						LocalNNCand(nearPtNow,6) = j
					ENDIF
					IF ( distNeighbor < LocalNNCand(nearPtNow,3) .AND. &
						distNeighbor > LocalNNCand(nearPtNow,2) .AND. i<sizeBubble-2) THEN
						! cascade old values down
						IF ( LocalNNCand(nearPtNow,3)<LocalNNCand(nearPtNow,4) ) THEN
							LocalNNCand(nearPtNow,4) = LocalNNCand(nearPtNow,3)
							LocalNNCand(nearPtNow,8) = LocalNNCand(nearPtNow,7)
						ENDIF

						! update with new, closer distances						
						LocalNNCand(nearPtNow,3) = distNeighbor
						LocalNNCand(nearPtNow,7) = j
					ELSEIF ( i>=sizeBubble-2 ) THEN
						LocalNNCand(nearPtNow,3) = outrangeVar
						LocalNNCand(nearPtNow,7) = j
					ENDIF
					IF ( distNeighbor < LocalNNCand(nearPtNow,4) .AND. &
						distNeighbor > LocalNNCand(nearPtNow,3) .AND. i<sizeBubble-3) THEN
						LocalNNCand(nearPtNow,4) = distNeighbor
						LocalNNCand(nearPtNow,8) = j
					ELSEIF ( i>=sizeBubble-3 ) THEN
						LocalNNCand(nearPtNow,4) = outrangeVar
						LocalNNCand(nearPtNow,8) = j
					ENDIF
				ENDIF
			ENDDO
			! note: LocalNNCand(i,1) is the shortest distance
			IF ( outrangeVar==LocalNNCand(nearPtNow,5) .OR. outrangeVar==LocalNNCand(nearPtNow,6) &
				.OR. outrangeVar==LocalNNCand(nearPtNow,7) .OR. outrangeVar==LocalNNCand(nearPtNow,8)) THEN
				PRINT *, 'Shading:SortWeighted err. npn,5,6,7,8=',INT(LocalNNCand(nearPtNow,5)), &
					INT(LocalNNCand(nearPtNow,6)), INT(LocalNNCand(nearPtNow,7)),INT(LocalNNCand(nearPtNow,8))
				PRINT *, 'i,nearPtNow,sizeBubble=',i,nearPtNow,sizeBubble
				STOP
			ENDIF

			LowestWeightPt=outrangeVar
			angNow=ATAN2(ynormNow, xnormNow)
			IF (angNow<0) THEN
				angNow=angNow+2*pi
			ENDIF
			nearPtNext=0

			DO j=5,8
				! calculate norm to pt
				delX = xOld(INT(LocalNNCand(nearPtNow,j))) - xOld(nearPtNow)
				delY = yOld(INT(LocalNNCand(nearPtNow,j))) - yOld(nearPtNow)

				PanLen2  = SQRT( delX*delX + delY*delY )
				xnormTemp= -delY/PanLen2
				ynormTemp= delX/PanLen2

				angNext=ATAN2(ynormTemp, xnormTemp)
				IF (angNext<0) THEN
					angNext=angNext+2*pi
				ENDIF
				LocalAngleCand(nearPtNow,j-4)=1.0*LocalNNCand(nearPtNow,j-4)/ &
					LocalNNCand(nearPtNow,1)
				LocalAngleCand(nearPtNow,j)=ABS( angNow-angNext )

				WeightPt=MultDist*LocalAngleCand(nearPtNow,j)+LocalAngleCand(nearPtNow,j-4)

				! now combined weighted score using somewhat arbitrary algorithm
				IF ( WeightPt < LowestWeightPt ) THEN
					LowestWeightPt = WeightPt
					nearPtNext=INT(LocalNNCand(nearPtNow,j))
					xnormNext=xnormTemp; ynormNext=ynormTemp
				ENDIF
			ENDDO 
			IF (0==nearPtNext) THEN
				PRINT *, 'Shading:SortAng err. nearPtNow,i=',nearPtNow,i
				PRINT *, 'NN 1,2,3,4=',LocalNNCand(nearPtNow,1),LocalNNCand(nearPtNow,2),LocalNNCand(nearPtNow,3),LocalNNCand(nearPtNow,4)
				PRINT *, 'NN 5,6,7,8=',LocalNNCand(nearPtNow,5),LocalNNCand(nearPtNow,6),LocalNNCand(nearPtNow,7),LocalNNCand(nearPtNow,8)
			ENDIF

			indexIn(nearPtNext)=0; 
			nearPtNow=nearPtNext; xnormNow=xnormNext; ynormNow=ynormNext
			ArrayNew(i+1,:)=ArrayOld(nearPtNext,:)
		ENDDO distangDo
	END SUBROUTINE SortWeighted_AngDist

! ***************
! SUBROUTINE BlobSeparate
! ***************
!
! DESCRIPTION: counts number of blobs and determines which node points are in each one
!
! INPUT: ls, ls blob number
!
! OUTPUT: updated 'insidemarker'

	RECURSIVE SUBROUTINE BlobSeparate ( LSMaxX, LSMaxY, NumBlob, WhichI, WhichJ )
		! incoming variables
		INTEGER, INTENT(in)			:: LSMaxX, LSMaxY, WhichI, WhichJ, NumBlob 

		CONTINUE
		insideMarker(WhichI,WhichJ)=NumBlob

		IF (WhichJ<LSMaxY) THEN
			IF ( phi(WhichI,WhichJ+1)>=0 .AND. -1==insideMarker(WhichI,WhichJ+1) ) THEN	! up liquid+
				CALL BlobSeparate( LSMaxX, LSMaxY, NumBlob, WhichI, WhichJ+1 )
			ENDIF
		ENDIF
		IF ( WhichI<LSMaxX) THEN
			IF ( phi(WhichI+1,WhichJ)>=0 .AND. -1==insideMarker(WhichI+1,WhichJ) ) THEN 	! right liquid+
				CALL BlobSeparate( LSMaxX, LSMaxY, NumBlob, WhichI+1, WhichJ )
			ENDIF
		ENDIF
		IF (WhichJ>1) THEN
			IF ( phi(WhichI,WhichJ-1)>=0 .AND. -1==insideMarker(WhichI,WhichJ-1) ) THEN	! down liquid+
				CALL BlobSeparate( LSMaxX, LSMaxY, NumBlob, WhichI, WhichJ-1 )
			ENDIF
		ENDIF
		IF (WhichI>1) THEN
			IF (phi(WhichI-1,WhichJ)>=0 .AND. -1==insideMarker(WhichI-1,WhichJ)) THEN	! left liquid+
				CALL BlobSeparate( LSMaxX, LSMaxY, NumBlob, WhichI-1, WhichJ )
			ENDIF
		ENDIF
	END SUBROUTINE BlobSeparate

! ***************
! SUBROUTINE Int2Str[x]
! ***************
! DESCRIPTION: convert integer to character without write()
!
! INPUT: string text, integer number
!
! OUTPUT: string of text+number
!
! from http://users.erols.com/dnagle/mppro/fthreads.f90
! accessed 05 July, 2005

	SUBROUTINE Int2Str(StrText, msg, code,exten,extenLen)
		! incoming
		INTEGER, INTENT(in)				:: code, extenLen
		CHARACTER(LEN=6), INTENT(in)	:: msg
		CHARACTER(LEN=extenLen), INTENT(in)	:: exten

		! outgoing
		CHARACTER(LEN=15), INTENT(out) 	:: StrText

		! internal
		INTEGER, parameter				:: cbuff_size=5 ! size of code string
		CHARACTER( len= cbuff_size) 	:: cbuff        ! construct string code
		CHARACTER( len= 1), PARAMETER 	:: blank = ' '    
		CHARACTER( len= 1), PARAMETER 	:: minus = '-'  ! if negative code
		INTEGER :: idigit, abs_code, lenStr             ! loop thru digits

		CONTINUE

		code_0: IF ( code /= 0 ) THEN                 	! special case zero
			!  process positive numbers only
			abs_code = abs(code)                    	! convert integer to ascii

			cbuff = blank
			idigit = cbuff_size

			each_digit: DO WHILE( abs_code > 0)    		! loop thru digits
				cbuff( idigit: idigit ) = achar( mod( abs_code, 10) + iachar( '0'))
				abs_code = abs_code / 10         		! next magnitude
				idigit = idigit - 1        				! next digit
			ENDDO each_digit                   			! loop thru digits

			!  if code was negative, add minus sign
			IF( code < 0 ) THEN
				cbuff( idigit: idigit) = minus     		! negative --> minus
			ENDIF

			lenStr=len_trim(ADJUSTL(cbuff))

			SELECT CASE (lenStr)							! assemble msg and code
			CASE (1)
				StrText = msg // '0000' // trim(ADJUSTL( cbuff)) // exten 
			CASE (2)
				StrText = msg // '000' // trim(ADJUSTL( cbuff)) // exten
			CASE (3)
				StrText = msg // '00' // trim(ADJUSTL( cbuff)) // exten
			CASE (4)
				StrText = msg // '0' // trim(ADJUSTL( cbuff)) // exten
			CASE (5)
				StrText = msg // '' // trim(ADJUSTL( cbuff)) // exten
			CASE DEFAULT
				PRINT *, 'Shading:Int2Str cbuff size err'
				STOP
			END SELECT
		ELSE                             				! do zero here
			StrText =trim(msg) // '00000' // exten      ! append zero
		ENDIF code_0                               		! code 0 or not
	END SUBROUTINE Int2Str

! ***************
! SUBROUTINE AddFileCommas
! ***************
! DESCRIPTION: reads a regular data file and writes a comma-separated one for ease of formatting.
!   Also, creates a histogram of the data
!
! INPUT: data file name
!
! OUTPUT: <none> new comma-ed and histogram data files written
!
! CALLING PROGRAM: Any
	SUBROUTINE AddFileCommas (IncomingFile, FileNameLen, exten, ExtenLen )
		! incoming
		INTEGER, INTENT(in)	:: FileNameLen, ExtenLen
		CHARACTER(len=FileNameLen), INTENT(in) 	:: IncomingFile
		CHARACTER(len=ExtenLen), INTENT(in) 	:: exten

		INTEGER, PARAMETER	:: Xunit=10, Yunit=11, Histunit=12, zoomHistunit=13, NumBin=25
		REAL, PARAMETER		:: zoomPercent=10.0

		CHARACTER			:: a
		CHARACTER(len=FileNameLen+6):: NewFileName
		INTEGER				:: i, j, sizelevel, zoomLevel, i1, i2
		INTEGER,DIMENSION(NumBin+1)	:: BinCount, zoomBinCount
		REAL				:: e1, e2, e3, e4, e5, e6, e7, e8, e9, maxMTCR, minMTCR, &
			binsize, zoombinsize
		REAL, DIMENSION(NumBin+1) 	:: HistBin, zoomHistBin

		CONTINUE
		PRINT *, '    Beginning AddFileCommas'

		maxMTCR=0; minMTCR=1e20; BinCount(:)=0; zoomBinCount(:)=0; zoomLevel=0

		IF (zoomPercent>100.0) THEN
			PRINT *, '  **Err: zoomPercent too high.'
			STOP
		ENDIF

		NewFileName= IncomingFile // '_comma' // exten
		OPEN(unit=Xunit, file=IncomingFile, status='old', action='read', err=990)
		OPEN(unit=Yunit, file=NewFileName, status='replace', action='write', err=991)
		OPEN(unit=Histunit, file='BlobHistogram.dat', status='replace', action='write', err=996)
		OPEN(unit=zoomHistunit, file='BlobHistogram_zoom.dat', status='replace', action='write', err=996)

		DO i=1,12	! skip beginning character description
			READ(UNIT=Xunit, FMT=10, ERR=992) a
		ENDDO
		sizelevel=1
		DO
			!Charge [C], Area [cm^2], Drop vol [m^3], Charge/Area [q/cm^2], Charge/Vol [C/m^3]
			!  Blob Xc [m], Blob Yc [m], Chi[#], MTCR[#atm/e-], Blob[#], Timestep
			READ(unit=Xunit, fmt=12, err=994) e1, e2, e3, e4, e5, e6, e7, e8, e9, i1, i2
			WRITE(unit=Yunit, fmt=13, err=995) e1, ',', e2, ',', e3, ',', e4, ',', &
				e5, ',', e6, ',', e7, ',', e8, ',', e9, ',', i1, ',', i2

			IF ( ABS(e9)>maxMTCR ) THEN
				maxMTCR=ABS(e9)
			ENDIF
			IF ( ABS(e9)<minMTCR ) THEN
				minMTCR=ABS(e9)
			ENDIF
			IF ( ABS(e9)<10 ) THEN
				PRINT *,' **Err: too small MTCR. i=',sizelevel
				PRINT *, e1, e2, e3, e4, e5, e6, e7, e8, e9, i1, i2
				STOP
			ENDIF

			sizelevel=sizelevel+1
		ENDDO
		994 sizelevel=sizelevel-2
		PRINT *, '   done with first loop reading. sizelevel=',sizelevel
		PRINT *, '   max/min MTCR=',maxMTCR,minMTCR
		CLOSE (unit=Yunit, status='keep', err=993)

		binsize=(maxMTCR-minMTCR)/(NumBin+1)
		zoombinsize=(zoomPercent*maxMTCR/100.0-minMTCR)/(NumBin+1)
		DO i=1,NumBin+1
			HistBin(i)=minMTCR+(i-1)*binsize
			zoomHistBin(i)=minMTCR+(i-1)*zoombinsize
		ENDDO

		REWIND(unit=Xunit, err=997)
		DO i=1,12	! skip beginning character description
			READ(UNIT=Xunit, FMT=10, ERR=992) a
		ENDDO
		PRINT *, '   beginning second histogram loop'
		DO
			READ(unit=Xunit, fmt=12, err=989) e1, e2, e3, e4, e5, e6, e7, e8, e9, i1, i2

			DO j=1,NumBin+1
				IF ( e9+0.1<HistBin(j) ) THEN
					BinCount(j-1)=BinCount(j-1)+1
					EXIT
				ENDIF
			ENDDO
			DO j=1,NumBin+1
				IF ( e9+0.1<zoomHistBin(j) ) THEN
					zoomBinCount(j-1)=zoomBinCount(j-1)+1
					EXIT
				ENDIF
			ENDDO
		ENDDO
		989 zoomLevel=SUM(zoomBinCount(1:NumBin))
		PRINT *, '   done with second loop reading. zoomlevel=',zoomLevel

		WRITE(unit=Histunit, fmt='(a)') 'VARIABLES="Mass to Charge Ratio [#atom / e<sup>-</sup> ]"'
		WRITE(unit=Histunit, fmt='(a)') '"<greek>h</greek><sub>total</sub> [%]"'
		WRITE(unit=Histunit, fmt='(a)') '"Bin Pop"'
		WRITE(unit=Histunit, fmt='(a)') 'ZONE T="main3"'
		WRITE(unit=zoomHistunit, fmt='(a)') 'VARIABLES="Mass to Charge Ratio [#atom / e<sup>-</sup> ]"'
		WRITE(unit=zoomHistunit, fmt='(a)') '"<greek>h</greek><sub>total</sub> [%]"'
		WRITE(unit=zoomHistunit, fmt='(a)') '"Bin Pop"'
		WRITE(unit=zoomHistunit, fmt='(a)') '"<greek>h</greek><sub>zoom</sub> [%]"'
		WRITE(unit=zoomHistunit, fmt='(a)') 'ZONE T="main3"'
		DO i=1,NumBin
			WRITE(unit=Histunit, fmt=14, err=988) HistBin(i), ',', &
				100.0*BinCount(i)/sizelevel, ',', BinCount(i)
			WRITE(unit=zoomHistunit, fmt=15, err=988) zoomHistBin(i), ',', &
				100.0*zoomBinCount(i)/sizelevel, ',', zoomBinCount(i), ',', &
				100.0*zoomBinCount(i)/zoomLevel
		ENDDO
		CLOSE (unit=Xunit, status='keep', err=998)
		CLOSE (unit=Histunit, status='keep', err=998)
		CLOSE (unit=zoomHistunit, status='keep', err=998)

		GOTO 999

		! --- formatting and error messages --
		10 FORMAT(a)
		12 FORMAT(9e13.5, i3, i4)
		13 FORMAT(e12.4, a, e12.4, a, e12.4, a, e12.4, a, e12.4, a, e12.4, a, e13.5, &
			a, e12.4, a, e12.4, a, i4, a, i5)
		14 FORMAT(e12.4, a, f8.3, a, i7)
		15 FORMAT(e12.4, a, f8.3, a, i7, a, f8.3)

		988 PRINT *, 'Err: Writing histogram bin populations'
		990 PRINT *, 'Err: Opening orig. file'
		991 PRINT *, 'Err: Opening replacement file'
		992 PRINT *, 'Err: Reading first dozen lines'
		993 PRINT *, 'Err: Closing 1st files'
		995 PRINT *, 'Err: Writing lines. i=',sizelevel
		996 PRINT *, 'Err: Opening histogram file'
		997 PRINT *, 'Err: Rewinding file'
		998 PRINT *, 'Err: Closing 2nd files'

	999 END SUBROUTINE AddFileCommas
END MODULE Shading
