MODULE MAT_Greens
	USE global
	USE Panels, ONLY   : AddBCandAround, CurveFit, Nexti, &
		polygon_area_2d_2, gmresIter, CalcDropletPot, CalcGridPot, GetPotGrid, &
		ForceToMeters, ZeroorCopyPanel, DoPtsCross

	IMPLICIT NONE

	INTERFACE MG
		SUBROUTINE adoublesetupwrap(a, b, c, d, e, f)
			USE global, ONLY : npoints, PANEL
			INTEGER, INTENT(in) :: a, f
			REAL*8, DIMENSION(a), INTENT(in)   :: c
			REAL*8, DIMENSION(npoints), INTENT(in)  :: e
			TYPE(PANEL), DIMENSION(a), INTENT(in)  :: b

			REAL*8, DIMENSION(a*a), INTENT(inout)  :: d
		END SUBROUTINE adoublesetupwrap

		SUBROUTINE calcEwrap(a, b, c, d, e, f, g)
			REAL*8, INTENT(in)  :: b, c, d, e, f, g

			REAL*8, INTENT(out) :: a
		END SUBROUTINE calcEwrap

		SUBROUTINE calcpotwrap(a, b, c, d, e, f)
			USE global, ONLY : npoints, PANEL
			INTEGER, INTENT(in) :: d
			REAL*8, INTENT(in)  :: b, c
			REAL*8, DIMENSION(npoints), INTENT(in)  :: e
			TYPE(PANEL), DIMENSION(d), INTENT(in)  :: f

			REAL*8, INTENT(out) :: a
		END SUBROUTINE calcpotwrap

		SUBROUTINE panelconstantswrap(a, b, c, d, e)
			USE global, ONLY : npoints, PANEL
			INTEGER, INTENT(in) :: a, e
			REAL*8, DIMENSION(npoints), INTENT(in)  :: d
			TYPE(PANEL), DIMENSION(a), INTENT(in)  :: b

			REAL*8, DIMENSION(a*a), INTENT(inout)  :: c
		END SUBROUTINE panelconstantswrap

		SUBROUTINE panelsolverealwrap(a, b, c)
			USE global, ONLY : PANEL
			INTEGER, INTENT(in) :: a
			TYPE(PANEL), DIMENSION(a), INTENT(in)  :: c

			REAL*8, DIMENSION(a*a), INTENT(inout)  :: b
		END SUBROUTINE panelsolverealwrap

		SUBROUTINE setconstantswrap(a, b, c, d)
			INTEGER, INTENT(in) :: a, c 
			REAL*8, INTENT(in)  :: b, d
		END SUBROUTINE setconstantswrap
	END INTERFACE MG

 	! used to restrict variables to prevent accidental calling outside
	PRIVATE      ! everything defaults to 'private' except those expressly labeled otherwise.
	PUBLIC  :: BoundaryForceCalc

CONTAINS 
! ***************
! Subroutine BoundaryForceCalc
! ***************
! DESCRIPTION: test program that defines the boundary points for the 2d axisymmetric shape
!
! INPUT: Coords, SelectedGrid
!
! OUTPUT: sigmas, potentials at SelectedGrid
!
! Primary author: Anton VanderWyst
! University of Michigan
! 1320 Beal Ave.
! Ann Arbor, MI 48109
! antonv@umich.edu 

	SUBROUTINE BoundaryForceCalc( BlobNum, blob_rc, blob_zc, NewNumPtsPan, MxPt_blob, radBlob, &
			thetaBlob, forceBlob, tempphi, LSMaxR, LSMaxZ, is_symmetric, rlo, zlo, &
			dr, dz, elecHeight, totalTime ) ! formal full version

  		! incoming variables
		INTEGER, INTENT(in) 	:: LSMaxR, LSMaxZ
		INTEGER, INTENT(inout)  :: is_symmetric 
		REAL*8, INTENT(in)     	:: rlo, zlo, dr, dz, elecHeight, totalTime
		REAL*8, DIMENSION(LSMaxR, LSMaxZ), INTENT(in)  :: tempphi

  		! outgoing variables
		INTEGER, INTENT(out)    :: BlobNum, MxPt_blob
		INTEGER, DIMENSION(SmallInitAry), INTENT(out)  :: NewNumPtsPan
		REAL*8, DIMENSION(SmallInitAry), INTENT(out)   :: blob_rc, blob_zc
		REAL*8, DIMENSION(SmallInitAry, InitArySize), INTENT(out) :: radBlob, thetaBlob, forceBlob

  		! program variables
    	! printing variables
		INTEGER  :: NewCoordsunit, Sigmaunit, tempLSunit, Phiunit, Phititleunit, &
			NewDropunit, NewDropunit2, PDFsingleunit, &
			NumBlobunit, PDFhighunit, BlobCentersunit, Aunit, Trackunit, Potunit, &
			BlobAreaunit, BlobVolunit, Allunit

    	! regular variables
		CHARACTER(LEN=15)		:: LSstep, cSurfCharge
		CHARACTER(LEN=12) 		:: PotCase
		LOGICAL        			:: RecalcPot

		INTEGER  :: i, i2, j, k, mm, m, z1, z2, NumCoords, totPanels, tempRange, &
			BlobPlus, NeedleNum, panelCount, pan, iPlus, NumBlobCoords, &
			numBorderCoord, numGaussPan, highblobNum, rSmoothEdge, panBlobNum, &
			MatSolveMethod, newdropNum, numBlobErr, filestat, filestat2
		INTEGER, DIMENSION(SmallInitAry)  :: sizeLevel
		INTEGER, DIMENSION(:), ALLOCATABLE :: NumTotNewPanels, tempNumPan, &
			swapblobNum, NumTotPanels, CumNumGaussPan

		REAL*8  :: r, z, delR, delZ, zp, rp, elecEdge, t1,t2, t3, timeCount, &
			maxR, minR, minR2, maxZ, minZ, electrodeWidth, a2, droprad, &
			MTCR, dropVol, chargeMassArea, chi, chargeMassVol, aa, ref, c, &
			ThreePtAngle, IgnoreImpact, CenterMove, ForceZeroR, DistMult, PotErr, LastPot, &
			maxDistStep_equalLS, AMR_testEnd, AMR_testStart
		REAL*8, DIMENSION(npoints)	:: wts, gpts
		REAL*8, DIMENSION(2)    	:: PtRef, PtC, PtA
		REAL*8, DIMENSION(3)     	:: PotPt, ddist
		REAL*8, DIMENSION(:), ALLOCATABLE :: En, temp1, temp2, temp3, temp7,  &
			panelarray, delRCoord, PanCharge, b, DCenter, rloc, &
			zloc, BlobCharge, AreaBlob, meshR, meshZ
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: temp4, temp5, Coords, PtsForce, &
			Minv, BlobTimeCharge
		TYPE(PANEL), DIMENSION(:), ALLOCATABLE :: GaussPan, tempSwapPan, orderedPan

		CONTINUE

    	!**********
    	! main program commands begin
    	!**********
  		! DEBUGGING LINE -------!
		is_symmetric=0; 
        !-- elecHeight=0.8 (set "zblob" in inputs file)  MS  --!
		electrodeWidth=0.005 ! 0.x cm = 0.00x m
		IF (counter<=0) THEN
			counter=0
			!PRINT *, '** starting wall clock **'
			!CALL cpu_time(AMR_testStart)
		ENDIF
		IF (11==counter) THEN
			!PRINT *, '** MAT_Greens:stopping on step 10 **'
			!CALL cpu_time(AMR_testEnd); timeCount=AMR_testEnd-AMR_testStart
			!WRITE(*, fmt='(a,f8.4,a)') 'Time for 10 steps ', timeCount, ' seconds.'
			!STOP
		ENDIF 

    	! defining parameters in main program
		PRINT *, '********************'
		CALL cpu_time(t3)
		WRITE(*, fmt='(a,2i6,i4,e11.4)') 'LSMaxR, LSMaxZ, is_symmetric, time=',LSMaxR, LSMaxZ, is_symmetric, totalTime
		WRITE(*, fmt='(a,i8,3e11.4)') ' step, dr, dz, elecHeight =', counter, dr, dz, elecHeight

		NewCoordsunit=14; Sigmaunit=11 ! file handles
		tempLSunit=16; Phiunit=17; Phititleunit=18; NewDropunit=19; NewDropunit2=20
		NumBlobunit=21; PDFsingleunit=22; PDFhighunit=23; BlobCentersunit=24; Aunit=25
		Trackunit=26; Potunit=27; BlobAreaunit=28;BlobVolunit=29; Allunit=30

		numBorderCoord=10 ! # of additional border coordinates are tacked on
		counter=counter+1 ! makes it possible to show how many times routine called
		MatSolveMethod=MatSolveMethodGoal

		! zeroing out all arrays initially
		ALLOCATE(phi(LSMaxR, LSMaxZ))
		ALLOCATE(NewCoords(MAX(4*LSMaxZ+10,InitArySize),8))
		ALLOCATE(alreadyChecked(LSMaxR,LSMaxZ)) ! nothing is initially checked
		ALLOCATE(insideMarker(LSMaxR,LSMaxZ)) ! nothing is initially marked inside
		ALLOCATE(NumTotPanels(SmallInitAry))
		ALLOCATE(Coords(InitArySize,8)); 
			! Coords:[r, z, pt val, pt type, blob#, blob_rc, blob_zc, flatpan] 

		DO i=1,LSMaxR
			DO j=1,LSMaxZ
				phi(i,j)=-1.
				alreadyChecked(i,j)=-1
				insideMarker(i,j)=-1
			ENDDO
		ENDDO
		k=4*LSMaxZ+10
		DO i=1,k
			NewCoords(i,1)=0.; NewCoords(i,2)=0.
			NewCoords(i,3)=0.; NewCoords(i,4)=0.
			NewCoords(i,5)=0.; NewCoords(i,6)=0.
			NewCoords(i,7)=0.; NewCoords(i,8)=0.
		ENDDO
		DO i=1,InitArySize
			Coords(i,1)=-1.; Coords(i,2)=-1.
			Coords(i,3)=-1.; Coords(i,4)=-1.
			Coords(i,5)=-1.; Coords(i,6)=-1.
			Coords(i,7)=-1.; Coords(i,8)=-1.
		ENDDO
		DO i=1,SmallInitAry
			sizeLevel(i)=-1
			ChokeCoord(i,1)=0.; ChokeCoord(i,2)=0.; ChokeCoord(i,3)=0.
			NumTotPanels(i)=0; NewNumPtsPan(i)=0
		ENDDO

  		! checking to see if inconsistent global params
		DistMult=1.
		CALL CheckErrs(DistMult);

		! create spaced points along the spline. Returns totPanels, NewCoords
  		!..... 
		CALL setconstantswrap(SkipFar, IgnoreImpact, RecastDist, HardwireCharge)
		IF (1==UseGaussQuad) THEN
			CALL gaussvals(gpts,wts)  	! gives Gaussian quadrature panel weights
		ELSE
			CALL linearvals(gpts,wts)	! gives linear quadrature panel weights
			PRINT *, '  *using linear panel placement'
		ENDIF
 		!..... 

  		!-----------------------------
		IF (1==LoadPastStep) THEN
   			! debugging step, loading latest number step
   			! writes a data file with the step number and electric field
			counter=1  ! makes it possible to show how many times routine called

   			! convert data number to string file name
   			!CALL Int2Str(LSstep, 'visual00000', counter,'',0)
			LSstep='Phi_pass.dat'

			OPEN(unit=tempLSunit, file=LSstep, status='old', IOSTAT=filestat)
			panelCount=0
			IF (0==filestat) THEN ! file exists
				DO i=1,LSMaxR
					DO j=1,LSMaxZ
						panelCount=panelCount+1
						READ(UNIT=tempLSunit, FMT=12, end=232) k, k, a2, phi(i,j), a2
					ENDDO
				ENDDO 
				232 CLOSE (UNIT=tempLSunit, STATUS='keep')
				PRINT *, '      ', LSstep, ' read.'
			ELSE
				PRINT *, 'File LSstep does not exist'
				STOP
			ENDIF
		ELSE
			!   phi=tempphi  ! use unaltered original data
			DO i=1,LSMaxR
				DO j=1,LSMaxZ
					phi(i,j)=tempphi(i,j)
				enddo
			enddo
		ENDIF
		ShowPhi:IF (1==ShowPassedPhi) THEN
   			! write inputted 'phi to a file'
			PRINT *, '  -Writing "Phi_pass.dat" '
			OPEN(unit=Phiunit, file='Phi_pass.dat', status='replace', IOSTAT=filestat)
			OPEN(unit=Phititleunit, file='Phi_pass2.dat', status='replace', IOSTAT=filestat2)  

			IF (0==filestat .AND. 0==filestat2) THEN
				WRITE(unit=Phititleunit, fmt='(a)') 'VARIABLES="R [cm]"'
				WRITE(unit=Phititleunit, fmt='(a)') '"Z [cm]"'
				WRITE(unit=Phititleunit, fmt='(a)') '"phi"'
				WRITE(unit=Phititleunit, fmt='(a)') 'ZONE T="main2"'
				WRITE(unit=Phititleunit, fmt='(a,i4,a,i4,a)') 'I=', LSMaxZ, ', J=', LSMaxR, ', K=1, F=POINT'

				DO i=1,LSMaxR
					r = (i-1)*dr + rlo
					DO j=1,LSMaxZ
						z = (j-1)*dz + zlo
						WRITE(unit=Phiunit, FMT=12) i, j, r, phi(i,j), z 
						WRITE(unit=Phititleunit, fmt='(3e14.6)') r, z, phi(i,j)
					ENDDO
				ENDDO
				CLOSE (unit=Phiunit, status='keep')
				CLOSE (unit=Phititleunit, status='keep')
			ELSE
				PRINT *, 'Phi_pass.dat file error'
				STOP
			ENDIF
		ENDIF ShowPhi

  		! determine all the surface points from passed levelset. 'UpdateBlobList' subfunc. assigns
  		!   default drop values
  		!.....
		CALL ShapeGrow( Coords, BlobNum, NeedleNum, sizeLevel, NumBlobCoords, &
			LSMaxR, LSMaxZ, rlo, zlo, dr, dz, is_symmetric, minR2 )
  		!.....

		IF (1==EqualLsBemNumCoord) THEN	! force roughly same number of points in NewCoord as LS from Coords
			maxDistStep_equalLS=0.7*(LSMaxR*dr + LSMaxR*dz+StepMult*(LSMaxR*dr + LSMaxR*dz))/NumBlobCoords
			PRINT *, '    Changing maxDistStep to,from=',maxDistStep_equalLS,maxDistStep
		ELSE
			maxDistStep_equalLS=maxDistStep
		ENDIF
		DEALLOCATE(insideMarker); DEALLOCATE(alreadyChecked)

		ForceZeroR=1e6 ! >= R value at which to set needle evolution e-field force to zero
		IF (1==ForceEdgeStable) THEN ! forces right edge stable
   			!.....
			CALL ForceRange(Coords, ForceZeroR, NeedleNum, sizeLevel, BlobNum)
			PRINT *, '    ZeroR= ', ForceZeroR
  	 		!.....
		ENDIF
		NumCoords= NumBlobCoords+numBorderCoord

  		!----------------------------- 
  		! add side N. BC and electrodes to all blobs
		BlobPlus=-1; elecEdge=0.1
		ALLOCATE(temp5(NumCoords,8)); 
		DO i=1,NumCoords
			temp5(i,1)=Coords(i,1); temp5(i,2)=Coords(i,2)
			temp5(i,3)=Coords(i,3); temp5(i,4)=Coords(i,4)
			temp5(i,5)=Coords(i,5); temp5(i,6)=Coords(i,6)
			temp5(i,7)=Coords(i,7); temp5(i,8)=Coords(i,8)
		ENDDO

		maxR=LSMaxR*dr+rlo
     	!..... 
		CALL AddBCandAround ( temp5, sizeLevel, NumBlobCoords, BlobPlus, BlobNum, &
			elecHeight, electrodeWidth, elecEdge, maxR, minR, maxZ, minZ, &
			numBorderCoord, NeedleNum )
  		!.....
		DO i=1,NumCoords
			Coords(i,1)=temp5(i,1); Coords(i,2)=temp5(i,2)
			Coords(i,3)=temp5(i,3); Coords(i,4)=temp5(i,4)
			Coords(i,5)=temp5(i,5); Coords(i,6)=temp5(i,6)
			Coords(i,7)=temp5(i,7); Coords(i,8)=temp5(i,8)
		ENDDO
		DEALLOCATE(temp5)

		pan=1; NumTotPanels(1)=1
		DO i=1,BlobPlus
			blob_rc(i)=Coords(pan,6)
			blob_zc(i)=Coords(pan,7)

			IF (i==NeedleNum) THEN
				WRITE(*, fmt='(a,e9.3,a)') '   Base electric field O(', -DistMult* &
					elecPotential/( (Coords(NumBlobCoords+NumBlobCoords+numBorderCoord,2) - &
					Coords(pan,2))*1e9 ), ') [V/nm]'

    			! sets 'IgnoreImpact' value
				IgnoreImpact=2.*SQRT( (elecHeight-Coords(pan,2))*(elecHeight-Coords(pan,2)) + &
					(elecEdge+electrodeWidth)*(elecEdge+electrodeWidth) )
				IF (1==SkipFar) THEN
					WRITE(*, fmt=17) '   Ignoring points further away then ', IgnoreImpact
				ENDIF
			ENDIF
			IF (i>1) THEN
				NumTotPanels(i)=NumTotPanels(i-1)+sizeLevel(i-1)
			ENDIF
			pan=pan+sizeLevel(i)
		ENDDO

		PRINT *, 'There are', BlobNum, ' blobs and a total', BlobPlus, ' shapes.'
  		! SWITCHES blob numbers to match those in past timesteps or leaves unchanged
		ALLOCATE(DCenter(BlobNum))
		ALLOCATE(swapblobNum(BlobNum))
		ALLOCATE(BlobCharge(BlobPlus))
		ALLOCATE(NumTotNewPanels(BlobPlus))
		ALLOCATE(GaussPan(1)); CALL ZeroorCopyPanel(GaussPan(1), 1)
		ALLOCATE(CumNumGaussPan(1)); CumNumGaussPan(1)=0
		DO i=1,BlobPlus
			BlobCharge(i)=-1.; NumTotNewPanels(i)=-1

			IF (i<=BlobNum) THEN
				DCenter(i)=-1.
				swapblobNum(i)=-1
			ENDIF
		ENDDO

		numBlobErr=1; RecalcPot=.FALSE.; PotErr=2.; 
		NewNumPtsPan(1)=1; LastPot=(elecPotential-fluidPot)/100.+fluidPot
		IF (1==ConnectDrops .AND. BlobNum>1) THEN
			ALLOCATE(tempNumPan(BlobNum))
			ALLOCATE(temp2(BlobNum))
			ALLOCATE(temp3(BlobNum))
			ALLOCATE(temp5(NumCoords,8))
			ALLOCATE(AreaBlob(BlobNum))

			DO i=1,NumCoords
				temp5(i,1)=Coords(i,1); temp5(i,2)=Coords(i,2)
				temp5(i,3)=Coords(i,3); temp5(i,4)=Coords(i,4)
				temp5(i,5)=Coords(i,5); temp5(i,6)=Coords(i,6)
				temp5(i,7)=Coords(i,7); temp5(i,8)=Coords(i,8)
			ENDDO

			pan=1
			DO i=1,BlobNum
				tempNumPan(i)=sizeLevel(i)
				temp2(i)=blob_rc(i)
				temp3(i)=blob_zc(i);
				AreaBlob(i)=-1.

				k=sizeLevel(i)
				ALLOCATE(temp4(k,2)); 
				DO j=1,k
					temp4(j,1)=Coords(pan+j-1,1); temp4(j,2)=Coords(pan+j-1,2)
				ENDDO

    			!..... 
				CALL polygon_area_2d_2 ( AreaBlob(i), temp4, k )
    			!..... 
				pan=pan+k

				DEALLOCATE(temp4)

				IF (ABS(AreaBlob(i))<smallepsil) THEN
					PRINT *, '  Zero areablob, giving small area. i,k,AreaBlob(i)=',i,k,AreaBlob(i)
					AreaBlob(i)=smallepsil
				ENDIF
			ENDDO

			CenterMove=0.; newdropNum=-1
   			!.....
			CALL TrackDrop( newdropNum, swapblobNum, BlobCharge, DCenter, CenterMove, &
				blob_rc, blob_zc, AreaBlob, BlobNum, NeedleNum )
   			!.....
   			! swapping the blob numbers around if needed. NeedleNum might /=1
			pan=1
			ALLOCATE(temp7(BlobNum)); temp7(:)=AreaBlob(:)
			DO i=1,BlobNum
				j=swapblobNum(i)
				k=tempNumPan(i)
				AreaBlob(j)=temp7(i)
				sizeLevel(i)=k

				PRINT *, '    oldblob#, newblob#, sizeblob, startpan=', &
					i, j, k, pan
				DO mm=1,k
					Coords(pan+mm-1,1)=temp5(NumTotPanels(i)+mm-1,1)
					Coords(pan+mm-1,2)=temp5(NumTotPanels(i)+mm-1,2)
					Coords(pan+mm-1,3)=temp5(NumTotPanels(i)+mm-1,3)
					Coords(pan+mm-1,4)=temp5(NumTotPanels(i)+mm-1,4)
					Coords(pan+mm-1,5)=i
					Coords(pan+mm-1,6)=temp5(NumTotPanels(i)+mm-1,6)
					Coords(pan+mm-1,7)=temp5(NumTotPanels(i)+mm-1,7)
					Coords(pan+mm-1,8)=temp5(NumTotPanels(i)+mm-1,8)
				ENDDO
				pan=pan+k

				IF (swapblobNum(i)/=swapblobNum(NeedleNum)) THEN
					numBlobErr=swapblobNum(i)
				ENDIF
			ENDDO
			DEALLOCATE(tempNumPan); DEALLOCATE(temp5)
			DEALLOCATE(temp2); DEALLOCATE(temp3); DEALLOCATE(temp7)
			NeedleNum=swapblobNum(NeedleNum)
			MxPt_blob=sizeLevel(NeedleNum)

			PRINT *, '    Beginning iteration for droplet potential'
			ALLOCATE(BlobTimeCharge(BlobNum, MaxNumInpotIter))
			DO i=1,BlobNum
				DO mm=1,MaxNumInpotIter
					BlobTimeCharge(i, mm)=0.
				ENDDO
			ENDDO

   			! *** MAIN iterative routine to create NewCoords, GaussPan and find sigma values!
			DO z2=1,MaxNumInpotIter
				IF (PotErr>0.01 .OR. z2<=2) THEN
					CALL InnerPotSet
				ELSE
					PRINT *, '   Done early, after ',z2-1,' iterations '
					EXIT !exit early
				ENDIF
			ENDDO
			DEALLOCATE(BlobTimeCharge)
		ELSE  ! don't change anything
			ALLOCATE(AreaBlob(BlobNum))
			pan=1
			DO z1=1,BlobNum
				k=sizeLevel(z1)
				ALLOCATE(temp4(k,2)); 
				DO j=1,k
					temp4(j,1)=Coords(pan+j-1,1); temp4(j,2)=Coords(pan+j-1,2)
				ENDDO

    			!..... 
				CALL polygon_area_2d_2 ( AreaBlob(z1), temp4, k )
				!..... 
				DEALLOCATE(temp4)
				pan=pan+k

				swapblobNum(z1)=z1
			ENDDO
			CenterMove=-1.
  			! *** MAIN iterative routine to create NewCoords, GaussPan and find sigma values!
			CALL InnerPotSet
   			! **********
		ENDIF
		IF (ALLOCATED(NumTotPanels)) DEALLOCATE(NumTotPanels)

  		! writing the data file of how many blobs existed. Used for seeing when another 
  		!   one breaks off
		PRINT *, '  -Writing "NumBlob.dat" '
		OPEN(unit=NumBlobunit, file='NumBlob.dat', status='old', IOSTAT=filestat)
		IF (0==filestat) THEN
			READ(unit=NumBlobunit, fmt='(i4)', end=494) z1
			494 REWIND(unit=NumBlobunit)
			WRITE(unit=NumBlobunit, fmt='(i4)') BlobNum     
			CLOSE(unit=NumBlobunit, status='keep')
		ELSE
			PRINT *, 'NumBlob.dat file error'
			STOP
		ENDIF
		IF (1==ShowCoords) THEN
			CALL PrintCoords(Coords, NumCoords)
		ENDIF

		ShowPanDetails: IF (1==ShowSubPan) THEN
			!.....
			CALL PrintGaussPan(GaussPan,numGaussPan)
			!.....
		ENDIF ShowPanDetails

		IF (1==RecastDist ) THEN  ! changes matrix to m-based system instead of centimeters.
			CALL ForceToMeters(GaussPan, numGaussPan, 1);
		ENDIF
		IF (1==ConnectDrops .AND. BlobNum>1) THEN
			PRINT *, '  **Done iterating for droplet potential'
		ENDIF

		IF (1==ShowAMat) THEN  ! write out the "A" matrix
			PRINT *, '  -Writing "Amat.dat" '
			OPEN(unit=Aunit, file='Amat.dat', status='replace', IOSTAT=filestat)
			IF (0==filestat) THEN
				WRITE(unit=Aunit, fmt='(a)') 'VARIABLES="Cell From [#]"'
				WRITE(unit=Aunit, fmt='(a)') '"Cell To [#]"'
				WRITE(unit=Aunit, fmt='(a)') '"A"'
				WRITE(unit=Aunit, fmt='(a)') 'ZONE T="main"'
				WRITE(unit=Aunit, fmt='(a,i3,a,i3,a)') 'I=', numGaussPan, ', J=', numGaussPan, &
					', K=1, F=POINT'

				OPEN(unit=Sigmaunit, file='Sigma.dat', status='replace')
				DO i=1,numGaussPan
					WRITE(unit=Sigmaunit, fmt='(e15.6)') GaussPan(i)%str
					DO j=1,numGaussPan
						WRITE(unit=Aunit, fmt='(2i8, e13.4)') i, j, panelarray((i-1)*numGaussPan+j)
					ENDDO
				ENDDO
				CLOSE (unit=Sigmaunit, status='keep')
				CLOSE (unit=Aunit, status='keep')
			ELSE
				PRINT *, 'Amat.dat file error'
				STOP
			ENDIF
		ENDIF

  		! Determines electric field by taking 3 points in the normal direction, calculating the potential
  		!   and taking the gradient
		ALLOCATE(En(totPanels))
		DO i=1,totPanels
			En(i)=-1.
		ENDDO
		DO i=1,panBlobNum
			CalcEnIf: IF (0==AformDouble) THEN ! differentiate potentials to get Efield
				PotPt(1)=GaussPan(i)%midPot; ddist(1)=0; j=1

				DO j=2,3
    				! calculate the points out from the surface
					rp=GaussPan(i)%rnrm*GaussPan(i)%length*0.5*(j-1)+GaussPan(i)%midr
					zp=GaussPan(i)%znrm*GaussPan(i)%length*0.5*(j-1)+GaussPan(i)%midz
					ddist(j)=SQRT((GaussPan(i)%midr-rp)*(GaussPan(i)%midr-rp) + &
						(GaussPan(i)%midz-zp)*(GaussPan(i)%midz-zp))

					IF (HardwireChargeYes/=1) THEN
						!..... 
						CALL calcpotwrap(PotPt(j), DistMult*zp, DistMult*rp, numGaussPan, wts, GaussPan)
					ELSE ! using 'PotPt' to calculate the electric field, not the potential at this point
      					!CALL calcEwrap(PotPt(j), zp, rp, (LSMaxZ*dz+zlo)/DistMult)
					ENDIF
				ENDDO

    			! Gear's method, O(h^2) -> dF/dx = 1/2delT * (3F(n) - 4F(n-1) + F(n-2))
    			! f(r)=z=r^3+a1r^2+a2r+a3; df/dr |_r=0 = 3r^2+2*a1*r^2+a2 = a2
    			!a3 = PotPt(1)
				a2 = ( PotPt(3)-PotPt(1) - ddist(3)**3 + ddist(2)*ddist(3)*ddist(3) - &
					(ddist(3)*ddist(3)/(ddist(2)*ddist(2)))*(PotPt(2)-PotPt(1)) ) / &
					(ddist(3)-ddist(3)*ddist(3)/ddist(2))
    			!a1 =(PotPt(2)-PotPt(1))/(ddist(2)*ddist(2)) - ddist(2) - a2/ddist(2)

				En(i) = a2
			ELSE ! get Efield directly, using a different D/N formulation
				En(i) = GaussPan(i)%str
			ENDIF CalcEnIf
		ENDDO
		IF (HardwireChargeYes/=1) PRINT *, '    En(4)=', En(4)
		PRINT *, 'Done after panel gradients'

  		!-----------------------------
  		! assigns all above values to output vector 'Pts_and_force' for the non-boundary pts
  		!   [xpt, ypt, xnormal, ynormal, radius to blob center, (0-2pi radians) theta to 
  		!   blob center, blob number, blob_x, blob_y, panel force] for each panel
		IF (ALLOCATED(PtsForce)) DEALLOCATE(PtsForce)
		ALLOCATE(PtsForce(numGaussPan,2)); 
		DO i=1,numGaussPan
			PtsForce(i,1)=-1.; PtsForce(i,2)=-1.
		ENDDO

		pan=1
		iBlobDo: DO i=1,BlobNum
			i2=NewNumPtsPan(i)  

  		 	! Create temp horizontal reference line
			aa=0.1; PtC(1)=blob_rc(i); PtC(2)=blob_zc(i)
			PtRef(1)=PtC(1)+0.1; PtRef(2)=PtC(2)

			jsizeBlobDo: DO j=1,i2
				IF (pan>numGaussPan) THEN
					PRINT *, 'Shading:BFC Pts_force panel err. panel,i,j=', pan, i,j
					STOP
				ENDIF
				PtsForce(pan,1)=GaussPan(pan)%r0; PtsForce(pan,2)=GaussPan(pan)%z0

				PtA(1)=GaussPan(pan)%r0; PtA(2)=GaussPan(pan)%z0

				ref = SQRT( (PtA(1)-PtC(1))*(PtA(1)-PtC(1)) + (PtA(2)-PtC(2))*(PtA(2)-PtC(2)) )
				c = SQRT( (PtA(1)-PtRef(1))*(PtA(1)-PtRef(1)) + (PtA(2)-PtRef(2))*(PtA(2)-PtRef(2)) )

    			! 0 < acos(x) < pi. Move to correct quadrant
				IF (0==ref) THEN
					ThreePtAngle=0
				ELSE
					ThreePtAngle=( (aa*aa) + (ref*ref) - (c*c) )/(2*aa*ref)
				ENDIF
				radBlob(i,j)=ref

				Flatcheck: IF ( ABS(ThreePtAngle+1)<eps ) THEN ! too close to -1
					thetaBlob(i,j)=pi
				ELSEIF ( ABS(ThreePtAngle-1)<eps ) THEN   ! too close to 1
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

				IF (1==HardwireChargeYes) THEN
					IF (50==pan) THEN
						PRINT *, 'debugging ewrap'
					ENDIF
					CALL calcEwrap( forceBlob(i,j), En(pan), GaussPan(pan)%z0, &
						GaussPan(pan)%z1, GaussPan(pan)%r0, &
						GaussPan(pan)%r1, maxZ/distMult )

					IF (1==i .AND. 2==j) PRINT *, '    i,j,forceBlob(i,j)=', i,j,forceBlob(i,j)
				ELSE
     				! Force/m^2 = -e_0 * En^2
					forceBlob(i,j)= -epsil0*En(pan)*En(pan)
					IF (1==i .AND. 2==j) PRINT *, '    i,j,forceBlob(i,j)=', i,j,forceBlob(i,j)
				ENDIF

				pan=pan+1
			ENDDO jsizeBlobDo
		ENDDO iBlobDo

		pan=1
		IF ( 1==ConnectDrops ) THEN 
   			! writes charges and locations of surface
			DO i=1,BlobNum
				IF ( NeedleNum==INT(NewCoords(pan,5)) ) THEN
					EXIT
				ENDIF
				IF (i>1) THEN
					pan=pan+NewNumPtsPan(i-1)
				ENDIF
			ENDDO

			IF (ShowEachSurfCharge/=1) THEN
				OPEN(unit=Trackunit, file='SurfaceCharge.dat', status='replace', IOSTAT=filestat)
				PRINT *, '  -Writing "SurfaceCharge.dat"'
			ELSE
				CALL Int2Str(cSurfCharge, 'SurfChr0000', counter-1,'.dat',4)
				OPEN(unit=Trackunit, file=cSurfCharge, status='replace', IOSTAT=filestat)
				PRINT *, '  -Writing ', cSurfCharge
			ENDIF

			IF (0==filestat) THEN
				WRITE(unit=Trackunit, fmt='(a,i5,a)') 'TITLE="Past surface location and charges for ', NewNumPtsPan(NeedleNum), ' panels"'
				WRITE(unit=Trackunit, fmt='(a)') 'VARIABLES="R [cm]"'
				WRITE(unit=Trackunit, fmt='(a)') '"Z [cm]"'
				WRITE(unit=Trackunit, fmt='(a)') '"Charge [C]"'
				WRITE(unit=Trackunit, fmt='(a)') '"Electric field [V/m]"'
				WRITE(unit=Trackunit, fmt='(a)') 'ZONE T="main"'
				WRITE(unit=Trackunit, fmt='(e15.4)') totalTime
				WRITE(unit=Trackunit, fmt='(2i5)') NewNumPtsPan(NeedleNum), NeedleNum

				j=NewNumPtsPan(NeedleNum)
				DO i=pan,pan+j-1
					WRITE(unit=Trackunit, fmt='(2e15.6, 2e13.4)') DistMult*GaussPan(i)%r0, &
						DistMult*GaussPan(i)%z0, ABS(epsil0*En(i)*GaussPan(i)%length), &
						En(i)
				ENDDO
				CLOSE (unit=Trackunit, status='keep')
			ELSE
				PRINT *, 'Surf[xxx].dat file error'
				STOP
			ENDIF
		ENDIF
		pan=0
		IF (1==ShowEachGaussPan .AND. 0==MOD(counter-1, SkipSteps)) THEN
			!IF (1==ShowEachGaussPan) THEN
			CALL Int2Str(cSurfCharge, 'AllSurf0000', counter-1,'.dat',4)
			OPEN(unit=Allunit, file=cSurfCharge, status='replace', IOSTAT=filestat)
			IF (0==filestat) THEN
				PRINT *, '  -Writing ', cSurfCharge
				WRITE(unit=Allunit, fmt='(a,i5,a)') 'TITLE="Past surface location"'
				WRITE(unit=Allunit, fmt='(a)') 'VARIABLES="R [cm]"'
				WRITE(unit=Allunit, fmt='(a)') '"Z [cm]"'
				WRITE(unit=Allunit, fmt='(a)') '"Blob num"'
				WRITE(unit=Allunit, fmt='(a)') 'ZONE T="main"'
				DO i=1,BlobNum
					j=NewNumPtsPan(i)
					DO k=1,j
						pan=pan+1
						WRITE(unit=Allunit, fmt='(2e15.6,i4)') DistMult*GaussPan(pan)%r0, &
							DistMult*GaussPan(pan)%z0, i
					ENDDO
				ENDDO
				CLOSE (unit=Allunit, status='keep')
			ELSE
				PRINT *, 'AllSurf[xxx].dat file error'
				STOP
			ENDIF
		ENDIF

  		!-----------------------------
  		! Smooth the edges out for the needle shape
		PRINT *, '    Forced smoothing for ',rSmoothEdge,' panels.'
		PRINT *, '  NeedleNum,NewNumPtsPan(NeedleNum)=',NeedleNum,NewNumPtsPan(NeedleNum)
		DO j=1,rSmoothEdge-1
			forceBlob(NeedleNum,NewNumPtsPan(NeedleNum)-j+1)=0
		ENDDO

  		!-----------------------------
  		! Write number of blobs and size/charge in scenario
		PRINT *, '  -Writing "BlobCenters.dat" '
		OPEN(UNIT=BlobCentersunit, FILE='BlobCenters.dat', STATUS='replace', &
			ACTION='write', IOSTAT=filestat)
		IF ( filestat/=0 ) THEN
			PRINT *, 'BlobCenters.dat file error, filestat=',filestat
			STOP
		ENDIF

		IF ( 1==ShowBlobArea ) THEN
			OPEN(UNIT=BlobAreaunit, FILE='BlobAreas.dat', STATUS='old', &
				ACCESS='sequential', POSITION='append', ERR=123)
			GOTO 124
			123 OPEN(UNIT=BlobAreaunit, FILE='BlobAreas.dat', STATUS='new', &
				ACCESS='sequential')
			WRITE(unit=BlobAreaunit, fmt='(a,i5,a,i5,a)') 'TITLE="Shape area using ', NewNumPtsPan(NeedleNum), &
				' panels and ', LSMaxR, ' points across"'
			WRITE(unit=BlobAreaunit, fmt='(a)') 'VARIABLES="R<sub>c</sub> [cm]"'
			WRITE(unit=BlobAreaunit, fmt='(a)') '"Z<sub>c</sub> [cm]"'
			WRITE(unit=BlobAreaunit, fmt='(a)') '"Area Err [%]"'
			WRITE(unit=BlobAreaunit, fmt='(a)') '"Area [cm<sup>2</sup>]"'
			WRITE(unit=BlobAreaunit, fmt='(a)') '"BlobNum [#]"'
			WRITE(unit=BlobAreaunit, fmt='(a)') '"Timestep [#]"'
		124 ENDIF

		ALLOCATE(PanCharge(BlobNum))
		tempRange=0
		DO i=1,BlobNum
			PanCharge(i)=BlobCharge(i)

   			! drop (r,z) new position, blob num, pancharge, panarea
			WRITE(UNIT=BlobCentersunit, FMT=15) blob_rc(i), blob_zc(i), i, &
				PanCharge(i), AreaBlob(i)

			IF (1==ShowBlobArea) THEN
				CALL cpu_time(t2); timeCount=t2-t3
				WRITE(*, fmt='(a,f8.4,a)') 'Time to blobArea ', timeCount, ' seconds.'

				WRITE(UNIT=BlobAreaunit, fmt='(3e12.4, e15.7, i3, i5)') blob_rc(i), &
					blob_zc(i), AreaBlob(i), AreaBlob(i), i, counter-1
			ENDIF
		ENDDO
		CLOSE(UNIT=BlobCentersunit, STATUS='keep')
		IF (1==ShowBlobArea) THEN
			CLOSE(UNIT=BlobAreaunit, STATUS='keep')
			CLOSE(UNIT=BlobVolunit, STATUS='keep')
		ENDIF

		DEALLOCATE(NewCoords); DEALLOCATE(PtsForce)

		MultBlobNumIf: IF ( BlobNum>1 ) THEN 
   			! Create snapoff PDF
			PRINT *, 'There are snapped off blobs! Number detached=', BlobNum-1

			PRINT *, '  -Writing "PDFBlob.dat"'
			OPEN(unit=NewDropunit, file='PDFBlob.dat', status='old', access='sequential', &
				position='append', IOSTAT=filestat)
			OPEN(unit=NewDropunit2, file='PDFBlob_comma.dat', status='old', &
				access='sequential', position='append', IOSTAT=filestat2)

			IF (0==filestat .AND. 0==filestat2) THEN
				DO i=1,BlobNum
					PrintCharge: IF (i /=NeedleNum) THEN
						droprad = SQRT( AreaBlob(i)/pi )/DistMult   ! radius [m]

     					!dropVol = (4/3.0)*pi*r2*droprad		! [m^3] sphere of charge
						!dropVol = pi*droprad*droprad*1.0    	! [m^3] cylinder of charge
						dropVol = (12.*pi*(blob_rc(i)/DistMult)*droprad*droprad)**(1./3) 	! [m^3] torus of charge

						chargeMassArea = PanCharge(i)/(AreaBlob(i)*q)! [#charges/cm^2]
						chargeMassVol = PanCharge(i)/dropVol   	! [C/m^3]
     					!chi = PanCharge(i)*PanCharge(i)/(64*pi*pi*epsil0*surften_in*aa) (sqrt?)
						chi = 58780*PanCharge(i)*PanCharge(i)/SQRT(droprad) ! 3.455e9 otherwise

     					!MTCR = 871080*SQRT( rho_in*rho_in*aa/(36*epsil0*surften_in) ) ! 
						!MTCR = dropVol*rho_in*5.068e24*q/PanCharge(i) ! [#in particle/1 short e-]
						MTCR= 1.216E11*(blob_rc(i)/DistMult)*droprad*droprad/PanCharge(i)

						WRITE(unit=NewDropunit, fmt=11) PanCharge(i), AreaBlob(i), dropVol, &
							chargeMassArea, chargeMassVol, blob_rc(i), blob_zc(i), chi, &
							MTCR, i, counter-1, totalTime
						WRITE(unit=NewDropunit2, fmt=20) PanCharge(i),',', AreaBlob(i),',', dropVol,',', &
							chargeMassArea,',', chargeMassVol,',', blob_rc(i),',', blob_zc(i),',', &
							chi,',', MTCR,',', i,',', counter-1,',', totalTime

     					! mark each blob # once as it passes about a global height
						IF (highPDFpt<GaussPan(i)%midz .AND. i>highblobNum) THEN
							PRINT *, '  -Writing "PDFhigh.dat"'
							highblobNum=highblobNum+1
							OPEN(unit=PDFhighunit, file='PDFhigh.dat', status='old', access='sequential', position='append')
							WRITE(unit=PDFhighunit, fmt=20) PanCharge(i),',', AreaBlob(i),',', dropVol,',', &
								chargeMassArea,',', chargeMassVol,',', blob_rc(i),',', blob_zc(i),',', &
								chi,',', MTCR,',', i,',', counter-1,',', totalTime
							CLOSE (unit=PDFhighunit, status='keep')
						ENDIF

     					! write to a file once as each new blob is created
						IF (z1<i) THEN
							PRINT *, '  -Writing "PDFsingle.dat"'
							OPEN(unit=PDFsingleunit, file='PDFsingle.dat', status='old', &
								access='sequential', position='append', IOSTAT=filestat)
							IF (0==filestat) THEN
								WRITE(unit=PDFsingleunit, fmt=20) PanCharge(i),',', AreaBlob(i),',', dropVol,',', &
									chargeMassArea,',', chargeMassVol,',', blob_rc(i),',', blob_zc(i),',', &
									chi,',', MTCR,',', i,',', counter-1,',', totalTime
								CLOSE (unit=PDFsingleunit, status='keep')
							ELSE
								PRINT *, 'PDFsingle.dat file writing error'
								STOP
							ENDIF
						ENDIF
					ENDIF PrintCharge
				ENDDO
				CLOSE (unit=NewDropunit, status='keep')
				CLOSE (unit=NewDropunit2, status='keep')

				DEALLOCATE(PanCharge); DEALLOCATE(BlobCharge)
			ELSE
				PRINT *, 'MultBlobNumIf file writing error'
				STOP
			ENDIF
		ENDIF MultBlobNumIf

  		!-----------------------------  
  		! run background 1-time potential field
		IF (1==counter .AND. 1==CalcPotBack) THEN
			CALL GetPotGrid(aa, PotCase, aa, aa, aa, 1)
   			!.....
			OPEN(UNIT=Potunit, file=PotCase, status='old', IOSTAT=filestat)
			IF (filestat/=0) THEN ! base background file does not exist
				CALL cpu_time(t1)
				ALLOCATE(meshZ(PotMeshSize))
				ALLOCATE(meshR(PotMeshSize))
				DO mm=1,PotMeshSize
					meshR(mm)=-1.; meshZ(mm)=-1.;
				ENDDO
				PRINT *, '  Creating baseline background electric field in MAT_Greens'

   				!.....
				CALL CalcGridPot(meshR, meshZ, BlobNum, NeedleNum, (minR+1e-4)/DistMult, &
					blob_rc(BlobNum+1)/DistMult, minZ/DistMult, maxZ/DistMult, numGaussPan, &
					wts, GaussPan, 0, NewNumPtsPan, BlobPlus, DistMult)
    			!.....

				CALL cpu_time(t2); timeCount=t2-t1
				WRITE(*, fmt=19) '     Solving for background potentials took ', timeCount, ' seconds.'
				DEALLOCATE(meshZ); DEALLOCATE(meshR)
			ENDIF
		ENDIF

  		!-----------------------------  
  		! run optional arbitrary point potential analysis to see phi plots
		OptionalCalcPot: IF ( 1==CalcPotentials .AND. 0==MOD(counter-1,PotRedoStep) ) THEN
			CALL cpu_time(t1)
			ALLOCATE(meshZ(PotMeshSize))
			ALLOCATE(meshR(PotMeshSize))
			DO mm=1,PotMeshSize
				meshR(mm)=-1.; meshZ(mm)=-1.;
			ENDDO
   			!.....
			CALL CalcGridPot(meshR, meshZ, BlobNum, NeedleNum, (minR+1e-4)/DistMult, &
				blob_rc(BlobNum+1)/DistMult, minZ/DistMult, maxZ/DistMult, numGaussPan, &
				wts, GaussPan, 1, NewNumPtsPan, BlobPlus, DistMult)
   			!.....

			DEALLOCATE(meshZ); DEALLOCATE(meshR)
			CALL cpu_time(t2); timeCount=t2-t1
			WRITE(*, fmt=19) '     Solving for potential field took ', timeCount, ' seconds.'
		ENDIF OptionalCalcPot

		ALLOCATE(rloc(panBlobNum))
		ALLOCATE(zloc(panBlobNum))
		DO i=1,panBlobNum
			rloc(i)=-1.; zloc(i)=-1.
		ENDDO
		i2=0
		DO i=1,BlobNum
			DO j=1,NewNumPtsPan(i)
				i2=i2+1
				rloc(i2) = blob_rc(i) + radBlob(i,j)*cos(thetaBlob(i,j))
				zloc(i2) = blob_zc(i) + radBlob(i,j)*sin(thetaBlob(i,j))
			ENDDO
		ENDDO
  		!-----------------------------
  		! write out last data file
		IF (1==ShowPt_Force .AND. counter>0 .AND. 0==MOD(counter-1, SkipSteps)) THEN
   			! convert data number to string file name
			CALL Int2Str(LSstep, 'PForce', counter-1,'.dat',4)
			PRINT *, '  -Writing ', LSstep
			OPEN(unit=Sigmaunit, file=LSstep, status='replace')
			WRITE(unit=Sigmaunit, fmt='(a)') 'VARIABLES="R [cm]"'
			WRITE(unit=Sigmaunit, fmt='(a)') '"Z [cm]"'
			WRITE(unit=Sigmaunit, fmt='(a)') '"Force [N/m]"'
			WRITE(unit=Sigmaunit, fmt='(a)') '"En [V/m]"'
			WRITE(unit=Sigmaunit, fmt='(a)') '"Panel number"'
			WRITE(unit=Sigmaunit, fmt='(a)') 'ZONE T="main"'

			i2=0
			DO i=1,BlobNum
				DO j=1,NewNumPtsPan(i)
					i2=i2+1
					WRITE(unit=Sigmaunit, fmt='(4e18.6, i8)'), rloc(i2), zloc(i2), &
						forceBlob(i,j), En(i2), i2
				ENDDO
			ENDDO
			CLOSE (unit=Sigmaunit, status='keep')
		ENDIF

		IF (1==RecastDist ) THEN  ! changes matrix back to a cm-based system instead of meters.
			CALL ForceToMeters(GaussPan, numGaussPan, 0);
		ENDIF
  		! prints out angle eta(r,t)
		IF (1==ShowChiAngle .AND. 0==MOD(counter-1, SkipSteps)) THEN
			CALL CalcChiAngle(rloc, zloc, panBlobNum)
		ENDIF

  		! deallocates everything - you're done
		DEALLOCATE(Coords); DEALLOCATE(NumTotNewPanels); DEALLOCATE(phi)
		IF (ALLOCATED(PtsForce)) DEALLOCATE(PtsForce)
		DEALLOCATE(DCenter); DEALLOCATE(swapblobNum); DEALLOCATE(GaussPan)
		DEALLOCATE(En); DEALLOCATE(panelarray)

		IF (ALLOCATED(rloc)) DEALLOCATE(rloc)
		IF (ALLOCATED(zloc)) DEALLOCATE(zloc)
		IF (ALLOCATED(AreaBlob)) DEALLOCATE(AreaBlob)
		IF (ALLOCATED(CumNumGaussPan)) DEALLOCATE(CumNumGaussPan)
		GOTO 999 ! actual program end, just lists formats and errors below

    	!-----------------------------
  		! FORMAT listings
		11 FORMAT(9e14.6, i4, i5, e14.6)
		12 FORMAT(2i5, 3f20.10)
		15 FORMAT(2e14.6, i6, 2e14.6)
		17 FORMAT(a,2e11.4)
		19 FORMAT(a,f8.2,a)
		20 FORMAT(e14.6, a, e14.6, a, e14.6, a, e14.6, a, e14.6, a, e14.6, a, e14.6, a, e14.6, a, e14.6, a, i4, a, i5, a, e14.6)

		999 CALL cpu_time(t2); timeCount=t2-t3
		WRITE(*, fmt=19) 'Entire UMich routine took ', timeCount, ' seconds.'
		PRINT *, '********************'

CONTAINS 
		SUBROUTINE InnerPotSet
   			!****************
   			! Potentially sets up loop to get correct detached droplet potential
  			!****************
			IF (1==ConnectDrops .AND. BlobNum>1) THEN
    			!---------------------------
    			! to calculate the detached blobs' dirichlet potentials
    			!---------------------------
				DO mm=BlobNum+1,BlobPlus
					BlobCharge(mm)=0.
				ENDDO
				IF (ALLOCATED(temp7)) DEALLOCATE(temp7)

				PRINT *, '       temp7(NumCoords)=',NumCoords
				ALLOCATE(temp7(NumCoords*3))
				DO m=1,NumCoords
					temp7(m)=Coords(m,3); temp7(m+NumCoords)=Coords(m,4); temp7(m+2*NumCoords)=Coords(m,5)
				ENDDO

    			!..... 
				CALL CalcDropletPot ( temp7, BlobCharge, PotErr, LastPot, numBlobErr, &
					blob_rc, blob_zc, sizeLevel, NeedleNum, BlobNum, BlobPlus, NumCoords, &
					newdropNum, wts, RecalcPot, GaussPan, numGaussPan, NewNumPtsPan, &
					CumNumGaussPan, swapblobNum, AreaBlob, BlobTimeCharge, z2)
    			!..... 

				DO m=1,NumCoords
					Coords(m,3)=temp7(m)
					Coords(m,4)=temp7(m+NumCoords)
					Coords(m,5)=temp7(m+2*NumCoords)
				ENDDO
				DEALLOCATE(temp7)

				!RecalcPot=.TRUE.
			ENDIF
			WRITE(*, fmt=17) '   The electrode middle edge is at r=', elecEdge

			ALLOCATE(delRCoord(NumCoords))
			DO i=1,NumCoords
				iPlus=Nexti( i, sizeLevel(INT(Coords(i,5))), INT(Coords(i,8)) )
				delRCoord(i)=Coords(iPlus,1)-Coords(i,1) 
			ENDDO
			PRINT *, '  Done adding points around boundary'

   			!----------------------------- 
     		! create spaced points along the spline. Returns NewCoords
			DO mm=1,SmallInitAry
				NewNumPtsPan(mm)=0
			ENDDO
			k=4*LSMaxZ+10
			DO mm=1,k
				NewCoords(mm,1)=-1.; NewCoords(mm,2)=-1.
				NewCoords(mm,3)=-1.; NewCoords(mm,4)=-1.
				NewCoords(mm,5)=-1.; NewCoords(mm,6)=-1.
				NewCoords(mm,7)=-1.; NewCoords(mm,8)=-1.
			ENDDO
			i2=0; panBlobNum=0;

			DO i=1,BlobPlus
				k=sizeLevel(i)

				IF (i>1) THEN
					NumTotNewPanels(i)=NumTotNewPanels(i-1)
				ELSE
					NumTotNewPanels(i)=0

				ENDIF
				i2=i2+sizeLevel(i); j=1

    			! get to starting number
				DO WHILE(INT(Coords(j,5))/=i)
					j=j+1
					IF(INT(Coords(j,5))==i) THEN
						EXIT
					ENDIF
				ENDDO

				IF (j>1) THEN
					j=j-1 ! find ending number
				ENDIF

				DO WHILE (INT(Coords(j,5))/=i)
					j=j+1
					IF (j>=NumCoords) THEN
						PRINT *, 'Shading:BoundaryForceCalc Call curvefit err. i=',i
						STOP
					ENDIF
				ENDDO

				ALLOCATE(temp1(k)); ALLOCATE(temp2(k)); ALLOCATE(temp3(k))
				ALLOCATE(temp4(k,2)); ALLOCATE(temp5(k,3)) 
				DO mm=1,k
					temp1(mm)=-1.; temp2(mm)=-1.; temp3(mm)=-1.; temp4(mm,1)=-1.; temp4(mm,2)=-1. 
					temp5(mm,1)=-1.; temp5(mm,2)=-1.; temp5(mm,3)=-1.;
				ENDDO

				tempRange=0
				DO mm=j,i2
					tempRange=tempRange+1
					temp1(tempRange)=Coords(mm,1); temp2(tempRange)=Coords(mm,2)
					temp3(tempRange)=delRCoord(mm); 
					temp4(tempRange,1)=Coords(mm,3)
					temp4(tempRange,2)=Coords(mm,4)
					temp5(tempRange,1)=Coords(mm,5)
					temp5(tempRange,2)=Coords(mm,6)
					temp5(tempRange,3)=Coords(mm,7)
				ENDDO

    			!-------------
    			! creates 'newcoords', # of points/panel and linear/cubic fit
    			!..... 
				CALL CurveFit ( NumTotNewPanels(i), NewNumPtsPan(i), temp1, &
					temp2, temp3, k, temp4, temp5, INT(Coords(j,8)), i, BlobNum, &
					maxDistStep_equalLS )
    			!..... 
				IF (i<=BlobNum) THEN
					panBlobNum=panBlobNum+NewNumPtsPan(i)
				ENDIF

				DEALLOCATE(temp1); DEALLOCATE(temp2); DEALLOCATE(temp3)
				DEALLOCATE(temp4); DEALLOCATE(temp5)
			ENDDO
			totPanels=NumTotNewPanels(BlobPlus); DEALLOCATE(delRCoord)
			PRINT *, 'Done after CurveFit_rz. # NewCoord nonedge panels=', panBlobNum

			SNC: IF (1==ShowNewCoords) THEN
      			! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatpan]
				PRINT *, '  -Writing "NewCoords.dat" '
				OPEN(unit=NewCoordsunit, file='NewCoords.dat', status='replace')
				WRITE(unit=NewCoordsunit, fmt='(a)') 'VARIABLES="R [cm]"'
				WRITE(unit=NewCoordsunit, fmt='(a)') '"Z [cm]"'
				WRITE(unit=NewCoordsunit, fmt='(a)') '"Pt val"'
				WRITE(unit=NewCoordsunit, fmt='(a)') '"Pt type"'
				WRITE(unit=NewCoordsunit, fmt='(a)') '"Blob#"'
				WRITE(unit=NewCoordsunit, fmt='(a)') '"Blob_rc"'
				WRITE(unit=NewCoordsunit, fmt='(a)') '"Blob_zc"'
				WRITE(unit=NewCoordsunit, fmt='(a)') '"OpenShape"'
				WRITE(unit=NewCoordsunit, fmt='(a)') '"Panel #"'

				DO panelCount=1,totPanels
					WRITE(unit=NewCoordsunit, fmt='(3e14.6, 2i4, 2e14.6, 2i5)') NewCoords(panelCount,1), &
						NewCoords(panelCount,2), NewCoords(panelCount,3), INT(NewCoords(panelCount,4)), &
						INT(NewCoords(panelCount,5)), NewCoords(panelCount,6), NewCoords(panelCount,7), &
						INT(NewCoords(panelCount,8)), panelCount
				ENDDO
				CLOSE (unit=NewCoordsunit, status='keep')
			ENDIF SNC

			IF (ALLOCATED(CumNumGaussPan)) DEALLOCATE(CumNumGaussPan)
			ALLOCATE(CumNumGaussPan(BlobPlus)); 
			DO mm=1,BlobPlus
				CumNumGaussPan(mm)=0
			ENDDO
			CumNumGaussPan(1)=1; z1=1

			IF (ALLOCATED(GaussPan)) DEALLOCATE(GaussPan)
			ALLOCATE(GaussPan(totPanels))
			DO i=1,totPanels ! set entire array to 1.s
				CALL ZeroorCopyPanel(GaussPan(i), 1)
			ENDDO
  		 	!----------------------------- 
   			! generate panels locations, normals and colocation pt. weights
			panelCount=0; rSmoothEdge=0
			PanDefineI: DO i = 1,totPanels
				panelCount=panelCount+1
				i2=Nexti( i, totPanels, 0 )

				IF ( INT(NewCoords(i,5)) == INT(NewCoords(i2,5)) ) THEN  ! middle of blob
					iPlus=i+1
				ELSE
					iPlus=i

					NewNumPtsPan(INT(NewCoords(i,5)))=NewNumPtsPan(INT(NewCoords(i,5)))-1
					panelCount=panelCount-1

					z1=z1+1  ! lists the overall 'GaussPan' beginning panel number at each blob
					IF (z1<=BlobPlus) CumNumGaussPan(z1)=panelCount+1
				ENDIF

				OpenEdge: IF ( iPlus /= i ) THEN
					GaussPan(panelCount)%midPot=NewCoords(i,3) 
					if (NewCoords(i,4)<0.9999) then
						GaussPan(panelCount)%type1=0.0;
					else
						GaussPan(panelCount)%type1=1.0;
					endif
					GaussPan(panelCount)%r0=NewCoords(i,1)
					GaussPan(panelCount)%z0=NewCoords(i,2)
					GaussPan(panelCount)%r1=NewCoords(i+1,1)
					GaussPan(panelCount)%z1=NewCoords(i+1,2)
					GaussPan(panelCount)%midr=0.5*(GaussPan(panelCount)%r0+GaussPan(panelCount)%r1)
					GaussPan(panelCount)%midz=0.5*(GaussPan(panelCount)%z0+GaussPan(panelCount)%z1)

					delR = NewCoords(iPlus,1) - NewCoords(i,1)
					delZ = NewCoords(iPlus,2) - NewCoords(i,2)

					GaussPan(panelCount)%length = SQRT(delR*delR + delZ*delZ)
					GaussPan(panelCount)%rnrm=-delZ/GaussPan(panelCount)%length ! inward normals
					GaussPan(panelCount)%znrm=delR/GaussPan(panelCount)%length
					GaussPan(panelCount)%str=0.

					IF (0==GaussPan(panelCount)%length) THEN
						PRINT *, 'MAT_Greens:BFC PanLen error. i,iPlus,k,panelCount=', i,iPlus,k,panelCount
						STOP 
					ENDIF

     				! Spread the collocation points IN EACH PANEL according to Gaussian quadrature
					PanDefineTo2: DO j = 1,npoints
						GaussPan(panelCount)%rpoints(j)=(delR+2.*GaussPan(panelCount)%r0)/2. + &
							gpts(j)*delR/2.
						GaussPan(panelCount)%zpoints(j)=(delZ+2.*GaussPan(panelCount)%z0)/2. + &
							gpts(j)*delZ/2.
						GaussPan(panelCount)%rscaled(j)=GaussPan(panelCount)%rpoints(j) / &
							GaussPan(panelCount)%midr
					ENDDO PanDefineTo2
				ENDIF OpenEdge

				IF ( NeedleNum==INT(NewCoords(i,5)) .AND. GaussPan(panelCount)%midr>ForceZeroR) THEN
     				! finds where the left boundary of arbitrary forced smoothing in the needle is located. Later run
					rSmoothEdge=rSmoothEdge+1
				ENDIF
			ENDDO PanDefineI  ! for i
			PRINT *, '  Done after full panel definition with panelCount=',panelCount
			numGaussPan=panelCount

   			!-----------------------------
   			! reorders nonneedle 'GaussPan' to be in a counter-clockwise order
			pan=1
			iBlobDo: DO i=1,BlobNum
				k=swapblobNum(i)
				i2=NewNumPtsPan(i)  

				NotNeedIf: IF (k/=NeedleNum) THEN
					ALLOCATE(temp1(i2)); ALLOCATE(temp2(i2))
					ALLOCATE(tempSwapPan(i2)); ALLOCATE(orderedPan(i2))
					DO m=1,i2
						temp1(m)=-1.; temp2(m)=-1.; 
						CALL ZeroorCopyPanel( tempSwapPan(m), 0, GaussPan(CumNumGaussPan(k)+m-1) )
						CALL ZeroorCopyPanel( orderedPan(m), 1 )
					ENDDO

     				! Create temp horizontal reference line
					aa=0.1; PtC(1)=blob_rc(i); PtC(2)=blob_zc(i)
					PtRef(1)=PtC(1)+0.1; PtRef(2)=PtC(2)

					jsizeBlobDo2: DO j=1,i2
						pan=CumNumGaussPan(k)+j-1
						PtA(1)=GaussPan(pan)%r0; PtA(2)=GaussPan(pan)%z0

						ref = SQRT( (PtA(1)-PtC(1))*(PtA(1)-PtC(1)) + (PtA(2)-PtC(2))*(PtA(2)-PtC(2)) )
						c = SQRT( (PtA(1)-PtRef(1))*(PtA(1)-PtRef(1)) + (PtA(2)-PtRef(2))*(PtA(2)-PtRef(2)) )

      					! 0 < acos(x) < pi. Move to correct quadrant
						IF (0==ref) THEN
							ThreePtAngle=0
						ELSE
							ThreePtAngle=( (aa*aa) + (ref*ref) - (c*c) )/(2*aa*ref)
						ENDIF

						Flatcheck: IF ( ABS(ThreePtAngle+1)<eps ) THEN ! too close to -1
							temp1(j)=pi
						ELSEIF ( ABS(ThreePtAngle-1)<eps ) THEN   ! too close to 1
							temp1(j)=0
						ELSE
							IF ( PtA(2)>=PtC(2) ) THEN
								temp1(j)=acos( ThreePtAngle )
							ELSE
								temp1(j)=-acos( ThreePtAngle )
							ENDIF

							IF ( temp1(j) < 0 ) THEN
								temp1(j)=temp1(j)+2.0*pi
							ENDIF
						ENDIF Flatcheck
						temp2(j)=j
					ENDDO jsizeBlobDo2

     				! now that the theta angle is calculated for this blob, sort by degrees
     				!.....
					CALL qsort(temp1, temp2, i2)
					! always sorts CCW, since largest to smallest angle
					DO j=1,i2
						CALL ZeroorCopyPanel( orderedPan(j), 0, tempSwapPan(INT(temp2(j))), 1 )
					ENDDO

     				! improve sort and then reassign points if need be to switch order
     				!.....
     				!CALL SortNN_NoLong( orderedPan, i2 )
					CALL SortNN_NoCross( orderedPan, i2, i )
					!CALL SortWeighted_AngDist( orderedPan, i2 ) 
					CALL FlipPanPtOrder(orderedPan, i2, PtC, maxR)
     				!.....

					DO j=1,i2 ! copy change back to the master list
						CALL ZeroorCopyPanel( GaussPan(CumNumGaussPan(k)+j-1), 0, orderedPan(j) )
					ENDDO
					DEALLOCATE(tempSwapPan); DEALLOCATE(orderedPan)
					DEALLOCATE(temp1); DEALLOCATE(temp2)
				ENDIF NotNeedIf
			ENDDO iBlobDo

			IF (1==RecastDist ) THEN ! changes matrix to a meters-based system instead of cm.
				CALL ForceToMeters(GaussPan, numGaussPan, 1);
			ENDIF

			CALL cpu_time(t1)

			IF (ALLOCATED(panelarray)) DEALLOCATE(panelarray)
			k=numGaussPan*numGaussPan
			ALLOCATE(panelarray(numGaussPan*numGaussPan)); panelarray(:)=0.

			IF (0==AformDouble) THEN
    			! Create matrix 'A' of elliptic integrals. Also solves for mu and nu as 'sigma', 
    			!   the matricies related to En and Phi. *Major* section that returns needed 
    			!   sigma value.

				!..... 
				IF (3==MatSolveMethod) THEN  ! Calculate panel array matrix
					CALL panelconstantswrap(numGaussPan, GaussPan, panelarray, wts, 1)
				ELSE

					CALL panelconstantswrap(numGaussPan, GaussPan, panelarray, wts, 0)
				ENDIF
   				!..... 
			ELSEIF (1==AformDouble) THEN
    			! Create matrix 'A' of elliptic integrals. Also solves for dphi, phi as 'sigma', 
   				!   Gets En directly
				ALLOCATE(b(numGaussPan))
				DO mm=1,numGaussPan
					b(mm)=-1.
				ENDDO
    			!..... 
				IF (3==MatSolveMethod) THEN  ! Calculate panel array matrix
					CALL adoublesetupwrap(numGaussPan, GaussPan, b, panelarray, wts, 1)
				ELSE
					CALL adoublesetupwrap(numGaussPan, GaussPan, b, panelarray, wts, 0)
				ENDIF
    			!..... 
			ENDIF

			CALL cpu_time(t2); timeCount=t2-t1
			IF (3==MatSolveMethod) THEN  ! Calculate panel array matrix
				WRITE(*, fmt=19) '   Setting up and inverting the A matrix took ', timeCount, ' seconds.'
			ELSE
				WRITE(*, fmt=19) '   Setting up the A matrix took ', timeCount, ' seconds.'
			ENDIF

   			! if preconditioning matrix
   			!---------------------------
			IF (1==UsePreCond .AND. 1==MatSolveMethod .AND. 0==AformDouble) THEN
				ALLOCATE(Minv(numGaussPan,numGaussPan)); ALLOCATE(b(numGaussPan))

				DO i = 1,numGaussPan
					Minv(i,i)=1.0/panelarray( (i-1)*numGaussPan+i )
					b(i)=Minv(i,i)*GaussPan(i)%midPot/chargeconstant
					DO j = 1,numGaussPan
						IF (i/=j) Minv(i,j)=0.0;
						pan = (i-1)*numGaussPan+j
						panelarray( pan ) = Minv(i,i)*panelarray( pan )
					ENDDO
				ENDDO

				DEALLOCATE(Minv)

				CALL cpu_time(t1); timeCount=t1-t2
				WRITE(*, fmt=19) '   Jacobinan preconditioning the A matrix took ', timeCount, ' seconds.'
			ELSEIF (1==MatSolveMethod .AND. 0==AformDouble) THEN
				PRINT *, '    No matrix preconditioner used.'
				ALLOCATE(b(numGaussPan));
				DO i = 1,numGaussPan
					b(i)=GaussPan(i)%midPot/chargeconstant
				ENDDO
				CALL cpu_time(t1)
			ELSE
				CALL cpu_time(t1)
			ENDIF

   			!-----------------------------
     		! solves for {X} in Ax=b
			WRITE(*, fmt=18), '   Done up to matrix Ax=b with ', numGaussPan, ' panels.'

			SELECT CASE( MatSolveMethod )
			CASE (1)
      			!---- gmres
				ALLOCATE(temp1(numGaussPan)); 
				DO i=1,numGaussPan
					temp1(i)=0.0
				ENDDO

				CALL gmresIter(temp1, panelarray, b, numGaussPan)

				DO i=1,numGaussPan
					GaussPan(i)%str=temp1(i)
				ENDDO
				DEALLOCATE(b); DEALLOCATE(temp1)
			CASE (2)  
       			!---- Gaussian inversion with partial pivoting
    			!ALLOCATE( indxLegs(numGaussPan) ); indxLegs(:)=-1
    			!ALLOCATE( sigmaGauss(numGaussPan) ); sigmaGauss(:)=0.
    			!ALLOCATE( sigma(numGaussPan)); sigma(:)=-1.
    			!CALL LEGS (A, numGaussPan, RHStemp, sigmaGauss, indxLegs)
    			!PRINT *, '    Done with inversion'

    			!DO j=1,numGaussPan ! correct the shuffled sigma locations
    			!sigma(indxLegs(j))=sigmaGauss(j)
    			!ENDDO
    			!DEALLOCATE(indxLegs); DEALLOCATE(sigmaGauss)
				PRINT *, '  Err - no code for full Gaussian inversion'
			CASE (3)
    			!---- Gaussian inversion WITHOUT partial pivoting
				CALL panelsolverealwrap(numGaussPan, panelarray, GaussPan);

				PRINT *, '  Done after regular Gaussian inversion'
			CASE DEFAULT
				STOP 'Wrong method for MatSolveMethod'
			END SELECT
			CALL cpu_time(t2); timeCount=t2-t1
			WRITE(*, fmt='(a,e10.3)') '     sigma(4)=', GaussPan(4)%str
			WRITE(*, fmt=19) '     Solving for sigma took ', timeCount, ' seconds.'

			IF ( GaussPan(1)%str>1e30 .AND.  1==MatSolveMethod ) THEN
				PRINT *, "WARNING: Iterative approach got NaN error. Re-running using full matrix inversion."
				CALL panelsolverealwrap(numGaussPan, panelarray, GaussPan);

				CALL cpu_time(t1); timeCount=t1-t2
				PRINT *, '    sigma(1)=', GaussPan(1)%str
				WRITE(*, fmt=19) '     Solving for sigma the 2nd time took ', timeCount, ' seconds.'
			ENDIF
			IF ( GaussPan(1)%str>1e30 ) THEN
				PRINT *, "Still getting a NaN error. Killing process"
				STOP
			ENDIF

			IF (1==RecastDist ) THEN  ! changes matrix back to a cm-based system instead of meters.
				CALL ForceToMeters(GaussPan, numGaussPan, 0);
			ENDIF

			GOTO 998 ! actual program end, just lists formats and errors below

     		!-----------------------------
   			! FORMAT listings
			17 FORMAT(a,2e11.4)
			18 FORMAT(a,i4,a)
			19 FORMAT(a,f8.2,a)
		998 END SUBROUTINE InnerPotSet
	END SUBROUTINE BoundaryForceCalc

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
! Primary author: Anton VanderWyst, University of Michigan 
!
! Version history:
!   01 Jul 2004 begun  
! 07 Oct   v3 removing the trapezoid component. Only Gaussian used
!  20 Jun 2005 using Jerry's automatic number generation code

	SUBROUTINE gaussvals(points, weights)
  		! incoming variables - only global number of g. points
  		! outgoing variables
		REAL*8, DIMENSION(npoints), INTENT(OUT) :: points, weights

  		! start regular variables
		INTEGER      :: i, j, ind, counter2
		REAL*8      :: error, oldval
		REAL*8, DIMENSION(npoints)  :: beta
		REAL*8, DIMENSION(npoints*npoints) :: z, q, r

		CONTINUE 
		error=1.; oldval=2.; z=0.

		DO i=2,npoints
			beta(i-1)=0.5/sqrt(1.0-1.0/(4.0*(i-1)*(i-1)))
		ENDDO

		DO i=1,npoints-1
			z((i-1)*npoints+i+1)=beta(i)
			z(i*npoints+i)=beta(i)
		ENDDO

		counter2=0
		DoErr: DO WHILE (error > 1e-14)  
			counter2=counter2+1
   			!..... 
			CALL mgs(z, q, r)   ! Compute the QR factorization of z
  			!..... 

			DO i=1,npoints    ! Multiply a and q to get z
				DO j=1,npoints
					z((i-1)*npoints+j)=0.0

					IF (i==1) THEN
						z((i-1)*npoints+j) = z((i-1)*npoints+j)+ beta(i)*q(i*npoints+j)
					ELSEIF (i==npoints) THEN
						z((i-1)*npoints+j) = z((i-1)*npoints+j) + beta(i-1)*q((i-2)*npoints+j)
					ELSE
						z((i-1)*npoints+j) = z((i-1)*npoints+j) + beta(i-1)*q((i-2)*npoints+j) + &
							beta(i)*q(i*npoints+j)
					ENDIF
				ENDDO
			ENDDO

   			! Compute change in first eigenvalue
			error=ABS((r(1)-oldval)/oldval)
			oldval=r(1)
		ENDDO DoErr

  		! Assign points from the diagonal of r (the eigenvalues)
		ind=1
		DO i=1,npoints,2 
			points(ind) = -r((i-1)*npoints+i)
			points(npoints-ind+1) = r((i-1)*npoints+i)
			ind = ind+1
		ENDDO  

  		! Assign weights from the square of the non-zero leading values in the 
  		!   columns of q (the eigenvectors)
		ind=1
		DO i=1,npoints
			IF (q(i) /= 0.0) THEN
				weights(ind)=q(i)*q(i)
				weights(npoints-ind+1)=q(i)*q(i)
				ind = ind+1
			ENDIF
		ENDDO

  		! Take care of the odd numbered case
  		!sum=0.0
		IF (MOD(npoints,2) == 1) THEN
			PRINT *, 'Err: AxiSym:gaussvals odd num'
			STOP
			!points(npoints/2)=0.0

			!DO i=1,npoints/2
			!sum = sum + weights(i)
			!ENDDO
			!weights(npoints/2)=2.0-2.0*sum
		ENDIF
	END SUBROUTINE gaussvals

! ***************
! SUBROUTINE linearvals
! ***************
! DESCRIPTION: Calculates the linear points and weights along panel
!
! INPUT: # points / panel, panel beginning and end points
!
! OUTPUT: location of linear points along panel
!
! CALLING PROGRAM: MAT_Greens

	SUBROUTINE linearvals(points, weights)
  		! incoming variables - only global number of points
  		! outgoing variables
		REAL*8, DIMENSION(npoints), INTENT(OUT) :: points, weights

  		! start regular variables
		INTEGER		:: i
		REAL*8		:: sizeStep
		CONTINUE 

		sizeStep=2./(npoints+1)
		DO i=1,npoints
			points(i)=i*sizeStep - 1
			weights(i)=1./npoints
		ENDDO
	END SUBROUTINE linearvals

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
		REAL*8, DIMENSION(npoints*npoints), INTENT(in)  :: a

  ! outgoing variables
		REAL*8, DIMENSION(npoints*npoints), INTENT(OUT) :: q, r

  ! start regular variables
		INTEGER        :: i, j, k, n
		REAL*8, DIMENSION(npoints*npoints)  :: v

		CONTINUE 
		n=npoints; v=a; r=0;

		ido: DO i=1,n
			r((i-1)*n+i)=0
			DO k=1,n
				r((i-1)*n+i) = r((i-1)*n+i) + v((k-1)*n+i)*v((k-1)*n+i)
			ENDDO
			r((i-1)*n+i)=sqrt(r((i-1)*n+i))  ! Calculate rjj

			DO k=1,n     ! Calculate q
				q((k-1)*n+i)=v((k-1)*n+i)/r((i-1)*n+i)
			ENDDO

			DO j=i+1,n
				DO k=1,n    ! Do the dot product...
					r((i-1)*n+j) = r((i-1)*n+j) + q((k-1)*n+i)*v((k-1)*n+j)
				ENDDO

				DO k=1,n    ! Calculate v
					v((k-1)*n+j) = v((k-1)*n+j) - r((i-1)*n+j)*q((k-1)*n+i)
				ENDDO
			ENDDO
		ENDDO iDo
	END SUBROUTINE mgs

! ***************
! Subroutine ShapeGrow
! ***************
! DESCRIPTION: for the Sussman level set visual file, finds and sorts the number and vertex members 
!   of each blob
!
! INPUT: visual data
!
! OUTPUT: members of each shape, number of shapes, ordered surface around
	SUBROUTINE ShapeGrow( SortedBlobPts, Countblob, NeedleBlob, sizeNeedleBlob, NumPt, &
			LSMaxX, LSMaxY, xLeft, yBottom, delX, delY, is_symmetric1, minX ) 
  		! incoming variables
		INTEGER, INTENT(in)   :: LSMaxX, LSMaxY, is_symmetric1
		REAL*8, INTENT(in)   :: xLeft, yBottom, delX, delY

  		! in/outgoing variables
		INTEGER, INTENT(inout)  :: NumPt, Countblob, NeedleBlob
		INTEGER, DIMENSION(SmallInitAry),INTENT(inout)  :: sizeNeedleBlob
		REAL*8, INTENT(inout)  :: minX
		REAL*8, DIMENSION(InitArySize,8), INTENT(inout) :: SortedBlobPts
       		! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatyesno]

  		! program variables
		INTEGER      :: i, j, k, aa, bb, cc, tempSize, CountUp, &
			Crossunit, Sortunit, startCase, prevNumPt, tempMaxSize, tempPrevStart, &
			Countblobbump
		INTEGER, DIMENSION(:), ALLOCATABLE :: InPts, temp3
		REAL*8, DIMENSION(InitArySize,4)   :: tempInPts, swapInPts
   			! [xpt, ypt, blob#, particle#]
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: temp1, temp2

		CONTINUE

  		! variable setup............
		Crossunit=12; Sortunit=13

  		! Countblob=number of discrete blobs. NumPt=total number of points in all blobs
		Countblob=0; NumPt=0; 
		Countblobbump=0
		DO i=1,InitArySize
			tempInPts(i,1)=-1.; tempInPts(i,2)=-1.
			tempInPts(i,3)=-1.; tempInPts(i,4)=-1.
			swapInPts(i,1)=0.; swapInPts(i,2)=0.
			swapInPts(i,3)=0.; swapInPts(i,4)=0.
		ENDDO

  		! begins computation............
  		! groups each point into a blob family
		PRINT *, '  Entering surface determination'

		minX=1e6
		DO i=StartNum,LSMaxX+StopNum-1
			DO j=LSMaxY+StopNum,StartNum,-1
				IF ( phi(i,j) >= 0 .AND. -1==insideMarker(i,j) ) THEN
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
						DO bb=aa+1,NumPt
							IF ( tempInPts(aa,1)==tempInPts(bb,1) .AND. &
								tempInPts(aa,2)==tempInPts(bb,2) ) THEN
								IF (bb<NumPt) THEN 
									DO cc=bb+1,NumPt
										tempInPts(cc-1,1)=tempInPts(cc,1)
										tempInPts(cc-1,2)=tempInPts(cc,2)
										tempInPts(cc-1,3)=tempInPts(cc,3)
										tempInPts(cc-1,4)=tempInPts(cc,4)
									ENDDO
								ELSE
									tempInPts(NumPt,1)=tempInPts(NumPt,1)
									tempInPts(NumPt,2)=tempInPts(NumPt,2)
									tempInPts(NumPt,3)=tempInPts(NumPt,3)
									tempInPts(NumPt,4)=tempInPts(NumPt,4)
								ENDIF
								NumPt=NumPt-1
								PRINT *, '    Remove Coord pt; aa,bb=',aa,bb
							ENDIF
						ENDDO
					ENDDO RmCoordDo
					NotCountIf: IF (NumPt<prevNumPt+3) THEN
						Countblob=Countblob-1
						NumPt=prevNumPt
					ELSE
      					! create mirror image of all blobs if symmetric
						symDo: IF ( 1 == is_symmetric1 ) THEN
							IF (Countblob>1) THEN
								PRINT *, '  Multiple blobs. Countblob=', Countblob
							ENDIF
							tempSize=NumPt-prevNumPt

							DO aa=prevNumPt+1,NumPt ! find x minimum
								IF ( tempInPts(aa,1)<minX ) THEN
									minX=tempInPts(aa,1)
								ENDIF
							ENDDO
							DO aa=prevNumPt+1,NumPt
								swapInPts(aa,1)=tempInPts(aa,1)
								swapInPts(aa,2)=tempInPts(aa,2)
								swapInPts(aa,3)=tempInPts(aa,3)
								swapInPts(aa,4)=tempInPts(aa,4)
							ENDDO
							k=0; tempMaxSize=tempSize*2+prevNumPt; tempPrevStart=prevNumPt+1
							DO aa=tempPrevStart,tempMaxSize ! now add mirrored points, {x,y, blob#, particle#}
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
							PRINT *, '  NumPt=',NumPt
						ENDIF symDo
					ENDIF NotCountIf
				ENDIF 
			ENDDO
		ENDDO
		PRINT *, '  After SurfaceTrack, there are ',NumPt,' coordinates.'
		Countblob=Countblob+Countblobbump

  		! =begin blob1, begin blob2
		ALLOCATE(InPts(Countblob)); 
		DO i=1,Countblob
			InPts(i)=-1
		ENDDO
		CountUp=1
		DO i=1,NumPt
			IF ( CountUp==INT(tempInPts(i,3)) ) THEN
				InPts(CountUp)=i
				CountUp=CountUp+1
			ENDIF
		ENDDO
		CountUp=CountUp-1
		PRINT *, '  Done grouping points into blobs. NumBlob=', CountUp
		IF (0==NumPt) THEN
			PRINT *, 'Not recognizing shape, terminating'
			STOP
		ENDIF

  		! ***************
		IF (1==ShowNewCross) THEN
			PRINT *, '  -Writing "BlobCross.dat" '
			OPEN(unit=Crossunit, file='BlobCross.dat', status='replace', err=990)
			WRITE(unit=Crossunit, fmt='(a)', err=992) 'VARIABLES="X [cm]"'
			WRITE(unit=Crossunit, fmt='(a)', err=992) '"Y [cm]"'
			WRITE(unit=Crossunit, fmt='(a)', err=992) '"BlobNum"'
			WRITE(unit=Crossunit, fmt='(a)', err=992) '"ParticleNum"'
			DO i=1,NumPt
				WRITE(unit=Crossunit, fmt='(2e13.5, 2i4)') tempInPts(i,1), tempInPts(i,2), &
					INT(tempInPts(i,3)), INT(tempInPts(i,4))
			ENDDO
			CLOSE (unit=Crossunit, status='keep', err=993)
		ENDIF

		ALLOCATE(temp1(NumPt,8)); ALLOCATE(temp2(NumPt,4)); ALLOCATE(temp3(CountBlob))
		DO i=1,NumPt
			temp1(i,1)=SortedBlobPts(i,1); temp1(i,2)=SortedBlobPts(i,2)
			temp1(i,3)=SortedBlobPts(i,3); temp1(i,4)=SortedBlobPts(i,4)
			temp1(i,5)=SortedBlobPts(i,5); temp1(i,6)=SortedBlobPts(i,6)
			temp1(i,7)=SortedBlobPts(i,7); temp1(i,8)=SortedBlobPts(i,8)
			temp2(i,1)=tempInPts(i,1); temp2(i,2)=tempInPts(i,2)
			temp2(i,3)=tempInPts(i,3); temp2(i,4)=tempInPts(i,4)
		ENDDO
		DO i=1,CountBlob
			temp3(i)=sizeNeedleBlob(i)
		ENDDO
  		!.....
		CALL UpdateBlobList( temp1, NeedleBlob, temp3, InPts, temp2, Countblob, NumPt )
  		!.....

		DO i=1,NumPt
			SortedBlobPts(i,1)=temp1(i,1); SortedBlobPts(i,2)=temp1(i,2)
			SortedBlobPts(i,3)=temp1(i,3); SortedBlobPts(i,4)=temp1(i,4)
			SortedBlobPts(i,5)=temp1(i,5); SortedBlobPts(i,6)=temp1(i,6)
			SortedBlobPts(i,7)=temp1(i,7); SortedBlobPts(i,8)=temp1(i,8)
		ENDDO
		DO i=1,CountBlob
			sizeNeedleBlob(i)=temp3(i)
		ENDDO
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

		PRINT *, 'Done with "ShapeGrow" in <MAT_Greens>'

	END SUBROUTINE ShapeGrow

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
		INTEGER, INTENT(in)   :: LSMaxX, LSMaxY, WhichI, WhichJ, NumBlob 

		CONTINUE
		insideMarker(WhichI,WhichJ)=NumBlob

		IF (WhichJ<LSMaxY) THEN
			IF ( phi(WhichI,WhichJ+1)>=0 .AND. -1==insideMarker(WhichI,WhichJ+1) ) THEN ! up liquid+
				CALL BlobSeparate( LSMaxX, LSMaxY, NumBlob, WhichI, WhichJ+1 )
			ENDIF
		ENDIF
		IF ( WhichI<LSMaxX) THEN
			IF ( phi(WhichI+1,WhichJ)>=0 .AND. -1==insideMarker(WhichI+1,WhichJ) ) THEN  ! right liquid+
				CALL BlobSeparate( LSMaxX, LSMaxY, NumBlob, WhichI+1, WhichJ )
			ENDIF
		ENDIF
		IF (WhichJ>1) THEN
			IF ( phi(WhichI,WhichJ-1)>=0 .AND. -1==insideMarker(WhichI,WhichJ-1) ) THEN ! down liquid+
				CALL BlobSeparate( LSMaxX, LSMaxY, NumBlob, WhichI, WhichJ-1 )
			ENDIF
		ENDIF
		IF (WhichI>1) THEN
			IF (phi(WhichI-1,WhichJ)>=0 .AND. -1==insideMarker(WhichI-1,WhichJ)) THEN ! left liquid+
				CALL BlobSeparate( LSMaxX, LSMaxY, NumBlob, WhichI-1, WhichJ )
			ENDIF
		ENDIF
	END SUBROUTINE BlobSeparate

! ***************
! SUBROUTINE SurfaceTrack
! ***************
!
! DESCRIPTION: tracks the blob from the left edge around the surface. Simplifies later
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
! CALLING ROUTINE: MAT_Greens:ShapeGrow

	RECURSIVE SUBROUTINE SurfaceTrack ( tempIn, NumPt, NumBlob, WhichI, &
			WhichJ, prevCase, delX, delY, xLeft, yBottom, LSMaxX, LSMaxY )
  		! incoming variables
		INTEGER, INTENT(in)   :: WhichI, WhichJ, NumBlob, LSMaxX, LSMaxY
		REAL*8, INTENT(in) :: delX, delY, xLeft, yBottom

  		! outgoing variable
		INTEGER, INTENT(inout)  :: NumPt, prevCase
		REAL*8, DIMENSION(InitArySize,4), INTENT(inout)  :: tempIn
    		! {x,y, blob#, particle#}

  		! internal variables
		LOGICAL      :: onEdge, negPhi, diagDir, LoopChoke
		INTEGER      :: i, iNext, jNext, caseplus2i, caseplus2j, &
			countLoop, plusPt, minusPt
		CONTINUE

		onEdge=.FALSE.; alreadyChecked(WhichI, WhichJ)=1
  		! check to make sure at least one of the 8points is phi<0 (i.e. you are on the boundary)
		DO i=1,8
			CALL CaseLoop ( iNext, jNext, caseplus2i, caseplus2j, i, &
				WhichI, WhichJ, LSMaxX, LSMaxY )
   			!..... 

			IF (phi(iNext, jNext)<=0) THEN
				onEdge=.TRUE.
				EXIT
			ENDIF
		ENDDO

		negPhi=.TRUE.; countLoop=0
		negLoop: DO WHILE (negPhi .AND. onEdge)
			countLoop=countLoop+1
			jMax: IF (iNext<LSMaxX-FLOOR(LSMaxX/200.)+1 .AND. jNext<LSMaxY-FLOOR(LSMaxY/200.)+1) THEN
				IF (1==mod(prevCase,2)) THEN
					diagDir=.TRUE.
				ELSE
					diagDir=.FALSE.
				ENDIF
    			!..... 
				CALL CaseLoop ( iNext, jNext, caseplus2i, caseplus2j, prevCase+1, &
					WhichI, WhichJ, LSMaxX, LSMaxY )
    			!..... 
				IF (NumPt>InitArySize-4) THEN
					PRINT *, 'Err: MAT_Greens:SurfaceTrack array size bigger than global:InitArySize'
					STOP
				ENDIF
				MoveMark: IF (1==WhichJ .AND. 7==iNext) THEN ! catches weird SW corner squish
					NumPt=NumPt+1
					tempIn(NumPt,1)=0
					tempIn(NumPt,2)=tempIn(NumPt-1,2)
					tempIn(NumPt,3)=1.0*NumBlob
					tempIn(NumPt,4)=1.0*NumPt

					prevCase = CaseNumShift(prevCase,-1)
				ELSEIF ( phi(iNext, jNext)<0 .AND. .NOT. (diagDir) .AND. &
					phi(WhichI,WhichJ)/=0 ) THEN ! * MARK *
					NumPt=NumPt+1

					IF ( 9==prevCase+1 .OR. 5==prevCase+1 ) THEN  ! if up/down
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
									alreadyChecked(WhichI,WhichJ)=0 ! allow one come-back node landing
									LoopChoke=.FALSE.
								ENDIF
							ENDDO
						ENDIF
					ELSEIF ( 3==prevCase+1 .OR. 7==prevCase+1 ) THEN ! if left/right
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
									alreadyChecked(WhichI,WhichJ)=0 ! allow one come-back node landing
									LoopChoke=.FALSE.
								ENDIF
							ENDDO
						ENDIF
					ELSE
						PRINT *, 'MAT_Greens:SurfaceTrack mark location err.'
						STOP
					ENDIF
					tempIn(NumPt,3)=1.0*NumBlob
					tempIn(NumPt,4)=1.0*NumPt

					prevCase = CaseNumShift(prevCase,1)
				ELSEIF (0==phi(iNext,jNext) .AND. alreadyChecked(iNext,jNext) /=1) THEN  
     				! catching unusual case where no marking would occur
					NumPt=NumPt+1

					tempIn(NumPt,1)=(iNext-1)*delX+xLeft
					tempIn(NumPt,2)=(jNext-1)*delY+yBottom
					tempIn(NumPt,3)=1.0*NumBlob
					tempIn(NumPt,4)=1.0*NumPt

					prevCase = CaseNumShift(prevCase,-3)

     				!..... 
					CALL SurfaceTrack ( tempIn, NumPt, NumBlob, iNext, jNext, prevCase, &
						delX, delY, xLeft, yBottom, LSMaxX, LSMaxY )
				ELSEIF ( (phi(iNext, jNext)>=0 .AND. alreadyChecked(iNext,jNext) /=1) .AND. &
					((phi(caseplus2i, caseplus2j)>=0 .AND. (diagDir)) .OR. .NOT. (diagDir)) ) THEN ! * MOVE *
					negPhi=.FALSE.
					prevCase = CaseNumShift(prevCase,-3)

     				!..... 
					CALL SurfaceTrack ( tempIn, NumPt, NumBlob, iNext, jNext, prevCase, &
						delX, delY, xLeft, yBottom, LSMaxX, LSMaxY )
				ELSE
					prevCase = CaseNumShift(prevCase,1)
				ENDIF MoveMark
			ELSE
				onEdge=.FALSE.
			ENDIF jMax
			IF (countLoop>8) THEN ! you're at the end of the shape
				EXIT
			ENDIF
		ENDDO negLoop
	END SUBROUTINE SurfaceTrack

! ***************
! FUNCTION LinFrac
! ***************
!
! DESCRIPTION: Gives the linear interpolation between two points where the levelset =0
!
! INPUT: levelset value at 2 points, position of 2 points
!
! OUTPUT: new intermediate crossing location
	FUNCTION LinFrac(LS1, LS2, Pos1, Pos2)
  		! incoming variables
		REAL*8, INTENT(in) :: LS1, LS2, Pos1, Pos2

  		! outgoing variable
		REAL*8     :: LinFrac

  		! program variables
		REAL*8    :: m, b
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
		INTEGER    :: CaseNumShift

  		! internal variables
		INTEGER    :: tempNum
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
		INTEGER, INTENT(in)  :: iPt, jPt, caseFrom, xMax, yMax

  		! outgoing variables
		INTEGER, INTENT(out) :: iNext, jNext, caseplus2i, caseplus2j
		CONTINUE

		SELECT CASE (caseFrom) 
		CASE (1) ! up
			iNext=iPt; jNext=jPt+1
			caseplus2i=iPt+1; caseplus2j=jPt+1
		CASE (2) ! up-right
			iNext=iPt+1; jNext=jPt+1
			caseplus2i=iPt+1; caseplus2j=jPt
		CASE (3) ! right
			iNext=iPt+1; jNext=jPt
			caseplus2i=iPt+1; caseplus2j=jPt-1
		CASE (4) ! down-right
			iNext=iPt+1; jNext=jPt-1
			caseplus2i=iPt; caseplus2j=jPt-1
		CASE (5) ! down
			iNext=iPt; jNext=jPt-1
			caseplus2i=iPt-1; caseplus2j=jPt-1
		CASE (6) ! down-left
			iNext=iPt-1; jNext=jPt-1
			caseplus2i=iPt-1; caseplus2j=jPt
		CASE (7) ! left
			iNext=iPt-1; jNext=jPt
			caseplus2i=iPt; caseplus2j=jPt+1
		CASE (8) ! left-up
			iNext=iPt-1; jNext=jPt+1
			caseplus2i=iPt; caseplus2j=jPt+1
		CASE (9) ! wrap-around, repeated up
			iNext=iPt; jNext=jPt+1
			caseplus2i=iPt+1; caseplus2j=jPt+1
		CASE (10) ! wrap-around, repeated up-right
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
! SUBROUTINE UpdateBlobList
! ***************
!
! DESCRIPTION: Sorts points in each blob clockwise. Finds circle center, sorts on decreasing
!   angles from an arbitrary horizontal point
!
! INPUT: array of n {x,y} pts in clump N, example list of clump pts z{3:3+n}
!
! OUTPUT: sorted clumped array

	RECURSIVE SUBROUTINE UpdateBlobList ( SortedBlobPts, NeedleBlob1, sizeBlob, InPts, &
			tempInPts, Countblob, NumPts )
  		! incoming variables
		INTEGER, INTENT(in)        :: Countblob, NumPts
		INTEGER, DIMENSION(Countblob), INTENT(in)  :: InPts
		REAL*8, DIMENSION(NumPts,4), INTENT(in)   :: tempInPts

  		! outgoing variable
		INTEGER, INTENT(inout)       :: NeedleBlob1
		INTEGER, DIMENSION(Countblob), INTENT(inout) :: sizeBlob
		REAL*8,DIMENSION(NumPts,8),INTENT(inout)   :: SortedBlobPts
        	! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, type panel]

  		! program variables
		INTEGER           :: i, k, CurrentPt, ptSum, flagOnAxis
		REAL*8, DIMENSION(2)       :: PtC

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
			IF(ABS(tempInPts(InPts(i),1))-eps<=0.) THEN ! an for onaxis detached blob 
				flagOnAxis=1
			ELSE
				flagOnAxis=0
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

				SortedBlobPts(CurrentPt,5)=tempInPts(k+ptSum,3)
				SortedBlobPts(CurrentPt,6)=PtC(1)
				SortedBlobPts(CurrentPt,7)=PtC(2)

				IF ( sizeBlob(i)<0.2*NumPts) THEN 
    			 	! unless global:ConnectDrops is used, default to a zero N. droplet cond.
					SortedBlobPts(CurrentPt,3)=0. 
					SortedBlobPts(CurrentPt,4)=1. !** == 1 Neumann flux **
				ELSEIF ( sizeBlob(i)>=0.2*NumPts ) THEN
     				! define the fluid coordinates with rarray, zarray and panels 
                	!   with potential -2.0 and dirichlet boundary conditions(==0)
					SortedBlobPts(CurrentPt,3)=fluidPot 
					SortedBlobPts(CurrentPt,4)=0.

					NeedleBlob1=i
				ELSE
					PRINT *, 'MAT_Greens:UpdateBlobList sizeBlob err'
					STOP
				ENDIF

				IF (1==flagOnAxis .OR. sizeBlob(i)>0.2*NumPts) THEN
					SortedBlobPts(CurrentPt,8)=1.
				ELSE
					SortedBlobPts(CurrentPt,8)=0.
				ENDIF
			ENDDO
		ENDDO iDo
	END SUBROUTINE UpdateBlobList

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
		INTEGER, INTENT(in)    :: code, extenLen
		CHARACTER(LEN=6), INTENT(in) :: msg
		CHARACTER(LEN=extenLen), INTENT(in) :: exten

  		! outgoing
		CHARACTER(LEN=15), INTENT(out)  :: StrText

  		! internal
		INTEGER, parameter    :: cbuff_size=5 ! size of code string
		CHARACTER( len= cbuff_size)  :: cbuff        ! construct string code
		CHARACTER( len= 1), PARAMETER  :: blank = ' '    
		CHARACTER( len= 1), PARAMETER  :: minus = '-'  ! if negative code
		INTEGER :: idigit, abs_code, lenStr             ! loop thru digits

		CONTINUE

		code_0: IF ( code /= 0 ) THEN                  ! special case zero
   			!  process positive numbers only
			abs_code = abs(code)                     ! convert integer to ascii

			cbuff = blank
			idigit = cbuff_size

			each_digit: DO WHILE( abs_code > 0)      ! loop thru digits
				cbuff( idigit: idigit ) = achar( mod( abs_code, 10) + iachar( '0'))
				abs_code = abs_code / 10           ! next magnitude
				idigit = idigit - 1            ! next digit
			ENDDO each_digit                      ! loop thru digits

   !  if code was negative, add minus sign
			IF( code < 0 ) THEN
				cbuff( idigit: idigit) = minus       ! negative --> minus
			ENDIF

			lenStr=len_trim(ADJUSTL(cbuff))

			SELECT CASE (lenStr)       ! assemble msg and code
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
		ELSE                                 ! do zero here
			StrText =trim(msg) // '00000' // exten      ! append zero
		ENDIF code_0                                 ! code 0 or not
	END SUBROUTINE Int2Str

! ***************
! SUBROUTINE TrackDrop
! ***************
! DESCRIPTION: tracks droplet creation/destruction
!
! INPUT: current droplet centers
!
! OUTPUT: the updated blob numbers, so it can track directions 
!
! ALGORITHM: compares old droplet (r,z) center from file to new (rc,zc) location.
!   Also assigns drops same charge as the previous timestep. New droplets are 
!   given a new blob number, lost ones leave an unused blob number.

	SUBROUTINE TrackDrop(newdropNum, DropLinks, ChargeDrop, centDist, AvgCenMove, &
			rc, zc, Areab, blobNum, NewNeedNum)
  		! incoming
		INTEGER, INTENT(in)      :: blobNum, NewNeedNum
		REAL*8, DIMENSION(blobNum), INTENT(in)  :: Areab

  		! outgoing
		INTEGER, INTENT(inout)     :: newdropNum
		INTEGER,DIMENSION(blobNum),INTENT(inout):: DropLinks
		REAL*8, INTENT(inout)     :: AvgCenMove
		REAL*8,DIMENSION(blobNum),INTENT(inout) :: ChargeDrop, rc, zc, centDist

  		! internal
		CHARACTER   :: chara
		CHARACTER(LEN=15) :: DLink2, cSurfCharge2
		INTEGER    :: BlobCentersunit2, DropLinkunit, Coordsunit2 ! file params

		INTEGER    :: i, ii, j, k, i2, filestat, sizeOld, sizeTot, &
			oldneedNum, m, mm, n, justSnapped

		INTEGER, DIMENSION(:), ALLOCATABLE :: OldblobNum, NeedleCloseBlob, UsedBlobNum
		INTEGER, DIMENSION(:,:),ALLOCATABLE :: CountShareBlob ! how many; their new blob numbers
		REAL*8    :: dist, distTemp, chargeTemp, distTempCharge, sumArea
		REAL*8, DIMENSION(:), ALLOCATABLE :: OldPCharge, OldSurfCharge, temprc, tempzc
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: OldDropCen

		CONTINUE
		BlobCentersunit2=34; DropLinkunit=25; Coordsunit2=26
		oldneedNum=-1; distTemp=-1.

		ALLOCATE(UsedBlobNum(blobNum))
		ALLOCATE(OldDropCen(InitArySize,2))
		ALLOCATE(OldPCharge(SmallInitAry))
		ALLOCATE(OldblobNum(SmallInitAry))
		ALLOCATE(NeedleCloseBlob(SmallInitAry))
		ALLOCATE(temprc(blobNum)); ALLOCATE(tempzc(blobNum))
		DO i=1,blobNum
			UsedBlobNum(i)=0
			temprc(i)=rc(i); tempzc(i)=zc(i)
		ENDDO
		DO i=1,InitArySize
			OldDropCen(i,1)=-1.; OldDropCen(i,2)=-1.
		ENDDO
		DO i=1,SmallInitAry
			OldPCharge(i)=-1.; OldblobNum(i)=-1; NeedleCloseBlob(i)=-1
		ENDDO

		j=0
		OPEN(UNIT=BlobCentersunit2, FILE='BlobCenters.dat', STATUS='old', &
			ACTION='read', IOSTAT=filestat)
  		! be careful of an existing but size zero file... Load past blob centers
		ReadBlobIf: IF (0==filestat) THEN ! the file exists
			DO WHILE (.TRUE.)     ! read until EOF
				j=j+1

    			! drop (r,z) new position, blob num, pancharge, panarea
				READ(unit=BlobCentersunit2, fmt=36, end=2424) OldDropCen(j,1), &
					OldDropCen(j,2), OldblobNum(j), OldPCharge(j), distTemp
			ENDDO
			2424 j=j-2
		ENDIF ReadBlobIf
		CLOSE(UNIT=BlobCentersunit2, STATUS='keep')

  		! load past surface shape
		IF (ShowEachSurfCharge/=1) THEN
			OPEN(UNIT=Coordsunit2, FILE='SurfaceCharge.dat', STATUS='old', ACTION='read')
		ELSE
			IF (counter <=1) THEN
				CALL Int2Str(cSurfCharge2, 'SurfChr0000', counter,'.dat',4)
			ELSE
				CALL Int2Str(cSurfCharge2, 'SurfChr0000', counter-2,'.dat',4)
			ENDIF

			OPEN(unit=Coordsunit2, file=cSurfCharge2, STATUS='old', ACTION='read', IOSTAT=filestat)
			IF (filestat/=0) THEN
				CLOSE(UNIT=Coordsunit2)
				OPEN(UNIT=Coordsunit2, FILE='SurfaceCharge.dat', STATUS='old', ACTION='read')
			ENDIF
		ENDIF

		DO i=1,7       ! skips header info
			READ(unit=Coordsunit2, fmt='(a)', END=2449) chara
		ENDDO
		2449 READ(unit=Coordsunit2, fmt='(2i5)', END=2451) sizeOld, oldneedNum
		2451 sizeTot=sizeOld+j
		ALLOCATE(OldSurfCharge(sizeTot)); 
		DO i=1,sizeTot
			OldSurfCharge(i)=-1.
		ENDDO

		DO i=j+1,sizeTot
			READ(unit=Coordsunit2, fmt='(2e15.6, 2e13.4)', END=2459) OldDropCen(i,1), OldDropCen(i,2), &
				OldSurfCharge(i), distTemp
		ENDDO
		2459 CLOSE(UNIT=Coordsunit2, STATUS='keep')

		justSnapped=0
		ALLOCATE(CountShareBlob(j, blobNum)); 
		DO i=1,j
			DO mm=1,blobNum
				CountShareBlob(i,mm)=0
			ENDDO
		ENDDO
		AvgCenMove=0.; m=0
		BlobCenDo: DO i=1,blobNum
			SkipNeed: IF ( i/=NewNeedNum ) THEN
				dist=1e6; i2=-1
    			!PRINT *, '    i, rc(i), zc(i)=',i, rc(i), zc(i)
				DO k=1,sizeTot
					distTemp = ( (rc(i)-OldDropCen(k,1))*(rc(i)-OldDropCen(k,1)) + &
						(zc(i)-OldDropCen(k,2))*(zc(i)-OldDropCen(k,2)) )

					IF ( distTemp<=dist ) THEN
						justSnapped=0
						dist=distTemp

						IF (k<=j) THEN
							i2=OldblobNum(k)
						ELSE   ! closest to old needle surface
							justSnapped=1
							IF (NewNeedNum==oldneedNum) THEN
								i2=i ! new blob is away from axis
							ELSE  ! new blob is on axis
								i2=NewNeedNum
							ENDIF
						ENDIF
					ENDIF
				ENDDO
				dist=SQRT(dist)
				DropLinks(i)=i2; centDist(i)=dist
				AvgCenMove=AvgCenMove+dist

    			! marking drops that just snapped off needle
				IF ( 1==justSnapped ) THEN
					m=m+1
					NeedleCloseBlob(m)=i
					newdropNum=i
				ELSEIF (j>0) THEN  ! linking non-needle blob counts to old blobs
					n=1
					DO WHILE (.TRUE.)
						n=n+1
						IF (CountShareBlob(i2,1) == 0) THEN
							CountShareBlob(i2,1)=1
							CountShareBlob(i2,2)=i
							EXIT
						ELSEIF (CountShareBlob(i2,n) == 0) THEN
							CountShareBlob(i2,1)=CountShareBlob(i2,1)+1
							CountShareBlob(i2,n)=i
							EXIT
						ENDIF
					ENDDO
				ENDIF
			ELSE
				DropLinks(i)=oldneedNum
    			!PRINT *, '    TrackDrop: oldneedNum,NewNeedNum=',oldneedNum,NewNeedNum
				ChargeDrop( oldneedNum ) = OldPCharge( oldneedNum )  
				rc(oldneedNum)=temprc(i); zc(oldneedNum)=tempzc(i)

				UsedBlobNum(oldneedNum)=1 ! have now used this blob number
			ENDIF SkipNeed
		ENDDO BlobCenDo
		AvgCenMove=1.*AvgCenMove/j

		DO i=1,j   ! going over all OLD blobs to update charges on NEW drops
			k=CountShareBlob(i,1)
			IF ( k>1 ) THEN ! a drop has broken into two, NOT from the needle
				sumArea=0.
				DO mm=2,k+1
					sumArea=sumArea+Areab(CountShareBlob(i,mm))
				ENDDO
				DO mm=2,k+1 ! (sec cur area/total cur area) *(old blob charge)
					n=CountShareBlob(i,mm)
					CALL ChooseNextDropNum(DropLinks(n), UsedBlobNum, blobNum)
					ii=DropLinks(n)
					ChargeDrop(ii) = Areab(ii) * OldPCharge(i) /sumArea
					rc(ii)=temprc(n); zc(ii)=tempzc(n)
     				!PRINT *, i, mm, n, ii, Areab(ii), OldPCharge(i), ChargeDrop(ii)
				ENDDO
			ELSEIF (1==k) THEN ! a drop has drifted along
				n=CountShareBlob(i,2)
				CALL ChooseNextDropNum(DropLinks(n), UsedBlobNum, blobNum)
				ii=DropLinks(n)
				ChargeDrop(ii) = OldPCharge(i)
				rc(ii)=temprc(n); zc(ii)=tempzc(n)
    			!PRINT *, '    i, ii, DropLinks(ii), ChargeDrop(ii)=',i, ii, DropLinks(ii), ChargeDrop(ii)
   			! ELSE, the needle to needle connection
			ENDIF
		ENDDO
		DEALLOCATE(CountShareBlob)

  		! -- just for new drops that have snapped off tip in last timestep
  		!   as a _very_ rough starting point, for the first timestep after the droplet
  		!   detaches, give it a total charge as the sum of panel charges as are in a 
  		!   circle 1.0x [arbitrary] distance to needle surface
		IF (m>0) THEN
			DO i=1,m
				n = NeedleCloseBlob(i)
				distTempCharge=1.2*centDist(n)

				chargeTemp=0
				DO k=j+1,sizeTot ! taking from last timestep's 'SurfaceCharge.dat' or numbered
					distTemp = SQRT( (temprc(n)-OldDropCen(k,1))*(temprc(n)-OldDropCen(k,1)) + &
						(tempzc(n)-OldDropCen(k,2))*(tempzc(n)-OldDropCen(k,2)) )
					IF (distTemp<distTempCharge) THEN
						chargeTemp=chargeTemp+OldSurfCharge(k)
					ENDIF
				ENDDO
				CALL ChooseNextDropNum(DropLinks(n), UsedBlobNum, blobNum)
				ii=DropLinks(n)
				ChargeDrop(ii)=chargeTemp
				rc(ii)=temprc(n); zc(ii)=tempzc(n)
    			!PRINT *, '   n, DropLinks(n)=',n, ii
			ENDDO
		ENDIF
		DEALLOCATE(NeedleCloseBlob); DEALLOCATE(OldSurfCharge)
		DEALLOCATE(UsedBlobNum)
		DEALLOCATE(temprc); DEALLOCATE(tempzc)

		IF ( 1==ShowDropLink ) THEN !.AND. 0==MOD(counter-1, SkipSteps)
   			! writes new blob numbers to file
			CALL Int2Str(DLink2, 'DropLn', counter-1,'.dat',4)
			PRINT *, '  -Writing ', DLink2
			OPEN(UNIT=DropLinkunit, FILE=DLink2, STATUS='replace')

			DO i=1,blobNum
				WRITE(UNIT=DropLinkunit, FMT='(2i8, e12.4)') DropLinks(i), i, centDist(i)
			ENDDO
			CLOSE(UNIT=DropLinkunit, STATUS='keep')
		ENDIF

		GOTO 993

		DLink2=chara ! doesn't do anything, but suppresses warning
  		!-----------------------------
  		! FORMAT listings
		36 FORMAT(2e14.6, i6, 2e14.6)

		993 DEALLOCATE(OldDropCen); DEALLOCATE(OldPCharge); DEALLOCATE(OldblobNum)

	CONTAINS 
		SUBROUTINE ChooseNextDropNum(NextNum, UsedBlobNmbr, totblobs)
   			! chooses the next available droplet number to assign the coord blobs to
			INTEGER, INTENT(in)  :: totblobs

			INTEGER, INTENT(inout) :: NextNum
			INTEGER, DIMENSION(totblobs), INTENT(inout) :: UsedBlobNmbr

			INTEGER :: ii

			CONTINUE 

			DO ii=1,totblobs
				IF (0==UsedBlobNmbr(ii)) THEN
					NextNum=ii
					UsedBlobNmbr(ii)=1
					EXIT
				ENDIF
				IF (totblobs==ii) THEN
					PRINT *, 'Err MAT_Greens:ChooseNextDropNum'
					STOP
				ENDIF
			ENDDO
		END SUBROUTINE ChooseNextDropNum
	END SUBROUTINE TrackDrop

! ***************
! SUBROUTINE ChangePanCharge
! ***************
! DESCRIPTION: recalculates the total charge and forces on each panel to keep electron
!   *density* (via area) constant between timesteps
!
! INPUT: New panel numbering matrix, array sizes, normal electric fields 
!
! OUTPUT: new panel charges and forces 

	SUBROUTINE ChangePanCharge(blobforce, newcharge, ptspan, efield, swapblobNum, &
			panareas, panChargeDens, dropmovedist, startpanNum, distToCen, totpts, &
			blobNum2, mxblobpts)

  			! incoming
		INTEGER, INTENT(in)       :: totpts, blobNum2, mxblobpts
		INTEGER, DIMENSION(blobNum2), INTENT(IN) :: ptspan, swapblobNum, startpanNum
		REAL*8, INTENT(in)       :: dropmovedist
		REAL*8, DIMENSION(blobNum2), INTENT(IN)  :: panChargeDens, panareas, distToCen

  			! outgoing
		REAL*8, DIMENSION(blobNum2), INTENT(inout)  :: newcharge
		REAL*8, DIMENSION(totpts), INTENT(inout) :: efield
		REAL*8, DIMENSION(blobNum2,mxblobpts), INTENT(out) :: blobforce

  			! internal
		INTEGER    :: i, j, k, pan
		REAL*8, DIMENSION(blobNum2) :: tempcharge, ChargeRatio
		REAL*8, DIMENSION(totpts) :: tempefield

		CONTINUE

  			! get new charges, assign efield values & recalculate forces
		tempcharge=newcharge; tempefield=efield

		DO i=1,blobNum2
			IF (distToCen(i)<2.*dropmovedist .AND. panChargeDens(i)>1e-10) THEN ! linked to an older dropler
				newcharge(i) = panareas(i)*panChargeDens(i)
			ELSE           ! a brand new drop. retain full charge
				newcharge(i) = tempcharge(i)
			ENDIF

			IF (ABS(tempcharge(i))>1e-10) THEN
				ChargeRatio(i) = newcharge(i)/tempcharge(i) ! (should be/is) droplet charge
			ELSE
				ChargeRatio(i)=0
			ENDIF

			k=ptspan(swapblobNum(i))
			pan=startpanNum(i)-1
			DO j=1,k
				pan=pan+1

				IF (ABS(ChargeRatio(i))>1e-10) THEN
					efield(pan)=tempefield(pan)*ChargeRatio(i)
				ELSE
					efield(pan)=tempefield(pan)
				ENDIF
				blobforce(i,j)= -epsil0*efield(pan)*efield(pan)
			ENDDO
		ENDDO
	END SUBROUTINE ChangePanCharge

	SUBROUTINE qsort(a, t, n)
  		!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
  		!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
  		!     ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
		INTEGER, INTENT(IN)    :: n
		REAL*8, DIMENSION(N), INTENT(INOUT)  :: a, t

  		!     Local Variables
		INTEGER                :: i, j, k, l, r, s, stackl(15), stackr(15)
		REAL*8                 :: w, ww, x

		s = 1
		stackl(1) = 1
		stackr(1) = n

  		!		KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.
		10 CONTINUE
		l = stackl(s)
		r = stackr(s)
		s = s - 1

  !     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.
		20 CONTINUE
		i = l
		j = r
		k = (l+r) / 2
		x = a(k)

  !     REPEAT UNTIL I > J.
		DO
			DO
				IF (a(i).LT.x) THEN  ! Search from lower end
					i = i + 1
					CYCLE
				ELSE
					EXIT
				END IF
			END DO

			DO
				IF (x.LT.a(j)) THEN  ! Search from upper end
					j = j - 1
					CYCLE
				ELSE
					EXIT
				END IF
			END DO

			IF (i.LE.j) THEN   ! Swap positions i & j
				w = a(i)
				ww = t(i)
				a(i) = a(j)
				t(i) = t(j)
				a(j) = w
				t(j) = ww
				i = i + 1
				j = j - 1
				IF (i.GT.j) EXIT
			ELSE
				EXIT
			END IF
		END DO

		IF (j-l.GE.r-i) THEN
			IF (l.LT.j) THEN
				s = s + 1
				stackl(s) = l
				stackr(s) = j
			END IF
			l = i
		ELSE
			IF (i.LT.r) THEN
				s = s + 1
				stackl(s) = i
				stackr(s) = r
			END IF
			r = j
		END IF

		IF (l.LT.r) GO TO 20
		IF (s.NE.0) GO TO 10

		RETURN
	END SUBROUTINE qsort

		!-----------------
		! checking bounds of most global variables
	SUBROUTINE CheckErrs(multyDist)
		REAL*8, INTENT(INOUT)  :: multyDist

		CONTINUE

		IF (SkipSteps<1) THEN
			PRINT *, 'Need to have a larger positive SkipSteps value in global.'
			STOP
		ENDIF

		IF (1==RecastDist) THEN
			PRINT *, '  *Expecting to get "phi" data in [cm]'
			multyDist = 100.
		ELSE
			PRINT *, '  *Expecting to get "phi" data in [m]'
			multyDist = 1.
		ENDIF

		IF (CaseShape<1 .OR. CaseShape>4) THEN
			PRINT *, 'global:CaseShape out of range'
			STOP
		ENDIF

  		! checking all global 1/0 variables
		CALL NotOneZero(AformDouble,  'AformDouble    ')
		CALL NotOneZero(CalcPotentials, 'CalcPotentials ')
		CALL NotOneZero(CalcPotBack,  'CalcPotBack    ') 
		CALL NotOneZero(ConnectDrops,  'ConnectDrops   ')
		CALL NotOneZero(CubicSplineFit, 'CubicSplineFit ')
		CALL NotOneZero(ForceEdgeStable,'ForceEdgeStable')
		CALL NotOneZero(LoadPastStep,  'LoadPastStep   ')
		CALL NotOneZero(RecastDist,  'RecastDist     ')
		CALL NotOneZero(ShowGrid,   'ShowGrid       ')
		CALL NotOneZero(ShowNewCross,  'ShowNewCross   ')
		CALL NotOneZero(ShowSort,   'ShowSort       ')
		CALL NotOneZero(ShowPt_Force,  'ShowPt_Force   ')
		CALL NotOneZero(ShowNewCoords,  'ShowNewCoords  ')
		CALL NotOneZero(ShowCoords,  'ShowCoords     ')
		CALL NotOneZero(ShowPassedPhi,  'ShowPassedPhi  ')
		CALL NotOneZero(ShowAngleSort,  'ShowAngleSort  ')
		CALL NotOneZero(ShowSubPan,  'ShowSubPan     ')
		CALL NotOneZero(ShowDropLink,  'ShowDropLink   ')
		CALL NotOneZero(ShowChiAngle,  'ShowChiAngle   ')
		CALL NotOneZero(ShowAMat,   'ShowAMat       ')
		CALL NotOneZero(ShowGMResConverge,  'ShowGMResConver')
		CALL NotOneZero(ShowEachSurfCharge, 'ShowEachSurfCha')
		CALL NotOneZero(SkipFar,   'SkipFar        ')
		CALL NotOneZero(UseGaussQuad,  'UseGaussQuad   ')
		CALL NotOneZero(UsePreCond,  'UsePreCond     ')
	END SUBROUTINE CheckErrs

! ***************
! SUBROUTINE CalcChiAngle
! ***************
! DESCRIPTION: optional routine that calculates & prints surface angles chi
!
! INPUT: PForce data (r,z)
!
! OUTPUT: angle, radial time chi=eta(r,t)
	SUBROUTINE CalcChiAngle(rpos, zpos, numPts)
  			! incoming
		INTEGER, INTENT(in)       :: numPts
		REAL*8, DIMENSION(numPts), INTENT(IN)  :: rpos, zpos

  			! outgoing -- nothing
  			! internal
		INTEGER, PARAMETER :: Angunit=10
		CHARACTER(LEN=15) :: angChar
		INTEGER    :: i
		REAL*8    :: surfang, diffr, changeunit

		CONTINUE

		IF (1==RecastDist) THEN
			changeunit=100. 
		ELSE 
			changeunit=1.
		ENDIF

		CALL Int2Str(angChar, 'ChiAng', counter-1,'.dat',4)
		PRINT *, '  -Writing ', angChar
		OPEN(unit=Angunit, file=angChar, status='replace')
		WRITE(unit=Angunit, fmt='(a)') 'VARIABLES="R [cm]"'
		WRITE(unit=Angunit, fmt='(a)') '"<greek>c</greek> [<sup>o</sup>]"'
		WRITE(unit=Angunit, fmt='(a)') '"Panel [#]"'

		DO i=1, numPts-2
			diffr=rpos(i+1)-rpos(i)
			IF ( ABS(diffr)-eps<0 ) THEN
				surfang=180.
			ELSE
				surfang=180.*ATAN((zpos(i+1)-zpos(i))/diffr)/pi+90.
			ENDIF
			WRITE(unit=Angunit, fmt='(2e12.4, i8)'), rpos(i)*changeunit, surfang, i
		ENDDO
		CLOSE(UNIT=Angunit, STATUS='keep')

	END SUBROUTINE CalcChiAngle

! ***************
! SUBROUTINE ForceRange
! ***************
! DESCRIPTION: optional routine that forces surface points near the right edge 
! to match edge and fit smooth function
!
! INPUT: Coordinates
!
! OUTPUT: smoothed Coordinates
	SUBROUTINE ForceRange(Coords, ZeroR, NeedleNum, sizeLevel, BlobNum)
  			! incoming
		INTEGER, INTENT(in)  :: BlobNum, NeedleNum
		INTEGER, DIMENSION(BlobNum), INTENT(in)   :: sizeLevel

  			! outgoing  
		REAL*8, INTENT(inout) :: ZeroR
		REAL*8,DIMENSION(InitArySize,8),INTENT(inout) :: Coords

  			! internal
		INTEGER     :: i, j, k, startCut, tempRange, &
			Needlestart, OldNumCut
		REAL*8     :: startHeight, maxR, minR

		CONTINUE

		tempRange=1
		DO i=1,BlobNum  ! finds the starting and ending point of needle section
			IF (i>1) THEN
				tempRange = tempRange+ sizeLevel(i-1)
			ENDIF
			IF (NeedleNum==i) THEN
				Needlestart=tempRange
			ENDIF
		ENDDO
		k = sizeLevel(NeedleNum)

		maxR=-1e6; minR=1e6
		DO i=1,k
			IF (Coords(tempRange+i-1,1)<minR) minR=Coords(tempRange+i-1,1)
			IF (Coords(tempRange+i-1,1)>maxR) maxR=Coords(tempRange+i-1,1)
		ENDDO
		ZeroR=0.9*(maxR-minR)+minR

		IF (0==ZeroR) THEN
			PRINT *, 'Err MAT_Greens:ForceRange ZeroR=0'
			STOP
		ENDIF

		i=1; j=0
		DO WHILE(.TRUE.)
			IF (Coords(i,1)<ZeroR) THEN
				j=j+1
			ELSE
				EXIT
			ENDIF
			i=i+1
		ENDDO
		OldNumCut=sizeLevel(NeedleNum)-j
		startCut=Needlestart+j
		k=startCut+OldNumCut-1

		startHeight=Coords(j,2)/EXP(-22.9*Coords(j,1)*Coords(j,1))
		DO i = startCut, k
			Coords(i,2) = startHeight*EXP(-22.9*Coords(i,1)*Coords(i,1))
		ENDDO

		PRINT *, '    Forcing ', OldNumCut, ' points to baseline values.'

	END SUBROUTINE ForceRange

! ***************
! SUBROUTINE NotOneZero
! ***************
! DESCRIPTION: base check to see if a variable is either a one (1) or zero (0)

	SUBROUTINE NotOneZero(checkvarI, namevar)
  			! incoming
		INTEGER, INTENT(in)  :: checkvarI

		CHARACTER(len=15), INTENT(in)  :: namevar
		CONTINUE

		IF ( .NOT.(0==checkvarI .OR. 1==checkvarI) ) THEN
			PRINT *, 'global variable', namevar, ' incorrect!'
			STOP
		ENDIF
	END SUBROUTINE NotOneZero

! ***************
! SUBROUTINE FlipPanPtOrder
! ***************
! DESCRIPTION: checks to see if the order of the panel points should be flipped. 
!   Uses farthest pt from blob center to see if the ordered list goes clockwise
!   or counterclockwise for an arbitrary shaped, non-crossing blob
!
! INPUT: panel points, blob centers, far point
!
! OUTPUT: updated panel points
	SUBROUTINE FlipPanPtOrder(PanPts, i2, PtC, maxRout)
  		! incoming variables
		INTEGER, INTENT(in)		:: i2
		REAL*8, INTENT(in)      :: maxRout
		REAL*8, DIMENSION(2), INTENT(in)   :: PtC

  		! outgoing function variable
		TYPE(PANEL), DIMENSION(i2),INTENT(inout):: PanPts

  		! internal variables 
		INTEGER		:: j, z1, iPlus, iMinus, reverseyes, connectDown
		REAL*8    	:: delR, delZ, dot1, dot2, dist, &
			ref, sgndot1, sgndot2, minZ, maxR
		REAL*8, DIMENSION(2):: s1norm, s2norm, scnorm, FarPt, PtA
		TYPE(PANEL), DIMENSION(i2) :: orderedPan

		CONTINUE
		dist=0.; minZ=1e10; maxR=0; reverseyes=5

		DO j=1,i2  ! makes local temporary copy
			CALL ZeroorCopyPanel( orderedPan(j), 0, PanPts(j) )
		ENDDO

  		! form sorted list and determine furthest point
		DO j=1,i2
			PtA(1)=PanPts(j)%r0; PtA(2)=PanPts(j)%z0
			IF (PtA(1)>maxR) maxR=PtA(1)
			IF (PtA(2)<minZ) minZ=PtA(2)

			ref = SQRT( (PtA(1)-PtC(1))*(PtA(1)-PtC(1)) + (PtA(2)-PtC(2))*(PtA(2)-PtC(2)) )
			IF (ref>dist) THEN
				FarPt(1)=PtA(1); FarPt(2)=PtA(2)
				dist=ref; z1=j
			ENDIF
		ENDDO

  		! look at point before and after
		iPlus=Nexti( z1, i2, 0 )
		connectDown=0
		IF (1==iPlus .AND. ABS(PanPts(i2)%r0/maxRout)>0.9) connectDown=1 
			! you are at the far edge
		IF (1==z1) THEN 
			iMinus=i2
		ELSE 
			iMinus=z1-1
		ENDIF

  		!------- center->farthest pt norm define
		delR = FarPt(1) - PtC(1)
		delZ = FarPt(2) - PtC(2)
		dist = SQRT(delR*delR + delZ*delZ)
		scnorm(1) = delR/dist
		scnorm(2) = delZ/dist

 		!------- left(usually) of farthest pt
		delR = PanPts(iMinus)%r0 - FarPt(1)
		delZ = PanPts(iMinus)%z0 - FarPt(2)

		dist = SQRT(delR*delR + delZ*delZ)
		s1norm(1) = -delZ/dist
		s1norm(2) = delR/dist

		dot1=DOT_PRODUCT(scnorm,s1norm)
		sgndot1=dot1/ABS(dot1)

  		!------- right (usually) of farthest pt
		IF (1==connectDown ) THEN
			delR = maxR - FarPt(1)
			delZ = minZ - FarPt(2)
		ELSE
			delR = PanPts(iPlus)%r0 - FarPt(1)
			delZ = PanPts(iPlus)%z0 - FarPt(2)
		ENDIF
		dist = SQRT(delR*delR + delZ*delZ)

		IF ( dist<eps ) THEN
			s2norm(1) = 0
			s2norm(2) = 0
			dot2=DOT_PRODUCT(scnorm,s2norm)
			sgndot2=0
		ELSE
			s2norm(1) = -delZ/dist
			s2norm(2) = delR/dist
			dot2=DOT_PRODUCT(scnorm,s2norm)
			sgndot2=dot2/ABS(dot2)
		ENDIF

		IF (sgndot1/=sgndot2 .AND. (dot1/=0 .AND. dot2/=0)) THEN !sign(0)=+
			IF (dot1<0) reverseyes=1 ! going CCW
			IF (dot2<0) reverseyes=0 ! going CW
		ELSEIF (dot1/=0 .AND. (0==dot2 .OR. sgndot1==sgndot2) ) THEN 
			reverseyes=INT( 0.5-sgndot1/2. ) ! 'L' opening left
		ELSEIF (dot2/=0 .AND. 0==dot1 ) THEN 
			reverseyes=INT( 0.5+sgndot2/2. ) ! 'L' opening right
		ELSE
			PRINT *, '  Err MAT_Greens:FlipPanPtOrder sgn1, sgn2, dot1,dot2=', &
				sgndot1, sgndot2, dot1,dot2
			reverseyes=1    ! degenerate baseline case
		ENDIF

		IF (1==reverseyes) THEN
			DO j=1,i2
    			!....
				CALL ZeroorCopyPanel(PanPts(j), 0, orderedPan(j), &
					reverseyes)
    			!....
			ENDDO
		ELSEIF(5==reverseyes) THEN
			PRINT *, 'Err, FlipPanPtOrder:reverseyes'
			STOP
		ENDIF
	END SUBROUTINE FlipPanPtOrder

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
! CALLING PROGRAM: MAT_Greens

	SUBROUTINE SortNN_NoLong( ArrayNew, sizeBubble )  
  		! incoming variables
		INTEGER, INTENT(in)        :: sizeBubble

  		! outgoing variables
		TYPE(PANEL), DIMENSION(sizeBubble), INTENT(inout) :: ArrayNew

  		! subroutine entirely variables
		LOGICAL        :: ReRun
		INTEGER, DIMENSION(sizeBubble)  :: indexIn, forceLink
		INTEGER        :: nearPtNow, nearPtNext, nearPtLast, &
			i, j, k, tempLeft, countRun, countRepeat
		REAL*8, PARAMETER    :: MultDist=4.0
		REAL*8       :: outrangeVar, dist, distNeighbor, refLength, &
			delX, delY, PanLen2, deltaNorm, maxNorm, minNorm, xnormNow, ynormNow, &
			xnormNext, ynormNext, xnormTemp, ynormTemp
		REAL*8, DIMENSION(sizeBubble) :: rOld, zOld
		TYPE(PANEL), DIMENSION(sizeBubble)  :: ArrayOld

		CONTINUE
 	 	!*************
  		! main program commands begin
		DO i=1,sizeBubble
   				!....
			CALL ZeroorCopyPanel(ArrayOld(i), 0, ArrayNew(i))
   				!....
			rOld(i)=ArrayOld(i)%r0
			zOld(i)=ArrayOld(i)%z0
			forceLink(i)=-1
		ENDDO

  		! left-most point gotten as furthest left x  
		outrangeVar=1e6; tempLeft=1
		DO i=1,sizeBubble
			IF (rOld(i)< rOld(tempLeft) ) THEN
				tempLeft=i
			ENDIF
		ENDDO
		nearPtNow=tempLeft; nearPtLast=nearPtNow

		ReRun=.TRUE.; countRun=0; countRepeat=0
		ReRunDo: DO WHILE (ReRun)
			102 ReRun=.FALSE.; countRun=countRun+1 ! hate doing this, but need to exit out 3 loops at once
			refLength=ABS(rOld(tempLeft)-rOld(tempLeft+1))

			DO i=1,sizeBubble
				indexIn(i)=1 
			ENDDO
			indexIn(tempLeft)=0; nearPtNow=tempLeft
   			!....
			CALL ZeroorCopyPanel(ArrayNew(1), 0, ArrayOld(tempLeft))
   			!....
   			! loops through and organizes 
			outer: DO i=1,sizeBubble-1
				forcelinkIF: IF ( forceLink(nearPtNow) /= -1 ) THEN ! you are on an iterated run and need to forcelink
     				! to a particular point to avoid overall error
					nearPtNext= forceLink(nearPtNow)

     				! calculate norm to pt
					delX = rOld(nearPtNext) - rOld(nearPtNow)
					delY = zOld(nearPtNext) - zOld(nearPtNow)

					PanLen2  = SQRT( delX*delX + delY*delY )
					xnormNext= -delY/PanLen2
					ynormNext= delX/PanLen2
				ELSE 
					dist=outrangeVar; nearPtNext=0; maxNorm=0; minNorm=1000
					NNNormDo: DO j=1,sizeBubble
						indexIf: IF ( 1==indexIn(j) ) THEN
							distNeighbor = (rOld(nearPtNow)-rOld(j))*(rOld(nearPtNow)-rOld(j)) + &
								(zOld(nearPtNow)-zOld(j))*(zOld(nearPtNow)-zOld(j))

       						! calculate norm to pt
							delX = rOld(j) - rOld(nearPtNow)
							delY = zOld(j) - zOld(nearPtNow)

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
						refLength=((i-1)*refLength+dist)/i ! compute a running average on the reference length
					ELSE
      					! do nothing - refLength is tentatively correct as-is
					ENDIF
					NN_noNormIF: IF ( 0==nearPtNext ) THEN ! if no connection, just do nearby NN
						DO j=1,sizeBubble
							IF ( 1==indexIn(j) ) then
								distNeighbor = (rOld(nearPtNow)-rOld(j))*(rOld(nearPtNow)-rOld(j)) + &
									(zOld(nearPtNow)-zOld(j))*(zOld(nearPtNow)-zOld(j))

								IF ( distNeighbor < dist .AND. distNeighbor<MultDist*refLength) THEN
									dist = distNeighbor
									nearPtNext=j

         							! calculate norm to pt
									delX = rOld(j) - rOld(nearPtNow)
									delY = zOld(j) - zOld(nearPtNow)

									PanLen2  = SQRT( delX*delX + delY*delY )
									xnormNext= -delY/PanLen2
									ynormNext= delX/PanLen2
								ENDIF
							ENDIF
						ENDDO
					ENDIF NN_noNormIF

					NN_noDistIF: IF ( 0==nearPtNext ) THEN 
      					! if still no connection, *major* miss. Next point is wrong, 
						!   because too far away; 1) find it. 2) it's NN already _taken_. 
						!   3) rerun
						dist=outrangeVar; countRepeat=countRepeat+1
						DO j=1,sizeBubble ! doing (1)
							IF ( 1==indexIn(j) ) THEN
								distNeighbor = (rOld(nearPtNow)-rOld(j))*(rOld(nearPtNow)-rOld(j)) + &
									(zOld(nearPtNow)-zOld(j))*(zOld(nearPtNow)-zOld(j))

								IF ( distNeighbor < dist ) THEN
									dist = distNeighbor
									nearPtNext=j
								ENDIF
							ENDIF
						ENDDO

						dist=outrangeVar
						DO j=1,sizeBubble ! doing (2)
							IF ( indexIn(j)==0 .AND. j/=nearPtNext ) THEN
								distNeighbor = (rOld(nearPtNext)-rOld(j))*(rOld(nearPtNext)-rOld(j)) + &
									(zOld(nearPtNext)-zOld(j))*(zOld(nearPtNext)-zOld(j))

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

						ReRun=.TRUE.  ! doing (3)
						forceLink(k)=nearPtNext

						IF (countRepeat>5) THEN ! you are stuck in an interative loop. Find
       							! nearest non-prior point and force connect
							dist=outrangeVar
							DO j=1,sizeBubble 
								IF ( j /= nearPtLast) THEN
									distNeighbor = (rOld(nearPtNow)-rOld(j))*(rOld(nearPtNow)-rOld(j)) + &
										(zOld(nearPtNow)-zOld(j))*(zOld(nearPtNow)-zOld(j))

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
    				!....
				CALL ZeroorCopyPanel(ArrayNew(i+1), 0, ArrayOld(nearPtNext))
    				!....
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
! CALLING PROGRAM: MAT_Greens

	SUBROUTINE SortWeighted_AngDist( ArrayNew, sizeBubble )  
  		! incoming variables
		INTEGER, INTENT(in)        :: sizeBubble

  		! outgoing variables
		TYPE(PANEL), DIMENSION(sizeBubble), INTENT(inout) :: ArrayNew

 		! subroutine entirely variables
		INTEGER       :: i, j, tempLeft, nearPtNow, nearPtNext
		INTEGER, DIMENSION(sizeBubble) :: indexIn
		REAL*8       :: WeightPt, LowestWeightPt, outrangeVar, xnormNow, &
			ynormNow, distNeighbor, delX, delY, PanLen2, xnormTemp, ynormTemp, xnormNext, &
			ynormNext, angNext, angNow
		REAL*8, PARAMETER    :: MultDist=1.0
		REAL*8, DIMENSION(sizeBubble) :: xOld, yOld
		REAL*8, DIMENSION(sizeBubble,8) :: LocalNNCand, LocalAngleCand
   			! 4 NN points, with ,1 being the closest and ,4 being the furthest
   			!   and 5=j of ,1 with 8=j of ,4. 'Angle(,1)'=weight dist. 'Angle(,5)'=weight angle norm
		TYPE(PANEL),DIMENSION(sizeBubble):: ArrayOld
    		! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, type panel]

		CONTINUE
  		! main program commands begin
		DO i=1,sizeBubble
   			!....
			CALL ZeroorCopyPanel(ArrayOld(i), 0, ArrayNew(i))
   			!....
			xOld(i)=ArrayOld(i)%r0
			yOld(i)=ArrayOld(i)%z0
			indexIn(i)=1
		ENDDO

  		! left-most point gotten as furthest left x  
		outrangeVar=1e6; tempLeft=1
		DO i=1,sizeBubble
			IF (xOld(i)< xOld(tempLeft) ) THEN
				tempLeft=i
			ENDIF
		ENDDO
		indexIn(tempLeft)=0; nearPtNow=tempLeft
  		!....
		CALL ZeroorCopyPanel(ArrayNew(1), 0, ArrayOld(tempLeft))
  		!....
		xnormNow=-1.0; ynormNow=0.0 
    	! b/c starts directly to right, assume a vertical panel going down

		DO i=1,sizeBubble
			LocalNNCand(i,1)=outrangeVar; LocalNNCand(i,5)=outrangeVar 
			LocalNNCand(i,2)=outrangeVar; LocalNNCand(i,6)=outrangeVar
			LocalNNCand(i,3)=outrangeVar; LocalNNCand(i,7)=outrangeVar
			LocalNNCand(i,4)=outrangeVar; LocalNNCand(i,8)=outrangeVar
			LocalAngleCand(i,1)=-1.; LocalAngleCand(i,5)=-1.
			LocalAngleCand(i,2)=-1.; LocalAngleCand(i,6)=-1.
			LocalAngleCand(i,3)=-1.; LocalAngleCand(i,7)=-1.
			LocalAngleCand(i,4)=-1.; LocalAngleCand(i,8)=-1.
		ENDDO
  		! find 4 nearest neighbors at each point
		distangDo: DO i=1,sizeBubble-1
			DO j=1,sizeBubble
				jIdx: IF ( 1==indexIn(j) ) THEN
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
				ENDIF jIdx
			ENDDO

   			! note: LocalNNCand(i,1) is the shortest distance
			IF ( outrangeVar==LocalNNCand(nearPtNow,5) .OR. outrangeVar==LocalNNCand(nearPtNow,6) &
				.OR. outrangeVar==LocalNNCand(nearPtNow,7) .OR. outrangeVar==LocalNNCand(nearPtNow,8)) THEN
				PRINT *, 'MAT_Greens:SortWeighted err. 1,2,3,4=', LocalNNCand(nearPtNow,1), &
					LocalNNCand(nearPtNow,2), LocalNNCand(nearPtNow,3), LocalNNCand(nearPtNow,4)
				PRINT *, 'MAT_Greens:SortWeighted err. 5,6,7,8=',INT(LocalNNCand(nearPtNow,5)), &
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
   			!....
			CALL ZeroorCopyPanel(ArrayNew(i+1), 0, ArrayOld(nearPtNext))
   			!....
		ENDDO distangDo
	END SUBROUTINE SortWeighted_AngDist

! ***************
! SUBROUTINE SortNN_NoCross
! ***************
! DESCRIPTION: takes theta sorted list of points and iteratively produces a new sorted list, based on 
! nearest neighbor AND not crossing over lines
!
! INPUT: Coords(x8)
!
! OUTPUT: sorted Coords(x8)
!
! CALLING PROGRAM: MAT_Greens

	SUBROUTINE SortNN_NoCross( ArrayNew, sizeBubble, WhichBlob )  
  		! incoming variables
		INTEGER, INTENT(in)    :: sizeBubble, WhichBlob

  		! outgoing variables
		TYPE(PANEL), DIMENSION(sizeBubble), INTENT(inout) :: ArrayNew

  		! subroutine entirely variables
		LOGICAL   	:: Crossed_Unknown, LoopCrossed, OkConct
		INTEGER   	:: nearPtNow, nearPtNext, i, j, m, tempLeft, countLoop
		INTEGER, DIMENSION(sizeBubble+1) 		:: NextPt, usedNum
		INTEGER, DIMENSION(sizeBubble,sizeBubble) :: indexNoConnect
		REAL*8   	:: outrangeVar, dist, distNeighbor
		REAL*8,DIMENSION(sizeBubble) 		:: xOld, yOld
		TYPE(PANEL), DIMENSION(sizeBubble) 	:: ArrayOld

		CONTINUE
  		!*************
  		! main program commands begin
		DO i=1,sizeBubble
   			!....
			CALL ZeroorCopyPanel(ArrayOld(i), 0, ArrayNew(i))
   			!....
			xOld(i)=ArrayOld(i)%r0
			yOld(i)=ArrayOld(i)%z0
		ENDDO

  		! left-most point gotten as furthest left x  
		tempLeft=1
		DO i=1,sizeBubble
			IF (xOld(i)< xOld(tempLeft) ) THEN
				tempLeft=i
			ENDIF
			DO j=1,sizeBubble-1
				indexNoConnect(i,j)=0 
			ENDDO
		ENDDO
		CALL ZeroorCopyPanel(ArrayNew(1), 0, ArrayOld(tempLeft))

		outrangeVar=1e6; Crossed_Unknown=.TRUE.
		countLoop=0
		NC_do: DO WHILE (Crossed_Unknown)
			103 nearPtNow=tempLeft; 
			countLoop=countLoop+1
			IF (0==MOD(countLoop,10000)) THEN
				PRINT *, '  ** MAT_Greens:SortNN_NoCross killing at countLoop=',countLoop
				PRINT *, '     Killed on blob, size=',WhichBlob,sizeBubble
				DO i=1,sizeBubble
   					!....
					CALL ZeroorCopyPanel(ArrayNew(i), 0, ArrayOld(i))
   					!....
				ENDDO
				EXIT
			ENDIF
			DO i=1,sizeBubble
				NextPt(i)=0; usedNum(i)=0
			ENDDO
			NextPt(1)=tempLeft; usedNum(tempLeft)=1

			LoopCrossed=.FALSE.

   			! loops through and organizes 
			outer: DO i=1,sizeBubble
				dist=outrangeVar; nearPtNext=0
				DO j=1,sizeBubble
					CALL OkConnect(OkConct, j, NextPt(i)) ! checks if ok to connect to pt(j)
					IF (countLoop>=10000) THEN
						PRINT *, 'i,j,NextPt(i),nearPtNext,OkConct=',i,j,NextPt(i),nearPtNext,OkConct
					ENDIF
					IF ( OkConct .AND. 0==usedNum(j) .AND. nearPtNow/=j) THEN
						distNeighbor = ( xOld(nearPtNow)-xOld(j))*(xOld(nearPtNow)-xOld(j)) + &
							(yOld(nearPtNow)-yOld(j))*(yOld(nearPtNow)-yOld(j) )

						IF ( distNeighbor < dist ) THEN
							dist = distNeighbor
							nearPtNext=j
						ENDIF
					ENDIF
				ENDDO
				IF ( outrangeVar == dist ) THEN ! last point, force connect to first
					nearPtNext=tempLeft
				ENDIF

				DO j=1,i-1
					CALL DoPtsCross( LoopCrossed, ArrayOld(nearPtNow)%r0,ArrayOld(nearPtNow)%z0, &
						ArrayOld(nearPtNext)%r0, ArrayOld(nearPtNext)%z0, ArrayOld(NextPt(j))%r0, &
						ArrayOld(NextPt(j))%z0, ArrayOld(NextPt(j+1))%r0, ArrayOld(NextPt(j+1))%z0 )

					IF (LoopCrossed) THEN  
      					! crossed, loop again. take (incorrect) final pt and force it to the nearest 
      					!   open neighbor 
						!PRINT *, '  Crossed link - try again', i, j, nearPtNow, nearPtNext, &
						!	NextPt(j),NextPt(j+1)
						CALL AddNoConnect(NextPt(j),NextPt(j+1))
						GOTO 103 ! hate to do this, but need to break out of only 2/3 loops
					ENDIF
				ENDDO

				NextPt(i+1)=nearPtNext; nearPtNow=nearPtNext
				usedNum(nearPtNext)=1
				IF (i<sizeBubble) THEN
    				!....
					CALL ZeroorCopyPanel(ArrayNew(i+1), 0, ArrayOld(nearPtNext))
    				!....
				ENDIF
			ENDDO outer 

			Crossed_Unknown=.FALSE. ! if you get this far, everything's good
		ENDDO NC_do
CONTAINS
  		! checks to see if there is an earlier prohibition against this point
		SUBROUTINE OkConnect(OkConct, LoopJ, LoopI) 
			INTEGER, INTENT(IN), OPTIONAL :: LoopJ, LoopI
			LOGICAL, INTENT(INOUT)   :: OkConct

			CONTINUE

			IF (LoopJ/=0 .AND. LoopI/=0) THEN
				m=1
				DO WHILE (m<sizeBubble)
					IF (LoopJ==indexNoConnect(LoopI,m)) THEN ! can't connect
						OkConct=.FALSE.
						EXIT
					ELSEIF (0==indexNoConnect(LoopI,m)) THEN  ! safe
						OkConct=.TRUE.
						EXIT
					ENDIF
					m=m+1
				END DO
			ENDIF
		END SUBROUTINE OkConnect

  		! adds a prohibition against this connection
		SUBROUTINE AddNoConnect(PtFromNoCon, PtToNoCon)
			INTEGER, INTENT(in) :: PtFromNoCon, PtToNoCon

			CONTINUE
			m=1
			DO WHILE (m<sizeBubble)
				IF (0==indexNoConnect(PtFromNoCon,m)) THEN  ! safe
					indexNoConnect(PtFromNoCon,m)=PtToNoCon
					EXIT
				ENDIF
				m=m+1
			END DO
		END SUBROUTINE AddNoConnect
	END SUBROUTINE SortNN_NoCross

! ***************
! SUBROUTINE PrintCoords
! ***************
! DESCRIPTION: 	prints out the file "Coords.dat"
	SUBROUTINE PrintCoords(Cord, NCord)
		INTEGER, INTENT(IN)	:: NCord
		REAL*8, DIMENSION(InitArySize,8), INTENT(IN)	:: Cord 

		INTEGER	:: i, Coordsunit

		CONTINUE
		Coordsunit=12

		! [r, z, pt val, pt type, blob#, blob_xc, blob_yc, open shape]
		PRINT *, '  -Writing "Coords.dat" '
		OPEN(unit=Coordsunit, file='Coords.dat', status='replace')
		WRITE(unit=Coordsunit, fmt='(a)') 'VARIABLES="R [cm]"'
		WRITE(unit=Coordsunit, fmt='(a)') '"Z [cm]"'
		WRITE(unit=Coordsunit, fmt='(a)') '"Pt val"'
		WRITE(unit=Coordsunit, fmt='(a)') '"Pt type"'
		WRITE(unit=Coordsunit, fmt='(a)') '"Blob#"'
		WRITE(unit=Coordsunit, fmt='(a)') '"Blob_rc"'
		WRITE(unit=Coordsunit, fmt='(a)') '"Blob_zc"'
		WRITE(unit=Coordsunit, fmt='(a)') '"OpenShape"'
		WRITE(unit=Coordsunit, fmt='(a)') '"Pan #"'
		DO i=1,NCord
			WRITE(unit=Coordsunit, fmt='(3e14.6, 2i4, 2e14.6, 2i5)') Cord(i,1), Cord(i,2), &
				Cord(i,3), INT(Cord(i,4)), INT(Cord(i,5)), Cord(i,6), Cord(i,7), &
				INT(Cord(i,8)), i
		ENDDO
		CLOSE (unit=Coordsunit, status='keep')
	END SUBROUTINE PrintCoords

! ***************
! SUBROUTINE PrintGaussPan
! ***************
! DESCRIPTION: 	prints out the file "PanDetail.dat"
	SUBROUTINE PrintGaussPan(GaussPan, numGaussPan)
		INTEGER, INTENT(IN)	:: numGaussPan
		TYPE(PANEL), DIMENSION(numGaussPan), INTENT(IN)	:: GaussPan
		INTEGER	:: i, k, tempPanunit, panelCount

		CONTINUE
		tempPanunit=13

		PRINT *, '  -Writing "PanDetail.dat" '
		OPEN(unit=tempPanunit, file='PanDetail.dat', status='replace')
		WRITE(unit=tempPanunit, fmt='(a)') 'VARIABLES="R [cm]"'
		WRITE(unit=tempPanunit, fmt='(a)') '"Z [cm]"'
		WRITE(unit=tempPanunit, fmt='(a)') '"Potential"'
		WRITE(unit=tempPanunit, fmt='(a)') '"Pt type"'
			!WRITE(unit=tempPanunit, fmt='(a)') '"Pt str"'		! only if no 'k' loop
			!WRITE(unit=tempPanunit, fmt='(a)') '"length [cm]"' !
			!WRITE(unit=tempPanunit, fmt='(a)') '"rnorm"'		!
			!WRITE(unit=tempPanunit, fmt='(a)') '"znorm"'		!
		WRITE(unit=tempPanunit, fmt='(a)') '"Pan #"'
		WRITE(unit=tempPanunit, fmt='(a)') 'ZONE T="main"'

		panelCount=0
		DO i = 1,numGaussPan
				!WRITE(unit=tempPanunit, fmt='(3e12.4, i2, 4e13.4, i5)') GaussPan(i)%midr, &
				!	GaussPan(i)%midz, GaussPan(i)%midPot, INT(GaussPan(i)%type1), &
				!	GaussPan(i)%str, GaussPan(i)%length, GaussPan(i)%rnrm, GaussPan(i)%znrm, i

			DO k=1,npoints
				panelCount=panelCount+1
				WRITE(unit=tempPanunit, fmt='(3e14.6, 2i6)') GaussPan(i)%rpoints(k), &
					GaussPan(i)%zpoints(k),  GaussPan(i)%midPot, INT(GaussPan(i)%type1), panelCount
			ENDDO
		ENDDO
		CLOSE (unit=tempPanunit, status='keep')
	END SUBROUTINE PrintGaussPan
END MODULE MAT_Greens
