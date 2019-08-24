! ***************
! Program Diagmatrix[x]
! ***************
! DESCRIPTION: main calling program to test call 2d axisymmetric BEM code
!
! INPUT: <na>
!
! OUTPUT: green's for given matrix
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   19 June 2005 begun

PROGRAM MAIN
	USE global, ONLY 	: eps, npoints, pi, SmallInitAry, InitArySize
	USE MAT_Greens,ONLY : BoundaryForceCalc

	IMPLICIT NONE

  	! PARAMETERS
	INTEGER, PARAMETER	:: Runit=10, Phiunit=11, CaseLoad=5, LSGridMaxZ=257, LSGridMaxR=129

	! regular variables
	INTEGER				:: i, j, NumBlobs, MxPt_blob, is_symmetric
	INTEGER, DIMENSION(SmallInitAry) 		:: NumPt_blob
	REAL*8				:: rlo, zlo, rhi, zhi, dr, dz, elecHeight, r, z, &
		rblob, zblob, radblob2, CosConst, baseExpZ, relsize, ptime
	REAL*8, DIMENSION(SmallInitAry)			:: blob_rc, blob_zc
	REAL*8,DIMENSION(LSGridMaxR,LSGridMaxZ) :: phi2
	REAL*8, DIMENSION(SmallInitAry, InitArySize) :: radBlob, thetaBlob, forceBlob

	CONTINUE
	! simulation values
	rlo=0.0; rhi=0.25; zlo=0.00; zhi=0.5
	rblob=0.5; zblob=0.5; radblob2=0.25
	elecHeight=1.5; is_symmetric=0; ptime=0.

	IF (1==mod(LSGridMaxZ,2)) THEN	! comparing to Mark's phi values
		dr=(rhi-rlo)/(LSGridMaxR-1.); dz=(zhi-zlo)/(LSGridMaxZ-1.)
		relsize=1.*(zhi-zlo)/(rhi-rlo) - (LSGridMaxZ-1.)/(LSGridMaxR-1.)
	ELSE							! running my test cases
		dr=(rhi-rlo)/(LSGridMaxR); dz=(zhi-zlo)/(LSGridMaxZ)
		relsize=1.*(zhi-zlo)/(rhi-rlo) - 1.*(LSGridMaxZ)/(LSGridMaxR)
	ENDIF

	IF ( ABS(relsize) > 1e-5) THEN
		PRINT *, 'AxiSym5: Grid axis relative size mismatch.', relsize
		STOP		
	ENDIF

	SELECT CASE(CaseLoad)
	CASE(1)	! flat line, for debugging. Use with CaseShape=1

		DO i =1,LSGridMaxR
			DO j=1,LSGridMaxZ
				IF (j>1 .OR. i>LSGridMaxR/2) THEN
					phi2(i,j)=-0.1
				ELSE
					phi2(i,j)=0.05
				ENDIF
			ENDDO
		ENDDO	
		CLOSE (UNIT=Runit, STATUS='keep', ERR=993)
	CASE(2)	! fixed data set, normal spread shape
		! set LSGridMaxY=400, LSGridMaxX=200
		! load data
		OPEN(unit=Runit, file='vis_20050705_624.dat', ERR=990)
		DO i =1,LSGridMaxR
			DO j=1,LSGridMaxZ
				READ(UNIT=Runit, FMT=11, ERR=991) NumBlobs, NumBlobs, r, phi2(i,j), r
			ENDDO
		ENDDO	
		CLOSE (UNIT=Runit, STATUS='keep', ERR=993)
	CASE(3)	! smooth data set, no bumps
		PRINT *, '  -Writing AxiSym5:"Phi.dat" '
		OPEN(unit=Phiunit, file='Phi.dat', status='replace')
		WRITE(unit=Phiunit, fmt='(a)') 'VARIABLES="R [cm]"'
		WRITE(unit=Phiunit, fmt='(a)') '"Z [cm]"'
		WRITE(unit=Phiunit, fmt='(a)') '"phi"'
		WRITE(unit=Phiunit, fmt='(a)') 'ZONE T="main2"'
		WRITE(unit=Phiunit, fmt='(a,i4,a,i4,a)') 'I=', LSGridMaxZ, ', J=', LSGridMaxR, ', K=1, F=POINT'

		CosConst = 2.0*pi/rblob
		DO i=1,LSGridMaxR
			r = (i+0.5)*dr + rlo
			DO j=1,LSGridMaxZ
				z = (j+0.5)*dz + zlo
				!phi2(i,j)=zblob+radblob2*cos(r*CosConst)-z
				phi2(i,j)=-( (r-0.)**2 + (z-0.4)**2 - 0.1**2)

				!WRITE(unit=Phiunit, fmt='(3e14.6)') r, z, phi2(i,j)
			ENDDO
		ENDDO
		CLOSE (unit=Phiunit, status='keep')
	CASE(4)	! no data set - load them in BFC
		! do nothing
	CASE(5)	! full Gaussian, to match paper 'Formation of the Taylor cone on
		!  the surface of liquid metal...'
		PRINT *, '  -Writing AxiSym5:"Phi.dat" '
		OPEN(unit=Phiunit, file='Phi.dat', status='replace')
		WRITE(unit=Phiunit, fmt='(a)') 'VARIABLES="R [cm]"'
		WRITE(unit=Phiunit, fmt='(a)') '"Z [cm]"'
		WRITE(unit=Phiunit, fmt='(a)') '"phi"'
		WRITE(unit=Phiunit, fmt='(a)') 'ZONE T="main3"'
		WRITE(unit=Phiunit, fmt='(a,i4,a,i4,a)') 'I=', LSGridMaxZ, ', J=', LSGridMaxR, ', K=1, F=POINT'

		dr=(0.4-0.0)/(LSGridMaxR-1)
		dz=(1.0-0.0)/(LSGridMaxZ-1)

		DO i=1,LSGridMaxR
			r = (i-1)*dr
			baseExpZ = 0.2*EXP(-22.9*r*r)	! for h=lamda test case
			!baseExpZ = 0.02*EXP(-17.33*r*r)	! for h=0.1 lamda test case
			DO j=1,LSGridMaxZ
				z = (j-1)*dz
				phi2(i,j)= baseExpZ - z

				!WRITE(unit=Phiunit, fmt='(3e14.6)') r, z, phi2(i,j)
			ENDDO
		ENDDO
		CLOSE (unit=Phiunit, status='keep')
	CASE DEFAULT
		STOP 'CallShading: err CaseLoad'
	END SELECT

	PRINT *, 'Done loading variables. Calling <BoundaryForceCalc>'

	CALL BoundaryForceCalc( NumBlobs, blob_rc, blob_zc, NumPt_blob, MxPt_blob, radBlob, &
		thetaBlob, forceBlob, phi2, LSGridMaxR, LSGridMaxZ, is_symmetric, rlo, zlo, &
		dr, dz, elecHeight, ptime )	! formal full version

	PRINT *, 'Done with AxiSym5'
	GOTO 500

	! FORMAT listings
	11 FORMAT(2i5, 3f20.10)

	! ERROR listings
	990 print *, 'Error open visarray'
	991 print *, 'Error read visarray'
	993 print *, 'Error close visarray'

500 END PROGRAM MAIN
