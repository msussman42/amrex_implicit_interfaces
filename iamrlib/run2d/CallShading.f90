! ***************
! Program CallShading[x]
! ***************
! DESCRIPTION: main calling program to test call my shading subroutines
!
! INPUT: <na>
!
! OUTPUT: inside/outside data points and ordered lines
!
! Primary author: Anton VanderWyst
! University of Michigan
! 1320 Beal Ave.
! Ann Arbor, MI 48109
! antonv@umich.edu 
!
! Version history:
!   24 Jul 2005 begun

PROGRAM CallShading

 USE ShadingGlobal, ONLY : r8, pi, InitArySize, SmallInitAry
 USE Shading, ONLY : BoundaryForceCalc

   ! PARAMETER
 INTEGER, PARAMETER   :: Xunit=10, Phiunit=11, CaseLoad=3, LSGridMaxY=400, LSGridMaxX=150, &
  is_symmetric=0, sizetemp=400
 INTEGER      :: i, j, iNode, jNode, NumBlobs, MxPt_blob
 INTEGER, DIMENSION(SmallInitAry)     :: NumPt_blob
 REAL(KIND=r8)    :: blah, xlo, ylo, xhi, yhi, dx, dy, elecHeight, x, y
 !REAL(KIND=r8), DIMENSION(:), ALLOCATABLE   :: tempa, tempb
 REAL(KIND=r8), DIMENSION(SmallInitAry)   :: blob_xc, blob_yc
 REAL(KIND=r8), DIMENSION(LSGridMaxX,LSGridMaxY) :: phi2, phiblob2
 REAL(KIND=r8), DIMENSION(SmallInitAry, InitArySize) :: radBlob, thetaBlob, forceBlob

 CONTINUE
 ! simulation values
 xlo=-0.375; xhi=0.375; ylo=0.0; yhi=1.5
 elecHeight=1.8
 dx=(xhi-xlo)/(LSGridMaxX-1); dy=(yhi-ylo)/(LSGridMaxY-1)

 SELECT CASE(CaseLoad)
 CASE(1) ! fixed data set, extreme spread shape
  ! splits shape incorrectly at edges
  ! set LSGridMaxY=300, LSGridMaxX=150
  ! load data
  OPEN(unit=Xunit, file='vis_20050706_1070.dat', ERR=990)
  DO i =1,LSGridMaxX
   DO j=1,LSGridMaxY
    READ(UNIT=Xunit, FMT=11, ERR=991) iNode, jNode, blah, phi2(i,j), blah
   ENDDO
  ENDDO 
  CLOSE (UNIT=Xunit, STATUS='keep', ERR=993)
 CASE(2) ! fixed data set, normal spread shape
  ! set LSGridMaxY=400, LSGridMaxX=200
  ! load data
  OPEN(unit=Xunit, file='vis_20050705_624.dat', ERR=990)
  DO i =1,LSGridMaxX
   DO j=1,LSGridMaxY
    READ(UNIT=Xunit, FMT=11, ERR=991) iNode, jNode, blah, phi2(i,j), blah
   ENDDO
  ENDDO 
  CLOSE (UNIT=Xunit, STATUS='keep', ERR=993)
 CASE(3) ! smooth data set, no bumps
  OPEN(unit=Phiunit, file='Phi.dat', status='replace')
  WRITE(unit=Phiunit, fmt='(a)') 'VARIABLES="X [cm]"'
  WRITE(unit=Phiunit, fmt='(a)') '"Y [cm]"'
  WRITE(unit=Phiunit, fmt='(a)') '"phi"'
  WRITE(unit=Phiunit, fmt='(a)') 'ZONE T="main2"'
  WRITE(unit=Phiunit, fmt='(a,i3,a,i3,a)') 'I=', LSGridMaxY, ', J=', LSGridMaxX, ', K=1, F=POINT'

  DO i=1,LSGridMaxX
   x = (i+0.5)*dx + xlo
   DO j=1,LSGridMaxY
    y = (j+0.5)*dy + ylo
    phi2(i,j)= -0.5+0.25*cos(2*pi*(x-0.5))+y

    WRITE(unit=Phiunit, fmt='(3e14.6)') x, y, phi2(i,j)
   ENDDO
  ENDDO
  CLOSE (unit=Phiunit, status='keep')
 CASE(4) ! no data set - load them in BFC
  ! set LSGridMaxY=400, LSGridMaxX=100; is_symmetric=1
  ! do nothing
 CASE DEFAULT
  STOP 'CallShading: err CaseLoad'
 END SELECT

 PRINT *, 'Done loading variables. Calling <BoundaryForceCalc>'

 IF ( CaseLoad /=4 ) THEN
  CALL BoundaryForceCalc( NumBlobs, blob_xc, blob_yc, NumPt_blob, MxPt_blob, radBlob, &
   thetaBlob, forceBlob, phi2, phiblob2, LSGridMaxX, LSGridMaxY, is_symmetric, xlo, ylo, &
   xhi, yhi, dx, dy, elecHeight) ! formal full version

  ! previous version, debugging work-around
  !ALLOCATE(tempa(sizetemp)); ALLOCATE(tempb(sizetemp))
  !CALL BoundaryForceCalc( xlo, xhi, yhi, xarray, yarray, sizetemp, qE )
  !DEALLOCATE(tempa); DEALLOCATE(tempb)
 ELSEIF ( 4==CaseLoad ) THEN
  !CALL BoundaryForceCalc( LSGridMaxX, LSGridMaxY, is_symmetric, xlo, ylo, dx, dy, elecHeight ) 
 ELSE
  PRINT *, 'CallShading Caseload BFC err'
  STOP
 ENDIF

 GOTO 999

 ! FORMAT listings
 11 FORMAT(2i5, 3f20.10)

 ! ERROR listings
 990 print *, 'Error open visarray'
 991 print *, 'Error read visarray'
 993 print *, 'Error close visarray'

999 end PROGRAM CallShading
