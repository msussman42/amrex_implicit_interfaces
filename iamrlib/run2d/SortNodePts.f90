MODULE SortNodePts

	use TreeCodeGlobal

	implicit none
! v6 14 Jan 05. Changing GradGradGreens to remove logical types
! v7 17 Jan	Much of the sorting will be done outside of the UMich code. Removing 
!			  duplicate functionality.
! v8 01 Feb Changes the electrodes from being the top boundary to a level 2x as 
!			  high to test how Neumann distance varies the answer
! v9 07 Mar v8 working. Debugging all Dirichlet for integration into GMRes

! used to restrict variables to prevent accidental calling outside
	private      ! everything defaults to 'private' except those expressly labeled otherwise.
	public			:: AddBCandAround, SortPts, GaussQuad_Data, TrapInit, &
		ExactPanelInteg, GaussInit, SortPtsNorm ! subroutines
	real(kind=r8), parameter :: Bigepsil=1e-3

contains 

! ***************
! SUBROUTINE SortPts[x]
! ***************
! DESCRIPTION: takes unsorted list of c0 surface points and produces
!   a sorted list, based on nearest neighbor
!
! INPUT: xarraySP, yarraySP
!
! OUTPUT: sorted xarraySP, yarraySP
!
! CALLING PROGRAM: CLSVOF[x]
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
!	19 Sep v2   changing to fit Sussman code. Using passed in grid value
!	05 Jan 2005 v4. Changing so that the left most point (domlo1,yleft) and the 
!				  right most point (domhi1,yright) are passed as parameters in 
!				  xarray, yarray. Adding 15 points to sizelevel to cover maxY, elecHeight
!	06 Jan v5	back to seven points, fewer variables passed
!	17 Jan v7	splitting sorting and filling arrays into separate sections

	SUBROUTINE SortPts( xarrayNew, yarrayNew, sizeLevel, xarrayOld, yarrayOld )

! incoming variables
		integer, intent(in)				    :: sizeLevel
		real(kind=r8), dimension(sizeLevel), intent(in) :: xarrayOld, yarrayOld

! outgoing variables
		real(kind=r8), dimension(sizeLevel), intent(out) :: xarrayNew, yarrayNew

! subroutine entirely variables
		integer, dimension(sizeLevel)		:: indexIn
		integer								:: nearPtNow, nearPtNext, i, j

		real(kind=r8)								:: outrangeVar, dist, distNeighbor

		continue
!*************
! main program commands begin

		outrangeVar=1e6
		do i=1,sizeLevel
			indexIn(i)=1	! forces the variable in the loop to be all points
		end do

! description of coordinate columns
! [1-2] coordinates to form a trapezoid around to test slanted panels
! [3] potential of boundary panel (if 1==[4]) or slope of panel (if 0==[4])
! [4] 1=dirch; 0=neumann boundary conditions

! left most point (domlo1,yleft) and the right most point (domhi1,yright) are 
! passed as parameters in xarray, yarray as the first and last points, respectively
		indexIn(1)=0; nearPtNow=1
		xarrayNew(1)=xarrayOld(1)
		yarrayNew(1)=yarrayOld(1)

! loops through and organizes 
		outer: do i=1,sizeLevel-1
			dist=outrangeVar; nearPtNext=0
			do j=1,sizeLevel
				if ( 1==indexIn(j) ) then
					distNeighbor = (xarrayOld(nearPtNow)-xarrayOld(j))**2 + &
						(yarrayOld(nearPtNow)-yarrayOld(j))**2

					if ( distNeighbor < dist ) then
						dist = distNeighbor
						nearPtNext=j
					endif
				endif
			end do

			xarrayNew(i+1)=xarrayOld(nearPtNext)
			yarrayNew(i+1)=yarrayOld(nearPtNext)
			indexIn(nearPtNext)=0; nearPtNow=nearPtNext
		end do outer 

		print *, 'Done with subroutine <SortPts>'

	end SUBROUTINE SortPts

! ***************
! SUBROUTINE SortPtsNorm[x]
! ***************
! DESCRIPTION: takes unsorted list of c0 surface points and produces
!   a sorted list, based on nearest neighbor and panel normals
!
! INPUT: xarraySP, yarraySP
!
! OUTPUT: sorted xarraySP, yarraySP
!
! CALLING PROGRAM: MAT_Greens:BoundaryForceCalc
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   22 Mar 2005 begun

	SUBROUTINE SortPtsNorm( xarrayNew, yarrayNew, sizeLevel, xarrayOld, &
			yarrayOld, xNorm, yNorm, TrackPts )

		! incoming variables
		integer, intent(in)							:: sizeLevel, TrackPts
		real(kind=r8), dimension(sizeLevel), intent(in)	:: xarrayOld, yarrayOld, xNorm, yNorm

		! outgoing variables
		real(kind=r8), dimension(sizeLevel), intent(out)	:: xarrayNew, yarrayNew

		! subroutine entirely variables
		integer										:: i, j, k
		integer, dimension(sizeLevel)				:: indxIn
		integer, dimension(sizeLevel, TrackPts)		:: NearNeighbors
		real(kind=r8)								:: outrangeVar, distNeighbor, AngleCutoff, &
			anglePre, anglePost, angleDiff
		real(kind=r8), dimension(sizeLevel)			:: TrackNewArray	
		real(kind=r8), dimension(sizeLevel, TrackPts) :: dist

		continue
		!*************
		! main program commands begin
		AngleCutoff=100	! max degrees that the normals can vary for sequential panel points
		outrangeVar=1e6; dist(:,:)=outrangeVar
		indxIn(:)=1; indxIn(1)=0

		! loops through and organizes 
		outerDo: do i=1,sizeLevel
			dist=outrangeVar
			do j=1,sizeLevel
				SelfSkipDo: if (i/=j) then
					distNeighbor = (xarrayOld(i)-xarrayOld(j))**2 + (yarrayOld(i)-yarrayOld(j))**2
					kDo: do k=1,TrackPts
						if (k>1) then
							if ( distNeighbor < dist(i,k) .AND. distNeighbor>dist(i,k-1) ) then 
								dist(i,k) = distNeighbor
								NearNeighbors(i,k) = j
							endif
						else
							if ( distNeighbor < dist(i,k) ) then 
								dist(i,k) = distNeighbor
								NearNeighbors(i,k) = j
							endif
						endif
					end do kDo
				endif SelfSkipDo
			end do
		end do outerDo 

		xarrayNew(1)=xarrayOld(1); yarrayNew(1)=yarrayOld(1); TrackNewArray(1)=1
		sortingDo: do i=2,sizeLevel
			if (64==i) then
				print *, 'testing'
			endif

			j=1
			do while ( .TRUE. )
				anglePre=atan2(yNorm(i),xNorm(i))
				anglePost=atan2(yNorm(NearNeighbors(i-1,j)),xNorm(NearNeighbors(i-1,j)))
				angleDiff=abs(anglePre-anglePost)*180.0/3.14159
				if ( (indxIn(NearNeighbors(i-1,j)) ==1 .AND. angleDiff<AngleCutoff) .OR. j==TrackPts ) then
					exit
				endif

				j=j+1
			enddo
			indxIn(NearNeighbors(i-1,j))=0
			xarrayNew(i)=xarrayOld( NearNeighbors(i-1,j) )
			yarrayNew(i)=yarrayOld( NearNeighbors(i-1,j) )
			TrackNewArray(i)=NearNeighbors(i-1,j)
		enddo sortingDo

		! temp command
		i=floor(TrackNewArray(1))
		print *, 'Done with subroutine <SortPtsNorm>'

	end SUBROUTINE SortPtsNorm

! ***************
! SUBROUTINE AddBCandAround[x]
! ***************
! DESCRIPTION: adds 3rd column (boundary type) and 4th column (boundary value) to
!	pre-sorted data
!
! INPUT: LevelPtsX, LevelPtsY 
!
! OUTPUT: fleshed out Coords
!
! CALLING PROGRAM: MAT_Greens[x]
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   17 Jan 05 split from SortPts
!	01 Feb v8 changes the electrodes from being the top boundary to a level 2x as 
!				high to test how Neumann distance varies the answer

	SUBROUTINE AddBCandAround( Coords, sizeLevel, xarraySP, yarraySP, elecHeight, elecWidth, &
			elecLedge, elecRedge)

		! incoming variables
		integer, intent(in)				    :: sizeLevel
		real(kind=r8), dimension(sizeLevel), intent(in) :: xarraySP, yarraySP
		real(kind=r8), intent(in)			:: elecHeight, elecWidth

		! outgoing variables
		real(kind=r8), intent(out)			:: elecLedge, elecRedge
		real(kind=r8), dimension(sizeLevel+18,4), intent(inout) :: Coords

		! subroutine entirely variables
		integer								:: i, CaseShape
		real(kind=r8)						:: maxX, minX, elecThick, JustAbove, elecPotential, &
			minY, maxY

		continue
		!*************
		! main program commands begin
		i=sizeLevel

		!maxX = xarraySP(sizeLevel); minX = xarraySP(1)
		maxX = 2.; minX=-0.5
		maxY = 2.; minY = yarraySP(1)

		elecThick = 0.1; JustAbove=0.01;
		elecPotential= 8e10; elecLedge=1.; elecRedge=0.5

  		! define the fluid coordinates with xarray, yarray and panels 
  		!   with potential -2.0 and dirichlet boundary conditions(==1)
		Coords(1:sizeLevel,1) = xarraySP(:)
		Coords(1:sizeLevel,2) = yarraySP(:)
		Coords(1:sizeLevel-1,3) = -2.0
		Coords(1:sizeLevel-1,4) = 1

		CaseShape=4			!=1, single box, simple pot; =2 single, real pot; 
							!=3, double high, real pot
		select case(CaseShape)
		case(1)				! test case around single box
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
			Coords(i+18,1)= minX;					Coords(i+18,2)= minY

			! debugging, test shape, force to all Dirichlet to check code bit by bit
			Coords(i,3)  = elecPotential;			Coords(i,4)  = 1
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

		case(2) 			!=2 single, real pot
			Coords(i+1,1)= maxX;					Coords(i+1,2)= minY
			Coords(i+2,1)= maxX;					Coords(i+2,2)= 0.5*elecHeight
			Coords(i+3,1)= maxX;					Coords(i+3,2)= 0.6*elecHeight
			Coords(i+4,1)= maxX;					Coords(i+4,2)= 0.7*elecHeight
			Coords(i+5,1)= maxX;					Coords(i+5,2)= elecHeight
			Coords(i+6,1)= elecLedge+elecWidth;		Coords(i+6,2)= elecHeight
			Coords(i+7,1)= elecLedge+elecWidth;		Coords(i+7,2)= elecHeight-elecThick
			Coords(i+8,1)= elecLedge;				Coords(i+8,2)= elecHeight-elecThick
			Coords(i+9,1)= elecLedge;				Coords(i+9,2)= elecHeight
			Coords(i+10,1)= elecRedge;				Coords(i+10,2)= elecHeight
			Coords(i+11,1)= elecRedge;				Coords(i+11,2)= elecHeight-elecThick
			Coords(i+12,1)= elecRedge-elecWidth;	Coords(i+12,2)= elecHeight-elecThick
			Coords(i+13,1)= elecRedge-elecWidth;	Coords(i+13,2)= elecHeight
			Coords(i+14,1)= minX;					Coords(i+14,2)= elecHeight
			Coords(i+15,1)= minX;					Coords(i+15,2)= 0.7*elecHeight
			Coords(i+16,1)= minX;					Coords(i+16,2)= 0.6*elecHeight
			Coords(i+17,1)= minX;					Coords(i+17,2)= 0.5*elecHeight
			Coords(i+18,1)= minX;					Coords(i+18,2)= minY

			! debugging, single test shape, real pot
			Coords(i,3)  = 0;						Coords(i,4)  = 0
			Coords(i+1,3)= 0;						Coords(i+1,4)= 0
			Coords(i+2,3)= 0;						Coords(i+2,4)= 0
			Coords(i+3,3)= 0;						Coords(i+3,4)= 0
			Coords(i+4,3)= 0;						Coords(i+4,4)= 0
			Coords(i+5,3)= 0;						Coords(i+5,4)= 0
			Coords(i+6,3)= elecPotential;			Coords(i+6,4)= 1
			Coords(i+7,3)= elecPotential;			Coords(i+7,4)= 1
			Coords(i+8,3)= elecPotential;			Coords(i+8,4)= 1
			Coords(i+9,3)= 0;						Coords(i+9,4)= 0
			Coords(i+10,3)= elecPotential;			Coords(i+10,4)= 1
			Coords(i+11,3)= elecPotential;			Coords(i+11,4)= 1
			Coords(i+12,3)= elecPotential;			Coords(i+12,4)= 1
			Coords(i+13,3)= 0;						Coords(i+13,4)= 0
			Coords(i+14,3)= 0;						Coords(i+14,4)= 0
			Coords(i+15,3)= 0;						Coords(i+15,4)= 0
			Coords(i+16,3)= 0;						Coords(i+16,4)= 0
			Coords(i+17,3)= 0;						Coords(i+17,4)= 0
			Coords(i+18,3)= 0;						Coords(i+18,4)= 0

		case(3)				! real case with curved electrodes and 2x higher full "T"
			Coords(i+1,1)= maxX;					Coords(i+1,2)= minY
			Coords(i+2,1)= maxX;					Coords(i+2,2)= elecHeight
			Coords(i+3,1)= elecLedge+elecWidth;		Coords(i+3,2)= elecHeight
			Coords(i+4,1)= elecLedge+elecWidth;		Coords(i+4,2)= elecHeight-elecThick
			Coords(i+5,1)= elecLedge;				Coords(i+5,2)= elecHeight-elecThick
			Coords(i+6,1)= elecLedge;				Coords(i+6,2)= elecHeight+JustAbove
			Coords(i+7,1)= elecLedge+elecWidth;		Coords(i+7,2)= elecHeight+JustAbove
			Coords(i+8,1)= maxX;					Coords(i+8,2)= elecHeight+JustAbove
			Coords(i+9,1)= maxX;					Coords(i+9,2)= maxY
			Coords(i+10,1)= minX;					Coords(i+10,2)= maxY
			Coords(i+11,1)= minX;					Coords(i+11,2)= elecHeight+JustAbove
			Coords(i+12,1)= elecRedge-elecWidth;	Coords(i+12,2)= elecHeight+JustAbove
			Coords(i+13,1)= elecRedge;				Coords(i+13,2)= elecHeight+JustAbove
			Coords(i+14,1)= elecRedge;				Coords(i+14,2)= elecHeight-elecThick
			Coords(i+15,1)= elecRedge-elecWidth;	Coords(i+15,2)= elecHeight-elecThick
			Coords(i+16,1)= elecRedge-elecWidth;	Coords(i+16,2)= elecHeight
			Coords(i+17,1)= minX;					Coords(i+17,2)= elecHeight
			Coords(i+18,1)= minX;					Coords(i+18,2)= minY

			Coords(i,3)  = 0;						Coords(i,4)  = 0
			Coords(i+1,3)= 0;						Coords(i+1,4)= 0
			Coords(i+2,3)= 0;						Coords(i+2,4)= 0
			Coords(i+3,3)= elecPotential;			Coords(i+3,4)= 1
			Coords(i+4,3)= elecPotential;			Coords(i+4,4)= 1
			Coords(i+5,3)= elecPotential;			Coords(i+5,4)= 1
			Coords(i+6,3)= -elecPotential;			Coords(i+6,4)= 1
			Coords(i+7,3)= 0;						Coords(i+7,4)= 0
			Coords(i+8,3)= 0;						Coords(i+8,4)= 0
			Coords(i+9,3)= 0;						Coords(i+9,4)= 0
			Coords(i+10,3)= 0;						Coords(i+10,4)= 0
			Coords(i+11,3)= 0;						Coords(i+11,4)= 0
			Coords(i+12,3)= -elecPotential;			Coords(i+12,4)= 1
			Coords(i+13,3)= elecPotential;			Coords(i+13,4)= 1
			Coords(i+14,3)= elecPotential;			Coords(i+14,4)= 1
			Coords(i+15,3)= elecPotential;			Coords(i+15,4)= 1
			Coords(i+16,3)= 0;						Coords(i+16,4)= 0
			Coords(i+17,3)= 0;						Coords(i+17,4)= 0
			Coords(i+18,3)= 0;						Coords(i+18,4)= 0

		case(4)	! real case with curved electrodes and 2x higher partial "T"
			Coords(i+1,1)= maxX;					Coords(i+1,2)= minY
			Coords(i+2,1)= maxX;					Coords(i+2,2)= 0.8*elecHeight
			Coords(i+3,1)= maxX;					Coords(i+3,2)= elecHeight
			Coords(i+4,1)= elecLedge+elecWidth;		Coords(i+4,2)= elecHeight
			Coords(i+5,1)= elecLedge+elecWidth;		Coords(i+5,2)= elecHeight-elecThick
			Coords(i+6,1)= elecLedge;				Coords(i+6,2)= elecHeight-elecThick
			Coords(i+7,1)= elecLedge;				Coords(i+7,2)= elecHeight+elecThick
			Coords(i+8,1)= elecLedge+elecWidth;		Coords(i+8,2)= elecHeight+elecThick
			Coords(i+9,1)= elecLedge+elecWidth;		Coords(i+9,2)= maxY
			Coords(i+10,1)= elecRedge-elecWidth;	Coords(i+10,2)= maxY
			Coords(i+11,1)= elecRedge-elecWidth;	Coords(i+11,2)= elecHeight+elecThick
			Coords(i+12,1)= elecRedge;				Coords(i+12,2)= elecHeight+elecThick
			Coords(i+13,1)= elecRedge;				Coords(i+13,2)= elecHeight-elecThick
			Coords(i+14,1)= elecRedge-elecWidth;	Coords(i+14,2)= elecHeight-elecThick
			Coords(i+15,1)= elecRedge-elecWidth;	Coords(i+15,2)= elecHeight
			Coords(i+16,1)= minX;					Coords(i+16,2)= elecHeight
			Coords(i+17,1)= minX;					Coords(i+17,2)= 0.8*elecHeight
			Coords(i+18,1)= minX;					Coords(i+18,2)= minY

			Coords(i,3)  = 0;						Coords(i,4)  = 0
			Coords(i+1,3)= 0;						Coords(i+1,4)= 0
			Coords(i+2,3)= 0;						Coords(i+2,4)= 0
			Coords(i+3,3)= 0;						Coords(i+3,4)= 0
			Coords(i+4,3)= elecPotential;			Coords(i+4,4)= 1
			Coords(i+5,3)= elecPotential;			Coords(i+5,4)= 1
			Coords(i+6,3)= elecPotential;			Coords(i+6,4)= 1
			Coords(i+7,3)= -elecPotential;			Coords(i+7,4)= 1
			Coords(i+8,3)= 0;						Coords(i+8,4)= 0
			Coords(i+9,3)= 0;						Coords(i+9,4)= 0
			Coords(i+10,3)= 0;						Coords(i+10,4)= 0
			Coords(i+11,3)= -elecPotential;			Coords(i+11,4)= 1
			Coords(i+12,3)= elecPotential;			Coords(i+12,4)= 1
			Coords(i+13,3)= elecPotential;			Coords(i+13,4)= 1
			Coords(i+14,3)= elecPotential;			Coords(i+14,4)= 1
			Coords(i+15,3)= 0;						Coords(i+15,4)= 0
			Coords(i+16,3)= 0;						Coords(i+16,4)= 0
			Coords(i+17,3)= 0;						Coords(i+17,4)= 0
			Coords(i+18,3)= 0;						Coords(i+18,4)= 0
		case default
			stop 'Err CaseShape'
		end select

		! now change the coordinates back to MKS units
		!Coords(:,1:2) = Coords(:,1:2)
		!elecLedge=elecLedge; elecRedge=elecRedge
		print *, 'Done with subroutine <AddBCandAround>'

	end SUBROUTINE AddBCandAround

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

		integer, intent(inout)	:: NumCoords
		integer					:: XDataPts
		real(kind=r8)			:: a

continue
		NumCoords=0; XDataPts=10
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

	SUBROUTINE GaussQuad_data ( PtlocationX, PtlocationY, Ptweigh, PtsPerPanel, &
			PanBeginX, PanLenX, PanBeginY, PanLenY, j )

! incoming variable
		integer, intent(in)				:: PtsPerPanel, j
		real(kind=r8), intent(in)		:: PanBeginX, PanLenX, PanBeginY, PanLenY

! outgoing variables
		real(kind=r8), intent(inout)	:: PtlocationX, PtlocationY, Ptweigh

! subroutine entirely variables
		real(kind=r8)					:: PanEndX, PanEndY
		real(kind=r8), dimension(PtsPerPanel) :: Mlocation, Mweight

continue

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
			Mlocation(10)= 0.97390653; Mweight(10)=0.06667134
			Mlocation(9)= 0.86506337; Mweight(9)= 0.14945135  
			Mlocation(8)= 0.67940957; Mweight(8)= 0.21908636
			Mlocation(7)= 0.43339539; Mweight(7)= 0.26926672  
			Mlocation(6)= 0.14887434; Mweight(6)= 0.29552422
			Mlocation(5)= -0.14887434; Mweight(5)= 0.29552422  
			Mlocation(4)= -0.43339539; Mweight(4)= 0.26926672
			Mlocation(3)= -0.67940957; Mweight(3)= 0.21908636  
			Mlocation(2)= -0.86506337; Mweight(2)= 0.14945135
			Mlocation(1)= -0.97390653; Mweight(1)= 0.06667134 
		case (12)
			Mlocation(12)= 0.981560634246732; Mweight(12)= 0.0471753363864754
			Mlocation(11)= 0.904117256370452; Mweight(11)= 0.1069393259953637
			Mlocation(10)= 0.769902674194317; Mweight(10)= 0.1600783285433586
			Mlocation(9)=  0.587317954286614; Mweight(9)=  0.2031674267230672
			Mlocation(8)=  0.367831498998180; Mweight(8)=  0.2334925365383534
			Mlocation(7)=  0.125233408511468; Mweight(7)=  0.2491470458134027
			Mlocation(6)= -0.125233408511468; Mweight(6)=  0.2491470458134027
			Mlocation(5)= -0.367831498998180; Mweight(5)=  0.2334925365383534
			Mlocation(4)= -0.587317954286614; Mweight(4)=  0.2031674267230672
			Mlocation(3)= -0.769902674194317; Mweight(3)=  0.1600783285433586
			Mlocation(2)= -0.904117256370452; Mweight(2)=  0.1069393259953637
			Mlocation(1)= -0.981560634246732; Mweight(1)=  0.0471753363864754
		case default
			stop 'Data not entered for this number of points per panel.'
		end select PtsPan

		PtlocationX=(PanEndX+PanBeginX)/2 + Mlocation(j)*(PanEndX-PanBeginX)/2
		PtlocationY=(PanEndY+PanBeginY)/2 + Mlocation(j)*(PanEndY-PanBeginY)/2
		Ptweigh=Mweight(j)

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
		integer, intent(in)				:: PanFromTo, flagTrap, k, Mpts
		real(kind=r8), intent(in)		:: PanLen, delX, delY, normXj, normYj, normXi, normYi, PanWeight

		! outgoing variables
		real(kind=r8), intent(inout)	:: GreenPart

		! subroutine entirely variables
		real(kind=r8)					:: HalfLen							

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

		real(kind=r8), intent(in)	:: zzz, aaa
		real(kind=r8)				:: Greens

		continue

		Greens = -log( (zzz)**2. + (aaa)**2. ) / (4.*pi)
	end FUNCTION Greens

    !-------------------------
	FUNCTION GradGreens( delX, delY, normX, normY )
    ! grad(G) subfunction

		real(kind=r8)				:: GradGreens
		real(kind=r8), intent(in)	:: delX, delY, normX, normY

		continue
		!   d/dn {ln (r_ij) } = grad_i {ln (r_ij) } . n_i
		! = (nx*dG/dx+ny*dG/dy)*dS

		GradGreens = ( delX*normX + delY*normY ) / ( ( delX**2 + delY**2 )*(2.*pi) )
	end FUNCTION GradGreens

!-------------------------
	FUNCTION GradGradGreens(delX, delY, normXj, normYj, normXi, normYi)
!neumann-dirichlet grad.grad subfunction
! d/dn_i (d/dn_j {ln (r_ij) }) = grad_i (grad_j {ln (r_ij) } . n_j) . n_i

		real(kind=r8), intent(in)	:: delX, delY, normXj, normYj, normXi, normYi
		real(kind=r8)				:: GradGradGreens, R, gamma, tempXY

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

	SUBROUTINE ExactPanelInteg ( GreenInt, xBeg, xEnd, a, AtanOk )
! incoming variables
		integer, intent(in)		:: AtanOk	! switch, =0 if a=0, =1 otherwise
		real(kind=r8), intent(in)		:: xBeg, xEnd, a
		real(kind=r8), intent(out)		:: GreenInt

		real(kind=r8)					:: tempxBeg, tempxEnd

continue
! make sure that the 'X' term is the one that you want to integrate with respect to.
! int(ln(x^2+a^2) dx = x ln(x^2+a^2) -2x +2a tan^-1(x/a)
		if (0==xBeg) then
			tempxBeg=1
		else
			tempxBeg=xBeg
		endif
		if (0==xEnd) then
			tempxEnd=1
		else
			tempxEnd=xEnd
		endif

		if (1==AtanOk) then
			GreenInt = (xEnd*log(xEnd*xEnd + a*a) - 2.*xEnd + 2.*a*atan2(a,xEnd))/(4.*pi) - &
				((xBeg*log(xBeg*xBeg + a*a) - 2.*xBeg + 2.*a*atan2(a,xBeg))/(4.*pi))
		else
			GreenInt = (xEnd*log(tempxEnd*tempxEnd) - 2.*xEnd - xBeg*log(tempxBeg*tempxBeg) + 2.*xBeg)/(4.*pi)
		endif

	END SUBROUTINE ExactPanelInteg

! ***************
! SUBROUTINE GaussInit
! ***************
! DESCRIPTION: integrates a function using Gaussian quadrature
!
! INPUT: panel length, weight, value
!
! OUTPUT: approximate numerically integrated value along panel
!
! CALLING PROGRAM: TreeOperations:TaylorCoeff
!
! Primary author: Anton VanderWyst, 
!	University of Michigan 
!
! Version history:
!   21 Mar 2005 begun
	SUBROUTINE GaussInit ( GaussPart, PanWeight, PanLen, k, Mpts, IntVariable )

! incoming variables
		integer, intent(in)		:: k, Mpts
		real(kind=r8), intent(in)		:: PanLen, PanWeight, IntVariable

! outgoing variables
		real(kind=r8), intent(out)		:: GaussPart

! subroutine entirely variables
		real(kind=r8)					:: HalfLen							

continue

		if ( 1==k .OR. Mpts==k ) then
			HalfLen = 0.5
		else
			HalfLen = 1.0 
		end if

		GaussPart = PanLen * HalfLen * PanWeight * IntVariable

	end SUBROUTINE GaussInit

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

	SUBROUTINE InSideYesNo( CrossCount, CrossYesNo, IntersecX, IntersecY, CrossDir, &
			Ax, Ay, Bx, By, Px, Py )

! incoming variables
		real(kind=r8), intent(in)		:: Ax, Ay, Bx, By, Px, Py
character, intent(in)	:: CrossDir

! outgoing variables
		integer, intent(inout)	:: CrossYesNo ! formerly logical
		integer, intent(inout)	:: CrossCount
		real(kind=r8), intent(inout)		:: IntersecX, IntersecY

! subroutine entirely variables
		real(kind=r8)					:: delBAx, delBAy, t

		continue

		CrossYesNo=0
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
			CrossCount = CrossCount+1; CrossYesNo=1
		endif

	end SUBROUTINE InSideYesNo
end MODULE SortNodePts
