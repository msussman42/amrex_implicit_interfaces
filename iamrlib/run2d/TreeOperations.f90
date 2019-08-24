MODULE TreeOperations
  ! 03 Feb 2005 v2 TaylorCoefficients
  ! 08 Feb		v3 Changing to pointers
  !				v4 Changing all sectors to a single pass
  ! 15 Feb 2005 v5 Pulling Cluster and Taylor together for a speedup
  ! 25 Feb 2005 v6 Adding individual charges qi to get overall force. 
  ! 01 Mar 2005 v7 Forming trees for Mat_Greens, removing qi
  ! 15 Mar 2005 v8 Adding vi to Cluster, removing it from Taylor. Integrating Greens function 
  !					along panels instead of just point charges
  ! 21 Mar 2005 v9 TaylorCoeff change to integrate Greens function for each cluster. Does Taylor
  !					expansion once outside loop.

	use directoryTree, only	: start, add_node, change_node, find2, current, node, finish
	use SortNodePts, only	: GaussInit
	use TreeCodeGlobal

	implicit none

! used to restrict variables to prevent accidental calling outside
	private      ! everything defaults to 'private' except those expressly labeled otherwise.
	public		:: Resizer, SplitAndCount, NewParticleList, SizeTrees, &
		TaylorCoeff, ClusterMoment, FormTrees
	real*8, parameter	:: tooStretched = 1.414213 ! 2^0.5

contains 

! ***************
! Subroutine Resizer
! ***************
! DESCRIPTION: shrinks the box largest and smallest dimensions around points
!
! INPUT: Coordates, size array of Coordinates
!
! OUTPUT: max and min x,y size
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   12 Dec 2004 begun
!	08 Feb 2005 v3 changed to pointers and using directoryTree

	SUBROUTINE Resizer( TreeBranchNum, SizePts, NumPts, XPts, YPts, tempPartList )

! incoming variables
		integer, intent(in)						:: SizePts, TreeBranchNum, NumPts
		integer, dimension(NumPts), intent(in)  :: tempPartList
		real*8, dimension(NumPts), intent(in)	:: XPts, YPts

! outgoing variables
		type(node), pointer						:: TreeBranch

! subroutine only vars
		integer				:: startBranch
		real*8				:: MaxX, MinX, MaxY, MinY

		continue

		call find2(TreeBranchNum); TreeBranch=>current
		startBranch = TreeBranch%StartPosition	! the point in tempPartList we begin at

		MaxX=MAXVAL(XPts( tempPartList(startBranch:startBranch+SizePts-1) ))
		MinX=MINVAL(XPts( tempPartList(startBranch:startBranch+SizePts-1) ))
		MaxY=MAXVAL(YPts( tempPartList(startBranch:startBranch+SizePts-1) ))
		MinY=MINVAL(YPts( tempPartList(startBranch:startBranch+SizePts-1) ))

  ! change the node pointes for min/max x/y, and xcen/ycen
  ! ymin=1, ymax=2, xmin=3, xmax=4, xcen=5, ycen=6, aspectRatio=7
  ! NumPtInclude=8, WhichPtInclude=9, WhichChildren=10, 
  ! StartPosition=11, NumBoxes=12

		call change_node(TreeBranch, MinX, 3)
		call change_node(TreeBranch, MinY, 1)
		call change_node(TreeBranch, MaxX, 4)
		call change_node(TreeBranch, MaxY, 2)
		call change_node(TreeBranch, (MaxX+MinX)/2., 5)
		call change_node(TreeBranch, (MaxY+MinY)/2., 6)

		if (MaxX-MinX /= 0) then
			call change_node(TreeBranch, abs( (MaxY-MinY)/(MaxX-MinX) ), 7)
		else
			call change_node(TreeBranch, 0._r8, 7)
		endif

	END SUBROUTINE Resizer

! ***************
! Subroutine SplitAndCount
! ***************
! DESCRIPTION: divides the tree area into smaller boxes and reallocates the included
!	points into one of those new boxes
!
! INPUT: Treebranch, NumPts
!
! OUTPUT: new, sub-Treebranch, sub-NumPts
!
! CALLING PROGRAM: TreeOperations:SizeTrees
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   14 Dec 2004 begun

	SUBROUTINE SplitAndCount( tempPartList, TreeBranch, TreeBranchNum, XPts, &
			YPts, SizePts, NumPts, whichgL)

! incoming variables
		integer									:: TreeBranchNum, SizePts, NumPts, whichgL
		real*8, dimension(NumPts), intent(in)	:: XPts, YPts
		type(node), pointer						:: TreeBranch

! outgoing variables
		integer, dimension(NumPts),intent(inout):: tempPartList

! subroutine only vars
		integer				:: i, j, i1, NPts1, NPts2, NPts3, NPts4, sizeQuad, startBranch
		integer, dimension(:), allocatable	:: Quad1, Quad2, Quad3, Quad4
		real*8				:: divX, divY, BXMax, BYMax, BXMin, BYMin, tempAspect
		type(node), pointer	:: NB1, NB2, NB3, NB4

		continue
		NPts1=0; NPts2=0; NPts3=0; NPts4=0

		BXMax = TreeBranch%xmax; BXMin = TreeBranch%xmin
		BYMax = TreeBranch%ymax; BYMin = TreeBranch%ymin
		sizeQuad = TreeBranch%NumPtInclude		! max size of the array
		startBranch = TreeBranch%StartPosition	! the point in tempPartList we begin at

! check to make sure boxes are close to rectangular
		tempAspect=TreeBranch%aspectRatio
		if ( tempAspect > tooStretched) then	! too tall, only divide on Y
			divX = BXMax
			divY = 1.*( BYMax+BYMin ) /2
		elseif(tempAspect<1./tooStretched) then ! too wide, only divide on X
			divX = 1.*( BXMax+BXMin ) /2
			divY = BYMax
		else
			divX = 1.*( BXMax+BXMin ) /2
			divY = 1.*( BYMax+BYMin ) /2
		endif

		allocate(Quad1(sizeQuad)); allocate(Quad2(sizeQuad))
		allocate(Quad3(sizeQuad)); allocate(Quad4(sizeQuad))
		Quad1(:)=0; Quad2(:)=0; Quad3(:)=0; Quad4(:)=0

		if (SizePts>PtsPerBoxMax) then
			do i=1,sizeQuad
				i1=tempPartList(startBranch)
				if (XPts(i1) <= divX .AND. YPts(i1) >= divY ) then		! upper left
					NPts1 = NPts1+1; Quad1(NPts1)=i1
					startBranch=startBranch+1
				elseif (XPts(i1) > divX .AND. YPts(i1) > divY ) then	! upper right
					NPts2 = NPts2+1; Quad2(NPts2)=i1
					startBranch=startBranch+1
				elseif (XPts(i1) < divX .AND. YPts(i1) < divY ) then	! lower left
					NPts3 = NPts3+1; Quad3(NPts3)=i1
					startBranch=startBranch+1
				elseif (XPts(i1) >= divX .AND. YPts(i1) <= divY) then	! lower right
					NPts4 = NPts4+1; Quad4(NPts4)=i1
					startBranch=startBranch+1
				else
					stop 'Error boxes'
				endif
			enddo
		endif

		j=0
		if (NPts1 /=0) then
			j=j+1
		endif
		if (NPts2 /=0) then
			j=j+1
		endif
		if (NPts3 /=0) then
			j=j+1
		endif
		if (NPts4 /=0) then
			j=j+1
		endif
		allocate(TreeBranch%WhichChildren(j))

  ! ymin=1, ymax=2, xmin=3, xmax=4, xcen=5, ycen=6, aspectRatio=7
  ! NumPtInclude=8, WhichPtInclude=9, WhichChildren=10, 
  ! StartPosition=11, NumBoxes=12
		j=0
		if (NPts1 > 0) then
			j=j+1; whichgL=whichgL+1
			call add_node(whichgL, TreeBranchNum)
			NB1 => current%child

			call change_node(NB1, divX, 4)
			call change_node(NB1, BYMax, 2)
			call change_node(NB1, BXMin, 3)
			call change_node(NB1, divY, 1)
			call change_node(NB1, divY, 5)
			call change_node(NB1, divY, 6)
			call change_node(NB1, abs( (BYMax-divY)/(divX-BXMin) ), 7)
			call change_node(NB1, real(NPts1,KIND=r8), 8)
			call change_node(TreeBranch, real(whichgL,KIND=r8), 10,j)

  ! updates the "WhichChildren" particle list
			call NewParticleList( tempPartList, Quad1(1:NPts1), TreeBranch%StartPosition, &
				NumPts, NPts1, NB1%Name)  
		endif
		if (NPts2 > 0) then
			j=j+1; whichgL=whichgL+1
			call add_node(whichgL, TreeBranchNum)
			NB2 => current%child

			call change_node(NB2, BXMax, 4)
			call change_node(NB2, BYMax, 2)
			call change_node(NB2, divX, 3)
			call change_node(NB2, divY, 1)
			call change_node(NB2, abs( (BYMax-divY)/(BXMax-divX) ), 7)
			call change_node(NB2, real(NPts2,KIND=r8), 8)
			call change_node(TreeBranch, real(whichgL,KIND=r8), 10,j)

  ! updates the "WhichChildren" particle list
			call NewParticleList( tempPartList, Quad2(1:NPts2), TreeBranch%StartPosition+NPts1, &
				NumPts, NPts2, NB2%Name)  
		endif
		if (NPts3 > 0) then
			j=j+1; whichgL=whichgL+1
			call add_node(whichgL, TreeBranchNum)
			NB3 => current%child

			call change_node(NB3, divX, 4)
			call change_node(NB3, divY, 2)
			call change_node(NB3, BXMin, 3)
			call change_node(NB3, BYMin, 1)
			call change_node(NB3, abs( (divY-BYMin)/(divX-BXMin) ), 7)
			call change_node(NB3, real(NPts3,KIND=r8), 8)
			call change_node(TreeBranch, real(whichgL,KIND=r8), 10,j)

  ! updates the "WhichChildren" particle list
			call NewParticleList( tempPartList, Quad3(1:NPts3), TreeBranch%StartPosition+NPts1+NPts2, &
				NumPts, NPts3, NB3%Name)  
		endif
		if (NPts4 > 0) then
			j=j+1; whichgL=whichgL+1
			call add_node(whichgL, TreeBranchNum)
			NB4 => current%child

			call change_node(NB4, BXMax, 4)
			call change_node(NB4, divY, 2)
			call change_node(NB4, divX, 3)
			call change_node(NB4, BYMin, 1)
			call change_node(NB4, abs( (divY-BYMin)/(BXMax-divX) ), 7)
			call change_node(NB4, real(NPts4,KIND=r8), 8)
			call change_node(TreeBranch, real(whichgL,KIND=r8), 10,j)

  ! updates the "WhichChildren" particle list
			call NewParticleList( tempPartList, Quad4(1:NPts4), TreeBranch%StartPosition+NPts1+NPts2+NPts3, &
				NumPts, NPts4, NB4%Name)  
		endif
		call change_node(TreeBranch, real(j,KIND=r8), 12)

		deallocate(Quad1); deallocate(Quad2)
		deallocate(Quad3); deallocate(Quad4)

	END SUBROUTINE SplitAndCount

! ***************
! Subroutine SortParticleList
! ***************
! DESCRIPTION: Sorts particles into <= 4 divisions in one master particle position listing.
!	Also assigns %WhichChildren the beginning and end particle _position_ integer
!
! INPUT: WhichChildren, ParticleList, tempList
!
! OUTPUT: Sorted tempList
!
! CALLING PROGRAM: TREECODE
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   15 Feb 2005 begun (v4)

	SUBROUTINE NewParticleList( tempList, KidsList, StartPos, NumPts, NumKids, TreeBranchNum)

! incoming variables
		integer, intent(in)						:: StartPos, NumPts, TreeBranchNum, NumKids
		integer, dimension(NumKids),intent(in)	:: KidsList

! outgoing variables
		integer,dimension(NumPts),intent(inout)	:: tempList

! subroutine only vars
		integer									:: i, i1
		real*8									:: End2
		type(node), pointer						:: TreeBranch

		continue

		call find2(TreeBranchNum); TreeBranch=>current
		i1=0
		do i=StartPos,StartPos+TreeBranch%NumPtInclude-1
			i1=i1+1
			tempList(i) = KidsList(i1)
		enddo

		End2=StartPos+TreeBranch%NumPtInclude-1

		call change_node(TreeBranch, real(StartPos,KIND=r8),9,1)
		call change_node(TreeBranch, End2,9,2)
		call change_node(TreeBranch, real(StartPos,KIND=r8),11)

	END SUBROUTINE NewParticleList

! ***************
! Subroutine SizeTrees
! ***************
! DESCRIPTION: Forms the pointer tree structure for dirichlet and neumann panels
!
! INPUT: size of trees, potential, x&y points
!
! OUTPUT: populated pointer tree 
!
! CALLING PROGRAM: TreeCode
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   01 Mar 2005 begun 

	SUBROUTINE SizeTrees (PartList, NewPartList, XPts2, YPts2, sizeNewTree, whichgL)
! incoming variables
		integer, intent(inout)						:: sizeNewTree, whichgL
		real*8, dimension(sizeNewTree), intent(in)	:: XPts2, YPts2

! outgoing variables
		integer, dimension(sizeNewTree), intent(inout) :: PartList, NewPartList

! subroutine only vars
		integer										:: i, LoopDivide
		type(node), pointer							:: TreePoint2

continue

! intializes and starts the tree
		call start()	
		call add_node(1,0)
		TreePoint2 => current%child
		do i=1,sizeNewTree ! create particle master list
			PartList(i)=i
		enddo
		call change_node(TreePoint2,real(sizeNewTree,KIND=r8),8)
		call change_node(TreePoint2,real(1.,KIND=r8),11)
		call NewParticleList( NewPartList, PartList, 1, sizeNewTree, sizeNewTree, 1)

!********
! construct tree
!********
		i=0; LoopDivide=1
		do while ( 1 == LoopDivide )
			i=i+1

			call find2(i)
			TreePoint2 => current

			if (TreePoint2%NumPtInclude > PtsPerBoxMax ) then
    			! find the size of the box
				call Resizer(  i, TreePoint2%NumPtInclude, sizeNewTree, XPts2, YPts2, NewPartList )

				call SplitAndCount ( NewPartList, TreePoint2, i, XPts2, YPts2, &
					TreePoint2%NumPtInclude, sizeNewTree, whichgL )
			elseif ( TreePoint2%NumPtInclude >= 1 ) then
				call Resizer(  i, TreePoint2%NumPtInclude, sizeNewTree, XPts2, YPts2, NewPartList )
			else ! no points in the box
				stop 'potential error in box'
			endif

  ! you are at the last child
			if (whichgL == i) then
				LoopDivide = 0
			endif
		enddo
	END SUBROUTINE SizeTrees

! ***************
! Subroutine FormTrees
! ***************
! DESCRIPTION: Takes pointers and makes regular array of tree and then erases pointers
!
! INPUT: size of trees, TreeStructure
!
! OUTPUT: populated TreeStructure 
!
! CALLING PROGRAM: TreeCode
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   06 Mar 2005 begun 

	SUBROUTINE FormTrees(NewTree, NumBranches)

! incoming variables
		integer, intent(in)							:: NumBranches

! outgoing variables
		type(node), dimension(NumBranches), intent(out)	:: NewTree

! subroutine only vars
		integer										:: i
		type(node), pointer 						:: TreePoint2

		continue

		do i=1,NumBranches
			call find2(i); TreePoint2 => current
			NewTree(i)= TreePoint2
			NewTree(i)%parent=TreePoint2%parent
		enddo
		call find2(0); call finish	! erase all pointers from above

	END SUBROUTINE FormTrees

! ***************
! Subroutine TaylorCoeff
! ***************
! DESCRIPTION: calculates ak(xi,yc) Taylor coefficiencts for ||k|| < p 
!	order approximation
!
! INPUT: Tree, p, pt
!
! OUTPUT: ak
!
! CALLING PROGRAM: TREECODE
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   16 Dec 2004 begun
!	03 Feb 2005 v2 Changing to a case structure so each component isn't the 
!				     sum of the previous components but a standalone
!	15 Feb 2005 v5 Changing to sum the entire thing in one loop and minimize computation
!					 with stand-in terms.
!	15 Mar 2005 v7 Removing vector influence coefficients vi. Should only be in moment calc
!					 Integrating Greens function along panels, not just a point. 
!	24 Mar 2005 v9 Removing Potential calc, outputting Taylor coeff so can compute only once

	SUBROUTINE TaylorCoeff ( Tlk, TreeBranch2, XPt, YPt, k, L  )

! incoming variables
		integer, intent(in)		:: k, L
		real*8, intent(in)		:: XPt, YPt					
		type(node), intent(in)	:: TreeBranch2

! outgoing variables
		real*8, intent(out)		:: Tlk

! subroutine only vars
		integer					:: casekL
		real*8					:: delX, delY, R, R2, delX2, delY2, TwopiR, TwopiR2
		continue

		delX = XPt - TreeBranch2%xcen; delY = YPt - TreeBranch2%ycen
		delX2 = delX*delX; delY2 = delY*delY;
		R = delX2 + delY2; R2 = R*R
		TwopiR = 2.*pi*R; TwopiR2 = TwopiR * R

		if (R /= 0) then
			casekL = 10*k+L
		else
			casekL = 1000
		endif

		select case(casekL)
		case (00)
			Tlk = log( R )/(4.*pi)
		case (10)
		!Ak = delY/(2.*pi*R)
			Tlk =  delY/(TwopiR)
		case (11)
		!Ak = delY/(2.*pi*R)
			Tlk =  delY/(TwopiR)
		case (20)
		!Ak = 0.5*( (R-2.*delY**2) )/(2.*pi*R2)
			Tlk =  0.5*( (R-2.*delY2) )/(TwopiR2)
		case (21)
		!Ak = 0.5*( -delX*delY )/(pi*R2)
			Tlk =  0.5*( -delX*delY )/(pi*R2)
		case (22)
		!Ak = 0.5*( (R-2.*delX**2) )/(2.*pi*R2)
			Tlk =  0.5*( (R-2.*delX2) )/(TwopiR2)
		case (30)
		!Ak = (8.*delY**3/R -6.*delY) /(12.*pi*R2) 
			Tlk =  (8.*delY2*delY/R -6.*delY) /(6*TwopiR2) 
		case (31)
		!Ak = (8.*delX*delY**2/R-2.*delX)/(4.*pi*R2)
			Tlk =  (8.*delX*delY2/R-2.*delX)/(2*TwopiR2)
		case (32)
		!Ak = (8.*delX**2*delY/R-2.*delY)/(4.*pi*R2)
			Tlk =  (8.*delX2*delY/R-2.*delY)/(2*TwopiR2)
		case (33)
		!Ak = (8.*delX**3/R -6.*delX) /(12.*pi*R2)
			Tlk =  (8.*delX2*delX/R -6.*delX) /(6*TwopiR2)
		case (40)
		!Ak = ( -48.*delY**4/R2 + 48.*delY**2/R - 6.)/(48.*pi*R2)
			Tlk =  ( -48.*delY2*delY2/R2 + 48.*delY2/R - 6.)/(24*TwopiR2)
		case (41)
		!Ak = ( 24.*delX*delY/R - 48.*delX*delY**3/R2 )/(12.*pi*R2) 
			Tlk =  ( 24.*delX*delY/R - 48.*delX*delY2*delY/R2 )/(6*TwopiR2)
		case (42)
		!Ak = ( -48.*(delX**2)*delY**2/R2 + 8.*delX**2/R + 8.*delY**2/R - 2 )/(4.*pi*R2)
			Tlk =  ( -48.*delX2*delY2/R2 + 8.*delX2/R + 8.*delY2/R - 2 )/(2*TwopiR2)
		case (43)
		!Ak = ( 24.*delX*delY/R - 48.*(delX**3)*delY/R2 )/(12.*pi*R2) 
			Tlk =  ( 24.*delX*delY/R - 48.*delX2*delX*delY/R2 )/(6*TwopiR2) 
		case (44)
		!Ak = (-48.*delX**4/R2 + 48*delX**2/R - 6)/(48.*pi*R2)
			Tlk =  (-48.*delX2*delX2/R2 + 48*delX2/R - 6)/(24*TwopiR2)
		case (1000)
		! do nothing, self point
		case default
			stop 'Err casekL'
		end select 

	END SUBROUTINE TaylorCoeff

! ***************
! Subroutine ClusterMoment
! ***************
! DESCRIPTION: calculates M(l,k), the moments around the cluster
!
! INPUT: Tree, l, k, coordinates
!
! OUTPUT: M(l,k)
!
! CALLING PROGRAM: TREECODE
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   15 Dec 2004 begun
!	17 Feb 2005 v5. Changing to do all calculations in here at once to reduce function calls

	SUBROUTINE ClusterMoment( M_lk4, XPts, YPts, PanGaussWeight, NumPts, SortedPositions2, &
			TreeBranch, whichgL, GaussPts, PanLen, vi )

		! incoming variables
		integer, intent(in)							:: NumPts, whichgL, GaussPts
		integer, dimension(NumPts), intent(in)		:: SortedPositions2
		real*8, dimension(NumPts), intent(in)		:: vi, PanLen
		real*8, dimension(GaussPts), intent(in)		:: PanGaussWeight
		real*8, dimension(NumPts*GaussPts), intent(in)	:: XPts, YPts
		type(node), dimension(whichgL), intent(in)	:: TreeBranch

		! outgoing variables
		real*8, dimension(whichgL, Pfac+1, Pfac+1), intent(out) :: M_lk4

		! subroutine only vars
		integer				:: i, k, L, m, o, i1, SizeBranch, startBranch, panNum
		real*8				:: xMom, yMom, Xpos, Ypos

		continue
		TreeDo: do m=1,whichgL	! loop over all Trees
			startBranch = TreeBranch(m)%StartPosition
			SizeBranch = TreeBranch(m)%NumPtInclude
			Xpos = TreeBranch(m)%xcen; Ypos = TreeBranch(m)%ycen
			PfacDo: do k=0,Pfac ! up to level of approximation Pfac
				do L=0,k  ! ..in 2 dimensions
					M_lk4(m,k+1,L+1)=0

					NumPtsDo: do i=1,SizeBranch
						i1=SortedPositions2(startBranch+i-1)
						ODo: do o=1,GaussPts
							panNum=(i1-1)*GaussPts+o
							xMom = ( XPts(panNum) - Xpos )**L 
							yMom = ( YPts(panNum) - Ypos )**(k - L)
							M_lk4(m,k+1,L+1) = M_lk4(m,k+1,L+1) - 0.5*PanLen(i1)*vi(i1) * xMom * yMom * PanGaussWeight(o)
						enddo ODo
					enddo NumPtsDo
				enddo
			enddo PfacDo
		enddo TreeDo

	END SUBROUTINE ClusterMoment

! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd
! http://www.fortran.com/fortran/qsort_c.f95

	recursive subroutine QsortC(A)
		real, intent(in out), dimension(:) :: A
		integer :: iq

		if(size(A) > 1) then
			call Partition(A, iq)
			call QsortC(A(:iq-1))
			call QsortC(A(iq:))
		endif
	end subroutine QsortC

	subroutine Partition(A, marker)
		real, intent(in out), dimension(:) :: A
		integer, intent(out) :: marker
		integer :: i, j
		real :: temp
		real :: x      ! pivot point
		x = A(1)
		i= 0
		j= size(A) + 1

		do
			j = j-1
			do
				if (A(j) <= x) exit
				j = j-1
			end do
			i = i+1
			do
				if (A(i) >= x) exit
				i = i+1
			end do
			if (i < j) then
        ! exchange A(i) and A(j)
				temp = A(i)
				A(i) = A(j)
				A(j) = temp
			elseif (i == j) then
				marker = i+1
				return
			else
				marker = i
				return
			endif
		end do

	end subroutine Partition
END MODULE TreeOperations
