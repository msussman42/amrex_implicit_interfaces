MODULE TreeCode
	use TreeOperations, only	: SizeTrees, FormTrees
	use TreeCodeGlobal

	implicit none

! used to restrict variables to prevent accidental calling outside
	private      ! everything defaults to 'private' except those expressly labeled otherwise.
	public		:: TreeCode_Matrix

contains 

! ***************
! Subroutine TreeCode_Matrix
! ***************
! DESCRIPTION: builds fluid and free space trees
!
! INPUT: xarray, yarray, sizeArray
!
! OUTPUT: built trees x 2
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   09 Dec 2004 begun
!	30 Jan 2005 v2. Precompiling moments and direct recursion potentials
!	07 Feb 2005 v3. Changing to subroutine to use with FSU code. Allocating
!					  %WhichPtInclude on a tree by tree basis. Removing tempTree
!	15 Feb 2005 v4. Current problem: extreme time calling 'ComputePotential'. Changing to 
!					  position allocation pointers to store smaller arrays. Increasing max size
!					  of trees to 50 members from 1.
!				v5. (v4 working). Changing to a listing of WhichChildren position beginning:end
!					  instead of listing all of them.
!   17 Feb 2005 v6. (v5 working, but very slow). Changing every array away from pointers for a speedup 
!	21 Feb 2005 v7. (v6 working!) Changing to be called from MAT_Greens
!	25 Feb 2005 v8. v7 not working. Changing to variable number of charges qi so can work on 
!					  force calculations
!	27 Feb 2005 v9. v8 not working. Removing variable qi, setting up to run only from MatGreens. Defining
!					  3 trees - surface, electrode, neumann. Not calculating potentials or forces, just 
!					  tree structure itself
!	07 Mar 2005 v10 v9 not working. Integrating TreeCode into a subfunction of GmRes. Passing out tree
!					  structure
!	15 Mar 2005 v11 v10 working. Changing back to only 2 cases, D and N

	SUBROUTINE TreeCode_Matrix (TreeGMRes, indxBounds, XBoundaryIn, YBoundaryIn, BoundType, NumPts)

	! incoming variables
		integer, intent(in)						:: NumPts
		integer,dimension(NumPts),intent(in)	:: BoundType
		real*8, dimension(NumPts), intent(in)	:: XBoundaryIn, YBoundaryIn

	! outgoing variables
		integer,dimension(NumPts,2), intent(out):: indxBounds
		type(TreePackage), dimension(2), intent(out) :: TreeGMRes

	! local variables
		integer									:: i, whichgL
		real*8									:: t1, t2, timeCount

		continue
	!-----------
		call cpu_time(t1)

	! There are neumann and dirichlet panels being passed to this subroutine (as evidenced by the
	!   values in 'BoundType'). Divide the XBoundary, YBoundary into D-Surface(1)  or Neumann(2) components
		TreeGMRes(1)%NumPtsTree=0; allocate(TreeGMRes(1)%XPtsTree(NumPts)); 
		allocate(TreeGMRes(1)%YPtsTree(NumPts))
		TreeGMRes(2)%NumPtsTree=0; allocate(TreeGMRes(2)%XPtsTree(NumPts)); 
		allocate(TreeGMRes(2)%YPtsTree(NumPts))

		do i=1,NumPts
			if ( 1==BoundType(i)) then	! it is a fixed-potential surface
				whichgL=1
			elseif ( 0==BoundType(i) ) then		! it is a fixed-flux
				whichgL=2
			else									! error
				print *, 'TreeCode11 divide PtsIn error'
				stop 
			endif
			TreeGMRes(whichgL)%NumPtsTree=TreeGMRes(whichgL)%NumPtsTree+1
			TreeGMRes(whichgL)%XPtsTree(TreeGMRes(whichgL)%NumPtsTree)=XBoundaryIn(i)
			TreeGMRes(whichgL)%YPtsTree(TreeGMRes(whichgL)%NumPtsTree)=YBoundaryIn(i)
			indxBounds(i,1)=whichgL; indxBounds(i,2)=TreeGMRes(whichgL)%NumPtsTree
		enddo

		! form both trees 
		allocate(TreeGMRes(1)%PtsTreeArray(TreeGMRes(1)%NumPtsTree))
		allocate(TreeGMRes(1)%NewPtsTreeArray(TreeGMRes(1)%NumPtsTree))
		allocate(TreeGMRes(2)%PtsTreeArray(TreeGMRes(2)%NumPtsTree))
		allocate(TreeGMRes(2)%NewPtsTreeArray(TreeGMRes(2)%NumPtsTree))

		!call SizeTrees (PartDElec, NewPartDElec, XBoundDElec, YBoundDElec, NumPtsDElec, glCountDElec)
		!allocate(TreeDElec(glCountDElec)); call FormTrees (TreeDElec, glCountDElec)
		if (TreeGMRes(1)%NumPtsTree >0) then
			call SizeTrees (TreeGMRes(1)%PtsTreeArray, TreeGMRes(1)%NewPtsTreeArray, &
				TreeGMRes(1)%XPtsTree, TreeGMRes(1)%YPtsTree, TreeGMRes(1)%NumPtsTree, glCountD)
			allocate(TreeGMRes(1)%TreeStructure(glCountD)); call FormTrees (TreeGMRes(1)%TreeStructure, glCountD)
		else
			glCountD=glCountD-1
		endif

		if (TreeGMRes(2)%NumPtsTree >0) then
			call SizeTrees (TreeGMRes(2)%PtsTreeArray, TreeGMRes(2)%NewPtsTreeArray, &
				TreeGMRes(2)%XPtsTree, TreeGMRes(2)%YPtsTree, TreeGMRes(2)%NumPtsTree, glCountNeu)
			allocate(TreeGMRes(2)%TreeStructure(glCountNeu)); call FormTrees (TreeGMRes(2)%TreeStructure, glCountNeu)
		else
			glCountNeu=glCountNeu-1
		endif

		call cpu_time(t2); timeCount = t2-t1
		write(*, fmt=13) '     BuildTree took ', timeCount, ' seconds. Theta= ', theta

		goto 999		! actual program end, just lists formats and errors below
! ***********

! FORMAT listings
		13 FORMAT(a,f8.3,a,f4.2)

		999 print *, '  Done with program <TREECODE>'

	end SUBROUTINE TreeCode_Matrix
end MODULE TreeCode
