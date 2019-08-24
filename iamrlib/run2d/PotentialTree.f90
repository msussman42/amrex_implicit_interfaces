MODULE PotentialTree
! 08 Feb 2005 v3 
! 15 Feb 2005 v4
! 17 Feb 2005 v5
! 01 Mar 2005 v6
! 15 Mar 2005 v7
! 24 Mar 2005 v8

	use directoryTree, only		: node
	use TreeOperations, only	: TaylorCoeff
	use SortNodePts,  only		: ExactPanelInteg, TrapInit
	use TreeCodeGlobal

	implicit none

! used to restrict variables to prevent accidental calling outside
	private      ! everything defaults to 'private' except those expressly labeled otherwise.
	public		:: ComputeAn, Greens, ExactIntegrateRoutine, CalcPot

	real*8, parameter	:: tooStretched = 1.414213 ! 2^0.5

contains 

! ***************
! Subroutine ComputeAn
! ***************
! DESCRIPTION: determines the largest acceptable cluster and calculates the 
!   potential along it
!
! INPUT: point, tree
!
! OUTPUT: potential phi
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
!   20 Dec 2004 begun
!   01 Feb 2005 v2. 1st shot at Taylor coefficients
!	08 Feb 2005 v3. Changing to pointers and calling from CallTree
!   15 Feb 2005 v4. Speeding up, changing to 2 %WhichPts
!	17 Feb 2005 v5. v4 working, changing to all array pointers to try and increase speed
!	01 Mar 2005 v6. v5 working, changing to multiple tree calls
!	15 Mar 2005 v7. v6 working, changing to integrate Greens function
!
!function ComputePotential(xBar; C) algorithm
! 1)is the Taylor approximation sufficiently accurate?
! 2) if yes, compute phi(xBar; C) by Taylor approximation
! 3) if no, does C have any sub-clusters?
! 4)   if yes, call ComputePotential(xBar; C') for each sub-cluster C' of C
! 5)   if no, compute phi(xBar; C) by direct summation

	recursive SUBROUTINE ComputeAn( XPts, YPts, iPt, Potential, TreeBranch, sizeTree, &
			NumPts, SortedParticleList, Mlk3, WhichPt, normV, RHSType2, normPanX, normPanY, &
			PanGaussWeight, PanLen, M, xPanCen, yPanCen, xPan, yPan, Tlk3 )

		! incoming variables
		integer, intent(in)						:: WhichPt, sizeTree, NumPts, iPt, M
		integer, dimension(NumPts), intent(in)	:: SortedParticleList 
		real*8, dimension(NumPts), intent(in)	:: XPts, YPts, normV, RHSType2, normPanX, &
			normPanY, PanLen, xPanCen, yPanCen
		real*8, dimension(M), intent(in)		:: PanGaussWeight 
		real*8, dimension(NumPts*M), intent(in) :: xPan, yPan
		real*8, dimension(sizeTree, Pfac+1,Pfac+1),intent(inout) :: Mlk3
		real*8, dimension(NumPts, sizeTree, Pfac+1,Pfac+1),intent(inout) :: Tlk3
		type(node),dimension(sizeTree),intent(in) :: TreeBranch

		! outgoing variables
		real*8, intent(inout)	:: Potential

		! subroutine only vars
		integer					:: i, j, k, L, i1, NoChildren, startBranch
		real*8					:: s, r, tempPot, delX, delY, testTlk

		continue
		if (associated(TreeBranch(WhichPt)%WhichChildren)) then
			NoChildren = size(TreeBranch(WhichPt)%WhichChildren)
		else
			NoChildren = 0
		endif
		delX = XPts(iPt)-TreeBranch(WhichPt)%xcen; delY = YPts(iPt)-TreeBranch(WhichPt)%ycen

		s = max((TreeBranch(WhichPt)%xmax-TreeBranch(WhichPt)%xmin)/2, &
			(TreeBranch(WhichPt)%ymax-TreeBranch(WhichPt)%ymin)/2)
		r = sqrt( delX*delX + delY*delY )

		! compares s/r to theta
		if ( 0==NoChildren ) then	! you are at an ending node (5)
			startBranch = TreeBranch(WhichPt)%StartPosition	! direct sum the points to get aN
			do j=1,TreeBranch(WhichPt)%NumPtInclude
				i1=SortedParticleList(startBranch+j-1) 
				delX = XPts(i1)-XPts(iPt); delY = YPts(i1)-YPts(iPt)

    			! attempt to Gaussian quadrature integrate panels from iPt panel centers.
				tempPot=GaussIntegrateRoutine( iPt, i1, RHSType2, xPanCen, yPanCen, xPan, yPan, normPanX, &
					normPanY, PanGaussWeight, PanLen, M, NumPts )

				Potential = Potential + tempPot*normV(i1)
			enddo
		elseif ( 1.*abs(s/r)<theta .AND. r/=0 ) then		! (2)
			! calculates the Taylor coefficient from each point to each cluster
			testTlk=Tlk3(iPt, WhichPt, 1, 1)
			do k=0,Pfac ! up to level of approximation Pfac
				do L=0,k  ! ..in 2 dimensions
					if (0==testTlk) then
						call TaylorCoeff ( Tlk3(iPt, WhichPt, k+1, L+1), TreeBranch(WhichPt), &
							XPanCen(iPt), YPanCen(iPt), k, L  )
					endif
					Potential=Potential-Tlk3(iPt,WhichPt,k+1,L+1)*Mlk3(WhichPt,k+1,L+1)
				enddo
			enddo
		elseif (1.*abs(s/r)>=theta) then !(4)
			do i=1,NoChildren
				call ComputeAn( XPts, YPts, iPt, Potential, TreeBranch, sizeTree, NumPts, SortedParticleList, &
					Mlk3, TreeBranch(WhichPt)%WhichChildren(i), normV, RHSType2, normPanX, normPanY, &
					PanGaussWeight, PanLen, M, xPanCen, yPanCen, xPan, yPan, Tlk3 )
			enddo
		else
			print *, 'Error ComputePotential'
			stop 
		endif

		goto 999

	999 END SUBROUTINE ComputeAn

!---------------------------
	FUNCTION Greens( zzz, aaa )
	! Greens subfunction

		real*8, intent(in)	:: zzz, aaa
		real*8				:: Greens

		continue

		Greens = -log( zzz*zzz + aaa*aaa ) / (4.*pi)
	end FUNCTION Greens

!---------------------------
	FUNCTION ExactIntegrateRoutine( xPanBegin, yPanBegin, xPanCen, yPanCen, n, FromPt, ToPt)
	! routine that integrates the panels from xCen,yCen through xPanBegin,yPanBegin with tan^-1(y/x)
		! incoming
		integer, intent(in)	:: n, FromPt, ToPt
		real*8, dimension(n), intent(in)	:: xPanBegin, yPanBegin, xPanCen, yPanCen

		! outgoing
		real*8				::  ExactIntegrateRoutine

		! internal
		integer				:: NextPan, AtanYes
		real*8				:: intBeg, intEnd, a, ei, fi, gi, sgnGi

		continue

		! FromPt==j, ToPt==i
		if (ToPt /=n) then
			NextPan=ToPt+1
		else
			NextPan=1
		endif

		! directly solve the Greens integral. dx component
		intBeg=xPanCen(FromPt)-xPanBegin(ToPt)
		intEnd=xPanCen(FromPt)-xPanBegin(NextPan)
		a=yPanCen(FromPt)-yPanCen(ToPt)
		if (0==a) then
			AtanYes=0
		else
			AtanYes=1
		endif
		call ExactPanelInteg(ei, intBeg,intEnd,a,AtanYes)

		! directly solve the Greens integral. dy component
		intBeg=yPanCen(FromPt)-yPanBegin(ToPt)
		intEnd=yPanCen(FromPt)-yPanBegin(NextPan)
		a=xPanCen(FromPt)-xPanCen(ToPt)
		if (0==a) then
			AtanYes=0
		else
			AtanYes=1
		endif
		call ExactPanelInteg(fi, intBeg,intEnd,a,AtanYes)

		gi=sqrt( fi*fi + ei*ei)
		if (abs(fi)>abs(ei)) then
			sgnGi=sign(1.0_r8, fi)
		elseif (abs(fi)<abs(ei)) then
			sgnGi=sign(1.0_r8, ei)
		else
			sgnGi=1
		endif
		ExactIntegrateRoutine = sgnGi*gi	  

	end FUNCTION ExactIntegrateRoutine

!---------------------------
	FUNCTION GaussIntegrateRoutine( FromPt, ToPt, RHSTypes, xPanCen, yPanCen, xPan, yPan, &
			normPanX, normPanY, PanGaussWeight, PanLen, M, n )
		! routine that integrates the panels from xCen,yCen through xPanBegin,yPanBegin with gaussian quadrature

		! incoming
		integer, intent(in)	:: FromPt, ToPt, M, n
		real*8, dimension(n), intent(in)	:: RHSTypes, normPanX, normPanY, PanLen, xPanCen, yPanCen
		real*8, dimension(M), intent(in)	:: PanGaussWeight
		real*8, dimension(n*M), intent(in)	:: xPan, yPan

		! outgoing
		real*8				::  GaussIntegrateRoutine

		! internal
		integer				:: k, PanFromTo, panel
		real*8				:: ai, bi, ci, di, delX, delY, GreenPart

		continue
		! j==FromPt, i==ToPt
		! normPans(:,1) = normPanX; normPans(:,2) = normPanY 
		PanFromTo = int( RHSTypes(FromPt)*10 + RHSTypes(ToPt) )
		ai=0; bi=0; ci=0; di=0

		do k = 1,M
			panel=(ToPt-1)*M + k
			delX=xPanCen(FromPt)-xPan(panel)
			delY=yPanCen(FromPt)-yPan(panel)

			select case (PanFromTo)		! 1=dirichlet, 0=neumann type BC
			case (11)
				GreenPart=0
				call TrapInit ( GreenPart, delX, delY, normPanX(ToPt), &
					normPanY(ToPt), normPanX(FromPt), normPanY(FromPt), PanGaussWeight(k), &
					PanLen(ToPt), 10, 0, k, M )

				ai = ai + GreenPart

			case (10)
				GreenPart=0
				call TrapInit ( GreenPart, delX, delY, normPanX(ToPt), &
					normPanY(ToPt), normPanX(FromPt), normPanY(FromPt), PanGaussWeight(k), &
					PanLen(ToPt), 11, 0, k, M )

				bi = bi - GreenPart 		  

			case (01)
				GreenPart=0
				call TrapInit ( GreenPart, delX, delY, normPanX(ToPt), &
					normPanY(ToPt), normPanX(FromPt), normPanY(FromPt), PanGaussWeight(k), &
					PanLen(ToPt), 10, 0, k, M )

				ci = ci + GreenPart 

			case (00)
				GreenPart=0

				if (ToPt /= FromPt) then
					call TrapInit ( GreenPart, delX, delY, normPanX(ToPt), &
						normPanY(ToPt), normPanX(FromPt), normPanY(FromPt), PanGaussWeight(k), &
						PanLen(ToPt), 00, 0, k, M )

					di = di + GreenPart
				elseif ( M==k ) then
					di = 0.5
				endif 
			case default
				stop 'Error on PanFromTo aLine'
			end select	! switch PanFromTo
		end do			! for k
		GaussIntegrateRoutine = ai + bi + ci + di  

	end FUNCTION GaussIntegrateRoutine

	!---------------------------
	! calculates whether something is inside or outside a shape
	SUBROUTINE CalcPot(InsideYes, meshi, meshj, NumGrid, NumPan, PanX, PanY)
		! incoming
		INTEGER, INTENT(in)	:: NumGrid, NumPan
		REAL(kind=r8), DIMENSION(NumPan),INTENT(in) :: PanX, PanY

		! outgoing
		INTEGER, DIMENSION(NumGrid,NumGrid),INTENT(out) :: InsideYes 
		REAL(kind=r8), DIMENSION(NumGrid), INTENT(out) 	:: meshi, meshj

		! internal
		INTEGER				:: z, CrossCount, CL, CrossBeforePx, CrossAfterPx, CrossBeforePy, &
			CrossAfterPy, countInt, countIntersecOut, sizeCross, NoDouble, xx, yy, zz, CrossYes, &
			countIntersecIn, countDeadXY
		REAL(kind=r8)		:: Ax, Bx, Ay, By, Px, Py, TriangleArea, IntersecX, IntersecY
		REAL(kind=r8), DIMENSION(4) 		:: XYPos ! xmin, xmax, ymin, ymax
		REAL(kind=r8), DIMENSION(NumPan,3) 	:: CrossLocation, NewCrossLocation


		CONTINUE
		XYPos(1)=0; XYPos(2)=1.5; XYPos(3)=0.; XYPos(4)=1.5

		XXLoop: do xx = 1,NumGrid
			meshi(xx) = 0.1*(XYPos(2)-XYPos(1))+(xx-1)*0.8*(XYPos(2)-XYPos(1))/(1.*NumGrid-1.)
			YYLoop: do yy = 1,NumGrid
				IF (1==xx) THEN
					meshj(yy) = 0.1*(XYPos(4)-XYPos(3))+(yy-1)*0.8*(XYPos(4)-XYPos(3))/(1.*NumGrid-1.)
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
			end do YYLoop		! for yy
		end do XXLoop           ! for xx
	END SUBROUTINE CalcPot

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

	SUBROUTINE InSideYesNo( CrossCount, CrossYes, IntersecX, IntersecY, CrossDir, &
			Ax, Ay, Bx, By, Px, Py )

	! incoming variables
		real(kind=r8), intent(in)	:: Ax, Ay, Bx, By, Px, Py
		character, intent(in)	:: CrossDir

 	! outgoing variables
		INTEGER, intent(out)	:: CrossYes
		integer, intent(inout)	:: CrossCount
		real(kind=r8), intent(out)	:: IntersecX, IntersecY

	! subroutine entirely variables
		real(kind=r8)			:: delBAx, delBAy, t

		continue

		CrossYes=0
		delBAx = Bx-Ax; delBAy = By-Ay

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
			CrossCount = CrossCount+1; CrossYes=1
		endif
	end SUBROUTINE InSideYesNo
END MODULE PotentialTree
