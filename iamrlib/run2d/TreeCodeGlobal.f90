! ***************
! Module TreeCodeGlobal[x]
! ***************
! DESCRIPTION: Global data, parameters, and structures  
!
! INPUT: <na>
!
! OUTPUT: Potential field at a given set of points
!
! Calling Program: TREECODE.f90
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   09 Dec 2004 begun
!   07 Feb 2005 v2. changing to use CallTree as main program
!	08 Feb		v3. removing data structures, leaving constants
!   01 Mar		v4. (v3 fine) Adding global integers to create 3 different trees
!	15 Mar		v5. (v4 fine) Reducing to 2 trees, N&D

module TreeCodeGlobal
	use directoryTree, only : node

	implicit none

!input user-specificed parameters
	integer, parameter, public	:: PtsPerBoxMax = 10, Pfac=4, &
		r8=SELECTED_REAL_KIND(14,16) ! type that has 14 digits of precision and e-16 coverage
	integer, SAVE, public 		:: counter

	real(kind=r8), parameter, public	:: pi=3.14159, & 
	epsil_0=8.8542e-12, &			! permittivity of free space, F/m
	q = 1.6022e-19, &				! charge of one electron, C
	q_epsil = 1.80954e-8, &			! q/epsil_0
	theta=0.6, &					! distance acceptance parameter
	epsil=1e-6, smallepsil=3e-15, toolarge=1e8

	real(kind=r8), dimension(:), allocatable :: sigma

!  Global data
	integer						:: glCountD=1, glCountNeu=1

! Define the basic tree package for GMRes
	type TreePackage
		integer								:: NumPtsTree
		integer, dimension(:), pointer  	:: PtsTreeArray, NewPtsTreeArray	
		real*8, dimension(:), pointer		:: XPtsTree, YPtsTree
		type(node), dimension(:), pointer	:: TreeStructure
	end type TreePackage

end module TreeCodeGlobal
