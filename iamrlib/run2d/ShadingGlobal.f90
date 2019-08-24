! ***************
! Module ShadingGlobal[x]
! ***************
! DESCRIPTION: Global data, parameters, and structures  
!
! Calling Program: Shading.f90
!
! Primary author: Anton VanderWyst
!	University of Michigan
!	1320 Beal Ave.
!	Ann Arbor, MI 48109
!	antonv@umich.edu 
!
! Version history:
!   24 Jul 2005 begun

MODULE ShadingGlobal
	IMPLICIT NONE

    !input user-specificed parameters
	LOGICAL, PARAMETER, PUBLIC	:: UsePerturb=.FALSE. ! do you perturb the initial surface with a sine wave?
	INTEGER, PUBLIC, SAVE		:: counter, highblobNum

	INTEGER, PARAMETER, PUBLIC	:: r8=SELECTED_REAL_KIND(14,19), &
    	! type that has 14 digits of precision and e-19 coverage
	i2=SELECTED_INT_KIND(2), & ! type that has covers -99:99. Use for *SHORT* integers
	StartNum=1,StopNum=0, & ! modifies reading the LS data to avoid spurious edge pts. 1,0=non-truncated
	flagNN=6, &				! =0, sort by increasing theta. =1, sort by nearest neighbor (NN)
		! =2, sort by NN+norm limiting, =3 sort by NN+no cross, =4 sort by NN+no long, =5 sort by weighted
		! angle/distance, =6 no sort (from 'surfaceTracker')
	InitArySize=4000, &		! how big is your initial array size guess
	SmallInitAry=100, &		! how big is your initial small array size guess
	M=8, &					! # pts / panel; # panels / shortest side
	PotMeshSize=120, &		! number of grid cells^2 for potential evaluation
	CubicSplineFit=0, & 	! if '0', linear spline fit. if '1', cubic spline fit.	
	MatSolveMethod=3, &		! how to solve Ax=b matrix? =1 GMRes; =2 Gaussian inversion w/ pp; =3 Gauss w/o pp
	MinNumPtsBlob=10, &		! what is the minimum number of points per blob?
	ShowPassedPhi=0, &		! flag for printing the passed phi array
	ShowNewCross=0, ShowSort=0, & ! flags to print data files about surface and sorted pts
	ShowPt_Force=1, &		! flag for printing out the [evenly spaced coords, 10] large data file
	ShowAngleSort=0, &		! flag for printing unsorted grouped blobs
	ShowSubPan=0, &			! flag for printing New Coordinates panel details
	ShowCoords=1, &			! flag for printing Coordinates
	ShowNewCoords=1, &		! flag for printing New Coordinates
	ShowPotentials=0, &		! if=1, run extra routine to calc. phi on arb. grid inside needle
	LoadPastStep=0			! temp fix to load last phi step and generate surface

	! physical values
	REAL(KIND=r8), PARAMETER, PUBLIC	:: pi=3.14159, & 
	epsil0=8.8542e-12, &	! permittivity of free space, F/m
	q = 1.6022e-19, &		! charge of one electron, C
	q_epsil = 1.80954e-8, &	! q/epsil_0
	temp_in = 450, &		! experiment temperature [K]
	surften_in = 0.555-(0.12/1000)*(temp_in-430), &		! indium's surface tension at temp_in [N/m]
	rho_in = 1000*(7.1295-0.67987e-3*(temp_in-273.15)), & ! indium's density [kg/m^3]
	epsil=1e-6, smallepsil=3e-15, toolarge=1e8

	! simulation values
	REAL(KIND=r8), PARAMETER, PUBLIC	:: tol=1e-7, & 	! relative tolerance for GMRes
	elecPotential=1e10, &	! electrode potential [V]
	fluidPot=-2.0, &		! base fluid potential [V]
	maxDistStep=0.05, &	! max dist. between evenly spaced points [MAJOR driver of total compute time]  originally 0.020
	StepMult=10.0, &		! multiply the number of interior panels by [* major time driver *]  originally 23.0
	toofarNorm=1.6, &		! how close do norms have to be to connect (in rad)
	maxPerturb=0.05, &		! what is the max perturbation value for the surface?
	highPDFpt=0.78			! at what height do you take the second PDF setup?

	! global arrays to pass between subroutines for old and new pts
	REAL(KIND=r8), DIMENSION(:), ALLOCATABLE, PUBLIC 	:: sigma
	REAL(KIND=r8), DIMENSION(InitArySize,8), PUBLIC		:: NewCoords
	  ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatpan]
	REAL(KIND=r8), DIMENSION(SmallInitAry,3), PUBLIC	:: ChokeCoord 	! lists choke points and how many times called

	INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC		:: alreadyChecked, insideMarker ! for 'surfacetrack' levelsets...

	REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE, PUBLIC	:: phi			! level set values and tying Coords->NewCoords
END MODULE ShadingGlobal
