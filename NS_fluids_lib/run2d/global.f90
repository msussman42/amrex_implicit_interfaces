
MODULE global
	! includes all constants and type declarations for the DiagMatrix
	IMPLICIT NONE
	PUBLIC
    !input user-specificed parameters
	!--------------------------------
	INTEGER, PUBLIC, SAVE		:: counter ! counts how many consecutive times BFC called

	! print debugging files/values
	!--------------------------------
	INTEGER, PARAMETER, PUBLIC	:: CalcPotentials=1, &	! if=1, calc. phi on arb. grid inside needle
	ShowGrid=0, ShowNewCross=0, ShowSort=0, & ! if=1, prints levelset, surface and sorted intermed pts
	ShowPt_Force=0, &		! if=1, prints out the [evenly spaced coords, *] large PForce* data file
	ShowNewCoords=0, &		! if=1, prints New Coordinates
	ShowCoords=1, &			! if=1, prints Coordinates
	ShowPassedPhi=0, &		! if=1, prints passed phi array
	ShowAngleSort=0, &		! if=1, prints unsorted grouped blobs
	ShowDropLink=0, &		! if=1, prints drop center distances
	ShowChiAngle=0, &		! if=1, calcs & prints surface angle every 'SkipSteps'
	ShowAMat=0, &			! if=1, prints A matrix
	ShowSubPan=0, &			! if=1, prints all GaussPan panel details
	ShowGMResConverge=0, &	! if=1, prints how GMRes converges over time. Needs final answer from 'sigma.dat'
	ShowEachSurfCharge=0, & ! if=1, prints out a different SurfCharge.dat for each timestep
	ShowEachGaussPan=1, &	! if=1, prints out AllSurf[xxxxx].dat for each timestep
	ShowBlobArea=0, &		! if=1, prints out a different BlobCen.dat for each timestep	
	ShowInLoopCoords=0, &	! if=1, prints out multiple copies of Coords, NewCoords within loop
	LoadPastStep=0			! if=1, loads last phi step and generates surface

	! simulation setup: adjustable
	!--------------------------------
	INTEGER, PARAMETER, PUBLIC	:: MatSolveMethodGoal=1, &		! how to solve Ax=b matrix? =1 GMRes; 
	    ! =2 Gaussian inversion w/ pp; =3 using cpp full inversion
	CaseShape=3				!Outer shape. =1, single box, simple pot; =2 ripple box, real pot; 
							!=3, Suvorov test case with electrode on axis; =4 is a 'T' top shape

	REAL*8, PARAMETER, PUBLIC	:: elecPotential=1e10, &	! electrode potential [V]
	fluidPot=-2., &			! base fluid potential [V]
	maxDistStep=7e-4, &	! max dist. between evenly spaced points [* major time driver *]
	StepMult=0.05			! multiply number of baseline exterior panels by [* major time driver *]

	! simulation setup: can be, but not meant to be adjusted
	!--------------------------------
	INTEGER, PARAMETER, PUBLIC	:: AformDouble=0, & ! if=1, use D/N A matrix setup. Otherwise, use sigma,mu
	CalcPotBack=0, &		! if=1, calcs base phi for elec. field background
	ConnectDrops=1, &		! if=1, find nearest neighbor for detached drops and link. Force Dirch. BC
	CubicSplineFit=0, & 	! if=0, linear spline fit. if '1', cubic spline fit.
	ForceEdgeStable=0, &	! if=1, forces the edge to be stable after x>ForceZeroR
	HardwireChargeYes=0, &	! if=1, hardwire a charge and see the shape evolve.
		! REQUIRES SKIPPING OVER BEM approach
	InitArySize=5000, &		! max # of points/blob. size of outgoing radblob var. 
	EqualLsBemNumCoord=0, &	! if=1, have the number of BEM coordinates around the entire shape = number
		! of LS points on surface. Keeps accuracy ~ even
	MaxNumInpotIter=1, &	! maximum number of 'InnerPotSet' iteration calls
	MinNumPtsBlob=8 , &		! what is the minimum number of points per blob?
	npoints=10, & ! # pts / panel; # panels / shortest side
		! NOTE: need to make sure 'constants.h:npoints' has the same value as npoints
	PotMeshSize=100, &		! number of grid cells^2 for optional potential evaluation
	PotRedoStep=100, &		! number of timesteps between recalculating the background E field
	RecastDist=1, &			! if=1, expects input data to be in cm
	SkipFar=0, &			! if=1, in elliptics:calcpot, ignore points 'IgnoreImpact' far away
	SkipSteps=25, &			! number of timesteps between printing .dat files
	SmallInitAry=50, &		! how big is your initial small array size guess
	StartNum=1,StopNum=-1, & ! modifies reading the LS data to avoid 
	    ! spurious edge pts. 1,0=non-truncated
	UseGaussQuad=1, &		! if=1, use Gaussian quadrature to space points within panel. if=0, use linear
	UsePreCond=1			! if=1, use Jacobi preconditioning (invert diagonal)

	REAL*8, PARAMETER, PUBLIC	:: tol=1e-7, & 	! relative tolerance for GMRes
	highPDFpt=0.33, &		! at what height do you take the second PDF setup?
	toofarNorm=1.6, &		! how close do norms have to be to connect (in rad)
	HardwireCharge=7.11902e-5, & ! hardwired charge [C] to give 1e9 E field at 8mm (7.1e-6)
	toolarge=1e8, smallepsil=3e-15, eps=1e-10

	! physical values
	!--------------------------------
	REAL*8, PARAMETER, PUBLIC	:: pi=3.14159265359, & 
	epsil0=8.8542e-12, &	! permittivity of free space, F/m
	q = 1.602178e-19, &		! charge of one electron, C
	q_epsil = 1.80954e-8, &	! q/epsil_0
	chargeConstant=q/(2.0*pi*pi*epsil0), &				! changing charges around theta 2d
	cmtom = 0.01, &			! factor to change centimeters to meters
	temp_in = 450, &		! experiment temperature [K]
	surften_in = 0.555-(0.12/1000)*(temp_in-430), &		! indium's surface tension at temp_in [N/m]
	rho_in = 1000*(7.1295-0.67987e-3*(temp_in-273.15))	! indium's density [kg/m^3]

	! global arrays to pass between subroutines for old and new pts
	!--------------------------------
	INTEGER, DIMENSION(:,:), ALLOCATABLE,PUBLIC :: alreadyChecked, insideMarker 
	  ! for 'surfacetrack' levelsets...
	REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC	:: NewCoords
	  ! [x, y, pt val, pt type, blob#, blob_xc, blob_yc, flatpan]
	REAL*8, DIMENSION(SmallInitAry,3), PUBLIC	:: ChokeCoord 	
	  ! lists choke points and how many times called
	REAL*8, DIMENSION(:,:), ALLOCATABLE, PUBLIC	:: phi			
	  ! level set values and tying Coords->NewCoords

	! Type description
	!--------------------------------
	TYPE PANEL
		! 1=Neumann, 0=Dirichlet
		REAL*8	:: type1
		! panel z start, panel z end, panel r start, panel r end, 
		!   amount of charge (mu, SOLVING for), znormal, rnormal, middle z, 
		!   middle r, panel length, potential at the middle of the panel
		REAL*8	:: z0, z1, r0, r1, str, znrm, rnrm, midz, midr, length, midPot
		! relative points along each panel in r and z and scaled 
		!   potential=array of rpt/r_pan_center and changing potential along the way
		REAL*8, DIMENSION(npoints) :: zpoints, rpoints, rscaled
	END TYPE PANEL
END MODULE global
