! --------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  Johns Hopkins University
!
!  VICAR3Dp, a parallelized version of VICAR3D.
!  VICAR3D, a viscous, Cartesian, 3D flow solver.
!
!  This is a contineously developing project.
!
!  Starting Developers:
!  Rajat Mittal
!  Fady Najjar
!
!  Other contributing programmers:
!     Haibo Dong
!     Haoxiang Luo
!     Meliha Bozkurttas
!     Qian Xue
!     Rupesh Babu K. A.
!     Xudong Zheng 
!     Reza Ghias 
!     S. A. Mohsen Karimian
!     Vijay Vedula
!     Jung-Hee Seo
!
!  Filename: AMODULE.PAR.F90
!  Last Modified : 2015.01.12
!  by JHSeo
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following modules:
!     global_parameters
!     flow_parameters
!     grid_arrays
!     flow_arrays
!     boundary_arrays
!     multiuse_arrays
!     pressure_arrays
!     nlold_arrays
!     solver_arrays
!     solver_ad_arrays
!     GCM_arrays
!     unstructured_surface_arrays
!     probe_parameters
!     stat_arrays
!     blasius_profile
!     mg_parameters
!     mg_arrays
!     usr_module
!     stat_vort_arrays
!     tahoe_parameters
! --------------------------------------------------------------------



MODULE global_parameters
   
    IMPLICIT NONE
 
    INTEGER, PARAMETER :: CGREAL = SELECTED_REAL_KIND(P=14,R=30) 
   
    INTEGER, PARAMETER :: VISCOUS_FLOW             = 1, &
                          POTENTIAL_FLOW           = 2
 
    INTEGER, PARAMETER :: UNIFORM_GRID             = 1, &
                          NONUNIFORM_GRID          = 2
    
    INTEGER, PARAMETER :: BC_TYPE_DIRICHLET        = 1, & 
                          BC_TYPE_ZERO_GRADIENT    = 2, &
                          BC_TYPE_PULSATILE_INFLOW = 3, &  
                          BC_TYPE_SYMMETRY         = 4, &
                          BC_TYPE_PERIODIC         = 5, &
                          BC_TYPE_USER_SPECIFIED   = 6, &
                          BC_TYPE_SHEAR            = 7, &
                          BC_TYPE_USER_FLOW        = 8, &        ! JHSeo2016
                          BC_TYPE_USER_PROFILE        = 9,  &
                          BC_TYPE_USER_PROFILE_RESCALE= 10, &
                          BC_TYPE_USER_FLOW_PROFILE   = 17, &
                          BC_TYPE_USER_RECTANGULAR = 19   ! Kourosh
      
    INTEGER, PARAMETER :: IT_SOLVER_TYPE_LSOR      = 1, &
                          IT_SOLVER_TYPE_MG        = 2, &
                          BICGSTAB_CPU             = 3            ! Vijay - 2.4.0 BICGSTAB !

    INTEGER, PARAMETER :: FIXED_BOUNDARY           = 1, & 
                          MOVING_BOUNDARY          = 2 
    INTEGER, PARAMETER :: UNCOMBINED               = 0, &  ! Chao
                          COMBINED                 = 1 

    INTEGER, PARAMETER :: STATIONARY               = 0, &
                          FORCED                   = 1, & 
                          FLOW_INDUCED             = 2, &
                          PRESCRIBED               = 3, &
                          FEA_FLOW_STRUC_INTERACTION = 4, & !SER_TO_PAR, QX, CH1     
                          EIGEN_MOTION             = 5, &
                          FEA_TAHOE                = 6, &      ! Added by Rajneesh
                          USER_MOTION              = 7, &      ! JHSeo2016
                          WAVE_MOTION              = 8, &      
                          USER_VEL                 = 9, &
                          ML_MOTION                =10
 

    INTEGER, PARAMETER :: PBC_DIRICHLET            = 1, &!added by H. Luo
                          PBC_NEUMANN              = 2

    INTEGER, PARAMETER :: ADAMS_BASHFORTH2         = 1, &!added by H. Luo
                          CRANK_NICOLSON1          = 2, &
                          CRANK_NICOLSON2          = 3   

    INTEGER, PARAMETER :: NONPOROUS_AND_NONSLIP    = 0, &
                          POROUS_OR_SLIP           = 1

    INTEGER, PARAMETER :: NONE                     = 0, &
                          GENERAL                  = 1, &
                          CANONICAL                = 2

    INTEGER, PARAMETER :: ELLIPTIC_CYLINDER        = 1, &
                          GENERAL_CYLINDER         = 2, &
                          ELLIPSOID                = 3, &
                          UNSTRUCTURED_SURFACE     = 4

    INTEGER, PARAMETER :: OPEN_MEMBRANE            = 1, &
                          CLOSED_MEMBRANE          = 2, &    ! JHSeo2016
                          DIFF_MEMBRANE            = 3

    INTEGER, PARAMETER :: SOLID_BODY               = 1, &
                          MEMBRANE                 = 2     
    
    INTEGER, PARAMETER :: BODY_DIM2                = 2, &
                          BODY_DIM3                = 3
    
    INTEGER, PARAMETER :: INTR_BOUND_NONE          = 0, &
                          INTR_BOUND_PRESENT       = 1

    INTEGER, PARAMETER :: NO_VAN_KAN               = 0, &
                          VAN_KAN                  = 1

!    INTEGER, PARAMETER :: TECPLOT                  = 1, & 
!                          FIELDVIEW                = 2 

!    INTEGER, PARAMETER :: IBLANK_READ              = 1

    INTEGER, PARAMETER :: DIM_2D                   = 2, &
                          DIM_3D                   = 3

!    INTEGER, PARAMETER :: IBLANK_USED              = 1

    INTEGER, PARAMETER :: IBLANK_SLOW              = 0, &
                          IBLANK_FAST              = 1
    
    INTEGER, PARAMETER :: NO_INTERNAL_BOUNDARY     = 0, &
                          SSM_METHOD               = 1, &
                          GCM_METHOD               = 2
    
!    INTEGER, PARAMETER :: INVISCID                 = 1                  

    INTEGER, PARAMETER :: ICOORD                   = 1, &
                          JCOORD                   = 2, &
                          KCOORD                   = 3

    INTEGER, PARAMETER :: STATS_NONE               = 0

    INTEGER, PARAMETER :: STDOUT                   = 6

    INTEGER, PARAMETER :: INACTIVE                 = 0, &
                          ACTIVE                   = 1, &
                          ERR_NONE                 = 0

!   File Unit Numbers
!   -----------------
    INTEGER, PARAMETER :: ifuMarkerTRI             = 145
                          
!--------Added For FEA---------------------------------  ! SER_TO_PAR, QX, CH2
    INTEGER, PARAMETER :: PLANE_STRESS             = 1,&
                          PLANE_STRAIN             = 2,&
                          GENERAL_3DBODY           = 3

    INTEGER, PARAMETER :: STATIC_ANALYSIS          = 1,&
                          DYNAMIC_RESPONSE_CD      = 2,&
                          DYNAMIC_RESPONSE_NEWMARK = 3,&
                          DYNAMIC_CHARACTER_INVERSE = 4,&
                          DYNAMIC_CHARACTER_SURFACE = 5
    INTEGER, PARAMETER :: PENALTY_CONTACT_MODEL     = 1,&
                          GEO_CONTACT_MODEL         = 2
                           
    INTEGER, PARAMETER :: SCALAR_ON                   = 1, &
                          SCALAR_OFF                  = 0
    INTEGER, PARAMETER :: ADIABATIC                = 1, &
                          ISOTHERMAL               = 2, &      ! JHSeo2016
                          DISTRIBUTE               = 3
    INTEGER, PARAMETER :: UNIFORM_INIT_T           = 0, &
                          NONUNIFRORM_INIT_T       = 1
    INTEGER, PARAMETER :: MOMENTUM           = 0, &
                          SCALAR             = 1               ! JHSeo2016		
    INTEGER, PARAMETER :: CONSERVATIVE       = 1, &
                          NON_CONSERVATIVE   = 2

   END MODULE global_parameters
!------------------------------------------------------

   MODULE flow_parameters 

    USE global_parameters
    
    IMPLICIT NONE
 
    LOGICAL :: monitorON, monitorIT, dumpIT
    LOGICAL :: ImtheBOSS

    REAL(KIND=CGREAL), PARAMETER :: sidw  = 2.0_CGREAL ! parameter used in IDW interpolation

    REAL(KIND=CGREAL) :: pi
    REAL(KIND=CGREAL) :: KillFor2D
    
#ifdef BL_USE_MPI 
    REAL(KIND=CGREAL) :: MainTimeStart, MainTimeEnd
#else
#ifdef BL_USE_MPI3 
    REAL(KIND=CGREAL) :: MainTimeStart, MainTimeEnd
#else
    INTEGER :: MainTimeStart, MainTimeEnd
#endif
#endif
    REAL(KIND=CGREAL) :: system_wall_time, exit_safe_time
    INTEGER :: last_stop_ntime, iSAFE_EXIT
    LOGICAL :: exit_safe

    INTEGER :: AD_SOLVER  ! added by JHSeo for the AD solver switching

!   Common domain and MPI variables for parallel and serial versions (SAMK)
!   -----------------------------------------------------------------------
    INTEGER :: myRank
    INTEGER :: Np(2)
    INTEGER :: myCoords(2)
    INTEGER :: myIs, myIe, myJs, myJe
    INTEGER :: myILL, myIUL, myJLL, myJUL
    INTEGER :: myImin, myImax, myJmin, myJmax
    INTEGER :: Ngl
!   -----------------------------------------------------------------------

    INTEGER  :: nElementCheckIBLANK      ! JHSeo2016
    INTEGER  :: IC_in, ifuVelIn, idirectforcing  ! JHSeo
    INTEGER  :: iacou, Nsa, iacou_rest, MEMB_CORONA_SIZE, ipart, iflow_solver,idmemb    ! additional input parameters JHSeo
    INTEGER  :: ifufullQ, ifumarkerfullQ, StackSize, StackStart, nStack, Stackcount, ndimfullQ, ifullQ, MarkerfullQ ! added for UTIL_FULLQ JHSeo
    INTEGER  :: iCC, iSSMP, iDragLift, iMergeType, MOMENTUM_CUTCELL    ! added for CutCell JHSeo
    INTEGER  :: nread
    INTEGER  :: ndim
    INTEGER  :: flow_type
    INTEGER  :: nx_GLBL,ny_GLBL
    INTEGER  :: nxc_GLBL,nyc_GLBL
    INTEGER  :: nx,ny,nz
    INTEGER  :: nxc,nyc,nzc
    INTEGER  :: xgrid_unif,ygrid_unif,zgrid_unif
    INTEGER  :: FDprof
    INTEGER  :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
    INTEGER  :: no_tsteps,nmonitor,ndump,nrestart,nstat,nmonitor_probe_liftdrag,ninit, qdump
    INTEGER  :: iverbose  ! Added By Rajneesh
!   INTEGER  :: format_dump
    INTEGER  :: ntime, ntime_start, ntime_skip
    INTEGER  :: it_solver_type, itermax
!   INTEGER  :: mlev,iwmg,mlw
    REAL(KIND=CGREAL) :: Fr
    INTEGER  :: restart_FIM
    INTEGER  :: nstart_FIM
    INTEGER  :: boundary_motion, nBody, nBody_Solid, nBody_Membrane, nElementCheckBI, nGroup_Combined
    INTEGER, DIMENSION(:), ALLOCATABLE  :: combined_Group_index
    INTEGER, DIMENSION(:), ALLOCATABLE  :: Surface_Integral    
!   INTEGER  :: ssm_Run

    INTEGER  :: boundary_formulation, extended_outflow
    INTEGER  :: iterMaxPoisson
    INTEGER  :: body_type

    character(50) :: ParLogName

    INTEGER  :: ifuInput, ifuParLog, ifuIblankIn, &
                ifuRstrtFlowIn, ifuRstrtFlowOut,  &
                ifuBodyIn, ifuBodyOut,            &
                ifuRstrtBodyIn, ifuRstrtBodyOut,  &
                ifuMarkerIn, ifuMarkerOut,        &
                ifuProbeIn, ifuProbeOut,          &
                ifuUnstrucSurfIn, ifuUnstrucSurfOut, &
                ifuOpenMembraneEdgeIn, ifurunningstatus, &
                ifuDyeIn, ifuScalarIn,ifuRstrtScalarIn,ifuRstrtScalarOut
                

    INTEGER  :: ifuStatOut
    INTEGER  :: ifuDragOut, ifuFreshCellOut

    INTEGER  :: non_inertial=0
    
    INTEGER  :: ifuFlux ! used for FEA 
    INTEGER  :: ifuVol  ! JHSeo2016	
                                                                                      ! SER_to_PAR. QX. CH3

    INTEGER  :: internal_boundary_present, nPtsMax, iblankFast
    INTEGER  :: idxRstrt, idxStat!!, indexStatVort
    INTEGER  :: frac_step_type
   
    INTEGER, DIMENSION(:), ALLOCATABLE  :: nPtsBodyMarkerOrig,canonical_body_type, &
                                           body_dim,boundary_motion_type,combined_type,wall_type, &
                                           membrane_type,unstruc_surface_type
 
    INTEGER, DIMENSION(:), ALLOCATABLE  :: nPtsBodyMarker

    INTEGER, DIMENSION(:), ALLOCATABLE  :: n_theta,n_phi

    INTEGER, DIMENSION(:), ALLOCATABLE  :: ntimePerCycle ! added by Haibo

!------------- new arrays -- added by Haibo
!    INTEGER, DIMENSION(:),        ALLOCATABLE :: imv,ipv,immv,ippv
!    INTEGER, DIMENSION(:),        ALLOCATABLE :: jmv,jpv,jmmv,jppv
!    INTEGER, DIMENSION(:),        ALLOCATABLE :: kmv,kpv,kmmv,kppv
!--------------------------------

!    LOGICAL            :: readIblankFlag

    REAL(KIND=CGREAL)  :: global_flux ! JHSeo2016
    REAL(KIND=CGREAL)  :: wave_ak, wave_c, wave_amp, wave_x_shift, wave_xc, wave_yc ! JHSeo
    REAL(KIND=CGREAL)  :: wave_CX1, wave_CX2, wave_CL
    REAL(KIND=CGREAL)  :: xout,yout,zout,Xext,DampFact
    REAL(KIND=CGREAL)  :: uinit,vinit,winit
    REAL(KIND=CGREAL)  :: ux1,ux2,vx1,vx2,wx1,wx2
    REAL(KIND=CGREAL)  :: uy1,uy2,vy1,vy2,wy1,wy2
    REAL(KIND=CGREAL)  :: uz1,uz2,vz1,vz2,wz1,wz2
    REAL(KIND=CGREAL)  :: freq_ux1,freq_vx1,freq_wx1,beta_x1,dir_x1  ! JHSeo2016
    REAL(KIND=CGREAL)  :: freq_ux2,freq_vx2,freq_wx2,beta_x2,dir_x2
    REAL(KIND=CGREAL)  :: freq_uy1,freq_vy1,freq_wy1,beta_y1,dir_y1
    REAL(KIND=CGREAL)  :: freq_uy2,freq_vy2,freq_wy2,beta_y2,dir_y2
    REAL(KIND=CGREAL)  :: freq_uz1,freq_vz1,freq_wz1,beta_z1,dir_z1
    REAL(KIND=CGREAL)  :: freq_uz2,freq_vz2,freq_wz2,beta_z2,dir_z2

    REAL(KIND=CGREAL)  :: freq_px1, beta_px1
    REAL(KIND=CGREAL)  :: freq_px2, beta_px2
    REAL(KIND=CGREAL)  :: freq_py1, beta_py1
    REAL(KIND=CGREAL)  :: freq_py2, beta_py2

  
    REAL(KIND=CGREAL)  :: re,dt,reinv,dtinv
    REAL(KIND=CGREAL)  :: restol, restolPoisson, omega, omega_adv
    REAL(KIND=CGREAL)  :: time
    REAL(KIND=CGREAL)  :: bodyInterceptWeight, imagePointWeight, probeLengthNormalized
    REAL(KIND=CGREAL)  :: gcmFlag

    REAL(KIND=CGREAL)  :: areax1,areax2,  &
                          areay1,areay2,  &
                          areaz1,areaz2

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: radiusx,radiusy,radiusz
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: alpha,cosalpha,sinalpha
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: phase,cosphase,sinphase
     
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: vxcent,vycent,vzcent,  &
                                                   angvx,angvy,angvz,     &
                                                   xcent,ycent,zcent

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: x_rot_cent, y_rot_cent, z_rot_cent
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: vxcentOld,vycentOld,vzcentOld
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: vxcentPre,vycentPre,vzcentPre ! non_inertial

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: dvxc,dvyc,dvzc          ! Aitken
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: dvxcOld,dvycOld,dvzcOld
    
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvxOld,angvyOld,angvzOld
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: vxcent_combinedOld,vycent_combinedOld,vzcent_combinedOld
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvx_combinedOld,angvy_combinedOld,angvz_combinedOld 
       
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: xcent_combined, ycent_combined, zcent_combined
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: vxcent_combined, vycent_combined, vzcent_combined
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvx_combined, angvy_combined, angvz_combined
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvx_combined_old, angvy_combined_old, angvz_combined_old             

! for new motion -- Added by Haibo
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angvx_old,angvy_old,angvz_old
     
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: xcentinit,ycentinit,zcentinit,        &
                                                   vxcentTrans,vycentTrans,vzcentTrans,  &
                                                   ampx,ampy,ampz,freqx,freqy,freqz
                                                        
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angvxinit,angvyinit,angvzinit,  &
                                                   ampangx,ampangy,ampangz,        &
                                                   freqangx,freqangy,freqangz
    INTEGER,DIMENSION(:), ALLOCATABLE :: i_fixed, fixed_mother_body                              

! wall velocity valiables
    INTEGER,           DIMENSION(:),     ALLOCATABLE  :: mMinWallVel, mMaxWallVel
    REAL(KIND=CGREAL), DIMENSION(:),     ALLOCATABLE  :: ampVelX, ampVelY, ampVelZ
    REAL(KIND=CGREAL), DIMENSION(:),     ALLOCATABLE  :: freqVelX, freqVelY, freqVelZ
    REAL(KIND=CGREAL), DIMENSION(:),     ALLOCATABLE  :: phaseVelX, phaseVelY, phaseVelZ

    REAL(KIND=CGREAL)  :: area_left,area_right,     &
                          area_bot,area_top,        &
                          area_back,area_front,     &
                          outflow_area, inflow_area, prim_left, outflow_area2  ! JHSeo2016

    REAL(KIND=CGREAL)  :: alfa       ! Weighting factor for hybrid scheme - added by Rupesh

    REAL(KIND=CGREAL)  :: vPert2Dto3D ! for 2D-3D initial random perturbations
  
! For Flow-Induced Motion. --- Added by veera

    REAL(KIND=CGREAL) ,DIMENSION(:),ALLOCATABLE  :: xcentConstr, ycentConstr, zcentConstr ! Centroid Constraint Flag
    
    REAL(KIND=CGREAL) :: density_fluid 
    
    REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: density_solid
   
    INTEGER            :: pbcx1,pbcx2, pbcy1,pbcy2, pbcz1,pbcz2  ! Added by H. Luo
	INTEGER            :: puserx1, puserx2, pusery1, pusery2, puserz1, puserz2  ! JHSeo2016
    REAL(KIND=CGREAL)  :: pppx1,pppx2, pppy1,pppy2, pppz1,pppz2  !
    INTEGER            :: advec_scheme                           ! added for implicit scheme

    INTEGER            :: ifluxplane, jfluxplane, kfluxplane     ! used for FEA
                                                                 ! SER_TO_PAR. QX. CH4

    INTEGER            :: is_scalar_on !xzheng flag for adding dyes into flow

    INTEGER            :: ppeTotIter             ! Added by Vijay - PAT 2.3.5 !

   END MODULE flow_parameters
!------------------------------------------------------

   MODULE grid_arrays

    USE global_parameters
   
    IMPLICIT NONE

!   Global Grid Variables
!   ---------------------
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: x,y,z
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: xc,yc,zc
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dx,dy,dz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxinv,dyinv,dzinv
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxc,dyc,dzc
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxcinv,dycinv,dzcinv
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: fx,fy,fz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: damper
   
   END MODULE grid_arrays
!------------------------------------------------------

   MODULE flow_arrays

    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: u,v,w
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: face_ue,face_vn,face_wf
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: face_uw,face_vs,face_wb
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxu,bcxv,bcxw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcyu,bcyv,bcyw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bczu,bczv,bczw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: viscTot
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxvisc,bcyvisc,bczvisc

    REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: uGhost, vGhost, wGhost
    REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: pGhost

    REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: div, temp1, temp2, temp3, vorQ, vdiss  ! JHSeo2016

   END MODULE flow_arrays
!------------------------------------------------------

   MODULE boundary_arrays

    USE global_parameters
    
    IMPLICIT NONE

    INTEGER  :: nFresh, num_fresh, num_dead

    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iblank, fresh_cell, ghostCellMark, boundCell, fresh_cell_temp, boundCellEdge
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ghostCellSolid, ghostCellMemb
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: bodyNum, bodyNumOld

    INTEGER, DIMENSION(:,:,:),   ALLOCATABLE :: iblank_solid
    INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: iblank_memb

    INTEGER, DIMENSION(:,:,:),   ALLOCATABLE :: INTENSITY ! JHSeo2016	

!    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: pot_flag
    REAL(KIND=CGREAL), DIMENSION(:),     ALLOCATABLE :: iBound, jBound, kBound
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: iup , ium , jup , jum , kup,  kum, exp_weight
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: iupp, iumm, jupp, jumm, kupp, kumm              ! For 2nd Upwinding - Added by Rupesh

    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarker,yBodyMarker,zBodyMarker, &
                                                        uBodyMarker,vBodyMarker,wBodyMarker, &
                                                        uBodyMarker_r,vBodyMarker_r,wBodyMarker_r, & ! non_inertial, Zhuoyu
                                                        sBodyMarker, pBodyMarker, pBodyMarker1, & !xudong added June 3rd 
                                                        ppBodyMarker, upBodyMarker, vpBodyMarker, wpBodyMarker  ! JHSeo2016

    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: theta_group
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: angv_roll,angv_yaw,angv_pitch      
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: R11,R12,R13,R21,R22,R23,R31,R32,R33 
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: C11,C12,C13,C21,C22,C23,C31,C32,C33                                                         
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyMarker_rel,vBodyMarker_rel,wBodyMarker_rel
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyMarker_rel_old,vBodyMarker_rel_old,wBodyMarker_rel_old

    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgradx1,pgradx2
    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgrady1,pgrady2 
    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgradz1,pgradz2 
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: boundPresSource

    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iblankUndecided, iblankTemp

   END MODULE boundary_arrays
!------------------------------------------------------

   MODULE multiuse_arrays

    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nlu,nlv,nlw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: uTilde,vTilde,wTilde
   
   END MODULE multiuse_arrays
!------------------------------------------------------

   MODULE pressure_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: p, pPrime
   
   END MODULE pressure_arrays
!------------------------------------------------------

   MODULE nlold_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nluold,nlvold,nlwold
   
   END MODULE nlold_arrays
!------------------------------------------------------

   MODULE solver_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amx,apx,acx
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amy,apy,acy
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amz,apz,acz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: rhs,dummy
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: face1, face2
   
   END MODULE solver_arrays
!------------------------------------------------------

   MODULE solver_ad_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amx_ad,apx_ad
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amy_ad,apy_ad
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amz_ad,apz_ad
   
   END MODULE solver_ad_arrays
!------------------------------------------------------

   MODULE GCM_arrays
   
    USE global_parameters
   
    IMPLICIT NONE
   
    INTEGER                                :: iRowMax, nGhost
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: incI, incJ, incK, iPvt
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestMarker,            &
                                              iGhost,jGhost,kGhost,     &
                                              iCellIndex,jCellIndex,kCellIndex
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: iFresh,jFresh,kFresh      &
                                             ,iFreshCellIndex,jFreshCellIndex,kFreshCellIndex
!    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestMarkerFresh
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestElementFresh
!    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: iBodyRank
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: iCellIndexS, jCellIndexS, kCellIndexS

!    REAL(KIND=CGREAL), DIMENSION(2)                :: det  
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: work
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: closestMarkerRatio, &
                                                      xBodyInterceptTang, &
                                                      yBodyInterceptTang, &
                                                      zBodyInterceptTang, &
                                                      xBodyInterceptNorm, &
                                                      yBodyInterceptNorm, &
                                                      zBodyInterceptNorm, &
                                                      xBodyIntercept,     &
                                                      yBodyIntercept,     &
                                                      zBodyIntercept,     &
                                                      xImagePoint,        &
                                                      yImagePoint,        &
                                                      zImagePoint,        &
                                                      probeLength,        &
                                                      uBodyIntercept,     &
                                                      vBodyIntercept,     &
                                                      wBodyIntercept,     &
                                                      pBodyIntercept,     &
                                                      dpdnBodyIntercept,  &
                                                      dpdtBodyIntercept

!    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: xBIG,yBIG,ZBIG

    REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: xBodyInterceptFresh,     &
                                                      yBodyInterceptFresh,     &
                                                      zBodyInterceptFresh,     &
                                                      uBodyInterceptFresh,     &
                                                      vBodyInterceptFresh,     &
                                                      wBodyInterceptFresh,     &
                                                      xBodyInterceptNormFresh, &
                                                      yBodyInterceptNormFresh, &
                                                      zBodyInterceptNormFresh

    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: coeffGCMD, coeffGCMN,         &
                                                      vanMatrixD, vanMatrixN,       &
                                                      coeffGCMFreshD, coeffGCMFreshN, coeffGCMFreshT !,dSFaceProject, &   
!                                                      xBodyCentroid,yBodyCentroid,  &
!                                                      sBodyCentroid,                &
!                                                      xCentroidTang,yCentroidTang, &
!                                                      xCentroidNorm,yCentroidNorm
 
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: probeLengthS,probeLengthNormalizedS, &
                                                      imagePointWeightS,                   &
                                                      xImagePointS,yImagePointS,zImagePointS

    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: coeffGCMDS, coeffGCMNS!, &
!                                                      vanMatrixDS, vanMatrixNS

   END MODULE GCM_arrays
!------------------------------------------------------

   MODULE unstructured_surface_arrays

    USE global_parameters
    
    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemNormX,      triElemNormY,      triElemNormZ,      triElemArea
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemCentX,      triElemCentY,      triElemCentZ
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemTang1X,     triElemTang1Y,     triElemTang1Z
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemTang2X,     triElemTang2Y,     triElemTang2Z
    REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: pointOutsideBodyX, pointOutsideBodyY, pointOutsideBodyZ, surfArea

    INTEGER, DIMENSION(:,:,:),         ALLOCATABLE :: triElemNeig
    INTEGER, DIMENSION(:),             ALLOCATABLE :: totNumTriElem
    INTEGER, DIMENSION(:),             ALLOCATABLE :: closestElementGC!,cElementG
    INTEGER, DIMENSION(:,:),           ALLOCATABLE :: edge_node

    REAL (KIND=CGREAL)                             :: normDirFlag

   END MODULE unstructured_surface_arrays
!------------------------------------------------------

   MODULE probe_parameters

    IMPLICIT NONE
    
    INTEGER                           :: nProbe, nProbe_Marker  ! JHSeo2016
    INTEGER, DIMENSION(:),ALLOCATABLE :: iProbe, jProbe, kProbe, Body_probe, Marker_probe
    
    LOGICAL, DIMENSION(:),ALLOCATABLE :: myProbe

   END MODULE probe_parameters
!------------------------------------------------------

   MODULE stat_arrays

    USE global_parameters

        IMPLICIT NONE

         INTEGER                                        :: statCtr
         REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uAv,vAv,wAv,pAv
         REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uvAv,vwAv,uwAv
         REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uuAv,vvAv,wwAv

   END MODULE stat_arrays
!------------------------------------------------------

    MODULE blasius_profile

        USE global_parameters

        IMPLICIT NONE

!        REAL(KIND=CGREAL)                          :: cavity_H,slot_H,ddratio,d,delta,uinf
!        REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:) :: eta,u_blasius
!        INTEGER                                    :: l,junk,i_start

    END MODULE blasius_profile
!------------------------------------------------------

    MODULE mg_parameters
 
    USE global_parameters
    USE flow_parameters
 
    IMPLICIT NONE
 
    INTEGER :: mgLevels_X, mgLevels_Y, mgLevels_Z
!    INTEGER :: mgcyclex, mgcycley, mgcyclez, infoconv, incrlev
    INTEGER :: mgcyclex, mgcycley, mgcyclez, infoconv
    INTEGER :: iterFinest, iterInter, iterCoarsest
!    INTEGER :: ittt1, nCount
 
    INTEGER :: iRedBlack, TNcolorX, TNcolorY, TNcolorZ, iStep, jStep, kStep    !new for Redblack LSOR
    
    integer :: nCoarseGrids

    integer :: maxIterInter                      ! Added by Vijay - PAT 2.3.5 ! 
   END MODULE mg_parameters
 
!------------------------------------------------------
   MODULE mg_arrays

    USE global_parameters
 
    IMPLICIT NONE
 
    TYPE :: MGXtype
!   ===============
      INTEGER :: nx_GLBL , nx
      INTEGER :: nxc_GLBL, nxc

      INTEGER :: myIs, myIe

      INTEGER, DIMENSION(:,:,:), POINTER :: iblank

      REAL(KIND=CGREAL), DIMENSION(:),     POINTER :: x, xc, dxinv, dxcinv
      REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: rhs, res, phi

#ifdef BL_USE_MPI 
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#else
#ifdef BL_USE_MPI3 
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#endif
#endif
    END TYPE MGXtype
    
    TYPE :: MGYtype
!   ===============
      INTEGER :: ny_GLBL , ny
      INTEGER :: nyc_GLBL, nyc

      INTEGER :: myJs, myJe

      INTEGER, DIMENSION(:,:,:),POINTER  :: iblank

      REAL(KIND=CGREAL), DIMENSION(:),     POINTER :: y, yc, dyinv, dycinv
      REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: rhs, res, phi

#ifdef BL_USE_MPI 
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#else
#ifdef BL_USE_MPI3 
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#endif
#endif
    END TYPE MGYtype
    
    TYPE :: MGZtype
!   ===============
      INTEGER :: nz, nzc

      INTEGER, DIMENSION(:,:,:), POINTER:: iblank

      REAL(KIND=CGREAL), DIMENSION(:),     POINTER:: z, zc, dzinv, dzcinv
      REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: rhs, res, phi

#ifdef BL_USE_MPI 
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#else
#ifdef BL_USE_MPI3 
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#endif
#endif
    END TYPE MGZtype
    
    TYPE(MGXtype), ALLOCATABLE :: MGX(:)
    TYPE(MGYtype), ALLOCATABLE :: MGY(:)
    TYPE(MGZtype), ALLOCATABLE :: MGZ(:)
    
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE              :: iblank_MG
!!    INTEGER, DIMENSION(:,:,:), ALLOCATABLE              :: ghostcellMark_MG, iblank_MG
 
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE    :: ium_mg, iup_mg, &
                                                           jum_mg, jup_mg, &
                                                           kum_mg, kup_mg

!------- new arrays ----------
!   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: iup_mg_Outer, jup_mg_Outer, kup_mg_Outer
!   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: ium_mg_Outer, jum_mg_Outer, kum_mg_Outer
 
   END MODULE mg_arrays
     
!------------------------------------------------------
MODULE usr_module

    USE global_parameters
     
   REAL(KIND=CGREAL) :: density_ratio
   REAL(KIND=CGREAL) :: lScale,vScale
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: I_XX_COMBINED,I_YY_COMBINED,I_ZZ_COMBINED,I_XY_COMBINED,I_YZ_COMBINED,I_XZ_COMBINED
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: non_dim_volume,volume
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: I_XX,I_YY,I_ZZ,I_XY,I_YZ,I_XZ
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: non_dim_mass
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: moment_x_combined,moment_y_combined,moment_z_combined
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: force_x_combined,force_y_combined,force_z_combined
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: non_dim_mass_combined
   REAL(KIND=CGREAL),DIMENSION(1:3,1:3) :: nonDimM_I_combined, invMI_combined
   REAL(KIND=CGREAL) :: vxcent_combined_prev,vycent_combined_prev,vzcent_combined_prev
   REAL(KIND=CGREAL) :: angvx_combined_prev, angvy_combined_prev, angvz_combined_prev   
   REAL(KIND=CGREAL) :: vxcent_prev,vycent_prev,vzcent_prev
   REAL(KIND=CGREAL) :: angvx_prev, angvy_prev, angvz_prev
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: force_x,force_y,force_z
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: moment_x,moment_y,moment_z   
   REAL(KIND=CGREAL),DIMENSION(1:3,1:3) :: invMI
   REAL(KIND=CGREAL),DIMENSION(:,:,:),ALLOCATABLE :: nonDimM_I
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: scx,scy,scz
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: scmx,scmy,scmz
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: shear_x,shear_y,shear_z

END MODULE usr_module 

!------------------------------------------------------
MODULE stat_vort_arrays

    USE global_parameters

    IMPLICIT NONE

    INTEGER                                        :: statCtrv
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: oxAv,oyAv,ozAv
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: oxoxAv,oyoyAv,ozozAv

END MODULE stat_vort_arrays
!------------------------------------------------------

!=====Added for FEA====================================

   MODULE finite_element_parameters  ! SER_TO_PAR. QX. CH5

    USE global_parameters

    IMPLICIT NONE

    INTEGER :: fea_nbody, fea_nnode, fea_nelementmax, fea_nptmax, fea_mbandmax, fea_nmat
    INTEGER :: fea_nfixptmax, fea_nloadmax, fea_neignvmax
    INTEGER :: fea_nndof, fea_nnstress, fea_solvertype, fea_probtype, fea_outputtype, fea_nva
    INTEGER :: fea_nedof, fea_ngdof
    INTEGER :: fea_ntime, fea_itime
    INTEGER :: fea_inputfile, fea_meshfile, fea_outputfile, fea_conloadfile, fea_solutionfile,&
               fea_staticoutput, fea_initdisvel, fea_correstablefile, fea_restartfilein, fea_restartfileout,&
               fea_stressfile, fea_probein, fea_ifuProbeOut, fea_Mfile, fea_keffectfile, fea_afile
    INTEGER :: fea_boundary_flag, fea_contactprobflag
    INTEGER :: fea_nd
    INTEGER :: fea_nfilterstart, fea_nfilterend, fea_dtratio
    INTEGER :: fea_readM, fea_readkeffect, fea_readinita
    INTEGER :: fea_contact_model,fea_trans_dir

    REAL(KIND=CGREAL)  :: fea_gravkey, fea_grav, fea_time
    REAL(KIND=CGREAL)  :: fea_freq, fea_cc1, fea_cc2, fea_dt, fea_beta, fea_gamma
    REAL(KIND=CGREAL)  :: fea_nmc0, fea_nmc1, fea_nmc2, fea_nmc3, fea_nmc4, fea_nmc5, fea_nmc6, fea_nmc7
    REAL(KIND=CGREAL)  :: fea_omegaload
    REAL(KIND=CGREAL)  :: fea_coeff_penalty

   END MODULE  finite_element_parameters

!------------------------------------------------------
   MODULE finite_element_arrays  ! SER_TO_PAR. QX. CH6

    USE global_parameters

    IMPLICIT NONE

    INTEGER, DIMENSION(:), ALLOCATABLE :: fea_nelement, fea_npt, fea_mband,fea_nfixpt, fea_nload, fea_muv, fea_neignv
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: fea_ielement, fea_ifixed, fea_iload
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: markertofeapointtable, featomarkerpointtable


    REAL(KIND=CGREAL),DIMENSION(:, :, :), ALLOCATABLE :: fea_cood, fea_vfixed, fea_vload
    REAL(KIND=CGREAL),DIMENSION(:, :, :), ALLOCATABLE :: fea_gmm, fea_gkm, fea_gcm, fea_keffect
    REAL(KIND=CGREAL),DIMENSION(:, :, :), ALLOCATABLE :: fea_veignv
    REAL(KIND=CGREAL),DIMENSION(:, :),    ALLOCATABLE :: fea_vmati, fea_ieignv
    REAL(KIND=CGREAL),DIMENSION(:, :),    ALLOCATABLE :: fea_d, fea_v, fea_a, fea_gp, fea_gp_old, fea_gu

    INTEGER, DIMENSION(:), ALLOCATABLE :: fea_icontactdir, fea_icontactplane, fea_icontactdirnorm
    INTEGER, DIMENSION(:), ALLOCATABLE :: fea_contactflag, fea_ncontactpoint, fea_ncontactsurf, fea_ncontactside
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: fea_icontactpoint, fea_icontactsurf, fea_icontactside

    REAL(KIND=CGREAL),DIMENSION(:, :),    ALLOCATABLE :: fea_original_d, fea_original_v, fea_original_a
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: fea_penaltycoeff
   
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: tempt,keffecttempt
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: tempta
 
   END MODULE finite_element_arrays

!------------------------------------------------------
   MODULE fea_probe_parameters  ! SER_TO_PAR. QX. CH7

    IMPLICIT NONE

    INTEGER                           :: fea_nProbe
    INTEGER, DIMENSION(:),ALLOCATABLE :: fea_iprobebody, fea_iprobenode

   END MODULE fea_probe_parameters

!------------------------------------------------------
   MODULE eigen_motion_parameters  !SER_TO_PAR. QX.CH8

    USE global_parameters

    IMPLICIT NONE

    INTEGER ::  neig_body, maxnmodes, fea_ifeigenin

   ENDMODULE eigen_motion_parameters

!-----------------------------------------------------
   MODULE eigen_motion_arrays  !SER_TO_PAR. QX. CH9

    USE global_parameters

    IMPLICIT NONE

    INTEGER, ALLOCATABLE, DIMENSION(:) :: ieig_body_table, nmodes_eig, ncontact_eig
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: icon_dir_eig, icon_plane_eig, icon_surnor_eig
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:, :, :) :: x_eig, y_eig, z_eig
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:, :) :: fre_eig, coe_eig, phase_eig
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:) :: x_coor_org, y_coor_org, z_coor_org

  ENDMODULE eigen_motion_arrays
 
!========End FEA====================================

    MODULE tahoe_parameters  ! Added by Rajneesh
 
   USE global_parameters
   IMPLICIT NONE
   
   INTEGER :: ndim_tahoe,elementType,ConstitutiveModel,TahoeSolver,nstart_tahoe,ndump_tahoe,nprobe_tahoe,&
   dtratio_tahoe,nmax,emax,FBCmax,KBCmax,Markmax, fea_boundary_flag_tahoe,itime_tahoe,KBC_flag,FBC_flag,&
   iflag_restart_tahoe,tahoe_inputfile,tahoe_fbcnodes,tahoe_kbcnodes,tahoe_markers2fem,tahoe_dataprobe,&
   tahoe_u,tahoe_Du,tahoe_DDu,tahoe_elem0,tahoe_rs,tahoe_xml,tahoe_bodydat,tahoe_geom,implicit_inputfile,&
   imconverg,time_data1,time_data2,restart_tahoebodyin,restart_tahoebodyout
   
   INTEGER, ALLOCATABLE, DIMENSION(:) :: numNodes,numElts,numMarkers,FBCnodes,KBCnodes
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: FBCnode2marker,Markers2femGridpt,nodesFBC,nodesKBC,BodyProbe,marker2elem
   
   REAL(KIND=CGREAL)  :: SolidDen,YmodulusEQ,YmodulusNEQ,PoissonR,etaDamping,muNEQ,muEQ,kappaNEQ,kappaEQ,tauBulk,tauShear,dtTahoe
   REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:) :: x0i,y0i,z0i
   REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: dispTahoeBody,veloTahoeBody,acclTahoeBody, dispTahoeBodyOld,&
   veloTahoeBodyOld,acclTahoeBodyOld
   
  ENDMODULE tahoe_parameters

!-----------------------------------------------------
MODULE cutcell_arrays
 USE global_parameters
 
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: AREA_W, AREA_E, AREA_S, AREA_N, AREA_B, AREA_F
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: VOL_CELL, CUT_AREA, cent_x,cent_y,cent_z
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: cc_nx, cc_ny, cc_nz, nlu_smallcell, div_smallcell
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: small_ue, small_uw, small_vn, small_vs, small_wf, small_wb
 INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: isolid

END MODULE  

!========End cutcell_arrays====================================
!-----------------------------------------------------
MODULE implicit_coupling_parameters  ! Added by Rajneesh
 
  USE global_parameters
   IMPLICIT NONE

   INTEGER :: implicit_coupling_flag,kimplicit,kimplicitMax,flag_underrelax,kimplicit_start,varyalphawithtime  ! JHSeo2016
   REAL(KIND=CGREAL) :: alpha_underrelax,slope_underrelax,ImplicitResidual,ImplicitResidualMax

   INTEGER :: implicit_coupling_flag_combined_motion, under_relaxation_method
   INTEGER :: convg_criterion
   INTEGER :: fim_inputfile
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarkerPrev,yBodyMarkerPrev,zBodyMarkerPrev 
   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: C11_prev,C12_prev,C13_prev,C21_prev,C22_prev,C23_prev,C31_prev,C32_prev,C33_prev  
   INTEGER :: kimplicit_FIM,kimplicitMax_FIM,flag_underrelax_FIM,kimplicit_start_FIM,varyalphawithtime_FIM
   REAL(KIND=CGREAL) :: alpha_underrelax_FIM,slope_underrelax_FIM,ImplicitResidual_FIM,ImplicitResidualMax_FIM
   REAL(KIND=CGREAL) :: timeimpl1_FIM,timeimpl2_FIM,alpha_underrelax2_FIM, alphaimplicit_FIM   

   REAL(KIND=CGREAL) :: timeimpl1,timeimpl2,alpha_underrelax2, alphaimplicit
   REAL(KIND=CGREAL) :: lambda1, lambda2 ! Aitken
   REAL(KIND=CGREAL) :: res1, res2
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iblankOld,iblank_solidOld
   INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: iblank_membOld
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: uOld,vOld,wOld,uGhostOld,vGhostOld,wGhostOld,pPrime_old      
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarkerOld,yBodyMarkerOld,zBodyMarkerOld
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: uBodyMarkerOld,vBodyMarkerOld,wBodyMarkerOld

   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: dxBodyMarker,dyBodyMarker,dzBodyMarker ! Aitken
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: dxBodyMarkerOld,dyBodyMarkerOld,dzBodyMarkerOld
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: duBodyMarker,dvBodyMarker,dwBodyMarker
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: duBodyMarkerOld,dvBodyMarkerOld,dwBodyMarkerOld

ENDMODULE implicit_coupling_parameters
   
!-----------------------------------------------------
MODULE diffusive_material_transport
 
 USE global_parameters
 IMPLICIT NONE
 
 INTEGER :: nscalars, T_discre_type, T_solver_type, reaction_ON, Restart_Scalar, T_user_src_flag  ! JHSeo2016
 INTEGER :: plt_model_ON, fibrin_model_ON
 
 INTEGER, DIMENSION(:), ALLOCATABLE :: T_inner_BCtype
 INTEGER, DIMENSION(:), ALLOCATABLE :: bcx1_scalar, bcy1_scalar, bcz1_scalar   !BOUNDARY CONDITION TYPE
 INTEGER, DIMENSION(:), ALLOCATABLE :: bcx2_scalar, bcy2_scalar, bcz2_scalar
 INTEGER, DIMENSION(:), ALLOCATABLE :: T_init_type, T_user_src_ON  ! JHSeo2016

 REAL(KIND=CGREAL) Total_T, Total_Told, T_ratio, T_ds, T_D, T_db, T_de, alfa_T, Tmax, Tmin  ! JHSeo2016
 
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: Tbcinner
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: sc
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: T_init_val
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: Tbcx1, Tbcx2, Tbcy1, Tbcy2, Tbcz1, Tbcz2

 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: T_scalar, nlT, nlTold, bcxT, bcyT, bczT
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: TGhost, Reaction_source  ! JHSeo2016

 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amx_ad_T,apx_ad_T
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amy_ad_T,apy_ad_T
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amz_ad_T,apz_ad_T

 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amx_ad_org,apx_ad_org
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amy_ad_org,apy_ad_org
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amz_ad_org,apz_ad_org 
 
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: TBodyMarker, T_user_src   ! JHSeo2016

ENDMODULE diffusive_material_transport

MODULE reaction_arrays    ! JHSeo2016

 USE global_parameters
 IMPLICIT NONE
 
 INTEGER  n_elements, n_reactions, n_coeff
 INTEGER, DIMENSION(:), ALLOCATABLE :: forward_order, backward_order, forward_coeff, backward_coeff
 INTEGER, DIMENSION(:,:), ALLOCATABLE :: forward_elements, backward_elements, Stoi_Matrix
 
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: Reaction_vector, coeff_k
  
END MODULE reaction_arrays


MODULE Plt_model      ! JHSeo2016
 USE global_parameters
  
 INTEGER iPTmu, iPTma, iPTb, iIIa, ifg, ifb
 REAL(kind=CGREAL) kadh1, kadh2, kcoh, kIIa, Cstar, kcat, Km, kp
 
END MODULE Plt_model

MODULE Rectangular_INFLOW  ! Kourosh
 USE global_parameters
 REAL(kind=CGREAL) ux_adjustCoef
END MODULE Rectangular_INFLOW

   
!******************************************************************************

   module subMGParameters                        ! Module added by Vijay - PAT 2.3.5 for subMG !
   use global_parameters
   use flow_parameters
   implicit none

   integer :: iSubMG, subMGLevelsX, subMGLevelsY
   double precision :: subMGTime, totSubMGTime, subMGCommTime, totSubMGCommTime
   logical :: iCheck

   end module subMGParameters

!******************************************************************************

   module subMGArrays                            ! Module added by Vijay - PAT 2.3.5 for subMG !
   use global_parameters
   implicit none

   type :: subMGXType
    integer :: nx, nxc, nxGLBL, nxcGLBL, myIs
    integer, dimension(:,:,:), allocatable :: iblank
    double precision, dimension(:), allocatable :: x, xc, dxInv, dxcInv
    double precision, dimension(:,:,:), allocatable :: phi, rhs, res
   end type subMGXType

   type :: subMGYType
    integer :: ny, nyc, nyGLBL, nycGLBL, myJs, myJe
    integer, dimension(:,:,:), allocatable :: iblank
    double precision, dimension(:), allocatable :: y, yc, dyInv, dycInv
    double precision, dimension(:,:,:), allocatable :: phi, rhs, res
   end type subMGYType

   type(subMGXType), dimension(:), allocatable :: subMGX
   type(subMGYType), dimension(:), allocatable :: subMGY

   integer, dimension(:,:,:), allocatable :: iblankSubMG  
   double precision, dimension(:,:,:), allocatable :: ium_subMG, iup_subMG
   double precision, dimension(:,:,:), allocatable :: jum_subMG, jup_subMG
   double precision, dimension(:,:,:), allocatable :: kum_subMG, kup_subMG

   double precision, dimension(:), allocatable :: amx_subMG,acx_subMG,apx_subMG
   double precision, dimension(:), allocatable :: amy_subMG,acy_subMG,apy_subMG
   double precision, dimension(:), allocatable :: amz_subMG,acz_subMG,apz_subMG
   double precision, dimension(:), allocatable :: rhs_subMG,dummy_subMG

   end module subMGArrays

!****************** Vijay - 2.4.0 BICGSTAB - Start ********************!

   module bicgVariables
   use flow_parameters
   implicit none

   real(kind=cgreal) :: bicgAlpha,bicgBeta,bicgOmega
   real(kind=cgreal) :: bicgRho,bicgRho1,bicgR0Vk,bicgTkSk,bicgTkTk
   real(kind=cgreal), dimension(:,:,:), allocatable :: bicgKinv
   real(kind=cgreal), dimension(:,:,:), allocatable :: bicgR0,bicgRk
   real(kind=cgreal), dimension(:,:,:), allocatable :: bicgPk,bicgYk,bicgYk_G,bicgVk
   real(kind=cgreal), dimension(:,:,:), allocatable :: bicgSk,bicgZk,bicgZk_G,bicgTk
   real(kind=cgreal), dimension(:    ), allocatable :: bicgR0_G,bicgRk_G
   real(kind=cgreal), dimension(:    ), allocatable :: bicgPk_G,bicgVk_G
   real(kind=cgreal), dimension(:    ), allocatable :: bicgSk_G,bicgTk_G

   real(kind=cgreal), dimension(:,:,:), allocatable :: bicgAMX,bicgAPX
   real(kind=cgreal), dimension(:,:,:), allocatable :: bicgAMY,bicgAPY
   real(kind=cgreal), dimension(:,:,:), allocatable :: bicgAMZ,bicgAPZ
   real(kind=cgreal), dimension(:,:,:), allocatable :: bicgAC

   end module bicgVariables

!****************** Vijay - 2.4.0 BICGSTAB - End ********************!
   
module directforcing_arrays     ! JHSeo2016
 use flow_parameters
 implicit NONE
 
 integer, dimension(:,:,:), allocatable :: mask_function
 

end module directforcing_arrays
   
   
MODULE scalar_phi ! Chao, FPM

 USE global_parameters
 IMPLICIT NONE
 
 
 INTEGER, PARAMETER :: XBC_DIRICHLET            = 1, &
                          XBC_NEUMANN              = 2 
                          
 INTEGER, PARAMETER :: PhiBC_DIRICHLET            = 1, &
                          PhiBC_NEUMANN              = 2                           
                          
                          
 INTEGER             :: VDVactive, PFactive
 INTEGER             :: n_visual, n_potential
 INTEGER             :: idirection
 INTEGER             :: iterMaxLaplace
 INTEGER             :: iterMaxPotential
 INTEGER             :: POISSON_INDEX
 INTEGER,PARAMETER   :: index_pressure   = 1
 INTEGER,PARAMETER   :: index_scalar     = 2
 INTEGER,PARAMETER   :: index_potential  = 3
 REAL(KIND=CGREAL)   :: restolLaplace
 REAL(KIND=CGREAL)   :: restolPotential
 INTEGER             :: it_solver_type_laplace
 INTEGER             :: Xbcx1, Xbcx2, Xbcy1, Xbcy2, Xbcz1, Xbcz2 
 INTEGER             :: Phibcx1, Phibcx2, Phibcy1, Phibcy2, Phibcz1, Phibcz2 
 REAL(KIND=CGREAL)   :: xxxx1,xxxx2, xxxy1,xxxy2, xxxz1,xxxz2
 REAL(KIND=CGREAL)   :: phix1,phix2, phiy1,phiy2, phiz1,phiz2
 REAL(KIND=CGREAL)   :: fx1, fx2, fy1, fy2, fz1, fz2 
 REAL(KIND=CGREAL)   :: botsurf_kinetic_energy, topsurf_kinetic_energy
 REAL(KIND=CGREAL)   :: kinetic_energy_difference
 REAL(KIND=CGREAL)   :: botsurf_pressure, topsurf_pressure
 REAL(KIND=CGREAL)   :: pressure_difference
 REAL(KIND=CGREAL)   :: outer_shearforce
 REAL(KIND=CGREAL)   :: Grad_u_square_dot_Grad_Dx_total
 REAL(KIND=CGREAL)   :: outer_Dx_dot_advection, outer_Dx_dot_gradu2
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: scalarX
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: potential
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Xgradx, Xgrady, Xgradz
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Potential_u, Potential_v, Potential_w
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Vorticity_u, Vorticity_v, Vorticity_w
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Xgradx1,Xgradx2
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Xgrady1,Xgrady2 
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Xgradz1,Xgradz2 
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Phigradx1,Phigradx2
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Phigrady1,Phigrady2 
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Phigradz1,Phigradz2 
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: boundScalarSource
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: boundPhiSource
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: zero
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: XGhost
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: PhiGhost 
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: omega_cross_u
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: XgradxBodyMarker, XgradyBodyMarker, XgradzBodyMarker
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: X_BodyMarker, X_BodyMarker_prev
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: advection_x_BodyMarker, advection_y_BodyMarker, advection_z_BodyMarker
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Lamb_x_BodyMarker, Lamb_y_BodyMarker, Lamb_z_BodyMarker
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: unsteady_term1, unsteady_term2, tau_Xgrad 
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: unsteady1, unsteady1_prev
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: unsteady_sum1, unsteady_sum2, unsteady_sum2p, unsteady_sum2v, advection_sum, Lamb_sum
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: X_Marker 
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: ssBodyMarker
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: dssBodyMarker 
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: x_force, y_force, z_force
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uTriElem, vTriElem, wTriElem
! REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: pTriElemCent0, pTriElemCent1
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: u_n, u_n_prev, u_n1_prev, u_n2_prev, u_n_dot, undot, udotn
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: advection, advection_prev, AB2_advection
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: F_i
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Advection_force, Unsteady_force, Shear_force, Vortex_force
 REAL(KIND=CGREAL) :: F_i_total, Advection_force_total, Unsteady_force_total, Shear_force_total, Vortex_force_total, Div_AB2_Xx_force_total, Div_Gradu2_Xx_force_total
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: diff_prev, diff_star, gradP, gradu2
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: PuPt, PvPt, PwPt
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyCent_prev, vBodyCent_prev, wBodyCent_prev  
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyCent_p_prev, vBodyCent_p_prev, wBodyCent_p_prev
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyCent_v_prev, vBodyCent_v_prev, wBodyCent_v_prev  
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: adpwl, sspwl, ppwl, totalpwl
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: adpw_l, sspw_l, ppw_l, totalpw_l
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: u_dot_n, up_dot_n, uv_dot_n
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: advection_n, Lamb_normal
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Div_Lamb
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: u_square, u_squareGhost
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: potential_u_square, vorticity_u_square
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: Grad_u_square
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: Grad_potential_u_square, Grad_vorticity_u_square
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_u_square_x, Grad_u_square_y, Grad_u_square_z
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_potential_u_square_x, Grad_potential_u_square_y, Grad_potential_u_square_z
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_vorticity_u_square_x, Grad_vorticity_u_square_y, Grad_vorticity_u_square_z
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_u_square_n
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_potential_u_square_n, Grad_vorticity_u_square_n
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: Grad_u_square_sum, Grad_potential_u_square_sum, Grad_vorticity_u_square_sum 
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: AB2_advection_Xx, Gradu2_Xx
 REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: Div_AB2_Xx, Div_Gradu2_Xx
 REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: Div_AB2_Xx_force, Div_Gradu2_Xx_force
 REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: Grad_u_square_dot_Grad_Dx
 REAL(KIND=CGREAL) :: outer_unsteady
 REAL(KIND=CGREAL) :: outer_shearforce_down
 REAL(KIND=CGREAL) :: outer_shearforce_up
 REAL(KIND=CGREAL) :: outer_shearforce_front
 REAL(KIND=CGREAL) :: outer_shearforce_back
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: unsteady_sum_check 
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: angle

ENDMODULE scalar_phi
