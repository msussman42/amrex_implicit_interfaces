! --------------------------------------------------------------------
!  Flow Physics and Computation Lab
!   
!
!  FSI module for
!  VICAR3D, a viscous, Cartesian, 3D flow solver.
!
!  This is a contineously developing project.
!
!  Starting Developers:
!  Kourosh Shoele
!
!
!  Final Filename: FiberShell26_20160629.F
!  Latest Modification: July, 28 2016 
!  by Kourosh Shoele
! --------------------------------------------------------------------

!@T
!\section{Introduction}
! This subroutine provides FSI module for {\bf Vicar3D} solver\\
! This version is written and modified by Kourosh Shoele 
!@q

!@T
!\section{Subroutines}
!@c
      MODULE HeaderFSI
      implicit none 
      include 'grid_def'
!@T
! Input parameters and runtime variables. Most variables are allocatable
! with the corresponding sizes are defined in [[grid_def]] file. \\
! The code needs to recompile as [[grid_def]] is modified.
!@q
      REAL*8  DELTA_T
     &       ,gx(2000)
     &       ,gy(2000)
     &       ,gz(2000)
     &       ,TIME,TIME2
     &       ,PI1
     &       ,EPS

      INTEGER TIME_STEP
     &       ,ndim
     &       ,nxgrid_ibm
     &       ,nygrid_ibm
     &       ,nzgrid_ibm
     &       ,Imp_vol_flag
     &       ,t0_restartfile
     &       ,dt_restartfile

      LOGICAL VARIABLE_DT
     &       ,FIRST_TIME
     &       ,CREATE_NEW_FLOW
     &       ,FLAG_3D
     &       ,PERIODIC_FLAG
     &       ,flag_restartfile
 
! IBM parameters
      real*8 ::
     &        MIN_GRID_X
     &       ,MIN_GRID_Y
     &       ,MIN_GRID_Z
     &       ,KS_IBM
     &       ,KB_IBM
     &       ,K_MASSIVE_IBM
     &       ,k_link	
     &       ,Kx_SPRING
     &       ,Ky_SPRING
     &       ,Ut0
     &       ,Vt0
     &       ,alpha_IBM
     &       ,beta_IBM
     &       ,alphabetaRise
     &       ,Fr
     &       ,the_grav
     &       ,phi_grav                                               
     &       ,temporary_sum
     &       ,betaGlobal
     &       ,gammaGlobal
     &       ,deltaGlobal
     &       ,alphaGlobal 
     &       ,ImpactPen(3) 


      real*8, allocatable, dimension(:) ::
     &        Density_Coef    
     &       ,Density_coefP    
     &       ,time_release    
     &       ,Ks_IBM_r    
     &       ,Kb_IBM_r    
     &       ,K_MASSIVE_IBM_r    
     &       ,cs_ibm    
     &       ,cs_ibm_target    
     &       ,Phi_ben_coef    
     &       ,Phi_nus_coef    
     &       ,wei    
     &       ,wei0    
     &       ,weiCont    
     &       ,wei_h    
     &       ,gau_h  
     &       ,Ccont    
     &       ,rot_seq
     &       ,alpham_Coef  
     &       ,alphaf_Coef  
     &       ,beta_Coef    
     &       ,gamma_Coef   


      real*8, allocatable, dimension(:,:) ::
     &        Gt0    
     &       ,Tramps    
     &       ,Trampx    
     &       ,ampx    
     &       ,freqx    
     &       ,phix    
     &       ,Trampt    
     &       ,ampt    
     &       ,freqt    
     &       ,phit    
     &       ,propDamp 

      real*8, allocatable, dimension(:,:) ::
     &        gau0    

      real*8, allocatable, dimension(:,:) ::
     &        Area    
     &       ,Phi_mem_coef    
     &       ,gau    
     &       ,gauCont    


      real*8, allocatable, dimension(:,:) ::
     &        contact_time_fsh 
     &       ,contact_time_esh 
     &       ,contact_time_fbc 


      real*8, allocatable, dimension(:,:) ::
     &        GX_IBM_fib    
     &       ,GY_IBM_fib    
     &       ,GZ_IBM_fib    
     &       ,GX_IBM_MASSIVE_fib    
     &       ,GY_IBM_MASSIVE_fib    
     &       ,GZ_IBM_MASSIVE_fib    
     &       ,GX_IBM_MASSIVEo_fib    
     &       ,GY_IBM_MASSIVEo_fib    
     &       ,GZ_IBM_MASSIVEo_fib    
     &       ,GX_IBMo1_fib    
     &       ,GY_IBMo1_fib    
     &       ,GZ_IBMo1_fib    
     &       ,GX_IBMpre_fib    
     &       ,GY_IBMpre_fib    
     &       ,GZ_IBMpre_fib    
     &       ,GX_IBM_MASSIVEpre_fib    
     &       ,GY_IBM_MASSIVEpre_fib    
     &       ,GZ_IBM_MASSIVEpre_fib    
     &       ,DSF_IBM_fib    
     &       ,DS_IBM_fib    
     &       ,FIBM1_fib    
     &       ,FIBM2_fib    
     &       ,FIBM3_fib    
     &       ,FS_1_IBM_fib    
     &       ,FS_2_IBM_fib    
     &       ,FS_1_IBMo_fib    
     &       ,FB_1_IBM_fib    
     &       ,FB_2_IBM_fib    
     &       ,FB_3_IBM_fib    
     &       ,FIner_1_IBM_fib    
     &       ,FIner_2_IBM_fib    
     &       ,FIner_3_IBM_fib    
     &       ,F_LINK1_fib    
     &       ,F_LINK2_fib    
     &       ,F_LINK3_fib    
     &       ,F_impuls1_fib    
     &       ,F_impuls2_fib    
     &       ,F_impuls3_fib    
     &       ,RIBM_fib    
     &       ,VIBM1_fib    
     &       ,VIBM2_fib    
     &       ,VIBM3_fib    
     &       ,VIBM1_pre_fib    
     &       ,VIBM2_pre_fib    
     &       ,VIBM3_pre_fib    
     &       ,MASS_IBM_fib    
     &       ,Q_MASSo_fib    
     &       ,Q_MASS_fib    
     &       ,GX_BP_fib    
     &       ,GY_BP_fib    
     &       ,GZ_BP_fib    
     &       ,GX_BP_fib0    
     &       ,GY_BP_fib0    
     &       ,GZ_BP_fib0    


      real*8, allocatable, dimension(:,:,:) ::
     &        GX_IBM_fsh    
     &       ,GY_IBM_fsh    
     &       ,GZ_IBM_fsh    
     &       ,GX_IBM_MASSIVE_fsh    
     &       ,GY_IBM_MASSIVE_fsh    
     &       ,GZ_IBM_MASSIVE_fsh    
     &       ,GX_IBM_MASSIVEo_fsh    
     &       ,GY_IBM_MASSIVEo_fsh    
     &       ,GZ_IBM_MASSIVEo_fsh    
     &       ,GX_IBMo1_fsh    
     &       ,GY_IBMo1_fsh    
     &       ,GZ_IBMo1_fsh    
     &       ,GX_IBMpre_fsh    
     &       ,GY_IBMpre_fsh    
     &       ,GZ_IBMpre_fsh    
     &       ,GX_IBM_MASSIVEpre_fsh    
     &       ,GY_IBM_MASSIVEpre_fsh    
     &       ,GZ_IBM_MASSIVEpre_fsh    
     &       ,DSF_IBM_fsh    
     &       ,DS_IBM_fsh    
     &       ,DSF2_IBM_fsh    
     &       ,DS2_IBM_fsh    
     &       ,DKF_IBM_fsh    
     &       ,DK_IBM_fsh    
     &       ,DKF2_IBM_fsh    
     &       ,DK2_IBM_fsh    
     &       ,FIBM1_fsh    
     &       ,FIBM2_fsh    
     &       ,FIBM3_fsh    
     &       ,FS_1_IBM_fsh    
     &       ,FS_2_IBM_fsh    
     &       ,FS_1_IBMo_fsh    
     &       ,FB_1_IBM_fsh    
     &       ,FB_2_IBM_fsh    
     &       ,FB_3_IBM_fsh    
     &       ,FIner_1_IBM_fsh    
     &       ,FIner_2_IBM_fsh    
     &       ,FIner_3_IBM_fsh    
     &       ,F_LINK1_fsh    
     &       ,F_LINK2_fsh    
     &       ,F_LINK3_fsh    
     &       ,F_impuls1_fsh    
     &       ,F_impuls2_fsh    
     &       ,F_impuls3_fsh    
     &       ,F_impulsHis1_fsh    
     &       ,F_impulsHis2_fsh    
     &       ,F_impulsHis3_fsh    
     &       ,RIBM_fsh    
     &       ,VIBM1_fsh    
     &       ,VIBM2_fsh    
     &       ,VIBM3_fsh    
     &       ,AIBM1_fsh    
     &       ,AIBM2_fsh    
     &       ,AIBM3_fsh    
     &       ,VIBM1_pre_fsh    
     &       ,VIBM2_pre_fsh    
     &       ,VIBM3_pre_fsh    
     &       ,MASS_IBM_fsh    
     &       ,Q_MASSo_fsh    
     &       ,Q_MASS_fsh    
     &       ,GX_BP_fsh    
     &       ,GY_BP_fsh    
     &       ,GZ_BP_fsh    
     &       ,GX_BP_fsh0    
     &       ,GY_BP_fsh0    
     &       ,GZ_BP_fsh0    



      real*8, allocatable, dimension(:,:) ::
     &        GX_IBM_esh    
     &       ,GY_IBM_esh    
     &       ,GZ_IBM_esh    
     &       ,GX_IBM_MASSIVE_esh    
     &       ,GY_IBM_MASSIVE_esh    
     &       ,GZ_IBM_MASSIVE_esh    
     &       ,GX_IBM_MASSIVEo_esh    
     &       ,GY_IBM_MASSIVEo_esh    
     &       ,GZ_IBM_MASSIVEo_esh    
     &       ,GX_IBMo1_esh    
     &       ,GY_IBMo1_esh    
     &       ,GZ_IBMo1_esh    
     &       ,GX_IBMpre_esh    
     &       ,GY_IBMpre_esh    
     &       ,GZ_IBMpre_esh    
     &       ,GX_IBM_MASSIVEpre_esh    
     &       ,GY_IBM_MASSIVEpre_esh    
     &       ,GZ_IBM_MASSIVEpre_esh    
     &       ,DS_IBM_esh    
     &       ,FIBM1_esh    
     &       ,FIBM2_esh    
     &       ,FIBM3_esh    
     &       ,FS_1_IBM_esh    
     &       ,FS_2_IBM_esh    
     &       ,FS_1_IBMo_esh    
     &       ,FB_1_IBM_esh    
     &       ,FB_2_IBM_esh    
     &       ,FB_3_IBM_esh    
     &       ,FIner_1_IBM_esh    
     &       ,FIner_2_IBM_esh    
     &       ,FIner_3_IBM_esh    
     &       ,F_LINK1_esh    
     &       ,F_LINK2_esh    
     &       ,F_LINK3_esh    
     &       ,F_impuls1_esh    
     &       ,F_impuls2_esh    
     &       ,F_impuls3_esh    
     &       ,RIBM_esh    
     &       ,VIBM1_esh    
     &       ,VIBM2_esh    
     &       ,VIBM3_esh    
     &       ,AIBM1_esh    
     &       ,AIBM2_esh    
     &       ,AIBM3_esh    
     &       ,VIBM1_pre_esh    
     &       ,VIBM2_pre_esh    
     &       ,VIBM3_pre_esh    
     &       ,MASS_IBM_esh    
     &       ,Q_MASSo_esh    
     &       ,Q_MASS_esh    
     &       ,GX_BP_esh    
     &       ,GY_BP_esh    
     &       ,GZ_BP_esh    
     &       ,GX_BP_esh0    
     &       ,GY_BP_esh0    
     &       ,GZ_BP_esh0    
     &       ,target_kvalue_esh 


      real*8, allocatable, dimension(:,:,:,:,:,:) ::
     &        gbase_con0SAVE_esh  
     &       ,gmetric_con0SAVE_esh  
      real*8, allocatable, dimension(:,:,:,:) ::
     &        detC_inplane0SAVE_esh 

      real*8, allocatable, dimension(:,:,:) ::
     &        Nn    
     &       ,dNndv    
     &       ,dNndw    
     &       ,ddNndvdv    
     &       ,ddNndwdw    
     &       ,ddNndvdw    
     &       ,contact_coef_esh 

      real*8, allocatable, dimension(:,:,:) ::
     &        NnCont    
     &       ,dNndvCont    
     &       ,dNndwCont    
     &       ,ddNndvdvCont    
     &       ,ddNndwdwCont    
     &       ,ddNndvdwCont    

      real*8, allocatable, dimension(:,:,:) ::
     &        Nn0    
     &       ,dNndv0    
     &       ,dNndw0    
     &       ,ddNndvdv0    
     &       ,ddNndwdw0    
     &       ,ddNndvdw0    


      real*8, allocatable, dimension(:,:,:) ::
     &        UIBM1_fib    
     &       ,UIBM2_fib    
     &       ,UIBM3_fib    
     &       ,FK_MASS1_fib    
     &       ,FK_MASS2_fib    
     &       ,FK_MASS3_fib    
     &       ,FK_MASS1o_fib    
     &       ,FK_MASS2o_fib    
     &       ,FK_MASS3o_fib    
     &       ,ffluidsum1_fib    
     &       ,ffluidsum2_fib    
     &       ,ffluidsum3_fib    


      real*8, allocatable, dimension(:,:,:) ::
     &        UIBM1_fsh    
     &       ,UIBM2_fsh    
     &       ,UIBM3_fsh    
     &       ,FK_MASS1_fsh    
     &       ,FK_MASS2_fsh    
     &       ,FK_MASS3_fsh    
     &       ,FK_MASS1o_fsh    
     &       ,FK_MASS2o_fsh    
     &       ,FK_MASS3o_fsh    
     &       ,ffluidsum1_fsh    
     &       ,ffluidsum2_fsh    
     &       ,ffluidsum3_fsh    
      real*8, allocatable, dimension(:,:,:,:) ::
     &        Mem_Coef_fsh    
     &       ,Ben_Coef_fsh    
     &       ,contact_coef_fsh 

      real*8, allocatable, dimension(:,:,:,:,:) ::
     &        Tzero_fsh    
     &       ,Bzero_fsh    



      real*8, allocatable, dimension(:,:) ::
     &        UIBM1_esh    
     &       ,UIBM2_esh    
     &       ,UIBM3_esh    
     &       ,FK_MASS1_esh    
     &       ,FK_MASS2_esh    
     &       ,FK_MASS3_esh    
     &       ,FK_MASS1o_esh    
     &       ,FK_MASS2o_esh    
     &       ,FK_MASS3o_esh    
     &       ,ffluidsum1_esh    
     &       ,ffluidsum2_esh    
     &       ,ffluidsum3_esh    

      real*8, allocatable, dimension(:,:,:,:) ::
     &        coorCurv_esh    

      real*8, allocatable, dimension(:,:,:) ::
     &        Jreal_2_imag_esh    
     &       ,crFibCrv_esh    
     &       ,xpcenter_esh    

      real*8, allocatable, dimension(:,:,:,:,:) ::
     &        Dmat0_esh    
     &       ,Kmat0_esh    



      integer ::   
     &        I_FIX  
     &       ,I_FIX2
     &       ,SPRING_START  
     &       ,SPRING_END    
     &       ,SAVE_FORCE_INT  
     &       ,Del_FORCE  
     &       ,CHECK_INT  
     &       ,Nr_IBM_f  
     &       ,Nr_IBM_fb
     &       ,n_jTypeAll
     &       ,SimulationType
     &       ,ngau
     &       ,ngau0
     &       ,ngauCont
     &       ,ngau_h
     &       ,contacttype
     &       ,s_START
     &       ,s_END
     &       ,imaster
     &       ,Nr_IBM_f_fib  
     &       ,Nr_IBM_fb_fib
     &       ,Nr_IBM_f_fsh  
     &       ,Nr_IBM_fb_fsh 
     &       ,Nr_IBM_f_esh   
     &       ,Nr_IBM_fb_esh 



      integer, allocatable, dimension(:) ::
     &        left_bndy    
     &       ,right_bndy    
     &       ,bottom_bndy    
     &       ,top_bndy    
     &       ,Nq_IBM_r    
     &       ,Ns_IBM_r    
     &       ,NS_IBM_i    
     &       ,Ns_IBM_rall    
     &       ,FiberBdry    
     &       ,BodyType    
     &       ,numelr    
     &       ,numeli    
     &       ,jTypeAll    
     &       , MasterBdy    
     &       ,iflaginext    
     &       ,BendingZero    
     &       ,MaterialTypeIBM    
     &       ,Nq_IBM_r_fib    
     &       ,Ns_IBM_r_fib    
     &       ,NS_IBM_i_fib    
     &       ,Ns_IBM_rall_fib    
     &       ,Nq_IBM_r_fsh    
     &       ,Ns_IBM_r_fsh    
     &       ,NS_IBM_i_fsh    
     &       ,Ns_IBM_rall_fsh    
     &       ,Nq_IBM_r_esh    
     &       ,Ns_IBM_r_esh    
     &       ,NS_IBM_i_esh    
     &       ,Ns_IBM_rall_esh    
     &       ,Iglbloc_fib  
     &       ,Iglbloc_esh  
     &       ,Iglbloc_fsh  
     &       ,Ilocglb_fib  
     &       ,Ilocglb_esh  
     &       ,Ilocglb_fsh  
     &       ,nMem_Coef_fsh
     &       ,n_matpara_esh
     &       ,Fibrous_esh  
     &       ,comprerssibleflag_esh  
     &       ,ShellModelType_esh  
     &       ,Genalpha_timesolver  
     &       ,Genalpha_niter 
     &       ,Piezo_Formulation_Flag

      integer, allocatable, dimension(:,:) ::
     &        BndyCnd_esh    
     &       ,neaType_esh    
     &       ,npos_esh    
     &       ,npos_fib    
     &       ,target_ktype_esh

      integer, allocatable, dimension(:,:,:) ::
     &        BoundryList_esh    
     &       ,nea_esh    
     &       ,nposele_esh    
     &       ,ContactShellFlag_fsh    
     &       ,ContactShellFlag_esh    
     &       ,BoundaryFibShell_fsh    
     &       ,npos_fsh    

      LOGICAL
     &        CLOSE_BD_flage
     &       ,para_coor_flag
     &       ,MASSIVE_IBM
     &       ,CREATE_IBM_flag
     &       ,STRU_SIMPLIFY_FLAG
     &       ,channelcontacty
     &       ,channelcontactz


      LOGICAL, allocatable, dimension(:) ::
     &        FluidForceFlag    
     &       ,ThermalForceFlag    
     &       ,ContactForceFlag    
     &       ,ExternalPressureFlag    
     &       ,Electromechanical
     &       ,ImplicitElecMech
     &       ,Piezo_Dist_Flag

C IBM parameters
      real*8  
     &        del_S
     &       ,alphaTEMP_IBM
     &       ,betaTEMP_IBM


      real*8, allocatable, dimension(:) ::
     &        TSIBMin    

      real*8, allocatable, dimension(:,:) ::
     &        QIBM_fib    
     &       ,TFIBM_fib    
     &       ,TSIBM_fib    
      real*8, allocatable, dimension(:,:,:) ::
     &        QIBM_fsh    
     &       ,TFIBM_fsh    
     &       ,TSIBM_fsh    
      real*8, allocatable, dimension(:,:) ::
     &        QIBM_esh    
     &       ,TFIBM_esh    
     &       ,TSIBM_esh    

C IBM parameters - Delta_type/width
      INTEGER 
     &        Del_X
     &       ,Del_Y
     &       ,Del_Z
     &       ,Delta_typeX
     &       ,Delta_typeY
     &       ,Delta_typeZ

C Target Points
       real*8  
     &        AROT_fin  
     &       ,AROT0_fin  
     &       ,maxErrorstruct  
     &       ,maxErrorfluid
     &       ,Solverabstol,Solverreltol

      INTEGER 
     &        SolverOutiter,SolverIniter,SolverPre,SolverILU

      real*8, allocatable, dimension(:,:) ::
     &        frequency_fin    

      real*8, allocatable, dimension(:,:) ::
     &        Target_K_link_fib    
     &       ,Target_t_link_fib    
     &       ,Target_K_link_fsh    
     &       ,Target_t_link_fsh    
     &       ,Target_K_link_esh    
     &       ,Target_t_link_esh    

      real*8, allocatable, dimension(:,:,:) ::
     &        phi_fin_fib    
     &       ,a_fin_fib    
     &       ,phi_fin_fsh    
     &       ,a_fin_fsh    
     &       ,phi_fin_esh    
     &       ,a_fin_esh    



      real*8, allocatable, dimension(:,:,:) ::
     &        Target_points_fib    
     &       ,force_points_fib    
     &       ,Target_points_v_fib    
     &       ,Target_points_fsh    
     &       ,force_points_fsh    
     &       ,Target_points_v_fsh    
     &       ,Target_points_esh    
     &       ,force_points_esh    
     &       ,Target_points_v_esh    


      real*8, allocatable, dimension(:,:,:) ::
     &       Materialpara_esh    

      integer maxiter,maxnstruct,maxnfluid


      integer, allocatable, dimension(:) ::
     &         Target_num    
     &        ,implicitflag    

      integer, allocatable, dimension(:,:) ::
     &         FlagFixed_fib    
     &        ,FlagFixed_fsh    
     &        ,FlagFixed_esh    

      integer, allocatable, dimension(:,:,:) ::
     &         Target_point_num_fsh    
      integer, allocatable, dimension(:,:) ::
     &         Target_point_num_fib    
     &        ,Target_point_num_esh    


       logical k_flag

C Bucket Search
       integer 
     &         NBucketnum
     &        ,IBucketnum(3)

      integer, allocatable, dimension(:,:) ::
     &         NBucket    
     &        ,PBucket    

      integer, allocatable, dimension(:,:,:) ::
     &         IBucket_fsh    
      integer, allocatable, dimension(:,:) ::
     &         IBucket_esh    
     &        ,IBucket_fib    

      integer, allocatable, dimension(:,:,:,:) ::
     &         SBucket_fsh    
      integer, allocatable, dimension(:,:,:) ::
     &         SBucket_esh    
     &        ,SBucket_fib   

       real*8 
     &         Bucketntinterval
     &        ,timebucket
     &        ,BucketXstart(3)
     &        ,BucketXEnd(3)


       real*8 
     &         p_inflation
     &        ,BucketdX(3)

      real*8, allocatable, dimension(:) ::
     &         Pnormal_Ext    

       integer,allocatable, dimension (:) ::
     &    n_TDpt_Pnormal_Ext,flag_TDpt_Pnormal_Ext 

      real*8, allocatable, dimension(:,:,:) ::
     &         data_TDpt_Pnormal_Ext 

C piezo paramters
      real*8, allocatable, dimension(:) ::
     &        piezo_beta
     &       ,piezo_gamma,piezo_damp,piezo_coef
     &       ,piezo_alpha11
     &       ,piezo_alpha22
     &       ,piezo_alpha12    
      real*8, allocatable, dimension(:,:,:,:) ::
     &        piezo_Coef_fsh  
      real*8, allocatable, dimension(:,:,:) ::
     &        Fpiezo1_fsh    
     &       ,Fpiezo2_fsh    
     &       ,Fpiezo3_fsh    
     &       ,ndot_ibm_fsh    
     &       ,ndotold_ibm_fsh    
     &       ,ndotpre_ibm_fsh    
     &       ,ndot_ibm_fsh0
     &       ,vPiezo_fsh    
     &       ,dvPiezo_fsh    
     &       ,ddvPiezo_fsh    
     &       ,vPiezoold_fsh    
     &       ,vPiezopre_fsh    

      real*8, allocatable, dimension(:,:,:,:) ::
     &        TzeroB_fsh    

       real*8, allocatable, dimension(:):: 
     &         e4coef 

       integer,allocatable, dimension (:) ::
     &    nibmptrecord  

       integer,allocatable, dimension (:,:,:) ::
     &    ibmptrecord  

       integer, allocatable, dimension (:) ::
     &    Ncontactplanes   
       LOGICAL, allocatable, dimension(:) ::
     &    contactplanesFlag  
       real*8, allocatable, dimension(:,:,:) ::
     &    contactplanesPara    


      real*8, allocatable, dimension(:,:) ::
     &        GX_IBM_fbc    
     &       ,GY_IBM_fbc    
     &       ,GZ_IBM_fbc    
     &       ,GX_IBM_MASSIVE_fbc    
     &       ,GY_IBM_MASSIVE_fbc    
     &       ,GZ_IBM_MASSIVE_fbc    
     &       ,GX_IBM_MASSIVEo_fbc    
     &       ,GY_IBM_MASSIVEo_fbc    
     &       ,GZ_IBM_MASSIVEo_fbc    
     &       ,GX_IBMo1_fbc    
     &       ,GY_IBMo1_fbc    
     &       ,GZ_IBMo1_fbc    
     &       ,GX_IBMpre_fbc    
     &       ,GY_IBMpre_fbc    
     &       ,GZ_IBMpre_fbc    
     &       ,GX_IBM_MASSIVEpre_fbc    
     &       ,GY_IBM_MASSIVEpre_fbc    
     &       ,GZ_IBM_MASSIVEpre_fbc    
     &       ,DS_IBM_fbc    
     &       ,FIBM1_fbc    
     &       ,FIBM2_fbc    
     &       ,FIBM3_fbc    
     &       ,FS_1_IBM_fbc    
     &       ,FS_2_IBM_fbc    
     &       ,FS_1_IBMo_fbc    
     &       ,FB_1_IBM_fbc    
     &       ,FB_2_IBM_fbc    
     &       ,FB_3_IBM_fbc    
     &       ,FIner_1_IBM_fbc    
     &       ,FIner_2_IBM_fbc    
     &       ,FIner_3_IBM_fbc    
     &       ,F_LINK1_fbc    
     &       ,F_LINK2_fbc    
     &       ,F_LINK3_fbc    
     &       ,F_impuls1_fbc    
     &       ,F_impuls2_fbc    
     &       ,F_impuls3_fbc    
     &       ,RIBM_fbc    
     &       ,VIBM1_fbc    
     &       ,VIBM2_fbc    
     &       ,VIBM3_fbc    
     &       ,AIBM1_fbc    
     &       ,AIBM2_fbc    
     &       ,AIBM3_fbc    
     &       ,VIBM1_pre_fbc    
     &       ,VIBM2_pre_fbc    
     &       ,VIBM3_pre_fbc    
     &       ,MASS_IBM_fbc    
     &       ,Q_MASSo_fbc    
     &       ,Q_MASS_fbc    
     &       ,GX_BP_fbc    
     &       ,GY_BP_fbc    
     &       ,GZ_BP_fbc    
     &       ,GX_BP_fbc0    
     &       ,GY_BP_fbc0    
     &       ,GZ_BP_fbc0    
     &       ,target_kvalue_fbc 



      real*8, allocatable, dimension(:,:,:) ::
     &        contact_coef_fbc 
     &       ,phi_fin_fbc    
     &       ,a_fin_fbc    


      real*8, allocatable, dimension(:,:) ::
     &        UIBM1_fbc    
     &       ,UIBM2_fbc    
     &       ,UIBM3_fbc    
     &       ,FK_MASS1_fbc    
     &       ,FK_MASS2_fbc    
     &       ,FK_MASS3_fbc    
     &       ,FK_MASS1o_fbc    
     &       ,FK_MASS2o_fbc    
     &       ,FK_MASS3o_fbc    
     &       ,ffluidsum1_fbc    
     &       ,ffluidsum2_fbc    
     &       ,ffluidsum3_fbc    

      real*8, allocatable, dimension(:,:,:) ::
     &        Jreal_2_imag_fbc    
     &       ,xpcenter_fbc    

      integer ::   
     &        Nr_IBM_f_fbc   
     &       ,Nr_IBM_fb_fbc 

      integer, allocatable, dimension(:) ::
     &        Nq_IBM_r_fbc    
     &       ,Ns_IBM_r_fbc    
     &       ,NS_IBM_i_fbc    
     &       ,Ns_IBM_rall_fbc    
     &       ,Iglbloc_fbc  
     &       ,Ilocglb_fbc  
     &       ,n_matpara_fbc

      integer, allocatable, dimension(:,:) ::
     &        BndyCnd_fbc    
     &       ,target_ktype_fbc 

      integer, allocatable, dimension(:,:,:) ::
     &        BoundryList_fbc    
     &       ,nea_fbc    
     &       ,ContactShellFlag_fbc    

      real*8, allocatable, dimension(:,:) ::
     &        QIBM_fbc    
     &       ,TFIBM_fbc    
     &       ,TSIBM_fbc    

      real*8, allocatable, dimension(:,:) ::
     &        Target_K_link_fbc    
     &       ,Target_t_link_fbc    

      real*8, allocatable, dimension(:,:) ::
     &        phi_fnr_ibm_fbc    
     &       ,a_fnr_ibm_fbc    

      real*8, allocatable, dimension(:,:,:) ::
     &        Target_points_fbc    
     &       ,force_points_fbc    
     &       ,Target_points_v_fbc    

      integer, allocatable, dimension(:,:) ::
     &         FlagFixed_fbc    

      integer, allocatable, dimension(:,:) ::
     &         Target_point_num_fbc    

C Bucket Search
      integer, allocatable, dimension(:,:) ::
     &         IBucket_fbc    

      integer, allocatable, dimension(:,:,:) ::
     &         SBucket_fbc    


      integer, allocatable, dimension(:) ::
     &         numedge    
     &        ,numeFacePair    
     &        ,inplane_type_fbc

      integer, allocatable, dimension(:,:) ::
     &         InplaneLinemode_fbc    

      integer, allocatable, dimension(:,:,:) ::
     &         EdgeList_fbc    
     &        ,facepairEle_fbc    
     &        ,facepairNode_fbc    

      real*8, allocatable, dimension(:) ::
     &         kappa_volume_fbc    
     &        ,KBT_fbc    
     &        ,volume_eq_fbc    
     &        ,kappa_Area_fbc    
     &        ,Area_eq_fbc    
     &        ,edge_length_p_fbc    
     &        ,edge_length_max_fbc    
     &        ,edge_length_eq_fbc    
     &        ,InplaneLinem_fbc    
     &        ,kappa_Bend_fbc    
     &        ,vertex_mass_fbc    
     &        ,inplane_E_fbc
     &        ,inplane_nu_fbc

      real*8, allocatable, dimension(:,:) ::
     &         EdgeLength0_fbc    
     &        ,EleArea0_fbc    
     &        ,areapnt_fbc    
     &        ,facepairtet0_fbc    
     &        ,inplane_Cq_fbc    
     &        ,kappa_inplane_fbc    

      real*8, allocatable, dimension(:,:,:) ::
     &         EleCenter0_fbc    
     &        ,EleNormal0_fbc    
     &        ,EleEdge0_fbc      
     &        ,EleInplane_Cq_fbc 
     &        ,EleInplaneParaK_fbc   
     &        ,EleInplaneParaC_fbc   

      logical, allocatable, dimension(:) ::
     &         volume_flag_fbc    
     &        ,Aera_flag_fbc    
     &        ,inplane_flag_fbc    
     &        ,Bending_flag_fbc    
     &        ,tet0_fbc            


       LOGICAL, allocatable, dimension(:) ::
     &        FabricMeshFlag_fbc     
     &       ,FabricLineFlag_fbc     

       LOGICAL ::
     &        FabricMesh_Presence
     &       ,FabricLine_Presence


      real*8, allocatable, dimension(:,:) ::
     &        FabricMesh_target_k_link   
     &       ,FabricMesh_target_t_link   
     &       ,FabricMesh_t   
     &       ,FabricLine_target_k_link   
     &       ,FabricLine_target_t_link   
     &       ,FabricLine_t   


      real*8, allocatable, dimension(:,:,:) ::
     &        FabricMesh_coord        
     &       ,FabricMesh_coordpre        
     &       ,FabricMesh_coordo1        
     &       ,FabricMesh_coordMass        
     &       ,FabricMesh_coordMasso        
     &       ,FabricMesh_coordbp        
     &       ,FabricMesh_coordbp0        
     &       ,FabricMesh_force        
     &       ,FabricMesh_v        
     &       ,FabricMesh_a        
     &       ,FabricMesh_SurfElAttr        
     &       ,FabricMesh_attr        
     &       ,FabricMeshL_attr        
     &       ,FabricMesh_a_fin        
     &       ,FabricMesh_phi_fin        
     &       ,FabricMesh_target_points        
     &       ,FabricMesh_target_points_v        
     &       ,FabricLine_coord        
     &       ,FabricLine_coordpre        
     &       ,FabricLine_coordo1        
     &       ,FabricLine_coordMass        
     &       ,FabricLine_coordMasso        
     &       ,FabricLine_coordbp        
     &       ,FabricLine_coordbp0        
     &       ,FabricLine_force        
     &       ,FabricLine_v        
     &       ,FabricLine_a        
     &       ,FabricLine_SurfElAttr        
     &       ,FabricLine_attr        
     &       ,FabricLine_a_fin        
     &       ,FabricLine_phi_fin        
     &       ,FabricLine_target_points        
     &       ,FabricLine_target_points_v        
     &       ,FabricLine_curv0                  

      integer ::
     &        numMeshpointRecordMax
     &       ,numMeshLinkRecordMax
     &       ,numMeshMax
     &       ,numMeshpointMax
     &       ,numMeshLinkMax
     &       ,FabricMesh_numattMax
     &       ,FabricMeshL_numattMax
     &       ,FabricMesh_numtargetMax
     &       ,numLinepointRecordMax
     &       ,numLineMax
     &       ,numLinepointMax
     &       ,FabricLine_numattMax
     &       ,FabricLine_numtargetMax

      integer, allocatable, dimension(:) ::
     &        FabricMesh_numrecord        
     &       ,FabricMeshL_numrecord     
     &       ,FabricMesh_nMesh        
     &       ,FabricMesh_numatt        
     &       ,FabricMeshL_numatt        
     &       ,FabricMesh_target_num        
     &       ,FabricLine_numrecord        
     &       ,FabricLine_nLine        
     &       ,FabricLine_numatt        
     &       ,FabricLine_target_num        

      integer, allocatable, dimension(:,:) ::
     &        FabricMesh_npoint        
     &       ,FabricMesh_nLine        
     &       ,FabricMesh_SurfElAddress        
     &       ,FabricMesh_target_point_num        
     &       ,FabricMesh_flagfixed        
     &       ,FabricLine_npoint        
     &       ,FabricLine_SurfElAddress        
     &       ,FabricLine_target_point_num        
     &       ,FabricLine_flagfixed        


      integer, allocatable, dimension(:,:,:) ::
     &        FabricMesh_address        
     &       ,FabricMeshL_con        
     &       ,FabricMeshL_address        
     &       ,FabricMeshL_addressRev        
     &       ,FabricLine_address        
     &       ,FabricLine_addressRev        

      integer :: 
     &        iMm0,iMc0,iMconR0,iMcon0,iMs0
     &       ,iMpD0_1,iMpD0_2,iMpD0_3
     &       ,n_iMnpara0,iMnpara0(100)
     &       ,iMks0,iMl0,n_iMlpara0,iMlpara0(100)

      integer :: 
     &        iLm0,iLc0,iLconR0,iLcon0,iLs0,iLks0,iLkb0,iLl0
     &       ,iLpD0_1,iLpD0_2,iLpD0_3,n_iLnpara0,iLnpara0(100)



      integer nFixcontact,nFixcontactL
      real*8, allocatable, dimension(:,:) ::
     &         xFixcontact       

      integer*8, allocatable, dimension(:,:) ::
     &         EleFixcontact     


       integer, allocatable, dimension(:) ::
     &    nifractal_esh          

       LOGICAL, allocatable, dimension(:) ::
     &    flag_fractal_esh          
     &    ,fractalinitFlag_esh      

       real*8, allocatable, dimension(:,:) ::
     &     fractalpara_esh          
     &    ,Dfractal_esh             
     &    ,Jfractal2_esh          

       real*8, allocatable, dimension(:,:,:) ::
     &     Dfractal2_esh          
     &    ,Hfractal2_esh          
     &    ,Ifractal2_esh          
     &    ,Bfractal_esh           
     &    ,Cfractal_esh           
     &    ,fibmsig_esh           

       real*8, allocatable, dimension(:,:,:,:) ::
     &     Afractal_esh           
     &    ,Afractal2_esh          
     &    ,Bfractal2_esh          
     &    ,Cfractal2_esh          
     &    ,Xfractal_esh           

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      END MODULE
