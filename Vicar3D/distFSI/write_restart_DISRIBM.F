c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE write_restart_DISRIBM(temp_ibm2,myoutput,input_Flag)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c write restart files
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      integer myoutput
      integer input_Flag,i,j,ijtemp
      integer  ifuRstrtDistIBMOut
      integer  ifuRstrtDistIBMIn
      real*8  temp_ibm2
      CHARACTER*80 fname1

      write(*,*) 'temp_ibm2,myoutput,input_Flag'
      write(*,*)  temp_ibm2,myoutput,input_Flag 

      ifuRstrtDistIBMOut=1310
      ifuRstrtDistIBMin =1311
      if((input_Flag .ne. 1).and.(myoutput.ne.zero)) then
!             print*,3
      WRITE(fname1,"('restart_distIBM_out.dat')") 
      OPEN(ifuRstrtDistIBMOut,FILE=fname1,FORM='UNFORMATTED')

      WRITE(ifuRstrtDistIBMOut)
     &        DELTA_T
     &       ,gx(2000)
     &       ,gy(2000)
     &       ,gz(2000)
     &       ,TIME,TIME2
     &       ,PI1,EPS,TIME_STEP,ndim,nxgrid_ibm
     &       ,nygrid_ibm,nzgrid_ibm,Imp_vol_flag
     &       ,t0_restartfile,dt_restartfile,VARIABLE_DT
     &       ,FIRST_TIME,CREATE_NEW_FLOW,FLAG_3D,PERIODIC_FLAG
     &       ,flag_restartfile
     &       ,MIN_GRID_X,MIN_GRID_Y,MIN_GRID_Z,KS_IBM,KB_IBM
     &       ,K_MASSIVE_IBM,k_link,Kx_SPRING,Ky_SPRING,Ut0,Vt0
     &       ,alpha_IBM,beta_IBM,alphabetaRise,Fr,the_grav,phi_grav                                               
     &       ,temporary_sum,betaGlobal,gammaGlobal,deltaGlobal
     &       ,alphaGlobal ,ImpactPen(1:3),dsecy_IBM,dsecz_IBM
     &             
     &       ,Density_Coef(1:Nr_IBM)
     &       ,Density_coefP(1:Nr_IBM)
     &       ,time_release(1:Nr_IBM)
     &       ,Ks_IBM_r(1:Nr_IBM)
     &       ,Kb_IBM_r(1:Nr_IBM)
     &       ,K_MASSIVE_IBM_r(1:Nr_IBM)
     &       ,cs_ibm(1:Nr_IBM)
     &       ,cs_ibm_target(1:Nr_IBM)                               
     &       ,Phi_ben_coef(1:Nr_IBM) 
     &       ,Phi_nus_coef(1:Nr_IBM)
     &       ,wei(1:ngaumax)
     &       ,weiCont(1:ngaumax)
     &       ,wei_h(1:ngaumaxh),gau_h(1:ngaumaxh)
     &       ,wei0(1:ngaumax)
     &       ,Ccont(0:20)
     &       ,rot_seq(1:Nr_IBM)
     &       ,alpham_Coef(1:Nr_IBM)
     &       ,alphaf_Coef(1:Nr_IBM)
     &       ,beta_Coef(1:Nr_IBM)
     &       ,gamma_Coef(1:Nr_IBM)
     &
     &       ,Gt0(1:Nr_IBM,1:3)
     &       ,Tramps(1:Nr_IBM,1:3)
     &       ,Trampx(1:Nr_IBM,1:3)
     &       ,ampx(1:Nr_IBM,1:3)
     &       ,freqx(1:Nr_IBM,1:3)
     &       ,phix(1:Nr_IBM,1:3)
     &       ,Trampt(1:Nr_IBM,1:3)
     &       ,ampt(1:Nr_IBM,1:3)
     &       ,freqt(1:Nr_IBM,1:3)
     &       ,phit(1:Nr_IBM,1:3)
     &       ,propDamp(1:Nr_IBM,1:2)
     &       ,Phi_mem_coef(1:Nr_IBM,1:25) 
     &       ,gau(1:2,1:ngaumax)
     &       ,gauCont(1:2,1:ngaumax)
     &       ,gau0(1:2,1:ngaumax)
     &       
     &       ,Nn(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,dNndv(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,dNndw(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,ddNndvdv(1:jTypeMax,1:ngaumax,1:maxordern)  
     &       ,ddNndwdw(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,ddNndvdw(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       
     &       ,NnCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,dNndvCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,dNndwCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,ddNndvdvCont(1:jTypeMax,1:ngaumax,1:maxordern)  
     &       ,ddNndwdwCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,ddNndvdwCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       
     &       ,Nn0(1:jTypeMax,1,1:maxordern) 
     &       ,dNndv0(1:jTypeMax,1,1:maxordern) 
     &       ,dNndw0(1:jTypeMax,1,1:maxordern) 
     &       ,ddNndvdv0(1:jTypeMax,1,1:maxordern)  
     &       ,ddNndwdw0(1:jTypeMax,1,1:maxordern) 
     &       ,ddNndvdw0(1:jTypeMax,1,1:maxordern) 
     &       
     &       ,I_FIX  
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
     &       ,ngauCont
     &       ,ngau0
     &       ,ngau_h
     &       ,Nsec_IBM
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
     &       
     &       ,left_bndy(1:Nr_IBM)
     &       ,right_bndy(1:Nr_IBM)
     &       ,bottom_bndy(1:Nr_IBM)
     &       ,top_bndy(1:Nr_IBM) 
     &       ,Nq_IBM_r(1:Nr_IBM)
     &       ,Ns_IBM_r(1:Nr_IBM)  
     &       ,Ns_IBM_i(1:Nr_IBM)  
     &       ,Ns_IBM_rall(1:Nr_IBM)  
     &       ,FiberBdry(1:Nr_IBM)
     &       ,BodyType(1:Nr_IBM)
     &       ,numelr(1:Nr_IBM)  
     &       ,numeli(1:Nr_IBM)  
     &       ,jTypeAll(1:jTypeMax)
     &       ,MasterBdy(1:Nr_IBM) 
     &       ,iflaginext(1:Nr_IBM)
     &       ,BendingZero(1:Nr_IBM)
     &       ,MaterialTypeIBM(1:Nr_IBM)
     &
     &       ,Iglbloc_fib(1:Nr_IBM)
     &       ,Iglbloc_esh(1:Nr_IBM)
     &       ,Iglbloc_fsh(1:Nr_IBM) 
     &       ,Ilocglb_fib(1:Nr_IBM)
     &       ,Ilocglb_esh(1:Nr_IBM)
     &       ,Ilocglb_fsh(1:Nr_IBM) 
     &       
     &       ,CLOSE_BD_flage
     &       ,para_coor_flag
     &       ,MASSIVE_IBM
     &       ,CREATE_IBM_flag
     &       ,STRU_SIMPLIFY_FLAG
     &       ,channelcontacty
     &       ,channelcontactz
     &       
     &       ,FluidForceFlag(1:Nr_IBM)
     &       ,ThermalForceFlag(1:Nr_IBM)
     &       ,ContactForceFlag(1:Nr_IBM)
     &       ,ExternalPressureFlag(1:Nr_IBM)
     &       ,Electromechanical(1:Nr_IBM)
     &       ,ImplicitElecMech(1:Nr_IBM)
     &       ,Piezo_Dist_Flag(1:Nr_IBM)
     &       
     &       ,del_S
     &       ,alphaTEMP_IBM
     &       ,betaTEMP_IBM
     &       
     &       ,TSIBMin(1:Nr_IBM)
     &       
     &       ,Del_X
     &       ,Del_Y
     &       ,Del_Z
     &       ,Delta_typeX
     &       ,Delta_typeY
     &       ,Delta_typeZ
     &       
     &       ,AROT_fin  
     &       ,AROT0_fin  
     &       ,maxErrorstruct  
     &       ,maxErrorfluid
     &       ,Solverabstol,Solverreltol,SolverOutiter
     &       ,SolverIniter,SolverPre, SolverILU
     &       ,frequency_fin(1:Nr_IBM,1:3)
     &       
     &       ,maxiter,maxnstruct,maxnfluid
     &       
     &       ,Target_num(1:Nr_IBM)
     &       ,implicitflag(1:Nr_IBM)
     &       
     &       ,k_flag
     &       
     &       ,NBucketnum
     &       ,IBucketnum(1:3)
     &       
     &       ,NBucket(1:Nr_IBM,1:NMaxnumbuck)
     &       ,PBucket(1:Nr_IBM,1:NMaxnumbuck)
     &       
     &       ,Bucketntinterval
     &       ,timebucket
     &       ,BucketXstart(1:3)
     &       ,BucketXEnd(1:3)
     &       
     &       ,p_inflation
     &       ,BucketdX(1:3)
     &       
     &       ,Pnormal_Ext(1:Nr_IBM)
     &       ,n_TDpt_Pnormal_Ext(1:Nr_IBM)
     &       ,flag_TDpt_Pnormal_Ext(1:Nr_IBM)
     &       ,data_TDpt_Pnormal_Ext(1:Nr_IBM,1:1000,1:2)
     &       
     &       ,piezo_beta(1:Nr_IBM)
     &       ,piezo_gamma(1:Nr_IBM)
     &       ,piezo_damp(1:Nr_IBM)
     &       ,piezo_coef(1:Nr_IBM)
     &       ,piezo_alpha11(1:Nr_IBM)
     &       ,piezo_alpha22(1:Nr_IBM)
     &       ,piezo_alpha12(1:Nr_IBM)
     &       ,e4coef(1:Nr_IBM)

         WRITE(ifuRstrtDistIBMOut)
     &       nibmptrecord(1:Nr_IBM)
     &      ,ibmptrecord(1:Nr_IBM,1:100,1:2)

         WRITE(ifuRstrtDistIBMOut)
     &        Genalpha_timesolver(1:Nr_IBM)
     &       ,Genalpha_niter(1:Nr_IBM)
     &       ,Piezo_Formulation_Flag(1:Nr_IBM)

       WRITE(ifuRstrtDistIBMOut)      
     &        contactplanesFlag(1:Nr_IBM)
     &       ,Ncontactplanes(1:Nr_IBM)
     &       ,contactplanesPara(1:Nr_IBM,1:NPlanesIBM,1:6)


      WRITE(ifuRstrtDistIBMOut)
     &        Nr_IBM_f_fbc   
     &       ,Nr_IBM_fb_fbc 
     &       ,Iglbloc_fbc(1:Nr_IBM)
     &       ,Ilocglb_fbc(1:Nr_IBM)

      if(nr_IBM_fib>=1) 
     &           call write_restart_DISRIBM_fib(ifuRstrtDistIBMOut)
      if(nr_IBM_fsh>=1) 
     &           call write_restart_DISRIBM_fsh(ifuRstrtDistIBMOut)
      if(nr_IBM_esh>=1) 
     &           call write_restart_DISRIBM_esh(ifuRstrtDistIBMOut)
      if(nr_IBM_fbc>=1) 
     &           call write_restart_DISRIBM_fbc(ifuRstrtDistIBMOut)

      WRITE(ifuRstrtDistIBMOut) nFixcontact,nFixcontactL
         IF(nFixcontact>0)
     &        WRITE(ifuRstrtDistIBMOut)  
     &        xFixcontact(1:nFixcontact ,1:3)
         IF(nFixcontactL>0)
     &        WRITE(ifuRstrtDistIBMOut)  
     &        ,EleFixcontact(1:nFixcontactL,1:3)
        CLOSE(ifuRstrtDistIBMout) 

      elseif(input_Flag .eq. 1) then

        WRITE(fname1,"('restart_distIBM_in.dat')") 
      !    print*,'K=',1
        OPEN(ifuRstrtDistIBMin,FILE=fname1,FORM='UNFORMATTED')
       !  print*,'K=',2
      READ (ifuRstrtDistIBMin)
     &        DELTA_T
     &       ,gx(2000)
     &       ,gy(2000)
     &       ,gz(2000)
     &       ,TIME,TIME2
     &       ,PI1,EPS,TIME_STEP,ndim,nxgrid_ibm
     &       ,nygrid_ibm,nzgrid_ibm,Imp_vol_flag
     &       ,t0_restartfile,dt_restartfile,VARIABLE_DT
     &       ,FIRST_TIME,CREATE_NEW_FLOW,FLAG_3D,PERIODIC_FLAG
     &       ,flag_restartfile
     &       ,MIN_GRID_X,MIN_GRID_Y,MIN_GRID_Z,KS_IBM,KB_IBM
     &       ,K_MASSIVE_IBM,k_link,Kx_SPRING,Ky_SPRING,Ut0,Vt0
     &       ,alpha_IBM,beta_IBM,alphabetaRise,Fr,the_grav,phi_grav                                               
     &       ,temporary_sum,betaGlobal,gammaGlobal,deltaGlobal
     &       ,alphaGlobal ,ImpactPen(1:3),dsecy_IBM,dsecz_IBM
     &       
     &       ,Density_Coef(1:Nr_IBM)
     &       ,Density_coefP(1:Nr_IBM)
     &       ,time_release(1:Nr_IBM)
     &       ,Ks_IBM_r(1:Nr_IBM)
     &       ,Kb_IBM_r(1:Nr_IBM)
     &       ,K_MASSIVE_IBM_r(1:Nr_IBM)
     &       ,cs_ibm(1:Nr_IBM)
     &       ,cs_ibm_target(1:Nr_IBM)                               
     &       ,Phi_ben_coef(1:Nr_IBM) 
     &       ,Phi_nus_coef(1:Nr_IBM)
     &       ,wei(1:ngaumax)
     &       ,weiCont(1:ngaumax)
     &       ,wei_h(1:ngaumaxh),gau_h(1:ngaumaxh)
     &       ,wei0(1:ngaumax)
     &       ,Ccont(0:20)
     &       ,rot_seq(1:Nr_IBM)
     &       ,alpham_Coef(1:Nr_IBM)
     &       ,alphaf_Coef(1:Nr_IBM)
     &       ,beta_Coef(1:Nr_IBM)
     &       ,gamma_Coef(1:Nr_IBM)
     &
     &       ,Gt0(1:Nr_IBM,1:3)
     &       ,Tramps(1:Nr_IBM,1:3)
     &       ,Trampx(1:Nr_IBM,1:3)
     &       ,ampx(1:Nr_IBM,1:3)
     &       ,freqx(1:Nr_IBM,1:3)
     &       ,phix(1:Nr_IBM,1:3)
     &       ,Trampt(1:Nr_IBM,1:3)
     &       ,ampt(1:Nr_IBM,1:3)
     &       ,freqt(1:Nr_IBM,1:3)
     &       ,phit(1:Nr_IBM,1:3)
     &       ,propDamp(1:Nr_IBM,1:2)
     &       ,Phi_mem_coef(1:Nr_IBM,1:25) 
     &       ,gau(1:2,1:ngaumax)
     &       ,gauCont(1:2,1:ngaumax)
     &       ,gau0(1:2,1:ngaumax)
     &       
     &       ,Nn(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,dNndv(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,dNndw(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,ddNndvdv(1:jTypeMax,1:ngaumax,1:maxordern)  
     &       ,ddNndwdw(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,ddNndvdw(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       
     &       ,NnCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,dNndvCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,dNndwCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,ddNndvdvCont(1:jTypeMax,1:ngaumax,1:maxordern)  
     &       ,ddNndwdwCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       ,ddNndvdwCont(1:jTypeMax,1:ngaumax,1:maxordern) 
     &       
     &       ,Nn0(1:jTypeMax,1,1:maxordern) 
     &       ,dNndv0(1:jTypeMax,1,1:maxordern) 
     &       ,dNndw0(1:jTypeMax,1,1:maxordern) 
     &       ,ddNndvdv0(1:jTypeMax,1,1:maxordern)  
     &       ,ddNndwdw0(1:jTypeMax,1,1:maxordern) 
     &       ,ddNndvdw0(1:jTypeMax,1,1:maxordern) 
     &       
     &       ,I_FIX  
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
     &       ,ngauCont
     &       ,ngau0
     &       ,ngau_h
     &       ,Nsec_IBM
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
     &       
     &       ,left_bndy(1:Nr_IBM)
     &       ,right_bndy(1:Nr_IBM)
     &       ,bottom_bndy(1:Nr_IBM)
     &       ,top_bndy(1:Nr_IBM) 
     &       ,Nq_IBM_r(1:Nr_IBM)
     &       ,Ns_IBM_r(1:Nr_IBM)  
     &       ,Ns_IBM_i(1:Nr_IBM)  
     &       ,Ns_IBM_rall(1:Nr_IBM)  
     &       ,FiberBdry(1:Nr_IBM)
     &       ,BodyType(1:Nr_IBM)
     &       ,numelr(1:Nr_IBM)  
     &       ,numeli(1:Nr_IBM)  
     &       ,jTypeAll(1:jTypeMax)
     &       ,MasterBdy(1:Nr_IBM) 
     &       ,iflaginext(1:Nr_IBM)
     &       ,BendingZero(1:Nr_IBM)
     &       ,MaterialTypeIBM(1:Nr_IBM)
     &
     &       ,Iglbloc_fib(1:Nr_IBM)
     &       ,Iglbloc_esh(1:Nr_IBM)
     &       ,Iglbloc_fsh(1:Nr_IBM) 
     &       ,Ilocglb_fib(1:Nr_IBM)
     &       ,Ilocglb_esh(1:Nr_IBM)
     &       ,Ilocglb_fsh(1:Nr_IBM) 
     &       
     &       ,CLOSE_BD_flage
     &       ,para_coor_flag
     &       ,MASSIVE_IBM
     &       ,CREATE_IBM_flag
     &       ,STRU_SIMPLIFY_FLAG
     &       ,channelcontacty
     &       ,channelcontactz
     &       
     &       ,FluidForceFlag(1:Nr_IBM)
     &       ,ThermalForceFlag(1:Nr_IBM)
     &       ,ContactForceFlag(1:Nr_IBM)
     &       ,ExternalPressureFlag(1:Nr_IBM)
     &       ,Electromechanical(1:Nr_IBM)
     &       ,ImplicitElecMech(1:Nr_IBM)
     &       ,Piezo_Dist_Flag(1:Nr_IBM)
     &       
     &       ,del_S
     &       ,alphaTEMP_IBM
     &       ,betaTEMP_IBM
     &       
     &       ,TSIBMin(1:Nr_IBM)
     &       
     &       ,Del_X
     &       ,Del_Y
     &       ,Del_Z
     &       ,Delta_typeX
     &       ,Delta_typeY
     &       ,Delta_typeZ
     &       
     &       ,AROT_fin  
     &       ,AROT0_fin  
     &       ,maxErrorstruct  
     &       ,maxErrorfluid
     &       ,Solverabstol,Solverreltol,SolverOutiter
     &       ,SolverIniter,SolverPre, SolverILU
     &       ,frequency_fin(1:Nr_IBM,1:3)
     &       
     &       ,maxiter,maxnstruct,maxnfluid
     &       
     &       ,Target_num(1:Nr_IBM)
     &       ,implicitflag(1:Nr_IBM)
     &       
     &       ,k_flag
     &       
     &       ,NBucketnum
     &       ,IBucketnum(1:3)
     &       
     &       ,NBucket(1:Nr_IBM,1:NMaxnumbuck)
     &       ,PBucket(1:Nr_IBM,1:NMaxnumbuck)
     &       
     &       ,Bucketntinterval
     &       ,timebucket
     &       ,BucketXstart(1:3)
     &       ,BucketXEnd(1:3)
     &       
     &       ,p_inflation
     &       ,BucketdX(1:3)
     &       
     &       ,Pnormal_Ext(1:Nr_IBM)
     &       ,n_TDpt_Pnormal_Ext(1:Nr_IBM)
     &       ,flag_TDpt_Pnormal_Ext(1:Nr_IBM)
     &       ,data_TDpt_Pnormal_Ext(1:Nr_IBM,1:1000,1:2)
     &       
     &       ,piezo_beta(1:Nr_IBM)
     &       ,piezo_gamma(1:Nr_IBM)
     &       ,piezo_damp(1:Nr_IBM)
     &       ,piezo_coef(1:Nr_IBM)
     &       ,piezo_alpha11(1:Nr_IBM)
     &       ,piezo_alpha22(1:Nr_IBM)
     &       ,piezo_alpha12(1:Nr_IBM)
     &       ,e4coef(1:Nr_IBM)

      READ (ifuRstrtDistIBMin)
     &       nibmptrecord(1:Nr_IBM)
     &      ,ibmptrecord(1:Nr_IBM,1:100,1:2)

      READ(ifuRstrtDistIBMin)
     &        Genalpha_timesolver(1:Nr_IBM)
     &       ,Genalpha_niter(1:Nr_IBM)
     &       ,Piezo_Formulation_Flag(1:Nr_IBM)

       READ(ifuRstrtDistIBMin)    
     &        contactplanesFlag(1:Nr_IBM)
     &       ,Ncontactplanes(1:Nr_IBM)
     &       ,contactplanesPara(1:Nr_IBM,1:NPlanesIBM,1:6)

      READ (ifuRstrtDistIBMin)
     &        Nr_IBM_f_fbc   
     &       ,Nr_IBM_fb_fbc 
     &       ,Iglbloc_fbc(1:Nr_IBM)
     &       ,Ilocglb_fbc(1:Nr_IBM)

      if(nr_IBM_fib>=1) 
     &           call read_restart_DISRIBM_fib(ifuRstrtDistIBMin)
      if(nr_IBM_fsh>=1) 
     &           call read_restart_DISRIBM_fsh(ifuRstrtDistIBMin)
      if(nr_IBM_esh>=1) 
     &           call read_restart_DISRIBM_esh(ifuRstrtDistIBMin)
      if(nr_IBM_fbc>=1) 
     &           call read_restart_DISRIBM_fbc(ifuRstrtDistIBMin)

         READ (ifuRstrtDistIBMin) nFixcontact,nFixcontactL

         allocate(xFixcontact(nFixcontact,3)
     &           ,EleFixcontact(nFixcontactL,3))

         IF(nFixcontact>0)
     &        READ (ifuRstrtDistIBMin)  
     &        xFixcontact(1:nFixcontact ,1:3)

         IF(nFixcontactL>0)
     &        READ (ifuRstrtDistIBMin)  
     &        ,EleFixcontact(1:nFixcontactL,1:3)

        temp_ibm2=time2


        CLOSE(ifuRstrtDistIBMin) 
       endif 

       end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine write_restart_DISRIBM_fib(ifuRstrtDistIBMOut)
      USE HeaderFSI
      integer  ifuRstrtDistIBMOut
      integer i,j,ijtemp
        WRITE(ifuRstrtDistIBMOut)
     &        GX_IBM_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBM_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBM_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBM_MASSIVE_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBM_MASSIVE_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBM_MASSIVE_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBM_MASSIVEo_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBM_MASSIVEo_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBM_MASSIVEo_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBMo1_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBMo1_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBMo1_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBMpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBMpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBMpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBM_MASSIVEpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBM_MASSIVEpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBM_MASSIVEpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,DSF_IBM_fib(1:Nr_IBM_fib
     &                          ,1-Ns_IBMB_fib:Ns_IBM_fib+Ns_IBMB_fib)    
     &       ,DS_IBM_fib(1:Nr_IBM_fib
     &                          ,1-Ns_IBMB_fib:Ns_IBM_fib+Ns_IBMB_fib)   
     &       ,FIBM1_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIBM2_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIBM3_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FS_1_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FS_2_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)  
     &       ,FS_1_IBMo_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FB_1_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FB_2_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FB_3_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIner_1_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIner_2_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIner_3_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_LINK1_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_LINK2_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_LINK3_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_impuls1_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_impuls2_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_impuls3_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,RIBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM1_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM2_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM3_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM1_pre_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM2_pre_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM3_pre_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,MASS_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,Q_MASSo_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,Q_MASS_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,GX_BP_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_BP_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_BP_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_BP_fib0(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_BP_fib0(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_BP_fib0(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,UIBM1_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,UIBM2_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,UIBM3_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS1_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS2_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS3_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS1o_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS2o_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS3o_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,ffluidsum1_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,ffluidsum2_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,ffluidsum3_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,npos_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,QIBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,TFIBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,TSIBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,Target_K_link_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,Target_t_link_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,phi_fin_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,a_fin_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,Target_points_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,force_points_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,Target_points_v_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,FlagFixed_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,Target_point_num_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,IBucket_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,SBucket_fib(1:Nr_IBM_fib,1:NMaxnumbuck,1:Ns_IBM_fib)
     &       ,Nq_IBM_r_fib(1:Nr_IBM_fib)
     &       ,Ns_IBM_r_fib(1:Nr_IBM_fib)  
     &       ,Ns_IBM_i_fib(1:Nr_IBM_fib)  
     &       ,Ns_IBM_rall_fib(1:Nr_IBM_fib)
       end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine write_restart_DISRIBM_fsh(ifuRstrtDistIBMOut)
      USE HeaderFSI
      integer  ifuRstrtDistIBMOut 
      integer i,j,ijtemp
          WRITE(ifuRstrtDistIBMOut)
     &        GX_IBM_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_IBM_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_IBM_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GX_IBM_MASSIVE_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GY_IBM_MASSIVE_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GZ_IBM_MASSIVE_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GX_IBM_MASSIVEo_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GY_IBM_MASSIVEo_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GZ_IBM_MASSIVEo_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GX_IBMo1_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_IBMo1_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_IBMo1_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GX_IBMpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_IBMpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_IBMpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GX_IBM_MASSIVEpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GY_IBM_MASSIVEpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GZ_IBM_MASSIVEpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,DSF_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)    
     &       ,DS_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)     
     &       ,DSF2_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)    
     &       ,DS2_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)    
     &       ,DKF_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)   
     &       ,DK_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)   
     &       ,DKF2_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)   
     &       ,DK2_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)    
     &       ,FIBM1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIBM2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIBM3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FS_1_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FS_2_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)  
     &       ,FS_1_IBMo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FB_1_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FB_2_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FB_3_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIner_1_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIner_2_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIner_3_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_LINK1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_LINK2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_LINK3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impuls1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impuls2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impuls3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impulsHis1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impulsHis2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impulsHis3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,RIBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,AIBM1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,AIBM2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,AIBM3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM1_pre_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM2_pre_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM3_pre_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,MASS_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Q_MASSo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Q_MASS_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,GX_BP_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_BP_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_BP_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GX_BP_fsh0(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_BP_fsh0(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_BP_fsh0(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       
     &       ,UIBM1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,UIBM2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,UIBM3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS1o_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS2o_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS3o_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ffluidsum1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ffluidsum2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ffluidsum3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       
     &       ,Tzero_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh,1:2,1:2)
     &       ,Bzero_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh,1:2,1:2)
     &       ,npos_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)      
     &       ,BoundaryFibShell_fsh(1:Nr_IBM_fsh,1:2,1:2)
     &       ,QIBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,TFIBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,TSIBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Target_K_link_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &       ,Target_t_link_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &       ,phi_fin_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,a_fin_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,Target_points_fsh
     &                 (1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,force_points_fsh
     &                 (1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,Target_points_v_fsh
     &                 (1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,contact_time_fsh(1:Nr_IBM_fsh,1:2)
     &       ,FlagFixed_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &       ,Target_point_num_fsh(1:Nr_IBM_fsh
     &                          ,1:2
     &                          ,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &       ,IBucket_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,SBucket_fsh(1:Nr_IBM_fsh,1:NMaxnumbuck
     &                          ,1:Nq_IBM_fsh*Ns_IBM_fsh,1:2)
     &       ,Nq_IBM_r_fsh(1:Nr_IBM_fsh)
     &       ,Ns_IBM_r_fsh(1:Nr_IBM_fsh)  
     &       ,Ns_IBM_i_fsh(1:Nr_IBM_fsh)  
     &       ,Ns_IBM_rall_fsh(1:Nr_IBM_fsh) 
     &       ,nMem_Coef_fsh(1:Nr_IBM_fsh)


         WRITE (ifuRstrtDistIBMout)
     &      piezo_Coef_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh,-1:5)
     &       ,Fpiezo1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Fpiezo2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Fpiezo3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ndot_ibm_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ndot_ibm_fsh0(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ndotold_ibm_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ndotpre_ibm_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,vPiezo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,dvPiezo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ddvPiezo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,vPiezoold_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,vPiezopre_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh) 
     &       ,TzeroB_fsh(1:Nr_IBM_fsh,1:4,1:Nq_IBM_fsh+Ns_IBM_fsh,1:3)
         
          do i=1,Nr_IBM_fsh  
         WRITE(ifuRstrtDistIBMOut)     
     &        Mem_Coef_fsh(i,1:Nq_IBM_fsh,1:Ns_IBM_fsh,1:Nmat_IBM_fsh)
     &       ,Ben_Coef_fsh(i,1:Nq_IBM_fsh,1:Ns_IBM_fsh,1:3)
     &       ,contact_coef_fsh(i,1:NPlanesIBM,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
          enddo
          do i=1,Nr_IBM_fsh
          do j=1,Ns_IBM_fsh*Nq_IBM_fsh
         WRITE(ifuRstrtDistIBMOut)    
     &        ContactShellFlag_fsh(i
     &                          ,j
     &                          ,1:Ns_IBM_fsh*Nq_IBM_fsh)
          enddo
          enddo
        
       end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine write_restart_DISRIBM_esh(ifuRstrtDistIBMOut)
      USE HeaderFSI
      integer  ifuRstrtDistIBMOut
      integer i,j,ijtemp
        WRITE(ifuRstrtDistIBMOut)
     &        GX_IBM_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBM_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBM_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBM_MASSIVE_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBM_MASSIVE_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBM_MASSIVE_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBM_MASSIVEo_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBM_MASSIVEo_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBM_MASSIVEo_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBMo1_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBMo1_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBMo1_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBMpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBMpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBMpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBM_MASSIVEpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBM_MASSIVEpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBM_MASSIVEpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,DS_IBM_esh(1:Nr_IBM_esh
     &                          ,1-Ns_IBMB_esh:Ns_IBM_esh+Ns_IBMB_esh)     
     &       ,FIBM1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIBM2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIBM3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FS_1_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FS_2_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)  
     &       ,FS_1_IBMo_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FB_1_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FB_2_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FB_3_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIner_1_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIner_2_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIner_3_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_LINK1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_LINK2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_LINK3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_impuls1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_impuls2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_impuls3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,RIBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,AIBM1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,AIBM2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,AIBM3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM1_pre_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM2_pre_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM3_pre_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,MASS_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,Q_MASSo_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,Q_MASS_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,GX_BP_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_BP_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_BP_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_BP_esh0(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_BP_esh0(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_BP_esh0(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,target_kvalue_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,contact_time_esh(1:Nr_IBM_esh,1:2)
     &       ,contact_coef_esh(1:Nr_IBM_esh,1:NPlanesIBM,0:Ns_IBM_esh+1)
     &       ,UIBM1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,UIBM2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,UIBM3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS1o_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS2o_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS3o_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,ffluidsum1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,ffluidsum2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,ffluidsum3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       
     &       ,coorCurv_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:maxordern,1:2)
     &       ,crFibCrv_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:6)
     &       ,xpcenter_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:3)
     &       ,Jreal_2_imag_esh(1:Nr_IBM_esh,1:Ns_IBM_esh,1:Ns_IBM_esh)  
     &       
     &       ,Dmat0_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:ngaumax,1:2,1:2)   
     &       ,Kmat0_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:ngaumax,1:2,1:2) 
     &       ,BndyCnd_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,neaType_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh)
     &       ,npos_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,target_ktype_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,BoundryList_esh(1:Nr_IBM_esh,1:Ns_IBM_img_esh,1:3)  
     &       ,nea_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:maxordern)
     &       ,nposele_esh(1:Nr_IBM_esh,1:Ns_IBM_esh,1:maxordern)
     &       ,QIBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,TFIBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,TSIBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,Target_K_link_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,Target_t_link_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,phi_fin_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,a_fin_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,Target_points_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)

     &       ,force_points_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,Target_points_v_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,Materialpara_esh(1:Nr_IBM_esh,
     &                         0:Ns_IBM_esh+1,1:Ns_IBM_Fbr_esh)
     &       ,FlagFixed_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,Target_point_num_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,IBucket_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,SBucket_esh(1:Nr_IBM_esh,1:NMaxnumbuck,1:Ns_IBM_esh)
     &       ,Area(1:Nr_IBM_esh,1:Ns_IBMEle_esh)   
     &       ,Nq_IBM_r_esh(1:Nr_IBM_esh)
     &       ,Ns_IBM_r_esh(1:Nr_IBM_esh)  
     &       ,Ns_IBM_i_esh(1:Nr_IBM_esh)  
     &       ,Ns_IBM_rall_esh(1:Nr_IBM_esh)  
     &       ,n_matpara_esh(1:Nr_IBM_esh)
     &       ,Fibrous_esh(1:Nr_IBM_esh)
     &       ,comprerssibleflag_esh(1:Nr_IBM_esh)
     &       ,ShellModelType_esh(1:Nr_IBM_esh)
     &       ,flag_fractal_esh(1:Nr_IBM_esh)

       ijtemp=0
       do i=1,Nr_IBM_esh
        do j=1,Ns_IBMmax_esh
         WRITE(ifuRstrtDistIBMOut)    
     &        ContactShellFlag_esh(i
     &                          ,j
     &                          ,1:Ns_IBMmax_esh)
        enddo
        ijtemp=max(ijtemp,abs(ShellModelType_esh(i)))
       enddo

       if(ijtemp .gt. 1) then
       do i=1,Nr_IBM_esh
        do j=1,Ns_IBMEle_esh
         WRITE(ifuRstrtDistIBMOut)    
     &        gmetric_con0SAVE_esh(i
     &                          ,j
     &                          ,1:ngaumax,1:ngaumaxh,1:2,1:2)
        enddo
       enddo
       do i=1,Nr_IBM_esh
        do j=1,Ns_IBMEle_esh
         WRITE(ifuRstrtDistIBMOut)    
     &        gbase_con0SAVE_esh(i
     &                          ,j
     &                          ,1:ngaumax,1:ngaumaxh,1:2,1:3)
        enddo
       enddo
       do i=1,Nr_IBM_esh
        do j=1,Ns_IBMEle_esh
         WRITE(ifuRstrtDistIBMOut)    
     &        detC_inplane0SAVE_esh(i
     &                          ,j
     &                          ,1:ngaumax,1:ngaumaxh)
        enddo
       enddo

       endif

       if(ANY(flag_fractal_esh(1:Nr_IBM_esh)) ) then
         WRITE(ifuRstrtDistIBMOut)  
     &     nifractal_esh(1:Nr_IBM_esh)
     &    ,fractalinitFlag_esh(1:Nr_IBM_esh)
     &    ,fractalpara_esh(1:Nr_IBM_esh,1:9)
     &    ,Dfractal_esh(1:3,1:Nr_IBM_esh)
     &    ,Jfractal2_esh(1:3,1:Nr_IBM_esh)
     &    ,Dfractal2_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,Hfractal2_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,Ifractal2_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,Bfractal_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,Cfractal_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,fibmsig_esh(1:Nr_IBM_esh,1:3,1:Ns_IBM_esh)
     &    ,Afractal_esh(1:12,1:12,1:3,1:Nr_IBM_esh)
     &    ,Afractal2_esh(1:12,1:12,1:3,1:Nr_IBM_esh)
     &    ,Bfractal2_esh(1:12,1:12,1:3,1:Nr_IBM_esh)
     &    ,Cfractal2_esh(1:12,1:12,1:3,1:Nr_IBM_esh)
     &    ,Xfractal_esh(1:72,1:ngaumax,1:Ns_IBMEle_esh,1:Nr_IBM_esh)
       endif

       end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine write_restart_DISRIBM_fbc(ifuRstrtDistIBMOut)
      USE HeaderFSI
      integer  ifuRstrtDistIBMOut
      integer FabricMeshGeneration,FabricLineGeneration
      integer i,j,ijtemp
      FabricMeshGeneration=0
      FabricLineGeneration=0



        WRITE(ifuRstrtDistIBMOut)
     &        GX_IBM_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBM_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBM_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBM_MASSIVE_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBM_MASSIVE_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBM_MASSIVE_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBM_MASSIVEo_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBM_MASSIVEo_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBM_MASSIVEo_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBMo1_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBMo1_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBMo1_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBMpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBMpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBMpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBM_MASSIVEpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBM_MASSIVEpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBM_MASSIVEpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,DS_IBM_fbc(1:Nr_IBM_fbc
     &                          ,1-Ns_IBMB_fbc:Ns_IBM_fbc+Ns_IBMB_fbc)     
     &       ,FIBM1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIBM2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIBM3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FS_1_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FS_2_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)  
     &       ,FS_1_IBMo_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FB_1_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FB_2_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FB_3_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIner_1_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIner_2_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIner_3_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_LINK1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_LINK2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_LINK3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_impuls1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_impuls2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_impuls3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,RIBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,AIBM1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,AIBM2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,AIBM3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM1_pre_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM2_pre_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM3_pre_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,MASS_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,Q_MASSo_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,Q_MASS_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,GX_BP_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_BP_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_BP_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_BP_fbc0(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_BP_fbc0(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_BP_fbc0(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,target_kvalue_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,contact_time_fbc(1:Nr_IBM_fbc,1:2)
     &       ,contact_coef_fbc(1:Nr_IBM_fbc,1:NPlanesIBM,0:Ns_IBM_fbc+1)
     &       ,UIBM1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,UIBM2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,UIBM3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS1o_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS2o_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS3o_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,ffluidsum1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,ffluidsum2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,ffluidsum3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       
     &       ,Jreal_2_imag_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc,1:Ns_IBM_fbc)  
     &       ,xpcenter_fbc(1:Nr_IBM_fbc,1:Ns_IBMEle_fbc,1:3)
     &       
     &       ,BndyCnd_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,target_ktype_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,BoundryList_fbc(1:Nr_IBM_fbc,1:Ns_IBM_img_fbc,1:3)  
     &       ,nea_fbc(1:Nr_IBM_fbc,1:Ns_IBMEle_fbc,1:3)
     &       ,QIBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,TFIBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,TSIBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,Target_K_link_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,Target_t_link_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,phi_fin_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,a_fin_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,Target_points_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,force_points_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,Target_points_v_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,FlagFixed_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,Target_point_num_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,IBucket_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,SBucket_fbc(1:Nr_IBM_fbc,1:NMaxnumbuck,1:Ns_IBM_fbc)
     &       ,Nq_IBM_r_fbc(1:Nr_IBM_fbc)
     &       ,Ns_IBM_r_fbc(1:Nr_IBM_fbc)  
     &       ,Ns_IBM_i_fbc(1:Nr_IBM_fbc)  
     &       ,Ns_IBM_rall_fbc(1:Nr_IBM_fbc)  
     &       ,n_matpara_fbc(1:Nr_IBM_fbc)

        WRITE(ifuRstrtDistIBMOut)
     &       ,numedge(1:nr_ibm_fbc)
     &       ,numeFacePair(1:nr_ibm_fbc)
     &       ,inplane_type_fbc(1:nr_ibm_fbc)
     &       ,InplaneLinemode_fbc(1:nr_ibm_fbc,1:3)
     &       ,EdgeList_fbc(1:nr_ibm_fbc,1:numedgeMAX_fbc,1:2)
     &       ,facepairEle_fbc(1:nr_ibm_fbc,1:numeFacePairMAX_fbc,1:2)
     &       ,facepairNode_fbc(1:nr_ibm_fbc,1:numeFacePairMAX_fbc,1:4)
     &       ,kappa_volume_fbc(1:nr_ibm_fbc)
     &       ,KBT_fbc(1:nr_ibm_fbc)
     &       ,volume_eq_fbc(1:nr_ibm_fbc)
     &       ,kappa_Area_fbc(1:nr_ibm_fbc)
     &       ,Area_eq_fbc(1:nr_ibm_fbc)
     &       ,kappa_inplane_fbc(1:nr_ibm_fbc,1:3)
     &       ,edge_length_p_fbc(1:nr_ibm_fbc)
     &       ,edge_length_eq_fbc(1:nr_ibm_fbc)
     &       ,edge_length_max_fbc(1:nr_ibm_fbc)
     &       ,InplaneLinem_fbc(1:nr_ibm_fbc) 
     &       ,inplane_Cq_fbc(1:nr_ibm_fbc,1:4)
     &       ,kappa_Bend_fbc(1:nr_ibm_fbc)
     &       ,inplane_E_fbc(1:nr_ibm_fbc)
     &       ,inplane_nu_fbc(1:nr_ibm_fbc)
     &       ,vertex_mass_fbc(1:nr_ibm_fbc)
     &       ,EdgeLength0_fbc(1:nr_ibm_fbc,1:numedgeMAX_fbc)
     &       ,EleArea0_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc)
     &       ,areapnt_fbc(1:nr_ibm_fbc,1:ns_ibm_fbc)        
     &       ,facepairtet0_fbc(1:nr_ibm_fbc,1:numeFacePairMAX_fbc)
     &       ,EleCenter0_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleNormal0_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleEdge0_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleInplaneParaK_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleInplaneParaC_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleInplane_Cq_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,volume_flag_fbc(1:nr_ibm_fbc)
     &       ,Aera_flag_fbc(1:nr_ibm_fbc)
     &       ,inplane_flag_fbc(1:nr_ibm_fbc)
     &       ,Bending_flag_fbc(1:nr_ibm_fbc)
     &       ,tet0_fbc(1:nr_ibm_fbc)

 
       do i=1,Nr_IBM_fbc
        do j=1,Ns_IBMmax_fbc
         WRITE(ifuRstrtDistIBMOut)    
     &        ContactShellFlag_fbc(i
     &                          ,j
     &                          ,1:Ns_IBMmax_fbc)
        enddo
       enddo
 
         WRITE(ifuRstrtDistIBMOut)    
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
     &       ,FabricMesh_Presence
     &       ,FabricLine_Presence
     &       ,FabricMeshFlag_fbc(1:nr_ibm_fbc)
     &       ,FabricLineFlag_fbc(1:nr_ibm_fbc)

         WRITE(ifuRstrtDistIBMOut)    
     &        iMm0,iMc0,iMconR0,iMcon0,iMs0
     &       ,iMpD0_1,iMpD0_1,iMpD0_3
     &       ,n_iMnpara0,iMnpara0(1:100)
     &       ,iMks0,iMl0,n_iMlpara0,iMlpara0(1:100)

         WRITE(ifuRstrtDistIBMOut)    
     &        iLm0,iLc0,iLconR0,iLcon0,iLs0,iLks0,iLkb0,iLl0
     &       ,iLpD0_1,iLpD0_1,iLpD0_3,n_iLnpara0,iLnpara0(1:100)
       
       FabricMeshGeneration=0
       if (FabricMesh_Presence)  then
            FabricMeshGeneration=1
       endif
       FabricLineGeneration=0
       if (FabricLine_Presence)  then
            FabricLineGeneration=1
       endif
    
       if( FabricMeshGeneration .eq. 1) then  
         WRITE(ifuRstrtDistIBMOut)   
     &        FabricMesh_target_k_link(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax)
     &       ,FabricMesh_target_t_link(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax)

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricMesh_coord(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordpre(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordo1(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordMass(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordMasso(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordbp(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordbp0(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_force(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_v(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_a(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_SurfElAttr(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:4)
     &       ,FabricMesh_attr(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax
     &                                 ,1:FabricMesh_numattMax)
     &       ,FabricMeshL_attr(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax
     &                                 ,1:FabricMeshL_numattMax)
     &       ,FabricMesh_a_fin(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax,1:3) 
     &       ,FabricMesh_phi_fin(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax,1:3) 
     &       ,FabricMesh_target_points(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax,1:3) 
     &       ,FabricMesh_target_points_v(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax,1:3) 

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricMesh_numrecord(1:Nr_IBM_fbc)
     &       ,FabricMeshL_numrecord(1:Nr_IBM_fbc)
     &       ,FabricMesh_nMesh(1:Nr_IBM_fbc)
     &       ,FabricMesh_numatt(1:Nr_IBM_fbc)
     &       ,FabricMeshL_numatt(1:Nr_IBM_fbc)
     &       ,FabricMesh_target_num(1:Nr_IBM_fbc)

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricMesh_npoint(1:Nr_IBM_fbc,1:numMeshMax)
     &       ,FabricMesh_nLine(1:Nr_IBM_fbc,1:numMeshMax)
     &       ,FabricMesh_SurfElAddress(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax)
     &       ,FabricMesh_target_point_num(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax)

     &       ,FabricMesh_flagfixed(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax)

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricMesh_address(1:Nr_IBM_fbc
     &                                 ,1:numMeshMax
     &                                 ,1:numMeshpointMax) 
     &       ,FabricMeshL_con(1:Nr_IBM_fbc
     &                                 ,1:numMeshLinkRecordMax,1:2)
     &       ,FabricMeshL_address(1:Nr_IBM_fbc
     &                                 ,1:numMeshMax
     &                                 ,1:numMeshLinkMax) 
     &       ,FabricMeshL_addressRev(1:Nr_IBM_fbc
     &                                 ,1:numMeshLinkRecordMax,1:2)
       endif 
      


       if( FabricLineGeneration .eq. 1) then  

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricLine_target_k_link(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax)
     &       ,FabricLine_target_t_link(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax)

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricLine_coord(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordpre(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordo1(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordMass(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordMasso(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordbp(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordbp0(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_force(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_v(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_a(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_SurfElAttr(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:4)
     &       ,FabricLine_attr(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax
     &                                 ,1:FabricLine_numattMax)  
     &       ,FabricLine_a_fin(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax,1:3)
     &       ,FabricLine_phi_fin(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax,1:3) 
     &       ,FabricLine_target_points(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax,1:3) 
     &       ,FabricLine_target_points_v(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax,1:3)
     &       ,FabricLine_curv0(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricLine_numrecord(1:Nr_IBM_fbc)
     &       ,FabricLine_nLine(1:Nr_IBM_fbc)
     &       ,FabricLine_numatt(1:Nr_IBM_fbc)
     &       ,FabricLine_target_num(1:Nr_IBM_fbc)

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricLine_npoint(1:Nr_IBM_fbc,1:numLineMax)
     &       ,FabricLine_SurfElAddress(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax)
     &       ,FabricLine_target_point_num(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax)
     &       ,FabricLine_flagfixed(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax)

         WRITE(ifuRstrtDistIBMOut) 
     &        FabricLine_address(1:Nr_IBM_fbc
     &                                 ,1:numLineMax
     &                                 ,1:numLinepointMax)  
     &       ,FabricLine_addressRev(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax ,1:2)  
        endif 
 

       end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine read_restart_DISRIBM_fib(ifuRstrtDistIBMin)
      USE HeaderFSI
      integer  ifuRstrtDistIBMin
      integer i,j,ijtemp
       READ (ifuRstrtDistIBMin)
     &        GX_IBM_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBM_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBM_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBM_MASSIVE_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBM_MASSIVE_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBM_MASSIVE_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBM_MASSIVEo_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBM_MASSIVEo_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBM_MASSIVEo_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBMo1_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBMo1_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBMo1_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBMpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBMpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBMpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_IBM_MASSIVEpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_IBM_MASSIVEpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_IBM_MASSIVEpre_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,DSF_IBM_fib(1:Nr_IBM_fib
     &                          ,1-Ns_IBMB_fib:Ns_IBM_fib+Ns_IBMB_fib)    
     &       ,DS_IBM_fib(1:Nr_IBM_fib
     &                          ,1-Ns_IBMB_fib:Ns_IBM_fib+Ns_IBMB_fib)   
     &       ,FIBM1_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIBM2_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIBM3_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FS_1_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FS_2_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)  
     &       ,FS_1_IBMo_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FB_1_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FB_2_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FB_3_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIner_1_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIner_2_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,FIner_3_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_LINK1_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_LINK2_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_LINK3_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_impuls1_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_impuls2_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,F_impuls3_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,RIBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM1_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM2_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM3_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM1_pre_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM2_pre_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,VIBM3_pre_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,MASS_IBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,Q_MASSo_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,Q_MASS_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,GX_BP_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_BP_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_BP_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GX_BP_fib0(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GY_BP_fib0(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,GZ_BP_fib0(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,UIBM1_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,UIBM2_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,UIBM3_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS1_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS2_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS3_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS1o_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS2o_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,FK_MASS3o_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,ffluidsum1_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,ffluidsum2_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,ffluidsum3_fib(1:Nr_IBM_fib,1:Nsec_IBMmax,1:Ns_IBM_fib)
     &       ,npos_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,QIBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,TFIBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,TSIBM_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,Target_K_link_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,Target_t_link_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,phi_fin_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,a_fin_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,Target_points_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,force_points_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,Target_points_v_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1,1:3)
     &       ,FlagFixed_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,Target_point_num_fib(1:Nr_IBM_fib,0:Ns_IBM_fib+1)
     &       ,IBucket_fib(1:Nr_IBM_fib,1:Ns_IBM_fib)
     &       ,SBucket_fib(1:Nr_IBM_fib,1:NMaxnumbuck,1:Ns_IBM_fib)
     &       ,Nq_IBM_r_fib(1:Nr_IBM_fib)
     &       ,Ns_IBM_r_fib(1:Nr_IBM_fib)  
     &       ,Ns_IBM_i_fib(1:Nr_IBM_fib)  
     &       ,Ns_IBM_rall_fib(1:Nr_IBM_fib) 
       end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine read_restart_DISRIBM_fsh(ifuRstrtDistIBMin)
      USE HeaderFSI
      integer  ifuRstrtDistIBMin
      integer i,j,ijtemp
        READ (ifuRstrtDistIBMin)
     &        GX_IBM_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_IBM_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_IBM_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GX_IBM_MASSIVE_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GY_IBM_MASSIVE_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GZ_IBM_MASSIVE_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GX_IBM_MASSIVEo_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GY_IBM_MASSIVEo_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GZ_IBM_MASSIVEo_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GX_IBMo1_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_IBMo1_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_IBMo1_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GX_IBMpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_IBMpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_IBMpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GX_IBM_MASSIVEpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GY_IBM_MASSIVEpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,GZ_IBM_MASSIVEpre_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1
     &                          ,0:Ns_IBM_fsh+1)
     &       ,DSF_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)    
     &       ,DS_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)     
     &       ,DSF2_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)    
     &       ,DS2_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)    
     &       ,DKF_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)   
     &       ,DK_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)   
     &       ,DKF2_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)   
     &       ,DK2_IBM_fsh(1:Nr_IBM_fsh
     &                          ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &                          ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)    
     &       ,FIBM1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIBM2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIBM3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FS_1_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FS_2_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)  
     &       ,FS_1_IBMo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FB_1_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FB_2_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FB_3_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIner_1_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIner_2_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FIner_3_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_LINK1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_LINK2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_LINK3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impuls1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impuls2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impuls3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impulsHis1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impulsHis2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,F_impulsHis3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,RIBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,AIBM1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,AIBM2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,AIBM3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM1_pre_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM2_pre_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,VIBM3_pre_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,MASS_IBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Q_MASSo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Q_MASS_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,GX_BP_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_BP_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_BP_fsh(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GX_BP_fsh0(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GY_BP_fsh0(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       ,GZ_BP_fsh0(1:Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &       
     &       ,UIBM1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,UIBM2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,UIBM3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS1o_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS2o_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,FK_MASS3o_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ffluidsum1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ffluidsum2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ffluidsum3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       
     &       ,Tzero_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh,1:2,1:2)
     &       ,Bzero_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh,1:2,1:2)
     &       ,npos_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)      
     &       ,BoundaryFibShell_fsh(1:Nr_IBM_fsh,1:2,1:2)
     &       ,QIBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,TFIBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,TSIBM_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Target_K_link_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &       ,Target_t_link_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &       ,phi_fin_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,a_fin_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,Target_points_fsh
     &                 (1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,force_points_fsh
     &                 (1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,Target_points_v_fsh
     &                 (1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,1:3)
     &       ,contact_time_fsh(1:Nr_IBM_fsh,1:2)
     &       ,FlagFixed_fsh(1:Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &       ,Target_point_num_fsh(1:Nr_IBM_fsh
     &                          ,1:2
     &                          ,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &       ,IBucket_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,SBucket_fsh(1:Nr_IBM_fsh,1:NMaxnumbuck
     &                          ,1:Nq_IBM_fsh*Ns_IBM_fsh,1:2)
     &       ,Nq_IBM_r_fsh(1:Nr_IBM_fsh)
     &       ,Ns_IBM_r_fsh(1:Nr_IBM_fsh)  
     &       ,Ns_IBM_i_fsh(1:Nr_IBM_fsh)  
     &       ,Ns_IBM_rall_fsh(1:Nr_IBM_fsh) 
     &       ,nMem_Coef_fsh(1:Nr_IBM_fsh)



        READ (ifuRstrtDistIBMin)
     &      piezo_Coef_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh,-1:5)
     &       ,Fpiezo1_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Fpiezo2_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,Fpiezo3_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ndot_ibm_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ndot_ibm_fsh0(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ndotold_ibm_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ndotpre_ibm_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,vPiezo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,dvPiezo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,ddvPiezo_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,vPiezoold_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,vPiezopre_fsh(1:Nr_IBM_fsh,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
     &       ,TzeroB_fsh(1:Nr_IBM_fsh,1:4,1:Nq_IBM_fsh+Ns_IBM_fsh,1:3)
 
 
          do i=1,Nr_IBM_fsh  
         READ(ifuRstrtDistIBMin)      
     &        Mem_Coef_fsh(i,1:Nq_IBM_fsh,1:Ns_IBM_fsh,1:Nmat_IBM_fsh)
     &       ,Ben_Coef_fsh(i,1:Nq_IBM_fsh,1:Ns_IBM_fsh,1:3)
     &       ,contact_coef_fsh(i,1:NPlanesIBM,1:Nq_IBM_fsh,1:Ns_IBM_fsh)
          enddo

          do i=1,Nr_IBM_fsh
          do j=1,Ns_IBM_fsh*Nq_IBM_fsh
         READ(ifuRstrtDistIBMin)    
     &        ContactShellFlag_fsh(i
     &                          ,j
     &                          ,1:Ns_IBM_fsh*Nq_IBM_fsh)
          enddo
          enddo
       end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine read_restart_DISRIBM_esh(ifuRstrtDistIBMin)
      USE HeaderFSI
      integer  ifuRstrtDistIBMin
      integer i,j,ijtemp
        READ (ifuRstrtDistIBMin)
     &        GX_IBM_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBM_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBM_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBM_MASSIVE_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBM_MASSIVE_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBM_MASSIVE_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBM_MASSIVEo_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBM_MASSIVEo_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBM_MASSIVEo_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBMo1_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBMo1_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBMo1_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBMpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBMpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBMpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_IBM_MASSIVEpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_IBM_MASSIVEpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_IBM_MASSIVEpre_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,DS_IBM_esh(1:Nr_IBM_esh
     &                          ,1-Ns_IBMB_esh:Ns_IBM_esh+Ns_IBMB_esh)     
     &       ,FIBM1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIBM2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIBM3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FS_1_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FS_2_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)  
     &       ,FS_1_IBMo_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FB_1_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FB_2_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FB_3_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIner_1_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIner_2_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FIner_3_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_LINK1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_LINK2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_LINK3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_impuls1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_impuls2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,F_impuls3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,RIBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,AIBM1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,AIBM2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,AIBM3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM1_pre_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM2_pre_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,VIBM3_pre_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,MASS_IBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,Q_MASSo_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,Q_MASS_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,GX_BP_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_BP_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_BP_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GX_BP_esh0(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GY_BP_esh0(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,GZ_BP_esh0(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,target_kvalue_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,contact_time_esh(1:Nr_IBM_esh,1:2)
     &       ,contact_coef_esh(1:Nr_IBM_esh,1:NPlanesIBM,0:Ns_IBM_esh+1)
     &       ,UIBM1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,UIBM2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,UIBM3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS1o_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS2o_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,FK_MASS3o_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,ffluidsum1_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,ffluidsum2_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,ffluidsum3_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       
     &       ,coorCurv_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:maxordern,1:2)
     &       ,crFibCrv_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:6)
     &       ,xpcenter_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:3)
     &       ,Jreal_2_imag_esh(1:Nr_IBM_esh,1:Ns_IBM_esh,1:Ns_IBM_esh)  
     &       
     &       ,Dmat0_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:ngaumax,1:2,1:2)   
     &       ,Kmat0_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:ngaumax,1:2,1:2) 
     &       ,BndyCnd_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,neaType_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh)
     &       ,npos_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,target_ktype_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,BoundryList_esh(1:Nr_IBM_esh,1:Ns_IBM_img_esh,1:3)  
     &       ,nea_esh(1:Nr_IBM_esh,1:Ns_IBMEle_esh,1:maxordern)
     &       ,nposele_esh(1:Nr_IBM_esh,1:Ns_IBM_esh,1:maxordern)
     &       ,QIBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,TFIBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,TSIBM_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,Target_K_link_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,Target_t_link_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,phi_fin_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,a_fin_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,Target_points_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,force_points_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,Target_points_v_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1,1:3)
     &       ,Materialpara_esh(1:Nr_IBM_esh,
     &                         0:Ns_IBM_esh+1,1:Ns_IBM_Fbr_esh)
     &       ,FlagFixed_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,Target_point_num_esh(1:Nr_IBM_esh,0:Ns_IBM_esh+1)
     &       ,IBucket_esh(1:Nr_IBM_esh,1:Ns_IBM_esh)
     &       ,SBucket_esh(1:Nr_IBM_esh,1:NMaxnumbuck,1:Ns_IBM_esh)
     &       ,Area(1:Nr_IBM_esh,1:Ns_IBMEle_esh)   
     &       ,Nq_IBM_r_esh(1:Nr_IBM_esh)
     &       ,Ns_IBM_r_esh(1:Nr_IBM_esh)  
     &       ,Ns_IBM_i_esh(1:Nr_IBM_esh)  
     &       ,Ns_IBM_rall_esh(1:Nr_IBM_esh)  
     &       ,n_matpara_esh(1:Nr_IBM_esh)
     &       ,Fibrous_esh(1:Nr_IBM_esh)
     &       ,comprerssibleflag_esh(1:Nr_IBM_esh)
     &       ,ShellModelType_esh(1:Nr_IBM_esh)
     &       ,flag_fractal_esh(1:Nr_IBM_esh)

       ijtemp=0
       do i=1,Nr_IBM_esh
        do j=1,Ns_IBMmax_esh
         READ(ifuRstrtDistIBMin)    
     &        ContactShellFlag_esh(i
     &                          ,j
     &                          ,1:Ns_IBMmax_esh)
        enddo
        ijtemp=max(ijtemp,abs(ShellModelType_esh(i)))
       enddo

       if(ijtemp .gt. 1) then
       do i=1,Nr_IBM_esh
        do j=1,Ns_IBMEle_esh
         READ(ifuRstrtDistIBMin)  
     &        gmetric_con0SAVE_esh(i
     &                          ,j
     &                          ,1:ngaumax,1:ngaumaxh,1:2,1:2)
        enddo
       enddo
       do i=1,Nr_IBM_esh
        do j=1,Ns_IBMEle_esh
         READ(ifuRstrtDistIBMin)    
     &        gbase_con0SAVE_esh(i
     &                          ,j
     &                          ,1:ngaumax,1:ngaumaxh,1:2,1:3)
        enddo
       enddo
       do i=1,Nr_IBM_esh
        do j=1,Ns_IBMEle_esh
         READ(ifuRstrtDistIBMin)    
     &        detC_inplane0SAVE_esh(i
     &                          ,j
     &                          ,1:ngaumax,1:ngaumaxh)
        enddo
       enddo

       endif


       if(ANY(flag_fractal_esh(1:Nr_IBM_esh)) ) then
         READ(ifuRstrtDistIBMin)  
     &     nifractal_esh(1:Nr_IBM_esh)
     &    ,fractalinitFlag_esh(1:Nr_IBM_esh)
     &    ,fractalpara_esh(1:Nr_IBM_esh,1:9)
     &    ,Dfractal_esh(1:3,1:Nr_IBM_esh)
     &    ,Jfractal2_esh(1:3,1:Nr_IBM_esh)
     &    ,Dfractal2_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,Hfractal2_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,Ifractal2_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,Bfractal_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,Cfractal_esh(1:12,1:3,1:Nr_IBM_esh)
     &    ,fibmsig_esh(1:Nr_IBM_esh,1:3,1:Ns_IBM_esh)
     &    ,Afractal_esh(1:12,1:12,1:3,1:Nr_IBM_esh)
     &    ,Afractal2_esh(1:12,1:12,1:3,1:Nr_IBM_esh)
     &    ,Bfractal2_esh(1:12,1:12,1:3,1:Nr_IBM_esh)
     &    ,Cfractal2_esh(1:12,1:12,1:3,1:Nr_IBM_esh)
     &    ,Xfractal_esh(1:72,1:ngaumax,1:Ns_IBMEle_esh,1:Nr_IBM_esh)
       endif



       end subroutine
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine read_restart_DISRIBM_fbc(ifuRstrtDistIBMin)
      USE HeaderFSI
      integer  ifuRstrtDistIBMin
      integer FabricMeshGeneration,FabricLineGeneration
      integer i,j,ijtemp
      FabricMeshGeneration=0
      FabricLineGeneration=0

       READ (ifuRstrtDistIBMin)
     &        GX_IBM_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBM_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBM_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBM_MASSIVE_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBM_MASSIVE_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBM_MASSIVE_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBM_MASSIVEo_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBM_MASSIVEo_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBM_MASSIVEo_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBMo1_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBMo1_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBMo1_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBMpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBMpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBMpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_IBM_MASSIVEpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_IBM_MASSIVEpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_IBM_MASSIVEpre_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,DS_IBM_fbc(1:Nr_IBM_fbc
     &                          ,1-Ns_IBMB_fbc:Ns_IBM_fbc+Ns_IBMB_fbc)     
     &       ,FIBM1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIBM2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIBM3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FS_1_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FS_2_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)  
     &       ,FS_1_IBMo_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FB_1_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FB_2_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FB_3_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIner_1_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIner_2_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FIner_3_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_LINK1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_LINK2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_LINK3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_impuls1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_impuls2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,F_impuls3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,RIBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,AIBM1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,AIBM2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,AIBM3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM1_pre_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM2_pre_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,VIBM3_pre_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,MASS_IBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,Q_MASSo_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,Q_MASS_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,GX_BP_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_BP_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_BP_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GX_BP_fbc0(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GY_BP_fbc0(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,GZ_BP_fbc0(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,target_kvalue_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,contact_time_fbc(1:Nr_IBM_fbc,1:2)
     &       ,contact_coef_fbc(1:Nr_IBM_fbc,1:NPlanesIBM,0:Ns_IBM_fbc+1)
     &       ,UIBM1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,UIBM2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,UIBM3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS1o_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS2o_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,FK_MASS3o_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,ffluidsum1_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,ffluidsum2_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,ffluidsum3_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       
     &       ,Jreal_2_imag_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc,1:Ns_IBM_fbc)  
     &       ,xpcenter_fbc(1:Nr_IBM_fbc,1:Ns_IBMEle_fbc,1:3)
     &       
     &       ,BndyCnd_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,target_ktype_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,BoundryList_fbc(1:Nr_IBM_fbc,1:Ns_IBM_img_fbc,1:3)  
     &       ,nea_fbc(1:Nr_IBM_fbc,1:Ns_IBMEle_fbc,1:3)
     &       ,QIBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,TFIBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,TSIBM_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,Target_K_link_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,Target_t_link_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,phi_fin_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,a_fin_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,Target_points_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,force_points_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,Target_points_v_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1,1:3)
     &       ,FlagFixed_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,Target_point_num_fbc(1:Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &       ,IBucket_fbc(1:Nr_IBM_fbc,1:Ns_IBM_fbc)
     &       ,SBucket_fbc(1:Nr_IBM_fbc,1:NMaxnumbuck,1:Ns_IBM_fbc)
     &       ,Nq_IBM_r_fbc(1:Nr_IBM_fbc)
     &       ,Ns_IBM_r_fbc(1:Nr_IBM_fbc)  
     &       ,Ns_IBM_i_fbc(1:Nr_IBM_fbc)  
     &       ,Ns_IBM_rall_fbc(1:Nr_IBM_fbc)  
     &       ,n_matpara_fbc(1:Nr_IBM_fbc)

       READ (ifuRstrtDistIBMin)
     &       ,numedge(1:nr_ibm_fbc)
     &       ,numeFacePair(1:nr_ibm_fbc)
     &       ,inplane_type_fbc(1:nr_ibm_fbc)
     &       ,InplaneLinemode_fbc(1:nr_ibm_fbc,1:3)
     &       ,EdgeList_fbc(1:nr_ibm_fbc,1:numedgeMAX_fbc,1:2)
     &       ,facepairEle_fbc(1:nr_ibm_fbc,1:numeFacePairMAX_fbc,1:2)
     &       ,facepairNode_fbc(1:nr_ibm_fbc,1:numeFacePairMAX_fbc,1:4)
     &       ,kappa_volume_fbc(1:nr_ibm_fbc)
     &       ,KBT_fbc(1:nr_ibm_fbc)
     &       ,volume_eq_fbc(1:nr_ibm_fbc)
     &       ,kappa_Area_fbc(1:nr_ibm_fbc)
     &       ,Area_eq_fbc(1:nr_ibm_fbc)
     &       ,kappa_inplane_fbc(1:nr_ibm_fbc,1:3)
     &       ,edge_length_p_fbc(1:nr_ibm_fbc)
     &       ,edge_length_eq_fbc(1:nr_ibm_fbc)
     &       ,edge_length_max_fbc(1:nr_ibm_fbc)
     &       ,InplaneLinem_fbc(1:nr_ibm_fbc) 
     &       ,inplane_Cq_fbc(1:nr_ibm_fbc,1:4)
     &       ,kappa_Bend_fbc(1:nr_ibm_fbc)
     &       ,inplane_E_fbc(1:nr_ibm_fbc)
     &       ,inplane_nu_fbc(1:nr_ibm_fbc)
     &       ,vertex_mass_fbc(1:nr_ibm_fbc)
     &       ,EdgeLength0_fbc(1:nr_ibm_fbc,1:numedgeMAX_fbc)
     &       ,EleArea0_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc)
     &       ,areapnt_fbc(1:nr_ibm_fbc,1:ns_ibm_fbc)        
     &       ,facepairtet0_fbc(1:nr_ibm_fbc,1:numeFacePairMAX_fbc)
     &       ,EleCenter0_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleNormal0_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleEdge0_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleInplaneParaK_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleInplaneParaC_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,EleInplane_Cq_fbc(1:nr_ibm_fbc,1:ns_IBMEle_fbc,1:3)
     &       ,volume_flag_fbc(1:nr_ibm_fbc)
     &       ,Aera_flag_fbc(1:nr_ibm_fbc)
     &       ,inplane_flag_fbc(1:nr_ibm_fbc)
     &       ,Bending_flag_fbc(1:nr_ibm_fbc)
     &       ,tet0_fbc(1:nr_ibm_fbc)


       do i=1,Nr_IBM_fbc
        do j=1,Ns_IBMmax_fbc
         READ (ifuRstrtDistIBMin)    
     &        ContactShellFlag_fbc(i
     &                          ,j
     &                          ,1:Ns_IBMmax_fbc)
        enddo
       enddo



         READ (ifuRstrtDistIBMin)    
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
     &       ,FabricMesh_Presence
     &       ,FabricLine_Presence
     &       ,FabricMeshFlag_fbc(1:nr_ibm_fbc)
     &       ,FabricLineFlag_fbc(1:nr_ibm_fbc)

         READ (ifuRstrtDistIBMin)    
     &        iMm0,iMc0,iMconR0,iMcon0,iMs0
     &       ,iMpD0_1,iMpD0_1,iMpD0_3
     &       ,n_iMnpara0,iMnpara0(1:100)
     &       ,iMks0,iMl0,n_iMlpara0,iMlpara0(1:100)

         READ (ifuRstrtDistIBMin)    
     &        iLm0,iLc0,iLconR0,iLcon0,iLs0,iLks0,iLkb0,iLl0
     &       ,iLpD0_1,iLpD0_1,iLpD0_3,n_iLnpara0,iLnpara0(1:100)
       
       FabricMeshGeneration=0
       do i=1,Nr_IBM_fbc 
          if (FabricMeshFlag_fbc(i) )  then
            FabricMeshGeneration=1
            exit
          endif
       enddo
       FabricLineGeneration=0
       do i=1,Nr_IBM_fbc 
          if (FabricLineFlag_fbc(i) )  then
            FabricLineGeneration=1
            exit
          endif
       enddo     
       if( FabricMeshGeneration .eq. 1) then  
                  allocate(
     &                    FabricMesh_target_k_link(Nr_IBM_fbc
     &                                       ,numMeshpointRecordMax)
     &                   ,FabricMesh_target_t_link(Nr_IBM_fbc
     &                                       ,numMeshpointRecordMax)
     &                     )

                  allocate(
     &                    FabricMesh_coord(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_coordpre(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_coordo1(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_coordMass(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_coordMasso(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_coordbp(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_coordbp0(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_force(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_v(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_a(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax
     &                                     ,3)
     &                   ,FabricMesh_t(Nr_IBM_fbc
     &                                     ,numMeshpointRecordMax)
     &                   ,FabricMesh_SurfElAttr(Nr_IBM_fbc
     &                                       ,numMeshpointRecordMax
     &                                       ,4)
     &                   ,FabricMesh_attr(Nr_IBM_fbc
     &                                       ,numMeshpointRecordMax
     &                                       ,FabricMesh_numattMax)
     &                   ,FabricMeshL_attr(Nr_IBM_fbc
     &                                       ,numMeshpointRecordMax
     &                                       ,FabricMeshL_numattMax)
     &                   ,FabricMesh_a_fin(Nr_IBM_fbc
     &                                       ,FabricMesh_numtargetMax
     &                                       ,3) 
     &                   ,FabricMesh_phi_fin(Nr_IBM_fbc
     &                                       ,FabricMesh_numtargetMax
     &                                       ,3) 
     &                   ,FabricMesh_target_points(Nr_IBM_fbc
     &                                       ,FabricMesh_numtargetMax
     &                                       ,3) 
     &                   ,FabricMesh_target_points_v(Nr_IBM_fbc
     &                                       ,FabricMesh_numtargetMax
     &                                       ,3) 
     &                     )


                  allocate(
     &                    FabricMesh_numrecord(Nr_IBM_fbc)
     &                   ,FabricMeshL_numrecord(Nr_IBM_fbc)
     &                   ,FabricMesh_nMesh(Nr_IBM_fbc)
     &                   ,FabricMesh_numatt(Nr_IBM_fbc)
     &                   ,FabricMeshL_numatt(Nr_IBM_fbc)
     &                   ,FabricMesh_target_num(Nr_IBM_fbc)
     &                     )


                  allocate(
     &                    FabricMesh_npoint(Nr_IBM_fbc
     &                                     ,numMeshMax)
     &                   ,FabricMesh_nLine(Nr_IBM_fbc
     &                                     ,numMeshMax)
     &                   ,FabricMesh_SurfElAddress(Nr_IBM_fbc
     &                                       ,numMeshpointRecordMax)
     &                   ,FabricMesh_target_point_num(Nr_IBM_fbc
     &                                       ,FabricMesh_numtargetMax)
     &                   ,FabricMesh_flagfixed(Nr_IBM_fbc
     &                                       ,FabricMesh_numtargetMax)
     &                     )


                  allocate(
     &                    FabricMesh_address(Nr_IBM_fbc
     &                                       ,numMeshMax
     &                                       ,numMeshpointMax) 
     &                   ,FabricMeshL_con(Nr_IBM_fbc
     &                                       ,numMeshLinkRecordMax
     &                                       ,2)
     &                   ,FabricMeshL_address(Nr_IBM_fbc
     &                                       ,numMeshMax
     &                                       ,numMeshLinkMax) 
     &                   ,FabricMeshL_addressRev(Nr_IBM_fbc
     &                                       ,numMeshLinkRecordMax
     &                                       ,2)
     &                     )

         READ (ifuRstrtDistIBMin)   
     &        FabricMesh_target_k_link(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax)
     &       ,FabricMesh_target_t_link(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax)

         READ (ifuRstrtDistIBMin) 
     &        FabricMesh_coord(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordpre(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordo1(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordMass(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordMasso(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordbp(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_coordbp0(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_force(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_v(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_a(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:3)
     &       ,FabricMesh_SurfElAttr(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax,1:4)
     &       ,FabricMesh_attr(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax
     &                                 ,1:FabricMesh_numattMax)
     &       ,FabricMeshL_attr(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax
     &                                 ,1:FabricMeshL_numattMax)
     &       ,FabricMesh_a_fin(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax,1:3) 
     &       ,FabricMesh_phi_fin(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax,1:3) 
     &       ,FabricMesh_target_points(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax,1:3) 
     &       ,FabricMesh_target_points_v(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax,1:3) 

         READ (ifuRstrtDistIBMin) 
     &        FabricMesh_numrecord(1:Nr_IBM_fbc)
     &       ,FabricMeshL_numrecord(1:Nr_IBM_fbc)
     &       ,FabricMesh_nMesh(1:Nr_IBM_fbc)
     &       ,FabricMesh_numatt(1:Nr_IBM_fbc)
     &       ,FabricMeshL_numatt(1:Nr_IBM_fbc)
     &       ,FabricMesh_target_num(1:Nr_IBM_fbc)

         READ (ifuRstrtDistIBMin) 
     &        FabricMesh_npoint(1:Nr_IBM_fbc,1:numMeshMax)
     &       ,FabricMesh_nLine(1:Nr_IBM_fbc,1:numMeshMax)
     &       ,FabricMesh_SurfElAddress(1:Nr_IBM_fbc
     &                                 ,1:numMeshpointRecordMax)
     &       ,FabricMesh_target_point_num(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax)
     &       ,FabricMesh_flagfixed(1:Nr_IBM_fbc
     &                                 ,1:FabricMesh_numtargetMax)

         READ (ifuRstrtDistIBMin) 
     &        FabricMesh_address(1:Nr_IBM_fbc
     &                                 ,1:numMeshMax
     &                                 ,1:numMeshpointMax) 
     &       ,FabricMeshL_con(1:Nr_IBM_fbc
     &                                 ,1:numMeshLinkRecordMax,1:2)
     &       ,FabricMeshL_address(1:Nr_IBM_fbc
     &                                 ,1:numMeshMax
     &                                 ,1:numMeshLinkMax) 
     &       ,FabricMeshL_addressRev(1:Nr_IBM_fbc
     &                                 ,1:numMeshLinkRecordMax,1:2)
       endif 

       if( FabricLineGeneration .eq. 1) then  
                  allocate(
     &                    FabricLine_target_k_link(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax)
     &                   ,FabricLine_target_t_link(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax)
     &                    )

                  allocate(
     &                    FabricLine_coord(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_coordpre(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_coordo1(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_coordMass(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_coordMasso(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_coordbp(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_coordbp0(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_force(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_v(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_a(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                   ,FabricLine_t(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax)
     &                   ,FabricLine_SurfElAttr(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,4)
     &                   ,FabricLine_attr(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,FabricLine_numattMax)  
     &                   ,FabricLine_a_fin(Nr_IBM_fbc
     &                                       ,FabricLine_numtargetMax
     &                                       ,3)
     &                   ,FabricLine_phi_fin(Nr_IBM_fbc
     &                                       ,FabricLine_numtargetMax
     &                                       ,3) 
     &                   ,FabricLine_target_points(Nr_IBM_fbc
     &                                       ,FabricLine_numtargetMax
     &                                       ,3) 
     &                   ,FabricLine_target_points_v(Nr_IBM_fbc
     &                                       ,FabricLine_numtargetMax
     &                                       ,3)
     &                   ,FabricLine_curv0(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax
     &                                       ,3)
     &                    )

                  allocate(
     &                    FabricLine_numrecord(Nr_IBM_fbc)
     &                   ,FabricLine_nLine(Nr_IBM_fbc)
     &                   ,FabricLine_numatt(Nr_IBM_fbc)
     &                   ,FabricLine_target_num(Nr_IBM_fbc)
     &                    )


                  allocate(
     &                    FabricLine_npoint(Nr_IBM_fbc
     &                                       ,numLineMax)
     &                   ,FabricLine_SurfElAddress(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax)
     &                   ,FabricLine_target_point_num(Nr_IBM_fbc
     &                                       ,FabricLine_numtargetMax)
     &                   ,FabricLine_flagfixed(Nr_IBM_fbc
     &                                       ,FabricLine_numtargetMax)
     &                    )

                  allocate(
     &                    FabricLine_address(Nr_IBM_fbc
     &                                       ,numLineMax
     &                                       ,numLinepointMax)  
     &                   ,FabricLine_addressRev(Nr_IBM_fbc
     &                                       ,numLinepointRecordMax 
     &                                       ,2)  
     &                    )

         READ (ifuRstrtDistIBMin) 
     &        FabricLine_target_k_link(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax)
     &       ,FabricLine_target_t_link(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax)

         READ (ifuRstrtDistIBMin) 
     &        FabricLine_coord(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordpre(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordo1(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordMass(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordMasso(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordbp(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_coordbp0(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_force(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_v(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_a(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)
     &       ,FabricLine_SurfElAttr(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:4)
     &       ,FabricLine_attr(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax
     &                                 ,1:FabricLine_numattMax)  
     &       ,FabricLine_a_fin(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax,1:3)
     &       ,FabricLine_phi_fin(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax,1:3) 
     &       ,FabricLine_target_points(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax,1:3) 
     &       ,FabricLine_target_points_v(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax,1:3)
     &       ,FabricLine_curv0(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax,1:3)

         READ (ifuRstrtDistIBMin) 
     &        FabricLine_numrecord(1:Nr_IBM_fbc)
     &       ,FabricLine_nLine(1:Nr_IBM_fbc)
     &       ,FabricLine_numatt(1:Nr_IBM_fbc)
     &       ,FabricLine_target_num(1:Nr_IBM_fbc)

         READ (ifuRstrtDistIBMin) 
     &        FabricLine_npoint(1:Nr_IBM_fbc,1:numLineMax)
     &       ,FabricLine_SurfElAddress(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax)
     &       ,FabricLine_target_point_num(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax)
     &       ,FabricLine_flagfixed(1:Nr_IBM_fbc
     &                                 ,1:FabricLine_numtargetMax)

         READ (ifuRstrtDistIBMin) 
     &        FabricLine_address(1:Nr_IBM_fbc
     &                                 ,1:numLineMax
     &                                 ,1:numLinepointMax)  
     &       ,FabricLine_addressRev(1:Nr_IBM_fbc
     &                                 ,1:numLinepointRecordMax ,1:2)  
        endif 
        end subroutine 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|


