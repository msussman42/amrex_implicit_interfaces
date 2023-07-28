c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!@c
      subroutine DistFSI_allocate_memory(iFlagin)
      USE HeaderFSI
      implicit none
      integer iFlagin
      logical, save :: FirstCall_Frac = .TRUE.
!@T Allocate global variables
!@q
    
      INTEGER, PARAMETER :: ERR_NONE = 0
      INTEGER  iERR

      if(iFlagin .eq. 0) then
! IBM parameters
      ALLOCATE(Density_Coef(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Density_Coef'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Density_coefP(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Density_coefP'
        STOP
      ENDIF ! iErr 

      ALLOCATE(time_release(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for time_release'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Ks_IBM_r(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ks_IBM_r'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Kb_IBM_r(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Kb_IBM_r'
        STOP
      ENDIF ! iErr 

      ALLOCATE(K_MASSIVE_IBM_r(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for K_MASSIVE_IBM_r'
        STOP
      ENDIF ! iErr 

      ALLOCATE(cs_ibm(nr_ibm),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for cs_ibm'
        STOP
      ENDIF ! iErr 

      ALLOCATE(cs_ibm_target(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for cs_ibm_target'
        STOP
      ENDIF ! iErr 
 
                               
      ALLOCATE(Phi_ben_coef(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Phi_ben_coef'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(Phi_nus_coef(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Phi_nus_coef'
        STOP
      ENDIF ! iErr 

      ALLOCATE(wei(ngaumax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for wei'
        STOP
      ENDIF ! iErr 

      ALLOCATE(weiCont(ngaumax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for weiCont'
        STOP
      ENDIF ! iErr 

      ALLOCATE(wei_h(ngaumaxh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for wei_h'
        STOP
      ENDIF ! iErr

      ALLOCATE(gau_h(ngaumaxh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for gau_h'
        STOP
      ENDIF ! iErr
 
      ALLOCATE(wei0(ngaumax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for wei0'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ccont(0:20),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ccont'
        STOP
      ENDIF ! iErr

      ALLOCATE(rot_seq(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for rot_seq'
        STOP
      ENDIF ! iErr    
      ALLOCATE(alpham_Coef(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for alpham_Coef'
        STOP
      ENDIF ! iErr

      ALLOCATE(alphaf_Coef(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for alphaf_Coef'
        STOP
      ENDIF ! iErr

      ALLOCATE(beta_Coef(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for beta_Coef'
        STOP
      ENDIF ! iErr

      ALLOCATE(gamma_Coef(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for gamma_Coef'
        STOP
      ENDIF ! iErr

      ALLOCATE(Gt0(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Gt0'
        STOP
      ENDIF ! iErr    
      ALLOCATE(Tramps(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Tramps'
        STOP
      ENDIF ! iErr

      ALLOCATE(Trampx(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Trampx'
        STOP
      ENDIF ! iErr

      ALLOCATE(ampx(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ampx'
        STOP
      ENDIF ! iErr

      ALLOCATE(freqx(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for freqx'
        STOP
      ENDIF ! iErr

      ALLOCATE(phix(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for phix'
        STOP
      ENDIF ! iErr

      ALLOCATE(Trampt(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Trampt'
        STOP
      ENDIF ! iErr

      ALLOCATE(ampt(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ampt'
        STOP
      ENDIF ! iErr

      ALLOCATE(freqt(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for freqt'
        STOP
      ENDIF ! iErr

      ALLOCATE(phit(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for phit'
        STOP
      ENDIF ! iErr

      ALLOCATE(propDamp(Nr_IBM,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for propDamp'
        STOP
      ENDIF ! iErr

      ALLOCATE(Area(Nr_IBM_esh,Ns_IBMEle_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Area'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Phi_mem_coef(Nr_IBM,25),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Phi_mem_coef'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(gau(2,ngaumax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for gau'
        STOP
      ENDIF ! iErr 

      ALLOCATE(gauCont(2,ngaumax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for gauCont'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(gau0(2,ngaumax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for gau0'
        STOP
      ENDIF ! iErr 


!%===============================================================
      ALLOCATE(GX_IBM_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBM_MASSIVE_fib(Nr_IBM_fib,0:Ns_IBM_fib+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVE_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_MASSIVE_fib(Nr_IBM_fib,0:Ns_IBM_fib+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVE_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_MASSIVE_fib(Nr_IBM_fib,0:Ns_IBM_fib+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVE_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBM_MASSIVEo_fib(Nr_IBM_fib,0:Ns_IBM_fib+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVEo_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_MASSIVEo_fib(Nr_IBM_fib,0:Ns_IBM_fib+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVEo_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_MASSIVEo_fib(Nr_IBM_fib,0:Ns_IBM_fib+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVEo_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBMo1_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBMo1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBMo1_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBMo1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBMo1_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBMo1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBMpre_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBMpre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBMpre_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBMpre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBMpre_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBMpre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GX_IBM_MASSIVEpre_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVEpre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GY_IBM_MASSIVEpre_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVEpre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GZ_IBM_MASSIVEpre_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVEpre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DSF_IBM_fib( 
     &Nr_IBM_fib,1-Ns_IBMB_fib:Ns_IBM_fib+Ns_IBMB_fib 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DSF_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DS_IBM_fib( 
     &Nr_IBM_fib,1-Ns_IBMB_fib:Ns_IBM_fib+Ns_IBMB_fib 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DS_IBM_fib'
        STOP
      ENDIF ! iErr 
 


      ALLOCATE(FIBM1_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM2_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM2_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM3_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM3_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FS_1_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_1_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FS_2_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_2_IBM_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(FS_1_IBMo_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_1_IBMo_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_1_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_1_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_2_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_2_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_3_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_3_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_1_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_1_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_2_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_2_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_3_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_3_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK1_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK2_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK2_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK3_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK3_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls1_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls2_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls2_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls3_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls3_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(RIBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for RIBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM1_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM2_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM2_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM3_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM3_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM1_pre_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM1_pre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM2_pre_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM2_pre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM3_pre_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM3_pre_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(MASS_IBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for MASS_IBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Q_MASSo_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Q_MASSo_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Q_MASS_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Q_MASS_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_BP_fib(Nr_ibm_fib,0:NS_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_BP_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_BP_fib(Nr_ibm_fib,0:NS_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_BP_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_BP_fib(Nr_ibm_fib,0:NS_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_BP_fib'
        STOP
      ENDIF ! iErr 


      ALLOCATE(GX_BP_fib0(Nr_ibm_fib,0:NS_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_BP_fib0'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_BP_fib0(Nr_ibm_fib,0:NS_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_BP_fib0'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_BP_fib0(Nr_ibm_fib,0:NS_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_BP_fib0'
        STOP
      ENDIF ! iErr 

!%===============================================================
      ALLOCATE(GX_IBM_fsh(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_fsh(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_fsh(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GX_IBM_MASSIVE_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVE_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GY_IBM_MASSIVE_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVE_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GZ_IBM_MASSIVE_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVE_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GX_IBM_MASSIVEo_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVEo_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GY_IBM_MASSIVEo_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVEo_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GZ_IBM_MASSIVEo_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVEo_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GX_IBMo1_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBMo1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GY_IBMo1_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBMo1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GZ_IBMo1_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBMo1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GX_IBMpre_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBMpre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GY_IBMpre_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBMpre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &         GZ_IBMpre_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBMpre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GX_IBM_MASSIVEpre_fsh(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVEpre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GY_IBM_MASSIVEpre_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVEpre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GZ_IBM_MASSIVEpre_fsh
     &         (Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVEpre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DSF_IBM_fsh( 
     &Nr_IBM_fsh
     &         ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &         ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DSF_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DS_IBM_fsh( 
     &Nr_IBM_fsh
     &         ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &         ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DS_IBM_fsh'
        STOP
      ENDIF ! iErr 
 



      ALLOCATE(DSF2_IBM_fsh( 
     &Nr_IBM_fsh
     &         ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &         ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DSF2_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DS2_IBM_fsh( 
     &Nr_IBM_fsh
     &         ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &         ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DS2_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DKF_IBM_fsh( 
     &Nr_IBM_fsh
     &         ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &         ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DKF_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DK_IBM_fsh( 
     &Nr_IBM_fsh
     &         ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &         ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DK_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DKF2_IBM_fsh( 
     &Nr_IBM_fsh
     &         ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &         ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DKF2_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DK2_IBM_fsh( 
     &Nr_IBM_fsh
     &         ,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &         ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DK2_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM3_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FS_1_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_1_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FS_2_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_2_IBM_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(FS_1_IBMo_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_1_IBMo_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_1_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_1_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_2_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_2_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_3_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_3_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_1_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_1_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_2_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_2_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_3_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_3_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK3_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impulsHis1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impulsHis1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impulsHis2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impulsHis2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impulsHis3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impulsHis3_fsh'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(F_impuls1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls3_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(RIBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for RIBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM3_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(AIBM1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(AIBM2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(AIBM3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM3_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM1_pre_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM1_pre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM2_pre_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM2_pre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM3_pre_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM3_pre_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(MASS_IBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for MASS_IBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Q_MASSo_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Q_MASSo_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Q_MASS_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Q_MASS_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_BP_fsh(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_BP_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_BP_fsh(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_BP_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_BP_fsh(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_BP_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(GX_BP_fsh0(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_BP_fsh0'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_BP_fsh0(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_BP_fsh0'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_BP_fsh0(Nr_IBM_fsh,0:Nq_IBM_fsh+1,0:Ns_IBM_fsh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_BP_fsh0'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(GX_IBM_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBM_MASSIVE_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVE_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_MASSIVE_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVE_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_MASSIVE_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVE_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBM_MASSIVEo_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVEo_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_MASSIVEo_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVEo_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_MASSIVEo_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVEo_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBMo1_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBMo1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBMo1_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBMo1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBMo1_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBMo1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBMpre_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBMpre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBMpre_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBMpre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBMpre_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBMpre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GX_IBM_MASSIVEpre_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVEpre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GY_IBM_MASSIVEpre_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVEpre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GZ_IBM_MASSIVEpre_esh(Nr_IBM_esh,0:Ns_IBM_esh+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVEpre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DS_IBM_esh( 
     &Nr_IBM_esh,1-Ns_IBMB_esh:Ns_IBM_esh+Ns_IBMB_esh 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DS_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM1_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM2_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM2_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM3_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM3_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FS_1_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_1_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FS_2_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_2_IBM_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(FS_1_IBMo_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_1_IBMo_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_1_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_1_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_2_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_2_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_3_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_3_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_1_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_1_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_2_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_2_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_3_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_3_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK1_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK2_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN

        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK2_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK3_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK3_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls1_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls2_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls2_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls3_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls3_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(RIBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for RIBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM1_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM2_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM2_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM3_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM3_esh'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(AIBM1_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(AIBM2_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM2_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(AIBM3_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM3_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(VIBM1_pre_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM1_pre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM2_pre_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM2_pre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM3_pre_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM3_pre_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(MASS_IBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for MASS_IBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Q_MASSo_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Q_MASSo_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Q_MASS_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Q_MASS_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_BP_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_BP_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_BP_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_BP_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_BP_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_BP_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(GX_BP_esh0(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_BP_esh0'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_BP_esh0(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_BP_esh0'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_BP_esh0(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_BP_esh0'
        STOP
      ENDIF ! iErr

      ALLOCATE(target_kvalue_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for target_kvalue_esh'
        STOP
      ENDIF ! iErr

!!?!!      ALLOCATE(contact_time_esh(Nr_IBM_esh,2),STAT=iErr)
!!?!!
!!?!!      IF ( iErr /= ERR_NONE ) THEN
!!?!!        WRITE(*,*) 
!!?!!     &  'Allocate_memory: Memory Allocation Error for contact_time_esh'
!!?!!        STOP
!!?!!      ENDIF ! iErr 

      ALLOCATE(
     & contact_coef_esh(Nr_IBM_esh,NPlanesIBM,0:Ns_IBM_esh+1),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for contact_coef_esh'
        STOP
      ENDIF ! iErr 
!%===============================================================
 
      ALLOCATE(Nn(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Nn'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(dNndv(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for dNndv'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(dNndw(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for dNndw'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(ddNndvdv(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndvdv'
        STOP
      ENDIF ! iErr

      ALLOCATE(ddNndwdw(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndwdw'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(ddNndvdw(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndvdw'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(NnCont(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for NnCont'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(dNndvCont(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for dNndvCont'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(dNndwCont(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for dNndwCont'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(ddNndvdvCont(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndvdvCont'
        STOP
      ENDIF ! iErr

      ALLOCATE(ddNndwdwCont(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndwdwCont'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(ddNndvdwCont(jTypeMax,ngaumax,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndvdwCont'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(Nn0(jTypeMax,1,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Nn0'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(dNndv0(jTypeMax,1,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for dNndv0'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(dNndw0(jTypeMax,1,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for dNndw0'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(ddNndvdv0(jTypeMax,1,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndvdv0'
        STOP
      ENDIF ! iErr

      ALLOCATE(ddNndwdw0(jTypeMax,1,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndwdw0'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(ddNndvdw0(jTypeMax,1,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddNndvdw0'
        STOP
      ENDIF ! iErr

      ALLOCATE(Jreal_2_imag_esh(Nr_IBM_esh,Ns_IBM_esh,Ns_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Jreal_2_imag_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Tzero_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh,2,2)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Tzero_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Bzero_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh,2,2)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Bzero_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Mem_Coef_fsh
     &           (Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh,Nmat_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Mem_Coef_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Ben_Coef_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh,3)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ben_Coef_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &     contact_coef_fsh(Nr_IBM_fsh,NPlanesIBM,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for contact_coef_fsh'
        STOP
      ENDIF ! iErr 
!==================
      ALLOCATE(UIBM1_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(UIBM2_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM2_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(UIBM3_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM3_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS1_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS1_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS2_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS2_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS3_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS3_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS1o_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS1o_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(FK_MASS2o_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS2o_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(FK_MASS3o_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS3o_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(ffluidsum1_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum1_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(ffluidsum2_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum2_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(ffluidsum3_fib(Nr_IBM_fib,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum3_fib'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(UIBM1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(UIBM2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(UIBM3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM3_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS1_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS3_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS1o_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS1o_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(FK_MASS2o_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS2o_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(FK_MASS3o_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS3o_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(ffluidsum1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum1_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(ffluidsum2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum2_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(ffluidsum3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum3_fsh'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(UIBM1_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(UIBM2_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM2_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(UIBM3_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM3_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS1_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS1_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS2_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS2_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS3_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS3_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS1o_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS1o_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(FK_MASS2o_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS2o_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(FK_MASS3o_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS3o_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(ffluidsum1_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum1_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(ffluidsum2_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum2_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(ffluidsum3_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum3_esh'
        STOP
      ENDIF ! iErr 
!%===============================================================


      ALLOCATE(coorCurv_esh(Nr_IBM_esh,Ns_IBMEle_esh,maxordern,2)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for coorCurv_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(crFibCrv_esh(Nr_IBM_esh,Ns_IBMEle_esh,6)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for crFibCrv_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(xpcenter_esh (Nr_IBM_esh,Ns_IBMEle_esh,3)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for xpcenter_esh'
        STOP
      ENDIF ! iErr 

 
      ALLOCATE(Dmat0_esh(Nr_IBM_esh,Ns_IBMEle_esh,ngaumax,2,2)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Dmat0_esh'
        STOP
      ENDIF ! iErr 
 
   
      ALLOCATE(Kmat0_esh(Nr_IBM_esh,Ns_IBMEle_esh,ngaumax,2,2)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Kmat0_esh'
        STOP
      ENDIF ! iErr 
!%===============================================================

      ALLOCATE(left_bndy(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for left_bndy'
        STOP
      ENDIF ! iErr 

      ALLOCATE(right_bndy(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 

     &  'Allocate_memory: Memory Allocation Error for right_bndy'
        STOP
      ENDIF ! iErr 

      ALLOCATE(bottom_bndy(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 


     &  'Allocate_memory: Memory Allocation Error for bottom_bndy'
        STOP
      ENDIF ! iErr 

      ALLOCATE(top_bndy(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for top_bndy'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(Nq_IBM_r(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Nq_IBM_r'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Ns_IBM_r(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_r'
        STOP
      ENDIF ! iErr

      ALLOCATE(NS_IBM_i(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for NS_IBM_i'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ns_IBM_rall(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_rall'
        STOP
      ENDIF ! iErr

      ALLOCATE(FiberBdry(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FiberBdry'
        STOP
      ENDIF ! iErr 

      ALLOCATE(BodyType(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for BodyType'
        STOP
      ENDIF ! iErr 

      ALLOCATE(numelr(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for numelr'
        STOP
      ENDIF ! iErr

      ALLOCATE(numeli(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for numeli'
        STOP
      ENDIF ! iErr

      ALLOCATE(jTypeAll(jTypeMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for jTypeAll'
        STOP
      ENDIF ! iErr 

      ALLOCATE( MasterBdy(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for  MasterBdy'
        STOP
      ENDIF ! iErr 
 
 
      ALLOCATE(iflaginext(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for iflaginext'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(BendingZero(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for BendingZero'
        STOP
      ENDIF ! iErr

      ALLOCATE(MaterialTypeIBM(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for MaterialTypeIBM'
        STOP
      ENDIF ! iErr 


      ALLOCATE(Nq_IBM_r_fib(Nr_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Nq_IBM_r_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Ns_IBM_r_fib(Nr_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_r_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(NS_IBM_i_fib(Nr_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for NS_IBM_i_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ns_IBM_rall_fib(Nr_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_rall_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(Nq_IBM_r_fsh(Nr_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Nq_IBM_r_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Ns_IBM_r_fsh(Nr_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_r_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(NS_IBM_i_fsh(Nr_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for NS_IBM_i_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ns_IBM_rall_fsh(Nr_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_rall_fsh'
        STOP
      ENDIF ! iErr 


      ALLOCATE(Nq_IBM_r_esh(Nr_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Nq_IBM_r_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Ns_IBM_r_esh(Nr_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_r_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(NS_IBM_i_esh(Nr_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for NS_IBM_i_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ns_IBM_rall_esh(Nr_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_rall_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Iglbloc_fib(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Iglbloc_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ilocglb_fib(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ilocglb_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(Iglbloc_esh(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Iglbloc_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ilocglb_esh(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ilocglb_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Iglbloc_fsh(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Iglbloc_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ilocglb_fsh(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ilocglb_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(nMem_Coef_fsh(Nr_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN

        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for nMem_Coef_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(n_matpara_esh(Nr_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN

        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for n_matpara_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Fibrous_esh(Nr_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Fibrous_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(comprerssibleflag_esh(Nr_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &'Allocate_memory: Memory Allocation Error',
     &' for comprerssibleflag_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(ShellModelType_esh(Nr_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for ShellModelType_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Genalpha_timesolver(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error',
     &  'for Genalpha_timesolver'
        STOP
      ENDIF ! iErr

      ALLOCATE(Genalpha_niter(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error',
     &  'for Genalpha_niter'
        STOP
      ENDIF ! iErr

      ALLOCATE(Piezo_Formulation_Flag(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error',
     &  'for Piezo_Formulation_Flag'
        STOP
      ENDIF ! iErr 



!%===============================================================

      ALLOCATE(BndyCnd_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for BndyCnd_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(neaType_esh(Nr_IBM_esh,Ns_IBMEle_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for neaType_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(BoundryList_esh(Nr_IBM_esh,Ns_IBM_img_esh,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for BoundryList_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(nea_esh(Nr_IBM_esh,Ns_IBMEle_esh,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for nea_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(npos_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for npos_esh'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(npos_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for npos_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(target_ktype_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for target_ktype_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(npos_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for npos_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(nposele_esh(Nr_IBM_esh,Ns_IBM_esh,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for nposele_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &     ContactShellFlag_fsh(
     &          Nr_IBM_fsh
     &          ,Ns_IBM_fsh*Nq_IBM_fsh
     &          ,Ns_IBM_fsh*Nq_IBM_fsh)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'ContactShellFlag_fsh'
        STOP
      ENDIF ! iErr 
 
      Ns_IBMmax_esh=max(Ns_IBMEle_esh,Ns_IBM_esh)
      ALLOCATE(
     &          ContactShellFlag_esh
     &          (Nr_IBM_esh,Ns_IBMmax_esh,Ns_IBMmax_esh)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'ContactShellFlag_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(BoundaryFibShell_fsh(Nr_IBM_fsh,2,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error ', 
     &          'for BoundaryFibShell_fsh'
        STOP
      ENDIF ! iErr 

!%===============================================================


      ALLOCATE(FluidForceFlag(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FluidForceFlag'
        STOP
      ENDIF ! iErr 

      ALLOCATE(ThermalForceFlag(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ThermalForceFlag'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(ContactForceFlag(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ContactForceFlag'
        STOP
      ENDIF ! iErr 

      ALLOCATE(ExternalPressureFlag(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &'Allocate_memory: Memory Allocation Error  ExternalPressureFlag'
        STOP
      ENDIF ! iErr

      ALLOCATE(Electromechanical(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Electromechanical'
        STOP
      ENDIF ! iErr

      ALLOCATE(ImplicitElecMech(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ImplicitElecMech'
        STOP
      ENDIF ! iErr


      ALLOCATE(Piezo_Dist_Flag(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &'Allocate_memory: Memory Allocation Error for Piezo_Dist_Flag'
        STOP
      ENDIF ! iErr 


C IBM parameters
      ALLOCATE(TSIBMin(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TSIBMin'
        STOP
      ENDIF ! iErr 


      ALLOCATE(QIBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for QIBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(TFIBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TFIBM_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(TSIBM_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TSIBM_fib'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(QIBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for QIBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(TFIBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TFIBM_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(TSIBM_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TSIBM_fsh'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(QIBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for QIBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(TFIBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TFIBM_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(TSIBM_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TSIBM_esh'
        STOP
      ENDIF ! iErr 
!%===============================================================


C Target Points
      ALLOCATE(frequency_fin(Nr_IBM,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for frequency_fin'
        STOP
      ENDIF ! iErr 


      ALLOCATE(Target_K_link_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_K_link_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Target_t_link_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_t_link_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(phi_fin_fib(Nr_IBM_fib,0:Ns_IBM_fib+1,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for phi_fin_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(a_fin_fib(Nr_IBM_fib,0:Ns_IBM_fib+1,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for a_fin_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(Target_points_fib(Nr_IBM_fib,0:Ns_IBM_fib+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_points_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(force_points_fib(Nr_IBM_fib,0:Ns_IBM_fib+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for force_points_fib'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Target_points_v_fib(Nr_IBM_fib,0:Ns_IBM_fib+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ',
     &  'Target_points_v_fib'
        STOP
      ENDIF ! iErr 
!%===============================================================

      ALLOCATE(Target_K_link_fsh(Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_K_link_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Target_t_link_fsh(Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_t_link_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(phi_fin_fsh(Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for phi_fin_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(a_fin_fsh(Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for a_fin_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(
     &      Target_points_fsh(Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_points_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &     force_points_fsh(Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for force_points_fsh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(
     &     Target_points_v_fsh(Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ',
     &  'Target_points_v_fsh'
        STOP
      ENDIF ! iErr 
!%===============================================================

      ALLOCATE(Target_K_link_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_K_link_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Target_t_link_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_t_link_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(phi_fin_esh(Nr_IBM_esh,0:Ns_IBM_esh+1,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for phi_fin_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(a_fin_esh(Nr_IBM_esh,0:Ns_IBM_esh+1,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for a_fin_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(Target_points_esh(Nr_IBM_esh,0:Ns_IBM_esh+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_points_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(force_points_esh(Nr_IBM_esh,0:Ns_IBM_esh+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for force_points_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Target_points_v_esh(Nr_IBM_esh,0:Ns_IBM_esh+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ',
     &  'Target_points_v_esh'
        STOP
      ENDIF ! iErr 
!%===============================================================

      ALLOCATE(Materialpara_esh
     &          (Nr_IBM_esh,0:Ns_IBM_esh+1,1:Ns_IBM_Fbr_esh)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Materialpara_esh'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Target_num(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_num'
        STOP
      ENDIF ! iErr 

      ALLOCATE(implicitflag(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for implicitflag'
        STOP
      ENDIF ! iErr

      ALLOCATE(FlagFixed_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FlagFixed_fib'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(FlagFixed_fsh(Nr_IBM_fsh,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FlagFixed_fsh'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(FlagFixed_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FlagFixed_esh'
        STOP
      ENDIF ! iErr 
!%===============================================================

      ALLOCATE( 
     &  Target_point_num_fib(Nr_IBM_fib,0:Ns_IBM_fib+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for '
     &  ,'Target_point_num_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE( 
     &  Target_point_num_fsh(Nr_IBM_fsh,2,0:Ns_IBM_fsh*Nq_IBM_fsh+1)
     &  ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for '
     &  ,'Target_point_num_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE( 
     &  Target_point_num_esh(Nr_IBM_esh,0:Ns_IBM_esh+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for '
     &  ,'Target_point_num_esh'
        STOP
      ENDIF ! iErr 

C Bucket Search
      ALLOCATE(NBucket(Nr_IBM,NMaxnumbuck),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for NBucket'
        STOP
      ENDIF ! iErr

      ALLOCATE(PBucket(Nr_IBM,NMaxnumbuck),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for PBucket'
        STOP
      ENDIF ! iErr

      ALLOCATE(IBucket_fib(Nr_IBM_fib,Ns_IBM_fib),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for IBucket_fib'
        STOP
      ENDIF ! iErr

      ALLOCATE(SBucket_fib(Nr_IBM_fib,NMaxnumbuck,Ns_IBM_fib)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for SBucket_fib'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(IBucket_esh(Nr_IBM_esh,Ns_IBM_esh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 

     &  'Allocate_memory: Memory Allocation Error for IBucket_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(SBucket_esh(Nr_IBM_esh,NMaxnumbuck,Ns_IBM_esh)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for SBucket_esh'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(IBucket_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for IBucket_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(
     & SBucket_fsh(Nr_IBM_fsh,NMaxnumbuck,Ns_IBM_fsh*Nq_IBM_fsh,2)
     &            ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for SBucket_fsh'
        STOP
      ENDIF ! iErr 
!%===============================================================

!++++++++++++++++++++ FOR  TEST  OF  AIRBAGE INFLATION : START ++++++++++++++++++++++++++++

      ALLOCATE(Pnormal_Ext(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Pnormal_Ext'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(n_TDpt_Pnormal_Ext(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for n_TDpt_Pnormal_Ext'
        STOP
      ENDIF ! iErr

      ALLOCATE(flag_TDpt_Pnormal_Ext(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &   'Allocate_memory: Memory Allocation Error for :CONT'
        WRITE(*,*) 
     &   'CONT: flag_TDpt_Pnormal_Ext '
        STOP
      ENDIF ! iErr

      ALLOCATE(data_TDpt_Pnormal_Ext(Nr_IBM,1000,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for :CONT'
        WRITE(*,*) 
     &   'CONT:data_TDpt_Pnormal_Ext'
        STOP
      ENDIF ! iErr 

C piezo paramters
      ALLOCATE(piezo_beta(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for piezo_beta'
        STOP
      ENDIF ! iErr

      ALLOCATE(piezo_gamma(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for piezo_gamma'
        STOP
      ENDIF ! iErr

      ALLOCATE(piezo_damp(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for piezo_damp'
        STOP
      ENDIF ! iErr

      ALLOCATE(piezo_Coef(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for piezo_Coef'
        STOP
      ENDIF ! iErr

      ALLOCATE(piezo_alpha11(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for piezo_alpha11'
        STOP
      ENDIF ! iErr

      ALLOCATE(piezo_alpha22(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for piezo_alpha22'
        STOP
      ENDIF ! iErr

      ALLOCATE(piezo_alpha12(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for piezo_alpha12'
        STOP
      ENDIF ! iErr   
  
      ALLOCATE(piezo_Coef_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh,-1:5)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for piezo_Coef_fsh'
        STOP
      ENDIF ! iErr   !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh,-1:5)

      ALLOCATE(Fpiezo1_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Fpiezo1_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(Fpiezo2_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Fpiezo2_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(Fpiezo3_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Fpiezo3_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(ndot_ibm_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ndot_ibm_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(ndot_ibm_fsh0(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ndot_ibm_fsh0'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(ndotold_ibm_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ndotold_ibm_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(ndotpre_ibm_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ndotpre_ibm_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(vPiezo_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for vPiezo_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(dvPiezo_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for dvPiezo_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(ddvPiezo_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ddvPiezo_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(vPiezoold_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for vPiezoold_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(vPiezopre_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for vPiezopre_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(TzeroB_fsh(Nr_IBM_fsh,4,Nq_IBM_fsh+Ns_IBM_fsh,3)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TzeroB_fsh'
        STOP
      ENDIF ! iErr     !(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)

      ALLOCATE(e4coef(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN

        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for e4coef'
        STOP
      ENDIF ! iErr     !(Nr_IBM)

      ALLOCATE(nibmptrecord(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for nibmptrecord'
        STOP
      ENDIF ! iErr     !(Nr_IBM)

      ALLOCATE(ibmptrecord(Nr_IBM,100,2)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ibmptrecord'
        STOP
      ENDIF ! iErr    !(Nr_IBM,100,2)


!%===============================================================

      ALLOCATE(Ncontactplanes(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ncontactplanes'
        STOP
      ENDIF ! iErr     !(Nr_IBM)

      ALLOCATE(contactplanesFlag(Nr_IBM)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for contactplanesFlag'
        STOP
      ENDIF ! iErr     !(Nr_IBM)

      ALLOCATE(contactplanesPara(Nr_IBM,NPlanesIBM,6)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for contactplanesPara'
        STOP
      ENDIF ! iErr    !(Nr_IBM,NPlanes,6)

      ALLOCATE(numedge(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for numedge '
        STOP
      ENDIF ! iErr

      ALLOCATE(EdgeList_fbc(nr_ibm_fbc,numedgeMAX_fbc,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for EdgeList_fbc '
        STOP
      ENDIF ! iErr

      ALLOCATE(InplaneLinemode_fbc(nr_ibm_fbc,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &'Allocate_memory: Memory Allocation Error for InplaneLinemode_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(numeFacePair(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for numeFacePair '
        STOP
      ENDIF ! iErr

      ALLOCATE(inplane_type_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for inplane_type_fbc '
        STOP
      ENDIF ! iErr

      ALLOCATE(facepairEle_fbc(nr_ibm_fbc,numeFacePairMAX_fbc,2)
     &   ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for facepairEle_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(facepairNode_fbc(nr_ibm_fbc,numeFacePairMAX_fbc,4)
     &   ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for facepairNode_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(EdgeLength0_fbc(nr_ibm_fbc,numedgeMAX_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for EdgeLength0_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(EleCenter0_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for EleCenter0_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(EleNormal0_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for EleNormal0_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(EleEdge0_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for EleEdge0_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(EleInplaneParaK_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3) 
     &  ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for' 
     &  ,' EleInplaneParaK_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(EleInplaneParaC_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3) 
     &  ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for' 
     &  ,' EleInplaneParaC_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(EleInplane_Cq_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for EleInplane_Cq_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(EleArea0_fbc(nr_ibm_fbc,Ns_IBMEle_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for EleArea0_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(areapnt_fbc(nr_ibm_fbc,ns_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for areapnt_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(kappa_volume_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for kappa_volume_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(KBT_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for kappa_volume_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(volume_eq_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for volume_eq_fbc '
        STOP
      ENDIF ! iErr

      ALLOCATE(kappa_Area_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for kappa_Area_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Area_eq_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Area_eq_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(kappa_inplane_fbc(nr_ibm_fbc,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for kappa_inplane_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(edge_length_p_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for edge_length_p_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(edge_length_eq_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &'Allocate_memory: Memory Allocation Error for edge_length_eq_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(edge_length_max_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &'Allocate_memory: Memory Allocation Error for edge_length_max_fbc'
        STOP
      ENDIF ! iErr


      ALLOCATE(InplaneLinem_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for InplaneLinem_fbc'
        STOP
      ENDIF ! iErr  
      ALLOCATE(inplane_Cq_fbc(nr_ibm_fbc,4),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for inplane_Cq_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(kappa_Bend_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for kappa_Bend_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(facepairtet0_fbc(nr_ibm_fbc,numeFacePairMAX_fbc)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for facepairtet0_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(vertex_mass_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for vertex_mass_fbc'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(inplane_E_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for inplane_E_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(inplane_nu_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for inplane_nu_fbc'
        STOP
      ENDIF ! iErr 
!logical:
       ALLOCATE(volume_flag_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for volume_flag_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Aera_flag_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Aera_flag_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(inplane_flag_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for inplane_flag_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Bending_flag_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Bending_flag_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(tet0_fbc(nr_ibm_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for tet0_fbc'
        STOP
      ENDIF ! iErr 
!somehow related variables
!%===============================================================
      ALLOCATE(GX_IBM_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBM_MASSIVE_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVE_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_MASSIVE_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVE_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_MASSIVE_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVE_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBM_MASSIVEo_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVEo_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBM_MASSIVEo_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVEo_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBM_MASSIVEo_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVEo_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBMo1_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBMo1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBMo1_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBMo1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBMo1_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBMo1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_IBMpre_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_IBMpre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_IBMpre_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_IBMpre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_IBMpre_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_IBMpre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GX_IBM_MASSIVEpre_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GX_IBM_MASSIVEpre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GY_IBM_MASSIVEpre_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GY_IBM_MASSIVEpre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &GZ_IBM_MASSIVEpre_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'GZ_IBM_MASSIVEpre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(DS_IBM_fbc( 
     &Nr_IBM_fbc,1-Ns_IBMB_fbc:Ns_IBM_fbc+Ns_IBMB_fbc 
     &),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for DS_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM1_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM2_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM2_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIBM3_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIBM3_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FS_1_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_1_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FS_2_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_2_IBM_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(FS_1_IBMo_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FS_1_IBMo_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_1_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_1_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_2_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_2_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FB_3_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FB_3_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_1_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_1_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_2_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_2_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FIner_3_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FIner_3_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK1_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK2_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN

        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK2_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_LINK3_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_LINK3_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls1_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls2_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls2_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(F_impuls3_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for F_impuls3_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(RIBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for RIBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM1_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM2_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM2_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM3_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM3_fbc'
        STOP
      ENDIF ! iErr 
 
      ALLOCATE(AIBM1_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(AIBM2_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM2_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(AIBM3_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for AIBM3_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(VIBM1_pre_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM1_pre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM2_pre_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM2_pre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(VIBM3_pre_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for VIBM3_pre_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(MASS_IBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for MASS_IBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Q_MASSo_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Q_MASSo_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Q_MASS_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Q_MASS_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GX_BP_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_BP_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_BP_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_BP_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_BP_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_BP_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(GX_BP_fbc0(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GX_BP_fbc0'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GY_BP_fbc0(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GY_BP_fbc0'
        STOP
      ENDIF ! iErr 

      ALLOCATE(GZ_BP_fbc0(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for GZ_BP_fbc0'
        STOP
      ENDIF ! iErr

      ALLOCATE(target_kvalue_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for target_kvalue_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(
     & contact_coef_fbc(Nr_IBM_fbc,NPlanesIBM,0:Ns_IBM_fbc+1),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for contact_coef_fbc'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(UIBM1_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(UIBM2_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM2_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(UIBM3_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for UIBM3_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS1_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS1_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS2_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS2_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS3_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS3_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FK_MASS1o_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS1o_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(FK_MASS2o_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS2o_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(FK_MASS3o_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FK_MASS3o_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(ffluidsum1_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum1_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(ffluidsum2_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum2_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(ffluidsum3_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ffluidsum3_fbc'
        STOP
      ENDIF ! iErr 
!%===============================================================
      ALLOCATE(Jreal_2_imag_fbc(Nr_IBM_fbc,Ns_IBM_fbc,Ns_IBM_fbc)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Jreal_2_imag_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(xpcenter_fbc (Nr_IBM_fbc,Ns_IBMEle_fbc,3)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for xpcenter_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Nq_IBM_r_fbc(Nr_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Nq_IBM_r_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ns_IBM_r_fbc(Nr_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_r_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(NS_IBM_i_fbc(Nr_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for NS_IBM_i_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ns_IBM_rall_fbc(Nr_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ns_IBM_rall_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Iglbloc_fbc(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Iglbloc_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Ilocglb_fbc(Nr_IBM),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Ilocglb_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(n_matpara_fbc(Nr_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN

        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for n_matpara_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(BndyCnd_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for BndyCnd_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(target_ktype_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for target_ktype_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(BoundryList_fbc(Nr_IBM_fbc,Ns_IBM_img_fbc,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for BoundryList_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(nea_fbc(Nr_IBM_fbc,Ns_IBMEle_fbc,maxordern),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for nea_fbc'
        STOP
      ENDIF ! iErr 

        Ns_IBMmax_fbc=max(Ns_IBMEle_fbc,Ns_IBM_fbc)

      ALLOCATE(
     &          ContactShellFlag_fbc
     &          (Nr_IBM_fbc,Ns_IBMmax_fbc,Ns_IBMmax_fbc)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ', 
     &          'ContactShellFlag_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(QIBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for QIBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(TFIBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TFIBM_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(TSIBM_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for TSIBM_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(Target_K_link_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_K_link_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Target_t_link_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_t_link_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(phi_fin_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for phi_fin_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(a_fin_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1,3),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for a_fin_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Target_points_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for Target_points_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(force_points_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for force_points_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(Target_points_v_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1,3)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for ',
     &  'Target_points_v_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FlagFixed_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for FlagFixed_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE( 
     &  Target_point_num_fbc(Nr_IBM_fbc,0:Ns_IBM_fbc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for '
     &  ,'Target_point_num_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(IBucket_fbc(Nr_IBM_fbc,Ns_IBM_fbc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 

     &  'Allocate_memory: Memory Allocation Error for IBucket_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(SBucket_fbc(Nr_IBM_fbc,NMaxnumbuck,Ns_IBM_fbc)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for SBucket_fbc'
        STOP
      ENDIF ! iErr 

      ALLOCATE(FabricMeshFlag_fbc(Nr_IBM_fbc)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for FabricMeshFlag_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(FabricLineFlag_fbc(Nr_IBM_fbc)
     &          ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for FabricLineFlag_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(contact_time_fsh(Nr_IBM_fsh,2),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for contact_time_fsh'
        STOP
      ENDIF ! iErr

      ALLOCATE(contact_time_esh(Nr_IBM_esh,2),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for contact_time_esh'
        STOP
      ENDIF ! iErr

      ALLOCATE(contact_time_fbc(Nr_IBM_fbc,2),STAT=iErr)

      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for contact_time_fbc'
        STOP
      ENDIF ! iErr

      ALLOCATE(flag_fractal_esh(Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for flag_fractal_esh'
        STOP
      ENDIF ! iErr    !(Nr_IBM_esh)
      endif
!%===============================================================

      if((iFlagin .eq. 1) .and.(FirstCall_Frac)) then
      ALLOCATE(nifractal_esh(Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     &  'Allocate_memory: Memory Allocation Error for nifractal_esh'
        STOP
      ENDIF ! iErr    !(Nr_IBM_esh)

      ALLOCATE(fractalinitFlag_esh(Nr_IBM_esh)
     &         ,STAT=iErr)
      fractalinitFlag_esh=.TRUE. !INITIALIZATION
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error ',
     & 'for fractalinitFlag_esh'
        STOP
      ENDIF ! iErr    !(Nr_IBM_esh)


      ALLOCATE(fractalpara_esh(Nr_IBM_esh,9)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for fractalpara_esh'
        STOP
      ENDIF ! iErr    !(Nr_IBM_esh,9)

      ALLOCATE(Dfractal_esh(3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Dfractal_esh'
        STOP
      ENDIF ! iErr    !(3,Nr_IBM_esh)

      ALLOCATE(Jfractal2_esh(3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Jfractal2_esh'
        STOP
      ENDIF ! iErr    !(3,Nr_IBM_esh)


      ALLOCATE(Dfractal2_esh(12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Dfractal2_esh'
        STOP
      ENDIF ! iErr    !(12,3,Nr_IBM_esh)

      ALLOCATE(Hfractal2_esh(12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Hfractal2_esh'
        STOP
      ENDIF ! iErr    !(12,3,Nr_IBM_esh)

      ALLOCATE(Ifractal2_esh(12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Ifractal2_esh'
        STOP
      ENDIF ! iErr    !(12,3,Nr_IBM_esh)

      ALLOCATE(Bfractal_esh(12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Bfractal_esh'
        STOP
      ENDIF ! iErr    !(12,3,Nr_IBM_esh)

      ALLOCATE(Cfractal_esh(12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Cfractal_esh'
        STOP
      ENDIF ! iErr    !(12,3,Nr_IBM_esh)

      ALLOCATE(fibmsig_esh(Nr_IBM_esh,3,Ns_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for fibmsig_esh'
        STOP
      ENDIF ! iErr    !(Nr_IBM_esh,3,Ns_IBM_esh)

      ALLOCATE(Afractal_esh(12,12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Afractal_esh'
        STOP
      ENDIF ! iErr    !(12,12,3,Nr_IBM_esh)

      ALLOCATE(Afractal2_esh(12,12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Afractal2_esh'
        STOP
      ENDIF ! iErr    !(12,12,3,Nr_IBM_esh)

      ALLOCATE(Bfractal2_esh(12,12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Bfractal2_esh'
        STOP
      ENDIF ! iErr    !(12,12,3,Nr_IBM_esh)

      ALLOCATE(Cfractal2_esh(12,12,3,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Cfractal2_esh'
        STOP
      ENDIF ! iErr    !(12,12,3,Nr_IBM_esh)

          !
      ALLOCATE(Xfractal_esh(72,ngaumax,Ns_IBMEle_esh,Nr_IBM_esh)
     &         ,STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(*,*) 
     & 'Allocate_memory: Memory Allocation Error for Xfractal_esh'
        STOP
      ENDIF ! iErr    !(72,ngaumax,Ns_IBMEle_esh,Nr_IBM_esh)

      FirstCall_Frac = .FALSE.
      endif
      END SUBROUTINE 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

