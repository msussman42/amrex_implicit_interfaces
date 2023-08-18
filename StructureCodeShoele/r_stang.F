c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine r_stang(ibm_solve,monitorON,isavesigma)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c solve for 2d and 3d structures
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      implicit none
      integer ibdy,np,npg,bodypoint(3),ne,jtype,itype,nos,ntem,noi,noj
     &     ,ntemi, ntemj

      real*8 wei_i,iner1,iner2,jner1,jner2
      integer ntype,info,isavesigma
      logical monitorON 
      integer ibm_solve
      real*8 mms(ns_ibm_esh,ns_ibm_esh)

      real*8 x(maxordern,3), sx(2,maxordern), mass(maxordern,maxordern)
     &    ,ks(maxordern,maxordern), kb(maxordern,maxordern)
     &    ,xpos(3),dxpos(2,3),d2xpos(3,3)
     &    ,ms(maxordern,maxordern)
     &    ,xj1,xj2,xj3,yj1,yj2,yj3,zj1,zj2,zj3

       real*8 s(1,maxordern), dsdv(2,maxordern),ddsdvdv(3,maxordern)
     &    , dmat(2,2),kmat(2,2)
     &    , kmatbar(2,2) , xi1, xi2, smat(2,2), bmat(3,3), jmat(2,2)
     &    , gmat(3,2), hmat(3,3)
     &    , amat1(2,maxordern), amat2(3,maxordern)
     &    , amat1T(maxordern,2), amat2T(maxordern,3)
     &    ,FIBM1pres_fib(Ns_IBM_fib,2)    
     &    ,FIBM2pres_fib(Ns_IBM_fib,2)    
     &    ,FIBM3pres_fib(Ns_IBM_fib,2)
     &    ,FIBM1pres_esh(Ns_IBM_esh)    
     &    ,FIBM2pres_esh(Ns_IBM_esh)    
     &    ,FIBM3pres_esh(Ns_IBM_esh)
     &  
     &    ,FIBM1cont_esh(Ns_IBM_esh)    
     &    ,FIBM2cont_esh(Ns_IBM_esh)    
     &    ,FIBM3cont_esh(Ns_IBM_esh)   
     &  
     &    ,normalpres_esh(Ns_IBM_esh,3) 
     &    ,FIBM1pres_fsh(nq_ibm_fsh,ns_ibm_fsh)    
     &    ,FIBM2pres_fsh(nq_ibm_fsh,ns_ibm_fsh)    
     &    ,FIBM3pres_fsh(nq_ibm_fsh,ns_ibm_fsh)  

      real*8 
     &     FIBM1pres_fbc(Ns_IBM_fbc)    
     &    ,FIBM2pres_fbc(Ns_IBM_fbc)    
     &    ,FIBM3pres_fbc(Ns_IBM_fbc)
     &  
     &    ,FIBM1cont_fbc(Ns_IBM_fbc)    
     &    ,FIBM2cont_fbc(Ns_IBM_fbc)    
     &    ,FIBM3cont_fbc(Ns_IBM_fbc)   

      integer i,j,k,k2,para,t_temp,n,j2,kcontact,jq,jq2
      integer ierrksh1,iksh1,iksh2,iksh3,i2,iksh1b,iksh1c

      real*8 aamat_fib(ns_ibm_fib,ns_ibm_fib)   !tobe changes
      real*8 aamat_total(nr_ibm_mas,ns_ibm_esh,ns_ibm_esh)    
     &    ,bbmat_total(nr_ibm_mas,ns_ibm_esh,3)   !tobe changes

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsr(:),iacsr(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsr(:)

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsrm(:),iacsrm(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsrm(:)

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsrfibershell(:),iacsrfibershell(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsrfibershell(:),vecRHSfibershell(:),
     &             vecXfibershell(:)

! contact parameters
      REAL*8 gau_i(3,ngaumax)
     &      ,XP1cont(3),XP2cont(3),xcont(maxordern,3),dxcont(3)
     &      ,ddcont,dncont,Fcont
     &      ,dXP1cont(2,3),Normalcont1(3),outcont(3),outcont2(3)
     &      ,Normalcont2(3),ContactForce
      REAL*8,SAVE, ALLOCATABLE :: 
     &             l_loc(:,:),iv_loc(:,:,:)
     &             ,phi_loc(:,:,:,:)
      INTEGER ntypeMAx_loc
     &        ,ItypecontELE1(3),ItypecontELE(3)
     &        ,ne2,ntype2,itype2,nos2,jtype2,ntem2

      INTEGER, SAVE :: ntermiMAX,ntermiMAXfibershell,ntermiMAXmaster 
      INTEGER, SAVE :: ntermiMAXm(nr_ibm_mas) 
!      real aamcsr(nr_ibm_mas,ns_ibm*20)
!      integer jamcsr(nr_ibm_mas,ns_ibm*20),iamcsr(nr_ibm_mas,ns_ibm)

      real temp1, temp2, temp3, temp4,ksh_tmp
      real temp_ibm1,temp_ibm2,temp_ibm3
     &      ,temp_ibm1b,temp_ibm2b,temp_ibm3b
     &      ,fn1,fn2,fn3
     &      ,temp_check1,temp_check2,temp_check3, temp2_ibm, temp3_ibm
      real   temp4_1_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)    
     &      , temp4_2_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)    
     &      , temp4_3_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)     
     &      ,a_fib(ns_ibm_fib,-3:3) 
     &      ,mat(4,4)
     &      ,tdx0_1_fib(nr_ibm_fib)    
     &      ,tdx0_2_fib(nr_ibm_fib)    
     &      ,tdx0_3_fib(nr_ibm_fib)    
     &      ,tdxnm1_1_fib(nr_ibm_fib)    
     &      ,tdxnm1_2_fib(nr_ibm_fib)    
     &      ,tdxnm1_3_fib(nr_ibm_fib)    

      real*8  s_function
      logical fix_bc,free_bc,fix_bcb,free_bcb 
      integer itmax,precond,success,itconverge
      real*8 tol,h_channel,KB_ibmt
    
      real*8 dt,dt2,time_alphaf 

      real*8 tmpv1(3),tmpv2(3),tmpv3(3),tmpv1x(3)
      real*8 areapnt_esh(nr_ibm_esh,ns_ibm_esh)
!     &       ,xpnt_esh(nr_ibm_esh,ns_ibm_esh)  !for inflation problem 

      integer jbndy,jqbndy

      real*8  X_i_j(3),X_ip_j(3),X_i_jp(3),X_ipf_j(3)
     &       ,X_imf_j(3),X_i_jpf(3),X_i_jmf(3)
     &       ,X_ipf_jp(3),X_imf_jp(3),X_ip_jpf(3)
     &       ,X_ip_jmf(3),X_ipf_jpf(3),X_imf_jpf(3)
     &       ,X_ipf_jmf(3),X_imf_jmf(3)
     &       ,DXX_i_j(3),DYY_i_j(3),DXY_i_j(3),DYX_i_j(3)
     &       ,DX_im_jm(3),DX_im_j(3),DX_im_jp(3)
     &       ,DX_i_jm(3),DX_i_j(3),DX_i_jp(3)
     &       ,DX_ip_jm(3),DX_ip_j(3),DX_ip_jp(3)
     &       ,DY_im_jm(3),DY_im_j(3),DY_im_jp(3)
     &       ,DY_i_jm(3),DY_i_j(3),DY_i_jp(3)
     &       ,DY_ip_jm(3),DY_ip_j(3),DY_ip_jp(3)
     &       ,kmem_i_j(50,3),kbend_i_j
     &       ,GXtmp(1:3,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &              ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh) 
     &       ,dGXtmp(1:3,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &              ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh) 
     &        ,Siga(nq_ibm_fsh,ns_ibm_fsh,2,2) 
     &        ,SigB(4,nq_ibm_fsh+ns_ibm_fsh,3) 
     &       ,gama(nq_ibm_fsh,ns_ibm_fsh,2,2) 

      integer ibkt,ibkt2,ibkt_i,ibkt2_i
     &       ,ibktmp1,ibktx,ibkty,ibktz,ibkts1,ibkts2,ibkts3
     &       ,ibktmin(3),ibktmax(3)

      real*8  dvcont(3),dnvcont,Vcont(maxordern,3),xvv(maxordern,3)
     &       ,VP1cont(3),VP2cont(3)
      integer i2l,i2g, ibdyg,i3g

      real*8  Stiff_points_fsh,target_timefn
      REAL*8,SAVE, ALLOCATABLE :: Stiff_link_fsh(:,:)
      REAL*8,SAVE, ALLOCATABLE :: vnPiezo_fsh(:,:,:),aamat_fsh(:,:,:)
      REAL*8,SAVE, ALLOCATABLE :: VECpiezo_tmp(:,:,:,:)
      REAL*8,SAVE, ALLOCATABLE :: ForceOut(:)
!temporal variables
      real*8  real_distibm(50),rtmpusrfn(50)
      integer int_distibm(50),itmpusrfn(50)
      integer jq_pm(2),j_pm(2)
      real*8  contactplane_n(4),contactplane_point(3)
      real*8  matFib_Esh(3,3)
      real*8  R_Fibmat(3,3),xposFib(3),x_FibVec(maxordern,3)
     &       ,n_FibVec(3),n3_FibVec(3),temp_piezo

!variable for genralized alpha
      real*8 X_n_esh(ns_ibm_esh,3)
     &       ,dX_n_esh(ns_ibm_esh,3)
     &       ,ddX_n_esh(ns_ibm_esh,3)
     &       ,X_np1_esh(ns_ibm_esh,3)
     &       ,dX_np1_esh(ns_ibm_esh,3)
     &       ,ddX_np1_esh(ns_ibm_esh,3)
     &       ,X_na_esh(ns_ibm_esh,3)
     &       ,dX_na_esh(ns_ibm_esh,3)
     &       ,ddX_na_esh(ns_ibm_esh,3)
     &       ,Xtar_na_esh(ns_ibm_esh,3)
     &       ,dXtar_na_esh(ns_ibm_esh,3)
     &       ,X0_np1_esh(ns_ibm_esh,3)
     &       ,ddU_np1_esh(ns_ibm_esh,3)
     &       ,dX_in_esh(ns_ibm_esh,3)
     &       ,ddX_in_esh(ns_ibm_esh,3)
     &       ,timeMarchcoef(20)
     &       ,d_smat(3,3)
     &       ,mattmp_eye(3,3)
     &       ,coef_iter

      real*8 X_n_fbc(ns_ibm_fbc,3)
     &       ,dX_n_fbc(ns_ibm_fbc,3)
     &       ,ddX_n_fbc(ns_ibm_fbc,3)
     &       ,X_np1_fbc(ns_ibm_fbc,3)
     &       ,dX_np1_fbc(ns_ibm_fbc,3)
     &       ,ddX_np1_fbc(ns_ibm_fbc,3)
     &       ,X_na_fbc(ns_ibm_fbc,3)
     &       ,dX_na_fbc(ns_ibm_fbc,3)
     &       ,ddX_na_fbc(ns_ibm_fbc,3)
     &       ,Xtar_na_fbc(ns_ibm_fbc,3)
     &       ,dXtar_na_fbc(ns_ibm_fbc,3)
     &       ,X0_np1_fbc(ns_ibm_fbc,3)
     &       ,ddU_np1_fbc(ns_ibm_fbc,3)
     &       ,dX_in_fbc(ns_ibm_fbc,3)
     &       ,ddX_in_fbc(ns_ibm_fbc,3)

      real*8 X_n_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,dX_n_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddX_n_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,X_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,dX_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddX_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,X_na_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,dX_na_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddX_na_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,Xtar_na_fsh(nq_ibm_fsh*ns_ibm_fsh,3)
     &       ,dXtar_na_fsh(nq_ibm_fsh*ns_ibm_fsh,3)
     &       ,X0_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddU_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,dX_in_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddX_in_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,Siga_d(nq_ibm_fsh,ns_ibm_fsh,5)
     &       ,SigB_d(4,nq_ibm_fsh+ns_ibm_fsh,3)
     &       ,gama_d(nq_ibm_fsh,ns_ibm_fsh,5)



      real*8 XvPiezo_n_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,XvPiezo_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,XvPiezo_na_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,XvPiezo_in_fsh(nq_ibm_fsh,ns_ibm_fsh,3)

      integer iter

      real*8, ALLOCATABLE ::
     &     aamat_esh(:,:)
!     (ns_ibm_esh or (3*ns_ibm_esh,
     &     ,aamat_fbc(:,:)
!     (ns_ibm_fbc or (3*ns_ibm_fbc,
     &     ,matdir1(:,:)
!     (3,maxordern)
     &     ,matdir2(:,:)
!     (3,maxordern)
     &     ,matdir3(:,:)
!     (3,maxordern)
     &     ,Dmat_Stiff1(:,:)
!     (maxordern,3)
     &     ,Dmat_Stiff2(:,:)
!     (maxordern,3)
     &     ,Dmat_Stiff3(:,:)
!     (maxordern,3)
     &     ,kb11(:,:)
!     (maxordern,maxordern)
     &     ,kb12(:,:)
!     (maxordern,maxordern)
     &     ,kb13(:,:)
!     (maxordern,maxordern)
     &     ,kb21(:,:)
!     (maxordern,maxordern)
     &     ,kb22(:,:)
!     (maxordern,maxordern)
     &     ,kb23(:,:)
!     (maxordern,maxordern)
     &     ,kb31(:,:)
!     (maxordern,maxordern)
     &     ,kb32(:,:)
!     (maxordern,maxordern)
     &     ,kb33(:,:)
!     (maxordern,maxordern)
     &     ,kks(:,:)
!     (ns_ibm_esh or (3*ns_ibm_esh,
     &     ,kkb(:,:)
!     (ns_ibm_esh or (3*ns_ibm_esh,
     &     ,kksf(:,:)
!     (ns_ibm_esh or
     &     ,kkbf(:,:)
!     (ns_ibm_esh or
     &     ,kkb_fbc(:,:)
!     (ns_ibm_fbc or (3*ns_ibm_fbc,
     &      ,mms_fbc(:,:) 

      real*8 smatsave_esh(3,ngaumax,Ns_IBMEle_esh)

      real*8 X22piezo(1:3),X2piezo(1:3),X11piezo(1:3),X1piezo(1:3)
     &      ,X12piezo(1:3),normalpiezo(1:3)
     &      ,Epiezo,Fpiezo,Gpiezo,Lpiezo,Mpiezo,Npiezo,Hnpiezo,Hdpiezo
     &      ,Ppiezo,Rpiezo
     &      ,VectempPiezo(3),Qmatpiezo(1:3,1:3),Pmatpiezo(1:3,1:3)

      logical ilucal_flag
      integer, save :: dist_esh(ns_ibm_esh)

      integer in_fbc
     &      ,iv0,iv1,iv2,iv3
     &      ,ie1,ie2

      real*8 EdgeLength_fbc(nr_ibm_fbc,numedgeMAX_fbc)
     &      ,EleCenter_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3)
     &      ,EleNormal_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3)
     &      ,EleEdge_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3)
     &      ,EleArea_fbc(nr_ibm_fbc,Ns_IBMEle_fbc)
     &      ,axtmp_fbc,aytmp_fbc,aztmp_fbc
     &      ,bxtmp_fbc,bytmp_fbc,bztmp_fbc
     &      ,cxtmp_fbc,cytmp_fbc,cztmp_fbc
     &      ,tmpcoord(3),tmpcoord2(3),tmpcoord3(3)
     &      ,total_volume_fbc(nr_ibm_fbc)
     &      ,total_Area_fbc(nr_ibm_fbc)
     &      ,cf_coef,cf_coef2,cf_coef1
     &      ,x0,x1,x2,x3
     &      ,y0,y1,y2,y3
     &      ,z0,z1,z2,z3
     &      ,Forcetmpfbc(4,3)
     &      ,Ktmp(4,3,4,3)
     &      ,L_over_Lmax
     &      ,sinsign
     &      ,costet0,sintet0

      INTEGER iDX0,iKS0,N_iPARA0,iPARA0(100),i_FABRICLINE0

      real*8 FabricML_dx(3)
     &      ,FabricML_att(FabricMeshL_numattMax)
     &      ,outcount2(3),outcount(3) 
     &      ,FabricML_tension,FabricML_Dtension
     &      ,FabricLine_curv(numLinepointRecordMax,3)

      real*8, dimension(numLinepointRecordMax,3) ::
     &       FabricLine_X_n,FabricLine_dX_n,FabricLine_ddX_n
     &      ,FabricLine_X_np1,FabricLine_dX_np1,FabricLine_ddX_np1
     &      ,FabricLine_ddX_in,FabricLine_dX_in
     &      ,FabricLine_X_na,FabricLine_dX_na,FabricLine_ddX_na
     &      ,FabricLine_force_points
     &      ,FabricLine_Xtar_np1,FabricLine_dXtar_np1,FabricLine_X0_np1
     &      ,FabricLine_Xtar_na,FabricLine_dXtar_na

      real*8, dimension(numMeshpointRecordMax,3) ::
     &       FabricMesh_X_n,FabricMesh_dX_n,FabricMesh_ddX_n
     &      ,FabricMesh_X_np1,FabricMesh_dX_np1,FabricMesh_ddX_np1
     &      ,FabricMesh_ddX_in,FabricMesh_dX_in
     &      ,FabricMesh_X_na,FabricMesh_dX_na,FabricMesh_ddX_na
     &      ,FabricMesh_force_points
     &      ,FabricMesh_Xtar_np1,FabricMesh_dXtar_np1,FabricMesh_X0_np1
     &      ,FabricMesh_Xtar_na,FabricMesh_dXtar_na

      REAL*8,  ALLOCATABLE :: 
     &             vecRHS(:,:)
     &            ,x0LHS(:,:)

      Logical tensionDerivFlag_fbc

      real*8 maxtmpx,normalpres_fbc(nr_ibm_fbc,3)
      Logical KtmpFlag
!=====newallocatables for FE shell Continuum model
      integer  ih
      real*8 
     &    weih_i,theta3
     &   ,ametric(3),bmetric(3),dnxpos(3)
     &   ,gmetric_cov(2,2),gmetric_con0(2,2),gmetric_con(2,2)
     &   ,CauchyGreen_C(3,3),CauchyGreen_Cbar(3,3),C33_old
     &   ,detC_inplane,tracC_inplane
     &   ,detC,tracC
     &   ,DiffPhi_C(3),DiffPhi_C33,Diff2Phi_CC33(3),Diff2Phi_C(3,3)
     &   ,Diff2Phi_C33,Sstress33 
     &   ,Smetric(3),CMAT(3,3),CMAT33,dC33_old,CMATi3(3)
     &   ,nmatMetric(3),mmatMetric(3),DmatMetric(3,3,3)
     &   ,abar3inv
     &   ,avectmp0(3),avectmp1(3),avectmp2(3),avectmp3(3)
     &   ,avectmp1b(3),avectmp2b(3),avectmp3b(3)
     &   ,xposFib_con(3),q_contrav_bases0(2,3),MatposFib_con(3,3)

      real*8, allocatable, dimension(:)::
     &    fibmRHS1_esh,fibmRHS2_esh,fibmRHS3_esh
     &   ,temp_ibmtype,a3bar_Dr

      real*8, allocatable, dimension(:,:)::
     &    Forcetemp1
     &   ,temp_ibmtype2,temp_ibmtype3
     &   ,a3bar_Drs,a3_Dr
     &   ,Da_ab_Dr,a3tilde_Dr,Db_ab_Dr

      real*8, allocatable, dimension(:,:,:) ::
     &    Da_ab_Drs,a3tilde_Drs,a3_Drs,Depsilon_ab_r,Dkappa_ab_r

      real*8, allocatable, dimension(:,:,:,:)::
     &    Db_ab_Drs,Depsilon_ab_rs

      real*8, allocatable, dimension(:,:,:,:,:)::
     &    Dkappa_ab_rs

      real*8   dxcont_piezo,outcont_piezo
     &   ,alpha_coefK(5),Elect_RHS 
      integer  n_dim_aamatPiezo_fsh,n_aamatPiezo_fsh
      real*8, allocatable, dimension(:) ::
     &    aamatPiezo_fsh

      integer  FlagForceCalculation(nq_ibm_fsh,ns_ibm_fsh)
!==========================================================

      parameter (tol = 1.0d-11, precond=1, itmax=200)
       
      iksh1=simulationtype
!     iksh1=0; only fiber or only shell?  iksh1 ~=0  fiber + shell 

c time step
      dt=delta_t
      dt2=dt**2.0
      eps=1.0d-8

      if(ibm_solve .ne. 2) then

#ifdef INCLUDE_FIB
       if(nr_ibm_fib.ge. 1) then
        fibm1_fib=0.0
        fibm2_fib=0.0
        fibm3_fib=0.0
        ffluidsum1_fib=0.0
        ffluidsum2_fib=0.0
        ffluidsum3_fib=0.0
       endif
c
c calculating feedback force for matching bcs
       do i=1,nr_ibm_fib
        i2g=Ilocglb_fib(i)
        if(FluidForceFlag(i2g)) then
         do j=1,ns_ibm_r_fib(i)
           fk_mass1_fib(i,j)=fk_mass1o_fib(i,j)
           fk_mass2_fib(i,j)=fk_mass2o_fib(i,j)

           ffluidsum1_fib(i,j)= fk_mass1_fib(i,j)
           ffluidsum2_fib(i,j)= fk_mass2_fib(i,j)

           if(ndim .gt. 2) then
            fk_mass3_fib(i,j)=fk_mass3o_fib(i,j)
            ffluidsum3_fib(i,j)= fk_mass3_fib(i,j)
           endif

         end do
         do j=1,ns_ibm_r_fib(i)

              fibm1_fib(i,j)=fibm1_fib(i,j)+ffluidsum1_fib(i,j)
              fibm2_fib(i,j)=fibm2_fib(i,j)+ffluidsum2_fib(i,j)
              if(ndim .gt. 2) 
     &           fibm3_fib(i,j)=fibm3_fib(i,j)+ffluidsum3_fib(i,j)

         end do
        endif  
  
        do j=1,ns_ibm_r_fib(i)
           vibm1_pre_fib(i,j)=vibm1_fib(i,j)
           vibm2_pre_fib(i,j)=vibm2_fib(i,j)
           if(ndim .gt. 2) 
     &         vibm3_pre_fib(i,j)=vibm3_fib(i,j)
        end do
      end do  !nr_ibm_fib
#endif      
#ifdef INCLUDE_ESH
      if(nr_ibm_esh.ge. 1) then
      fibm1_esh=0.0
      fibm2_esh=0.0
      fibm3_esh=0.0
      ffluidsum1_esh=0.0
      ffluidsum2_esh=0.0
      ffluidsum3_esh=0.0
      endif
c
c calculating feedback force for matching bcs

      do i=1,nr_ibm_esh
       i2g=Ilocglb_esh(i)
       if(FluidForceFlag(i2g)) then
        do j=1,ns_ibm_r_esh(i)
           fk_mass1_esh(i,j)=fk_mass1o_esh(i,j)
           fk_mass2_esh(i,j)=fk_mass2o_esh(i,j)
           ffluidsum1_esh(i,j)= fk_mass1_esh(i,j)
           ffluidsum2_esh(i,j)= fk_mass2_esh(i,j)
           if(ndim .gt. 2) then
            fk_mass3_esh(i,j)=fk_mass3o_esh(i,j)
            ffluidsum3_esh(i,j)= fk_mass3_esh(i,j)
           endif
         end do

         do j=1,ns_ibm_r_esh(i)
              fibm1_esh(i,j)=fibm1_esh(i,j)+ffluidsum1_esh(i,j)
              fibm2_esh(i,j)=fibm2_esh(i,j)+ffluidsum2_esh(i,j)
              if(ndim .gt. 2) 
     &           fibm3_esh(i,j)=fibm3_esh(i,j)+ffluidsum3_esh(i,j)
          end do
        endif  
      end do  !nr_ibm_esh


      do i=1,nr_ibm_esh
        do j=1,ns_ibm_r_esh(i)
           vibm1_pre_esh(i,j)=vibm1_esh(i,j)
           vibm2_pre_esh(i,j)=vibm2_esh(i,j)
           if(ndim .gt. 2) 
     &         vibm3_pre_esh(i,j)=vibm3_esh(i,j)
        end do
      end do  !nr_ibm_esh

#endif
#ifdef INCLUDE_FSH
      if(nr_ibm_fsh.ge. 1) then
      fibm1_fsh=0.0
      fibm2_fsh=0.0
      fibm3_fsh=0.0
      ffluidsum1_fsh=0.0
      ffluidsum2_fsh=0.0
      ffluidsum3_fsh=0.0
      endif
c
c calculating feedback force for matching bcs
      do i=1,nr_ibm_fsh
       i2g=Ilocglb_fsh(i)
       if(FluidForceFlag(i2g)) then
       do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
           fk_mass1_fsh(i,jq,j)=fk_mass1o_fsh(i,jq,j)
           fk_mass2_fsh(i,jq,j)=fk_mass2o_fsh(i,jq,j)

           ffluidsum1_fsh(i,jq,j)= fk_mass1_fsh(i,jq,j)
           ffluidsum2_fsh(i,jq,j)= fk_mass2_fsh(i,jq,j)

           if(ndim .gt. 2) then
            fk_mass3_fsh(i,jq,j)=fk_mass3o_fsh(i,jq,j)
            ffluidsum3_fsh(i,jq,j)= fk_mass3_fsh(i,jq,j)
           endif

        end do
        end do

       do jq=1,nq_ibm_r_fsh(i)
         do j=1,ns_ibm_r_fsh(i)

              fibm1_fsh(i,jq,j)=fibm1_fsh(i,jq,j)+ffluidsum1_fsh(i,jq,j)
              fibm2_fsh(i,jq,j)=fibm2_fsh(i,jq,j)+ffluidsum2_fsh(i,jq,j)
              if(ndim .gt. 2) 
     &        fibm3_fsh(i,jq,j)=fibm3_fsh(i,jq,j)+ffluidsum3_fsh(i,jq,j)

          end do
        end do
        endif  
      end do  !nr_ibm_fsh
  
      do i=1,nr_ibm_fsh
       do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
           vibm1_pre_fsh(i,jq,j)=vibm1_fsh(i,jq,j)
           vibm2_pre_fsh(i,jq,j)=vibm2_fsh(i,jq,j)
           if(ndim .gt. 2) 
     &         vibm3_pre_fsh(i,jq,j)=vibm3_fsh(i,jq,j)
        end do
       enddo
      end do  !nr_ibm_fsh

#endif
#ifdef INCLUDE_FBC
      if(nr_ibm_fbc.ge. 1) then
      fibm1_fbc=0.0
      fibm2_fbc=0.0
      fibm3_fbc=0.0
      ffluidsum1_fbc=0.0
      ffluidsum2_fbc=0.0
      ffluidsum3_fbc=0.0
      endif
c
c calculating feedback force for matching bcs

      do i=1,nr_ibm_fbc
       i2g=Ilocglb_fbc(i)
       if(FluidForceFlag(i2g)) then
        do j=1,ns_ibm_r_fbc(i)
           fk_mass1_fbc(i,j)=fk_mass1o_fbc(i,j)
           fk_mass2_fbc(i,j)=fk_mass2o_fbc(i,j)
           ffluidsum1_fbc(i,j)= fk_mass1_fbc(i,j)
           ffluidsum2_fbc(i,j)= fk_mass2_fbc(i,j)
           if(ndim .gt. 2) then
            fk_mass3_fbc(i,j)=fk_mass3o_fbc(i,j)
            ffluidsum3_fbc(i,j)= fk_mass3_fbc(i,j)
           endif
         end do

         do j=1,ns_ibm_r_fbc(i)
              fibm1_fbc(i,j)=fibm1_fbc(i,j)+ffluidsum1_fbc(i,j)
              fibm2_fbc(i,j)=fibm2_fbc(i,j)+ffluidsum2_fbc(i,j)
              if(ndim .gt. 2) 
     &           fibm3_fbc(i,j)=fibm3_fbc(i,j)+ffluidsum3_fbc(i,j)
          end do
        endif  
      end do  !nr_ibm_fbc


      do i=1,nr_ibm_fbc
        do j=1,ns_ibm_r_fbc(i)
           vibm1_pre_fbc(i,j)=vibm1_fbc(i,j)
           vibm2_pre_fbc(i,j)=vibm2_fbc(i,j)
           if(ndim .gt. 2) 
     &         vibm3_pre_fbc(i,j)=vibm3_fbc(i,j)
        end do
      end do  !nr_ibm_fbc

#endif

401   format(3(i8,' '),20(f20.9,' '))
      else  ! ibm_solve ==2 and only structure
      if(nr_ibm_fib.ge. 1) then
      fibm1_fib=0.0
      fibm2_fib=0.0
      fibm3_fib=0.0
      ffluidsum1_fib=0.0
      ffluidsum2_fib=0.0
      ffluidsum3_fib=0.0
      endif

      if(nr_ibm_esh.ge. 1) then
      fibm1_esh=0.0
      fibm2_esh=0.0
      fibm3_esh=0.0
      ffluidsum1_esh=0.0
      ffluidsum2_esh=0.0
      ffluidsum3_esh=0.0
      endif

      if(nr_ibm_fsh.ge. 1) then
      fibm1_fsh=0.0
      fibm2_fsh=0.0
      fibm3_fsh=0.0
      ffluidsum1_fsh=0.0
      ffluidsum2_fsh=0.0
      ffluidsum3_fsh=0.0
      endif

      if(nr_ibm_fbc.ge. 1) then
      fibm1_fbc=0.0
      fibm2_fbc=0.0
      fibm3_fbc=0.0
      ffluidsum1_fbc=0.0
      ffluidsum2_fbc=0.0
      ffluidsum3_fbc=0.0
      if(FabricMesh_Presence) then
           FabricMesh_force=0.0
      endif
      if(FabricLine_Presence) then
           FabricLine_force=0.0
      endif
      endif

      endif  !ibm_solve

      n=1

c imposing bcs
      call boundary_ibm

      if(nr_ibm_fib.ge. 1) then
         f_impuls1_fib=0.d0
         f_impuls2_fib=0.d0
         f_impuls3_fib=0.d0
         force_points_fib=0.d0
         f_link1_fib=0.0
         f_link2_fib=0.0
         f_link3_fib=0.0
       endif

      if(nr_ibm_esh.ge. 1) then
         f_impuls1_esh=0.d0
         f_impuls2_esh=0.d0
         f_impuls3_esh=0.d0
         areapnt_esh=0.0
!         xpnt_esh=0.0
         force_points_esh=0.d0
         f_link1_esh=0.0
         f_link2_esh=0.0
         f_link3_esh=0.0
      endif

      if(nr_ibm_fsh.ge. 1) then
         f_impuls1_fsh=0.d0
         f_impuls2_fsh=0.d0
         f_impuls3_fsh=0.d0
         force_points_fsh=0.d0
         f_link1_fsh=0.0
         f_link2_fsh=0.0
         f_link3_fsh=0.0

         f_impulsHis1_fsh=0.d0
         f_impulsHis2_fsh=0.d0
         f_impulsHis3_fsh=0.d0
       endif

      if(nr_ibm_fbc.ge. 1) then
         f_impuls1_fbc=0.d0
         f_impuls2_fbc=0.d0
         f_impuls3_fbc=0.d0
         areapnt_fbc=0.0

         force_points_fbc=0.d0
         f_link1_fbc=0.0
         f_link2_fbc=0.0
         f_link3_fbc=0.0
      endif

!=============================  MAIN PART of the SOLVER ========================
      do ibdyg=1,nr_ibm
#ifdef INCLUDE_FIB
        if(bodytype(ibdyg) .eq. fiber) then
          include "r_stang_Fib.Fi"
        endif
#endif         
#ifdef INCLUDE_ESH
        if (bodytype(ibdyg) .eq. shell) then
          include "r_stang_Esh.Fi"
        endif
#endif         
#ifdef INCLUDE_FSH
        if (bodytype(ibdyg) .eq. fibrousshell) then
          include "r_stang_Fsh.Fi"
        endif
#endif
#ifdef INCLUDE_FBC
        if (bodytype(ibdyg) .eq. fabric) then
          include "r_stang_Fbc.Fi"
        endif
#endif
     
!===============================================================================

!CAUTION:need to be rewrite for fibrous shell
      if(iksh1 .ne. 0) then
        imaster=masterbdy(ibdyg)
        if(BodyType(imaster) .ne. shell) then 
           write(*,*)'The master should be shell'
          stop
        endif
        if(BodyType(ibdyg) .eq. fiber) then 
        do noi=1,ns_ibm_r_fib(ibdy)
           ntemi=npos_fib(ibdy,noi)
           do noj=1,ns_ibm_r_fib(ibdy)
              ntemj=npos_fib(ibdy,noj) 
              aamat_total(imaster,ntemi,ntemj)=
     &        aamat_total(imaster,ntemi,ntemj)
     &                   +aamat_fib(noi,noj)
           end do
           bbmat_total(imaster,ntemi,1)=bbmat_total(imaster,ntemi,1)
     &                                +fibm1_fib(i,noi)
           bbmat_total(imaster,ntemi,2)=bbmat_total(imaster,ntemi,2)
     &                                +fibm2_fib(i,noi)
           if(ndim .gt. 2) 
     &     bbmat_total(imaster,ntemi,3)=bbmat_total(imaster,ntemi,3)
     &                                +fibm3_fib(i,noi)
        end do
        endif

        if(BodyType(ibdyg) .eq. shell) then 
        do noi=1,ns_ibm_r_esh(ibdy)
           ntemi=npos_esh(ibdy,noi)
           do noj=1,ns_ibm_r_esh(ibdy)
              ntemj=npos_esh(ibdy,noj) 
              aamat_total(imaster,ntemi,ntemj)=
     &        aamat_total(imaster,ntemi,ntemj)
     &                   +aamat_esh(noi,noj)
           end do
           bbmat_total(imaster,ntemi,1)=bbmat_total(imaster,ntemi,1)
     &                                +fibm1_esh(i,noi)
           bbmat_total(imaster,ntemi,2)=bbmat_total(imaster,ntemi,2)
     &                                +fibm2_esh(i,noi)
           if(ndim .gt. 2) 
     &     bbmat_total(imaster,ntemi,3)=bbmat_total(imaster,ntemi,3)
     &                                +fibm3_esh(i,noi)
        end do
        endif
      end if
      

      end do  !ibdy


!CAUTION:need to be rewrite for fibrous shell

      if(iksh1 .ne. 0) then
         jq=1
         do imaster=1,Nr_IBM_mas
         i=imaster
        if(BodyType(i) .ne. shell) then 
           write(*,*)'The master should be shell'
          stop
        endif
         if(ntermiMAXm(imaster) .le. 1)  then      
              ntemi=0
              do noj=1,ns_ibm_rall(i)
                 iacsr(noj)=ntemi+1
                 do noi=1,ns_ibm_rall(i)
                    ksh_tmp=aamat_total(imaster,noj,noi)
                    if(abs(ksh_tmp) .gt. 1d-12) then
                       ntemi=ntemi+1
                    endif
                 enddo
              enddo
              ntermiMAXm(imaster)=ntemi+1 
!         k_flagMAT=.false.
         allocate(aacsrm(ntermiMAXm(imaster))
     &            ,jacsrm(ntermiMAXm(imaster)))
         aacsrm=0.0
         jacsrm=0
         iacsrm=0
         ntemi=0
         do noj=1,ns_ibm_rall(i)
           iacsrm(noj)=ntemi+1
           do noi=1,ns_ibm_rall(i)
             ksh_tmp=aamat_total(imaster,noj,noi)
             if(abs(ksh_tmp) .gt. 1d-12) then
               ntemi=ntemi+1
               aacsrm(ntemi)=ksh_tmp
               jacsrm(ntemi)=noi
             endif
           enddo
          enddo
          iacsrm(noj)=ntemi+1 
          ntermiMAXmaster=ntermiMAXm(imaster)
         endif

         i2l=Iglbloc_esh(imaster)
          call gmres_csr(
     &      aacsrm(1: ntermiMAXmaster)
     &     ,jacsrm(1: ntermiMAXmaster)
     &     ,iacsrm(1: ns_ibm_rall(i)+1)
     &     ,bbmat_total(imaster,1: ns_ibm_rall(i),1)
     &     ,gx_ibm_massive_esh(i2l,1: ns_ibm_rall(i))
     &     ,tol,ns_ibm_rall(i),ntermiMAXmaster)

          call gmres_csr(
     &      aacsrm(1: ntermiMAXmaster)
     &     ,jacsrm(1: ntermiMAXmaster)
     &     ,iacsrm(1: ns_ibm_rall(i)+1)
     &     ,bbmat_total(imaster,1: ns_ibm_rall(i),2)
     &     ,gy_ibm_massive_esh(i2l,1: ns_ibm_rall(i))
     &     ,tol,ns_ibm_rall(i),ntermiMAXmaster)
         if(ndim .gt. 2) 
     &      call gmres_csr(
     &      aacsrm(1: ntermiMAXmaster)
     &     ,jacsrm(1: ntermiMAXmaster)
     &     ,iacsrm(1: ns_ibm_rall(i)+1)
     &     ,bbmat_total(imaster,1: ns_ibm_rall(i),3)
     &     ,gz_ibm_massive_esh(i2l,1:ns_ibm_rall(i))
     &     ,tol,ns_ibm_rall(i),ntermiMAXmaster)
c initialized master matrices for later use
         aamat_total(imaster,1: ns_ibm_rall(i),1: ns_ibm_rall(i))=0.0
         bbmat_total(imaster,1: ns_ibm_rall(i),1:3)=0.0
      enddo
      do ibdy=1,nr_ibm
        i=ibdy
        imaster=masterbdy(ibdy)
        if(BodyType(i) .eq. fiber) then 
        i2l=Iglbloc_fib(i)
        do noi=1,ns_ibm_r_fib(i2l)
           ntemi=npos_fib(i2l,noi)
           gx_ibm_massive_fib(i2l,noi)=gx_ibm_massive_esh(imaster,ntemi)
           gy_ibm_massive_fib(i2l,noi)=gy_ibm_massive_esh(imaster,ntemi)
           if(ndim .gt. 2) 
     &    gz_ibm_massive_fib(i2l,noi)=gz_ibm_massive_esh(imaster,ntemi)
        end do
        elseif(BodyType(i) .eq. shell) then 
        i2l=Iglbloc_esh(i)
        do noi=1,ns_ibm_r_esh(i2l)
           ntemi=npos_esh(i2l,noi)
           gx_ibm_massive_esh(i2l,noi)=gx_ibm_massive_esh(imaster,ntemi)
           gy_ibm_massive_esh(i2l,noi)=gy_ibm_massive_esh(imaster,ntemi)
           if(ndim .gt. 2) 
     &     gz_ibm_massive_esh(i2l,noi)=gz_ibm_massive_esh(imaster,ntemi)
        end do
        endif
      end do ! ibdy
      end if !iksh1

      do ibdy=1,nr_ibm
        i=ibdy
        if(BodyType(i) .eq. fiber) then 
        i2l=Iglbloc_fib(i)
        do j=1,ns_ibm_r_fib(i2l)
          vibm1_fib(i2l,j)= 
     &        (gx_ibm_massive_fib(i2l,j)-gx_ibm_fib(i2l,j))/dt
          vibm2_fib(i2l,j)= 
     &        (gy_ibm_massive_fib(i2l,j)-gy_ibm_fib(i2l,j))/dt
          if(ndim .gt. 2) 
     &     vibm3_fib(i2l,j)= 
     &        (gz_ibm_massive_fib(i2l,j)-gz_ibm_fib(i2l,j))/dt
        end do

        elseif(BodyType(i) .eq. shell) then
        i2l=Iglbloc_esh(i)
        do j=1,ns_ibm_r_esh(i2l)
          vibm1_esh(i2l,j)= 
     &        (gx_ibm_massive_esh(i2l,j)-gx_ibm_esh(i2l,j))/dt
          vibm2_esh(i2l,j)= 
     &        (gy_ibm_massive_esh(i2l,j)-gy_ibm_esh(i2l,j))/dt
          if(ndim .gt. 2) 
     &     vibm3_esh(i2l,j)= 
     &        (gz_ibm_massive_esh(i2l,j)-gz_ibm_esh(i2l,j))/dt
        end do
        elseif(BodyType(i) .eq. fabric) then
        i2l=Iglbloc_fbc(i)
        do j=1,ns_ibm_r_fbc(i2l)
          vibm1_fbc(i2l,j)= 
     &        (gx_ibm_massive_fbc(i2l,j)-gx_ibm_fbc(i2l,j))/dt
          vibm2_fbc(i2l,j)= 
     &        (gy_ibm_massive_fbc(i2l,j)-gy_ibm_fbc(i2l,j))/dt
          if(ndim .gt. 2) 
     &     vibm3_fbc(i2l,j)= 
     &        (gz_ibm_massive_fbc(i2l,j)-gz_ibm_fbc(i2l,j))/dt
        end do
        if(FabricMeshFlag_fbc(i2l)) then
          do j=1,FabricMesh_numrecord(i2l)
             FabricMesh_v(i2l,j,1:3)= 
     &       ( FabricMesh_coordMass(i2l,j,1:3) 
     &       - FabricMesh_coord(i2l,j,1:3) )/dt
          enddo
        endif
        if(FabricLineFlag_fbc(i2l)) then
          do j=1,FabricLine_numrecord(i2l)
             FabricLine_v(i2l,j,1:3)= 
     &       ( FabricLine_coordMass(i2l,j,1:3) 
     &       - FabricLine_coord(i2l,j,1:3) )/dt
          enddo
        endif
        elseif(BodyType(i) .eq. fibrousshell) then
        i2l=Iglbloc_fsh(i)
        do jq=1,nq_ibm_r_fsh(i2l)
        do j=1,ns_ibm_r_fsh(i2l)
          vibm1_fsh(i2l,jq,j)= 
     &     (gx_ibm_massive_fsh(i2l,jq,j)-gx_ibm_fsh(i2l,jq,j))/dt
          vibm2_fsh(i2l,jq,j)= 
     &     (gy_ibm_massive_fsh(i2l,jq,j)-gy_ibm_fsh(i2l,jq,j))/dt
          if(ndim .gt. 2) 
     &     vibm3_fsh(i2l,jq,j)= 
     &     (gz_ibm_massive_fsh(i2l,jq,j)-gz_ibm_fsh(i2l,jq,j))/dt

        end do
        end do 
        endif 

      end do ! ibdy

      do i=1,nr_ibm_esh
       i2g=Ilocglb_esh(i)
      if((Genalpha_timesolver(i2g) .gt. 1)) then
       fibm1_esh(i,1:Ns_IBM_esh)=0.0
       fibm2_esh(i,1:Ns_IBM_esh)=0.0
       fibm3_esh(i,1:Ns_IBM_esh)=0.0

       if(FluidForceFlag(i2g)) then
         do j=1,ns_ibm_r_esh(i)
           fk_mass1_esh(i,j)=fk_mass1o_esh(i,j)
           fk_mass2_esh(i,j)=fk_mass2o_esh(i,j)
           ffluidsum1_esh(i,j)= fk_mass1_esh(i,j)
           ffluidsum2_esh(i,j)= fk_mass2_esh(i,j)
           if(ndim .gt. 2) then
            fk_mass3_esh(i,j)=fk_mass3o_esh(i,j)
            ffluidsum3_esh(i,j)= fk_mass3_esh(i,j)
           endif
         end do

         do j=1,ns_ibm_r_esh(i)
              fibm1_esh(i,j)=fibm1_esh(i,j)+ffluidsum1_esh(i,j)
              fibm2_esh(i,j)=fibm2_esh(i,j)+ffluidsum2_esh(i,j)
              if(ndim .gt. 2) 
     &           fibm3_esh(i,j)=fibm3_esh(i,j)+ffluidsum3_esh(i,j)
          end do
        endif  
      endif  !Genalpha_timesolver(i2g)
      end do  !nr_ibm_esh
      do i=1,nr_ibm_fbc
       i2g=Ilocglb_fbc(i)
      if((Genalpha_timesolver(i2g) .gt. 1)) then
       fibm1_fbc(i,1:Ns_IBM_fbc)=0.0
       fibm2_fbc(i,1:Ns_IBM_fbc)=0.0
       fibm3_fbc(i,1:Ns_IBM_fbc)=0.0

       if(FluidForceFlag(i2g)) then
        do j=1,ns_ibm_r_fbc(i)
           fk_mass1_fbc(i,j)=fk_mass1o_fbc(i,j)
           fk_mass2_fbc(i,j)=fk_mass2o_fbc(i,j)

            ffluidsum1_fbc(i,j)= fk_mass1_fbc(i,j)
            ffluidsum2_fbc(i,j)= fk_mass2_fbc(i,j)
            if(ndim .gt. 2) then
             fk_mass3_fbc(i,j)=fk_mass3o_fbc(i,j)
             ffluidsum3_fbc(i,j)= fk_mass3_fbc(i,j)
            endif
         end do

         do j=1,ns_ibm_r_fbc(i)
              fibm1_fbc(i,j)=fibm1_fbc(i,j)+ffluidsum1_fbc(i,j)
              fibm2_fbc(i,j)=fibm2_fbc(i,j)+ffluidsum2_fbc(i,j)
              if(ndim .gt. 2) 
     &           fibm3_fbc(i,j)=fibm3_fbc(i,j)+ffluidsum3_fbc(i,j)
          end do
        endif  
      endif  !Genalpha_timesolver(i2g)
      end do  !nr_ibm_fbc

      do i=1,nr_ibm_fsh
       i2g=Ilocglb_fsh(i)

       if(FluidForceFlag(i2g)) then
       do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
           fk_mass1_fsh(i,jq,j)=fk_mass1o_fsh(i,jq,j)
           fk_mass2_fsh(i,jq,j)=fk_mass2o_fsh(i,jq,j)
           ffluidsum1_fsh(i,jq,j)= fk_mass1_fsh(i,jq,j)
           ffluidsum2_fsh(i,jq,j)= fk_mass2_fsh(i,jq,j)
           if(ndim .gt. 2) then
            fk_mass3_fsh(i,jq,j)=fk_mass3o_fsh(i,jq,j)
            ffluidsum3_fsh(i,jq,j)= fk_mass3_fsh(i,jq,j)
           endif
        end do
       end do

       do jq=1,nq_ibm_r_fsh(i)
         do j=1,ns_ibm_r_fsh(i)

              fibm1_fsh(i,jq,j)=fibm1_fsh(i,jq,j)+ffluidsum1_fsh(i,jq,j)
              fibm2_fsh(i,jq,j)=fibm2_fsh(i,jq,j)+ffluidsum2_fsh(i,jq,j)
              if(ndim .gt. 2) 
     &        fibm3_fsh(i,jq,j)=fibm3_fsh(i,jq,j)+ffluidsum3_fsh(i,jq,j)

          end do
        end do
        endif
      end do  !nr_ibm_fsh


3000    format(1x,24(f8.4,2x))
3001    format(1x,1(i5,2x),24(f8.4,2x))
3002    format(1x,2(i5,2x),24(f8.4,2x))
3003    format(1x,3(i5,2x),24(f8.4,2x))
3004    format(1x,4(i5,2x),24(f8.4,2x))
3005    format(1x,5(i5,2x),24(f8.4,2x))
3006    format(1x,6(i5,2x),24(f8.4,2x))
3007    format(1x,7(i5,2x),24(f8.4,2x))
3008    format(1x,8(i5,2x),24(f8.4,2x))
3009    format(1x,9(i5,2x),24(f8.4,2x))
3010    format(1x,10(i5,2x),24(f8.4,2x))
3011    format(1x,11(i5,2x))

      return 
      end

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine r_stang_fib(ibm_solve,monitorON,isavesigma)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c solve for 2d and 3d structures
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      implicit none

      integer ibdy,np,npg,bodypoint(3),ne,jtype,itype,nos,ntem,noi,noj
     &     ,ntemi, ntemj

      real*8 wei_i,iner1,iner2,jner1,jner2
      integer ntype,info,isavesigma
      logical monitorON 
      integer ibm_solve
      real*8 mms(ns_ibm_esh,ns_ibm_esh)

      real*8 x(maxordern,3), sx(2,maxordern), mass(maxordern,maxordern)
     &    ,ks(maxordern,maxordern), kb(maxordern,maxordern)
     &    ,xpos(3),dxpos(2,3),d2xpos(3,3)
     &    ,ms(maxordern,maxordern)

       real*8 s(1,maxordern), dsdv(2,maxordern),ddsdvdv(3,maxordern)
     &    , dmat(2,2),kmat(2,2)
     &    , kmatbar(2,2) , smat(2,2), bmat(3,3), jmat(2,2)
     &    , gmat(3,2), hmat(3,3)
     &    , amat1(2,maxordern), amat2(3,maxordern)
     &    , amat1T(maxordern,2), amat2T(maxordern,3)

      integer i,j,k,k2,para,t_temp,n,j2,kcontact,jq,jq2
      integer ierrksh1,iksh1,iksh2,iksh3,i2,iksh1b,iksh1c

      real*8 aamat_fib(ns_ibm_fib,ns_ibm_fib)   !tobe changes
      real*8 aamat_total(nr_ibm_mas,ns_ibm_esh,ns_ibm_esh)    
     &    ,bbmat_total(nr_ibm_mas,ns_ibm_esh,3)   !tobe changes

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsr(:),iacsr(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsr(:)

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsrm(:),iacsrm(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsrm(:)

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsrfibershell(:),iacsrfibershell(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsrfibershell(:),vecRHSfibershell(:),
     &             vecXfibershell(:)

! contact parameters
      REAL*8 gau_i(3,ngaumax)
     &      ,xcont(maxordern,3),dxcont(3)
     &      ,ddcont,dncont,Fcont
     &      ,dXP1cont(2,3),Normalcont1(3),outcont(3),outcont2(3)
     &      ,Normalcont2(3),ContactForce
      REAL*8,SAVE, ALLOCATABLE :: 
     &             l_loc(:,:),iv_loc(:,:,:)
     &             ,phi_loc(:,:,:,:)
      INTEGER ntypeMAx_loc
     &        ,ItypecontELE1(3),ItypecontELE(3)
     &        ,ne2,ntype2,itype2,nos2,jtype2,ntem2

      INTEGER, SAVE :: ntermiMAX,ntermiMAXfibershell,ntermiMAXmaster 
      INTEGER, SAVE :: ntermiMAXm(nr_ibm_mas) 
!      real aamcsr(nr_ibm_mas,ns_ibm*20)
!      integer jamcsr(nr_ibm_mas,ns_ibm*20),iamcsr(nr_ibm_mas,ns_ibm)

      real temp1, temp2, temp3, temp4,ksh_tmp
      real temp_ibm1,temp_ibm2,temp_ibm3
     &      ,temp_ibm1b,temp_ibm2b,temp_ibm3b
     &      ,fn1,fn2,fn3
     &      ,temp_check1,temp_check2,temp_check3, temp2_ibm, temp3_ibm
      real   temp4_1_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)    
     &      , temp4_2_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)    
     &      , temp4_3_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)     
     &      ,a_fib(ns_ibm_fib,-3:3) 
     &      ,mat(4,4)
     &      ,tdx0_1_fib(nr_ibm_fib)    
     &      ,tdx0_2_fib(nr_ibm_fib)    
     &      ,tdx0_3_fib(nr_ibm_fib)    
     &      ,tdxnm1_1_fib(nr_ibm_fib)    
     &      ,tdxnm1_2_fib(nr_ibm_fib)    
     &      ,tdxnm1_3_fib(nr_ibm_fib)    

      real*8  s_function
      logical fix_bc,free_bc,fix_bcb,free_bcb 
      integer itmax,precond,success,itconverge
      real*8 tol,h_channel,KB_ibmt
    
      real*8 dt,dt2,time_alphaf 

      real*8 tmpv1(3),tmpv2(3),tmpv3(3),tmpv1x(3)
      real*8 areapnt_esh(nr_ibm_esh,ns_ibm_esh)
!     &       ,xpnt_esh(nr_ibm_esh,ns_ibm_esh)  !for inflation problem 

      integer jbndy,jqbndy

      real*8  X_i_j(3),X_ip_j(3),X_i_jp(3)
     &       ,X_imf_j(3),X_i_jpf(3),X_i_jmf(3)
     &       ,X_imf_jp(3)
     &       ,X_ip_jmf(3),X_imf_jpf(3)
     &       ,X_imf_jmf(3)
     &       ,DXX_i_j(3),DYY_i_j(3),DXY_i_j(3),DYX_i_j(3)
     &       ,DX_im_jm(3),DX_im_j(3),DX_im_jp(3)
     &       ,DX_i_jm(3),DX_i_j(3),DX_i_jp(3)
     &       ,DX_ip_jm(3),DX_ip_j(3),DX_ip_jp(3)
     &       ,DY_im_jm(3),DY_im_j(3),DY_im_jp(3)
     &       ,DY_i_jm(3),DY_i_j(3),DY_i_jp(3)
     &       ,DY_ip_jm(3),DY_ip_j(3),DY_ip_jp(3)
     &       ,kmem_i_j(50,3),kbend_i_j

      integer ibkt,ibkt2,ibkt_i,ibkt2_i
     &       ,ibktmp1,ibktx,ibkty,ibktz,ibkts1,ibkts2,ibkts3
     &       ,ibktmin(3),ibktmax(3)

      real*8  dvcont(3),dnvcont,Vcont(maxordern,3)
     &       ,VP1cont(3),VP2cont(3)
      integer i2l,i2g, ibdyg,i3g

      real*8 target_timefn
      REAL*8,SAVE, ALLOCATABLE :: ForceOut(:)
!temporal variables
      real*8  real_distibm(50),rtmpusrfn(50)
      integer int_distibm(50),itmpusrfn(50)
      integer jq_pm(2),j_pm(2)
      real*8  contactplane_n(4),contactplane_point(3)
      real*8  matFib_Esh(3,3)
      real*8  R_Fibmat(3,3),xposFib(3),x_FibVec(maxordern,3)
     &       ,n_FibVec(3),n3_FibVec(3)

!variable for genralized alpha
      real*8  timeMarchcoef(20)
     &       ,d_smat(3,3)
     &       ,mattmp_eye(3,3)
     &       ,coef_iter

      integer iter

      real*8, ALLOCATABLE ::
!     (ns_ibm_fbc or (3*ns_ibm_fbc,
     &     matdir1(:,:)
!     (3,maxordern)
     &     ,matdir2(:,:)
!     (3,maxordern)
     &     ,matdir3(:,:)
!     (3,maxordern)
     &     ,Dmat_Stiff1(:,:)
!     (maxordern,3)
     &     ,Dmat_Stiff2(:,:)
!     (maxordern,3)
     &     ,Dmat_Stiff3(:,:)
!     (maxordern,3)
     &     ,kb11(:,:)
!     (maxordern,maxordern)
     &     ,kb12(:,:)
!     (maxordern,maxordern)
     &     ,kb13(:,:)
!     (maxordern,maxordern)
     &     ,kb21(:,:)
!     (maxordern,maxordern)
     &     ,kb22(:,:)
!     (maxordern,maxordern)
     &     ,kb23(:,:)
!     (maxordern,maxordern)
     &     ,kb31(:,:)
!     (maxordern,maxordern)
     &     ,kb32(:,:)
!     (maxordern,maxordern)
     &     ,kb33(:,:)
!     (maxordern,maxordern)
     &     ,kks(:,:)
!     (ns_ibm_esh or (3*ns_ibm_esh,
     &     ,kkb(:,:)
!     (ns_ibm_esh or (3*ns_ibm_esh,
     &     ,kksf(:,:)
!     (ns_ibm_esh or
     &     ,kkbf(:,:)

      logical ilucal_flag

      integer iv0,iv1,iv2,iv3,ie1,ie2

      real*8 tmpcoord(3),tmpcoord2(3),tmpcoord3(3)
     &      ,cf_coef,cf_coef2,cf_coef1
     &      ,x0,x1,x2,x3
     &      ,Ktmp(4,3,4,3)
     &      ,L_over_Lmax
     &      ,sinsign
     &      ,costet0,sintet0

      INTEGER iDX0,iKS0,N_iPARA0,iPARA0(100),i_FABRICLINE0

      REAL*8,  ALLOCATABLE :: 
     &             vecRHS(:,:)
     &            ,x0LHS(:,:)

      real*8 maxtmpx
      Logical KtmpFlag
!=====newallocatables for FE shell Continuum model
      integer  ih
      real*8 
     &    weih_i,theta3
     &   ,ametric(3),bmetric(3),dnxpos(3)
     &   ,gmetric_cov(2,2),gmetric_con0(2,2),gmetric_con(2,2)
     &   ,CauchyGreen_C(3,3),CauchyGreen_Cbar(3,3),C33_old
     &   ,detC_inplane,tracC_inplane
     &   ,detC,tracC
     &   ,DiffPhi_C(3),DiffPhi_C33,Diff2Phi_CC33(3),Diff2Phi_C(3,3)
     &   ,Diff2Phi_C33,Sstress33 
     &   ,Smetric(3),CMAT(3,3),CMAT33,dC33_old,CMATi3(3)
     &   ,nmatMetric(3),mmatMetric(3),DmatMetric(3,3,3)
     &   ,abar3inv
     &   ,avectmp0(3),avectmp1(3),avectmp2(3),avectmp3(3)
     &   ,avectmp1b(3),avectmp2b(3),avectmp3b(3)
     &   ,q_contrav_bases0(2,3),MatposFib_con(3,3)

      real*8, allocatable, dimension(:,:)::
     &    Forcetemp1
     &   ,temp_ibmtype2,temp_ibmtype3
     &   ,a3bar_Drs,a3_Dr
     &   ,Da_ab_Dr,a3tilde_Dr,Db_ab_Dr

      real*8, allocatable, dimension(:,:,:) ::
     &    Da_ab_Drs,a3tilde_Drs,a3_Drs,Depsilon_ab_r,Dkappa_ab_r

      real*8, allocatable, dimension(:,:,:,:)::
     &    Db_ab_Drs,Depsilon_ab_rs

      real*8, allocatable, dimension(:,:,:,:,:)::
     &    Dkappa_ab_rs

      real*8   alpha_coefK(5),Elect_RHS 



!==========================================================

      parameter (tol = 1.0d-11, precond=1, itmax=200)
       
#ifdef INCLUDE_FIB
      iksh1=simulationtype
      if (iksh1.eq.0) then
       ! do nothing
      else
       print *,"expecting iksh1.eq.0"
       stop
      endif
!     iksh1=0; only fiber or only shell   iksh1 ~=0  fiber + shell 

c time step
      dt=delta_t
      dt2=dt**2.0
      eps=1.0d-8

      if(ibm_solve .ne. 2) then
       ! do nothing
      else
       print *,"expecting ibm_solve.ne.2"
       stop
      endif

c
      do i=1,nr_ibm_fib

       do j=1,ns_ibm_r_fib(i)

        fk_mass1_fib(i,j)=fk_mass1o_fib(i,j)
        fk_mass2_fib(i,j)=fk_mass2o_fib(i,j)
        fk_mass3_fib(i,j)=0.0
        ffluidsum1_fib(i,j)= fk_mass1o_fib(i,j)
        ffluidsum2_fib(i,j)= fk_mass2o_fib(i,j)
        ffluidsum3_fib(i,j)=0.0 
        fibm1_fib(i,j)=fk_mass1o_fib(i,j)
        fibm2_fib(i,j)=fk_mass2o_fib(i,j)
        fibm3_fib(i,j)=0.0
        if(ndim .gt. 2) then
         print *,"3d fibers not supported, try FSH"
         stop
         fk_mass3_fib(i,j)=fk_mass3o_fib(i,j)
         ffluidsum3_fib(i,j)= fk_mass3o_fib(i,j)
         fibm3_fib(i,j)=fk_mass3o_fib(i,j)
        endif
! j           
       end do

       do j=1,ns_ibm_r_fib(i)
        vibm1_pre_fib(i,j)=vibm1_fib(i,j)
        vibm2_pre_fib(i,j)=vibm2_fib(i,j)
        if(ndim .gt. 2) then
         vibm3_pre_fib(i,j)=vibm3_fib(i,j)
        endif
       end do

      end do  !i=1,nr_ibm_fib

401   format(3(i8,' '),20(f20.9,' '))

      n=1

c imposing bcs
      call boundary_ibm

      if(nr_ibm_fib.ge. 1) then
       f_impuls1_fib=0.d0
       f_impuls2_fib=0.d0
       f_impuls3_fib=0.d0
       force_points_fib=0.d0
       f_link1_fib=0.0
       f_link2_fib=0.0
       f_link3_fib=0.0
      endif


!=============================  MAIN PART of the SOLVER ========================
      do ibdyg=1,nr_ibm
       if(bodytype(ibdyg) .eq. fiber) then
        include "r_stang_Fib.Fi"
       else
        print *,"expecting (bodytype(ibdyg) .eq. fiber)"
        stop
       endif
!===============================================================================

!CAUTION:need to be rewrite for fibrous shell
       if(iksh1 .ne. 0) then
        print *,"expecting iksh1=0"
        stop
       end if

      end do  
!ibdy

!CAUTION:need to be rewrite for fiber + shell
      if(iksh1 .ne. 0) then
       print *,"expecting iksh1=0"
       stop
      end if 
!iksh1

      do ibdy=1,nr_ibm
       i=ibdy
       if(BodyType(i) .eq. fiber) then 
        i2l=Iglbloc_fib(i)
        do j=1,ns_ibm_r_fib(i2l)
         vibm1_fib(i2l,j)= 
     &    (gx_ibm_massive_fib(i2l,j)-gx_ibm_fib(i2l,j))/dt
         vibm2_fib(i2l,j)= 
     &    (gy_ibm_massive_fib(i2l,j)-gy_ibm_fib(i2l,j))/dt
         if(ndim .gt. 2) then
          vibm3_fib(i2l,j)= 
     &     (gz_ibm_massive_fib(i2l,j)-gz_ibm_fib(i2l,j))/dt
         endif
        end do

       elseif(BodyType(i) .eq. shell) then
        print *,"this code should not be called"
        stop
       elseif(BodyType(i) .eq. fabric) then
        print *,"this code should not be called"
        stop
       elseif(BodyType(i) .eq. fibrousshell)then
        print *,"this code should not be called"
        stop
       else
        print *,"expecting BodyType(i)=fiber"
        stop
       endif 

      end do ! ibdy

      do i=1,nr_ibm_esh
       print *,"this code should not be called"
       stop
      end do  !nr_ibm_esh
      do i=1,nr_ibm_fbc
       print *,"this code should not be called"
       stop
      end do  !nr_ibm_fbc

      do i=1,nr_ibm_fsh
       print *,"this code should not be called"
       stop
      end do  !nr_ibm_fsh

3000  format(1x,24(f8.4,2x))
3001  format(1x,1(i5,2x),24(f8.4,2x))
3002  format(1x,2(i5,2x),24(f8.4,2x))
3003  format(1x,3(i5,2x),24(f8.4,2x))
3004  format(1x,4(i5,2x),24(f8.4,2x))
3005  format(1x,5(i5,2x),24(f8.4,2x))
3006  format(1x,6(i5,2x),24(f8.4,2x))
3007  format(1x,7(i5,2x),24(f8.4,2x))
3008  format(1x,8(i5,2x),24(f8.4,2x))
3009  format(1x,9(i5,2x),24(f8.4,2x))
3010  format(1x,10(i5,2x),24(f8.4,2x))
3011  format(1x,11(i5,2x))

#else 
      print *,"include_fib not defined"
      stop
#endif

      return 
      end


      subroutine r_stang_fsh(ibm_solve,monitorON,isavesigma)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c solve for 2d and 3d structures
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      implicit none
      integer ibdy,np,npg,bodypoint(3),ne,jtype,itype,nos,ntem,noi,noj
     &     ,ntemi, ntemj

      real*8 wei_i,iner1,iner2,jner1,jner2
      integer ntype,info,isavesigma
      logical monitorON 
      integer ibm_solve
      real*8 mms(ns_ibm_esh,ns_ibm_esh)

      real*8 x(maxordern,3), sx(2,maxordern), mass(maxordern,maxordern)
     &    ,ks(maxordern,maxordern), kb(maxordern,maxordern)
     &    ,xpos(3),dxpos(2,3),d2xpos(3,3)
     &    ,ms(maxordern,maxordern)
     &    ,xj1,xj2,xj3,yj1,yj2,yj3,zj1,zj2,zj3

       real*8 s(1,maxordern), dsdv(2,maxordern),ddsdvdv(3,maxordern)
     &    , dmat(2,2),kmat(2,2)
     &    , kmatbar(2,2) , xi1, xi2, smat(2,2), bmat(3,3), jmat(2,2)
     &    , gmat(3,2), hmat(3,3)
     &    , amat1(2,maxordern), amat2(3,maxordern)
     &    , amat1T(maxordern,2), amat2T(maxordern,3)
     &    ,FIBM1pres_fib(Ns_IBM_fib,2)    
     &    ,FIBM2pres_fib(Ns_IBM_fib,2)    
     &    ,FIBM3pres_fib(Ns_IBM_fib,2)
     &    ,FIBM1pres_esh(Ns_IBM_esh)    
     &    ,FIBM2pres_esh(Ns_IBM_esh)    
     &    ,FIBM3pres_esh(Ns_IBM_esh)
     &  
     &    ,FIBM1cont_esh(Ns_IBM_esh)    
     &    ,FIBM2cont_esh(Ns_IBM_esh)    
     &    ,FIBM3cont_esh(Ns_IBM_esh)   
     &  
     &    ,normalpres_esh(Ns_IBM_esh,3) 
     &    ,FIBM1pres_fsh(nq_ibm_fsh,ns_ibm_fsh)    
     &    ,FIBM2pres_fsh(nq_ibm_fsh,ns_ibm_fsh)    
     &    ,FIBM3pres_fsh(nq_ibm_fsh,ns_ibm_fsh)  

      real*8 
     &     FIBM1pres_fbc(Ns_IBM_fbc)    
     &    ,FIBM2pres_fbc(Ns_IBM_fbc)    
     &    ,FIBM3pres_fbc(Ns_IBM_fbc)
     &  
     &    ,FIBM1cont_fbc(Ns_IBM_fbc)    
     &    ,FIBM2cont_fbc(Ns_IBM_fbc)    
     &    ,FIBM3cont_fbc(Ns_IBM_fbc)   

      integer i,j,k,k2,para,t_temp,n,j2,kcontact,jq,jq2
      integer ierrksh1,iksh1,iksh2,iksh3,i2,iksh1b,iksh1c

      real*8 aamat_fib(ns_ibm_fib,ns_ibm_fib)   !tobe changes
      real*8 aamat_total(nr_ibm_mas,ns_ibm_esh,ns_ibm_esh)    
     &    ,bbmat_total(nr_ibm_mas,ns_ibm_esh,3)   !tobe changes

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsr(:),iacsr(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsr(:)

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsrm(:),iacsrm(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsrm(:)

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsrfibershell(:),iacsrfibershell(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsrfibershell(:),vecRHSfibershell(:),
     &             vecXfibershell(:)

! contact parameters
      REAL*8 gau_i(3,ngaumax)
     &      ,XP1cont(3),XP2cont(3),xcont(maxordern,3),dxcont(3)
     &      ,ddcont,dncont,Fcont
     &      ,dXP1cont(2,3),Normalcont1(3),outcont(3),outcont2(3)
     &      ,Normalcont2(3),ContactForce
      REAL*8,SAVE, ALLOCATABLE :: 
     &             l_loc(:,:),iv_loc(:,:,:)
     &             ,phi_loc(:,:,:,:)
      INTEGER ntypeMAx_loc
     &        ,ItypecontELE1(3),ItypecontELE(3)
     &        ,ne2,ntype2,itype2,nos2,jtype2,ntem2

      INTEGER, SAVE :: ntermiMAX,ntermiMAXfibershell,ntermiMAXmaster 
      INTEGER, SAVE :: ntermiMAXm(nr_ibm_mas) 
!      real aamcsr(nr_ibm_mas,ns_ibm*20)
!      integer jamcsr(nr_ibm_mas,ns_ibm*20),iamcsr(nr_ibm_mas,ns_ibm)

      real temp1, temp2, temp3, temp4,ksh_tmp
      real temp_ibm1,temp_ibm2,temp_ibm3
     &      ,temp_ibm1b,temp_ibm2b,temp_ibm3b
     &      ,fn1,fn2,fn3
     &      ,temp_check1,temp_check2,temp_check3, temp2_ibm, temp3_ibm
      real   temp4_1_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)    
     &      , temp4_2_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)    
     &      , temp4_3_ibm_fib(nr_ibm_fib,0:ns_ibm_fib)     
     &      ,a_fib(ns_ibm_fib,-3:3) 
     &      ,mat(4,4)
     &      ,tdx0_1_fib(nr_ibm_fib)    
     &      ,tdx0_2_fib(nr_ibm_fib)    
     &      ,tdx0_3_fib(nr_ibm_fib)    
     &      ,tdxnm1_1_fib(nr_ibm_fib)    
     &      ,tdxnm1_2_fib(nr_ibm_fib)    
     &      ,tdxnm1_3_fib(nr_ibm_fib)    

      real*8  s_function
      logical fix_bc,free_bc,fix_bcb,free_bcb 
      integer itmax,precond,success,itconverge
      real*8 tol,h_channel,KB_ibmt
    
      real*8 dt,dt2,time_alphaf 

      real*8 tmpv1(3),tmpv2(3),tmpv3(3),tmpv1x(3)
      real*8 areapnt_esh(nr_ibm_esh,ns_ibm_esh)
!     &       ,xpnt_esh(nr_ibm_esh,ns_ibm_esh)  !for inflation problem 

      integer jbndy,jqbndy

      real*8  X_i_j(3),X_ip_j(3),X_i_jp(3),X_ipf_j(3)
     &       ,X_imf_j(3),X_i_jpf(3),X_i_jmf(3)
     &       ,X_ipf_jp(3),X_imf_jp(3),X_ip_jpf(3)
     &       ,X_ip_jmf(3),X_ipf_jpf(3),X_imf_jpf(3)
     &       ,X_ipf_jmf(3),X_imf_jmf(3)
     &       ,DXX_i_j(3),DYY_i_j(3),DXY_i_j(3),DYX_i_j(3)
     &       ,DX_im_jm(3),DX_im_j(3),DX_im_jp(3)
     &       ,DX_i_jm(3),DX_i_j(3),DX_i_jp(3)
     &       ,DX_ip_jm(3),DX_ip_j(3),DX_ip_jp(3)
     &       ,DY_im_jm(3),DY_im_j(3),DY_im_jp(3)
     &       ,DY_i_jm(3),DY_i_j(3),DY_i_jp(3)
     &       ,DY_ip_jm(3),DY_ip_j(3),DY_ip_jp(3)
     &       ,kmem_i_j(50,3),kbend_i_j
     &       ,GXtmp(1:3,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &              ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh) 
     &       ,dGXtmp(1:3,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh
     &              ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh) 
     &        ,Siga(nq_ibm_fsh,ns_ibm_fsh,2,2) 
     &        ,SigB(4,nq_ibm_fsh+ns_ibm_fsh,3) 
     &       ,gama(nq_ibm_fsh,ns_ibm_fsh,2,2) 

      integer ibkt,ibkt2,ibkt_i,ibkt2_i
     &       ,ibktmp1,ibktx,ibkty,ibktz,ibkts1,ibkts2,ibkts3
     &       ,ibktmin(3),ibktmax(3)

      real*8  dvcont(3),dnvcont,Vcont(maxordern,3),xvv(maxordern,3)
     &       ,VP1cont(3),VP2cont(3)
      integer i2l,i2g, ibdyg,i3g

      real*8  Stiff_points_fsh,target_timefn
      REAL*8,SAVE, ALLOCATABLE :: Stiff_link_fsh(:,:)
      REAL*8,SAVE, ALLOCATABLE :: vnPiezo_fsh(:,:,:),aamat_fsh(:,:,:)
      REAL*8,SAVE, ALLOCATABLE :: VECpiezo_tmp(:,:,:,:)
      REAL*8,SAVE, ALLOCATABLE :: ForceOut(:)
!temporal variables
      real*8  real_distibm(50),rtmpusrfn(50)
      integer int_distibm(50),itmpusrfn(50)
      integer jq_pm(2),j_pm(2)
      real*8  contactplane_n(4),contactplane_point(3)
      real*8  matFib_Esh(3,3)
      real*8  R_Fibmat(3,3),xposFib(3),x_FibVec(maxordern,3)
     &       ,n_FibVec(3),n3_FibVec(3),temp_piezo

!variable for genralized alpha
      real*8 X_n_esh(ns_ibm_esh,3)
     &       ,dX_n_esh(ns_ibm_esh,3)
     &       ,ddX_n_esh(ns_ibm_esh,3)
     &       ,X_np1_esh(ns_ibm_esh,3)
     &       ,dX_np1_esh(ns_ibm_esh,3)
     &       ,ddX_np1_esh(ns_ibm_esh,3)
     &       ,X_na_esh(ns_ibm_esh,3)
     &       ,dX_na_esh(ns_ibm_esh,3)
     &       ,ddX_na_esh(ns_ibm_esh,3)
     &       ,Xtar_na_esh(ns_ibm_esh,3)
     &       ,dXtar_na_esh(ns_ibm_esh,3)
     &       ,X0_np1_esh(ns_ibm_esh,3)
     &       ,ddU_np1_esh(ns_ibm_esh,3)
     &       ,dX_in_esh(ns_ibm_esh,3)
     &       ,ddX_in_esh(ns_ibm_esh,3)
     &       ,timeMarchcoef(20)
     &       ,d_smat(3,3)
     &       ,mattmp_eye(3,3)
     &       ,coef_iter

      real*8 X_n_fbc(ns_ibm_fbc,3)
     &       ,dX_n_fbc(ns_ibm_fbc,3)
     &       ,ddX_n_fbc(ns_ibm_fbc,3)
     &       ,X_np1_fbc(ns_ibm_fbc,3)
     &       ,dX_np1_fbc(ns_ibm_fbc,3)
     &       ,ddX_np1_fbc(ns_ibm_fbc,3)
     &       ,X_na_fbc(ns_ibm_fbc,3)
     &       ,dX_na_fbc(ns_ibm_fbc,3)
     &       ,ddX_na_fbc(ns_ibm_fbc,3)
     &       ,Xtar_na_fbc(ns_ibm_fbc,3)
     &       ,dXtar_na_fbc(ns_ibm_fbc,3)
     &       ,X0_np1_fbc(ns_ibm_fbc,3)
     &       ,ddU_np1_fbc(ns_ibm_fbc,3)
     &       ,dX_in_fbc(ns_ibm_fbc,3)
     &       ,ddX_in_fbc(ns_ibm_fbc,3)

      real*8 X_n_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,dX_n_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddX_n_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,X_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,dX_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddX_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,X_na_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,dX_na_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddX_na_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,Xtar_na_fsh(nq_ibm_fsh*ns_ibm_fsh,3)
     &       ,dXtar_na_fsh(nq_ibm_fsh*ns_ibm_fsh,3)
     &       ,X0_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddU_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,dX_in_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,ddX_in_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,Siga_d(nq_ibm_fsh,ns_ibm_fsh,5)
     &       ,SigB_d(4,nq_ibm_fsh+ns_ibm_fsh,3)
     &       ,gama_d(nq_ibm_fsh,ns_ibm_fsh,5)



      real*8 XvPiezo_n_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,XvPiezo_np1_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,XvPiezo_na_fsh(nq_ibm_fsh,ns_ibm_fsh,3)
     &       ,XvPiezo_in_fsh(nq_ibm_fsh,ns_ibm_fsh,3)

      integer iter

      real*8, ALLOCATABLE ::
     &     aamat_esh(:,:)
!     (ns_ibm_esh or (3*ns_ibm_esh,
     &     ,aamat_fbc(:,:)
!     (ns_ibm_fbc or (3*ns_ibm_fbc,
     &     ,matdir1(:,:)
!     (3,maxordern)
     &     ,matdir2(:,:)
!     (3,maxordern)
     &     ,matdir3(:,:)
!     (3,maxordern)
     &     ,Dmat_Stiff1(:,:)
!     (maxordern,3)
     &     ,Dmat_Stiff2(:,:)
!     (maxordern,3)
     &     ,Dmat_Stiff3(:,:)
!     (maxordern,3)
     &     ,kb11(:,:)
!     (maxordern,maxordern)
     &     ,kb12(:,:)
!     (maxordern,maxordern)
     &     ,kb13(:,:)
!     (maxordern,maxordern)
     &     ,kb21(:,:)
!     (maxordern,maxordern)
     &     ,kb22(:,:)
!     (maxordern,maxordern)
     &     ,kb23(:,:)
!     (maxordern,maxordern)
     &     ,kb31(:,:)
!     (maxordern,maxordern)
     &     ,kb32(:,:)
!     (maxordern,maxordern)
     &     ,kb33(:,:)
!     (maxordern,maxordern)
     &     ,kks(:,:)
!     (ns_ibm_esh or (3*ns_ibm_esh,
     &     ,kkb(:,:)
!     (ns_ibm_esh or (3*ns_ibm_esh,
     &     ,kksf(:,:)
!     (ns_ibm_esh or
     &     ,kkbf(:,:)
!     (ns_ibm_esh or
     &     ,kkb_fbc(:,:)
!     (ns_ibm_fbc or (3*ns_ibm_fbc,
     &      ,mms_fbc(:,:) 

      real*8 smatsave_esh(3,ngaumax,Ns_IBMEle_esh)

      real*8 X22piezo(1:3),X2piezo(1:3),X11piezo(1:3),X1piezo(1:3)
     &      ,X12piezo(1:3),normalpiezo(1:3)
     &      ,Epiezo,Fpiezo,Gpiezo,Lpiezo,Mpiezo,Npiezo,Hnpiezo,Hdpiezo
     &      ,Ppiezo,Rpiezo
     &      ,VectempPiezo(3),Qmatpiezo(1:3,1:3),Pmatpiezo(1:3,1:3)

      logical ilucal_flag
      integer, save :: dist_esh(ns_ibm_esh)

      integer in_fbc
     &      ,iv0,iv1,iv2,iv3
     &      ,ie1,ie2

      real*8 EdgeLength_fbc(nr_ibm_fbc,numedgeMAX_fbc)
     &      ,EleCenter_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3)
     &      ,EleNormal_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3)
     &      ,EleEdge_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3)
     &      ,EleArea_fbc(nr_ibm_fbc,Ns_IBMEle_fbc)
     &      ,axtmp_fbc,aytmp_fbc,aztmp_fbc
     &      ,bxtmp_fbc,bytmp_fbc,bztmp_fbc
     &      ,cxtmp_fbc,cytmp_fbc,cztmp_fbc
     &      ,tmpcoord(3),tmpcoord2(3),tmpcoord3(3)
     &      ,total_volume_fbc(nr_ibm_fbc)
     &      ,total_Area_fbc(nr_ibm_fbc)
     &      ,cf_coef,cf_coef2,cf_coef1
     &      ,x0,x1,x2,x3
     &      ,y0,y1,y2,y3
     &      ,z0,z1,z2,z3
     &      ,Forcetmpfbc(4,3)
     &      ,Ktmp(4,3,4,3)
     &      ,L_over_Lmax
     &      ,sinsign
     &      ,costet0,sintet0

      INTEGER iDX0,iKS0,N_iPARA0,iPARA0(100),i_FABRICLINE0

      real*8 FabricML_dx(3)
     &      ,FabricML_att(FabricMeshL_numattMax)
     &      ,outcount2(3),outcount(3) 
     &      ,FabricML_tension,FabricML_Dtension
     &      ,FabricLine_curv(numLinepointRecordMax,3)

      real*8, dimension(numLinepointRecordMax,3) ::
     &       FabricLine_X_n,FabricLine_dX_n,FabricLine_ddX_n
     &      ,FabricLine_X_np1,FabricLine_dX_np1,FabricLine_ddX_np1
     &      ,FabricLine_ddX_in,FabricLine_dX_in
     &      ,FabricLine_X_na,FabricLine_dX_na,FabricLine_ddX_na
     &      ,FabricLine_force_points
     &      ,FabricLine_Xtar_np1,FabricLine_dXtar_np1,FabricLine_X0_np1
     &      ,FabricLine_Xtar_na,FabricLine_dXtar_na

      real*8, dimension(numMeshpointRecordMax,3) ::
     &       FabricMesh_X_n,FabricMesh_dX_n,FabricMesh_ddX_n
     &      ,FabricMesh_X_np1,FabricMesh_dX_np1,FabricMesh_ddX_np1
     &      ,FabricMesh_ddX_in,FabricMesh_dX_in
     &      ,FabricMesh_X_na,FabricMesh_dX_na,FabricMesh_ddX_na
     &      ,FabricMesh_force_points
     &      ,FabricMesh_Xtar_np1,FabricMesh_dXtar_np1,FabricMesh_X0_np1
     &      ,FabricMesh_Xtar_na,FabricMesh_dXtar_na

      REAL*8,  ALLOCATABLE :: 
     &             vecRHS(:,:)
     &            ,x0LHS(:,:)

      Logical tensionDerivFlag_fbc

      real*8 maxtmpx,normalpres_fbc(nr_ibm_fbc,3)
      Logical KtmpFlag
!=====newallocatables for FE shell Continuum model
      integer  ih
      real*8 
     &    weih_i,theta3
     &   ,ametric(3),bmetric(3),dnxpos(3)
     &   ,gmetric_cov(2,2),gmetric_con0(2,2),gmetric_con(2,2)
     &   ,CauchyGreen_C(3,3),CauchyGreen_Cbar(3,3),C33_old
     &   ,detC_inplane,tracC_inplane
     &   ,detC,tracC
     &   ,DiffPhi_C(3),DiffPhi_C33,Diff2Phi_CC33(3),Diff2Phi_C(3,3)
     &   ,Diff2Phi_C33,Sstress33 
     &   ,Smetric(3),CMAT(3,3),CMAT33,dC33_old,CMATi3(3)
     &   ,nmatMetric(3),mmatMetric(3),DmatMetric(3,3,3)
     &   ,abar3inv
     &   ,avectmp0(3),avectmp1(3),avectmp2(3),avectmp3(3)
     &   ,avectmp1b(3),avectmp2b(3),avectmp3b(3)
     &   ,xposFib_con(3),q_contrav_bases0(2,3),MatposFib_con(3,3)

      real*8, allocatable, dimension(:)::
     &    fibmRHS1_esh,fibmRHS2_esh,fibmRHS3_esh
     &   ,temp_ibmtype,a3bar_Dr

      real*8, allocatable, dimension(:,:)::
     &    Forcetemp1
     &   ,temp_ibmtype2,temp_ibmtype3
     &   ,a3bar_Drs,a3_Dr
     &   ,Da_ab_Dr,a3tilde_Dr,Db_ab_Dr

      real*8, allocatable, dimension(:,:,:) ::
     &    Da_ab_Drs,a3tilde_Drs,a3_Drs,Depsilon_ab_r,Dkappa_ab_r

      real*8, allocatable, dimension(:,:,:,:)::
     &    Db_ab_Drs,Depsilon_ab_rs

      real*8, allocatable, dimension(:,:,:,:,:)::
     &    Dkappa_ab_rs

      real*8   dxcont_piezo,outcont_piezo
     &   ,alpha_coefK(5),Elect_RHS 
      integer  n_dim_aamatPiezo_fsh,n_aamatPiezo_fsh
      real*8, allocatable, dimension(:) ::
     &    aamatPiezo_fsh

      integer  FlagForceCalculation(nq_ibm_fsh,ns_ibm_fsh)
!==========================================================

      parameter (tol = 1.0d-11, precond=1, itmax=200)
       
#ifdef INCLUDE_FSH
      iksh1=simulationtype
      if (iksh1.eq.0) then
       ! do nothing
      else
       print *,"expecting iksh1.eq.0"
       stop
      endif
!     iksh1=0; only fiber or only shell   iksh1 ~=0  fiber + shell 

c time step
      dt=delta_t
      dt2=dt**2.0
      eps=1.0d-8

      if(ibm_solve .ne. 2) then

       do i=1,nr_ibm_fsh
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)

         fk_mass1_fsh(i,jq,j)=fk_mass1o_fsh(i,jq,j)
         fk_mass2_fsh(i,jq,j)=fk_mass2o_fsh(i,jq,j)
         fk_mass3_fsh(i,jq,j)=0.0
         ffluidsum1_fsh(i,jq,j)= fk_mass1o_fsh(i,jq,j)
         ffluidsum2_fsh(i,jq,j)= fk_mass2o_fsh(i,jq,j)
         ffluidsum3_fsh(i,jq,j)=0.0 
         fibm1_fsh(i,jq,j)=fk_mass1o_fsh(i,jq,j)
         fibm2_fsh(i,jq,j)=fk_mass2o_fsh(i,jq,j)
         fibm3_fsh(i,jq,j)=0.0
         if(ndim .gt. 2) then
          fk_mass3_fsh(i,jq,j)=fk_mass3o_fsh(i,jq,j)
          ffluidsum3_fsh(i,jq,j)= fk_mass3o_fsh(i,jq,j)
          fibm3_fsh(i,jq,j)=fk_mass3o_fsh(i,jq,j)
         endif
        end do
        end do
!nr_ibm_fsh
       end do  
  
       do i=1,nr_ibm_fsh
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
         vibm1_pre_fsh(i,jq,j)=vibm1_fsh(i,jq,j)
         vibm2_pre_fsh(i,jq,j)=vibm2_fsh(i,jq,j)
         if(ndim .gt. 2) 
     &    vibm3_pre_fsh(i,jq,j)=vibm3_fsh(i,jq,j)
        end do
        enddo
       end do  !nr_ibm_fsh

401    format(3(i8,' '),20(f20.9,' '))

      else  ! ibm_solve ==2 and only structure

       if(nr_ibm_fsh.ge. 1) then
        fibm1_fsh=0.0
        fibm2_fsh=0.0
        fibm3_fsh=0.0
        ffluidsum1_fsh=0.0
        ffluidsum2_fsh=0.0
        ffluidsum3_fsh=0.0
       else
        print *,"expecting nr_ibm_fsh.ge.1"
        stop
       endif

      endif  !ibm_solve

      n=1

c imposing bcs
      call boundary_ibm

      if(nr_ibm_fsh.ge. 1) then
       f_impuls1_fsh=0.d0
       f_impuls2_fsh=0.d0
       f_impuls3_fsh=0.d0
       force_points_fsh=0.d0
       f_link1_fsh=0.0
       f_link2_fsh=0.0
       f_link3_fsh=0.0

       f_impulsHis1_fsh=0.d0
       f_impulsHis2_fsh=0.d0
       f_impulsHis3_fsh=0.d0
      else
       print *,"expecting nr_ibm_fsh.ge.1"
       stop
      endif

!=============================  MAIN PART of the SOLVER ========================
      do ibdyg=1,nr_ibm
       if (bodytype(ibdyg) .eq. fibrousshell) then
        include "r_stang_Fsh.Fi"
       else
        print *,"expecting (bodytype(ibdyg) .eq. fibrousshell)"
        stop
       endif
     
!===============================================================================

!CAUTION:need to be rewrite for fiber + shell
       if(iksh1 .ne. 0) then
        print *,"expecting iksh1=0"
        stop
       end if 
!iksh1
      end do  
!ibdy


!CAUTION:need to be rewrite for fiber + shell

      if(iksh1 .ne. 0) then
       print *,"expecting iksh1=0"
       stop
      end if 

      do ibdy=1,nr_ibm
       i=ibdy
       if(BodyType(i) .eq. fiber) then 
        print *,"this code should not be called"
        stop
       elseif(BodyType(i) .eq. shell) then
        print *,"this code should not be called"
        stop
       elseif(BodyType(i) .eq. fabric) then
        print *,"this code should not be called"
        stop
       elseif(BodyType(i) .eq. fibrousshell) then
        i2l=Iglbloc_fsh(i)
        do jq=1,nq_ibm_r_fsh(i2l)
        do j=1,ns_ibm_r_fsh(i2l)
         vibm1_fsh(i2l,jq,j)= 
     &    (gx_ibm_massive_fsh(i2l,jq,j)-gx_ibm_fsh(i2l,jq,j))/dt
         vibm2_fsh(i2l,jq,j)= 
     &    (gy_ibm_massive_fsh(i2l,jq,j)-gy_ibm_fsh(i2l,jq,j))/dt
         if(ndim .gt. 2) then
          vibm3_fsh(i2l,jq,j)= 
     &     (gz_ibm_massive_fsh(i2l,jq,j)-gz_ibm_fsh(i2l,jq,j))/dt
         endif
        end do
        end do 
       endif 

      end do ! ibdy

      do i=1,nr_ibm_esh
       print *,"this code should not be called"
       stop
      end do  !nr_ibm_esh
      do i=1,nr_ibm_fbc
       print *,"this code should not be called"
       stop
      end do  !nr_ibm_fbc

3000    format(1x,24(f8.4,2x))
3001    format(1x,1(i5,2x),24(f8.4,2x))
3002    format(1x,2(i5,2x),24(f8.4,2x))
3003    format(1x,3(i5,2x),24(f8.4,2x))
3004    format(1x,4(i5,2x),24(f8.4,2x))
3005    format(1x,5(i5,2x),24(f8.4,2x))
3006    format(1x,6(i5,2x),24(f8.4,2x))
3007    format(1x,7(i5,2x),24(f8.4,2x))
3008    format(1x,8(i5,2x),24(f8.4,2x))
3009    format(1x,9(i5,2x),24(f8.4,2x))
3010    format(1x,10(i5,2x),24(f8.4,2x))
3011    format(1x,11(i5,2x))

#else 
      print *,"include_fsh not defined"
      stop
#endif

      return 
      end





c----*|--.---------.---------.---------.---------.---------.---------.-|-------|



