c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!@c
      subroutine initialize_ibm(nIBM_rq,nIBM_rq_fsh
     &                         ,nIBM_r,nIBM_r_fib,nIBM_r_fsh
     &                         ,nIBM_r_esh,nIBM_r_fbc,dtypeDelta 
     &                         ,min_gridx
     &                         ,min_gridy
     &                         ,min_gridz
     &                         ,nx1,ny1,nz1
     &                         ,gx1,gy1,gz1
     &                         ,idimin
     &                         ,nRead
     &                         ,theboss)
!@q
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!@T 
!This subroutine sets the initial position of immersed boundaries and their mass
!@q
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      implicit none

      integer nIBM_r(nr_IBM)
      integer nIBM_r_fib(nr_IBM_fib)
      integer nIBM_r_fsh(nr_IBM_fsh)
      integer nIBM_r_esh(nr_IBM_esh)
      integer nIBM_r_fbc(nr_IBM_fbc)
     
      integer dtypeDelta(3) 

      integer nIBM_rq(nr_IBM),nIBM_rq_fsh(nr_IBM_fsh)

      integer i,j,k,ig,jq
! Eulerian grid spacing
      real*8 min_gridx,min_gridy,min_gridz

      integer nx1,ny1,nz1

      real*8 gx1(nx1),gy1(ny1),gz1(nz1)

      integer ibdy,np,bodypoint(3),npg,ne
     &    ,jtype,itype,nos,ntype,ntem,noi,noj,n 

       real*8 s(1,maxordern), dsdv(2,maxordern),ddsdvdv(3,maxordern)

      real*8 x(maxordern,3)
     &    ,xpos(3),dxpos(2,3),d2xpos(3,3), kmatbar_esh(2,2), jmat(2,2)
     &    , gmat(3,2), hmat(3,3),hmat0(3,3)
     &    , amat1(2,maxordern), amat2(3,maxordern),sx(2,maxordern)
     &    ,x_FibVec(maxordern,3),xFibvecLoc(3)

      real*8 :: ds_ibm_esh_hold

      real*8 das,xj1,xj2,xj3,yj1,yj2,yj3,zj1,zj2,zj3  

      integer idimin
      integer i2

      real*8 tmpv1val,tmpv2val,tmpv1(3),tmpv2(3),tmpv3(3),tmpv1x(3)
     &      ,tmpvec(3),tmpvecval
     &      ,tmpvecvalahat(2)
     &      ,matFib_Esh(3,3)
     &      ,kmem_i_j(50,3),kbend_i_j

      integer jq_pm(2),j_pm(2)
      integer itmpnode(3)

      integer nRead
      logical temprestart
      logical theboss

      real*8 
     &     GXtmp(1:3,1-Nq_IBMB_fsh:Nq_IBM_fsh+Nq_IBMB_fsh 
     &              ,1-Ns_IBMB_fsh:Ns_IBM_fsh+Ns_IBMB_fsh)


      real*8  coor1(3),coor2(3)
     &       ,X_i_j(3),X_ip_j(3),X_i_jp(3),X_ipf_j(3)
     &       ,X_imf_j(3),X_i_jpf(3),X_i_jmf(3)
     &       ,X_ipf_jp(3),X_imf_jp(3),X_ip_jpf(3)
     &       ,X_ip_jmf(3),X_ipf_jpf(3),X_imf_jpf(3)
     &       ,X_ipf_jmf(3),X_imf_jmf(3)
     &       ,DX_im_jp(3),DY_im_jp(3),DY_ip_jm(3)
     &       ,DX_ip_j(3),DY_ip_j(3),DX_i_jp(3),DY_i_jp(3)
     &       ,DX_im_j(3),DY_im_j(3),DX_i_jm(3),DY_i_jm(3)
     &       ,DXX_i_j(3),DYY_i_j(3),DXY_i_j(3),DYX_i_j(3)
     &       ,X22piezo(3),X2piezo(3),X11piezo(3)
     &       ,X1piezo(3),X12piezo(3),normalpiezo(3)
     &       ,Epiezo,Fpiezo,Gpiezo,Lpiezo,Mpiezo,Npiezo,Hnpiezo,Hdpiezo
     &       ,dist,ds_ibm0,dr_ibm0,dt_ibm0

      integer jbndy,jqbndy
       
! contact parameters
      integer ne2,ntype2,itype2,nos2,ItypecontELE(3),jtype2,ntem2

      real XP1cont(3),XP2cont(3),xcont2(3,3),xcont(3,3)
     &    ,dxcont(3),ddcont,contactlength

      INTEGER,SAVE, ALLOCATABLE :: 
     &             jacsrfibershell(:),iacsrfibershell(:)
      REAL*8,SAVE, ALLOCATABLE :: 
     &             aacsrfibershell(:),vecRHSfibershell(:),
     &             vecXfibershell(:)

      real*8 
     &        Siga(nq_ibm_fsh,ns_ibm_fsh,2,2) 
     &        ,SigB(4,nq_ibm_fsh+ns_ibm_fsh,3) 
     &       ,gama(nq_ibm_fsh,ns_ibm_fsh,2,2) 
     &       ,dXP1cont(2,3),Normalcont1(3)

      integer ntermiMAXfibershell,j2,i2g
      REAL*8,SAVE, ALLOCATABLE :: Stiff_link_fsh(:,:)
      real*8 temp_ibm1,temp_ibm3

      integer ibkt2
      real*8 contactplane_point(3),contactplane_n(4),dncont,ksh_tmp

      integer in_fbc
     &      ,iv0,iv1,iv2,iv3
     &      ,ie1,ie2

      real   
     &       axtmp_fbc,aytmp_fbc,aztmp_fbc
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

      Logical KtmpFlag
      integer nQPoints,FF(1),jq2
      real QPoints(1,3),DD(1),BB(1,3),PP(1,3),outcount2(3),temp_ibm2
!=================== new variables for continuum shell ===============
      real*8 dnxpos(3),ametric(3),bmetric(3)
     &   ,avectmp0(3),avectmp1(3),avectmp2(3)
     &   ,avectmp3(3),avectmp0b(2,3),avectmp0c(2,3),theta3
     &   ,gmetric_cov(2,2),gmetric_con(2,2)
      integer ih
      integer  FlagForceCalculation(nq_ibm_fsh,ns_ibm_fsh)

      call DistFSI_allocate_memory(0)

      if (nRead .eq. 1) then
       temprestart=.true.
       call write_restart_DISRIBM(temp_ibm1,1,1)
      else
        
      min_grid_x=min_gridx
      min_grid_y=min_gridy
      min_grid_z=min_gridz

      nxgrid_ibm=nx1
      nygrid_ibm=ny1
      nzgrid_ibm=nz1

      gx(1:nx1)=gx1(1:nx1)
      gy(1:ny1)=gy1(1:ny1) 
      gz(1:nz1)=gz1(1:nz1)

      ndim=idimin
      if(theboss) then
       write(*,*) 'ibm grid size: nr_ibm =',nr_ibm
       write (*,*) 'coordinate of ibm points in r'
       write(*,*) 'ibm grid size: ns_ibm_fib =',ns_ibm_fib
       write(*,*) 'ibm grid size: ns_ibm_fsh =',ns_ibm_fsh
       write(*,*) 'ibm grid size: nq_ibm_fsh =',nq_ibm_fsh
       write(*,*) 'ibm grid size: ns_ibm_esh =',ns_ibm_esh
       write(*,*) 'ibm grid size: ns_ibm_fbc =',ns_ibm_fbc
      endif
! initialized ibm grid or reading saved grid file
      ds_ibm0=0.0
      dr_ibm0=0.0
      dt_ibm0=0.0
      create_ibm_flag=.true.
      create_new_flow=.true.
      massive_ibm =.true.
      para_coor_flag= .false.

      if(theboss)  print*, 'start with IBM input read'
! read input files and create structures
       call create_grid2(ds_ibm0,dr_ibm0,contactlength,theboss)
      if(theboss)  print*, 'done with  IBM input read'

!write section information for ourtput
 
      if (create_ibm_flag) then
       if (create_new_flow) then
        do i=1,nr_ibm_fib
          do j=1,ns_ibm_fib
             fk_mass1_fib(i,j)=0.0d0
             fk_mass2_fib(i,j)=0.0d0
             fk_mass3_fib(i,j)=0.0d0
          enddo
          do j=1,ns_ibm_fib
             if(  j .le. ns_ibm_r_fib(i)
     &          ) then
!             write(14,*)gx_ibm_fib(i,j),gy_ibm_fib(i,j),gz_ibm_fib(i,j)
             end if
          end do
        end do 

        do i=1,nr_ibm_fsh
           do jq=1,nq_ibm_fsh
           do j=1,ns_ibm_fsh
             fk_mass1_fsh(i,jq,j)=0.0d0
             fk_mass2_fsh(i,jq,j)=0.0d0
             fk_mass3_fsh(i,jq,j)=0.0d0
            enddo
          enddo
          do jq=1,nq_ibm_fsh
          do j=1,ns_ibm_fsh
             if(   ( j .le. ns_ibm_r_fsh(i))
     &       .and. (jq .le. nq_ibm_r_fsh(i))) then
!             write(14,*)
!     &      gx_ibm_fsh(i,jq,j),gy_ibm_fsh(i,jq,j),gz_ibm_fsh(i,jq,j)
             end if
          end do
          end do
        end do 
        jq=0
        do i=1,nr_ibm_esh
          jq=max(jq,abs(ShellModelType_esh(i)))
          do j=1,ns_ibm_esh
             fk_mass1_esh(i,j)=0.0d0
             fk_mass2_esh(i,j)=0.0d0
             fk_mass3_esh(i,j)=0.0d0
          enddo
          do j=1,ns_ibm_esh
             if(  j .le. ns_ibm_r_esh(i)+ns_ibm_i_esh(i)
     &          ) then
!             write(14,*)gx_ibm_esh(i,j),gy_ibm_esh(i,j),gz_ibm_esh(i,j)
             end if
          end do
        end do 
!        close(14)

       if(jq .gt. 1) then
         if(.not. (allocated(gmetric_con0SAVE_esh)))  
     &      allocate(gmetric_con0SAVE_esh(Nr_IBM_esh 
     &                                   ,Ns_IBMEle_esh 
     &                                   ,ngaumax 
     &                                   ,ngaumaxh 
     &                                   ,2
     &                                   ,2) 
     &                                   )
         if(.not. (allocated(gbase_con0SAVE_esh)))  
     &      allocate(gbase_con0SAVE_esh(Nr_IBM_esh 
     &                                   ,Ns_IBMEle_esh 
     &                                   ,ngaumax 
     &                                   ,ngaumaxh 
     &                                   ,2
     &                                   ,3)
     &                                   )

         if(.not. (allocated(detC_inplane0SAVE_esh)))  
     &      allocate(detC_inplane0SAVE_esh(Nr_IBM_esh 
     &                                   ,Ns_IBMEle_esh 
     &                                   ,ngaumax 
     &                                   ,ngaumaxh)
     &                                   )
       endif


        do i=1,nr_ibm_fbc
          do j=1,ns_ibm_fbc
             fk_mass1_fbc(i,j)=0.0d0
             fk_mass2_fbc(i,j)=0.0d0
             fk_mass3_fbc(i,j)=0.0d0
          enddo
          do j=1,ns_ibm_fbc
             if(  j .le. ns_ibm_r_fbc(i)+ns_ibm_i_fbc(i)
     &          ) then
!             write(14,*)gx_ibm_esh(i,j),gy_ibm_esh(i,j),gz_ibm_esh(i,j)
             end if
          end do
        end do 
!        close(14)

        if (massive_ibm) then 
!        open (15,file='massibm.txt',form='formatted',status='unknown')
!         write(15,*) nr_ibm_fib,ns_ibm_fib
         do i=1,nr_ibm_fib
           do j=1,ns_ibm_r_fib(i)
!             write(15,*) mass_ibm_fib(i,j)
           end do
         end do
         do i=1,nr_ibm_fsh
           do jq=1,nq_ibm_r_fsh(i)
           do j=1,ns_ibm_r_fsh(i)
!             write(15,*) mass_ibm_fsh(i,jq,j)
           end do
           end do
         end do
         do i=1,nr_ibm_esh
           do j=1,ns_ibm_r_esh(i)+ns_ibm_i_esh(i)
!             write(15,*) mass_ibm_esh(i,j)
           end do
         end do
!         close(15)
         do i=1,nr_ibm_fbc
           do j=1,ns_ibm_r_fbc(i)+ns_ibm_i_fbc(i)
!             write(15,*) mass_ibm_esh(i,j)
           end do
         end do
!         close(15)
        end if

       end if !(create_new_flow)
      end if  !(create_ibm_flag)    
   

! initialize lagrangian coord. massive part
         if (create_new_flow) then
            do i=1,nr_ibm_fib
              do j=1,ns_ibm_r_fib(i)+ns_ibm_i_fib(i)
                gx_ibm_massive_fib(i,j)=gx_ibm_fib(i,j)
                gy_ibm_massive_fib(i,j)=gy_ibm_fib(i,j)
                gz_ibm_massive_fib(i,j)=gz_ibm_fib(i,j)
                gx_ibmpre_fib(i,j)=gx_ibm_fib(i,j)
                gy_ibmpre_fib(i,j)=gy_ibm_fib(i,j)
                gz_ibmpre_fib(i,j)=gz_ibm_fib(i,j)
              end do
            end do
            do i=1,nr_ibm_fsh
              do jq=1,nq_ibm_r_fsh(i)
              do j=1,ns_ibm_r_fsh(i)
                gx_ibm_massive_fsh(i,jq,j)=gx_ibm_fsh(i,jq,j)
                gy_ibm_massive_fsh(i,jq,j)=gy_ibm_fsh(i,jq,j)
                gz_ibm_massive_fsh(i,jq,j)=gz_ibm_fsh(i,jq,j)
                gx_ibmpre_fsh(i,jq,j)=gx_ibm_fsh(i,jq,j)
                gy_ibmpre_fsh(i,jq,j)=gy_ibm_fsh(i,jq,j)
                gz_ibmpre_fsh(i,jq,j)=gz_ibm_fsh(i,jq,j)
              end do
              end do
            end do
            do i=1,nr_ibm_esh
              do j=1,ns_ibm_r_esh(i)+ns_ibm_i_esh(i)
                gx_ibm_massive_esh(i,j)=gx_ibm_esh(i,j)
                gy_ibm_massive_esh(i,j)=gy_ibm_esh(i,j)
                gz_ibm_massive_esh(i,j)=gz_ibm_esh(i,j)
                gx_ibmpre_esh(i,j)=gx_ibm_esh(i,j)
                gy_ibmpre_esh(i,j)=gy_ibm_esh(i,j)
                gz_ibmpre_esh(i,j)=gz_ibm_esh(i,j)
              end do
            end do
            do i=1,nr_ibm_fbc
              do j=1,ns_ibm_r_fbc(i)+ns_ibm_i_fbc(i)
                gx_ibm_massive_fbc(i,j)=gx_ibm_fbc(i,j)
                gy_ibm_massive_fbc(i,j)=gy_ibm_fbc(i,j)
                gz_ibm_massive_fbc(i,j)=gz_ibm_fbc(i,j)
                gx_ibmpre_fbc(i,j)=gx_ibm_fbc(i,j)
                gy_ibmpre_fbc(i,j)=gy_ibm_fbc(i,j)
                gz_ibmpre_fbc(i,j)=gz_ibm_fbc(i,j)
              end do
             if(FabricMeshFlag_fbc(i)) then
                do j=1,FabricMesh_numrecord(i)
                  FabricMesh_coordMass(i,j,1:3)= 
     &              FabricMesh_coord(i,j,1:3)
                  FabricMesh_coordpre(i,j,1:3)= 
     &              FabricMesh_coord(i,j,1:3)  
                enddo
              endif
              if(FabricLineFlag_fbc(i)) then
                do j=1,FabricLine_numrecord(i)
                  FabricLine_coordMass(i,j,1:3)= 
     &              FabricLine_coord(i,j,1:3)
                  FabricLine_coordpre(i,j,1:3)= 
     &              FabricLine_coord(i,j,1:3)   
                enddo
                do j=1,FabricLine_numrecord(i)
                 i2=FabricLine_addressRev(i,j,1) 
                 jq=FabricLine_addressRev(i,j,2) 
                 if(jq .eq. FabricLine_npoint(i,i2)) then
                    outcount2(1:3)=0.0
                 elseif(jq .eq. 1) then
                    outcount2(1:3)=0.0
                 else
                   jq2=FabricLine_address(i,i2,jq+1)
                   k  =FabricLine_address(i,i2,jq-1)
                   tmpcoord(1:3) =FabricLine_coordMass(i,k,1:3) 
                   dxcont(1:3)   =FabricLine_coordMass(i,j,1:3)
                   tmpcoord2(1:3)=FabricLine_coordMass(i,jq2,1:3)

                   temp_ibm1=FabricLine_attr(i,k,iLl0)
                   temp_ibm2=FabricLine_attr(i,j,iLl0)
                   dncont=0.5*(temp_ibm1+temp_ibm1)

                   outcount2(1:3)
     &                   =( (dxcont(1:3)   -tmpcoord(1:3))/temp_ibm1
     &                     -(tmpcoord2(1:3)-dxcont(1:3)  )/temp_ibm2)
     &                    /(dncont+1d-9)
                 endif
                 FabricLine_curv0(i,j,1:3)
     &                   =outcount2(1:3)
                 enddo
              endif
            end do
         end if !(create_new_flow)

! find regular and staggered distance of fiber points 
!j'=j+1/2
!dsf_{j'}= | X_{j+1} - X_{j}    |
!ds_{j}  = | X_{j'}  - X_{j'-1} | 
      if (create_new_flow) then
        if (para_coor_flag .eqv. .false.) then
          do i=1,nr_ibm_fib
           s_start=1 
           s_end=ns_ibm_r_fib(i)+ns_ibm_i_fib(i)
           do j=1,s_end-1
             coor1(1)=gx_bp_fib(i,j)
             coor1(2)=gy_bp_fib(i,j)
             coor1(3)=gz_bp_fib(i,j)
             coor2(1)=gx_bp_fib(i,j+1)
             coor2(2)=gy_bp_fib(i,j+1)
             coor2(3)=gz_bp_fib(i,j+1)
             dsf_ibm_fib(i,j)=dist(coor1,coor2)
           end do
           dsf_ibm_fib(i,0)=dsf_ibm_fib(i,1)
           dsf_ibm_fib(i,s_end)=dsf_ibm_fib(i,s_end-1)
          end do
          do i=1,nr_ibm_fib
           do j=2,ns_ibm_r_fib(i)+ns_ibm_i_fib(i)
            ds_ibm_fib(i,j)=
     &            (dsf_ibm_fib(i,j-1) + dsf_ibm_fib(i,j))/ 2.0d0
           end do
            ds_ibm_fib(i,1)=ds_ibm_fib(i,2)
          end do

! do more calculation for the fibrous shell and find the 
! extended position of the shell
          do i=1,nr_ibm_fsh
           do jq=1,nq_ibm_r_fsh(i)
           do j=1,ns_ibm_r_fsh(i)
               GXtmp(1,jq,j)=gx_bp_fsh(i,jq,j)
               GXtmp(2,jq,j)=gy_bp_fsh(i,jq,j)
               GXtmp(3,jq,j)=gz_bp_fsh(i,jq,j)
           enddo
           enddo

           do jq=1,nq_ibm_r_fsh(i)
            if (BoundaryFibShell_fsh(i,1,1) .eq. 3) then
             do j=1-(Ns_IBMB_fsh+1),0,1
               jbndy=ns_ibm_r_fsh(i)+j
               GXtmp(1:3,jq,j)=GXtmp(1:3,jq,jbndy)
             enddo
            else
             do j=0,1-(1+Ns_IBMB_fsh),-1
               jbndy=j+1
               GXtmp(1:3,jq,j)=
     &            2.0*GXtmp(1:3,jq,jbndy)-GXtmp(1:3,jq,jbndy+1)
             enddo 
            endif
            if (BoundaryFibShell_fsh(i,1,2) .eq. 3) then
             do j=ns_ibm_r_fsh(i)+1,ns_ibm_r_fsh(i)+(1+Ns_IBMB_fsh),1
               jbndy=-ns_ibm_r_fsh(i)+j
               GXtmp(1:3,jq,j)=GXtmp(1:3,jq,jbndy)
             enddo
            else
             do j=ns_ibm_r_fsh(i)+1,ns_ibm_r_fsh(i)+(1+Ns_IBMB_fsh),1
               jbndy=j-1
               GXtmp(1:3,jq,j)=
     &            2.0*GXtmp(1:3,jq,jbndy)-GXtmp(1:3,jq,jbndy-1)
             enddo 
            endif
          enddo

           do j=1-(1+Ns_IBMB_fsh),ns_ibm_r_fsh(i)+(1+Ns_IBMB_fsh)
            if (BoundaryFibShell_fsh(i,2,1) .eq. 3) then
             do jq=1-(Nq_IBMB_fsh+1),0,1
               jqbndy=nq_ibm_r_fsh(i)+jq
               GXtmp(1:3,jq,j)=GXtmp(1:3,jqbndy,j)
             enddo
            else
             do jq=0,1-(Nq_IBMB_fsh+1),-1
               jqbndy=jq+1
               GXtmp(1:3,jq,j)=
     &            2.0*GXtmp(1:3,jqbndy,j)-GXtmp(1:3,jqbndy+1,j)
             enddo 
            endif
            if (BoundaryFibShell_fsh(i,2,2) .eq. 3) then
             do jq=nq_ibm_r_fsh(i)+1,nq_ibm_r_fsh(i)+(Nq_IBMB_fsh+1),1
               jqbndy=-nq_ibm_r_fsh(i)+jq
               GXtmp(1:3,jq,j)=GXtmp(1:3,jqbndy,j)
             enddo
            else
             do jq=nq_ibm_r_fsh(i)+1,nq_ibm_r_fsh(i)+(Nq_IBMB_fsh+1),1
               jqbndy=jq-1
               GXtmp(1:3,jq,j)=
     &            2.0*GXtmp(1:3,jqbndy,j)-GXtmp(1:3,jqbndy-1,j)
             enddo 
            endif
          enddo

! define non-uniform Lagrangian grid system:
!
! j'= j+ 1/2    i'=i+1/2
! dsf_ibm_{i',j}  = | X_{i+1,j} - X{i,j} |
! dsf2_ibm_{i,j'} = | X_{i,j+1} - X{i,j} |

! ds_ibm_{i,j}   = | X_{i',j} - X{i'-1,j} | 
! ds2_ibm_{i,j}  = | X_{i,j'} - X{i,j'-1} |

! dk_ibm_{i',j'}   = | X_{i+1,j'} - X{i,j'} | 
! dk2_ibm_{i',j'}  = | X_{i',j+1} - X{i',j} |

! dkf_ibm_{i,j'}   = | X_{i',j'} - X{i'-1,j'} | 
! dkf2_ibm_{i',j}  = | X_{i',j'} - X{i',j'-1} |
!
           do jq=0,nq_ibm_r_fsh(i)+1
           do j=0,ns_ibm_r_fsh(i)+1
                X_i_j (1:3)=GXtmp(1:3,jq,j)
                X_i_jp(1:3)=GXtmp(1:3,jq,j+1)
                X_ip_j(1:3)=GXtmp(1:3,jq+1,j)

                X_i_jpf(1:3)=0.5*(GXtmp(1:3,jq,j+1)+GXtmp(1:3,jq,j))
                X_i_jmf(1:3)=0.5*(GXtmp(1:3,jq,j-1)+GXtmp(1:3,jq,j))
                X_ipf_j(1:3)=0.5*(GXtmp(1:3,jq+1,j)+GXtmp(1:3,jq,j))
                X_imf_j(1:3)=0.5*(GXtmp(1:3,jq-1,j)+GXtmp(1:3,jq,j))

                X_ip_jpf(1:3)=
     &                    0.5*(GXtmp(1:3,jq+1,j+1)+GXtmp(1:3,jq+1,j))
                X_ip_jmf(1:3)=
     &                    0.5*(GXtmp(1:3,jq+1,j-1)+GXtmp(1:3,jq+1,j))
                X_ipf_jp(1:3)=
     &                    0.5*(GXtmp(1:3,jq+1,j+1)+GXtmp(1:3,jq,j+1))
                X_imf_jp(1:3)=
     &                    0.5*(GXtmp(1:3,jq-1,j+1)+GXtmp(1:3,jq,j+1))

                X_ipf_jpf(1:3)=
     &                    0.25*(GXtmp(1:3,jq+1,j+1)+GXtmp(1:3,jq+1,j)
     &                         +GXtmp(1:3,jq,j+1)  +GXtmp(1:3,jq,j))
                X_ipf_jmf(1:3)=
     &                    0.25*(GXtmp(1:3,jq+1,j-1)+GXtmp(1:3,jq+1,j)
     &                         +GXtmp(1:3,jq,j-1)  +GXtmp(1:3,jq,j))
                X_imf_jpf(1:3)=
     &                    0.25*(GXtmp(1:3,jq-1,j+1)+GXtmp(1:3,jq-1,j)
     &                         +GXtmp(1:3,jq,j+1)  +GXtmp(1:3,jq,j))
                X_imf_jmf(1:3)=
     &                    0.25*(GXtmp(1:3,jq-1,j-1)+GXtmp(1:3,jq-1,j)
     &                         +GXtmp(1:3,jq,j-1)  +GXtmp(1:3,jq,j))

                dsf_ibm_fsh(i,jq,j) = dist(X_ip_j,X_i_j)
                dsf2_ibm_fsh(i,jq,j)= dist(X_i_jp,X_i_j)

                ds_ibm_fsh(i,jq,j) = dist(X_ipf_j,X_imf_j)
                ds2_ibm_fsh(i,jq,j)= dist(X_i_jpf,X_i_jmf)

                dk_ibm_fsh(i,jq,j) = dist(X_ip_jpf,X_i_jpf)
                dk2_ibm_fsh(i,jq,j)= dist(X_ipf_jp,X_ipf_j)

                dkf_ibm_fsh(i,jq,j) = dist(X_ipf_jpf,X_imf_jpf)
                dkf2_ibm_fsh(i,jq,j)= dist(X_ipf_jpf,X_ipf_jmf)
           enddo
           enddo
          end do
        else   !(para_coor_flag .eq. .false.) then
        end if !(para_coor_flag .eq. .false.) then

! here force for feedback algorithm is calculated 
        do i=1,nr_ibm_fib
           do j=1,ns_ibm_r_fib(i)+ns_ibm_i_fib(i)-1
! computing (dx.dx)^{n}
              ds_ibm0=(gx_ibm_fib(i,j+1)-gx_ibm_fib(i,j))
     &                /dsf_ibm_fib(i,j)
              dr_ibm0=(gy_ibm_fib(i,j+1)-gy_ibm_fib(i,j))
     &                /dsf_ibm_fib(i,j) 
              dt_ibm0=(gz_ibm_fib(i,j+1)-gz_ibm_fib(i,j))
     &                /dsf_ibm_fib(i,j) 
! saving computed term for using later time
              fs_1_ibm_fib(i,j)=ds_ibm0**2.+dr_ibm0**2.+dt_ibm0**2.
           end do
        end do

! initializing velocities of structures
        if(theboss)  print*, 'initialization of vibm1'
        if (nr_ibm_fib>=1) then
           vibm1_fib=0.0d0;vibm1_pre_fib=0.0d0
           vibm2_fib=0.0d0;vibm2_pre_fib=0.0d0
           vibm3_fib=0.0d0;vibm3_pre_fib=0.0d0
        endif
        if (nr_ibm_fsh>=1) then
           vibm1_fsh=0.0d0;vibm1_pre_fsh=0.0d0;aibm1_fsh=0.0d0
           vibm2_fsh=0.0d0;vibm2_pre_fsh=0.0d0;aibm2_fsh=0.0d0
           vibm3_fsh=0.0d0;vibm3_pre_fsh=0.0d0;aibm3_fsh=0.0d0
        endif
        if (nr_ibm_esh>=1) then
           vibm1_esh=0.0d0;vibm1_pre_esh=0.0d0;aibm1_esh=0.0d0
           vibm2_esh=0.0d0;vibm2_pre_esh=0.0d0;aibm2_esh=0.0d0
           vibm3_esh=0.0d0;vibm3_pre_esh=0.0d0;aibm3_esh=0.0d0
        endif
        if (nr_ibm_fbc>=1) then
           vibm1_fbc=0.0d0;vibm1_pre_fbc=0.0d0;aibm1_fbc=0.0d0
           vibm2_fbc=0.0d0;vibm2_pre_fbc=0.0d0;aibm2_fbc=0.0d0
           vibm3_fbc=0.0d0;vibm3_pre_fbc=0.0d0;aibm3_fbc=0.0d0
        endif
      endif !(create_new_flow) then

! determine the width of delta function
      call delta_width(del_x,del_y,del_z
     &                ,delta_typex,delta_typey,delta_typez)
      close(12)
      if(nr_ibm_fsh>=1) ContactShellFlag_fsh=0
      if(nr_ibm_esh>=1) ContactShellFlag_esh=0
      if(nr_ibm_fbc>=1) ContactShellFlag_fbc=0

#ifdef INCLUDE_ESH
!!======== Perform intial calculations for FE shells
        do i=1,nr_ibm_esh
           ibdy=i
           i2g=Ilocglb_esh(ibdy)
! for free boundary condition
           jreal_2_imag_esh=0.d0
           do np=1,ns_ibm_r_esh(ibdy)  
              jreal_2_imag_esh(ibdy,np,np)=1.0
           end do       
           do np=1,ns_ibm_i_esh(ibdy)  
              npg=ns_ibm_r_esh(ibdy)+np
              bodypoint(1)=boundrylist_esh(ibdy,np,1)
              bodypoint(2)=boundrylist_esh(ibdy,np,2)
              bodypoint(3)=boundrylist_esh(ibdy,np,3)
              gx_ibm_massive_esh(ibdy,npg)=
     &                     gx_ibm_massive_esh(ibdy,bodypoint(1))
     &                    +gx_ibm_massive_esh(ibdy,bodypoint(2))       
     &                    -gx_ibm_massive_esh(ibdy,bodypoint(3))       

              gy_ibm_massive_esh(ibdy,npg)=
     &                     gy_ibm_massive_esh(ibdy,bodypoint(1))
     &                    +gy_ibm_massive_esh(ibdy,bodypoint(2))       
     &                    -gy_ibm_massive_esh(ibdy,bodypoint(3))       
    
              gz_ibm_massive_esh(ibdy,npg)=
     &                     gz_ibm_massive_esh(ibdy,bodypoint(1))
     &                    +gz_ibm_massive_esh(ibdy,bodypoint(2))       
     &                    -gz_ibm_massive_esh(ibdy,bodypoint(3))

              jreal_2_imag_esh(ibdy,npg,bodypoint(1))=1.0
              jreal_2_imag_esh(ibdy,npg,bodypoint(2))=1.0
              jreal_2_imag_esh(ibdy,npg,bodypoint(3))=-1.0


              gx_bp_esh0(ibdy,npg)=
     &                     gx_bp_esh0(ibdy,bodypoint(1))
     &                    +gx_bp_esh0(ibdy,bodypoint(2))       
     &                    -gx_bp_esh0(ibdy,bodypoint(3))       

              gy_bp_esh0(ibdy,npg)=
     &                     gy_bp_esh0(ibdy,bodypoint(1))
     &                    +gy_bp_esh0(ibdy,bodypoint(2))       
     &                    -gy_bp_esh0(ibdy,bodypoint(3))       
    
              gz_bp_esh0(ibdy,npg)=
     &                     gz_bp_esh0(ibdy,bodypoint(1))
     &                    +gz_bp_esh0(ibdy,bodypoint(2))       
     &                    -gz_bp_esh0(ibdy,bodypoint(3))
           end do       


           do ne=1,numelr(ibdy)
              jtype=neatype_esh(ibdy,ne)
              itype=jtypeall(jtype)
              if(itype .eq. 0) then
                ntype=12  
              else
                ntype=6+itype
              end if

              do nos=1,ntype
                     ntem=nea_esh(ibdy,ne,nos)
                     noj=1
                     x(nos,noj)=gx_ibm_massive_esh(ibdy,ntem)
                     noj=2
                     x(nos,noj)=gy_ibm_massive_esh(ibdy,ntem)
                     noj=3
                     x(nos,noj)=gz_ibm_massive_esh(ibdy,ntem)
              enddo
              if(itype .eq. 0) then
                     itmpnode(1)=4; itmpnode(2)=7; itmpnode(3)=8
              else
                     itmpnode(1)=1; itmpnode(2)=2; itmpnode(3)=itype+1
              endif
      
              tmpv1val=0.d0
              tmpv2val=0.d0
              do noj=1,3   
               tmpv1(noj)=x(itmpnode(2),noj)-x(itmpnode(1),noj)
               tmpv1val=tmpv1val+tmpv1(noj)**2
               tmpv2(noj)=x(itmpnode(3),noj)-x(itmpnode(1),noj)
               tmpv2val=tmpv2val+tmpv2(noj)**2
              enddo
              tmpv1val=sqrt(tmpv1val)
              tmpv2val=sqrt(tmpv2val)
              do noj=1,3   
               tmpv1(noj)=tmpv1(noj)/tmpv1val
               tmpv2(noj)=tmpv2(noj)/tmpv2val
              enddo
              call crossvec(tmpv3,tmpv1,tmpv2)

              call crossvec(tmpv2,tmpv3,tmpv1)

              crFibCrv_esh(ibdy,ne,1:3)=tmpv1(1:3)
              crFibCrv_esh(ibdy,ne,4:6)=tmpv3(1:3)

              do nos=1,ntype
                 do noj=1,3 
                    tmpvec(noj)=x(nos,noj)-x(itmpnode(1),noj)
                 enddo
                 tmpvecval=dot_product(tmpvec, tmpv3)
                 do noj=1,3 
                    tmpvec(noj)=tmpvec(noj)-tmpvecval*tmpv3(noj)
                 enddo
                 coorcurv_esh(ibdy,ne,nos,1)=dot_product(tmpvec, tmpv1)
                 coorcurv_esh(ibdy,ne,nos,2)=dot_product(tmpvec, tmpv2)

                 do noj=1,2
                   sx(noj,nos)=coorcurv_esh(ibdy,ne,nos,noj)
                 enddo
              enddo

!+++++++++++  SKIPED START+++++++++++++++++++++++++++++++++++++++
           goto 2121
!     for now skip over new implementation of Orthotropic shell
           if(Fibrous_Esh(ibdy) .eq. 1) then
            do nos=1,ntype
            ntem=nea_esh(ibdy,ne,nos)
            noj=1
            x_FibVec(nos,noj)=Materialpara_esh(ibdy,ntem,noj)
            noj=2
            x_FibVec(nos,noj)=Materialpara_esh(ibdy,ntem,noj)
            noj=3
            x_FibVec(nos,noj)=Materialpara_esh(ibdy,ntem,noj)
            enddo
            do ig=1,ngau0
             s(1,1:ntype)=nn0(jtype,ig,1:ntype)
             dsdv(1,1:ntype)=dnndv0(jtype,ig,1:ntype)
             dsdv(2,1:ntype)=dnndw0(jtype,ig,1:ntype)
             ddsdvdv(1,1:ntype)=ddnndvdv0(jtype,ig,1:ntype)
             ddsdvdv(2,1:ntype)=ddnndvdw0(jtype,ig,1:ntype)
             ddsdvdv(3,1:ntype)=ddnndwdw0(jtype,ig,1:ntype)

             tmpvecvalahat=0.d0  ! to calculate the norm of dxpos
             tmpv1val=0.d0
             do noj=1,3
                xFibvecLoc(noj)=
     &          dot_product(s(1,1:ntype),x_FibVec(1:ntype,noj))
                tmpv1val=tmpv1val+xFibvecLoc(noj)**2
                do noi=1,2
                   dxpos(noi,noj)=
     &             dot_product(dsdv(noi,1:ntype),x(1:ntype,noj))
                   tmpvecvalahat(noi)
     &                   =tmpvecvalahat(noi)+dxpos(noi,noj)**2
                end do
             enddo
             tmpvecvalahat(1)=sqrt(tmpvecvalahat(1))
             tmpvecvalahat(2)=sqrt(tmpvecvalahat(2))

             xFibvecLoc(1:3)=xFibvecLoc(1:3)/sqrt(tmpv1val)

             call crossvec(tmpv1(1:3),dxpos(1,1:3),dxpos(2,1:3))    



             hmat0(1,1)= dxpos(1,1)
             hmat0(2,1)= dxpos(1,2)
             hmat0(3,1)= dxpos(1,3)
 
             hmat0(1,2)= dxpos(2,1)
             hmat0(2,2)= dxpos(2,2)
             hmat0(3,2)= dxpos(2,3)

             hmat0(1,3)= tmpv1(1)
             hmat0(2,3)= tmpv1(2)
             hmat0(3,3)= tmpv1(3)  

             call invert(hmat0,3)


             tmpvecval=dot_product(xFibvecLoc, tmpv1)
             xFibvecLoc(1:3)=xFibvecLoc(1:3)-tmpvecval*tmpv1(1:3)

             call crossvec_noscale
     &                    (tmpv2(1:3),tmpv1(1:3),xFibvecLoc(1:3)) 

             call crossvec_noscale
     &                    (tmpvec(1:3),xFibvecLoc(1:3),dxpos(1,1:3))

            tmpv1val=0.d0
            do noj=1,3   
             tmpv1val=tmpv1val+tmpvec(noj)**2
            enddo
            tmpv1val=sqrt(tmpv1val)

  
             tmpvecval=dot_product(xFibvecLoc(1:3),dxpos(1,1:3))
             if(abs(tmpvecval).lt. 1.0d-8) then
                if (tmpvecval .lt. 0.0) then
                      tmpvecval=-1.d-8
                else
                      tmpvecval= 1.d-8
                endif
             endif
                           
             tmpv2val=atan(tmpv1val/tmpvecval)
             tmpv2val=sin(tmpv2val)
             xFibvecLoc(1:3)=xFibvecLoc(1:3)/(tmpv2val*tmpvecvalahat(2))
             tmpv2(1:3)     =tmpv2(1:3)     /(tmpv2val*tmpvecvalahat(1))

             hmat(1,1)= xFibvecLoc(1)
             hmat(2,1)= xFibvecLoc(2)
             hmat(3,1)= xFibvecLoc(3)
 
             hmat(1,2)= tmpv2(1)
             hmat(2,2)= tmpv2(2)
             hmat(3,2)= tmpv2(3)

             hmat(1,3)= tmpv1(1)
             hmat(2,3)= tmpv1(2)
             hmat(3,3)= tmpv1(3)

             hmat=matmul(hmat0,hmat)  
            crFibCrv_esh(ibdy,ne,1)=hmat(1,1)  !eta_1p
            crFibCrv_esh(ibdy,ne,2)=hmat(1,2)  !eta_2p
            crFibCrv_esh(ibdy,ne,3)=hmat(2,1)  !xi_1p
            crFibCrv_esh(ibdy,ne,4)=hmat(2,2)  !xi_2p
            enddo
            endif !Fibrous_Esh(ibdy) .eq. 1
 2121   continue            
!+++++++++++  SKIPED END+++++++++++++++++++++++++++++++++++++++
          enddo  !ne=1,numelr(ibdy)

          do ne=1,numelr(ibdy)
             jtype=neatype_esh(ibdy,ne)
             itype=jtypeall(jtype)
             if(itype .eq. 0) then
               ntype=12  
             else
               ntype=6+itype
             end if

             do nos=1,ntype
                ntem=nea_esh(ibdy,ne,nos)
                noj=1
                x(nos,noj)=gx_ibm_massive_esh(ibdy,ntem)
                noj=2
                x(nos,noj)=gy_ibm_massive_esh(ibdy,ntem)
                noj=3
                x(nos,noj)=gz_ibm_massive_esh(ibdy,ntem)

                do noj=1,2
                   sx(noj,nos)=coorcurv_esh(ibdy,ne,nos,noj)
                enddo
            enddo

             if(itype .eq. 0) then
                 xj1=x(4,1);xj2=x(7,1);xj3=x(8,1)
                 yj1=x(4,2);yj2=x(7,2);yj3=x(8,2)
                 zj1=x(4,3);zj2=x(7,3);zj3=x(8,3)
             else
                 xj1=x(1,1);xj2=x(2,1);xj3=x(itype+1,1)
                 yj1=x(1,2);yj2=x(2,2);yj3=x(itype+1,2)
                 zj1=x(1,3);zj2=x(2,3);zj3=x(itype+1,3)
             endif


             area(ibdy,ne)=
     &       (xj1*yj2 - xj2*yj1 - xj1*yj3 
     &      + xj3*yj1 + xj2*yj3 - xj3*yj2)**2
     &      +(xj1*zj2 - xj2*zj1 - xj1*zj3 
     &      + xj3*zj1 + xj2*zj3 - xj3*zj2)**2
     &      +(yj1*zj2 - yj2*zj1 - yj1*zj3 
     &      + yj3*zj1 + yj2*zj3 - yj3*zj2)**2
             area(ibdy,ne)=0.5*sqrt(area(ibdy,ne))

             tmpv1x(1)=0.333333d0*(xj1+xj2+xj3)
             tmpv1x(2)=0.333333d0*(yj1+yj2+yj3)
             tmpv1x(3)=0.333333d0*(zj1+zj2+zj3)

             xpcenter_esh(ibdy,ne,1)=sqrt((xj1-tmpv1x(1))**2
     &               +     (yj1-tmpv1x(2))**2 
     &               +     (zj1-tmpv1x(3))**2 )
            
             xpcenter_esh(ibdy,ne,2)=sqrt((xj2-tmpv1x(1))**2
     &               +     (yj2-tmpv1x(2))**2 
     &               +     (zj2-tmpv1x(3))**2 )

             xpcenter_esh(ibdy,ne,3)=sqrt((xj3-tmpv1x(1))**2
     &               +     (yj3-tmpv1x(2))**2 
     &               +     (zj3-tmpv1x(3))**2 )
            
             do ig=1,ngau

               s(1,1:ntype)=nn(jtype,ig,1:ntype)
               dsdv(1,1:ntype)=dnndv(jtype,ig,1:ntype)
               dsdv(2,1:ntype)=dnndw(jtype,ig,1:ntype)
               ddsdvdv(1,1:ntype)=ddnndvdv(jtype,ig,1:ntype)
               ddsdvdv(2,1:ntype)=ddnndvdw(jtype,ig,1:ntype)
               ddsdvdv(3,1:ntype)=ddnndwdw(jtype,ig,1:ntype)

!+++++++++++  SKIPED START+++++++++++++++++++++++++++++++++++++++
               goto 2122
!     for now skip over new implementation of Orthotropic shell
                if(Fibrous_Esh(ibdy) .eq. 1) then
                matFib_Esh(1,1)=crFibCrv_esh(ibdy,ne,1)
                matFib_Esh(1,2)=crFibCrv_esh(ibdy,ne,2)
                matFib_Esh(2,1)=crFibCrv_esh(ibdy,ne,3)
                matFib_Esh(2,2)=crFibCrv_esh(ibdy,ne,4)

                dsdv(1:2,1:ntype)=
     &              matmul(matFib_Esh(1:2,1:2),dsdv(1:2,1:ntype))


                matFib_Esh(1,1)=crFibCrv_esh(ibdy,ne,1)**2
                matFib_Esh(1,2)=2.0*crFibCrv_esh(ibdy,ne,1)
     &                           *crFibCrv_esh(ibdy,ne,2)
                matFib_Esh(1,3)=crFibCrv_esh(ibdy,ne,2)**2

                matFib_Esh(2,1)=    crFibCrv_esh(ibdy,ne,1)
     &                           *crFibCrv_esh(ibdy,ne,3)
                matFib_Esh(2,2)=   (crFibCrv_esh(ibdy,ne,1)
     &                           *crFibCrv_esh(ibdy,ne,4)
     &                           +crFibCrv_esh(ibdy,ne,2)
     &                           *crFibCrv_esh(ibdy,ne,3))
                matFib_Esh(2,3)=    crFibCrv_esh(ibdy,ne,2)
     &                           *crFibCrv_esh(ibdy,ne,4)

                matFib_Esh(3,1)=crFibCrv_esh(ibdy,ne,3)**2
                matFib_Esh(3,2)=2.0*crFibCrv_esh(ibdy,ne,3)
     &                           *crFibCrv_esh(ibdy,ne,4)
                matFib_Esh(3,3)=crFibCrv_esh(ibdy,ne,4)**2

                ddsdvdv(1:3,1:ntype)=
     &              matmul(matFib_Esh(1:3,1:3),ddsdvdv(1:3,1:ntype))
                endif
2122            continue
!+++++++++++  SKIPED END+++++++++++++++++++++++++++++++++++++++

                 do noi=1,2
                   do noj=1,2
                      jmat(noj,noi)=
     &                  dot_product(dsdv(noj,1:ntype),sx(noi,1:ntype))
                   end do
                   do noj=1,3
                      gmat(noj,noi)=
     &                dot_product(ddsdvdv(noj,1:ntype),sx(noi,1:ntype))
                   end do
                end do 
                hmat(1,1)= jmat(1,1)**2
                hmat(2,1)= jmat(1,1)*jmat(2,1)
                hmat(3,1)= jmat(2,1)**2  
                hmat(1,2)= 2.0*jmat(1,1)*jmat(1,2)
                hmat(2,2)= jmat(1,1)*jmat(2,2)+jmat(1,2)*jmat(2,1)
                hmat(3,2)= 2.0*jmat(2,1)*jmat(2,2) 
                hmat(1,3)= jmat(1,2)**2
                hmat(2,3)= jmat(1,2)*jmat(2,2)
                hmat(3,3)= jmat(2,2)**2  

                call invert(jmat,2)
                call invert(hmat,3)
  
                amat1(1:2,1:ntype)=matmul(jmat,dsdv(1:2,1:ntype))
                amat2(1:3,1:ntype)=matmul(gmat,amat1(1:2,1:ntype))
!     shape for memebrane effect
                amat2(1:3,1:ntype)
     &                    =ddsdvdv(1:3,1:ntype)-amat2(1:3,1:ntype)
                amat2(1:3,1:ntype)=matmul(hmat,amat2(1:3,1:ntype))
!     shape for bending part

                do noj=1,3
                   xpos(noj)=dot_product(s(1,1:ntype),x(1:ntype,noj))
                   do noi=1,2
                      dxpos(noi,noj)=
     &                  dot_product(amat1(noi,1:ntype),x(1:ntype,noj))
                   end do
                   do noi=1,3
                      d2xpos(noi,noj)=
     &                 dot_product(amat2(noi,1:ntype),x(1:ntype,noj))
                   end do
                end do


            if(abs(ShellModelType_esh(ibdy)) .gt. 2) then

            call crossvec(dnxpos(1:3),dxpos(1,1:3),dxpos(2,1:3))
            call crossvec_noscale(
     &            avectmp0(1:3),dxpos(1,1:3),dxpos(2,1:3))
            dncont=sqrt(dot_product( avectmp0(1:3),avectmp0(1:3)))    
            
! a_ab
            ametric(1)=dot_product(dxpos(1,1:3),dxpos(1,1:3))
            ametric(2)=dot_product(dxpos(1,1:3)  ,dxpos(2,1:3))
            ametric(3)=dot_product(dxpos(2,1:3),dxpos(2,1:3))

! b_ab
            bmetric(1)=dot_product(d2xpos(1,1:3),dnxpos(1:3))
            bmetric(2)=dot_product(d2xpos(2,1:3),dnxpos(1:3))
            bmetric(3)=dot_product(d2xpos(3,1:3),dnxpos(1:3))

            call crossvec_noscale(
     &            avectmp1(1:3),d2xpos(1,1:3),dxpos(2,1:3))
            call crossvec_noscale(
     &            avectmp2(1:3),dxpos(1,1:3),d2xpos(2,1:3))             
            avectmp3(1:3)=avectmp1(1:3)+avectmp2(1:3)
            temp_ibm3=dot_product(avectmp0(1:3),dnxpos(1:3))
            avectmp0b(1,1:3)=
     &            (avectmp3(1:3)-temp_ibm3*dnxpos(1:3))/dncont

            call crossvec_noscale(
     &            avectmp1(1:3),d2xpos(2,1:3),dxpos(2,1:3))
            call crossvec_noscale(
     &            avectmp2(1:3),dxpos(1,1:3),d2xpos(3,1:3))             
            avectmp3(1:3)=avectmp1(1:3)+avectmp2(1:3)
            temp_ibm3=dot_product(avectmp0(1:3),dnxpos(1:3))
            avectmp0b(2,1:3)=
     &            (avectmp3(1:3)-temp_ibm3*dnxpos(1:3))/dncont



            do ih=1,ngau_h
             theta3=Phi_mem_coef(i2g,1)/2.0*gau_h(ih)
             gmetric_cov(1,1)=ametric(1)-2.0*theta3*bmetric(1)
             gmetric_cov(2,2)=ametric(3)-2.0*theta3*bmetric(3)
             gmetric_cov(1,2)=ametric(2)-2.0*theta3*bmetric(2)
             gmetric_cov(2,1)=gmetric_cov(2,2)
             gmetric_con(1:2,1:2)= gmetric_cov(1:2,1:2)
             call invert(gmetric_con,2)
             gmetric_con0SAVE_esh(ibdy,ne,ig,ih,1:2,1:2)=
     &          gmetric_con(1:2,1:2)

             detC_inplane0SAVE_esh(ibdy,ne,ig,ih)=
     &        sqrt(gmetric_cov(1,1)*gmetric_cov(2,2)
     &        -gmetric_cov(1,2)**2)


             avectmp0c(1,1:3)=dxpos(1,1:3)+theta3*avectmp0b(1,1:3)
             avectmp0c(2,1:3)=dxpos(2,1:3)+theta3*avectmp0b(2,1:3)

             gbase_con0SAVE_esh(ibdy,ne,i,ih,1,1:3)
     &                       =gmetric_con(1,1)*avectmp0c(1,1:3)
     &                       +gmetric_con(1,2)*avectmp0c(2,1:3)

             gbase_con0SAVE_esh(ibdy,ne,i,ih,2,1:3)
     &                       =gmetric_con(2,1)*avectmp0c(1,1:3)
     &                       +gmetric_con(2,2)*avectmp0c(2,1:3)
            enddo

            else !othe shell models

             noi=1
             noj=1
             dmat0_esh(ibdy,ne,ig,noi,noj)=1/2.0*
     &        (dot_product(dxpos(noi,1:3),dxpos(noj,1:3))
     &        )
             kmatbar_esh(noi,noj)=sqrt(
     &        dot_product(d2xpos(noi,1:3),d2xpos(noj,1:3)))
                
             noi=1
             noj=2
             dmat0_esh(ibdy,ne,ig,noi,noj)=1/2.0*
     &        (dot_product(dxpos(noi,1:3),dxpos(noj,1:3))
     &        )
             kmatbar_esh(noi,noj)=sqrt(
     &        dot_product(d2xpos(noj,1:3),d2xpos(noj,1:3)))
             kmatbar_esh(noj,noi)=kmatbar_esh(noi,noj)
             dmat0_esh(ibdy,ne,ig,noj,noi)=dmat0_esh(ibdy,ne,ig,noi,noj)
             noi=2
             noj=2
             dmat0_esh(ibdy,ne,ig,noi,noj)=1/2.0*
     &        (dot_product(dxpos(noi,1:3),dxpos(noj,1:3))
     &         )
             kmatbar_esh(noi,noj)=sqrt(dot_product(d2xpos(noj+1,1:3)
     &        ,d2xpos(noj+1,1:3)))



             kmat0_esh(ibdy,ne,ig,1,1)=kmatbar_esh(1,1)
             kmat0_esh(ibdy,ne,ig,1,2)=kmatbar_esh(1,2)
             kmat0_esh(ibdy,ne,ig,2,1)=kmatbar_esh(2,1)
             kmat0_esh(ibdy,ne,ig,2,2)=kmatbar_esh(2,2)

            endif
          end do  ! gauss
      end do


      if(contacttype .le. 1) then
!     old contact model => to find wether two points are initially close
      do np=1,ns_ibm_r_esh(ibdy)  
           do noj=np,ns_ibm_r_esh(ibdy)  
                           ContactShellFlag_esh(ibdy,np,noj)=0
                           ContactShellFlag_esh(ibdy,noj,np)=0

                if((
     &              abs((gy_ibm_massive_esh(ibdy,np)
     &                  -gy_ibm_massive_esh(ibdy,noj))/min_grid_y)
     &              .le.13.0) .and. (
     &              abs((gx_ibm_massive_esh(ibdy,np)
     &                  -gx_ibm_massive_esh(ibdy,noj))/min_grid_x)
     &              .le.13.0) .and.  (
     &              abs((gz_ibm_massive_esh(ibdy,np)
     &                  -gz_ibm_massive_esh(ibdy,noj))/min_grid_z)
     &              .le.13.0) ) then
                           ContactShellFlag_esh(ibdy,np,noj)=1
                           ContactShellFlag_esh(ibdy,noj,np)=1
                 endif
            enddo
       enddo
       elseif(contacttype.gt. 1) then !nonlinear/simpler contact model
       jq=1
        do np=1,numelr(ibdy) 
               jtype2=neatype_esh(ibdy,np)
               itype2=jtypeall(jtype2)
               if(itype2 .eq. 0) then
                   ntype2=12  
               else
                   ntype2=6+itype2
               end if   
               if(itype2 .eq. 0) then
                  ItypecontELE(1)=4
                  ItypecontELE(2)=7
                  ItypecontELE(3)=8
               else              
                  ItypecontELE(1)=1
                  ItypecontELE(2)=2
                  ItypecontELE(3)=itype2+1
               endif
               do nos2=1,3
                   ntem2=nea_esh(ibdy,np,ItypecontELE(nos2))
                    noj=1
                    Xcont(noj,nos2)=gx_ibm_massive_esh(ibdy,ntem2)
                    noj=2
                    xcont(noj,nos2)=gy_ibm_massive_esh(ibdy,ntem2)
                    noj=3
                    xcont(noj,nos2)=gz_ibm_massive_esh(ibdy,ntem2)
               enddo
               do j=1,3
                    XP1cont(j)=(xcont(j,1)+xcont(j,2)+xcont(j,3))/3.0d0
               enddo

           do noj=np,numelr(ibdy)  
               ContactShellFlag_esh(ibdy,np,noj)=0
               ContactShellFlag_esh(ibdy,noj,np)=0

               jtype2=neatype_esh(ibdy,noj)
               itype2=jtypeall(jtype2)
               if(itype2 .eq. 0) then
                   ntype2=12  
               else
                   ntype2=6+itype2
               end if   
               if(itype2 .eq. 0) then
                  ItypecontELE(1)=1


                  ItypecontELE(2)=2
                  ItypecontELE(3)=3
               else              
                  ItypecontELE(1)=1
                  ItypecontELE(2)=2
                  ItypecontELE(3)=itype2+1
               endif
               do nos2=1,3
                   ntem2=nea_esh(ibdy,noj,ItypecontELE(nos2))
                    j=1
                    Xcont2(j,nos2)=gx_ibm_massive_esh(ibdy,ntem2)
                    j=2
                    xcont2(j,nos2)=gy_ibm_massive_esh(ibdy,ntem2)
                    j=3
                    xcont2(j,nos2)=gz_ibm_massive_esh(ibdy,ntem2)
               enddo
               do j=1,3
                    XP2cont(j)=
     &                     (xcont2(j,1)+xcont2(j,2)+xcont2(j,3))/3.0d0
               enddo

              dxcont(1:3)=XP2cont(1:3)-XP1cont(1:3)
              ddcont= sqrt(dot_product(dxcont(1:3),dxcont(1:3)))

            
              if(ddcont .lt. max(Ccont(9)*Ccont(1),contactlength))then
                    ContactShellFlag_esh(ibdy,np,noj)=1
                    ContactShellFlag_esh(ibdy,noj,np)=1
              endif
            enddo
       enddo
       else
           stop 'there is no more if here ???'
       endif

        if(ndim .gt. 2) then
          ddcont=(min_grid_x*min_grid_y*min_grid_z)**(1.0/3.0)
        else
          ddcont=sqrt(min_grid_x*min_grid_y)
        endif
       
        do ibkt2=1,Ncontactplanes(i2g) 
                contactplane_point(1:3)= 
     &                       contactplanesPara(i2g,ibkt2,1:3)
                contactplane_n(1:3)=
     &                       contactplanesPara(i2g,ibkt2,4:6)

                dncont=0.0
                do j=1,3
                   dncont=dncont+contactplane_n(j) 
     &                          *contactplane_point(j)
                enddo
                contactplane_n(4)=-dncont
!     plane defined as ax+by+cz+d=0, here we save d

                do j=1,ns_ibm_r_esh(ibdy)
                  ksh_tmp=contactplane_n(1) 
     &                   *gx_ibm_massive_esh(ibdy,j) 
     &                   +contactplane_n(2) 
     &                   *gy_ibm_massive_esh(ibdy,j) 
     &                   +contactplane_n(3) 
     &                   *gz_ibm_massive_esh(ibdy,j)
     &                   +contactplane_n(4) 

                  if(ksh_tmp .ge. impactpen(3)*ddcont) then
                     contact_coef_esh(ibdy,ibkt2,j)=impactpen(3)*ddcont
                  elseif(ksh_tmp .ge. 0.05*impactpen(3)*ddcont) then
                     contact_coef_esh(ibdy,ibkt2,j)=ksh_tmp
                  else
                     contact_coef_esh(ibdy,ibkt2,j)= 
     &                0.05*impactpen(3)*ddcont
                  endif
                enddo 
  
         enddo
      enddo !nr_ibm_esh
#endif      
#ifdef INCLUDE_FBC
!!======== Perform intial calculations for fabric shells
        do i=1,nr_ibm_fbc
           ibdy=i
           in_fbc=ibdy
           i2g=Ilocglb_fbc(ibdy)
! for free boundary condition
           jreal_2_imag_fbc=0.d0
           do np=1,ns_ibm_r_fbc(ibdy)  
              jreal_2_imag_fbc(ibdy,np,np)=1.0
           end do       
           do np=1,ns_ibm_i_fbc(ibdy)  
              npg=ns_ibm_r_fbc(ibdy)+np
              bodypoint(1)=boundrylist_fbc(ibdy,np,1)
              bodypoint(2)=boundrylist_fbc(ibdy,np,2)
              bodypoint(3)=boundrylist_fbc(ibdy,np,3)
           gx_ibm_massive_fbc(ibdy,npg)=
     &                     gx_ibm_massive_fbc(ibdy,bodypoint(1))
     &                    +gx_ibm_massive_fbc(ibdy,bodypoint(2))       
     &                    -gx_ibm_massive_fbc(ibdy,bodypoint(3))       




           gy_ibm_massive_fbc(ibdy,npg)=
     &                     gy_ibm_massive_fbc(ibdy,bodypoint(1))
     &                    +gy_ibm_massive_fbc(ibdy,bodypoint(2))       
     &                    -gy_ibm_massive_fbc(ibdy,bodypoint(3))       
    
           gz_ibm_massive_fbc(ibdy,npg)=
     &                     gz_ibm_massive_fbc(ibdy,bodypoint(1))
     &                    +gz_ibm_massive_fbc(ibdy,bodypoint(2))       
     &                    -gz_ibm_massive_fbc(ibdy,bodypoint(3))

              jreal_2_imag_fbc(ibdy,npg,bodypoint(1))=1.0
              jreal_2_imag_fbc(ibdy,npg,bodypoint(2))=1.0
              jreal_2_imag_fbc(ibdy,npg,bodypoint(3))=-1.0


           gx_bp_fbc0(ibdy,npg)=
     &                     gx_bp_fbc0(ibdy,bodypoint(1))
     &                    +gx_bp_fbc0(ibdy,bodypoint(2))       
     &                    -gx_bp_fbc0(ibdy,bodypoint(3))       

           gy_bp_fbc0(ibdy,npg)=
     &                     gy_bp_fbc0(ibdy,bodypoint(1))
     &                    +gy_bp_fbc0(ibdy,bodypoint(2))       
     &                    -gy_bp_fbc0(ibdy,bodypoint(3))       
    
           gz_bp_fbc0(ibdy,npg)=
     &                     gz_bp_fbc0(ibdy,bodypoint(1))
     &                    +gz_bp_fbc0(ibdy,bodypoint(2))       
     &                    -gz_bp_fbc0(ibdy,bodypoint(3))
           end do       

!update center of elements
               do j=1,numelr(i2)   
                  EleCenter0_fbc(in_fbc,j,1:3)=0.0           
                  do j2=1,3  !nodes
                     k=nea_fbc(in_fbc,j,j2)
                     EleCenter0_fbc(in_fbc,j,1)=
     &                   EleCenter0_fbc(in_fbc,j,1)
     &                   +gx_ibm_fbc(in_fbc,k)/3.0
                     EleCenter0_fbc(in_fbc,j,2)=
     &                   EleCenter0_fbc(in_fbc,j,2)
     &                   +gy_ibm_fbc(in_fbc,k)/3.0
                     EleCenter0_fbc(in_fbc,j,3)=
     &                   EleCenter0_fbc(in_fbc,j,3)
     &                   +gz_ibm_fbc(in_fbc,k)/3.0
                  enddo
                 enddo
                  do j=1,numelr(i2) 
                  do j2=1,3  !nodes
                     k=nea_fbc(in_fbc,j,j2)
                     xpcenter_fbc(in_fbc,j,j2)=sqrt(
     &              (gx_ibm_fbc(in_fbc,k)-EleCenter0_fbc(in_fbc,j,1))**2
     &             +(gy_ibm_fbc(in_fbc,k)-EleCenter0_fbc(in_fbc,j,2))**2
     &             +(gz_ibm_fbc(in_fbc,k)-EleCenter0_fbc(in_fbc,j,3))**2
     &               )
                  enddo
                 enddo            
!update area of element
               do j=1,numelr(i2)   
                 iv0=nea_fbc(in_fbc,j,1)
                 iv1=nea_fbc(in_fbc,j,2)
                 iv2=nea_fbc(in_fbc,j,3)
                 axtmp_fbc=gx_ibm_fbc(in_fbc,iv1)-gx_ibm_fbc(in_fbc,iv0)
                 aytmp_fbc=gy_ibm_fbc(in_fbc,iv1)-gy_ibm_fbc(in_fbc,iv0)
                 aztmp_fbc=gz_ibm_fbc(in_fbc,iv1)-gz_ibm_fbc(in_fbc,iv0)
                 bxtmp_fbc=gx_ibm_fbc(in_fbc,iv2)-gx_ibm_fbc(in_fbc,iv0)
                 bytmp_fbc=gy_ibm_fbc(in_fbc,iv2)-gy_ibm_fbc(in_fbc,iv0)
                 bztmp_fbc=gz_ibm_fbc(in_fbc,iv2)-gz_ibm_fbc(in_fbc,iv0)
                 cxtmp_fbc=aytmp_fbc*bztmp_fbc-aztmp_fbc*bytmp_fbc
                 cytmp_fbc=aztmp_fbc*bxtmp_fbc-axtmp_fbc*bztmp_fbc
                 cztmp_fbc=axtmp_fbc*bytmp_fbc-aytmp_fbc*bxtmp_fbc
                 EleArea0_fbc(in_fbc,j)=
     &                0.5*sqrt(cxtmp_fbc**2+cytmp_fbc**2+cztmp_fbc**2)
                 areapnt_fbc(ibdy,iv0)=
     &                  areapnt_fbc(ibdy,iv0)
     &                                   +EleArea0_fbc(in_fbc,j)/3.0
                 areapnt_fbc(ibdy,iv1)=
     &                  areapnt_fbc(ibdy,iv1)
     &                                   +EleArea0_fbc(in_fbc,j)/3.0
                 areapnt_fbc(ibdy,iv2)=
     &                  areapnt_fbc(ibdy,iv2)
     &                                   +EleArea0_fbc(in_fbc,j)/3.0

!compute the normal
                 axtmp_fbc=gx_ibm_fbc(in_fbc,iv0)-gx_ibm_fbc(in_fbc,iv2)
                 aytmp_fbc=gy_ibm_fbc(in_fbc,iv0)-gy_ibm_fbc(in_fbc,iv2)
                 aztmp_fbc=gz_ibm_fbc(in_fbc,iv0)-gz_ibm_fbc(in_fbc,iv2)
                 bxtmp_fbc=gx_ibm_fbc(in_fbc,iv0)-gx_ibm_fbc(in_fbc,iv1)
                 bytmp_fbc=gy_ibm_fbc(in_fbc,iv0)-gy_ibm_fbc(in_fbc,iv1)
                 bztmp_fbc=gz_ibm_fbc(in_fbc,iv0)-gz_ibm_fbc(in_fbc,iv1)
                 tmpcoord(1)=aytmp_fbc*bztmp_fbc-aztmp_fbc*bytmp_fbc
                 tmpcoord(2)=aztmp_fbc*bxtmp_fbc-axtmp_fbc*bztmp_fbc
                 tmpcoord(3)=axtmp_fbc*bytmp_fbc-aytmp_fbc*bxtmp_fbc
                 cxtmp_fbc
     &                =sqrt(dot_product(tmpcoord(1:3),tmpcoord(1:3)))
                 tmpcoord(1)=tmpcoord(1)/cxtmp_fbc
                 tmpcoord(2)=tmpcoord(2)/cxtmp_fbc
                 tmpcoord(3)=tmpcoord(3)/cxtmp_fbc 
                 EleNormal0_fbc(in_fbc,j,1:3)= tmpcoord(1:3) 

!compute the inner element angles and edge sizes
               tmpcoord(1)=gx_ibm_fbc(in_fbc,iv2)-gx_ibm_fbc(in_fbc,iv1)
               tmpcoord(2)=gy_ibm_fbc(in_fbc,iv2)-gy_ibm_fbc(in_fbc,iv1)
               tmpcoord(3)=gz_ibm_fbc(in_fbc,iv2)-gz_ibm_fbc(in_fbc,iv1)

                 EleEdge0_fbc(in_fbc,j,1) 
     &                =sqrt(dot_product(tmpcoord(1:3),tmpcoord(1:3)))

              tmpcoord2(1)=gx_ibm_fbc(in_fbc,iv2)-gx_ibm_fbc(in_fbc,iv0)
              tmpcoord2(2)=gy_ibm_fbc(in_fbc,iv2)-gy_ibm_fbc(in_fbc,iv0)
              tmpcoord2(3)=gz_ibm_fbc(in_fbc,iv2)-gz_ibm_fbc(in_fbc,iv0)

                 EleEdge0_fbc(in_fbc,j,2) 
     &                =sqrt(dot_product(tmpcoord2(1:3),tmpcoord2(1:3)))

              tmpcoord3(1)=gx_ibm_fbc(in_fbc,iv1)-gx_ibm_fbc(in_fbc,iv0)
              tmpcoord3(2)=gy_ibm_fbc(in_fbc,iv1)-gy_ibm_fbc(in_fbc,iv0)
              tmpcoord3(3)=gz_ibm_fbc(in_fbc,iv1)-gz_ibm_fbc(in_fbc,iv0)

                 EleEdge0_fbc(in_fbc,j,3) 
     &                =sqrt(dot_product(tmpcoord3(1:3),tmpcoord3(1:3)))

                 if(inplane_type_fbc(in_fbc) .ne. 1) then

                 axtmp_fbc
     &             =acos(dot_product(tmpcoord3(1:3),tmpcoord2(1:3))
     &             /(EleEdge0_fbc(in_fbc,j,3)*EleEdge0_fbc(in_fbc,j,2)))

                 aytmp_fbc
     &             =acos(dot_product(tmpcoord(1:3),-tmpcoord3(1:3))
     &             /(EleEdge0_fbc(in_fbc,j,1)*EleEdge0_fbc(in_fbc,j,3)))

                 aztmp_fbc
     &             =acos(dot_product(-tmpcoord(1:3),-tmpcoord2(1:3))
     &             /(EleEdge0_fbc(in_fbc,j,1)*EleEdge0_fbc(in_fbc,j,2)))

                 bxtmp_fbc
     &             =inplane_E_fbc(in_fbc)/(
     &               16.0*(1.0-inplane_nu_fbc(in_fbc)**2)
     &               *EleArea0_fbc(in_fbc,j)
     &               )

                    EleInplaneParaK_fbc(in_fbc,j,1)
     &              =bxtmp_fbc
     &              *(2*(
     &              1.0/tan(axtmp_fbc))**2+1.0-inplane_nu_fbc(in_fbc))
                                        
                    EleInplaneParaK_fbc(in_fbc,j,2)
     &              =bxtmp_fbc
     &              *(2*(
     &              1.0/tan(aytmp_fbc))**2+1.0-inplane_nu_fbc(in_fbc))

                    EleInplaneParaK_fbc(in_fbc,j,3)
     &              =bxtmp_fbc
     &              *(2*(
     &              1.0/tan(aztmp_fbc))**2+1.0-inplane_nu_fbc(in_fbc))

                    EleInplaneParaC_fbc(in_fbc,j,1)
     &              =bxtmp_fbc
     &              *(2*
     &              1.0/tan(aytmp_fbc)*
     &              1.0/tan(aztmp_fbc)+1.0
     &                -inplane_nu_fbc(in_fbc))
                                        
                    EleInplaneParaC_fbc(in_fbc,j,2)
     &              =bxtmp_fbc
     &              *(2*
     &              1.0/tan(axtmp_fbc)*
     &              1.0/tan(aztmp_fbc)+1.0
     &                -inplane_nu_fbc(in_fbc))

                    EleInplaneParaC_fbc(in_fbc,j,3)
     &              =bxtmp_fbc
     &              *(2*
     &              1.0/tan(axtmp_fbc)*
     &              1.0/tan(aytmp_fbc)+1.0
     &                -inplane_nu_fbc(in_fbc))
                 else
                   axtmp_fbc=(  EleEdge0_fbc(in_fbc,j,1)
     &                         +EleEdge0_fbc(in_fbc,j,2)
     &                         +EleEdge0_fbc(in_fbc,j,3)) /3.0

                   if(edge_length_max_fbc(in_fbc) .lt. 0) then
                      bxtmp_fbc=abs(edge_length_max_fbc(in_fbc))
     &                          *axtmp_fbc
                   else
                      bxtmp_fbc=edge_length_max_fbc(in_fbc)
                   endif
                   cxtmp_fbc=axtmp_fbc/bxtmp_fbc
                   cytmp_fbc=EleArea0_fbc(in_fbc,j)
                   if(inplane_Cq_fbc(in_fbc,1)<0.0) then
                      EleInplane_Cq_fbc(in_fbc,j,1)=
     &                         sqrt(3.0)
     &                        *cytmp_fbc**(inplane_Cq_fbc(in_fbc,3)+1)
     &                        *kappa_inplane_fbc(in_fbc,1)
     &                        *KBT_fbc(in_fbc)
     &                        *(4.0*cxtmp_fbc**2-9.0*cxtmp_fbc+6.0)
     &                        /(
     &                          4.0*edge_length_p_fbc(in_fbc)
     &                          *inplane_Cq_fbc(in_fbc,3)
     &                          *bxtmp_fbc
     &                          *(1.0-cxtmp_fbc)**2
     &                          )
                    else
                      EleInplane_Cq_fbc(in_fbc,j,1)
     &                         =inplane_Cq_fbc(in_fbc,1)
                    endif
                   if(inplane_Cq_fbc(in_fbc,2)<0.0) then
                      EleInplane_Cq_fbc(in_fbc,j,2)=
     &                         sqrt(3.0)
     &                        *cytmp_fbc**(inplane_Cq_fbc(in_fbc,4)+1)
     &                        *kappa_inplane_fbc(in_fbc,2)
     &                        /(
     &                           inplane_Cq_fbc(in_fbc,4)
     &                          *(1.0-cxtmp_fbc**2)
     &                          )
                    else
                      EleInplane_Cq_fbc(in_fbc,j,2)
     &                         =inplane_Cq_fbc(in_fbc,2)
                    endif
                 endif
               enddo
             if(bending_flag_fbc(in_fbc).and.tet0_fbc(in_fbc))then
               do j=1,numeFacePair(i2)
                  ie1=facepairEle_fbc(in_fbc,j,1)
                  ie2=facepairEle_fbc(in_fbc,j,2) 
                  tmpcoord(1:3) =EleNormal0_fbc(in_fbc,ie1,1:3) 
     &                          -EleNormal0_fbc(in_fbc,ie2,1:3)
                  tmpcoord2(1:3)=EleCenter0_fbc(in_fbc,ie1,1:3) 
     &                          -EleCenter0_fbc(in_fbc,ie2,1:3)
                  axtmp_fbc=dot_product(tmpcoord(1:3),tmpcoord2(1:3))
                  if(axtmp_fbc>0) then
	                    sinsign=1.0
                  else
	                    sinsign=-1.0
                  endif                         
                  aytmp_fbc=dot_product(EleNormal0_fbc(in_fbc,ie1,1:3) 
     &                                 ,EleNormal0_fbc(in_fbc,ie2,1:3)) 
                  aztmp_fbc=acos(aytmp_fbc)
                  if(sinsign<0) aztmp_fbc=2.0*pi1-aztmp_fbc  
                  facepairtet0_fbc(in_fbc,j)=aztmp_fbc
               enddo
             else
                if (tet0_fbc(in_fbc)) then
                   facepairtet0_fbc(in_fbc,1:numeFacePair(i2))=1.d0
                else
                   facepairtet0_fbc(in_fbc,1:numeFacePair(i2))=0.d0
                end if
             endif               

             if(volume_flag_fbc(in_fbc) 
     &          .and. (volume_eq_fbc(in_fbc).lt. 0)) then
                volume_eq_fbc(in_fbc)=0.0
!calculate total volume
                do j=1,numelr(i2) 
                   volume_eq_fbc(in_fbc)=volume_eq_fbc(in_fbc)
     &                   +EleArea0_fbc(in_fbc,j)*(
     &                            EleNormal0_fbc(in_fbc,j,1)           
     &                           *EleCenter0_fbc(in_fbc,j,1) 
     &                           +EleNormal0_fbc(in_fbc,j,2)           
     &                           *EleCenter0_fbc(in_fbc,j,2) 
     &                           +EleNormal0_fbc(in_fbc,j,3)           
     &                           *EleCenter0_fbc(in_fbc,j,3) 
     &                           )
               enddo
               volume_eq_fbc(in_fbc)=volume_eq_fbc(in_fbc)/3.0
            endif 
             if(Aera_flag_fbc(in_fbc) 
     &          .and. (Area_eq_fbc(in_fbc).lt. 0)) then
                Area_eq_fbc(in_fbc)=0.0
                do j=1,numelr(i2)           
                   Area_eq_fbc(in_fbc)=Area_eq_fbc(in_fbc)
     &                   +EleArea0_fbc(in_fbc,j)
               enddo
            endif 
            if(contacttype .le. 1) then
!     old contact model => to find wether two points are initially close
      do np=1,ns_ibm_r_fbc(ibdy)  
           do noj=np,ns_ibm_r_fbc(ibdy)  
                           ContactShellFlag_fbc(ibdy,np,noj)=0
                           ContactShellFlag_fbc(ibdy,noj,np)=0

                if((
     &              abs((gy_ibm_massive_fbc(ibdy,np)
     &                  -gy_ibm_massive_fbc(ibdy,noj))/min_grid_y)
     &              .le.13.0) .and. (
     &              abs((gx_ibm_massive_fbc(ibdy,np)
     &                  -gx_ibm_massive_fbc(ibdy,noj))/min_grid_x)
     &              .le.13.0) .and.  (
     &              abs((gz_ibm_massive_fbc(ibdy,np)
     &                  -gz_ibm_massive_fbc(ibdy,noj))/min_grid_z)
     &              .le.13.0) ) then
                           ContactShellFlag_fbc(ibdy,np,noj)=1
                           ContactShellFlag_fbc(ibdy,noj,np)=1
                 endif
            enddo
       enddo
       elseif(contacttype.gt. 1) then !nonlinear/simpler contact model
       jq=1
        do np=1,numelr(ibdy) 
               do nos2=1,3
                   ntem2=nea_fbc(ibdy,np,nos2)
                    noj=1
                    Xcont(noj,nos2)=gx_ibm_massive_fbc(ibdy,ntem2)
                    noj=2
                    xcont(noj,nos2)=gy_ibm_massive_fbc(ibdy,ntem2)
                    noj=3
                    xcont(noj,nos2)=gz_ibm_massive_fbc(ibdy,ntem2)
               enddo
               do j=1,3
                    XP1cont(j)=(xcont(j,1)+xcont(j,2)+xcont(j,3))/3.0d0
               enddo

           do noj=np,numelr(ibdy)  
               ContactShellFlag_fbc(ibdy,np,noj)=0
               ContactShellFlag_fbc(ibdy,noj,np)=0

               do nos2=1,3
                   ntem2=nea_fbc(ibdy,noj,nos2)
                    j=1
                    Xcont2(j,nos2)=gx_ibm_massive_fbc(ibdy,ntem2)
                    j=2
                    xcont2(j,nos2)=gy_ibm_massive_fbc(ibdy,ntem2)
                    j=3
                    xcont2(j,nos2)=gz_ibm_massive_fbc(ibdy,ntem2)
               enddo
               do j=1,3
                    XP2cont(j)=
     &                     (xcont2(j,1)+xcont2(j,2)+xcont2(j,3))/3.0d0
               enddo

              dxcont(1:3)=XP2cont(1:3)-XP1cont(1:3)
              ddcont= sqrt(dot_product(dxcont(1:3),dxcont(1:3)))

            
              if(ddcont .lt. max(Ccont(9)*Ccont(1),contactlength))then
                    ContactShellFlag_fbc(ibdy,np,noj)=1
                    ContactShellFlag_fbc(ibdy,noj,np)=1
              endif
            enddo
       enddo
       else
           stop 'there is no more if here ???'
       endif

        if(ndim .gt. 2) then
          ddcont=(min_grid_x*min_grid_y*min_grid_z)**(1.0/3.0)
        else
          ddcont=sqrt(min_grid_x*min_grid_y)
        endif
       
        do ibkt2=1,Ncontactplanes(i2g) 
                contactplane_point(1:3)= 
     &                       contactplanesPara(i2g,ibkt2,1:3)
                contactplane_n(1:3)=
     &                       contactplanesPara(i2g,ibkt2,4:6)

                dncont=0.0
                do j=1,3
                   dncont=dncont+contactplane_n(j) 
     &                          *contactplane_point(j)
                enddo
                contactplane_n(4)=-dncont
!     plane defined as ax+by+cz+d=0, here we save d

                do j=1,ns_ibm_r_fbc(ibdy)
                  ksh_tmp=contactplane_n(1) 
     &                   *gx_ibm_massive_fbc(ibdy,j) 
     &                   +contactplane_n(2) 
     &                   *gy_ibm_massive_fbc(ibdy,j) 
     &                   +contactplane_n(3) 
     &                   *gz_ibm_massive_fbc(ibdy,j)
     &                   +contactplane_n(4) 

                  if(ksh_tmp .ge. impactpen(3)*ddcont) then
                     contact_coef_fbc(ibdy,ibkt2,j)=impactpen(3)*ddcont
                  elseif(ksh_tmp .ge. 0.05*impactpen(3)*ddcont) then
                     contact_coef_fbc(ibdy,ibkt2,j)=ksh_tmp
                  else
                     contact_coef_fbc(ibdy,ibkt2,j)= 
     &                0.05*impactpen(3)*ddcont
                  endif
                enddo 
         enddo
      enddo !nr_ibm_fbc

!!======== Calculation of intersection between line/mesh and fbc membrane
      do i=1,nr_ibm_fbc
           ibdy=i
           in_fbc=ibdy
           i2g=Ilocglb_fbc(ibdy)
           if(FabricMeshFlag_fbc(i)) then
             nQPoints=1

             do k=1,FabricMesh_nMesh(i)              
               do j=1,FabricMesh_npoint(i,k)
                  jq=FabricMesh_address(i,k,j)
                  qPoints(1,1:3)=FabricMesh_coord(i,jq,1:3)
                  call  point2trimesh(
     &                 DD
     &                ,PP
     &                ,BB
     &                ,FF
     &                ,numelr(i2g),numedge(in_fbc)
     &                ,ns_ibm_r_fbc(i2g),nQPoints
     &                ,nea_fbc(in_fbc,1:numelr(i2g),1:3)
     &                ,EdgeList_fbc(in_fbc,1:numedge(in_fbc),1:2)
     &                ,gx_ibm_massive_fbc(in_fbc,1:ns_ibm_r_fbc(i2g))
     &                ,gy_ibm_massive_fbc(in_fbc,1:ns_ibm_r_fbc(i2g))
     &                ,gz_ibm_massive_fbc(in_fbc,1:ns_ibm_r_fbc(i2g))
     &                ,EleNormal0_fbc(i,1:numelr(i2g),1:3),qPoints) 

                  FabricMesh_SurfElAddress(i,jq) =FF(1)
                  FabricMesh_SurfElAttr(i,jq,1  )=abs(DD(1))
                  FabricMesh_SurfElAttr(i,jq,2:4)=BB(1,1:3) 
               enddo
             enddo
          endif                 
          if(FabricLineFlag_fbc(i)) then
             nQPoints=1

             do k=1,FabricLine_nLine(i)              
               do j=1,FabricLine_npoint(i,k)
                  jq=FabricLine_address(i,k,j)
                  qPoints(1,1:3)=FabricLine_coord(i,jq,1:3)
                  call  point2trimesh(
     &                 DD
     &                ,PP
     &                ,BB
     &                ,FF
     &                ,numelr(i2g),numedge(in_fbc)
     &                ,ns_ibm_r_fbc(i2g),nQPoints
     &                ,nea_fbc(in_fbc,1:numelr(i2g),1:3)
     &                ,EdgeList_fbc(in_fbc,1:numedge(in_fbc),1:2)
     &                ,gx_ibm_massive_fbc(in_fbc,1:ns_ibm_r_fbc(i2g))
     &                ,gy_ibm_massive_fbc(in_fbc,1:ns_ibm_r_fbc(i2g))
     &                ,gz_ibm_massive_fbc(in_fbc,1:ns_ibm_r_fbc(i2g))
     &                ,EleNormal0_fbc(i,1:numelr(i2g),1:3),qPoints) 

                  FabricLine_SurfElAddress(i,jq) =FF(1)
                  FabricLine_SurfElAttr(i,jq,1  )=abs(DD(1))
                  FabricLine_SurfElAttr(i,jq,2:4)=BB(1,1:3) 
               enddo
             enddo
          endif                 
      enddo !nr_ibm_fbc
#endif
#ifdef INCLUDE_FSH
!!====== calculate initial stretch in the structure for later use.
      do i=1,nr_ibm_fsh
           i2g=Ilocglb_fsh(i)
           ibdy=i
           do jq=1,nq_ibm_r_fsh(i)
           do j=1,ns_ibm_r_fsh(i)
               GXtmp(1,jq,j)=gx_ibm_massive_fsh(i,jq,j)
               GXtmp(2,jq,j)=gy_ibm_massive_fsh(i,jq,j)
               GXtmp(3,jq,j)=gz_ibm_massive_fsh(i,jq,j)
           enddo
           enddo

           do jq=1,nq_ibm_r_fsh(i)
            if (BoundaryFibShell_fsh(i,1,1) .eq. 3) then
             do j=1-(Ns_IBMB_fsh+1),0,1
               jbndy=ns_ibm_r_fsh(i)+j
               GXtmp(1:3,jq,j)=GXtmp(1:3,jq,jbndy)
             enddo
            else
             do j=0,1-(1+Ns_IBMB_fsh),-1

               jbndy=j+1
               GXtmp(1:3,jq,j)=
     &            2.0*GXtmp(1:3,jq,jbndy)-GXtmp(1:3,jq,jbndy+1)
             enddo 
            endif
            if (BoundaryFibShell_fsh(i,1,2) .eq. 3) then
             do j=ns_ibm_r_fsh(i)+1,ns_ibm_r_fsh(i)+(1+Ns_IBMB_fsh),1
               jbndy=-ns_ibm_r_fsh(i)+j
               GXtmp(1:3,jq,j)=GXtmp(1:3,jq,jbndy)
             enddo
            else
             do j=ns_ibm_r_fsh(i)+1,ns_ibm_r_fsh(i)+(1+Ns_IBMB_fsh),1
               jbndy=j-1
               GXtmp(1:3,jq,j)=
     &            2.0*GXtmp(1:3,jq,jbndy)-GXtmp(1:3,jq,jbndy-1)
             enddo 
            endif
          enddo


           do j=1-(1+Ns_IBMB_fsh),ns_ibm_r_fsh(i)+(1+Ns_IBMB_fsh)
            if (BoundaryFibShell_fsh(i,2,1) .eq. 3) then
             do jq=1-(Nq_IBMB_fsh+1),0,1
               jqbndy=nq_ibm_r_fsh(i)+jq
               GXtmp(1:3,jq,j)=GXtmp(1:3,jqbndy,j)
             enddo
            else
             do jq=0,1-(Nq_IBMB_fsh+1),-1
               jqbndy=jq+1
               GXtmp(1:3,jq,j)=
     &            2.0*GXtmp(1:3,jqbndy,j)-GXtmp(1:3,jqbndy+1,j)
             enddo 
            endif
            if (BoundaryFibShell_fsh(i,2,2) .eq. 3) then
             do jq=nq_ibm_r_fsh(i)+1,nq_ibm_r_fsh(i)+(Nq_IBMB_fsh+1),1
               jqbndy=-nq_ibm_r_fsh(i)+jq
               GXtmp(1:3,jq,j)=GXtmp(1:3,jqbndy,j)
             enddo
            else
             do jq=nq_ibm_r_fsh(i)+1,nq_ibm_r_fsh(i)+(Nq_IBMB_fsh+1),1
               jqbndy=jq-1
               GXtmp(1:3,jq,j)=
     &            2.0*GXtmp(1:3,jqbndy,j)-GXtmp(1:3,jqbndy-1,j)
             enddo 
            endif
          enddo


           do jq=1,nq_ibm_r_fsh(i)
           do j=1,ns_ibm_r_fsh(i)

             if(e4coef(i2g) .ge. 0) then
              DX_ip_j(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))/dsf_IBM_fsh(i,jq,j)
!x+1/2,y

              DX_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))/dsf_IBM_fsh(i,jq-1,j)
!     x-1/2,y

              DX_i_jp(1:3)=
     &         (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq,j+1))
     &             /dsf_IBM_fsh(i,jq,j+1)
!     x+1/2,y+1

              DX_im_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq-1,j+1))
     &             /dsf_IBM_fsh(i,jq-1,j+1)
!     x-1/2,y+1

              DX_i_jm(1:3)=0.25*(
     &             DX_i_jp(1:3)+DX_im_jp(1:3)+DX_ip_j(1:3)+DX_im_j(1:3))
!     x,y+1/2


              DY_i_jp(1:3)=
     &             (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &             /dsf2_IBM_fsh(i,jq,j)
!     x,y+1/2

              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &             /dsf2_IBM_fsh(i,jq,j-1)
!     x,y-1/2

              DY_ip_j(1:3)=
     &        (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq+1,j))
     &             /dsf2_IBM_fsh(i,jq+1,j)
!     x+1,y+1/2

              DY_ip_jm(1:3)=
     &        (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq+1,j-1))
     &             /dsf2_IBM_fsh(i,jq+1,j-1)
!     x+1,y-1/2 

              DY_im_j(1:3)=0.25*(
     &             DY_i_jp(1:3)+DY_i_jm(1:3)+DY_ip_j(1:3)+DY_ip_jm(1:3))
!     x+1/2,y



               if(e4coef(i2g) .lt. 50.0) then 
                 Tzero_fsh(i,jq,j,1,1)=
     &          (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &           +e4coef(i2g)
     &           *(1.0-1.0/sqrt(dot_product(DY_im_j(1:3),DY_im_j(1:3))))

                 Tzero_fsh(i,jq,j,2,2)=
     &           (1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))
     &           +e4coef(i2g)
     &           *(1.0-1.0/sqrt(dot_product(DX_i_jm(1:3),DX_i_jm(1:3))))        
                elseif(e4coef(i2g) .le. 100.0) then 
                  temp_ibm3= e4coef(i2g)-50.0 

                  Tzero_fsh(i,jq,j,1,1)=
     &            (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &             +temp_ibm3
     &           *( (sqrt(dot_product(DY_im_j(1:3),DY_im_j(1:3)))-1.0)        
     &            /(sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))) )

                  Tzero_fsh(i,jq,j,2,2)=
     &            (1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))
     &             +temp_ibm3
     &           *( (sqrt(dot_product(DX_i_jm(1:3),DX_i_jm(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))) )
                else
                  Tzero_fsh(i,jq,j,1,1)=0.0
                  Tzero_fsh(i,jq,j,2,2)=0.0
                endif

              endif

              DX_ip_j(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &         /dsf_IBM_fsh(i,jq,j) 
              DX_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))
     &         /dsf_IBM_fsh(i,jq-1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &         /ds_ibm_fsh(i,jq,j)

              if(e4coef(i2g) .lt. 0) then
              Tzero_fsh(i,jq,j,1,1)=
     &         dot_product(DX_ip_j(1:3),DX_ip_j(1:3))
              endif

              Bzero_fsh(i,jq,j,1,1)=
     &         dot_product(DXX_i_j(1:3),DXX_i_j(1:3))

 
              DY_i_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &         /dsf2_IBM_fsh(i,jq,j) 
              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &         /dsf2_IBM_fsh(i,jq,j-1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &         /ds2_ibm_fsh(i,jq,j)

              if(e4coef(i2g) .lt. 0) then
              Tzero_fsh(i,jq,j,2,2)=
     &         dot_product(DY_i_jp(1:3),DY_i_jp(1:3))
              endif

              Bzero_fsh(i,jq,j,2,2)=
     &         dot_product(DYY_i_j(1:3),DYY_i_j(1:3))


              DY_ip_j(1:3)=
     &        (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq+1,j))
     &         /dsf2_IBM_fsh(i,jq+1,j) 
              DY_im_j(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &         /dsf2_IBM_fsh(i,jq,j) 
              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &         /dk_ibm_fsh(i,jq,j)
              Bzero_fsh(i,jq,j,1,2)=
     &         dot_product(DXY_i_j(1:3),DXY_i_j(1:3))


              DX_i_jp(1:3)=
     &         (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq,j+1))
     &         /dsf_IBM_fsh(i,jq,j+1) 
              DX_i_jm(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &         /dsf_IBM_fsh(i,jq,j) 
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &         /dk2_ibm_fsh(i,jq,j)
              Bzero_fsh(i,jq,j,2,1)=
     &         dot_product(DYX_i_j(1:3),DYX_i_j(1:3))

              temp_ibm1=0.0
              do noj=1,3
                temp_ibm1=temp_ibm1+0.25d0*
     &         (DX_i_jp(noj)+DX_i_jm(noj))*(DY_ip_j(noj)+DY_im_j(noj))
              enddo

              Tzero_fsh(i,jq,j,2,1)=
     &         temp_ibm1
              Tzero_fsh(i,jq,j,1,2)=Tzero_fsh(i,jq,j,2,1)
           enddo
           enddo

           if(e4coef(i2g) .ge. 0) then
           do j=1,ns_ibm_r_fsh(i)
! boundary 1
           jq=1
              DX_im_j(1:3)=
     &         (   -GXtmp(1:3,jq+2,j)
     &         +4.0*GXtmp(1:3,jq+1,j)
     &         -3.0*GXtmp(1:3,jq,j)   )
     &          /(3.0*dsf_IBM_fsh(i,jq,j)-dsf_IBM_fsh(i,jq+1,j) )
!     x,y

           jq=1
              DX_ip_j(1:3)=
     &         (   -GXtmp(1:3,jq+2,j+1)
     &         +4.0*GXtmp(1:3,jq+1,j+1)
     &         -3.0*GXtmp(1:3,jq,j+1)   )
     &          /(3.0*dsf_IBM_fsh(i,jq,j+1)-dsf_IBM_fsh(i,jq+1,j+1) )
!     x,y+1

              DX_ip_j(1:3)=0.5*(DX_ip_j(1:3)+DX_im_j(1:3))

              DY_ip_j(1:3)=
     &             (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &             /dsf2_IBM_fsh(i,jq,j)
!     x,y+1/2

             if(e4coef(i2g) .lt. 50.0) then
             TzeroB_fsh(i,1,j,2)=
     &            (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &        +e4coef(i2g)
     &           *(1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))        
              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0
             TzeroB_fsh(i,1,j,2)=
     &           (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))) )
              else
             TzeroB_fsh(i,1,j,2)=0.0
              endif  
             TzeroB_fsh(i,1,j,1)=0.0
             TzeroB_fsh(i,1,j,3)=
     &              dot_product(DX_ip_j(1:3),DY_ip_j(1:3))

! boundary 2
           jq=nq_ibm_r_fsh(i)
              DX_im_j(1:3)=
     &         (    GXtmp(1:3,jq-2,j)
     &         -4.0*GXtmp(1:3,jq-1,j)
     &         +3.0*GXtmp(1:3,jq,j)   )
     &         /(3.0*dsf_IBM_fsh(i,jq-1,j)- dsf_IBM_fsh(i,jq-2,j) )

           jq=nq_ibm_r_fsh(i)
              DX_ip_j(1:3)=
     &         (    GXtmp(1:3,jq-2,j+1)
     &         -4.0*GXtmp(1:3,jq-1,j+1)
     &         +3.0*GXtmp(1:3,jq,j+1)   )
     &          /(3.0*dsf_IBM_fsh(i,jq-1,j+1)-dsf_IBM_fsh(i,jq-2,j+1) )
!     x,y+1

              DX_ip_j(1:3)=0.5*(DX_ip_j(1:3)+DX_im_j(1:3))

              DY_ip_j(1:3)=
     &             (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &         /dsf2_IBM_fsh(i,jq,j)
!     x,y+1/2


             if(e4coef(i2g) .lt. 50.0) then 
             TzeroB_fsh(i,2,j,2)=
     &              ( 
     &            (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &        +e4coef(i2g)
     &           *(1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))) 
     &              )
              elseif(e4coef(i2g) .le. 100.0) then   
              temp_ibm3= e4coef(i2g)-50.0
             TzeroB_fsh(i,2,j,2)=
     &              ( 
     &           (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))) )
     &              )
              else
             TzeroB_fsh(i,2,j,2)=0.0
              endif  
             TzeroB_fsh(i,2,j,1)=0.0
             TzeroB_fsh(i,1,j,3)=
     &              dot_product(DX_ip_j(1:3),DY_ip_j(1:3))
           enddo


           do jq=1,nq_ibm_r_fsh(i)
! boundary 3
           j=1
              DY_i_jm(1:3)=
     &         (   -GXtmp(1:3,jq,j+2)
     &         +4.0*GXtmp(1:3,jq,j+1)
     &         -3.0*GXtmp(1:3,jq,j)   )
     &          /(3.0*dsf2_IBM_fsh(i,jq,j)-dsf2_IBM_fsh(i,jq,j+1))
!     x+1/2,y

              DY_i_jp(1:3)=
     &         (   -GXtmp(1:3,jq+1,j+2)
     &         +4.0*GXtmp(1:3,jq+1,j+1)
     &         -3.0*GXtmp(1:3,jq+1,j)   )
     &             /(3.0*dsf2_IBM_fsh(i,jq+1,j)
     &         -dsf2_IBM_fsh(i,jq+1,j+1))
!     x+1/2,y


              DY_i_jp(1:3)=
     &           0.5*(DY_i_jm(1:3)+DY_i_jp(1:3))
              DX_ip_j(1:3)=
     &             (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &             /dsf_IBM_fsh(i,jq,j)
!     x,y+1/2


              TzeroB_fsh(i,3,jq,2)=0.0
              if(e4coef(i2g) .lt. 50.0) then
              TzeroB_fsh(i,3,jq,1)=
     &         (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &        +e4coef(i2g)
     &        *(1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))        
              elseif(e4coef(i2g) .le. 100.0) then   
              temp_ibm3= e4coef(i2g)-50.0
             TzeroB_fsh(i,3,jq,1)=
     &         (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))-1.0)        
     &            /(sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))) )
              else
             TzeroB_fsh(i,3,jq,1)=0.0
              endif  

              TzeroB_fsh(i,3,jq,3)=
     &              dot_product(DX_ip_j(1:3),DY_i_jp(1:3))
! boundary 4
           j=ns_ibm_r_fsh(i)
              DY_i_jm(1:3)=
     &         (    GXtmp(1:3,jq,j-2)
     &         -4.0*GXtmp(1:3,jq,j-1)
     &         +3.0*GXtmp(1:3,jq,j)   )
     &          /(3.0*dsf2_IBM_fsh(i,jq,j-1)-dsf2_IBM_fsh(i,jq,j-2))
!     x+1/2,y

              DY_i_jp(1:3)=
     &         (    GXtmp(1:3,jq+1,j-2)
     &         -4.0*GXtmp(1:3,jq+1,j-1)
     &         +3.0*GXtmp(1:3,jq+1,j)   )
     &             /(3.0*dsf2_IBM_fsh(i,jq+1,j-1)
     &             -dsf2_IBM_fsh(i,jq+1,j-2))
!     x+1/2,y

              DY_i_jp(1:3)=
     &           0.5*(DY_i_jm(1:3)+DY_i_jp(1:3))
              DX_ip_j(1:3)=
     &             (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &             /dsf_IBM_fsh(i,jq,j)
!     x,y+1/2


              TzeroB_fsh(i,4,jq,2)=0.0

              if(e4coef(i2g) .lt. 50.0) then 
              TzeroB_fsh(i,4,jq,1)=
     &              ( 
     &         (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &        +e4coef(i2g)
     &        *(1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))        
     &                  )      
              elseif(e4coef(i2g) .le. 100.0) then   
              temp_ibm3= e4coef(i2g)-50.0
             TzeroB_fsh(i,4,jq,1)=
     &              ( 
     &           (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))-1.0)        
     &            /(sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))) )
     &                 )
              else
             TzeroB_fsh(i,4,jq,1)=0.0
              endif  
              TzeroB_fsh(i,4,jq,3)=
     &              dot_product(DX_ip_j(1:3),DY_i_jp(1:3))
           enddo
           endif

           do jq=1,nq_ibm_r_fsh(i)
              jq_pm(2)=jq+1
              jq_pm(1)=jq-1
              if(jq.eq. nq_ibm_r_fsh(i)) then
                 if (BoundaryFibShell_fsh(i,2,2) .eq. 3) then
                    jqbndy=jq
                 else
                    jqbndy=1
                 endif 
                 jq_pm(2)=jqbndy
              elseif(jq.eq. 1) then
                 if (BoundaryFibShell_fsh(i,2,1) .eq. 3) then
                    jqbndy=nq_ibm_r_fsh(i)
                 else
                    jqbndy=jq
                 endif 
                 jq_pm(1)=jqbndy
              endif
           do j=1,ns_ibm_r_fsh(i)
              j_pm(2)=j+1
              j_pm(1)=j-1
              if(j.eq. ns_ibm_r_fsh(i)) then
                 if (BoundaryFibShell_fsh(i,1,2) .eq. 3) then
                    jbndy=j
                 else
                    jbndy=1
                 endif 
                 j_pm(2)=jbndy
              elseif(j.eq. 1) then
                 if (BoundaryFibShell_fsh(i,1,1) .eq. 3) then
                    jbndy=ns_ibm_r_fsh(i)
                 else
                    jbndy=j
                 endif
                 j_pm(1)=jbndy 
              endif

              do noi=1,nMem_Coef_fsh(i)
              kmem_i_j(noi,1)=
     &         (Mem_Coef_fsh(i,jq      ,j,noi)
     &         +Mem_Coef_fsh(i,jq_pm(2),j,noi))/2.0
              kmem_i_j(noi,2)=
     &         (Mem_Coef_fsh(i,jq,j      ,noi)
     &         +Mem_Coef_fsh(i,jq,j_pm(2),noi))/2.0
              kmem_i_j(noi,3)=
     &         (Mem_Coef_fsh(i,jq      ,j      ,noi)
     &         +Mem_Coef_fsh(i,jq_pm(2),j      ,noi)
     &         +Mem_Coef_fsh(i,jq      ,j_pm(2),noi)
     &         +Mem_Coef_fsh(i,jq_pm(2),j_pm(2),noi))/4.0
              enddo

              do noi=1,1
               kbend_i_j=
     &         (Ben_Coef_fsh(i,jq      ,j      ,3)
     &         +Ben_Coef_fsh(i,jq_pm(2),j      ,3)
     &         +Ben_Coef_fsh(i,jq      ,j_pm(2),3)
     &         +Ben_Coef_fsh(i,jq_pm(2),j_pm(2),3))/4.0
              enddo



              DX_ip_j(1:3)=
     &             (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &             /dsf_IBM_fsh(i,jq,j)
!     x+1/2,y

              if(e4coef(i2g) .ge. 0) then
              DX_im_j(1:3)=
     &                (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))
     &                /dsf_IBM_fsh(i,jq-1,j)
!     x-1/2,y

              DX_i_jp(1:3)=
     &         (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq,j+1))
     &             /dsf_IBM_fsh(i,jq,j+1)
!     x+1/2,y+1

              DX_im_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq-1,j+1))
     &             /dsf_IBM_fsh(i,jq-1,j+1)
!     x-1/2,y+1


              DX_i_jm(1:3)=0.25*(
     &             DX_ip_j(1:3)+DX_im_j(1:3)+DX_i_jp(1:3)+DX_im_jp(1:3))
!     x,y+1/2


              DY_i_jp(1:3)=
     &             (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &             /dsf2_IBM_fsh(i,jq,j)
!     x,y+1/2

              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &             /dsf2_IBM_fsh(i,jq,j-1)
!     x,y-1/2

              DY_ip_j(1:3)=
     &        (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq+1,j))
     &             /dsf2_IBM_fsh(i,jq+1,j)
!     x+1,y+1/2

              DY_ip_jm(1:3)=
     &        (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq+1,j-1))
     &             /dsf2_IBM_fsh(i,jq+1,j-1)
!     x+1,y-1/2 

              DY_im_j(1:3)=0.25*(
     &             DY_i_jp(1:3)+DY_i_jm(1:3)+DY_ip_j(1:3)+DY_ip_jm(1:3))
!     x+1/2,y

               if(e4coef(i2g) .lt. 50.0) then  
              Siga(jq,j,1,1)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &         (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &        +e4coef(i2g)
     &        *(1.0-1.0/sqrt(dot_product(DY_im_j(1:3),DY_im_j(1:3))))        
     &                  -Tzero_fsh(i,jq,j,1,1)  )
              elseif(e4coef(i2g) .le. 100.0) then 
                     temp_ibm3= e4coef(i2g)-50.0  
                Siga(jq,j,1,1)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DY_im_j(1:3),DY_im_j(1:3)))-1.0)        
     &            /(sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))) )
     &                  -Tzero_fsh(i,jq,j,1,1)  )
              else
                     temp_ibm1=(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))
     &                         -1.0)/2.0  
                     temp_ibm2=(dot_product(DY_im_j(1:3),DY_im_j(1:3))
     &                         -1.0)/2.0  
                     ksh_tmp=kmem_i_j(1,1)*temp_ibm1**2
     &                      +2.0*kmem_i_j(4,1)*temp_ibm1*temp_ibm2
     &                      +kmem_i_j(2,1)*temp_ibm2**2
                  Siga(jq,j,1,1)=kmem_i_j(5,1)*exp(ksh_tmp)
     &                             *( kmem_i_j(1,1)*temp_ibm1
     &                               +kmem_i_j(4,1)*temp_ibm2)
              endif

              else
              Siga(jq,j,1,1)=Mem_Coef_fsh(i,jq,j,1)
     &                       *(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))
     &                         -Tzero_fsh(i,jq,j,1,1))

              endif
 
              if(MaterialTypeIBM(i2g).eq. 1) then
              DX_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))/dsf_IBM_fsh(i,jq-1,j) 
              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /ds_ibm_fsh(i,jq,j)

              gama(jq,j,1,1)=Ben_Coef_fsh(i,jq,j,1)
     &                       *(dot_product(DXX_i_j(1:3),DXX_i_j(1:3))
     &                         -Bzero_fsh(i,jq,j,1,1))
              else
              gama(jq,j,1,1)=Ben_Coef_fsh(i,jq,j,1) 
              endif
              
              if(e4coef(i2g) .ge. 0.0) then
!     calculate for continuum fiber shell model (INCOMPLETE)
               if(e4coef(i2g) .lt. 50.0) then 
              Siga(jq,j,2,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &         (1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))
     &        +e4coef(i2g)
     &        *(1.0-1.0/sqrt(dot_product(DX_i_jm(1:3),DX_i_jm(1:3))))        
     &                  -Tzero_fsh(i,jq,j,2,2)   )
              elseif(e4coef(i2g) .le. 100.0) then 
                     temp_ibm3= e4coef(i2g)-50.0  
                Siga(jq,j,2,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DX_i_jm(1:3),DX_i_jm(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))) )
     &                  -Tzero_fsh(i,jq,j,1,1)  )
              else
                     temp_ibm1=(dot_product(DX_i_jm(1:3),DX_i_jm(1:3))
     &                         -1.0)/2.0  
                     temp_ibm2=(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))
     &                         -1.0)/2.0  
                  ksh_tmp=kmem_i_j(1,2)*temp_ibm1**2
     &                      +2.0*kmem_i_j(4,2)*temp_ibm1*temp_ibm2
     &                      +kmem_i_j(2,2)*temp_ibm2**2
               Siga(jq,j,2,2)=kmem_i_j(5,2)*exp(ksh_tmp)
     &                             *( kmem_i_j(1,2)*temp_ibm1
     &                               +kmem_i_j(4,2)*temp_ibm2)
               endif
              else
              DY_i_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))/dsf2_IBM_fsh(i,jq,j) 

              Siga(jq,j,2,2)=Mem_Coef_fsh(i,jq,j,2)
     &                       *(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))
     &                         -Tzero_fsh(i,jq,j,2,2))
              endif

! calculation for different bending energy density function 
              if(MaterialTypeIBM(i2g).eq. 1) then
              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &                     /dsf2_IBM_fsh(i,jq,j-1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /ds2_ibm_fsh(i,jq,j)

              gama(jq,j,2,2)=Ben_Coef_fsh(i,jq,j,2)
     &                       *(dot_product(DYY_i_j(1:3),DYY_i_j(1:3))
     &                         -Bzero_fsh(i,jq,j,2,2))
              else
              gama(jq,j,2,2)=Ben_Coef_fsh(i,jq,j,2)
              endif


              if(MaterialTypeIBM(i2g).eq. 1) then
              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &                     /dsf2_IBM_fsh(i,jq,j-1) 
              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /ds2_ibm_fsh(i,jq,j)

              gama(jq,j,2,2)=Ben_Coef_fsh(i,jq,j,2)
     &                       *(dot_product(DYY_i_j(1:3),DYY_i_j(1:3))
     &                         -Bzero_fsh(i,jq,j,2,2))
              else
              gama(jq,j,2,2)=Ben_Coef_fsh(i,jq,j,2)
              endif

              DY_ip_j(1:3)=
     &        (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq+1,j))
     &                     /dsf2_IBM_fsh(i,jq+1,j) 
              DY_im_j(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))/dsf2_IBM_fsh(i,jq,j) 

              DX_i_jp(1:3)=
     &         (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq,j+1))
     &                     /dsf_IBM_fsh(i,jq,j+1) 
              DX_i_jm(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &                     /dsf_IBM_fsh(i,jq,j) 

              if(MaterialTypeIBM(i2g).eq. 1) then
              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dk_ibm_fsh(i,jq,j)

              gama(jq,j,1,2)=Ben_Coef_fsh(i,jq,j,3)
     &                       *(dot_product(DXY_i_j(1:3),DXY_i_j(1:3))
     &                         -Bzero_fsh(i,jq,j,1,2))

              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dk2_ibm_fsh(i,jq,j)

              gama(jq,j,2,1)=Ben_Coef_fsh(i,jq,j,3)
     &                       *(dot_product(DYX_i_j(1:3),DYX_i_j(1:3))
     &                         -Bzero_fsh(i,jq,j,2,1))

              else
              gama(jq,j,1,2)=Ben_Coef_fsh(i,jq,j,3)
              gama(jq,j,2,1)=Ben_Coef_fsh(i,jq,j,3)
              endif

              temp_ibm1=0.0
              do noj=1,3
                temp_ibm1=temp_ibm1+0.25d0*
     &         (DX_i_jp(noj)+DX_i_jm(noj))*(DY_ip_j(noj)+DY_im_j(noj))
              enddo

              Siga(jq,j,2,1)=Mem_Coef_fsh(i,jq,j,3)*(
     &          temp_ibm1
     &         -Tzero_fsh(i,jq,j,2,1))

              Siga(jq,j,1,2)=Siga(jq,j,2,1)
           enddo !jq=1,nq_ibm_r_fsh(i)
           enddo !j =1,ns_ibm_r_fsh(i)

              if(Electromechanical(i2g)) then
              if(Piezo_Formulation_Flag(i2g) .eq. 0) then
              do jq=1,nq_ibm_r_fsh(i)
              do j=1,ns_ibm_r_fsh(i)
              DX_ip_j(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))/dsf_IBM_fsh(i,jq,j) 
              DX_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))/dsf_IBM_fsh(i,jq-1,j) 

              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /ds_ibm_fsh(i,jq,j)
 
              DY_i_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))/dsf2_IBM_fsh(i,jq,j) 

              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &                     /dsf2_IBM_fsh(i,jq,j-1) 

              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /ds2_ibm_fsh(i,jq,j)


              dXP1cont(1,1:3)=0.5*(DX_im_j(1:3)+DX_ip_j(1:3))
              dXP1cont(2,1:3)=0.5*(DY_i_jm(1:3)+DY_i_jp(1:3))
              call crossvec( 
     &             Normalcont1(1:3),dXP1cont(1,1:3),dXP1cont(2,1:3))

              do k=1,3
                ndot_ibm_fsh0(i,jq,j)=
     &             piezo_Coef_fsh(i,jq,j,3)*DXX_i_j(k)*Normalcont1(k)
     &            +piezo_Coef_fsh(i,jq,j,4)*DYY_i_j(k)*Normalcont1(k)
              enddo
              

              DY_ip_j(1:3)=0.25*( 
     &           GXtmp(1:3,jq,j)
     &          +GXtmp(1:3,jq,j+1)
     &          +GXtmp(1:3,jq+1,j)
     &          +GXtmp(1:3,jq+1,j+1))            
              DY_im_j(1:3)=0.25*( 
     &           GXtmp(1:3,jq,j) 
     &          +GXtmp(1:3,jq,j+1) 
     &          +GXtmp(1:3,jq-1,j) 
     &          +GXtmp(1:3,jq-1,j+1))

              DX_i_jp(1:3)=0.25*( 
     &           GXtmp(1:3,jq,j) 
     &          +GXtmp(1:3,jq,j-1) 
     &          +GXtmp(1:3,jq+1,j) 
     &          +GXtmp(1:3,jq+1,j-1))     
              DX_i_jm(1:3)=0.25*( 
     &           GXtmp(1:3,jq,j) 
     &          +GXtmp(1:3,jq,j-1) 
     &          +GXtmp(1:3,jq-1,j) 
     &          +GXtmp(1:3,jq-1,j-1))



              DXY_i_j(1:3)=(DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dk_ibm_fsh(i,jq,j)
              DYX_i_j(1:3)=(DX_i_jp(1:3)-DX_i_jm(1:3))
     &                     /dk_ibm_fsh(i,jq,j-1)
              DYY_i_j(1:3)=(DXY_i_j(1:3)-DYX_i_j(1:3)) 
     &                      /ds2_ibm_fsh(i,jq,j)

              DXY_i_j(1:3)=(DY_ip_j(1:3)-DX_i_jp(1:3))
     &                     /dk2_ibm_fsh(i,jq,j)
              DYX_i_j(1:3)=(DY_im_j(1:3)-DX_i_jm(1:3))
     &                     /dk2_ibm_fsh(i,jq-1,j)
              DXX_i_j(1:3)=(DXY_i_j(1:3)-DYX_i_j(1:3)) 
     &                      /ds_ibm_fsh(i,jq,j)

              do k=1,3
                ndot_ibm_fsh0(i,jq,j)=ndot_ibm_fsh0(i,jq,j)
     &            +piezo_Coef_fsh(i,jq,j,5)
     &                *(DXX_i_j(k)+DYY_i_j(k))*Normalcont1(k)
              enddo

              enddo
              enddo

              elseif(Piezo_Formulation_Flag(i2g) .eq. 1) then
              do jq=1,nq_ibm_r_fsh(i)
              do j=1,ns_ibm_r_fsh(i)

              DX_ip_j(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))/dsf_IBM_fsh(i,jq,j) 
              DX_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))/dsf_IBM_fsh(i,jq-1,j) 

              DXX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /ds_ibm_fsh(i,jq,j)
 
              DY_i_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))/dsf2_IBM_fsh(i,jq,j) 

              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &                     /dsf2_IBM_fsh(i,jq,j-1) 

              DYY_i_j(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /ds2_ibm_fsh(i,jq,j)


              dXP1cont(1,1:3)=0.5*(DX_im_j(1:3)+DX_ip_j(1:3))
              dXP1cont(2,1:3)=0.5*(DY_i_jm(1:3)+DY_i_jp(1:3))
              call crossvec( 
     &             Normalcont1(1:3),dXP1cont(1,1:3),dXP1cont(2,1:3))

!central canculation
              DY_i_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &                     /dsf2_IBM_fsh(i,jq,j) 
              DY_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &                     /dsf2_IBM_fsh(i,jq,j-1) 
              X22piezo(1:3)=( DY_i_jp(1:3)-DY_i_jm(1:3))
     &                     /ds2_ibm_fsh(i,jq,j) 
              X2piezo(1:3)=( DY_i_jp(1:3)+DY_i_jm(1:3))
     &                     /2.0

              DX_ip_j(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &                     /dsf_IBM_fsh(i,jq,j) 
              DX_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))

     &                     /dsf_IBM_fsh(i,jq-1,j) 
              X11piezo(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /ds_ibm_fsh(i,jq,j) 
              X1piezo(1:3)=( DX_ip_j(1:3)+DX_im_j(1:3))
     &                     /2.0



              DY_ip_j(1:3)=
     &        (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq+1,j))
     &                     /dsf2_IBM_fsh(i,jq+1,j) 
              DY_im_j(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &                     /dsf2_IBM_fsh(i,jq,j) 

              DX_i_jp(1:3)=
     &         (GXtmp(1:3,jq+1,j+1)-GXtmp(1:3,jq,j+1))
     &                     /dsf_IBM_fsh(i,jq,j+1) 
              DX_i_jm(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &                     /dsf_IBM_fsh(i,jq,j) 


              DXY_i_j(1:3)=( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dk_ibm_fsh(i,jq,j)
              DYX_i_j(1:3)=( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dk2_ibm_fsh(i,jq,j)


              DY_ip_j(1:3)=
     &        (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &                     /dsf2_IBM_fsh(i,jq,j) 
              DY_im_j(1:3)=
     &         (GXtmp(1:3,jq-1,j+1)-GXtmp(1:3,jq-1,j))
     &                     /dsf2_IBM_fsh(i,jq-1,j) 

              DX_i_jp(1:3)=
     &         (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq-1,j+1))
     &                     /dsf_IBM_fsh(i,jq-1,j+1) 
              DX_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))
     &                     /dsf_IBM_fsh(i,jq-1,j) 


              DXY_i_j(1:3)=
     &                     DXY_i_j(1:3)
     &                     +( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dk_ibm_fsh(i,jq-1,j)
              DYX_i_j(1:3)=
     &                     DYX_i_j(1:3)
     &                     +( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dk2_ibm_fsh(i,jq-1,j)


              DY_ip_j(1:3)=
     &        (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq+1,j-1))
     &                     /dsf2_IBM_fsh(i,jq+1,j-1) 
              DY_im_j(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &                     /dsf2_IBM_fsh(i,jq,j-1) 

              DX_i_jp(1:3)=
     &         (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &                     /dsf_IBM_fsh(i,jq,j) 
              DX_i_jm(1:3)=
     &         (GXtmp(1:3,jq+1,j-1)-GXtmp(1:3,jq,j-1))
     &                     /dsf_IBM_fsh(i,jq,j-1) 


              DXY_i_j(1:3)=
     &                     DXY_i_j(1:3)
     &                     +( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dk_ibm_fsh(i,jq,j-1)
              DYX_i_j(1:3)=
     &                     DYX_i_j(1:3)
     &                     +( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dk2_ibm_fsh(i,jq,j-1)




              DY_ip_j(1:3)=
     &        (GXtmp(1:3,jq,j)-GXtmp(1:3,jq,j-1))
     &                     /dsf2_IBM_fsh(i,jq,j-1) 
              DY_im_j(1:3)=
     &         (GXtmp(1:3,jq-1,j)-GXtmp(1:3,jq-1,j-1))
     &                     /dsf2_IBM_fsh(i,jq-1,j-1) 

              DX_i_jp(1:3)=
     &         (GXtmp(1:3,jq,j)-GXtmp(1:3,jq-1,j))
     &                     /dsf_IBM_fsh(i,jq-1,j) 
              DX_i_jm(1:3)=
     &         (GXtmp(1:3,jq,j-1)-GXtmp(1:3,jq-1,j-1))
     &                     /dsf_IBM_fsh(i,jq-1,j-1) 


              DXY_i_j(1:3)=
     &                     DXY_i_j(1:3)
     &                     +( DY_ip_j(1:3)-DY_im_j(1:3))
     &                     /dk_ibm_fsh(i,jq-1,j-1)
              DYX_i_j(1:3)=
     &                     DYX_i_j(1:3)
     &                     +( DX_ip_j(1:3)-DX_im_j(1:3))
     &                     /dk2_ibm_fsh(i,jq-1,j-1)

              X12piezo(1:3)=0.25*( (DXY_i_j(1:3)+DYX_i_j(1:3))/2.0)


              normalpiezo(1:3)=Normalcont1(1:3)

              Epiezo=dot_product(X1piezo(1:3),X1piezo(1:3))
!     E=X_1 . X_1
              Fpiezo=dot_product(X1piezo(1:3),X2piezo(1:3))
!     F=X_1 . X_2
              Gpiezo=dot_product(X2piezo(1:3),X2piezo(1:3))
!     G=X_2 . X_2

              Lpiezo=dot_product(X11piezo(1:3),normalpiezo(1:3))
!     L=X_11 . n
              Mpiezo=dot_product(X12piezo(1:3),normalpiezo(1:3))
!     M=X_12 . n
              Npiezo=dot_product(X22piezo(1:3),normalpiezo(1:3))
!     N=X_22 . n
           
              Hnpiezo=Epiezo*Npiezo-2.0*Fpiezo*Mpiezo+Gpiezo*Lpiezo
!     Hn=EN-2FM+GL
              Hdpiezo=2.0*(Epiezo*Gpiezo-Fpiezo**2)
!     Hd=2(EG-F^2)

              ndot_ibm_fsh0(i,jq,j)=Hnpiezo/Hdpiezo
              enddo
              enddo
              else  !Piezo_Formulation_Flag(i2g)
                 write(*,*) 'wrong Piezo_Formulation_Flag for ',i2g
                 stop
              endif !Piezo_Formulation_Flag(i2g)
               do jq=1,nq_ibm_r_fsh(i)
               do j=1,ns_ibm_r_fsh(i)
                  ndot_ibm_fsh(i,jq,j)   =ndot_ibm_fsh0(i,jq,j)
                  ndotold_ibm_fsh(i,jq,j)=ndot_ibm_fsh0(i,jq,j)
                  ndotpre_ibm_fsh(i,jq,j)=ndot_ibm_fsh0(i,jq,j)
               enddo
               enddo
              endif


           if(e4coef(i2g) .ge. 0.0)then
           do j=1,ns_ibm_r_fsh(i)
! boundary 1
           jq=1
              DX_im_j(1:3)=
     &         (   -GXtmp(1:3,jq+2,j)
     &         +4.0*GXtmp(1:3,jq+1,j)
     &         -3.0*GXtmp(1:3,jq,j)   )
     &          /(3.0*dsf_IBM_fsh(i,jq,j)-dsf_IBM_fsh(i,jq+1,j) )
!     x,y

           jq=1
              DX_ip_j(1:3)=
     &         (   -GXtmp(1:3,jq+2,j+1)
     &         +4.0*GXtmp(1:3,jq+1,j+1)
     &         -3.0*GXtmp(1:3,jq,j+1)   )
     &          /(3.0*dsf_IBM_fsh(i,jq,j+1)-dsf_IBM_fsh(i,jq+1,j+1) )
!     x,y+1

              DX_ip_j(1:3)=0.5*(DX_ip_j(1:3)+DX_im_j(1:3))

              DY_ip_j(1:3)=
     &             (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &             /dsf2_IBM_fsh(i,jq,j)
!     x,y+1/2

              if(e4coef(i2g) .lt. 50.0) then   
              sigB(1,j,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &            (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &        +e4coef(i2g)
     &           *(1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))        
     &                  -TzeroB_fsh(i,1,j,2)  )

              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0          
                     sigB(1,j,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))) )
     &                  -TzeroB_fsh(i,1,j,2)  )

              endif
              sigB(1,j,1)=0.0
              sigB(1,j,3)=Mem_Coef_fsh(i,jq,j,3)
     &              *(dot_product(DX_ip_j(1:3),DY_ip_j(1:3))
     &                  -TzeroB_fsh(i,1,j,3)  )
! boundary 2
           jq=nq_ibm_r_fsh(i)
              DX_im_j(1:3)=
     &         (    GXtmp(1:3,jq-2,j)
     &         -4.0*GXtmp(1:3,jq-1,j)
     &         +3.0*GXtmp(1:3,jq,j)   )
     &         /(3.0*dsf_IBM_fsh(i,jq-1,j)- dsf_IBM_fsh(i,jq-2,j) )

           jq=nq_ibm_r_fsh(i)
              DX_ip_j(1:3)=
     &         (    GXtmp(1:3,jq-2,j+1)
     &         -4.0*GXtmp(1:3,jq-1,j+1)
     &         +3.0*GXtmp(1:3,jq,j+1)   )
     &          /(3.0*dsf_IBM_fsh(i,jq-1,j+1)-dsf_IBM_fsh(i,jq-2,j+1))
!     x,y+1

              DX_ip_j(1:3)=0.5*(DX_ip_j(1:3)+DX_im_j(1:3))

              DY_ip_j(1:3)=
     &             (GXtmp(1:3,jq,j+1)-GXtmp(1:3,jq,j))
     &             /dsf2_IBM_fsh(i,jq,j)
!     x,y+1/2

              sigB(2,j,1)=0.0
              sigB(2,j,3)=Mem_Coef_fsh(i,jq,j,3)
     &              *(dot_product(DX_ip_j(1:3),DY_ip_j(1:3))
     &                  -TzeroB_fsh(i,2,j,3)  )
              if(e4coef(i2g) .lt. 50.0) then

              sigB(2,j,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &            (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &        +e4coef(i2g)
     &           *(1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))) 
     &                  -TzeroB_fsh(i,2,j,2)  )

              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0          
                     sigB(2,j,2)=Mem_Coef_fsh(i,jq,j,2)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))-1.0)        
     &            /(sqrt(dot_product(DY_ip_j(1:3),DY_ip_j(1:3)))) )
     &                  -TzeroB_fsh(i,2,j,2)  )
              endif
           enddo


           do jq=1,nq_ibm_r_fsh(i)
! boundary 3
           j=1
              DY_i_jm(1:3)=
     &         (   -GXtmp(1:3,jq,j+2)
     &         +4.0*GXtmp(1:3,jq,j+1)
     &         -3.0*GXtmp(1:3,jq,j)   )
     &          /(3.0*dsf2_IBM_fsh(i,jq,j)-dsf2_IBM_fsh(i,jq,j+1))
!     x+1/2,y

              DY_i_jp(1:3)=
     &         (   -GXtmp(1:3,jq+1,j+2)
     &         +4.0*GXtmp(1:3,jq+1,j+1)
     &         -3.0*GXtmp(1:3,jq+1,j)   )
     &             /(3.0*dsf2_IBM_fsh(i,jq+1,j)
     &             -dsf2_IBM_fsh(i,jq+1,j+1))
!     x+1/2,y


              DY_i_jp(1:3)=
     &           0.5*(DY_i_jm(1:3)+DY_i_jp(1:3))
              DX_ip_j(1:3)=
     &             (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &             /dsf_IBM_fsh(i,jq,j)
!     x,y+1/2



              sigB(3,jq,2)=0.0
              if(e4coef(i2g) .lt. 50.0) then 
              sigB(3,jq,1)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &         (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &        +e4coef(i2g)
     &        *(1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))        
     &                  -TzeroB_fsh(i,3,jq,1)  )
              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0          
                     sigB(3,jq,1)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))-1.0)        
     &            /(sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))) )
     &                  -TzeroB_fsh(i,3,jq,1)  )
              endif
              sigB(3,jq,3)=Mem_Coef_fsh(i,jq,j,3)
     &              *(dot_product(DX_ip_j(1:3),DY_i_jp(1:3))
     &                  -TzeroB_fsh(i,3,jq,3)  )
! boundary 4
           j=ns_ibm_r_fsh(i)
              DY_i_jm(1:3)=
     &         (    GXtmp(1:3,jq,j-2)
     &         -4.0*GXtmp(1:3,jq,j-1)
     &         +3.0*GXtmp(1:3,jq,j)   )
     &          /(3.0*dsf2_IBM_fsh(i,jq,j-1)-dsf2_IBM_fsh(i,jq,j-2))
!     x+1/2,y

              DY_i_jp(1:3)=
     &         (    GXtmp(1:3,jq+1,j-2)
     &         -4.0*GXtmp(1:3,jq+1,j-1)
     &         +3.0*GXtmp(1:3,jq+1,j)   )
     &             /(3.0*dsf2_IBM_fsh(i,jq+1,j-1)
     &             -dsf2_IBM_fsh(i,jq+1,j-2))
!     x+1/2,y

              DY_i_jp(1:3)=
     &           0.5*(DY_i_jm(1:3)+DY_i_jp(1:3))
              DX_ip_j(1:3)=
     &             (GXtmp(1:3,jq+1,j)-GXtmp(1:3,jq,j))
     &             /dsf_IBM_fsh(i,jq,j)
!     x,y+1/2


              sigB(4,jq,2)=0.0
              sigB(4,jq,3)=Mem_Coef_fsh(i,jq,j,3)
     &              *(dot_product(DX_ip_j(1:3),DY_i_jp(1:3))
     &                  -TzeroB_fsh(i,4,jq,3)  )

              if(e4coef(i2g) .lt. 50.0) then 
              sigB(4,jq,1)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &         (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &        +e4coef(i2g)
     &        *(1.0-1.0/sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3))))        
     &                  -TzeroB_fsh(i,4,jq,1)  )
              elseif(e4coef(i2g) .le. 100.0) then   
                     temp_ibm3= e4coef(i2g)-50.0          
                     sigB(4,jq,1)=Mem_Coef_fsh(i,jq,j,1)
     &              *( 
     &           (1.0-1.0/sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3))))
     &             +temp_ibm3
     &          *( (sqrt(dot_product(DY_i_jp(1:3),DY_i_jp(1:3)))-1.0)        
     &            /(sqrt(dot_product(DX_ip_j(1:3),DX_ip_j(1:3)))) )
     &                  -TzeroB_fsh(i,4,jq,1)  )
              endif
           enddo
        endif

! if initial bending is not zero calculate initial bending 
      if(BendingZero(i2g) .ne. 0) then
      if(.not. allocated(aacsrfibershell))  then   
        allocate (iacsrfibershell(ns_ibm_fsh*nq_ibm_fsh+1))   
        allocate(vecRHSfibershell(ns_ibm_fsh*nq_ibm_fsh),
     &           vecXfibershell(ns_ibm_fsh*nq_ibm_fsh) )
        iacsrfibershell=0
        
        ntermiMAXfibershell=ns_ibm_fsh*nq_ibm_fsh*25
        allocate(aacsrfibershell(ntermiMAXfibershell)
     &          ,jacsrfibershell(ntermiMAXfibershell))
        jacsrfibershell=0;aacsrfibershell=0.d0
      endif
      if(.not. allocated(Stiff_link_fsh)) 
     &     allocate (Stiff_link_fsh(ns_ibm_fsh,nq_ibm_fsh)) 
         Stiff_link_fsh=0.0
         FlagForceCalculation=1
            if(e4coef(i2g) .ge. 0 ) then
       Call  matfibrous_Yue(
     &         ns_ibm_fsh*nq_ibm_fsh,ntermiMAXfibershell
     &        ,iacsrfibershell,jacsrfibershell,aacsrfibershell
     &        ,nq_ibm_r_fsh(i),ns_ibm_r_fsh(i)
     &        ,1.0d0
     &        ,density_coef(i2g),cs_ibm(i2g),density_coefp(i2g)
     &        ,BoundaryFibShell_fsh(i,1:2,1:2)
     &        ,Stiff_link_fsh(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))
     &        ,Siga(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),1,1)   
     &        ,Siga(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),1,2) 
     &        ,Siga(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),2,1) 
     &        ,Siga(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),2,2)            
     &        ,gama(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),1,1)   
     &        ,gama(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),1,2) 
     &        ,gama(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),2,1) 
     &        ,gama(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),2,2)
     &        ,ds_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))   
     &        ,dsf_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i)) 
     &        ,dk_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i)) 
     &        ,dkf_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))    
     &        ,ds2_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))   
     &        ,dsf2_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i)) 
     &        ,dk2_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i)) 
     &        ,dkf2_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))    
     &        ,sigB(1:4,1:nq_ibm_r_fsh(i)+ns_ibm_r_fsh(i),1)   
     &        ,sigB(1:4,1:nq_ibm_r_fsh(i)+ns_ibm_r_fsh(i),2) 
     &        ,sigB(1:4,1:nq_ibm_r_fsh(i)+ns_ibm_r_fsh(i),3) 
     &        ,e4coef(i2g)            
     &        ,0
     &        ,FlagForceCalculation(1:nq_ibm_r_fsh(i)
     &                             ,1:ns_ibm_r_fsh(i) ) )
            else
       Call  matfibrous_Huang(
     &         ns_ibm_fsh*nq_ibm_fsh,ntermiMAXfibershell
     &        ,iacsrfibershell,jacsrfibershell,aacsrfibershell
     &        ,nq_ibm_r_fsh(i),ns_ibm_r_fsh(i)
     &        ,1.0d0
     &        ,density_coef(i2g),cs_ibm(i2g),density_coefp(i2g)
     &        ,BoundaryFibShell_fsh(i,1:2,1:2)
     &        ,Stiff_link_fsh(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))
     &        ,Siga(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),1,1)   
     &        ,Siga(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),1,2) 
     &        ,Siga(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),2,1) 
     &        ,Siga(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),2,2)            
     &        ,gama(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),1,1)   
     &        ,gama(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),1,2) 
     &        ,gama(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),2,1) 
     &        ,gama(1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i),2,2)
     &        ,ds_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))   
     &        ,dsf_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i)) 
     &        ,dk_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i)) 
     &        ,dkf_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))    
     &        ,ds2_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))   
     &        ,dsf2_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i)) 
     &        ,dk2_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i)) 
     &        ,dkf2_ibm_fsh(i,1:nq_ibm_r_fsh(i),1:ns_ibm_r_fsh(i))    
     &        ,0)
            endif
!Solve x direction
       j2=0
       do jq=1,nq_ibm_r_fsh(i) 
       do j=1,ns_ibm_r_fsh(i) 
            j2=j2+1
            vecXfibershell(j2) = gx_ibm_massive_fsh(i,jq,j)
       end do
       end do
       call ax_cr ( ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i)
     & ,ntermiMAXfibershell
     & ,iacsrfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i)+1)
     & , jacsrfibershell(1: ntermiMAXfibershell)
     & , aacsrfibershell(1: ntermiMAXfibershell)
     & , vecXfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i))
     & , vecRHSfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i))
     &   )
       j2=0
       do jq=1,nq_ibm_r_fsh(i) 
       do j=1,ns_ibm_r_fsh(i) 
            j2=j2+1
            FB_1_IBM_fsh(i,jq,j)=vecRHSfibershell(j2)
            if(FlagForceCalculation(jq,j) .eq. 0) then
              FB_1_IBM_fsh(i,jq,j)=0.0
            endif
       end do
       end do
       vecXfibershell=0.0
       vecRHSfibershell=0.0
!Solve y direction
       j2=0
       do jq=1,nq_ibm_r_fsh(i) 
       do j=1,ns_ibm_r_fsh(i) 
            j2=j2+1
            vecXfibershell(j2) = gy_ibm_massive_fsh(i,jq,j)
       end do
       end do
       call ax_cr ( ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i)
     & ,ntermiMAXfibershell
     & ,iacsrfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i)+1)
     & , jacsrfibershell(1: ntermiMAXfibershell)
     & , aacsrfibershell(1: ntermiMAXfibershell)
     & , vecXfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i))
     & , vecRHSfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i))
     &   )
       j2=0
       do jq=1,nq_ibm_r_fsh(i) 
       do j=1,ns_ibm_r_fsh(i) 
            j2=j2+1
            FB_2_IBM_fsh(i,jq,j)=vecRHSfibershell(j2)
            if(FlagForceCalculation(jq,j) .eq. 0) then
              FB_2_IBM_fsh(i,jq,j)=0.0
            endif
       end do
       end do
       vecXfibershell=0.0
       vecRHSfibershell=0.0
!Solve z direction
        j2=0
        do jq=1,nq_ibm_r_fsh(i) 
        do j=1,ns_ibm_r_fsh(i) 
            j2=j2+1
            vecXfibershell(j2) = gz_ibm_massive_fsh(i,jq,j)
        end do
        end do
        call ax_cr ( ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i)
     &  ,ntermiMAXfibershell
     &  ,iacsrfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i)+1)
     &  , jacsrfibershell(1: ntermiMAXfibershell)
     &  , aacsrfibershell(1: ntermiMAXfibershell)
     &  , vecXfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i))
     &  , vecRHSfibershell(1:ns_ibm_r_fsh(i)*nq_ibm_r_fsh(i))
     &   )
        j2=0
        do jq=1,nq_ibm_r_fsh(i) 
        do j=1,ns_ibm_r_fsh(i) 
            j2=j2+1
            FB_3_IBM_fsh(i,jq,j)=vecRHSfibershell(j2)
            if(FlagForceCalculation(jq,j) .eq. 0) then
              FB_3_IBM_fsh(i,jq,j)=0.0
            endif

        end do
        end do
        vecXfibershell=0.0
        vecRHSfibershell=0.0

        else
         do jq=1,nq_ibm_r_fsh(i) 
         do j=1,ns_ibm_r_fsh(i) 
            FB_1_IBM_fsh(i,jq,j)=0.0
            FB_2_IBM_fsh(i,jq,j)=0.0
            FB_3_IBM_fsh(i,jq,j)=0.0
         end do
         end do
        endif

        do ibkt2=1,Ncontactplanes(i2g) 
                contactplane_point(1:3)= 
     &                       contactplanesPara(i2g,ibkt2,1:3)
                contactplane_n(1:3)=
     &                       contactplanesPara(i2g,ibkt2,4:6)

                dncont=0.0
                do j=1,3
                   dncont=dncont+contactplane_n(j) 
     &                          *contactplane_point(j)
                enddo
                contactplane_n(4)=-dncont
!     plane defined as ax+by+cz+d=0, here we save d

                do jq=1,nq_ibm_r_fsh(ibdy) 
                do j=1,ns_ibm_r_fsh(ibdy) 
                  ksh_tmp=contactplane_n(1) 
     &                   *GXtmp(1,jq,j) 
     &                   +contactplane_n(2) 
     &                   *GXtmp(2,jq,j)
     &                   +contactplane_n(3) 
     &                   *GXtmp(3,jq,j)
     &                   +contactplane_n(4) 

                  if(ksh_tmp .ge. impactpen(3)*ddcont) then
                     contact_coef_fsh(ibdy,ibkt2,jq,j)=impactpen(3) 
     &                                                *ddcont
                  elseif(ksh_tmp .ge. 0.05*impactpen(3)*ddcont) then
                     contact_coef_fsh(ibdy,ibkt2,jq,j)=ksh_tmp
                  else
                     contact_coef_fsh(ibdy,ibkt2,jq,j)= 
     &                0.05*impactpen(3)*ddcont
                  endif
                enddo 
                enddo
         enddo
      enddo ! nr_ibm_fsh
#endif
      
      endif  !nRead

! passing to outputs

      do i=1,nr_ibm
       nIBM_r(i) =ns_IBM_r(i)
       nIBM_rq(i)=nq_IBM_r(i)
      enddo

      do i=1,nr_ibm_fib
       nIBM_r_fib(i) =ns_IBM_r_fib(i)
      enddo

      do i=1,nr_ibm_fsh
       nIBM_r_fsh(i) =ns_IBM_r_fsh(i)
       nIBM_rq_fsh(i)=nq_IBM_r_fsh(i)
      enddo

      do i=1,nr_ibm_esh
       nIBM_r_esh(i) =ns_IBM_r_esh(i)
       do j=1,ns_ibm_r_esh(i)
          das=0.d0
          do i2=1,nposele_esh(i,j,1)
            das=das+area(i,nposele_esh(i,j,1+i2))
          end do
          ds_ibm_esh_hold=das/dble(nposele_esh(i,j,1))
          ds_ibm_esh(i,j)=ds_ibm_esh_hold
       end do
      enddo

      do i=1,nr_ibm_fbc
       nIBM_r_fbc(i) =ns_IBM_r_fbc(i)
       do j=1,ns_ibm_r_fbc(i)
          ds_ibm_fbc(i,j)=areapnt_fbc(i,j)
       end do
      enddo

      dtypeDelta(1)=Delta_typeX
      dtypeDelta(2)=Delta_typeY
      dtypeDelta(3)=Delta_typeZ

      call  BucketSorting
      timebucket=0.0
      return
      end 

      subroutine copy_nIBM( 
     &  nIBM_rq,nIBM_rq_fsh
     & ,nIBM_r,nIBM_r_fib,nIBM_r_fsh
     & ,nIBM_r_esh,nIBM_r_fbc
     & ,Ns_IBM_fib_out
     & ,Ns_IBM_fsh_out
     & ,Nq_IBM_fsh_out
     & ,Ns_IBM_esh_out
     & ,Ns_IBM_fbc_out)
      USE HeaderFSI
      implicit none

      integer nIBM_r(nr_IBM)
      integer nIBM_r_fib(nr_IBM_fib)
      integer nIBM_r_fsh(nr_IBM_fsh)
      integer nIBM_r_esh(nr_IBM_esh)
      integer nIBM_r_fbc(nr_IBM_fbc)
      integer nIBM_rq(nr_IBM),nIBM_rq_fsh(nr_IBM_fsh)
      integer Ns_IBM_fib_out
      integer Ns_IBM_fsh_out
      integer Nq_IBM_fsh_out
      integer Ns_IBM_esh_out
      integer Ns_IBM_fbc_out

      integer i

      Ns_IBM_fib_out = Ns_IBM_fib
      Ns_IBM_fsh_out = Ns_IBM_fsh
      Nq_IBM_fsh_out = Nq_IBM_fsh
      Ns_IBM_esh_out = Ns_IBM_esh
      Ns_IBM_fbc_out = Ns_IBM_fbc

      do i=1,nr_ibm
       nIBM_r(i) =ns_IBM_r(i)
       nIBM_rq(i)=nq_IBM_r(i)
      enddo
      do i=1,nr_ibm_fib
       nIBM_r_fib(i) =ns_IBM_r_fib(i)
      enddo
      do i=1,nr_ibm_fsh
       nIBM_r_fsh(i) =ns_IBM_r_fsh(i)
       nIBM_rq_fsh(i)=nq_IBM_r_fsh(i)
      enddo
      do i=1,nr_ibm_esh
       nIBM_r_esh(i) =ns_IBM_r_esh(i)
      enddo
      do i=1,nr_ibm_fbc
       nIBM_r_fbc(i) =ns_IBM_r_fbc(i)
      enddo

      return
      end 

      subroutine copy_nmaxIBM( 
     &  nr_IBM_out
     & ,CTML_num_scalars_out,
     & ,nr_IBM_fib_out
     & ,nr_IBM_fsh_out
     & ,nr_IBM_esh_out
     & ,nr_IBM_fbc_out
     & ,Ns_IBM_fib_out
     & ,Ns_IBM_fsh_out
     & ,Nq_IBM_fsh_out
     & ,Ns_IBM_esh_out
     & ,Ns_IBM_fbc_out)
      USE HeaderFSI
      implicit none

      integer Nr_IBM_out
      integer CTML_num_scalars_out
      integer Nr_IBM_fib_out
      integer Nr_IBM_fsh_out
      integer Nr_IBM_esh_out
      integer Nr_IBM_fbc_out
      integer Ns_IBM_fib_out
      integer Ns_IBM_fsh_out
      integer Nq_IBM_fsh_out
      integer Ns_IBM_esh_out
      integer Ns_IBM_fbc_out

      Nr_IBM_out = Nr_IBM
      CTML_num_scalars_out=0
      if (Nr_IBM_fsh.eq.0) then
! do nothing
      else if (Nr_IBM_fsh.gt.0) then
! ndot, vPiezo
       CTML_num_scalars_out=2
      else
       print *,"Nr_IBM_fsh invalid"
       stop
      endif

      Nr_IBM_fib_out = Nr_IBM_fib
      Nr_IBM_fsh_out = Nr_IBM_fsh
      Nr_IBM_esh_out = Nr_IBM_esh
      Nr_IBM_fbc_out = Nr_IBM_fbc

      Ns_IBM_fib_out = Ns_IBM_fib
      Ns_IBM_fsh_out = Ns_IBM_fsh
      Nq_IBM_fsh_out = Nq_IBM_fsh
      Ns_IBM_esh_out = Ns_IBM_esh
      Ns_IBM_fbc_out = Ns_IBM_fbc

      return
      end 

      subroutine copy_ibm_fib(
     &   ilo,ihi,jlo,jhi,klo,khi
     &  ,ds_fib
     &  ,gxyz_ibm_massiveout_fib)
      USE HeaderFSI
      implicit none

      integer ilo,ihi,jlo,jhi,klo,khi
      real*8 GXYZ_IBM_MASSIVEout_fib(Nr_IBM_fib
     &   ,ilo:ihi,jlo:jhi,klo:khi,3)
      real*8 ds_fib(Nr_IBM_fib
     &   ,ilo:ihi,jlo:jhi,klo:khi)

      integer i,j

      do i=1,nr_ibm_fib

       if ((ilo.le.1-Ns_IBMB_fib).and. 
     &     (ilo.le.0).and.
     &     (ihi.ge.ns_ibm_r_fib(i)+1).and.
     &     (ihi.ge.ns_ibm_r_fib(i)+Ns_IBMB_fib).and.
     &     (jlo.eq.1).and.
     &     (jhi.eq.1).and. 
     &     (klo.eq.1).and. 
     &     (khi.eq.1)) then
!       do nothing
       else
        print *,"ilo,ihi,jlo,jhi,klo,khi invalid"
        stop
       endif

       do j=0,ns_ibm_r_fib(i)+1
        gxyz_ibm_massiveout_fib(i,j,jlo,klo,1)=gx_ibm_massive_fib(i,j)
        gxyz_ibm_massiveout_fib(i,j,jlo,klo,2)=gy_ibm_massive_fib(i,j)
        gxyz_ibm_massiveout_fib(i,j,jlo,klo,3)=gz_ibm_massive_fib(i,j)
       end do
       do j=1,ns_ibm_r_fib(i)
        ds_fib(i,j,jlo,klo)=ds_ibm_fib(i,j)
       end do
      enddo

      return
      end 
     
      subroutine copy_ibm_fsh(
     &   ilo,ihi,jlo,jhi,klo,khi,nscalar
     &  ,ds_fsh
     &  ,scalar_fsh
     &  ,gxyz_ibm_massiveout_fsh)
      USE HeaderFSI
      implicit none

      integer ilo,ihi,jlo,jhi,klo,khi,nscalar
      real*8 GXYZ_IBM_MASSIVEout_fsh(Nr_IBM_fsh
     &   ,ilo:ihi,jlo:jhi,klo:khi,3)
      real*8 ds_fsh(Nr_IBM_fsh
     &   ,ilo:ihi,jlo:jhi,klo:khi)
      real*8 scalar_fsh(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,nscalar)

      integer i,j,jq,i_glb
      real*8 tmpv1val
      real*8 tmpv2val

      if (nscalar.eq.2) then
! do nothing
      else
       print *,"expecting nscalar.eq.2"
       stop
      endif

      do i=1,nr_ibm_fsh

       if ((ilo.le.1-Nq_IBMB_fsh).and. 
     &     (ilo.le.0).and.
     &     (ihi.ge.Nq_IBM_fsh+1).and.
     &     (ihi.ge.Nq_IBM_fsh+Nq_IBMB_fsh).and. &
     &     (jlo.le.1-Ns_IBMB_fsh).and. 
     &     (jlo.le.0).and.
     &     (jhi.ge.Ns_IBM_fsh+1).and. 
     &     (jhi.ge.Ns_IBM_fsh+Ns_IBMB_fsh).and.
     &     (klo.eq.1).and.
     &     (khi.eq.1)) then
!       do nothing
       else
        print *,"ilo,ihi,jlo,jhi,klo,khi invalid"
        stop
       endif

       i_glb=Ilocglb_fsh(i);

       if(Electromechanical(i_glb)) then

        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
         scalar_fsh(i,jq,j,klo,1)=
     &    ndot_ibm_fsh(i,jq,j)
        end do
        end do

        if(Piezo_Dist_Flag(i_glb)) then
         do jq=1,nq_ibm_r_fsh(i)
         do j=1,ns_ibm_r_fsh(i)
          scalar_fsh(i,jq,j,klo,2)=
     &     vPiezo_fsh(i,jq,j)
         enddo
         enddo
        else
         do jq=1,1
         do j=1,1
          scalar_fsh_out(i,jq,j,klo,2)=
     &     vPiezo_fsh(i,jq,j)
         enddo
         enddo
        endif
       endif

       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
        gxyz_ibm_massiveout_fsh(i,jq,j,klo,1)=
     &     gx_ibm_massive_fsh(i,jq,j)
        gxyz_ibm_massiveout_fsh(i,jq,j,klo,2)=
     &     gy_ibm_massive_fsh(i,jq,j)
        gxyz_ibm_massiveout_fsh(i,jq,j,klo,3)=
     &     gz_ibm_massive_fsh(i,jq,j)
        tmpv1val=1.0d0
        tmpv2val=1.0d0
        if(j .eq.1 .or. j .eq. ns_ibm_r_fsh(i)) tmpv1val=0.5d0
        if(jq.eq.1 .or. jq.eq. nq_ibm_r_fsh(i)) tmpv2val=0.5d0
        ds_fsh(i,jq,j,klo)=ds_ibm_fsh(i,jq,j)*ds2_ibm_fsh(i,jq,j)
     &         *tmpv1val*tmpv2val
       enddo
       enddo
      enddo

      return
      end 


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
