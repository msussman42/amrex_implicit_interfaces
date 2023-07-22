c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine create_grid2(ds_ibm0,dr_ibm0,contactlength,theboss)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c this subroutine reads input fuile and create the structures
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      USE HeaderFSI
      implicit none
      integer i,j,j2,i_fin,type_tmp1,i2,jq,k
      integer ixfin(5000)
      real*8 a_r,r,r0,s,s0,alpha,beta,gama,ds_ibm0,dr_ibm0,temp
     &     ,temp2,temp3,temp4,c_prime,deltatmp,ksh_tmp
     &     ,delta_tmp1,delta_tmp2,dx_tmp1
     &     ,delta_tmpx(2),delta_tmpy(2)
     &     ,delta_tmp3,delta_tmp4,lmax_f
     &     ,delta_tmp2z,delta_tmp4z,temp3z,contactlength,contactdx
     &     ,tmpdata(10)

      integer  itmpdata(10),ialphafile,itausfile
      integer  target_ktype
     &                 (nq_ibm_fsh+1,ns_ibm_fbc+ns_ibm_esh+ns_ibm_fsh)
      real*8   target_kvalue
     &                 (nq_ibm_fsh+1,ns_ibm_fbc+ns_ibm_esh+ns_ibm_fsh)
      real*8  target_k_linktmp
     &                 (nq_ibm_fsh+1,ns_ibm_fbc+ns_ibm_esh+ns_ibm_fsh)
     &     ,target_t_linktmp
     &                 (nq_ibm_fsh+1,ns_ibm_fbc+ns_ibm_esh+ns_ibm_fsh)
     &     ,a_fintmp(1,ns_ibm_fbc+ns_ibm_esh+ns_ibm_fsh*nq_ibm_fsh,3)
     &     ,phi_fintmp(1,ns_ibm_fbc+ns_ibm_esh+ns_ibm_fsh*nq_ibm_fsh,3)
     &     ,flagfixedtmp(1,ns_ibm_fbc+ns_ibm_esh+ns_ibm_fsh*nq_ibm_fsh)
      integer i_fintype(ns_ibm_esh+ns_ibm_fsh+ns_ibm_fib),numi,numj
      real*8 maxtmpx,maxtmpy,maxtmpz
     &      ,coefmaxtmpx,coefmaxtmpy,coefmaxtmpz,phi_ben2_coef
      character*34 infile,infilepres,infileFabricMesh,infileFabricLine
      character*150 namebody
      character*10 numalphafile ,numtausfile 
      integer i_phi_mem_coefFLG,i_phi_mem_coef
      real*8  tmpcoord(10),tmpcoordshift(3),tmpcoordR(3),tmpangle(3)
     &       ,tmpcoordscale(6)
      real*8  tmpRmat(3,3)
      real*8   stmp1,stmp2,stmp3,ctmp1,ctmp2,ctmp3,scalegriddx
      integer  in_fib,in_fsh,in_esh,in_fbc
      real*8  Rmattmp(1:3,1:3,1:3)
      integer iseqtmp, i0tmp, i1tmp, i2tmp, i3tmp, tmpseq
      real*8   S_Alpha,coef_j,coef_jq
      logical  fix_bc1(2),fix_bc2(2)
      integer  istartS(2),iraiseS(2)
      logical rampshear,infileflag,theboss
      integer Rot_matPara(3)
      real*8  rhoinf_Coef
      integer method_rhoinf_Coef

      real*8  scale_fbc,scale_fbc2,scale_fbc3

      integer iv0,iv1,iv2,iv3
     &       ,ie0,ie1,ie2,ie3
     &       ,ifound
     &       ,iedgelist
     &       ,ifacepair

      integer conn_face_edge_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3)  
     &       ,conn_face_face_fbc(nr_ibm_fbc,Ns_IBMEle_fbc,3)  
     &       ,faces_vertisin_fbc(nr_ibm_fbc,ns_ibm_fbc,50)  
     &       ,faces_nvertisin_fbc(nr_ibm_fbc,ns_ibm_fbc)

      real*8  maxlength_fbc(ns_ibm_fbc), avlength_fbc(ns_ibm_fbc)
      integer FabricMeshGeneration,FabricLineGeneration
     &       ,numrecord,numrecordL,iMeshmax,ilinemax,tmpint(0:2)

      FabricMeshGeneration=0
      FabricLineGeneration=0
      FabricMesh_Presence=.false.
      FabricLine_Presence=.false.

      numMeshpointRecordMax=0
      numMeshLinkRecordMax=0
      numMeshMax=0
      numMeshpointMax=0
      numMeshLinkMax=0
      FabricMesh_numattMax=0
      FabricMeshL_numattMax=0
      FabricMesh_numtargetMax=0

      numLinepointRecordMax=0
      numLineMax=0
      numLinepointMax=0
      FabricLine_numattMax=0
      FabricLine_numtargetMax=0
     
      open (9997,file='./inputdistibm.dat',status="unknown")


      pi1 = 4. * atan(1.0)
      temp4=0.01
!      alpha=alphaglobal; gama= gammaglobal !gamma here length of body 
!      deltatmp=deltaglobal
!      beta=betaglobal 
!      dx_tmp1=gammaglobal/dble(ns_ibm_fib-1)

      dr_ibm0=0.0
      Genalpha_timesolver=0  !initialization of timemarching method

      if (nr_ibm_fib .gt. 1) dr_ibm0=1./dble(nr_ibm_fib-1)
      ds_ibm0=1.d0/dble(ns_ibm_fib-1)
      r0= dr_ibm0/2. ; s0=0.0*ds_ibm0/2.

      if(theboss) print*, ' ==> Start reading global parameters'

      read(9997,*) k_link, simulationtype, impactpen(1:3), contacttype
     &            ,contactlength,contactdx
      read(9997,*) Ccont(0),Ccont(1),Ccont(2),Ccont(3),Ccont(4)
     &            ,Ccont(5),Ccont(6),Ccont(7),Ccont(8),Ccont(9)
     &            ,Ccont(10)
      read(9997,*) BucketXstart(1:3)
      read(9997,*) BucketXEnd(1:3)
      read(9997,*) BucketdX(1:3)
      read(9997,*) Bucketntinterval
      read(9997,*) delta_typex, delta_typey, delta_typez
      read(9997,*) alpha_ibm , beta_ibm, alphabetaRise
      read(9997,*) save_force_int, del_force,  check_int, fr
      read(9997,*) phi_grav,the_grav
      read(9997,*) gammaglobal
      read(9997,*) maxnstruct,maxErrorstruct,maxnfluid,maxErrorfluid
     &            ,maxiter
      read(9997,*) Solverabstol,Solverreltol,SolverOutiter
     &             ,SolverIniter,SolverILU,SolverPre
      read(9997,*) del_s,alphatemp_ibm,betatemp_ibm
      read(9997,*) channelcontacty,channelcontactz
      read(9997,*) flag_restartfile, t0_restartfile,dt_restartfile
      read(9997,*) 

      contactdx=min(contactdx,min_grid_x)
      contactdx=min(contactdx,min_grid_y)
      contactdx=min(contactdx,min_grid_z)

      min_grid_x=min(min_grid_x,contactlength)
      min_grid_y=min(min_grid_y,contactlength)
      min_grid_z=min(min_grid_z,contactlength)


! ====== intial Contact calculations /Scaling ===================
! scale Ccont for the contact model
!Ccont(0) : angle separation
!Ccont(1) :initial distance check  (scale dx)
!Ccont(2) :false contact calculation   (scale dx)
!Ccont(3) :h distance  (scale dx)
!Ccont(4) :a distance  (scale dx)
!Ccont(5) :b distance  (scale dx)
!Ccont(6) :ntrial
!Ccont(7) :tolx
!Ccont(8) :tolf
!Ccont(9) :initial ratio to ignore contact  for point < Ccont(1)*Ccont(9) apart
!Ccont(10):viscoelastic damping
!Ccont(15) :Contact penalty (scale by 1/dx)

          scalegriddx=contactdx
          Ccont(1)=Ccont(1)*scalegriddx
          Ccont(2)=Ccont(2)*scalegriddx
          Ccont(3)=Ccont(3)*scalegriddx
          Ccont(4)=Ccont(4)*scalegriddx
          Ccont(5)=Ccont(5)*scalegriddx
          Ccont(15)=impactpen(1)/scalegriddx
!bucketCalculation
       do j=1,3
       IBucketnum(j)=int((BucketXEnd(j)-BucketXStart(j))/BucketdX(j))+1
       BucketdX(j)=(BucketXEnd(j)-BucketXStart(j))/dble(IBucketnum(j))
       enddo

       NBucketnum=IBucketnum(1)*IBucketnum(2)*IBucketnum(3)
       if (NBucketnum .gt. NMaxnumbuck) then
              write(*,*)'More than Max Buck=',NMaxnumbuck
     &         ,' have beenused=',IBucketnum
              stop 
       endif


      k_massive_ibm=0.0
      the_grav=the_grav/180.d0*4.d0*atan(1.d0)
      phi_grav=phi_grav/180.d0*4.d0*atan(1.d0)

      read(9997,*) nr_ibm_f,nr_ibm_fb
      if(theboss) print*, 'num of flex             ',nr_ibm_f
      if(theboss) print*, 'num of flex for contact ',nr_ibm_fb
      read(9997,'(a)') infile 
      if(theboss) print*, 'trangulated contact file ',infile
      open(9995,file=infile,status='old',form='formatted')
      read(9997,*) nFixcontact,nFixcontactL

      allocate(xFixcontact(nFixcontact,3),EleFixcontact(nFixcontactL,3))

      do k=1,nFixcontact
             read(9995,*)xFixcontact(k,1:3)
      enddo
      do k=1,nFixcontactL
             read(9995,*)EleFixcontact(k,1:3)
      enddo
      close(9995)
!===============================================================
      lmax_f=gammaglobal

      in_fib=0
      in_fsh=0
      in_esh=0
      in_fbc=0
      do i=0,nr_ibm-1
          i_fix=1

          read(9997,*) 
          read(9997,*)  bodytype(i+1)
!     bodytype =0 fiber, =1 shell, =2 fibrous Shell,  =3 fabric

          if(bodytype(i+1).eq. fiber) then
             in_fib=in_fib+1
             Ilocglb_fib(in_fib)=i+1
             Iglbloc_fib(i+1)=in_fib
          elseif(bodytype(i+1).eq. shell) then
             in_esh=in_esh+1
             Ilocglb_esh(in_esh)=i+1
             Iglbloc_esh(i+1)=in_esh
          elseif(bodytype(i+1).eq. fibrousshell) then
             in_fsh=in_fsh+1
             Ilocglb_fsh(in_fsh)=i+1
             Iglbloc_fsh(i+1)=in_fsh
          elseif(bodytype(i+1).eq. fabric) then
             in_fbc=in_fbc+1
             Ilocglb_fbc(in_fbc)=i+1
             Iglbloc_fbc(i+1)=in_fbc
          endif

          if(theboss) print*,  '==>=>start reading bdy',i+1
          if(theboss) print*, 'BODYtype=',bodytype(i+1)

          read(9997,*) infileflag
          read(9997,'(a)') infile
          if(infileflag) then
          if(theboss) print*,  '==>=>start reading ',infile

          open(9995,file=infile,status='old',form='formatted')
          read(9995,*) nibmptrecord(i+1)
          do k=1,nibmptrecord(i+1)
             read(9995,*)ibmptrecord(i+1,k,1:2)
          enddo
          close(9995)
          else
           nibmptrecord(i+1)=0
          endif

          read(9997,*)  tsibmin(i+1)    !temprature of the body
          if(theboss) print*, 'tsibmin=',tsibmin(i+1)
          read(9997,*)  cs_ibm(i+1),cs_ibm_target(i+1)
     &                 ,density_coef(i+1),density_coefp(i+1)
          if(theboss) print*, 'cs_ibm=',cs_ibm(i+1)
          if(theboss) print*, 'cs_ibm_target=',cs_ibm_target(i+1)
          if(theboss) print*, 'density_coef=',density_coef(i+1)
          if(theboss) print*, 'Virtual Addedmass=',density_coefp(i+1)

          read(9997,*)  FluidForceFlag(i+1),ThermalForceFlag(i+1)
     &       ,ContactForceFlag(i+1),ExternalPressureFlag(i+1)
     &       ,Electromechanical(i+1),ImplicitElecMech(i+1)

          if(theboss) print*, 'FluidForceFlag=',FluidForceFlag(i+1)
          if(theboss) print*, 'ThermalFlag=',ThermalForceFlag(i+1)
          if(theboss) print*, 'ContactFlag=',ContactForceFlag(i+1)
          if(theboss) print*, 'ExtPresFlag=',ExternalPressureFlag(i+1)
          if(theboss) print*, 'EleMechFlag=',Electromechanical(i+1)
          if(theboss) print*, 'Implicit EleMechFlag='
     &                        ,ImplicitElecMech(i+1)

          read(9997,*) contactplanesFlag(i+1)
          if(theboss) print*,'contactplanesFlag=',contactplanesFlag(i+1) 
          if( contactplanesFlag(i+1)) then
             read(9997,*) Ncontactplanes(i+1)
          if(theboss)print*,'Num of Contact planes=',Ncontactplanes(i+1)
             do j=1,Ncontactplanes(i+1)
               if(j .eq. 1) then
                if(theboss)print*,'================================='
               endif
               read(9997,*) contactplanesPara(i+1,j,1:6)
!     X0 and n of plate   
               ctmp1=sqrt(contactplanesPara(i+1,j,4)**2
!     normalized the normal vector
     &                    +contactplanesPara(i+1,j,5)**2 
     &                    +contactplanesPara(i+1,j,6)**2)  
                contactplanesPara(i+1,j,4:6)= 
     &                    contactplanesPara(i+1,j,4:6)/ctmp1 
                if(theboss)print*,contactplanesPara(i+1,j,1:6)        
             enddo
          endif


          read(9997,*) tmpcoordshift(1:3)
          if(theboss)print*,'coordshift=',tmpcoordshift(1:3)
          read(9997,*) tmpcoordscale(1:6)
          if(theboss)print*,'coordscale=',tmpcoordscale(1:6)
          read(9997,*) tmpcoordR(1:3)
          if(theboss)print*,'Rotation center=',tmpcoordR(1:3)
          read(9997,*) tmpangle(1:3),tmpseq
          if(theboss)print*,'Rot angle(deg)=',tmpangle(1:3)
          if(theboss)print*,'Rot sequence=',tmpseq      
!for angle convention look at http://en.wikipedia.org/wiki/Euler_angles
       tmpangle(1:3)=tmpangle(1:3)/180.0d0*pi1
       ctmp1=cos(tmpangle(1))
       ctmp2=cos(tmpangle(2))
       ctmp3=cos(tmpangle(3))
       stmp1=sin(tmpangle(1))
       stmp2=sin(tmpangle(2))
       stmp3=sin(tmpangle(3))


       Rmattmp(1,1,1:3)=(/1.0d+0,   0.0d+0,   0.0d+0/)
       Rmattmp(1,2,1:3)=(/0.0d+0, ctmp1,-stmp1/)              
       Rmattmp(1,3,1:3)=(/0.0d+0, stmp1, ctmp1/)


       Rmattmp(2,1,1:3)=(/ ctmp2,0.0d+0, stmp2/)
       Rmattmp(2,2,1:3)=(/0.0d+0,1.0d+0, 0.0d+0/)              
       Rmattmp(2,3,1:3)=(/-stmp2,0.0d+0, ctmp2/)


       Rmattmp(3,1,1:3)=(/ ctmp3,-stmp3,0.0d+0/)
       Rmattmp(3,2,1:3)=(/ stmp3, ctmp3,0.0d+0/)             
       Rmattmp(3,3,1:3)=(/0.0d+0,0.0d+0,1.0d+0/)

       iseqtmp=mod(tmpseq,10)

       i0tmp=int((tmpseq+1.0e-5)/10.0)
       i1tmp=int(i0tmp/100)
       i3tmp=mod(i0tmp,10)
       i2tmp=mod(int(i0tmp/10),10)

       if(iseqtmp.eq. 0)  then !current coordinate
       tmpRmat(1:3,1:3)=matmul(Rmattmp(i2tmp,1:3,1:3)
     &               ,Rmattmp(i3tmp,1:3,1:3))
       tmpRmat(1:3,1:3)=matmul(Rmattmp(i1tmp,1:3,1:3)
     &               ,tmpRmat(1:3,1:3))
       else !fixed coordinate
       tmpRmat(1:3,1:3)=matmul(Rmattmp(i2tmp,1:3,1:3)
     &               ,Rmattmp(i1tmp,1:3,1:3))
       tmpRmat(1:3,1:3)=matmul(Rmattmp(i3tmp,1:3,1:3)
     &               ,tmpRmat(1:3,1:3))
       endif
       read(9997,*) 
!==================== FIBER ============================================
      if (bodytype(i+1) .eq. fiber) then   !fiber
         read(9997,*)  ks_ibm,  kb_ibm, iflaginext(i+1)
     &               , implicitflag(i+1)
          if(theboss)print*,'ks_ibm,  kb_ibm, iflaginext(i+1)'
     &               , 'implicitflag(i+1)'
          if(theboss)print*, ks_ibm,  kb_ibm, iflaginext(i+1)
     &               , implicitflag(i+1)
          read(9997,*)  type_tmp1
!     1=create fiber mesh 0= read fiber mesh
          if(theboss)print*,' type of fiber creation=',type_tmp1
         read(9997,*)  fiberbdry(i+1)
          if(theboss)print*,' Boundary code=',fiberbdry(i+1)

         if (real(type_tmp1) .gt. 0.5) then
             read(9997,*)  delta_tmp1,delta_tmp2,delta_tmp2z, 
     &                     delta_tmp3,delta_tmp4, 
     &                     delta_tmp4z
              if(theboss)print*,' coordinate of the first point='
     &                          ,delta_tmp1,delta_tmp2,delta_tmp2z
              if(theboss)print*,' coordinate of the last point='
     &                          ,delta_tmp3,delta_tmp4,delta_tmp4z
             if(ndim .lt. 3) then
               delta_tmp2z=0.0
               delta_tmp4z=0.0
             endif
 
              masterbdy(i+1)=0

              a_r  =sqrt( (delta_tmp1-delta_tmp3  )**2+
     &                    (delta_tmp2-delta_tmp4  )**2+
     &                    (delta_tmp2z-delta_tmp4z)**2)

              ns_ibm_r(i+1)=int(dble(ns_ibm_fib)/lmax_f*a_r)+0
              ns_ibm_r_fib(in_fib)= ns_ibm_r(i+1)
              ns_ibm_i_fib(in_fib)=0

              nq_ibm_r(i+1)=1
              ds_ibm0=1.d0/dble(ns_ibm_r(i+1)-1)
              s0=0.0*ds_ibm0/2.
              r=r0+(i)*dr_ibm0
              a_r=1.0   ! +dsin(2*pi1*r-pi1/2.)
c variable stiffness
              ks_ibm_r(i+1)=ks_ibm*a_r
              kb_ibm_r(i+1)=kb_ibm*a_r
              k_massive_ibm_r(i+1)=k_massive_ibm*a_r
              temp=0.5*(alpha+beta)
            if (create_new_flow) then
              gama= gammaglobal !gamma here length of body 

              gx_ibm_fib(in_fib,1)=delta_tmp1
              gy_ibm_fib(in_fib,1)=delta_tmp2
              gz_ibm_fib(in_fib,1)=delta_tmp2z


              read(9997,*) frequency_fin(i+1,1:3)
              read(9997,*) i_fin
              read(9997,*) ixfin(1:i_fin)    
              do j=1,i_fin
                    if(ixfin(j) .le. 0) then
                        target_point_num_fib(in_fib,j)
     &                                =ns_ibm_r_fib(in_fib)+ixfin(j)+0
                    else
                        target_point_num_fib(in_fib,j)=ixfin(j)+0
                    endif
                    i_fix=i_fix+1
              end do

              do j=0,ns_ibm_r_fib(in_fib)

                 s=s0+dble(j)*ds_ibm0
                 gx_ibm_fib(in_fib,j+1)=delta_tmp1
     &                     +(delta_tmp3-delta_tmp1)*s
                 gy_ibm_fib(in_fib,j+1)=delta_tmp2
     &                +(delta_tmp4-delta_tmp2)*s
!     dsin(2.0d0*pi1*s)
                 gz_ibm_fib(in_fib,j+1)=delta_tmp2z
     &                +(delta_tmp4z-delta_tmp2z)*s
!     dsin(2.0d0*pi1*s)

                 tmpcoord(1:3)=
     &           (/gx_ibm_fib(in_fib,j+1)
     &           , gy_ibm_fib(in_fib,j+1)
     &           , gz_ibm_fib(in_fib,j+1) /)

                 do jq=1,3
                     tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                 enddo

                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)

                 gx_ibm_fib(in_fib,j+1)=tmpcoord(1)
                 gy_ibm_fib(in_fib,j+1)=tmpcoord(2)
                 gz_ibm_fib(in_fib,j+1)=tmpcoord(3)

                 gx_bp_fib(in_fib,j+1)=gx_ibm_fib(in_fib,j+1)
                 gy_bp_fib(in_fib,j+1)=gy_ibm_fib(in_fib,j+1)
                 gz_bp_fib(in_fib,j+1)=gz_ibm_fib(in_fib,j+1)
                if (j .eq. 0) then 
                   temp2=gx_bp_fib(in_fib,1)-gx_ibm_fib(in_fib,1)
                   temp3=gy_bp_fib(in_fib,1)-gy_ibm_fib(in_fib,1)
                   temp3z=gz_bp_fib(in_fib,1)-gz_ibm_fib(in_fib,1)
                end if
                 gx_ibm_fib(in_fib,j+1)=temp2+gx_ibm_fib(in_fib,j+1)
                 gy_ibm_fib(in_fib,j+1)=temp3+gy_ibm_fib(in_fib,j+1)
                 gz_ibm_fib(in_fib,j+1)=temp3z+gz_ibm_fib(in_fib,j+1)

                 gx_bp_fib(in_fib,j+1)=gx_ibm_fib(in_fib,j+1)
                 gy_bp_fib(in_fib,j+1)=gy_ibm_fib(in_fib,j+1)
                 gz_bp_fib(in_fib,j+1)=gz_ibm_fib(in_fib,j+1)
                 gx_bp_fib0(in_fib,j+1)=gx_ibm_fib(in_fib,j+1)
                 gy_bp_fib0(in_fib,j+1)=gy_ibm_fib(in_fib,j+1)
                 gz_bp_fib0(in_fib,j+1)=gz_ibm_fib(in_fib,j+1)
                 gx_ibmpre_fib(in_fib,j+1)=gx_ibm_fib(in_fib,j+1)
                 gy_ibmpre_fib(in_fib,j+1)=gy_ibm_fib(in_fib,j+1)
                 gz_ibmpre_fib(in_fib,j+1)=gz_ibm_fib(in_fib,j+1)
                 if (massive_ibm) mass_ibm_fib(in_fib,j+1)=1.0*ds_ibm0
              end do
            end if
               target_num(i+1)=i_fix-1
               read(9997,*) target_k_link_fib(in_fib,1:target_num(i+1)) 
               read(9997,*) target_t_link_fib(in_fib,1:target_num(i+1))   
               read(9997,*) a_fin_fib(in_fib,1:target_num(i+1),1)    
               read(9997,*) phi_fin_fib(in_fib,1:target_num(i+1),1)    
               read(9997,*) a_fin_fib(in_fib,1:target_num(i+1),2)    
               read(9997,*) phi_fin_fib(in_fib,1:target_num(i+1),2)
               read(9997,*) a_fin_fib(in_fib,1:target_num(i+1),3)    
               read(9997,*) phi_fin_fib(in_fib,1:target_num(i+1),3)  
               read(9997,*) flagfixed_fib(in_fib,1:target_num(i+1))   
            else
                a_r=1.0   ! +dsin(2*pi1*r-pi1/2.)
c variable stiffness
                ks_ibm_r(i+1)=ks_ibm*a_r
                kb_ibm_r(i+1)=kb_ibm*a_r
                k_massive_ibm_r(i+1)=k_massive_ibm*a_r
                i2=i+1

                read(9997,*)  fiberbdry(i+1)
                read(9997,'(a)') infile
                read(9997,*) frequency_fin(i+1,1:3)
                read(9997,*) i_fin
                read(9997,*) ixfin(1:i_fin) 
                i_fix=0;   
                do j=1,i_fin
                   i_fix=i_fix+1
                end do
               target_num(i+1)=i_fix-1
               read(9997,*) target_k_link_fib(in_fib,1:target_num(i+1)) 
               read(9997,*) target_t_link_fib(in_fib,1:target_num(i+1))   
               read(9997,*) a_fin_fib(in_fib,1:target_num(i+1),1)    
               read(9997,*) phi_fin_fib(in_fib,1:target_num(i+1),1)    
               read(9997,*) a_fin_fib(in_fib,1:target_num(i+1),2)    
               read(9997,*) phi_fin_fib(in_fib,1:target_num(i+1),2)
               read(9997,*) a_fin_fib(in_fib,1:target_num(i+1),3)    
               read(9997,*) phi_fin_fib(in_fib,1:target_num(i+1),3)
               read(9997,*) flagfixed_fib(in_fib,1:target_num(i+1))


               open(9998,file=infile,status='old',form='formatted')
               read(9998,'(a)') namebody(i*10+1:(i2)*10)
               read(9998,*) ns_ibm_r(i2),ns_ibm_i(i2) 
     &                     ,ns_ibm_rall(i2)
     &                     ,masterbdy(i2)

                ns_ibm_r_fib(in_fib)=ns_ibm_r(i2)
                ns_ibm_rall_fib(in_fib)=ns_ibm_rall(i2)
                ns_ibm_i_fib(in_fib)=ns_ibm_i(i2)

                do j=1,i_fin
                    if(ixfin(j) .le. 0) then
                        target_point_num_fib(in_fib,j)=
     &                                ns_ibm_r_fib(in_fib)+ixfin(j)+0
                    else
                        target_point_num_fib(in_fib,j)=ixfin(j)+0
                    endif
                end do
 
           
                 do j=1,ns_ibm_rall_fib(in_fib)
                      read(9998,*)   
     &                     numi   
     &                    ,gx_ibm_fib(in_fib,j)   
     &                    ,gy_ibm_fib(in_fib,j)   
     &                    ,gz_ibm_fib(in_fib,j) 

                 tmpcoord(1:3)=
     &           (/gx_ibm_fib(in_fib,j)
     &           , gy_ibm_fib(in_fib,j)
     &           , gz_ibm_fib(in_fib,j) /)

                 do jq=1,3
                     tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                 enddo


                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)
                 gx_ibm_fib(in_fib,j)=tmpcoord(1)
                 gy_ibm_fib(in_fib,j)=tmpcoord(2)
                 gz_ibm_fib(in_fib,j)=tmpcoord(3)

                 end do

                if (masterbdy(i2) .gt. 0) then
                   do j=1,ns_ibm_r_fib(in_fib)
                      read(9998,*)  
     &                    numi, npos_fib(in_fib,j)
                   end do

                 end if
                 do j=1,ns_ibm_r_fib(in_fib)
                    gx_bp_fib(in_fib,j+1)=gx_ibm_fib(in_fib,j+1)
                    gy_bp_fib(in_fib,j+1)=gy_ibm_fib(in_fib,j+1)
                    gz_bp_fib(in_fib,j+1)=gz_ibm_fib(in_fib,j+1)
                    if (j .eq. 0) then 
                       temp2=gx_bp_fib(in_fib,1)-gx_ibm_fib(in_fib,1)
                       temp3=gy_bp_fib(in_fib,1)-gy_ibm_fib(in_fib,1)
                       temp3z=gz_bp_fib(in_fib,1)-gz_ibm_fib(in_fib,1)
                    end if
                    gx_ibm_fib(in_fib,j+1)=temp2+gx_ibm_fib(in_fib,j+1)
                    gy_ibm_fib(in_fib,j+1)=temp3+gy_ibm_fib(in_fib,j+1)
                    gz_ibm_fib(in_fib,j+1)=temp3z+gz_ibm_fib(in_fib,j+1)
                    gx_bp_fib(in_fib,j+1)=gx_ibm_fib(in_fib,j+1)
                    gy_bp_fib(in_fib,j+1)=gy_ibm_fib(in_fib,j+1)
                    gz_bp_fib(in_fib,j+1)=gz_ibm_fib(in_fib,j+1)

                    gx_bp_fib0(in_fib,j+1)=gx_ibm_fib(in_fib,j+1)
                    gy_bp_fib0(in_fib,j+1)=gy_ibm_fib(in_fib,j+1)
                    gz_bp_fib0(in_fib,j+1)=gz_ibm_fib(in_fib,j+1)

                    gx_ibmpre_fib(in_fib,j+1)=gx_ibm_fib(in_fib,j+1)
                    gy_ibmpre_fib(in_fib,j+1)=gy_ibm_fib(in_fib,j+1)
                    gz_ibmpre_fib(in_fib,j+1)=gz_ibm_fib(in_fib,j+1)
                    if (massive_ibm) mass_ibm_fib(in_fib,j+1)=0.0
                 enddo

                close (9998)

             endif
!     end of whether to make the mesh or read that for fibers

               read(9997,*) flag_TDpt_Pnormal_Ext(i+1),Pnormal_Ext(i+1)
               read(9997,'(a)') infilepres
               if (flag_TDpt_Pnormal_Ext(i+1) .eq. 1) then
                  open(9998,
     &                 file=infilepres,status='old',form='formatted' )
                  read(9998,*) n_TDpt_Pnormal_Ext(i+1)
                  if (n_TDpt_Pnormal_Ext(i+1) .gt. 1000) then
                     write(*,*) 'size of press is larger than allowed'
                     write(*,*) 'For body', i+1


                     stop
                  endif
                  do j=1,n_TDpt_Pnormal_Ext(i+1)
                     read(9998,*) data_TDpt_Pnormal_Ext(i+1,j,1)
     &                           ,data_TDpt_Pnormal_Ext(i+1,j,2)
                  enddo
                  close (9998)
               endif
!============= SHELL FE ================================================
          elseif (bodytype(i+1) .eq. shell) then    ! shell input read 
                i2=i+1
                read(9997,'(a)') infile
                read(9997,*) ShellModelType_esh(in_esh)   
     &                      ,MaterialTypeIBM(i2)   
     &                      ,i_phi_mem_coef,Fibrous_esh(in_esh)
     &               ,comprerssibleflag_esh(in_esh)
!     !NEW INPUT
     &               ,flag_fractal_esh(in_esh)
!     !====INSERT FRACTAL
                if(flag_fractal_esh(in_esh)) then
                      call DistFSI_allocate_memory(1)
                      read(9997,*) fractalpara_esh(in_esh,1:9)
!     alphaS11,tausS11,taueS11,alphaS22,tausS22,taueS22,alphaS12,tausS12,taueS12
!     !====INSERT FRACTAL
                endif
                read(9997,*) phi_mem_coef(i2,1:i_phi_mem_coef)   
                read(9997,*) phi_ben_coef(i2)   
     &                      ,phi_nus_coef(i2)
               read(9997,*) Genalpha_timesolver(i2)
     &                     ,Genalpha_niter(i2)
     &                     ,propDamp(i2,1:2)
               if(Genalpha_timesolver(i2) .gt. 5) then
                  write(*,*) 'Genalpha_timesolver for', i2,'is >5?'
                  stop
               endif
               read(9997,*) rhoinf_Coef,method_rhoinf_Coef
     &                     ,alpham_Coef(i2),alphaf_Coef(i2)

               if(method_rhoinf_Coef .eq. 2) then
                  !rhoinf_Coef=1.0/rhoinf_Coef
                  alpham_coef(i2)=(2.0*rhoinf_Coef-1)/(1.0+rhoinf_Coef)
!     2nd order system
                  alphaf_coef(i2)=rhoinf_Coef/(1.0+rhoinf_Coef)
               elseif(method_rhoinf_Coef .eq. 1) then
                  alpham_coef(i2)=
     &                 0.5*(3.0*rhoinf_Coef-1)/(1.0+rhoinf_Coef)
!     1st order system
                  alphaf_coef(i2)=rhoinf_Coef/(1.0+rhoinf_Coef)
               endif
               beta_Coef(i2)=0.25*
     &                        (1.0-alpham_coef(i2)+alphaf_coef(i2))**2
               gamma_Coef(i2)=0.5-alpham_coef(i2)+alphaf_coef(i2)


              if(flag_fractal_esh(in_esh)) then
                  do j=1,3
                   tmpdata(1:3)=
     &           (/fractalpara_esh(in_esh,(j-1)*3+1)
     &           , fractalpara_esh(in_esh,(j-1)*3+2)
     &           , fractalpara_esh(in_esh,(j-1)*3+3) /)     


!                   j2=int((tmpcoord(1)+0.0001)/0.1)
!                   if((j2 .lt. 1) .or. (j2 .gt. 9) ) 
!     &         stop 'there is no more fractal database for your choice'
                  itmpdata(1)=int(tmpdata(1)*10+0.0001)
                  itmpdata(2)=int(tmpdata(2)*100+0.0001)

                  if(itmpdata(1) .lt. 10) then
	            write(numalphafile,'(i1.1)') itmpdata(1)
                  ialphafile=1
                  elseif(itmpdata(1) .lt. 100) then
	            write(numalphafile,'(i2.2)') itmpdata(1)   
                  ialphafile=2
                  endif

                  if(itmpdata(2) .lt. 10) then
	            write(numtausfile,'(i1.1)') itmpdata(2)
                  itausfile=1
                  elseif(itmpdata(2) .lt. 100) then
	            write(numtausfile,'(i2.2)') itmpdata(2)  
                  itausfile=2 
                  elseif(itmpdata(2) .lt. 1000) then
	            write(numtausfile,'(i3.3)') itmpdata(2)   
                  itausfile=3
                  elseif(itmpdata(2) .lt. 10000) then
	            write(numtausfile,'(i4.4)') itmpdata(2)  
                  itausfile=4 
                  endif
                  print*,numtausfile,itausfile, itmpdata(2)
                  infilepres='./FractalDatabase/'
     &                      //numalphafile(1:ialphafile)
     &                      //'_'
     &                      //numtausfile(1:itausfile)//'.txt'
                  open(9998,
     &                 file=infilepres,status='old',form='formatted' ) 
                  read(9998,*) !alpha
                  read(9998,*) !taus
                  nifractal_esh(in_esh)=12
                  if (itmpdata(1).eq. 10) nifractal_esh(in_esh)=1
                  jq=nifractal_esh(in_esh)
                  do j2=1,jq
                     read(9998,*) Afractal_esh(j2,1:jq,j,in_esh)   
                  enddo
                  do j2=1,jq
                     read(9998,*) Bfractal_esh(j2,j,in_esh)   
                  enddo                                
                  read(9998,*) Cfractal_esh(1:jq,j,in_esh) 
                  read(9998,*) Dfractal_esh(j,in_esh)  
                  close(9998)

                  fractalpara_esh(in_esh,(j-1)*3+2)=
     &                   ( tmpdata(3)/tmpdata(2))**tmpdata(1)

                  fractalpara_esh(in_esh,(j-1)*3+3)=
     &                   1.0/tmpdata(2)**tmpdata(1)*
     &                   ( 1.0-fractalpara_esh(in_esh,(j-1)*3+2))
                 enddo
               endif

               read(9997,*) frequency_fin(i+1,1:3)
               read(9997,*) i_fintype(i+1)
               read(9997,*) target_k_linktmp(1,1:i_fintype(i+1)) 
               read(9997,*) target_t_linktmp(1,1:i_fintype(i+1))  
 
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),1)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),1)   
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),2)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),2)    
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),3)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),3)    
 
               read(9997,*) flagfixedtmp(1,1:i_fintype(i+1)) 
               read(9997,*) flag_TDpt_Pnormal_Ext(i+1),Pnormal_Ext(i+1)
               read(9997,'(a)') infilepres
               if (flag_TDpt_Pnormal_Ext(i+1) .eq. 1) then
                  open(9998,
     &                 file=infilepres,status='old',form='formatted' )
                  read(9998,*) n_TDpt_Pnormal_Ext(i+1)
                  if (n_TDpt_Pnormal_Ext(i+1) .gt. 1000) then
                     write(*,*) 'size of press is larger than allowed'
                     write(*,*) 'For body', i+1
                     stop
                  endif
                  do j=1,n_TDpt_Pnormal_Ext(i+1)
                     read(9998,*) data_TDpt_Pnormal_Ext(i+1,j,1)
     &                           ,data_TDpt_Pnormal_Ext(i+1,j,2)
                  enddo
                  close (9998)
               endif

               open(9998,file=infile,status='old',form='formatted' )
               read(9998,'(a)' ) namebody(i*10+1:(i2)*10)
               read(9998,*)   
     &                    ns_ibm_r(i2),ns_ibm_i(i2), ns_ibm_rall(i2) 
     &                   ,numelr(i2), numeli(i2) 
     &                   ,masterbdy(i2)
               read(9998,*) Rot_matPara(1:3), n_matpara_esh(in_esh)

                ns_ibm_r_esh(in_esh)=ns_ibm_r(i2)
                ns_ibm_rall_esh(in_esh)=ns_ibm_rall(i2)
                ns_ibm_i_esh(in_esh)=ns_ibm_i(i2)

               nq_ibm_r(i2)=0
                read(9998,*) jtypeall(1:jtypemax)
               n_jtypeall=0
               do j=1,jtypemax
                  if(jtypeall(j)>-1) then 
                      n_jtypeall=n_jtypeall+1
                  end if
               end do

               j2=0
               maxtmpx=0.0
               maxtmpy=0.0
               maxtmpz=0.0
               do j=1,ns_ibm_r_esh(in_esh)
                  j2=j2+1
                  if(n_matpara_esh(in_esh) .gt. 0) then
                    read(9998,*) numi   
     &                    ,gx_ibm_esh(in_esh,j)   
     &                    ,gy_ibm_esh(in_esh,j)   
     &                    ,gz_ibm_esh(in_esh,j)   
     &                    ,target_ktype(1,j)
     &                    ,target_kvalue(1,j)
     &             ,Materialpara_esh(in_esh,j,1:n_matpara_esh(in_esh))
                  else
                    read(9998,*) numi   
     &                    ,gx_ibm_esh(in_esh,j)   
     &                    ,gy_ibm_esh(in_esh,j)   
     &                    ,gz_ibm_esh(in_esh,j)   
     &                    ,target_ktype(1,j)
     &                    ,target_kvalue(1,j)
                  endif

                  target_ktype_esh(in_esh,j)=target_ktype(1,j)
                  target_kvalue_esh(in_esh,j)=target_kvalue(1,j)

                 tmpcoord(1:3)=
     &           (/gx_ibm_esh(in_esh,j)
     &           , gy_ibm_esh(in_esh,j)
     &           , gz_ibm_esh(in_esh,j) /)

                 do jq=1,3
                     tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                 enddo
                 jq=1

                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)

                 gx_ibm_esh(in_esh,j)=tmpcoord(1)
                 gy_ibm_esh(in_esh,j)=tmpcoord(2)
                 gz_ibm_esh(in_esh,j)=tmpcoord(3)

                 if(Rot_matPara(1) .gt. 0) then
                 tmpcoord(1:3)=
     &           (/Materialpara_esh(in_esh,j,Rot_matPara(1))
     &           , Materialpara_esh(in_esh,j,Rot_matPara(2))
     &           , Materialpara_esh(in_esh,j,Rot_matPara(3)) /)                     

                 tmpcoord(1:3)=matmul(tmpRmat,tmpcoord(1:3))

                 Materialpara_esh(in_esh,j,Rot_matPara(1))=tmpcoord(1)
                 Materialpara_esh(in_esh,j,Rot_matPara(2))=tmpcoord(2)
                 Materialpara_esh(in_esh,j,Rot_matPara(3))=tmpcoord(3)
                 endif
                  maxtmpx=maxtmpx+abs(gx_ibm_esh(in_esh,j))
                  maxtmpy=maxtmpy+abs(gy_ibm_esh(in_esh,j))
                  maxtmpz=maxtmpz+abs(gz_ibm_esh(in_esh,j))
               end do
               do j=1,ns_ibm_i_esh(in_esh)
                      j2=j2+1
                      read(9998,*) numi,   
     &                    gx_ibm_esh(in_esh,j2),   
     &                    gy_ibm_esh(in_esh,j2),   
     &                    gz_ibm_esh(in_esh,j2), 
     &                    bndycnd_esh(in_esh,j2),
     &                    boundrylist_esh(in_esh,j,1),   
     &                    boundrylist_esh(in_esh,j,2),   
     &                    boundrylist_esh(in_esh,j,3)

                 tmpcoord(1:3)=
     &           (/gx_ibm_esh(in_esh,j2)
     &           , gy_ibm_esh(in_esh,j2)
     &           , gz_ibm_esh(in_esh,j2) /)

                 do jq=1,3
                     tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                 enddo
                 jq=1
               
                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)
                 gx_ibm_esh(in_esh,j2)=tmpcoord(1)
                 gy_ibm_esh(in_esh,j2)=tmpcoord(2)
                 gz_ibm_esh(in_esh,j2)=tmpcoord(3)

                 do jq=1,n_matpara_esh(in_esh)
                Materialpara_esh(in_esh,j2,jq)=
     &          Materialpara_esh(in_esh,boundrylist_esh(in_esh,j,1),jq)        
     &         +Materialpara_esh(in_esh,boundrylist_esh(in_esh,j,2),jq) 
     &         -Materialpara_esh(in_esh,boundrylist_esh(in_esh,j,3),jq) 
                enddo
               enddo
               do j=ns_ibm_r_esh(in_esh)+ns_ibm_i_esh(in_esh)+1
     &             ,ns_ibm_rall_esh(in_esh)
                      j2=j2+1
                      read(9998,*) numi,   
     &                    gx_ibm_esh(in_esh,j2),   
     &                    gy_ibm_esh(in_esh,j2),   
     &                    gz_ibm_esh(in_esh,j2) 

                 tmpcoord(1:3)=
     &           (/gx_ibm_esh(in_esh,j2)
     &           , gy_ibm_esh(in_esh,j2)
     &           , gz_ibm_esh(in_esh,j2) /)

                 do jq=1,3
                     tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                 enddo
                 jq=1

                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)
                 gx_ibm_esh(in_esh,j2)=tmpcoord(1)
                 gy_ibm_esh(in_esh,j2)=tmpcoord(2)
                 gz_ibm_esh(in_esh,j2)=tmpcoord(3)

                end do
               if (masterbdy(i2) .gt. 0) then
                  do j=1,ns_ibm_r_esh(in_esh)+ns_ibm_i_esh(in_esh)
                      read(9998,*) numi, npos_esh(in_esh,j) 
                  end do
               end if
                i_fix=0             
               do j=1,ns_ibm_r_esh(i2)
                  if(target_ktype(1,j) .ne. 0) then
                    i_fix=i_fix+1
                    target_point_num_esh(in_esh,i_fix)=j
                    target_k_link_esh(in_esh,i_fix)=
     &                            target_k_linktmp(1,target_ktype(1,j))
     &                           *target_kvalue(1,j)
                    target_t_link_esh(in_esh,i_fix)=
     &                            target_t_linktmp(1,target_ktype(1,j))
                    a_fin_esh(in_esh,i_fix,1:3)=
     &                            a_fintmp(1,target_ktype(1,j),1:3)
                    phi_fin_esh(in_esh,i_fix,1:3)=
     &                            phi_fintmp(1,target_ktype(1,j),1:3)
                    flagfixed_esh(in_esh,i_fix)=
     &                            flagfixedtmp(1,target_ktype(1,j))
                  end if     
               end do
               target_num(i+1)=i_fix             
               do j=1,ns_ibm_r_esh(in_esh)
                  read(9998,*) nposele_esh(in_esh,j,1:maxordern-6+2)
               end do

               do j=1,numelr(i2)
                  read(9998,*) neatype_esh(in_esh,j)
     &                        , nea_esh(in_esh,j,1:maxordern)
               end do

               coefmaxtmpx=0.0
               if(maxtmpx .lt. 1d-6) then 
                  coefmaxtmpx=1.0d-6
               endif

               coefmaxtmpy=0.0
               if(maxtmpy .lt. 1d-6) then 
                  coefmaxtmpy=1.0d-6
               endif

               coefmaxtmpz=0.0
               if(maxtmpz .lt. 1d-6) then 
                  coefmaxtmpz=1.0d-6
               endif

               do j=1,ns_ibm_rall_esh(in_esh)
                 gx_ibm_esh(in_esh,j)=gx_ibm_esh(in_esh,j)+ coefmaxtmpx
                 gy_ibm_esh(in_esh,j)=gy_ibm_esh(in_esh,j)+ coefmaxtmpy
                 gz_ibm_esh(in_esh,j)=gz_ibm_esh(in_esh,j)+ coefmaxtmpz

                 gx_bp_esh(in_esh,j)=gx_ibm_esh(in_esh,j)
                 gy_bp_esh(in_esh,j)=gy_ibm_esh(in_esh,j)
                 gz_bp_esh(in_esh,j)=gz_ibm_esh(in_esh,j)
                 gx_bp_esh0(in_esh,j)=gx_ibm_esh(in_esh,j)
                 gy_bp_esh0(in_esh,j)=gy_ibm_esh(in_esh,j)
                 gz_bp_esh0(in_esh,j)=gz_ibm_esh(in_esh,j)

                 gx_ibmpre_esh(in_esh,j)=gx_ibm_esh(in_esh,j)
                 gy_ibmpre_esh(in_esh,j)=gy_ibm_esh(in_esh,j)
                 gz_ibmpre_esh(in_esh,j)=gz_ibm_esh(in_esh,j)
               end do
                close (9998)
!============= Fabric FE ===============================================
          elseif (bodytype(i+1) .eq. fabric) then    ! shell input read 
                i2=i+1
                scale_fbc
     &               =(tmpcoordscale(4)
     &                 *tmpcoordscale(5)
     &                 *tmpcoordscale(6))**(1/3)
                scale_fbc2=scale_fbc*scale_fbc
                scale_fbc3=scale_fbc*scale_fbc*scale_fbc 

                read(9997,'(a)') infile
         
                read(9997,*) inplane_flag_fbc(in_fbc)
     &                      ,inplane_type_fbc(in_fbc)

                if (inplane_type_fbc(in_fbc) .eq. 1) then
                read(9997,*) kappa_inplane_fbc(in_fbc,1:3)
     &                      ,inplane_Cq_fbc(in_fbc,1:4)
     &                      ,InplaneLinemode_fbc(in_fbc,1:3)
     &                      ,InplaneLinem_fbc(in_fbc)
                else
                read(9997,*) kappa_inplane_fbc(in_fbc,1)
     &                      ,InplaneLinemode_fbc(in_fbc,1)
     &                      ,inplane_E_fbc(in_fbc)
     &                      ,inplane_nu_fbc(in_fbc)
                endif


                
                read(9997,*) volume_flag_fbc(in_fbc)
     &                      ,volume_eq_fbc(in_fbc)
     &                      ,kappa_volume_fbc(in_fbc)
                volume_eq_fbc(in_fbc)=volume_eq_fbc(in_fbc)*scale_fbc3
          
                read(9997,*) Aera_flag_fbc(in_fbc)
     &                      ,area_eq_fbc(in_fbc)
     &                      ,kappa_area_fbc(in_fbc)
                area_eq_fbc(in_fbc)=area_eq_fbc(in_fbc)*scale_fbc2
          
                read(9997,*) Bending_flag_fbc(in_fbc)
     &                      ,tet0_fbc(in_fbc)
     &                      ,kappa_bend_fbc(in_fbc)
          
                read(9997,*) edge_length_eq_fbc(in_fbc)
     &                      ,edge_length_max_fbc(in_fbc)
     &                      ,edge_length_p_fbc(in_fbc)
                edge_length_eq_fbc(in_fbc)
     &               =edge_length_eq_fbc(in_fbc)*scale_fbc
                edge_length_max_fbc(in_fbc)
     &               =edge_length_max_fbc(in_fbc)*scale_fbc
                edge_length_p_fbc(in_fbc)
     &               =edge_length_p_fbc(in_fbc)*scale_fbc
          
                read(9997,*) vertex_mass_fbc(in_fbc)
     &                      , KBT_fbc(in_fbc)

               read(9997,*) Genalpha_timesolver(i2)
     &                     ,Genalpha_niter(i2) 
     &                     ,propDamp(i2,1:2)

               if(Genalpha_timesolver(i2) .gt. 5) then
                  write(*,*) 'Genalpha_timesolver for', i2,'is >5?'
                  stop
               endif          
               read(9997,*) rhoinf_Coef,method_rhoinf_Coef
     &                     ,alpham_Coef(i2),alphaf_Coef(i2)

               if(method_rhoinf_Coef .eq. 2) then
                  rhoinf_Coef=1.0/rhoinf_Coef
                  alpham_coef(i2)=(2.0-rhoinf_Coef)/(1.0+rhoinf_Coef)
!     2nd order system
                  alphaf_coef(i2)=1.0/(1.0+rhoinf_Coef)
               elseif(method_rhoinf_Coef .eq. 1) then
                  alpham_coef(i2)=
     &                 0.5*(3.0-rhoinf_Coef)/(1.0+rhoinf_Coef)
!     1st order system
                  alphaf_coef(i2)=1.0/(1.0+rhoinf_Coef)
               endif
               beta_Coef(i2)=0.25*
     &                        (1.0-alpham_coef(i2)+alphaf_coef(i2))**2
               gamma_Coef(i2)=0.5-alpham_coef(i2)+alphaf_coef(i2)

               read(9997,*) FabricMeshFlag_fbc(in_fbc)
                  if(FabricMeshFlag_fbc(in_fbc)) then
                  if(FabricMeshGeneration .eq. 0) then
                  
                  open(9998
     &                ,file='FabricMeshLineSize.txt'
     &                ,status='old',form='formatted' )
                  read(9998,*)
                  read(9998,*)
     &                       numMeshpointRecordMax
     &                      ,numMeshLinkRecordMax
     &                      ,numMeshMax
     &                      ,numMeshpointMax
     &                      ,numMeshLinkMax
     &                      ,FabricMesh_numattMax
     &                      ,FabricMeshL_numattMax
     &                      ,FabricMesh_numtargetMax
                  close (9998)
                  FabricMesh_numattMax=FabricMesh_numattMax+1
!     append zero column to end
                  FabricMeshL_numattMax=FabricMeshL_numattMax+1
!     append zero column to end
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

                  FabricMeshGeneration=1
                  FabricMesh_Presence=.true.
                  endif
               
                  read(9997,'(a)') infileFabricMesh
                  
                  open(9998,
     &           file=infileFabricMesh,status='old',form='formatted' )

                  read(9998,*) numrecord, FabricMesh_numatt(in_fbc)
     &                      ,numrecordL, FabricMeshL_numatt(in_fbc)
                  FabricMesh_npoint(in_fbc,1:numMeshMax)=0
                  FabricMesh_nLine(in_fbc,1:numMeshMax)=0
                  FabricMesh_numrecord(in_fbc)=numrecord
                  FabricMeshL_numrecord(in_fbc)=numrecordL

                  if(numrecord .gt. numMeshpointRecordMax) then
                      write(*,*) 'not enough numMeshpointRecordMax'
                      stop
                  endif
                  if(numrecordL .gt. numMeshLinkRecordMax) then
                      write(*,*) 'not enough numMeshLinkRecordMax'
                      stop
                  endif

                  if(FabricMesh_numatt(in_fbc) 
     &                     .gt. FabricMesh_numattMax) then
                      write(*,*) 'not enough FabricMesh_numattMax'
                      stop
                  endif
                  if(FabricMeshL_numatt(in_fbc) 
     &                     .gt. FabricMeshL_numattMax) then
                      write(*,*) 'not enough FabricMeshL_numattMax'
                      stop
                  endif
                  read(9998,*) 
     &                      iMm0,iMc0,iMconR0,iMcon0,iMs0
     &                     ,iMpD0_3,iMpD0_1,iMpD0_2
                  read(9998,*) n_iMnpara0
                  if( n_iMnpara0 >0) read(9998,*) iMnpara0(1:n_iMnpara0)
                  if(iMm0<=0) iMm0=FabricMesh_numattMax
                  if(iMc0<=0) iMc0=FabricMesh_numattMax
                  if(iMconR0<=0) iMconR0=FabricMesh_numattMax
                  if(iMcon0<=0) iMcon0=FabricMesh_numattMax
                  if(iMs0<=0) iMs0=FabricMesh_numattMax
                  if(iMpD0_1<=0) iMpD0_1=FabricMesh_numattMax
                  if(iMpD0_2<=0) iMpD0_2=FabricMesh_numattMax
                  if(iMpD0_3<=0) iMpD0_3=FabricMesh_numattMax
                  do j=1,n_iMnpara0     
                    if(iMnpara0(j)<0) iMnpara0(j)=FabricMesh_numattMax
                  enddo
                  k=FabricMesh_numattMax
                  if(iMm0>k) then
                        write(*,*) 'not enough attr record for iMm0'
                        stop 
                  endif                 
                  if(iMc0>k) then
                        write(*,*) 'not enough attr record for iMc0'
                        stop 
                  endif                 
                  if(iMconR0>k) then
                        write(*,*) 'not enough attr record for iMconR0'
                        stop 
                  endif                 
                  if(iMcon0>k) then
                        write(*,*) 'not enough attr record for iMcon0'
                        stop 
                  endif                 
                  if(iMs0>k) then
                        write(*,*) 'not enough attr record for iMs0'
                        stop 
                  endif 
                  if(iMpD0_1>k) then
                        write(*,*) 'not enough attr record for iMpD0_1'
                        stop 
                  endif 
                  if(iMpD0_2>k) then
                        write(*,*) 'not enough attr record for iMpD0_2'
                        stop 
                  endif 
                  if(iMpD0_3>k) then
                        write(*,*) 'not enough attr record for iMpD0_3'
                        stop 
                  endif 
                  if(n_iMnpara0 > 100) then
                        write(*,*) 'very large n_iMnpara0'
                        write(*,*) 'please change its hard-coded limit'
                        write(*,*) 'which is 100'
                        stop 
                  endif 
                  do j=1,n_iMnpara0     
                    if(iMnpara0(j)>k) then
                       write(*,*) 'not enough attr record for iMnpara0('
     &                      ,j,')'
                        stop 
                    endif
                  enddo
! read link records
                  read(9998,*) 
     &                      iMks0,iMl0
                  read(9998,*) n_iMlpara0
                  if( n_iMlpara0 >0) read(9998,*) iMlpara0(1:n_iMlpara0)

                  if(iMks0<=0) iMks0=FabricMeshL_numattMax
                  if(iMl0<=0) iMl0=FabricMeshL_numattMax
                  do j=1,n_iMlpara0     
                    if(iMlpara0(j)<0) iMlpara0(j)=FabricMeshL_numattMax
                  enddo
                  k=FabricMeshL_numattMax
                  if(iMm0>k) then
                   write(*,*) 'not enough link attr record for iMks0'
                        stop 
                  endif                 
                  if(iMl0>k) then
                   write(*,*) 'not enough link attr record for iMl0'
                        stop 
                  endif                 
                  if(n_iMlpara0 > 100) then
                        write(*,*) 'very large n_iMlpara0'
                        write(*,*) 'please change its hard-coded limit'
                        write(*,*) 'which is 100'
                        stop 
                  endif 
                  do j=1,n_iMlpara0     
                    if(iMlpara0(j)>k) then
                  write(*,*) 'not enough link attr record for iMlpara0('
     &                      ,j,')'
                        stop 
                    endif
                  enddo

                  iMeshmax=0
                  do k=1,numrecord
                     read(9998,*) tmpint(0:2)
     &                       ,FabricMesh_coord(in_fbc,k,1:3)
     &                       ,FabricMesh_attr(
     &                                in_fbc
     &                               ,k
     &                               ,1:FabricMesh_numatt(in_fbc))
                     iMeshmax=max(tmpint(1),iMeshmax)

                     if(tmpint(2) .gt. numMeshpointMax) then
                        write(*,*) 'not enough numMeshpointMax'
                        stop
                     endif

                     if(iMeshmax .gt. numMeshMax) then
                        write(*,*) 'not enough numMeshMax'
                        stop
                     endif
                     FabricMesh_npoint(in_fbc,tmpint(1))
     &                       =max(
     &                        FabricMesh_npoint(in_fbc,tmpint(1))
     &                       ,tmpint(2)
     &                           )
                     FabricMesh_address(in_fbc,tmpint(1),tmpint(2))=k
                  enddo
                  FabricMesh_nMesh(in_fbc)= iMeshmax  

                  iMeshmax=0
                  do k=1,numrecordL
                     read(9998,*) tmpint(0:2)
     &                           ,FabricMeshL_con(in_fbc,k,1:2)
     &                           ,FabricMeshL_attr(
     &                                    in_fbc
     &                                   ,k
     &                                   ,1:FabricMeshL_numatt(in_fbc))
                     iMeshmax=max(tmpint(1),iMeshmax)
                     if(iMeshmax .gt. numMeshMax) then
                        write(*,*) 'not enough numMeshMax'
                        stop
                     endif
                     FabricMesh_nline(in_fbc,tmpint(1))
     &                       =max(
     &                           FabricMesh_nline(in_fbc,tmpint(1))
     &                          ,tmpint(2)
     &                            )
                     if(tmpint(2) .gt. numMeshLinkMax) then
                        write(*,*) 'not enough numMeshLinkMax'
                        stop
                     endif

                     FabricMeshL_address(in_fbc,tmpint(1),tmpint(2))=k
                     FabricMeshL_addressRev(in_fbc,k,1:2)=tmpint(1:2)
                  enddo
                  close (9998)
               endif


!line records
               read(9997,*) FabricLineFlag_fbc(in_fbc)
               if(FabricLineFlag_fbc(in_fbc)) then
                  if(FabricLineGeneration .eq. 0) then
                  open(9998
     &                ,file='FabricMeshLineSize.txt'
     &                ,status='old',form='formatted' )
                  read(9998,*)!meshname
                  read(9998,*)!meshsize data
                  read(9998,*)!line name
                  read(9998,*)
     &                        numLinepointRecordMax
     &                       ,numLineMax
     &                       ,numLinepointMax
     &                       ,FabricLine_numattMax
     &                       ,FabricLine_numtargetMax
                  close (9998)
                  FabricMesh_numattMax=FabricMesh_numattMax+1
!     append zero column to end

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

                  FabricLineGeneration=1
                  FabricLine_Presence=.true.
                  endif


               read(9997,'(a)') infileFabricLine 
                  open(9998,
     &           file=infileFabricLine,status='old',form='formatted' )

                  read(9998,*) numrecord, FabricLine_numatt(in_fbc)

                  if(numrecord .gt. numLinepointRecordMax) then
                      write(*,*) 'not enough numLinepointRecordMax'
                      stop
                  endif
                  FabricLine_numrecord(in_fbc)=numrecord
                  if(FabricLine_numatt(in_fbc) 
     &                     .gt. FabricLine_numattMax) then
                      write(*,*) 'not enough FabricLine_numattMax'
                      stop
                  endif

                  read(9998,*) 
     &                      iLm0,iLc0,iLconR0,iLcon0,iLs0
     &                     ,iLks0,iLkb0,iLl0
     &                     ,iLpD0_3,iLpD0_1,iLpD0_2
                  read(9998,*) n_iLnpara0
                  if( n_iLnpara0 >0) read(9998,*) iLnpara0(1:n_iMnpara0)
                  if(iLm0<=0) iMm0=FabricLine_numattMax
                  if(iLc0<=0) iLc0=FabricLine_numattMax
                  if(iLconR0<=0) iLconR0=FabricLine_numattMax
                  if(iLcon0<=0) iLcon0=FabricLine_numattMax
                  if(iLs0<=0) iLs0=FabricLine_numattMax
                  if(iLks0<=0) iLks0=FabricLine_numattMax
                  if(iLkb0<=0) iLkb0=FabricLine_numattMax
                  if(iLl0<=0) iLl0=FabricLine_numattMax
                  if(iLpD0_1<=0) iLpD0_1=FabricLine_numattMax
                  if(iLpD0_2<=0) iLpD0_2=FabricLine_numattMax
                  if(iLpD0_3<=0) iLpD0_3=FabricLine_numattMax
                  do j=1,n_iLnpara0     
                    if(iLnpara0(j)<0) iLnpara0(j)=FabricLine_numattMax
                  enddo
                  k=FabricLine_numattMax
                  if(iLm0>k) then
                        write(*,*) 'not enough attr record for iLm0'
                        stop 
                  endif                 
                  if(iLc0>k) then
                        write(*,*) 'not enough attr record for iLc0'
                        stop 
                  endif                 
                  if(iLconR0>k) then
                        write(*,*) 'not enough attr record for iLconR0'
                        stop 
                  endif                 
                  if(iLcon0>k) then
                        write(*,*) 'not enough attr record for iLcon0'
                        stop 
                  endif                 
                  if(iLs0>k) then
                        write(*,*) 'not enough attr record for iLs0'
                        stop 
                  endif 
                  if(iLks0>k) then
                        write(*,*) 'not enough attr record for iLks0'
                        stop 
                  endif 
                  if(iLkb0>k) then
                        write(*,*) 'not enough attr record for iLkb0'
                        stop 
                  endif 
                  if(iLl0>k) then
                        write(*,*) 'not enough attr record for iLl0'
                        stop 
                  endif 
                  if(iLpD0_1>k) then
                        write(*,*) 'not enough attr record for iLpD0_1'
                        stop 
                  endif 
                  if(iLpD0_2>k) then
                        write(*,*) 'not enough attr record for iLpD0_2'
                        stop 
                  endif 
                  if(iLpD0_3>k) then
                        write(*,*) 'not enough attr record for iLpD0_3'
                        stop 
                  endif 
                  if(n_iMnpara0 > 100) then
                        write(*,*) 'very large n_iMnpara0'
                        write(*,*) 'please change its hard-coded limit'
                        write(*,*) 'which is 100'
                        stop 
                  endif 
                  do j=1,n_iLnpara0     
                    if(iLnpara0(j)>k) then
                       write(*,*) 'not enough attr record for iLnpara0('
     &                      ,j,')'

                        stop 
                    endif
                  enddo

                  FabricLine_npoint(in_fbc,1:numLineMax)=0.0
                  ilinemax=0
                  do k=1,numrecord
                     read(9998,*) tmpint(1:2)
     &                           ,FabricLine_coord(in_fbc,k,1:3)
     &                           ,FabricLine_attr(
     &                                    in_fbc
     &                                   ,k
     &                                   ,1:FabricLine_numatt(in_fbc))
                     ilinemax=max(tmpint(1),ilinemax)
                     if(ilinemax .gt. numlineMax) then
                        write(*,*) 'not enough numlineMax'
                        stop
                     endif
                     FabricLine_npoint(in_fbc,tmpint(1))
     &                       =max(
     &                            FabricLine_npoint(in_fbc,tmpint(1))
     &                           ,tmpint(2)
     &                            )
                     FabricLine_address(in_fbc,tmpint(1),tmpint(2))=k
                     FabricLine_addressRev(in_fbc,k,1:2)=tmpint(1:2)
                  enddo
                  FabricLine_nLine(in_fbc)= ilinemax  
                  close (9998)
               endif

               read(9997,*) frequency_fin(i+1,1:3)
               read(9997,*) i_fintype(i+1)
               read(9997,*) target_k_linktmp(1,1:i_fintype(i+1)) 
               read(9997,*) target_t_linktmp(1,1:i_fintype(i+1))  
 
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),1)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),1)   
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),2)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),2)    
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),3)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),3)    
 
               read(9997,*) flagfixedtmp(1,1:i_fintype(i+1)) 
               read(9997,*) flag_TDpt_Pnormal_Ext(i+1),Pnormal_Ext(i+1)
               read(9997,'(a)') infilepres
               if (flag_TDpt_Pnormal_Ext(i+1) .eq. 1) then
                  open(9998,
     &                 file=infilepres,status='old',form='formatted' )
                  read(9998,*) n_TDpt_Pnormal_Ext(i+1)
                  if (n_TDpt_Pnormal_Ext(i+1) .gt. 1000) then
                     write(*,*) 'size of press is larger than allowed'
                     write(*,*) 'For body', i+1
                     stop
                  endif
                  do j=1,n_TDpt_Pnormal_Ext(i+1)
                     read(9998,*) data_TDpt_Pnormal_Ext(i+1,j,1)
     &                           ,data_TDpt_Pnormal_Ext(i+1,j,2)
                  enddo
                  close (9998)
               endif



                read(9998,'(a)' ) namebody(i*10+1:(i2)*10)
                open(9998,file=infile,status='old',form='formatted' )
                read(9998,*) ns_ibm_r(i2),ns_ibm_i(i2), ns_ibm_rall(i2) 
     &                       ,numelr(i2), numeli(i2) 

               ns_ibm_r_fbc(in_fbc)=ns_ibm_r(i2)
               ns_ibm_rall_fbc(in_fbc)=ns_ibm_rall(i2)
               ns_ibm_i_fbc(in_fbc)=ns_ibm_i(i2)

               j2=0
               maxtmpx=0.0
               maxtmpy=0.0
               maxtmpz=0.0
               do j=1,ns_ibm_r_fbc(in_fbc)
                  j2=j2+1
                  read(9998,*) numi   
     &                    ,gx_ibm_fbc(in_fbc,j)   
     &                    ,gy_ibm_fbc(in_fbc,j)   
     &                    ,gz_ibm_fbc(in_fbc,j)   
     &                    ,target_ktype(1,j)
     &                    ,target_kvalue(1,j)
    
                  target_ktype_fbc(in_fbc,j)=target_ktype(1,j)
                  target_kvalue_fbc(in_fbc,j)=target_kvalue(1,j)


                 tmpcoord(1:3)=
     &           (/gx_ibm_fbc(in_fbc,j)
     &           , gy_ibm_fbc(in_fbc,j)
     &           , gz_ibm_fbc(in_fbc,j) /)

                 do jq=1,3
                     tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                 enddo
                 jq=1

                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)

                 gx_ibm_fbc(in_fbc,j)=tmpcoord(1)
                 gy_ibm_fbc(in_fbc,j)=tmpcoord(2)
                 gz_ibm_fbc(in_fbc,j)=tmpcoord(3)
                  maxtmpx=maxtmpx+abs(gx_ibm_fbc(in_fbc,j))
                  maxtmpy=maxtmpy+abs(gy_ibm_fbc(in_fbc,j))
                  maxtmpz=maxtmpz+abs(gz_ibm_fbc(in_fbc,j))
               end do
               do j=1,ns_ibm_i_fbc(in_fbc)
                      j2=j2+1
                      read(9998,*) numi,   
     &                    gx_ibm_fbc(in_fbc,j2),   
     &                    gy_ibm_fbc(in_fbc,j2),   
     &                    gz_ibm_fbc(in_fbc,j2), 
     &                    bndycnd_fbc(in_fbc,j2),
     &                    boundrylist_fbc(in_fbc,j,1),   
     &                    boundrylist_fbc(in_fbc,j,2),   
     &                    boundrylist_fbc(in_fbc,j,3)

                 tmpcoord(1:3)=
     &           (/gx_ibm_fbc(in_fbc,j2)
     &           , gy_ibm_fbc(in_fbc,j2)
     &           , gz_ibm_fbc(in_fbc,j2) /)

                 do jq=1,3
                     tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                 enddo
                 jq=1
               
                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)
                 gx_ibm_fbc(in_fbc,j2)=tmpcoord(1)
                 gy_ibm_fbc(in_fbc,j2)=tmpcoord(2)
                 gz_ibm_fbc(in_fbc,j2)=tmpcoord(3)
               enddo

               if(FabricMeshFlag_fbc(in_fbc)) then
                  do j=1,FabricMesh_numrecord(in_fbc)
                     tmpcoord(1:3)=FabricMesh_coord(in_fbc,j,1:3)
                     do jq=1,3
                        tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                     enddo
                     jq=1
               
                     tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                     tmpcoord(1:3)=
     &                 tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &                +tmpcoordR(1:3)
                     FabricMesh_coord(in_fbc,j,1:3)=tmpcoord(1:3)
                  enddo
               endif
               if(FabricLineFlag_fbc(in_fbc)) then
                  do j=1,FabricLine_numrecord(in_fbc)
                     tmpcoord(1:3)=FabricLine_coord(in_fbc,j,1:3)
                     do jq=1,3
                        tmpcoord(jq)=tmpcoordscale(jq)
     &                           +(tmpcoord(jq)-tmpcoordscale(jq))
     &                            *tmpcoordscale(3+jq)
                     enddo
                     jq=1
               
                     tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                     tmpcoord(1:3)=
     &                 tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &                +tmpcoordR(1:3)
                     FabricLine_coord(in_fbc,j,1:3)=tmpcoord(1:3)
                  enddo
               endif

                i_fix=0             
               do j=1,ns_ibm_r_fbc(i2)
                  if(target_ktype(1,j) .ne. 0) then
                    i_fix=i_fix+1
                    target_point_num_fbc(in_fbc,i_fix)=j
                    target_k_link_fbc(in_fbc,i_fix)=
     &                            target_k_linktmp(1,target_ktype(1,j))
     &                           *target_kvalue(1,j)
                    target_t_link_fbc(in_fbc,i_fix)=
     &                            target_t_linktmp(1,target_ktype(1,j))
                    a_fin_fbc(in_fbc,i_fix,1:3)=
     &                            a_fintmp(1,target_ktype(1,j),1:3)
                    phi_fin_fbc(in_fbc,i_fix,1:3)=
     &                            phi_fintmp(1,target_ktype(1,j),1:3)
                    flagfixed_fbc(in_fbc,i_fix)=
     &                            flagfixedtmp(1,target_ktype(1,j))
                  end if     
               end do
               target_num(i+1)=i_fix   

               if(FabricMeshFlag_fbc(in_fbc)) then
                  i_fix=0 
                  do j=1,FabricMesh_numrecord(in_fbc)
                     jq=int(FabricMesh_attr(i,j,iMconR0)+1d-4)
                     if(jq .ne. 0) then
                       i_fix=i_fix+1
                       FabricMesh_target_point_num(in_fbc,i_fix)=j
                       FabricMesh_target_k_link(in_fbc,i_fix)=
     &                            target_k_linktmp(1,jq)
     &                           *FabricMesh_attr(i,j,iMconR0)
                       FabricMesh_target_t_link(in_fbc,i_fix)=
     &                            target_t_linktmp(1,jq)
                       FabricMesh_a_fin(in_fbc,i_fix,1:3)=
     &                            a_fintmp(1,jq,1:3)
                       FabricMesh_phi_fin(in_fbc,i_fix,1:3)=
     &                            phi_fintmp(1,jq,1:3)
                       FabricMesh_flagfixed(in_fbc,i_fix)=
     &                            flagfixedtmp(1,jq)
                    end if 
                  enddo 
                  FabricMesh_target_num(in_fbc)=i_fix 
               endif
               if(FabricLineFlag_fbc(in_fbc)) then
                  i_fix=0 
                  do j=1,FabricLine_numrecord(in_fbc)
                     jq=int(FabricLine_attr(i,j,iLconR0)+1d-4)
                     if(jq .ne. 0) then
                       i_fix=i_fix+1
                       FabricLine_target_point_num(in_fbc,i_fix)=j
                       FabricLine_target_k_link(in_fbc,i_fix)=
     &                            target_k_linktmp(1,jq)
     &                           *FabricLine_attr(i,j,iLconR0)
                       FabricLine_target_t_link(in_fbc,i_fix)=
     &                            target_t_linktmp(1,jq)
                       FabricLine_a_fin(in_fbc,i_fix,1:3)=
     &                            a_fintmp(1,jq,1:3)
                       FabricLine_phi_fin(in_fbc,i_fix,1:3)=
     &                            phi_fintmp(1,jq,1:3)
                       FabricLine_flagfixed(in_fbc,i_fix)=
     &                            flagfixedtmp(1,jq)
                    end if 
                  enddo 
                  FabricLine_target_num(in_fbc)=i_fix
                endif

          
               do j=1,numelr(i2)+numeli(i2)
                  read(9998,*) nea_fbc(in_fbc,j,1:3)
               end do

               coefmaxtmpx=0.0
               if(maxtmpx .lt. 1d-6) then 
                  coefmaxtmpx=1.0d-6
               endif

               coefmaxtmpy=0.0
               if(maxtmpy .lt. 1d-6) then 
                  coefmaxtmpy=1.0d-6
               endif

               coefmaxtmpz=0.0
               if(maxtmpz .lt. 1d-6) then 
                  coefmaxtmpz=1.0d-6
               endif

               do j=1,ns_ibm_rall_fbc(in_fbc)
                 gx_ibm_fbc(in_fbc,j)=gx_ibm_fbc(in_fbc,j)+ coefmaxtmpx
                 gy_ibm_fbc(in_fbc,j)=gy_ibm_fbc(in_fbc,j)+ coefmaxtmpy
                 gz_ibm_fbc(in_fbc,j)=gz_ibm_fbc(in_fbc,j)+ coefmaxtmpz

                 gx_bp_fbc(in_fbc,j)=gx_ibm_fbc(in_fbc,j)
                 gy_bp_fbc(in_fbc,j)=gy_ibm_fbc(in_fbc,j)
                 gz_bp_fbc(in_fbc,j)=gz_ibm_fbc(in_fbc,j)
                 gx_bp_fbc0(in_fbc,j)=gx_ibm_fbc(in_fbc,j)
                 gy_bp_fbc0(in_fbc,j)=gy_ibm_fbc(in_fbc,j)
                 gz_bp_fbc0(in_fbc,j)=gz_ibm_fbc(in_fbc,j)

                 gx_ibmpre_fbc(in_fbc,j)=gx_ibm_fbc(in_fbc,j)
                 gy_ibmpre_fbc(in_fbc,j)=gy_ibm_fbc(in_fbc,j)
                 gz_ibmpre_fbc(in_fbc,j)=gz_ibm_fbc(in_fbc,j)
               end do

               if(FabricMeshFlag_fbc(i)) then
                do j=1,FabricMesh_numrecord(i)
                  FabricMesh_coordbp(i,j,1:3)= 
     &              FabricMesh_coord(i,j,1:3) 
                  FabricMesh_coordbp0(i,j,1:3)= 
     &              FabricMesh_coord(i,j,1:3) 
                  FabricMesh_coordpre(i,j,1:3)= 
     &              FabricMesh_coord(i,j,1:3) 
                enddo
               endif
               if(FabricLineFlag_fbc(i)) then
                do j=1,FabricLine_numrecord(i)
                  FabricLine_coordbp(i,j,1:3)= 
     &              FabricLine_coord(i,j,1:3) 
                  FabricLine_coordbp0(i,j,1:3)= 
     &              FabricLine_coord(i,j,1:3) 
                  FabricLine_coordpre(i,j,1:3)= 
     &              FabricLine_coord(i,j,1:3) 
                enddo
               endif

               close (9998)
!some pre processing
! form connectivities
               !two vector 
               !conn_face_edge_fbc(nr_ibm_fbc,elem,3)  INT
               !conn_face_face_fbc(nr_ibm_fbc,elem,3)  INT
               !faces_vertisin_fbc(nr_ibm_fbc,ns_ibm_fbc,50)  INT

               !faces_nvertisin_fbc(nr_ibm_fbc,ns_ibm_fbc)
               faces_nvertisin_fbc=0          
               do j=1,numelr(i2)
                  do j2=1,3
                     jq=nea_fbc(in_fbc,j,j2)
                     faces_nvertisin_fbc(in_fbc,jq)
     &                            =faces_nvertisin_fbc(in_fbc,jq)+1
                     faces_vertisin_fbc(in_fbc
     &                                 ,jq
     &                                 ,faces_nvertisin_fbc(in_fbc,jq)
     &                                  )=j
                  enddo
               enddo
               ie2=-1
               do j=1,numelr(i2)
                  do j2=1,3
                     ifound=0
                     iv0=nea_fbc(in_fbc,j,mod(j2-1,3)+1)
                     iv1=nea_fbc(in_fbc,j,mod(j2  ,3)+1)                  
                     do k=1,faces_nvertisin_fbc(in_fbc,iv0)
                        ie0=faces_vertisin_fbc(in_fbc,iv0,k)
                        if(ie0 .ne. j) then
                          ifound=0
                          do jq=1,faces_nvertisin_fbc(in_fbc,iv1)
                             ie1=faces_vertisin_fbc(in_fbc,iv1,jq)
                             if(ie2 .ne. j)   ifound=1       
                          enddo
                        endif
                        if(ifound==1) ie2=j
                     enddo
                     if(ie2==-1) then
                        write(*,*) 'Error: Found boundary in real part'
                        write(*,*) 'FBC',j,j2,ie2
                        STOP
                     else
                        jq=-1
                        do k=1,3
                           if((iv0+iv1) .eq. 
     &                      (nea_fbc(in_fbc,ie2,mod(k-1,3)+1)
     &                      +nea_fbc(in_fbc,ie2,mod(k  ,3)+1))) then
                               jq=k
                               exit
                           endif
                        enddo
                        if(jq ==-1) then
                         write(*,*) 'Error: sth worng if coonectivity'
                         write(*,*) 'FBC',j,j2,ie2
                         STOP 
                        endif
                     endif                       
                     conn_face_edge_fbc(in_fbc,j,j2)=jq
                     conn_face_face_fbc(in_fbc,j,j2)=ie2
                  enddo
               enddo
               do j=1,numelr(i2)
                  do j2=1,3
                     iv1=conn_face_face_fbc(in_fbc,j,j2)   
                     iv0=conn_face_edge_fbc(in_fbc,j,j2) 
                     if(j .ne. conn_face_face_fbc(in_fbc,iv1,iv0)) then
                       write(*,*)'Error: conn info: face mismatch'  
                       write(*,*) 'FBC',j,j2
                       STOP
                     endif
                     if(j2 .ne. conn_face_edge_fbc(in_fbc,iv1,iv0)) then
                       write(*,*)'Error: conn info: edge mismatch'  
                       write(*,*) 'FBC',j,j2   
                       STOP
                     endif
                     jq=nea_fbc(in_fbc,j,mod(j2-1,3)+1)
                     k  =nea_fbc(in_fbc,j,mod(j2  ,3)+1) 
                     ie0=nea_fbc(in_fbc,j,mod(j2-1,3)+1)    
                     ie1=nea_fbc(in_fbc,j,mod(j2-1,3)+1)            
                     if(dabs(gx_ibm_fbc(in_fbc,jq)
     &                      +gx_ibm_fbc(in_fbc,k)
     &                      -gx_ibm_fbc(in_fbc,ie0)
     &                      -gx_ibm_fbc(in_fbc,ie1)) .ge. 1d-4) then
                       write(*,*)'Error: conn info: edge mismatch'  
                       write(*,*) 'FBC',j,j2   
                       write(*,*) jq,k,ie0,ie1
                       STOP
                     endif
                  enddo
                enddo
! generate edge list
               iedgelist=0
               do j=1,numelr(i2)             
                  do j2=1,3
                     ie1=conn_face_face_fbc(in_fbc,j,j2)
                     if(j .lt. ie1) then
                     iv0=nea_fbc(in_fbc,j,mod(j2-1,3)+1)
                     iv1=nea_fbc(in_fbc,j,mod(j2  ,3)+1)
                     iedgelist=iedgelist+1
                     EdgeList_fbc(in_fbc,iedgelist,1:2)=(/iv0,iv1/)
                     endif
                  enddo
               enddo
               numedge(in_fbc)=iedgelist
               maxlength_fbc(in_fbc)=-1.0
               avlength_fbc(in_fbc)=0
               do j=1,numedge(in_fbc)
                  iv0=EdgeList_fbc(in_fbc,j,1) 
                  iv1=EdgeList_fbc(in_fbc,j,2)
                  maxtmpx= sqrt(
     &             (gx_ibm_fbc(in_fbc,iv0)-gx_ibm_fbc(in_fbc,iv1))**2   
     &            +(gy_ibm_fbc(in_fbc,iv0)-gy_ibm_fbc(in_fbc,iv1))**2   
     &            +(gz_ibm_fbc(in_fbc,iv0)-gz_ibm_fbc(in_fbc,iv1))**2)  
                  EdgeLength0_fbc(in_fbc,j)=maxtmpx                     
                  avlength_fbc(in_fbc)=avlength_fbc(in_fbc)+ maxtmpx
                  maxlength_fbc(in_fbc)=
     &                    max(maxlength_fbc(in_fbc),maxtmpx)
               enddo
                if(theboss)  then
                   write(*,*) 'For Fabric Body ',i2
                   write(*,*) 'There are ',numedge(in_fbc),' edges'
                   write(*,*) 'with average legth of '
     &                         ,avlength_fbc(in_fbc)
                   write(*,*) 'and maximum legth of '
     &                         , maxlength_fbc(in_fbc)    
                endif
! genrate face pair
               ifacepair=0
               do j=1,numelr(i2)             
                  do j2=1,3
                     ie1=conn_face_face_fbc(in_fbc,j,j2)
                     if(j .lt. ie1) then
                        jq=mod(j2+1,3)  !face1 nodes
                        iv0=nea_fbc(in_fbc,j,mod(jq-1,3)+1)
                        iv1=nea_fbc(in_fbc,j,mod(jq-0,3)+1)   
                        iv2=nea_fbc(in_fbc,j,mod(jq+1,3)+1) 

                        k=conn_face_edge_fbc(in_fbc,j,j2)
                        iv3=nea_fbc(in_fbc,ie1,mod(k+1,3)+1)
                        ifacepair= ifacepair+1
                        facepairNode_fbc(in_fbc,ifacepair,1:4)=
     &                                  (/iv0,iv1,iv2,iv3/)
                        facepairEle_fbc(in_fbc,ifacepair,1:2)=
     &                                  (/j,ie1/)
                     endif
                  enddo
               enddo
               numeFacePair(i2)=ifacepair
!=================== Fibrous Shell ====================================
            elseif (bodytype(i+1) .eq. fibrousshell) then
!     fibrous shell input read

                i2=i+1

                if(theboss) 
     &              write(*,*) '====Body ',i2,'===='

                read(9997,'(a)') infile

                if(theboss) 
     &              write(*,*) 'name of the inputfile=',infile

                read(9997,*) MaterialTypeIBM(i2),BendingZero(i2)   
     &                      ,i_phi_mem_coefFLG,i_phi_mem_coef
                if(theboss) then
                  write(*,*) 'Bending energy type ',MaterialTypeIBM(i2)
                  write(*,*) 'BendingZero='
                  write(*,*) 'Have non-zero intial Bending '   
     &                      ,BendingZero(i2)
                  write(*,*) 'i_phi_mem_coefFLG= ' 
                  write(*,*) 'read from fsh inp file~0; multiplier>0; '   
     &                      ,i_phi_mem_coefFLG
                endif

                read(9997,*) e4coef(i2)   
     &                      ,phi_mem_coef(i2,1:i_phi_mem_coef)   

                if(theboss) 
     &                     write(*,*) phi_mem_coef(i2,1:3)   

!====check for consistency == START
                if(e4coef(i+1) .ge. 0.0) then
                if(theboss) write(*,*) 'Option of e4coef >=0.0 is for',
     &          'elastic mem. with near orthogonal curvilinear grid '
                 temp3=e4coef(i+1)
                 if((e4coef(i+1).gt. 50.0).and.(e4coef(i+1).lt. 100.0)) 
     &            temp3=e4coef(i+1)-50.0
                phi_mem_coef(i2,3) =phi_mem_coef(i2,3)*(1.0-temp3)/2.0
               endif
!====check for consistency == END

                read(9997,*) phi_ben_coef(i2),phi_ben2_coef   
     &                      ,phi_nus_coef(i2)
                if(theboss) 
     &                     write(*,*) phi_ben_coef(i2),phi_ben2_coef   
     &                      ,phi_nus_coef(i2)

                read(9997,*) Piezo_Formulation_Flag(i2)
     &                      ,Piezo_Dist_Flag(i2)
     &                      ,piezo_beta(i2)
     &                      ,piezo_gamma(i2)
     &                      ,piezo_damp(i2)
     &                      ,piezo_Coef(i2)
     &                      ,piezo_alpha11(i2)   
     &                      ,piezo_alpha22(i2)
     &                      ,piezo_alpha12(i2) 
                read(9997,*) rampshear
                read(9997,*)
     &                      istartS(1),iraiseS(1)
     &                     ,istartS(2),iraiseS(2)


! read BCs  (0:fixed, 1:free, 2:pinned, 3:periodic)
                read(9997,*) BoundaryFibShell_fsh(in_fsh,1,1)
!     BC for smin   
     &               ,BoundaryFibShell_fsh(in_fsh,1,2)
!     BC for smax   
     &               ,BoundaryFibShell_fsh(in_fsh,2,1)
!     BC for qmin   
     &               ,BoundaryFibShell_fsh(in_fsh,2,2)
!     BC for qmax

                if(theboss) 
     &               write(*,*)   BoundaryFibShell_fsh(in_fsh,1,1)
!     BC for smin   
     &               ,BoundaryFibShell_fsh(in_fsh,1,2)
!     BC for smax   
     &               ,BoundaryFibShell_fsh(in_fsh,2,1)
!     BC for qmin   
     &               ,BoundaryFibShell_fsh(in_fsh,2,2)
!     BC for qmax

               fix_bc1(1)=.true.  !,1:j  ,2:jq 
               fix_bc1(2)=.true.
               fix_bc2(1)=.true.
               fix_bc2(2)=.true.
 
               if(BoundaryFibShell_fsh(in_fsh,1,1).eq. 1)
     &               fix_bc1(1)=.false.
               if(BoundaryFibShell_fsh(in_fsh,1,2).eq. 1)
     &               fix_bc1(2)=.false.
               if(BoundaryFibShell_fsh(in_fsh,2,1).eq. 1)
     &               fix_bc2(1)=.false.
               if(BoundaryFibShell_fsh(in_fsh,2,2).eq. 1)
     &               fix_bc2(2)=.false.


! correct for potential periodic boundary conditions
               if (BoundaryFibShell_fsh(in_fsh,1,1).eq.3)   
     &             BoundaryFibShell_fsh(in_fsh,1,2)=3
               if (BoundaryFibShell_fsh(in_fsh,1,2).eq.3)   
     &             BoundaryFibShell_fsh(in_fsh,1,1)=3
               if (BoundaryFibShell_fsh(in_fsh,2,1).eq.3)   
     &             BoundaryFibShell_fsh(in_fsh,2,2)=3
               if (BoundaryFibShell_fsh(in_fsh,2,2).eq.3)   
     &             BoundaryFibShell_fsh(in_fsh,2,1)=3


               read(9997,*) Genalpha_timesolver(i2)
     &                     ,Genalpha_niter(i2)
     &                     ,propDamp(i2,1:2)
               if(Genalpha_timesolver(i2) .gt. 5) then
                  write(*,*) 'Genalpha_timesolver for', i2,'is >5?'
                  stop
               endif
               read(9997,*) rhoinf_Coef,method_rhoinf_Coef
     &                     ,alpham_Coef(i2),alphaf_Coef(i2)

               if(method_rhoinf_Coef .eq. 2) then
                  !rhoinf_Coef=1.0/rhoinf_Coef
                  alpham_coef(i2)=(2.0*rhoinf_Coef-1)/(1.0+rhoinf_Coef)
!     2nd order system
                  alphaf_coef(i2)=rhoinf_Coef/(1.0+rhoinf_Coef)
               elseif(method_rhoinf_Coef .eq. 1) then
                  alpham_coef(i2)=
     &                 0.5*(3.0*rhoinf_Coef-1)/(1.0+rhoinf_Coef)
!     1st order system
                  alphaf_coef(i2)=rhoinf_Coef/(1.0+rhoinf_Coef)
               endif
               beta_Coef(i2)=0.25*
     &                        (1.0-alpham_coef(i2)+alphaf_coef(i2))**2
               gamma_Coef(i2)=0.5-alpham_coef(i2)+alphaf_coef(i2)



               read(9997,*) frequency_fin(i+1,1:3)
               read(9997,*) i_fintype(i+1)
               read(9997,*) target_k_linktmp(1,1:i_fintype(i+1)) 
               read(9997,*) target_t_linktmp(1,1:i_fintype(i+1))   
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),1)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),1)    
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),2)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),2)   
               read(9997,*) a_fintmp(1,1:i_fintype(i+1),3)    
               read(9997,*) phi_fintmp(1,1:i_fintype(i+1),3)   
               read(9997,*) flagfixedtmp(1,1:i_fintype(i+1)) 
               read(9997,*) flag_TDpt_Pnormal_Ext(i+1),Pnormal_Ext(i+1)
               read(9997,'(a)') infilepres
               if (flag_TDpt_Pnormal_Ext(i+1) .eq. 1) then
                  open(9998,
     &                 file=infilepres,status='old',form='formatted' )
                  read(9998,*) n_TDpt_Pnormal_Ext(i+1)
                  if (n_TDpt_Pnormal_Ext(i+1) .gt. 1000) then
                     write(*,*) 'size of press is larger than allowed'
                     write(*,*) 'For body', i+1
                     stop
                  endif
                  do j=1,n_TDpt_Pnormal_Ext(i+1)
                     read(9998,*) data_TDpt_Pnormal_Ext(i+1,j,1)
     &                           ,data_TDpt_Pnormal_Ext(i+1,j,2)
                  enddo
                  close (9998)
               endif

               if(theboss) then
               write(*,*) frequency_fin(i+1,1:3)
               write(*,*) i_fintype(i+1)
               write(*,*) target_k_linktmp(1,1:i_fintype(i+1)) 
               write(*,*) target_t_linktmp(1,1:i_fintype(i+1))   
               write(*,*) a_fintmp(1,1:i_fintype(i+1),1)    
               write(*,*) phi_fintmp(1,1:i_fintype(i+1),1)    
               write(*,*) a_fintmp(1,1:i_fintype(i+1),2)    
               write(*,*) phi_fintmp(1,1:i_fintype(i+1),2)   
               write(*,*) a_fintmp(1,1:i_fintype(i+1),3)    
               write(*,*) phi_fintmp(1,1:i_fintype(i+1),3)   
               write(*,*) flagfixedtmp(1,1:i_fintype(i+1)) 
               write(*,*) flag_TDpt_Pnormal_Ext(i+1),Pnormal_Ext(i+1)
               endif


               open(9998,file=infile,status='old',form='formatted' )
               read(9998,'(a)' ) namebody(i*10+1:(i2)*10)
               read(9998,*)   
     &                    nq_ibm_r(i2),ns_ibm_r(i2), masterbdy(i2)   
     &                   ,nMem_Coef_fsh(in_fsh)
                if(i_phi_mem_coef .ne. nMem_Coef_fsh(in_fsh)) then
                   write(*,*)    
     &   'number of membrane coefs are not correct in inputdist file'
                   stop
                endif
                ns_ibm_rall(i2)=ns_ibm_r(i2)
                ns_ibm_i(i2)=0

                ns_ibm_r_fsh(in_fsh)=ns_ibm_r(i2)
                nq_ibm_r_fsh(in_fsh)=nq_ibm_r(i2)
                ns_ibm_rall_fsh(in_fsh)=ns_ibm_rall(i2)
                ns_ibm_i_fsh(in_fsh)=ns_ibm_i(i2)

               do jq=1,nq_ibm_r_fsh(in_fsh)
               do j =1,ns_ibm_r_fsh(in_fsh)
                if (i_phi_mem_coefFLG .ne. 0) then
                  read(9998,*) numi,numj   
     &                    ,gx_ibm_fsh(in_fsh,jq,j)  
     &                    ,gy_ibm_fsh(in_fsh,jq,j)  
     &                    ,gz_ibm_fsh(in_fsh,jq,j)   
     &                    ,target_ktype(jq,j)
     &                    ,Mem_Coef_fsh(in_fsh,jq,j,
     &                          1:nMem_Coef_fsh(in_fsh))
     &                    ,Ben_Coef_fsh(in_fsh,jq,j,1)
     &                    ,Ben_Coef_fsh(in_fsh,jq,j,2)
     &                    ,Ben_Coef_fsh(in_fsh,jq,j,3)

                else
                  read(9998,*) numi,numj   
     &                    ,gx_ibm_fsh(in_fsh,jq,j)  
     &                    ,gy_ibm_fsh(in_fsh,jq,j)  
     &                    ,gz_ibm_fsh(in_fsh,jq,j)   
     &                    ,target_ktype(jq,j)
                  Mem_Coef_fsh(in_fsh,jq,j,1:nMem_Coef_fsh(in_fsh))=1.d0

                  Ben_Coef_fsh(in_fsh,jq,j,1)=1.d0
                  Ben_Coef_fsh(in_fsh,jq,j,2)=1.d0
                  Ben_Coef_fsh(in_fsh,jq,j,3)=1.d0
                endif
                if (i_phi_mem_coefFLG .ge. 0) then
                  do k=1,nMem_Coef_fsh(in_fsh)
                    Mem_Coef_fsh(in_fsh,jq,j,k)=   
     &                Mem_Coef_fsh(in_fsh,jq,j,k)*phi_mem_coef(i2,k)
                  enddo
                  Ben_Coef_fsh(in_fsh,jq,j,1)=   
     &              Ben_Coef_fsh(in_fsh,jq,j,1)*phi_ben_coef(i2)
                  Ben_Coef_fsh(in_fsh,jq,j,2)=   
     &              Ben_Coef_fsh(in_fsh,jq,j,2)*phi_ben2_coef
                  Ben_Coef_fsh(in_fsh,jq,j,3)=   
     &              Ben_Coef_fsh(in_fsh,jq,j,3)*phi_nus_coef(i2)
                endif
                if(istartS(1) .gt. 0) then
                  coef_j= S_Alpha(j,ns_ibm_r_fsh(in_fsh)
     &                 ,fix_bc1(1),fix_bc2(1),istartS(1),iraiseS(1))
                else
                  coef_j=1.0
                endif
                if(istartS(2) .gt. 0) then
                  coef_jq= S_Alpha(jq,nq_ibm_r_fsh(in_fsh)
     &                  ,fix_bc1(2),fix_bc2(2),istartS(2),iraiseS(2))
                else
                  coef_jq=1.0
                endif
                 if(rampshear) then
                 Mem_Coef_fsh(in_fsh,jq,j,3)= 
     &                      coef_j*coef_jq*Mem_Coef_fsh(in_fsh,jq,j,3)
                 endif

                 piezo_Coef_fsh(in_fsh,jq,j,-1)= 
     &                      piezo_Coef(i2)
                 piezo_Coef_fsh(in_fsh,jq,j,0)= 
     &                      piezo_damp(i2)
                 piezo_Coef_fsh(in_fsh,jq,j,1)= 
     &                      piezo_beta(i2)
                 piezo_Coef_fsh(in_fsh,jq,j,2)= 
     &                      piezo_gamma(i2)
                 piezo_Coef_fsh(in_fsh,jq,j,3)= 
     &                      coef_j*coef_jq*piezo_alpha11(i2)
                 piezo_Coef_fsh(in_fsh,jq,j,4)= 
     &                      coef_j*coef_jq*piezo_alpha22(i2)
                 piezo_Coef_fsh(in_fsh,jq,j,5)= 
     &                      coef_j*coef_jq*piezo_alpha12(i2)

                 tmpcoord(1:3)=
     &           (/gx_ibm_fsh(in_fsh,jq,j)
     &           , gy_ibm_fsh(in_fsh,jq,j)
     &           , gz_ibm_fsh(in_fsh,jq,j) /)

                 do k=1,3
                     tmpcoord(k)=tmpcoordscale(k)
     &                           +(tmpcoord(k)-tmpcoordscale(k))
     &                            *tmpcoordscale(3+k)
                 enddo

                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)
                 gx_ibm_fsh(in_fsh,jq,j)=tmpcoord(1)
                 gy_ibm_fsh(in_fsh,jq,j)=tmpcoord(2)
                 gz_ibm_fsh(in_fsh,jq,j)=tmpcoord(3)

               end do
               end do              
               ns_ibm_rall_fsh(in_fsh)=ns_ibm_r_fsh(in_fsh)

               if (masterbdy(i2) .gt. 0) then
                  do jq=1,nq_ibm_r_fsh(in_fsh)
                  do j=1,ns_ibm_r_fsh(in_fsh)
                      read(9998,*) numi, npos_fsh(in_fsh,jq,j) 
                  end do
                  end do

               end if
                i_fix=0   
               do jq=1,nq_ibm_r_fsh(in_fsh)          
               do j=1,ns_ibm_r_fsh(in_fsh)
                  if(target_ktype(jq,j) .ne. 0) then
                    i_fix=i_fix+1
                    target_point_num_fsh(in_fsh,1,i_fix)=j
                    target_point_num_fsh(in_fsh,2,i_fix)=jq
                    target_k_link_fsh(in_fsh,i_fix)=
     &                            target_k_linktmp(1,target_ktype(jq,j))
                    target_t_link_fsh(in_fsh,i_fix)=
     &                            target_t_linktmp(1,target_ktype(jq,j))
                    a_fin_fsh(in_fsh,i_fix,1:3)=
     &                            a_fintmp(1,target_ktype(jq,j),1:3)
                    phi_fin_fsh(in_fsh,i_fix,1:3)=
     &                            phi_fintmp(1,target_ktype(jq,j),1:3)
                    flagfixed_fsh(in_fsh,i_fix)=
     &                            flagfixedtmp(1,target_ktype(jq,j))
                  end if     
               end do
               end do
               target_num(i+1)=i_fix             
               i2=i+1
               do jq=1,nq_ibm_r_fsh(in_fsh) 
               do j=1,ns_ibm_rall_fsh(in_fsh)
                 gx_ibm_fsh(in_fsh,jq,j)=gx_ibm_fsh(in_fsh,jq,j)
                 gy_ibm_fsh(in_fsh,jq,j)=gy_ibm_fsh(in_fsh,jq,j)
                 gz_ibm_fsh(in_fsh,jq,j)=gz_ibm_fsh(in_fsh,jq,j)

                 gx_bp_fsh(in_fsh,jq,j)=gx_ibm_fsh(in_fsh,jq,j)
                 gy_bp_fsh(in_fsh,jq,j)=gy_ibm_fsh(in_fsh,jq,j)
                 gz_bp_fsh(in_fsh,jq,j)=gz_ibm_fsh(in_fsh,jq,j)
                 gx_bp_fsh0(in_fsh,jq,j)=gx_ibm_fsh(in_fsh,jq,j)
                 gy_bp_fsh0(in_fsh,jq,j)=gy_ibm_fsh(in_fsh,jq,j)
                 gz_bp_fsh0(in_fsh,jq,j)=gz_ibm_fsh(in_fsh,jq,j)

                 gx_ibmpre_fsh(in_fsh,jq,j)=gx_ibm_fsh(in_fsh,jq,j)
                 gy_ibmpre_fsh(in_fsh,jq,j)=gy_ibm_fsh(in_fsh,jq,j)
                 gz_ibmpre_fsh(in_fsh,jq,j)=gz_ibm_fsh(in_fsh,jq,j)
                  gx_ibm_massive_fsh(in_fsh,jq,j)=
     &                            gx_ibm_fsh(in_fsh,jq,j)
                  gy_ibm_massive_fsh(in_fsh,jq,j)=
     &                            gy_ibm_fsh(in_fsh,jq,j)
                  gz_ibm_massive_fsh(in_fsh,jq,j)=
     &                            gz_ibm_fsh(in_fsh,jq,j)

                  vPiezo_fsh(in_fsh,jq,j)=0.0
                  vPiezoold_fsh(in_fsh,jq,j)=0.0

                end do
               end do
                close (9998)
          else
               stop 'this option of structure is not available!'              
          end if

!++++++++++++++++++++ Rigid Body motions of the ith body : START ++++++++++++++++++++++++++
               read(9997,*) 
               read(9997,*) Gt0(i+1,1:3)
               read(9997,*) Tramps(i+1,1:3)
               read(9997,*) Trampx(i+1,1:3)
               read(9997,*) ampx(i+1,1:3)
               read(9997,*) freqx(i+1,1:3)
               read(9997,*) phix(i+1,1:3)
               read(9997,*) Trampt(i+1,1:3)
               read(9997,*) ampt(i+1,1:3),rot_seq(i+1)
               read(9997,*) freqt(i+1,1:3)
               read(9997,*) phit(i+1,1:3)


               ampt(i+1,1:3)=ampt(i+1,1:3)/180.0d0*pi1
               phix(i+1,1:3)=phix(i+1,1:3)/180.0d0*pi1
               phit(i+1,1:3)=phit(i+1,1:3)/180.0d0*pi1    
           
                 tmpcoord(1:3)=
     &           (/Gt0(i+1,1)
     &           , Gt0(i+1,2)
     &           , Gt0(i+1,3) /)

                 do k=1,3
                     tmpcoord(k)=tmpcoordscale(k)
     &                           +(tmpcoord(k)-tmpcoordscale(k))
     &                            *tmpcoordscale(3+k)
                 enddo


                 tmpcoord(1:3)=tmpcoord(1:3)-tmpcoordR(1:3)
                 tmpcoord(1:3)=
     &                tmpcoordshift(1:3)+matmul(tmpRmat,tmpcoord(1:3))
     &               +tmpcoordR(1:3)

                 Gt0(i+1,1)=tmpcoord(1)
                 Gt0(i+1,2)=tmpcoord(2)
                 Gt0(i+1,3)=tmpcoord(3)

      enddo !end body numbers i

               read(9997,*) 
               read(9997,*) p_inflation,dsecy_IBM,dsecz_IBM 

               if(theboss) print*,  'done with inputdist.dat '

               dsecy_IBM=1.0
               dsecz_IBM=0.0

      close(9997)  


         call shapeinitialize(jtypeall
     *   ,nn,dnndv,dnndw,ddnndvdv,ddnndwdw,ddnndvdw
     *   ,nnCont,dnndvCont,dnndwCont,
     *    ddnndvdvCont,ddnndwdwCont,ddnndvdwCont
     *   ,nn0,dnndv0,dnndw0,ddnndvdv0,ddnndwdw0,ddnndvdw0
     *   ,gau,ngau,wei
     *   ,gauCont,ngauCont,weiCont
     *   ,gau0   ,ngau0   ,wei0
     *   ,gau_h   ,ngau_h   ,wei_h
     *   ,n_jtypeall)

      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
