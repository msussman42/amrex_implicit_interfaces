c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!@c
      subroutine tick(rtime,dt,nt,monitorON,ntsave
     &     ,uibm1in_fib,uibm2in_fib,uibm3in_fib
     &     ,uibm1in_fsh,uibm2in_fsh,uibm3in_fsh
     &     ,uibm1in_esh,uibm2in_esh,uibm3in_esh
     &     ,uibm1in_fbc,uibm2in_fbc,uibm3in_fbc
     &     ,fibm1out_fib,fibm2out_fib,fibm3out_fib
     &     ,fibm1out_fsh,fibm2out_fsh,fibm3out_fsh
     &     ,fibm1out_esh,fibm2out_esh,fibm3out_esh
     &     ,fibm1out_fbc,fibm2out_fbc,fibm3out_fbc
     &     ,gx_ibm_massiveout_fib
     &     ,gy_ibm_massiveout_fib
     &     ,gz_ibm_massiveout_fib
     &     ,gx_ibm_massiveout_fsh
     &     ,gy_ibm_massiveout_fsh
     &     ,gz_ibm_massiveout_fsh
     &     ,gx_ibm_massiveout_esh
     &     ,gy_ibm_massiveout_esh
     &     ,gz_ibm_massiveout_esh
     &     ,gx_ibm_massiveout_fbc
     &     ,gy_ibm_massiveout_fbc
     &     ,gz_ibm_massiveout_fbc
     &     ,theboss
     &       )
      USE HeaderFSI
!@T This subroutine calls the structural solver and pass {\bf Vicar3D}
!variables from subroutine [[solve_membranetest]] or [[solve_membrane]] \\

! this subroutines are called from [[UTIL_BOUNDARY_FORCE_FSI.F90]] file.\\
! Depending on the mass ratio of the structure $m_s$ different initial values
! is assumed for the calculation at time $t (n)$ based on the variables at time
! $t-\Delta t (n-1)$. \\
!$ x^*=2 * x^{n} - x^{n-1}\/\/ \textrm{if}\/\/  m_s\/\/  \ge\/\/  1.05$ \\
!$ x^*= x^{n}\/\/\/\/\/\/  \textrm{if}\/\/   m_s\/\/  \ge\/\/  1.05$ \\
!@q
      implicit none
      real*8 UIBM1in_fib(Nr_IBM_fib,nsec_IBMmax,Ns_IBM_fib)
     &      ,UIBM2in_fib(Nr_IBM_fib,nsec_IBMmax,Ns_IBM_fib)
     &      ,UIBM3in_fib(Nr_IBM_fib,nsec_IBMmax,Ns_IBM_fib) 
     &      ,UIBM1in_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,UIBM2in_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,UIBM3in_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh) 
     &      ,UIBM1in_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,UIBM2in_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,UIBM3in_esh(Nr_IBM_esh,Ns_IBM_esh) 
     &      ,UIBM1in_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
     &      ,UIBM2in_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
     &      ,UIBM3in_fbc(Nr_IBM_fbc,Ns_IBM_fbc) 
      real*8 FIBM1out_fib(Nr_IBM_fib,nsec_IBMmax,Ns_IBM_fib)
     &      ,FIBM2out_fib(Nr_IBM_fib,nsec_IBMmax,Ns_IBM_fib)
     &      ,FIBM3out_fib(Nr_IBM_fib,nsec_IBMmax,Ns_IBM_fib)
     &      ,FIBM1out_fsh(Nr_IBM,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,FIBM2out_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,FIBM3out_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,FIBM1out_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,FIBM2out_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,FIBM3out_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,FIBM1out_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
     &      ,FIBM2out_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
     &      ,FIBM3out_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
      real*8 GX_IBM_MASSIVEout_fib(Nr_IBM_fib,Ns_IBM_fib)
     &      ,Gy_IBM_MASSIVEout_fib(Nr_IBM_fib,Ns_IBM_fib)
     &      ,Gz_IBM_MASSIVEout_fib(Nr_IBM_fib,Ns_IBM_fib)
     &      ,GX_IBM_MASSIVEout_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,Gy_IBM_MASSIVEout_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,Gz_IBM_MASSIVEout_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,GX_IBM_MASSIVEout_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,Gy_IBM_MASSIVEout_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,Gz_IBM_MASSIVEout_esh(Nr_IBM_esh,Ns_IBM_esh)  
     &      ,GX_IBM_MASSIVEout_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
     &      ,Gy_IBM_MASSIVEout_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
     &      ,Gz_IBM_MASSIVEout_fbc(Nr_IBM_fbc,Ns_IBM_fbc) 
      real*8 rtime, dt
      logical monitorON,theboss
    
      integer nt, ntsave,i,j,iiter2,jq,jq2
      integer ibm_solve,isavesigma
      integer isec,i_glb, ne1,i1
      real*8 x(3,3) ,xcont(3,3), xp1cont(3)
     &          ,xp2cont(3), outcont(3),outcont2(3)

      ibm_solve=0
! if converged
! check to see we need resorting of buckets "for contact algorithms"
      timebucket=timebucket+dt
      if(timebucket .gt. Bucketntinterval) then
            call  BucketSorting
            timebucket=0.0
      endif
#ifdef INCLUDE_FIB
! passing old fiber values to pre variables and initialize new variables 
      do i=1,nr_ibm_fib
       do j=1,ns_ibm_r_fib(i)
            ! gx_ibm_fib data at t^{n-1}
          gx_ibmpre_fib(i,j)=gx_ibm_fib(i,j)
          gy_ibmpre_fib(i,j)=gy_ibm_fib(i,j)
            ! gx_ibm_massive_fib data at t^{n}
          gx_ibm_fib(i,j)=gx_ibm_massive_fib(i,j)
          gy_ibm_fib(i,j)=gy_ibm_massive_fib(i,j)
            ! gx_ibm_massive_fib data at t^{n}
          gx_ibmo1_fib(i,j)=gx_ibm_massive_fib(i,j)
          gy_ibmo1_fib(i,j)=gy_ibm_massive_fib(i,j)
       end do
       if(ndim .gt.2) then

       do j=1,ns_ibm_r_fib(i)
          gz_ibmpre_fib(i,j)=gz_ibm_fib(i,j)
          gz_ibm_fib(i,j)=gz_ibm_massive_fib(i,j)
          gz_ibmo1_fib(i,j)=gz_ibm_massive_fib(i,j)
       end do
       endif
       do isec=1,nsec_ibm
       do j=1,ns_ibm_r_fib(i)
          UIBM1_fib(i,isec,j)=UIBM1in_fib(i,isec,j)
          UIBM2_fib(i,isec,j)=UIBM2in_fib(i,isec,j)
          fk_mass1o_fib(i,isec,j)=fk_mass1_fib(i,isec,j)
          fk_mass2o_fib(i,isec,j)=fk_mass2_fib(i,isec,j)
       end do
       if(ndim .gt.2) then
       do j=1,ns_ibm_r_fib(i)
          UIBM3_fib(i,isec,j)=UIBM3in_fib(i,isec,j)
          fk_mass3o_fib(i,isec,j)=fk_mass3_fib(i,isec,j)
       end do
       endif
      end do
      end do
#endif
#ifdef INCLUDE_FSH
! passing old fibrous shells values to pre variables and initialize new variables
      do i=1,nr_ibm_fsh
       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
          gx_ibmpre_fsh(i,jq,j)=gx_ibm_fsh(i,jq,j)
          gy_ibmpre_fsh(i,jq,j)=gy_ibm_fsh(i,jq,j)
          gx_ibm_fsh(i,jq,j)=gx_ibm_massive_fsh(i,jq,j)
          gy_ibm_fsh(i,jq,j)=gy_ibm_massive_fsh(i,jq,j)
          gx_ibmo1_fsh(i,jq,j)=gx_ibm_massive_fsh(i,jq,j)
          gy_ibmo1_fsh(i,jq,j)=gy_ibm_massive_fsh(i,jq,j)
       end do
       end do
       if(ndim .gt.2) then
       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
          gz_ibmpre_fsh(i,jq,j)=gz_ibm_fsh(i,jq,j)
          gz_ibm_fsh(i,jq,j)=gz_ibm_massive_fsh(i,jq,j)
          gz_ibmo1_fsh(i,jq,j)=gz_ibm_massive_fsh(i,jq,j)
       end do
       end do
       endif

       i_glb=Ilocglb_fsh(i);
       if(Electromechanical(i_glb)) then
         do jq=1,nq_ibm_r_fsh(i)
         do j=1,ns_ibm_r_fsh(i)
            ndotpre_ibm_fsh(i,jq,j)=ndotold_ibm_fsh(i,jq,j)
            ndotold_ibm_fsh(i,jq,j)=ndot_ibm_fsh(i,jq,j)
         end do
         end do

         if(Piezo_Dist_Flag(i_glb)) then
           do jq=1,nq_ibm_r_fsh(i)
           do j=1,ns_ibm_r_fsh(i)
              vPiezopre_fsh(i,jq,j)=vPiezoold_fsh(i,jq,j)
              vPiezoold_fsh(i,jq,j)=vPiezo_fsh(i,jq,j)
           enddo
           enddo
          else
           do jq=1,1
           do j=1,1
              vPiezopre_fsh(i,jq,j)=vPiezoold_fsh(i,jq,j)
              vPiezoold_fsh(i,jq,j)=vPiezo_fsh(i,jq,j)
           enddo
           enddo
          endif
       endif

       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
          UIBM1_fsh(i,jq,j)=UIBM1in_fsh(i,jq,j)
          UIBM2_fsh(i,jq,j)=UIBM2in_fsh(i,jq,j)
          fk_mass1o_fsh(i,jq,j)=fk_mass1_fsh(i,jq,j)
          fk_mass2o_fsh(i,jq,j)=fk_mass2_fsh(i,jq,j)
       end do
       end do
       if(ndim .gt.2) then
       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
          UIBM3_fsh(i,jq,j)=UIBM3in_fsh(i,jq,j)
          fk_mass3o_fsh(i,jq,j)=fk_mass3_fsh(i,jq,j)
       end do
       end do
       endif
      end do
#endif
#ifdef INCLUDE_ESH
! passing old FE shells values to pre variables and initialize new variables
      do i=1,nr_ibm_esh
       i_glb=Ilocglb_esh(i)
       do j=1,ns_ibm_r_esh(i)
          gx_ibmpre_esh(i,j)=gx_ibm_esh(i,j)
          gy_ibmpre_esh(i,j)=gy_ibm_esh(i,j)
          gx_ibm_esh(i,j)=gx_ibm_massive_esh(i,j)
          gy_ibm_esh(i,j)=gy_ibm_massive_esh(i,j)
          gx_ibmo1_esh(i,j)=gx_ibm_massive_esh(i,j)
          gy_ibmo1_esh(i,j)=gy_ibm_massive_esh(i,j)
       end do

       if(ndim .gt.2) then
       do j=1,ns_ibm_r_esh(i)
          gz_ibmpre_esh(i,j)=gz_ibm_esh(i,j)
          gz_ibm_esh(i,j)=gz_ibm_massive_esh(i,j)
          gz_ibmo1_esh(i,j)=gz_ibm_massive_esh(i,j)
       end do
       endif

       do j=1,ns_ibm_r_esh(i)
          UIBM1_esh(i,j)=UIBM1in_esh(i,j)
          UIBM2_esh(i,j)=UIBM2in_esh(i,j)
          fk_mass1o_esh(i,j)=fk_mass1_esh(i,j)
          fk_mass2o_esh(i,j)=fk_mass2_esh(i,j)
       end do
       if(ndim .gt.2) then
       do j=1,ns_ibm_r_esh(i)
          UIBM3_esh(i,j)=UIBM3in_esh(i,j)
          fk_mass3o_esh(i,j)=fk_mass3_esh(i,j)
       end do
       endif
      if(flag_fractal_esh(i)) then
         jq2=nifractal_esh(i)
       if(Genalpha_timesolver(i_glb) .eq. 0) then
       if(fractalinitFlag_esh(i)) then
         do j=1,3

         do jq=1,jq2
            Afractal2_esh(1:jq2,jq,j,i)=-Afractal_esh(1:jq2,jq,j,i)*dt
         enddo
         do jq=1,jq2
            Afractal2_esh(jq,jq,j,i)=Afractal2_esh(jq,jq,j,i)
     &                              +1.0
         enddo
! call inverse( Afractal2_esh(1:jq2,1:jq2,j,i), Afractal2_esh(1:jq2,1:jq2,j,i),jq2)
         print *,"inverse is commneted out to compile with gnu"
         stop
       end do
       fractalinitFlag_esh(i)=.false.
      endif
         jq=nifractal_esh(i)
         jq2=3*jq
         do ne1=1,numelr(i_glb)
           do i1=1,ngau 
              Xfractal_esh(1:jq,i1,ne1,i)=
     &        Xfractal_esh(jq2+1:jq2+jq,i1,ne1,i)

              Xfractal_esh(jq+1:2*jq,i1,ne1,i)=
     &        Xfractal_esh(jq2+jq+1:jq2+2*jq,i1,ne1,i)

              Xfractal_esh(2*jq+1:3*jq,i1,ne1,i)=
     &        Xfractal_esh(jq2+2*jq+1:jq2+3*jq,i1,ne1,i)
           enddo
         enddo  

      else
       if(fractalinitFlag_esh(i)) then
         jq2=nifractal_esh(i)
         do j=1,3

         do jq=1,jq2
            Afractal2_esh(1:jq2,jq,j,i)=
     &                              -(1.0-alphaf_Coef(i_glb))
     &                              *Afractal_esh(1:jq2,jq,j,i)
         enddo
         do jq=1,jq2
            Afractal2_esh(jq,jq,j,i)=Afractal2_esh(jq,jq,j,i)
     &               +(1.0-alpham_Coef(i_glb))/(dt*gamma_Coef(i_glb))
         enddo

! call inverse( Afractal2_esh(1:jq2,1:jq2,j,i),Afractal2_esh(1:jq2,1:jq2,j,i),jq2)
         print *,"inverse is commneted out to compile with gnu"
         stop

         do jq=1,jq2
            Bfractal2_esh(1:jq2,jq,j,i)
     &                =Afractal_esh(1:jq2,jq,j,i)
            Cfractal2_esh(1:jq2,jq,j,i)=0.0
         enddo

         do jq=1,jq2         
            Bfractal2_esh(jq,jq,j,i)=Bfractal2_esh(jq,jq,j,i)
     &              +(1.0-alpham_Coef(i_glb))/(dt*gamma_Coef(i_glb))
     &              -(1.0-alphaf_Coef(i_glb))
            Cfractal2_esh(jq,jq,j,i)=
     &                -1.0
     &                -(1.0-alpham_Coef(i_glb))/gamma_Coef(i_glb)
         enddo                

         Hfractal2_esh(1:jq2,j,i)=(1.0-alphaf_Coef(i_glb))
     &        *matmul(Cfractal_esh(1:jq2,j,i),
     &            matmul(Afractal2_esh(1:jq2,1:jq2,j,i)
     &                  ,Bfractal2_esh(1:jq2,1:jq2,j,i) ) )
     &        +alphaf_Coef(i_glb)*Cfractal_esh(1:jq2,j,i)

         Ifractal2_esh(1:jq2,j,i)=(1.0-alphaf_Coef(i_glb))
     &        *matmul(Cfractal_esh(1:jq2,j,i),
     &            matmul(Afractal2_esh(1:jq2,1:jq2,j,i)
     &                  ,Cfractal2_esh(1:jq2,1:jq2,j,i) ) )

         Jfractal2_esh(j,i)=(1.0-alphaf_Coef(i_glb))
     &        *dot_product(Cfractal_esh(1:jq2,j,i),
     &            matmul(Afractal2_esh(1:jq2,1:jq2,j,i)
     &                  ,Bfractal_esh(1:jq2,j,i) ) )
     &        +Dfractal_esh(j,i)


          Bfractal2_esh(1:jq2,1:jq2,j,i)
     &          =matmul(Afractal2_esh(1:jq2,1:jq2,j,i)
     &                 ,Bfractal2_esh(1:jq2,1:jq2,j,i) )

          Cfractal2_esh(1:jq2,1:jq2,j,i)
     &          =matmul(Afractal2_esh(1:jq2,1:jq2,j,i)
     &                 ,Cfractal2_esh(1:jq2,1:jq2,j,i) )

          Dfractal2_esh(1:jq2,j,i)
     &          =matmul(Afractal2_esh(1:jq2,1:jq2,j,i)
     &                 ,Bfractal_esh(1:jq2,j,i))
       end do
       fractalinitFlag_esh(i)=.false.
      endif
      endif  !Genalpha_timesolver(i_glb) 
      endif  !flag_fractal_esh(i)
      end do
#endif
#ifdef INCLUDE_FBC
! passing old fabric shells values to pre variables and initialize new variables
      do i=1,nr_ibm_fbc
       do j=1,ns_ibm_r_fbc(i)
          gx_ibmpre_fbc(i,j)=gx_ibm_fbc(i,j)
          gy_ibmpre_fbc(i,j)=gy_ibm_fbc(i,j)
          gx_ibm_fbc(i,j)=gx_ibm_massive_fbc(i,j)
          gy_ibm_fbc(i,j)=gy_ibm_massive_fbc(i,j)
          gx_ibmo1_fbc(i,j)=gx_ibm_massive_fbc(i,j)
          gy_ibmo1_fbc(i,j)=gy_ibm_massive_fbc(i,j)
       end do

       if(ndim .gt.2) then
       do j=1,ns_ibm_r_fbc(i)
          gz_ibmpre_fbc(i,j)=gz_ibm_fbc(i,j)
          gz_ibm_fbc(i,j)=gz_ibm_massive_fbc(i,j)
          gz_ibmo1_fbc(i,j)=gz_ibm_massive_fbc(i,j)
       end do
       endif

       if(FabricMeshFlag_fbc(i)) then
        do j=1,FabricMesh_numrecord(i)
          FabricMesh_coordpre(i,j,1:3)= FabricMesh_coord(i,j,1:3)         
          FabricMesh_coord(i,j,1:3)= FabricMesh_coordMass(i,j,1:3) 
          FabricMesh_coordo1(i,j,1:3)= FabricMesh_coordMass(i,j,1:3) 
        enddo
       endif
       if(FabricLineFlag_fbc(i)) then
        do j=1,FabricLine_numrecord(i)
          FabricLine_coordpre(i,j,1:3)= FabricLine_coord(i,j,1:3)         
          FabricLine_coord(i,j,1:3)= FabricLine_coordMass(i,j,1:3) 
          FabricLine_coordo1(i,j,1:3)= FabricLine_coordMass(i,j,1:3) 
        enddo
       endif

       do j=1,ns_ibm_r_fbc(i)
          UIBM1_fbc(i,j)=UIBM1in_fbc(i,j)
          UIBM2_fbc(i,j)=UIBM2in_fbc(i,j)
          fk_mass1o_fbc(i,j)=fk_mass1_fbc(i,j)
          fk_mass2o_fbc(i,j)=fk_mass2_fbc(i,j)
       end do
       if(ndim .gt.2) then
       do j=1,ns_ibm_r_fbc(i)
          UIBM3_fbc(i,j)=UIBM3in_fbc(i,j)
          fk_mass3o_fbc(i,j)=fk_mass3_fbc(i,j)
       end do
       endif
      end do
#endif
      
#ifdef INCLUDE_FIB
! computing predicted x^*=2 * x^{n} - x^{n-1}
      do i=1,nr_ibm_fib
       i_glb=Ilocglb_fib(i);
       do j=1,ns_ibm_r_fib(i)
          gx_ibm_massiveo_fib(i,j)=
     &                 2.0d0*gx_ibm_fib(i,j)-gx_ibmpre_fib(i,j)
          gy_ibm_massiveo_fib(i,j)=
     &                 2.0d0*gy_ibm_fib(i,j)-gy_ibmpre_fib(i,j)
          if(ndim .gt.2) gz_ibm_massiveo_fib(i,j)=
     &                 2.0d0*gz_ibm_fib(i,j)-gz_ibmpre_fib(i,j)

        if(density_coef(i_glb) .ge. 1.05) then
          gx_ibm_massive_fib(i,j)=gx_ibm_massiveo_fib(i,j)
          gy_ibm_massive_fib(i,j)=gy_ibm_massiveo_fib(i,j)
          if(ndim .gt.2) gz_ibm_massive_fib(i,j)=
     &                   gz_ibm_massiveo_fib(i,j)
         else
          gx_ibm_massive_fib(i,j)=gx_ibm_fib(i,j)
          gy_ibm_massive_fib(i,j)=gy_ibm_fib(i,j)
          if(ndim .gt.2) gz_ibm_massive_fib(i,j)=gz_ibm_fib(i,j)
         endif

          fs_1_ibmo_fib(i,j)=fs_1_ibm_fib(i,j)
          gx_ibmo1_fib(i,j)=gx_ibm_massiveo_fib(i,j)
          gy_ibmo1_fib(i,j)=gy_ibm_massiveo_fib(i,j)
          if(ndim .gt.2) 
     &       gz_ibmo1_fib(i,j)=gz_ibm_massiveo_fib(i,j)
  
        end do
       end do
#endif
#ifdef INCLUDE_FSH
      do i=1,nr_ibm_fsh
       i_glb=Ilocglb_fsh(i);
       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
          gx_ibm_massiveo_fsh(i,jq,j)=
     &                 2.0d0*gx_ibm_fsh(i,jq,j)-gx_ibmpre_fsh(i,jq,j)
          gy_ibm_massiveo_fsh(i,jq,j)=
     &                 2.0d0*gy_ibm_fsh(i,jq,j)-gy_ibmpre_fsh(i,jq,j)
          if(ndim .gt.2) gz_ibm_massiveo_fsh(i,jq,j)=
     &                 2.0d0*gz_ibm_fsh(i,jq,j)-gz_ibmpre_fsh(i,jq,j)
        if(density_coef(i_glb) .ge. 1.05) then
          gx_ibm_massive_fsh(i,jq,j)=gx_ibm_massiveo_fsh(i,jq,j)
          gy_ibm_massive_fsh(i,jq,j)=gy_ibm_massiveo_fsh(i,jq,j)
          if(ndim .gt.2) gz_ibm_massive_fsh(i,jq,j)=
     &                   gz_ibm_massiveo_fsh(i,jq,j)
         else
          gx_ibm_massive_fsh(i,jq,j)=gx_ibm_fsh(i,jq,j)
          gy_ibm_massive_fsh(i,jq,j)=gy_ibm_fsh(i,jq,j)
          if(ndim .gt.2) gz_ibm_massive_fsh(i,jq,j)=gz_ibm_fsh(i,jq,j)
         endif

          fs_1_ibmo_fsh(i,jq,j)=fs_1_ibm_fsh(i,jq,j)
          gx_ibmo1_fsh(i,jq,j)=gx_ibm_massiveo_fsh(i,jq,j)
          gy_ibmo1_fsh(i,jq,j)=gy_ibm_massiveo_fsh(i,jq,j)
          if(ndim .gt.2) 
     &       gz_ibmo1_fsh(i,jq,j)=gz_ibm_massiveo_fsh(i,jq,j)
        end do
        end do
       end do
#endif
#ifdef INCLUDE_ESH
      do i=1,nr_ibm_esh
       i_glb=Ilocglb_esh(i);
       do j=1,ns_ibm_r_esh(i)
          gx_ibm_massiveo_esh(i,j)=
     &                 2.0d0*gx_ibm_esh(i,j)-gx_ibmpre_esh(i,j)
          gy_ibm_massiveo_esh(i,j)=
     &                 2.0d0*gy_ibm_esh(i,j)-gy_ibmpre_esh(i,j)
          if(ndim .gt.2) gz_ibm_massiveo_esh(i,j)=
     &                 2.0d0*gz_ibm_esh(i,j)-gz_ibmpre_esh(i,j)
        if(density_coef(i_glb) .ge. 1.05) then
          gx_ibm_massive_esh(i,j)=gx_ibm_massiveo_esh(i,j)
          gy_ibm_massive_esh(i,j)=gy_ibm_massiveo_esh(i,j)
          if(ndim .gt.2) gz_ibm_massive_esh(i,j)=
     &                   gz_ibm_massiveo_esh(i,j)
         else
          gx_ibm_massive_esh(i,j)=gx_ibm_esh(i,j)
          gy_ibm_massive_esh(i,j)=gy_ibm_esh(i,j)
          if(ndim .gt.2) gz_ibm_massive_esh(i,j)=gz_ibm_esh(i,j)
         endif

          fs_1_ibmo_esh(i,j)=fs_1_ibm_esh(i,j)
          gx_ibmo1_esh(i,j)=gx_ibm_massiveo_esh(i,j)
          gy_ibmo1_esh(i,j)=gy_ibm_massiveo_esh(i,j)
          if(ndim .gt.2) 
     &       gz_ibmo1_esh(i,j)=gz_ibm_massiveo_esh(i,j)
  
        end do
       end do
       isavesigma=0 
#endif

#ifdef INCLUDE_FBC
      do i=1,nr_ibm_fbc
       i_glb=Ilocglb_fbc(i)
       do j=1,ns_ibm_r_fbc(i)

          gx_ibm_massiveo_fbc(i,j)=
     &                 2.0d0*gx_ibm_fbc(i,j)-gx_ibmpre_fbc(i,j)
          gy_ibm_massiveo_fbc(i,j)=
     &                 2.0d0*gy_ibm_fbc(i,j)-gy_ibmpre_fbc(i,j)
          if(ndim .gt.2) gz_ibm_massiveo_fbc(i,j)=
     &                 2.0d0*gz_ibm_fbc(i,j)-gz_ibmpre_fbc(i,j)
        if(density_coef(i_glb) .ge. 1.05) then
          gx_ibm_massive_fbc(i,j)=
     &                 2.0d0*gx_ibm_fbc(i,j)-gx_ibmpre_fbc(i,j)
          gy_ibm_massive_fbc(i,j)=
     &                 2.0d0*gy_ibm_fbc(i,j)-gy_ibmpre_fbc(i,j)
          if(ndim .gt.2) gz_ibm_massive_fbc(i,j)=
     &                 2.0d0*gz_ibm_fbc(i,j)-gz_ibmpre_fbc(i,j)
         else
          gx_ibm_massive_fbc(i,j)=gx_ibm_fbc(i,j)
          gy_ibm_massive_fbc(i,j)=gy_ibm_fbc(i,j)
          if(ndim .gt.2) gz_ibm_massive_fbc(i,j)=gz_ibm_fbc(i,j)
         endif

          fs_1_ibmo_fbc(i,j)=fs_1_ibm_fbc(i,j)
          gx_ibmo1_fbc(i,j)=gx_ibm_massiveo_fbc(i,j)
          gy_ibmo1_fbc(i,j)=gy_ibm_massiveo_fbc(i,j)
          if(ndim .gt.2) 
     &       gz_ibmo1_fbc(i,j)=gz_ibm_massiveo_fbc(i,j)
  
        end do

       if(FabricMeshFlag_fbc(i)) then
        do j=1,FabricMesh_numrecord(i)
          FabricMesh_coordMasso(i,j,1:3)= 
     &      2.0*FabricMesh_coord(i,j,1:3) 
     &      -FabricMesh_coordpre(i,j,1:3)         
        enddo
       endif
       if(FabricLineFlag_fbc(i)) then
        do j=1,FabricLine_numrecord(i)
          FabricLine_coordMasso(i,j,1:3)= 
     &      2.0*FabricLine_coord(i,j,1:3) 
     &      -FabricLine_coordpre(i,j,1:3)   
        enddo
       endif

       if(density_coef(i_glb) .ge. 1.05) then
       if(FabricMeshFlag_fbc(i)) then
        do j=1,FabricMesh_numrecord(i)
          FabricMesh_coordMass(i,j,1:3)=FabricMesh_coordMasso(i,j,1:3)       
        enddo
       endif
       if(FabricLineFlag_fbc(i)) then
        do j=1,FabricLine_numrecord(i)
          FabricLine_coordMass(i,j,1:3)=FabricLine_coordMasso(i,j,1:3)  
        enddo
       endif
       else  !FORif(density_coef(i_glb) .ge. 1.05) then
       if(FabricMeshFlag_fbc(i)) then
        do j=1,FabricMesh_numrecord(i)
          FabricMesh_coordMass(i,j,1:3)= 
     &      FabricMesh_coord(i,j,1:3) 
        enddo
       endif
       if(FabricLineFlag_fbc(i)) then
        do j=1,FabricLine_numrecord(i)
          FabricLine_coordMass(i,j,1:3)= 
     &      FabricLine_coord(i,j,1:3) 
  
        enddo
       endif
       endif !FORif(density_coef(i_glb) .ge. 1.05) then
       end do
#endif

       delta_t=dt
       time=rtime
       time2=rtime
!perfrom iteration over structural solver
         do iiter2=1,maxiter

             call r_stang(ibm_solve,monitorON,isavesigma)
#ifdef INCLUDE_FIB
             do i=1,nr_ibm_fib
              do j=1,ns_ibm_r_fib(i)
               gx_ibmo1_fib(i,j)=gx_ibm_massive_fib(i,j)
               gy_ibmo1_fib(i,j)=gy_ibm_massive_fib(i,j)
               if(ndim .gt.2) 
     &            gz_ibmo1_fib(i,j)=gz_ibm_massive_fib(i,j)
              end do
             end do
#endif
#ifdef INCLUDE_FSH
             do i=1,nr_ibm_fsh
              do jq=1,nq_ibm_r_fsh(i)
              do j=1,ns_ibm_r_fsh(i)
               gx_ibmo1_fsh(i,jq,j)=gx_ibm_massive_fsh(i,jq,j)
               gy_ibmo1_fsh(i,jq,j)=gy_ibm_massive_fsh(i,jq,j)
               if(ndim .gt.2) 
     &            gz_ibmo1_fsh(i,jq,j)=gz_ibm_massive_fsh(i,jq,j)
              end do
              end do
             end do

#endif
#ifdef INCLUDE_ESH
             do i=1,nr_ibm_esh
              do j=1,ns_ibm_r_esh(i)
               gx_ibmo1_esh(i,j)=gx_ibm_massive_esh(i,j)
               gy_ibmo1_esh(i,j)=gy_ibm_massive_esh(i,j)
               if(ndim .gt.2) 
     &            gz_ibmo1_esh(i,j)=gz_ibm_massive_esh(i,j)
              end do
             end do

#endif
#ifdef INCLUDE_FBC
             do i=1,nr_ibm_fbc
              do j=1,ns_ibm_r_fbc(i)
               gx_ibmo1_fbc(i,j)=gx_ibm_massive_fbc(i,j)
               gy_ibmo1_fbc(i,j)=gy_ibm_massive_fbc(i,j)
               if(ndim .gt.2) 
     &            gz_ibmo1_fbc(i,j)=gz_ibm_massive_fbc(i,j)
              end do

             if(FabricMeshFlag_fbc(i)) then
                do j=1,FabricMesh_numrecord(i)
                  FabricMesh_coordo1(i,j,1:3)= 
     &              FabricMesh_coordMass(i,j,1:3) 
                enddo
              endif
              if(FabricLineFlag_fbc(i)) then
                do j=1,FabricLine_numrecord(i)
                  FabricLine_coordo1(i,j,1:3)= 
     &              FabricLine_coordMass(i,j,1:3) 
                enddo
              endif
             end do
#endif
          end do
          
! call one more time at the end of the iteration
         time_step=nt
         if(mod(time_step,ntsave).eq. 0 ) isavesigma=1         
         call r_stang(ibm_solve,monitorON,isavesigma)
!!?!!         call r_stang(ibm_solve,monitorON)
!!?!!         time_step=nt
!computing output forces

         if(theboss) then
         if (mod(time_step,save_force_int).lt. del_force) then
           if(mod(time_step,save_force_int)==0) then
               call compute_force(.true.)
           else
               call compute_force(.false.)
           end if
         endif

! write tecplot files for structural components
           if(mod(time_step,ntsave).eq. 0 ) then
              call outputstr(nt,ntsave)
           endif

! write geometrical restart file for Fe shell  structures
           if(     (flag_restartfile)  
     &        .and.(time_step .ge.t0_restartfile) 
     &        .and.(mod(time_step,dt_restartfile).eq. 0 ) 
     &       ) then
            call restartInput_esh(nt,nt-t0_restartfile,dt_restartfile)
           endif
         endif

! initialize and assign output variables
         if(nr_ibm_fib>=1) FIBM1out_fib=0.0
         if(nr_ibm_fsh>=1) FIBM1out_fsh=0.0
         if(nr_ibm_esh>=1) FIBM1out_esh=0.0
         if(nr_ibm_fbc>=1) FIBM1out_fbc=0.0

         if(nr_ibm_fib>=1) FIBM2out_fib=0.0
         if(nr_ibm_fsh>=1) FIBM2out_fsh=0.0
         if(nr_ibm_esh>=1) FIBM2out_esh=0.0
         if(nr_ibm_fbc>=1) FIBM2out_fbc=0.0

         if(nr_ibm_fib>=1) FIBM3out_fib=0.0
         if(nr_ibm_fsh>=1) FIBM3out_fsh=0.0
         if(nr_ibm_esh>=1) FIBM3out_esh=0.0
         if(nr_ibm_fbc>=1) FIBM3out_fbc=0.0

#ifdef INCLUDE_FIB
      do i=1,nr_ibm_fib
        do isec=1,nsec_ibm
        do j=1,ns_ibm_r_fib(i)
          FIBM1out_fib(i,isec,j)=ffluidsum1_fib(i,isec,j)
          FIBM2out_fib(i,isec,j)=ffluidsum2_fib(i,isec,j)
          if(ndim .gt.2) 
     &       FIBM3out_fib(i,isec,j)=ffluidsum3_fib(i,isec,j)
        end do
        end do
        do j=1,ns_ibm_r_fib(i)
          gx_ibm_massiveout_fib(i,j)=gx_ibm_massive_fib(i,j)
          gy_ibm_massiveout_fib(i,j)=gy_ibm_massive_fib(i,j)
          if(ndim .gt.2) 
     &       gz_ibm_massiveout_fib(i,j)=gz_ibm_massive_fib(i,j)
        end do
       end do
#endif

#ifdef INCLUDE_FSH
       do i=1,nr_ibm_fsh
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
          FIBM1out_fsh(i,jq,j)=ffluidsum1_fsh(i,jq,j)
          FIBM2out_fsh(i,jq,j)=ffluidsum2_fsh(i,jq,j)
          if(ndim .gt.2) 
     &       FIBM3out_fsh(i,jq,j)=ffluidsum3_fsh(i,jq,j)
        end do
        end do
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
          gx_ibm_massiveout_fsh(i,jq,j)=gx_ibm_massive_fsh(i,jq,j)
          gy_ibm_massiveout_fsh(i,jq,j)=gy_ibm_massive_fsh(i,jq,j)
          if(ndim .gt.2) 
     &       gz_ibm_massiveout_fsh(i,jq,j)=gz_ibm_massive_fsh(i,jq,j)
        end do
        end do
      end do
#endif

#ifdef INCLUDE_ESH
      do i=1,nr_ibm_esh
        do j=1,ns_ibm_r_esh(i)
          FIBM1out_esh(i,j)=ffluidsum1_esh(i,j)
          FIBM2out_esh(i,j)=ffluidsum2_esh(i,j)
          if(ndim .gt.2) 
     &       FIBM3out_esh(i,j)=ffluidsum3_esh(i,j)
        end do
        do j=1,ns_ibm_r_esh(i)
          gx_ibm_massiveout_esh(i,j)=gx_ibm_massive_esh(i,j)
          gy_ibm_massiveout_esh(i,j)=gy_ibm_massive_esh(i,j)
          if(ndim .gt.2) 
     &       gz_ibm_massiveout_esh(i,j)=gz_ibm_massive_esh(i,j)
        end do
      end do
#endif
#ifdef INCLUDE_FBC
      do i=1,nr_ibm_fbc
        do j=1,ns_ibm_r_fbc(i)
          FIBM1out_fbc(i,j)=ffluidsum1_fbc(i,j)
          FIBM2out_fbc(i,j)=ffluidsum2_fbc(i,j)
          if(ndim .gt.2) 
     &       FIBM3out_fbc(i,j)=ffluidsum3_fbc(i,j)
        end do
        do j=1,ns_ibm_r_fbc(i)
          gx_ibm_massiveout_fbc(i,j)=gx_ibm_massive_fbc(i,j)
          gy_ibm_massiveout_fbc(i,j)=gy_ibm_massive_fbc(i,j)
          if(ndim .gt.2) 
     &       gz_ibm_massiveout_fbc(i,j)=gz_ibm_massive_fbc(i,j)
        end do
      end do
#endif

      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
