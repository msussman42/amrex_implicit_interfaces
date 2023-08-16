c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c rtime = t^{n+1}
c
      subroutine tick(rtime
     &     ,dt,nt,monitorON,ntsave
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
      implicit none
      real*8 UIBM1in_fib(Nr_IBM_fib,Ns_IBM_fib)
     &      ,UIBM2in_fib(Nr_IBM_fib,Ns_IBM_fib)
     &      ,UIBM3in_fib(Nr_IBM_fib,Ns_IBM_fib) 
     &      ,UIBM1in_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,UIBM2in_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh)
     &      ,UIBM3in_fsh(Nr_IBM_fsh,Nq_IBM_fsh,Ns_IBM_fsh) 
     &      ,UIBM1in_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,UIBM2in_esh(Nr_IBM_esh,Ns_IBM_esh)
     &      ,UIBM3in_esh(Nr_IBM_esh,Ns_IBM_esh) 
     &      ,UIBM1in_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
     &      ,UIBM2in_fbc(Nr_IBM_fbc,Ns_IBM_fbc)
     &      ,UIBM3in_fbc(Nr_IBM_fbc,Ns_IBM_fbc) 
      real*8 FIBM1out_fib(Nr_IBM_fib,Ns_IBM_fib)
     &      ,FIBM2out_fib(Nr_IBM_fib,Ns_IBM_fib)
     &      ,FIBM3out_fib(Nr_IBM_fib,Ns_IBM_fib)
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
      integer i_glb, ne1,i1
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
      i1=0
      ne1=0
      jq=0
      jq2=0
      outcont=0.0d0
      outcont2=0.0d0
      x=0.0d0
      xcont=0.0d0
      xp1cont=0.0d0
      xp2cont=0.0d0
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
       do j=1,ns_ibm_r_fib(i)
          UIBM1_fib(i,j)=UIBM1in_fib(i,j)
          UIBM2_fib(i,j)=UIBM2in_fib(i,j)
          fk_mass1o_fib(i,j)=fk_mass1_fib(i,j)
          fk_mass2o_fib(i,j)=fk_mass2_fib(i,j)
       end do
       if(ndim .gt.2) then
        do j=1,ns_ibm_r_fib(i)
          UIBM3_fib(i,j)=UIBM3in_fib(i,j)
          fk_mass3o_fib(i,j)=fk_mass3_fib(i,j)
        end do
       endif
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
        do j=1,ns_ibm_r_fib(i)
          FIBM1out_fib(i,j)=ffluidsum1_fib(i,j)
          FIBM2out_fib(i,j)=ffluidsum2_fib(i,j)
          if(ndim .gt.2) 
     &       FIBM3out_fib(i,j)=ffluidsum3_fib(i,j)
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


      subroutine tick_fsh(
     &      ilo,ihi,jlo,jhi,klo,khi
     &     ,nscalar
     &     ,ilo_dom,ihi_dom
     &     ,jlo_dom,jhi_dom
     &     ,klo_dom,khi_dom
     &     ,rtime
     &     ,dt
     &     ,nt
     &     ,monitorON
     &     ,ntsave
     &     ,prev_uibm123in_fsh
     &     ,uibm123in_fsh
     &     ,uibm123out_fsh
     &     ,fibm123in_fsh
     &     ,prev_gxyz_ibm_massivein_fsh
     &     ,gxyz_ibm_massivein_fsh
     &     ,gxyz_ibm_massiveout_fsh
     &     ,ds_fsh_in
     &     ,ds_fsh_out
     &     ,prev_scalar_fsh_in
     &     ,scalar_fsh_in
     &     ,scalar_fsh_out
     &     ,theboss
     &       )
      USE HeaderFSI
      implicit none
      integer ilo,ihi,jlo,jhi,klo,khi
      integer nscalar
      integer ilo_dom,ihi_dom,jlo_dom,jhi_dom,klo_dom,khi_dom

      real*8 UIBM123out_fsh(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 prev_uibm123in_fsh(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 uibm123in_fsh(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 FIBM123in_fsh(
     &  Nr_IBM_fsh,ilo_dom:ihi_dom,jlo_dom:jhi_dom
     & ,klo_dom:khi_dom,1:3)
      real*8 GXYZ_IBM_MASSIVEout_fsh(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 GXYZ_IBM_MASSIVEin_fsh(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 prev_GXYZ_IBM_MASSIVEin_fsh(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,1:3)

      real*8 ds_fsh_out(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi)
      real*8 ds_fsh_in(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi)

      real*8 prev_scalar_fsh_in(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,nscalar)
      real*8 scalar_fsh_in(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,nscalar)
      real*8 scalar_fsh_out(
     &  Nr_IBM_fsh,ilo:ihi,jlo:jhi,klo:khi,nscalar)

      real*8 rtime, dt
      logical monitorON,theboss
    
      integer nt, ntsave
      integer i,j,iiter2,jq,jq2
      integer ibm_solve,isavesigma
      integer i_glb, ne1,i1
      real*8 x(3,3) ,xcont(3,3), xp1cont(3)
     &          ,xp2cont(3), outcont(3),outcont2(3)
      real*8 tmpv1val
      real*8 tmpv2val
      integer debug_fsh
      integer debug_user
      integer ntsave_debug

#ifdef INCLUDE_FSH

      debug_fsh=0
      ntsave_debug=1

      if (nscalar.eq.2) then
! do nothing
      else
       print *,"expecting nscalar.eq.2"
       stop
      endif

      ibm_solve=0
! if converged
! check to see we need resorting of buckets "for contact algorithms"
      timebucket=timebucket+dt
      if(timebucket .gt. Bucketntinterval) then
            call  BucketSorting
            timebucket=0.0
      endif
      do i=1,nr_ibm_fsh

       if ((ilo.le.1-Nq_IBMB_fsh).and. 
     &     (ilo.le.0).and.
     &     (ihi.ge.nq_ibm_r_fsh(i)+1).and.
     &     (ihi.ge.nq_ibm_r_fsh(i)+Nq_IBMB_fsh).and. 
     &     (jlo.le.1-Ns_IBMB_fsh).and.
     &     (jlo.le.0).and.
     &     (jhi.ge.ns_ibm_r_fsh(i)+1).and.
     &     (jhi.ge.ns_ibm_r_fsh(i)+Ns_IBMB_fsh).and. 
     &     (klo.eq.1).and. 
     &     (khi.eq.1)
     &     ) then
!       do nothing
       else
        print *,"ilo,ihi,jlo,jhi,klo, or khi invalid"
        stop
       endif

       if ((ilo_dom.eq.1).and. 
     &     (ihi_dom.ge.nq_ibm_r_fsh(i)).and.
     &     (jlo_dom.eq.1).and.
     &     (jhi_dom.ge.ns_ibm_r_fsh(i)).and.
     &     (klo_dom.eq.1).and.(khi_dom.eq.1)) then
!       do nothing
       else
        print *,"ilo_dom,ihi_dom,..., or khi_dom invalid"
        stop
       endif

       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
        gx_ibmpre_fsh(i,jq,j)=
     &   prev_gxyz_ibm_massivein_fsh(i,jq,j,klo,1)

        gy_ibmpre_fsh(i,jq,j)=
     &   prev_gxyz_ibm_massivein_fsh(i,jq,j,klo,2)

        gx_ibm_fsh(i,jq,j)= 
     &   gxyz_ibm_massivein_fsh(i,jq,j,klo,1)

        gy_ibm_fsh(i,jq,j)= 
     &   gxyz_ibm_massivein_fsh(i,jq,j,klo,2)

        gx_ibmo1_fsh(i,jq,j)= 
     &   gxyz_ibm_massivein_fsh(i,jq,j,klo,1)
        gy_ibmo1_fsh(i,jq,j)= 
     &   gxyz_ibm_massivein_fsh(i,jq,j,klo,2)

        gx_ibm_massive_fsh(i,jq,j)= 
     &   gxyz_ibm_massivein_fsh(i,jq,j,klo,1)
        gy_ibm_massive_fsh(i,jq,j)= 
     &   gxyz_ibm_massivein_fsh(i,jq,j,klo,2)

       end do
       end do

       if(ndim .gt.2) then
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)

         gz_ibmpre_fsh(i,jq,j)=
     &    prev_gxyz_ibm_massivein_fsh(i,jq,j,klo,3)

         gz_ibm_fsh(i,jq,j)= 
     &    gxyz_ibm_massivein_fsh(i,jq,j,klo,3)

         gz_ibmo1_fsh(i,jq,j)= 
     &    gxyz_ibm_massivein_fsh(i,jq,j,klo,3)

         gz_ibm_massive_fsh(i,jq,j)= 
     &    gxyz_ibm_massivein_fsh(i,jq,j,klo,3)
        end do
        end do
       endif

       i_glb=Ilocglb_fsh(i);

       if(Electromechanical(i_glb)) then
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
          ndotpre_ibm_fsh(i,jq,j)=
     &     prev_scalar_fsh_in(i,jq,j,klo,1)
          ndotold_ibm_fsh(i,jq,j)=
     &     scalar_fsh_in(i,jq,j,klo,1)
        end do
        end do

        if(Piezo_Dist_Flag(i_glb)) then
         do jq=1,nq_ibm_r_fsh(i)
         do j=1,ns_ibm_r_fsh(i)
          vPiezopre_fsh(i,jq,j)=
     &     prev_scalar_fsh_in(i,jq,j,klo,2)
          vPiezoold_fsh(i,jq,j)=
     &     scalar_fsh_in(i,jq,j,klo,2)
         enddo
         enddo
        else
         do jq=1,1
         do j=1,1
          vPiezopre_fsh(i,jq,j)=
     &     prev_scalar_fsh_in(i,jq,j,klo,2)
          vPiezoold_fsh(i,jq,j)=
     &     scalar_fsh_in(i,jq,j,klo,2)
         enddo
         enddo
        endif
       endif

       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
        UIBM1_fsh(i,jq,j)=0.0d0
        UIBM2_fsh(i,jq,j)=0.0d0
        vibm1_fsh(i,jq,j)=uibm123in_fsh(i,jq,j,klo,1)
        vibm2_fsh(i,jq,j)=uibm123in_fsh(i,jq,j,klo,2)
        fk_mass1o_fsh(i,jq,j)=FIBM123in_fsh(i,jq,j,1,1)
        if (abs(fk_mass1o_fsh(i,jq,j)).le.1.0D+20) then
!do nothing
        else
         print *,"fk_mass1o_fsh overflow tick.F"
         stop
        endif
        fk_mass2o_fsh(i,jq,j)=FIBM123in_fsh(i,jq,j,1,2)
        if (abs(fk_mass2o_fsh(i,jq,j)).le.1.0D+20) then
!do nothing
        else
         print *,"fk_mass2o_fsh overflow tick.F"
         stop
        endif
       end do
       end do

       do jq=1-Nq_IBMB_fsh,nq_ibm_r_fsh(i)+Nq_IBMB_fsh
       do j=1-Ns_IBMB_fsh,ns_ibm_r_fsh(i)+Ns_IBMB_fsh
! do nothing DS_IBM_fsh(i,jq,j) and DS2_IBM_fsh(i,jq,j) never change.
       enddo
       enddo

       if(ndim .gt.2) then
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
         UIBM3_fsh(i,jq,j)=0.0d0
         vibm3_fsh(i,jq,j)=uibm123in_fsh(i,jq,j,klo,3)
         fk_mass3o_fsh(i,jq,j)=FIBM123in_fsh(i,jq,j,1,3)
         if (abs(fk_mass3o_fsh(i,jq,j)).le.1.0D+20) then
!do nothing
         else
          print *,"fk_mass3o_fsh overflow tick.F"
          stop
         endif
        end do
        end do
       endif
! i
      end do

      do i=1,nr_ibm_fsh
       i_glb=Ilocglb_fsh(i);
       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
        gx_ibm_massiveo_fsh(i,jq,j)=
     &   2.0d0*gx_ibm_fsh(i,jq,j)-gx_ibmpre_fsh(i,jq,j)
        gy_ibm_massiveo_fsh(i,jq,j)=
     &   2.0d0*gy_ibm_fsh(i,jq,j)-gy_ibmpre_fsh(i,jq,j)
        if(ndim .gt.2) then
         gz_ibm_massiveo_fsh(i,jq,j)=
     &     2.0d0*gz_ibm_fsh(i,jq,j)-gz_ibmpre_fsh(i,jq,j)
        endif
        if(density_coef(i_glb) .ge. 1.05) then
         gx_ibm_massive_fsh(i,jq,j)=gx_ibm_massiveo_fsh(i,jq,j)
         gy_ibm_massive_fsh(i,jq,j)=gy_ibm_massiveo_fsh(i,jq,j)
         if(ndim .gt.2) then
          gz_ibm_massive_fsh(i,jq,j)=
     &     gz_ibm_massiveo_fsh(i,jq,j)
         endif
        else
         gx_ibm_massive_fsh(i,jq,j)=gx_ibm_fsh(i,jq,j)
         gy_ibm_massive_fsh(i,jq,j)=gy_ibm_fsh(i,jq,j)
         if(ndim .gt.2) then
          gz_ibm_massive_fsh(i,jq,j)=
     &      gz_ibm_fsh(i,jq,j)
         endif
        endif

        fs_1_ibmo_fsh(i,jq,j)=fs_1_ibm_fsh(i,jq,j)
        gx_ibmo1_fsh(i,jq,j)=gx_ibm_massiveo_fsh(i,jq,j)
        gy_ibmo1_fsh(i,jq,j)=gy_ibm_massiveo_fsh(i,jq,j)
        if(ndim .gt.2) then
         gz_ibmo1_fsh(i,jq,j)=gz_ibm_massiveo_fsh(i,jq,j)
        endif
       end do
       end do
      end do

      delta_t=dt
      time=rtime
      time2=rtime

      if(theboss) then

       if (debug_fsh.eq.0) then
! do nothing
       else if (debug_fsh.eq.1) then
        call outputstr(nt,ntsave_debug)
        print *,"debug_fsh before r_stang"
        call SOLIDFLUID_F90_KEYBOARD()
       else
        print *,"debug_fsh invalid"
        stop
       endif

      endif

!perfrom iteration over structural solver
      do iiter2=1,maxiter

       call r_stang_fsh(ibm_solve,monitorON,isavesigma)

       do i=1,nr_ibm_fsh
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
         gx_ibmo1_fsh(i,jq,j)=gx_ibm_massive_fsh(i,jq,j)
         gy_ibmo1_fsh(i,jq,j)=gy_ibm_massive_fsh(i,jq,j)
         if(ndim .gt.2) then
          gz_ibmo1_fsh(i,jq,j)=gz_ibm_massive_fsh(i,jq,j)
         endif
        end do
        end do
       end do

      end do
          
! call one more time at the end of the iteration
      time_step=nt
      if(mod(time_step,ntsave).eq. 0 ) then
       isavesigma=1         
      endif

      call r_stang_fsh(ibm_solve,monitorON,isavesigma)

      if(theboss) then

       if (debug_fsh.eq.0) then
! do nothing
       else if (debug_fsh.eq.1) then
        call outputstr(nt,ntsave_debug)
        print *,"debug_fsh after r_stang"
        call SOLIDFLUID_F90_KEYBOARD()
       else
        print *,"debug_fsh invalid"
        stop
       endif

      endif

      if(theboss) then

! write tecplot files for structural components
       if(mod(time_step,ntsave).eq. 0 ) then
        call outputstr(nt,ntsave)
       endif

      endif

      do i=1,nr_ibm_fsh

       i_glb=Ilocglb_fsh(i);

       if(Electromechanical(i_glb)) then
        do jq=1,nq_ibm_r_fsh(i)
        do j=1,ns_ibm_r_fsh(i)
          scalar_fsh_out(i,jq,j,klo,1)=
     &     ndot_ibm_fsh(i,jq,j)
        end do
        end do

        if(Piezo_Dist_Flag(i_glb)) then
         do jq=1,nq_ibm_r_fsh(i)
         do j=1,ns_ibm_r_fsh(i)
          scalar_fsh_out(i,jq,j,klo,2)=
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
        uIBM123out_fsh(i,jq,j,klo,1)=vibm1_fsh(i,jq,j)
        uIBM123out_fsh(i,jq,j,klo,2)=vibm2_fsh(i,jq,j)
        if(ndim .gt.2) then
         uIBM123out_fsh(i,jq,j,klo,3)=vibm3_fsh(i,jq,j)
        endif
       end do
       end do

       do jq=0,nq_ibm_r_fsh(i)+1
       do j=0,ns_ibm_r_fsh(i)+1
! position at t^{n+1}
        gxyz_ibm_massiveout_fsh(i,jq,j,klo,1)=
     &     gx_ibm_massive_fsh(i,jq,j)
        gxyz_ibm_massiveout_fsh(i,jq,j,klo,2)=
     &     gy_ibm_massive_fsh(i,jq,j)
        if(ndim .gt.2) then
         gxyz_ibm_massiveout_fsh(i,jq,j,klo,3)= 
     &     gz_ibm_massive_fsh(i,jq,j)
        endif
       end do
       end do

       do jq=1,nq_ibm_r_fsh(i)
       do j=1,ns_ibm_r_fsh(i)
        tmpv1val=1.0d0
        tmpv2val=1.0d0
        if(j .eq.1 .or. j .eq. ns_ibm_r_fsh(i)) tmpv1val=0.5d0
        if(jq.eq.1 .or. jq.eq. nq_ibm_r_fsh(i)) tmpv2val=0.5d0
        ds_fsh_out(i,jq,j,klo)=ds_ibm_fsh(i,jq,j)*ds2_ibm_fsh(i,jq,j)
     &           *tmpv1val*tmpv2val

        if (abs(ds_fsh_in(i,jq,j,klo)-
     &         ds_fsh_out(i,jq,j,klo)).le.1.0D-10) then
!        do nothing
        else
         print *,"expecting ds_fsh not to change"
         stop
        endif
       enddo
       enddo

      end do

#else 
      print *,"include_fsh not defined"
      stop
#endif

      return
      end

c rtime = t^{n+1}
      subroutine tick_fib(
     &      ilo,ihi,jlo,jhi,klo,khi
     &     ,ilo_dom,ihi_dom
     &     ,jlo_dom,jhi_dom
     &     ,klo_dom,khi_dom
! t^{n+1}              
     &     ,rtime
     &     ,dt
! step              
     &     ,nt
     &     ,monitorON
! plot_int
     &     ,ntsave
     &     ,prev_uibm123in_fib
     &     ,uibm123in_fib
     &     ,uibm123out_fib
     &     ,fibm123in_fib
     &     ,prev_gxyz_ibm_massivein_fib
     &     ,gxyz_ibm_massivein_fib
     &     ,gxyz_ibm_massiveout_fib
     &     ,ds_fib_in
     &     ,ds_fib_out
     &     ,theboss
     &       )
      USE HeaderFSI
      implicit none
      integer ilo,ihi,jlo,jhi,klo,khi
      integer ilo_dom,ihi_dom,jlo_dom,jhi_dom,klo_dom,khi_dom
      real*8 UIBM123out_fib(
     &  Nr_IBM_fib,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 prev_uibm123in_fib(
     &  Nr_IBM_fib,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 uibm123in_fib(
     &  Nr_IBM_fib,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 FIBM123in_fib(
     &  Nr_IBM_fib,ilo_dom:ihi_dom,jlo_dom:jhi_dom
     & ,klo_dom:khi_dom,1:3)
      real*8 GXYZ_IBM_MASSIVEout_fib(
     &  Nr_IBM_fib,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 GXYZ_IBM_MASSIVEin_fib(
     &  Nr_IBM_fib,ilo:ihi,jlo:jhi,klo:khi,1:3)
      real*8 prev_GXYZ_IBM_MASSIVEin_fib(
     &  Nr_IBM_fib,ilo:ihi,jlo:jhi,klo:khi,1:3)

      real*8 ds_fib_out(
     &  Nr_IBM_fib,ilo:ihi,jlo:jhi,klo:khi)
      real*8 ds_fib_in(
     &  Nr_IBM_fib,ilo:ihi,jlo:jhi,klo:khi)

      real*8 rtime, dt
      logical monitorON,theboss
    
      integer nt, ntsave
      integer i,j,iiter2
      integer ibm_solve,isavesigma
      integer i_glb
      integer debug_fib
      integer debug_user
      integer ntsave_debug

#ifdef INCLUDE_FIB
 
      debug_fib=0
      ntsave_debug=1

      ibm_solve=0
! if converged
! check to see we need resorting of buckets "for contact algorithms"
      timebucket=timebucket+dt
      if(timebucket .gt. Bucketntinterval) then
       call  BucketSorting
       timebucket=0.0
      endif
! passing old fiber values to pre variables and initialize new variables 
      do i=1,nr_ibm_fib
       if ((ilo.le.1-Ns_IBMB_fib).and. 
     &     (ilo.le.0).and.
     &     (ihi.ge.ns_ibm_r_fib(i)+1).and.
     &     (ihi.ge.ns_ibm_r_fib(i)+Ns_IBMB_fib).and. 
     &     (jlo.eq.1).and.(jhi.eq.1).and.
     &     (klo.eq.1).and.(khi.eq.1)) then
!       do nothing
       else
        print *,"ilo,ihi,jlo,jhi,klo, or khi invalid"
        stop
       endif

       if ((ilo_dom.eq.1).and. 
     &     (ihi_dom.ge.ns_ibm_r_fib(i)).and.
     &     (jlo_dom.eq.1).and.(jhi_dom.eq.1).and.
     &     (klo_dom.eq.1).and.(khi_dom.eq.1)) then
!       do nothing
       else
        print *,"ilo_dom,ihi_dom,..., or khi_dom invalid"
        stop
       endif

       do j=0,ns_ibm_r_fib(i)+1
          ! gx_ibm_fib data (position) at t^{n-1}
        gx_ibmpre_fib(i,j)=prev_gxyz_ibm_massivein_fib(i,j,jlo,klo,1)
        gy_ibmpre_fib(i,j)=prev_gxyz_ibm_massivein_fib(i,j,jlo,klo,2)
          ! gx_ibm_massive_fib data (position) at t^{n}
        gx_ibm_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,1)
        gy_ibm_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,2)
          ! gx_ibm_massive_fib data (position) at t^{n}
        gx_ibmo1_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,1)
        gy_ibmo1_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,2)
          ! gx_ibm_massive_fib data (position) at t^{n}
        gx_ibm_massive_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,1)
        gy_ibm_massive_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,2)
! j
       end do

       if(ndim .gt.2) then

        do j=0,ns_ibm_r_fib(i)+1
         gz_ibmpre_fib(i,j)=prev_gxyz_ibm_massivein_fib(i,j,jlo,klo,3)
         gz_ibm_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,3)
         gz_ibmo1_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,3)
         gz_ibm_massive_fib(i,j)=gxyz_ibm_massivein_fib(i,j,jlo,klo,3)
        end do

       endif

       do j=1,ns_ibm_r_fib(i)
        UIBM1_fib(i,j)=0.0d0
        UIBM2_fib(i,j)=0.0d0
        vibm1_fib(i,j)=uibm123in_fib(i,j,jlo,klo,1)
        vibm2_fib(i,j)=uibm123in_fib(i,j,jlo,klo,2)
! force
        fk_mass1o_fib(i,j)=FIBM123in_fib(i,j,1,1,1)
        if (abs(fk_mass1o_fib(i,j)).le.1.0D+20) then
!do nothing
        else
         print *,"fk_mass1 overflow tick.F 1232"
         stop
        endif
        fk_mass2o_fib(i,j)=FIBM123in_fib(i,j,1,1,2)
        if (abs(fk_mass2o_fib(i,j)).le.1.0D+20) then
!do nothing
        else
         print *,"fk_mass2 overflow tick.F 1239"
         stop
        endif
       end do

       do j=1-Ns_IBMB_fib,ns_ibm_r_fib(i)+Ns_IBMB_fib
! do nothing ds_ibm_fib not expected to change.
       enddo

       if(ndim .gt.2) then
        do j=1,ns_ibm_r_fib(i)
         UIBM3_fib(i,j)=0.0d0
         vibm3_fib(i,j)=uibm123in_fib(i,j,jlo,klo,3)
         fk_mass3o_fib(i,j)=FIBM123in_fib(i,j,1,1,3)
         if (abs(fk_mass3o_fib(i,j)).le.1.0D+20) then
!do nothing
         else
          print *,"fk_mass3 overflow tick.F 1239"
          stop
         endif
        end do
       endif

! i
      end do

! computing predicted x^*=2 * x^{n} - x^{n-1}
      do i=1,nr_ibm_fib
       i_glb=Ilocglb_fib(i);
       do j=1,ns_ibm_r_fib(i)
        gx_ibm_massiveo_fib(i,j)=
     &          2.0d0*gx_ibm_fib(i,j)-gx_ibmpre_fib(i,j)
        gy_ibm_massiveo_fib(i,j)=
     &          2.0d0*gy_ibm_fib(i,j)-gy_ibmpre_fib(i,j)
        if(ndim .gt.2) then
         gz_ibm_massiveo_fib(i,j)=
     &       2.0d0*gz_ibm_fib(i,j)-gz_ibmpre_fib(i,j)
        endif

        if(density_coef(i_glb) .ge. 1.05) then
         gx_ibm_massive_fib(i,j)=gx_ibm_massiveo_fib(i,j)
         gy_ibm_massive_fib(i,j)=gy_ibm_massiveo_fib(i,j)
         if(ndim .gt.2) then
          gz_ibm_massive_fib(i,j)=gz_ibm_massiveo_fib(i,j)
         endif
        else
         gx_ibm_massive_fib(i,j)=gx_ibm_fib(i,j)
         gy_ibm_massive_fib(i,j)=gy_ibm_fib(i,j)
         if(ndim .gt.2) then
          gz_ibm_massive_fib(i,j)=gz_ibm_fib(i,j)
         endif
        endif

        fs_1_ibmo_fib(i,j)=0.0d0
        gx_ibmo1_fib(i,j)=gx_ibm_massiveo_fib(i,j)
        gy_ibmo1_fib(i,j)=gy_ibm_massiveo_fib(i,j)
        if(ndim .gt.2) then
         gz_ibmo1_fib(i,j)=gz_ibm_massiveo_fib(i,j)
        endif
  
       end do
      end do

      delta_t=dt
      time=rtime
      time2=rtime

      if(theboss) then

       if (debug_fib.eq.0) then
! do nothing
       else if (debug_fib.eq.1) then
        call outputstr(nt,ntsave_debug)
        print *,"debug_fib before r_stang"
        call SOLIDFLUID_F90_KEYBOARD()
       else
        print *,"debug_fib invalid"
        stop
       endif

      endif

!perfrom iteration over structural solver
      do iiter2=1,maxiter

       call r_stang_fib(ibm_solve,monitorON,isavesigma)
       do i=1,nr_ibm_fib
        do j=1,ns_ibm_r_fib(i)
         gx_ibmo1_fib(i,j)=gx_ibm_massive_fib(i,j)
         gy_ibmo1_fib(i,j)=gy_ibm_massive_fib(i,j)
         if(ndim .gt.2) then
          gz_ibmo1_fib(i,j)=gz_ibm_massive_fib(i,j)
         endif
        end do
       end do
! iiter2
      end do
          
! call one more time at the end of the iteration
      time_step=nt
      if(mod(time_step,ntsave).eq. 0 ) then
       isavesigma=1  
      endif

      call r_stang_fib(ibm_solve,monitorON,isavesigma)

      if(theboss) then

       if (debug_fib.eq.0) then
! do nothing
       else if (debug_fib.eq.1) then
        call outputstr(nt,ntsave_debug)
        print *,"debug_fib after r_stang"
        call SOLIDFLUID_F90_KEYBOARD()
       else
        print *,"debug_fib invalid"
        stop
       endif

      endif


      if(theboss) then

! write tecplot files for structural components
! g[xyz]_ibm_massive_fib(i,j)
! vibm[1-3]_fib(i,j)
! filename='ibm-2dFib-B ... .dat'
! filename='ibm-3dFib-B ... .dat'
       if(mod(time_step,ntsave).eq. 0 ) then
        call outputstr(nt,ntsave)
       endif

      endif

      do i=1,nr_ibm_fib

       do j=1,ns_ibm_r_fib(i)
        uIBM123out_fib(i,j,jlo,klo,1)=vibm1_fib(i,j)
        uIBM123out_fib(i,j,jlo,klo,2)=vibm2_fib(i,j)
        if(ndim .gt.2) then
         uIBM123out_fib(i,j,jlo,klo,3)=vibm3_fib(i,j)
        endif
       end do

       do j=0,ns_ibm_r_fib(i)+1
! position at t^{n+1}
        gxyz_ibm_massiveout_fib(i,j,jlo,klo,1)=gx_ibm_massive_fib(i,j)
        gxyz_ibm_massiveout_fib(i,j,jlo,klo,2)=gy_ibm_massive_fib(i,j)
        if(ndim .gt.2) then
         gxyz_ibm_massiveout_fib(i,j,jlo,klo,3)=gz_ibm_massive_fib(i,j)
        endif
       end do

       do j=1,ns_ibm_r_fib(i)
        DS_fib_out(i,j,jlo,klo)=ds_IBM_fib(i,j)

        if (abs(ds_fib_in(i,j,jlo,klo)-
     &         ds_fib_out(i,j,jlo,klo)).le.1.0D-10) then
!        do nothing
        else
         print *,"expecting ds_ibm_fib not to change"
         stop
        endif
       enddo

      end do

#else 
      print *,"include_fib not defined"
      stop
#endif

      return
      end
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
