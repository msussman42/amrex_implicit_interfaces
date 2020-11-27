#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 0

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "PLIC_F.H"


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

#if (STANDALONE==1)
      module plic_cpp_module
      contains
#endif

! mask=0 at coarse/fine border cells and physical boundaries.
! mask=1 at fine/fine and periodic border cells.
! mask=1 at symmetric border cells
! update_flag=0 do not update error
! update_flag=1 update error
! vof,ref centroid,order,slope,intercept  x nmat
! last comp. of solid <0 in solid
! vof is inputs, slopes is output.

      ! masknbr:
      ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
      ! (2) =1 interior  =0 otherwise
      ! (3) =1 interior+ngrow-1  =0 otherwise
      ! (4) =1 interior+ngrow    =0 otherwise
      subroutine FORT_SLOPERECON( &
        tid, &
        gridno, &
        level, &
        finest_level, &
        max_level, &
        ngrow, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        xlo,dx, &
        masknbr,DIMS(masknbr), &
        snew,DIMS(snew), &
        vof,DIMS(vof), &
        LS,DIMS(LS), &
        slopes,DIMS(slopes), &
        nsteps, &
        time, &
        nmat,nten, &
        latent_heat, &
        update_flag, &
        total_calls, &
        total_iterations, &
        continuous_mof, &
        radius_cutoff)
#if (STANDALONE==0)
      use probf90_module
#elif (STANDALONE==1)
      use probcommon_module
#endif
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: gridno
      INTEGER_T, intent(in) :: level,finest_level,max_level
      INTEGER_T, intent(in) :: nsteps

      INTEGER_T, intent(in) :: nmat

      INTEGER_T, intent(in) :: radius_cutoff(nmat)

      INTEGER_T, intent(in) :: continuous_mof
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: update_flag
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: DIMDEC(masknbr)
      INTEGER_T, intent(in) :: DIMDEC(snew)
      INTEGER_T, intent(in) :: DIMDEC(vof)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(slopes)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
     
      REAL_T, intent(in) :: masknbr(DIMV(masknbr),4) 
      REAL_T, intent(in) :: vof(DIMV(vof),nmat*ngeom_raw) 
      REAL_T, intent(in) :: LS(DIMV(LS),nmat) 
      REAL_T, intent(out) :: slopes(DIMV(slopes),nmat*ngeom_recon) 
      REAL_T, intent(inout) :: snew(DIMV(snew),nmat*ngeom_raw+1) 
      REAL_T, intent(in) :: latent_heat(2*nten)
      
      INTEGER_T i,j,k,dir
      INTEGER_T igridlo(3),igridhi(3)

      INTEGER_T cmoflo(SDIM),cmofhi(SDIM)

      INTEGER_T im
      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T mofdata_super(nmat*ngeom_recon)
      REAL_T vof_super(nmat)
      REAL_T mofsten(nmat*ngeom_recon)
      REAL_T multi_centroidA(nmat,SDIM)

      INTEGER_T vofcomprecon
      INTEGER_T vofcompraw

#if (STANDALONE==0)
      REAL_T err,errsave
      INTEGER_T local_mask
      INTEGER_T stencil_valid
#elif (STANDALONE==1)
      ! do nothing
#else
      print *,"stand alone bust"
      stop
#endif
      REAL_T orderflag
      INTEGER_T total_calls(nmat)
      INTEGER_T total_iterations(nmat)
      INTEGER_T i1,j1,k1
      REAL_T LS_stencil(D_DECL(-1:1,-1:1,-1:1),nmat)

      INTEGER_T nmax

      INTEGER_T continuous_mof_parm
     
      INTEGER_T klosten,khisten
      INTEGER_T nhalf
      REAL_T xsten(-3:3,SDIM)
      REAL_T xstenbox(-1:1,SDIM)
      INTEGER_T nmat_in_cell
      INTEGER_T nmat_in_stencil
      INTEGER_T nten_test
      REAL_T volume_super
      REAL_T volume_super_mofdata
      REAL_T volsten,volmat
      REAL_T censten(SDIM)
      REAL_T cen_super(SDIM)
      REAL_T voflist(nmat)
      REAL_T voflist_stencil(nmat)
      REAL_T voflist_test
      INTEGER_T mof_verbose,use_ls_data,nhalfbox_sten
      REAL_T dxmaxLS
      INTEGER_T debugslope
      INTEGER_T tessellate
      INTEGER_T nhalf_box

#include "mofdata.H"

      nhalf_box=1

      tessellate=0

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      debugslope=0
      nhalf=3
      nhalfbox_sten=1

      if (bfact.lt.1) then
       print *,"bfact invalid170"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten inv. sloperecon nten, nten_test ",nten,nten_test
       stop
      endif

      if (nsteps.lt.0) then
       print *,"nsteps invalid in sloperecon, nsteps=",nsteps
       stop
      endif
      if ((gridno.lt.0).or. &
          (level.lt.0).or. &
          (level.gt.finest_level)) then
       print *,"grid or level bust"
       stop
      endif
      if (max_level.lt.finest_level) then
       print *,"max_level invalid sloperecon"
       stop
      endif 
     
      if (debugslope.eq.1) then
       print *,"BEFORE SLOPERECON --------------------------------"
       print *,"grid,level,finest ",gridno,level,finest_level
       print *,"STEP,TIME ",nsteps,time
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: SLOPE_RECON

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngrow.lt.1) then
       print *,"ngrow invalid in slope recon"
       stop
      endif
      if ((continuous_mof.ne.0).and. &
          (continuous_mof.ne.2).and. &
          (continuous_mof.ne.3).and. &
          (continuous_mof.ne.4).and. &
          (continuous_mof.ne.5)) then
       print *,"continuous_mof invalid"
       stop
      endif

      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid"
       stop
      endif

      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.3) then
       ! do nothing
      else
       print *,"levelrz invalid slope recon"
       stop
      endif

      if ((update_flag.eq.0).or.(update_flag.eq.1)) then
       ! do nothing
      else
       print *,"update_flag invalid1: ",update_flag
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(masknbr),1,-1,12)
      call checkbound(fablo,fabhi,DIMS(snew),1,-1,12)
      call checkbound(fablo,fabhi,DIMS(vof),1,-1,12)
      call checkbound(fablo,fabhi,DIMS(LS),1,-1,12)
      call checkbound(fablo,fabhi,DIMS(slopes),ngrow,-1,12)

      if (SDIM.eq.3) then
       klosten=-1
       khisten=1
      else if (SDIM.eq.2) then
       klosten=0
       khisten=0
      else
       print *,"dimension bust"
       stop
      endif

      call get_dxmaxLS(dx,bfact,dxmaxLS)

      call growntilebox(tilelo,tilehi, &
        fablo,fabhi,igridlo,igridhi,0)

      do i = igridlo(1),igridhi(1)
      do j = igridlo(2),igridhi(2)
      do k = igridlo(3),igridhi(3)

       if ((update_flag.eq.0).or.(update_flag.eq.1)) then
        ! do nothing
       else
        print *,"update_flag invalid loop ",update_flag
        print *,"i,j,k ",i,j,k
        print *,"igridlo ",igridlo(1),igridlo(2),igridlo(3)
        print *,"igridhi ",igridhi(1),igridhi(2),igridhi(3)
        stop
       endif
       call gridsten_level(xsten,i,j,k,level,nhalf)

       use_ls_data=1

       do dir=1,SDIM
        cmoflo(dir)=-1
        cmofhi(dir)=1
       enddo

       do im=1,nmat

        vofcomprecon=(im-1)*ngeom_recon+1
        vofcompraw=(im-1)*ngeom_raw+1
        do dir=0,SDIM
         mofdata(vofcomprecon+dir)=vof(D_DECL(i,j,k),vofcompraw+dir)
        enddo
        if ((mofdata(vofcomprecon).lt.-0.1).or. &
            (mofdata(vofcomprecon).gt.1.1)) then
         print *,"mofdata(vofcomprecon) invalid"
         print *,"mofdata(vofcomprecon)=",mofdata(vofcomprecon)
         stop
        endif

        orderflag=zero
        mofdata(vofcomprecon+SDIM+1)=orderflag

        do dir=SDIM+3,ngeom_recon
         mofdata(vofcomprecon+dir-1)=zero
        enddo

       enddo  ! im=1..nmat

        ! sum of F_fluid=1
        ! sum of F_rigid<=1
       nhalf_box=1
       call make_vfrac_sum_ok_base( &
         xsten,nhalf,nhalf_box, &
         bfact,dx, &
         tessellate,mofdata,nmat,SDIM,6)

       do im=1,nmat
        vofcomprecon=(im-1)*ngeom_recon+1
        voflist(im)=mofdata(vofcomprecon)
        voflist_stencil(im)=zero
       enddo ! im=1..nmat

       if ((level.ge.0).and. &
           (level.le.finest_level)) then

        do i1=-1,1
        do j1=-1,1
        do k1=klosten,khisten
         do im=1,nmat
          LS_stencil(D_DECL(i1,j1,k1),im)= &
            LS(D_DECL(i+i1,j+j1,k+k1),im)
          vofcompraw=(im-1)*ngeom_raw+1
          voflist_test=vof(D_DECL(i+i1,j+j1,k+k1),vofcompraw)
          if ((voflist_test.lt.-0.1).or. &
              (voflist_test.gt.1.1)) then
           print *,"voflist_test invalid"
           print *,"im,voflist_test= ",im,voflist_test
           print *,"i1,j1,k1 ",i1,j1,k1
           print *,"i,j,k ",i,j,k
           print *,"igridlo ",igridlo(1),igridlo(2),igridlo(3)
           print *,"igridhi ",igridhi(1),igridhi(2),igridhi(3)
           stop
          else if ((voflist_test.ge.-0.1).and. &
                   (voflist_test.le.1.1)) then
           if (voflist_test.gt.voflist_stencil(im)) then
            voflist_stencil(im)=voflist_test
           endif
          else
           print *,"voflist_test bust"
           stop
          endif
         enddo ! im
        enddo
        enddo
        enddo  ! i1,j1,k1 

        nmat_in_cell=0
        nmat_in_stencil=0
        do im=1,nmat
         if (is_rigid(nmat,im).eq.0) then
          if (voflist_stencil(im).gt.VOFTOL) then
           nmat_in_stencil=nmat_in_stencil+1
          endif
          if (voflist(im).gt.VOFTOL) then
           nmat_in_cell=nmat_in_cell+1
          endif
         else if (is_rigid(nmat,im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif
        enddo ! im=1..nmat
 
        if (nmat_in_cell.gt.nmat_in_stencil) then
         print *,"nmat_in_cell invalid"
         stop
        endif

        do im=1,nmat*ngeom_recon
         mofdata_super(im)=mofdata(im)
        enddo

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volume_super,cen_super,SDIM)

        if (nmat_in_cell.eq.1) then
         continuous_mof_parm=0
        else if ((nmat_in_cell.ge.2).and.(nmat_in_cell.le.nmat)) then

         continuous_mof_parm=continuous_mof

         if ((nmat_in_stencil.ge.3).and. &
             (nmat_in_stencil.le.nmat)) then
          if (continuous_mof.eq.3) then ! CLSVOF 2 materials, MOF > 2mat
           continuous_mof_parm=0 ! MOF
          else if ((continuous_mof.eq.4).or. & ! CLSVOF 2 mat., CMOF > 2 mat
                   (continuous_mof.eq.2)) then ! CMOF
           if (nmat_in_cell.eq.nmat_in_stencil) then
            continuous_mof_parm=2 ! CMOF
           else if (nmat_in_cell.lt.nmat_in_stencil) then
            continuous_mof_parm=2 ! CMOF
           else
            print *,"nmat_in_cell invalid"
            stop
           endif
          else if ((continuous_mof.eq.0).or. & ! mof
                   (continuous_mof.eq.5)) then ! clsvof
           ! do nothing
          else
           print *,"continuous_mof invalid"
           stop
          endif
         else if (nmat_in_stencil.eq.2) then
          if (nmat_in_cell.ne.2) then
           print *,"nmat_in_cell invalid"
           stop
          endif
          if ((continuous_mof.eq.3).or. & ! CLSVOF if 2 materials in stencil
              (continuous_mof.eq.4).or. & ! CLSVOF if 2 materials in stencil
              (continuous_mof.eq.5)) then ! CLSVOF if 2 materials in stencil
           continuous_mof_parm=5
          else if ((continuous_mof.eq.0).or. & ! MOF
                   (continuous_mof.eq.2)) then ! CMOF
           if (continuous_mof_parm.ne.continuous_mof) then
            print *,"continuous_mof_parm.ne.continuous_mof"
            stop
           endif
          else
           print *,"continuous_mof invalid"
           stop
          endif
         else
          print *,"nmat_in_cell or nmat_in_stencil invalid"
          stop
         endif
        
        else
         print *,"nmat_in_cell invalid"
         stop
        endif
           
         ! supercell for centroid cost function.
         ! center cell for volume constraint.
        if (continuous_mof_parm.eq.2) then
  
         volume_super=zero
         volume_super_mofdata=zero
         do dir=1,SDIM
          cen_super(dir)=zero
         enddo
         do im=1,nmat
          vofcomprecon=(im-1)*ngeom_recon+1
          do dir=1,SDIM
           mofdata_super(vofcomprecon+dir)=zero
          enddo
          vof_super(im)=zero
         enddo
         do i1=-1,1
         do j1=-1,1
         do k1=klosten,khisten
          call CISBOX(xstenbox,nhalfbox_sten, &
           xlo,dx,i+i1,j+j1,k+k1, &
           bfact,level, &
           volsten,censten,SDIM)
          volume_super=volume_super+volsten
          do dir=1,SDIM
           cen_super(dir)=cen_super(dir)+volsten*censten(dir)
          enddo

          do im=1,nmat
           vofcomprecon=(im-1)*ngeom_recon+1
           vofcompraw=(im-1)*ngeom_raw+1
           do dir=0,SDIM
            mofsten(vofcomprecon+dir)= &
             vof(D_DECL(i+i1,j+j1,k+k1),vofcompraw+dir)
           enddo
           orderflag=zero
           mofsten(vofcomprecon+SDIM+1)=orderflag
           do dir=SDIM+3,ngeom_recon
            mofsten(vofcomprecon+dir-1)=zero
           enddo
          enddo  ! im=1..nmat

           ! sum of F_fluid=1
           ! sum of F_rigid<=1
          nhalf_box=1
          call make_vfrac_sum_ok_base( &
            xstenbox,nhalfbox_sten,nhalf_box, &
            bfact,dx,tessellate,mofsten,nmat,SDIM,6)

          do im=1,nmat
           vofcomprecon=(im-1)*ngeom_recon+1
           volmat=volsten*mofsten(vofcomprecon)
           vof_super(im)=vof_super(im)+volmat
           do dir=1,SDIM
            mofdata_super(vofcomprecon+dir)= &
             mofdata_super(vofcomprecon+dir)+ &
             volmat*(censten(dir)+mofsten(vofcomprecon+dir))
           enddo ! dir
           if (is_rigid(nmat,im).eq.0) then
            volume_super_mofdata=volume_super_mofdata+volmat
           else if (is_rigid(nmat,im).eq.1) then
            ! do nothing
           else
            print *,"is_rigid(nmat,im) invalid"
            stop
           endif
          enddo ! im=1..nmat
         enddo
         enddo
         enddo ! i1,j1,k1

         if (volume_super.le.zero) then
          print *,"volume_super invalid"
          stop
         endif
         if (volume_super_mofdata.le.zero) then
          print *,"volume_super_mofdata invalid"
          stop
         endif
         do dir=1,SDIM
          cen_super(dir)=cen_super(dir)/volume_super
         enddo

         do im=1,nmat
          vofcomprecon=(im-1)*ngeom_recon+1
          if (vof_super(im).gt.zero) then
           do dir=1,SDIM
            mofdata_super(vofcomprecon+dir)= &
             mofdata_super(vofcomprecon+dir)/ &
             vof_super(im)- &
             cen_super(dir)
           enddo
           vof_super(im)=vof_super(im)/volume_super_mofdata
          else if (vof_super(im).eq.zero) then
           do dir=1,SDIM
            mofdata_super(vofcomprecon+dir)=zero
           enddo
          else
           print *,"vof_super(im) invalid"
           stop
          endif

           ! always standard MOF for the rigid materials.
          if (is_rigid(nmat,im).eq.1) then
           do dir=1,SDIM
            mofdata_super(vofcomprecon+dir)=mofdata(vofcomprecon+dir)
           enddo
          else if (is_rigid(nmat,im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif

         enddo ! im=1..nmat

        else if ((continuous_mof_parm.eq.0).or. &
                 (continuous_mof_parm.eq.5)) then
         ! do nothing
        else
         print *,"continuous_mof_parm invalid"
         stop
        endif

        mof_verbose=0

        call multimaterial_MOF( &
          bfact,dx,xsten,nhalf, &
          mof_verbose, &
          use_ls_data, &
          LS_stencil, &
          geom_xtetlist(1,1,1,tid+1), &
          geom_xtetlist_old(1,1,1,tid+1), &
          nmax, &
          nmax, &
          mofdata_super, &
          multi_centroidA, &
          continuous_mof_parm, &
          cmoflo,cmofhi, &
          nmat,SDIM,2)

        if (continuous_mof_parm.eq.2) then
           ! center cell centroids.
         do im=1,nmat
          vofcomprecon=(im-1)*ngeom_recon+1
          do dir=1,SDIM
           mofdata_super(vofcomprecon+dir)=mofdata(vofcomprecon+dir)
          enddo
         enddo ! im
        else if ((continuous_mof_parm.eq.0).or. &
                 (continuous_mof_parm.eq.5)) then
         ! do nothing
        else
         print *,"continuous_mof_parm invalid"
         stop
        endif

        do im=1,nmat
         total_calls(im)=total_calls(im)+mof_calls(tid+1,im)
         total_iterations(im)= &
          total_iterations(im)+mof_iterations(tid+1,im)
        enddo  ! im

       else
        print *,"level invalid slope recon"
        stop
       endif

       do dir=1,nmat*ngeom_recon
        slopes(D_DECL(i,j,k),dir)=mofdata_super(dir)
       enddo

       if (update_flag.eq.1) then

#if (STANDALONE==0)
        if (use_ls_data.ne.1) then
         print *,"use_ls_data invalid"
         stop
        endif

        if ((level.ge.0).and.(level.le.finest_level)) then

         stencil_valid=1
         do i1=-1,1
         do j1=-1,1
         do k1=klosten,khisten
          local_mask=NINT(masknbr(D_DECL(i+i1,j+j1,k+k1),1))
          if (local_mask.eq.1) then ! fine-fine ghost in domain or interior.
           ! do nothing
          else if (local_mask.eq.0) then ! ghost value is low order accurate.
           stencil_valid=0
          else
           print *,"local_mask invalid"
           stop
          endif
         enddo 
         enddo 
         enddo 

         call calc_error_indicator( &
          stencil_valid, &
          level,max_level, &
          xsten,nhalf,dx,bfact, &
          voflist, &
          LS_stencil, &
          nmat, &
          nten, &
          latent_heat, &
          radius_cutoff, &
          err,time)
         errsave=snew(D_DECL(i,j,k),nmat*ngeom_raw+1)
         if (errsave.lt.err) then
          snew(D_DECL(i,j,k),nmat*ngeom_raw+1)=err
         endif
         if ((errsave.lt.zero).or.(err.lt.zero)) then
          print *,"err bust"
          stop
         endif

        else
         print *,"level invalid slope recon 2"
         stop
        endif
#elif (STANDALONE==1)
        print *,"update_flag invalid for stand alone SLOPERECON"
        stop
#else
        print *,"bust compiling PLIC_3D.F90"
        stop
#endif
       else if (update_flag.eq.0) then
        ! do nothing
       else 
        print *,"update_flag invalid2: ",update_flag
        stop
       endif

      enddo
      enddo
      enddo

      if (debugslope.eq.1) then
       print *,"AFTER SLOPERECON --------------------------------"
       print *,"grid,level,finest ",gridno,level,finest_level
       print *,"STEP,TIME ",nsteps,time
      endif

      return
      end subroutine FORT_SLOPERECON



#if (STANDALONE==1)
      end module plic_cpp_module
#endif

#undef STANDALONE

