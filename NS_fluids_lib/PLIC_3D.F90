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

      module plic_cpp_module
      contains

! mask=0 at coarse/fine border cells and physical boundaries.
! mask=1 at fine/fine and periodic border cells.
! mask=1 at symmetric border cells
! update_flag=0 do not update error
! update_flag=1 update error
! vof,ref centroid,order,slope,intercept  x num_materials
! last comp. of solid <0 in solid
! vof is inputs, slopes is output.

      ! masknbr:
      ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
      ! (2) =1 interior  =0 otherwise
      ! (3) =1 interior+ngrow-1  =0 otherwise
      ! (4) =1 interior+ngrow    =0 otherwise
      subroutine fort_sloperecon( &
        tid, &
        gridno, &
        level, &
        finest_level, &
        max_level, &
        ngrow, &
        vofbc, &
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
        update_flag, &
        total_calls, &
        total_iterations, &
        continuous_mof, &
        force_cmof_at_triple_junctions, &
        partial_cmof_stencil_at_walls) &
      bind(c,name='fort_sloperecon')

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

      INTEGER_T, intent(in) :: continuous_mof
      INTEGER_T, intent(in) :: force_cmof_at_triple_junctions
      INTEGER_T, intent(in) :: partial_cmof_stencil_at_walls
      INTEGER_T, intent(in) :: update_flag
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: vofbc(SDIM,2)
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
     
      REAL_T, intent(in), target :: masknbr(DIMV(masknbr),4) 
      REAL_T, intent(in), target :: vof(DIMV(vof),num_materials*ngeom_raw) 
      REAL_T, intent(in), target :: LS(DIMV(LS),num_materials) 
      REAL_T, intent(out), target :: slopes(DIMV(slopes),num_materials*ngeom_recon) 
      REAL_T, pointer :: slopes_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(inout), target :: snew(DIMV(snew),num_materials*ngeom_raw+1) 
      REAL_T, pointer :: snew_ptr(D_DECL(:,:,:),:)
      
      INTEGER_T i,j,k,dir
      INTEGER_T igridlo(3),igridhi(3)

      INTEGER_T cmofsten(D_DECL(-1:1,-1:1,-1:1))

      INTEGER_T im
      REAL_T mofdata(num_materials*ngeom_recon)
      REAL_T mofdata_super(num_materials*ngeom_recon)
      REAL_T vof_super(num_materials)
      REAL_T mofsten(num_materials*ngeom_recon)
      REAL_T multi_centroidA(num_materials,SDIM)

      INTEGER_T vofcomprecon
      INTEGER_T vofcompraw

#if (STANDALONE==0)
      REAL_T err,errsave
      INTEGER_T local_mask
#elif (STANDALONE==1)
      ! do nothing
#else
      print *,"stand alone bust"
      stop
#endif
      REAL_T orderflag
      INTEGER_T total_calls(num_materials)
      INTEGER_T total_iterations(num_materials)
      INTEGER_T i1,j1,k1
      REAL_T LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)

      INTEGER_T nmax

      INTEGER_T continuous_mof_parm
      INTEGER_T continuous_mof_base
     
      INTEGER_T klosten,khisten
      INTEGER_T nhalf
      REAL_T xsten(-3:3,SDIM)
      REAL_T xstenbox(-1:1,SDIM)
      INTEGER_T num_materials_in_cell
      INTEGER_T num_materials_in_stencil
      INTEGER_T nten_test
      REAL_T volume_super
      REAL_T volume_super_mofdata
      REAL_T volsten
      REAL_T volmat
      REAL_T censten(SDIM)
      REAL_T cen_super(SDIM)
      REAL_T voflist_center(num_materials)
      REAL_T voflist_stencil(num_materials)
      REAL_T voflist_test
      INTEGER_T mof_verbose,use_ls_data,nhalfbox_sten
      REAL_T dxmaxLS
      INTEGER_T debugslope
      INTEGER_T tessellate
      INTEGER_T nhalf_box

      REAL_T vfrac_fluid_sum
      REAL_T vfrac_solid_sum
      REAL_T vfrac_solid_sum_center
      REAL_T vfrac_raster_solid
      REAL_T vfrac_local(num_materials)
      INTEGER_T im_raster_solid
      INTEGER_T mod_cmofsten
      INTEGER_T local_mod_cmofsten
      INTEGER_T loc_indx

#include "mofdata.H"

      nhalf_box=1

      slopes_ptr=>slopes
      snew_ptr=>snew

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
       print *,"BEFORE fort_sloperecon --------------------------------"
       print *,"grid,level,finest ",gridno,level,finest_level
       print *,"STEP,TIME ",nsteps,time
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: SLOPE_RECON

      if (ngrow.lt.1) then
       print *,"ngrow invalid in slope recon"
       stop
      endif
      if ((continuous_mof.eq.0).or. & ! MOF
          (continuous_mof.eq.2)) then ! CMOF
       ! do nothing
      else
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

      call checkbound_array(fablo,fabhi,masknbr,1,-1,12)
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1,12)
      call checkbound_array(fablo,fabhi,vof,1,-1,12)
      call checkbound_array(fablo,fabhi,LS,1,-1,12)
      call checkbound_array(fablo,fabhi,slopes_ptr,ngrow,-1,12)

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

       do i1=-1,1
       do j1=-1,1
       do k1=klosten,khisten
        cmofsten(D_DECL(i1,j1,k1))=1
       enddo
       enddo
       enddo

       do im=1,num_materials

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

       enddo  ! im=1..num_materials

        ! sum of F_fluid=1
        ! sum of F_rigid<=1
       nhalf_box=1
       call make_vfrac_sum_ok_base( &
         cmofsten, &
         xsten,nhalf,nhalf_box, &
         bfact,dx, &
         tessellate, & ! =0
         mofdata,SDIM,6)

       vfrac_fluid_sum=zero
       vfrac_solid_sum=zero
       im_raster_solid=0
       vfrac_raster_solid=zero
       mod_cmofsten=0
       local_mod_cmofsten=0

       do im=1,num_materials
        vofcomprecon=(im-1)*ngeom_recon+1
        voflist_center(im)=mofdata(vofcomprecon)
        voflist_stencil(im)=zero

        if (is_rigid(im).eq.0) then
         vfrac_fluid_sum=vfrac_fluid_sum+voflist_center(im)
        else if (is_rigid(im).eq.1) then
         if (im_raster_solid.eq.0) then
          im_raster_solid=im
          vfrac_raster_solid=voflist_center(im)
         else if ((im_raster_solid.ge.1).and. &
                  (im_raster_solid.le.num_materials).and. &
                  (is_rigid(im_raster_solid).eq.1)) then
          if (vfrac_raster_solid.lt.voflist_center(im)) then
           im_raster_solid=im
           vfrac_raster_solid=voflist_center(im)
          endif
         else
          print *,"im_raster_solid invalid"
          stop
         endif
      
         vfrac_solid_sum=vfrac_solid_sum+voflist_center(im)
        else
         print *,"is_rigid(im) invalid"
         stop
        endif

       enddo ! im=1..num_materials

       vfrac_solid_sum_center=vfrac_solid_sum

       if (abs(vfrac_fluid_sum-one).le.VOFTOL) then
        ! do nothing
       else
        print *,"vfrac_fluid_sum invalid"
        stop
       endif

       if ((level.ge.0).and. &
           (level.le.finest_level)) then

        do i1=-1,1
        do j1=-1,1
        do k1=klosten,khisten
         do im=1,num_materials
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

        num_materials_in_cell=0
        num_materials_in_stencil=0
        do im=1,num_materials
         if (is_rigid(im).eq.0) then
          if (voflist_stencil(im).gt.VOFTOL) then
           num_materials_in_stencil=num_materials_in_stencil+1
          endif
          if (voflist_center(im).gt.VOFTOL) then
           num_materials_in_cell=num_materials_in_cell+1
          endif
         else if (is_rigid(im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid PLIC_3D.F90"
          stop
         endif
        enddo ! im=1..num_materials
 
        if (num_materials_in_cell.gt.num_materials_in_stencil) then
         print *,"num_materials_in_cell invalid"
         stop
        endif

        do im=1,num_materials*ngeom_recon
         mofdata_super(im)=mofdata(im)
        enddo

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volume_super,cen_super,SDIM)

        mod_cmofsten=0

         ! check if fluid cell near a wall.
        if (vfrac_solid_sum_center.lt.half) then
         do i1=-1,1
         do j1=-1,1
         do k1=klosten,khisten

          do dir=1,SDIM
           if (dir.eq.1) then
            loc_indx=i+i1
           else if (dir.eq.2) then
            loc_indx=j+j1
           else if ((dir.eq.3).and.(SDIM.eq.3)) then
            loc_indx=k+k1
           else
            print *,"dir invalid"
            stop
           endif
           if (loc_indx.lt.fablo(dir)) then
            if (vofbc(dir,1).ne.INT_DIR) then
             mod_cmofsten=1
            endif
           else if (loc_indx.gt.fabhi(dir)) then
            if (vofbc(dir,2).ne.INT_DIR) then
             mod_cmofsten=1
            endif
           else if ((loc_indx.ge.fablo(dir)).and. &
                    (loc_indx.le.fabhi(dir))) then
            ! do nothing
           else
            print *,"loc_indx invalid"
            stop
           endif
          enddo ! dir=1..sdim

          if (mod_cmofsten.eq.0) then

           do im=1,num_materials
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
           enddo  ! im=1..num_materials

           call CISBOX(xstenbox,nhalfbox_sten, &
            xlo,dx,i+i1,j+j1,k+k1, &
            bfact,level, &
            volsten,censten,SDIM)

           ! sum of F_fluid=1
           ! sum of F_rigid<=1
           nhalf_box=1
           call make_vfrac_sum_ok_base( &
            cmofsten, &
            xstenbox,nhalfbox_sten,nhalf_box, &
            bfact,dx, &
            tessellate, & ! =0
            mofsten,SDIM,6)

           vfrac_fluid_sum=zero
           vfrac_solid_sum=zero

           do im=1,num_materials
            vofcomprecon=(im-1)*ngeom_recon+1
            vfrac_local(im)=mofsten(vofcomprecon)

            if (is_rigid(im).eq.0) then
             vfrac_fluid_sum=vfrac_fluid_sum+vfrac_local(im)
            else if (is_rigid(im).eq.1) then
             vfrac_solid_sum=vfrac_solid_sum+vfrac_local(im)
            else
             print *,"is_rigid(im) invalid"
             stop
            endif
           enddo ! im=1..num_materials

           if (abs(vfrac_fluid_sum-one).le.VOFTOL) then
            ! do nothing
           else
            print *,"vfrac_fluid_sum invalid"
            stop
           endif

           if (vfrac_solid_sum.ge.half) then
            if ((i1.eq.0).and.(j1.eq.0).and.(k1.eq.0)) then
             print *,"expecting i1 or j1 or k1 not 0"
             stop
            endif
            mod_cmofsten=1
           else if (vfrac_solid_sum.lt.half) then
            ! do nothing
           else
            print *,"vfrac_solid_sum invalid"
            stop
           endif

          else if (mod_cmofsten.eq.1) then
           ! do nothing
          else
           print *,"mod_cmofsten invalid"
           stop
          endif
         enddo !k1
         enddo !j1
         enddo !i1
        else if (vfrac_solid_sum_center.ge.half) then
         ! do nothing (ok to use MOF)
        else
         print *,"vfrac_solid_sum_center invalid"
         stop
        endif

        continuous_mof_base=continuous_mof

        if (mod_cmofsten.eq.1) then
         if (force_cmof_at_triple_junctions.eq.0) then
          ! do nothing
         else if ((force_cmof_at_triple_junctions.eq.1).or. &
                  (force_cmof_at_triple_junctions.eq.2)) then
          continuous_mof_base=2
         else
          print *,"force_cmof_at_triple_junctions invalid"
          stop
         endif
        else if (mod_cmofsten.eq.0) then
         ! do nothing
        else
         print *,"mod_cmofsten invalid"
         stop
        endif

        if (num_materials_in_cell.eq.1) then
         continuous_mof_parm=0
        else if ((num_materials_in_cell.ge.2).and. &
                 (num_materials_in_cell.le.num_materials)) then

         continuous_mof_parm=continuous_mof_base

         if ((num_materials_in_stencil.ge.3).and. &
             (num_materials_in_stencil.le.num_materials)) then

          if (force_cmof_at_triple_junctions.eq.2) then
           continuous_mof_base=2
          else if ((force_cmof_at_triple_junctions.eq.0).or. &
                   (force_cmof_at_triple_junctions.eq.1)) then
           ! do nothing
          else
           print *,"force_cmof_at_triple_junctions invalid"
           stop
          endif

          if (continuous_mof_base.eq.2) then ! CMOF
           if (num_materials_in_cell.eq.num_materials_in_stencil) then
            continuous_mof_parm=2 ! CMOF
           else if (num_materials_in_cell.lt.num_materials_in_stencil) then
            continuous_mof_parm=2 ! CMOF
           else
            print *,"num_materials_in_cell invalid"
            stop
           endif
          else if (continuous_mof_base.eq.0) then
           ! do nothing
          else
           print *,"continuous_mof_base invalid"
           stop
          endif
         else if (num_materials_in_stencil.eq.2) then
          if (num_materials_in_cell.ne.2) then
           print *,"num_materials_in_cell invalid"
           stop
          endif
          if ((continuous_mof_base.eq.0).or. & ! MOF
              (continuous_mof_base.eq.2)) then ! CMOF
           if (continuous_mof_parm.ne.continuous_mof_base) then
            print *,"continuous_mof_parm.ne.continuous_mof_base"
            stop
           endif
          else
           print *,"continuous_mof_base invalid"
           stop
          endif
         else
          print *,"num_materials_in_cell or num_materials_in_stencil invalid"
          stop
         endif
        
        else
         print *,"num_materials_in_cell invalid"
         stop
        endif

        mod_cmofsten=0
        local_mod_cmofsten=0

         ! supercell for centroid cost function.
         ! center cell for volume constraint.
        if (continuous_mof_parm.eq.2) then
  
         volume_super=zero ! volume of the extended region
         volume_super_mofdata=zero !same as volume_super, except by im_fluid.

         do dir=1,SDIM
          cen_super(dir)=zero
         enddo
         do im=1,num_materials
          vofcomprecon=(im-1)*ngeom_recon+1
          do dir=1,SDIM
           mofdata_super(vofcomprecon+dir)=zero
          enddo
          vof_super(im)=zero
         enddo ! im=1..num_materials

         do i1=-1,1
         do j1=-1,1
         do k1=klosten,khisten

          call CISBOX(xstenbox,nhalfbox_sten, &
            xlo,dx,i+i1,j+j1,k+k1, &
            bfact,level, &
            volsten,censten,SDIM)

          do im=1,num_materials
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
          enddo  ! im=1..num_materials

           ! sum of F_fluid=1
           ! sum of F_rigid<=1
          nhalf_box=1
          call make_vfrac_sum_ok_base( &
            cmofsten, &
            xstenbox,nhalfbox_sten,nhalf_box, &
            bfact,dx, &
            tessellate, & ! =0
            mofsten,SDIM,6)

          vfrac_fluid_sum=zero
          vfrac_solid_sum=zero
          im_raster_solid=0
          vfrac_raster_solid=zero

          do im=1,num_materials
           vofcomprecon=(im-1)*ngeom_recon+1
           vfrac_local(im)=mofsten(vofcomprecon)

           if (is_rigid(im).eq.0) then
            vfrac_fluid_sum=vfrac_fluid_sum+vfrac_local(im)
           else if (is_rigid(im).eq.1) then
            if (im_raster_solid.eq.0) then
             im_raster_solid=im
             vfrac_raster_solid=vfrac_local(im)
            else if ((im_raster_solid.ge.1).and. &
                     (im_raster_solid.le.num_materials).and. &
                     (is_rigid(im_raster_solid).eq.1)) then
             if (vfrac_raster_solid.lt.vfrac_local(im)) then
              im_raster_solid=im
              vfrac_raster_solid=vfrac_local(im)
             endif
            else
             print *,"im_raster_solid invalid"
             stop
            endif
     
            vfrac_solid_sum=vfrac_solid_sum+vfrac_local(im)
           else
            print *,"is_rigid(im) invalid"
            stop
           endif
          enddo ! im=1..num_materials

          if (abs(vfrac_fluid_sum-one).le.VOFTOL) then
           ! do nothing
          else
           print *,"vfrac_fluid_sum invalid"
           stop
          endif

          local_mod_cmofsten=0

          if (vfrac_solid_sum_center.ge.half) then
           ! do nothing, we can do the full cmof stencil in masked
           ! off is_rigid=1 cells
          else if (vfrac_solid_sum_center.lt.half) then

           if (vfrac_solid_sum.ge.half) then
            if ((i1.eq.0).and.(j1.eq.0).and.(k1.eq.0)) then
             print *,"expecting i1 or j1 or k1 not 0"
             stop
            endif

            if (partial_cmof_stencil_at_walls.eq.1) then
             cmofsten(D_DECL(i1,j1,k1))=0
             mod_cmofsten=1
             local_mod_cmofsten=1
            else if (partial_cmof_stencil_at_walls.eq.0) then
             ! do nothing
            else
             print *,"partial_cmof_stencil_at_walls invalid"
             stop
            endif
           else if (vfrac_solid_sum.lt.half) then
            ! do nothing
           else
            print *,"vfrac_solid_sum invalid"
            stop
           endif

          else
           print *,"vfrac_solid_sum_center invalid"
           stop
          endif

          if (local_mod_cmofsten.eq.0) then

           volume_super=volume_super+volsten
           do dir=1,SDIM
            cen_super(dir)=cen_super(dir)+volsten*censten(dir)
           enddo

           do im=1,num_materials

            vofcomprecon=(im-1)*ngeom_recon+1
            volmat=volsten*vfrac_local(im)
            vof_super(im)=vof_super(im)+volmat
            do dir=1,SDIM
             mofdata_super(vofcomprecon+dir)= &
                mofdata_super(vofcomprecon+dir)+ &
                volmat*(censten(dir)+mofsten(vofcomprecon+dir))
            enddo ! dir
            if (is_rigid(im).eq.0) then
             volume_super_mofdata=volume_super_mofdata+volmat
            else if (is_rigid(im).eq.1) then
             ! do nothing
            else
             print *,"is_rigid(im) invalid"
             stop
            endif

           enddo ! im=1..num_materials

          else if (local_mod_cmofsten.eq.1) then
           ! do nothing
          else
           print *,"local_mod_cmofsten invalid"
           stop
          endif

         enddo
         enddo
         enddo ! i1,j1,k1

         if (volume_super.gt.zero) then
          ! do nothing
         else
          print *,"volume_super invalid"
          stop
         endif
         if (volume_super_mofdata.gt.zero) then
          ! do nothing
         else
          print *,"volume_super_mofdata invalid"
          stop
         endif
         do dir=1,SDIM
          cen_super(dir)=cen_super(dir)/volume_super
         enddo

         do im=1,num_materials
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

           ! always standard MOF centroid for the rigid materials.
          if (is_rigid(im).eq.1) then
           do dir=1,SDIM
            mofdata_super(vofcomprecon+dir)=mofdata(vofcomprecon+dir)
           enddo
          else if (is_rigid(im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid invalid PLIC_3D.F90"
           stop
          endif

         enddo ! im=1..num_materials

         if (1.eq.0) then
          if (mod_cmofsten.eq.1) then
           print *,"mod_cmofsten=1:level,finest_level ", &
                   level,finest_level
           print *,"mod_cmofsten=1:i,j,k ",i,j,k
           print *,"x,y,z ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
           do i1=-1,1
           do j1=-1,1
           do k1=klosten,khisten
            print *,"i1,j1,k1,cmofsten(i1,j1,k1) ",i1,j1,k1, &
                    cmofsten(D_DECL(i1,j1,k1))
           enddo
           enddo
           enddo
          else if (mod_cmofsten.eq.0) then
           ! do nothing
          else
           print *,"mod_cmofsten invalid"
           stop
          endif 
         endif
        else if (continuous_mof_parm.eq.0) then
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
          cmofsten, &
          SDIM, &
          2)

        if (continuous_mof_parm.eq.2) then
           ! center cell centroids.
         do im=1,num_materials
          vofcomprecon=(im-1)*ngeom_recon+1
          do dir=1,SDIM
           mofdata_super(vofcomprecon+dir)=mofdata(vofcomprecon+dir)
          enddo
         enddo ! im
        else if (continuous_mof_parm.eq.0) then
         ! do nothing
        else
         print *,"continuous_mof_parm invalid"
         stop
        endif

        do im=1,num_materials
         total_calls(im)=total_calls(im)+mof_calls(tid+1,im)
         total_iterations(im)= &
          total_iterations(im)+mof_iterations(tid+1,im)
        enddo  ! im

       else
        print *,"level invalid slope recon"
        stop
       endif

       do dir=1,num_materials*ngeom_recon
        slopes(D_DECL(i,j,k),dir)=mofdata_super(dir)
       enddo

       if (update_flag.eq.1) then

#if (STANDALONE==0)
        if (use_ls_data.ne.1) then
         print *,"use_ls_data invalid"
         stop
        endif

        if ((level.ge.0).and.(level.le.finest_level)) then

         do i1=-1,1
         do j1=-1,1
         do k1=klosten,khisten
          local_mask=NINT(masknbr(D_DECL(i+i1,j+j1,k+k1),1))
          if (local_mask.eq.1) then ! fine-fine ghost in domain or interior.
           ! do nothing
          else if (local_mask.eq.0) then ! ghost value is low order accurate.
           ! do nothing
          else
           print *,"local_mask invalid"
           stop
          endif
         enddo 
         enddo 
         enddo 

         call calc_error_indicator( &
          level,max_level, &
          xsten,nhalf,dx,bfact, &
          voflist_center, &
          LS_stencil, &
          err,time)
         errsave=snew(D_DECL(i,j,k),num_materials*ngeom_raw+1)
         if (errsave.lt.err) then
          snew(D_DECL(i,j,k),num_materials*ngeom_raw+1)=err
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
        print *,"update_flag invalid for stand alone fort_sloperecon"
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
       print *,"AFTER fort_sloperecon --------------------------------"
       print *,"grid,level,finest ",gridno,level,finest_level
       print *,"STEP,TIME ",nsteps,time
      endif

      return
      end subroutine fort_sloperecon

      end module plic_cpp_module

#undef STANDALONE

