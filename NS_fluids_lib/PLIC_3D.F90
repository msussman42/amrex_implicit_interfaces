#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "EXTRAP_COMP.H"

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
      use amrex_fort_module, only : amrex_real
      contains

! mask=0 at coarse/fine border cells and physical boundaries.
! mask=1 at fine/fine and periodic border cells.
! mask=1 at symmetric border cells
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
        vofbc, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        xlo,dx, &
        maskcov,DIMS(maskcov), &
        masknbr,DIMS(masknbr), &
        snew,DIMS(snew), &
        vof,DIMS(vof), &
        LS,DIMS(LS), &
        slopes,DIMS(slopes), &
        nsteps, &
        time, &
        update_flag, & !RECON_UPDATE_ ...
        number_centroid_per_core, &
        delta_centroid_per_core, &
        total_calls, &
        total_iterations, &
        total_errors, &
        continuous_mof, &
        partial_cmof_stencil_at_walls) &
      bind(c,name='fort_sloperecon')

      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: gridno
      integer, INTENT(in) :: level,finest_level,max_level
      integer, INTENT(in) :: nsteps

      integer, INTENT(in) :: continuous_mof
      integer, INTENT(in) :: partial_cmof_stencil_at_walls
      integer, INTENT(in) :: update_flag
      integer, INTENT(out) :: number_centroid_per_core
      real(amrex_real), INTENT(out) :: delta_centroid_per_core
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: vofbc(SDIM,2)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: bfact
      integer, parameter :: ngrow=1
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(masknbr)
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(vof)
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(slopes)
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
    
      real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov)) 
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: masknbr(DIMV(masknbr),4) 
      real(amrex_real), pointer :: masknbr_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
        vof(DIMV(vof),num_materials*ngeom_raw) 
      real(amrex_real), pointer :: vof_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: LS(DIMV(LS),num_materials) 
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              slopes(DIMV(slopes),num_materials*ngeom_recon) 
      real(amrex_real), pointer :: slopes_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              snew(DIMV(snew),num_materials*ngeom_raw+1) 
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)
      
      integer :: i,j,k
      integer :: iside,jside,kside
      integer :: dir,side
      integer :: dir_local
      integer :: imask
      integer :: igridlo(3),igridhi(3)

      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))

      integer :: grid_index(SDIM)
      integer :: grid_side(SDIM)
      integer :: grid_level=-1

      integer im
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdata_super(num_materials*ngeom_recon)
      real(amrex_real) mofdata_super_vfrac(num_materials*ngeom_recon)
      real(amrex_real) vof_super(num_materials)
      real(amrex_real) mofsten(num_materials*ngeom_recon)
      real(amrex_real) multi_centroidA(num_materials,SDIM)

      integer vofcomprecon
      integer vofcompraw

      real(amrex_real) err,errsave
      integer local_mask

      real(amrex_real) orderflag
      integer total_calls(num_materials)
      integer total_iterations(num_materials)
      real(amrex_real) total_errors(num_materials)
      integer i1,j1,k1
      real(amrex_real) LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)

      integer, PARAMETER :: nmax=POLYGON_LIST_MAX !in: fort_sloperecon

      integer continuous_mof_parm
      integer continuous_mof_parm_super
     
      integer klosten,khisten
      integer, parameter :: nhalf=3

      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstencil_point(SDIM)

      real(amrex_real) xstenbox(-1:1,SDIM)
      integer num_fluid_materials_in_cell
      integer num_fluid_materials_in_stencil
      real(amrex_real) volume_super
      real(amrex_real) volume_super_mofdata
      real(amrex_real) volsten
      real(amrex_real) volmat
      real(amrex_real) censten(SDIM)
      real(amrex_real) cen_super(SDIM)
      real(amrex_real) voflist_center(num_materials)
      real(amrex_real) voflist_stencil(num_materials)
      real(amrex_real) voflist_test
      integer mof_verbose
      integer, parameter :: use_ls_data=1
      integer, parameter :: nhalfbox_sten=1
      real(amrex_real) dxmaxLS
      integer debugslope
      integer, parameter :: tessellate=0
      integer, parameter :: shapeflag=0
      integer, parameter :: continuous_mof_standard=STANDARD_MOF

      real(amrex_real) vfrac_fluid_sum
      real(amrex_real) vfrac_solid_sum
      real(amrex_real) vfrac_solid_sum_center
      real(amrex_real) vfrac_raster_solid
      real(amrex_real) vfrac_local(num_materials)
      integer im_raster_solid
      integer mod_cmofsten
      integer :: interior_flag,outside_fab

      real(amrex_real) :: xtet(SDIM+1,SDIM)
      real(amrex_real) :: cmof_centroid
      real(amrex_real) :: delta_centroid
      real(amrex_real) :: multi_area(num_materials)
      real(amrex_real) :: multi_volume(num_materials)
      real(amrex_real) :: multi_cen(SDIM,num_materials)

      integer :: local_maskcov
      integer :: fluid_obscured
      real(amrex_real) LS_clamped
      real(amrex_real) VEL_clamped(SDIM)
      real(amrex_real) temperature_clamped
      integer :: prescribed_flag

#include "mofdata.H"

      maskcov_ptr=>maskcov
      masknbr_ptr=>masknbr
      vof_ptr=>vof
      LS_ptr=>LS

      slopes_ptr=>slopes
      snew_ptr=>snew

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      debugslope=0

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
       print *,"gridno=",gridno
       print *,"level=",level
       print *,"finest_level=",finest_level
       print *,"max_level=",max_level
       stop
      endif
      if (max_level.lt.finest_level) then
       print *,"max_level invalid sloperecon: ",max_level
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

      if (ngrow.ne.1) then
       print *,"ngrow invalid in fort_sloperecon: ",ngrow
       stop
      endif
      if ((continuous_mof.eq.STANDARD_MOF).or. & !MOF
          (continuous_mof.eq.CMOF_X)) then !CMOF
       ! do nothing
      else
       print *,"continuous_mof invalid (fort_sloperecon): ",continuous_mof
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

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"levelrz invalid fort_sloperecon"
       stop
      endif

      if ((update_flag.eq.RECON_UPDATE_NULL).or. &
          (update_flag.eq.RECON_UPDATE_STATE_ERR).or. &
          (update_flag.eq.RECON_UPDATE_STATE_CENTROID).or. &
          (update_flag.eq.RECON_UPDATE_STATE_ERR_AND_CENTROID)) then
       ! do nothing
      else
       print *,"update_flag invalid1: ",update_flag
       stop
      endif

      number_centroid_per_core=0
      delta_centroid_per_core=zero

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,masknbr_ptr,1,-1)
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,vof_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      call checkbound_array(fablo,fabhi,slopes_ptr,ngrow,-1)

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

      do k = igridlo(3),igridhi(3)
      do j = igridlo(2),igridhi(2)
      do i = igridlo(1),igridhi(1)

        ! we must visit all of the covered cells too since
        ! the AMR error indicator is needed on the coarse levels.
        ! Also, the covered cell reconstruction(s) are needed in order
        ! reinitialize coarse grid level set function.
       local_maskcov=NINT(maskcov(D_DECL(i,j,k)))
       if ((local_maskcov.eq.1).or. &
           (local_maskcov.eq.0)) then
        ! do nothing
       else
        print *,"local_maskcov invalid: ",local_maskcov
        stop
       endif

       grid_index(1)=i
       grid_index(2)=j
       if (SDIM.eq.3) then
        grid_index(SDIM)=k
       endif

       interior_flag=1

       do dir=1,SDIM
        do side=1,2

         do dir_local=1,SDIM
          grid_side(dir_local)=grid_index(dir_local)
         enddo

         outside_fab=0
      
         if (side.eq.1) then
          grid_side(dir)=grid_side(dir)-1
          if (grid_side(dir).lt.fablo(dir)) then
           outside_fab=1
          else if (grid_side(dir).lt.fabhi(dir)) then
           !do nothing
          else
           print *,"grid_side(dir) invalid"
           stop
          endif
         else if (side.eq.2) then
          grid_side(dir)=grid_side(dir)+1
          if (grid_side(dir).gt.fabhi(dir)) then
           outside_fab=1
          else if (grid_side(dir).gt.fablo(dir)) then
           !do nothing
          else
           print *,"grid_side(dir) invalid"
           stop
          endif
         else
          print *,"side invalid: ",side
          stop
         endif

         iside=grid_side(1)
         jside=grid_side(2)
         kside=grid_side(SDIM)
         imask=NINT(masknbr(D_DECL(iside,jside,kside),1))

         if (outside_fab.eq.1) then

          if (vofbc(dir,side).eq.REFLECT_EVEN) then
           interior_flag=0
           if (imask.eq.0) then
            ! do nothing
           else
            print *,"expecting imask=0"
            stop
           endif
          else if (vofbc(dir,side).eq.FOEXTRAP) then
           interior_flag=0
           if (imask.eq.0) then
            ! do nothing
           else
            print *,"expecting imask=0"
            stop
           endif
          else if (vofbc(dir,side).eq.EXT_DIR) then
           interior_flag=0
           if (imask.eq.0) then
            ! do nothing
           else
            print *,"expecting imask=0:",imask
            stop
           endif
          else if (vofbc(dir,side).eq.INT_DIR) then

           if (imask.eq.1) then
            ! do nothing
           else if (imask.eq.0) then
            interior_flag=0
           else
            print *,"imask invalid: ",imask
            stop
           endif

          else
           print *,"vofbc invalid; dir,side,vofbc: ", &
            dir,side,vofbc(dir,side)
           stop
          endif

         else if (outside_fab.eq.0) then
          if (imask.eq.1) then
           ! do nothing
          else
           print *,"expecting imask=1"
           stop
          endif
         else
          print *,"outside_fab invalid: ",outside_fab
          stop
         endif
        enddo !side=1,2
       enddo !dir=1,SDIM

       if ((update_flag.eq.RECON_UPDATE_NULL).or. &
           (update_flag.eq.RECON_UPDATE_STATE_ERR).or. &
           (update_flag.eq.RECON_UPDATE_STATE_CENTROID).or. &
           (update_flag.eq.RECON_UPDATE_STATE_ERR_AND_CENTROID)) then
        ! do nothing
       else
        print *,"update_flag invalid loop ",update_flag
        print *,"i,j,k ",i,j,k
        print *,"igridlo ",igridlo(1),igridlo(2),igridlo(3)
        print *,"igridhi ",igridhi(1),igridhi(2),igridhi(3)
        stop
       endif

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir=1,SDIM
        xstencil_point(dir)=xsten(0,dir)
       enddo

       do k1=klosten,khisten
       do j1=-1,1
       do i1=-1,1
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

         !vof,cenref,order,slope,intercept
        do dir=1,SDIM
         mofdata(vofcomprecon+SDIM+1+dir)= &
             slopes(D_DECL(i,j,k),vofcomprecon+SDIM+1+dir)
        enddo

        if ((mofdata(vofcomprecon).ge.-0.1d0).and. &
            (mofdata(vofcomprecon).le.1.1d0)) then
         ! do nothing
        else if ((mofdata(vofcomprecon).lt.-0.1d0).or. &
                 (mofdata(vofcomprecon).gt.1.1d0)) then
         print *,"mofdata(vofcomprecon) out of range"
         print *,"mofdata(vofcomprecon)=",mofdata(vofcomprecon)
         stop
        else
         print *,"mofdata(vofcomprecon) is NaN"
         print *,"mofdata(vofcomprecon)=",mofdata(vofcomprecon)
         stop
        endif

        orderflag=zero
        mofdata(vofcomprecon+SDIM+1)=orderflag

         !initialize the intercept to be zero
        mofdata(vofcomprecon+ngeom_recon-1)=zero

       enddo  ! im=1..num_materials

        ! sum of F_fluid=1
        ! sum of F_rigid<=1
       call make_vfrac_sum_ok_base( &
         cmofsten, &  ! INTENT(in)
         xsten, &
         nhalf, &
         continuous_mof_standard, &
         bfact,dx, &
         tessellate, & ! =0
         mofdata, &  ! INTENT(inout)
         SDIM)

       vfrac_fluid_sum=zero
       vfrac_solid_sum=zero
       im_raster_solid=0

        ! i,j,k
        ! vfrac_raster_solid=max_{is_rigid(im)==1} vfrac(im)
       vfrac_raster_solid=zero

       do im=1,num_materials
        vofcomprecon=(im-1)*ngeom_recon+1
        voflist_center(im)=mofdata(vofcomprecon)
        vof_super(im)=voflist_center(im)

         ! voflist_stencil(im)=max_{3x3x3 stencil} F(im,stencil)
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
          else if (vfrac_raster_solid.ge.voflist_center(im)) then
           ! do nothing
          else
           print *,"vfrac_raster_solid or voflist_center is NaN:", &
             im,vfrac_raster_solid,voflist_center(im)
           stop
          endif
         else
          print *,"im_raster_solid invalid:",im_raster_solid
          stop
         endif
      
         vfrac_solid_sum=vfrac_solid_sum+voflist_center(im)
        else
         print *,"is_rigid(im) invalid: ",im,is_rigid(im)
         stop
        endif

       enddo ! im=1..num_materials

       vfrac_solid_sum_center=vfrac_solid_sum

       if (abs(vfrac_fluid_sum-one).le.EPS1) then
        ! do nothing
       else
        print *,"vfrac_fluid_sum invalid: ",vfrac_fluid_sum
        stop
       endif

       if ((level.ge.0).and. &
           (level.le.finest_level)) then

        do k1=klosten,khisten
        do j1=-1,1
        do i1=-1,1
         do im=1,num_materials
          LS_stencil(D_DECL(i1,j1,k1),im)= &
            LS(D_DECL(i+i1,j+j1,k+k1),im)
          vofcompraw=(im-1)*ngeom_raw+1
          voflist_test=vof(D_DECL(i+i1,j+j1,k+k1),vofcompraw)
          if ((voflist_test.lt.-0.1d0).or. &
              (voflist_test.gt.1.1d0)) then
           print *,"voflist_test invalid"
           print *,"im,voflist_test= ",im,voflist_test
           print *,"i1,j1,k1 ",i1,j1,k1
           print *,"i,j,k ",i,j,k
           print *,"igridlo ",igridlo(1),igridlo(2),igridlo(3)
           print *,"igridhi ",igridhi(1),igridhi(2),igridhi(3)
           stop
          else if ((voflist_test.ge.-0.1d0).and. &
                   (voflist_test.le.1.1d0)) then
            ! voflist_stencil(im)=max_{3x3x3 stencil} F(im,stencil)
           if (voflist_test.gt.voflist_stencil(im)) then
            voflist_stencil(im)=voflist_test
           endif
          else
           print *,"voflist_test bust"
           stop
          endif
         enddo ! im=1..num_materials
        enddo
        enddo
        enddo  ! i1,j1,k1  = -1,1

        num_fluid_materials_in_cell=0
        num_fluid_materials_in_stencil=0
        do im=1,num_materials
         if (is_rigid(im).eq.0) then
           ! voflist_stencil(im)=max_{3x3x3 stencil} F(im,stencil)
          if (voflist_stencil(im).gt.VOFTOL) then
           num_fluid_materials_in_stencil=num_fluid_materials_in_stencil+1
          endif
          if (voflist_center(im).gt.VOFTOL) then
           num_fluid_materials_in_cell=num_fluid_materials_in_cell+1
          endif
         else if (is_rigid(im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid PLIC_3D.F90"
          stop
         endif
        enddo ! im=1..num_materials
 
        if (num_fluid_materials_in_cell.gt.num_fluid_materials_in_stencil) then
         print *,"num_fluid_materials_in_cell invalid"
         print *,"num_fluid_materials_in_cell: ",num_fluid_materials_in_cell
         print *,"num_fluid_materials_in_stencil: ", &
           num_fluid_materials_in_stencil
         stop
        endif

        call SUB_clamped_LS(xstencil_point,time,LS_clamped, &
          VEL_clamped,temperature_clamped,prescribed_flag,dx)

        fluid_obscured=0
        if (vfrac_solid_sum.ge.half) then
         fluid_obscured=1
        else if (vfrac_solid_sum.le.half) then
         !do nothing
        else
         print *,"vfrac_solid_sum invalid:",vfrac_solid_sum
         stop
        endif
        if (LS_clamped.ge.zero) then
         fluid_obscured=1
        else if (LS_clamped.le.zero) then
         !do nothing
        else
         print *,"LS_clamped invalid:",LS_clamped
         stop
        endif 

        do im=1,num_materials*ngeom_recon
         mofdata_super(im)=mofdata(im)
        enddo

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volume_super,cen_super,SDIM)

        if ((level.lt.max_level).or. &
            (level.ne.decision_tree_max_level).or. &
            (fluid_obscured.eq.1)) then

         ! always use MOF on the coarser levels or if decision tree data
         ! is not available or if the fluids in a cell are ``mostly''
         ! obscured by "is_rigid" materials.
         continuous_mof_parm=STANDARD_MOF

        else if ((level.eq.max_level).and. &
                 (level.eq.decision_tree_max_level).and. &
                 (fluid_obscured.eq.0)) then

         if (num_fluid_materials_in_cell.eq.1) then
          continuous_mof_parm=STANDARD_MOF
         else if ((num_fluid_materials_in_cell.ge.2).and. &
                  (num_fluid_materials_in_cell.le.num_materials)) then

          continuous_mof_parm=continuous_mof

         else
          print *,"num_fluid_materials_in_cell invalid: ", &
           num_fluid_materials_in_cell
          stop
         endif

        else
         print *,"level or decision_tree_max_level or fluid_obscured invalid"
         print *,"level=",level
         print *,"max_level=",max_level
         print *,"decision_tree_max_level=",decision_tree_max_level
         print *,"fluid_obscured=",fluid_obscured
         stop
        endif

         ! supercell for centroid cost function.
         ! center cell for volume constraint.
        if (continuous_mof_parm.eq.CMOF_X) then
  
         volume_super=zero ! volume of the extended region
         volume_super_mofdata=zero !same as volume_super, except by im_fluid.

         do dir=1,SDIM
          cen_super(dir)=zero
         enddo

         do im=1,num_materials

          if (is_rigid(im).eq.1) then
           ! do nothing
          else if (is_rigid(im).eq.0) then
           vofcomprecon=(im-1)*ngeom_recon+1
           do dir=1,SDIM
            mofdata_super(vofcomprecon+dir)=zero
           enddo
           vof_super(im)=zero
          else
           print *,"is_rigid(im) invalid"
           stop
          endif

         enddo ! im=1..num_materials

         do k1=klosten,khisten
         do j1=-1,1
         do i1=-1,1

          call CISBOX(xstenbox, &
            nhalfbox_sten, & ! =1
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

           !vof,cenref,order,slope,intercept
           do dir=1,SDIM
            mofsten(vofcomprecon+SDIM+1+dir)= &
             slopes(D_DECL(i+i1,j+j1,k+k1),vofcomprecon+SDIM+1+dir)
           enddo

           orderflag=zero
           mofsten(vofcomprecon+SDIM+1)=orderflag

           !initialize the intercept to be zero
           mofsten(vofcomprecon+ngeom_recon-1)=zero

          enddo  ! im=1..num_materials

           ! sum of F_fluid=1
           ! sum of F_rigid<=1
          call make_vfrac_sum_ok_base( &
            cmofsten, & !intent(in)
            xstenbox, &
            nhalfbox_sten, & ! =1
            continuous_mof_standard, & ! =STANDARD_MOF
            bfact,dx, &
            tessellate, & ! =0
            mofsten, &
            SDIM)

          vfrac_fluid_sum=zero
          vfrac_solid_sum=zero
          im_raster_solid=0
           ! i+i',j+j',k+k'
           ! vfrac_raster_solid=max_{is_rigid(im)==1} vfrac(im)
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
             print *,"im_raster_solid invalid:",im_raster_solid
             stop
            endif
     
            vfrac_solid_sum=vfrac_solid_sum+vfrac_local(im)
           else
            print *,"is_rigid(im) invalid"
            stop
           endif
          enddo ! im=1..num_materials

          if (abs(vfrac_fluid_sum-one).le.EPS1) then
           ! do nothing
          else
           print *,"vfrac_fluid_sum invalid: ",vfrac_fluid_sum
           stop
          endif

          mod_cmofsten=0

          if (vfrac_solid_sum_center.ge.half) then

           ! do nothing, we can do the full cmof stencil in masked
           ! off is_rigid=1 cells (for reconstructing the fluid
           ! interfaces)

          else if (vfrac_solid_sum_center.lt.half) then

           if (vfrac_solid_sum.ge.half) then

            if ((i1.eq.0).and.(j1.eq.0).and.(k1.eq.0)) then
             print *,"expecting i1 or j1 or k1 not 0"
             stop
            endif

            if (partial_cmof_stencil_at_walls.eq.1) then
             cmofsten(D_DECL(i1,j1,k1))=0
             mod_cmofsten=1
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

          if (mod_cmofsten.eq.0) then

           volume_super=volume_super+volsten
           do dir=1,SDIM
            cen_super(dir)=cen_super(dir)+volsten*censten(dir)
           enddo

           do im=1,num_materials

            if (is_rigid(im).eq.0) then

             vofcomprecon=(im-1)*ngeom_recon+1
             volmat=volsten*vfrac_local(im)
             vof_super(im)=vof_super(im)+volmat
             do dir=1,SDIM
              mofdata_super(vofcomprecon+dir)= &
                mofdata_super(vofcomprecon+dir)+ &
                volmat*(censten(dir)+mofsten(vofcomprecon+dir))
             enddo ! dir
             volume_super_mofdata=volume_super_mofdata+volmat

            else if (is_rigid(im).eq.1) then
             ! do nothing
            else
             print *,"is_rigid(im) invalid"
             stop
            endif

           enddo ! im=1..num_materials

          else if (mod_cmofsten.eq.1) then
           ! do nothing
          else
           print *,"mod_cmofsten invalid"
           stop
          endif

         enddo
         enddo
         enddo ! i1,j1,k1

         if (volume_super.gt.zero) then
          ! do nothing
         else
          print *,"volume_super invalid: ",volume_super
          stop
         endif

         if (volume_super_mofdata.gt.zero) then
          ! do nothing
         else
          print *,"volume_super_mofdata invalid: ",volume_super_mofdata
          stop
         endif

         do dir=1,SDIM
          cen_super(dir)=cen_super(dir)/volume_super
         enddo

         do im=1,num_materials
          vofcomprecon=(im-1)*ngeom_recon+1

           ! always standard MOF centroid for the rigid materials.
          if (is_rigid(im).eq.1) then

           do dir=1,SDIM
            mofdata_super(vofcomprecon+dir)=mofdata(vofcomprecon+dir)
           enddo

          else if (is_rigid(im).eq.0) then

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
            print *,"vof_super(im) invalid: ",im,vof_super(im)
            stop
           endif

          else
           print *,"is_rigid invalid PLIC_3D.F90"
           stop
          endif

         enddo ! im=1..num_materials

        else if (continuous_mof_parm.eq.STANDARD_MOF) then
         ! do nothing
        else
         print *,"continuous_mof_parm invalid: ",continuous_mof_parm
         stop
        endif

        mof_verbose=0

        grid_level=-1

        if (interior_flag.eq.1) then

         if ((level.eq.training_max_level).or. &
             (level.eq.decision_tree_max_level))  then
          if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
              (levelrz.eq.COORDSYS_RZ)) then
           grid_level=level
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           grid_level=level
          else
           print *,"levelrz invalid: ",levelrz
           stop
          endif
         else if ((level.ge.0).and. &
                  (level.le.finest_level)) then
          ! do nothing
         else
          print *,"level invalid: ",level
          print *,"finest_level: ",finest_level
          print *,"max_level: ",max_level
          print *,"training_max_level: ",training_max_level
          print *,"decision_tree_max_level: ",decision_tree_max_level
          stop
         endif

        else if (interior_flag.eq.0) then
         ! do nothing
        else
         print *,"interior_flag invalid: ",interior_flag
         stop
        endif
 
        if ((continuous_mof_parm.eq.STANDARD_MOF).or. &
            (num_fluid_materials_in_stencil.eq.1).or. &
            (num_fluid_materials_in_stencil.eq.2).or. &
            ((num_fluid_materials_in_stencil.ge.3).and. &
             (num_fluid_materials_in_cell.eq. &
              num_fluid_materials_in_stencil))) then

         call multimaterial_MOF( &
          bfact,dx, &
          xsten, &
          nhalf, &
          mof_verbose, &
          use_ls_data, & ! use_ls_data=1
          LS_stencil, &
          geom_xtetlist(1,1,1,tid+1), &
          geom_xtetlist_old(1,1,1,tid+1), &
          nmax, &
          nmax, &
          mofdata_super, & !intent(inout)
          vof_super, &
          multi_centroidA, & ! (num_materials,sdim) relative to supercell
          continuous_mof_parm, &
          cmofsten, & !intent(in)
          grid_index, &
          grid_level, &
          SDIM)

         if (continuous_mof_parm.eq.CMOF_X) then
           ! center cell centroids.
          do im=1,num_materials
           vofcomprecon=(im-1)*ngeom_recon+1
           do dir=1,SDIM
            mofdata_super(vofcomprecon+dir)=mofdata(vofcomprecon+dir)
           enddo
          enddo ! im
         else if (continuous_mof_parm.eq.STANDARD_MOF) then
          ! do nothing
         else
          print *,"continuous_mof_parm invalid(2) ",continuous_mof_parm
          stop
         endif

         ! mof_calls, mof_iterations, mof_errors are init. in 
         ! multimaterial_MOF
         do im=1,num_materials
          total_calls(im)=total_calls(im)+mof_calls(tid+1,im)
          total_iterations(im)= &
           total_iterations(im)+mof_iterations(tid+1,im)
          total_errors(im)= &
           total_errors(im)+mof_errors(tid+1,im)
         enddo  ! im=1..num_materials

         ! for general CMOF reconstruction (>=3 materials):
         ! 1. find CMOF reconstruction using both CMOF centroids and
         !    CMOF volumes
         ! 2. Find the centroids of the resulting reconstruction in the
         !    center cell.
         ! 3. If both the original center cell volume fractions and
         !    center cell reconstructed volume fraction satisfy 0<F<1,
         !    then replace center cell centroids with those from step 3.
         ! 4. Do a standard MOF reconstruction.
         ! Future: if same number of materials 0<F<1 in center cell as in
         !  the CMOF stencil, just do one step:
         ! 1. find CMOF reconstruction using CMOF centroids and MOF 
         !    (cell center) volumes.
        else if ((continuous_mof_parm.eq.CMOF_X).and. &
                 (num_fluid_materials_in_stencil.ge.3).and. &
                 (num_fluid_materials_in_stencil.gt. &
                  num_fluid_materials_in_cell)) then

         continuous_mof_parm_super=CMOF_F_AND_X

         do dir=1,num_materials*ngeom_recon
          mofdata_super_vfrac(dir)=mofdata_super(dir)
         enddo
         do im=1,num_materials
          vofcomprecon=(im-1)*ngeom_recon+1
          mofdata_super_vfrac(vofcomprecon)=vof_super(im)
         enddo

          !mofdata_super: centroid super cell, vfrac regular cell
          !mofdata_super_vfrac: centroid super cell, vfrac super cell

         call multimaterial_MOF( &
          bfact,dx, &
          xsten, &
          nhalf, &
          mof_verbose, &
          use_ls_data, & ! use_ls_data=1
          LS_stencil, &
          geom_xtetlist(1,1,1,tid+1), &
          geom_xtetlist_old(1,1,1,tid+1), &
          nmax, &
          nmax, &
          mofdata_super_vfrac, & !intent(inout)
          vof_super, &
          multi_centroidA, & ! (num_materials,sdim) relative to supercell
          continuous_mof_parm_super, &
          cmofsten, & !intent(in)
          grid_index, &
          grid_level, &
          SDIM)

         ! mof_calls, mof_iterations, mof_errors are init. in 
         ! multimaterial_MOF
         do im=1,num_materials
          total_calls(im)=total_calls(im)+mof_calls(tid+1,im)
          total_iterations(im)= &
           total_iterations(im)+mof_iterations(tid+1,im)
          total_errors(im)= &
           total_errors(im)+mof_errors(tid+1,im)
         enddo  ! im=1..num_materials

         do im=1,num_materials
          vofcomprecon=(im-1)*ngeom_recon+1
          do dir=1,SDIM
           mofdata_super(vofcomprecon+dir)=mofdata(vofcomprecon+dir)
          enddo
           ! order,slope,intercept (SDIM+2)
          do dir=SDIM+1,ngeom_recon-1
           mofdata_super(vofcomprecon+dir)= &
               mofdata_super_vfrac(vofcomprecon+dir)
          enddo
         enddo ! im=1..num_materials

          !mofdata_super: centroid regular cell, vfrac regular cell,
          !               reconstruction super cell.

         call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volume_super,cen_super,SDIM)

         call multi_get_volume_grid( &
          EPS_8_4, &
          tessellate, & ! =0
          bfact,dx, &
          xsten,nhalf, & ! phi = n dot (x-x0) + intercept
          mofdata_super, &
          xsten,nhalf, & ! find volumes within xsten (cell i,j,k)
          xtet, &        ! not within xtet
          multi_volume, &
          multi_cen, & !(sdim,num_materials) absolute frame of ref.
          multi_area, & !(num_materials)
          geom_xtetlist_uncapt(1,1,1,tid+1), &
          nmax, &
          nmax, &
          SDIM, &
          shapeflag) !shapeflag=0

         vfrac_fluid_sum=zero
         do im=1,num_materials
          if (is_rigid(im).eq.0) then
           vfrac_fluid_sum=vfrac_fluid_sum+multi_volume(im)
          else if (is_rigid(im).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif
         enddo ! im=1,..,num_materials

         if (vfrac_fluid_sum.gt.zero) then
          do im=1,num_materials
           if (is_rigid(im).eq.0) then
            multi_volume(im)=multi_volume(im)/vfrac_fluid_sum
           else if (is_rigid(im).eq.1) then
            ! do nothing
           else
            print *,"is_rigid invalid"
            stop
           endif
          enddo ! im=1,..,num_materials
         else
          print *,"vfrac_fluid_sum invalid"
          stop
         endif

         do im=1,num_materials
          if (is_rigid(im).eq.0) then

           vofcomprecon=(im-1)*ngeom_recon+1
           vfrac_local(im)=mofdata_super(vofcomprecon)
           vof_super(im)=vfrac_local(im)

           if ((vfrac_local(im).ge.0.01D0).and. &
               (vfrac_local(im).le.0.99D0).and. &
               (multi_volume(im).ge.0.01D0).and. &
               (multi_volume(im).le.0.99D0)) then

            do dir=1,SDIM
             cmof_centroid=multi_cen(dir,im)-cen_super(dir)
             mofdata_super(vofcomprecon+dir)=cmof_centroid
            enddo
           else if (abs(vfrac_local(im)).le.0.01D0) then
            ! do nothing
           else if (abs(one-vfrac_local(im)).le.0.01D0) then
            ! do nothing
           else if (abs(multi_volume(im)).le.0.01D0) then
            ! do nothing
           else if (abs(one-multi_volume(im)).le.0.01D0) then
            ! do nothing
           else
            print *,"vfrac_local or multi_volumeinvalid"
            stop
           endif

          else if (is_rigid(im).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif
         enddo ! im=1,num_materials

          ! replace the slopes with the particle container derived slopes.
         do im=1,num_materials
          vofcomprecon=(im-1)*ngeom_recon+1
          !vof,cenref,order,slope,intercept
          do dir=1,SDIM
           mofdata_super(vofcomprecon+SDIM+1+dir)= &
             slopes(D_DECL(i,j,k),vofcomprecon+SDIM+1+dir)
          enddo
         enddo ! im=1,num_materials

         continuous_mof_parm=STANDARD_MOF

         call multimaterial_MOF( &
          bfact,dx, &
          xsten, &
          nhalf, &
          mof_verbose, &
          use_ls_data, & ! use_ls_data=1
          LS_stencil, &
          geom_xtetlist(1,1,1,tid+1), &
          geom_xtetlist_old(1,1,1,tid+1), &
          nmax, &
          nmax, &
          mofdata_super, & !intent(inout)
          vof_super, &
          multi_centroidA, & ! (num_materials,sdim) relative to supercell
          continuous_mof_parm, &
          cmofsten, & !intent(in)
          grid_index, &
          grid_level, &
          SDIM)

         ! mof_calls, mof_iterations, mof_errors are init. in 
         ! multimaterial_MOF
         do im=1,num_materials
          total_calls(im)=total_calls(im)+mof_calls(tid+1,im)
          total_iterations(im)= &
           total_iterations(im)+mof_iterations(tid+1,im)
          total_errors(im)= &
           total_errors(im)+mof_errors(tid+1,im)
         enddo  ! im=1..num_materials
 
        else
         print *,"continuous_mof_parm or num_fluid_materials bad"
         print *,"continuous_mof_parm=",continuous_mof_parm
         print *,"num_fluid_materials_in_cell=",num_fluid_materials_in_cell
         print *,"num_fluid_materials_in_stencil=", &
             num_fluid_materials_in_stencil
         stop
        endif

       else
        print *,"level invalid fort_sloperecon: ",level
        stop
       endif

       if ((update_flag.eq.RECON_UPDATE_STATE_CENTROID).or. &
           (update_flag.eq.RECON_UPDATE_STATE_ERR_AND_CENTROID)) then

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volume_super,cen_super,SDIM)

        call multi_get_volume_grid( &
         EPS_8_4, &
         tessellate, & ! =0
         bfact,dx, &
         xsten,nhalf, & ! phi = n dot (x-x0) + intercept
         mofdata_super, &
         xsten,nhalf, & ! find volumes within xsten
         xtet, &        ! not within xtet
         multi_volume, &
         multi_cen, & !(sdim,num_materials) absolute frame of ref.
         multi_area, & !(num_materials)
         geom_xtetlist_uncapt(1,1,1,tid+1), &
         nmax, &
         nmax, &
         SDIM, &
         shapeflag)

        do im=1,num_materials
         if (is_rigid(im).eq.0) then

          if (local_maskcov.eq.0) then
           ! do nothing
          else if (local_maskcov.eq.1) then
           
           vofcomprecon=(im-1)*ngeom_recon+1
           vofcompraw=(im-1)*ngeom_raw+1

           vfrac_local(im)=mofdata_super(vofcomprecon)

           if ((vfrac_local(im).ge.0.01D0).and. &
               (vfrac_local(im).le.0.99D0)) then

            number_centroid_per_core=number_centroid_per_core+1
            delta_centroid=zero
            do dir=1,SDIM
             cmof_centroid=multi_cen(dir,im)-cen_super(dir)
             delta_centroid=delta_centroid+ &
              ((mofdata_super(vofcomprecon+dir)-cmof_centroid)/dx(dir))**2
             mofdata_super(vofcomprecon+dir)=cmof_centroid
             snew(D_DECL(i,j,k),vofcompraw+dir)=cmof_centroid
            enddo
            delta_centroid=sqrt(delta_centroid)
            delta_centroid_per_core=delta_centroid_per_core+delta_centroid
           else if (abs(vfrac_local(im)).le.0.01D0) then
            ! do nothing
           else if (abs(one-vfrac_local(im)).le.0.01D0) then
            ! do nothing
           else
            print *,"vfrac_local invalid"
            stop
           endif

          else
           print *,"local_maskcov invalid: ",local_maskcov
           stop
          endif

         else if (is_rigid(im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif
        enddo ! im=1,num_materials

       else if (update_flag.eq.RECON_UPDATE_NULL) then
        ! do nothing
       else if (update_flag.eq.RECON_UPDATE_STATE_ERR) then
        ! do nothing
       else 
        print *,"update_flag invalid2: ",update_flag
        stop
       endif

       do dir=1,num_materials*ngeom_recon
        slopes(D_DECL(i,j,k),dir)=mofdata_super(dir)
       enddo
       
       if ((update_flag.eq.RECON_UPDATE_STATE_ERR).or. &
           (update_flag.eq.RECON_UPDATE_STATE_ERR_AND_CENTROID)) then

        if ((level.ge.0).and.(level.le.finest_level)) then

         do k1=klosten,khisten
         do j1=-1,1
         do i1=-1,1
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

          !calc_error_indicator is declared in PROB.F90
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
         if ((errsave.ge.zero).and. &
             (err.ge.zero)) then 
          !do nothing
         else
          print *,"err bust"
          stop
         endif

        else
         print *,"level invalid fort_sloperecon 2: ",level
         stop
        endif
       else if (update_flag.eq.RECON_UPDATE_NULL) then
        ! do nothing
       else if (update_flag.eq.RECON_UPDATE_STATE_CENTROID) then
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


      subroutine fort_MOF_training( &
        num_samples, &
        op_training, & !0=alloc 1=create data,python proc 2=read network data
        cpp_training_lo, &
        cpp_training_hi, &
        i,j,k, &
        finest_level, &
        max_level, &
        bfact, &
        domlo,domhi, &
        dx, &
        continuous_mof) &
      bind(c,name='fort_MOF_training')

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: num_samples
      integer, INTENT(in) :: op_training
      integer, INTENT(inout) :: cpp_training_lo(SDIM)
      integer, INTENT(inout) :: cpp_training_hi(SDIM)
      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: max_level
      integer, INTENT(in) :: continuous_mof !=STANDARD_MOF or CMOF_X
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(SDIM)

      integer :: nmax=POLYGON_LIST_MAX

      real(amrex_real) :: vof_training(num_samples)
      real(amrex_real) :: phi_training(num_samples)
      real(amrex_real) :: theta_training(num_samples)
       ! centroid, VOF, angle_exact, angle_init
      real(amrex_real) :: data_training(NTRAINING,num_samples)
      real(amrex_real) :: xc0(SDIM)

      real(amrex_real) :: angle_exact_db(SDIM-1)
      real(amrex_real) :: angle_init_db(SDIM-1)
      real(amrex_real) :: angle_and_vfrac(SDIM)
      real(amrex_real) :: refvfrac(1)
      real(amrex_real) :: vof_single 
      real(amrex_real) :: refcen(SDIM)
      integer :: try_new_vfrac

      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      integer dir
      integer i1,j1,k1
      integer i_training

      real(amrex_real) :: npredict(3,SDIM)
      real(amrex_real) :: vof_null
      real(amrex_real) :: centroid_null(SDIM)
      real(amrex_real) :: mag_centroid(3)
      integer :: critical_material

      real(amrex_real) :: DT_cost,NN_cost,RF_cost
      real(amrex_real) :: angle_exact_db_data(SDIM-1)

      integer sysret
      integer cmof_idx
      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))
      integer klosten,khisten
      integer tid
  
      tid=0

      if (num_samples.gt.0) then
       ! do nothing
      else
       print *,"num_samples invalid: ",num_samples
       stop
      endif

      if (max_level.eq.fort_max_level) then
       training_max_level=max_level
      else
       print *,"max_level and fort_max_level mismatch"
       stop
      endif

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
      do k1=klosten,khisten
      do j1=-1,1
      do i1=-1,1
       cmofsten(D_DECL(i1,j1,k1))=1
      enddo
      enddo
      enddo

      if (bfact.lt.1) then
       print *,"bfact invalid170"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((continuous_mof.eq.STANDARD_MOF).or. & ! MOF
          (continuous_mof.eq.CMOF_X)) then ! CMOF
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

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"levelrz invalid fort_MOF_training"
       stop
      endif

      if (op_training.eq.0) then

       do dir=1,SDIM
        training_lo(dir)=0
        training_hi(dir)=0
       enddo 

       do dir=1,SDIM
        if (domlo(dir).eq.0) then
         ! do nothing
        else
         print *,"expecting domlo=0"
         stop
        endif
        training_hi(dir)=domlo(dir)+bfact-1
        if (training_hi(dir).le.domhi(dir)) then
         ! do nothing 
        else
         print *,"training_hi invalid"
         stop
        endif
       enddo ! dir=1..sdim

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        dir=1
        training_hi(dir)=domhi(dir)
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        dir=1
        training_hi(dir)=domhi(dir)
        dir=2
        training_hi(dir)=domhi(dir)
       else
        print *,"levelrz invalid fort_MOF_training"
        stop
       endif

       allocate(training_array( &
         D_DECL(training_lo(1):training_hi(1),training_lo(2):training_hi(2),training_lo(SDIM):training_hi(SDIM) ), &
                  0:1))

       do dir=1,SDIM
        cpp_training_lo(dir)=training_lo(dir)
        cpp_training_hi(dir)=training_hi(dir)
       enddo

      else if (op_training.eq.1) then

       call gridsten_level(xsten,i,j,k,max_level,nhalf)

       do dir=1,SDIM
        if (cpp_training_lo(dir).eq.training_lo(dir)) then
         ! do nothing
        else
         print *,"mismatch cpp_training_lo and training_lo"
         stop
        endif
        if (cpp_training_hi(dir).eq.training_hi(dir)) then
         ! do nothing
        else
         print *,"mismatch cpp_training_hi and training_hi"
         stop
        endif
       enddo ! do dir=1,sdim

       print *,"training data num_samples,i,j,k,continuous_mof ", &
          num_samples,i,j,k,continuous_mof

       Call random_number(vof_training) ! 0<=vof_training<1

       Call random_number(phi_training) ! 0<=phi_training<1

       if (SDIM.eq.2) then
        phi_training=(phi_training-half) * Pi * two
        Do i_training = 1, num_samples
         theta_training(i_training)=zero
        enddo
       else if (SDIM.eq.3) then
        Call random_number(theta_training) ! 0<=theta_training<1
        phi_training=(phi_training-half) * Pi * two
        theta_training=theta_training * Pi 
       else
        print *,"sdim invalid"
        stop
       endif

       Do i_training = 1, num_samples

        if (SDIM.eq.2) then
         angle_exact_db(1)=phi_training(i_training)
        else if (SDIM.eq.3) then
         angle_exact_db(1)=phi_training(i_training)
         angle_exact_db(SDIM-1)=theta_training(i_training)
        else
         print *,"dimension bust"
         stop
        endif

        try_new_vfrac=1

        do while (try_new_vfrac.eq.1)

         refvfrac(1)=vof_training(i_training)

         if ((refvfrac(1).ge.zero).and. &
             (refvfrac(1).lt.VOFTOL)) then
          ! do nothing
         else if ((refvfrac(1).gt.one-VOFTOL).and. &
                  (refvfrac(1).le.one)) then
          ! do nothing
         else if ((refvfrac(1).ge.VOFTOL).and. &
                  (refvfrac(1).le.one-VOFTOL)) then

           ! given the slope, find the centroid.
          call angle_init_from_angle_recon_and_F( &
           bfact,dx, &
           xsten,nhalf, &
           refvfrac, & 
           continuous_mof, & !=STANDARD_MOF or CMOF_X
           cmofsten, & 
           geom_xtetlist(1,1,1,tid+1), &
           geom_xtetlist_old(1,1,1,tid+1), &
           nmax, &
           angle_init_db, & ! INTENT(out)
           refcen, &  ! INTENT(out)
           angle_exact_db, & ! INTENT(in)
           nmax, &
           SDIM)

          try_new_vfrac=0

          do dir=1,SDIM-1

           if ((angle_init_db(dir).lt.-Pi).and. &
               (angle_init_db(dir).ge.-Pi-EPS_8_4)) then
            try_new_vfrac=1
           else if ((angle_init_db(dir).gt.Pi).and. &
                    (angle_init_db(dir).le.Pi+EPS_8_4)) then
            try_new_vfrac=1
           else if ((angle_init_db(dir).ge.-Pi).and. &
                    (angle_init_db(dir).le.Pi)) then
            ! do nothing
           else
            print *,"angle_init_db invalid"
            print *,"dir=",dir
            print *,"angle_init_db(dir)=",angle_init_db(dir)
            stop
           endif

           if ((angle_exact_db(dir).lt.-Pi).and. &
               (angle_exact_db(dir).ge.-Pi-EPS_8_4)) then
            try_new_vfrac=1
           else if ((angle_exact_db(dir).gt.Pi).and. &
                    (angle_exact_db(dir).le.Pi+EPS_8_4)) then
            try_new_vfrac=1
           else if ((angle_exact_db(dir).ge.-Pi).and. &
                    (angle_exact_db(dir).le.Pi)) then
            ! do nothing
           else
            print *,"angle_exact_db invalid"
            print *,"dir=",dir
            print *,"angle_exact_db(dir)=",angle_exact_db(dir)
            stop
           endif

          enddo !dir=1,sdim-1

          if (SDIM.eq.2) then
           ! do nothing
          else if (SDIM.eq.3) then
           dir=SDIM-1

           if ((angle_init_db(dir).lt.zero).and. &
               (angle_init_db(dir).ge.-EPS_8_4)) then
            try_new_vfrac=1
           else if ((angle_init_db(dir).gt.Pi).and. &
                    (angle_init_db(dir).le.Pi+EPS_8_4)) then
            try_new_vfrac=1
           else if ((angle_init_db(dir).ge.zero).and. &
                    (angle_init_db(dir).le.Pi)) then
            ! do nothing
           else
            print *,"angle_init_db invalid"
            print *,"dir=",dir
            print *,"angle_init_db(dir)=",angle_init_db(dir)
            stop
           endif

           if ((angle_exact_db(dir).lt.zero).and. &
               (angle_exact_db(dir).ge.-EPS_8_4)) then
            try_new_vfrac=1
           else if ((angle_exact_db(dir).gt.Pi).and. &
                    (angle_exact_db(dir).le.Pi+EPS_8_4)) then
            try_new_vfrac=1
           else if ((angle_exact_db(dir).ge.zero).and. &
                    (angle_exact_db(dir).le.Pi)) then
            ! do nothing
           else
            print *,"angle_exact_db invalid"
            print *,"dir=",dir
            print *,"angle_exact_db(dir)=",angle_exact_db(dir)
            stop
           endif

          else
           print *,"dimension bust"
           stop
          endif

          if (try_new_vfrac.eq.0) then

           do dir=1,SDIM
            centroid_null(dir)=zero
           enddo
           vof_null=one

           critical_material=1
           call find_predict_slope( &
            npredict, & ! INTENT(out)
            mag_centroid, & ! INTENT(out)
            vof_null, &
            vof_null, &
            centroid_null, & ! centroid of uncaptured region
            refvfrac(1), &
             ! relative to cell centroid of the super cell; INTENT(in)
            refcen, & 
            bfact,dx,xsten,nhalf,SDIM)

           try_new_vfrac=0

           if (mag_centroid(1).gt.VOFTOL*dx(1)) then
            ! do nothing
           else if (mag_centroid(1).le.VOFTOL*dx(1)) then
            try_new_vfrac=1
           else
            print *,"mag_centroid bust: ",mag_centroid(1)
            stop
           endif

          else if (try_new_vfrac.eq.1) then
           !do nothing
          else
           print *,"try_new_vfrac invalid"
           stop
          endif

         else 
          print *,"refvfrac(1) out of range: ",refvfrac(1)
          stop
         endif

         if (try_new_vfrac.eq.1) then
          Call random_number(vof_single)
          vof_training(i_training)=vof_single
         else if (try_new_vfrac.eq.0) then
          ! do nothing
         else
          print *,"try_new_vfrac invalid"
          stop
         endif

        enddo ! do while (try_new_vfrac.eq.1)

        do dir=1,SDIM
         xc0(dir)=refcen(dir)
        enddo

        if (NTRAINING.eq.ANGLE_INIT_TRAIN2) then
         !do nothing
        else
         print *,"ntraining or angle_init_train2 invalid"
         stop
        endif

        do dir=1,SDIM
         data_training(dir,i_training) = xc0(dir)
        enddo
        data_training(VOFTRAIN,i_training) = vof_training(i_training)
        do dir=1,SDIM-1
         data_training(VOFTRAIN+dir,i_training) = angle_exact_db(dir)
         data_training(ANGLE_EXACT_TRAIN2+dir,i_training) = angle_init_db(dir)
        enddo
       End Do ! i_training = 1, num_samples

! unit number 5: standard input
! unit number 6: standard output
   
       call execute_command_line('rm exact_centroid.dat',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm exact_f.dat',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm exact_angle.dat',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm initial_angle.dat',wait=.true., &
               exitstat=sysret)

       call execute_command_line('rm exact_centroid.npy',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm exact_f.npy',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm exact_angle.npy',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm initial_angle.npy',wait=.true., &
               exitstat=sysret)

       call execute_command_line('rm nn_coef.dat',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm dt_coef.dat',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm rf_coef.dat',wait=.true., &
               exitstat=sysret)

       open(10,file='exact_centroid.dat',status='unknown')
       open(11,file='exact_f.dat',status='unknown')
       open(12,file='exact_angle.dat',status='unknown')
       open(13,file='initial_angle.dat',status='unknown')
        ! previous: F16.12
        !      now: E25.16
       Do i_training = 1, num_samples

        if (SDIM.eq.3) then
         Write(10,'(3E25.16)')data_training(1:SDIM,i_training)
        else if (SDIM.eq.2) then
         Write(10,'(2E25.16)')data_training(1:SDIM,i_training)
        else
         print *,"dimension bust"
         stop
        endif

        Write(11,'(E25.16)')data_training(VOFTRAIN,i_training)

        if (SDIM.eq.3) then
         Write(12,'(2E25.16)') &
           data_training(ANGLE_EXACT_TRAIN1:ANGLE_EXACT_TRAIN2,i_training)
         Write(13,'(2E25.16)') &
           data_training(ANGLE_INIT_TRAIN1:ANGLE_INIT_TRAIN2,i_training)
        else if (SDIM.eq.2) then
         Write(12,'(E25.16)') &
           data_training(ANGLE_EXACT_TRAIN1:ANGLE_EXACT_TRAIN2,i_training)
         Write(13,'(E25.16)') &
           data_training(ANGLE_INIT_TRAIN1:ANGLE_INIT_TRAIN2,i_training)
        else
         print *,"dimension bust"
         stop
        endif

       End Do ! Do i_training = 1, num_samples

       close(10)
       close(11)
       close(12)
       close(13)

       call execute_command_line('python convert2binary.py',wait=.true., &
               exitstat=sysret)

       call execute_command_line('python training.py',wait=.true., &
               exitstat=sysret)

       call execute_command_line('rm exact_centroid.dat',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm exact_f.dat',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm exact_angle.dat',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm initial_angle.dat',wait=.true., &
               exitstat=sysret)

       call execute_command_line('rm exact_centroid.npy',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm exact_f.npy',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm exact_angle.npy',wait=.true., &
               exitstat=sysret)
       call execute_command_line('rm initial_angle.npy',wait=.true., &
               exitstat=sysret)

        ! sanity test for sk2f.py
       if (1.eq.0) then
        if (continuous_mof.eq.STANDARD_MOF) then
         cmof_idx=STANDARD_MOF
        else if (continuous_mof.eq.CMOF_X) then
         cmof_idx=CMOF_X
        else
         print *,"continuous_mof invalid: ",continuous_mof
         stop
        endif
        call training_array(D_DECL(i,j,k),cmof_idx)% &
         NN_ZHOUTENG_LOCAL%Initialization()
        call training_array(D_DECL(i,j,k),cmof_idx)% &
         DT_ZHOUTENG_LOCAL%Initialization()
        call training_array(D_DECL(i,j,k),cmof_idx)% &
         RF_ZHOUTENG_LOCAL%Initialization()

        NN_cost=zero
        DT_cost=zero
        RF_cost=zero
        Do i_training = 1, num_samples

         do dir=1,SDIM-1
          angle_init_db(dir)=data_training(ANGLE_EXACT_TRAIN2+dir,i_training)
          angle_exact_db_data(dir)=data_training(VOFTRAIN+dir,i_training)
          angle_and_vfrac(dir)=angle_init_db(dir)
         enddo
         refvfrac(1)=data_training(VOFTRAIN,i_training)
         angle_and_vfrac(SDIM)=refvfrac(1)

         angle_exact_db= &
           training_array(D_DECL(i,j,k),cmof_idx)%DT_ZHOUTENG_LOCAL% &
             predict(angle_and_vfrac)
         do dir=1,SDIM-1
          DT_cost=DT_cost+(angle_exact_db(dir)-angle_exact_db_data(dir))**2
         enddo

         if (1.eq.0) then
          print *,"DT; i_training ",i_training
          print *,"angle_init_db ",angle_init_db
          print *,"angle_exact_db_data ",angle_exact_db_data
          print *,"angle_exact_db ",angle_exact_db
          print *,"angle_and_vfrac ",angle_and_vfrac
         endif

         angle_exact_db= &
           training_array(D_DECL(i,j,k),cmof_idx)%NN_ZHOUTENG_LOCAL% &
             predict(angle_and_vfrac)
         do dir=1,SDIM-1
          NN_cost=NN_cost+(angle_exact_db(dir)-angle_exact_db_data(dir))**2
         enddo

         if (1.eq.0) then
          print *,"NN; i_training ",i_training
          print *,"angle_init_db ",angle_init_db
          print *,"angle_exact_db_data ",angle_exact_db_data
          print *,"angle_exact_db ",angle_exact_db
          print *,"angle_and_vfrac ",angle_and_vfrac
         endif

         angle_exact_db= &
           training_array(D_DECL(i,j,k),cmof_idx)%RF_ZHOUTENG_LOCAL% &
             predict(angle_and_vfrac)
         do dir=1,SDIM-1
          RF_cost=RF_cost+(angle_exact_db(dir)-angle_exact_db_data(dir))**2
         enddo

         if (1.eq.0) then
          print *,"RF; i_training ",i_training
          print *,"angle_init_db ",angle_init_db
          print *,"angle_exact_db_data ",angle_exact_db_data
          print *,"angle_exact_db ",angle_exact_db

          stop
         endif

        enddo ! i_training = 1, num_samples

        print *,"FORTRAN: DT cost ",DT_cost
        print *,"FORTRAN: NN cost ",NN_cost
        print *,"FORTRAN: RF cost ",RF_cost

        stop
       endif

      else if (op_training.eq.2) then

       do dir=1,SDIM
        if (cpp_training_lo(dir).eq.training_lo(dir)) then
         ! do nothing
        else
         print *,"mismatch cpp_training_lo and training_lo"
         stop
        endif
        if (cpp_training_hi(dir).eq.training_hi(dir)) then
         ! do nothing
        else
         print *,"mismatch cpp_training_hi and training_hi"
         stop
        endif
       enddo ! do dir=1,sdim

       if (continuous_mof.eq.STANDARD_MOF) then
        cmof_idx=STANDARD_MOF
       else if (continuous_mof.eq.CMOF_X) then
        cmof_idx=CMOF_X
       else
        print *,"continuous_mof invalid: ",continuous_mof
        stop
       endif
       call training_array(D_DECL(i,j,k),cmof_idx)% &
         NN_ZHOUTENG_LOCAL%Initialization()
       call training_array(D_DECL(i,j,k),cmof_idx)% &
         DT_ZHOUTENG_LOCAL%Initialization()
       call training_array(D_DECL(i,j,k),cmof_idx)% &
         RF_ZHOUTENG_LOCAL%Initialization()

      else
       print *,"op_training invalid"
       stop
      endif

      return
      end subroutine fort_MOF_training

      !NavierStokes::post_restart calls NavierStokes::prepare_post_process
      !NavierStokes::post_init_state calls NavierStokes::prepare_post_process
      !NavierStokes::prepare_post_process calls NavierStokes::MOF_training
      !NavierStokes::MOF_training calls fort_MOF_DT_training
      subroutine fort_MOF_DT_training( &
        num_samples, & !mof_decision_tree_learning
        finest_level, &
        max_level, &
        bfact, &
        domlo,domhi, &
        dx) &
      bind(c,name='fort_MOF_DT_training')

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: num_samples
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: max_level
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(SDIM)

      integer :: i,j,k
      integer :: local_continuous_mof
      integer :: nmax=POLYGON_LIST_MAX
      real(amrex_real) :: vof_training(num_samples)
      real(amrex_real) :: phi_training(num_samples)
      real(amrex_real) :: theta_training(num_samples)
      real(amrex_real) :: &
        data_decisions(num_samples,MOF_TRAINING_NDIM_DECISIONS)
      real(amrex_real) :: data_classify(num_samples,MOF_TRAINING_NDIM_CLASSIFY)

      real(amrex_real) :: angle_exact_db(SDIM-1)
      real(amrex_real) :: angle_init_db(SDIM-1)
      real(amrex_real) :: angle_and_vfrac(MOF_TRAINING_NDIM_DECISIONS)
      real(amrex_real) :: refvfrac(1)
      real(amrex_real) :: vof_single 
      real(amrex_real) :: refcen(SDIM)
      integer :: try_new_vfrac

      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      integer dir
      integer i1,j1,k1
      integer i_training

      integer :: loc_lo(SDIM)
      integer :: loc_hi(SDIM)

      real(amrex_real) :: npredict(3,SDIM)
      real(amrex_real) :: vof_null
      real(amrex_real) :: centroid_null(SDIM)
      real(amrex_real) :: mag_centroid(3)
      integer :: critical_material

      real(amrex_real) :: DT_cost
      real(amrex_real) :: angle_exact_db_data(SDIM-1)

      integer cmof_idx
      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))
      integer klosten,khisten
      integer tid
   
      tid=0

      if (num_samples.eq.0) then

       decision_tree_max_level=-1

      else if (num_samples.gt.0) then

       if (max_level.eq.fort_max_level) then
        decision_tree_max_level=max_level
       else
        print *,"max_level and fort_max_level mismatch"
        print *,"max_level: ",max_level
        print *,"fort_max_level: ",fort_max_level
        print *,"decision_tree_max_level: ",decision_tree_max_level
        stop
       endif

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
       do k1=klosten,khisten
       do j1=-1,1
       do i1=-1,1
        cmofsten(D_DECL(i1,j1,k1))=1
       enddo
       enddo
       enddo

       if (bfact.lt.1) then
        print *,"bfact invalid170"
        stop
       endif

       if (num_state_base.ne.2) then
        print *,"num_state_base invalid"
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

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        ! do nothing
       else
        print *,"levelrz invalid fort_MOF_DT_training"
        stop
       endif

       do dir=1,3
        decision_tree_lo(dir)=0
        decision_tree_hi(dir)=0
       enddo 

       do dir=1,SDIM
        if (domlo(dir).eq.0) then
         ! do nothing
        else
         print *,"expecting domlo=0"
         stop
        endif
        decision_tree_hi(dir)=domlo(dir)+bfact-1
        if (decision_tree_hi(dir).le.domhi(dir)) then
         ! do nothing 
        else
         print *,"decision_tree_hi invalid"
         stop
        endif
       enddo ! dir=1..sdim

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        dir=1
        decision_tree_hi(dir)=domhi(dir)
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        dir=1
        decision_tree_hi(dir)=domhi(dir)
        dir=2
        decision_tree_hi(dir)=domhi(dir)
       else
        print *,"levelrz invalid fort_MOF_DT_training"
        stop
       endif

       do dir=1,SDIM
        loc_lo(dir)=decision_tree_lo(dir)
        loc_hi(dir)=decision_tree_hi(dir)
       enddo

       allocate(decision_tree_array( &
        D_DECL(loc_lo(1):loc_hi(1),loc_lo(2):loc_hi(2),loc_lo(SDIM):loc_hi(SDIM) ), &
        0:1))

       do k=decision_tree_lo(3),decision_tree_hi(3)
       do j=decision_tree_lo(2),decision_tree_hi(2)
       do i=decision_tree_lo(1),decision_tree_hi(1)
       do cmof_idx=0,1

        local_continuous_mof=cmof_idx

        call gridsten_level(xsten,i,j,k,max_level,nhalf)

        print *,"DT: training data num_samples,i,j,k,continuous_mof ", &
           num_samples,i,j,k,local_continuous_mof

        Call random_number(vof_training) ! 0<=vof_training<1

        Call random_number(phi_training) ! 0<=phi_training<1

        if (SDIM.eq.2) then
         phi_training=(phi_training-half) * Pi * two
         Do i_training = 1, num_samples
          theta_training(i_training)=zero
         enddo
        else if (SDIM.eq.3) then
         Call random_number(theta_training) ! 0<=theta_training<1
         phi_training=(phi_training-half) * Pi * two
         theta_training=theta_training * Pi
        else
         print *,"sdim invalid"
         stop
        endif

        Do i_training = 1, num_samples

         if (SDIM.eq.2) then
          angle_exact_db(1)=phi_training(i_training)
         else if (SDIM.eq.3) then
          angle_exact_db(1)=phi_training(i_training)
          angle_exact_db(SDIM-1)=theta_training(i_training)
         else
          print *,"dimension bust"
          stop
         endif

         try_new_vfrac=1

         do while (try_new_vfrac.eq.1)

          refvfrac(1)=vof_training(i_training)

          if ((refvfrac(1).ge.zero).and. &
              (refvfrac(1).lt.VOFTOL)) then
           ! do nothing
          else if ((refvfrac(1).gt.one-VOFTOL).and. &
                   (refvfrac(1).le.one)) then
           ! do nothing
          else if ((refvfrac(1).ge.VOFTOL).and. &
                   (refvfrac(1).le.one-VOFTOL)) then

            ! given the slope, find the centroid.
           call angle_init_from_angle_recon_and_F( &
            bfact,dx,xsten,nhalf, &
            refvfrac, & 
            local_continuous_mof, & 
            cmofsten, & 
            geom_xtetlist(1,1,1,tid+1), &
            geom_xtetlist_old(1,1,1,tid+1), &
            nmax, &
            angle_init_db, & ! INTENT(out)
            refcen, &  ! INTENT(out)
            angle_exact_db, & ! INTENT(in)
            nmax, &
            SDIM)

           try_new_vfrac=0

           do dir=1,SDIM-1

            if ((angle_init_db(dir).lt.-Pi).and. &
                (angle_init_db(dir).ge.-Pi-EPS_8_4)) then
             try_new_vfrac=1
            else if ((angle_init_db(dir).gt.Pi).and. &
                     (angle_init_db(dir).le.Pi+EPS_8_4)) then
             try_new_vfrac=1
            else if ((angle_init_db(dir).ge.-Pi).and. &
                     (angle_init_db(dir).le.Pi)) then
             ! do nothing
            else
             print *,"angle_init_db invalid"
             print *,"dir=",dir
             print *,"angle_init_db(dir)=",angle_init_db(dir)
             stop
            endif

            if ((angle_exact_db(dir).lt.-Pi).and. &
                (angle_exact_db(dir).ge.-Pi-EPS_8_4)) then
             try_new_vfrac=1
            else if ((angle_exact_db(dir).gt.Pi).and. &
                     (angle_exact_db(dir).le.Pi+EPS_8_4)) then
             try_new_vfrac=1
            else if ((angle_exact_db(dir).ge.-Pi).and. &
                     (angle_exact_db(dir).le.Pi)) then
             ! do nothing
            else
             print *,"angle_exact_db invalid"
             print *,"dir=",dir
             print *,"angle_exact_db(dir)=",angle_exact_db(dir)
             stop
            endif

           enddo !dir=1,sdim-1

           if (SDIM.eq.2) then
            ! do nothing
           else if (SDIM.eq.3) then
            dir=SDIM-1
 
            if ((angle_init_db(dir).lt.zero).and. &
                (angle_init_db(dir).ge.-EPS_8_4)) then
             try_new_vfrac=1
            else if ((angle_init_db(dir).gt.Pi).and. &
                     (angle_init_db(dir).le.Pi+EPS_8_4)) then
             try_new_vfrac=1
            else if ((angle_init_db(dir).ge.zero).and. &
                     (angle_init_db(dir).le.Pi)) then
             ! do nothing
            else
             print *,"angle_init_db invalid"
             print *,"dir=",dir
             print *,"angle_init_db(dir)=",angle_init_db(dir)
             stop
            endif

            if ((angle_exact_db(dir).lt.zero).and. &
                (angle_exact_db(dir).ge.-EPS_8_4)) then
             try_new_vfrac=1
            else if ((angle_exact_db(dir).gt.Pi).and. &
                     (angle_exact_db(dir).le.Pi+EPS_8_4)) then
             try_new_vfrac=1
            else if ((angle_exact_db(dir).ge.zero).and. &
                     (angle_exact_db(dir).le.Pi)) then
             ! do nothing
            else
             print *,"angle_exact_db invalid"
             print *,"dir=",dir
             print *,"angle_exact_db(dir)=",angle_exact_db(dir)
             stop
            endif

           else
            print *,"dimension bust"
            stop
           endif

           if (try_new_vfrac.eq.0) then

            do dir=1,SDIM
             centroid_null(dir)=zero
            enddo
            vof_null=one

            critical_material=1
            call find_predict_slope( &
             npredict, & ! INTENT(out)
             mag_centroid, & ! INTENT(out)
             vof_null, &
             vof_null, &
             centroid_null, & ! centroid of uncaptured region
             refvfrac(1), &
              ! relative to cell centroid of the super cell; INTENT(in)
             refcen, & 
             bfact,dx,xsten,nhalf,SDIM)

            try_new_vfrac=0

            if (mag_centroid(1).gt.VOFTOL*dx(1)) then
             ! do nothing
            else if (mag_centroid(1).le.VOFTOL*dx(1)) then
             try_new_vfrac=1
            else
             print *,"mag_centroid bust"
             stop
            endif

           else if (try_new_vfrac.eq.1) then
            !do nothing
           else
            print *,"try_new_vfrac invalid"
            stop
           endif

          else 
           print *,"refvfrac(1) out of range: ",refvfrac(1)
           stop
          endif

          if (try_new_vfrac.eq.1) then
           Call random_number(vof_single)
           vof_training(i_training)=vof_single
          else if (try_new_vfrac.eq.0) then
           ! do nothing
          else
           print *,"try_new_vfrac invalid"
           stop
          endif

         enddo ! do while (try_new_vfrac.eq.1)

         data_decisions(i_training,MOF_TRAINING_NDIM_DECISIONS) = &
            vof_training(i_training)
         do dir=1,SDIM-1
          data_classify(i_training,dir) = angle_exact_db(dir)
          data_decisions(i_training,dir) = angle_init_db(dir)
         enddo
        End Do ! i_training = 1, num_samples

        call initialize_decision_tree(data_decisions,data_classify, &
                num_samples,MOF_TRAINING_NDIM_DECISIONS, &
                MOF_TRAINING_NDIM_CLASSIFY, &
                decision_tree_array(D_DECL(i,j,k),cmof_idx))

         ! sanity test
        if (1.eq.1) then

         DT_cost=zero
         Do i_training = 1, num_samples

          do dir=1,SDIM-1
           angle_init_db(dir)=data_decisions(i_training,dir)
           angle_exact_db_data(dir)=data_classify(i_training,dir)
           angle_and_vfrac(dir)=angle_init_db(dir)
          enddo
          refvfrac(1)=data_decisions(i_training,MOF_TRAINING_NDIM_DECISIONS)
          angle_and_vfrac(MOF_TRAINING_NDIM_DECISIONS)=refvfrac(1)

          call decision_tree_predict(angle_and_vfrac,angle_exact_db, &
            MOF_TRAINING_NDIM_DECISIONS, &
            MOF_TRAINING_NDIM_CLASSIFY, &
            decision_tree_array(D_DECL(i,j,k),cmof_idx))

          do dir=1,SDIM-1
           DT_cost=DT_cost+(angle_exact_db(dir)-angle_exact_db_data(dir))**2
          enddo

          if (1.eq.0) then
           print *,"DT(fortran); i_training ",i_training
           print *,"angle_init_db ",angle_init_db
           print *,"angle_exact_db_data ",angle_exact_db_data
           print *,"angle_exact_db ",angle_exact_db
           print *,"angle_and_vfrac ",angle_and_vfrac
          endif

         enddo ! i_training = 1, num_samples

         print *,"FORTRAN: DT cost ",DT_cost

        endif

       enddo !local_continuous_mof=STANDARD_MOF,CMOF_X
       enddo !i
       enddo !j
       enddo !k

      else
       print *,"num_samples invalid: ",num_samples
       stop
      endif

      return
      end subroutine fort_MOF_DT_training


      end module plic_cpp_module

