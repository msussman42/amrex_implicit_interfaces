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

!update_flag=RECON_UPDATE_STATE_ERR_AND_CENTROID or
!update_flag=RECON_UPDATE_STATE_CENTROID if called from
! (a) post_init_state
! (b) post_restart
! (c) level_phase_change_convertALL
! (d) nonlinear_advection (just before split_scalar_advectionALL)
! (e) phase_change_code_segment
! (f) no mass_transfer_code_segment
! (g) nucleation_code_segment
!and update_centroid_after_recon=1

      ! masknbr:
      ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
      ! (2) =1 interior  =0 otherwise
      ! (3) =1 interior+ngrow-1  =0 otherwise
      ! (4) =1 interior+ngrow    =0 otherwise
      subroutine fort_sloperecon( &
        tid_in, &
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
        total_calls, &
        total_iterations, &
        total_errors, &
        ngrow_slope_recon, &
        ngrow_recon, &
        continuous_mof, &
        continuous_mof_radius) &
      bind(c,name='fort_sloperecon')

      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_in
      integer, INTENT(in) :: gridno
      integer, INTENT(in) :: level,finest_level,max_level
      integer, INTENT(in) :: nsteps

      integer, INTENT(in) :: continuous_mof
      integer, INTENT(in) :: continuous_mof_radius
      integer :: local_continuous_mof_radius
      integer, INTENT(in) :: update_flag
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: vofbc(SDIM,2)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: ngrow_slope_recon
      integer, INTENT(in) :: ngrow_recon
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(masknbr)
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(vof)
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(slopes)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) :: dx_extended(SDIM)
    
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
      integer :: dir
      integer :: igridlo(3),igridhi(3)

      integer im
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdata_extended(num_materials*ngeom_recon)
      real(amrex_real) mofdata_super(num_materials*ngeom_recon)
      real(amrex_real) vof_super(num_materials)
      real(amrex_real) vof_extended(num_materials)
      real(amrex_real) mofsten(num_materials*ngeom_recon)
      real(amrex_real) multi_centroidA(num_materials,SDIM)

      integer vofcomprecon
      integer vofcompraw

      real(amrex_real) err,errsave
      integer local_mask

      integer total_calls(num_materials)
      integer total_iterations(num_materials)
      real(amrex_real) total_errors(num_materials)
      integer i1,j1,k1
      real(amrex_real) LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)

      integer, PARAMETER :: nmax=POLYGON_LIST_MAX !in: fort_sloperecon

      integer continuous_mof_parm
     
      integer klosten,khisten
      integer klostenLS,khistenLS
      integer, parameter :: nhalf=3
      integer :: nhalf_extend

      integer ihalf
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_extended(-nhalf:nhalf,SDIM)
      real(amrex_real), dimension(:,:), allocatable :: xsten_temp
      real(amrex_real) xstencil_point(SDIM)

      real(amrex_real) xstenbox(-1:1,SDIM)
      integer num_fluid_materials_in_cell
      real(amrex_real) volume_super
      real(amrex_real) volsten
      real(amrex_real) volmat
      real(amrex_real) censten(SDIM)
      real(amrex_real) cen_super(SDIM)
      real(amrex_real) voflist_center(num_materials)
      real(amrex_real) voflist_test
      integer mof_verbose
      integer, parameter :: use_ls_data=1
      integer, parameter :: nhalfbox_sten=1
      real(amrex_real) dxmaxLS
      integer debugslope
      integer, parameter :: tessellate=TESSELLATE_FLUIDS
      integer, parameter :: shapeflag=0
      integer, parameter :: continuous_mof_standard=STANDARD_MOF

      real(amrex_real) vfrac_fluid_sum
      real(amrex_real) vfrac_solid_sum
      real(amrex_real) vfrac_elastic_sum
      real(amrex_real) vfrac_local(num_materials)

      real(amrex_real) :: xtet(SDIM+1,SDIM)
      real(amrex_real) :: cmof_centroid
      real(amrex_real) :: multi_area(num_materials)
      real(amrex_real) :: multi_volume(num_materials)
      real(amrex_real) :: multi_cen(SDIM,num_materials)

      integer :: local_maskcov
      integer :: fluid_obscured
      real(amrex_real) LS_clamped
      real(amrex_real) VEL_clamped(SDIM)
      real(amrex_real) temperature_clamped
      integer :: prescribed_flag
      integer :: verification_flag
      integer, parameter :: caller_id=7
      integer, parameter :: caller_id2=8

#include "mofdata.H"


      maskcov_ptr=>maskcov
      masknbr_ptr=>masknbr
      vof_ptr=>vof
      LS_ptr=>LS

      slopes_ptr=>slopes
      snew_ptr=>snew

      if ((tid_in.lt.0).or.(tid_in.ge.geom_nthreads)) then
       print *,"tid_in invalid: ",tid_in
       stop
      endif

      debugslope=0

      if (bfact.lt.1) then
       print *,"bfact invalid fort_sloperecon ",bfact
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

      if (ngrow_slope_recon.ne.1) then
       print *,"ngrow_slope_recon invalid in fort_sloperecon: ", &
          ngrow_slope_recon
       stop
      endif

      if (continuous_mof.eq.STANDARD_MOF) then
       if (ngrow_recon.ge.ngrow_slope_recon) then
        !do nothing
       else
        print *,"ngrow_recon invalid ",ngrow_recon
        stop
       endif
      else if (continuous_mof.eq.CMOF_X) then !CMOF
       if (ngrow_recon.ge.ngrow_slope_recon+continuous_mof_radius) then
        !do nothing
       else
        print *,"ngrow_recon invalid ",ngrow_recon
        stop
       endif
      else
       print *,"continuous_mof invalid (fort_sloperecon): ",continuous_mof
       stop
      endif
      if (continuous_mof_radius.ge.1) then
       !do nothing
      else
       print *,"continuous_mof_radius invalid ",continuous_mof_radius
       stop
      endif

      nhalf_extend=2*continuous_mof_radius+1
      allocate(xsten_temp(-nhalf_extend:nhalf_extend,SDIM))

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

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,masknbr_ptr,1,-1)
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,vof_ptr,ngrow_recon,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      call checkbound_array(fablo,fabhi,slopes_ptr,ngrow_slope_recon,-1)

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

       local_continuous_mof_radius=continuous_mof_radius

       call gridsten_level(xsten_temp,i,j,k,level,nhalf_extend)
       do dir=1,SDIM
        xsten_extended(1,dir)=xsten_temp(nhalf_extend,dir)
        xsten_extended(-1,dir)=xsten_temp(-nhalf_extend,dir)
        xsten_extended(0,dir)=xsten_temp(0,dir)

        if (dir.eq.1) then
         if (xsten_extended(-1,dir).lt.-dx(1)*EPS1) then
          if ((levelrz.eq.COORDSYS_RZ).or. &
              (levelrz.eq.COORDSYS_CYLINDRICAL)) then
           ihalf=-nhalf_extend
           do while ((xsten_temp(ihalf,dir).lt.-dx(1)*EPS1).and. &
                     (ihalf.lt.0))
            ihalf=ihalf+1
           enddo
           if (ihalf.ge.0) then
            print *,"ihalf became corrupt ",ihalf
            stop
           endif
           if (ihalf.eq.-1) then
            xsten_extended(1,dir)=xsten_temp(1,dir)
            xsten_extended(-1,dir)=xsten_temp(-1,dir)
            xsten_extended(0,dir)=xsten_temp(0,dir)
            local_continuous_mof_radius=1
           else if ((ihalf.lt.-1).and.(ihalf.gt.-nhalf_extend)) then
            local_continuous_mof_radius=(-ihalf-1)/2
            if (local_continuous_mof_radius*2.ne.-ihalf-1) then
             print *,"ihalf invalid ",ihalf
             stop
            endif
            xsten_extended(1,dir)= &
              xsten_temp(2*local_continuous_mof_radius+1,dir)
            xsten_extended(-1,dir)= &
              xsten_temp(-(2*local_continuous_mof_radius+1),dir)
            xsten_extended(0,dir)=xsten_temp(0,dir)
           else
            print *,"ihalf invalid ",ihalf
            stop
           endif
          else if (levelrz.eq.COORDSYS_CARTESIAN) then
           !do nothing
          else
           print *,"levelrz invalid fort_sloperecon ",levelrz
           stop
          endif
         endif !xsten_extended(-1,dir)<-dx(1)*eps1 ??
        else if ((dir.eq.2).or.(dir.eq.SDIM)) then
         !do nothing
        else
         print *,"dir invalid fort_sloperecon ",dir
         stop
        endif

        dx_extended(dir)=xsten_extended(1,dir)-xsten_extended(-1,dir)
        if (dx_extended(dir).gt.zero) then
         !do nothing
        else
         print *,"dx_extended invalid: ",dx_extended
         stop
        endif
        xsten_extended(3,dir)=xsten_extended(1,dir)+dx_extended(dir)
        xsten_extended(-3,dir)=xsten_extended(-1,dir)-dx_extended(dir)
        xsten_extended(2,dir)=xsten_extended(1,dir)+half*dx_extended(dir)
        xsten_extended(-2,dir)=xsten_extended(-1,dir)-half*dx_extended(dir)
       enddo !dir=1,sdim

       if (SDIM.eq.3) then
        klosten=-local_continuous_mof_radius
        khisten=local_continuous_mof_radius
        klostenLS=-1
        khistenLS=1
       else if (SDIM.eq.2) then
        klosten=0
        khisten=0
        klostenLS=0
        khistenLS=0
       else
        print *,"dimension bust"
        stop
       endif

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir=1,SDIM
        xstencil_point(dir)=xsten(0,dir)
       enddo

       do im=1,num_materials

        vofcomprecon=(im-1)*ngeom_recon+1
        vofcompraw=(im-1)*ngeom_raw+1
        do dir=0,SDIM
         mofdata(vofcomprecon+dir)=vof(D_DECL(i,j,k),vofcompraw+dir)
        enddo

         !vof,cenref,order,slope,intercept
        do dir=1,SDIM
         mofdata(vofcomprecon+SDIM+1+dir)=zero
        enddo
        mofdata(vofcomprecon+SDIM+1)=zero !order

        if ((mofdata(vofcomprecon).ge.-0.1d0).and. &
            (mofdata(vofcomprecon).le.1.1d0)) then
         ! do nothing
        else if ((mofdata(vofcomprecon).lt.-0.1d0).or. &
                 (mofdata(vofcomprecon).gt.1.1d0)) then
         print *,"mofdata(vofcomprecon) out of range"
         print *,"mofdata=",mofdata
         stop
        else
         print *,"mofdata(vofcomprecon) is NaN"
         print *,"mofdata=",mofdata
         stop
        endif

         !initialize the intercept to be zero
        mofdata(vofcomprecon+ngeom_recon-1)=zero

       enddo  ! im=1..num_materials

        ! sum of F_fluid=1
        ! sum of F_rigid<=1
       call make_vfrac_sum_ok_base( &
         xsten, &
         nhalf, &
         bfact,dx, &
         tessellate, & ! =TESSELLATE_FLUIDS
         mofdata, &  ! INTENT(inout)
         SDIM)

       vfrac_fluid_sum=zero
       vfrac_solid_sum=zero
       vfrac_elastic_sum=zero

       do im=1,num_materials
        vofcomprecon=(im-1)*ngeom_recon+1
        voflist_center(im)=mofdata(vofcomprecon)
        vof_super(im)=voflist_center(im)

        if (is_rigid(im).eq.1) then
         vfrac_solid_sum=vfrac_solid_sum+voflist_center(im)
        else if (is_elastic(im).eq.1) then
         vfrac_elastic_sum=vfrac_elastic_sum+voflist_center(im)
        else if ((is_rigid(im).eq.0).and. &
                 (is_elastic(im).eq.0)) then
         vfrac_fluid_sum=vfrac_fluid_sum+voflist_center(im)
        else
         print *,"is_rigid(im) invalid: ",im,is_rigid(im)
         print *,"or is_elastic(im) invalid: ",im,is_elastic(im)
         stop
        endif

       enddo ! im=1..num_materials

       if (abs(vfrac_fluid_sum-one).le.EPS1) then
        ! do nothing
       else
        print *,"vfrac_fluid_sum invalid: ",vfrac_fluid_sum
        stop
       endif

       if ((level.ge.0).and. &
           (level.le.finest_level)) then
        !do nothing
       else
        print *,"level invalid ",level
        stop
       endif

       do k1=klostenLS,khistenLS
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
          !do nothing
         else
          print *,"voflist_test bust: ",voflist_test
          stop
         endif
        enddo ! im=1..num_materials
       enddo
       enddo
       enddo  ! i1,j1,k1  = -1,1

       num_fluid_materials_in_cell=0
       do im=1,num_materials
        if ((is_rigid(im).eq.0).and.(is_elastic(im).eq.0)) then
         if (voflist_center(im).gt.VOFTOL) then
          num_fluid_materials_in_cell=num_fluid_materials_in_cell+1
         endif
        else if ((is_rigid(im).eq.1).or. &
                 (is_elastic(im).eq.1)) then
         ! do nothing
        else
         print *,"is_rigid invalid PLIC_3D.F90: ",im,is_rigid(im)
         print *,"or is_elastic invalid PLIC_3D.F90: ",im,is_elastic(im)
         stop
        endif
       enddo ! im=1..num_materials
 
       call SUB_clamped_LS(xstencil_point,time,LS_clamped, &
         VEL_clamped,temperature_clamped,prescribed_flag,dx)

       call SUB_verification_flag(verification_flag)

       fluid_obscured=0
       if ((vfrac_solid_sum.ge.half).or. &
           (vfrac_elastic_sum.ge.half)) then
        fluid_obscured=1
       else if ((vfrac_solid_sum.le.half).and. &
                (vfrac_elastic_sum.le.half)) then
        !do nothing
       else
        print *,"vfrac_solid_sum invalid:",vfrac_solid_sum
        print *,"or vfrac_elastic_sum invalid:",vfrac_elastic_sum
        stop
       endif

       if ((LS_clamped.ge.zero).and. &
           (verification_flag.eq.0)) then
        fluid_obscured=1
       else if ((LS_clamped.le.zero).or. &
                (verification_flag.eq.1)) then
        !do nothing
       else
        print *,"LS_clamped invalid:",LS_clamped
        print *,"or verification_flag invalid:",verification_flag
        stop
       endif 

       do im=1,num_materials*ngeom_recon
        mofdata_super(im)=mofdata(im)
        mofdata_extended(im)=mofdata(im)
       enddo

        !center cell
       call Box_volumeFAST(bfact,dx,xsten,nhalf, &
         volume_super,cen_super,SDIM)

       if ((level.lt.finest_level).or. &
           (fluid_obscured.eq.1)) then

        continuous_mof_parm=STANDARD_MOF

       else if ((level.ge.finest_level).and. &
                 (fluid_obscured.eq.0)) then

        if (num_fluid_materials_in_cell.eq.1) then
         continuous_mof_parm=STANDARD_MOF
        else if (num_fluid_materials_in_cell.eq.2) then
         continuous_mof_parm=continuous_mof
        else if ((num_fluid_materials_in_cell.ge.3).and. &
                 (num_fluid_materials_in_cell.le.num_materials)) then
         continuous_mof_parm=continuous_mof
        else
         print *,"num_fluid_materials_in_cell invalid: ", &
          num_fluid_materials_in_cell
         stop
        endif

       else
        print *,"level or fluid_obscured invalid"
        print *,"level=",level
        print *,"max_level=",max_level
        print *,"fluid_obscured=",fluid_obscured
        stop
       endif

        ! supercell for centroid cost function.
        ! center cell for volume constraint.
       if (continuous_mof_parm.eq.CMOF_X) then
  
        volume_super=zero ! volume of the extended region

        do dir=1,SDIM
         cen_super(dir)=zero
        enddo

        do im=1,num_materials

         vofcomprecon=(im-1)*ngeom_recon+1

         do dir=0,SDIM
          mofdata_extended(vofcomprecon+dir)=zero
         enddo

         vof_extended(im)=zero

        enddo ! im=1..num_materials

        do k1=klosten,khisten
        do j1=-local_continuous_mof_radius,local_continuous_mof_radius
        do i1=-local_continuous_mof_radius,local_continuous_mof_radius

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
           mofsten(vofcomprecon+SDIM+1+dir)=zero
          enddo

          mofsten(vofcomprecon+SDIM+1)=zero !order

          !initialize the intercept to be zero
          mofsten(vofcomprecon+ngeom_recon-1)=zero

         enddo  ! im=1..num_materials

          ! sum of F_fluid=1
          ! sum of F_rigid<=1
         call make_vfrac_sum_ok_base( &
           xstenbox, &
           nhalfbox_sten, & ! =1
           bfact,dx, &
           tessellate, & ! =TESSELLATE_FLUIDS
           mofsten, &
           SDIM)

         do im=1,num_materials
          vofcomprecon=(im-1)*ngeom_recon+1
          vfrac_local(im)=mofsten(vofcomprecon)
          volmat=volsten*vfrac_local(im)
          do dir=1,SDIM
           mofdata_extended(vofcomprecon+dir)= &
               mofdata_extended(vofcomprecon+dir)+ &
               volmat*(censten(dir)+mofsten(vofcomprecon+dir))
          enddo ! dir
          mofdata_extended(vofcomprecon)= &
            mofdata_extended(vofcomprecon)+volmat
          vof_extended(im)=vof_extended(im)+volmat
         enddo ! im=1..num_materials

         volume_super=volume_super+volsten
         do dir=1,SDIM
          cen_super(dir)=cen_super(dir)+volsten*censten(dir)
         enddo

        enddo
        enddo
        enddo ! i1,j1,k1

        if (volume_super.gt.zero) then
         ! do nothing
        else
         print *,"volume_super invalid: ",volume_super
         stop
        endif

        do dir=1,SDIM
         cen_super(dir)=cen_super(dir)/volume_super
        enddo

        do im=1,num_materials
         vofcomprecon=(im-1)*ngeom_recon+1

         if (vof_extended(im).gt.zero) then
          do dir=1,SDIM
           mofdata_extended(vofcomprecon+dir)= &
             mofdata_extended(vofcomprecon+dir)/ &
             vof_extended(im)- &
             cen_super(dir)
          enddo
          vof_extended(im)=vof_extended(im)/volume_super
          mofdata_extended(vofcomprecon)= &
             mofdata_extended(vofcomprecon)/volume_super
         else if (vof_extended(im).eq.zero) then
          do dir=0,SDIM
           mofdata_extended(vofcomprecon+dir)=zero
          enddo
          vof_extended(im)=zero
         else
          print *,"vof_extended invalid ",vof_extended
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

       if (continuous_mof_parm.eq.STANDARD_MOF) then

        call multimaterial_MOF( &
         tessellate, & ! =TESSELLATE_FLUIDS
         tid_in, &
         bfact,dx, &
         xsten, &
         nhalf, &
         mof_verbose, &
         use_ls_data, & ! use_ls_data=1
         LS_stencil, &
         geom_xtetlist(1,1,1,tid_in+1), &
         nmax, &
         nmax, &
         mofdata_super, & !intent(inout)
         vof_super, &
         multi_centroidA, & ! (num_materials,sdim) relative to supercell
         SDIM)

        ! mof_calls, mof_iterations, mof_errors are init. in 
        ! multimaterial_MOF
        do im=1,num_materials
         total_calls(im)=total_calls(im)+mof_calls(tid_in+1,im)
         total_iterations(im)= &
          total_iterations(im)+mof_iterations(tid_in+1,im)
         total_errors(im)= &
          total_errors(im)+mof_errors(tid_in+1,im)
        enddo  ! im=1..num_materials

       else if (continuous_mof_parm.eq.CMOF_X) then

        call multimaterial_MOF( &
          tessellate, & ! =TESSELLATE_FLUIDS
          tid_in, &
          bfact, &
          dx_extended, &
          xsten_extended, &
          nhalf, &
          mof_verbose, &
          use_ls_data, & ! use_ls_data=1
          LS_stencil, &
          geom_xtetlist(1,1,1,tid_in+1), &
          nmax, &
          nmax, &
          mofdata_extended, & !intent(inout)
          vof_extended, & !intent(in)
          multi_centroidA, & ! (num_materials,sdim) relative to supercell
          SDIM)

        ! mof_calls, mof_iterations, mof_errors are init. in 
        ! multimaterial_MOF
        do im=1,num_materials
          total_calls(im)=total_calls(im)+mof_calls(tid_in+1,im)
          total_iterations(im)= &
           total_iterations(im)+mof_iterations(tid_in+1,im)
          total_errors(im)= &
           total_errors(im)+mof_errors(tid_in+1,im)
        enddo  ! im=1..num_materials

        do im=1,num_materials

         vofcomprecon=(im-1)*ngeom_recon+1

         do dir=1,SDIM
          mofdata_super(vofcomprecon+SDIM+1+dir)= &
           mofdata_extended(vofcomprecon+SDIM+1+dir)
         enddo
         mofdata_super(vofcomprecon+SDIM+1)= &
           mofdata_extended(vofcomprecon+SDIM+1)

        enddo ! im=1,num_materials

        call multimaterial_MOF( &
          tessellate, & ! =TESSELLATE_FLUIDS
          tid_in, &
          bfact,dx, &
          xsten, &
          nhalf, &
          mof_verbose, &
          use_ls_data, & ! use_ls_data=1
          LS_stencil, &
          geom_xtetlist(1,1,1,tid_in+1), &
          nmax, &
          nmax, &
          mofdata_super, & !intent(inout)
          vof_super, &
          multi_centroidA, & ! (num_materials,sdim) relative to supercell
          SDIM)

        ! mof_calls, mof_iterations, mof_errors are init. in 
        ! multimaterial_MOF
        do im=1,num_materials
          total_calls(im)=total_calls(im)+mof_calls(tid_in+1,im)
          total_iterations(im)= &
           total_iterations(im)+mof_iterations(tid_in+1,im)
          total_errors(im)= &
           total_errors(im)+mof_errors(tid_in+1,im)
        enddo  ! im=1..num_materials
 
       else
        print *,"continuous_mof_parm invalid"
        print *,"continuous_mof_parm=",continuous_mof_parm
        stop
       endif

       if ((update_flag.eq.RECON_UPDATE_STATE_CENTROID).or. &
           (update_flag.eq.RECON_UPDATE_STATE_ERR_AND_CENTROID)) then

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volume_super,cen_super,SDIM)

        call multi_get_volume_grid( &
         caller_id2, &
         tid_in, &
         EPS_11_4, &
         tessellate, & ! =TESSELLATE_FLUIDS
         tessellate, & ! =TESSELLATE_FLUIDS
         bfact,dx, &
         xsten,nhalf, & ! phi = n dot (x-x0) + intercept
         mofdata_super, &
         xsten,nhalf, & ! find volumes within xsten
         xtet, &        ! not within xtet
         multi_volume, &
         multi_cen, & !(sdim,num_materials) absolute frame of ref.
         multi_area, & !(num_materials)
         geom_xtetlist_uncapt(1,1,1,tid_in+1), &
         nmax, &
         nmax, &
         SDIM, &
         shapeflag)

        do im=1,num_materials

         if (local_maskcov.eq.0) then
          ! do nothing
         else if (local_maskcov.eq.1) then
           
          vofcomprecon=(im-1)*ngeom_recon+1
          vofcompraw=(im-1)*ngeom_raw+1

          vfrac_local(im)=mofdata_super(vofcomprecon)

          if ((vfrac_local(im).ge.0.01D0).and. &
              (vfrac_local(im).le.0.99D0)) then

           do dir=1,SDIM
            cmof_centroid=multi_cen(dir,im)-cen_super(dir)
            mofdata_super(vofcomprecon+dir)=cmof_centroid
            snew(D_DECL(i,j,k),vofcompraw+dir)=cmof_centroid
           enddo

          else if (abs(vfrac_local(im)).le.0.01D0) then
           ! do nothing
          else if (abs(one-vfrac_local(im)).le.0.01D0) then
           ! do nothing
          else
           print *,"vfrac_local invalid ",vfrac_local
           stop
          endif

         else
          print *,"local_maskcov invalid: ",local_maskcov
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
          print *,"err bust: ",errsave,err
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

      deallocate(xsten_temp)

      return
      end subroutine fort_sloperecon

      end module plic_cpp_module

