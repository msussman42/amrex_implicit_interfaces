!get rid of autoindent   :setl noai nocin nosi inde=
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "INTEGRATED_QUANTITY.H"
#include "EXTRAP_COMP.H"
#include "GODUNOV_F.H"


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif


      module godunov_module
      use amrex_fort_module, only : amrex_real
      use probf90_module

      contains

      subroutine departure_node_split( &
        xstenMAC,nhalf,dx,bfact, &
        delta,passive_veltime, &
        normdir,dt,map_forward)
      use mass_transfer_module
      IMPLICIT NONE

      integer, INTENT(in) :: normdir,nhalf,bfact
      integer :: dir
      real(amrex_real), INTENT(in) :: xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(inout) :: delta
      real(amrex_real), INTENT(in) :: dx(SDIM) 
      real(amrex_real) vel0_test(SDIM)
      real(amrex_real) x0(SDIM)
      integer, INTENT(in) :: map_forward
      real(amrex_real), INTENT(in) :: dt,passive_veltime
      real(amrex_real) :: RR
      real(amrex_real) LS_clamped
      real(amrex_real) temperature_clamped
      integer prescribed_flag

      if (nhalf.lt.1) then
       print *,"nhalf invalid departure node split"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid41"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif
      if ((map_forward.ne.0).and.(map_forward.ne.1)) then
       print *,"map_forward invalid"
       stop
      endif

      RR=one
      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
       if (normdir.eq.1) then
        RR=xstenMAC(0,1)
       endif
      else
       print *,"levelrz invalid departure node split"
       stop
      endif

      do dir=1,SDIM
       x0(dir)=xstenMAC(0,dir)
      enddo

      call SUB_clamped_LS(x0,passive_veltime,LS_clamped, &
         vel0_test,temperature_clamped,prescribed_flag,dx)
       
      if (LS_clamped.ge.zero) then 
       delta=dt*vel0_test(normdir+1)/RR
      else if (LS_clamped.lt.zero) then
       ! do nothing
      else
       print *,"LS_clamped invalid"
       stop
      endif

      if ((levelrz.eq.COORDSYS_RZ).or. &
          (levelrz.eq.COORDSYS_CYLINDRICAL)) then
       if (normdir.eq.0) then
        if (x0(normdir+1).le.EPS2*dx(normdir+1)) then
         delta=zero
        else
         call adjust_du(delta,normdir,x0(normdir+1),map_forward)
        endif
       else if ((normdir.eq.1).or.(normdir.eq.SDIM-1)) then
        if (x0(1).le.zero) then
         delta=zero
        endif
       else
        print *,"normdir invalid"
        stop
       endif
      endif  ! levelrz=1 or 3

      return
      end subroutine departure_node_split

      subroutine derive_density( &
       voldepart,voltarget,voltotal_depart, & !intent(in)
       constant_density_all_time, & !intent(in)
       massdepart, & !intent(in)
       im, & !intent(in)
       density)  !intent(out)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: im
      integer, INTENT(in) :: constant_density_all_time(num_materials)
      real(amrex_real), INTENT(in) :: voldepart,voltarget,voltotal_depart
      real(amrex_real), INTENT(in) :: massdepart
      real(amrex_real), INTENT(out) :: density
     
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid22"
       stop
      endif 

       ! sanity check
      if (is_rigid(im).eq.1) then
       if (fort_material_type(im).ne.999) then
        print *,"fort_material_type(im) invalid"
        stop
       endif
      else if (is_rigid(im).eq.0) then
       ! do nothing
      else
       print *,"is_rigid(im) invalid"
       stop
      endif

      if ((is_rigid(im).eq.1).or. &
          (fort_material_type(im).eq.999)) then
       if (constant_density_all_time(im).eq.1) then
        density=fort_denconst(im)
       else
        print *,"constant_density_all_time invalid"
        stop
       endif
      else if ((voldepart.le.VOFTOL*voltotal_depart).or. &
               (voltarget.le.VOFTOL*voltotal_depart)) then
       density=fort_denconst(im)
      else if ((fort_material_type(im).ge.0).and. &
               (fort_material_type(im).le.MAX_NUM_EOS)) then

       if (fort_material_type(im).eq.0) then

         if (constant_density_all_time(im).eq.1) then
          density=fort_denconst(im)
         else if (constant_density_all_time(im).eq.0) then 
          density=massdepart/voldepart
         else
          print *,"constant_density_all_time invalid"
          stop
         endif

       else if (fort_material_type(im).gt.0) then
        if (constant_density_all_time(im).eq.0) then 
         density=massdepart/voltarget
        else
         print *,"constant_density_all_time invalid"
         stop
        endif
       else 
        print *,"fort_material_type is corrupt"
        stop
       endif

      else
       print *,"fort_material_type(im) invalid"
       stop
      endif

      return
      end subroutine derive_density



        ! i,j,k is upper adjoining cell
        ! i-ii,j-jj,k-kk is lower adjoining cell
        ! is,ie,js,je,ks,ke are bounds for face stencil in the
        ! TANGENTIAL (dirtan) direction.
        ! levelpc(im).
      subroutine slopecrossterm( &
        massfrac, &
        total_mass, &
        levelpc, &
        faceLS, & ! =0 if imL<>imR or coarse/fine bc or ext. bc
        mdata, &
        tdata, &
        ii,jj,kk, &
        i,j,k,dir, &
        dirtan, &
        tcomp, &
        is,ie,js,je,ks,ke, &
        slopeterm)
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE


       ! pointers are always INTENT(in) but the data itself inherits
       ! its INTENT property from the target.
      real(amrex_real), INTENT(in), pointer :: faceLS(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: mdata(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: tdata(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: levelpc(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: massfrac(num_materials)
      real(amrex_real), INTENT(in) :: total_mass

      integer, INTENT(in) :: ii,jj,kk
      integer, INTENT(in) :: i,j,k,dir,dirtan
      integer, INTENT(in) :: tcomp
      integer, INTENT(in) :: is,ie,js,je,ks,ke
      real(amrex_real), INTENT(inout) :: slopeterm

      integer im1,jm1,km1,i1,j1,k1
      integer im,im_primary,im_face,imL,imR
      real(amrex_real) weight_total,sumslope
      real(amrex_real) testslope
      integer try_stencil
      real(amrex_real) LSLEFT(num_materials)
      real(amrex_real) LSRIGHT(num_materials)

      if (dir.eq.dirtan) then
       print *,"dir or dirtan invalid"
       stop
      endif
      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid slopecrossterm"
       stop
      endif
      if ((dirtan.lt.1).or.(dirtan.gt.SDIM)) then
       print *,"dirtan invalid"
       stop
      endif

      if (SDIM.eq.3) then
       if ((i.lt.is).or.(i.gt.ie).or. &
           (j.lt.js).or.(j.gt.je).or. &
           (k.lt.ks).or.(k.gt.ke)) then
        print *,"i,j,k invalid"
        stop
       endif
      else if (SDIM.eq.2) then
       if ((i.lt.is).or.(i.gt.ie).or. &
           (j.lt.js).or.(j.gt.je)) then
        print *,"i,j invalid"
        stop
       endif
      else
       print *,"dimension bust"
       stop
      endif

      im1=i-ii
      jm1=j-jj
      km1=k-kk
   
      if (total_mass.gt.zero) then ! not in RZ or r>0.

        ! massfrac derived from tessellating interfaces.
       im_primary=0
       do im=1,num_materials
        if (im_primary.eq.0) then
         im_primary=im
        else if ((im_primary.ge.1).and.(im_primary.le.num_materials)) then
         if (massfrac(im).gt.massfrac(im_primary)) then
          im_primary=im
         endif
        else
         print *,"im_primary invalid"
         stop
        endif
       enddo ! im
    
       if ((im_primary.ge.1).and.(im_primary.le.num_materials).and. &
           (massfrac(im_primary).gt.zero)) then

        do im=1,num_materials
         LSLEFT(im)=levelpc(D_DECL(im1,jm1,km1),im)
         LSRIGHT(im)=levelpc(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(LSLEFT,imL)
        call get_primary_material(LSRIGHT,imR)
 
        if ((imL.ne.im_primary).and. &
            (imR.ne.im_primary)) then 
         ! do not increment slopeterm
        else if ((imL.eq.im_primary).or. &
                 (imR.eq.im_primary)) then

         if (mdata(D_DECL(i,j,k),dir).eq.zero) then
          ! do not increment slopeterm if both adjoining cells are 
          ! solid or the cell pair is outside the grid.
         else if (mdata(D_DECL(i,j,k),dir).eq.one) then 
          ! one of adjoining cells is fluid
          weight_total=zero
          sumslope=zero

          do k1=ks,ke 
          do j1=js,je 
          do i1=is,ie 

            !im_face=0 if, wrt i1,j1,k1, imL_1<>imR_1 or 
            !coarse/fine or exterior BC or is_clamped_face>=1
           im_face=NINT(faceLS(D_DECL(i1,j1,k1),dirtan))
           if ((im_face.ge.0).and.(im_face.le.num_materials)) then
            ! do nothing
           else
            print *,"im_face invalid"
            stop
           endif

            ! the solid is represented as a rasterized boundary.
            ! if n=(1 0 0) then 2 mu D dot n=
            ! ( 2 u_x   u_y + v_x  u_z + w_x
            !   v_x+u_y   2 v_y    v_z + w_y
            !   u_z + w_x  v_z+w_y   2 w_z   ) dot (1 0 0) =
            ! (2 u_x  v_x+u_y  u_z+w_x)
           if (im_face.eq.0) then
            try_stencil=0
           else if (is_rigid(im_face).eq.1) then 
            try_stencil=0  
           else if (is_rigid(im_face).eq.0) then ! im_face=fluid
            if ((im_face.ge.1).and.(im_face.le.num_materials)) then
             if (im_face.eq.im_primary) then
              try_stencil=1
             else if (im_face.ne.im_primary) then
              if (is_rigid(im_primary).eq.1) then ! im_primary=solid
               try_stencil=1
              else if (is_rigid(im_primary).eq.0) then
               try_stencil=0
              else
               print *,"is_rigid(im_primary) invalid"
               stop
              endif
             else
              print *,"im_face invalid"
              stop
             endif
            else
             print *,"im_face invalid"
             stop
            endif
           else
            print *,"is_rigid(im_face) invalid"
            stop
           endif
    
           if (try_stencil.eq.0) then
            ! do nothing - solid face or face not in stencil.
           else if (try_stencil.eq.1) then
    
            weight_total=weight_total+one
            testslope=tdata(D_DECL(i1,j1,k1),tcomp)
            sumslope=sumslope+testslope

           else
            print *,"try_stencil invalid"
            stop
           endif
          enddo
          enddo
          enddo ! i1,j1,k1

          if (weight_total.gt.zero) then

           sumslope=sumslope/weight_total
           slopeterm=slopeterm+sumslope

          else if (weight_total.eq.zero) then
           ! do nothing
          else
           print *,"weight_total invalid"
           stop
          endif 
         else
          print *,"mdata invalid"
          stop
         endif

        else
         print *,"imL or imR bust"
         stop 
        endif 
 
       else
        print *,"im_primary invalid"
        stop
       endif

      else if (total_mass.eq.zero) then
       ! do nothing
      else
       print *,"total_mass invalid"
       stop
      endif

      return
      end subroutine slopecrossterm



      !xsten_accept=stencil for ipart,jpart,kpart,iside_part
      !xsten_recon=stencil for idonate,jdonate,kdonate
      !xsten_donate=stencil for idonate,jdonate,kdonate,isidedonate
      !xsten_depart default value: xsten_accept
      !xsten_target default value: xsten_accept
      subroutine derive_mappings( &
       xsten_accept, &
       xsten_donate, &
       xsten_target, &
       xsten_depart, &
       usten_accept, &
       usten_donate, &
       xdepartsize, &
       xtargetsize, &
       xloint, &
       xhiint, &
       volint, &
       coeff, &
       bfact,dx,map_forward,normdir)
      IMPLICIT NONE

      integer map_forward,bfact,normdir
      real(amrex_real) dx(SDIM)
      real(amrex_real) xsten_accept(-1:1,SDIM)
      real(amrex_real) xsten_donate(-1:1,SDIM)
      real(amrex_real) xsten_target(-1:1,SDIM)
      real(amrex_real) xsten_depart(-1:1,SDIM)
      real(amrex_real) usten_accept(-1:1)
      real(amrex_real) usten_donate(-1:1)
      real(amrex_real) xdepartsize,xtargetsize,xloint,xhiint
      real(amrex_real) volint
      real(amrex_real) coeff(2)
      real(amrex_real) coeffINV(2)
      integer ihalf

      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid42"
       stop
      endif
      do ihalf=-1,1
       xsten_target(ihalf,normdir+1)=xsten_accept(ihalf,normdir+1)
       xsten_depart(ihalf,normdir+1)=xsten_donate(ihalf,normdir+1)
      enddo

      if (map_forward.eq.0) then ! backwards tracing

       do ihalf=-1,1
        xsten_depart(ihalf,normdir+1)=xsten_accept(ihalf,normdir+1)- &
         usten_accept(ihalf)
       enddo

      else if (map_forward.eq.1) then ! forwards tracing

       do ihalf=-1,1
        xsten_target(ihalf,normdir+1)=xsten_donate(ihalf,normdir+1)+ &
         usten_donate(ihalf)
       enddo

      else
       print *,"map_forward invalid"
       stop
      endif

      xdepartsize=xsten_depart(1,normdir+1)-xsten_depart(-1,normdir+1)
      xtargetsize=xsten_target(1,normdir+1)-xsten_target(-1,normdir+1)

      if ((xdepartsize.lt.VOFTOL*dx(normdir+1)).or. &
          (xdepartsize.gt.two*bfact*dx(normdir+1)).or. &
          (xtargetsize.lt.VOFTOL*dx(normdir+1)).or. &
          (xtargetsize.gt.two*bfact*dx(normdir+1)).or. &
          (xdepartsize.le.VOFTOL*xtargetsize).or. &
          (xtargetsize.le.VOFTOL*xdepartsize)) then
       print *,"WARNING xtarget or xdepart invalid size"
       print *,"xdepartsize ",xdepartsize
       print *,"xtargetsize ",xtargetsize
       print *,"dx= ",dx(normdir+1)
       print *,"map_forward= ",map_forward
       print *,"try lowering cfl<1/2 "
      endif

      if (map_forward.eq.0) then 
        ! xdonate is box that is intersected with the departure region.
       xloint=max(xsten_donate(-1,normdir+1),xsten_depart(-1,normdir+1))
       xhiint=min(xsten_donate(1,normdir+1),xsten_depart(1,normdir+1))
      else if (map_forward.eq.1) then
       xloint=max(xsten_accept(-1,normdir+1),xsten_target(-1,normdir+1))
       xhiint=min(xsten_accept(1,normdir+1),xsten_target(1,normdir+1))
      else
       print *,"map_forward invalid"
       stop
      endif
      volint=xhiint-xloint

      if (volint.gt.VOFTOL*dx(normdir+1)) then
       coeff(1)=xtargetsize/xdepartsize
       coeff(2)=xsten_target(-1,normdir+1)- &
         coeff(1)*xsten_depart(-1,normdir+1)
       coeffINV(1)=xdepartsize/xtargetsize
       coeffINV(2)=xsten_depart(-1,normdir+1)- &
         coeffINV(1)*xsten_target(-1,normdir+1)

       if (map_forward.eq.0) then ! xloint,xhiint in depart region
        xsten_target(-1,normdir+1)=coeff(1)*xloint+coeff(2)
        xsten_target(1,normdir+1)=coeff(1)*xhiint+coeff(2)
        xsten_depart(-1,normdir+1)=xloint
        xsten_depart(1,normdir+1)=xhiint
       else if (map_forward.eq.1) then ! xloint,xhiint in target region
        xsten_depart(-1,normdir+1)=coeffINV(1)*xloint+coeffINV(2)
        xsten_depart(1,normdir+1)=coeffINV(1)*xhiint+coeffINV(2)
        xsten_target(-1,normdir+1)=xloint
        xsten_target(1,normdir+1)=xhiint
       else
        print *,"map_forward invalid"
        stop
       endif

      else if ((volint.ge.zero).and. &
               (volint.le.VOFTOL*dx(normdir+1))) then
       volint=zero
      else if (volint.lt.zero) then
       volint=zero
      else
       print *,"volint invalid"
       stop
      endif

      return
      end subroutine derive_mappings

      subroutine get_default_scalar_diffusion(project_option, &
                     thermal_k, &
                     LS1,im_source,im_dest, &
                     den, &
                     heatcoeff)
      use global_utility_module
      IMPLICIT NONE
      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: im_source
      integer, INTENT(in) :: im_dest
      real(amrex_real), INTENT(in) :: thermal_k(num_materials)
      real(amrex_real), INTENT(in) :: den
      real(amrex_real), INTENT(in) :: LS1
      real(amrex_real), INTENT(out) :: heatcoeff
      integer :: ispec

      if (den.gt.zero) then
       ! do nothing
      else
       print *,"den invalid"
       stop
      endif

      if (project_option.eq.SOLVETYPE_HEAT) then ! thermal diffusion
       if (LS1.ge.zero) then ! center cell owned by im_source
        heatcoeff=thermal_k(im_source)
       else  ! center cell owned by im_dest
        heatcoeff=thermal_k(im_dest)
       endif
      else if ((project_option.ge.SOLVETYPE_SPEC).and. & ! species diffusion
               (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
       ispec=project_option-SOLVETYPE_SPEC
       if (LS1.ge.zero) then ! center cell owned by im_source
        heatcoeff=fort_speciesviscconst(ispec*num_materials+im_source)*den
       else ! center cell owned by im_dest
        heatcoeff=fort_speciesviscconst(ispec*num_materials+im_dest)*den
       endif
      else
       print *,"project_option invalid: get_default_scalar_diffusion"
       stop
      endif
      if (heatcoeff.ge.zero) then
       ! do nothing
      else
       print *,"heatcoeff invalid"
       stop
      endif

      end subroutine get_default_scalar_diffusion


! Prior to calling this routine:
!  a) init_gradu_tensor(...,LOCAL_CELLTENSOR_MF,LOCAL_FACETENSOR_MF)
!     i)  doit_gradu_tensor  spectral_loop==0  itensor_iter==0
!     ii) doit_gradu_tensor  spectral_loop==1  itensor_iter==0
!     ii) doit_gradu_tensor  spectral_loop==0  itensor_iter==1
!     in doit_gradu_tensor:
!         fort_face_gradients, tileloop==0 (low order), tileloop==1 (SEM)
!  b) spectral_loop=0,1
!     dir=1..sdim
!     tileloop=0...3
!       fort_crossterm
!
! in fort_crossterm:
!  if tileloop==0  spectral_loop==0 
!   low order fluxes (slopecrossterm)
!  if tileloop==0  spectral_loop==1
!   do nothing 
!  if tileloop==1  spectral_loop==0
!   high order fluxes
!  if tileloop==1  spectral_loop==1
!   do nothing 
!  if tileloop==2  spectral_loop==0
!   xflux=visc_coef*xface(facevisc_index+1)*(xflux+divterm)
!  if tileloop==2  spectral_loop==1
!   do nothing 
!  if tileloop==3  spectral_loop==0
!   semflux=xflux
!  if tileloop==3  spectral_loop==1
!   xflux=(semflux_in + semflux_out)/2
!
! radial velocity is negated if r<0
! -dt * visc_coef * viscface * (grad U + grad U^T)
! fluxes found on "dir" faces
      subroutine fort_crossterm( &
       nsolve, &
       tileloop, &
       dir, &  ! dir=1..sdim
       operation_flag, & ! OP_UGRAD_COUPLING_MAC
       enable_spectral, &
       spectral_loop, &
       ncfluxreg, &
       semflux,semfluxlo,semfluxhi, & 
       mask,masklo,maskhi, &  ! 1=fine/fine 0=coarse/fine
       maskcoef,maskcoeflo,maskcoefhi, & ! 1=not cov by level+1 or outside.
       faceLS,faceLSlo,faceLShi, & 
       mdata,mdatalo,mdatahi, & 
       tdata,tdatalo,tdatahi, & 
       c_tdata,c_tdatalo,c_tdatahi, & 
       maskSEM,maskSEMlo,maskSEMhi, &
       xlo,dx, &
       dt, &
       cur_time, &
       vel,vello,velhi, &
       levelpc,levelpclo,levelpchi, &
       xflux,xfluxlo,xfluxhi, &
       xface,xfacelo,xfacehi, & !FACE_VAR_MF
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       rzflag, &
       velbc, &
       visc_coef, &
       nden, &
       uncoupled_viscosity, &
       homflag) &
      bind(c,name='fort_crossterm')

      use probcommon_module
      use global_utility_module
      use MOF_routines_module
 
      IMPLICIT NONE

      integer, INTENT(in) :: nsolve
      integer, INTENT(in) :: tileloop
      integer, INTENT(in) :: spectral_loop
      integer, INTENT(in) :: ncfluxreg
      integer, INTENT(in) :: operation_flag
      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: level
      integer, INTENT(in) :: homflag
      integer :: nc
      integer, INTENT(in) :: uncoupled_viscosity
      integer, INTENT(in) :: nden
      integer, INTENT(in) :: velbc(SDIM,2,SDIM) 
      integer, INTENT(in) :: rzflag 
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer :: indexlo(SDIM),indexhi(SDIM)
      integer :: sideidx(SDIM)
      integer :: indexmid(SDIM)
      integer :: index_flux(SDIM)
      integer :: index_edge(SDIM)
      integer :: index_opp(SDIM)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: semfluxlo(SDIM),semfluxhi(SDIM)
      integer, INTENT(in) :: masklo(SDIM),maskhi(SDIM)
      integer, INTENT(in) :: maskcoeflo(SDIM),maskcoefhi(SDIM)
      integer, INTENT(in) :: maskSEMlo(SDIM),maskSEMhi(SDIM)
      integer, INTENT(in) :: faceLSlo(SDIM),faceLShi(SDIM)
      integer, INTENT(in) :: mdatalo(SDIM),mdatahi(SDIM)
      integer, INTENT(in) :: tdatalo(SDIM),tdatahi(SDIM)
      integer, INTENT(in) :: c_tdatalo(SDIM),c_tdatahi(SDIM)
      integer, INTENT(in) :: vello(SDIM),velhi(SDIM)
      integer, INTENT(in) :: levelpclo(SDIM),levelpchi(SDIM)
      integer, INTENT(in) :: xfluxlo(SDIM),xfluxhi(SDIM)
      integer, INTENT(in) :: xfacelo(SDIM),xfacehi(SDIM)
  
      real(amrex_real), INTENT(in) :: dt 
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM) 
      real(amrex_real), INTENT(in), target :: mask( &
       D_DECL(masklo(1):maskhi(1),masklo(2):maskhi(2),masklo(3):maskhi(3)))
      real(amrex_real), INTENT(in), target :: maskcoef( &
       D_DECL(maskcoeflo(1):maskcoefhi(1),maskcoeflo(2):maskcoefhi(2),maskcoeflo(3):maskcoefhi(3)))
      real(amrex_real), INTENT(inout), target :: semflux( &
       D_DECL(semfluxlo(1):semfluxhi(1),semfluxlo(2):semfluxhi(2),semfluxlo(3):semfluxhi(3)),ncfluxreg)
      real(amrex_real), pointer :: semflux_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: faceLS( &
       D_DECL(faceLSlo(1):faceLShi(1),faceLSlo(2):faceLShi(2),faceLSlo(3):faceLShi(3)),SDIM)

      real(amrex_real), pointer :: faceLS_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: mdata( &
       D_DECL(mdatalo(1):mdatahi(1),mdatalo(2):mdatahi(2),mdatalo(3):mdatahi(3)),SDIM)

      real(amrex_real), pointer :: mdata_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: tdata( &
       D_DECL(tdatalo(1):tdatahi(1),tdatalo(2):tdatahi(2),tdatalo(3):tdatahi(3)),AMREX_SPACEDIM_SQR)

      real(amrex_real), pointer :: tdata_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: c_tdata( &
       D_DECL(c_tdatalo(1):c_tdatahi(1),c_tdatalo(2):c_tdatahi(2),c_tdatalo(3):c_tdatahi(3)),AMREX_SPACEDIM_SQR)

      real(amrex_real), INTENT(in), target :: maskSEM( &
       D_DECL(maskSEMlo(1):maskSEMhi(1),maskSEMlo(2):maskSEMhi(2),maskSEMlo(3):maskSEMhi(3)))

      real(amrex_real), INTENT(in), target :: vel( &
       D_DECL(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3)),STATE_NCOMP_VEL)

      real(amrex_real), INTENT(in), target :: levelpc( &
       D_DECL(levelpclo(1):levelpchi(1),levelpclo(2):levelpchi(2),levelpclo(3):levelpchi(3)),num_materials)

      real(amrex_real), pointer :: levelpc_ptr(D_DECL(:,:,:),:)

       !u
      real(amrex_real), INTENT(out), target :: xflux( &
       D_DECL(xfluxlo(1):xfluxhi(1),xfluxlo(2):xfluxhi(2),xfluxlo(3):xfluxhi(3)),nsolve)

      real(amrex_real), pointer :: xflux_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: xface( &
       D_DECL(xfacelo(1):xfacehi(1),xfacelo(2):xfacehi(2),xfacelo(3):xfacehi(3)),FACECOMP_NCOMP)

      real(amrex_real), INTENT(in) :: visc_coef

      integer, INTENT(in) :: dir
 
      integer ilo,ihi 
      integer jlo,jhi 
      integer klo,khi 

      integer i,j,k
      integer dir2
      integer ic,jc,kc
      integer dirtan(2)
      integer coupling(2)
      integer ii,jj,kk,im1,jm1,km1
      real(amrex_real) gradterm,alpha
      integer side
      integer nbase

      integer im

      real(amrex_real) LSleft(num_materials)
      real(amrex_real) LSright(num_materials)
      real(amrex_real) local_LS

      real(amrex_real) divterm
      integer compressible_face
      real(amrex_real) uxterm,vyterm,wzterm
      real(amrex_real) visc_constant
      real(amrex_real) diff_flux(SDIM)
      integer imL,imR
      real(amrex_real) total_mass,DMface
      real(amrex_real) massfrac(num_materials)
      real(amrex_real) massF(2*num_materials)
      integer, PARAMETER :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) xclamped_minus_sten(-nhalf:nhalf,SDIM)
      real(amrex_real) xclamped_plus_sten(-nhalf:nhalf,SDIM)
      real(amrex_real) xclamped_minus(SDIM)
      real(amrex_real) xclamped_plus(SDIM)

      real(amrex_real) LS_clamped_minus
      real(amrex_real) LS_clamped_plus
      real(amrex_real) vel_clamped_minus(SDIM)
      real(amrex_real) vel_clamped_plus(SDIM)
      real(amrex_real) vel_clamped_face(SDIM)
      real(amrex_real) temperature_clamped_minus
      real(amrex_real) temperature_clamped_plus
      integer prescribed_flag
      integer is_clamped_face
      integer nbr_outside_domain_flag(2)
      integer nbr_covered_flag ! 0=covered 1=not covered
      integer isten

      integer local_bctype(2)
      integer local_maskSEM
      real(amrex_real) x_sep(2)
      real(amrex_real) local_bcval(2)
      real(amrex_real) local_interp(0:bfact)
      real(amrex_real) local_vel(0:bfact)
      real(amrex_real) RRface(0:bfact)
      real(amrex_real) lineflux(0:bfact,SDIM)
      real(amrex_real) local_data(1:bfact)
      real(amrex_real) local_data_side(2)
      real(amrex_real) local_grad(0:bfact)

      integer maskcov
      integer mask_out
      integer shared_face ! in: fort_crossterm
      integer test_maskSEM
      integer stripstat
      integer elemlo(3),elemhi(3)
      integer ielem,jelem,kelem
      real(amrex_real) avgflux(SDIM)
      integer i_in,j_in,k_in
      integer i_out,j_out,k_out
      integer iflux,jflux,kflux
      integer velcomp
      integer tcomp
      real(amrex_real) xflux_temp
      integer uncoupled_viscosity_override
      integer side_cell,side_face
      integer velcomp_alt
      integer inorm
      integer inorm_elem
      integer local_bc

      real(amrex_real) dxmaxLS

      real(amrex_real) local_flux_val
      real(amrex_real) local_flux_val_in
      real(amrex_real) local_flux_val_out
      integer project_option

      semflux_ptr=>semflux
      xflux_ptr=>xflux
      levelpc_ptr=>levelpc
      faceLS_ptr=>faceLS
      mdata_ptr=>mdata
      tdata_ptr=>tdata

      project_option=SOLVETYPE_VISC

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid79"
       stop
      endif

      if (nsolve.eq.SDIM) then
       ! do nothing
      else
       print *,"nsolve invalid in fort_crossterm"
       stop
      endif

      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.1)) then
       print *,"enable_spectral invalid crossterm"
       stop
      endif
      if (operation_flag.ne.OP_UGRAD_COUPLING_MAC) then
       print *,"operation_flag invalid5"
       stop
      endif
      if (ncfluxreg.ne.AMREX_SPACEDIM_SQR) then
       print *,"ncfluxreg invalid18 ",ncfluxreg
       stop
      endif
      if ((spectral_loop.ne.0).and. &
          (spectral_loop.ne.1)) then
       print *,"spectral_loop invalid"
       stop
      endif

      if ((homflag.eq.0).or.(homflag.eq.1)) then
       ! do nothing
      else
       print *,"homflag invalid 2"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt must be positive in crossterm: ",dt
       stop
      endif
      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time must be nonneg in crossterm: ",cur_time
       stop
      endif
      if ((uncoupled_viscosity.ne.0).and. &
          (uncoupled_viscosity.ne.1)) then
       print *,"uncoupled_viscosity invalid"
       stop
      endif

      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid: ",visc_coef
       stop
      endif

      do im=1,num_materials
       if (fort_denconst(im).gt.zero) then
        !do nothing
       else
        print *,"denconst invalid: ",im,fort_denconst(im)
        stop
       endif
      enddo
      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      if (rzflag.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (rzflag.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"rzflag invalid"
       stop
      endif

      call checkbound_array(fablo,fabhi,faceLS_ptr,1,-1)
      call checkbound_array(fablo,fabhi,mdata_ptr,1,-1)
      call checkbound_array(fablo,fabhi,tdata_ptr,1,-1)
      call checkbound_array(fablo,fabhi,c_tdata,1,-1)
      call checkbound_array(fablo,fabhi,vel,1,-1)
      call checkbound_array(fablo,fabhi,levelpc_ptr,2,-1)
      call checkbound_array1(fablo,fabhi,maskSEM,1,-1)
      call checkbound_array(fablo,fabhi,semflux_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,mask,1,-1)
      call checkbound_array1(fablo,fabhi,maskcoef,1,-1)
      call checkbound_array(fablo,fabhi,xflux_ptr,0,dir-1)
      call checkbound_array(fablo,fabhi,xface,0,dir-1)

      ! mdata(i,j,k,dir)=1 if at least one adjoining cell is a fluid cell.

      ii=0
      jj=0
      kk=0

      if (dir.eq.1) then
       ii=1
       nbase=TENSOR_TRANSPOSE_UX-1
      else if (dir.eq.2) then
       jj=1
       nbase=TENSOR_TRANSPOSE_UY-1
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       kk=1
       nbase=TENSOR_TRANSPOSE_UZ-1
      else
       print *,"dir invalid crossterm"
       stop
      endif

      if (dir.eq.1) then  ! fluxes on x-face
       coupling(1)=TENSOR_TRANSPOSE_UY  !VX+UY term
       coupling(2)=TENSOR_TRANSPOSE_UZ  !WX+UZ term
       dirtan(1)=2
       dirtan(2)=SDIM
      else if (dir.eq.2) then  ! fluxes on y-face
       coupling(1)=TENSOR_TRANSPOSE_VX !UY+VX term
       coupling(2)=TENSOR_TRANSPOSE_VZ !WY+VZ term
       dirtan(1)=1
       dirtan(2)=SDIM
      else if ((dir.eq.3).and.(SDIM.eq.3)) then ! fluxes on z-face
       coupling(1)=TENSOR_TRANSPOSE_WX !UZ+WX term
       coupling(2)=TENSOR_TRANSPOSE_WY !VZ+WY term
       dirtan(1)=1
       dirtan(2)=2
      else
       print *,"dir invalid crossterm 2"
       stop
      endif

      call get_dxmaxLS(dx,bfact,dxmaxLS)
      if (dxmaxLS.gt.zero) then
       !do nothing
      else
       print *,"dxmaxLS invalid: ",dxmaxLS
       stop
      endif

      if (tileloop.eq.0) then ! low order (grad U + grad U^T)

       if (spectral_loop.eq.0) then

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
          growlo,growhi,0,dir-1)

        do k=growlo(3),growhi(3)
        do j=growlo(2),growhi(2)
        do i=growlo(1),growhi(1)

         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir-1)

         if (dir.eq.1) then
          inorm=i
         else if (dir.eq.2) then
          inorm=j
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          inorm=k
         else
          print *,"dir invalid"
          stop
         endif

         uncoupled_viscosity_override=0

         side_face=0
         if (inorm.eq.fablo(dir)) then
          side_face=1
         else if (inorm.eq.fabhi(dir)+1) then
          side_face=2
         else if ((inorm.gt.fablo(dir)).and. &
                  (inorm.lt.fabhi(dir)+1)) then
          ! do nothing
         else
          print *,"inorm invalid"
          stop
         endif

         do velcomp_alt=1,SDIM
          if (side_face.eq.0) then
           ! do nothing
          else if ((side_face.eq.1).or.(side_face.eq.2)) then
           local_bc=velbc(dir,side_face,velcomp_alt)
           if ((local_bc.eq.EXT_DIR).or. &
               (local_bc.eq.REFLECT_EVEN).or. &
               (local_bc.eq.REFLECT_ODD).or. &
               (local_bc.eq.FOEXTRAP)) then
            uncoupled_viscosity_override=1
           else if (local_bc.eq.INT_DIR) then
            ! do nothing
           else
            print *,"local_bc invalid"
            stop
           endif
          else
           print *,"side_face invalid"
           stop
          endif
         enddo ! velcomp_alt=1..sdim

         im1=i-ii
         jm1=j-jj
         km1=k-kk

         do im=1,2*num_materials
          massF(im)=xface(D_DECL(i,j,k),FACECOMP_MASSFACE+im)
         enddo
         do im=1,num_materials
          massfrac(im)=zero
         enddo
         total_mass=zero
         do side=1,2
          do im=1,num_materials
           DMface=massF(2*(im-1)+side)
           if (DMface.lt.zero) then
            print *,"DMface bust"
            stop
           endif
           total_mass=total_mass+DMface
           massfrac(im)=massfrac(im)+DMface
          enddo ! im
         enddo ! side=1,2
         if (total_mass.gt.zero) then
          do im=1,num_materials
           massfrac(im)=massfrac(im)/total_mass
          enddo
         else if (total_mass.eq.zero) then
          ! do nothing
         else
          print *,"total_mass invalid"
          stop
         endif

         do velcomp=1,SDIM
          diff_flux(velcomp)=zero
         enddo  ! velcomp

         if ((uncoupled_viscosity.eq.0).and. &
             (uncoupled_viscosity_override.eq.0)) then

          do nc=1,SDIM-1

            ! find face stencil
           ilo=i
           ihi=i
           jlo=j
           jhi=j
           klo=k
           khi=k

           if (dir.eq.1) then ! u component  x-face
            ilo=i-1
            if (nc.eq.1) then  ! du/dy
             jhi=j+1
            else if (nc.eq.2) then  ! du/dz
             khi=k+1
            else
             print *,"nc invalid"
             stop
            endif
           else if (dir.eq.2) then  ! v component y-face
            jlo=j-1
            if (nc.eq.1) then ! dv/dx
             ihi=i+1
            else if (nc.eq.2) then ! dv/dz
             khi=k+1
            else
             print *,"nc invalid"
             stop
            endif
           else if ((dir.eq.3).and.(SDIM.eq.3)) then ! w component z-face
            klo=k-1
            if (nc.eq.1) then ! dw/dx
             ihi=i+1
            else if (nc.eq.2) then  ! dw/dy
             jhi=j+1
            else
             print *,"nc invalid"
             stop
            endif
           else
            print *,"dir invalid crossterm 3"
            stop
           endif

           ! mdata=0 if both adjoining cells to a face are solid cells or
           ! a cell pair is outside the grid.
           call slopecrossterm( &
             massfrac, &
             total_mass, &
             levelpc_ptr, &
             faceLS_ptr, &
             mdata_ptr, &
             tdata_ptr, &
             ii,jj,kk, &
             i,j,k,dir, &
             dirtan(nc), &
             coupling(nc), &
             ilo,ihi, &
             jlo,jhi, &
             klo,khi, &
             diff_flux(dirtan(nc)))

          enddo ! nc=1..sdim-1

         else if ((uncoupled_viscosity.eq.1).or. &
                  (uncoupled_viscosity_override.eq.1)) then
          ! do nothing
         else
          print *,"uncoupled_viscosity or uncoupled_viscosity_override invalid"
          stop
         endif
    

! 2 u_x, v_x, w_x or
! u_y, 2 v_y, w_y or
! u_z, v_z, 2 w_z 

         do nc=1,SDIM

          if (uncoupled_viscosity.eq.0) then
           if (nc.eq.dir) then
            alpha=two
           else
            alpha=one
           endif
          else if (uncoupled_viscosity.eq.1) then
           alpha=one
          else
           print *,"uncoupled_viscosity invalid"
           stop
          endif

           ! ux,vx, wx  or
           ! uy,vy, wy  or
           ! uz,vz, wz  
          tcomp=nbase+nc

          local_flux_val=tdata(D_DECL(i,j,k),tcomp)
          if (abs(local_flux_val).le.OVERFLOW_CUTOFF) then
           ! do nothing
          else
           print *,"tdata overflow: i,j,k,nbase,nc,tdata ", &
            i,j,k,nbase,nc,local_flux_val
           print *,"dt,cur_time,level ",dt,cur_time,level
           stop
          endif

          gradterm=alpha*local_flux_val

          if (side_face.eq.0) then
           ! do nothing
          else if ((side_face.eq.1).or.(side_face.eq.2)) then
           local_bc=velbc(dir,side_face,nc)
           if ((local_bc.eq.INT_DIR).or. &
               (local_bc.eq.REFLECT_ODD).or. &
               (local_bc.eq.EXT_DIR)) then
            ! do nothing
           else if ((local_bc.eq.REFLECT_EVEN).or. &
                    (local_bc.eq.FOEXTRAP)) then
            gradterm=zero
           else
            print *,"local_bc invalid"
            stop
           endif
          else
           print *,"side_face invalid"
           stop
          endif
    
          if (uncoupled_viscosity.eq.0) then 
           diff_flux(nc)=diff_flux(nc)+gradterm
          else if (uncoupled_viscosity.eq.1) then
           if (diff_flux(nc).eq.zero) then
            diff_flux(nc)=gradterm
           else
            print *,"diff_flux should be zero"
            stop
           endif
          else
           print *,"uncoupled_viscosity invalid"
           stop
          endif

         enddo  ! nc=1..sdim

         do velcomp=1,SDIM
          xflux(D_DECL(i,j,k),velcomp)=diff_flux(velcomp)
         enddo  ! velcomp

        enddo !i
        enddo !j
        enddo !k faces

       else if (spectral_loop.eq.1) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

      else if (tileloop.eq.1) then ! high order (grad U + grad U^T)

       if (spectral_loop.eq.0) then

         ! overwrite fluxes in spectral elements 
         ! CROSSTERM -> TO_MAC  (interpolate CELL_TO_MAC)
        if (enable_spectral.eq.1) then

         if (bfact.ge.2) then

          call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
          do k=growlo(3),growhi(3)
          do j=growlo(2),growhi(2)
          do i=growlo(1),growhi(1)

           if ((dir.ge.1).and.(dir.le.SDIM)) then
            ! do nothing
           else
            print *,"dir invalid crossterm 4"
            stop
           endif

           if (dir.eq.1) then
            inorm=i
           else if (dir.eq.2) then
            inorm=j
           else if ((dir.eq.3).and.(SDIM.eq.3)) then
            inorm=k
           else
            print *,"dir invalid"
            stop
           endif

           call strip_status(i,j,k,bfact,stripstat)

           if (stripstat.eq.1) then

            local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
            maskcov=NINT(maskcoef(D_DECL(i,j,k)))

            if ((local_maskSEM.ge.1).and. &
                (local_maskSEM.le.num_materials).and. &
                (maskcov.eq.1)) then

              ! elemhi(dir)=elemlo(dir)
             call elementbox(i,j,k,bfact,dir-1,elemlo,elemhi)
             do kelem=elemlo(3),elemhi(3)
             do jelem=elemlo(2),elemhi(2)
             do ielem=elemlo(1),elemhi(1)

              if (dir.eq.1) then ! x-fluxes
               inorm_elem=ielem
              else if (dir.eq.2) then ! y-fluxes
               inorm_elem=jelem
              else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then ! z-fluxes
               inorm_elem=kelem
              else
               print *,"dir invalid crossterm 5"
               stop
              endif

              if ((inorm_elem/bfact)*bfact.ne.inorm_elem) then
               print *,"inorm_elem invalid"
               stop
              endif
              if (inorm_elem.lt.0) then
               print *,"inorm_elem invalid"
               stop
              endif

              do dir2=1,SDIM
               if (fablo(dir2).lt.0) then
                print *,"fablo invalid"
                stop
               endif
              enddo ! dir2

              indexlo(1)=ielem
              indexlo(2)=jelem
              if (SDIM.eq.3) then
               indexlo(SDIM)=kelem
              endif
              do dir2=1,SDIM
               indexhi(dir2)=indexlo(dir2)
              enddo ! dir2
              indexhi(dir)=indexlo(dir)+bfact-1

              do isten=0,bfact
               do velcomp=1,SDIM
                lineflux(isten,velcomp)=zero
               enddo
              enddo

              if (uncoupled_viscosity.eq.0) then

               do nc=1,SDIM-1

                do side=1,2

                 nbr_outside_domain_flag(side)=0

                 do dir2=1,SDIM
                  sideidx(dir2)=indexlo(dir2)
                 enddo

                 if (side.eq.1) then

                  sideidx(dir)=indexlo(dir)-1

                  do dir2=1,SDIM
                   if (sideidx(dir2).lt.fablo(dir2)) then
                    if (velbc(dir2,side,dir2).ne.INT_DIR) then
                     nbr_outside_domain_flag(side)=1
                    endif
                   endif
                  enddo ! dir2

                  i_out=indexlo(1)-ii
                  j_out=indexlo(2)-jj
                  k_out=indexlo(SDIM)-kk

                 else if (side.eq.2) then

                  sideidx(dir)=indexhi(dir)+1
                  do dir2=1,SDIM
                   if (sideidx(dir2).gt.fabhi(dir2)) then
                    if (velbc(dir2,side,dir2).ne.INT_DIR) then
                     nbr_outside_domain_flag(side)=1
                    endif
                   endif
                  enddo ! dir2

                  i_out=indexhi(1)+ii
                  j_out=indexhi(2)+jj
                  k_out=indexhi(SDIM)+kk

                 else 
                  print *,"side invalid"
                  stop
                 endif

                 if (nbr_outside_domain_flag(side).eq.0) then

                  ic=i_out
                  jc=j_out
                  kc=k_out

                   ! 0=covered 1=not covered
                  nbr_covered_flag=NINT(maskcoef(D_DECL(ic,jc,kc)))

                  call gridsten(xsten,xlo, &
                   ic,jc,kc, &
                   fablo,bfact,dx,nhalf)

                  local_bctype(side)=SEM_INTERIOR

                  local_bcval(side)=zero

                  tcomp=coupling(nc)
  
                  if (nbr_covered_flag.eq.1) then 
                   local_data_side(side)=c_tdata(D_DECL(ic,jc,kc),tcomp)
                  else if (nbr_covered_flag.eq.0) then
                   local_data_side(side)=zero
                   local_bctype(side)=SEM_EXTRAP
                   local_bcval(side)=zero
                  else
                   print *,"nbr_covered_flag invalid"
                   stop
                  endif

                 else if (nbr_outside_domain_flag(side).eq.1) then 

                  call gridsten(xsten,xlo, &
                   i_out,j_out,k_out, &
                   fablo,bfact,dx,nhalf)

                  local_data_side(side)=zero

                  if (velbc(dir,side,dir).eq.INT_DIR) then
                   print *,"velbc bust "
                   print *,"side=",side
                   print *,"nbr_outside_domain_flag= ", &
                     nbr_outside_domain_flag(side)
                   stop
                  endif
                   ! currently in: fort_crossterm
                  if ((velbc(dir,side,dir).eq.REFLECT_EVEN).or. &
                      (velbc(dir,side,dir).eq.FOEXTRAP).or. &
                      (velbc(dir,side,dir).eq.REFLECT_ODD).or. &
                      (velbc(dir,side,dir).eq.EXT_DIR)) then
                   local_bctype(side)=SEM_EXTRAP
                   local_bcval(side)=zero
                  else
                   print *,"velbc is corrupt"
                   stop
                  endif
                 else
                  print *,"nbr_outside_domain_flag invalid 1"
                  stop
                 endif

                enddo ! side=1..2

                do isten=0,bfact-1
                 do dir2=1,SDIM
                  indexmid(dir2)=indexlo(dir2)
                 enddo
                 indexmid(dir)=indexlo(dir)+isten
               
                 ic=indexmid(1)
                 jc=indexmid(2)
                 kc=indexmid(SDIM)

                 call gridsten(xsten,xlo, &
                  ic,jc,kc, &
                  fablo,bfact,dx,nhalf)

                 tcomp=coupling(nc)

                 local_data(isten+1)=c_tdata(D_DECL(ic,jc,kc),tcomp)
                enddo !isten=0..bfact-1

                do isten=0,bfact
                 indexmid(dir)=indexlo(dir)+isten
                 ic=indexmid(1)
                 jc=indexmid(2)
                 kc=indexmid(SDIM)

                 call gridstenMAC(xstenMAC,xlo, &
                  ic,jc,kc, &
                  fablo,bfact,dx,nhalf,dir-1)

                 RRface(isten)=xstenMAC(0,1)

                 local_vel(isten)=zero
                enddo ! isten=0..bfact

                call lineGRAD( &
                 levelrz, &
                 dir, &
                 nc, &
                 RRface, &
                 local_bctype, &
                 local_bcval, &
                 local_vel, &
                 local_data, &
                 local_data_side, &
                 local_grad, &
                 local_interp, &
                 bfact, &
                 dx(dir), &
                 x_sep, &
                 operation_flag) ! operation_flag==OP_UGRAD_COUPLING_MAC

                do isten=0,bfact

                 side_cell=0
                 if (inorm+isten.eq.fablo(dir)) then
                  side_cell=1
                 else if (inorm+isten.eq.fabhi(dir)+1) then
                  side_cell=2
                 else if ((inorm+isten.gt.fablo(dir)).and. &
                          (inorm+isten.lt.fabhi(dir)+1)) then
                  ! do nothing
                 else
                  print *,"inorm or isten invalid"
                  stop
                 endif

                 uncoupled_viscosity_override=0
                
                 do velcomp_alt=1,SDIM
                  if (side_cell.eq.0) then
                   ! do nothing
                  else if (((side_cell.eq.1).and.(isten.eq.0)).or. &
                           ((side_cell.eq.2).and.(isten.eq.bfact))) then
                   local_bc=velbc(dir,side_cell,velcomp_alt)
                   if ((local_bc.eq.EXT_DIR).or. &
                       (local_bc.eq.REFLECT_EVEN).or. &
                       (local_bc.eq.REFLECT_ODD).or. &
                       (local_bc.eq.FOEXTRAP)) then
                    uncoupled_viscosity_override=1
                   else if (local_bc.eq.INT_DIR) then
                    ! do nothing
                   else
                    print *,"local_bc invalid"
                    stop
                   endif
                  else if ((isten.gt.0).and.(isten.lt.bfact)) then
                   ! do nothing
                  else
                   print *,"side_cell invalid1"
                   stop
                  endif
                 enddo ! velcomp_alt=1..sdim

                 if (uncoupled_viscosity_override.eq.0) then
                  lineflux(isten,dirtan(nc))= & 
                   lineflux(isten,dirtan(nc))+local_interp(isten)
                 else if (uncoupled_viscosity_override.eq.1) then
                  ! do nothing
                 else
                  print *,"uncoupled_viscosity_override invalid"
                  stop
                 endif
                enddo ! isten=0..bfact

               enddo ! nc=1..sdim-1

              else if (uncoupled_viscosity.eq.1) then
               ! do nothing
              else
               print *,"uncoupled_viscosity invalid"
               stop
              endif
                 
              do isten=0,bfact

               side_cell=0
               if (inorm+isten.eq.fablo(dir)) then
                side_cell=1
               else if (inorm+isten.eq.fabhi(dir)+1) then
                side_cell=2
               else if ((inorm+isten.gt.fablo(dir)).and. &
                        (inorm+isten.lt.fabhi(dir)+1)) then
                ! do nothing
               else
                print *,"inorm or isten invalid"
                stop
               endif

               ! in: crossterm; prevent race condition when doing tiling.
               ! shared_face==1 => two threads would compete for same
               ! point without intervention.
               shared_face=0 

               do dir2=1,SDIM
                indexmid(dir2)=indexlo(dir2)
               enddo
               indexmid(dir)=indexlo(dir)+isten

               ic=indexmid(1)
               jc=indexmid(2)
               kc=indexmid(SDIM)

               call gridstenMAC_level(xstenMAC, &
                ic,jc,kc,level,nhalf,dir-1)

               if (isten.eq.bfact) then ! right side of element

                test_maskSEM=NINT(maskSEM(D_DECL(ic,jc,kc)))

                 ! maskcov=1 => cell is on finest possible level.(not cov)
                 ! maskcov=0 => cell is not on finest possible level.(cov)

                maskcov=NINT(maskcoef(D_DECL(ic,jc,kc)))

                if ((indexmid(dir).gt.fablo(dir)).and. &
                    (indexmid(dir).le.fabhi(dir)).and. &
                    (test_maskSEM.eq.local_maskSEM).and. &
                    (maskcov.eq.1)) then
                 shared_face=1
                else if ((indexmid(dir).eq.fabhi(dir)+1).or. &
                         (test_maskSEM.ne.local_maskSEM).or. &
                         (maskcov.eq.0)) then
                 ! do nothing
                else
                 print *,"indexmid,test_maskSEM, or maskcov invalid"
                 stop
                endif

               else if ((isten.ge.0).and.(isten.lt.bfact)) then
                ! do nothing
               else
                print *,"isten invalid"
                stop
               endif

               do nc=1,SDIM

                if (uncoupled_viscosity.eq.0) then
                 if (nc.eq.dir) then
                  alpha=two
                 else
                  alpha=one
                 endif
                else if (uncoupled_viscosity.eq.1) then
                 alpha=one
                else 
                 print *,"uncoupled_viscosity invalid"
                 stop
                endif

                ! ux,vx, wx  or
                ! uy,vy, wy  or
                ! uz,vz, wz  
                tcomp=nbase+nc
        
                local_flux_val=tdata(D_DECL(ic,jc,kc),tcomp)

                gradterm=alpha*local_flux_val

                if (side_cell.eq.0) then
                 ! do nothing
                else if (((side_cell.eq.1).and.(isten.eq.0)).or. &
                         ((side_cell.eq.2).and.(isten.eq.bfact))) then
                 local_bc=velbc(dir,side_cell,nc)
                 if ((local_bc.eq.INT_DIR).or. &
                     (local_bc.eq.REFLECT_ODD).or. &
                     (local_bc.eq.EXT_DIR)) then
                  ! do nothing
                 else if ((local_bc.eq.REFLECT_EVEN).or. &
                          (local_bc.eq.FOEXTRAP)) then
                  gradterm=zero
                 else
                  print *,"local_bc invalid"
                  stop
                 endif
                else if ((isten.gt.0).and.(isten.lt.bfact)) then
                 ! do nothing
                else
                 print *,"side_cell invalid2"
                 print *,"side_cell=",side_cell
                 print *,"isten,bfact= ",isten,bfact
                 print *,"dir=",dir
                 print *,"nc=",nc
                 print *,"local_bc=",local_bc
                 print *,"int dir= ",INT_DIR
                 print *,"ref odd= ",REFLECT_ODD
                 print *,"ext dir= ",EXT_DIR
                 print *,"ref evn= ",REFLECT_EVEN
                 print *,"fo ext=",FOEXTRAP
                 print *,"i,j,k=",i,j,k
                 print *,"inorm=",inorm
                 print *,"growlo=",growlo(1),growlo(2),growlo(3)
                 print *,"growhi=",growhi(1),growhi(2),growhi(3)
                 print *,"fablo(dir),fabhi(dir) ",fablo(dir),fabhi(dir)
                 print *,"ielem,jelem,kelem ",ielem,jelem,kelem
                 print *,"inorm_elem ",inorm_elem
                 stop
                endif
               
                if (uncoupled_viscosity.eq.0) then 
                 lineflux(isten,nc)=lineflux(isten,nc)+gradterm
                else if (uncoupled_viscosity.eq.1) then
                 if (lineflux(isten,nc).eq.zero) then
                  lineflux(isten,nc)=gradterm
                 else
                  print *,"lineflux should be zero"
                  stop
                 endif
                else
                 print *,"uncoupled_viscosity invalid"
                 stop
                endif

               enddo  ! nc=1..sdim
            
               ! shared_face=1 if right element edge and not the right
               ! edge of the grid or next to a low order or covered
               ! element.
               if (shared_face.eq.0) then

                do velcomp=1,SDIM
                 xflux(D_DECL(ic,jc,kc),velcomp)=lineflux(isten,velcomp)
                enddo  ! velcomp

               else if (shared_face.eq.1) then
                ! do nothing
               else
                print *,"shared_face invalid"
                stop
               endif

              enddo ! isten=0..bfact

             enddo
             enddo
             enddo !ielem,jelem,kelem

            else if ((local_maskSEM.eq.0).or. &
                     (maskcov.eq.0)) then
             ! do nothing
            else
             print *,"local_maskSEM or maskcov invalid"
             stop
            endif 

           else if (stripstat.eq.0) then
            ! do nothing
           else
            print *,"stripstat invalid"
            stop
           endif

          enddo !i
          enddo !j
          enddo !k

         else if (bfact.eq.1) then
          ! do nothing
         else
          print *,"bfact invalid80"
          stop
         endif

        else if (enable_spectral.eq.0) then
         ! do nothing
        else
         print *,"enable_spectral invalid"
         stop
        endif

       else if (spectral_loop.eq.1) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

      ! (1) fluxes+=divu  (2) fluxes*=visccoef
      else if (tileloop.eq.2) then 

       if (spectral_loop.eq.0) then 

          ! does not include the right face of the tile unless
          ! tilehi==fabhi
        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
         growlo,growhi,0,dir-1)

         ! u_x+v_y+w_z on flux face
         ! multiply by visc_constant
        do k=growlo(3),growhi(3)
        do j=growlo(2),growhi(2)
        do i=growlo(1),growhi(1)

          ! dir=1..sdim
         call gridstenMAC_level(xstenMAC, &
           i,j,k,level,nhalf,dir-1)

         im1=i-ii
         jm1=j-jj
         km1=k-kk

         if (dir.eq.1) then
          inorm=i
         else if (dir.eq.2) then
          inorm=j
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          inorm=k
         else
          print *,"dir invalid"
          stop
         endif

         side_face=0
         if (inorm.eq.fablo(dir)) then
          side_face=1
         else if (inorm.eq.fabhi(dir)+1) then
          side_face=2
         else if ((inorm.gt.fablo(dir)).and. &
                  (inorm.lt.fabhi(dir)+1)) then
          ! do nothing
         else
          print *,"inorm invalid"
          stop
         endif

         do im=1,num_materials
          LSleft(im)=levelpc(D_DECL(im1,jm1,km1),im)
          LSright(im)=levelpc(D_DECL(i,j,k),im)
         enddo
         call get_primary_material(LSleft,imL)
         call get_primary_material(LSright,imR)

         if ((imL.lt.1).or.(imL.gt.num_materials)) then
          print *,"imL invalid"
          stop
         endif
         if ((imR.lt.1).or.(imR.gt.num_materials)) then
          print *,"imR invalid"
          stop
         endif

         call gridsten_level(xclamped_minus_sten,im1,jm1,km1,level,nhalf)
         call gridsten_level(xclamped_plus_sten,i,j,k,level,nhalf)
         do dir2=1,SDIM
          xclamped_minus(dir2)=xclamped_minus_sten(0,dir2)
          xclamped_plus(dir2)=xclamped_plus_sten(0,dir2)
         enddo

          ! LS>0 if clamped
         call SUB_clamped_LS(xclamped_minus,cur_time,LS_clamped_minus, &
          vel_clamped_minus,temperature_clamped_minus,prescribed_flag,dx)
         call SUB_clamped_LS(xclamped_plus,cur_time,LS_clamped_plus, &
          vel_clamped_plus,temperature_clamped_plus,prescribed_flag,dx)
         if ((LS_clamped_minus.ge.zero).or. &
             (LS_clamped_plus.ge.zero)) then
          is_clamped_face=1
          do dir2=1,SDIM
           if (LS_clamped_minus.lt.zero) then
            vel_clamped_face(dir2)=vel_clamped_plus(dir2)
            is_clamped_face=2
           else if (LS_clamped_plus.lt.zero) then
            vel_clamped_face(dir2)=vel_clamped_minus(dir2)
            is_clamped_face=3
           else
            vel_clamped_face(dir2)=half*(vel_clamped_plus(dir2)+ &
              vel_clamped_minus(dir2))
           endif
          enddo
         else if ((LS_clamped_minus.lt.zero).and. &
                  (LS_clamped_plus.lt.zero)) then
          is_clamped_face=0
         else
          print *,"LS_clamped plus or minus is NaN"
          stop
         endif

         compressible_face=1

         if ((is_clamped_face.eq.1).or. &
             (is_clamped_face.eq.2).or. &
             (is_clamped_face.eq.3)) then
          compressible_face=0
         else if (is_clamped_face.eq.0) then

          if ((is_rigid(imL).eq.1).or. &
              (is_rigid(imR).eq.1)) then
           compressible_face=0
          else if ((is_rigid(imL).eq.0).and. &
                   (is_rigid(imR).eq.0)) then
           ! do nothing
          else
           print *,"is_rigid invalid GODUNOV_3D.F90"
           stop
          endif
         else
          print *,"is_clamped_face invalid"
          stop
         endif

         if ((fort_material_type(imL).eq.0).or. &
             (fort_material_type(imR).eq.0)) then
          compressible_face=0
         endif

         do velcomp_alt=1,SDIM
          if (side_face.eq.0) then
           ! do nothing
          else if ((side_face.eq.1).or.(side_face.eq.2)) then
           local_bc=velbc(dir,side_face,velcomp_alt)
           if ((local_bc.eq.EXT_DIR).or. &
               (local_bc.eq.REFLECT_EVEN).or. &
               (local_bc.eq.REFLECT_ODD).or. &
               (local_bc.eq.FOEXTRAP)) then
            compressible_face=0
           else if (local_bc.eq.INT_DIR) then
            ! do nothing
           else
            print *,"local_bc invalid: ",local_bc
            stop
           endif
          else
           print *,"side_face invalid: ",side_face
           stop
          endif
         enddo ! velcomp_alt=1..sdim

         if (compressible_face.eq.0) then
          ! do nothing
         else if (compressible_face.eq.1) then

          do im=1,num_materials

           do side=1,2
            if (side.eq.1) then
             local_LS=LSleft(im)
            else if (side.eq.2) then
             local_LS=LSright(im)
            else
             print *,"side invalid (im) : ",im,side
             stop
            endif
            if (local_LS.ge.-incomp_thickness*dxmaxLS) then
             if (is_rigid(im).eq.1) then
              compressible_face=0
             else if (is_rigid(im).eq.0) then
              if (is_compressible_mat(im).eq.0) then
               compressible_face=0
              else if (is_compressible_mat(im).eq.1) then
               !do nothing
              else
               print *,"is_compressible_mat(im) invalid: ", &
                im,is_compressible_mat(im)
               stop
              endif
             else 
              print *,"is_rigid(im) invalid: ",im,is_rigid(im)
              stop
             endif
            else if (local_LS.le.-incomp_thickness*dxmaxLS) then
             ! do nothing
            else
             print *,"local_LS corrupt,fort_crossterm"
             print *,"im,local_LS: ",im,local_LS
             stop
            endif 
           enddo !side=1,2

          enddo !im=1,num_materials

         else 
          print *,"compressible_face invalid: ",compressible_face
          stop
         endif

         wzterm=zero

         if (dir.eq.1) then

          uxterm=tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_UX)

          vyterm=(tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_VY)+ &
                  tdata(D_DECL(i,j+1,k),TENSOR_TRANSPOSE_VY)+ &
                  tdata(D_DECL(i-1,j,k),TENSOR_TRANSPOSE_VY)+ &
                  tdata(D_DECL(i-1,j+1,k),TENSOR_TRANSPOSE_VY))/four
          if (SDIM.eq.3) then
           wzterm= &
            (tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_WZ)+ &
             tdata(D_DECL(i,j,k+1),TENSOR_TRANSPOSE_WZ)+ &
             tdata(D_DECL(i-1,j,k),TENSOR_TRANSPOSE_WZ)+ &
             tdata(D_DECL(i-1,j,k+1),TENSOR_TRANSPOSE_WZ))/four
          endif

         else if (dir.eq.2) then

          uxterm=(tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_UX)+ &
                  tdata(D_DECL(i+1,j,k),TENSOR_TRANSPOSE_UX)+ &
                  tdata(D_DECL(i,j-1,k),TENSOR_TRANSPOSE_UX)+ &
                  tdata(D_DECL(i+1,j-1,k),TENSOR_TRANSPOSE_UX))/four

          vyterm=tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_VY)

          if (SDIM.eq.3) then
           wzterm= &
            (tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_WZ)+ &
             tdata(D_DECL(i,j,k+1),TENSOR_TRANSPOSE_WZ)+ &
             tdata(D_DECL(i,j-1,k),TENSOR_TRANSPOSE_WZ)+ &
             tdata(D_DECL(i,j-1,k+1),TENSOR_TRANSPOSE_WZ))/four
          endif

         else if ((dir.eq.3).and.(SDIM.eq.3)) then

          uxterm=(tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_UX)+ &
                  tdata(D_DECL(i+1,j,k),TENSOR_TRANSPOSE_UX)+ &
                  tdata(D_DECL(i,j,k-1),TENSOR_TRANSPOSE_UX)+ &
                  tdata(D_DECL(i+1,j,k-1),TENSOR_TRANSPOSE_UX))/four

          vyterm=(tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_VY)+ &
                  tdata(D_DECL(i,j+1,k),TENSOR_TRANSPOSE_VY)+ &
                  tdata(D_DECL(i,j,k-1),TENSOR_TRANSPOSE_VY)+ &
                  tdata(D_DECL(i,j+1,k-1),TENSOR_TRANSPOSE_VY))/four

          wzterm=tdata(D_DECL(i,j,k),TENSOR_TRANSPOSE_WZ)

         else
          print *,"dir invalid crossterm 6: ",dir
          stop
         endif

         if (compressible_face.eq.0) then
          divterm=zero
         else if (mdata(D_DECL(i,j,k),dir).eq.zero) then
          divterm=zero
         else if (compressible_face.eq.1) then
          divterm=-(two/SDIM)*(uxterm+vyterm+wzterm)
         else
          print *,"compressible_face bust: ",compressible_face
          stop
         endif

         visc_constant=visc_coef*xface(D_DECL(i,j,k),FACECOMP_FACEVISC+1)
         if (visc_constant.lt.zero) then
          print *,"visc_constant cannot be negative"
          stop
         endif
         visc_constant=-dt*visc_constant

         do velcomp=1,SDIM

          local_flux_val=xflux(D_DECL(i,j,k),velcomp)+divterm

          xflux(D_DECL(i,j,k),velcomp)=visc_constant*local_flux_val

         enddo  ! velcomp

        enddo !i 
        enddo !j
        enddo !k add divterm to fluxes and multiply by visc_constant

       else if (spectral_loop.eq.1) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

       ! average "left" and "right" fluxes at faces separating
       ! 2 spectral elements.
      else if (tileloop.eq.3) then 

       if (enable_spectral.eq.1) then

        if (bfact.ge.2) then

         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
         if ((dir.lt.1).or.(dir.gt.SDIM)) then
          print *,"dir invalid crossterm"
          stop
         endif
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          call strip_status(i,j,k,bfact,stripstat)

          if (stripstat.eq.1) then

           local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
           maskcov=NINT(maskcoef(D_DECL(i,j,k)))

           if ((local_maskSEM.ge.1).and. &
               (local_maskSEM.le.num_materials).and. &
               (maskcov.eq.1)) then

            call elementbox(i,j,k,bfact,dir-1,elemlo,elemhi)
            do kelem=elemlo(3),elemhi(3)
            do jelem=elemlo(2),elemhi(2)
            do ielem=elemlo(1),elemhi(1)

             if (dir.eq.1) then ! x-fluxes
              inorm_elem=ielem
             else if (dir.eq.2) then ! y-fluxes
              inorm_elem=jelem
             else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then ! z-fluxes
              inorm_elem=kelem
             else
              print *,"dir invalid crossterm 7"
              stop
             endif

             if ((inorm_elem/bfact)*bfact.ne.inorm_elem) then
              print *,"inorm_elem invalid"
              stop
             endif
             if (inorm_elem.lt.0) then
              print *,"inorm_elem invalid"
              stop
             endif

             indexlo(1)=ielem
             indexlo(2)=jelem
             if (SDIM.eq.3) then
              indexlo(SDIM)=kelem
             endif
             do dir2=1,SDIM
              index_edge(dir2)=indexlo(dir2)
              index_opp(dir2)=indexlo(dir2)
              index_flux(dir2)=indexlo(dir2)
             enddo ! dir2

             do side=1,2

              ! in: crossterm  avoid race conditions when doing tiling.
              ! shared_face==1 => two threads would compete for same
              ! point without intervention.
              shared_face=0 

              if (side.eq.1) then
               index_flux(dir)=indexlo(dir)
               index_edge(dir)=indexlo(dir)
               index_opp(dir)=index_edge(dir)-1
              else if (side.eq.2) then
               index_flux(dir)=indexlo(dir)+bfact
               index_edge(dir)=indexlo(dir)+bfact-1
               index_opp(dir)=index_edge(dir)+1
              else
               print *,"side invalid"
               stop
              endif

              i_in=index_edge(1)
              j_in=index_edge(2)
              k_in=index_edge(SDIM)

              i_out=index_opp(1)
              j_out=index_opp(2)
              k_out=index_opp(SDIM)

              iflux=index_flux(1)
              jflux=index_flux(2)
              kflux=index_flux(SDIM)

              call gridstenMAC_level(xstenMAC, &
               iflux,jflux,kflux,level,nhalf,dir-1)

              test_maskSEM=NINT(maskSEM(D_DECL(i_out,j_out,k_out)))
              maskcov=NINT(maskcoef(D_DECL(i_out,j_out,k_out)))

              if (side.eq.2) then

               if ((index_edge(dir).ge.fablo(dir)).and. &
                   (index_edge(dir).lt.fabhi(dir)).and. &
                   (test_maskSEM.eq.local_maskSEM).and. &
                   (maskcov.eq.1)) then
                shared_face=1
               else if ((index_edge(dir).eq.fabhi(dir)).or. &
                        (test_maskSEM.ne.local_maskSEM).or. &
                        (maskcov.eq.0)) then
                ! do nothing
               else
                print *,"index_edge invalid"
                stop
               endif

              else if (side.eq.1) then
               ! do nothing
              else
               print *,"side invalid"
               stop
              endif

               ! set mask_out=0 at coarse-fine grid border.
               ! use the (already init) flux from the fine side only.
              mask_out=1
              do dir2=1,SDIM
               if ((index_opp(dir2).lt.fablo(dir2)).or. &
                   (index_opp(dir2).gt.fabhi(dir2))) then
                mask_out=NINT(mask(D_DECL(i_out,j_out,k_out)))
               endif
              enddo

               ! do not average with flux from 
               ! an element that is covered.
              if (maskcov.eq.1) then
               ! do nothing
              else if (maskcov.eq.0) then
               mask_out=0
              else
               print *,"maskcov invalid"
               stop
              endif

               ! set mask_out=0 at high/low order interface.
               ! use the (already init) flux from the high order side only.
              if (test_maskSEM.ne.local_maskSEM) then
               mask_out=0
              endif

               ! shared_face=1 for faces on right side of elements and not
               !  touching the right side of the grid, not touching
               !  a maskSEM==0 element, and not touching a covered element.
               ! shared_face=0 for faces on the left side of elements or for
               !  the face touching the right side of the grid, touching
               !  a maskSEM!=cen_maskSEM element, or touching a covered element.
               ! no need to average fluxes at a "shared_face"
               ! since iface=0 for following element corresponds to
               ! iface=bfact of the previous.
              if (shared_face.eq.1) then
               mask_out=0
              else if (shared_face.eq.0) then
               ! do nothing
              else
               print *,"shared_face invalid"
               stop
              endif

              if (spectral_loop.eq.0) then

               do nc=1,SDIM
                xflux_temp=xflux(D_DECL(iflux,jflux,kflux),nc)
              
                semflux(D_DECL(i_in,j_in,k_in),nbase+nc)=xflux_temp
               enddo ! nc

              else if (spectral_loop.eq.1) then

               if (mask_out.eq.1) then

                do nc=1,SDIM

                 local_flux_val_in=semflux(D_DECL(i_in,j_in,k_in),nbase+nc)
                 local_flux_val_out=semflux(D_DECL(i_out,j_out,k_out),nbase+nc)

                 avgflux(nc)=half*(local_flux_val_in+local_flux_val_out)

                enddo ! nc=1..sdim

                do nc=1,SDIM
                 xflux(D_DECL(iflux,jflux,kflux),nc)=avgflux(nc)
                enddo ! nc=1..sdim
               else if (mask_out.eq.0) then
                ! do nothing
               else
                print *,"mask_out invalid"
                stop
               endif

              else
               print *,"spectral_loop invalid"
               stop
              endif

             enddo ! side=1,2

            enddo
            enddo
            enddo !elem,jelem,kelem

           else if ((local_maskSEM.eq.0).or. &
                    (maskcov.eq.0)) then
            ! do nothing
           else
            print *,"loca_maskSEM_cen or maskcov invalid"
            stop
           endif
          else if (stripstat.eq.0) then
           ! do nothing
          else
           print *,"stripstat invalid"
           stop
          endif

         enddo !i
         enddo !j
         enddo !k

        else if (bfact.eq.1) then
         ! do nothing
        else
         print *,"bfact invalid81"
         stop
        endif

       else if (enable_spectral.eq.0) then
        ! do nothing
       else
        print *,"enable_spectral invalid"
        stop
       endif

      else 
       print *,"tileloop invalid"
       stop
      endif

      return 
      end subroutine fort_crossterm


      ! uu_estdt_max(1..sdim) is max vel in dir.
      ! uu_estdt_max(sdim+1) is max c^2
      ! called from: void NavierStokes::MaxAdvectSpeed when static_flag==0
      subroutine fort_estdt( &
        interface_mass_transfer_model, &
        tid, &
        enable_spectral, &
        AMR_min_phase_change_rate, &
        AMR_max_phase_change_rate, &
        elastic_time, &
        microlayer_substrate, &
        microlayer_angle, &
        microlayer_size, &
        macrolayer_size, &
        reaction_rate, &
        freezing_model, &
        Tanasawa_or_Schrage_or_Kassemi, &
        distribute_from_target, &
        saturation_temp, &
        mass_fraction_id, &
        molar_mass, &
        species_molar_mass, &
        denconst_interface, & !fort_estdt
        denconst_interface_min, & !fort_estdt
        velmac,DIMS(velmac), &
        velcell,DIMS(velcell), &
        solidfab,DIMS(solidfab), &
        den,DIMS(den), &
        vof,DIMS(vof), &
        dist,DIMS(dist), &
        xlo,dx, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        min_stefan_velocity_for_dt, &
        cap_wave_speed, &
        uu_estdt_max, & ! fort_estdt
        dt_min, &
        rzflag, &
        denconst, &
        visc_coef, &
        gravity_reference_wavelen_in, &
        dirnormal, &
        nparts, &
        nparts_def, &
        im_solid_map, &
        material_type, &
        time, &
        shock_timestep, &
        cfl, &
        EILE_flag, &
        level,finest_level) &
      bind(c,name='fort_estdt')

      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use hydrateReactor_module
      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def)
      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: level,finest_level
      real(amrex_real), INTENT(in) :: cfl
      integer, INTENT(in) :: EILE_flag
      real(amrex_real), INTENT(in) :: AMR_min_phase_change_rate(SDIM)
      real(amrex_real), INTENT(in) :: AMR_max_phase_change_rate(SDIM)
      real(amrex_real), INTENT(in) :: elastic_time(num_materials)
      integer, INTENT(in) :: shock_timestep(num_materials)
      integer, INTENT(in) :: material_type(num_materials)
      integer, INTENT(in) :: microlayer_substrate(num_materials)
      real(amrex_real), INTENT(in) :: microlayer_angle(num_materials)
      real(amrex_real), INTENT(in) :: microlayer_size(num_materials)
      real(amrex_real), INTENT(in) :: macrolayer_size(num_materials)
      integer, INTENT(in) :: interface_mass_transfer_model(2*num_interfaces)
      real(amrex_real), INTENT(in) :: reaction_rate(2*num_interfaces)
      real(amrex_real) :: K_f
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      integer, INTENT(in) :: Tanasawa_or_Schrage_or_Kassemi(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
      integer, INTENT(in) :: mass_fraction_id(2*num_interfaces)
      real(amrex_real), INTENT(in) :: molar_mass(num_materials)
      real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var+1)
      real(amrex_real), INTENT(in) :: denconst_interface(num_interfaces)
      real(amrex_real), INTENT(in) :: denconst_interface_min(num_interfaces)
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real) uu_estdt
      real(amrex_real) uu_estdt_core
      real(amrex_real) uu_estdt_phase_change
      real(amrex_real) c_core
      real(amrex_real) cc,cleft,cright
      real(amrex_real) cc_diag,cleft_diag,cright_diag
      integer i,j,k
      integer icell,jcell,kcell
      integer ialt,jalt,kalt
      integer, INTENT(in) :: rzflag
      integer, INTENT(in) :: dirnormal
      integer side,dir2
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(inout) :: uu_estdt_max(SDIM+1) ! fort_estdt
      real(amrex_real), INTENT(inout) :: dt_min
      real(amrex_real) user_tension(num_interfaces)
      real(amrex_real), INTENT(in) :: denconst(num_materials)
      real(amrex_real), INTENT(in) :: visc_coef
      real(amrex_real), INTENT(in) :: gravity_reference_wavelen_in
      real(amrex_real) :: gravity_reference_wavelen

      integer, INTENT(in) :: DIMDEC(velmac)
      integer, INTENT(in) :: DIMDEC(velcell)
      integer, INTENT(in) :: DIMDEC(vof)
      integer, INTENT(in) :: DIMDEC(dist)
      integer, INTENT(in) :: DIMDEC(solidfab)
      integer, INTENT(in) :: DIMDEC(den)
      real(amrex_real), target, INTENT(in) :: velmac(DIMV(velmac))
      real(amrex_real), pointer :: velmac_ptr(D_DECL(:,:,:))
      real(amrex_real), target, INTENT(in) :: &
        velcell(DIMV(velcell),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: velcell_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: &
        solidfab(DIMV(solidfab),nparts_def*SDIM) 
      real(amrex_real), pointer :: solidfab_ptr(D_DECL(:,:,:),:)
       ! den,temp,species
      real(amrex_real), target, INTENT(in) :: &
              den(DIMV(den),num_state_material*num_materials)  
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: &
        vof(DIMV(vof),num_materials*ngeom_raw)
      real(amrex_real), pointer :: vof_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: dist(DIMV(dist),num_materials)
      real(amrex_real), pointer :: dist_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: min_stefan_velocity_for_dt
      real(amrex_real), INTENT(inout) :: &
        cap_wave_speed(num_interfaces) !fort_estdt
      real(amrex_real) hx,hxmac
      real(amrex_real) dthold
      integer ii,jj,kk
      integer im,im_primaryL,im_primaryR
      integer ibase
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) LSleft(num_materials)
      real(amrex_real) LSright(num_materials)
      real(amrex_real) LSsub(num_materials)
      integer im_left_interface,im_right_interface
      integer im_opp
      integer iten
      integer im_source,im_dest
      real(amrex_real) temperature_left,temperature_right
      real(amrex_real) density_left,density_right
      real(amrex_real) internal_energy_left,internal_energy_right
      real(amrex_real) massfrac_parm_left(num_species_var+1)
      real(amrex_real) massfrac_parm_right(num_species_var+1)
      real(amrex_real) gradh
      real(amrex_real) weymouth_factor,weymouth_cfl
      real(amrex_real) dxmin,dxmax,dxmaxLS
      real(amrex_real) den1,den2,visc1,visc2
      integer recompute_wave_speed
      real(amrex_real) uulocal
      real(amrex_real) denjump_scale
      real(amrex_real) denjump_scale_temp
      real(amrex_real) denmax
      real(amrex_real) USTEFAN,USTEFAN_hold
      real(amrex_real) LS1,LS2,Tsrc,Tdst,Dsrc,Ddst,Csrc,Cdst,delta
      real(amrex_real) VOFsrc,VOFdst
      real(amrex_real) LL
      integer velcomp
      integer dcompsrc,dcompdst
      integer tcompsrc,tcompdst
      integer ireverse
      integer ifaceR,jfaceR,kfaceR
      real(amrex_real) uleft,uright
      real(amrex_real) C_w0,PHYDWATER,Cmethane_in_hydrate
      integer local_freezing_model
      integer local_Tanasawa_or_Schrage_or_Kassemi
      integer distribute_from_targ
      integer vofcompsrc,vofcompdst
      real(amrex_real) TSAT,Tsrcalt,Tdstalt
      real(amrex_real) uleftcell,urightcell,udiffcell,umaxcell
      real(amrex_real) velsum
      real(amrex_real) RR
      real(amrex_real) level_cap_wave_speed(num_interfaces) !fort_estdt
      real(amrex_real) ksource_predict,kdest_predict
      real(amrex_real) ksource_physical,kdest_physical
      real(amrex_real) alpha,beta,dt_heat
      integer for_estdt
      real(amrex_real) xI(SDIM)
      real(amrex_real) mu
      integer partid
      integer ispec
      real(amrex_real) vapor_den
      real(amrex_real) elastic_wave_speed
      real(amrex_real) source_perim_factor
      real(amrex_real) dest_perim_factor
      real(amrex_real) effective_velocity
      real(amrex_real) local_elastic_time
      real(amrex_real) ugrav
      real(amrex_real) local_gravity_mag
      real(amrex_real) local_temperature

      integer istenlo(3),istenhi(3)
      integer ivec(3)
      integer triple_flag
      integer im_sub

      if ((tid.lt.0).or. &
          (tid.ge.geom_nthreads)) then
       print *,"tid invalid in fort_estdt"
       stop
      endif

      velmac_ptr=>velmac
      velcell_ptr=>velcell
      solidfab_ptr=>solidfab
      den_ptr=>den
      vof_ptr=>vof
      dist_ptr=>dist

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_estdt"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_estdt"
       stop
      endif
      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.1)) then
       print *,"enable_spectral invalid fort_estdt"
       stop
      endif

      vapor_den=one

      denjump_scale=zero

       ! dxmin=min_d min_i dxsub_{gridtype,d,i} d=1..sdim  i=0..bfact-1
       ! gridtype=MAC or CELL
       ! if cylindrical coordinates, then dx_{\theta}*=problox
      call get_dxmin(dx,bfact,dxmin)
      if (dxmin.gt.zero) then
       ! do nothing
      else
       print *,"dxmin must be positive: ",dxmin
       stop
      endif

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif

      if ((finest_level.lt.0).or. &
          (level.lt.0).or. &
          (level.gt.finest_level)) then
       print *,"level or finest level invalid estdt"
       stop
      endif

      if ((EILE_flag.eq.-1).or. & ! Weymouth and Yue
          (EILE_flag.eq.1).or.  & ! EILE
          (EILE_flag.eq.2).or.  & ! always EI
          (EILE_flag.eq.3)) then  ! always LE
       ! do nothing
      else 
       print *,"EILE flag invalid"
       stop
      endif

      if ((cfl.gt.zero).and. &
          (cfl.le.0.95d0)) then
       ! do nothing
      else
       print *,"cfl invalid (fort_estdt): ",cfl
       stop
      endif

      if ((EILE_flag.eq.1).or. & ! EI-LE
          (EILE_flag.eq.2).or. & ! always EI
          (EILE_flag.eq.3)) then ! always LE
       weymouth_cfl=half  ! we advect half cells.
      else if (EILE_flag.eq.-1) then ! Weymouth and Yue
       weymouth_cfl=one/(two*SDIM)
      else
       print *,"EILE_flag invalid"
       stop
      endif

      if ((cfl.gt.zero).and. &
          (cfl.le.0.95d0)) then
       weymouth_factor=weymouth_cfl/cfl
      else
       print *,"cfl invalid (fort_estdt): ",cfl
       print *,"weymouth_cfl ",weymouth_cfl
       stop
      endif

      do im=1,num_materials 

       if ((shock_timestep(im).ne.0).and. &
           (shock_timestep(im).ne.1).and. &
           (shock_timestep(im).ne.2)) then
        print *,"shock_timestep invalid"
        stop
       endif
       if (denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"denconst invalid"
        stop
       endif
       mu=get_user_viscconst(im,fort_denconst(im),fort_tempconst(im))
       if (mu.ge.zero) then
        ! do nothing
       else
        print *,"viscconst invalid"
        stop
       endif

       if (is_rigid(im).eq.0) then
        ! do nothing
       else if (is_rigid(im).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid GODUNOV_3D.F90"
        stop
       endif

      enddo  ! im=1..num_materials

      call get_max_user_tension(fort_tension,user_tension)

         ! finest_level is first level tried.
      recompute_wave_speed=0
      if (level.eq.finest_level) then
       recompute_wave_speed=1
      endif
      
      if (recompute_wave_speed.eq.1) then

       do im=1,num_materials-1
        do im_opp=im+1,num_materials
         if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
          print *,"im or im_opp bust 6"
          stop
         endif
         call get_iten(im,im_opp,iten)
         if ((is_rigid(im).eq.1).or. &
             (is_rigid(im_opp).eq.1)) then
          cap_wave_speed(iten)=zero
         else if ((is_ice_or_FSI_rigid_material(im).eq.1).or. &
                  (is_ice_or_FSI_rigid_material(im_opp).eq.1)) then
          cap_wave_speed(iten)=zero
         else if (user_tension(iten).eq.zero) then
          cap_wave_speed(iten)=zero
         else if (user_tension(iten).gt.zero) then
          den1=denconst(im)
          den2=denconst(im_opp)
          mu=get_user_viscconst(im,den1,fort_tempconst(im))
          visc1=visc_coef*mu+EPS10
          mu=get_user_viscconst(im_opp,den2,fort_tempconst(im_opp))
          visc2=visc_coef*mu+EPS10

          if (denconst_interface(iten).eq.zero) then
           ! do nothing
          else if (denconst_interface(iten).gt.zero) then
           den1=denconst_interface(iten)
           den2=den1
          else
           print *,"denconst_interface invalid fort_estdt"
           stop
          endif

          if (denconst_interface_min(iten).eq.zero) then
           ! do nothing
          else if (denconst_interface_min(iten).gt.zero) then
           if (min(den1,den2).lt.denconst_interface_min(iten)) then
            den1=denconst_interface_min(iten)
            den2=den1
           endif
          else
           print *,"denconst_interface_min invalid fort_estdt"
           stop
          endif


           ! capillary_wave_speed declared in PROB.F90
           ! theory:
           ! wavespeed=sqrt(2\pi tension/((den1+den2)*dx)
           ! practical:
           ! wavespeed=sqrt(2\pi tension/(min(den1,den2)*dx)
           ! =>
           ! sqrt( (N/m) (m^3/kg)(1/m) )=sqrt( (kg/s^2)(m^2/kg) )=m/s
           !
           ! alternate derivation:
           ! ( dt tension kappa grad H/rho )dt < dx
           ! dt^2 < rho dx^3/tension
           ! dt < \sqrt(rho/tension)dx^{3/2} \equiv dx/U
           ! U=sqrt(tension/(rho dx))

          call capillary_wave_speed( &
           dxmin, & !wavelen
           den1,den2, & 
           visc1,visc2, &
           user_tension(iten), &
           cap_wave_speed(iten)) !INTENT(out)

         else
          print *,"user_tension invalid fort_estdt"
          stop
         endif
         level_cap_wave_speed(iten)=cap_wave_speed(iten)
        enddo ! im_opp=im+1..num_materials
       enddo ! im=1..num_materials-1
      else if (recompute_wave_speed.eq.0) then
       do im=1,num_materials-1
        do im_opp=im+1,num_materials
         call get_iten(im,im_opp,iten)
         level_cap_wave_speed(iten)=zero
        enddo
       enddo
      else
       print *,"recompute wave speed invalid fort_estdt"
       print *,"num_materials=",num_materials
       print *,"level=",level
       print *,"finest_level=",finest_level
       print *,"dxmin= ",dxmin
       print *,"cfl= ",cfl
       print *,"denconst(1)= ",denconst(1)
       print *,"get_user_viscconst(1)= ", &
        get_user_viscconst(1,fort_denconst(1),fort_tempconst(1))
       print *,"visc_coef= ",visc_coef
       print *,"user_tension(1)= ",user_tension(1)
       print *,"EILE_flag= ",EILE_flag
       print *,"num_interfaces=",num_interfaces
       print *,"dirnormal=",dirnormal
       print *,"rzflag=",rzflag
       stop
      endif

      if (rzflag.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (rzflag.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust fort_estdt"
        stop
       endif
      else if (rzflag.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"rzflag invalid fort_estdt"
       stop
      endif

      if ((dirnormal.lt.0).or.(dirnormal.ge.SDIM)) then
       print *,"dirnormal invalid fort_estdt"
       stop
      endif
      if ((adv_dir.lt.1).or.(adv_dir.gt.2*SDIM+1)) then
       print *,"adv_dir invalid fort_estdt"
       stop
      endif

      call checkbound_array1(fablo,fabhi,velmac_ptr,0,dirnormal)
      call checkbound_array(fablo,fabhi,velcell_ptr,1,-1)
      call checkbound_array(fablo,fabhi,solidfab_ptr,0,dirnormal)
      call checkbound_array(fablo,fabhi,den_ptr,1,-1)
      call checkbound_array(fablo,fabhi,vof_ptr,1,-1)
       ! need enough ghost cells for the calls to derive_dist.
      call checkbound_array(fablo,fabhi,dist_ptr,2,-1)

      if (rzflag.ne.levelrz) then
       print *,"rzflag invalid fort_estdt"
       stop
      endif

      local_gravity_mag=gravity_vector(1)**2+ &
        gravity_vector(2)**2
      if (SDIM.eq.3) then
       local_gravity_mag=local_gravity_mag+gravity_vector(SDIM)**2
      endif
      local_gravity_mag=sqrt(local_gravity_mag)

       ! get_max_denjump_scale is declared in: PROB.F90
       ! denjump_scale=max_{im,im_opp (fluids)} 
       !   |den(im)-den(im_opp)|/max(den(im),den(im_opp),den_added)
       ! denjump_scale in [0,1]
      call get_max_denjump_scale(denjump_scale,denconst_interface)

      if ((denjump_scale.ge.zero).and.(denjump_scale.le.one)) then
       ! do nothing
      else
       print *,"denjump_scale invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dirnormal.eq.0) then
        ii=1
      else if (dirnormal.eq.1) then
        jj=1
      else if ((dirnormal.eq.2).and.(SDIM.eq.3)) then
        kk=1
      else
       print *,"dirnormal invalid estdt 2"
       stop
      endif

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi, &
              0,dirnormal)

      if (1.eq.0) then
       print *,"dt_min before estdt loop: ",dt_min
       print *,"weymouth_factor= ",weymouth_factor
      endif

      uu_estdt_core = zero
      c_core = zero  ! max of c^2

! if rz and dirnormal=0 and u>0, need u dt r/(r-u dt) < dx
! u dt r < dx(r-u dt)
! u dt dx + udt r < dx r
! u dt (dx/r+1) < dx
!
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       ivec(1)=i
       ivec(2)=j
       if (SDIM.eq.3) then
        ivec(SDIM)=k
       endif

       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dirnormal)
       hx=xstenMAC(1,dirnormal+1)-xstenMAC(-1,dirnormal+1)

       RR=one
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        if (dirnormal.eq.1) then ! theta direction
         RR=xstenMAC(0,1)
        endif
       else
        print *,"levelrz invalid"
        stop
       endif
       if (RR.gt.zero) then
        ! do nothing
       else
        print *,"RR invalid fort_estdt"
        stop
       endif
 
       hx=hx*RR
       hxmac=hx

       if (hx.gt.(one-EPS2)*dxmin) then
        ! do nothing
       else
        print *,"expecting hx>(1-EPS2)*dxmin"
        print *,"xstenMAC invalid estdt"
        print *,"hx= ",hx
        print *,"dxmin= ",dxmin
        print *,"i,j,k ",i,j,k
        print *,"dirnormal=",dirnormal
        print *,"xright ",xstenMAC(1,dirnormal+1)
        print *,"xleft ",xstenMAC(-1,dirnormal+1)
        stop
       endif

       if (enable_spectral.eq.1) then
        if (bfact.ge.2) then
         hx=dxmin
        else if (bfact.eq.1) then
         ! do nothing
        else
         print *,"bfact invalid63"
         stop
        endif
       else if (enable_spectral.eq.0) then
        ! do nothing
       else
        print *,"enable_spectral invalid"
        stop
       endif

       uu_estdt=zero
       uu_estdt_phase_change=zero

        ! first check that characteristics do not collide or that
        ! a cell does not become a vacuum.
       ifaceR=i+ii
       jfaceR=j+jj
       kfaceR=k+kk
       if ((ifaceR.le.growhi(1)).and. &
           (jfaceR.le.growhi(2)).and. &
           (kfaceR.le.growhi(3))) then
        uleft=velmac(D_DECL(i,j,k))
        uright=velmac(D_DECL(ifaceR,jfaceR,kfaceR))
        if (uleft*uright.le.zero) then
         uu_estdt=max(uu_estdt,abs(uleft-uright))
        else if (uleft*uright.gt.zero) then
         !do nothing
        else
         print *,"uleft*uright invalid"
         stop
        endif
       endif
      
       uleftcell=velcell(D_DECL(i-ii,j-jj,k-kk),dirnormal+1)
       urightcell=velcell(D_DECL(i,j,k),dirnormal+1)
       uu_estdt=max(uu_estdt,abs(uleftcell))
       uu_estdt=max(uu_estdt,abs(urightcell))

       if (enable_spectral.eq.1) then
        if (bfact.ge.2) then
         velsum=zero
         do dir2=1,SDIM
          uleftcell=velcell(D_DECL(i-ii,j-jj,k-kk),dir2)
          urightcell=velcell(D_DECL(i,j,k),dir2)
          udiffcell=abs(uleftcell-urightcell)
          umaxcell=abs(uleftcell)
          if (abs(urightcell).gt.umaxcell) then
           umaxcell=abs(urightcell)
          endif
          if (abs(udiffcell).gt.umaxcell) then
           umaxcell=abs(udiffcell)
          endif
          velsum=velsum+umaxcell
         enddo ! dir2=1..sdim
         uu_estdt=max(uu_estdt,velsum)
        else if (bfact.eq.1) then
         ! do nothing
        else
         print *,"bfact invalid64"
         stop
        endif
       else if (enable_spectral.eq.0) then
        ! do nothing
       else
        print *,"enable_spectral invalid"
        stop
       endif

       do partid=0,nparts-1
        velcomp=partid*SDIM+dirnormal+1
        uu_estdt=max(uu_estdt,abs(solidfab(D_DECL(i,j,k),velcomp)))
       enddo ! partid=0..nparts-1

       uulocal=abs(velmac(D_DECL(i,j,k)))

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if ((levelrz.eq.COORDSYS_RZ).or. &
                (levelrz.eq.COORDSYS_CYLINDRICAL)) then
        if ((dirnormal.eq.0).and. &
            (xstenMAC(0,1).gt.EPS2*dx(1))) then
         if (xstenMAC(0,1).ge.hxmac) then
          uulocal=uulocal/( one-three*hxmac/(four*xstenMAC(0,1)) )  
         else if ((xstenMAC(0,1).gt.zero).and. &
                  (xstenMAC(-1,1).gt.zero)) then
          uulocal=four*uulocal
         else
          print *,"xstenMAC invalid estdt 2"
          print *,"i,j,k,dirnormal ",i,j,k,dirnormal
          print *,"hx=",hx
          print *,"hxmac=",hxmac
          print *,"xstenMAC(0,1)= ",xstenMAC(0,1)
          stop
         endif
        else if ((dirnormal.gt.0).or. &
                 (xstenMAC(0,1).le.EPS2*dx(1))) then
         !do nothing
        else
         print *,"dirnormal or xstenMAC invalid"
         stop
        endif
       else
        print *,"levelrz invalid estdt"
        stop
       endif
       uu_estdt=max(uu_estdt,uulocal)

       cleft=zero
       cright=zero
       cc=zero

       cleft_diag=zero
       cright_diag=zero
       cc_diag=zero

       do im=1,num_materials
        if (fort_denconst(im).gt.zero) then
         if (fort_viscosity_state_model(im).ge.0) then
          if (visc_coef.ge.zero) then
           if (elastic_time(im).ge.zero) then
            if (elastic_time(im).eq.zero) then
             local_elastic_time=one
            else if (elastic_time(im).ge.one) then
             local_elastic_time=one
            else
             local_elastic_time=elastic_time(im)
            endif

             ! rho u_t = div beta (grad X + grad X^T)/2
             ! kg/m^3  m/s^2  = (1/m^2) beta m
             ! kg/(m^2 s^2) = (1/m) beta
             ! beta = kg/(m s^2)
             ! beta/rho = kg/(m s^2)   / (kg/m^3) = m^2/s^2
            if (fort_elastic_viscosity(im).ge.zero) then
             elastic_wave_speed=visc_coef*fort_elastic_viscosity(im)/ &
                (local_elastic_time*fort_denconst(im))
            else
             print *,"fort_elastic_viscosity(im) invalid"
             stop
            endif
            if (elastic_wave_speed.gt.zero) then
             elastic_wave_speed=sqrt(elastic_wave_speed)
             dthold=hx/elastic_wave_speed
             dt_min=min(dt_min,dthold)
            else if (elastic_wave_speed.eq.zero) then
             ! do nothing
            else
             print *,"elastic_wave_speed invalid"
             stop
            endif
           else
            print *,"elastic_time invalid"
            print *,"im= ",im
            print *,"elastic_time(im)=",elastic_time(im)
            stop
           endif
          else
           print *,"visc_coef invalid"
           stop
          endif
         else
          print *,"fort_viscosity_state_model invalid"
          stop
         endif
        else
         print *,"fort_denconst(im) invalid"
         stop
        endif
       enddo ! im=1..num_materials

       do im=1,num_materials
        LSleft(im)=dist(D_DECL(i-ii,j-jj,k-kk),im)
        LSright(im)=dist(D_DECL(i,j,k),im)
       enddo
       call get_primary_material(LSleft,im_primaryL)
       call get_primary_material(LSright,im_primaryR)

       USTEFAN=zero

        ! fluid region
       if ((is_rigid(im_primaryL).eq.0).and. &
           (is_rigid(im_primaryR).eq.0)) then 

        do ireverse=0,1
        do im=1,num_materials-1
        do im_opp=im+1,num_materials

         if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
          print *,"im or im_opp bust 7"
          stop
         endif
         call get_iten(im,im_opp,iten)

         tcompsrc=(im_primaryL-1)*num_state_material+1+ENUM_TEMPERATUREVAR
         tcompdst=(im_primaryR-1)*num_state_material+1+ENUM_TEMPERATUREVAR
         Tsrc=den(D_DECL(i-ii,j-jj,k-kk),tcompsrc)
         Tdst=den(D_DECL(i,j,k),tcompdst)

         local_temperature=half*(Tsrc+Tdst)

         LL=get_user_latent_heat(iten+ireverse*num_interfaces, &
                 local_temperature,0)

         K_f=reaction_rate(iten+ireverse*num_interfaces)
         local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
         local_Tanasawa_or_Schrage_or_Kassemi= &
           Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces)
         distribute_from_targ= &
             distribute_from_target(iten+ireverse*num_interfaces)
         TSAT=saturation_temp(iten+ireverse*num_interfaces)

         if (is_hydrate_freezing_modelF(local_freezing_model).eq.1) then
          if (num_species_var.eq.1) then
           ! do nothing
          else
           print *,"num_species_var invalid for hydrate"
           stop
          endif 
         else if (is_valid_freezing_modelF(local_freezing_model).eq.1) then
          ! do nothing
         else
          print *,"local_freezing_model invalid"
          stop
         endif 

         if (ireverse.eq.0) then
          im_source=im
          im_dest=im_opp
         else if (ireverse.eq.1) then
          im_source=im_opp
          im_dest=im
         else
          print *,"ireverse invalid"
          stop
         endif
         dcompsrc=(im_source-1)*num_state_material+1+ENUM_DENVAR
         tcompsrc=(im_source-1)*num_state_material+1+ENUM_TEMPERATUREVAR
         vofcompsrc=(im_source-1)*ngeom_raw+1 
         dcompdst=(im_dest-1)*num_state_material+1+ENUM_DENVAR
         tcompdst=(im_dest-1)*num_state_material+1+ENUM_TEMPERATUREVAR
         vofcompdst=(im_dest-1)*ngeom_raw+1 

         ispec=mass_fraction_id(iten+ireverse*num_interfaces)

         if (is_multi_component_evapF(local_freezing_model, &
              local_Tanasawa_or_Schrage_or_Kassemi,LL).eq.1) then 
          if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
           ! do nothing
          else if (ispec.eq.0) then
           print *,"ispec invalid"
           stop
          endif
         else if (is_multi_component_evapF(local_freezing_model, &
              local_Tanasawa_or_Schrage_or_Kassemi,LL).eq.0) then 
          if (ispec.eq.0) then
           ! do nothing
          else
           print *,"expecting ispec=0"
           stop
          endif
         else
          print *,"is_multi_component_evapF invalid"
          stop
         endif

         if ((is_rigid(im).eq.1).or. &
             (is_rigid(im_opp).eq.1)) then
          ! do nothing
         else if (LL.ne.zero) then

          do side=1,2
           if (side.eq.1) then
            icell=i-ii
            jcell=j-jj
            kcell=k-kk
            ialt=i
            jalt=j
            kalt=k
           else if (side.eq.2) then
            ialt=i-ii
            jalt=j-jj
            kalt=k-kk
            icell=i
            jcell=j
            kcell=k
           else
            print *,"side invalid"
            stop
           endif

           VOFsrc=vof(D_DECL(icell,jcell,kcell),vofcompsrc)
           VOFdst=vof(D_DECL(icell,jcell,kcell),vofcompdst)

           call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)

            ! adjust LS1 if R-theta.
           call derive_dist(xsten,nhalf, &
            dist_ptr,icell,jcell,kcell,im_source,LS1)
            ! adjust LS2 if R-theta.
           call derive_dist(xsten,nhalf, &
            dist_ptr,icell,jcell,kcell,im_dest,LS2)
           
           Tsrc=den(D_DECL(icell,jcell,kcell),tcompsrc)
           Tdst=den(D_DECL(icell,jcell,kcell),tcompdst)
           Tsrcalt=den(D_DECL(ialt,jalt,kalt),tcompsrc)
           Tdstalt=den(D_DECL(ialt,jalt,kalt),tcompdst)
           Dsrc=den(D_DECL(icell,jcell,kcell),dcompsrc)
           Ddst=den(D_DECL(icell,jcell,kcell),dcompdst)

           if (LL.gt.zero) then ! evaporation
            vapor_den=Ddst
           else if (LL.lt.zero) then ! condensation
            vapor_den=Dsrc
           else
            print *,"LL invalid"
            stop
           endif

           if (vapor_den.gt.zero) then
            ! do nothing
           else
            print *,"vapor_den invalid"
            stop
           endif  

           Csrc=zero
           Cdst=zero

           if (is_hydrate_freezing_modelF(local_freezing_model).eq.1) then 
            if (distribute_from_targ.ne.0) then
             print *,"distribute_from_targ invalid"
             stop
            endif
            Csrc=den(D_DECL(icell,jcell,kcell),tcompsrc+1)
            Cdst=den(D_DECL(icell,jcell,kcell),tcompdst+1)
           else if (is_hydrate_freezing_modelF(local_freezing_model).eq.0) then 
            ! do nothing
           else
            print *,"is_hydrate_freezing_modelF invalid"
            stop
           endif

           if ((Dsrc.gt.zero).and.(Ddst.gt.zero)) then
            ! do nothing
           else
            print *,"density must be positive estdt "
            print *,"im,im_opp ",im,im_opp
            print *,"im_source,im_dest ",im_source,im_dest
            print *,"Dsrc(source),Ddst(dest) ",Dsrc,Ddst
            stop
           endif

           if ((abs(LS1).le.two*dxmaxLS).and. &
               (abs(LS2).le.two*dxmaxLS)) then 

            delta=dxmin

            if ((LS1.ge.zero).or. &
                (LS2.ge.zero)) then

              ! coordinate of center of cell adjacent to face (i,j,k)
             do dir2=1,SDIM
              xI(dir2)=xsten(0,dir2)
             enddo
              ! either: 1-den_dst/den_src
              !     or: 1-den_src/den_dst
             if (fort_expansion_factor(iten+ireverse*num_interfaces).ge. &
                 one) then
              print *,"fort_expansion_factor(iten+ireverse*num_interfaces) bad"
              stop
             endif

             ksource_predict=get_user_heatviscconst(im_source)* &
                     fort_heatflux_factor(im_source)
             ksource_physical=get_user_heatviscconst(im_source)

             if ((ksource_predict.ge.zero).and. &
                 (ksource_physical.ge.zero)) then
              ! do nothing
             else
              print *,"ksource_predict or ksource_physical invalid"
              stop
             endif

             kdest_predict=get_user_heatviscconst(im_dest)* &
                     fort_heatflux_factor(im_dest)
             kdest_physical=get_user_heatviscconst(im_dest)

             if ((kdest_predict.ge.zero).and. &
                 (kdest_physical.ge.zero)) then
              ! do nothing
             else
              print *,"kdest_predict or kdest_physical invalid"
              stop
             endif

             alpha=fort_alpha(iten+ireverse*num_interfaces)
             beta=fort_beta(iten+ireverse*num_interfaces)
             if ((alpha.lt.zero).or.(beta.lt.zero)) then
              print *,"alpha or beta are invalid"
              stop
             endif
             dt_heat=one

             C_w0=fort_denconst(1)  ! density of water
             PHYDWATER=2.0D+19
             Cmethane_in_hydrate=Cdst ! hydrate is destination material.
       
             for_estdt=1

             source_perim_factor=one
             dest_perim_factor=one

             call get_vel_phasechange( &
              interface_mass_transfer_model(iten+ireverse*num_interfaces), &
              for_estdt, &
              xI, &
              ispec, &
              molar_mass, &
              species_molar_mass, &
              local_freezing_model, &
              local_Tanasawa_or_Schrage_or_Kassemi, &
              distribute_from_targ, &
              USTEFAN_hold, &
              Dsrc,Ddst, &
              Dsrc,Ddst, &
              ksource_predict, &
              kdest_predict, &
              ksource_physical, &
              kdest_physical, &
              Tsrc,Tdst, &
              TSAT, &
              Tsrcalt,Tdstalt, &
              LL, &
              source_perim_factor, &
              dest_perim_factor, &
              microlayer_substrate(im_source), &
              microlayer_angle(im_source), &
              microlayer_size(im_source), &
              macrolayer_size(im_source), &
              microlayer_substrate(im_dest), &
              microlayer_angle(im_dest), &
              microlayer_size(im_dest), &
              macrolayer_size(im_dest), &
              delta, &  ! dxprobe_source
              delta, &  ! dxprobe_dest
              im_source,im_dest, &
              time,dt_heat, &
              alpha, &
              beta, &
              fort_expansion_factor(iten+ireverse*num_interfaces), &
              K_f, &
              Cmethane_in_hydrate, &
              C_w0, &
              PHYDWATER, &
              VOFsrc,VOFdst)

             if (USTEFAN_hold.ge.zero) then
              ! do nothing
             else
              print *,"USTEFAN_hold.lt.zero"
              stop
             endif

             USTEFAN=USTEFAN+USTEFAN_hold

            else if ((LS1.lt.zero).and.(LS2.lt.zero)) then
             ! do nothing
            else
             print *,"LS1 or LS2 invalid"
             stop
            endif

           endif ! LS1,LS2 both close to 0
 
          enddo ! side 
         endif ! latent_heat<>0
        enddo ! im_opp
        enddo ! im
        enddo ! ireverse=0,1

       else if ((is_rigid(im_primaryL).eq.1).or. &
                (is_rigid(im_primaryR).eq.1)) then
        ! do nothing
       else
        print *,"im_primaryL, or im_primaryR invalid"
        stop 
       endif

       if (time.eq.zero) then
        ! do nothing
       else if (time.gt.zero) then
        USTEFAN=zero
        do dir2=1,SDIM
         if (USTEFAN.lt.abs(AMR_min_phase_change_rate(dir2))) then
          USTEFAN=abs(AMR_min_phase_change_rate(dir2))
         endif
         if (USTEFAN.lt.abs(AMR_max_phase_change_rate(dir2))) then
          USTEFAN=abs(AMR_max_phase_change_rate(dir2))
         endif
        enddo ! dir2=1..sdim
       else
        print *,"time invalid in fort_estdt"
        stop
       endif

        ! factor of 4 in order to guarantee that characteristics do not
        ! collide.
        ! also, the factor of 4 should guarantee that a swept cell is not
        ! full at the end of fort_convertmaterial.
       if (min_stefan_velocity_for_dt.ge.zero) then
        if (USTEFAN.lt.min_stefan_velocity_for_dt) then
         USTEFAN=min_stefan_velocity_for_dt
        else if (USTEFAN.ge.min_stefan_velocity_for_dt) then
         ! do nothing
        else
         print *,"USTEFAN bust"
         stop
        endif
       else
        print *,"min_stefan_velocity_for_dt invalid"
        stop
       endif

       uu_estdt_phase_change=abs(uu_estdt_phase_change)+two*USTEFAN

       if (is_rigid(im_primaryL).eq.0) then
        ibase=(im_primaryL-1)*num_state_material
        density_left= &
           den(D_DECL(i-ii,j-jj,k-kk),ibase+ENUM_DENVAR+1)

        if (density_left.gt.zero) then
         ! do nothing
        else
         print *,"density_left invalid"
         stop
        endif

        if (material_type(im_primaryL).gt.0) then

         call init_massfrac_parm(density_left,massfrac_parm_left,im_primaryL)
         do ispec=1,num_species_var
          massfrac_parm_left(ispec)= &
           den(D_DECL(i-ii,j-jj,k-kk),ibase+ENUM_SPECIESVAR+ispec)
         enddo

         temperature_left=den(D_DECL(i-ii,j-jj,k-kk), &
                 ibase+ENUM_TEMPERATUREVAR+1)
         call INTERNAL_material(density_left,massfrac_parm_left, &
          temperature_left, &
          internal_energy_left, &
          material_type(im_primaryL),im_primaryL)
         call SOUNDSQR_material(density_left,massfrac_parm_left, &
          internal_energy_left, &
          cleft_diag, &
          material_type(im_primaryL),im_primaryL)

         if ((shock_timestep(im_primaryL).eq.1).or. &
             ((shock_timestep(im_primaryL).eq.0).and.(time.eq.zero))) then
          cleft=cleft_diag
         endif
        else if (material_type(im_primaryL).eq.0) then
         ! do nothing
        else
         print *,"material_type(im_primaryL) invalid"
         stop
        endif 
       else if (is_rigid(im_primaryL).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid GODUNOV_3D.F90"
        stop
       endif  

       if (is_rigid(im_primaryR).eq.0) then
        ibase=(im_primaryR-1)*num_state_material
        density_right=den(D_DECL(i,j,k),ibase+ENUM_DENVAR+1)

        if (density_right.gt.zero) then
         ! do nothing
        else
         print *,"density_right invalid"
         stop
        endif

        if (material_type(im_primaryR).gt.0) then

         call init_massfrac_parm(density_right,massfrac_parm_right,im_primaryR)
         do ispec=1,num_species_var
          massfrac_parm_right(ispec)= &
             den(D_DECL(i,j,k),ibase+ENUM_SPECIESVAR+ispec)
         enddo

         temperature_right=den(D_DECL(i,j,k),ibase+ENUM_TEMPERATUREVAR+1)
         call INTERNAL_material(density_right,massfrac_parm_right, &
          temperature_right, &
          internal_energy_right, &
          material_type(im_primaryR),im_primaryR)
         call SOUNDSQR_material(density_right,massfrac_parm_right, &
          internal_energy_right, &
          cright_diag, &
          material_type(im_primaryR),im_primaryR)

         if ((shock_timestep(im_primaryR).eq.1).or. &
             ((shock_timestep(im_primaryR).eq.0).and.(time.eq.zero))) then
          cright=cright_diag
         endif

         if (im_primaryR.eq.im_primaryL) then
          denmax=max(density_left,density_right)
          if (denmax.gt.zero) then
           denjump_scale_temp=abs(density_left-density_right)/denmax

           if ((denjump_scale_temp.ge.zero).and. &
               (denjump_scale_temp.le.one)) then
            ! do nothing
           else
            print *,"denjump_scale_temp invalid"
            stop
           endif

           if (denjump_scale_temp.gt.denjump_scale) then
            denjump_scale=denjump_scale_temp
           endif

          else
           print *,"denmax invalid"
           stop
          endif

         endif
        else if (material_type(im_primaryR).eq.0) then
         ! do nothing
        else
         print *,"material_type(im_primaryR) invalid"
         stop
        endif 
       else if (is_rigid(im_primaryR).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid GODUNOV_3D.F90"
        stop
       endif

       cc_diag=max(cleft_diag,cright_diag)  ! c^2
       cc=max(cleft,cright)  ! c^2

        ! check_user_defined_velbc is declared in: PROB.F90
       call check_user_defined_velbc(time,dirnormal, &
             uu_estdt, & ! INTENT(inout)
             dx)

       uu_estdt_core = max(uu_estdt_core, &
           abs(uu_estdt)+abs(uu_estdt_phase_change))
       c_core = max(c_core,abs(cc_diag))  ! c^2

        !weymouth_factor=weymouth_cfl/cfl
       effective_velocity=abs(uu_estdt_phase_change/weymouth_factor)
       if (effective_velocity.gt.zero) then
         !u_eff=u*cfl/weymouth_cfl
         !dt=dx/u_eff=dx * weymouth_cfl / (u cfl)
         !cfl * dt = dx * weymouth_cfl/u
        dt_min=min(dt_min,hx/effective_velocity)
       else if (effective_velocity.eq.zero) then
        ! do nothing
       else
        print *,"effective_velocity invalid"
        stop
       endif

       effective_velocity=abs(uu_estdt/weymouth_factor)+sqrt(cc)
       if (effective_velocity.gt.zero) then
        dt_min=min(dt_min,hx/effective_velocity)
       else if (effective_velocity.eq.zero) then
        ! do nothing
       else
        print *,"effective_velocity invalid"
        stop
       endif

        ! fluid region
       if ((is_rigid(im_primaryL).eq.0).and. &
           (is_rigid(im_primaryR).eq.0)) then 

        call fluid_interface(LSleft,LSright,gradh, &
         im_opp,im, &
         im_left_interface,im_right_interface)

        if (gradh.ne.zero) then
         if ((im.gt.num_materials).or. &
             (im_opp.gt.num_materials).or. &
             (im.lt.1).or. &
             (im_opp.lt.1)) then
          print *,"im or im_opp bust 8"
          stop
         endif
         call get_iten(im,im_opp,iten)
         if ((iten.ge.1).and.(iten.le.num_interfaces)) then
          ! do nothing
         else
          print *,"iten invalid"
          stop
         endif
        else if (gradh.eq.zero) then
         im=0
         im_opp=0
         iten=0
        else
         print *,"gradh invalid"
         stop
        endif

       else if ((is_rigid(im_primaryL).eq.1).or. &
                (is_rigid(im_primaryR).eq.1)) then
        gradh=zero
        im=0
        im_opp=0
        iten=0
       else
        print *,"im_primaryL, or im_primaryR invalid"
        stop 
       endif

        ! gradh<>0 if the levelset function changes sign across a MAC
        ! face.
       if (gradh.ne.zero) then

        if (level_cap_wave_speed(iten).lt.zero) then
         print *,"level_cap wave speed not initialized fort_estdt"
         stop
        else if (level_cap_wave_speed(iten).eq.zero) then
         ! do nothing
        else if (level_cap_wave_speed(iten).gt.zero) then

         istenlo(3)=0
         istenhi(3)=0
         do dir2=1,SDIM
          istenlo(dir2)=ivec(dir2)-1
          istenhi(dir2)=ivec(dir2)+1
          if (dirnormal+1.eq.dir2) then
           istenlo(dir2)=ivec(dir2)-2
          endif
         enddo !dir2=1..sdim

         triple_flag=0
         do kalt=istenlo(3),istenhi(3)
         do jalt=istenlo(2),istenhi(2)
         do ialt=istenlo(1),istenhi(1)
          do im_sub=1,num_materials
           LSsub(im_sub)=dist(D_DECL(ialt,jalt,kalt),im_sub)
          enddo
          call get_primary_material(LSsub,im_sub)
          if ((im_sub.ne.im).and. &
              (im_sub.ne.im_opp)) then
           triple_flag=1
          endif
         enddo ! kalt
         enddo ! jalt
         enddo ! ialt

         dthold=hx/level_cap_wave_speed(iten)
         if (triple_flag.eq.1) then
          dthold=dthold/three
         endif

         dt_min=min(dt_min,dthold)
        else
         print *,"level_cap_wave_speed is NaN in fort_estdt"
         stop
        endif

       else if (gradh.eq.zero) then
        ! do nothing
       else
        print *,"gradh invalid fort_estdt"
        stop
       endif 

      enddo !i
      enddo !j
      enddo !k

      if (local_gravity_mag.gt.zero) then

       if (denjump_scale.gt.zero) then

        if (denjump_scale.le.one) then
         gravity_reference_wavelen=gravity_reference_wavelen_in
         call SUB_reference_wavelen(gravity_reference_wavelen)
         if (gravity_reference_wavelen.gt.zero) then
           ! gravity_wave_speed is declared in PROB.F90
          call gravity_wave_speed(gravity_reference_wavelen, &
            local_gravity_mag*denjump_scale,ugrav)
         else
          print *,"gravity_reference_wavelen invalid"
          stop
         endif
        else 
         print *,"require denjump_scale in [0,1]"
         stop
        endif

        if (ugrav.gt.zero) then
         dthold=dxmin/ugrav 
         dt_min=min(dt_min,dthold)
        else
         print *,"ugrav invalid 1 fort_estdt"
         print *,"uu_estdt ",uu_estdt
         print *,"denjump_scale ",denjump_scale
         print *,"local_gravity_mag ",local_gravity_mag
         print *,"dxmin ",dxmin
         stop
        endif
       else if (denjump_scale.eq.zero) then
        ! do nothing
       else
        print *,"denjump_scale cannot be negative, denjump = ",denjump_scale
        stop
       endif

      else if (local_gravity_mag.eq.zero) then
       ! do nothing
      else
       print *,"local_gravity_mag is NaN fort_estdt"
       stop
      endif 

       ! dirnormal=0..sdim-1
      if (uu_estdt_max(dirnormal+1).lt.uu_estdt_core) then
       uu_estdt_max(dirnormal+1)=uu_estdt_core
      endif
      if (uu_estdt_max(SDIM+1).lt.c_core) then
       uu_estdt_max(SDIM+1)=c_core  ! c^2
      endif

      return
      end subroutine fort_estdt


       ! mask=1 if cell not covered.
       ! masknbr=1 fine-fine border cells and interior cells.
       ! masknbr=0 coarse-fine cells and cells outside domain.
       ! called from getStateMOM_DEN
      subroutine fort_derive_mom_den( &
       im_parm, & ! 1..num_materials
       ngrow, &
       constant_density_all_time, & ! 1..num_materials
       spec_material_id_AMBIENT, &  ! 1..num_species_var
       presbc_arr, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       dt, &
       mask,DIMS(mask), &
       masknbr,DIMS(masknbr), &
       vol,DIMS(vol), &
       eosdata,DIMS(eosdata), &
       momden,DIMS(momden), &
       recon,DIMS(recon), &
       xlo,dx, &
       DrhoDT, &
       override_density, &
       level,finest_level) &
      bind(c,name='fort_derive_mom_den')

      use probf90_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: im_parm ! 1..num_materials
      integer, INTENT(in) :: ngrow
      integer, INTENT(in) :: constant_density_all_time(num_materials)
      integer, INTENT(in) :: spec_material_id_AMBIENT(num_species_var+1)
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact

      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: DIMDEC(vol)
      integer, INTENT(in) :: DIMDEC(eosdata)
      integer, INTENT(in) :: DIMDEC(momden)
      integer, INTENT(in) :: DIMDEC(recon)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: DIMDEC(masknbr)
     
      real(amrex_real), INTENT(in), target :: mask(DIMV(mask)) 
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: masknbr(DIMV(masknbr)) 
      real(amrex_real), pointer :: masknbr_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: vol(DIMV(vol)) 
      real(amrex_real), pointer :: vol_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: eosdata(DIMV(eosdata), &
               num_state_material*num_materials)
      real(amrex_real), pointer :: eosdata_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out), target :: &
              momden(DIMV(momden),num_materials)
      real(amrex_real), pointer :: momden_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
              recon(DIMV(recon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: recon_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: presbc_arr(SDIM,2)

      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: override_density(num_materials)
      real(amrex_real), INTENT(in) :: DrhoDT(num_materials)
     
      integer i,j,k
      integer dir

      integer dencomp
      real(amrex_real) xpos(SDIM)
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) temperature
      real(amrex_real) density_of_TZ
      real(amrex_real) rho_base
      integer vofcomp
      real(amrex_real) local_vfrac

      mask_ptr=>mask
      masknbr_ptr=>masknbr
      vol_ptr=>vol
      eosdata_ptr=>eosdata
      momden_ptr=>momden
      recon_ptr=>recon

      if (bfact.lt.1) then
       print *,"bfact invalid66"
       stop
      endif
      if (ngrow.ge.1) then
       ! do nothing
      else
       print *,"ngrow>=1 required"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid dencor"
       stop
      endif
      if ((im_parm.ge.1).and.(im_parm.le.num_materials)) then
       ! do nothing
      else
       print *,"fort_derive_mom_den: im_parm invalid, im_parm=",im_parm
       stop
      endif
      vofcomp=(im_parm-1)*ngeom_recon+1

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
       print *,"levelrz invalid dencor"
       stop
      endif

      call checkbound_array1(fablo,fabhi,vol_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,eosdata_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,momden_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,recon_ptr,ngrow,-1)
      call checkbound_array1(fablo,fabhi,mask_ptr,ngrow,-1)
      call checkbound_array1(fablo,fabhi,masknbr_ptr,ngrow,-1)

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir=1,SDIM
        xpos(dir)=xsten(0,dir)
       enddo

       if (mask(D_DECL(i,j,k)).eq.zero) then
          ! a default for coarse grid cells covered by finer levels.
        momden(D_DECL(i,j,k),im_parm)=fort_denconst(im_parm)
       else if (mask(D_DECL(i,j,k)).eq.one) then
     
        dencomp=(im_parm-1)*num_state_material+ENUM_DENVAR+1

        if (constant_density_all_time(im_parm).eq.1) then
         rho_base=fort_denconst(im_parm)
        else if (constant_density_all_time(im_parm).eq.0) then
         rho_base=eosdata(D_DECL(i,j,k),dencomp)
        else
         print *,"constant_density_all_time(im_parm) invalid"
         stop
        endif

        ! rho=rho(T)
        if (override_density(im_parm).eq.1) then

         if (fort_material_type(im_parm).eq.0) then
          ! do nothing
         else
          print *,"override_density==1 for incomp material only"
          stop
         endif

         local_vfrac=recon(D_DECL(i,j,k),vofcomp)

         if ((local_vfrac.ge.-EPS1).and. &
             (local_vfrac.le.VOFTOL)) then
          momden(D_DECL(i,j,k),im_parm)=rho_base
         else if ((local_vfrac.ge.VOFTOL).and. &
                  (local_vfrac.le.one+EPS1)) then

           ! den,T
          temperature=eosdata(D_DECL(i,j,k),dencomp+1)

          if (DrhoDT(im_parm).le.zero) then
           ! do nothing
          else
           print *,"DrhoDT invalid"
           stop
          endif

            ! default: fort_DrhoDT(im)*(T-T0)
            ! fort_DrhoDT units: 1/Temperature
          call SUB_UNITLESS_EXPANSION_FACTOR( &
            im_parm,temperature,fort_tempconst(im_parm),density_of_TZ)
          density_of_TZ=rho_base*(one+density_of_TZ)

          if ((temperature.ge.zero).and. &
              (fort_tempconst(im_parm).ge.zero).and. &
              (fort_denconst(im_parm).gt.zero).and. &
              (rho_base.gt.zero)) then 
           ! do nothing
          else
           print *,"invalid parameters to get the density"
           print *,"im_parm=",im_parm
           print *,"temperature=",temperature
           print *,"density_of_TZ=",density_of_TZ
           print *,"rho_base=",rho_base
           print *,"fort_tempconst(im_parm)=",fort_tempconst(im_parm)
           stop
          endif

          if (density_of_TZ.gt.zero) then
           ! do nothing
          else if (density_of_TZ.le.zero) then
           print *,"WARNING density_of_TZ.le.zero"
           print *,"im_parm=",im_parm
           print *,"temperature=",temperature
           print *,"density_of_TZ=",density_of_TZ
           print *,"rho_base=",rho_base
           print *,"fort_tempconst(im_parm)=",fort_tempconst(im_parm)
           print *,"fort_tempcutoffmax(im_parm)=",fort_tempcutoffmax(im_parm)
          
           temperature=fort_tempcutoffmax(im_parm)

           call SUB_UNITLESS_EXPANSION_FACTOR( &
             im_parm,temperature,fort_tempconst(im_parm),density_of_TZ)
           density_of_TZ=rho_base*(one+density_of_TZ)
  
           if (density_of_TZ.gt.zero) then
            ! do nothing
           else
            print *,"density_of_TZ.le.zero (STILL)"
            stop
           endif

          else
           print *,"density_of_TZ is NaN"
           stop
          endif

          momden(D_DECL(i,j,k),im_parm)=density_of_TZ

         else
          print *,"local_vfrac invalid: ",local_vfrac
          stop
         endif

        else if ((override_density(im_parm).eq.0).or. &
                 (override_density(im_parm).eq.2)) then
         momden(D_DECL(i,j,k),im_parm)=rho_base
        else
         print *,"override_density invalid"
         stop
        endif

       else
        print *,"mask invalid"
        stop
       endif

      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_derive_mom_den


       ! (vel/RR) if R-THETA
       ! vel=vel*dt , call adjust_du if RZ, override if passive advect.
      subroutine fort_velmac_override( &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       velbc, &
       dt,time, &
       passive_veltime, &
       vel_time, &
       dir_absolute_direct_split, &
       normdir, &
       umactemp,DIMS(umactemp), &
       umac_displace,DIMS(umac_displace), &
       xlo,dx, &
       mac_grow, &
       map_forward, &
       level, &
       finest_level, &
       SDC_outer_sweeps, &
       ns_time_order, &
       divu_outer_sweeps, &
       num_divu_outer_sweeps) &
      bind(c,name='fort_velmac_override')

      use probf90_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: SDC_outer_sweeps
      integer, INTENT(in) :: ns_time_order
      integer, INTENT(in) :: divu_outer_sweeps
      integer, INTENT(in) :: num_divu_outer_sweeps
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: dir_absolute_direct_split
      integer, INTENT(in) :: normdir
      integer, INTENT(in) :: mac_grow,map_forward
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dt,time,vel_time,passive_veltime
      integer, INTENT(in) :: DIMDEC(umactemp)
      integer, INTENT(in) :: DIMDEC(umac_displace)
     
      real(amrex_real), INTENT(in), target :: umactemp(DIMV(umactemp)) 
      real(amrex_real), pointer :: umactemp_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: umac_displace(DIMV(umac_displace)) 
      real(amrex_real), pointer :: umac_displace_ptr(D_DECL(:,:,:))
      integer, INTENT(in) :: velbc(SDIM,2)

      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
     
      integer i,j,k
      integer ii,jj,kk
      integer idx,side
      real(amrex_real) delta
      real(amrex_real) hx
      real(amrex_real) RR
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      integer localbc

      umactemp_ptr=>umactemp
      umac_displace_ptr=>umac_displace

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid65"
       stop
      endif
 
      if (mac_grow.ne.2) then
       print *,"mac_grow invalid mac_grow=",mac_grow
       stop
      endif
      if ((map_forward.ne.0).and.(map_forward.ne.1)) then
       print *,"map_forward invalid"
       stop
      endif
      if ((normdir.ge.0).and.(normdir.lt.SDIM)) then
       ! do nothing
      else
       print *,"normdir invalid"
       stop
      endif
      if ((dir_absolute_direct_split.ge.0).and. &
          (dir_absolute_direct_split.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir_absolute_direct_split invalid"
       stop
      endif
      if ((level.ge.0).and.(level.lt.finest_level)) then
       if (fabhi(normdir+1)-fablo(normdir+1)+1.lt.4) then
        print *,"blocking factor should be at least 4"
        print *,"level,finest_level ",level,finest_level
        stop
       endif
      else if (level.eq.finest_level) then
       if (fabhi(normdir+1)-fablo(normdir+1)+1.lt.2) then
        print *,"blocking factor should be at least 2"
        print *,"level,finest_level ",level,finest_level
        stop
       endif
      else
       print *,"level or finest_level invalid"
       stop
      endif

      if (level.gt.finest_level) then
       print *,"finest_level invalid velmac override"
       stop
      else if (level.lt.0) then
       print *,"level invalid velmac override"
       stop
      endif
      if ((SDC_outer_sweeps.ge.0).and. &
          (SDC_outer_sweeps.lt.ns_time_order)) then
       ! do nothing
      else
       print *,"SDC_outer_sweeps invalid"
       stop
      endif
      if ((divu_outer_sweeps.lt.0).or. &
          (divu_outer_sweeps.ge.num_divu_outer_sweeps)) then
       print *,"divu_outer_sweeps invalid fort_velmac_override"
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
       print *,"levelrz invalid velmac override"
       stop
      endif

      call checkbound_array1(fablo,fabhi,umactemp_ptr,mac_grow,normdir)
      call checkbound_array1(fablo,fabhi,umac_displace_ptr,mac_grow,normdir)

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (normdir.eq.0) then
       ii=1
      else if (normdir.eq.1) then
       jj=1
      else if ((normdir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif

        ! 1. multiply velocity by dt.
        ! 2. adjust velocity if RZ.
        ! 3. override velocity if it is a passive advection problem.
        ! 4. copy into mac_velocity

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,mac_grow,normdir)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

        if (normdir.eq.0) then
         idx=i
        else if (normdir.eq.1) then
         idx=j
        else if ((normdir.eq.2).and.(SDIM.eq.3)) then
         idx=k
        else
         print *,"normdir invalid"
         stop
        endif

        call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,normdir)
        hx=xstenMAC(1,normdir+1)-xstenMAC(-1,normdir+1)
        if (hx.gt.zero) then
         ! do nothing
        else
         print *,"xstenMAC bust"
         stop
        endif

        RR=one
        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
        else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
         if (normdir.eq.1) then
          RR=xstenMAC(0,1)
         endif
        else
         print *,"levelrz invalid velmac override"
         stop
        endif

        delta=umactemp(D_DECL(i,j,k))

        side=0

        if (idx.le.fablo(normdir+1)) then
         side=1
        else if (idx.ge.fabhi(normdir+1)+1) then
         side=2
        else if ((idx.gt.fablo(normdir+1)).and. &
                 (idx.lt.fabhi(normdir+1)+1)) then
         ! do nothing
        else
         print *,"idx invalid"
         stop
        endif

        if ((side.eq.1).or.(side.eq.2)) then
         localbc=velbc(normdir+1,side)
         if (localbc.eq.REFLECT_ODD) then
          delta=zero
         else if (localbc.eq.EXT_DIR) then
          call velbc_override(vel_time,normdir+1,side,normdir+1, &
           delta, &
           xstenMAC,nhalf,dx,bfact)
         else if (localbc.eq.INT_DIR) then
          ! do nothing
         else if (localbc.eq.REFLECT_EVEN) then
          ! do nothing
         else if (localbc.eq.FOEXTRAP) then
          ! do nothing
         else
          print *,"localbc invalid"
          stop
         endif  ! cases for localbc 
        else if (side.eq.0) then
         ! do nothing
        else
         print *,"side invalid"
         stop
        endif  

        delta=dt*delta/RR

        ! modifies "delta" if "passive_advect_flag=1" or RZ. 
        call departure_node_split( &
          xstenMAC,nhalf,dx,bfact, &
          delta,passive_veltime, &
          normdir,dt,map_forward)

        if (abs(delta).ge.(one-0.001)*hx) then
         print *,"in: velmac_override"
         print *,"MAC: displacement exceeds grid cell"
         print *,"parameters to check:"
         print *,"ns.cfl"
         print *,"ns.fixed_dt"
         print *,"ns.fixed_dt_init"
         print *,"ns.min_velocity_for_dt"
         print *,"ns.init_shrink"
         print *,"ns.change_max"
         print *,"initial velocity"
         print *,"umactemp ",umactemp(D_DECL(i,j,k))
         print *,"delta (u dt) = ",delta
         print *,"hx=    ",hx
         print *,"dt=    ",dt
         print *,"dir_absolute_direct_split= ",dir_absolute_direct_split
         print *,"normdir= ",normdir
         print *,"i,j,k ",i,j,k 
         print *,"level,finest_level ",level,finest_level
         print *,"SDC_outer_sweeps,ns_time_order ", &
            SDC_outer_sweeps,ns_time_order
         print *,"levelrz=",levelrz
         stop
        endif
      
         ! find displacements 
        umac_displace(D_DECL(i,j,k))=delta

      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_velmac_override

        ! recon:
        ! vof,ref centroid,order,slope,intercept  x num_materials
        !
        ! FACECOMP_ICEMASK and FACECOMP_ICEFACECUT components are 
        !   initialized to one
        !   in "fort_init_physics_vars"
        ! if num_materials=2, num_interfaces=1
        ! if num_materials=3, num_interfaces=3    12 13 23
        ! if num_materials=4, num_interfaces=6    12 13 14 23 24 34
      subroutine fort_initjumpterm( &
       mdotplus, &
       mdotminus, &
       mdotcount, &
       time, &
       level, &
       finest_level, &
       saturation_temp, &
       freezing_model, &
       distribute_from_target, &
       constant_volume_mdot, &
       constant_density_all_time, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx,dt, &
       maskcov,DIMS(maskcov), &
       JUMPFAB,DIMS(JUMPFAB), &
       mdot,DIMS(mdot), &
       LSnew,DIMS(LSnew), &
       recon,DIMS(recon) ) &
      bind(c,name='fort_initjumpterm')

      use probf90_module
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      real(amrex_real), INTENT(inout) :: mdotplus
      real(amrex_real), INTENT(inout) :: mdotminus
      real(amrex_real), INTENT(inout) :: mdotcount
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: level,finest_level
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      integer, INTENT(in) :: constant_volume_mdot(2*num_interfaces)
      integer, INTENT(in) :: constant_density_all_time(num_materials)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(JUMPFAB)
      integer, INTENT(in) :: DIMDEC(mdot)
      integer, INTENT(in) :: DIMDEC(LSnew)
      integer, INTENT(in) :: DIMDEC(recon)
      real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: &
              JUMPFAB(DIMV(JUMPFAB),2*num_interfaces)
      real(amrex_real), pointer :: JUMPFAB_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: mdot(DIMV(mdot))
      real(amrex_real), pointer :: mdot_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: LSnew(DIMV(LSnew),num_materials)
      real(amrex_real), pointer :: LSnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target ::  &
        recon(DIMV(recon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: recon_ptr(D_DECL(:,:,:),:)

      integer i,j,k
      integer im,im_opp,ireverse,iten
      integer iten_shift
      integer im_source,im_dest
      real(amrex_real) jump_strength

      real(amrex_real) LL
      real(amrex_real) divu_material
      integer local_mask
      real(amrex_real) F_solid_sum
      integer imlocal,vofcomp,im_primary
      real(amrex_real) LS_local(num_materials)

      mdot_ptr=>mdot

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in convert_material"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      do im=1,num_materials-1
       do im_opp=im+1,num_materials
        do ireverse=0,1
         call get_iten(im,im_opp,iten)
         iten_shift=iten+ireverse*num_interfaces
         if (is_valid_freezing_modelF(freezing_model(iten_shift)).eq.1) then
          ! do nothing 
         else
          print *,"freezing_model invalid init jump term"
          print *,"iten,ireverse,num_interfaces ",iten,ireverse,num_interfaces
          stop
         endif
         if ((distribute_from_target(iten_shift).lt.0).or. &
             (distribute_from_target(iten_shift).gt.1)) then
          print *,"distribute_from_target invalid init jump term"
          print *,"iten,ireverse,num_interfaces ",iten,ireverse,num_interfaces
          stop
         endif
         if (constant_volume_mdot(iten_shift).eq.0) then 
          ! do nothing
         else if (constant_volume_mdot(iten_shift).eq.1) then

          ! distribute -sum mdot to the source:

          if ((fort_material_type(im).eq.0).and. &
              (fort_material_type(im_opp).eq.0)) then
           !do nothing
          else
           print *,"expecting constant_volume_mdot==0"
           print *,"im,fort_material_type ",im,fort_material_type(im)
           print *,"im_opp,fort_material_type ", &
                   im_opp,fort_material_type(im_opp)
           stop
          endif

         else if (constant_volume_mdot(iten_shift).eq.-1) then

          ! distribute -sum mdot to the dest:

          if ((fort_material_type(im).eq.0).and. &
              (fort_material_type(im_opp).eq.0)) then
           !do nothing
          else
           print *,"expecting constant_volume_mdot==0"
           print *,"im,fort_material_type ",im,fort_material_type(im)
           print *,"im_opp,fort_material_type ", &
                   im_opp,fort_material_type(im_opp)
           stop
          endif

         else
          print *,"constant_volume_mdot(iten_shift) invalid"
          stop
         endif
        enddo ! ireverse
       enddo ! im_opp
      enddo ! im

      maskcov_ptr=>maskcov
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      JUMPFAB_ptr=>JUMPFAB
      call checkbound_array(fablo,fabhi,JUMPFAB_ptr,ngrow_distance,-1)
      call checkbound_array1(fablo,fabhi,mdot_ptr,0,-1)
      LSnew_ptr=>LSnew
      call checkbound_array(fablo,fabhi,LSnew_ptr,1,-1)
      recon_ptr=>recon
      call checkbound_array(fablo,fabhi,recon_ptr,1,-1)
 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        do im=1,num_materials-1
         do im_opp=im+1,num_materials
          do ireverse=0,1
           if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
            print *,"im or im_opp bust 10"
            stop
           endif

           call get_iten(im,im_opp,iten)
           iten_shift=iten+ireverse*num_interfaces

           LL=get_user_latent_heat(iten_shift,room_temperature,1)

           jump_strength=JUMPFAB(D_DECL(i,j,k),iten_shift)
  
           if ((is_rigid(im).eq.1).or. &
               (is_rigid(im_opp).eq.1)) then 
            ! do nothing
           else if (LL.ne.zero) then
            if (ireverse.eq.0) then
             im_source=im
             im_dest=im_opp
            else if (ireverse.eq.1) then
             im_source=im_opp
             im_dest=im
            else
             print *,"ireverse invalid"
             stop
            endif
      
            F_solid_sum=zero 
            do imlocal=1,num_materials
             LS_local(imlocal)=LSnew(D_DECL(i,j,k),imlocal)
             if (is_rigid(imlocal).eq.1) then
              vofcomp=(imlocal-1)*ngeom_recon+1
              F_solid_sum=F_solid_sum+recon(D_DECL(i,j,k),vofcomp)
             else if (is_rigid(imlocal).eq.0) then
              ! do nothing
             else
              print *,"is_rigid invalid"
              stop
             endif
            enddo ! imlocal=1..num_materials 

            call get_primary_material(LS_local,im_primary)

            if (F_solid_sum.ge.zero) then
             if (F_solid_sum.lt.half) then
              if (is_rigid(im_primary).eq.0) then
             
               ! jump_strength units: cm^3/s^2
               divu_material=jump_strength

               if (divu_material.gt.zero) then
                mdotplus=mdotplus+divu_material
                mdotcount=mdotcount+one
               else if (divu_material.lt.zero) then
                mdotminus=mdotminus+divu_material
                mdotcount=mdotcount+one
               else if (divu_material.eq.zero) then
                ! do nothing
               else
                print *,"divu_material bust"
                stop
               endif
            
               mdot(D_DECL(i,j,k))=mdot(D_DECL(i,j,k))+divu_material
              else if (is_rigid(im_primary).eq.1) then
               ! do nothing
              else
               print *,"is_rigid invalid"
               stop
              endif
             else if ((F_solid_sum.ge.half).and. &
                      (F_solid_sum.le.one+EPS1)) then
              ! do nothing
             else
              print *,"F_solid_sum invalid"
              stop
             endif 
            else
             print *,"F_solid_sum invalid"
             stop
            endif 

           else if (LL.eq.zero) then
            ! do nothing
           else
            print *,"LL bust"
            stop
           endif  
          enddo ! ireverse=0...1
         enddo ! im_opp=im+1...num_materials
        enddo ! im=1...num_materials-1

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_initjumpterm

        ! recon:
        ! vof,ref centroid,order,slope,intercept  x num_materials
        !
        ! xface(FACECOMP_ICEMASK+1)=one and 
        ! xface(FACECOMP_ICEFACECUT+1) = one in
        ! fluid regions (not ice and not FSI_RIGID).
        !
        ! local_face(FACECOMP_ICEMASK+1) and 
        ! local_face(FACECOMP_ICEFACECUT+1) 
        ! are initialized to one in "fort_init_physics_vars."
        !
        ! if num_materials=2, num_interfaces=1
        ! if num_materials=3, num_interfaces=3    12 13 23
        ! if num_materials=4, num_interfaces=6    12 13 14 23 24 34
        ! This routine is called from NavierStokes.cpp: 
        !  NavierStokes::level_init_icemask_and_icefacecut()
        !   which is called from
        !     NavierStokes::make_physics_varsALL
      subroutine fort_init_icemask_and_icefacecut( &
       nden, &
       time, &
       level,finest_level, &
       saturation_temp, &
       freezing_model, &
       distribute_from_target, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       dt, &
       maskcov,DIMS(maskcov), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       denstate,DIMS(denstate), &
       LSnew,DIMS(LSnew), &
       recon,DIMS(recon) ) &
      bind(c,name='fort_init_icemask_and_icefacecut')
      use probf90_module
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: nden
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: level,finest_level
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growloMAC(3),growhiMAC(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(xface)
      integer, INTENT(in) :: DIMDEC(yface)
      integer, INTENT(in) :: DIMDEC(zface)
      integer, INTENT(in) :: DIMDEC(denstate)
      integer, INTENT(in) :: DIMDEC(LSnew)
      integer, INTENT(in) :: DIMDEC(recon)
      real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: &
              xface(DIMV(xface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: xface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              yface(DIMV(yface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: yface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              zface(DIMV(zface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: zface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: denstate(DIMV(denstate),nden)
      real(amrex_real), pointer :: denstate_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: LSnew(DIMV(LSnew),num_materials)
      real(amrex_real), pointer :: LSnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              recon(DIMV(recon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: recon_ptr(D_DECL(:,:,:),:)

      integer complement_flag

      integer i,j,k
      integer dir
      integer dir2
      integer ii,jj,kk
      integer ireverse
      integer iten
      integer im,im_opp
      integer im_left,im_opp_left,im_primary_left
      integer im_right,im_opp_right,im_primary_right
      integer ireverse_left,ireverse_right
      real(amrex_real) denstateleft(nden)
      real(amrex_real) denstateright(nden)
      real(amrex_real) LSleft(num_materials)
      real(amrex_real) LSright(num_materials)
      real(amrex_real) VOFleft(num_materials)
      real(amrex_real) VOFright(num_materials)
      integer vofcomp

      real(amrex_real) ice_test,cut_test
      real(amrex_real) icemask_left
      real(amrex_real) icemask_right
      real(amrex_real) icefacecut_left
      real(amrex_real) icefacecut_right
      real(amrex_real) icemask
      real(amrex_real) icefacecut
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) xmac(SDIM)
      integer local_mask_right
      integer local_mask_left

      maskcov_ptr=>maskcov
      xface_ptr=>xface
      yface_ptr=>yface
      zface_ptr=>zface
      denstate_ptr=>denstate
      LSnew_ptr=>LSnew
      recon_ptr=>recon


      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid: ",time
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in fort_init_icemask_and_icefacecut"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nden.eq.num_materials*num_state_material) then
       ! do nothing
      else
       print *,"nden invalid"
       stop
      endif

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid: ",dt
       stop
      endif

      do im=1,num_materials-1
       do im_opp=im+1,num_materials
        do ireverse=0,1
         call get_iten(im,im_opp,iten)
         if (is_valid_freezing_modelF( &
              freezing_model(iten+ireverse*num_interfaces)).eq.1) then
          ! do nothing 
         else
          print *,"freezing_model invalid fort_init_icemask_and_icefacecut"
          print *,"iten,ireverse,num_interfaces ",iten,ireverse,num_interfaces
          stop
         endif
         if ((distribute_from_target(iten+ireverse*num_interfaces).lt.0).or. &
             (distribute_from_target(iten+ireverse*num_interfaces).gt.1)) then
          print *,"distribute_from_target err fort_init_icemask_and_icefacecut"
          print *,"iten,ireverse,num_interfaces ",iten,ireverse,num_interfaces
          stop
         endif
        enddo ! ireverse
       enddo ! im_opp
      enddo ! im

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,xface_ptr,0,0)
      call checkbound_array(fablo,fabhi,yface_ptr,0,1)
      call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1)
      call checkbound_array(fablo,fabhi,denstate_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LSnew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,recon_ptr,1,-1)
 
      do dir=0,SDIM-1
       ii=0
       jj=0
       kk=0
       if (dir.eq.0) then
        ii=1
       else if (dir.eq.1) then
        jj=1
       else if ((dir.eq.2).and.(SDIM.eq.3)) then
        kk=1
       else
        print *,"dir invalid fort_init_icemask_and_icefacecut"
        stop
       endif

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growloMAC,growhiMAC,0,dir) 
       do k=growloMAC(3),growhiMAC(3)
       do j=growloMAC(2),growhiMAC(2)
       do i=growloMAC(1),growhiMAC(1)

        local_mask_right=NINT(maskcov(D_DECL(i,j,k))) 
        local_mask_left=NINT(maskcov(D_DECL(i-ii,j-jj,k-kk))) 

         ! check if this is an uncovered face.
        if ((local_mask_right.eq.1).or. &
            (local_mask_left.eq.1)) then
         
         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir)
         do dir2=1,SDIM
          xmac(dir2)=xstenMAC(0,dir2)
         enddo 
         do im=1,num_materials
          LSleft(im)=LSnew(D_DECL(i-ii,j-jj,k-kk),im)
          LSright(im)=LSnew(D_DECL(i,j,k),im)
          vofcomp=(im-1)*ngeom_recon+1
          VOFleft(im)=recon(D_DECL(i-ii,j-jj,k-kk),vofcomp)
          VOFright(im)=recon(D_DECL(i,j,k),vofcomp)
         enddo
         do im=1,nden
          denstateleft(im)=denstate(D_DECL(i-ii,j-jj,k-kk),im)
          denstateright(im)=denstate(D_DECL(i,j,k),im)
         enddo

         complement_flag=0

          ! get_icemask_and_icefacecut defined in PROB.F90
          ! get_icemask_and_icefacecut is "triggered" for both ice 
          !  materials and "is_FSI_rigid" materials.
          ! this routine: fort_init_icemask_and_icefacecut
         call get_icemask_and_icefacecut( &
          nden, &
          xmac, &
          time, &
          dx,bfact, &
          icemask_left, &  ! 0 or 1
          icefacecut_left, & ! 0<=f<=1
          im_left, &
          im_opp_left, &
          im_primary_left, &
          ireverse_left, &
          denstateleft, &
          LSleft, &
          VOFleft, &
          distribute_from_target, &
          complement_flag)

          ! get_icemask_and_icefacecut defined in PROB.F90
          ! get_icemask_and_icefacecut is "triggered" for both ice 
          !  materials and "is_FSI_rigid" materials.
          ! this routine: fort_init_icemask_and_icefacecut
         call get_icemask_and_icefacecut( &
          nden, &
          xmac, &
          time, &
          dx,bfact, &
          icemask_right, &  ! 0 or 1
          icefacecut_right, & ! 0<=f<=1
          im_right, &
          im_opp_right, &
          im_primary_right, &
          ireverse_right, &
          denstateright, &
          LSright, &
          VOFright, &
          distribute_from_target, &
          complement_flag)

         icemask=min(icemask_left,icemask_right)

         if (icemask.eq.one) then
          ! do nothing
         else if (icemask.eq.zero) then
          ! do nothing
         else
          print *,"icemask invalid"
          stop
         endif

         if ((icefacecut_left.ge.zero).and. &
             (icefacecut_right.ge.zero).and. &
             (icefacecut_left.le.one).and. &
             (icefacecut_right.le.one)) then

          if ((icefacecut_left.eq.zero).and. &
              (icefacecut_right.eq.zero)) then
           icefacecut=zero
          else if ((icefacecut_left.eq.zero).or. &
                   (icefacecut_right.eq.zero)) then
           icefacecut=zero
          else if ((icefacecut_left.eq.one).and. &
                   (icefacecut_right.eq.one)) then
           icefacecut=one
          else if ((icefacecut_left.gt.zero).and. &
                   (icefacecut_left.le.one).and. &
                   (icefacecut_right.gt.zero).and. &
                   (icefacecut_right.le.one)) then
           icefacecut=min(icefacecut_left,icefacecut_right)
          else
           print *,"icefacecut_left or icefacecut_right invalid"
           print *,"icefacecut_left ",icefacecut_left
           print *,"icefacecut_right ",icefacecut_right
           stop
          endif

         else
          print *,"icefacecut_left or icefacecut_right invalid"
          print *,"icefacecut_left ",icefacecut_left
          print *,"icefacecut_right ",icefacecut_right
          stop
         endif

         if (icemask.eq.one) then
          if (icefacecut.eq.one) then
           ! do nothing
          else
           print *,"icefacecut invalid(1)"
           print *,"icefacecut ",icefacecut
           print *,"icemask ",icemask
           stop
          endif
         else if (icemask.eq.zero) then
          if ((icefacecut.ge.zero).and. &
              (icefacecut.lt.one)) then
           ! do nothing
          else
           print *,"icefacecut invalid(2)"
           print *,"icefacecut ",icefacecut
           print *,"icefacecut_left ",icefacecut_left
           print *,"icefacecut_right ",icefacecut_right
           print *,"icemask ",icemask
           print *,"icemask_left ",icemask_left
           print *,"icemask_right ",icemask_right
           stop
          endif
         else
          print *,"icemask invalid icemask=",icemask
          print *,"icefacecut=",icefacecut
          stop
         endif

         if ((icefacecut.ge.zero).and. &
             (icefacecut.le.one)) then
          ! do nothing
         else
          print *,"icefacecut invalid icefacecut=",icefacecut 
          stop
         endif

         if (icefacecut.ge.icemask) then
          ! do nothing
         else
          print *,"expecting icefacecut>=icemask"
          print *,"icefacecut=",icefacecut
          print *,"icemask=",icemask
          stop
         endif

         if (dir.eq.0) then
          ice_test=xface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)
          cut_test=xface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
         else if (dir.eq.1) then
          ice_test=yface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)
          cut_test=yface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
         else if ((dir.eq.2).and.(SDIM.eq.3)) then
          ice_test=zface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)
          cut_test=zface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
         else
          print *,"dir invalid fort_init_icemask_and_icefacecut 2"
          stop
         endif

         if (ice_test.eq.zero) then
          ! do nothing
         else if (ice_test.eq.one) then
          ! do nothing
         else
          print *,"ice_test invalid ice_test=",ice_test
          print *,"cut_test=",cut_test
          stop
         endif
         if ((cut_test.ge.zero).and. &
             (cut_test.le.one)) then
          ! do nothing
         else
          print *,"cut_test invalid" 
          print *,"cut_test ",cut_test
          stop
         endif
  
         if (dir.eq.0) then
          xface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)=icefacecut
          xface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)=icemask
         else if (dir.eq.1) then
          yface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)=icefacecut
          yface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)=icemask
         else if ((dir.eq.2).and.(SDIM.eq.3)) then
          zface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)=icefacecut
          zface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)=icemask
         else
          print *,"dir invalid fort_init_icemask_and_icefacecut 3"
          stop
         endif

        else if ((local_mask_right.eq.0).and. &
                 (local_mask_left.eq.0)) then
         ! do nothing
        else
         print *,"local_mask_right or local_mask_left invalid"
         stop
        endif

       enddo !i
       enddo !j
       enddo !k
      enddo ! dir=0..sdim-1

      return
      end subroutine fort_init_icemask_and_icefacecut

      subroutine check_for_closest_UMAC( &
       xtarget, &
       fablo,fabhi, &
       level,finest_level, &
       i1,j1,k1, &
       mac_cell_index, &
       normdir, & ! 0..sdim-1
       im_cpp, & ! 0..num_materials-1
       umac_ptr, & 
       LS_ptr, & 
       mindist, &
       umac_trial)

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: normdir
      integer, INTENT(in) :: im_cpp
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: mac_cell_index(SDIM)
      integer, INTENT(in) :: i1,j1,k1
      real(amrex_real), INTENT(inout) :: mindist
      real(amrex_real), INTENT(inout) :: umac_trial

      real(amrex_real), INTENT(in), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: umac_ptr(D_DECL(:,:,:))

      integer, PARAMETER :: nhalf=3
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) local_umac
      real(amrex_real) local_dist
      integer local_dir
      integer ii,jj,kk
      integer i,j,k
      real(amrex_real) LSLEFT(num_materials)
      real(amrex_real) LSRIGHT(num_materials)
      integer imL,imR
      integer imlocal

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid check_for_closest_UMAC"
       stop
      endif

      if ((normdir.ge.0).and.(normdir.lt.SDIM)) then
       ! do nothing
      else
       print *,"normdir invalid check_for_closest_UMAC"
       stop
      endif
      if (ngrow_distance.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance invalid"
       stop
      endif
      if ((im_cpp.ge.0).and.(im_cpp.lt.num_materials)) then
       ! do nothing
      else
       print *,"im_cpp invalid"
       stop
      endif
      if (is_rigid(im_cpp+1).eq.0) then
       ! do nothing
      else
       print *,"is_rigid(im_cpp+1) invalid"
       stop
      endif

      call checkbound_array1(fablo,fabhi,umac_ptr, &
        ngrow_distance,normdir)
      call checkbound_array(fablo,fabhi,LS_ptr,ngrow_distance,-1)

      ii=0
      jj=0
      kk=0
      if (normdir.eq.0) then
       ii=1
      else if (normdir.eq.1) then
       jj=1
      else if ((normdir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif

      i=mac_cell_index(1)+i1
      j=mac_cell_index(2)+j1
      k=0
      if (SDIM.eq.3) then
       k=mac_cell_index(SDIM)+k1
      endif

       ! normdir=0..sdim-1
      call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,normdir)

      do imlocal=1,num_materials
       call safe_data(i-ii,j-jj,k-kk,imlocal,LS_ptr,LSLEFT(imlocal))
       call safe_data(i,j,k,imlocal,LS_ptr,LSRIGHT(imlocal))
      enddo
      call get_primary_material(LSLEFT,imL)
      call get_primary_material(LSRIGHT,imR)

      if ((imL.ge.1).and.(imL.le.num_materials).and. &
          (imR.ge.1).and.(imR.le.num_materials)) then
       ! do nothing
      else
       print *,"imL or imR invalid"
       stop
      endif
   
      if ((LSLEFT(im_cpp+1).ge.zero).or. &
          (LSRIGHT(im_cpp+1).ge.zero)) then
       call safe_data_single(i,j,k,umac_ptr,local_umac)
       local_dist=zero
       do local_dir=1,SDIM
        local_dist=local_dist+(xstenMAC(0,local_dir)-xtarget(local_dir))**2
       enddo
       local_dist=sqrt(local_dist)
       if (mindist.eq.-one) then
        mindist=local_dist
        umac_trial=local_umac
       else if (mindist.ge.zero) then
        if (local_dist.lt.mindist) then
         mindist=local_dist
         umac_trial=local_umac
        else if (local_dist.ge.mindist) then
         ! do nothing
        else
         print *,"local_dist or mindist invalid"
         stop
        endif
       else
        print *,"mindist invalid"
        stop
       endif
      else if ((LSLEFT(im_cpp+1).lt.zero).and. &
               (LSRIGHT(im_cpp+1).lt.zero)) then
       ! do nothing
      else
       print *,"LSLEFT or LSRIGHT invalid"
       stop
      endif

      return
      end subroutine check_for_closest_UMAC


      subroutine fort_extend_mac_vel( &
       tid_current, &
       level, &
       finest_level, &
       normdir, & ! 0..sdim-1
       im_cpp, & ! 0..num_materials-1
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       time, &
       dt, &
       velbc, &
       mask, &  !mask=1 if not covered by level+1 or outside the domain
       DIMS(mask), &
       umac,DIMS(umac), & 
       umac_mask,DIMS(umac_mask), & 
       scalar_mask,DIMS(scalar_mask), & 
       divu_mask,DIMS(divu_mask), & 
       LS,DIMS(LS)) &
      bind(c,name='fort_extend_mac_vel')

      use probcommon_module
      use geometry_intersect_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: tid_current
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: normdir
      integer, INTENT(in) :: im_cpp
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: velbc(SDIM,2,SDIM)
      integer, INTENT(in) :: DIMDEC(mask) 
      integer, INTENT(in) :: DIMDEC(umac) 
      integer, INTENT(in) :: DIMDEC(umac_mask) 
      integer, INTENT(in) :: DIMDEC(scalar_mask) 
      integer, INTENT(in) :: DIMDEC(divu_mask) 
      integer, INTENT(in) :: DIMDEC(LS) 

      real(amrex_real), INTENT(in), target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: umac(DIMV(umac))
      real(amrex_real), pointer :: umac_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: umac_mask(DIMV(umac_mask))
      real(amrex_real), pointer :: umac_mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: scalar_mask(DIMV(scalar_mask))
      real(amrex_real), pointer :: scalar_mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: divu_mask(DIMV(divu_mask))
      real(amrex_real), pointer :: divu_mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: LS(DIMV(LS),num_materials*(SDIM+1))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) xclamped_minus_sten(-nhalf:nhalf,SDIM)
      real(amrex_real) xclamped_plus_sten(-nhalf:nhalf,SDIM)
      real(amrex_real) xclamped_minus(SDIM)
      real(amrex_real) xclamped_plus(SDIM)
      real(amrex_real) xtarget(SDIM)
      real(amrex_real) LSLEFT(num_materials)
      real(amrex_real) LSRIGHT(num_materials)

      integer bctest
      integer vel_boundary_fixed
      integer i,j,k
      integer ii,jj,kk
      integer i1,j1,k1
      integer k1low,k1high
      integer mask_left,mask_right
      integer imL,imR
      integer :: ivec(3)
      integer :: growlo(3),growhi(3)
      integer :: growlo_cell(3),growhi_cell(3)
      real(amrex_real) :: LS_clamped_minus,LS_clamped_plus
      real(amrex_real) :: temperature_clamped_minus,temperature_clamped_plus
      integer :: prescribed_flag
      real(amrex_real) :: vel_clamped_minus(SDIM),vel_clamped_plus(SDIM)
      integer :: need_closest_point
      real(amrex_real) :: nrmCP_LEFT(SDIM)
      real(amrex_real) :: nrmCP_RIGHT(SDIM)
      real(amrex_real) :: xI_LEFT(SDIM)
      real(amrex_real) :: xI_RIGHT(SDIM)
      integer :: mac_cell_index_LEFT(SDIM)
      integer :: mac_cell_index_RIGHT(SDIM)
      integer :: imlocal
      integer :: local_dir
      real(amrex_real) :: mindist
      real(amrex_real) :: dxmin
      real(amrex_real) :: EXTEND_BAND_WIDTH

      if ((tid_current.lt.0).or. &
          (tid_current.ge.geom_nthreads)) then
       print *,"tid_current invalid in fort_extend_mac_vel"
       stop
      endif
      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid fort_extend_mac_vel"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid46"
       stop
      endif

      if ((normdir.ge.0).and.(normdir.lt.SDIM)) then
       ! do nothing
      else
       print *,"normdir invalid fort_extend_mac_vel"
       stop
      endif
      if (ngrow_distance.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance invalid"
       stop
      endif
      if ((im_cpp.ge.0).and.(im_cpp.lt.num_materials)) then
       ! do nothing
      else
       print *,"im_cpp invalid"
       stop
      endif
      if (is_rigid(im_cpp+1).eq.0) then
       ! do nothing
      else
       print *,"is_rigid(im_cpp+1) invalid"
       stop
      endif

      call get_dxmin(dx,bfact,dxmin)
      EXTEND_BAND_WIDTH=two*dxmin

      mask_ptr=>mask
      umac_ptr=>umac
      umac_mask_ptr=>umac_mask
      scalar_mask_ptr=>scalar_mask
      divu_mask_ptr=>divu_mask
      LS_ptr=>LS

      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,umac_ptr, &
        ngrow_distance,normdir)
      call checkbound_array1(fablo,fabhi,umac_mask_ptr, &
        ngrow_distance,normdir)
      call checkbound_array1(fablo,fabhi,scalar_mask_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,divu_mask_ptr,0,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,ngrow_distance,-1)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0,normdir)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo_cell,growhi_cell,0) 

      ii=0
      jj=0
      kk=0
      if (normdir.eq.0) then
       ii=1
      else if (normdir.eq.1) then
       jj=1
      else if ((normdir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif

      k1low=-1
      k1high=1
      if (SDIM.eq.2) then
       k1low=0
       k1high=0
      endif

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       ivec(1)=i
       ivec(2)=j
       ivec(3)=k

       mask_left=NINT(mask(D_DECL(i-ii,j-jj,k-kk)))
       mask_right=NINT(mask(D_DECL(i,j,k)))

       if ((mask_left.eq.1).or.(mask_right.eq.1)) then

         ! normdir=0..sdim-1
        call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,normdir)
        call gridsten_level(xclamped_minus_sten,i-ii,j-jj,k-kk,level,nhalf)
        call gridsten_level(xclamped_plus_sten,i,j,k,level,nhalf)

        if (ivec(normdir+1).eq.fablo(normdir+1)) then
         bctest=velbc(normdir+1,1,normdir+1)
         if ((bctest.eq.REFLECT_ODD).or. &
             (bctest.eq.EXT_DIR)) then
          vel_boundary_fixed=1
         else if ((bctest.eq.INT_DIR).or. &
                  (bctest.eq.FOEXTRAP)) then
          ! do nothing
         else
          print *,"bctest invalid"
          stop
         endif
        else if (ivec(normdir+1).eq.fabhi(normdir+1)+1) then
         bctest=velbc(normdir+1,2,normdir+1)
         if ((bctest.eq.REFLECT_ODD).or. &
             (bctest.eq.EXT_DIR)) then
          vel_boundary_fixed=2
         else if ((bctest.eq.INT_DIR).or. &
                  (bctest.eq.FOEXTRAP)) then
          ! do nothing
         else
          print *,"bctest invalid"
          stop
         endif
        else if ((ivec(normdir+1).gt.fablo(normdir+1)).and. &
                 (ivec(normdir+1).lt.fabhi(normdir+1)+1)) then
         vel_boundary_fixed=0
        else
         print *,"ivec invalid"
         stop
        endif
      
        do local_dir=1,SDIM
         xclamped_minus(local_dir)=xclamped_minus_sten(0,local_dir)
         xclamped_plus(local_dir)=xclamped_plus_sten(0,local_dir)
         xtarget(local_dir)=xstenMAC(0,local_dir)
        enddo

        do imlocal=1,num_materials
         LSLEFT(imlocal)=LS(D_DECL(i-ii,j-jj,k-kk),imlocal)
         LSRIGHT(imlocal)=LS(D_DECL(i,j,k),imlocal)
        enddo
        call get_primary_material(LSLEFT,imL)
        call get_primary_material(LSRIGHT,imR)

        if ((imL.ge.1).and.(imL.le.num_materials).and. &
            (imR.ge.1).and.(imR.le.num_materials)) then
         ! do nothing
        else
         print *,"imL or imR invalid"
         stop
        endif
    
        if ((vel_boundary_fixed.eq.1).or. &
            (vel_boundary_fixed.eq.2)) then
         ! do nothing
        else if (vel_boundary_fixed.eq.0) then

         if ((imL.eq.im_cpp+1).or.(imR.eq.im_cpp+1)) then
          ! do nothing
         else if ((is_rigid(imL).eq.1).or. &
                  (is_rigid(imR).eq.1)) then
          ! do nothing
         else if ((is_rigid(imL).eq.0).and. &
                  (is_rigid(imR).eq.0)) then
           ! LS>0 if clamped
          call SUB_clamped_LS(xclamped_minus,time,LS_clamped_minus, &
           vel_clamped_minus,temperature_clamped_minus,prescribed_flag,dx)
          call SUB_clamped_LS(xclamped_plus,time,LS_clamped_plus, &
           vel_clamped_plus,temperature_clamped_plus,prescribed_flag,dx)

          if ((LS_clamped_minus.ge.zero).or. &
              (LS_clamped_plus.ge.zero)) then
           ! do nothing
          else if ((LS_clamped_minus.lt.zero).and. &
                   (LS_clamped_plus.lt.zero)) then

           need_closest_point=0

           if ((LSLEFT(im_cpp+1).ge.zero).or. &
               (LSRIGHT(im_cpp+1).ge.zero)) then
            ! do nothing
           else if ((LSLEFT(im_cpp+1).lt.-EXTEND_BAND_WIDTH).and. &
                    (LSRIGHT(im_cpp+1).lt.-EXTEND_BAND_WIDTH)) then
            if (ivec(normdir+1).le.growhi_cell(normdir+1)) then
             scalar_mask(D_DECL(i,j,k))=scalar_mask(D_DECL(i,j,k))+one
            endif
            if (ivec(normdir+1)-1.ge.growlo_cell(normdir+1)) then
             scalar_mask(D_DECL(i-ii,j-jj,k-kk))= &
               scalar_mask(D_DECL(i-ii,j-jj,k-kk))+one
            endif
           else if ((LSLEFT(im_cpp+1).lt.-EXTEND_BAND_WIDTH).and. &
                    (LSRIGHT(im_cpp+1).lt.zero).and. &
                    (LSRIGHT(im_cpp+1).ge.-EXTEND_BAND_WIDTH)) then
            if (ivec(normdir+1)-1.ge.growlo_cell(normdir+1)) then
             scalar_mask(D_DECL(i-ii,j-jj,k-kk))= &
                 scalar_mask(D_DECL(i-ii,j-jj,k-kk))+one
            endif
            need_closest_point=1
           else if ((LSRIGHT(im_cpp+1).lt.-EXTEND_BAND_WIDTH).and. &
                    (LSLEFT(im_cpp+1).lt.zero).and. &
                    (LSLEFT(im_cpp+1).ge.-EXTEND_BAND_WIDTH)) then
            if (ivec(normdir+1).le.growhi_cell(normdir+1)) then
             scalar_mask(D_DECL(i,j,k))=scalar_mask(D_DECL(i,j,k))+one
            endif
            need_closest_point=1
           else if ((LSLEFT(im_cpp+1).ge.-EXTEND_BAND_WIDTH).and. &
                    (LSRIGHT(im_cpp+1).ge.-EXTEND_BAND_WIDTH).and. &
                    (LSLEFT(im_cpp+1).lt.zero).and. &
                    (LSRIGHT(im_cpp+1).lt.zero)) then
            need_closest_point=1
           else
            print *,"LSLEFT or LSRIGHT invalid"
            stop
           endif

           if (need_closest_point.eq.1) then
            do local_dir=1,SDIM
             nrmCP_LEFT(local_dir)=LS(D_DECL(i-ii,j-jj,k-kk), &
                num_materials+im_cpp*SDIM+local_dir)
             xI_LEFT(local_dir)=xclamped_minus(local_dir)- &
                LSLEFT(im_cpp+1)*nrmCP_LEFT(local_dir)
             nrmCP_RIGHT(local_dir)=LS(D_DECL(i,j,k), &
                num_materials+im_cpp*SDIM+local_dir)
             xI_RIGHT(local_dir)=xclamped_plus(local_dir)- &
                LSRIGHT(im_cpp+1)*nrmCP_RIGHT(local_dir)
            enddo !local_dir=1..sdim
            call containing_MACcell(bfact,dx,xlo,fablo,xI_LEFT, &
             normdir,mac_cell_index_LEFT)
            call containing_MACcell(bfact,dx,xlo,fablo,xI_RIGHT, &
             normdir,mac_cell_index_RIGHT)

            mindist=-one
            do k1=k1low,k1high
            do j1=-1,1
            do i1=-1,1
             call check_for_closest_UMAC( &
               xtarget, &
               fablo,fabhi, &
               level,finest_level, &
               i1,j1,k1, &
               mac_cell_index_LEFT, &
               normdir, &
               im_cpp, &
               umac_ptr, &
               LS_ptr, &
               mindist, &
               umac(D_DECL(i,j,k)))
             call check_for_closest_UMAC( &
               xtarget, &
               fablo,fabhi, &
               level,finest_level, &
               i1,j1,k1, &
               mac_cell_index_RIGHT, &
               normdir, &
               im_cpp, &
               umac_ptr, &
               LS_ptr, &
               mindist, &
               umac(D_DECL(i,j,k)))
            enddo
            enddo
            enddo

            umac_mask(D_DECL(i,j,k))=one

           else if (need_closest_point.eq.0) then
            ! do nothing
           else
            print *,"need_closest_point invaid"
            stop
           endif

          else
           print *,"LS_clamped_minus or LS_clamped_plus invalid"
           stop
          endif
         else 
          print *,"(is_rigid(imL or imR) invalid"
          stop
         endif

        else
         print *,"vel_boundary_fixed invalid"
         stop
        endif

       else if ((mask_left.eq.0).and.(mask_right.eq.0)) then
        ! do nothing
       else
        print *,"mask_left or mask_right invalid"
        stop
       endif

      enddo !i 
      enddo !j
      enddo !k (face center "conserved" variables) 

      return
      end subroutine fort_extend_mac_vel


         ! 1=T11 2=T12 3=T22 4=T33 5=T13 6=T23
         ! rhoinverse is 1/den
      subroutine fort_tensorheat( &
       nstate, &
       xlo,dx,  &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       lsfab,DIMS(lsfab), &
       DeDTinverse, &
       DIMS(DeDTinverse), &
       vischeat,DIMS(vischeat), &
       tensor,DIMS(tensor), &
       gradu,DIMS(gradu), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       dt, &
       irz, &
       im_parm, &
       nden) &
      bind(c,name='fort_tensorheat')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: im_parm
      integer, INTENT(in) :: nden,nstate,level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(xface)
      integer, INTENT(in) :: DIMDEC(yface)
      integer, INTENT(in) :: DIMDEC(zface)
      integer, INTENT(in) :: DIMDEC(lsfab)
      integer, INTENT(in) :: DIMDEC(DeDTinverse)
      integer, INTENT(in) :: DIMDEC(vischeat)
      integer, INTENT(in) :: DIMDEC(tensor)
      integer, INTENT(in) :: DIMDEC(gradu)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in),target :: xface(DIMV(xface),FACECOMP_NCOMP)
      real(amrex_real), INTENT(in),target :: yface(DIMV(yface),FACECOMP_NCOMP)
      real(amrex_real), INTENT(in),target :: zface(DIMV(zface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: xface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: yface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: zface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
        lsfab(DIMV(lsfab),num_materials*(1+SDIM))
      real(amrex_real), pointer :: lsfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: DeDTinverse(DIMV(DeDTinverse))
      real(amrex_real), pointer :: DeDTinverse_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout),target :: vischeat(DIMV(vischeat))
      real(amrex_real), pointer :: vischeat_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: &
        tensor(DIMV(tensor),ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: tensor_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
         gradu(DIMV(gradu),AMREX_SPACEDIM_SQR)
      real(amrex_real), pointer :: gradu_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: irz
      integer :: i,j,k
      integer veldir,dir
      integer local_mask
      real(amrex_real) IEforce,Tforce
      real(amrex_real) one_over_DeDT
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) local_gradu(3,3)
      integer nbase
      integer grdcomp
      integer imlocal
      real(amrex_real) LScen(num_materials)
      real(amrex_real) Q(3,3)
      integer ii,jj

      if (bfact.lt.1) then
       print *,"bfact invalid50"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (irz.ne.levelrz) then
       print *,"irz invalid"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if ((im_parm.lt.0).or.(im_parm.ge.num_materials)) then
       print *,"im_parm invalid24"
       stop
      endif

      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      xface_ptr=>xface
      yface_ptr=>yface
      zface_ptr=>zface
      call checkbound_array(fablo,fabhi,xface_ptr,0,0)
      call checkbound_array(fablo,fabhi,yface_ptr,0,1)
      call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1)

      lsfab_ptr=>lsfab
      call checkbound_array(fablo,fabhi,lsfab_ptr,1,-1)
      DeDTinverse_ptr=>DeDTinverse
      call checkbound_array1(fablo,fabhi,DeDTinverse_ptr,0,-1)
      vischeat_ptr=>vischeat
      call checkbound_array1(fablo,fabhi,vischeat_ptr,0,-1)
      tensor_ptr=>tensor
      call checkbound_array(fablo,fabhi,tensor_ptr,0,-1)
      gradu_ptr=>gradu
      call checkbound_array(fablo,fabhi,gradu_ptr,0,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       do veldir=1,3
       do dir=1,3
        local_gradu(veldir,dir)=zero
       enddo 
       enddo 
       do dir=1,SDIM
        if (dir.eq.1) then
         nbase=TENSOR_TRANSPOSE_UX-1
        else if (dir.eq.2) then
         nbase=TENSOR_TRANSPOSE_UY-1
        else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
         nbase=TENSOR_TRANSPOSE_UZ-1
        else
         print *,"dir invalid viscoelastic heating"
         stop
        endif
        do veldir=1,SDIM
         grdcomp=nbase+veldir
         local_gradu(veldir,dir)=gradu(D_DECL(i,j,k),grdcomp)
        enddo ! veldir
       enddo ! dir

       if (irz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (irz.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (xsten(0,1).le.zero) then
         print *,"no neg domain in r-z"
         stop
        endif
       else if (irz.eq.COORDSYS_CYLINDRICAL) then
        if (xsten(-2,1).le.zero) then
         print *,"no neg domain in r-T"
         stop
        endif
       else
        print *,"irz invalid"
        stop
       endif

       do imlocal=1,num_materials
        LScen(imlocal)=lsfab(D_DECL(i,j,k),imlocal)
       enddo
       call get_primary_material(LScen,local_mask)

       if ((local_mask.eq.im_parm+1).and. &
           (LScen(im_parm+1).gt.zero)) then
        local_mask=1
       else if ((local_mask.ge.1).and.(local_mask.le.num_materials)) then
        local_mask=0
       else
        print *,"local_mask invalid"
        stop
       endif

       do ii=1,3
       do jj=1,3
        Q(ii,jj)=zero
       enddo
       enddo
       do dir=1,ENUM_NUM_TENSOR_TYPE
        call stress_index(dir,ii,jj)
        Q(ii,jj)=tensor(D_DECL(i,j,k),dir)
       enddo
       Q(2,1)=Q(1,2)
       Q(3,1)=Q(1,3)
       Q(3,2)=Q(2,3)

        ! E: div( u dot tau )=(u tau_11)_x+(u tau_12)_y+(u tau_13)_z+...
        ! e: grad u : tau

       if (local_mask.eq.0) then
        IEforce=zero
       else if (local_mask.eq.1) then

        IEforce=zero
        do veldir=1,SDIM
        do dir=1,SDIM
         IEforce=IEforce+local_gradu(veldir,dir)*Q(veldir,dir)
        enddo
        enddo

        one_over_DeDT=DeDTinverse(D_DECL(i,j,k))  ! 1/(rho cv)

        if (one_over_DeDT.gt.zero) then
         ! do nothing
        else
         print *,"one_over_DeDT invalid"
         stop
        endif

        Tforce=IEforce*dt*one_over_DeDT
    
        vischeat(D_DECL(i,j,k))=vischeat(D_DECL(i,j,k))+Tforce

       else
        print *,"local_mask invalid"
        stop
       endif

      enddo !i
      enddo !j
      enddo !k
 
      return
      end subroutine fort_tensorheat


         ! rhoinverse is 1/den
      subroutine fort_visctensorheat( &
       nsolve, &
       nstate, &
       xlo,dx,  &
       lsfab,DIMS(lsfab), &
       DeDTinverse, &
       DIMS(DeDTinverse), &
       vischeat,DIMS(vischeat), &
       xstress,DIMS(xstress), &
       ystress,DIMS(ystress), &
       zstress,DIMS(zstress), &
       gradu,DIMS(gradu), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       dt,irz, &
       nden) &
      bind(c,name='fort_visctensorheat')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: nsolve
      integer, INTENT(in) :: nden,nstate,level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(lsfab)
      integer, INTENT(in) :: DIMDEC(DeDTinverse)
      integer, INTENT(in) :: DIMDEC(vischeat)
      integer, INTENT(in) :: DIMDEC(xstress)
      integer, INTENT(in) :: DIMDEC(ystress)
      integer, INTENT(in) :: DIMDEC(zstress)
      integer, INTENT(in) :: DIMDEC(gradu)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in),target :: lsfab(DIMV(lsfab),num_materials)
      real(amrex_real), pointer :: lsfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: DeDTinverse(DIMV(DeDTinverse))
      real(amrex_real), pointer :: DeDTinverse_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout),target :: vischeat(DIMV(vischeat))
      real(amrex_real), pointer :: vischeat_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: xstress(DIMV(xstress),nsolve)
      real(amrex_real), INTENT(in),target :: ystress(DIMV(ystress),nsolve)
      real(amrex_real), INTENT(in),target :: zstress(DIMV(zstress),nsolve)
      real(amrex_real), pointer :: xstress_ptr(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: ystress_ptr(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: zstress_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in),target :: gradu(DIMV(gradu),AMREX_SPACEDIM_SQR)
      real(amrex_real), pointer :: gradu_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: irz

      integer :: i,j,k
      integer veldir,dir
      real(amrex_real) IEforce,Tforce
      real(amrex_real) one_over_DeDT
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) local_gradu(3,3)
      real(amrex_real) tensor(SDIM,SDIM)
      integer nbase
      integer grdcomp

      if (bfact.lt.1) then
       print *,"bfact invalid51"
       stop
      endif
      if (nsolve.ne.SDIM) then
       print *,"nsolve invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (irz.ne.levelrz) then
       print *,"irz invalid"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif

      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      lsfab_ptr=>lsfab
      call checkbound_array(fablo,fabhi,lsfab_ptr,1,-1)
      DeDTinverse_ptr=>DeDTinverse
      call checkbound_array1(fablo,fabhi,DeDTinverse_ptr,0,-1)
      vischeat_ptr=>vischeat
      call checkbound_array1(fablo,fabhi,vischeat_ptr,0,-1)
      xstress_ptr=>xstress
      ystress_ptr=>ystress
      zstress_ptr=>zstress
      call checkbound_array(fablo,fabhi,xstress_ptr,0,0)
      call checkbound_array(fablo,fabhi,ystress_ptr,0,1)
      call checkbound_array(fablo,fabhi,zstress_ptr,0,SDIM-1)
      gradu_ptr=>gradu
      call checkbound_array(fablo,fabhi,gradu_ptr,0,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       if (irz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (irz.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (xsten(0,1).gt.zero) then
         ! do nothing
        else
         print *,"no neg domain in r-z"
         stop
        endif
       else if (irz.eq.COORDSYS_CYLINDRICAL) then
        if (xsten(-2,1).gt.zero) then
         ! do nothing
        else
         print *,"no neg domain in r-T"
         stop
        endif
       else
        print *,"irz invalid"
        stop
       endif

        ! E: div( u dot tau )=(u tau_11)_x+(u tau_12)_y+(u tau_13)_z+...
        ! e: grad u : tau

       do veldir=1,3
       do dir=1,3
        local_gradu(veldir,dir)=zero
       enddo 
       enddo 
       do dir=1,SDIM
        if (dir.eq.1) then
         nbase=TENSOR_TRANSPOSE_UX-1
        else if (dir.eq.2) then
         nbase=TENSOR_TRANSPOSE_UY-1
        else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
         nbase=TENSOR_TRANSPOSE_UZ-1
        else
         print *,"dir invalid viscoelastic heating"
         stop
        endif
        do veldir=1,SDIM
         grdcomp=nbase+veldir
         local_gradu(veldir,dir)=gradu(D_DECL(i,j,k),grdcomp)
        enddo ! veldir
       enddo ! dir=1..sdim

       do veldir=1,SDIM
        tensor(veldir,1)=half*(xstress(D_DECL(i,j,k),veldir)+ &
               xstress(D_DECL(i+1,j,k),veldir))
        tensor(veldir,2)=half*(ystress(D_DECL(i,j,k),veldir)+ &
               ystress(D_DECL(i,j+1,k),veldir))
        if (SDIM.eq.3) then
         tensor(veldir,SDIM)=half*(zstress(D_DECL(i,j,k),veldir)+ &
                zstress(D_DECL(i,j,k+1),veldir))
        endif
       enddo ! veldir

       IEforce=zero
       do veldir=1,SDIM
       do dir=1,SDIM
        IEforce=IEforce+local_gradu(veldir,dir)*tensor(veldir,dir)
       enddo
       enddo

       one_over_DeDT=DeDTinverse(D_DECL(i,j,k))  ! 1/(rho cv)

       if (one_over_DeDT.gt.zero) then
        ! do nothing
       else
        print *,"one_over_DeDT invalid"
        stop
       endif
       Tforce=-IEforce*dt*one_over_DeDT
       vischeat(D_DECL(i,j,k))= &
        vischeat(D_DECL(i,j,k))+Tforce

      enddo !i
      enddo !j
      enddo !k
 
      return
      end subroutine fort_visctensorheat

      ! rhoinverse is 1/den
      ! curv((iten-1)*CURVCOMP_NCOMP+CURVCOMP_MARANGONI+dir)=
      !  mgoni_force(dir)=
      !  (I-nn^T)(grad sigma)delta 
      subroutine fort_marangoniforce( &
       nstate, &
       num_curv, &
       xlo,dx,  &
       ls,DIMS(ls), &
       rhoinverse, &
       DIMS(rhoinverse), &
       curv,DIMS(curv), &
       velnew,DIMS(velnew), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       bfact_grid, &
       level, &
       finest_level, &
       dt, &
       cur_time) &
      bind(c,name='fort_marangoniforce')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: num_curv
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(ls)
      integer, INTENT(in) :: DIMDEC(rhoinverse)
      integer, INTENT(in) :: DIMDEC(curv)
      integer, INTENT(in) :: DIMDEC(velnew)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: bfact_grid
      real(amrex_real), INTENT(in), target :: ls(DIMV(ls),num_materials*(SDIM+1))
      real(amrex_real), pointer :: ls_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: rhoinverse(DIMV(rhoinverse))
      real(amrex_real), pointer :: rhoinverse_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: curv(DIMV(curv),num_curv)
      real(amrex_real), pointer :: curv_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
        velnew(DIMV(velnew),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: velnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: dt,cur_time

      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) LScen(num_materials)
      integer im
      integer iten
      integer iforce
      integer dirloc
      integer i,j,k
      real(amrex_real) surface_tension_force(SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid52"
       stop
      endif
      if ((level.ge.0).and.(level.lt.finest_level)) then
       if (bfact_grid.lt.4) then
        print *,"bfact_grid invalid in fort_marangoniforce: ",bfact_grid
        stop
       endif
      else if (level.eq.finest_level) then
       if (bfact_grid.lt.2) then
        print *,"bfact_grid invalid in fort_marangoniforce: ",bfact_grid
        stop
       endif
      else
       print *,"level invalid in fort_marangoniforce"
       stop
      endif

      do dirloc=1,SDIM
       if ((fablo(dirloc)/bfact_grid)*bfact_grid.ne.fablo(dirloc)) then
        print *,"fablo mod bfact_grid not 0"
        stop
       endif
       if (((fabhi(dirloc)+1)/bfact_grid)*bfact_grid.ne.fabhi(dirloc)+1) then
        print *,"fabhi+1 mod bfact_grid not 0"
        stop
       endif
      enddo ! dirloc=1..sdim

      if (num_curv.ne.num_interfaces*CURVCOMP_NCOMP) then
       print *,"num_curv invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid marangoni force"
       stop
      endif
      if (finest_level.ne.fort_finest_level) then
       print *,"finest_level invalid marangoni force"
       stop
      endif

      ls_ptr=>ls
      call checkbound_array(fablo,fabhi,ls_ptr,2,-1)
      rhoinverse_ptr=>rhoinverse
      call checkbound_array1(fablo,fabhi,rhoinverse_ptr,1,-1)
      curv_ptr=>curv
      call checkbound_array(fablo,fabhi,curv_ptr,1,-1)
      velnew_ptr=>velnew
      call checkbound_array(fablo,fabhi,velnew_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       do dirloc=1,SDIM
        surface_tension_force(dirloc)=zero
       enddo

       do im=1,num_materials
        LScen(im)=ls(D_DECL(i,j,k),im)
       enddo
       call get_primary_material(LScen,im)

       if (is_rigid(im).eq.0) then ! fluid region

        do iten=1,num_interfaces

         do dirloc=1,SDIM
          iforce=(iten-1)*CURVCOMP_NCOMP+CURVCOMP_MARANGONI+dirloc
          surface_tension_force(dirloc)= &
           surface_tension_force(dirloc)+ &
           curv(D_DECL(i,j,k),iforce)*dt*rhoinverse(D_DECL(i,j,k))
         enddo  ! dirloc

        enddo !iten=1,num_interfaces
  
       else if (is_rigid(im).eq.1) then 
        ! do nothing
       else
        print *,"is_rigid(im) invalid"
        stop
       endif 
      
       do dirloc=1,SDIM
        velnew(D_DECL(i,j,k),dirloc)= &
         velnew(D_DECL(i,j,k),dirloc)+surface_tension_force(dirloc)
       enddo ! dirloc

      enddo!i
      enddo!j
      enddo!k

      return
      end subroutine fort_marangoniforce

       ! MEHDI VAHAB HEAT SOURCE
       ! T^new=T^* + dt * Q/(rho cv)
       ! Q units: J/(m^3 s)
       ! called from: make_heat_source
       ! make_heat_source is called from veldiffuseALL
      subroutine fort_heatsource( &
       nstate, &
       nden, &
       xlo,dx,  &
       temperature_source, &
       temperature_source_cen, &
       temperature_source_rad, &
       DeDTinverse, &
       DIMS(DeDTinverse), &
       Tnew,DIMS(Tnew), &
       lsfab,DIMS(lsfab), &
       recon,DIMS(recon), &
       vol,DIMS(vol), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       dt,time) &
      bind(c,name='fort_heatsource')
      use probf90_module
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: nden
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(DeDTinverse)
      integer, INTENT(in) :: DIMDEC(Tnew)
      integer, INTENT(in) :: DIMDEC(lsfab)
      integer, INTENT(in) :: DIMDEC(recon)
      integer, INTENT(in) :: DIMDEC(vol)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: temperature_source
      real(amrex_real), INTENT(in) :: temperature_source_cen(SDIM)
      real(amrex_real), INTENT(in) :: temperature_source_rad(SDIM)
      real(amrex_real), INTENT(in),target :: DeDTinverse(DIMV(DeDTinverse)) ! 1/(rho cv)
      real(amrex_real), pointer :: DeDTinverse_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout),target :: Tnew(DIMV(Tnew),nden)
      real(amrex_real), pointer :: Tnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: lsfab(DIMV(lsfab),num_materials)
      real(amrex_real), pointer :: lsfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: recon(DIMV(recon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: recon_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: vol(DIMV(vol))
      real(amrex_real), pointer :: vol_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in) :: dt,time
                         
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_cell(SDIM)
      real(amrex_real) LS(num_materials)
      real(amrex_real) VFRAC(num_materials)
      real(amrex_real) heat_source_local(num_materials)
      real(amrex_real) T_local(num_materials)
      real(amrex_real) den_local(num_materials)
      integer i,j,k
      integer im
      integer vofcomp,dencomp
      integer dirloc
      integer imattype
      real(amrex_real) heat_source_total,vfrac_total
      real(amrex_real) DeDT_local(num_materials)
      integer ispec
      real(amrex_real) massfrac_parm(num_species_var+1)

      if (bfact.lt.1) then
       print *,"bfact invalid53"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid heat source"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif
      if (temperature_source.ge.zero) then
       ! do nothing
      else
       print *,"temperature_source invalid"
       stop
      endif

      DeDTinverse_ptr=>DeDTinverse
      call checkbound_array1(fablo,fabhi,DeDTinverse_ptr,1,-1)
      Tnew_ptr=>Tnew
      call checkbound_array(fablo,fabhi,Tnew_ptr,1,-1)
      lsfab_ptr=>lsfab
      call checkbound_array(fablo,fabhi,lsfab_ptr,1,-1)
      recon_ptr=>recon
      call checkbound_array(fablo,fabhi,recon_ptr,1,-1)
      vol_ptr=>vol
      call checkbound_array1(fablo,fabhi,vol_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       do im=1,num_materials
        LS(im)=lsfab(D_DECL(i,j,k),im)
        vofcomp=(im-1)*ngeom_recon+1
        VFRAC(im)=recon(D_DECL(i,j,k),vofcomp) 
        dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
        den_local(im)=Tnew(D_DECL(i,j,k),dencomp)
        T_local(im)=Tnew(D_DECL(i,j,k),dencomp+1)
        if ((VFRAC(im).ge.-EPS1).and. &
            (VFRAC(im).le.one+EPS1)) then
         ! do nothing
        else
         print *,"VFRAC invalid: ",im,VFRAC(im)
         stop
        endif 
        if (den_local(im).gt.zero) then
         ! do nothing
        else
         print *,"den_local must be positive"
         stop
        endif
        if (T_local(im).gt.zero) then
         ! do nothing
        else
         print *,"T_local must be positive"
         stop
        endif
        call init_massfrac_parm(den_local(im),massfrac_parm,im)
        do ispec=1,num_species_var
         massfrac_parm(ispec)=Tnew(D_DECL(i,j,k),dencomp+1+ispec)
         if (massfrac_parm(ispec).ge.zero) then
          ! do nothing
         else
          print *,"massfrac_parm(ispec) invalid"
          stop
         endif
        enddo

        imattype=fort_material_type(im) 
          ! in otherwords, find c_v for material im
        call DeDT_material(den_local(im),massfrac_parm, &
         T_local(im), &
         DeDT_local(im),imattype,im)
       enddo ! im=1..num_materials

       do dirloc=1,SDIM
        xsten_cell(dirloc)=xsten(0,dirloc)
       enddo

       ! MEHDI VAHAB HEAT SOURCE
       ! for right flux boundary condition:  
       !   T_new = T_old + dt (-q_right + qleft)/dx
       ! T^new=T^* + dt * Q/(rho cv)
       ! Q units: J/(m^3 s)
       ! get_local_heat_source in PROB.F90
       call get_local_heat_source( &
         time,dt, &
         xsten_cell, &
         xsten, &  ! xsten(-nhalf:nhalf,SDIM) xsten(0,dir)=cell center
                   ! xsten(1,dir)=xcell + dx/2   xsten(-1,dir)=xcell-dx/2
                   ! xsten(2,dir)=xcell + dx     xsten(-2,dir)=xcell-dx
         nhalf, &
         temperature_source, &
         temperature_source_cen, &
         temperature_source_rad, &
         LS,VFRAC,T_local,den_local, &
         DeDT_local, &
         heat_source_local)

       heat_source_total=zero
       vfrac_total=zero 
       do im=1,num_materials
        heat_source_total=heat_source_total+VFRAC(im)*heat_source_local(im)
        vfrac_total=vfrac_total+VFRAC(im)
       enddo

       if (vfrac_total.gt.zero) then
        ! do nothing
       else
        print *,"vfrac_total invalid"
        stop
       endif

       heat_source_total=heat_source_total/vfrac_total
 
         ! DeDTinverse = 1/(rho cv)
       do im=1,num_materials
        T_local(im)=T_local(im)+ &
          dt*DeDTinverse(D_DECL(i,j,k))*heat_source_total
        if (1.eq.0) then
         if (heat_source_local(im).ne.zero) then
          print *,"x,im,heat_source_local ",xsten(0,1),xsten(0,2), &
           im,heat_source_local(im)
         endif
        endif

        dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
        Tnew(D_DECL(i,j,k),dencomp+1)=T_local(im)
       enddo ! im=1..num_materials

      enddo!i
      enddo!j
      enddo!k
 
      return
      end subroutine fort_heatsource


         ! rhoinverse is 1/den
      subroutine fort_semdeltaforce( &
       nstate, &
       project_option, &
       xlo,dx,  &
       deltafab,DIMS(deltafab), &
       maskSEM,DIMS(maskSEM), &
       rhoinverse, &
       DIMS(rhoinverse), &
       DeDTinverse, &
       DIMS(DeDTinverse), &
       velnew,DIMS(velnew), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,level, &
       dt) &
      bind(c,name='fort_semdeltaforce')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE


      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer :: i,j,k
      integer :: veldir
      integer, INTENT(in) :: DIMDEC(deltafab)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: DIMDEC(rhoinverse)
      integer, INTENT(in) :: DIMDEC(DeDTinverse)
      integer, INTENT(in) :: DIMDEC(velnew)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in),target :: deltafab(DIMV(deltafab),NSTATE_SDC)
      real(amrex_real), pointer :: deltafab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: rhoinverse(DIMV(rhoinverse))
      real(amrex_real), pointer :: rhoinverse_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: DeDTinverse(DIMV(DeDTinverse))
      real(amrex_real), pointer :: DeDTinverse_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout),target :: velnew(DIMV(velnew),nstate)
      real(amrex_real), pointer :: velnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: dt
      integer idst,isrc,im
      integer local_maskSEM


      if (bfact.lt.1) then
       print *,"bfact invalid54"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
      if ((project_option.eq.SOLVETYPE_HEAT).or. &
          (project_option.eq.SOLVETYPE_VISC)) then
       ! do nothing
      else
       print *,"project_option invalid; fort_semdeltaforce"
       stop
      endif

      deltafab_ptr=>deltafab
      call checkbound_array(fablo,fabhi,deltafab_ptr,0,-1)
      rhoinverse_ptr=>rhoinverse
      call checkbound_array1(fablo,fabhi,rhoinverse_ptr,1,-1)
      DeDTinverse_ptr=>DeDTinverse
      call checkbound_array1(fablo,fabhi,DeDTinverse_ptr,1,-1)
      velnew_ptr=>velnew
      call checkbound_array(fablo,fabhi,velnew_ptr,1,-1)
      maskSEM_ptr=>maskSEM
      call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
       if ((local_maskSEM.ge.1).and. &
           (local_maskSEM.le.num_materials)) then

        if (project_option.eq.SOLVETYPE_VISC) then ! viscosity

         do veldir=1,SDIM
          velnew(D_DECL(i,j,k),STATECOMP_VEL+veldir)= &
            velnew(D_DECL(i,j,k),STATECOMP_VEL+veldir)- &
            rhoinverse(D_DECL(i,j,k))* &
            deltafab(D_DECL(i,j,k),veldir)
         enddo ! veldir=1..sdim

        else if (project_option.eq.SOLVETYPE_HEAT) then ! thermal conduction

         isrc=1

         do im=1,num_materials 
          idst=STATECOMP_STATES+ &
            (im-1)*num_state_material+ENUM_TEMPERATUREVAR+1
          velnew(D_DECL(i,j,k),idst)= &
           velnew(D_DECL(i,j,k),idst)- &
           DeDTinverse(D_DECL(i,j,k))*deltafab(D_DECL(i,j,k),isrc)
         enddo ! im=1..num_materials

        else
         print *,"project_option invalid; fort_semdeltaforce"
         stop
        endif
       else if (local_maskSEM.eq.0) then
        ! do nothing
       else
        print *,"local_maskSEM invalid"
        stop
       endif

      enddo!i
      enddo!j
      enddo!k

      return
      end subroutine fort_semdeltaforce

      subroutine fort_semdeltaforce_face( &
       dir, &
       xlo,dx,  &
       deltafab,DIMS(deltafab), &
       maskSEM,DIMS(maskSEM), &
       xface,DIMS(xface), &
       xmac,DIMS(xmac), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,level, &
       dt) &
      bind(c,name='fort_semdeltaforce_face')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: dir
      integer, INTENT(in) :: level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer :: i,j,k
      integer :: ii,jj,kk
      integer, INTENT(in) :: DIMDEC(deltafab)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: DIMDEC(xface)
      integer, INTENT(in) :: DIMDEC(xmac)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in),target :: deltafab(DIMV(deltafab))
      real(amrex_real), pointer :: deltafab_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: xface(DIMV(xface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: xface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: xmac(DIMV(xmac))
      real(amrex_real), pointer :: xmac_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in) :: dt
      integer maskleft,maskright

      if (bfact.lt.1) then
       print *,"bfact invalid55"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid sem delta force face"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid sem delta force face 2"
       stop
      endif

      deltafab_ptr=>deltafab
      call checkbound_array1(fablo,fabhi,deltafab_ptr,0,dir)
      xface_ptr=>xface
      call checkbound_array(fablo,fabhi,xface_ptr,0,dir)
      xmac_ptr=>xmac
      call checkbound_array1(fablo,fabhi,xmac_ptr,0,dir)
      maskSEM_ptr=>maskSEM
      call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0,dir)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       maskleft=NINT(maskSEM(D_DECL(i-ii,j-jj,k-kk)))
       maskright=NINT(maskSEM(D_DECL(i,j,k)))
       if ((maskleft.lt.0).or.(maskleft.gt.num_materials)) then
        print *,"maskleft invalid"
        stop
       endif
       if ((maskright.lt.0).or.(maskright.gt.num_materials)) then
        print *,"maskright invalid"
        stop
       endif
       if (((maskleft.ge.1).and.(maskleft.le.num_materials)).or. &
           ((maskright.ge.1).and.(maskright.le.num_materials))) then

        if ((maskleft.ne.0).and.(maskright.ne.0)) then
         if (maskleft.ne.maskright) then
          print *,"cannot have maskleft.ne.maskright"
          stop
         endif
        endif

        xmac(D_DECL(i,j,k))= &
          xmac(D_DECL(i,j,k))- &
          xface(D_DECL(i,j,k),FACECOMP_FACEDEN+1)* &
          deltafab(D_DECL(i,j,k))

       else if ((maskleft.eq.0).and.(maskright.eq.0)) then
        ! do nothing
       else
        print *,"maskleft or right invalid"
       endif

      enddo!i
      enddo!j
      enddo!k

      return
      end subroutine fort_semdeltaforce_face

       ! called from: NavierStokes::update_SEM_delta_force (NavierStokes.cpp)
       !   which is
       ! called from: NavierStokes::update_SEM_forces (MacProj.cpp)
       !   which is 
       ! called from: NavierStokes::update_SEM_forcesALL (MacProj.cpp)
       !   which is 
       ! called from: NavierStokes::veldiffuseALL (NavierStokes3.cpp)
       !   or 
       ! called from: NavierStokes::do_the_advance (NavierStokes3.cpp)
      subroutine fort_updatesemforce( &
       ns_time_order, &
       slab_step, &
       nsolve, &
       update_spectral, &
       update_stable, &
       nstate, &
       project_option, &
       xlo,dx,  &
       divfab,DIMS(divfab), &
       hoopfab,DIMS(hoopfab), &
       HOfab,DIMS(HOfab), &
       LOfab,DIMS(LOfab), &
       maskSEM,DIMS(maskSEM), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,level, &
       dt) &
      bind(c,name='fort_updatesemforce')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: ns_time_order
      integer, INTENT(in) :: slab_step
      integer, INTENT(in) :: nsolve
      integer, INTENT(in) :: update_spectral
      integer, INTENT(in) :: update_stable
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: project_option,level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(divfab)
      integer, INTENT(in) :: DIMDEC(hoopfab)
      integer, INTENT(in) :: DIMDEC(HOfab)
      integer, INTENT(in) :: DIMDEC(LOfab)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in),target :: divfab(DIMV(divfab),nsolve)
      real(amrex_real), pointer :: divfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: hoopfab(DIMV(hoopfab),nsolve)
      real(amrex_real), pointer :: hoopfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: HOfab(DIMV(HOfab),NSTATE_SDC)
      real(amrex_real), pointer :: HOfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: LOfab(DIMV(LOfab),NSTATE_SDC)
      real(amrex_real), pointer :: LOfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in) :: dt
      integer local_maskSEM
      integer velcomp
      integer :: veldir
      integer :: i,j,k


       ! in: fort_updatesemforce
      if ((ns_time_order.ge.2).and.(ns_time_order.le.32)) then
       ! do nothing
      else
       print *,"ns_time_order invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid56"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.eq.one) then
       ! do nothing
      else
       print *,"dt invalid; fort_updatesemforce"
       stop
      endif
      if ((project_option.eq.SOLVETYPE_HEAT).or. & !thermal conductivity
          (project_option.eq.SOLVETYPE_VISC)) then !viscosity
       ! do nothing
      else
       print *,"project_option invalid; fort_updatesemforce"
       stop
      endif
      if ((nsolve.ne.1).and. &
          (nsolve.ne.SDIM)) then
       print *,"nsolve invalid 1"
       stop
      endif

      divfab_ptr=>divfab
      call checkbound_array(fablo,fabhi,divfab_ptr,0,-1)
      hoopfab_ptr=>hoopfab
      call checkbound_array(fablo,fabhi,hoopfab_ptr,0,-1)
      HOfab_ptr=>HOfab
      call checkbound_array(fablo,fabhi,HOfab_ptr,0,-1)
      LOfab_ptr=>LOfab
      call checkbound_array(fablo,fabhi,LOfab_ptr,0,-1)
      maskSEM_ptr=>maskSEM
      call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
       if ((local_maskSEM.ge.1).and. &
           (local_maskSEM.le.num_materials)) then

        if (project_option.eq.SOLVETYPE_VISC) then !viscosity

         do veldir=1,SDIM

          if (update_spectral.eq.1) then
           if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
            print *,"slab_step invalid"
            stop
           endif
           ! I-scheme,thermal conduction,viscosity (-div(2 mu D)-force)
           ! HOfab=-div 2 mu D-HOOP_FORCE_MARK_MF
           HOfab(D_DECL(i,j,k),SEMDIFFUSE_U+veldir)= &
            divfab(D_DECL(i,j,k),veldir)- &
            hoopfab(D_DECL(i,j,k),veldir)
          else if (update_spectral.eq.0) then
           ! do nothing
          else
           print *,"update_spectral invalid"
           stop
          endif

          if (update_stable.eq.1) then
           if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
            print *,"slab_step invalid"
            stop
           endif
           ! I-scheme,thermal conduction,viscosity (-div(2 mu D)-force)
           ! LOfab=-div 2 mu D-HOOP_FORCE_MARK_MF
           LOfab(D_DECL(i,j,k),SEMDIFFUSE_U+veldir)= &
            divfab(D_DECL(i,j,k),veldir)- &
            hoopfab(D_DECL(i,j,k),veldir)
          else if (update_stable.eq.0) then
           ! do nothing
          else
           print *,"update_stable invalid"
           stop
          endif
         enddo ! veldir=1..sdim

        else if (project_option.eq.SOLVETYPE_HEAT) then  ! temperature diffusion

         velcomp=1

         if (update_spectral.eq.1) then
          if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
           print *,"slab_step invalid"
           stop
          endif

          ! I-scheme,thermal conduction,viscosity
          ! HOfab=-div k grad T
          HOfab(D_DECL(i,j,k),SEMDIFFUSE_T+1)= &
           divfab(D_DECL(i,j,k),velcomp)
         else if (update_spectral.eq.0) then
          ! do nothing
         else
          print *,"update_spectral invalid"
          stop
         endif

         if (update_stable.eq.1) then
          if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
           print *,"slab_step invalid"
           stop
          endif

          ! I-scheme,thermal conduction,viscosity
          ! LOfab=-div k grad T
          LOfab(D_DECL(i,j,k),SEMDIFFUSE_T+1)= &
           divfab(D_DECL(i,j,k),velcomp)
         else if (update_stable.eq.0) then
          ! do nothing
         else
          print *,"update_stable invalid"
          stop
         endif

        else
         print *,"project_option invalid; fort_updatesemforce"
         stop
        endif

       else if (local_maskSEM.eq.0) then
        ! do nothing
       else
        print *,"local_maskSEM invalid"
        stop
       endif
      enddo!i
      enddo!j
      enddo!k

      return
      end subroutine fort_updatesemforce

       ! delta_MF, stableF_MF, spectralF_MF are initialized
       ! to zero in NavierStokes::init_delta_SDC()
      subroutine fort_updatesemforce_face( &
       project_option, &
       ns_time_order, &
       dir, &
       slab_step, &
       update_spectral, &
       update_stable, &
       xlo,dx,  &
       gpfab,DIMS(gpfab), &
       HOfab,DIMS(HOfab), &
       LOfab,DIMS(LOfab), &
       maskSEM,DIMS(maskSEM), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       dt) &
      bind(c,name='fort_updatesemforce_face')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: ns_time_order
      integer, INTENT(in) :: dir
      integer, INTENT(in) :: slab_step
      integer, INTENT(in) :: update_spectral
      integer, INTENT(in) :: update_stable
      integer, INTENT(in) :: level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(gpfab)
      integer, INTENT(in) :: DIMDEC(HOfab)
      integer, INTENT(in) :: DIMDEC(LOfab)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in),target :: gpfab(DIMV(gpfab))
      real(amrex_real), pointer :: gpfab_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out),target :: HOfab(DIMV(HOfab))
      real(amrex_real), pointer :: HOfab_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out),target :: LOfab(DIMV(LOfab))
      real(amrex_real), pointer :: LOfab_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in) :: dt
      integer maskleft,maskright
      integer :: i,j,k
      integer ii,jj,kk
      integer im_crit

      if (project_option.eq.SOLVETYPE_PRES) then
       ! do nothing
      else if (project_option.eq.SOLVETYPE_HEAT) then ! thermal diffusion
       print *,"thermal diffusion force is only cell centered"
       stop
      else if (project_option.eq.SOLVETYPE_VISC) then ! viscosity
       print *,"viscosity force is only cell centered"
       stop
      else
       print *,"project_option invalid fort_updatesemforce_face"
       stop
      endif

      if ((ns_time_order.ge.2).and.(ns_time_order.le.32)) then
       ! do nothing
      else
       print *,"ns_time_order invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid57"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (dt.ne.one) then
       print *,"dt invalid"
       stop
      endif
      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid update sem force face"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid update sem force face 2"
       stop
      endif

      gpfab_ptr=>gpfab
      call checkbound_array1(fablo,fabhi,gpfab_ptr,0,dir)
      HOfab_ptr=>HOfab
      call checkbound_array1(fablo,fabhi,HOfab_ptr,0,dir)
      LOfab_ptr=>LOfab
      call checkbound_array1(fablo,fabhi,LOfab_ptr,0,dir)
      maskSEM_ptr=>maskSEM
      call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0,dir)
 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       maskleft=NINT(maskSEM(D_DECL(i-ii,j-jj,k-kk)))
       maskright=NINT(maskSEM(D_DECL(i,j,k)))
       if ((maskleft.lt.0).or.(maskleft.gt.num_materials)) then
        print *,"maskleft invalid"
        stop
       endif
       if ((maskright.lt.0).or.(maskright.gt.num_materials)) then
        print *,"maskright invalid"
        stop
       endif
       im_crit=0
       if ((maskleft.eq.0).and.(maskright.eq.0)) then
        ! do nothing
       else if (maskleft.eq.0) then
        im_crit=maskright
       else if (maskright.eq.0) then
        im_crit=maskleft
       else if (maskleft.eq.maskright) then
        im_crit=maskleft
       else
        print *,"maskleft or maskright invalid"
        stop
       endif

       if ((im_crit.ge.1).and.(im_crit.le.num_materials)) then

        if (update_spectral.eq.1) then
         if ((slab_step.lt.-1).or.(slab_step.ge.bfact_time_order)) then
          print *,"slab_step invalid"
          stop
         endif 
         HOfab(D_DECL(i,j,k))=gpfab(D_DECL(i,j,k))
        else if (update_spectral.eq.0) then
         ! do nothing
        else
         print *,"update_spectral invalid"
         stop
        endif

        if (update_stable.eq.1) then
         if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
          print *,"slab_step invalid"
          stop
         endif

         LOfab(D_DECL(i,j,k))=gpfab(D_DECL(i,j,k))
        else if (update_stable.eq.0) then
         ! do nothing
        else
         print *,"update_stable invalid"
         stop
        endif

       else if (im_crit.eq.0) then
        ! do nothing
       else
        print *,"im_crit invalid"
        stop
       endif

      enddo!i
      enddo!j
      enddo!k

      return
      end subroutine fort_updatesemforce_face


      subroutine fort_sdc_time_quad( &
       HOncomp, &
       LOncomp, &
       delta_ncomp, &
       nstate, &
       xlo,dx,  &
       deltafab,DIMS(deltafab), &
       HOfab,DIMS(HOfab), &
       LOfab,DIMS(LOfab), &
       maskSEM,DIMS(maskSEM), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       finest_level, &
       dt) &
      bind(c,name='fort_sdc_time_quad')
      use LegendreNodes
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: HOncomp
      integer, INTENT(in) :: LOncomp
      integer, INTENT(in) :: delta_ncomp
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(deltafab)
      integer, INTENT(in) :: DIMDEC(HOfab)
      integer, INTENT(in) :: DIMDEC(LOfab)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(inout), target :: deltafab(DIMV(deltafab),delta_ncomp)
      real(amrex_real), pointer :: deltafab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: HOfab(DIMV(HOfab),HOncomp)
      real(amrex_real), pointer :: HOfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: LOfab(DIMV(LOfab),LOncomp)
      real(amrex_real), pointer :: LOfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in) :: dt
      integer i,j,k
      integer local_maskSEM
      integer slab_step
      integer ibase,ibase2,icomp
      integer kinterval,jstencil,iquad,i1
      real(amrex_real) force_integral,dt_sub
      real(amrex_real) sanity_sum
      real(amrex_real) GQws(0:bfact_time_order,0:bfact_time_order-1, &
                  1:bfact_time_order)
      real(amrex_real) GQwsQUAD(0:bfact_time_order,1:bfact_time_order)
      real(amrex_real) yGL(0:bfact_time_order)
      real(amrex_real) ydiff(1:bfact_time_order)

      if (bfact.lt.1) then
       print *,"bfact invalid58"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
      if (HOncomp.ne.NSTATE_SDC*(bfact_time_order+1)) then
       print *,"HOncomp invalid"
       stop
      endif
      if (LOncomp.ne.NSTATE_SDC*bfact_time_order) then
       print *,"LOncomp invalid"
       stop
      endif
      if (delta_ncomp.ne.NSTATE_SDC*bfact_time_order) then
       print *,"delta_ncomp invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid sdc time quad"
       stop
      endif

      HOfab_ptr=>HOfab
      call checkbound_array(fablo,fabhi,HOfab_ptr,0,-1)
      LOfab_ptr=>LOfab
      call checkbound_array(fablo,fabhi,LOfab_ptr,0,-1)
      deltafab_ptr=>deltafab
      call checkbound_array(fablo,fabhi,deltafab_ptr,0,-1)
      maskSEM_ptr=>maskSEM
      call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      call SDC_GQweights(bfact_time_order,bfact_time_order,GQws)

      do kinterval=1,bfact_time_order
       sanity_sum=zero
       do jstencil=0,bfact_time_order
        GQwsQUAD(jstencil,kinterval)=zero
        do iquad=0,bfact_time_order-1
         GQwsQUAD(jstencil,kinterval)=GQwsQUAD(jstencil,kinterval)+ &
          cache_gauss_w(bfact_time_order,iquad,TMTYPE)* &
          GQws(jstencil,iquad,kinterval)
        enddo ! iquad
        sanity_sum=sanity_sum+GQwsQUAD(jstencil,kinterval)
       enddo  ! jstencil
       if (abs(sanity_sum-two).gt.EPS12) then
        print *,"SDC sanity check failed"
        stop
       endif
      enddo  ! kinterval

      do i1=0,bfact_time_order
       yGL(i1)=cache_gauss_lobatto(bfact_time_order,i1,TMTYPE)
      enddo
      do i1=1,bfact_time_order
       ydiff(i1)=yGL(i1)-yGL(i1-1)
       if (ydiff(i1).le.zero) then
        print *,"yGL bust"
        stop
       endif
      enddo ! i1

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))

       if (local_maskSEM.eq.0) then
        ! do nothing
       else if ((local_maskSEM.ge.1).and.(local_maskSEM.le.num_materials)) then
        ! do nothing
       else
        print *,"local_maskSEM invalid"
        stop
       endif

       do slab_step=1,bfact_time_order
        ibase=(slab_step-1)*NSTATE_SDC
        do icomp=1,NSTATE_SDC
         deltafab(D_DECL(i,j,k),ibase+icomp)=zero
        enddo
       enddo

       if (local_maskSEM.eq.0) then
        ! do nothing
       else if ((local_maskSEM.ge.1).and. &
                (local_maskSEM.le.num_materials)) then

        do slab_step=1,bfact_time_order

         ibase=(slab_step-1)*NSTATE_SDC

         do icomp=1,NSTATE_SDC
           force_integral=zero
           do jstencil=0,bfact_time_order
            ibase2=jstencil*NSTATE_SDC+icomp
            force_integral=force_integral+ &
             GQwsQUAD(jstencil,slab_step)* &
             HOfab(D_DECL(i,j,k),ibase2)
           enddo ! jstencil
           dt_sub=dt*ydiff(slab_step)/two

           deltafab(D_DECL(i,j,k),ibase+icomp)= &
            force_integral*dt_sub/two- &
            dt_sub*LOfab(D_DECL(i,j,k),ibase+icomp)
         enddo ! icomp=1..NSTATE_SDC

        enddo ! slab_step

       else
        print *,"local_maskSEM invalid"
        stop
       endif
 
      enddo!i
      enddo!j
      enddo!k

      return
      end subroutine fort_sdc_time_quad

      subroutine fort_sdc_time_quad_face( &
       dir, &
       HOncomp, &
       LOncomp, &
       delta_ncomp, &
       nstate, &
       xlo,dx,  &
       deltafab,DIMS(deltafab), &
       HOfab,DIMS(HOfab), &
       LOfab,DIMS(LOfab), &
       maskSEM,DIMS(maskSEM), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       finest_level, &
       dt) &
      bind(c,name='fort_sdc_time_quad_face')
      use LegendreNodes
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: dir
      integer, INTENT(in) :: HOncomp
      integer, INTENT(in) :: LOncomp
      integer, INTENT(in) :: delta_ncomp
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(deltafab)
      integer, INTENT(in) :: DIMDEC(HOfab)
      integer, INTENT(in) :: DIMDEC(LOfab)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(inout), target :: deltafab(DIMV(deltafab),delta_ncomp)
      real(amrex_real), pointer :: deltafab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: HOfab(DIMV(HOfab),HOncomp)
      real(amrex_real), pointer :: HOfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: LOfab(DIMV(LOfab),LOncomp)
      real(amrex_real), pointer :: LOfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in) :: dt
      integer i,j,k
      integer ii,jj,kk
      integer slab_step
      integer kinterval,jstencil,iquad,i1
      real(amrex_real) force_integral,dt_sub
      real(amrex_real) sanity_sum
      real(amrex_real) GQws(0:bfact_time_order,0:bfact_time_order-1, &
                  1:bfact_time_order)
      real(amrex_real) GQwsQUAD(0:bfact_time_order,1:bfact_time_order)
      real(amrex_real) yGL(0:bfact_time_order)
      real(amrex_real) ydiff(1:bfact_time_order)
      integer maskleft,maskright

      if (bfact.lt.1) then
       print *,"bfact invalid59"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
      if (HOncomp.ne.(bfact_time_order+1)) then
       print *,"HOncomp invalid"
       stop
      endif
      if (LOncomp.ne.bfact_time_order) then
       print *,"LOncomp invalid"
       stop
      endif
      if (delta_ncomp.ne.bfact_time_order) then
       print *,"delta_ncomp invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid sdc time quad face"
       stop
      endif

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid sdc time quad face"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid sdc time quad face 2"
       stop
      endif

      HOfab_ptr=>HOfab
      call checkbound_array(fablo,fabhi,HOfab_ptr,0,dir)
      LOfab_ptr=>LOfab
      call checkbound_array(fablo,fabhi,LOfab_ptr,0,dir)
      deltafab_ptr=>deltafab
      call checkbound_array(fablo,fabhi,deltafab_ptr,0,dir)
      maskSEM_ptr=>maskSEM
      call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0,dir)

      call SDC_GQweights(bfact_time_order,bfact_time_order,GQws)

      do kinterval=1,bfact_time_order
       sanity_sum=zero
       do jstencil=0,bfact_time_order
        GQwsQUAD(jstencil,kinterval)=zero
        do iquad=0,bfact_time_order-1
         GQwsQUAD(jstencil,kinterval)=GQwsQUAD(jstencil,kinterval)+ &
          cache_gauss_w(bfact_time_order,iquad,TMTYPE)* &
          GQws(jstencil,iquad,kinterval)
        enddo ! iquad
        sanity_sum=sanity_sum+GQwsQUAD(jstencil,kinterval)
       enddo  ! jstencil
       if (abs(sanity_sum-two).le.EPS12) then
        !do nothing
       else
        print *,"SDC sanity check failed: ",sanity_sum
        stop
       endif
      enddo  ! kinterval

      do i1=0,bfact_time_order
       yGL(i1)=cache_gauss_lobatto(bfact_time_order,i1,TMTYPE)
      enddo
      do i1=1,bfact_time_order
       ydiff(i1)=yGL(i1)-yGL(i1-1)
       if (ydiff(i1).le.zero) then
        print *,"yGL bust"
        stop
       endif
      enddo ! i1

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       maskleft=NINT(maskSEM(D_DECL(i-ii,j-jj,k-kk)))
       maskright=NINT(maskSEM(D_DECL(i,j,k)))
       if ((maskleft.lt.0).or.(maskleft.gt.num_materials)) then
        print *,"maskleft invalid"
        stop
       endif
       if ((maskright.lt.0).or.(maskright.gt.num_materials)) then
        print *,"maskright invalid"
        stop
       endif

       do slab_step=1,bfact_time_order
        deltafab(D_DECL(i,j,k),slab_step)=zero
       enddo
       if (((maskleft.ge.1).and.(maskleft.le.num_materials)).or. &
           ((maskright.ge.1).and.(maskright.le.num_materials))) then

        do slab_step=1,bfact_time_order

         force_integral=zero
         do jstencil=0,bfact_time_order
          force_integral=force_integral+ &
           GQwsQUAD(jstencil,slab_step)* &
           HOfab(D_DECL(i,j,k),jstencil+1)
         enddo ! jstencil
         dt_sub=dt*ydiff(slab_step)/two

         deltafab(D_DECL(i,j,k),slab_step)= &
          force_integral*dt_sub/two- &
          dt_sub*LOfab(D_DECL(i,j,k),slab_step)

        enddo ! slab_step

       endif
 
      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_sdc_time_quad_face

       ! viscoelastic force = div ( H(phi) mu_p (1/lambda)(f(A)A-I)) if
       !  FENE-P
       ! Pushpendra Singh (NJIT)
      subroutine fort_maketensor( &
       partid, & ! 0..num_materials_viscoelastic-1
       level, &
       finest_level, &
       ncomp_visc, &
       im_parm, & ! 0..num_materials-1
       xlo,dx, &
       visc,DIMS(visc), &
       tensor,DIMS(tensor), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       elastic_viscosity, &
       etaS, &
       elastic_time, &
       viscoelastic_model, &
       polymer_factor, &
       irz) &
      bind(c,name='fort_maketensor')

      use probcommon_module
      use global_utility_module
      use mass_transfer_module
      IMPLICIT NONE

      integer, INTENT(in) :: partid !0..num_materials_viscoelastic-1
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: ncomp_visc
      integer, INTENT(in) :: im_parm
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(visc)
      integer, INTENT(in) :: DIMDEC(tensor)
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact

      real(amrex_real), INTENT(in), target :: visc(DIMV(visc),ncomp_visc)
      real(amrex_real), pointer :: visc_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(inout), target :: tensor(DIMV(tensor), &
              ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: tensor_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in) :: elastic_viscosity,etaS
      real(amrex_real), INTENT(in) :: elastic_time,polymer_factor
      integer, INTENT(in) :: viscoelastic_model
      integer, INTENT(in) :: irz

      integer ii,jj
      real(amrex_real) Q(3,3),TQ(3,3),Q_plus_I(3,3)
      integer i,j,k
      integer dir_local
      integer im_elastic_p1
      real(amrex_real) xcenter(SDIM)

       ! Q=A-I
      real(amrex_real) trace_A
      real(amrex_real) determinant_factor
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      tensor_ptr=>tensor
      visc_ptr=>visc

      if (ENUM_NUM_TENSOR_TYPE.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ENUM_NUM_TENSOR_TYPE invalid"
       stop
      endif

      if ((partid.ge.0).and. &
          (partid.lt.num_materials_viscoelastic)) then
       ! do nothing
      else
       print *,"partid invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid fort_maketensor"
       stop
      endif
      if (finest_level.ne.fort_finest_level) then
       print *,"finest_level invalid fort_maketensor"
       stop
      endif

      if ((im_parm.lt.0).or. &
          (im_parm.ge.num_materials).or. &
          (is_rigid(im_parm+1).eq.1)) then
       print *,"im_parm invalid26"
       stop
      endif

      im_elastic_p1=im_parm+1

      if (ncomp_visc.ne.3*num_materials) then
       print *,"ncomp_visc invalid"
       stop
      endif
      if (polymer_factor.ge.zero) then !1/L
       ! do nothing
      else
       print *,"polymer_factor out of range"
       stop
      endif

      call checkbound_array(fablo,fabhi,visc_ptr,1,-1)
      call checkbound_array(fablo,fabhi,tensor_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir_local=1,SDIM
        xcenter(dir_local)=xsten(0,dir_local)
       enddo

       do ii=1,3
       do jj=1,3
        Q(ii,jj)=zero
       enddo
       enddo
       do dir_local=1,ENUM_NUM_TENSOR_TYPE
        call stress_index(dir_local,ii,jj)
        Q(ii,jj)=tensor(D_DECL(i,j,k),dir_local)
       enddo
       Q(2,1)=Q(1,2)
       Q(3,1)=Q(1,3)
       Q(3,2)=Q(2,3)

       if (viscoelastic_model.eq.0) then ! FENE-CR
        ! coeff=(visc-etaS)/(modtime+dt)
        ! modtime=max(0.0,elastic_time*(1-Tr(A)/L^2))
       else if (viscoelastic_model.eq.1) then ! Oldroyd-B
        ! coeff=(visc-etaS)/(modtime+dt)
        ! modtime=elastic_time
       else if (viscoelastic_model.eq.5) then ! FENE-P
        ! coeff=(visc-etaS)/(modtime+dt)
        ! modtime=max(0.0,elastic_time*(1-Tr(A)/L^2))
       else if (viscoelastic_model.eq.6) then ! linear PTT
        ! coeff=(visc-etaS)/(modtime+dt)
        ! modtime=elastic_time
       else if (viscoelastic_model.eq.3) then ! incremental model
        ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
        ! coeff=elastic_viscosity
       else if (viscoelastic_model.eq.7) then ! incremental model
        ! Xia, Lu, Tryggvason 2018
        ! coeff=elastic_viscosity*|Q+I|^{-5/6}
       else
        print *,"viscoelastic_model invalid: ",viscoelastic_model
        stop
       endif

       trace_A=zero
       do ii=1,3
        trace_A=trace_A+Q(ii,ii)+one
        if ((viscoelastic_model.eq.0).or. & !FENE-CR
            (viscoelastic_model.eq.1).or. & !OLDROYD-B
            (viscoelastic_model.eq.5).or. & !FENE-P
            (viscoelastic_model.eq.6)) then !linear PTT
         if (Q(ii,ii)+one.gt.zero) then
          ! do nothing
         else
          print *,"A=Q+I should be positive definite"
          stop
         endif
        else if (viscoelastic_model.eq.3) then ! incremental model, plastic
         ! do nothing
        else if (viscoelastic_model.eq.7) then ! incremental model, Neo-Hookian 
         if (Q(ii,ii)+one.gt.zero) then
          ! do nothing
         else
          print *,"A=Q+I should be positive definite"
          stop
         endif
        else
         print *,"viscoelastic_model invalid: ",viscoelastic_model
         stop
        endif
       enddo ! ii=1,3

       determinant_factor=one

       if (viscoelastic_model.eq.5) then !FENE-P          
        if (trace_A.gt.zero) then
         if (polymer_factor.gt.zero) then !1/L
             ! (f(A)A-I)/lambda  f(A)=1/(1-trac(A)/L^2)
             ! ((f(A)(Q+I)-I)/lambda
             ! (f(A)Q + f(A)I-I)/lambda=
             ! (f(A)/lambda)(Q+I-I/f(A))=
             ! (f(A)/lambda)(Q+(f(A)I-I)/f(A))
             ! (f(A)-1)/f(A)=1-1/f(A)=1-(1-trac(A)/L^2)=trace(A)/L^2
             ! (f(A)A-I)/lambda = (f(A)/lambda)*(Q+trace(A)I/L^{2})
             ! 
          do ii=1,3
           Q(ii,ii)=Q(ii,ii)+min(trace_A*(polymer_factor**2),one)
          enddo
         else
          print *,"polymer_factor out of range for FENE-P"
          stop
         endif
        else
         print *,"trace_A should be positive for FENE-P"
         stop
        endif
       else if (viscoelastic_model.eq.0) then !FENE-CR
        ! do nothing 
       else if (viscoelastic_model.eq.1) then !OLDROYD-B
        ! do nothing 
       else if (viscoelastic_model.eq.3) then ! incremental model,plastic
        ! do nothing
       else if (viscoelastic_model.eq.7) then ! incremental model,Neo-Hookean
        do ii=1,3
        do jj=1,3
         Q_plus_I(ii,jj)=Q(ii,jj)
        enddo
        enddo
        do ii=1,3
         Q_plus_I(ii,ii)=Q_plus_I(ii,ii)+one
        enddo
        call abs_value_determinant(Q_plus_I,3,determinant_factor)

        if (determinant_factor.gt.zero) then
!         determinant_factor=one/(determinant_factor**(five/six))
         determinant_factor=one  !*incompressible* Neo-Hookean model.
        else if (determinant_factor.eq.zero) then
         print *,"determinant_factor must be positive"
         do ii=1,3
         do jj=1,3
          print *,"ii,jj,Q_plus_I ",ii,jj,Q_plus_I(ii,jj)
         enddo
         enddo
         print *,"determinant_factor: ",determinant_factor
         stop
        else
         print *,"determinant_factor invalid"
         stop
        endif
        
       else if (viscoelastic_model.eq.6) then !linearPTT
        ! do nothing
       else
        print *,"viscoelastic_model invalid: ",viscoelastic_model
        stop
       endif

        ! in fort_derviscosity:
        !  etaS=etaL-etaP=viscconst-elastic_viscosity 
        !  viscoelastic_coeff= &
        !   (visc(D_DECL(i,j,k),im_parm+1)-etaS)/(modtime+dt)
        !  visc(D_DECL(i,j,k),num_materials+im_parm+1)=
        !     viscoelastic_coeff*visc_coef
        !  visc(D_DECL(i,j,k),2*num_materials+im_parm+1)=modtime
       do ii=1,3
       do jj=1,3
         ! FENE-P or FENE-CP visc(num_materials+im+1)=f(A)/lambda
        TQ(ii,jj)=Q(ii,jj)*determinant_factor* &
             visc(D_DECL(i,j,k),num_materials+im_parm+1)
       enddo
       enddo

       do dir_local=1,ENUM_NUM_TENSOR_TYPE
        call stress_index(dir_local,ii,jj)
        tensor(D_DECL(i,j,k),dir_local)=TQ(ii,jj)
       enddo
      enddo!i
      enddo!j
      enddo!k

      return
      end subroutine fort_maketensor

      subroutine fort_maketensor_mac( &
       flux_grid_type, &
       partid, & ! 0..num_materials_viscoelastic-1
       level, &
       finest_level, &
       ncomp_visc, &
       im_parm, & ! 0..num_materials-1
       xlo,dx, &
       visc,DIMS(visc), &
       tensor,DIMS(tensor), &
       tensorMAC,DIMS(tensorMAC), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       elastic_viscosity, &
       etaS, &
       elastic_time, &
       viscoelastic_model, &
       polymer_factor, &
       irz) &
      bind(c,name='fort_maketensor_mac')

      use probcommon_module
      use global_utility_module
      use mass_transfer_module
      IMPLICIT NONE

      integer, INTENT(in) :: flux_grid_type
      integer, INTENT(in) :: partid !0..num_materials_viscoelastic-1
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: ncomp_visc
      integer, INTENT(in) :: im_parm  ! 0..num_materials-1
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(visc)
      integer, INTENT(in) :: DIMDEC(tensor)
      integer, INTENT(in) :: DIMDEC(tensorMAC)
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact

      real(amrex_real), INTENT(in), target :: visc(DIMV(visc),ncomp_visc)
      real(amrex_real), pointer :: visc_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: tensor(DIMV(tensor), &
              ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: tensor_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out), target :: tensorMAC(DIMV(tensorMAC), &
              ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: tensorMAC_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in) :: elastic_viscosity,etaS
      real(amrex_real), INTENT(in) :: elastic_time,polymer_factor
      integer, INTENT(in) :: viscoelastic_model
      integer, INTENT(in) :: irz

      integer i,j,k
      integer dir_local
      integer im_elastic_p1
      real(amrex_real) xcenter(SDIM)

      type(deriv_from_grid_parm_type) :: data_in
      type(interp_from_grid_out_parm_type) :: data_out

      real(amrex_real), target :: cell_data_deriv(1)

      integer itensor
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      tensorMAC_ptr=>tensorMAC

      if (ENUM_NUM_TENSOR_TYPE.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ENUM_NUM_TENSOR_TYPE invalid"
       stop
      endif

      if ((partid.ge.0).and. &
          (partid.lt.num_materials_viscoelastic)) then
       ! do nothing
      else
       print *,"partid invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid MAKETENSOR"
       stop
      endif
      if (finest_level.ne.fort_finest_level) then
       print *,"finest_level invalid MAKETENSOR"
       stop
      endif
      if (polymer_factor.ge.zero) then !1/L
       ! do nothing
      else
       print *,"polymer_factor out of range"
       stop
      endif

      if ((im_parm.lt.0).or. &
          (im_parm.ge.num_materials).or. &
          (is_rigid(im_parm+1).eq.1)) then
       print *,"im_parm invalid26"
       stop
      endif

      if ((flux_grid_type.eq.0).or. &
          (flux_grid_type.eq.1).or. &
          ((flux_grid_type.eq.2).and.(SDIM.eq.3))) then
       ! do nothing
      else
       print *,"flux_grid_type invalid"
       stop
      endif

      im_elastic_p1=im_parm+1

      if (ncomp_visc.ne.3*num_materials) then
       print *,"ncomp_visc invalid"
       stop
      endif


      visc_ptr=>visc
      tensor_ptr=>tensor
      call checkbound_array(fablo,fabhi,visc_ptr,1,-1)
      call checkbound_array(fablo,fabhi,tensor_ptr,1,-1)
      call checkbound_array(fablo,fabhi,tensorMAC_ptr,0,flux_grid_type)

      !type(deriv_from_grid_parm_type) :: data_in
      data_in%level=level
      data_in%finest_level=finest_level
      data_in%bfact=bfact ! bfact=kind of spectral element grid 
      data_in%dx=dx
      data_in%xlo=xlo
      data_in%fablo=fablo
      data_in%fabhi=fabhi

      data_in%grid_type_flux=flux_grid_type
      call grid_type_to_box_type(flux_grid_type,data_in%box_type_flux)

      data_out%data_interp=>cell_data_deriv

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0, &
         flux_grid_type) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridstenMAC_level(xsten,i,j,k,level,nhalf,flux_grid_type)
       do dir_local=1,SDIM
        xcenter(dir_local)=xsten(0,dir_local)
       enddo

       data_in%index_flux(1)=i
       data_in%index_flux(2)=j
       if (SDIM.eq.3) then
        data_in%index_flux(SDIM)=k
       else if (SDIM.eq.2) then
        !do nothing
       else
        print *,"dimension bust"
        stop
       endif

       do itensor=1,ENUM_NUM_TENSOR_TYPE

         data_in%dir_deriv=-1 !interpolation
         data_in%grid_type_data=-1 !cell centered
         do dir_local=1,SDIM
          data_in%box_type_data(dir_local)=0
         enddo
         data_in%ncomp=1
         data_in%scomp=itensor

         call deriv_from_grid_util(data_in,tensor_ptr,data_out)

         tensorMAC(D_DECL(i,j,k),itensor)=cell_data_deriv(1)
       enddo ! itensor=1..ENUM_NUM_TENSOR_TYPE

      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_maketensor_mac


      subroutine fort_copy_vel_on_sign( &
       im_part, &
       nparts, &
       partid, &
       ngrow_make_distance_in, &
       nFSI, &
       xlo,dx, &
       snew,DIMS(snew), &
       fsi,DIMS(fsi), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       nstate) &
      bind(c,name='fort_copy_vel_on_sign')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: im_part
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: partid
      integer, INTENT(in) :: ngrow_make_distance_in
      integer, INTENT(in) :: nFSI
      integer, INTENT(in) :: nstate
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(fsi)
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact

      real(amrex_real), INTENT(inout),target :: snew(DIMV(snew),nstate)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: fsi(DIMV(fsi),nFSI)
      real(amrex_real), pointer :: fsi_ptr(D_DECL(:,:,:),:)

      integer i,j,k,dir,ibase
      real(amrex_real) LS

      snew_ptr=>snew
      fsi_ptr=>fsi

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (nFSI.ne.nparts*NCOMP_FSI) then
       print *,"nFSI invalid"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"ngrow_make_distance.ne.3"
       stop
      endif
      if (ngrow_make_distance_in.ne.3) then
       print *,"ngrow_make_distance_in.ne.3"
       stop
      endif
      if ((nparts.lt.1).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_copy_vel_on_sign"
       stop
      endif
      if ((partid.lt.0).or.(partid.ge.nparts)) then
       print *,"partid invalid"
       stop
      endif
      if ((im_part.lt.0).or. &
          (im_part.ge.num_materials).or. &
          (is_lag_part(im_part+1).ne.1)) then
       print *,"im_part invalid"
       stop
      endif

      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,fsi_ptr,ngrow_make_distance,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      ibase=partid*NCOMP_FSI

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       LS=fsi(D_DECL(i,j,k),ibase+FSI_LEVELSET+1)
       if (LS.ge.zero) then
        do dir=1,SDIM
         snew(D_DECL(i,j,k),dir)=fsi(D_DECL(i,j,k),ibase+FSI_VELOCITY+dir)
        enddo
       else if (LS.le.zero) then
        ! do nothing
       else
        print *,"LS bust"
        stop
       endif

      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_copy_vel_on_sign

      subroutine fort_build_moment( &
       cur_time, &
       level, &
       finest_level, &
       nFSI, &
       nparts, &
       ngrow_make_distance_in, &
       im_solid_map, &
       xlo,dx, &
       snew,DIMS(snew), &
       lsnew,DIMS(lsnew), &
       fsi,DIMS(fsi), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       nstate) &
      bind(c,name='fort_build_moment')

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: cur_time
      integer, INTENT(in) :: level 
      integer, INTENT(in) :: finest_level 
      integer, INTENT(in) :: nFSI
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: ngrow_make_distance_in
      integer, INTENT(in) :: im_solid_map(nparts)
      integer, INTENT(in) :: nstate
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(lsnew)
      integer, INTENT(in) :: DIMDEC(fsi)
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact

      real(amrex_real), INTENT(inout), target :: snew(DIMV(snew),nstate)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: lsnew(DIMV(lsnew),num_materials*(1+SDIM))
      real(amrex_real), pointer :: lsnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: fsi(DIMV(fsi),nFSI)
      real(amrex_real), pointer :: fsi_ptr(D_DECL(:,:,:),:)

      integer i,j,k
      integer i1,j1,k1
      integer im_part
      integer partid
      integer k1lo,k1hi
      integer dir,ibase
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) ldata(D_DECL(3,3,3))
      real(amrex_real) volume_frac,facearea
      real(amrex_real) centroid(SDIM)
      real(amrex_real) cencell(SDIM)
      real(amrex_real) volcell
      real(amrex_real) LS_center,mag,delta
      real(amrex_real) nrm(SDIM)
      integer vofcomp
      integer ok_to_modify_EUL

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (nFSI.ne.nparts*NCOMP_FSI) then
       print *,"nFSI invalid"
       stop
      endif
      if ((nparts.lt.1).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_build_moment"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"ngrow_make_distance invalid"
       stop
      endif
      if (ngrow_make_distance_in.ne.3) then
       print *,"ngrow_make_distance_in invalid"
       print *,"fort_build_moment"
       print *,"ngrow_make_distance_in=",ngrow_make_distance_in
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 33"
       stop
      endif
      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid"
       stop
      endif

      lsnew_ptr=>lsnew
      call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)
      snew_ptr=>snew
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      fsi_ptr=>fsi
      call checkbound_array(fablo,fabhi,fsi_ptr,ngrow_make_distance,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

       do partid=1,nparts
        im_part=im_solid_map(partid)+1
        if ((im_part.lt.1).or.(im_part.gt.num_materials)) then
         print *,"im_part invalid fort_build_moment"
         stop
        endif

        if (fort_read_from_CAD(FSI_flag(im_part)).eq.1) then

         if (is_lag_part(im_part).ne.1) then
          print *,"is_lag_part(im_part).ne.1"
          stop
         endif

         ok_to_modify_EUL=1
         if ((FSI_flag(im_part).eq.FSI_ICE_NODES_INIT).or. & 
             (FSI_flag(im_part).eq.FSI_FLUID_NODES_INIT)) then 
          if (cur_time.eq.zero) then
           ok_to_modify_EUL=1
          else if (cur_time.gt.zero) then
           ok_to_modify_EUL=0
          else
           print *,"cur_time invalid: ",cur_time
           stop
          endif
         else if ((FSI_flag(im_part).eq.FSI_PRESCRIBED_NODES).or. & 
                  (FSI_flag(im_part).eq.FSI_SHOELE_CTML)) then 
          ok_to_modify_EUL=1
         else
          print *,"FSI_flag invalid in fort_build_moment"
          print *,"im_part,FSI_flag(im_part): ",im_part,FSI_flag(im_part)
          stop
         endif

         if (ok_to_modify_EUL.eq.1) then

          ibase=(partid-1)*NCOMP_FSI
          do k1=k1lo,k1hi
          do j1=-1,1
          do i1=-1,1
           ldata(D_DECL(i1+2,j1+2,k1+2))= &
             fsi(D_DECL(i+i1,j+j1,k+k1),ibase+FSI_LEVELSET+1)
          enddo 
          enddo 
          enddo  ! i1,j1,k1=-1..1

          mag=zero
          do dir=1,SDIM
           i1=0
           j1=0
           k1=0
           if (dir.eq.1) then
            i1=1
           else if (dir.eq.2) then
            j1=1
           else if ((dir.eq.3).and.(SDIM.eq.3)) then
            k1=1
           else
            print *,"dir invalid"
            stop
           endif
           delta=xsten(2,dir)-xsten(-2,dir)
           if (delta.le.zero) then
            print *,"delta invalid delta=",delta
            print *,"dir=",dir
            print *,"xsten(2,dir)=",xsten(2,dir)
            print *,"xsten(-2,dir)=",xsten(-2,dir)
            print *,"i,j,k ",i,j,k
            print *,"xlo=",xlo(1),xlo(2),xlo(SDIM)
            print *,"fablo=",fablo(1),fablo(2),fablo(SDIM)
            print *,"bfact=",bfact
            print *,"nhalf=",nhalf
            print *,"dx=",dx(1),dx(2),dx(SDIM)
            stop
           endif
           nrm(dir)=(ldata(D_DECL(2+i1,2+j1,2+k1))- &
                     ldata(D_DECL(2-i1,2-j1,2-k1)))/delta
           mag=mag+nrm(dir)**2
          enddo ! dir=1..sdim
          mag=sqrt(mag)
          if (mag.eq.zero) then
           nrm(1)=one
          else if (mag.gt.zero) then
           do dir=1,SDIM
            nrm(dir)=nrm(dir)/mag
           enddo
          else
           print *,"mag invalid 6"
           stop
          endif

          do dir=1,SDIM
           lsnew(D_DECL(i,j,k),num_materials+(im_part-1)*SDIM+dir)=nrm(dir)
          enddo
  
          LS_center=fsi(D_DECL(i,j,k),ibase+FSI_LEVELSET+1)

          call getvolume( &
           bfact,dx,xsten,nhalf, &
           ldata,volume_frac,facearea, &
           centroid,VOFTOL,SDIM)
          call CISBOX(xsten,nhalf, &
           xlo,dx,i,j,k, &
           bfact,level, &
           volcell,cencell,SDIM)
          vofcomp=STATECOMP_MOF+(im_part-1)*ngeom_raw

          if (volume_frac.lt.zero) then
           print *,"volume_frac.lt.zero"
           stop
          else if (volume_frac.le.VOFTOL) then
           if (LS_center.ge.-VOFTOL*dx(1)) then
            volume_frac=VOFTOL_SLOPES
           endif
          else if (volume_frac.gt.one) then
           print *,"volume_frac.gt.one"
           stop
          else if (volume_frac.ge.one-VOFTOL) then
           if (LS_center.le.VOFTOL*dx(1)) then
            volume_frac=one-VOFTOL_SLOPES
           endif
          else if ((volume_frac.ge.VOFTOL).and. &
                   (volume_frac.le.one-VOFTOL)) then
           ! do nothing
          else
           print *,"volume_frac invalid: ",volume_frac
           stop
          endif
  
          snew(D_DECL(i,j,k),vofcomp+1)=volume_frac 
          do dir=1,SDIM
           snew(D_DECL(i,j,k),vofcomp+1+dir)=centroid(dir)-cencell(dir)
          enddo 

         else if (ok_to_modify_EUL.eq.0) then
          ! do nothing
         else
          print *,"ok_to_modify_EUL invalid"
          stop
         endif

        else if (FSI_flag(im_part).eq.FSI_PRESCRIBED_PROBF90) then 
         ! do nothing
        else
         print *,"FSI_flag invalid in fort_build_moment"
         print *,"im_part,FSI_flag(im_part): ",im_part,FSI_flag(im_part)
         stop
        endif
       enddo ! partid=1..nparts

      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_build_moment


        ! called from tensor_advecton_update() in NavierStokes.cpp
        ! vel is the advective velocity
      subroutine fort_updatetensor( &
       level, &
       finest_level, &
       im_critical, &  ! 0<=im_critical<=num_materials-1
       ncomp_visc, & 
       visc,DIMS(visc), &
       one_over_den, &
       DIMS(one_over_den), &
       tendata,DIMS(tendata), & !tendata:fort_getshear,only_scalar=0
       dx,xlo, &
       vel,DIMS(vel), &
       tnew,DIMS(tnew), &
       told,DIMS(told), &
       tilelo, tilehi,  &
       fablo, fabhi, &
       bfact,  &
       dt, &
       elastic_time, &
       viscoelastic_model, &
       polymer_factor, &
       elastic_viscosity, &
       irz, &
       bc, &
       transposegradu) &
      bind(c,name='fort_updatetensor')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: im_critical
      integer, INTENT(in) :: ncomp_visc
      integer, INTENT(in) :: DIMDEC(visc)
      integer, INTENT(in) :: DIMDEC(one_over_den)
      integer, INTENT(in) :: DIMDEC(tendata)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(tnew)
      integer, INTENT(in) :: DIMDEC(told)
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(SDIM),xlo(SDIM)

      real(amrex_real), INTENT(in), target :: visc(DIMV(visc),ncomp_visc)
      real(amrex_real), pointer :: visc_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: one_over_den(DIMV(one_over_den))
      real(amrex_real), pointer :: one_over_den_ptr(D_DECL(:,:,:))

      ! D=(1/2)(gradU + gradU^Transpose)
      ! DERIVE_TENSOR_MAG+1: sqrt(2 * D : D)
      ! DERIVE_TENSOR_RATE_DEFORM+1: D11,D12,D13,D21,D22,D23,D31,D32,D33
      ! DERIVE_TENSOR_GRAD_VEL+1: ux,uy,uz,vx,vy,vz,wx,wy,wz
      real(amrex_real), INTENT(in), target :: tendata(DIMV(tendata),DERIVE_TENSOR_NCOMP)
      real(amrex_real), pointer :: tendata_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: vel(DIMV(vel),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(out), target :: tnew(DIMV(tnew),ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: tnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real) :: point_tnew(ENUM_NUM_TENSOR_TYPE)

      real(amrex_real), INTENT(in), target :: told(DIMV(told),ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: told_ptr(D_DECL(:,:,:),:)
      real(amrex_real) :: point_told(ENUM_NUM_TENSOR_TYPE)

      integer :: i,j,k
      real(amrex_real), INTENT(in) :: dt,elastic_time
      integer, INTENT(in) :: viscoelastic_model
      real(amrex_real), INTENT(in) :: polymer_factor
      real(amrex_real), INTENT(in) :: elastic_viscosity
      integer, INTENT(in) :: transposegradu
      integer, INTENT(in) :: bc(SDIM,2,SDIM)
      integer, INTENT(in) :: irz
      integer :: dir_local

      tnew_ptr=>tnew

      if (irz.ne.levelrz) then
       print *,"irz invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid60"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 34"
       stop
      endif

      if (polymer_factor.ge.zero) then !1/L
       ! do nothing
      else
       print *,"polymer_factor out of range"
       stop
      endif
      if (elastic_viscosity.gt.zero) then
       ! do nothing
      else
       print *,"elastic_viscosity invalid"
       stop
      endif

      if (viscoelastic_model.eq.0) then ! FENE-CR
       ! coeff=(visc-etaS)/(modtime+dt)
       ! modtime=max(0.0,elastic_time*(1-Tr(A)/L^2))
      else if (viscoelastic_model.eq.1) then ! Oldroyd-B
       ! coeff=(visc-etaS)/(modtime+dt)
       ! modtime=elastic_time
      else if (viscoelastic_model.eq.5) then ! FENE-P
       ! coeff=(visc-etaS)/(modtime+dt)
       ! modtime=max(0.0,elastic_time*(1-Tr(A)/L^2))
      else if (viscoelastic_model.eq.6) then ! linear PTT
       ! coeff=(visc-etaS)/(modtime+dt)
       ! modtime=elastic_time
      else if (viscoelastic_model.eq.3) then ! incremental model
       ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
       ! coeff=elastic_viscosity
      else if (viscoelastic_model.eq.7) then ! incremental model
       ! Xia, Lu, Tryggvason 2018
       ! coeff=elastic_viscosity
      else
       print *,"viscoelastic_model invalid: ",viscoelastic_model
       stop
      endif

      if ((im_critical.lt.0).or.(im_critical.ge.num_materials)) then
       print *,"im_critical invalid27"
       stop
      endif
      if (ncomp_visc.ne.3*num_materials) then
       print *,"ncomp visc invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      visc_ptr=>visc
      call checkbound_array(fablo,fabhi,visc_ptr,0,-1)
      one_over_den_ptr=>one_over_den
      call checkbound_array1(fablo,fabhi,one_over_den_ptr,0,-1)
      tendata_ptr=>tendata
      call checkbound_array(fablo,fabhi,tendata_ptr,0,-1)
      vel_ptr=>vel
      call checkbound_array(fablo,fabhi,vel_ptr,1,-1)
      call checkbound_array(fablo,fabhi,tnew_ptr,0,-1)

      told_ptr=>told
      call checkbound_array(fablo,fabhi,told_ptr,0,-1)

      if ((transposegradu.ne.0).and. &
          (transposegradu.ne.1)) then
       print *,"transposegradu invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       do dir_local=1,ENUM_NUM_TENSOR_TYPE
        point_told(dir_local)=told(D_DECL(i,j,k),dir_local)
       enddo

        !point_updatetensor is declared in: GLOBALUTIL.F90
       call point_updatetensor( &
        i,j,k, &
        level, &
        finest_level, &
        im_critical, &  ! 0<=im_critical<=num_materials-1
        ncomp_visc, & 
        visc_ptr, &
        one_over_den_ptr, &
        tendata_ptr, & !tendata:fort_getshear,only_scalar=0
        dx,xlo, &
        vel_ptr, &
        point_tnew, &
        point_told, &
        tilelo, tilehi,  &
        fablo, fabhi, &
        bfact,  &
        dt, &
        elastic_time, &
        viscoelastic_model, &
        polymer_factor, &
        elastic_viscosity, &
        irz, &
        bc, &
        transposegradu) 

       do dir_local=1,ENUM_NUM_TENSOR_TYPE
        tnew(D_DECL(i,j,k),dir_local)=point_tnew(dir_local)
       enddo

      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_updatetensor

        ! called from tensor_advecton_update() in NavierStokes.cpp
      subroutine fort_extrapolate_tensor( &
       level, &
       finest_level, &
       im_critical, &  ! 0<=im_critical<=num_materials-1
       dx,xlo, &
       LS,DIMS(LS), &
       tnew,DIMS(tnew), &
       told,DIMS(told), &
       tilelo, tilehi,  &
       fablo, fabhi, &
       bfact) &
      bind(c,name='fort_extrapolate_tensor')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: im_critical
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(tnew)
      integer, INTENT(in) :: DIMDEC(told)
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(SDIM),xlo(SDIM)

      real(amrex_real), INTENT(in), target :: &
            LS(DIMV(LS),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(out), target :: tnew(DIMV(tnew),ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: tnew_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: told(DIMV(told),ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: told_ptr(D_DECL(:,:,:),:)

      integer :: i,j,k
      integer :: i1,j1,k1
      integer :: k1low,k1high

      integer :: dir_local
      integer, PARAMETER :: nhalf=3
      integer :: im
      integer :: im_local
      integer :: im_sten

      real(amrex_real) x_sten(-nhalf:nhalf,SDIM)
      real(amrex_real) x_extrap(-nhalf:nhalf,SDIM)
      real(amrex_real) LS_local(num_materials)
      real(amrex_real) LS_sten(num_materials)
      real(amrex_real) Q_extrap(ENUM_NUM_TENSOR_TYPE)
      real(amrex_real) wtsum
      real(amrex_real) wt_local

      k1low=0
      k1high=0
      if (SDIM.eq.3) then
       k1low=-2
       k1high=2
      else if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif

      tnew_ptr=>tnew

      if (ENUM_NUM_TENSOR_TYPE.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ENUM_NUM_TENSOR_TYPE INVALID"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid60"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 34"
       stop
      endif

      if ((im_critical.lt.0).or.(im_critical.ge.num_materials)) then
       print *,"im_critical invalid27"
       stop
      endif

      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,2,-1)

      call checkbound_array(fablo,fabhi,tnew_ptr,0,-1)

      told_ptr=>told
      call checkbound_array(fablo,fabhi,told_ptr,2,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten_level(x_sten,i,j,k,level,nhalf)

       do im=1,num_materials
        LS_local(im)=LS(D_DECL(i,j,k),im)
       enddo
       call get_primary_material(LS_local,im_local)

       if ((im_local.eq.im_critical+1).and. &
           (LS_local(im_critical+1).ge.zero)) then
        ! do nothing
       else if ((im_local.ne.im_critical+1).or. &
                (LS_local(im_critical+1).lt.zero)) then

        if ((im_local.ge.1).and.(im_local.le.num_materials)) then

         do dir_local=1,ENUM_NUM_TENSOR_TYPE
          Q_extrap(dir_local)=zero
         enddo
         wtsum=zero

         do k1=k1low,k1high
         do j1=-2,2
         do i1=-2,2
          do im=1,num_materials
           LS_sten(im)=LS(D_DECL(i+i1,j+j1,k+k1),im)
          enddo
          call get_primary_material(LS_sten,im_sten)
          if ((im_sten.eq.im_critical+1).and. &
              (LS_sten(im_critical+1).ge.zero)) then
           call gridsten_level(x_extrap,i+i1,j+j1,k+k1,level,nhalf)
           wt_local=zero
           do dir_local=1,SDIM
            wt_local=wt_local+(x_extrap(0,dir_local)-x_sten(0,dir_local))**2
           enddo
           if (wt_local.gt.zero) then
            wt_local=one/wt_local
           else
            print *,"wt_local invalid"
            stop
           endif
           wtsum=wtsum+wt_local
           do dir_local=1,ENUM_NUM_TENSOR_TYPE
            Q_extrap(dir_local)=Q_extrap(dir_local)+ &
               wt_local*Told(D_DECL(i+i1,j+j1,k+k1),dir_local)
           enddo
          else if ((im_sten.ne.im_critical+1).or. &
                   (LS_sten(im_critical+1).lt.zero)) then
           ! do nothing
          else
           print *,"im_sten or LS_sten invalid"
           stop
          endif
         enddo !k1
         enddo !j1
         enddo !i1
           
         if (wtsum.eq.zero) then
          ! do nothing
         else if (wtsum.gt.zero) then
          do dir_local=1,ENUM_NUM_TENSOR_TYPE
           Q_extrap(dir_local)=Q_extrap(dir_local)/wtsum
          enddo
         else
          print *,"wtsum invalid"
          stop
         endif

         do dir_local=1,ENUM_NUM_TENSOR_TYPE
          tnew(D_DECL(i,j,k),dir_local)=Q_extrap(dir_local)
         enddo

        else
         print *,"im_local invalid"
         stop
        endif
       else
        print *,"im_local or LS_local invalid"
        stop
       endif

      enddo !i
      enddo !j
      enddo !k

      return
      end subroutine fort_extrapolate_tensor

       ! adjust_temperature==1  modify temperature (Snew and coeff)
       ! adjust_temperature==0  modify coefficient (coeff)
       ! adjust_temperature==-1 modify heatx,heaty,heatz
       ! if project_option==SOLVETYPE_HEAT:
       !  heatxyz correspond to thermal diffusivity
       ! else if project_option>=SOLVETYPE_SPEC:
       !  heatxyz correspond to rho D
       !
      subroutine fort_stefansolver( &
       project_option, & 
       solidheat_flag, & ! 0=diffuse in solid 1=dirichlet 2=Neumann
       microlayer_size, & ! 1..num_materials
       microlayer_substrate, & ! 1..num_materials
       microlayer_temperature_substrate, & ! 1..num_materials
       adjust_temperature, &
       nstate, &
       ntsat, &
       nden, &
       freezing_model, &
       distribute_from_target, &
       saturation_temp, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       xlo,dx, &
       dt, &
       maskfab,DIMS(maskfab), &
       conductstate,DIMS(conductstate), & ! num_materials components
       STATEFAB,DIMS(STATEFAB), &
       TgammaFAB,DIMS(TgammaFAB), &
       swept,DIMS(swept), &
       LS,DIMS(LS),  &
       T_fab,DIMS(T_fab),  &
       TorY_fab,DIMS(TorY_fab),  &
       Snew,DIMS(Snew), & 
       DeDT,DIMS(DeDT), &
       den,DIMS(den), &
       coeff,DIMS(coeff), &
       vol,DIMS(vol), &
       heatx,DIMS(heatx), &
       heaty,DIMS(heaty), &
       heatz,DIMS(heatz), &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz) ) &
      bind(c,name='fort_stefansolver')

      use probf90_module
      use global_utility_module
      use MOF_routines_module
      use mass_transfer_module
      IMPLICIT NONE

      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: ntsat
      integer, INTENT(in) :: nden

      integer, INTENT(in) :: solidheat_flag
      real(amrex_real), INTENT(in) :: microlayer_size(num_materials)
      integer, INTENT(in) :: microlayer_substrate(num_materials)
      real(amrex_real), INTENT(in) :: microlayer_temperature_substrate(num_materials)
      
      integer, INTENT(in) :: adjust_temperature
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: DIMDEC(conductstate)
      integer, INTENT(in) :: DIMDEC(maskfab)
      integer, INTENT(in) :: DIMDEC(STATEFAB)
      integer, INTENT(in) :: DIMDEC(TgammaFAB)
      integer, INTENT(in) :: DIMDEC(swept)
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(T_fab)
      integer, INTENT(in) :: DIMDEC(TorY_fab)
      integer, INTENT(in) :: DIMDEC(Snew)
      integer, INTENT(in) :: DIMDEC(DeDT)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(coeff)
      integer, INTENT(in) :: DIMDEC(vol)
      integer, INTENT(in) :: DIMDEC(heatx)
      integer, INTENT(in) :: DIMDEC(heaty)
      integer, INTENT(in) :: DIMDEC(heatz)
      integer, INTENT(in) :: DIMDEC(areax)
      integer, INTENT(in) :: DIMDEC(areay)
      integer, INTENT(in) :: DIMDEC(areaz)
      real(amrex_real), target, INTENT(in) :: &
           conductstate(DIMV(conductstate),num_materials)
      real(amrex_real), pointer :: conductstate_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: maskfab(DIMV(maskfab)) 
      real(amrex_real), pointer :: maskfab_ptr(D_DECL(:,:,:))
      real(amrex_real), target, INTENT(in) :: STATEFAB(DIMV(STATEFAB),nden) 
      real(amrex_real), pointer :: STATEFAB_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: TgammaFAB(DIMV(TgammaFAB),ntsat) 
      real(amrex_real), pointer :: TgammaFAB_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: swept(DIMV(swept),num_materials)
      real(amrex_real), pointer :: swept_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: LS(DIMV(LS),num_materials*(SDIM+1))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: T_fab(DIMV(T_fab),num_materials)
      real(amrex_real), pointer :: T_fab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: TorY_fab(DIMV(TorY_fab),num_materials)
      real(amrex_real), pointer :: TorY_fab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(out) :: Snew(DIMV(Snew),nstate)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)

      ! 1/(rho cv) (cv=DeDT)
      real(amrex_real), target, INTENT(in) :: DeDT(DIMV(DeDT))  
      real(amrex_real), pointer :: DeDT_ptr(D_DECL(:,:,:))

      ! 1/den (i.e. den actually stores 1/den)
      real(amrex_real), target, INTENT(in) :: den(DIMV(den))  
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:))

       ! alphanovolume or outer_iter_pressure
      real(amrex_real), target, INTENT(out) :: coeff(DIMV(coeff))  
      real(amrex_real), pointer :: coeff_ptr(D_DECL(:,:,:))

      real(amrex_real), target, INTENT(in) :: vol(DIMV(vol))
      real(amrex_real), pointer :: vol_ptr(D_DECL(:,:,:))
       ! thermal conductivity
      real(amrex_real), target, INTENT(out) :: heatx(DIMV(heatx))
      real(amrex_real), pointer :: heatx_ptr(D_DECL(:,:,:))
      real(amrex_real), target, INTENT(out) :: heaty(DIMV(heaty))
      real(amrex_real), pointer :: heaty_ptr(D_DECL(:,:,:))
      real(amrex_real), target, INTENT(out) :: heatz(DIMV(heatz))
      real(amrex_real), pointer :: heatz_ptr(D_DECL(:,:,:))

      real(amrex_real), target, INTENT(in) :: areax(DIMV(areax))
      real(amrex_real), target, INTENT(in) :: areay(DIMV(areay))
      real(amrex_real), target, INTENT(in) :: areaz(DIMV(areaz))
      real(amrex_real), pointer :: areax_ptr(D_DECL(:,:,:))
      real(amrex_real), pointer :: areay_ptr(D_DECL(:,:,:))
      real(amrex_real), pointer :: areaz_ptr(D_DECL(:,:,:))

      integer i,j,k
      integer i1,j1,k1
      integer k1lo,k1hi
      integer dir,side
      integer ii,jj,kk
      integer ic,jc,kc
      integer iface,jface,kface
      integer at_interface
      integer im_loop
      integer im_primary  ! LS_center(im_primary)=max_im LS_center(im)
      integer im_side_primary ! LS_side(im_side_primary)=max_im LS_side(im)
      integer im_crit
      integer im_crit_save
      integer im_adjust
      integer im
      integer im_opp
      integer ireverse
      integer iten
      integer im_source
      integer im_dest
      integer im_source_substrate
      integer im_dest_substrate
      integer ireverse_crit
      integer iten_crit
      integer im_source_crit
      integer im_dest_crit
      integer im_source_substrate_crit
      integer im_dest_substrate_crit
      integer local_freezing_model
      integer distribute_from_targ
      real(amrex_real) TGRAD_MAX,TGRAD_test,TorY_test
      real(amrex_real) LL
      real(amrex_real) Tgamma
      real(amrex_real) TorYgamma_BC
      integer Tgamma_STATUS
      integer tsat_comp
      integer ngrow_tsat
      real(amrex_real) T_MIN(num_materials)
      real(amrex_real) T_MAX(num_materials)
      real(amrex_real) TorY_MIN(num_materials)
      real(amrex_real) TorY_MAX(num_materials)
      integer T_STATUS(num_materials)
      integer TorY_STATUS(num_materials)

      real(amrex_real) over_den,over_cv
      real(amrex_real) single_material_den
      real(amrex_real) local_vol
      real(amrex_real) original_coeff,delta_coeff,coeff_Tgamma
      real(amrex_real) aface
      real(amrex_real) LS_center(num_materials)
      real(amrex_real) LS_side(num_materials)
      real(amrex_real) LS_no_tess(num_materials)
      real(amrex_real) LS1,LS2
      real(amrex_real) theta
      real(amrex_real) theta_cutoff
      real(amrex_real) heatcoeff
      real(amrex_real) side_coeff,T_adjust
      integer tcomp
      real(amrex_real) SWEPTFACTOR,hx
      integer ncomp_per_tsat
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_side(-nhalf:nhalf,SDIM)
      real(amrex_real) x_interface(SDIM)
      integer dir_inner
      real(amrex_real) T_or_Y_min_sanity
      real(amrex_real) thermal_k(num_materials)
      integer maskcell

      snew_ptr=>snew
      coeff_ptr=>coeff
      heatx_ptr=>heatx
      heaty_ptr=>heaty
      heatz_ptr=>heatz
      TgammaFAB_ptr=>TgammaFAB
      conductstate_ptr=>conductstate

      theta_cutoff=0.001

      ncomp_per_tsat=EXTRAP_PER_TSAT
      if (ntsat.eq.EXTRAP_NCOMP_TSAT) then
       ! do nothing
      else
       print *,"nstat invalid"
       stop
      endif
      if (nden.eq.num_materials*num_state_material) then
       ! do nothing
      else
       print *,"nden invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid67"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base must be 2"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level or finest_level invalid stefan solver"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
       ! solidheat_flag==0 diffuse in solid
       ! solidheat_flag==1 dirichlet solid/fluid
       ! solidheat_flag==2 Neumann solid/fluid
      if ((solidheat_flag.lt.0).or. &
          (solidheat_flag.gt.2)) then
       print *,"solidheat_flag invalid"
       stop
      endif
      do im=1,num_materials
       im_crit=microlayer_substrate(im)
       if (im_crit.eq.0) then
        ! do nothing
       else if ((im_crit.ge.1).and.(im_crit.le.num_materials)) then
        if (is_rigid(im_crit).ne.1) then
         print *,"is_rigid(im_crit) invalid"
         stop
        endif
       else
        print *,"microlayer_substrate invalid"
        stop
       endif
      enddo ! im=1..num_materials

      if (project_option.eq.SOLVETYPE_HEAT) then ! thermal diffusion
       T_or_Y_min_sanity=zero
      else if ((project_option.ge.SOLVETYPE_SPEC).and. & ! species diffusion
               (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
       T_or_Y_min_sanity=zero
      else
       print *,"project_option invalid; fort_stefansolver"
       stop
      endif

      if (SDIM.eq.2) then 
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      call checkbound_array(fablo,fabhi,conductstate_ptr,1,-1)
      STATEFAB_ptr=>STATEFAB
      call checkbound_array(fablo,fabhi,STATEFAB_ptr,1,-1)
      call checkbound_array(fablo,fabhi,TgammaFAB_ptr,1,-1)

      maskfab_ptr=>maskfab
      call checkbound_array1(fablo,fabhi,maskfab_ptr,1,-1)

      swept_ptr=>swept
      call checkbound_array(fablo,fabhi,swept_ptr,0,-1)

      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      T_fab_ptr=>T_fab
      call checkbound_array(fablo,fabhi,T_fab_ptr,1,-1)
      TorY_fab_ptr=>TorY_fab
      call checkbound_array(fablo,fabhi,TorY_fab_ptr,1,-1)
      call checkbound_array(fablo,fabhi,Snew_ptr,1,-1)

      DeDT_ptr=>DeDT
      call checkbound_array1(fablo,fabhi,DeDT_ptr,1,-1) !1/(density * cv)
      den_ptr=>den
      call checkbound_array1(fablo,fabhi,den_ptr,1,-1)  !1/(density)
      call checkbound_array1(fablo,fabhi,coeff_ptr,0,-1)
      vol_ptr=>vol
      call checkbound_array1(fablo,fabhi,vol_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,heatx_ptr,0,0)
      call checkbound_array1(fablo,fabhi,heaty_ptr,0,1)
      call checkbound_array1(fablo,fabhi,heatz_ptr,0,SDIM-1)
      areax_ptr=>areax
      areay_ptr=>areay
      areaz_ptr=>areaz
      call checkbound_array1(fablo,fabhi,areax_ptr,0,0)
      call checkbound_array1(fablo,fabhi,areay_ptr,0,1)
      call checkbound_array1(fablo,fabhi,areaz_ptr,0,SDIM-1)
 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       maskcell=NINT(maskfab(D_DECL(i,j,k)))

       if (maskcell.eq.1) then

        call gridsten_level(xsten,i,j,k,level,nhalf)

        im_crit_save=-1 ! the dominant material in the center cell
        im_dest_crit=-1
        im_source_crit=-1
        im_dest_substrate_crit=-1
        im_source_substrate_crit=-1
        iten_crit=-1
        ireverse_crit=-1
        TGRAD_MAX=-1.0e+10

        do im=1,num_materials
         LS_no_tess(im)=LS(D_DECL(i,j,k),im)
         thermal_k(im)=conductstate(D_DECL(i,j,k),im)
        enddo
        call LS_tessellate(LS_no_tess,LS_center)
        im_primary=1
        do im=2,num_materials
         if (LS_center(im).gt.LS_center(im_primary)) then
          im_primary=im
         endif
        enddo

        if (LS_center(im_primary).ge.zero) then

         do im=1,num_materials
          T_MIN(im)=zero
          T_MAX(im)=zero
          TorY_MIN(im)=zero
          TorY_MAX(im)=zero
          T_STATUS(im)=0
          TorY_STATUS(im)=0
         enddo ! im=1..num_materials

         do k1=k1lo,k1hi
         do j1=-1,1 
         do i1=-1,1 
          do im=1,num_materials
           LS_no_tess(im)=LS(D_DECL(i+i1,j+j1,k+k1),im)
          enddo
          call LS_tessellate(LS_no_tess,LS_side)
          im_side_primary=1
          do im=2,num_materials
           if (LS_side(im).gt.LS_side(im_side_primary)) then
            im_side_primary=im
           endif
          enddo
          TorY_test=T_fab(D_DECL(i+i1,j+j1,k+k1),im_side_primary)
          if (TorY_test.ge.T_or_Y_min_sanity) then
           ! do nothing
          else
           print *,"TorY_test= ",TorY_test
           print *,"adjust_temperature=",adjust_temperature
           print *,"i,j,k ",i,j,k
           print *,"i1,j1,k1 ",i1,j1,k1
           print *,"im_side_primary ",im_side_primary
           print *,"im_primary ",im_primary
           do im=1,num_materials
            print *,"im,T_fab ",im, &
             T_fab(D_DECL(i+i1,j+j1,k+k1),im) 
           enddo
           print *,"level, finest_level ",level,finest_level
           print *,"TorY_test.le.zero fort_stefansolver"
           stop
          endif
          if (T_STATUS(im_side_primary).eq.0) then
           T_MIN(im_side_primary)=TorY_test
           T_MAX(im_side_primary)=TorY_test
          else if (T_STATUS(im_side_primary).eq.1) then
           if (TorY_test.gt.T_MAX(im_side_primary)) then
            T_MAX(im_side_primary)=TorY_test
           endif
           if (TorY_test.lt.T_MIN(im_side_primary)) then
            T_MIN(im_side_primary)=TorY_test
           endif
          else
           print *,"T_STATUS(im_side_primary) invalid"
           stop
          endif
          T_STATUS(im_side_primary)=1

          TorY_test=TorY_fab(D_DECL(i+i1,j+j1,k+k1),im_side_primary)
          if (TorY_test.ge.T_or_Y_min_sanity) then
           ! do nothing
          else
           print *,"TorY_test= ",TorY_test
           print *,"adjust_temperature=",adjust_temperature
           print *,"i,j,k ",i,j,k
           print *,"i1,j1,k1 ",i1,j1,k1
           print *,"im_side_primary ",im_side_primary
           print *,"im_primary ",im_primary
           do im=1,num_materials
            print *,"im,TorY_fab ",im, &
             TorY_fab(D_DECL(i+i1,j+j1,k+k1),im) 
           enddo
           print *,"level, finest_level ",level,finest_level
           print *,"TorY_test.le.zero fort_stefansolver"
           stop
          endif
          if (TorY_STATUS(im_side_primary).eq.0) then
           TorY_MIN(im_side_primary)=TorY_test
           TorY_MAX(im_side_primary)=TorY_test
          else if (TorY_STATUS(im_side_primary).eq.1) then
           if (TorY_test.gt.TorY_MAX(im_side_primary)) then
            TorY_MAX(im_side_primary)=TorY_test
           endif
           if (TorY_test.lt.TorY_MIN(im_side_primary)) then
            TorY_MIN(im_side_primary)=TorY_test
           endif
          else
           print *,"TorY_STATUS(im_side_primary) invalid"
           stop
          endif
          TorY_STATUS(im_side_primary)=1

         enddo
         enddo
         enddo ! i1,j1,k1

         do im=1,num_materials-1
          do im_opp=im+1,num_materials
           do ireverse=0,1
            if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
             print *,"im or im_opp bust 8"
             stop
            endif

            call get_iten(im,im_opp,iten)
            LL=get_user_latent_heat(iten+ireverse*num_interfaces, &
                    room_temperature,1)

            Tgamma_STATUS=NINT(TgammaFAB(D_DECL(i,j,k),iten))

            if (ireverse.eq.0) then
             ! do nothing
            else if (ireverse.eq.1) then
             Tgamma_STATUS=-Tgamma_STATUS
            else
             print *,"ireverse invalid"
             stop
            endif

            if ((Tgamma_STATUS.eq.1).or.(Tgamma_STATUS.eq.2)) then

             local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
             distribute_from_targ=distribute_from_target(iten+ireverse*num_interfaces)

             if ((distribute_from_targ.ne.0).and. &
                 (distribute_from_targ.ne.1)) then
              print *,"distribute_from_targ invalid"
              stop
             endif

             if ((is_rigid(im).eq.0).and. &
                 (is_rigid(im_opp).eq.0)) then

              if (LL.ne.zero) then

               if (ireverse.eq.0) then
                im_source=im
                im_dest=im_opp
               else if (ireverse.eq.1) then
                im_source=im_opp
                im_dest=im
               else
                print *,"ireverse invalid"
                stop
               endif

                ! local_freezing_model=0 (sharp interface stefan model)
                ! local_freezing_model=1 (source term model)
                ! local_freezing_model=2 (hydrate model)
                ! local_freezing_model=3 (wildfire)
                ! local_freezing_model=5 (evaporation/condensation)
                ! local_freezing_model=6 (Palmore Desjardins)
                ! local_freezing_model=7 (Cavitation)
               if (is_GFM_freezing_modelF(local_freezing_model).eq.1) then 

                if (project_option.eq.SOLVETYPE_HEAT) then
                   ! default Tgamma
                 Tgamma=saturation_temp(iten+ireverse*num_interfaces)
                 TorYgamma_BC=Tgamma
                 if (Tgamma.gt.zero) then
                  tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+1
                  Tgamma=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                  TorYgamma_BC=Tgamma
                  if (Tgamma.gt.zero) then
                   ! do nothing
                  else
                   print *,"Tgamma must be positive1"
                   stop
                  endif
                 else
                  print *,"saturation temperature must be positive2"
                  stop
                 endif
                else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                      (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
                 Tgamma=saturation_temp(iten+ireverse*num_interfaces)
                 TorYgamma_BC=one
                 if (Tgamma.gt.zero) then
                  tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+1
                  Tgamma=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                  tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+2
                  TorYgamma_BC=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                  if (Tgamma.gt.zero) then
                   ! do nothing
                  else
                   print *,"Tgamma must be positive22"
                   stop
                  endif
                  if ((TorYgamma_BC.ge.zero).and.(TorYgamma_BC.le.one)) then
                   ! do nothing
                  else
                   print *,"TorYgamma_BC (aka Y) must be >= 0 and <=1"
                   stop
                  endif
                 else
                  print *,"saturation temperature must be positive33"
                  stop
                 endif
                else
                 print *,"project_option invalid; fort_stefansolver"
                 stop
                endif

                im_source_substrate=im_source
                im_dest_substrate=im_dest

                if (local_freezing_model.eq.0) then ! stefan

                  ! Tgamma BC at thin filament interface.
                  if (solidheat_flag.eq.0) then ! diffuse in solid
                   if (microlayer_substrate(im_source).ne.0) then
                    im_source_substrate=microlayer_substrate(im_source)
                    if (is_rigid(im_source_substrate).ne.1) then
                     print *,"is_rigid(im_source_substrate) invalid"
                     stop
                    endif
                   endif
                   if (microlayer_substrate(im_dest).ne.0) then
                    im_dest_substrate=microlayer_substrate(im_dest)
                    if (is_rigid(im_dest_substrate).ne.1) then
                     print *,"is_rigid(im_dest_substrate) invalid"
                     stop
                    endif
                   endif
                  else if ((solidheat_flag.eq.1).or. & ! dirichlet 
                           (solidheat_flag.eq.2)) then ! neumann 
                   ! do nothing
                  else
                   print *,"solidheat_flag invalid"
                   stop
                  endif

                else if (local_freezing_model.eq.5) then ! stefan evap/diff
                  ! do nothing
                else if (local_freezing_model.eq.6) then !Palmore/Desjardins
                  ! do nothing
                else
                  print *,"local_freezing_model invalid 15"
                  stop
                endif

                if ((im_source.eq.im_primary).or. &
                    (im_source_substrate.eq.im_primary).or. &
                    (im_dest.eq.im_primary).or. &
                    (im_dest_substrate.eq.im_primary)) then

                  im_crit=im_primary

                   ! im_primary is found in the stencil
                  if (TorY_STATUS(im_crit).eq.1) then

                   if (LL.lt.zero) then ! freezing or condensation
                    TGRAD_test=Tgamma-T_MIN(im_crit)
                   else if (LL.gt.zero) then  ! melting or boiling
                    TGRAD_test=T_MAX(im_crit)-Tgamma
                   else
                    print *,"LL invalid"
                    stop
                   endif

                   if (DEBUG_EVAPORATION.eq.1) then
                    print *,"DEBUG_EVAPORATION STATEMENT 1"
                    print *,"i,j,k,x,y,z ",i,j,k, &
                        xsten(0,1),xsten(0,2),xsten(0,SDIM)
                    print *,"im,im_opp,ireverse ",im,im_opp,ireverse
                    print *,"im_source,im_dest ",im_source,im_dest
                    print *,"LL ",LL
                    print *,"project_option=",project_option
                    print *,"Tgamma,TorYgamma_BC ",Tgamma,TorYgamma_BC
                    print *,"im_crit=",im_crit
                    print *,"T_MIN(im_crit) ",T_MIN(im_crit)
                    print *,"T_MAX(im_crit) ",T_MAX(im_crit)
                    print *,"TorY_MIN(im_crit) ",TorY_MIN(im_crit)
                    print *,"TorY_MAX(im_crit) ",TorY_MAX(im_crit)
                    print *,"TGRAD_test=",TGRAD_TEST
                   endif

                    ! local_freezing_model=0 (sharp interface stefan model)
                    ! local_freezing_model=1 (source term model)
                    ! local_freezing_model=2 (hydrate model)
                    ! local_freezing_model=3 (wildfire)
                    ! local_freezing_model=5 (stefan evaporation/condensation)
                    ! local_freezing_model=6 (Palmore/Desjardins)
                    ! local_freezing_model=7 (Cavitation)
                   if ((local_freezing_model.eq.0).or. & !stefan model
                       ((local_freezing_model.eq.5).and. & !stefan:evap or cond
                        (TGRAD_test.gt.zero)).or. &
                       ((local_freezing_model.eq.6).and. & !Palmore/Desjardins
                        (TGRAD_test.gt.zero))) then 

                    if (im_dest_crit.eq.-1) then
                     im_crit_save=im_crit
                     im_dest_crit=im_dest
                     im_source_crit=im_source
                     im_dest_substrate_crit=im_dest_substrate
                     im_source_substrate_crit=im_source_substrate
                     iten_crit=iten
                     ireverse_crit=ireverse
                     TGRAD_MAX=TGRAD_test
                    else if ((im_dest_crit.ge.1).and. &
                             (im_dest_crit.le.num_materials)) then
                     if (TGRAD_test.gt.TGRAD_MAX) then
                      im_crit_save=im_crit
                      im_dest_crit=im_dest
                      im_source_crit=im_source
                      im_dest_substrate_crit=im_dest_substrate
                      im_source_substrate_crit=im_source_substrate
                      iten_crit=iten
                      ireverse_crit=ireverse
                      TGRAD_MAX=TGRAD_test
                     endif
                    else
                     print *,"im_dest_crit invalid"
                     stop
                    endif

                   else if ((local_freezing_model.eq.5).and. &
                            (TGRAD_test.le.zero)) then
                    ! do nothing
                   else if ((local_freezing_model.eq.6).and. &
                            (TGRAD_test.le.zero)) then
                    ! do nothing
                   else
                    print *,"local_freezing_model invalid 16"
                    stop
                   endif

                  else if (TorY_STATUS(im_crit).eq.0) then
                   ! do nothing
                  else
                   print *,"TorY_STATUS(im_crit) invalid"
                   stop
                  endif

                else if ((im_source.ne.im_primary).and. &
                         (im_source_substrate.ne.im_primary).and. &
                         (im_dest.ne.im_primary).and. &
                         (im_dest_substrate.ne.im_primary)) then
                  ! do nothing
                else
                 print *,"LS_center invalid"
                 stop
                endif

               else if (is_GFM_freezing_modelF(local_freezing_model).eq.0) then 
                ! do nothing
               else
                print *,"freezing_model invalid in stefansolver"
                print *,"local_freezing_model= ",local_freezing_model
                print *,"iten,ireverse,num_interfaces ", &
                        iten,ireverse,num_interfaces
                stop
               endif

              else if (LL.eq.zero) then
               ! do nothing
              else
               print *,"LL invalid"
               stop
              endif

             else if ((is_rigid(im).eq.1).or. &
                      (is_rigid(im_opp).eq.1)) then
              ! do nothing
             else
              print *,"is_rigid(im) or is_rigid(im_opp) invalid"
              stop
             endif

            else if ((Tgamma_STATUS.eq.-1).or.(Tgamma_STATUS.eq.-2)) then
             ! do nothing
            else if (Tgamma_STATUS.eq.0) then
             ! do nothing
            else
             print *,"Tgamma_STATUS invalid"
             stop
            endif

           enddo !ireverse=0,1
          enddo ! im_opp=im+1..num_materials
         enddo ! im=1..num_materials-1

         if (im_dest_crit.eq.-1) then
          ! do nothing
         else if ((im_dest_crit.ge.1).and. &
                  (im_dest_crit.le.num_materials)) then

          SWEPTFACTOR=swept(D_DECL(i,j,k),im_dest_crit) !default:SWEPTFACTOR==1
          if ((SWEPTFACTOR.ge.EPS2).and. &
              (SWEPTFACTOR.le.one)) then
           !do nothing
          else
           print *,"SWEPTFACTOR INVALID: ",SWEPTFACTOR
           stop
          endif
          over_den=den(D_DECL(i,j,k))  ! 1/(rho)
          over_cv=DeDT(D_DECL(i,j,k))  ! 1/(rho cv)
          local_vol=vol(D_DECL(i,j,k))
          single_material_den=STATEFAB(D_DECL(i,j,k), &
            (im_primary-1)*num_state_material+1)

          if ((over_den.gt.zero).and. &
              (over_cv.gt.zero).and. &
              (local_vol.gt.zero).and. &
              (single_material_den.gt.zero)) then
           ! do nothing
          else
           print *,"over_den, over_cv, local_vol, or single_mat_den invalid"
           stop
          endif

          original_coeff=one/(dt*SWEPTFACTOR)
          if (project_option.eq.SOLVETYPE_HEAT) then
           original_coeff=original_coeff/over_cv
          else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                   (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
           original_coeff=original_coeff/over_den
          else
           print *,"project_option invalid; fort_stefansolver"
           stop
          endif

          delta_coeff=zero
          coeff_Tgamma=zero

          im_crit=im_crit_save

          iten=iten_crit
          ireverse=ireverse_crit 

          im_source=im_source_crit
          im_source_substrate=im_source_substrate_crit
          im_dest=im_dest_crit
          im_dest_substrate=im_dest_substrate_crit

          LL=get_user_latent_heat(iten+ireverse*num_interfaces, &
                  room_temperature,1)
          local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
          distribute_from_targ= &
                distribute_from_target(iten+ireverse*num_interfaces)

          if (project_option.eq.SOLVETYPE_HEAT) then
           TorYgamma_BC=saturation_temp(iten+ireverse*num_interfaces)
          else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                   (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
           TorYgamma_BC=one
          else
           print *,"project_option invalid; subroutine fort_stefansolver"
           stop
          endif


          do dir=1,SDIM
            ii=0
            jj=0
            kk=0
            if (dir.eq.1) then
             ii=1
            else if (dir.eq.2) then
             jj=1
            else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
             kk=1
            else
             print *,"dir invalid stefansolver"
             stop
            endif 
            do side=-1,1,2
             if (side.eq.-1) then
              ic=i-ii
              jc=j-jj
              kc=k-kk
              iface=i
              jface=j
              kface=k
             else if (side.eq.1) then
              ic=i+ii
              jc=j+jj
              kc=k+kk
              iface=i+ii
              jface=j+jj
              kface=k+kk
             else
              print *,"side invalid"
              stop
             endif

             call gridsten_level(xsten_side,ic,jc,kc,level,nhalf)

             if (dir.eq.1) then
              aface=areax(D_DECL(iface,jface,kface))
             else if (dir.eq.2) then
              aface=areay(D_DECL(iface,jface,kface))
             else if ((dir.eq.3).and.(SDIM.eq.3)) then
              aface=areaz(D_DECL(iface,jface,kface))
             else
              print *,"dir invalid stefansolver 2"
              stop
             endif

             do im_loop=1,num_materials
              LS_no_tess(im_loop)=LS(D_DECL(ic,jc,kc),im_loop)
             enddo
             call LS_tessellate(LS_no_tess,LS_side)

             ! thermal diffusivity==0 where LS changes sign
             ! and latent_heat <> 0.
             at_interface=0
             if ((LS_center(im_source)*LS_side(im_source).le.zero).and. &
                 (LS_center(im_dest)*LS_side(im_dest).le.zero)) then
              at_interface=1
              LS1=LS_center(im_source)-LS_center(im_dest)
              LS2=LS_side(im_source)-LS_side(im_dest)
              call get_default_scalar_diffusion(project_option, &
                      thermal_k, &
                      LS1,im_source,im_dest, &
                      single_material_den, &
                      heatcoeff)
             else if ((LS_center(im_source_substrate)* &
                       LS_side(im_source_substrate).le.zero).and. &
                      (LS_center(im_dest_substrate)* &
                       LS_side(im_dest_substrate).le.zero)) then
              at_interface=1
              LS1=LS_center(im_source_substrate)-LS_center(im_dest_substrate)
              LS2=LS_side(im_source_substrate)-LS_side(im_dest_substrate)
              call get_default_scalar_diffusion(project_option, &
                      thermal_k, &
                      LS1,im_source_substrate,im_dest_substrate, &
                      single_material_den, &
                      heatcoeff)
             endif

             if (at_interface.eq.1) then

               ! cannot do tiling here.
              if (adjust_temperature.eq.-1) then ! modify heatx,heaty,heatz

               if (dir.eq.1) then
                heatx(D_DECL(iface,jface,kface))=zero
               else if (dir.eq.2) then
                heaty(D_DECL(iface,jface,kface))=zero
               else if ((dir.eq.3).and.(SDIM.eq.3)) then
                heatz(D_DECL(iface,jface,kface))=zero
               else
                print *,"dir invalid stefansolver 2B"
                stop
               endif
              
              else if ((adjust_temperature.eq.0).or. & ! modify coeff
                       (adjust_temperature.eq.1)) then ! modify Snew and coeff
               ! do nothing
              else
               print *,"adjust_temperature invalid"
               stop
              endif
         
              if (heatcoeff.lt.zero) then
               print *,"heatcoeff invalid"
               stop
              endif

              if (LS1*LS2.le.zero) then

               if ((LS1.eq.zero).and.(LS2.eq.zero)) then
                theta=half
               else if ((LS1.ne.zero).or.(LS2.ne.zero)) then
                theta=LS1/(LS1-LS2)
               else
                print *,"LS1 or LS2 invalid"
                stop
               endif
               if ((theta.ge.zero).and.(theta.le.one+EPS2)) then
                ! do nothing
               else
                print *,"theta invalid"
                stop
               endif
               if (theta.le.theta_cutoff) then
                theta=theta_cutoff
               endif

               do dir_inner=1,SDIM
                x_interface(dir_inner)=theta*xsten_side(0,dir_inner)+ &
                        (one-theta)*xsten(0,dir_inner)
               enddo

               if (project_option.eq.SOLVETYPE_HEAT) then
                tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+1
               else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                        (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
                tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+2
               else
                print *,"project_option invalid; fort_stefansolver"
                stop
               endif

               ngrow_tsat=1
               call interpfab_tsat( &
                i,j,k, &
                ireverse, &
                iten, &
                ntsat, &
                bfact, &
                level, &
                finest_level, &
                dx,xlo, &
                x_interface, &
                tsat_comp, &
                fablo,fabhi, &
                TgammaFAB_ptr, &
                TorYgamma_BC)  ! TorYgamma_BC here is an output

               hx=abs(xsten(0,dir)-xsten(2*side,dir))
               if ((levelrz.eq.COORDSYS_CYLINDRICAL).and.(dir.eq.2)) then
                hx=hx*xsten(side,1)
               endif
               side_coeff=aface*heatcoeff/(theta*hx)
               delta_coeff=delta_coeff+side_coeff
               coeff_Tgamma=coeff_Tgamma+TorYgamma_BC*side_coeff

               if (DEBUG_EVAPORATION.eq.1) then
                print *,"DEBUG_EVAPORATION STATEMENT 2"
                print *,"project_option,i,j,k,dir,side ", &
                        project_option,i,j,k,dir,side
                print *,"im_source,im_dest,TorYgamma_BC ", &
                        im_source,im_dest,TorYgamma_BC
               endif
              endif ! LS1 * LS2 <=0

             else if (at_interface.eq.0) then
              ! do nothing
             else
              print *,"at_interface invalid in fort_stefansolver"
              print *,"project_option=",project_option
              print *,"solidheat_flag=",solidheat_flag
              stop
             endif

            enddo ! side=-1,1,2
          enddo ! dir=1..sdim
        
          if (adjust_temperature.eq.-1) then ! modify heatxyz
            ! do nothing
          else if ((adjust_temperature.eq.0).or. & ! modify coeff
                   (adjust_temperature.eq.1)) then ! modify Snew and coeff
  
            if (delta_coeff.gt.zero) then

              ! im_crit is material that dominates center cell.
             if ((im_crit.lt.1).or. &
                 (im_crit.gt.num_materials)) then
              print *,"im_crit invalid"
              stop
             endif

             delta_coeff=delta_coeff/local_vol
             coeff_Tgamma=coeff_Tgamma/local_vol
    
             if (adjust_temperature.eq.1) then

              TorY_test=TorY_fab(D_DECL(i,j,k),im_crit)
              if (TorY_test.ge.T_or_Y_min_sanity) then
               ! do nothing
              else
               print *,"TorY_test<T_or_Y_min_sanity"
               stop
              endif

              T_adjust=(original_coeff*TorY_test+coeff_Tgamma)/ &
                       (original_coeff+delta_coeff) 

              do im_adjust=1,num_materials
               if (project_option.eq.SOLVETYPE_HEAT) then
                tcomp=STATECOMP_STATES+ &
                 (im_adjust-1)*num_state_material+ &
                 ENUM_TEMPERATUREVAR+1
               else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                        (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
                tcomp=STATECOMP_STATES+ &
                  (im_adjust-1)*num_state_material+ &
                  ENUM_SPECIESVAR+1+ &
                  project_option-SOLVETYPE_SPEC
               else
                print *,"project_option invalid; fort_stefansolver"
                stop
               endif
               Snew(D_DECL(i,j,k),tcomp)=T_adjust
              enddo  ! im_adjust=1..num_materials

              coeff(D_DECL(i,j,k))=T_adjust

             else if (adjust_temperature.eq.0) then

              coeff(D_DECL(i,j,k))=original_coeff+delta_coeff

             else
              print *,"adjust_temperature invalid"
              stop
             endif  

            else if (delta_coeff.eq.zero) then
             ! do nothing
            else
             print *,"delta_coeff invalid"
             stop
            endif

          else
            print *,"adjust_temperature invalid"
            stop
          endif

         else 
          print *,"im_dest_crit invalid"
          stop
         endif

        else if (LS_center(im_primary).lt.zero) then
         ! do nothing (vacuum should not be found)
        else
         print *,"LS_center(im_primary) invalid"
         stop
        endif

       else if (maskcell.eq.0) then
        ! do nothing
       else
        print *,"maskcell invalid"
        stop
       endif

      enddo ! i=growlo(1),growhi(1)
      enddo ! j=growlo(2),growhi(2)
      enddo ! k=growlo(3),growhi(3)

      return
      end subroutine fort_stefansolver

! MEHDI VAHAB HEAT SOURCE
! T^new=T^* + dt A Q/(rho cv V) 
! Q units: J/(m^2 s)
      subroutine fort_heatsource_face( &
       nstate, &
       saturation_temp, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       dt, &
       time, &
       level, &
       finest_level, &
       LS,DIMS(LS),  &
       Snew,DIMS(Snew), & 
       DeDT,DIMS(DeDT), &
       den,DIMS(den), &
       vol,DIMS(vol), &
       heatx,DIMS(heatx), &
       heaty,DIMS(heaty), &
       heatz,DIMS(heatz), &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz) ) &
      bind(c,name='fort_heatsource_face')

      use probf90_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nstate
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(Snew)
      integer, INTENT(in) :: DIMDEC(DeDT)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(vol)
      integer, INTENT(in) :: DIMDEC(heatx)
      integer, INTENT(in) :: DIMDEC(heaty)
      integer, INTENT(in) :: DIMDEC(heatz)
      integer, INTENT(in) :: DIMDEC(areax)
      integer, INTENT(in) :: DIMDEC(areay)
      integer, INTENT(in) :: DIMDEC(areaz)
      real(amrex_real), INTENT(in),target :: LS(DIMV(LS),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: Snew(DIMV(Snew),nstate)
      real(amrex_real), pointer :: Snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: DeDT(DIMV(DeDT))  ! 1/(rho cv) (cv=DeDT)
      real(amrex_real), pointer :: DeDT_ptr(D_DECL(:,:,:))
       ! 1/den (i.e. den actually stores 1/den)
      real(amrex_real), INTENT(in),target :: den(DIMV(den)) 
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: vol(DIMV(vol))
      real(amrex_real), pointer :: vol_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: heatx(DIMV(heatx))
      real(amrex_real), INTENT(in),target :: heaty(DIMV(heaty))
      real(amrex_real), INTENT(in),target :: heatz(DIMV(heatz))
      real(amrex_real), pointer :: heatx_ptr(D_DECL(:,:,:))
      real(amrex_real), pointer :: heaty_ptr(D_DECL(:,:,:))
      real(amrex_real), pointer :: heatz_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: areax(DIMV(areax))
      real(amrex_real), INTENT(in),target :: areay(DIMV(areay))
      real(amrex_real), INTENT(in),target :: areaz(DIMV(areaz))
      real(amrex_real), pointer :: areax_ptr(D_DECL(:,:,:))
      real(amrex_real), pointer :: areay_ptr(D_DECL(:,:,:))
      real(amrex_real), pointer :: areaz_ptr(D_DECL(:,:,:))

      integer i,j,k
      integer heat_dir
      integer heat_side
      integer ii,jj,kk
      integer iface,jface,kface
      integer icell,jcell,kcell
      integer im
      integer im_primary
      integer im_primary_cell
      real(amrex_real) over_den,over_cv,local_vol
      real(amrex_real) aface,hface
      integer tcomp
      real(amrex_real) heat_source_term
      real(amrex_real) heat_flux
      real(amrex_real) flux_sign
      real(amrex_real) ls_cell_or_face(num_materials)
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid68"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base must be 2"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid heat source face"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif

      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      Snew_ptr=>Snew
      call checkbound_array(fablo,fabhi,Snew_ptr,1,-1)
      DeDT_ptr=>DeDT
      call checkbound_array1(fablo,fabhi,DeDT_ptr,0,-1)
      den_ptr=>den
      call checkbound_array1(fablo,fabhi,den_ptr,0,-1)
      vol_ptr=>vol
      call checkbound_array1(fablo,fabhi,vol_ptr,0,-1)

       ! thermal conductivity
      heatx_ptr=>heatx
      heaty_ptr=>heaty
      heatz_ptr=>heatz
      call checkbound_array1(fablo,fabhi,heatx_ptr,0,0)
      call checkbound_array1(fablo,fabhi,heaty_ptr,0,1)
      call checkbound_array1(fablo,fabhi,heatz_ptr,0,SDIM-1)

      areax_ptr=>areax
      areay_ptr=>areay
      areaz_ptr=>areaz
      call checkbound_array1(fablo,fabhi,areax_ptr,0,0)
      call checkbound_array1(fablo,fabhi,areay_ptr,0,1)
      call checkbound_array1(fablo,fabhi,areaz_ptr,0,SDIM-1)
 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       ! MEHDI VAHAB HEAT SOURCE
       ! heat_dir=1,2,3
       ! heat_side=1,2
       if (is_in_probtype_list().eq.1) then
        call SUB_EB_heat_source(time,dt,xsten,nhalf, &
               heat_flux,heat_dir,heat_side)
       else
        heat_flux=zero
       endif

       if (heat_flux.gt.zero) then
 
        do im=1,num_materials
         ls_cell_or_face(im)=LS(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(ls_cell_or_face,im_primary)

        if (is_rigid(im_primary).eq.1) then

         ii=0
         jj=0
         kk=0
         if (heat_dir.eq.1) then
          ii=1
         else if (heat_dir.eq.2) then
          jj=1
         else if ((heat_dir.eq.3).and.(SDIM.eq.3)) then
          kk=1
         else
          print *,"heat_dir invalid heatsource"
          stop
         endif
         if (heat_side.eq.2) then
          iface=i+ii
          jface=j+jj
          kface=k+kk
          icell=i+ii
          jcell=j+jj
          kcell=k+kk
          flux_sign=one
         else if (heat_side.eq.1) then
          iface=i
          jface=j
          kface=k
          icell=i-ii
          jcell=j-jj
          kcell=k-kk
          flux_sign=-one
         else
          print *,"heat_side invalid"
          stop
         endif

         do im=1,num_materials
          ls_cell_or_face(im)=LS(D_DECL(icell,jcell,kcell),im)
         enddo
         call get_primary_material(ls_cell_or_face,im_primary_cell)

         if (is_rigid(im_primary_cell).eq.0) then

           ! in: subroutine fort_heatsource_face
          if (heat_dir.eq.1) then
           aface=areax(D_DECL(iface,jface,kface))
           hface=heatx(D_DECL(iface,jface,kface))
          else if (heat_dir.eq.2) then
           aface=areay(D_DECL(iface,jface,kface))
           hface=heaty(D_DECL(iface,jface,kface))
          else if ((heat_dir.eq.3).and.(SDIM.eq.3)) then
           aface=areaz(D_DECL(iface,jface,kface))
           hface=heatz(D_DECL(iface,jface,kface))
          else
           print *,"heat_dir invalid heatsource 2"
           stop
          endif

          if ((aface.ge.zero).and.(hface.ge.zero)) then
           ! do nothing
          else
           print *,"aface or hface (thermal conductivity) invalid"
           stop
          endif

          over_den=den(D_DECL(i,j,k))
          over_cv=DeDT(D_DECL(i,j,k))  ! 1/(rho cv)
          local_vol=vol(D_DECL(i,j,k))

          if ((over_den.gt.zero).and. &
              (over_cv.gt.zero).and. &
              (local_vol.gt.zero)) then
           ! do nothing
          else
           print *,"over_den, over_cv, or local_vol invalid"
           stop
          endif

          heat_source_term=flux_sign*dt*over_cv*aface*heat_flux/local_vol

          tcomp=STATECOMP_STATES+ENUM_TEMPERATUREVAR+1

          Snew(D_DECL(i,j,k),tcomp)= &
           Snew(D_DECL(i,j,k),tcomp)+heat_source_term

         else if (is_rigid(im_primary_cell).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(im_primary_cell) invalid"
          stop
         endif 
        else if (is_rigid(im_primary).eq.0) then
         ! do nothing
        else
         print *,"is_rigid(im_primary) invalid"
         stop
        endif 
       else if (heat_flux.eq.zero) then
        ! do nothing
       else
        print *,"heat_flux invalid"
        stop
       endif 

      enddo 
      enddo 
      enddo 

      return
      end subroutine fort_heatsource_face

       ! called from:NavierStokes::init_FSI_GHOST_MAC_MF(int ngrow) 
       ! (in NavierStokes.cpp)
       ! called when "law_of_the_wall=0,1,2"
       ! if nparts==0 => interpolate state cell velocity to MAC grid.
       ! if nparts>0 and law_of_the_wall==0 => interpolate solid cell velocity
       ! to MAC grid.
      subroutine fort_wallfunction( &
       data_dir, &
       law_of_the_wall, &
       NS_sumdata_size, &
       NS_sumdata, &
       ncomp_sum_int_user1, &
       ncomp_sum_int_user2, &
       wall_model_velocity, &
       im_solid_map, &
       level, &
       finest_level, &
       ngrow_distance_in, &
       nparts, &
       nparts_ghost, &
       nden, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx, &
       dt, &
       time, &
       LSCP,DIMS(LSCP),  &
       LSFD,DIMS(LSFD),  &
       state,DIMS(state), &
       ufluid,DIMS(ufluid), &
       usolid,DIMS(usolid), &
       ughost,DIMS(ughost), &
       history_dat, &
       DIMS(history_dat), &
       nhistory, &
       visc_coef) &
      bind(c,name='fort_wallfunction')
      use probf90_module
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: ncomp_sum_int_user1
      integer, INTENT(in) :: ncomp_sum_int_user2
      integer :: ncomp_sum_int_user12
      integer, INTENT(in) :: data_dir
      integer, INTENT(in) :: nhistory
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: ngrow_distance_in
      integer, INTENT(in) :: NS_sumdata_size
      integer, INTENT(in) :: law_of_the_wall(num_materials)
      real(amrex_real), INTENT(in) :: wall_model_velocity(num_materials)
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_ghost
      integer, INTENT(in) :: nden
      integer, INTENT(in) :: im_solid_map(nparts_ghost)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in), target :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: NS_sumdata(NS_sumdata_size)
      real(amrex_real), INTENT(in), target :: xlo(SDIM)
      real(amrex_real), INTENT(in), target :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: visc_coef
       ! DIMDEC is defined in ArrayLim.H in the BoxLib/Src/C_BaseLib
      integer, INTENT(in) :: DIMDEC(LSCP)
      integer, INTENT(in) :: DIMDEC(LSFD)
      integer, INTENT(in) :: DIMDEC(state)
      integer, INTENT(in) :: DIMDEC(ufluid) ! declare x,y,z dimensions of LS
      integer, INTENT(in) :: DIMDEC(usolid)
      integer, INTENT(in) :: DIMDEC(ughost)
      integer, INTENT(in) :: DIMDEC(history_dat)

        ! LS1,LS2,..,LSn,normal1,normal2,...normal_n 
        ! normal points from negative to positive
        !DIMV(LS)=x,y,z  
      !CP=Closest Point
      real(amrex_real), INTENT(in), target :: &
             LSCP(DIMV(LSCP),num_materials*(SDIM+1)) 
      real(amrex_real), pointer :: LSCP_ptr(D_DECL(:,:,:),:)

      ! FD=Finite Difference
      real(amrex_real), INTENT(in), target :: &
              LSFD(DIMV(LSFD),num_materials*SDIM)  
      real(amrex_real), pointer :: LSFD_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: state(DIMV(state),nden)
      real(amrex_real), pointer :: state_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: &
          ufluid(DIMV(ufluid),STATE_NCOMP_VEL+STATE_NCOMP_PRES) ! u,v,w,p
      real(amrex_real), pointer :: ufluid_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: &
              usolid(DIMV(usolid),nparts_ghost*SDIM) 
      real(amrex_real), pointer :: usolid_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(out),target :: &
              ughost(DIMV(ughost),nparts_ghost*SDIM) 
      real(amrex_real), pointer :: ughost_ptr(D_DECL(:,:,:),:)

       ! nhistory=nparts_ghost * (usolid_law_of_wall,uimage,usolid,angle)
      real(amrex_real), INTENT(out),target :: &
              history_dat(DIMV(history_dat),nhistory) 
      real(amrex_real), pointer :: history_dat_ptr(D_DECL(:,:,:),:)

      integer i,j,k
      integer ii,jj,kk
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten_solid(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_probe(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_fluid(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_MAC(-nhalf:nhalf,SDIM)
      real(amrex_real) LS_left(num_materials)
      real(amrex_real) LS_left_probe(num_materials)
      real(amrex_real) LS_right(num_materials)
      real(amrex_real) LS_right_probe(num_materials)
      integer side_solid,side_image
      integer partid
      integer im_solid
      integer im_fluid
      integer im_primary_left
      integer im_primary_left_probe
      integer im_primary_right
      integer im_primary_right_probe
      integer im
      integer dir
      integer isideSOLID,jsideSOLID,ksideSOLID
      integer isideFLUID,jsideFLUID,ksideFLUID
      integer iside_probe,jside_probe,kside_probe
      real(amrex_real), target :: n_raster(SDIM) ! points to solid
      real(amrex_real), target :: x_projection_raster(SDIM)
      real(amrex_real), target :: x_image_raster(SDIM)
      real(amrex_real), target :: x_probe_raster(SDIM)
      real(amrex_real) usolid_law_of_wall(SDIM)
      real(amrex_real) uimage_raster(SDIM)
      real(amrex_real) temperature_image
      real(amrex_real) temperature_wall
      real(amrex_real) temperature_wall_max
      real(amrex_real) dist_probe
      real(amrex_real) dist_fluid
      real(amrex_real), target :: usolid_raster(SDIM)
      real(amrex_real) angle_ACT_cell
      integer nhistory_sub
      type(law_of_wall_parm_type) :: law_of_wall_parm

      ncomp_sum_int_user12=ncomp_sum_int_user1+ncomp_sum_int_user2

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in wallfunction"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (ngrow_distance.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance invalid"
       stop
      endif
      if (ngrow_distance_in.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance_in invalid"
       stop
      endif
      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid"
       stop
      endif
      if ((nparts.ge.0).and.(nparts.le.num_materials)) then 
       ! do nothing
      else
       print *,"nparts invalid fort_wallfunction"
       stop
      endif
      if ((nparts_ghost.ge.1).and. &
          (nparts_ghost.le.num_materials).and. &
          (nparts_ghost.ge.nparts)) then 
       ! do nothing
      else
       print *,"nparts_ghost invalid fort_wallfunction"
       stop
      endif

      if ((nparts_ghost.eq.nparts).or.(nparts_ghost.eq.1)) then
       ! do nothing
      else
       print *,"nparts_ghost invalid"
       stop
      endif
      if (NS_sumdata_size.ne.IQ_TOTAL_SUM_COMP) then
       print *,"mismatch between NS_sumdata_size and IQ_TOTAL_SUM_COMP"
       stop
      endif

       ! ughost,imgVR,solVR,angle
      nhistory_sub=3*SDIM+1

      if (nhistory.eq.nparts_ghost*nhistory_sub) then
       ! do nothing
      else
       print *,"nhistory invalid"
       stop
      endif

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif 
      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif 
      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif 
      if (visc_coef.eq.fort_visc_coef) then
       ! do nothing
      else
       print *,"visc_coef.eq.fort_visc_coef is false"
       stop
      endif

      do im=1,num_materials

       if (abs(wall_model_velocity(im)).le.1.0D+20) then
        ! do nothing
       else
        print *,"wall_model_velocity(im) is corrupt"
        stop
       endif

       if ((law_of_the_wall(im).eq.0).or. &
           (law_of_the_wall(im).eq.1).or. &  ! CODY (turbulent)
           (law_of_the_wall(im).eq.2)) then  ! ZEYU (GNBC) 
        ! do nothing
       else
        print *,"law_of_the_wall invalid"
        stop
       endif
      enddo ! im=1..num_materials

      if ((data_dir.ge.0).and.(data_dir.le.SDIM-1)) then
       ! do nothing
      else
       print *,"data_dir invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (data_dir.eq.0) then
       ii=1
      else if (data_dir.eq.1) then
       jj=1
      else if ((data_dir.eq.SDIM-1).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"data_dir invalid"
       stop
      endif

       ! check that LS has ngrow_distance border cells.
       ! the "-1" means that LS is a cell centered variable
       ! instead of a face centered (staggared) variable.
       ! valid values for position type for staggared variables
       ! are 0,1,..,SDIM-1. 
      LSCP_ptr=>LSCP
      call checkbound_array(fablo,fabhi,LSCP_ptr,ngrow_distance,-1)
      LSFD_ptr=>LSFD
      call checkbound_array(fablo,fabhi,LSFD_ptr,ngrow_distance,-1)
      state_ptr=>state
      call checkbound_array(fablo,fabhi,state_ptr,ngrow_distance,-1)
      ufluid_ptr=>ufluid
      call checkbound_array(fablo,fabhi,ufluid_ptr,ngrow_distance,-1)
      usolid_ptr=>usolid
      call checkbound_array(fablo,fabhi,usolid_ptr,ngrow_distance,-1)
      ughost_ptr=>ughost
      call checkbound_array(fablo,fabhi,ughost_ptr,0,data_dir)
      history_dat_ptr=>history_dat
      call checkbound_array(fablo,fabhi,history_dat_ptr,0,data_dir)

      law_of_wall_parm%visc_coef=visc_coef
      law_of_wall_parm%time=time
      law_of_wall_parm%dt=dt
      law_of_wall_parm%level=level
      law_of_wall_parm%finest_level=finest_level
      law_of_wall_parm%bfact=bfact
      law_of_wall_parm%dx=>dx
      law_of_wall_parm%xlo=>xlo
      law_of_wall_parm%fablo=>fablo
      law_of_wall_parm%fabhi=>fabhi
      law_of_wall_parm%LSCP=>LSCP ! LS + gradient from closest point.
      law_of_wall_parm%LSFD=>LSFD ! gradient from finite differences
      law_of_wall_parm%state=>state ! density, temperature, species
      law_of_wall_parm%ufluid=>ufluid ! velocity and pressure
      law_of_wall_parm%usolid=>usolid ! solid velocity in solid bulk

      law_of_wall_parm%dxmin=dx(1)
      if (dx(2).lt.law_of_wall_parm%dxmin) then
       law_of_wall_parm%dxmin=dx(2)
      endif
      if (dx(SDIM).lt.law_of_wall_parm%dxmin) then
       law_of_wall_parm%dxmin=dx(SDIM)
      endif
      
       ! data_dir=0,1, or 2. 
      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi, &
              0,data_dir) 

       ! A FAB (fortran array box) is tessellated into tiles.
       ! i.e. a single FAB can contain multiple tiles.
       ! BOXLIB stores the FAB.   FAB's store data on a rectangular grid
       ! traverse faces of a given tile.
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       if (nparts.eq.0) then
         ! no solids, as placeholder put fictitious fluid velocity
         ! on the face.
        if (nparts_ghost.eq.1) then 
         do dir=1,SDIM
          ughost(D_DECL(i,j,k),dir)= &
            half*(ufluid(D_DECL(i,j,k),dir)+ &
                  ufluid(D_DECL(i-ii,j-jj,k-kk),dir))
         enddo
        else
         print *,"nparts_ghost invalid"
         stop
        endif
       else if ((nparts.ge.1).and.(nparts.le.num_materials)) then

        do partid=1,nparts

         do dir=1,SDIM
          usolid_raster(dir)= &
            half*(usolid(D_DECL(i,j,k),(partid-1)*SDIM+dir)+ &
                  usolid(D_DECL(i-ii,j-jj,k-kk),(partid-1)*SDIM+dir))
          ughost(D_DECL(i,j,k),(partid-1)*SDIM+dir)=usolid_raster(dir)
         enddo  

         do im=1,num_materials
          LS_right(im)=LSCP(D_DECL(i,j,k),im)
          LS_right_probe(im)=LSCP(D_DECL(i+ii,j+jj,k+kk),im)
          LS_left(im)=LSCP(D_DECL(i-ii,j-jj,k-kk),im)
          LS_left_probe(im)=LSCP(D_DECL(i-2*ii,j-2*jj,k-2*kk),im)
         enddo
         call get_primary_material(LS_right,im_primary_right)
         call get_primary_material(LS_right_probe,im_primary_right_probe)
         call get_primary_material(LS_left,im_primary_left)
         call get_primary_material(LS_left_probe,im_primary_left_probe)
         if ((im_primary_right.ge.1).and. &
             (im_primary_right.le.num_materials).and. &
             (im_primary_right_probe.ge.1).and. &
             (im_primary_right_probe.le.num_materials).and. &
             (im_primary_left.ge.1).and. &
             (im_primary_left.le.num_materials).and. &
             (im_primary_left_probe.ge.1).and. &
             (im_primary_left_probe.le.num_materials)) then
          ! do nothing
         else
          print *,"im_primary_left, im_primary_right, "
          print *,"im_primary_left_probe, or im_primary_right_probe "
          print *,"invalid"
          stop
         endif

         im_solid=im_solid_map(partid)+1  ! type integer: material id
         if ((im_solid.ge.1).and.(im_solid.le.num_materials)) then
          ! do nothing
         else
          print *,"im_solid invalid fort_wallfunction"
          stop
         endif

         side_solid=-1
         side_image=-1
         im_fluid=-1

         if (is_lag_part(im_solid).eq.1) then

          ! law of wall or dynamic contact angle treatment
          ! only for rigid substrates; not flexible substrates.
          if (is_rigid(im_solid).eq.1) then
           ! im_solid=material id of a rigid solid.
           ! Here, we test if cell center is in the solid.
           ! zero is defined in CONSTANTS.H
           ! CONSTANTS.H is defined in: ./BoxLib/Src/C_BaseLib/CONSTANTS.H

           do dir=1,SDIM
            n_raster(dir)=zero ! points to solid
           enddo

           if ((LS_right(im_solid).ge.zero).and. &
               (im_primary_right.eq.im_solid)) then

            side_solid=1 ! right side
            isideSOLID=i
            jsideSOLID=j
            ksideSOLID=k
            isideFLUID=i-ii
            jsideFLUID=j-jj
            ksideFLUID=k-kk
            iside_probe=i-2*ii
            jside_probe=j-2*jj
            kside_probe=k-2*kk

            if (1.eq.0) then
             print *,"sideSOLID ",side_solid,isideSOLID,jsideSOLID,ksideSOLID
             print *,"sideFLUID ",side_solid,isideFLUID,jsideFLUID,ksideFLUID
            endif

            if ((is_rigid(im_primary_left).eq.0).and. &
                (is_rigid(im_primary_left_probe).eq.0).and. &
                (im_primary_left.eq.im_primary_left_probe)) then
             side_image=0  ! left side
             im_fluid=im_primary_left
             n_raster(data_dir+1)=one ! points to solid
             do dir=1,SDIM
              uimage_raster(dir)= &
                ufluid(D_DECL(iside_probe,jside_probe,kside_probe),dir)
             enddo  
             dist_probe=LSCP(D_DECL(iside_probe,jside_probe,kside_probe), &
              im_fluid) 
             dist_fluid=LSCP(D_DECL(isideFLUID,jsideFLUID,ksideFLUID), &
              im_fluid) 
             temperature_image= &
                state(D_DECL(iside_probe,jside_probe,kside_probe), &
                     (im_primary_left-1)*num_state_material+ &
                     ENUM_TEMPERATUREVAR+1)
             temperature_wall= &
                state(D_DECL(isideSOLID,jsideSOLID,ksideSOLID), &
                     (im_solid-1)*num_state_material+ &
                     ENUM_TEMPERATUREVAR+1)
             temperature_wall_max= &
               NS_sumdata(IQ_MAXSTATE_SUM_COMP+2*(im_solid-1)+2)
            else if ((is_rigid(im_primary_left).eq.1).or. &
                     (is_rigid(im_primary_left_probe).eq.1).or. &
                     (im_primary_left.ne.im_primary_left_probe)) then
             ! do nothing
            else
             print *,"is_rigid(im_primary_left), or "
             print *,"is_rigid(im_primary_left_probe), "
             print *,"invalid"
             stop
            endif

           else if ((LS_left(im_solid).ge.zero).and. &
                    (im_primary_left.eq.im_solid)) then 
            side_solid=0  ! left side
            isideSOLID=i-ii
            jsideSOLID=j-jj
            ksideSOLID=k-kk
            isideFLUID=i
            jsideFLUID=j
            ksideFLUID=k
            iside_probe=i+ii
            jside_probe=j+jj
            kside_probe=k+kk

            if (1.eq.0) then
             print *,"sideSOLID ",side_solid,isideSOLID,jsideSOLID,ksideSOLID
             print *,"sideFLUID ",side_solid,isideFLUID,jsideFLUID,ksideFLUID
            endif

            if ((is_rigid(im_primary_right).eq.0).and. &
                (is_rigid(im_primary_right_probe).eq.0).and. &
                (im_primary_right.eq.im_primary_right_probe)) then
             side_image=1 ! right side
             im_fluid=im_primary_right
             n_raster(data_dir+1)=-one ! points to solid
             do dir=1,SDIM
              uimage_raster(dir)= &
                ufluid(D_DECL(iside_probe,jside_probe,kside_probe),dir)
             enddo  
             dist_probe=LSCP(D_DECL(iside_probe,jside_probe,kside_probe), &
              im_fluid) 
             dist_fluid=LSCP(D_DECL(isideFLUID,jsideFLUID,ksideFLUID), &
              im_fluid) 
             temperature_image= &
                state(D_DECL(iside_probe,jside_probe,kside_probe), &
                     (im_primary_right-1)*num_state_material+ &
                     ENUM_TEMPERATUREVAR+1)
             temperature_wall= &
                state(D_DECL(isideSOLID,jsideSOLID,ksideSOLID), &
                     (im_solid-1)*num_state_material+ &
                     ENUM_TEMPERATUREVAR+1)
             temperature_wall_max= &
               NS_sumdata(IQ_MAXSTATE_SUM_COMP+2*(im_solid-1)+2)
            else if ((is_rigid(im_primary_right).eq.1).or. &
                     (is_rigid(im_primary_right_probe).eq.1).or. &
                     (im_primary_right.ne.im_primary_right_probe)) then
             ! do nothing
            else
             print *,"is_rigid(im_primary_right), or "
             print *,"is_rigid(im_primary_right_probe), "
             print *,"invalid"
             stop
            endif

           else
            side_solid=-1
            side_image=-1
           endif

           if (((side_solid.eq.0).and.(side_image.eq.1)).or. &
               ((side_solid.eq.1).and.(side_image.eq.0))) then

            if ((im_fluid.ge.1).and.(im_fluid.le.num_materials)) then

             if (law_of_the_wall(im_fluid).eq.0) then
                     ! do nothing
             else if ((law_of_the_wall(im_fluid).eq.1).or. & !CODY(turbulent)
                      (law_of_the_wall(im_fluid).eq.2)) then !ZEYU(GNBC)

              !xsten(0,dir) gives dir'th component of coordinate of the storage
              !location of cell (i,j,k)
              !e.g. 1D:
              !      xsten(-2,1)  xsten(0,1)  xsten(2,1)
              !     |     .     |      .     |     .    |
              !          i-1           i          i+1
              ! xsten(-3,1)   xsten(-1,1) xsten(1,1)  xsten(3,1)
              call gridsten_level(xsten_solid, &
                isideSOLID,jsideSOLID,ksideSOLID, &
                level,nhalf)

              call gridsten_level(xsten_probe, &
                iside_probe,jside_probe,kside_probe, &
                level,nhalf)

              call gridsten_level(xsten_fluid, &
                isideFLUID,jsideFLUID,ksideFLUID, &
                level,nhalf)

              ! data_dir=0,1, or 2. 
              call gridstenMAC_level(xsten_MAC,i,j,k, &
                level,nhalf,data_dir)

              do dir=1,SDIM
               x_projection_raster(dir)=xsten_MAC(0,dir)
               x_image_raster(dir)=xsten_fluid(0,dir)
               x_probe_raster(dir)=xsten_probe(0,dir)
              enddo

              law_of_wall_parm%x_image_raster=>x_image_raster
              law_of_wall_parm%x_probe_raster=>x_probe_raster
              law_of_wall_parm%x_projection_raster=>x_projection_raster
              law_of_wall_parm%usolid_raster=>usolid_raster
              law_of_wall_parm%n_raster=>n_raster ! points to solid

              ! call CODY ESTEBEs routine here 
              ! (getGhostVel is declared in: GLOBALUTIL.F90)
              call getGhostVel( &
               law_of_wall_parm, & ! INTENT(in)
               law_of_the_wall(im_fluid), &
               isideSOLID, &
               jsideSOLID, &
               ksideSOLID, &
               isideFLUID, &
               jsideFLUID, &
               ksideFLUID, &
               iside_probe, &
               jside_probe, &
               kside_probe, &
               side_solid, &
               side_image, &
               data_dir, & ! data_dir=0,1, or 2
               uimage_raster, & ! INTENT(in)
               wall_model_velocity(im_fluid), &
               dist_probe, & ! INTENT(in)
               dist_fluid, & ! INTENT(in)
               temperature_image, & ! INTENT(in)
               temperature_wall, & ! INTENT(in)
               temperature_wall_max, & ! INTENT(in)
               usolid_law_of_wall, & ! INTENT(out)
               angle_ACT_cell, & !INTENT(out) dyn. contact angle at image point
               im_fluid, &
               im_solid)

               ! solid "ghost" velocity in the solid regions.
              do dir=1,SDIM
               !ughost is an output.
               ughost(D_DECL(i,j,k),(partid-1)*SDIM+dir)= &
                usolid_law_of_wall(dir)

               history_dat(D_DECL(i,j,k),(partid-1)*nhistory_sub+dir)= &
                usolid_law_of_wall(dir)
               history_dat(D_DECL(i,j,k),(partid-1)*nhistory_sub+SDIM+dir)= &
                uimage_raster(dir)
               history_dat(D_DECL(i,j,k),(partid-1)*nhistory_sub+2*SDIM+dir)= &
                usolid_raster(dir)
              enddo  ! dir=1..sdim

              history_dat(D_DECL(i,j,k), &
                (partid-1)*nhistory_sub+nhistory_sub)=angle_ACT_cell
             else
              print *,"law_of_the_wall(im_fluid) invalid"
              stop
             endif
            else
             print *,"im_fluid invalid"
             stop
            endif
           
           else if ((side_solid.eq.-1).or.(side_image.eq.-1)) then
            ! do nothing
           else
            print *,"side_solid or side_image invalid"
            stop
           endif

          else if (is_rigid(im_solid).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(im_solid) invalid"
           stop
          endif
         else 
          print *,"is_lag_part(im_solid) invalid"
          stop
         endif
        enddo ! partid=1..nparts
       else
        print *,"nparts invalid"
        stop
       endif
      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine fort_wallfunction


       ! called from:NavierStokes::init_FSI_GHOST_MAC_MF_predict(int ngrow) 
       ! (in NavierStokes.cpp)
       ! if nparts==0 => interpolate state cell velocity to MAC grid.
       ! if nparts>0 => interpolate solid cell velocity
       ! to MAC grid.
      subroutine fort_wallfunction_predict( &
       data_dir, &
       im_solid_map, &
       level, &
       finest_level, &
       ngrow_distance_in, &
       nparts, &
       nparts_ghost, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx, &
       dt, &
       time, &
       ufluid,DIMS(ufluid), &
       usolid,DIMS(usolid), &
       ughost,DIMS(ughost)) &
      bind(c,name='fort_wallfunction_predict')
      use probf90_module
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: data_dir
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: ngrow_distance_in
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_ghost
      integer, INTENT(in) :: im_solid_map(nparts_ghost)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in), target :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in), target :: xlo(SDIM)
      real(amrex_real), INTENT(in), target :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: DIMDEC(ufluid) ! declare x,y,z dimensions of LS
      integer, INTENT(in) :: DIMDEC(usolid)
      integer, INTENT(in) :: DIMDEC(ughost)

      real(amrex_real), INTENT(in), target :: &
           ufluid(DIMV(ufluid),STATE_NCOMP_VEL+STATE_NCOMP_PRES) ! u,v,w,p
      real(amrex_real), pointer :: ufluid_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: usolid(DIMV(usolid),nparts_ghost*SDIM) 
      real(amrex_real), pointer :: usolid_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(out),target :: ughost(DIMV(ughost),nparts_ghost*SDIM) 
      real(amrex_real), pointer :: ughost_ptr(D_DECL(:,:,:),:)

      integer i,j,k
      integer ii,jj,kk
      integer dir
      integer partid
      real(amrex_real) :: usolid_raster(SDIM)


      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in wallfunction_predict"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (ngrow_distance_in.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance_in invalid"
       stop
      endif
      if (ngrow_distance.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance invalid"
       stop
      endif

      if ((nparts.ge.0).and.(nparts.le.num_materials)) then 
       ! do nothing
      else
       print *,"nparts invalid fort_wallfunction_predict"
       stop
      endif
      if ((nparts_ghost.ge.1).and. &
          (nparts_ghost.le.num_materials).and. &
          (nparts_ghost.ge.nparts)) then 
       ! do nothing
      else
       print *,"nparts_ghost invalid fort_wallfunction"
       stop
      endif

      if ((nparts_ghost.eq.nparts).or.(nparts_ghost.eq.1)) then
       ! do nothing
      else
       print *,"nparts_ghost invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif 
      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif 

      if ((data_dir.ge.0).and.(data_dir.le.SDIM-1)) then
       ! do nothing
      else
       print *,"data_dir invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (data_dir.eq.0) then
       ii=1
      else if (data_dir.eq.1) then
       jj=1
      else if ((data_dir.eq.SDIM-1).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"data_dir invalid"
       stop
      endif

      ufluid_ptr=>ufluid
      call checkbound_array(fablo,fabhi,ufluid_ptr,ngrow_distance,-1)
      usolid_ptr=>usolid
      call checkbound_array(fablo,fabhi,usolid_ptr,ngrow_distance,-1)
      ughost_ptr=>ughost
      call checkbound_array(fablo,fabhi,ughost_ptr,0,data_dir)

       ! data_dir=0,1, or 2. 
      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi, &
              0,data_dir) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       if (nparts.eq.0) then
         ! no solids, as placeholder put fictitious fluid velocity
         ! on the face.
        if (nparts_ghost.eq.1) then 
         do dir=1,SDIM
          ughost(D_DECL(i,j,k),dir)= &
            half*(ufluid(D_DECL(i,j,k),dir)+ &
                  ufluid(D_DECL(i-ii,j-jj,k-kk),dir))
         enddo
        else
         print *,"nparts_ghost invalid"
         stop
        endif
       else if ((nparts.ge.1).and.(nparts.le.num_materials)) then

        do partid=1,nparts

         do dir=1,SDIM
          usolid_raster(dir)= &
            half*(usolid(D_DECL(i,j,k),(partid-1)*SDIM+dir)+ &
                  usolid(D_DECL(i-ii,j-jj,k-kk),(partid-1)*SDIM+dir))
          ughost(D_DECL(i,j,k),(partid-1)*SDIM+dir)=usolid_raster(dir)
         enddo  

        enddo ! partid=1..nparts
       else
        print *,"nparts invalid"
        stop
       endif
      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine fort_wallfunction_predict

      ! tag = 1 -> donor cell
      ! tag = 2 -> receving cell
      ! tag = 0 -> none of above
      subroutine fort_tagexpansion(&
       nden, &
       freezing_model, &
       distribute_from_target, &
       time, &
       vofbc, &
       expect_mdot_sign, &
       mdot_sum, &
       mdot_sum_comp, &
       im_source, & !intent(in)
       im_dest, & !intent(in)
       indexEXP, & !intent(in) 0<=indexEXP<2*num_interfaces
       level,finest_level, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx,dt, &
       maskcov,DIMS(maskcov), &
       tag, &
       DIMS(tag), &
       tag_comp, &
       DIMS(tag_comp), &
       expan,DIMS(expan), &
       expan_comp,DIMS(expan_comp), &
       denstate,DIMS(denstate), &
       LS,DIMS(LS), &  ! newdistfab=(*localMF[LSNEW_MF])[mfi]
       recon,DIMS(recon)) &
      bind(c,name='fort_tagexpansion')

      use probf90_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: nden
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: mdot_sum
      real(amrex_real), INTENT(inout) :: mdot_sum_comp
      real(amrex_real), INTENT(in) :: expect_mdot_sign
      integer, INTENT(in) :: im_source,im_dest
      integer :: im_ice
      integer, INTENT(in) :: indexEXP
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(tag)
      integer, INTENT(in) :: DIMDEC(tag_comp)
      integer, INTENT(in) :: DIMDEC(expan)
      integer, INTENT(in) :: DIMDEC(expan_comp)
      integer, INTENT(in) :: DIMDEC(denstate)
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(recon)
      real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out), target :: tag(DIMV(tag))
      real(amrex_real), pointer :: tag_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out), target :: tag_comp(DIMV(tag_comp))
      real(amrex_real), pointer :: tag_comp_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: expan(DIMV(expan),2*num_interfaces)
      real(amrex_real), pointer :: expan_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
           expan_comp(DIMV(expan_comp),2*num_interfaces)
      real(amrex_real), pointer :: expan_comp_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: denstate(DIMV(denstate),nden)
      real(amrex_real), pointer :: denstate_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: LS(DIMV(LS),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
           recon(DIMV(recon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: recon_ptr(D_DECL(:,:,:),:)

      integer local_freezing_model
      integer vofbc(SDIM,2)
      integer i,j,k
      real(amrex_real) VFRAC(num_materials)
      real(amrex_real) VDOT
      real(amrex_real) local_denstate(nden)
      real(amrex_real) LSCELL(num_materials)
      real(amrex_real) ICEMASK
      real(amrex_real) icefacecut
      integer im,im_opp
      integer ireverse
      integer iten
      integer im_primary
      integer im_primary_icemask
      integer vofcomp
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_center(SDIM)
      integer local_mask
      integer dir
      integer index_compare
      integer complement_flag
      real(amrex_real) LL

      tag_ptr=>tag
      tag_comp_ptr=>tag_comp

      denstate_ptr=>denstate

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      !! Sanity checks

      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif

      if ((im_source.ge.1).and.(im_source.le.num_materials)) then
       ! do nothing
      else
       print *,"im_source invalid"
       stop
      endif
      if ((im_dest.ge.1).and.(im_dest.le.num_materials)) then
       ! do nothing
      else
       print *,"im_dest invalid"
       stop
      endif
      if (im_dest.eq.im_source) then
       print *,"im_dest or im_source invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in distribute_expansion"
       stop
      endif
      if ((indexEXP.lt.0).or.(indexEXP.ge.2*num_interfaces)) then
       print *,"indexEXP invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      if (ngrow_distance.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance invalid"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nden.eq.num_materials*num_state_material) then
       ! do nothing
      else
       print *,"nden invalid"
       stop
      endif

      maskcov_ptr=>maskcov
      LS_ptr=>LS
      recon_ptr=>recon
      expan_ptr=>expan
      expan_comp_ptr=>expan_comp
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,denstate_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      call checkbound_array(fablo,fabhi,recon_ptr,1,-1)
      call checkbound_array(fablo,fabhi,expan_ptr,ngrow_distance,-1)
      call checkbound_array(fablo,fabhi,expan_comp_ptr,ngrow_distance,-1)
      call checkbound_array1(fablo,fabhi,tag_ptr,ngrow_distance,-1)
      call checkbound_array1(fablo,fabhi,tag_comp_ptr,ngrow_distance,-1)

      local_freezing_model=freezing_model(indexEXP+1)
      if (is_valid_freezing_modelF(local_freezing_model).eq.1) then
       ! do nothing
      else
       print *,"local_freezing_model invalid 17"
       stop
      endif
      if ((distribute_from_target(indexEXP+1).lt.0).or. &
          (distribute_from_target(indexEXP+1).gt.1)) then
       print *,"distribute_from_target invalid"
       stop
      endif
      LL=get_user_latent_heat(indexEXP+1,room_temperature,1)
      if (LL.eq.zero) then
       ! check nothing
      else if (LL.ne.zero) then
       im_ice=0
       if (is_ice(im_dest).eq.1) then
        im_ice=im_dest
       else if (is_ice(im_source).eq.1) then
        im_ice=im_source
       else if ((is_ice(im_dest).eq.0).and. &
                (is_ice(im_source).eq.0)) then
        im_ice=0
       else
        print *,"is_ice invalid"
        stop
       endif
       if (im_ice.eq.0) then
        ! do nothing
       else if (im_ice.eq.im_dest) then ! freezing
        if (distribute_from_target(indexEXP+1).eq.1) then
         ! distribute_from_target=true
         ! distribute_to_target=false
         ! source=melt  dest (target)=ice
        else 
         print *,"distribute_from_target should be 1(freezing)"
         stop
        endif
       else if (im_ice.eq.im_source) then ! melting
        if (distribute_from_target(indexEXP+1).eq.0) then
         ! distribute_from_target=false
         ! distribute_to_target=true
         ! source=ice  dest (target)=melt
        else 
         print *,"distribute_from_target should be 0(melting)"
         stop
        endif
       else
        print *,"im_ice invalid"
        stop
       endif
      else
       print *,"LL invalid: ",LL
       stop
      endif

      ! Iterate over the box
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        call gridsten_level(xsten,i,j,k,level,nhalf)
        do dir=1,SDIM
         xsten_center(dir)=xsten(0,dir)
        enddo

        !! water freezing to ice
        !! im_source -> index for water material
        !! im_dest   -> index for ice material
            
         ! check for being a donor cell:
         ! 1. non-zero expansion term
         ! 2. (F_ice < 0.5) or
         !    (icemask=0.0)
         ! 
         ! check for being a receiving cell:
         ! (F_ice_tessellate > 0.5)&&(icemask>0.0)
         ! 
         ! weight for redistribution: equal weight for all cells in the
         ! redistribution stencil.
              
        tag(D_DECL(i,j,k)) = zero
        tag_comp(D_DECL(i,j,k)) = zero

        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         VFRAC(im)=recon(D_DECL(i,j,k),vofcomp)
         LSCELL(im)=LS(D_DECL(i,j,k),im)
        enddo
        do im=1,nden
         local_denstate(im)=denstate(D_DECL(i,j,k),im)
        enddo

         ! "get_primary_material"
         ! first checks the rigid materials for a positive LS; if none
         ! exist, then "get_primary_material" checks the fluid materials.
        call get_primary_material(LSCELL,im_primary)

        VDOT=expan(D_DECL(i,j,k),indexEXP+1)
        if (expect_mdot_sign.eq.one) then
         if (VDOT.lt.zero) then
          print *,"expecting VDOT>=0"
          print *,"sign+ VDOT invalid i,j,k,vdot ",i,j,k,VDOT
          stop
         else if (VDOT.ge.zero) then
          ! do nothing
         else
          print *,"VDOT bust"
          stop
         endif
        else if (expect_mdot_sign.eq.-one) then
         if (VDOT.gt.zero) then
          print *,"expecting VDOT<=0"
          print *,"sign- VDOT invalid i,j,k,vdot ",i,j,k,VDOT
          stop
         else if (VDOT.le.zero) then
          ! do nothing
         else
          print *,"VDOT bust"
          stop
         endif
        else
         print *,"expect_mdot_sign invalid: ",expect_mdot_sign
         stop
        endif

        do complement_flag=0,1

         im_ice=0
         if (is_ice(im_dest).eq.1) then
          im_ice=im_dest
         else if (is_ice(im_source).eq.1) then
          im_ice=im_source
         else if ((is_ice(im_dest).eq.0).and. &
                  (is_ice(im_source).eq.0)) then
          im_ice=0
         else
          print *,"is_ice invalid"
          stop
         endif

         if (im_ice.eq.0) then!both source and dest can be distributed to.
          ICEMASK=one
         else if ((im_ice.ge.1).and. &
                  (im_ice.le.num_materials)) then

          ! in: fort_tagexpansion
          ! ICEMASK=0 => mask off this cell.
          ! ICEMASK=1 => do nothing
          ! get_icemask_and_icefacecut declared in PROB.F90
          call get_icemask_and_icefacecut( &
           nden, &
           xsten_center, &
           time, &
           dx,bfact, &
           ICEMASK, &
           icefacecut, &
           im, &
           im_opp, &
           im_primary_icemask, &
           ireverse, &
           local_denstate, &
           LSCELL, &
           VFRAC, &
           distribute_from_target, &
           complement_flag)

          if (ireverse.eq.-1) then!both source and dest can be distributed to.
           ICEMASK=one
          else if ((ireverse.eq.0).or.(ireverse.eq.1)) then
           call get_iten(im,im_opp,iten)
           index_compare=iten+ireverse*num_interfaces-1
           if ((index_compare.ge.0).and. &
               (index_compare.lt.2*num_interfaces)) then
            if (index_compare.eq.indexEXP) then
             ! do nothing
            else
             ICEMASK=one !assume both source and dest can be distributed to.
            endif
           else
            print *,"index_compare invalid"
            stop
           endif
          else
           print *,"ireverse invalid"
           stop
          endif
 
         else 
          print *,"im_ice invalid"
          stop
         endif

         if (complement_flag.eq.0) then
          mdot_sum=mdot_sum+VDOT
         else if (complement_flag.eq.1) then
          mdot_sum_comp=mdot_sum_comp+VDOT
         else
          print *,"complement_flag invalid"
          stop
         endif

         if ((is_rigid(im_primary).eq.0).and. &
             (is_FSI_rigid(im_primary).eq.0)) then

           ! first tag donor cells (tag=one)
          if (VDOT.ne.zero) then ! nonzero source

           if ( &
             swap1_0(distribute_from_target(indexEXP+1),complement_flag) &
             .eq.0) then

            if ((VFRAC(im_dest).lt.half).or. &
                (ICEMASK.eq.zero)) then
             if (complement_flag.eq.0) then
              tag(D_DECL(i,j,k)) = one ! donor cell
             else if (complement_flag.eq.1) then
              tag_comp(D_DECL(i,j,k)) = one ! donor cell
             else
              print *,"complement_flag invalid"
              stop
             endif
            else if ((VFRAC(im_dest).ge.half).and. &
                     (ICEMASK.eq.one)) then
             ! do nothing - acceptor cell
            else
             print *,"VFRAC or ICEMASK bust"
             stop      
            endif

           else if ( &
             swap1_0(distribute_from_target(indexEXP+1),complement_flag) &
             .eq.1) then

            if ((VFRAC(im_source).lt.half).or. &
                (ICEMASK.eq.zero)) then
             if (complement_flag.eq.0) then
              tag(D_DECL(i,j,k)) = one ! donor cell
             else if (complement_flag.eq.1) then
              tag_comp(D_DECL(i,j,k)) = one ! donor cell
             else
              print *,"complement_flag invalid"
              stop
             endif
            else if ((VFRAC(im_source).ge.half).and. &
                     (ICEMASK.eq.one)) then
             ! do nothing - acceptor cell
            else
             print *,"VFRAC or ICEMASK bust"     
             stop
            endif

           else
            print *,"distribute_from_target(indexEXP+1) invalid"
            stop
           endif

          else if (VDOT.eq.zero) then
           ! do nothing
          else
           print *,"VDOT became corrupt"
           stop
          endif 

           ! now we tag receiver cells (tag=two)
          if ( &
            swap1_0(distribute_from_target(indexEXP+1),complement_flag) &
            .eq.0) then

           if ((VFRAC(im_dest).ge.half).and. &
               (ICEMASK.eq.one)) then
            if (complement_flag.eq.0) then
             tag(D_DECL(i,j,k)) = two ! receiver
            else if (complement_flag.eq.1) then
             tag_comp(D_DECL(i,j,k)) = two ! receiver
            else
             print *,"complement_flag invalid"
             stop
            endif
           else if ((VFRAC(im_dest).lt.half).or. &
                    (ICEMASK.eq.zero)) then
            ! do nothing - donor cell if VDOT<>0
           else
            print *,"VFRAC or ICEMASK bust"
            stop      
           endif
 
          else if ( &
              swap1_0(distribute_from_target(indexEXP+1),complement_flag) &
              .eq.1) then

           if ((VFRAC(im_source).ge.half).and. &
               (ICEMASK.eq.one)) then
            if (complement_flag.eq.0) then
             tag(D_DECL(i,j,k)) = two ! receiver
            else if (complement_flag.eq.1) then
             tag_comp(D_DECL(i,j,k)) = two ! receiver
            else
             print *,"complement_flag invalid"
             stop
            endif
           else if ((VFRAC(im_source).lt.half).or. &
                    (ICEMASK.eq.zero)) then
            ! do nothing - donor cell if VDOT<>0
           else
            print *,"VFRAC or ICEMASK bust"     
            stop
           endif
 
          else
           print *,"distribute_from_target(indexEXP+1) invalid"
           stop
          endif

         !in the prescribed solid.
         else if ((is_rigid(im_primary).eq.1).or. &
                  (is_FSI_rigid(im_primary).eq.1)) then 
          ! do nothing (tag initialized to 0, neither donor nor receiver)
         else
          print *,"is_rigid(im_primary) invalid or"
          print *,"is_FSI_rigid(im_primary) invalid"
          stop
         endif 

        enddo ! complement_flag=0..1

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine fort_tagexpansion

      subroutine redistribute_weight(xmain,xside,crit_weight)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: xmain(SDIM)
      real(amrex_real), INTENT(in) :: xside(SDIM)
      real(amrex_real), INTENT(out) :: crit_weight

      crit_weight=one

      return
      end subroutine redistribute_weight


      ! recon( num_materials * ngeom_recon )
      ! ngeom_recon=2 * SDIM + 3
      ! volume fraction, reference centroid, order, slope, intercept
      ! num_interfaces = 
      ! (num_materials * (num_materials-1))/2  -> number of possible surface 
      !                                 tension coefficients
      ! if num_materials=2, num_interfaces=1
      ! if num_materials=3, num_interfaces=3    12 13 23
      ! if num_materials=4, num_interfaces=6    12 13 14 23 24 34
      
      ! This will be called before fort_initjumpterm and after
      ! fort_convertmaterial
      ! tag = 1 -> donor cell
      ! tag = 2 -> receving cell
      ! tag = 0 -> non of above
      subroutine fort_distributeexpansion(&
       mdot_sum, &
       mdot_lost, &
       mdot_sum_comp, &
       mdot_lost_comp, &
       im_source, &
       im_dest, &
       indexEXP, &
       level,finest_level, &
       domlo,domhi, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx,dt, &
       maskcov,DIMS(maskcov),&
       LS,DIMS(LS),&
       tag, &
       DIMS(tag),&
       tag_comp, &
       DIMS(tag_comp),&
       weightfab, &
       DIMS(weightfab),&
       weight_comp, &
       DIMS(weight_comp),&
       expan, &
       DIMS(expan), &
       expan_comp, &
       DIMS(expan_comp) ) &
       bind(c,name='fort_distributeexpansion')

       use probf90_module
       use global_utility_module
       use geometry_intersect_module

       IMPLICIT NONE

       real(amrex_real), INTENT(inout) :: mdot_sum,mdot_lost
       real(amrex_real), INTENT(inout) :: mdot_sum_comp,mdot_lost_comp
       integer, INTENT(in) :: im_source,im_dest,indexEXP
       integer, INTENT(in) :: level,finest_level
       integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer :: growlo(3),growhi(3)
       integer :: stenlo(3),stenhi(3)
       integer, INTENT(in) :: bfact
       real(amrex_real), INTENT(in) :: xlo(SDIM)
       real(amrex_real), INTENT(in) :: dx(SDIM)
       real(amrex_real), INTENT(in) :: dt
       integer, INTENT(in) :: DIMDEC(maskcov)
       integer, INTENT(in) :: DIMDEC(LS)
       integer, INTENT(in) :: DIMDEC(tag)
       integer, INTENT(in) :: DIMDEC(tag_comp)
       integer, INTENT(in) :: DIMDEC(weightfab)
       integer, INTENT(in) :: DIMDEC(weight_comp)
       integer, INTENT(in) :: DIMDEC(expan)
       integer, INTENT(in) :: DIMDEC(expan_comp)
       real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov))
       real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in), target :: LS(DIMV(LS),num_materials)
       real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in), target :: tag(DIMV(tag))
       real(amrex_real), pointer :: tag_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in), target :: tag_comp(DIMV(tag_comp))
       real(amrex_real), pointer :: tag_comp_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in), target :: weightfab(DIMV(weightfab))
       real(amrex_real), pointer :: weightfab_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in), target :: weight_comp(DIMV(weight_comp))
       real(amrex_real), pointer :: weight_comp_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(inout), target :: expan(DIMV(expan),2*num_interfaces)
       real(amrex_real), pointer :: expan_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(inout), target ::  &
         expan_comp(DIMV(expan_comp),2*num_interfaces)
       real(amrex_real), pointer :: expan_comp_ptr(D_DECL(:,:,:),:)

       real(amrex_real), dimension(D_DECL(:,:,:),:), allocatable,target :: expan_new
       real(amrex_real), pointer :: expan_new_ptr(D_DECL(:,:,:),:)
       real(amrex_real) local_expan_new
       real(amrex_real) local_expan_old

       integer i,j,k
       integer dir
       integer i_n,j_n,k_n
       real(amrex_real) total_weight,crit_weight
       real(amrex_real) TAGLOC,TAGSIDE
       integer, PARAMETER :: nhalf=1
       real(amrex_real) xsten(-nhalf:nhalf,SDIM)
       real(amrex_real) xsten_n(-nhalf:nhalf,SDIM)
       real(amrex_real) xmain(SDIM)
       real(amrex_real) xside(SDIM)
       integer local_mask

       expan_ptr=>expan
       expan_comp_ptr=>expan_comp

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       !! Sanity checks

       if ((im_source.ge.1).and.(im_source.le.num_materials)) then
        ! do nothing
       else
        print *,"im_source invalid"
        stop
       endif
       if ((im_dest.ge.1).and.(im_dest.le.num_materials)) then
        ! do nothing
       else
        print *,"im_dest invalid"
        stop
       endif
       if (im_dest.eq.im_source) then
        print *,"im_dest or im_source invalid"
        stop
       endif

       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid in distribute_expansion"
        stop
       end if
       if ((indexEXP.lt.0).or.(indexEXP.ge.2*num_interfaces)) then
        print *,"indexEXP invalid"
        stop
       endif
       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if (ngrow_distance.eq.4) then
        ! do nothing
       else
        print *,"ngrow_distance invalid"
        stop
       endif

       allocate(expan_new(DIMV(expan),2))
       expan_new_ptr=>expan_new

       maskcov_ptr=>maskcov
       LS_ptr=>LS
       tag_ptr=>tag
       tag_comp_ptr=>tag_comp
       weightfab_ptr=>weightfab
       weight_comp_ptr=>weight_comp

       call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
       call checkbound_array(fablo,fabhi,LS_ptr,ngrow_distance,-1)

       call checkbound_array(fablo,fabhi,expan_ptr,ngrow_distance,-1)
       call checkbound_array(fablo,fabhi,expan_comp_ptr,ngrow_distance,-1)
       call checkbound_array(fablo,fabhi,expan_new_ptr,ngrow_distance,-1)

       call checkbound_array1(fablo,fabhi,tag_ptr,ngrow_distance,-1)
       call checkbound_array1(fablo,fabhi,tag_comp_ptr,ngrow_distance,-1)
       call checkbound_array1(fablo,fabhi,weightfab_ptr,ngrow_distance,-1)
       call checkbound_array1(fablo,fabhi,weight_comp_ptr,ngrow_distance,-1)

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        expan_new(D_DECL(i,j,k),1)=zero
        expan_new(D_DECL(i,j,k),2)=zero
       enddo
       enddo
       enddo

       ! Iterate over the box
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        local_mask=NINT(maskcov(D_DECL(i,j,k)))

        if (local_mask.eq.1) then

         ! if a receiving cell
         TAGLOC=tag(D_DECL(i,j,k))
         if(TAGLOC.eq.two) then ! receiver
          call stencilbox(i,j,k,fablo,fabhi,stenlo,stenhi,ngrow_distance)

          do dir=1,SDIM
           if (stenlo(dir).lt.domlo(dir)) then
            stenlo(dir)=domlo(dir)
           endif
           if (stenhi(dir).gt.domhi(dir)) then
            stenhi(dir)=domhi(dir)
           endif
          enddo

          call gridsten_level(xsten,i,j,k,level,nhalf)
          do dir=1,SDIM
           xmain(dir)=xsten(0,dir)
          enddo
            
          do k_n=stenlo(3),stenhi(3)
          do j_n=stenlo(2),stenhi(2)
          do i_n=stenlo(1),stenhi(1)

           TAGSIDE=tag(D_DECL(i_n,j_n,k_n))
           if ((TAGSIDE.eq.one).or. & ! doner
               (TAGSIDE.eq.two)) then ! receiver

            call gridsten_level(xsten_n,i_n,j_n,k_n,level,nhalf)
            do dir=1,SDIM
             xside(dir)=xsten_n(0,dir)
            enddo
            call redistribute_weight(xmain,xside,crit_weight)
            total_weight=weightfab(D_DECL(i_n,j_n,k_n))
            if (total_weight.gt.zero) then
             call redistribute_weight(xmain,xside,crit_weight)
             if (crit_weight.le.total_weight) then
              expan_new(D_DECL(i,j,k),1) = &
               expan_new(D_DECL(i,j,k),1) + &
               expan(D_DECL(i_n,j_n,k_n),indexEXP+1)*crit_weight/total_weight
             else
              print *,"crit_weight invalid"
              stop
             endif
            else
             print *,"total_weight invalid"
             stop
            endif

           else if (TAGSIDE.eq.zero) then
            ! do nothing
           else
            print *,"TAGSIDE invalid"
            stop
           endif 

          end do ! k_n
          end do ! j_n
          end do ! i_n

         else if ((TAGLOC.eq.one).or. & ! doner
                  (TAGLOC.eq.zero)) then
          ! do nothing
         else
          print *,"TAGLOC invalid"
          stop
         endif 

          ! ---------------- DISTRIBUTE FOR COMPLEMENT ----------------

         ! if a receiving cell
         TAGLOC=tag_comp(D_DECL(i,j,k))
         if(TAGLOC.eq.two) then ! receiver
          call stencilbox(i,j,k,fablo,fabhi,stenlo,stenhi,ngrow_distance)

          do dir=1,SDIM
           if (stenlo(dir).lt.domlo(dir)) then
            stenlo(dir)=domlo(dir)
           endif
           if (stenhi(dir).gt.domhi(dir)) then
            stenhi(dir)=domhi(dir)
           endif
          enddo

          call gridsten_level(xsten,i,j,k,level,nhalf)
          do dir=1,SDIM
           xmain(dir)=xsten(0,dir)
          enddo
            
          do k_n=stenlo(3),stenhi(3)
          do j_n=stenlo(2),stenhi(2)
          do i_n=stenlo(1),stenhi(1)

           TAGSIDE=tag_comp(D_DECL(i_n,j_n,k_n))
           if ((TAGSIDE.eq.one).or. & ! doner
               (TAGSIDE.eq.two)) then ! receiver

            call gridsten_level(xsten_n,i_n,j_n,k_n,level,nhalf)
            do dir=1,SDIM
             xside(dir)=xsten_n(0,dir)
            enddo
            call redistribute_weight(xmain,xside,crit_weight)
            total_weight=weight_comp(D_DECL(i_n,j_n,k_n))
            if (total_weight.gt.zero) then
             call redistribute_weight(xmain,xside,crit_weight)
             if (crit_weight.le.total_weight) then
              expan_new(D_DECL(i,j,k),2) = &
               expan_new(D_DECL(i,j,k),2) + &
               expan_comp(D_DECL(i_n,j_n,k_n),indexEXP+1)* &
               crit_weight/total_weight
             else
              print *,"crit_weight invalid"
              stop
             endif
            else
             print *,"total_weight invalid"
             stop
            endif

           else if (TAGSIDE.eq.zero) then
            ! do nothing
           else
            print *,"TAGSIDE invalid"
            stop
           endif 

          end do ! k_n
          end do ! j_n
          end do ! i_n

         else if ((TAGLOC.eq.one).or. & ! doner
                  (TAGLOC.eq.zero)) then
          ! do nothing
         else
          print *,"TAGLOC invalid"
          stop
         endif 

        else if (local_mask.eq.0) then
         ! do nothing
        else
         print *,"local_mask invalid"
         stop
        endif

       end do ! k
       end do ! j
       end do ! i

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        local_mask=NINT(maskcov(D_DECL(i,j,k)))
        if (local_mask.eq.1) then
         local_expan_new=expan_new(D_DECL(i,j,k),1)
         local_expan_old=expan(D_DECL(i,j,k),indexEXP+1)
         expan(D_DECL(i,j,k),indexEXP+1)=local_expan_new
         mdot_sum=mdot_sum+local_expan_new
         mdot_lost=mdot_lost+local_expan_old-local_expan_new

         local_expan_new=expan_new(D_DECL(i,j,k),2)
         local_expan_old=expan_comp(D_DECL(i,j,k),indexEXP+1)
         expan_comp(D_DECL(i,j,k),indexEXP+1)=local_expan_new
         mdot_sum_comp=mdot_sum_comp+local_expan_new
         mdot_lost_comp=mdot_lost_comp+local_expan_old-local_expan_new

        else if (local_mask.eq.0) then
         ! do nothing
        else
         print *,"local_mask invalid"
         stop
        endif
       end do ! k
       end do ! j
       end do ! i

       deallocate(expan_new)

      end subroutine fort_distributeexpansion

      subroutine fort_accept_weight(&
       im_source, &
       im_dest, &
       indexEXP, &
       level,finest_level, &
       domlo,domhi, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx,dt, &
       maskcov,DIMS(maskcov),&
       LS,DIMS(LS),&
       tag, &
       DIMS(tag),&
       tag_comp, &
       DIMS(tag_comp),&
       weightfab, &
       DIMS(weightfab),&
       weight_comp, &
       DIMS(weight_comp),&
       expan, &
       DIMS(expan), &
       expan_comp, &
       DIMS(expan_comp) ) &
       bind(c,name='fort_accept_weight')

       use probf90_module
       use global_utility_module
       use geometry_intersect_module

       IMPLICIT NONE

       integer, INTENT(in) :: im_source,im_dest,indexEXP
       integer, INTENT(in) :: level,finest_level
       integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer :: growlo(3),growhi(3)
       integer :: stenlo(3),stenhi(3)
       integer, INTENT(in) :: bfact
       real(amrex_real), INTENT(in) :: xlo(SDIM)
       real(amrex_real), INTENT(in) :: dx(SDIM)
       real(amrex_real), INTENT(in) :: dt
       integer, INTENT(in) :: DIMDEC(maskcov)
       integer, INTENT(in) :: DIMDEC(LS)
       integer, INTENT(in) :: DIMDEC(tag)
       integer, INTENT(in) :: DIMDEC(tag_comp)
       integer, INTENT(in) :: DIMDEC(weightfab)
       integer, INTENT(in) :: DIMDEC(weight_comp)
       integer, INTENT(in) :: DIMDEC(expan)
       integer, INTENT(in) :: DIMDEC(expan_comp)
       real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov))
       real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in), target :: LS(DIMV(LS),num_materials)
       real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in), target :: tag(DIMV(tag))
       real(amrex_real), pointer :: tag_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in), target :: tag_comp(DIMV(tag_comp))
       real(amrex_real), pointer :: tag_comp_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(out), target :: weightfab(DIMV(weightfab))
       real(amrex_real), pointer :: weightfab_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(out), target :: weight_comp(DIMV(weight_comp))
       real(amrex_real), pointer :: weight_comp_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(inout), target :: expan(DIMV(expan),2*num_interfaces)
       real(amrex_real), pointer :: expan_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(inout), target ::  &
         expan_comp(DIMV(expan_comp),2*num_interfaces)
       real(amrex_real), pointer :: expan_comp_ptr(D_DECL(:,:,:),:)

       integer i,j,k
       integer dir
       integer i_n,j_n,k_n
       real(amrex_real) total_weight,crit_weight
       real(amrex_real) TAGLOC,TAGSIDE
       integer, PARAMETER :: nhalf=1
       real(amrex_real) xsten(-nhalf:nhalf,SDIM)
       real(amrex_real) xsten_n(-nhalf:nhalf,SDIM)
       real(amrex_real) xmain(SDIM)
       real(amrex_real) xside(SDIM)
       integer local_mask

       expan_ptr=>expan
       expan_comp_ptr=>expan_comp

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       !! Sanity checks

       if ((im_source.ge.1).and.(im_source.le.num_materials)) then
        ! do nothing
       else
        print *,"im_source invalid"
        stop
       endif
       if ((im_dest.ge.1).and.(im_dest.le.num_materials)) then
        ! do nothing
       else
        print *,"im_dest invalid"
        stop
       endif
       if (im_dest.eq.im_source) then
        print *,"im_dest or im_source invalid"
        stop
       endif

       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid in distribute_expansion"
        stop
       end if
       if ((indexEXP.lt.0).or.(indexEXP.ge.2*num_interfaces)) then
        print *,"indexEXP invalid"
        stop
       endif
       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if (ngrow_distance.eq.4) then
        ! do nothing
       else
        print *,"ngrow_distance invalid"
        stop
       endif

       maskcov_ptr=>maskcov
       LS_ptr=>LS
       tag_ptr=>tag
       tag_comp_ptr=>tag_comp
       weightfab_ptr=>weightfab
       weight_comp_ptr=>weight_comp

       call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
       call checkbound_array(fablo,fabhi,LS_ptr,ngrow_distance,-1)
       call checkbound_array(fablo,fabhi,expan_ptr,ngrow_distance,-1)
       call checkbound_array(fablo,fabhi,expan_comp_ptr,ngrow_distance,-1)
       call checkbound_array1(fablo,fabhi,tag_ptr,ngrow_distance,-1)
       call checkbound_array1(fablo,fabhi,tag_comp_ptr,ngrow_distance,-1)
       call checkbound_array1(fablo,fabhi,weightfab_ptr,ngrow_distance,-1)
       call checkbound_array1(fablo,fabhi,weight_comp_ptr,ngrow_distance,-1)

       ! Iterate over the box
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        local_mask=NINT(maskcov(D_DECL(i,j,k)))

        if (local_mask.eq.1) then

         TAGLOC=tag(D_DECL(i,j,k))
         if ((TAGLOC.eq.one).or. & ! doner cell
             (TAGLOC.eq.two)) then ! receiver cell
          call stencilbox(i,j,k,fablo,fabhi,stenlo,stenhi,ngrow_distance)

          do dir=1,SDIM
           if (stenlo(dir).lt.domlo(dir)) then
            stenlo(dir)=domlo(dir)
           endif
           if (stenhi(dir).gt.domhi(dir)) then
            stenhi(dir)=domhi(dir)
           endif
          enddo

          call gridsten_level(xsten,i,j,k,level,nhalf)
          do dir=1,SDIM
           xmain(dir)=xsten(0,dir)
          enddo
            
          total_weight=zero
          do k_n=stenlo(3),stenhi(3)
          do j_n=stenlo(2),stenhi(2)
          do i_n=stenlo(1),stenhi(1)

           ! if there is receiver neighbor cell
           TAGSIDE=tag(D_DECL(i_n,j_n,k_n))
           if (TAGSIDE.eq.two) then ! receiver cell

            call gridsten_level(xsten_n,i_n,j_n,k_n,level,nhalf)
            do dir=1,SDIM
             xside(dir)=xsten_n(0,dir)
            enddo
            call redistribute_weight(xmain,xside,crit_weight)
            total_weight=total_weight+crit_weight

           else if (TAGSIDE.eq.one) then ! donate
            ! do nothing
           else if (TAGSIDE.eq.zero) then
            ! do nothing
           else
            print *,"TAGSIDE invalid"
            stop
           endif

          enddo
          enddo
          enddo ! i_n,j_n,k_n

          if (total_weight.ge.zero) then
           ! do nothing
          else
           print *,"total_weight invalid"
           stop
          endif

          weightfab(D_DECL(i,j,k))=total_weight

         else if (TAGLOC.eq.zero) then
          ! do nothing
         else
          print *,"TAGLOC invalid"
          stop
         endif 

          ! ---------------- DISTRIBUTE FOR COMPLEMENT ----------------

         TAGLOC=tag_comp(D_DECL(i,j,k))
         if ((TAGLOC.eq.one).or. & ! doner cell
             (TAGLOC.eq.two)) then ! receiver cell
          call stencilbox(i,j,k,fablo,fabhi,stenlo,stenhi,ngrow_distance)

          do dir=1,SDIM
           if (stenlo(dir).lt.domlo(dir)) then
            stenlo(dir)=domlo(dir)
           endif
           if (stenhi(dir).gt.domhi(dir)) then
            stenhi(dir)=domhi(dir)
           endif
          enddo

          call gridsten_level(xsten,i,j,k,level,nhalf)
          do dir=1,SDIM
           xmain(dir)=xsten(0,dir)
          enddo
            
          total_weight=zero
          do k_n=stenlo(3),stenhi(3)
          do j_n=stenlo(2),stenhi(2)
          do i_n=stenlo(1),stenhi(1)

           ! if there is receiver neighbor cell
           TAGSIDE=tag_comp(D_DECL(i_n,j_n,k_n))
           if(TAGSIDE.eq.two) then ! receiver cell

            call gridsten_level(xsten_n,i_n,j_n,k_n,level,nhalf)
            do dir=1,SDIM
             xside(dir)=xsten_n(0,dir)
            enddo
            call redistribute_weight(xmain,xside,crit_weight)
            total_weight=total_weight+crit_weight
           else if (TAGSIDE.eq.one) then ! donate
            ! do nothing
           else if (TAGSIDE.eq.zero) then
            ! do nothing
           else
            print *,"TAGSIDE invalid"
            stop
           endif

          enddo
          enddo
          enddo ! i_n,j_n,k_n

          if (total_weight.ge.zero) then
           ! do nothing
          else
           print *,"total_weight invalid"
           stop
          endif

          weight_comp(D_DECL(i,j,k))=total_weight

         else if (TAGLOC.eq.zero) then
          ! do nothing
         else
          print *,"TAGLOC invalid"
          stop
         end if ! receiving cell

        else if (local_mask.eq.0) then
         ! do nothing
        else
         print *,"local_mask invalid"
         stop
        endif

       end do ! k
       end do ! j
       end do ! i
      end subroutine fort_accept_weight

      subroutine fort_aggressive( &
       datatype, &
       warning_cutoff, &
       tilelo,tilehi, &
       fablo,fabhi, &
       growlo,growhi, &
       bfact, &
       dx, &
       scomp, &
       ncomp, &
       ndefined, &
       ngrow, &
       dir, &
       verbose, &
       force_check, &
       gridno, &
       ngrid,level,finest_level, &
       mf,DIMS(mf)) &
      bind(c,name='fort_aggressive')

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: datatype
      real(amrex_real), INTENT(in) :: warning_cutoff
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: growlo(SDIM),growhi(SDIM)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: scomp,ncomp,ndefined
      integer, INTENT(in) :: ngrow
      integer, INTENT(in) :: dir
      integer, INTENT(in) :: verbose
      integer, INTENT(in) :: force_check
      integer, INTENT(in) :: gridno,ngrid,level,finest_level
      integer, INTENT(in) :: DIMDEC(mf)
      real(amrex_real), INTENT(in), target :: mf(DIMV(mf),ndefined)
      real(amrex_real), pointer :: mf_ptr(D_DECL(:,:,:),:)
      real(amrex_real) :: critical_cutoff_low
      real(amrex_real) :: critical_cutoff_high

      mf_ptr=>mf

      critical_cutoff_low=-1.0D+30
      critical_cutoff_high=1.0D+30

      if (bfact.lt.1) then
       print *,"bfact invalid69"
       stop
      endif

      if (((verbose.eq.0).or.(verbose.eq.1)).and.(force_check.eq.0)) then
       ! do nothing
      else if ((verbose.eq.2).or.(force_check.eq.1)) then
       call aggressive_worker( &
        datatype, &
        warning_cutoff, &
        tilelo,tilehi, &
        fablo,fabhi, &
        growlo,growhi, &
        bfact, &
        dx, &
        scomp, &
        ncomp, &
        ndefined, &
        ngrow, &
        dir, &
        verbose, &
        force_check, &
        gridno,ngrid,level,finest_level, &
        mf_ptr, &
        critical_cutoff_low, &
        critical_cutoff_high)
      else
       print *,"verbose invalid"
       stop
      endif

      return
      end subroutine fort_aggressive

       ! "coarray fortran"  (MPI functionality built in)
       ! masknbr=1.0 in the interior
       !        =1.0 fine-fine ghost cells
       !        =0.0 coarse-fine ghost cells and outside the domain.
       ! mask=tag if not covered by level+1 or outside the domain.
      subroutine fort_vfrac_split( &
       nprocessed, &
       tid, &
       density_floor, &
       density_ceiling, &
       solidheat_flag, &
       freezing_model, &
       distribute_from_target, &
       constant_density_all_time, &
       velbc, &
       divu_outer_sweeps, &
       num_divu_outer_sweeps, &
       EILE_flag, &
       dir_counter, &
       normdir, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       bfact_f, &
       dt, &
       time, &
       cur_time, &
       passive_veltime, &
       LS,DIMS(LS), &  ! original data; ngrow=2
       den, &
       DIMS(den), &
       mom_den, &
       DIMS(mom_den), &
       tensor,DIMS(tensor), &
       refineden,DIMS(refineden), &
       velfab,DIMS(velfab), & !VELADVECT_MF
       PLICSLP,DIMS(PLICSLP), &  ! slope data
       snew,DIMS(snew), &  ! this is the result
       tennew,DIMS(tennew), & 
       refinedennew,DIMS(refinedennew), & 
       LSnew,DIMS(LSnew), &
       vof0,DIMS(vof0), &  
       mask,DIMS(mask), & !mask=1 if not covered by level+1 or outside domain
       masknbr,DIMS(masknbr), &
       umac_displace, & ! vel*dt
       DIMS(umac_displace), & 
       xlo,dx, &
       conserve,DIMS(conserve), & ! local variables
       xmomside,DIMS(xmomside), & ! 1..2
       ymomside,DIMS(ymomside), &
       zmomside,DIMS(zmomside), &
       xmassside,DIMS(xmassside), & ! 1..2
       ymassside,DIMS(ymassside), &
       zmassside,DIMS(zmassside), &
       xmac_new,DIMS(xmac_new), &
       ymac_new,DIMS(ymac_new), &
       zmac_new,DIMS(zmac_new), &
       xmac_old,DIMS(xmac_old), &
       ymac_old,DIMS(ymac_old), &
       zmac_old,DIMS(zmac_old), &
       stokes_flow, &
       denconst_interface, & !unused: fort_vfrac_split
       nc_conserve, &
       map_forward, &
       recon_ncomp, &
       den_recon_ncomp, &
       ncomp_state, &
       nc_bucket, &
       verbose, &
       gridno,ngrid, &
       level, &
       finest_level, &
       dombc, &
       domlo,domhi) &
      bind(c,name='fort_vfrac_split')

      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, PARAMETER :: nhalf=1
      integer, INTENT(inout) :: nprocessed
      integer, INTENT(in) :: tid

      integer, INTENT(in) :: nc_conserve
      integer, PARAMETER :: ngrow=2
      integer, INTENT(in) :: stokes_flow
      real(amrex_real), INTENT(in) :: denconst_interface(num_interfaces)
      integer, INTENT(in) :: solidheat_flag
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)

      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      integer, INTENT(in) :: dombc(SDIM,2)
      integer, INTENT(in) :: divu_outer_sweeps
      integer, INTENT(in) :: num_divu_outer_sweeps
      integer, INTENT(in) :: EILE_flag
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(in) :: passive_veltime
      integer, INTENT(in) :: dir_counter
      integer, INTENT(in) :: normdir
      integer, INTENT(in) :: verbose
      integer :: force_check
      integer, INTENT(in) :: gridno,ngrid
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: recon_ncomp
      integer, INTENT(in) :: den_recon_ncomp
      integer, INTENT(in) :: ncomp_state
      integer, INTENT(in) :: nc_bucket
      integer :: nc_bucket_test
      integer, INTENT(in) :: map_forward
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: bfact_f
      integer, INTENT(in) :: constant_density_all_time(num_materials)
      real(amrex_real), INTENT(in) :: dt,time
       ! original data
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(mom_den)
      integer, INTENT(in) :: DIMDEC(tensor)
      integer, INTENT(in) :: DIMDEC(refineden)
      integer, INTENT(in) :: DIMDEC(velfab)
       ! slope data
      integer, INTENT(in) :: DIMDEC(PLICSLP)
       ! new data
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(tennew)
      integer, INTENT(in) :: DIMDEC(refinedennew)
      integer, INTENT(in) :: DIMDEC(LSnew)
       ! other vars
      integer, INTENT(in) :: DIMDEC(vof0)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: DIMDEC(masknbr)
      integer, INTENT(in) :: DIMDEC(umac_displace)
       ! local variables
      integer, INTENT(in) :: DIMDEC(conserve)

      integer, INTENT(in) :: DIMDEC(xmomside) 
      integer, INTENT(in) :: DIMDEC(ymomside) 
      integer, INTENT(in) :: DIMDEC(zmomside) 
      integer, INTENT(in) :: DIMDEC(xmassside) 
      integer, INTENT(in) :: DIMDEC(ymassside) 
      integer, INTENT(in) :: DIMDEC(zmassside) 
      integer, INTENT(in) :: DIMDEC(xmac_new) 
      integer, INTENT(in) :: DIMDEC(ymac_new) 
      integer, INTENT(in) :: DIMDEC(zmac_new) 
      integer, INTENT(in) :: DIMDEC(xmac_old) 
      integer, INTENT(in) :: DIMDEC(ymac_old) 
      integer, INTENT(in) :: DIMDEC(zmac_old) 

       ! FABS
       ! original data
      real(amrex_real), INTENT(in), target :: LS(DIMV(LS),num_materials)
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: den(DIMV(den),den_recon_ncomp)
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              mom_den(DIMV(mom_den),num_materials)
      real(amrex_real), pointer :: mom_den_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              tensor(DIMV(tensor),NUM_CELL_ELASTIC)
      real(amrex_real), pointer :: tensor_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              refineden(DIMV(refineden),NUM_CELL_REFINE_DENSITY)
      real(amrex_real), pointer :: refineden_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: & !VELADVECT_MF
           velfab(DIMV(velfab),STATE_NCOMP_VEL+STATE_NCOMP_PRES)
      real(amrex_real), pointer :: velfab_ptr(D_DECL(:,:,:),:)
       ! slope data
      real(amrex_real), INTENT(in), target :: PLICSLP(DIMV(PLICSLP),recon_ncomp)
      real(amrex_real), pointer :: PLICSLP_ptr(D_DECL(:,:,:),:)
       ! new data
      real(amrex_real), INTENT(inout), target :: snew(DIMV(snew),ncomp_state)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              tennew(DIMV(tennew),NUM_CELL_ELASTIC)
      real(amrex_real), pointer :: tennew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              refinedennew(DIMV(refinedennew),NUM_CELL_REFINE_DENSITY)
      real(amrex_real), pointer :: refinedennew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              LSnew(DIMV(LSnew),num_materials)
      real(amrex_real), pointer :: LSnew_ptr(D_DECL(:,:,:),:)
       ! other vars
      real(amrex_real), INTENT(in), target :: vof0(DIMV(vof0),num_materials)
      real(amrex_real), pointer :: vof0_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))
      ! =1 int. =1 fine-fine in domain =0 o.t.
      real(amrex_real), INTENT(in), target :: masknbr(DIMV(masknbr)) 
      real(amrex_real), pointer :: masknbr_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: umac_displace(DIMV(umac_displace))
      real(amrex_real), pointer :: umac_displace_ptr(D_DECL(:,:,:))
       ! local variables
      real(amrex_real), INTENT(inout), target ::  &
            conserve(DIMV(conserve),nc_conserve)
      real(amrex_real), pointer :: conserve_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(inout), target :: xmomside(DIMV(xmomside),2)
      real(amrex_real), pointer :: xmomside_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: ymomside(DIMV(ymomside),2)
      real(amrex_real), pointer :: ymomside_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: zmomside(DIMV(zmomside),2)
      real(amrex_real), pointer :: zmomside_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(inout), target :: xmassside(DIMV(xmassside),2)
      real(amrex_real), pointer :: xmassside_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: ymassside(DIMV(ymassside),2)
      real(amrex_real), pointer :: ymassside_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: zmassside(DIMV(zmassside),2)
      real(amrex_real), pointer :: zmassside_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(out), target :: xmac_new(DIMV(xmac_new))
      real(amrex_real), pointer :: xmac_new_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out), target :: ymac_new(DIMV(ymac_new))
      real(amrex_real), pointer :: ymac_new_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out), target :: zmac_new(DIMV(zmac_new))
      real(amrex_real), pointer :: zmac_new_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in), target :: xmac_old(DIMV(xmac_old))
      real(amrex_real), pointer :: xmac_old_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: ymac_old(DIMV(ymac_old))
      real(amrex_real), pointer :: ymac_old_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: zmac_old(DIMV(zmac_old))
      real(amrex_real), pointer :: zmac_old_ptr(D_DECL(:,:,:))
    
      real(amrex_real), INTENT(in) :: density_floor(num_materials)
      real(amrex_real), INTENT(in) :: density_ceiling(num_materials)

      integer, INTENT(in) :: velbc(SDIM,2)

      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      integer ibucket

      real(amrex_real) refine_den_bucket(num_materials_compressible)
      real(amrex_real) refine_vol_bucket(num_materials_compressible)

      real(amrex_real) xsten_crse(-nhalf:nhalf,SDIM)
      integer dir2
      integer iside
      integer vofcomp
      integer im
      integer incompressible_interface_flag
      real(amrex_real) dxmaxLS

      real(amrex_real) mom2(SDIM)
      real(amrex_real) xsten_MAC(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_accept(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_donate(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_target(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_depart(-nhalf:nhalf,SDIM)
      real(amrex_real) usten_accept(-nhalf:nhalf)
      real(amrex_real) usten_donate(-nhalf:nhalf)

      real(amrex_real) xdepartsize,xtargetsize,xloint,xhiint
      real(amrex_real) volint
      real(amrex_real) coeff(2)
      integer nmax
      integer ii,jj,kk
     
      integer veldir

      real(amrex_real) totalmass_depart
     
      integer istate,ispecies,igeom
    
      real(amrex_real) KE,local_internal
      real(amrex_real) :: vel_coarse(SDIM)
      real(amrex_real) :: vel_fine(SDIM)

      integer no_material_flag

      integer dencomp_data,statecomp_data,tempcomp_data,speccomp_data
      integer imap

      real(amrex_real) volcell_recon
      real(amrex_real) cencell_recon(SDIM)
      real(amrex_real) volcell_accept
      real(amrex_real) cencell_accept(SDIM)
      real(amrex_real) volcell_donate
      real(amrex_real) cencell_donate(SDIM)
      integer iii,jjj,kkk

      real(amrex_real) massdepart
      real(amrex_real) massdepart_mom

      integer idonate,jdonate,kdonate
      integer growlo(3),growhi(3)
      integer datatype

      integer istencil
      integer maskleft,maskright
      real(amrex_real) donate_data
      real(amrex_real) donate_density
      real(amrex_real) donate_mom_density
      real(amrex_real) ETcore

      integer icrse,jcrse,kcrse
      integer ifine,jfine,kfine
      integer nfine
      integer ifine_stencil,jfine_stencil,kfine_stencil
      integer ifine_stencil_lo,jfine_stencil_lo,kfine_stencil_lo
      integer ifine_stencil_hi,jfine_stencil_hi,kfine_stencil_hi
      integer nfine_stencil
      integer fine_offset
      integer im_comp,im_refine_density

      integer idonatelow
      integer idonatehigh

      real(amrex_real) voltotal_target
      real(amrex_real) voltotal_depart
      real(amrex_real) LS_voltotal_depart

      real(amrex_real) mofdata_grid(recon_ncomp)
      real(amrex_real) snew_hold(ncomp_state)
      real(amrex_real) tennew_hold(NUM_CELL_ELASTIC)
      real(amrex_real) mom_dencore(num_materials)
      real(amrex_real) dencore(num_materials)
      real(amrex_real) oldLS(num_materials)
      real(amrex_real) newLS(num_materials)
      real(amrex_real) newvfrac_weymouth(num_materials)
      real(amrex_real) newvfrac_cor(num_materials)
      real(amrex_real) newvfrac(num_materials)
      real(amrex_real) volmat_depart(num_materials)
      real(amrex_real) volmat_target(num_materials)
      real(amrex_real) volmat_depart_cor(num_materials)
      real(amrex_real) volmat_target_cor(num_materials)
      real(amrex_real) multi_volume(num_materials)
      real(amrex_real) multi_volume_grid(num_materials)
      real(amrex_real) multi_cen(SDIM,num_materials)
      real(amrex_real) multi_cen_grid(SDIM,num_materials)
      real(amrex_real) newcen(SDIM,num_materials)
      real(amrex_real) veldata(nc_bucket)

      integer ihalf
      integer check_intersection
      real(amrex_real) xsten_recon(-nhalf:nhalf,SDIM)

      real(amrex_real) warning_cutoff

      integer all_incomp
      real(amrex_real) vol_target_local
      integer k1lo,k1hi

      real(amrex_real) massfrac_parm(num_species_var+1)

      real(amrex_real) :: critical_cutoff_low
      real(amrex_real) :: critical_cutoff_high

      integer :: dir_local
      real(amrex_real) :: wt_oldvel
      integer :: zapvel
      integer :: iright,jright,kright
      integer :: ileft,jleft,kleft
      integer :: icell,jcell,kcell
      real(amrex_real) :: momface_total
      real(amrex_real) :: massface_total
      real(amrex_real) :: massquarter
      real(amrex_real) :: momquarter

      real(amrex_real) :: xclamped(SDIM)
      real(amrex_real) :: LS_clamped
      real(amrex_real) :: vel_clamped(SDIM)
      real(amrex_real) :: temperature_clamped
      real(amrex_real) :: local_temperature
      integer :: prescribed_flag

! fort_vfrac_split code starts here

      LS_ptr=>LS
      den_ptr=>den
      mom_den_ptr=>mom_den
      velfab_ptr=>velfab
      PLICSLP_ptr=>PLICSLP

      vof0_ptr=>vof0
      mask_ptr=>mask
      masknbr_ptr=>masknbr
      umac_displace_ptr=>umac_displace
      conserve_ptr=>conserve

      snew_ptr=>snew
      LSnew_ptr=>LSnew

      xmomside_ptr=>xmomside
      ymomside_ptr=>ymomside
      zmomside_ptr=>zmomside

      xmassside_ptr=>xmassside
      ymassside_ptr=>ymassside
      zmassside_ptr=>zmassside

      xmac_new_ptr=>xmac_new
      ymac_new_ptr=>ymac_new
      zmac_new_ptr=>zmac_new

      xmac_old_ptr=>xmac_old
      ymac_old_ptr=>ymac_old
      zmac_old_ptr=>zmac_old

      critical_cutoff_low=-1.0D+30
      critical_cutoff_high=1.0D+30

      if ((tid.lt.0).or. &
          (tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if ((num_materials_viscoelastic.ge.1).and. &
          (num_materials_viscoelastic.le.num_materials)) then
       tennew_ptr=>tennew
       tensor_ptr=>tensor
      else if (num_materials_viscoelastic.eq.0) then
       tennew_ptr=>snew
       tensor_ptr=>den
      else
       print *,"num_materials_viscoelastic invalid:fort_vfrac_split"
       stop
      endif


      if ((num_materials_compressible.ge.1).and. &
          (num_materials_compressible.le.num_materials)) then
       refinedennew_ptr=>refinedennew
       refineden_ptr=>refineden
      else if (num_materials_compressible.eq.0) then
       refinedennew_ptr=>snew
       refineden_ptr=>den
      else
       print *,"num_materials_compressible invalid:fort_vfrac_split"
       stop
      endif


      nmax=POLYGON_LIST_MAX ! in: fort_vfrac_split

      k1lo=0
      k1hi=0
      if (SDIM.eq.2) then
       ! do nothing
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid70"
       stop
      endif
      if (bfact_f.lt.1) then
       print *,"bfact_f invalid"
       stop
      endif
      if ((bfact.ne.bfact_f).and.(bfact.ne.2*bfact_f)) then
       print *,"bfact invalid71"
       stop
      endif
      if ((solidheat_flag.lt.0).or. &
          (solidheat_flag.gt.2)) then
       print *,"solidheat_flag invalid"
       stop
      endif

      if ((level.lt.0).or. &
          (level.gt.finest_level)) then
       print *,"level invalid vfrac split"
       stop
      endif
      if ((verbose.lt.0).or.(verbose.gt.2)) then
       print *,"verbose invalid"
       stop
      endif
      if ((gridno.lt.0).or.(gridno.ge.ngrid)) then
       print *,"gridno invalid in vfrac split"
       stop
      endif

      if (ngrow.eq.2) then
       ! do nothing
      else
       print *,"ngrow invalid: ",ngrow
       stop
      endif

      if (ncomp_state.ne.STATECOMP_STATES+ &
          num_materials*(num_state_material+ngeom_raw)+1) then
       print *,"ncomp_state invalid"
       stop
      endif

      do im=1,num_interfaces
       if (denconst_interface(im).gt.zero) then
        ! do nothing
       else if (denconst_interface(im).eq.zero) then
        ! do nothing
       else
        print *,"denconst_interface invalid fort_vfrac_split"
        print *,"index,value ",im,denconst_interface(im)
        stop
       endif
      enddo !im=1,num_interfaces

      all_incomp=1

      do im=1,num_materials

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).ge.1).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
        all_incomp=0
       else
        print *,"fort_material_type invalid"
        stop
       endif

       if ((density_floor(im).ge.zero).and. &
           (density_floor(im).lt.fort_denconst(im))) then
        !do nothing
       else
        print *,"density_floor invalid: ",im,density_floor(im);
        stop
       endif
       if ((density_ceiling(im).gt.zero).and. &
           (density_ceiling(im).ge.fort_denconst(im))) then
        !do nothing
       else
        print *,"density_ceiling invalid: ",im,density_ceiling(im)
        stop
       endif

      enddo  ! im=1..num_materials

      if ((all_incomp.eq.0).or. &
          (all_incomp.eq.1)) then
       ! do nothing
      else
       print *,"all_incomp invalid: ",all_incomp
       stop
      endif

      if (num_state_material.ne. &
          num_state_base+num_species_var) then
       print *,"num_state_material invalid"
       stop
      endif

      if (NUM_CELL_ELASTIC.eq. &
          num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) then
       ! do nothing
      else
       print *,"NUM_CELL_ELASTIC invalid"
       stop
      endif
      if (ENUM_NUM_TENSOR_TYPE.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ENUM_NUM_TENSOR_TYPE invalid"
       stop
      endif 

      if ((num_materials_viscoelastic.ge.0).and. &
          (num_materials_viscoelastic.le.num_materials)) then
       ! do nothing
      else
       print *,"num_materials_viscoelastic invalid:fort_vfrac_split"
       stop
      endif

      if (NUM_CELL_REFINE_DENSITY.eq. &
          num_materials_compressible*ENUM_NUM_REFINE_DENSITY_TYPE) then
       ! do nothing
      else
       print *,"NUM_CELL_REFINE_DENSITY invalid"
       stop
      endif
      if (ENUM_NUM_REFINE_DENSITY_TYPE.eq.4*(SDIM-1)) then
       ! do nothing
      else
       print *,"ENUM_NUM_REFINE_DENSITY_TYPE invalid"
       stop
      endif 

      if ((num_materials_compressible.ge.0).and. &
          (num_materials_compressible.le.num_materials)) then
       ! do nothing
      else
       print *,"num_materials_compressible invalid:fort_vfrac_split"
       stop
      endif



      if ((divu_outer_sweeps.ge.0).and. &
          (divu_outer_sweeps.lt.num_divu_outer_sweeps)) then
       ! do nothing
      else
       print *,"divu_outer_sweeps invalid: ",divu_outer_sweeps
       stop
      endif

      if ((EILE_flag.eq.-1).or. & ! Weymouth and Yue
          (EILE_flag.eq.1).or.  & ! EILE
          (EILE_flag.eq.2).or.  & ! always EI
          (EILE_flag.eq.3)) then  ! always LE
       ! do nothing
      else 
       print *,"EILE flag invalid"
       stop
      endif

      if ((dir_counter.lt.0).or.(dir_counter.ge.SDIM)) then
       print *,"dir_counter invalid"
       stop
      endif

      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid: ",normdir
       stop
      endif

      if (den_recon_ncomp.ne.num_materials*num_state_material) then
       print *,"den_recon_ncomp invalid: ",den_recon_ncomp
       stop
      endif

      if (recon_ncomp.ne.num_materials*ngeom_recon) then
       print *,"recon_ncomp invalid"
       stop
      endif

      if ((map_forward.ne.0).and.(map_forward.ne.1)) then
       print *,"map_forward invalid"
       stop
      endif

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension crash"
        stop
       endif
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"levelrz invalid vfrac split"
       stop
      endif

      call get_dxmaxLS(dx,bfact,dxmaxLS)
      if (dxmaxLS.gt.zero) then
       !do nothing
      else
       print *,"dxmaxLS invalid: ",dxmaxLS
       stop
      endif

       ! SANITY CHECKS TO MAKE SURE THAT input FABs have the expected number of
       ! ghost cells.

       ! original data
      call checkbound_array(fablo,fabhi,LS_ptr,ngrow,-1)
       ! ngrow=2
      call checkbound_array(fablo,fabhi,den_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,mom_den_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,tensor_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,refineden_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,velfab_ptr,ngrow,-1)
       ! slope data
      call checkbound_array(fablo,fabhi,PLICSLP_ptr,ngrow,-1)
       ! new data
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,tennew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,refinedennew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LSnew_ptr,1,-1)
       ! other vars
      call checkbound_array(fablo,fabhi,vof0_ptr,ngrow,-1)
      call checkbound_array1(fablo,fabhi,mask_ptr,ngrow,-1)
      call checkbound_array1(fablo,fabhi,masknbr_ptr,ngrow,-1)
     
       ! example: imac=0, left side; then the parcel at imac=-1, left side
       ! might be advected: imac=-1, left side is in icell=-2.  icell=-2 is
       ! advected using imac=-2 and imac=-1.
      call checkbound_array1(fablo,fabhi,umac_displace_ptr, &
              ngrow,normdir)

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (normdir.eq.0) then
       ii=1
      else if (normdir.eq.1) then
       jj=1
      else if ((normdir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"normdir invalid"
       stop
      endif

      growlo(3)=0
      growhi(3)=0

      call checkbound_array(fablo,fabhi,xmomside_ptr,1,-1)
      call checkbound_array(fablo,fabhi,ymomside_ptr,1,-1)
      call checkbound_array(fablo,fabhi,zmomside_ptr,1,-1)

      call checkbound_array(fablo,fabhi,xmassside_ptr,1,-1)
      call checkbound_array(fablo,fabhi,ymassside_ptr,1,-1)
      call checkbound_array(fablo,fabhi,zmassside_ptr,1,-1)

      call checkbound_array1(fablo,fabhi,xmac_new_ptr,0,0)
      call checkbound_array1(fablo,fabhi,ymac_new_ptr,0,1)
      call checkbound_array1(fablo,fabhi,zmac_new_ptr,0,SDIM-1)

      call checkbound_array1(fablo,fabhi,xmac_old_ptr,ngrow,0)
      call checkbound_array1(fablo,fabhi,ymac_old_ptr,ngrow,1)
      call checkbound_array1(fablo,fabhi,zmac_old_ptr,ngrow,SDIM-1)

      if (nc_conserve.ne.CISLCOMP_CONS_NCOMP*ENUM_NUM_REFINE_DENSITY_TYPE) then
       print *,"nc_conserve invalid: ",nc_conserve
       stop
      endif

      call checkbound_array(fablo,fabhi,conserve_ptr,ngrow,-1)
     
      force_check=0
      datatype=0 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow)
 
      warning_cutoff=two
      call aggressive_worker( &
       datatype, &
       warning_cutoff, &
       tilelo,tilehi, &
       fablo,fabhi, &
       growlo,growhi, &
       bfact, &
       dx, &
       0, & !scomp=0
       num_materials*ngeom_recon, &
       num_materials*ngeom_recon, &
       ngrow, &
       -1, & ! dir
       verbose, &
       force_check, &
       gridno,ngrid, &
       level,finest_level, &
       PLICSLP_ptr, &
       critical_cutoff_low, &
       critical_cutoff_high)

      warning_cutoff=1.0e+15
      call aggressive_worker( &
       datatype, & 
       warning_cutoff, &
       tilelo,tilehi, &
       fablo,fabhi, &
       growlo,growhi, &
       bfact, &
       dx, &
       0, &  !scomp=0
       nc_conserve, &
       nc_conserve, &
       ngrow, &
       -1, & !dir
       verbose, &
       force_check, &
       gridno,ngrid, &
       level, &
       finest_level, &
       conserve_ptr, &
       critical_cutoff_low, &
       critical_cutoff_high)

      nc_bucket_test=CISLCOMP_NCOMP
      if (nc_bucket_test.ne.nc_bucket) then
       print *,"nc_bucket invalid: ",nc_bucket_test,nc_bucket
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow)

      do kcrse=growlo(3),growhi(3)
      do jcrse=growlo(2),growhi(2)
      do icrse=growlo(1),growhi(1)

       call gridsten_level(xsten_crse,icrse,jcrse,kcrse,level,nhalf)

       do veldir=1,SDIM
        vel_coarse(veldir)=velfab(D_DECL(icrse,jcrse,kcrse),veldir)
       enddo

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        if (xsten_crse(0,1).lt.zero) then
         do veldir=1,SDIM
          vel_coarse(veldir)=zero
         enddo
         if (icrse.lt.0) then
          !do nothing
         else
          print *,"expecting icrse<0"
          stop
         endif
        endif
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        if (xsten_crse(0,1).lt.zero) then
         do veldir=1,SDIM
          vel_coarse(veldir)=zero
         enddo
         if (icrse.lt.0) then
          !do nothing
         else
          print *,"expecting icrse<0"
          stop
         endif
        endif
       else
        print *,"levelrz invalid fort_vfrac_split (vel_coarse)"
        stop
       endif

       kfine=0
#if (AMREX_SPACEDIM==3)
       do kfine=0,1
#endif
       do jfine=0,1
       do ifine=0,1
        nfine=4*kfine+2*jfine+ifine+1
        fine_offset=CISLCOMP_CONS_NCOMP*(nfine-1)

        veldir=1
        if (ifine.eq.0) then
         vel_fine(veldir)=xmac_old(D_DECL(icrse,jcrse,kcrse))
        else if (ifine.eq.1) then
         vel_fine(veldir)=xmac_old(D_DECL(icrse+1,jcrse,kcrse))
        else
         print *,"ifine invalid"
         stop
        endif
        veldir=2
        if (jfine.eq.0) then
         vel_fine(veldir)=ymac_old(D_DECL(icrse,jcrse,kcrse))
        else if (jfine.eq.1) then
         vel_fine(veldir)=ymac_old(D_DECL(icrse,jcrse+1,kcrse))
        else
         print *,"jfine invalid"
         stop
        endif
        if (SDIM.eq.3) then
         veldir=SDIM
         if (kfine.eq.0) then
          vel_fine(veldir)=zmac_old(D_DECL(icrse,jcrse,kcrse))
         else if (kfine.eq.1) then
          vel_fine(veldir)=zmac_old(D_DECL(icrse,jcrse,kcrse+1))
         else
          print *,"kfine invalid: ",kfine
          stop
         endif
        endif

        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then
         if ((xsten_crse(-1,1).le.EPS2*dx(1)).and.(ifine.eq.0)) then
          vel_fine(1)=zero
          if (icrse.le.0) then
           !do nothing
          else
           print *,"expecting icrse<=0"
           stop
          endif
         endif
         if (xsten_crse(0,1).lt.zero) then
          do veldir=1,SDIM
           vel_fine(veldir)=zero
          enddo
          if (icrse.lt.0) then
           !do nothing
          else
           print *,"expecting icrse<0"
           stop
          endif
         endif
        else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
         if ((xsten_crse(-1,1).le.EPS2*dx(1)).and.(ifine.eq.0)) then
          vel_fine(1)=zero
          if (icrse.le.0) then
           !do nothing
          else
           print *,"expecting icrse<=0"
           stop
          endif
         endif
         if (xsten_crse(0,1).lt.zero) then
          do veldir=1,SDIM
           vel_fine(veldir)=zero
          enddo
          if (icrse.lt.0) then
           !do nothing
          else
           print *,"expecting icrse<0"
           stop
          endif
         endif
        else
         print *,"levelrz invalid fort_vfrac_split (vel_fine)"
         stop
        endif

        ! KE=u dot u/2
        KE=zero
        do veldir=1,SDIM
         conserve(D_DECL(icrse,jcrse,kcrse),fine_offset+veldir)= &
              vel_fine(veldir)
         KE=KE+vel_coarse(veldir)**2
        enddo
        KE=half*KE

        if (KE.ge.zero) then
         ! do nothing
        else
         print *,"KE invalid: ",KE
         stop
        endif

        im_refine_density=0

        do im=1,num_materials
         istate=ENUM_DENVAR+1
         dencomp_data=(im-1)*num_state_material+istate
         if (is_compressible_mat(im).eq.0) then
          dencore(im)=den(D_DECL(icrse,jcrse,kcrse),dencomp_data)
          mom_dencore(im)=mom_den(D_DECL(icrse,jcrse,kcrse),im)
          if (constant_density_all_time(im).eq.1) then
           if (abs(dencore(im)-fort_denconst(im)).le. &
               fort_denconst(im)*EPS_8_4) then
            ! do nothing
           else
            print *,"dencore(im) invalid"
            print *,"im,icrse,jcrse,kcrse,den ",im,icrse,jcrse,kcrse,dencore(im)
            print *,"fort_denconst(im) ",fort_denconst(im)
            print *,"dencomp_data=",dencomp_data
            print *,"normdir=",normdir
            stop
           endif
          else if (constant_density_all_time(im).eq.0) then 
           ! do nothing
          else
           print *,"constant_density_all_time invalid"
           stop
          endif
         else if (is_compressible_mat(im).eq.1) then

          im_refine_density=im_refine_density+1
          if (fort_im_refine_density_map(im_refine_density).eq.im-1) then
           !do nothing
          else
           print *,"fort_im_refine_density_map invalid"
           stop
          endif

          dencore(im)=refineden(D_DECL(icrse,jcrse,kcrse), &
             (im_refine_density-1)*ENUM_NUM_REFINE_DENSITY_TYPE+nfine)
          mom_dencore(im)=dencore(im)

          if (constant_density_all_time(im).eq.0) then
           !do nothing
          else
           print *,"expecting constant_density_all_time=0"
           stop
          endif

         else
          print *,"is_compressible(im) invalid"
          stop
         endif
         if (dencore(im).gt.zero) then
          ! do nothing
         else
          print *,"density must be positive fort_vfrac_split"
          print *,"im,dencore(im) ",im,dencore(im)
          print *,"im,fort_denconst(im) ",im,fort_denconst(im)
          stop
         endif  

         if (mom_dencore(im).gt.zero) then
          ! do nothing
         else
          print *,"mom_density must be positive fort_vfrac_split"
          print *,"im,mom_dencore(im) ",im,mom_dencore(im)
          print *,"im,fort_denconst(im) ",im,fort_denconst(im)
          stop
         endif  

        enddo ! im=1..num_materials

        do im=1,num_materials

         ! in: fort_vfrac_split
         istate=1
         do while (istate.le.num_state_material)

          if (istate.eq.ENUM_DENVAR+1) then ! Density
           dencomp_data=(im-1)*num_state_material+ENUM_DENVAR+1
           conserve(D_DECL(icrse,jcrse,kcrse), &
             fine_offset+CISLCOMP_STATES+dencomp_data)=dencore(im)
           istate=istate+1
          else if (istate.eq.ENUM_TEMPERATUREVAR+1) then ! Temperature
           tempcomp_data=(im-1)*num_state_material+ENUM_TEMPERATUREVAR+1
           local_temperature=den(D_DECL(icrse,jcrse,kcrse),tempcomp_data)
           if ((is_compressible_mat(im).eq.0).or. &
               (fort_material_conservation_form(im).eq.0)) then
             ! den * T
            conserve(D_DECL(icrse,jcrse,kcrse), &
              fine_offset+CISLCOMP_STATES+tempcomp_data)= &
                   dencore(im)*local_temperature

           else if ((is_compressible_mat(im).eq.1).and. &
                    (fort_material_conservation_form(im).eq.1)) then

            call init_massfrac_parm(dencore(im),massfrac_parm,im)
            do ispecies=1,num_species_var
             massfrac_parm(ispecies)= &
               den(D_DECL(icrse,jcrse,kcrse),tempcomp_data+ispecies)
             if (massfrac_parm(ispecies).ge.zero) then
              ! do nothing
             else
              print *,"massfrac_parm(ispecies) invalid: ", &
                ispecies,massfrac_parm(ispecies)
              stop
             endif
            enddo ! do ispecies=1,num_species_var

            ! den * (u dot u/2 + cv T)
            call INTERNAL_material(dencore(im),massfrac_parm, &
              local_temperature,local_internal, &
              fort_material_type(im),im)
            conserve(D_DECL(icrse,jcrse,kcrse), &
              fine_offset+CISLCOMP_STATES+tempcomp_data)= &
                   dencore(im)*(KE+local_internal) 

           else
            print *, &
             "is_compressible_mat or fort_material_conservation_form invalid"
            stop
           endif
           istate=istate+1

          else if ((istate.eq.num_state_base+1).and. &
                   (num_species_var.gt.0)) then 
           ! den * Y
           do ispecies=1,num_species_var
            speccomp_data=(im-1)*num_state_material+num_state_base+ispecies
            conserve(D_DECL(icrse,jcrse,kcrse), &
             fine_offset+CISLCOMP_STATES+speccomp_data)= &
               dencore(im)*den(D_DECL(icrse,jcrse,kcrse),speccomp_data)
            istate=istate+1
           enddo ! ispecies=1..num_species_var
          else 
           print *,"istate invalid"
           stop
          endif

         enddo ! do while (istate.le.num_state_material)

         if (dencore(im).gt.zero) then 
          ! do nothing
         else
          print *,"dencore must be positive"
          stop
         endif

        enddo ! im=1..num_materials
       enddo !ifine
       enddo !jfine
#if (AMREX_SPACEDIM==3)
       enddo !kfine
#endif

      enddo !icrse 
      enddo !jcrse
      enddo !kcrse (cell center "conserved" variables) 

       ! xmomside,ymomside,zmomside already init to 0 
       ! xmassside,ymassside,zmassside already init to 0 

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,1)

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       growlo(1)=max(0,growlo(1))
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
       growlo(1)=max(0,growlo(1))
      else
       print *,"levelrz invalid fort_vfrac_split (growlo(1))"
       stop
      endif

      do kcrse=growlo(3),growhi(3)
      do jcrse=growlo(2),growhi(2)
      do icrse=growlo(1),growhi(1)
       nprocessed=nprocessed+1

       call gridsten_level(xsten_crse,icrse,jcrse,kcrse,level,nhalf)

       voltotal_depart=zero

       do istate=1,nc_bucket
        veldata(istate)=zero
       enddo

       incompressible_interface_flag=0

       do im=1,num_materials
        oldLS(im)=LS(D_DECL(icrse,jcrse,kcrse),im)
       enddo

       do im=1,num_materials
        if (oldLS(im).ge.-incomp_thickness*dxmaxLS) then
         if (is_rigid(im).eq.1) then
          incompressible_interface_flag=1
         else if (is_rigid(im).eq.0) then 
          !do nothing
         else 
          print *,"is_rigid(im) invalid: ",im,is_rigid(im)
          stop
         endif
        else if (oldLS(im).le.-incomp_thickness*dxmaxLS) then
         ! do nothing
        else
         print *,"oldLS(im) corrupt,fort_vfrac_split"
         print *,"im,oldLS(im): ",im,oldLS(im)
         stop
        endif 
       enddo !im=1,num_materials

       idonate=icrse
       jdonate=jcrse
       kdonate=kcrse

       kfine=0
#if (AMREX_SPACEDIM==3)
       do kfine=0,1
#endif
       do jfine=0,1
       do ifine=0,1

        do im_comp=1,num_materials_compressible
         refine_den_bucket(im_comp)=zero
         refine_vol_bucket(im_comp)=zero
        enddo

        nfine=4*kfine+2*jfine+ifine+1

        call CISBOXFINE(xsten_accept,nhalf, &
         xlo,dx, &
         icrse,jcrse,kcrse, &
         ifine,jfine,kfine, &
         bfact,level, &
         volcell_accept,cencell_accept,SDIM)

        if (volcell_accept.gt.zero) then
         ! do nothing
        else
         print *,"volcell_accept invalid: ",volcell_accept
         stop
        endif

        do ihalf=-1,1
        do dir2=1,SDIM
         xsten_target(ihalf,dir2)=xsten_accept(ihalf,dir2)
         xsten_depart(ihalf,dir2)=xsten_accept(ihalf,dir2)
        enddo
        enddo

        usten_accept(-1)=umac_displace(D_DECL(icrse,jcrse,kcrse))
        usten_accept(1)=umac_displace(D_DECL(icrse+ii,jcrse+jj,kcrse+kk))

        if (normdir.eq.0) then

         if (levelrz.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if ((levelrz.eq.COORDSYS_RZ).or. &
                  (levelrz.eq.COORDSYS_CYLINDRICAL)) then
          if (icrse.le.0) then
           usten_accept(-1)=zero
          endif
          if (icrse.lt.0) then
           usten_accept(1)=zero
           print *,"expecting icrse>=0"
           stop
          endif
         else
          print *,"levelrz invalid"
          stop
         endif

        endif

        usten_accept(0)=half*(usten_accept(-1)+usten_accept(1))

        if (usten_accept(-1).gt.zero) then
         idonatelow=-1
        else if (usten_accept(-1).le.zero) then
         idonatelow=0
        else
         print *,"usten_accept(-1) invalid: ",usten_accept(-1)
         stop
        endif

        if (usten_accept(1).lt.zero) then
         idonatehigh=1
        else if (usten_accept(1).ge.zero) then
         idonatehigh=0
        else
         print *,"usten_accept(1) invalid: ",usten_accept(1)
         stop
        endif

        if (normdir.eq.0) then
         if (ifine.eq.0) then
          usten_accept(1)=usten_accept(0)
         else if (ifine.eq.1) then
          usten_accept(-1)=usten_accept(0)
         else
          print *,"ifine invalid"
          stop
         endif
        else if (normdir.eq.1) then
         if (jfine.eq.0) then
          usten_accept(1)=usten_accept(0)
         else if (jfine.eq.1) then
          usten_accept(-1)=usten_accept(0)
         else
          print *,"jfine invalid"
          stop
         endif
        else if ((normdir.eq.2).and.(SDIM.eq.3)) then
         if (kfine.eq.0) then
          usten_accept(1)=usten_accept(0)
         else if (kfine.eq.1) then
          usten_accept(-1)=usten_accept(0)
         else
          print *,"kfine invalid: ",kfine
          stop
         endif
        else
         print *,"normdir invalid: ",normdir
         stop
        endif
        usten_accept(0)=half*(usten_accept(-1)+usten_accept(1))

        do istencil=idonatelow,idonatehigh
 
         if (normdir.eq.0) then
          idonate=icrse+istencil
         else if (normdir.eq.1) then 
          jdonate=jcrse+istencil
         else if ((normdir.eq.2).and.(SDIM.eq.3)) then
          kdonate=kcrse+istencil
         else
          print *,"normdir invalid"
          stop
         endif
         kfine_stencil_lo=kfine
         jfine_stencil_lo=jfine
         ifine_stencil_lo=ifine
         kfine_stencil_hi=kfine
         jfine_stencil_hi=jfine
         ifine_stencil_hi=ifine
         if (normdir.eq.0) then
          ifine_stencil_lo=0 
          ifine_stencil_hi=1
         else if (normdir.eq.1) then 
          jfine_stencil_lo=0 
          jfine_stencil_hi=1
         else if ((normdir.eq.2).and.(SDIM.eq.3)) then
          kfine_stencil_lo=0 
          kfine_stencil_hi=1
         else
          print *,"normdir invalid: ",normdir
          stop
         endif
         kfine_stencil=0
#if (AMREX_SPACEDIM==3)
         do kfine_stencil=kfine_stencil_lo,kfine_stencil_hi
#endif
         do jfine_stencil=jfine_stencil_lo,jfine_stencil_hi
         do ifine_stencil=ifine_stencil_lo,ifine_stencil_hi

          nfine_stencil=4*kfine_stencil+2*jfine_stencil+ifine_stencil+1
          fine_offset=CISLCOMP_CONS_NCOMP*(nfine_stencil-1)

          call CISBOX(xsten_recon,1, &
           xlo,dx, &
           idonate,jdonate,kdonate, &
           bfact,level, &
           volcell_recon,cencell_recon,SDIM)
 
          call CISBOXFINE(xsten_donate,1, &
           xlo,dx, &
           idonate,jdonate,kdonate, &
           ifine_stencil,jfine_stencil,kfine_stencil, &
           bfact,level, &
           volcell_donate,cencell_donate,SDIM)

          check_intersection=1

          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if (levelrz.eq.COORDSYS_RZ) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if (xsten_recon(0,1).le.EPS2*dx(1)) then
            check_intersection=0
            print *,"idonate invalid (RZ): ",idonate
            print *,"[ijk]crse ",icrse,jcrse,kcrse
            print *,"normdir, [ijk]donate ",normdir,idonate,jdonate,kdonate
            stop
           endif
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           if (xsten_recon(0,1).le.EPS2*dx(1)) then
            check_intersection=0
            print *,"idonate invalid (RTZ): ",idonate
            stop
           endif
          else
           print *,"levelrz invalid add to bucket 2"
           stop
          endif

          if (check_intersection.eq.1) then 

           usten_donate(-1)=umac_displace(D_DECL(idonate,jdonate,kdonate))
           usten_donate(1)= &
              umac_displace(D_DECL(idonate+ii,jdonate+jj,kdonate+kk))

           if (normdir.eq.0) then

            if (levelrz.eq.COORDSYS_CARTESIAN) then
             ! do nothing
            else if ((levelrz.eq.COORDSYS_RZ).or. &
                     (levelrz.eq.COORDSYS_CYLINDRICAL)) then
             if (idonate.le.0) then
              usten_donate(-1)=zero
             endif
             if (idonate.lt.0) then
              usten_donate(1)=zero
              print *,"idonate invalid: ",idonate
              print *,"[ijk]donate: ",idonate,jdonate,kdonate
              stop
             endif
            else
             print *,"levelrz invalid"
             stop
            endif

           endif

           usten_donate(0)=half*(usten_donate(-1)+usten_donate(1))

           if (normdir.eq.0) then
            if (ifine_stencil.eq.0) then
             usten_donate(1)=usten_donate(0)
            else if (ifine_stencil.eq.1) then
             usten_donate(-1)=usten_donate(0)
            else
             print *,"ifine_stencil invalid"
             stop
            endif
           else if (normdir.eq.1) then
            if (jfine_stencil.eq.0) then
             usten_donate(1)=usten_donate(0)
            else if (jfine_stencil.eq.1) then
             usten_donate(-1)=usten_donate(0)
            else
             print *,"jfine_stencil invalid"
             stop
            endif
           else if ((normdir.eq.2).and.(SDIM.eq.3)) then
            if (kfine_stencil.eq.0) then
             usten_donate(1)=usten_donate(0)
            else if (kfine_stencil.eq.1) then
             usten_donate(-1)=usten_donate(0)
            else
             print *,"kfine_stencil invalid: ",kfine_stencil
             stop
            endif
           else
            print *,"normdir invalid: ",normdir
            stop
           endif
           usten_donate(0)=half*(usten_donate(-1)+usten_donate(1))

             ! normdir=0..sdim-1
           call derive_mappings( &
            xsten_accept, &
            xsten_donate, &
            xsten_target, &
            xsten_depart, &
            usten_accept, &
            usten_donate, &
            xdepartsize, &
            xtargetsize, &
            xloint, &
            xhiint, &
            volint, &
            coeff, &
            bfact, & !only used for sanity checks
            dx, &
            map_forward, &
            normdir)

           if (volint.gt.zero) then  

             ! we are inside the istencil loop.
            LS_voltotal_depart=zero

            do dir2=1,num_materials*ngeom_recon
             mofdata_grid(dir2)= &
              PLICSLP(D_DECL(idonate,jdonate,kdonate),dir2)
            enddo

              ! the volumes and centroids are tessellating for the fluid
              ! materials, but not the solid materials.  Solid materials are
              ! immersed into the domain.

            call multi_get_volume_grid_and_map( &
              normdir, & ! normdir=0..sdim-1
              coeff, &
              bfact,dx, &
              xsten_recon,1, &
              mofdata_grid, &
              xsten_depart,1, &
              multi_volume_grid, & ! intersection of departure with grid.
              multi_cen_grid, &
              multi_volume, & ! intersection of target with grid.
              multi_cen, &
              geom_xtetlist_uncapt(1,1,1,tid+1), &
              nmax, &
              nmax, &
              SDIM)

             ! normdir=0..sdim-1
            do im=1,num_materials
             vofcomp=(im-1)*ngeom_recon+1
  
              ! fluid materials tessellate the domain. 
             if (is_rigid(im).eq.0) then 
              LS_voltotal_depart=LS_voltotal_depart+ &
               multi_volume_grid(im)
             else if (is_rigid(im).eq.1) then
              ! do nothing
             else
              print *,"is_rigid invalid GODUNOV_3D.F90"
              stop
             endif

            enddo  ! im=1,..,num_materials

            if (LS_voltotal_depart.gt.zero) then
             ! do nothing
            else
             print *,"EPS_8_4= ",EPS_8_4
             print *,"LS_voltotal_depart bust (multi_get_volume_grid_and_map)"
             print *,"LS_voltotal_depart ",LS_voltotal_depart
             print *,"map_forward,volint ",map_forward,volint
             print *,"istencil ",istencil
             print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
             print *,"level,finest_level ",level,finest_level
             do im=1,num_materials 
              vofcomp=(im-1)*ngeom_recon+1
              print *,"im,multi_volume_grid ",im,multi_volume_grid(im)
              print *,"im,multi_volume ",im,multi_volume(im)
              print *,"im,vfrac ",im,mofdata_grid(vofcomp)
              print *,"im,flag ",im,mofdata_grid(vofcomp+SDIM+1)
             enddo
             stop
            endif
            voltotal_depart=voltotal_depart+LS_voltotal_depart

            im_refine_density=0

            do im=1,num_materials
             ! density
             dencomp_data=(im-1)*num_state_material+ENUM_DENVAR+1
              ! conserve is initialized in the beginning of this routine.
              ! donate_density is equal to the density that is stored in the
              ! old state variable.
             donate_density= &
              conserve(D_DECL(idonate,jdonate,kdonate), &
                fine_offset+CISLCOMP_STATES+dencomp_data) 
             if (is_compressible_mat(im).eq.0) then
              donate_mom_density= &
               mom_den(D_DECL(idonate,jdonate,kdonate),im) 
             else if (is_compressible_mat(im).eq.1) then
              im_refine_density=im_refine_density+1
              if (fort_im_refine_density_map(im_refine_density).eq.im-1) then
               !do nothing
              else
               print *,"fort_im_refine_density_map invalid"
               stop
              endif
              donate_mom_density=refineden(D_DECL(idonate,jdonate,kdonate), &
                (im_refine_density-1)*ENUM_NUM_REFINE_DENSITY_TYPE+ &
                nfine_stencil)
             else
              print *,"is_compressible_mat(im) invalid"
              stop
             endif

             if (donate_density.gt.zero) then
              ! do nothing
             else
              print *,"donate_density must be positive"
              stop
             endif
             if (donate_mom_density.gt.zero) then
              ! do nothing
             else
              print *,"donate_mom_density must be positive"
              stop
             endif
             massdepart=donate_density
             massdepart_mom=donate_mom_density

             if (massdepart.gt.zero) then
              ! do nothing
             else
              print *,"density invalid in vfrac split 2"
              print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
              print *,"idonate,jdonate,kdonate ",idonate,jdonate,kdonate
              print *,"im ",im
              print *,"donate_density ",donate_density
              print *,"multi_volume_grid(im) ",multi_volume_grid(im)
              print *,"massdepart ",massdepart
              stop
             endif

             if (massdepart_mom.gt.zero) then
              ! do nothing
             else
              print *,"density (for mom) invalid in vfrac split 2"
              print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
              print *,"idonate,jdonate,kdonate ",idonate,jdonate,kdonate
              print *,"im ",im
              print *,"donate_density ",donate_density
              print *,"multi_volume_grid(im) ",multi_volume_grid(im)
              print *,"massdepart_mom ",massdepart_mom
              stop
             endif

             do veldir=1,SDIM
              donate_data= &
               conserve(D_DECL(idonate,jdonate,kdonate), &
                fine_offset+veldir)
              mom2(veldir)=multi_volume_grid(im)*massdepart_mom*donate_data
             enddo  ! veldir=1..sdim (velocity)

             massdepart=massdepart*multi_volume_grid(im)
             massdepart_mom=massdepart_mom*multi_volume_grid(im)

             veldata(CISLCOMP_DEN_MOM+im)= &
              veldata(CISLCOMP_DEN_MOM+im)+massdepart_mom
             veldata(CISLCOMP_STATES+dencomp_data)= &
              veldata(CISLCOMP_STATES+dencomp_data)+massdepart

             if (is_compressible_mat(im).eq.0) then
              !do nothing
             else if (is_compressible_mat(im).eq.1) then
              refine_den_bucket(im_refine_density)= &
               refine_den_bucket(im_refine_density)+massdepart
             else
              print *,"is_compressible_mat(im) invalid"
              stop
             endif

             ! skip density,then do energy,scalars,Q, ...
             ! for temperature:
             ! if incompressible: conserve=den * temp
             ! if compressible  : conserve=0.5 den |u|^2 + den * temp
             ! this is energy in the departure region.
             istate=2
             do while (istate.le.num_state_material)
              statecomp_data=(im-1)*num_state_material+istate

               ! conserve initialized in the beginning of this routine.
               ! Temperature and species variables are multiplied by 
               ! dencore(im) in BUILD_CONSERVE.  (dencore(im) is the
               ! value of density stored in the state variable)
              donate_data= &
               conserve(D_DECL(idonate,jdonate,kdonate), &
                 fine_offset+CISLCOMP_STATES+statecomp_data) 
              
              veldata(CISLCOMP_STATES+statecomp_data)= &
               veldata(CISLCOMP_STATES+statecomp_data)+ & 
               multi_volume_grid(im)*donate_data

              if (istate.eq.ENUM_TEMPERATUREVAR+1) then
               if (veldata(CISLCOMP_STATES+statecomp_data).ge.zero) then
                ! do nothing
               else
                print *,"energy became negative "
                print *,"im,comp2 ",im,CISLCOMP_STATES+statecomp_data
                print *,"current donated value ", &
                 veldata(CISLCOMP_STATES+statecomp_data)
                print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse 
                print *,"idonate,jdonate,kdonate ", &
                 idonate,jdonate,kdonate
                print *,"normdir,dir_counter ",normdir,dir_counter
                print *,"num_materials ",num_materials
                stop
               endif
              endif

              istate=istate+1
             enddo  !do while (istate.le.num_state_material)

             if ((num_materials_viscoelastic.ge.1).and. &
                 (num_materials_viscoelastic.le.num_materials)) then

              if (fort_store_elastic_data(im).eq.1) then
               imap=1
               do while ((fort_im_elastic_map(imap)+1.ne.im).and. &
                         (imap.le.num_materials_viscoelastic))
                imap=imap+1
               enddo

               if (imap.le.num_materials_viscoelastic) then

                do istate=1,ENUM_NUM_TENSOR_TYPE
                 statecomp_data=(imap-1)*ENUM_NUM_TENSOR_TYPE+istate
                 ! configuration tensor is stored at the cell centers, not the
                 ! corresponding material centroid.
                 donate_data= &
                  tensor(D_DECL(idonate,jdonate,kdonate),statecomp_data)
                 veldata(CISLCOMP_TENSOR+statecomp_data)= &
                  veldata(CISLCOMP_TENSOR+statecomp_data)+ &
                  LS_voltotal_depart*donate_data 
                enddo !istate=1..ENUM_NUM_TENSOR_TYPE

               else 
                print *,"imap invalid"
                stop
               endif
              else if (fort_store_elastic_data(im).eq.0) then
               ! do nothing
              else
               print *,"fort_store_elastic_data(im) invalid"
               stop
              endif

             else if (num_materials_viscoelastic.eq.0) then
              ! do nothing
             else
              print *,"num_materials_viscoelastic invalid:fort_vfrac_split"
              stop
             endif

             ! level set function for im material.
             ! level set function is stored at the cell centers, not the
             ! corresponding material centroid.
             donate_data=LS(D_DECL(idonate,jdonate,kdonate),im) 
             veldata(CISLCOMP_LS+im)=veldata(CISLCOMP_LS+im)+ &
              LS_voltotal_depart*donate_data

             vofcomp=(im-1)*ngeom_raw+1
             ! material volume from departure (donating) region
             veldata(CISLCOMP_MOF+vofcomp)= &
              veldata(CISLCOMP_MOF+vofcomp)+multi_volume_grid(im)
             ! material volume from target (accepting) region
             veldata(CISLCOMP_FTARGET+im)= &
              veldata(CISLCOMP_FTARGET+im)+multi_volume(im)

             if (is_compressible_mat(im).eq.0) then
              !do nothing
             else if (is_compressible_mat(im).eq.1) then
              refine_vol_bucket(im_refine_density)= &
               refine_vol_bucket(im_refine_density)+multi_volume(im)
             else
              print *,"is_compressible_mat(im) invalid"
              stop
             endif


             ! material centroid from target (accepting) region
             do dir2=1,SDIM
              veldata(CISLCOMP_MOF+vofcomp+dir2)= &
               veldata(CISLCOMP_MOF+vofcomp+dir2)+ &
               multi_volume(im)*multi_cen(dir2,im)
             enddo 

             do veldir=1,SDIM
               ! fluid materials tessellate the domain.
              if (is_rigid(im).eq.0) then
               veldata(veldir)=veldata(veldir)+mom2(veldir) 

               if (veldir.eq.1) then
                xmomside(D_DECL(icrse,jcrse,kcrse),ifine+1)= &
                 xmomside(D_DECL(icrse,jcrse,kcrse),ifine+1)+ &
                 mom2(veldir)
                xmassside(D_DECL(icrse,jcrse,kcrse),ifine+1)= &
                 xmassside(D_DECL(icrse,jcrse,kcrse),ifine+1)+ &
                 massdepart_mom
               else if (veldir.eq.2) then
                ymomside(D_DECL(icrse,jcrse,kcrse),jfine+1)= &
                 ymomside(D_DECL(icrse,jcrse,kcrse),jfine+1)+ &
                 mom2(veldir)
                ymassside(D_DECL(icrse,jcrse,kcrse),jfine+1)= &
                 ymassside(D_DECL(icrse,jcrse,kcrse),jfine+1)+ &
                 massdepart_mom
               else if ((veldir.eq.3).and.(SDIM.eq.3)) then
                zmomside(D_DECL(icrse,jcrse,kcrse),kfine+1)= &
                 zmomside(D_DECL(icrse,jcrse,kcrse),kfine+1)+ &
                 mom2(veldir)
                zmassside(D_DECL(icrse,jcrse,kcrse),kfine+1)= &
                 zmassside(D_DECL(icrse,jcrse,kcrse),kfine+1)+ &
                 massdepart_mom
               else
                print *,"veldir invalid: ",veldir
                stop
               endif

              else if (is_rigid(im).eq.1) then
               ! do nothing
              else
               print *,"is_rigid invalid GODUNOV_3D.F90"
               stop
              endif
             enddo ! veldir=1..sdim
    
            enddo ! im=1,..,num_materials (state variables, geometry, velocity)

           endif ! volint>0

          else if (check_intersection.eq.0) then
           ! do nothing
          else
           print *,"check_intersection invalid"
           stop
          endif

         enddo  ! ifine_stencil
         enddo  ! jfine_stencil
#if (AMREX_SPACEDIM==3)
         enddo  ! kfine_stencil
#endif

        enddo !istencil=idonatelow,idonaatehigh

        im_refine_density=0
        do im=1,num_materials
         if (is_compressible_mat(im).eq.0) then
          !do nothing
         else if (is_compressible_mat(im).eq.1) then
          im_refine_density=im_refine_density+1
          if (fort_im_refine_density_map(im_refine_density).eq.im-1) then
           !do nothing
          else
           print *,"fort_im_refine_density_map invalid"
           stop
          endif
          if (refine_vol_bucket(im_refine_density).gt.zero) then
           refinedennew(D_DECL(icrse,jcrse,kcrse), &
            (im_refine_density-1)*ENUM_NUM_REFINE_DENSITY_TYPE+nfine)= &
             refine_den_bucket(im_refine_density)/ &
             refine_vol_bucket(im_refine_density)
          else if (refine_vol_bucket(im_refine_density).eq.zero) then
           refinedennew(D_DECL(icrse,jcrse,kcrse), &
            (im_refine_density-1)*ENUM_NUM_REFINE_DENSITY_TYPE+nfine)= &
             fort_denconst(im) 
          else
           print *,"refine_vol_bucket invalid: ",im_refine_density, &
             refine_vol_bucket(im_refine_density)
           stop
          endif
         else
          print *,"is_compressible_mat(im) invalid"
          stop
         endif
        enddo ! im=1..num_materials

       enddo  ! ifine
       enddo  ! jfine
#if (AMREX_SPACEDIM==3)
       enddo  ! kfine
#endif

       voltotal_depart=zero
       voltotal_target=zero
       do im=1,num_materials
        vofcomp=(im-1)*ngeom_raw+1 
        volmat_target(im)=veldata(CISLCOMP_FTARGET+im)
        volmat_depart(im)=veldata(CISLCOMP_MOF+vofcomp)

        volmat_target_cor(im)=volmat_target(im)
        volmat_depart_cor(im)=volmat_depart(im)

         ! fluid materials tessellate the domain.
        if (is_rigid(im).eq.0) then
         voltotal_target=voltotal_target+volmat_target(im)
         voltotal_depart=voltotal_depart+volmat_depart(im)
        else if (is_rigid(im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid GODUNOV_3D.F90"
         stop
        endif
       enddo ! im=1..num_materials
       
       if (voltotal_depart.gt.zero) then
        ! do nothing
       else
        print *,"voltotal_depart bust "
        stop
       endif
       if (voltotal_target.gt.zero) then
        ! do nothing
       else
        print *,"voltotal_target bust "
        stop
       endif

       do im=1,num_materials
        vofcomp=(im-1)*ngeom_raw+1 

        newvfrac(im)=volmat_target(im)/voltotal_target
        newvfrac_cor(im)=volmat_target_cor(im)/voltotal_target

        if (vof0(D_DECL(icrse,jcrse,kcrse),im).le.half) then
         newvfrac_weymouth(im)=volmat_depart_cor(im)/voltotal_target
         if (newvfrac_weymouth(im).gt.one) then
          newvfrac_weymouth(im)=one
         endif
        else if (vof0(D_DECL(icrse,jcrse,kcrse),im).ge.half) then
         newvfrac_weymouth(im)=one- &
          (voltotal_depart-volmat_depart_cor(im))/voltotal_target
         if (newvfrac_weymouth(im).lt.zero) then
          newvfrac_weymouth(im)=zero
         endif
        else
         print *,"vof0 bust"
         stop
        endif

        if (newvfrac(im).le.VOFTOL) then
         newvfrac_weymouth(im)=newvfrac(im)
         newvfrac_cor(im)=newvfrac(im)
        endif
     
        call CISBOX(xsten_accept,nhalf, &
         xlo,dx,icrse,jcrse,kcrse, &
         bfact,level, &
         volcell_accept,cencell_accept,SDIM)
  
        do dir2=1,SDIM
         if (newvfrac(im).gt.VOFTOL) then
          newcen(dir2,im)= &
           veldata(CISLCOMP_MOF+vofcomp+dir2)/ &
           volmat_target(im)- &
           cencell_accept(dir2)
         else
          newcen(dir2,im)=zero
         endif
        enddo ! dir2

       enddo  ! im=1..num_materials (geometry)

       call consistent_materials(newvfrac_cor,newcen)

       if ((EILE_flag.eq.1).or. & ! EILE
           (EILE_flag.eq.2).or. & ! EI
           (EILE_flag.eq.3)) then ! LE
        ! do nothing
       else if (EILE_flag.eq.-1) then ! weymouth and Yue
        do im=1,num_materials
         newvfrac_cor(im)=newvfrac_weymouth(im)
        enddo
        call consistent_materials(newvfrac_cor,newcen)
       else
        print *,"EILE_flag invalid"
        stop
       endif
   
       ! pressure
       statecomp_data=STATECOMP_PRES+1

       if (divu_outer_sweeps.eq.0) then
        snew_hold(statecomp_data)= &
          velfab(D_DECL(icrse,jcrse,kcrse),statecomp_data)
       else if ((divu_outer_sweeps.ge.1).and. &
                (divu_outer_sweeps.lt.num_divu_outer_sweeps)) then
        snew_hold(statecomp_data)= &
          snew(D_DECL(icrse,jcrse,kcrse),statecomp_data)
       else
        print *,"divu_outer_sweeps invalid: ",divu_outer_sweeps
        stop
       endif

       ! density
       do im=1,num_materials

        dencomp_data=(im-1)*num_state_material+ENUM_DENVAR+1
        massdepart=veldata(CISLCOMP_STATES+dencomp_data)
        if (massdepart.ge.zero) then
         ! do nothing
        else
         print *,"new mass cannot be negative"
         print *,"im= ",im
         print *,"new mass= ",massdepart
         stop
        endif
         ! if is_rigid==1 or voldepart<eps or voltarget<eps then
         !  den=fort_denconst(im)
         ! else if mat_type==0 then
         !  if override==0 or 2 then
         !   den=fort_denconst(im)
         !  else if override==1 then
         !   den=massdepart/voldepart
         !  endif
         ! else if mat_type>0 then
         !  den=massdepart/voltarget
         ! endif
        if (all_incomp.eq.1) then
         vol_target_local=volmat_depart_cor(im)
        else if (all_incomp.eq.0) then
         vol_target_local=volmat_target_cor(im)
        else
         print *,"all_incomp invalid: ",all_incomp
         stop
        endif

        if ((is_compressible_mat(im).eq.0).or. &
            (fort_material_conservation_form(im).eq.0)) then
         vol_target_local=volmat_depart_cor(im)
        else if ((is_compressible_mat(im).eq.1).and. &
                 (fort_material_conservation_form(im).eq.1)) then
         if (incompressible_interface_flag.eq.0) then
          ! do nothing
         else if (incompressible_interface_flag.eq.1) then
          vol_target_local=volmat_depart_cor(im)
         else 
          print *,"incompressible_interface_flag invalid"
          stop
         endif
        else
         print *, &
           "is_compressible_mat or fort_material_conservation_form invalid"
         print *,"im,is_compressible_mat: ", &
              im,is_compressible_mat(im)
         print *,"im,fort_material_conservation_form: ", &
              im,fort_material_conservation_form(im)
         stop
        endif

        ! if is_rigid(im), density=fort_denconst(im)
        ! if incompressible,
        !   if constant_density_all_time==1 then density=fort_denconst(im)
        !   if constant_density_all_time==0 then 
        !                                  density=mass_depart/vol_depart
        ! if compressible,
        !  if constant_density_all_time==0 then
        !   density=massdepart/voltarget
        !  else
        !   return error.
        ! subroutine derive_density declared in GODUNOV_3D.F90 (this file)
        call derive_density(volmat_depart_cor(im), &
         vol_target_local,voltotal_depart, &
         constant_density_all_time, &
         massdepart, & !intent(in)
         im, &
         dencore(im)) ! intent(out)

        istate=STATECOMP_STATES+(im-1)*num_state_material+ENUM_DENVAR+1
        if (dencore(im).gt.zero) then
         ! do nothing
        else
         print *,"density must be positive vfrac_split 2"
         print *,"im,dencore(im) ",im,dencore(im)
         stop
        endif
        if (dencore(im).lt.density_floor(im)) then
         dencore(im)=density_floor(im)
        endif
        if (density_ceiling(im).gt.zero) then
         if (dencore(im).gt.density_ceiling(im)) then
          dencore(im)=density_ceiling(im)
         endif
        else
         print *,"density_ceiling(im) invalid"
         stop
        endif
        snew_hold(istate)=dencore(im)

       enddo ! im, updating density

       ! levelset function
       ! voltotal_depart=sum_{fluid mat} volmat_depart(im)
       do im=1,num_materials
        if (voltotal_depart.gt.zero) then
         newLS(im)=veldata(CISLCOMP_LS+im)/voltotal_depart
        else
         print *,"voltotal_depart invalid"
         stop
        endif 
       enddo  ! im=1..num_materials (updating levelset vars)

       do im=1,num_materials

        if ((num_materials_viscoelastic.ge.1).and. &
            (num_materials_viscoelastic.le.num_materials)) then

         if (fort_store_elastic_data(im).eq.1) then
          imap=1
          do while ((fort_im_elastic_map(imap)+1.ne.im).and. &
                    (imap.le.num_materials_viscoelastic))
           imap=imap+1
          enddo

          if (imap.le.num_materials_viscoelastic) then

           do istate=1,ENUM_NUM_TENSOR_TYPE
            statecomp_data=(imap-1)*ENUM_NUM_TENSOR_TYPE+istate
            if (voltotal_depart.gt.zero) then
             tennew_hold(statecomp_data)= &
               veldata(CISLCOMP_TENSOR+statecomp_data)/voltotal_depart
            else
             print *,"voltotal_depart invalid"
             stop
            endif 
 
           enddo !istate=1..ENUM_NUM_TENSOR_TYPE

          else 
           print *,"imap invalid"
           stop
          endif
         else if (fort_store_elastic_data(im).eq.0) then
          ! do nothing
         else
          print *,"fort_store_elastic_data(im) invalid"
          stop
         endif

        else if (num_materials_viscoelastic.eq.0) then
         ! do nothing
        else
         print *,"num_materials_viscoelastic invalid:fort_vfrac_split"
         stop
        endif
 
       enddo ! im=1..num_materials (updating viscoelastic vars)

       ! velocity=mom/mass
       ! fluid materials tessellate the domain.
       totalmass_depart=zero
       do im=1,num_materials
        if (is_rigid(im).eq.0) then
         massdepart_mom=veldata(CISLCOMP_DEN_MOM+im)
         totalmass_depart=totalmass_depart+massdepart_mom
        else if (is_rigid(im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid GODUNOV_3D.F90"
         stop
        endif
       enddo ! im=1..num_materials

       if (totalmass_depart.gt.zero) then
        ! do nothing
       else
        print *,"totalmass_depart bust totalmass_depart=",totalmass_depart
        do dir2=1,SDIM
         print *,"dir,fablo,fabhi ",dir2,fablo(dir2),fabhi(dir2)
        enddo
        print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
        print *,"num_state_material=",num_state_material
        print *,"normdir=",normdir
        print *,"dir_counter=",dir_counter
        print *,"num_materials,map_forward,level,finest_level ", &
         num_materials,map_forward,level,finest_level
        stop
       endif

       do veldir=1,SDIM
        snew_hold(veldir)=veldata(veldir)/totalmass_depart
       enddo

        ! make sure 0<=F<=1 and sum F_i = 1.
        ! also truncation 1.0e-8 to 0 and 1-1.0e-8 to 1.
       call consistent_materials(newvfrac_cor,newcen)

       do im=1,num_materials

        vofcomp=(im-1)*ngeom_raw+1

        KE=zero
        do veldir=1,SDIM
         KE=KE+snew_hold(veldir)**2
        enddo ! veldir
        KE=half*KE

        if (ngeom_raw.eq.SDIM+1) then
         snew_hold(STATECOMP_MOF+vofcomp)=newvfrac_cor(im)
         do dir2=1,SDIM
          snew_hold(STATECOMP_MOF+vofcomp+dir2)=newcen(dir2,im)
         enddo
        else
         print *,"ngeom_raw invalid in vfrac split"
         print *,"ngeom_raw= ",ngeom_raw
         stop
        endif

        no_material_flag=0

        if ( (volmat_depart(im).le. &
              VOFTOL*voltotal_depart).or. &
             (volmat_depart_cor(im).le. &
              VOFTOL*voltotal_depart).or. &
             (newvfrac_cor(im).le.VOFTOL).or. &
             (volmat_target(im).le. &
              VOFTOL*voltotal_depart).or. &
             (volmat_target_cor(im).le. &
              VOFTOL*voltotal_depart) ) then

         no_material_flag=1

        endif

         ! in: fort_vfrac_split
        dencomp_data=(im-1)*num_state_material+ENUM_DENVAR+1

        istate=1
        do while (istate.le.num_state_material)

         if (istate.eq.ENUM_DENVAR+1) then
          ! do nothing, density updated above
          istate=istate+1
         else if (istate.eq.ENUM_TEMPERATUREVAR+1) then 

          do ispecies=1,num_species_var
           speccomp_data=(im-1)*num_state_material+num_state_base+ &
             ispecies
           if (no_material_flag.eq.1) then ! no material (im)
            snew_hold(STATECOMP_STATES+speccomp_data)=zero
           else if (no_material_flag.eq.0) then
            if (is_rigid(im).eq.1) then ! mass fraction=0 in solids.
             snew_hold(STATECOMP_STATES+speccomp_data)=zero
            else if (is_rigid(im).eq.0) then
             massdepart=veldata(CISLCOMP_STATES+dencomp_data)
             if (massdepart.gt.zero) then
              snew_hold(STATECOMP_STATES+speccomp_data)= &
               veldata(CISLCOMP_STATES+speccomp_data)/massdepart
             else
              print *,"massdepart invalid"
              stop
             endif 
            else
             print *,"is_rigid invalid GODUNOV_3D.F90"
             stop
            endif
           else 
            print *,"no_material_flag invalid"
            stop
           endif

          enddo ! ispecies=1..num_species_var

          tempcomp_data=(im-1)*num_state_material+ENUM_TEMPERATUREVAR+1

          if (no_material_flag.eq.1) then
           snew_hold(STATECOMP_STATES+tempcomp_data)=fort_tempconst(im)
          else if (no_material_flag.eq.0) then
           if (is_rigid(im).eq.1) then
            if (fort_material_type(im).ne.999) then
             print *,"fort_material_type(im).ne.999"
             stop
            endif

            ! solidheat_flag==0 diffuse in solid
            ! solidheat_flag==1 dirichlet bc at solid-fluid
            ! solidheat_flag==2 insulating bc at solid-fluid
            if (solidheat_flag.eq.0) then ! diffuse in solid

             massdepart=veldata(CISLCOMP_STATES+dencomp_data)
             if (massdepart.gt.zero) then
              !do nothing
             else
              print *,"massdepart invalid: ",massdepart
              stop
             endif 
             ETcore=veldata(CISLCOMP_STATES+tempcomp_data)/massdepart

            else if (solidheat_flag.eq.2) then ! neumann

             ! placeholder
             ETcore=fort_tempconst(im)

            else if (solidheat_flag.eq.1) then ! dirichlet

             ! placeholder
             ETcore=fort_tempconst(im)

            else
             print *,"solidheat_flag invalid: ",solidheat_flag
             stop
            endif

            if (ETcore.lt.fort_tempcutoff(im)) then
             ETcore=fort_tempcutoff(im)
            endif
            if (ETcore.gt.fort_tempcutoffmax(im)) then
             ETcore=fort_tempcutoffmax(im)
            endif

            if (ETcore.gt.zero) then
             ! do nothing
            else
             print *,"Energy (ETcore) went negative: ",ETcore
             stop
            endif

            snew_hold(STATECOMP_STATES+tempcomp_data)=ETcore
           else if (is_rigid(im).eq.0) then
            if ((fort_material_type(im).ge.0).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
             ! do nothing
            else
             print *,"fort_material_type invalid"
             stop
            endif
            massdepart=veldata(CISLCOMP_STATES+dencomp_data)
            if (massdepart.gt.zero) then
             ! do nothing
            else
             print *,"massdepart invalid: ",massdepart
             stop
            endif 
            ! integral_omega_depart rho T F_m /
            ! integral_omega_depart rho F_m
            if ((is_compressible_mat(im).eq.0).or. &
                (fort_material_conservation_form(im).eq.0)) then
             ETcore=veldata(CISLCOMP_STATES+tempcomp_data)/massdepart
            else if ((is_compressible_mat(im).eq.1).and. &
                     (fort_material_conservation_form(im).eq.1)) then
             ! integral_omega_depart rho (u dot u/2 + c_v T) F_m /
             ! integral_omega_depart rho F_m
             ETcore=veldata(CISLCOMP_STATES+tempcomp_data)/massdepart

             local_internal=ETcore-KE
             if (local_internal.gt.zero) then

              call init_massfrac_parm(dencore(im),massfrac_parm,im)
              do ispecies=1,num_species_var
               speccomp_data=(im-1)*num_state_material+num_state_base+ &
                 ispecies
               massfrac_parm(ispecies)= &
                 snew_hold(STATECOMP_STATES+speccomp_data)
               if (massfrac_parm(ispecies).ge.zero) then
                ! do nothing
               else
                print *,"massfrac_parm(ispecies) invalid"
                stop
               endif
              enddo ! ispecies=1..num_species_var

              call TEMPERATURE_material(dencore(im),massfrac_parm, &
                ETcore,local_internal,fort_material_type(im),im) 
             else if (local_internal.le.zero) then
              ETcore=fort_tempcutoff(im)
             else
              print *,"local_internal invalid:",im,local_internal
              stop
             endif
            else
             print *, &
              "is_compressible_mat or fort_material_conservation_form invalid"
             stop
            endif
            if (ETcore.lt.fort_tempcutoff(im)) then
             ETcore=fort_tempcutoff(im)
            endif
            if (ETcore.gt.fort_tempcutoffmax(im)) then
             ETcore=fort_tempcutoffmax(im)
            endif
            if (ETcore.gt.zero) then
             ! do nothing
            else
             print *,"Energy (ETcore) went negative(2)"
             stop
            endif
            snew_hold(STATECOMP_STATES+tempcomp_data)=ETcore
           else
            print *,"is_rigid invalid GODUNOV_3D.F90"
            stop
           endif
          else 
           print *,"no_material_flag invalid"
           stop
          endif

          istate=istate+1+num_species_var
         else
          print *,"istate invalid"
          stop
         endif

        enddo ! do while (istate.le.num_state_material)

       enddo  ! im=1..num_materials

       do istate=1,STATECOMP_STATES
        snew(D_DECL(icrse,jcrse,kcrse),istate)=snew_hold(istate)
       enddo

       if (stokes_flow.eq.1) then
        wt_oldvel=one
       else if (stokes_flow.eq.0) then
        wt_oldvel=zero
       else
        print *,"stokes_flow invalid"
        stop
       endif

       do istate=1,SDIM
        snew(D_DECL(icrse,jcrse,kcrse),istate)= &
           (one-wt_oldvel)*snew_hold(istate)+ &
           wt_oldvel*velfab(D_DECL(icrse,jcrse,kcrse),istate)
       enddo

        ! density, temperature, other scalars
        ! volume fractions, centroids
       do im=1,num_materials

        if (is_rigid(im).eq.0) then

         do istate=1,num_state_material
          statecomp_data=STATECOMP_STATES+(im-1)*num_state_material+istate
          snew(D_DECL(icrse,jcrse,kcrse),statecomp_data)= &
            snew_hold(statecomp_data)
         enddo ! istate=1..num_state_material

         do igeom=1,ngeom_raw
          statecomp_data=STATECOMP_STATES+ &
            num_materials*num_state_material+ &
            (im-1)*ngeom_raw+igeom
          snew(D_DECL(icrse,jcrse,kcrse),statecomp_data)= &
            snew_hold(statecomp_data)
         enddo

        else if (is_rigid(im).eq.1) then

         if (solidheat_flag.eq.0) then ! diffuse in solid
          tempcomp_data=STATECOMP_STATES+(im-1)*num_state_material+ &
            ENUM_TEMPERATUREVAR+1
          snew(D_DECL(icrse,jcrse,kcrse),tempcomp_data)= &
            snew_hold(tempcomp_data)
         else if (solidheat_flag.eq.2) then ! neumann
          ! do nothing
         else if (solidheat_flag.eq.1) then ! dirichlet
          ! do nothing
         else
          print *,"solidheat_flag invalid: ",solidheat_flag
          stop
         endif

        else
         print *,"is_rigid invalid GODUNOV_3D.F90: ",im,is_rigid(im)
         stop
        endif

        if ((num_materials_viscoelastic.ge.1).and. &
            (num_materials_viscoelastic.le.num_materials)) then

         if (fort_store_elastic_data(im).eq.1) then
          imap=1
          do while ((fort_im_elastic_map(imap)+1.ne.im).and. &
                    (imap.le.num_materials_viscoelastic))
           imap=imap+1
          enddo
          if (imap.le.num_materials_viscoelastic) then

           do istate=1,ENUM_NUM_TENSOR_TYPE
            statecomp_data=(imap-1)*ENUM_NUM_TENSOR_TYPE+istate
            tennew(D_DECL(icrse,jcrse,kcrse),statecomp_data)= &
             tennew_hold(statecomp_data)
           enddo !istate=1..ENUM_NUM_TENSOR_TYPE

          else 
           print *,"imap invalid"
           stop
          endif
         else if (fort_store_elastic_data(im).eq.0) then
          ! do nothing
         else
          print *,"fort_store_elastic_data(im) invalid"
          stop
         endif

        else if (num_materials_viscoelastic.eq.0) then
         ! do nothing
        else
         print *,"num_materials_viscoelastic invalid:fort_vfrac_split"
         stop
        endif

       enddo ! im=1..num_materials

       do im=1,num_materials
        if (is_rigid(im).eq.0) then
         LSnew(D_DECL(icrse,jcrse,kcrse),im)=newLS(im)

         ! level set comes from Lagrangian representation.
        else if (is_rigid(im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(im) invalid"
         stop
        endif
       enddo ! im=1..num_materials

      enddo
      enddo
      enddo ! icrse,jcrse,kcrse -> growntilebox(0 ghost cells)

      do veldir=1,SDIM

        iii=0
        jjj=0
        kkk=0 

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
          growlo,growhi,0,veldir-1)

        if (veldir.eq.1) then
         iii=1
        else if (veldir.eq.2) then
         jjj=1
        else if ((veldir.eq.3).and.(SDIM.eq.3)) then
         kkk=1
        else
         print *,"veldir invalid"
         stop
        endif

        do kcrse=growlo(3),growhi(3)
        do jcrse=growlo(2),growhi(2)
        do icrse=growlo(1),growhi(1)

          ! veldir=1..sdim
         call gridstenMAC_level(xsten_MAC, &
            icrse,jcrse,kcrse, &
            level,nhalf,veldir-1)

         do dir_local=1,SDIM
          xclamped(dir_local)=xsten_MAC(0,dir_local)
         enddo

         zapvel=0
         if (levelrz.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if (levelrz.eq.COORDSYS_RZ) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          if ((xsten_MAC(0,1).le.EPS2*dx(1)).and. &
              (veldir.eq.1)) then
           zapvel=1
           if (icrse.eq.0) then
            !do nothing
           else
            print *,"icrse invalid"
            stop
           endif 
          endif
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
          if ((xsten_MAC(0,1).le.EPS2*dx(1)).and. &
              (veldir.eq.1)) then
           zapvel=1
           if (icrse.eq.0) then
            !do nothing
           else
            print *,"icrse invalid"
            stop
           endif 
          endif
         else
          print *,"levelrz invalid"
          stop
         endif

         if (zapvel.eq.1) then

          if (veldir.eq.1) then
           xmac_new(D_DECL(icrse,jcrse,kcrse))=zero
           wt_oldvel=zero
          else
           print *,"veldir invalid"
           stop
          endif

         else if (zapvel.eq.0) then

          iright=icrse
          jright=jcrse
          kright=kcrse
          ileft=iright-iii
          jleft=jright-jjj
          kleft=kright-kkk
           ! mask=1 if cell is not covered by level+1 or cell is 
           ! outside the domain.
          maskleft=mask(D_DECL(ileft,jleft,kleft))
          maskright=mask(D_DECL(iright,jright,kright))

          if ((maskleft.eq.one).and.(maskright.eq.one)) then

            momface_total=zero
            massface_total=zero

             ! iside=-1 (left of face)
             ! iside=1  (right of face)
            do iside=-1,1,2

              ! iside=-1 (left of face)
             if (iside.eq.-1) then
              icell=ileft
              jcell=jleft
              kcell=kleft
              ibucket=2  ! right side of cell ileft,jleft,kleft

               ! iside=1  (right of face)
             else if (iside.eq.1) then
              icell=iright
              jcell=jright
              kcell=kright
              ibucket=1  ! left side of cell iright,jright,kright
             else 
              print *,"iside invalid"
              stop
             endif

             if (veldir.eq.1) then
              massquarter=xmassside(D_DECL(icell,jcell,kcell),ibucket)
              momquarter=xmomside(D_DECL(icell,jcell,kcell),ibucket)
             else if (veldir.eq.2) then
              massquarter=ymassside(D_DECL(icell,jcell,kcell),ibucket)
              momquarter=ymomside(D_DECL(icell,jcell,kcell),ibucket)
             else if ((veldir.eq.3).and.(SDIM.eq.3)) then
              massquarter=zmassside(D_DECL(icell,jcell,kcell),ibucket)
              momquarter=zmomside(D_DECL(icell,jcell,kcell),ibucket)
             else
              print *,"veldir invalid"
              stop
             endif
             if (massquarter.ge.zero) then
              ! do nothing
             else
              print *,"massquarter cannot be negative"
              print *,"icell,jcell,kcell,ibucket ", &
                icell,jcell,kcell,ibucket 
              print *,"icrse,jcrse,kcrse ",icrse,jcrse,kcrse
              print *,"massquarter ",massquarter
              print *,"momquarter ",momquarter
              stop
             endif

             massface_total=massface_total+massquarter
             momface_total=momface_total+momquarter

            enddo ! iside: do iside=-1,1,2

            if (stokes_flow.eq.1) then
             wt_oldvel=one
            else if (stokes_flow.eq.0) then
             wt_oldvel=zero
            else
             print *,"stokes_flow invalid"
             stop
            endif

            if (massface_total.gt.zero) then

             momface_total=momface_total/massface_total
             ! LS>0 if clamped
             call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
              vel_clamped,temperature_clamped,prescribed_flag,dx)
             if (LS_clamped.ge.zero) then
              momface_total=vel_clamped(veldir)
              wt_oldvel=zero
             else if (LS_clamped.lt.zero) then
              ! do nothing
             else
              print *,"LS_clamped is NaN"
              stop
             endif

            else
             print *,"massface_total invalid"
             stop
            endif
          
            if (veldir.eq.1) then
             xmac_new(D_DECL(icrse,jcrse,kcrse))= &
              (one-wt_oldvel)*momface_total+ &
              wt_oldvel*xmac_old(D_DECL(icrse,jcrse,kcrse))
            else if (veldir.eq.2) then
             ymac_new(D_DECL(icrse,jcrse,kcrse))= &
              (one-wt_oldvel)*momface_total+ &
              wt_oldvel*ymac_old(D_DECL(icrse,jcrse,kcrse))
            else if ((veldir.eq.3).and.(SDIM.eq.3)) then
             zmac_new(D_DECL(icrse,jcrse,kcrse))= &
              (one-wt_oldvel)*momface_total+ &
              wt_oldvel*zmac_old(D_DECL(icrse,jcrse,kcrse))
            else
             print *,"veldir invalid"
             stop
            endif

          else if ((maskleft.eq.zero).or.(maskright.eq.zero)) then
            ! do nothing
          else 
           print *,"maskleft or maskright invalid"
           stop
          endif

         else
          print *,"zapvel invalid fort_vfrac_split"
          stop
         endif

        enddo
        enddo
        enddo  ! i,j,k

      enddo ! veldir=1..sdim

      return
      end subroutine fort_vfrac_split


      ! combine_flag==0 (FVM -> GFM)
      ! combine_flag==1 (GFM -> FVM)
      ! combine_flag==2 (combine if vfrac<VOFTOL)
      ! project_option==SOLVETYPE_VISC (cell centered velocity)
      ! project_option==SOLVETYPE_PRES (MAC velocity - COMBINEVELFACE is called)
      subroutine fort_combinevel( &
       tid, &
       hflag, &
       num_materials_combine, &
       mass_fraction_id, &
       freezing_model, &
       Tanasawa_or_Schrage_or_Kassemi, &
       distribute_from_target, &
       saturation_temp, &
       hydrate_flag, & ! scalar
       nparts, &
       nparts_def, &
       im_solid_map, &
       nsolve, &
       project_option, &
       combine_idx, &
       combine_flag, &
       interface_cond_avail, &
       nstate_main, &
       ncomp_cell, &
       scomp, &
       ncomp, &
       scomp_size, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       ntsat, &
       TgammaFAB,DIMS(TgammaFAB), &
       maskcov,DIMS(maskcov), &
       solxfab,DIMS(solxfab), &
       solyfab,DIMS(solyfab), &
       solzfab,DIMS(solzfab), &
       LSNEW,DIMS(LSNEW), &
       LS,DIMS(LS), &
       vof,DIMS(vof), &
       cellfab,DIMS(cellfab), &
       newcell,DIMS(newcell), &
       state,DIMS(state), &  !Snew
       velbc, &
       listbc, &
       xlo,dx, &
       cur_time) &
      bind(c,name='fort_combinevel')
      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use mass_transfer_module
 
      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: num_materials_combine
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def)
      integer, INTENT(in) :: hflag
      integer, INTENT(in) :: mass_fraction_id(2*num_interfaces)
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      integer, INTENT(in) :: Tanasawa_or_Schrage_or_Kassemi(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
      integer, INTENT(in) :: hydrate_flag
      integer, INTENT(in) :: nsolve
      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: combine_idx
      integer, INTENT(in) :: combine_flag
      integer, INTENT(in) :: interface_cond_avail
      integer, INTENT(in) :: nstate_main
      integer, INTENT(in) :: ncomp_cell
      integer, INTENT(in) :: scomp_size
      integer, INTENT(in) :: scomp(scomp_size)
      integer, INTENT(in) :: ncomp(scomp_size)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: ntsat

      integer, INTENT(in) :: DIMDEC(TgammaFAB)
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(solxfab)
      integer, INTENT(in) :: DIMDEC(solyfab)
      integer, INTENT(in) :: DIMDEC(solzfab)
      integer, INTENT(in) :: DIMDEC(LSNEW)
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(vof)
      integer, INTENT(in) :: DIMDEC(cellfab)
      integer, INTENT(in) :: DIMDEC(newcell)
      integer, INTENT(in) :: DIMDEC(state)

      integer, INTENT(in) :: velbc(SDIM,2,SDIM)
      integer, INTENT(in) :: listbc(SDIM,2, &
             nsolve*num_materials_combine)

      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM) 
      real(amrex_real) :: dxmaxLS

      real(amrex_real), INTENT(in) :: cur_time

      real(amrex_real), INTENT(in),target :: TgammaFAB(DIMV(TgammaFAB),ntsat)
      real(amrex_real), pointer :: TgammaFAB_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: &
              solxfab(DIMV(solxfab),nparts_def*SDIM)
      real(amrex_real), pointer :: solxfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
              solyfab(DIMV(solyfab),nparts_def*SDIM)
      real(amrex_real), pointer :: solyfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
              solzfab(DIMV(solzfab),nparts_def*SDIM)
      real(amrex_real), pointer :: solzfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
              LSNEW(DIMV(LSNEW),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LSNEW_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: LS(DIMV(LS),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
        vof(DIMV(vof),num_materials*ngeom_recon)
      real(amrex_real), pointer :: vof_ptr(D_DECL(:,:,:),:)
       ! output if,
       ! (1) project_option==SOLVETYPE_VISC and combine_flag==2
       ! (2) project_option==SOLVETYPE_HEAT, SOLVETYPE_SPEC, and 
       !     combine_flag==2
      real(amrex_real), INTENT(inout),target :: &
       cellfab(DIMV(cellfab),ncomp_cell)
      real(amrex_real), pointer :: cellfab_ptr(D_DECL(:,:,:),:)
       ! output if,
       !  project_option==SOLVETYPE_HEAT, SOLVETYPE_SPEC, and 
       !     combine_flag==0 or combine_flag==1
      real(amrex_real), INTENT(out),target :: newcell(DIMV(newcell), &
              nsolve*num_materials_combine)
      real(amrex_real), pointer :: newcell_ptr(D_DECL(:,:,:),:)
       !Snew
      real(amrex_real), INTENT(in),target :: state(DIMV(state),nstate_main) 
      real(amrex_real), pointer :: state_ptr(D_DECL(:,:,:),:)

      real(amrex_real) DATA_FLOOR
 
      integer i,j,k
      integer dir
      integer side
      integer i1,j1,k1
      integer k1lo,k1hi
      integer im,im_opp
      integer im_crit
      integer im_primary
      integer im_secondary
      integer ireverse,iten

      integer im_source
      integer im_source_master
      integer im_dest
      integer im_dest_master
      integer tsat_flag

      integer vofcomp
      integer, PARAMETER :: nhalf=3

      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_ofs(-nhalf:nhalf,SDIM)

      real(amrex_real) total_vol_cell
      real(amrex_real) mass_sum
      real(amrex_real) DeDT_total
      real(amrex_real) weight_sum

      real(amrex_real) DeDT
      real(amrex_real) local_volume
      real(amrex_real) local_mass

      real(amrex_real) volcell
      real(amrex_real) cencell(SDIM)

      real(amrex_real) solid_dist

      real(amrex_real) cell_LS(num_materials)
      real(amrex_real) cell_vfrac(num_materials)
      real(amrex_real) cell_species_mfrac(num_materials)
      real(amrex_real) cell_DeDT_mfrac(num_materials)
      real(amrex_real) local_diffusion_coeff

      integer imattype

       !=0 no solid 1<=is_solid_cell<=num_materials+1 otherwise
      integer is_solid_cell 

      integer dencomp
      integer tempcomp
      integer cellcomp

      real(amrex_real) vsol(nsolve)
      real(amrex_real) ucombine(nsolve)

      real(amrex_real) state_mass_average

      real(amrex_real) test_density
      real(amrex_real) test_temperature
      real(amrex_real) massfrac_parm(num_species_var+1)

      real(amrex_real) cell_temperature(num_materials)
      real(amrex_real) new_temperature(num_materials)

      real(amrex_real) LS_source,LS_dest

      real(amrex_real) LL
      integer local_freezing_model
      integer distribute_from_targ
      real(amrex_real) TSAT_master
      real(amrex_real) TDIFF
      real(amrex_real) TDIFF_master
      real(amrex_real) T_out(1)
      real(amrex_real) Tcenter(num_materials)
      real(amrex_real) thermal_state(num_materials)

      real(amrex_real) xtarget(SDIM)
      real(amrex_real) xI(SDIM)
      real(amrex_real) nrm(SDIM)

      real(amrex_real) mofdata(num_materials*ngeom_recon)
      integer nmax

      integer mask_test

      real(amrex_real) :: local_VOF(num_materials)
      integer :: im_local
      integer :: local_vofcomp
      integer im_primary_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      real(amrex_real) VF_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) LS_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) T_sten(D_DECL(-1:1,-1:1,-1:1))

      integer partid
      integer partid_vel_plus
      integer partid_vel
      integer im_solid_vel_plus
      integer im_solid_vel
      integer im_solid_vel_max
      real(amrex_real) LSCRIT_solid_plus
      real(amrex_real) LSCRIT_solid
      real(amrex_real) LSTEST
      integer ncomp_per_tsat
      integer Tgamma_STATUS
      integer ispec
      real(amrex_real) Tgamma
      real(amrex_real) TorYgamma_BC
      integer tsat_comp
      integer local_tessellate
      real(amrex_real) xclamped(SDIM)
      real(amrex_real) LS_clamped
      real(amrex_real) vel_clamped(SDIM)
      real(amrex_real) temperature_clamped
      integer prescribed_flag

      DATA_FLOOR=zero

      nmax=POLYGON_LIST_MAX ! in: fort_combinevel

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if ((hflag.ne.0).and.(hflag.ne.1)) then
       print *,"hflag invalid1  hflag=",hflag
       stop
      endif 
      if (combine_idx.lt.-1) then
       print *,"combine_idx invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid75"
       stop
      endif

      ncomp_per_tsat=EXTRAP_PER_TSAT
      if (ntsat.eq.EXTRAP_NCOMP_TSAT) then
       ! do nothing
      else
       print *,"nstat invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid combine vel"
       stop
      endif
      if ((hydrate_flag.ne.0).and.(hydrate_flag.ne.1)) then
       print *,"hydrate_flag invalid"
       stop
      endif
      if ((nsolve.ne.1).and. &
          (nsolve.ne.STATE_NCOMP_VEL)) then
       print *,"nsolve invalid"
       stop
      endif

      if (nstate_main.ne.STATECOMP_STATES+ &
          num_materials*(num_state_material+ngeom_raw)+1) then
       print *,"nstate_main invalid (1): ",nstate_main
       stop
      endif
      if (nstate_main.ne.STATE_NCOMP) then
       print *,"nstate_main invalid (2): ",nstate_main
       stop
      endif

      if ((ncomp_cell.ne.num_materials_combine).and. &
          (ncomp_cell.ne.nstate_main)) then
       print *,"ncomp_cell invalid"
       stop
      endif
   
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_combinevel: ",nparts
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_combinevel: ",nparts_def
       stop
      endif

      if (combine_idx.eq.-1) then

       if (project_option.eq.SOLVETYPE_PRES) then
        print *,"project_option==SOLVETYPE_PRES not allowed here"
        print *,"combine_idx= ",combine_idx
        stop
       else if (project_option.eq.SOLVETYPE_HEAT) then ! thermal conduction
        if (scomp_size.ne.num_materials) then
         print *,"scomp_size invalid: ",scomp_size
         stop
        endif
        do im=1,num_materials
         if (ncomp(im).ne.1) then
          print *,"ncomp(im) invalid37: ",im,ncomp(im)
          stop
         endif
         if (scomp(im).ne.STATECOMP_STATES+ &
             (im-1)*num_state_material+ENUM_TEMPERATUREVAR) then
          print *,"scomp(im) invalid(1): ",im,scomp(im)
          stop
         endif
        enddo ! im=1..num_materials
       else if (project_option.eq.SOLVETYPE_VISC) then  ! viscosity
        if (scomp(1).ne.STATECOMP_VEL) then
         print *,"scomp(1) invalid: ",scomp(1)
         stop
        endif
        if (ncomp(1).ne.STATE_NCOMP_VEL) then
         print *,"ncomp(1) invalid38: ",ncomp(1)
         stop
        endif
        if (scomp_size.ne.1) then
         print *,"scomp_size invalid(1): ",scomp_size
         stop
        endif
       else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
        if (scomp_size.ne.num_materials) then
         print *,"scomp_size invalid(2): ",scomp_size
         stop
        endif
        do im=1,num_materials
         if (ncomp(im).ne.1) then
          print *,"ncomp(im) invalid39: ",im,ncomp(im)
          stop
         endif
         if (scomp(im).ne.STATECOMP_STATES+ &
             (im-1)*num_state_material+ &
             ENUM_SPECIESVAR+ &
             project_option-SOLVETYPE_SPEC) then
          print *,"scomp(im) invalid(2): ",im,scomp(im)
          stop
         endif
        enddo ! im=1..num_materials
       else
        print *,"project_option invalid fort_combinevel: ",project_option
        stop
       endif

      else if (combine_idx.ge.0) then
       ! do nothing (data in localMF[combine_idx])
      else
       print *,"combine_idx invalid"
       stop
      endif

      if ((interface_cond_avail.eq.0).or. &
          (interface_cond_avail.eq.1)) then
       ! do nothing
      else
       print *,"interface_cond_avail invalid: ",interface_cond_avail
       stop
      endif

      maskcov_ptr=>maskcov
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      solxfab_ptr=>solxfab
      call checkbound_array(fablo,fabhi,solxfab_ptr,0,0)
      solyfab_ptr=>solyfab
      call checkbound_array(fablo,fabhi,solyfab_ptr,0,1)
      solzfab_ptr=>solzfab
      call checkbound_array(fablo,fabhi,solzfab_ptr,0,SDIM-1)
      LSNEW_ptr=>LSNEW
      call checkbound_array(fablo,fabhi,LSNEW_ptr,1,-1)

      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)

      vof_ptr=>vof
      call checkbound_array(fablo,fabhi,vof_ptr,1,-1)
      cellfab_ptr=>cellfab
      call checkbound_array(fablo,fabhi,cellfab_ptr,1,-1)
      newcell_ptr=>newcell
      call checkbound_array(fablo,fabhi,newcell_ptr,0,-1)
      state_ptr=>state
      call checkbound_array(fablo,fabhi,state_ptr,1,-1)
      TgammaFAB_ptr=>TgammaFAB
      call checkbound_array(fablo,fabhi,TgammaFAB_ptr,1,-1)

      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       mask_test=NINT(maskcov(D_DECL(i,j,k)))

       !mask=tag if not covered by level+1 or outside the domain.
       if (mask_test.eq.1) then

        call gridsten_level(xsten,i,j,k,level,nhalf)
        do dir=1,SDIM
         xclamped(dir)=xsten(0,dir)
        enddo

        do im=1,num_materials*ngeom_recon
         mofdata(im)=vof(D_DECL(i,j,k),im)
        enddo
        !  if solid material(s) dominate the cell, then F_solid_raster=1
        !  and F_fluid=0.
        !  if fluid material(s) dominate the cell, then F_solid=0,
        !  sum F_fluid=1
        local_tessellate=3
         !EPS2 tolerance
        call multi_get_volume_tessellate( &
         local_tessellate, &  ! =3
         bfact, &
         dx,xsten, &
         nhalf, & !nhalf=3
         mofdata, &
         geom_xtetlist(1,1,1,tid+1), &
         nmax, &
         nmax, &
         SDIM)

        mass_sum=zero
        total_vol_cell=zero
        DeDT_total=zero

        do im=1,num_materials

         vofcomp=(im-1)*ngeom_recon+1
         local_volume=mofdata(vofcomp)
         
         if ((local_volume.ge.-EPS1).and. &
             (local_volume.le.one+EPS1)) then
          if (local_volume.lt.zero) then
           local_volume=zero
          endif
          if (local_volume.gt.one) then
           local_volume=one
          endif
         else
          print *,"local_volume invalid: ",local_volume
          stop
         endif
         cell_vfrac(im)=local_volume

         if (project_option.eq.SOLVETYPE_HEAT) then
          local_diffusion_coeff=fort_heatviscconst(im)
         else if ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
                  (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
          local_diffusion_coeff= &
           fort_speciesviscconst((project_option-SOLVETYPE_SPEC)* &
              num_materials+im)
         else if (project_option.eq.SOLVETYPE_VISC) then
          local_diffusion_coeff=one
         else
          print *,"project_option invalid: ",project_option
          stop
         endif

         if (local_diffusion_coeff.ge.zero) then 
          if (local_diffusion_coeff.eq.zero) then
           local_diffusion_coeff=EPS5
          else if (local_diffusion_coeff.gt.zero) then
           local_diffusion_coeff=one
          else
           print *,"local_diffusion_coeff invalidA: ",local_diffusion_coeff
           stop
          endif
         else
          print *,"local_diffusion_coeff invalidB: ",local_diffusion_coeff
          stop
         endif

         dencomp=STATECOMP_STATES+(im-1)*num_state_material+ENUM_DENVAR+1
         test_density=state(D_DECL(i,j,k),dencomp)
         if (test_density.gt.zero) then
          ! do nothing
         else
          print *,"test_density invalid: ",test_density
          stop
         endif
         local_mass=test_density*local_volume ! local_volume is a volume frac.
         total_vol_cell=total_vol_cell+local_volume

         mass_sum=mass_sum+local_mass

         cell_species_mfrac(im)=local_mass*local_diffusion_coeff

         imattype=fort_material_type(im)

         test_temperature=state(D_DECL(i,j,k),dencomp+1)
         if (test_temperature.gt.zero) then
          ! do nothing
         else
          print *,"test_temperature invalid: ",test_temperature
          stop
         endif

         call init_massfrac_parm(test_density,massfrac_parm,im)
         do ispec=1,num_species_var
          massfrac_parm(ispec)=state(D_DECL(i,j,k),dencomp+1+ispec)
          if (massfrac_parm(ispec).ge.zero) then
           ! do nothing
          else
           print *,"massfrac_parm(ispec) invalid: ",massfrac_parm(ispec)
           stop
          endif
         enddo ! ispec=1,num_species_var

          !DeDT=cv
          !DeDT_material is declared in GLOBALUTIL.F90
         call DeDT_material(test_density, & !intent(in)
           massfrac_parm, & !intent(in)
           test_temperature, & !intent(in)
           DeDT, & !intent(out)
           imattype,im) !intent(in)
         if (DeDT.gt.zero) then
          ! do nothing
         else
          print *,"DeDT must be positive: ",DeDT
          stop
         endif

         cell_DeDT_mfrac(im)=local_mass*DeDT*local_diffusion_coeff
         DeDT_total=DeDT_total+local_mass*DeDT

        enddo ! im=1,num_materials

        if (total_vol_cell.gt.zero) then
         ! do nothing
        else
         print *,"total_vol_cell invalid: ",total_vol_cell
         stop
        endif

        if (mass_sum.gt.zero) then
         ! do nothing
        else
         print *,"mass_sum invalid: ",mass_sum
         stop
        endif

        if (DeDT_total.gt.zero) then
         ! do nothing
        else
         print *,"DeDT_total invalid: ",DeDT_total
         stop
        endif

        do im=1,num_materials
         cell_vfrac(im)=cell_vfrac(im)/total_vol_cell
         cell_LS(im)=cell_vfrac(im)-half
        enddo ! im

        partid=0
        partid_vel_plus=0
        partid_vel=0
        im_solid_vel_plus=0
        im_solid_vel=0
        LSCRIT_solid_plus=-1.0D+30 !only investigate LSNEW(solid)>=0
        LSCRIT_solid=-1.0D+30      !investigate LSNEW>=0 and LSNEW<0
        
        do im=1,num_materials
         if (is_lag_part(im).eq.1) then
          if (is_rigid(im).eq.1) then
           LSTEST=LSNEW(D_DECL(i,j,k),im)

           if (LSTEST.ge.zero) then
            if (im_solid_vel_plus.eq.0) then
             partid_vel_plus=partid
             im_solid_vel_plus=im
             LSCRIT_solid_plus=LSTEST
            else if ((im_solid_vel_plus.ge.1).and. &
                     (im_solid_vel_plus.le.num_materials)) then
             if (LSTEST.ge.LSCRIT_solid_plus) then
              partid_vel_plus=partid
              im_solid_vel_plus=im
              LSCRIT_solid_plus=LSTEST
             endif
            else
             print *,"im_solid_vel_plus invalid combinevel:",im_solid_vel_plus
             stop
            endif
           else if (LSTEST.lt.zero) then
            ! do nothing
           else
            print *,"LSTEST invalid: ",LSTEST
            stop
           endif

           if ((LSTEST.ge.zero).or. &
               (LSTEST.le.zero)) then

            if (im_solid_vel.eq.0) then
             partid_vel=partid
             im_solid_vel=im
             LSCRIT_solid=LSTEST
            else if ((im_solid_vel.ge.1).and. &
                     (im_solid_vel.le.num_materials)) then
             if (LSTEST.ge.LSCRIT_solid) then
              partid_vel=partid
              im_solid_vel=im
              LSCRIT_solid=LSTEST
             endif
            else
             print *,"im_solid_vel invalid: ",im_solid_vel
             stop
            endif

           else
            print *,"LSTEST invalid: ",LSTEST
            stop
           endif

          else if (is_rigid(im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(im) invalid"
           stop
          endif
          partid=partid+1
         else if (is_lag_part(im).eq.0) then
          if (is_rigid(im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(im) invalid: ",im,is_rigid(im)
           stop
          endif
         else
          print *,"is_lag_part(im) invalid: ",im,is_lag_part(im)
          stop
         endif
        enddo ! im=1..num_materials

        if (partid.ne.nparts) then
         print *,"partid invalid: ",partid,nparts
         stop
        endif

        is_solid_cell=0

        if ((im_solid_vel_plus.ge.1).and. &
            (im_solid_vel_plus.le.num_materials)) then
         if (is_rigid(im_solid_vel_plus).eq.1) then
          is_solid_cell=im_solid_vel_plus
          if (im_solid_map(partid_vel_plus+1)+1.ne.im_solid_vel_plus) then
           print *,"im_solid_map(partid_vel_plus+1)+1.ne.im_solid_vel_plus"
           print *,"partid_vel_plus=",partid_vel_plus
           print *,"im_solid_vel_plus=",im_solid_vel_plus
           print *,"im_solid_map(partid_vel_plus+1)+1=",partid_vel_plus, &
             im_solid_map(partid_vel_plus+1)+1
           stop
          endif
         else if (is_rigid(im_solid_vel_plus).eq.0) then
          ! do nothing
         else
          print *,"is_rigid(im_solid_vel_plus) invalid: ", &
           is_rigid(im_solid_vel_plus) 
          print *,"im_solid_vel_plus=",im_solid_vel_plus
          stop
         endif
        else if (im_solid_vel_plus.eq.0) then
         ! do nothing
        else
         print *,"im_solid_vel_plus invalid combinevel:2",im_solid_vel_plus
         stop
        endif

         ! solid_dist>0 in the solid
         ! returns: im_solid_vel_max=max_{is_rigid(im)==1} cell_LS(im)
        call combine_solid_LS(cell_LS,solid_dist,im_solid_vel_max)

        if (solid_dist.ge.zero) then
         if ((im_solid_vel_max.ge.1).and. &
             (im_solid_vel_max.le.num_materials)) then
          if (is_rigid(im_solid_vel_max).eq.1) then
           is_solid_cell=im_solid_vel_max
          else if (is_rigid(im_solid_vel_max).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(im_solid_vel_max) invalid"
           stop
          endif
         else
          print *,"im_solid_vel_max invalid"
          stop
         endif
        else if (solid_dist.lt.zero) then
         ! do nothing
        else
         print *,"solid_dist invalid; fort_combinevel: ",solid_dist
         stop
        endif

         ! first checks the rigid materials for a positive LS; if none
         ! exist, then the subroutine checks the fluid materials.
        call get_primary_material(cell_LS,im_primary)

        if (is_rigid(im_primary).eq.1) then
         is_solid_cell=im_primary
        else if (is_rigid(im_primary).eq.0) then
         ! do nothing
        else
         print *,"is_rigid invalid; fort_combinevel"
         stop
        endif

        if (project_option.eq.SOLVETYPE_VISC) then ! viscosity

         if (combine_flag.eq.2) then !combine if vfrac<VOFTOL

          if (nsolve.ne.STATE_NCOMP_VEL) then
           print *,"nsolve invalid: ",nsolve
           stop
          endif
          if (scomp_size.ne.1) then
           print *,"scomp_size invalid: ",scomp_size
           stop
          endif
          if (scomp(1).ne.STATECOMP_VEL) then
           print *,"scomp invalid: ",scomp(1)
           stop
          endif
          if (ncomp(1).ne.nsolve) then
           print *,"ncomp(1) invalid: ",ncomp(1)
           print *,"nsolve=",nsolve
           stop
          endif

          do cellcomp=1,SDIM
           if (hflag.eq.0) then !inhomogeneous solid velocity
            if (im_solid_vel.eq.0) then
             vsol(cellcomp)=zero
            else if ((im_solid_vel.ge.1).and. &
                     (im_solid_vel.le.num_materials)) then
             if (im_solid_map(partid_vel+1)+1.ne.im_solid_vel) then
              print *,"im_solid_map(partid_vel+1)+1.ne.im_solid_vel"
              stop
             endif
             if (cellcomp.eq.1) then
              vsol(cellcomp)=half* &
                (solxfab(D_DECL(i,j,k),partid_vel*SDIM+cellcomp)+ &
                 solxfab(D_DECL(i+1,j,k),partid_vel*SDIM+cellcomp))
             else if (cellcomp.eq.2) then
              vsol(cellcomp)=half* &
                (solyfab(D_DECL(i,j,k),partid_vel*SDIM+cellcomp)+ &
                 solyfab(D_DECL(i,j+1,k),partid_vel*SDIM+cellcomp))
             else if ((cellcomp.eq.SDIM).and.(SDIM.eq.3)) then
              vsol(cellcomp)=half* &
                (solzfab(D_DECL(i,j,k),partid_vel*SDIM+cellcomp)+ &
                 solzfab(D_DECL(i,j,k+1),partid_vel*SDIM+cellcomp))
             else
              print *,"cellcomp invalid"
              stop
             endif

            else
             print *,"im_solid_vel invalid"
             stop
            endif
           else if (hflag.eq.1) then !homogeneous solid velocity
            vsol(cellcomp)=zero
           else
            print *,"hflag invalid2 hflag=",hflag
            stop
           endif
          enddo ! cellcomp=1..sdim
          
           ! LS>0 if clamped
          call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
            vel_clamped,temperature_clamped,prescribed_flag,dx)

          if (LS_clamped.ge.zero) then
           is_solid_cell=num_materials+1
           if (hflag.eq.0) then ! inhomogeneous
            do dir=1,SDIM
             vsol(dir)=vel_clamped(dir)
            enddo
           else if (hflag.eq.1) then !homogeneous
            do dir=1,SDIM
             vsol(dir)=zero
            enddo
           else
            print *,"hflag invalid3 hflag=",hflag
            stop
           endif
          else if (LS_clamped.lt.zero) then
           ! do nothing
          else
           print *,"LS_clamped is NaN"
           stop
          endif

          if ((is_solid_cell.ge.1).and. &
              (is_solid_cell.le.num_materials+1)) then
   
           do cellcomp=1,SDIM
            ucombine(cellcomp)=vsol(cellcomp)
           enddo ! cellcomp

          else if (is_solid_cell.eq.0) then

           do cellcomp=1,SDIM
            ucombine(cellcomp)=cellfab(D_DECL(i,j,k),cellcomp)
           enddo ! cellcomp

          else
           print *,"is_solid_cell invalid: ",is_solid_cell
           stop
          endif

          do cellcomp=1,SDIM
           cellfab(D_DECL(i,j,k),cellcomp)=ucombine(cellcomp)
          enddo ! cellcomp

         else
          print *,"combine_flag!=2: ",combine_flag
          print *,"project_option: ",project_option
          stop
         endif
  
        else if ((project_option.eq.SOLVETYPE_HEAT).or. &   ! temperature
                 ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
                  (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then

         state_mass_average=zero
         weight_sum=zero

         do im_crit=1,num_materials

          if (project_option.eq.SOLVETYPE_HEAT) then
           weight_sum=weight_sum+cell_DeDT_mfrac(im_crit)
          else if ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
                   (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
           weight_sum=weight_sum+cell_species_mfrac(im_crit)
          else
           print *,"project_option invalid"
           stop
          endif

          if (combine_idx.eq.-1) then !State_Type (combine_flag=0,1)
           cellcomp=scomp(im_crit)+1
          else if (combine_idx.ge.0) then !localMF[combine_idx] input
           cellcomp=im_crit
          else
           print *,"combine_idx invalid: ",combine_idx
           stop
          endif

          cell_temperature(im_crit)=cellfab(D_DECL(i,j,k),cellcomp)

          if (hflag.eq.0) then !inhomogeneous solid velocity
           if (cell_temperature(im_crit).ge.zero) then
            ! do nothing
           else
            print *,"cell_temperature(im_crit) must be positive: combinevel"
            print *,"cell_temperature(im_crit)=",cell_temperature(im_crit)
            print *,"im_crit=",im_crit
            print *,"cellcomp=",cellcomp
            stop
           endif
          else if (hflag.eq.1) then !homogeneous solid velocity
           if ((cell_temperature(im_crit).ge.zero).or. &
               (cell_temperature(im_crit).le.zero)) then
            ! do nothing
           else
            print *,"cell_temperature(im_crit) is NaN: ", &
               cell_temperature(im_crit)
            stop
           endif
          else
           print *,"hflag invalid3 hflag=",hflag
           stop
          endif

          if (combine_flag.eq.1) then !center->centroid

           new_temperature(im_crit)=newcell(D_DECL(i,j,k),im_crit)

           if (abs(cell_temperature(1)).le.VOFTOL) then

            if (abs(cell_temperature(1)- &
                    cell_temperature(im_crit)).le.VOFTOL) then
             !do nothing
            else
             print *,"cell_temperature invalid:", &
               im_crit, &
               cell_temperature(1), &
               cell_temperature(im_crit)
             stop
            endif

            if (abs(cell_temperature(1)- &
                    new_temperature(im_crit)).le.VOFTOL) then
             !do nothing
            else
             print *,"new_temperature invalid:", &
               im_crit, &
               cell_temperature(1), &
               new_temperature(im_crit)
             stop
            endif
           else if (abs(cell_temperature(1)).ge.VOFTOL) then
            if (abs(cell_temperature(1)- &
                    new_temperature(im_crit)).le. &
                VOFTOL*abs(cell_temperature(1))) then
             !do nothing
            else
             print *,"new_temperature invalid:", &
               im_crit, &
               cell_temperature(1), &
               new_temperature(im_crit)
             stop
            endif
           else
            print *,"new_temperature(im_crit) invalid:", &
              new_temperature(im_crit)
            stop
           endif

          else if ((combine_flag.eq.0).or. & !centroid -> center
                   (combine_flag.eq.2)) then !VOF=0 cells updated.
           ! do nothing
          else
           print *,"combine_flag invalid: ",combine_flag
           stop
          endif
 
          if (project_option.eq.SOLVETYPE_HEAT) then
           state_mass_average=state_mass_average+ &
             cell_DeDT_mfrac(im_crit)*cell_temperature(im_crit)
          else if ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
                   (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
           state_mass_average=state_mass_average+ &
             cell_species_mfrac(im_crit)*cell_temperature(im_crit)
          else
           print *,"project_option invalid: ",project_option
           stop
          endif

         enddo ! im_crit=1 .. num_materials

         if (weight_sum.gt.zero) then
          state_mass_average=state_mass_average/weight_sum
         else
          print *,"weight_sum invalid line 16060: ",weight_sum
          stop
         endif 

         if (hflag.eq.0) then !inhomogeneous solid velocity
          if (state_mass_average.ge.zero) then
           ! do nothing
          else
           print *,"state_mass_average must be nonneg: combinevel:", &
             state_mass_average
           stop
          endif
         else if (hflag.eq.1) then !homogeneous solid velocity
          if ((state_mass_average.ge.zero).or. &
              (state_mass_average.le.zero)) then
           ! do nothing
          else 
           print *,"state_mass_average is NaN: combinevel:", &
             state_mass_average
           stop
          endif
         else
          print *,"hflag invalid4 hflag=",hflag
          stop
         endif

         do im=1,num_materials

          if ((cell_vfrac(im).le.one-VOFTOL).and. &
              (cell_vfrac(im).ge.-EPS1)) then

           vofcomp=(im-1)*ngeom_recon+1
          
           if (((im_primary.eq.im).and. &
               (combine_flag.eq.0)).or. &  ! FVM -> GFM
              ((cell_vfrac(im).ge.VOFTOL).and. &
               (combine_flag.eq.1))) then  ! GFM -> FVM

            if (combine_idx.eq.-1) then
             ! do nothing
            else
             print *,"combine_idx invalid: ",combine_idx
             stop
            endif

            if (num_materials_combine.ne.num_materials) then
             print *,"num_materials_combine invalid"
             stop
            endif

            do k1=k1lo,k1hi
            do j1=-1,1
            do i1=-1,1
             call gridsten_level(xsten_ofs,i+i1,j+j1,k+k1,level,nhalf)
             call Box_volumeFAST(bfact,dx,xsten_ofs,nhalf, &
              volcell,cencell,SDIM)
             do dir=1,SDIM
              XC_sten(D_DECL(i1,j1,k1),dir)= &
               vof(D_DECL(i+i1,j+j1,k+k1),vofcomp+dir)+cencell(dir)
             enddo

             do im_local=1,num_materials
              local_vofcomp=(im_local-1)*ngeom_recon+1
              local_VOF(im_local)= &
                vof(D_DECL(i+i1,j+j1,k+k1),local_vofcomp)
             enddo !im_local=1..num_materials

             call get_primary_material_VFRAC(local_VOF, &
               im_primary_sten(D_DECL(i1,j1,k1)))

             VF_sten(D_DECL(i1,j1,k1))=local_VOF(im)
             LS_sten(D_DECL(i1,j1,k1))= &
              LS(D_DECL(i+i1,j+j1,k+k1),im)
            enddo
            enddo
            enddo ! i1,j1,k1

            if (combine_flag.eq.0) then ! centroid -> center
             do dir=1,SDIM
              xtarget(dir)=xsten(0,dir)
             enddo
            else if (combine_flag.eq.1) then ! center -> centroid
             do dir=1,SDIM
              xtarget(dir)=XC_sten(D_DECL(0,0,0),dir)
             enddo
            else
             print *,"combine_flag invalid: ",combine_flag
             stop
            endif

            do dir=1,SDIM
             xI(dir)=xtarget(dir)
            enddo

            ! check for TSAT BC.
            tsat_flag=0 !default: no TSAT available and only use cells w/F>0
            im_source_master=0
            im_dest_master=0
            TSAT_master=273.0d0

             ! check for Tgamma or Ygamma boundary condition.

            do im_crit=1,num_materials
             dencomp=STATECOMP_STATES+ &
               (im_crit-1)*num_state_material+ENUM_DENVAR+1
             tempcomp=STATECOMP_STATES+ &
               (im_crit-1)*num_state_material+ENUM_TEMPERATUREVAR+1
             Tcenter(im_crit)=cellfab(D_DECL(i,j,k),scomp(im_crit)+1)
             if (Tcenter(im_crit).ge.zero) then
              ! do nothing
             else
              print *,"Tcenter(im_crit) invalid"
              stop
             endif
             if (project_option.eq.SOLVETYPE_HEAT) then
              thermal_state(im_crit)=Tcenter(im_crit)
             else if ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
                      (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
              if (tempcomp.eq.dencomp+1) then
               thermal_state(im_crit)= &
                 state(D_DECL(i,j,k),dencomp+1)
              else
               print *,"dencomp or tempcomp invalid"
               print *,"dencomp: ",dencomp
               print *,"tempcomp: ",tempcomp
               stop
              endif
             else
              print *,"project_option invalid: ",project_option
              stop
             endif

            enddo ! im_crit=1..num_materials

            do ireverse=0,1
             do im_opp=1,num_materials
              if (im_opp.ne.im) then

               call get_iten(im,im_opp,iten)
               LL=get_user_latent_heat(iten+ireverse*num_interfaces, &
                       room_temperature,1)

               if (interface_cond_avail.eq.1) then

                if (LL.ne.zero) then

                 if (((ireverse.eq.0).and.(im.lt.im_opp)).or. &
                     ((ireverse.eq.1).and.(im.gt.im_opp))) then
                  im_source=im
                  im_dest=im_opp
                 else if (((ireverse.eq.0).and.(im.gt.im_opp)).or. &
                          ((ireverse.eq.1).and.(im.lt.im_opp))) then
                  im_source=im_opp
                  im_dest=im
                 else
                  print *,"ireverse invalid"
                  stop
                 endif

                 local_freezing_model= &
                   freezing_model(iten+ireverse*num_interfaces)
                 distribute_from_targ= &
                   distribute_from_target(iten+ireverse*num_interfaces)
                 if ((distribute_from_targ.lt.0).or. &
                     (distribute_from_targ.gt.1)) then
                  print *,"distribute_from_targ invalid"
                  stop
                 endif

                 if (is_GFM_freezing_modelF(local_freezing_model).eq.1) then 

                  if ((im_primary.eq.im).or.(im_primary.eq.im_opp)) then

                   call get_secondary_material(cell_LS,im_primary,im_secondary)

                   if (im_primary.eq.im_secondary) then
                    print *,"cannot have im_primary.eq.im_secondary"
                    stop
                   endif

                   if ((im_secondary.eq.im).or. &
                       (im_secondary.eq.im_opp)) then

                    if ((cell_vfrac(im).ge.VOFTOL).and. &
                        (cell_vfrac(im_opp).ge.VOFTOL)) then

                     Tgamma_STATUS=NINT(TgammaFAB(D_DECL(i,j,k),iten))
                     if (ireverse.eq.0) then
                      ! do nothing
                     else if (ireverse.eq.1) then
                      Tgamma_STATUS=-Tgamma_STATUS
                     else
                      print *,"ireverse invalid"
                      stop
                     endif

                     if (project_option.eq.SOLVETYPE_HEAT) then
                      ! do nothing (no need to reset Tgamma_STATUS to 0)
                     else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                              (project_option.lt. &
                               SOLVETYPE_SPEC+num_species_var)) then
                      if ((Tgamma_STATUS.eq.1).or. &
                          (Tgamma_STATUS.eq.2)) then
                       if (is_multi_component_evapF(local_freezing_model, &
                        Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces),&
                        LL).eq.0) then
                        Tgamma_STATUS=0
                       else if &
                        (is_multi_component_evapF(local_freezing_model, &
                        Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces),&
                        LL).eq.1) then
                        ispec=mass_fraction_id(iten+ireverse*num_interfaces)
                        if (ispec.eq.project_option-SOLVETYPE_SPEC+1) then
                         ! do nothing
                        else if ((ispec.ge.1).and. &
                                 (ispec.le.num_species_var)) then
                         Tgamma_STATUS=0
                        else
                         print *,"ispec invalid: ",ispec
                         stop
                        endif
                       else
                        print *,"is_multi_component_evapF invalid"
                        stop
                       endif
                      else if (Tgamma_STATUS.eq.0) then
                       ! do nothing
                      else
                       print *,"Tgamma_STATUS invalid"
                       stop
                      endif
                     else
                      print *,"project_option invalid; fort_combinevel"
                      stop
                     endif

                     if ((Tgamma_STATUS.eq.1).or. &
                         (Tgamma_STATUS.eq.2)) then

                      if (project_option.eq.SOLVETYPE_HEAT) then
                       ! default Tgamma
                       Tgamma=saturation_temp(iten+ireverse*num_interfaces)
                       TorYgamma_BC=Tgamma
                       if (Tgamma.gt.zero) then
                        tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+1
                        Tgamma=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                        TorYgamma_BC=Tgamma
                        if (Tgamma.gt.zero) then
                         ! do nothing
                        else
                         print *,"Tgamma must be positive1"
                         stop
                        endif
                       else
                        print *,"saturation temperature must be positive2"
                        stop
                       endif
                      else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                       (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
                       Tgamma=saturation_temp(iten+ireverse*num_interfaces)
                       TorYgamma_BC=one
                       if (Tgamma.gt.zero) then
                        tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+1
                        Tgamma=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                        tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+2
                        TorYgamma_BC=TgammaFAB(D_DECL(i,j,k),tsat_comp)
                        if (Tgamma.gt.zero) then
                         ! do nothing
                        else
                         print *,"Tgamma must be positive22"
                         stop
                        endif
                        if ((TorYgamma_BC.ge.zero).and. &
                            (TorYgamma_BC.le.one)) then
                         ! do nothing
                        else
                         print *,"TorYgamma_BC (aka Y) must be >= 0 and <=1"
                         stop
                        endif
                       else
                        print *,"saturation temperature must be positive33"
                        stop
                       endif
                      else
                       print *,"project_option invalid; fort_combinevel"
                       stop
                      endif

                      if (LL.lt.zero) then ! freezing
                       TDIFF=max(Tgamma-thermal_state(im), &
                                 Tgamma-thermal_state(im_opp))
                      else if (LL.gt.zero) then ! melting
                       TDIFF=max(thermal_state(im)-Tgamma, &
                                 thermal_state(im_opp)-Tgamma)
                      else
                       print *,"LL invalid; fort_combinevel"
                       stop
                      endif
                      if (tsat_flag.eq.0) then
                       tsat_flag=1
                       TSAT_master=TorYgamma_BC
                       TDIFF_master=TDIFF
                       im_source_master=im_source
                       im_dest_master=im_dest
                      else if (tsat_flag.eq.1) then
                       if (TDIFF.gt.TDIFF_master) then
                        TSAT_master=TorYgamma_BC
                        TDIFF_master=TDIFF
                        im_source_master=im_source
                        im_dest_master=im_dest
                       endif
                      else
                       print *,"tsat_flag invalid"
                       stop
                      endif 

                     else if (Tgamma_STATUS.eq.0) then
                      ! do nothing
                     else
                      print *,"Tgamma_STATUS invalid"
                      stop
                     endif

                    else if ((abs(cell_vfrac(im)).le.VOFTOL).or. &
                             (abs(cell_vfrac(im_opp)).le.VOFTOL)) then
                     ! do nothing
                    else if ((abs(cell_vfrac(im)).le.EPS1).or. &
                             (abs(cell_vfrac(im_opp)).le.EPS1)) then
                     ! do nothing
                    else
                     print *,"cell_vfrac invalid"
                     stop
                    endif
                   endif ! im_secondary==im or im_opp
                  endif ! im_primary=im or im_opp

                 else if (is_GFM_freezing_modelF( &
                           local_freezing_model).eq.0) then 
                  ! do nothing
                 else
                  print *,"local_freezing_model invalid:",local_freezing_model
                  stop
                 endif

                else if (LL.eq.zero) then
                 ! do nothing
                else
                 print *,"LL invalid"
                 stop
                endif

               else if (interface_cond_avail.eq.0) then
                ! do nothing
               else
                print *,"interface_cond_avail invalid"
                stop
               endif

              else if (im_opp.eq.im) then
               ! do nothing
              else
               print *,"im_opp invalid"
               stop
              endif
         
             enddo ! im_opp
            enddo ! ireverse

             ! combine_flag==0 (FVM->GFM) or
             ! combine_flag==1 (GFM->FVM) 
            do k1=k1lo,k1hi
            do j1=-1,1
            do i1=-1,1
             test_temperature=cellfab(D_DECL(i+i1,j+j1,k+k1),scomp(im)+1)

             if ((project_option.eq.SOLVETYPE_HEAT).or. & ! thermal combine
                 ((project_option.ge.SOLVETYPE_SPEC).and. &
                  (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then

              if (hflag.eq.0) then !inhomogeneous solid velocity

               if (test_temperature.ge.zero) then
                ! do nothing
               else
                print *,"project_option ",project_option
                print *,"test_temperature invalid:",test_temperature
                print *,"combine_flag=",combine_flag
                print *,"level,finest_level ",level,finest_level
                print *,"i,j,k ",i,j,k
                print *,"i1,j1,k1 ",i1,j1,k1
                do dir=1,SDIM
                 print *,"dir,growlo,growhi ",dir,growlo(dir),growhi(dir)
                enddo
                print *,"x= ",xsten(2*i1,1)
                print *,"y= ",xsten(2*j1,2)
                print *,"im= ",im
                print *,"INT_DIR= ",INT_DIR
                print *,"EXT_DIR= ",EXT_DIR
                print *,"REFLECT_EVEN= ",REFLECT_EVEN
                print *,"REFLECT_ODD= ",REFLECT_ODD
                print *,"FOEXTRAP= ",FOEXTRAP
                print *,"scomp_size=",scomp_size
                do im_opp=1,num_materials
                 print *,"im_opp,F ",im_opp,cell_vfrac(im_opp)
                 print *,"im_opp,scomp,ncomp ", &
                  im_opp,scomp(im_opp),ncomp(im_opp)
                 do dir=1,SDIM
                  do side=1,2
                   print *,"dir,side,im_opp,listbc ",dir,side,im_opp, &
                     listbc(dir,side,im_opp)
                  enddo
                 enddo
                enddo ! im_opp
                 
                stop
               endif

              else if (hflag.eq.1) then !homogeneous solid velocity

               if ((test_temperature.ge.zero).or. &
                   (test_temperature.le.zero)) then
                ! do nothing
               else
                print *,"test_temperature bad:",test_temperature
                stop
               endif
              else
               print *,"hflag invalid: ",hflag
               stop
              endif

             else
              print *,"project_option invalid; fort_combinevel: ", &
               project_option
              stop
             endif
        
             T_sten(D_DECL(i1,j1,k1))=test_temperature
            enddo
            enddo
            enddo ! i1,j1,k1

            if (nsolve.ne.1) then
             print *,"nsolve invalid"
             stop
            endif

            if (tsat_flag.eq.1) then

             if (is_rigid(im_primary).eq.0) then

              im_source=im_source_master
              im_dest=im_dest_master
              LS_source=LS(D_DECL(i,j,k),im_source)
              LS_dest=LS(D_DECL(i,j,k),im_dest)

              if ((abs(LS_source).le.two*dxmaxLS).and. &
                  (abs(LS_dest).le.two*dxmaxLS)) then
               if (LS_dest.ge.zero) then
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),num_materials+(im_source-1)*SDIM+dir)
                 xI(dir)=xsten(0,dir)-LS_source*nrm(dir)
                enddo
               else if (LS_source.ge.zero) then
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),num_materials+(im_dest-1)*SDIM+dir)
                 xI(dir)=xsten(0,dir)-LS_dest*nrm(dir)
                enddo
               else if ((LS_dest.le.zero).and. &
                        (LS_source.le.zero)) then
                tsat_flag=0
               else
                print *,"LS_dest or LS_source invalid"
                stop
               endif
              else if ((abs(LS_source).ge.two*dxmaxLS).or. &
                       (abs(LS_dest).ge.two*dxmaxLS)) then
               tsat_flag=0
              else
               print *,"LS_dest or LS_source invalid"
               stop
              endif

             else if (is_rigid(im_primary).eq.1) then
              tsat_flag=0
             else
              print *,"im_primary invalid"
              stop
             endif

            else if (tsat_flag.eq.0) then
             ! do nothing
            else
             print *,"tsat_flag invalid: ",tsat_flag
             stop
            endif

             ! this subroutine: fort_combinevel
             ! subroutine center_centroid_interchange is declared in:
             ! MASS_TRANSFER_3D.F90
            call center_centroid_interchange( &
             DATA_FLOOR, &
             nsolve, &
             combine_flag,  & ! 0=>centroid -> center   1=>center->centroid
             tsat_flag, & !-1=>use all cells in stencil;=0=>no tsat;=1=>tsat
             bfact, &
             level, &
             finest_level, &
             dx,xlo, &
             xsten,nhalf, &
             T_sten, &  
             XC_sten, &  
             xI, &
             xtarget, &
             im, &
             im_primary_sten, &
             VF_sten, &
             LS_sten, &
             TSAT_master, &
             T_out)

            if (combine_flag.eq.0) then !centroid -> center (FVM->GFM)

             if (im_primary.eq.im) then
              !do nothing
             else
              print *,"expecting im=im_primary:",im,im_primary
              stop
             endif

             if (VF_sten(D_DECL(0,0,0)).gt.VOFTOL) then
              ! do nothing
             else
              print *,"expecting VF_sten(0,0,0)>VOFTOL"
              stop
             endif
             if (cell_vfrac(im).gt.VOFTOL) then
              ! do nothing
             else
              print *,"expecting cell_vfrac(im)>VOFTOL"
              stop
             endif

             if (tsat_flag.eq.0) then

              !  T_out(1)=T_sten(D_DECL(0,0,0))
              T_out(1)=state_mass_average

             else if (tsat_flag.eq.1) then
              ! do nothing
             else
              print *,"tsat_flag invalid"
              stop
             endif

             do im_opp=1,num_materials
              newcell(D_DECL(i,j,k),im_opp)=T_out(1)
             enddo

            else if (combine_flag.eq.1) then  !center -> centroid (GFM->FVM)

             if (tsat_flag.eq.0) then

              ! T_out(1)=cellfab(D_DECL(i,j,k),scomp(im_primary)+1)
              T_out(1)=state_mass_average

             else if (tsat_flag.eq.1) then
              ! do nothing
             else
              print *,"tsat_flag invalid"
              stop
             endif

             newcell(D_DECL(i,j,k),im)=T_out(1)

            else
             print *,"combine_flag invalid: ",combine_flag
             stop
            endif

           else if (combine_flag.eq.0) then !centroid->center (FVM->GFM)

            if (im_primary.ne.im) then
             ! do nothing
            else
             print *,"combine_flag,im, or im_primary invalid: ", &
                combine_flag,im,im_primary
             stop
            endif

           else if ((combine_flag.eq.1).or. & ! GFM->FVM
                    (combine_flag.eq.2)) then ! combine if F==0

            if (abs(cell_vfrac(im)).le.VOFTOL) then

             if (nsolve.ne.1) then
              print *,"nsolve invalid: ",nsolve
              stop
             endif
             if (num_materials_combine.ne.num_materials) then
              print *,"num_materials_combine invalid: ", &
                      num_materials_combine
              stop
             endif

             if (combine_idx.eq.-1) then
              cellcomp=scomp(im)+1
             else if (combine_idx.ge.0) then
              cellcomp=im
             else
              print *,"combine_idx invalid"
              stop
             endif

             if (combine_flag.eq.1) then ! GFM -> FVM (after diffusion)
              newcell(D_DECL(i,j,k),im)=state_mass_average
             else if (combine_flag.eq.2) then !combine if F==0
              cellfab(D_DECL(i,j,k),cellcomp)=state_mass_average
             else
              print *,"combine_flag invalid:",combine_flag
              stop
             endif
 
            else if ((cell_vfrac(im).ge.VOFTOL).and. &
                     (cell_vfrac(im).le.one+EPS1).and. &
                     (combine_flag.eq.2)) then !combine if vfrac<VOFTOL
             ! do nothing
            else
             print *,"cell_vfrac or combine_flag bust"
             print *,"im=",im
             print *,"cell_vfrac=",cell_vfrac(im)
             print *,"combine_flag=",combine_flag
             stop
            endif

           else
            print *,"combine_flag invalid: ",combine_flag
            stop
           endif
   
          else if (cell_vfrac(im).ge.one-VOFTOL) then
 
           if (im_primary.eq.im) then
            !do nothing
           else
            print *,"im,im_primary: ",im,im_primary
            stop
           endif

           if (combine_flag.eq.0) then ! centroid -> center

            T_out(1)=cellfab(D_DECL(i,j,k),scomp(im)+1)

            if (abs(state_mass_average).le.VOFTOL) then
             if (abs(state_mass_average-T_out(1)).le.VOFTOL) then
              !do nothing
             else
              print *,"T_out(1),state_mass_average invalid:", &
                T_out(1),state_mass_average
              stop
             endif
            else if (abs(state_mass_average).ge.VOFTOL) then
             if (abs(state_mass_average-T_out(1)).le. &
                 VOFTOL*abs(state_mass_average)) then
              !do nothing
             else
              print *,"T_out(1),state_mass_average invalid:", &
                T_out(1),state_mass_average
              stop
             endif
            else
             print *,"state_mass_average invalid:",state_mass_average
             stop
            endif
                             
             ! since cell_vfrac(im)==1 => cell_vfrac(im_opp)==0.0 (im_opp!=im)
            do im_opp=1,num_materials
             newcell(D_DECL(i,j,k),im_opp)=T_out(1)
            enddo

           else if (combine_flag.eq.1) then ! center->centroid
            ! do nothing
           else if (combine_flag.eq.2) then ! extrap where vfrac=0
            ! do nothing
           else
            print *,"combine_flag invalid: ",combine_flag
            stop
           endif
          
          else
           print *,"cell_vfrac or combine_flag bust"
           print *,"im,cell_vfrac=",im,cell_vfrac(im)
           print *,"combine_flag " ,combine_flag
           stop
          endif

         enddo ! im=1,num_materials

        else
         print *,"project_option invalid: ",project_option
         stop
        endif

       else if (mask_test.eq.0) then
        ! do nothing
       else
        print *,"mask_test invalid: ",mask_test
        stop
       endif
 
      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine fort_combinevel

       ! combine_flag==2 (only overwrite if F=0)
      subroutine fort_combinevelface( &
       tid, &
       hflag, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       combine_idx, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       velbc, &
       vof,DIMS(vof), &
       mac,DIMS(mac), &
       xface,DIMS(xface), &
       LS,DIMS(LS), &
       solfab,DIMS(solfab), &
       xlo,dx, &
       dir, &
       cur_time) &
      bind(c,name='fort_combinevelface')
      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
 
      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: hflag

      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def)
      integer, INTENT(in) :: combine_idx

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: bfact
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: velbc(SDIM,2,SDIM)

      integer, INTENT(in) :: DIMDEC(vof)
      integer, INTENT(in) :: DIMDEC(mac)
      integer, INTENT(in) :: DIMDEC(xface)
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(solfab)
  
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM) 
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: cur_time

      real(amrex_real), INTENT(in),target :: vof(DIMV(vof),num_materials*ngeom_recon)
      real(amrex_real), pointer :: vof_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: mac(DIMV(mac))
      real(amrex_real), pointer :: mac_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: xface(DIMV(xface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: xface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: LS(DIMV(LS),num_materials*(SDIM+1))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: solfab(DIMV(solfab),nparts_def*SDIM)
      real(amrex_real), pointer :: solfab_ptr(D_DECL(:,:,:),:)

      integer ii,jj,kk 
      integer i,j,k
      integer icell,jcell,kcell
      integer side
      integer im
      integer idx
      integer, PARAMETER :: nhalf=3

      real(amrex_real) ucombine
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer is_solid_face
      real(amrex_real) vsol
      integer iboundary
      integer at_RZ_face

      real(amrex_real) face_vfrac_cell(num_materials)
      real(amrex_real) face_vfrac(2,num_materials)
      real(amrex_real) face_mfrac(2,num_materials)
      real(amrex_real) fluid_vfrac_face
      real(amrex_real) total_vol_face
      real(amrex_real) total_vol_face_cell
      real(amrex_real) local_volume
      real(amrex_real) local_mass
      integer vofcomp
      integer velcomp
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      integer nmax
      integer im_solid_crit,imL,imR
      integer partid
      integer partid_crit
      integer partidL,partidR
      real(amrex_real) lsleft(num_materials)
      real(amrex_real) lsright(num_materials)
      real(amrex_real) xface_local
      integer dir_local
      real(amrex_real) xclamped_face(SDIM)
      real(amrex_real) xclamped_minus(SDIM)
      real(amrex_real) xclamped_plus(SDIM)
      real(amrex_real) LS_clamped_minus
      real(amrex_real) LS_clamped_plus
      real(amrex_real) vel_clamped_minus(SDIM)
      real(amrex_real) vel_clamped_plus(SDIM)
      real(amrex_real) vel_clamped_face(SDIM)
      real(amrex_real) temperature_clamped_minus
      real(amrex_real) temperature_clamped_plus
      integer prescribed_flag
      integer is_clamped_face
      integer local_tessellate

      nmax=POLYGON_LIST_MAX ! in: COMBINEVELFACE

      if ((tid.lt.0).or. &
          (tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid"
       stop
      endif

      if ((hflag.ne.0).and.(hflag.ne.1)) then
       print *,"hflag invalid5 hflag=",hflag
       stop
      endif 
      if (combine_idx.lt.-1) then
       print *,"combine_idx invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid76"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid combinevel face"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_combinevelface"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_combinevelface"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir out of range in fort_combinevelface"
       stop
      endif

      vof_ptr=>vof
      call checkbound_array(fablo,fabhi,vof_ptr,1,-1)
      mac_ptr=>mac
      call checkbound_array1(fablo,fabhi,mac_ptr,0,dir)
      xface_ptr=>xface
      call checkbound_array(fablo,fabhi,xface_ptr,0,dir)
      solfab_ptr=>solfab
      call checkbound_array(fablo,fabhi,solfab_ptr,0,dir)
      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,0,dir)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir)
       do dir_local=1,SDIM
        xclamped_face(dir_local)=xstenMAC(0,dir_local)
       enddo

       do side=1,2
        if (side.eq.1) then
         icell=i-ii
         jcell=j-jj
         kcell=k-kk
        else if (side.eq.2) then
         icell=i
         jcell=j
         kcell=k
        else
         print *,"side invalid"
         stop
        endif

        call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)
        do dir_local=1,SDIM
         if (side.eq.1) then
          xclamped_minus(dir_local)=xsten(0,dir_local)
         else if (side.eq.2) then
          xclamped_plus(dir_local)=xsten(0,dir_local)
         else
          print *,"side invalid"
          stop
         endif
        enddo ! dir_local=1..sdim
       enddo ! side=1,2

       if (dir.eq.0) then
        idx=i
       else if (dir.eq.1) then
        idx=j
       else if ((dir.eq.2).and.(SDIM.eq.3)) then
        idx=k
       else
        print *,"dir invalid combinevel"
        stop
       endif

       do im=1,num_materials
        face_vfrac_cell(im)=zero
       enddo
       total_vol_face_cell=zero

       fluid_vfrac_face=zero
       total_vol_face=zero

       at_RZ_face=0
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if ((levelrz.eq.COORDSYS_RZ).or. &
                (levelrz.eq.COORDSYS_CYLINDRICAL)) then
        if ((dir.eq.0).and. &
            (xstenMAC(0,1).le.EPS2*dx(1))) then
         at_RZ_face=1
        endif
       else
        print *,"levelrz invalid combine vel"
        stop
       endif

       is_solid_face=0

       vsol=zero

       if (at_RZ_face.eq.1) then

        is_solid_face=1
        vsol=zero

       else if (at_RZ_face.eq.0) then

        do side=1,2
 
         if (side.eq.1) then
          icell=i-ii
          jcell=j-jj
          kcell=k-kk
         else if (side.eq.2) then
          icell=i
          jcell=j
          kcell=k
         else
          print *,"side invalid"
          stop
         endif

         call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)

         do im=1,num_materials*ngeom_recon
          mofdata(im)=vof(D_DECL(icell,jcell,kcell),im)
         enddo
         local_tessellate=3
          !EPS2 tolerance
         call multi_get_volume_tessellate( &
          local_tessellate, & !  =3
          bfact, &
          dx,xsten, &
          nhalf, & !=3
          mofdata, &
          geom_xtetlist(1,1,1,tid+1), &
          nmax, &
          nmax, &
          SDIM)

         do im=1,num_materials

          vofcomp=(im-1)*ngeom_recon+1
          local_volume=mofdata(vofcomp)

          if ((local_volume.ge.-EPS1).and. &
              (local_volume.le.one+EPS1)) then
           if (local_volume.lt.zero) then
            local_volume=zero
           endif
           if (local_volume.gt.one) then
            local_volume=one
           endif
          else
           print *,"local_volume invalid"
           stop
          endif

          face_vfrac_cell(im)=face_vfrac_cell(im)+local_volume
          total_vol_face_cell=total_vol_face_cell+local_volume

           !side=1,2
          local_volume=xface(D_DECL(i,j,k),FACECOMP_VOFFACE+2*(im-1)+side)
          local_mass=xface(D_DECL(i,j,k),FACECOMP_MASSFACE+2*(im-1)+side)

          total_vol_face=total_vol_face+local_volume
          if (is_rigid(im).eq.0) then
           fluid_vfrac_face=fluid_vfrac_face+local_volume
          else if (is_rigid(im).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif

          face_vfrac(side,im)=local_volume
          face_mfrac(side,im)=local_mass

         enddo ! im=1,num_materials

        enddo ! side=1,2

        if (total_vol_face.gt.zero) then
         ! do nothing
        else
         print *,"fluid_vfrac_face ",fluid_vfrac_face
         print *,"total_vol_face ",total_vol_face
         print *,"total_vol_face bust combinevel1" 
         stop
        endif
        if (total_vol_face_cell.gt.zero) then
         ! do nothing
        else
         print *,"total_vol_face_cell invalid"
         stop
        endif

        fluid_vfrac_face=fluid_vfrac_face/total_vol_face
        if ((fluid_vfrac_face.ge.zero).and. &
            (fluid_vfrac_face.le.one+EPS1)) then
         ! do nothing
        else
         print *,"fluid_vfrac_face invalid"
         stop
        endif

        im_solid_crit=0
        partid_crit=0
        partid=0

        do im=1,num_materials
         face_vfrac_cell(im)=face_vfrac_cell(im)/total_vol_face_cell
         lsleft(im)=LS(D_DECL(i-ii,j-jj,k-kk),im)
         lsright(im)=LS(D_DECL(i,j,k),im)
         if (is_rigid(im).eq.1) then
          if (im_solid_crit.eq.0) then
           im_solid_crit=im
           partid_crit=partid
          else if ((im_solid_crit.ge.1).and. &
                   (im_solid_crit.le.num_materials)) then
           if (lsleft(im)+lsright(im).ge. &
               lsleft(im_solid_crit)+lsright(im_solid_crit)) then
            im_solid_crit=im
            partid_crit=partid
           endif
          else
           print *,"im_solid_crit invalid combinevelface:",im_solid_crit
           stop
          endif
         else if (is_rigid(im).eq.0) then
          ! do nothing
         else
          print *,"is_rigid(im) invalid"
          stop
         endif

         if (is_lag_part(im).eq.1) then
          partid=partid+1
         else if (is_lag_part(im).eq.0) then
          ! do nothing
         else
          print *,"is_lag_part invalid"
          stop
         endif

        enddo ! im=1..num_materials

        if (partid.ne.nparts) then
         print *,"partid.ne.nparts"
         stop
        endif

        call get_primary_material(lsleft,imL)
        call get_primary_material(lsright,imR)

        partidL=0
        partidR=0

        if (is_rigid(imL).eq.1) then 
         do im=1,imL-1
          if (is_lag_part(im).eq.1) then
           partidL=partidL+1
          else if (is_lag_part(im).eq.0) then
           ! do nothing
          else
           print *,"is_lag_part(im) invalid"
           stop
          endif
         enddo ! im=1..imL-1
         if (im_solid_map(partidL+1)+1.ne.imL) then
          print *,"im_solid_map(partidL+1)+1.ne.imL"
          stop
         endif
        else if (is_rigid(imL).eq.0) then
         ! do nothing 
        else
         print *,"is_rigid(imL) invalid"
         stop
        endif

        if (is_rigid(imR).eq.1) then
         do im=1,imR-1
          if (is_lag_part(im).eq.1) then
           partidR=partidR+1
          else if (is_lag_part(im).eq.0) then
           ! do nothing
          else
           print *,"is_lag_part(im) invalid"
           stop
          endif
         enddo ! im=1..imR-1
         if (im_solid_map(partidR+1)+1.ne.imR) then
          print *,"im_solid_map(partidR+1)+1.ne.imR"
          stop
         endif
        else if (is_rigid(imR).eq.0) then
         ! do nothing 
        else
         print *,"is_rigid(imR) invalid"
         stop
        endif

        if (fluid_vfrac_face.le.VOFTOL_AREAFRAC) then
         is_solid_face=2
        else if ((fluid_vfrac_face.ge.VOFTOL_AREAFRAC).and. &
                 (fluid_vfrac_face.le.one+EPS1)) then
         xface_local=xface(D_DECL(i,j,k),FACECOMP_FACECUT+1)
         if ((xface_local.ge.zero).and.(xface_local.le.half)) then
          xface_local=zero
          is_solid_face=3
         else if ((xface_local.ge.half).and.(xface_local.le.one)) then
          if ((is_rigid(imR).eq.1).or. &
              (is_rigid(imL).eq.1)) then
           is_solid_face=3
          else if ((is_rigid(imR).eq.0).and. &
                   (is_rigid(imL).eq.0)) then
           ! do nothing
          else
           print *,"imR or imL invalid"
           stop
          endif
         else
          print *,"xface_local invalid"
          stop
         endif
        else
         print *,"fluid_vfrac_face invalid"
         stop
        endif

        if (hflag.eq.0) then
         if ((im_solid_crit.ge.1).and. &
             (im_solid_crit.le.num_materials)) then
          if (im_solid_map(partid_crit+1)+1.ne.im_solid_crit) then
           print *,"im_solid_map(partid_crit+1)+1.ne.im_solid_crit"
           stop
          endif
          velcomp=partid_crit*SDIM+dir+1
          vsol=solfab(D_DECL(i,j,k),velcomp)
         else if (im_solid_crit.eq.0) then
          vsol=zero
         else
          print *,"im_solid_crit invalid combinevelface2:",im_solid_crit
          stop
         endif
        else if (hflag.eq.1) then
         vsol=zero
        else
         print *,"hflag invalid6 hflag=",hflag
         stop
        endif

       else
        print *,"at_RZ_face invalid"
        stop
       endif

       iboundary=0
       if (idx.eq.fablo(dir+1)) then
        iboundary=1
        side=1
       else if (idx.eq.fabhi(dir+1)+1) then
        iboundary=1
        side=2
       else if ((idx.gt.fablo(dir+1)).and. &
                (idx.lt.fabhi(dir+1)+1)) then
        ! do nothing
       else
        print *,"idx invalid"
        stop
       endif

       if (iboundary.eq.1) then
        if (velbc(dir+1,side,dir+1).eq.REFLECT_ODD) then
         vsol=zero
         is_solid_face=4
        else if (velbc(dir+1,side,dir+1).eq.EXT_DIR) then
         if (hflag.eq.0) then
          call velbc_override(cur_time,dir+1,side,dir+1, &
           vsol, &
           xstenMAC,nhalf,dx,bfact)
         else if (hflag.eq.1) then
          vsol=zero
         else
          print *,"hflag invalid7 hflag=",hflag
          stop
         endif
         is_solid_face=5
        else if ((velbc(dir+1,side,dir+1).eq.INT_DIR).or. &
                 (velbc(dir+1,side,dir+1).eq.REFLECT_EVEN).or. &
                 (velbc(dir+1,side,dir+1).eq.FOEXTRAP)) then
         ! do nothing
        else
         print *,"velbc invalid"
         stop
        endif
       else if (iboundary.eq.0) then
        ! do nothing
       else
        print *,"iboundary invalid"
        stop
       endif

         ! LS>0 if clamped
       call SUB_clamped_LS(xclamped_minus,cur_time,LS_clamped_minus, &
        vel_clamped_minus,temperature_clamped_minus,prescribed_flag,dx)
       call SUB_clamped_LS(xclamped_plus,cur_time,LS_clamped_plus, &
        vel_clamped_plus,temperature_clamped_plus,prescribed_flag,dx)
       if ((LS_clamped_minus.ge.zero).or. &
           (LS_clamped_plus.ge.zero)) then
        is_clamped_face=1
        do dir_local=1,SDIM
         if (LS_clamped_minus.lt.zero) then
          vel_clamped_face(dir_local)=vel_clamped_plus(dir_local)
          is_clamped_face=2
         else if (LS_clamped_plus.lt.zero) then
          vel_clamped_face(dir_local)=vel_clamped_minus(dir_local)
          is_clamped_face=3
         else
          vel_clamped_face(dir_local)=half*(vel_clamped_plus(dir_local)+ &
            vel_clamped_minus(dir_local))
         endif
        enddo ! dir_local=1..sdim
       else if ((LS_clamped_minus.lt.zero).and. &
                (LS_clamped_plus.lt.zero)) then
        is_clamped_face=0
       else
        print *,"LS_clamped plus or minus is NaN"
        stop
       endif

       if ((is_clamped_face.eq.1).or. &
           (is_clamped_face.eq.2).or. &
           (is_clamped_face.eq.3)) then
        if (is_solid_face.eq.0) then
         is_solid_face=2
         if (hflag.eq.1) then
          vsol=zero
         else if (hflag.eq.0) then
          vsol=vel_clamped_face(dir+1)
         else
          print *,"hflag invalid"
          stop
         endif
        else if ((is_solid_face.ge.1).and. &
                 (is_solid_face.le.5)) then
         ! do nothing
        else
         print *,"is_solid_face invalid"
         stop
        endif
       else if (is_clamped_face.eq.0) then
        ! do nothing
       else
        print *,"is_clamped_face invalid"
        stop
       endif

       if (is_solid_face.eq.1) then ! RZ face
        ucombine=vsol
       else if (is_solid_face.eq.4) then ! REFLECT_ODD
        ucombine=vsol
       else if (is_solid_face.eq.5) then ! EXT_DIR
        ucombine=vsol
       else if ((is_solid_face.eq.2).or. &
                (is_solid_face.eq.3)) then
        ucombine=vsol
       else if (is_solid_face.eq.0) then
        ucombine=mac(D_DECL(i,j,k))
       else
        print *,"is_solid_face invalid 1 ",is_solid_face
        stop
       endif
       mac(D_DECL(i,j,k))=ucombine

      enddo
      enddo
      enddo

      return
      end subroutine fort_combinevelface

       ! PART I (low order):
       !
       ! itensor_iter==0 => face grad u
       !    tileloop==0  => low order
       !       spectral_loop==0 => low order
       !       spectral_loop==1 => do nothing
       !    tileloop==1  => do nothing
       ! itensor_iter==1 => cell grad u
       !    tileloop==0  => low order
       !
       ! PART II (high order)
       !
       ! tileloop==0 => do nothing
       !
       ! tileloop==1
       !   itensor_iter==0 => face grad u
       !     spectral_loop==0 initialize semflux
       !     spectral_loop==1 flux=average of values shared at face.
       !   itensor_iter==1 => cell grad u  
       !   (only called when spectral_loop==0)
       !     
      subroutine fort_face_gradients( &
       ns_time_order, &
       divu_outer_sweeps, &
       num_divu_outer_sweeps, &
       SDC_outer_sweeps, &
       tileloop, &
       dir, &
       slab_step, &
       itensor_iter, & ! 0=>gradient face  1=>gradient at cell
       time, &
       enable_spectral, &
       velbc, &
       spectral_loop, &
       ncfluxreg, &
       semflux,DIMS(semflux), &
       amrsync,DIMS(amrsync), &
       mask0,DIMS(mask0), &  ! mask0=1 if not cov. by finer level or outside. 
       mask3,DIMS(mask3), &  ! 1 at fine-fine ghost cells 0 at other ghost.
       maskSEM,DIMS(maskSEM), &
       faceLS,DIMS(faceLS), & 
       mdata,DIMS(mdata), & 
       tdata,DIMS(tdata), & 
       c_tdata,DIMS(c_tdata), & 
       vel,DIMS(vel), &
       solidx,DIMS(solidx), &
       solidy,DIMS(solidy), &
       solidz,DIMS(solidz), &
       levelpc,DIMS(levelpc), &
       xlo,dx, &
       rzflag, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact,bfact_c,bfact_f, &
       level, &
       finest_level, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       homflag, &
       simple_AMR_BC_flag_viscosity) &
      bind(c,name='fort_face_gradients')

      use probcommon_module
      use global_utility_module
      use probf90_module
      use MOF_routines_module
 
      IMPLICIT NONE

      real(amrex_real) def_dt

      integer, INTENT(in) :: ns_time_order
      integer, INTENT(in) :: divu_outer_sweeps
      integer, INTENT(in) :: num_divu_outer_sweeps
      integer, INTENT(in) :: SDC_outer_sweeps
      integer, INTENT(in) :: tileloop
      integer, INTENT(in) :: dir
      integer, INTENT(in) :: itensor_iter
      integer, INTENT(in) :: spectral_loop
      integer, INTENT(in) :: slab_step
      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level

      integer, INTENT(in) :: simple_AMR_BC_flag_viscosity

      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def)
      integer, INTENT(in) :: ncfluxreg
      integer, INTENT(in) :: rzflag 
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact,bfact_c,bfact_f
      integer, INTENT(in) :: DIMDEC(semflux)
      integer, INTENT(in) :: DIMDEC(amrsync)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: DIMDEC(mask0)
      !1 at fine-fine ghost cells 0 at other ghost.
      integer, INTENT(in) :: DIMDEC(mask3) 
      integer, INTENT(in) :: DIMDEC(faceLS)
      integer, INTENT(in) :: DIMDEC(mdata)
      integer, INTENT(in) :: DIMDEC(tdata)
      integer, INTENT(in) :: DIMDEC(c_tdata)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(solidx)
      integer, INTENT(in) :: DIMDEC(solidy)
      integer, INTENT(in) :: DIMDEC(solidz)
      integer, INTENT(in) :: DIMDEC(levelpc)
 
      integer, INTENT(in) :: velbc(SDIM,2,SDIM)
      real(amrex_real), INTENT(in) :: time
 
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM) 
      real(amrex_real), INTENT(inout), target :: semflux(DIMV(semflux),ncfluxreg)
      real(amrex_real), pointer :: semflux_ptr(D_DECL(:,:,:),:)

       ! INTENT(inout) instead of INTENT(in) since
       ! this parameter doubles as "xp" in SEM_CELL_TO_MAC
      real(amrex_real), INTENT(inout), target :: amrsync(DIMV(amrsync),SDIM)
      real(amrex_real), pointer :: amrsync_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
       ! mask0=tag if not covered by level+1 or outside the domain.
      real(amrex_real), INTENT(in), target :: mask0(DIMV(mask0))
      real(amrex_real), pointer :: mask0_ptr(D_DECL(:,:,:))
       ! mask3=tag at exterior fine/fine border.
       ! mask3=1-tag at other exterior boundaries.
      real(amrex_real), INTENT(in), target :: mask3(DIMV(mask3))
      real(amrex_real), pointer :: mask3_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: faceLS(DIMV(faceLS),SDIM)
      real(amrex_real), pointer :: faceLS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out), target :: mdata(DIMV(mdata),SDIM)
      real(amrex_real), pointer :: mdata_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              tdata(DIMV(tdata),AMREX_SPACEDIM_SQR)
      real(amrex_real), pointer :: tdata_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out), target :: &
              c_tdata(DIMV(c_tdata),AMREX_SPACEDIM_SQR)
      real(amrex_real), pointer :: c_tdata_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: vel(DIMV(vel),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              solidx(DIMV(solidx),nparts_def*SDIM)
      real(amrex_real), pointer :: solidx_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              solidy(DIMV(solidy),nparts_def*SDIM)
      real(amrex_real), pointer :: solidy_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              solidz(DIMV(solidz),nparts_def*SDIM)
      real(amrex_real), pointer :: solidz_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              levelpc(DIMV(levelpc),num_materials)
      real(amrex_real), pointer :: levelpc_ptr(D_DECL(:,:,:),:)
  
      integer i,j,k
      integer dir2
      integer dirtan
      integer mask_boundary
      integer local_maskSEM
      integer maskcov
      integer index_adjoin(2,3)
      integer masktest,bctest,icrit
      integer ii,jj,kk
      integer ishift,jshift,kshift
      integer shift_flag
      integer im1,jm1,km1
      integer ip1,jp1,kp1
      integer nc
      real(amrex_real) delta
      integer side
      integer nbase
      real(amrex_real) vplus(SDIM)
      real(amrex_real) vminus(SDIM)

      real(amrex_real) lsleft(num_materials)
      real(amrex_real) lsright(num_materials)

      real(amrex_real) solidvelleft(SDIM)
      real(amrex_real) solidvelright(SDIM)
      integer indexleft,indexright
      integer im,imL,imR
      integer homflag

      real(amrex_real) hold_grad
      real(amrex_real) RR
      integer, PARAMETER :: nhalf=1
      integer, PARAMETER :: nhalfcell=3
      integer, PARAMETER :: nhalfclamped=1

      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten(-nhalfcell:nhalfcell,SDIM)
      integer scomp,scomp_bc,dcomp
      integer ncomp_source
      integer ncomp_dest
      integer ncomp_xvel
      integer ncomp_denold
      integer ncomp_veldest
      integer ncomp_dendest
      integer ncomp_cterm
      integer operation_flag
      integer energyflag 
      integer project_option_vel
      real(amrex_real) lspoint(num_materials)
      integer im_primary
      real(amrex_real) gradu(SDIM)
      real(amrex_real) int_xlo,int_xhi
      real(amrex_real) leftwt,rightwt
      real(amrex_real) slopeLT,slopeRT
      real(amrex_real) theta_factor
      integer stripstat
      integer ielem,jelem,kelem
      integer velcomp
      integer tensorcomponent
      integer ncomp_xgp
      integer ncomp_xp
      integer left_rigid,right_rigid
      integer partidL,partidR

      real(amrex_real) xclamped_minus_sten(-nhalfclamped:nhalfclamped,SDIM)
      real(amrex_real) xclamped_plus_sten(-nhalfclamped:nhalfclamped,SDIM)
      real(amrex_real) xclamped_minus(SDIM)
      real(amrex_real) xclamped_plus(SDIM)
      real(amrex_real) xclamped_cen(SDIM)
      real(amrex_real) LS_clamped_cen
      real(amrex_real) LS_clamped_minus
      real(amrex_real) LS_clamped_plus
      real(amrex_real) vel_clamped_cen(SDIM)
      real(amrex_real) vel_clamped_minus(SDIM)
      real(amrex_real) vel_clamped_plus(SDIM)
      real(amrex_real) vel_clamped_face(SDIM)
      real(amrex_real) temperature_clamped_cen
      real(amrex_real) temperature_clamped_minus
      real(amrex_real) temperature_clamped_plus
      integer prescribed_flag
      integer is_clamped_face

      semflux_ptr=>semflux
      amrsync_ptr=>amrsync
      faceLS_ptr=>faceLS
      mdata_ptr=>mdata
      tdata_ptr=>tdata
      c_tdata_ptr=>c_tdata
      mask3_ptr=>mask3
      mask0_ptr=>mask0
      maskSEM_ptr=>maskSEM
      vel_ptr=>vel

      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid face gradients"
       stop
      endif

      if ((slab_step.lt.-1).or.(slab_step.gt.bfact_time_order)) then
       print *,"slab_step invalid face gradients "
       stop
      endif
      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid: ",time
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid77"
       stop
      endif
      if ((bfact_c.ne.bfact).and.(bfact_c.ne.2*bfact)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_f.ne.bfact).and.(bfact.ne.2*bfact_f)) then
       print *,"bfact_f invalid"
       stop
      endif
      if ((spectral_loop.ne.0).and. &
          (spectral_loop.ne.1)) then
       print *,"spectral_loop invalid"
       stop
      endif
      if (itensor_iter.eq.0) then ! face grad U
       if ((ncfluxreg/SDIM)*SDIM.ne.ncfluxreg) then
        print *,"ncfluxreg invalid15 ",ncfluxreg
        stop
       endif
       if (ncfluxreg.ne.AMREX_SPACEDIM_SQR) then
        print *,"ncfluxreg invalid16 ",ncfluxreg
        stop
       endif
      else if (itensor_iter.eq.1) then ! cell grad U
       if (ncfluxreg.lt.SDIM) then
        print *,"ncfluxreg invalid17 ",ncfluxreg
        stop
       endif
       if (spectral_loop.eq.0) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        print *,"unnecessary to average cell based gradients"
        stop
       endif
      else
       print *,"itensor_iter invalid"
       stop
      endif

      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.1)) then
       print *,"enable_spectral invalid face gradients"
       stop
      endif
 
      if ((simple_AMR_BC_flag_viscosity.eq.0).or. &
          (simple_AMR_BC_flag_viscosity.eq.1)) then
       ! do nothing
      else
       print *,"simple_AMR_BC_flag_viscosity invalid"
       stop
      endif

      do im=1,num_materials

       if (fort_denconst(im).gt.zero) then
        !do nothing
       else
        print *,"fort_denconst invalid: ",im,fort_denconst(im)
        stop
       endif

      enddo ! im=1..num_materials

       ! in: fort_face_gradients
      if ((ns_time_order.ge.1).and.(ns_time_order.le.32)) then
       ! do nothing
      else
       print *,"ns_time_order invalid"
       stop
      endif
      if ((divu_outer_sweeps.ge.0).and. &
          (divu_outer_sweeps.lt.num_divu_outer_sweeps)) then
       ! do nothing
      else
       print *,"divu_outer_sweeps invalid fort_face_gradients"
       stop
      endif
       ! in: fort_face_gradients
      if ((SDC_outer_sweeps.ge.0).and. &
          (SDC_outer_sweeps.lt.ns_time_order)) then
       ! do nothing
      else
       print *,"SDC_outer_sweeps invalid in fort_face_gradients"
       stop
      endif

      if (rzflag.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (rzflag.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"rzflag invalid"
       stop
      endif

      call checkbound_array(fablo,fabhi,amrsync_ptr,0,dir-1)

      call checkbound_array(fablo,fabhi,semflux_ptr,1,-1)
      call checkbound_array(fablo,fabhi,faceLS_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,mask0_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,mask3_ptr,1,-1)
      call checkbound_array(fablo,fabhi,mdata_ptr,1,-1)
      call checkbound_array(fablo,fabhi,tdata_ptr,1,-1)
      call checkbound_array(fablo,fabhi,c_tdata_ptr,1,-1)
      call checkbound_array(fablo,fabhi,vel_ptr,1,-1)

      solidx_ptr=>solidx
      solidy_ptr=>solidy
      solidz_ptr=>solidz
      call checkbound_array(fablo,fabhi,solidx_ptr,0,0)
      call checkbound_array(fablo,fabhi,solidy_ptr,0,1)
      call checkbound_array(fablo,fabhi,solidz_ptr,0,SDIM-1)

      levelpc_ptr=>levelpc
      call checkbound_array(fablo,fabhi,levelpc_ptr,2,-1)

      call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)

      ii=0
      jj=0
      kk=0
      if (dir.eq.1) then
       ii=1
       nbase=TENSOR_TRANSPOSE_UX-1
      else if (dir.eq.2) then
       jj=1
       nbase=TENSOR_TRANSPOSE_UY-1
      else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
       kk=1
       nbase=TENSOR_TRANSPOSE_UZ-1
      else
       print *,"dir invalid face gradients 2"
       stop
      endif

      if (itensor_iter.eq.0) then  ! face grad U

       if (tileloop.eq.0) then ! low order

        if (spectral_loop.eq.0) then

          ! same as growntileboxMAC, except includes one layer of ghost
          ! cells in the tangential directions.
         call growntileboxTENSOR(tilelo,tilehi,fablo,fabhi, &
           growlo,growhi,dir-1)

         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          !dir=1..sdim
          call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir-1)

          im1=i-ii
          jm1=j-jj
          km1=k-kk

           ! kluge to prevent reading out of the array
           ! bounds for coupling terms that have solid cells
           ! in the stencil.  
           ! The corner values will never be used if cell pair is
           ! owned by a different material (or not owned at all)
           ! from the material owned on the face.
           ! If the cell pair is owned by a prescribed solid, then
           ! it is not used in the coupling stencil.
          shift_flag=0
          ishift=i
          jshift=j
          kshift=k
          if ((dir.eq.1).or. &
              (dir.eq.3)) then
           dirtan=2
           if (j.eq.fablo(dirtan)-1) then
            jshift=j+1
            shift_flag=1
           else if (j.eq.fabhi(dirtan)+1) then
            jshift=j-1
            shift_flag=1
           else if ((j.ge.fablo(dirtan)).and. &
                    (j.le.fabhi(dirtan))) then
            ! do nothing
           else
            print *,"j invalid"
            stop
           endif
          endif

          if ((dir.eq.2).or. &
              (dir.eq.3)) then
           dirtan=1
           if (i.eq.fablo(dirtan)-1) then
            ishift=i+1
            shift_flag=1
           else if (i.eq.fabhi(dirtan)+1) then
            ishift=i-1
            shift_flag=1
           else if ((i.ge.fablo(dirtan)).and. &
                    (i.le.fabhi(dirtan))) then
            ! do nothing
           else
            print *,"i invalid"
            stop
           endif
          endif

          if (((dir.eq.1).or. &
               (dir.eq.2)).and.  &
              (SDIM.eq.3)) then
           dirtan=SDIM
           if (k.eq.fablo(dirtan)-1) then
            kshift=k+1
            shift_flag=1
           else if (k.eq.fabhi(dirtan)+1) then
            kshift=k-1
            shift_flag=1
           else if ((k.ge.fablo(dirtan)).and. &
                    (k.le.fabhi(dirtan))) then
            ! do nothing
           else
            print *,"k invalid fort_face_gradients"
            stop
           endif
          endif

          side=1
          index_adjoin(side,1)=im1
          index_adjoin(side,2)=jm1
          index_adjoin(side,3)=km1
          side=2
          index_adjoin(side,1)=i
          index_adjoin(side,2)=j
          index_adjoin(side,3)=k

          mask_boundary=0

          do dir2=1,SDIM
           do side=1,2
            if (side.eq.1) then
             ! mask3=tag at exterior fine/fine border.
             ! mask3=1-tag at other exterior boundaries.
             ! mask0=tag if not covered by level+1 or outside the domain.
             masktest=NINT(mask3(D_DECL(im1,jm1,km1)))
             maskcov=NINT(mask0(D_DECL(im1,jm1,km1)))
            else if (side.eq.2) then
             masktest=NINT(mask3(D_DECL(i,j,k)))
             maskcov=NINT(mask0(D_DECL(i,j,k)))
            else
             print *,"side invalid"
             stop
            endif
            if (maskcov.eq.1) then
             ! do nothing
            else if (maskcov.eq.0) then
             mask_boundary=1 ! do not include covered cells in the stencil.
            else
             print *,"maskcov invalid"
             stop
            endif

            icrit=index_adjoin(side,dir2)

            if (icrit.lt.fablo(dir2)) then
             if (masktest.eq.0) then ! not fine/fine
              do nc=1,SDIM
               bctest=velbc(dir2,1,nc)
               if ((bctest.eq.INT_DIR).or. &!coarse/fine not fine/fine,periodic
                   (bctest.eq.FOEXTRAP).or. &
                   (bctest.eq.EXT_DIR).or. &
                   (bctest.eq.REFLECT_ODD).or. &
                   (bctest.eq.REFLECT_EVEN)) then
                mask_boundary=1
               else
                print *,"bctest invalid"
                stop
               endif
              enddo ! nc
             else if (masktest.eq.1) then
              ! do nothing fine-fine border
             else
              print *,"masktest invalid"
              stop
             endif
            endif
            if (icrit.gt.fabhi(dir2)) then
             if (masktest.eq.0) then
              do nc=1,SDIM
               bctest=velbc(dir2,2,nc)
               if ((bctest.eq.INT_DIR).or. &
                   (bctest.eq.FOEXTRAP).or. &
                   (bctest.eq.EXT_DIR).or. &
                   (bctest.eq.REFLECT_ODD).or. &
                   (bctest.eq.REFLECT_EVEN)) then
                mask_boundary=1
               else
                print *,"bctest invalid"
                stop
               endif
              enddo ! nc
             else if (masktest.eq.1) then
              ! do nothing fine-fine border
             else
              print *,"masktest invalid"
              stop
             endif
            endif
           enddo ! side
          enddo ! dir2
            
          mdata(D_DECL(i,j,k),dir)=one

          do im=1,num_materials
           lsleft(im)=levelpc(D_DECL(im1,jm1,km1),im)
           lsright(im)=levelpc(D_DECL(i,j,k),im)
          enddo ! im=1..num_materials

          call get_primary_material(lsleft,imL)
          call get_primary_material(lsright,imR)

          call gridsten_level(xclamped_minus_sten,im1,jm1,km1,level, &
                  nhalfclamped)
          call gridsten_level(xclamped_plus_sten,i,j,k,level, &
                  nhalfclamped)
          do dir2=1,SDIM
           xclamped_minus(dir2)=xclamped_minus_sten(0,dir2)
           xclamped_plus(dir2)=xclamped_plus_sten(0,dir2)
          enddo

           ! LS>0 if clamped
          call SUB_clamped_LS(xclamped_minus,time,LS_clamped_minus, &
           vel_clamped_minus,temperature_clamped_minus,prescribed_flag,dx)
          call SUB_clamped_LS(xclamped_plus,time,LS_clamped_plus, &
           vel_clamped_plus,temperature_clamped_plus,prescribed_flag,dx)
          if ((LS_clamped_minus.ge.zero).or. &
              (LS_clamped_plus.ge.zero)) then
           is_clamped_face=1
           do dir2=1,SDIM
            if (LS_clamped_minus.lt.zero) then
             vel_clamped_face(dir2)=vel_clamped_plus(dir2)
             is_clamped_face=2
            else if (LS_clamped_plus.lt.zero) then
             vel_clamped_face(dir2)=vel_clamped_minus(dir2)
             is_clamped_face=3
            else
             vel_clamped_face(dir2)=half*(vel_clamped_plus(dir2)+ &
               vel_clamped_minus(dir2))
            endif
           enddo
          else if ((LS_clamped_minus.lt.zero).and. &
                   (LS_clamped_plus.lt.zero)) then
           is_clamped_face=0
          else
           print *,"LS_clamped plus or minus is NaN"
           stop
          endif


          if ((is_clamped_face.eq.1).or. &
              (is_clamped_face.eq.2).or. &
              (is_clamped_face.eq.3)) then
           faceLS(D_DECL(i,j,k),dir)=zero ! mask off this gradient
          else if (is_clamped_face.eq.0) then
           if (mask_boundary.eq.0) then
            if (imL.eq.imR) then
             faceLS(D_DECL(i,j,k),dir)=imL
            else
             faceLS(D_DECL(i,j,k),dir)=zero ! mask off this gradient
            endif
           else if (mask_boundary.eq.1) then
            faceLS(D_DECL(i,j,k),dir)=zero ! mask off this gradient
           else
            print *,"mask_boundary invalid"
            stop
           endif
          else
           print *,"is_clamped_face invalid"
           stop
          endif

          indexleft=index_adjoin(1,dir)
          indexright=index_adjoin(2,dir)

          partidL=0
          partidR=0

          if (is_rigid(imL).eq.1) then
           do im=1,imL-1
            if (is_lag_part(im).eq.1) then
             partidL=partidL+1
            else if (is_lag_part(im).eq.0) then
             ! do nothing
            else
             print *,"is_lag_part(im) invalid"
             stop
            endif
           enddo ! im=1..imL-1
           if (im_solid_map(partidL+1)+1.ne.imL) then
            print *,"im_solid_map(partidL+1)+1.ne.imL"
            stop
           endif
          else if (is_rigid(imL).eq.0) then
           ! do nothing 
          else
           print *,"is_rigid(imL) invalid"
           stop
          endif

          if (is_rigid(imR).eq.1) then
           do im=1,imR-1
            if (is_lag_part(im).eq.1) then
             partidR=partidR+1
            else if (is_lag_part(im).eq.0) then
             ! do nothing
            else
             print *,"is_lag_part(im) invalid"
             stop
            endif
           enddo ! im=1..imR-1
           if (im_solid_map(partidR+1)+1.ne.imR) then
            print *,"im_solid_map(partidR+1)+1.ne.imR"
            stop
           endif
          else if (is_rigid(imR).eq.0) then
           ! do nothing 
          else
           print *,"is_rigid(imR) invalid"
           stop
          endif
 
          do nc=1,SDIM

           velcomp=nc
           left_rigid=0
           right_rigid=0

           if (is_rigid(imL).eq.1) then
            left_rigid=1

            if (dir.eq.1) then
             solidvelleft(nc)= &
               solidx(D_DECL(ishift,jshift,kshift),partidL*SDIM+nc)
            else if (dir.eq.2) then
             solidvelleft(nc)= &
               solidy(D_DECL(ishift,jshift,kshift),partidL*SDIM+nc)
            else if ((dir.eq.3).and.(SDIM.eq.3)) then
             solidvelleft(nc)= &
               solidz(D_DECL(ishift,jshift,kshift),partidL*SDIM+nc)
            else
             print *,"dir invalid"
             stop
            endif
           else if (is_rigid(imL).eq.0) then
            solidvelleft(nc)=zero
           else
            print *,"is_rigid(imL) invalid"
            stop
           endif

           if (is_rigid(imR).eq.1) then
            right_rigid=1

            if (dir.eq.1) then
             solidvelright(nc)= &
               solidx(D_DECL(ishift,jshift,kshift),partidR*SDIM+nc)
            else if (dir.eq.2) then
             solidvelright(nc)= &
               solidy(D_DECL(ishift,jshift,kshift),partidR*SDIM+nc)
            else if ((dir.eq.3).and.(SDIM.eq.3)) then
             solidvelright(nc)= &
               solidz(D_DECL(ishift,jshift,kshift),partidR*SDIM+nc)
            else
             print *,"dir invalid"
             stop
            endif

           else if (is_rigid(imR).eq.0) then
            solidvelright(nc)=zero
           else
            print *,"is_rigid(imR) invalid"
            stop
           endif

           if ((indexleft.lt.fablo(dir)).and. &
               (velbc(dir,1,nc).eq.EXT_DIR).and. &
               (is_rigid(imL).eq.0)) then
            left_rigid=1
            solidvelleft(nc)=vel(D_DECL(im1,jm1,km1),velcomp)
           endif

           if ((indexright.gt.fabhi(dir)).and. &
               (velbc(dir,2,nc).eq.EXT_DIR).and. &
               (is_rigid(imR).eq.0)) then
            right_rigid=1
            solidvelright(nc)=vel(D_DECL(i,j,k),velcomp)
           endif

           if (is_clamped_face.eq.1) then
            left_rigid=1
            right_rigid=1
            solidvelleft(nc)=vel_clamped_face(nc)
            solidvelright(nc)=vel_clamped_face(nc)
           else if (is_clamped_face.eq.2) then
            right_rigid=1
            solidvelright(nc)=vel_clamped_face(nc)
           else if (is_clamped_face.eq.3) then
            left_rigid=1
            solidvelleft(nc)=vel_clamped_face(nc)
           else if (is_clamped_face.eq.0) then
            ! do nothing
           else
            print *,"is_clamped_face invalid"
            stop
           endif

           if (homflag.eq.1) then
            solidvelleft(nc)=zero
            solidvelright(nc)=zero
           else if (homflag.eq.0) then
            ! do nothing
           else
            print *,"homflag invalid 1"
            stop
           endif

           ! lsleft_solid and lsright_solid < 0 in the solid.
           if ((left_rigid.eq.1).and.(right_rigid.eq.1)) then
            mdata(D_DECL(i,j,k),dir)=zero
           else if ((left_rigid.eq.1).or.(right_rigid.eq.1)) then
            if (shift_flag.eq.0) then
             ! do nothing

             !prevent reading solidvel outside the array bounds
            else if (shift_flag.eq.1) then 
             mdata(D_DECL(i,j,k),dir)=zero
            else
             print *,"shift_flag invalid"
             stop
            endif
           else if ((left_rigid.eq.0).and.(right_rigid.eq.0)) then
            ! do nothing
           else
            print *,"left_rigid or right_rigid invalid"
            stop
           endif

           ! extra factor of r for theta gradient in cylindrical coordinates.
           RR=one
           if (rzflag.eq.COORDSYS_CARTESIAN) then
            ! do nothing
           else if (rzflag.eq.COORDSYS_RZ) then
            if (SDIM.ne.2) then
             print *,"dimension bust"
             stop
            endif
            ! do nothing
           else if (rzflag.eq.COORDSYS_CYLINDRICAL) then
            if (dir.eq.2) then ! theta direction  s_theta/r
             RR=xstenMAC(0,1)
            endif
           else
            print *,"rzflag invalid"
            stop
           endif
           if (RR.le.zero) then
            print *,"RR invalid"
            stop
           else if (RR.gt.zero) then
            delta=RR*(xstenMAC(1,dir)-xstenMAC(-1,dir))
           else
            print *,"RR bust"
            stop
           endif
           if (delta.le.zero) then
            print *,"delta invalid"
            stop
           endif

           theta_factor=one

           if ((left_rigid.eq.0).and.(right_rigid.eq.0)) then
            vplus(nc)=vel(D_DECL(i,j,k),velcomp)
            vminus(nc)=vel(D_DECL(im1,jm1,km1),velcomp)
           else if ((left_rigid.eq.1).and.(right_rigid.eq.1)) then
            vplus(nc)=solidvelright(nc)
            vminus(nc)=solidvelleft(nc)
            theta_factor=zero ! zero out the gradient.
           else if ((left_rigid.eq.0).and.(right_rigid.eq.1)) then
            vplus(nc)=solidvelright(nc)
            vminus(nc)=vel(D_DECL(im1,jm1,km1),velcomp)
            theta_factor=two
           else if ((left_rigid.eq.1).and.(right_rigid.eq.0)) then
            vplus(nc)=vel(D_DECL(i,j,k),velcomp)
            vminus(nc)=solidvelleft(nc)
            theta_factor=two
           else
            print *,"left_rigid or right_rigid invalid"
            stop
           endif

 ! grad u in cylindrical coordinates:
 !
 ! S= (grad u + grad u^T)/2 
 !
 ! grad u=| u_r  u_t/r-v/r  u_z  |
 !        | v_r  v_t/r+u/r  v_z  |
 !        | w_r  w_t/r      w_z  |
 !
 ! S=
 !   |u_r     (u_t/r+v_r-v/r)/2   (u_z+w_r)/2   |
 !   | .      v_t/r+u/r           (v_z+w_t/r)/2 |
 !   | .      .                   w_z           |
 !
 ! note: v_r-v/r=r(v/r)_r
 !
 ! 2S=
 !
 !   |2u_r     (u_t/r+v_r-v/r)   (u_z+w_r)   |
 !   | .       2v_t/r+2u/r       (v_z+w_t/r) |
 !   | .      .                    2w_z      |
 !
 ! 
 ! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
 !         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
 !         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |
 ! 
 ! ur =  costheta u + sintheta v
 ! ut = -sintheta u + costheta v
 ! 
 ! u = costheta ur - sintheta ut
 ! v = sintheta ur + costheta ut
 !
 ! e.g. theta=pi/2  ur=0 
 !   u=-utheta  v=0
 ! if constant viscosity:
 ! div u = (ru)_r/r + v_t/r + w_z= u_r + u/r +v_t/r + w_z=0
 ! (div u)_r=u_rr+u_r/r-u/r^2+v_tr/r-v_t/r^2+w_zr=0
 ! (div u)_t/r=u_rt/r+u_t/r^2+v_tt/r^2+w_zt/r=0
 ! (div u)_z=u_rz+u_z/r+v_tz/r+w_zz=0
 !
 ! div(2 S)=
 ! |2u_rr+2u_r/r+u_tt/r^2+v_tr/r-v_t/r^2-2v_t/r^2-2u/r^2+u_zz+w_rz |
 ! |u_tr/r+v_rr+v_r/r-v_r/r+2v_tt/r^2+2u_t/r^2+u_t/r^2+v_r/r-v/r^2+v_zz+w_tz/r|
 ! |u_zr+u_z/r+w_rr+w_r/r+v_zt/r+w_tt/r^2 + 2w_zz |=
 !
 ! |u_rr+u_r/r-u/r^2+u_tt/r^2-2v_t/r^2+u_zz |    
 ! |v_rr+v_r/r+v_tt/r^2+2u_t/r^2-v/r^2+v_zz |
 ! |w_rr+w_r/r+w_tt/r^2+w_zz                |
 !
 ! compromise: 
 !
 ! GU=| u_r       u_t/r  u_z  |
 !    | v_r       v_t/r  v_z  |
 !    | w_r       w_t/r  w_z  |
 !
 ! hoop term 1st component:  -3 v_t/r^2 - 2 u/r^2
 ! hoop term 2nd component:   3 u_t/r^2 - v/r^2
 ! 
 ! If uncoupled_viscosity==true:
 ! hoop term 1st component:  -2 v_t/r^2 - u/r^2
 ! hoop term 2nd component:   2 u_t/r^2 - v/r^2
 ! No coupling terms.
 ! Diagonal terms not multiplied by 2.

           hold_grad=theta_factor*(vplus(nc)-vminus(nc))/delta

           tensorcomponent=nbase+nc

           tdata(D_DECL(i,j,k),tensorcomponent)=hold_grad

          enddo ! nc=1..sdim

         enddo
         enddo
         enddo  ! i,j,k faces

        else if (spectral_loop.eq.1) then
         ! do nothing
        else
         print *,"spectral_loop invalid"
         stop
        endif

       else if (tileloop.eq.1) then ! high order (see below)
        ! do nothing
       else
        print *,"tileloop invalid"
        stop
       endif

      else if (itensor_iter.eq.1) then  ! cell grad U

       if (spectral_loop.eq.0) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

       if (tileloop.eq.0) then ! low order gradients at cells

         ! same as growntilebox, except includes one layer of ghost
         ! cells in the tangential directions.
        call growntileboxTENSOR_SEM( &
         tilelo,tilehi,fablo,fabhi,growlo,growhi,dir-1) 

        do k=growlo(3),growhi(3)
        do j=growlo(2),growhi(2)
        do i=growlo(1),growhi(1)
         call gridsten_level(xsten,i,j,k,level,nhalfcell)
         do dir2=1,SDIM
          xclamped_cen(dir2)=xsten(0,dir2)
         enddo

          ! LS>0 if clamped
         call SUB_clamped_LS(xclamped_cen,time,LS_clamped_cen, &
           vel_clamped_cen,temperature_clamped_cen,prescribed_flag,dx)

         im1=i-ii
         jm1=j-jj
         km1=k-kk
         ip1=i+ii
         jp1=j+jj
         kp1=k+kk
         do im=1,num_materials
          lspoint(im)=levelpc(D_DECL(i,j,k),im)
         enddo
         call get_primary_material(lspoint,im_primary)
         if ((im_primary.lt.1).or.(im_primary.gt.num_materials)) then
          print *,"im_primary invalid"
          stop
         endif

         if (LS_clamped_cen.ge.zero) then
          do nc=1,SDIM
           gradu(nc)=zero
          enddo
         else if (LS_clamped_cen.lt.zero) then

          if (is_rigid(im_primary).eq.1) then
           do nc=1,SDIM
            gradu(nc)=zero
           enddo
          else if (is_rigid(im_primary).eq.0) then
           int_xlo=max(xsten(-2,dir),xsten(-1,dir))
           int_xhi=min(xsten(0,dir),xsten(1,dir))
           if (int_xhi.gt.int_xlo) then
            leftwt=int_xhi-int_xlo
           else
            leftwt=zero
           endif
           int_xlo=max(xsten(0,dir),xsten(-1,dir))
           int_xhi=min(xsten(2,dir),xsten(1,dir))
           if (int_xhi.gt.int_xlo) then
            rightwt=int_xhi-int_xlo
           else
            rightwt=zero
           endif
           if ((leftwt.le.zero).or.(rightwt.le.zero)) then
            print *,"weights invalid"
            stop
           endif
           lsleft(im_primary)=levelpc(D_DECL(im1,jm1,km1),im_primary)
           lsright(im_primary)=levelpc(D_DECL(ip1,jp1,kp1),im_primary)
            
           do nc=1,SDIM
            tensorcomponent=nbase+nc
            slopeLT=tdata(D_DECL(i,j,k),tensorcomponent)
            slopeRT=tdata(D_DECL(ip1,jp1,kp1),tensorcomponent)

            if (((lsleft(im_primary).ge.zero).and. &
                 (lsright(im_primary).ge.zero)).or. &
                ((lsleft(im_primary).lt.zero).and. &
                 (lsright(im_primary).lt.zero))) then
             gradu(nc)=(leftwt*slopeLT+rightwt*slopeRT)/(leftwt+rightwt)
            else if (lsleft(im_primary).ge.zero) then
             gradu(nc)=slopeLT
            else if (lsright(im_primary).ge.zero) then
             gradu(nc)=slopeRT
            else
             print *,"lsleft or lsright invalid"
             stop
            endif

           enddo ! nc=1..sdim
          else
           print *,"is_rigid(im_primary) invalid"
           stop
          endif

         else
          print *,"LS_clamped_cen is NaN"
          stop
         endif
            
         do nc=1,SDIM
          tensorcomponent=nbase+nc
          c_tdata(D_DECL(i,j,k),tensorcomponent)=gradu(nc)
         enddo ! nc=1..sdim

        enddo ! k
        enddo ! j
        enddo ! i

       else if (tileloop.eq.1) then ! high order below
        ! do nothing
       else
        print *,"tileloop invalid"
        stop
       endif

      else
       print *,"itensor_iter invalid"
       stop
      endif

       ! in: fort_face_gradients
      if (enable_spectral.eq.1) then

       if (bfact.ge.2) then

        if (tileloop.eq.1) then ! tileloop==1 => high order

         if (itensor_iter.eq.0) then ! tensor on MAC grid.

          if ((dir.lt.1).or.(dir.gt.SDIM)) then
           print *,"dir invalid face gradients 3"
           stop
          endif

          ! same as growntilebox, except includes one layer of ghost
          ! cells in the tangential directions.
          call growntileboxTENSOR_SEM(tilelo,tilehi,fablo,fabhi, &
           growlo,growhi,dir-1)

          do k=growlo(3),growhi(3)
          do j=growlo(2),growhi(2)
          do i=growlo(1),growhi(1)
           local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
           maskcov=NINT(mask0(D_DECL(i,j,k)))

           if ((local_maskSEM.ge.1).and. &
               (local_maskSEM.le.num_materials).and. &
               (maskcov.eq.1)) then

            call strip_status_dir(i,j,k,bfact,dir-1,stripstat)

             ! stripstat==1 if (idx_dir/bfact)*bfact == idx_dir
             ! In otherwords, one can have stripstat==1 even if
             ! one of the tangential indexes is outside of the box.
            if (stripstat.eq.1) then

              ielem=i
              jelem=j
              kelem=k

              scomp=1  

              dcomp=nbase+1

              ncomp_source=SDIM
              ncomp_dest=SDIM
              ncomp_xgp=AMREX_SPACEDIM_SQR
              ncomp_xp=SDIM ! number of amrsync components
              scomp_bc=1
              operation_flag=OP_UGRAD_MAC  ! evaluate tensor values
              energyflag=SUB_OP_DEFAULT
              project_option_vel=SOLVETYPE_VISC
              def_dt=one

               ! in: fort_face_gradients
               ! the boundary conditions for "vel" are already set in the
               ! ghost cell.  e.g. homogeneous versus inhomogeneous.
              call SEM_CELL_TO_MAC( &
               ncomp_xp, &  ! number of amrsync components
               simple_AMR_BC_flag_viscosity, &
               level, &
               finest_level, &
               operation_flag, &  ! OP_UGRAD_MAC
               energyflag, &
               project_option_vel, &
               time, &  ! beta
               time, &  ! visc_coef
               time, &
               def_dt, &  ! dt
               ielem,jelem,kelem, &
               tilelo,tilehi, &
               fablo,fabhi, &
               xlo, &
               dx, &
               dir, &
               bfact,bfact_c,bfact_f, &
               velbc, &  ! presbc
               velbc, &
               scomp, &
               scomp_bc, &
               dcomp, &
               ncomp_dest, &
               ncomp_source, &
               ncomp_xgp, &     ! =AMREX_SPACEDIM_SQR
               ncomp_dest, &    ! ncphys
               spectral_loop, &
               ncfluxreg, &
               semflux_ptr, &
               mask3_ptr, &
               mask0_ptr, & !mask0=1 if not cov. by finer or outside.
               vel_ptr, &
               vel_ptr, &  ! pres
               vel_ptr, &  ! den
               tdata_ptr, &  !xface
               tdata_ptr, &  !xgp (destination)
               tdata_ptr, &  !xcut
               amrsync_ptr, & !xp
               tdata_ptr, &  !xvel
               maskSEM_ptr)

            else if (stripstat.eq.0) then
             ! do nothing
            else
             print *,"stripstat invalid"
             stop
            endif
         
           else if ((local_maskSEM.eq.0).or. &
                    (maskcov.eq.0)) then

            ! do nothing

           else
            print *,"local_maskSEM invalid"
            stop
           endif 

          enddo
          enddo
          enddo ! i,j,k

         else if (itensor_iter.eq.1) then ! cell grad U

          if (spectral_loop.eq.0) then
           ! do nothing
          else
           print *,"spectral_loop invalid"
           stop
          endif

          if ((dir.lt.1).or.(dir.gt.SDIM)) then
           print *,"dir invalid face gradients 4"
           stop
          endif

           ! same as growntilebox, except includes one layer of ghost
           ! cells in the tangential directions.
          call growntileboxTENSOR_SEM(tilelo,tilehi,fablo,fabhi, &
           growlo,growhi,dir-1)

          do k=growlo(3),growhi(3)
          do j=growlo(2),growhi(2)
          do i=growlo(1),growhi(1)

           local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
           maskcov=NINT(mask0(D_DECL(i,j,k)))

           if ((local_maskSEM.ge.1).and. &
               (local_maskSEM.le.num_materials).and. &
               (maskcov.eq.1)) then

            call strip_status_dir(i,j,k,bfact,dir-1,stripstat)

             ! stripstat==1 if (idx_dir/bfact)*bfact == idx_dir
            if (stripstat.eq.1) then

              ielem=i
              jelem=j
              kelem=k

              dcomp=nbase+1

              scomp=dcomp

              scomp_bc=1
              ncomp_dest=SDIM ! ncomp
              ncomp_xvel=AMREX_SPACEDIM_SQR
              ncomp_denold=AMREX_SPACEDIM_SQR
              ncomp_veldest=AMREX_SPACEDIM_SQR
              ncomp_dendest=AMREX_SPACEDIM_SQR
              ncomp_cterm=AMREX_SPACEDIM_SQR
              operation_flag=OP_GRADU_MAC_TO_CELL 
              project_option_vel=SOLVETYPE_VISC
              energyflag=SUB_OP_DEFAULT
              def_dt=one

               ! 1<=dir<=sdim
              call SEM_MAC_TO_CELL( &
               ncomp_denold, &
               ncomp_veldest, &
               ncomp_dendest, &
               ns_time_order, &
               divu_outer_sweeps, &
               num_divu_outer_sweeps, &
               SDC_outer_sweeps, &
               level, &
               finest_level, &
               operation_flag, & ! operation_flag==OP_GRADU_MAC_TO_CELL
               project_option_vel, &
               energyflag, &
               energyflag, & ! homflag
               local_maskSEM, &
               time, &
               slab_step, &
               def_dt, & ! dt
               ielem,jelem,kelem, &
               tilelo,tilehi, &
               fablo,fabhi, &
               xlo,dx,dir,bfact, &
               velbc, &
               velbc, & ! presbc
               scomp, &
               scomp_bc, &
               dcomp, &
               ncomp_dest, & ! ncomp
               ncomp_xvel, &
               ncomp_cterm, &
               tdata_ptr, & ! vol
               tdata_ptr, & ! xface
               tdata_ptr, & ! xvel
               tdata_ptr, & ! maskcoef
               tdata_ptr, & ! cterm
               tdata_ptr, & ! mdotcell
               tdata_ptr, & ! maskdivres
               tdata_ptr, & ! pold
               tdata_ptr, & ! denold
               tdata_ptr, & ! ustar
               c_tdata_ptr, & ! veldest
               c_tdata_ptr, & ! dendest
               c_tdata_ptr) !divdest

            else if (stripstat.eq.0) then
             ! do nothing
            else
             print *,"stripststat invalid"
             stop
            endif
  
           else if ((local_maskSEM.eq.0).or. &
                    (maskcov.eq.0)) then
            ! do nothing
           else
            print *,"local_maskSEM invalid"
            stop
           endif

          enddo
          enddo
          enddo

         else
          print *,"itensor_iter invalid"
          stop
         endif 

        else if (tileloop.eq.0) then
         ! do nothing
        else
         print *,"tileloop invalid"
         stop
        endif

       else if (bfact.eq.1) then
        ! do nothing
       else
        print *,"bfact invalid78"
        stop
       endif

      else if (enable_spectral.eq.0) then
       ! do nothing
      else
       print *,"enable_spectral invalid"
       stop
      endif

      return 
      end subroutine fort_face_gradients


! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
!         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
!         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |
       ! The continuum model for FSI:
       !  du/dt + u dot grad u = -grad p/rho + div(H tau_elastic)/rho
       !  H=1 in the elastic solid
       !  H=0 in the fluid
       ! fort_elastic_force is called from NavierStokes2.cpp:
       ! void NavierStokes::MAC_GRID_ELASTIC_FORCE
      subroutine fort_elastic_force( &
       im_elastic, & ! 0..num_materials-1
       partid, & ! 0..num_materials_viscoelastic-1
       force_dir, & ! 0..sdim-1
       ncomp_visc, &
       visc_coef, &
       velbc, &
       dt, &
       cur_time, &
       xlo,dx, &
       grid_type_CC, &
       MACFLUX_CC, &
       DIMS(MACFLUX_CC), &
       grid_type_X, &
       MACFLUX_X, &
       DIMS(MACFLUX_X), &
       grid_type_Y, &
       MACFLUX_Y, &
       DIMS(MACFLUX_Y), &
       grid_type_Z, &
       MACFLUX_Z, &
       DIMS(MACFLUX_Z), &
       visc,DIMS(visc), &
       mask,DIMS(mask), &  ! 1=fine/fine 0=coarse/fine
       maskcoef,DIMS(maskcoef), & ! 1=not cov by level+1 or outside.
       levelpc,DIMS(levelpc), &
       rhoinvfab, &
       DIMS(rhoinvfab), &
       SNEW, &
       DIMS(SNEW), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       rzflag, &
       domlo,domhi) &
      bind(c,name='fort_elastic_force')

      use probcommon_module
      use global_utility_module
      use MOF_routines_module
      use mass_transfer_module
 
      IMPLICIT NONE

      integer, INTENT(in) :: im_elastic !0..num_materials-1
      integer, INTENT(in) :: partid !0..num_materials_viscoelastic-1

       ! MAC force component, force_dir=0..sdim-1
      integer, INTENT(in) :: force_dir  
      integer, INTENT(in) :: ncomp_visc
      real(amrex_real), INTENT(in) :: visc_coef
      integer, INTENT(in) :: velbc(SDIM,2,SDIM) 
      real(amrex_real), INTENT(in) :: dt 
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM) 
      integer, INTENT(in) :: grid_type_CC
      integer, INTENT(in) :: DIMDEC(MACFLUX_CC)
      integer, INTENT(in) :: grid_type_X
      integer, INTENT(in) :: DIMDEC(MACFLUX_X)
      integer, INTENT(in) :: grid_type_Y
      integer, INTENT(in) :: DIMDEC(MACFLUX_Y)
      integer, INTENT(in) :: grid_type_Z
      integer, INTENT(in) :: DIMDEC(MACFLUX_Z)
      integer, INTENT(in) :: DIMDEC(visc)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: DIMDEC(maskcoef)
      integer, INTENT(in) :: DIMDEC(levelpc)
      integer, INTENT(in) :: DIMDEC(rhoinvfab)
      integer, INTENT(in) :: DIMDEC(SNEW)

      real(amrex_real), INTENT(in), target :: MACFLUX_CC( &
             DIMV(MACFLUX_CC), &
             ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: MACFLUX_CC_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: MACFLUX_X( &
             DIMV(MACFLUX_X), &
             ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: MACFLUX_X_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: MACFLUX_Y( &
             DIMV(MACFLUX_Y), &
             ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: MACFLUX_Y_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: MACFLUX_Z( &
             DIMV(MACFLUX_Z), &
             ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: MACFLUX_Z_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: visc(DIMV(visc),ncomp_visc)
      real(amrex_real), pointer :: visc_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in), target :: maskcoef(DIMV(maskcoef))
      real(amrex_real), pointer :: maskcoef_ptr(D_DECL(:,:,:))

      real(amrex_real), target, INTENT(in) :: &
              levelpc(DIMV(levelpc),num_materials*(1+SDIM))
      real(amrex_real), pointer :: levelpc_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: rhoinvfab(DIMV(rhoinvfab))
      real(amrex_real), pointer :: rhoinvfab_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(inout), target :: SNEW(DIMV(SNEW))
      real(amrex_real), pointer :: SNEW_ptr(D_DECL(:,:,:))

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: rzflag 
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      integer :: i,j,k
      integer :: ii,jj,kk
      integer :: dir_flux,side_flux
      integer, PARAMETER :: nhalf=3
      real(amrex_real) :: xstenCELL(-nhalf:nhalf,SDIM)
      real(amrex_real) :: xsten_flux(-nhalf:nhalf,SDIM)
      integer dircomp
      integer dir_row,dir_column
      real(amrex_real), target :: x_CELL_control_volume(SDIM)
      real(amrex_real) DISP_TEN(3,3)
      real(amrex_real) dxmin
      integer im_elastic_p1
      integer im_LS

      integer dir_local
      real(amrex_real) :: LS_control_volume(num_materials)
      real(amrex_real) :: LS_outside(num_materials)

      real(amrex_real) xflux_local(0:1,3,3)
      real(amrex_real) yflux_local(0:1,3,3)
      real(amrex_real) zflux_local(0:1,3,3)
      real(amrex_real) center_flux(3,3)
      real(amrex_real) center_hoop_22

      integer mask_flux_point(2,SDIM)
      real(amrex_real) x_at_flux_point(2,SDIM,SDIM)
      real(amrex_real) x_at_flux_point_local(SDIM)
      real(amrex_real) hx,hy,hz,rplus,rminus,rval
      real(amrex_real) LS_clamped
      real(amrex_real) vel_clamped(SDIM)
      real(amrex_real) temperature_clamped
      integer prescribed_flag
      integer local_mask
      integer mask_control_volume
      real(amrex_real) force(SDIM)
      real(amrex_real) bodyforce
      real(amrex_real) deninv
      real(amrex_real) XFORCE_local
      integer box_type_flux(SDIM)
      integer grid_type_flux
      integer grid_type_sanity
      integer iflux_array(SDIM)
      integer iflux,jflux,kflux
      integer itensor
      real(amrex_real) local_tensor_data(ENUM_NUM_TENSOR_TYPE)

      SNEW_ptr=>SNEW

      im_elastic_p1=im_elastic+1

      if (ENUM_NUM_TENSOR_TYPE.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ENUM_NUM_TENSOR_TYPE invalid"
       stop
      endif

      MACFLUX_CC_ptr=>MACFLUX_CC
      call checkbound_array(fablo,fabhi,MACFLUX_CC_ptr,1,grid_type_CC)

      MACFLUX_X_ptr=>MACFLUX_X
      call checkbound_array(fablo,fabhi,MACFLUX_X_ptr,0,grid_type_X)

      MACFLUX_Y_ptr=>MACFLUX_Y
      call checkbound_array(fablo,fabhi,MACFLUX_Y_ptr,0,grid_type_Y)

      MACFLUX_Z_ptr=>MACFLUX_Z
      call checkbound_array(fablo,fabhi,MACFLUX_Z_ptr,0,grid_type_Z)

      visc_ptr=>visc
      call checkbound_array(fablo,fabhi,visc_ptr,1,-1)
      mask_ptr=>mask
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      maskcoef_ptr=>maskcoef
      call checkbound_array1(fablo,fabhi,maskcoef_ptr,1,-1)
      levelpc_ptr=>levelpc
      call checkbound_array(fablo,fabhi,levelpc_ptr,2,-1)

      rhoinvfab_ptr=>rhoinvfab
      call checkbound_array1(fablo,fabhi,rhoinvfab_ptr,0,-1)

      call checkbound_array1(fablo,fabhi,SNEW_ptr,1,-1)

      call get_dxmin(dx,bfact,dxmin)
      if (dxmin.gt.zero) then
       ! do nothing
      else
       print *,"dxmin invalid"
       stop
      endif

      if ((im_elastic.ge.0).and. &
          (im_elastic.lt.num_materials)) then
       ! do nothing
      else
       print *,"im_elastic invalid"
       stop
      endif
      if ((partid.ge.0).and. &
          (partid.lt.num_materials_viscoelastic)) then
       ! do nothing
      else
       print *,"partid invalid"
       stop
      endif
      if ((level.ge.0).and. &
          (level.le.finest_level)) then
       ! do nothing
      else
       print *,"level invalid"
       stop
      endif
      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid"
       stop
      endif
      if ((force_dir.ge.0).and.(force_dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"force_dir invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

       ! traverse the cell centered grid.
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten_level(xstenCELL,i,j,k,level,nhalf)

       do dircomp=1,SDIM
        x_CELL_control_volume(dircomp)=xstenCELL(0,dircomp)
       enddo

        !rzflag=0 => volume and area independent of r.
       if (rzflag.eq.COORDSYS_CARTESIAN) then
        !volume=dx dy dz
        rval=one
       else if ((rzflag.eq.COORDSYS_RZ).or. &
                (rzflag.eq.COORDSYS_CYLINDRICAL)) then
        !volume=(dtheta/2)*(r_{2}^2 - r_{1}^{2}) dz=dtheta*dz*dr r
        rval=x_CELL_control_volume(1)
       else
        print *,"rzflag invalid"
        stop
       endif

       if (rval.gt.zero) then
        ! do nothing
       else
        print *,"rval invalid"
        stop
       endif

       do im_LS=1,num_materials
        LS_control_volume(im_LS)=levelpc(D_DECL(i,j,k),im_LS)
       enddo

       call get_primary_material(LS_control_volume,local_mask)

       if (local_mask.eq.im_elastic_p1) then
        local_mask=1
       else if ((local_mask.ge.1).and. &
                (local_mask.le.num_materials)) then
        local_mask=0
       else
        print *,"local_mask invalid"
        stop
       endif

        ! LS>0 if clamped
       call SUB_clamped_LS(x_CELL_control_volume,cur_time,LS_clamped, &
         vel_clamped,temperature_clamped,prescribed_flag,dx)
       if (LS_clamped.ge.zero) then
        local_mask=0
       else if (LS_clamped.lt.zero) then
        ! do nothing
       else
        print *,"LS_clamped invalid"
        stop
       endif

       XFORCE_local=zero
           
       mask_control_volume=local_mask

        ! im_elastic_p1 dominates the center of the MAC control volume.
       if (mask_control_volume.eq.1) then

        do itensor=1,ENUM_NUM_TENSOR_TYPE
         local_tensor_data(itensor)=MACFLUX_CC(D_DECL(i,j,k),itensor)
        enddo

        do ii=1,3
        do jj=1,3
         center_flux(ii,jj)=zero
        enddo
        enddo

        do itensor=1,ENUM_NUM_TENSOR_TYPE
         call stress_index(itensor,ii,jj)
         center_flux(ii,jj)=local_tensor_data(itensor)
        enddo
        center_flux(2,1)=center_flux(1,2)
        center_flux(3,1)=center_flux(1,3)
        center_flux(3,2)=center_flux(2,3)

! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
!         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
!         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |
        center_hoop_22=center_flux(3,3)

         ! traverse all the flux face centroids associated with the
         ! "force_dir" MAC grid control volume (i,j,k)
        do dir_flux=0,SDIM-1
        do side_flux=0,1

         iflux_array(1)=i
         iflux_array(2)=j
         if (SDIM.eq.3) then
          iflux_array(SDIM)=k
         endif

         do dir_local=1,SDIM
          box_type_flux(dir_local)=0
         enddo
         box_type_flux(dir_flux+1)=1

         if (side_flux.eq.0) then
          ! do nothing
         else if (side_flux.eq.1) then
          iflux_array(dir_flux+1)=iflux_array(dir_flux+1)+1
         else
          print *,"side_flux invalid"
          stop
         endif
         iflux=iflux_array(1)
         jflux=iflux_array(2)
         if (SDIM.eq.3) then
          kflux=iflux_array(SDIM)
         else
          kflux=0
         endif

         call box_type_to_grid_type(grid_type_flux,box_type_flux)

          ! grid_type_flux=0,1, or 2
         call gridstenMAC_level(xsten_flux,iflux,jflux,kflux, &
                 level,nhalf,grid_type_flux)
         do dir_local=1,SDIM
          x_at_flux_point(side_flux+1,dir_flux+1,dir_local)= &
                  xsten_flux(0,dir_local)
          x_at_flux_point_local(dir_local)= &
                  xsten_flux(0,dir_local)
         enddo ! dir_local=1..sdim
       
         do itensor=1,ENUM_NUM_TENSOR_TYPE
          if (grid_type_flux.eq.0) then
           grid_type_sanity=grid_type_X
           local_tensor_data(itensor)= &
             MACFLUX_X(D_DECL(iflux,jflux,kflux),itensor) 
          else if (grid_type_flux.eq.1) then
           grid_type_sanity=grid_type_Y
           local_tensor_data(itensor)= &
             MACFLUX_Y(D_DECL(iflux,jflux,kflux),itensor) 
          else if ((grid_type_flux.eq.SDIM-1).and.(SDIM.eq.3)) then
           grid_type_sanity=grid_type_Z
           local_tensor_data(itensor)= &
             MACFLUX_Z(D_DECL(iflux,jflux,kflux),itensor) 
          else
           print *,"grid_type_flux invalid"
           stop
          endif
         enddo ! itensor=1..ENUM_NUM_TENSOR_TYPE

         if (grid_type_flux.eq.grid_type_sanity) then
          ! do nothing
         else
          print *,"grid_type_sanity failed"
          stop
         endif

         do ii=1,3
         do jj=1,3
          DISP_TEN(ii,jj)=zero
         enddo
         enddo

! grad u=| u_r  u_t/r-v/r  u_z  |
!        | v_r  v_t/r+u/r  v_z  |
!        | w_r  w_t/r      w_z  |
! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
!         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
!         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |
         do itensor=1,ENUM_NUM_TENSOR_TYPE
          call stress_index(itensor,ii,jj)
          DISP_TEN(ii,jj)=local_tensor_data(itensor)
         enddo
         DISP_TEN(2,1)=DISP_TEN(1,2)
         DISP_TEN(3,1)=DISP_TEN(1,3)
         DISP_TEN(3,2)=DISP_TEN(2,3)

         do dir_row=1,3
          do dir_column=1,3
           if (dir_flux.eq.0) then
            xflux_local(side_flux,dir_row,dir_column)= &
              DISP_TEN(dir_row,dir_column)
           else if (dir_flux.eq.1) then
            yflux_local(side_flux,dir_row,dir_column)= &
              DISP_TEN(dir_row,dir_column)
           else if ((dir_flux.eq.2).and.(SDIM.eq.3)) then
            zflux_local(side_flux,dir_row,dir_column)= &
              DISP_TEN(dir_row,dir_column)
           else
            print *,"dir_flux invalid"
            stop
           endif
          enddo ! dir_column
         enddo ! dir_row

         ii=0
         jj=0
         kk=0
         if (dir_flux.eq.0) then
          ii=1
         else if (dir_flux.eq.1) then
          jj=1
         else if ((dir_flux.eq.2).and.(SDIM.eq.3)) then
          kk=1
         else
          print *,"dir_flux invalid"
          stop
         endif
         do im_LS=1,num_materials
          if (side_flux.eq.0) then
           LS_outside(im_LS)=levelpc(D_DECL(i-ii,j-jj,k-kk),im_LS)
          else if (side_flux.eq.1) then
           LS_outside(im_LS)=levelpc(D_DECL(i+ii,j+jj,k+kk),im_LS)
          else
           print *,"side_flux invalid"
           stop
          endif
         enddo !im_LS=1...num_materials

         call get_primary_material(LS_outside,local_mask)
         if (local_mask.eq.im_elastic_p1) then
          local_mask=1
         else if ((local_mask.ge.1).and. &
                  (local_mask.le.num_materials)) then
          local_mask=0
         else
          print *,"local_mask invalid"
          stop
         endif
         ! LS>0 if clamped
         call SUB_clamped_LS(x_at_flux_point_local,cur_time,LS_clamped, &
          vel_clamped,temperature_clamped,prescribed_flag,dx)
         if (LS_clamped.ge.zero) then
          local_mask=0
         else if (LS_clamped.lt.zero) then
          ! do nothing
         else
          print *,"LS_clamped invalid"
          stop
         endif

         mask_flux_point(side_flux+1,dir_flux+1)=local_mask

        enddo ! side_flux=0..1
        enddo ! dir_flux=0..sdim-1

         ! divergence of fluxes goes here.

         ! [n dot tau dot n] = - sigma kappa
         ! [n dot tau dot tj] = 0
    
         !x_at_flux_point(side_flux+1,dir_flux+1,dir_local)

        dir_local=1
        hx=x_at_flux_point(2,dir_local,dir_local)- &
           x_at_flux_point(1,dir_local,dir_local)
        dir_local=2
        hy=x_at_flux_point(2,dir_local,dir_local)- &
           x_at_flux_point(1,dir_local,dir_local)
        dir_local=SDIM
        hz=x_at_flux_point(2,dir_local,dir_local)- &
           x_at_flux_point(1,dir_local,dir_local)

        if ((hx.gt.zero).and.(hy.gt.zero).and.(hz.gt.zero)) then

         if (rzflag.eq.COORDSYS_CARTESIAN) then
          !areax=dz dy
          !areay=dx dz
          !areaz=dx dy
          !volume=dx dy dz
          rplus=one
          rminus=one
         else if ((rzflag.eq.COORDSYS_RZ).or. &
                  (rzflag.eq.COORDSYS_CYLINDRICAL)) then
          ! areax=dz * dtheta * r
          ! areaz=(dtheta/2)*(r_{2}^2 - r_{1}^{2})=dtheta*dr*r
          ! areay=dr * dz
          ! volume=dtheta*dz*dr r
          rplus=x_at_flux_point(2,1,1)
          if (i.eq.0) then
           rminus=zero
          else if (i.gt.0) then
           rminus=x_at_flux_point(1,1,1)
          else
           print *,"i invalid"
           stop
          endif
          if (rzflag.eq.COORDSYS_CYLINDRICAL) then
           hy=hy*rval
          endif
         else
          print *,"rzflag invalid"
          stop
         endif

         do dir_row=1,SDIM

          force(dir_row)=zero

          dir_local=1

          if ((rplus.gt.zero).and. &
              (rminus.ge.zero).and. &
              (rval.gt.zero).and. &
              (hx.gt.zero)) then
           force(dir_row)=force(dir_row)+ &
            (rplus*mask_flux_point(2,dir_local)* &
                   xflux_local(1,dir_row,dir_local)- &
             rminus*mask_flux_point(1,dir_local)* &
                    xflux_local(0,dir_row,dir_local))/ &
            (rval*hx)
          else
           print *,"rplus=",rplus
           print *,"rminus=",rminus
           print *,"rval=",rval
           print *,"hx=",hx
           print *,"rplus, rminus, rval, or hx invalid"
           stop
          endif

          if (hy.gt.zero) then
           dir_local=2
           force(dir_row)=force(dir_row)+ &
            (mask_flux_point(2,dir_local)*yflux_local(1,dir_row,dir_local)- &
             mask_flux_point(1,dir_local)*yflux_local(0,dir_row,dir_local))/hy
          else
           print *,"hy=",hy
           print *,"hy invalid"
           stop
          endif

          if (SDIM.eq.3) then
           if (hz.gt.zero) then
            dir_local=SDIM
            force(dir_row)=force(dir_row)+ &
             (mask_flux_point(2,dir_local)*zflux_local(1,dir_row,dir_local)- &
              mask_flux_point(1,dir_local)*zflux_local(0,dir_row,dir_local))/hz
           else
            print *,"hz=",hz
            print *,"hz invalid"
            stop
           endif
          endif

         enddo ! dir_row=1..sdim
                   
         if (rzflag.eq.COORDSYS_CARTESIAN) then
          bodyforce=zero
         else if (rzflag.eq.COORDSYS_RZ) then

          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
! div S = | (r S_11)_r/r + (S_12)_t/r - S_22/r  + (S_13)_z |
!         | (r S_21)_r/r + (S_22)_t/r + S_12/r  + (S_23)_z |
!         | (r S_31)_r/r + (S_32)_t/r +           (S_33)_z |
           ! -T33/r
          dir_row=1
          if (rval.gt.zero) then
           bodyforce=-center_hoop_22/rval
          else
           print *,"rval invalid"
           stop
          endif

          if (abs(bodyforce).lt.OVERFLOW_CUTOFF) then
           ! do nothing
          else
           print *,"bodyforce overflow bodyforce,rval:",bodyforce,rval
           stop
          endif
          force(dir_row)=force(dir_row)+bodyforce

         else if (rzflag.eq.COORDSYS_CYLINDRICAL) then
          ! -T22/r
          dir_row=1

          if (rval.gt.zero) then
           bodyforce=-center_flux(2,2)/rval
          else
           print *,"rval invalid"
           stop
          endif

          force(dir_row)=force(dir_row)+bodyforce
         else
          print *,"rzflag invalid"
          stop
         endif 

         if (rzflag.eq.COORDSYS_CARTESIAN) then
          bodyforce=zero
         else if (rzflag.eq.COORDSYS_RZ) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          ! do nothing
         else if (rzflag.eq.COORDSYS_CYLINDRICAL) then ! T12/r
          dir_row=2

          if (rval.gt.zero) then
           bodyforce=center_flux(1,2)/rval
          else
           print *,"rval invalid"
           stop
          endif

          force(dir_row)=force(dir_row)+bodyforce
         else
          print *,"rzflag invalid"
          stop
         endif 

         if (is_rigid(im_elastic_p1).eq.1) then
          print *,"im_elastic should not be an is_rigid material"
          stop
         else if (is_rigid(im_elastic_p1).eq.0) then
          do dir_row=1,SDIM
           force(dir_row)=force(dir_row)*dt
          enddo
         else
          print *,"is_rigid invalid"
          stop
         endif

         do dir_row=1,SDIM
          deninv=rhoinvfab(D_DECL(i,j,k))

          if (deninv.ge.zero) then 
           if (abs(force(dir_row)).lt.OVERFLOW_CUTOFF) then
            ! do nothing
           else
            print *,"elastic overflow dir_row,force ",dir_row,force(dir_row)
            print *,"i,j,k,deninv ",i,j,k,deninv
            stop
           endif
 
           if (dir_row.eq.force_dir+1) then
            XFORCE_local=force(dir_row)*deninv
            SNEW(D_DECL(i,j,k))=SNEW(D_DECL(i,j,k))+XFORCE_local
           endif

          else
           print *,"deninv invalid"
           stop
          endif
         enddo ! dir_row = 1 ..sdim

        else
         print *,"hx,hy, or hz invalid"
         stop
        endif

       else if (mask_control_volume.eq.0) then
        ! do nothing
       else
        print *,"mask_control_volume invalid"
        stop
       endif

      enddo ! k
      enddo ! j
      enddo ! i
      
      return 
      end subroutine fort_elastic_force

        ! NavierStokes::veldiffuseALL() (NavierStokes3.cpp)
        ! NavierStokes::assimilate_state_data() (NavierStokes.cpp)
        !     NavierStokes::init_boundary() 
        !     fort_assimilate_statedata(...)
        ! put ns.wall_slip_weight=0.5 for example in the inputs file.
        ! ns.wall_slip_weight=0.0 (default) => do not strengthen the slip BC
        ! ns.wall_slip_weight=1.0 => strongest imposition of slip BC
      subroutine fort_assimilate_statedata( &
       isweep, &
       law_of_the_wall, &
       wall_slip_weight, &
       static_damping_coefficient, &
       im_solid_map, &
       level, &
       finest_level, &
       nstate, &
       nparts, &
       nparts_ghost, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx, &
       dt, &
       cur_time, & ! cur_time
       LS_state,DIMS(LS_state), & ! state data
       state,DIMS(state), &       ! state data
       macx,DIMS(macx), &
       macy,DIMS(macy), &
       macz,DIMS(macz), &
       ughostx,DIMS(ughostx), &  ! stores the slip velocity
       ughosty,DIMS(ughosty), &  
       ughostz,DIMS(ughostz)) &
      bind(c,name='fort_assimilate_statedata')
      use probf90_module
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: isweep
      real(amrex_real), INTENT(in) :: wall_slip_weight
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: law_of_the_wall(num_materials)
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_ghost
      real(amrex_real), INTENT(in) :: static_damping_coefficient(num_materials)
      integer, INTENT(in), target :: im_solid_map(nparts_ghost)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in), target :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in), target :: xlo(SDIM)
      real(amrex_real), INTENT(in), target :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: cur_time
      integer, INTENT(in) :: DIMDEC(LS_state)
      integer, INTENT(in) :: DIMDEC(state)
      integer, INTENT(in) :: DIMDEC(macx)
      integer, INTENT(in) :: DIMDEC(macy)
      integer, INTENT(in) :: DIMDEC(macz)
      integer, INTENT(in) :: DIMDEC(ughostx)
      integer, INTENT(in) :: DIMDEC(ughosty)
      integer, INTENT(in) :: DIMDEC(ughostz)

      real(amrex_real), INTENT(inout), target :: LS_state(DIMV(LS_state), &
           num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_state_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: state(DIMV(state),nstate)
      real(amrex_real), pointer :: state_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: macx(DIMV(macx))
      real(amrex_real), pointer :: macx_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: macy(DIMV(macy))
      real(amrex_real), pointer :: macy_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: macz(DIMV(macz))
      real(amrex_real), pointer :: macz_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: ughostx(DIMV(ughostx),nparts_ghost*SDIM) 
      real(amrex_real), pointer :: ughostx_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: ughosty(DIMV(ughosty),nparts_ghost*SDIM) 
      real(amrex_real), pointer :: ughosty_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: ughostz(DIMV(ughostz),nparts_ghost*SDIM) 
      real(amrex_real), pointer :: ughostz_ptr(D_DECL(:,:,:),:)
      integer i,j,k
      integer iface,jface,kface
      integer icell,jcell,kcell
      integer ileft,jleft,kleft
      integer iright,jright,kright
      integer i_nbr,j_nbr,k_nbr
      integer ii,jj,kk
      integer iii,jjj,kkk
      integer im
      integer im_primary
      integer im_stencil
      integer im_stencil_left
      integer im_stencil_right
      real(amrex_real) :: local_damping
      real(amrex_real) :: LS_local(num_materials)
      real(amrex_real) :: LS_stencil(num_materials)
      real(amrex_real) :: LS_LEFT(num_materials)
      real(amrex_real) :: LS_RIGHT(num_materials)
      integer, PARAMETER :: nhalf=3
      real(amrex_real), target :: xsten(-nhalf:nhalf,SDIM)
      type(assimilate_parm_type) :: assimilate_parm
      type(assimilate_out_parm_type) :: assimilate_out_parm
      integer cell_flag
      integer veldir
      integer dir
      integer dirtan
      integer side
      integer side_nbr
      integer partid
      real(amrex_real) velsum(SDIM)
      real(amrex_real) wtsum
      real(amrex_real) velface

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in ratemasschange"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid in GODUNOV_3D.F90 "
       print *,"nstate=",nstate
       stop
      endif
      if ((nparts.ge.0).and.(nparts.le.num_materials)) then 
       ! do nothing
      else
       print *,"nparts invalid fort_assimilate_statedata"
       stop
      endif
      if ((nparts_ghost.ge.1).and. &
          (nparts_ghost.le.num_materials).and. &
          (nparts_ghost.ge.nparts)) then 
       ! do nothing
      else
       print *,"nparts_ghost invalid fort_assimilate_statedata"
       stop
      endif

      if ((nparts_ghost.eq.nparts).or.(nparts_ghost.eq.1)) then
       ! do nothing
      else
       print *,"nparts_ghost invalid"
       stop
      endif

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif 
      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid"
       stop
      endif 
      do im=1,num_materials
       if ((law_of_the_wall(im).eq.0).or. &
           (law_of_the_wall(im).eq.1).or. &
           (law_of_the_wall(im).eq.2)) then
        ! do nothing
       else
        print *,"law_of_the_wall invalid"
        stop
       endif
      enddo ! im=1..num_materials

      LS_state_ptr=>LS_state
      call checkbound_array(fablo,fabhi,LS_state_ptr,1,-1)
      state_ptr=>state
      call checkbound_array(fablo,fabhi,state_ptr,1,-1)
      macx_ptr=>macx
      call checkbound_array1(fablo,fabhi,macx_ptr,0,0)
      macy_ptr=>macy
      call checkbound_array1(fablo,fabhi,macy_ptr,0,1)
      macz_ptr=>macz
      call checkbound_array1(fablo,fabhi,macz_ptr,0,SDIM-1)

      ughostx_ptr=>ughostx
      call checkbound_array(fablo,fabhi,ughostx_ptr,0,0)
      ughosty_ptr=>ughosty
      call checkbound_array(fablo,fabhi,ughosty_ptr,0,1)
      ughostz_ptr=>ughostz
      call checkbound_array(fablo,fabhi,ughostz_ptr,0,SDIM-1)

      assimilate_parm%cur_time=cur_time
      assimilate_parm%dt=dt
      assimilate_parm%nhalf=nhalf
      assimilate_parm%nstate=nstate
      assimilate_parm%nparts=nparts
      assimilate_parm%nparts_ghost=nparts_ghost
      assimilate_parm%im_solid_map=>im_solid_map
      assimilate_parm%level=level
      assimilate_parm%finest_level=finest_level
      assimilate_parm%bfact=bfact
      assimilate_parm%dx=>dx
      assimilate_parm%xlo=>xlo
      assimilate_parm%fablo=>fablo
      assimilate_parm%fabhi=>fabhi
      assimilate_parm%ughostx=>ughostx
      assimilate_parm%ughosty=>ughosty
      assimilate_parm%ughostz=>ughostz

      assimilate_parm%dxmin=dx(1)
      if (dx(2).lt.assimilate_parm%dxmin) then
       assimilate_parm%dxmin=dx(2)
      endif
      if (dx(SDIM).lt.assimilate_parm%dxmin) then
       assimilate_parm%dxmin=dx(SDIM)
      endif
       
      assimilate_out_parm%LS_state=>LS_state
      assimilate_out_parm%state=>state
      assimilate_out_parm%macx=>macx
      assimilate_out_parm%macy=>macy
      assimilate_out_parm%macz=>macz

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      cell_flag=-1

       ! wall_slip_weight=zero is the default.
      if ((wall_slip_weight.ge.zero).and. &
          (wall_slip_weight.le.one)) then
       ! do nothing
      else
       print *,"wall_slip_weight invalid"
       stop
      endif

      if (isweep.eq.0) then

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        call gridsten_level(xsten,i,j,k,level,nhalf)

        assimilate_parm%xsten=>xsten
        if (is_in_probtype_list().eq.1) then
         call SUB_ASSIMILATE(assimilate_parm,assimilate_out_parm, &
          i,j,k,cell_flag)
        else
         ! do nothing
        endif

         ! check if a fluid cell neighbors a solid cell
        do im=1,num_materials
         LS_local(im)=LS_state(D_DECL(i,j,k),im)
        enddo
        ! first checks the rigid materials for a positive LS; if none
        ! exist, then the subroutine checks the fluid materials.
        call get_primary_material(LS_local,im_primary) 

        if (is_rigid(im_primary).eq.0) then
         ! check all neighbors in "star stencil" for solid cells.
         ! for each solid cell, update the present cell center velocity
         ! with a weighted average of the slip velocity and the present
         ! velocity.
         do veldir=1,SDIM
          velsum(veldir)=zero
         enddo
         wtsum=zero

         do dir=1,SDIM
          ii=0
          jj=0
          kk=0
          if (dir.eq.1) then
           ii=1
          else if (dir.eq.2) then
           jj=1
          else if ((dir.eq.3).and.(SDIM.eq.3)) then
           kk=1
          else
           print *,"dir invalid"
           stop
          endif
          do side=1,2
           if (side.eq.1) then
            iface=i
            jface=j
            kface=k
            icell=i-ii
            jcell=j-jj
            kcell=k-kk
           else if (side.eq.2) then
            iface=i+ii
            jface=j+jj
            kface=k+kk
            icell=i+ii
            jcell=j+jj
            kcell=k+kk
           else
            print *,"side invalid"
            stop
           endif

           do im=1,num_materials
            LS_stencil(im)=LS_state(D_DECL(icell,jcell,kcell),im)
           enddo
           call get_primary_material(LS_stencil,im_stencil) 
           if (is_rigid(im_stencil).eq.0) then
            ! do nothing
           else if (is_rigid(im_stencil).eq.1) then
            partid=1
            do while ((im_solid_map(partid)+1.ne.im_stencil).and. &
                      (partid.lt.nparts_ghost))
             partid=partid+1
            enddo
            if (im_solid_map(partid)+1.eq.im_stencil) then
             do veldir=1,SDIM
              if (dir.eq.1) then
               velface=ughostx(D_DECL(iface,jface,kface), &
                       (partid-1)*SDIM+veldir)
              else if (dir.eq.2) then
               velface=ughosty(D_DECL(iface,jface,kface), &
                       (partid-1)*SDIM+veldir)
              else if ((dir.eq.3).and.(SDIM.eq.3)) then
               velface=ughostz(D_DECL(iface,jface,kface), &
                       (partid-1)*SDIM+veldir)
              else
               print *,"dir invalid"
               stop
              endif

              velsum(veldir)=velsum(veldir)+velface
             enddo ! veldir=1..sdim
             wtsum=wtsum+one
            else
             print *,"im_solid_map(partid) invalid"
             stop
            endif
           else
            print *,"is_rigid(im_stencil) invalid"
            stop
           endif 
          enddo ! side=1,2
         enddo ! dir=1..SDIM
         if (wtsum.eq.zero) then
          ! do nothing
         else if (wtsum.gt.zero) then
          do veldir=1,SDIM
           state(D_DECL(i,j,k),veldir)= &
                  (one-wall_slip_weight)*state(D_DECL(i,j,k),veldir)+ &
                  wall_slip_weight*velsum(veldir)/wtsum
          enddo
         else
          print *,"wtsum invalid"
          stop
         endif 

        else if (is_rigid(im_primary).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(im_primary) invalid"
         stop
        endif

         ! static_damping_coefficient=zero is the default.
         ! cell centered velocity
        if (static_damping_coefficient(im_primary).gt.zero) then
         ! v_t = -mu v =>  v^{n+1} - v^{n} = -mu dt v^{n+1}
         ! v^{n+1}=v^{n}/(1+mu dt)
         do veldir=1,SDIM
          state(D_DECL(i,j,k),veldir)=state(D_DECL(i,j,k),veldir)/ &
                  (one+static_damping_coefficient(im_primary)*dt)
         enddo
        else if (static_damping_coefficient(im_primary).eq.zero) then
         ! do nothing
        else
         print *,"static_damping_coefficient(im_primary) invalid"
         stop
        endif
       enddo ! k
       enddo ! j
       enddo ! i

      else if (isweep.eq.1) then

       do cell_flag=0,SDIM-1

        ii=0
        jj=0
        kk=0
        if (cell_flag.eq.0) then
         ii=1
        else if (cell_flag.eq.1) then
         jj=1
        else if ((cell_flag.eq.2).and.(SDIM.eq.3)) then
         kk=1
        else
         print *,"cell_flag invalid"
         stop
        endif
        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
             growlo,growhi,0,cell_flag) 

        do k=growlo(3),growhi(3)
        do j=growlo(2),growhi(2)
        do i=growlo(1),growhi(1)

          !cell_flag=0..sdim-1
         call gridstenMAC_level(xsten,i,j,k,level,nhalf,cell_flag)

         assimilate_parm%xsten=>xsten
         if (is_in_probtype_list().eq.1) then
          call SUB_ASSIMILATE(assimilate_parm,assimilate_out_parm, &
           i,j,k,cell_flag)
         else
          ! do nothing
         endif

          ! check that both adjoining cells are fluid cells
         ileft=i-ii
         jleft=j-jj
         kleft=k-kk
         iright=i
         jright=j
         kright=k
         do im=1,num_materials
          LS_LEFT(im)=LS_state(D_DECL(ileft,jleft,kleft),im)
         enddo
         call get_primary_material(LS_LEFT,im_stencil_left)
         do im=1,num_materials
          LS_RIGHT(im)=LS_state(D_DECL(iright,jright,kright),im)
         enddo
         call get_primary_material(LS_RIGHT,im_stencil_right)

         if (is_rigid(im_stencil_left).eq.0) then
          if (is_rigid(im_stencil_right).eq.0) then
           ! check if any neighbor solid faces oriented perpendicular to the
           ! given face.
           !   side_nbr=1   side_nbr=2
           !     --------------------
           !     |    s?   |  s?    |
           !     |         |        |
           !     --------------------
           !     |    f   f|f  f    |
           !     |        f|f       |
           !     --------------------
           !     |    s?   |  s?    |
           !     |         |        |
           !     --------------------
           ! loop through the "s?" cells, wtsum will contain the count
           wtsum=zero
           do dirtan=1,SDIM
            if (dirtan.ne.cell_flag+1) then
             do side_nbr=1,2
              if (side_nbr.eq.1) then
               i_nbr=ileft
               j_nbr=jleft
               k_nbr=kleft
              else if (side_nbr.eq.2) then
               i_nbr=iright
               j_nbr=jright
               k_nbr=kright
              else
               print *,"side_nbr invalid"
               stop
              endif
              iii=0
              jjj=0
              kkk=0
              if (dirtan.eq.1) then
               iii=1
              else if (dirtan.eq.2) then
               jjj=1
              else if ((dirtan.eq.3).and.(SDIM.eq.3)) then
               kkk=1
              else
               print *,"dirtan invalid"
               stop
              endif
            
              do side=1,2
               if (side.eq.1) then
                icell=i_nbr-iii
                jcell=j_nbr-jjj
                kcell=k_nbr-kkk
               else if (side.eq.2) then
                icell=i_nbr+iii
                jcell=j_nbr+jjj
                kcell=k_nbr+kkk
               else
                print *,"side invalid"
                stop
               endif

               do im=1,num_materials
                LS_stencil(im)=LS_state(D_DECL(icell,jcell,kcell),im)
               enddo
               call get_primary_material(LS_stencil,im_stencil) 
               if (is_rigid(im_stencil).eq.0) then
                ! do nothing
               else if (is_rigid(im_stencil).eq.1) then
                partid=1
                do while ((im_solid_map(partid)+1.ne.im_stencil).and. &
                          (partid.lt.nparts_ghost))
                 partid=partid+1
                enddo
                if (im_solid_map(partid)+1.eq.im_stencil) then
                 wtsum=wtsum+one
                else
                 print *,"im_solid_map(partid) invalid"
                 stop
                endif
               else
                print *,"is_rigid(im_stencil) invalid"
                stop
               endif 
              enddo ! side=1,2
             enddo ! side_nbr=1..2
            else if (dirtan.eq.cell_flag+1) then
             ! do nothing
            else
             print *,"dirtan invalid"
             stop
            endif
           enddo ! dirtan=1..SDIM
           if (wtsum.eq.zero) then
            ! do nothing
           else if (wtsum.gt.zero) then
            do veldir=1,SDIM

             velsum(veldir)= &
                half*(state(D_DECL(ileft,jleft,kleft),veldir)+ &
                      state(D_DECL(iright,jright,kright),veldir))

              !wall_slip_weight=0.0 is the default.
             if (cell_flag+1.eq.veldir) then
              if (veldir.eq.1) then
               macx(D_DECL(i,j,k))= &
                   (one-wall_slip_weight)*macx(D_DECL(i,j,k))+ &
                   wall_slip_weight*velsum(veldir)
              else if (veldir.eq.2) then
               macy(D_DECL(i,j,k))= &
                   (one-wall_slip_weight)*macy(D_DECL(i,j,k))+ &
                   wall_slip_weight*velsum(veldir)
              else if ((veldir.eq.3).and.(SDIM.eq.3)) then
               macz(D_DECL(i,j,k))= &
                   (one-wall_slip_weight)*macz(D_DECL(i,j,k))+ &
                   wall_slip_weight*velsum(veldir)
              else
               print *,"veldir invalid"
               stop
              endif
             else if ((cell_flag+1.ge.1).and. &
                      (cell_flag+1.le.SDIM)) then
              ! do nothing
             else
              print *,"cell_flag invalid"
              stop
             endif
            enddo ! veldir=1..sdim
           else
            print *,"wtsum invalid"
            stop
           endif 
          else if (is_rigid(im_stencil_right).eq.1) then
           ! do nothing
          else
           print *,"is_rigid(im_stencil_right) (right) invalid"
           stop
          endif
         else if (is_rigid(im_stencil_left).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(im_stencil_left) (left) invalid"
          stop
         endif
         do im=1,num_materials
          LS_stencil(im)=half*(LS_LEFT(im)+LS_RIGHT(im))
         enddo
         call get_primary_material(LS_stencil,im_stencil)

         ! static_damping_coefficient=zero is the default.
         local_damping=static_damping_coefficient(im_stencil)
         if (local_damping.gt.zero) then
          if (cell_flag.eq.0) then
           macx(D_DECL(i,j,k))=macx(D_DECL(i,j,k))/(one+dt*local_damping)
          else if (cell_flag.eq.1) then
           macy(D_DECL(i,j,k))=macy(D_DECL(i,j,k))/(one+dt*local_damping)
          else if ((cell_flag.eq.2).and.(SDIM.eq.3)) then
           macz(D_DECL(i,j,k))=macz(D_DECL(i,j,k))/(one+dt*local_damping)
          else
           print *,"cell_flag invalid"
           stop
          endif
         else if (local_damping.eq.zero) then
          ! do nothing
         else
          print *,"local_damping invalid"
          stop
         endif

        enddo ! k
        enddo ! j
        enddo ! i

       enddo ! cell_flag=0...sdim-1

      else
       print *,"isweep invalid"
       stop
      endif


      return
      end subroutine fort_assimilate_statedata


      ! enable_spectral:
      ! 0 - low order
      ! 1 - space/time spectral
      subroutine fort_build_masksem( &
       dx, &
       spectral_cells_level, &
       mask_sweep, &
       level, &
       finest_level, &
       cur_time, &
       enable_spectral, &
       domlo,domhi, &
       vofbc, &
       maskcov,DIMS(maskcov), &
       masknbr,DIMS(masknbr), &
       mask,DIMS(mask), &
       oldmask,DIMS(oldmask), &
       vfrac,DIMS(vfrac), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       bfact_fine) &
      bind(c,name='fort_build_masksem')

      use probcommon_module
      use global_utility_module
      use MOF_routines_module
      use probf90_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(inout) :: spectral_cells_level(num_materials)
      integer, INTENT(in) :: mask_sweep
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: cur_time
      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: bfact_fine
      integer, INTENT(in) :: vofbc(SDIM,2)
      integer, INTENT(in) :: DIMDEC(maskcov) 
      integer, INTENT(in) :: DIMDEC(masknbr) 
      integer, INTENT(in) :: DIMDEC(mask) 
      integer, INTENT(in) :: DIMDEC(oldmask) 
      integer, INTENT(in) :: DIMDEC(vfrac) 
      real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: masknbr(DIMV(masknbr))
      real(amrex_real), pointer :: masknbr_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out), target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: oldmask(DIMV(oldmask))
      real(amrex_real), pointer :: oldmask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: vfrac(DIMV(vfrac),num_materials)
      real(amrex_real), pointer :: vfrac_ptr(D_DECL(:,:,:),:)

      integer im,imcrit,im_max
      integer sumtag
      integer tag(num_materials)
      integer local_maskSEM
      integer old_maskSEM
      integer i,j,k
      integer dir
      integer iofs,jofs,kofs,kofs_hi
      integer stripstat
      integer test_maskcov
      integer test_maskSEM
      integer covered_count
      integer uncovered_count
      real(amrex_real) vfracmax
      real(amrex_real) vfractest(num_materials)
      real(amrex_real) vfrac_fluid
      real(amrex_real) vfrac_sum_fluid,vfrac_sum_solid
      integer clamped_cell_in_element
      real(amrex_real) xclamped(SDIM)
      integer iregions
      real(amrex_real) LS_clamped,charfn
      real(amrex_real) vel_clamped(SDIM)
      real(amrex_real) temperature_clamped
      integer prescribed_flag
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid build masksem"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid43"
       stop
      endif
      if ((bfact_fine.gt.bfact).or.(bfact_fine.lt.1)) then
       print *,"bfact_fine invalid"
       stop
      endif

      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.1)) then
       print *,"enable_spectral invalid build masksem"
       stop
      endif

      maskcov_ptr=>maskcov
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      masknbr_ptr=>masknbr
      call checkbound_array1(fablo,fabhi,masknbr_ptr,1,-1)
      mask_ptr=>mask
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      oldmask_ptr=>oldmask
      call checkbound_array1(fablo,fabhi,oldmask_ptr,1,-1)
      vfrac_ptr=>vfrac
      call checkbound_array(fablo,fabhi,vfrac_ptr,0,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call strip_status(i,j,k,bfact,stripstat)

       if (stripstat.eq.1) then

        do im=1,num_materials
         tag(im)=0
        enddo
        imcrit=0

        covered_count=0
        uncovered_count=0

        im_max=0
        vfracmax=zero

        if (SDIM.eq.2) then
         kofs_hi=0
        else if (SDIM.eq.3) then
         kofs_hi=bfact-1
        else
         print *,"dimension bust"
         stop
        endif

        clamped_cell_in_element=0

        do kofs=0,kofs_hi
        do jofs=0,bfact-1
        do iofs=0,bfact-1

         call gridsten_level(xsten,i+iofs,j+jofs,k+kofs,level,nhalf)
         do dir=1,SDIM
          xclamped(dir)=xsten(0,dir)
         enddo
          ! LS>0 if clamped
         call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
           vel_clamped,temperature_clamped,prescribed_flag,dx)

         if (LS_clamped.ge.zero) then
          clamped_cell_in_element=1
         else if (LS_clamped.lt.zero) then
          ! do nothing
         else
          print *,"LS_clamped is NaN"
          stop
         endif

         do iregions=1,number_of_source_regions
          call SUB_CHARFN_REGION(iregions,xclamped,cur_time,charfn)
          if (charfn.eq.one) then
           clamped_cell_in_element=1
          else if (charfn.eq.zero) then
           ! do nothing
          else
           print *,"charfn invalid"
           stop
          endif
         enddo !iregions=1,number_of_source_regions

         vfrac_sum_fluid=zero
         vfrac_sum_solid=zero

         do im=1,num_materials
          vfractest(im)=vfrac(D_DECL(i+iofs,j+jofs,k+kofs),im)
          if (is_rigid(im).eq.1) then
           vfrac_sum_solid=vfrac_sum_solid+vfractest(im)

           if (vfractest(im).gt.vfracmax) then
            im_max=im
            vfracmax=vfractest(im)
           else if (vfractest(im).le.vfracmax) then
            ! do nothing
           else
            print *,"vfractest(im) is NaN vfractest,vfracmax=", &
              vfractest(im),vfracmax
            stop
           endif
           if (vfractest(im).gt.VOFTOL) then
            imcrit=im
            tag(im)=1
           else if (vfractest(im).le.VOFTOL) then
            ! do nothing
           else
            print *,"vfractest(im) is NaN vfractest,VOFTOL=", &
              vfractest(im),VOFTOL
            stop
           endif

          else if (is_rigid(im).eq.0) then

           vfrac_sum_fluid=vfrac_sum_fluid+vfractest(im)

          else
           print *,"is_rigid(im) invalid"
           stop
          endif
         enddo ! im=1..num_materials

         if ((vfrac_sum_solid.lt.-EPS1).or. &
             (vfrac_sum_solid.gt.one+EPS1)) then
          print *,"vfrac_sum_solid out of range"
          print *,"vfrac_sum_solid=",vfrac_sum_solid
          stop
         else if ((vfrac_sum_solid.ge.-EPS1).and. &
                  (vfrac_sum_solid.le.one+EPS1)) then
          ! do nothing
         else
          print *,"vfrac_sum_solid is NaN: ",vfrac_sum_solid
          stop
         endif

         if (abs(one-vfrac_sum_fluid).gt.0.01) then
          print *,"vfrac_sum_fluid out of range in build mask sem"
          print *,"vfrac_sum_solid=",vfrac_sum_solid
          print *,"vfrac_sum_fluid=",vfrac_sum_fluid
          print *,"i,j,k ",i,j,k
          print *,"iofs,jofs,kofs ",iofs,jofs,kofs
          print *,"level,finest_level ",level,finest_level
          stop
         endif
 
         do im=1,num_materials
          if (is_rigid(im).eq.1) then
           ! do nothing
          else if (is_rigid(im).eq.0) then
           vfrac_fluid=(one-vfrac_sum_solid)*vfractest(im)
           if (vfrac_fluid.gt.vfracmax) then
            im_max=im
            vfracmax=vfrac_fluid
           endif  
           if (vfrac_fluid.gt.VOFTOL) then
            imcrit=im
            tag(im)=1
           endif
          else
           print *,"is_rigid(im) invalid"
           stop
          endif
         enddo ! im=1..num_materials

         test_maskcov=NINT(maskcov(D_DECL(i+iofs,j+jofs,k+kofs))) 
         if (test_maskcov.eq.1) then
          uncovered_count=uncovered_count+1
         else if (test_maskcov.eq.0) then
          covered_count=covered_count+1
         else
          print *,"test_maskcov invalid"
          stop
         endif
           
        enddo
        enddo
        enddo ! iofs,jofs,kofs

        if ((imcrit.lt.1).or.(imcrit.gt.num_materials)) then
         print *,"imcrit invalid"
         stop
        endif
        if ((im_max.lt.1).or.(im_max.gt.num_materials)) then
         print *,"im_max invalid"
         stop
        endif

        if (vfracmax.gt.VOFTOL) then
         ! do nothing
        else
         print *,"vfracmax invalid:",vfracmax
         stop
        endif

        sumtag=0
        do im=1,num_materials
         sumtag=sumtag+tag(im)
        enddo

        if ((sumtag.le.0).or.(sumtag.gt.num_materials)) then
         print *,"sumtag invalid"
         stop
        else if (sumtag.eq.1) then

         if (im_max.eq.imcrit) then
          ! do nothing
         else
          print *,"if sumtag==1, need im_max=imcrit"
          print *,"im_max=",im_max
          print *,"imcrit=",imcrit
          stop
         endif

         if (clamped_cell_in_element.eq.1) then
          local_maskSEM=0
         else if (clamped_cell_in_element.eq.0) then
          if (is_rigid(imcrit).eq.1) then
           local_maskSEM=0
          else if (is_ice(imcrit).eq.1) then
           local_maskSEM=0
          else if (is_FSI_rigid(imcrit).eq.1) then
           local_maskSEM=0
          else if ((FSI_flag(imcrit).eq.FSI_FLUID).or. & 
                   (FSI_flag(imcrit).eq.FSI_FLUID_NODES_INIT)) then
           local_maskSEM=imcrit
          else
           print *,"FSI_flag invalid"
           stop
          endif
         else
          print *,"clamped_cell_in_element invalid"
          stop
         endif

        else if (sumtag.gt.1) then
         local_maskSEM=0
        else
         print *,"sumtag bust"
         stop
        endif

        if (covered_count.gt.0) then
         if (uncovered_count.ne.0) then
          print *,"cannot have an element partially covered"
          print *,"must have bfact_{l+1}>=2 * order_{l}"
          stop
         endif
        else if (uncovered_count.gt.0) then
         if (covered_count.ne.0) then
          print *,"cannot have an element partially covered"
          print *,"must have bfact_{l+1}>=2 * order_{l}"
          stop
         endif
        else
         print *,"covered_count or uncovered_count invalid"
         stop
        endif 

        test_maskcov=NINT(maskcov(D_DECL(i,j,k))) 
        if ((test_maskcov.eq.0).or.(test_maskcov.eq.1)) then
         ! do nothing
        else
         print *,"test_maskcov invalid"
         stop
        endif

        if (mask_sweep.eq.0) then
         ! do nothing
        else if (mask_sweep.eq.1) then ! check neighboring elements

         old_maskSEM=NINT(oldmask(D_DECL(i,j,k))) 

         if (old_maskSEM.ne.local_maskSEM) then
          print *,"old_maskSEM.ne.local_maskSEM"
          stop
         endif

         if (test_maskcov.eq.1) then
          ! do nothing

          ! high order stencil never includes covered values.
         else if (test_maskcov.eq.0) then 
          local_maskSEM=0
         else
          print *,"test_maskcov invalid"
          stop
         endif

         if ((local_maskSEM.ge.1).and.(local_maskSEM.le.num_materials)) then

          if (SDIM.eq.2) then
           kofs_hi=-1
          else if (SDIM.eq.3) then
           kofs_hi=bfact
          else
           print *,"dimension bust"
           stop
          endif

          do kofs=-1,kofs_hi
          do jofs=-1,bfact
          do iofs=-1,bfact
           test_maskSEM=NINT(oldmask(D_DECL(i+iofs,j+jofs,k+kofs))) 
           if (test_maskSEM.ne.local_maskSEM) then
            local_maskSEM=0
           endif
          enddo
          enddo
          enddo ! iofs,jofs,kofs

         else if (local_maskSEM.eq.0) then
          ! do nothing
         else
          print *,"local_maskSEM invalid"
          stop
         endif

        else
         print *,"mask_sweep invalid"
         stop
        endif


        if (SDIM.eq.2) then
         kofs_hi=0
        else if (SDIM.eq.3) then
         kofs_hi=bfact-1
        else
         print *,"dimension bust"
         stop
        endif

        do kofs=0,kofs_hi
        do jofs=0,bfact-1
        do iofs=0,bfact-1

         mask(D_DECL(i+iofs,j+jofs,k+kofs))=local_maskSEM

         if ((local_maskSEM.ge.1).and.(local_maskSEM.le.num_materials)) then
          spectral_cells_level(local_maskSEM)= &
            spectral_cells_level(local_maskSEM)+one
         else if (local_maskSEM.eq.0) then
          ! do nothing
         else
          print *,"local maskSEM invalid"
          stop
         endif
  
        enddo
        enddo
        enddo ! iofs,jofs,kofs

        if (bfact.eq.1) then
         ! do nothing
        else if ((bfact.ge.2).and.(bfact.le.16)) then
         if (mask_sweep.eq.0) then
          ! do nothing
         else if (mask_sweep.eq.1) then
          if ((local_maskSEM.ge.1).and.(local_maskSEM.le.num_materials)) then
           ! do nothing
          else if (local_maskSEM.eq.0) then
           if (level.eq.finest_level) then
            if (test_maskcov.eq.1) then
             ! do nothing
            else 
             print *,"test_maskcov = 1 on the finest level"
             stop
            endif
           else if ((level.ge.0).and.(level.lt.finest_level)) then
            if (test_maskcov.eq.1) then
             ! do nothing
            else if (test_maskcov.eq.0) then
             ! do nothing
            else
             print *,"test_maskcov invalid"
             stop
            endif
           else
            print *,"level invalid"
            stop
           endif
          else
           print *,"local_maskSEM invalid"
           stop
          endif
         else
          print *,"mask_sweep invalid"
          stop
         endif

        else
         print *,"bfact out of range"
         stop
        endif

       else if (stripstat.eq.0) then
        ! do nothing
       else
        print *,"stripstat invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i,j,k (only interior cells)

      return
      end subroutine fort_build_masksem


       ! called from: NavierStokes3.cpp
      subroutine fort_heatadvance( &
       level, &
       finest_level, &
       cur_time, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       nsolve, &
       nstate, &
       xlo,dx, &
       solxfab,DIMS(solxfab), &
       solyfab,DIMS(solyfab), &
       solzfab,DIMS(solzfab), &
       snew,DIMS(snew), &
       lsnew,DIMS(lsnew), &
       du,DIMS(du), &
       tilelo,tilehi, &
       fablo,fabhi,bfact) &
      bind(c,name='fort_heatadvance')

      use probcommon_module
      use global_utility_module
      use probf90_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: cur_time
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def)
      integer, INTENT(in) :: nsolve
      integer, INTENT(in) :: nstate
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: bfact
      integer :: growlo(3),growhi(3)

      integer, INTENT(in) :: DIMDEC(solxfab)
      integer, INTENT(in) :: DIMDEC(solyfab)
      integer, INTENT(in) :: DIMDEC(solzfab)
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(lsnew)
      integer, INTENT(in) :: DIMDEC(du)
      integer :: i,j,k
      integer :: dir
      integer :: im
      integer :: ibase
      integer :: velcomp

      real(amrex_real), INTENT(in),target :: &
              solxfab(DIMV(solxfab),nparts_def*SDIM)
      real(amrex_real), INTENT(in),target :: &
              solyfab(DIMV(solyfab),nparts_def*SDIM)
      real(amrex_real), INTENT(in),target :: &
              solzfab(DIMV(solzfab),nparts_def*SDIM)
      real(amrex_real), pointer :: solxfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: solyfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: solzfab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(inout),target :: snew(DIMV(snew),nstate)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
              lsnew(DIMV(lsnew),num_materials*(SDIM+1))
      real(amrex_real), pointer :: lsnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: du(DIMV(du),nsolve)
      real(amrex_real), pointer :: du_ptr(D_DECL(:,:,:),:)

      real(amrex_real) Tforce,new_temperature,TEMPERATURE
      real(amrex_real) xclamped(SDIM)
      real(amrex_real) LS_clamped
      real(amrex_real) vel_clamped(SDIM)
      real(amrex_real) temperature_clamped
      integer prescribed_flag
      integer, PARAMETER :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid veladvance"
       stop
      endif
      if (finest_level.ne.fort_finest_level) then
       print *,"finest_level invalid veladvance"
       stop
      endif
      if (nsolve.ne.1) then
       print *,"nsolve invalid"
       stop
      endif

      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid FORT_HEATADVANCE"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid FORT_HEATADVANCE"
       stop
      endif

      solxfab_ptr=>solxfab
      solyfab_ptr=>solyfab
      solzfab_ptr=>solzfab
      snew_ptr=>snew
      lsnew_ptr=>lsnew
      du_ptr=>du

      call checkbound_array(fablo,fabhi,solxfab_ptr,0,0)
      call checkbound_array(fablo,fabhi,solyfab_ptr,0,1)
      call checkbound_array(fablo,fabhi,solzfab_ptr,0,SDIM-1)
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,du_ptr,0,-1)
     
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir=1,SDIM
        xclamped(dir)=xsten(0,dir)
       enddo

        ! LS>0 if clamped
       call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
         vel_clamped,temperature_clamped,prescribed_flag,dx)

       do im=1,num_materials
        ibase=STATECOMP_STATES+(im-1)*num_state_material
        TEMPERATURE=snew(D_DECL(i,j,k),ibase+ENUM_TEMPERATUREVAR+1)
        if (TEMPERATURE.gt.zero) then
         ! do nothing
        else
          print *,"HEATADVANCE: temperature must be positive"
          print *,"i,j,k,im ",i,j,k,im
          print *,"TEMPERATURE= ",TEMPERATURE
          stop
        endif
         ! viscous heating term.
        velcomp=1
        Tforce=du(D_DECL(i,j,k),velcomp)
        new_temperature=TEMPERATURE+Tforce
        if (new_temperature.le.zero) then
         new_temperature=TEMPERATURE
        endif
        snew(D_DECL(i,j,k),ibase+ENUM_TEMPERATUREVAR+1)=new_temperature
       enddo ! im = 1..num_materials

      enddo
      enddo
      enddo
  
      return 
      end subroutine fort_heatadvance

      end module godunov_module

