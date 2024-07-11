#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "EXTRAP_COMP.H"
#include "PROB_F.H"
#include "LEVEL_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

#define maxnumline 10

#define GRADHTOL (0.01)
#define MARCO 0

#define NR_MAX 60
#define NS_MAX NR_MAX*NR_MAX


! new random number generator
      module randomNG
        IMPLICIT NONE
        private
        public::rng_t, rng_seed, rng_uniform

        ! Dimension of the state
        integer, parameter :: ns = 4

        ! Default seed vector
        integer, parameter, dimension(ns) :: default_seed &
             = (/ 521288629, 362436069, 16163801, 1131199299 /)

        ! A data type for storing the state of the RNG
        type :: rng_t
           integer, dimension(ns) :: state =default_seed
        end type rng_t

      contains

        ! Seeds the RNG using a single
        ! integer and a default seed
        ! vector.
        subroutine rng_seed(self,seed)
          type(rng_t), INTENT(inout) :: self
          integer, INTENT(in) :: seed
          self%state(1) = seed
          self%state(2:ns) = default_seed(2:ns)
        end subroutine rng_seed


        ! Draws a uniform real number on [0,1].
        function rng_uniform(self) result(u)
          type(rng_t), INTENT(inout) :: self
          real :: u
          integer :: imz

          imz = self%state(1) - self%state(3)

          if (imz < 0) imz = imz + 2147483579

          self%state(1) = self%state(2)
          self%state(2) = self%state(3)
          self%state(3) = imz
          self%state(4) = 69069 * self%state(4) + 101390424             
          imz = imz + self%state(4)
          u = 0.5d0 + 0.23283064d-9 * imz
        end function rng_uniform

      end module randomNG



      module probf90_module
      use amrex_fort_module, only : amrex_real
      use probcommon_module_types
      use probcommon_module

      real(amrex_real) lowerdiag_initdata(1000),upperdiag_initdata(1000)
      real(amrex_real) diag_initdata(1000),soln_initdata(1000)
      real(amrex_real) rhs_initdata(1000)
      real(amrex_real) xlodiss_initdata,xhidiss_initdata
      real(amrex_real) dtdiss_initdata,posdiss_initdata
      integer ndiss_initdata,ispace_initdata

      contains


      subroutine get_mach_number(tessellate, &
        vel,den,vof,mach)
      use MOF_routines_module
      use global_utility_module
      IMPLICIT NONE
 
      integer, INTENT(in) :: tessellate
      real(amrex_real), INTENT(in) :: vel(SDIM+1)
      real(amrex_real), INTENT(in) :: den(num_materials*num_state_material)
      real(amrex_real), INTENT(in) :: vof(num_materials)
      real(amrex_real), INTENT(out) :: mach
      integer dir
      real(amrex_real) UMACH_local
      real(amrex_real) test_vel
      real(amrex_real) USOUND_local
      real(amrex_real) den_local
      real(amrex_real) temp_local
      real(amrex_real) energy_local
      integer im_primary,im
      integer ispec
      real(amrex_real) :: massfrac_parm(num_species_var+1)

      if (tessellate.eq.0) then
       call get_primary_material_VFRAC(vof,im_primary)
      else if ((tessellate.eq.1).or. &
               (tessellate.eq.3)) then
       im_primary=0
       do im=1,num_materials
        if (im_primary.eq.0) then
         im_primary=im
        else if (vof(im).gt.vof(im_primary)) then
         im_primary=im
        endif
       enddo ! im 
      else
       print *,"tessellate invalid63"
       stop
      endif
     
      if ((im_primary.lt.1).or.(im_primary.gt.num_materials)) then
       print *,"im_primary invalid"
       stop
      endif
 
      UMACH_local=zero
      do dir=1,SDIM
       test_vel=vel(dir)
       UMACH_local=max(UMACH_local,abs(test_vel))
      enddo ! dir

      if (is_rigid(im_primary).eq.1) then

       USOUND_local=1.0D+20  ! these materials do not compress

      else if (is_rigid(im_primary).eq.0) then

       if (fort_material_type(im_primary).eq.0) then
        USOUND_local=1.0D+20  ! these materials do not compress
       else if ((fort_material_type(im_primary).ge.1).and. &
                (fort_material_type(im_primary).le.MAX_NUM_EOS)) then
        den_local=den((im_primary-1)*num_state_material+1)
        temp_local=den((im_primary-1)*num_state_material+2)

        call init_massfrac_parm(den_local,massfrac_parm,im_primary)
        do ispec=1,num_species_var
         massfrac_parm(ispec)=den((im_primary-1)*num_state_material+2+ispec)
        enddo

        call INTERNAL_material(den_local,massfrac_parm, &
         temp_local,energy_local, &
         fort_material_type(im_primary),im_primary)
        call SOUNDSQR_material(den_local,massfrac_parm, &
         energy_local, &
         USOUND_local, &
         fort_material_type(im_primary),im_primary)
        USOUND_local=sqrt(USOUND_local)
        if (USOUND_local.le.zero) then
         print *,"USOUND_local invalid"
         stop
        endif
       else 
        print *,"fort_material_type(im_primary) invalid"
        stop
       endif
      else
       print *,"is_rigid invalid PROB.F90"
       stop
      endif
      mach=UMACH_local/USOUND_local

      return
      end subroutine get_mach_number 


       ! called from fort_init_icemask_and_icefacecut (GODUNOV_3D.F90) and
       !             tagexpansion (GODUNOV_3D.F90)
       ! input: LS, latent_heat, distribute_from_target
       ! output: icemask,icefacecut,im,im_opp,ireverse
      subroutine get_icemask_and_icefacecut( &
        nden, &
        xtarget, &
        time, &
        dx,bfact, &
        icemask, &
        icefacecut, &
        im, &
        im_opp, &
        im_primary, &
        ireverse, &
        denstate, &
        LS, &
        VOF, &
        distribute_from_target, &
        complement_flag)
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      integer, INTENT(in) :: nden
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(out) :: im
      integer, INTENT(out) :: im_opp
      integer, INTENT(out) :: im_primary
      integer, INTENT(out) :: ireverse
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(out) :: icemask
      real(amrex_real), INTENT(out) :: icefacecut
      real(amrex_real), INTENT(in) :: denstate(nden)
      real(amrex_real), INTENT(in) :: LS(num_materials)
      real(amrex_real), INTENT(in) :: VOF(num_materials)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      integer, INTENT(in) :: complement_flag

      real(amrex_real) dist_mask_override
      integer im_primary_vof
      integer im_secondary
      integer im_tertiary
      integer im_ice
      integer im_FSI_rigid
      integer im_dest,im_source
      integer iten
      real(amrex_real) LL(0:1)
      real(amrex_real) dxmax
      real(amrex_real) :: def_thermal(num_materials)
      integer :: im_local

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      do im_local=1,num_materials
       def_thermal(im_local)=room_temperature
      enddo

      if ((complement_flag.eq.0).or. &
          (complement_flag.eq.1)) then
       ! do nothing
      else
       print *,"complement_flag invalid"
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

      call get_dxmax(dx,bfact,dxmax)

        ! we are in "get_icemask_and_icefacecut"

      call get_primary_material(LS,im_primary)
      call get_primary_material_VFRAC(VOF,im_primary_vof)
      call get_secondary_material(LS,im_primary,im_secondary)

      if ((im_primary.ge.1).and. &
          (im_primary.le.num_materials)) then
       ! do nothing
      else
       print *,"im_primary invalid : ",im_primary
       stop
      endif

      if ((im_secondary.ge.1).and. &
          (im_secondary.le.num_materials)) then
        ! get_tertiary_material is declared in MOF.F90
        ! is_rigid(im_tertiary)=0
       call get_tertiary_material(LS,im_primary,im_secondary,im_tertiary)
      else
       print *,"im_secondary invalid : ",im_secondary
       stop
      endif

      if ((im_tertiary.ge.0).and.(im_tertiary.le.num_materials)) then
       ! do nothing
      else
       print *,"im_tertiary invalid: ",im_tertiary
       stop
      endif

      if (im_primary.eq.im_secondary) then
       print *,"im_primary.eq.im_secondary"
       stop
      else if (im_primary.eq.im_tertiary) then
       print *,"im_primary.eq.im_tertiary"
       stop
      else if (im_primary.lt.im_secondary) then
       im=im_primary
       im_opp=im_secondary
      else if (im_primary.gt.im_secondary) then
       im_opp=im_primary
       im=im_secondary
      else
       print *,"im_primary or im_secondary invalid"
       stop
      endif

      call get_iten(im,im_opp,iten)

      do ireverse=0,1
       LL(ireverse)= &
        get_user_latent_heat(iten+ireverse*num_interfaces,room_temperature,1)
      enddo

      if ((is_ice(im).eq.0).and. &
          (is_ice(im_opp).eq.0)) then
       im_ice=0
      else if ((is_ice(im).eq.1).and. &
               (is_ice(im_opp).eq.0)) then
       im_ice=im
      else if ((is_ice(im).eq.0).and. &
               (is_ice(im_opp).eq.1)) then
       im_ice=im_opp
      else if ((is_ice(im).eq.1).and. &
               (is_ice(im_opp).eq.1)) then
       im_ice=im_primary
      else
       print *,"is_ice invalid"
       stop
      endif
     
      if ((is_FSI_rigid(im).eq.0).and. &
          (is_FSI_rigid(im_opp).eq.0)) then
       im_FSI_rigid=0
      else if ((is_FSI_rigid(im).eq.1).and. &
               (is_FSI_rigid(im_opp).eq.0)) then
       im_FSI_rigid=im
      else if ((is_FSI_rigid(im).eq.0).and. &
               (is_FSI_rigid(im_opp).eq.1)) then
       im_FSI_rigid=im_opp
      else if ((is_FSI_rigid(im).eq.1).and. &
               (is_FSI_rigid(im_opp).eq.1)) then
       im_FSI_rigid=im_primary
      else
       print *,"is_FSI_rigid invalid"
       stop
      endif

      if (im_FSI_rigid.eq.im_primary) then

       ireverse=-1
       icemask=zero
       icefacecut=zero

      else if ((im_FSI_rigid.ge.0).and. &
               (im_FSI_rigid.le.num_materials).and. &
               (im_FSI_rigid.ne.im_primary)) then

        ! either the primary or secondary material is "ice"
       if ((im_ice.ge.1).and. &
           (im_ice.le.num_materials)) then 

         ! if the associated "melt" material
         ! is not the primary or secondary material,
         ! then check if it is the tertiary material.
        if ((LL(0).eq.zero).and.(LL(1).eq.zero)) then
         if ((im_tertiary.ge.1).and. &
             (im_tertiary.le.num_materials)) then
          if (is_rigid(im_tertiary).eq.0) then
           if (is_FSI_rigid(im_tertiary).eq.0) then
            if (is_ice(im_tertiary).eq.0) then
             if (im_ice.lt.im_tertiary) then
              im=im_ice
              im_opp=im_tertiary
             else if (im_ice.gt.im_tertiary) then
              im_opp=im_ice
              im=im_tertiary
             else
              print *,"im_ice or im_tertiary invalid"
              stop
             endif
             call get_iten(im,im_opp,iten)
             do ireverse=0,1
              LL(ireverse)= &
               get_user_latent_heat(iten+ireverse*num_interfaces, &
                   room_temperature,1)
             enddo
            else if (is_ice(im_tertiary).eq.1) then
             ! do nothing
            else
             print *,"is_ice(im_tertiary) invalid"
             stop
            endif
           else if (is_FSI_rigid(im_tertiary).eq.1) then
            ! do nothing
           else
            print *,"is_FSI_rigid(im_tertiary) invalid"
            stop
           endif
          else
           print *,"is_rigid(im_tertiary) invalid"
           print *,"contradiction with: get_tertiary_material"
           print *,"is_rigid(im_tertiary): ",is_rigid(im_tertiary)
           stop
          endif
         else if (im_tertiary.eq.0) then
          print *,"expecting im_tertiary>=1 and <=num_materials"
          print *,"im_tertiary: ",im_tertiary
          print *,"im_ice: ",im_ice
          print *,"im_FSI_rigid: ",im_FSI_rigid
          stop
         else
          print *,"im_tertiary invalid: ",im_tertiary
          stop
         endif
        else if ((LL(0).ne.zero).or.(LL(1).ne.zero)) then
         ! do nothing
        else
         print *,"LL invalid"
         stop
        endif

         ! an associated melt material was not found:
        if ((LL(0).eq.zero).and.(LL(1).eq.zero)) then

         ireverse=-1

         if (is_ice(im_primary).eq.1) then ! in ice bulk region

          icemask=zero
          icefacecut=zero

         else if (is_ice(im_primary).eq.0) then
          icemask=one
          icefacecut=one
         else
          print *,"is_ice(im_primary) invalid"
          print *,"im_primary ",im_primary
          print *,"is_ice(im_primary) ",is_ice(im_primary)
          stop
         endif

        else if ((LL(0).ne.zero).and.(LL(1).eq.zero)) then
         ireverse=0
         im_source=im
         im_dest=im_opp
        else if ((LL(0).eq.zero).and.(LL(1).ne.zero)) then
         ireverse=1
         im_source=im_opp
         im_dest=im
        else if ((LL(0).ne.zero).and.(LL(1).ne.zero)) then
         if (LS(im).ge.LS(im_opp)) then
          ireverse=0
          im_source=im
          im_dest=im_opp
         else if (LS(im_opp).ge.LS(im)) then
          ireverse=1
          im_source=im_opp
          im_dest=im
         else
          print *,"LS(im) or LS(im_opp) invalid"
          print *,"im, im_opp ",im,im_opp
          print *,"LS(im),LS(im_opp) ",LS(im),LS(im_opp)
          stop
         endif
        else
         print *,"LL invalid"
         print *,"LL(0),LL(1) ",LL(0),LL(1)
         stop
        endif
   
        if (ireverse.eq.-1) then
         ! do nothing
        else if ((ireverse.eq.0).or. &
                 (ireverse.eq.1)) then

         if (im_ice.eq.im_dest) then  ! freezing

          if (distribute_from_target(iten+num_interfaces*ireverse).eq.1) then
           ! do nothing
          else
           print *,"required freezing: "
           print *,"distribute_from_target(iten+num_interfaces*ireverse)=1"
           stop
          endif

          ! dist_mask_override>0 in the substrate. (dest is ice)  this routine
          ! tells one whether to force mask=0
          call icemask_override(xtarget,im_source,im_dest,dist_mask_override)

          if (dist_mask_override.ge.zero) then
           icemask=zero  ! mask off this cell.
          else if (dist_mask_override.le.zero) then
           if ((LS(im_ice).ge.zero).or. &
               (im_primary_vof.eq.im_ice).or. &
               (im_primary.eq.im_ice).or. &
               (VOF(im_ice).ge.0.5d0)) then
            icemask=zero
           else if ((LS(im_ice).le.zero).and. &
                    (im_primary_vof.ne.im_ice).and. &
                    (im_primary.ne.im_ice).and. &
                    (VOF(im_ice).le.0.5d0)) then
            icemask=one
           else
            print *,"LS(im_ice) bust"
            stop
           endif
          else
           print *,"dist_mask_override bust"
           stop
          endif

         else if (im_ice.eq.im_source) then ! melting

          if (distribute_from_target(iten+num_interfaces*ireverse).eq.0) then
           ! dist. from the ice to the liquid
          else
           print *,"required melting:"
           print *,"distribute_from_target(iten+num_interfaces*ireverse)=0"
           stop
          endif

          if ((LS(im_ice).ge.zero).or. &
              (im_primary_vof.eq.im_ice).or. &
              (im_primary.eq.im_ice).or. &
              (VOF(im_ice).ge.0.5d0)) then
           icemask=zero
          else if ((LS(im_ice).le.zero).and. &
                   (im_primary_vof.ne.im_ice).and. &
                   (im_primary.ne.im_ice).and. &
                   (VOF(im_ice).le.0.5d0)) then
           icemask=one
          else
           print *,"LS(im_ice) bust"
           stop
          endif

         else
          print *,"im_ice invalid"
          stop
         endif

         if ((icemask.eq.zero).or.(icemask.eq.one)) then
          ! do nothing
         else
          print *,"icemask invalid: ",icemask
          stop
         endif

         if (im_ice.eq.im_dest) then ! freezing

          if (dist_mask_override.ge.zero) then ! in a substrate
           icefacecut=zero
          else if (icemask.eq.zero) then
           icefacecut=zero
          else if (icemask.eq.one) then
           icefacecut=one
          else
           print *,"icemask invalid: ",icemask
           stop
          endif

         else if (im_ice.eq.im_source) then ! melting

          if (icemask.eq.zero) then

           icefacecut=zero

          else if (icemask.eq.one) then
           icefacecut=one
          else
           print *,"icemask invalid: ",icemask
           stop
          endif

         else
          print *,"im_ice invalid"
          stop
         endif

        else
         print *,"ireverse invalid"
         stop
        endif
 
       else if (im_ice.eq.0) then

        ireverse=-1
        icemask=one
        icefacecut=one

       else
        print *,"im_ice invalid:",im_ice
        stop
       endif

      else
       print *,"im_FSI_rigid invalid: ",im_FSI_rigid
       print *,"num_materials: ",num_materials
       print *,"im_primary: ",im_primary
       stop
      endif

      if ((icemask.eq.zero).or.(icemask.eq.one)) then
       ! do nothing
      else
       print *,"icemask invalid: ",icemask
       stop
      endif
   
      if ((icefacecut.ge.zero).and. &
          (icefacecut.le.one)) then
       ! do nothing
      else
       print *,"icefacecut invalid(6): ",icefacecut
       print *,"icemask=",icemask
       stop
      endif
      if (icemask.le.icefacecut) then
       ! do nothing
      else
       print *,"expecting icemask<=icefacecut"
       stop
      endif
 
      return
      end subroutine get_icemask_and_icefacecut

       ! MEHDI VAHAB HEAT SOURCE
       ! T^new=T^* + dt * (Q)/(rho cv)
       ! Q units: J/(m^3 s)
       ! called from: GODUNOV_3D.F90, subroutine fort_heatsource
       ! in fort_heatsource:
       ! T_local(im)=T_local(im)+ &
       !   dt*DeDTinverse(D_DECL(i,j,k),1)*heat_source_total im=1..num_materials

      subroutine get_local_heat_source( &
       time,dt, &
       x, &
       xsten, &  ! xsten(-nhalf:nhalf,SDIM)
       nhalf, &
       temperature_source, &  ! inputs file variable; default = 0.0
       temperature_source_cen, &
       temperature_source_rad, &
       LS,VFRAC,TEMPERATURE,DENSITY, &
       CV, &
       HEAT_SOURCE_OUT)
      use global_utility_module
      use global_distance_module
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use rigid_FSI_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time,dt
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: temperature_source
      real(amrex_real), INTENT(in) :: temperature_source_cen(SDIM)
      real(amrex_real), INTENT(in) :: temperature_source_rad(SDIM)
      real(amrex_real), INTENT(in) :: LS(num_materials)
      real(amrex_real), INTENT(in) :: VFRAC(num_materials)
      real(amrex_real), INTENT(in) :: TEMPERATURE(num_materials)
      real(amrex_real), INTENT(in) :: DENSITY(num_materials)
      real(amrex_real), INTENT(in) :: CV(num_materials)
      real(amrex_real), INTENT(out) :: HEAT_SOURCE_OUT(num_materials)
      real(amrex_real) dist
      real(amrex_real) dist_gas
      real(amrex_real) eta,depth,xs,ys,zs,phiE
      real(amrex_real) alpha_absorp,p_laser,r0,v_laser,xi,sigma_rad,eps_rad,hc
      real(amrex_real) T_chill
      real(amrex_real) time_unit
      integer im 
      real(amrex_real) vfrac_cutoff
      real(amrex_real) dist_cutoff
      integer localdir,inbox
      real(amrex_real) xlo,xhi

      ! eta: absorption efficiency
      ! phiE: electron beam diameter
      ! xs,ys,zs: center of electronic beam 

      vfrac_cutoff=0.1  ! was 0.99
      dist_cutoff=-1.0e+10 ! was 0.0

      depth=0.0062
      eta=0.9
      phiE=0.055 ! cm
      xs=0.5+7.0*sin(2.0*Pi*5.0*time)**3
      ys=0.5*floor(time*5.0)
      zs=zero

      alpha_absorp=3.0  ! was 0.3
      p_laser=100.0
      v_laser=1.25
      r0=3.5e-5
      xi=5.0e-6
      sigma_rad=5.67e-8
      eps_rad=1.0
      hc=100.0

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
      if (temperature_source.ge.zero) then
       ! do nothing
      else
       print *,"temperature_source invalid"
       stop
      endif
      do im=1,num_materials
       if ((VFRAC(im).ge.-EPS1).and. &
           (VFRAC(im).le.one+EPS1)) then
        ! do nothing
       else
        print *,"VFRAC invalid: ",VFRAC(im)
        stop
       endif 
       if (DENSITY(im).gt.zero) then
        ! do nothing
       else
        print *,"DENSITY must be positive"
        stop
       endif
       if (TEMPERATURE(im).gt.zero) then
        ! do nothing
       else
        print *,"TEMPERATURE must be positive"
        stop
       endif
       if (CV(im).gt.zero) then
        ! do nothing
       else
        print *,"CV must be positive"
        stop
       endif
       HEAT_SOURCE_OUT(im)=zero

       if (is_in_probtype_list().eq.1) then

        call SUB_HEATSOURCE( &
          im,VFRAC,time, &
          x, &
          xsten, & ! xsten(-nhalf:nhalf,SDIM)
          nhalf, &
          TEMPERATURE, &
          HEAT_SOURCE_OUT(im),DENSITY,CV,dt, &
          num_materials)

       else if (probtype.eq.401) then

        call HELIX_HEATSOURCE(im,VFRAC,time,x,TEMPERATURE, &
               HEAT_SOURCE_OUT(im),DENSITY,CV,dt)

       else if (probtype.eq.402) then
      
        call TSPRAY_HEATSOURCE(im,VFRAC,time,x,TEMPERATURE, &
               HEAT_SOURCE_OUT(im),DENSITY,CV,dt)

       else if (probtype.eq.533) then

        call rigid_FSI_HEATSOURCE(im,VFRAC,time,x,TEMPERATURE, &
            HEAT_SOURCE_OUT(im),DENSITY,CV,dt)

       else if (probtype.eq.311) then ! sample user defined call

        call USERDEF_HEATSOURCE(im,VFRAC,time,x,TEMPERATURE, &
               HEAT_SOURCE_OUT(im),DENSITY,CV,dt)

       else if ((probtype.eq.299).or. &
                (probtype.eq.301)) then ! melting or additive manufacturing
        if (im.eq.3) then ! metal
         if ((radblob3.lt.zero).or.(radblob3.gt.radblob)) then
          print *,"radblob3 invalid"
          stop
         endif
         if (probtype.eq.299) then
          call INIT_LS_SOLID_MELT(x(1),x(2),x(SDIM),time,dist)  !>0 in solid
          call INIT_LS_GAS_MELT(x(1),x(2),x(SDIM),time,dist_gas)  !>0 in gas
         else if (probtype.eq.301) then
          call INIT_LS_SOLID_AM(x(1),x(2),x(SDIM),time,dist)  !>0 in solid
          call INIT_LS_GAS_AM(x(1),x(2),x(SDIM),time,dist_gas)  !>0 in gas
         else
          print *,"probtype invalid"
          stop
         endif

         if (axis_dir.eq.0) then

          if ((dist.gt.radblob-radblob3).and. &
              (vfrac(im).ge.vfrac_cutoff)) then
           if (TEMPERATURE(im).lt.fort_tempconst(im)) then
            HEAT_SOURCE_OUT(im)=(fort_tempconst(im)-TEMPERATURE(im))* &
             DENSITY(im)*CV(im)/dt 
           endif 
          endif ! dist>radblob-radblob3

         else if (axis_dir.eq.1) then

          if ((dist.gt.dist_cutoff).and. &
              (VFRAC(im).ge.vfrac_cutoff)) then 
!           time_unit=dt
           time_unit=one
           HEAT_SOURCE_OUT(im)=(alpha_absorp*p_laser/(Pi*r0*r0)* & 
            exp(-2.0*(x(1)-50.0e-6)**2/(r0*r0)))
           !- & 
           !sigma_rad*eps_rad*(TEMPERATURE(im)**4-300.0**4)- &
           !hc*(TEMPERATURE(im)-300.0))*DENSITY(im)*CV(im)/dt

           !INTENSITY=max(0.0, (-2.25*dist*dist/(depth*depth) & 
           !           +1.5*dist/depth+0.75)/0.75)
           !HS=exp(-2.0*(x-xs)*(x-xs)/(phiE*phiE))           
           !HEAT_SOURCE_OUT(im)=2.0*60.0*6.7/(Pi*phiE*phiE)* &
           !     DENSITY(im)*CV(im)*eta*HS*INTENSITY/depth/dt
          endif ! dist>dist_cutoff and vfrac>vfrac_cutoff?
         else if (axis_dir.eq.2) then

          if (VFRAC(im).ge.vfrac_cutoff) then
!           time_unit=dt
           time_unit=one
           HEAT_SOURCE_OUT(im)=alpha_absorp*p_laser/(Pi*r0*r0*r0)* &
            exp(-2.0*(x(1)-time*v_laser-xi)**2/(r0*r0)) - &
            sigma_rad*eps_rad*(TEMPERATURE(im)**4-300.0**4)- &
            hc*(TEMPERATURE(im)-300.0)
           if (HEAT_SOURCE_OUT(im).lt.zero) then
            HEAT_SOURCE_OUT(im)=zero
           endif
          endif ! vfrac>vfrac_cutoff?

         else if (axis_dir.eq.3) then

          if ((dist.gt.radblob-radblob3).and. &
              (vfrac(im).ge.vfrac_cutoff)) then
           if (TEMPERATURE(im).lt.fort_tempconst(im)) then
            HEAT_SOURCE_OUT(im)=(fort_tempconst(im)-TEMPERATURE(im))* &
             DENSITY(im)*CV(im)/dt
           endif
          endif ! dist>radblob-radblob3

         else
          print *,"axis_dir invalid"
          stop
         endif

        else if (im.eq.1) then ! melt

         if (axis_dir.eq.0) then
          ! do nothing
         else if (axis_dir.eq.1) then

          if ((radblob3.lt.zero).or.(radblob3.gt.radblob)) then
           print *,"radblob3 invalid"
           stop
          endif
          if (probtype.eq.299) then
           call INIT_LS_LIQUID_MELT(x(1),x(2),x(SDIM),time,dist)  !>0 in melt
          else if (probtype.eq.301) then
           call INIT_LS_LIQUID_AM(x(1),x(2),x(SDIM),time,dist)  !>0 in melt
          else
           print *,"probtype invalid"
           stop
          endif

           ! inside the liquid metal
          if ((dist.gt.dist_cutoff).and. &
              (vfrac(im).ge.vfrac_cutoff)) then 
           HEAT_SOURCE_OUT(im)=(alpha_absorp*p_laser/(Pi*r0*r0)* &
              exp(-2.0*(x(1)-50.0e-6)**2/(r0*r0)))
          endif ! dist>dist_cutoff and F>vfrac_cutoff in melt 

         else if (axis_dir.eq.2) then

          if ((radblob3.lt.zero).or.(radblob3.gt.radblob)) then
           print *,"radblob3 invalid"
           stop
          endif
          if (probtype.eq.299) then
           call INIT_LS_LIQUID_MELT(x(1),x(2),x(SDIM),time,dist)  !>0 in melt
          else if (probtype.eq.301) then
           call INIT_LS_LIQUID_AM(x(1),x(2),x(SDIM),time,dist)  !>0 in melt
          else
           print *,"probtype invalid"
           stop
          endif

           ! inside the liquid metal
          if ((dist.gt.dist_cutoff).and. &
              (vfrac(im).ge.vfrac_cutoff)) then

           HEAT_SOURCE_OUT(im)=alpha_absorp*p_laser/(Pi*r0*r0*r0)* &
            exp(-2.0*(x(1)-time*v_laser-xi)**2/(r0*r0)) - &
            sigma_rad*eps_rad*(TEMPERATURE(im)**4-300.0**4)- &
            hc*(TEMPERATURE(im)-300.0)
           if (HEAT_SOURCE_OUT(im).lt.zero) then
            HEAT_SOURCE_OUT(im)=zero
           endif
          endif ! dist>dist_cutoff and F>vfrac_cutoff in melt 

         else if (axis_dir.eq.3) then
          ! do nothing
         else
          print *,"axis_dir invalid"
          stop
         endif

        else if (im.eq.2) then ! air

         if (axis_dir.eq.0) then
          ! do nothing
         else if (axis_dir.eq.1) then
          ! do nothing
         else if (axis_dir.eq.2) then
          ! do nothing
         else if (axis_dir.eq.3) then
          
          if ((dist_gas.gt.radblob4).and. &
              (vfrac(im).ge.vfrac_cutoff)) then
           T_chill=250.0
           if (TEMPERATURE(im).gt.T_chill) then
            HEAT_SOURCE_OUT(im)=(T_chill-TEMPERATURE(im))* &
             DENSITY(im)*CV(im)/dt 
           endif 
          endif ! dist_gas>radblob4

         else
          print *,"axis_dir invalid"
          stop
         endif

        else
         print *,"im invalid61"
         stop
        endif 
       endif ! probtype==299 or probtype==301

       if (temperature_source.eq.zero) then
        ! do nothing
       else if (temperature_source.gt.zero) then
        inbox=1
        do localdir=1,SDIM
         xlo=temperature_source_cen(localdir)-temperature_source_rad(localdir)
         xhi=temperature_source_cen(localdir)+temperature_source_rad(localdir)
         if ((x(localdir).lt.xlo).or.(x(localdir).gt.xhi)) then
          inbox=0
         endif
        enddo

        if (inbox.eq.1) then
         if (TEMPERATURE(im).lt.temperature_source) then
          HEAT_SOURCE_OUT(im)=(temperature_source-TEMPERATURE(im))* &
             DENSITY(im)*CV(im)/dt 
         else if (TEMPERATURE(im).ge.temperature_source) then
          ! do nothing
         else
          print *,"TEMPERATURE(im) invalid"
          stop
         endif
        else if (inbox.eq.0) then
         ! do nothing
        else
         print *,"inbox invalid"
         stop
        endif
     
       else
        print *,"temperature_source invalid"
        stop
       endif

       if (1.eq.0) then
        if (HEAT_SOURCE_OUT(im).ne.zero) then
         print *,"x,im,heat_source_out ",xsten(0,1),xsten(0,2), &
           im,HEAT_SOURCE_OUT(im)
        endif
       endif
      enddo ! im=1..num_materials

      return
      end subroutine get_local_heat_source 


       !check_user_defined_velbc is called from fort_estdt (GODUNOV_3D.F90)
      subroutine check_user_defined_velbc(time,dir,uu,dx)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: uu
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) utest,uscale

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid"
       stop
      endif

      if (is_in_probtype_list().eq.1) then

       call SUB_CFL_HELPER(time,dir,uu,dx)

      else

       if (dir.eq.adv_dir-1) then
        if (probtype.eq.32) then ! flow past cylinder
         uu=max(abs(uu),abs(adv_vel))
         uu=max(abs(uu),abs(advbot))
        endif
        if (probtype.eq.41) then ! pipe flow
         uu=max(abs(uu),abs(adv_vel))
         uu=max(abs(uu),abs(advbot))
         uu=max(abs(uu),abs(vinletgas))
        endif
       endif
       
       if (dir.eq.0) then
        if ((probtype.eq.1).and.(axis_dir.eq.15)) then
         ! u = U t
         ! U(t+k)k=h 
         ! k^2+k t = h/U
         ! k=(-t+sqrt(t^2+4h/U))/2=(2h/U)/(t+sqrt(t^2+4h/U))
         ! h/k=(t+sqrt(t^2+4h/U))U/2
         uscale=abs(adv_vel)/two
         if (time.lt.two) then
          utest=(time+sqrt(time**2+four*dx(1)/uscale))
          utest=utest*uscale/two
         else
          utest=abs(adv_vel)
         endif
         uu=max(abs(uu),utest)
        else if ((probtype.eq.1).and.(axis_dir.lt.150).and. &
                 (axis_dir.ge.0)) then
         uu=max(abs(uu),abs(vinletgas))
         if (SDIM.eq.2) then
          uu=max(abs(uu),abs(adv_vel))
         endif
        endif
        if (probtype.eq.5700) then
         uu=max(abs(uu),abs(vinletgas))
        endif
       else if (dir.eq.1) then
 
         ! advbot is inflow vel at base 
        if (SDIM.eq.2) then
         if ((probtype.eq.531).or. &
             (probtype.eq.538).or. &  ! inputs.injA
             (probtype.eq.53).or. &
             (probtype.eq.541)) then
          uu=max(abs(uu),abs(advbot))
         endif
        endif
        if (probtype.eq.5700) then
         uu=max(abs(uu),abs(advbot))
        endif
        if (probtype.eq.5700) then
         uu=max(abs(uu),abs(vinletgas))
        endif

       else if ((dir.eq.2).and.(SDIM.eq.3)) then
        if (probtype.eq.36) then
         uu=max(abs(uu),abs(yblob9))
         uu=max(abs(uu),abs(yblob10))
        endif
         ! advbot is inflow vel at base  
        if ((probtype.eq.537).or. &
            (probtype.eq.538).or. &
            (probtype.eq.53).or. &
            (probtype.eq.541)) then
         uu=max(abs(uu),abs(advbot))
        endif
       else
        print *,"dir invalid check user defined velbc"
        stop
       endif

       ! vinletgas is a prescribed boiling rate
       if (probtype.eq.801) then
        uu=max(abs(uu),abs(vinletgas)*fort_denconst(1)/fort_denconst(2))
       endif
       if ((probtype.eq.701).and.(axis_dir.eq.2)) then
        uu=max(abs(uu),abs(adv_vel))
        uu=max(abs(uu),yblob*two*Pi)
       endif 
       if (probtype.eq.201) then
        uu=max(abs(uu),abs(advbot))
       endif 
      endif

      return
      end subroutine check_user_defined_velbc

      subroutine Zuzio_velbc(time,dir,side,veldir, &
       vel_out,vel_in)
      IMPLICIT NONE

      real(amrex_real) time,vel_out,vel_in
      integer dir,side,veldir

      vel_out=vel_in
      if ((dir.eq.1).and.(veldir.eq.1)) then
       vel_out=min(adv_vel,time*adv_vel/two)
      else if (dir.eq.1) then
       vel_out=zero
      else if (dir.eq.veldir) then
       vel_out=zero
      endif

      return
      end subroutine Zuzio_velbc


        ! fractional error always less than 1.2E-7
      real(amrex_real) function erf_recipe(x)
      IMPLICIT NONE

      real(amrex_real) dumerfc, x
      real(amrex_real) t, z

      z = abs(x)
      t = one / ( one + half * z )

      dumerfc = t * exp(-z * z - 1.26551223d0 + t *   &
        ( 1.00002368d0 + t * ( 0.37409196d0 + t *     &
        ( 0.09678418d0 + t * (-0.18628806d0 + t *     &
        ( 0.27886807d0 + t * (-1.13520398d0 + t *     &
        ( 1.48851587d0 + t * (-0.82215223d0 + t * 0.17087277d0 )))))))))

      if ( x.lt.zero ) dumerfc = two - dumerfc
     
      erf_recipe = one - dumerfc

      return
      end function erf_recipe

      real(amrex_real) function transf(x,cp)
      IMPLICIT NONE

      real(amrex_real) x,pi,cp

      if (cp.le.zero) then
       print *,"cp invalid"
       stop
      endif

      pi=four*atan(one)
      transf=x*exp(x*x)*erf_recipe(x)-cp/sqrt(pi)

      return
      end function transf

       ! ice is supercooled, liquid is saturation temperature
      subroutine stefan1D(SUBLEN,TSUBSTRATE,TSAT,CP,L,K,den,time,rad)
      IMPLICIT NONE

      real(amrex_real) SUBLEN,tcrit
      real(amrex_real) TSUBSTRATE,TSAT,CP,L,K,den,time,rad
      real(amrex_real) a,b,fb,fa,c,fc,lambda,alpha,fact
      integer i

      if (TSAT.le.TSUBSTRATE) then
       print *,"TSAT invalid in stefan1D"
       stop
      endif
      if ((L.le.zero).or.(den.le.zero).or.(CP.le.zero).or.(K.le.zero).or. &
          (TSUBSTRATE.le.zero).or.(time.le.zero)) then
       print *,"stefan1D parameters invalid"
       stop
      endif
      fact=CP*(TSAT-TSUBSTRATE)/L
      alpha=K/(den*CP)

      a=zero
      b=two
      fb=transf(b,fact)
      do while (fb.le.zero) 
       b=2.0*b
       fb=transf(b,fact)
      enddo 
      fa=transf(a,fact)
      fb=transf(b,fact)
      if (fa*fb.gt.zero) then
       print *,"a,b invalid"
      endif
      do i=1,20
       c=(a+b)/two  
       fc=transf(c,fact)
       if (fa*fc.gt.zero) then
        a=c
        fa=fc
       else
        b=c
        fb=fc
       endif
      enddo
      lambda=c
       ! R=2 lam sqrt(alpha t)
       ! (R/(2 lam))^2 = alpha t
       ! t=(R/(2 lam))^2 / alpha
      if ((lambda.le.zero).or.(alpha.le.zero).or. &
          (SUBLEN.lt.zero)) then
       print *,"parameters invalid"
       stop
      endif
      tcrit=(SUBLEN/(two*lambda))**2/alpha 
      rad=two*lambda*sqrt(alpha*(tcrit+time))-SUBLEN
      if (rad.le.zero) then
       print *,"rad invalid"
       stop
      endif

      return
      end subroutine stefan1D


! for bering straight, radblob=sillrad, yblob=width, xblob=length
!  zblob=height  zblob-twall=free surface height
! height= height of sill
      subroutine beringfree(x,y,z,dist)
      IMPLICIT NONE


      real(amrex_real) x,y,z,dist

      dist=zblob-twall-z

      return
      end subroutine beringfree
      
      subroutine bering(x,y,z,height,width,length,sillrad,dist)
      use global_utility_module
      use global_distance_module

      IMPLICIT NONE

      real(amrex_real) x,y,z,dist
      real(amrex_real) height,width,length,sillrad
      real(amrex_real) midx,midy,yy,maxx,dist1
      real(amrex_real) hy,hz,ylo,yhi,zlo,zhi

      midx=half*length
      midy=half*width
      maxx=midx+height

      yy=y
      if (y.gt.midy) then
       yy=two*midy-y
      endif

      if (x.le.midx) then
       if ((z.le.height).or.(yy.le.midy-sillrad)) then
        dist=midx-x
       else 
        dist=sqrt( (midx-x)**2+(yy-midy+sillrad)**2 )
       endif
      else 
       hz=height-(x-midx)
       hy=midy-sillrad-(x-midx)  
       ylo=hy
       yhi=midy+sillrad+(x-midx)
       zlo=hz
       zhi=99999.0
       call squaredist(yy,z,ylo,yhi,zlo,zhi,dist1)
       dist=-dist1
       if (dist.lt.zero) then
        dist1=midx-x
        if (dist1.gt.dist) then
         dist=dist1
        endif
       endif
      endif

      return
      end subroutine bering



! XIAOYI LI: Adding turbulent boundary layer for gas flow      
      subroutine BL_inlet(vel,z,time)
      IMPLICIT NONE
      real(amrex_real) z,vel,time
      real(amrex_real) TBLT
      real(amrex_real) freq,amp,time0,velmax

      time0=0.0
      freq=0.5
      amp=0.1*adv_vel

!      velmax=adv_vel+amp*sin(2.0*3.1415926*(time-time0)*freq)
      velmax=adv_vel

      TBLT=1.0

      if(z < TBLT) then
       vel = velmax*(z/TBLT)**(1/7.0)
      else
       vel = velmax
      endif

!      if (z>2.55) print *, "Imposed gas flow = ", velmax                                            
      return
      end subroutine BL_inlet



! unit are in meters

      subroutine boatparms(depth,waterdepth,height,length,botwidth, &
         xback,yback,zback)
      IMPLICIT NONE
      real(amrex_real) depth,waterdepth,height,length,botwidth, &
             xback,yback,zback


      if (probtype.eq.9) then  ! ship wave with solid parms and
                               ! free surface parms prescribed here ?

      if (axis_dir.eq.0) then
       depth=1.0
       waterdepth=13.0
       height=7.0
       length=21.0
       botwidth=4.0
       xback=40.0
       yback=15.0
       zback=waterdepth-depth
         ! solid comes from .cas file
      else if ((axis_dir.eq.1).or.(axis_dir.eq.2).or. &  
               (axis_dir.eq.3)) then
       depth=-0.1
       waterdepth=0.0
       height=0.1
       length=1.0
       botwidth=0.1
       xback=1.0
       yback=0.1
       zback=waterdepth-depth
      else
       print *,"invalid boat selection"
       stop
      endif 

      else
       print *,"probtype should be 9"
       stop
      endif

      return
      end subroutine boatparms

      subroutine boatdist(x,y,z,dist)
      IMPLICIT NONE
      real(amrex_real) x,y,z,dist
      real(amrex_real) depth,waterdepth,height,length,botwidth, &
             xback,yback,zback

      call boatparms(depth,waterdepth,height,length,botwidth, &
             xback,yback,zback) 
      dist=waterdepth-z ! waterdepth=0 if axis_dir=1
      if (axis_dir .eq. 3) then
       if (x .gt. -0.65 .and. x .lt. -0.49 .and. &
           y .gt. -0.02 .and. y .lt. 0.02         ) then
        dist = -0.03 - z
       endif
      endif
 
      return
      end subroutine boatdist

       ! polymer_factor=1/maximum_extensibility=1/L
      subroutine get_mod_elastic_time(elastic_time,traceA, &
          polymer_factor,modtime)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: elastic_time
      real(amrex_real), INTENT(in) :: traceA
      real(amrex_real), INTENT(in) :: polymer_factor ! 1/max_extend
      real(amrex_real), INTENT(out) :: modtime

      if (polymer_factor.ge.zero) then
       ! do nothing
      else
       print *,"polymer_factor invalid"
       stop
      endif
      if (elastic_time.ge.zero) then
       ! do nothing
      else
       print *,"elastic_time invalid  elastic_time=",elastic_time
       stop
      endif
      if (traceA.gt.zero) then
       ! do nothing
      else
       print *,"trace A must be postive since A is SPD, traceA=",traceA
       stop
      endif

      modtime=elastic_time*(one-traceA*(polymer_factor**2))
      if (modtime.ge.zero) then
       ! do nothing
      else if (modtime.lt.zero) then
       modtime=zero
      else
       print *,"modtime invalid"
       stop
      endif

      return
      end subroutine get_mod_elastic_time

      subroutine viscosity(system,vis,shear)
      IMPLICIT NONE
!
      integer   system
!
      real(amrex_real)    zeroviscosity,  contviscosity
      real(amrex_real)    parameterA,  parameterB,  parameterN
      real(amrex_real)    vis, shear, shearmax,etacutoff
      real(amrex_real)    term1,powerterm
!
!     * This is a subroutine for setting pseudoplastic viscosity
!
!     * Generalized Cross-Carreau model parameters (4 parameters)
!     * Zero sheare rate viscosity (cgs)
!     * theree parameters a, B, n (-),(s),(-)

      shearmax=1.0D+20
      etacutoff=zero

      if (system.eq.1) then
        zeroviscosity = 0.5
        parameterA  = 2.0
        parameterB  = 1.0
        parameterN  = 0.5
      else if (system.eq.2) then
        zeroviscosity = 0.5
        parameterA  = 2.0
        parameterB  = 0.01
        parameterN  = 0.5
      else if (system.eq.3) then
        zeroviscosity = 0.5
        parameterA  = 2.0
        parameterB  = 1.0
        parameterN  = 0.8
      else if (system.eq.4) then
        zeroviscosity = 0.5
        parameterA  = 2.0
        parameterB  = 0.01
        parameterN  = 0.8
      else if (system.eq.5) then
        zeroviscosity = 0.31
        parameterA  = 5.0
        parameterB  = 0.167
        parameterN  = 0.875
      else if (system.eq.6) then
        zeroviscosity = 0.5
        parameterA  = 2.5
        parameterB  = 0.9
        parameterN  = 0.25
        shearmax=153.6
        etacutoff=1.24D-2
      else if (system.eq.7) then
        zeroviscosity = 0.26
        parameterA  = 4.0
        parameterB  = 1.429
        parameterN  = 0.4
        shearmax=103.1
        etacutoff=1.30D-2
      else
        print *,"system invalid"
        stop
      endif

       if (parameterN.gt.one) then
        print *,"parameterN invalid"
        stop
       endif

       if (shear.ge.shearmax) then
        contviscosity=etacutoff
       else
        term1=one+(parameterB*shear)**parameterA
        powerterm=(parameterN-one)/parameterA
        contviscosity=zeroviscosity*(term1**powerterm)
       endif
       vis=contviscosity

      return
      end subroutine viscosity


      subroutine get_max_user_tension(tension,new_tension)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: tension(num_interfaces)
      real(amrex_real), INTENT(out) :: new_tension(num_interfaces)
      integer iten

      do iten=1,num_interfaces
       new_tension(iten)=tension(iten)
      enddo

      return
      end subroutine get_max_user_tension



 
       subroutine get_scaled_tension(tension_in,tension_out)
       IMPLICIT NONE

       real(amrex_real), INTENT(in) :: tension_in
       real(amrex_real), INTENT(out) :: tension_out


       tension_out=tension_in/global_pressure_scale

       return 
       end subroutine get_scaled_tension


!assume x_proj locates in centeral cell ((ni+1)/2,(nj+1)/2,(nk+1)/2)
!The stencil is better to be 7x7x7 (or larger), 
!because the range of delta function in GNBC model is 2*dx, 
!if cannot find contact line in stencil, or closest_distance > 2*dx, or 
!the cell cloest to CL is (i = 1 or ni, or j = ..., or k = ...(dim=3)), 
!then return?
!LS1_xp: fluid level set value of projection point
!x(ni, nj, nk, dim): positions of stencil
!x_proj(dim): position of projection point
!dx: cell scale
!if dim = 2, then nk = 1
subroutine closest_distance_to_CL( &
     LS_stencil, &
     LS_xp, &! bilinear interp of LS to x_projection
     x_stencil, & ! stencil of x values
     x_proj, &
     n_rad, &
     actual_angle, &
     closest_distance, &
     im_primary, & ! im_primary owns adjoining cell to given structure cell
     im_secondary, & ! im_secondary is a fluid too.
     im_solid, &  ! structure material id.
     prob_dim)
use global_utility_module
implicit none

integer, INTENT(in) :: n_rad, prob_dim
integer, INTENT(in) :: im_primary,im_secondary,im_solid
real(amrex_real), INTENT(in) :: &
   LS_stencil(-n_rad:n_rad,-n_rad:n_rad,-n_rad:n_rad,num_materials)
real(amrex_real), INTENT(in) :: &
   x_stencil(-n_rad:n_rad,-n_rad:n_rad,-n_rad:n_rad,prob_dim)
real(amrex_real), INTENT(in) :: x_proj(prob_dim)
real(amrex_real), INTENT(in) :: LS_xp(num_materials)

real(amrex_real), INTENT(out) :: actual_angle, closest_distance

real(amrex_real) :: dx

integer dir
integer i, j, k, d, i_method, icl, jcl, kcl
integer find_cl ! estimate if contact line exits
real(amrex_real) costheta, eps, dis, mag, phimin, tmp(3), tmp1(3), &
                 nphi(3), npsi(3), nphi_xp(3), tpsi(3), nalpha(3), &
                 x_inf_proj(prob_dim), x_inf_proj_alpha(prob_dim), &
                 x_psi_proj(prob_dim), x_contact_point(prob_dim)
!nphi: normal vector of contact line, grad(LS1)/|grad(LS1)|
!npsi: normal vector of substrate, grad(LS3)/|grad(LS3)|
!nphi_xp: unit vector in gradient direction of projection point, 
!    grad(LS1)/|grad(LS1)|, assume it equal the value at centeral cell
!tpsi: unit vector which is tangent to substrate, on plane_alpha 
!    (constructed by npsi and nphi_xp), point outward to interface
!nalpha: normal vector of plane_alpha
!x_inf_proj: projection point of (icl,jcl,kcl) on interface
!x_inf_proj_alpha: projection point of x_inf_proj_alpha on plane_alpha
!x_psi_proj: projection point of x_inf_proj_alpha on substrate
!x_contact_point: location of closest contact point

    if (n_rad.eq.3) then
     ! do nothing
    else
     print *,"n_rad invalid"
     stop
    endif

    dir=1
    dx=0.5d0*(x_stencil(1,0,0,dir)-x_stencil(-1,0,0,dir))
    i_method = 2
    find_cl = 0
    eps = 1.1d0*dx
    closest_distance = 1.d10
    actual_angle = 0.d0

!calculate normal vector of substrate
!this used to be LS3
    dir=1
    tmp(dir) = (LS_stencil(1,0,0,im_solid)-LS_stencil(-1,0,0,im_solid))/ &
               (x_stencil(1,0,0,dir)-x_stencil(-1,0,0,dir))
    dir=2
    tmp(dir) = (LS_stencil(0,1,0,im_solid)-LS_stencil(0,-1,0,im_solid))/ &
               (x_stencil(0,1,0,dir)-x_stencil(0,-1,0,dir))
    if (prob_dim .eq. 3) then
     dir=prob_dim
     tmp(dir) = (LS_stencil(0,0,1,im_solid)-LS_stencil(0,0,-1,im_solid))/ &
                (x_stencil(0,0,1,dir)-x_stencil(0,0,-1,dir))
    else
       tmp(3) = 0.d0
    endif
    mag = 0.d0
    do d = 1, prob_dim
       mag = mag+tmp(d)*tmp(d)
    enddo
    mag = sqrt(mag)
    do d = 1, prob_dim
       npsi(d) = tmp(d)/mag
    enddo

!find closest distance and actual angle

!method 1: approximate value, find closest cell center (with different LS 
!          sign and close to substrate)
    if (i_method .eq. 1) then
       closest_distance = 1.d10
       do k = -n_rad,n_rad
       do j = -n_rad,n_rad
       do i = -n_rad,n_rad
        if ((LS_stencil(i,j,k,im_primary)*LS_xp(im_primary).le.0.d0).and. &
            (abs(LS_stencil(i,j,k,im_solid)).lt.eps)) then
         find_cl = 1
         dis = 0.d0
         do d = 1, prob_dim
          dis = dis+(x_stencil(i,j,k,d)-x_proj(d))**2.d0
         enddo
         dis = sqrt(dis)
         if (dis .lt. closest_distance) then
          closest_distance = dis
          if (LS_stencil(i,j,k,im_solid) .lt. 0.d0) then
           icl = i
           jcl = j
           kcl = k
          endif
         endif
        endif
       enddo
       enddo
       enddo
       if (find_cl .ne. 1) then
          print *, "Cannot find contact line!"
          return
       endif
       if (icl .eq. -n_rad .or. icl .eq. n_rad .or. &
           jcl .eq. -n_rad .or. jcl .eq. n_rad .or. &
           (prob_dim .eq. 3 .and. &
           (kcl .eq. -n_rad .or. kcl .eq. n_rad))) then
          print *, "Too far from contact line!"
          return
       endif
!calculate normal vector of contact line
!this used to be LS1
       dir=1
       tmp(dir) = (LS_stencil(icl+1,jcl,kcl,im_primary)- &
                   LS_stencil(-1+icl,jcl,kcl,im_primary))/ &
         (x_stencil(icl+1,jcl,kcl,dir)-x_stencil(-1+icl,jcl,kcl,dir))
       dir=2
       tmp(dir) = (LS_stencil(icl,jcl+1,kcl,im_primary)- &
                   LS_stencil(icl,-1+jcl,kcl,im_primary))/ &
         (x_stencil(icl,jcl+1,kcl,dir)-x_stencil(icl,-1+jcl,kcl,dir))
       if (prob_dim .eq. 3) then
        dir=prob_dim
        tmp(dir) = (LS_stencil(icl,jcl,kcl+1,im_primary)- &
                    LS_stencil(icl,jcl,-1+kcl,im_primary))/ &
         (x_stencil(icl,jcl,kcl+1,dir)-x_stencil(icl,jcl,-1+kcl,dir))
       else
          tmp(3) = 0.d0
       endif
       mag = 0.d0
       do d = 1, prob_dim
          mag = mag+tmp(d)*tmp(d)
       enddo
       mag = sqrt(mag)
!calculate actual contact angle
       costheta = 0.d0
       do d = 1, prob_dim
          nphi(d) = tmp(d)/mag
          costheta = costheta + nphi(d) * npsi(d)
       enddo
       actual_angle = acos(costheta)

!method 2:
!  define: phi = LS1, psi = LS3
!  step1: calculate nphi_xp based on LS1 value in centeral cell, 
!         with nphi_xp and npsi, determine plane-alpha (which is 
!         perpendicular to substrate and interface)
!  step2: calculate tpsi on substrate and plane-alpha, with point x_proj
!         and vector tpsi, we can determine line L (contact point is on L)
!  step3: when -eps<psi<0 and |phi|<eps and |(x-x_proj) X tpsi|<eps,
!         find (icl,jcl,kcl) makes |phi| minimum
!  step4: calculate nphi based on LS1 value in (icl,jcl,kcl)
!  step5: calculate actual_angle based on nphi and npsi
!  step6: project (icl,jcl,kcl) to interface based on nphi and 
!         LS1(icl,jcl,kcl), get point x_inf_proj
!  step7: project x_inf_proj to plane_alpha, get x_inf_proj_alpha
!  step8: get distance between x_inf_proj_alpha and line L, dis, and 
!         intersection x_psi_proj
!  step9:get x_contact_point = x_psi_proj + dis/tan(actual_angle) * tpsi
!  step10:get closest_distance = |x_proj - x_contact_point|

    else if (i_method .eq. 2) then

!calculate unit vector in gradient direction of projection point, nphi_xp
       icl=0
       jcl=0
       kcl=0
       dir=1
         ! used to be LS1
       tmp(dir) = (LS_stencil(icl+1,jcl,kcl,im_primary)- &
                   LS_stencil(-1+icl,jcl,kcl,im_primary))/ &
         (x_stencil(icl+1,jcl,kcl,dir)-x_stencil(-1+icl,jcl,kcl,dir))
       dir=2
       tmp(dir) = (LS_stencil(icl,jcl+1,kcl,im_primary)- &
                   LS_stencil(icl,-1+jcl,kcl,im_primary))/ &
         (x_stencil(icl,jcl+1,kcl,dir)-x_stencil(icl,-1+jcl,kcl,dir))
       if (prob_dim .eq. 3) then
        dir=prob_dim
        tmp(dir) = (LS_stencil(icl,jcl,kcl+1,im_primary)- &
                    LS_stencil(icl,jcl,-1+kcl,im_primary))/ &
         (x_stencil(icl,jcl,kcl+1,dir)-x_stencil(icl,jcl,-1+kcl,dir))
       else
          tmp(3) = 0.d0
       endif

       mag = 0.d0
       do d = 1, prob_dim
          mag = mag+tmp(d)*tmp(d)
       enddo
       mag = sqrt(mag)
       do d = 1, prob_dim
          nphi_xp(d) = tmp(d)/mag
       enddo
!calculate tpsi
!  nalpha = nphi_xp X npsi (X: cross product)
!  nalpha = nalpha / |nalpha|
!  tpsi = nalpha X npsi
!  tpsi = tpsi / |tpsi|
       call crossprod(nphi_xp, npsi, nalpha)
       mag = 0.d0
       do d = 1, prob_dim
          mag = mag+nalpha(d)*nalpha(d)
       enddo
       mag = sqrt(mag)
       do d = 1, prob_dim
          nalpha(d) = nalpha(d)/mag
       enddo
       call crossprod(nalpha, npsi, tpsi)
       mag = 0.d0
       do d = 1, prob_dim
          mag = mag+tpsi(d)*tpsi(d)
       enddo
       mag = sqrt(mag)
       do d = 1, prob_dim
          tpsi(d) = tpsi(d)/mag
       enddo
!find icl, jcl, kcl
       phimin = 1.d10
       do k = -n_rad,n_rad
       do j = -n_rad,n_rad
       do i = -n_rad,n_rad
        do d = 1, 2
         tmp(d) = x_stencil(i,j,k,d) - x_proj(d)
        enddo
        if (prob_dim .eq. 3) then
         tmp(3) = x_stencil(i,j,k,prob_dim) - x_proj(prob_dim)
        else
         tmp(3) = 0.d0
        endif
        call crossprod(tmp, tpsi, tmp1)
        dis = 0.d0
        do d = 1, prob_dim
         dis = dis+tmp1(d)*tmp1(d)
        enddo
        dis = sqrt(dis)
         ! used to be LS3
        if (LS_stencil(i,j,k,im_solid) .le. 0.d0 .and. &
            LS_stencil(i,j,k,im_solid) .gt. -eps .and. &
            abs(LS_stencil(i,j,k,im_primary)) .lt. eps .and. &
            dis .lt. eps) then
         find_cl = 1
         if (abs(LS_stencil(i,j,k,im_primary)) .lt. phimin) then
          phimin = abs(LS_stencil(i,j,k,im_primary))
          icl = i
          jcl = j
          kcl = k
          !print *, "phimin = ", phimin
          !print *, "icl, jcl, kcl = ", icl, jcl, kcl
         endif
        endif
       enddo
       enddo
       enddo
       if (find_cl .ne. 1) then
          print *, "Cannot find contact line!"
          return
       endif
       if (icl .eq. -n_rad .or. icl .eq. n_rad .or. &
           jcl .eq. -n_rad .or. jcl .eq. n_rad .or. &
           (prob_dim .eq. 3 .and. &
           (kcl .eq. -n_rad .or. kcl .eq. n_rad))) then
          print *, "Too far from contact line!"
          return
       endif
!calculate nphi and actual angle

        ! used to be LS1
       dir=1
       tmp(dir) = (LS_stencil(icl+1,jcl,kcl,im_primary)- &
                   LS_stencil(-1+icl,jcl,kcl,im_primary))/ &
         (x_stencil(icl+1,jcl,kcl,dir)-x_stencil(-1+icl,jcl,kcl,dir))
       dir=2
       tmp(dir) = (LS_stencil(icl,jcl+1,kcl,im_primary)- &
                   LS_stencil(icl,-1+jcl,kcl,im_primary))/ &
         (x_stencil(icl,jcl+1,kcl,dir)-x_stencil(icl,-1+jcl,kcl,dir))
       if (prob_dim .eq. 3) then
        dir=prob_dim
        tmp(dir) = (LS_stencil(icl,jcl,kcl+1,im_primary)- &
                    LS_stencil(icl,jcl,-1+kcl,im_primary))/ &
         (x_stencil(icl,jcl,kcl+1,dir)-x_stencil(icl,jcl,-1+kcl,dir))
       else
          tmp(3) = 0.d0
       endif


       mag = 0.d0
       do d = 1, prob_dim
          mag = mag+tmp(d)*tmp(d)
       enddo
       mag = sqrt(mag)

       costheta = 0.d0
       do d = 1, prob_dim
          nphi(d) = tmp(d)/mag
          costheta = costheta + nphi(d) * npsi(d)
          !print *, "d, nphi(d), npsi(d) = ", d, nphi(d), npsi(d)
       enddo
       actual_angle = acos(costheta)
!get point x_inf_proj
       do d = 1, prob_dim
           ! used to be LS1
          x_inf_proj(d) = x_stencil(icl,jcl,kcl,d)- &
                  LS_stencil(icl,jcl,kcl,im_primary)*nphi(d)
          !print *, "d = ", d, "x_inf_proj = ", x_inf_proj(d)
       enddo
!get point x_inf_proj_alpha
       if (prob_dim .eq. 2) then
          do d = 1, prob_dim
             x_inf_proj_alpha(d) = x_inf_proj(d)
          enddo
       else if (prob_dim .eq. 3) then
          dis = 0.d0
          do d = 1, prob_dim
             tmp(d) = x_inf_proj(d) - x_proj(d)
             dis = dis + tmp(d) * nalpha(d)
          enddo
          !if x_inf_proj locates in left side of plane_alpha, dis < 0
          do d = 1, prob_dim
             x_inf_proj_alpha(d) = x_inf_proj(d)-dis*nalpha(d)
             !print *, "d = ", d, "x_inf_proj_alpha = ", x_inf_proj_alpha(d)
          enddo
       else
          print *, "prob_dim invalid"
          stop
       endif
!get point x_psi_proj
       dis = 0.d0
       do d = 1, prob_dim
          tmp(d) = x_inf_proj_alpha(d) - x_proj(d)
          dis = dis + tmp(d) * npsi(d)
       enddo
       !if x_inf_proj_alpha locates in fluid zone, dis < 0
       do d = 1, prob_dim
          x_psi_proj(d) = x_inf_proj_alpha(d)-dis*npsi(d)
          !print *, "d = ", d, "x_psi_proj = ", x_psi_proj(d)
       enddo
!get point x_contact_point
       if (costheta .eq. 1.d0) then
          print *, "contact angle equal zero"
          stop
       else if (costheta .eq. 0.d0) then
          do d = 1, prob_dim
             x_contact_point(d) = x_psi_proj(d)
          enddo
       else
          do d = 1, prob_dim
             x_contact_point(d) = x_psi_proj(d)-dis/tan(actual_angle) &
                                  *tpsi(d)
             !print *, "d = ", d, "x_contact_point = ", x_contact_point(d)
          enddo
       endif
       closest_distance = 0.d0
       do d = 1, prob_dim
          closest_distance = closest_distance+(x_proj(d)- &
                             x_contact_point(d))**2.d0
       enddo
       closest_distance = sqrt(closest_distance)

    else
       print *, "i_method invalid"
       stop
    endif

  end subroutine closest_distance_to_CL

      subroutine get_vortex_info(x,time, &
            neg_force,vel,vort,energy_moment,temperature)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real) :: xprime(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(out) :: neg_force(SDIM)
      real(amrex_real), INTENT(out) :: vel(SDIM)
      real(amrex_real), INTENT(out) :: vort
      real(amrex_real), INTENT(out) :: energy_moment
      real(amrex_real), INTENT(out) :: temperature
      real(amrex_real) problo(SDIM)
      real(amrex_real) probhi(SDIM)
      real(amrex_real) problen(SDIM)
      real(amrex_real) rr,alpha,alpha_t,alpha_r,r_x,r_y,uy,vx
      real(amrex_real) fx,gy,fy,gx,fxp,gyp,fyp,gxp
      integer dir,dir_x,dir_y
      real(amrex_real) decay_term

      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif

      temperature=zero

      if ((probtype.eq.26).and. &
          ((axis_dir.eq.2).or. &  ! vortex confinement, no interface
           (axis_dir.eq.3))) then ! vortex confinement, interface

       if (SDIM.eq.2) then

        if ((adv_dir.eq.1).or. &
            (adv_dir.eq.2).or. &
            (adv_dir.eq.3)) then

         do dir=1,SDIM
          xprime(dir)=x(dir)
         enddo
         if ((adv_dir.eq.1).or.(adv_dir.eq.3)) then
          xprime(1)=xprime(1)-adv_vel*time
         endif
         if ((adv_dir.eq.2).or.(adv_dir.eq.3)) then
          xprime(2)=xprime(2)-adv_vel*time
         endif
         problo(1)=problox
         problo(2)=probloy
         probhi(1)=probhix
         probhi(2)=probhiy

         do dir=1,SDIM
          problen(dir)=probhi(dir)-problo(dir)
          if (problen(dir).gt.zero) then
           do while (xprime(dir).lt.problo(dir))
            xprime(dir)=xprime(dir)+problen(dir)
           enddo
           do while (xprime(dir).gt.probhi(dir))
            xprime(dir)=xprime(dir)-problen(dir)
           enddo
          else
           print *,"problen invalid"
           stop
          endif
         enddo ! dir=1..sdim

         if ((radblob2.lt.zero).or. &
             (radblob3.le.zero)) then
          print *,"radblob2 or radblob3 invalid"
          stop
         endif

         rr=sqrt((xprime(1)-xblob)**2+(xprime(2)-yblob)**2)-radblob
         energy_moment=rr+radblob

         if (rr+radblob.gt.1.0E-14) then
          alpha=half*vinletgas*(one-tanh(30.0*rr))* &
                (one+radblob2*sin(two*Pi*time/radblob3))
          vel(1)=alpha*(xprime(2)-yblob)       
          vel(2)=-alpha*(xprime(1)-xblob)       
          if ((adv_dir.eq.1).or.(adv_dir.eq.3)) then
           vel(1)=vel(1)+adv_vel
          endif
          if ((adv_dir.eq.2).or.(adv_dir.eq.3)) then
           vel(2)=vel(2)+adv_vel
          endif
          alpha_t=half*vinletgas*(one-tanh(30.0*rr))* &    
                radblob2*cos(two*Pi*time/radblob3)* &
                two*Pi/radblob3
          neg_force(1)=-alpha_t*(xprime(2)-yblob)
          neg_force(2)=alpha_t*(xprime(1)-xblob)       
          alpha_r=half*vinletgas*(tanh(30.0*rr)**2-one)* &    
                (one+radblob2*sin(two*Pi*time/radblob3))*30.0
          r_x=(xprime(1)-xblob)/(rr+radblob)
          r_y=(xprime(2)-yblob)/(rr+radblob)
          uy=alpha_r*r_y*(xprime(2)-yblob)+alpha
          vx=-alpha_r*r_x*(xprime(1)-xblob)-alpha
          vort=abs(uy-vx)
         else if (abs(rr+radblob).le.1.0E-14) then
          vel(1)=zero
          vel(2)=zero
          neg_force(1)=zero
          neg_force(2)=zero
          vort=two*vinletgas*(one+radblob2*sin(two*Pi*time/radblob3))
         else
          print *,"rr invalid"
          stop
         endif

        else
         print *,"adv_dir invalid probtype==26"
         stop
        endif

       else
        print *,"dimension other than 2 is not supported here"
        stop
       endif

      else if ((probtype.eq.26).and. &
               (axis_dir.eq.11)) then ! BCG smooth test, periodic BC.

       do dir=1,SDIM
        xprime(dir)=x(dir)
       enddo
       problo(1)=problox
       problo(2)=probloy
       probhi(1)=probhix
       probhi(2)=probhiy
       if (SDIM.eq.3) then
        problo(SDIM)=probloz
        probhi(SDIM)=probhiz
       endif

       if (SDIM.eq.2) then

        if ((adv_dir.ge.1).and.(adv_dir.le.3)) then
         ! do nothing
        else
         print *,"adv_dir invalid probtype==26"
         stop
        endif

        if ((adv_dir.eq.1).or.(adv_dir.eq.3)) then
         xprime(1)=xprime(1)-adv_vel*time
        endif
        if ((adv_dir.eq.2).or.(adv_dir.eq.3)) then
         xprime(2)=xprime(2)-adv_vel*time
        endif

       else if (SDIM.eq.3) then

        if ((adv_dir.ge.1).and.(adv_dir.le.7)) then
         ! do nothing
        else
         print *,"adv_dir invalid"
         stop
        endif

        if ((adv_dir.eq.1).or.(adv_dir.eq.4).or. &
            (adv_dir.eq.5).or.(adv_dir.eq.7)) then
         xprime(1)=xprime(1)-adv_vel*time
        endif
        if ((adv_dir.eq.2).or.(adv_dir.eq.4).or. &
            (adv_dir.eq.6).or.(adv_dir.eq.7)) then
         xprime(2)=xprime(2)-adv_vel*time
        endif
        if ((adv_dir.eq.3).or.(adv_dir.eq.5).or. &
            (adv_dir.eq.6).or.(adv_dir.eq.7)) then
         xprime(SDIM)=xprime(SDIM)-adv_vel*time
        endif

       else
        print *,"dimension bust"
        stop
       endif

       do dir=1,SDIM
        problen(dir)=probhi(dir)-problo(dir)
        if (problen(dir).gt.zero) then
         do while (xprime(dir).lt.problo(dir))
          xprime(dir)=xprime(dir)+problen(dir)
         enddo
         do while (xprime(dir).gt.probhi(dir))
          xprime(dir)=xprime(dir)-problen(dir)
         enddo
        else
         print *,"problen invalid"
         stop
        endif
       enddo ! dir=1..sdim

       energy_moment=zero
       decay_term=exp(-fort_viscconst(1)*eight*Pi*Pi*time)
       alpha=cos(two*Pi*radblob*time)*decay_term

       if ((radblob.eq.zero).and.(fort_viscconst(1).eq.zero)) then
        ! do nothing, ok set of parameters
       else if ((radblob.gt.zero).and. &
                (fort_viscconst(1).eq.zero)) then
        ! do nothing, ok set of parameters
       else if ((radblob.eq.zero).and. &
                (fort_viscconst(1).gt.zero)) then
        ! do nothing, ok set of parameters
       else
        print *,"radblob and/or fort_viscconst(1) must be changed"
        stop
       endif

       if (SDIM.eq.2) then
        dir_x=1
        dir_y=2
        if ((probhix.eq.one).and.(probhiy.eq.one)) then
         ! do nothing
        else
         print *,"probhix or probhiy invalid"
         stop
        endif
       else if (SDIM.eq.3) then

        if ((probhix.eq.one).and. &
            (probhiy.eq.one).and. &
            (probhiz.eq.half)) then
         dir_x=1
         dir_y=2
        else if ((probhix.eq.one).and. &
                 (probhiy.eq.half).and. &
                 (probhiz.eq.one)) then
         dir_x=1
         dir_y=SDIM
        else if ((probhix.eq.half).and. &
                 (probhiy.eq.one).and. &
                 (probhiz.eq.one)) then
         dir_x=2
         dir_y=SDIM
        else
         print *,"probhi x,y, or z invalid"
         stop
        endif

       else
        print *,"dimension bust"
        stop
       endif

       do dir=1,SDIM
        vel(dir)=zero
        neg_force(dir)=zero
       enddo
       fx=sin(two*Pi*xprime(dir_x))
       gy=cos(two*Pi*xprime(dir_y))
       fy=sin(two*Pi*xprime(dir_y))
       gx=cos(two*Pi*xprime(dir_x))
       vel(dir_x)=-alpha*fx*gy
       vel(dir_y)=alpha*fy*gx
       temperature=alpha*fx*fy+fort_tempconst(1)
       if (SDIM.eq.2) then
        if ((adv_dir.eq.1).or.(adv_dir.eq.3)) then
         vel(1)=vel(1)+adv_vel
        endif
        if ((adv_dir.eq.2).or.(adv_dir.eq.3)) then
         vel(2)=vel(2)+adv_vel
        endif
       else if (SDIM.eq.3) then

        if ((adv_dir.eq.1).or.(adv_dir.eq.4).or. &
            (adv_dir.eq.5).or.(adv_dir.eq.7)) then
         vel(1)=vel(1)+adv_vel
        endif
        if ((adv_dir.eq.2).or.(adv_dir.eq.4).or. &
            (adv_dir.eq.6).or.(adv_dir.eq.7)) then
         vel(2)=vel(2)+adv_vel
        endif
        if ((adv_dir.eq.3).or.(adv_dir.eq.5).or. &
            (adv_dir.eq.6).or.(adv_dir.eq.7)) then
         vel(SDIM)=vel(SDIM)+adv_vel
        endif
       else
        print *,"dimension bust"
        stop
       endif

       alpha_t=-sin(two*Pi*radblob*time)*two*Pi*radblob*decay_term
         
       fxp=gx*two*Pi
       gyp=-fy*two*Pi 
       fyp=gy*two*Pi
       gxp=-fx*two*Pi 
      
       neg_force(dir_x)=alpha_t*fx*gy
       neg_force(dir_y)=-alpha_t*fy*gx

       uy=-alpha*fx*gyp
       vx=alpha*fy*gxp
       vort=abs(uy-vx)
      else if ((probtype.eq.26).and. &
               (axis_dir.eq.12)) then ! buoyancy test
       print *,"get_vortex_info should not have been called:buoyancy"
       stop
      else
       print *,"get_vortex_info should not have been called"
       stop
      endif

      end subroutine get_vortex_info

      subroutine get_vort_vel_error(time,x,vort,vel,temperature, &
          vort_err,vel_err,temperature_err, &
          vort_expect,vel_expect,temperature_expect, &
          energy_moment)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: vort
      real(amrex_real), INTENT(in) :: vel(SDIM)
      real(amrex_real), INTENT(in) :: temperature
      real(amrex_real), INTENT(out) :: vort_err
      real(amrex_real), INTENT(out) :: vel_err
      real(amrex_real), INTENT(out) :: temperature_err
      real(amrex_real), INTENT(out) :: vort_expect
      real(amrex_real), INTENT(out) :: temperature_expect
      real(amrex_real), INTENT(out) :: energy_moment
      real(amrex_real), INTENT(out) :: vel_expect(SDIM)
      real(amrex_real) :: neg_force(SDIM)
      integer dir

      do dir=1,SDIM
       vel_expect(dir)=zero
      enddo
      vort_expect=zero
      temperature_expect=zero
      vort_err=zero
      temperature_err=zero
      vel_err=zero
      energy_moment=zero

      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
       stop
      endif
      
      if (probtype.eq.26) then
       if ((axis_dir.eq.2).or. & ! vortex confinement, no interface.
           (axis_dir.eq.3)) then ! vortex confinement, interface.
        if ((adv_dir.eq.1).or. &
            (adv_dir.eq.2).or. &
            (adv_dir.eq.3)) then
         if (SDIM.eq.2) then
          call get_vortex_info(x,time,neg_force,vel_expect,vort_expect, &
           energy_moment,temperature_expect)
          vort_err=abs(vort_expect-vort)
          do dir=1,SDIM
           vel_err=vel_err+(vel_expect(dir)-vel(dir))**2
          enddo
          vel_err=sqrt(vel_err)
         else
          print *,"sdim invalid"
          stop
         endif
        else
         print *,"adv_dir invalid (2) probtype==26"
         stop
        endif
       else if ((axis_dir.eq.0).or.(axis_dir.eq.1)) then
        ! do nothing
       else if (axis_dir.eq.10) then
        ! do nothing
       else if (axis_dir.eq.11) then ! BCG periodic
        call get_vortex_info(x,time,neg_force,vel_expect,vort_expect, &
          energy_moment,temperature_expect)
        vort_err=abs(vort_expect-vort)
        temperature_err=abs(temperature_expect-temperature)
        do dir=1,SDIM
         vel_err=vel_err+(vel_expect(dir)-vel(dir))**2
        enddo
        vel_err=sqrt(vel_err)
       else if (axis_dir.eq.12) then !buoyancy test.
        ! do nothing
       else
        print *,"axis_dir invalid"
        stop
       endif
      endif ! probtype.eq.26

      return
      end subroutine get_vort_vel_error

       ! called by fort_estdt: determine maximum force due to buoyancy. 
       ! if denconst_interface==0.0 (default), then
       !  denjump_scale=(rhoA - rhoB)/max(rhoA,rhoB)
      subroutine get_max_denjump_scale( &
              denjump_scale, &
              denconst_interface)
      use global_utility_module

      IMPLICIT NONE

      integer im,im_opp
      integer iten
      real(amrex_real), INTENT(in) :: denconst_interface(num_interfaces)
      real(amrex_real), INTENT(out) :: denjump_scale
      real(amrex_real) denjump_scale_temp
      real(amrex_real) max_den_interface
      real(amrex_real) den_interface
      integer internal_wave_exists

      denjump_scale=zero

      do im=1,num_materials 
       if (is_rigid(im).eq.0) then
        do im_opp=im+1,num_materials
         if (is_rigid(im_opp).eq.0) then
          call get_iten(im,im_opp,iten)
          if ((iten.ge.1).and.(iten.le.num_interfaces)) then
           ! do nothing
          else
           print *,"iten invalid"
           stop
          endif
          max_den_interface=max(fort_denconst(im),fort_denconst(im_opp))
          if (max_den_interface.gt.zero) then
      
           den_interface=denconst_interface(iten)
           if (den_interface.eq.zero) then
            ! do nothing
           else if (den_interface.gt.zero) then 
            if (den_interface.gt.max_den_interface) then
             max_den_interface=den_interface
            else
             print *,"need den_interface.gt.max_den_interface"
             stop
            endif
           else
            print *,"den_interface invalid"
            stop
           endif

           denjump_scale_temp=abs(fort_denconst(im)-fort_denconst(im_opp))/ &
             max_den_interface

           if ((denjump_scale_temp.ge.zero).and. &
               (denjump_scale_temp.le.one)) then

            if (denjump_scale_temp.gt.denjump_scale) then
             denjump_scale=denjump_scale_temp
            endif

           else
            print *,"require denjump_scale_temp in [0,1]"
            stop
           endif

          else
           print *,"max_den_interface invalid: ",max_den_interface
           stop
          endif
         else if (is_rigid(im_opp).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid PROB.F90"
          stop
         endif
        enddo ! im_opp=im+1,num_materials
       else if (is_rigid(im).eq.1) then
        ! do nothing
       else
        print *,"is_rigid invalid PROB.F90"
        stop
       endif
      enddo ! im=1..num_materials

      call SUB_INTERNAL_GRAVITY_WAVE_FLAG(internal_wave_exists)

      if (internal_wave_exists.eq.1) then
       if (twall.ge.fort_tempconst(1)) then
        ! do nothing
       else
        print *,"twall invalid get_max_denjump_scale"
        stop
       endif
       max_den_interface=max(fort_denconst(1),fort_denconst(2))
       if (max_den_interface.gt.zero) then
        ! density(T) = density_base * (1+expansion_factor(T))
        ! |density(T_outer)-density(T_inner)|=
        ! |density_base*(expansion_factor(T_outer,T_base)-
        !                expansion_factor(T_inner,T_base))| = 
        ! |density_base*fort_DrhoDT(1)*(T_outer-T_inner)|=
        ! fort_denconst(1)*|expansion_factor(T_inner,T_outer)|
        ! im=1
        ! temperature=twall (outer wall temperature)
        ! temperature_base=fort_tempconst(1) (inner wall temperature)
        ! default: expansion_factor=denjump_scale temp=
        !   fort_DrhoDT(1)*(temperature-temperature_base)
        call SUB_UNITLESS_EXPANSION_FACTOR(1,twall, &
          fort_tempconst(1),denjump_scale_temp)
        denjump_scale_temp=abs(denjump_scale_temp)*fort_denconst(1)/ &
           max_den_interface
        if (denjump_scale_temp.gt.denjump_scale) then
         denjump_scale=denjump_scale_temp
        endif
       else
        print *,"max_den_interface invalid"
        stop
       endif
      else if (internal_wave_exists.eq.0) then
       ! do nothing
      else
       print *,"internal_wave_exits invalid"
       stop
      endif

      return
      end subroutine get_max_denjump_scale
   

       ! Du/Dt=-grad p/rho - omega cross (omega cross r)- 2 omega cross u +
       !       \vec{g} 
       ! omega z^hat cross r r^hat=omega r theta^hat
       ! omega z^hat cross (omega z^hat cross r)=
       ! omega z^hat cross omega r theta^hat = -omega^2 r r^hat
       !  
       ! p=dt( \vec{g}\cdot\vec{x} + (1/2)Omega^2 r^2 )
       ! force=grad p=dt( \vec{g} + Omega^2 r r^hat )
       !
       ! called from fort_init_potential. (NAVIERSTOKES_3D.F90)
       ! called from EOS_error_ind (PROB.F90)
       ! called from presBDRYCOND (PROB.F90)
       ! called from fort_initdata (PROB.F90)
      subroutine general_hydrostatic_pressure_density( &
        i,j,k,level, &
        angular_velocity, &!intent(in) general_hydrostatic_pressure_density
        centrifugal_force_factor, &!intent(in) "     "
        dt, &
        rho_hydrostatic, &
        pres_hydrostatic, &
        state_ptr)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: i,j,k,level
       !general_hydrostatic_pressure_density
      real(amrex_real), INTENT(in) :: angular_velocity 
      real(amrex_real), INTENT(in) :: centrifugal_force_factor
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(out) :: rho_hydrostatic
      real(amrex_real), INTENT(out) :: pres_hydrostatic
      real(amrex_real), INTENT(in),pointer :: state_ptr(D_DECL(:,:,:),:)
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xcell(SDIM)
      integer :: local_dir
      integer :: gravity_dir
      integer, PARAMETER :: from_boundary_hydrostatic=0

      call fort_derive_gravity_dir(gravity_vector,gravity_dir)

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt must be positive"
       stop
      endif

      if (gravity_vector(gravity_dir).eq.zero) then
       ! do nothing
      else if (gravity_vector(gravity_dir).ne.zero) then
       ! do nothing
      else
       print *,"gravity_vector is NaN"
       stop
      endif

      call gridsten_level(xsten,i,j,k,level,nhalf)
      do local_dir=1,SDIM
       xcell(local_dir)=xsten(0,local_dir)
      enddo 

       ! the force is grad p^hydrostatic/rho^hydrostatic
      rho_hydrostatic=fort_denconst(1) 
      if (rho_hydrostatic.gt.zero) then
       pres_hydrostatic=zero
       do local_dir=1,SDIM
        pres_hydrostatic=pres_hydrostatic+ &
          gravity_vector(local_dir)*rho_hydrostatic*xcell(local_dir)
       enddo
    
       if (angular_velocity.ge.zero) then
        ! do nothing
       else
        print *,"angular_velocity should be nonneg (counter-clockwise): ", &
           angular_velocity
        stop
       endif

       if ((centrifugal_force_factor.ge.zero).and. &
           (centrifugal_force_factor.le.one)) then
        ! do nothing
       else
        print *,"expecting 0<=centrifugal_force_factor<=1"
        stop
       endif

       if (levelrz.eq.COORDSYS_CARTESIAN) then

        if (centrifugal_force_factor.eq.zero) then
         ! do nothing
        else if ((centrifugal_force_factor.gt.zero).and. &
                 (centrifugal_force_factor.le.one)) then
         pres_hydrostatic=pres_hydrostatic+ &
           half*rho_hydrostatic*centrifugal_force_factor* &
           (angular_velocity**2)*(xcell(1)**2+xcell(2)**2)
        else
         print *,"expecting 0<=centrifugal_force_factor<=1"
         stop
        endif
     
       else if (levelrz.eq.COORDSYS_RZ) then

        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (angular_velocity.eq.zero) then
         ! do nothing
        else
         print *,"angular_velocity must be 0 for RZ"
         stop
        endif

       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then

        pres_hydrostatic=pres_hydrostatic+ &
           half*rho_hydrostatic*centrifugal_force_factor* &
           (angular_velocity**2)*(xcell(1)**2)

       else
        print *,"levelrz invalid general hydrostatic pressure density"
        print *,"levelrz=",levelrz
        stop
       endif

      else
       print *,"rho_hydrostatic invalid: ",rho_hydrostatic
       stop
      endif

      if (is_in_probtype_list().eq.1) then
       call SUB_correct_pres_rho_hydrostatic( &
        i,j,k,level, &
        angular_velocity, &
        centrifugal_force_factor, &
        dt, &
        rho_hydrostatic, &
        pres_hydrostatic, &
        state_ptr)
      else if (fort_material_type(1).eq.1) then !TAIT EOS
       call tait_hydrostatic_pressure_density(xcell, &
         rho_hydrostatic, &
         pres_hydrostatic, &
         from_boundary_hydrostatic)
      endif

      if (1.eq.0) then
       print *,"before: rho,pres,dt,global_pressure_scale ", &
         rho_hydrostatic,pres_hydrostatic,dt,global_pressure_scale
      endif
 
        ! dt multiplied by velocity scale.
      pres_hydrostatic=pres_hydrostatic*dt/global_pressure_scale

      return
      end subroutine general_hydrostatic_pressure_density

       ! called from presBDRYCOND when material_type(1)==13 TAIT EOS
      subroutine boundary_hydrostatic( &
        xpos,rho,pres)
      use global_utility_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: xpos(SDIM)
      real(amrex_real), INTENT(inout) :: rho
      real(amrex_real), INTENT(inout) :: pres
      integer, PARAMETER :: from_boundary_hydrostatic=1

       ! first material obeys TAIT EOS
      if (fort_material_type(1).eq.13) then

        if ((probtype.eq.36).or. &   ! bubble in liquid
            (probtype.eq.601).or. &  ! cooling disk
            (probtype.eq.602)) then  ! Rayleigh-Taylor (and checkerboard test)
         if ((probtype.eq.36).and. &
             (axis_dir.eq.10)) then
          print *,"cannot use hydrostatic pressure for heat pipe experiment"
          print *,"modify: tait_hydrostatic_pressure_density for microgravity"
          stop
         endif
          !calling from boundary_hydrostatic
          !boundary_hydrostatic is called from presBDRYCOND
         call tait_hydrostatic_pressure_density(xpos,rho,pres, &
                 from_boundary_hydrostatic)
        else
         print *,"expecting probtype=36,601, or 602"
         stop
        endif

      else
       print *,"expecting liquid to be compressible (Tait EOS)"
       stop
      endif
 
      return
      end subroutine boundary_hydrostatic

! err is initialized already with the slope error
! pres_in comes from the pressure computed from the compressible
! projection method, not the equation of state.
! err=max(err,1.0d0)
! this routine called from: fort_pressure_indicator
      subroutine EOS_error_ind( &
       pressure_error_flag, &
       xsten,nhalf,bfact, &
       vort, &
       pres_in, &
       temp_in, &
       err, &
       imattype, &
       vorterr, &
       pressure_cutoff, &
       temperature_cutoff)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: imattype,pressure_error_flag,nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: vort
      real(amrex_real), INTENT(in) :: pres_in(D_DECL(3,3,3))
      real(amrex_real), INTENT(in) :: temp_in(D_DECL(3,3,3))
      real(amrex_real), INTENT(inout) :: err
      real(amrex_real), INTENT(in) :: vorterr
      real(amrex_real), INTENT(in) :: pressure_cutoff
      real(amrex_real), INTENT(in) :: temperature_cutoff

      real(amrex_real) atmos_pres,atmos_den
      real(amrex_real) pres_scale
      real(amrex_real) pres_test
      real(amrex_real) xpos(SDIM)
      real(amrex_real) temp_variation
      real(amrex_real) pres_variation
      integer dir,side,ii,jj,kk
      integer, PARAMETER :: from_boundary_hydrostatic=0

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid eos error ind"
       stop
      endif
      do dir=1,SDIM
       xpos(dir)=xsten(0,dir)
      enddo
      if ((pressure_error_flag.ne.0).and. &
          (pressure_error_flag.ne.1)) then
       print *,"pressure_error_flag invalid: ",pressure_error_flag
       stop
      endif
      if (vorterr.eq.zero) then
       ! do nothing
      else if (vorterr.gt.zero) then
       if (vort*global_velocity_scale.gt.vorterr) then
        if (err.eq.zero) then
         err=one
        else if (err.eq.one) then
         ! do nothing
        else
         print *,"err invalid"
         stop
        endif
       endif
      else
       print *,"vorterr invalid"
       stop
      endif

      if (temperature_cutoff.eq.zero) then
       ! do nothing
      else if (temperature_cutoff.gt.zero) then

       temp_variation=zero

       do dir=1,SDIM
        do side=0,1
         ii=0
         jj=0
         kk=0
         if (dir.eq.1) then
          ii=2*side-1
         else if (dir.eq.2) then
          jj=2*side-1
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          kk=2*side-1
         else
          print *,"dir invalid eos error ind"
          stop
         endif

         temp_variation=temp_variation+ &
          abs(temp_in(D_DECL(2,2,2))- &
              temp_in(D_DECL(2+ii,2+jj,2+kk)))
        enddo
       enddo

       temp_variation=temp_variation/(two*SDIM)

       if (temp_variation.gt.temperature_cutoff) then

        if (err.eq.zero) then
         err=one
        else if (err.eq.one) then
         ! do nothing
        else
         print *,"err invalid"
         stop
        endif

       endif

      else
       print *,"temperature cutoff invalid"
       stop
      endif

      pres_scale=pres_in(D_DECL(2,2,2))*global_pressure_scale
      pres_variation=zero

      if (pressure_cutoff.eq.zero) then
       ! do nothing
      else if (pressure_cutoff.gt.zero) then

       do dir=1,SDIM
        do side=0,1
         ii=0
         jj=0
         kk=0
         if (dir.eq.1) then
          ii=2*side-1
         else if (dir.eq.2) then
          jj=2*side-1
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          kk=2*side-1
         else
          print *,"dir invalid eos error ind 2"
          stop
         endif

         pres_variation=pres_variation+  &
          abs(pres_in(D_DECL(2,2,2))- &
              pres_in(D_DECL(2+ii,2+jj,2+kk)))
        enddo ! side
       enddo ! dir
       pres_variation=pres_variation/(two*SDIM)
       pres_variation=pres_variation*global_pressure_scale

       if (pressure_error_flag.eq.1) then
        pres_test=pres_variation
       else if (pressure_error_flag.eq.0) then
        pres_test=abs(pres_scale)
       else
        print *,"pressure_error_flag invalid"
        stop
       endif

       if (imattype.eq.0) then
        print *,"pressure cutoff should be 0 for incomp materials"
        stop
       else if ((imattype.ge.1).and. &
                (imattype.le.MAX_NUM_EOS)) then

         ! pressure comes from the compressible projection method, not
         ! the equation of state.
        if ((pres_scale.le.zero).and.(1.eq.0)) then
         print *,"pressure underflow"
         print *,"imattype=",imattype
         print *,"pres_scale= ",pres_scale
         stop
        endif

        if (imattype.eq.1) then  ! tait EOS 
           !calling from EOS_error_ind
         call tait_hydrostatic_pressure_density(xpos,atmos_den,atmos_pres, &
                 from_boundary_hydrostatic)
        else if ((imattype.ge.2).and. &
                 (imattype.le.MAX_NUM_EOS)) then
         call general_hydrostatic_pressure(atmos_pres)
        else 
         print *,"imattype invalid EOS_error_ind"
         stop
        endif

        if (pres_test/atmos_pres.gt.pressure_cutoff) then

         if (err.eq.zero) then
          err=one
         else if (err.eq.one) then
          ! do nothing
         else
          print *,"err invalid: ",err
          stop
         endif

        endif

       else
        print *,"imattype invalid EOS_error_ind"
        stop
       endif

      else
       print *,"pressure cutoff invalid"
       stop
      endif

      return
      end subroutine EOS_error_ind


       ! ADIABATIC_EOS_FLAG==0 if pressure does depend on internal energy.
       ! ADIABATIC_EOS_FLAG==1 if pressure does NOT depend on internal energy.
       ! NOT USED ANYMORE, BUT SHOULD STILL BE MAINTAINED FOR REINFORCING
       ! WHAT EACH EOS ASSUMES.
      subroutine ADIABATIC_EOS_FLAG(imattype,flag)
      IMPLICIT NONE

      integer imattype,flag


      if (imattype.eq.999) then
       flag=1
      else if (imattype.eq.0) then
       flag=1
      else if (imattype.eq.1) then  ! EOS_tait
       flag=1
      else if (imattype.eq.7) then  ! EOS_tait_rho
       flag=1
      else if (imattype.eq.13) then  ! EOS_tait_rhohydro
       flag=1
      else if (imattype.eq.2) then  ! jwl (adiabatic)
       flag=1
      else if (imattype.eq.8) then  ! air adiabatic
       flag=1
      else if (imattype.eq.9) then  ! EOS_tait_rho3
       flag=1
      else if (imattype.eq.10) then  ! EOS_tait_rho2
       flag=1
      else if (imattype.eq.11) then  ! EOS_koren_rho1
       flag=1
      else if (imattype.eq.12) then  ! EOS_koren_rho2
       flag=1
      else if (imattype.eq.15) then  ! EOS_dodecane
       flag=0
      else if (imattype.eq.19) then  ! EOS_vacuum
       flag=1
      else if (imattype.eq.20) then  ! EOS_tait_vacuum
       flag=1
      else if ((imattype.gt.0).and.(imattype.le.MAX_NUM_EOS)) then
       flag=0
      else
       print *,"material type invalid"
       stop
      endif

      return
      end subroutine
      

        ! returns e/scale
      subroutine INTERNAL_ENTROPY_material(rho,entropy,internal_energy, &
        imattype,im)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(in) :: entropy
      real(amrex_real), INTENT(out) :: internal_energy


      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (entropy.le.zero) then
       print *,"entropy invalid"
       stop
      endif

      if (imattype.eq.999) then
       call INTERNAL_default(rho,entropy,internal_energy,imattype,im)
      else if (imattype.eq.0) then
       call INTERNAL_default(rho,entropy,internal_energy,imattype,im)
      else if (imattype.eq.1) then
       call INTERNAL_tait(rho,entropy,internal_energy)
      else if (imattype.eq.2) then  ! this is adiabatic JWL
       call INTERNAL_jwl(rho,entropy,internal_energy)
      else if (imattype.eq.3) then  ! this is non-adiabatic JWL
       call INTERNAL_ENTROPY_jwl(rho,entropy,internal_energy)
      else if (imattype.eq.4) then
       call INTERNAL_ENTROPY_SF6(rho,entropy,internal_energy)
      else if (imattype.eq.5) then
       call INTERNAL_ENTROPY_air(rho,entropy,internal_energy)
      else if (imattype.eq.14) then
       call INTERNAL_ENTROPY_air_rho2(rho,entropy,internal_energy)
      else if (imattype.eq.6) then
       print *,"define INTERNAL_ENTROPY_Marquina"
       stop
       call INTERNAL_Marquina(rho,entropy,internal_energy)
      else if (imattype.eq.7) then
       call INTERNAL_tait_rho(rho,entropy,internal_energy)
      else if (imattype.eq.13) then
       call INTERNAL_tait_rhohydro(rho,entropy,internal_energy)
      else if (imattype.eq.8) then
       call INTERNAL_airADIABAT(rho,entropy,internal_energy)
      else if (imattype.eq.9) then
       call INTERNAL_tait_rho3(rho,entropy,internal_energy)
      else if (imattype.eq.10) then
       call INTERNAL_tait_rho2(rho,entropy,internal_energy)
      else if (imattype.eq.11) then
       call INTERNAL_koren_rho1(rho,entropy,internal_energy)
      else if (imattype.eq.12) then
       call INTERNAL_koren_rho2(rho,entropy,internal_energy)
      else if (imattype.eq.15) then
       print *,"define INTERNAL_ENTROPY_dodecane"
       stop
       call INTERNAL_dodecane(rho,entropy,internal_energy)
      else if (imattype.eq.16) then
       call INTERNAL_SF6ADIABAT(rho,entropy,internal_energy)
      else if (imattype.eq.17) then
       print *,"define INTERNAL_ENTROPY_stiffened"
       stop
       call INTERNAL_stiffened(rho,entropy,internal_energy,im)
      else if (imattype.eq.19) then
       call INTERNAL_vacuum(rho,entropy,internal_energy)
      else if (imattype.eq.20) then
       call INTERNAL_tait_vacuum(rho,entropy,internal_energy)
      else
       print *,"imattype invalid INTERNAL_ENTROPY_material"
       stop
      endif

      internal_energy=internal_energy/global_pressure_scale

      return
      end subroutine INTERNAL_ENTROPY_material



      subroutine TEMPERATURE_ENTROPY_material(rho,massfrac_parm, &
        entropy,temperature, &
        imattype,im)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)
      real(amrex_real), INTENT(in) :: entropy
      real(amrex_real) :: internal_energy
      real(amrex_real), INTENT(out) :: temperature


      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (entropy.le.zero) then
       print *,"entropy invalid"
       stop
      endif

      if (imattype.eq.4) then
       call TEMPERATURE_ENTROPY_SF6(rho,entropy,temperature)
      else
       call INTERNAL_ENTROPY_material(rho,entropy,internal_energy, &
        imattype,im)
       call TEMPERATURE_material(rho,massfrac_parm, &
        temperature,internal_energy, &
        imattype,im)
      endif

      return
      end subroutine TEMPERATURE_ENTROPY_material


       ! extracts temperature from density and internal energy
       ! returns T(e*scale)
      subroutine TEMPERATURE_material(rho,massfrac_parm, &
        temperature,internal_energy_in, &
        imattype,im)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho,internal_energy_in
      real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)
      real(amrex_real) internal_energy
      real(amrex_real), INTENT(out) :: temperature
      integer :: ispec

      internal_energy=internal_energy_in*global_pressure_scale

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid71"
       stop
      endif
      if (rho.gt.zero) then
       ! do nothing
      else
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.gt.zero) then
       ! do nothing
      else
       print *,"internal energy invalid in temperature material"
       print *,"rho,energy,imat ",rho,internal_energy,imattype 
       stop
      endif
      do ispec=1,num_species_var
       if (massfrac_parm(ispec).ge.zero) then
        ! do nothing
       else
        print *,"massfrac_parm(ispec) invalid"
        stop
       endif
      enddo

      if (is_in_probtype_list().eq.1) then
       call SUB_TEMPERATURE(rho,massfrac_parm, &
         temperature,internal_energy, &
         imattype,im,num_species_var)
      else 
       call TEMPERATURE_material_CORE(rho,massfrac_parm, &
        temperature,internal_energy, &
        imattype,im)
      endif

      return
      end subroutine TEMPERATURE_material



       ! extracts entropy from density and internal energy
       ! returns E(e*scale)
      subroutine ENTROPY_material(rho,massfrac_parm, &
        entropy,internal_energy_in, &
        imattype,im)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)
      real(amrex_real), INTENT(out) :: entropy
      real(amrex_real) :: internal_energy
      real(amrex_real), INTENT(in) :: internal_energy_in


      internal_energy=internal_energy_in*global_pressure_scale

      if ((im.le.0).or.(im.gt.num_materials)) then
       print *,"im invalid72"
       stop
      endif
      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy invalid in entropy material"
       print *,"rho,energy,imat ",rho,internal_energy,imattype 
       stop
      endif

      if (imattype.eq.999) then 
       call TEMPERATURE_default(rho,entropy,internal_energy, &
         imattype,im)
      else if (imattype.eq.0) then
       call TEMPERATURE_default(rho,entropy,internal_energy, &
         imattype,im)
      else if (imattype.eq.1) then
       call TEMPERATURE_tait(rho,entropy,internal_energy)
      else if (imattype.eq.2) then  ! adiabatic JWL
       call TEMPERATURE_jwl(rho,entropy,internal_energy)
      else if (imattype.eq.3) then  ! non-adiabatic JWL
       call ENTROPY_jwl(rho,internal_energy,entropy)
      else if (imattype.eq.4) then
       call ENTROPY_SF6(rho,internal_energy,entropy)
      else if (imattype.eq.5) then
       call ENTROPY_air(rho,internal_energy,entropy)
      else if (imattype.eq.14) then
       call ENTROPY_air_rho2(rho,internal_energy,entropy)
      else if (imattype.eq.6) then
       print *,"define ENTROPY_Marquina"
       stop
       call TEMPERATURE_Marquina(rho,entropy,internal_energy)
      else if (imattype.eq.7) then
       call TEMPERATURE_tait_rho(rho,entropy,internal_energy)
      else if (imattype.eq.13) then
       call TEMPERATURE_tait_rhohydro(rho,entropy,internal_energy)
      else if (imattype.eq.8) then
       call TEMPERATURE_airADIABAT(rho,entropy,internal_energy)
      else if (imattype.eq.9) then
       call TEMPERATURE_tait_rho3(rho,entropy,internal_energy)
      else if (imattype.eq.10) then
       call TEMPERATURE_tait_rho2(rho,entropy,internal_energy)
      else if (imattype.eq.11) then
       call TEMPERATURE_koren_rho1(rho,entropy,internal_energy)
      else if (imattype.eq.12) then
       call TEMPERATURE_koren_rho2(rho,entropy,internal_energy)
      else if (imattype.eq.15) then
       print *,"define ENTROPY_dodecane"
       stop
       call TEMPERATURE_dodecane(rho,entropy,internal_energy)
      else if (imattype.eq.16) then
       call TEMPERATURE_SF6ADIABAT(rho,entropy,internal_energy)
      else if (imattype.eq.17) then
       print *,"define ENTROPY_stiffened"
       stop
       call TEMPERATURE_stiffened(rho,entropy,internal_energy,im)
      else if (imattype.eq.19) then
       call TEMPERATURE_vacuum(rho,entropy,internal_energy)
      else if (imattype.eq.20) then
       call TEMPERATURE_tait_vacuum(rho,entropy,internal_energy)
      else
       print *,"imattype invalid ENTROPY_material"
       print *,"imattype= ",imattype
       stop
      endif

      return
      end subroutine ENTROPY_material



      subroutine ENTROPY_TEMPERATURE_material(rho,massfrac_parm, &
        entropy,temperature_in, &
        imattype,im)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: imattype,im
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(out) :: entropy
      real(amrex_real) :: internal_energy
      real(amrex_real), INTENT(in) :: temperature_in
      real(amrex_real), INTENT(in) :: massfrac_parm(num_species_var+1)

      if (imattype.eq.4) then
       call ENTROPY_TEMPERATURE_SF6(rho,temperature_in,entropy)
      else
       call INTERNAL_material(rho,massfrac_parm, &
         temperature_in,internal_energy, &
         imattype,im)
       call ENTROPY_material(rho,massfrac_parm, &
        entropy,internal_energy, &
        imattype,im)
      endif

      return
      end subroutine ENTROPY_TEMPERATURE_material

      function CLS(phi,eps)
      IMPLICIT NONE
      real(amrex_real) CLS,phi,eps,temp

      temp=phi/(two*eps)
      CLS=half*(sinh(temp)/cosh(temp)+one)

      return
      end function CLS

      function myfact(n)
      integer myfact,i,n

      if (n.eq.0) then
       myfact=1
      else
       myfact=1
       do i=2,n
        myfact=myfact*i 
       enddo
      endif

      return
      end function myfact


      subroutine get_pipe_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)
      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer bfact,nhalf
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      integer, parameter :: nhalf2=1
      real(amrex_real) xsten2(-nhalf2:nhalf2,SDIM)

      real(amrex_real) dx(SDIM)
      real(amrex_real) cenbc(num_materials,SDIM)
      real(amrex_real) vfrac(num_materials)
      integer im

      integer dir2,i1,j1,k1,k1lo,k1hi
      real(amrex_real) centroid(num_materials,SDIM)
      real(amrex_real) lsgrid(D_DECL(3,3,3),num_materials)
      real(amrex_real) distbatch(num_materials)
      real(amrex_real) facearea(num_materials)
      real(amrex_real) pipexlo,pipexhi,vfrac_sum
      real(amrex_real) EBVOFTOL
      integer isten

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid get pipe vfrac"
       stop
      endif

      if ((probtype.eq.41).and.(axis_dir.eq.4)) then
       if (rigid_count().ne.0) then
        print *,"solid not allowed for comparison with LSA"
        stop
       endif
      else
       if (rigid_count().eq.0) then
        print *,"solid missing for pipe"
        stop
       endif
      endif

      if (probtype.eq.41) then
       ! do nothing
      else
       print *,"probtype invalid in get pipe vfrac"
       stop
      endif

      pipexlo=problox
      pipexhi=probhix
      if ((axis_dir.eq.1).or. &
          (axis_dir.eq.2).or. &
          (axis_dir.eq.3)) then
       pipexlo=zero
       pipexhi=two*radblob3
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

      do k1=k1lo,k1hi
      do j1=-1,1
      do i1=-1,1
       do isten=-1,1
        dir2=1
        xsten2(isten,dir2)=xsten(isten+2*i1,dir2)
        dir2=2
        xsten2(isten,dir2)=xsten(isten+2*j1,dir2)
        if (SDIM.eq.3) then
         dir2=SDIM
         xsten2(isten,dir2)=xsten(isten+2*k1,dir2)
        endif
       enddo ! isten

       call inletpipedist(xsten2(0,1),xsten2(0,2),xsten2(0,SDIM),distbatch)
       do im=1,num_materials
        lsgrid(D_DECL(i1+2,j1+2,k1+2),im)=distbatch(im) 
       enddo
      enddo
      enddo
      enddo
      EBVOFTOL=VOFTOL
      call getvolumebatch(bfact,dx,xsten,nhalf, &
        lsgrid,vfrac, &
        facearea,centroid,EBVOFTOL,SDIM)
      do im=1,num_materials
       do dir2=1,SDIM
        cenbc(im,dir2)=centroid(im,dir2)-xsten(0,dir2)
       enddo
      enddo

      if (rigid_exists().eq.1) then

       if ((axis_dir.ge.0).and.(axis_dir.le.4)) then

        if ((xsten(0,1).lt.pipexlo).or.(xsten(0,1).gt.pipexhi)) then
         do im=1,num_materials
          vfrac(im)=zero
          do dir2=1,SDIM
           cenbc(im,dir2)=zero
          enddo
          if (is_rigid(im).eq.1) then
           vfrac(im)=one
          endif
         enddo !im=1..nat
        endif
       else if (axis_dir.eq.5) then
        ! do nothing
       else
        print *,"axis_dir invalid"
        stop
       endif
      else if (rigid_exists().eq.0) then
       ! do nothing
      else
       print *,"rigid exists bust"
       stop
      endif

       ! kluge
      vfrac_sum=zero
      do im=1,num_materials
       if (is_rigid(im).eq.0) then
        vfrac_sum=vfrac_sum+vfrac(im)
       else if (is_rigid(im).eq.1) then
        ! do nothing
       else
        print *,"is_rigid(im) invalid"
        stop
       endif
      enddo ! im=1..num_materials
      if (vfrac_sum.lt.VOFTOL) then
       vfrac(1)=one
       vfrac(2)=zero
      endif

      return
      end subroutine get_pipe_vfrac

         ! nslope points into the solid
         ! materialdistsolid: dist>0 in solid
      subroutine find_LS_stencil_slope( &
       bfact,dx,xsten0,nhalf0, &
       nslope,time,im)
      use global_utility_module
      use global_distance_module

      IMPLICIT NONE

      integer im
      integer bfact,nhalf0
      real(amrex_real) xsten0(-nhalf0:nhalf0,SDIM)
      real(amrex_real) dx(SDIM)
      real(amrex_real) nslope(SDIM)
      real(amrex_real) xn,yn,zn,time,mag
      integer dir
      integer inode,jnode,knode
      integer idx_array(3)
      integer knodelo,knodehi
      real(amrex_real) dist

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid73"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in find_ls_stencil_slope"
      else if (time.lt.zero) then
       print *,"time invalid in find_ls_stencil_slope"
       stop
      else
       print *,"time bust in find_ls_stencil_slope"
       stop
      endif

      if (is_rigid(im).eq.1) then

       do dir=1,SDIM
        nslope(dir)=zero
       enddo

       if (SDIM.eq.2) then
        knodelo=0
        knodehi=0
       else if (SDIM.eq.3) then
        knodelo=-1
        knodehi=1
       else
        print *,"dimension bust"
        stop
       endif

       do knode=knodelo,knodehi,2
       do jnode=-1,1,2
       do inode=-1,1,2

        idx_array(1)=inode
        idx_array(2)=jnode
        idx_array(3)=knode

        xn=xsten0(inode,1)
        yn=xsten0(jnode,2)
        if (SDIM.eq.2) then
         zn=yn
        else if (SDIM.eq.3) then
         zn=xsten0(knode,SDIM)
        else
         print *,"dimension bust prior to materialdistsolid"
         stop
        endif
        call materialdistsolid(xn,yn,zn,dist,time,im)

        do dir=1,SDIM
         nslope(dir)=nslope(dir)+idx_array(dir)*dist
        enddo

       enddo
       enddo
       enddo ! inode,jnode,knode

       mag=zero
       do dir=1,SDIM
        nslope(dir)=nslope(dir)/ &
          (xsten0(1,dir)-xsten0(-1,dir))
        mag=mag+nslope(dir)**2
       enddo
       mag=sqrt(mag)
       if (mag.eq.zero) then
        do dir=1,SDIM
         nslope(dir)=zero
        enddo
        nslope(SDIM)=one
       else if (mag.gt.zero) then
        do dir=1,SDIM
         nslope(dir)=nslope(dir)/mag
        enddo
       else
        print *,"mag invalid find_LS_stencil_slope"
        stop
       endif

      else
       print *,"expecting is_rigid(im)==1"
       stop
      endif

      return
      end subroutine find_LS_stencil_slope

       ! cen in absolute coordinates
       ! returns a volume fraction
      subroutine find_LS_stencil_volume_coarse( &
       bfact,dx,xsten,nhalf, &
       time,vfrac,cen, &
       volbox,cenbox,im)
      use global_utility_module
      use global_distance_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer im,bfact,nhalf
      real(amrex_real) time
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) dx(SDIM)
      real(amrex_real) vfrac
      real(amrex_real) cen(SDIM)
      integer dir
      real(amrex_real) lnode(4*(SDIM-1))
      integer inode,jnode,knode,klo,khi
      integer isynth
      real(amrex_real) xn,yn,zn,facearea
      real(amrex_real) volbox
      real(amrex_real) cenbox(SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid find override vfrac coarse"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in find_ls_stencil_volume_coarse"
      else if (time.lt.zero) then
       print *,"time invalid in find_ls_stencil_volume_coarse"
       stop
      else
       print *,"time bust in find_ls_stencil_volume_coarse"
       stop
      endif

      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid74"
       stop
      endif

      if (is_rigid(im).eq.1) then

       isynth=1
       if (SDIM.eq.2) then
        klo=0
        khi=0
       else if (SDIM.eq.3) then
        klo=-1
        khi=1
       else
        print *,"SDIM invalid"
        stop
       endif

       do knode=klo,khi,2
       do jnode=-1,1,2
       do inode=-1,1,2
        xn=xsten(inode,1)
        yn=xsten(jnode,2)
        if (SDIM.eq.2) then
         zn=yn
        else if (SDIM.eq.3) then
         zn=xsten(knode,SDIM)
        else
         print *,"dimension bust calling materialdistsolid"
         stop
        endif
        call materialdistsolid(xn,yn,zn,lnode(isynth),time,im)
        isynth=isynth+1
       enddo
       enddo
       enddo

        ! returns centroid in absolute coordinate system
       call fast_cell_intersection_grid( &
        bfact,dx,xsten,nhalf, &
        lnode, &
        vfrac,cen,facearea, &
        volbox,cenbox,SDIM)
       if (volbox.le.zero) then
        vfrac=zero
       else
        vfrac=vfrac/volbox
       endif

       if (vfrac.le.VOFTOL) then
        vfrac=zero
        do dir=1,SDIM
         cen(dir)=cenbox(dir)
        enddo
       endif
       if (vfrac.ge.one-VOFTOL) then
        vfrac=one
        do dir=1,SDIM
         cen(dir)=cenbox(dir)
        enddo
       endif

      else
       print *,"is_rigid invalid PROB.F90"
       stop
      endif

      return
      end subroutine find_LS_stencil_volume_coarse

       ! cen in absolute coordinates
      subroutine find_LS_stencil_volume( &
       bfact, &
       dx, &
       xsten,nhalf, &
       nrefine, &
       time, &
       vfrac, &
       cen, &
       im)
      use global_utility_module

      IMPLICIT NONE

      integer im
      integer nrefine,bfact,nhalf
      real(amrex_real) time
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer, parameter :: nhalf2=1
      real(amrex_real) xsten2(-nhalf2:nhalf2,SDIM)
      real(amrex_real) dx(SDIM)
      real(amrex_real) vfrac
      real(amrex_real) cen(SDIM)
      integer dir
      integer inode,jnode,knode,khi
      integer nside,iside
      real(amrex_real) dxrefine(SDIM)
      real(amrex_real) vfrac_node
      real(amrex_real) cen_node(SDIM)
      real(amrex_real) volbox,volbox_node
      real(amrex_real) cenbox(SDIM)
      real(amrex_real) cenbox_node(SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid find override vfrac"
       stop
      endif
      if ((nrefine.lt.0).or.(nrefine.gt.2)) then
       print *,"nrefine out of range"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid75"
       stop
      endif

      if (is_rigid(im).eq.1) then

       nside=1
       do iside=1,nrefine
        nside=nside*2
       enddo
       if (SDIM.eq.2) then
        khi=1
       else if (SDIM.eq.3) then
        khi=nside
       else
        print *,"SDIM invalid"
        stop
       endif

       do dir=1,SDIM
        cen(dir)=zero
        cenbox(dir)=zero
       enddo
       vfrac=zero
       volbox=zero

       do dir=1,SDIM
        dxrefine(dir)=(xsten(1,dir)-xsten(-1,dir))/nside
        if (dxrefine(dir).le.zero) then
         print *,"dxrefine invalid"
         stop
        endif
       enddo ! dir=1..sdim

       do knode=1,khi
       do jnode=1,nside
       do inode=1,nside
        dir=1
        xsten2(-1,dir)=xsten(-1,dir)+(inode-1)*dxrefine(dir)
        dir=2
        xsten2(-1,dir)=xsten(-1,dir)+(jnode-1)*dxrefine(dir)
        if (SDIM.eq.3) then
         dir=SDIM
         xsten2(-1,dir)=xsten(-1,dir)+(knode-1)*dxrefine(dir)
        endif
        do dir=1,SDIM
         xsten2(1,dir)=xsten2(-1,dir)+dxrefine(dir)
         xsten2(0,dir)=(xsten2(-1,dir)+xsten2(1,dir))/two
        enddo

        call find_LS_stencil_volume_coarse( &
         bfact,dx,xsten2,nhalf2, &
         time,vfrac_node,cen_node, &
         volbox_node,cenbox_node,im)

        vfrac=vfrac+vfrac_node*volbox_node
        volbox=volbox+volbox_node
        do dir=1,SDIM
         cenbox(dir)=cenbox(dir)+cenbox_node(dir)*volbox_node
         cen(dir)=cen(dir)+cen_node(dir)*vfrac_node*volbox_node
        enddo 
       enddo
       enddo
       enddo

       if (volbox.le.zero) then
        vfrac=zero
        do dir=1,SDIM
         cenbox(dir)=xsten(0,dir)
        enddo
       else
        vfrac=vfrac/volbox
        do dir=1,SDIM
         cenbox(dir)=cenbox(dir)/volbox
         if (vfrac.ge.VOFTOL) then
          cen(dir)=cen(dir)/(vfrac*volbox)
         endif
        enddo ! dir
       endif

       if (vfrac.le.VOFTOL) then
        vfrac=zero
        do dir=1,SDIM
         cen(dir)=cenbox(dir)
        enddo
       endif
       if (vfrac.ge.one-VOFTOL) then
        vfrac=one
        do dir=1,SDIM
         cen(dir)=cenbox(dir)
        enddo
       endif

      else
       print *,"is_rigid invalid PROB.F90"
       stop
      endif

      return
      end subroutine find_LS_stencil_volume


      subroutine inletpipedist(x,y,z,dist)
      use global_utility_module
      use global_distance_module

      IMPLICIT NONE

      integer im
      real(amrex_real), INTENT(out) :: dist(num_materials)
      real(amrex_real), INTENT(in) :: x,y,z
      real(amrex_real) ht,rr,initial_time
      integer im_solid_pipe

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y in inletpipedist"
        stop
       endif
      endif

      im_solid_pipe=im_solid_primary()

      initial_time=zero

      if (probtype.eq.41) then
       ! do nothing
      else
       print *,"probtype invalid for inlet pipedist"
       stop
      endif

      if ((axis_dir.lt.0).or.(axis_dir.gt.5)) then
       print *,"axis dir invalid for pipe problem"
       stop
      endif

      if (axis_dir.eq.0) then

       if (SDIM.eq.2) then
        dist(1)=xblob+radblob*cos(two*Pi*y/yblob)-x
       else if (SDIM.eq.3) then
        rr=sqrt(y**2+z**2)
        dist(1)=rr-zblob-radblob*sin(two*Pi*x/xblob)
       else
        print *,"dimension bust"
        stop
       endif

      else if ((axis_dir.eq.1).or. &
               (axis_dir.eq.2).or. &
               (axis_dir.eq.3)) then

       if (SDIM.eq.3) then
        dist(1)=z-zblob-radblob*sin(two*Pi*x/xblob) 
       else if (SDIM.eq.2) then
        dist(1)=abs(x)-xblob-radblob*sin(two*Pi*y/yblob)
        if (axis_dir.eq.2) then ! gas in the middle
         ht=radblob2+radblob*sin(two*Pi*y/yblob)
         if (x.gt.xblob) then
          dist(1)=x-(xblob+ht)
         else
          dist(1)=(xblob-ht)-x
         endif
        endif
        if (axis_dir.eq.3) then ! blowout problem
         dist(1)=y
        endif
       else
        print *,"dimension bust"
        stop
       endif

      else if (axis_dir.eq.4) then

       dist(1)=abs(x)-xblob-radblob*sin(two*Pi*y/yblob)
       dist(1)=-dist(1)

      else if (axis_dir.eq.5) then

       dist(1)=-99999.0

      else 
       print *,"axis dir invalid"
       stop
      endif

      dist(2)=-dist(1)
      do im=3,num_materials
       dist(im)=-99999.0
      enddo

        ! axis_dir=4 is the comparison with Linear Stability Analysis
      if (axis_dir.ne.4) then
       if ((im_solid_pipe.lt.1).or. &
           (im_solid_pipe.gt.num_materials)) then
        print *,"im_solid_pipe invalid 2"
        stop
       endif
         ! in inlet_pipe_dist; positive in solid
       call materialdistsolid(x,y,z,dist(im_solid_pipe),initial_time, &
        im_solid_pipe)
      endif

      return
      end subroutine inletpipedist

      subroutine get_pipe_velocity(xsten,nhalf,dx,bfact,vel,time)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) x,y,z,r
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) cenbc(num_materials,SDIM)
      real(amrex_real), INTENT(out) :: vel(SDIM)
      real(amrex_real) VOF(num_materials)
      integer dir2
      integer im_solid_pipe
      real(amrex_real) x_vel,y_vel,z_vel

      im_solid_pipe=im_solid_primary()

      if (nhalf.lt.3) then
       print *,"nhalf invalid get pipe velocity"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      if (probtype.eq.41) then
       ! do nothing
      else 
       print *,"probtype invalid in get pipe velocity"
       stop
      endif

      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y line 4073"
        stop
       endif
      endif

      if ((axis_dir.lt.0).or.(axis_dir.gt.5)) then
       print *,"get_pipe_velocity: axis dir invalid for pipe problem"
       stop
      endif

      do dir2=1,SDIM
       vel(dir2)=zero
      enddo

      call get_pipe_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc)  

        ! axis_dir=4 is LSA comparison
      if (axis_dir.eq.4) then

       if (VOF(1).gt.zero) then
        if (time.gt.zero) then
         vel(SDIM)=vinletgas
        else if (time.eq.zero) then
         vel(SDIM)=yblob4
        else
         print *,"time invalid in get pipe velocity"
         print *,"time= ",time
         print *,"vof1,vof2 ",VOF(1),VOF(2)
         print *,"axis_dir ",axis_dir
         stop
        endif
       else if ((VOF(1).eq.zero).and.(VOF(2).gt.zero)) then
        if (time.gt.zero) then
         vel(SDIM)=advbot
        else if (time.eq.zero) then
         vel(SDIM)=xblob4
        else
         print *,"time invalid in get pipe velocity"
         print *,"time= ",time
         print *,"vof1,vof2 ",VOF(1),VOF(2)
         print *,"axis_dir ",axis_dir
         stop
        endif
       endif
      else if (axis_dir.eq.3) then
       if ((im_solid_pipe.lt.1).or. &
           (im_solid_pipe.gt.num_materials)) then
        print *,"im_solid_pipe invalid 2.9"
        stop
       endif
       if (VOF(im_solid_pipe).gt.zero) then
        vel(SDIM)=zero
       else if (z.ge.VOFTOL*dx(SDIM)) then !y=z if 2D
        if (time.gt.zero) then
         vel(SDIM)=vinletgas
        else if (time.eq.zero) then
         vel(SDIM)=yblob4
        else
         print *,"time invalid in get pipe velocity"
         print *,"time= ",time
         print *,"vof1,vof2 ",VOF(1),VOF(2)
         print *,"axis_dir ",axis_dir
         stop
        endif
       else if (z.le.VOFTOL*dx(SDIM)) then ! y=z if 2D
        if (time.gt.zero) then
         vel(SDIM)=advbot
        else if (time.eq.zero) then
         vel(SDIM)=xblob4
        else
         print *,"time invalid in get pipe velocity"
         print *,"time= ",time
         print *,"vof1,vof2 ",VOF(1),VOF(2)
         print *,"axis_dir ",axis_dir
         stop
        endif
       endif

      else if ((axis_dir.eq.0).or. &
               (axis_dir.eq.1).or. &
               (axis_dir.eq.2)) then

       if ((im_solid_pipe.lt.1).or. &
           (im_solid_pipe.gt.num_materials)) then
        print *,"im_solid_pipe invalid 3"
        stop
       endif
       if (VOF(im_solid_pipe).gt.zero) then
        vel(SDIM)=zero
       else if (VOF(1).gt.zero) then
        if (time.gt.zero) then
         vel(SDIM)=vinletgas
        else if (time.eq.zero) then
         vel(SDIM)=yblob4
        else
         print *,"time invalid in get pipe velocity"
         print *,"time= ",time
         print *,"vof1,vof2 ",VOF(1),VOF(2)
         print *,"axis_dir ",axis_dir
         stop
        endif
       else if ((VOF(1).eq.zero).and.(VOF(2).gt.zero)) then
        if (time.gt.zero) then
         vel(SDIM)=advbot
        else if (time.eq.zero) then
         vel(SDIM)=xblob4
        else
         print *,"time invalid in get pipe velocity"
         print *,"time= ",time
         print *,"vof1,vof2 ",VOF(1),VOF(2)
         print *,"axis_dir ",axis_dir
         stop
        endif
       endif
 
       ! above: axis_dir=0,1,2 (and above that 3,4)

      else if (axis_dir.eq.5) then

       if (SDIM.eq.3) then
        ! do nothing
       else if (SDIM.eq.2) then
        z=zero
       else
        print *,"dimension bust"
        stop
       endif
       r = sqrt(y*y+z*z)
     
       if (r.gt.radblob) then
        x_vel=zero
       else 
        x_vel=adv_vel*r*(radblob-r)
       endif
       y_vel=zero
       z_vel=zero
       if((z-0.3)**2+(y-0.3)**2+(x-0.5)**2.le.0.4)then
         x_vel = .3d0
         y_vel = .4d0
         z_vel = -.2d0
       else if((z+0.3)**2+(y-1.7)**2+(x-5.5)**2.le.0.4)then
         x_vel = -.3d0
         y_vel = -.4d0
         z_vel = .2d0
       else if((z-0.9)**2+(y)**2+(x-1.5)**2.le.1.)then
         x_vel = -x_vel
         y_vel = -y_vel
         z_vel = -z_vel 
       else if((z+0.9)**2+(y)**2+(x-3.5)**2.le.1.)then
         x_vel = -y_vel
         y_vel = -x_vel
         z_vel = -z_vel 
       else if((z)**2+(y)**2+(x-2.5)**2.le.0.1)then
         x_vel = -z_vel
         y_vel = -y_vel
         z_vel = -x_vel 
       endif
       vel(1)=x_vel
       vel(2)=y_vel
       if (SDIM.eq.3) then
        vel(SDIM)=z_vel
       endif
      else
       print *,"axis_dir invalid"
       stop
      endif

      return
      end subroutine get_pipe_velocity


! time is in microseconds, LL is in microns,PTERM is in atmospheres
! 1atm=1.013x10^6 dyne/cm^2

      subroutine microfabpressure(LL,PTERM,time)

      use global_distance_module

      IMPLICIT NONE
      real(amrex_real) LL,PTERM,time
      real(amrex_real) realtime,realpress
      integer error


      real(amrex_real) NOD,NID,NPT,CHH,CHW,JLEN

      call microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)

      LL=half*(JLEN-NPT)
! simple pressure 2
! assume flat meniscus, so quiescent pressure is 1atm
      if (axis_dir.eq.13) then
       PTERM=zero
      else
       if (time.le.6.0) then
        PTERM=0.7-1.0
       else if (time.le.14.0) then
        PTERM=1.8-1.0
       else if (time.le.21.0) then
        PTERM=0.2-1.0
       else if (time.le.28.0) then
        PTERM=1.2-1.0
       else 
        PTERM=zero
       endif

       realtime=time*EPS6
       call pressure_bc(realpress,realtime,error)
       PTERM=realpress/1.00D+06 - one 
      endif
            
      return
      end subroutine

      subroutine microfabdist(rr,z,dist)

      use global_distance_module

      IMPLICIT NONE
      real(amrex_real) x,y,z,rr,dist


      real(amrex_real) NOD,NID,NPT,CHH,CHW,JLEN,incline

      call microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)

      if (SDIM.eq.2) then

      incline=JLEN 
      if (rr.le.half*NOD) then
       if (axis_dir.eq.12) then
        incline=JLEN+sqrt(16.9**2-(half*NOD)**2)-sqrt(16.9**2-rr**2)
       endif
       dist=incline-z
      else
       if (rr.le.half*(NPT+NOD)) then
        dist=incline-z-(rr-half*NOD)
       else
        dist=incline-z-half*NPT
       endif
      endif

      else if (SDIM.eq.3) then

      if (axis_dir.eq.13) then
       xblob=zero
       yblob=zero
      endif

      rr=sqrt( (x-xblob)**2 + (y-yblob)**2 )

      incline=JLEN
      if (rr.le.half*NOD) then
       if (axis_dir.eq.7) then
        incline=JLEN+sqrt(16.9**2-(half*NOD)**2)-sqrt(16.9**2-rr**2)
       endif
       dist=incline-z
      else
       if (rr.le.half*NID) then
        dist=incline-z-(rr-half*NOD)
       else
        dist=incline-z-half*(NID-NOD)
       endif
      endif

      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine
 
      subroutine pressure_bc( press, t    , error )
      IMPLICIT NONE

!     * Return the pressure boundary condition at time t provided by MicroFab
!     *  for the Okidata problems.

!     * CAUTION! Pressure is in dynes/ cm^2 and time is in seconds!
      
!     * This routine assumes that time is measured in seconds, the time values
!     * are equally spaced with spacing dt and that the pressure values are 
!     * known up to time t = 7.00D-05. For the MicroFab 
!     * test problems j= 0 ... = 70.  However,  to allow for future use of this
!     * type of BC, the size of the arrays time and pressbc may have to be 
!     * increased or passed as an argument to this subroutine.  Right now the
!     * number of data points in press_file is hard coded to be 70.

!     * Variables passed in ...

      integer   error  
      
      real(amrex_real)    press       , t
      
!     * The array pressbc(i,j).  The first array contains the time t, the
!     * second contains the value of the pressure on the inflow boundary at
!     * time t.

      real(amrex_real)  sigma    

!     integer  1234567890, 1234567890, 1234567890, 1234567890, 1234567890

      integer  j         
      
      error = 0
      j = int(t / dt_pressure_bcs)
      
      if (t .gt. 7.00D-05) then

!       * From 70 microseconds on, the pressure BC is 1 atm
        
        press = 1.00D+06
        
      else if ((time_pressure_bcs(j) .le. t) .and.  &
               (t .le. time_pressure_bcs(j+1))) then
        
!       * Lineraly interpolate between the given time values.
        
        sigma = (t - time_pressure_bcs(j)) /dt_pressure_bcs
        press = (1 - sigma) *  &
          pressbc_pressure_bcs(j,selectpress) +  &
          sigma * pressbc_pressure_bcs(j+1,selectpress) 
        
      else
       error = -1
       print *,"error in pressure_bc, t= ",t
       stop
      end if

      return
      end subroutine

!
!     * This is a subroutine for setting inflow B.C.(Poiseuille flow)
!     * Nozzle radius (cm)   
!     probtype.eq.25 and axis_dir.gt.0
!
      subroutine bubble_formation_inflow_bc(xsten,nhalf,x,phi)
      IMPLICIT NONE

      integer nhalf
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) delta
      real(amrex_real) aveQ,aveV,radius,diameter,x,phi
      real(amrex_real) density,  viscosity,  sigma
      real(amrex_real) Weber,Reynolds,Froude
      integer old_flag,Bright_flag

      if (nhalf.lt.1) then
       print *,"nhalf invalid bubble formation inflow bc"
       stop
      endif
      if (probtype.eq.25) then

       old_flag=0
       Bright_flag=1

       if (old_flag.eq.1) then
        radius=0.127
       else if (old_flag.eq.0) then
        if (Bright_flag.eq.0) then
         radius=0.05
        else if (Bright_flag.eq.1) then
         radius=0.08
        else
         print *,"Bright_flag invalid"
         stop
        endif
       else
        print *,"old_flag invalid"
        stop
       endif

       diameter=two*radius

       if (old_flag.eq.1) then

        if (axis_dir.eq.1) then
         aveQ = 1.3D-2
        else if (axis_dir.eq.2) then
         aveQ = 4.8D-2
        else if (axis_dir.eq.3) then
         aveQ = 2.0D-1
        else if (axis_dir.eq.4) then
         aveQ = 5.0D-1
        else if (axis_dir.eq.5) then
         aveQ = 1.1
        else if (axis_dir.eq.6) then
         aveQ = two
        else if (axis_dir.eq.7) then
         aveQ = five
        else if (axis_dir.eq.8) then
         aveQ = 6.8
        else if (axis_dir.eq.9) then
         aveQ = 7.2
        else if (axis_dir.eq.10) then
         aveQ = 15.0
        else if (axis_dir.eq.11) then
         aveQ = 20.0
        else
         print *,"axis_dir invalid for bubble formation"
         stop
        endif
        aveV=aveQ/(Pi*(radius**2))

        density=0.996
        sigma=51.1
        viscosity=0.00958

       else if (old_flag.eq.0) then

        if (axis_dir.eq.1) then
         aveV = 10.0D0
        else if (axis_dir.eq.2) then
         aveV = 20.0D0
        else if (axis_dir.eq.3) then
         aveV = 30.0D0
        else if (axis_dir.eq.4) then
         aveV = 40.0D0
        else if (axis_dir.eq.5) then
         aveV = 50.0D0
        else if (axis_dir.eq.6) then
         aveV = 60.0D0
        else if (axis_dir.eq.7) then
         aveV = 70.0D0
        else if (axis_dir.eq.8) then
         aveV = 80.0D0
        else if (axis_dir.eq.9) then
         aveV = 90.0D0
        else if (axis_dir.eq.10) then
         if (Bright_flag.eq.0) then
          aveV = 100.0D0
         else if (Bright_flag.eq.1) then
          aveV = 44.0D0
         else
          print *,"Bright_flag invalid"
          stop
         endif
        else
         print *,"axis_dir invalid for bubble formation"
         stop
        endif

        if (Bright_flag.eq.0) then
         density=0.9
         viscosity=10.0
         sigma=25.0
        else if (Bright_flag.eq.1) then
         density=0.732
         viscosity=0.0099
         sigma=22.5
        else
         print *,"Bright_flag invalid"
         stop
        endif

       else
        print *,"old flag invalid"
        stop
       endif

       Weber=(aveV**2)*radius*density/sigma
       Reynolds=density*radius*aveV/viscosity
       Froude=(aveV**2)/(radius*980.0)

       if (1.eq.0) then
        print *,"Weber,Reynolds,Froude ",Weber,Reynolds,Froude
        print *,"1/Weber,1/Reynolds,1/Froude ",one/Weber, &
          one/Reynolds,one/Froude
        stop
       endif
       phi=zero
       aveV=one
       radius=radblob
       if (abs(x).le.radius) then
        delta=xsten(1,1)-xsten(-1,1)
        if (delta.le.zero) then
         print *,"delta invalid bubble_formation_inflow_bc"
         stop
        endif
        phi=two*aveV*(1.0-(abs(x)/radius)**2-((delta/radius)**2)/four) 
       endif

      else
       print *,"need probtype.eq.25 in bubble_formation_inflow_bc"
       stop
      endif

      return
      end subroutine bubble_formation_inflow_bc

      subroutine vbc( velocity, t , yval,zval, error )

      use global_distance_module

      IMPLICIT NONE
! velocity is in m/s (or microns/microseconds), t is in microseconds
! zval in microns

      integer   error  
      
      real(amrex_real)    velocity, t, yval,zval
      
      integer  i,j,istar,jstar 
      real(amrex_real) zdiff,tdiff,tlocal
      
      real(amrex_real) NOD,NID,NPT,CHH,CHW,JLEN

      tlocal=t
      if (tlocal.lt.timehist_velbc(0)) then
       tlocal=timehist_velbc(0)
      endif

      error = 0
      call microfabparm(NOD,NID,NPT,CHH,CHW,JLEN)
!     print *,"in vbc yval,jlen-npt-193,tlocal,period,t(it) ",yval,
!    &  JLEN-NPT-193.0,tlocal,period_velbc,timehist(itime_velbc)
      if ((yval.le.half).or.(yval.ge.JLEN-NPT-193.0)) then
       velocity=zero
      else

      if (tlocal.gt.period_velbc) then
       velocity=zero
      else
       if (tlocal.ge.timehist_velbc(itime_velbc)) then
        velocity=zero
       else
        do j=0,itime_velbc-1
         if ((timehist_velbc(j).le.tlocal).and. &
             (timehist_velbc(j+1).ge.tlocal)) then
          jstar=j
         endif
        enddo
        if (AMREX_SPACEDIM.eq.2) then
         zval=half*zpos_velbc(ipos_velbc)
        endif
        if (zval.ge.zpos_velbc(ipos_velbc)) then
         velocity=zero
        else if (zval.le.zpos_velbc(1)) then
         velocity=zero
        else
         do i=1,ipos_velbc-1
          if ((zpos_velbc(i).le.zval).and. &
              (zpos_velbc(i+1).ge.zval)) then
           istar=i
          endif
         enddo
         tdiff=(timehist_velbc(jstar+1)-tlocal)/ &
           (timehist_velbc(jstar+1)-timehist_velbc(jstar))
         zdiff=(zpos_velbc(istar+1)-zval)/ &
           (zpos_velbc(istar+1)-zpos_velbc(istar))
         velocity=zdiff*tdiff*velbc_velbc(jstar,istar)+ &
           zdiff*(one-tdiff)*velbc_velbc(jstar+1,istar)+ &
           (one-zdiff)*tdiff*velbc_velbc(jstar,istar+1)+ &
           (one-zdiff)*(one-tdiff)*velbc_velbc(jstar+1,istar+1)
        endif
       endif
      endif
! tlocal in range
      endif
! yval in range          
      return
      end subroutine

      subroutine readpress( press_file, error )
      IMPLICIT NONE

!     * Read in the pressure boundary conditions provided by MicroFab for the 
!     * Okidata problems and store in the array pressbc. 

!     * This routine assumes that time is measured in seconds, the time values
!     * are equally spaced with spacing dt and that the pressure values are 
!     * known up to time t = 7.00D-05. For the MicroFab 
!     * test problems j= 0 ... = 70.  However,  to allow for future use of this
!     * type of BC, the size of the arrays time and pressbc may have to be 
!     * increased or passed as an argument to this subroutine.  Right now the
!     * number of data points in press_file is hard coded to be 71.

!     * Variables passed in ...

      character press_file*20

      integer   error  
      
! element 1-complicated 2-simple1 3-simple2


!     integer  1234567890, 1234567890, 1234567890, 1234567890, 1234567890

      integer  j         
      
!     * Open the data file containing the pressure BCs.  The name of the file
!     * is passed as an argument to the subroutine.

!      open(7, file=press_file, form="formatted", status="old", err=900)

      selectpress=1
      print *,"selectpress = ",selectpress

      do j = 0, 70
         
!        read(7,100) time_pressure_bcs(j), &
!         pressbc_pressure_bcs(j,1),  &
!         pressbc_pressure_bcs(j,2),  &
!         pressbc_pressure_bcs(j,3)
        
      end do
 
!      close(7)
 
      do j=0,70
!      print *,"time,pressure ", &
!        time_pressure_bcs(j), &
!        pressbc_pressure_bcs(j,selectpress)
      enddo

      dt_pressure_bcs = time_pressure_bcs(1) -  &
        time_pressure_bcs(0)
      error = 0
      
      return
      
! 100   format(' ', e10.6, e9.6, e9.6, e9.6 )
! 900   write(6,910) press_file
! 910   format(" Can't open file = ",a," stopping ...")

      error = -1
      
      return
      end subroutine

      subroutine readvel( vel_file, error )
      IMPLICIT NONE

      character vel_file*20

      integer   error  
! time is in microseconds, velocity in meter/s, zpos in microns 
    
      integer  i,j         

      period_velbc=28.0
      rigidwall_velbc=180.0
      itime_velbc=11
      ipos_velbc=11
      
!      open(7, file=vel_file, form="formatted", status="old", err=901)

!      read(7,101) timehist_velbc(0),zpos_velbc(1),zpos_velbc(2), &
!       zpos_velbc(3),zpos_velbc(4),zpos_velbc(5), &
!       zpos_velbc(6),zpos_velbc(7),zpos_velbc(8), &
!       zpos_velbc(9),zpos_velbc(10)
      zpos_velbc(11)=rigidwall_velbc

      do j = 0, itime_velbc
!       read(7,102) timehist_velbc(j),velbc_velbc(j,1),velbc_velbc(j,2), &
!        velbc_velbc(j,3),velbc_velbc(j,4), &
!        velbc_velbc(j,5),velbc_velbc(j,6),velbc_velbc(j,7), &
!        velbc_velbc(j,8),velbc_velbc(j,9),velbc_velbc(j,10)
       velbc_velbc(j,11)=zero
      end do
!  101 format(f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2)  
!  102 format(f7.2,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3,f7.3)  
      close(7)
 
      print *,"period,rigidwall ",period_velbc,rigidwall_velbc
      do j=1,ipos_velbc
       print *,"j,zpos ",j,zpos_velbc(j)
      enddo
      do j=0,itime_velbc
       print *,"j,timehist ",j,timehist_velbc(j)
      enddo
      do i=0,itime_velbc
      do j=1,ipos_velbc
       print *,"i,j,velbc ",i,j,velbc_velbc(i,j)
      enddo
      enddo
      
      error = 0
      
      return
      
! 901   write(6,911) vel_file
! 911   format(" Can't open file = ",a," stopping ...")

      error = -1
      
      return
      end subroutine


      subroutine localparam(Dbdry,grainbdry)
      IMPLICIT NONE
      real(amrex_real) Dbdry,grainbdry

      Dbdry=zero
      grainbdry=zero

      return
      end subroutine


      subroutine rtdist(x,y,dist)
      IMPLICIT NONE
      real(amrex_real) x,y,dist,rholevel

      dist=radblob*cos(xblob*Pi*x)-y
      rholevel=dist


      return
      end subroutine

        ! currently not used.  April 8, 2016
      subroutine get_ternary(im,im_opp,im_3,iten,iten_13, &
        iten_23,sin1,sin2,sin3,th1,th2,th3,tension)
      IMPLICIT NONE

      integer, INTENT(in) :: im,im_opp,im_3,iten,iten_13,iten_23
      real(amrex_real), INTENT(out) :: sin1,sin2,sin3
      real(amrex_real), INTENT(in) :: tension(num_interfaces)
      real(amrex_real) t1,t2,t3
      real(amrex_real) th1,th2,th3
      real(amrex_real) t1old,t2old,err
      integer iter,maxiter
      integer i12,i13,i23

      if ((num_materials.lt.3).or. &
          (im.le.0).or.(im.gt.num_materials).or. &
          (im_opp.le.0).or.(im_opp.gt.num_materials).or. &
          (im_3.le.0).or.(im_3.gt.num_materials).or. &
          (iten.le.0).or.(iten.gt.num_interfaces).or. & 
          (iten_13.le.0).or.(iten_13.gt.num_interfaces).or. & 
          (iten_23.le.0).or.(iten_23.gt.num_interfaces).or. &
          (tension(iten).lt.zero).or. &
          (tension(iten_13).lt.zero).or. &
          (tension(iten_23).lt.zero)) then
       print *,"parameter bust in get_ternary"
       stop
      endif

      maxiter=200
      iter=1

      if ((tension(iten).ge.tension(iten_13)).and. &
          (tension(iten).ge.tension(iten_23))) then
       i12=iten
       i13=iten_13
       i23=iten_23
      else if (tension(iten_13).ge.tension(iten_23)) then
       i12=iten_13
       i13=iten
       i23=iten_23
      else
       i12=iten_23
       i13=iten
       i23=iten_13
      endif

      t1=two*Pi/three 
      t2=two*Pi/three
      t3=two*Pi-t1-t2
      sin1=sin(t1)
      sin2=sin(t2)
      sin3=sin(t3)

      if (tension(i12).gt.zero) then

       err=1.0D+20
       do while ((err.gt.VOFTOL).and.(iter.lt.maxiter))
        t3=two*Pi-t1-t2
        t1old=t1
        t2old=t2
        t1=asin(tension(i23)*sin(t3)/tension(i12))
        t2=asin(tension(i13)*sin(t3)/tension(i12))
        if (t1.lt.zero) then
         t1=t1+two*Pi
        endif
        if (t1.gt.two*Pi) then
         t1=t1-two*Pi
        endif
        if (t2.lt.zero) then
         t2=t2+two*Pi
        endif
        if (t2.gt.two*Pi) then
         t2=t2-two*Pi
        endif
        err=abs(t1-t1old)+abs(t2-t2old)
        iter=iter+1
       enddo
       if (iter.ge.maxiter) then
        print *,"iter exceeds maxiter in get_ternary"
        stop
       endif

       t3=two*Pi-t1-t2

       if (i12.eq.iten) then
        th3=t3  ! t3 opposite i12
        th2=t2  ! t2 opposite i13
        th1=t1  ! t1 opposite i23
       else if (i12.eq.iten_13) then
        th2=t3  ! th2 opposite iten_13, t3 opposite i12
        th3=t2  ! th3 opposite iten, t2 opposite i13
        th1=t1  ! th1 opposite iten_23, t1 opposite i23
       else if (i12.eq.iten_23) then
        th1=t3  ! th1 opposite iten_23, t3 opposite i12
        th3=t2  ! th3 opposite iten, t2 opposite i13
        th2=t1  ! th2 opposite iten_13, t1 opposite i23
       else
        print *,"i12 invalid" 
        stop
       endif
       sin1=sin(th1)   
       sin2=sin(th2)   
       sin3=sin(th3)
 
      else if (tension(i12).eq.zero) then
       th1=t1
       th2=t2
       th3=t3
       sin1=sin(th1)   
       sin2=sin(th2)   
       sin3=sin(th3)
      else
       print *,"tension invalid"
       stop
      endif

      if (1.eq.1) then
       print *,"iter,im,im_opp,im_3,th1,th2,th3 ", &
        iter,im,im_opp,im_3,th1,th2,th3
      endif

      return
      end subroutine get_ternary



        ! input: LSleft,LSright,num_materials
        ! output: gradh,im,im_opp
        ! gradh=H(LSright)-H(LSleft) 
        !
        ! contact line treatment is done if "im" and "im_opp"
        ! are majority fluids in adjoining cells.  Then "im_3"
        ! is a majority fluid in a 3rd cell of the 3x3x3 stencil.
        ! If there are more than one "im_3", pick the one with the
        ! largest value.
        !
      subroutine fluid_interface( &
        LSleft,LSright, &
        gradh, & !intent(out)
        im_opp,im, & !intent(out)
        imL,imR) !intent(out)

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(out) :: im_opp,im
      real(amrex_real), INTENT(in) :: LSleft(num_materials)
      real(amrex_real), INTENT(in) :: LSright(num_materials)
      real(amrex_real), INTENT(out) :: gradh
      integer, INTENT(out) :: imL,imR

      im=0
      im_opp=0
      gradh=zero

      call get_primary_material(LSleft,imL)
      call get_primary_material(LSright,imR)

      if ((imL.lt.1).or.(imL.gt.num_materials).or. &
          (imR.lt.1).or.(imR.gt.num_materials)) then
       print *,"imL or imR invalid"
       stop
      endif

      if (is_rigid(imL).eq.1) then
       ! do nothing
      else if (is_rigid(imR).eq.1) then
       ! do nothing
      else if (imL.eq.imR) then
       ! do nothing
      else if ((is_rigid(imL).eq.0).and. &
               (is_rigid(imR).eq.0).and. &
               (imL.ne.imR)) then

       if (imL.lt.imR) then 
        gradh=-one
        im=imL
        im_opp=imR
       else if (imL.gt.imR) then 
        gradh=one
        im=imR
        im_opp=imL
       else
        print *,"imL or imR bust"
        stop
       endif

      else
       print *,"is_rigid, imL, or imR invalid PROB.F90"
       stop
      endif

      return
      end subroutine fluid_interface

       ! partid=0..nparts-1
       ! im_solid=0..num_materials
       ! is_solid_face==1 if:
       !   0.0<=facecut_solid<=VOFTOL_AREAFRAC  or
       !   max(LSleft(im_solid),LSright(im_solid))>=0.0
      subroutine fixed_face( &
       facecut_solid, &   !intent(in)
       LSleft,LSright, &  !intent(in)
       is_solid_face, &
       im_solid, &
       im_solid_valid, &
       partid_solid)
      use global_utility_module
      IMPLICIT NONE

        !surface tension coefficient is zero
      real(amrex_real), INTENT(in) :: facecut_solid 
      real(amrex_real), INTENT(in) :: LSleft(num_materials)
      real(amrex_real), INTENT(in) :: LSright(num_materials)
      real(amrex_real) LScrit_solid,LStest
      integer, INTENT(out) :: is_solid_face
      integer im
      integer, INTENT(out) :: im_solid
      integer, INTENT(out) :: im_solid_valid
      integer, INTENT(out) :: partid_solid
      integer nparts

      is_solid_face=0
      im_solid=0
      im_solid_valid=0

      nparts=0
      partid_solid=-1

      do im=1,num_materials

        ! FSI_PRESCRIBED_PROBF90
        ! FSI_PRESCRIBED_NODES
        ! FSI_SHOELE_CTML
        ! FSI_ICE_NODES_INIT
        ! FSI_FLUID_NODES_INIT
       if (is_lag_part(im).eq.1) then

        ! if is_rigid==1:
        ! FSI_PRESCRIBED_PROBF90
        ! FSI_PRESCRIBED_NODES
        ! FSI_SHOELE_CTML
        if (is_rigid(im).eq.0) then

         ! FSI_PRESCRIBED_PROBF90
         ! FSI_PRESCRIBED_NODES
         ! FSI_SHOELE_CTML
        else if (is_rigid(im).eq.1) then

         if (im_solid.eq.0) then
          im_solid=im
          partid_solid=nparts
          LScrit_solid=max(LSleft(im),LSright(im))
         else if ((im_solid.ge.1).and. &
                  (im_solid.le.num_materials)) then
          LStest=max(LSleft(im),LSright(im))
          if (LStest.gt.LScrit_solid) then
           LScrit_solid=LStest
           im_solid=im
           partid_solid=nparts
          else if (LStest.le.LScrit_solid) then
           ! do nothing
          else
           print *,"LStest or LScrit_solid bust"
           stop
          endif
         else
          print *,"im_solid invalid 7"
          stop
         endif

        else
         print *,"is_rigid invalid PROB.F90"
         stop
        endif

        nparts=nparts+1

       else if (is_lag_part(im).eq.0) then

        if (is_rigid(im).eq.0) then
         ! do nothing
        else 
         print *,"is_rigid(im) invalid"
         stop
        endif
       
       else
        print *,"is_lag_part invalid"
        stop
       endif

      enddo ! im=1..num_materials

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fixed_face"
       stop
      endif

      if ((facecut_solid.le.VOFTOL_AREAFRAC).and. &
          (facecut_solid.ge.zero)) then
       is_solid_face=1
       if (im_solid.eq.0) then
        ! do nothing
       else if ((im_solid.ge.1).and.(im_solid.le.num_materials)) then
        im_solid_valid=1 
        if ((partid_solid.lt.0).or. &
            (partid_solid.ge.nparts)) then
         print *,"partid_solid invalid"
         stop
        endif 
       else
        print *,"im_solid invalid 20"
        stop
       endif
      else if ((facecut_solid.ge.VOFTOL_AREAFRAC).and. &
               (facecut_solid.le.one)) then
       if (im_solid.eq.0) then
        ! do nothing
       else if ((im_solid.ge.1).and. &
                (im_solid.le.num_materials)) then
        im_solid_valid=1
        if (LScrit_solid.ge.zero) then
         is_solid_face=1
        else if (LScrit_solid.lt.zero) then
         ! do nothing
        else
         print *,"LScrit_solid invalid: ",LScrit_solid
         stop
        endif
        if ((partid_solid.lt.0).or. &
            (partid_solid.ge.nparts)) then
         print *,"partid_solid invalid: ",partid_solid
         stop
        endif 
       else
        print *,"im_solid invalid 9"
        stop
       endif
      else
       print *,"facecut_solid invalid: ",facecut_solid
       stop
      endif
  
      return
      end subroutine fixed_face

       !sin(theta1)/sigma23 = sin(theta2)/sigma13 = sin(theta3)/sigma12
       !if theta_air=Pi => sigma_ice_melt=0.0
       !In general: if sigma_{ij}=0.0, then merge materials i and j.
       !For contact line dynamics:
       !sigma_{jk}-sigma_{ik}=sigma_{ij}cos(theta_{i})
       !if sigma_{jk}=0 and sigma_{ik}=sigma_{ij} => theta_{i}=180 deg.
       ! => merge materials j and k.
       !if sigma_{ik}=0 and sigma_{jk}=sigma_{ij} => theta_{i}=0 deg.
       ! => merge materials i and k.
      subroutine merge_levelset(xpos,time,LS,LS_merge)
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: xpos(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: LS(num_materials)
      real(amrex_real), INTENT(out) :: LS_merge(num_materials)
      integer im,im_opp
      integer iten
      integer im_primary
      integer im_secondary
      integer im_tertiary
      real(amrex_real) :: user_tension(num_interfaces)
      real(amrex_real) :: def_thermal(num_materials)

      do im=1,num_materials
       def_thermal(im)=room_temperature
       LS_merge(im)=LS(im)
      enddo

      do im=1,num_materials
       if ((is_ice(im).eq.1).or. &
           (is_rigid(im).eq.1).or. &
           (is_rigid_CL(im).eq.1)) then
        do im_opp=1,num_materials
         if (im_opp.ne.im) then
          if ((is_rigid(im_opp).eq.0).and. &
              (is_ice(im_opp).eq.0).and. &
              (is_rigid_CL(im_opp).eq.0)) then
           call get_iten(im,im_opp,iten)
           if ((iten.ge.1).and.(iten.le.num_interfaces)) then
            ! do nothing
           else
            print *,"iten invalid: ",iten
            stop
           endif

           call get_user_tension( &
            xpos,time,fort_tension,user_tension,def_thermal)

              ! "merge_levelset" is for the algorithm described by:
              ! Lyu, Wang, Zhang, Pedrono, Sun, Legendre JCP 2021
              ! sigma_ice_melt=0 => theta_ambient=0 (=>growth_angle=0)
           if (user_tension(iten).eq.zero) then

             !get_primary_material is declared in GLOBALUTIL.F90
            call get_primary_material(LS,im_primary)
             !get_secondary_material is declared in MOF.F90
            call get_secondary_material(LS,im_primary,im_secondary)
             !get_tertiary_material is declared in MOF.F90
             !is_rigid(im_tertiary)=0
            call get_tertiary_material(LS, &
              im_primary,im_secondary,im_tertiary)
            if (im_tertiary.eq.0) then
             ! do nothing (there is no surface tension in this case)
            else if ((im_tertiary.ge.1).and. &
                     (im_tertiary.le.num_materials)) then

             LS_merge(im)=-99999.0d0 ! delete the ice/rigid/rigid_CL material.

             if (im_primary.eq.im) then ! ice/rigid/rigid_CL was primary

              if (im_secondary.eq.im_opp) then !"water" was secondary
               LS_merge(im_opp)=-LS(im_tertiary) !im_tertiary="air"
              else if (im_secondary.ne.im_opp) then
               LS_merge(im_opp)=-LS(im_secondary) !im_secondary="air"
              else
               print *,"im_secondary bust: ",im_secondary
               stop
              endif

             else if (im_primary.eq.im_opp) then !"water" is primary

              if (im_secondary.eq.im) then ! ice/rigid/rigid_CL is secondary
               LS_merge(im_opp)=-LS(im_tertiary) !im_tertiary="air"
              else
               ! do nothing (im_secondary is "air")
              endif

             else if ((im_primary.ne.im).and. &
                      (im_primary.ne.im_opp)) then

              ! "water" replaces ice/rigid/rigid_CL
              if (im_secondary.eq.im) then
               LS_merge(im_opp)=LS(im) 
              else if (im_secondary.eq.im_opp) then
               ! do nothing
              else if ((im_secondary.ne.im).and. &
                       (im_secondary.ne.im_opp)) then
               if (LS(im).gt.LS(im_opp)) then
                LS_merge(im_opp)=LS(im) 
               else if (LS(im).le.LS(im_opp)) then
                ! do nothing
               else
                print *,"LS(im) invalid: ",im,LS(im)
                stop 
               endif
              else
               print *,"im_secondary invalid: ",im_secondary
               stop
              endif

             else
              print *,"im_primary invalid: ",im_primary
              stop
             endif

            else
             print *,"im_tertiary invalid: ",im_tertiary
             stop
            endif

           else if (user_tension(iten).gt.zero) then
            !do nothing
           else
            print *,"user_tension invalid: ",user_tension(iten)
            stop
           endif
          else if ((is_rigid(im_opp).eq.1).or. &
                   (is_ice(im_opp).eq.1).or. &
                   (is_rigid_CL(im_opp).eq.1)) then
           ! do nothing
          else
           print *,"im_opp inconsistency: ",im_opp
           print *,"is_rigid(im_opp): ",is_rigid(im_opp)
           print *,"is_ice(im_opp): ",is_ice(im_opp)
           print *,"is_rigid_CL(im_opp): ",is_rigid_CL(im_opp)
           stop
          endif
         else if (im_opp.eq.im) then
          ! do nothing
         else
          print *,"im_opp or im bust: ",im,im_opp
          stop
         endif
        enddo !im_opp=1,num_materials
       else if ((is_ice(im).eq.0).and. &
                (is_rigid(im).eq.0).and. &
                (is_rigid_CL(im).eq.0)) then
        ! do nothing
       else
        print *,"im inconsistency: ",im
        print *,"is_rigid(im): ",is_rigid(im)
        print *,"is_ice(im): ",is_ice(im)
        print *,"is_rigid_CL(im): ",is_rigid_CL(im)
        stop
       endif
      enddo !im=1,num_materials

      end subroutine merge_levelset


       !sin(theta1)/sigma23 = sin(theta2)/sigma13 = sin(theta3)/sigma12
       !if theta_air=Pi => sigma_ice_melt=0
       !In general: if sigma_{ij}=0.0, then merge materials i and j.
       !For contact line dynamics:
       !sigma_{jk}-sigma_{ik}=sigma_{ij}cos(theta_{i})
       !if sigma_{jk}=0 and sigma_{ik}=sigma_{ij} => theta_{i}=180 deg.
       ! => merge materials j and k.
       !if sigma_{ik}=0 and sigma_{jk}=sigma_{ij} => theta_{i}=0 deg.
       ! => merge materials i and k.
      subroutine merge_normal(xpos,time,LS,nrm,nrm_merge)
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: xpos(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: LS(num_materials)
      real(amrex_real), INTENT(in) :: nrm(num_materials*SDIM)
      real(amrex_real), INTENT(out) :: nrm_merge(num_materials*SDIM)
      integer im,im_opp
      integer iten
      integer dir
      integer im_primary
      integer im_secondary
      integer im_tertiary
      real(amrex_real) :: user_tension(num_interfaces)
      real(amrex_real) :: def_thermal(num_materials)

      do im=1,num_materials
       def_thermal(im)=room_temperature
      enddo
      do im=1,num_materials*SDIM
       nrm_merge(im)=nrm(im)
      enddo
      do im=1,num_materials
       if ((is_ice(im).eq.1).or. &
           (is_rigid(im).eq.1).or. &
           (is_rigid_CL(im).eq.1)) then
        do im_opp=1,num_materials
         if (im_opp.ne.im) then
          if ((is_rigid(im_opp).eq.0).and. &
              (is_ice(im_opp).eq.0).and. &
              (is_rigid_CL(im_opp).eq.0)) then
           call get_iten(im,im_opp,iten)
           if ((iten.ge.1).and.(iten.le.num_interfaces)) then
            ! do nothing
           else
            print *,"iten invalid: ",iten
            stop
           endif

           call get_user_tension( &
             xpos,time,fort_tension,user_tension,def_thermal)

              ! "merge_normal" is for the algorithm described by:
              ! Lyu, Wang, Zhang, Pedrono, Sun, Legendre JCP 2021
              ! sigma_ice_melt=0 => theta_ambient=0 (=>growth_angle=0)
           if (user_tension(iten).eq.zero) then

             !get_primary_material is declared in GLOBALUTIL.F90
            call get_primary_material(LS,im_primary)
             !get_secondary_material is declared in MOF.F90
            call get_secondary_material(LS,im_primary,im_secondary)
             !get_tertiary_material is declared in MOF.F90
             !is_rigid(im_tertiary)=0
            call get_tertiary_material(LS, &
              im_primary,im_secondary,im_tertiary)
            if (im_tertiary.eq.0) then
             ! do nothing (there is no surface tension in this case)
            else if ((im_tertiary.ge.1).and. &
                     (im_tertiary.le.num_materials)) then

             if (im_primary.eq.im) then ! ice/rigid/rigid_CL was primary

              if (im_secondary.eq.im_opp) then !"water" was secondary
               do dir=1,SDIM
                nrm_merge((im_opp-1)*SDIM+dir)= &
                      -nrm((im_tertiary-1)*SDIM+dir)
               enddo
              else if (im_secondary.ne.im_opp) then
               do dir=1,SDIM
                nrm_merge((im_opp-1)*SDIM+dir)= &
                      -nrm((im_secondary-1)*SDIM+dir)
               enddo
              else
               print *,"im_secondary bust: ",im_secondary
               stop
              endif

             else if (im_primary.eq.im_opp) then !"water" is primary

              if (im_secondary.eq.im) then ! ice/rigid/rigid_CL is secondary
               do dir=1,SDIM
                nrm_merge((im_opp-1)*SDIM+dir)= &
                     -nrm((im_tertiary-1)*SDIM+dir)
               enddo
              else
               ! do nothing (im_secondary is "air")
              endif

             else if ((im_primary.ne.im).and. &
                      (im_primary.ne.im_opp)) then

              ! "water" replaces ice/rigid/rigid_CL
              if (im_secondary.eq.im) then
               do dir=1,SDIM
                nrm_merge((im_opp-1)*SDIM+dir)= &
                      nrm((im-1)*SDIM+dir)
               enddo
              else if (im_secondary.eq.im_opp) then
               ! do nothing
              else if ((im_secondary.ne.im).and. &
                       (im_secondary.ne.im_opp)) then
               if (LS(im).gt.LS(im_opp)) then
                do dir=1,SDIM
                 nrm_merge((im_opp-1)*SDIM+dir)= &
                      nrm((im-1)*SDIM+dir)
                enddo
               else if (LS(im).le.LS(im_opp)) then
                ! do nothing
               else
                print *,"LS(im) invalid: ",im,LS(im)
                stop 
               endif
              else
               print *,"im_secondary invalid: ",im_secondary
               stop
              endif

             else
              print *,"im_primary invalid: ",im_primary
              stop
             endif

            else
             print *,"im_tertiary invalid: ",im_tertiary
             stop
            endif

           else if (user_tension(iten).gt.zero) then
            ! do nothing
           else
            print *,"user_tension invalid: ",user_tension(iten)
            stop
           endif
          else if ((is_rigid(im_opp).eq.1).or. &
                   (is_ice(im_opp).eq.1).or. &
                   (is_rigid_CL(im_opp).eq.1)) then
           ! do nothing
          else
           print *,"im_opp inconsistency: ",im_opp
           print *,"is_rigid(im_opp): ",is_rigid(im_opp)
           print *,"is_ice(im_opp): ",is_ice(im_opp)
           print *,"is_rigid_CL(im_opp): ",is_rigid_CL(im_opp)
           stop
          endif
         else if (im_opp.eq.im) then
          ! do nothing
         else
          print *,"im_opp or im bust: ",im,im_opp
          stop
         endif
        enddo !im_opp=1,num_materials
       else if ((is_ice(im).eq.0).and. &
                (is_rigid(im).eq.0).and. &
                (is_rigid_CL(im).eq.0)) then
        ! do nothing
       else
        print *,"im inconsistency: ",im
        print *,"is_rigid(im): ",is_rigid(im)
        print *,"is_ice(im): ",is_ice(im)
        print *,"is_rigid_CL(im): ",is_rigid_CL(im)
        stop
       endif
      enddo !im=1,num_materials

      end subroutine merge_normal


       !sin(theta1)/sigma23 = sin(theta2)/sigma13 = sin(theta3)/sigma12
       !if theta_air=Pi => sigma_ice_melt=0
       !In general: if sigma_{ij}=0.0, then merge materials i and j.
       !For contact line dynamics:
       !sigma_{jk}-sigma_{ik}=sigma_{ij}cos(theta_{i})
       !if sigma_{jk}=0 and sigma_{ik}=sigma_{ij} => theta_{i}=180 deg.
       ! => merge materials j and k.
       !if sigma_{ik}=0 and sigma_{jk}=sigma_{ij} => theta_{i}=0 deg.
       ! => merge materials i and k.
      subroutine merge_vof(xpos,time,vof,vof_merge)
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE
      real(amrex_real), INTENT(in) :: xpos(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: vof(num_materials)
      real(amrex_real), INTENT(out) :: vof_merge(num_materials)
      integer im,im_opp
      integer iten
      real(amrex_real) :: user_tension(num_interfaces)
      real(amrex_real) :: def_thermal(num_materials)

      do im=1,num_materials
       def_thermal(im)=room_temperature
       vof_merge(im)=vof(im)
      enddo
      do im=1,num_materials
       if ((is_ice(im).eq.1).or. &
           (is_rigid(im).eq.1).or. &
           (is_rigid_CL(im).eq.1)) then
        do im_opp=1,num_materials
         if (im_opp.ne.im) then
          if ((is_rigid(im_opp).eq.0).and. &
              (is_ice(im_opp).eq.0).and. &
              (is_rigid_CL(im_opp).eq.0)) then
           call get_iten(im,im_opp,iten)
           if ((iten.ge.1).and.(iten.le.num_interfaces)) then
            ! do nothing
           else
            print *,"iten invalid: ",iten
            stop
           endif

           call get_user_tension( &
             xpos,time,fort_tension,user_tension,def_thermal)

             ! "merge_vof" code is for the algorithm described by:
             ! Lyu, Wang, Zhang, Pedrono, Sun, Legendre JCP 2021
             ! sigma_ice_melt=0 => theta_ambient=0 (=>growth_angle=0)
           if (user_tension(iten).eq.zero) then

            vof_merge(im)=zero
            vof_merge(im_opp)=vof(im_opp)+vof(im)

           else if (user_tension(iten).gt.zero) then
            ! do nothing
           else
            print *,"user_tension invalid: ",user_tension(iten)
            stop
           endif
          else if ((is_rigid(im_opp).eq.1).or. &
                   (is_ice(im_opp).eq.1).or. &
                   (is_rigid_CL(im_opp).eq.1)) then
           ! do nothing
          else
           print *,"im_opp inconsistency: ",im_opp
           print *,"is_rigid(im_opp): ",is_rigid(im_opp)
           print *,"is_ice(im_opp): ",is_ice(im_opp)
           print *,"is_rigid_CL(im_opp): ",is_rigid_CL(im_opp)
           stop
          endif

         else if (im_opp.eq.im) then
          ! do nothing
         else
          print *,"im_opp or im bust: ",im,im_opp
          stop
         endif
        enddo !im_opp=1,num_materials
       else if ((is_ice(im).eq.0).and. &
                (is_rigid(im).eq.0).and. &
                (is_rigid_CL(im).eq.0)) then
        ! do nothing
       else
        print *,"im inconsistency: ",im
        print *,"is_rigid(im): ",is_rigid(im)
        print *,"is_ice(im): ",is_ice(im)
        print *,"is_rigid_CL(im): ",is_rigid_CL(im)
        stop
       endif
      enddo !im=1,num_materials

      end subroutine merge_vof

       ! called from:
       !   fort_init_physics_vars
       !   fort_cell_to_mac (surface tension force on MAC grid)
      subroutine fluid_interface_tension( &
         xpos,time, &
         LSleft,LSright, &
         gradh, & !INTENT(out)
         im_opp,im, & !INTENT(out)
         imL,imR)     !INTENT(out)

      use global_utility_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: xpos(SDIM)
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(out) :: im_opp,im
      real(amrex_real), INTENT(in) :: LSleft(num_materials)
      real(amrex_real), INTENT(in) :: LSright(num_materials)
      real(amrex_real), INTENT(out) :: gradh
      real(amrex_real) :: LSleft_merge(num_materials)
      real(amrex_real) :: LSright_merge(num_materials)
      integer, INTENT(out) :: imL,imR

      im=0
      im_opp=0
      gradh=zero

      call merge_levelset(xpos,time,LSleft,LSleft_merge)
      call merge_levelset(xpos,time,LSright,LSright_merge)

      call get_primary_material(LSleft_merge,imL)
      call get_primary_material(LSright_merge,imR)

      if ((imL.lt.1).or.(imL.gt.num_materials).or. &
          (imR.lt.1).or.(imR.gt.num_materials)) then
       print *,"imL or imR invalid (fluid_interface_tension): ",imL,imR
       stop
      endif

      if (is_rigid_CL(imL).eq.1) then
       ! do nothing
      else if (is_rigid_CL(imR).eq.1) then
       ! do nothing
      else if (imL.eq.imR) then
       ! do nothing
      else if ((is_rigid_CL(imL).eq.0).and. &
               (is_rigid_CL(imR).eq.0).and. &
               (imL.ne.imR)) then

       if (imL.lt.imR) then
        gradh=-one
        im=imL
        im_opp=imR
       else if (imL.gt.imR) then
        gradh=one
        im=imR
        im_opp=imL
       else
        print *,"imL or imR bust"
        stop
       endif

      else
       print *,"is_rigid_CL, imL, or imR invalid PROB.F90"
       stop
      endif 

      return
      end subroutine fluid_interface_tension


      subroutine adapt_at_nozzle(adapt_nozzle_flag)
      IMPLICIT NONE

      integer adapt_nozzle_flag


      adapt_nozzle_flag=0


      if ((probtype.eq.532).or. &
          (probtype.eq.538).or. &  ! inputs.injA
          (probtype.eq.537).or. &
          (probtype.eq.541)) then
       adapt_nozzle_flag=1

      else if (SDIM.eq.2) then

       ! do nothing

      else if (SDIM.eq.3) then

        ! for impinging jet with like materials, or diesel jet,
        ! always adapt at the nozzle. (unless otherwise prescribed)
       if ((probtype.eq.53).or.(probtype.eq.536).or. &
           (probtype.eq.530)) then
        adapt_nozzle_flag=1
       endif

      else
       print *,"dimension bust"
       stop
      endif
       
      return
      end subroutine adapt_at_nozzle


       ! called from fort_sloperecon and fort_initdata
       ! (note see NavierStokes::init_pressure_error_indicator() too,
       !  which calls fort_pressure_indicator, which calls
       !  EOS_error_ind; init_pressure_error_indicator() is called from 
       !  multiphase_project with project_option=SOLVETYPE_PRES)
      subroutine calc_error_indicator( &
        level, &
        max_level, &
        xsten, &
        nhalf, &
        dx, &
        bfact, &
        voflist, &
        LS_stencil, &
        err, &
        time)

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: max_level
      integer, INTENT(in) :: nhalf
      integer, INTENT(in) :: bfact

      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: &
            LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)
      real(amrex_real), INTENT(in) :: voflist(num_materials)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(out) :: err

      integer inear
      integer im,im_primary

      real(amrex_real) dxmin
      real(amrex_real) dist,dist3
      integer i1,j1,k1
      integer material_count
      integer adapt_nozzle_flag
      integer material_present_flag(num_materials)
      integer k1lo,k1hi
      real(amrex_real) vfrac_rigid_sum
      real(amrex_real) LS_temp(num_materials)

      if ((level.lt.0).or.(level.gt.max_level)) then
       print *,"level invalid calc_error_indicator"
       stop
      endif

      if (nhalf.lt.2) then
       print *,"nhalf invalid in calc_error_indicator"
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
 
      call get_dxmin(dx,bfact,dxmin)
      if (dxmin.gt.zero) then
       ! do nothing
      else
       print *,"dxmin invalid"
       stop
      endif

      inear=0 

      call adapt_at_nozzle(adapt_nozzle_flag)

      if (adapt_nozzle_flag.eq.1) then

       im=1
       call materialdist(xsten,nhalf,dx,bfact,dist,im,time)
       if ((FSI_flag(im).eq.FSI_FLUID).or. & 
           (FSI_flag(im).eq.FSI_FLUID_NODES_INIT)) then 
        ! do nothing
       else
        print *,"FSI_flag(im) invalid"
        print *,"im=",im
        print *,"FSI_flag(im)=",FSI_flag(im)
        stop
       endif
       if (abs(dist).le.two*dxmin) then
        inear=2
       else if (abs(dist).ge.two*dxmin) then
        ! do nothing
       else
        print *,"dist is NaN"
        stop
       endif

       if (num_materials.eq.2) then
        ! do nothing
       else if (num_materials.gt.2) then
        im=3
        if ((FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90).or. &
            (FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. &
            (FSI_flag(im).eq.FSI_ICE_PROBF90).or. &
            (FSI_flag(im).eq.FSI_ICE_STATIC).or. &
            (FSI_flag(im).eq.FSI_ICE_NODES_INIT).or. &
            (FSI_flag(im).eq.FSI_SHOELE_CTML).or. &
            (FSI_flag(im).eq.FSI_RIGID_NOTPRESCRIBED)) then
         ! do nothing
        else if ((FSI_flag(im).eq.FSI_FLUID).or. &
                 (FSI_flag(im).eq.FSI_FLUID_NODES_INIT)) then  
         call materialdist(xsten,nhalf,dx,bfact,dist3,im,time)
         if (abs(dist3).le.two*dxmin) then
          inear=2
         else if (abs(dist3).ge.two*dxmin) then
          ! do nothing
         else
          print *,"dist3 is NaN"
          stop
         endif
        else
         print *,"FSI_flag(im) invalid in calc_error_indicator"
         print *,"im=",im
         print *,"FSI_flag(im)=",FSI_flag(im)
         stop
        endif
       else
        print *,"num_materials invalid"
        stop
       endif

      else if (adapt_nozzle_flag.eq.0) then
       ! do nothing
      else
       print *,"adapt_nozzle_flag invalid"
       stop
      endif  

      do im=1,num_materials 
       material_present_flag(im)=0
      enddo

      vfrac_rigid_sum=zero
      do im=1,num_materials 
       LS_temp(im)=LS_stencil(D_DECL(0,0,0),im)
       if (is_rigid(im).eq.0) then
        ! do nothing
       else if (is_rigid(im).eq.1) then
        vfrac_rigid_sum=vfrac_rigid_sum+voflist(im)
       else
        print *,"is_rigid(im) invalid"
        stop
       endif
      enddo ! im=1..num_materials

      if ((vfrac_rigid_sum.ge.one-VOFTOL).and. &
          (vfrac_rigid_sum.le.one+EPS1)) then
       vfrac_rigid_sum=one
      else if ((vfrac_rigid_sum.ge.-EPS1).and. &
               (vfrac_rigid_sum.le.VOFTOL)) then
       vfrac_rigid_sum=zero
      else if ((vfrac_rigid_sum.gt.zero).and. &
               (vfrac_rigid_sum.lt.one)) then
       ! do nothing
      else
       print *,"vfrac_rigid_sum invalid: ",vfrac_rigid_sum
       stop
      endif

      call get_primary_material(LS_temp,im_primary)

      material_present_flag(im_primary)=1

      if ((is_rigid(im_primary).eq.0).or. &   ! fluid primary
          ((is_rigid(im_primary).eq.1).and. & ! solid primary, but fluids
           (vfrac_rigid_sum.le.one-VOFTOL))) then  ! in the cell.

       do im=1,num_materials 
        if ((voflist(im).ge.VOFTOL).and. &
            (voflist(im).le.one+EPS1)) then
         material_present_flag(im)=1
        else if ((voflist(im).ge.-EPS1).and. &
                 (voflist(im).le.VOFTOL)) then
         ! do nothing
        else
         print *,"voflist invalid: ",voflist(im)
         stop
        endif
       enddo ! im=1..num_materials

       do k1=k1lo,k1hi
       do j1=-1,1
       do i1=-1,1
        do im=1,num_materials 
         LS_temp(im)=LS_stencil(D_DECL(i1,j1,k1),im)
        enddo
        call get_primary_material(LS_temp,im_primary)
        material_present_flag(im_primary)=1
        if (is_rigid(im_primary).eq.0) then
         if (probtype.eq.46) then ! cavitation
          if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
           im=2 ! jwl
           if (LS_temp(im).gt.-dxmin) then
            material_present_flag(im)=1
           endif
          else if (axis_dir.eq.10) then
           ! do nothing (sphere impact)
          else if (axis_dir.eq.20) then
           ! do nothing (CODY ESTEBE created test problem)
          else
           print *,"axis_dir invalid"
           stop
          endif
         endif 
        else if (is_rigid(im_primary).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid PROB.F90"
         stop
        endif
       enddo
       enddo
       enddo ! i1,j1,k1

      else if ((is_rigid(im_primary).eq.1).and. &
               (vfrac_rigid_sum.gt.one-VOFTOL)) then
       ! do nothing
      else
       print *,"is_rigid or vfrac_rigid_sum invalid"
       stop
      endif

      material_count=0

      do im=1,num_materials 

       if (material_present_flag(im).eq.1) then

        material_count=material_count+1

       else if (material_present_flag(im).eq.0) then
        ! do nothing
       else 
        print *,"material_present_flag invalid"
        stop
       endif
        
      enddo ! im=1..num_materials

      if (material_count.gt.num_materials) then
       print *,"material_count is corrupt"
       stop
      endif

      if (material_count.ge.2) then
       inear=max(inear,1)
      endif
      if (material_count.ge.3) then
       inear=max(inear,2)
      endif

       ! if err>0.0 => adapt
      if (inear.eq.1) then
       err=one
      else if (inear.eq.2) then
       err=one
      else if (inear.eq.0) then
       err=zero
      else
       print *,"inear invalid"
       stop
      endif 
 
      return
      end subroutine calc_error_indicator


       ! called from get_symmetric_error
      subroutine exactdist( &
        xsten, &
        nhalf, &
        bfact, &
        dx, &
        dist, &
        imaterial, &
        time)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(out) :: dist
      integer, INTENT(in) :: imaterial
      real(amrex_real), INTENT(in) :: time
      integer dir
      real(amrex_real) x,y,z
      real(amrex_real) xstar,ystar,zstar
      real(amrex_real) xvec(SDIM)
      real(amrex_real) distsolid,distgas,dist_liquid,dist_ice
      integer im_solid_exactdist
      real(amrex_real) LS(num_materials)
      real(amrex_real) diamblob

      im_solid_exactdist=im_solid_primary()

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid exact dist"
       stop
      endif
 
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      do dir=1,SDIM
       xvec(dir)=xsten(0,dir)
      enddo

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        !do nothing
       else
        print *,"expecting z=y 5990"
        stop
       endif
      endif

      xstar=x
      ystar=y
      zstar=z

      if (fort_is_passive_advect_test().eq.1) then

       call SUB_LS(xvec,time,LS,num_materials)
       dist=LS(imaterial)

       ! drop on slope (exactdist)
      else if (probtype.eq.55) then
     
       if (SDIM.eq.2) then

        if ((num_materials.eq.3).and. &
            (im_solid_exactdist.eq.3).and. &
            (axis_dir.eq.0).and. &
            (radblob3.eq.zero).and. &
            (radblob5.eq.zero).and. &
            (radblob6.eq.zero).and. &
            (radblob7.eq.zero).and. &
            (abs(xblob-xblob2).lt.1.0E-6).and. &
            (abs(yblob-yblob2).lt.1.0E-6)) then

         if ((imaterial.eq.1).or.(imaterial.eq.2)) then
          ! distsolid>0 in solid
          call materialdist(xsten,nhalf,dx,bfact,distsolid,3,time)
          ! in: exactdist (maxtall=2 * radblob => dist_ice=dist_liquid)
          diamblob=two*radblob
          call drop_slope_dist(xstar,ystar,zstar, &
           time,diamblob,dist_ice,dist_liquid)
          distgas=-dist_liquid

          if (imaterial.eq.1) then
           dist=dist_liquid
          else if (imaterial.eq.2) then
           dist=distgas
          else
           print *,"imaterial invalid in exactdist"
           print *,"imaterial= ",imaterial
           stop
          endif
         endif  ! imaterial=1,2

         ! drop falling on ice (exactdist)
        else if ((num_materials.eq.3).and.(axis_dir.eq.1)) then
         dist=-9999.0d0
        else if ((num_materials.eq.4).and.(axis_dir.eq.1)) then
         dist=-9999.0d0
        else
         dist=-9999.0d0
        endif  ! drop on slope problem

       else if (SDIM.eq.3) then

        if ((num_materials.eq.3).and. &
            (im_solid_exactdist.eq.3).and. &
            (axis_dir.eq.0).and. &
            (radblob5.eq.zero).and. &
            (radblob6.eq.zero).and. &
            (radblob7.eq.zero).and. &
            (abs(xblob-xblob2).lt.1.0E-6).and. &
            (abs(yblob-yblob2).lt.1.0E-6).and. &
            (abs(zblob-zblob2).lt.1.0E-6)) then

         if ((imaterial.eq.1).or.(imaterial.eq.2)) then
          ! distsolid>0 in solid
          call materialdist(xsten,nhalf,dx,bfact,distsolid,3,time)
          ! in: exactdist (maxtall=2 * radblob => dist_ice=dist_liquid)
          diamblob=two*radblob
          call drop_slope_dist(xstar,ystar,zstar, &
           time,diamblob,dist_ice,dist_liquid)
          distgas=-dist_liquid

          if (imaterial.eq.1) then
           dist=dist_liquid
          else if (imaterial.eq.2) then
           dist=distgas
          else
           print *,"imaterial invalid"
           stop
          endif
         endif  ! imaterial=1,2
        else
         dist=-9999.0d0
        endif  ! drop on slope problem

       else
        print *,"dimension bust for probtype==55 case"
        stop
       endif

      else 
       call materialdist(xsten,nhalf,dx,bfact,dist,imaterial,time)
      endif

      return
      end subroutine exactdist
        
        ! xsten0,nhalf0 corresponds to top level cell
        ! xsten,nhalf center of cell in question
        ! error: f1,e1,f2,e2,f3,e3,f4,e4, ....
        ! called from: stackerror
      subroutine get_symmetric_error( &
        xtrilist, &
        xsten0,nhalf0, & ! top level cell
        dx,bfact, &
        xsten,nhalf, &   ! refined cell
        mofdata, &
        mofdata_tess, &
        errorparm,cutflag,time)

      use MOF_routines_module
      use geometry_intersect_module

      IMPLICIT NONE


      integer, INTENT(in) :: nhalf0,nhalf,bfact
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: mofdata(num_materials*ngeom_recon)
      real(amrex_real), INTENT(in) :: mofdata_tess(num_materials*ngeom_recon)
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,SDIM)  ! top level
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)     ! refined
      integer, parameter :: nhalf_test=1
      real(amrex_real) :: xsten_test(-nhalf_test:nhalf_test,SDIM)
      integer isten
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) dxlevel(SDIM)
      real(amrex_real), INTENT(out) :: errorparm(2*num_materials)
      integer, INTENT(out) :: cutflag
      integer imaterial
      integer minusflag(2)
      integer plusflag(2)
      integer i1,j1,k1,k1lo,k1hi,ii,dir
      real(amrex_real) volbox
      real(amrex_real) cenbox(SDIM)
      real(amrex_real) cenallagain(SDIM)
      real(amrex_real) ltest(D_DECL(3,3,3))
      real(amrex_real) lnode(4*(SDIM-1))
      real(amrex_real) volallagain
      real(amrex_real) volcut(num_materials)
      real(amrex_real) cencut(SDIM)
      real(amrex_real) facearea
      ! get_symmetric_error
      real(amrex_real), INTENT(inout) :: xtrilist(SDIM+1,SDIM,POLYGON_LIST_MAX) 
      integer nmax
      integer tessellate
      integer vofcomp
      real(amrex_real) multi_volume(num_materials)
      real(amrex_real) multi_cen(SDIM,num_materials)
      integer combine_materials,imat1a,imat1b,imaterial_temp
      real(amrex_real) vfrac,dxmaxREFINE
      real(amrex_real) local_LS
      integer inode
     
      if (nhalf.lt.3) then
       print *,"nhalf invalid get symmetric error"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      combine_materials=0
      imat1a=0
      imat1b=0

      nmax=POLYGON_LIST_MAX ! get_symmetric_error

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
       print *,"levelrz invalid get_symmetric error"
       stop
      endif

      call Box_volumeFAST(bfact,dx,xsten,nhalf,volbox,cenbox,SDIM)

      cutflag=0
      do ii=1,2*num_materials
       errorparm(ii)=zero
      enddo

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

      dxmaxREFINE=zero
      do dir=1,SDIM
       dxlevel(dir)=xsten(1,dir)-xsten(-1,dir)
       if (dxlevel(dir).gt.zero) then
        ! do nothing
       else
        print *,"dxlevel invalid"
        stop
       endif
       if (dxlevel(dir).gt.dxmaxREFINE) then
        dxmaxREFINE=dxlevel(dir)
       endif
      enddo ! dir
      if (dxmaxREFINE.gt.zero) then
       ! do nothing
      else
       print *,"dxmaxREFINE invalid"
       stop
      endif

      do imaterial=1,num_materials

       do ii=1,2
        minusflag(ii)=0
        plusflag(ii)=0
       enddo

       do k1=k1lo,k1hi
       do j1=-1,1
       do i1=-1,1
        do isten=-1,1
         dir=1
         xsten_test(isten,dir)=xsten(2*i1+isten,dir)
         dir=2
         xsten_test(isten,dir)=xsten(2*j1+isten,dir)
         if (SDIM.eq.3) then
          dir=SDIM
          xsten_test(isten,dir)=xsten(2*k1+isten,dir)
         endif
        enddo ! isten

        imaterial_temp=imaterial
        if (combine_materials.eq.1) then
         if (imaterial.eq.imat1b) then
          imaterial_temp=imat1a
         endif
        else if (combine_materials.ne.0) then
         print *,"bust"
         stop
        endif

        call exactdist(xsten_test,nhalf_test,bfact, &
         dx,ltest(D_DECL(i1+2,j1+2,k1+2)), &
         imaterial_temp,time)

        local_LS=ltest(D_DECL(i1+2,j1+2,k1+2))

        if ((local_LS.ge.zero).or.(local_LS.le.zero)) then

         if (local_LS.le.zero) then
          minusflag(1)=1
         endif
         if (local_LS.ge.zero) then
          plusflag(1)=1
         endif
         if (abs(local_LS).le.dxmaxREFINE) then
          minusflag(1)=1
          plusflag(1)=1
         endif
       
        else
         print *,"local_LS is NaN : ",local_LS
         stop
        endif
     
        vofcomp=(imaterial-1)*ngeom_recon+1
        vfrac=mofdata_tess(vofcomp)
        if (combine_materials.eq.1) then
         if (imaterial.eq.imat1a) then 
          vofcomp=(imat1b-1)*ngeom_recon+1
          vfrac=vfrac+mofdata_tess(vofcomp)
         else if (imaterial.eq.imat1b) then
          vofcomp=(imat1a-1)*ngeom_recon+1
          vfrac=vfrac+mofdata_tess(vofcomp)
         endif
        else if (combine_materials.eq.0) then
         ! do nothing
        else 
         print *,"combine_materials invalid"
         stop
        endif

        if (vfrac.lt.VOFTOL) then
         minusflag(2)=1
        else if (vfrac.gt.one-VOFTOL) then
         plusflag(2)=1
        else if ((vfrac.ge.VOFTOL).and.(vfrac.le.one-VOFTOL)) then
         minusflag(2)=1
         plusflag(2)=1
        else
         print *,"vfrac out of range:",vfrac
         stop
        endif
         
       enddo
       enddo
       enddo

       if ((minusflag(1).eq.1).and. &
           (minusflag(2).eq.1).and. &
           (plusflag(1).eq.0).and. &
           (plusflag(2).eq.0)) then
        errorparm(2*imaterial-1)=zero
        errorparm(2*imaterial)=zero
       else if ((minusflag(1).eq.0).and. &
                (minusflag(2).eq.0).and. &
                (plusflag(1).eq.1).and. &
                (plusflag(2).eq.1)) then
        errorparm(2*imaterial-1)=volbox
        errorparm(2*imaterial)=zero
       else if ((minusflag(1).eq.1).and. &
                (minusflag(2).eq.0).and. &
                (plusflag(1).eq.0).and. &
                (plusflag(2).eq.1)) then
        errorparm(2*imaterial-1)=volbox
        errorparm(2*imaterial)=volbox
       else if ((minusflag(1).eq.0).and. &
                (minusflag(2).eq.1).and. &
                (plusflag(1).eq.1).and. &
                (plusflag(2).eq.0)) then
        errorparm(2*imaterial-1)=zero
        errorparm(2*imaterial)=volbox
       else 
        cutflag=1
       endif

      enddo ! imaterial
      
      if (cutflag.eq.1) then

        ! first initialize the volumes for the expected solution
       do imaterial=1,num_materials

        inode=1
        do k1=k1lo,k1hi,2
        do j1=-1,1,2
        do i1=-1,1,2

         do isten=-1,1
          dir=1
          xsten_test(isten,dir)=xsten(i1+isten,dir)
          dir=2
          xsten_test(isten,dir)=xsten(j1+isten,dir)
          if (SDIM.eq.3) then
           dir=SDIM
           xsten_test(isten,dir)=xsten(k1+isten,dir)
          endif
         enddo ! isten

         call exactdist(xsten_test,nhalf_test,bfact, &
          dx,lnode(inode), &
          imaterial,time)

         inode=inode+1
        enddo
        enddo
        enddo
 
        call fast_cell_intersection_grid( &
         bfact,dx,xsten,nhalf, &
         lnode, &
         volcut(imaterial), &
         cencut,facearea, &
         volallagain,cenallagain,SDIM)

       enddo ! imaterial

       tessellate=1
       if (nmax.lt.10) then
        print *,"nmax bust 1"
        stop
       endif
        ! in: get_symmetric_error
        ! EPS2
       call multi_get_volume_grid_simple( &
        EPS2, &
        tessellate, &  ! =1
        bfact,dx,xsten0,nhalf0, &
        mofdata, &
        xsten,nhalf, &
        multi_volume,multi_cen, &
        xtrilist, &
        nmax, &
        nmax, &
        SDIM)

       if (combine_materials.eq.1) then
        volcut(imat1b)=volcut(imat1a)
        multi_volume(imat1a)=multi_volume(imat1a)+multi_volume(imat1b)
        multi_volume(imat1b)=multi_volume(imat1a)
       endif 

       do imaterial=1,num_materials
        errorparm(2*imaterial-1)=multi_volume(imaterial)
        errorparm(2*imaterial)=abs(volcut(imaterial)-multi_volume(imaterial))
       enddo

      endif ! cutflag=1
 
      return
      end subroutine get_symmetric_error

        ! xsten0,nhalf0 corresponds to top level cell
        ! xsten,nhalf center of cell in question
        ! error: f1,e1,f2,e2,f3,e3,f4,e4
        ! called from: SUMMASS (NAVIERSTOKES_3D.F90)
      recursive subroutine stackerror( &
       xtrilist, &
       xsten0,nhalf0, &
       dxin,bfact, &
       xsten,nhalf, &
       mofdata, &
       mofdata_tess, &
       errorparm,level, &
       max_level,time)
      IMPLICIT NONE


      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: nhalf0,nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten0(-nhalf0:nhalf0,SDIM)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: dxin(SDIM)
      real(amrex_real), INTENT(inout) :: errorparm(2*num_materials)
      integer, INTENT(inout) :: max_level
      integer, INTENT(in) :: level
      real(amrex_real), INTENT(inout) :: mofdata(num_materials*ngeom_recon)
      real(amrex_real), INTENT(inout) :: mofdata_tess(num_materials*ngeom_recon)
       ! in: stackerror
      real(amrex_real), INTENT(inout) :: xtrilist(SDIM+1,SDIM,POLYGON_LIST_MAX) 

      integer i1,j1,k1,k1lo,k1hi,dir,cutflag,im,isten
      real(amrex_real), allocatable, dimension(:) :: localerror
      real(amrex_real), allocatable, dimension(:,:) :: xstensub
      real(amrex_real), allocatable, dimension(:) :: dxsub
      real(amrex_real), allocatable, dimension(:) :: xmid

      allocate(localerror(2*num_materials)) 
      allocate(dxsub(SDIM))
      allocate(xmid(SDIM))
      allocate(xstensub(-nhalf:nhalf,SDIM))

      if (nhalf.lt.1) then
       print *,"nhalf invalid stackerror"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=0
       k1hi=1
      else
       print *,"dimension bust"
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
       print *,"levelrz invalid stackerror"
       stop
      endif


      if (level.eq.0) then
       do im=1,2*num_materials
        errorparm(im)=zero
       enddo
      endif

      call get_symmetric_error( &
        xtrilist, &
        xsten0,nhalf0, &
        dxin,bfact, &
        xsten,nhalf, &
        mofdata, &
        mofdata_tess, &
        localerror,cutflag,time)

      if ((level.eq.max_level).or.(cutflag.eq.0)) then
       do im=1,2*num_materials
        errorparm(im)=errorparm(im)+localerror(im)
       enddo
      else
       do dir=1,SDIM
        dxsub(dir)=half*(xsten(1,dir)-xsten(-1,dir))
        xmid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        if (dxsub(dir).gt.zero) then
         ! do nothing
        else
         print *,"dxsub invalid"
         stop
        endif
       enddo
       do k1=k1lo,k1hi
       do j1=0,1
       do i1=0,1
        do isten=-nhalf,nhalf
         dir=1
         xstensub(isten,dir)=xmid(dir)+(i1-half)*dxsub(dir)+ &
          half*isten*dxsub(dir)
         dir=2
         xstensub(isten,dir)=xmid(dir)+(j1-half)*dxsub(dir)+ &
          half*isten*dxsub(dir)
         if (SDIM.eq.3) then
          dir=SDIM
          xstensub(isten,dir)=xmid(dir)+(k1-half)*dxsub(dir)+ &
            half*isten*dxsub(dir)
         endif
        enddo !isten
    
        call stackerror( &
         xtrilist, &
         xsten0,nhalf0,dxin,bfact, &
         xstensub,nhalf, &
         mofdata, &
         mofdata_tess, &
         errorparm,level+1,max_level,time)
       enddo
       enddo 
       enddo 
      endif

      deallocate(localerror)
      deallocate(dxsub)
      deallocate(xmid)
      deallocate(xstensub)

      return
      end subroutine stackerror

      subroutine rfile_dist( x,y,dist,rfilemax )
      IMPLICIT NONE


!     * j= 0 ... = rfilemax.

      real(amrex_real) x,y,dist

      integer  j,distset,rfilemax 
      real(amrex_real) deltaz,hval,rval,zval,z1,z2,r1,r2,rmid
      real(amrex_real) distmin,rcoeff,zcoeff,ccoeff,dist1,zint,rint
      real(amrex_real) dist2,dist3,rmin,rmax

      hval=zstatic(rfilemax)
      deltaz=zstatic(1)-zstatic(0)
      rval=x
      zval=yblob-y

      distset=0
      distmin=zero
      do j=1,rfilemax
       z1=zstatic(j-1)
       z2=zstatic(j)
       r1=rstatic(j-1)
       r2=rstatic(j)
       if (r1.lt.r2) then
        rmin=r1
        rmax=r2
       else
        rmin=r2
        rmax=r1
       endif
       if (abs(r2-r1).gt.abs(z2-z1)) then
        rcoeff=(z2-z1)/(r2-r1)
        zcoeff=-one
       else
        rcoeff=-one
        zcoeff=(r2-r1)/(z2-z1)
       endif
       dist=sqrt(rcoeff**2 + zcoeff**2)
       rcoeff=rcoeff/dist
       zcoeff=zcoeff/dist
       ccoeff=-rcoeff*r1-zcoeff*z1
       dist1=rcoeff*rval+zcoeff*zval+ccoeff
       rint=rval-dist1*rcoeff
       zint=zval-dist1*zcoeff
       if ((zint.lt.z1-EPS10).or.(zint.gt.z2+EPS10).or. &
           (rint.lt.rmin-EPS10).or.(rint.gt.rmax+EPS10)) then
        dist2=sqrt( (rval-r1)**2 + (zval-z1)**2 )
        dist3=sqrt( (rval-r2)**2 + (zval-z2)**2 )
        if (dist2.le.dist3) then
         dist1=dist2
        else
         dist1=dist3
        endif
       endif
       if ((abs(dist1).lt.distmin).or.(distset.eq.0)) then
        distmin=abs(dist1)
        distset=1
       endif
      enddo
      
      if (zval.ge.hval) then
       dist=distmin  
      else if (zval.le.zero) then
       dist=rval-rstatic(0)
      else
       j=int( zval/deltaz )
       z1=j*deltaz
       r1=rstatic(j)
       r2=rstatic(j+1) 
       rmid=r1+(zval-z1)*(r2-r1)/deltaz 
       dist=rval-rmid
       if (dist.lt.zero) then
        dist=-distmin
       else
        dist=distmin
       endif
      endif
        
      return
      end subroutine

      subroutine readrfile( rfile, error, rfilemax )
      IMPLICIT NONE


!     *  0..rfilemax

!     * Variables passed in ...

      character rfile*20
      integer   error,rfilemax 
      
      integer  j         

      print *,"will assume rfilemax=",rfilemax
 
!     open(7, file=rfile, form="formatted", status="old", err=900)

      do j = 0,rfilemax 
         
!       read(7,100) zstatic(j),rstatic(j)
        
      end do
 
!     close(7)
 
      do j=0,rfilemax
       print *,"zstatic,rstatic ",zstatic(j),rstatic(j)
      enddo

      error = 0
      
      return
      
! 100  format(' ', e20.14, e21.14 )
! 900  write(6,910) rfile
! 910  format(" Can't open file = ",a," stopping ...")

      error = -1
      
      return
      end subroutine readrfile


       ! called when temperature is prescribed in the solid. (solid distance
       ! function has the appropriate sign)
       ! this routine is called from:
       ! 1. "initsolidtemp"
       !    "initsolidtemp" is called from "solid_temperature()"
       !    when solidheat_flag=1 or 2.
       ! 2. renormalize routine (when solidheat_flag<>0)
       ! 3. vfrac_split (when solidheat_flag<>0)
      subroutine tempsolid(x,y,z,temp,time,im)
      use global_utility_module
      use probcommon_module
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module

      IMPLICIT NONE
     
      real(amrex_real), INTENT(in) :: x,y,z
      real(amrex_real), INTENT(out) :: temp
      real(amrex_real), INTENT(in) :: time 
      integer im
      integer im_solid_temp
      real(amrex_real) LS(num_materials)
      real(amrex_real) xvec(SDIM)
      real(amrex_real) STATE(num_materials*num_state_material)
      integer ibase
      integer bcflag

      bcflag=0

      im_solid_temp=im_solid_primary()

      if (is_rigid(im).ne.1) then
       print *,"is_rigid invalid PROB.F90 in tempsolid"
       stop
      endif
      xvec(1)=x      
      xvec(2)=y      
      if (SDIM.eq.3) then
       xvec(SDIM)=z
      endif 

      if ((im_solid_temp.ge.1).and. &
          (im_solid_temp.le.num_materials)) then

       if (is_in_probtype_list().eq.1) then

        call SUB_LS(xvec,time,LS,num_materials)
        call SUB_STATE(xvec,time,LS,STATE, &
                bcflag,num_materials,num_state_material)
        ibase=(im-1)*num_state_material
        temp=STATE(ibase+ENUM_TEMPERATUREVAR+1) 

       else if (probtype.eq.401) then ! helix user defined
        call HELIX_LS(xvec,time,LS)
        call HELIX_STATE(xvec,time,LS,STATE)
        ibase=(im-1)*num_state_material
        temp=STATE(ibase+ENUM_TEMPERATUREVAR+1) 
       else if (probtype.eq.402) then  ! thermal spray
        call TSPRAY_LS(xvec,time,LS)
        call TSPRAY_STATE(xvec,time,LS,STATE)
        ibase=(im-1)*num_state_material
        temp=STATE(ibase+ENUM_TEMPERATUREVAR+1)
       else if (probtype.eq.412) then ! user defined
        call CAV2Dstep_LS(xvec,time,LS)
        call CAV2Dstep_STATE(xvec,time,LS,STATE)
        ibase=(im-1)*num_state_material
        temp=STATE(ibase+ENUM_TEMPERATUREVAR+1)
       else if (probtype.eq.413) then ! Zeyu's gnbc validation case
        call ZEYU_droplet_impact_LS(xvec,time,LS)
        call ZEYU_droplet_impact_STATE(xvec,time,LS,STATE)
        ibase=(im-1)*num_state_material
        temp=STATE(ibase+ENUM_TEMPERATUREVAR+1)
       else if (probtype.eq.311) then ! user defined
        call USERDEF_LS(xvec,time,LS)
        call USERDEF_STATE(xvec,time,LS,STATE)
        ibase=(im-1)*num_state_material
        temp=STATE(ibase+ENUM_TEMPERATUREVAR+1) 
       else
        temp=fort_tempconst(im)
       endif
      else
       print *,"im_solid_temp invalid 4"
       stop
      endif

      return
      end subroutine tempsolid


       ! return Q=-k dT/dx_dir * nstaircase_dir
       ! nstaircase points from fluid to solid
       ! Q>0 if boundary condition is cooling (decreases energy)
       ! dir=0,1,2
      subroutine tempfluxsolid(x,y,z,tempflux,time,dir)
      use global_utility_module

      IMPLICIT NONE

      real(amrex_real) x,y,z,tempflux,time
      integer dir
      integer im_solid_tempflux

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid tempflux solid"
       stop
      endif

      im_solid_tempflux=im_solid_primary()
 
      if ((im_solid_tempflux.ge.1).and. &
          (im_solid_tempflux.le.num_materials)) then
       tempflux=zero
       if ((probtype.eq.32).and.(SDIM.eq.2)) then
        tempflux= &
         -(fort_tempconst(im_solid_tempflux)-fort_tempconst(1))* &
         get_user_heatviscconst(im_solid_tempflux)
       endif
      else
       print *,"im_solid_tempflux invalid 5"
       stop
      endif

      return
      end subroutine tempfluxsolid

      subroutine dumbbelldist(x,y,z,xc,yc,zc,R,wave,eps,dist)
      IMPLICIT NONE


      real(amrex_real) x,y,z,dist
      real(amrex_real) xc,yc,zc,R,wave,eps
      real(amrex_real) z1,r1,cos1

      if(z.le.(zc-wave)) then
        z1 = z-(zc-wave)
        dist = sqrt((x-xc)**2+(y-yc)**2+z1*z1)-R
      else if(z.ge.(zc+wave)) then
        z1 = z-(zc+wave)
        dist = sqrt((x-xc)**2+(y-yc)**2+z1*z1)-R
      else
        cos1 = cos(3.1415926*(z-zc)/wave)
        r1 = R-eps*(1+cos1)
        dist = abs(cos1)*(sqrt((x-xc)**2+(y-yc)**2)-r1)
      endif

      return
      end subroutine dumbbelldist

       ! override_tagflag is called from fort_vfracerror (PROB.F90) if 
       ! level<max_level_for_use.
      subroutine override_tagflag( &
        i,j,k, &
        level,max_level, &
        snew_ptr,lsnew_ptr, &
        xsten,nhalf,time, &
        rflag,tagflag)
      use probcommon_module
      use global_utility_module
      use global_distance_module
      IMPLICIT NONE

      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: level,max_level
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: rflag
      integer, INTENT(inout) :: tagflag
      real(amrex_real), INTENT(in),pointer :: snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),pointer :: lsnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real) radx,radshrink,dist
      real(amrex_real) x,y,z
      real(amrex_real) local_delta(SDIM)
      integer dir

      if (nhalf.lt.3) then
       print *,"nhalf invalid override_tagflag"
       stop
      endif

      if ((level.ge.0).and.(level.lt.max_level)) then
       ! do nothing
      else
       print *,"level and/or max_level invalid"
       print *,"level=",level
       print *,"max_level=",max_level
       stop
      endif

      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)
      do dir=1,SDIM
       local_delta(dir)=xsten(1,dir)-xsten(-1,dir)
       if (local_delta(dir).gt.zero) then
        ! do nothing
       else
        print *,"local_delta invalid override_tagflag"
        stop
       endif
      enddo ! dir=1..sdim

      if (is_in_probtype_list().eq.1) then

       call SUB_OVERRIDE_TAGFLAG( &
         i,j,k, &
         level,max_level, &
         snew_ptr,lsnew_ptr, &
         xsten,nhalf,time, &
         rflag,tagflag)

       ! bubble formation
      else if (probtype.eq.25) then

       if ((axis_dir.gt.0).and.(SDIM.eq.2)) then
        if ( (abs(x-xblob).le.radblob).and.(y.le.zblob) ) then
         tagflag=1
        endif
        if (rflag.eq.one) then
         tagflag=1
        else if (rflag.eq.zero) then
         ! do nothing
        else
         print *,"rflag invalid: ",rflag
         print *,"probtype: ",probtype
         print *,"axis_dir: ",axis_dir
         stop
        endif
       else if ((axis_dir.eq.0).or.(SDIM.eq.3)) then
        ! do nothing
       else
        print *,"invalid axis_dir or SDIM"
        stop
       endif

      else if (probtype.eq.102) then

       if (SDIM.eq.2) then
        radx=radblob3+y*(radblob4-radblob3)/yblob3  
        if ((y.le.yblob3).and. &
            (x.ge.radx-local_delta(1)).and. &
            (x.le.radx+radblob6+local_delta(1))) then
         tagflag=1
        endif
        if ((y.le.yblob+yblob2).and. &
            (x.le.radblob5+local_delta(1))) then
         tagflag=1
        endif
        radshrink=radblob7**2-radblob5**2
        radshrink=sqrt(radshrink)
        radx=radblob4+(y-yblob3)*(radshrink-radblob4)/(probhiy-yblob3)
        if ((y.ge.yblob3).and. &
            (x.ge.radx-local_delta(1)).and. &
            (x.le.radx+radblob6+local_delta(1))) then
         tagflag=1
        endif
       else if (SDIM.eq.3) then
        ! do nothing
       else
        print *,"dimension bust"
        stop
       endif

       ! in override_tagflag 
      else if (probtype.eq.701) then

        ! dist>0 in the airfoil
       call naca_dist(x,y,z,time,dist)
       if (abs(dist).le.local_delta(1)) then
        tagflag=1
       else if (abs(dist).ge.local_delta(1)) then
        ! do nothing
       else
        print *,"dist or local_delta is NaN"
        stop
       endif
       if ((abs(y+0.1d0).le.0.4d0).and. &
           (x.le.3.0d0)) then
        tagflag=1
       endif

      else if (probtype.eq.66) then

       if (SDIM.eq.2) then
        if (xblob2.lt.xblob3) then
         if ((x.lt.xblob2).or. &
             (x.gt.xblob3)) then
          tagflag=0
          rflag=0.0d0
         endif
        endif
       else if (SDIM.eq.3) then
        ! do nothing
       else
        print *,"dimension bust"
        stop
       endif

      else if (probtype.ge.0) then

       ! do nothing

      else
       print *,"expecting probtype>=0: ",probtype
       stop
      endif 

      return
      end subroutine override_tagflag

       ! called from:
       !  subroutine mask_velocity
       !  subroutine fort_initdatasolid
      subroutine velsolid(x,y,z,vel,time,im,dx)
      use global_utility_module
      use global_distance_module
      use probcommon_module
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      IMPLICIT NONE

      integer, INTENT(in) :: im
      integer dir
      real(amrex_real), INTENT(in) :: x,y,z,time
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(out) :: vel(SDIM)
      real(amrex_real) areacross
      real(amrex_real) tadv,dist
      real(amrex_real) LS(num_materials)
      real(amrex_real) xvec(SDIM)
      integer velsolid_flag

      velsolid_flag=1

      xvec(1)=x      
      xvec(2)=y      
      if (SDIM.eq.3) then
       xvec(SDIM)=z
      endif      

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in velsolid"
      else if (time.lt.zero) then
       print *,"time invalid in velsolid"
       stop
      else
       print *,"time bust in velsolid"
       stop
      endif

      if ((adv_dir.lt.1).or.(adv_dir.gt.2*SDIM+1)) then
       print *,"adv_dir invalid velsolid (6)"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid76"
       stop
      endif

      do dir=1,SDIM
       vel(dir)=zero
      enddo

      if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then

       if (is_in_probtype_list().eq.1) then

         ! pass dx
        call SUB_LS(xvec,time,LS,num_materials)
        call SUB_VEL(xvec,time,LS,vel,velsolid_flag,dx, &
                num_materials)

       else if (probtype.eq.401) then
        call HELIX_LS(xvec,time,LS)
        call HELIX_VEL(xvec,time,LS,vel,velsolid_flag)

       else if (probtype.eq.402) then ! thermal spray
        call TSPRAY_LS(xvec,time,LS)
        call TSPRAY_VEL(xvec,time,LS,vel,velsolid_flag)

       else if (probtype.eq.412) then ! step
        call CAV2Dstep_LS(xvec,time,LS)
        call CAV2Dstep_VEL(xvec,time,LS,vel,velsolid_flag)

       else if (probtype.eq.413) then ! ZEYU droplet impact
         ! pass dx
        call ZEYU_droplet_impact_LS(xvec,time,LS)
        call ZEYU_droplet_impact_LS_VEL(xvec,time,LS,vel,velsolid_flag,dx)
       else if (probtype.eq.311) then ! user defined
        call USERDEF_LS(xvec,time,LS)
        call USERDEF_VEL(xvec,time,LS,vel,velsolid_flag)

        ! cavitation
       else if ((probtype.eq.46).and.(axis_dir.eq.10)) then
          ! dist>0 in the steel sphere
        call stainless_steel_dist_rate(x,y,z,time, &
         dist,vel(SDIM))
          ! CODY ESTEBE created test problem.
       else if ((probtype.eq.46).and.(axis_dir.eq.20)) then
        do dir=1,SDIM
         vel(dir)=zero
        enddo
       else if (probtype.eq.701) then 
        call naca_velocity(x,y,z,time,vel)
       else if (probtype.eq.5700) then ! microfluidics (velsolid)
        if ((axis_dir.eq.4).or.(axis_dir.eq.5)) then
         ! do nothing
        else
         print *,"must have FSI_flag=FSI_PRESCRIBED_PROBF90"
         stop
        endif
        do dir=1,SDIM
         vel(dir)=zero
        enddo
       else if ((probtype.eq.102).and.(SDIM.eq.2)) then
        ! Roper "choked flow" problem
        ! advbot=flow rate L/s=1000 cm^3/s
        ! radblob3= entry gas nozzle radius
        ! radblob5= liquid nozzle outer wall radius

        areacross=Pi*(radblob3**2-radblob5**2)
        if (areacross.le.zero) then
         print *,"cross section area bust"
         stop
        endif
        print *,"radblob7 because what used to be hy here"
        stop
        if ((y.le.radblob7).and.(x.ge.radblob5).and.(x.le.radblob3)) then
         vel(SDIM)=advbot*1000.0/areacross
        endif
       else if ((probtype.eq.32).and.(SDIM.eq.3)) then

        if (advbot.ne.zero) then
         vel(adv_dir)=advbot
        else if ((xblob4.ne.zero).and. &
                 (xblob3.ne.zero)) then
         tadv=time
         vel(adv_dir)= &
           two*Pi*xblob3*sin(two*Pi*tadv/xblob4)/xblob4
        endif

         ! velsolid: pipe velocity=0 at walls
       else if (probtype.eq.41) then

        do dir=1,SDIM
         vel(dir)=zero
        enddo

       else if (probtype.eq.50) then
        print *,"obsolete"
        stop
       else if (probtype.eq.54) then
        vel(1)=zero
        vel(2)=advbot
       else if (probtype.eq.52) then
        print *,"this option obsolete"
        stop
       else if (probtype.eq.56) then
        print *,"obsolete"
        stop
       else if ((probtype.eq.531).and.(SDIM.eq.2)) then ! velsolid
        if (axis_dir.eq.0) then
         vel(SDIM)=advbot
        else if ((axis_dir.eq.1).or.(axis_dir.eq.2).or. &
                 (axis_dir.eq.3)) then
         ! do nothing
        else 
         print *,"axis_dir invalid probtype=531"
         stop
        endif

        ! melting block of ice
       else if (probtype.eq.59) then
        do dir=1,SDIM
         vel(dir)=zero
        enddo
       else if (probtype.eq.32) then  ! flow past moving cylinder
        if (advbot.ne.zero) then
         vel(adv_dir)=advbot
        else if ((xblob4.ne.zero).and. &
                 (xblob3.ne.zero)) then
         tadv=time
         vel(adv_dir)= &
           two*Pi*xblob3*sin(two*Pi*tadv/xblob4)/xblob4
        endif
       else if (probtype.eq.bubbleInPackedColumn) then 
        do dir=1,SDIM
         vel(dir)=zero
        enddo
       endif
      else if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
               (FSI_flag(im).eq.FSI_SHOELE_CTML)) then 

! as simulation progresses: FSI_MF multifab copied to fortran.
! closest value(s) on same processor are used.

       do dir=1,SDIM
        vel(dir)=zero
       enddo

      else
       print *,"FSI_flag invalid in velsolid"
       print *,"im,FSI_flag(im) ",im,FSI_flag(im)
       stop
      endif

      return
      end subroutine velsolid



      subroutine get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)
      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real),  INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, parameter :: nhalf2=1
      real(amrex_real) xsten2(-nhalf2:nhalf2,SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(out) :: cenbc(num_materials,SDIM)
      real(amrex_real), INTENT(out) :: vfrac(num_materials)
      integer im

      integer dir2,i1,j1,k1,k1lo,k1hi,isten
      real(amrex_real) centroid(num_materials,SDIM)
      real(amrex_real) lsgrid(D_DECL(3,3,3),num_materials)
      real(amrex_real), dimension(:), allocatable :: distbatch
      real(amrex_real) facearea(num_materials)
      real(amrex_real) EBVOFTOL
      integer im_solid_microfluidic
      real(amrex_real) initial_time

      initial_time=zero

      im_solid_microfluidic=im_solid_primary()

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid get microfluidic vfrac"
       stop
      endif

      allocate(distbatch(num_materials))

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

      if (probtype.ne.5700) then  ! get_microfluidic_vfrac
       print *,"probtype invalid get_microfluidic_vfrac"
       stop
      endif

      if ( ((xsten(0,SDIM).gt.probhiz+VOFTOL*dx(SDIM)).or. &
            (xsten(0,SDIM).lt.probloz-VOFTOL*dx(SDIM))).and. &
           (SDIM.eq.3)) then

       do im=1,num_materials
        do dir2=1,SDIM
         cenbc(im,dir2)=zero
        enddo
        vfrac(im)=zero
       enddo ! im
       if ((im_solid_microfluidic.lt.1).or. &
           (im_solid_microfluidic.gt.num_materials)) then
        print *,"im_solid_microfluidic invalid 7"
        stop
       endif
       vfrac(im_solid_microfluidic)=one

      else

       do k1=k1lo,k1hi
       do j1=-1,1
       do i1=-1,1
        do isten=-1,1
         dir2=1
         xsten2(isten,dir2)=xsten(isten+2*i1,dir2)
         dir2=2
         xsten2(isten,dir2)=xsten(isten+2*j1,dir2)
         if (SDIM.eq.3) then
          dir2=SDIM
          xsten2(isten,dir2)=xsten(isten+2*k1,dir2)
         endif
        enddo ! isten

        call materialdist_batch( &
         xsten2,nhalf2,dx,bfact, &
         distbatch,initial_time)
        do im=1,num_materials
         lsgrid(D_DECL(i1+2,j1+2,k1+2),im)=distbatch(im) 
        enddo
       enddo
       enddo
       enddo

       EBVOFTOL=VOFTOL
       call getvolumebatch(bfact,dx,xsten,nhalf, &
        lsgrid,vfrac,facearea, &
        centroid,EBVOFTOL,SDIM)

       do im=1,num_materials
        do dir2=1,SDIM
         cenbc(im,dir2)=centroid(im,dir2)-xsten(0,dir2)
        enddo
       enddo ! im

      endif ! not above or below the channel.

      deallocate(distbatch)

      return
      end subroutine get_microfluidic_vfrac


      subroutine microfluidic_velbc(xsten,nhalf,dir,side,vel)
      IMPLICIT NONE

      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) dx(SDIM)
      real(amrex_real) x,y,z
      real(amrex_real), INTENT(inout) :: vel(SDIM)
      integer, INTENT(in) :: dir,side
      integer :: dir2,veldir
      real(amrex_real) LX,LY,LZ,xscale,yscale,zscale
      real(amrex_real) dxscale,dyscale,dzscale
      real(amrex_real) xterm,yterm,zterm
      integer iplug

      if (nhalf.lt.1) then
       print *,"nhalf invalid microfluidic velbc"
       stop
      endif 
      x=xsten(0,1) 
      y=xsten(0,2) 
      z=xsten(0,SDIM)
      do dir2=1,SDIM
       dx(dir2)=xsten(1,dir2)-xsten(-1,dir2)
       if (dx(dir2).le.zero) then
        print *,"dx(dir2) invalid"
        stop
       endif
      enddo
 
      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y 7352"
        stop
       endif
      endif
      if (probtype.eq.5700) then  ! microfluidic_velbc

       iplug=0

       do dir2=1,SDIM
        vel(dir2)=zero
       enddo

       if ((dir.eq.1).and.(side.eq.1)) then ! xlo

        veldir=1

        if (SDIM.eq.2) then

        if (axis_dir.eq.2) then
         print *,"axis_dir=2 obsolete"
         stop
        else if (axis_dir.eq.5) then
         ! do nothing, no inflow from xlo for squeezing geom.
        else if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
                 (axis_dir.eq.3)) then
         vel(veldir)=abs(vinletgas)   ! plug flow
        else if (axis_dir.eq.4) then
         if (zblob2.ge.zblob3) then
          print *,"zblob2 or zblob3 invalid"
          stop
         endif
          ! microfluidic probtype=5700
         if ((yblob3.eq.zero).and.(yblob2.eq.zero)) then
          vel(veldir)=abs(vinletgas)   ! plug flow
         else if (yblob2.lt.yblob3) then
          if ((y.le.yblob2).or.(y.ge.yblob3)) then
           ! do nothing
          else
           if (iplug.eq.1) then
            vel(veldir)=abs(vinletgas)   ! plug flow
           else if (iplug.eq.0) then
            LY=yblob3-yblob2
            yscale=two*((y-yblob2)/LY-half)
            dyscale=dx(2)/LY
            yterm=one-yscale**2-dyscale**2/3.0
            if (yterm.lt.zero) then
             yterm=zero
            endif
            vel(veldir)=abs(vinletgas)*(3.0/2.0)*yterm
           else
            print *,"iplug invalid"
            stop
           endif
          endif
         else
          print *,"inflow dimensions invalid"
          stop
         endif
        else
         print *,"axis_dir invalid microfluidic velbc"
         stop
        endif

        else if (SDIM.eq.3) then

        if (axis_dir.eq.2) then
         print *,"axis_dir=2 obsolete"
         stop
        else if (axis_dir.eq.5) then
         ! do nothing, no inflow from xlo for squeezing geom.
        else if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
                 (axis_dir.eq.3)) then
         vel(veldir)=abs(vinletgas)   ! plug flow
        else if (axis_dir.eq.4) then
         if (zblob2.ge.zblob3) then
          print *,"zblob2 or zblob3 invalid"
          stop
         endif
         if ((zblob3.eq.zero).and.(zblob2.eq.zero)) then
          vel(veldir)=abs(vinletgas)   ! plug flow
         else if ((yblob2.lt.yblob3).and.(zblob2.lt.zblob3)) then
          if ((y.le.yblob2).or.(y.ge.yblob3).or. &
              (z.le.zblob2).or.(z.ge.zblob3)) then
           ! do nothing
          else
           if (iplug.eq.1) then
            vel(veldir)=abs(vinletgas)   ! plug flow
           else if (iplug.eq.0) then
            LY=yblob3-yblob2
            LZ=zblob3-zblob2
            yscale=two*((y-yblob2)/LY-half)
            zscale=two*((z-zblob2)/LZ-half)
            dyscale=dx(2)/LY
            dzscale=dx(SDIM)/LZ
            yterm=one-yscale**2-dyscale**2/3.0
            if (yterm.lt.zero) then
             yterm=zero
            endif
            zterm=one-zscale**2-dzscale**2/3.0
            if (zterm.lt.zero) then
             zterm=zero
            endif
            vel(veldir)=abs(vinletgas)*((3.0/2.0)**2)*yterm*zterm
           else
            print *,"iplug invalid"
            stop
           endif
          endif
         else
          print *,"inflow dimensions invalid"
          stop
         endif
        else
         print *,"axis_dir invalid microfluidic velbc"
         stop
        endif

        else
         print *,"dimension bust"
         stop
        endif

       else if ((dir.eq.1).and.(side.eq.2)) then  ! xhi
        ! do nothing

       else if ((dir.eq.2).and.(side.eq.1)) then  ! ylo

        if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
            (axis_dir.eq.3).or.(axis_dir.eq.4)) then
         ! do nothing
        else if (axis_dir.eq.2) then 
         print *,"this option obsolete"
         stop
        else if (axis_dir.eq.5) then ! squeeze geom, ylo
         veldir=2
         if ((x.ge.xblob2).and. &
             (x.le.xblob3)) then
          vel(veldir)=abs(vinletgas)  ! plug flow
         else
          vel(veldir)=zero
         endif
        else
         print *,"axis_dir invalid ylo velbc"
         stop
        endif

       else if ((dir.eq.2).and.(side.eq.2)) then  ! yhi
        veldir=2

        if (SDIM.eq.2) then

        if (axis_dir.eq.2) then 
         print *,"axis_dir = 2 obsolete"
         stop
        else if (axis_dir.eq.5) then  ! squeezing channel
         if ((x.ge.xblob2).and.(x.le.xblob3)) then
          vel(veldir)=-abs(advbot)  ! plug flow
         else
          vel(veldir)=zero
         endif
        else if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
                 (axis_dir.eq.3)) then
         vel(veldir)=-abs(advbot)  !plug flow
        else if (axis_dir.eq.4) then
         if ((yblob3.eq.zero).and.(yblob2.eq.zero)) then
          vel(veldir)=-abs(advbot)  ! plug flow
         else if (xblob2.lt.xblob3) then

          if ((x.le.xblob2).or.(x.ge.xblob3)) then
           ! do nothing
          else
           if (iplug.eq.1) then
            vel(veldir)=-abs(advbot)  ! plug flow
           else if (iplug.eq.0) then
            LX=xblob3-xblob2
            xscale=two*((x-xblob2)/LX-half)
            dxscale=dx(1)/LX
            xterm=one-xscale**2-dxscale**2/3.0
            if (xterm.lt.zero) then
             xterm=zero
            endif
            vel(veldir)=-abs(advbot)*(3.0/2.0)*xterm
           else
            print *,"iplug invalid"
            stop
           endif
          endif
         else
          print *,"inflow dimensions invalid"
          stop
         endif
        else
         print *,"axis_dir invalid microfluidic velbc"
         stop
        endif

        else if (SDIM.eq.3) then

        if (axis_dir.eq.2) then
         print *,"axis_dir = 2 obsolete"
         stop
        else if (axis_dir.eq.5) then  ! squeezing channel
         if ((x.ge.xblob2).and.(x.le.xblob3)) then
          vel(veldir)=-abs(advbot)  ! plug flow
         else
          vel(veldir)=zero
         endif
        else if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
                 (axis_dir.eq.3)) then
         vel(veldir)=-abs(advbot)  ! plug flow
        else if (axis_dir.eq.4) then
         if ((zblob3.eq.zero).and.(zblob2.eq.zero)) then
          vel(veldir)=-abs(advbot)  ! plug flow
         else if ((xblob2.lt.xblob3).and.(zblob2.lt.zblob3)) then

          if ((x.le.xblob2).or.(x.ge.xblob3).or. &
              (z.le.zblob2).or.(z.ge.zblob3)) then
           ! do nothing
          else
           if (iplug.eq.1) then
            vel(veldir)=-abs(advbot)  ! plug flow
           else if (iplug.eq.0) then
            LX=xblob3-xblob2
            LZ=zblob3-zblob2
            xscale=two*((x-xblob2)/LX-half)
            zscale=two*((z-zblob2)/LZ-half)
            dxscale=dx(1)/LX
            dzscale=dx(SDIM)/LZ
            xterm=one-xscale**2-dxscale**2/3.0
            if (xterm.lt.zero) then
             xterm=zero
            endif
            zterm=one-zscale**2-dzscale**2/3.0
            if (zterm.lt.zero) then
             zterm=zero
            endif
            vel(veldir)=-abs(advbot)*((3.0/2.0)**2)*xterm*zterm 
           else
            print *,"iplug invalid"
            stop
           endif
          endif
         else
          print *,"inflow dimensions invalid"
          stop
         endif
        else
         print *,"axis_dir invalid microfluidic velbc"
         stop
        endif

        else
         print *,"dimension bust"
         stop
        endif

       else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then
        ! do nothing
       else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then
        ! do nothing
       else
        print *,"dir or side invalid"
        stop
       endif

      else
       print *,"probtype invalid in microfluidic velbc"
       stop
      endif

      return
      end subroutine microfluidic_velbc


        ! ice behaves like rigid solid where dist>0
        ! called from "get_icemask_and_icefacecut" (PROB.F90) and 
      subroutine icemask_override(xtarget,im_source,im_dest,dist)
      use global_utility_module
      use global_distance_module

      IMPLICIT NONE
     
      real(amrex_real), INTENT(in) :: xtarget(SDIM) 
      real(amrex_real), INTENT(out) :: dist
      integer, INTENT(in) :: im_source,im_dest

      if ((im_source.lt.1).or.(im_source.gt.num_materials)) then
       print *,"im_source invalid"
       stop
      endif
      if ((im_dest.lt.1).or.(im_dest.gt.num_materials)) then
       print *,"im_dest invalid"
       stop
      endif
      if (is_ice(im_dest).ne.1) then
       print *,"is_ice invalid"
       stop
      endif
      dist=-9999.0d0

      end subroutine icemask_override

      subroutine cavitation_bubble_dist(xsten,nhalf,dist,dx,bfact)
      use global_utility_module
      use global_distance_module
      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) x,y,z
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) xmin,xmax,ymin,ymax,temprad

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid cavitation bubble dist"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        !do nothing
       else
        print *,"expecting z=y 7684"
        stop
       endif
      endif
      if (probtype.ne.46) then
       print *,"probtype invalid"
       stop
      endif

      if (1.eq.0) then

       dist=radblob-sqrt( (x-xblob)**2 + (y-yblob)**2 )

      else

       if (radblob2.gt.zero) then
        xmin=xblob-radblob3-radblob2+half*radblob6
        xmax=xblob+radblob3+radblob2-half*radblob6
        ymin=yblob-radblob4-radblob5+half*radblob6
        ymax=yblob+radblob4+radblob5-half*radblob6
          ! negative in the square
        call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
        dist=-dist
       else
        temprad=radblob
        xmin=-temprad
        xmax=temprad
        ymin=yblob-temprad
        ymax=yblob+temprad
          ! negative in the square
        call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
        dist=-dist
       endif

      endif

      end subroutine cavitation_bubble_dist

        ! imaterial = 1..num_materials
        ! liquid,gas,alt,solid
      subroutine materialdist_batch(xsten,nhalf,dx,bfact,dist,time)
      use global_utility_module
      use global_distance_module
      use hydrateReactor_module
      use unimaterialChannel_module
      use River
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use rigid_FSI_module
      use sinking_particle_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: dx(SDIM) 
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(out) :: dist(:)
      real(amrex_real) x,y,z
      integer imaterial
      real(amrex_real) distline
      real(amrex_real) distcircle
      real(amrex_real) distleft,distright
      real(amrex_real) raddist,distfilament,distfilm
      real(amrex_real) distsolid
      real(amrex_real) drat,veltop,velbot,ytop,ybot
      integer im_solid_materialdist
      integer dir
      real(amrex_real) x_in(SDIM)
      real(amrex_real) maxdx
      real(amrex_real) :: box_xlo,box_xhi
      real(amrex_real) :: box_ylo,box_yhi


      im_solid_materialdist=im_solid_primary()

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid materialdistbatch"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in materialdistbatch"
      else if (time.lt.zero) then
       print *,"time invalid in materialdistbatch"
       stop
      else
       print *,"time bust in materialdistbatch"
       stop
      endif

      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)
      do dir=1,SDIM
       x_in(dir)=xsten(0,dir)
      enddo

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y  7792"
        stop
       endif
      endif

      do imaterial=1,num_materials
       dist(imaterial)=-9999.0
      enddo

      distsolid=-9999.0
      do imaterial=1,num_materials
       if (is_rigid(imaterial).eq.1) then
         ! pos in solid, calling from materialdist_batch
        call materialdistsolid(x,y,z,dist(imaterial),time,imaterial)  
        if (dist(imaterial).gt.distsolid) then
         distsolid=dist(imaterial)
        endif
       else if (is_rigid(imaterial).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid PROB.F90"
        stop
       endif
      enddo ! imaterial=1..num_materials

      if (distsolid.ge.zero) then
       if ((im_solid_materialdist.lt.1).or. &
           (im_solid_materialdist.gt.num_materials)) then
        print *,"im_solid_materialdist invalid: ",im_solid_materialdist
        print *,"probtype= ",probtype
        print *,"axis_dir=",axis_dir
        print *,"num_materials=",num_materials
        stop
       endif
      endif
      if (im_solid_materialdist.eq.1) then
       print *,"im_solid_materialdist=1 invalid: ",im_solid_materialdist
       stop
      endif

      if (is_in_probtype_list().eq.1) then
       call SUB_LS(x_in,time,dist,num_materials)
      else if (probtype.eq.401) then
       call HELIX_LS(x_in,time,dist)
      else if (probtype.eq.402) then
       call TSPRAY_LS(x_in,time,dist)
      else if (probtype.eq.412) then ! step
       call CAV2Dstep_LS(x_in,time,dist)
      else if (probtype.eq.413) then ! zeyu
       call ZEYU_droplet_impact_LS(x_in,time,dist)
      else if (probtype.eq.533) then
       call rigid_FSI_LS(x_in,time,dist)
      else if (probtype.eq.534) then
       call sinking_FSI_LS(x_in,time,dist)
      else if (probtype.eq.311) then ! user defined problem
       call USERDEF_LS(x_in,time,dist)

       ! HYDRATE (materialdist_batch)
      else if (probtype.eq.199) then
       if (num_materials.ne.3) then
        print *,"num_materials invalid for hydrate problem"
        stop
       endif
       call INIT_LS_WATER(x,y,z,time,dist(1))
       call INIT_LS_GAS(x,y,z,time,dist(2))
       call INIT_LS_HYDRATE(x,y,z,time,dist(3))
       do imaterial=1,3
        if (is_rigid(imaterial).ne.0) then
         print *,"all hydrate problem materials should be fluids"
         stop
        endif
       enddo

       !in: materialdist_batch
      else if (probtype.eq.220) then
       if (num_materials.ne.3) then
        print *,"num_materials invalid for unimaterial problem"
        stop
       endif
       maxdx=max(dx(1),dx(2))
       call UNIMAT_INIT_LS_MAT(maxdx,dist(1))
       call UNIMAT_INIT_LS_GST(maxdx,dist(2))

       ! melting (materialdist_batch) (initial level set functions)
      else if (probtype.eq.299) then
       call INIT_LS_LIQUID_MELT(x,y,z,time,dist(1))
       call INIT_LS_GAS_MELT(x,y,z,time,dist(2))
       call INIT_LS_SOLID_MELT(x,y,z,time,dist(3))
       do imaterial=1,3
        if (is_rigid(imaterial).ne.0) then
         print *,"all additive manufacturing materials should be fluids"
         stop
        endif
       enddo

       ! melting block of ice (materialdist_batch)
       ! fluids must tessellate the whole domain.
      else if (probtype.eq.59) then

       if (SDIM.eq.2) then
        if (abs(yblob2-(yblob-half*radblob)).gt.EPS2) then
         print *,"bottom of original ice block must coincide w/substrate"
         print *,"probtype=",probtype
         print *,"substrate_height (yblob2) =",yblob2
         print *,"ice_vertical (yblob if 2d)=",yblob
         print *,"radblob=",radblob
         stop
        endif

        !dist<0 inside the square
        !water below the ice
        box_xlo=xblob-half*radblob
        box_xhi=xblob+half*radblob
        box_ylo=-yblob2-radblob3
        box_yhi=yblob2+radblob3
        call squaredist(x,y,box_xlo,box_xhi, &
         box_ylo,box_yhi,dist(1))
        dist(1)=-dist(1)
        !ice
        box_ylo=yblob2+radblob3
        box_yhi=yblob2+radblob
        call squaredist(x,y,box_xlo,box_xhi, &
         box_ylo,box_yhi,dist(3))
        dist(3)=-dist(3)

        !air; dist<0 inside the square
        box_ylo=-yblob2-radblob
        box_yhi=yblob2+radblob
        call squaredist(x,y,box_xlo,box_xhi, &
         box_ylo,box_yhi,dist(2))

       else if (SDIM.eq.3) then
        print *,"3D case not implemented yet"
        stop
       endif

      else if (probtype.eq.301) then
       call INIT_LS_LIQUID_AM(x,y,z,time,dist(1))
       call INIT_LS_GAS_AM(x,y,z,time,dist(2))
       call INIT_LS_SOLID_AM(x,y,z,time,dist(3))
       do imaterial=1,3
        if (is_rigid(imaterial).ne.0) then
         print *,"all additive manufacturing materials should be fluids"
         stop
        endif
       enddo

       ! in: materialdist_batch
       ! cavitation
      else if ((probtype.eq.46).and.(SDIM.eq.2)) then

         ! water, jwl, air, vacuum
       if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
        if (num_materials.lt.3) then
         print *,"num_materials invalid"
         stop
        endif
        call vapordist(xsten,nhalf,dx,bfact,dist(1))  ! water
        call cavitation_bubble_dist(xsten,nhalf,dist(2),dx,bfact) ! jwl
        dist(3)=y-(zblob+yblob) ! air
        if (num_materials.eq.4) then
         dist(num_materials)=-99999.0 ! vacuum
        endif
       else if (axis_dir.eq.10) then ! cavitation due to falling steel ball
        if (num_materials.ne.3) then
         print *,"num_materials invalid: num_materials=",num_materials
         stop
        endif
         ! dist(num_materials)>0 in the steel sphere
        call vapordist(xsten,nhalf,dx,bfact,dist(num_materials)) 
        dist(1)=99999.0 ! water
        dist(2)=-99999.0 ! ambient
        if (im_solid_materialdist.ne.num_materials) then
         print *,"im_solid_materialdist invalid: ",im_solid_materialdist
         stop
        endif
        ! CODY ESTEBE created test problem
       else if (axis_dir.eq.20) then
        dist(1)=99999.0
        dist(2)=-99999.0 ! ambient
       else
        print *,"axis_dir invalid"
        stop
       endif

       ! River
      else if (probtype.eq.209) then
       if (num_materials.ne.2) then
        print *,"num_materials invalid for River problem"
        stop
       endif
       call RiverHeight(x,y,dist(1),axis_dir)
       dist(1)=dist(1)-z
       dist(2)=-dist(1)

      else

       call vapordist(xsten,nhalf,dx,bfact,dist(1)) 
       if (im_solid_materialdist.ne.2) then
        dist(2)=-dist(1)
       endif

         ! pos in solid, calling from materialdist_batch
       if (probtype.eq.531) then
        if (((axis_dir.eq.3).and.(SDIM.eq.2)).or.  & ! barbell swimmer
            ((axis_dir.eq.2).and.(SDIM.eq.2)).or.  & ! ice crystal in drop
            ((axis_dir.eq.1).and.(SDIM.eq.2))) then ! falling solid on pool
         if (num_materials.ne.3) then
          print *,"expecting num_materials=3 if probtype=531"
          stop
         endif
         print *,"FSI algorithm not implemented yet"
         stop
        endif
       endif ! probtype.eq.531 ?

        ! freezing disk: ice, water, air
       if ((probtype.eq.801).and. &
           (num_materials.eq.3).and.(radblob2.gt.zero)) then
        if (axis_dir.eq.3) then
         dist(3)=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob2
        else if (axis_dir.eq.0) then
         dist(3)=sqrt( (x-xblob)**2 )-radblob2
        else
         print *,"axis_dir invalid probtype=801"
         stop
        endif
        if (dist(1).gt.zero) then
         if (dist(3).ge.zero) then
          dist(1)=-dist(3)
         else
          dist(1)=min(dist(1),-dist(3))
         endif
        endif
       endif

        ! material 1 = top half  material 2=inside circle
        ! material 3 = bottom half
        ! liquid lens (materialdistbatch)
       if ((probtype.eq.202).and.(SDIM.eq.2)) then 
        if (num_materials.ne.3) then
         print *,"num_materials should be 3"
         stop
        endif
        distline=zblob2-y  ! negative at top half
        distcircle=dist(1)  ! positive outside circle
        distright=sqrt( (x-xblob-radblob)**2 + (y-yblob)**2 )
        distleft=sqrt( (x-xblob+radblob)**2 + (y-yblob)**2 )

        dist(2)=-distcircle  ! positive in the circle

        if (distline.le.zero) then  ! top half
         if (distcircle.gt.zero) then
          if (distcircle.gt.abs(distline)) then
           dist(1)=abs(distline)
          else
           dist(1)=distcircle
          endif
         else
          dist(1)=distcircle
         endif
         if (distcircle.le.zero) then
          if (distleft.le.distright) then
           dist(3)=-distleft
          else
           dist(3)=-distright
          endif
         else
          dist(3)=distline
         endif
        else if (distline.ge.zero) then ! bottom half
         if (distcircle.gt.zero) then
          if (distcircle.gt.abs(distline)) then
           dist(3)=abs(distline)
          else
           dist(3)=distcircle
          endif
         else
          dist(3)=distcircle
         endif
         if (distcircle.le.zero) then
          if (distleft.le.distright) then
           dist(1)=-distleft
          else
           dist(1)=-distright
          endif
         else
          dist(1)=-distline
         endif
        endif
       endif ! probtype.eq.202

        ! material 2 is the bubble or drop
        ! material 1 is the lower stratified fluid
        ! material 3 is the upper stratified fluid
       if (probtype.eq.201) then ! bubble-stratified (materialdistbatch)
        if (num_materials.ne.3) then
         print *,"num_materials should be 3"
         stop
        endif
        if (SDIM.eq.2) then
         if (yblob2.eq.zblob2) then
          ! do nothing
         else 
          print *,"expecting yblob2==zblob2 in 2D"
          stop
         endif
        else if (SDIM.eq.3) then
         ! check nothing
        else
         print *,"dimension bust"
         stop
        endif
        distline=zblob2-z
        dist(3)=dist(1)
        if (distline.lt.dist(1)) then
         dist(1)=distline
        endif
        if (-distline.lt.dist(3)) then
         dist(3)=-distline
        endif
       endif

        ! Rieber problem (materialdistbatch)
       if ((probtype.eq.540).and.(SDIM.eq.3)) then  

        if ((radblob3.ne.zero).and.(radblob4.ne.zero)) then
         print *,"conflict of parameters 540"
         stop
        else if (radblob4.ne.zero) then
         if (num_materials.lt.3) then
          print *,"num_materials invalid"
          stop
         endif
         raddist=sqrt( (x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2 )
         if (raddist.le.radblob+half*radblob4) then
          distfilament=raddist-radblob
         else if (raddist.ge.radblob+half*radblob4) then
          distfilament=radblob+radblob4-raddist
         else
          print *,"bust"
          stop
         endif
         dist(3)=distfilament

         distfilm=radblob2-z  ! radblob2-z in 3D  radblob2-y in 2D
         if (raddist.le.radblob+radblob4) then
          dist(2)=raddist-radblob-radblob4
         else if (distfilm.ge.zero) then
          dist(2)=-distfilm
         else
          dist(2)=raddist-radblob-radblob4
          if (dist(2).gt.-distfilm) then
           dist(2)=-distfilm
          endif
         endif
        else if (radblob3.ne.zero) then
         if (num_materials.lt.3) then
          print *,"num_materials invalid"
          stop
         endif
         if (radblob2.ne.zero) then
          print *,"cannot have both radblob2 and radblob3 <>0"
          stop
         endif
         ! if radblob3<>0 then: dist(1)>0 in drop, dist(2)<0 in drop
         ! summary for radblob3<>0:
         !   1. material 1: drop
         !   2. material 2: outside drop and above z=radblob3
         !   3. material 3: outside drop and below z=radblob3
         dist(3)=radblob3-z
         if (dist(2).gt.-dist(3)) then
           dist(2)=-dist(3)
         endif
        endif
       endif  ! 3D rieber problem  probtype.eq.540


        ! 2D Rieber problem (materialdistbatch)
       if ((probtype.eq.540).and.(SDIM.eq.2)) then  

        if ((radblob3.ne.zero).and.(radblob4.ne.zero)) then
         print *,"conflict of parameters 540"
         stop
        else if (radblob4.ne.zero) then
         if (num_materials.lt.3) then
          print *,"num_materials invalid"
          stop
         endif
         raddist=sqrt( (x-xblob)**2 + (y-yblob)**2 )
         if (raddist.le.radblob+half*radblob4) then
          distfilament=raddist-radblob
         else if (raddist.ge.radblob+half*radblob4) then
          distfilament=radblob+radblob4-raddist
         else
          print *,"bust"
          stop
         endif
         dist(3)=distfilament

         distfilm=radblob2-y  ! radblob2-z in 3D  radblob2-y in 2D
         if (raddist.le.radblob+radblob4) then
          dist(2)=raddist-radblob-radblob4
         else if (distfilm.ge.zero) then
          dist(2)=-distfilm
         else
          dist(2)=raddist-radblob-radblob4
          if (dist(2).gt.-distfilm) then
           dist(2)=-distfilm
          endif
         endif
        else if (radblob3.ne.zero) then
         if (num_materials.lt.3) then
          print *,"num_materials invalid"
          stop
         endif
         if (radblob2.ne.zero) then
          print *,"cannot have both radblob2 and radblob3 <>0"
          stop
         endif
          ! if radblob3<>0 then dist(2) is the signed distance
          ! to the falling drop.
         dist(3)=radblob3-y
         if (dist(2).gt.-dist(3)) then
           dist(2)=-dist(3)
         endif
        endif
       endif  ! 2D rieber problem probtype.eq.540

       if ((probtype.eq.530).and.(SDIM.eq.3)) then ! impinging jets
        if (axis_dir.eq.1) then  ! impinging unlike jets
         call get_jet_dist(x,y,z,dist)
        else if (axis_dir.eq.0) then
         ! do nothing
        else
         print *,"axis_dir invalid"
         stop
        endif
       endif !probtype=530


! ysl 05/12/14
        ! water (material 3) on top
        ! diesel (material 1) on bottom
       if ((probtype.eq.17).and.(SDIM.eq.3)) then
        if (num_materials.ne.3) then
         print *,"num_materials must be 3"
         stop
        endif
        if (im_solid_materialdist.ne.0) then
         print *,"no solid in this problem"
         stop
        endif
! dist(i)
! if >0: is material i, the nearest distance to another material
! if <0: the shortest distance to material i
! dist(2): air
! oil: 1
! air: 2
! water: 3
! materialdistbatch
        drat=fort_denconst(3)/fort_denconst(1)  ! dentop/denbot
        veltop=-one/(drat+one)
        velbot=drat/(drat+one)
        ytop=yblob+1.25*radblob
        ybot=yblob-radblob+(ytop-yblob-radblob)*velbot/veltop

        if (y.ge.zero) then
         dist(3)=-dist(2)  ! distance to water
         dist(1)=radblob-sqrt((x-xblob)**2+(y-ybot)**2+(z-zblob)**2)
        else
         dist(1)=-dist(2)
         dist(3)=radblob-sqrt((x-xblob)**2+(y-ytop)**2+(z-zblob)**2)
        endif
       endif  ! probtype=17

        ! water (material 3) on top
        ! diesel (material 1) on bottom
       if ((probtype.eq.17).and.(SDIM.eq.2)) then
        if (num_materials.ne.3) then
         print *,"num_materials must be 3"
         stop
        endif
        if (im_solid_materialdist.ne.0) then
         print *,"no solid in this problem"
         stop
        endif
        drat=fort_denconst(3)/fort_denconst(1)  ! dentop/denbot
        veltop=-one/(drat+one)
        velbot=drat/(drat+one)
        ytop=yblob+1.25*radblob
        ybot=yblob-radblob+(ytop-yblob-radblob)*velbot/veltop

        if (y.ge.zero) then
         dist(3)=-dist(2)  
         dist(1)=radblob-sqrt((x-xblob)**2+(y-ybot)**2)
        else
         dist(1)=-dist(2)
         dist(3)=radblob-sqrt((x-xblob)**2+(y-ytop)**2)
        endif
       endif

      endif ! not custom distance functions

      return
      end subroutine materialdist_batch

        ! imaterial = 1,2,3,4
        ! liquid,gas,alt,solid
      subroutine materialdist(xsten,nhalf,dx,bfact,dist,imaterial,time)
      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(out) :: dist
      integer, INTENT(in) :: imaterial
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), dimension(:), allocatable :: distbatch

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      allocate(distbatch(num_materials))
      if ((imaterial.lt.1).or.(imaterial.gt.num_materials)) then
       print *,"imaterial invalid in materialdist"
       print *,"imaterial = ",imaterial
       stop
      endif

      call materialdist_batch(xsten,nhalf,dx,bfact,distbatch,time)
      dist=distbatch(imaterial)

      deallocate(distbatch)

      return
      end subroutine materialdist



      subroutine pulseheight(x,t,ht)
      IMPLICIT NONE
      real(amrex_real) x,t,ht
      real(amrex_real) xprime,cc,xx


      cc=12.0  ! (118.9453125+0.1953125)/10
      xx=x-cc*t
      xprime=sqrt(three*radblob/(four*zblob*zblob*zblob))
      ht=zblob+radblob/(cosh(xx*xprime)**2)

      return
      end subroutine

      subroutine min_jetdist(x,y,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,dist


      if (levelrz.ne.COORDSYS_RZ) then
       print *,"levelrz invalid min jetdist"
       stop
      endif

      if ((xblob3.ge.radblob).and. &
          (xblob3.le.xblob2)) then
       !do nothing
      else
       print *,"xblob3 out of range"
       stop
      endif 
      if (y.le.yblob3) then
       dist=xblob3-abs(x)
      else
       dist=xblob3-sqrt(x*x+(y-yblob3)*(y-yblob3))
      endif

      return
      end subroutine min_jetdist


      subroutine get_jet_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)
      use geometry_intersect_module
      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, parameter :: nhalf2=1
      real(amrex_real) xsten2(-nhalf2:nhalf2,SDIM)

      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(out) :: cenbc(num_materials,SDIM)
      real(amrex_real), INTENT(out) :: vfrac(num_materials)
      integer im

      integer dir2,i1,j1,k1,k1lo,k1hi,isten
      real(amrex_real) centroid(num_materials,SDIM)
      real(amrex_real) lsgrid(D_DECL(3,3,3),num_materials)
      real(amrex_real) facearea(num_materials)
      real(amrex_real) distbatch(num_materials)
      real(amrex_real) EBVOFTOL

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid get jet vfrac"
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

      do k1=k1lo,k1hi
      do j1=-1,1
      do i1=-1,1
       do isten=-1,1
        dir2=1
        xsten2(isten,dir2)=xsten(isten+2*i1,dir2)
        dir2=2
        xsten2(isten,dir2)=xsten(isten+2*j1,dir2)
        if (SDIM.eq.3) then
         dir2=SDIM
         xsten2(isten,dir2)=xsten(isten+2*k1,dir2)
        endif
       enddo ! isten

       call get_jet_dist(xsten2(0,1),xsten2(0,2), &
         xsten2(0,SDIM),distbatch)
       do im=1,num_materials
        lsgrid(D_DECL(i1+2,j1+2,k1+2),im)=distbatch(im) 
       enddo
      enddo
      enddo
      enddo
      EBVOFTOL=VOFTOL
      call getvolumebatch(bfact,dx,xsten,nhalf, &
        lsgrid,vfrac, &
        facearea,centroid,EBVOFTOL,SDIM)
      do im=1,num_materials
       do dir2=1,SDIM
        cenbc(im,dir2)=centroid(im,dir2)-xsten(0,dir2)
       enddo
      enddo

      return
      end subroutine get_jet_vfrac


      subroutine get_initial_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)
      use geometry_intersect_module
      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, parameter :: nhalf2=1
      real(amrex_real) xsten2(-nhalf2:nhalf2,SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(out) :: cenbc(num_materials,SDIM)
      real(amrex_real), INTENT(out) :: vfrac(num_materials)
      real(amrex_real) :: initial_time
      integer im

      integer dir2,i1,j1,k1,k1lo,k1hi,isten
      real(amrex_real) centroid(num_materials,SDIM)
      real(amrex_real) lsgrid(D_DECL(3,3,3),num_materials)
      real(amrex_real) facearea(num_materials)
      real(amrex_real), dimension(:), allocatable :: distbatch
      real(amrex_real) EBVOFTOL
      real(amrex_real) LS_center

      initial_time=zero

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid get initial vfrac"
       stop
      endif
      allocate(distbatch(num_materials))

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

      do k1=k1lo,k1hi
      do j1=-1,1
      do i1=-1,1
       do isten=-1,1
        dir2=1
        xsten2(isten,dir2)=xsten(isten+2*i1,dir2)
        dir2=2
        xsten2(isten,dir2)=xsten(isten+2*j1,dir2)
        if (SDIM.eq.3) then
         dir2=SDIM
         xsten2(isten,dir2)=xsten(isten+2*k1,dir2)
        endif
       enddo ! isten

       call materialdist_batch( &
        xsten2,nhalf2,dx,bfact, &
        distbatch,initial_time)
       do im=1,num_materials
        lsgrid(D_DECL(i1+2,j1+2,k1+2),im)=distbatch(im) 
       enddo
      enddo
      enddo
      enddo ! i1,j1,k1

      EBVOFTOL=VOFTOL
        ! in: MOF.F90
      call getvolumebatch(bfact,dx,xsten,nhalf, &
        lsgrid,vfrac,facearea, &
        centroid,EBVOFTOL,SDIM)
      do im=1,num_materials

       if (vfrac(im).lt.zero) then
        print *,"vfrac invalid in get_initial_vfrac 1"
        stop
       else if (vfrac(im).le.VOFTOL) then
         ! if the interface is linear and LS_center>=0, then
         !  F>=1/2.  If LS_center>=0, and F<1/2 => nonlinear
         !  interface.
        LS_center=lsgrid(D_DECL(2,2,2),im)
        if (LS_center.ge.-VOFTOL*dx(1)) then
         vfrac(im)=VOFTOL_SLOPES
        endif
       else if ((vfrac(im).gt.zero).and. &
                (vfrac(im).le.one+EPS1)) then
        ! do nothing
       else
        print *,"vfrac invalid in get_initial_vfrac 2: ",vfrac(im)
        stop
       endif

       do dir2=1,SDIM
        cenbc(im,dir2)=centroid(im,dir2)-xsten(0,dir2)
       enddo
      enddo
      deallocate(distbatch)

      return
      end subroutine get_initial_vfrac


! called if probtype.eq.53, probtype.eq.532, probtype.eq.538,
! probtype.eq.541 (2D)
! called if probtype.eq.53, probtype.eq.536, probtype.eq.537,
! probtype.eq.530,probtype.eq.532,probtype.eq.538,
! probtype.eq.541 (3D)
      subroutine get_jetbend_velocity(xsten,nhalf,dx,bfact,vel)
      IMPLICIT NONE

      integer, INTENT(in) :: nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) x,y,z
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) cenbc(num_materials,SDIM)
      real(amrex_real), INTENT(inout) :: vel(SDIM)

      integer dir2
      real(amrex_real) vfrac(num_materials)
      real(amrex_real) angle
      real(amrex_real) xrot,yrot,xrot2,yrot2,radrot,xcen,ycen
      real(amrex_real) dnode1,dnode2

      if (nhalf.lt.3) then
       print *,"nhalf invalid get jetbend velocity"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      do dir2=1,SDIM
       vel(dir2)=zero
      enddo

      if (probtype.eq.537) then  ! get_jetbend_velocity

       call get_jet_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)
       if (vfrac(1).gt.zero) then
        vel(SDIM)=advbot
       else
        vel(1)=adv_vel
       endif

      else if (SDIM.eq.2) then

       if ((probtype.eq.53).and.(axis_dir.ne.2)) then
       
        call get_jet_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)  
        if (vfrac(1).gt.zero) then
          ! this is plug flow
         vel(SDIM)=advbot
         if (1.eq.0) then
          ! u=A(1-(r/r0)^2)  uavg=A(r^2/2-r^4/4 (1/r0^2)) |_0^r0 /
          !                       r^2/2 |_0^r0 =
          ! A(r0^2/2-r0^2/4)/r0^2/2 = A/2
          if (abs(x).lt.radblob) then
           vel(SDIM)=two*advbot*(one-(x/radblob)**2) 
          else
           vel(SDIM)=zero
          endif
         endif
        else
         vel(1)=adv_vel
        endif

       else if ((probtype.eq.53).and.(axis_dir.eq.2)) then
         ! do nothing
       else if (probtype.eq.532) then !impinge from the sides get_jetbend_vel

        angle=0.5235988

        call get_jet_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc) 

        if (vfrac(1).gt.zero) then  
         vel(SDIM)=advbot*cos(angle)
         if (x.lt.zero) then
          vel(SDIM-1)=advbot*sin(angle)
         else
          vel(SDIM-1)=-advbot*sin(angle)
         endif
        else
         vel(SDIM)=zero
        endif

        ! 2D diesel injector w/needle
       else if ((probtype.eq.538).or. & ! inputs.injA
                (probtype.eq.541)) then 

        call get_initial_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)
        if (vfrac(1).gt.zero) then
         vel(SDIM)=advbot
        endif

       else if (probtype.eq.539) then ! supnozz - jetbend_vel

        call get_initial_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)
        if (vfrac(1).gt.VOFTOL) then !initial velocity in the gap
         angle=1.07961377
         vel(1)=advbot*cos(angle)
         vel(2)=advbot*sin(angle)
        endif
  
       else
        print *,"probtype invalid in get jetbend vel 2D"
        stop
       endif

      else if (SDIM.eq.3) then

       if (probtype.eq.532) then ! impinge from the sides get_jetbend_vel

        angle=0.5235988

        call get_jet_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc) 

        if (vfrac(1).gt.zero) then  
         vel(SDIM)=advbot*cos(angle)
         if (y.lt.zero) then
          vel(SDIM-1)=advbot*sin(angle)
         else
          vel(SDIM-1)=-advbot*sin(angle)
         endif
        else
         vel(SDIM)=zero
        endif

       else if (probtype.eq.530) then  ! impinge jets

        angle=0.5235988

        if (axis_dir.eq.0) then ! impinge like jets
         radrot=0.2
         xcen=xblob
         ycen=yblob+radrot
         xrot2=sin(xblob2)*radrot+xcen
         yrot2=cos(xblob2)*radrot+ycen
         xrot=xcen-sin(xblob2)*radrot
         yrot=ycen-cos(xblob2)*radrot
         dnode1=sqrt( (x-xrot)**2 + (y-yrot)**2 )
         dnode2=sqrt( (x-xrot2)**2 + (y-yrot2)**2 )

         call get_jet_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc) 

         if (vfrac(1).gt.zero) then  
          vel(SDIM)=advbot*cos(angle)

          if (dnode1.lt.dnode2) then
           vel(SDIM-1)=advbot*sin(angle)*cos(xblob2)
           vel(1)=advbot*sin(angle)*sin(xblob2)
          else
           vel(SDIM-1)=-advbot*sin(angle)*cos(xblob2)
           vel(1)=-advbot*sin(angle)*sin(xblob2)
          endif
         else
          vel(SDIM)=zero
         endif

        else if (axis_dir.eq.1) then  ! impinge unlike jets

         call get_jet_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc) 
         if ((vfrac(1).gt.zero).or.(vfrac(3).gt.zero)) then
          vel(SDIM)=advbot*cos(angle)
          if (vfrac(1).gt.zero) then
           vel(SDIM-1)=advbot*sin(angle)
          endif 
          if (vfrac(3).gt.zero) then
           vel(SDIM-1)=-advbot*sin(angle)
          endif 
         else
          vel(SDIM)=zero
         endif
        else
         print *,"axis_dir invalid"
         stop
        endif

         ! 3D get_jetbend_velocity - not impinging jets options
         ! 530 and 532 should not appear here.
       else if ( ((probtype.eq.53).and.(axis_dir.ne.2)).or. &  
                 (probtype.eq.536)) then

        call get_jet_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc) 
        if (vfrac(1).gt.zero) then
         vel(SDIM)=advbot
        else
         vel(1)=adv_vel
        endif

        ! 3D diesel injector w/needle
       else if ((probtype.eq.538).or. & ! inputs.injA
                (probtype.eq.541)) then 

        call get_initial_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc) 
        if (vfrac(1).gt.zero) then
         vel(SDIM)=advbot
        endif

        ! jetbend with nozzle and pressure bc.
       else if ((probtype.eq.53).and.(axis_dir.eq.2)) then
         ! do nothing
       else
        print *,"probtype invalid in get jetbend vel"
        stop
       endif

      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine get_jetbend_velocity


      ! called when: 
      ! probtype.eq.53, probtype.eq.532, 
      ! probtype.eq.530, probtype.eq.538,
      ! probtype.eq.539, probtype.eq.529,
      ! probtype.eq.531, probtype.eq.536,
      ! probtype.eq.537 
      subroutine get_jet_dist(x,y,z,dist)
      use global_utility_module
      use global_distance_module

      IMPLICIT NONE

      real(amrex_real) dist(num_materials)
      real(amrex_real) x,y,z,zmin,zmax,dist1,dist2
      real(amrex_real) xmin,xmax
      real(amrex_real) distsolid
      integer im
      real(amrex_real) xrot,yrot,xrot2,yrot2,radrot
      real(amrex_real) xcen,ycen,zcen
      real(amrex_real) angle
      real(amrex_real) xcrit,ycrit,zcrit
      real(amrex_real) xcenbase,xcritbase
      real(amrex_real) ycenbase,ycritbase
      real(amrex_real) y2d
      real(amrex_real) initial_time

      integer im_solid_jet


      im_solid_jet=im_solid_primary()

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y 8810"
        stop
       endif
      endif

      if (SDIM.eq.2) then
       y2d=zero
      else
       y2d=y
      endif

      initial_time=zero

      do im=1,num_materials
       dist(im)=-99999.0
      enddo


      distsolid=-9999.0
      do im=1,num_materials
       if (is_rigid(im).eq.1) then
         ! in get_jet_dist; positive in solid.
        call materialdistsolid(x,y2d,z,dist(im),initial_time,im)
        if (dist(im).gt.distsolid) then
         distsolid=dist(im)
        endif
       else if (is_rigid(im).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid PROB.F90"
        stop
       endif
      enddo ! im=1..num_materials
   
      if (probtype.eq.538) then ! get_jet_dist (inputs.injA)

        ! 2D get_jet_dist: diesel injector w/needle 
       if (SDIM.eq.2) then
        zmin=zblob-0.1
        zmax=zblob+zblob2
        xmin=xblob-radblob
        xmax=xblob+radblob
        call squaredist(x,y,xmin,xmax,zmin,zmax,dist(1))
        dist(1)=-dist(1)
        dist(2)=-dist(1)
       else if (SDIM.eq.3) then
        if (1.eq.0) then ! old
         zmin=zblob-0.1
         zmax=zblob+zblob2
        else   ! new, January 2018
         zmin=zblob
         zmax=zblob2
        endif
        call cylinderdist(x,y,z,xblob,yblob,radblob,zmin,zmax,dist(1))
        dist(1)=-dist(1)
        dist(2)=-dist(1)
       else
        print *,"dimension bust"
        stop
       endif
 
      else if (probtype.eq.537) then ! get_jet_dist

       zmin=-1.0e+10
       zmax=zblob

        ! dist>0 outside the cylinder in cylinderdist.
       call cylinderdist(x,y2d,z,xblob,yblob,radblob,zmin,zmax,dist(1))
       dist(1)=-dist(1)
       dist(2)=-dist(1)
       if ((im_solid_jet.lt.1).or. &
           (im_solid_jet.gt.num_materials)) then
        print *,"im_solid_jet invalid 8"
        stop
       endif

      else if (SDIM.eq.2) then

       if (probtype.eq.541) then
        zmin=zblob-0.5
        zmax=zblob+zblob2
        xmin=xblob-radblob
        xmax=xblob+radblob
        call squaredist(x,y,xmin,xmax,zmin,zmax,dist(1))
        dist(1)=-dist(1)
        dist(2)=-dist(1)
       else if (probtype.eq.539) then ! oblique gap - get_jet_dist
         call gapdist(x,y,xblob,radblob,yblob,yblob2,xblob2,dist(1))
         dist(1)=-dist(1)
         dist(2)=-dist(1)

          ! 2D atomization problem
       else if (probtype.eq.53) then

        zmin=zblob-1.0e+10
        zmax=zblob+radblob
        xmin=xblob-radblob
        xmax=xblob+radblob
         ! dist<0 in the square.
        call squaredist(x,y,xmin,xmax,zmin,zmax,dist(1))
        dist(1)=-dist(1)
        dist(2)=-dist(1)

        if ((axis_dir.eq.1).or.(axis_dir.eq.2)) then
         if ((im_solid_jet.lt.1).or. &
             (im_solid_jet.gt.num_materials)) then
          print *,"im_solid_jet invalid 9"
          stop
         endif
        else if (axis_dir.ne.0) then
         print *,"axis_dir invalid probtype=53"
         stop
        endif

       else if (probtype.eq.532) then ! impinge from the sides get_jet_dist

         ! find distance by tracing back along characteristic

        angle=0.5235988

         ! ycrit=y coordinate of point to find distance
        xcrit=abs(x)
        if (xcrit.lt.probhix-radblob) then
         xcrit=probhix-radblob
        endif
        zcrit=y

          ! given point on center axis
        xcen=probhix+radblob
        zcen=zblob

        xcenbase=tan(angle)*zcen+xcen
        xcritbase=tan(angle)*zcrit+xcrit 

        dist1=radblob-sqrt( (xcritbase-xcenbase)**2 )
        if (abs(x).lt.xcrit) then
         if (dist1.ge.zero) then
          dist1=-sqrt( (xcrit-abs(x))**2 )
         else
          dist1=-sqrt( (xcrit-abs(x))**2+dist1**2 )
         endif
        endif
        dist(1)=dist1
        dist(2)=-dist(1)
       else
        print *,"probtype invalid get_jet_dist probtype=",probtype
        stop
       endif

      else if (SDIM.eq.3) then

        ! 3d diesel injector w/needle
       if (probtype.eq.541) then 
        zmin=zblob-0.1
        zmax=zblob+zblob2
        call cylinderdist(x,y,z,xblob,yblob,radblob,zmin,zmax,dist(1))
        dist(1)=-dist(1)
        dist(2)=-dist(1)
       else if (probtype.eq.529) then  ! 3d airblast

        zmin=xblob-5*radblob
        zmax=xblob+2*radblob
        if (levelrz.eq.COORDSYS_CARTESIAN) then
         call cylinderdist(z,y,x,zblob,yblob,radblob,zmin,zmax,dist(1))
         dist(1)=-dist(1)
        else
         print *,"levelrz invalid get jet dist"
         stop
        endif
        dist(2)=-dist(1)
        
          ! 3D atomization
       else if (probtype.eq.53) then  ! get_jet_dist

           ! dist>0 outside the cylinder in cylinderdist
        zmin=zblob-1.0e+10
        zmax=zblob+radblob
        if (levelrz.eq.COORDSYS_CARTESIAN) then
         call cylinderdist(x,y,z,xblob,yblob,radblob,zmin,zmax,dist(1))
         dist(1)=-dist(1)
        else
         print *,"levelrz invalid get jet dist 2"
         stop
        endif 
        dist(2)=-dist(1)

        if (axis_dir.eq.100) then

         ! do nothing; this routine called from materialdist_batch 
         ! which takes into account the nozzle.  
         ! Set default values here.

        else if ((axis_dir.eq.1).or.(axis_dir.eq.2)) then
           ! cylindrical solid nozzle if axis_dir=1
         if ((im_solid_jet.lt.1).or. &
             (im_solid_jet.gt.num_materials)) then
          print *,"im_solid_jet invalid 10"
          stop
         endif
        else if (axis_dir.ne.0) then
         print *,"axis_dir invalid"
         stop
        endif

       else if (probtype.eq.532) then  ! impinge from the sides get_jet_dist
         ! find distance by tracing back along characteristic

        angle=0.5235988

         ! ycrit=y coordinate of point to find distance
        ycrit=abs(y)
        if (ycrit.lt.probhiy-radblob) then
         ycrit=probhiy-radblob
        endif
        zcrit=z 

          ! given point on center axis
        ycen=probhiy+radblob
        zcen=zblob

          ! v=advbot * sin(angle)
          ! w=advbot * cos(angle)
          ! slope of center axis is dz/dy=-cos(angle)/sin(angle)
          ! or dy/dz=-tan(angle)
          ! y=-tan(angle) (z-zcen) + ycen  is equation for center axis.
          ! the equation of the characteristic going through (y*,z*)
          ! is:
          ! y=-tan(angle) (z-z*) + y*
          ! if base is at z=0, then
          ! ybase=-tan(angle)*(-z*)+y*=tan(angle)*(z*)+y*
        ycenbase=tan(angle)*zcen+ycen
        ycritbase=tan(angle)*zcrit+ycrit 

        dist1=radblob-sqrt( x**2 + (ycritbase-ycenbase)**2 )
        if (abs(y).lt.ycrit) then
         if (dist1.ge.zero) then
          dist1=-sqrt( (ycrit-abs(y))**2 )
         else
          dist1=-sqrt( (ycrit-abs(y))**2+dist1**2 )
         endif
        endif
        dist(1)=dist1
        dist(2)=-dist(1)
        
       else if (probtype.eq.530) then  ! impinging jets

        if (axis_dir.eq.0) then
         radrot=0.2
         xcen=xblob
         ycen=yblob+radrot
         xrot2=sin(xblob2)*radrot+xcen
         yrot2=cos(xblob2)*radrot+ycen
         xrot=xcen-sin(xblob2)*radrot
         yrot=ycen-cos(xblob2)*radrot
   
         zmin=zblob-1.0e+10
         zmax=zblob+radblob
         call cylinderdist(x,y,z,xrot,yrot,radblob,zmin,zmax,dist1)
         call cylinderdist(x,y,z,xrot2,yrot2,radblob,zmin,zmax,dist2)
         dist(1) = max(-dist1,-dist2)
         dist(2)=-dist(1)
        else if (axis_dir.eq.1) then ! impinging jets two materials

         if (num_materials.lt.3) then
          print *,"num_materials too small"
          stop
         endif

         zmin=zblob-1.0e+10
         zmax=zblob+radblob

           ! dist>0 outside the cylinder, in cylinderdist
         call cylinderdist(x,y,z,xblob,yblob,radblob,zmin,zmax,dist1)
         yblob2=yblob+0.4
         call cylinderdist(x,y,z,xblob,yblob2,radblob,zmin,zmax,dist2)
         dist(2) = -max(-dist1,-dist2)
         dist(1)=-dist1
         dist(3)=-dist2
        else
         print *,"axis_dir invalid"
         stop
        endif

       else
        print *,"probtype invalid get_jet_dist 3d probtype=",probtype
        stop
       endif

      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine get_jet_dist

      subroutine is_dissolution(iflag)
      IMPLICIT NONE

      integer iflag


      if (SDIM.eq.3) then
       iflag=0
      else if (SDIM.eq.2) then
       if (probtype.eq.802) then
        iflag=1
       else 
        iflag=0
       endif
      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine is_dissolution


      subroutine blob_array_dist(x,y,z,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,z,dist
      real(amrex_real) hugedist
      integer icomp,dir
      real(amrex_real) distarr(10)
      real(amrex_real) rr
      real(amrex_real) xx(SDIM)

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        !do nothing
       else
        print *,"z=y expected if 2D (in blob_array_dist)"
        stop
       endif
      endif
      hugedist=99999.0

      dist=hugedist

      do icomp=1,10
       distarr(icomp)=hugedist
       rr=radblobarr(icomp)
       if (rr.gt.zero) then
        xx(1)=x-xblobarr(icomp)
        xx(2)=y-yblobarr(icomp)
        if (SDIM.eq.3) then
         xx(SDIM)=z-zblobarr(icomp)
        endif
        distarr(icomp)=zero
        do dir=1,SDIM
         distarr(icomp)=distarr(icomp)+xx(dir)**2
        enddo
        distarr(icomp)=sqrt(distarr(icomp))-rr
      
        dist=min(dist,distarr(icomp))
       endif ! rr>0 
      enddo ! icomp

      return
      end subroutine blob_array_dist

! called from:
!  materialdist_batch
!  velbc_override
!  denBC
!
! dist>0 in material 1
! dist<0 in material 2
      subroutine vapordist(xsten,nhalf,dx,bfact,dist)
      use global_utility_module
      use global_distance_module
      use marangoni

      IMPLICIT NONE
      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) x,y,z
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real), INTENT(in) :: dx(SDIM)

      real(amrex_real) NPT,HSB,NID,NOD,CHH,scaleCHH,VRAD,dist1,dist2
      real(amrex_real) xmin,xmax,ymin,ymax,zmin,zmax,zz,temprad
      real(amrex_real) m,b
      integer igeom
      real(amrex_real) costheta,sintheta,xprime,yprime,zprime,delta
      real(amrex_real) h1
      real(amrex_real) hugedist
      real(amrex_real) distbatch(num_materials)
      real(amrex_real) initial_time
      real(amrex_real) ypretend
      real(amrex_real) wave_number
      real(amrex_real) box_ylo
      real(amrex_real) cylinder_zlo
      real(amrex_real) ktermx,velperturbx
      real(amrex_real) drat,veltop,velbot,ytop,ybot,y2d
      real(amrex_real), parameter :: stub_zero=zero

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        !do nothing
       else
        print *,"z=y expected if 2D (in vapordist)"
        stop
       endif
      endif

      if (nhalf.lt.1) then
       print *,"nhalf invalid vapordist"
       stop
      endif

      initial_time=zero

      hugedist=99999.0

      igeom=0

      dist=hugedist

       ! vapordist
      if ((probtype.eq.1).and.(axis_dir.eq.15)) then
       if (SDIM.eq.2) then
        dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
       else
        dist=-sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)+radblob
       endif
       ! vapordist 2d or 3d
       ! dist>0 in material 1
       ! dist<0 in material 2
      else if (probtype.eq.537) then ! vapordist

        zmin=-1.0e+10
        zmax=zblob
        if (SDIM.eq.2) then
         y2d=zero
        else
         y2d=y
        endif
         ! dist>0 outside the cylinder in cylinderdist.
        call cylinderdist(x,y2d,z,xblob,yblob,radblob,zmin,zmax,dist)
        dist=-dist

      else if (probtype.eq.801) then

       ! vapordist:
       ! test problems from Welch and Wilson 2000
       ! vapor on left, water on right.

       if (axis_dir.eq.0) then
        dist=x-xblob
        if (num_materials.eq.3) then
         dist=sqrt( (x-xblob)**2 )-radblob
        endif
       else if (axis_dir.eq.1) then
        dist=y-yblob
       else if ((axis_dir.eq.2).and.(SDIM.eq.3)) then
        dist=z-zblob
       else if (axis_dir.eq.3) then
        if (SDIM.eq.2) then
         dist=sqrt( (x-xblob)**2+(y-yblob)**2 )-radblob
        else if (SDIM.eq.3) then
         dist=sqrt( (x-xblob)**2+(y-yblob)**2+(z-zblob)**2 )-radblob
        else
         print *,"dimension bust"
         stop
        endif
       else
        print *,"axis_dir invalid"
        stop
       endif

      else if (probtype.eq.92) then ! shock tube
       ! do nothing
      else if (probtype.eq.93) then ! shock tube with interface

       dist=xblob-x

      else if (SDIM.eq.2) then

        ! dissolution (vapordist)
       if (probtype.eq.802) then
        ktermx=two*Pi*yblob3*x/(two*radblob)
        velperturbx=one+radblob3*cos(ktermx)
        if (y.gt.yblob) then
         dist=yblob+radblob*velperturbx-y
        else 
         dist=y-(yblob-radblob*velperturbx)
        endif

         ! Rayleigh Taylor and checkerboard test (vapordist)
       else if (probtype.eq.602) then
         ! yblob is the average height of the interface 
         ! between water and air
         ! radblob is the amplitude of the perturbation
         ! problox <= x <= probhix
         ! y=yblob+radblob*cos(2.0*pi*(x-problox)/(probhix-problox))
        if (xblob.lt.zero) then
         print *,"xblob invalid for rayleigh taylor test"
         stop
        else if (xblob.eq.zero) then
         wave_number=one
        else
         wave_number=xblob
        endif
        if ((xblob.ge.1.0D+10).and.(radblob.le.(EPS10))) then
         dist=y-(yblob+half*(1.0d0/16.0d0))
        else 
         dist=y-(yblob+ &
          radblob*cos(two*Pi*wave_number*(x-problox)/(probhix-problox)))
        endif

       else if (probtype.eq.603) then  ! Benard advection (vapordist)
        dist=zblob-y

        ! driven wave problem, water on bottom, air on top.
       else if (probtype.eq.90) then
        dist=yblob-y
       else if (probtype.eq.540) then ! Rieber simulation vapordist 2D
        dist=radblob-sqrt(x**2 + (y-yblob)**2)
       
        if (radblob2.ne.zero) then
         if (radblob3.ne.zero) then
          print *,"cannot have both radblob2 and radblob3 <>0"
          stop
         endif
         dist2=radblob2-y
         if (dist2.gt.dist) then
          dist=dist2
         endif
        endif

! Roper nozzle problem
       else if (probtype.eq.102) then
        dist=half*(radblob5+radblob2)-abs(x)
        dist2=yblob+yblob2-y
        if (dist2.lt.dist) then
         dist=dist2
        endif 
       else if (probtype.eq.3) then
        dist=xblob+radblob*cos(two*Pi*y/yblob)-x
       else if (probtype.eq.4) then
        call rtdist(x,y,dist)
       else if (probtype.eq.7) then
        dist = radblob-sqrt((x-xblob)**2 + (y-yblob)**2)
        dist1 = yblob-two-y
        if (dist1.gt.dist) then
         dist=dist1
        endif
       else if (probtype.eq.8) then
        dist = -radblob+sqrt((x-xblob)**2 + (y-yblob)**2)
        dist1 = yblob+radblob-y+0.2
        if (dist1.lt.dist) then
         dist=dist1
        endif
       else if (probtype.eq.531) then  ! rigid body problem - vapordist
        if ((axis_dir.eq.0).or.(axis_dir.eq.1)) then
         dist=yblob2-y
        else if (axis_dir.eq.2) then
         dist=yblob2-y
         dist2=radblob2-sqrt((x-xblob)**2 + (y-yblob)**2)
         if (dist.lt.dist2) then
          dist=dist2
         endif 
        else if (axis_dir.eq.3) then
         dist=99999.0
        else
         print *,"axis_dir invalid probtype=531"
         stop
        endif

! microfluidics channel -- 0,3,4 Roper, 1 Comsol, 2,5 squeeze vapordist
       else if (probtype.eq.5700) then
        if (axis_dir.eq.2) then
         print *,"axis_dir=2 obsolete"
         stop
        else if (axis_dir.eq.5) then
         dist=y-yblob
        else if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
                 (axis_dir.eq.3).or.(axis_dir.eq.4)) then
         dist=y-yblob
         if (xblob4.gt.zero) then
          if ((y.ge.yblob).and.(x.le.half*xblob4)) then
           if (y-yblob.le.x-half*xblob4) then
            dist=y-yblob
           else
            dist=half*xblob4-x
           endif
          else if ((y.le.yblob).and.(x.le.half*xblob4)) then
           dist=y-yblob
          else if (y.ge.yblob) then
           dist=half*xblob4-x
          else
           dist=-sqrt( (half*xblob4-x)**2+(y-yblob)**2 )
          endif
         endif
        else
         print *,"axis_dir invalid probtype=5700"
         stop
        endif

       else if (probtype.eq.101) then
! dist<0 inside the square
        if (axis_dir.eq.0) then
         call squaredist(x,y,xblob-radblob,xblob+radblob,yblob-radblob, &
          yblob+radblob,dist)
        else if (axis_dir.eq.1) then
         h1=1.10*radblob
         call squaredist(x,y,xblob-radblob-h1,xblob+radblob-h1,yblob-radblob, &
          yblob+radblob,dist)
         call squaredist(x,y,xblob-radblob+h1,xblob+radblob+h1,yblob-radblob, &
          yblob+radblob,dist2)
         if (dist2.lt.dist) then
          dist=dist2
         endif
        endif
       else if (probtype.eq.58) then
        dist=yblob+radblob*cos(two*Pi*x/xblob)-y
       elseif ((probtype.eq.63).or.(probtype.eq.64)) then
         if (y.ge.two*xblob10) then
          dist=half*xblob10-x
         else
          dist=half*xblob10-x+(two*xblob10-y)
         endif
       else if ((probtype.eq.1).and.(axis_dir.eq.13)) then

           ! dist<0 in the blobs.
        call blob_array_dist(x,y,z,dist)

       else if ((probtype.eq.1).and.(axis_dir.eq.11)) then ! vapordist 2d
        dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
       else if ((probtype.eq.1).and.(axis_dir.eq.12)) then ! vapordist 2d
        dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
       else if ((probtype.eq.1).and.(axis_dir.eq.14)) then
        dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
       else if ((probtype.eq.1).and.(axis_dir.eq.140)) then
        dist=max(-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob, &
                 -sqrt( (x-xblob)**2 + (y-yblob2)**2 )+radblob)
       else if ((probtype.eq.1).and.(axis_dir.eq.141)) then
        dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob
       else if ((probtype.eq.1).and.(axis_dir.lt.150)) then
        dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
       else if (probtype.eq.5) then
        call ellipsedist(x,y,radblob,zblob,xblob,yblob,dist)
       else if (probtype.eq.12) then
        call legenddist(x,y,dist)
       else if (probtype.eq.13) then
        call legenddist(x,y,dist)
       else if (probtype.eq.14) then
        dist = sqrt((x-xblob)**2 + (y-yblob)**2)-radblob
       else if (probtype.eq.16) then
         if (x.le.xblob+radblob) then
          dist=yblob-y-EPS3
         else
          dist=sqrt((x-xblob-radblob)*(x-xblob-radblob)+ &
                    (yblob-y)*(yblob-y))
         endif
         ! vapordist
         ! water (material 3) on top
         ! diesel (material 1) on bottom
         ! for probtype=18, same liquid top or bottom
       else if ((probtype.eq.17).or.(probtype.eq.18)) then
        if ((xblob.ne.zero).or.(yblob.ne.zero).or. &
            (radblob.ne.half)) then
         print *,"set xblob=yblob=0 radblob=1/2 drop collide"
         stop
        endif
        if (probtype.eq.18) then
         drat=one
        else
         drat=fort_denconst(3)/fort_denconst(1)  ! dentop/denbot
        endif
        veltop=-one/(drat+one)
        velbot=drat/(drat+one)
        ytop=yblob+1.25*radblob
        ybot=yblob-radblob+(ytop-yblob-radblob)*velbot/veltop

        dist = radblob-sqrt((x-xblob)**2+(y-ytop)**2)  ! water (3,top)
        dist1 = radblob-sqrt((x-xblob)**2+(y-ybot)**2) ! oil (1,bot)
        if (dist1.gt.dist) then
         dist=dist1
        endif
       else if (probtype.eq.23) then
        dist=yblob+radblob*cos(two*Pi*x/xblob)-y
       else if ((probtype.eq.24).or.(probtype.eq.27)) then
        dist=half-y
       else if ((probtype.eq.25).and.(axis_dir.eq.0)) then ! vapordist
        dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
        dist=-dist
        ! vapordist - bubble formation
       else if ((probtype.eq.25).and.(axis_dir.ge.1).and. &
                (axis_dir.le.11)) then
        if (zblob.le.zero) then
         print *,"zblob (nozzle ht) must be positive"
         stop
        else  
          ! dist<0 in the square.
         box_ylo=-1000.0d0*zblob
         call squaredist(x,y,-radblob,radblob,box_ylo,zblob,dist)
        endif
       else if (probtype.eq.51) then  ! vapordist oscillating column
        dist=radblob-abs(x-xblob) 
       else if (probtype.eq.43) then  ! vapordist 
        xmin=-xblob-radblob
        xmax=xblob+radblob
        ymin=-yblob/two
        ymax=yblob/two
        if (zblob.gt.zero) then
         ymin=ymin-zblob/two
         ymax=ymax+zblob/two
        endif
        call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
       else if (probtype.eq.48) then
        dist1=y-yblob
        m=-one
        b=yblob-m*xblob
        dist2=y-(m*x+b)
        dist=dist1
        if (dist.lt.dist2) then
         dist=dist2
        endif
       else if (probtype.eq.11) then
        print *,"cavitation with top wall outflow cond. is deleted"
        stop
       else if (probtype.eq.47) then
        if (radblob.ne.zero) then
          dist=radblob-abs(x-half*xblob)
        else
          print *,"radblob cannot equal zero!"
          stop
        endif
       else if (probtype.eq.35) then
        dist=half*yblob-y
       else if (probtype.eq.2) then
        dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
       else if (probtype.eq.201) then  ! vapordist stratified 2D
        dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
       else if (probtype.eq.202) then  ! vapordist liquid lens
        dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
       else if (probtype.eq.36) then ! vapordist 2D: bubble

         !dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob
        call spheredist(x,y,z,dist)  ! dist<0 in the sphere

        if ((probtype.eq.36).and.(axis_dir.eq.100)) then
         dist=-dist
        endif

        if ((probtype.eq.36).and.(axis_dir.eq.210)) then
         dist=(y+radblob*(x-xblob)-yblob)/sqrt(one+radblob**2)
        endif

        if ((probtype.eq.36).and.(axis_dir.eq.7)) then
         xmin=xblob-radblob
         xmax=xblob+radblob
         ymin=yblob-radblob
         ymax=yblob+radblob
         call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
        endif
         ! marangoni (heat pipe) problem dist<0 in the bubble
        if ((probtype.eq.36).and.(axis_dir.eq.10)) then
         call dist_long_bubble(radblob,radblob2,x,y,z,dist)
        endif

        ! vapordist bubble jetting 2D - dist>0 in the water.
       else if (probtype.eq.42) then 

        dist=sqrt( (x-xblob)**2 + (y-yblob)**2 )-radblob

         ! if free surface is in the domain
        dist1=zblob+yblob-y
        if (dist1.lt.dist) then
         dist=dist1
        endif

        ! vapordist cavitation 2D - dist>0 in the water for jwl case,
        ! dist>0 in the sphere for the sphere impact case.
       else if (probtype.eq.46) then

        if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then 
         call cavitation_bubble_dist(xsten,nhalf,dist,dx,bfact)
         dist=-dist

          ! if free surface is in the domain
         dist1=zblob+yblob-y
         if (dist1.lt.dist) then
          dist=dist1
         endif
        else if (axis_dir.eq.10) then
         dist=radblob-sqrt( (x-xblob)**2 + (y-yblob)**2 )

         ! CODY ESTEBE created test problem
        else if (axis_dir.eq.20) then
         dist=99999.0
        else
         print *,"axis_dir invalid"
         stop
        endif

       else if (probtype.eq.9) then  ! ship wave free surface (vapordist) 2D
        ypretend=zero
        call boatdist(x,ypretend,y,dist)
       else if (probtype.eq.37) then
        dist=radblob-sqrt( (x-xblob)**2 + (y-yblob)**2 )
       else if ((probtype.eq.22).and.(axis_dir.eq.14)) then
        call initjetparms(HSB,NOD,NPT,NID,CHH,scaleCHH)
        VRAD=scaleCHH/eight
        dist=sqrt( (x-zero)**2 + (y-zero)**2 ) - VRAD  
       else if (probtype.eq.39) then
        dist=yblob+radblob*cos(two*Pi*x/xblob)-y
! Suvorov Formation of Taylor cone...
        if (radblob3.gt.zero) then
         dist=0.2*exp(-radblob3*x*x)-y+yblob3
        endif
! pipe - vapordist 2D
       else if (probtype.eq.41) then
        call inletpipedist(x,y,z,distbatch) 
        dist=distbatch(1)
       else if (probtype.eq.44) then
         ! in 2d, y=z
        call damdist(x,z,dist,igeom)
       else if (probtype.eq.45) then
        dist=yblob+1/(2*Pi)*(radblob*cos(2*Pi*(x-0.1)/xblob)+ &
           1/2*radblob**2*cos(4*Pi*(x-0.1)/xblob)+ &
           3/8*radblob**3*cos(6*Pi*(x-0.1)/xblob))-y
       else if (probtype.eq.21) then
        dist=yblob/two-y
       else if (probtype.eq.49) then
        dist=sqrt((x-xblob)**2+(y-yblob)**2) - radblob
        dist1= 1 - (y-yblob)
        if (dist1.lt.dist) then
         dist=dist1
        endif
       else if (probtype.eq.53) then  ! 2D JICF vapordist
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
       else if (probtype.eq.532) then ! impinge from sides (vapordist)
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
         ! 2d diesel injector w/needle (vapordist)
       else if ((probtype.eq.538).or. &  ! inputs.injA
                (probtype.eq.541)) then
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
       else if (probtype.eq.701) then ! flapping wing (vapordist 2D)
        if ((axis_dir.eq.0).or.(axis_dir.eq.1)) then
         dist=zblob-y
        else if (axis_dir.eq.2) then
         dist=hugedist
        else
         print *,"axis_dir invalid"
         stop
        endif
       else if (probtype.eq.539) then ! jet w/ supersonic crossflow - vapordist
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
       else if(probtype.eq.72) then
        xmin=xblob-radblob
        xmax=xblob+radblob
        ymin=zblob-1.0e+10
        ymax=zblob+radblob
        call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
        dist=-dist
       else if (probtype.eq.54) then
        xmin=-xblob
        xmax=xblob
        ymin=-yblob
        ymax=yblob+half*yblob
        call squaredist(x,y,xmin,xmax,ymin,ymax,dist)
        dist=-dist
       else if (probtype.eq.61) then
        dist=-(sqrt((x-xblob)**2+(y-yblob)**2) - radblob)
        if (axis_dir.eq.0) then
         dist1=-(3 + (y-yblob))
        else if (axis_dir.eq.1) then
         dist1=-(radblob + (y-yblob))
        else
         print *,"axis_dir invalid probtype=61"
         stop
        endif
  
        if (dist1.gt.dist) then
         dist=dist1
        endif
       else if (probtype.eq.62) then
        dist=zblob3-y
        costheta=cos(xblob2)
        sintheta=sin(xblob2)
        xprime=costheta*(x-xblob)-sintheta*(y-yblob)
        yprime=sintheta*(x-xblob)+costheta*(y-yblob)
        delta=half*(radblob2-radblob)
        temprad=1.0D+10
        box_ylo=-half*zblob2-delta
        call squaredist(xprime,yprime,-radblob-delta,radblob+delta, &
          box_ylo,temprad,dist1)
        if (dist1.lt.dist) then
         dist=dist1
        endif 
       else if (probtype.eq.50) then
        xmin=-50.0
        xmax=half*twall
        ymin=-50.0
        ymax=50.0
        zmin=-50.0
        zmax=3.8
        zz=yblob
        call cubedist(xmin,xmax,ymin,ymax,zmin,zmax, &
                      x,zz,y,dist)
        dist=-dist

        xmin=-50.0
        xmax=50.0
        ymin=-50.0
        ymax=50.0
        zmin=-50.0
        zmax=0.2
        call cubedist(xmin,xmax,ymin,ymax,zmin,zmax, &
                      x,zz,y,dist1)
        dist1=-dist1

        if (dist1.gt.dist) then
         dist=dist1
        endif

       else if (probtype.eq.52) then
        dist=zblob2-y
       else if (probtype.eq.56) then
        dist=zblob2-y
       else if (probtype.eq.66) then
        call pulseheight(x,stub_zero,dist)
        dist=dist-y
       else if (probtype.eq.38) then
        dist=yblob+radblob*cos(two*Pi*x/xblob)-y
! 2d: vapordist
       else if (probtype.eq.26) then
        if (axis_dir.eq.1) then !swirl with interface
         dist=abs(y-half)-one/four
        else if (axis_dir.eq.0) then ! swirl, no interface
         ! do nothing
        else if (axis_dir.eq.2) then ! inputs.vortex_confine (no interface)
         ! do nothing
        else if (axis_dir.eq.3) then ! vortex_confinement w/interface
         dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
        else if (axis_dir.eq.10) then ! inputs.BCG, Dirichlet BC.
         ! do nothing
        else if (axis_dir.eq.11) then ! inputs.BCG_periodic
         ! do nothing
        else if (axis_dir.eq.12) then ! inputs.buoyancy
         ! do nothing
        else
         print *,"axis_dir 0, 1,2,3,10, 11 or 12 expected"
         stop
        endif
! natural convection in triangular enclosure
       else if (probtype.eq.81) then
        dist=yblob2-y
       else if (probtype.eq.110) then
!       dist=yblob-y
        call local_shallow_water_elevation(initial_time,x,dist)
        dist=dist-y
       endif

      else if (SDIM.eq.3) then

       if (probtype.eq.540) then ! Rieber simulation vapordist 3D
        dist=radblob-sqrt((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)

        if (radblob2.ne.zero) then
         if (radblob3.ne.zero) then
          print *,"cannot have both radblob2 and radblob3 <>0"
          stop
         endif
         dist2=radblob2-z
         if (dist2.gt.dist) then
          dist=dist2
         endif
        endif
       else if (probtype.eq.5501) then  ! rough surface problem.
        if (axis_dir.ne.0) then
         print *,"axis_dir invalid"
         stop
        endif
        dist=radblob-sqrt( (x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2 )

       else if (probtype.eq.4) then
        dist=sqrt((x-(xblob+0.5))**2+(y-yblob)**2+(z-zblob)**2 )-radblob
        dist1=sqrt((x-(xblob-0.5))**2+(y-yblob)**2+(z-(zblob+2.3))**2)- &
            radblob
        if (dist.gt.dist1) then
         dist=dist1
        endif
       else if (probtype.eq.5) then
        dist=sqrt((x-(xblob+radblob/two))**2+(y-yblob)**2+(z-zblob)**2)- &
           one
        dist1=sqrt((x-(xblob-radblob/two))**2+(y-yblob)**2+(z-zblob)**2)- &
           one
        if (dist.gt.dist1) then
         dist=dist1
        endif
       else if (probtype.eq.603) then  ! Benard advection (vapordist)
        dist=zblob-z
       else if (probtype.eq.101) then
        call cubedist(xblob-radblob,xblob+radblob,yblob-radblob,yblob+radblob, &
         zblob-radblob,zblob+radblob,x,y,z,dist)
       else if ((probtype.eq.1).and.(axis_dir.le.10)) then ! vapordist 3d
        dist=sqrt( (x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2 )-radblob
       else if ((probtype.eq.1).and.(axis_dir.eq.11)) then ! vapordist 3d
        dist=sqrt( (x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2 )-radblob
       else if ((probtype.eq.1).and.(axis_dir.eq.13)) then ! vapordist 3d
        ! dist<0 in the blobs.
        call blob_array_dist(x,y,z,dist)
       else if ((probtype.eq.1).and.(axis_dir.eq.12)) then
        dist=-sqrt( (x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2 )+radblob
       else if (probtype.eq.35) then
        dist=half*zblob-z 
       else if (probtype.eq.2) then
        dist=sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)-radblob
       else if (probtype.eq.201) then  ! vapordist stratified 3D
        dist=sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)-radblob
       else if (probtype.eq.390) then ! dumbbell shape (vapordist)
        call dumbbelldist(x,y,z,xblob,yblob,zblob, &
          radblob,zblob2,radblob2,dist)
       else if (probtype.eq.36) then ! vapordist 3D

        call spheredist(x,y,z,dist)  ! dist < 0 in the sphere

        if ((probtype.eq.36).and.(axis_dir.eq.100)) then
         dist=-dist
        endif

        if ((probtype.eq.36).and.(axis_dir.eq.200)) then
         dist=sqrt((y-yblob)**2+(z-zblob)**2)-radblob
        else if ((probtype.eq.36).and.(axis_dir.eq.201)) then
         dist=sqrt((x-xblob)**2+(z-zblob)**2)-radblob
        else if ((probtype.eq.36).and.(axis_dir.eq.202)) then
         dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
        endif

        if ((probtype.eq.36).and.(axis_dir.eq.210)) then
         dist=(z+radblob*(x-xblob)+radblob2*(y-yblob)-zblob)/ &
              sqrt(one+radblob**2+radblob2**2)
        endif

         ! for inputs3d.bubbleburst,
         !  gas bubble submerged in liquid with a horizontal
         !  gas/liquid interface above the bubble.
         !  zblob=vertical coordinate of center of bubble
         !  zblob10=vertical coordinate of the horizontal interface
         !  be careful, zblob10 defaults to 0.
        if ((probtype.eq.36).and. &
            (zblob.gt.zero).and. &
            (zblob10.gt.zblob).and. &
            (probloz.ge.zero)) then
         dist2=zblob10-z
         if (dist2.lt.dist) then
          dist=dist2
         endif
        endif 

       else if (probtype.eq.42) then  ! vapordist 3D
        dist=sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)-radblob
! zblob2 (3d)  zblob (2d)
        dist1=zblob2+zblob-z
        if (dist1.lt.dist) then
         dist=dist1
        endif
       else if (probtype.eq.46) then
        print *,"vapordist 3D: this problem not ready in 3d"
        stop
       else if (probtype.eq.37) then
        dist=radblob-sqrt( (x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2 )
       else if (probtype.eq.9) then  ! ship wave free surface (vapordist) 3D
        call boatdist(x,y,z,dist)
!
! ysl 05/12/14  ! calculate dist(2): the vapor distance to its nearest boundary
!
       else if ((probtype.eq.17).or.(probtype.eq.18)) then
        if ((xblob.ne.zero).or.(yblob.ne.zero).or. &
            (radblob.ne.half)) then
         print *,"set xblob=yblob=0 radblob=1/2 drop collide"
         stop
        endif
        if (probtype.eq.18) then
         drat=one
        else
         drat=fort_denconst(3)/fort_denconst(1)  ! dentop/denbot
        endif
        veltop=-one/(drat+one)
        velbot=drat/(drat+one)
        ytop=yblob+1.25*radblob
        ybot=yblob-radblob+(ytop-yblob-radblob)*velbot/veltop

        dist = radblob- &
         sqrt((x-xblob)**2+(y-ytop)**2+(z-zblob)**2)  ! water (3)
        dist1 = radblob- &
         sqrt((x-xblob)**2+(y-ybot)**2+(z-zblob)**2) ! oil (1)
! to dist is the nearest distance to a neighbouring fluid
        if (dist1.gt.dist) then
         dist=dist1
        endif

       else if (probtype.eq.13) then
        print *,"this option obsolete"
        stop
       else if (probtype.eq.14) then
        print *,"this option obsolete"
        stop
       else if (probtype.eq.50) then
        xmin=-50.0
        xmax=twall
        ymin=yblob-radblob
        ymax=yblob+radblob
        zmin=zblob-radblob
        zmax=zblob+radblob
        call cubedist(xmin,xmax,ymin,ymax,zmin,zmax, &
                      x,y,z,dist)
        dist=-dist

        xmin=-50.0
        xmax=50.0
        ymin=-50.0
        ymax=50.0
        zmin=-50.0
        zmax=zblob2
        call cubedist(xmin,xmax,ymin,ymax,zmin,zmax, &
                      x,y,z,dist1)
        dist1=-dist1

        if (dist1.gt.dist) then
         dist=dist1
        endif

       else if (probtype.eq.51) then
        print *,"option obsolete"
        stop
       else if (probtype.eq.52) then
        print *,"option obsolete"
        stop
! vapor level set for gear test problem. (denliquid,visunburn)
! one can have a pool of liquid and/or an incoming jet.
       else if (probtype.eq.563) then
        if (levelrz.eq.COORDSYS_CARTESIAN) then
         dist=zblob2-z
         if (radblob.gt.zero) then
          call cylinderdist(y,z,x,yblob,zblob,radblob,xblob-radblob, &
           xblob+radblob,dist2)
          dist2=-dist2
          if (dist2.gt.dist) then
           dist=dist2
          endif
         endif 
        else 
         print *,"levelrz invalid vapordist"
         stop
        endif
       else if ((probtype.eq.56).or.(probtype.eq.562)) then
        dist=zblob2-z
! dog
       else if (probtype.eq.5600) then
        dist=zblob2-z
! vessel (viorel sphere)
       else if (probtype.eq.5601) then
        dist=zblob2-z
       else if (probtype.eq.5602) then ! internal inflow
        dist=zblob2-z
! microfluidics channel -- 0,3,4 Roper, 1 Comsol, 2,5 squeeze vapordist
       else if (probtype.eq.5700) then
        if (axis_dir.eq.2) then
         print *,"axis_dir=2 obsolete"
         stop
        else if (axis_dir.eq.5) then
         dist=y-yblob
        else if ((axis_dir.eq.0).or.(axis_dir.eq.1).or. &
                 (axis_dir.eq.3).or.(axis_dir.eq.4)) then
         dist=y-yblob
         if (xblob4.gt.zero) then
          if (x.ge.xblob4) then  ! oil initially in the curly cue channel
           dist=-99999.0
          endif
         endif
        else
         print *,"axis_dir invalid"
         stop
        endif
       else if (probtype.eq.57) then
        dist=zblob2-z
       else if (probtype.eq.58) then
        dist=zblob2-z
       else if (probtype.eq.529) then ! airblast with coaxial injection
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
         ! 3D jetbend problem
       else if (probtype.eq.53) then ! JICF (vapordist)
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
       else if (probtype.eq.530) then ! impinge (vapordist)
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
        ! 3d diesel injector w/needle (vapordist)
       else if ((probtype.eq.538).or. &  ! inputs.injA
                (probtype.eq.541)) then
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
       else if (probtype.eq.701) then ! flapping wing (vapordist 3D)
        dist=zblob-z
       else if (probtype.eq.532) then ! impinge from sides (vapordist)
        call get_jet_dist(x,y,z,distbatch)
        dist=distbatch(1)
       else if (probtype.eq.536) then
        zmin=zblob-1.0e+10
        zmax=zblob+1.3
        radblob2=0.85
        call cylinderdist(x,y,z,xblob,yblob,radblob2,zmin,zmax,dist)
        dist=-dist
       else if (probtype.eq.41) then ! pipe vapordist 3D
        call inletpipedist(x,y,z,distbatch) 
        dist=distbatch(1)
       else if (probtype.eq.54) then
        zmin=-1.0e+10
        zmax=zblob+zblob/four
        call cylinderdist(x,y,z,xblob,yblob,radblob,zmin,zmax,dist)
        dist=-dist
       else if (probtype.eq.62) then
        dist=zblob3-z
        costheta=cos(xblob2)
        sintheta=sin(xblob2)
        xprime=costheta*(x-xblob)-sintheta*(z-zblob)
        yprime=y-yblob
        zprime=sintheta*(x-xblob)+costheta*(z-zblob)
        delta=half*(radblob2-radblob)
        temprad=1.0D+10
        cylinder_zlo=-half*zblob2-delta
        call cylinderdist(xprime,yprime,zprime, &
          stub_zero,stub_zero,radblob+delta, &
          cylinder_zlo,temprad,dist1) 
        if (dist1.lt.dist) then
         dist=dist1
        endif 
       else if (probtype.eq.63) then
         if (z.ge.two*xblob10) then
          dist=half*xblob10-sqrt(x**2+y**2)
         else
          dist=half*xblob10-sqrt(x**2+y**2)+(two*xblob10-z)
         endif
       else if (probtype.eq.44) then
        call damdist(x,z,dist,igeom)
       else if (probtype.eq.65) then
        dist=zblob3-z
        costheta=cos(xblob2)
        sintheta=sin(xblob2)
        xprime=costheta*(x-xblob)-sintheta*(z-zblob)
        yprime=y-yblob
        zprime=sintheta*(x-xblob)+costheta*(z-zblob)
        delta=half*(radblob2-radblob)
        temprad=1.0D+10
        cylinder_zlo=-half*zblob2-delta
        call tcylinderdist(xprime,yprime,zprime, &
          stub_zero,stub_zero,radblob+delta, &
          cylinder_zlo,temprad,dist1) 
        dist1 = -dist1
        if ((dist1.gt.dist).and.(z > zblob3)) then
         dist=dist1
        endif 
       else if (probtype.eq.66) then
        xprime=sqrt(three*radblob/(four*zblob*zblob*zblob))
        dist=zblob+radblob/(cosh(x*xprime)**2)-z
       else if (probtype.eq.61) then
        dist=-(sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2) - radblob)
        dist1=-(3 + (z-zblob))
        if (axis_dir.eq.2) then
         dist1=-(z-zblob2) 
        endif
        if (dist1.gt.dist) then
         dist=dist1
        endif
       else if (probtype.eq.26) then ! 3D: vapordist
! x-y plane
        if (axis_dir.eq.3) then !3D vortex confinement, with interface
         dist=abs(y-half)-one/four
! x-z plane
        else if (axis_dir.eq.2) then !3D vortex confinement, with interface
         dist=abs(z-half)-one/four
! y-z plane
        else if (axis_dir.eq.1) then !3D vortex confinement, with interface
         dist=abs(z-half)-one/four
        else if (axis_dir.eq.0) then ! 3D swirl, no interface
         ! do nothing
        else if (axis_dir.eq.4) then ! 3D vortex confinement, no interface
         ! do nothing
        else if (axis_dir.eq.5) then ! 3D vortex confinement, w/interface
         dist=sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)-radblob
        else if (axis_dir.eq.11) then ! BCG smooth test, periodic
         ! do nothing
        else if (axis_dir.eq.12) then ! buoyancy test, periodic
         ! do nothing
        else
         print *,"axis_dir 0,1,2,3,4,5,11 or 12 expected"
         stop
        endif
       endif

      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine vapordist


      subroutine airgunsolid(x,y,xcen,ycen,xhole,yhole, &
                             height,width,gunthick,dist)
      use global_utility_module
      use global_distance_module

      IMPLICIT NONE

!        width/2
! ____    ____  
!  ___|  |__  |  
!  |        | | height/2
! _|        |_|          
! _          _  2yhole
!  |        | | 
!  |__    __| |
!  ___|  |____| gunthick  
!     2xhole
! Negative on the inside of the L-shaped solids!
! The subroutine below is intended to take advantage
! of axisymmetry so we only draw the right portion.
!
      real(amrex_real) x,y,xcen,ycen,xhole,yhole
      real(amrex_real) dist
      real(amrex_real) height,width,gunthick
      real(amrex_real) xlo1,xhi1,ylo1,yhi1
      real(amrex_real) xlo2,xhi2,ylo2,yhi2
      real(amrex_real) xlo3,xhi3,ylo3,yhi3
      real(amrex_real) xlo4,xhi4,ylo4,yhi4
      
      if ( (xhole.lt.zero).or.(yhole.lt.zero).or. &
          (height.le.zero).or.(width.le.zero).or. &
          (gunthick.le.zero) ) then 
       print *,"These parameters must be positive ",xhole, &
                yhole,height,width,gunthick
       stop
      endif
! bottom horizontal rectangle
      xlo1=xcen+xhole
      xhi1=xcen+xhole+half*width
      ylo1=ycen-half*height-yhole
      yhi1=ycen-half*height-yhole+gunthick
! top horizontal rectangle
      xlo2=xcen+xhole
      xhi2=xcen+xhole+half*width
      ylo2=ycen+half*height+yhole-gunthick
      yhi2=ycen+half*height+yhole
! bottom-right vertical rectangle
      xlo3=xcen+xhole+half*width-gunthick
      xhi3=xcen+xhole+half*width
      ylo3=ycen-half*height-yhole+gunthick
      yhi3=ycen-yhole
! top-right vertical rectangle
      xlo4=xcen+xhole+half*width-gunthick
      xhi4=xcen+xhole+half*width
      ylo4=ycen+yhole
      yhi4=ycen+half*height+yhole-gunthick
      
      if (y.le.yhi1) then
       call squaredist(x,y,xlo1,xhi1,ylo1,yhi1,dist)       
      endif
      if (y.gt.ylo2) then 
       call squaredist(x,y,xlo2,xhi2,ylo2,yhi2,dist)       
      endif
      if ((y.ge.ylo3).and.(y.le.(yhi3+yhole))) then
       call squaredist(x,y,xlo3,xhi3,ylo3,yhi3,dist)          
      endif
      if ((y.le.yhi4).and.(y.ge.(ylo4-yhole))) then
       call squaredist(x,y,xlo4,xhi4,ylo4,yhi4,dist)       
      endif
      
      return
      end subroutine



! dist>0 inside of gap
      subroutine gapdist(x,y,xcen,gap,ymin,ymax,tana,dist)
      IMPLICIT NONE


      real(amrex_real) x,y,xcen,gap,dist
      real(amrex_real) ymin,ymax,tana ! tangent of direction angle

      if (ymin.ge.ymax-1.0E-10) then
       print *,"invalid parameters ",ymin,ymax
       stop
      endif
      dist=abs(x-xcen-tana*(y-ymin))-gap
      if (y.ge.ymax) then
       if (dist.le.zero) then
        dist=y-ymax
       else
        dist=sqrt(dist**2+(y-ymax)**2)
       endif
      else if (y.le.ymin) then
       if (dist.le.zero) then
        dist=ymin-y
       else
        dist=sqrt(dist**2+(ymin-y)**2)
       endif
      endif

      return
      end subroutine





! dist>0 inside of droplet
      subroutine legenddist(x,y,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,dist,mag,costheta,sintheta
      real(amrex_real) cos2theta,sin2theta

      mag=sqrt((x-xblob)**2+(y-yblob)**2)
      if (mag.lt.EPS5) then
       dist=radblob+zblob
      else
       costheta=(y-yblob)/mag
       sintheta=(x-xblob)/mag
       cos2theta=costheta**2-sintheta**2
       sin2theta=two*costheta*sintheta
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        dist=radblob+zblob*(two*(costheta**2)-one)/two - mag
       else if (levelrz.eq.COORDSYS_RZ) then 
        dist=radblob+zblob*(three*(costheta**2)-one)/two - mag
!       dist=radblob+zblob*(three*(cos2theta**2)-one)/two - mag
       else
        print *,"levelrz invalid legenddist"
        stop
       endif
      endif

      return
      end subroutine

! dist>0 inside of the ellipse
      subroutine ellipsedist(x,y,a,b,xc,yc,dist)
      IMPLICIT NONE

      real(amrex_real) x,y,a,b,xc,yc,dist
      real(amrex_real) xprime,yprime,xcritical,ycritical
      real(amrex_real) factor

      xprime=x-xc
      yprime=y-yc
      if (xprime.lt.zero) then
       xprime=-xprime
      endif
      if (yprime.lt.zero) then
       yprime=-yprime
      endif
      if ((xprime.eq.zero).and.(yprime.eq.zero)) then
       dist=a
       if (dist.gt.b) then
        dist=b
       endif
      else if (xprime.gt.yprime) then
       factor=one/a**2 + (yprime/(xprime*b))**2
       xcritical=one/sqrt(factor)
       ycritical=yprime*xcritical/xprime
      else 
       factor=one/b**2 + (xprime/(yprime*a))**2
       ycritical=one/sqrt(factor)
       xcritical=xprime*ycritical/yprime
      endif
      dist=sqrt(xcritical**2+ycritical**2)- &
           sqrt(xprime**2+yprime**2)

      return
      end subroutine
     
      subroutine get_bump_dist(x,y,z,dist,time)
      use global_utility_module
      use global_distance_module

      IMPLICIT NONE

      real(amrex_real), INTENT(out) :: dist(num_materials)
      real(amrex_real), INTENT(in) :: x,y,z,time
      real(amrex_real) distsolid
      real(amrex_real) distleft,distright
      real(amrex_real) xleft,xright,elev
      integer im_solid_bump

      im_solid_bump=im_solid_primary()

      if ((im_solid_bump.lt.1).or. &
          (im_solid_bump.gt.num_materials)) then
       print *,"im_solid_bump invalid in get_bump_dist"
       stop
      endif
      if (probtype.ne.110) then
       print *,"probtype invalid"
       stop
      endif
      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y in get_bump_dist"
        stop
       endif
      endif

      xleft=-594.36
      xright=502.92
      call get_left_elevationIOWA(time,distleft)
      call get_right_elevationIOWA(time,distright)
      call local_shallow_water_elevation(time,x,elev)
      dist(1)=elev-y
      dist(2)=-dist(1)

        ! pos in solid (get_bump_dist)
      call materialdistsolid(x,y,z,distsolid,time,im_solid_bump)
      dist(im_solid_bump)=distsolid

      return
      end subroutine get_bump_dist


 
      subroutine get_bump_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc,time)
      use geometry_intersect_module
      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, parameter :: nhalf2=1
      real(amrex_real) xsten2(-nhalf2:nhalf2,SDIM)

      real(amrex_real), INTENT(in) :: time
      integer im
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(out) :: vfrac(num_materials)
      real(amrex_real), INTENT(out) :: cenbc(num_materials,SDIM)

      integer dir2,i1,j1,k1,k1lo,k1hi,isten
      real(amrex_real) centroid(num_materials,SDIM)
      real(amrex_real) lsgrid(D_DECL(3,3,3),num_materials)
      real(amrex_real) facearea(num_materials)
      real(amrex_real) distbatch(num_materials)
      real(amrex_real) EBVOFTOL

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid get bump vfrac"
       stop
      endif


      if (probtype.ne.110) then
       print *,"probtype invalid"
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

      do k1=k1lo,k1hi
      do j1=-1,1
      do i1=-1,1
       do isten=-1,1
        dir2=1
        xsten2(isten,dir2)=xsten(isten+2*i1,dir2)
        dir2=2
        xsten2(isten,dir2)=xsten(isten+2*j1,dir2)
        if (SDIM.eq.3) then
         dir2=SDIM
         xsten2(isten,dir2)=xsten(isten+2*k1,dir2)
        endif
       enddo ! isten

       call get_bump_dist(xsten2(0,1),xsten2(0,2), &
        xsten2(0,SDIM),distbatch,time)
       do im=1,num_materials
        lsgrid(D_DECL(i1+2,j1+2,k1+2),im)=distbatch(im) 
       enddo
 
      enddo
      enddo
      enddo

      EBVOFTOL=VOFTOL
      call getvolumebatch(bfact,dx,xsten,nhalf, &
        lsgrid,vfrac, &
        facearea,centroid,EBVOFTOL,SDIM)

      do im=1,num_materials
       do dir2=1,SDIM
        cenbc(im,dir2)=centroid(im,dir2)-xsten(0,dir2)
       enddo
      enddo

      return
      end subroutine get_bump_vfrac


      subroutine get_bump_velocity(xsten,nhalf,dx,bfact,vel,time)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(out) :: vel

      real(amrex_real) VOF(num_materials)
      real(amrex_real) cenbc(num_materials,SDIM)
      real(amrex_real) xleft,xright
      real(amrex_real) velleft,velright

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid get bump velocity"
       stop
      endif

      if (probtype.ne.110) then
       print *,"probtype invalid"
       stop
      endif

      xleft=-594.36
      xright=502.92
      call get_left_velocityIOWA(time,velleft)
      call get_right_velocityIOWA(time,velright)

      call get_bump_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc,time)  

      if (VOF(1).gt.zero) then

       call local_shallow_water_velocity(time,xsten(0,1),vel)

      else
       vel=zero
      endif

      return
      end subroutine get_bump_velocity

      subroutine mask_velocity(xsten,nhalf,dx,bfact,velcell,time)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(inout) :: velcell(SDIM)

      real(amrex_real) VOF(num_materials)
      real(amrex_real) cenbc(num_materials,SDIM)
      integer im

      if (nhalf.lt.3) then
       print *,"nhalf invalid mask velocity"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      call get_initial_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc)

      do im=1,num_materials
       if (is_rigid(im).eq.1) then
        if (VOF(im).gt.zero) then
         call velsolid(xsten(0,1),xsten(0,2),xsten(0,SDIM),velcell,time,im,dx)
        endif
       else if (is_rigid(im).eq.0) then
        ! do nothing
       else
        print *,"is_rigid(im) invalid"
        stop
       endif
      enddo ! im=1..num_materials

      return
      end subroutine mask_velocity

      subroutine mofBC(time,dir,side,VOF,VOFwall, &
       xsten,nhalf,dx,bfact,im)
      IMPLICIT NONE

      integer, INTENT(in) :: dir,side,im,bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(out) :: VOF
      real(amrex_real) :: xwall
      real(amrex_real), INTENT(in) :: VOFwall
      real(amrex_real) :: x,y,z
      real(amrex_real), INTENT(in) :: dx(SDIM)

      if (nhalf.lt.1) then
       print *,"nhalf invalid mofbc"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)
      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else
        xwall=probhix
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else
        xwall=probhiy
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else
        xwall=probhiz
       endif
      else
       print *,"dir invalid mofbc"
       stop
      endif

      if ((im.lt.0).or.(im.ge.num_materials)) then
       print *,"im invalid77"
       stop
      endif
      if ((dir.eq.1).and.(side.eq.1)) then
       if (xwall.lt.x) then
        print *,"xwall,x invalid"
        stop
       endif
       VOF=VOFwall
      else if ((dir.eq.1).and.(side.eq.2)) then
       if (xwall.gt.x) then
        print *,"xwall,x invalid"
        stop
       endif
       VOF=VOFwall
      else if ((dir.eq.2).and.(side.eq.1)) then
       if (xwall.lt.y) then
        print *,"xwall,y invalid"
        stop
       endif
       VOF=VOFwall
      else if ((dir.eq.2).and.(side.eq.2)) then
       if (xwall.gt.y) then
        print *,"xwall,y invalid"
        stop
       endif
       VOF=VOFwall
      else if ((dir.eq.3).and.(side.eq.1)) then
       if (xwall.lt.z) then
        print *,"xwall,z invalid"
        stop
       endif
       VOF=VOFwall
      else if ((dir.eq.3).and.(side.eq.2)) then
       if (xwall.gt.z) then
        print *,"xwall,z invalid"
        stop
       endif
       VOF=VOFwall
      else
       print *,"dir side invalid"
       stop
      endif

      return
      end subroutine mofBC

      subroutine copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
      IMPLICIT NONE

      real(amrex_real), INTENT(out) :: VOF(num_materials*ngeom_raw)
      real(amrex_real), INTENT(in) :: VOFwall(num_materials*ngeom_raw)
      real(amrex_real), INTENT(in) :: cenbc(num_materials,SDIM)
      real(amrex_real), INTENT(in) :: vofarray(num_materials)
      integer im,vofcomp,dir2

      do im=1,num_materials
       vofcomp=(im-1)*(ngeom_raw)+1
       if ((FSI_flag(im).eq.FSI_FLUID).or. & 
           (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90).or. & 
           (FSI_flag(im).eq.FSI_ICE_PROBF90).or. & 
           (FSI_flag(im).eq.FSI_ICE_STATIC).or. & 
           (FSI_flag(im).eq.FSI_RIGID_NOTPRESCRIBED)) then 
        VOF(vofcomp)=vofarray(im)
        do dir2=1,SDIM
         VOF(vofcomp+dir2)=cenbc(im,dir2)
        enddo
       else if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
                (FSI_flag(im).eq.FSI_ICE_NODES_INIT).or. & 
                (FSI_flag(im).eq.FSI_FLUID_NODES_INIT).or. & 
                (FSI_flag(im).eq.FSI_SHOELE_CTML)) then 
        VOF(vofcomp)=VOFwall(vofcomp)
        do dir2=1,SDIM
         VOF(vofcomp+dir2)=VOFwall(vofcomp+dir2)
        enddo
       else
        print *,"FSI_flag invalid in copy_mofbc_to_result"
        print *,"im,FSI_flag(im): ",im,FSI_flag(im)
        stop
       endif
      enddo ! im=1..num_materials

      return
      end subroutine copy_mofbc_to_result

      subroutine check_lsbc_extrap(LS,LSwall)
      IMPLICIT NONE

      real(amrex_real), INTENT(out) :: LS(num_materials)
      real(amrex_real), INTENT(in) :: LSwall(num_materials)
      integer im

      do im=1,num_materials
       if ((FSI_flag(im).eq.FSI_FLUID).or. & 
           (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90).or. & 
           (FSI_flag(im).eq.FSI_ICE_PROBF90).or. & 
           (FSI_flag(im).eq.FSI_ICE_STATIC).or. & 
           (FSI_flag(im).eq.FSI_RIGID_NOTPRESCRIBED)) then 
        ! do nothing
       else if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
                (FSI_flag(im).eq.FSI_ICE_NODES_INIT).or. & 
                (FSI_flag(im).eq.FSI_FLUID_NODES_INIT).or. & 
                (FSI_flag(im).eq.FSI_SHOELE_CTML)) then 
        LS(im)=LSwall(im)
       else
        print *,"FSI_flag invalid in check_lsbc_extrap"
        print *,"im,FSI_flag(im): ",im,FSI_flag(im)
        stop
       endif
      enddo ! im=1..num_materials

      return
      end subroutine check_lsbc_extrap


! inflow, outflow, slipwall, noslipwall bcs have ext dir; i.e.
! this routine called for inflow, outflow, slipwall, noslipwall

      subroutine groupmofBC(time,dir,side,VOF,VOFwall, &
        xsten,nhalf,dx,bfact)
      use global_utility_module
      use randomNG
      use rainControl_module, only : get_rain_vfrac,get_rain_velocity
      use bubbleControl_module, only : get_pack_vfrac
      use hydrateReactor_module
      use unimaterialChannel_module
      use River
      IMPLICIT NONE

      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: dir,side
      real(amrex_real), INTENT(out) :: VOF(num_materials*ngeom_raw)
      real(amrex_real), INTENT(in) :: VOFwall(num_materials*ngeom_raw)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real) :: xwall,x,y,z
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) cenbc(num_materials,SDIM)
      real(amrex_real) LS
      real(amrex_real) vofarray(num_materials)
      integer im
      integer vofcomp
      integer dir2

      real(amrex_real) vof_sum_check
      integer im_solid_mofbc

      im_solid_mofbc=im_solid_primary()

      if (nhalf.lt.1) then
       print *,"nhalf invalid group mof bc"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else
        xwall=probhix
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else
        xwall=probhiy
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else
        xwall=probhiz
       endif
      else
       print *,"dir invalid group mof bc 0"
       stop
      endif

      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid in group mof bc 1"
       stop
      endif

      do im=1,num_materials*ngeom_raw
       VOF(im)=zero
      enddo

      vof_sum_check=zero

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_raw+1
       VOF(vofcomp)=VOFwall(vofcomp)
       do dir2=1,SDIM
        VOF(vofcomp+dir2)=VOFwall(vofcomp+dir2)
        cenbc(im,dir2)=zero
       enddo
       vofarray(im)=VOFwall(vofcomp)

       vof_sum_check=vof_sum_check+vofarray(im)
      enddo ! im=1..num_materials

      if (vof_sum_check.ge.half) then
       ! do nothing
      else
       print *,"(breakpoint) break point and gdb: "
       print *,"(1) compile with the -g option"
       print *,"(2) break PROB.F90:10814"
       print *,"vof_sum_check failed"
       print *,"vof_sum_check=",vof_sum_check
       print *,"time=",time
       print *,"dir,side ",dir,side
       print *,"xsten(0) ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
       print *,"nhalf ",nhalf
       print *,"dx= ",dx(1),dx(2),dx(SDIM)
       print *,"bfact=",bfact
       print *,"num_materials=",num_materials
       print *,"VOFwall= ",VOFwall
       print *,"vofarray= ",vofarray
       stop
      endif

      if (is_in_probtype_list().eq.1) then
       call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
       call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)

      else if (probtype.eq.401) then
       call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
       call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)

      else if (probtype.eq.412) then
       call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
       call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)

      else if (probtype.eq.311) then ! user defined

       call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
       call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)

      else if (probtype.eq.199) then !hydrates (groupmofBC)
       print *,"this code must be upgrated"
       stop 

       do im=1,num_materials
        if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. &
            (FSI_flag(im).eq.FSI_SHOELE_CTML)) then
         !do nothing
        else
         vofcomp=(im-1)*ngeom_raw+1
         call HYD_VOLF_BC(time,dir,side,VOF(vofcomp), &
          xwall,VOFwall(vofcomp),x,y,z,dx,im)
         call HYD_MOMF_BC(time,dir,side,VOF(vofcomp+1), &
          xwall,VOFwall(vofcomp+1),x,y,z,dx,im)
        endif
       enddo ! im=1..num_materials

       ! in: groupmofBC
      else if (probtype.eq.220) then !UNIMATERIAL problem

       do im=1,num_materials-1
        vofcomp=(im-1)*ngeom_raw+1
        call UNIMAT_VOLF_BC(time,dir,side,VOF(vofcomp), &
         xwall,VOFwall(vofcomp),x,y,z,dx,im)
        call UNIMAT_MOMF_BC(time,dir,side,VOF(vofcomp+1), &
         xwall,VOFwall(vofcomp+1),x,y,z,dx,im)
       enddo

       im=num_materials
       vofcomp=(im-1)*ngeom_raw+1
       if (FSI_flag(im).eq.FSI_SHOELE_CTML) then
        VOF(vofcomp)=zero
        do dir2=1,SDIM
         VOF(vofcomp+dir2)=zero
        enddo
       else
        print *,"expecting last material to be a CTML material"
        stop
       endif
        
      else if ((probtype.eq.299).or. &
               (probtype.eq.301)) then !melting (boundary conditions for F,X)
       call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
       call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
      else if (probtype.eq.209) then ! River
       call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
       call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        ! marangoni (heat pipe) problem
      else if ((probtype.eq.36).and.(axis_dir.eq.10)) then
       call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
       call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)

       ! sanity check: curvature for a plane should be 0
      else if ((probtype.eq.36).and.(axis_dir.eq.210)) then
       call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
       call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
      else
 
       if ((dir.eq.1).and.(side.eq.1)) then  ! xlo
        if (xwall.lt.x) then
         print *,"xwall,x invalid"
         stop
        endif

        if (SDIM.eq.2) then

        if (probtype.eq.58) then
         LS=yblob2-y
         if (LS.lt.zero) then
          vofarray(1)=zero
         else
          vofarray(1)=one
         endif
         vofarray(2)=one-vofarray(1)
         do im=3,num_materials
          vofarray(im)=zero
         enddo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.6) then
         vofarray(1)=zero
         vofarray(2)=one-vofarray(1)
         do im=3,num_materials
          vofarray(im)=zero
         enddo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.110) then
           ! xlo
         call get_bump_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc,time) 
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.59) then  ! xlo dir=1 side=1 groupmofBC 2d
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
         ! xlo groupmofBC 2d
        else if ((probtype.eq.801).and.(axis_dir.eq.dir-1)) then 
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.802) then ! xlo - dissolution
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5700) then ! xlo
         call get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then  ! xlo, groupmofBC, 2D
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.41) then
         call get_pipe_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)  ! xlo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.532) then ! xlo
         call get_jet_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)   
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.539) then ! xlo, sup injector - groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.202) then ! xlo, liquidlens - groupmofbc
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.701) then  ! xlo dir=1 side=1 groupmofBC rain
         if ((axis_dir.eq.0).or. &
             (axis_dir.eq.1)) then
          call get_rain_vfrac(x,y,z,dx,vofarray,cenbc,time,dir)
          call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
         else if (axis_dir.eq.2) then

          do im=1,num_materials
           if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. &
               (FSI_flag(im).eq.FSI_SHOELE_CTML)) then
            ! do nothing
           else
            vofcomp=(im-1)*(ngeom_raw)+1
            if (im.eq.1) then
             VOF(vofcomp)=one
            else 
             VOF(vofcomp)=zero
            endif
            do dir2=1,SDIM
             VOF(vofcomp+dir2)=zero
            enddo
           endif
          enddo ! im=1..num_materials
         else
          print *,"axis_dir invalid"
          stop
         endif
        endif

        else if (SDIM.eq.3) then

        if (probtype.eq.50) then
         vofarray(1)=one
         vofarray(2)=one-vofarray(1)
         do im=3,num_materials
          vofarray(im)=zero
         enddo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.58) then
         LS=zblob2-z
         if (LS.lt.zero) then
          vofarray(1)=zero
         else
          vofarray(1)=one
         endif
         vofarray(2)=one-vofarray(1)
         do im=3,num_materials
          vofarray(im)=zero
         enddo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.41) then ! xlo 3D
         call get_pipe_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)  ! xlo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5700) then ! xlo
         call get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc) 
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
          ! airblast
        else if (probtype.eq.529) then  ! dir=1 side=1 groupmofBC
         call get_jet_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.59) then  ! xlo dir=1 side=1 groupmofBC 3d
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5501) then  ! xlo, groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then  ! xlo, groupmofBC, 3D
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

        else
         print *,"dimension bust"
         stop
        endif

       else if ((dir.eq.1).and.(side.eq.2)) then  ! xhi
        if (xwall.gt.x) then
         print *,"xwall,x invalid"
         stop
        endif

        if (SDIM.eq.2) then

        if (probtype.eq.110) then
            ! xhi
         call get_bump_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc,time)  
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.59) then  ! dir=1 side=2 (xhi) 2d
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)

         ! outflow or inflow=EXT_DIR; for microfluid channels, 
         ! want extrap bc either case.
         ! so this code is disabled
        else if ((probtype.eq.5700).and.(1.eq.0)) then ! xhi
         call get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then ! xhi, groupmofbc, 2D
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.41) then
         call get_pipe_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)  ! xhi
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.532) then ! xhi
         call get_jet_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc) 
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.539) then ! xhi, sup injector
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.202) then ! xhi, liquidlens, groupmofbc
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if ((probtype.eq.25).and.(axis_dir.gt.0)) then ! xhi, bubble frm.
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

        else if (SDIM.eq.3) then

         ! outflow or inflow=EXT_DIR; for microfluid channels, 
         ! want extrap bc either case.
         ! so this code is disabled
        if ((probtype.eq.5700).and.(1.eq.0)) then ! xhi
         call get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc) 
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.59) then  ! dir=1 side=2 (xhi) 3d groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5501) then  ! xhi groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then  ! xhi, groupmofBC, 3D
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

        else 
         print *,"dimension bust"
         stop
        endif

       else if ((dir.eq.2).and.(side.eq.1)) then ! ylo
        if (xwall.lt.y) then
         print *,"xwall,y invalid"
         stop
        endif

        if (SDIM.eq.2) then

        if ((probtype.eq.25).and.(axis_dir.gt.0)) then  ! ylo2D bubble formation

         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)

         ! ylo dir=2 side=1 groupmofBC 
        else if(probtype.eq.bubbleInPackedColumn)then  
         call get_pack_vfrac(x,y,z,dx,vofarray,cenbc,time,dir)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
         ! ylo 2D groupmofBC
        else if ((probtype.eq.801).and.(axis_dir.eq.dir-1)) then 
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5700) then  ! ylo 2D
         call get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.41) then
         call get_pipe_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)  ! ylo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.72) then
         if (abs(x-xblob).gt.radblob) then
          vofarray(1)=zero
         else
          vofarray(1)=one
         endif
         vofarray(2)=one-vofarray(1)
         do im=3,num_materials
          vofarray(im)=zero
         enddo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if ((probtype.eq.53).or.(probtype.eq.537)) then !ylo 2D,groupmofBC
         call get_jet_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc) 
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if ((probtype.eq.538).or. &
                 (probtype.eq.539).or. &
                 (probtype.eq.541)) then ! ylo 2D, diesel injector, groupmofBC

         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)

         ! 540 Rieber problem
         ! ylo dir=2 side=1 2d groupmofBC
        else if ((probtype.eq.59).or. &
                 (probtype.eq.540)) then  
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then ! ylo, groupmofbc, 2D
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

        else if (SDIM.eq.3) then

        if (probtype.eq.5700) then  ! ylo
         call get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc) 
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.59) then  ! ylo dir=2 side=1 3d groupmofbc
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5501) then  ! ylo groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then  ! ylo, groupmofBC, 3D
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.532) then  ! impinge from side, ylo
         call get_jet_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

        else
         print *,"dimension bust"
         stop
        endif

       else if ((dir.eq.2).and.(side.eq.2)) then ! yhi
        if (xwall.gt.y) then
         print *,"xwall,y invalid"
         stop
        endif

        if (SDIM.eq.2) then

        if (probtype.eq.62) then
         vofarray(1)=zero
         vofarray(2)=one
         do im=3,num_materials
          vofarray(im)=zero
         enddo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5700) then ! yhi
         call get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if ((probtype.eq.14).or.(probtype.eq.16).or. &
          ((probtype.eq.25).and.(axis_dir.eq.0)) ) then
         if (probtype.eq.25) then
          vofarray(1)=one
          vofarray(2)=zero
         else if (probtype.eq.16) then
          vofarray(1)=zero
          vofarray(2)=one
         endif
         do im=3,num_materials
          vofarray(im)=zero
         enddo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
         ! top (y=yhi)  y=ymax wall
        else if (probtype.eq.41) then
         call get_pipe_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc) 
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.540) then  ! Rieber problem y=yhi
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.539) then ! yhi, sup injector
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then ! yhi, groupmofBC, 2D
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if ((probtype.eq.701).and.(1.eq.0)) then ! yhi, groupmofBC 2D
         call get_rain_vfrac(x,y,z,dx,vofarray,cenbc,time,dir)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

        else if (SDIM.eq.3) then

        if (probtype.eq.5700) then ! yhi
         call get_microfluidic_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc) 
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5501) then  ! yhi groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.59) then  ! yhi dir=2 side=2 3d groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then ! yhi, groupmofBC, 3D
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.532) then  ! impinge from side, yhi
         call get_jet_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

        else 
         print *,"dimension bust"
         stop
        endif

       else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then  ! zlo
        if (xwall.lt.z) then
         print *,"xwall,z invalid"
         stop
        endif

         ! called at inflow and outflow bc

        if ((probtype.eq.538).or. &  ! inputs.injA
            (probtype.eq.541)) then  ! zlo, diesel injector

         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)   
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
         
        else if ((probtype.eq.53).or.(probtype.eq.531).or. &  ! zlo 3D
                 (probtype.eq.530).or.(probtype.eq.536).or. &
                 (probtype.eq.537).or.(probtype.eq.532)) then
         call get_jet_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)   
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
         ! 540 Rieber problem
         ! zlo dir=3 side=1 groupmofBC (3D)
        else if ((probtype.eq.59).or. &
                 (probtype.eq.540)) then  
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5501) then  ! zlo groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then  ! zlo, groupmofBC
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5700) then  ! microfluidics zlo
         do im=1,num_materials
          vofarray(im)=zero
         enddo
         if ((im_solid_mofbc.le.0).or. &
             (im_solid_mofbc.gt.num_materials)) then
          print *,"im_solid_mofbc invalid 12"
          stop
         endif
         vofarray(im_solid_mofbc)=one
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

       else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then ! zhi
        if (xwall.gt.z) then
         print *,"xwall,z invalid"
         stop
        endif
        if ((probtype.eq.62).or.(probtype.eq.65)) then
         vofarray(1)=zero
         vofarray(2)=one-vofarray(1)
         do im=3,num_materials
          vofarray(im)=zero
         enddo
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.540) then  ! Rieber problem z=zhi
         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.9) then  ! zhi, groupmofBC

         call get_initial_vfrac(xsten,nhalf,dx,bfact,vofarray,cenbc)
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        else if (probtype.eq.5700) then  ! microfluidics zhi
         do im=1,num_materials
          vofarray(im)=zero
         enddo
         if ((im_solid_mofbc.lt.1).or. &
             (im_solid_mofbc.gt.num_materials)) then
          print *,"im_solid_mofbc invalid 13"
          stop
         endif
         vofarray(im_solid_mofbc)=one
         call copy_mofbc_to_result(VOF,vofarray,cenbc,VOFwall)
        endif

       else
        print *,"dir side invalid"
        stop
       endif

      endif ! non-hydrate materials

      return
      end subroutine groupmofBC


      subroutine lsBC(time,dir,side,VOF,VOFwall, &
        xsten,nhalf,dx,bfact,im)
      use hydrateReactor_module
      use unimaterialChannel_module
      use River

      IMPLICIT NONE

      integer, INTENT(in) :: dir,side,im,nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: VOF
      real(amrex_real) :: xwall
      real(amrex_real), INTENT(inout) :: VOFwall
      real(amrex_real) :: x,y,z
      real(amrex_real), INTENT(in) :: dx(SDIM)

      if (nhalf.lt.1) then
       print *,"nhalf invalid ls bc"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((im.lt.0).or.(im.ge.num_materials)) then
       print *,"im invalid78"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else
        xwall=probhix
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else
        xwall=probhiy
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else
        xwall=probhiz
       endif
      else
       print *,"dir invalid lsbc"
       stop
      endif

      print *,"lsBC is a placeholder, call GROUP_LSFILL instead"
      print *,"GROUP_LSFILL calls grouplsBC"
      stop

      return
      end subroutine lsBC

! inflow, outflow, slipwall, noslipwall bcs have ext dir; i.e.
! this routine called for inflow, outflow, slipwall, noslipwall
! input: LSwall (LS just inside the domain)
! output: LS
      subroutine grouplsBC(time,dir,side,LS,LSwall, &
        xsten,nhalf,dx,bfact)
      use global_utility_module
      use rainControl_module
      use bubbleControl_module
      use hydrateReactor_module
      use unimaterialChannel_module
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use rigid_FSI_module
      use sinking_particle_module

      IMPLICIT NONE

      integer, INTENT(in) :: nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      integer, INTENT(in) :: dir,side
      real(amrex_real), INTENT(inout) :: LS(num_materials)
      real(amrex_real), INTENT(in) :: LSwall(num_materials)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real) xwall,x,y,z
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) local_LS
      integer im,im_opp,imls
      real(amrex_real) bigdist
      integer im_solid_lsbc
      real(amrex_real) xvec(SDIM)
      integer local_dir
      integer iprob

      im_solid_lsbc=im_solid_primary()

       ! coordinate of cell center that is outside the domain.
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)
      do local_dir=1,SDIM
       xvec(local_dir)=xsten(0,local_dir)
      enddo

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y in 11443"
        stop
       endif
      endif

      bigdist=three*dx(1)

      if (nhalf.lt.1) then
       print *,"nhalf invalid group ls bc"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else
        xwall=probhix
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else
        xwall=probhiy
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else
        xwall=probhiz
       endif
      else
       print *,"dir invalid grouplsbc"
       stop
      endif

      do imls=1,num_materials
       LS(imls)=LSwall(imls)
      enddo

      if (is_in_probtype_list().eq.1) then
       if (1.eq.0) then
        print *,"probtype= ",probtype
        print *,"probtype_list_size= ",probtype_list_size
        do iprob=1,probtype_list_size
         print *,"iprob,used_probtypes(iprob) ",iprob, &
                 used_probtypes(iprob)
        enddo
        stop
       endif
       call SUB_LS_BC(xwall,xvec,time,LS,LSwall,dir,side,dx, &
        num_materials)
       call check_lsbc_extrap(LS,LSWALL)
      else if (probtype.eq.401) then
       call HELIX_LS_BC(xwall,xvec,time,LS,LSwall,dir,side,dx)
       call check_lsbc_extrap(LS,LSWALL)
      else if (probtype.eq.402) then
       call TSPRAY_LS_BC(xwall,xvec,time,LS,LSwall,dir,side,dx)
       call check_lsbc_extrap(LS,LSWALL)
      else if (probtype.eq.412) then ! step
       call CAV2Dstep_LS_BC(xwall,xvec,time,LS,LSwall,dir,side,dx)
       call check_lsbc_extrap(LS,LSWALL)
      else if (probtype.eq.413) then ! zeyu
       call ZEYU_droplet_impact_LS_BC(xwall,xvec,time,LS,LSwall,dir,side,dx)
       call check_lsbc_extrap(LS,LSWALL)
      else if (probtype.eq.533) then
       call rigid_FSI_LS_BC(xwall,xvec,time,LS,LSwall,dir,side,dx)
       call check_lsbc_extrap(LS,LSWALL)
      else if (probtype.eq.534) then
       call sinking_FSI_LS_BC(xwall,xvec,time,LS,LSwall,dir,side,dx)
       call check_lsbc_extrap(LS,LSWALL)

      else if (probtype.eq.311) then ! user defined problem

       call USERDEF_LS_BC(xwall,xvec,time,LS,LSwall,dir,side,dx)
       call check_lsbc_extrap(LS,LSWALL)

       ! curvature sanity check (line in 2D, plane in 3D)
      else if ((probtype.eq.36).and.(axis_dir.eq.210)) then

       call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
       call check_lsbc_extrap(LS,LSWALL)

      else if (probtype.eq.199) then  ! in grouplsBC

       do imls=1,num_materials
        im=imls
        call HYD_LVLS_BC(time,dir,side,LS(imls), &
         xwall,LSwall(imls),x,y,z,dx,im)
       enddo
       call check_lsbc_extrap(LS,LSWALL)

       ! in: grouplsBC
      else if (probtype.eq.220) then  ! UNIMATERIAL problem

       do imls=1,num_materials-1
        im=imls
        call UNIMAT_LVLS_BC(time,dir,side,LS(imls), &
         xwall,LSwall(imls),x,y,z,dx,im)
       enddo
       imls=num_materials
       if (FSI_flag(imls).eq.FSI_SHOELE_CTML) then
        LS(imls)=-ten
       else
        print *,"expecting FSI_flag equal to FSI_SHOELE_CTML"
        stop
       endif

      else if ((probtype.eq.299).or. &
               (probtype.eq.301)) then !melting (boundary condition LS)

       call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
       call check_lsbc_extrap(LS,LSWALL)

      else if (probtype.eq.209) then  ! River
       call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
       call check_lsbc_extrap(LS,LSWALL)

       ! marangoni (heat pipe) problem
      else if ((probtype.eq.36).and.(axis_dir.eq.10)) then
       call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
       call check_lsbc_extrap(LS,LSWALL)
      else

       if ((dir.eq.1).and.(side.eq.1)) then  ! xlo
         if (xwall.lt.x) then
          print *,"xwall,x invalid"
          stop
         endif

         if (SDIM.eq.2) then

         if (probtype.eq.58) then
          local_LS=yblob2-y
          LS(1)=local_LS
          LS(2)=-local_LS
          do im=3,num_materials
           LS(im)=-bigdist
          enddo
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.6) then
          LS(1)=-bigdist
          LS(2)=bigdist
          do im=3,num_materials
           LS(im)=-bigdist
          enddo
          call check_lsbc_extrap(LS,LSWALL)

          ! ysl add level set function for BC 
         else if (probtype.eq.701) then ! xlo grouplsBC
          if ((axis_dir.eq.0).or.(axis_dir.eq.1)) then
           call xloLS_rain(x,y,z,LS,time,bigdist)
          else if (axis_dir.eq.2) then
           LS(1)=bigdist
           LS(2)=-bigdist
          else
           print *,"axis_dir invalid"
           stop
          endif
          call check_lsbc_extrap(LS,LSWALL)
         else if ((probtype.eq.bubbleInPackedColumn).and.(1.eq.2)) then ! xlo 
          call xloLS_pack(x,y,z,LS,time,bigdist)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.110) then
          call get_bump_dist(x,y,z,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
          ! xlo 2d - grouplsBC
         else if ((probtype.eq.801).and.(axis_dir.eq.dir-1)) then 
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.802) then ! xlo: dissolution grouplsBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.59) then  !  xlo: dir=1 side=1 grouplsBC 2d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5700) then ! xlo
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then  ! xlo, groupLSBC, 2d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.41) then ! xlo, groupLSBC, 2d
          call inletpipedist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.532) then ! xlo
          call get_jet_dist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.539) then ! xlo, sup injector - grouplsBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.202) then ! xlo, liquidlens - grouplsbc
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         endif

         else if (SDIM.eq.3) then

         if (probtype.eq.50) then
          LS(1)=bigdist
          LS(2)=-bigdist
          do im=3,num_materials
           LS(im)=-bigdist
          enddo
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.58) then
          local_LS=zblob2-z
          LS(1)=local_LS
          LS(2)=-local_LS
          do im=3,num_materials
           LS(im)=-bigdist
          enddo
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.41) then ! xlo 3D, grouplsbc
          call inletpipedist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5700) then ! xlo
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)

           ! airblast
         else if (probtype.eq.529) then  ! dir=1 side=1 grouplsBC
          call get_jet_dist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.59) then  ! xlo dir=1 side=1 grouplsBC 3d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5501) then  ! xlo, grouplsBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then  ! xlo, grouplsBC, 3D
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         endif

         else
          print *,"dimension bust"
          stop
         endif

       else if ((dir.eq.1).and.(side.eq.2)) then ! xhi
         if (xwall.gt.x) then
          print *,"xwall,x invalid"
          stop
         endif
          ! outflow or inflow=EXT_DIR; for microfluid channels, 
          ! want extrap bc either case.
          ! so this code is disabled

         if (SDIM.eq.2) then

         if (probtype.eq.110) then
          call get_bump_dist(x,y,z,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.59) then  ! dir=1 side=2 (xhi) grouplsBC 2d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if ((probtype.eq.5700).and.(1.eq.0)) then ! xhi
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then  ! xhi, groupLSBC, 2D
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.532) then ! xhi
          call get_jet_dist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.539) then ! xhi, sup injector
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.41) then  ! xhi
          call inletpipedist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.202) then ! xhi, liquidlens, grouplsbc
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if ((probtype.eq.25).and.(axis_dir.gt.0)) then ! xhi, bubble frm.
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         endif

         else if (SDIM.eq.3) then

         if ((probtype.eq.5700).and.(1.eq.0)) then ! xhi
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.59) then  ! xhi dir=1 side=2 3d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5501) then  ! xhi grouplsBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then  ! xhi grouplsBC, 3D
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         endif

         else
          print *,"dimension bust"
          stop
         endif

       else if ((dir.eq.2).and.(side.eq.1)) then  ! ylo
         if (xwall.lt.y) then
          print *,"xwall,y invalid"
          stop
         endif

         if (SDIM.eq.2) then

          !ylo 2D bubble formation
         if ((probtype.eq.25).and.(axis_dir.gt.0)) then 
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
          ! ylo 2D grouplsBC
         else if ((probtype.eq.801).and.(axis_dir.eq.dir-1)) then
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5700) then  ! ylo
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.41) then ! ylo
          call inletpipedist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.72) then
          local_LS=abs(x-xblob)-radblob
          LS(1)=-local_LS
          LS(2)=local_LS
          do im=3,num_materials
           LS(im)=-bigdist
          enddo
          call check_lsbc_extrap(LS,LSWALL)
         else if ((probtype.eq.53).or.(probtype.eq.537)) then ! ylo 2D
          call get_jet_dist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
           ! 2D ylo, diesel injector, grouplsBC
         else if ((probtype.eq.538).or. &
                  (probtype.eq.539).or. &
                  (probtype.eq.541)) then
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
          ! 540 Rieber problem
         else if (probtype.eq.540) then !dir=2 side=1,2d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.59) then !ylo dir=2 side=1,2d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then ! ylo, groupLSBC, 2D
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.bubbleInPackedColumn) then ! ylo grouplsBC
          call yloLS_pack(x,y,z,LS,time,bigdist)
          call check_lsbc_extrap(LS,LSWALL)
         endif

         else if (SDIM.eq.3) then

         if (probtype.eq.5700) then  ! ylo
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.59) then  ! ylo dir=2 side=1 3d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5501) then  ! ylo grouplsBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then  ! ylo, groupLSBC, 3D
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.532) then  ! impinge from side, ylo
          call get_jet_dist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         endif

         else
          print *,"dimension bust"
          stop
         endif

       else if ((dir.eq.2).and.(side.eq.2)) then  ! yhi
         if (xwall.gt.y) then
          print *,"xwall,y invalid"
          stop
         endif

         if (SDIM.eq.2) then

         if (probtype.eq.62) then
          LS(1)=-bigdist
          LS(2)=bigdist
          do im=3,num_materials
           LS(im)=-bigdist
          enddo
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5700) then ! yhi 2D
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if ((probtype.eq.14).or.(probtype.eq.16).or. &
               ((probtype.eq.25).and.(axis_dir.eq.0)) ) then
          if (probtype.eq.25) then
           LS(1)=bigdist
           LS(2)=-bigdist
          else if (probtype.eq.16) then
           LS(1)=-bigdist
           LS(2)=bigdist
          endif
          do im=3,num_materials
           LS(im)=-bigdist
          enddo
          call check_lsbc_extrap(LS,LSWALL)
          ! top y=ymax wall
         else if (probtype.eq.41) then
          call inletpipedist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.540) then  ! Rieber problem y=yhi
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.539) then ! yhi, sup injector
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then ! yhi, groupLSBC, 2D
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.701) then ! yhi, groupLSBC, 2D
          if ((axis_dir.eq.0).or.(axis_dir.eq.1)) then
           call yhiLS_rain(x,y,z,LS,time,bigdist)
           call check_lsbc_extrap(LS,LSWALL)
          else if (axis_dir.eq.2) then
           ! do nothing
          else
           print *,"axis_dir invalid"
           stop
          endif
         endif

         else if (SDIM.eq.3) then

         if (probtype.eq.5700) then ! yhi
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5501) then  ! yhi grouplsBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.59) then  ! yhi dir=2 side=2 3d
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then ! yhi, groupLSBC, 3D
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.532) then  ! impinge from side, yhi
          call get_jet_dist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)
         endif

         else
          print *,"dimension bust"
          stop
         endif

       else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then  ! zlo
         if (xwall.lt.z) then
          print *,"xwall,z invalid"
          stop
         endif

          ! called at inflow and outflow bc

         if ((probtype.eq.538).or. &
             (probtype.eq.541)) then  ! zlo, diesel injector

          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         
         else if ((probtype.eq.53).or.(probtype.eq.531).or. &  ! zlo 3D
                  (probtype.eq.530).or.(probtype.eq.536).or. &
                  (probtype.eq.537).or.(probtype.eq.532)) then
          call get_jet_dist(x,y,z,LS)
          call check_lsbc_extrap(LS,LSWALL)

          ! 540 Rieber problem
         else if (probtype.eq.540) then  ! dir=3 side=1
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.59) then  ! zlo dir=3 side=1 groupLSBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5501) then  ! zlo grouplsBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then  ! zlo grouplsBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5700) then  ! microfluidics zlo
          do im=1,num_materials
           LS(im)=-bigdist
          enddo
          if ((im_solid_lsbc.le.0).or. &
              (im_solid_lsbc.gt.num_materials)) then
           print *,"im_solid_lsbc invalid 15"
           stop
          endif
          LS(im_solid_lsbc)=bigdist
          call check_lsbc_extrap(LS,LSWALL)
         endif

       else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then ! zhi
         if (xwall.gt.z) then
          print *,"xwall,z invalid"
          stop
         endif
         if ((probtype.eq.62).or.(probtype.eq.65)) then
          LS(1)=-bigdist
          LS(2)=bigdist
          do im=3,num_materials
           LS(im)=-bigdist
          enddo
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.540) then  ! Rieber problem z=zhi
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.9) then  ! zhi, groupLSBC
          call materialdist_batch(xsten,nhalf,dx,bfact,LS,time)
          call check_lsbc_extrap(LS,LSWALL)
         else if (probtype.eq.5700) then  ! microfluidics zhi
          do im=1,num_materials
           LS(im)=-bigdist
          enddo
          if ((im_solid_lsbc.lt.1).or. &
              (im_solid_lsbc.gt.num_materials)) then
           print *,"im_solid_lsbc invalid 16"
           stop
          endif
          LS(im_solid_lsbc)=bigdist
          call check_lsbc_extrap(LS,LSWALL)
         endif

       else
         print *,"dir side invalid"
         stop
       endif

      endif ! not hydrates

      if (ls_homflag.eq.0) then
       ! do nothing
      else if (ls_homflag.eq.1) then
       do im=1,num_materials
        if (is_rigid(im).eq.1) then
         if (LS(im).ge.zero) then ! point is in the solid
          do im_opp=1,num_materials
           if (im.ne.im_opp) then
            if (is_rigid(im_opp).eq.0) then
             LS(im_opp)=LSwall(im_opp)  ! default 90 degree contact angle.
            else if (is_rigid(im_opp).eq.1) then
             ! do nothing
            else
             print *,"is_rigid invalid PROB.F90"
             stop
            endif
           else if (im.eq.im_opp) then
            ! do nothing
           else
            print *,"im_opp invalid"
            stop
           endif
          enddo ! im_opp=1..num_materials
         else if (LS(im).lt.zero) then
          ! do nothing
         else
          print *,"in: subroutine grouplsBC"
          print *,"time,dir,side,bfact,num_materials ", &
            time,dir,side,bfact,num_materials
          print *,"probtype,axis_dir ",probtype,axis_dir
          print *,"ls_homflag=",ls_homflag
          print *,"dimension= ",SDIM
          print *,"LS(im) invalid: im, LS(im)= ",im,LS(im)
          do imls=1,num_materials
           print *,"imls,LSwall(imls) ",imls,LSwall(imls)
          enddo
          print *,"x,y,z ",x,y,z
          print *,"bigdist=",bigdist
          stop
         endif
        else if (is_rigid(im).eq.0) then
         ! do nothing
        else
         print *,"is_rigid invalid PROB.F90"
         stop
        endif
       enddo ! im=1..num_materials

      else
       print *,"ls_homflag invalid"
       stop
      endif

      return
      end subroutine grouplsBC

       ! called from:
       !   fort_buildfacewt (LEVELSET_3D.F90)
       !   fort_fluidsolidcor (NAVIERSTOKES_3D.F90)
      subroutine eval_face_coeff( &
       xsten,nhalf, &
       level,finest_level, &
       cc, &
       cc_ice, & !intent(in)
       cc_ice_mask, & !intent(in)
       cc_group, &  ! intent(out)
       dd, &
       dd_group, & ! intent(out)
       visc_coef, &
       nsolve, &
       dir, &
       veldir, & ! veldir=1..nsolve
       project_option, &
       uncoupled_viscosity, &
       side, &
       local_presbc, &
       local_wt)  ! intent(out)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: cc,cc_ice,cc_ice_mask
      real(amrex_real), INTENT(out) :: cc_group
      real(amrex_real), INTENT(in) :: dd
      real(amrex_real) :: ddfactor
      real(amrex_real), INTENT(out) :: dd_group
      real(amrex_real), INTENT(in) :: visc_coef
      integer, INTENT(in) :: nsolve,dir,veldir,project_option
      integer, INTENT(in) :: uncoupled_viscosity
      integer, INTENT(in) :: side
      integer, INTENT(in) :: local_presbc
      real(amrex_real), INTENT(out) :: local_wt(nsolve)
      integer :: at_RZ_boundary

      if ((nhalf.ge.3).and. &
          (cc.ge.zero).and. &
          (cc.le.one).and. &
          (cc_ice.ge.zero).and. &
          (cc_ice.le.one).and. &
          ((cc_ice_mask.eq.zero).or.(cc_ice_mask.eq.one)).and. &
          (dd.ge.zero).and. &
          (visc_coef.ge.zero).and. &
          (nsolve.ge.1).and. &
          (veldir.ge.1).and. &
          (veldir.le.nsolve).and. &
          (dir.ge.0).and. &
          (dir.lt.SDIM).and. &
          (project_option_is_validF(project_option).eq.1).and. &
          ((uncoupled_viscosity.eq.0).or. &
           (uncoupled_viscosity.eq.1)).and. &
          ((side.eq.0).or.(side.eq.1).or.(side.eq.2)).and. &
          (level.ge.0).and. &
          (level.le.finest_level)) then

       at_RZ_boundary=0

       if (SDIM.eq.2) then
        if (dir.eq.0) then
         if (side.eq.1) then !inorm=fablo(dir+1)

          if (levelrz.eq.COORDSYS_RZ) then
           if (local_presbc.eq.REFLECT_EVEN) then
            if ((xsten(-1,dir+1).lt.zero).and. &
                (xsten(1,dir+1).gt.zero)) then
             at_RZ_boundary=1
            else if ((xsten(-1,dir+1).gt.zero).and. &
                     (xsten(1,dir+1).gt.zero)) then
             ! do nothing
            else
             print *,"xsten bust eval_face_coeff"
             stop
            endif
           else if ((local_presbc.eq.INT_DIR).or. &
                    (local_presbc.eq.EXT_DIR).or. &
                    (local_presbc.eq.REFLECT_ODD).or. &
                    (local_presbc.eq.FOEXTRAP)) then
            ! do nothing
           else
            print *,"local_presbc invalid"
            stop
           endif
          else if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           ! do nothing
          else
           print *,"levelrz invalid"
           stop
          endif

         else if ((side.eq.0).or. & !fablo<inorm<fabhi+1
                  (side.eq.2)) then !inorm=fabhi(dir+1)+1
          ! do nothing
         else
          print *,"side invalid; side=",side
          stop
         endif
        else if (dir.eq.1) then
         ! do nothing
        else
         print *,"dir invalid"
         stop
        endif
       else if (SDIM.eq.3) then
        ! do nothing
       else
        print *,"dimension bust"
        stop
       endif

       local_wt(veldir)=zero
       cc_group=cc
       dd_group=dd

        ! SOLVETYPE_PRES, 
        ! SOLVETYPE_PRESGRAVITY, 
        ! SOLVETYPE_INITPROJ
       if (project_option_projectionF(project_option).eq.1) then

        if (project_option.eq.SOLVETYPE_PRES) then!regular pressure projection
         cc_group=cc*cc_ice
        else if (project_option.eq.SOLVETYPE_INITPROJ) then!initial projection
         cc_group=cc*cc_ice_mask
        else if (project_option.eq.SOLVETYPE_PRESGRAVITY) then!grav projection
         cc_group=cc ! we do not mask off the ice or "FSI is rigid" regions
        else
         print *,"project_option invalid eval_face_coeff"
         stop
        endif

        if (nsolve.ne.1) then
         print *,"nsolve invalid for projection"
         stop
        endif

        if (at_RZ_boundary.eq.1) then
         local_wt(veldir)=zero
        else if (at_RZ_boundary.eq.0) then
         if ((dd_group.gt.zero).and. &
             (cc_group.ge.zero).and. &
             (cc_group.le.one)) then
          local_wt(veldir)=dd_group*cc_group
          if (side.eq.0) then
           ! do nothing
          else if ((side.eq.1).or.(side.eq.2)) then
           if ((local_presbc.eq.REFLECT_EVEN).or. &
               (local_presbc.eq.FOEXTRAP)) then
            local_wt(veldir)=zero
           else if ((local_presbc.eq.INT_DIR).or. &
                    (local_presbc.eq.EXT_DIR)) then
            ! do nothing
           else
            print *,"local_presbc invalid"
            stop
           endif
          else
           print *,"side invalid"
           stop
          endif
         else
          print *,"put breakpoint here to see the caller"
          print *,"dd_group or cc_group invalid1"
          print *,"dd_group= ",dd_group
          print *,"cc_group= ",cc_group
          print *,"cc,cc_ice,cc_ice_mask,dd,nsolve,dir,side ", &
             cc,cc_ice,cc_ice_mask,dd,nsolve,dir,side
          print *,"level,finest_level ", &
            level,finest_level
          print *,"at_RZ_boundary ",at_RZ_boundary
          print *,"local_presbc ",local_presbc
          print *,"project_option= ",project_option
          stop
         endif
        else
         print *,"at_RZ_boundary invalid;eval_face_coeff"
         stop
        endif

       else if (project_option.eq.SOLVETYPE_PRESEXTRAP) then 

        if (nsolve.ne.1) then
         print *,"nsolve invalid for pressure extrapolation"
         stop
        endif
        cc_group=cc*cc_ice

        if (at_RZ_boundary.eq.1) then
         local_wt(veldir)=zero
        else if (at_RZ_boundary.eq.0) then
         if ((dd_group.gt.zero).and. &
             (cc_group.ge.zero).and. &
             (cc_group.le.one)) then
          if ((cc_group.gt.zero).and. &
              (cc_group.le.one)) then
           local_wt(veldir)=zero
          else if (cc_group.eq.zero) then
           local_wt(veldir)=one
          else
           print *,"cc_group invalid SOLVETYPE_PRESEXTRAP: ",cc_group
           stop
          endif
          if (side.eq.0) then
           ! do nothing
          else if ((side.eq.1).or.(side.eq.2)) then
           if (local_presbc.eq.INT_DIR) then
            ! do nothing (periodic BC)
           else if (local_presbc.eq.EXT_DIR) then !pressure extrap case
            ! do nothing (same BC as regular pressue)
           else if ((local_presbc.eq.REFLECT_EVEN).or. &
                    (local_presbc.eq.FOEXTRAP)) then
            local_wt(veldir)=zero ! pressure extrap case
           else
            print *,"local_presbc invalid"
            stop
           endif
          else
           print *,"side invalid"
           stop
          endif
         else
          print *,"dd_group or cc_group invalid2 SOLVETYPE_PRESEXTRAP"
          print *,"dd_group= ",dd_group
          print *,"cc_group= ",cc_group
          print *,"project_option=",project_option
          print *,"subroutine eval_face_coeff"
          stop
         endif
        else
         print *,"at_RZ_boundary invalid"
         stop
        endif

       else if ((project_option.eq.SOLVETYPE_HEAT).or. & ! temperature
                ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
                 (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then

        if (nsolve.ne.1) then
         print *,"nsolve invalid for temperature/species"
         stop
        endif
        cc_group=one  ! rely on solid coefficient at solid
        if ((dd_group.ge.zero).and. &
            (cc_group.eq.one)) then
         local_wt(veldir)=dd_group*cc_group
         if (side.eq.0) then
          ! do nothing
         else if ((side.eq.1).or.(side.eq.2)) then
          if (local_presbc.eq.INT_DIR) then
           ! do nothing
          else if (local_presbc.eq.EXT_DIR) then
           ! do nothing (first order Dirichlet for low order BC)
           ! 2nd order or higher for spectral element method.
          else if (local_presbc.eq.REFLECT_EVEN) then
           local_wt(veldir)=zero
          else if (local_presbc.eq.FOEXTRAP) then
           local_wt(veldir)=zero
          else
           print *,"local_presbc invalid"
           stop
          endif
         else
          print *,"side invalid"
          stop
         endif
        else
         print *,"dd_group or cc_group invalid3"
         print *,"dd_group= ",dd_group
         print *,"cc_group= ",cc_group
         print *,"project_option=",project_option
         stop
        endif

       else if (project_option.eq.SOLVETYPE_VISC) then ! viscosity

        if (nsolve.ne.SDIM) then
         print *,"nsolve invalid for viscosity"
         stop
        endif
        dd_group=dd*visc_coef
        cc_group=one
        ddfactor=one
        if ((dd_group.ge.zero).and. &
            (cc_group.eq.one)) then
         if (uncoupled_viscosity.eq.0) then
          if (dir+1.eq.veldir) then
           ddfactor=two
          endif
         else if (uncoupled_viscosity.eq.1) then
          ! do nothing
         else
          print *,"uncoupled_viscosity invalid"
          stop
         endif

         local_wt(veldir)=dd_group*cc_group*ddfactor

         if (side.eq.0) then
          ! do nothing
         else if ((side.eq.1).or.(side.eq.2)) then
          if (local_presbc.eq.INT_DIR) then
           ! do nothing
          else if (local_presbc.eq.EXT_DIR) then
           ! do nothing 
           ! a) 2nd order BC if standard FVM 
           !    discretization (fort_face_gadients)
           ! b) 2nd order or higher BC if SEM discretization.
          else if (local_presbc.eq.REFLECT_ODD) then
           ! do nothing
          else if (local_presbc.eq.REFLECT_EVEN) then
           local_wt(veldir)=zero
          else if (local_presbc.eq.FOEXTRAP) then
           local_wt(veldir)=zero
          else
           print *,"local_presbc invalid"
           stop
          endif
         else
          print *,"side invalid"
          stop
         endif
        else
         print *,"dd_group or cc_group invalid4"
         print *,"dd_group= ",dd_group
         print *,"cc_group= ",cc_group
         print *,"project_option=",project_option
         stop
        endif

       else
        print *,"project_option invalid eval_face_coeff"
        stop
       endif

      else 
       print *,"coefficients bust"
       print *,"cc=",cc
       print *,"cc_ice=",cc_ice
       print *,"cc_ice_mask=",cc_ice_mask
       print *,"dd=",dd
       print *,"visc_coef=",visc_coef
       print *,"nsolve=",nsolve
       print *,"veldir=",veldir
       print *,"dir=",dir
       print *,"side=",side
       print *,"project_option=",project_option
       print *,"uncoupled_viscosity=",uncoupled_viscosity
       stop
      endif

      return
      end subroutine eval_face_coeff

! OP_PRESGRAD_MAC (0)
! operation_flag=0  pressure gradient on MAC grid
! OP_POTGRAD_TO_MAC (2)
! operation_flag=2  potential gradient on MAC grid
! OP_UNEW_CELL_TO_MAC (3)
! operation_flag=3  unew^MAC=unew^CELL->MAC
! OP_UNEW_USOL_MAC_TO_MAC (4)
! operation_flag=4  unew^MAC=uSOLID^MAC or uFLUID^MAC (not used here)
! OP_UMAC_PLUS_VISC_CELL_TO_MAC (5)
! operation_flag=5  unew^MAC=unew^MAC+beta diffuse_ref^CELL->MAC
! OP_UGRAD_MAC (6)
! operation_flag=6  evaluate tensor values.
!   (called from FACE_GRADIENTS)
! OP_ISCHEME_MAC (7)
! operation_flag=7  advection.
! OP_UGRAD_COUPLING_MAC (8)
! operation_flag=8  reserved for coupling terms in CROSSTERM
!   (SEM_CELL_TO_MAC not called with 
!    operation_flag==OP_UGRAD_COUPLING_MAC)
! OP_U_SEM_CELL_MAC_TO_MAC (11)
! operation_flag=11 unew^MAC=uold^MAC +(unew^cell-uold^cell)^{cell->MAC}
! SEM="spectral element method"
! CELL=Gauss-Legendre,Gauss-Legendre,Gauss-Legendre input data
! MAC=Gauss Lobatto Legendre,Gauss-Legendre, Gauss-Legendre
!  or
! MAC=Gauss-Legendre,Gauss Lobatto Legendre,Gauss-Legendre
!  or
! MAC=Gauss-Legendre,Gauss-Legendre,Gauss Lobatto Legendre
      subroutine SEM_CELL_TO_MAC( &
       ncomp_xp, & !number of amrsync components if op=0,3,5,6,7,9,10,11
       simple_AMR_BC_flag_in, &
       level, &
       finest_level, &
       operation_flag, &
       energyflag, &
       project_option, &
       beta, &
       visc_coef, &
       time, &
       dt, &
       i,j,k, &
       tilelo,tilehi, &
       fablo,fabhi, &
       xlo, &
       dx, &
       dir, &
       bfact,bfact_c,bfact_f, &
       presbc_in, &
       velbc_in, &
       scomp, &
       scomp_bc, &
       dcomp, &
       ncomp_dest, &
       ncomp_source, &
       ncomp_xgp, &
       ncphys, &
       spectral_loop, &
       ncfluxreg, &
       semflux, &
       maskCF, & !maskCF=1.0 at interior fine bc ghost cells
       maskcov, & ! 1=not cov. or outside domain  
       vel, & 
       pres, & 
       den, & 
       xface, & 
       ! Umac_old if OP_UMAC_PLUS_VISC_CELL_TO_MAC or 
       !             OP_U_SEM_CELL_MAC_TO_MAC 
       xgp, &  
       xcut, & 
       xp, &  ! holds amrsync if op==0,3,5,6,7,9,10,11
       xvel, & 
       maskSEM )
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: ncomp_xp
      integer, INTENT(in) :: simple_AMR_BC_flag_in
      integer :: simple_AMR_BC_flag
      integer :: local_AMR_BC_flag
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: beta,visc_coef
      integer, INTENT(in) :: operation_flag
      integer, INTENT(in) :: energyflag
      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: dir
      integer, INTENT(in) :: bfact,bfact_c,bfact_f
      integer, INTENT(in) :: scomp,scomp_bc
      integer, INTENT(in) :: dcomp
      integer, INTENT(in) :: ncomp_dest
      integer, INTENT(in) :: ncomp_source
      integer, INTENT(in) :: ncomp_xgp
      integer, INTENT(in) :: ncphys
      integer, INTENT(in) :: spectral_loop
      integer, INTENT(in) :: ncfluxreg
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: presbc_in(SDIM,2,num_materials*num_state_material)
      integer, INTENT(in) :: velbc_in(SDIM,2,SDIM)
        ! INTENT(in) means the pointer cannot be reassigned.
        ! The data itself inherits the INTENT attribute from the
        ! target.
      real(amrex_real), INTENT(in), pointer :: semflux(D_DECL(:,:,:),:)
      ! 1=fine-fine  0=coarse-fine
      real(amrex_real), INTENT(in), pointer :: maskCF(D_DECL(:,:,:)) 
      real(amrex_real), INTENT(in), pointer :: maskcov(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), pointer :: maskSEM(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), pointer :: vel(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: pres(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: den(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: xface(D_DECL(:,:,:),:)
       ! INTENT(in) means the pointer cannot be reassigned.
       ! The data itself inherits the INTENT attribute from the
       ! target.
       !xgp is usually a dest variable except that it holds umac_old if 
       !op==11 or op 5
      real(amrex_real), INTENT(in), pointer :: xgp(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: xcut(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: xp(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: xvel(D_DECL(:,:,:),:)

      integer local_bctype(2)
      integer local_bctype_den(2)
      real(amrex_real) x_sep(2)
      real(amrex_real) local_bcval(2)
      real(amrex_real) local_bcval_den(2)
      integer cen_maskSEM
      real(amrex_real) local_vel(0:bfact)
      real(amrex_real) RRface(0:bfact)
      real(amrex_real) local_data(bfact)
      real(amrex_real) local_data_side(2)
      real(amrex_real) local_data_den(bfact)
      real(amrex_real) local_data_side_den(2)
      real(amrex_real) local_grad(bfact+1)
      real(amrex_real) local_grad_den(bfact+1)
      real(amrex_real) local_interp(bfact+1)
      real(amrex_real) local_interp_den(bfact+1)
      integer ii,jj,kk
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_coarse(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_fine(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_face(-nhalf:nhalf,SDIM)
      integer nc
      integer side
       ! 0=fine-fine neighbor  1=fine-wall neighbor  
       ! -1=fine (current) -coarse neighbor -2=coarse (current) -fine neighbor
      integer nbr_outside_domain_flag(2)
      integer cen_outside_domain_flag
      integer strip_outside_fab_flag
      integer isten
      integer dir2
      integer sideidx(SDIM)
      integer i_out,j_out,k_out
      integer mask_out
      integer local_maskcov
      integer local_maskCF
      integer shared_face ! in: SEM_CELL_TO_MAC
      real(amrex_real) shared_face_value
      integer i_in,j_in,k_in
      integer ic,jc,kc
      integer icoarse,jcoarse,kcoarse
      integer iface_out,jface_out,kface_out
      integer ifine,jfine,kfine
      integer ibase
      integer ngroup,fluxbase
      integer indexlo(SDIM)
      integer indexhi(SDIM)
      integer indexmid(SDIM)
      integer index_edge(SDIM)
      integer index_opp(SDIM)
      real(amrex_real) templocal
      integer test_maskSEM
      real(amrex_real) shared_xcut
      integer nbase
      integer testbc
      real(amrex_real) problo(SDIM),probhi(SDIM),problen(SDIM)
      real(amrex_real) dx_c(SDIM),dx_f(SDIM)
      integer domlo(SDIM)
      integer bctype_tag
      real(amrex_real) inside_flux,outside_flux
      real(amrex_real) udotn_boundary

      simple_AMR_BC_flag=simple_AMR_BC_flag_in
      if (1.eq.0) then
       simple_AMR_BC_flag=1
      endif

      problo(1)=problox 
      problo(2)=probloy
      probhi(1)=probhix 
      probhi(2)=probhiy
      if (SDIM.eq.3) then
       problo(SDIM)=probloz
       probhi(SDIM)=probhiz
      endif

      do dir2=1,SDIM
       problen(dir2)=probhi(dir2)-problo(dir2)
       if (problen(dir2).gt.zero) then
        ! do nothing
       else
        print *,"problen invalid"
        stop
       endif
       dx_c(dir2)=dx(dir2)
       dx_f(dir2)=dx(dir2)
       if (level.gt.0) then
        dx_c(dir2)=two*dx(dir2)
       endif
       if (level.lt.finest_level) then
        dx_f(dir2)=half*dx(dir2)
       endif
       domlo(dir2)=0
      enddo ! dir2=1..sdim

      if ((simple_AMR_BC_flag.eq.0).or. &
          (simple_AMR_BC_flag.eq.1)) then
       ! do nothing
      else
       print *,"simple_AMR_BC_flag invalid"
       stop
      endif

      if (ncomp_xp.lt.1) then
       print *,"ncomp_xp invalid(10) ",ncomp_xp
       stop
      endif
      if (ncomp_xgp.lt.1) then
       print *,"ncomp_xgp invalid"
       stop
      endif
 
      if (ncphys.lt.1) then
       print *,"ncphys invalid"
       stop
      endif
      if (bfact.lt.2) then
       print *,"bfact invalid200"
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

      if ((level.lt.0).or. &
          (level.gt.finest_level)) then
       print *,"level invalid sem cell to mac"
       stop
      endif

      if ((dir.lt.1).or.(dir.gt.SDIM)) then
       print *,"dir invalid sem cell to mac1"
       stop
      endif

      ngroup=(ncfluxreg/SDIM)
      if (ngroup*SDIM.ne.ncfluxreg) then
       print *,"ncfluxreg invalid1 ",ncfluxreg
       stop
      endif
      fluxbase=(dir-1)*ngroup

      cen_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
      local_maskcov=NINT(maskcov(D_DECL(i,j,k)))
      if (local_maskcov.ne.1) then
       print *,"local_maskcov invalid in SEM_CELL_TO_MAC"
       stop
      endif
      if ((cen_maskSEM.lt.0).or.(cen_maskSEM.gt.num_materials)) then
       print *,"cen_maskSEM invalid"
       stop
      endif

      if (operation_flag.eq.OP_UGRAD_MAC) then ! tensor

       if (ncomp_xgp.ne.AMREX_SPACEDIM_SQR) then
        print *,"ncomp_xgp invalid1"
        stop
       endif
       if (ncomp_xp.ne.SDIM) then
        print *,"ncomp_xp invalid (it is supposed to be sdim)",ncomp_xp
        print *,"operation_flag==OP_UGRAD_MAC"
        stop
       endif
       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_UGRAD_MAC"
        stop
       endif
       if (project_option.ne.SOLVETYPE_VISC) then
        print *,"project_option.ne.SOLVETYPE_VISC; SEM_CELL_TO_MAC"
        stop
       endif
        ! number of components for flux synchronization.
       if (ncfluxreg.ne.AMREX_SPACEDIM_SQR) then
        print *,"ncfluxreg invalid2 ",ncfluxreg
        stop
       endif
! compute u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  

       nbase=dcomp-1

       if ((nbase.ne.0).and.(nbase.ne.SDIM).and. &
           (nbase.ne.(SDIM-1)*SDIM)) then
        print *,"nbase invalid"
        stop
       endif

       if (scomp.ne.1) then
        print *,"scomp invalid"
        stop
       endif
       if (ncomp_source.ne.SDIM) then
        print *,"ncomp_source invalid"
        stop
       endif
       if (ncomp_dest.ne.SDIM) then
        print *,"ncomp_dest invalid"
        stop
       endif
       if (scomp_bc.ne.1) then
        print *,"scomp_bc invalid"
        stop
       endif

      else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

       if (ncomp_xp.ne.NFLUXSEM) then
        print *,"ncomp_xp invalid(11) ",ncomp_xp
        print *,"operation_flag==OP_ISCHEME_MAC"
        stop
       endif
       if (ncomp_xgp.ne.NFLUXSEM) then
        print *,"expecting ncomp_xgp==NFLUXSEM invalid2"
        print *,"operation_flag==OP_ISCHEME_MAC"
        stop
       endif
       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_ISCHEME_MAC"
        stop
       endif
       if (ncfluxreg.ne.SDIM*NFLUXSEM) then
        print *,"ncfluxreg invalid operation_flag.eq.OP_ISCHEME_MAC"
        stop
       endif
       if ((scomp.ne.1).or. &
           (dcomp.ne.1).or. &
           (ncomp_dest.ne.ncphys).or. &
           (ncphys.ne.NFLUXSEM).or. &
           (ncomp_source.ne.SDIM).or. &
           (scomp_bc.ne.1)) then
        print *,"parameters invalid for op=7"
        print *,"operation_flag==OP_ISCHEME_MAC"
        stop
       endif
       if ((cen_maskSEM.ge.1).and.(cen_maskSEM.le.num_materials)) then
        ! do nothing
       else
        print *,"cen_maskSEM invalid"
        stop
       endif

      else if (operation_flag.eq.OP_PRESGRAD_MAC) then ! grad p

       if (ncomp_xp.eq.1) then
        ! do nothing
       else
        print *,"expecting ncomp_xp=1 if OP_PRESGRAD_MAC: ",ncomp_xp
        stop
       endif

       if (ncomp_xgp.ne.1) then
        print *,"ncomp_xgp invalid3"
        print *,"expecting ncomp_xgp=1 if OP_PRESGRAD_MAC: ",ncomp_xgp
        stop
       endif
       if ((energyflag.ne.SUB_OP_FOR_MAIN).and. & ! regular solver
           (energyflag.ne.SUB_OP_FOR_SDC)) then  ! for SDC
        print *,"energyflag invalid OP_PRESGRAD_MAC: ",energyflag
        stop
       endif
       if ((scomp.ne.1).or. &
           (dcomp.ne.1)) then
        print *,"parameters invalid for op=0=OP_PRESGRAD_MAC"
        print *,"scomp=",scomp
        print *,"dcomp=",dcomp
        stop
       endif
       if ((ncomp_dest.ne.1).or. &
           (ncomp_source.ne.1).or. &
           (scomp_bc.ne.1)) then
        print *,"parameters invalid for op=0=OP_PRESGRAD_MAC"
        print *,"ncomp_dest=",ncomp_dest
        print *,"ncomp_source=",ncomp_source
        print *,"scomp_bc=",scomp_bc
        stop
       endif

      else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then ! grad ppot

       if (ncomp_xp.ne.1) then
        print *,"ncomp_xp invalid5 OP_POTGRAD_TO_MAC ",ncomp_xp
        print *,"operation_flag.eq.OP_POTGRAD_TO_MAC"
        stop
       endif

       if (ncomp_xgp.ne.1) then
        print *,"ncomp_xgp invalid OP_POTGRAD_TO_MAC ",ncomp_xgp
        stop
       endif

       if (energyflag.ne.SUB_OP_FORCE_MASK_BASE+3) then
        print *,"energyflag invalid OP_POTGRAD_TO_MAC"
        stop
       endif

       if ((ncomp_dest.ne.1).or. &
           (ncomp_source.ne.1).or. &
           (scomp.ne.1).or. &
           (dcomp.ne.1).or. &
           (scomp_bc.ne.1)) then
        print *,"parameters invalid for op=2 OP_POTGRAD_TO_MAC"
        stop
       endif

      else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & !unew^CELL->MAC
               (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC).or. &
               (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC)) then 

       if (ncomp_xgp.ne.1) then
        print *,"expecting ncomp_xgp=1 OP_UNEW,OP_U_COMP,OP_UMAC"
        stop
       endif
       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_U ETC CELL/MAC to MAC"
        stop
       endif

       if ((scomp.ne.dir).or. &
           (dcomp.ne.1)) then
        print *,"parameters invalid op=3,10,11,5, OP_UNEW,OP_U_COMP,OP_UMAC"
        stop
       endif
       if ((ncomp_dest.ne.1).or. &
           (ncomp_source.ne.1).or. &
           (scomp_bc.ne.dir)) then
        print *,"parameters invalid for op=3 OP_UNEW,OP_U_COMP,OP_UMAC"
        stop
       endif

      else
       print *,"operation_flag invalid20"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.1) then
       ii=1
       if ((i/bfact)*bfact.ne.i) then
        print *,"i invalid"
        stop
       endif
       if (i.lt.0) then
        print *,"i invalid"
        stop
       endif
      else if (dir.eq.2) then
       jj=1
       if ((j/bfact)*bfact.ne.j) then
        print *,"j invalid"
        stop
       endif
       if (j.lt.0) then
        print *,"j invalid"
        stop
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       kk=1
       if ((k/bfact)*bfact.ne.k) then
        print *,"k invalid SEM_CELL_TO_MAC"
        stop
       endif
       if (k.lt.0) then
        print *,"k invalid SEM_CELL_TO_MAC"
        stop
       endif
      else
       print *,"dir invalid sem cell to mac2"
       stop
      endif

      do dir2=1,SDIM
       if (fablo(dir2).lt.0) then
        print *,"fablo invalid"
        stop
       endif
      enddo ! dir2

      if ((dir.lt.1).or.(dir.gt.SDIM).or. &
          (bfact.lt.2)) then
       print *,"SEM_CELL_TO_MAC parms corrupt"
       stop
      endif

       ! do nothing unless the strip is a spectral element strip.
      if ((cen_maskSEM.ge.1).and. &
          (cen_maskSEM.le.num_materials)) then
       ! do nothing
      else if (cen_maskSEM.eq.0) then
       ! do nothing
      else
       print *,"cen_maskSEM invalid SEM_CELL_TO_MAC: ",cen_maskSEM
       print *,"i,j,k : ",i,j,k
       stop
      endif

       ! ncomp_dest=SDIM if OP_UGRAD_MAC  (du/dn, dv/dn, dw/dn)
       ! ncomp_dest=1 if OP_PRESGRAD_MAC
      do nc=1,ncomp_dest

        ! do nothing if element is not a spectral element
       if ((cen_maskSEM.ge.1).and. &
           (cen_maskSEM.le.num_materials)) then

        indexlo(1)=i
        indexlo(2)=j
        if (SDIM.eq.3) then
         indexlo(SDIM)=k
        endif

        do dir2=1,SDIM
         indexhi(dir2)=indexlo(dir2)
        enddo ! dir2
        indexhi(dir)=indexlo(dir)+bfact-1

        strip_outside_fab_flag=0

        cen_outside_domain_flag=0

        do dir2=1,SDIM
 
         side=1
         if (indexlo(dir2).eq.fablo(dir2)-1) then
          if (dir2.eq.dir) then
           print *,"indexlo(dir) invalid; strip must at least"
           print *,"be alongside a strip that is in the domain"
           stop
          endif

          strip_outside_fab_flag=1

          if (velbc_in(dir2,side,dir2).ne.INT_DIR) then
           cen_outside_domain_flag=1
          endif
         else if (indexlo(dir2).ge.fablo(dir2)) then
          ! do nothing
         else
          print *,"indexlo(dir2)<fablo(dir2)-1"
          print *,"dir,dir2,indexlo(dir2),fablo(dir2) ", &
               dir,dir2,indexlo(dir2),fablo(dir2)
          stop
         endif

         side=2
         if (indexhi(dir2).eq.fabhi(dir2)+1) then
          if (dir2.eq.dir) then
           print *,"indexhi(dir) invalid; strip must at least"
           print *,"be alongside a strip that is in the domain"
           stop
          endif

          strip_outside_fab_flag=1

          if (velbc_in(dir2,side,dir2).ne.INT_DIR) then
           cen_outside_domain_flag=1
          endif
         else if (indexhi(dir2).le.fabhi(dir2)) then
          ! do nothing
         else
          print *,"indexhi(dir2)>fabhi(dir2)+1"
          print *,"dir,dir2,indexhi(dir2),fabhi(dir2) ", &
               dir,dir2,indexhi(dir2),fabhi(dir2)
          stop
         endif

        enddo ! dir2=1..sdim

        if (cen_outside_domain_flag.eq.1) then

         ! do nothing

        else if (cen_outside_domain_flag.eq.0) then

         ibase=num_state_material*(cen_maskSEM-1)
     
         do side=1,2

          nbr_outside_domain_flag(side)=0

          do dir2=1,SDIM
           sideidx(dir2)=indexlo(dir2)
          enddo

          if (side.eq.1) then

           sideidx(dir)=indexlo(dir)-1

           local_maskcov= &
            NINT(maskcov(D_DECL(sideidx(1),sideidx(2),sideidx(SDIM))))
           
           do dir2=1,SDIM
            if (sideidx(dir2).eq.fablo(dir2)-1) then
             testbc=velbc_in(dir2,side,dir2)
             if ((testbc.eq.INT_DIR).and.(dir2.eq.dir)) then
              local_maskCF= &
               NINT(maskCF(D_DECL(sideidx(1),sideidx(2),sideidx(SDIM))))
              if (local_maskCF.eq.1) then ! fine-fine
               ! do nothing
              else if (local_maskCF.eq.0) then ! fine (current) -coarse
               nbr_outside_domain_flag(side)=-1
              else
               print *,"local_maskCF invalid"
               stop
              endif
             else if ((testbc.eq.INT_DIR).and.(dir2.ne.dir)) then
              nbr_outside_domain_flag(side)=0
              if (operation_flag.eq.OP_UGRAD_MAC) then
               ! do nothing (we expect this case for grad U)
              else
               print *,"operation_flag invalid"
               stop
              endif
             else if ((testbc.eq.EXT_DIR).or. &
                      (testbc.eq.REFLECT_EVEN).or. &
                      (testbc.eq.REFLECT_ODD).or. &
                      (testbc.eq.FOEXTRAP)) then
              nbr_outside_domain_flag(side)=1
             else
              print *,"testbc invalid"
              stop 
             endif
            else if (sideidx(dir2).ge.fablo(dir2)) then
             if (local_maskcov.eq.1) then
              ! do nothing
             else if (local_maskcov.eq.0) then ! coarse(current) - fine
              nbr_outside_domain_flag(side)=-2
             else
              print *,"local_maskcov invalid"
              stop
             endif
            else
             print *,"sideidx(dir2) < fablo(dir2)-1"
             print *,"dir,dir2,sideidx(dir2),fablo(dir2) ", &
                  dir,dir2,sideidx(dir2),fablo(dir2)
             stop 
            endif
           enddo ! dir2=1..sdim

           i_out=indexlo(1)-ii
           j_out=indexlo(2)-jj
           k_out=indexlo(SDIM)-kk

           iface_out=indexlo(1)
           jface_out=indexlo(2)
           kface_out=indexlo(SDIM)

          else if (side.eq.2) then

           sideidx(dir)=indexhi(dir)+1

           local_maskcov= &
            NINT(maskcov(D_DECL(sideidx(1),sideidx(2),sideidx(SDIM))))

           do dir2=1,SDIM
            if (sideidx(dir2).eq.fabhi(dir2)+1) then
             testbc=velbc_in(dir2,side,dir2)
             if ((testbc.eq.INT_DIR).and.(dir2.eq.dir)) then
              local_maskCF= &
               NINT(maskCF(D_DECL(sideidx(1),sideidx(2),sideidx(SDIM))))
              if (local_maskCF.eq.1) then ! fine-fine
               ! do nothing
              else if (local_maskCF.eq.0) then ! fine (current) -coarse
               nbr_outside_domain_flag(side)=-1
              else
               print *,"local_maskCF invalid"
               stop
              endif
             else if ((testbc.eq.INT_DIR).and.(dir2.ne.dir)) then
              nbr_outside_domain_flag(side)=0
              if (operation_flag.eq.OP_UGRAD_MAC) then
               ! do nothing (we expect this case for grad U)
              else
               print *,"operation_flag invalid"
               stop
              endif
             else if ((testbc.eq.EXT_DIR).or. &
                      (testbc.eq.REFLECT_EVEN).or. &
                      (testbc.eq.REFLECT_ODD).or. &
                      (testbc.eq.FOEXTRAP)) then
              nbr_outside_domain_flag(side)=1
             else
              print *,"testbc invalid"
              stop 
             endif
            else if (sideidx(dir2).le.fabhi(dir2)) then
             if (local_maskcov.eq.1) then
              ! do nothing
             else if (local_maskcov.eq.0) then ! coarse(current) - fine
              nbr_outside_domain_flag(side)=-2
             else
              print *,"local_maskcov invalid"
              stop
             endif
            else
             print *,"sideidx(dir2) > fabhi(dir2)+1"
             print *,"dir,dir2,sideidx(dir2),fabhi(dir2) ", &
                  dir,dir2,sideidx(dir2),fabhi(dir2)
             stop 
            endif
           enddo ! dir2=1..sdim

           i_out=indexhi(1)+ii
           j_out=indexhi(2)+jj
           k_out=indexhi(SDIM)+kk

           iface_out=i_out
           jface_out=j_out
           kface_out=k_out
          else 
           print *,"side invalid"
           stop
          endif

          call gridsten(xsten,xlo, &
           i_out,j_out,k_out, &
           fablo,bfact,dx,nhalf)

           ! dir=1..sdim
          call gridstenMAC(xsten_face,xlo, &
           iface_out,jface_out,kface_out, &
           fablo,bfact,dx,nhalf,dir-1)

           ! fine-fine boundary or element properly contained in grid.
          if (nbr_outside_domain_flag(side).eq.0) then

           local_bctype(side)=SEM_INTERIOR 

           local_bcval(side)=zero
            
           ic=i_out
           jc=j_out
           kc=k_out

             !advection(values outside elem)
           if (operation_flag.eq.OP_ISCHEME_MAC) then

            templocal=den(D_DECL(ic,jc,kc),ibase+ENUM_TEMPERATUREVAR+1)

            if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then
             ! u_{i+1/2}S_{i+1/2}-u_{i-1/2}S_{i-1/2}=
             ! (S_{i+1/2}+S_{i-1/2})(u_{i+1/2}-u_{i-1/2})/2 +
             ! (u_{i+1/2}+u_{i-1/2})(S_{i+1/2}-S_{i-1/2})/2
             ! u dot grad u = div(umac u)-I(umac) div umac
             local_data_side(side)= &
               vel(D_DECL(ic,jc,kc),nc)   ! velocity

             ! I(umac) dot grad T
            else if (nc.eq.SEM_T+1) then 
             local_data_side(side)=templocal ! temperature
            else
             print *,"nc invalid"
             stop
            endif

           else if (operation_flag.eq.OP_UGRAD_MAC) then ! tensor derivatives

            if ((nc.ge.1).and.(nc.le.SDIM)) then
             local_data_side(side)=vel(D_DECL(ic,jc,kc),nc)
            else
             print *,"nc invalid"
             stop
            endif

           else if (operation_flag.eq.OP_PRESGRAD_MAC) then!MAC pres. gradient
            if (nc.eq.1) then
             if (scomp.eq.1) then
              local_data_side(side)=pres(D_DECL(ic,jc,kc),1)
             else
              print *,"scomp invalid OP_PRESGRAD_MAC"
              stop
             endif
            else
             print *,"nc invalid OP_PRESGRAD_MAC"
             stop
            endif
           else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 
            if (nc.eq.1) then
             if (scomp.eq.1) then
              local_data_side(side)=pres(D_DECL(ic,jc,kc),1)
              local_data_side_den(side)=den(D_DECL(ic,jc,kc),1)
              if ((abs(local_data_side(side)).lt.1.0D+20).and. &
                  (local_data_side_den(side).gt.zero).and. &
                  (abs(local_data_side_den(side)).lt.1.0D+20)) then
               ! do nothing
              else
               print *,"local_data_side or local_data_side_den invalid"
               stop
              endif
             else
              print *,"scomp invalid OP_POTGRAD_TO_MAC"
              stop
             endif
            else
             print *,"nc invalid OP_POTGRAD_TO_MAC"
             stop
            endif
           else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
                    (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC).or. &
                    (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC)) then 

             !if OP_UNEW_CELL_TO_MAC:
             !primary_vel_data="vel"=CURRENT_CELL_VEL_MF; 
             !if OP_UMAC_PLUS_VISC_CELL_TO_MAC:
             !primary_vel_data="vel"=idx_velcell;  (increment)
             !if OP_U_SEM_CELL_MAC_TO_MAC:
             !primary_vel_data="vel"=DELTA_CELL_VEL_MF; 

            if (nc.eq.1) then
             if (scomp.eq.dir) then
              local_data_side(side)=vel(D_DECL(ic,jc,kc),scomp)
             else
              print *,"scomp invalid(OP_UNEW, or OP_U_COMP, or OP_UMAC)"
              stop
             endif
            else
             print *,"nc invalid(OP_UNEW, or OP_U_COMP, or OP_UMAC)"
             stop
            endif
           else
            print *,"operation_flag invalid21"
            stop
           endif

           ! element touches the domain wall
          else if (nbr_outside_domain_flag(side).eq.1) then 

           local_data_side(side)=zero

           if (velbc_in(dir,side,dir).eq.INT_DIR) then
            print *,"velbc_in bust "
            print *,"cen_outside_domain_flag= ",cen_outside_domain_flag
            print *,"nbr_outside_domain_flag= ",nbr_outside_domain_flag(side)
            stop
           endif

           if (operation_flag.eq.OP_UGRAD_MAC) then ! grad U

            if (velbc_in(dir,side,nc).eq.REFLECT_EVEN) then
             local_bctype(side)=SEM_REFLECT_EVEN
             local_bcval(side)=zero
            else if (velbc_in(dir,side,nc).eq.FOEXTRAP) then
             local_bctype(side)=SEM_NEUMANN
             local_bcval(side)=zero
            else if (velbc_in(dir,side,nc).eq.REFLECT_ODD) then
             local_bctype(side)=SEM_REFLECT_ODD
             local_bcval(side)=zero
            else if (velbc_in(dir,side,nc).eq.EXT_DIR) then
             local_bctype(side)=SEM_DIRICHLET

             local_bcval(side)=vel(D_DECL(i_out,j_out,k_out),nc)
            else
             print *,"velbc_in is corrupt"
             stop
            endif

           else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
                    (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC).or. &
                    (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC)) then 

            if (velbc_in(dir,side,scomp_bc).eq.REFLECT_EVEN) then
             local_bctype(side)=SEM_REFLECT_EVEN
             local_bcval(side)=zero
            else if (velbc_in(dir,side,scomp_bc).eq.FOEXTRAP) then
             local_bctype(side)=SEM_NEUMANN
             local_bcval(side)=zero
            else if (velbc_in(dir,side,scomp_bc).eq.REFLECT_ODD) then
             local_bctype(side)=SEM_REFLECT_ODD
             local_bcval(side)=zero
            else if (velbc_in(dir,side,scomp_bc).eq.EXT_DIR) then
             local_bctype(side)=SEM_DIRICHLET
             local_bcval(side)=vel(D_DECL(i_out,j_out,k_out),scomp)
            else
             print *,"velbc_in is corrupt"
             stop
            endif

           else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection (bc's)

             ! normal points out of the computational domain.
            if (side.eq.1) then
             udotn_boundary=-vel(D_DECL(i_out,j_out,k_out),dir)
            else if (side.eq.2) then
             udotn_boundary=vel(D_DECL(i_out,j_out,k_out),dir)
            else
             print *,"side invalid"
             stop
            endif

            if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then ! velocity

             if (velbc_in(dir,side,nc).eq.REFLECT_EVEN) then
              local_bctype(side)=SEM_REFLECT_EVEN
              local_bcval(side)=zero
             else if (velbc_in(dir,side,nc).eq.FOEXTRAP) then
              local_bctype(side)=SEM_NEUMANN
              local_bcval(side)=zero
             else if (velbc_in(dir,side,nc).eq.REFLECT_ODD) then
              local_bctype(side)=SEM_REFLECT_ODD
              local_bcval(side)=zero
             else if (velbc_in(dir,side,nc).eq.EXT_DIR) then

               ! normal points out of the computational domain.
               ! udotn<0 => characteristics enter into the domain.
              if (udotn_boundary.lt.zero) then
               local_bctype(side)=SEM_DIRICHLET

               local_bcval(side)= &
                 vel(D_DECL(i_out,j_out,k_out),nc) 
              else if (udotn_boundary.ge.zero) then
               local_bctype(side)=SEM_NEUMANN
               local_bcval(side)=zero
              else
               print *,"udotn_boundary invalid"
               stop
              endif

             else
              print *,"velbc_in is corrupt"
              stop
             endif

            else if (nc.eq.SEM_T+1) then ! temperature

             if (presbc_in(dir,side,ibase+ENUM_TEMPERATUREVAR+1).eq. &
                 REFLECT_EVEN) then
              local_bctype(side)=SEM_REFLECT_EVEN
              local_bcval(side)=zero 
             else if (presbc_in(dir,side,ibase+ENUM_TEMPERATUREVAR+1).eq. &
                      FOEXTRAP) then
              local_bctype(side)=SEM_NEUMANN
              local_bcval(side)=zero 
             else if (presbc_in(dir,side,ibase+ENUM_TEMPERATUREVAR+1).eq. &
                      REFLECT_ODD) then
              print *,"cannot have reflect odd BC for temperature"
              stop
              local_bctype(side)=SEM_REFLECT_ODD
              local_bcval(side)=zero 
             else if (presbc_in(dir,side,ibase+ENUM_TEMPERATUREVAR+1).eq. &
                      EXT_DIR) then
              if (udotn_boundary.lt.zero) then !INFLOW BOUNDARY CONDITION
               templocal=den(D_DECL(i_out,j_out,k_out), &
                 ibase+ENUM_TEMPERATUREVAR+1)
               local_bctype(side)=SEM_DIRICHLET
               local_bcval(side)=templocal
              else if (udotn_boundary.ge.zero) then
               local_bctype(side)=SEM_NEUMANN !NUMERICAL BOUNDARY CONDITION
               local_bcval(side)=zero
              else
               print *,"udotn_boundary invalid"
               stop
              endif
             else
              print *,"presbc_in is corrupt"
              stop
             endif
            else
             print *,"nc invalid"
             stop
            endif

           else if ((operation_flag.eq.OP_PRESGRAD_MAC).or. & 
                    (operation_flag.eq.OP_POTGRAD_TO_MAC)) then 

            if (presbc_in(dir,side,1).eq.REFLECT_EVEN) then
             local_bctype(side)=SEM_REFLECT_EVEN
             local_bcval(side)=zero
            else if (presbc_in(dir,side,1).eq.FOEXTRAP) then
             local_bctype(side)=SEM_NEUMANN
             local_bcval(side)=zero
            else if (presbc_in(dir,side,1).eq.EXT_DIR) then
             local_bctype(side)=SEM_DIRICHLET
             local_bcval(side)=pres(D_DECL(i_out,j_out,k_out),1)
            else
             print *,"presbc_in corrupt OP_PRESGRAD_MAC or OP_POTGRAD_TO_MAC"
             stop
            endif
           else
            print *,"operation_flag invalid22"
            stop
           endif

           ! -1=fine(current) next to coarse
           ! -2=coarse(current) next to fine
          else if ((nbr_outside_domain_flag(side).eq.-1).or. &
                   (nbr_outside_domain_flag(side).eq.-2)) then

           local_bctype(side)=SEM_INTERIOR
           local_bcval(side)=zero
           ic=i_out
           jc=j_out
           kc=k_out
                
            ! fine(current) next to coarse 
           if (nbr_outside_domain_flag(side).eq.-1) then 
            bctype_tag=SEM_FINE_NEXT_TO_COARSE
            if (level.ge.1) then
             if (side.eq.1) then
              icoarse=(sideidx(1)-1)/2
              jcoarse=(sideidx(2)-1)/2
              kcoarse=(sideidx(SDIM)-1)/2
             else if (side.eq.2) then
              icoarse=sideidx(1)/2
              jcoarse=sideidx(2)/2
              kcoarse=sideidx(SDIM)/2
             else
              print *,"side invalid"
              stop
             endif
            
             call gridsten(xsten_coarse,problo, &
              icoarse,jcoarse,kcoarse, &
              domlo,bfact_c,dx_c,nhalf)

             if (side.eq.1) then
              x_sep(side)=two*(xsten_face(0,dir)-xsten_coarse(0,dir))/ &
               (dx(dir)*bfact)
              if ((x_sep(side).gt.zero).and. &
                  (x_sep(side).lt.two)) then
               ! do nothing
              else
               print *,"x_sep(side) invalid"
               stop
              endif
             else if (side.eq.2) then
              x_sep(side)=two*(xsten_coarse(0,dir)-xsten_face(0,dir))/ &
               (dx(dir)*bfact)
              if ((x_sep(side).gt.zero).and. &
                  (x_sep(side).lt.two)) then
               ! do nothing
              else
               print *,"x_sep(side) invalid"
               stop
              endif
             else
              print *,"side invalid"
              stop
             endif
            else
             print *,"level invalid 44"
             print *,"level=",level
             print *,"side=",side
             print *,"operation_flag=",operation_flag
             print *,"finest_level=",finest_level
             print *,"dir=",dir
             stop
            endif
             
             ! coarse(current) next to fine
           else if (nbr_outside_domain_flag(side).eq.-2) then
            bctype_tag=SEM_COARSE_NEXT_TO_FINE
            if (level.lt.finest_level) then
             if (side.eq.1) then
              ifine=2*sideidx(1)+1
              jfine=2*sideidx(2)+1
              kfine=2*sideidx(SDIM)+1
             else if (side.eq.2) then
              ifine=2*sideidx(1)
              jfine=2*sideidx(2)
              kfine=2*sideidx(SDIM)
             else
              print *,"side invalid"
              stop
             endif
            
             call gridsten(xsten_fine,problo, &
              ifine,jfine,kfine, &
              domlo,bfact_f,dx_f,nhalf)

             if (side.eq.1) then
              x_sep(side)=two*(xsten_face(0,dir)-xsten_fine(0,dir))/ &
               (dx(dir)*bfact)
              if ((x_sep(side).gt.zero).and. &
                  (x_sep(side).lt.half)) then
               ! do nothing
              else
               print *,"x_sep(side) invalid"
               stop
              endif
             else if (side.eq.2) then
              x_sep(side)=two*(xsten_fine(0,dir)-xsten_face(0,dir))/ &
               (dx(dir)*bfact)
              if ((x_sep(side).gt.zero).and. &
                  (x_sep(side).lt.half)) then
               ! do nothing
              else
               print *,"x_sep(side) invalid"
               stop
              endif
             else
              print *,"side invalid"
              stop
             endif
            else
             print *,"level invalid 45"
             stop
            endif

           else
            print *,"nbr_outside_domain_flag invalid 2"
            stop
           endif

           if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

            if (simple_AMR_BC_flag.eq.0) then
             local_bctype(side)=bctype_tag
             templocal=xp(D_DECL(iface_out,jface_out,kface_out),SEM_T+1)

             if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then
              local_data_side(side)= &
               xp(D_DECL(iface_out,jface_out,kface_out),nc)
             else if (nc.eq.SEM_T+1) then ! temperature
              local_data_side(side)=templocal
             else
              print *,"nc invalid"
              stop
             endif
            else if (simple_AMR_BC_flag.eq.1) then

             templocal=den(D_DECL(ic,jc,kc),ibase+ENUM_TEMPERATUREVAR+1)

             if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then
              local_data_side(side)= &
               vel(D_DECL(ic,jc,kc),nc) 
             else if (nc.eq.SEM_T+1) then ! temperature
              local_data_side(side)=templocal
             else
              print *,"nc invalid"
              stop
             endif

            else
             print *,"simple_AMR_BC_flag invalid"
             stop
            endif

            if ((abs(local_data_side(side)).lt.1.0D+20).and. &
                (abs(templocal).lt.1.0D+20)) then
             ! do nothing
            else
             print *,"data overflow SEM nc=",nc
             print *,"local_data_side(side)=",local_data_side(side)
             print *,"templocal=",templocal
             print *,"side=",side
             print *,"simple_AMR_BC_flag=",simple_AMR_BC_flag
             print *,"local_bctype(side) ",local_bctype(side)
             print *,"iface_out,jface_out,kface_out ", &
                iface_out,jface_out,kface_out
             print *,"ic,jc,kc ",ic,jc,kc
             print *,"level ",level
             print *,"finest_level ",finest_level
             print *,"bfact,bfact_c,bfact_f ",bfact,bfact_c,bfact_f
             print *,"den(D_DECL(ic,jc,kc),ibase+ENUM_TEMPERATUREVAR+1) ", &
                den(D_DECL(ic,jc,kc),ibase+ENUM_TEMPERATUREVAR+1)
             print *,"ibase= ",ibase
             print *,"cen_maskSEM= ",cen_maskSEM
             stop
            endif

           else if (operation_flag.eq.OP_UGRAD_MAC) then ! tensor derivatives

            local_AMR_BC_flag=simple_AMR_BC_flag

            if (cen_outside_domain_flag.eq.0) then

             if (strip_outside_fab_flag.eq.0) then
              ! do nothing
             else if (strip_outside_fab_flag.eq.1) then
              local_AMR_BC_flag=1
             else
              print *,"strip_outside_fab_flag invalid"
              stop
             endif

            else 
             print *,"cen_outside_domain_flag invalid"
             stop
            endif

            if (local_AMR_BC_flag.eq.0) then
             local_bctype(side)=bctype_tag
             if ((nc.ge.1).and.(nc.le.SDIM)) then
              local_data_side(side)= &
               xp(D_DECL(iface_out,jface_out,kface_out),nc)
             else
              print *,"nc invalid OP_UGRAD_MAC"
              stop
             endif
            else if (local_AMR_BC_flag.eq.1) then
             if ((nc.ge.1).and.(nc.le.SDIM)) then
              local_data_side(side)=vel(D_DECL(ic,jc,kc),nc)
             else
              print *,"nc invalid OP_UGRAD_MAC"
              stop
             endif
            else
             print *,"local_AMR_BC_flag invalid"
             stop
            endif

           else if (operation_flag.eq.OP_PRESGRAD_MAC) then!MAC pres gradient

            if (simple_AMR_BC_flag.eq.0) then
             local_bctype(side)=bctype_tag
             if (nc.eq.1) then
              if (scomp.eq.1) then
               local_data_side(side)= &
                xp(D_DECL(iface_out,jface_out,kface_out),scomp)
              else
               print *,"scomp invalid OP_PRESGRAD_MAC"
               stop
              endif
             else
              print *,"nc invalid OP_PRESGRAD_MAC"
              stop
             endif
            else if (simple_AMR_BC_flag.eq.1) then
             if (nc.eq.1) then
              if (scomp.eq.1) then
               local_data_side(side)=pres(D_DECL(ic,jc,kc),1)
              else
               print *,"scomp invalid OP_PRESGRAD_MAC"
               stop
              endif
             else
              print *,"nc invalid OP_PRESGRAD_MAC"
              stop
             endif
            else
             print *,"simple_AMR_BC_flag invalid: ",simple_AMR_BC_flag
             stop
            endif

           else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 

            if (simple_AMR_BC_flag.eq.1) then

             if (nc.eq.1) then
              if (scomp.eq.1) then
               local_data_side(side)=pres(D_DECL(ic,jc,kc),1)
               local_data_side_den(side)=den(D_DECL(ic,jc,kc),1)
              else
               print *,"scomp invalid OP_POTGRAD_TO_MAC"
               stop
              endif
             else
              print *,"nc invalid OP_POTGRAD_TO_MAC"
              stop
             endif

            else
             print *,"simple_AMR_BC_flag invalid OP_POTGRAD_TO_MAC"
             stop
            endif

            if ((abs(local_data_side(side)).lt.1.0D+20).and. &
                (local_data_side_den(side).gt.zero).and. &
                (abs(local_data_side_den(side)).lt.1.0D+20)) then
             ! do nothing
            else
             print *,"local_data_side or local_data_side_den invalid"
             stop
            endif

           else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
                    (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC).or. &
                    (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC)) then 
            if (simple_AMR_BC_flag.eq.0) then
             local_bctype(side)=bctype_tag
             if (nc.eq.1) then
              if (scomp.eq.dir) then
               ! AMRSYNC_VEL_MF passed in as xp
               local_data_side(side)= &
                xp(D_DECL(iface_out,jface_out,kface_out),1)
              else
               print *,"scomp invalid OP_UNEW or OP_U_COMP or OP_UMAC"
               stop
              endif
             else
              print *,"nc invalid OP_UNEW or OP_U_COMP or OP_UMAC"
              stop
             endif
            else if (simple_AMR_BC_flag.eq.1) then
             if (nc.eq.1) then
              if (scomp.eq.dir) then
               local_data_side(side)=vel(D_DECL(ic,jc,kc),scomp)
              else
               print *,"scomp invalid OP_UNEW or OP_U_COMP or OP_UMAC"
               stop
              endif
             else
              print *,"nc invalid OP_UNEW or OP_U_COMP or OP_UMAC"
              stop
             endif
            else
             print *,"simple_AMR_BC_flag invalid"
             stop
            endif
            if (abs(local_data_side(side)).lt.1.0D+20) then
             ! do nothing
            else
             print *,"SEM data overflow"
             stop
            endif

           else
            print *,"operation_flag invalid21"
            stop
           endif

          else
           print *,"nbr_outside_domain_flag invalid 3"
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

          if (operation_flag.eq.OP_UGRAD_MAC) then ! face grad U

           if ((nc.ge.1).and.(nc.le.SDIM)) then
            local_data(isten+1)=vel(D_DECL(ic,jc,kc),nc)
           else
            print *,"nc invalid(OP_UGRAD_MAC)"
            stop
           endif

          else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
                   (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC).or. &
                   (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC)) then

           if (nc.eq.1) then
            if (scomp.eq.dir) then
              !if OP_U_SEM_CELL_MAC_TO_MAC:
              !vel=primary_velfab=DELTA_CELL_VEL_MF
             local_data(isten+1)=vel(D_DECL(ic,jc,kc),scomp)
            else
             print *,"scomp invalid(OP_UNEW, or OP_U_COMP, or OP_UMAC)"
             stop
            endif
           else
            print *,"nc invalid(OP_UNEW, or OP_U_COMP, or OP_UMAC)"
            stop
           endif
         
          else if (operation_flag.eq.OP_PRESGRAD_MAC) then!pres grad on MAC

           if (nc.eq.1) then
            local_data(isten+1)=pres(D_DECL(ic,jc,kc),1)
           else
            print *,"nc invalid OP_PRESGRAD_MAC: ",nc
            stop
           endif

          else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 

           if (nc.eq.1) then
            local_data(isten+1)=pres(D_DECL(ic,jc,kc),1)
            local_data_den(isten+1)=den(D_DECL(ic,jc,kc),1)
           else
            print *,"nc invalid OP_POTGRAD_TO_MAC"
            stop
           endif

           !advection(values inside elem)
          else if (operation_flag.eq.OP_ISCHEME_MAC) then

           templocal=den(D_DECL(ic,jc,kc),ibase+ENUM_TEMPERATUREVAR+1)

            ! u dot grad u = div(umac u)-u div umac
           if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then ! velocity

            local_data(isten+1)= &
              vel(D_DECL(ic,jc,kc),nc) 

           else if (nc.eq.SEM_T+1) then ! temperature

            local_data(isten+1)=templocal

           else
            print *,"nc invalid"
            stop
           endif

          else
           print *,"operation_flag invalid23"
           stop
          endif
         enddo ! isten=0..bfact-1

         do isten=0,bfact
          indexmid(dir)=indexlo(dir)+isten
          ic=indexmid(1)
          jc=indexmid(2)
          kc=indexmid(SDIM)

           ! dir=1..sdim
          call gridstenMAC(xsten,xlo, &
           ic,jc,kc, &
           fablo,bfact,dx,nhalf,dir-1)

          RRface(isten)=xsten(0,1)

          local_vel(isten)=xvel(D_DECL(ic,jc,kc),1)
         enddo ! isten=0..bfact

         if (spectral_loop.eq.0) then

          ! if operation_flag.eq.OP_ISCHEME_MAC (advection), then
          ! velocity/temperature flux 
          ! will be multiplied by umac (local_vel) in
          ! lineGRAD.
          ! u u_x + v u_y + w u_z = (u umac)_x + (u vmac)_y + (u wmac)_z -
          !                         u umac_x - u vmac_y - u wmac_z
          ! calling from SEM_CELL_TO_MAC
          ! SEM_IMAGE_BC_ALG.eq.0 => spectral accuracy
          ! SEM_IMAGE_BC_ALG.eq.1 => 2nd order reflection
          call lineGRAD( &
           levelrz, &
           dir, &
           nc, &
           RRface, &
           local_bctype, &
           local_bcval, &
           local_vel, &
           local_data,local_data_side,  &
           local_grad,local_interp, &
           bfact, &
           dx(dir), &
           x_sep, &
           operation_flag)

          if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 

           do side=1,2
            local_bcval_den(side)=one  ! will not be used since "extrap" bc.
            local_bctype_den(side)=local_bctype(side)
            if (local_bctype_den(side).ne.SEM_INTERIOR) then
             local_bctype_den(side)=SEM_EXTRAP
            endif
           enddo
           ! calling from SEM_CELL_TO_MAC
           ! SEM_IMAGE_BC_ALG.eq.0 => spectral accuracy
           ! SEM_IMAGE_BC_ALG.eq.1 => 2nd order reflection
           call lineGRAD( &
            levelrz, &
            dir, &
            nc, &
            RRface, &
            local_bctype_den, &
            local_bcval_den, &
            local_vel, &
            local_data_den,local_data_side_den, &
            local_grad_den,local_interp_den, &
            bfact, &
            dx(dir), &
            x_sep, &
            operation_flag)
          else if ((operation_flag.eq.OP_PRESGRAD_MAC).or. & !grad p_MAC
                   (operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & !u^{c->mac}
                   (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & !u^mac+du^{c->mac}
                   (operation_flag.eq.OP_UGRAD_MAC).or. & !rate of strain
                   (operation_flag.eq.OP_ISCHEME_MAC).or. & !advection
                   (operation_flag.eq.OP_UGRAD_COUPLING_MAC).or. & !coupling
                   (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then
           ! do nothing
          else
           print *,"operation_flag invalid24"
           stop
          endif

         else if (spectral_loop.eq.1) then
          ! do nothing
         else
          print *,"spectral_loop invalid"
          stop
         endif

         do isten=0,bfact

          ! prevent race condition if tiling
          ! if two "pure" elements are neighbors and separate tiles, then we
          ! update the left most flux of each element.  The right
          ! most flux on the right most tile can be updated too.
          shared_face=0

          indexlo(1)=i
          indexlo(2)=j
          if (SDIM.eq.3) then
           indexlo(SDIM)=k
          endif

          do dir2=1,SDIM
           indexmid(dir2)=indexlo(dir2)
          enddo ! dir2
          indexmid(dir)=indexlo(dir)+isten ! isten=0..bfact

          ic=indexmid(1)
          jc=indexmid(2)
          kc=indexmid(SDIM)

          do dir2=1,SDIM
           index_edge(dir2)=indexmid(dir2)
           index_opp(dir2)=indexmid(dir2)
          enddo

          side=0

          if (isten.eq.0) then ! left most GL node

           side=1
           index_edge(dir)=indexlo(dir) ! left most G node
           index_opp(dir)=index_edge(dir)-1 ! right most G node of left elem.

          else if (isten.eq.bfact) then  ! right most GL node

           side=2
           index_edge(dir)=indexlo(dir)+bfact-1 ! right most G node
           index_opp(dir)=index_edge(dir)+1 ! left most G node of right elem.

          else if ((isten.ge.1).and.(isten.lt.bfact)) then
           ! do nothing
          else
           print *,"isten invalid isten=",isten
           stop
          endif

          i_in=index_edge(1)
          j_in=index_edge(2)
          k_in=index_edge(SDIM)
   
          i_out=index_opp(1)
          j_out=index_opp(2)
          k_out=index_opp(SDIM)

          test_maskSEM=NINT(maskSEM(D_DECL(i_out,j_out,k_out)))
          local_maskcov=NINT(maskcov(D_DECL(i_out,j_out,k_out)))

          if (side.eq.2) then ! right side of the element

           if ((index_edge(dir).ge.fablo(dir)).and. & !can conflict w/rt nbr.
               (index_edge(dir).lt.fabhi(dir)).and. &
               (test_maskSEM.eq.cen_maskSEM).and. &
               (local_maskcov.eq.1)) then
            shared_face=1
           else if ((index_edge(dir).eq.fabhi(dir)).or. &
                    (test_maskSEM.ne.cen_maskSEM).or. & !neighbor elem low ord?
                    (local_maskcov.eq.0)) then ! neighbor elem covered?
            ! do nothing
           else
            print *,"index_edge,test_maskSEM, or local_maskcov invalid"
            stop
           endif

          else if (side.eq.1) then
           ! do nothing (left side of element)
          else if (side.eq.0) then
           ! do nothing
          else
           print *,"side invalid"
           stop
          endif

          mask_out=1

          if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection
           do dir2=1,SDIM
            if (dir2.ne.dir) then
             if ((index_opp(dir2).lt.fablo(dir2)).or. &
                 (index_opp(dir2).gt.fabhi(dir2))) then
              print *,"index_opp invalid"
              print *,"dir2=",dir2
              print *,"index_opp(dir2)=",index_opp(dir2)
              print *,"fablo(dir2)=",fablo(dir2)
              print *,"fabhi(dir2)=",fabhi(dir2)
              stop
             endif
            else if (dir2.eq.dir) then
             if ((index_opp(dir2).lt.fablo(dir2)-1).or. &
                 (index_opp(dir2).gt.fabhi(dir2)+1)) then
              print *,"index_opp invalid"
              stop
             endif
            else
             print *,"dir2 invalid"
             stop
            endif
           enddo ! dir2=1..sdim
           if (index_opp(dir).eq.fablo(dir)-1) then
            if (side.eq.1) then
             if (nbr_outside_domain_flag(side).eq.1) then
              mask_out=0
             else if (nbr_outside_domain_flag(side).eq.0) then
              ! do nothing (fine next to fine)
             else if (nbr_outside_domain_flag(side).eq.-1) then
              ! do nothing (fine(current) next to coarse)
             else if (nbr_outside_domain_flag(side).eq.-2) then
              ! do nothing (coarse(current) next to fine)
              ! (this case will be masked off below)
             else
              print *,"nbr_outside_domain_flag invalid4"
              print *,"nbr_outside_domain_flag(side)=", &
                nbr_outside_domain_flag(side)
              print *,"dir,side ",dir,side
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           else if (index_opp(dir).eq.fabhi(dir)+1) then
            if (side.eq.2) then
             if (nbr_outside_domain_flag(side).eq.1) then
              mask_out=0
             else if (nbr_outside_domain_flag(side).eq.0) then
              ! do nothing (fine next to fine)
             else if (nbr_outside_domain_flag(side).eq.-1) then
              ! do nothing (fine(current) next to coarse)
             else if (nbr_outside_domain_flag(side).eq.-2) then
              ! do nothing (coarse(current) next to fine)
              ! (this case will be masked off below)
             else
              print *,"nbr_outside_domain_flag invalid5"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           else if ((index_opp(dir).ge.fablo(dir)).and. &
                    (index_opp(dir).le.fabhi(dir))) then
            ! do nothing
           else
            print *,"index_opp invalid"
            stop
           endif

           ! except for advection:
           ! do not average with flux from neighboring coarse grid or
           ! from outside the domain.
           ! maskCF==0 at coarse/fine   maskCF==1 at fine/fine

          else if ((operation_flag.eq.OP_PRESGRAD_MAC).or. & !grad p_MAC
                   (operation_flag.eq.OP_POTGRAD_TO_MAC).or. & 
                   (operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & !u^{c->mac}
                   (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                   (operation_flag.eq.OP_UGRAD_MAC).or. & !rate of strain
                   (operation_flag.eq.OP_UGRAD_COUPLING_MAC).or. & !coupling
                   (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then

           do dir2=1,SDIM
            if ((index_opp(dir2).lt.fablo(dir2)).or. &
                (index_opp(dir2).gt.fabhi(dir2))) then
             mask_out=NINT(maskCF(D_DECL(i_out,j_out,k_out)))
            endif
           enddo ! dir2

          else
           print *,"operation_flag invalid"
           stop
          endif

           ! do not average with flux from 
           ! an element that is covered. (avgDownEdge will handle this case)
          if (local_maskcov.eq.1) then
           ! do nothing
          else if (local_maskcov.eq.0) then
           mask_out=0
          else
           print *,"local_maskcov invalid"
           stop
          endif

           ! do not average with flux from 
           ! a low order  element.
          if (test_maskSEM.ne.cen_maskSEM) then
           mask_out=0
          endif

           ! shared_face=1 for faces on right side of elements and not
           !  touching the right side of the grid, not touching
           !  a maskSEM==0 element, and not touching a covered element.
           ! shared_face=0 for faces on the left side of elements and for
           !  the face touching the right side of the grid or touching
           !  a maskSEM!=cen_maskSEM element, or touching a covered element.
          if (shared_face.eq.1) then
           mask_out=0 
          else if (shared_face.eq.0) then
           ! do nothing
          else
           print *,"shared_face invalid"
           stop
          endif

          if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

           if (spectral_loop.eq.0) then

            if (shared_face.eq.0) then

             if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then ! u * velocity 

              xface(D_DECL(ic,jc,kc),nc)=local_interp(isten+1)

               ! temperature
             else if (nc.eq.SEM_T+1) then ! u * temperature

              xface(D_DECL(ic,jc,kc),nc)=local_interp(isten+1)

             else
              print *,"nc invalid"
              stop
             endif

            else if (shared_face.eq.1) then
             ! do nothing
            else
             print *,"shared_face invalid"
             stop
            endif

             ! semflux is a cell centered FAB with 1 ghost cell.
             ! fluxbase=(dir-1)*ngroup
            if ((side.eq.1).or.(side.eq.2)) then
             semflux(D_DECL(i_in,j_in,k_in),fluxbase+nc)=local_interp(isten+1)
            else if (side.eq.0) then
             ! do nothing
            else
             print *,"side invalid"
             stop
            endif

           else if (spectral_loop.eq.1) then

             ! 1. neighbor element is not a low order element
             ! 2. shared_face==0
            if (mask_out.eq.1) then

             if (side.eq.0) then
              ! do nothing
             else if ((side.eq.1).or.(side.eq.2)) then
               ! fine(current) next to coarse
              if (nbr_outside_domain_flag(side).eq.-1) then
               outside_flux=xgp(D_DECL(ic,jc,kc),nc) ! interpolated flux.
              else if (nbr_outside_domain_flag(side).eq.0) then
               outside_flux=semflux(D_DECL(i_out,j_out,k_out),fluxbase+nc)
              else
               print *,"nbr_outside_domain_flag invalid6"
               stop
              endif
              inside_flux=semflux(D_DECL(i_in,j_in,k_in),fluxbase+nc)

              if ((abs(inside_flux).lt.1.0D+20).and. &
                  (abs(outside_flux).lt.1.0D+20)) then

                if (((side.eq.1).and.(local_vel(0).ge.zero)).or. &
                    ((side.eq.2).and.(local_vel(bfact).le.zero))) then
                 xface(D_DECL(ic,jc,kc),nc)=outside_flux
                else if (((side.eq.1).and.(local_vel(0).le.zero)).or. &
                         ((side.eq.2).and.(local_vel(bfact).ge.zero))) then
                 xface(D_DECL(ic,jc,kc),nc)=inside_flux
                else
                 print *,"side invalid"
                 stop
                endif

              else
               print *,"inside_flux or outside_flux overflow"
               stop
              endif

             else
              print *,"side invalid"
              stop
             endif

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

           ! tensor derivatives
          else if (operation_flag.eq.OP_UGRAD_MAC) then 

           if (spectral_loop.eq.0) then

            if (shared_face.eq.0) then
             xgp(D_DECL(ic,jc,kc),dcomp+nc-1)=local_grad(isten+1)
            else if (shared_face.eq.1) then
             ! do nothing
            else
             print *,"shared_face invalid"
             stop
            endif

            if ((side.eq.1).or.(side.eq.2)) then
             semflux(D_DECL(i_in,j_in,k_in),nbase+nc)=local_grad(isten+1)
            else if (side.eq.0) then
             ! do nothing
            else
             print *,"side invalid"
             stop
            endif

           else if (spectral_loop.eq.1) then

            if (mask_out.eq.1) then

             if ((side.eq.1).or.(side.eq.2)) then
              xgp(D_DECL(ic,jc,kc),dcomp+nc-1)=half*( &
               semflux(D_DECL(i_in,j_in,k_in),nbase+nc)+ &
               semflux(D_DECL(i_out,j_out,k_out),nbase+nc))
             else if (side.eq.0) then
              ! do nothing
             else
              print *,"side invalid"
              stop
             endif

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

          else if (operation_flag.eq.OP_PRESGRAD_MAC) then 

           if (ncfluxreg.ne.SDIM) then
            print *,"ncfluxreg invalid5 ",ncfluxreg
            stop
           endif

           shared_xcut=xcut(D_DECL(ic,jc,kc),1)

           if (side.eq.0) then
            ! do nothing
           else if ((side.eq.1).or.(side.eq.2)) then
            if (nbr_outside_domain_flag(side).eq.1) then
             if ((local_bctype(side).eq.SEM_REFLECT_EVEN).or. & 
                 (local_bctype(side).eq.SEM_NEUMANN)) then ! neumann
              if (shared_xcut.eq.zero) then
               ! do nothing
              else
               print *,"shared_xcut invalid: ",shared_xcut
               stop
              endif
             else if ((local_bctype(side).eq.SEM_DIRICHLET).or. &
                      (local_bctype(side).eq.SEM_INTERIOR)) then
              ! do nothing
             else
              print *,"local_bctype invalid"
              stop
             endif
            else if (nbr_outside_domain_flag(side).eq.0) then
             ! do nothing
            else if (nbr_outside_domain_flag(side).eq.-2) then
             ! do nothing (coarse(current) next to fine)
            else if (nbr_outside_domain_flag(side).eq.-1) then 
             ! do nothing (fine(current) next to coarse)
            else
             print *,"nbr_outside_domain_flag invalid 4"
             stop
            endif 
           else
            print *,"side invalid"
            stop
           endif

           if (spectral_loop.eq.0) then
         
            ! regular solver or SDC viscosity or thermal flux.
            if (energyflag.eq.SUB_OP_FOR_MAIN) then 
             shared_face_value=-dt*local_grad(isten+1)*shared_xcut

            ! for SDC pressure gradient
            else if (energyflag.eq.SUB_OP_FOR_SDC) then 
             shared_face_value=local_grad(isten+1)

            else
             print *,"energyflag invalid OP_PRESGRAD_MAC"
             stop
            endif

            if (shared_face.eq.0) then
             xgp(D_DECL(ic,jc,kc),dcomp)=shared_face_value
            else if (shared_face.eq.1) then
             ! do nothing
            else
             print *,"shared_face invalid"
             stop
            endif

            if ((side.eq.1).or. &
                (side.eq.2)) then
             semflux(D_DECL(i_in,j_in,k_in),fluxbase+nc)=shared_face_value
            else if (side.eq.0) then
             ! do nothing
            else
             print *,"side invalid"
             stop
            endif

           else if (spectral_loop.eq.1) then

            if (mask_out.eq.1) then
             if ((side.eq.1).or.(side.eq.2)) then
              xgp(D_DECL(ic,jc,kc),dcomp)=half*( &
               semflux(D_DECL(i_in,j_in,k_in),fluxbase+nc)+ &
               semflux(D_DECL(i_out,j_out,k_out),fluxbase+nc))
             else if (side.eq.0) then
              ! do nothing
             else
              print *,"side invalid"
              stop
             endif
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

          else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 

           ! potential pressure gradient/den 
           if (ncfluxreg.ne.SDIM) then
            print *,"ncfluxreg invalid7 OP_POTGRAD_TO_MAC ",ncfluxreg
            stop
           endif

           if (spectral_loop.eq.0) then

            if (local_interp_den(isten+1).gt.zero) then
             ! do nothing
            else
             print *,"local_interp_den underflow"
             stop
            endif

            shared_face_value=local_grad(isten+1)/local_interp_den(isten+1)

            if (shared_face.eq.0) then
             if (nc.eq.1) then
              ! grad ppot/rho_pot at face.
              xgp(D_DECL(ic,jc,kc),nc)=shared_face_value
             else
              print *,"nc invalid OP_POTGRAD_TO_MAC"
              stop
             endif
            else if (shared_face.eq.1) then
             ! do nothing
            else
             print *,"shared_face invalid OP_POTGRAD_TO_MAC"
             stop
            endif

            if ((side.eq.1).or.(side.eq.2)) then
              ! grad ppot/rho_pot at face.
             semflux(D_DECL(i_in,j_in,k_in),fluxbase+1)=shared_face_value
            else if (side.eq.0) then
             ! do nothing
            else
             print *,"side invalid"
             stop
            endif
           else if (spectral_loop.eq.1) then
            if (mask_out.eq.1) then
             if ((side.eq.1).or.(side.eq.2)) then
              ! grad ppot/rho_pot at face.
              xgp(D_DECL(ic,jc,kc),nc)=half*( &
               semflux(D_DECL(i_in,j_in,k_in),fluxbase+1)+ &
               semflux(D_DECL(i_out,j_out,k_out),fluxbase+1))
             else if (side.eq.0) then
              ! do nothing
             else
              print *,"side invalid"
              stop
             endif
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

          else if (operation_flag.eq.OP_UNEW_CELL_TO_MAC) then   

           if (ncfluxreg.ne.SDIM) then
            print *,"ncfluxreg invalid8 ",ncfluxreg
            stop
           endif

           if (spectral_loop.eq.0) then

            if (shared_face.eq.0) then
             xvel(D_DECL(ic,jc,kc),1)=local_interp(isten+1)
            else if (shared_face.eq.1) then
             ! do nothing
            else
             print *,"shared_face invalid"
             stop
            endif

            if ((side.eq.1).or.(side.eq.2)) then
             semflux(D_DECL(i_in,j_in,k_in),fluxbase+1)=local_interp(isten+1)
            else if (side.eq.0) then
             ! do nothing
            else
             print *,"side invalid"
             stop
            endif

           else if (spectral_loop.eq.1) then

            if (mask_out.eq.1) then
             if ((side.eq.1).or.(side.eq.2)) then
              xvel(D_DECL(ic,jc,kc),1)=half*( &
               semflux(D_DECL(i_in,j_in,k_in),fluxbase+1)+ &
               semflux(D_DECL(i_out,j_out,k_out),fluxbase+1))
             else if (side.eq.0) then
              ! do nothing
             else
              print *,"side invalid"
              stop
             endif
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

          else if (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC) then 

           if (ncfluxreg.ne.SDIM) then
            print *,"ncfluxreg invalid8 ",ncfluxreg
            stop
           endif
           if (spectral_loop.eq.0) then

             !xgp=Umac_old if OP_U_SEM_CELL_MAC_TO_MAC
            shared_face_value=xgp(D_DECL(ic,jc,kc),dcomp)+ &
              local_interp(isten+1)

            if (shared_face.eq.0) then
             xvel(D_DECL(ic,jc,kc),1)=shared_face_value
            else if (shared_face.eq.1) then
             ! do nothing
            else
             print *,"shared_face invalid"
             stop
            endif

            if ((side.eq.1).or.(side.eq.2)) then
             semflux(D_DECL(i_in,j_in,k_in),fluxbase+1)=shared_face_value
            else if (side.eq.0) then
             ! do nothing
            else
             print *,"side invalid"
             stop
            endif

           else if (spectral_loop.eq.1) then

            if (mask_out.eq.1) then
             if ((side.eq.1).or.(side.eq.2)) then
              xvel(D_DECL(ic,jc,kc),1)=half*( &
               semflux(D_DECL(i_in,j_in,k_in),fluxbase+1)+ &
               semflux(D_DECL(i_out,j_out,k_out),fluxbase+1))
             else if (side.eq.0) then
              ! do nothing
             else
              print *,"side invalid"
              stop
             endif
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

          else if (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC) then 

           if (ncfluxreg.ne.SDIM) then
            print *,"ncfluxreg invalid9 ",ncfluxreg
            stop
           endif
           if (spectral_loop.eq.0) then

            shared_face_value=xgp(D_DECL(ic,jc,kc),dcomp)+ &
              beta*local_interp(isten+1)

            if (shared_face.eq.0) then
             xvel(D_DECL(ic,jc,kc),1)=shared_face_value
            else if (shared_face.eq.1) then
             ! do nothing
            else
             print *,"shared_face invalid"
             stop
            endif

            if ((side.eq.1).or.(side.eq.2)) then
             semflux(D_DECL(i_in,j_in,k_in),fluxbase+1)=shared_face_value
            else if (side.eq.0) then
             ! do nothing
            else
             print *,"side invalid"
             stop
            endif
           else if (spectral_loop.eq.1) then
            if (mask_out.eq.1) then
             if ((side.eq.1).or.(side.eq.2)) then
              xvel(D_DECL(ic,jc,kc),1)=half*( &
               semflux(D_DECL(i_in,j_in,k_in),fluxbase+1)+ &
               semflux(D_DECL(i_out,j_out,k_out),fluxbase+1))
             else if (side.eq.0) then
              ! do nothing
             else
              print *,"side invalid"
              stop
             endif
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

          else
           print *,"operation_flag invalid25"
           stop
          endif

         enddo ! isten=0..bfact
 
        else
         print *,"cen_outside_domain_flag: ",cen_outside_domain_flag
         print *,"cen_maskSEM: ",cen_maskSEM
         stop
        endif

       else if (cen_maskSEM.eq.0) then
        ! do nothing
       else
        print *,"cen_maskSEM invalid SEM_CELL_TO_MAC: ",cen_maskSEM
        print *,"i,j,k : ",i,j,k
        stop
       endif

      enddo ! nc=1..ncomp_dest

      return
      end subroutine SEM_CELL_TO_MAC
     
       ! OP_RHS_CELL 
       ! operation_flag=100 (right hand side for solver)
       ! OP_DIV_CELL 
       ! operation_flag=110 (divergence)
       ! OP_VEL_MAC_TO_CELL 
       ! operation_flag=103 (mac -> cell velocity in solver)
       ! OP_GRADU_MAC_TO_CELL 
       ! operation_flag=105 interpolate grad u from MAC to CELL.
       ! OP_ISCHEME_CELL 
       ! operation_flag=107 advection
       ! mask>0 SEM; mask=0 piecewise finite volume
      subroutine SEM_MAC_TO_CELL( &
       ncomp_denold, &
       ncomp_veldest, &
       ncomp_dendest, &
       ns_time_order, &
       divu_outer_sweeps, &
       num_divu_outer_sweeps, &
       SDC_outer_sweeps, &
       level, &
       finest_level, &
       operation_flag, &
       project_option, &
       energyflag, &
       homflag, &
       maskSEM, &
       time, &
       slab_step, &
       dt, &
       i,j,k, &
       tilelo,tilehi, &
       fablo,fabhi, &
       xlo,dx, &
       dir_main, &
       bfact, &
       velbc_in, &
       presbc_in, &
       scomp, &
       scomp_bc, &
       dcomp, &
       ncomp, &
       ncomp_xvel, &
       ncomp_cterm, &
       vol, & 
       xface, & 
       xvel, & 
       maskcoef, & 
       cterm, & 
       mdotcell, &  ! VELADVECT_MF, OP_ISCHEME_CELL
       maskdivres, & ! DEN_RECON_MF, OP_ISCHEME_CELL
       pold, & 
       denold, & 
       vel_old_fab, & 
       veldest, &
       dendest, &
       divdest)
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: ncomp_denold
      integer, INTENT(in) :: ncomp_veldest
      integer, INTENT(in) :: ncomp_dendest
      integer, INTENT(in) :: ns_time_order
      integer, INTENT(in) :: divu_outer_sweeps
      integer, INTENT(in) :: num_divu_outer_sweeps
      integer, INTENT(in) :: SDC_outer_sweeps
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: slab_step
      integer, INTENT(in) :: operation_flag
      integer, INTENT(in) :: homflag
      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: energyflag
      integer :: advect_iter
      integer :: source_term
      integer, INTENT(in) :: maskSEM
      real(amrex_real), INTENT(in) :: time,dt
      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: dir_main
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: scomp
      integer, INTENT(in) :: scomp_bc
      integer, INTENT(in) :: dcomp
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: ncomp_xvel
      integer, INTENT(in) :: ncomp_cterm
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: velbc_in(SDIM,2,SDIM)
      integer, INTENT(in) :: presbc_in(SDIM,2)
      real(amrex_real), INTENT(in), pointer :: vol(D_DECL(:,:,:),:)
       ! flux data for I-scheme
      real(amrex_real), INTENT(in), pointer :: xface(D_DECL(:,:,:),:)  
      real(amrex_real), INTENT(in), pointer :: xvel(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: maskcoef(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: cterm(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: &
          mdotcell(D_DECL(:,:,:),:) !VELADVECT_MF,if OP_ISCHEME_CELL
      real(amrex_real), INTENT(in), pointer :: &
          maskdivres(D_DECL(:,:,:),:) !DEN_RECON_MF,if OP_ISCHEME_CELL
      real(amrex_real), INTENT(in), pointer :: pold(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: denold(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: dendest(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: vel_old_fab(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: veldest(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), pointer :: divdest(D_DECL(:,:,:),:)
      real(amrex_real) local_data(bfact+1)
      real(amrex_real) local_vel_data(bfact+1)
      real(amrex_real) local_vel_data_div(bfact+1)
      real(amrex_real) local_cell(bfact)
      real(amrex_real) local_vel_cell(bfact)
      real(amrex_real) local_vel_cell_div(bfact)
      real(amrex_real) local_vel_div(bfact)
      real(amrex_real) local_vel_div_div(bfact)
      real(amrex_real) local_div(bfact)
      integer ii,jj,kk
      integer ic,jc,kc
      integer i_out,j_out,k_out
      integer nc,side
      integer nbr_outside_domain_flag(2)
      integer cen_outside_domain_flag
      integer isten
      integer dir2
      integer sideidx(SDIM)
      integer indexlo(SDIM)
      integer indexhi(SDIM)
      integer indexmid(SDIM)
      integer, parameter :: nhalf=3
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) RR,RRTHETA,RR_DIVIDE
      integer ibase,nc2
      real(amrex_real) divu,CC,CC_DUAL,MDOT,RHS
      real(amrex_real) local_POLD
      real(amrex_real) local_POLD_DUAL
      real(amrex_real) divflux(ncomp)
      real(amrex_real) vel_old,vel_update,T_old,T_new
      real(amrex_real) vel_new(SDIM)
      real(amrex_real) VOLTERM
      integer nbase
      real(amrex_real) local_div_val

      if ((dir_main.ge.1).and.(dir_main.le.SDIM)) then
       ! do nothing
      else
       print *,"dir_main invalid sem mac to cell"
       stop
      endif

      if ((slab_step.lt.-1).or.(slab_step.gt.bfact_time_order)) then
       print *,"slab_step invalid sem mac to cell "
       stop
      endif
      if ((level.lt.0).or. &
          (level.gt.finest_level)) then
       print *,"level or finest_level invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid200 mac to cell"
       stop
      endif

      if ((ns_time_order.ge.1).and.(ns_time_order.le.32)) then
       ! do nothing
      else
       print *,"ns_time_order invalid"
       stop
      endif
      if ((SDC_outer_sweeps.ge.0).and. &
          (SDC_outer_sweeps.lt.ns_time_order)) then
       ! do nothing
      else
       print *,"SDC_outer_sweeps invalid"
       stop
      endif
      if (num_divu_outer_sweeps.lt.1) then
       print *,"num_divu_outer_sweeps invalid SEM_MAC_TO_CELL"
       stop
      endif
      if ((divu_outer_sweeps.ge.0).and. &
          (divu_outer_sweeps.lt.num_divu_outer_sweeps)) then
       ! do nothing
      else
       print *,"divu_outer_sweeps invalid SEM_MAC_TO_CELL"
       stop
      endif

      advect_iter=energyflag
      source_term=homflag

      if ((maskSEM.ge.1).and. &
          (maskSEM.le.num_materials)) then

       ! do nothing

      else if (maskSEM.eq.0) then
       ! do nothing
      else
       print *,"maskSEM invalid sem_mac_to_cell: ",maskSEM
       print *,"i,j,k : ",i,j,k
       stop
      endif

      if (operation_flag.eq.OP_GRADU_MAC_TO_CELL) then ! grad U: MAC -> CELL

       if ((maskSEM.lt.1).or.(maskSEM.gt.num_materials)) then
        print *,"maskSEM invalid"
        stop
       endif 
       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_GRADU_MAC_TO_CELL"
        stop
       endif
       if (project_option.ne.SOLVETYPE_VISC) then
        print *,"project_option.ne.SOLVETYPE_VISC; SEM_MAC_TO_CELL"
        stop
       endif

       ! u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
       nbase=dcomp-1

       if ((nbase.ne.0).and.(nbase.ne.SDIM).and. &
           (nbase.ne.(SDIM-1)*SDIM)) then
        print *,"nbase invalid"
        stop
       endif

       if (scomp.ne.dcomp) then
        print *,"scomp invalid"
        stop
       endif
       if (ncomp.ne.SDIM) then
        print *,"ncomp invalid1"
        stop
       endif
       if ((ncomp_xvel.eq.AMREX_SPACEDIM_SQR).and. &
           (ncomp_veldest.eq.AMREX_SPACEDIM_SQR).and. &
           (ncomp_dendest.eq.AMREX_SPACEDIM_SQR).and. &
           (ncomp_denold.eq.AMREX_SPACEDIM_SQR).and. &
           (ncomp_cterm.eq.AMREX_SPACEDIM_SQR)) then
        ! do nothing
       else
        print *,"nc_xvel,nc_veldest,nc_dendest,nc_denold or nc_cterm invalid"
        stop
       endif

      else if (operation_flag.eq.OP_RHS_CELL) then ! RHS for solver

        ! SOLVETYPE_PRES, PRESGRAVITY, INITPROJ
       if (project_option_projectionF(project_option).eq.1) then
        if (ncomp.ne.1) then
         print *,"ncomp invalid2 SEM_MAC_TO_CELL"
         stop
        endif
        if (ncomp_xvel.ne.1) then
         print *,"ncomp_xvel invalid2 project_option_projectionF"
         stop
        endif
       else if (project_option.eq.SOLVETYPE_PRESEXTRAP) then 
        print *,"extension project should be low order"
        stop
       else if (project_option.eq.SOLVETYPE_VISC) then ! viscosity
        if (ncomp.ne.SDIM) then
         print *,"ncomp invalid3"
         stop
        endif
        if (ncomp_xvel.ne.SDIM) then
         print *,"ncomp_xvel invalid3 SOLVETYPE_VISC"
         stop
        endif
       else if (project_option.eq.SOLVETYPE_HEAT) then ! thermal conduction
        if (ncomp.ne.1) then
         print *,"ncomp invalid4"
         stop
        endif
        if (ncomp_xvel.ne.1) then
         print *,"ncomp_xvel invalid4 SOLVETYPE_HEAT"
         stop
        endif
       else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
        if (ncomp.ne.1) then
         print *,"ncomp invalid5"
         stop
        endif
        if (ncomp_xvel.ne.1) then
         print *,"ncomp_xvel invalid5 SOLVETYPE_SPEC"
         stop
        endif
       else
        print *,"project_option invalid SEM_MAC_TO_CELL"
        stop
       endif
       if ((maskSEM.lt.1).or. &
           (maskSEM.gt.num_materials)) then!operation_flag=OP_RHS_CELL
        print *,"maskSEM invalid SEM_MAC_TO_CELL"
        stop
       endif 
       if ((scomp.ne.1).or. &
           (dcomp.ne.1)) then
        print *,"scomp or dcomp invalid"
        stop
       endif
       if ((ncomp.ne.1).and. &
           (ncomp.ne.SDIM)) then
        print *,"ncomp invalid6"
        stop
       endif
       if ((scomp_bc.lt.1).or.(scomp_bc.gt.SDIM)) then
        print *,"scomp_bc invalid"
        stop
       endif
       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_RHS_CELL"
        stop
       endif
       if ((ncomp_xvel.ne.1).and. &
           (ncomp_xvel.ne.SDIM)) then
        print *,"ncomp_xvel invalid OP_RHS_CELL"
        stop
       endif
       if (ncomp_veldest.eq.ncomp) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.eq.ncomp) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncomp_denold.eq.ncomp) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

       if (ncomp_cterm.ne.ncomp) then
        print *,"ncomp_cterm invalid"
        stop
       endif

      else if (operation_flag.eq.OP_DIV_CELL) then ! divergence

       if ((maskSEM.lt.1).or.(maskSEM.gt.num_materials)) then
        print *,"maskSEM invalid"
        stop
       endif 
       if ((scomp_bc.lt.1).or.(scomp_bc.gt.SDIM)) then
        print *,"scomp_bc invalid"
        stop
       endif
       if ((scomp.ne.1).or. &
           (dcomp.ne.1).or. &
           (ncomp.ne.1)) then
        print *,"scomp, dcomp, or ncomp invalid"
        stop
       endif
       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_DIV_CELL"
        stop
       endif
       if (ncomp_xvel.ne.1) then
        print *,"ncomp_xvel invalid OP_DIV_CELL"
        stop
       endif

       if (ncomp_veldest.eq.1) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.eq.1) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncomp_denold.eq.1) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif
       if (ncomp_cterm.ne.1) then
        print *,"ncomp_cterm invalid"
        stop
       endif

      else if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then

       if ((maskSEM.lt.1).or.(maskSEM.gt.num_materials)) then
        print *,"maskSEM invalid"
        stop
       endif 
       if ((scomp_bc.lt.1).or.(scomp_bc.gt.SDIM)) then
        print *,"scomp_bc invalid"
        stop
       endif
       if ((dcomp.ge.1).and.(dcomp.le.SDIM)) then
        ! do nothing
       else
        print *,"dcomp invalid"
        stop
       endif
       if (dcomp.ne.dir_main) then
        print *,"dcomp invalid"
        stop
       endif
       if ((scomp.ne.1).or. &
           (ncomp.ne.1)) then
        print *,"scomp or ncomp invalid"
        stop
       endif
       if ((energyflag.ne.SUB_OP_THERMAL_DIVUP_NULL).and. &
           (energyflag.ne.SUB_OP_THERMAL_DIVUP_OK)) then
        print *,"energyflag invalid OP_VEL_MAC_TO_CELL"
        stop
       endif
       if (ncomp_xvel.ne.1) then
        print *,"ncomp_xvel invalid OP_VEL_MAC_TO_CELL"
        stop
       endif

       if (ncomp_veldest.ge. &
           SDIM+num_state_material*num_materials) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.ge.num_state_material*num_materials) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if ((ncomp_denold.eq.ncomp).or. &
           (ncomp_denold.eq.1)) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

       if (ncomp_cterm.ne.1) then
        print *,"ncomp_cterm invalid"
        stop
       endif

      else if (operation_flag.eq.OP_ISCHEME_CELL) then ! advection

       if ((maskSEM.ge.1).and.(maskSEM.le.num_materials)) then
        ! do nothing
       else
        print *,"maskSEM invalid"
        stop
       endif 

       if (ncomp_veldest.eq.STATE_NCOMP) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.eq.ncomp_veldest-STATECOMP_STATES) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncomp_denold.eq.num_materials*num_state_material) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

       if (ncomp.eq.NFLUXSEM) then
        ! do nothing
       else
        print *,"ncomp invalid7"
        stop
       endif
       if (source_term.eq.SUB_OP_SDC_LOW_TIME) then 
        if (advect_iter.ne.SUB_OP_ISCHEME_PREDICT) then 
         print *,"advect_iter invalid"
         stop
        endif
       else if (source_term.eq.SUB_OP_SDC_ISCHEME) then 
        if ((advect_iter.ne.SUB_OP_ISCHEME_PREDICT).and. &
            (advect_iter.ne.SUB_OP_ISCHEME_CORRECT)) then
         print *,"advect_iter invalid"
         stop
        endif
       else
        print *,"source_term invalid"
        stop
       endif
       if ((scomp.ne.1).or. &
           (dcomp.ne.1).or. &
           (scomp_bc.ne.1)) then
        print *,"scomp, dcomp, or scomp_bc invalid"
        stop
       endif
       if (ncomp_xvel.ne.NFLUXSEM) then
        print *,"ncomp_xvel invalid OP_ISCHEME_CELL"
        stop
       endif
       if (ncomp_cterm.ne.ncomp) then
        print *,"ncomp_cterm invalid"
        stop
       endif

      else
       print *,"operation_flag invalid26:",operation_flag
       stop
      endif

      ii=0
      jj=0
      kk=0

      if (dir_main.eq.1) then
       ii=1
       if ((i/bfact)*bfact.ne.i) then
        print *,"i invalid"
        stop
       endif
       if (i.lt.0) then
        print *,"i invalid"
        stop
       endif
      else if (dir_main.eq.2) then
       jj=1
       if ((j/bfact)*bfact.ne.j) then
        print *,"j invalid"
        stop
       endif
       if (j.lt.0) then
        print *,"j invalid"
        stop
       endif
      else if ((dir_main.eq.3).and.(SDIM.eq.3)) then
       kk=1
       if ((k/bfact)*bfact.ne.k) then
        print *,"k invalid SEM_MAC_TO_CELL"
        stop
       endif
       if (k.lt.0) then
        print *,"k invalid SEM_MAC_TO_CELL"
        stop
       endif
      else
       print *,"dir_main invalid sem mac to cell2"
       stop
      endif

      do dir2=1,SDIM
       if (fablo(dir2).ge.0) then
        ! do nothing
       else
        print *,"fablo invalid"
        stop
       endif
      enddo ! dir2=1..sdim

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid dt=",dt
       print *,"operation_flag= ",operation_flag
       stop
      endif

      if ((maskSEM.lt.1).or.(maskSEM.gt.num_materials)) then
       print *,"maskSEM invalid"
       stop
      endif 

      ibase=(maskSEM-1)*num_state_material

      if ((dir_main.ge.1).and. &
          (dir_main.le.SDIM).and. &
          (bfact.ge.1)) then
       ! do nothing
      else
       print *,"dir_main or bfact invalid200 in SEM_MAC_TO_CELL"
       stop
      endif

      indexlo(1)=i
      indexlo(2)=j
      if (SDIM.eq.3) then
       indexlo(SDIM)=k
      endif
 
      do dir2=1,SDIM
       indexhi(dir2)=indexlo(dir2)
      enddo
      indexhi(dir_main)=indexlo(dir_main)+bfact-1

      cen_outside_domain_flag=0

      do dir2=1,SDIM
       if (indexlo(dir2).lt.fablo(dir2)) then
        if (dir2.eq.dir_main) then
         print *,"indexlo(dir_main) invalid"
         stop
        endif
        if (velbc_in(dir2,1,dir2).ne.INT_DIR) then
         cen_outside_domain_flag=1
        endif
       endif 
       if (indexhi(dir2).gt.fabhi(dir2)) then
        if (dir2.eq.dir_main) then
         print *,"indexhi(dir_main) invalid"
         stop
        endif
        if (velbc_in(dir2,2,dir2).ne.INT_DIR) then
         cen_outside_domain_flag=1
        endif
       endif
      enddo ! dir2=1..sdim

      do nc=1,ncomp

       if (cen_outside_domain_flag.eq.1) then
        ! do nothing
       else if (cen_outside_domain_flag.eq.0) then

        do isten=0,bfact

         do dir2=1,SDIM
          indexmid(dir2)=indexlo(dir2)
         enddo
         indexmid(dir_main)=indexlo(dir_main)+isten
 
         ic=indexmid(1)
         jc=indexmid(2)
         kc=indexmid(SDIM)

         if (operation_flag.eq.OP_GRADU_MAC_TO_CELL) then 
          local_data(isten+1)=xvel(D_DECL(ic,jc,kc),scomp+nc-1)
         else if (operation_flag.eq.OP_ISCHEME_CELL) then ! advection
           ! u dot grad S=div(uS)-S div(u)
          if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then
           local_data(isten+1)=xface(D_DECL(ic,jc,kc),nc) ! u * umac 
          else if (nc.eq.SEM_T+1) then
           local_data(isten+1)=xface(D_DECL(ic,jc,kc),nc) ! temperature * umac
          else
           print *,"nc invalid"
           stop
          endif
          local_vel_data(isten+1)=xvel(D_DECL(ic,jc,kc),1) ! umac
          local_vel_data_div(isten+1)=xvel(D_DECL(ic,jc,kc),1) ! umac

         else if (operation_flag.eq.OP_RHS_CELL) then ! RHS

          local_data(isten+1)=xvel(D_DECL(ic,jc,kc),nc)

         else if (operation_flag.eq.OP_DIV_CELL) then ! divergence

          local_data(isten+1)=xvel(D_DECL(ic,jc,kc),1)

         else if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then

          local_data(isten+1)=xvel(D_DECL(ic,jc,kc),1)

         else
          print *,"operation_flag invalid27:",operation_flag
          stop
         endif

         RR=one

         if (dir_main.eq.1) then  ! r direction

           ! dir_main=1..sdim
          call gridstenMAC(xstenMAC,xlo, &
           ic,jc,kc, &
           fablo,bfact,dx,nhalf,dir_main-1)

          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if (levelrz.eq.COORDSYS_RZ) then
           RR=xstenMAC(0,1)
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           RR=xstenMAC(0,1)
          else
           print *,"levelrz invalid"
           stop
          endif

         else if ((dir_main.eq.2).or.(dir_main.eq.SDIM)) then
          ! do nothing
         else
          print *,"dir_main invalid sem mac to cell 3"
          stop
         endif 

         if (operation_flag.eq.OP_GRADU_MAC_TO_CELL) then ! interp grad U^T
          ! do nothing
         else if (operation_flag.eq.OP_RHS_CELL) then !RHS
          local_data(isten+1)=local_data(isten+1)*RR
         else if (operation_flag.eq.OP_ISCHEME_CELL) then !advection
          local_data(isten+1)=local_data(isten+1)*RR
          local_vel_data_div(isten+1)=local_vel_data_div(isten+1)*RR
         else if (operation_flag.eq.OP_DIV_CELL) then !divergence
          local_data(isten+1)=local_data(isten+1)*RR
         else if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then
          local_data(isten+1)=local_data(isten+1)*RR
         else
          print *,"operation_flag invalid28"
          stop
         endif

        enddo ! isten=0..bfact

        do side=1,2

         nbr_outside_domain_flag(side)=0

         do dir2=1,SDIM
          sideidx(dir2)=indexlo(dir2)
         enddo

         if (side.eq.1) then
          sideidx(dir_main)=indexlo(dir_main)-1
          do dir2=1,SDIM
           if (sideidx(dir2).lt.fablo(dir2)) then
            if (velbc_in(dir2,side,dir2).ne.INT_DIR) then
             nbr_outside_domain_flag(side)=1
            endif
           endif
          enddo ! dir2=1..sdim

          i_out=indexlo(1)-ii
          j_out=indexlo(2)-jj
          k_out=indexlo(SDIM)-kk

         else if (side.eq.2) then

          sideidx(dir_main)=indexhi(dir_main)+1
          do dir2=1,SDIM
           if (sideidx(dir2).gt.fabhi(dir2)) then
            if (velbc_in(dir2,side,dir2).ne.INT_DIR) then
             nbr_outside_domain_flag(side)=1
            endif
           endif
          enddo ! dir2=1..sdim

          i_out=indexhi(1)+ii
          j_out=indexhi(2)+jj
          k_out=indexhi(SDIM)+kk

         else 
          print *,"side invalid"
          stop
         endif

         if (nbr_outside_domain_flag(side).eq.0) then

          ! do nothing
 
         else if (nbr_outside_domain_flag(side).eq.1) then 

          if (velbc_in(dir_main,side,dir_main).eq.INT_DIR) then
           print *,"velbc_in bust "
           print *,"cen_outside_domain_flag= ",cen_outside_domain_flag
           print *,"nbr_outside_domain_flag= ",nbr_outside_domain_flag(side)
           stop
          endif

          if (operation_flag.eq.OP_GRADU_MAC_TO_CELL) then!grad u^T: MAC->CELL
         
             ! u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
           if ((velbc_in(dir_main,side,scomp_bc+nc-1).eq.REFLECT_EVEN).or. &
               (velbc_in(dir_main,side,scomp_bc+nc-1).eq.FOEXTRAP)) then
            if (side.eq.1) then
             local_data(1)=zero
            else if (side.eq.2) then
             local_data(bfact+1)=zero
            else
             print *,"side invalid"
             stop
            endif
           else if (velbc_in(dir_main,side,scomp_bc+nc-1).eq.REFLECT_ODD) then
            ! do nothing
           else if (velbc_in(dir_main,side,scomp_bc+nc-1).eq.EXT_DIR) then
            ! do nothing
           else
            print *,"velbc_in is corrupt"
            stop
           endif
          else if (operation_flag.eq.OP_RHS_CELL) then!RHS
           ! do nothing
          else if (operation_flag.eq.OP_ISCHEME_CELL) then!advection
           ! do nothing
          else if (operation_flag.eq.OP_DIV_CELL) then!divergence
           ! do nothing
          else if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then
           ! do nothing
          else
           print *,"operation_flag invalid29:",operation_flag
           stop
          endif

         else
          print *,"nbr_outside_domain_flag invalid 5"
          stop
         endif

        enddo ! side=1..2

         ! mask=0 piecewise FV  mask>0 SEM
        call line_MAC_TO_CELL(local_data,local_cell,local_div, &
         bfact,maskSEM,dx(dir_main))

        if (operation_flag.eq.OP_ISCHEME_CELL) then !advection
         if ((nc.ge.SEM_U+1).and.(nc.le.SEM_T+1)) then
          call line_MAC_TO_CELL(local_vel_data, &
           local_vel_cell,local_vel_div, &
           bfact,maskSEM,dx(dir_main))
          call line_MAC_TO_CELL(local_vel_data_div, &
           local_vel_cell_div,local_vel_div_div, &
           bfact,maskSEM,dx(dir_main))
         else
          print *,"nc invalid"
          stop
         endif
        else if ((operation_flag.eq.OP_RHS_CELL).or. & ! RHS
                 (operation_flag.eq.OP_DIV_CELL).or. & ! div u
                 (operation_flag.eq.OP_VEL_MAC_TO_CELL).or. & ! VEL MAC to CELL
                 (operation_flag.eq.OP_GRADU_MAC_TO_CELL)) then 
         ! do nothing
        else
         print *,"operation_flag invalid16991: ",operation_flag
         stop
        endif

        do isten=0,bfact-1

         indexlo(1)=i
         indexlo(2)=j
         if (SDIM.eq.3) then
          indexlo(SDIM)=k
         endif
         do dir2=1,SDIM
          indexmid(dir2)=indexlo(dir2)
         enddo ! dir2

         indexmid(dir_main)=indexlo(dir_main)+isten

         ic=indexmid(1)
         jc=indexmid(2)
         kc=indexmid(SDIM)

         RR=one
         RRTHETA=one

         call gridsten(xsten,xlo, &
          ic,jc,kc, &
          fablo,bfact,dx,nhalf)

         if (dir_main.eq.1) then ! r direction

          if (levelrz.eq.COORDSYS_CARTESIAN) then
           RR=one
          else if (levelrz.eq.COORDSYS_RZ) then
           if (SDIM.eq.2) then
            RR=xsten(0,1)
           else
            print *,"dimension bust"
            stop
           endif
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           RR=xsten(0,1)
          else
           print *,"levelrz invalid"
           stop
          endif

         else if ((dir_main.eq.2).or.(dir_main.eq.SDIM)) then
          RR=one
         else
          print *,"dir_main invalid sem mac to cell 4"
          stop
         endif 

         if ((dir_main.eq.2).and.(levelrz.eq.COORDSYS_CYLINDRICAL)) then ! theta direction
          RRTHETA=xsten(0,1)
         endif

         if (operation_flag.eq.OP_GRADU_MAC_TO_CELL) then ! interp grad U

          veldest(D_DECL(ic,jc,kc),dcomp+nc-1)=local_cell(isten+1)

         else if (operation_flag.eq.OP_ISCHEME_CELL) then ! advection

          if ((nc.ge.1).and.(nc.le.ncomp)) then

           if (dir_main.eq.1) then
            RR_divide=RR
           else if ((dir_main.eq.2).or.(dir_main.eq.SDIM)) then
            RR_divide=RRTHETA
           else
            print *,"dir_main invalid"
            stop
           endif

            ! local_data (the MAC data) is multiplied by RR above.

            ! velocity = div (umac u) - u div(u)
            ! temperature = div (umac T) - T div(u)
            ! mdotcell corresponds to VELADVECT_MF in the caller.
            ! maskdivres corresponds to DEN_RECON_MF in the caller.
           if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then ! velocity
            local_div(isten+1)=local_div(isten+1)- &
             mdotcell(D_DECL(ic,jc,kc),nc)*local_vel_div_div(isten+1)
           else if (nc.eq.SEM_T+1) then ! temperature
            local_div(isten+1)=local_div(isten+1)- & 
             maskdivres(D_DECL(ic,jc,kc),ibase+ENUM_TEMPERATUREVAR+1)* &
             local_vel_div_div(isten+1)
           else
            print *,"nc invalid"
            stop
           endif
           local_div(isten+1)=local_div(isten+1)/RR_divide

            ! divdest corresponds to rhs
            ! rhs corresponds to (*rhs)[mfi]
           if (dir_main.eq.1) then ! r direction
            divdest(D_DECL(ic,jc,kc),nc)=local_div(isten+1)
           else if ((dir_main.eq.2).or.(dir_main.eq.SDIM)) then
            divdest(D_DECL(ic,jc,kc),nc)= &
             divdest(D_DECL(ic,jc,kc),nc)+local_div(isten+1)
           else
            print *,"dir_main invalid mac to cell 5:",dir_main
            stop
           endif

           if (dir_main.eq.SDIM) then

            if (nc.eq.ncomp) then
      
             if (source_term.eq.SUB_OP_SDC_LOW_TIME) then 

              if ((slab_step.lt.0).or.(slab_step.gt.bfact_time_order)) then
               print *,"slab_step invalid"
               stop
              endif 

             else if (source_term.eq.SUB_OP_SDC_ISCHEME) then

              if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
               print *,"slab_step invalid"
               stop
              endif 

              do nc2=1,ncomp
               divflux(nc2)=divdest(D_DECL(ic,jc,kc),nc2)
              enddo ! nc2

              do nc2=1,SDIM
               vel_old=vel_old_fab(D_DECL(ic,jc,kc),nc2) 
               vel_update=vel_old-dt*divflux(nc2)

                ! SDC correction term: momentum
               if ((ns_time_order.ge.2).and. &
                   (ns_time_order.le.32).and. &
                   (advect_iter.eq.SUB_OP_ISCHEME_CORRECT).and. &
                   (SDC_outer_sweeps.gt.0).and. &
                   (divu_outer_sweeps+1.eq.num_divu_outer_sweeps)) then
                vel_update=vel_update-cterm(D_DECL(ic,jc,kc),nc2)
               else if ((ns_time_order.eq.1).or. &
                        (advect_iter.eq.SUB_OP_ISCHEME_PREDICT).or. &
                        (SDC_outer_sweeps.eq.0).or. &
                        (divu_outer_sweeps+1.lt.num_divu_outer_sweeps)) then
                ! do nothing
               else
                print *,"ns_time_order, SDC_outer_sweeps, or divu_outer.. bad"
                stop
               endif

               vel_new(nc2)=vel_update
              enddo ! nc2=1..sdim

              T_old=denold(D_DECL(ic,jc,kc),ibase+ENUM_TEMPERATUREVAR+1)
              T_new=T_old-dt*divflux(SEM_T+1)

                ! SDC correction term: temperature
              if ((ns_time_order.ge.2).and. &
                  (ns_time_order.le.32).and. &
                  (advect_iter.eq.SUB_OP_ISCHEME_CORRECT).and. &
                  (SDC_outer_sweeps.gt.0).and. &
                  (divu_outer_sweeps+1.eq.num_divu_outer_sweeps)) then
                 ! cterm corresponds to (*localMF[delta_MF])[mfi]
               T_new=T_new-cterm(D_DECL(ic,jc,kc),SEM_T+1)
              else if ((ns_time_order.eq.1).or. &
                       (advect_iter.eq.SUB_OP_ISCHEME_PREDICT).or. &
                       (SDC_outer_sweeps.eq.0).or. &
                       (divu_outer_sweeps+1.lt.num_divu_outer_sweeps)) then
               ! do nothing
              else
               print *,"ns_time_order, SDC_outer_sweeps, or divu_outer.. bad"
               stop
              endif

              if (T_new.lt.fort_tempcutoff(maskSEM)) then
               T_new=fort_tempcutoff(maskSEM)
              endif
              if (T_new.gt.fort_tempcutoffmax(maskSEM)) then
               T_new=fort_tempcutoffmax(maskSEM)
              endif 

              if (T_new.gt.zero) then
               ! do nothing
              else
               print *,"T_new underflow"
               print *,"dt=",dt
               print *,"divflux (u dot grad T)= ",divflux(SEM_T+1)
               print *,"cterm= ",cterm(D_DECL(ic,jc,kc),SEM_T+1)
               print *,"energyflag (advect_iter) =",energyflag
               print *,"ic,jc,kc ",ic,jc,kc
               stop
              endif

              do nc2=1,SDIM
               veldest(D_DECL(ic,jc,kc),nc2)=vel_new(nc2)
              enddo ! nc2
              dendest(D_DECL(ic,jc,kc),ibase+ENUM_TEMPERATUREVAR+1)=T_new

             else
              print *,"source_term invalid"
              stop
             endif

            else if ((nc.ge.1).and.(nc.lt.ncomp)) then
             ! do nothing
            else
             print *,"nc invalid"
             stop    
            endif 

           else if ((dir_main.ge.1).and.(dir_main.lt.SDIM)) then
            ! do nothing
           else
            print *,"dir_main invalid"
            stop    
           endif 

          else
           print *,"nc out of range"
           stop
          endif 

         else if (operation_flag.eq.OP_RHS_CELL) then ! RHS

          if (maskcoef(D_DECL(ic,jc,kc),1).eq.one) then ! not covered

           if (dir_main.eq.1) then
            divdest(D_DECL(ic,jc,kc),nc)= &
             local_div(isten+1)/RR
           else if ((dir_main.ge.2).and. &
                    (dir_main.le.SDIM)) then
            divdest(D_DECL(ic,jc,kc),nc)= &
             divdest(D_DECL(ic,jc,kc),nc)+ &
             local_div(isten+1)/RRTHETA
           else
            print *,"dir_main invalid mac to cell 6"
            stop
           endif

           if (dir_main.eq.SDIM) then
            VOLTERM=vol(D_DECL(ic,jc,kc),1)
            if (VOLTERM.gt.zero) then
             ! do nothing
            else
             print *,"VOLTERM invalid"
             stop
            endif
            divu=divdest(D_DECL(ic,jc,kc),nc)*VOLTERM
            CC=cterm(D_DECL(ic,jc,kc),nc) ! already x VOLTERM
            CC_DUAL=veldest(D_DECL(ic,jc,kc),nc) ! already x VOLTERM
            if (CC_DUAL.eq.CC) then
             ! do nothing
            else
             print *,"CC_DUAL invalid"
             stop
            endif

            MDOT=mdotcell(D_DECL(ic,jc,kc),nc)

            if (CC.ge.zero) then
             ! do nothing
            else
             print *,"CC invalid in SEM_MAC_TO_CELL"
             stop
            endif

            if (1.eq.0) then
             if (CC.eq.zero) then
              ! do nothing
             else
              print *,"CC should be 0 for incompressible fluids"
              stop
             endif
             if (MDOT.ne.zero) then
              print *,"expecting MDOT=0   MDOT=",MDOT
              print *,"ic,jc,kc,nc ",ic,jc,kc,nc
              stop
             endif
            endif

             ! divu=-dt VOLTERM * div(k grad T)  project_option==SOLVETYPE_HEAT
             ! divu=-dt VOLTERM * visc_coef div(2 mu D) 
             !   project_option==SOLVETYPE_VISC
             ! use_dt=1 dir=-1
             ! use_HO=1
             ! enable_spectral=1
             ! uncoupled_viscosity=1
            local_div_val=divu/VOLTERM

            local_POLD=pold(D_DECL(ic,jc,kc),nc)
            local_POLD_DUAL=dendest(D_DECL(ic,jc,kc),nc)

            if (homflag.eq.0) then
             if (local_POLD.eq.local_POLD_DUAL) then
              RHS=local_POLD*CC
             else
              print *,"local_POLD invalid"
              stop
             endif
            else if (homflag.eq.1) then
             if (local_POLD.eq.local_POLD_DUAL) then
              RHS=local_POLD*CC
             else
              print *,"local_POLD invalid"
              stop
             endif
            else if (homflag.eq.2) then
             if (local_POLD.eq.local_POLD_DUAL) then
              RHS=-local_POLD*CC
             else
              print *,"local_POLD invalid"
              stop
             endif
            else if (homflag.eq.3) then
             if ((local_POLD.eq.zero).and.(local_POLD_DUAL.eq.zero)) then
              RHS=zero
             else
              print *,"local_POLD or local_POLD_DUAL invalid"
              stop
             endif
            else if (homflag.eq.4) then
             if (local_POLD.eq.local_POLD_DUAL) then
              RHS=zero
             else
              print *,"local_POLD invalid"
              stop
             endif
            else
             print *,"homflag invalid 10"
             stop
            endif

            if (homflag.eq.0) then
             RHS=RHS-divu/dt+MDOT
            else if (homflag.eq.1) then
             RHS=RHS+divu/dt
            else if (homflag.eq.2) then
             RHS=RHS-divu/dt+MDOT
            else if (homflag.eq.3) then
             RHS=-divu/dt+MDOT
             if (level.eq.finest_level) then
              if (divu.eq.zero) then
               ! do nothing
              else 
               print *,"divu invalid"
               stop
              endif
             else if ((level.ge.0).and.(level.lt.finest_level)) then
              ! do nothing
             else
              print *,"level invalid"
              stop
             endif
            else if (homflag.eq.4) then
             if (RHS.eq.zero) then
              RHS=divu/VOLTERM
             else
              print *,"RHS invalid"
              stop
             endif
            else
             print *,"homflag invalid"
             stop
            endif

            divdest(D_DECL(ic,jc,kc),nc)=RHS
           else if ((dir_main.ge.1).and.(dir_main.lt.SDIM)) then
            ! do nothing
           else
            print *,"dir_main invalid mac_to_cell 7"
            stop
           endif 

          else if (maskcoef(D_DECL(ic,jc,kc),1).eq.zero) then
           ! do nothing (covered)
          else 
           print *,"maskcoef invalid"
           stop
          endif

         else if (operation_flag.eq.OP_DIV_CELL) then ! div

          if (dir_main.eq.1) then
           divdest(D_DECL(ic,jc,kc),1)= &
            local_div(isten+1)/RR
          else if ((dir_main.eq.2).or. &
                   (dir_main.eq.SDIM)) then
           divdest(D_DECL(ic,jc,kc),1)= &
            divdest(D_DECL(ic,jc,kc),1)+ &
            local_div(isten+1)/RRTHETA
          else
           print *,"dir_main invalid sem mac to cell 9"
           stop
          endif

         else if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then
           ! local_data (the MAC data) is multiplied by RR above.
          veldest(D_DECL(ic,jc,kc),dir_main)= &
           local_cell(isten+1)/RR
         else
          print *,"operation_flag invalid30:",operation_flag
          stop
         endif
         
        enddo ! isten=0...bfact-1

       else
        print *,"cen_outside_domain_flag invalid"
        stop
       endif 

      enddo ! nc=1..ncomp

      return
      end subroutine SEM_MAC_TO_CELL

      subroutine SUB_clamped_LS(x,t,LS,vel,temperature,prescribed_flag,dx)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(out) :: LS
      real(amrex_real), INTENT(out) :: vel(SDIM)
      real(amrex_real), INTENT(out) :: temperature
      integer, INTENT(out) :: prescribed_flag

      integer dir

      call SUB_clamped_LS_no_scale(x,t,LS,vel,temperature,prescribed_flag,dx)
      do dir=1,SDIM
       if (vel_homflag.eq.1) then
        vel(dir)=zero
       else if (vel_homflag.eq.0) then
        if (global_velocity_scale.gt.zero) then
         vel(dir)=vel(dir)/global_velocity_scale
        else
         print *,"global_velocity_scale invalid"
         stop
        endif
       else
        print *,"vel_homflag invalid"
        stop
       endif
      enddo ! dir=1..sdim

      return
      end subroutine SUB_clamped_LS 

       ! 1<=dir<=sdim  1<=side<=2
       ! 1<=veldir<=sdim
       ! input: vel_in
       ! output: vel_in
      subroutine velbc_override(time,dir,side,veldir,vel_in, &
        xsten,nhalf,dx,bfact)
      use global_utility_module
      use rainControl_module
      use bubbleControl_module
      use hydrateReactor_module
      use unimaterialChannel_module
      use River
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use rigid_FSI_module
      use sinking_particle_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: dir,side,veldir,bfact
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) vel
      real(amrex_real), INTENT(inout) :: vel_in
      real(amrex_real) x,y,z,r
      real(amrex_real), INTENT(in) ::time
      real(amrex_real) x_vel,y_vel,z_vel
      real(amrex_real) dist
      integer error
      real(amrex_real) costheta,sintheta
      real(amrex_real) xprime,yprime,zprime
      real(amrex_real) velcell(SDIM)
      integer dir2
      real(amrex_real) velx_rain,vely_rain
      real(amrex_real) xvec(SDIM)
      integer local_dir
      real(amrex_real) xwall
      real(amrex_real) local_LS(num_materials)
      real(amrex_real), parameter :: stub_zero=zero
      
      velx_rain=adv_vel
      vely_rain=uDrop

      if (nhalf.lt.3) then
       print *,"nhalf invalid velbc override"
       print *,"nhalf= ",nhalf
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)
      do local_dir=1,SDIM
       xvec(local_dir)=xsten(0,local_dir)
      enddo
      xwall=xvec(dir)

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y in 30916"
        stop
       endif
      endif

      vel=vel_in*global_velocity_scale

      if (vel_homflag.eq.1) then
       vel_in=zero
      else if (vel_homflag.eq.0) then

       do dir2=1,SDIM
        velcell(dir2)=zero
       enddo

       if (is_in_probtype_list().eq.1) then
        call SUB_LS(xvec,time,local_LS,num_materials)
        call SUB_VEL_BC(xwall,xvec,time,local_LS, &
         velcell(veldir),vel,veldir,dir,side,dx,num_materials)
       else if (probtype.eq.401) then
        call HELIX_LS(xvec,time,local_LS)
        call HELIX_VEL_BC(xwall,xvec,time,local_LS, &
         velcell(veldir),vel,veldir,dir,side,dx)
       else if (probtype.eq.402) then
        call TSPRAY_LS(xvec,time,local_LS)
        call TSPRAY_VEL_BC(xwall,xvec,time,local_LS, &
         velcell(veldir),vel,veldir,dir,side,dx)
       else if (probtype.eq.412) then ! step
        call CAV2Dstep_LS(xvec,time,local_LS)
        call CAV2Dstep_VEL_BC(xwall,xvec,time,local_LS, &
         velcell(veldir),vel,veldir,dir,side,dx)

       else if (probtype.eq.413) then ! ZEYU droplet impact
        call ZEYU_droplet_impact_LS(xvec,time,local_LS)
        ! pass dx
        call ZEYU_droplet_impact_VEL_BC(xwall,xvec,time,local_LS, &
         velcell(veldir),vel,veldir,dir,side,dx)

       else if (probtype.eq.533) then
        call rigid_FSI_LS(xvec,time,local_LS)
        call rigid_FSI_VEL_BC(xwall,xvec,time,local_LS, &
         velcell(veldir),vel,veldir,dir,side,dx)
       else if (probtype.eq.534) then
        call sinking_FSI_LS(xvec,time,local_LS)
        call sinking_FSI_VEL_BC(xwall,xvec,time,local_LS, &
         velcell(veldir),vel,veldir,dir,side,dx)

       else if (probtype.eq.311) then ! user defined

        call USERDEF_LS(xvec,time,local_LS)
        call USERDEF_VEL_BC(xwall,xvec,time,local_LS, &
         velcell(veldir),vel,veldir,dir,side,dx)

        ! velbc_override
       else if ((probtype.eq.26).and. &
                ((axis_dir.eq.10).or. & ! BCG test
                 (axis_dir.eq.11).or. & !BCG_periodic test (shouldn't be here)
                 (axis_dir.eq.12))) then!buoyancy test 
        do dir2=1,SDIM
         velcell(dir2)=zero
        enddo

        ! in: velbc_override
       else if (probtype.eq.199) then ! hydrates

        call HYD_VELO_BC(time,dir,side,veldir, &
         velcell(veldir),vel,x,y,z,dx)

        ! in: velbc_override
       else if (probtype.eq.220) then ! UNIMAT problem

        call UNIMAT_VELO_BC(time,dir,side,veldir, &
         velcell(veldir),vel,x,y,z,dx)

       else if ((probtype.eq.299).or. &
                (probtype.eq.301)) then !melting (boundary condition velocity)
        do dir2=1,SDIM
         velcell(dir2)=zero
        enddo
       else if (probtype.eq.209) then ! River
        if (dir.eq.SDIM) then
         call RiverVelocity(x,y,z,velcell,axis_dir,probloz,probhiz)
        else
         print *,"dir invalid in velbc_override"
         stop
        endif
       else if ((probtype.eq.1).and.(axis_dir.eq.15)) then
        call Zuzio_velbc(time,dir,side,veldir, &
         velcell(veldir),vel)
       else 
 
        if (side.eq.1) then
  
           ! xlo
         if (dir.eq.1) then
   
           ! xlo,xvel
          if (veldir.eq.1) then
  
           if (adv_dir.eq.1)then
            call rampvel(time,x_vel)  ! default is adv_vel
           else  
            x_vel = zero
           endif
  
           velcell(veldir)=x_vel  ! xlo, velx

           if (probtype.eq.bubbleInPackedColumn) then ! xlo,velx
            call Pack_velbc(dir,side,veldir,velcell)
           else if (SDIM.eq.2) then
  
            if ((probtype.eq.1).and.(axis_dir.eq.11)) then
             velcell(veldir)=vinletgas*(y/yblob-one) 
            else if (probtype.eq.58) then
             velcell(veldir)=zero
            else if (probtype.eq.802) then ! dissolution, xlo,velx
             print *,"802 obsolete"
             stop
            else if (probtype.eq.110) then
             call get_bump_velocity(xsten,nhalf,dx,bfact,velcell(veldir),time)  ! xvel,xlo 
            else if (probtype.eq.59) then ! xvel,xlo,velbc_override 2d
             ! vel=solidvel in solid
             call mask_velocity(xsten,nhalf,dx,bfact,velcell,time)
            else if (probtype.eq.532) then
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell) ! xvel,xlo
            else if (probtype.eq.5700) then  ! xvel,xlo
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            endif

           else if (SDIM.eq.3) then

            if (probtype.eq.59) then ! xvel,xlo,velbc_override 3d
               ! vel=solvel in sol
             call mask_velocity(xsten,nhalf,dx,bfact,velcell,time) 
            else if (probtype.eq.5501) then  ! xvel,xlo,velbc_override
               ! vel=solvel in sol
             call mask_velocity(xsten,nhalf,dx,bfact,velcell,time) 
   
! XIAOYI LI : Adding turbulent boundary layer for gas flow          
            else if (probtype.eq.53) then
!          call BL_inlet(velcell(veldir),z,time)  
            else if (probtype.eq.529) then ! airblast inlet
             r = sqrt((y-yblob)**2+(z-zblob)**2)
             if( (r.le.radblob2) .and. (r.ge.radblob+twall) ) then ! gas co-flow
              velcell(veldir) = adv_vel
             else if (r.le.radblob) then ! liquid injection
              velcell(veldir) = advbot
             else
              velcell(veldir) = zero
             endif
   
! gear problem
            else if (probtype.eq.563) then
             if (levelrz.eq.COORDSYS_CARTESIAN) then
              if (radblob.gt.zero) then
               if (sqrt((y-yblob)**2+(z-zblob)**2).le.radblob) then
                velcell(veldir)=advbot
               endif
              endif
             else
              print *,"levelrz invalid velbc override"
              stop
             endif
            else if (probtype.eq.58) then
              velcell(veldir)=zero
            else if (probtype.eq.51) then
              velcell(veldir)=zero
            else if (probtype.eq.50) then
              if ((z.lt.zblob-radblob).or.(z.gt.zblob+radblob).or. &
                  (y.lt.yblob-radblob).or.(y.gt.yblob+radblob)) then
               velcell(veldir)=zero
              endif
! xlo, velx 3D
            else if (probtype.eq.41) then 
             velcell(veldir)=zero
            endif
! microfluidics problem horizontal vel xlo inflow
            if (probtype.eq.5700) then  ! xvel,xlo
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            endif

           else
            print *,"dimension bust"
            stop
           endif
  
          else if (veldir.eq.2) then ! xlo, vely
  
           if (adv_dir .eq. 2) then
            y_vel = adv_vel
           else
            y_vel = zero
           endif
           if ((probtype.eq.53).and.(SDIM.eq.2)) then  
            ! y velocity, left wall
            y_vel = zero
           endif
  
           velcell(veldir)=y_vel  ! xlo, vely

           if (probtype.eq.bubbleInPackedColumn) then ! xlo,vely
            call Pack_velbc(dir,side,veldir,velcell)
           else if (SDIM.eq.2) then
  
! xlo wall vertical pipe velocity
            if (probtype.eq.41) then
             velcell(veldir)=zero
            else if (probtype.eq.532) then ! xlo,yvel
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            else if (probtype.eq.802) then ! dissolution, xlo,vely
             print *,"802 obsolete"
             stop
            endif

           else if (SDIM.eq.3) then

            ! do nothing

           else
            print *,"dimension bust"
            stop
           endif

! dir=1, side=1  (xlo face of domain)
          else if ((veldir.eq.3).and.(SDIM.eq.3)) then
   
           if (adv_dir.eq.3) then
            z_vel = adv_vel
           else  
            z_vel = zero
           endif
           if ((probtype.eq.53).or.(probtype.eq.530).or. &
                (probtype.eq.531).or. &
                (probtype.eq.536).or.(probtype.eq.537).or. &
                (probtype.eq.538).or.(probtype.eq.541)) then
            z_vel=zero
           endif
           velcell(veldir)=z_vel  ! xlo,velz
  
! velz, xlo
           if ((probtype.eq.36).and. &
               ((yblob9.ne.zero).or.(yblob10.ne.zero))) then
            velcell(veldir)=yblob9
           else if (probtype.eq.bubbleInPackedColumn) then ! xlo,velz
            call Pack_velbc(dir,side,veldir,velcell)
           endif
          else
           print *,"veldir invalid"
           stop
          endif
  
         else if (dir.eq.2) then ! ylo
  
          if (veldir.eq.1) then
  
           if (adv_dir.eq.1)then
            call rampvel(time,x_vel)
           else  
            x_vel = zero
           endif
  
           velcell(veldir)=x_vel  ! ylo, velx

           if (probtype.eq.bubbleInPackedColumn) then ! ylo,velx
            call Pack_velbc(dir,side,veldir,velcell)
           else if (SDIM.eq.2) then

              ! ylo, velx, 2D
            if ((probtype.eq.701).and.(axis_dir.eq.2)) then
             velcell(veldir)=adv_vel 
            else if ((probtype.eq.1).and.(axis_dir.eq.13)) then
             velcell(veldir)=zero
            else if ((probtype.eq.1).and.(axis_dir.ne.14)) then
             velcell(veldir)=-vinletgas
            else if (probtype.eq.21) then
             velcell(veldir)=adv_vel
            else if (probtype.eq.58) then
             velcell(veldir)=zero
            else if (probtype.eq.53) then  ! x vel bottom wall
             velcell(veldir)=zero
              ! 2d, ylo, velx bottom wall velbc_override
            else if ((probtype.eq.538).or.(probtype.eq.541)) then
             velcell(veldir)=zero
            else if (probtype.eq.10) then
             if ((x.ge.xblob-radblob).and.(x.le.xblob+radblob)) then
              velcell(veldir)=zero
             endif
            endif

           else if (SDIM.eq.3) then

            if (probtype.eq.58) then
             velcell(veldir)=zero
            endif

           else
            print *,"dimension bust"
            stop
           endif
  
          else if (veldir.eq.2) then  ! ylo, vely
  
           if (adv_dir .eq. 2) then
            y_vel = adv_vel
           else
            y_vel = zero
           endif
  
           velcell(veldir)=y_vel   ! y velocity bottom wall (ylo, vely)

           if (probtype.eq.bubbleInPackedColumn) then ! ylo,vely
            call Pack_velbc(dir,side,veldir,velcell)
           else if (SDIM.eq.2) then
  
             ! side=1, dir=2,veldir=2  (vely at ylo)
             ! bubble formation
             ! bubble_formation_inflow_bc sets of Poiseuille flow
             ! quantity in cell is average velocity for Poiseuille flow
            if ((probtype.eq.25).and.(axis_dir.gt.0)) then
             call bubble_formation_inflow_bc(xsten,nhalf,x,velcell(veldir))
            else if (probtype.eq.22) then
             velcell(veldir)=y_vel*(xblob**2-x**2)*three/ &
                   (two*xblob*xblob*xblob)
             if (x.gt.xblob) then
              velcell(veldir)=zero
             endif
            else if ((probtype.eq.53).or.(probtype.eq.537)) then  ! ylo, vely 2D
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            else if ((probtype.eq.538).or.(probtype.eq.541)) then  ! ylo, vely
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            else if (probtype.eq.539) then  ! ylo, vely - velbc_override
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            else if ((probtype.eq.63).or.(probtype.eq.64)) then
             velcell(veldir)=y_vel
            else if (probtype.eq.10) then
             velcell(veldir)=zero
             if ((x.ge.xblob-radblob).and. &
                 (x.le.xblob+radblob)) then
              velcell(veldir)=advbot
             endif
            else if ((probtype.eq.36).and.(xblob10.gt.zero).and. &
                     (yblob10.ne.zero)) then
             velcell(veldir)=x*yblob10/xblob10
! ylo vertical velocity
            else if (probtype.eq.41) then
             call get_pipe_velocity(xsten,nhalf,dx,bfact,velcell,time)
            else if (probtype.eq.72) then
             call vapordist(xsten,nhalf,dx,bfact,dist)
             if (dist.ge.zero) then
              velcell(veldir)=advbot
             else
              velcell(veldir)=adv_vel
             endif
! microfluidics squeeze vertical vel ylo inflow
            else if (probtype.eq.5700) then  ! yvel,ylo
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            endif

           else if (SDIM.eq.3) then

            ! ylo, yvel
            if (probtype.eq.532) then
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            else if (probtype.eq.59) then
               ! vel=solvel in sol
             call mask_velocity(xsten,nhalf,dx,bfact,velcell,time) 
            else if (probtype.eq.5501) then
               ! vel=solvel in sol
             call mask_velocity(xsten,nhalf,dx,bfact,velcell,time) 
! microfluidics squeeze vertical vel ylo inflow
            else if (probtype.eq.5700) then ! yvel,ylo
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            endif

           else
            print *,"dimension bust"
            stop
           endif

! dir=2, side=1 (ylo face) zvel yloface
          else if ((veldir.eq.3).and.(SDIM.eq.3)) then
  
           if (adv_dir.eq.3) then
            z_vel = adv_vel
           else  
            z_vel = zero
           endif
            ! ylo,zvel
           if (probtype.eq.bubbleInPackedColumn) then ! ylo,velz
            call Pack_velbc(dir,side,veldir,velcell)
           else if (probtype.eq.532) then
            call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
           else if ((probtype.eq.53).or.(probtype.eq.530).or. &
               (probtype.eq.531).or. &
               (probtype.eq.536).or.(probtype.eq.537).or. &
               (probtype.eq.538).or. &
               (probtype.eq.541)) then
            z_vel=zero
            velcell(veldir)=z_vel  ! ylo,velz
           else
            velcell(veldir)=z_vel  ! ylo,velz
           endif
  
          else
           print *,"veldir invalid"
           stop
          endif

         else if ((dir.eq.3).and.(SDIM.eq.3)) then
  
! zlo face, x velocity
          if (veldir.eq.1) then
  
           if (adv_dir.eq.1)then
            call rampvel(time,x_vel)
           else  
            x_vel = zero
           endif
           velcell(veldir)=x_vel  ! zlo, velx
  
           if (probtype.eq.bubbleInPackedColumn) then ! zlo,velx
            call Pack_velbc(dir,side,veldir,velcell)
           else if ((probtype.eq.1).and.(axis_dir.eq.13)) then
            velcell(veldir)=zero
           else if ((probtype.eq.1).and.(axis_dir.lt.150)) then
            if ((axis_dir.eq.11).or.(axis_dir.eq.12)) then
             if (xblob10.eq.one) then
              velcell(veldir)=radblob10
             else
              velcell(veldir)=-vinletgas
             endif
            else
             velcell(veldir)=-vinletgas
            endif
! zlo face, x velocity
           else if (probtype.eq.530) then
            call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
           else if ((probtype.eq.58).or.(probtype.eq.53).or. &
                    (probtype.eq.531).or. &
                    (probtype.eq.536).or.(probtype.eq.537).or. &
                    (probtype.eq.538).or. &
                    (probtype.eq.541)) then
            velcell(veldir)=zero
            if ((probtype.eq.541).and.(axis_dir.eq.0)) then
               !  zlo,velx
             call velread_bc_point(x,y,z,dir,veldir,0,time,velcell(dir))
            endif
! wall gas pipe horizontal velocity (velx) z=zlo 
           else if (probtype.eq.41) then
            print *,"modify for pipe 3d"
            stop
           endif
  
! zlo face, y velocity
  
          else if (veldir.eq.2) then
  
           if (adv_dir .eq. 2) then
            y_vel = adv_vel
           else
            y_vel = zero
           endif
           velcell(veldir)=y_vel  ! zlo, vely
  
           if (probtype.eq.bubbleInPackedColumn) then ! zlo,vely
            call Pack_velbc(dir,side,veldir,velcell)
! zlo face, y velocity
           else if (probtype.eq.530) then
            call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
           else if ((probtype.eq.541).and.(axis_dir.eq.0)) then
            call velread_bc_point(x,y,z,dir,veldir,0,time,velcell(dir))
           endif
  
          else if ((veldir.eq.3).and.(SDIM.eq.3)) then
  
! zlo face, z velocity
           if (adv_dir.eq.3) then
            z_vel = adv_vel
           else  
            z_vel = zero
           endif
           velcell(veldir)=z_vel  ! zlo, velz
  
! zlo, velz 3D, velbc_override
           if (probtype.eq.bubbleInPackedColumn) then ! zlo,velz
            call Pack_velbc(dir,side,veldir,velcell)
           else if ((probtype.eq.53).or.(probtype.eq.536).or. &
               (probtype.eq.537).or.(probtype.eq.530).or. &
               (probtype.eq.538).or.(probtype.eq.541)) then

            call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            if ((probtype.eq.541).and.(axis_dir.eq.0)) then
             call velread_bc_point(x,y,z,dir,veldir,0,time,velcell(dir))
            endif


           else if (probtype.eq.531) then

            z_vel=zero
            velcell(veldir)=z_vel  ! zlo, velz

           else if ((probtype.eq.63).or.(probtype.eq.64)) then
            velcell(veldir)=z_vel
! velz, zlo
           else if ((probtype.eq.36).and.(xblob10.gt.zero).and. &
                    ((yblob9.ne.zero).or.(yblob10.ne.zero))) then
            if (probhix-problox.le.zero) then
             print *,"probhix or problox invalid"
             stop
            endif
            velcell(veldir)=yblob9+(x-problox)*(yblob10-yblob9)/  &
              (probhix-problox)
           endif
  
          else
           print *,"veldir invalid"
           stop
          endif
  
         else
          print *,"dir invalid velbcoverride 2"
          stop
         endif
  
! -------------------------------side=2 ----------------------------
        else if (side.eq.2) then
  
         if (dir.eq.1) then
  
! uvel at xhi
          if (veldir.eq.1) then
  
           if (adv_dir.eq.1)then
            call rampvel(time,x_vel)
           else  
            x_vel = zero
           endif
  
           velcell(veldir)=x_vel  ! xhi, velx

           if (probtype.eq.bubbleInPackedColumn) then ! xhi,velx
            call Pack_velbc(dir,side,veldir,velcell)
           else if (SDIM.eq.2) then
  
              ! xhi,xvel
            if (probtype.eq.532) then
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            else if ((probtype.eq.1).and.(axis_dir.eq.11)) then
             velcell(veldir)=vinletgas*(y/yblob-one) 
            else if (probtype.eq.58) then
             velcell(veldir)=zero
            else if (probtype.eq.59) then ! xvel,xhi,velbc_override 2d
             ! vel=solidvel in solid
             call mask_velocity(xsten,nhalf,dx,bfact,velcell,time)
            else if ((probtype.eq.22).and.(axis_dir.eq.13)) then
             call vbc(velcell(veldir),time,y,stub_zero,error)
            else if (probtype.eq.32) then
             if (levelrz.ne.COORDSYS_CARTESIAN) then
              print *,"levelrz invalid for probtype=32"
              stop
             endif
            else if (probtype.eq.110) then
             call get_bump_velocity(xsten,nhalf,dx,bfact,velcell(veldir),time) ! xvel,xhi
            endif
            if (probtype.eq.5700) then  ! xvel,xhi
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            endif

           else if (SDIM.eq.3) then

            if (probtype.eq.58) then
             velcell(veldir)=zero
            else if ((axis_dir.eq.13).and.(probtype.eq.22)) then
             call vbc(velcell(veldir),time,z,y,error)
            endif
            if (probtype.eq.5700) then  ! xhi,velx
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            endif

           else
            print *,"dimension bust"
            stop
           endif
  
! yvel at xhi
          else if (veldir.eq.2) then
  
           if (adv_dir .eq. 2) then
            y_vel = adv_vel
           else
            y_vel = zero
           endif
            ! side=2  veldir=2 dir=1 (yvel, right wall)
           if ((probtype.eq.53).and.(SDIM.eq.2)) then 
            y_vel = zero
           endif
  
           velcell(veldir)=y_vel  ! xhi, vely

           if (probtype.eq.bubbleInPackedColumn) then ! xhi,vely
            call Pack_velbc(dir,side,veldir,velcell)
           else if (SDIM.eq.2) then
  
            if ((probtype.eq.36).and.(yblob10.ne.zero)) then
             velcell(veldir)=yblob10
! xhi wall vertical pipe velocity
            else if (probtype.eq.41) then
             velcell(veldir)=zero
            else if (probtype.eq.532) then
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell) ! yvel,xhi
            else if (probtype.eq.32) then 
             ! do nothing ?
            endif
            if (probtype.eq.5700) then  ! yvel,xhi
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            endif

           else if (SDIM.eq.3) then

            if (probtype.eq.5700) then  ! xhi,vely
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            endif

           else
            print *,"dimension bust"
            stop
           endif

! xhi face, z vel, velbc_override
          else if ((veldir.eq.3).and.(SDIM.eq.3)) then
  
           if (adv_dir.eq.3) then
            z_vel = adv_vel
           else  
            z_vel = zero
           endif
           if (probtype.eq.bubbleInPackedColumn) then ! xhi,velz
            call Pack_velbc(dir,side,veldir,velcell)
           else if ((probtype.eq.53).or.(probtype.eq.530).or. &
               (probtype.eq.531).or. &
               (probtype.eq.536).or.(probtype.eq.537).or. &
               (probtype.eq.538).or.(probtype.eq.541)) then
            z_vel=zero
           endif
           velcell(veldir)=z_vel  ! xhi, velz
  
! velz,xhi 
           if ((probtype.eq.36).and. &
               ((yblob9.ne.zero).or.(yblob10.ne.zero))) then
            velcell(veldir)=yblob10
           endif

          else
           print *,"veldir invalid"
           stop
          endif

          ! yhi face 
         else if (dir.eq.2) then
 
           ! yhi, velx
          if (veldir.eq.1) then
  
           if (adv_dir.eq.1)then
            call rampvel(time,x_vel)
           else  
            x_vel = zero
           endif
  
           velcell(veldir)=x_vel ! yhi, velx

           if (probtype.eq.bubbleInPackedColumn) then ! yhi,velx
            call Pack_velbc(dir,side,veldir,velcell)
           else if (SDIM.eq.2) then

             ! yhi, velx, 2D
            if ((probtype.eq.701).and.(axis_dir.eq.2)) then 
             velcell(veldir)=adv_vel
            else if ((probtype.eq.1).and.(axis_dir.ne.14)) then
             velcell(veldir)=vinletgas
            else if (probtype.eq.21) then
             velcell(veldir)=adv_vel
            else if (probtype.eq.58) then
             velcell(veldir)=zero
            endif

           else if (SDIM.eq.3) then

            if (probtype.eq.58) then
             velcell(veldir)=zero
            endif

           else
            print *,"dimension bust"
            stop
           endif

            ! yvel, yhi face 
          else if (veldir.eq.2) then
  
           if (adv_dir .eq. 2) then
            y_vel = adv_vel
           else
            y_vel = zero
           endif
           if (((probtype.eq.53).or. &
                (probtype.eq.538).or. &
                (probtype.eq.541)).and. &
               (SDIM.eq.2)) then 
            ! side=2 dir=2 veldir=2 (yvel top wall)
            y_vel = zero
           endif
  
           velcell(veldir)=y_vel  ! yhi, vely

           if (probtype.eq.bubbleInPackedColumn) then ! yhi,vely
            call Pack_velbc(dir,side,veldir,velcell)
           else if (SDIM.eq.2) then
  
            if ((probtype.eq.16).or. &
                ((probtype.eq.25).and.(axis_dir.eq.0)) ) then
             velcell(veldir)=zero
             if ((x.ge.xblob-radblob).and. &
                 (x.le.xblob+radblob)) then
              velcell(veldir)=-abs(advbot)
             endif
! yhi vertical velocity
            else if (probtype.eq.41) then
             call get_pipe_velocity(xsten,nhalf,dx,bfact,velcell,time)
            else if (probtype.eq.72) then
   
            else if (probtype.eq.62) then
             costheta=cos(xblob2)
             sintheta=sin(xblob2)
             xprime=costheta*(x-xblob)-sintheta*(y-yblob)
             yprime=sintheta*(x-xblob)+costheta*(y-yblob)
             if (xprime**2<radblob**2) then
              velcell(veldir)=-radblob3
             endif
            else if ((probtype.eq.63).or.(probtype.eq.64)) then
             velcell(veldir)=four*y_vel
            else if ((probtype.eq.36).and.(xblob10.gt.zero).and. &
                     (yblob10.ne.zero)) then
             velcell(veldir)=x*yblob10/xblob10
! microfluidics problem yvel yhi 
! microfluidics channel -- 0,3,4 Roper, 1 Comsol, 2 squeeze
            else if (probtype.eq.5700) then  ! yvel,yhi
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)
            else if (probtype.eq.701) then
             if ((axis_dir.eq.0).or.(axis_dir.eq.1)) then
              call get_rain_velocity(x,y,z,dx,velcell(veldir), &
                      vely_rain,time,dir)  ! yhi veldir=2
             else if (axis_dir.eq.2) then
              velcell(veldir)=zero  ! yhi veldir=2
             else
              print *,"axis_dir invalid"
              stop
             endif

            endif

           else if (SDIM.eq.3) then

             ! yhi,yvel,3D
            if (probtype.eq.532) then
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            else if (probtype.eq.531) then
             call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
            endif
! microfluidics problem yvel yhi,3D 
! microfluidics channel -- 0,3,4 Roper, 1 Comsol, 2,5 squeeze
            if (probtype.eq.5700) then  ! yvel,yhi
             call microfluidic_velbc(xsten,nhalf,dir,side,velcell)     
            endif

           else
            print *,"dimension bust"
            stop
           endif

! yhi, zvel
! z vel, yhi face
          else if ((veldir.eq.3).and.(SDIM.eq.3)) then
  
           if (adv_dir.eq.3) then
            z_vel = adv_vel
           else  
            z_vel = zero
           endif 
             ! yhi,zvel
           if (probtype.eq.bubbleInPackedColumn) then ! yhi,zvel
            call Pack_velbc(dir,side,veldir,velcell)
           else if (probtype.eq.532) then
            call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
           else if ((probtype.eq.53).or.(probtype.eq.530).or. &
               (probtype.eq.531).or. &
               (probtype.eq.536).or.(probtype.eq.537).or. &
               (probtype.eq.538).or.(probtype.eq.541)) then
            z_vel=zero
            velcell(veldir)=z_vel  ! yhi, velz
           else
            velcell(veldir)=z_vel  ! yhi, velz
           endif
  
          else
           print *,"veldir invalid"
           stop
          endif

         else if ((dir.eq.3).and.(SDIM.eq.3)) then
  
          if (veldir.eq.1) then
  
           if (adv_dir.eq.1)then
            call rampvel(time,x_vel)
           else  
            x_vel = zero
           endif
           velcell(veldir)=x_vel  ! zhi, velx
  
           if (probtype.eq.bubbleInPackedColumn) then ! zhi,velx
            call Pack_velbc(dir,side,veldir,velcell)
           else if ((probtype.eq.1).and.(axis_dir.lt.150)) then
            velcell(veldir)=vinletgas
           else if (probtype.eq.58) then
            velcell(veldir)=zero
           endif
  
          else if (veldir.eq.2) then
  
           if (adv_dir .eq. 2) then
            y_vel = adv_vel
           else
            y_vel = zero
           endif
           velcell(veldir)=y_vel  ! zhi, vely
  
           if (probtype.eq.bubbleInPackedColumn) then ! zhi,vely
            call Pack_velbc(dir,side,veldir,velcell)
           endif
  
! zhi, velz
          else if ((veldir.eq.3).and.(SDIM.eq.3)) then
  
           if (adv_dir.eq.3) then
            z_vel = adv_vel
           else  
            z_vel = zero
           endif
           if ((probtype.eq.53).or.(probtype.eq.530).or. &
               (probtype.eq.531).or. &
               (probtype.eq.536).or.(probtype.eq.537).or. &
               (probtype.eq.538).or.(probtype.eq.541)) then
            z_vel=zero
           endif
           velcell(veldir)=z_vel  ! zhi, velz
  
           if (probtype.eq.bubbleInPackedColumn) then ! zhi,velz
            call Pack_velbc(dir,side,veldir,velcell)
           else if ((probtype.eq.62).or.(probtype.eq.65)) then
            costheta=cos(xblob2)
            sintheta=sin(xblob2)
            xprime=costheta*(x-xblob)-sintheta*(z-zblob)
            yprime=y-yblob
            zprime=sintheta*(x-xblob)+costheta*(z-zblob)
            if (xprime**2+yprime**2<radblob**2) then
             velcell(veldir)=-radblob3
            endif
           else if ((probtype.eq.63).or.(probtype.eq.64)) then
            velcell(veldir)=four*z_vel
! velz, zhi
           else if ((probtype.eq.36).and.(xblob10.gt.zero).and. &
                    ((yblob9.ne.zero).or.(yblob10.ne.zero))) then
            if (probhix-problox.le.zero) then
             print *,"probhix or problox invalid"
             stop
            endif
            velcell(veldir)=yblob9+(x-problox)*(yblob10-yblob9)/ &
               (probhix-problox)

           endif
  
          else
           print *,"veldir invalid"
           stop
          endif
  
         else
          print *,"dir invalid velbc_override 3"
          stop
         endif
  
        else 
         print *,"side invalid"
         stop
        endif
  
       endif ! non-hydrate problems

       vel=velcell(veldir)

       if (global_velocity_scale.gt.zero) then
        vel_in=vel/global_velocity_scale
       else
        print *,"globa_velocity_scale invalid"
        stop
       endif

      else
       print *,"vel_homflag invalid"
       stop
      endif

      return
      end subroutine velbc_override


 
      subroutine scalarBC(time,dir,side,ADV,ADVwall, &
       xsten,nhalf,dx,bfact)
      IMPLICIT NONE

      integer, INTENT(in) :: dir,side,nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: ADV
      real(amrex_real) :: xwall
      real(amrex_real), INTENT(inout) :: ADVwall
      real(amrex_real) :: x,y,z
      real(amrex_real), INTENT(in) :: dx(SDIM)

      if (nhalf.lt.1) then
       print *,"nhalf invalid scalar bc"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else
        xwall=probhix
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else
        xwall=probhiy
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else
        xwall=probhiz
       endif
      else
       print *,"dir invalid scalarbc"
       stop
      endif

      if ((dir.eq.1).and.(side.eq.1)) then
       if (xwall.lt.x) then
        print *,"xwall,x invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.1).and.(side.eq.2)) then
       if (xwall.gt.x) then
        print *,"xwall,x invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.2).and.(side.eq.1)) then
       if (xwall.lt.y) then
        print *,"xwall,y invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.2).and.(side.eq.2)) then
       if (xwall.gt.y) then
        print *,"xwall,y invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then
       if (xwall.lt.z) then
        print *,"xwall,z invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then
       if (xwall.gt.z) then
        print *,"xwall,z invalid"
        stop
       endif
       ADV=ADVwall
      else
       print *,"dir side invalid"
       stop
      endif

      return
      end subroutine scalarBC

       ! sets all the physical BCs to the interior value.
      subroutine extrapBC(time,dir,side,ADV,ADVwall, &
       xsten,nhalf,dx,bfact)
      IMPLICIT NONE

      integer, INTENT(in) :: dir,side,nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: ADV
      real(amrex_real) :: xwall
      real(amrex_real), INTENT(inout) :: ADVwall
      real(amrex_real) :: x,y,z
      real(amrex_real), INTENT(in) :: dx(SDIM)

      if (nhalf.lt.1) then
       print *,"nhalf invalid extrap bc"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else if (side.eq.2) then
        xwall=probhix
       else
        print *,"side invalid"
        stop
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else if (side.eq.2) then
        xwall=probhiy
       else
        print *,"side invalid"
        stop
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else if (side.eq.2) then
        xwall=probhiz
       else
        print *,"side invalid"
        stop
       endif
      else
       print *,"dir invalid scalarbc"
       stop
      endif

      if ((dir.eq.1).and.(side.eq.1)) then
       if (xwall+dx(dir)*VOFTOL.ge.x) then
        ! do nothing
       else
        print *,"xwall,x invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.1).and.(side.eq.2)) then
       if (xwall-dx(dir)*VOFTOL.le.x) then
        ! do nothing
       else
        print *,"xwall,x invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.2).and.(side.eq.1)) then
       if (xwall+dx(dir)*VOFTOL.ge.y) then
        ! do nothing
       else
        print *,"xwall,y invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.2).and.(side.eq.2)) then
       if (xwall-dx(dir)*VOFTOL.le.y) then
        ! do nothing
       else
        print *,"xwall,y invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then
       if (xwall+dx(dir)*VOFTOL.ge.z) then
        ! do nothing
       else
        print *,"xwall,z invalid"
        stop
       endif
       ADV=ADVwall
      else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then
       if (xwall-dx(dir)*VOFTOL.le.z) then
        ! do nothing
       else
        print *,"xwall,z invalid"
        stop
       endif
       ADV=ADVwall
      else
       print *,"dir side invalid"
       stop
      endif

      return
      end subroutine extrapBC

       !called from fort_pressurefill
      subroutine presBDRYCOND(time,dir,side,ADV,ADVwall_in, &
        xsten,nhalf,dx,bfact)
      use global_utility_module
      use hydrateReactor_module
      use unimaterialChannel_module
      use River
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use rigid_FSI_module
      use sinking_particle_module

      IMPLICIT NONE


      integer, INTENT(in) :: dir,side,nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: ADV
      real(amrex_real) :: xwall
      real(amrex_real) :: ADVwall
      real(amrex_real) :: x,y,z
      real(amrex_real), INTENT(inout) :: ADVwall_in
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) pipexlo,pipexhi
      real(amrex_real) rhohydro
      real(amrex_real) base_pres
      real(amrex_real) VOF(num_materials)
      real(amrex_real) cenbc(num_materials,SDIM)
      real(amrex_real) waterdepth
      real(amrex_real) xpos(SDIM)
      integer local_dir
      real(amrex_real) local_LS(num_materials)
      integer, PARAMETER :: from_boundary_hydrostatic=0

      integer :: gravity_dir

      call fort_derive_gravity_dir(gravity_vector,gravity_dir)

      if (nhalf.lt.1) then
       print *,"nhalf invalid presBDRYCOND"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)
      do local_dir=1,SDIM
       xpos(local_dir)=xsten(0,local_dir)
      enddo

      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else
        xwall=probhix
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else
        xwall=probhiy
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else
        xwall=probhiz
       endif
      else
       print *,"dir invalid presbc"
       stop
      endif

      if (SDIM.eq.2) then
       if (abs(z-y).le.EPS2) then
        ! do nothing
       else
        print *,"z<>y in 17219"
        stop
       endif
      endif

      ADVwall=ADVwall_in*global_pressure_scale

      waterdepth=zero

      if (pres_homflag.eq.1) then
       ADV=zero
      else if (pres_homflag.eq.0) then

       if (is_in_probtype_list().eq.1) then
        call SUB_LS(xpos,time,local_LS,num_materials)
        call SUB_PRES_BC(xwall,xpos,time,local_LS, &
         ADV,ADVwall,dir,side,dx,num_materials)

       else if (probtype.eq.401) then

        call HELIX_LS(xpos,time,local_LS)
        call HELIX_PRES_BC(xwall,xpos,time,local_LS, &
          ADV,ADVwall,dir,side,dx)
       else if (probtype.eq.402) then

        call TSPRAY_LS(xpos,time,local_LS)
        call TSPRAY_PRES_BC(xwall,xpos,time,local_LS, &
          ADV,ADVwall,dir,side,dx)
       else if (probtype.eq.412) then ! step
        call CAV2Dstep_LS(xpos,time,local_LS)
        call CAV2Dstep_PRES_BC(xwall,xpos,time,local_LS, &
          ADV,ADVwall,dir,side,dx)
       else if (probtype.eq.413) then ! zeyu
        call ZEYU_droplet_impact_LS(xpos,time,local_LS)
        call ZEYU_droplet_impact_PRES_BC(xwall,xpos,time,local_LS, &
          ADV,ADVwall,dir,side,dx)

       else if (probtype.eq.533) then
        call rigid_FSI_LS(xpos,time,local_LS)
        call rigid_FSI_PRES_BC(xwall,xpos,time,local_LS, &
          ADV,ADVwall,dir,side,dx)
       else if (probtype.eq.534) then
        call sinking_FSI_LS(xpos,time,local_LS)
        call sinking_FSI_PRES_BC(xwall,xpos,time,local_LS, &
          ADV,ADVwall,dir,side,dx)

       else if (probtype.eq.311) then ! user defined

        call USERDEF_LS(xpos,time,local_LS)
        call USERDEF_PRES_BC(xwall,xpos,time,local_LS, &
          ADV,ADVwall,dir,side,dx)

       else if (probtype.eq.199) then ! hydrates

        call HYD_PRES_BC(time,dir,side,ADV,xwall,ADVwall, &
          x,y,z,dx)

        ! in: presBDRYCOND
       else if (probtype.eq.220) then ! UNIMATERIAL problem
        call UNIMAT_PRES_BC(time,dir,side,ADV,xwall,ADVwall, &
          x,y,z,dx)

       else if ((probtype.eq.299).or. &
                (probtype.eq.301)) then !melting (boundary condition pressure)
        ADV=zero
       else if (probtype.eq.209) then ! River
        call RiverPressure(x,y,z,time,ADV, &
         fort_denconst(2),fort_denconst(1),axis_dir)
       else  
 
        base_pres=zero
        if ((probtype.eq.53).and.(axis_dir.eq.2)) then ! injection shock
         call general_hydrostatic_pressure(base_pres)
        else if ((probtype.eq.53).and. &  ! compressible JIC
                 (fort_material_type(2).gt.0)) then
         call general_hydrostatic_pressure(base_pres)
        else if ((probtype.eq.538).and.(SDIM.eq.2)) then 
         ! diesel injection w/needle
!        call general_hydrostatic_pressure(base_pres)
          base_pres=outflow_pressure
        else if ((probtype.eq.541).and.(SDIM.eq.2)) then
          base_pres=outflow_pressure
        else if (probtype.eq.539) then  ! sup injector
         call general_hydrostatic_pressure(base_pres)
        else if ((probtype.eq.530).and. &  ! impinging jets
                 (axis_dir.eq.1).and. &
                 (fort_material_type(2).gt.0).and. &
                 (SDIM.eq.3)) then
         call general_hydrostatic_pressure(base_pres)
        else if ((probtype.eq.538).and.(SDIM.eq.3)) then 
          ! diesel injection w/needle
!        call general_hydrostatic_pressure(base_pres)
          base_pres=outflow_pressure
        else if ((probtype.eq.541).and.(SDIM.eq.3)) then
          base_pres=outflow_pressure

         ! CODY ESTEBE created test problem
        else if ((probtype.eq.46).and.(axis_dir.eq.20)) then
          base_pres=outflow_pressure
        endif

        if ((dir.eq.1).and.(side.eq.1)) then  ! xlo
         if (xwall.lt.x) then
          print *,"xwall,x invalid"
          stop
         endif
         ADV=base_pres

          ! xlo,presBDRYCOND, stratified fluid
         if (probtype.eq.201) then

          if (z.le.zblob2) then
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-zblob2)
          else
           ADV=-fort_denconst(3)*abs(gravity_vector(gravity_dir))*(z-zblob2)
          endif
  
         else if (SDIM.eq.2) then

          if ((probtype.eq.92).or. & !contact captured
              (probtype.eq.93)) then !contact tracked
           if (axis_dir.eq.2) then ! shock-turbulence pressure xlo
            ADV=10.333333
           else if (axis_dir.eq.3) then ! mach>4
            ADV=10.0*(1.4-1.0)
           else if ((axis_dir.eq.0).or. &
                    (axis_dir.eq.1).or. &
                    (axis_dir.eq.4)) then
            ! do nothing
           else
            print *,"axis_dir invalid"
            stop
           endif
          else if (fort_material_type(1).eq.13) then ! xlo, 2D, Tait EOS
           call boundary_hydrostatic(xpos,rhohydro,ADV)

          ! xlo, presBDRYCOND, 2D
          else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  

           if (y.le.waterdepth) then
            ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(y-waterdepth)
           else
            ADV=-fort_denconst(2)*abs(gravity_vector(gravity_dir))*(y-waterdepth)
           endif

          else if (probtype.eq.539) then  ! supnozzle xlo - presBDRYCOND
           ADV=inflow_pressure ! inlet side is all gas

          endif

         else if (SDIM.eq.3) then

          if (fort_material_type(1).eq.13) then !xlo 3D(1st material Tait EOS)
           call boundary_hydrostatic(xpos,rhohydro,ADV)
          else if ((probtype.eq.36).and.(axis_dir.eq.0)) then ! xlo
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*z

           ! xlo, presBDRYCOND, 3D
          else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  
           if (z.le.waterdepth) then
            ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
           else
            ADV=-fort_denconst(2)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
           endif
          endif

         else
          print *,"dimension bust"
          stop
         endif

        else if ((dir.eq.1).and.(side.eq.2)) then  ! xhi, presBDRYCOND
         if (xwall.gt.x) then
          print *,"xwall,x invalid"
          stop
         endif
         ADV=base_pres
   
          ! xhi, presBDRYCOND, stratified
         if (probtype.eq.201) then

          if (z.le.zblob2) then
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-zblob2)
          else
           ADV=-fort_denconst(3)*abs(gravity_vector(gravity_dir))*(z-zblob2)
          endif

         else if (SDIM.eq.2) then

          if ((probtype.eq.92).or. & !contact captured
              (probtype.eq.93)) then !contact tracked
           if (axis_dir.eq.2) then ! shock-turbulence pressure xhi
            ADV=one
           else if (axis_dir.eq.3) then ! mach>4
            ADV=(1.4-1.0)
           else if ((axis_dir.eq.0).or. &
                    (axis_dir.eq.1).or. &
                    (axis_dir.eq.4)) then
            ! do nothing
           else
            print *,"axis_dir invalid"
            stop
           endif
          else if (fort_material_type(1).eq.13) then ! xhi, 2D, Tait EOS
           call boundary_hydrostatic(xpos,rhohydro,ADV)
          else if ((probtype.eq.36).and.(axis_dir.eq.0)) then ! xhi
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*y
          else if ((probtype.eq.36).and.(axis_dir.eq.2)) then
            !calling from presBDRYCOND
           call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                   from_boundary_hydrostatic)
          else if (probtype.eq.42) then ! bubble jetting 2D
            !calling from presBDRYCOND
           call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                   from_boundary_hydrostatic)
            ! in presBDRYCOND: 2d, xhi
          else if (probtype.eq.46) then ! cavitation 2D
           if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
             !calling from presBDRYCOND
            call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                    from_boundary_hydrostatic)
           else if (axis_dir.eq.10) then
             !calling from presBDRYCOND
            call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                    from_boundary_hydrostatic)
           else if (axis_dir.eq.20) then
            ! do nothing (ADV=base_pres, base_pres=outflow_pressure above) 
           else
            print *,"axis_dir invalid"
            stop
           endif

           ! xhi, presBDRYCOND, 2D
          else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  
           if (y.le.waterdepth) then
            ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(y-waterdepth)
           else
            ADV=-fort_denconst(2)*abs(gravity_vector(gravity_dir))*(y-waterdepth)
           endif
          else if (probtype.eq.539) then  ! supnozzle xhi - presBDRYCOND
           ADV=outflow_pressure

          endif

         else if (SDIM.eq.3) then

          if (fort_material_type(1).eq.13) then ! xhi, presBDRYCOND, TAIT EOS
           call boundary_hydrostatic(xpos,rhohydro,ADV)
          else if ((probtype.eq.36).and.(axis_dir.eq.0)) then ! xhi
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*z

           ! xhi, presBDRYCOND, 3D
          else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  
           if (z.le.waterdepth) then
            ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
           else
            ADV=-fort_denconst(2)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
           endif
          endif

         else
          print *,"dimension bust"
          stop
         endif

        else if ((dir.eq.2).and.(side.eq.1)) then  ! ylo

         if (xwall.lt.y) then
          print *,"xwall,y invalid"
          stop
         endif
         ADV=base_pres

          ! ylo,presBDRYCOND, stratified fluid, material 1 is bottom
         if (probtype.eq.201) then

          if (z.le.zblob2) then
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-zblob2)
          else
           ADV=-fort_denconst(3)*abs(gravity_vector(gravity_dir))*(z-zblob2)
          endif

         else if (SDIM.eq.2) then

          if (fort_material_type(1).eq.13) then  !ylo 2D (TAIT EOS)
           call boundary_hydrostatic(xpos,rhohydro,ADV)
          else if ((probtype.eq.36).and.(axis_dir.eq.0)) then ! ylo
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*y
          else if ((probtype.eq.36).and.(axis_dir.eq.2)) then
             !calling from presBDRYCOND
           call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                   from_boundary_hydrostatic)
          else if (probtype.eq.42) then ! bubble jetting 2D
             !calling from presBDRYCOND
           call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                   from_boundary_hydrostatic)
          else if (probtype.eq.46) then ! cavitation ylo 2D
           if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
             !calling from presBDRYCOND
            call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                    from_boundary_hydrostatic)
           else if (axis_dir.eq.10) then
             !calling from presBDRYCOND
            call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                    from_boundary_hydrostatic)
           else if (axis_dir.eq.20) then
            ! do nothing (ADV=base_pres, base_pres=outflow_pressure above) 
           else
            print *,"axis_dir invalid"
            stop
           endif

           ! ylo, presBDRYCOND, 2D
          else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  
           if (y.le.waterdepth) then
            ADV= &
             -fort_denconst(1)*abs(gravity_vector(gravity_dir))*(y-waterdepth)
           else
            ADV=-fort_denconst(2)*abs(gravity_vector(gravity_dir))*(y-waterdepth)
           endif
          endif

           ! nozzle + pressure bc
          if ((probtype.eq.53).and.(axis_dir.eq.2)) then ! injection+shock
           call get_jet_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc) 
           if (VOF(1).gt.zero) then
            ADV=inflow_pressure
           endif 
           ! 2d diesel injection w/needle ylo presBDRYCOND
          else if (probtype.eq.538) then 
           call get_initial_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc)
           if (VOF(1).gt.zero) then
            ADV=inflow_pressure
           endif
           ! 2d diesel injection w/needle ylo presBDRYCOND
          else if ((probtype.eq.541).and.(axis_dir.eq.2)) then
           call get_initial_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc)
           if (VOF(1).gt.zero) then
            ADV=inflow_pressure
           endif
          else if (probtype.eq.539) then  ! supnozzle ylo - presBDRYCOND
           ADV=outflow_pressure
           call get_jet_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc) 
           if (VOF(1).gt.zero) then
! Hack: should be input, but "inflow_pressure" used at xlo
            ADV=1.6e7
           endif
          endif

         else if (SDIM.eq.3) then

          if (fort_material_type(1).eq.13) then ! ylo 3D, Tait EOS
           call boundary_hydrostatic(xpos,rhohydro,ADV)
          else if ((probtype.eq.36).and.(axis_dir.eq.0)) then  ! ylo
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*z

           ! ylo, presBDRYCOND, 3D
          else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  
           if (z.le.waterdepth) then
            ADV= &
             -fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
           else
            ADV= &
             -fort_denconst(2)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
           endif
          endif

         else
          print *,"dimension bust"
          stop
         endif

        else if ((dir.eq.2).and.(side.eq.2)) then  ! yhi
         if (xwall.gt.y) then
          print *,"xwall,y invalid"
          stop
         endif
         ADV=base_pres

          ! yhi,presBDRYCOND, stratified
         if (probtype.eq.201) then

          if (z.le.zblob2) then
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-zblob2)
          else
           ADV=-fort_denconst(3)*abs(gravity_vector(gravity_dir))*(z-zblob2)
          endif

         else if (SDIM.eq.2) then

          if (fort_material_type(1).eq.13) then ! yhi,2D, Tait EOS
           call boundary_hydrostatic(xpos,rhohydro,ADV)
          else if ((probtype.eq.36).and.(axis_dir.eq.0)) then
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*y
          else if ((probtype.eq.36).and.(axis_dir.eq.2)) then  ! yhi
            !calling from presBDRYCOND
           call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                   from_boundary_hydrostatic)
           ! yhi, presBDRYCOND, 2D
          else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  
           if (y.le.waterdepth) then
            ADV= &
             -fort_denconst(1)*abs(gravity_vector(gravity_dir))*(y-waterdepth)
           else
            ADV= &
             -fort_denconst(2)*abs(gravity_vector(gravity_dir))*(y-waterdepth)
           endif
          else if (probtype.eq.42) then ! bubble jetting 2D
            !calling from presBDRYCOND
           call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                   from_boundary_hydrostatic)
          else if (probtype.eq.46) then ! cavitation yhi 2D
           if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
             !calling from presBDRYCOND
            call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                    from_boundary_hydrostatic)
           else if (axis_dir.eq.10) then
             !calling from presBDRYCOND
            call tait_hydrostatic_pressure_density(xpos,rhohydro,ADV, &
                    from_boundary_hydrostatic)
           else if (axis_dir.eq.20) then
            ! do nothing (ADV=base_pres, base_pres=outflow_pressure above) 
           else
            print *,"axis_dir invalid"
            stop
           endif

          endif

          if (probtype.eq.41) then ! presBDRYCOND 2D yhi
           pipexlo=problox
           pipexhi=probhix
           if ((axis_dir.eq.1).or.(axis_dir.eq.2).or.(axis_dir.eq.3)) then
            pipexlo=zero
            pipexhi=two*radblob3
           endif

           if (axis_dir.eq.0) then
            ! do nothing
           else if (axis_dir.eq.1) then

            if (x.lt.pipexlo) then
             ADV=zero
            else if (x.lt.xblob) then
             ADV=fort_denconst(2)*abs(gravity_vector(gravity_dir))*(x-pipexlo)
            else
             ADV= &
              fort_denconst(2)*abs(gravity_vector(gravity_dir))* &
              (xblob-pipexlo)+ &
              fort_denconst(1)*abs(gravity_vector(gravity_dir))*(x-xblob) 
            endif

           else if (axis_dir.eq.2) then

            if (x.lt.pipexlo) then
             ADV=zero
            else if (x.lt.xblob-radblob2) then
             ADV=fort_denconst(1)*abs(gravity_vector(gravity_dir))*(x-pipexlo)
            else if (x.lt.xblob+radblob2) then
             ADV=fort_denconst(1)*abs(gravity_vector(gravity_dir))*(xblob-radblob2-pipexlo)+ &
                 fort_denconst(2)*abs(gravity_vector(gravity_dir))*(x-xblob+radblob2) 
            else
             ADV=fort_denconst(1)*abs(gravity_vector(gravity_dir))*(xblob-radblob2-pipexlo)+ &
                 fort_denconst(2)*abs(gravity_vector(gravity_dir))*(two*radblob2)+ & 
                 fort_denconst(1)*abs(gravity_vector(gravity_dir))*(x-xblob-radblob2) 
            endif

           else if (axis_dir.eq.3) then

            if (x.lt.pipexlo) then
             ADV=zero
            else 
             ADV=fort_denconst(1)*abs(gravity_vector(gravity_dir))*(x-pipexlo)
            endif

           else if (axis_dir.eq.4) then
            ! do nothing 
           else
            print *,"axis_dir invalid for pipe problem"
            stop
           endif  
  
           ! above, probtype=41 below supnozzle yhi - presBDRYCOND
          else if (probtype.eq.539) then
           ADV=outflow_pressure
          endif  

         else if (SDIM.eq.3) then

          if (fort_material_type(1).eq.13) then ! yhi, 3D, Tait EOS
           call boundary_hydrostatic(xpos,rhohydro,ADV)
          else if ((probtype.eq.36).and.(axis_dir.eq.0)) then ! yhi
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*z
           ! yhi, presBDRYCOND, 3D
          else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  
           if (z.le.waterdepth) then
            ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
           else
            ADV=-fort_denconst(2)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
           endif
          endif

         else 
          print *,"dimension bust"
          stop
         endif

        else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then  ! zlo
         if (xwall.lt.z) then
          print *,"xwall,z invalid"
          stop
         endif
         ADV=base_pres
         if (fort_material_type(1).eq.13) then ! zlo, 3D, Tait EOS
          call boundary_hydrostatic(xpos,rhohydro,ADV)
         else if (probtype.eq.201) then ! zlo stratified material 1 at bottom
          ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-zblob2)
         else if ((probtype.eq.36).and.(axis_dir.eq.0)) then  ! zlo
          ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*z
         else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  ! zlo, presBDRYCOND
          if (z.le.waterdepth) then
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
          else
           ADV=-fort_denconst(2)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
          endif
         endif

          ! pressure bc at inflow
         if ((probtype.eq.53).and.(axis_dir.eq.2)) then ! injection shock
          call get_jet_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc)   
          if (VOF(1).gt.zero) then
           ADV=inflow_pressure
          endif 
           ! 3D zlo presBDRYCOND diesel injection w/needle
         else if (probtype.eq.538) then
          call get_initial_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc)
          if (VOF(1).gt.zero) then
           ADV=inflow_pressure
          endif
           ! 3D zlo presBDRYCOND diesel injection w/needle
         else if ((probtype.eq.541).and.(axis_dir.eq.2)) then
          call get_initial_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc)
          if (VOF(1).gt.zero) then
           ADV=inflow_pressure
          endif
         endif

        else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then  ! zhi
         if (xwall.gt.z) then
          print *,"xwall,z invalid"
          stop
         endif
         ADV=base_pres
         if (fort_material_type(1).eq.13) then ! zhi, 3D, Tait EOS
          call boundary_hydrostatic(xpos,rhohydro,ADV)
         else if (probtype.eq.201) then ! zhi, stratified
          ADV=-fort_denconst(3)*abs(gravity_vector(gravity_dir))*(z-zblob2)
         else if ((probtype.eq.36).and.(axis_dir.eq.0)) then ! zhi
          ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*z
         else if ((probtype.eq.9).and.(axis_dir.eq.1)) then  ! zhi, presBDRYCOND
          if (z.le.waterdepth) then
           ADV=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
          else
           ADV=-fort_denconst(2)*abs(gravity_vector(gravity_dir))*(z-waterdepth)
          endif
         endif

        else
         print *,"dir side invalid"
         stop
        endif

       endif ! non-hydrate problems

      else
       print *,"pres_homflag invalid"
       stop
      endif

      ADV=ADV/global_pressure_scale

      return
      end subroutine presBDRYCOND

      ! called from FORT_STATEFILL, FORT_GROUP_STATEFILL if EXT_DIR
      subroutine denBC(time,dir,side,ADV,ADVwall_in, &
         xsten,nhalf,dx,bfact,istate,im)
      use global_utility_module
      use hydrateReactor_module
      use unimaterialChannel_module
      use marangoni
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use rigid_FSI_module
      use sinking_particle_module

      IMPLICIT NONE

      integer, INTENT(in) :: dir,side,istate,im,nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: ADV
      real(amrex_real) xwall
      real(amrex_real) ADVwall
      real(amrex_real) x,y,z
      real(amrex_real), INTENT(inout) :: ADVwall_in
      real(amrex_real) ADV_merge
      real(amrex_real) dx(SDIM)
      real(amrex_real) dist,water_temp
      integer species_flag,local_homflag
      real(amrex_real) xvec(SDIM)
      integer local_dir
      real(amrex_real) local_LS(num_materials)
      integer try_merge
      real(amrex_real) :: massfrac_parm(num_species_var+1)

      if (nhalf.lt.1) then
       print *,"nhalf invalid denbc"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((im.le.0).or.(im.gt.num_materials)) then
       print *,"im invalid80: ",im
       stop
      endif
      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)
      do local_dir=1,SDIM
       xvec(local_dir)=xsten(0,local_dir)
      enddo

      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else
        xwall=probhix
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else
        xwall=probhiy
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else
        xwall=probhiz
       endif
      else
       print *,"dir invalid denbc"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid in denBC"
       stop
      endif

      if ((istate.ge.1).and. &
          (istate.le.num_state_base)) then
       species_flag=0
       if (istate.eq.ENUM_DENVAR+1) then ! density
        try_merge=0
       else if (istate.eq.ENUM_TEMPERATUREVAR+1) then ! temperature
        try_merge=1
       else
        print *,"istate invalid"
        stop
       endif
       if ((temp_homflag.eq.1).and. &
           (istate.eq.2)) then
        local_homflag=1
       else
        local_homflag=0
       endif
      else if ((istate.gt.num_state_base).and. & ! species
               (istate.le.num_state_material)) then
       species_flag=1
       try_merge=1
       if (species_homflag.eq.1) then
        local_homflag=1
       else
        local_homflag=0
       endif
      else
       print *,"istate invalid"
       stop
      endif

        ! homogeneous BC for temperature solver
      if ((temp_homflag.eq.1).and. &
          (species_flag.eq.0).and. &
          (istate.eq.2).and. &  ! temperature
          (local_homflag.eq.1)) then

       ADV=zero
       ADV_merge=zero

        ! homogeneous BC for species solver
      else if ((species_homflag.eq.1).and. &
               (species_flag.eq.1).and. &
               (istate.gt.num_state_base).and. &
               (istate.le.num_state_material).and. &
               (local_homflag.eq.1)) then

       ADV=zero
       ADV_merge=zero

      else if (local_homflag.eq.0) then

       if (istate.eq.ENUM_DENVAR+1) then
        ADVwall=ADVwall_in   ! den
       else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
        ADVwall=ADVwall_in   ! T
       else if ((istate.gt.2).and. &
                (istate.le.num_state_material)) then
        ADVwall=ADVwall_in   ! species
       else
        print *,"istate invalid"
        stop
       endif

       if ((im.lt.1).or.(im.gt.num_materials)) then
        print *,"im invalid79"
        stop
       endif

       if (is_in_probtype_list().eq.1) then

        call SUB_LS(xvec,time,local_LS,num_materials)
        call SUB_STATE_BC(xwall,xvec,time,local_LS, &
         ADV,ADV_merge,ADVwall,im,istate,dir,side,dx,num_materials)

       else if (probtype.eq.401) then

        call HELIX_LS(xvec,time,local_LS)
        call HELIX_STATE_BC(xwall,xvec,time,local_LS, &
          ADV,ADV_merge,ADVwall,im,istate,dir,side,dx) 

       else if (probtype.eq.402) then

        call TSPRAY_LS(xvec,time,local_LS)
        call TSPRAY_STATE_BC(xwall,xvec,time,local_LS, &
          ADV,ADV_merge,ADVwall,im,istate,dir,side,dx) 

       else if (probtype.eq.412) then ! step

        call CAV2Dstep_LS(xvec,time,local_LS)
        call CAV2Dstep_STATE_BC(xwall,xvec,time,local_LS, &
          ADV,ADV_merge,ADVwall,im,istate,dir,side,dx)

       else if (probtype.eq.413) then ! zeyu

        call ZEYU_droplet_impact_LS(xvec,time,local_LS)
        call ZEYU_droplet_impact_STATE_BC(xwall,xvec,time,local_LS, &
          ADV,ADV_merge,ADVwall,im,istate,dir,side,dx)

       else if (probtype.eq.533) then

        call rigid_FSI_LS(xvec,time,local_LS)
        call rigid_FSI_STATE_BC(xwall,xvec,time,local_LS, &
          ADV,ADV_merge,ADVwall,im,istate,dir,side,dx) 
       else if (probtype.eq.534) then

        call sinking_FSI_LS(xvec,time,local_LS)
        call sinking_FSI_STATE_BC(xwall,xvec,time,local_LS, &
          ADV,ADV_merge,ADVwall,im,istate,dir,side,dx) 

       else if (probtype.eq.311) then ! user defined

        call USERDEF_LS(xvec,time,local_LS)
        call USERDEF_STATE_BC(xwall,xvec,time,local_LS, &
          ADV,ADV_merge,ADVwall,im,istate,dir,side,dx) 

       else if (species_flag.eq.0) then

         ! in: subroutine denBC
         ! HYDRATES
        if (probtype.eq.199) then

         if (istate.eq.ENUM_DENVAR+1) then
          call HYD_DENS_BC(time,dir,side,ADV,xwall,ADVwall, &
            x,y,z,dx,im)
         else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
          call HYD_TEMP_BC(time,dir,side,ADV,xwall,ADVwall, &
            x,y,z,dx,im)
         else
          print *,"istate invalid"
          stop
         endif

         ! in: subroutine denBC
         ! UNIMATERIAL problem
        else if (probtype.eq.220) then

         if (istate.eq.ENUM_DENVAR+1) then
          call UNIMAT_DENS_BC(time,dir,side,ADV,xwall,ADVwall, &
            x,y,z,dx,im)
         else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
          call UNIMAT_TEMP_BC(time,dir,side,ADV,xwall,ADVwall, &
            x,y,z,dx,im)
         else
          print *,"istate invalid"
          stop
         endif

        else if ((probtype.eq.299).or. &
                 (probtype.eq.301)) then !melting

         if (istate.eq.ENUM_DENVAR+1) then ! density bdry cond
          ADV=ADVwall
         else if (istate.eq.ENUM_TEMPERATUREVAR+1) then ! temperature bdry cond
          ADV=fort_tempconst(2)  ! gas temperature
         else
          print *,"istate invalid"
          stop
         endif

         ! marangoni (heat pipe) problem
        else if ((probtype.eq.36).and.(axis_dir.eq.10)) then

         if (istate.eq.ENUM_DENVAR+1) then
          print *,"density bc should be FOEXTRAP"
          stop
         else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
           ! flag=0 (no buffer, do not start w/room temp in the bulk)
           ! second parameter would be the lateral buffer size. (delta)
           ! third parameter would be the front/back buffer size. (delta2)
          call position_Temp(0,radblob,radblob2,x,y,z,ADV)
         else 
          print *,"istate invalid"
          stop
         endif
  
        else

          ! xlo
         if ((dir.eq.1).and.(side.eq.1)) then
          if (xwall.lt.x) then
           print *,"xwall,x invalid"
           stop
          endif
          ADV=ADVwall

          if (SDIM.eq.2) then

            ! xlo states (freezing singularity problem, boiling problem)
           if (probtype.eq.59) then

            if (istate.eq.ENUM_DENVAR+1) then
             ! do nothing (density)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ! bcflag=1 (calling from denBC - boundary conditions
             ! for density, temperature and species variables)
             call outside_temperature(time,x,y,z,ADV,im,1) 
            else
             print *,"istate invalid"
             stop
            endif

            ! xlo states
            ! denBC: compare w/Welch and Wilson 2000
           else if (probtype.eq.801) then 

            if (num_materials.lt.2) then
             print *,"num_materials invalid probtype=801"
             stop
            endif
             ! temperature (probtype==801)
            if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             call vapordist(xsten,nhalf,dx,bfact,dist)
             if (dist.lt.zero) then
              ADV=fort_tempconst(2) ! vapor/ice temperature
             else 
              ADV=fort_tempconst(1) ! liquid temperature 
             endif
            endif
           else if (probtype.eq.802) then ! dissolution, xlo
            print *,"802 obsolete"
            stop
           else if (probtype.eq.532) then
            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
            else
             print *,"istate invalid"
             stop
            endif 
           else if (probtype.eq.539) then ! xlo, sup, denBC

            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
! Hack: inlet temperature (should not be called since outflow BC?)
             ADV=513.0
            else
             print *,"istate invalid"
             stop
            endif 

           else if (probtype.eq.53) then ! xlo, denBC

            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
            else
             print *,"istate invalid"
             stop
            endif 
           endif  ! probtype=53

          else if (SDIM.eq.3) then

           if (probtype.eq.53) then
            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
            else
             print *,"istate invalid"
             stop
            endif
           endif  ! probtype=53

          else
           print *,"dimension bust"
           stop
          endif

          ! xhi
         else if ((dir.eq.1).and.(side.eq.2)) then
          if (xwall.gt.x) then
           print *,"xwall,x invalid"
           stop
          endif
          ADV=ADVwall

          if (SDIM.eq.2) then

            ! xhi states (freezing singularity problem, boiling problem)
           if (probtype.eq.59) then

            if (istate.eq.ENUM_DENVAR+1) then
             ! do nothing (density)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ! bcflag=1 (calling from denBC)
             call outside_temperature(time,x,y,z,ADV,im,1) 
            else
             print *,"istate invalid"
             stop
            endif

            ! xhi states
           else if (probtype.eq.801) then !denBC:compare w/Welch and Wilson 2000
            if (num_materials.lt.2) then
             print *,"num_materials invalid probtype=801"
             stop
            endif
            if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             call vapordist(xsten,nhalf,dx,bfact,dist)
             if (dist.lt.zero) then
              ADV=fort_tempconst(2) !vapor/ice temperature
             else 
              ADV=fort_tempconst(1) !liquid temperature (air temp in far field)
             endif
            endif
           else if (probtype.eq.532) then
            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
            else
             print *,"istate invalid"
             stop
            endif 
           else if (probtype.eq.539) then  ! xhi, sup, denBC

            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
            else
             print *,"istate invalid"
             stop
            endif 
           endif

          else if (SDIM.eq.3) then

          else
           print *,"dimension bust"
           stop
          endif

          ! ylo (subroutine denBC)
         else if ((dir.eq.2).and.(side.eq.1)) then
          if (xwall.lt.y) then
           print *,"xwall,y invalid"
           stop
          endif
          ADV=ADVwall

            ! ylo states

          if (SDIM.eq.2) then

            ! melting ice block on substrate (subroutine denBC) 
            ! ylo
            ! prescribe_temperature_outflow:
            ! 0=dirichlet at inflow
            ! 1=dirichlet at inflow and outflow
            ! 2=dirichlet at inflow and walls.
            ! 3=dirichlet at inflow, outflow, and walls.
           if ((prescribe_temperature_outflow.eq.3).and. &
               (probtype.eq.59)) then

            if (istate.eq.ENUM_DENVAR+1) then
             ! do nothing (density)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ! bcflag=1 (calling from denBC)
             call outside_temperature(time,x,y,z,ADV,im,1) 
            else
             print *,"istate invalid"
             stop
            endif
  
            ! denBC: compare w/Welch and Wilson 2000
           else if (probtype.eq.801) then 

            if (num_materials.lt.2) then
             print *,"num_materials invalid probtype=801"
             stop
            endif
             ! ylo
            if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             call vapordist(xsten,nhalf,dx,bfact,dist)
             if (dist.lt.zero) then
              ADV=fort_tempconst(2) ! vapor/ice temperature
             else 
              ADV=fort_tempconst(1) ! liquid temperature + far field
             endif
            endif
           else if ((probtype.eq.538).or.(probtype.eq.53).or. &
                    (probtype.eq.541)) then
            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
            else
             print *,"istate invalid"
             stop
            endif 
            ! above denBC - probtype=53 or 538 ylo face
            ! below Benard convection
           else if (probtype.eq.603) then
            water_temp=radblob2+fort_tempconst(1)
            call init_massfrac_parm(fort_denconst(1),massfrac_parm,1)
            call INTERNAL_material(fort_denconst(1),massfrac_parm, &
             water_temp, &
             fort_energyconst(1),fort_material_type(1),1)
            if (fort_energyconst(1).gt.zero) then
             ! do nothing
            else
             print *,"fort_energyconst(1) invalid in denBC ylo"
             print *,"fort_energyconst(1)=",fort_energyconst(1)
             stop
            endif

            ! this density will never be used at the wall.
            if (istate.eq.ENUM_DENVAR+1) then 
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=water_temp
            else
             print *,"istate invalid"
             stop
            endif 

           endif   ! probtype.eq.603

          else if (SDIM.eq.3) then

           ! ylo states
           if (probtype.eq.532) then
            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
            else
             print *,"istate invalid"
             stop
            endif
           endif  ! probtype=532

          else
           print *,"dimension bust"
           stop
          endif

          ! yhi
         else if ((dir.eq.2).and.(side.eq.2)) then
          if (xwall.gt.y) then
           print *,"xwall,y invalid"
           stop
          endif
          ADV=ADVwall

          if (SDIM.eq.2) then

            ! in: subroutine denBC
            ! yhi states 2d (freezing singularity problem, boiling)
           if (probtype.eq.59) then

            if (istate.eq.ENUM_DENVAR+1) then
             ! do nothing
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ! bcflag=1 (calling from denBC)
             call outside_temperature(time,x,y,z,ADV,im,1) 
            else
             print *,"istate invalid"
             stop
            endif

            ! yhi states
           else if (probtype.eq.801) then !denBC:compare w/Welch and Wilson 2000
            if (num_materials.lt.2) then
             print *,"num_materials invalid probtype=801"
             stop
            endif
            if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             call vapordist(xsten,nhalf,dx,bfact,dist)
             if (dist.lt.zero) then
              ADV=fort_tempconst(2) ! vapor/ice temperature
             else 
              ADV=fort_tempconst(1) ! liquid temperature then far field.
             endif
            endif

           ! Benard convection yhi
           else if (probtype.eq.603) then
            water_temp=fort_tempconst(1)
            call init_massfrac_parm(fort_denconst(1),massfrac_parm,1)
            call INTERNAL_material(fort_denconst(1),massfrac_parm, &
             water_temp, &
             fort_energyconst(1),fort_material_type(1),1)
            if (fort_energyconst(1).gt.zero) then
             ! do nothing
            else
             print *,"fort_energyconst(1) invalid in denBC yhi"
             print *,"fort_energyconst(1)=",fort_energyconst(1)
             stop
            endif

            ! this density will never be used at the wall.
            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=water_temp
            else
             print *,"istate invalid"
             stop
            endif 

           endif  

          else if (SDIM.eq.3) then

           ! 3d yhi states denBC
           if ((probtype.eq.532).or.(probtype.eq.541)) then
            if (istate.eq.ENUM_DENVAR+1) then
             ADV=fort_denconst(im)
            else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ADV=fort_tempconst(im)
            else
             print *,"istate invalid"
             stop
            endif
           endif  ! probtype=532

          else
           print *,"dimension bust"
           stop
          endif

           ! zlo
         else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then
          if (xwall.lt.z) then
           print *,"xwall,z invalid"
           stop
          endif
          ADV=ADVwall

           ! melting ice block: zlo
          if ((prescribe_temperature_outflow.eq.3).and. &
              (probtype.eq.59)) then

           if (istate.eq.ENUM_DENVAR+1) then
            ! do nothing (density)
           else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
             ! bcflag=1 (calling from denBC)
            call outside_temperature(time,x,y,z,ADV,im,1)
           else
            print *,"istate invalid"
            stop
           endif

          else if (probtype.eq.53) then

           if (istate.eq.ENUM_DENVAR+1) then
            ADV=fort_denconst(im)
           else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
            ADV=fort_tempconst(im)
           else
            print *,"istate invalid"
            stop
           endif 
           ! above probtype.eq.53
           ! zlo states
          else if ((probtype.eq.530).or.(probtype.eq.538).or. &
                   (probtype.eq.541)) then
           if (istate.eq.ENUM_DENVAR+1) then
            ADV=fort_denconst(im)
           else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
            ADV=fort_tempconst(im)
           else
            print *,"istate invalid"
            stop
           endif 
           ! above probtype=530 or 538
          else if (probtype.eq.603) then ! Benard convection zlo denBC
           water_temp=radblob2+fort_tempconst(1)
           call init_massfrac_parm(fort_denconst(1),massfrac_parm,1)
           call INTERNAL_material(fort_denconst(1),massfrac_parm, &
            water_temp, &
            fort_energyconst(1),fort_material_type(1),1)
           if (fort_energyconst(1).gt.zero) then
            ! do nothing
           else
            print *,"fort_energyconst(1) invalid in denBC zlo"
            print *,"fort_energyconst(1)=",fort_energyconst(1)
            stop
           endif

            ! this density will never be used at the wall.
           if (istate.eq.ENUM_DENVAR+1) then 
            ADV=fort_denconst(im) 
           else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
            ADV=water_temp
           else
            print *,"istate invalid"
            stop
           endif
   
          endif  ! probtype.eq.603

          ! zhi
         else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then
          if (xwall.gt.z) then
           print *,"xwall,z invalid"
           stop
          endif
          ADV=ADVwall

           ! Benard convection zhi
          if (probtype.eq.603) then
           water_temp=fort_tempconst(1)
           call init_massfrac_parm(fort_denconst(1),massfrac_parm,1)
           call INTERNAL_material(fort_denconst(1),massfrac_parm,  &
            water_temp, &
            fort_energyconst(1),fort_material_type(1),1)
           if (fort_energyconst(1).gt.zero) then
            ! do nothing
           else
            print *,"fort_energyconst(1) invalid in denBC zhi"
            print *,"fort_energyconst(1)=",fort_energyconst(1)
            stop
           endif

            ! this density will never be used at the wall.
           if (istate.eq.ENUM_DENVAR+1) then
            ADV=fort_denconst(im)
           else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
            ADV=water_temp
           else
            print *,"istate invalid"
            stop
           endif 

          endif  

         else
          print *,"dir side invalid"
          stop
         endif

         if (istate.eq.ENUM_DENVAR+1) then
          ! do nothing: den
         else if (istate.eq.ENUM_TEMPERATUREVAR+1) then
          ! do nothing: T
         else
          print *,"istate invalid"
          stop
         endif

        endif ! non-hydrate problems

        ADV_merge=ADV

       else if (species_flag.eq.1) then ! not density and not temperature

        ! CODY ESTEBE created test problem (species=0 at inflow)
        if ((probtype.eq.46).and.(axis_dir.eq.20)) then

         ADV=zero

        else if (probtype.eq.199) then

         call HYD_CCNT_BC(time,dir,side,ADV,xwall,ADVwall, &
          x,y,z,dx,im)

        else 

         if ((dir.eq.1).and.(side.eq.1)) then
          if (xwall.lt.x) then
           print *,"xwall,x invalid"
           stop
          endif
          ADV=ADVwall
         else if ((dir.eq.1).and.(side.eq.2)) then
          if (xwall.gt.x) then
           print *,"xwall,x invalid"
           stop
          endif
          ADV=ADVwall
         else if ((dir.eq.2).and.(side.eq.1)) then
          if (xwall.lt.y) then
           print *,"xwall,y invalid"
           stop
          endif
          ADV=ADVwall
         else if ((dir.eq.2).and.(side.eq.2)) then
          if (xwall.gt.y) then
           print *,"xwall,y invalid"
           stop
          endif
          ADV=ADVwall
         else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then
          if (xwall.lt.z) then
           print *,"xwall,z invalid"
           stop
          endif
          ADV=ADVwall
         else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then
          if (xwall.gt.z) then
           print *,"xwall,z invalid"
           stop
          endif
          ADV=ADVwall
         else
          print *,"dir side invalid"
          stop
         endif

        endif ! not hydrate problem

        ADV_merge=ADV

       else
        print *,"species_flag invalid"
        stop
       endif

       if (try_merge.eq.1) then ! temperature or species
        ADV=ADV_merge
       else if (try_merge.eq.0) then
        ! do nothing
       else
        print *,"try_merge invalid"
        stop
       endif
      else
       print *,"temp_homflag,species_homflag, or local_homflag invalid"
       stop
      endif  

      return
      end subroutine denBC

           ! sets all the EXT_DIR bcs to 0.0
      subroutine tensorBC(time,dir,side,ADV,ADVwall, &
        xsten,nhalf,dx,bfact,ipart,im)
      IMPLICIT NONE

      integer, INTENT(in) :: ipart,im
      integer, INTENT(in) :: dir,side,bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(out) :: ADV
      real(amrex_real), INTENT(in) :: ADVwall
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) :: x,y,z
      real(amrex_real) :: xwall

      if (nhalf.lt.1) then
       print *,"nhalf invalid tensorbc"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((ipart.lt.1).or. &
          (ipart.gt.num_materials_viscoelastic)) then
       print *,"ipart invalid:tensorBC"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid81:tensorBC"
       stop
      endif
      if (fort_im_elastic_map(ipart)+1.ne.im) then
       print *,"fort_im_elastic_map(ipart)+1.ne.im"
       stop
      endif

      x=xsten(0,1)
      y=xsten(0,2)
      z=xsten(0,SDIM)

      if (dir.eq.1) then
       if (side.eq.1) then
        xwall=problox
       else
        xwall=probhix
       endif
      else if (dir.eq.2) then
       if (side.eq.1) then
        xwall=probloy
       else
        xwall=probhiy
       endif
      else if ((dir.eq.3).and.(SDIM.eq.3)) then
       if (side.eq.1) then
        xwall=probloz
       else
        xwall=probhiz
       endif
      else
       print *,"dir invalid tensorbc"
       stop
      endif

      if ((dir.eq.1).and.(side.eq.1)) then
       if (xwall.lt.x) then
        print *,"xwall,x invalid"
        stop
       endif
!      ADV=ADVwall
       ADV=zero
      else if ((dir.eq.1).and.(side.eq.2)) then
       if (xwall.gt.x) then
        print *,"xwall,x invalid"
        stop
       endif
!      ADV=ADVwall
       ADV=zero
      else if ((dir.eq.2).and.(side.eq.1)) then
       if (xwall.lt.y) then
        print *,"xwall,y invalid"
        stop
       endif
!      ADV=ADVwall
       ADV=zero
      else if ((dir.eq.2).and.(side.eq.2)) then
       if (xwall.gt.y) then
        print *,"xwall,y invalid"
        stop
       endif
!      ADV=ADVwall
       ADV=zero
      else if ((dir.eq.3).and.(side.eq.1).and.(SDIM.eq.3)) then
       if (xwall.lt.z) then
        print *,"xwall,z invalid"
        stop
       endif
!      ADV=ADVwall
       ADV=zero
      else if ((dir.eq.3).and.(side.eq.2).and.(SDIM.eq.3)) then
       if (xwall.gt.z) then
        print *,"xwall,z invalid"
        stop
       endif
!      ADV=ADVwall
       ADV=zero
      else
       print *,"dir side invalid"
       stop
      endif

      return
      end subroutine tensorBC

      subroutine velread_bc_point(x,y,z,dir,velcomp,isann,time,vbc_point)

      IMPLICIT NONE
      integer dir,velcomp,isann
      real(amrex_real)    x,y,z,time,vbc_point


      integer nr,na
      real(amrex_real)    dt1,dt2,r,a0,a,daa,twopi
      real(amrex_real)    rbc(NR_MAX),vbc1(NS_MAX),vbc2(NS_MAX)

      ! interpolation on polar grid
      twopi=8.0*ATAN(1.0)
      if(dir.eq.1) then
       r = sqrt((y-yblob)**2+(z-zblob)**2)
       a = atan2((z-zblob),(y-yblob))
      else if(dir.eq.2) then
       r = sqrt((z-zblob)**2+(x-xblob)**2)
       a = atan2((x-xblob),(z-zblob))
      else if(dir.eq.3) then
       r = sqrt((x-xblob)**2+(y-yblob)**2)
       a = atan2((y-yblob),(x-xblob))
      endif

      ! Hack: this is taking too long
      if(velcomp.eq.3.and.r.le.radblob) then
       vbc_point = advbot
      else
       vbc_point = 0.
      endif
      return
        
      ! Unfortunately, the full file needs to be read for this point
      call velread_polar(dir,velcomp,isann,time, &
                         dt1,dt2,a0,rbc,vbc1,vbc2,nr,na)
      a = a-a0
      if (a.lt.0) then
       a = a+twopi
      endif
      daa = twopi/(na+1.0)    ! assume regular angular spacing
      if (isann.eq.0) then
       vbc_point = zero  ! initialize to zero only once
       if(r.lt.rbc(nr)) then
        call interp_polar(r,a,daa,dt1,dt2,rbc,vbc1,vbc2,nr,na,vbc_point)
       endif
      else if ((isann.eq.1).and.(r.lt.rbc(nr)).and.(r.gt.rbc(1))) then
       call interp_polar(r,a,daa,dt1,dt2,rbc,vbc1,vbc2,nr,na,vbc_point)
      endif

      return
      end subroutine


      subroutine interp_polar(r,a,daa,dt1,dt2,rbc,vbc1,vbc2,nr,na,vinterp)

      IMPLICIT NONE
      integer    nr,na
      real(amrex_real)       r,a,daa,dt1,dt2,vinterp
      real(amrex_real)       rbc(NR_MAX),vbc1(NS_MAX),vbc2(NS_MAX)

      integer l1,l2,m,m1,m2
      real(amrex_real) r1,r2,da1,da2,dr1,dr2,v1,v2
      real(amrex_real) vel111,vel112,vel121,vel122,vel211,vel212,vel221,vel222

       dr1 = 0.
       dr2 = 0.
       if(r.le.rbc(1)) then    ! extrapolate to the centerline i
                               ! (special case that happens for isann =
                               ! 0 only)
        dr2 = rbc(1)
        m1 = 1
        m2 = 1
       endif
       if(r.lt.rbc(nr).and.r.gt.rbc(1)) then  ! radial interpolation 
        do m = 2,nr
         if((r.gt.rbc(m-1)).and.(r.le.rbc(m))) then
          m1 = m-1
          m2 = m
          goto 30
         endif
        enddo
  30    r1 = rbc(m1)
        r2 = rbc(m2)
        dr1 = r-r1
        dr2 = r2-r
       endif

       if(dr2.gt.0) then          ! angle interpolation
        if ((m1.lt.1).or.(m2.lt.1).or.(m1.gt.nr).or.(m2.gt.nr)) then
         print *,"incorrect radial interpolation point(s) ",m1,m2
         stop
        endif
        l1 = ceiling(a/daa)
        if(l1.ge.na) then
         l1 = na
         l2 = 1
        else
         l2 = l1+1
        endif
        da1 = a-(l1-1)*daa
        da2 = daa-da1
        if ((l1.lt.1).or.(l2.lt.1).or.(l1.gt.na).or.(l2.gt.na)) then
         print *,"incorrect angular interpolation point(s) ",l1,l2
         stop
        endif
        vel111 = vbc1(m1+nr*(l1-1))
        vel121 = vbc1(m2+nr*(l1-1))
        vel112 = vbc1(m1+nr*(l2-1))
        vel122 = vbc1(m2+nr*(l2-1))
        v1 = vel111*dr2*da2+vel121*dr1*da2+vel112*dr2*da1+vel122*dr1*da1
        vel211 = vbc2(m1+nr*(l1-1))
        vel221 = vbc2(m2+nr*(l1-1))
        vel212 = vbc2(m1+nr*(l2-1))
        vel222 = vbc2(m2+nr*(l2-1))
        v2 = vel211*dr2*da2+vel221*dr1*da2+vel212*dr2*da1+vel222*dr1*da1
        vinterp = (dt2*v1+dt1*v2)/(dr1+dr2)/(da1+da2)/(dt1+dt2)
       else
        vinterp = 0.
       endif ! if(dr2.gt.0.)

      return
      end subroutine


      subroutine velread_polar(dir,velcomp,isann,time, &
                               dt1,dt2,a0,r_store,v1_store,v2_store,nr,na)

      ! Reads two consecutive velocity fields from a sequence of files
      ! spaced by 
      ! time "tinterv" and covering a total time "period"
      ! Two polar geometries allowed: filled circle and ring (isann =
      ! 0,1)

      IMPLICIT NONE

      integer    isann ! to distinguish between tube (0) and annular (1) geometry
      integer    dir,velcomp ! face direction,velocity component
      integer    nr,na ! number of radial and tangential points
      real(amrex_real) time,a0,r_store(NR_MAX),v1_store(NS_MAX),v2_store(NS_MAX)
      ! reference angle, radial positions, velocity fields


      character vel_file1*80, vel_file2*80
      integer i,j,nr2,na2,interv1,interv2
      real(amrex_real) rscale,rinner,router,dt1,dt2
      real(amrex_real) lref,vref,ptime1,ptime2,period,tinterv
      real(amrex_real), dimension(:), allocatable :: rbc
      real(amrex_real), dimension(:,:), allocatable :: vbc1
      real(amrex_real), dimension(:,:), allocatable :: vbc2

      if(isann.eq.0) then
       rinner = 0.
       router = radblob
       vref = advbot
      else
       rinner = radblob2
       router = radblob3
       vref = radblob4
      endif
      lref = router-rinner
      period = radblob5*lref/vref      ! overall duration of the records
      tinterv = radblob6*lref/vref     ! time spacing between records

      ! find correct time and open file
      ptime1 = mod(time,period)
      interv1 = floor(ptime1/tinterv)+1
      ! for interpolation in time
      dt1 = ptime1-(interv1-1)*tinterv
      ptime2 = mod(time+tinterv,period)
      interv2 = floor(ptime2/tinterv)+1
      dt2 = tinterv-dt1
      if ((dt1.lt.0.).or.(dt2.le.0)) then
       print *,"incorrect interpolation interval(s) in velread_polar",dt1,dt2
       stop
      endif

      ! Generate file name
      ! Note: synthetic turbulent inflow was generated along the x axis 
      if(isann.eq.0) then

       if(dir.eq.1) then
        if(velcomp.eq.1) then
         WRITE(vel_file1,111) interv1
         WRITE(vel_file2,111) interv2
        else if(velcomp.eq.2) then
         WRITE(vel_file1,222) interv1
         WRITE(vel_file2,222) interv2
        else if(velcomp.eq.3) then
         WRITE(vel_file1,333) interv1
         WRITE(vel_file2,333) interv2
        else
         print *,"invalid velocity component in velread_polar "
         stop
        endif
       else if(dir.eq.2) then
        if(velcomp.eq.1) then
         WRITE(vel_file1,333) interv1
         WRITE(vel_file2,333) interv2
        else if(velcomp.eq.2) then
         WRITE(vel_file1,111) interv1
         WRITE(vel_file2,111) interv2
        else if(velcomp.eq.3) then
         WRITE(vel_file1,222) interv1
         WRITE(vel_file2,222) interv2
        else
         print *,"invalid velocity component in velread_polar "
         stop
        endif
       else if(dir.eq.3) then
        if(velcomp.eq.1) then
         WRITE(vel_file1,222) interv1
         WRITE(vel_file2,222) interv2
        else if(velcomp.eq.2) then
         WRITE(vel_file1,333) interv1
         WRITE(vel_file2,333) interv2
        else if(velcomp.eq.3) then
         WRITE(vel_file1,111) interv1
         WRITE(vel_file2,111) interv2
        else
         print *,"invalid velocity component in velread_polar "
         stop
        endif
       endif

      else if(isann.eq.1) then

       if(dir.eq.1) then
        if(velcomp.eq.1) then
         WRITE(vel_file1,444) interv1
         WRITE(vel_file2,444) interv2
        else if(velcomp.eq.2) then
         WRITE(vel_file1,555) interv1
         WRITE(vel_file2,555) interv2
        else if(velcomp.eq.3) then
         WRITE(vel_file1,666) interv1
         WRITE(vel_file2,666) interv2
        else
         print *,"invalid velocity component in velread_polar "
         stop
        endif
       else if(dir.eq.2) then
        if(velcomp.eq.1) then
         WRITE(vel_file1,666) interv1
         WRITE(vel_file2,666) interv2
        else if(velcomp.eq.2) then
         WRITE(vel_file1,444) interv1
         WRITE(vel_file2,444) interv2
        else if(velcomp.eq.3) then
         WRITE(vel_file1,555) interv1
         WRITE(vel_file2,555) interv2
        else
         print *,"invalid velocity component in velread_polar "
         stop
        endif
       else if(dir.eq.3) then
        if(velcomp.eq.1) then
         WRITE(vel_file1,555) interv1
         WRITE(vel_file2,555) interv2
        else if(velcomp.eq.2) then
         WRITE(vel_file1,666) interv1
         WRITE(vel_file2,666) interv2
        else if(velcomp.eq.3) then
         WRITE(vel_file1,444) interv1
         WRITE(vel_file2,444) interv2
        else
         print *,"invalid velocity component in velread_polar "
         stop
        endif
       endif

      endif

      OPEN(unit=25,file=vel_file1,access='sequential', &
      form="formatted",status='old',err=901)
      READ(25,*) nr,na  ! number of radial and tangential components
      if ((nr.le.0).or.(na.le.0)) then
       print *,"Parameters ",nr,na, " are inconsistent ",nr,na
       stop
      endif
      if ((nr.gt.NR_MAX).or.(nr*na.gt.NS_MAX)) then
       print *,"Increase the buffer size: nr = ",nr," and na = ",nr,na
       stop
      endif
      allocate(rbc(nr))
      allocate(vbc1(nr,na))
      allocate(vbc2(nr,na))
      READ(25,*) a0  ! angle of the first coordinate
      READ(25,*) (rbc(i),i=1,nr)
      do j = 1,na
       READ(25,*) (vbc1(i,j),i=1,nr)
      enddo
      CLOSE(25)
      goto 10
 901  WRITE(6,911) vel_file1
      stop

  10  OPEN(unit=26,file=vel_file2,access='sequential', &
      form="formatted",status='old',err=902)
      READ(26,*) nr2,na2 ! should be the same as above
      if ((nr2.ne.nr).or.(na2.ne.na)) then
       print *,"the size ",nr2,na2, &
        " of the velocity tables does not match existing ",nr,na
       stop
      endif
      READ(26,*) a0
      READ(26,*) (rbc(i),i=1,nr)
      do j = 1,na
       READ(26,*) (vbc2(i,j),i=1,nr)
      enddo
      CLOSE(26)
      goto 20
 902  WRITE(6,911) vel_file2
      stop

 111  FORMAT ('./turbinlet/uinn',I3.3,'.dat')
 222  FORMAT ('./turbinlet/vinn',I3.3,'.dat')
 333  FORMAT ('./turbinlet/winn',I3.3,'.dat')
 444  FORMAT ('./turbinlet/uout',I3.3,'.dat')
 555  FORMAT ('./turbinlet/vout',I3.3,'.dat')
 666  FORMAT ('./turbinlet/wout',I3.3,'.dat')
 911  FORMAT(" Can't open file ",a," stopping ...")

  20  rscale = 2.*router      ! the normalized radius is 0.5
      do i = 1,nr
       r_store(i)=rscale*rbc(i)
       do j = 1,na
        v1_store(i+nr*(j-1)) = vref*vbc1(i,j)
        v2_store(i+nr*(j-1)) = vref*vbc2(i,j)
       enddo
      enddo

      deallocate(rbc)
      deallocate(vbc1)
      deallocate(vbc2)

      return
      end subroutine








! L rho= rho^* dx
! rho^*=L rho/dx
! xlo,xhi is bounding box for "side" neighbor of face
! side=0,1  dir=0,1
      subroutine override_mass_side(xlo,xhi,mass_side,delta,side,dir)
      IMPLICIT NONE

      integer side,dir
      real(amrex_real) xlo,xhi,mass_side,delta

     
      if ((side.ne.0).and.(side.ne.1)) then
       print *,"side invalid"
       stop
      endif
      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid override mass side"
       stop
      endif

      if (probtype.eq.102) then
       if (yblob4.le.zero) then
        print *,"yblob4 invalid"
        stop
       endif
       if (yblob5.le.zero) then
        print *,"yblob5 invalid"
        stop
       endif
       if (dir.eq.1) then

        if (side.eq.0) then
         if (half*(xlo+xhi).lt.zero) then
          mass_side=mass_side*yblob4/delta
         endif
        else if (side.eq.1) then
         if (half*(xlo+xhi).gt.probhiy) then
          mass_side=mass_side*yblob5/delta
         endif
        else
         print *,"side invalid"
         stop
        endif
 
       else if (dir.eq.0) then
 
       else
        print *,"dir invalid override mass side 2"
        stop
       endif

      endif

      return
      end subroutine override_mass_side

      subroutine override_facecut(facecut,xface,dir)
      IMPLICIT NONE

      real(amrex_real) facecut
      real(amrex_real) xface(SDIM)
      integer dir


      if (probtype.eq.102) then
       if (dir.eq.1) then
        if ((xface(2).le.EPS10).and.(xface(1).ge.radblob5).and. &
            (xface(1).le.radblob3)) then
         facecut=zero  ! override gas inflow
        endif
       else if (dir.eq.0) then

       else
        print *,"dir invalid override_facecut"
        stop
       endif
      endif

      return
      end subroutine override_facecut


! called if probtype=540
      subroutine get_Rieber_velocity(xsten,nhalf,bfact,dx,vel)
      IMPLICIT NONE

      integer, INTENT(in) :: nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real) cenbc(num_materials,SDIM)
      real(amrex_real), INTENT(out) :: vel(SDIM)

      integer dir2
      real(amrex_real) vfrac(num_materials)
      real(amrex_real) midz

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid get_Rieber_velocity"
       stop
      endif

      do dir2=1,SDIM
       vel(dir2)=zero
      enddo

      if (probtype.eq.540) then

       call get_initial_vfrac(xsten,nhalf,dx,bfact,vfrac,cenbc)  
       if (SDIM.eq.2) then
        midz=half*((yblob-radblob)+radblob2)
        if ((vfrac(1).gt.zero).and.(xsten(0,SDIM).ge.midz)) then  
         vel(SDIM)=-abs(advbot)
        endif
       else if (SDIM.eq.3) then
        midz=half*((zblob-radblob)+radblob2)
        if ((vfrac(1).gt.zero).and.(xsten(0,SDIM).ge.midz)) then  
         vel(SDIM)=-abs(advbot)
        endif
       else
        print *,"dimension bust"
        stop
       endif
 
      else
       print *,"probtype invalid"
       stop
      endif

      return
      end subroutine get_Rieber_velocity

        ! RGASRWATER is the ratio of the total domain width to the
        ! liquid domain width.
        ! 0<r<1 liquid
        ! 1<r<RGASRWATER  gas
      subroutine get_surface_tension_wave_speed(dxmin,Uscale,Lscale, &
        Re,We,RGASRWATER,wave_speed,density_ratio,viscosity_ratio, &
        wave_speed_target)
      IMPLICIT NONE


      real(amrex_real) wave_speed_target
      real(amrex_real) density_ratio,viscosity_ratio
      real(amrex_real) wave_speed,local_wave_speed
      real(amrex_real) Uscale,Lscale
      real(amrex_real) dxmin
      real(amrex_real) Re,We,RGASRWATER
      integer N1parm,N2parm
      real(amrex_real) :: alpha_real  ! set to 2 pi/dxmin
      real(amrex_real) :: alpha_imag  ! set to 0.0
      real(amrex_real) :: PI_LSA

      integer :: M, j, k
      complex*16 :: alpha
      complex*16, dimension(:,:), allocatable :: eigenv
      complex*16, dimension(:), allocatable :: omega1

      real(amrex_real), dimension(:), allocatable :: r1
      real(amrex_real), dimension(:), allocatable :: r2

      real(amrex_real)  :: omega_real,omega_imag,ki,kr
      real(amrex_real)  :: omega_before_real,omega_before_imag
      real(amrex_real)  :: err

      PI_LSA=4.0*ATAN(1.0)
      N1parm=64
      N2parm=N1parm
      alpha_imag=0.0
      alpha_real=2.0*PI_LSA*Lscale/dxmin

      print *,"surface tension alpha_real,alpha_imag ", &
         alpha_real,alpha_imag
      alpha=DCMPLX(alpha_real, alpha_imag)
      
     
        ! u, w, p, f  (v is radial velocity; not used in 2d)
        ! noslip and interface continuity conditions are hardwired.
      M=(N1parm+N2parm+2)*3+1-6

      allocate(r1(0:N1parm))
      allocate(r2(0:N2parm))

      allocate(omega1(M))
      allocate(eigenv(M,M))

      print *,"calling tlsa2d M=",M
      call tlsa2d(M, alpha, omega1, eigenv, N1parm, N2parm,  &
        r1, r2, Re,We,RGASRWATER,density_ratio,viscosity_ratio)
      print *,"done with tlsa2d"

      wave_speed=zero

      k=0
      do j=1,M
       kr=real(omega1(j))
       ki=aimag(omega1(j))
       omega_real=kr
       omega_imag=ki
       local_wave_speed=omega_real*Uscale/alpha_real

       if (1.eq.0) then

        print *,"**********************************************"
        print *,"the growthrate is: j,omega real, omega imag: ",   &
         j,omega_real,omega_imag

        print *,"dimensional speed: j,c ",j,local_wave_speed
        print *,"**********************************************"
     
       endif

       if (abs(local_wave_speed).le.wave_speed_target) then

        if (abs(local_wave_speed).ge.wave_speed_target*EPS6) then

         if (omega_imag.lt.zero) then

          do k=1,M
           if (k.ne.j) then 
          
            omega_before_real=real(omega1(k))
            omega_before_imag=aimag(omega1(k))
         
            if (omega_before_real*omega_real.lt.zero) then

             err=(omega_before_imag-omega_imag)/omega_imag
             if (abs(err).lt.EPS3) then
              err=(omega_before_real+omega_real)/omega_real
              if (abs(err).lt.EPS3) then
               if (abs(local_wave_speed).gt.wave_speed) then
                wave_speed=abs(local_wave_speed)
               endif
              endif
             endif

            endif ! real part opposite signs

           endif ! k<>j
          enddo ! k

         endif ! omega_imag<0

        endif ! computed wave speed > small fraction of predicted wave speed

       endif ! computed wave speed < than predicted wave speed
       
      enddo

! mode is e^{i(alpha z + n theta - omega t)}
! let n=0
! =e^{i (alpha z- omega t) } 
! the imaginary part of alpha is 0, so we have,
! =e^{i alpha (z-omega/alpha t)}
! let c=omega/alpha
! =e^{i alpha (z-ct)}=
! =e^{alpha imag(c)}e^{i alpha (z- real(c) t)}
! note: imag(c) should always be <=0 (no growing waves in 2d)
! dimensionless speed of waves is real(c).

      deallocate(r1)
      deallocate(r2)

      deallocate(omega1)
      deallocate(eigenv)

      return
      end subroutine


      subroutine gravity_wave_speed(wavelen,gravity_mag,wavespeed)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: wavelen
      real(amrex_real), INTENT(in) :: gravity_mag
      real(amrex_real), INTENT(out) :: wavespeed

      if ((wavelen.gt.zero).and. &
          (gravity_mag.gt.zero)) then
       ! do nothing
      else
       print *,"parameters invalid in gravity wave speed"
       stop
      endif

       ! Denner and van Wachem, JCP 285 (2015) 24-40 (27)
       ! omega/k = sqrt(g/k) 
       ! for SHALLOW (h<<1) waves: omega/k = sqrt(g h)  where h is the depth.
       ! Denner and van Wachem's formula is derived from:
       ! c=omega/k=sqrt((g/k) tanh(kh))  h=depth k=2 pi/lambda
       ! since 0<=tanh(kh)<=1 =>
       ! c<=sqrt(g/k)=sqrt(g * lambda/(2 pi))
      wavespeed=sqrt(gravity_mag * wavelen/(two*Pi))

      return
      end subroutine gravity_wave_speed


      subroutine capillary_wave_speed(wavelen,den1,den2,visc1,visc2, &
       tension,wavespeed)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: wavelen,den1,den2,visc1,visc2,tension
      real(amrex_real), INTENT(out) :: wavespeed
      real(amrex_real) omega,k

      if ((wavelen.gt.zero).and.(den1.gt.zero).and. &
          (den2.gt.zero).and.(visc1.gt.zero).and. &
          (visc2.gt.zero).and.(tension.gt.zero)) then
       ! do nothing
      else
       print *,"parameters invalid in capillary wave speed"
       stop
      endif

       ! Denner and van Wachem, JCP 285 (2015) 24-40 (1)
       ! omega/k = k^{1/2} (sigma/(den_liquid+den_gas))^{1/2}
       ! wavelen = dxmin
       ! wavespeed * dt < dx
       ! dt < dx * k / omega =
       ! dx / (sqrt(k * tension / (2 den ) ) ) =
       ! dx ^ {3/2}/sqrt( pi * tension / den ) =
       ! dx^{3/2} sqrt{den/(tension * pi)}
      k=two*Pi/wavelen
      omega=(k**(1.5d0))*sqrt(tension/min(den1,den2))
       ! wavespeed=omega/k=
       ! sqrt(k)*sqrt(tension/(den1+den2))=
       ! sqrt(2\pi tension/((den1+den2)*dx)=
       ! (den1+den2) replaced by min(den1,den2)
      wavespeed=omega/k 

      return
      end subroutine capillary_wave_speed


! distance to star with center at origin

      subroutine stardist(x,y,z,radstar,radthick,dist)
      use global_utility_module
      use global_distance_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x,y,z
      real(amrex_real), INTENT(out) :: dist
      real(amrex_real), INTENT(in) :: radstar,radthick
      real(amrex_real) :: dist1,dist2

      call cubedist(-radstar,radstar,-radstar,radstar, &
                    -radthick,radthick,x,y,z,dist1)
      call cubedist(-radthick,radthick,-radstar,radstar, &
                    -radstar,radstar,x,y,z,dist2)

      if ((dist1.le.zero).and.(dist2.le.zero)) then
       dist=-sqrt(dist1**2+dist2**2)
      else if (dist1.le.zero) then
       dist=dist1
      else if (dist2.le.zero) then
       dist=dist2
      else
       dist=dist1
       if (dist.gt.dist2) then
        dist=dist2
       endif
      endif

      return
      end subroutine

! rotate is in radians, positive angle rotates clockwise.
      subroutine paddlegeom(x,y,z,dist,onlypaddle)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: x,y,z
      real(amrex_real), INTENT(out) :: dist
      integer, INTENT(in) :: onlypaddle


      print *,"obsolete"
      stop

      return
      end subroutine paddlegeom


! u=2s_yz/pi^2
! v=-s_xz/pi^2
! w=-s_xy/pi^2


      subroutine jettingdist(x,y,dist,HSB,NOD,NPT)

      use global_distance_module

      IMPLICIT NONE
      real(amrex_real) x,y,dist,HSB,NOD,NPT,slope
      real(amrex_real) xx1(maxnumline),yy1(maxnumline)
      real(amrex_real) xx2(maxnumline),yy2(maxnumline)
      integer dd(maxnumline),lessflag(maxnumline),numline


      xx1(1)=zero
      xx2(1)=half*NOD
      yy1(1)=HSB+NPT
      yy2(1)=HSB+NPT
      dd(1)=0
      lessflag(1)=0 

      slope = -0.4
      xx1(2)=half*NOD
      yy1(2)=HSB+NPT
      xx2(2)=1.0D+3
      yy2(2)=slope*(xx2(2)-xx1(2))+yy1(2)
      dd(2)=0
      lessflag(2)=0

      numline=2

      call construct(x,y,xx1,yy1,xx2,yy2,dd,lessflag,numline,dist)
      if (y.ge.HSB+NPT) then
       dist=-abs(dist)
      endif

      if ((axis_dir.ge.8).and.(axis_dir.le.10)) then
       print *,"this option disabled"
       stop
      else if ((axis_dir.eq.11).or.(axis_dir.eq.12).or.(axis_dir.eq.13)) then
       call microfabdist(x,y,dist)
      endif
    
      return
      end subroutine jettingdist


! -------------------
! CODY ESTEBE
! -----------------
! representitive parameter values:
! values from (Singhal 2002)
! rho_l = 1.0 !density liquid (g/cm^3) (Brusiani: 0.845 g/cm^3 IDEAL fuel)
! rho_v = 0.00002558 ! density gas (g/cm^3)
! sigma = 71.7 !surface tension g/s^2 (Singhal), or can calculate
! mu = 0.00853 dynamic viscosity g/(cm s) 
! (water: 0.00853 at 300K)
! saturation pressure: 35400 g/(cm s^2)
! See also Bicer and Sou, 2016
! P=liquid pressure
! f_v=vapor mass fraction
! mu=liquid viscosity
! sigma=liquid/vapor surface tension.
subroutine RatePhaseChange(P,f_v,saturation_pressure, &
  rho_l,rho_v,mu,sigma,R_e,R_c,D_alpha_D_t)
 implicit none
 real(amrex_real), INTENT(in) :: P, f_v !pressure, vapor mass fraction
 real(amrex_real), INTENT(in) :: rho_l, rho_v ! density liquid, density vapor
 real(amrex_real), INTENT(in) :: sigma ! surface tension
 real(amrex_real), INTENT(in) :: mu ! dynamic viscosity of liquid
 real(amrex_real), INTENT(in) :: saturation_pressure ! saturation vapor pressure
 real(amrex_real) :: alpha_v
 real(amrex_real) :: k !k: turbulent kinetic energy
 real(amrex_real) :: C_e, C_c !empirical constants
 real(amrex_real) :: I, U, L, Re !initial turbulence intensity %, initial velocity magnitude, dynamic viscosity, length scale, Reynolds number (for obtaining k)
 real(amrex_real) :: rho
 real(amrex_real) :: P_v_Singhal  ! (18) in Singhal 2002
 real(amrex_real) :: saturation_pressure_local
 real(amrex_real) :: n_0,R_b,R_0
 real(amrex_real) :: P_injector
 real(amrex_real) :: R_crit
 real(amrex_real) :: P_crit
 real(amrex_real) :: h_l,h_v,h
 real(amrex_real) :: f_bar
 real(amrex_real) :: a,b,theta_0,psi,theta
 real(amrex_real) :: DfDt
 integer :: lowpressure
 real(amrex_real), INTENT(out) :: R_e, R_c !evaporation, condensation rates
 real(amrex_real), INTENT(out) :: D_alpha_D_t !volume fraction rate

 if ((rho_l.le.zero).or.(rho_v.le.zero)) then
  print *,"rho_l or rho_v invalid"
  stop
 endif
 if ((mu.le.zero).or.(sigma.le.zero)) then
  print *,"mu or sigma invalid"
  stop
 endif
 if ((f_v.lt.zero).or.(f_v.gt.one)) then
  print *,"f_v invalid"
  stop
 endif
 
 !T=300K (Singhal)
 U = 11000.0 !cm/s
 L = 0.03*0.07! length scale cm (7% nozzle diameter) (Brusiani)

 !(https://en.wikipedia.org/wiki/Turbulence_kinetic_energy)
 Re = rho_l*U*L/mu !Reynolds number
 I = 0.16*Re**(-1.0/8.0) !turbulent intensity
 k = (3.0/2.0)*((U*I)**2.0) !turbulent kinetic energy 

  ! note:
  ! rho=rho_l + alpha_v (rho_v - rho_l)
  ! alpha_v=(rho-rho_l)/(rho_v-rho_l)
  ! alpha_v rho_v = f_v rho
  ! rho=(1-alpha_v)rho_l + alpha_v rho_v=
  !     (1-f_v rho/rho_v)rho_l + f_v rho=
  ! rho rho_v =rho_l rho_v - f_v rho rho_l + f_v rho rho_v
  ! rho ((1-f_v)rho_v+f_v rho_l)= rho_l rho_v
  ! rho=rho_l rho_v/((1-f_v)rho_v+f_v rho_l)

 rho=rho_l*rho_v/((one-f_v)*rho_v+f_v*rho_l)
 if (rho.le.zero) then
  print *,"rho.le.zero"
  stop
 endif
 alpha_v=f_v*rho/rho_v ! alpha_v=vapor volume fraction f_v=vapor mass fraction
 if ((alpha_v.lt.zero).or.(alpha_v.gt.one)) then
  print *,"(alpha_v.lt.zero).or.(alpha_v.gt.one)"
  stop
 endif

 saturation_pressure_local=saturation_pressure ! 35400 (g/(cm s^2))
 if (saturation_pressure_local.le.zero) then
  print *,"saturation_pressure_local.le.zero"
  stop
 endif
 P_v_Singhal=saturation_pressure_local+0.39*rho*k/2.0 
 if (P_v_Singhal.le.zero) then
  print *,"P_v_Singhal.le.zero"
  print *,"P_v_Singhal=",P_v_Singhal
  print *,"saturation_pressure_local=",saturation_pressure_local
  print *,"k=",k
  print *,"rho=",rho
  stop
 endif

 if (1.eq.0) then
  print *,"saturation_pressure ",saturation_pressure
  print *,"P_v_Singhal ",P_v_Singhal
  print *,"P ",P
  print *,"rho_l=",rho_l
  print *,"rho_v=",rho_v
  print *,"f_v=",f_v
  print *,"k=",k
  print *,"sigma=",sigma
  stop
 endif
 
 C_e = 0.02
 C_c = 0.01

 !(Singhal 2002), (Brusiani et al 2013), Bicer and Sou 2016
 R_e = zero
 R_c = zero

  ! Singhal 2002
 if (1.eq.1) then

  if (P.lt.P_v_Singhal) then
   R_e = C_e*(sqrt(k)/sigma)*rho_l*rho_v* &
    sqrt((2.0/3.0)*(P_v_Singhal-P)/rho_l)*(one-f_v)
  else if (P.ge.P_v_Singhal) then
   R_c = C_c*(sqrt(k)/sigma)*rho_l*rho_v* &
    sqrt((2.0/3.0)*(P-P_v_Singhal)/rho_l)*f_v
  else
   print *,"P or  P_v_Singhal invalid"
   stop
  endif

 else if (1.eq.0) then ! Bicer and Sou
  n_0=10.0e+8  ! 1/cm^3
  R_b=(3.0*alpha_v/(4.0*Pi*(one-alpha_v)*n_0))**(1.0/3.0)
  R_0=1.0e-4
  P_injector=0.22e+7 ! gm/(cm s^2)  dyne/cm^2=gm cm/s^2  / cm^2
  if (P_injector.le.saturation_pressure_local) then
   print *,"P_injector.le.saturation_pressure_local"
   stop
  endif
  R_crit=R_0*sqrt(3.0*(R_0*(P_injector-saturation_pressure_local)/ &
         (2.0*sigma)+1.0))
  P_crit = saturation_pressure_local-(4.0*sigma)/(3.0*R_crit)
  if (P.lt.P_crit) then
   R_e = C_e*3.0*(rho_l*rho_v/rho)*alpha_v*(1.0-alpha_v)/R_b* &
         sqrt(2.0*abs(saturation_pressure_local-P)/(3.0*rho_l))
  else if (P.gt.saturation_pressure_local) then
   R_c = C_c*3.0*(rho_l*rho_v/rho)*alpha_v*(1.0-alpha_v)/R_b* &
         sqrt(1.27*abs(saturation_pressure_local-P)/rho_l)
  else if ((P.ge.P_crit).and.(P.le.saturation_pressure_local)) then
   R_e=zero
   R_c=zero
  else
   print *,"P, P_crit, or P_saturation_pressure_local invalid"
   stop
  endif

 else if (1.eq.0) then ! HRM
  n_0=10.0e+8  ! 1/cm^3
  R_b=(3.0*alpha_v/(4.0*Pi*(one-alpha_v)*n_0))**(1.0/3.0)
  R_0=1.0e-4
  P_injector=0.22e+7 ! gm/(cm s^2)  dyne/cm^2=gm cm/s^2  / cm^2
  if (P_injector.le.saturation_pressure_local) then
   print *,"P_injector.le.saturation_pressure_local"
   stop
  endif
  R_crit=R_0*sqrt(3.0*(R_0*(P_injector-saturation_pressure_local)/ &
         (2.0*sigma)+1.0))
  P_crit = saturation_pressure_local-(4.0*sigma)/(3.0*R_crit)
  h_l=P_injector/rho_l
  h_v=P_crit/rho_v
  h=P/rho
  if ((h.le.h_l).and.(h_v.lt.h_l)) then
   f_bar=(h_l-h)/(h_l-h_v)
  else
   print *,"h>h_l or h_v >= h_l"
   stop
  endif

  lowpressure = 0 !flag low or high pressure (below 10bar)
  if (lowpressure.eq.1) then
   a=-0.257
   b=-2.24
   theta_0=6.51e-4
   psi = abs((saturation_pressure_local-P)/saturation_pressure_local)
   theta = theta_0*(alpha_v**a)*(psi**b)
  else if (lowpressure.eq.0) then
   a=-0.54
   b=-1.76
   theta_0=6.51e-4
   psi = abs((saturation_pressure_local-P)/ &
             (P_crit-saturation_pressure_local))
   theta = theta_0*(alpha_v**a)*(psi**b)
  else
   print *,"lowpressure invalid"
   stop
  endif
   ! DfDt = rate of change of the vapor mass fraction
   ! rho DfDt = volume D (rho f)/Dt = rate of change of the vapor mass
  DfDt = (f_bar-f_v)/theta
   ! the effective velocity (assuming unidirectional normal) is:
   ! velocity=dx * D_alpha_D_t
  D_alpha_D_t=rho*DfDt/rho_v
  if (DfDt.eq.zero) then
   ! do nothing
  else if (DfDt.gt.zero) then
   R_e=DfDt*rho
  else if (DfDt.lt.zero) then
   R_c=DfDt*rho
  else
   print *,"DfDt invalid"
   stop
  endif

 else
  print *,"option does not exist"
  stop
 endif

 if ((R_e.lt.zero).or.(R_c.lt.zero)) then
  print *,"(R_e.lt.zero).or.(R_c.lt.zero)"
  stop
 endif
 
end subroutine RatePhaseChange


       ! vapor_mass_source<0 if condensation
       ! vapor_mass_source>0 if evaporation
      subroutine get_vapor_mass_source( &
       liq_viscosity, &
       liq_vap_tension, &
       total_density, &
       liquid_density,vapor_density, &
       liquid_pressure, &
       saturation_pressure, &
       vapor_mass_frac, &
       vapor_vol_frac, &
       vapor_mass_source, &
       vapor_vfrac_source)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: liq_viscosity
      real(amrex_real), INTENT(in) :: liq_vap_tension
      real(amrex_real), INTENT(in) :: total_density
      real(amrex_real), INTENT(in) :: liquid_density,vapor_density
      real(amrex_real), INTENT(in) :: liquid_pressure,saturation_pressure
      real(amrex_real), INTENT(in) :: vapor_mass_frac
      real(amrex_real), INTENT(in) :: vapor_vol_frac
      real(amrex_real), INTENT(out) :: vapor_mass_source
      real(amrex_real), INTENT(out) :: vapor_vfrac_source

      real(amrex_real) R_e,R_c

      if ((liquid_density.gt.zero).and.(vapor_density.gt.zero)) then
       ! do nothing
      else
       print *,"liquid_density or vapor_density invalid"
       stop
      endif
      if (liq_vap_tension.ge.zero) then
       ! do nothing
      else
       print *,"liq_vap_tension invalid"
       stop
      endif

      if ((vapor_mass_frac.ge.zero).and.(vapor_mass_frac.le.one)) then

       if ((vapor_vol_frac.ge.zero).and.(vapor_vol_frac.le.one)) then

         ! R_e>=0  R_c>=0
        call RatePhaseChange( &
         liquid_pressure, &
         vapor_mass_frac, &
         saturation_pressure, &
         liquid_density,vapor_density, &
         liq_viscosity,liq_vap_tension, &
         R_e,R_c, &
         vapor_vfrac_source)

        if ((R_e.ge.zero).and.(R_c.ge.zero)) then
         vapor_mass_source=R_e-R_c
        else
         print *,"R_e or R_c invalid"
         stop
        endif 

       else
        print *,"vapor_vol_frac invalid"
        stop
       endif
 
      else
       print *,"vapor_mass_frac invalid"
       stop
      endif

      return
      end subroutine get_vapor_mass_source

! -------------------
! end CODY ESTEBE
! -----------------


      subroutine get_vel_phasechange_NUCLEATE( &
                      nucleate_in,nucleate_out)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

       ! nucleation_parm_type_input is declared in PROBCOMMON.F90
      type(nucleation_parm_type_input), INTENT(in) :: nucleate_in
      type(nucleation_parm_type_inout), INTENT(inout) :: nucleate_out
      real(amrex_real) VOFTOL_NUCLEATE
      integer i,j,k
      integer denbase
      integer mofbase
      integer vofcomp
      integer local_freezing_model
      integer im_local
      integer im_dest
      integer im_source
      integer im_vapor
      integer im_liquid
      real(amrex_real) LL
      integer make_seed ! 0 no seed, 1 yes to seed
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) prev_time,cur_time,dt
      real(amrex_real) saturation_pres
      real(amrex_real) total_density
      real(amrex_real) test_pressure
      real(amrex_real) test_den
      real(amrex_real) vapor_density
      real(amrex_real) test_temp
      real(amrex_real) test_visc
      real(amrex_real) vapor_vol_frac
      real(amrex_real) vapor_mass_frac
      real(amrex_real) liquid_mass_frac
      real(amrex_real) vapor_mass_source
      real(amrex_real) vapor_vfrac_source
      integer nmax
      integer, parameter :: use_ls_data=0
      integer, parameter :: mof_verbose=0
      integer, parameter :: continuous_mof=STANDARD_MOF
      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))

      integer :: grid_index(SDIM)
      integer, parameter :: grid_level=-1

      integer tessellate
      integer ibasesrc,ibasedst
      integer ibase_raw,ibase_recon

      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) vof_super(num_materials)

      real(amrex_real) :: LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)
      real(amrex_real) :: multi_centroidA(num_materials,SDIM)
      real(amrex_real) :: volcell
      real(amrex_real) :: cencell(SDIM)
      integer :: dir
      real(amrex_real) cen_src(SDIM)
      real(amrex_real) cen_dst(SDIM)
      real(amrex_real) vfluid_sum
      real(amrex_real) VOF_source,VOF_dest
      integer local_tessellate

      nmax=POLYGON_LIST_MAX 

      call checkbound_array(nucleate_in%fablo,nucleate_in%fabhi, &
        nucleate_in%EOS,1,-1)
      call checkbound_array1(nucleate_in%fablo,nucleate_in%fabhi, &
        nucleate_in%pres,1,-1)
      call checkbound_array1(nucleate_in%fablo,nucleate_in%fabhi, &
        nucleate_in%pres_eos,1,-1)
      call checkbound_array(nucleate_in%fablo,nucleate_in%fabhi, &
        nucleate_in%Snew,1,-1)
      call checkbound_array(nucleate_in%fablo,nucleate_in%fabhi, &
        nucleate_in%LSnew,1,-1)

      if (nucleate_in%nstate.eq.STATE_NCOMP) then
       ! do nothing
      else
       print *,"nstate invalid"
       stop
      endif
      if (n_sites.gt.0) then
       if (nucleate_in%nucleate_pos_size.ne.n_sites*4) then
        print *,"nucleate_pos_size invalid"
        stop
       endif
      endif

       ! redistancing for phase change sees materials in which
       ! F>EPS3
      VOFTOL_NUCLEATE=EPS3*two

      denbase=STATECOMP_STATES
      mofbase=STATECOMP_MOF
      im_dest=nucleate_in%im_dest
      im_source=nucleate_in%im_source
      vofcomp=mofbase+(im_dest-1)*ngeom_raw+1

      local_freezing_model=nucleate_in%local_freezing_model

      i=nucleate_in%i
      j=nucleate_in%j
      k=nucleate_in%k
      if (nucleate_out%Snew(D_DECL(i,j,k),vofcomp).ge.VOFTOL_NUCLEATE) then
       ! do nothing (already destination material in the cell)
      else if (nucleate_out%Snew(D_DECL(i,j,k),vofcomp).le.VOFTOL_NUCLEATE) then

       LL=nucleate_in%LL
       make_seed=0

       grid_index(1)=i
       grid_index(2)=j
       if (SDIM.eq.3) then
        grid_index(SDIM)=k
       endif

       call gridsten_level(xsten,i,j,k,nucleate_in%level,nhalf)

       if ((local_freezing_model.eq.0).or. & !stefan model freezing/boiling
           (local_freezing_model.eq.5).or. & !stefan model evap/condensation
           (local_freezing_model.eq.6)) then !TSAT variable evap/condensation

        if (is_in_probtype_list().eq.1) then
         call SUB_nucleation(nucleate_in,xsten,nhalf,make_seed)
        else
         ! do nothing
        endif
        
       else if (local_freezing_model.eq.7) then ! cavitation

        prev_time=nucleate_in%prev_time
        cur_time=nucleate_in%cur_time
        dt=nucleate_in%dt
        if (prev_time.gt.zero) then
         if (cur_time.gt.prev_time) then
          if (dt.gt.zero) then
           test_pressure=nucleate_in%pres_eos(D_DECL(i,j,k))
           if (LL.gt.zero) then
            im_vapor=im_dest ! evaporation
            im_liquid=im_source
            vapor_mass_frac=zero
            liquid_mass_frac=one
            vapor_density=fort_denconst(im_vapor)
            test_den=nucleate_in%EOS(D_DECL(i,j,k), &
                   (im_liquid-1)*num_state_material+1)
            test_temp=nucleate_in%EOS(D_DECL(i,j,k), &
                   (im_liquid-1)*num_state_material+2)
           else if (LL.lt.zero) then
            im_vapor=im_source ! condensation
            im_liquid=im_dest
            vapor_mass_frac=one
            liquid_mass_frac=zero
            vapor_density=nucleate_in%EOS(D_DECL(i,j,k), &
                   (im_vapor-1)*num_state_material+1)
            test_den=fort_denconst(im_liquid)
            test_temp=nucleate_in%EOS(D_DECL(i,j,k), &
                   (im_vapor-1)*num_state_material+2)
           else
            print *,"LL invalid"
            stop
           endif
           saturation_pres=nucleate_in%cavitation_pressure(im_liquid)
           if (vapor_density.gt.zero) then
            ! let Y=vapor_mass_frac dV=vapor density dL=liquid den
            ! alpha*dV=Y*d
            ! d=alpha dV + (1-alpha) dL=alpha(dV-dL)+dL=
            ! Y*d(dV-dL)/dV+dL
            ! d(1-Y(dV-dL)/dV)=dL
            ! d=dL dV/(dV-Y(dV-dL))=dL dV/((1-Y)dV+Y dL)
            ! 
            total_density=test_den*vapor_density/ &
               (vapor_mass_frac*test_den+liquid_mass_frac*vapor_density)
            vapor_vol_frac=total_density*vapor_mass_frac/vapor_density

            if ((vapor_vol_frac.ge.zero).and. &
                (vapor_vol_frac.le.one)) then

             if (nucleate_in%cavitation_tension(im_liquid).ge.zero) then
              test_visc=get_user_viscconst(im_liquid,test_den,test_temp)

              if (test_visc.ge.zero) then
               call get_vapor_mass_source( &
                 test_visc, &
                 nucleate_in%cavitation_tension(im_liquid), &
                 total_density, &
                 test_den,vapor_density, &
                 test_pressure, &
                 saturation_pres, &
                 vapor_mass_frac, &
                 vapor_vol_frac, &
                 vapor_mass_source, &
                 vapor_vfrac_source)
              else
               print *,"test_visc invalid"
               stop
              endif

              if (vapor_vfrac_source*LL.gt.zero) then
               make_seed=1
              else if (vapor_vfrac_source*LL.le.zero) then
               ! do nothing
              else
               print *,"vapor_vfrac_source invalid"
               stop
              endif
             else
              print *,"cavitation_tension invalid"
              stop
             endif
            else
             print *,"vapor_vol_frac invalid"
             stop
            endif
           else
            print *,"vapor_density invalid"
            stop
           endif
          else
           print *,"dt invalid"
           stop
          endif
         else
          print *,"cur_time invalid"
          stop
         endif
        else if (prev_time.eq.zero) then
         ! do nothing
        else
         print *,"prev_time invalid"
         stop
        endif
       else if (is_valid_freezing_modelF(local_freezing_model).eq.1) then
        ! do nothing
       else
        print *,"local_freezing_model invalid 10"
        stop
       endif

       if (make_seed.eq.1) then
        ! create unidirectional seed
        ! initialize volume fraction, level set function, centroid,
        ! density, temperature, and vapor mass fraction.
        tessellate=0

        do im_local=1,num_materials
         ibase_raw=(im_local-1)*ngeom_raw+1
         ibase_recon=(im_local-1)*ngeom_recon+1
         do dir=0,SDIM
          mofdata(ibase_recon+dir)= &
            nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibase_raw+dir)
         enddo
          ! order=0
         mofdata(ibase_recon+SDIM+1)=zero
         do dir=SDIM+2,ngeom_recon-1
          mofdata(ibase_recon+dir)=zero  ! slope, intercept
         enddo
        enddo  ! im_local=1..num_materials

        call make_vfrac_sum_ok_base( &
          cmofsten, &
          xsten,nhalf, &
          continuous_mof, &
          nucleate_in%bfact, &
          nucleate_in%dx, &
          tessellate, &  ! =0
          mofdata, &
          SDIM)

        do im_local=1,num_materials
         ibase_recon=(im_local-1)*ngeom_recon+1
         vof_super(im_local)=mofdata(ibase_recon)
        enddo

        call multimaterial_MOF( &
         nucleate_in%bfact, &
         nucleate_in%dx, &
         xsten,nhalf, &
         mof_verbose, & ! =0
         use_ls_data, & ! =0
         LS_stencil, &
         geom_xtetlist(1,1,1,nucleate_in%tid+1), &
         geom_xtetlist(1,1,1,nucleate_in%tid+1), &
         nmax, &
         nmax, &
         mofdata, & !intent(inout)
         vof_super, &
         multi_centroidA, &
         continuous_mof, & ! =STANDARD_MOF
         cmofsten, &
         grid_index, &
         grid_level, &
         SDIM)

         ! if local_tessellate==3: (rasterized reconstruction for solids)
         !  if solid material(s) dominate the cell, then F_solid_raster=1
         !  and F_fluid=0.
         !  if fluid material(s) dominate the cell, then F_solid=0,
         !  sum F_fluid=1
        local_tessellate=3
         !EPS2
        call multi_get_volume_tessellate( &
         local_tessellate, & ! =3
         nucleate_in%bfact, &
         nucleate_in%dx, &
         xsten,nhalf, &
         mofdata, &
         geom_xtetlist(1,1,1,nucleate_in%tid+1), &
         nmax, &
         nmax, &
         SDIM)

        call CISBOX( &
         xsten,nhalf, &
         nucleate_in%xlo, &
         nucleate_in%dx, &
         i,j,k, &
         nucleate_in%bfact, &
         nucleate_in%level, &
         volcell,cencell,SDIM)
 
        ibasesrc=(im_source-1)*ngeom_recon+1

         ! is there enough "source material" in the cell to be 
         ! converted to "destination material"?
        if (mofdata(ibasesrc).gt.VOFTOL_NUCLEATE) then

         nucleate_out%LSnew(D_DECL(i,j,k),im_source)=-nucleate_in%dx(1)
         nucleate_out%LSnew(D_DECL(i,j,k),im_dest)=nucleate_in%dx(1)
          !levelset gradient=0 in the source and dest materials.
         do dir=1,SDIM
          nucleate_out% &
            LSnew(D_DECL(i,j,k),num_materials+SDIM*(im_source-1)+dir)=zero
          nucleate_out% &
            LSnew(D_DECL(i,j,k),num_materials+SDIM*(im_dest-1)+dir)=zero
         enddo

         vfluid_sum=zero
         do im_local=1,num_materials
          if (is_rigid(im_local).eq.0) then
           ibase_raw=(im_local-1)*ngeom_raw+1
           vfluid_sum=vfluid_sum+ &
            nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibase_raw)
          else if (is_rigid(im_local).eq.1) then
           ! do nothing
          else
           print *,"is_rigid(im_local) invalid"
           stop
          endif
         enddo ! im_local=1..num_materials
         if (vfluid_sum.gt.VOFTOL_NUCLEATE) then
          ibasesrc=(im_source-1)*ngeom_raw+1
          VOF_source=nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibasesrc)
          if (VOF_source.gt.VOFTOL_NUCLEATE) then
           nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibasesrc)=zero
           ibasedst=(im_dest-1)*ngeom_raw+1
           VOF_dest=nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibasedst)
           nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibasedst)= &
              VOF_dest+VOF_source
           if (VOF_dest+VOF_source.gt.zero) then
            do dir=1,SDIM
             cen_dst(dir)=nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibasedst+dir)
             cen_src(dir)=nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibasesrc+dir)
             nucleate_out%Snew(D_DECL(i,j,k),mofbase+ibasedst+dir)= &
               (VOF_dest*cen_dst(dir)+VOF_source*cen_src(dir))/ &
               (VOF_dest+VOF_source)
            enddo ! dir=1..sdim
           else
            print *,"VOF_dest+VOF_source invalid"
            stop
           endif
          else
           print *,"VOF_source invalid"
           stop
          endif
         else
          print *,"vfluid_sum invalid; fluid volume fractions should"
          print *,"tessellate the domain.  vfluid_sum=",vfluid_sum
          stop
         endif
        else
         print *,"mofdata(ibasesrc) invalid: ",mofdata(ibasesrc)
         print *,"(expecting mofdata(ibasesrc)>VOFTOL_NUCLEATE)"
         stop
        endif
           
       else if (make_seed.eq.0) then
        ! do nothing
       else
        print *,"make_seed invalid"
        stop
       endif
      else
       print *,"Snew(vofcomp) invalid"
       stop
      endif

      end subroutine get_vel_phasechange_NUCLEATE

       ! distribute_from_target=0 => distribute to the destination material
       ! => velocity in the source will be divergence free (continuous) => 
       ! velocity on the interface will correspond to the source velocity.
       !  V=u_src dot n + mdot/rho_src
       ! distribute_from_target=1 => distribute to the source material
       ! => velocity in the destination will be divergence free (continuous) => 
       ! velocity on the interface will correspond to the destination velocity.
       !  V=u_dst dot n + mdot/rho_dst
       !
       ! expansion_fact=
       ! either: 1-den_dst/den_src (distribute_from_target==0)
       !     or: 1-den_src/den_dst (distribute_from_target==1)
      subroutine get_vel_phasechange( &
       interface_mass_transfer_model, &
       for_estdt, &
       xI, &
       ispec, &
       molar_mass, &
       species_molar_mass, &
       local_freezing_model, &
         ! 1=Tanasawa  2=Schrage 3=Kassemi 
       local_Tanasawa_or_Schrage_or_Kassemi, & 
       distribute_from_target, &
       vel, &
       densrc_I,dendst_I, & ! replaced with vapor_den if freezing_model=5,6
       densrc_probe,dendst_probe, &
       ksrc_derived, &
       kdst_derived, &
       ksrc_physical, &
       kdst_physical, &
       Tsrc_probe, &
       Tdst_probe, &
       Tsat, &
       Tsrc_INT,Tdst_INT, &
       LL, &
       source_perim_factor, &
       dest_perim_factor, &
       microlayer_substrate_source, &
       microlayer_angle_source, &
       microlayer_size_source, &
       macrolayer_size_source, &
       microlayer_substrate_dest, &
       microlayer_angle_dest, &
       microlayer_size_dest, &
       macrolayer_size_dest, &
       dxprobe_source, &
       dxprobe_dest, &
       im_source,im_dest, &
       time,dt, &
       alpha, &
       beta, &
       expansion_fact, &
       K_f, &
       Cmethane_in_hydrate, &
       C_w0, &
       PHYDWATER, &
       VOFsrc,VOFdst) ! VOFsrc,VOFdst used for the hydrate model.
      use hydrateReactor_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: interface_mass_transfer_model
      integer, INTENT(in) :: for_estdt
      real(amrex_real), INTENT(in) :: xI(SDIM)
      integer, INTENT(in) :: local_freezing_model
!1=Tanasawa 2=Schrage 3=Kassemi 
      integer, INTENT(in) :: local_Tanasawa_or_Schrage_or_Kassemi
       ! MEHDI EVAPORATION
      integer, INTENT(in) :: ispec ! 0 if no species  1..num_species_var
       ! MEHDI EVAPORATION
      real(amrex_real), INTENT(in) :: molar_mass(num_materials)
      real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var+1)
      integer, INTENT(in) :: distribute_from_target
       ! MEHDI EVAPORATION  im_source,im_dest = 1..num_materials
      integer, INTENT(in) :: im_source,im_dest
      real(amrex_real), INTENT(out) :: vel
      real(amrex_real), INTENT(in) :: densrc_I,dendst_I
      real(amrex_real), INTENT(in) :: densrc_probe,dendst_probe
      real(amrex_real), INTENT(in) :: time,dt,alpha,beta
      real(amrex_real), INTENT(in) :: expansion_fact
      real(amrex_real), INTENT(in) :: ksrc_derived,kdst_derived
      real(amrex_real), INTENT(in) :: ksrc_physical,kdst_physical
      real(amrex_real), INTENT(in) :: Tsrc_probe
      real(amrex_real), INTENT(in) :: Tdst_probe
      real(amrex_real), INTENT(in) :: Tsat
      real(amrex_real), INTENT(in) :: LL
      real(amrex_real), INTENT(in) :: source_perim_factor,dest_perim_factor
      integer, INTENT(in) :: microlayer_substrate_source
      integer, INTENT(in) :: microlayer_substrate_dest
      real(amrex_real), INTENT(in) :: microlayer_angle_source,microlayer_angle_dest
      real(amrex_real), INTENT(in) :: microlayer_size_source,microlayer_size_dest
      real(amrex_real), INTENT(in) :: macrolayer_size_source,macrolayer_size_dest
      real(amrex_real), INTENT(in) :: dxprobe_source
      real(amrex_real), INTENT(in) :: dxprobe_dest
      real(amrex_real), INTENT(in) :: Tsrc_INT,Tdst_INT
      real(amrex_real), INTENT(in) :: K_f
      real(amrex_real), INTENT(in) :: Cmethane_in_hydrate,C_w0,PHYDWATER
      real(amrex_real), INTENT(in) :: VOFsrc,VOFdst

      real(amrex_real) DTsrc,DTdst
      real(amrex_real) velsrc,veldst,velsum
      real(amrex_real) velsrc_micro,veldst_micro
      real(amrex_real) psi_upper,psi_lower,micro_slope
      integer verb_hydrate
      integer mdot_override
      real(amrex_real) mdot

      if ((im_source.lt.1).or.(im_source.gt.num_materials)) then
       print *,"im_source invalid"
       stop
      endif
      if ((im_dest.lt.1).or.(im_dest.gt.num_materials)) then
       print *,"im_dest invalid"
       stop
      endif
      if ((densrc_I.gt.zero).and.(dendst_I.gt.zero)) then
       ! do nothing
      else
       print *,"density I must be positive"
       stop
      endif
      if ((densrc_probe.gt.zero).and.(dendst_probe.gt.zero)) then
       ! do nothing
      else
       print *,"density probe be positive"
       stop
      endif
      if ((distribute_from_target.ne.0).and. &
          (distribute_from_target.ne.1)) then
       print *,"distribute_from_target invalid"
       stop
      endif
      if (is_valid_freezing_modelF(local_freezing_model).eq.1) then
       ! do nothing
      else
       print *,"local_freezing_model invalid 11"
       stop
      endif
      if (is_hydrate_freezing_modelF(local_freezing_model).eq.1) then 
       if (num_species_var.eq.1) then
        ! do nothing
       else
        print *,"must define species var if hydrate model"
        stop
       endif
      else if (is_hydrate_freezing_modelF(local_freezing_model).eq.0) then 
       ! do nothing
      else
       print *,"is_hydrate_freezing_modelF invalid"
       stop
      endif
      if (is_multi_component_evapF(local_freezing_model, &
           local_Tanasawa_or_Schrage_or_Kassemi,LL).eq.1) then
       if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
        ! do nothing
       else
        print *,"ispec invalid"
        stop
       endif
      else if (is_multi_component_evapF(local_freezing_model, &
           local_Tanasawa_or_Schrage_or_Kassemi,LL).eq.0) then
       if (ispec.eq.0) then
        ! do nothing
       else
        print *,"ispec invalid"
        stop
       endif
      else
       print *,"is_multi_component_evapF invalid"
       stop
      endif

      if ((VOFsrc.ge.-EPS1).and. &
          (VOFdst.ge.-EPS1).and. &
          (VOFsrc.le.one+EPS1).and. &
          (VOFdst.le.one+EPS1)) then
       ! do nothing
      else
       print *,"VOFsrc or VOFdst invalid: ",VOFsrc,VOFdst
       stop
      endif

      if ((Tsrc_probe.ge.zero).and.(Tdst_probe.ge.zero).and. &
          (Tsat.ge.zero).and.(Tsrc_INT.ge.zero).and. &
          (Tdst_INT.ge.zero)) then
       ! do nothing
      else
       print *,"temperature cannot be negative in get_vel_phasechange"
       print *,"Tsrc_probe,Tdst_probe,Tsat,TsrcI,TdstI ", &
         Tsrc_probe, Tdst_probe, Tsat, Tsrc_INT, Tdst_INT
       print *,"for_estdt= ",for_estdt
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
      if (alpha.ge.zero) then
       ! do nothing
      else
       print *,"alpha invalid"
       stop
      endif
      if (beta.ge.zero) then
       ! do nothing
      else
       print *,"beta invalid"
       stop
      endif
       ! expansion_fact=
       ! either: 1-den_dst/den_src
       !     or: 1-den_src/den_dst
      if (expansion_fact.lt.one) then
       ! do nothing
      else
       print *,"expansion_fact invalid get_vel_phasechange:",expansion_fact
       stop
      endif
      if ((for_estdt.ne.0).and.(for_estdt.ne.1)) then
       print *,"for_estdt invalid"
       stop
      endif
      if ((source_perim_factor.gt.zero).and. &
          (source_perim_factor.le.one).and. &
          (dest_perim_factor.gt.zero).and. &
          (dest_perim_factor.le.one)) then
       ! do nothing
      else
       print *,"source_perim_factor or dest_perim_factor invalid"
       stop
      endif

      if ((probtype.eq.801).and.(vinletgas.gt.zero)) then
         ! in: get_vel_phasechange
         ! fixed rate of phase change (sanity check)
        vel=vinletgas  ! MDOT/rho_src

        ! local_freezing_model=0 (sharp interface stefan model)
        ! local_freezing_model=1 (source term model)
        ! local_freezing_model=2 (hydrate model)
        ! local_freezing_model=3 (wildfire)
        ! local_freezing_model=5 (Stefan model evaporation/condensation)
        ! local_freezing_model=6 (Palmore and Desjardins)
        ! local_freezing_model=7 (Cavitation)
      else if ((is_GFM_freezing_modelF(local_freezing_model).eq.1).or. &
                (local_freezing_model.eq.1)) then
         ! LL<0 if freezing

        if ((LL.ne.zero).and. &
            (dxprobe_source.gt.zero).and. &
            (dxprobe_dest.gt.zero)) then
         ! do nothing
        else
         print *,"LL, dxprobe_source, or dxprobe_dest invalid"
         stop
        endif

         ! local_freezing_model==5, if Tsrc_probe > Tsat then
         ! evaporation will occur.  (source==water destination==vapor within
         !                           the air)
         ! if local_freezing_model==5, and Tsrc_probe < Tsat then
         !  velsrc<0 => then the "rate of mass transfer is negative" which
         ! is disallowed; i.e. no evaporation occurs.
        ! Tsrc_probe is the probe temperature in the source
        DTsrc=Tsrc_probe-Tsat 
    
        ! Tdst_probe is the probe temperature in the destination
        DTdst=Tdst_probe-Tsat  

        velsrc=ksrc_derived*DTsrc/(LL*dxprobe_source)
        veldst=kdst_derived*DTdst/(LL*dxprobe_dest)

        if (local_freezing_model.eq.0) then

         velsrc_micro=zero
         veldst_micro=zero
       
         if (microlayer_size_source.gt.zero) then
          if ((microlayer_substrate_source.ge.1).and. &
              (microlayer_substrate_source.le.num_materials)) then
           if (is_rigid(microlayer_substrate_source).ne.1) then
            print *,"is_rigid(microlayer_substrate_source).ne.1"
            stop
           endif
           if (macrolayer_size_source.gt.microlayer_size_source) then
            if ((microlayer_angle_source.gt.zero).and. &
                (microlayer_angle_source.lt.Pi)) then
             micro_slope=tan(half*microlayer_angle_source)
             psi_upper=half*macrolayer_size_source/micro_slope
             psi_lower=half*microlayer_size_source/micro_slope
             velsrc_micro=( &
              half*(Tsat+fort_tempconst(microlayer_substrate_source))-Tsat)* &
              log(psi_upper/psi_lower)* &
              sqrt(one+micro_slope**2)/ &
              (micro_slope*(psi_upper-psi_lower))
             velsrc_micro=source_perim_factor*abs(ksrc_derived*velsrc_micro/LL)
             velsrc=velsrc_micro
            else
             print *,"microlayer_angle_source invalid"
             stop
            endif
           else
            print *,"macrolayer_size_source invalid"
            stop
           endif
          else if (microlayer_substrate_source.eq.0) then
           ! do nothing
          else
           print *,"microlayer_substrate_source invalid"
           stop
          endif
         else if (microlayer_size_source.eq.zero) then
          ! do nothing
         else
          print *,"microlayer_size_source invalid"
          stop
         endif
  
         if (microlayer_size_dest.gt.zero) then
          if ((microlayer_substrate_dest.ge.1).and. &
              (microlayer_substrate_dest.le.num_materials)) then
           if (is_rigid(microlayer_substrate_dest).ne.1) then
            print *,"is_rigid(microlayer_substrate_dest).ne.1"
            stop
           endif
           if (macrolayer_size_dest.gt.microlayer_size_dest) then
            if ((microlayer_angle_dest.gt.zero).and. &
                (microlayer_angle_dest.lt.Pi)) then
             micro_slope=tan(half*microlayer_angle_dest)
             psi_upper=half*macrolayer_size_dest/micro_slope
             psi_lower=half*microlayer_size_dest/micro_slope
             veldst_micro=( &
              half*(Tsat+fort_tempconst(microlayer_substrate_dest))-Tsat)* &
              log(psi_upper/psi_lower)* &
              sqrt(one+micro_slope**2)/ &
              (micro_slope*(psi_upper-psi_lower))
             veldst_micro=dest_perim_factor*abs(kdst_derived*veldst_micro/LL)
             veldst=veldst_micro
            else
             print *,"microlayer_angle_dest invalid"
             stop
            endif
           else
            print *,"macrolayer_size_dest invalid"
            stop
           endif
          else if (microlayer_substrate_dest.eq.0) then
           ! do nothing
          else
           print *,"microlayer_substrate_dest invalid"
           stop
          endif
         else if (microlayer_size_dest.eq.zero) then
          ! do nothing
         else
          print *,"microlayer_size_dest invalid"
          stop
         endif

         ! above: local_freezing_model==0
        else if (local_freezing_model.eq.1) then
         ! do nothing

         ! our evaporation model corresponds to "boiling"
         ! in which the temperature gradient is assumed constant
         ! in the vapor and the vapor is captured instead of tracked.
         ! The saturation temperature is currently constant, but later
         ! on, the saturation temperature should make use of the
         ! Clausius Clapeyron condition.
         ! 
         ! Interface Temperature used for solving the heat equation:
         !  1. advection (rho Y)_t + div (rho u Y) = div (rho D grad Y)
         !      rho_t Y + rho Y_t + div (rho u) Y + rho u grad Y =
         !      div (rho D grad Y)
         !      Y_t + u dot grad Y = div (rho D grad Y)/rho 
         !  2. rate of mass transfer, input: probe Temperature, mass fraction,
         !     density; output: rate of mass transfer, interface temperature,
         !     interface species
         !  3. diffusion species mass fraction, and temperature
         !  Supermesh for Temperature and species is good.
         !  Supermesh for viscous solver and pressure projection???
        else if (local_freezing_model.eq.5) then

         if (LL.gt.zero) then ! evaporation
          veldst=zero ! ignore temperature gradient in the air.
         else if (LL.lt.zero) then ! condensation
          velsrc=zero
         else if (LL.eq.zero) then
          print *,"LL invalid"
          stop
         else
          print *,"LL invalid"
          stop
         endif 
        else if (local_freezing_model.eq.6) then ! Palmore and Desjardins
         ! do nothing
        else if (local_freezing_model.eq.7) then ! cavitation
         print *,"cavitation model in development"
         stop
        else
         print *,"local_freezing_model invalid 12"
         stop
        endif

        if (distribute_from_target.eq.0) then ! default
         velsrc=velsrc/densrc_I
         veldst=veldst/densrc_I
        else if (distribute_from_target.eq.1) then
         velsrc=velsrc/dendst_I
         veldst=veldst/dendst_I
        else
         print *,"distribute_from_target invalid"
         stop
        endif

        if (for_estdt.eq.1) then
         velsrc=abs(velsrc)
         veldst=abs(veldst)
         velsrc=max(velsrc,velsrc/(one-expansion_fact))
         veldst=max(veldst,veldst/(one-expansion_fact))
        else if (for_estdt.eq.0) then
         ! do nothing
        else
         print *,"for_estdt invalid"
         stop
        endif

        velsum=velsrc+veldst
        if (velsum.gt.zero) then
         ! do nothing
        else if (velsum.le.zero) then
         velsum=zero
        else
         print *,"velsum invalid in get_vel_phasechange"
         print *,"velsum=",velsum
         print *,"velsrc=",velsrc
         print *,"veldst=",veldst
         print *,"for_estdt=",for_estdt
         print *,"expansion_fact=",expansion_fact
         print *,"densrc_I=",densrc_I
         print *,"densrc_I=",dendst_I
         print *,"distribute_from_target=",distribute_from_target
         print *,"LL=",LL
         print *,"local_freezing_model=",local_freezing_model
         print *,"ksrc_derived,kdst_derived ",ksrc_derived,kdst_derived
         print *,"ksrc_physical,kdst_physical ",ksrc_physical,kdst_physical
         print *,"DTsrc,DTdst ",DTsrc,DTdst
         print *,"dxprobe_source, dxprobe_dest ", &
                 dxprobe_source, dxprobe_dest
         stop
        endif 

        vel=velsum

      else if (local_freezing_model.eq.2) then

        if (distribute_from_target.ne.0) then
         print *,"distribute_from_target invalid"
         stop
        endif

        verb_hydrate=0
        if (((VOFsrc.ge.VOFTOL).and.(VOFdst.ge.VOFTOL)).or. &
            (for_estdt.eq.0)) then
         call HYDRATE_FORMATION_RATE(time,Cmethane_in_hydrate, &
          C_w0,Tsrc_probe,PHYDWATER,vel,Tsat,K_f,verb_hydrate)

         if (for_estdt.eq.1) then
          vel=max(vel,vel/(one-expansion_fact))
         else if (for_estdt.eq.0) then
          ! do nothing
         else
          print *,"for_estdt invalid"
          stop
         endif

        else if (((VOFsrc.le.VOFTOL).or.(VOFdst.le.VOFTOL)).and. &
                 (for_estdt.eq.1)) then
         vel=zero
        else
         print *,"VOFsrc, VOFdst, or for_estdt invalid"
         stop
        endif

      else
        print *,"local_freezing_model invalid 13"
        stop
      endif

      mdot_override=0

       ! user defined mdot=[k grad T dot n]/L
       ! k units=W/(m K)
       ! W=J/s
       ! L=J/kg
       ! J/(m s K) * (1/m) (K)/(J/kg)=J/(m^2 s)  * kg/J = kg/(m^2 s)
       ! velocity=mdot/rho=kg/(m^2 s)  / (kg/m^3)= m/s
      if (is_in_probtype_list().eq.1) then
       ! compute mdot here
       call SUB_MDOT( &
        num_materials, &
        num_species_var, &
        interface_mass_transfer_model, &
        xI, &
        ispec, &
        molar_mass, & ! 1..num_materials
        species_molar_mass, & ! 1..num_species_var+1
        im_source, &
        im_dest, &
        mdot, & ! INTENT(out)
        mdot_override, & ! INTENT(inout)
        ksrc_derived, &
        kdst_derived, &
        ksrc_physical, &
        kdst_physical, &
        Tsrc_probe, &
        Tdst_probe, &
        Tsat, &
        LL, &
        dxprobe_source, &
        dxprobe_dest)
      endif 

      if (mdot_override.eq.1) then

       if (distribute_from_target.eq.0) then 
        vel=mdot/densrc_I
       else if (distribute_from_target.eq.1) then
        vel=mdot/dendst_I
       else
        print *,"distribute_from_target invalid"
        stop
       endif

       if (for_estdt.eq.1) then
        vel=max(abs(vel),abs(vel)/(one-expansion_fact))
       else if (for_estdt.eq.0) then
        ! do nothing
       else
        print *,"for_estdt invalid"
        stop
       endif

      else if (mdot_override.eq.0) then
       ! do nothing
      else
       print *,"mdot_override invalid"
       stop
      endif

      return
      end subroutine get_vel_phasechange


      subroutine length1(z,zout)
      IMPLICIT NONE
      real(amrex_real) z(2),zout
      
      zout = sqrt(z(1)*z(1)+z(2)*z(2))
      return 

      end subroutine length1
      
      subroutine DOT1(a,b,c)
      IMPLICIT NONE
      real(amrex_real) a(2),b(2),c
      
      c = a(1)*b(1)+a(2)*b(2)
      return 

      end subroutine DOT1

      subroutine distsub(p,a,b,dist)
      IMPLICIT NONE
      real(amrex_real) t,t1,t2,t3,t4,p(2),b(2),a(2)
      real(amrex_real) dist,pb(2),pa(2),ab(2),temp(2)
      integer i
       
      do i=1,2
       ab(i)=b(i)-a(i)
       pb(i)=b(i)-p(i)
       pa(i)=a(i)-p(i)
      enddo
! one should now check if dot(pa,pb)/dot(pa,pa) is between 0 and 1
      call DOT1(pb,ab,t1) 
      call DOT1(ab,ab,t2)
      t = t1/t2       
! t =DOT(pa,pb)/DOT(ab,ab)
      if ((t.gt.0).and.(t.lt.1))  then
       call DOT1(pa,ab,t3)
       do i=1,2
! the second part of temp is the projection of p onto ab	  
        temp(i) = t2*p(i)-t1*a(i)+t3*b(i)
! temp(i) =DOT(ab,ab)*p(i)-DOT(pb,ab)*a(i)+DOT(pa,ab)*b(i)
       enddo
       call length1(temp,t4)
       dist = t4/t2
! dist = 1/DOT(ab,ab)*length(temp)
      else
       call length1(pa,t1)
       call length1(pb,t2)
       dist = min(t1,t2)
      endif

      return 
      end subroutine distsub
      
      
      subroutine setupwave(zed,t)
      IMPLICIT NONE
      real(amrex_real) o(101),kk,t,ar(101),aa,bb
      integer i
      COMPLEX*16 zed(199),A,B,A3,A4,A5,II
      COMPLEX*16 P0(101),P1(101),P2(101),Q0(101),Q1(101)
      COMPLEX*16 Q2(101),Q3(101),Q4(101),Q5(101), tang,d

      kk=20
      II=DCMPLX(cos(Pi/2),sin(Pi/2))
!     plunging breaker
      A=2.4*DCMPLX(cos(Pi/3),sin(Pi/3))
      B=2.4*DCMPLX(cos(Pi/3),sin(Pi/3))
      A3=DCMPLX(0,0.6)
      A4=DCMPLX(0,-0.6)
      A5=DCMPLX(-0.35,0)

!     spilling breaker        
!     A=DCMPLX(cos(pi/4),sin(pi/4))
!     B=DCMPLX(cos(pi/4),sin(pi/4))
!     A3=DCMPLX(0,0.5)
!     A4=DCMPLX(0,-0.5)
!     A5=DCMPLX(-0.2,0)

      do i=1,101
       ar(i)=1
!      o(i)=-2+(i-1)/(23.0+2*t)
       o(i)=-2.2+(i-1)/(25.0)
      enddo
      do i=1,101
       P0(i) = t*ar(i)
       P1(i) = t*o(i)-sqrt(kk)*(1/(2*t)-1/(6*t**3))*II*ar(i)
       P2(i) = t*o(i)**2-kk*(1/(12*t**3)-7/(90*t**5)+1/(84*t**7))*ar(i)- &
                sqrt(kk)*(1/t-1/(3*t**3))*II*o(i)

       Q0(i) = 1*ar(i)
       Q1(i) = o(i)-sqrt(kk)/2.*(1/(3*t**2)-1/(5*t**4))*II*ar(i)
       Q2(i) = o(i)**2-kk*(1/(60*t**4)-13/(630*t**6)+1/(180*t**8))*ar(i)- &
                sqrt(kk)*(1/(3*t**2)-1/(5*t**4))*II*o(i)
       Q3(i) = o(i)**3-1.5*sqrt(kk)*(1/(3*t**2)-1/(5*t**4))*II*o(i)**2+ &
                3*kk*(1/(60*t**4)-13/(630*t**6)+1/(180*t**8))*o(i)+ &
                II*kk**1.5*(1/(840*t**6)-17/(7560*t**8))*ar(i)
       Q4(i) = o(i)**4-2*sqrt(kk)*(1/(3*t**2)-1/(5*t**4))*II*o(i)**3+ &
                6*kk*(1/(60*t**4)-13/(630*t**6)+1/(180*t**8))*o(i)**2+ &
                4*II*kk**1.5*(1/(840*t**6)-17/(7560*t**8))*o(i)+ &
                2*kk**2/(30240*t**8)*ar(i) 
       Q5(i) = o(i)**5-2.5*sqrt(kk)*(1/(3*t**2)-1/(5*t**4))*II*o(i)**4+ &
                10*kk*(1/(60*t**4)-13/(630*t**6)+1/(180*t**8))*o(i)**3+ &
                10*II*kk**1.5*(1/(840*t**6)-17/(7560*t**8))*o(i)**2+ &
                10*kk**2/(30240*t**8)*o(i)

       zed(i+49) = A*(P0(i)+P2(i)/2.)*II-B*(Q0(i)+Q2(i)/2.)+A3*Q3(i)+A4*Q4(i)+ &
                A5*Q5(i)-0.03*II*Q5(i)+8.1*t*ar(i)-t**1.5*II*ar(i)

      enddo


      tang = zed(50)-zed(51)
      tang = t*DREAL(tang)+II*DIMAG(tang)
      do i = 1,49
       zed(50-i) = zed(50-i+1)+tang
       tang = DREAL(tang)+II*DIMAG(tang)/i**(0.5)
      enddo


!     aa and bb are as in aa*exp(-t^2/c)+bb
      d=zed(150)-zed(50)
      aa=DIMAG(d)
      bb=DIMAG(zed(50))
      do i=1,49
       zed(i+150) = DREAL(zed(150))-i+II*(aa*exp(-REAL(i)**2/800)+bb)
      enddo
      return 
      end subroutine setupwave

 
       ! vel is initialized with the current velocity.
      subroutine vel_freestream(x,dir,vel,time, &
       presbc_array, &
       outflow_velocity_buffer_size, & !(1,1),(2,1),(3,1),(1,2),(2,2),(3,2)
       problo,probhi)
      use global_utility_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: problo(SDIM),probhi(SDIM)
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: vel
      integer, INTENT(in) :: dir
      integer dirloc
      integer, INTENT(in) :: presbc_array(SDIM,2)
      real(amrex_real), INTENT(in) :: outflow_velocity_buffer_size(2*SDIM)
      real(amrex_real) local_buffer(2*SDIM)
      real(amrex_real) buf,dist
      integer dirbc,side,ibuf
      real(amrex_real) problen(SDIM)


      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid in vel_freestream"
       stop
      endif

      do dirloc=1,SDIM
       problen(dirloc)=probhi(dirloc)-problo(dirloc)
       if (problen(dirloc).le.zero) then
        print *,"problen(dirloc).le.zero"
        stop
       else if (problen(dirloc).gt.zero) then
        ! do nothing
       else
        print *,"problen bust"
        stop
       endif
      enddo

      do ibuf=1,2*SDIM
       local_buffer(ibuf)=outflow_velocity_buffer_size(ibuf)
      enddo

      buf=problen(1)/128.0

      if (is_in_probtype_list().eq.1) then

       call SUB_velfreestream(problen,local_buffer)

      else if (probtype.eq.201) then
       if (SDIM.eq.3) then
        dirbc=2
        side=1 
        ibuf=(side-1)*SDIM+dirbc
        if (local_buffer(ibuf).eq.zero) then
         local_buffer(ibuf)=buf
        endif
        side=2
        ibuf=(side-1)*SDIM+dirbc
        if (local_buffer(ibuf).eq.zero) then
         local_buffer(ibuf)=buf
        endif
       endif ! sdim==3
       dirbc=1
       side=1
       ibuf=(side-1)*SDIM+dirbc
       if (local_buffer(ibuf).eq.zero) then
        local_buffer(ibuf)=buf
       endif
       side=2
       ibuf=(side-1)*SDIM+dirbc
       if (local_buffer(ibuf).eq.zero) then
        local_buffer(ibuf)=buf
       endif
      else if (probtype.eq.5700) then
       buf=problen(1)/32.0
       dirbc=1
       side=2
       ibuf=(side-1)*SDIM+dirbc
       if (local_buffer(ibuf).eq.zero) then
        local_buffer(ibuf)=buf
       endif
! 3D jetbend outflow conditions
! axis_dir=2 => pressure inflow BC
      else if ((probtype.eq.53).and.(axis_dir.ne.2)) then
       buf=problen(1)/32.0
       dirbc=2
       side=1
       ibuf=(side-1)*SDIM+dirbc
       if (local_buffer(ibuf).eq.zero) then
        local_buffer(ibuf)=buf
       endif
       dirbc=1
       side=2
       ibuf=(side-1)*SDIM+dirbc
       if (local_buffer(ibuf).eq.zero) then
        local_buffer(ibuf)=buf
       endif
       dirbc=2
       ibuf=(side-1)*SDIM+dirbc
       if (local_buffer(ibuf).eq.zero) then
        local_buffer(ibuf)=buf
       endif
       dirbc=SDIM
       ibuf=(side-1)*SDIM+dirbc
       if (local_buffer(ibuf).eq.zero) then
        local_buffer(ibuf)=buf
       endif
      endif

      do dirbc=1,SDIM
       do side=1,2
        ibuf=(side-1)*SDIM+dirbc
        if (local_buffer(ibuf).eq.zero) then
         ! do nothing
        else if ((local_buffer(ibuf).gt.zero).and. &
                 (local_buffer(ibuf).le.problen(dirbc)*(one+VOFTOL))) then
      
         if (presbc_array(dirbc,side).eq.EXT_DIR) then
          if (side.eq.1) then
           dist=(x(dirbc)-problo(dirbc))
           if (dist.lt.-VOFTOL*problen(dirbc)) then
            print *,"dist.lt.-VOFTOL*problen"
            stop
           else if (dist.le.local_buffer(ibuf)) then
            if (dir+1.eq.dirbc) then
             if (vel.gt.zero) then
              vel=zero
             else if (vel.le.zero) then
              ! do nothing
             else
              print *,"vel invalid"
              stop
             endif
            else if ((dir.ge.0).and.(dir.lt.SDIM)) then
             vel=zero
            else
             print *,"dir invalid"
             stop
            endif
           else if ((dist.ge.local_buffer(ibuf)).and. &
                    (dist.le.problen(dirbc)*(one+VOFTOL))) then
            ! do nothing
           else
            print *,"dist invalid"
            stop
           endif

          else if (side.eq.2) then

           dist=(probhi(dirbc)-x(dirbc))
           if (dist.lt.-VOFTOL*problen(dirbc)) then
            print *,"dist.lt.-VOFTOL*problen"
            stop
           else if (dist.le.local_buffer(ibuf)) then
            if (dir+1.eq.dirbc) then
             if (vel.lt.zero) then
              vel=zero
             else if (vel.ge.zero) then
              ! do nothing
             else
              print *,"vel invalid"
              stop
             endif
            else if ((dir.ge.0).and.(dir.lt.SDIM)) then
             vel=zero
            else
             print *,"dir invalid"
             stop
            endif
           else if ((dist.ge.local_buffer(ibuf)).and. &
                    (dist.le.problen(dirbc)*(one+VOFTOL))) then
            ! do nothing
           else
            print *,"dist invalid"
            stop
           endif

          else
           print *,"side invalid"
           stop
          endif
         else if ((presbc_array(dirbc,side).eq.INT_DIR).or. &
                  (presbc_array(dirbc,side).eq.REFLECT_EVEN).or. &
                  (presbc_array(dirbc,side).eq.FOEXTRAP)) then
          ! do nothing
         else
          print *,"presbc_array invalid"
          print *,"probtype,dirbc,side,presbc_array ",probtype,dirbc,side, &
            presbc_array(dirbc,side)
          stop
         endif
        else 
         print *,"local_buffer invalid"
         print *,"ibuf,dirbc,side,problen,local_buffer ", &
          ibuf,dirbc,side,problen(dirbc),local_buffer(ibuf)
         print *,"probtype= ",probtype
         stop
        endif
       enddo ! side=1..2
      enddo ! dirbc=1..SDIM

      return
      end subroutine vel_freestream

subroutine initialvel(N1parm,N2parm,vel_lr,     &
 vel_lz,vel_gr, vel_gz, r1,r2,W1bar,W2bar,Re,We,RGASRWATER)

implicit none

real(amrex_real) Re,We,RGASRWATER
integer, INTENT(in)          :: N1parm
integer, INTENT(in)          :: N2parm
integer, parameter          :: nn=0
real(amrex_real) :: alpha_real  ! xblob9
real(amrex_real) :: alpha_imag  ! yblob9

integer                     :: M, j, k, ii,jj
complex*16                  :: alpha
complex*16, dimension(:,:), allocatable :: eigenv
complex*16, dimension(:), allocatable :: omega1

complex*16, dimension(0:N1parm) :: vel_lr
complex*16, dimension(0:N1parm) :: vel_lz
complex*16, dimension(0:N2parm) :: vel_gr
complex*16, dimension(0:N2parm) :: vel_gz

real(amrex_real), dimension(0:N1parm) :: r1
real(amrex_real), dimension(0:N2parm) :: r2
real(amrex_real) :: fr,fi,f
real(amrex_real)  :: omega_real,omega_imag,ki,kr
real(amrex_real), dimension(0:N1parm)          :: W1bar
real(amrex_real), dimension(0:N2parm)          :: W2bar


alpha_real=xblob9
alpha_imag=yblob9
print *,"alpha_real,alpha_imag ",alpha_real,alpha_imag
alpha=DCMPLX(alpha_real, alpha_imag)

! u, v, w, p, f
M=(N1parm+N2parm+2)*4+1

allocate(omega1(M))
allocate(eigenv(M,M))

print *,"calling tlsa M=",M
call tlsa(M, alpha, omega1, eigenv, N1parm, N2parm, nn, r1, r2,W1bar,W2bar, &
 Re,We,RGASRWATER)
print *,"done with tlsa"

k=0
do j=1,M
  kr=real(omega1(j))
  ki=aimag(omega1(j))
!  if ((abs(kr).lt.EPS8).and.(ki.gt.EPS2).and.(ki.lt.one)) then
  if ((abs(kr).lt.2.0).and.(ki.gt.zero).and.(ki.lt.one)) then
   k=k+1
   omega_real=kr
   omega_imag=ki
   fr=real(eigenv(M,j))
   fi=aimag(eigenv(M,j))
   f=fr*fr+fi*fi
   write(*,*) "the growthrate is (omega real, omega imag: ",   &
     omega_real,omega_imag

   do ii=0,N1parm
     if (ii+1.gt.M) then
      print *,"ii too big"
      stop
     endif
     kr=real(eigenv(ii+1,j))
     ki=aimag(eigenv(ii+1,j))
     vel_lr(ii)=DCMPLX((kr*fr+ki*fi)/f,(ki*fr-kr*fi)/f)
     if (ii+2*N1parm+3.gt.M) then
      print *,"ii too big"
      stop
     endif
     kr=real(eigenv(ii+2*N1parm+3,j))
     ki=aimag(eigenv(ii+2*N1parm+3,j))
     vel_lz(ii)=DCMPLX((kr*fr+ki*fi)/f,(ki*fr-kr*fi)/f)
   enddo
   jj=4*N1parm+5
   do ii=0,N2parm
     if (ii+jj.gt.M) then
      print *,"ii+jj too big"
      stop
     endif
     kr=real(eigenv(jj+ii,j))
     ki=aimag(eigenv(jj+ii,j))
     vel_gr(ii)=DCMPLX((kr*fr+ki*fi)/f,(ki*fr-kr*fi)/f)
     if (jj+2*N2parm+2+ii.gt.M) then
      print *,"index too big"
      stop
     endif
     kr=real(eigenv(jj+2*N2parm+2+ii,j))
     ki=aimag(eigenv(jj+2*N2parm+2+ii,j))
     vel_gz(ii)=DCMPLX((kr*fr+ki*fi)/f,(ki*fr-kr*fi)/f)
   enddo
  endif
enddo
if(k==0) then
  write(*,*) "no suitable eigenvalues"
  stop
endif
deallocate(omega1)
deallocate(eigenv)
return

end subroutine initialvel

subroutine normcomplexmat(A,M)
implicit none

integer M,argi,argj,i,j
complex*16 A(M,M)
real(amrex_real) magmax,minreal,maxreal,mincomplex,maxcomplex,tr,ti,magtemp

 magmax=0.0
 minreal=1.0D+20
 maxreal=-1.0D+20
 mincomplex=1.0D+20
 maxcomplex=-1.0D+20
 do i=1,M
 do j=1,M
  tr=real(A(i,j))
  ti=aimag(A(i,j))
  magtemp=sqrt(tr**2+ti**2)
  if (magtemp.gt.magmax) then
   magmax=magtemp
   argi=i
   argj=j
  endif
  if (tr.gt.maxreal) then
   maxreal=tr
  endif
  if (tr.lt.minreal) then
   minreal=tr
  endif
  if (ti.gt.maxcomplex) then
   maxcomplex=ti
  endif
  if (ti.lt.mincomplex) then
   mincomplex=ti
  endif
 enddo
 enddo
 print *,"M=",M
 print *,"magmax,argi,argj ",magmax,argi,argj
 print *,"maxreal=",maxreal
 print *,"minreal=",minreal
 print *,"maxcomplex=",maxcomplex
 print *,"mincomplex=",mincomplex

return
end subroutine


subroutine tlsa(M, alpha, omega, eigenv,N1parm,N2parm,nn,r1,r2,W1bar,W2bar, &
  Re,We,RGASRWATER)
implicit none

real(amrex_real) Re,We,RGASRWATER
integer, INTENT(in)                      :: M,N1parm,N2parm,nn
complex*16, INTENT(in)                   :: alpha
complex*16, dimension(1:M)               :: omega
complex*16, dimension(1:M,1:M)           :: eigenv
real(amrex_real), dimension(0:N1parm)              :: r1
real(amrex_real), dimension(0:N2parm)              :: r2


integer                            :: i, j
integer                            :: LWORK
integer                            :: LDA
integer                            :: LDB
integer                            :: ldvl
integer                            :: ldvr
complex*16, dimension(M)           :: alpha1
complex*16, dimension(M)           :: beta1
complex*16, dimension(M,M)         :: vr
complex*16, dimension(M,M)         :: A
complex*16, dimension(M,M)         :: B
complex*16, dimension(:),allocatable      :: WORK
real(amrex_real),dimension(:),allocatable      :: RWORK

real(amrex_real)                   :: temp,maxbeta
real(amrex_real), dimension(0:N1parm)          :: W1bar
real(amrex_real), dimension(0:N2parm)          :: W2bar

LWORK=20*M
LDA=M
LDB=M
ldvl=M
ldvr=M

print *,"before initialize N1,N2,nn,alpha ",N1parm,N2parm,nn,alpha
call initialize(M,A,B,alpha,N1parm,N2parm,nn,r1,r2,W1bar,W2bar, &
  Re,We,RGASRWATER)
print *,"after initialize"

allocate(RWORK(20*M))
allocate(WORK(LWORK))

#ifdef HAS_BLAS
print *,"before zggev M,LDA,LDB ",M,LDA,LDB
print *,"norm of A"
call normcomplexmat(A,M)
print *,"norm of B"
call normcomplexmat(B,M)
call zggev(jobvl,jobvr,M,A,LDA,B,LDB,alpha1,beta1,vl,ldvl,vr,ldvr,WORK, &
          LWORK,RWORK,INFO)
print *,"after zggev"
#endif
#ifndef HAS_BLAS
print *,"this routine can only be called if USEBLAS=TRUE"
stop
#endif

maxbeta=0.0d0
do i=1,M
 temp=sqrt(aimag(beta1(i))**2+real(beta1(i))**2)
 if (temp.gt.maxbeta) then
  maxbeta=temp
 endif
enddo
do i=1,M
  temp=sqrt(aimag(beta1(i))**2+real(beta1(i))**2)
  if(temp<=(EPS10)*maxbeta) then
     omega(i)=1.0D+30
  else
     omega(i)=alpha1(i)/beta1(i)
  endif
  do j=1,M
     eigenv(i,j)=vr(i,j)
  enddo
enddo
deallocate(WORK)
deallocate(RWORK)
end subroutine tlsa

subroutine initialize(M,A,B,alpha,N1parm,N2parm,nn,r1,r2,W1bar,W2bar, &
  ReIN,WeIN,RGASRWATER)
use global_utility_module
IMPLICIT NONE

real(amrex_real) ReIN,WeIN,RGASRWATER
integer, INTENT(in)                :: M,N1parm,N2parm,nn
complex*16, INTENT(in)             :: alpha 
complex*16, dimension(0:M-1,0:M-1) :: A
complex*16, dimension(0:M-1,0:M-1) :: B

real(amrex_real), parameter :: PI_lsa = 3.1415926535898d0
!lsaproblem_type=0 Rayleigh capillary test
!lsaproblem_type=1 Co-flowing jet test with Lin's basic flow
integer :: lsaproblem_type

real(amrex_real)  ::  N, Q, l, Re, We, Fr

integer                            :: i, j, k
real(amrex_real)                             :: alphai,alphar,alphai2,alphar2,alphari
real(amrex_real), dimension(0:N1parm)         :: C1
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: E1
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: ER1
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: F1
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: ERE1
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: E1alt
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: F1alt
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: ERE1alt
real(amrex_real), dimension(0:N2parm)         :: C2
real(amrex_real), dimension(0:N2parm,0:N2parm)    :: E2
real(amrex_real), dimension(0:N2parm,0:N2parm)    :: ER2
real(amrex_real), dimension(0:N2parm,0:N2parm)    :: F2
real(amrex_real), dimension(0:N2parm,0:N2parm)    :: ERE2
real(amrex_real), dimension(0:M-1,0:M-1)  :: Ar
real(amrex_real), dimension(0:M-1,0:M-1)  :: Ai
real(amrex_real), dimension(0:M-1,0:M-1)  :: Br
real(amrex_real), dimension(0:M-1,0:M-1)  :: Bi
 
real(amrex_real)          :: NQ,Q1,NQ2,Q2
real(amrex_real), dimension(0:N1parm)          :: eta
real(amrex_real), dimension(0:N2parm)          :: xi
real(amrex_real), dimension(0:N1parm)          :: r1
real(amrex_real), dimension(0:N2parm)          :: r2
real(amrex_real), dimension(0:N1parm)          :: W1bar,DW1bar
real(amrex_real), dimension(0:N1parm)          :: DDW1bar
real(amrex_real), dimension(0:N2parm)          :: W2bar,DW2bar
real(amrex_real), dimension(0:N2parm)          :: DDW2bar

real(amrex_real) :: rj2,rj
  
integer IU,IV,IW,IP
integer IU2,IV2,IW2,IP2
integer IFVAR
integer ICONT,IUMOM,IVMOM,IWMOM
integer ICONT2,IUMOM2,IVMOM2,IWMOM2
integer IUNOSLIP,IVNOSLIP,IWNOSLIP
integer IUJUMP,IVJUMP,IWJUMP,IPJUMP,IFEQN,ITAN1parmEQN,ITAN2parmEQN
integer IEXTRAEQN
real(amrex_real) :: coeff,coeff2,coeff3,coeff4,coeff5,coeff6,coeff7,RR


lsaproblem_type=NINT(radblob9)
print *,"lsaproblem_type=",lsaproblem_type


N=get_user_viscconst(2,fort_denconst(2),fort_tempconst(2))/ &
  get_user_viscconst(1,fort_denconst(1),fort_tempconst(1))

Q=fort_denconst(2)/fort_denconst(1)
l=RGASRWATER
Re=ReIN
We=WeIN
Fr=1.0E12

print *,"N,Q,l,Re,We,Fr ",N,Q,l,Re,We,Fr

do k=0,N1parm
   eta(k)=-dcos(PI_lsa*DBLE(k)/DBLE(N1parm))
   r1(k)=(1.0d0+eta(k))/2.0d0
   if(lsaproblem_type==0) then
      W1bar(k)=0.0d0
      DW1bar(k)=-0.0d0
      DDW1bar(k)=-0.0d0
   elseif(lsaproblem_type==1) then
      RR=Re/Fr
      coeff=(1.0d0-Q)*RR*(2.0d0*log(l)+1.0d0-l*l)/4.0d0
      coeff2=(N-coeff)/(N-1.0d0+l*l)
      W1bar(k)=1.0d0-coeff2*(r1(k)**2)
      DW1bar(k)=-2.0d0*coeff2*r1(k)
      DDW1bar(k)=-2.0d0*coeff2
      W1bar(k)=-W1bar(k)
      DW1bar(k)=-DW1bar(k)
      DDW1bar(k)=-DDW1bar(k)
   endif

enddo
do k=0,N2parm
   xi(k)=-dcos(PI_lsa*DBLE(k)/DBLE(N2parm))
   r2(k)=1.0d0+(xi(k)+1.0d0)*(l-1.0d0)/2.0d0

   if(lsaproblem_type==0) then
      W2bar(k)=0.0d0
      DW2bar(k)=-0.0d0
      DDW2bar(k)=-0.0d0
   elseif(lsaproblem_type==1) then
      RR=Re/Fr
      coeff=(1.0d0-Q)*RR*(2.0d0*log(l)+1.0d0-l*l)/4.0d0
      coeff2=(N-coeff)/(N*(N-1.0d0+l*l))
      coeff3=(1.0d0-Q)*RR/(4.0d0*N)
      coeff4=l*l*(coeff3-coeff2)
      coeff5=(coeff2-coeff3)
      coeff6=-2.0d0*coeff3
      coeff7=coeff4+coeff6*log(l)
      W2bar(k)=-coeff5*(r2(k))**2+coeff6*log(r2(k))-coeff7
      DW2bar(k)=-2.0d0*coeff5*r2(k)+coeff6/r2(k)
      DDW2bar(k)=-2.0d0*coeff5-coeff6/(r2(k)**2)
      W2bar(k)=-W2bar(k)
      DW2bar(k)=-DW2bar(k)
      DDW2bar(k)=-DDW2bar(k)
   endif
enddo

do i=1,N1parm-1
   C1(i)=1.0d0
   enddo
C1(0)=2.0d0
C1(N1parm)=2.0d0
do j=0,N1parm
   do k=0,N1parm
    if (j.ne.k) then
     E1(j,k)=C1(j)/(C1(k)*(eta(j)-eta(k)))
     if(mod(j+k,2).NE.0) then
      E1(j,k)=-E1(j,k)
     endif
    endif
   enddo
enddo
do j=1,N1parm-1
   E1(j,j)=-eta(j)/(2.0d0*(1.0d0-eta(j)*eta(j)))
enddo
E1(0,0)=-(2.0d0*N1parm*N1parm+1.0d0)/6.0d0
E1(N1parm,N1parm)=-E1(0,0)

do j=0,N1parm
do k=0,N1parm
 E1(j,k)=E1(j,k)*2.0d0
enddo
enddo

do j=0,N1parm
 do k=0,N1parm
  F1(j,k)=0.0d0
  do i=0,N1parm
    F1(j,k)=F1(j,k)+E1(j,i)*E1(i,k)
  enddo
 enddo
enddo

! multiply the columns by r1
do j=0,N1parm
do k=0,N1parm
 ER1(j,k)=E1(j,k)*r1(k)
enddo
enddo
    
do j=0,N1parm
 do k=0,N1parm
  ERE1(j,k)=0.0d0
  do i=0,N1parm
    ERE1(j,k)=ERE1(j,k)+ER1(j,i)*E1(i,k)
  enddo
 enddo
enddo
   
do j=0,N1parm
 do k=0,N1parm
  if (j.eq.0) then
   E1alt(j,k)=E1(j,k)
  else
   E1alt(j,k)=E1(j,k)*r1(k)/r1(j)
  endif
 enddo
 if (j.ne.0) then
  E1alt(j,j)=E1alt(j,j)-1.0d0/r1(j)
 endif
enddo
    
do j=0,N1parm
 do k=0,N1parm
  if (j.eq.0) then
   F1alt(j,k)=F1(j,k)
  else
   F1alt(j,k)=(F1(j,k)*r1(k)-2.0d0*E1alt(j,k))/r1(j)
  endif
 enddo
enddo

do j=0,N1parm
 do k=0,N1parm
  ERE1alt(j,k)=E1alt(j,k)+r1(j)*F1alt(j,k)
 enddo
enddo
    
do i=1,N2parm-1
   C2(i)=1.0d0
enddo
C2(0)=2.0d0
C2(N2parm)=2.0d0
do j=0,N2parm
 do k=0,N2parm
  if (j.ne.k) then
   E2(j,k)=C2(j)/(C2(k)*(xi(j)-xi(k)))
   if(mod(j+k,2).NE.0) then
    E2(j,k)=-E2(j,k)
   endif
  endif
 enddo
enddo
do j=1,N2parm-1
   E2(j,j)=-xi(j)/(2.0d0*(1.0d0-xi(j)*xi(j)))
enddo
E2(0,0)=-(2.0d0*N2parm*N2parm+1.0d0)/6.0d0
E2(N2parm,N2parm)=-E2(0,0)

do j=0,N2parm
do k=0,N2parm
 E2(j,k)=E2(j,k)*2.0d0/(l-1.0d0)
enddo
enddo

do j=0,N2parm
 do k=0,N2parm
  F2(j,k)=0.0d0
  do i=0,N2parm
   F2(j,k)=F2(j,k)+E2(j,i)*E2(i,k)
  enddo
 enddo
enddo

! multiply the columns by r1
do j=0,N2parm
do k=0,N2parm
 ER2(j,k)=E2(j,k)*r2(k)
enddo
enddo

do j=0,N2parm
 do k=0,N2parm
  ERE2(j,k)=0.0d0
  do i=0,N2parm
    ERE2(j,k)=ERE2(j,k)+ER2(j,i)*E2(i,k)
  enddo
 enddo
enddo


alphai=aimag(alpha)
alphar=real(alpha)
alphai2=alphai*alphai
alphar2=alphar*alphar
alphari=alphai*alphar

!write(*,*) alphar, alphai, alphari, alphar2, alphai2

! initialize Ar and Ai
do j=0, M-1
   do k=0,M-1
      Ar(j,k)=0.0d0
      Ai(j,k)=0.0d0
   enddo
enddo

IU=0
IV=IU+N1parm+1
IW=IV+N1parm+1
IP=IW+N1parm+1

IU2=IP+N1parm+1
IV2=IU2+N2parm+1
IW2=IV2+N2parm+1
IP2=IW2+N2parm+1

IFVAR=IP2+N2parm+1

if (IFVAR+1.ne.M) then
 print *,"number of variables invalid"
 stop
endif

ICONT=0
IUMOM=ICONT+N1parm
IVMOM=IUMOM+N1parm
IWMOM=IVMOM+N1parm

ICONT2=IWMOM+N1parm
IUMOM2=ICONT2+N2parm-1
IVMOM2=IUMOM2+N2parm
IWMOM2=IVMOM2+N2parm-1
   
IUNOSLIP=IWMOM2+N2parm-1
IVNOSLIP=IUNOSLIP+1
IWNOSLIP=IVNOSLIP+1

IUJUMP=IWNOSLIP+1
IVJUMP=IUJUMP+1
IWJUMP=IVJUMP+1
IPJUMP=IWJUMP+1
IFEQN=IPJUMP+1
ITAN1parmEQN=IFEQN+1
ITAN2parmEQN=ITAN1parmEQN+1
   
IEXTRAEQN=ITAN2parmEQN+1
   
if (IEXTRAEQN+2.ne.M) then
 print *,"number of equations invalid"
 stop
endif
   
NQ=1.0
Q1=1.0
  
NQ2=N/Q
Q2=Q
! at centerline:
! if nn=0, u=v=0  p and w are finite
! if nn=1, 2u'+iv'=0, u+iv=0  w=p=0
! o.t. u=v=w=p=0
j=0
Ar(j+ICONT,j+IP)=1.0d0

do j=1,N1parm-1
 rj=r1(j)
 do k=0,N1parm
  Ar(j+ICONT,k+IU)=ER1(j,k)/rj
 enddo
 Ai(j+ICONT,j+IV)=nn/rj
 Ar(j+ICONT,j+IW)=-alphai
 Ai(j+ICONT,j+IW)=alphar

 do k=0,N1parm
   Ar(j+ICONT,k+IU)=Ar(j+ICONT,k+IU)*rj
 enddo
 Ai(j+ICONT,j+IV)=Ai(j+ICONT,j+IV)*rj
 Ar(j+ICONT,j+IW)=Ar(j+ICONT,j+IW)*rj
 Ai(j+ICONT,j+IW)=Ai(j+ICONT,j+IW)*rj
enddo

j=0
if (nn.ne.1) then 
 Ar(j+IUMOM,j+IU)=1.0d0
else if (nn.eq.1) then
 do k=0,N1parm
  Ar(j+IUMOM,k+IU)=2.0d0*E1(0,k)
  Ai(j+IUMOM,k+IV)=E1(0,k)
 enddo
endif

do j=1,N1parm-1
 rj=r1(j)
 rj2=rj**2
 do k=0,N1parm
   Ar(j+IUMOM,k+IU)=Q1*NQ*ERE1(j,k)/rj
 enddo
 do k=0,N1parm
  Ar(j+IUMOM,k+IP)=-Re*E1alt(j,k)
 enddo
 Ar(j+IUMOM,j+IU)=Ar(j+IUMOM,j+IU)-Q1*(NQ*((nn**2+1.0d0)/rj2+     &
            alphar2-alphai2)-Re*alphai*W1bar(j))
 Ai(j+IUMOM,j+IU)=-Q1*(Re*alphar*W1bar(j)+2.0d0*NQ*alphari)
 Ai(j+IUMOM,j+IV)=-2.0*nn*NQ*Q1/rj2
 do k=0,N1parm
   Ar(j+IUMOM,k+IU)=Ar(j+IUMOM,k+IU)*rj2
   Ar(j+IUMOM,k+IP)=Ar(j+IUMOM,k+IP)*rj2
 enddo
 Ai(j+IUMOM,j+IU)=Ai(j+IUMOM,j+IU)*rj2
 Ai(j+IUMOM,j+IV)=Ai(j+IUMOM,j+IV)*rj2
enddo

j=0
if (nn.ne.1) then
 Ar(j+IVMOM,j+IV)=1.0d0
else
 Ar(j+IVMOM,j+IU)=1.0d0
 Ai(j+IVMOM,j+IV)=1.0d0
endif

do j=1,N1parm-1
 rj=r1(j)
 rj2=rj**2
 do k=0,N1parm
   Ar(j+IVMOM,k+IV)=NQ*Q1*ERE1(j,k)/rj
 enddo
 Ar(j+IVMOM,j+IV)=Ar(j+IVMOM,j+IV)-Q1*(NQ*((nn**2+1.0d0)/rj2+     &
            alphar2-alphai2)-Re*alphai*W1bar(j))
 Ai(j+IVMOM,j+IV)=-Q1*(Re*alphar*W1bar(j)+2.0d0*NQ*alphari)
 Ai(j+IVMOM,j+IU)=Q1*2.0d0*nn*NQ/rj2
 Ai(j+IVMOM,j+IP)=-Re*nn/rj

 do k=0,N1parm
   Ar(j+IVMOM,k+IV)=Ar(j+IVMOM,k+IV)*rj2
 enddo
 Ai(j+IVMOM,j+IV)=Ai(j+IVMOM,j+IV)*rj2
 Ai(j+IVMOM,j+IU)=Ai(j+IVMOM,j+IU)*rj2
 Ai(j+IVMOM,j+IP)=Ai(j+IVMOM,j+IP)*rj2

enddo

j=0
Ar(j+IWMOM,j+IW)=1.0d0

do j=1,N1parm-1
 rj=r1(j)
 rj2=rj**2
 do k=0,N1parm
   Ar(j+IWMOM,k+IW)=NQ*Q1*ERE1alt(j,k)/rj
 enddo
 Ar(j+IWMOM,j+IW)=Ar(j+IWMOM,j+IW)-Q1*(NQ*((nn**2)/rj2+     &
            alphar2-alphai2)-Re*alphai*W1bar(j))
 Ai(j+IWMOM,j+IW)=-Q1*(Re*alphar*W1bar(j)+2.0d0*NQ*alphari)
 Ar(j+IWMOM,j+IU)=-Re*Q1*DW1bar(j)
 Ar(j+IWMOM,j+IP)=Re*alphai
 Ai(j+IWMOM,j+IP)=-Re*alphar
 
 do k=0,N1parm
   Ar(j+IWMOM,k+IW)=Ar(j+IWMOM,k+IW)*rj2
 enddo
 Ai(j+IWMOM,j+IW)=Ai(j+IWMOM,j+IW)*rj2
 Ar(j+IWMOM,j+IU)=Ar(j+IWMOM,j+IU)*rj2
 Ar(j+IWMOM,j+IP)=Ar(j+IWMOM,j+IP)*rj2
 Ai(j+IWMOM,j+IP)=Ai(j+IWMOM,j+IP)*rj2
 
enddo
  
 
do j=1,N2parm-1
 rj=r2(j)
 do k=0,N2parm
   Ar(j+ICONT2-1,k+IU2)=ER2(j,k)/rj
 enddo
 Ai(j+ICONT2-1,j+IV2)=nn/rj
 Ar(j+ICONT2-1,j+IW2)=-alphai
 Ai(j+ICONT2-1,j+IW2)=alphar
 
 do k=0,N2parm
   Ar(j+ICONT2-1,k+IU2)=Ar(j+ICONT2-1,k+IU2)*rj
 enddo
 Ai(j+ICONT2-1,j+IV2)=Ai(j+ICONT2-1,j+IV2)*rj
 Ar(j+ICONT2-1,j+IW2)=Ar(j+ICONT2-1,j+IW2)*rj
 Ai(j+ICONT2-1,j+IW2)=Ai(j+ICONT2-1,j+IW2)*rj

enddo
 
do j=1,N2parm-1
 rj=r2(j)
 rj2=rj**2
 do k=0,N2parm
   Ar(j+IUMOM2-1,k+IU2)=NQ2*Q2*ERE2(j,k)/rj
 enddo
 do k=0,N2parm
  Ar(j+IUMOM2-1,k+IP2)=-Re*E2(j,k)
 enddo
 Ar(j+IUMOM2-1,j+IU2)=Ar(j+IUMOM2-1,j+IU2)-Q2*(NQ2*((nn**2+1.0d0)/rj2+     &
            alphar2-alphai2)-Re*alphai*W2bar(j))
 Ai(j+IUMOM2-1,j+IU2)=-Q2*(Re*alphar*W2bar(j)+2.0d0*NQ2*alphari)
 Ai(j+IUMOM2-1,j+IV2)=-2.0*nn*Q2*NQ2/rj2
 
 do k=0,N2parm
   Ar(j+IUMOM2-1,k+IU2)=Ar(j+IUMOM2-1,k+IU2)*rj2
   Ar(j+IUMOM2-1,k+IP2)=Ar(j+IUMOM2-1,k+IP2)*rj2
 enddo
 Ai(j+IUMOM2-1,j+IU2)=Ai(j+IUMOM2-1,j+IU2)*rj2
 Ai(j+IUMOM2-1,j+IV2)=Ai(j+IUMOM2-1,j+IV2)*rj2
enddo

j=N2parm
rj=r2(j)
rj2=rj**2
do k=0,N2parm
  Ar(j+IUMOM2-1,k+IU2)=NQ2*Q2*ERE2(j,k)/rj
enddo
do k=0,N2parm
 Ar(j+IUMOM2-1,k+IP2)=-Re*E2(j,k)
enddo
Ar(j+IUMOM2-1,j+IU2)=Ar(j+IUMOM2-1,j+IU2)-Q2*(NQ2*((nn**2+1.0d0)/rj2+     &
           alphar2-alphai2)-Re*alphai*W2bar(j))
Ai(j+IUMOM2-1,j+IU2)=-Q2*(Re*alphar*W2bar(j)+2.0d0*NQ2*alphari)
Ai(j+IUMOM2-1,j+IV2)=-2.0*nn*Q2*NQ2/rj2

do k=0,N2parm
  Ar(j+IUMOM2-1,k+IU2)=Ar(j+IUMOM2-1,k+IU2)*rj2
  Ar(j+IUMOM2-1,k+IP2)=Ar(j+IUMOM2-1,k+IP2)*rj2
enddo
Ai(j+IUMOM2-1,j+IU2)=Ai(j+IUMOM2-1,j+IU2)*rj2
Ai(j+IUMOM2-1,j+IV2)=Ai(j+IUMOM2-1,j+IV2)*rj2

do j=1,N2parm-1
 rj=r2(j)
 rj2=rj**2
 do k=0,N2parm
   Ar(j+IVMOM2-1,k+IV2)=NQ2*Q2*ERE2(j,k)/rj
 enddo
 Ar(j+IVMOM2-1,j+IV2)=Ar(j+IVMOM2-1,j+IV2)-Q2*(NQ2*((nn**2+1.0d0)/rj2+     &
            alphar2-alphai2)-Re*alphai*W2bar(j))
 Ai(j+IVMOM2-1,j+IV2)=-Q2*(Re*alphar*W2bar(j)+2.0d0*NQ2*alphari)
 Ai(j+IVMOM2-1,j+IU2)=2.0*nn*NQ2*Q2/rj2
 Ai(j+IVMOM2-1,j+IP2)=-Re*nn/rj

 do k=0,N2parm
   Ar(j+IVMOM2-1,k+IV2)=Ar(j+IVMOM2-1,k+IV2)*rj2
 enddo
 Ai(j+IVMOM2-1,j+IV2)=Ai(j+IVMOM2-1,j+IV2)*rj2
 Ai(j+IVMOM2-1,j+IU2)=Ai(j+IVMOM2-1,j+IU2)*rj2
 Ai(j+IVMOM2-1,j+IP2)=Ai(j+IVMOM2-1,j+IP2)*rj2
enddo

do j=1,N2parm-1
 rj=r2(j)
 rj2=rj**2
 do k=0,N2parm
   Ar(j+IWMOM2-1,k+IW2)=NQ2*Q2*ERE2(j,k)/rj
 enddo
 Ar(j+IWMOM2-1,j+IW2)=Ar(j+IWMOM2-1,j+IW2)-Q2*(NQ2*((nn**2)/rj2+     &
            alphar2-alphai2)-Re*alphai*W2bar(j))
 Ai(j+IWMOM2-1,j+IW2)=-Q2*(Re*alphar*W2bar(j)+2.0d0*NQ2*alphari)
 Ar(j+IWMOM2-1,j+IU2)=-Re*Q2*DW2bar(j)
 Ar(j+IWMOM2-1,j+IP2)=Re*alphai
 Ai(j+IWMOM2-1,j+IP2)=-Re*alphar

 do k=0,N2parm
   Ar(j+IWMOM2-1,k+IW2)=Ar(j+IWMOM2-1,k+IW2)*rj2
 enddo
 Ai(j+IWMOM2-1,j+IW2)=Ai(j+IWMOM2-1,j+IW2)*rj2
 Ar(j+IWMOM2-1,j+IU2)=Ar(j+IWMOM2-1,j+IU2)*rj2
 Ar(j+IWMOM2-1,j+IP2)=Ar(j+IWMOM2-1,j+IP2)*rj2
 Ai(j+IWMOM2-1,j+IP2)=Ai(j+IWMOM2-1,j+IP2)*rj2

enddo


Ar(IUNOSLIP,IU2+N2parm)=1.0d0
Ar(IVNOSLIP,IV2+N2parm)=1.0d0
Ar(IWNOSLIP,IW2+N2parm)=1.0d0

Ar(IUJUMP,IU+N1parm)=1.0
Ar(IUJUMP,IU2)=-1.0

Ar(IVJUMP,IV+N1parm)=1.0
Ar(IVJUMP,IV2)=-1.0
Ar(IVJUMP,IFVAR)=0.0d0

Ar(IWJUMP,IW+N1parm)=1.0
Ar(IWJUMP,IW2)=-1.0
Ar(IWJUMP,IFVAR)=DW1bar(N1parm)-DW2bar(0)
  
Ar(IFEQN,IFVAR)=-alphai*W1bar(N1parm)
Ai(IFEQN,IFVAR)=alphar*W1bar(N1parm)
Ar(IFEQN,IU+N1parm)=-1.0

Ar(IPJUMP,IP+N1parm)=1.0
Ar(IPJUMP,IP2)=-1.0
do k=0,N1parm
 Ar(IPJUMP,IU+k)=-(2.0/Re)*E1(N1parm,k)
enddo
do k=0,N2parm
 Ar(IPJUMP,IU2+k)=(2.0/Re)*N*E2(0,k)
enddo
Ar(IPJUMP,IFVAR)=(1.0/We)*(1.0-nn*nn-alphar2+alphai2)
Ai(IPJUMP,IFVAR)=(1.0/We)*(-2.0d0*alphari)
 
do k=0,N1parm
 Ar(ITAN1parmEQN,IV+k)=E1(N1parm,k)
enddo
do k=0,N2parm
 Ar(ITAN1parmEQN,IV2+k)=-N*E2(0,k)
enddo
Ar(ITAN1parmEQN,IV+N1parm)=Ar(ITAN1parmEQN,IV+N1parm)-1.0d0
Ar(ITAN1parmEQN,IV2)=Ar(ITAN1parmEQN,IV2)+N*1.0d0
Ai(ITAN1parmEQN,IU+N1parm)=nn
Ai(ITAN1parmEQN,IU2)=-N*nn
Ar(ITAN1parmEQN,IFVAR)=0.0d0
 
Ar(ITAN2parmEQN,IU+N1parm)=-alphai
Ai(ITAN2parmEQN,IU+N1parm)=alphar
Ar(ITAN2parmEQN,IU2)=N*alphai
Ai(ITAN2parmEQN,IU2)=-N*alphar
do k=0,N1parm
 Ar(ITAN2parmEQN,IW+k)=E1alt(N1parm,k)
enddo
do k=0,N2parm
 Ar(ITAN2parmEQN,IW2+k)=-N*E2(0,k)
enddo
Ar(ITAN2parmEQN,IFVAR)=(DDW1bar(N1parm)-N*DDW2bar(0))
 
j=N1parm
rj=r1(j)
do k=0,N1parm
   Ar(IEXTRAEQN,k+IU)=ER1(j,k)/rj
enddo
Ai(IEXTRAEQN,j+IV)=nn/rj
Ai(IEXTRAEQN,j+IW)=alphar
Ar(IEXTRAEQN,j+IW)=-alphai
do k=0,N1parm
   Ar(IEXTRAEQN,k+IU)=Ar(IEXTRAEQN,k+IU)*rj
enddo
Ai(IEXTRAEQN,j+IV)=Ai(IEXTRAEQN,j+IV)*rj
Ar(IEXTRAEQN,j+IW)=Ar(IEXTRAEQN,j+IW)*rj
Ai(IEXTRAEQN,j+IW)=Ai(IEXTRAEQN,j+IW)*rj
j=0
rj=r2(j)
do k=0,N2parm
   Ar(IEXTRAEQN+1,k+IU2)=ER2(j,k)/rj
enddo
Ai(IEXTRAEQN+1,j+IV2)=nn/rj
Ai(IEXTRAEQN+1,j+IW2)=alphar
Ar(IEXTRAEQN+1,j+IW2)=-alphai

do k=0,N2parm
   Ar(IEXTRAEQN+1,k+IU2)=Ar(IEXTRAEQN+1,k+IU2)*rj
enddo
Ai(IEXTRAEQN+1,j+IV2)=Ai(IEXTRAEQN+1,j+IV2)*rj
Ar(IEXTRAEQN+1,j+IW2)=Ar(IEXTRAEQN+1,j+IW2)*rj
Ai(IEXTRAEQN+1,j+IW2)=Ai(IEXTRAEQN+1,j+IW2)*rj

do i=0,M-1
   do j=0,M-1
      A(j,i)=DCMPLX(Ar(j,i),Ai(j,i))
   enddo
enddo

do j=0, M-1
   do i=0, M-1
      Br(i,j)=0.0d0
      Bi(i,j)=0.0d0
   enddo
enddo

! opposite signs since B represents terms moved to RHS
do j=1,N1parm-1
 rj=r1(j)
 rj2=rj**2
 Bi(IUMOM+j,IU+j)=-Re*Q1
 Bi(IUMOM+j,IU+j)=Bi(IUMOM+j,IU+j)*rj2
enddo
do j=1,N1parm-1
 rj=r1(j)
 rj2=rj**2
 Bi(IVMOM+j,IV+j)=-Re*Q1
 Bi(IVMOM+j,IV+j)=Bi(IVMOM+j,IV+j)*rj2
enddo
do j=1,N1parm-1
 rj=r1(j)
 rj2=rj**2
 Bi(IWMOM+j,IW+j)=-Re*Q1
 Bi(IWMOM+j,IW+j)=Bi(IWMOM+j,IW+j)*rj2
enddo

do j=1,N2parm
 rj=r2(j)
 rj2=rj**2
 Bi(IUMOM2+j-1,IU2+j)=-Re*Q2
 Bi(IUMOM2+j-1,IU2+j)=Bi(IUMOM2+j-1,IU2+j)*rj2
enddo
do j=1,N2parm-1
 rj=r2(j)
 rj2=rj**2
 Bi(IVMOM2+j-1,IV2+j)=-Re*Q2
 Bi(IVMOM2+j-1,IV2+j)=Bi(IVMOM2+j-1,IV2+j)*rj2
enddo
do j=1,N2parm-1
 rj=r2(j)
 rj2=rj**2
 Bi(IWMOM2+j-1,IW2+j)=-Re*Q2
 Bi(IWMOM2+j-1,IW2+j)=Bi(IWMOM2+j-1,IW2+j)*rj2
enddo
Bi(IFEQN,IFVAR)=1.0d0

do i=0,M-1
   do j=0,M-1
      B(j,i)=DCMPLX(Br(j,i),Bi(j,i))
   enddo
enddo
   
return
end subroutine initialize

! RGASRWATER=ratio of domain width to liquid domain width
subroutine tlsa2d(M, alpha, omega, eigenv,N1parm,N2parm, &
  r1,r2,Re,We,RGASRWATER,density_ratio,viscosity_ratio)
implicit none

real(amrex_real) density_ratio,viscosity_ratio
real(amrex_real) Re,We,RGASRWATER
integer, INTENT(in)                      :: M,N1parm,N2parm
complex*16, INTENT(in)                   :: alpha
complex*16, dimension(1:M)               :: omega
complex*16, dimension(1:M,1:M)           :: eigenv
real(amrex_real), dimension(0:N1parm)              :: r1
real(amrex_real), dimension(0:N2parm)              :: r2


integer                            :: i, j
integer                            :: LWORK
integer                            :: LDA
integer                            :: LDB
integer                            :: ldvl
integer                            :: ldvr
complex*16, dimension(M)           :: alpha1
complex*16, dimension(M)           :: beta1
complex*16, dimension(M,M)         :: vr
complex*16, dimension(M,M)         :: A
complex*16, dimension(M,M)         :: B
complex*16, dimension(:),allocatable      :: WORK
real(amrex_real),dimension(:),allocatable      :: RWORK

real(amrex_real)                   :: temp,maxbeta

LWORK=20*M
LDA=M
LDB=M
ldvl=M
ldvr=M

print *,"before initialize2d N1,N2,alpha ",N1parm,N2parm,alpha
if ((N1parm+N2parm+2)*3+1-6.ne.M) then
 print *,"M invalid in tlsa2d"
 stop
endif

call initialize2d(M,A,B,alpha,N1parm,N2parm,r1,r2, &
  Re,We,RGASRWATER,density_ratio,viscosity_ratio)
print *,"after initialize2d"

allocate(RWORK(20*M))
allocate(WORK(LWORK))

#ifdef HAS_BLAS
print *,"before zggev M,LDA,LDB ",M,LDA,LDB
print *,"norm of A"
call normcomplexmat(A,M)
print *,"norm of B"
call normcomplexmat(B,M)
call zggev(jobvl,jobvr,M,A,LDA,B,LDB,alpha1,beta1,vl,ldvl,vr,ldvr,WORK, &
          LWORK,RWORK,INFO)
print *,"after zggev"
#endif
#ifndef HAS_BLAS
print *,"this routine can only be called if USEBLAS=TRUE"
stop
#endif

maxbeta=0.0d0
do i=1,M
 temp=sqrt(aimag(beta1(i))**2+real(beta1(i))**2)
 if (temp.gt.maxbeta) then
  maxbeta=temp
 endif
enddo
do i=1,M
  temp=sqrt(aimag(beta1(i))**2+real(beta1(i))**2)
  if(temp<=(EPS10)*maxbeta) then
     omega(i)=1.0D+30
  else
     omega(i)=alpha1(i)/beta1(i)
  endif
  do j=1,M
     eigenv(i,j)=vr(i,j)
  enddo
enddo
deallocate(WORK)
deallocate(RWORK)
end subroutine tlsa2d

subroutine initialize2d(M,A,B,alpha,N1parm,N2parm,r1,r2, &
  ReIN,WeIN,RGASRWATER,density_ratio,viscosity_ratio)
IMPLICIT NONE

real(amrex_real) density_ratio,viscosity_ratio
real(amrex_real) ReIN,WeIN,RGASRWATER
integer, INTENT(in)                :: M,N1parm,N2parm
complex*16, INTENT(in)             :: alpha 
complex*16, dimension(0:M-1,0:M-1) :: A
complex*16, dimension(0:M-1,0:M-1) :: B

real(amrex_real), parameter :: PI_lsa = 3.1415926535898d0

real(amrex_real)  ::  N, Q, l, Re, We

integer :: i, j, k
real(amrex_real)  :: alphai,alphar,alphai2,alphar2,alphari
real(amrex_real), dimension(0:N1parm)         :: C1
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: E1
real(amrex_real), dimension(0:N1parm,0:N1parm)    :: F1
real(amrex_real), dimension(0:N2parm)         :: C2
real(amrex_real), dimension(0:N2parm,0:N2parm)    :: E2
real(amrex_real), dimension(0:N2parm,0:N2parm)    :: F2
real(amrex_real), dimension(0:M-1,0:M-1)  :: Ar
real(amrex_real), dimension(0:M-1,0:M-1)  :: Ai
real(amrex_real), dimension(0:M-1,0:M-1)  :: Br
real(amrex_real), dimension(0:M-1,0:M-1)  :: Bi
 
real(amrex_real)          :: NQ,Q1,NQ2,Q2
real(amrex_real), dimension(0:N1parm)          :: eta
real(amrex_real), dimension(0:N2parm)          :: xi
real(amrex_real), dimension(0:N1parm)          :: r1
real(amrex_real), dimension(0:N2parm)          :: r2

! no more IV and IV2 
integer IU,IW,IP
integer IU2,IW2,IP2
integer IFVAR
! no more IVMOM and IVMOM2
integer ICONT,IUMOM,IWMOM
integer ICONT2,IUMOM2,IWMOM2
! no more NOSLIP,JUMP for u,v,w
! no more ITAN2parmEQN
integer IPJUMP,IFEQN,ITAN1parmEQN
integer IEXTRAEQN


! u,w,p,F
if ((N1parm+N2parm+2)*3+1-6.ne.M) then
 print *,"M invalid in initialize2d"
 stop
endif

! N=fort_viscconst(2)/fort_viscconst(1)
! Q=fort_denconst(2)/fort_denconst(1)
N=viscosity_ratio
Q=density_ratio
l=RGASRWATER
Re=ReIN
We=WeIN

print *,"N is mugas/muwater"
print *,"Q is rhogas/rhowater"
print *,"N,Q,l,Re,We ",N,Q,l,Re,We

! liquid for 0<r<1
do k=0,N1parm
   eta(k)=-dcos(PI_lsa*DBLE(k)/DBLE(N1parm))
   r1(k)=(1.0d0+eta(k))/2.0d0
enddo
! gas for 1<r<l
do k=0,N2parm
   xi(k)=-dcos(PI_lsa*DBLE(k)/DBLE(N2parm))
   r2(k)=1.0d0+(xi(k)+1.0d0)*(l-1.0d0)/2.0d0
enddo

! see pages 34 and 35 of Yaohong's thesis
do i=1,N1parm-1
   C1(i)=1.0d0
enddo
C1(0)=2.0d0
C1(N1parm)=2.0d0
do j=0,N1parm
   do k=0,N1parm
    if (j.ne.k) then
     E1(j,k)=C1(j)/(C1(k)*(eta(j)-eta(k)))
     if(mod(j+k,2).NE.0) then
      E1(j,k)=-E1(j,k)
     endif
    endif
   enddo
enddo
do j=1,N1parm-1
   E1(j,j)=-eta(j)/(2.0d0*(1.0d0-eta(j)*eta(j)))
enddo
E1(0,0)=-(2.0d0*N1parm*N1parm+1.0d0)/6.0d0
E1(N1parm,N1parm)=-E1(0,0)

! account for S_11 and S_12 terms (page 36)
do j=0,N1parm
do k=0,N1parm
 E1(j,k)=E1(j,k)*2.0d0
enddo
enddo

do j=0,N1parm
 do k=0,N1parm
  F1(j,k)=0.0d0
  do i=0,N1parm
    F1(j,k)=F1(j,k)+E1(j,i)*E1(i,k)
  enddo
 enddo
enddo

! treatment at r=0 (page 36 of Yaohong's thesis) is not used in 2D.

do i=1,N2parm-1
   C2(i)=1.0d0
enddo
C2(0)=2.0d0
C2(N2parm)=2.0d0
do j=0,N2parm
 do k=0,N2parm
  if (j.ne.k) then
   E2(j,k)=C2(j)/(C2(k)*(xi(j)-xi(k)))
   if(mod(j+k,2).NE.0) then
    E2(j,k)=-E2(j,k)
   endif
  endif
 enddo
enddo
do j=1,N2parm-1
   E2(j,j)=-xi(j)/(2.0d0*(1.0d0-xi(j)*xi(j)))
enddo
E2(0,0)=-(2.0d0*N2parm*N2parm+1.0d0)/6.0d0
E2(N2parm,N2parm)=-E2(0,0)

! account for S_21 and S_22 terms (page 36)
do j=0,N2parm
do k=0,N2parm
 E2(j,k)=E2(j,k)*2.0d0/(l-1.0d0)
enddo
enddo

do j=0,N2parm
 do k=0,N2parm
  F2(j,k)=0.0d0
  do i=0,N2parm
   F2(j,k)=F2(j,k)+E2(j,i)*E2(i,k)
  enddo
 enddo
enddo

! treatment for "r" is not needed in 2D X-Y.


alphai=aimag(alpha)
alphar=real(alpha)
alphai2=alphai*alphai
alphar2=alphar*alphar
alphari=alphai*alphar

!write(*,*) alphar, alphai, alphari, alphar2, alphai2

! initialize Ar and Ai
do j=0, M-1
   do k=0,M-1
      Ar(j,k)=0.0d0
      Ai(j,k)=0.0d0
   enddo
enddo

! these are the unknowns (correspond to matrix columns)
! u, w, p, f  (v is radial velocity; not used in 2d)
! implicit noslip condition at r=0
! implicit continuity at r=1
IU=0
IW=IU+N1parm-1
IP=IW+N1parm-1

! implicit noslip condition at r=l
IU2=IP+N1parm+1
IW2=IU2+N2parm
IP2=IW2+N2parm

IFVAR=IP2+N2parm+1

if (IFVAR+1.ne.M) then
 print *,"number of variables invalid for 2D problem"
 stop
endif

! these are the various equations
! no continuity equation at the origin, instead dp/dr-(1/Re)(u_rr+u_zz)=0
! no equations for liquid or gas at r=0 or r=1
ICONT=0
IUMOM=ICONT+N1parm  
IWMOM=IUMOM+N1parm-1  
ICONT2=IWMOM+N1parm-1  ! continuity equation at r=1 comes later
IUMOM2=ICONT2+N2parm-1  ! last comp. for IUMOM2 is dp/dr-(1/Re)u_rr=0 at r=l
IWMOM2=IUMOM2+N2parm
IPJUMP=IWMOM2+N2parm-1
IFEQN=IPJUMP+1
ITAN1parmEQN=IFEQN+1

! 2 extra equations corresponding to continuity equations for 
! liquid and gas at r=1
IEXTRAEQN=ITAN1parmEQN+1
   
if (IEXTRAEQN+2.ne.M) then
 print *,"number of equations invalid"
 stop
endif
   
NQ=1.0  ! dimensionless kinematic viscosity in water.
Q1=1.0  ! dimensionless density in water
  
NQ2=N/Q  ! dimensionless kinematic viscosity in gas
Q2=Q     ! dimensionless density in gas

! at r=0,l, dp/dr=nu (u_rr+u_zz)
! at r=0,l, u=w=0 
do k=0,N1parm
 Ar(ICONT,k+IP)=-Re*E1(0,k)  ! -Re * dp/dr
enddo
! N=mu_gas/mu_liquid
! Q=den_gas/den_liquid
do k=0,N1parm-2
 Ar(ICONT,k+IU)=NQ*Q1*F1(0,k+1) ! nu u_rr,  NQ=1 Q1=1 NQ2=N/Q Q2=Q
enddo
Ar(ICONT,IU2)=NQ*Q1*F1(0,N1parm)

do j=0,N1parm-2
 do k=0,N1parm-2
  Ar(j+ICONT+1,k+IU)=E1(j+1,k+1)  ! du/dr
 enddo
 Ar(j+ICONT+1,IU2)=E1(j+1,N1parm)

   ! dw/dz
 Ar(j+ICONT+1,j+IW)=-alphai  
 Ai(j+ICONT+1,j+IW)=alphar
enddo

do j=0,N1parm-2
 do k=0,N1parm-2
   Ar(j+IUMOM,k+IU)=Q1*NQ*F1(j+1,k+1)  ! nu u_rr
 enddo
 Ar(j+IUMOM,IU2)=Q1*NQ*F1(j+1,N1parm)  ! nu u_rr

 do k=0,N1parm
  Ar(j+IUMOM,k+IP)=-Re*E1(j+1,k)  ! -Re dp/dr
 enddo
 Ar(j+IUMOM,j+IU)=Ar(j+IUMOM,j+IU)-Q1*NQ*(alphar2-alphai2) ! nu u_zz
 Ai(j+IUMOM,j+IU)=-Q1*2.0d0*NQ*alphari  ! nu u_zz
enddo

do j=0,N1parm-2
 do k=0,N1parm-2
   Ar(j+IWMOM,k+IW)=NQ*Q1*F1(j+1,k+1)  ! nu w_rr
 enddo
 Ar(j+IWMOM,IW2)=Q1*NQ*F1(j+1,N1parm)  ! nu w_rr


 Ar(j+IWMOM,j+IW)=Ar(j+IWMOM,j+IW)-Q1*NQ*(alphar2-alphai2) ! nu w_zz
 Ai(j+IWMOM,j+IW)=-Q1*2.0d0*NQ*alphari  ! nu w_zz
 Ar(j+IWMOM,j+IP+1)=Re*alphai  ! -Re dp/dr
 Ai(j+IWMOM,j+IP+1)=-Re*alphar ! -Re dp/dr
 
enddo
  
do j=0,N2parm-2
 do k=0,N2parm-1
   Ar(j+ICONT2,k+IU2)=E2(j+1,k)  !  du/dr
 enddo
 Ar(j+ICONT2,j+IW2+1)=-alphai  ! dw/dz
 Ai(j+ICONT2,j+IW2+1)=alphar   ! dw/dz
enddo
 
do j=0,N2parm-2
 do k=0,N2parm-1
   Ar(j+IUMOM2,k+IU2)=NQ2*Q2*F2(j+1,k)  ! nu u_rr
 enddo
 do k=0,N2parm
  Ar(j+IUMOM2,k+IP2)=-Re*E2(j+1,k)  ! -Re dp/dr
 enddo
   ! nu u_zz
 Ar(j+IUMOM2,j+IU2+1)=Ar(j+IUMOM2,j+IU2+1)-Q2*NQ2*(alphar2-alphai2)
 Ai(j+IUMOM2,j+IU2+1)=-Q2*2.0d0*NQ2*alphari
enddo

! last comp. of IUMOM2 is dp/dr-(1/Re)(u_rr+u_zz)=0 at r=l
 ! nu u_rr
do k=0,N2parm-1
  Ar(N2parm+IUMOM2-1,k+IU2)=NQ2*Q2*F2(N2parm,k)  !NQ=1 Q1=1  NQ2=N/Q  Q2=Q
enddo
 ! -Re dp/dr
do k=0,N2parm
 Ar(N2parm+IUMOM2-1,k+IP2)=-Re*E2(N2parm,k)
enddo

do j=0,N2parm-2
 do k=0,N2parm-1
   Ar(j+IWMOM2,k+IW2)=NQ2*Q2*F2(j+1,k)  ! nu w_rr
 enddo
  ! nu w_zz
 Ar(j+IWMOM2,j+IW2+1)=Ar(j+IWMOM2,j+IW2+1)-Q2*NQ2*(alphar2-alphai2)
 Ai(j+IWMOM2,j+IW2+1)=-Q2*2.0d0*NQ2*alphari
  ! -Re dp/dr
 Ar(j+IWMOM2,j+IP2+1)=Re*alphai
 Ai(j+IWMOM2,j+IP2+1)=-Re*alphar
enddo

Ar(IFEQN,IU2)=-1.0

Ar(IPJUMP,IP+N1parm)=1.0
Ar(IPJUMP,IP2)=-1.0
do k=0,N1parm-2
 Ar(IPJUMP,IU+k)=-(2.0/Re)*E1(N1parm,k+1)
enddo

Ar(IPJUMP,IU2)=-(2.0/Re)*E1(N1parm,N1parm)+(2.0/Re)*N*E2(0,0)

do k=1,N2parm-1
 Ar(IPJUMP,IU2+k)=(2.0/Re)*N*E2(0,k)
enddo
Ar(IPJUMP,IFVAR)=(1.0/We)*(-alphar2+alphai2)
Ai(IPJUMP,IFVAR)=(1.0/We)*(-2.0d0*alphari)
 
Ar(ITAN1parmEQN,IU2)=N*alphai-alphai
Ai(ITAN1parmEQN,IU2)=-N*alphar+alphar
do k=0,N1parm-2
 Ar(ITAN1parmEQN,IW+k)=E1(N1parm,k+1)
enddo
Ar(ITAN1parmEQN,IW2)=E1(N1parm,N1parm)-N*E2(0,0)
do k=1,N2parm-1
 Ar(ITAN1parmEQN,IW2+k)=-N*E2(0,k)
enddo
 
! du1/dr at r=1
do k=0,N1parm-2
   Ar(IEXTRAEQN,k+IU)=E1(N1parm,k+1)
enddo
Ar(IEXTRAEQN,IU2)=E1(N1parm,N1parm)

! dw1/dz at r=1
Ai(IEXTRAEQN,IW2)=alphar
Ar(IEXTRAEQN,IW2)=-alphai

! du2/dr at r=1
do k=0,N2parm-1
   Ar(IEXTRAEQN+1,k+IU2)=E2(0,k)
enddo
! dw2/dz at r=1
Ai(IEXTRAEQN+1,IW2)=alphar
Ar(IEXTRAEQN+1,IW2)=-alphai

do i=0,M-1
   do j=0,M-1
      A(j,i)=DCMPLX(Ar(j,i),Ai(j,i))
   enddo
enddo

do j=0, M-1
   do i=0, M-1
      Br(i,j)=0.0d0
      Bi(i,j)=0.0d0
   enddo
enddo

! opposite signs since B represents terms moved to RHS
do j=0,N1parm-2
 Bi(IUMOM+j,IU+j)=-Re*Q1
enddo
do j=0,N1parm-2
 Bi(IWMOM+j,IW+j)=-Re*Q1
enddo

do j=0,N2parm-2
 Bi(IUMOM2+j,IU2+j+1)=-Re*Q2
enddo
do j=0,N2parm-2
 Bi(IWMOM2+j,IW2+j+1)=-Re*Q2
enddo
Bi(IFEQN,IFVAR)=1.0d0

do i=0,M-1
   do j=0,M-1
      B(j,i)=DCMPLX(Br(j,i),Bi(j,i))
   enddo
enddo
   
return
end subroutine initialize2d





       subroutine velinterpolation(N1parm,N2parm,r1,r2,vel_lr,   &
         vel_lz,vel_gr, vel_gz, x, y, dist, x_vel, y_vel,W1bar,W2bar)
       implicit none
       integer :: N1parm, N2parm
       real(amrex_real), dimension(0:N1parm) :: r1
       real(amrex_real), dimension(0:N2parm) :: r2
       real(amrex_real), dimension(0:N1parm) :: W1bar
       real(amrex_real), dimension(0:N2parm) :: W2bar

       complex*16, dimension(0:N1parm) :: vel_lr, vel_lz
       complex*16, dimension(0:N2parm) :: vel_gr, vel_gz
       real(amrex_real) :: x, y, dist
       real(amrex_real) :: x_vel, y_vel,W

       real(amrex_real)    :: Pi_lsa
       real(amrex_real)    :: lint,xr, xi, yr, yi
       real(amrex_real)    :: txr, txi, tyr, tyi, kr,ki, theta
       integer :: i, j


       Pi_lsa=4.0*ATAN(1.0)
       if (dist.ge.0) then
        xr=0.0d0
        xi=0.0d0
        yr=0.0d0
        yi=0.0d0
        W=0.0d0
        do i=1, N1parm
         lint=1.0d0
         txr=real(vel_lr(i))
         txi=aimag(vel_lr(i))
         tyr=real(vel_lz(i))
         tyi=aimag(vel_lz(i))
         do j=1, i-1
          lint=lint*(x-r1(j))/(r1(i)-r1(j))
         enddo
         do j=i+1,N1parm
          lint=lint*(x-r1(j))/(r1(i)-r1(j))
         enddo
         xr=xr+txr*lint
         xi=xi+txi*lint
         yr=yr+tyr*lint
         yi=yi+tyi*lint
         W=W+W1bar(i)*lint
        enddo
       else
        xr=0.0d0
        xi=0.0d0
        yr=0.0d0
        yi=0.0d0
        W=0.0d0
        do i=0, N2parm
         lint=1.0d0
         txr=real(vel_gr(i))
         txi=aimag(vel_gr(i))
         tyr=real(vel_gz(i))
         tyi=aimag(vel_gz(i))
         do j=0, i-1
          lint=lint*(x-r2(j))/(r2(i)-r2(j))
         enddo
         do j=i+1,N2parm
          lint=lint*(x-r2(j))/(r2(i)-r2(j))
         enddo
         xr=xr+txr*lint
         xi=xi+txi*lint
         yr=yr+tyr*lint
         yi=yi+tyi*lint
         W=W+W2bar(i)*lint
        enddo
       endif
       theta=2.0d0*Pi_lsa*y/yblob
       kr=cos(theta)
       ki=sin(theta)
       x_vel=radblob*(xr*ki+xi*kr)
       y_vel=radblob*(yr*ki+yi*kr)+W
       return
       end subroutine velinterpolation

       subroutine integrand_function(x,den_ratio,beta,ff)
       IMPLICIT NONE

       real(amrex_real) x,den_ratio,beta,ff,eps

       eps=one-one/den_ratio
       if (den_ratio.lt.one) then
        print *,"den_ratio invalid"
        stop
       endif
       if ((eps.lt.zero).or.(eps.gt.one)) then
        print *,"eps invalid in integrand_function"
        stop
       endif
       ff=two*(beta**3)*exp((beta-x)*(beta+x-two*eps*beta*beta/x))/(x*x)

       return
       end subroutine integrand_function

       subroutine exp_integral(lobound,hibound,beta,den_ratio,N,ff)
       IMPLICIT NONE

       integer N,i
       real(amrex_real) lobound,hibound,beta,den_ratio,h,ff,midpt,x

       h=(hibound-lobound)/N
       ff=zero
       do i=1,N
        x=lobound+(i-half)*h
        call integrand_function(x,den_ratio,beta,midpt)
        ff=ff+midpt*h
       enddo

       return
       end subroutine exp_integral

       subroutine bisection_function(beta,den_ratio,JA,ff_bisect)
       IMPLICIT NONE

       real(amrex_real) beta,den_ratio,JA,hibound,ff_bisect,ff_int
       integer N

       if (beta.le.zero) then
        print *,"beta cannot be <= 0"
        stop
       endif
  
         ! beta is dimensionless:
         ! r_interface=2 beta sqrt(alpha t)
       hibound=ten+beta
       N=1000
      
       call exp_integral(beta,hibound,beta,den_ratio,N,ff_int)
       ff_bisect=JA-ff_int

       return
       end subroutine bisection_function

       subroutine find_beta(beta,den_ratio,Jacob_number) 
       IMPLICIT NONE

       real(amrex_real) Jacob_number,den_ratio
       real(amrex_real) beta
       real(amrex_real) ff_a,ff_b,ff_c
       real(amrex_real) a,b,c
       integer niter

       beta=one
       call bisection_function(beta,den_ratio,Jacob_number,ff_a)
       niter=0
       do while (ff_a.lt.zero)
        beta=beta/two
        call bisection_function(beta,den_ratio,Jacob_number,ff_a)
        niter=niter+1
       enddo
       a=beta

       b=beta
       call bisection_function(beta,den_ratio,Jacob_number,ff_b)
       niter=0
       do while (ff_b.gt.0.0) 
        beta=2.0*beta
        call bisection_function(beta,den_ratio,Jacob_number,ff_b)
        b=beta
        niter=niter+1
       enddo
       b=beta
 
       do niter=1,64
        c=(a+b)/2.0
        beta=c
        call bisection_function(beta,den_ratio,Jacob_number,ff_c)
        if (ff_a*ff_c.le.0.0) then
         b=c
         ff_b=ff_c
        else if (ff_c*ff_b.le.0.0) then
         a=c
         ff_a=ff_c
        else
         print *,"bisection failure"
         print *,"niter=",niter
         print *,"a,b,c ",a,b,c
         print *,"fa,fb,fc ",ff_a,ff_b,ff_c
         stop
        endif
       enddo ! niter
       return

       end subroutine find_beta

       subroutine superheat_temperature( &
         alpha,beta,den_ratio,T_superheat, &
         Jacob_number,T_saturation,rstefan,time,temperature)
       IMPLICIT NONE

       real(amrex_real) alpha,den_ratio
       real(amrex_real) beta,T_superheat,Jacob_number,T_saturation
       real(amrex_real) rstefan,time,temperature,s,hibound,ff_int
       integer N

       if (time.le.zero) then
        print *,"time invalid"
        stop
       endif
       if (beta.le.zero) then
        print *,"beta cannot be <= 0"
        stop
       endif
       s=half*rstefan/sqrt(alpha*time)

       hibound=ten+s
       N=1000
       call exp_integral(s,hibound,beta,den_ratio,N,ff_int)
       temperature=ff_int*(T_saturation-T_superheat)/Jacob_number+ &
         T_superheat

       return
       end subroutine superheat_temperature

       ! before this routine:
       ! calc_error_indicator
       ! EOS_error_ind
      subroutine fort_vfracerror(  &
       tid_current, &
       error_set_count, &
       tag, &
       DIMS(tag), &
       set,clear, &
       errfab, &
       DIMS(errfab), &
       snew, &
       DIMS(snew), &
       lsnew, &
       DIMS(lsnew), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       domlo,domhi, &
       dx,xlo, &
       problo, &
       time, &
       level, &
       max_level, &
       max_level_for_use, &
       nblocks, &
       xblocks,yblocks,zblocks, &
       rxblocks,ryblocks,rzblocks, &
       ncoarseblocks, &
       xcoarseblocks,ycoarseblocks,zcoarseblocks, &
       rxcoarseblocks,rycoarseblocks,rzcoarseblocks) &
      bind(c,name='fort_vfracerror')

      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_current
      integer, INTENT(inout) :: error_set_count

      integer, INTENT(in) :: nblocks,ncoarseblocks
      real(amrex_real), INTENT(in) :: xblocks(10),yblocks(10),zblocks(10)
      real(amrex_real), INTENT(in) :: rxblocks(10),ryblocks(10),rzblocks(10)
      real(amrex_real), INTENT(in) :: xcoarseblocks(10)
      real(amrex_real), INTENT(in) :: ycoarseblocks(10)
      real(amrex_real), INTENT(in) :: zcoarseblocks(10)
      
      real(amrex_real), INTENT(in) :: rxcoarseblocks(10)
      real(amrex_real), INTENT(in) :: rycoarseblocks(10)
      real(amrex_real), INTENT(in) :: rzcoarseblocks(10)

      integer, INTENT(in) :: DIMDEC(tag)
      integer, INTENT(in) :: DIMDEC(errfab)
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(lsnew)
      integer, INTENT(in) :: set, clear
      integer, INTENT(in) :: level
      integer, INTENT(in) :: max_level
      integer, INTENT(in) :: max_level_for_use
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer   growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: domlo(SDIM), domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer, INTENT(out), target :: tag(DIMV(tag))
      integer, pointer :: tag_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: errfab(DIMV(errfab))
      real(amrex_real), pointer :: errfab_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: &
           lsnew(DIMV(lsnew),num_materials*(1+SDIM))
      real(amrex_real), pointer :: lsnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: snew(DIMV(snew),STATE_NCOMP)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real) :: x, y, z
      real(amrex_real) :: rflag
      integer i,j,k
      integer np
      integer iregions
      integer dir_local
      real(amrex_real) x_local(SDIM)
      real(amrex_real) LS_clamped,charfn
      real(amrex_real) vel_clamped(SDIM)
      real(amrex_real) temperature_clamped
      integer prescribed_flag
      integer tagflag
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer coarseblocks_available

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      if (level.lt.0) then
       print *,"level invalid vfrac error: ",level
       stop
      endif
      if (max_level.le.level) then
       print *,"max_level invalid: ",max_level
       stop
      endif
      if ((max_level_for_use.lt.0).or. &
          (max_level_for_use.gt.max_level)) then
       print *,"max_level_for_use invalid: ",max_level_for_use
       stop
      endif
      if ((nblocks.lt.0).or.(ncoarseblocks.lt.0).or. &
          (nblocks.ge.10).or.(ncoarseblocks.ge.10)) then
       print *,"nblocks or ncoarseblocks out of range"
       stop
      endif
      if ((tid_current.lt.0).or. &
          (tid_current.ge.geom_nthreads)) then
       print *,"tid_current invalid"
       stop
      endif

      tag_ptr=>tag
      errfab_ptr=>errfab
      lsnew_ptr=>lsnew
      snew_ptr=>snew

      call checkbound_int_array1(tilelo,tilehi,tag_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,errfab_ptr,1,-1)
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
       x=xsten(0,1)
       y=xsten(0,2)
       z=xsten(0,SDIM)
       do dir_local=1,SDIM
        x_local(dir_local)=xsten(0,dir_local)
       enddo

        ! GRIDS ARE CREATED AS FOLLOWS:
        ! 1. cells are tagged for refinement (tagflag==1) or not
        !    tagged (tagflag==0)
        ! 2. depending on "amr.n_error_buf", neighboring cells of tagged
        !    cells are also tagged.
        !    For each cell (i,j,k) that is tagged with tagflag==1, the
        !    following cells are also tagged:
        !    i-n_error_buf<= i* <=i+n_error_buf
        !    j-n_error_buf<= j* <=j+n_error_buf
        !    k-n_error_buf<= k* <=k+n_error_buf
        ! 3. The Berger and Rigoustos clustering algorithm is invoked which
        !    forms minimal boxes surrounding all the tagged cells, making
        !    sure that appropriate proper nesting and blocking factor
        !    conditions are satisfied.
        ! NOTE: errfab(D_DECL(i,j,k)) is initialized ultimately 
        !   from calc_error_indicator and EOS_error_ind
       rflag=errfab(D_DECL(i,j,k))
       tagflag=0

       if (level.lt.max_level_for_use) then

         ! updates "tagflag" and/or "rflag"
        call override_tagflag( &
         i,j,k, &
         level,max_level, &
         snew_ptr,lsnew_ptr, &
         xsten,nhalf,time, &
         rflag,tagflag)

        call SUB_clamped_LS(x_local,time,LS_clamped,vel_clamped, &
          temperature_clamped,prescribed_flag,dx)
        if (LS_clamped.ge.zero) then
         if (prescribed_flag.eq.0) then
          rflag=one
          tagflag=1
         else if (prescribed_flag.eq.1) then
          ! do nothing
         else
          print *,"prescribed_flag invalid"
          stop
         endif
        else if (LS_clamped.lt.zero) then
         ! do nothing
        else
         print *,"LS_clamped invalid: ",LS_clamped
         stop
        endif
        do iregions=1,number_of_source_regions

         call SUB_CHARFN_REGION(iregions,x_local,time,charfn)

         if (charfn.eq.one) then

          rflag=one
          tagflag=1

         else if (charfn.eq.zero) then
          ! do nothing
         else
          print *,"charfn invalid"
          stop
         endif
        enddo !iregions=1,number_of_source_regions

        if (rflag.eq.one) then
         tagflag=1
        else if (rflag.eq.zero) then
         ! do nothing
        else
         print *,"rflag invalid in fort_vfracerror: ",rflag
         print *,"probtype=",probtype
         print *,"axis_dir=",axis_dir
         stop
        endif

        if (ractivex.gt.zero) then
         if ((abs(x-xactive).gt.ractivex).or. &
#if (AMREX_SPACEDIM==3)
             (abs(z-zactive).gt.ractivez).or. &
#endif
             (abs(y-yactive).gt.ractivey)) then
          tagflag=0
         endif
        endif

        if (nblocks.gt.0) then
         do np=1,nblocks
          if ((abs(x-xblocks(np)).le.rxblocks(np)).and. &
#if (AMREX_SPACEDIM==3)
              (abs(z-zblocks(np)).le.rzblocks(np)).and. &
#endif
              (abs(y-yblocks(np)).le.ryblocks(np))) then
           tagflag=1
          endif
         enddo
        endif

        coarseblocks_available=1
        if ((probtype.eq.541).and.(level.le.3)) then
         coarseblocks_available=0
        else
         ! do nothing
        endif

        if (coarseblocks_available.eq.1) then

         if (ncoarseblocks.gt.0) then
          do np=1,ncoarseblocks
           if ((abs(x-xcoarseblocks(np)).ge.rxcoarseblocks(np)).or. &
#if (AMREX_SPACEDIM==3)
               (abs(z-zcoarseblocks(np)).ge.rzcoarseblocks(np)).or. &
#endif
               (abs(y-ycoarseblocks(np)).ge.rycoarseblocks(np))) then
            tagflag=0
           endif
          enddo ! do np=1,ncoarseblocks
         else if (ncoarseblocks.eq.0) then
          ! do nothing
         else
          print *,"ncoarseblocks invalid"
          stop
         endif

        else if (coarseblocks_available.eq.0) then
         ! do nothing
        else
         print *,"coarseblocks_available invalid"
         stop
        endif
        
       else if (level.ge.max_level_for_use) then
        ! do nothing
       else
        print *,"level invalid"
        stop
       endif

       if (tagflag.eq.1) then
        error_set_count=error_set_count+1
        tag(D_DECL(i,j,k))=set
       endif

      enddo
      enddo
      enddo

      return
      end subroutine fort_vfracerror

      subroutine fort_group_extrapfill ( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_group_extrapfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer icomp
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer extrecon_scomp
      integer mask_scomp
      integer burnvel_scomp
      integer tsat_scomp
      integer ncomp_per
      integer ncomp_per_burning
      integer ncomp_per_tsat
      integer ncomp_expect

      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      ncomp_per_burning=EXTRAP_PER_BURNING
      ncomp_per_tsat=EXTRAP_PER_TSAT  ! interface temperature, mass fraction

       ! extrap, velx, vely, velz
      extrecon_scomp=EXTRAPCOMP_MOF
       ! extrap, velx, vely, velz, mof recon
      mask_scomp=EXTRAPCOMP_MASK
      burnvel_scomp=EXTRAPCOMP_BURNVEL
      tsat_scomp=EXTRAPCOMP_TSAT

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 12"
       stop
      endif

       ! c++ index
      if (scomp.eq.burnvel_scomp) then
       ncomp_per=ncomp_per_burning
       ncomp_expect=num_interfaces*(1+ncomp_per)
      else if (scomp.eq.tsat_scomp) then
       ncomp_per=ncomp_per_tsat ! interface temperature, mass fraction
       ncomp_expect=num_interfaces*(1+ncomp_per)
      else if (scomp.eq.EXTRAPCOMP_DRAG) then
       ncomp_per=0
       ncomp_expect=N_DRAG
      else
       print *,"scomp invalid group extrapfill"
       print *,"scomp, ncomp, bfact, time ",scomp,ncomp,bfact,time
       print *,"num_materials,extrecon_scomp,mask_scomp,burnvel_scomp ", &
        num_materials,extrecon_scomp,mask_scomp,burnvel_scomp
       print *,"EXTRAPCOMP_DRAG=",EXTRAPCOMP_DRAG
       stop
      endif

      if (ncomp.eq.ncomp_expect) then
       ! do nothing
      else
       print *,"ncomp invalid14"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      u_ptr=>u
      do icomp=1,ncomp
       call local_filcc4D(bfact, &
        u_ptr,icomp, &
        domlo,domhi,bc)
      enddo

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

      do icomp=1,ncomp

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side,icomp)

       if (test_bc.eq.EXT_DIR) then

        print *,"exterior dirichlet BC not allowed in fort_group_extrapfill"
        print *,"test_bc=",test_bc
        stop

        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

          ! low order extrap
         call extrapBC(time,dir2,side, &
           u(D_DECL(i,j,k),icomp), &
           u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),icomp), &
           xsten,nhalf,dx,bfact)
        enddo
        enddo
        enddo
       endif            
      enddo ! icomp=1,ncomp

      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_group_extrapfill

      subroutine fort_extrapfill ( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_extrapfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u))
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:))

      integer, INTENT(in) :: bc(SDIM,2)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer box_type(SDIM)
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 13"
       stop
      endif

      if (ncomp.lt.1) then
       print *,"ncomp invalid in extrap fill"
       stop
      endif

       ! box_type(dir)=0 => CELL
       ! box_type(dir)=1 => NODE
      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (box_type(dir2).eq.0) then
        if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
         print *,"domhi+1 not divisible by bfact"
         stop
        endif
       else if (box_type(dir2).eq.1) then
        if ((domhi(dir2)/bfact)*bfact.ne.domhi(dir2)) then
         print *,"domhi not divisible by bfact"
         stop
        endif
       else
        print *,"box_type(dir2) invalid"
        stop
       endif
      enddo  ! dir2=1..sdim

      u_ptr=>u
      call efilcc(bfact, &
       u_ptr, &
       domlo,domhi,bc, &
       grid_type)

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side)

       if (test_bc.eq.EXT_DIR) then

        print *,"exterior dirichlet BC not allowed in fort_extrapfill"
        print *,"test_bc=",test_bc
        stop

        if (box_type(dir2).eq.1) then

         if (side.eq.1) then
          if (fablo(dir2).le.domlo(dir2)) then
           ext_dir_flag=1
           borderhi(dir2)=domlo(dir2)
           inside_index=domlo(dir2)
          endif
         else if (side.eq.2) then
          if (fabhi(dir2).ge.domhi(dir2)) then
           ext_dir_flag=1
           borderlo(dir2)=domhi(dir2)
           inside_index=domhi(dir2)
          endif
         else
          print *,"side invalid"
          stop 
         endif

        else if (box_type(dir2).eq.0) then

         if (side.eq.1) then
          if (fablo(dir2).lt.domlo(dir2)) then
           ext_dir_flag=1
           borderhi(dir2)=domlo(dir2)-1
           inside_index=domlo(dir2)
          endif
         else if (side.eq.2) then
          if (fabhi(dir2).gt.domhi(dir2)) then
           ext_dir_flag=1
           borderlo(dir2)=domhi(dir2)+1
           inside_index=domhi(dir2)
          endif
         else
          print *,"side invalid"
          stop
         endif

        else
         print *,"dir2 bust"
         stop
        endif

       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  
    
       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)
 
         call gridstenMAC(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf,grid_type)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

          ! sets all the physical BCs to the interior value.
         call extrapBC(time,dir2,side, &
           u(D_DECL(i,j,k)), &
           u(D_DECL(IWALL(1),IWALL(2),IWALL(3))), &
           xsten,nhalf,dx,bfact)

        enddo
        enddo
        enddo
       endif 
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_extrapfill


       subroutine fort_setfortscales(pressure_scale, &
        velocity_scale) &
       bind(c,name='fort_setfortscales')

       use probcommon_module

       IMPLICIT NONE

       real(amrex_real), INTENT(in) :: pressure_scale,velocity_scale

       if ((pressure_scale.gt.zero).and. &
           (velocity_scale.gt.zero)) then
        ! do nothing
       else
        print *,"scales invalid in set fort scales "
        print *,"pressure scale = ",pressure_scale
        print *,"velocity scale = ",velocity_scale
        stop
       endif

       if (abs(one-pressure_scale/velocity_scale**2).gt.EPS5) then
        print *,"scales invalid in set fort scales "
        print *,"pressure scale = ",pressure_scale
        print *,"velocity scale = ",velocity_scale
        stop
       endif

       global_pressure_scale=pressure_scale
       global_velocity_scale=velocity_scale

       return
       end subroutine fort_setfortscales

       subroutine fort_overridelsbc(homflag) &
       bind(c,name='fort_overridelsbc')

       use probcommon_module
       IMPLICIT NONE

       integer, INTENT(in) :: homflag


       if (homflag.eq.0) then
         ls_homflag=0
       else if (homflag.eq.1) then
         ls_homflag=1
       else
        print *,"homflag invalid in override ls bc"
        stop
       endif

       return
       end subroutine fort_overridelsbc

       subroutine fort_overridepbc(homflag,project_option) &
       bind(c,name='fort_overridepbc')

       use probcommon_module
       use global_utility_module
       IMPLICIT NONE

       integer, INTENT(in) :: homflag,project_option

        ! project_option_singular_possibleF is declared in: GLOBALUTIL.F90
       if (project_option_singular_possibleF(project_option).eq.1) then
        if (homflag.eq.0) then
         pres_homflag=0
        else if (homflag.eq.1) then
         pres_homflag=1
        else
         print *,"homflag invalid in override pbc"
         stop
        endif
       else if (project_option.eq.SOLVETYPE_VISC) then ! viscosity
        if (homflag.eq.0) then
         vel_homflag=0
        else if (homflag.eq.1) then
         vel_homflag=1
        else
         print *,"homflag invalid in override pbc 2"
         stop
        endif
       else if (project_option.eq.SOLVETYPE_HEAT) then  ! temperature
        if (homflag.eq.0) then
         temp_homflag=0
        else if (homflag.eq.1) then
         temp_homflag=1
        else
         print *,"homflag invalid in override pbc 3"
         stop
        endif
       else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
        if (homflag.eq.0) then
         species_homflag=0
        else if (homflag.eq.1) then
         species_homflag=1
        else
         print *,"homflag invalid in override pbc 4"
         stop
        endif
       else 
        print *,"project_option invalid in override pbc"
        stop
       endif

       return
       end subroutine fort_overridepbc

       subroutine fort_flush_fortran() &
       bind(c,name='fort_flush_fortran')
       IMPLICIT NONE

       call FLUSH(6)   ! unit=6 screen

       end subroutine fort_flush_fortran

       subroutine fort_set_periodic_var(periodic_in) &
       bind(c,name='fort_set_periodic_var')

       IMPLICIT NONE

       integer, INTENT(in) :: periodic_in(SDIM)
       integer dir

       do dir=1,SDIM
        fort_is_periodic(dir)=periodic_in(dir)
        if ((fort_is_periodic(dir).ne.0).and. &
            (fort_is_periodic(dir).ne.1)) then
         print *,"fort_is_periodic invalid"
         stop
        endif
       enddo

       end subroutine fort_set_periodic_var
     

      subroutine fort_initdatasolid( &
       nparts, &
       ngrow_make_distance_in, &
       im_solid_map, &
       time, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       solid,DIMS(solid), &
       LS,DIMS(LS), &
       SNEW,DIMS(SNEW), &
       dx,xlo,xhi) &
      bind(c,name='fort_initdatasolid')

      use global_distance_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: ngrow_make_distance_in
      integer, INTENT(in) :: im_solid_map(nparts)
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: bfact
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: DIMDEC(solid)
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(SNEW)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM), xhi(SDIM)

      real(amrex_real), INTENT(out), target :: solid(DIMV(solid),nparts*SDIM)
      real(amrex_real), pointer :: solid_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out), target :: &
           LS(DIMV(LS),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out), target :: SNEW(DIMV(SNEW),STATE_NCOMP)
      real(amrex_real), pointer :: SNEW_ptr(D_DECL(:,:,:),:)

      integer i,j,k,dir
      real(amrex_real) distsolid
      real(amrex_real) temp_solid_mat
      real(amrex_real) vel(SDIM)
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer im
      integer partid

      solid_ptr=>solid
      LS_ptr=>LS
      SNEW_ptr=>SNEW

      if ((nparts.lt.1).or.(nparts.ge.num_materials)) then
       print *,"nparts invalid"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"ngrow_make_distance invalid"
       stop
      endif
      if (ngrow_make_distance_in.ne.3) then
       print *,"ngrow_make_distance_in invalid"
       print *,"fort_initdatasolid"
       print *,"ngrow_make_distance_in=",ngrow_make_distance_in
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in initdata_solid"
      else if (time.lt.zero) then
       print *,"time invalid in initdata_solid"
       stop
      else
       print *,"time bust in initdata_solid"
       stop
      endif

      call checkbound_array(fablo,fabhi,solid_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      call checkbound_array(fablo,fabhi,SNEW_ptr,1,-1)

      if ((adv_dir.lt.1).or.(adv_dir.gt.2*SDIM+1)) then
       print *,"adv_dir invalid initdatasolid (10)"
       stop
      endif
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

       do partid=1,nparts

        im=im_solid_map(partid)+1
        if ((im.lt.1).or.(im.gt.num_materials)) then
         print *,"im invalid82"
         stop
        endif

        if (is_lag_part(im).eq.1) then

         if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
          call materialdistsolid(xsten(0,1),xsten(0,2),xsten(0,SDIM), &
            distsolid,time,im) 
          call velsolid(xsten(0,1),xsten(0,2),xsten(0,SDIM),vel,time,im,dx)
          call tempsolid(xsten(0,1),xsten(0,2),xsten(0,SDIM), &
           temp_solid_mat,time,im)
! LSnew,Snew (MOF data) modified in fort_renormalize_prescribe
!          LS(D_DECL(i,j,k),im)=distsolid
          do dir=1,SDIM
           solid(D_DECL(i,j,k),(partid-1)*SDIM+dir)=vel(dir)
          enddo
         else if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
                  (FSI_flag(im).eq.FSI_ICE_NODES_INIT).or. & 
                  (FSI_flag(im).eq.FSI_FLUID_NODES_INIT).or. & 
                  (FSI_flag(im).eq.FSI_SHOELE_CTML)) then 
          ! do nothing
         else
          print *,"FSI_flag invalid in fort_initdatasolid"
          print *,"im,FSI_flag(im): ",im,FSI_flag(im)
          stop
         endif

        else
         print *,"is_lag_part invalid fort_initdatasolid"
         stop
        endif

       enddo ! partid=1..nparts

      enddo
      enddo
      enddo

      return
      end subroutine fort_initdatasolid

       ! called when solidheat_flag=1,2  (not =0)
      subroutine fort_initsolidtemp( &
       nden, &
       time, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       snew,DIMS(snew), &
       lsnew,DIMS(lsnew), &
       dx,xlo) &
      bind(c,name='fort_initsolidtemp')

      use global_utility_module
      use global_distance_module

      IMPLICIT NONE

      integer, INTENT(in) :: nden
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(lsnew)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)

      real(amrex_real), INTENT(inout), target :: snew(DIMV(snew),nden)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: lsnew(DIMV(lsnew),num_materials*(SDIM+1))
      real(amrex_real), pointer :: lsnew_ptr(D_DECL(:,:,:),:)

      integer i,j,k
      integer im
      integer im_solid_crit
      integer tcomp
      real(amrex_real) distsolid
      real(amrex_real) disttest
      real(amrex_real) temp_solid_mat
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer im_solid_thermal

      snew_ptr=>snew
      lsnew_ptr=>lsnew

      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact out of range"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in initsolidtemp"
      else if (time.lt.zero) then
       print *,"time invalid in initsolidtemp"
       stop
      else
       print *,"time bust in initsolidtemp"
       stop
      endif

      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)

      im_solid_thermal=im_solid_primary()

      if ((im_solid_thermal.lt.1).or. &
          (im_solid_thermal.gt.num_materials)) then
       print *,"im_solid_thermal invalid 18"
       stop
      endif
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

       distsolid=-99999.0
       disttest=-99999.0
       im_solid_crit=0

       do im=1,num_materials
        if (is_rigid(im).eq.1) then
          ! in: fort_initsolidtemp
         call materialdistsolid(xsten(0,1),xsten(0,2), &
           xsten(0,SDIM),disttest,time,im)

         if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. &
             (FSI_flag(im).eq.FSI_SHOELE_CTML)) then
          disttest=lsnew(D_DECL(i,j,k),im)
         else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
          ! do nothing
         else
          print *,"FSI_flag invalid in fort_initsolidtemp"
          print *,"im,FSI_flag(im): ",im,FSI_flag(im)
          stop
         endif
         if (disttest.gt.distsolid) then
          distsolid=disttest
          im_solid_crit=im
         endif
        else if (is_rigid(im).eq.0) then
         ! do nothing
        else
         print *,"is_rigid invalid PROB.F90 (fort_initsolidtemp)"
         stop
        endif
       enddo ! im=1..num_materials
     
       if ((im_solid_crit.lt.1).or. &
           (im_solid_crit.gt.num_materials)) then
        print *,"im_solid_crit invalid in initsolidtemp"
        stop
       endif
 
       call tempsolid(xsten(0,1),xsten(0,2),xsten(0,SDIM), &
         temp_solid_mat,time,im_solid_crit)

       if ((FSI_flag(im_solid_crit).eq.FSI_PRESCRIBED_NODES).or. &
           (FSI_flag(im_solid_crit).eq.FSI_SHOELE_CTML)) then
        tcomp=(im_solid_crit-1)*num_state_material+ENUM_TEMPERATUREVAR+1
        temp_solid_mat=snew(D_DECL(i,j,k),tcomp)
       else if (FSI_flag(im_solid_crit).eq.FSI_PRESCRIBED_PROBF90) then 
        ! do nothing
       else
        print *,"FSI_flag(im_solid_crit) invalid in initsolidtemp"
        print *,"im_solid_crit,FSI_flag(im_solid_crit) ", &
         im_solid_crit,FSI_flag(im_solid_crit)
        stop
       endif

       do im=1,num_materials
        tcomp=(im-1)*num_state_material+ENUM_TEMPERATUREVAR+1
        if ((distsolid.ge.zero).or. &
            (im.eq.im_solid_crit)) then
         snew(D_DECL(i,j,k),tcomp)=temp_solid_mat
        endif
       enddo ! im=1..num_materials

      enddo
      enddo
      enddo

      return
      end subroutine fort_initsolidtemp


       ! grad= -dt k grad S 
       ! called from: NavierStokes::viscous_boundary_fluxes
       ! which is called from: NavierStokes::apply_pressure_grad
      subroutine fort_viscfluxfill( &
       macrolayer_size, &
       microlayer_substrate, &
       microlayer_temperature_substrate, &
       freezing_model, &
       saturation_temp, &
       nsolve, &
       dir, &
       xlo,dx, &
       velbc, &
       tempbc, &
       temp_dombc, &
       LS,DIMS(LS), &
       area,DIMS(area), &
       xflux,DIMS(xflux), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       domlo,domhi, &
       dt, &
       solidheat_flag, &
       project_option, &
       time) &
      bind(c,name='fort_viscfluxfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: nsolve
      integer, INTENT(in) :: dir
      integer dir2
      integer, INTENT(in) :: solidheat_flag
      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(area)
      integer, INTENT(in) :: DIMDEC(xflux)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer growloMAC(3),growhiMAC(3)
      integer growlo_strip(3),growhi_strip(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in), target :: LS(DIMV(LS),num_materials*(SDIM+1))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: area(DIMV(area))
      real(amrex_real), pointer :: area_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: xflux(DIMV(xflux),nsolve)
      real(amrex_real), pointer :: xflux_ptr(D_DECL(:,:,:),:)
      integer, INTENT(in) :: velbc(SDIM,2,SDIM)
      integer, INTENT(in) :: tempbc(SDIM,2)
      integer, INTENT(in) :: temp_dombc(SDIM,2)
      real(amrex_real), INTENT(in) :: macrolayer_size(num_materials)
      integer, INTENT(in) :: microlayer_substrate(num_materials)
      real(amrex_real), INTENT(in) :: microlayer_temperature_substrate(num_materials)
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)

      integer i,j,k,ii,jj,kk
      integer side
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) dist
      real(amrex_real) LSleft,LSright
      real(amrex_real) nn
      real(amrex_real) tempflux
      real(amrex_real) xflux_local
      integer im,im1,im2,ireverse
      integer im_solid_tempflux
      real(amrex_real) LL,TSAT,TSUPER,thermal_layer
      integer local_freezing_model,heat_flux_model,iten

      im_solid_tempflux=im_solid_primary()
 
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
 
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif

      if ((solidheat_flag.lt.0).or. &
          (solidheat_flag.gt.2)) then
       print *,"solidheat_flag invalid"
       stop
      endif

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid viscfluxfill"
       stop
      endif

      if (dt.lt.zero) then
       print *,"dt invalid"
       stop
      endif

      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      area_ptr=>area
      call checkbound_array1(fablo,fabhi,area_ptr,0,dir)
      xflux_ptr=>xflux
      call checkbound_array(fablo,fabhi,xflux_ptr,0,dir)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
          growloMAC,growhiMAC,0,dir) 

      if (project_option.eq.SOLVETYPE_VISC) then ! viscosity
       if (nsolve.ne.AMREX_SPACEDIM) then
        print *,"nsolve invalid"
        stop
       endif
      else if (project_option.eq.SOLVETYPE_HEAT) then !thermal conduction
       if (nsolve.ne.1) then
        print *,"nsolve invalid"
        stop
       endif
      else
       print *,"project_option not supported; fort_viscfluxfill"
       stop
      endif
      
      if (project_option.eq.SOLVETYPE_VISC) then ! viscosity
       ! do nothing
      else if (project_option.eq.SOLVETYPE_HEAT) then ! thermal conduction

       do im=1,num_materials

        if (is_rigid(im).eq.0) then
         ! do nothing
        else if (is_rigid(im).eq.1) then

         if ((im_solid_tempflux.lt.1).or. &
             (im_solid_tempflux.gt.num_materials)) then
          print *,"im_solid_tempflux invalid"
          stop
         endif
      
         ! 0=diffuse in solid
         ! 1=dirichlet
         ! 2=neumann (insulating or heat flux source/sink)
         if (solidheat_flag.eq.2) then

          heat_flux_model=0

          do ireverse=0,1
          do im1=1,num_materials-1
          do im2=im1+1,num_materials
           call get_iten(im1,im2,iten)
           LL=get_user_latent_heat(iten+ireverse*num_interfaces,room_temperature,1)
           local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
           TSAT=saturation_temp(iten+ireverse*num_interfaces)

           if ((LL.ne.zero).and. &
               (local_freezing_model.eq.0).and. &
               ((microlayer_substrate(im1).eq.im).or. &
                (microlayer_substrate(im2).eq.im))) then
            heat_flux_model=1
            if (microlayer_substrate(im1).eq.im) then
             TSUPER=microlayer_temperature_substrate(im1)
             thermal_layer=macrolayer_size(im1) 
             if ((TSUPER.le.zero).or.(thermal_layer.le.zero)) then
              print *,"TSUPER or thermal_layer invalid"
              stop
             endif
              ! -k dT/dx * n
              ! case 1: solid left, fluid right  n=-1
              !   dT/dx=TSAT-TSUPER  
              ! case 2: solid right, fluid left, n=1
              !   dT/dx=TSUPER-TSAT 
             tempflux=-fort_heatviscconst(im1)*(TSUPER-TSAT)/thermal_layer
            else if (microlayer_substrate(im2).eq.im) then
             TSUPER=microlayer_temperature_substrate(im2)
             thermal_layer=macrolayer_size(im2)
             if ((TSUPER.le.zero).or.(thermal_layer.le.zero)) then
              print *,"TSUPER or thermal_layer invalid"
              stop
             endif
              ! -k dT/dx * n
              ! case 1: solid left, fluid right  n=-1
              !   dT/dx=TSAT-TSUPER  
              ! case 2: solid right, fluid left, n=1
              !   dT/dx=TSUPER-TSAT 
             tempflux=-fort_heatviscconst(im2)*(TSUPER-TSAT)/thermal_layer
            else
             print *,"microlayer_substrate invalid"
             stop
            endif
           else if ((LL.eq.zero).or. &
                    (local_freezing_model.gt.0).or. &
                    ((microlayer_substrate(im1).ne.im).and. &
                     (microlayer_substrate(im2).ne.im))) then
            ! do nothing
           else
            print *,"LL, local_freezing_model, or microlayer_substrate invalid"
            stop
           endif
          enddo ! im2
          enddo ! im1
          enddo ! ireverse

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
           print *,"dir invalid viscfluxfill 2"
           stop
          endif

          do k=growloMAC(3),growhiMAC(3)
          do j=growloMAC(2),growhiMAC(2)
          do i=growloMAC(1),growhiMAC(1)
           ! dir=0..sdim-1
           call gridstenMAC(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf,dir)

           LSleft=LS(D_DECL(i-ii,j-jj,k-kk),im)
           LSright=LS(D_DECL(i,j,k),im)
            ! n points from from fluid to solid
            ! Neumann BC approximation assumes that the
            ! interface has a staircase (rasterized) reconstruction.
           if ((LSright.ge.zero).and.(LSleft.lt.zero)) then
            nn=one
           else if ((LSleft.ge.zero).and.(LSright.lt.zero)) then
            nn=-one
           else
            nn=zero
           endif

             ! flux is -k dt grad T
             ! tempfluxsolid returns Q=-k dT/dx * n  (Q>0 => cooling)
             ! n points from fluid to solid
             ! solid on left and fluid on right => n=-1
             ! code should return Q dt n
             ! units of Q are (erg/(cm s Kelvin)) kelvin/cm=
             ! erg/(cm^2 s)
           if (LSleft*LSright.le.zero) then

            if (heat_flux_model.eq.0) then

             call tempfluxsolid(xsten(0,1), &
              xsten(0,2), &
              xsten(0,SDIM), &
              tempflux,time,dir)

             xflux_local=dt*tempflux*nn

            else if (heat_flux_model.eq.1) then
             ! 
             ! tempflux=-k (TSUPER-TSAT)/thermal_layer=-k dT/dx n
             xflux_local=dt*tempflux*nn
            else
             print *,"heat_flux_model invalid"
             stop
            endif

            xflux(D_DECL(i,j,k),1)=xflux_local
 
           else if (LSleft*LSright.ge.zero) then
            ! do nothing
           else
            print *,"LSleft or LSright bust"
            stop
           endif

          enddo
          enddo
          enddo

         ! 0=diffuse in solid
         ! 1=dirichlet
         ! 2=neumann 
         else if ((solidheat_flag.eq.0).or. &
                  (solidheat_flag.eq.1)) then
          ! do nothing
         else
          print *,"solidheat_flag invalid"
          stop
         endif

        else
         print *,"is_rigid invalid PROB.F90"
         stop
        endif
    
       enddo ! im=1..num_materials
 
       do dir2=1,3
        growlo_strip(dir2)=growloMAC(dir2) 
        growhi_strip(dir2)=growhiMAC(dir2) 
       enddo

       side=1
       if ((temp_dombc(dir+1,side).eq.FOEXTRAP).and. &
           (tilelo(dir+1).eq.domlo(dir+1))) then

        growhi_strip(dir+1)=growlo_strip(dir+1)

        do k=growlo_strip(3),growhi_strip(3) 
        do j=growlo_strip(2),growhi_strip(2) 
        do i=growlo_strip(1),growhi_strip(1) 

         ! do nothing - no problems yet with a prescribed heat flux at
         ! the bottom

        enddo
        enddo
        enddo
       endif  ! side==1 case (left side)

       do dir2=1,3
        growlo_strip(dir2)=growloMAC(dir2) 
        growhi_strip(dir2)=growhiMAC(dir2) 
       enddo

       side=2
       if ((temp_dombc(dir+1,side).eq.FOEXTRAP).and. &
           (tilehi(dir+1).eq.domhi(dir+1))) then

        growlo_strip(dir+1)=growhi_strip(dir+1)

        do k=growlo_strip(3),growhi_strip(3) 
        do j=growlo_strip(2),growhi_strip(2) 
        do i=growlo_strip(1),growhi_strip(1) 

         ! dir=0..sdim-1
         call gridstenMAC(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf,dir)

         ! cooling disk (top wall side==2)
         ! in order to cool from the top, q=-k grad T>0
         ! q=-k (TCOLD-T_ROOM)>0
         if ((probtype.eq.601).and. &
             (side.eq.2).and. &
             (dir.eq.SDIM-1)) then

          if (SDIM.eq.3) then
           dist=sqrt((xsten(0,1)-xblob)**2+(xsten(0,2)-yblob)**2)-radblob
          else if (SDIM.eq.2) then
           dist=abs(xsten(0,1))-radblob
          else
           print *,"dimension bust"
           stop
          endif

          ! Q=- k grad T 
          ! k grad T<0 => Q>0
          if (dist.le.zero) then
           if (radblob2.le.zero) then
            print *,"cooling Q should be positive"
            stop
           endif
           xflux_local=dt*radblob2
          else
           xflux_local=zero
          endif

          xflux(D_DECL(i,j,k),1)=xflux_local

         endif ! cooling disk

        enddo
        enddo
        enddo
       endif  ! side==2 case (right side)

      else
       print *,"project_option not supported; fort_viscfluxfill"
       stop
      endif  

      return
      end subroutine fort_viscfluxfill


       subroutine init_initdata(nc, &
         freezing_model, &
         distribute_from_target, &
         saturation_temp, &
         dx)
       use global_utility_module
       use supercooled_exact_sol
       IMPLICIT NONE

       integer, INTENT(in) :: nc
       integer nc_expect
       real(amrex_real), INTENT(in) :: dx(SDIM)
       integer, INTENT(in) :: freezing_model(2*num_interfaces)
       integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
       real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
       real(amrex_real) :: lmSt

       integer ireverse,im,im_opp,iten,local_freezing_model
       integer im_source,im_dest
       real(amrex_real) LL,TSAT
       real(amrex_real) cp_source,k_source,TDIFF,rho_source,rho_dest
       real(amrex_real) TDIFF_source
       real(amrex_real) TDIFF_dest
       real(amrex_real) TDIFF_solid
       real(amrex_real) den_ratio
       integer im_solid_initdata

       im_solid_initdata=im_solid_primary()

       nc_expect=STATE_NCOMP
       if (nc.ne.nc_expect) then
        print *,"fort: nc invalid"
        stop
       endif

       if (probtype.eq.802) then ! dissolution
        xlodiss_initdata=zero
        xhidiss_initdata=radblob
        ndiss_initdata=NINT(radblob/dx(SDIM))-1
        if (ndiss_initdata.ge.1000) then
         print *,"ndiss_initdata too big"
         stop
        endif
        xhidiss_initdata=ndiss_initdata*dx(SDIM)
        dtdiss_initdata=one

        do ispace_initdata=1,ndiss_initdata-1 
         lowerdiag_initdata(ispace_initdata)=denfact/(dx(SDIM)**2)
         upperdiag_initdata(ispace_initdata)=denfact/(dx(SDIM)**2)
         diag_initdata(ispace_initdata)= &
                 -lowerdiag_initdata(ispace_initdata)- &
                 upperdiag_initdata(ispace_initdata)-1.0/dtdiss_initdata
         rhs_initdata(ispace_initdata)=zero
         if (ispace_initdata.eq.1) then
          rhs_initdata(ispace_initdata)=rhs_initdata(ispace_initdata)- &
                  lowerdiag_initdata(ispace_initdata)
         endif
        enddo
         ! in: GLOBALUTIL.F90
        call tridiag_solve(lowerdiag_initdata,upperdiag_initdata, &
                diag_initdata,ndiss_initdata-1,rhs_initdata,soln_initdata)
       endif  ! probtype=802

       do ireverse=0,1
       do im=1,num_materials-1
       do im_opp=im+1,num_materials
        if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
         print *,"im or im_opp bust 8"
         stop
        endif
         ! 1<=iten<=num_interfaces
        call get_iten(im,im_opp,iten)
        LL=get_user_latent_heat(iten+ireverse*num_interfaces,room_temperature,1)
        local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
        TSAT=saturation_temp(iten+ireverse*num_interfaces)
        if (is_hydrate_freezing_modelF(local_freezing_model).eq.1) then 
         if (num_species_var.eq.1) then
          ! do nothing
         else
          print *,"must define species var if hydrate model"
          stop
         endif
        else if (is_hydrate_freezing_modelF(local_freezing_model).eq.0) then 
         ! do nothing
        else
         print *,"is_hydrate_freezing_modelF invalid"
         stop
        endif

        fort_alpha(iten+ireverse*num_interfaces)=zero
        fort_beta(iten+ireverse*num_interfaces)=zero
        fort_expansion_factor(iten+ireverse*num_interfaces)=zero
        fort_stefan_number(iten+ireverse*num_interfaces)=zero
        fort_jacob_number(iten+ireverse*num_interfaces)=zero
        fort_time_radblob(iten+ireverse*num_interfaces)=zero

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
          ! we find the rate for the "sucking" problem which
          ! has a thin thermal layer (interface moves towards the source).
         cp_source=get_user_stiffCP(im_source) ! J/Kelvin
         k_source=get_user_heatviscconst(im_source) ! W/(m Kelvin)

          ! in: init_initdata
         if ((im_source.ge.1).and.(im_source.le.num_materials).and. &
             (im_dest.ge.1).and.(im_dest.le.num_materials)) then
          TDIFF_source=fort_tempconst(im_source)-TSAT
          TDIFF_dest=fort_tempconst(im_dest)-TSAT
          TDIFF_solid=TDIFF_dest
          if ((im_solid_initdata.ge.1).and. &
              (im_solid_initdata.le.num_materials)) then
           TDIFF_solid=fort_tempconst(im_solid_initdata)-TSAT
          else if (im_solid_initdata.eq.0) then
           ! do nothing
          else
           print *,"im_solid_initdata invalid"
           stop
          endif
         else
          print *,"im_source or im_dest invalid"
          stop
         endif

         if (LL.lt.zero) then
          TDIFF_source=-TDIFF_source
          TDIFF_dest=-TDIFF_dest
          TDIFF_solid=-TDIFF_solid
         else if (LL.gt.zero) then
          ! do nothing
         else
          print *,"LL invalid"
          stop
         endif 
         TDIFF=max(TDIFF_source,TDIFF_dest)
         TDIFF=max(TDIFF,TDIFF_solid)
         TDIFF=abs(TDIFF)

         rho_source=fort_denconst(im_source) ! kg/m^3 
         rho_dest=fort_denconst(im_dest) ! kg/m^3 
         den_ratio=max(rho_dest,rho_source)/min(rho_dest,rho_source)
         if ((den_ratio.gt.1.0D+5).or.(den_ratio.lt.one)) then
          print *,"den_ratio invalid"
          stop
         endif

         if (distribute_from_target(iten+ireverse*num_interfaces).eq.0) then
          fort_expansion_factor(iten+ireverse*num_interfaces)=one-rho_dest/rho_source
         else if (distribute_from_target(iten+ireverse*num_interfaces).eq.1) then
          fort_expansion_factor(iten+ireverse*num_interfaces)=one-rho_source/rho_dest
         else
          print *,"distribute_from_target invalid"
          stop
         endif


         ! (W/(m Kelvin))/((kg/m^3)(J/(kg Kelvin)))
         ! (J/(s m Kelvin))/((1/m^3)(J/Kelvin))
         ! (1/(s m))/(1/m^3)
         ! (1/(s))/(1/m^2)=m^2/s
         fort_alpha(iten+ireverse*num_interfaces)=k_source/(rho_source*cp_source) 

         if (1.eq.1) then
          print *,"TDIFF=",TDIFF
         endif

         ! (J/(kg Kelvin)) Kelvin/(J/kg)=1
         fort_stefan_number(iten+ireverse*num_interfaces)= &
                 cp_source*TDIFF/abs(LL)
        
         fort_jacob_number(iten+ireverse*num_interfaces)= &
            (rho_source/rho_dest)* &
            fort_stefan_number(iten+ireverse*num_interfaces)

         ! solidification
         ! circular freeze.
         if (den_ratio.lt.10.0) then
          call find_lambda(lmSt, &
            fort_stefan_number(iten+ireverse*num_interfaces))

         ! spherical boiling
         else if (den_ratio.ge.10.0) then

          call find_beta(lmSt, &
            den_ratio,fort_jacob_number(iten+ireverse*num_interfaces))

         else
          print *,"den_ratio bust"
          stop
         endif

         fort_beta(iten+ireverse*num_interfaces)=lmSt

          ! sqrt(alpha time_radblob)*two*lmSt=radblob 
         call solidification_front_time(lmSt, &
          fort_alpha(iten+ireverse*num_interfaces), &
          fort_time_radblob(iten+ireverse*num_interfaces), &
          radblob)

         print *,"iten,ireverse ",iten,ireverse
         print *,"im_source,im_dest ",im_source,im_dest 
         print *,"fort_expansion_factor= ", &
          fort_expansion_factor(iten+ireverse*num_interfaces)
         print *,"distribute_from_target= ", &
          distribute_from_target(iten+ireverse*num_interfaces)
         print *,"den_ratio= ",den_ratio
         print *,"Stefan_number= ",fort_stefan_number(iten+ireverse*num_interfaces)
         print *,"Jacob_number= ",fort_jacob_number(iten+ireverse*num_interfaces)
         print *,"lmSt= ",lmSt
         print *,"alpha= ",fort_alpha(iten+ireverse*num_interfaces)
         print *,"beta= ",fort_beta(iten+ireverse*num_interfaces)
         print *,"time_radblob is the time for the front to grow from"
         print *,"r=0 to r=radblob"
         print *,"time_radblob=",fort_time_radblob(iten+ireverse*num_interfaces)
         print *,"radius doubling time=4 * time_radblob - time_radblob=", &
          three*fort_time_radblob(iten+ireverse*num_interfaces)
         print *,"front location: 2 beta sqrt(alpha t) "
        else if (LL.eq.zero) then
         ! do nothing
        else
         print *,"LL invalid"
         stop
        endif
       enddo ! im_opp
       enddo ! im
       enddo ! ireverse

       return
       end subroutine init_initdata

        ! fort_initgridmap is called from:
        !   NavierStokes::post_restart() (called from AmrCore::restart)
        !   NavierStokes::initData() 
       subroutine fort_initgridmap( &
        verbose, &
        ioproc, &
        max_level, &
        bfact_space_level, & 
        bfact_grid_level, & 
        domlo,domhi, &
        dx, &
        problo,probhi) &
       bind(c,name='fort_initgridmap')

       use global_utility_module
       use supercooled_exact_sol

       IMPLICIT NONE

       integer, INTENT(in) :: verbose
       integer, INTENT(in) :: ioproc
       integer, INTENT(in) :: max_level
       integer, INTENT(in) :: bfact_space_level(0:max_level)
       integer, INTENT(in) :: bfact_grid_level(0:max_level)
       integer, INTENT(in) :: domlo(SDIM)
       integer, INTENT(in) :: domhi(SDIM)
       real(amrex_real), INTENT(in) :: dx(SDIM)
       real(amrex_real), INTENT(in) :: problo(SDIM)
       real(amrex_real), INTENT(in) :: probhi(SDIM)
       integer, parameter :: nhalf=1
       real(amrex_real) xsten(-nhalf:nhalf)
       integer bfactmax
       integer ilev,max_ncell,dir,inode,i
       integer ncell(SDIM)
       real(amrex_real) dxlevel(SDIM)
       integer domlo_level(SDIM)
       integer domhi_level(SDIM)

       if (max_level.lt.0) then
        print *,"max_level invalid"
        stop
       endif

       bfactmax=0
       do ilev=0,max_level
        if ((bfact_space_level(ilev).lt.1).or. &
            (bfact_grid_level(ilev).lt.2)) then
         print *,"bfact invalid200 in initgridmap"
         stop
        endif
        if (bfact_space_level(ilev).gt.bfactmax) then
         bfactmax=bfact_space_level(ilev)
        endif
        if (bfact_grid_level(ilev).gt.bfactmax) then
         bfactmax=bfact_grid_level(ilev)
        endif
       enddo ! ilev=0,max_level

       max_ncell=0
       do dir=1,SDIM
        if (domlo(dir).ne.0) then
         print *,"domlo invalid"
         stop
        endif
        ncell(dir)=domhi(dir)-domlo(dir)+1
        if (ncell(dir).lt.2) then
         print *,"ncell invalid"
         stop
        endif
        if (ncell(dir).gt.max_ncell) then
         max_ncell=ncell(dir)
        endif
       enddo ! dir

       if (max_ncell.lt.2) then
        print *,"max_ncell invalid"
        stop
       endif

       do ilev=1,max_level
        max_ncell=max_ncell*2
       enddo

       if (bfactmax.lt.8) then
        bfactmax=8
       endif

       if (grid_cache_allocated.eq.0) then

        cache_index_low=-4*bfactmax
        cache_index_high=2*max_ncell+4*bfactmax
        cache_max_level=max_level

        print *,"allocate grid_cache"
        do dir=1,SDIM
         print *,"dir,domlo,domhi ",dir,domlo(dir),domhi(dir)
         print *,"dir,problo,probhi ",dir,problo(dir),probhi(dir)
         print *,"dir,dx ",dir,dx(dir)
        enddo
        print *,"cache_index_low ",cache_index_low
        print *,"cache_index_high ",cache_index_high
        print *,"bfactmax,max_ncell ",bfactmax,max_ncell
        print *,"cache_max_level ",cache_max_level

        allocate(grid_cache(0:cache_max_level, &
         cache_index_low:cache_index_high,SDIM))

        do dir=1,SDIM
         dxlevel(dir)=dx(dir)
         domlo_level(dir)=domlo(dir)
         domhi_level(dir)=domhi(dir)
        enddo

        do ilev=0,cache_max_level

         if (verbose.gt.0) then
          if (ioproc.eq.1) then
           print *,"fort_initgridmap: ilev,bfact,dxlevel(1) ", &
              ilev,bfact_space_level(ilev),dxlevel(1) 
          else if (ioproc.eq.0) then
           ! do nothing
          else
           print *,"ioproc invalid"
           stop
          endif
         else if (verbose.le.0) then
          ! do nothing
         else
          print *,"verbose invalid"
          stop
         endif

         do dir=1,SDIM
          if (domlo_level(dir).ne.0) then
           print *,"domlo_level invalid"
           stop
          endif
          do i=domlo_level(dir)-2*bfactmax,domhi_level(dir)+2*bfactmax
           inode=2*i
           if ((inode.lt.cache_index_low).or. &
               (inode+1.gt.cache_index_high)) then
            print *,"icell outside of cache range"
            stop
           endif
           call gridsten1D(xsten,problo,i,domlo, &
            bfact_space_level(ilev),dxlevel,dir,nhalf) 
           grid_cache(ilev,inode,dir)=xsten(0)
           grid_cache(ilev,inode+1,dir)=xsten(1)
           if (1.eq.0) then
            print *,"lev,inode,dir,x,xnxt ",ilev,inode,dir,xsten(0),xsten(1)
           endif
          enddo !i=domlo_level(dir)-2*bfactmax,domhi_level(dir)+2*bfactmax
         enddo ! dir=1,SDIM
         do dir=1,SDIM
          dxlevel(dir)=half*dxlevel(dir)
          domlo_level(dir)=2*domlo_level(dir)
          domhi_level(dir)=2*(domhi_level(dir)+1)-1
         enddo
        enddo !ilev=0,cache_max_level

       else if (grid_cache_allocated.eq.1) then
        ! do nothing
       else
        print *,"grid_cache_allocated invalid"
        stop
       endif

       grid_cache_allocated=1

       return
       end subroutine fort_initgridmap

       subroutine fort_initdata_alloc( &
        nc, &
        freezing_model, &
        distribute_from_target, &
        saturation_temp, &
        dx) &
       bind(c,name='fort_initdata_alloc')

       use supercooled_exact_sol
       use global_utility_module
       IMPLICIT NONE

       integer, INTENT(in) :: nc
       real(amrex_real), INTENT(in) :: dx(SDIM)
       integer, INTENT(in) :: freezing_model(2*num_interfaces)
       integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
       real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)

       call init_initdata(nc, &
        freezing_model, &
        distribute_from_target, &
        saturation_temp, &
        dx)

       return
       end subroutine fort_initdata_alloc

       subroutine fort_initdata( &
        tid, &
        adapt_quad_depth, &
        level, &
        max_level, &
        time, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        nc, &
        saturation_temp, &
        scal, &
        DIMS(scal), &
        refineden, &
        DIMS(refineden), &
        LS, &
        DIMS(LS), &
        dx,xlo,xhi, &
        centroid_noise_factor) &
       bind(c,name='fort_initdata')

       use MOF_routines_module
       use geometry_intersect_module
       use hydrateReactor_module
       use supercooled_exact_sol
       use unimaterialChannel_module
       use global_utility_module
       use global_distance_module
       use marangoni
       use USERDEF_module
       use HELIX_module
       use TSPRAY_module
       use CAV2Dstep_module
       use ZEYU_droplet_impact_module
       use rigid_FSI_module
       use sinking_particle_module
       use stackvolume_module

       IMPLICIT NONE

       integer, INTENT(in) :: adapt_quad_depth,tid
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer growlo(3),growhi(3)
       integer, INTENT(in) :: bfact
       integer, INTENT(in) :: level
       integer, INTENT(in) :: max_level
       integer, INTENT(in) :: nc
       integer imls
       real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
       real(amrex_real), INTENT(in) :: time
       integer, INTENT(in) :: DIMDEC(scal)
       integer, INTENT(in) :: DIMDEC(refineden)
       integer, INTENT(in) :: DIMDEC(LS)

       real(amrex_real), INTENT(inout), target :: scal(DIMV(scal),nc)
       real(amrex_real), pointer :: scal_ptr(D_DECL(:,:,:),:)

       real(amrex_real), INTENT(inout), target :: &
         refineden(DIMV(refineden),NUM_CELL_REFINE_DENSITY)
       real(amrex_real), pointer :: refineden_ptr(D_DECL(:,:,:),:)

       real(amrex_real), INTENT(inout), target :: &
               LS(DIMV(LS),num_materials*(1+SDIM))
       real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)

       real(amrex_real), INTENT(in) :: dx(SDIM)
       real(amrex_real), INTENT(in) :: xlo(SDIM), xhi(SDIM)
       real(amrex_real), INTENT(in) :: centroid_noise_factor(num_materials)

       real(amrex_real) centroid_noise,noise_amplitude
       integer ibase
       integer refine_comp
       integer ic,jc,kc
       integer ifine,jfine,kfine
       integer nfine
       integer n,im
       integer dir

       integer, parameter :: num_particles=0
       real(amrex_real) :: particle_list(1,SDIM+1)

       real(amrex_real) vfracsum_test

       real(amrex_real) fluiddata(num_materials,2*SDIM+2)
       real(amrex_real) mofdata(num_materials*ngeom_recon)
       real(amrex_real) distbatch(num_materials)
       real(amrex_real) LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)
       real(amrex_real) err
       real(amrex_real) vofdark(num_materials)
       real(amrex_real) voflight(num_materials)
       real(amrex_real) cendark(num_materials,SDIM)
       real(amrex_real) cenlight(num_materials,SDIM)

       integer, parameter :: nhalf=3
       integer, parameter :: nhalf2=1
       integer :: nmax
       real(amrex_real) xsten(-nhalf:nhalf,SDIM)
       real(amrex_real) xsten2(-nhalf2:nhalf2,SDIM)

       integer i1,j1,k1,k1lo,k1hi
       real(amrex_real) xpos(SDIM)
       real(amrex_real) scalc(nc)
       real(amrex_real) LSc(num_materials*(1+SDIM))
       real(amrex_real) x,y,z,rr
       real(amrex_real) volcell
       real(amrex_real) cencell(SDIM)
       real(amrex_real) centroid_absolute(SDIM)
       
       real(amrex_real) debug_vfrac_sum
       real(amrex_real) vel(SDIM)
       real(amrex_real) temp,dens,ccnt
       real(amrex_real) distsolid

       real(amrex_real) den_jwl_left,den_jwl_right
       real(amrex_real) temp_jwl_left,temp_jwl_right
       real(amrex_real) e_jwl_left,e_jwl_right
       real(amrex_real) p_jwl_left,p_jwl_right
       real(amrex_real) u_jwl_left,u_jwl_right
       real(amrex_real) xshock
       real(amrex_real) den_jwl,denroom,e_jwl,e_room,eps_benard
       real(amrex_real) gamma_jwl
       real(amrex_real) p_hyd,p_jwl,p_room,preshydro,rhohydro
       real(amrex_real) temp_jwl,temp_slope,temproom
       real(amrex_real) water_temp
       integer imattype,isten
       integer max_levelstack
       integer vofcomp_raw
       integer vofcomp_recon
       real(amrex_real) jumpval
       real(amrex_real) voflist(num_materials)

       integer im_source,im_dest,ireverse,iten
       integer im_refine_density
       real(amrex_real) L_ice_melt,TSAT,T_EXTREME,cp_melt,k_melt,rstefan
       real(amrex_real) T_FIELD
       real(amrex_real) den_ratio
       real(amrex_real) dxmaxLS
       integer im_solid_initdata
       real(amrex_real) lsnormal(num_materials,SDIM)
       integer lsnormal_valid(num_materials)
       real(amrex_real) ls_intercept(num_materials)
       integer doubly_flag
       real(amrex_real) local_state(num_materials*num_state_material)
       real(amrex_real) massfrac_parm(num_species_var+1)
       integer local_ibase
       integer tessellate
       integer bcflag
       integer, PARAMETER :: from_boundary_hydrostatic=0
       integer, parameter :: continuous_mof=STANDARD_MOF
       integer cmofsten(D_DECL(-1:1,-1:1,-1:1))
       real(amrex_real) theta_initdata
       real(amrex_real) concentration_initdata
       real(amrex_real) concen1_initdata,concen2_initdata

       real(amrex_real) local_neg_force(SDIM)
       real(amrex_real) local_vel(SDIM)
       real(amrex_real) local_vort
       real(amrex_real) local_energy_moment

       real(amrex_real) zcrit
       real(amrex_real) z_extrema
       real(amrex_real) a1,a2,D2
       real(amrex_real) pz,pz_sanity,fpz,gpz
       real(amrex_real) T_HOT,T_COLD
       real(amrex_real), parameter :: slope_checker=0.25d0

       if ((level.ge.0).and. &
           (level.le.max_level)) then
        ! do nothing
       else
        print *,"level invalid"
        stop
       endif

       scal_ptr=>scal
       refineden_ptr=>refineden
       LS_ptr=>LS

       tessellate=0

       bcflag=0

       if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
        print *,"tid invalid"
        stop
       endif

       if ((time.ge.zero).and.(time.le.1.0D+20)) then
        ! do nothing
       else if (time.ge.1.0D+20) then
        print *,"WARNING time.ge.1.0D+20 in initdata"
       else if (time.lt.zero) then
        print *,"time invalid in initdata"
        stop
       else
        print *,"time bust in initdata"
        stop
       endif

       call get_dxmaxLS(dx,bfact,dxmaxLS)

       im_solid_initdata=im_solid_primary()

       nmax=POLYGON_LIST_MAX ! in: fort_initdata

       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if ((adapt_quad_depth.lt.1).or.(adapt_quad_depth.gt.10)) then
        print *,"adapt_quad_depth invalid"
        stop
       endif
       max_levelstack=adapt_quad_depth

      if (NUM_CELL_REFINE_DENSITY.eq. &
          num_materials_compressible*ENUM_NUM_REFINE_DENSITY_TYPE) then
       ! do nothing
      else
       print *,"NUM_CELL_REFINE_DENSITY invalid"
       stop
      endif

       if (nc.ne.STATE_NCOMP) then
        print *,"nc invalid: ",nc
        print *,"STATE_NCOMP: ",STATE_NCOMP
        stop
       endif
      
       if (num_state_base.ne.2) then
        print *,"num_state_base invalid: ",num_state_base
        stop
       endif

       call checkbound_array(fablo,fabhi,scal_ptr,1,-1)
       call checkbound_array(fablo,fabhi,refineden_ptr,1,-1)
       call checkbound_array(fablo,fabhi,LS_ptr,1,-1)

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

       do kc=growlo(3),growhi(3)
       do jc=growlo(2),growhi(2)
       do ic=growlo(1),growhi(1)
        
        do n=1,nc
         scalc(n)=zero
        enddo
        do imls=1,num_materials*(1+SDIM)
         LSc(imls)=zero
        enddo
        volcell=zero
        do dir=1,SDIM
         cencell(dir)=zero
        enddo

        call gridsten(xsten,xlo,ic,jc,kc,fablo,bfact,dx,nhalf)

        x=xsten(0,1)
        y=xsten(0,2)
        z=xsten(0,SDIM)

        do dir=1,SDIM
         xpos(dir)=xsten(0,dir)
        enddo

        call Box_volumeFAST(bfact,dx,xsten,nhalf,volcell,cencell,SDIM)

        scalc(STATECOMP_PRES+1)=zero

        if (is_in_probtype_list().eq.1) then

         call SUB_LS(xpos,time,distbatch,num_materials)

             ! bcflag=0 (calling from fort_initdata)
         call SUB_STATE(xpos,time,distbatch,local_state, &
                 bcflag,num_materials,num_state_material)
         do im=1,num_materials
          ibase=STATECOMP_STATES+(im-1)*num_state_material
          local_ibase=(im-1)*num_state_material
          scalc(ibase+ENUM_DENVAR+1)= &
            local_state(local_ibase+ENUM_DENVAR+1) ! density
          scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
            local_state(local_ibase+ENUM_TEMPERATUREVAR+1) 
          ! species
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            local_state(local_ibase+num_state_base+n)
           if (scalc(ibase+num_state_base+n).ge.zero) then
            ! do nothing
           else
            print *,"mass fraction cannot be negative"
            stop
           endif
          enddo

          if (scalc(ibase+ENUM_DENVAR+1).gt.zero) then
           ! do nothing
          else
           print *,"density invalid probtype==421 "
           print *,"im,ibase,num_materials ",im,ibase,num_materials
           print *,"density=",scalc(ibase+ENUM_DENVAR+1)
           stop
          endif

          if (scalc(ibase+ENUM_TEMPERATUREVAR+1).gt.zero) then
           ! do nothing
          else
           print *,"temperature invalid probtype==421 "
           print *,"im,ibase,num_materials ",im,ibase,num_materials
           print *,"temperature=",scalc(ibase+ENUM_TEMPERATUREVAR+1)
           stop
          endif

         enddo ! im=1..num_materials
         call SUB_PRES(xpos,time,distbatch,p_hyd,num_materials)
         scalc(STATECOMP_PRES+1)=p_hyd

         if (p_hyd.ge.zero) then
          ! do nothing
         else
          print *,"p_hyd invalid"
          print *,"probtype=",probtype
          print *,"p_hyd=",p_hyd
          stop
         endif

        else if (probtype.eq.401) then

         call HELIX_LS(xpos,time,distbatch)
         call HELIX_STATE(xpos,time,distbatch,local_state)
         do im=1,num_materials
          ibase=STATECOMP_STATES+(im-1)*num_state_material
          local_ibase=(im-1)*num_state_material
          scalc(ibase+ENUM_DENVAR+1)=local_state(local_ibase+ENUM_DENVAR+1) 
          scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
               local_state(local_ibase+ENUM_TEMPERATUREVAR+1) 
           ! species
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            local_state(local_ibase+num_state_base+n)
          enddo
         enddo ! im=1..num_materials
         call HELIX_PRES(xpos,time,distbatch,p_hyd)
         scalc(STATECOMP_PRES+1)=p_hyd

        else if (probtype.eq.402) then

         call TSPRAY_LS(xpos,time,distbatch)
         call TSPRAY_STATE(xpos,time,distbatch,local_state)
         do im=1,num_materials
          ibase=STATECOMP_STATES+(im-1)*num_state_material
          local_ibase=(im-1)*num_state_material
          scalc(ibase+ENUM_DENVAR+1)=local_state(local_ibase+ENUM_DENVAR+1) 
          scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
             local_state(local_ibase+ENUM_TEMPERATUREVAR+1) 
           ! species
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            local_state(local_ibase+num_state_base+n)
          enddo
         enddo ! im=1..num_materials
         call TSPRAY_PRES(xpos,time,distbatch,p_hyd)
         scalc(STATECOMP_PRES+1)=p_hyd

        else if (probtype.eq.412) then ! step

         call CAV2Dstep_LS(xpos,time,distbatch)
         call CAV2Dstep_STATE(xpos,time,distbatch,local_state)
         do im=1,num_materials
          ibase=STATECOMP_STATES+(im-1)*num_state_material
          local_ibase=(im-1)*num_state_material
          scalc(ibase+ENUM_DENVAR+1)=local_state(local_ibase+ENUM_DENVAR+1) 
          scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
             local_state(local_ibase+ENUM_TEMPERATUREVAR+1) 
          ! species
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            local_state(local_ibase+num_state_base+n)
          enddo
         enddo ! im=1..num_materials
         call CAV2Dstep_PRES(xpos,time,distbatch,p_hyd)
         scalc(STATECOMP_PRES+1)=p_hyd

        else if (probtype.eq.413) then ! zeyu

         call ZEYU_droplet_impact_LS(xpos,time,distbatch)
         call ZEYU_droplet_impact_STATE(xpos,time,distbatch,local_state)
         do im=1,num_materials
          ibase=STATECOMP_STATES+(im-1)*num_state_material
          local_ibase=(im-1)*num_state_material
          scalc(ibase+ENUM_DENVAR+1)=local_state(local_ibase+ENUM_DENVAR+1) 
          scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
             local_state(local_ibase+ENUM_TEMPERATUREVAR+1) ! temperature
          ! species
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            local_state(local_ibase+num_state_base+n)
          enddo
         enddo ! im=1..num_materials
         call ZEYU_droplet_impact_PRES(xpos,time,distbatch,p_hyd)
         scalc(STATECOMP_PRES+1)=p_hyd

        else if (probtype.eq.533) then

         call rigid_FSI_LS(xpos,time,distbatch)
         call rigid_FSI_STATE(xpos,time,distbatch,local_state)
         do im=1,num_materials
          ibase=STATECOMP_STATES+(im-1)*num_state_material
          local_ibase=(im-1)*num_state_material
          scalc(ibase+ENUM_DENVAR+1)=local_state(local_ibase+ENUM_DENVAR+1) 
          scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
             local_state(local_ibase+ENUM_TEMPERATUREVAR+1) ! temperature
           ! species
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            local_state(local_ibase+num_state_base+n)
          enddo
         enddo ! im=1..num_materials
         call rigid_FSI_PRES(xpos,time,distbatch,p_hyd)
         scalc(STATECOMP_PRES+1)=p_hyd

        else if (probtype.eq.534) then

         call sinking_FSI_LS(xpos,time,distbatch)
         call sinking_FSI_STATE(xpos,time,distbatch,local_state)
         do im=1,num_materials
          ibase=STATECOMP_STATES+(im-1)*num_state_material
          local_ibase=(im-1)*num_state_material
          scalc(ibase+ENUM_DENVAR+1)=local_state(local_ibase+ENUM_DENVAR+1) 
          scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
              local_state(local_ibase+ENUM_TEMPERATUREVAR+1) ! temperature
           ! species
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            local_state(local_ibase+num_state_base+n)
          enddo
         enddo ! im=1..num_materials
         call sinking_FSI_PRES(xpos,time,distbatch,p_hyd)
         scalc(STATECOMP_PRES+1)=p_hyd

        else if (probtype.eq.311) then ! user defined

         call USERDEF_LS(xpos,time,distbatch)
         call USERDEF_STATE(xpos,time,distbatch,local_state)
         do im=1,num_materials
          ibase=STATECOMP_STATES+(im-1)*num_state_material
          local_ibase=(im-1)*num_state_material
          scalc(ibase+ENUM_DENVAR+1)=local_state(local_ibase+ENUM_DENVAR+1) 
          scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
              local_state(local_ibase+ENUM_TEMPERATUREVAR+1) ! temperature
           ! species
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            local_state(local_ibase+num_state_base+n)
          enddo
         enddo ! im=1..num_materials
         call USERDEF_PRES(xpos,time,distbatch,p_hyd)
         scalc(STATECOMP_PRES+1)=p_hyd

        else

         do im=1,num_materials

          ibase=STATECOMP_STATES+(im-1)*num_state_material

          scalc(ibase+ENUM_DENVAR+1)=fort_denconst(im)  ! den
          scalc(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im) 
     
          do n=1,num_species_var
           scalc(ibase+num_state_base+n)= &
            fort_speciesconst((n-1)*num_materials+im)
          enddo

           ! in: INITDATA
          if (probtype.eq.26) then ! swirl if axis_dir=0 or 1.

           if (axis_dir.eq.10) then ! BCG test
            scalc(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(1)
           else if (axis_dir.eq.11) then ! BCG periodic test
            call get_vortex_info(xpos,time, &
              local_neg_force,local_vel,local_vort, &
              local_energy_moment, &
              scalc(ibase+ENUM_TEMPERATUREVAR+1))
           else if (axis_dir.eq.12) then ! buoyancy test
            T_HOT=fort_initial_temperature(1)
            T_COLD=T_HOT-abs(radblob2)
            if ((T_HOT.gt.T_COLD).and. &
                (T_COLD.gt.zero)) then
             ! do nothing
            else
             print *,"T_HOT or T_COLD invalid"
             stop
            endif

            if (xpos(SDIM).le.problo_array(SDIM)) then
             scalc(ibase+ENUM_TEMPERATUREVAR+1)=T_HOT
            else if (xpos(SDIM).ge.probhi_array(SDIM)) then
             scalc(ibase+ENUM_TEMPERATUREVAR+1)=T_COLD
            else if ((xpos(SDIM).gt.problo_array(SDIM)).and. &
                     (xpos(SDIM).lt.probhi_array(SDIM))) then
             ! the "(T_HOT+T_COLD)/2" temperature contour lives on the
             ! curve z=(problo+probhi)/2+radblob3*cos(beta (x-problo(1)) - pi)
             ! beta=2 pi/(probhi(1)-problo(1))
             ! tanh(x/eps)=1 as x->inf
             ! tanh(x/eps)=-1 as x->-inf
             ! tanh(0)=0
             ! quadratic mapping p(z):
             ! p(zcrit)=0
             ! p(problo)=-1
             ! p(probhi)=1
             ! problo -1 a1=1/(zcrit-problo)      (a2-a1)
             ! zcrit  0  a2=1/(probhi-zcrit) b=-----------------
             ! probhi 1                         probhi-problo
             ! p(z)=-1+a1(z-problo)+b (z-problo)(z-zcrit)
             ! p'(z)=a1+b(2z-problo-zcrit)=a1+2b(z-(problo+zcrit)/2)=0
             ! z=-a1/(2b)+(problo+zcrit)/2
             ! 
             ! f(p(z))=tanh(p(z)/eps)/tanh(1/eps)
             ! g(z)=(f(p(z))+1)/2
             ! T(z)=T_HOT * (1-g(z)) + T_COLD * g(z)
             ! 
             if ((problen_array(1).gt.zero).and. &
                 (problen_array(2).gt.zero).and. &
                 (problen_array(SDIM).gt.zero)) then
              ! do nothing
             else
              print *,"problen_array invalid"
              stop
             endif

             if (levelrz.eq.COORDSYS_CARTESIAN) then
              zcrit=half*(problo_array(SDIM)+probhi_array(SDIM))+ &
               abs(radblob3)* &
               cos(two*Pi*(xpos(1)-problo_array(1))/problen_array(1)-Pi)
             else if ((levelrz.eq.COORDSYS_RZ).and.(SDIM.eq.2)) then
              if (problo_array(1).eq.zero) then
               zcrit=half*(problo_array(SDIM)+probhi_array(SDIM))+ &
                abs(radblob3)* &
                cos(Pi*(xpos(1)+probhi_array(1))/problen_array(1)-Pi)
              else
               print *,"problo_array(1) invalid"
               stop
              endif
             else
              print *,"levelrz invalid probtype=26 and axis_dir=12"
              stop
             endif

             if ((zcrit.gt.problo_array(SDIM)).and. &
                 (zcrit.lt.probhi_array(SDIM))) then
              ! do nothing
             else
              print *,"zcrit invalid"
              stop
             endif
             a1=one/(zcrit-problo_array(SDIM))
             a2=one/(probhi_array(SDIM)-zcrit)
             D2=(a2-a1)/problen_array(SDIM)

             if ((a1.gt.zero).and.(a2.gt.zero)) then
              ! do nothing
             else
              print *,"a1 or a2 invalid"
              stop
             endif

             if (D2.eq.zero) then
              ! do nothing
             else if (D2.ne.zero) then
              z_extrema=half*(-a1/D2+problo_array(SDIM)+zcrit)
              if ((z_extrema.lt.problo_array(SDIM)).or. &
                  (z_extrema.gt.probhi_array(SDIM))) then
               ! do nothing
              else
               print *,"z_extrema invalid: ",z_extrema
               stop
              endif
             else
              print *,"D2 is NaN"
              stop
             endif
             pz=-one+a1*(xpos(SDIM)-problo_array(SDIM))+ &
                D2*(xpos(SDIM)-problo_array(SDIM))* &
                   (xpos(SDIM)-zcrit)
             if ((pz.ge.-one).and.(pz.le.one)) then
              ! do nothing
             else
              print *,"pz out of range"
              stop
             endif
             pz_sanity=-one+a1*(probhi_array(SDIM)-problo_array(SDIM))+ &
                D2*(probhi_array(SDIM)-problo_array(SDIM))* &
                   (probhi_array(SDIM)-zcrit)
             if (abs(pz_sanity-one).le.VOFTOL) then
              ! do nothing
             else
              print *,"pz_sanity invalid"
              stop
             endif
             pz_sanity=-one+a1*(zcrit-problo_array(SDIM))
             if (abs(pz_sanity).le.VOFTOL) then
              ! do nothing
             else
              print *,"pz_sanity invalid"
              stop
             endif


             fpz=tanh(pz/radblob4)/tanh(one/radblob4)
             gpz=half*(fpz+one)
             
             scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
               T_HOT*(one-gpz)+T_COLD*gpz
            else
             print *,"xpos is NaN"
             stop
            endif
           else if ((axis_dir.ge.0).and. & !swirl,vortex confinement
                    (axis_dir.le.5)) then
            doubly_flag=1
            if (SDIM.eq.2) then
             if ((axis_dir.eq.0).or. & ! swirl
                 (axis_dir.eq.1)) then
              rr=y
             else if ((axis_dir.eq.2).or. & ! vortex confinement
                      (axis_dir.eq.3)) then
              rr=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
              doubly_flag=0
             else
              print *,"axis_dir invalid"
              stop
             endif
            else if (SDIM.eq.3) then
             if ((axis_dir.ge.0).and. & !swirl or vortex confinement
                 (axis_dir.le.3)) then
              if (adv_dir.eq.3) then ! vortex confinement
               rr=y
              else if (adv_dir.eq.2) then ! vortex confinement
               rr=z
              else if (adv_dir.eq.1) then ! swirl
               rr=z
              else
               print *,"adv_dir invalid probtype==26 (11)"
               stop
              endif
             else if ((axis_dir.eq.4).or. & ! vortex confinement.
                      (axis_dir.eq.5)) then
              rr=sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)-radblob
              doubly_flag=0
             else
              print *,"axis_dir invalid"
              stop
             endif
            else
             print *,"dimension bust"
             stop
            endif
            if (doubly_flag.eq.1) then
             if (rr.le.half) then
              jumpval=tanh( (rr-one/four)*30.0 )
             else
              jumpval=tanh( (three/four-rr)*30.0 )
             endif
            else if (doubly_flag.eq.0) then
             jumpval=tanh(30.0*rr)
            else
             print *,"doubly_flag invalid"
             stop
            endif

            jumpval=(jumpval+one)/two
            scalc(ibase+ENUM_TEMPERATUREVAR+1)= &
             jumpval*fort_initial_temperature(1)+ &
             (one-jumpval)*fort_initial_temperature(2)
  
           else
            print *,"axis_dir invalid"
            stop
           endif

          endif ! probtype==26

            ! 2B from Wardlaws list
          if ((probtype.eq.36).and.(axis_dir.eq.2).and. &
              (SDIM.eq.2)) then
           if (im.eq.1) then
             !calling from fort_initdata
            call tait_hydrostatic_pressure_density(xpos, &
             rhohydro,preshydro,from_boundary_hydrostatic)
            scalc(ibase+ENUM_DENVAR+1)=rhohydro
           else if (im.eq.2) then
            e_jwl=4.2814D+10
            den_jwl=1.63D0
            call init_massfrac_parm(den_jwl,massfrac_parm,im)
            call TEMPERATURE_material(den_jwl,massfrac_parm, &
             temp_jwl,e_jwl, &
             fort_material_type(im),im)
            scalc(ibase+ENUM_DENVAR+1)=den_jwl
            scalc(ibase+ENUM_TEMPERATUREVAR+1)=temp_jwl
           endif
          endif  ! probtype=36

           ! initial temperature for melting ice block on a substrate.
          if (probtype.eq.59) then

           if (num_materials.lt.4) then
            print *,"num_materials too small for melting ice on substrate"
            stop
           endif
           ! substrate (initial temperature)
           if (im.eq.4) then
            ! bcflag=0 (calling from fort_initdata)
            call outside_temperature(time,x,y,z,water_temp,im,0)
            scalc(ibase+ENUM_TEMPERATUREVAR+1)=water_temp
           endif

          endif ! probtype.eq.59

           ! Benard instability problem initdata
           ! density above is a filler; density will be
           ! replaced by rho(T,z) after "nonlinear_advection"
           ! (correct_density)
          if (probtype.eq.603) then
           if (z.gt.yblob) then
            water_temp=fort_initial_temperature(1) 
           else if (z.gt.zero) then 
            temp_slope=-radblob2/yblob
            if (SDIM.eq.2) then
             eps_benard=radblob*cos(two*Pi*x/xblob)*four*y*(yblob-y)/(yblob**2) 
            else if (SDIM.eq.3) then
             eps_benard=cos(two*Pi*(x-problox)/xblob)* &
                cos(two*Pi*(y-probloy)/xblob)* &
                radblob*four*z*(yblob-z)/(yblob**2)
            else
             print *,"dimension bust"
             stop
            endif

            water_temp=radblob2+fort_initial_temperature(1)+ &
              temp_slope*(z+eps_benard)
           else
            water_temp=radblob2+fort_initial_temperature(1)
           endif
           scalc(ibase+ENUM_TEMPERATUREVAR+1)=water_temp  ! temperature
          endif  ! probtype=603

           ! shock tube problems
           ! do not use material_type=5 (EOS_air),
           ! use material_type=18 instead.
           ! (results should be similar though)
          if ((probtype.eq.92).or. & !contact captured
              (probtype.eq.93)) then !contact tracked
           gamma_jwl=1.4

           if (axis_dir.eq.0) then  ! Sod shock tube
            den_jwl_left=one
            den_jwl_right=0.125d0
            p_jwl_left=one
            p_jwl_right=0.1d0
            u_jwl_left=zero
            u_jwl_right=zero
            xshock=half
           else if (axis_dir.eq.1) then ! strong shock tube
            den_jwl_left=one
            den_jwl_right=0.125d0
            p_jwl_left=1.0D+10
            p_jwl_right=0.1d0
            u_jwl_left=zero
            u_jwl_right=zero
            xshock=half
           else if (axis_dir.eq.2) then ! shock turbulence interaction
            den_jwl_left=3.857148d0
            den_jwl_right=one+0.2d0*sin(five*x-five)
            p_jwl_left=10.333333d0
            p_jwl_right=one
            u_jwl_left=2.629369d0
            u_jwl_right=zero
            xshock=one
           else if (axis_dir.eq.3) then ! mach>4
            den_jwl_left=10.0d0
            den_jwl_right=1.0d0
            p_jwl_left=10.0d0*(1.4d0-1.0d0)
            p_jwl_right=(1.4d0-1.0d0)
            u_jwl_left=5.0d0
            u_jwl_right=zero
            xshock=one
            ! Kadioglu, Sussman, Osher, Wright, Kang (smooth test problem)
           else if (axis_dir.eq.4) then 
            u_jwl_left=zero
            u_jwl_right=zero
            if (adv_dir.eq.1) then
             rr=x
            else if (adv_dir.eq.2) then
             rr=y
            else if ((adv_dir.eq.3).and.(SDIM.eq.3)) then
             rr=z
            else
             print *,"adv_dir invalid probtype==92,93 (12)"
             stop
            endif
            p_jwl_left=(1.0D+6)+60.0d0*cos(two*Pi*rr)+100.0d0*sin(four*Pi*rr)
            p_jwl_right=p_jwl_left
            den_jwl_left=fort_denconst(2)*((p_jwl/1.0D+6)**(one/gamma_jwl))
            den_jwl_right=den_jwl_left
            xshock=one
           else 
            print *,"axis_dir invalid probtype=92 or 93: ",axis_dir
            stop
           endif
           if ((axis_dir.eq.0).or. &
               (axis_dir.eq.1).or. &
               (axis_dir.eq.2).or. &
               (axis_dir.eq.3)) then
            e_jwl_left=p_jwl_left/((gamma_jwl-one)*den_jwl_left)
            e_jwl_right=p_jwl_right/((gamma_jwl-one)*den_jwl_right)
            call init_massfrac_parm(den_jwl_left,massfrac_parm,im)
            call TEMPERATURE_material(den_jwl_left,massfrac_parm, &
             temp_jwl_left,e_jwl_left, &
             fort_material_type(im),im)
            call init_massfrac_parm(den_jwl_right,massfrac_parm,im)
            call TEMPERATURE_material(den_jwl_right,massfrac_parm, &
             temp_jwl_right, &
             e_jwl_right,fort_material_type(im),im)
           else if (axis_dir.eq.4) then
            temp_jwl_left=fort_initial_temperature(1)
            temp_jwl_right=temp_jwl_left
           else
            print *,"axis_dir invalid"
            stop
           endif
           if (probtype.eq.92) then!the contact discontinuity is captured.
            if (x.le.xshock) then
             den_jwl=den_jwl_left 
             temp_jwl=temp_jwl_left 
            else if (x.gt.xshock) then
             den_jwl=den_jwl_right
             temp_jwl=temp_jwl_right
            else
             print *,"x invalid"
             stop
            endif
           else if (probtype.eq.93) then!the contact discontinuity is tracked.
            if (xblob.lt.xshock) then
             if (im.eq.1) then
              den_jwl=den_jwl_left
              temp_jwl=temp_jwl_left 
             else if (im.eq.2) then
              if (x.le.xshock) then
               den_jwl=den_jwl_left
               temp_jwl=temp_jwl_left 
              else if (x.gt.xshock) then
               den_jwl=den_jwl_right
               temp_jwl=temp_jwl_right
              else
               print *,"x invalid"
               stop
              endif
             else
              print *,"im invalid84"
              stop
             endif
            else if (xblob.gt.xshock) then
             if (im.eq.1) then
              if (x.le.xshock) then
               den_jwl=den_jwl_left
               temp_jwl=temp_jwl_left 
              else if (x.gt.xshock) then
               den_jwl=den_jwl_right
               temp_jwl=temp_jwl_right
              else
               print *,"x invalid"
               stop
              endif
             else if (im.eq.2) then
              den_jwl=den_jwl_right
              temp_jwl=temp_jwl_right
             else
              print *,"im invalid85"
              stop
             endif
            else if (xblob.eq.xshock) then 
             if (im.eq.1) then
              den_jwl=den_jwl_left
             else if (im.eq.2) then
              den_jwl=den_jwl_right
             else
              print *,"im invalid86"
              stop
             endif
             if (x.le.xshock) then
              temp_jwl=temp_jwl_left 
             else if (x.gt.xshock) then
              temp_jwl=temp_jwl_right
             else
              print *,"x invalid"
              stop
             endif
            else
             print *,"xblob or xshock invalid"
             stop
            endif
           else 
            print *,"probtype invalid"
            stop
           endif

           scalc(ibase+ENUM_DENVAR+1)=den_jwl
           scalc(ibase+ENUM_TEMPERATUREVAR+1)=temp_jwl

          endif ! probtype=92 or 93 (shock tube test problems)

           ! in: fort_initdata
           ! material=1 is TAIT EOS
          if (fort_material_type(1).eq.13) then
           if (im.eq.1) then
             !calling from fort_initdata
            call tait_hydrostatic_pressure_density(xpos, &
             rhohydro,preshydro,from_boundary_hydrostatic)
            scalc(ibase+ENUM_DENVAR+1)=rhohydro

           endif  ! im=1
          endif ! material_type(1)=13

           ! 2E from Wardlaw's list (bubble jetting) 
           ! in: initdata
          if ((probtype.eq.42).and.(axis_dir.eq.1)) then
           if (im.eq.1) then

             !calling from fort_initdata
            call tait_hydrostatic_pressure_density(xpos,rhohydro,preshydro, &
                    from_boundary_hydrostatic)
            scalc(ibase+ENUM_DENVAR+1)=rhohydro
   
           else if (im.eq.2) then
            e_jwl=4.2945D+10
            den_jwl=1.63D0
            call init_massfrac_parm(den_jwl,massfrac_parm,im)
            call TEMPERATURE_material(den_jwl,massfrac_parm, &
             temp_jwl,e_jwl, &
             fort_material_type(im),im)

            scalc(ibase+ENUM_DENVAR+1)=den_jwl     ! density
            scalc(ibase+ENUM_TEMPERATUREVAR+1)=temp_jwl
           endif
          endif  ! probtype=42 axis_dir=1 (bubble jetting)

          ! Marangoni (heat pipe) test problem
          ! flag=0 
          if ((probtype.eq.36).and.(axis_dir.eq.10)) then
           call position_Temp(0,radblob,radblob2,x,y,z,temp_jwl)
           scalc(ibase+ENUM_TEMPERATUREVAR+1)=temp_jwl
          endif

           ! in: initdata
          if (probtype.eq.46) then ! cavitation

           if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
            if (im.eq.1) then ! water

              !calling from fort_initdata
             call tait_hydrostatic_pressure_density(xpos,rhohydro,preshydro, &
                     from_boundary_hydrostatic)
             scalc(ibase+ENUM_DENVAR+1)=rhohydro

            else if (im.eq.2) then ! jwl
             e_jwl=4.2945D+10
             den_jwl=1.63D0
             call init_massfrac_parm(den_jwl,massfrac_parm,im)
             call TEMPERATURE_material(den_jwl,massfrac_parm, &
              temp_jwl,e_jwl, &
              fort_material_type(im),im)

             scalc(ibase+ENUM_DENVAR+1)=den_jwl     ! density
             scalc(ibase+ENUM_TEMPERATUREVAR+1)=temp_jwl
            else if (im.eq.3) then  ! air
             call general_hydrostatic_pressure(p_hyd)
             den_jwl=fort_denconst(im)
             temp_jwl=fort_initial_temperature(im)
             call init_massfrac_parm(den_jwl,massfrac_parm,im)
             call INTERNAL_material(den_jwl,massfrac_parm, &
              temp_jwl,e_jwl, &
              fort_material_type(im),im)
             call EOS_material(den_jwl,massfrac_parm, &
              e_jwl,p_jwl, &
              fort_material_type(im),im)
             temp_jwl=temp_jwl*p_hyd/p_jwl
        
             scalc(ibase+ENUM_DENVAR+1)=den_jwl
             scalc(ibase+ENUM_TEMPERATUREVAR+1)=temp_jwl

            endif
           else if (axis_dir.eq.10) then
            if (im.eq.1) then ! water
              !calling from fort_initdata
             call tait_hydrostatic_pressure_density(xpos,rhohydro,preshydro, &
                     from_boundary_hydrostatic)
             scalc(ibase+ENUM_DENVAR+1)=rhohydro
            endif
           else if (axis_dir.eq.20) then
            !do nothing (CODY ESTEBE created test problem,fort_denconst(im) ok)
           else
            print *,"axis_dir invalid"
            stop
           endif
          endif  ! probtype=46 (cavitation)


           ! (initdata) shock injection with nozzle and pressure BC
          if ((probtype.eq.53).and.(axis_dir.eq.2)) then

           denroom=fort_denconst(im)
           temproom=fort_initial_temperature(im)  ! room temp

           if (im.eq.1) then  ! liquid
            scalc(ibase+ENUM_DENVAR+1)=denroom
           else if (im.eq.2) then  ! gas
            scalc(ibase+ENUM_DENVAR+1)=denroom
            imattype=fort_material_type(im)

            if (imattype.eq.0) then
             ! do nothing
            else if (imattype.gt.0) then

             call init_massfrac_parm(denroom,massfrac_parm,im)
             call INTERNAL_material(denroom,massfrac_parm, &
              temproom,e_room, &
              imattype,im)
             call EOS_material(denroom,massfrac_parm, &
               e_room,p_room,imattype,im)
             call general_hydrostatic_pressure(p_hyd)
             temproom=temproom*p_hyd/p_room
             e_room=e_room*p_hyd/p_room
             scalc(ibase+ENUM_TEMPERATUREVAR+1)=temproom  ! temperature
            else 
             print *,"imattype invalid fort_initdata"
             stop
            endif

           endif

           ! shock injection JICF with compressible gas.
          else if ((probtype.eq.53).and.(fort_material_type(2).gt.0)) then

           call general_hydrostatic_pressure(p_hyd)
           scalc(STATECOMP_PRES+1)=p_hyd

          else if ((probtype.eq.530).and. &
                   (axis_dir.eq.1).and. &
                   (fort_material_type(2).gt.0).and. &
                   (SDIM.eq.3)) then
           call general_hydrostatic_pressure(p_hyd)
           scalc(STATECOMP_PRES+1)=p_hyd
          endif

          ! circular freezing disk
          ! or spherical boiling
          if ((probtype.eq.801).and.(axis_dir.eq.3)) then 
   
           im_source=1
           im_dest=2
           ireverse=0
           call get_iten(im_source,im_dest,iten)
           L_ice_melt= &
            abs(get_user_latent_heat(iten+ireverse*num_interfaces, &
                 room_temperature,1))
           TSAT=saturation_temp(iten+ireverse*num_interfaces)
           T_EXTREME=fort_initial_temperature(im_source)
           cp_melt=get_user_stiffCP(im_source) ! J/Kelvin
           k_melt=get_user_heatviscconst(im_source) ! W/(m Kelvin)
  
           if (SDIM.eq.2) then
            rstefan=sqrt((x-xblob)**2+(y-yblob)**2)
           else
            rstefan=sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)
           endif
           if (rstefan.le.radblob) then
            T_FIELD=TSAT
           else
            den_ratio=max(fort_denconst(im_dest),fort_denconst(im_source))/ &
                      min(fort_denconst(im_dest),fort_denconst(im_source))

            if (den_ratio.lt.10.0d0) then 
             call liquid_temperature( &
              fort_beta(iten+ireverse*num_interfaces), & ! lmSt
              T_EXTREME, &
              L_ice_melt, &
              cp_melt, &
              fort_stefan_number(iten+ireverse*num_interfaces), &
              rstefan, &
              fort_time_radblob(iten+ireverse*num_interfaces), &
              k_melt, &
              T_FIELD)
            else if (den_ratio.ge.10.0) then
             call superheat_temperature( &
              fort_alpha(iten+ireverse*num_interfaces), &
              fort_beta(iten+ireverse*num_interfaces), &
              den_ratio, &
              T_EXTREME, &
              fort_jacob_number(iten+ireverse*num_interfaces), &
              TSAT, &
              rstefan, &
              fort_time_radblob(iten+ireverse*num_interfaces), &
              T_FIELD)
            else
             print *,"for_expansion_factor invalid"
             stop
            endif
           endif
           scalc(ibase+ENUM_TEMPERATUREVAR+1)=T_FIELD

          endif ! (probtype.eq.801).and.(axis_dir.eq.3)

          if (probtype.eq.802) then ! dissolution
           scalc(ibase+ENUM_TEMPERATUREVAR+1)=two ! T (concentration_initdata)
           call vapordist(xsten,nhalf,dx,bfact,posdiss_initdata)
           if (posdiss_initdata.le.zero) then
            concentration_initdata=two
           else if (posdiss_initdata.ge.xhidiss_initdata-two*dx(SDIM)) then
            concentration_initdata=one
           else
            ispace_initdata=NINT(posdiss_initdata/dx(SDIM)-half)
            if (ispace_initdata.eq.0) then
             concen1_initdata=two
             concen2_initdata=soln_initdata(ispace_initdata+1)+one
            else if (ispace_initdata.ge.ndiss_initdata-2) then
             concen1_initdata=soln_initdata(ispace_initdata)+one
             concen2_initdata=one
            else
             concen1_initdata=soln_initdata(ispace_initdata)+one
             concen2_initdata=soln_initdata(ispace_initdata+1)+one
            endif
            theta_initdata=(posdiss_initdata-ispace_initdata*dx(SDIM))/dx(SDIM)
            concentration_initdata= &
               (one-theta_initdata)*concen1_initdata+ &
               theta_initdata*concen2_initdata
           endif
            ! T (concentration_initdata)
           scalc(ibase+ENUM_TEMPERATUREVAR+1)=concentration_initdata   
          endif ! 802 (dissolution)

           ! in: subroutine fort_initdata
           ! hydrates
          if (probtype.eq.199) then
           if (num_materials.ne.3) then
            print *,"3 materials for hydrate reactor"
            stop
           endif
           if (num_species_var.ne.1) then
            print *,"num_species_var should be 1"
            stop
           endif
           if (im.eq.1) then
            call INIT_STATE_WATER(x,y,z,time,vel,temp,dens,ccnt) 
           else if (im.eq.2) then
            call INIT_STATE_GAS(x,y,z,time,vel,temp,dens,ccnt) 
           else if (im.eq.3) then
            call INIT_STATE_HYDRATE(x,y,z,time,vel,temp,dens,ccnt) 
           else
            print *,"im invalid87"
            stop
           endif

            ! density comes from the inputs file.
           scalc(ibase+ENUM_TEMPERATUREVAR+1)=temp
           scalc(ibase+3)=ccnt

           ! in: subroutine fort_initdata
          else if (probtype.eq.220) then
           ! do nothing, density and temperature are set from the input file
           ! in the beginning of the loop on im. NO INIT_STATE_*** is called

           ! in: subroutine fort_initdata
          else if ((probtype.eq.299).or. &
                   (probtype.eq.301)) then !melting (initial temperature field)

           temp=fort_initial_temperature(im)
           scalc(ibase+ENUM_TEMPERATUREVAR+1)=temp

          endif  

         enddo  ! im=1..num_materials

        endif ! if (probtype.eq.user_def_probtype) then ... else ... endif

        call materialdist_batch(xsten,nhalf,dx,bfact,distbatch,time)
        do im=1,num_materials
          !FSI_PRESCRIBED_PROBF90
          !FSI_PRESCRIBED_NODES
          !FSI_SHOELE_CTML
         if (is_rigid(im).eq.1) then

          if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. &
              (FSI_flag(im).eq.FSI_SHOELE_CTML)) then
           distbatch(im)=LS(D_DECL(ic,jc,kc),im)
          else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
           ! do nothing
          else
           print *,"FSI_flag(im) invalid in fort_initdata"
           print *,"im,FSI_flag(im) ",im,FSI_flag(im)
           stop
          endif

         else if (is_rigid(im).eq.0) then
          ! do nothing
         else
          print *,"is_rigid(im) invalid"
          stop
         endif
        enddo ! im=1..num_materials

         ! in: fort_initdata
         ! "stackvolume_batch" is declared in STACKVOLUME.F90
         ! "level" = 0
        call stackvolume_batch(xsten,nhalf,dx,bfact,fluiddata, &
         0,max_levelstack,materialdist_batch,time)
        call extract_vof_cen_batch(fluiddata,vofdark,voflight, &
         cendark,cenlight)

         ! Rayleigh-Taylor, checkerboard test
        if (probtype.eq.602) then
          ! in subroutine vapordist: dist=y-(yblob+(1/16)/2)
          ! material 1 is on top (dist>0)
          ! material 2 is on bottom (dist<0)
         if ((xblob.ge.1.0D+10).and.(radblob.le.EPS10)) then
          if (num_materials.ne.2) then
           print *,"num_materials invalid"
           stop
          endif

          if (1.eq.0) then
           print *,"ic,jc,x,y ",ic,jc,x,y
           do im=1,2
            print *,"im,dist,vof,cen ",im,distbatch(im),vofdark(im), &
              cendark(im,1),cendark(im,2)
           enddo
          endif

          if (abs(1.0d0/16.0d0-dx(2)).le.EPS14) then
           ! do nothing
          else
           print *,"checkerboard test distance function assumes dy=1/16"
           stop
          endif
          if (abs(1.0d0/16.0d0-dx(1)).le.EPS14) then
           ! do nothing
          else
           print *,"checkerboard test distance function assumes dx=1/16"
           stop
          endif

           ! if 1x1 centroid is (alpha,-1/4+beta) then: 
           !  for saw tooth:
           !  3x3 volume=3(3/2)=9/2
           !  3x3 bot centroidx:(1/2)((-1+alpha)+alpha+(1+alpha))/(9/2)=alpha/3
           !  3x3 bot centroidy: ((3/2)(-1/4+beta)-3)/(9/2)= 
           !  (-3/8)/(9/2)-2/3+beta/3=-1/12-8/12+beta/3=-3/4+beta/3
           !  for continuation as line:
           !  3x3 volume=9/2
           !  3x3 bot centroidx: 3 alpha
           !  3x3 bot centroidy: -3/4+3 beta
           !  3 beta' = beta/3   beta'=beta/9=(m/3)^2/12
           !  3 alpha' = alpha/3   alpha'=alpha/9=(m/9)/6
          dir=2
          im=1 ! material 1 on top
          cendark(im,dir)=(0.25d0-(slope_checker**2)/12.0d0)*dx(2)
          im=2 ! material 2 on bottom
          cendark(im,dir)=-cendark(1,dir)
          dir=1
          im=1 ! material 1 on top
          cendark(im,dir)=-slope_checker*dx(1)/6.0d0 
          im=2 ! material 2 on bottom
          cendark(im,dir)=-cendark(1,dir)
         endif
        endif ! Rayleigh-Taylor checkerboard test

        debug_vfrac_sum=zero

        do imls=1,num_materials
         LSc(imls)=distbatch(imls)
        enddo

        do im=1,num_materials

         vofcomp_raw=STATECOMP_MOF+(im-1)*ngeom_raw+1

         scalc(vofcomp_raw)=vofdark(im)
          ! centroid relative to centroid of cell; not cell center.
         do dir=1,SDIM
          scalc(vofcomp_raw+dir)=cendark(im,dir)
         enddo

         if (is_rigid(im).eq.0) then
          debug_vfrac_sum=debug_vfrac_sum+vofdark(im)
         else if (is_rigid(im).eq.1) then

          if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. &
              (FSI_flag(im).eq.FSI_SHOELE_CTML)) then
           scalc(vofcomp_raw)=scal(D_DECL(ic,jc,kc),vofcomp_raw)
           do dir=1,SDIM 
            scalc(vofcomp_raw+dir)=scal(D_DECL(ic,jc,kc),vofcomp_raw+dir)
           enddo
          else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
           ! do nothing
          else
           print *,"FSI_flag invalid in fort_initdata"
           print *,"im,FSI_flag(im): ",im,FSI_flag(im)
           stop
          endif
         else
          print *,"is_rigid(im) invalid"
          stop
         endif

        enddo  ! im=1..num_materials
        
        if (debug_vfrac_sum.le.half) then
         print *,"WARNING in 'process_initdata'"
         print *,"debug_vfrac_sum= ",debug_vfrac_sum
         print *,"ic,jc,kc,x,y,z ",ic,jc,kc,x,y,z
         print *,"time=",time
         do im=1,num_materials
          print *,"im,vofdark ",im,vofdark(im)
          print *,"im,distbatch ",im,distbatch(im)
         enddo
         if ((im_solid_initdata.ge.1).and. &
             (im_solid_initdata.le.num_materials)) then
          call materialdistsolid(x,y,z,distsolid,time,im_solid_initdata)
          if ((FSI_flag(im_solid_initdata).eq.FSI_PRESCRIBED_NODES).or. &
              (FSI_flag(im_solid_initdata).eq.FSI_SHOELE_CTML)) then
           distsolid=LS(D_DECL(ic,jc,kc),im_solid_initdata)
          else if (FSI_flag(im_solid_initdata).eq.FSI_PRESCRIBED_PROBF90) then
           ! do nothing
          else
           print *,"FSI_flag invalid in fort_initdata"
           print *,"im_solid_initdata,FSI_flag(im_solid_initdata) ", &
            im_solid_initdata,FSI_flag(im_solid_initdata)
           stop
          endif
         else
          print *,"im_solid_initdata invalid"
          stop
         endif

         print *,"result of materialdistsolid: distsolid=",distsolid

         print *,"adjusting the volume fraction of the 1st material"
         im=1
         if (is_rigid(im).ne.0) then
          print *,"is_rigid(im).ne.0"
          stop
         endif
         vofdark(im)=vofdark(im)+one-debug_vfrac_sum
         vofcomp_raw=STATECOMP_MOF+(im-1)*ngeom_raw+1
         scalc(vofcomp_raw)=vofdark(im)
        endif

        do k1=k1lo,k1hi
        do j1=-1,1
        do i1=-1,1

         do isten=-nhalf2,nhalf2
          dir=1
          xsten2(isten,dir)=xsten(isten+2*i1,dir)
          dir=2
          xsten2(isten,dir)=xsten(isten+2*j1,dir)
          if (SDIM.eq.3) then
           dir=SDIM
           xsten2(isten,dir)=xsten(isten+2*k1,dir)
          endif
         enddo ! isten
         call materialdist_batch(xsten2,nhalf2,dx,bfact,distbatch,time)
         do im=1,num_materials
          if (is_rigid(im).eq.1) then

           if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. &
               (FSI_flag(im).eq.FSI_SHOELE_CTML)) then
            distbatch(im)=LS(D_DECL(ic+i1,jc+j1,kc+k1),im)
           else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
            ! do nothing
           else
            print *,"FSI_flag(im) invalid in fort_initdata"
            print *,"im,FSI_flag(im) ",im,FSI_flag(im)
            stop
           endif
          else if (is_rigid(im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid invalid PROB.F90"
           stop
          endif 
          LS_stencil(D_DECL(i1,j1,k1),im)=distbatch(im)
         enddo

        enddo
        enddo
        enddo ! i1,j1,k1 = -1..1

        do im=1,num_materials
         vofcomp_recon=(im-1)*ngeom_recon+1
         vofcomp_raw=STATECOMP_MOF+(im-1)*ngeom_raw+1
         mofdata(vofcomp_recon)=scalc(vofcomp_raw)
         mofdata(vofcomp_recon+SDIM+1)=zero ! order
         mofdata(vofcomp_recon+2*SDIM+2)=zero ! intercept

         do dir=1,SDIM

          if (level.eq.max_level) then
           if (centroid_noise_factor(im).eq.zero) then
            centroid_noise=zero
           else if (centroid_noise_factor(im).gt.zero) then
             ! 0.0<=noise_amplitude<=1.0
            Call random_number(noise_amplitude)
            centroid_noise=(two*noise_amplitude-one)* &
               centroid_noise_factor(im)*dx(1)
           else
            print *,"centroid_noise_factor invalid"
            stop
           endif
          else if ((level.lt.max_level).and. &
                   (level.ge.0)) then
           centroid_noise=zero
          else
           print *,"level invalid"
           stop
          endif
          
           ! centroid relative to centroid of cell; not cell center.
          mofdata(vofcomp_recon+dir)=scalc(vofcomp_raw+dir)+centroid_noise
          mofdata(vofcomp_recon+SDIM+dir+1)=zero ! slope
         enddo  !dir=1..SDIM
        enddo  ! im=1..num_materials

        ! sum F_fluid=1  sum F_solid <= 1
        ! centroids are projected to the cell in question.
        call make_vfrac_sum_ok_base( &
          cmofsten, &
          xsten,nhalf, &
          continuous_mof, &
          bfact,dx, &
          tessellate, & ! =0
          mofdata,SDIM)

        do im=1,num_materials
         vofcomp_recon=(im-1)*ngeom_recon+1
         voflist(im)=mofdata(vofcomp_recon)
        enddo

        call calc_error_indicator( &
         level,max_level, &
         xsten,nhalf,dx,bfact, &
         voflist, &
         LS_stencil, &
         err,time)

        scalc(nc)=err

        if (volcell.gt.zero) then
         ! do nothing
        else
         print *,"volcell invalid in INITDATA: ",volcell
         stop
        endif 

        vfracsum_test=zero
        do im=1,num_materials
         vofcomp_raw=STATECOMP_MOF+(im-1)*ngeom_raw+1
         if (is_rigid(im).eq.0) then
          vfracsum_test=vfracsum_test+scalc(vofcomp_raw)
         else if (is_rigid(im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(im) invalid"
          stop
         endif
        enddo !im=1..num_materials

        if ((vfracsum_test.gt.half).and. &
            (vfracsum_test.le.1.5d0)) then
         ! do nothing
        else
         print *,"FAILED: vfracsum_test= ",vfracsum_test
         do im=1,num_materials
          vofcomp_raw=STATECOMP_MOF+(im-1)*ngeom_raw+1
          print *,"im,vfrac ",im,scalc(vofcomp_raw)
          print *,"im,LS ",im,LSc(im)
         enddo
         print *,"level,max_level ",level,max_level
         print *,"x,y,z= ",x,y,z
         print *,"dx,dy,dz= ",dx(1),dx(2),dx(SDIM)
         print *,"ic,jc,kc= ",ic,jc,kc
         print *,"radblob,radblob2,radblob3,radblob4 ", &
           radblob,radblob2,radblob3,radblob4
         print *,"xblob,yblob,zblob ",xblob,yblob,zblob
         print *,"probtype ",probtype
         stop
        endif

        do imls=1,num_materials 

         vofcomp_recon=(imls-1)*ngeom_recon+1
         do dir=1,SDIM
          centroid_absolute(dir)=cencell(dir)+ &
             mofdata(vofcomp_recon+dir)
         enddo

         call find_cut_geom_slope_CLSVOF( &
          continuous_mof, & !STANDARD_MOF
          LS_stencil, &
          particle_list, &
          num_particles, &
          lsnormal, &
          lsnormal_valid, &
          ls_intercept, &
          bfact,dx, &
          xsten,nhalf, &
          centroid_absolute, &
          imls, &
          dxmaxLS, &
          SDIM)

         if (lsnormal_valid(imls).eq.1) then
          do dir=1,SDIM
           LSc(num_materials+SDIM*(imls-1)+dir)=lsnormal(imls,dir)
          enddo
         else if (lsnormal_valid(imls).eq.0) then
          ! do nothing
         else
          print *,"lsnormal_valid invalid"
          stop
         endif
        enddo !imls=1..num_materials 

        do im=1,num_materials

         vofcomp_raw=STATECOMP_MOF+(im-1)*ngeom_raw+1
         vofcomp_recon=(im-1)*ngeom_recon+1

         if (abs(voflist(im)-mofdata(vofcomp_recon)).le.1.0E-12) then
          ! do nothing
         else
          print *,"voflist and mofdata mismatch"
          stop
         endif

         scalc(vofcomp_raw)=voflist(im)
         do dir=1,SDIM
          scalc(vofcomp_raw+dir)=mofdata(vofcomp_recon+dir)
         enddo

        enddo ! im=1..num_materials

         ! nc=STATE_NCOMP
        do n=1,nc
         scal(D_DECL(ic,jc,kc),n)=scalc(n)
        enddo

        if ((num_materials_compressible.ge.1).and. &
            (num_materials_compressible.le.num_materials)) then
         kfine=0
#if (AMREX_SPACEDIM==3)
         do kfine=0,1
#endif
         do jfine=0,1
         do ifine=0,1
          nfine=4*kfine+2*jfine+ifine+1

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
            refine_comp=ENUM_NUM_REFINE_DENSITY_TYPE* &
                   (im_refine_density-1)+nfine
            ibase=STATECOMP_STATES+(im-1)*num_state_material+ &
                    ENUM_DENVAR+1
            refineden(D_DECL(ic,jc,kc),refine_comp)=scalc(ibase)
           else
            print *,"is_compressible_mat(im) invalid"
            stop
           endif
          enddo ! im=1..num_materials

         enddo !ifine
         enddo !jfine
#if (AMREX_SPACEDIM==3)
         enddo !kfine
#endif

        else if (num_materials_compressible.eq.0) then
         !do nothing
        else
         print *,"num_materials_compressible invalid"
         stop
        endif

        do imls=1,num_materials*(1+SDIM)
         LS(D_DECL(ic,jc,kc),imls)=LSc(imls)
        enddo

       enddo
       enddo
       enddo ! ic,jc,kc

       return
       end subroutine fort_initdata

       subroutine fort_addnoise( &
        dir, &
        angular_velocity, & !INTENT(in): fort_addnoise
        perturbation_mode, &
        perturbation_eps_temp, &
        perturbation_eps_vel, &
        nstate, &
        xlo,dx,  &
        Snew,DIMS(Snew), &
        LSnew,DIMS(LSnew), &
        MAC,DIMS(MAC), &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        level, &
        finest_level) &
       bind(c,name='fort_addnoise')

       use global_utility_module

       IMPLICIT NONE

      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: angular_velocity !fort_addnoise
      integer, INTENT(in) :: perturbation_mode
      real(amrex_real), INTENT(in) :: perturbation_eps_temp
      real(amrex_real), INTENT(in) :: perturbation_eps_vel
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: DIMDEC(Snew)
      integer, INTENT(in) :: DIMDEC(LSnew)
      integer, INTENT(in) :: DIMDEC(MAC)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(inout),target :: Snew(DIMV(Snew),nstate)
      real(amrex_real), pointer :: Snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: LSnew(DIMV(LSnew),num_materials)
      real(amrex_real), pointer :: LSnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: MAC(DIMV(MAC))
      real(amrex_real), pointer :: MAC_ptr(D_DECL(:,:,:))


      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer i,j,k,ii,jj,kk,dir2
      real(amrex_real) problo_arr(SDIM)
      real(amrex_real) probhi_arr(SDIM)
      real(amrex_real) sinprod

      if (bfact.lt.1) then
       print *,"bfact invalid200"
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
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid add noise"
       stop
      endif

      problo_arr(1)=problox
      problo_arr(2)=probloy
      probhi_arr(1)=probhix
      probhi_arr(2)=probhiy
      if (SDIM.eq.3) then
       problo_arr(SDIM)=probloz
       probhi_arr(SDIM)=probhiz
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
       print *,"dir invalid in addnoise "
       stop
      endif

      Snew_ptr=>Snew
      LSnew_ptr=>LSnew
      MAC_ptr=>MAC

      call checkbound_array(fablo,fabhi,Snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LSnew_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,MAC_ptr,0,dir)

      if (perturbation_mode.le.0) then
       print *,"perturbation_mode invalid"
       stop
      endif 
      if (perturbation_mode.gt.1024) then
       print *,"perturbation_mode too large"
       stop
      endif 
      if (perturbation_eps_temp.lt.zero) then
       print *,"perturbation_eps_temp invalid"
       stop
      endif 
      if (perturbation_eps_vel.lt.zero) then
       print *,"perturbation_eps_vel invalid"
       stop
      endif 

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       sinprod=one
       do dir2=1,SDIM
        if (probhi_arr(dir2).le.problo_arr(dir2)) then
         print *,"probhi_arr invalid"
         stop
        endif
        sinprod=sinprod*sin(two*Pi*perturbation_mode* &
          (xsten(0,dir2)-problo_arr(dir2))/ &
          (probhi_arr(dir2)-problo_arr(dir2)))
       enddo ! dir2

      enddo
      enddo
      enddo ! i,j,k

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
        ! dir=0..sdim-1
       call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir)

       sinprod=one
       do dir2=1,SDIM
        if (probhi_arr(dir2).le.problo_arr(dir2)) then
         print *,"probhi_arr invalid"
         stop
        endif
        sinprod=sinprod*sin(two*Pi*perturbation_mode* &
          (xsten(0,dir2)-problo_arr(dir2))/ &
          (probhi_arr(dir2)-problo_arr(dir2)))
       enddo

      enddo
      enddo
      enddo ! i,j,k

      return
      end subroutine fort_addnoise

      subroutine fort_init_regions_list( &
       constant_density_all_time, &
       num_threads_in) &
      bind(c,name='fort_init_regions_list')

      use probcommon_module
      use geometry_intersect_module

      IMPLICIT NONE

      integer, INTENT(in) :: num_threads_in
      integer, INTENT(in) :: constant_density_all_time(num_materials)
      integer :: im
      integer :: iregion
      integer :: ithread
      integer :: dir

      if (num_threads_in.eq.geom_nthreads) then
       ! do nothing
      else
       print *,"num_threads_in invalid"
       stop
      endif
      do im=1,num_materials
       if ((constant_density_all_time(im).eq.0).or. &
           (constant_density_all_time(im).eq.1)) then
        ! do nothing
       else
        print *,"constant_density_all_time(im) invalid"
        stop
       endif
      enddo ! im=1..num_materials

      call SUB_INIT_REGIONS_LIST( &
       constant_density_all_time, &
       num_materials, &
       num_threads_in)

      if (number_of_source_regions.eq.0) then
       ! do nothing
      else if (number_of_source_regions.gt.0) then
       do ithread=1,num_threads_in

        do iregion=1,number_of_source_regions
         regions_list(iregion,ithread)%region_material_id= &
           regions_list(iregion,0)%region_material_id

         regions_list(iregion,ithread)%region_dt= &
           regions_list(iregion,0)%region_dt

         regions_list(iregion,ithread)%region_mass_flux= &
           regions_list(iregion,0)%region_mass_flux

         regions_list(iregion,ithread)%region_volume_flux= &
           regions_list(iregion,0)%region_volume_flux

         regions_list(iregion,ithread)%region_temperature_prescribe= &
           regions_list(iregion,0)%region_temperature_prescribe

         do dir=1,SDIM
          regions_list(iregion,ithread)%region_velocity_prescribe(dir)= &
           regions_list(iregion,0)%region_velocity_prescribe(dir)
         enddo

         regions_list(iregion,ithread)%region_energy_flux= &
           regions_list(iregion,0)%region_energy_flux

         regions_list(iregion,ithread)%region_volume_raster= &
           regions_list(iregion,0)%region_volume_raster

         regions_list(iregion,ithread)%region_volume= &
           regions_list(iregion,0)%region_volume
         regions_list(iregion,ithread)%region_mass= &
           regions_list(iregion,0)%region_mass
         regions_list(iregion,ithread)%region_energy= &
           regions_list(iregion,0)%region_energy
         regions_list(iregion,ithread)%region_energy_per_kelvin= &
           regions_list(iregion,0)%region_energy_per_kelvin
         regions_list(iregion,ithread)%region_volume_after= &
           regions_list(iregion,0)%region_volume_after
         regions_list(iregion,ithread)%region_mass_after= &
           regions_list(iregion,0)%region_mass_after
         regions_list(iregion,ithread)%region_energy_after= &
           regions_list(iregion,0)%region_energy_after

        enddo ! iregion=1,number_of_source_regions

       enddo ! ithread=1..num_threads_in
      else
       print *,"number_of_source_regions invalid"
       stop
      endif

      end subroutine fort_init_regions_list

      subroutine fort_delete_regions_list(ioproc) &
      bind(c,name='fort_delete_regions_list')

      use probcommon_module
      use geometry_intersect_module
      IMPLICIT NONE
      integer, INTENT(in) :: ioproc
      integer lower_bound(2)
      integer upper_bound(2)
      integer iregions
      integer dir

      call SUB_DELETE_REGIONS_LIST()

      if (number_of_source_regions.eq.0) then
       ! do nothing
      else if (number_of_source_regions.gt.0) then
       lower_bound=LBOUND(regions_list)
       upper_bound=UBOUND(regions_list)
       if ((lower_bound(1).eq.1).and. &
           (lower_bound(2).eq.0).and. &
           (upper_bound(1).eq.number_of_source_regions).and. &
           (upper_bound(2).eq.number_of_threads_regions)) then

        if (ioproc.eq.1) then
         do iregions=1,number_of_source_regions

          if (regions_list(iregions,0)%region_dt.gt.zero) then
           ! do nothing
          else
           print *,"region_dt must be positive"
           stop
          endif

          print *,"iregions=",iregions
          print *,"regions_list(iregions,0)%region_material_id ", &
            regions_list(iregions,0)%region_material_id
          print *,"regions_list(iregions,0)%region_dt ", &
            regions_list(iregions,0)%region_dt
          print *,"regions_list(iregions,0)%region_mass_after ", &
            regions_list(iregions,0)%region_mass_after
          print *,"regions_list(iregions,0)%region_volume_after ", &
            regions_list(iregions,0)%region_volume_after
          print *,"regions_list(iregions,0)%region_energy_after ", &
            regions_list(iregions,0)%region_energy_after

          print *,"regions_list(iregions,0)%region_mass_flux ", &
            regions_list(iregions,0)%region_mass_flux
          print *,"regions_list(iregions,0)%region_volume_flux ", &
            regions_list(iregions,0)%region_volume_flux

          print *,"regions_list(iregions,0)%region_temperature_prescribe ", &
            regions_list(iregions,0)%region_temperature_prescribe
          do dir=1,SDIM
           print *,"dir,region_velocity_prescribe ", &
            dir,regions_list(iregions,0)%region_velocity_prescribe(dir)
          enddo

          print *,"regions_list(iregions,0)%region_energy_flux ", &
            regions_list(iregions,0)%region_energy_flux

          print *,"regions_list(iregions,0)%region_mass_flux measured: ", &
            (regions_list(iregions,0)%region_mass_after- &
             regions_list(iregions,0)%region_mass)/ &
            regions_list(iregions,0)%region_dt

          print *,"regions_list(iregions,0)%region_volume_flux measured: ", &
            (regions_list(iregions,0)%region_volume_after- &
             regions_list(iregions,0)%region_volume)/ &
            regions_list(iregions,0)%region_dt

          print *,"regions_list(iregions,0)%region_energy_flux measured: ", &
            (regions_list(iregions,0)%region_energy_after- &
             regions_list(iregions,0)%region_energy)/ &
            regions_list(iregions,0)%region_dt

         enddo ! do iregions=1,number_of_source_regions
        else if (ioproc.eq.0) then
         ! do nothing
        else
         print *,"ioproc invalid"
         stop
        endif

        deallocate(regions_list)
       else
        print *,"lower_bound or upper_bound invalid"
        stop
       endif
      else
       print *,"number_of_source_regions invalid"
       stop
      endif

      end subroutine fort_delete_regions_list

      subroutine fort_initvelocity( &
        level,time, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        vel,DIMS(vel), &
        dx,xlo,xhi) &
      bind(c,name='fort_initvelocity')

      use global_distance_module
      use global_utility_module
      use hydrateReactor_module
      use unimaterialChannel_module
      use River
      use USERDEF_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use rigid_FSI_module
      use sinking_particle_module

      IMPLICIT NONE

      real(amrex_real) kterm,velperturb
      real(amrex_real) ktermx,velperturbx

      integer, INTENT(in) :: level
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: DIMDEC(vel)
      real(amrex_real), INTENT(in) :: time, dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM), xhi(SDIM)

      real(amrex_real), INTENT(out), target :: vel(DIMV(vel),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)

!     ::::: local variables
      integer i,j,k
      real(amrex_real) x,y,z
      real(amrex_real) x_vel,y_vel,z_vel,dist
      real(amrex_real) xx_vel,yy_vel,zz_vel
      real(amrex_real) xtemp,ytemp,ztemp
      real(amrex_real) ytop,radcross,rtest
      real(amrex_real) outer_rad,areacross,radshrink
      real(amrex_real) velcell(SDIM)
      real(amrex_real) cenbc(num_materials,SDIM)
      real(amrex_real) vfracbatch(num_materials)
      real(amrex_real) drat
      real(amrex_real) temp,dens,ccnt
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xvec(SDIM)
      integer dir
      real(amrex_real) jumpval,alpha
      real(amrex_real), parameter :: stub_zero=zero

      real(amrex_real), allocatable, dimension(:) :: distbatch
      integer velsolid_flag
 
      vel_ptr=>vel

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif 

      velsolid_flag=0

      call checkbound_array(fablo,fabhi,vel_ptr,1,-1)

      allocate(distbatch(num_materials))

      if (time.eq.zero) then
       ! do nothing
      else
       print *,"time should be zero in initvelocity"
       stop
      endif

      call default_rampvel(time,xx_vel,yy_vel,zz_vel)

      if (adv_vel.ne.zero) then
        print *,"adv_dir,adv_vel ",adv_dir,adv_vel
      endif

      if (SDIM.eq.2) then

! shear (initvelocity)
      if (probtype.eq.1) then
        if (axis_dir.eq.0) then
         print *,"newtonian liquid"
        else if ((axis_dir.gt.0).and.(axis_dir.le.7)) then
         print *,"shear thinning liquid"
        else if (axis_dir.eq.11) then
         print *,"viscoelastic outer fluid"
        else if (axis_dir.eq.12) then
         print *,"viscoelastic drop"
        else if (axis_dir.eq.140) then
         print *,"droplets head on problem vinletgas=initial velocity"
        else if (axis_dir.eq.141) then
         print *,"diff. droplets head on problem vinletgas=initial velocity"
        else if (axis_dir.eq.13) then
         print *,"middle earth flow"
        else if (axis_dir.eq.14) then
         print *,"droplet collision problem vinletgas=initial velocity"
        else if (axis_dir.eq.15) then
         print *,"test problem from Zuzio et al"
        else if (axis_dir.eq.150) then
         print *,"shock drop interaction problem"
        else if (axis_dir.eq.151) then
         print *,"shock column interaction problem"
        else
         print *,"axis_dir invalid probtype=1"
         stop
        endif
! bubble
      else if (probtype.eq.2) then
       if ((axis_dir.lt.0).or.(axis_dir.gt.7)) then
        print *,"axis_dir out of range in initbubble"
        stop
       else if (axis_dir.eq.0) then
        print *,"Newtonian liquid being computed...."
       else
        print *,"non-newtonian generalized cross carreau model liquid"
        print *,"axis_dir=",axis_dir
       endif
! capillary
      else if ((probtype.eq.3).or. &
               (probtype.eq.41)) then
       print *,"2D pipe problem or rayleigh capillary break up test problem"
      else if (probtype.eq.4) then
       print *,"wavenumber is xblob : ",xblob
       print *,"y=radblob*cos(xblob*pi*x), xblob=2 for rt"
       print *,"xx_vel,yy_vel ",xx_vel,yy_vel
! gas burst
      else if (probtype.eq.8) then
       print *,"INITIALIZING RZ (axisym) GAS BURST PROBLEM "
      else if (probtype.eq.14) then
       print *,"this probtype obsolete"
       stop
! jetting 
      else if (probtype.eq.22) then
       print *,"jetting obselete"
! standing wave problem
      else if (probtype.eq.23) then
       print *,"standing wave problem (NOT r-z)"
       print *,"wavelen is xblob : ",xblob
       print *,"perturbation is radblob ",radblob
       print *,"base amplitude is yblob ",yblob
       print *,"y=yblob+radblob*cos(2pi x/xblob)"
       print *,"levelset < 0 in gas, levelset >0 in liquid"
       print *,"vfrac = 0 in gas, vfrac =1 in liquid"
! hanging
      else if (probtype.eq.25) then
       print *,"hanging drop problem or bubble column problem"
       print *,"axis_dir=1..11 if bubble column, axis_dir: ",axis_dir
       print *,"radius of orifice is radblob: ",radblob
       print *,"if axis_dir>0, zblob = height of column=",zblob
       print *,"if axis_dir=0 zblob = radius preejected fluid=",zblob
       print *,"do not set zblob<0"
       print *,"advbot = rate water poured in =",advbot
       print *,"xblob should be 0, xblob=",xblob
       if (xblob.ne.0.0) then
        stop
       endif
       print *,"yblob=y value of inflow, yblob=",yblob
    
       if ((axis_dir.gt.11).or.(axis_dir.lt.0)) then
        print *,"axis_dir out of range for probtype=25"
       endif
       if ((axis_dir.gt.0).and.(zblob.le.zero)) then
        print *,"zblob should be positive for bubble column problem"
        stop
       endif 
! shed
      else if ((probtype.eq.30).or.(probtype.eq.32).or. &
               (probtype.eq.33).or.(probtype.eq.34) ) then
       print *,"probtype=30 means half circle, probtype=32 means full"
       print *,"probtype=33 means drop on a slope"
       print *,"probtype=34 means capillary tube"
       print *,"probtype=",probtype
! meniscus
      else if (probtype.eq.35) then
       print *,"radblob is NID/2 radblob= ",radblob
       print *,"yblob is NPT  yblob= ",yblob
       print *,"in 3d, xblob is domain base size xblob= ",xblob
      else if (probtype.eq.39) then
       print *,"standing wave problem (NOT r-z)"
       print *,"wavelen is xblob : ",xblob
       print *,"perturbation is radblob ",radblob
       print *,"base amplitude is yblob ",yblob
       print *,"y=yblob+radblob*cos(2pi x/xblob)"
       print *,"levelset < 0 in gas, levelset >0 in liquid"
       print *,"vfrac = 0 in gas, vfrac =1 in liquid"
      else if (probtype.eq.40) then
       if (adv_dir .eq. 1) then
         print *,"translation in x-direction with adv_vel=",adv_vel
       else if (adv_dir .eq. 2) then
         print *,"translation in y-direction with adv_vel=",adv_vel
       else if (adv_dir.eq.3) then
         print *,"translation in x and y-direction with adv_vel=",adv_vel
       else if (adv_dir.eq.4) then
         print *,"solid body rotation with adv_vel=",adv_vel
       else if (adv_dir.eq.5) then
         print *,"stretching with adv_vel=",adv_vel
       else
         write(6,*) "error: initvortpatch: adv_dir = ",adv_dir
         stop
       endif
! overturn
      else if (probtype.eq.45) then
       ytop=0.5
       print *,"using ytop=.5"
      endif


      else if (SDIM.eq.3) then

! shear (initvelocity)
      if (probtype.eq.1) then
        if (axis_dir.eq.0) then
         print *,"newtonian liquid"
        else if ((axis_dir.gt.0).and.(axis_dir.le.7)) then
         print *,"shear thinning liquid"
        else if (axis_dir.eq.11) then
         print *,"LS<0 inside of drop (gas) and LS>0 outside drop (liquid)" 
        else if (axis_dir.eq.12) then
         print *,"viscoelastic drop (LS>0 inside, LS<0 outside)"
        else if (axis_dir.eq.13) then
         print *,"middle earth flow"
        else if (axis_dir.eq.14) then
         print *,"droplet collision problem vinletgas=initial velocity"
        else if (axis_dir.eq.15) then
         print *,"Zuzio test problem"
        else if (axis_dir.eq.150) then
         print *,"shock drop interaction problem"
        else if (axis_dir.eq.151) then
         print *,"shock column interaction problem"
        else
         print *,"axis_dir invalid shear probtype axis_dir ", &
          probtype,axis_dir
         stop
        endif
! bubble
      else if (probtype.eq.2) then
       if ((axis_dir.lt.0).or.(axis_dir.gt.7)) then
        print *,"axis_dir out of range in initbubble"
        stop
       else if (axis_dir.eq.0) then
        print *,"Newtonian liquid being computed...."
       else
        print *,"non-newtonian generalized cross carreau model liquid"
        print *,"axis_dir=",axis_dir
       endif
! pipe
      else if (probtype.eq.41) then

       if (axis_dir.eq.5) then
        ! do nothing
       else
        print *,"pipe problem setup should be modified in 3d"
        stop
       endif

      else if (probtype.eq.4) then
       print *,"wavenumber is xblob : ",xblob
       print *,"y=radblob*cos(xblob*pi*x), xblob=2 for rt"
       print *,"xx_vel,yy_vel ",xx_vel,yy_vel
! splash
      else if (probtype.eq.7) then
       print *,"this probtype obsolete"
       stop
! gas burst
      else if (probtype.eq.18) then
       print *,"not a 3d problem"
       stop
! jetting 
      else if (probtype.eq.22) then
       print *,"jetting obselete"
! standing wave problem
      else if (probtype.eq.23) then
       print *,"standing wave problem (NOT r-z)"
       print *,"wavelen is xblob : ",xblob
       print *,"perturbation is radblob ",radblob
       print *,"base amplitude is yblob ",yblob
       print *,"y=yblob+radblob*cos(2pi x/xblob)"
       print *,"levelset < 0 in gas, levelset >0 in liquid"
       print *,"vfrac = 0 in gas, vfrac =1 in liquid"
! hanging
      else if (probtype.eq.25) then
       print *,"hanging drop problem or bubble column problem"
       print *,"axis_dir=1..11 if bubble column, axis_dir: ",axis_dir
       print *,"radius of orifice is radblob: ",radblob
       print *,"is axis_dir>0, zblob = height of column=",zblob
       print *,"otherwise zblob = radius preejected fluid=",zblob
       print *,"do not set zblob<0"
       print *,"advbot = rate water poured in =",advbot
       print *,"xblob should be 0, xblob=",xblob
       if (xblob.ne.0.0) then
        stop
       endif
       print *,"yblob=y value of inflow, yblob=",yblob
    
       if ((axis_dir.gt.11).or.(axis_dir.lt.0)) then
        print *,"axis_dir out of range in inithanging"
       endif
       if ((axis_dir.gt.0).and.(zblob.lt.zero)) then
        print *,"zblob should be non-negative for bubble column problem"
        stop
       endif 
! shed
      else if ((probtype.eq.30).or.(probtype.eq.32).or. &
               (probtype.eq.33).or.(probtype.eq.34) ) then
       print *,"probtype=30 means half circle, probtype=32 means full"
       print *,"probtype=33 means drop on a slope"
       print *,"probtype=34 means capillary tube"
       print *,"probtype=",probtype
! meniscus
      else if (probtype.eq.35) then
       print *,"radblob is NID/2 radblob= ",radblob
       print *,"yblob is NPT  yblob= ",yblob
       print *,"in 3d, xblob is domain base size xblob= ",xblob
      else if (probtype.eq.39) then
       print *,"standing wave problem (NOT r-z)"
       print *,"wavelen is xblob : ",xblob
       print *,"perturbation is radblob ",radblob
       print *,"base amplitude is yblob ",yblob
       print *,"y=yblob+radblob*cos(2pi x/xblob)"
       print *,"levelset < 0 in gas, levelset >0 in liquid"
       print *,"vfrac = 0 in gas, vfrac =1 in liquid"
! overturn
      else if (probtype.eq.45) then
       ytop=0.5
       print *,"using ytop=.5"
      endif

      else
       print *,"dimension bust"
       stop
      endif

      if (1.eq.0) then
       print *,"xx_vel= ",xx_vel 
       print *,"yy_vel= ",yy_vel 
       print *,"zz_vel= ",zz_vel 
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

        x_vel=xx_vel
        y_vel=yy_vel
        z_vel=zz_vel

        call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

        x=xsten(0,1)
        y=xsten(0,2)
        z=xsten(0,SDIM)
        do dir=1,SDIM
         xvec(dir)=xsten(0,dir)
        enddo

        if (is_in_probtype_list().eq.1) then

         call SUB_LS(xvec,time,distbatch,num_materials)
          ! pass dx
         call SUB_VEL(xvec,time,distbatch,velcell, &
          velsolid_flag,dx,num_materials)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

        else if (probtype.eq.401) then
         call HELIX_LS(xvec,time,distbatch)
         call HELIX_VEL(xvec,time,distbatch,velcell,velsolid_flag)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

        else if (probtype.eq.402) then
         call TSPRAY_LS(xvec,time,distbatch)
         call TSPRAY_VEL(xvec,time,distbatch,velcell,velsolid_flag)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

        else if (probtype.eq.412) then ! step
         call CAV2Dstep_LS(xvec,time,distbatch)
         call CAV2Dstep_VEL(xvec,time,distbatch,velcell,velsolid_flag)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

        else if (probtype.eq.413) then ! ZEYU droplet impact
         call ZEYU_droplet_impact_LS(xvec,time,distbatch)
          ! pass dx
         call ZEYU_droplet_impact_LS_VEL(xvec,time,distbatch,velcell, &
          velsolid_flag,dx)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

        else if (probtype.eq.533) then
         call rigid_FSI_LS(xvec,time,distbatch)
         call rigid_FSI_VEL(xvec,time,distbatch,velcell,velsolid_flag)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)
        else if (probtype.eq.534) then
         call sinking_FSI_LS(xvec,time,distbatch)
         call sinking_FSI_VEL(xvec,time,distbatch,velcell,velsolid_flag)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

        else if (probtype.eq.311) then ! user defined example
         call USERDEF_LS(xvec,time,distbatch)
         call USERDEF_VEL(xvec,time,distbatch,velcell,velsolid_flag)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

         ! HYDRATE  (in initvelocity)
        else if (probtype.eq.199) then
         call materialdist_batch(xsten,nhalf,dx,bfact,distbatch,time)
         if (distbatch(1).ge.zero) then
          call INIT_STATE_WATER(x,y,z,time,velcell,temp,dens,ccnt)
         else if (distbatch(2).ge.zero) then
          call INIT_STATE_GAS(x,y,z,time,velcell,temp,dens,ccnt)
         else
          call INIT_STATE_HYDRATE(x,y,z,time,velcell,temp,dens,ccnt)
         endif
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

         ! in: fort_initvelocity
        else if (probtype.eq.220) then
         call UNIMAT_INIT_VEL(x,y,z,velcell)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

        else if ((probtype.eq.299).or. &
                 (probtype.eq.301)) then ! melting, initial velocity

         x_vel=zero
         y_vel=zero
         z_vel=zero

        else if (probtype.eq.209) then  ! river

         call RiverVelocity(x,y,z,velcell,axis_dir,probloz,probhiz)
         x_vel=velcell(1)
         y_vel=velcell(2)
         z_vel=velcell(SDIM)

          ! Zuzio, initvelocity
        else if ((probtype.eq.1).and.(axis_dir.eq.15)) then 
         x_vel=zero
         y_vel=zero
         z_vel=zero
         ! Marioff injector
        else if (probtype.eq.537) then
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
          z_vel=velcell(SDIM)

          ! in "initvelocity":
          ! melting ice block on substrate.
        else if (probtype.eq.59) then
          x_vel=zero
          y_vel=zero
          z_vel=zero

        else if (probtype.eq.201) then ! stratified bubble (initvelocity)

         if (advbot.eq.zero) then
          ! do nothing
         else
          x_vel=zero
          y_vel=zero
          z_vel=zero
          call get_initial_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc)
          if (vfracbatch(2).gt.zero) then
           if (SDIM.eq.2) then
            y_vel=-abs(advbot)
           else if (SDIM.eq.3) then
            z_vel=-abs(advbot)
           else
            print *,"dimension bust"
            stop
           endif
          endif
         endif ! advbot <> 0

         ! shock tube problems (initvelocity)
        else if ((probtype.eq.92).or. & !contact is captured
                 (probtype.eq.93)) then !contact is tracked
         x_vel=zero
         y_vel=zero
         z_vel=zero
         if (axis_dir.eq.0) then ! Sod shock tube
          if (x.le.half) then
           x_vel=zero
          else
           x_vel=zero
          endif
         else if (axis_dir.eq.1) then ! strong shock tube
          if (x.le.half) then
           x_vel=zero
          else
           x_vel=zero
          endif
         else if (axis_dir.eq.2) then ! shock-turbulence
          if (x.le.one) then
           x_vel=2.629369
          else
           x_vel=zero
          endif
         else if (axis_dir.eq.3) then ! mach>4
          if (x.le.one) then
           x_vel=5.0
          else
           x_vel=zero
          endif
         else if (axis_dir.eq.4)  then ! smooth problem
          x_vel=zero
          y_vel=zero
          z_vel=zero
         else
          print *,"axis_dir invalid probtype=92,93: ",axis_dir
          stop
         endif

        else if (SDIM.eq.2) then


         if (probtype.eq.801.and.axis_dir.eq.3) then ! convective evaporation
          if(sqrt( (x-xblob)**2+(y-yblob)**2 ).lt.radblob) then
            x_vel=zero
          endif

           ! dissolution initial velocity
         else if (probtype.eq.802) then
          x_vel=zero
          y_vel=zero
          call vapordist(xsten,nhalf,dx,bfact,dist)
          if (dist.gt.radblob) then
           dist=radblob
          endif
          if (dist.ge.half*dx(SDIM)) then
           kterm=two*Pi*yblob2*y/(two*radblob)
           velperturb=one+radblob2*cos(kterm)
           ktermx=two*Pi*yblob2*x/(two*radblob)
           velperturbx=one+radblob2*cos(ktermx)
           if (1.eq.0) then
            x_vel=1.5*adv_vel*velperturb*velperturbx* &
               (one-(dist/radblob)**2)
           else
            x_vel=adv_vel*velperturb*velperturbx  ! plug flow
           endif
           y_vel=adv_vel*radblob2*cos(kterm)*cos(ktermx)
          endif
         else if (probtype.eq.602) then  ! Rayleigh Taylor and checkerboard
          x_vel=zero
          y_vel=zero
         else if ((probtype.eq.1).and.(axis_dir.lt.150)) then
          if (axis_dir.eq.11) then
           x_vel=vinletgas*(y/yblob-one)
          else if (axis_dir.eq.14) then
           dist=radblob-sqrt((x-xblob)**2+(y-yblob)**2)
           if (dist.ge.zero) then
             x_vel=zero
             y_vel=vinletgas
           endif
          else if ((axis_dir.eq.140).or.(axis_dir.eq.141)) then
           dist=max(-sqrt( (x-xblob)**2 + (y-yblob)**2 )+radblob, &
                    -sqrt( (x-xblob)**2 + (y-yblob2)**2 )+radblob)

           if (dist.ge.zero) then
             x_vel=zero
             if( y > half*(yblob+yblob2))then
              y_vel=-abs(vinletgas)
             else
              y_vel=abs(vinletgas)
             endif
           endif
          else if ((axis_dir.eq.12).and.(adv_dir.eq.1).and. &
                   (vinletgas.eq.zero)) then
           x_vel=zero  ! no velocity in the droplet at t=0
          endif
         else if (probtype.eq.531) then ! falling sphere - INIT_VELOCITY
          x_vel=zero
          y_vel=zero
! Reiber problem
         else if (probtype.eq.540) then
          call get_Rieber_velocity(xsten,nhalf,bfact,dx,velcell)
          y_vel=velcell(SDIM)
         else if (probtype.eq.17) then  ! drop collide of diesel and water
! vb-vt=1
! vb db + vt dt=0
! -vt dt/db - vt = 1
! vt=-1/(dt/db + 1)
! vb=1+vt=dt/db / (1+dt/db)
          call get_initial_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc)
          drat=fort_denconst(3)/fort_denconst(1)  ! dt/db
          if (vfracbatch(1).gt.zero) then  ! diesel on bottom
           y_vel=drat/(one+drat)
          else if (vfracbatch(3).gt.zero) then  ! water on top
           y_vel=-one/(one+drat)
          else
           y_vel=zero
          endif  
         else if (probtype.eq.18) then  ! drop collide same material
          call get_initial_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc)
          if (vfracbatch(1).gt.zero) then
           if (y.gt.zero) then
            y_vel=-half
           else
            y_vel=half
           endif
          else
           y_vel=zero
          endif
         else if (probtype.eq.51) then ! oscillating column
          call get_initial_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc)
          if (vfracbatch(1).gt.zero) then
           x_vel=adv_vel
          else
           x_vel=zero
          endif
         else if (probtype.eq.102) then ! nozzle
          if (yblob3.le.zero) then
           print *,"yblob3 invalid"
           stop
          endif
          x_vel=zero
          y_vel=zero 
           ! liquid nozzle: 0<y<yblob+yblob2
          if (y.le.yblob+yblob2) then
           outer_rad=radblob3-y*(radblob3-radblob4)/yblob3 
           areacross=Pi*(outer_rad**2-radblob5**2)
           if ((x.ge.radblob5).and.(x.le.outer_rad)) then
            y_vel=advbot*1000.0d0/areacross
           endif 
           ! gas nozzle: 0<y<yblob3
          else if (y.le.yblob3) then
           outer_rad=radblob3-y*(radblob3-radblob4)/yblob3
           areacross=Pi*outer_rad**2
           if (x.le.outer_rad) then
            y_vel=advbot*1000.0d0/areacross
           endif
          else  ! expansion region
           radshrink=radblob7**2-radblob5**2
           if (radshrink.le.zero) then
            print *,"radshrink invalid"
            stop
           endif
           radshrink=sqrt(radshrink)
           outer_rad=radblob4+ &
             (y-yblob3)*(radshrink-radblob4)/(probhiy-yblob3)
           areacross=Pi*outer_rad**2
           if (x.le.outer_rad) then
            y_vel=advbot*1000.0d0/areacross
           endif
          endif 
! microfluidics problem initial velocity at t=0
         else if (probtype.eq.5700) then
          x_vel=zero
          y_vel=zero
!         z_vel=zero

           ! in: fort_initvelocity (2D)
         else if ((probtype.eq.3).or. &
                  (probtype.eq.41)) then

          if ((axis_dir.eq.0).or. &
              (axis_dir.eq.1).or. &
              (axis_dir.eq.2).or. &
              (axis_dir.eq.3)) then
           call get_pipe_velocity(xsten,nhalf,dx,bfact,velcell,stub_zero)
           x_vel=velcell(1)
           y_vel=velcell(2)

           ! in: fort_initvelocity (2D)
          else if (axis_dir.eq.5) then
           call get_pipe_velocity(xsten,nhalf,dx,bfact,velcell,stub_zero)
           x_vel=velcell(1)
           y_vel=velcell(2)
!          z_vel=zero
          else if (axis_dir.eq.4) then

           if (1.eq.0) then
            call get_pipe_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc) 
            if (vfracbatch(1).gt.zero) then
             dist=half
            else
             dist=-half
            endif
           else
            call inletpipedist(x,y,z,distbatch)   
            dist=distbatch(1)
           endif
  
           call get_pipe_velocity(xsten,nhalf,dx,bfact,velcell,stub_zero)!time=0
           y_vel=velcell(2)
           x_vel=zero
          else
           print *,"axis_dir invalid initvel axis_dir=",axis_dir
           stop
          endif

! rotate
         else if (probtype.eq.5) then
          print *,"this problem obsolete"
          stop
! oilexpel
         else if (probtype.eq.16) then
          if ((y.gt.yblob).and.(x.le.xblob+radblob)) then
           y_vel = -abs(advbot)
          endif
         else if (probtype.eq.23) then
          call vapordist(xsten,nhalf,dx,bfact,dist)
          if (dist.ge.zero) then
           x_vel=zero
           y_vel=zero
          endif
! validate
         else if (probtype.eq.24) then
          x_vel=-sin(Pi*x)*sin(Pi*x)*sin(two*Pi*y)
          y_vel=sin(Pi*y)*sin(Pi*y)*sin(two*Pi*x)
         else if (probtype.eq.25) then  ! in initvelocity
          if (axis_dir.eq.0) then
           if ((y.gt.yblob).and.(x.le.xblob+radblob)) then
            y_vel = -abs(advbot)
           endif
          endif
! swirl 2D, in: fort_initvelocity
         else if (probtype.eq.26) then

          if ((axis_dir.eq.0).or. & !swirl
              (axis_dir.eq.1)) then

           if (y.le.half) then
            x_vel=tanh( (y-one/four)*30.0d0 )
           else
            x_vel=tanh( (three/four-y)*30.0d0 )
           endif
           y_vel=0.05d0*sin(two*Pi*x)

          else if ((axis_dir.eq.2).or. & !vortex confinement
                   (axis_dir.eq.3)) then

           dist=sqrt((x-xblob)**2+(y-yblob)**2)-radblob
           jumpval=tanh(30.0d0*dist)
           jumpval=(jumpval+one)/two
           alpha=(one-jumpval)*vinletgas
           x_vel=alpha*(y-yblob)
           y_vel=-alpha*(x-xblob)
           if ((adv_dir.eq.1).or.(adv_dir.eq.3)) then
            x_vel=x_vel+adv_vel
           endif
           if ((adv_dir.eq.2).or.(adv_dir.eq.3)) then
            y_vel=y_vel+adv_vel
           endif
          else if (axis_dir.eq.10) then !BCG homogeneous bc
           x_vel=-(sin(Pi*x)**2)*sin(two*Pi*y)
           y_vel=sin(two*Pi*x)*(sin(Pi*y)**2)
          else if (axis_dir.eq.11) then  ! 2D BCG periodic
           x_vel=-sin(two*Pi*x)*cos(two*Pi*y)
           y_vel=cos(two*Pi*x)*sin(two*Pi*y)
           if ((adv_dir.eq.1).or.(adv_dir.eq.3)) then
            x_vel=x_vel+adv_vel
           endif
           if ((adv_dir.eq.2).or.(adv_dir.eq.3)) then
            y_vel=y_vel+adv_vel
           endif
          else if (axis_dir.eq.12) then  ! buoyancy
           x_vel=zero
           y_vel=zero
          else
           print *,"axis_dir invalid"
           stop
          endif
             
         else if (probtype.eq.202) then  ! liquid lens
          ! do nothing (adv_vel used above if prescribed)
         else if (probtype.eq.36) then ! bubble 2D
          if ((axis_dir.eq.2).or.(axis_dir.eq.4)) then
           x_vel=zero
           y_vel=zero
          endif
          if ((xblob10.gt.zero).and. &
              (yblob10.ne.zero)) then
           y_vel=x*yblob10/xblob10
          endif
         else if (probtype.eq.37) then
          ! do nothing
         else if (probtype.eq.11) then
          print *,"cavitation with outflow top is deleted"
          stop

          ! in: fort_initvelocity (2D section)
         else if (probtype.eq.42) then
          ! do nothing - bubble jetting 2D
         else if (probtype.eq.46) then
          if ((axis_dir.ge.0).and.(axis_dir.lt.10)) then
           ! do nothing cavitation 2D, jwl
          else if (axis_dir.eq.10) then
           if (1.eq.0) then
            call get_initial_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc)
            if (vfracbatch(num_materials).gt.zero) then ! sphere
             y_vel=advbot
            else
             y_vel=zero
            endif
           endif
          else if (axis_dir.eq.20) then
           x_vel=adv_vel
           y_vel=zero
          else
           print *,"axis_dir invalid"
           stop
          endif
! kh
         else if (probtype.eq.38) then
          call vapordist(xsten,nhalf,dx,bfact,dist)
          if (dist.lt.zero) then
           y_vel=-x_vel*radblob*two*Pi*sin(two*Pi*x/xblob)/xblob
          else
           x_vel=zero
           y_vel=zero
          endif
! vstanding
         else if (probtype.eq.39) then
          if (axis_dir.eq.0) then
           x_vel=zero
           y_vel=zero
          else
           print *,"this problem obsolete"
           stop
          endif
! vortpatch
         else if (probtype.eq.40) then
          x_vel=zero
          y_vel=zero
! overturn
         else if (probtype.eq.45) then
          print *,"this probtype obsolete"
          stop
! paddle
         else if (probtype.eq.50) then
          if (y.lt.zblob) then
           x_vel=zero
          endif
         else if (probtype.eq.58) then
          x_vel=zero
! jetbend
         else if (probtype.eq.53) then
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
! 2d diesel injector w/needle
         else if ((probtype.eq.538).or.(probtype.eq.541)) then
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
! supersonic nozzle: fort_initvelocity
         else if (probtype.eq.539) then
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
         else if (probtype.eq.532) then ! imp jets from sides (initvelocity)
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)

! 3D jet coaxial
         else if (probtype.eq.72) then
          call vapordist(xsten,nhalf,dx,bfact,dist)
          if (dist.ge.zero) then
           x_vel=advbot
          else
           x_vel=adv_vel
          endif
! milkdrop
         else if (probtype.eq.61) then
          if (sqrt( (x-xblob)**2+(y-yblob)**2 ).le.radblob) then
           if (axis_dir.eq.1) then
            y_vel=-one
           endif
          endif
! nozzle
         else if ((probtype.eq.63).or.(probtype.eq.64)) then
          call nozzlerad(z,radcross,stub_zero)
          if (x.gt.radcross) then
           y_vel=zero
          else
           y_vel=y_vel*(xblob10**2/radcross**2)
          endif
! pulse
         else if (probtype.eq.66) then
          xtemp=sqrt(three*radblob/(four*zblob*zblob*zblob))
          x_vel=sqrt(9.8*zblob)*(radblob/zblob)/(cosh(xtemp*x)**2)
          y_vel=sqrt(three*9.8*zblob)*((radblob/zblob)**(1.5))* &
            (y/zblob)*tanh(xtemp*x)/(cosh(xtemp*x)**2)
         else if (probtype.eq.110) then
          call get_bump_velocity(xsten,nhalf,dx,bfact,x_vel,time)
          y_vel=zero
         else if (probtype.eq.701) then
          ! flapping wing, init_velocity, do nothing: x_vel=adv_vel
         endif

        else if (SDIM.eq.3) then

         if ((probtype.eq.1).and.(axis_dir.lt.150)) then
          if ((axis_dir.eq.11).or.(axis_dir.eq.12)) then
           if (zblob.gt.zero) then
            x_vel=vinletgas*(z/zblob-one)
           else if (zblob.eq.zero) then
            if (probhiz.le.zero) then
             print *,"probhiz invalid"
             stop
            endif
            x_vel=vinletgas*z/probhiz
           else
            print *,"parameters invalid for shear problem"
            stop
           endif

           if (xblob10.eq.one) then
            if (zblob.gt.zero) then
             x_vel=radblob10+(vinletgas-radblob10)*z/(two*zblob)
            else
             if (radblob10.ne.zero) then
              print *,"zero velocity at axis of symmetry"
             endif
             if (zblob10.eq.zero) then
              print *,"zblob10 should be domain height"
              stop
             endif
             x_vel=vinletgas*z/zblob10
            endif
           endif
          else if (axis_dir.eq.14) then
           dist=radblob-sqrt((x-xblob)**2+(y-yblob)**2)
           if (dist.ge.zero) then
            x_vel=zero
            y_vel=vinletgas
           endif
          endif

! pipe setup at t=0
! in: fort_initvelocity, 3D
         else if (probtype.eq.41) then
          if (axis_dir.eq.5) then
           call get_pipe_velocity(xsten,nhalf,dx,bfact,velcell,stub_zero)
           x_vel=velcell(1)
           y_vel=velcell(2)
           z_vel=velcell(SDIM)
          else
           print *,"this problem not ready for 3d yet"
           stop
          endif

! wave
         else if (probtype.eq.13) then
          print *,"option obsolete"
          stop
         else if (probtype.eq.14) then
          print *,"option obsolete"
          stop
         else if (probtype.eq.16) then
          if ((y.gt.yblob).and.(x.le.xblob+radblob)) then
           y_vel = -abs(advbot)
          endif
         else if (probtype.eq.23) then
          print *,"this option called in error"
          stop
! validate
         else if (probtype.eq.24) then
          x_vel=-sin(Pi*x)*sin(Pi*x)*sin(two*Pi*y)
          y_vel=sin(Pi*y)*sin(Pi*y)*sin(two*Pi*x)
         else if (probtype.eq.25) then
          if (axis_dir.eq.0) then
           if ((y.gt.yblob).and.(x.le.xblob+radblob)) then
            y_vel = -abs(advbot)
           endif
          endif
! swirl 3D
         else if (probtype.eq.26) then 

          if ((axis_dir.eq.0).or. & !swirl
              (axis_dir.eq.1)) then
           ! x-y
           if (adv_dir.eq.3) then
            if (y.le.half) then
             x_vel=tanh( (y-one/four)*30.0 )
            else
             x_vel=tanh( (three/four-y)*30.0 )
            endif
            y_vel=0.05*sin(two*Pi*x)
            z_vel=zero
           ! x-z
           else if (adv_dir.eq.2) then
            if (z.le.half) then
             x_vel=tanh( (z-one/four)*30.0 )
            else
             x_vel=tanh( (three/four-z)*30.0 )
            endif
            z_vel=0.05*sin(two*Pi*x)
            y_vel=zero
           ! y-z
           else if (adv_dir.eq.1) then
            if (z.le.half) then
             y_vel=tanh( (z-one/four)*30.0 )
            else
             y_vel=tanh( (three/four-z)*30.0 )
            endif
            z_vel=0.05*sin(two*Pi*y)
            x_vel=zero
           else
            print *,"adv_dir invalid probtype==26 (13)"
            stop
           endif

          else if ((axis_dir.eq.2).or. & !vortex confinement 3D
                   (axis_dir.eq.3)) then

           dist=sqrt((x-xblob)**2+(y-yblob)**2+(z-zblob)**2)-radblob
           jumpval=tanh(30.0*dist)
           jumpval=(jumpval+one)/two
           alpha=(one-jumpval)*vinletgas
           x_vel=alpha*(y-yblob)
           y_vel=-alpha*(x-xblob)
           if ((adv_dir.eq.1).or.(adv_dir.eq.4)) then
            x_vel=x_vel+adv_vel
           endif
           if ((adv_dir.eq.2).or.(adv_dir.eq.4)) then
            y_vel=y_vel+adv_vel
           endif
           if ((adv_dir.eq.3).or.(adv_dir.eq.4)) then
            z_vel=z_vel+adv_vel
           endif

          else if (axis_dir.eq.11) then ! 3D BCG periodic

           if ((probhix.eq.one).and. &
               (probhiy.eq.one).and. &
               (probhiz.eq.half)) then
            x_vel=-sin(two*Pi*x)*cos(two*Pi*y)
            y_vel=cos(two*Pi*x)*sin(two*Pi*y)
            z_vel=zero
           else if ((probhix.eq.one).and. &
                    (probhiy.eq.half).and. &
                    (probhiz.eq.one)) then
            x_vel=-sin(two*Pi*x)*cos(two*Pi*z)
            z_vel=cos(two*Pi*x)*sin(two*Pi*z)
            y_vel=zero
           else if ((probhix.eq.half).and. &
                    (probhiy.eq.one).and. &
                    (probhiz.eq.one)) then
            y_vel=-sin(two*Pi*y)*cos(two*Pi*z)
            z_vel=cos(two*Pi*y)*sin(two*Pi*z)
            x_vel=zero
           else
            print *,"probhi x,y, or z invalid"
            stop
           endif

           if ((adv_dir.eq.1).or.(adv_dir.eq.4).or. &
               (adv_dir.eq.5).or.(adv_dir.eq.7)) then
            x_vel=x_vel+adv_vel
           endif
           if ((adv_dir.eq.2).or.(adv_dir.eq.4).or. &
               (adv_dir.eq.6).or.(adv_dir.eq.7)) then
            y_vel=y_vel+adv_vel
           endif
           if ((adv_dir.eq.3).or.(adv_dir.eq.5).or. &
               (adv_dir.eq.6).or.(adv_dir.eq.7)) then
            z_vel=z_vel+adv_vel
           endif

          else if (axis_dir.eq.12) then ! buoyancy
           x_vel=zero
           y_vel=zero
           z_vel=zero
          else
           print *,"axis_dir invalid"
           stop
          endif

! vbubble - this routine is initvelocity
         else if (probtype.eq.36) then ! bubble 3D
          if ((axis_dir.eq.2).or.(axis_dir.eq.4)) then
           x_vel=zero
           y_vel=zero
           z_vel=zero
          endif

          if ((xblob10.gt.zero).and. &
              ((yblob9.ne.zero).or.(yblob10.ne.zero))) then
           if (probhix-problox.le.zero) then
            print *,"probhix or problox invalid"
            stop
           endif
           z_vel=yblob9+(x-problox)*(yblob10-yblob9)/(probhix-problox)
          endif
         else if (probtype.eq.37) then
          ! do nothing

          ! in: fort_initvelocity (3D section)
         else if (probtype.eq.42) then
          ! do nothing: bubble jetting 3D
         else if (probtype.eq.46) then
          ! do nothing: cavitation 3D
! kh
         else if (probtype.eq.38) then
          call vapordist(xsten,nhalf,dx,bfact,dist)
          if (dist.lt.zero) then
           y_vel=-x_vel*radblob*two*Pi*sin(two*Pi*x/xblob)/xblob
          else
           x_vel=zero
           y_vel=zero
          endif
! vstanding
         else if (probtype.eq.39) then
          print *,"this problem obsolete"
          stop
! vortpatch
         else if (probtype.eq.40) then
          x_vel=zero
          y_vel=zero
! overturn
         else if (probtype.eq.45) then
          print *,"this problem obsolete" 
          stop
! paddle
         else if (probtype.eq.50) then
          if ((z.lt.zblob-radblob).or.(z.gt.zblob+radblob).or. &
              (y.lt.yblob-radblob).or.(y.gt.yblob+radblob)) then
           x_vel=zero
          endif
! bering
         else if (probtype.eq.51) then
          x_vel=zero
         else if (probtype.eq.58) then
          x_vel=zero
         else if (probtype.eq.5501) then ! drop hitting rough surface
          if (axis_dir.ne.0) then
           print *,"axis_dir invalid"
           stop
          endif
          x_vel=zero
          y_vel=zero
          z_vel=zero
           ! initialize velocity of droplet only
          call get_initial_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc)
          if (vfracbatch(1).gt.zero) then
           z_vel=-abs(advbot)
          endif
! microfluidics problem initial velocity at t=0
         else if (probtype.eq.5700) then
          x_vel=zero
          y_vel=zero
          z_vel=zero
! airblast with coaxial air flow  at t=0
         else if (probtype.eq.529) then
          x_vel=zero
          y_vel=zero
          z_vel=zero
! jetbend at t=0
         else if (probtype.eq.53) then
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
          z_vel=velcell(SDIM)
! impinging jets - AIAA 2008-4847
         else if (probtype.eq.530) then
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
          z_vel=velcell(SDIM)
         else if (probtype.eq.532) then ! imp jets from sides (initvelocity)
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
          z_vel=velcell(SDIM)
         else if (probtype.eq.536) then
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
          z_vel=velcell(SDIM)

          ! initvelocity
          ! 3D diesel injector w/needle
         else if ((probtype.eq.538).or.(probtype.eq.541)) then
          call get_jetbend_velocity(xsten,nhalf,dx,bfact,velcell)
          x_vel=velcell(1)
          y_vel=velcell(2)
          z_vel=velcell(SDIM)

! Reiber problem
         else if (probtype.eq.540) then
          call get_Rieber_velocity(xsten,nhalf,bfact,dx,velcell)
          z_vel=velcell(SDIM)

! ysl 05/12/14
         else if (probtype.eq.17) then  ! drop collide of diesel and water
          call get_initial_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc)
          drat=fort_denconst(3)/fort_denconst(1)  ! dt/db
          if (vfracbatch(1).gt.zero) then  ! diesel on bottom
           y_vel=drat/(one+drat)
          else if (vfracbatch(3).gt.zero) then  ! water on top
           y_vel=-one/(one+drat)
          else
           y_vel=zero
          endif
         else if (probtype.eq.18) then  ! drop collide same material
          call get_initial_vfrac(xsten,nhalf,dx,bfact,vfracbatch,cenbc)
          if (vfracbatch(1).gt.zero) then
           if (y.gt.zero) then
            y_vel=-half
           else
            y_vel=half
           endif
          else
           y_vel=zero
          endif

! milkdrop (initvelocity)
         else if ((probtype.eq.61).or.(probtype.eq.64)) then
          dist=sqrt( (x-xblob)**2+(y-yblob)**2+(z-zblob)**2 )-radblob
          if (dist.le.half*dx(SDIM)) then
           if (axis_dir.eq.1) then
            z_vel=-one
           else if (axis_dir.eq.2) then
            z_vel=vinletgas
           else
            print *,"axis_dir invalid"
            stop
           endif
          endif
! nozzle
         else if (probtype.eq.63) then
          call nozzlerad(z,radcross,stub_zero)
          rtest=sqrt(x**2+y**2)
          if (rtest.gt.radcross) then
           z_vel=zero
          else
           z_vel = z_vel*(xblob10**2/radcross**2)
          endif
! coffee
         else if (probtype.eq.65) then
! center of whirlpool is (xtemp,ytemp) = (5,5); strength = radblob5
          xtemp = (xlo(1)+xhi(1))/2.0
          ytemp = (xlo(2)+xhi(2))/2.0
          ztemp = (ytemp-y)*(ytemp-y)+(x-xtemp)*(x-xtemp)+one
          x_vel = radblob5*(ytemp-y)/ztemp
          y_vel = radblob5*(x-xtemp)/ztemp
          if (z.gt.zblob3) then
           x_vel = zero
           y_vel = zero
          endif

! pulse
         else if (probtype.eq.66) then
          xtemp=sqrt(three*radblob/(four*zblob*zblob*zblob))
          x_vel=sqrt(9.8*zblob)*(radblob/zblob)/(cosh(xtemp*x)**2)
          y_vel=zero
          z_vel=sqrt(three*9.8*zblob)*((radblob/zblob)**(1.5))* &
            (z/zblob)*tanh(xtemp*x)/(cosh(xtemp*x)**2)
! gear
         else if (probtype.eq.563) then
          if (levelrz.eq.COORDSYS_CARTESIAN) then
           if (radblob.gt.zero) then
            call cylinderdist(y,z,x,yblob,zblob,radblob,xblob-radblob, &
              xblob+radblob,dist)
            dist=-dist
            if (dist.ge.zero) then
             x_vel=advbot
            endif
           endif
          else 
           print *,"levelrz invalid probtype = 563"
           stop
          endif
         else if (probtype.eq.701) then
          ! flapping wing, init_velocity, do nothing: x_vel=adv_vel
         endif

        else
         print *,"dimension bust"
         stop
        endif

        vel(D_DECL(i,j,k),1) = x_vel
        vel(D_DECL(i,j,k),2) = y_vel
        if (SDIM.eq.3) then
         vel(D_DECL(i,j,k),SDIM) = z_vel
        endif

      enddo
      enddo
      enddo

      deallocate(distbatch)

      return
      end subroutine fort_initvelocity

      subroutine fort_forcevelocity( &
        problo,probhi, &
        vel,DIMS(vel), &
        velmac,DIMS(velmac), &
        dir, &
        xlo,dx, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        time, &
        presbc_array, &
        outflow_velocity_buffer_size) & !(1,1),(2,1),(3,1),(1,2),(2,2),(3,2)
      bind(c,name='fort_forcevelocity')

      use global_utility_module

      IMPLICIT NONE
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(velmac)
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: problo(SDIM),probhi(SDIM)
      real(amrex_real), INTENT(inout),target :: vel(DIMV(vel),SDIM)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: velmac(DIMV(velmac))
      real(amrex_real), pointer :: velmac_ptr(D_DECL(:,:,:))
      integer, INTENT(in) :: presbc_array(SDIM,2)
      real(amrex_real), INTENT(in) :: outflow_velocity_buffer_size(2*SDIM)
      real(amrex_real) vel_in

      integer i,j,k,ii,jj,kk
      integer dirloc
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_cell(SDIM)
      integer velcomp

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid forcevelocity"
       stop
      endif
      if (bfact.lt.1) then
       print *,"fact invalid"
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
       print *,"dir invalid forcevelocity 2"
       stop
      endif

      vel_ptr=>vel
      velmac_ptr=>velmac
      call checkbound_array(fablo,fabhi,vel_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,velmac_ptr,0,dir)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
       do dirloc=1,SDIM
        xsten_cell(dirloc)=xsten(0,dirloc)
       enddo
       velcomp=dir+1
       vel_in=vel(D_DECL(i,j,k),velcomp)*global_velocity_scale
        ! vel_freestream is declared in: PROB.F90
       call vel_freestream( &
        xsten_cell, &
        dir,vel_in,time, &
        presbc_array, &
        outflow_velocity_buffer_size, &
        problo,probhi)
       vel(D_DECL(i,j,k),velcomp)=vel_in/global_velocity_scale
      enddo
      enddo
      enddo

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
         growlo,growhi,0,dir) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridstenMAC(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf,dir)
       do dirloc=1,SDIM
        xsten_cell(dirloc)=xsten(0,dirloc)
       enddo
       vel_in=velmac(D_DECL(i,j,k))*global_velocity_scale
       call vel_freestream( &
         xsten_cell, &
         dir,vel_in,time, &
         presbc_array, &
         outflow_velocity_buffer_size, &
         problo,probhi)
       velmac(D_DECL(i,j,k))=vel_in/global_velocity_scale
      enddo
      enddo
      enddo

      return
      end subroutine fort_forcevelocity

      subroutine fort_velfill( &
       grid_type, &
       level, &
       u,DIMS(u), &
       domlo,domhi,dx, &
       xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_velfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time

      real(amrex_real), INTENT(inout), target :: u(DIMV(u))
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:))

      integer, INTENT(in) :: bc(SDIM,2)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer velcomp
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (ncomp.ne.1) then
       print *,"ncomp invalid9"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 1"
       stop
      endif

      if ((scomp.ge.0).and.(scomp.lt.SDIM)) then
       ! do nothing
      else
       print *,"scomp invalid"
       stop
      endif
      velcomp=scomp+1

      if ((velcomp.lt.1).or.(velcomp.gt.SDIM)) then
       print *,"velcomp invalid"
       stop
      endif

      u_ptr=>u
      call local_filcc(bfact, &
       u_ptr, &
       domlo,domhi,bc)

      fablo(1)=LBOUND(u,1) ! ulox
      fablo(2)=LBOUND(u,2) ! uloy
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM) ! uloz
#endif
      fabhi(1)=UBOUND(u,1) ! uhix
      fabhi(2)=UBOUND(u,2) ! uhiy
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM) ! uhiz
#endif
      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side)

       if (test_bc.eq.EXT_DIR) then
        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif

       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        if (MARCO.eq.1) then
         print *,"Marco Arienti's code needs to be migrated"
         stop
        endif
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         call velbc_override(time,dir2,side,velcomp, &
          u(D_DECL(i,j,k)), &
          xsten,nhalf,dx,bfact)
        enddo
        enddo
        enddo
       endif            
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_velfill

       ! gets all components of the velocity at once. 
       ! (for just a single material)
      subroutine fort_group_velfill( &
       grid_type, &
       level, &
       u,DIMS(u), &
       domlo,domhi,dx, &
       xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_group_velfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time

      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer velcomp
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((scomp.ge.0).and.(scomp.lt.SDIM)) then
       ! do nothing
      else
       print *,"scomp invalid"
       stop
      endif

      if (ncomp.ne.SDIM) then
       print *,"ncomp invalid10"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 2"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif
      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      u_ptr=>u
      do velcomp=1,SDIM
       call local_filcc4D(bfact, &
        u_ptr,velcomp, &
        domlo,domhi, &
        bc)
      enddo ! velcomp

      do velcomp=1,SDIM
      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side,velcomp)

       if (test_bc.eq.EXT_DIR) then

        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        if (MARCO.eq.1) then
         print *,"Marco Arienti's code needs to be migrated"
         stop
        endif
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         call velbc_override(time,dir2,side,velcomp, &
          u(D_DECL(i,j,k),velcomp), &
          xsten,nhalf,dx,bfact)
        enddo
        enddo
        enddo
       endif            
      enddo
      enddo
      enddo ! dir2,side,velcomp

      return
      end subroutine fort_group_velfill


       ! associated with Solid_State_Type
      subroutine fort_solvfill( &
       grid_type, &
       level, &
       u,DIMS(u), &
       domlo,domhi,dx, &
       xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_solvfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u))
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:))
      integer, INTENT(in) :: bc(SDIM,2)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer velcomp
      integer im_vel
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer nparts

      u_ptr=>u

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (ncomp.ne.1) then
       print *,"ncomp invalid11"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      nparts=0
      do im_vel=1,num_materials
       if (is_lag_part(im_vel).eq.1) then
        nparts=nparts+1
       else if (is_lag_part(im_vel).eq.0) then
        ! do nothing
       else
        print *,"is_lag_part(im_vel) invalid"
        stop
       endif
      enddo
      if ((nparts.lt.1).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid SOLVFILL"
       stop
      endif

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 1"
       stop
      endif
 
      im_vel=scomp/SDIM
      velcomp=scomp-im_vel*SDIM+1

      if ((im_vel.lt.0).or.(im_vel.ge.nparts)) then
       print *,"scomp out of range in solv fill"
       stop
      endif
      if ((velcomp.lt.1).or.(velcomp.gt.SDIM)) then
       print *,"velcomp invalid"
       stop
      endif

      call local_filcc(bfact, &
       u_ptr, &
       domlo,domhi,bc)

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif
      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side)

       if (test_bc.eq.EXT_DIR) then

        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        if (MARCO.eq.1) then
         print *,"Marco Arienti's code needs to be migrated"
         stop
        endif
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         call velbc_override(time,dir2,side,velcomp, &
          u(D_DECL(i,j,k)), &
          xsten,nhalf,dx,bfact)
        enddo
        enddo
        enddo
       endif            
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_solvfill

       ! associated with Solid_State_Type
       ! gets all components of the velocity at once. 
       ! (for just a single material)
      subroutine fort_group_solvfill( &
       grid_type, &
       level, &
       u,DIMS(u), &
       domlo,domhi,dx, &
       xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_group_solvfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)
      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer velcomp
      integer im_vel
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer nparts

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      nparts=0
      do im_vel=1,num_materials
       if (is_lag_part(im_vel).eq.1) then
        nparts=nparts+1
       else if (is_lag_part(im_vel).eq.0) then
        ! do nothing
       else
        print *,"is_lag_part(im_vel) invalid"
        stop
       endif
      enddo
      if ((nparts.lt.1).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_group_solvfill"
       stop
      endif

      im_vel=scomp/SDIM
      if ((im_vel.lt.0).or.(im_vel.ge.nparts).or. &
          (im_vel*SDIM.ne.scomp)) then
       print *,"scomp invalid in fort_group_solvfill"
       stop
      endif
      if (ncomp.ne.SDIM) then
       print *,"ncomp invalid12"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 2"
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif
      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      u_ptr=>u
      do velcomp=1,SDIM
       call local_filcc4D(bfact, &
        u_ptr,velcomp, &
        domlo,domhi, &
        bc)
      enddo ! velcomp

      do velcomp=1,SDIM
      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side,velcomp)

       if (test_bc.eq.EXT_DIR) then

        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        if (MARCO.eq.1) then
         print *,"Marco Arienti's code needs to be migrated"
         stop
        endif
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         call velbc_override(time,dir2,side,velcomp, &
          u(D_DECL(i,j,k),velcomp), &
          xsten,nhalf,dx,bfact)
        enddo
        enddo
        enddo
       endif            
      enddo
      enddo
      enddo ! dir2,side,velcomp

      return
      end subroutine fort_group_solvfill


      subroutine fort_umacfill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_umacfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u))
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:))
      integer, INTENT(in) :: bc(SDIM,2)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer velcomp,veldir
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      u_ptr=>u

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 3"
       stop
      endif

      if ((scomp.lt.0).or.(scomp+ncomp.gt.1)) then
       print *,"scomp invalid fort_umacfill"
       stop
      endif

      if (ncomp.ne.1) then
       print *,"ncomp invalid fort_umacfill 13, grid_type=",grid_type
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      if ((grid_type.ge.0).and.(grid_type.lt.SDIM)) then
       velcomp=grid_type+1
       veldir=grid_type
      else
       print *,"grid_type invalid in umacfill"
       stop
      endif

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (dir2.ne.velcomp) then
        if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
         print *,"domhi+1 not divisible by bfact"
         stop
        endif
       else
        if ((domhi(dir2)/bfact)*bfact.ne.domhi(dir2)) then
         print *,"domhi not divisible by bfact"
         stop
        endif
       endif
      enddo  ! dir2

      call efilcc(bfact, &
       u_ptr, &
       domlo,domhi,bc,veldir) ! veldir=grid_type

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side)

       if (test_bc.eq.EXT_DIR) then

        if (dir2.eq.velcomp) then

         if (side.eq.1) then
          if (fablo(dir2).le.domlo(dir2)) then
           ext_dir_flag=1
           borderhi(dir2)=domlo(dir2)
           inside_index=domlo(dir2)
          endif
         else if (side.eq.2) then
          if (fabhi(dir2).ge.domhi(dir2)) then
           ext_dir_flag=1
           borderlo(dir2)=domhi(dir2)
           inside_index=domhi(dir2)
          endif
         else
          print *,"side invalid"
          stop 
         endif

        else if (dir2.ne.velcomp) then

         if (side.eq.1) then
          if (fablo(dir2).lt.domlo(dir2)) then
           ext_dir_flag=1
           borderhi(dir2)=domlo(dir2)-1
           inside_index=domlo(dir2)
          endif
         else if (side.eq.2) then
          if (fabhi(dir2).gt.domhi(dir2)) then
           ext_dir_flag=1
           borderlo(dir2)=domhi(dir2)+1
           inside_index=domhi(dir2)
          endif
         else
          print *,"side invalid"
          stop
         endif

        else
         print *,"dir2 bust"
         stop
        endif

       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  
    
       if (ext_dir_flag.eq.1) then
        if (MARCO.eq.1) then
         print *,"Marco Arienti's code needs to be migrated"
         stop
        endif
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)
 
          ! velcomp=1..sdim
         call gridstenMAC(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf,velcomp-1)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         call velbc_override(time,dir2,side,velcomp, &
          u(D_DECL(i,j,k)), &
          xsten,nhalf,dx,bfact)
        enddo
        enddo
        enddo
       endif 
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_umacfill

      subroutine fort_moffill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_moffill')

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout) :: u(DIMV(u))
      integer, INTENT(in) :: bc(SDIM,2)

      print *,"fort_moffill should never be called"
      stop

      return
      end subroutine fort_moffill


      subroutine fort_refine_densityfill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_refine_densityfill')

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout) :: u(DIMV(u))
      integer, INTENT(in) :: bc(SDIM,2)

      print *,"fort_refine_densityfill should never be called"
      stop

      return
      end subroutine fort_refine_densityfill


      subroutine fort_extmoffill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_extmoffill')

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout) :: u(DIMV(u))
      integer, INTENT(in) :: bc(SDIM,2)

      print *,"fort_extmoffill should never be called"
      stop

      return
      end subroutine fort_extmoffill

      subroutine fort_group_moffill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_group_moffill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)
      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3
      integer side
      integer side_debug
      integer ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer im
      integer im_debug
      real(amrex_real) uwall(num_materials*ngeom_raw)
      real(amrex_real) uboundary(num_materials*ngeom_raw)
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 6"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      if (num_materials*ngeom_raw.ne.ncomp) then
       print *,"ncomp invalid mof group fill"
       stop
      endif
      if (scomp.ne.STATECOMP_MOF) then
       print *,"scomp invalid mof group fill: ",scomp
       stop
      endif

      u_ptr=>u
      do im=1,num_materials*ngeom_raw
       call local_filcc4D(bfact, &
        u_ptr,im, &
        domlo,domhi,bc)
      enddo

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side,1)

       if (test_bc.eq.EXT_DIR) then

        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         do im=1,num_materials*ngeom_raw
          uwall(im)=u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),im)
          if ((uwall(im).ge.zero).or. &
              (uwall(im).le.zero)) then
           ! do nothing
          else
           print *,"(breakpoint) break point and gdb: "
           print *,"(1) compile with the -g option"
           print *,"(2) break PROB.F90:29153"

           print *,"u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),im)=", &
             u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),im)

           if (1.eq.0) then
            print *,"u(D_DECL(IWALL(1)-1,IWALL(2),IWALL(3)),im)=", &
             u(D_DECL(IWALL(1)-1,IWALL(2),IWALL(3)),im)
           endif

           print *,"uwall(im) is NaN: ",im,uwall(im)
           print *,"IWALL= ",IWALL
           print *,"inside_index= ",inside_index
           print *,"i,j,k ",i,j,k
           print *,"borderlo=",borderlo
           print *,"borderhi=",borderhi
           print *,"fablo=",fablo
           print *,"fabhi=",fabhi
           print *,"domlo=",domlo
           print *,"domhi=",domhi
           print *,"dir2,side ",dir2,side
           print *,"EXT_DIR=",EXT_DIR
           print *,"INT_DIR=",INT_DIR
           print *,"FOEXTRAP=",FOEXTRAP
           print *,"REFLECT_EVEN=",REFLECT_EVEN
           print *,"REFLECT_ODD=",REFLECT_ODD
           do dir3=1,SDIM
           do side_debug=1,2
           do im_debug=1,ncomp
            print *,"dir,side,nc,bc(dir,side,nc) ",dir3,side_debug,im_debug, &
              bc(dir3,side_debug,im_debug)
           enddo 
           enddo 
           enddo 

           stop
          endif
         enddo
         call groupmofBC(time,dir2,side, &
          uboundary, &
          uwall, &
          xsten,nhalf,dx,bfact)
         do im=1,num_materials*ngeom_raw
          u(D_DECL(i,j,k),im)=uboundary(im)
         enddo
        enddo
        enddo
        enddo
       endif 
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_group_moffill


      subroutine fort_group_refine_densityfill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_group_refine_densityfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)
      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer irefine,jrefine,krefine
      integer nrefine
      integer dir2
      integer dir3
      integer side
      integer ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer im
      integer im_compressible
      real(amrex_real) uwall_avg
      real(amrex_real) ughost
      integer, parameter :: istate=1
      integer, parameter :: scomp_data=1
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in group_refine_densityfill"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid group_refine_densityfill"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid group_refine_densityfill"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid group_refine_densityfill"
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      if (4*(SDIM-1).ne.ncomp) then
       print *,"ncomp invalid group_refine_densityfill"
       stop
      endif

      if ((num_materials_compressible.ge.1).and. &
          (num_materials_compressible.le.num_materials)) then
       ! do nothing
      else
       print *,"num_materials_compressible invalid group_refine_densityfill"
       stop
      endif

      if ((scomp.ge.0).and. &
          (scomp.le.(num_materials_compressible-1)*ncomp)) then
       !do nothing
      else
       print *,"scomp invalid group_refine_densityfill: ",scomp
       stop
      endif

      im_compressible=scomp/ncomp
      im_compressible=im_compressible+1
      if (ncomp*(im_compressible-1).ne.scomp) then
       print *,"scomp invalid group_refine_densityfill(2): ",scomp
       stop
      endif
      if ((im_compressible.ge.1).and. &
          (im_compressible.le.num_materials_compressible)) then
       ! do nothing
      else
       print *,"im_compressible invalid: ",im_compressible
       stop
      endif

      u_ptr=>u
      call local_filcc4D_refine(bfact, &
       u_ptr, &
       scomp_data,ncomp, &
       domlo,domhi,bc)

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side,1)

       if (test_bc.eq.EXT_DIR) then

        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid fort_group_refine_densityfill"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid fort_group_refine_densityfill: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         im=fort_im_refine_density_map(im_compressible)+1
         if ((im.ge.1).and.(im.le.num_materials)) then
          !do nothing
         else
          print *,"im invalid (im,im_compressible): ",im,im_compressible
          stop
         endif

         uwall_avg=zero 
         krefine=0
#if (AMREX_SPACEDIM==3)
         do krefine=0,1
#endif
         do jrefine=0,1
         do irefine=0,1
          nrefine=4*krefine+2*jrefine+irefine+1
          uwall_avg=uwall_avg+u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),nrefine)
         enddo
         enddo
#if (AMREX_SPACEDIM==3)
         enddo
#endif
         uwall_avg=uwall_avg/ncomp
         call denBC(time,dir2,side, &
           ughost, &
           uwall_avg, &
           xsten,nhalf,dx,bfact,istate,im)

         krefine=0
#if (AMREX_SPACEDIM==3)
         do krefine=0,1
#endif
         do jrefine=0,1
         do irefine=0,1
          nrefine=4*krefine+2*jrefine+irefine+1
          u(D_DECL(i,j,k),nrefine)=ughost
         enddo
         enddo
#if (AMREX_SPACEDIM==3)
         enddo
#endif

        enddo
        enddo
        enddo
       endif            

      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_group_refine_densityfill


      subroutine fort_group_extmoffill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_group_extmoffill')

      use filcc_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer isub
      integer ibasesrc
      integer ibasedst
      integer, parameter :: use_ls_data=0
      integer, parameter :: mof_verbose=0
      integer, parameter :: continuous_mof=STANDARD_MOF
      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))

      integer :: grid_index(SDIM)
      integer, parameter :: grid_level=-1

      real(amrex_real) LS_stencil(D_DECL(-1:1,-1:1,-1:1),1)  ! not used
      real(amrex_real) multi_centroidA(num_materials,SDIM)
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) vof_super(num_materials)
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer im
      real(amrex_real) uwall(num_materials*ngeom_raw)
      real(amrex_real) uboundary(num_materials*ngeom_raw)
      integer vofcomp
      real(amrex_real) voffluid_wall,vofsolid_wall
      real(amrex_real) voffluid_bound,vofsolid_bound
      real(amrex_real) voftest_wall,voftest_bound

      integer tessellate

      integer, parameter :: nhalf=3

      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer tid,nmax
#ifdef _OPENMP
      integer omp_get_thread_num
#endif

      nmax=POLYGON_LIST_MAX ! in: fort_group_extmoffill

      tid=0       
#ifdef _OPENMP
      tid=omp_get_thread_num()
#endif
      if ((tid.ge.geom_nthreads).or.(tid.lt.0)) then
       print *,"tid invalid"
       stop
      endif 
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 7"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      tessellate=0

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      if (num_materials*ngeom_recon.ne.ncomp) then
       print *,"ncomp invalid in fort_group_extmoffill"
       stop
      endif
      if (scomp.ne.EXTRAPCOMP_MOF) then
       print *,"scomp invalid in fort_group_extmoffill"
       stop
      endif

      u_ptr=>u
      do im=1,num_materials*ngeom_recon
       call local_filcc4D(bfact, &
        u_ptr,im, &
        domlo,domhi,bc)
      enddo

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side,1)

       if (test_bc.eq.EXT_DIR) then
        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         grid_index(1)=i
         grid_index(2)=j
         if (SDIM.eq.3) then
          grid_index(SDIM)=k
         endif

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         do im=1,num_materials
          ibasesrc=(im-1)*ngeom_recon
          ibasedst=(im-1)*ngeom_raw
          do isub=1,ngeom_raw
           uwall(ibasedst+isub)= &
             u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),ibasesrc+isub)
          enddo 
         enddo ! im

         call groupmofBC(time,dir2,side, &
          uboundary, &
          uwall, &
          xsten,nhalf,dx,bfact)

         voffluid_wall=zero
         vofsolid_wall=zero
         voffluid_bound=zero
         vofsolid_bound=zero
         do im=1,num_materials
          vofcomp=(im-1)*ngeom_raw+1
          voftest_wall=uwall(vofcomp)
          voftest_bound=uboundary(vofcomp)
          if (is_rigid(im).eq.0) then
           voffluid_wall=voffluid_wall+voftest_wall
           voffluid_bound=voffluid_bound+voftest_bound
          else if (is_rigid(im).eq.1) then
           vofsolid_wall=vofsolid_wall+voftest_wall
           vofsolid_bound=vofsolid_bound+voftest_bound
          else
           print *,"is_rigid invalid PROB.F90"
           stop
          endif
         enddo ! im=1..num_materials

         if ((voffluid_wall.le.zero).or.(voffluid_bound.le.zero)) then
          print *,"fluid disappeared in GROUP_EXTMOFFILL"
          print *,"i,j,k ",i,j,k
          print *,"IWALL ",IWALL(1),IWALL(2),IWALL(3)
          print *,"voffluid_wall,vofsolid_wall ",voffluid_wall,vofsolid_wall
          print *,"voffluid_bound,vofsolid_bound ",voffluid_bound, &
            vofsolid_bound
          stop
         else if ((voffluid_wall.le.two).and.(voffluid_bound.le.two)) then
          ! do nothing
         else
          print *,"voffluid_wall or voffluid_bound invalid"
          stop
         endif

         do im=1,num_materials

          ibasesrc=(im-1)*ngeom_raw+1
          ibasedst=(im-1)*ngeom_recon+1
          do dir3=0,SDIM
           mofdata(ibasedst+dir3)=uboundary(ibasesrc+dir3)
          enddo

           !slope=0
          do dir3=1,SDIM
           mofdata(ibasedst+SDIM+1+dir3)=zero
          enddo

           ! order=0
          mofdata(ibasedst+SDIM+1)=zero

         enddo  ! im=1,num_materials

         call make_vfrac_sum_ok_base( &
           cmofsten, &
           xsten,nhalf, &
           continuous_mof, &
           bfact,dx, &
           tessellate, &  ! =0
           mofdata, &
           SDIM)

         do im=1,num_materials
          ibasedst=(im-1)*ngeom_recon+1
          vof_super(im)=mofdata(ibasedst)
         enddo

         call multimaterial_MOF( &
          bfact,dx,xsten,nhalf, &
          mof_verbose, & !=0
          use_ls_data, & !=0
          LS_stencil, & 
          geom_xtetlist(1,1,1,tid+1), &
          geom_xtetlist(1,1,1,tid+1), &
          nmax, &
          nmax, &
          mofdata, & !intent(inout)
          vof_super, &
          multi_centroidA, &
          continuous_mof, & !=STANDARD_MOF
          cmofsten, &
          grid_index, &
          grid_level, &
          SDIM)

         do dir3=1,num_materials*ngeom_recon
          u(D_DECL(i,j,k),dir3)=mofdata(dir3)
         enddo

        enddo
        enddo
        enddo
       endif            
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_group_extmoffill

      subroutine fort_ls_fill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_ls_fill')

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout) :: u(DIMV(u))
      integer, INTENT(in) :: bc(SDIM,2)

      print *,"fort_ls_fill should never be called"
      print *,"fort_group_ls_fill should be called instead"
      stop

      return
      end subroutine fort_ls_fill

      subroutine fort_group_ls_fill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_group_ls_fill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM) 
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer imls
      
      real(amrex_real) uwall(ncomp)
      real(amrex_real) uboundary(ncomp)
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer ncomp_ho,icomp

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 10"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      ncomp_ho=(SDIM+1)*num_materials

      if (ncomp_ho.ne.ncomp) then
       print *,"ncomp invalid fort_group_ls_fill"
       stop
      endif
      if ((scomp.ne.SDIM*num_materials).and. &  ! called from ghost fill
          (scomp.ne.0)) then           ! called from main fill
       print *,"scomp invalid fort_group_ls_fill"
       print *,"scomp= ",scomp
       stop
      endif
      if ((ls_homflag.ne.0).and.(ls_homflag.ne.1)) then
       print *,"ls_homflag invalid"
       stop
      endif

      u_ptr=>u
      do imls=1,ncomp_ho
       call local_filcc4D( &
        bfact, &
        u_ptr,imls, &
        domlo,domhi, &
        bc)
      enddo

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

        ! 1..num_materials  level set functions
       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo

       ext_dir_flag=0

       test_bc=bc(dir2,side,1)

       if (test_bc.eq.EXT_DIR) then

        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index
         do imls=1,num_materials
          uwall(imls)=u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),imls)
         enddo
         call grouplsBC(time,dir2,side, &
          uboundary, &
          uwall, &
          xsten,nhalf,dx,bfact)
         do imls=1,num_materials
          u(D_DECL(i,j,k),imls)=uboundary(imls)
         enddo
        enddo
        enddo
        enddo
       else if (ext_dir_flag.eq.0) then
        ! do nothing
       else
        print *,"ext_dir_flag invalid"
        stop
       endif 

        ! num_materials+1 ... ncomp_ho: levelset normals
       do icomp=num_materials+1,ncomp_ho

        borderlo(3)=0
        borderhi(3)=0
        do dir3=1,SDIM
         borderlo(dir3)=fablo(dir3)
         borderhi(dir3)=fabhi(dir3)
        enddo
        ext_dir_flag=0

        test_bc=bc(dir2,side,icomp)

        if (test_bc.eq.EXT_DIR) then

         if (side.eq.1) then
          if (fablo(dir2).lt.domlo(dir2)) then
           ext_dir_flag=1
           borderhi(dir2)=domlo(dir2)-1
           inside_index=domlo(dir2)
          endif
         else if (side.eq.2) then
          if (fabhi(dir2).gt.domhi(dir2)) then
           ext_dir_flag=1
           borderlo(dir2)=domhi(dir2)+1
           inside_index=domhi(dir2)
          endif
         else
          print *,"side invalid"
          stop
         endif
        else if ((test_bc.eq.FOEXTRAP).or. &
                 (test_bc.eq.HOEXTRAP).or. &
                 (test_bc.eq.REFLECT_EVEN).or. &
                 (test_bc.eq.REFLECT_ODD).or. &
                 (test_bc.eq.INT_DIR)) then
         ! do nothing
        else
         print *,"test_bc invalid: ",test_bc
         stop
        endif  

        if (ext_dir_flag.eq.1) then
         do k=borderlo(3),borderhi(3)
         do j=borderlo(2),borderhi(2)
         do i=borderlo(1),borderhi(1)

          call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
          IWALL(1)=i
          IWALL(2)=j
          IWALL(3)=k
          IWALL(dir2)=inside_index

          call extrapBC(time,dir2,side, &
            u(D_DECL(i,j,k),icomp), &
            u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),icomp), &
            xsten,nhalf,dx,bfact)
         enddo
         enddo
         enddo
        else if (ext_dir_flag.eq.0) then
         ! do nothing
        else
         print *,"ext_dir_flag invalid"
         stop
        endif 

       enddo ! icomp=num_materials+1 ... ncomp_ho

      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_group_ls_fill

       ! this is for the "errorind" variable.
      subroutine fort_scalarfill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_scalarfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u))
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:))

      integer, INTENT(in) :: bc(SDIM,2)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer nc
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 11"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      if (ncomp.ne.1) then
       print *,"ncomp invalid in fort_scalarfill"
       stop
      endif
       ! "errorind" variable.
      nc=STATECOMP_ERR
      if (scomp.ne.nc) then
       print *,"scomp invalid in fort_scalarfill"
       stop
      endif

      u_ptr=>u
      call local_filcc(bfact, &
       u_ptr, &
       domlo,domhi,bc)

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo

       ext_dir_flag=0

       test_bc=bc(dir2,side)

       if (test_bc.eq.EXT_DIR) then
        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         call scalarBC(time,dir2,side, &
           u(D_DECL(i,j,k)), &
           u(D_DECL(IWALL(1),IWALL(2),IWALL(3))), &
           xsten,nhalf,dx,bfact)
        enddo
        enddo
        enddo
       endif            
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_scalarfill

      subroutine fort_statefill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_statefill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u))
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:))

      integer, INTENT(in) :: bc(SDIM,2)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer im,istate
      integer icomplo,icomphi
      integer scomp_spec,num_state_material_test
      integer dencomp
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 14"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      if (ncomp.ne.1) then
       print *,"ncomp invalid15"
       stop
      endif
       ! c++ index
      icomplo=STATECOMP_STATES
      icomphi=icomplo+num_materials*num_state_material
      if ((scomp.lt.icomplo).or.(scomp.ge.icomphi)) then
       print *,"scomp out of range in fort_statefill"
       stop
      endif
      im=(scomp-icomplo)/num_state_material+1
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im out of range"
       stop
      endif
       ! c++ convention
      dencomp=icomplo+(im-1)*num_state_material+ENUM_DENVAR
  
       ! fortran convention
      istate=scomp-dencomp+1
      if ((istate.lt.1).or.(istate.gt.num_state_material)) then
       print *,"istate invalid"
       stop
      endif

        ! fortran convention
      scomp_spec=num_state_base+1
      num_state_material_test=scomp_spec+num_species_var-1

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
   
      if (num_state_material_test.ne.num_state_material) then
       print *,"num_state_material invalid"
       stop
      endif

      u_ptr=>u
      call local_filcc(bfact, &
       u_ptr, &
       domlo,domhi,bc)

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side)

       if (test_bc.eq.EXT_DIR) then
        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         if ((istate.ge.1).and. &
             (istate.le.num_state_material)) then
          if ((im.ge.1).and.(im.le.num_materials)) then
           call denBC(time,dir2,side, &
            u(D_DECL(i,j,k)), &
            u(D_DECL(IWALL(1),IWALL(2),IWALL(3))), &
            xsten,nhalf,dx,bfact,istate,im)
          else
           print *,"im invalid (call to denBC): ",im
           stop
          endif
         else
          print *,"istate not supported"
          stop
         endif
        enddo
        enddo
        enddo
       endif            
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_statefill

      subroutine fort_tensorfill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_tensorfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u))
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:))

      integer, INTENT(in) :: bc(SDIM,2)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer im,ipart
      integer icomplo,icomphi
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 14"
       stop
      endif

      if (ENUM_NUM_TENSOR_TYPE.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ENUM_NUM_TENSOR_TYPE invalid"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      if (ncomp.ne.1) then
       print *,"ncomp invalid16"
       stop
      endif
       ! c++ index
      icomplo=0
      icomphi=num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE
      if ((scomp.lt.icomplo).or.(scomp.ge.icomphi)) then
       print *,"scomp out of range in fort_tensorfill"
       stop
      endif

      ipart=(scomp-icomplo)/ENUM_NUM_TENSOR_TYPE+1
      if ((ipart.lt.1).or. &
          (ipart.gt.num_materials_viscoelastic)) then
       print *,"ipart out of range:FORT_TENSORFILL"
       stop
      endif
      im=fort_im_elastic_map(ipart)+1
      if ((im.ge.1).and.(im.le.num_materials)) then

       u_ptr=>u
       call local_filcc(bfact, &
        u_ptr, &
        domlo,domhi,bc)

       fablo(1)=LBOUND(u,1)
       fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
       fablo(SDIM)=LBOUND(u,SDIM)
#endif
       fabhi(1)=UBOUND(u,1)
       fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
       fabhi(SDIM)=UBOUND(u,SDIM)
#endif

       do dir2=1,SDIM
        if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
         print *,"domlo not divisible by bfact"
         stop
        endif
        if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
         print *,"domhi+1 not divisible by bfact"
         stop
        endif
       enddo  ! dir2

       do dir2=1,SDIM
       do side=1,2

        borderlo(3)=0
        borderhi(3)=0
        do dir3=1,SDIM
         borderlo(dir3)=fablo(dir3)
         borderhi(dir3)=fabhi(dir3)
        enddo
        ext_dir_flag=0

        test_bc=bc(dir2,side)

        if (test_bc.eq.EXT_DIR) then
         if (side.eq.1) then
          if (fablo(dir2).lt.domlo(dir2)) then
           ext_dir_flag=1
           borderhi(dir2)=domlo(dir2)-1
           inside_index=domlo(dir2)
          endif
         else if (side.eq.2) then
          if (fabhi(dir2).gt.domhi(dir2)) then
           ext_dir_flag=1
           borderlo(dir2)=domhi(dir2)+1
           inside_index=domhi(dir2)
          endif
         else
          print *,"side invalid"
          stop
         endif
        else if ((test_bc.eq.FOEXTRAP).or. &
                 (test_bc.eq.HOEXTRAP).or. &
                 (test_bc.eq.REFLECT_EVEN).or. &
                 (test_bc.eq.REFLECT_ODD).or. &
                 (test_bc.eq.INT_DIR)) then
         ! do nothing
        else
         print *,"test_bc invalid: ",test_bc
         stop
        endif  

        if (ext_dir_flag.eq.1) then

         print *,"expecting all BCs to be reflect even or odd"
         stop

         do k=borderlo(3),borderhi(3)
         do j=borderlo(2),borderhi(2)
         do i=borderlo(1),borderhi(1)

          call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
          IWALL(1)=i
          IWALL(2)=j
          IWALL(3)=k
          IWALL(dir2)=inside_index

           ! sets all the EXT_DIR bcs to 0.0
          call tensorBC(time,dir2,side, &
           u(D_DECL(i,j,k)), &
           u(D_DECL(IWALL(1),IWALL(2),IWALL(3))), &
           xsten,nhalf,dx,bfact,ipart,im)
         enddo
         enddo
         enddo
        endif            
       enddo ! side
       enddo ! dir2
      else
       print *,"im invalid in fort_tensorfill"
       stop
      endif

      return
      end subroutine fort_tensorfill

      subroutine fort_pressurefill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_pressurefill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u))
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:))

      integer, INTENT(in) :: bc(SDIM,2)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif
      if (ncomp.ne.1) then
       print *,"ncomp invalid17"
       stop
      endif
      if ((scomp.ne.STATECOMP_PRES).and.(scomp.ne.0)) then 
       print *,"scomp invalid fort_pressurefill"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 15"
       stop
      endif

      u_ptr=>u
      call local_filcc(bfact, &
       u_ptr, &
       domlo,domhi,bc)

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      do dir2=1,SDIM
       if ((domlo(dir2)/bfact)*bfact.ne.domlo(dir2)) then
        print *,"domlo not divisible by bfact"
        stop
       endif
       if (((domhi(dir2)+1)/bfact)*bfact.ne.domhi(dir2)+1) then
        print *,"domhi+1 not divisible by bfact"
        stop
       endif
      enddo  ! dir2

      do dir2=1,SDIM
      do side=1,2
       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side)

       if (test_bc.eq.EXT_DIR) then

        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         call presBDRYCOND(time,dir2,side, &
           u(D_DECL(i,j,k)), &
           u(D_DECL(IWALL(1),IWALL(2),IWALL(3))), &
           xsten,nhalf,dx,bfact)
        enddo
        enddo
        enddo
       endif            
      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_pressurefill

      subroutine fort_group_statefill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc,scomp,ncomp,bfact) &
      bind(c,name='fort_group_statefill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer im,icomp,istate
      integer scomp_spec,num_state_material_test
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

       ! c++ index
      if (scomp.ne.STATECOMP_STATES) then
       print *,"scomp invalid fort_group_statefill"
       stop
      endif
      if (ncomp.ne.num_state_material*num_materials) then
       print *,"ncomp invalid18"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 16"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

        ! fortran convention
      scomp_spec=num_state_base+1
      num_state_material_test=scomp_spec+num_species_var-1

      if (num_state_material_test.ne.num_state_material) then
       print *,"num_state_material invalid"
       stop
      endif
      

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      u_ptr=>u
      do icomp=1,num_state_material*num_materials
       call local_filcc4D(bfact, &
        u_ptr,icomp, &
        domlo,domhi,bc)
      enddo

      do dir2=1,SDIM
      do side=1,2

      icomp=0 
      do im=1,num_materials
      do istate=1,num_state_material
       icomp=icomp+1

       borderlo(3)=0
       borderhi(3)=0
       do dir3=1,SDIM
        borderlo(dir3)=fablo(dir3)
        borderhi(dir3)=fabhi(dir3)
       enddo
       ext_dir_flag=0

       test_bc=bc(dir2,side,icomp)

       if (test_bc.eq.EXT_DIR) then
        if (side.eq.1) then
         if (fablo(dir2).lt.domlo(dir2)) then
          ext_dir_flag=1
          borderhi(dir2)=domlo(dir2)-1
          inside_index=domlo(dir2)
         endif
        else if (side.eq.2) then
         if (fabhi(dir2).gt.domhi(dir2)) then
          ext_dir_flag=1
          borderlo(dir2)=domhi(dir2)+1
          inside_index=domhi(dir2)
         endif
        else
         print *,"side invalid"
         stop
        endif
       else if ((test_bc.eq.FOEXTRAP).or. &
                (test_bc.eq.HOEXTRAP).or. &
                (test_bc.eq.REFLECT_EVEN).or. &
                (test_bc.eq.REFLECT_ODD).or. &
                (test_bc.eq.INT_DIR)) then
        ! do nothing
       else
        print *,"test_bc invalid: ",test_bc
        stop
       endif  

       if (ext_dir_flag.eq.1) then
        do k=borderlo(3),borderhi(3)
        do j=borderlo(2),borderhi(2)
        do i=borderlo(1),borderhi(1)

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
         IWALL(1)=i
         IWALL(2)=j
         IWALL(3)=k
         IWALL(dir2)=inside_index

         if ((istate.ge.1).and. &
             (istate.le.num_state_material)) then
          if ((im.ge.1).and.(im.le.num_materials)) then
           call denBC(time,dir2,side, &
            u(D_DECL(i,j,k),icomp), &
            u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),icomp), &
            xsten,nhalf,dx,bfact,istate,im)
          else
           print *,"im invalid (call to denBC2): ",im
           stop
          endif
         else
          print *,"istate not supported"
          stop
         endif
        enddo
        enddo
        enddo
       endif            
      enddo ! istate
      enddo ! im

      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_group_statefill

      subroutine fort_group_tensorfill( &
      grid_type, &
      level, &
      u,DIMS(u), &
      domlo,domhi,dx, &
      xlo,time,bc, &
      scomp,ncomp,bfact) &
      bind(c,name='fort_group_tensorfill')

      use filcc_module
      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: grid_type
      integer, INTENT(in) :: scomp,ncomp,bfact,level
      integer, INTENT(in) :: DIMDEC(u)  ! ulox,uloy,uloz,uhix,uhiy,uhiz
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM), time
      real(amrex_real), INTENT(inout), target :: u(DIMV(u),ncomp)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: bc(SDIM,2,ncomp)
      integer :: test_bc

      integer i,j,k
      integer dir2,dir3,side,ext_dir_flag,inside_index
      integer fablo(SDIM)
      integer fabhi(SDIM)
      integer borderlo(3)
      integer borderhi(3)
      integer IWALL(3)
      integer ipart,im,istate
      integer icomp_total
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      integer check_scomp,check_ncomp,max_ncomp
      
      if (ENUM_NUM_TENSOR_TYPE.eq.2*SDIM) then
       ! do nothing
      else
       print *,"ENUM_NUM_TENSOR_TYPE invalid"
       stop
      endif
      if (grid_type.eq.-1) then
       ! do nothing
      else
       print *,"grid_type invalid"
       stop
      endif

      max_ncomp=num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE

      check_scomp=0
      check_ncomp=0
      if (scomp.eq.0) then
       check_scomp=1
      endif
      if (ncomp.eq.max_ncomp) then
       check_ncomp=1
      endif

       ! c++ index
      if (check_scomp.eq.1) then
       ! do nothing
      else
       print *,"scomp invalid fort_group_tensorfill"
       print *,"scomp=",scomp
       print *,"ncomp=",ncomp
       print *,"num_materials_viscoelastic=",num_materials_viscoelastic
       print *,"ENUM_NUM_TENSOR_TYPE= ",ENUM_NUM_TENSOR_TYPE
       print *,"level=",level
       print *,"fort_finest_level=",fort_finest_level
       stop
      endif
      if (check_ncomp.eq.1) then
       ! do nothing
      else
       print *,"ncomp invalid19 ncomp=",ncomp
       print *,"num_materials_viscoelastic=",num_materials_viscoelastic
       print *,"scomp=",scomp
       stop
      endif
      if (scomp+ncomp.le.max_ncomp) then
       ! do nothing
      else
       print *,"scomp+ncomp invalid19"
       stop
      endif

      if ((level.lt.0).or.(level.gt.fort_finest_level)) then
       print *,"level invalid in fill 16"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid200"
       stop
      endif

      fablo(1)=LBOUND(u,1)
      fablo(2)=LBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fablo(SDIM)=LBOUND(u,SDIM)
#endif
      fabhi(1)=UBOUND(u,1)
      fabhi(2)=UBOUND(u,2)
#if (AMREX_SPACEDIM==3)
      fabhi(SDIM)=UBOUND(u,SDIM)
#endif

      u_ptr=>u
      do icomp_total=scomp+1,scomp+ncomp
       call local_filcc4D(bfact, &
        u_ptr,icomp_total-scomp, &
        domlo,domhi,bc)
      enddo

      do dir2=1,SDIM
      do side=1,2

       icomp_total=0

       if (scomp.eq.0) then
 
        do ipart=1,num_materials_viscoelastic
        do istate=1,ENUM_NUM_TENSOR_TYPE

         icomp_total=icomp_total+1

         if ((icomp_total.ge.1).and. &
             (icomp_total.le.ncomp).and. &
             (icomp_total.le. &
              num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE)) then

          im=fort_im_elastic_map(ipart)+1

          borderlo(3)=0
          borderhi(3)=0
          do dir3=1,SDIM
           borderlo(dir3)=fablo(dir3)
           borderhi(dir3)=fabhi(dir3)
          enddo
          ext_dir_flag=0

          test_bc=bc(dir2,side,icomp_total)

          if (test_bc.eq.EXT_DIR) then

           if (side.eq.1) then
            if (fablo(dir2).lt.domlo(dir2)) then
             ext_dir_flag=1
             borderhi(dir2)=domlo(dir2)-1
             inside_index=domlo(dir2)
            endif
           else if (side.eq.2) then
            if (fabhi(dir2).gt.domhi(dir2)) then
             ext_dir_flag=1
             borderlo(dir2)=domhi(dir2)+1
             inside_index=domhi(dir2)
            endif
           else
            print *,"side invalid"
            stop
           endif
          else if ((test_bc.eq.FOEXTRAP).or. &
                   (test_bc.eq.HOEXTRAP).or. &
                   (test_bc.eq.REFLECT_EVEN).or. &
                   (test_bc.eq.REFLECT_ODD).or. &
                   (test_bc.eq.INT_DIR)) then
           ! do nothing
          else
           print *,"test_bc invalid: ",test_bc
           stop
          endif  

          if (ext_dir_flag.eq.1) then

           print *,"in: fort_group_tensorfill:"
           print *,"expecting all BCs to be reflect even or odd"

           stop

           do k=borderlo(3),borderhi(3)
           do j=borderlo(2),borderhi(2)
           do i=borderlo(1),borderhi(1)
  
            call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
            IWALL(1)=i
            IWALL(2)=j
            IWALL(3)=k
            IWALL(dir2)=inside_index

              ! sets all the physical BCs to 0.0
            call tensorBC(time,dir2,side, &
             u(D_DECL(i,j,k),icomp_total), &
             u(D_DECL(IWALL(1),IWALL(2),IWALL(3)),icomp_total), &
             xsten,nhalf,dx,bfact,ipart,im)
           enddo
           enddo
           enddo
          else if (ext_dir_flag.eq.0) then
           ! do nothing
          else
           print *,"ext_dir_flag invalid"
           stop
          endif  

         else
          print *,"icomp_total invalid"
          stop
         endif

        enddo ! istate
        enddo ! ipart
       else
        print *,"scomp invalid for tensor bc"
        stop
       endif

      enddo ! side
      enddo ! dir2

      return
      end subroutine fort_group_tensorfill

      end module probf90_module

