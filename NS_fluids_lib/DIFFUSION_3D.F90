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
#include "DIFFUSION_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif
!
! 1. DA/Dt=0  or DQ/Dt=0
! 2. A^** = S A^* S^T   S=I+dt grad u
! 3. dA/dt = -(A-I)/lambda  or DQ/Dt=-Q/lambda
!
! force term: div(mu H Q)/rho or div(mu H A)/rho  Q=A-I
!
!ux,vx,wx,uy,vy,wy,uz,vz,wz
! grad u in cylindrical coordinates:
!
! S= (grad u + grad u^T)/2 
!
! grad u=| u_r  u_t/r-v/r  u_z  |
!        | v_r  v_t/r+u/r  v_z  |
!        | w_r  w_t/r      w_z  |
! in RZ:
! grad u=| u_r  0  u_z  |
!        | 0   u/r  0   |
!        | w_r  0  w_z  |
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
! ur =     costheta u + sintheta v
! utheta = -sintheta u + costheta v
! 
! u = costheta ur - sintheta utheta
! v = sintheta ur + costheta utheta
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
! the above (constant viscosity) matches mathematica website.
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
!
! Rotating frame of reference: (X-Y-Z coordinates)
!
! Eady (page 35) and Lappa (page 20):
!  \vec{u}_{t} = - 2 \vec{Omega} x \vec{u} - grad p/rho - |g|\vec{z}
!   i       j        k
!   0       0        Omega
!   u       v        w
!
! u_t = 2 Omega v - p_{x}/rho0 
! v_t = -2 Omega u - p_{y}/rho0
! w_t = -p_{z}/rho0 + |g beta|(T-T0)
! T_t + u T_x + v T_y +w T_z=0
!
! K=2 Omega
! Base state: ubase=\Gamma z 
!             Tbase=T0 + A y + B z  
!             \Gamma=-A |beta g|/(2 Omega)
!             pbase=rho0 |beta g| B z^{2}/2 + rho0 |beta g| A y z 
!
! check: pbase_{x}=vbase=0
!        pbase_{y}/rho0=|beta g|A z 
!        -2 Omega ubase=-2 Omega \Gamma z=A |beta g| z
!        pbase_{z}/rho0=|beta g| B z + |beta g|A y 
! 
! Lewis and Nagata (ignoring viscosity and nonlinear terms):
! \vec{u}=u e_r + v e_{theta} + w e_{z}
! e_{z} x \vec{u}= e_r     e_theta     e_z
!                   0        0          1
!                   u        v          w  = -v e_r + u e_theta
! u_t= 2 Omega v - p_r/rho0
! v_t=-2 Omega u - p_theta/rho0
! w_t=           - p_z/rho0 + |g beta| (T-T0)
! T_t + u T_r + v T_theta +w T_z=0
!
! Base state: vbase=\Gamma z 
!             Tbase=T0 + A r + B z  
!             \Gamma=A |beta g|/(2 Omega)
!             pbase=rho0 |beta g| B z^{2}/2 + rho0 |beta g| A r z 
!
! check: pbase_{theta}=ubase=0
!        pbase_{r}/rho0=A |beta g| z 
!        2 Omega vbase=2 Omega \Gamma z=A |beta g| z
!        pbase_{z}/rho0=|beta g| B z + |beta g| A r 
!
! sanity check for cylindrical coordinates: suppose particle has
! velocity \vec{u} = (-1, 0, 0), then particle will deflect counter clockwise
! (same direction as Omega if Omega>0)
! i.e. new velocity will be (-1, 2 Omega dt,   0)
!
! Let p'=p-pbase  
! u_t= 2 Omega (v-vbase+vbase) - (p'+pbase)_r/rho0
! v_t=-2 Omega (u-ubase+ubase) - (p'+pbase)_theta/rho0
! w_t=           - (p'+pbase)_z/rho0 + |g beta| (T-Tbase+Tbase-T0)
!
! u_t= 2 Omega (v-vbase) - (p')_r/rho0
! v_t=-2 Omega (u-ubase) - (p')_theta/rho0
! w_t=           - (p')_z/rho0 + |g beta| (T-Tbase)
! p'=0 T=Tbase u=ubase v=vbase w=0 is a solution.

       module hoop_module
       use amrex_fort_module, only : amrex_real

       contains

       subroutine fort_hoopimplicit( &
         im_elastic_map, &
         num_FSI_outer_sweeps, &
         FSI_outer_sweeps, &
         override_density, &
         constant_density_all_time, & ! 1..num_materials
         force,DIMS(force), &
         tensor,DIMS(tensor), &
         thermal,DIMS(thermal), &
         recon,DIMS(recon), &
         solxfab,DIMS(solxfab), &
         solyfab,DIMS(solyfab), &
         solzfab,DIMS(solzfab), &
         xlo,dx, &
         uold,DIMS(uold), &
         unew,DIMS(unew), &
         lsnew,DIMS(lsnew), &
         den,DIMS(den), &  ! 1/density
         mu,DIMS(mu), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         finest_level, &
         visc_coef, &
         angular_velocity, & !intent(in): fort_hoopimplicit
         centrifugal_force_factor, & !intent(in): fort_hoopimplicit
         uncoupled_viscosity, &
         update_state, &
         dt, &
         cur_time_slab, &
         rzflag, &
         nparts, &
         nparts_def, &
         im_solid_map) &
       bind(c,name='fort_hoopimplicit')

       use probf90_module
       use global_utility_module 

       IMPLICIT NONE

       integer, INTENT(in) :: num_FSI_outer_sweeps
       integer, INTENT(in) :: FSI_outer_sweeps
       integer, INTENT(in) :: im_elastic_map(num_FSI_outer_sweeps-1)
       integer, INTENT(in) :: override_density(num_materials)
       integer, INTENT(in) :: constant_density_all_time(num_materials)
 
       integer, INTENT(in) :: nparts
       integer, INTENT(in) :: nparts_def
       integer, INTENT(in) :: im_solid_map(nparts_def)
       integer, INTENT(in) :: level
       integer, INTENT(in) :: finest_level
       integer, INTENT(in) :: rzflag
       real(amrex_real), INTENT(in) :: angular_velocity
       real(amrex_real), INTENT(in) :: centrifugal_force_factor
       real(amrex_real), INTENT(in) :: visc_coef
       integer, INTENT(in) :: uncoupled_viscosity
       integer, INTENT(in) :: update_state
       real(amrex_real), INTENT(in) :: dt
       real(amrex_real), INTENT(in) :: cur_time_slab
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer :: growlo(3),growhi(3)
       integer, INTENT(in) :: bfact
    
       integer, INTENT(in) :: DIMDEC(force)
       integer, INTENT(in) :: DIMDEC(tensor)
       integer, INTENT(in) :: DIMDEC(thermal)
       integer, INTENT(in) :: DIMDEC(recon)
       integer, INTENT(in) :: DIMDEC(solxfab)
       integer, INTENT(in) :: DIMDEC(solyfab)
       integer, INTENT(in) :: DIMDEC(solzfab)
       integer, INTENT(in) :: DIMDEC(uold)
       integer, INTENT(in) :: DIMDEC(unew)
       integer, INTENT(in) :: DIMDEC(lsnew)
       integer, INTENT(in) :: DIMDEC(den)
       integer, INTENT(in) :: DIMDEC(mu)

       real(amrex_real), INTENT(out),target :: &
               force(DIMV(force),AMREX_SPACEDIM)
       real(amrex_real), pointer :: force_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in),target :: &
               tensor(DIMV(tensor),AMREX_SPACEDIM_SQR)
       real(amrex_real), pointer :: tensor_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in),target :: thermal(DIMV(thermal))
       real(amrex_real), pointer :: thermal_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in),target :: &
               recon(DIMV(recon),num_materials*ngeom_recon)
       real(amrex_real), pointer :: recon_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in),target :: &
               solxfab(DIMV(solxfab),nparts_def*SDIM)
       real(amrex_real), INTENT(in),target :: &
               solyfab(DIMV(solyfab),nparts_def*SDIM)
       real(amrex_real), INTENT(in),target :: &
               solzfab(DIMV(solzfab),nparts_def*SDIM)
       real(amrex_real), pointer :: solxfab_ptr(D_DECL(:,:,:),:)
       real(amrex_real), pointer :: solyfab_ptr(D_DECL(:,:,:),:)
       real(amrex_real), pointer :: solzfab_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in),target :: uold(DIMV(uold),AMREX_SPACEDIM)
       real(amrex_real), pointer :: uold_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(inout),target :: unew(DIMV(unew),STATE_NCOMP)
       real(amrex_real), pointer :: unew_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in),target :: &
               lsnew(DIMV(lsnew),num_materials*(SDIM+1))
       real(amrex_real), pointer :: lsnew_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in),target :: den(DIMV(den))
       real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in),target :: mu(DIMV(mu))
       real(amrex_real), pointer :: mu_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in) ::  xlo(SDIM)
       integer, PARAMETER :: nhalf=3
       real(amrex_real) :: xsten(-nhalf:nhalf,SDIM)
       real(amrex_real) :: xpoint(SDIM)
       real(amrex_real), INTENT(in) ::  dx(SDIM)

       integer i,j,k
       integer dir
       integer im
       real(amrex_real) un(AMREX_SPACEDIM)
       real(amrex_real) unp1(AMREX_SPACEDIM)
       real(amrex_real) RCEN
       real(amrex_real) inverseden
       real(amrex_real) mu_cell
       integer vofcomp
       integer dencomp
       real(amrex_real) localF
       real(amrex_real) rho_base
       real(amrex_real) rho_factor
       real(amrex_real) local_temp
       real(amrex_real) DTEMP
       real(amrex_real) T_BASE
       real(amrex_real) V_BASE(SDIM)
       real(amrex_real) cell_density_denom
       real(amrex_real) Fsolid
       real(amrex_real) vt_over_r,ut_over_r
       real(amrex_real) param1,param2,hoop_force_coef
       integer ut_comp
       integer partid,im_solid,partid_crit,im_FSI
       real(amrex_real) LStest,LScrit
       integer im_rigid_CL

       if (FSI_outer_sweeps.eq.0) then
        im_rigid_CL=num_materials
       else if (FSI_outer_sweeps.ge.1) then
        im_rigid_CL=im_elastic_map(FSI_outer_sweeps)+1
       else
        print *,"FSI_outer_sweeps invalid"
        stop
       endif

       if (num_FSI_outer_sweeps.ge.1) then
        !do nothing
       else
        print *,"num_FSI_outer_sweeps invalid: ",num_FSI_outer_sweeps
        stop
       endif
       if ((FSI_outer_sweeps.ge.0).and. &
           (FSI_outer_sweeps.le.num_FSI_outer_sweeps-1)) then
        !do nothing
       else
        print *,"FSI_outer_sweeps invalid: ",FSI_outer_sweeps
        stop
       endif

       if (dt.gt.zero) then
        ! do nothing
       else
        print *,"expecting dt>0: ",dt
        stop
       endif
       if (cur_time_slab.ge.zero) then
        ! do nothing
       else
        print *,"expecting cur_time_slab>=0.0d0: ",cur_time_slab
        stop
       endif

       if (abs(angular_velocity-fort_angular_velocity).le.EPS2) then
        ! do nothing
       else
        print *,"angular_velocity or fort_angular_velocity invalid"
        print *,"angular_velocity=",angular_velocity
        print *,"fort_angular_velocity=",fort_angular_velocity
        stop
       endif

       if ((uncoupled_viscosity.ne.0).and. &
           (uncoupled_viscosity.ne.1)) then
        print *,"uncoupled_viscosity invalid"
        stop
       endif
       if ((update_state.ne.OP_HOOP_BOUSSINESQ_EXPLICIT).and. &
           (update_state.ne.OP_HOOP_BOUSSINESQ_IMPLICIT)) then
        print *,"update_state invalid"
        stop
       endif
       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid hoop implicit"
        stop
       endif
       if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
        print *,"nparts invalid fort_hoopimplicit"
        stop
       endif
       if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
        print *,"nparts_def invalid fort_hoopimplicit"
        stop
       endif

       if (bfact.lt.1) then
        print *,"bfact invalid8"
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
       if (num_state_base.ne.2) then
        print *,"num_state_base invalid"
        stop
       endif

       if (angular_velocity.ge.zero) then
        ! do nothing
       else
        print *,"angular_velocity cannot be negative"
        print *,"expecting counterclockwise"
        print *,"angular_velocity=",angular_velocity
        stop
       endif
       if ((centrifugal_force_factor.le.one).and. &
           (centrifugal_force_factor.ge.zero)) then
        ! do nothing
       else
        print *,"expecting 0<=centrifugal_force_factor<=1"
        stop
       endif

       force_ptr=>force
       call checkbound_array(fablo,fabhi,force_ptr,1,-1)
       tensor_ptr=>tensor
       call checkbound_array(fablo,fabhi,tensor_ptr,0,-1)
       thermal_ptr=>thermal
       call checkbound_array1(fablo,fabhi,thermal_ptr,1,-1)
       recon_ptr=>recon
       call checkbound_array(fablo,fabhi,recon_ptr,1,-1)
       solxfab_ptr=>solxfab
       solyfab_ptr=>solyfab
       solzfab_ptr=>solzfab
       call checkbound_array(fablo,fabhi,solxfab_ptr,0,0)
       call checkbound_array(fablo,fabhi,solyfab_ptr,0,1)
       call checkbound_array(fablo,fabhi,solzfab_ptr,0,SDIM-1)
       uold_ptr=>uold
       call checkbound_array(fablo,fabhi,uold_ptr,1,-1)
       unew_ptr=>unew
       call checkbound_array(fablo,fabhi,unew_ptr,1,-1)
       lsnew_ptr=>lsnew
       call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)

       den_ptr=>den
       call checkbound_array1(fablo,fabhi,den_ptr,1,-1)

       mu_ptr=>mu
       call checkbound_array1(fablo,fabhi,mu_ptr,1,-1)

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do dir=1,SDIM
         xpoint(dir)=xsten(0,dir)
        enddo

        do dir=1,SDIM
         un(dir)=uold(D_DECL(i,j,k),dir)
         unp1(dir)=uold(D_DECL(i,j,k),dir)
        enddo ! dir=1..sdim

        partid=0
        im_solid=0
        im_FSI=0
        partid_crit=0

        do im=1,num_materials
         LStest=lsnew(D_DECL(i,j,k),im)

         if (is_lag_part(im).eq.1) then
          if (is_rigid(im).eq.1) then
           if (LStest.ge.zero) then
            if (im_solid.eq.0) then
             im_solid=im
             partid_crit=partid
             LScrit=LStest
            else if ((im_solid.ge.1).and. &
                     (im_solid.le.num_materials)) then
             if (LStest.ge.LScrit) then
              im_solid=im
              partid_crit=partid
              LScrit=LStest
             endif
            else
             print *,"im_solid invalid 1"
             stop
            endif
           else if (LStest.lt.zero) then
            ! do nothing
           else
            print *,"LStest invalid: ",LStest
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

          if (FSI_outer_sweeps.eq.0) then
           !do nothing
          else if (FSI_outer_sweeps.ge.1) then
           if (LStest.ge.zero) then
            if ((is_rigid_CL(im).eq.1).and. &
                (im.le.im_rigid_CL)) then
             im_FSI=im
            else if (is_rigid_CL(im).eq.0) then
             !do nothing
            else
             print *,"is_rigid_CL(im) invalid: ",im,is_rigid_CL(im)
             stop
            endif 
           else if (LStest.lt.zero) then
            ! do nothing
           else
            print *,"LStest invalid: ",im,LStest
            stop
           endif

          else
           print *,"FSI_outer_sweeps invalid: ",FSI_outer_sweeps
           stop
          endif

          if (is_rigid(im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid invalid DIFFUSION_3D.F90"
           stop
          endif
         else
          print *,"is_lag_part invalid"
          stop
         endif
        enddo ! im=1..num_materials

        if (partid.ne.nparts) then
         print *,"partid invalid"
         stop
        endif

        inverseden=den(D_DECL(i,j,k))
        mu_cell=mu(D_DECL(i,j,k))

        if (inverseden.gt.zero) then
         ! do nothing
        else
         print *,"inverseden invalid: ",inverseden
         stop
        endif

        if (mu_cell.ge.zero) then
         ! do nothing
        else
         print *,"mu_cell invalid: ",mu_cell
         stop
        endif
 
        if ((im_solid.ge.1).and.(im_solid.le.num_materials)) then

         if (im_solid_map(partid_crit+1)+1.ne.im_solid) then
          print *,"im_solid_map(partid_crit+1)+1.ne.im_solid"
          stop
         endif

         dir=1
         unp1(dir)=half*( &
            solxfab(D_DECL(i,j,k),partid_crit*SDIM+dir)+ &
            solxfab(D_DECL(i+1,j,k),partid_crit*SDIM+dir))
         dir=2
         unp1(dir)=half*( &
            solyfab(D_DECL(i,j,k),partid_crit*SDIM+dir)+ &
            solyfab(D_DECL(i,j+1,k),partid_crit*SDIM+dir))
         if (SDIM.eq.3) then
          dir=SDIM
          unp1(dir)=half*( &
            solzfab(D_DECL(i,j,k),partid_crit*SDIM+dir)+ &
            solzfab(D_DECL(i,j,k+1),partid_crit*SDIM+dir))
         endif

        else if ((im_FSI.ge.1).and.(im_FSI.le.num_materials)) then
         !do nothing
        else if ((im_solid.eq.0).and.(im_FSI.eq.0))  then ! in the fluid

         RCEN=xsten(0,1)

         DTEMP=zero
         cell_density_denom=zero
         Fsolid=zero

         do im=1,num_materials

          vofcomp=(im-1)*ngeom_recon+1
          localF=recon(D_DECL(i,j,k),vofcomp)
          if ((localF.ge.-EPS1).and. &
              (localF.le.VOFTOL)) then
           localF=zero
          else if ((localF.ge.one-VOFTOL).and. &
                   (localF.le.one+EPS1)) then
           localF=one
          else if ((localF.gt.zero).and.(localF.lt.one)) then
           ! do nothing
          else
           print *,"localF out of range: ",localF
           stop
          endif
          dencomp=STATECOMP_STATES+(im-1)*num_state_material+1+ENUM_DENVAR
          if (constant_density_all_time(im).eq.1) then
           rho_base=fort_denconst(im)
          else if (constant_density_all_time(im).eq.0) then
           rho_base=unew(D_DECL(i,j,k),dencomp)
          else
           print *,"constant_density_all_time(im) invalid"
           stop
          endif
          if (rho_base.gt.zero) then
           ! do nothing
          else
           print *,"rho_base invalid: ",rho_base
           stop
          endif
          if (is_rigid(im).eq.1) then
           Fsolid=Fsolid+localF
          else if (is_rigid(im).eq.0) then
           cell_density_denom=cell_density_denom+localF*rho_base

           if ((override_density(im).eq.0).or. & ! rho_t + div (rho u) = 0
               (override_density(im).eq.1)) then ! rho=rho(T,Y)
            ! do nothing
            !
            ! NOTE: DrhoDT<=0.0.
            !
            ! if override_density(im)==1,
            ! GODUNOV_3D.F90: fort_derive_mom_den
            ! rho=rho0*(1+fort_DrhoDT(im)*(T-T0))
            ! units of fort_DrhoDT: 1/temperature

            !Boussinesq approximation:
            !rho Du/Dt=-grad p - rho |g| zhat +
            !  rho Omega^2 r rhat
            !rho0 Du/Dt=-grad p - rho |g| zhat +
            !  rho Omega^2 r rhat
            !rho/rho0=(1-|beta|(T-T0))
            !Du/Dt=-grad p/rho0-|g| zhat + |g beta| zhat (T-T0)+
            !  Omega^{2} r rhat -
            !  Omega^{2} r rhat |beta| (T-T0)
            !  (beta<0)

           else if (override_density(im).eq.2) then 

            local_temp=thermal(D_DECL(i,j,k))

            if (local_temp.gt.zero) then
             ! do nothing
            else
             print *,"local_temp cannot be <= 0 fort_hoopimplicit(1)"
             stop
            endif

            if (fort_DrhoDT(im).le.zero) then
             ! do nothing
            else
             print *,"fort_DrhoDT has invalid sign: fort_hoopimplicit"
             print *,"im=",im
             print *,"fort_DrhoDT(im)=",fort_DrhoDT(im)
             stop
            endif

             !STUB_PROCS.F90: T_BASE=fort_tempconst(im)
            call SUB_T0_Boussinesq(xpoint,dx,cur_time_slab,im,T_BASE)

            ! units of DrhoDT are 1/(degrees Kelvin)
            ! DTEMP will have no units after dividing by total density.
            ! fort_tempconst is the temperature of the inner boundary
            ! for the differentially heated rotating annulus problem.
            ! example (default): rho_factor=fort_DrhoDT(im)*(T-Tbase)
            call SUB_UNITLESS_EXPANSION_FACTOR( &
              im, & !intent(in)
              local_temp, & !intent(in)
              T_BASE, & !intent(in)
              rho_factor) !intent(out) (unitless)

             !rho_base=density of material "im"
             !rho_factor=fort_DrhoDT(im)*(T-Tbase)
            DTEMP=DTEMP+localF*rho_base*rho_factor

           else
            print *,"override_density invalid"
            stop
           endif
          else
           print *,"is_rigid(im) invalid"
           stop
          endif
         enddo ! im=1..num_materials

         if (cell_density_denom.gt.zero) then
          ! do nothing
         else
          print *,"cell_density_denom invalid: fort_hoopimplicit"
          print *,"cell_density_denom=",cell_density_denom
          stop
         endif

         if ((Fsolid.ge.zero).and.(Fsolid.le.one+EPS1)) then
          ! do nothing
         else
          print *,"Fsolid invalid: fort_hoopimplicit: ",Fsolid
          stop
         endif

         if (Fsolid.ge.half) then
          DTEMP=zero
         else if (Fsolid.le.half) then
          DTEMP=DTEMP/cell_density_denom !DTEMP will be unitless after this.
         else
          print *,"Fsolid is NaN"
          stop
         endif

          ! gravity force (temperature dependence)
          ! units of gravity: m/s^2
          ! DTEMP has no units.
          ! DTEMP=beta(T-T0)  (beta<0)
          ! usually gravity_vector(SDIM)<0
         if (abs(DTEMP).ge.zero) then
          do dir=1,SDIM
           unp1(dir)=unp1(dir)+ &
              dt*gravity_vector(dir)*DTEMP
          enddo
         else
          print *,"DTEMP is NaN: ",DTEMP
          stop
         endif

         call SUB_V0_Coriolis(xpoint,dx,cur_time_slab,V_BASE)

          ! polar coordinates: coriolis force (temperature dependence)
          !                    centrifugal force (temperature dependence).
          ! angular_velocity>0 => counter clockwise
          ! angular_velocity<0 => clockwise (not allowed)
          ! in PROB.F90:
          ! R-Theta-Z 
          ! pres=pres+half*rho*(angular_velocity**2)*(xpos(1)**2)
          ! X-Y-Z 
          ! pres=pres+half*rho*(angular_velocity**2)*(xpos(1)**2+xpos(2)**2)

         if (rzflag.eq.COORDSYS_CYLINDRICAL) then

          if (RCEN.gt.zero) then
           ! do nothing
          else
           print *,"RCEN invalid: ",RCEN
           stop
          endif
           ! Coriolis "force"
           ! Lewis and Nagata 2004:
           ! -2 Omega e_{z} \Times \vec{u}
          unp1(1)=unp1(1)+dt*(centrifugal_force_factor*(un(2)**2)/RCEN+ &
             two*angular_velocity*(un(2)-V_BASE(2)))
          unp1(2)=unp1(2)-dt*(centrifugal_force_factor*(un(1)*un(2))/RCEN+ &
             two*angular_velocity*(un(1)-V_BASE(1)))

           ! DTEMP has no units.
           ! Lewis and Nagata 2004:
           ! -Omega^{2} r \rhat DrhoDT*(T-Tbase)
           ! DTEMP=beta*(T-T0)  beta<0
          unp1(1)=unp1(1)+ &
           dt*DTEMP*centrifugal_force_factor*(angular_velocity**2)*RCEN

         else if (rzflag.eq.COORDSYS_CARTESIAN) then
          ! assume that RCEN > eps > 0 ?
          ! coriolis force:
          ! -2 omega cross v =
          !  rhat  theta_hat   zhat 
          !  0        0      angular_velocity
          !  u        v         w
          ! = -2(-angular_vel. v,angular_velocity u)
          unp1(1)=unp1(1)+ &
               dt*( two*angular_velocity*(un(2)-V_BASE(2)) )
          unp1(2)=unp1(2)- &
               dt*( two*angular_velocity*(un(1)-V_BASE(1)) )

          if ((DTEMP.eq.zero).or. &
              (angular_velocity.eq.zero)) then
           ! do nothing
          else if ((DTEMP.ne.zero).and. &
                   (angular_velocity.gt.zero)) then
           ! DTEMP has no units.
           ! Lewis and Nagata 2004:
           ! -Omega^{2} r \rhat DrhoDT*(T-Tbase)
           ! DTEMP=beta*(T-T0)  beta<0
           unp1(1)=unp1(1)+ &
            dt*DTEMP*centrifugal_force_factor*(angular_velocity**2)*xsten(0,1)
           unp1(2)=unp1(2)+ &
            dt*DTEMP*centrifugal_force_factor*(angular_velocity**2)*xsten(0,2)
          else
           print *,"DTEMP or angular_velocity invalid"
           print *,"DTEMP=",DTEMP
           print *,"angular_velocity=",angular_velocity
           stop
          endif

         else if (rzflag.eq.COORDSYS_RZ) then

          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif

          if (angular_velocity.eq.zero) then
           ! do nothing
          else
           print *,"angular_velocity<>0 not implemented RZ:", &
            angular_velocity
           stop
          endif

         else
          print *,"rzflag invalid"
          stop
         endif
            
         if (uncoupled_viscosity.eq.0) then
          param1=three
          param2=two
         else if (uncoupled_viscosity.eq.1) then
          param1=two
          param2=one
         else
          print *,"uncoupled_viscosity invalid"
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

          ut_comp=SDIM+1 ! u_theta/r
          ut_over_r=tensor(D_DECL(i,j,k),ut_comp)
          vt_over_r=tensor(D_DECL(i,j,k),ut_comp+1)

           ! units of viscosity: kg/(m s)
           ! units of update:
           ! s kg/(m s) m^3/kg (m/s) (1/m^2)=m/s
          unp1(1)=unp1(1)-param1* &
            dt*visc_coef*mu_cell*inverseden*vt_over_r/RCEN
          unp1(2)=unp1(2)+param1* &
            dt*visc_coef*mu_cell*inverseden*ut_over_r/RCEN

         else
          print *,"rzflag invalid"
          stop
         endif
 
         if (rzflag.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if (rzflag.eq.COORDSYS_RZ) then

          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          if (RCEN.gt.zero) then
           ! do nothing
          else
           print *,"RCEN invalid"
           stop
          endif

           ! units of viscosity: kg/(m s)
           ! units of hoop_force_coef: s kg/(m s) m^3/kg (1/m^2)=1
          hoop_force_coef=dt*visc_coef*mu_cell*inverseden/(RCEN**2)
          if (hoop_force_coef.lt.zero) then
           print *,"hoop_force_coef invalid"
           stop
          endif

          if (update_state.eq.OP_HOOP_BOUSSINESQ_IMPLICIT) then
           unp1(1)=unp1(1)/(one+param2*hoop_force_coef)
          else if (update_state.eq.OP_HOOP_BOUSSINESQ_EXPLICIT) then
           unp1(1)=unp1(1)-param2*hoop_force_coef*un(1)
          else
           print *,"update_state invalid"
           stop
          endif

         else if (rzflag.eq.COORDSYS_CYLINDRICAL) then

          if (RCEN.gt.zero) then
           ! do nothing
          else
           print *,"RCEN invalid"
           stop
          endif

          hoop_force_coef=dt*visc_coef*mu_cell*inverseden/(RCEN**2)
          if (hoop_force_coef.lt.zero) then
           print *,"hoop_force_coef invalid"
           stop
          endif

          if (update_state.eq.OP_HOOP_BOUSSINESQ_IMPLICIT) then
           unp1(1)=unp1(1)/(one+param2*hoop_force_coef)
           unp1(2)=unp1(2)/(one+hoop_force_coef)
          else if (update_state.eq.OP_HOOP_BOUSSINESQ_EXPLICIT) then
           unp1(1)=unp1(1)-param2*hoop_force_coef*un(1)
           unp1(2)=unp1(2)-hoop_force_coef*un(2)
          else
           print *,"update_state invalid"
           stop
          endif

         else
          print *,"rzflag invalid"
          stop
         endif

        else
         print *,"im_solid or im_fsi invalid 2: ",im_FSI,im_solid
         stop
        endif

        if (dt.gt.zero) then
         ! do nothing
        else
         print *,"dt invalid: ",dt
         stop
        endif

        ! viscosity force=-div(2 mu D)-HOOP_FORCE_MARK_MF
        do dir=1,SDIM
         force(D_DECL(i,j,k),dir)=(unp1(dir)-un(dir))/(inverseden*dt)
         if (update_state.eq.OP_HOOP_BOUSSINESQ_EXPLICIT) then
          ! do nothing
         else if (update_state.eq.OP_HOOP_BOUSSINESQ_IMPLICIT) then
          unew(D_DECL(i,j,k),dir)=unp1(dir)
         else
          print *,"update_state invalid"
          stop
         endif

        enddo ! dir=1..sdim

       enddo
       enddo
       enddo

       return
       end subroutine fort_hoopimplicit

       subroutine fort_user_defined_momentum_force( &
         thermal,DIMS(thermal), &
         xlo,dx, &
         uold,DIMS(uold), &
         unew,DIMS(unew), &
         lsnew,DIMS(lsnew), &
         den,DIMS(den), &  ! 1/density
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         finest_level, &
         dt, &
         cur_time_slab, &
         rzflag, &
         nparts, &
         nparts_def, &
         im_solid_map) &
       bind(c,name='fort_user_defined_momentum_force')

       use probf90_module
       use global_utility_module 

       IMPLICIT NONE

       integer, INTENT(in) :: nparts
       integer, INTENT(in) :: nparts_def
       integer, INTENT(in) :: im_solid_map(nparts_def)
       integer, INTENT(in) :: level
       integer, INTENT(in) :: finest_level
       integer, INTENT(in) :: rzflag
       real(amrex_real), INTENT(in) :: dt
       real(amrex_real), INTENT(in) :: cur_time_slab
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, target, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer :: growlo(3),growhi(3)
       integer, INTENT(in) :: bfact
    
       integer, INTENT(in) :: DIMDEC(thermal)
       integer, INTENT(in) :: DIMDEC(uold)
       integer, INTENT(in) :: DIMDEC(unew)
       integer, INTENT(in) :: DIMDEC(lsnew)
       integer, INTENT(in) :: DIMDEC(den)

       real(amrex_real), INTENT(in),target :: thermal(DIMV(thermal))
       real(amrex_real), pointer :: thermal_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in),target :: uold(DIMV(uold),AMREX_SPACEDIM)
       real(amrex_real), pointer :: uold_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(inout),target :: unew(DIMV(unew),STATE_NCOMP)
       real(amrex_real), pointer :: unew_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in),target :: lsnew(DIMV(lsnew),num_materials*(SDIM+1))
       real(amrex_real), pointer :: lsnew_ptr(D_DECL(:,:,:),:)
       real(amrex_real), INTENT(in),target :: den(DIMV(den))
       real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:))
       real(amrex_real), target, INTENT(in) ::  xlo(SDIM)
       integer, PARAMETER :: nhalf=3
       real(amrex_real) :: xsten(-nhalf:nhalf,SDIM)
       real(amrex_real) :: xpoint(SDIM)
       real(amrex_real), target, INTENT(in) :: dx(SDIM)

       integer :: i,j,k
       integer :: dir
       real(amrex_real) :: inverseden
       real(amrex_real) :: output_force(SDIM)

       type(user_defined_force_parm_type_input) :: force_input

       force_input%thermal=>thermal
       force_input%uold=>uold
       force_input%lsnew=>lsnew
       force_input%one_over_den=>den

       force_input%dx=>dx
       force_input%xlo=>xlo
       force_input%fablo=>fablo
       force_input%fabhi=>fabhi
       
       force_input%dt=dt
       force_input%cur_time=cur_time_slab

       if (dt.gt.zero) then
        ! do nothing
       else
        print *,"expecting dt>0"
        stop
       endif
       if (cur_time_slab.ge.zero) then
        ! do nothing
       else
        print *,"expecting cur_time_slab>=0.0d0"
        stop
       endif

       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid hoop implicit"
        stop
       endif
       if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
        print *,"nparts invalid fort_hoopimplicit"
        stop
       endif
       if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
        print *,"nparts_def invalid fort_hoopimplicit"
        stop
       endif

       if (bfact.lt.1) then
        print *,"bfact invalid8"
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
       if (num_state_base.ne.2) then
        print *,"num_state_base invalid"
        stop
       endif

       thermal_ptr=>thermal
       call checkbound_array1(fablo,fabhi,thermal_ptr,1,-1)
       uold_ptr=>uold
       call checkbound_array(fablo,fabhi,uold_ptr,1,-1)
       unew_ptr=>unew
       call checkbound_array(fablo,fabhi,unew_ptr,1,-1)
       lsnew_ptr=>lsnew
       call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)

       den_ptr=>den
       call checkbound_array1(fablo,fabhi,den_ptr,1,-1)

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do dir=1,SDIM
         xpoint(dir)=xsten(0,dir)
        enddo

        inverseden=den(D_DECL(i,j,k))

        if (inverseden.gt.zero) then
         ! do nothing
        else
         print *,"inverseden invalid"
         stop
        endif

        force_input%i=i
        force_input%j=j
        force_input%k=k
        call SUB_USER_DEFINED_FORCE(xpoint,output_force,force_input)
        do dir=1,SDIM
         unew(D_DECL(i,j,k),dir)=uold(D_DECL(i,j,k),dir)+ &
                dt*output_force(dir)*inverseden
        enddo
       enddo
       enddo
       enddo

       return
       end subroutine fort_user_defined_momentum_force

      end module hoop_module

