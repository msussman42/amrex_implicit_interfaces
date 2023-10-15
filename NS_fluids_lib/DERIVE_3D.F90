#undef BL_LANG_CC
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
#include "DERIVE_F.H"
#include "PROB_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

      module derive_module
      use amrex_fort_module, only : amrex_real
      use probcommon_module

      implicit none

      contains


       ! WALE model
       ! Wall Adapting Local Eddy-viscosity models for simulations in complex
       ! geometries.   Ducros, Nicoud, Poinsot 1998
       ! called from getStateVISC
      subroutine fort_derturbvisc( &
       les_model, &
       level, &
       im, &  ! im=1..num_materials
       dt, &
       denstate,DIMS(denstate), &
       vof,DIMS(vof), &
       vel,DIMS(vel), &
       visc,DIMS(visc), &
       cellten,DIMS(cellten), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       cur_time, &
       dx, &
       xlo, &
       dx_coarsest, &
       ngrow, &
       ncompvisc) &
      bind(c,name='fort_derturbvisc')

      use global_utility_module
      use probf90_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: les_model
      integer, INTENT(in) :: level
      integer, INTENT(in) :: im ! 1..num_materials
      integer, INTENT(in) :: ncompvisc
      integer, INTENT(in) :: ngrow
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: cur_time

      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx_coarsest(SDIM)

      integer, INTENT(in) :: DIMDEC(denstate)
      integer, INTENT(in) :: DIMDEC(vof)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(visc)
      integer, INTENT(in) :: DIMDEC(cellten)
  
      real(amrex_real), INTENT(in), target :: &
              denstate(DIMV(denstate),num_materials*num_state_material)
      real(amrex_real), pointer :: denstate_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: vof(DIMV(vof),num_materials*ngeom_recon)
      real(amrex_real), pointer :: vof_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: vel(DIMV(vel),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: visc(DIMV(visc),ncompvisc)
      real(amrex_real), pointer :: visc_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: cellten(DIMV(cellten),AMREX_SPACEDIM_SQR)
      real(amrex_real), pointer :: cellten_ptr(D_DECL(:,:,:),:)

      real(amrex_real) g(3,3),s(3,3),sd(3,3),g2(3),g2tr,ss,sdsd
      real(amrex_real) g2_12,g2_13,g2_21,g2_23,g2_31,g2_32,turb_visc,density
      real(amrex_real) rr
      integer i,j,k
      integer i1,j1,k1
      integer k1lo,k1hi
      integer veldir,dir
      integer flagcomp
      integer vofcomp
      integer nbase
      real(amrex_real) xsten(-1:1,SDIM)
      integer nhalf
      integer near_solid
      integer im_primary
      integer im_local
      real(amrex_real) VFRAC(num_materials)
     
      nhalf=1

      if (ngrow.eq.1) then
       ! do nothing
      else
       print *,"expecting ngrow==1"
       stop
      endif

      if ((ncompvisc.ne.num_materials).and.(ncompvisc.ne.3*num_materials)) then
       print *,"ncompvisc invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid1"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid1"
       stop
      endif
      if (dt.ge.zero) then 
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
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (level.lt.0) then
       print *,"level invalid 30 "
       stop
      endif

      vof_ptr=>vof 
      call checkbound_array(fablo,fabhi,vof_ptr,ngrow+1,-1)
      visc_ptr=>visc
      call checkbound_array(fablo,fabhi,visc_ptr,ngrow,-1)
      vel_ptr=>vel
      call checkbound_array(fablo,fabhi,vel_ptr,ngrow,-1)
      denstate_ptr=>denstate
      call checkbound_array(fablo,fabhi,denstate_ptr,ngrow,-1)
      cellten_ptr=>cellten
      call checkbound_array(fablo,fabhi,cellten_ptr,ngrow,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      if (is_rigid(im).eq.1) then
       print *,"cannot have eddy wall viscosity and (is_rigid(im).eq.1)"
       stop
      else if (is_rigid(im).eq.0) then
       ! do nothing
      else
       print *,"is_rigid invalid DERIVE_3D.F90 1"
       stop
      endif

       ! density
      flagcomp=(im-1)*num_state_material+1+ENUM_DENVAR

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       do im_local=1,num_materials
        vofcomp=(im_local-1)*ngeom_recon+1
        VFRAC(im_local)=vof(D_DECL(i,j,k),vofcomp)
       enddo
       call get_primary_material_VFRAC(VFRAC,im_primary)

       do veldir=1,3
       do dir=1,3
        g(veldir,dir)=zero
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
         print *,"dir invalid get shear"
         stop
        endif
        do veldir=1,SDIM
         g(veldir,dir)=cellten(D_DECL(i,j,k),nbase+veldir) 
        enddo
       enddo ! dir

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        rr=xsten(0,1)
        g(3,3)=vel(D_DECL(i,j,k),1)/abs(rr)
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        rr=xsten(0,1)
        g(2,2)=g(2,2)+vel(D_DECL(i,j,k),1)/abs(rr)
        g(1,2)=g(1,2)-vel(D_DECL(i,j,k),2)/abs(rr)
       else
        print *,"levelrz invalid getturbvisc"
        stop
       endif
 
       ! Build the strain tensor
       do veldir=1,3
       do dir=1,3
        s(veldir,dir)=half*(g(veldir,dir)+g(dir,veldir))
       enddo
       enddo

        ! ss=s : s = s_{ij}s_{ij}
       ss=s(1,1)*s(1,1)+s(1,2)*s(1,2)+s(1,3)*s(1,3)+&
          s(2,1)*s(2,1)+s(2,2)*s(2,2)+s(2,3)*s(2,3)+&
          s(3,1)*s(3,1)+s(3,2)*s(3,2)+s(3,3)*s(3,3)

       ! Builds the traceless symmetric part of the square of the
       ! velocity tensor
       ! g2 (the vector)=diagonal(g^T g)
       ! g2_{i}=g_{ik}g_{ki}=g^{2}_{ii} (no summation on last)
       do veldir=1,3
        g2(veldir)=zero
        do dir=1,3
         g2(veldir)=g2(veldir)+g(veldir,dir)*g(dir,veldir)
        enddo
       enddo
        ! g2tr=g^{2}_{kk}/3
       g2tr=(g2(1)+g2(2)+g2(3))/3.0d0 ! tensor trace

        ! the off diagonals of g2 (the matrix) are 
        ! g^{2}_{ij} \equiv g_{ik}g_{kj}\equiv (G * G)_{ij}
       g2_12=g(1,1)*g(1,2)+g(1,2)*g(2,2)+g(1,3)*g(3,2)
       g2_13=g(1,1)*g(1,3)+g(1,2)*g(2,3)+g(1,3)*g(3,3)
       g2_21=g(2,1)*g(1,1)+g(2,2)*g(2,1)+g(2,3)*g(3,1)
       g2_23=g(2,1)*g(1,3)+g(2,2)*g(2,3)+g(2,3)*g(3,3)
       g2_31=g(3,1)*g(1,1)+g(3,2)*g(2,1)+g(3,3)*g(3,1)
       g2_32=g(3,1)*g(1,2)+g(3,2)*g(2,2)+g(3,3)*g(3,2)

        ! S^{d}_{ij}=(1/2)(g_{ij}^{2}+g_{ji}^{2})-
        !  (1/3)delta_{ij}g_{kk}^{2} 
       sd(1,1)=g2(1)-g2tr  !g^T g - (1/3)Tr(g^T g)
       sd(1,2)=0.5d0*(g2_12+g2_21) !(g g + (g g)^T)/2
       sd(1,3)=0.5d0*(g2_13+g2_31)

       sd(2,1)=sd(1,2)
       sd(2,2)=g2(2)-g2tr
       sd(2,3)=0.5*(g2_23+g2_32)

       sd(3,1)=sd(1,3)
       sd(3,2)=sd(2,3)
       sd(3,3)=g2(3)-g2tr

       ! tensor squared
       ! sdsd=S^{d}_{ij}S^{d}_{ij}
       sdsd=sd(1,1)*sd(1,1)+sd(1,2)*sd(1,2)+sd(1,3)*sd(1,3)+&
            sd(2,1)*sd(2,1)+sd(2,2)*sd(2,2)+sd(2,3)*sd(2,3)+&
            sd(3,1)*sd(3,1)+sd(3,2)*sd(3,2)+sd(3,3)*sd(3,3)

       if (les_model.eq.1) then

         ! units of viscosity: Pa s = Newton sec/m^2 = (kg m/sec^2)(sec/m^2) =
         !  kg/(m sec) 
         ! units of g: (1/sec)
         ! units of s: (1/sec)
         ! units of ss: (1/sec^2)
         ! units of g2: (1/sec^2)
         ! units of sd: (1/sec^2)
         ! units of sdsd: (1/sec^4)
         ! units of sdsd^(3/2): (1/sec^6)
         ! units of sdsd^(5/4): (1/sec^5)
         ! units of ss^(5/2)=(1/sec^5)
         ! units of turb_visc: m^2 (1/sec^6)/(1/sec^5)=m^2/s
         ! units of turb_visc * density=kg/(m s)
        if(sdsd.gt.zero) then
         turb_visc=((0.5d0*dx_coarsest(1))**2)*(sdsd**1.5d0)/ &
                 (ss**2.5d0+sdsd**1.25d0)
         density=denstate(D_DECL(i,j,k),flagcomp)
         ! limiter: arbitrary at this point
         turb_visc=min(turb_visc, &
           two*(visc(D_DECL(i,j,k),im)+fort_viscconst_eddy_bulk(im))/density)

         if (im_primary.eq.im) then
          if (VFRAC(im).ge.one-VOFTOL) then
           visc(D_DECL(i,j,k),im)=visc(D_DECL(i,j,k),im)+turb_visc*density
          else if ((VFRAC(im).ge.-VOFTOL).and. &
                   (VFRAC(im).le.one-VOFTOL)) then
           ! do nothing
          else
           print *,"VFRAC invalid"
           stop
          endif
         else if ((im_primary.ne.im).and. &
                  (im_primary.ge.1).and. &
                  (im_primary.le.num_materials)) then
          ! do nothing
         else
          print *,"im_primary invalid"
          stop
         endif
        else if (sdsd.eq.zero) then
         turb_visc=zero
        else
         print *,"sdsd nan"
         stop
        endif

       else if (les_model.eq.0) then
        ! do nothing
       else
        print *,"les_model invalid"
        stop
       endif

       if (fort_viscconst_eddy_bulk(im).gt.zero) then
        visc(D_DECL(i,j,k),im)=visc(D_DECL(i,j,k),im)+ &
               fort_viscconst_eddy_bulk(im)
       else if (fort_viscconst_eddy_bulk(im).eq.zero) then
        ! do nothing
       else
        print *,"fort_viscconst_eddy_bulk invalid"
        stop
       endif

       if (fort_viscconst_eddy_wall(im).gt.zero) then
        k1lo=0
        k1hi=0
        if (SDIM.eq.3) then
         k1lo=-1
         k1hi=1
        else if (SDIM.eq.2) then
         ! do nothing
        else
         print *,"dimension problem"
         stop
        endif
        near_solid=0
        do i1=-1,1
        do j1=-1,1
        do k1=k1lo,k1hi
         do im_local=1,num_materials
          vofcomp=(im_local-1)*ngeom_recon+1
          VFRAC(im_local)=vof(D_DECL(i+i1,j+j1,k+k1),vofcomp)
         enddo
         call get_primary_material_VFRAC(VFRAC,im_primary)
         if (is_rigid(im_primary).eq.1) then
          near_solid=1
         else if (is_rigid(im_primary).eq.0) then
          ! do nothing
         else
          print *,"is_rigid invalid DERIVE_3D.F90 2"
          stop
         endif
        enddo !k1
        enddo !j1
        enddo !i1

        if (near_solid.eq.1) then
         visc(D_DECL(i,j,k),im)=visc(D_DECL(i,j,k),im)+ &
                 fort_viscconst_eddy_wall(im)
        else if (near_solid.eq.0) then
         ! do nothing
        else
         print *,"near_solid invalid"
         stop
        endif

       else if (fort_viscconst_eddy_wall(im).eq.zero) then
        ! do nothing
       else
        print *,"fort_viscconst_eddy_wall invalid"
        stop
       endif
  
      enddo 
      enddo 
      enddo 

      return
      end subroutine fort_derturbvisc

      subroutine fort_getshear( &
       cellten,DIMS(cellten), &
       vel,DIMS(vel), &
       dx,xlo, &
       tensordata, &
       DIMS(tensordata), &
       only_scalar, &
       time, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       bc, &
       ngrow) &
      bind(c,name='fort_getshear')

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real) xsten(-1:1,SDIM)
      integer nhalf
      integer i,j,k,n,dir,veldir
      integer, INTENT(in) :: ngrow,only_scalar
      integer i1,j1
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: bc(SDIM,2,SDIM)
      integer, INTENT(in) :: DIMDEC(cellten)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(tensordata)
  
      real(amrex_real), INTENT(in), target :: cellten(DIMV(cellten),AMREX_SPACEDIM_SQR)
      real(amrex_real), pointer :: cellten_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: vel(DIMV(vel),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(out), target :: &
           tensordata(DIMV(tensordata),DERIVE_TENSOR_NCOMP)
      real(amrex_real), pointer :: tensordata_ptr(D_DECL(:,:,:),:)

      real(amrex_real) visctensor(3,3),gradu(3,3)
      real(amrex_real) shear
      real(amrex_real) rr
      real(amrex_real) vort(3)
      integer nbase

      nhalf=1

      tensordata_ptr=>tensordata
      cellten_ptr=>cellten
      vel_ptr=>vel

      if (bfact.lt.1) then
       print *,"bfact invalid2"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      call checkbound_array(fablo,fabhi,cellten_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,vel_ptr,ngrow+1,-1)
      call checkbound_array(fablo,fabhi,tensordata_ptr,ngrow,-1)

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
       print *,"levelrz invalid getshear"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       do veldir=1,3
       do dir=1,3
        gradu(veldir,dir)=zero
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
         print *,"dir invalid get shear"
         stop
        endif
        do veldir=1,SDIM
         gradu(veldir,dir)=cellten(D_DECL(i,j,k),nbase+veldir) 
        enddo
       enddo ! dir=1..sdim

! grad u=| u_r  u_t/r-v/r  u_z  |
!        | v_r  v_t/r+u/r  v_z  |
!        | w_r  w_t/r      w_z  |

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        rr=abs(xsten(0,1))
        if (rr.gt.zero) then
         gradu(3,3)=vel(D_DECL(i,j,k),1)/rr
        else
         print *,"rr invalid"
         stop
        endif
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        rr=abs(xsten(0,1))
        if (rr.gt.zero) then
         gradu(2,2)=gradu(2,2)+vel(D_DECL(i,j,k),1)/rr
         gradu(1,2)=gradu(1,2)-vel(D_DECL(i,j,k),2)/rr
        else
         print *,"rr invalid"
         stop
        endif
       else
        print *,"levelrz invalid getshear 2"
        stop
       endif

       do i1=1,3
       do j1=1,3
        visctensor(i1,j1)=half*(gradu(i1,j1)+gradu(j1,i1))
       enddo
       enddo

#if (AMREX_SPACEDIM==3)
       shear=visctensor(1,1)**2+visctensor(2,2)**2+ &
          visctensor(3,3)**2+two*(visctensor(1,2)**2)+ &
          two*(visctensor(1,3)**2)+two*(visctensor(2,3)**2)
#elif (AMREX_SPACEDIM==2)
       shear=visctensor(1,1)**2+visctensor(2,2)**2+ &
          visctensor(3,3)**2+two*(visctensor(1,2)**2)
#else
       print *,"dimension bust"
       stop
#endif
       shear=sqrt(two*shear)

       if (only_scalar.eq.1) then
        tensordata(D_DECL(i,j,k),DERIVE_TENSOR_MAG+1)=shear
       else if (only_scalar.eq.2) then
        vort(1)=gradu(3,2)-gradu(2,3)
        vort(2)=gradu(1,3)-gradu(3,1)
        vort(3)=gradu(2,1)-gradu(1,2)
        tensordata(D_DECL(i,j,k),DERIVE_TENSOR_MAG+1)= &
          sqrt(vort(1)**2+vort(2)**2+vort(3)**2)
       else if (only_scalar.eq.0) then
        tensordata(D_DECL(i,j,k),DERIVE_TENSOR_MAG+1)=shear
        n=DERIVE_TENSOR_RATE_DEFORM+1
        do i1=1,3
        do j1=1,3
         tensordata(D_DECL(i,j,k),n)=visctensor(i1,j1)
         n=n+1
        enddo
        enddo
        if (n.eq.DERIVE_TENSOR_GRAD_VEL+1) then
         ! do nothing
        else
         print *,"n.eq.DERIVE_TENSOR_GRAD_VEL+1 failed"
         stop
        endif

        do i1=1,3
        do j1=1,3
         tensordata(D_DECL(i,j,k),n)=gradu(i1,j1) !gradu(veldir,dir)
         n=n+1
        enddo
        enddo

        if (n.eq.DERIVE_TENSOR_NCOMP+1) then
         ! do nothing
        else
         print *,"n.eq.DERIVE_TENSOR_NCOMP+1 failed"
         stop
        endif
       else 
        print *,"only_scalar invalid"
        stop
       endif

      enddo
      enddo
      enddo

      return
      end subroutine fort_getshear

       ! called from NavierStokes::getStateVISC() (NavierStokes2.cpp) 
      subroutine fort_derviscosity( &
        level, &
        finest_level, &
        visc_coef, &
        im_parm, & !1<=im_parm<=num_materials
        dt, & ! used for the viscoelastic coefficient.
        viscosity_coefficient, & ! viscconst(im_parm)
        shear_thinning_fluid, &
        Carreau_alpha, &
        Carreau_beta, &
        Carreau_n, &
        Carreau_mu_inf, &
        concentration, &
        elastic_time, &
        viscosity_state_model, &
        viscoelastic_model, &
        elastic_viscosity, &
        etaL,etaP,etaS, &
        polymer_factor, &
        visc,DIMS(visc), &
        vel,DIMS(vel), &
        eosdata,DIMS(eosdata), &
        tensor,DIMS(tensor), &
        gammadot,DIMS(gammadot), &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        time, &
        dx,xlo, &
        bc,ngrow, &
        ncompvisc) &
      bind(c,name='fort_derviscosity')

      use global_utility_module
      use probf90_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: visc_coef
      integer, INTENT(in) :: im_parm !1<=im_parm<=num_materials
      integer, INTENT(in) :: ncompvisc
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: viscosity_coefficient
      integer, INTENT(in) :: shear_thinning_fluid
      real(amrex_real), INTENT(in) :: Carreau_alpha
      real(amrex_real), INTENT(in) :: Carreau_beta
      real(amrex_real), INTENT(in) :: Carreau_n
      real(amrex_real), INTENT(in) :: Carreau_mu_inf
      real(amrex_real), INTENT(in) :: concentration,etaL,etaP,etaS
      real(amrex_real), INTENT(in) :: elastic_time
      real(amrex_real), INTENT(in) :: elastic_viscosity
      integer, INTENT(in) :: viscosity_state_model
      integer, INTENT(in) :: viscoelastic_model
      real(amrex_real), INTENT(in) :: polymer_factor
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: DIMDEC(visc)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(eosdata)
      integer, INTENT(in) :: DIMDEC(tensor)
      integer, INTENT(in) :: DIMDEC(gammadot)
      integer, INTENT(in) :: bc(SDIM,2,SDIM)
      integer, INTENT(in) :: ngrow
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM)
       !ncompvisc=3*num_materials
      real(amrex_real), INTENT(out), target :: visc(DIMV(visc),ncompvisc) 
      real(amrex_real), pointer :: visc_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: vel(DIMV(vel),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: gammadot(DIMV(gammadot))
      real(amrex_real), pointer :: gammadot_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: eosdata(DIMV(eosdata), &
              num_materials*num_state_material)
      real(amrex_real), pointer :: eosdata_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: tensor(DIMV(tensor),ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), pointer :: tensor_ptr(D_DECL(:,:,:),:)

      integer i,j,k
      real(amrex_real)    shear,density,temperature,mu
      integer flagcomp
      real(amrex_real) bterm,pterm
      integer numstatetest,ii,jj
      integer dir_local
      real(amrex_real) Q(3,3)
      real(amrex_real) traceA,modtime,viscoelastic_coeff

      visc_ptr=>visc

      if (bfact.lt.1) then
       print *,"bfact invalid3"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 31"
       stop
      endif

      if ((im_parm.lt.1).or.(im_parm.gt.num_materials)) then
       print *,"im_parm invalid3"
       stop
      endif

      if (shear_thinning_fluid.eq.fort_shear_thinning_fluid(im_parm)) then
       ! do nothing
      else
       print *,"shear_thinning_fluid invalid"
       stop
      endif

      if (Carreau_alpha.eq.fort_Carreau_alpha(im_parm)) then
       ! do nothing
      else
       print *,"fort_Carreau_alpha(im_parm) invalid"
       stop
      endif

      if (Carreau_beta.eq.fort_Carreau_beta(im_parm)) then
       ! do nothing
      else
       print *,"fort_Carreau_beta(im_parm) invalid"
       stop
      endif

      if (Carreau_n.eq.fort_Carreau_n(im_parm)) then
       ! do nothing
      else
       print *,"fort_Carreau_n(im_parm) invalid"
       stop
      endif

      if (Carreau_mu_inf.eq.fort_Carreau_mu_inf(im_parm)) then
       ! do nothing
      else
       print *,"fort_Carreau_mu_inf(im_parm) invalid"
       stop
      endif

      if (concentration.eq.fort_concentration(im_parm)) then
       ! do nothing
      else
       print *,"fort_concentration(im_parm) invalid"
       stop
      endif
      if (etaL.eq.fort_etaL(im_parm)) then
       ! do nothing
      else
       print *,"fort_etaL(im_parm) invalid"
       stop
      endif

      if (etaP.eq.fort_etaP(im_parm)) then
       ! do nothing
      else
       print *,"fort_etaP(im_parm) invalid"
       stop
      endif

      if (etaS.eq.fort_etaS(im_parm)) then
       ! do nothing
      else
       print *,"fort_etaS(im_parm) invalid"
       stop
      endif

      if (polymer_factor.eq.fort_polymer_factor(im_parm)) then
       ! do nothing
      else
       print *,"fort_polymer_factor(im_parm) invalid"
       stop
      endif

      if (visc_coef.eq.fort_visc_coef) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif

      if (elastic_viscosity.eq.fort_elastic_viscosity(im_parm)) then
       ! do nothing
      else
       print *,"fort_elastic_viscosity(im_parm) invalid"
       stop
      endif
      if (elastic_time.eq.fort_elastic_time(im_parm)) then
       ! do nothing
      else
       print *,"fort_elastic_time(im_parm) invalid"
       stop
      endif
      if (viscoelastic_model.eq.fort_viscoelastic_model(im_parm)) then
       ! do nothing
      else
       print *,"fort_viscoelastic_model(im_parm) invalid"
       stop
      endif

      if (viscosity_state_model.ne. &
          fort_viscosity_state_model(im_parm)) then
       print *,"viscosity_state_model invalid"
       stop
      endif
      if (viscosity_state_model.ge.0) then
       ! do nothing
      else
       print *,"viscosity_state_model invalid"
       stop
      endif

      if (ENUM_NUM_TENSOR_TYPE.eq.2*AMREX_SPACEDIM) then
       ! do nothing
      else
       print *,"expecting ENUM_NUM_TENSOR_TYPE.eq.2*AMREX_SPACEDIM"
       stop
      endif

      if (fort_built_in_elastic_model(elastic_viscosity, &
              viscoelastic_model).eq.1) then
       ! do nothing
      else if (fort_built_in_elastic_model(elastic_viscosity, &
                 viscoelastic_model).eq.0) then
       ! do nothing
      else
       print *,"fort_built_in_elastic_model invalid"
       stop
      endif

      if (dt.gt.zero) then 
       ! do nothing
      else
       print *,"dt invalid in fort_derviscosity"
       stop
      endif

      if (polymer_factor.ge.zero) then
       ! do nothing
      else
       print *,"polymer_factor invalid"
       stop
      endif
      if ((shear_thinning_fluid.ne.0).and. &
          (shear_thinning_fluid.ne.1)) then
       print *,"shear_thinning_fluid invalid"
       stop
      endif
      if (viscosity_coefficient.ge.zero) then
       ! do nothing
      else
       print *,"viscosity_coefficient invalid"
       stop
      endif
      if (abs(viscosity_coefficient-fort_viscconst(im_parm)).le.1.0E-14) then
       ! do nothing
      else
       print *,"viscosity_coefficient invalid"
       stop
      endif

      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif

      if (ncompvisc.ne.3*num_materials) then
       print *,"ncompvisc invalid"
       stop
      endif

      call checkbound_array(fablo,fabhi,visc_ptr,ngrow,-1)
      gammadot_ptr=>gammadot
      call checkbound_array1(fablo,fabhi,gammadot_ptr,ngrow,-1)
      eosdata_ptr=>eosdata
      call checkbound_array(fablo,fabhi,eosdata_ptr,ngrow,-1)
      tensor_ptr=>tensor
      call checkbound_array(fablo,fabhi,tensor_ptr,ngrow,-1)
      vel_ptr=>vel
      call checkbound_array(fablo,fabhi,vel_ptr,ngrow,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      if (is_rigid(im_parm).eq.1) then
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        visc(D_DECL(i,j,k),im_parm)=viscosity_coefficient 
       enddo
       enddo
       enddo
      else if (is_rigid(im_parm).eq.0) then

       vel_ptr=>vel
       call checkbound_array(fablo,fabhi,vel_ptr,ngrow+1,-1)

       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

         ! den,T
        if (viscosity_state_model.ge.1) then
         flagcomp=(im_parm-1)*num_state_material+1+ENUM_DENVAR
         density=eosdata(D_DECL(i,j,k),flagcomp)
         temperature=eosdata(D_DECL(i,j,k),flagcomp+1)
        else if (viscosity_state_model.eq.0) then
         density=fort_denconst(im_parm)
         temperature=fort_tempconst(im_parm)
        else
         print *,"viscosity_state_model invalid"
         stop
        endif

         !viscconst(im_parm) in default case.
        mu=get_user_viscconst(im_parm,density,temperature)

        if ((viscoelastic_model.eq.0).or. & !FENE-CR
            (viscoelastic_model.eq.1).or. & ! Oldroyd B
            (viscoelastic_model.eq.5).or. & ! FENE-P
            (Viscoelastic_model.eq.6)) then ! linear PTT
         ! do nothing
        else if (viscoelastic_model.eq.3) then !incremental
         ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
         if (elastic_viscosity.ge.zero) then
          ! do nothing
         else
          print *,"elastic_viscosity invalid"
          stop
         endif
        else if (viscoelastic_model.eq.7) then !incremental
         ! Xia, Lu, Tryggvason 2018
         if (elastic_viscosity.ge.zero) then
          ! do nothing
         else
          print *,"elastic_viscosity invalid"
          stop
         endif
        else
         print *,"viscoelastic_model invalid: ",viscoelastic_model
         stop
        endif

        visc(D_DECL(i,j,k),im_parm) = mu

        if (shear_thinning_fluid.eq.1) then
         shear=gammadot(D_DECL(i,j,k))
         if ((im_parm.eq.1).and.(Carreau_beta.eq.zero).and. &
             (probtype.eq.2).and.(axis_dir.gt.0)) then
          call viscosity(axis_dir,visc(D_DECL(i,j,k),im_parm),shear)
         else if (Carreau_beta.gt.zero) then
          if (Carreau_n.gt.one) then
           print *,"Carreau_n invalid"
           stop
          endif
          if (Carreau_mu_inf.lt.zero) then
           print *,"Carreau_mu_inf invalid"
          endif
          pterm=(Carreau_n-one)/Carreau_alpha
          bterm=one+(Carreau_beta*shear)**Carreau_alpha
          if ((elastic_time.gt.zero).and. &
              (elastic_viscosity.gt.zero)) then
           visc(D_DECL(i,j,k),im_parm)=etaS+etaP*(bterm**pterm)
          else if ((elastic_time.eq.zero).and. &
                   (elastic_viscosity.eq.zero)) then
           visc(D_DECL(i,j,k),im_parm)=Carreau_mu_inf+ &
              (etaL-Carreau_mu_inf)*(bterm**pterm)
          else
           print *,"elastic_time or elastic_viscosity invalid"
           stop
          endif
         else
          print *,"carreau beta invalid"
          stop
         endif
        else if (shear_thinning_fluid.eq.0) then
         ! do nothing
        else
         print *,"shear_thinning_fluid invalid"
         stop
        endif

        if ((elastic_time.gt.zero).and. &
            (elastic_viscosity.gt.zero)) then

          if (ncompvisc.ne.3*num_materials) then
           print *,"ncompvisc invalid"
           stop
          endif

          numstatetest=num_state_base+num_species_var
          if (num_state_material.ne.numstatetest) then
           print *,"num_state_material invalid"
           stop
          endif
          if ((num_materials_viscoelastic.lt.1).or. &
              (num_materials_viscoelastic.gt.num_materials)) then
           print *,"num_materials_viscoelastic invalid:fort_derviscosity"
           stop
          endif

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

          if (SDIM.eq.3) then
           ! do nothing
          else if (SDIM.eq.2) then
           if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
               (levelrz.eq.COORDSYS_CYLINDRICAL)) then
            ! do nothing
           else if (levelrz.eq.COORDSYS_RZ) then
            ! do nothing
           else
            print *,"levelrz invalid"
            stop
           endif
          else
           print *,"dimension bust"
           stop
          endif

           ! A^{n+1}=(I-dt grad U)^T A^{adv} (I - dt grad U) =
           !  W^T A^{adv} W; a similarity transformation, in which
           !  W is non-singular, preserves the eigenvalues of A^{adv}.
           !  initially A is the identity matrix.
           ! For elastic materials, Q=grad X + (grad X)^{T}
           ! and sigma_ij= lambda delta_ij trace(grad X) +
           !               mu Q_ij
           ! lambda=Lame's first parameter.
           ! for elastic materials, Q does not need to be SPD.
          traceA=zero

          if ((viscoelastic_model.eq.0).or. & ! FENE-CR
              (viscoelastic_model.eq.1).or. & ! Oldroyd B
              (viscoelastic_model.eq.5).or. & ! FENE-P
              (viscoelastic_model.eq.6)) then ! linear PTT

           traceA=zero
           do ii=1,3
            traceA=traceA+Q(ii,ii)+one
            if (Q(ii,ii)+one.gt.zero) then
             ! do nothing
            else
             print *,"diagonal entries of SPD matrix cannot be negative"
             print *,"ii,Q(ii,ii)= ",ii,Q(ii,ii)
             print *,"level,finest_level ",level,finest_level
             print *,"im_parm,ngrow,i,j,k = ",im_parm,ngrow,i,j,k
             print *,"dt= ",dt
             stop
            endif
           enddo ! ii=1..3

           if ((viscoelastic_model.eq.0).or. & !FENE-CR
               (viscoelastic_model.eq.5)) then ! FENE-P

            ! declared in PROB.F90
            ! modtime=max(0.0,elastic_time*(1-Tr(A)/L^2))=max(0,lambda/f(A))
            ! polymer_factor=1/L
            call get_mod_elastic_time(elastic_time,traceA, &
             polymer_factor,modtime)

           else if (viscoelastic_model.eq.1) then ! Oldroyd B

            modtime=elastic_time

           else if (viscoelastic_model.eq.6) then ! linear PTT

            modtime=elastic_time

           else
            print *,"viscoelastic_model invalid"
            stop
           endif

           ! modtime=elastic_time >> 1 for 2,3,7
          else if ((viscoelastic_model.eq.3).or. & !incremental (plastic)
                   (viscoelastic_model.eq.7)) then !incremental (Neo-Hookean)
           modtime=elastic_time
          else
           print *,"viscoelastic_model invalid: ",viscoelastic_model
           stop
          endif

          if (modtime.ge.zero) then
           ! do nothing
          else
           print *,"modtime invalid"
           stop
          endif
         
           ! etaS=etaL-etaP=viscconst-elastic_viscosity 
          if (modtime+dt.gt.zero) then
           if ((viscoelastic_model.eq.0).or. &!FENE-CR(lambda_tilde=f(A)/lambda)
               (viscoelastic_model.eq.1).or. &!Oldroyd-B(modtime=elastic_time)
               (viscoelastic_model.eq.5).or. &!FENE-P (lambda_tilde=f(A)/lambda
               (viscoelastic_model.eq.6)) then!linearPTT(modtime=elastic_time)
             ! etaS=etaL-etaP=viscconst-elastic_viscosity 
            viscoelastic_coeff= &
             (visc(D_DECL(i,j,k),im_parm)-etaS)/(modtime+dt)
           else if (viscoelastic_model.eq.3) then !incremental
            ! Maire et al
            viscoelastic_coeff=elastic_viscosity
           else if (viscoelastic_model.eq.7) then !incremental
            ! Xia, Lu, Tryggvason
            viscoelastic_coeff=elastic_viscosity
           else
            print *,"viscoelastic_model invalid: ",viscoelastic_model
            stop
           endif
          else
           print *,"expecting modtime+dt to be positive"
           print *,"modtime=",modtime
           print *,"dt=",dt
           stop
          endif

          if (viscoelastic_coeff.ge.zero) then
           ! do nothing
          else
           print *,"viscoelastic_coef invalid"
           stop
          endif

          visc(D_DECL(i,j,k),num_materials+im_parm)= &
             viscoelastic_coeff*visc_coef
          visc(D_DECL(i,j,k),2*num_materials+im_parm)=modtime

        else if ((elastic_time.eq.zero).and. &
                 (elastic_viscosity.eq.zero)) then
          ! do nothing
        else
         print *,"elastic_time/elastic_viscosity invalid"
         stop
        endif 

       enddo
       enddo
       enddo

      else
       print *,"im_parm invalid4"
       stop
      endif

      return
      end subroutine fort_derviscosity


       ! called from getStateCONDUCTIVITY
      subroutine fort_derconductivity( &
       NS_sumdata_size, &
       NS_sumdata, &
       ncomp_sum_int_user1, &
       ncomp_sum_int_user2, &
       level, &
       finest_level, &
       im_parm, & ! =1..num_materials
       dt, &
       conduct,DIMS(conduct), &
       eosdata,DIMS(eosdata), &
       vof,DIMS(vof), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       time, &
       dx,xlo, &
       ngrow) &
      bind(c,name='fort_derconductivity')

      use global_utility_module
      use probf90_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: ncomp_sum_int_user1
      integer, INTENT(in) :: ncomp_sum_int_user2
      integer :: ncomp_sum_int_user12
      integer, INTENT(in) :: NS_sumdata_size
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: im_parm !=1..num_materials
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: DIMDEC(conduct)
      integer, INTENT(in) :: DIMDEC(eosdata)
      integer, INTENT(in) :: DIMDEC(vof)
      integer, INTENT(in) :: ngrow
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM)

      real(amrex_real), INTENT(in) :: NS_sumdata(NS_sumdata_size)

      real(amrex_real), INTENT(out), target :: conduct(DIMV(conduct),num_materials) 
      real(amrex_real), pointer :: conduct_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: eosdata(DIMV(eosdata), &
              num_materials*num_state_material)
      real(amrex_real), pointer :: eosdata_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: vof(DIMV(vof),num_materials*ngeom_recon)
      real(amrex_real), pointer :: vof_ptr(D_DECL(:,:,:),:)

      integer i,j,k
      integer ii,jj,kk
      integer isolid,jsolid,ksolid
      integer iprobe,jprobe,kprobe
      integer dir,side
      integer dir_local
      integer flagcomp
      integer vofcomp
      real(amrex_real) density
      real(amrex_real) temperature
      real(amrex_real) temperature_wall
      real(amrex_real) temperature_wall_max
      real(amrex_real) temperature_probe
      real(amrex_real) thermal_k
      real(amrex_real) thermal_k_max
      real(amrex_real) xvec(SDIM)
      real(amrex_real) xsten(-1:1,SDIM)
      integer nhalf
      integer im_local
      integer im_solid_crit
      integer im_primary_center
      integer im_primary_side
      integer im_primary_probe
      real(amrex_real) VFRAC(num_materials)
      integer near_interface
      real(amrex_real) nrm(SDIM)

      ncomp_sum_int_user12=ncomp_sum_int_user1+ncomp_sum_int_user2
      if (NS_sumdata_size.ne.IQ_TOTAL_SUM_COMP) then
       print *,"mismatch between NS_sumdata_size and IQ_TOTAL_SUM_COMP"
       stop
      endif

      nhalf=1

      conduct_ptr=>conduct

      if (bfact.lt.1) then
       print *,"bfact invalid3"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 31"
       stop
      endif

      if ((im_parm.lt.1).or.(im_parm.gt.num_materials)) then
       print *,"im_parm invalid3"
       stop
      endif

      if (dt.gt.zero) then 
       ! do nothing
      else
       print *,"dt invalid in fort_derconductivity"
       stop
      endif
      if (time.ge.zero) then 
       ! do nothing
      else
       print *,"time invalid in fort_derconductivity"
       stop
      endif
      if (fort_heatviscconst_eddy_bulk(im_parm).ge.zero) then
       ! do nothing
      else
       print *,"fort_heatviscconst_eddy_bulk invalid"
       stop
      endif
      if (fort_heatviscconst_eddy_wall(im_parm).ge.zero) then
       ! do nothing
      else
       print *,"fort_heatviscconst_eddy_wall invalid"
       stop
      endif

      call checkbound_array(fablo,fabhi,conduct_ptr,ngrow,-1)
      eosdata_ptr=>eosdata
      call checkbound_array(fablo,fabhi,eosdata_ptr,ngrow+2,-1)

      vof_ptr=>vof 
      call checkbound_array(fablo,fabhi,vof_ptr,ngrow+2,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir=1,SDIM
        xvec(dir)=xsten(0,dir)
       enddo

       flagcomp=(im_parm-1)*num_state_material+1+ENUM_DENVAR
       density=eosdata(D_DECL(i,j,k),flagcomp)
       temperature=eosdata(D_DECL(i,j,k),flagcomp+1)
       temperature_wall=temperature
       temperature_wall_max=temperature
       temperature_probe=temperature

       thermal_k=get_user_heatviscconst(im_parm)

       near_interface=0
       do dir_local=1,SDIM
        nrm(dir_local)=zero
       enddo
       im_solid_crit=0

       if (is_in_probtype_list().eq.1) then
        call SUB_THERMAL_K( &
          xvec,dx,time, &
          density, &
          temperature, &
          thermal_k, & ! INTENT(inout)
          im_parm, &
          near_interface, &
          im_solid_crit, &
          temperature_wall, &
          temperature_wall_max, &
          temperature_probe, &
          nrm) ! nrm points from solid to fluid
       endif

       thermal_k_max=thermal_k

       do im_local=1,num_materials
        vofcomp=(im_local-1)*ngeom_recon+1
        VFRAC(im_local)=vof(D_DECL(i,j,k),vofcomp)
       enddo
       call get_primary_material_VFRAC(VFRAC,im_primary_center)

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
        do side=-1,1,2
         isolid=i
         jsolid=j
         ksolid=k
         iprobe=i 
         jprobe=j 
         kprobe=k 
         im_solid_crit=0
         do im_local=1,num_materials
          vofcomp=(im_local-1)*ngeom_recon+1
          VFRAC(im_local)=vof(D_DECL(i+side*ii,j+side*jj,k+side*kk),vofcomp)
         enddo
         call get_primary_material_VFRAC(VFRAC,im_primary_side)

         near_interface=0
         do dir_local=1,SDIM
          nrm(dir_local)=zero
         enddo

         if ((is_rigid(im_primary_side).eq.1).and. &
             (is_rigid(im_primary_center).eq.1)) then
          ! do nothing 
         else if ((is_rigid(im_primary_side).eq.1).and. &
                  (is_rigid(im_primary_center).eq.0).and. &
                  (im_primary_center.eq.im_parm)) then
          im_solid_crit=im_primary_side
          near_interface=1
          isolid=i+side*ii
          jsolid=j+side*jj
          ksolid=k+side*kk
          iprobe=isolid-2*side*ii
          jprobe=jsolid-2*side*jj
          kprobe=ksolid-2*side*kk
          nrm(dir)=-side  ! points from solid to fluid
         else if ((is_rigid(im_primary_center).eq.1).and. &
                  (is_rigid(im_primary_side).eq.0).and. &
                  (im_primary_side.eq.im_parm)) then
          im_solid_crit=im_primary_center
          near_interface=1
          isolid=i
          jsolid=j
          ksolid=k
          iprobe=isolid+2*side*ii
          jprobe=jsolid+2*side*jj
          kprobe=ksolid+2*side*kk
          nrm(dir)=side ! points from solid to fluid
         else if ((is_rigid(im_primary_center).eq.0).and. &
                  (is_rigid(im_primary_side).eq.0)) then
          ! do nothing
         else if ((is_rigid(im_primary_side).eq.1).and. &
                  (is_rigid(im_primary_center).eq.0).and. &
                  (im_primary_center.ne.im_parm)) then
          ! do nothing
         else if ((is_rigid(im_primary_center).eq.1).and. &
                  (is_rigid(im_primary_side).eq.0).and. &
                  (im_primary_side.ne.im_parm)) then
          ! do nothing
         else
          print *,"is_rigid invalid DERIVE_3D.F90 3"
          stop
         endif
         if (near_interface.eq.1) then
          do im_local=1,num_materials
           vofcomp=(im_local-1)*ngeom_recon+1
           VFRAC(im_local)=vof(D_DECL(iprobe,jprobe,kprobe),vofcomp)
          enddo
          call get_primary_material_VFRAC(VFRAC,im_primary_probe)
          if ((is_rigid(im_primary_probe).eq.1).or. &
              (im_primary_probe.ne.im_parm)) then
           near_interface=0
          else if ((is_rigid(im_primary_probe).eq.0).and. &
                   (im_primary_probe.eq.im_parm)) then
           ! do nothing
          else
           print *,"is_rigid or im_primary_probe invalid"
           stop
          endif
         else if (near_interface.eq.0) then
          ! do nothing
         else
          print *,"near_interface invalid"
          stop
         endif

         if (near_interface.eq.1) then
          flagcomp=(im_solid_crit-1)*num_state_material+1+ENUM_DENVAR
          temperature_wall=eosdata(D_DECL(isolid,jsolid,ksolid),flagcomp+1)
          temperature_wall_max= &
            NS_sumdata(IQ_MAXSTATE_SUM_COMP+2*(im_solid_crit-1)+2)
          flagcomp=(im_parm-1)*num_state_material+1+ENUM_DENVAR
          temperature_probe=eosdata(D_DECL(iprobe,jprobe,kprobe),flagcomp+1)
          thermal_k=get_user_heatviscconst(im_parm)+ &
               fort_heatviscconst_eddy_wall(im_parm)
         else if (near_interface.eq.0) then
          temperature_wall=temperature
          temperature_wall_max=temperature
          temperature_probe=temperature
          thermal_k=get_user_heatviscconst(im_parm)
         else
          print *,"near_interface invalid"
          stop
         endif

         if (is_in_probtype_list().eq.1) then
          call SUB_THERMAL_K( &
            xvec,dx,time, &
            density, &
            temperature, &
            thermal_k, & ! INTENT(inout)
            im_parm, &
            near_interface, &
            im_solid_crit, &
            temperature_wall, &
            temperature_wall_max, &
            temperature_probe, &
            nrm) ! nrm points from solid to fluid
         endif

         if (thermal_k.ge.zero) then
          ! do nothing
         else
          print *,"thermal_k became corrupt"
          stop
         endif
         if (thermal_k.ge.thermal_k_max) then
          thermal_k_max=thermal_k
         else if (thermal_k.le.thermal_k_max) then
          ! do nothing
         else
          print *,"thermal_k or thermal_k_max corrupt"
          stop
         endif
        enddo ! side=-1,1,2
       enddo ! dir=1,SDIM

       if (1.eq.0) then
        print *,"i,j,k,thermal_k,thermal_k_max,im_parm ", &
            i,j,k,thermal_k,thermal_k_max,im_parm
       endif

       conduct(D_DECL(i,j,k),im_parm) = thermal_k_max+ &
          fort_heatviscconst_eddy_bulk(im_parm)

      enddo
      enddo
      enddo

      return
      end subroutine fort_derconductivity

! if FENE-CR+Carreau,
! liquid viscosity=etaS+etaP ( 1+ (beta gamma_dot)^alpha )^((n-1)/alpha)
!
! for each material, there are 5 components:
! 1. \dot{gamma}
! 2. Tr(A) if viscoelastic
!    \dot{gamma} o.t.
! 3. Tr(A) (liquid viscosity - etaS)/etaP  if FENE-CR+Carreau
!    Tr(A) if FENE-CR
!    \dot{gamma} o.t.
! 4. (3) * f(A)  if viscoelastic
!    \dot{gamma} o.t.
! 5. magnitude of vorticity

      subroutine fort_dermagtrace( &
       level, &
       finest_level, &
       im, &   ! im=0..num_materials-1
       cellten,DIMS(cellten), &
       dest,DIMS(dest), &
       den,DIMS(den), &
       tensor,DIMS(tensor), &
       vel,DIMS(vel), &
       visc,DIMS(visc), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       ngrow, &
       dx,xlo,time, &
       bc, &
       ncomp_den, &
       ncomp_visc, &
       n_trace, &
       polymer_factor, &
       etaS, &
       etaP, &
       Carreau_beta, &
       elastic_time, &
       viscoelastic_model, &
       elastic_viscosity) &
      bind(c,name='fort_dermagtrace')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: ncomp_den
      integer, INTENT(in) :: ncomp_visc
      integer, INTENT(in) :: n_trace

      real(amrex_real), INTENT(in) :: polymer_factor(num_materials)
      real(amrex_real), INTENT(in) :: etaS(num_materials)
      real(amrex_real), INTENT(in) :: etaP(num_materials)
      real(amrex_real), INTENT(in) :: Carreau_beta(num_materials)
      real(amrex_real), INTENT(in) :: elastic_time(num_materials)
      integer, INTENT(in) :: viscoelastic_model(num_materials)
      real(amrex_real), INTENT(in) :: elastic_viscosity(num_materials)

      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: DIMDEC(dest)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(tensor)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(visc)
      integer, INTENT(in) :: DIMDEC(cellten)
      integer, INTENT(in) :: bc(SDIM,2,SDIM)
      integer, INTENT(in) :: ngrow
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dx(SDIM), xlo(SDIM)
      real(amrex_real), INTENT(in), target :: cellten(DIMV(cellten),AMREX_SPACEDIM_SQR)

      real(amrex_real), INTENT(inout), target :: dest(DIMV(dest),5)
      real(amrex_real), pointer :: dest_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: den(DIMV(den),ncomp_den)
      real(amrex_real), INTENT(in), target :: tensor(DIMV(tensor),ENUM_NUM_TENSOR_TYPE)
      real(amrex_real), INTENT(in), target :: vel(DIMV(vel),STATE_NCOMP_VEL)
      real(amrex_real), INTENT(in), target :: visc(DIMV(visc),ncomp_visc)

      integer im  ! im=0..num_materials-1
      integer i,j,k
      real(amrex_real) T11,T22,T33,traceA,modtime
      integer nbase
      integer dir,veldir
      real(amrex_real) gradu(3,3)
      real(amrex_real) xsten(-1:1,SDIM)
      integer nhalf
      real(amrex_real) rr
      real(amrex_real) vort(3)

      nhalf=1

      dest_ptr=>dest

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid dermagtrace"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid4"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((im.lt.0).or.(im.ge.num_materials)) then
       print *,"im invalid5"
       stop
      endif
      if (n_trace.ne.num_materials*5) then
       print *,"n_trace invalid"
       stop
      endif
      if (ncomp_den.ne.num_materials*num_state_material) then
       print *,"ncomp_den invalid"
       stop
      endif
      if (ncomp_visc.lt.num_materials) then
       print *,"ncomp_visc invalid"
       stop
      endif  
      if (fort_built_in_elastic_model(elastic_viscosity(im+1), &
             viscoelastic_model(im+1)).eq.1) then 
       ! do nothing
      else if (fort_built_in_elastic_model(elastic_viscosity(im+1), &
             viscoelastic_model(im+1)).eq.0) then 
       ! do nothing
      else
       print *,"fort_built_in_elastic_model invalid"
       stop
      endif
      if (fort_elastic_viscosity(im+1).ne.elastic_viscosity(im+1)) then
       print *,"elastic_viscosity(im+1) invalid"
       stop
      endif
      if (fort_built_in_elastic_model(elastic_viscosity(im+1), &
             fort_viscoelastic_model(im+1)).eq.1) then 
       if (num_materials_viscoelastic.le.0) then
        print *,"num_materials_viscoelastic.le.0:fort_dermagtrace"
        stop
       endif
      else if (fort_built_in_elastic_model(elastic_viscosity(im+1), &
                 fort_viscoelastic_model(im+1)).eq.0) then 
       ! do nothing
      else
       print *,"fort_built_in_elastic_model invalid"
       stop
      endif  

      call checkbound_array(fablo,fabhi,cellten,ngrow,-1)
      call checkbound_array(fablo,fabhi,dest_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,den,ngrow,-1)
      call checkbound_array(fablo,fabhi,tensor,ngrow,-1)
      call checkbound_array(fablo,fabhi,vel,ngrow+1,-1)
      call checkbound_array(fablo,fabhi,visc,ngrow,-1)

       ! in: fort_dermagtrace
       ! visc=sqrt(2*(a11**2+a22**2+a33**2+2*a12**2+2*a13**2+2*a23**2))
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 
   
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
 
       call gridsten_level(xsten,i,j,k,level,nhalf)

       do veldir=1,3
       do dir=1,3
        gradu(veldir,dir)=zero
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
         print *,"dir invalid get shear"
         stop
        endif
        do veldir=1,SDIM
         gradu(veldir,dir)=cellten(D_DECL(i,j,k),nbase+veldir) 
        enddo
       enddo ! dir=1..sdim

! in RZ:
! grad u=| u_r  0  u_z  |
!        | 0   u/r  0   |
!        | w_r  0  w_z  |
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        rr=xsten(0,1)
        gradu(3,3)=vel(D_DECL(i,j,k),1)/abs(rr)
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        rr=xsten(0,1)
        gradu(2,2)=gradu(2,2)+vel(D_DECL(i,j,k),1)/abs(rr)
        gradu(1,2)=gradu(1,2)-vel(D_DECL(i,j,k),2)/abs(rr)
       else
        print *,"levelrz invalid dermagtrace 2"
        stop
       endif

       vort(1)=gradu(3,2)-gradu(2,3)
       vort(2)=gradu(1,3)-gradu(3,1)
       vort(3)=gradu(2,1)-gradu(1,2)
       dest(D_DECL(i,j,k),5)= &
         sqrt(vort(1)**2+vort(2)**2+vort(3)**2)

        ! dest(1) calculated in fort_getshear.
       if (fort_built_in_elastic_model(elastic_viscosity(im+1), &
             fort_viscoelastic_model(im+1)).eq.0) then 
        dest(D_DECL(i,j,k),2)=dest(D_DECL(i,j,k),1)
        dest(D_DECL(i,j,k),3)=dest(D_DECL(i,j,k),1)
        dest(D_DECL(i,j,k),4)=dest(D_DECL(i,j,k),1)
       else if (fort_built_in_elastic_model(elastic_viscosity(im+1), &
                 fort_viscoelastic_model(im+1)).eq.1) then 
        if (ncomp_visc.ne.3*num_materials) then
         print *,"ncomp_visc invalid"
         stop
        endif

        T11=tensor(D_DECL(i,j,k),1)+one
        T22=tensor(D_DECL(i,j,k),3)+one
        T33=tensor(D_DECL(i,j,k),4)+one
        traceA=T11+T22+T33

        if ((T11.gt.zero).and. &
            (T22.gt.zero).and. &
            (T33.gt.zero)) then
         ! do nothing
        else if ((T11.le.zero).or. &
                 (T22.le.zero).or. &
                 (T33.le.zero)) then
         if ((fort_viscoelastic_model(im+1).eq.0).or. &!FENE-CR
             (fort_viscoelastic_model(im+1).eq.1).or. &!Oldroyd-B 
             (fort_viscoelastic_model(im+1).eq.5).or. &!FENE-P 
             (fort_viscoelastic_model(im+1).eq.7).or. &!incremental,neohookean
             (fort_viscoelastic_model(im+1).eq.6)) then!linear PTT
          print *,"T11, T22, T33 must be positive"
          stop
         else if ((fort_viscoelastic_model(im+1).eq.3).or. & !incremental,plst
                  (fort_viscoelastic_model(im+1).eq.4)) then !pres vel coupling
          ! check nothing
         else
          print *,"fort_viscoelastic_model(im+1) invalid"
          stop
         endif
        else
         print *,"T11,T22, or T33 is NaN"
         stop
        endif
        dest(D_DECL(i,j,k),2)=traceA

        modtime=visc(D_DECL(i,j,k),2*num_materials+im+1)

        if (elastic_time(im+1).eq.zero) then
         modtime=zero
        else if (elastic_time(im+1).gt.zero) then
         modtime=modtime/elastic_time(im+1)
        else
         print *,"elastic time invalid"
         stop
        endif

        if (is_rigid(im+1).eq.1) then
         dest(D_DECL(i,j,k),3)=traceA
         dest(D_DECL(i,j,k),4)=zero
        else if (Carreau_beta(im+1).eq.zero) then
         dest(D_DECL(i,j,k),3)=traceA
         dest(D_DECL(i,j,k),4)=traceA*modtime
        else if (Carreau_beta(im+1).gt.zero) then
         if (etaP(im+1).le.zero) then
          print *,"etaP invalid"
          stop
         endif
         dest(D_DECL(i,j,k),3)=traceA* &
          (visc(D_DECL(i,j,k),im+1)-etaS(im+1))/etaP(im+1)
         dest(D_DECL(i,j,k),4)= &
          dest(D_DECL(i,j,k),3)*modtime
        else
         print *,"carreau beta invalid"
         stop
        endif

       else
        print *,"fort_built_in_elastic_model invalid"
        stop
       endif

      enddo
      enddo
      enddo

      return
      end subroutine fort_dermagtrace

       ! decompose: 2 (mu0-mu0+mu(grad U) D = 2 mu0 D + 2(mu(grad U)-mu0)D
       ! NavierStokes::GetDrag is called from
       !  NavierStokes::volWgtSumALL
       ! fort_getdrag is called from NavierStokes::GetDrag
       ! see DRAG_COMP.H
      subroutine fort_getdrag( &
       tid, &
       level, &
       finest_level, &
       isweep, &
       globalsum, &
       globalsum_sweep, &
       localsum, &
       tdata,DIMS(tdata), &  ! grad U (CELLTENSOR_MF)
       viscoten,DIMS(viscoten), &  ! viscoelastic configuration tensor
       den,DIMS(den), & 
       mask,DIMS(mask), & 
       slrecon,DIMS(slrecon), &
       levelpc,DIMS(levelpc), &
       vol,DIMS(vol), &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       cvisc,DIMS(cvisc), &
       c_mat_visc, &
       DIMS(c_mat_visc), &
       xlo,dx, &
       solxfab,DIMS(solxfab), &
       solyfab,DIMS(solyfab), &
       solzfab,DIMS(solzfab), &
       pres,DIMS(pres), &
       vel,DIMS(vel), &
       drag,DIMS(drag), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       rzflag,bc, &
       time, &
       visc_coef, &
       nparts, &
       nparts_def, &
       im_solid_map) &
      bind(c,name='fort_getdrag')

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probf90_module
      use godunov_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level

      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def)

      integer, INTENT(in) :: isweep
      real(amrex_real), INTENT(in) :: globalsum(N_DRAG_IQ)
      integer, INTENT(in) :: globalsum_sweep(N_DRAG_IQ)
      real(amrex_real), INTENT(inout) :: localsum(N_DRAG_IQ)
      real(amrex_real), INTENT(in) :: visc_coef
      integer, INTENT(in) :: rzflag
      real(amrex_real), INTENT(in) :: time
      real(amrex_real) presmag
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: bc(SDIM,2,SDIM)
      integer, INTENT(in) :: DIMDEC(tdata)
      integer, INTENT(in) :: DIMDEC(viscoten)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: DIMDEC(slrecon)
      integer, INTENT(in) :: DIMDEC(levelpc)
      integer, INTENT(in) :: DIMDEC(vol)
      integer, INTENT(in) :: DIMDEC(areax)
      integer, INTENT(in) :: DIMDEC(areay)
      integer, INTENT(in) :: DIMDEC(areaz)

      integer, INTENT(in) :: DIMDEC(xface)
      integer, INTENT(in) :: DIMDEC(yface)
      integer, INTENT(in) :: DIMDEC(zface)

      integer, INTENT(in) :: DIMDEC(cvisc)
      integer, INTENT(in) :: DIMDEC(c_mat_visc)

      integer, INTENT(in) :: DIMDEC(solxfab)
      integer, INTENT(in) :: DIMDEC(solyfab)
      integer, INTENT(in) :: DIMDEC(solzfab)

      integer, INTENT(in) :: DIMDEC(pres)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(drag)

      real(amrex_real), INTENT(in),target :: tdata(DIMV(tdata),AMREX_SPACEDIM_SQR)
      real(amrex_real), pointer :: tdata_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: viscoten(DIMV(viscoten),NUM_CELL_ELASTIC)
      real(amrex_real), pointer :: viscoten_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: den(DIMV(den),num_materials*num_state_material)
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: slrecon(DIMV(slrecon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: slrecon_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: levelpc(DIMV(levelpc),num_materials*(SDIM+1))
      real(amrex_real), pointer :: levelpc_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: vol(DIMV(vol))
      real(amrex_real), pointer :: vol_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: areax(DIMV(areax))
      real(amrex_real), pointer :: areax_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: areay(DIMV(areay))
      real(amrex_real), pointer :: areay_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: areaz(DIMV(areaz))
      real(amrex_real), pointer :: areaz_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in),target :: xface(DIMV(xface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: xface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: yface(DIMV(yface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: yface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: zface(DIMV(zface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: zface_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in),target :: cvisc(DIMV(cvisc))
      real(amrex_real), pointer :: cvisc_ptr(D_DECL(:,:,:))

        ! c_mat_visc initialized in fort_derviscosity
        ! 1..num_materials viscosity
        ! num_materials+1 ... 2*num_materials viscoelastic coeff.
        ! 2*num_materials+1 ... 3*num_materials modtime
      real(amrex_real), INTENT(in),target :: c_mat_visc(DIMV(c_mat_visc),3*num_materials)
      real(amrex_real), pointer :: c_mat_visc_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in),target :: solxfab(DIMV(solxfab),nparts_def*SDIM)
      real(amrex_real), pointer :: solxfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: solyfab(DIMV(solyfab),nparts_def*SDIM)
      real(amrex_real), pointer :: solyfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: solzfab(DIMV(solzfab),nparts_def*SDIM)
      real(amrex_real), pointer :: solzfab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in),target :: pres(DIMV(pres))
      real(amrex_real), pointer :: pres_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: vel(DIMV(vel),STATE_NCOMP_VEL)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: drag(DIMV(drag),N_DRAG)
      real(amrex_real), pointer :: drag_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      integer ii,jj,kk
      integer imac,jmac,kmac
      integer i_face,j_face,k_face
      integer i_side,j_side,k_side
      integer i_side_visc,j_side_visc,k_side_visc
      integer ii_visc,jj_visc,kk_visc
      integer dir
      integer veldir
      integer facedir
      integer side_cell  ! 0 or 1
      integer side_visc
      integer dirend
      integer dir_visc
      integer i1,j1
      integer vofcomp,dencomp,viscbase
      integer icell,jcell,kcell
      real(amrex_real) vel6point(SDIM,2,SDIM)
      real(amrex_real) ls_visc(num_materials)
      real(amrex_real) ls_side(num_materials)
      real(amrex_real) pressure_load(3)
      real(amrex_real) viscous_stress_load(3)
      real(amrex_real) viscous0_stress_load(3)
      real(amrex_real) visco_stress_load(3)

      real(amrex_real) pressure_stress_tensor(3,3)
      real(amrex_real) viscous_stress_tensor(3,3)
      real(amrex_real) viscous0_stress_tensor(3,3)
      real(amrex_real) visco_stress_tensor(3,3)

      real(amrex_real) xsten(-3:3,SDIM)
      real(amrex_real) xsten_face(-3:3,SDIM)
      integer nhalf

      real(amrex_real) volume
      real(amrex_real) gradu(3,3)

      real(amrex_real) facearea
      integer im
      integer im_visc
      integer im_test
      integer im_primary
      integer im_side
      real(amrex_real) volgrid,mass,local_density
      real(amrex_real) cengrid(SDIM)
      real(amrex_real) global_centroid(SDIM)
      real(amrex_real) rvec(3)
      real(amrex_real) gravvector(3)
      real(amrex_real) rcross(3)
      real(amrex_real) pressure_rcross(3)
      real(amrex_real) viscous0_rcross(3)
      real(amrex_real) viscous_rcross(3)
      real(amrex_real) visco_rcross(3)
      real(amrex_real) grav_localtorque(SDIM)
      real(amrex_real) pressure_localtorque(SDIM)
      real(amrex_real) viscous_localtorque(SDIM)
      real(amrex_real) viscous0_localtorque(SDIM)
      real(amrex_real) visco_localtorque(SDIM)
      real(amrex_real) Q(3,3)
      real(amrex_real) nsolid(3)
      integer mask_cell
      real(amrex_real) ls_sort(num_materials)
      real(amrex_real) mu_0 ! viscosity for ambient case.
      real(amrex_real) mu_non_ambient ! viscosity for non ambient case.
      real(amrex_real) delx
      integer partid
      integer ibase
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdata_tess(num_materials*ngeom_recon)
      integer local_tessellate
      integer nmax

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid fort_getdrag"
       stop
      endif

      nhalf=3
      nmax=POLYGON_LIST_MAX  ! in: fort_getdrag

      tdata_ptr=>tdata
      viscoten_ptr=>viscoten
      den_ptr=>den
      mask_ptr=>mask
      slrecon_ptr=>slrecon
      levelpc_ptr=>levelpc
      vol_ptr=>vol
      areax_ptr=>areax
      areay_ptr=>areay
      areaz_ptr=>areaz
      xface_ptr=>xface
      yface_ptr=>yface
      zface_ptr=>zface

      cvisc_ptr=>cvisc
      c_mat_visc_ptr=>c_mat_visc

      solxfab_ptr=>solxfab
      solyfab_ptr=>solyfab
      solzfab_ptr=>solzfab
      pres_ptr=>pres
      vel_ptr=>vel
      drag_ptr=>drag

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in fort_getdrag"
      else if (time.lt.zero) then
       print *,"time invalid in fort_getdrag"
       stop
      else
       print *,"time bust in fort_getdrag"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_getdrag"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_getdrag"
       stop
      endif

        ! cell centered grad U
      call checkbound_array(fablo,fabhi,tdata_ptr,0,-1)
      call checkbound_array(fablo,fabhi,viscoten_ptr,1,-1)
      call checkbound_array(fablo,fabhi,den_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      call checkbound_array(fablo,fabhi,slrecon_ptr,1,-1)
      call checkbound_array(fablo,fabhi,levelpc_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,vol_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,areax_ptr,0,0)
      call checkbound_array1(fablo,fabhi,areay_ptr,0,1)
      call checkbound_array1(fablo,fabhi,areaz_ptr,0,SDIM-1)

      call checkbound_array(fablo,fabhi,xface_ptr,0,0)
      call checkbound_array(fablo,fabhi,yface_ptr,0,1)
      call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1)

      call checkbound_array1(fablo,fabhi,cvisc_ptr,0,-1)
      call checkbound_array(fablo,fabhi,c_mat_visc_ptr,1,-1)

      call checkbound_array(fablo,fabhi,solxfab_ptr,0,0)
      call checkbound_array(fablo,fabhi,solyfab_ptr,0,1)
      call checkbound_array(fablo,fabhi,solzfab_ptr,0,SDIM-1)

      call checkbound_array1(fablo,fabhi,pres_ptr,1,-1)
      call checkbound_array(fablo,fabhi,vel_ptr,1,-1)
       ! DRAG_MF has ngrow_make_distance=3 ghost cells
      call checkbound_array(fablo,fabhi,drag_ptr,3,-1)

      if (bfact.lt.1) then
       print *,"bfact invalid5"
       stop
      endif
      if (visc_coef.eq.fort_visc_coef) then
       ! do nothing
      else
       print *,"fort_visc_coef invalid"
       stop
      endif
      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif
      do im=1,num_materials
       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"denconst invalid"
        stop
       endif
       mu_0=get_user_viscconst(im,fort_denconst(im),fort_tempconst(im))
       if (mu_0.ge.zero) then
        ! do nothing
       else
        print *,"mu_0=get_user_viscconst invalid"
        stop
       endif
      enddo
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

      if ((isweep.lt.0).or.(isweep.gt.1)) then
       print *,"isweep invalid"
       stop
      endif
      if (ENUM_NUM_TENSOR_TYPE.ne.2*SDIM) then
       print *,"ENUM_NUM_TENSOR_TYPE invalid"
       stop
      endif
      if ((num_materials_viscoelastic.ge.1).and. &
          (num_materials_viscoelastic.le.num_materials)) then
       if (ENUM_NUM_TENSOR_TYPE*num_materials_viscoelastic.ne. &
           NUM_CELL_ELASTIC) then
        print *,"NUM_CELL_ELASTIC invalid1"
        stop
       endif
      else if (num_materials_viscoelastic.eq.0) then
       ! do nothing
      else
       print *,"num_materials_viscoelastic invalid: fort_getdrag"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       ! first sweep - find the mass and centroid of materials
      if (isweep.eq.0) then
       ! do nothing
      else if (isweep.eq.1) then 
        ! im_test is the material on which a force/torque is applied.
       do im_test=1,num_materials
        mass=globalsum(DRAGCOMP_IQ_MASS+im_test)
        if (mass.lt.zero) then
         print *,"mass cannot be negative  im_test,mass=",im_test,mass
         stop
        else if (mass.eq.zero) then
         print *,"WARNING: mass=0"
         print *,"im_test=",im_test
         print *,"FSI_flag(im_test)=",FSI_flag(im_test)
         do im=1,num_materials
          print *,"im,mass ",im,globalsum(DRAGCOMP_IQ_MASS+im)
         enddo
         do dir=1,SDIM
          global_centroid(dir)=zero
         enddo
        else if (mass.gt.zero) then
         do dir=1,SDIM
          global_centroid(dir)= &
           globalsum(DRAGCOMP_IQ_COM+3*(im_test-1)+dir)/mass
         enddo
         if (levelrz.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if (levelrz.eq.COORDSYS_RZ) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          global_centroid(1)=zero
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
          ! do nothing
         else
          print *,"levelrz invalid fort_getdrag"
          stop
         endif
        else
         print *,"mass=NaN"
         stop
        endif
       enddo ! im_test=1..num_materials
      else
       print *,"isweep invalid"
       stop
      endif

      do icell=growlo(1),growhi(1)
      do jcell=growlo(2),growhi(2)
      do kcell=growlo(3),growhi(3)

       do dir=1,num_materials*ngeom_recon
        mofdata(dir)=slrecon(D_DECL(icell,jcell,kcell),dir)
        mofdata_tess(dir)=mofdata(dir)
       enddo
       call gridsten(xsten,xlo,icell,jcell,kcell,fablo,bfact,dx,nhalf)
       call Box_volumeFAST(bfact,dx,xsten,nhalf,volgrid,cengrid,SDIM)
       if (volgrid.gt.zero) then
        ! do nothing
       else
        print *,"volgrid invalid"
        stop
       endif

        ! before (mofdata): fluids tessellate
        ! after  (mofdata): fluids and solids tessellate
       local_tessellate=1
       call multi_get_volume_tessellate( &
         local_tessellate, &
         bfact, &
         dx,xsten,nhalf, &
         mofdata_tess, &
         geom_xtetlist(1,1,1,tid+1), &
         nmax, &
         nmax, &
         SDIM)

       ! first sweep - find the mass and centroid of materials
       if (isweep.eq.0) then
 
        do dir=1,N_DRAG
         drag(D_DECL(icell,jcell,kcell),dir)=zero
        enddo

        mask_cell=NINT(mask(D_DECL(icell,jcell,kcell)))
        if (mask_cell.eq.1) then

         ! im_test is the material on which a force/torque is applied.
         do im_test=1,num_materials
          vofcomp=(im_test-1)*ngeom_recon+1
          dencomp=(im_test-1)*num_state_material+1+ENUM_DENVAR
          local_density=den(D_DECL(icell,jcell,kcell),dencomp)
          if (local_density.gt.zero) then
           ! do nothing
          else
           print *,"local_density invalid"
           stop
          endif
          mass=local_density*volgrid*mofdata_tess(vofcomp)
          if (globalsum_sweep(DRAGCOMP_IQ_MASS+im_test).eq.0) then
           localsum(DRAGCOMP_IQ_MASS+im_test)= &
             localsum(DRAGCOMP_IQ_MASS+im_test)+mass
          else
           print *,"globalsum_sweep invalid"
           stop
          endif
          do dir=1,SDIM
           if (globalsum_sweep(DRAGCOMP_IQ_COM+3*(im_test-1)+dir).eq.0) then
            localsum(DRAGCOMP_IQ_COM+3*(im_test-1)+dir)= &
             localsum(DRAGCOMP_IQ_COM+3*(im_test-1)+dir)+ &
              mass*(mofdata_tess(vofcomp+dir)+cengrid(dir))
           else
            print *,"globalsum_sweep invalid"
            stop
           endif
          enddo ! dir=1..sdim

         enddo ! im_test=1..num_materials

        else if (mask_cell.eq.0) then
         ! do nothing
        else
         print *,"mask_cell invalid"
         stop
        endif 

       else if (isweep.eq.1) then ! above, mass and centroid

        mask_cell=NINT(mask(D_DECL(icell,jcell,kcell)))
        if (mask_cell.eq.1) then

         ! force=integral body forces + integral_boundary tau dot n dA
         !  tau=-pI + 2 mu D + mu_p f(A)/lambda  \tilde{Q} 
         ! buoyancy force (body forces within the materials).
         ! Also, update the moment of inertia integral.

         ! DRAGCOMP_FLAG not changed when accumulating the body forces
         ! within materials.
 
         ! calculate the forces exerted on material "im_test"
         ! i.e. im_test is the material on which a force/torque is applied.
         do im_test=1,num_materials
          vofcomp=(im_test-1)*ngeom_recon+1
          dencomp=(im_test-1)*num_state_material+1+ENUM_DENVAR
          local_density=den(D_DECL(icell,jcell,kcell),dencomp)
          if (local_density.gt.zero) then
           ! do nothing
          else
           print *,"local_density invalid"
           stop
          endif
          mass=local_density*volgrid*mofdata_tess(vofcomp)
          do dir=1,SDIM
           gravvector(dir)=mass*gravity_vector(dir)
           rvec(dir)=mofdata_tess(vofcomp+dir)+ &
            cengrid(dir)-global_centroid(dir)
          enddo
          if (SDIM.eq.2) then
           rvec(3)=zero
           gravvector(3)=zero
          endif 

          do dir=1,SDIM
           ibase=DRAGCOMP_IQ_BODYFORCE+3*(im_test-1)+dir
           localsum(ibase)=localsum(ibase)+gravvector(dir)
          enddo

          ! torque=r cross F
          ! a x b=(a2 b3 - a3 b2) hat(i)+
          !       (-a1 b3 + b1 a3) hat(j)+
          !       (a1 b2 - a2 b1) hat(k)
          call crossprod(rvec,gravvector,rcross)
          do dir=1,SDIM
           grav_localtorque(dir)=zero
          enddo
          grav_localtorque(1)=rcross(3) ! x-y plane of rotation
          if (SDIM.eq.3) then
           grav_localtorque(2)=rcross(1) ! y-z plane of rotation
           grav_localtorque(SDIM)=rcross(2) ! x-z plane of rotation
          endif

          ibase=DRAGCOMP_IQ_MOMINERTIA+3*(im_test-1)
          localsum(ibase+1)=localsum(ibase+1)+ &
            (rvec(1)**2+rvec(2)**2)*mass ! x-y
          dirend=1
          if (SDIM.eq.3) then
           localsum(ibase+2)=localsum(ibase+2)+ &
              (rvec(2)**2+rvec(3)**2)*mass ! y-z
           localsum(ibase+3)=localsum(ibase+3)+ &
              (rvec(1)**2+rvec(3)**2)*mass ! x-z
           dirend=3
          endif
 
          ibase=DRAGCOMP_IQ_BODYTORQUE+3*(im_test-1)
          do dir=1,dirend
           localsum(ibase+dir)=localsum(ibase+dir)+grav_localtorque(dir)
          enddo

          ! cells - pressure, viscosity, viscoelastic
          ! if im_primary<>im_test,
          !  then a force is applied from (icell,jcell,kcell) to a 
          !  neighbor "im_test" cell.

          do im=1,num_materials
           ls_sort(im)=levelpc(D_DECL(icell,jcell,kcell),im)
          enddo
           ! im_primary is the material applying a given force/torque onto
           ! "im_test"
           ! "get_primary_material" declared in GLOBALUTIL.F90
          call get_primary_material(ls_sort,im_primary)

          if (im_primary.eq.im_test) then
           ! do nothing
          else if ((im_primary.ne.im_test).and. &
                   (im_primary.ge.1).and. &
                   (im_primary.le.num_materials)) then

           volume=vol(D_DECL(icell,jcell,kcell))
           if (volume.gt.zero) then
            ! do nothing
           else
            print *,"volume invalid"
            stop
           endif
           presmag=pres(D_DECL(icell,jcell,kcell))

            ! check if im_test is a neighbor.
           do facedir=1,SDIM
            ii=0
            jj=0
            kk=0
            if (facedir.eq.1) then
             ii=1
            else if (facedir.eq.2) then
             jj=1
            else if ((facedir.eq.3).and.(SDIM.eq.3)) then
             kk=1
            else
             print *,"facedir invalid"
             stop
            endif

            do side_cell=0,1

             if (side_cell.eq.0) then
              i_side=icell-ii
              j_side=jcell-jj
              k_side=kcell-kk
             else if (side_cell.eq.1) then
              i_side=icell+ii
              j_side=jcell+jj
              k_side=kcell+kk
             else
              print *,"side_cell invalid"
              stop
             endif
             imac=icell+side_cell*ii
             jmac=jcell+side_cell*jj
             kmac=kcell+side_cell*kk

             do im=1,num_materials
              ls_side(im)=levelpc(D_DECL(i_side,j_side,k_side),im)
             enddo
             call get_primary_material(ls_side,im_side)

             if (im_side.eq.im_test) then

              if (facedir.eq.1) then
               facearea=areax(D_DECL(imac,jmac,kmac))
              else if (facedir.eq.2) then
               facearea=areay(D_DECL(imac,jmac,kmac))
              else if ((facedir.eq.3).and.(SDIM.eq.3)) then
               facearea=areaz(D_DECL(imac,jmac,kmac))
              else
               print *,"facedir invalid"
               stop
              endif
              if (facearea.ge.zero) then
               ! do nothing
              else
               print *,"facearea invalid"
               stop
              endif

              ! nsolid points into the forced material
              do dir=1,3
               nsolid(dir)=zero
              enddo
              if (side_cell.eq.0) then
               nsolid(facedir)=-one
              else if (side_cell.eq.1) then
               nsolid(facedir)=one
              else 
               print *,"side_cell invalid"
               stop
              endif

              drag(D_DECL(icell,jcell,kcell),DRAGCOMP_FLAG+im_test)=1

              !facedir=1..sdim
              call gridstenMAC(xsten_face,xlo,imac,jmac,kmac,fablo,bfact, &
               dx,nhalf,facedir-1)

              ! im_primary is the forcing fluid at cell (icell,jcell,kcell) 
              ! which is applying a force/torque to im_test.
              do dir_visc=1,SDIM
               ii_visc=0
               jj_visc=0
               kk_visc=0
               if (dir_visc.eq.1) then
                ii_visc=1
               else if (dir_visc.eq.2) then
                jj_visc=1
               else if ((dir_visc.eq.3).and.(SDIM.eq.3)) then
                kk_visc=1
               else
                print *,"dir_visc invalid"
                stop
               endif

               do side_visc=1,2

                if (side_visc.eq.1) then
                 i_side_visc=icell-ii_visc
                 j_side_visc=jcell-jj_visc
                 k_side_visc=kcell-kk_visc
                 i_face=icell
                 j_face=jcell
                 k_face=kcell
                else if (side_visc.eq.2) then
                 i_side_visc=icell+ii_visc
                 j_side_visc=jcell+jj_visc
                 k_side_visc=kcell+kk_visc
                 i_face=icell+ii_visc
                 j_face=jcell+jj_visc
                 k_face=kcell+kk_visc
                else
                 print *,"side_visc invalid"
                 stop
                endif
              
                do im=1,num_materials
                 ls_visc(im)= &
                    levelpc(D_DECL(i_side_visc,j_side_visc,k_side_visc),im)
                enddo
                call get_primary_material(ls_visc,im_visc)
                if (is_rigid(im_visc).eq.1) then

                 partid=0
                 do im=1,im_visc-1
                  if (is_lag_part(im).eq.1) then
                   partid=partid+1
                  endif
                 enddo
                 if (im_solid_map(partid+1)+1.ne.im_visc) then
                  print *,"im_solid_map(partid+1)+1.ne.im_visc"
                  stop
                 endif
                 ibase=partid*SDIM
                 do dir=1,SDIM
                  if (dir_visc.eq.1) then
                   vel6point(dir_visc,side_visc,dir)= &
                    solxfab(D_DECL(i_face,j_face,k_face),ibase+dir)
                  else if (dir_visc.eq.2) then
                   vel6point(dir_visc,side_visc,dir)= &
                    solyfab(D_DECL(i_face,j_face,k_face),ibase+dir)
                  else if ((dir_visc.eq.3).and.(SDIM.eq.3)) then
                   vel6point(dir_visc,side_visc,dir)= &
                    solzfab(D_DECL(i_face,j_face,k_face),ibase+dir)
                  else
                   print *,"dir_visc invalid"
                   stop
                  endif
                 enddo ! dir=1..sdim

                else if (is_rigid(im_visc).eq.0) then
                 do dir=1,SDIM
                  vel6point(dir_visc,side_visc,dir)= &
                   half*(vel(D_DECL(icell,jcell,kcell),dir)+ &
                         vel(D_DECL(i_side_visc,j_side_visc,k_side_visc),dir))
                 enddo
                else
                 print *,"is_rigid(im_visc) invalid"
                 stop
                endif

               enddo ! side_visc=1,2
              enddo ! dir_visc=1..sdim
                      
              do veldir=1,3
               do dir=1,3
                gradu(veldir,dir)=zero
               enddo
              enddo

              do veldir=1,SDIM
               do dir=1,SDIM
                delx=xsten(1,dir)-xsten(-1,dir)
                if (delx.gt.zero) then
                 gradu(veldir,dir)= &
                  (vel6point(dir,2,veldir)-vel6point(dir,1,veldir))/delx
                else
                 print *,"delx invalid"
                 stop
                endif
               enddo ! dir=1..sdim
              enddo ! veldir=1..sdim 

              do j1=1,3
              do i1=1,3
               Q(i1,j1)=zero
              enddo
              enddo

              ! if the forcing fluid "im_primary" is a viscoelastic material,
              ! then the viscoelastic force onto "im_test" must be considered.
              if (fort_built_in_elastic_model( &
                   fort_elastic_viscosity(im_primary), &
                   fort_viscoelastic_model(im_primary)).eq.1) then 
               partid=1
               do while ((fort_im_elastic_map(partid)+1.ne.im_primary).and. &
                         (partid.le.num_materials_viscoelastic))
                partid=partid+1
               enddo
               if (partid.le.num_materials_viscoelastic) then
                viscbase=(partid-1)*ENUM_NUM_TENSOR_TYPE
                do dir=1,ENUM_NUM_TENSOR_TYPE
                 call stress_index(dir,i1,j1)
                 Q(i1,j1)=viscoten(D_DECL(icell,jcell,kcell),viscbase+dir)
                enddo
                Q(2,1)=Q(1,2)
                Q(3,1)=Q(1,3)
                Q(3,2)=Q(2,3)
               else
                print *,"partid invalid in fort_getdrag"
                stop
               endif
              else if (fort_built_in_elastic_model( &
                        fort_elastic_viscosity(im_primary), &
                        fort_viscoelastic_model(im_primary)).eq.0) then 
               ! do nothing
              else
               print *,"fort_built_in_elastic_model invalid"
               stop
              endif

              mu_0=get_user_viscconst(im_primary, &
               fort_denconst(im_primary),fort_tempconst(im_primary))
              ! c_mat_visc is initialized in NavierStokes2.cpp: getStateVISC,
              ! getStateVISC_ALL.  "c_mat_visc" includes WALE model and
              ! "viscconst_eddy_wall" effects.
              mu_non_ambient= &
               c_mat_visc(D_DECL(icell,jcell,kcell),im_primary)
              if ((mu_0.ge.zero).and.(mu_non_ambient.ge.zero)) then
               ! do nothing
              else
               print *,"mu_0 or mu_non_ambient invalid"
               stop
              endif

              do j1=1,3
              do i1=1,3

               viscous0_stress_tensor(i1,j1)=mu_0*visc_coef* &
                 (gradu(i1,j1)+gradu(j1,i1))
               viscous_stress_tensor(i1,j1)=mu_non_ambient*visc_coef* &
                 (gradu(i1,j1)+gradu(j1,i1))
               visco_stress_tensor(i1,j1)=Q(i1,j1)

               pressure_stress_tensor(i1,j1)=zero
               if (i1.eq.j1) then
                pressure_stress_tensor(i1,j1)=presmag
               endif

              enddo !i1
              enddo !j1

              do j1=1,3
               viscous_stress_load(j1)=zero
               viscous0_stress_load(j1)=zero
               visco_stress_load(j1)=zero
               pressure_load(j1)=zero

               do i1=1,3
                viscous0_stress_load(j1)= &
                 viscous0_stress_load(j1)- &
                   viscous0_stress_tensor(i1,j1)*nsolid(i1)
                viscous_stress_load(j1)= &
                 viscous_stress_load(j1)- &
                   viscous_stress_tensor(i1,j1)*nsolid(i1)
                visco_stress_load(j1)= &
                 visco_stress_load(j1)- &
                  visco_stress_tensor(i1,j1)*nsolid(i1)
               enddo ! i1=1,3

               pressure_load(j1)=presmag*nsolid(j1)*facearea
               viscous0_stress_load(j1)=viscous0_stress_load(j1)*facearea
               viscous_stress_load(j1)=viscous_stress_load(j1)*facearea
               visco_stress_load(j1)=visco_stress_load(j1)*facearea
              enddo ! j1=1,3

! global_centroid will be incorrect if the solid geometry is reflected
! across a domain boundary.

              do dir=1,SDIM
               rvec(dir)=xsten(0,dir)-global_centroid(dir)
              enddo
              if (SDIM.eq.2) then
               rvec(3)=zero
               viscous_stress_load(3)=zero
               viscous0_stress_load(3)=zero
               visco_stress_load(3)=zero
               pressure_load(3)=zero
              endif

              ibase=DRAGCOMP_TORQUE_ARM+3*(im_test-1)
              do dir=1,SDIM
               drag(D_DECL(icell,jcell,kcell),ibase+dir)=rvec(dir)
              enddo

              call crossprod(rvec,visco_stress_load,visco_rcross)
              call crossprod(rvec,viscous0_stress_load,viscous0_rcross)
              call crossprod(rvec,viscous_stress_load,viscous_rcross)
              call crossprod(rvec,pressure_load,pressure_rcross)

              do dir=1,SDIM
               viscous_localtorque(dir)=zero
               viscous0_localtorque(dir)=zero
               visco_localtorque(dir)=zero
               pressure_localtorque(dir)=zero
              enddo

              viscous_localtorque(1)=viscous_rcross(3) ! x-y plane of rotation
              viscous0_localtorque(1)=viscous0_rcross(3) ! x-y plane of rotation
              visco_localtorque(1)=visco_rcross(3) ! x-y plane of rotation
              pressure_localtorque(1)=pressure_rcross(3) ! x-y plane of rotation

              if (SDIM.eq.3) then

               viscous_localtorque(2)=viscous_rcross(1) !y-z plane of rotation
               viscous0_localtorque(2)=viscous0_rcross(1)!y-z plane of rotation
               visco_localtorque(2)=visco_rcross(1) ! y-z plane of rotation
               pressure_localtorque(2)=pressure_rcross(1) !y-z plane of rotation

               viscous_localtorque(SDIM)=viscous_rcross(2)!x-z plane of rotation
               viscous0_localtorque(SDIM)= &
                      viscous0_rcross(2)!x-z plane of rotation
               visco_localtorque(SDIM)=visco_rcross(2) ! x-z plane rotation
               pressure_localtorque(SDIM)=pressure_rcross(2)!x-z plane rotation
              endif

! units of watts: W=J/s=N m/s
! r=(x,0,z)
! f=(f1,f2,f3)
! torque=r x f = (y f3 - z f2)ihat+(z f1 - x f3)jhat+(x f2 - y f1)khat
! power= torque dot w
! rotation about y axis
! vinletgas RPM
! w=2 pi vinletgas/60 radians/second
! so for gear problem, scale torque by 2 pi abs(vinletgas)/60

              ! im_test is the material on which a force/torque is applied.

              do dir=1,6

               call stress_index(dir,i1,j1)
               ibase=DRAGCOMP_STRESS+6*(im_test-1)+dir
               drag(D_DECL(icell,jcell,kcell),ibase)= &
                     pressure_stress_tensor(i1,j1)+ &
                     viscous_stress_tensor(i1,j1)+ &
                     visco_stress_tensor(i1,j1)

               ibase=DRAGCOMP_PSTRESS+6*(im_test-1)+dir
               drag(D_DECL(icell,jcell,kcell),ibase)= &
                 pressure_stress_tensor(i1,j1)

               ibase=DRAGCOMP_VISCOUSSTRESS+6*(im_test-1)+dir
               drag(D_DECL(icell,jcell,kcell),ibase)= &
                 viscous_stress_tensor(i1,j1)

               ibase=DRAGCOMP_VISCOUS0STRESS+6*(im_test-1)+dir
               drag(D_DECL(icell,jcell,kcell),ibase)= &
                 viscous0_stress_tensor(i1,j1)

               ibase=DRAGCOMP_VISCOSTRESS+6*(im_test-1)+dir
               drag(D_DECL(icell,jcell,kcell),ibase)= &
                 visco_stress_tensor(i1,j1)

              enddo ! dir=1..6

              do dir=1,SDIM
               ibase=DRAGCOMP_IQ_FORCE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                pressure_load(dir)+ &
                viscous_stress_load(dir)+ &
                visco_stress_load(dir)

               ibase=DRAGCOMP_IQ_PFORCE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                pressure_load(dir)

               ibase=DRAGCOMP_IQ_VISCOUSFORCE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                viscous_stress_load(dir)

               ibase=DRAGCOMP_IQ_VISCOUS0FORCE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                viscous0_stress_load(dir)

               ibase=DRAGCOMP_IQ_VISCOFORCE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                visco_stress_load(dir)

              enddo ! dir=1..sdim

              dirend=1
              if (SDIM.eq.3) then
               dirend=3
              endif

              do dir=1,dirend

               ibase=DRAGCOMP_IQ_TORQUE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                pressure_localtorque(dir)+ &
                viscous_localtorque(dir)+ &
                visco_localtorque(dir)

               ibase=DRAGCOMP_IQ_PTORQUE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                pressure_localtorque(dir)

               ibase=DRAGCOMP_IQ_VISCOUSTORQUE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                viscous_localtorque(dir)

               ibase=DRAGCOMP_IQ_VISCOUS0TORQUE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                viscous0_localtorque(dir)

               ibase=DRAGCOMP_IQ_VISCOTORQUE+3*(im_test-1)+dir

               localsum(ibase)=localsum(ibase)+ &
                visco_localtorque(dir)

              enddo ! dir=1..dirend

              ibase=DRAGCOMP_IQ_PERIM+im_test
              localsum(ibase)=localsum(ibase)+facearea

             else if ((im_side.ne.im_test).and. &
                      (im_side.ge.1).and. &
                      (im_side.le.num_materials)) then
              ! do nothing
             else
              print *,"im_side invalid"
              stop
             endif

            enddo ! side_cell=0,1
           enddo ! facedir=1..sdim

          else
           print *,"im_primary or im_test corruption"
           stop
          endif

         enddo ! im_test=1..num_materials

        else if (mask_cell.eq.0) then
         ! do nothing
        else
         print *,"mask_cell invalid"
         stop
        endif

       else
        print *,"isweep invalid"
        stop
       endif

      enddo
      enddo
      enddo  ! icell,jcell,kcell

      return
      end subroutine fort_getdrag



      subroutine fort_maxpresvel( &
        minpres, &
        maxpres, &
        maxvel, &
        maxvel_collide, &
        xlo,dx, &
        mask,DIMS(mask), &
        vel,DIMS(vel), &
        velx,DIMS(velx), &
        vely,DIMS(vely), &
        velz,DIMS(velz), &
        tilelo,tilehi, &
        fablo,fabhi,bfact) &
      bind(c,name='fort_maxpresvel')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      real(amrex_real), INTENT(inout) :: minpres,maxpres,maxvel,maxvel_collide

      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer, INTENT(in) :: bfact
      integer growlo(3), growhi(3)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(velx)
      integer, INTENT(in) :: DIMDEC(vely)
      integer, INTENT(in) :: DIMDEC(velz)

      real(amrex_real), INTENT(in),target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in),target :: &
          vel(DIMV(vel),STATE_NCOMP_VEL+STATE_NCOMP_PRES)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: velx(DIMV(velx))
      real(amrex_real), pointer :: velx_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: vely(DIMV(vely))
      real(amrex_real), pointer :: vely_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: velz(DIMV(velz))
      real(amrex_real), pointer :: velz_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      integer dir
      real(amrex_real) magvel,magvel_mac,magvel_collide
      real(amrex_real) vello,velhi,vellohi
      integer i,j,k

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      mask_ptr=>mask
      vel_ptr=>vel
      velx_ptr=>velx
      vely_ptr=>vely
      velz_ptr=>velz
      call checkbound_array1(fablo,fabhi,mask_ptr,0,-1)
      call checkbound_array(fablo,fabhi,vel_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,velx_ptr,0,0)
      call checkbound_array1(fablo,fabhi,vely_ptr,0,1)
      call checkbound_array1(fablo,fabhi,velz_ptr,0,SDIM-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
    
       if (mask(D_DECL(i,j,k)).eq.one) then

        if (vel(D_DECL(i,j,k),STATECOMP_PRES+1).gt.maxpres) then
          maxpres=vel(D_DECL(i,j,k),STATECOMP_PRES+1)
        endif
        if (vel(D_DECL(i,j,k),STATECOMP_PRES+1).lt.minpres) then
          minpres=vel(D_DECL(i,j,k),STATECOMP_PRES+1)
        endif
        magvel=zero
        magvel_mac=zero
        magvel_collide=zero
        do dir=1,SDIM
         magvel=magvel+vel(D_DECL(i,j,k),dir)**2
         if (dir.eq.1) then
          vello=velx(D_DECL(i,j,k))
          velhi=velx(D_DECL(i+1,j,k))
         else if (dir.eq.2) then
          vello=vely(D_DECL(i,j,k))
          velhi=vely(D_DECL(i,j+1,k))
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          vello=velz(D_DECL(i,j,k))
          velhi=velz(D_DECL(i,j,k+1))
         else
          print *,"dir invalid"
          stop
         endif
         magvel_collide=magvel_collide+(velhi-vello)**2
         vellohi=max(abs(vello),abs(velhi))
         magvel_mac=magvel_mac+vellohi**2
        enddo ! dir=1..sdim
        magvel=sqrt(magvel)
        magvel_mac=sqrt(magvel_mac)
        magvel_collide=sqrt(magvel_collide)
        if (magvel.gt.maxvel) then
         maxvel=magvel
        endif
        if (magvel_mac.gt.maxvel) then
         maxvel=magvel_mac
        endif
        if (magvel_collide.gt.maxvel_collide) then
         maxvel_collide=magvel_collide
        endif
       else if (mask(D_DECL(i,j,k)).ne.zero) then
        print *,"mask invalid"
        stop 
       endif
! mask=1
      enddo
      enddo
      enddo

      return
      end subroutine fort_maxpresvel

      end module derive_module

