#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "DRAG_COMP.H"
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
       im, &
       nmat, &
       dt, &
       ntensor, &  ! for declaring cellten
       denstate,DIMS(denstate), &
       vof,DIMS(vof), &
       vel,DIMS(vel), &
       visc,DIMS(visc), &
       cellten,DIMS(cellten), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       cur_time, &
       dx,xlo, &
       ngrow, &
       ncompvisc) &
      bind(c,name='fort_derturbvisc')

      use global_utility_module
      use probf90_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: les_model
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: im
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: ncompvisc
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: ntensor
      REAL_T, intent(in) :: dt
      REAL_T, intent(in) :: cur_time

      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: dx(SDIM),xlo(SDIM)

      INTEGER_T, intent(in) :: DIMDEC(denstate)
      INTEGER_T, intent(in) :: DIMDEC(vof)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(visc)
      INTEGER_T, intent(in) :: DIMDEC(cellten)
  
      REAL_T, intent(in), target :: &
              denstate(DIMV(denstate),nmat*num_state_material)
      REAL_T, pointer :: denstate_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T, pointer :: vof_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: vel(DIMV(vel),SDIM)
      REAL_T, pointer :: vel_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(inout), target :: visc(DIMV(visc),ncompvisc)
      REAL_T, pointer :: visc_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: cellten(DIMV(cellten),ntensor)
      REAL_T, pointer :: cellten_ptr(D_DECL(:,:,:),:)

      REAL_T g(3,3),s(3,3),sd(3,3),g2(3),g2tr,ss,sdsd
      REAL_T g2_12,g2_13,g2_21,g2_23,g2_31,g2_32,turb_visc,density
      REAL_T rr
      INTEGER_T i,j,k
      INTEGER_T i1,j1,k1
      INTEGER_T k1lo,k1hi
      INTEGER_T veldir,dir
      INTEGER_T flagcomp
      INTEGER_T vofcomp
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz,nbase
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf
      INTEGER_T near_solid
      INTEGER_T caller_id
      INTEGER_T im_primary
      INTEGER_T im_local
      REAL_T VFRAC(nmat)
     
      nhalf=1

      caller_id=113

      if (ngrow.eq.1) then
       ! do nothing
      else
       print *,"expecting ngrow==1"
       stop
      endif

      if ((ncompvisc.ne.nmat).and.(ncompvisc.ne.3*nmat)) then
       print *,"ncompvisc invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid1"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid1"
       stop
      endif
      if (dt.lt.zero) then 
       print *,"dt invalid"
       stop
      endif
      if (cur_time.lt.zero) then 
       print *,"cur_time invalid"
       stop
      endif
      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
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

      ! compute u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)
    
      vof_ptr=>vof 
      call checkbound_array(fablo,fabhi,vof_ptr,ngrow+1,-1,310)
      visc_ptr=>visc
      call checkbound_array(fablo,fabhi,visc_ptr,ngrow,-1,311)
      vel_ptr=>vel
      call checkbound_array(fablo,fabhi,vel_ptr,ngrow,-1,312)
      denstate_ptr=>denstate
      call checkbound_array(fablo,fabhi,denstate_ptr,ngrow,-1,313)
      cellten_ptr=>cellten
      call checkbound_array(fablo,fabhi,cellten_ptr,ngrow,-1,314)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      if (is_rigid(nmat,im).eq.1) then
       print *,"cannot have eddy viscosity and (is_rigid(nmat,im).eq.1)"
       stop
      endif

       ! density
      flagcomp=(im-1)*num_state_material+1

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       do veldir=1,3
       do dir=1,3
        g(veldir,dir)=zero
       enddo
       enddo

       do dir=1,SDIM
        if (dir.eq.1) then
         nbase=ux-1
        else if (dir.eq.2) then
         nbase=uy-1
        else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
         nbase=uz-1
        else
         print *,"dir invalid get shear"
         stop
        endif
        do veldir=1,SDIM
         g(veldir,dir)=cellten(D_DECL(i,j,k),nbase+veldir) 
        enddo
       enddo ! dir

       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        rr=xsten(0,1)
        g(3,3)=vel(D_DECL(i,j,k),1)/abs(rr)
       else if (levelrz.eq.3) then
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

       ss=s(1,1)*s(1,1)+s(1,2)*s(1,2)+s(1,3)*s(1,3)+&
          s(2,1)*s(2,1)+s(2,2)*s(2,2)+s(2,3)*s(2,3)+&
          s(3,1)*s(3,1)+s(3,2)*s(3,2)+s(3,3)*s(3,3)

       ! Builds the traceless symmetric part of the square of the
       ! velocity tensor
       do veldir=1,3
        g2(veldir)=zero
        do dir=1,3
         g2(veldir)=g2(veldir)+g(veldir,dir)*g(dir,veldir)
        enddo
       enddo
       g2tr=(g2(1)+g2(2)+g2(3))/3.0d0 ! tensor trace
       g2_12=g(1,1)*g(1,2)+g(1,2)*g(2,2)+g(1,3)*g(3,2)
       g2_13=g(1,1)*g(1,3)+g(1,2)*g(2,3)+g(1,3)*g(3,3)
       g2_21=g(2,1)*g(1,1)+g(2,2)*g(2,1)+g(2,3)*g(3,1)
       g2_23=g(2,1)*g(1,3)+g(2,2)*g(2,3)+g(2,3)*g(3,3)
       g2_31=g(3,1)*g(1,1)+g(3,2)*g(2,1)+g(3,3)*g(3,1)
       g2_32=g(3,1)*g(1,2)+g(3,2)*g(2,2)+g(3,3)*g(3,2)
 
       sd(1,1)=g2(1)-g2tr
       sd(1,2)=0.5d0*(g2_12+g2_21)
       sd(1,3)=0.5d0*(g2_13+g2_31)

       sd(2,1)=sd(1,2)
       sd(2,2)=g2(2)-g2tr
       sd(2,3)=0.5*(g2_23+g2_23)

       sd(3,1)=sd(1,3)
       sd(3,2)=sd(2,3)
       sd(3,3)=g2(3)-g2tr

       ! tensor squared
       sdsd=sd(1,1)*sd(1,1)+sd(1,2)*sd(1,2)+sd(1,3)*sd(1,3)+&
            sd(2,1)*sd(2,1)+sd(2,2)*sd(2,2)+sd(2,3)*sd(2,3)+&
            sd(3,1)*sd(3,1)+sd(3,2)*sd(3,2)+sd(3,3)*sd(3,3)

       if (les_model.eq.1) then

        if(sdsd.ne.zero) then
         turb_visc=(0.5d0*dx(1))**2*sdsd**1.5d0/(ss**2.5d0+sdsd**1.25d0)
         density=denstate(D_DECL(i,j,k),flagcomp)
         ! limiter: arbitrary at this point
         turb_visc=min(turb_visc,two*visc(D_DECL(i,j,k),im)/density)
         visc(D_DECL(i,j,k),im)=visc(D_DECL(i,j,k),im)+turb_visc*density
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

       if (fort_viscconst_eddy(im).gt.zero) then
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
         do im_local=1,nmat
          vofcomp=(im_local-1)*ngeom_recon+1
          VFRAC(im_local)=vof(D_DECL(i+i1,j+j1,k+k1),vofcomp)
         enddo
         call get_primary_material_VFRAC(VFRAC,nmat,im_primary,caller_id)
         if (is_rigid(nmat,im_primary).eq.1) then
          near_solid=1
         else if (is_rigid(nmat,im_primary).eq.0) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif
        enddo !k1
        enddo !j1
        enddo !i1

        if (near_solid.eq.1) then
         visc(D_DECL(i,j,k),im)=visc(D_DECL(i,j,k),im)+ &
                 fort_viscconst_eddy(im)
        else if (near_solid.eq.0) then
         ! do nothing
        else
         print *,"near_solid invalid"
         stop
        endif

       else if (fort_viscconst_eddy(im).eq.zero) then
        ! do nothing
       else
        print *,"fort_viscconst_eddy invalid"
        stop
       endif
  
      enddo 
      enddo 
      enddo 

      return
      end subroutine fort_derturbvisc

      subroutine fort_getshear( &
       ntensor, &
       cellten,DIMS(cellten), &
       vel,DIMS(vel), &
       dx,xlo, &
       tensordata, &
       DIMS(tensordata), &
       iproject,onlyscalar, &
       time, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       bc, &
       ngrow,nmat) &
      bind(c,name='fort_getshear')

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level,ntensor
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf
      INTEGER_T i,j,k,n,dir,veldir
      INTEGER_T, intent(in) :: iproject,ngrow,onlyscalar
      INTEGER_T i1,j1
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: bc(SDIM,2,SDIM)
      INTEGER_T, intent(in) :: DIMDEC(cellten)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(tensordata)
  
      REAL_T, intent(in), target :: cellten(DIMV(cellten),ntensor)
      REAL_T, pointer :: cellten_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(in), target :: vel(DIMV(vel),SDIM)
      REAL_T, pointer :: vel_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(out), target :: tensordata(DIMV(tensordata),20)
      REAL_T, pointer :: tensordata_ptr(D_DECL(:,:,:),:)

      REAL_T visctensor(3,3),gradu(3,3)
      REAL_T shear
      REAL_T a,b,c
      REAL_T rr
      REAL_T vort(3)
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz,nbase

      nhalf=1

      tensordata_ptr=>tensordata
      cellten_ptr=>cellten
      vel_ptr=>vel

      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid2"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      ! compute u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)

      call checkbound_array(fablo,fabhi,cellten_ptr,ngrow,-1,64)
      call checkbound_array(fablo,fabhi,vel_ptr,ngrow+1,-1,64)
      call checkbound_array(fablo,fabhi, &
       tensordata_ptr, &
       ngrow,-1,65)

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

        !cellten:
        !u_x,v_x,w_x,u_y,v_y,w_y,u_z,v_z,w_z
       do dir=1,SDIM
        if (dir.eq.1) then
         nbase=ux-1
        else if (dir.eq.2) then
         nbase=uy-1
        else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
         nbase=uz-1
        else
         print *,"dir invalid get shear"
         stop
        endif
        do veldir=1,SDIM
         gradu(veldir,dir)=cellten(D_DECL(i,j,k),nbase+veldir) 
        enddo
       enddo ! dir

       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        rr=xsten(0,1)
        gradu(3,3)=vel(D_DECL(i,j,k),1)/abs(rr)
       else if (levelrz.eq.3) then
        rr=xsten(0,1)
        gradu(2,2)=gradu(2,2)+vel(D_DECL(i,j,k),1)/abs(rr)
        gradu(1,2)=gradu(1,2)-vel(D_DECL(i,j,k),2)/abs(rr)
       else
        print *,"levelrz invalid getshear 2"
        stop
       endif

! project trace of gradu to zero (valid for incompressible flow)
! ( iproject==0 if called from NavierStokes::tensor_advection_update() )
       if (iproject.eq.1) then

#if (AMREX_SPACEDIM==3)
        a=(two*gradu(1,1)-gradu(2,2)-gradu(3,3))/three
        b=(two*gradu(2,2)-gradu(1,1)-gradu(3,3))/three
        c=(two*gradu(3,3)-gradu(1,1)-gradu(2,2))/three
#elif (AMREX_SPACEDIM==2)
        a=(gradu(1,1)-gradu(2,2)-gradu(3,3))/two
        b=(gradu(2,2)-gradu(1,1)-gradu(3,3))/two
        c=gradu(3,3)
#else
        print *,"dimension bust"
        stop
#endif
        gradu(1,1)=a
        gradu(2,2)=b
        gradu(3,3)=c
       else if (iproject.eq.0) then
        ! do nothing
       else 
        print *,"iproject invalid"
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

       if (onlyscalar.eq.1) then
        tensordata(D_DECL(i,j,k),1)=shear
       else if (onlyscalar.eq.2) then
        vort(1)=gradu(3,2)-gradu(2,3)
        vort(2)=gradu(1,3)-gradu(3,1)
        vort(3)=gradu(2,1)-gradu(1,2)
        tensordata(D_DECL(i,j,k),1)= &
          sqrt(vort(1)**2+vort(2)**2+vort(3)**2)
       else if (onlyscalar.eq.3) then
        tensordata(D_DECL(i,j,k),1)=abs(gradu(1,1)+gradu(2,2)+gradu(3,3))
       else if (onlyscalar.eq.0) then
        tensordata(D_DECL(i,j,k),1)=shear
        n=2
        do i1=1,3
        do j1=1,3
         tensordata(D_DECL(i,j,k),n)=visctensor(i1,j1)
         n=n+1
        enddo
        enddo
        do i1=1,3
        do j1=1,3
         tensordata(D_DECL(i,j,k),n)=gradu(i1,j1) !gradu(veldir,dir)
         n=n+1
        enddo
        enddo
       else 
        print *,"onlyscalar invalid"
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
        im_parm, &
        nmat, &
        dt, &
        viscosity_coefficient, & ! viscconst
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
        elastic_regularization, &
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

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: visc_coef
      INTEGER_T, intent(in) :: im_parm
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: ncompvisc
      REAL_T, intent(in) :: dt
      REAL_T, intent(in) :: viscosity_coefficient
      INTEGER_T, intent(in) :: shear_thinning_fluid
      REAL_T, intent(in) :: Carreau_alpha
      REAL_T, intent(in) :: Carreau_beta
      REAL_T, intent(in) :: Carreau_n
      REAL_T, intent(in) :: Carreau_mu_inf
      REAL_T, intent(in) :: concentration,etaL,etaP,etaS
      REAL_T, intent(in) :: elastic_time
      REAL_T, intent(in) :: elastic_viscosity
      REAL_T, intent(in) :: elastic_regularization
      INTEGER_T, intent(in) :: viscosity_state_model
      INTEGER_T, intent(in) :: viscoelastic_model
      REAL_T, intent(in) :: polymer_factor
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(visc)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(eosdata)
      INTEGER_T, intent(in) :: DIMDEC(tensor)
      INTEGER_T, intent(in) :: DIMDEC(gammadot)
      INTEGER_T, intent(in) :: bc(SDIM,2,SDIM)
      INTEGER_T, intent(in) :: ngrow
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: dx(SDIM), xlo(SDIM)
       !ncompvisc=3*nmat
      REAL_T, intent(out), target :: visc(DIMV(visc),ncompvisc) 
      REAL_T, pointer :: visc_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(in), target :: vel(DIMV(vel),SDIM)
      REAL_T, pointer :: vel_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: gammadot(DIMV(gammadot))
      REAL_T, pointer :: gammadot_ptr(D_DECL(:,:,:))
      REAL_T, intent(in), target :: eosdata(DIMV(eosdata), &
              nmat*num_state_material)
      REAL_T, pointer :: eosdata_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: tensor(DIMV(tensor),FORT_NUM_TENSOR_TYPE)
      REAL_T, pointer :: tensor_ptr(D_DECL(:,:,:),:)

      INTEGER_T i,j,k
      REAL_T    shear,density,temperature,mu
      INTEGER_T flagcomp
      REAL_T bterm,pterm
      INTEGER_T numstatetest,ii,jj
      REAL_T Q(3,3)
      REAL_T traceA,modtime,viscoelastic_coeff
      REAL_T bulk_modulus

      visc_ptr=>visc

      if (bfact.lt.1) then
       print *,"bfact invalid3"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 31"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((im_parm.lt.1).or.(im_parm.gt.nmat)) then
       print *,"im_parm invalid3"
       stop
      endif

      if (elastic_regularization.ge.zero) then
       ! do nothing
      else
       print *,"elastic_regularization must be >= 0.0"
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

      if (fort_is_eulerian_elastic_model(elastic_viscosity, &
              viscoelastic_model).eq.1) then
       ! do nothing
      else if (fort_is_eulerian_elastic_model(elastic_viscosity, &
                 viscoelastic_model).eq.0) then
       ! do nothing
      else
       print *,"fort_is_eulerian_elastic_model invalid"
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


      if (ncompvisc.ne.3*nmat) then
       print *,"ncompvisc invalid"
       stop
      endif

      call checkbound_array(fablo,fabhi,visc_ptr,ngrow,-1,316)
      gammadot_ptr=>gammadot
      call checkbound_array1(fablo,fabhi,gammadot_ptr,ngrow,-1,317)
      eosdata_ptr=>eosdata
      call checkbound_array(fablo,fabhi,eosdata_ptr,ngrow,-1,318)
      tensor_ptr=>tensor
      call checkbound_array(fablo,fabhi,tensor_ptr,ngrow,-1,319)
      vel_ptr=>vel
      call checkbound_array(fablo,fabhi,vel_ptr,ngrow,-1,320)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      if (is_rigid(nmat,im_parm).eq.1) then
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        visc(D_DECL(i,j,k),im_parm)=viscosity_coefficient 
       enddo
       enddo
       enddo
      else if (is_rigid(nmat,im_parm).eq.0) then

       vel_ptr=>vel
       call checkbound_array(fablo,fabhi,vel_ptr,ngrow+1,-1,321)

       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

         ! den,T
        if (viscosity_state_model.ge.1) then
         flagcomp=(im_parm-1)*num_state_material+1
         density=eosdata(D_DECL(i,j,k),flagcomp)
         temperature=eosdata(D_DECL(i,j,k),flagcomp+1)
        else if (viscosity_state_model.eq.0) then
         density=fort_denconst(im_parm)
         temperature=fort_tempconst(im_parm)
        else
         print *,"viscosity_state_model invalid"
         stop
        endif

        mu=get_user_viscconst(im_parm,density,temperature)

        if ((viscoelastic_model.eq.0).or. & !FENE-CR
            (viscoelastic_model.eq.1).or. & ! Oldroyd B
            (viscoelastic_model.eq.5).or. & ! FENE-P
            (Viscoelastic_model.eq.6)) then ! linear PTT
         ! do nothing
        else if (viscoelastic_model.eq.4) then ! pressure velocity coupling
         ! do nothing
        else if ((viscoelastic_model.eq.2).or. & !displacement gradient
                 (viscoelastic_model.eq.3)) then !incremental
         bulk_modulus=elastic_viscosity
         if (bulk_modulus.gt.zero) then
          if (visc_coef.gt.zero) then
            ! notes: 
            !  viscoelastic_coeff*visc_coef down below.
            !  dd_group=dd*visc_coef in PROB.F90 
            !  xflux*=-dt * visc_coef * facevisc_index in CROSSTERM
           if (mu.ge.zero) then
            if (elastic_regularization.ge.zero) then
             mu=mu+dt*elastic_regularization*bulk_modulus  
            else
             print *,"elastic_regularization invalid"
             stop
            endif
           else
            print *,"mu invalid"
            stop
           endif
          else if (visc_coef.eq.zero) then
           ! do nothing
          else
           print *,"visc_coef invalid"
           stop
          endif
         else if (bulk_modulus.eq.zero) then
          ! do nothing
         else
          print *,"bulk_modulus invalid"
          stop
         endif
        else
         print *,"viscoelastic_model invalid"
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

          if (ncompvisc.ne.3*nmat) then
           print *,"ncompvisc invalid"
           stop
          endif

          numstatetest=num_state_base+num_species_var
          if (num_state_material.ne.numstatetest) then
           print *,"num_state_material invalid"
           stop
          endif
          if ((num_materials_viscoelastic.lt.1).or. &
              (num_materials_viscoelastic.gt.nmat)) then
           print *,"num_materials_viscoelastic invalid:fort_derviscosity"
           stop
          endif

          do ii=1,3
          do jj=1,3
           Q(ii,jj)=zero
          enddo
          enddo
          Q(1,1)=tensor(D_DECL(i,j,k),1)
          Q(1,2)=tensor(D_DECL(i,j,k),2)
          Q(2,2)=tensor(D_DECL(i,j,k),3)
          Q(3,3)=tensor(D_DECL(i,j,k),4)
#if (AMREX_SPACEDIM==3)
          Q(1,3)=tensor(D_DECL(i,j,k),5)
          Q(2,3)=tensor(D_DECL(i,j,k),6)
#endif
          Q(2,1)=Q(1,2)
          Q(3,1)=Q(1,3)
          Q(3,2)=Q(2,3)

          if (SDIM.eq.3) then
           ! do nothing
          else if (SDIM.eq.2) then
           if ((levelrz.eq.0).or.(levelrz.eq.3)) then
            ! do nothing
           else if (levelrz.eq.1) then

            if (i.lt.0) then
             if (Q(3,3).le.-half) then
              Q(3,3)=-half
             else if (Q(3,3).ge.-half) then
              ! do nothing
             else
              print *,"Q(3,3) bust"
              stop
             endif
            else if (i.eq.0) then
             if (Q(3,3).ge.half) then
              Q(3,3)=half
             else if (Q(3,3).le.half) then
              ! do nothing
             else
              print *,"Q(3,3) bust"
              stop
             endif
            else if (i.gt.0) then
             ! do nothing
            else
             print *,"i invalid"
             stop
            endif

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

           ! modtime=elastic_time >> 1 for 2,3
          else if ((viscoelastic_model.eq.2).or. & !displacement gradient
                   (viscoelastic_model.eq.4).or. & !pressure velocity coupling
                   (viscoelastic_model.eq.3)) then !incremental
           modtime=elastic_time
          else
           print *,"viscoelastic_model invalid"
           stop
          endif

          if (modtime.ge.zero) then
           ! do nothing
          else
           print *,"modtime invalid"
           stop
          endif
         
           ! etaS=etaL-etaP=viscconst-elastic_viscosity 
          if (modtime+dt.le.zero) then
           viscoelastic_coeff=zero
          else
           if ((viscoelastic_model.eq.0).or. &!FENE-CR(lambda_tilde=f(A)/lambda)
               (viscoelastic_model.eq.1).or. &!Oldroyd-B(modtime=elastic_time)
               (viscoelastic_model.eq.5).or. &!FENE-P (lambda_tilde=f(A)/lambda
               (viscoelastic_model.eq.6)) then!linearPTT(modtime=elastic_time)
             ! etaS=etaL-etaP=viscconst-elastic_viscosity 
            viscoelastic_coeff= &
             (visc(D_DECL(i,j,k),im_parm)-etaS)/(modtime+dt)
           else if (viscoelastic_model.eq.2) then !displacement gradient
            viscoelastic_coeff=elastic_viscosity
           else if (viscoelastic_model.eq.3) then !incremental
            viscoelastic_coeff=elastic_viscosity
           else if (viscoelastic_model.eq.4) then !pressure velocity coupling
            viscoelastic_coeff=elastic_viscosity
           else
            print *,"viscoelastic_model invalid"
            stop
           endif
          endif
          if (viscoelastic_coeff.ge.zero) then
           ! do nothing
          else
           print *,"viscoelastic_coef invalid"
           stop
          endif

          visc(D_DECL(i,j,k),nmat+im_parm)=viscoelastic_coeff*visc_coef
          visc(D_DECL(i,j,k),2*nmat+im_parm)=modtime

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
        level, &
        finest_level, &
        im_parm, & ! =1..nmat
        nmat, &
        dt, &
        conduct,DIMS(conduct), &
        eosdata,DIMS(eosdata), &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        time, &
        dx,xlo, &
        ngrow) &
      bind(c,name='fort_derconductivity')

      use global_utility_module
      use probf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: im_parm !=1..nmat
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(conduct)
      INTEGER_T, intent(in) :: DIMDEC(eosdata)
      INTEGER_T, intent(in) :: ngrow
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: dx(SDIM), xlo(SDIM)
      REAL_T, intent(out), target :: conduct(DIMV(conduct),nmat) 
      REAL_T, pointer :: conduct_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(in), target :: eosdata(DIMV(eosdata), &
              nmat*num_state_material)
      REAL_T, pointer :: eosdata_ptr(D_DECL(:,:,:),:)

      INTEGER_T i,j,k
      INTEGER_T dir
      INTEGER_T flagcomp
      REAL_T density,temperature
      REAL_T thermal_k
      REAL_T xvec(SDIM)
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf

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

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((im_parm.lt.1).or.(im_parm.gt.nmat)) then
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

      call checkbound_array(fablo,fabhi,conduct_ptr,ngrow,-1,316)
      eosdata_ptr=>eosdata
      call checkbound_array(fablo,fabhi,eosdata_ptr,ngrow,-1,318)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)
       do dir=1,SDIM
        xvec(dir)=xsten(0,dir)
       enddo

       flagcomp=(im_parm-1)*num_state_material+1
       density=eosdata(D_DECL(i,j,k),flagcomp)
       temperature=eosdata(D_DECL(i,j,k),flagcomp+1)

       thermal_k=get_user_heatviscconst(im_parm)

       if (is_in_probtype_list().eq.1) then
        call SUB_THERMAL_K(xvec,dx,time,density,temperature, &
                thermal_k,im_parm)
       endif

       conduct(D_DECL(i,j,k),im_parm) = thermal_k

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
       im, &   ! im=0..nmat-1
       ntensor, &
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
       ncomp_tensor, &
       ncomp_visc, &
       n_trace, &
       nmat, &
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

      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: ntensor
      INTEGER_T, intent(in) :: ncomp_den
      INTEGER_T, intent(in) :: ncomp_tensor
      INTEGER_T, intent(in) :: ncomp_visc
      INTEGER_T, intent(in) :: n_trace
      INTEGER_T, intent(in) :: nmat

      REAL_T, intent(in) :: polymer_factor(nmat)
      REAL_T, intent(in) :: etaS(nmat)
      REAL_T, intent(in) :: etaP(nmat)
      REAL_T, intent(in) :: Carreau_beta(nmat)
      REAL_T, intent(in) :: elastic_time(nmat)
      INTEGER_T, intent(in) :: viscoelastic_model(nmat)
      REAL_T, intent(in) :: elastic_viscosity(nmat)

      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(dest)
      INTEGER_T, intent(in) :: DIMDEC(den)
      INTEGER_T, intent(in) :: DIMDEC(tensor)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(visc)
      INTEGER_T, intent(in) :: DIMDEC(cellten)
      INTEGER_T, intent(in) :: bc(SDIM,2,SDIM)
      INTEGER_T, intent(in) :: ngrow
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: dx(SDIM), xlo(SDIM)
      REAL_T, intent(in), target :: cellten(DIMV(cellten),ntensor)

      REAL_T, intent(inout), target :: dest(DIMV(dest),5)
      REAL_T, pointer :: dest_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(in), target :: den(DIMV(den),ncomp_den)
      REAL_T, intent(in), target :: tensor(DIMV(tensor),ncomp_tensor)
      REAL_T, intent(in), target :: vel(DIMV(vel),SDIM)
      REAL_T, intent(in), target :: visc(DIMV(visc),ncomp_visc)

      INTEGER_T im  ! im=0..nmat-1
      INTEGER_T i,j,k
      INTEGER_T iproject,onlyscalar
      REAL_T T11,T22,T33,traceA,modtime
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz
      INTEGER_T nbase
      INTEGER_T dir,veldir
      REAL_T gradu(3,3)
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T rr
      REAL_T vort(3)

      nhalf=1

      dest_ptr=>dest

      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif
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
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (n_trace.ne.nmat*5) then
       print *,"n_trace invalid"
       stop
      endif
      if (ncomp_den.ne.nmat*num_state_material) then
       print *,"ncomp_den invalid"
       stop
      endif
      if (ncomp_visc.lt.nmat) then
       print *,"ncomp_visc invalid"
       stop
      endif  
      if (fort_is_eulerian_elastic_model(elastic_viscosity(im+1), &
             viscoelastic_model(im+1)).eq.1) then 
       ! do nothing
      else if (fort_is_eulerian_elastic_model(elastic_viscosity(im+1), &
             viscoelastic_model(im+1)).eq.0) then 
       ! do nothing
      else
       print *,"fort_is_eulerian_elastic_model invalid"
       stop
      endif
      if (fort_elastic_viscosity(im+1).ne.elastic_viscosity(im+1)) then
       print *,"elastic_viscosity(im+1) invalid"
       stop
      endif
      if (fort_is_eulerian_elastic_model(elastic_viscosity(im+1), &
             fort_viscoelastic_model(im+1)).eq.1) then 
       if (num_materials_viscoelastic.le.0) then
        print *,"num_materials_viscoelastic.le.0:fort_dermagtrace"
        stop
       endif
       if (ncomp_tensor.ne.FORT_NUM_TENSOR_TYPE) then
        print *,"ncomp_tensor.ne.FORT_NUM_TENSOR_TYPE"
        stop
       endif
      else if (fort_is_eulerian_elastic_model(elastic_viscosity(im+1), &
                 fort_viscoelastic_model(im+1)).eq.0) then 
       if (ncomp_tensor.ne.ncomp_den) then
        print *,"ncomp_tensor.ne.ncomp_den"
        stop
       endif
      else
       print *,"fort_is_eulerian_elastic_model invalid"
       stop
      endif  

      ! compute u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)
 
      call checkbound_array(fablo,fabhi,cellten,ngrow,-1,64)
      call checkbound_array(fablo,fabhi,dest_ptr,ngrow,-1,323)
      call checkbound_array(fablo,fabhi,den,ngrow,-1,324)
      call checkbound_array(fablo,fabhi,tensor,ngrow,-1,325)
      call checkbound_array(fablo,fabhi,vel,ngrow+1,-1,326)
      call checkbound_array(fablo,fabhi,visc,ngrow,-1,327)

      iproject=0
      onlyscalar=1  ! mag(trace gradu)
        
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
         nbase=ux-1
        else if (dir.eq.2) then
         nbase=uy-1
        else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
         nbase=uz-1
        else
         print *,"dir invalid get shear"
         stop
        endif
        do veldir=1,SDIM
         gradu(veldir,dir)=cellten(D_DECL(i,j,k),nbase+veldir) 
        enddo
       enddo ! dir

       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        rr=xsten(0,1)
        gradu(3,3)=vel(D_DECL(i,j,k),1)/abs(rr)
       else if (levelrz.eq.3) then
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
       if (fort_is_eulerian_elastic_model(elastic_viscosity(im+1), &
             fort_viscoelastic_model(im+1)).eq.0) then 
        dest(D_DECL(i,j,k),2)=dest(D_DECL(i,j,k),1)
        dest(D_DECL(i,j,k),3)=dest(D_DECL(i,j,k),1)
        dest(D_DECL(i,j,k),4)=dest(D_DECL(i,j,k),1)
       else if (fort_is_eulerian_elastic_model(elastic_viscosity(im+1), &
                 fort_viscoelastic_model(im+1)).eq.1) then 
        if (ncomp_visc.ne.3*nmat) then
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
         if ((fort_viscoelastic_model(im+1).eq.0).or. & !FENE-CR
             (fort_viscoelastic_model(im+1).eq.1).or. & !Oldroyd-B 
             (fort_viscoelastic_model(im+1).eq.5).or. & !FENE-P 
             (fort_viscoelastic_model(im+1).eq.6)) then !linear PTT
          print *,"T11, T22, T33 must be positive"
          stop
         else if ((fort_viscoelastic_model(im+1).eq.2).or. & !displacement grad
                  (fort_viscoelastic_model(im+1).eq.3).or. & !incremental
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

        modtime=visc(D_DECL(i,j,k),2*nmat+im+1)

        if (elastic_time(im+1).eq.zero) then
         modtime=zero
        else if (elastic_time(im+1).gt.zero) then
         modtime=modtime/elastic_time(im+1)
        else
         print *,"elastic time invalid"
         stop
        endif

        if (is_rigid(nmat,im+1).eq.1) then
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
        print *,"fort_is_eulerian_elastic_model invalid"
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
       ! gravity_normalized>0 if pointing downwards
       ! 1<=gravity_dir<=dim
       ! see DRAG_COMP.H
      subroutine fort_getdrag( &
       isweep, &
       globalsum, &
       localsum, &
       gravity_normalized, &
       gravdir, &
       ntenvisco, &
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
       facevisc_index, &
       faceheat_index, &
       ncphys, &
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
       ntensor, &
       nmat, &
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

      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)

      INTEGER_T, intent(in) :: ntensor

      INTEGER_T, intent(in) :: facevisc_index
      INTEGER_T, intent(in) :: faceheat_index
      INTEGER_T, intent(in) :: ncphys

      INTEGER_T, intent(in) :: isweep
      REAL_T, intent(in) :: globalsum(N_DRAG)
      REAL_T, intent(inout) :: localsum(N_DRAG)
      REAL_T, intent(in) :: gravity_normalized
      INTEGER_T, intent(in) :: gravdir
      INTEGER_T, intent(in) :: ntenvisco
      REAL_T, intent(in) :: visc_coef
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: rzflag
      REAL_T, intent(in) :: time
      REAL_T presmag
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: bc(SDIM,2,SDIM)
      INTEGER_T, intent(in) :: DIMDEC(tdata)
      INTEGER_T, intent(in) :: DIMDEC(viscoten)
      INTEGER_T, intent(in) :: DIMDEC(den)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(slrecon)
      INTEGER_T, intent(in) :: DIMDEC(levelpc)
      INTEGER_T, intent(in) :: DIMDEC(vol)
      INTEGER_T, intent(in) :: DIMDEC(areax)
      INTEGER_T, intent(in) :: DIMDEC(areay)
      INTEGER_T, intent(in) :: DIMDEC(areaz)

      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)

      INTEGER_T, intent(in) :: DIMDEC(cvisc)
      INTEGER_T, intent(in) :: DIMDEC(c_mat_visc)

      INTEGER_T, intent(in) :: DIMDEC(solxfab)
      INTEGER_T, intent(in) :: DIMDEC(solyfab)
      INTEGER_T, intent(in) :: DIMDEC(solzfab)

      INTEGER_T, intent(in) :: DIMDEC(pres)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(drag)

      REAL_T, intent(in),target :: tdata(DIMV(tdata),ntensor)
      REAL_T, pointer :: tdata_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: viscoten(DIMV(viscoten),ntenvisco)
      REAL_T, pointer :: viscoten_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: den(DIMV(den),nmat*num_state_material)
      REAL_T, pointer :: den_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      REAL_T, intent(in),target :: slrecon(DIMV(slrecon),nmat*ngeom_recon)
      REAL_T, pointer :: slrecon_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: levelpc(DIMV(levelpc),nmat*(SDIM+1))
      REAL_T, pointer :: levelpc_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: vol(DIMV(vol))
      REAL_T, pointer :: vol_ptr(D_DECL(:,:,:))
      REAL_T, intent(in),target :: areax(DIMV(areax))
      REAL_T, pointer :: areax_ptr(D_DECL(:,:,:))
      REAL_T, intent(in),target :: areay(DIMV(areay))
      REAL_T, pointer :: areay_ptr(D_DECL(:,:,:))
      REAL_T, intent(in),target :: areaz(DIMV(areaz))
      REAL_T, pointer :: areaz_ptr(D_DECL(:,:,:))

      REAL_T, intent(in),target :: xface(DIMV(xface),ncphys)
      REAL_T, pointer :: xface_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: yface(DIMV(yface),ncphys)
      REAL_T, pointer :: yface_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: zface(DIMV(zface),ncphys)
      REAL_T, pointer :: zface_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(in),target :: cvisc(DIMV(cvisc))
      REAL_T, pointer :: cvisc_ptr(D_DECL(:,:,:))

        ! c_mat_visc initialized in fort_derviscosity
        ! 1..nmat viscosity
        ! nmat+1 ... 2*nmat viscoelastic coeff.
        ! 2*nmat+1 ... 3*nmat modtime
      REAL_T, intent(in),target :: c_mat_visc(DIMV(c_mat_visc),3*nmat)
      REAL_T, pointer :: c_mat_visc_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(in),target :: solxfab(DIMV(solxfab),nparts_def*SDIM)
      REAL_T, pointer :: solxfab_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: solyfab(DIMV(solyfab),nparts_def*SDIM)
      REAL_T, pointer :: solyfab_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in),target :: solzfab(DIMV(solzfab),nparts_def*SDIM)
      REAL_T, pointer :: solzfab_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(in),target :: pres(DIMV(pres))
      REAL_T, pointer :: pres_ptr(D_DECL(:,:,:))
      REAL_T, intent(in),target :: vel(DIMV(vel),SDIM)
      REAL_T, pointer :: vel_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(inout), target :: drag(DIMV(drag),N_DRAG)
      REAL_T, pointer :: drag_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T ii,jj,kk
      INTEGER_T imac,jmac,kmac
      INTEGER_T i_face,j_face,k_face
      INTEGER_T i_side,j_side,k_side
      INTEGER_T ii_visc,jj_visc,kk_visc
      INTEGER_T dir
      INTEGER_T veldir
      INTEGER_T facedir
      INTEGER_T side_cell  ! 0 or 1
      INTEGER_T side_visc
      INTEGER_T dirend
      INTEGER_T dir_visc
      INTEGER_T i1,j1
      INTEGER_T vofcomp,dencomp,viscbase
      INTEGER_T icell,jcell,kcell
      REAL_T vel6point(SDIM,2,SDIM)
      REAL_T ls_visc(nmat)
      REAL_T ls_side(nmat)
      REAL_T pressure_load(3)
      REAL_T viscous_stress_load(3)
      REAL_T viscous0_stress_load(3)
      REAL_T visco_stress_load(3)

      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_face(-3:3,SDIM)
      INTEGER_T nhalf

      REAL_T volume
      REAL_T gradu(3,3)

      REAL_T facearea
      INTEGER_T im
      INTEGER_T im_visc
      INTEGER_T im_test
      INTEGER_T im_primary
      INTEGER_T im_side(2)
      INTEGER_T equal_opposite_force
      REAL_T volgrid,mass
      REAL_T cengrid(SDIM)
      REAL_T global_centroid(SDIM)
      REAL_T rvec(3),gravvector(3),rcross(3)
      REAL_T pressure_rcross(3)
      REAL_T viscous0_rcross(3)
      REAL_T viscous_rcross(3)
      REAL_T visco_rcross(3)
      REAL_T grav_localtorque(SDIM)
      REAL_T pressure_localtorque(SDIM)
      REAL_T viscous_localtorque(SDIM)
      REAL_T viscous0_localtorque(SDIM)
      REAL_T visco_localtorque(SDIM)
      REAL_T Q(3,3)
      REAL_T nsolid(3)
      INTEGER_T mask_cell
      REAL_T ls_sort(nmat)
      REAL_T mu_0 ! viscosity for ambient case.
      REAL_T mu_non_ambient ! viscosity for non ambient case.
      REAL_T delx
      INTEGER_T partid
      INTEGER_T ibase

      nhalf=3

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

      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid fort_getdrag"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid fort_getdrag"
       stop
      endif

        ! cell centered grad U
      call checkbound_array(fablo,fabhi,tdata_ptr,0,-1,1252)
      call checkbound_array(fablo,fabhi,viscoten_ptr,1,-1,1253)
      call checkbound_array(fablo,fabhi,den_ptr,1,-1,1254)
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1,1255)
      call checkbound_array(fablo,fabhi,slrecon_ptr,1,-1,12560)
      call checkbound_array(fablo,fabhi,levelpc_ptr,1,-1,1257)
      call checkbound_array1(fablo,fabhi,vol_ptr,0,-1,6600)
      call checkbound_array1(fablo,fabhi,areax_ptr,0,0,6601)
      call checkbound_array1(fablo,fabhi,areay_ptr,0,1,6602)
      call checkbound_array1(fablo,fabhi,areaz_ptr,0,SDIM-1,6603)

      call checkbound_array(fablo,fabhi,xface_ptr,0,0,1258)
      call checkbound_array(fablo,fabhi,yface_ptr,0,1,1259)
      call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1,1261)

      call checkbound_array1(fablo,fabhi,cvisc_ptr,0,-1,1262)
      call checkbound_array(fablo,fabhi,c_mat_visc_ptr,1,-1,1262)

      call checkbound_array(fablo,fabhi,solxfab_ptr,0,0,6604)
      call checkbound_array(fablo,fabhi,solyfab_ptr,0,1,6604)
      call checkbound_array(fablo,fabhi,solzfab_ptr,0,SDIM-1,6604)

      call checkbound_array1(fablo,fabhi,pres_ptr,1,-1,6605)
      call checkbound_array(fablo,fabhi,vel_ptr,1,-1,6606)
      call checkbound_array(fablo,fabhi,drag_ptr,0,-1,6607)

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
      do im=1,nmat
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
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (rzflag.eq.0) then
       ! do nothing
      else if (rzflag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rzflag.eq.3) then
       ! do nothing
      else
       print *,"rzflag invalid"
       stop
      endif

      if (ntensor.ne.SDIM*SDIM) then
       print *,"ntensor invalid"
       stop
      endif

      if ((isweep.lt.0).or.(isweep.gt.1)) then
       print *,"isweep invalid"
       stop
      endif
      if ((gravdir.lt.1).or.(gravdir.gt.SDIM)) then
       print *,"gravdir invalid"
       stop
      endif
      if (FORT_NUM_TENSOR_TYPE.ne.2*SDIM) then
       print *,"FORT_NUM_TENSOR_TYPE invalid"
       stop
      endif
      if ((num_materials_viscoelastic.ge.1).and. &
          (num_materials_viscoelastic.le.nmat)) then
       if (FORT_NUM_TENSOR_TYPE*num_materials_viscoelastic.ne.ntenvisco) then
        print *,"ntenvisco invalid1"
        stop
       endif
      else if (num_materials_viscoelastic.eq.0) then
       ! do nothing
      else
       print *,"num_materials_viscoelastic invalid: fort_getdrag"
       stop
      endif
       ! indexes start at 0
      if ((facevisc_index.ne.6).or. &
          (faceheat_index.ne.7)) then
       print *,"face_index bust 10"
       stop
      endif
      if (ncphys.lt.8) then
       print *,"ncphys invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       ! first sweep - find the mass and centroid of materials
      if (isweep.eq.0) then
 
       do icell=growlo(1),growhi(1)
       do jcell=growlo(2),growhi(2)
       do kcell=growlo(3),growhi(3)
        do dir=1,N_DRAG
         drag(D_DECL(icell,jcell,kcell),dir)=zero
        enddo
        do im_test=1,nmat
         mask_cell=NINT(mask(D_DECL(icell,jcell,kcell)))
         if (mask_cell.eq.1) then
          call gridsten(xsten,xlo,icell,jcell,kcell,fablo,bfact,dx,nhalf)
          call Box_volumeFAST(bfact,dx,xsten,nhalf, &
           volgrid,cengrid,SDIM)
          vofcomp=(im_test-1)*ngeom_recon+1
          dencomp=(im_test-1)*num_state_material+1
          mass=den(D_DECL(icell,jcell,kcell),dencomp)*volgrid* &
           slrecon(D_DECL(icell,jcell,kcell),vofcomp)
          localsum(DRAGCOMP_MASS+im_test)= &
             localsum(DRAGCOMP_MASS+im_test)+mass
          do dir=1,SDIM
           localsum(DRAGCOMP_COM+3*(im_test-1)+dir)= &
            localsum(DRAGCOMP_COM+3*(im_test-1)+dir)+ &
             mass* &
             (slrecon(D_DECL(icell,jcell,kcell),vofcomp+dir)+cengrid(dir))
          enddo ! dir=1..sdim
         else if (mask_cell.eq.0) then
          ! do nothing
         else
          print *,"mask_cell invalid"
          stop
         endif 
        enddo ! im_test=1..nmat

       enddo
       enddo
       enddo  ! isweep=0 centroids

      else if (isweep.eq.1) then ! above, mass and centroid

       do im_test=1,nmat

        mass=globalsum(DRAGCOMP_MASS+im_test)
        if (mass.lt.zero) then
         print *,"mass cannot be negative  im_test,mass=",im_test,mass
         stop
        else if (mass.eq.zero) then
         print *,"WARNING: mass=0"
         print *,"im_test=",im_test
         print *,"FSI_flag(im_test)=",FSI_flag(im_test)
         do dir=1,SDIM
          global_centroid(dir)=zero
         enddo
        else if (mass.gt.zero) then
         do dir=1,SDIM
          global_centroid(dir)=globalsum(DRAGCOMP_COM+3*(im_test-1)+dir)/mass
         enddo
         if (levelrz.eq.0) then
          ! do nothing
         else if (levelrz.eq.1) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          global_centroid(1)=zero
         else if (levelrz.eq.3) then
          ! do nothing
         else
          print *,"levelrz invalid fort_getdrag"
          stop
         endif
        else
         print *,"mass bust"
         stop
        endif
       enddo ! im_test=1..nmat
    
        ! force=integral body forces + integral_boundary tau dot n dA
        !  tau=-pI + 2 mu D + mu_p f(A)/lambda  \tilde{Q} 
        ! buoyancy force (body forces within the materials).
        ! Also, update the moment of inertia integral.
       do icell=growlo(1),growhi(1)
       do jcell=growlo(2),growhi(2)
       do kcell=growlo(3),growhi(3)

         ! calculate the forces exerted on material "im_test"
        do im_test=1,nmat
         mask_cell=NINT(mask(D_DECL(icell,jcell,kcell)))
         if (mask_cell.eq.1) then
          call gridsten(xsten,xlo,icell,jcell,kcell,fablo,bfact,dx,nhalf)
          call Box_volumeFAST(bfact,dx,xsten,nhalf, &
           volgrid,cengrid,SDIM)
          vofcomp=(im_test-1)*ngeom_recon+1
          dencomp=(im_test-1)*num_state_material+1
          mass=den(D_DECL(icell,jcell,kcell),dencomp)*volgrid* &
           slrecon(D_DECL(icell,jcell,kcell),vofcomp)
          do dir=1,SDIM
           gravvector(dir)=zero
           rvec(dir)=slrecon(D_DECL(icell,jcell,kcell),vofcomp+dir)+ &
            cengrid(dir)-global_centroid(dir)
          enddo
          gravvector(gravdir)=-mass*gravity_normalized
          if (SDIM.eq.2) then
           rvec(3)=zero
           gravvector(3)=zero
          endif 

          do dir=1,SDIM
           ibase=DRAGCOMP_FORCE+3*(im_test-1)+dir
           localsum(ibase)=localsum(ibase)+gravvector(dir)
           drag(D_DECL(icell,jcell,kcell),ibase)= &
             drag(D_DECL(icell,jcell,kcell),ibase)+gravvector(dir)
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

          ibase=DRAGCOMP_MOMINERTIA+3*(im_test-1)
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
 
          ibase=DRAGCOMP_TORQUE+3*(im_test-1)
          do dir=1,dirend
           localsum(ibase+dir)=localsum(ibase+dir)+grav_localtorque(dir)
           drag(D_DECL(icell,jcell,kcell),ibase+dir)= &
             drag(D_DECL(icell,jcell,kcell),ibase+dir)+grav_localtorque(dir)
          enddo

         else if (mask_cell.eq.0) then
          ! do nothing
         else
          print *,"mask_cell invalid"
          stop
         endif 

        enddo ! im_test=1..nmat

       enddo
       enddo
       enddo  ! icell,jcell,kcell (cells - body forces and torques (buoancy) )

        ! cells - pressure, viscosity, viscoelastic
       do icell=growlo(1),growhi(1)
       do jcell=growlo(2),growhi(2)
       do kcell=growlo(3),growhi(3)

        mask_cell=NINT(mask(D_DECL(icell,jcell,kcell)))

        if (mask_cell.eq.1) then

          ! im_test is the material for which a given force is applied.
         do im_test=1,nmat

          do im=1,nmat
           ls_sort(im)=levelpc(D_DECL(icell,jcell,kcell),im)
          enddo
           ! declared in GLOBALUTIL.F90
          call get_primary_material(ls_sort,nmat,im_primary)

          if (im_primary.eq.im_test) then
           ! do nothing
          else if (im_primary.ne.im_test) then

           call gridsten(xsten,xlo, &
            icell,jcell,kcell,fablo,bfact,dx,nhalf)

           volume=vol(D_DECL(icell,jcell,kcell))
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
             do im=1,nmat
              ls_side(im)=levelpc(D_DECL(i_side,j_side,k_side),im)
             enddo
             call get_primary_material(ls_side,nmat,im_side(side_cell+1))
            enddo !side_cell=0..1
           
            equal_opposite_force=1

            if ((im_side(1).eq.im_test).and. &
                (im_side(2).eq.im_test)) then
             ! do nothing - both adjoining cells in the forced material
            else if ((im_side(1).ne.im_test).and. &
                     (im_side(2).ne.im_test)) then
             ! do nothing - both adjoining cells not in forced material.
            else if ( ((im_side(1).eq.im_test).and. &
                       (im_side(2).ne.im_test)).or. &
                      ((im_side(1).ne.im_test).and. &
                       (im_side(2).eq.im_test)) ) then
             equal_opposite_force=0
            else
             print *,"im_side became corrupt"
             stop
            endif

            if (equal_opposite_force.eq.0) then

             if (im_side(1).eq.im_test) then
              side_cell=0
             else if (im_side(2).eq.im_test) then
              side_cell=1
             else
              print *,"im_side became corrupt"
              stop
             endif

             imac=icell+side_cell*ii
             jmac=jcell+side_cell*jj
             kmac=kcell+side_cell*kk

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

             ibase=DRAGCOMP_PERIM_VECTOR+3*(im_test-1)
             drag(D_DECL(icell,jcell,kcell),ibase+facedir)= &
                     facearea*nsolid(facedir)

             localsum(ibase+facedir)=localsum(ibase+facedir)+ &
               facearea*nsolid(facedir)

              !facedir=1..sdim
             call gridstenMAC(xsten_face,xlo,imac,jmac,kmac,fablo,bfact, &
               dx,nhalf,facedir-1,91)

             ! im_primary is the forcing fluid at cell (icell,jcell,kcell)
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
                i_side=icell-ii_visc
                j_side=jcell-jj_visc
                k_side=kcell-kk_visc
                i_face=icell
                j_face=jcell
                k_face=kcell
               else if (side_visc.eq.2) then
                i_side=icell+ii_visc
                j_side=jcell+jj_visc
                k_side=kcell+kk_visc
                i_face=icell+ii_visc
                j_face=jcell+jj_visc
                k_face=kcell+kk_visc
               else
                print *,"side_visc invalid"
                stop
               endif
              
               do im=1,nmat
                ls_visc(im)=levelpc(D_DECL(i_side,j_side,k_side),im)
               enddo
               call get_primary_material(ls_visc,nmat,im_visc)
               if (is_rigid(nmat,im_visc).eq.1) then

                partid=0
                do im=1,im_visc-1
                 if (is_lag_part(nmat,im).eq.1) then
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

               else if (is_rigid(nmat,im_visc).eq.0) then
                do dir=1,SDIM
                 vel6point(dir_visc,side_visc,dir)= &
                   half*(vel(D_DECL(icell,jcell,kcell),dir)+ &
                         vel(D_DECL(i_side,j_side,k_side),dir))
                enddo
               else
                print *,"is_rigid(nmat,im_visc) invalid"
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

             if (fort_is_eulerian_elastic_model( &
                   fort_elastic_viscosity(im_primary), &
                   fort_viscoelastic_model(im_primary)).eq.1) then 
              partid=1
              do while ((fort_im_elastic_map(partid)+1.ne.im_primary).and. &
                        (partid.le.num_materials_viscoelastic))
               partid=partid+1
              enddo
              if (partid.le.num_materials_viscoelastic) then
               viscbase=(partid-1)*FORT_NUM_TENSOR_TYPE
               Q(1,1)=viscoten(D_DECL(icell,jcell,kcell),viscbase+1)
               Q(1,2)=viscoten(D_DECL(icell,jcell,kcell),viscbase+2)
               Q(2,2)=viscoten(D_DECL(icell,jcell,kcell),viscbase+3)
               Q(3,3)=viscoten(D_DECL(icell,jcell,kcell),viscbase+4)
               Q(1,3)=zero
               Q(2,3)=zero
#if (AMREX_SPACEDIM==3)
               Q(1,3)=viscoten(D_DECL(icell,jcell,kcell),viscbase+5)
               Q(2,3)=viscoten(D_DECL(icell,jcell,kcell),viscbase+6)
#endif
               Q(2,1)=Q(1,2)
               Q(3,1)=Q(1,3)
               Q(3,2)=Q(2,3)
              else
               print *,"partid invalid in fort_getdrag"
               stop
              endif
             else if (fort_is_eulerian_elastic_model( &
                       fort_elastic_viscosity(im_primary), &
                       fort_viscoelastic_model(im_primary)).eq.0) then 
              ! do nothing
             else
              print *,"fort_is_eulerian_elastic_model invalid"
              stop
             endif

             mu_0=get_user_viscconst(im_primary, &
               fort_denconst(im_primary),fort_tempconst(im_primary))
              ! c_mat_visc is initialized in NavierStokes2.cpp: getStateVISC,
              ! getStateVISC_ALL.  "c_mat_visc" includes WALE model and
              ! "viscconst_eddy" effects.
             mu_non_ambient= &
               c_mat_visc(D_DECL(icell,jcell,kcell),im_primary)
             if ((mu_0.ge.zero).and.(mu_non_ambient.ge.zero)) then
              ! do nothing
             else
              print *,"mu_0 or mu_non_ambient invalid"
              stop
             endif

             do j1=1,3
              viscous_stress_load(j1)=zero
              viscous0_stress_load(j1)=zero
              visco_stress_load(j1)=zero
              pressure_load(j1)=zero

              do i1=1,3
               viscous0_stress_load(j1)= &
                 viscous0_stress_load(j1)- &
                   mu_0*visc_coef* &
                   (gradu(i1,j1)+gradu(j1,i1))*nsolid(i1)
               viscous_stress_load(j1)= &
                 viscous_stress_load(j1)- &
                   mu_non_ambient*visc_coef* &
                   (gradu(i1,j1)+gradu(j1,i1))*nsolid(i1)
               visco_stress_load(j1)= &
                 visco_stress_load(j1)- &
                  Q(i1,j1)*nsolid(i1)
              enddo ! i1=1,3

              pressure_load(j1)=presmag*nsolid(j1)*facearea
              viscous0_stress_load(j1)=viscous0_stress_load(j1)*facearea
              viscous_stress_load(j1)=viscous_stress_load(j1)*facearea
              visco_stress_load(j1)=visco_stress_load(j1)*facearea
             enddo ! j1=1,3

! global_centroid will be incorrect if the solid geometry is reflected
! across a domain boundary.

             do dir=1,SDIM
              rvec(dir)=xsten_face(0,dir)-global_centroid(dir)
             enddo
             if (SDIM.eq.2) then
              rvec(3)=zero
              viscous_stress_load(3)=zero
              viscous0_stress_load(3)=zero
              visco_stress_load(3)=zero
              pressure_load(3)=zero
             endif

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
              visco_localtorque(SDIM)=visco_rcross(2) ! x-z plane of rotation
              pressure_localtorque(SDIM)=pressure_rcross(2) !x-z plane of rotation
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

             do dir=1,SDIM
              ibase=DRAGCOMP_FORCE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               pressure_load(dir)+ &
               viscous_stress_load(dir)+ &
               visco_stress_load(dir)

              localsum(ibase)=localsum(ibase)+ &
               pressure_load(dir)+ &
               viscous_stress_load(dir)+ &
               visco_stress_load(dir)

              ibase=DRAGCOMP_PFORCE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               pressure_load(dir)

              localsum(ibase)=localsum(ibase)+ &
               pressure_load(dir)

              ibase=DRAGCOMP_VISCOUSFORCE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               viscous_stress_load(dir)

              localsum(ibase)=localsum(ibase)+ &
               viscous_stress_load(dir)

              ibase=DRAGCOMP_VISCOUS0FORCE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               viscous0_stress_load(dir)

              localsum(ibase)=localsum(ibase)+ &
               viscous0_stress_load(dir)

              ibase=DRAGCOMP_VISCOFORCE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               visco_stress_load(dir)

              localsum(ibase)=localsum(ibase)+ &
               visco_stress_load(dir)

             enddo ! dir=1..sdim

             dirend=1
             if (SDIM.eq.3) then
              dirend=3
             endif

             do dir=1,dirend

              ibase=DRAGCOMP_TORQUE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               pressure_localtorque(dir)+ &
               viscous_localtorque(dir)+ &
               visco_localtorque(dir)

              localsum(ibase)=localsum(ibase)+ &
               pressure_localtorque(dir)+ &
               viscous_localtorque(dir)+ &
               visco_localtorque(dir)

              ibase=DRAGCOMP_PTORQUE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               pressure_localtorque(dir)

              localsum(ibase)=localsum(ibase)+ &
               pressure_localtorque(dir)

              ibase=DRAGCOMP_VISCOUSTORQUE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               viscous_localtorque(dir)

              localsum(ibase)=localsum(ibase)+ &
               viscous_localtorque(dir)

              ibase=DRAGCOMP_VISCOUS0TORQUE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               viscous0_localtorque(dir)

              localsum(ibase)=localsum(ibase)+ &
               viscous0_localtorque(dir)

              ibase=DRAGCOMP_VISCOTORQUE+3*(im_test-1)+dir

              drag(D_DECL(icell,jcell,kcell),ibase)= &
               drag(D_DECL(icell,jcell,kcell),ibase)+ &
               visco_localtorque(dir)

              localsum(ibase)=localsum(ibase)+ &
               visco_localtorque(dir)

             enddo ! dir=1..dirend

             ibase=DRAGCOMP_PERIM+im_test
             drag(D_DECL(icell,jcell,kcell),ibase)= &
              drag(D_DECL(icell,jcell,kcell),ibase)+facearea

             localsum(ibase)=localsum(ibase)+facearea

            else if (equal_opposite_force.eq.1) then
             ! do nothing
            else
             print *,"equal_opposite_force became corrupt"
              stop
            endif

           enddo ! facedir=1..sdim
          else
           print *,"im_primary or im_test corruption"
           stop
          endif

         enddo ! im_test=1..nmat

        else if (mask_cell.eq.0) then
         ! do nothing
        else
         print *,"mask_cell invalid"
         stop
        endif

       enddo
       enddo
       enddo  ! icell,jcell,kcell

      else
       print *,"isweep invalid"
       stop
      endif

      return
      end subroutine fort_getdrag

      end module derive_module

      subroutine FORT_INTEGRATE_RECALESCE( &
       isweep, &
       globalsum, &
       localsum, &
       recalesce_material_in, &
       snew,DIMS(snew), & 
       mask,DIMS(mask), & 
       vol,DIMS(vol), &
       xlo,dx, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       velbc, &
       time, &
       num_integrate,nmat,ncomp_state, &
       level,finest_level)

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probf90_module
      use godunov_module


      IMPLICIT NONE

      INTEGER_T isweep,num_integrate,nmat
      INTEGER_T ncomp_state,level,finest_level
      INTEGER_T ncomp_state_test
      INTEGER_T recalesce_material_in(nmat)
      REAL_T globalsum(num_integrate*nmat)
      REAL_T localsum(num_integrate*nmat)
      REAL_T time
      INTEGER_T tilelo(SDIM), tilehi(SDIM)
      INTEGER_T fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact
      INTEGER_T velbc(SDIM,2,SDIM)
      INTEGER_T DIMDEC(mask)
      INTEGER_T DIMDEC(snew)
      INTEGER_T DIMDEC(vol)

      REAL_T mask(DIMV(mask))
      REAL_T snew(DIMV(snew),ncomp_state)
      REAL_T vol(DIMV(vol))
      REAL_T xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k,im
      INTEGER_T ibase
      INTEGER_T dir
      INTEGER_T dencomp,vofcomp,nhalf
      REAL_T xsten(-1:1,SDIM)
      REAL_T xsten_center(SDIM)
      REAL_T cengrid(SDIM)
      INTEGER_T mask_cell
      REAL_T volgrid,mass,centroid_term,vfrac
      REAL_T dxmax,dxmaxLS,dist_substrate,temperature
      REAL_T mu
      INTEGER_T im_ice

      call checkbound(fablo,fabhi,DIMS(mask),0,-1,1255)
      call checkbound(fablo,fabhi,DIMS(snew),0,-1,12561)
      call checkbound(fablo,fabhi,DIMS(vol),0,-1,6608)

      nhalf=1

      if (num_integrate.ne.5+SDIM+1) then
       print *,"num_integrate invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      ncomp_state_test=(SDIM+1)+ &
       nmat*(num_state_material+ngeom_raw)+1
      if (ncomp_state.ne.ncomp_state_test) then
       print *,"ncomp_state invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid integrate recalesce"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid6"
       stop
      endif
      do im=1,nmat
       if (fort_denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif
       mu=get_user_viscconst(im,fort_denconst(im),fort_tempconst(im))
       if (mu.lt.zero) then
        print *,"viscconst invalid"
        stop
       endif
      enddo !im=1..nmat

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if ((isweep.lt.0).or.(isweep.gt.1)) then
       print *,"isweep invalid"
       stop
      endif

      im_ice=0
      do im=1,nmat
       if (is_ice(nmat,im).eq.1) then
        im_ice=im
       else if (is_ice(nmat,im).eq.0) then
        ! do nothing
       else
        print *,"is_ice invalid"
        stop
       endif
      enddo

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       ! first sweep - find integral of temperature, dist_farthest,
       ! dist_closest, touch_flag, centroid, volume
       ! 2nd sweep = find temperature_farthest
       ! all the vars:
       ! average temperature, temperature_farthest, dist_farthest,
       ! dist_closest, touch_flag, xcentroid, mass

      if (isweep.eq.0) then
 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
       
        mask_cell=NINT(mask(D_DECL(i,j,k)))
        if (mask_cell.eq.1) then
         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
         do dir=1,SDIM
          xsten_center(dir)=xsten(0,dir)
         enddo
         call Box_volumeFAST(bfact,dx,xsten,nhalf, &
           volgrid,cengrid,SDIM)

         do im=1,nmat
          if ((recalesce_material_in(im).eq.1).or. &
              (recalesce_material_in(im).eq.2)) then
           ibase=(im-1)*num_integrate

             ! dist_substrate>0 in the substrate
           call icemask_override( &
             xsten_center, &
             im,im_ice,dist_substrate)
             ! make dist_substrate>0 outside the substrate
           dist_substrate=-dist_substrate

           vofcomp=(SDIM+1)+ &
             nmat*num_state_material+(im-1)*ngeom_raw+1
           dencomp=(SDIM+1)+ &
             (im-1)*num_state_material+1
           temperature=snew(D_DECL(i,j,k),dencomp+1) 
           vfrac=snew(D_DECL(i,j,k),vofcomp)
           if (vfrac.ge.half) then
            if (temperature.lt.zero) then
             print *,"temperature invalid"
             stop
            endif
            if (dist_substrate.le.0.75*dxmaxLS) then
             localsum(ibase+5)=one
            endif
             ! find furthest distance from the substrate
            if (dist_substrate.gt.localsum(ibase+3)) then
             localsum(ibase+3)=dist_substrate
            endif
             ! find closest distance to the substrate
            if (dist_substrate.lt.localsum(ibase+4)) then
             localsum(ibase+4)=dist_substrate
            endif
           endif

           mass=snew(D_DECL(i,j,k),dencomp)*volgrid*vfrac
           localsum(ibase+1)=localsum(ibase+1)+mass*temperature
           localsum(ibase+SDIM+6)=localsum(ibase+SDIM+6)+mass
           do dir=1,SDIM
            centroid_term=mass* &
             (snew(D_DECL(i,j,k),vofcomp+dir)+cengrid(dir))
            if (dir.eq.1) then
             if (levelrz.eq.0) then
              ! do nothing
             else if (levelrz.eq.1) then
              if (SDIM.ne.2) then
               print *,"dimension bust"
               stop
              endif
              centroid_term=zero
             else if (levelrz.eq.3) then
              ! do nothing
             else
              print *,"levelrz invalid integrate recalesce"
              stop
             endif
            endif ! dir=1?

            localsum(ibase+5+dir)=localsum(ibase+5+dir)+centroid_term
           enddo ! dir=1 ... sdim
          else if (recalesce_material_in(im).eq.0) then
           ! do nothing
          else
           print *,"recalesce_material_in invalid"
           stop
          endif
         enddo ! im 
        else if (mask_cell.eq.0) then
         ! do nothing
        else
         print *,"mask_cell invalid"
         stop
        endif 

       enddo
       enddo
       enddo  ! isweep=0 everything except for temperature farthest.

      else if (isweep.eq.1) then

       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
          
        mask_cell=NINT(mask(D_DECL(i,j,k)))
        if (mask_cell.eq.1) then
         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
         do dir=1,SDIM
          xsten_center(dir)=xsten(0,dir)
         enddo
         call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volgrid,cengrid,SDIM)
             
         do im=1,nmat
          if ((recalesce_material_in(im).eq.1).or. &
              (recalesce_material_in(im).eq.2)) then
           ibase=(im-1)*num_integrate
            
            ! dist_substrate>0 in the substrate
           call icemask_override( &
             xsten_center, &
             im,im_ice,dist_substrate)
            ! make dist_substrate>0 outside the substrate
           dist_substrate=-dist_substrate

           vofcomp=(SDIM+1)+ &
             nmat*num_state_material+(im-1)*ngeom_raw+1
           dencomp=(SDIM+1)+ &
             (im-1)*num_state_material+1
           temperature=snew(D_DECL(i,j,k),dencomp+1) 
           vfrac=snew(D_DECL(i,j,k),vofcomp)
           if (vfrac.ge.half) then
            if (temperature.lt.zero) then
             print *,"temperature invalid"
             stop
            endif
             ! temperature furthest from the substrate.
            if (dist_substrate.ge.globalsum(ibase+3)*(one-VOFTOL)) then
             localsum(ibase+2)=temperature
            endif
           endif
          else if (recalesce_material_in(im).eq.0) then
           ! do nothing
          else
           print *,"recalesce_material_in invalid"
           stop
          endif
         enddo ! im 
        else if (mask_cell.eq.0) then
         ! do nothing
        else
         print *,"mask_cell invalid"
         stop
        endif 

       enddo
       enddo
       enddo  ! isweep=1 temperature farthest.

      else
       print *,"isweep invalid"
       stop
      endif

      return
      end subroutine FORT_INTEGRATE_RECALESCE



      subroutine FORT_RESET_TEMPERATURE( &
       im_source, &
       TSAT, &
       snew,DIMS(snew), & 
       xlo,dx, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       velbc, &
       time, &
       nmat,ncomp_state, &
       level,finest_level)

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probf90_module
      use godunov_module


      IMPLICIT NONE

      INTEGER_T im_source,nmat
      INTEGER_T ncomp_state,level,finest_level
      INTEGER_T ncomp_state_test
      REAL_T time,TSAT
      INTEGER_T tilelo(SDIM), tilehi(SDIM)
      INTEGER_T fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact
      INTEGER_T velbc(SDIM,2,SDIM)
      INTEGER_T DIMDEC(snew)

      REAL_T snew(DIMV(snew),ncomp_state)
      REAL_T xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k,im
      INTEGER_T dencomp,vofcomp
      REAL_T vfrac
      REAL_T mu

      call checkbound(fablo,fabhi,DIMS(snew),0,-1,12562)

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      ncomp_state_test=(SDIM+1)+ &
       nmat*(num_state_material+ngeom_raw)+1
      if (ncomp_state.ne.ncomp_state_test) then
       print *,"ncomp_state invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid reset_temperature"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid7"
       stop
      endif
      do im=1,nmat
       if (fort_denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif
       mu=get_user_viscconst(im,fort_denconst(im),fort_tempconst(im))
       if (mu.lt.zero) then
        print *,"viscconst invalid"
        stop
       endif
      enddo

      if ((im_source.lt.0).or.(im_source.ge.nmat)) then
       print *,"im_source invalid"
       stop
      endif
      if (TSAT.le.zero) then
       print *,"TSAT invalid in fort_reset_temperature"
       stop
      endif
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       
       vofcomp=(SDIM+1)+ &
          nmat*num_state_material+im_source*ngeom_raw+1
       vfrac=snew(D_DECL(i,j,k),vofcomp)
       if (vfrac.ge.half) then
        do im=1,nmat
         dencomp=(SDIM+1)+ &
           (im-1)*num_state_material+1
         snew(D_DECL(i,j,k),dencomp+1)=TSAT
        enddo
       endif
      enddo
      enddo
      enddo  

      return
      end subroutine FORT_RESET_TEMPERATURE


      subroutine FORT_MAXPRESVEL( &
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
        fablo,fabhi,bfact)

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      REAL_T, intent(inout) :: minpres,maxpres,maxvel,maxvel_collide

      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(velx)
      INTEGER_T, intent(in) :: DIMDEC(vely)
      INTEGER_T, intent(in) :: DIMDEC(velz)

      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: vel(DIMV(vel),(SDIM+1))
      REAL_T, intent(in) :: velx(DIMV(velx))
      REAL_T, intent(in) :: vely(DIMV(vely))
      REAL_T, intent(in) :: velz(DIMV(velz))
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T dir
      REAL_T magvel,magvel_mac,magvel_collide
      REAL_T vello,velhi,vellohi
      INTEGER_T i,j,k

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(mask),0,-1,6609)
      call checkbound(fablo,fabhi,DIMS(vel),0,-1,6610)
      call checkbound(fablo,fabhi,DIMS(velx),0,0,6610)
      call checkbound(fablo,fabhi,DIMS(vely),0,1,6610)
      call checkbound(fablo,fabhi,DIMS(velz),0,SDIM-1,6610)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
    
       if (mask(D_DECL(i,j,k)).eq.one) then

        if (vel(D_DECL(i,j,k),SDIM+1).gt.maxpres) then
          maxpres=vel(D_DECL(i,j,k),SDIM+1)
        endif
        if (vel(D_DECL(i,j,k),SDIM+1).lt.minpres) then
          minpres=vel(D_DECL(i,j,k),SDIM+1)
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
      end subroutine FORT_MAXPRESVEL


