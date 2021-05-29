#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

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


      end module derive_module

       ! WALE model
       ! Wall Adapting Local Eddy-viscosity models for simulations in complex
       ! geometries.   Ducros, Nicoud, Poinsot 1998
       ! called from getStateVISC
      subroutine FORT_DERTURBVISC( &
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
       ncompvisc)

      use global_utility_module
      use probf90_module
      use derive_module

      IMPLICIT NONE

      INTEGER_T level
      INTEGER_T im
      INTEGER_T nmat
      INTEGER_T ncompvisc
      INTEGER_T ngrow
      INTEGER_T ntensor
      REAL_T dt
      REAL_T cur_time

      INTEGER_T tilelo(SDIM), tilehi(SDIM)
      INTEGER_T fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact
      REAL_T dx(SDIM),xlo(SDIM)

      INTEGER_T DIMDEC(denstate)
      INTEGER_T DIMDEC(vof)
      INTEGER_T DIMDEC(vel)
      INTEGER_T DIMDEC(visc)
      INTEGER_T DIMDEC(cellten)
  
      REAL_T denstate(DIMV(denstate),nmat*num_state_material)
      REAL_T vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T vel(DIMV(vel),SDIM)
      REAL_T visc(DIMV(visc),ncompvisc)
      REAL_T cellten(DIMV(cellten),ntensor)

      REAL_T g(3,3),s(3,3),sd(3,3),g2(3),g2tr,ss,sdsd
      REAL_T g2_12,g2_13,g2_21,g2_23,g2_31,g2_32,turb_visc,density
      REAL_T rr
      INTEGER_T i,j,k,veldir,dir,flagcomp
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz,nbase
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf
     
      nhalf=1

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
     
      call checkbound(fablo,fabhi,DIMS(vof),ngrow+1,-1,310)
      call checkbound(fablo,fabhi,DIMS(visc),ngrow,-1,311)
      call checkbound(fablo,fabhi,DIMS(vel),ngrow,-1,312)
      call checkbound(fablo,fabhi,DIMS(denstate),ngrow,-1,313)
      call checkbound(fablo,fabhi,DIMS(cellten),ngrow,-1,314)

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
       g2tr=(g2(1)+g2(2)+g2(3))/3. ! tensor trace
       g2_12=g(1,1)*g(1,2)+g(1,2)*g(2,2)+g(1,3)*g(3,2)
       g2_13=g(1,1)*g(1,3)+g(1,2)*g(2,3)+g(1,3)*g(3,3)
       g2_21=g(2,1)*g(1,1)+g(2,2)*g(2,1)+g(2,3)*g(3,1)
       g2_23=g(2,1)*g(1,3)+g(2,2)*g(2,3)+g(2,3)*g(3,3)
       g2_31=g(3,1)*g(1,1)+g(3,2)*g(2,1)+g(3,3)*g(3,1)
       g2_32=g(3,1)*g(1,2)+g(3,2)*g(2,2)+g(3,3)*g(3,2)
 
       sd(1,1)=g2(1)-g2tr
       sd(1,2)=0.5*(g2_12+g2_21)
       sd(1,3)=0.5*(g2_13+g2_31)

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

       if(sdsd.ne.zero) then

        turb_visc=(0.5*dx(1))**2*sdsd**1.5/(ss**2.5+sdsd**1.25)
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
  
      enddo 
      enddo 
      enddo 

      return
      end subroutine DERTURBVISC

      subroutine FORT_GETSHEAR( &
       im, &
       ntensor, &
       cellten,DIMS(cellten), &
       vof,DIMS(vof), &
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
       ngrow,nmat)

      use global_utility_module
      use derive_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: im,level,ntensor
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
      INTEGER_T, intent(in) :: DIMDEC(vof)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(tensordata)
  
      REAL_T, intent(in) :: cellten(DIMV(cellten),ntensor)
      REAL_T, intent(in) :: vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T, intent(in) :: vel(DIMV(vel),SDIM)
      REAL_T, intent(out) :: tensordata(DIMV(tensordata),20)
      REAL_T visctensor(3,3),gradu(3,3)
      REAL_T shear
      REAL_T a,b,c
      REAL_T rr
      REAL_T vort(3)
      INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz,nbase

      nhalf=1

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
      if ((im.lt.0).or.(im.ge.num_materials)) then
       print *,"im invalid2"
       stop
      endif

      ! compute u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)

      call checkbound(fablo,fabhi,DIMS(cellten),ngrow,-1,64)
      call checkbound(fablo,fabhi,DIMS(vof),ngrow+1,-1,64)
      call checkbound(fablo,fabhi,DIMS(vel),ngrow+1,-1,64)
      call checkbound(fablo,fabhi, &
       DIMS(tensordata), &
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
      end subroutine FORT_GETSHEAR

       ! called from getStateVISC
      subroutine FORT_DERVISCOSITY( &
        level, &
        finest_level, &
        visc_coef, &
        im_parm, &
        nmat, &
        dt, &
        viscosity_coefficient, &
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
        vof,DIMS(vof), &
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
        ncompvisc)

      use global_utility_module
      use probf90_module
      use derive_module

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
      INTEGER_T, intent(in) :: DIMDEC(vof)
      INTEGER_T, intent(in) :: DIMDEC(visc)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(eosdata)
      INTEGER_T, intent(in) :: DIMDEC(tensor)
      INTEGER_T, intent(in) :: DIMDEC(gammadot)
      INTEGER_T, intent(in) :: bc(SDIM,2,SDIM)
      INTEGER_T, intent(in) :: ngrow
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: dx(SDIM), xlo(SDIM)
      REAL_T, intent(in) :: vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T, intent(out) :: visc(DIMV(visc),ncompvisc) !ncompvisc=3*nmat
      REAL_T, intent(in) :: vel(DIMV(vel),SDIM)
      REAL_T, intent(in) :: gammadot(DIMV(gammadot))
      REAL_T, intent(in) :: eosdata(DIMV(eosdata),nmat*num_state_material)
      REAL_T, intent(in) :: tensor(DIMV(tensor),FORT_NUM_TENSOR_TYPE)

      INTEGER_T i,j,k
      REAL_T    shear,density,temperature,mu
      INTEGER_T flagcomp
      REAL_T bterm,pterm
      INTEGER_T numstatetest,ii,jj
      REAL_T Q(3,3)
      REAL_T traceA,modtime,viscoelastic_coeff
      REAL_T bulk_modulus

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

      if ((viscoelastic_model.ge.0).and. &
          (viscoelastic_model.le.3)) then
       ! do nothing
      else
       print *,"viscoelastic_model invalid"
       stop
      endif

      if (dt.gt.zero) then 
       ! do nothing
      else
       print *,"dt invalid in FORT_DERVISCOSITY"
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
     
      call checkbound(fablo,fabhi,DIMS(vof),ngrow+1,-1,315)
      call checkbound(fablo,fabhi,DIMS(visc),ngrow,-1,316)
      call checkbound(fablo,fabhi,DIMS(gammadot),ngrow,-1,317)
      call checkbound(fablo,fabhi,DIMS(eosdata),ngrow,-1,318)
      call checkbound(fablo,fabhi,DIMS(tensor),ngrow,-1,319)
      call checkbound(fablo,fabhi,DIMS(vel),ngrow,-1,320)

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

       call checkbound(fablo,fabhi,DIMS(vel),ngrow+1,-1,321)

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

        if ((viscoelastic_model.eq.0).or. &
            (Viscoelastic_model.eq.1)) then
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
           print *,"num_materials_viscoelastic invalid"
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

          if ((viscoelastic_model.eq.0).or. &
              (viscoelastic_model.eq.1)) then

           do ii=1,3
            traceA=traceA+Q(ii,ii)+one
            if (Q(ii,ii)+one.le.zero) then
             print *,"WARNING: A violates pos. def"
             print *,"ii,Q(ii,ii)= ",ii,Q(ii,ii)
             print *,"level,finest_level ",level,finest_level
             print *,"im_parm,ngrow,i,j,k = ",im_parm,ngrow,i,j,k
             print *,"dt= ",dt
            endif
           enddo

            ! elastic_time*(1-Tr(A)/L^2)
           call get_mod_elastic_time(elastic_time,traceA, &
            polymer_factor,modtime)

           ! modtime=elastic_time >> 1
          else if ((viscoelastic_model.eq.2).or. & !displacement gradient
                   (viscoelastic_model.eq.3)) then !incremental
           modtime=elastic_time
          else
           print *,"viscoelastic_model invalid"
           stop
          endif

          if (modtime.lt.zero) then
           print *,"modtime invalid"
           stop
          endif
         
           ! etaS=etaL-etaP=viscconst-elastic_viscosity 
          if (modtime+dt.le.zero) then
           viscoelastic_coeff=zero
          else
           if (viscoelastic_model.eq.0) then
            viscoelastic_coeff= &
             (visc(D_DECL(i,j,k),im_parm)-etaS)/(modtime+dt)
           else if (viscoelastic_model.eq.1) then
            viscoelastic_coeff= &
             (visc(D_DECL(i,j,k),im_parm)-etaS)
           else if (viscoelastic_model.eq.2) then !displacement gradient
            viscoelastic_coeff=elastic_viscosity
           else if (viscoelastic_model.eq.3) then !incremental
            viscoelastic_coeff=elastic_viscosity
           else
            print *,"viscoelastic_model invalid"
            stop
           endif
          endif
          if (viscoelastic_coeff.lt.zero) then
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
      end subroutine FORT_DERVISCOSITY



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

      subroutine FORT_DERMAGTRACE( &
       level, &
       finest_level, &
       im, &   ! im=0..nmat-1
       ntensor, &
       cellten,DIMS(cellten), &
       vof,DIMS(vof), &
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
       elastic_viscosity)

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T level,finest_level
      INTEGER_T ntensor
      INTEGER_T ncomp_den
      INTEGER_T ncomp_tensor
      INTEGER_T ncomp_visc
      INTEGER_T n_trace
      INTEGER_T nmat

      REAL_T polymer_factor(nmat)
      REAL_T etaS(nmat)
      REAL_T etaP(nmat)
      REAL_T Carreau_beta(nmat)
      REAL_T elastic_time(nmat)
      INTEGER_T viscoelastic_model(nmat)
      REAL_T elastic_viscosity(nmat)

      INTEGER_T tilelo(SDIM), tilehi(SDIM)
      INTEGER_T fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact
      INTEGER_T DIMDEC(vof)
      INTEGER_T DIMDEC(dest)
      INTEGER_T DIMDEC(den)
      INTEGER_T DIMDEC(tensor)
      INTEGER_T DIMDEC(vel)
      INTEGER_T DIMDEC(visc)
      INTEGER_T DIMDEC(cellten)
      INTEGER_T bc(SDIM,2,SDIM)
      INTEGER_T ngrow
      REAL_T time
      REAL_T dx(SDIM), xlo(SDIM)
      REAL_T cellten(DIMV(cellten),ntensor)
      REAL_T vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T dest(DIMV(dest),5)
      REAL_T den(DIMV(den),ncomp_den)
      REAL_T tensor(DIMV(tensor),ncomp_tensor)
      REAL_T vel(DIMV(vel),SDIM)
      REAL_T visc(DIMV(visc),ncomp_visc)

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
      if ((viscoelastic_model(im+1).ge.0).and. &
          (viscoelastic_model(im+1).le.3)) then
       ! do nothing
      else
       print *,"viscoelastic_model(im+1) invalid"
       stop
      endif
      if (fort_elastic_viscosity(im+1).ne.elastic_viscosity(im+1)) then
       print *,"elastic_viscosity(im+1) invalid"
       stop
      endif
      if (elastic_viscosity(im+1).gt.zero) then
       if (num_materials_viscoelastic.le.0) then
        print *,"num_materials_viscoelastic.le.0"
        stop
       endif
       if (ncomp_tensor.ne.FORT_NUM_TENSOR_TYPE) then
        print *,"ncomp_tensor.ne.FORT_NUM_TENSOR_TYPE"
        stop
       endif
      else if (elastic_viscosity(im+1).eq.zero) then
       if (ncomp_tensor.ne.ncomp_den) then
        print *,"ncomp_tensor.ne.ncomp_den"
        stop
       endif
      else
       print *,"elastic_viscosity invalid"
       stop
      endif  

      ! compute u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
      call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)
 
      call checkbound(fablo,fabhi,DIMS(cellten),ngrow,-1,64)
      call checkbound(fablo,fabhi,DIMS(vof),ngrow+1,-1,322)
      call checkbound(fablo,fabhi,DIMS(dest),ngrow,-1,323)
      call checkbound(fablo,fabhi,DIMS(den),ngrow,-1,324)
      call checkbound(fablo,fabhi,DIMS(tensor),ngrow,-1,325)
      call checkbound(fablo,fabhi,DIMS(vel),ngrow+1,-1,326)
      call checkbound(fablo,fabhi,DIMS(visc),ngrow,-1,327)

      iproject=0
      onlyscalar=1  ! mag(trace gradu)
        
       ! in: DERMAGTRACE
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
        print *,"levelrz invalid getshear 2"
        stop
       endif

       vort(1)=gradu(3,2)-gradu(2,3)
       vort(2)=gradu(1,3)-gradu(3,1)
       vort(3)=gradu(2,1)-gradu(1,2)
       dest(D_DECL(i,j,k),5)= &
         sqrt(vort(1)**2+vort(2)**2+vort(3)**2)

       if (elastic_viscosity(im+1).eq.zero) then
        dest(D_DECL(i,j,k),2)=dest(D_DECL(i,j,k),1)
        dest(D_DECL(i,j,k),3)=dest(D_DECL(i,j,k),1)
        dest(D_DECL(i,j,k),4)=dest(D_DECL(i,j,k),1)
       else if (elastic_viscosity(im+1).gt.zero) then
        if (ncomp_visc.ne.3*nmat) then
         print *,"ncomp_visc invalid"
         stop
        endif

        T11=tensor(D_DECL(i,j,k),1)
        T22=tensor(D_DECL(i,j,k),3)
        T33=tensor(D_DECL(i,j,k),4)
        traceA=T11+T22+T33+three
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
        print *,"elastic_viscosity invalid"
        stop
       endif

      enddo
      enddo
      enddo

      return
      end subroutine FORT_DERMAGTRACE

       ! gravity_normalized>0 if pointing downwards
       ! 1<=gravity_dir<=dim
       ! 1-3 Force, 4-6 Torque, 7-9 moment of inertia, 10-12 center of mass,
       ! 13 solid mass
       ! 1-dim force (drag)
       ! dim+1 - 2 dim pforce (drag)
       ! 2 dim+1 - 3 dim torque
       ! 3 dim+1 - 4 dim ptorque
       ! 4 dim+1     perimeter
      subroutine FORT_GETDRAG( &
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
       ntensorMM, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map)

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
      INTEGER_T, intent(in) :: ntensorMM

      INTEGER_T, intent(in) :: facevisc_index
      INTEGER_T, intent(in) :: faceheat_index
      INTEGER_T, intent(in) :: ncphys

      INTEGER_T, intent(in) :: isweep
      REAL_T, intent(in) :: globalsum(13)
      REAL_T, intent(inout) :: localsum(13)
      REAL_T, intent(in) :: gravity_normalized
      INTEGER_T, intent(in) :: gravdir
      INTEGER_T, intent(in) :: ntenvisco
      REAL_T, intent(in) :: visc_coef
      INTEGER_T, intent(in) :: nmat,rzflag
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

      INTEGER_T, intent(in) :: DIMDEC(solxfab)
      INTEGER_T, intent(in) :: DIMDEC(solyfab)
      INTEGER_T, intent(in) :: DIMDEC(solzfab)

      INTEGER_T, intent(in) :: DIMDEC(pres)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(drag)

      REAL_T, intent(in) :: tdata(DIMV(tdata),ntensorMM)
      REAL_T, intent(in) :: viscoten(DIMV(viscoten),ntenvisco)
      REAL_T, intent(in) :: den(DIMV(den),nmat*num_state_material)
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: slrecon(DIMV(slrecon),nmat*ngeom_recon)
      REAL_T, intent(in) :: levelpc(DIMV(levelpc),nmat*(SDIM+1))
      REAL_T, intent(in) :: vol(DIMV(vol))
      REAL_T, intent(in) :: areax(DIMV(areax))
      REAL_T, intent(in) :: areay(DIMV(areay))
      REAL_T, intent(in) :: areaz(DIMV(areaz))

      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: yface(DIMV(yface),ncphys)
      REAL_T, intent(in) :: zface(DIMV(zface),ncphys)
      REAL_T, intent(in) :: cvisc(DIMV(cvisc))

      REAL_T, intent(in) :: solxfab(DIMV(solxfab),nparts_def*SDIM)
      REAL_T, intent(in) :: solyfab(DIMV(solyfab),nparts_def*SDIM)
      REAL_T, intent(in) :: solzfab(DIMV(solzfab),nparts_def*SDIM)

      REAL_T, intent(in) :: pres(DIMV(pres),num_materials_vel)
      REAL_T, intent(in) :: vel(DIMV(vel),SDIM*num_materials_vel)
      REAL_T, intent(inout) :: drag(DIMV(drag),4*SDIM+1)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T ii,jj,kk
      INTEGER_T i,j,k
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
      INTEGER_T i1,j1,n
      INTEGER_T vofcomp,dencomp,viscbase
      INTEGER_T icell,jcell,kcell
      REAL_T vel6point(SDIM,2,SDIM)
      REAL_T ls_visc(nmat)
      REAL_T lsleft(nmat)
      REAL_T lsright(nmat)
      REAL_T solid_dist_primary
      INTEGER_T im_solid_crit
      REAL_T stress(3),pressure_load(3)

      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_face(-3:3,SDIM)
      INTEGER_T nhalf

      REAL_T volume
      REAL_T gradu(SDIM,SDIM)

      REAL_T facearea
      INTEGER_T im
      INTEGER_T im_visc
      INTEGER_T im_fluid
      INTEGER_T im_test
      INTEGER_T im_primary
      INTEGER_T im_left
      INTEGER_T im_right
      INTEGER_T FSI_exclude
      REAL_T volgrid,mass
      REAL_T cengrid(SDIM)
      REAL_T global_centroid(SDIM)
      REAL_T rvec(3),gravvector(3),rcross(3),prcross(3)
      REAL_T localtorque(SDIM)
      REAL_T plocaltorque(SDIM)
      REAL_T Q(3,3)
      REAL_T viscous_stress,elastic_stress
      REAL_T nsolid(SDIM)
      INTEGER_T mask_cell
      REAL_T ls_sort(nmat)
      INTEGER_T sorted_list(nmat)
      REAL_T mu
      REAL_T delx
      INTEGER_T partid
      INTEGER_T ibase

      nhalf=3

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in getdrag"
      else if (time.lt.zero) then
       print *,"time invalid in getdrag"
       stop
      else
       print *,"time bust in getdrag"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_GETDRAG"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_GETDRAG"
       stop
      endif

        ! cell centered grad U
      call checkbound(fablo,fabhi,DIMS(tdata),0,-1,1252)
      call checkbound(fablo,fabhi,DIMS(viscoten),1,-1,1253)
      call checkbound(fablo,fabhi,DIMS(den),1,-1,1254)
      call checkbound(fablo,fabhi,DIMS(mask),1,-1,1255)
      call checkbound(fablo,fabhi,DIMS(slrecon),1,-1,12560)
      call checkbound(fablo,fabhi,DIMS(levelpc),1,-1,1257)
      call checkbound(fablo,fabhi,DIMS(vol),0,-1,6600)
      call checkbound(fablo,fabhi,DIMS(areax),0,0,6601)
      call checkbound(fablo,fabhi,DIMS(areay),0,1,6602)
      call checkbound(fablo,fabhi,DIMS(areaz),0,SDIM-1,6603)

      call checkbound(fablo,fabhi,DIMS(xface),0,0,1258)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,1259)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,1261)
      call checkbound(fablo,fabhi,DIMS(cvisc),0,-1,1262)

      call checkbound(fablo,fabhi,DIMS(solxfab),0,0,6604)
      call checkbound(fablo,fabhi,DIMS(solyfab),0,1,6604)
      call checkbound(fablo,fabhi,DIMS(solzfab),0,SDIM-1,6604)

      call checkbound(fablo,fabhi,DIMS(pres),1,-1,6605)
      call checkbound(fablo,fabhi,DIMS(vel),1,-1,6606)
      call checkbound(fablo,fabhi,DIMS(drag),0,-1,6607)

      if (bfact.lt.1) then
       print *,"bfact invalid5"
       stop
      endif
      if (visc_coef.lt.zero) then
       print *,"visc_coef invalid"
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
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
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
      if (ntensorMM.ne.ntensor*num_materials_vel) then
       print *,"ntensorMM invalid"
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
      if ((num_materials_viscoelastic.ge.1).and. &
          (num_materials_viscoelastic.le.nmat)) then
       if (FORT_NUM_TENSOR_TYPE*num_materials_viscoelastic.ne.ntenvisco) then
        print *,"ntenvisco invalid1"
        stop
       endif
      else
       print *,"num_materials_viscoelastic invalid"
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

       ! first sweep - find the mass and centroid of the solid 
      if (isweep.eq.0) then
 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        do n=1,4*SDIM+1
         drag(D_DECL(i,j,k),n)=zero
        enddo
        do im_test=1,nmat
         if (is_rigid(nmat,im_test).eq.0) then
          ! do nothing
         else if (is_rigid(nmat,im_test).eq.1) then
          mask_cell=NINT(mask(D_DECL(i,j,k)))
          if (mask_cell.eq.1) then
           call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
           call Box_volumeFAST(bfact,dx,xsten,nhalf, &
            volgrid,cengrid,SDIM)
           vofcomp=(im_test-1)*ngeom_recon+1
           dencomp=(im_test-1)*num_state_material+1
           mass=den(D_DECL(i,j,k),dencomp)*volgrid* &
            slrecon(D_DECL(i,j,k),vofcomp)
           localsum(13)=localsum(13)+mass
           ! 1-3 force, 4-6 torque, 7-9 moment of inertia
           do dir=1,SDIM
            localsum(9+dir)=localsum(9+dir)+ &
             mass* &
             (slrecon(D_DECL(i,j,k),vofcomp+dir)+cengrid(dir))
           enddo
          else if (mask_cell.eq.0) then
           ! do nothing
          else
           print *,"mask_cell invalid"
           stop
          endif 
         else
          print *,"is_rigid invalid 90"
          stop
         endif

        enddo ! im_test=1..nmat

       enddo
       enddo
       enddo  ! isweep=0 centroids

      else if (isweep.eq.1) then ! above, mass and centroid

       do im_test=1,nmat
        if (is_rigid(nmat,im_test).eq.0) then
         ! do nothing
        else if (is_rigid(nmat,im_test).eq.1) then

         mass=globalsum(13)
         if (mass.lt.zero) then
          print *,"mass cannot be negative  im_test,mass=",im_test,mass
          stop
         else if (mass.eq.zero) then
          print *,"WARNING: FSI_flag(im_test)=1,2,4 mass=0"
          print *,"im_test=",im_test
          print *,"FSI_flag(im_test)=",FSI_flag(im_test)
          do dir=1,SDIM
           global_centroid(dir)=zero
          enddo
         else if (mass.gt.zero) then
          do dir=1,SDIM
           global_centroid(dir)=globalsum(9+dir)/mass
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
           print *,"levelrz invalid getdrag"
           stop
          endif
         else
          print *,"mass bust"
          stop
         endif
        else
         print *,"im_rigid invalid 91"
         stop
        endif
       enddo ! im_test=1..nmat
     
        ! buoyancy force (body forces in the solid) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        do im_test=1,nmat
         if (is_rigid(nmat,im_test).eq.0) then
          ! do nothing
         else if (is_rigid(nmat,im_test).eq.1) then
          mask_cell=NINT(mask(D_DECL(i,j,k)))
          if (mask_cell.eq.1) then
           call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
           call Box_volumeFAST(bfact,dx,xsten,nhalf, &
            volgrid,cengrid,SDIM)
           vofcomp=(im_test-1)*ngeom_recon+1
           dencomp=(im_test-1)*num_state_material+1
           mass=den(D_DECL(i,j,k),dencomp)*volgrid* &
            slrecon(D_DECL(i,j,k),vofcomp)
           do dir=1,SDIM
            gravvector(dir)=zero
            rvec(dir)=slrecon(D_DECL(i,j,k),vofcomp+dir)+ &
             cengrid(dir)-global_centroid(dir)
           enddo
           gravvector(gravdir)=-mass*gravity_normalized
           if (SDIM.eq.2) then
            rvec(3)=zero
            gravvector(3)=zero
           endif 

           do dir=1,SDIM
            localsum(dir)=localsum(dir)+gravvector(dir)
            drag(D_DECL(i,j,k),dir)= &
             drag(D_DECL(i,j,k),dir)+gravvector(dir)
           enddo

           ! torque=r cross F
           ! a x b=(a2 b3 - a3 b2) hat(i)+
           !       (-a1 b3 + b1 a3) hat(j)+
           !       (a1 b2 - a2 b1) hat(k)
           call crossprod(rvec,gravvector,rcross)
           localtorque(1)=rcross(3) ! x-y plane of rotation
           if (SDIM.eq.3) then
            localtorque(2)=rcross(1) ! y-z plane of rotation
            localtorque(SDIM)=rcross(2) ! x-z plane of rotation
           endif

           localsum(7)=localsum(7)+(rvec(1)**2+rvec(2)**2)*mass ! x-y
           dirend=1
           if (SDIM.eq.3) then
            localsum(8)=localsum(8)+(rvec(2)**2+rvec(3)**2)*mass ! y-z
            localsum(9)=localsum(9)+(rvec(1)**2+rvec(3)**2)*mass ! x-z
            dirend=3
           endif
 
           do dir=1,dirend
            localsum(3+dir)=localsum(3+dir)+localtorque(dir)
            drag(D_DECL(i,j,k),2*SDIM+dir)= &
             drag(D_DECL(i,j,k),2*SDIM+dir)+localtorque(dir)
           enddo

          else if (mask_cell.eq.0) then
           ! do nothing
          else
           print *,"mask_cell invalid"
           stop
          endif 
         else
          print *,"is_rigid invalid 92"
          stop
         endif

        enddo ! im_test=1..nmat

       enddo
       enddo
       enddo  ! i,j,k (cells - body forces and torques (buoancy) )

        ! cells - pressure and viscosity
       do icell=growlo(1),growhi(1)
       do jcell=growlo(2),growhi(2)
       do kcell=growlo(3),growhi(3)

        mask_cell=NINT(mask(D_DECL(icell,jcell,kcell)))

        if (mask_cell.eq.1) then

         do im=1,nmat
          ls_sort(im)=levelpc(D_DECL(icell,jcell,kcell),im)
         enddo
         call get_primary_material(ls_sort,nmat,im_primary)

          ! is center cell in the fluid?
         if (is_rigid(nmat,im_primary).eq.0) then

          call gridsten(xsten,xlo, &
           icell,jcell,kcell,fablo,bfact,dx,nhalf)

          volume=vol(D_DECL(icell,jcell,kcell))
          presmag=pres(D_DECL(icell,jcell,kcell),1)

           ! returns: im_solid_crit=max_{is_rigid(im)==1} ls_sort(im)
          call combine_solid_LS(ls_sort,nmat,solid_dist_primary,im_solid_crit)
          if ((im_solid_crit.ge.0).and.(im_solid_crit.le.nmat)) then
           ! do nothing
          else
           print *,"im_solid_crit invalid in getdrag:",im_solid_crit
           stop
          endif

          FSI_exclude=1
          call sort_volume_fraction(ls_sort,FSI_exclude,sorted_list,nmat)
          im_fluid=sorted_list(1)
          if ((im_fluid.ge.1).and.(im_fluid.le.nmat)) then
           ! do nothing
          else
           print *,"im_fluid invalid in GETDRAG"
           stop
          endif

           ! first check for change of sign of the solid LS.
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

            i=icell+side_cell*ii
            j=jcell+side_cell*jj
            k=kcell+side_cell*kk

            do im=1,nmat
              ! cell to left of face
             lsleft(im)=levelpc(D_DECL(i-ii,j-jj,k-kk),im)
              ! cell to right of face
             lsright(im)=levelpc(D_DECL(i,j,k),im)
            enddo
            call get_primary_material(lsleft,nmat,im_left)
            call get_primary_material(lsright,nmat,im_right)
            if ((is_rigid(nmat,im_left).eq.1).and. &
                (is_rigid(nmat,im_right).eq.1)) then
             ! do nothing - both adjoining cells in the solid
            else if ((is_rigid(nmat,im_left).eq.0).and. &
                     (is_rigid(nmat,im_right).eq.0)) then
             ! do nothing - both adjoining cells in the fluid
            else if ( ((is_rigid(nmat,im_left).eq.1).and. &
                       (is_rigid(nmat,im_right).eq.0)).or. &
                      ((is_rigid(nmat,im_left).eq.0).and. &
                       (is_rigid(nmat,im_right).eq.1)) ) then

             if (mask_cell.eq.1) then
              
              if (facedir.eq.1) then
               facearea=areax(D_DECL(i,j,k))
              else if (facedir.eq.2) then
               facearea=areay(D_DECL(i,j,k))
              else if ((facedir.eq.3).and.(SDIM.eq.3)) then
               facearea=areaz(D_DECL(i,j,k))
              else
               print *,"facedir invalid"
               stop
              endif

               ! nsolid points into the solid
              do dir=1,SDIM 
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

               !facedir=1..sdim
              call gridstenMAC(xsten_face,xlo,i,j,k,fablo,bfact, &
                dx,nhalf,facedir-1)

              ! im_primary is a fluid at cell (icell,jcell,kcell)
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

              do j1=1,SDIM
              do i1=1,SDIM
               Q(i1,j1)=zero
              enddo
              enddo

              if (fort_elastic_viscosity(im_fluid).gt.zero) then
               partid=1
               do while ((fort_im_elastic_map(partid)+1.ne.im_fluid).and. &
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
                print *,"partid invalid in GETDRAG"
                stop
               endif
              else if (fort_elastic_viscosity(im_fluid).eq.zero) then
               ! do nothing
              else
               print *,"fort_elastic_viscosity(im_fluid) invalid"
               stop
              endif

              mu=get_user_viscconst(im_fluid, &
                fort_denconst(im_fluid),fort_tempconst(im_fluid))

              do j1=1,SDIM
               stress(j1)=zero
               do i1=1,SDIM
                viscous_stress= &
                 mu*visc_coef* &
                  (gradu(i1,j1)+gradu(j1,i1))
                elastic_stress=Q(i1,j1)
                stress(j1)=stress(j1)- &
                  (viscous_stress+elastic_stress)*nsolid(i1)
               enddo
               pressure_load(j1)=presmag*nsolid(j1)
               stress(j1)=stress(j1)*facearea
               pressure_load(j1)=pressure_load(j1)*facearea
              enddo  ! j1

! global_centroid will be incorrect if the solid geometry is reflected
! across a domain boundary.

              do dir=1,SDIM
               rvec(dir)=xsten_face(0,dir)-global_centroid(dir)
              enddo
              if (SDIM.eq.2) then
               rvec(3)=zero
               stress(3)=zero
               pressure_load(3)=zero
              endif
              call crossprod(rvec,stress,rcross)
              call crossprod(rvec,pressure_load,prcross)
              localtorque(1)=rcross(3) ! x-y plane of rotation
              plocaltorque(1)=prcross(3) ! x-y plane of rotation
              if (SDIM.eq.3) then
               localtorque(2)=rcross(1) ! y-z plane of rotation
               plocaltorque(2)=prcross(1) ! y-z plane of rotation
               localtorque(SDIM)=rcross(2) ! x-z plane of rotation
               plocaltorque(SDIM)=prcross(2) ! x-z plane of rotation
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

              do n=1,SDIM
               drag(D_DECL(icell,jcell,kcell),n)= &
                drag(D_DECL(icell,jcell,kcell),n)+stress(n)
               drag(D_DECL(icell,jcell,kcell),n)= &
                drag(D_DECL(icell,jcell,kcell),n)+pressure_load(n)

               localsum(n)=localsum(n)+stress(n)
               localsum(n)=localsum(n)+pressure_load(n)

               drag(D_DECL(icell,jcell,kcell),n+SDIM)= &
                drag(D_DECL(icell,jcell,kcell),n+SDIM)+pressure_load(n)
              enddo

              dirend=1
              if (SDIM.eq.3) then
               dirend=3
              endif

              do dir=1,dirend
               localsum(3+dir)=localsum(3+dir)+localtorque(dir)
               localsum(3+dir)=localsum(3+dir)+plocaltorque(dir)

               drag(D_DECL(icell,jcell,kcell),2*SDIM+dir)= &
                drag(D_DECL(icell,jcell,kcell),2*SDIM+dir)+ &
                localtorque(dir)
               drag(D_DECL(icell,jcell,kcell),2*SDIM+dir)= &
                drag(D_DECL(icell,jcell,kcell),2*SDIM+dir)+ &
                plocaltorque(dir)
               drag(D_DECL(icell,jcell,kcell),3*SDIM+dir)= &
                drag(D_DECL(icell,jcell,kcell),3*SDIM+dir)+ &
                plocaltorque(dir)
              enddo ! dir=1..dirend

              drag(D_DECL(icell,jcell,kcell),4*SDIM+1)= &
               drag(D_DECL(icell,jcell,kcell),4*SDIM+1)+facearea

             else if (mask_cell.eq.0) then
              ! do nothing
             else
              print *,"mask_cell invalid"
              stop
             endif

            else
             print *,"is_rigid(nmat,im_left) or is_rigid(nmat,im_right) bad"
             stop
            endif

           enddo ! side_cell=0,1
          enddo ! facedir=1..sdim 

         else if (is_rigid(nmat,im_primary).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(nmat,im_primary) invalid"
          stop
         endif

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
      end subroutine FORT_GETDRAG


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
      ncomp_state_test=num_materials_vel*(SDIM+1)+ &
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

           vofcomp=(SDIM+1)*num_materials_vel+ &
             nmat*num_state_material+(im-1)*ngeom_raw+1
           dencomp=(SDIM+1)*num_materials_vel+ &
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

           vofcomp=(SDIM+1)*num_materials_vel+ &
             nmat*num_state_material+(im-1)*ngeom_raw+1
           dencomp=(SDIM+1)*num_materials_vel+ &
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
      ncomp_state_test=num_materials_vel*(SDIM+1)+ &
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
       
       vofcomp=(SDIM+1)*num_materials_vel+ &
          nmat*num_state_material+im_source*ngeom_raw+1
       vfrac=snew(D_DECL(i,j,k),vofcomp)
       if (vfrac.ge.half) then
        do im=1,nmat
         dencomp=(SDIM+1)*num_materials_vel+ &
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
      REAL_T, intent(in) :: vel(DIMV(vel),num_materials_vel*(SDIM+1))
      REAL_T, intent(in) :: velx(DIMV(velx))
      REAL_T, intent(in) :: vely(DIMV(vely))
      REAL_T, intent(in) :: velz(DIMV(velz))
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T dir,im,ibase
      REAL_T magvel,magvel_mac,magvel_collide
      REAL_T vello,velhi,vellohi
      INTEGER_T i,j,k

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if (num_materials_vel.eq.1) then
       ! do nothing
      else
       print *,"num_materials_vel invalid"
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

        ibase=num_materials_vel*SDIM
        do im=1,num_materials_vel
         if (vel(D_DECL(i,j,k),ibase+im).gt.maxpres) then
          maxpres=vel(D_DECL(i,j,k),ibase+im)
         endif
         if (vel(D_DECL(i,j,k),ibase+im).lt.minpres) then
          minpres=vel(D_DECL(i,j,k),ibase+im)
         endif
        enddo
        do im=1,num_materials_vel
         ibase=(im-1)*SDIM
         magvel=zero
         magvel_mac=zero
         magvel_collide=zero
         do dir=1,SDIM
          magvel=magvel+vel(D_DECL(i,j,k),ibase+dir)**2
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
         enddo
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
        enddo ! im
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


