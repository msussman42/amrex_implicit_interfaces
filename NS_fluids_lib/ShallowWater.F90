#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "ShallowWater_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

      module shallowwater_module
      use global_utility_module
      use probf90_module

      contains

      subroutine shallow_water_pressure(density,internal_energy,pressure)
      IMPLICIT NONE

      REAL_T density,internal_energy,pressure
      REAL_T R,cp,cv,gamma_constant,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
      if (density.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=omega*density*internal_energy

      return
      end subroutine shallow_water_pressure


      subroutine shallow_water_sound_speed_sqr(density,internal_energy, &
        sound_speed_sqr)
      IMPLICIT NONE

      REAL_T density,internal_energy,pressure,sound_speed_sqr
      REAL_T R,cp,cv,gamma_constant,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
      if (density.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal energy cannot be negative"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif
      pressure=omega*density*internal_energy
      sound_speed_sqr=gamma_constant*pressure/density

      return
      end subroutine shallow_water_sound_speed_sqr


      subroutine shallow_water_internal(density,temperature,internal_energy)
      IMPLICIT NONE

      REAL_T density,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (density.le.zero) then
       print *,"density negative"
       stop
      endif
      if (temperature.le.zero) then
       print *,"temperature cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      internal_energy=cv*temperature

      return
      end subroutine shallow_water_internal


      subroutine shallow_water_temperature(density,temperature,internal_energy)
      IMPLICIT NONE

      REAL_T density,temperature,internal_energy,gamma_constant
      REAL_T cp,cv,R,omega

      call air_parms(R,cp,cv,gamma_constant,omega)
    
      if (density.le.zero) then
       print *,"density negative"
       stop
      endif
      if (internal_energy.le.zero) then
       print *,"internal_energy cannot be <=0"
       stop
      endif
      if (cv.le.zero) then
       print *,"cv error"
       stop
      endif
      if (cp.le.zero) then
       print *,"cp error"
       stop
      endif

      temperature=internal_energy/cv

      return
      end subroutine shallow_water_temperature


      subroutine convert_to_conserve(primitive,conserve,total_mass, &
        nstate,nflux,xcell_cen,nmat,normdir)
      IMPLICIT NONE

      INTEGER_T normdir
      INTEGER_T nmat
      INTEGER_T nstate
      INTEGER_T nstate_test
      INTEGER_T nflux
      INTEGER_T nflux_test
      REAL_T xcell_cen(SDIM)
      REAL_T primitive(nflux) 
      REAL_T conserve(nflux) 

      INTEGER_T dir
      INTEGER_T im
      INTEGER_T igeom,ithermal,ils,icomp,irest
      INTEGER_T vofcomp,lscomp
      REAL_T total_mass
      REAL_T total_temperature
      REAL_T total_internal_energy
      REAL_T total_pressure
      REAL_T vfrac(nmat)
      REAL_T density(nmat)
      REAL_T temperature(nmat)
      REAL_T internal_energy(nmat)
      REAL_T total_kinetic_energy
      REAL_T total_energy
   
      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif 

      nstate_test=(SDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1
      nflux_test=nstate_test+nmat*(SDIM+1)

      if (nflux.ne.nflux_test) then
       print *,"nflux invalid"
       stop
      endif
      if (nstate.ne.nstate_test) then
       print *,"nstate invalid"
       stop
      endif 
      igeom=(AMREX_SPACEDIM+1)+nmat*num_state_material
      ithermal=(AMREX_SPACEDIM+1)
      ils=(AMREX_SPACEDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1

      total_mass=zero
      total_temperature=zero
      do im=1,nmat
       vfrac(im)=primitive(igeom+(im-1)*ngeom_raw+1)
       density(im)=primitive(ithermal+(im-1)*num_state_material+1)
       if (density(im).le.zero) then
        print *,"density(im).le.zero"
        stop
       endif
       temperature(im)=primitive(ithermal+(im-1)*num_state_material+2)
       call shallow_water_internal(density(im),temperature(im), &
        internal_energy(im))
       total_mass=total_mass+vfrac(im)*density(im)
       total_temperature=total_temperature+vfrac(im)*temperature(im)
      enddo

      if (total_mass.le.zero) then
       print *,"total_mass invalid"
       stop
      endif

      call shallow_water_internal(total_mass,total_temperature, &
        total_internal_energy)
      call shallow_water_pressure(total_mass,total_internal_energy, &
        total_pressure)

      total_kinetic_energy=zero
      do dir=1,SDIM
       total_kinetic_energy=total_kinetic_energy+ &
          half*total_mass*(primitive(dir)**2)
      enddo

      do dir=1,SDIM
       conserve(dir)=primitive(dir)*total_mass
      enddo
      conserve(nstate)=primitive(normdir+1)

      do im=1,nmat
       icomp=ithermal+(im-1)*num_state_material+1
       conserve(icomp)=total_mass ! density
       icomp=icomp+1
       total_energy=total_mass*total_internal_energy+total_kinetic_energy
        ! energy
       conserve(icomp)=total_energy
       do irest=1,num_state_material-2
        icomp=ithermal+(im-1)*num_state_material+2+irest
        conserve(icomp)=primitive(icomp)*total_mass
       enddo
      enddo ! im=1..nmat

      ! for centroids: xcen=int (x-xcellcen) H/int H=int xH/int H - xcellcen
      ! xcen+xcellcen=int xH/int H
      ! (xcen+xcellcen)int H = int xH
      ! (xcen+xcellcen)F = int xH/vol
      do im=1,nmat
       vofcomp=igeom+(im-1)*ngeom_raw+1
       conserve(vofcomp)=primitive(vofcomp)  ! volume fraction
       do dir=1,SDIM
        conserve(vofcomp+dir)= &
         (primitive(vofcomp+dir)+xcell_cen(dir))*primitive(vofcomp)
       enddo
       lscomp=ils+(im-1)*(AMREX_SPACEDIM+1)+1
       conserve(lscomp)=primitive(lscomp)
       do dir=1,SDIM
        conserve(lscomp+dir)=primitive(lscomp+dir)
       enddo
      enddo ! im=1..nmat

      return
      end subroutine convert_to_conserve


      subroutine convert_to_primitive(primitive,conserve,nstate,nflux, &
       xcell_cen,nmat,normdir)
      IMPLICIT NONE

      INTEGER_T normdir
      INTEGER_T nmat
      INTEGER_T nflux
      INTEGER_T nflux_test
      INTEGER_T nstate
      INTEGER_T nstate_test
      REAL_T xcell_cen(SDIM)
      REAL_T primitive(nflux) 
      REAL_T conserve(nflux) 

      INTEGER_T dir
      INTEGER_T im
      INTEGER_T igeom,ithermal,ils,icomp,irest
      INTEGER_T vofcomp,lscomp
      REAL_T total_mass
      REAL_T total_temperature
      REAL_T total_internal_energy
      REAL_T total_pressure
      REAL_T vfrac(nmat)
      REAL_T density(nmat)
      REAL_T temperature(nmat)
      REAL_T internal_energy(nmat)
      REAL_T energy(nmat)
      REAL_T total_kinetic_energy
    

      nstate_test=(SDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1
      nflux_test=nstate_test+nmat*(SDIM+1)

      if (nflux.ne.nflux_test) then
       print *,"nflux invalid"
       stop
      endif
      if (nstate.ne.nstate_test) then
       print *,"nstate invalid"
       stop
      endif
      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif
 
      igeom=(AMREX_SPACEDIM+1)+nmat*num_state_material
      ithermal=(AMREX_SPACEDIM+1)
      ils=(AMREX_SPACEDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1

      total_mass=zero
      total_temperature=zero
      do im=1,nmat
       vfrac(im)=conserve(igeom+(im-1)*ngeom_raw+1)
       density(im)=conserve(ithermal+(im-1)*num_state_material+1)
       if (density(im).le.zero) then
        print *,"density(im).le.zero"
        stop
       endif
       energy(im)=conserve(ithermal+(im-1)*num_state_material+2)
       total_mass=total_mass+vfrac(im)*density(im)
      enddo ! im=1..nmat

      if (total_mass.le.zero) then
       print *,"total_mass invalid"
       stop
      endif

       ! velocity
      do dir=1,SDIM
       primitive(dir)=conserve(dir)/total_mass
      enddo
      primitive(nstate)=primitive(normdir+1)

      total_kinetic_energy=zero
      do dir=1,SDIM
       total_kinetic_energy=total_kinetic_energy+ &
          half*total_mass*(primitive(dir)**2)
      enddo

      do im=1,nmat
       internal_energy(im)=(energy(im)-total_kinetic_energy)/total_mass
       call shallow_water_temperature(density(im),temperature(im), &
        internal_energy(im))
       total_temperature=total_temperature+vfrac(im)*temperature(im)
      enddo ! im=1..nmat

      call shallow_water_internal(total_mass,total_temperature, &
        total_internal_energy)
      call shallow_water_pressure(total_mass,total_internal_energy, &
        total_pressure)

      do im=1,nmat
       icomp=ithermal+(im-1)*num_state_material+1
       primitive(icomp)=total_mass ! density
       icomp=icomp+1
        ! energy
       primitive(icomp)=total_temperature
       do irest=1,num_state_material-2
        icomp=ithermal+(im-1)*num_state_material+2+irest
        primitive(icomp)=conserve(icomp)/total_mass
       enddo
      enddo ! im=1..nmat

      ! for centroids: xcen=int (x-xcellcen) H/int H=int xH/int H - xcellcen
      ! xcen+xcellcen=int xH/int H
      ! (xcen+xcellcen)int H = int xH
      ! (xcen+xcellcen)F = int xH/vol
      do im=1,nmat
       vofcomp=igeom+(im-1)*ngeom_raw+1
       primitive(vofcomp)=conserve(vofcomp)  ! volume fraction
       do dir=1,SDIM
        if (primitive(vofcomp).eq.zero) then
         primitive(vofcomp+dir)=zero
        else if (primitive(vofcomp).gt.zero) then
         primitive(vofcomp+dir)= &
          conserve(vofcomp+dir)/primitive(vofcomp)-xcell_cen(dir)
        else
         print *,"primitive(vofcomp) invalid"
         stop
        endif
       enddo ! dir=1..sdim
       lscomp=ils+(im-1)*(AMREX_SPACEDIM+1)+1
       primitive(lscomp)=conserve(lscomp)
       do dir=1,SDIM
        primitive(lscomp+dir)=conserve(lscomp+dir)
       enddo
      enddo ! im=1..nmat

      return
      end subroutine convert_to_primitive


      subroutine shallow_water_flux(state,flux,soundsqr, &
        max_speed,nstate,nflux,normal_vel,nmat,normdir,xcell_cen)
      IMPLICIT NONE

      INTEGER_T nmat
      INTEGER_T normdir
      INTEGER_T nflux
      INTEGER_T nflux_test
      INTEGER_T nstate
      INTEGER_T nstate_test
      REAL_T normal_vel
      REAL_T state(nflux)
      REAL_T conserve(nflux)
      REAL_T flux(nflux)
      REAL_T xcell_cen(SDIM)
      INTEGER_T dir
      INTEGER_T im
      INTEGER_T igeom,ithermal,ils,icomp,irest
      INTEGER_T vofcomp,lscomp
      REAL_T total_mass
      REAL_T total_temperature
      REAL_T total_internal_energy
      REAL_T total_pressure
      REAL_T soundsqr
      REAL_T max_speed
      REAL_T vfrac(nmat)
      REAL_T density(nmat)
      REAL_T temperature(nmat)
      REAL_T internal_energy(nmat)
      REAL_T total_kinetic_energy
      REAL_T total_energy
    
      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif

      nstate_test=(SDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1
      nflux_test=nstate_test+nmat*(SDIM+1)

      if (nflux.ne.nflux_test) then
       print *,"nflux invalid"
       stop
      endif
      if (nstate.ne.nstate_test) then
       print *,"nstate invalid"
       stop
      endif
 
      do dir=1,nflux
       flux(dir)=zero
      enddo

      igeom=(AMREX_SPACEDIM+1)+nmat*num_state_material
      ithermal=(AMREX_SPACEDIM+1)
      ils=(AMREX_SPACEDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1

      call convert_to_conserve(state,conserve,total_mass, &
         nstate,nflux,xcell_cen,nmat,normdir)

      total_temperature=zero
      do im=1,nmat
       vfrac(im)=state(igeom+(im-1)*ngeom_raw+1)
       density(im)=state(ithermal+(im-1)*num_state_material+1)
       temperature(im)=state(ithermal+(im-1)*num_state_material+2)
       call shallow_water_internal(density(im),temperature(im), &
        internal_energy(im))
       total_temperature=total_temperature+vfrac(im)*temperature(im)
      enddo ! im=1..nmat

      call shallow_water_internal(total_mass,total_temperature, &
        total_internal_energy)
      call shallow_water_pressure(total_mass,total_internal_energy, &
        total_pressure)
      call shallow_water_sound_speed_sqr(total_mass,total_internal_energy, &
        soundsqr)

      total_kinetic_energy=zero
      max_speed=zero
      do dir=1,SDIM
       total_kinetic_energy=total_kinetic_energy+ &
          half*total_mass*(state(dir)**2)
       max_speed=max(max_speed,abs(state(dir))+sqrt(soundsqr))
      enddo

      do dir=1,SDIM
       flux(dir)=conserve(dir)*normal_vel
      enddo
      flux(normdir+1)=flux(normdir+1)+total_pressure
      flux(nstate)=normal_vel

      do im=1,nmat
       icomp=ithermal+(im-1)*num_state_material+1
       flux(icomp)=total_mass*normal_vel ! density
       icomp=icomp+1
       total_energy=total_mass*total_internal_energy+total_kinetic_energy
        ! energy
       flux(icomp)=total_energy*normal_vel+total_pressure*normal_vel
       do irest=1,num_state_material-2
        icomp=ithermal+(im-1)*num_state_material+2+irest
        flux(icomp)=conserve(icomp)*normal_vel
       enddo
      enddo ! im=1..nmat

      do im=1,nmat
       vofcomp=igeom+(im-1)*ngeom_raw+1
       flux(vofcomp)=conserve(vofcomp)*normal_vel  ! volume fraction
       do dir=1,SDIM
        flux(vofcomp+dir)=conserve(vofcomp+dir)*normal_vel ! centroid
       enddo
       lscomp=ils+(im-1)*(AMREX_SPACEDIM+1)+1
       flux(lscomp)=conserve(lscomp)*normal_vel
       do dir=1,SDIM
        flux(lscomp+dir)=conserve(lscomp+dir)*normal_vel ! slope
       enddo
      enddo ! im=1..nmat
       
      return
      end subroutine shallow_water_flux
 
      subroutine shallow_water_LLF(uLeft,uRight,normal_vel, &
        dir,nstate,nflux,LLF_flux,nmat,xcell_cen_left,xcell_cen_right)
      IMPLICIT NONE

      INTEGER_T dir,nmat
      INTEGER_T nflux
      INTEGER_T nflux_test
      INTEGER_T nstate
      INTEGER_T nstate_test
      REAL_T xcell_cen_left(SDIM)
      REAL_T xcell_cen_right(SDIM)
      REAL_T uLeft(nflux)
      REAL_T uRight(nflux)
      REAL_T LLF_flux(nflux)
      REAL_T flux_left(nflux)
      REAL_T flux_right(nflux)
      REAL_T normal_vel
      REAL_T normal_velLeft
      REAL_T normal_velRight
      REAL_T soundsqr_left
      REAL_T soundsqr_right
      REAL_T max_speed_left
      REAL_T max_speed_right
      REAL_T max_speed
      INTEGER_T n

    
      nstate_test=(AMREX_SPACEDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1
      nflux_test=nstate_test+nmat*(SDIM+1)
      if (nflux.ne.nflux_test) then
       print *,"nflux invalid"
       stop
      endif
      if (nstate.ne.nstate_test) then
       print *,"nstate invalid"
       stop
      endif
      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid"
       stop
      endif
     
      normal_velLeft=uLeft(dir+1) 
      normal_velRight=uRight(dir+1) 
      call shallow_water_flux(uLeft,flux_left,soundsqr_left, &
        max_speed_left,nstate,nflux,normal_vel,nmat, &
        dir,xcell_cen_left)
      call shallow_water_flux(uRight,flux_right,soundsqr_right, &
        max_speed_right,nstate,nflux,normal_vel,nmat, &
        dir,xcell_cen_right)
      max_speed=max(max_speed_left,max_speed_right)

       ! dt max = dx
       ! dx/dt = max
      do n=1,nflux
       LLF_flux(n)=half*(flux_left(n)+flux_right(n))- &
        half*max_speed*(uRight(n)-uLeft(n))
      enddo 
      LLF_flux(nstate)=normal_vel

      return
      end subroutine shallow_water_LLF

      end module shallowwater_module

 
      subroutine FORT_FEED_FORWARD_ADVECT( &
       init_fluxes, &
       normdir, & ! 0..sdim-1
       nflux_feed_forward, & 
       nstate, &
       velbc_in, &
       slab_step, &
       dt, &
       time, &
       xlo, &
       dx, &
       mask,DIMS(mask), & ! 1=fine/fine  0=coarse/fine
       maskcoef,DIMS(maskcoef), & ! 1=not cov. or outside domain  0=covered
       ax,DIMS(ax), &   
       vol,DIMS(vol), &  
       xflux,DIMS(xflux), &  
       velfab,DIMS(velfab), &
       lsfab,DIMS(lsfab), &
       snew,DIMS(snew), &
       sold,DIMS(sold), &
       lsnew,DIMS(lsnew), &
       lsold,DIMS(lsold), &
       xmacnew,DIMS(xmacnew), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       domlo,domhi, &
       nmat)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probf90_module
      use probcommon_module
      use shallowwater_module
      IMPLICIT NONE

      INTEGER_T nmat
      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T init_fluxes
      INTEGER_T normdir
      INTEGER_T nflux_feed_forward
      INTEGER_T nflux_test
      INTEGER_T nstate
      INTEGER_T nstate_test
      INTEGER_T velbc_in(nstate)
      INTEGER_T slab_step
      REAL_T dt
      REAL_T time
      REAL_T xlo(SDIM)
      REAL_T dx(SDIM)
      INTEGER_T DIMDEC(mask)
      INTEGER_T DIMDEC(maskcoef)
      INTEGER_T DIMDEC(ax)
      INTEGER_T DIMDEC(vol)
      INTEGER_T DIMDEC(xflux)
      INTEGER_T DIMDEC(velfab)
      INTEGER_T DIMDEC(lsfab)
      INTEGER_T DIMDEC(snew)
      INTEGER_T DIMDEC(sold)
      INTEGER_T DIMDEC(lsnew)
      INTEGER_T DIMDEC(lsold)
      INTEGER_T DIMDEC(xmacnew)

      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      INTEGER_T domlo(SDIM),domhi(SDIM)

      REAL_T mask(DIMV(mask))
      REAL_T maskcoef(DIMV(maskcoef))

      REAL_T ax(DIMV(ax))
      REAL_T vol(DIMV(vol))
      REAL_T xflux(DIMV(xflux),nflux_feed_forward)
      REAL_T velfab(DIMV(velfab),nstate)
      REAL_T lsfab(DIMV(lsfab),(SDIM+1)*nmat)
      REAL_T snew(DIMV(snew),nstate)
      REAL_T sold(DIMV(sold),nstate)
      REAL_T lsnew(DIMV(lsnew),(SDIM+1)*nmat)
      REAL_T lsold(DIMV(lsold),(SDIM+1)*nmat)
      REAL_T xmacnew(DIMV(xmacnew))

      INTEGER_T i,j,k,ii,jj,kk,veldir
      INTEGER_T nhalf
      REAL_T xstenMAC(-3:3,SDIM)
      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_left(-3:3,SDIM)
      REAL_T xsten_right(-3:3,SDIM)
      REAL_T xcell_cen(SDIM)
      REAL_T xcell_cen_left(SDIM)
      REAL_T xcell_cen_right(SDIM)
      REAL_T vol_left,vol_right,vol_cen
      REAL_T primitive(nflux_feed_forward)
      REAL_T conserve(nflux_feed_forward)
      REAL_T conserve_Left(nflux_feed_forward)
      REAL_T conserve_Right(nflux_feed_forward)
      REAL_T uLeft(nflux_feed_forward)
      REAL_T uRight(nflux_feed_forward)
      REAL_T fluxdiff(nflux_feed_forward)
      REAL_T LLF_flux(nflux_feed_forward)
      REAL_T divu
      REAL_T normal_vel
      INTEGER_T icomp,ithermal,irest,ils,igeom,lscomp,vofcomp
      INTEGER_T nc,im
      INTEGER_T iii,jjj,kkk
      REAL_T total_mass_left,total_mass_right,total_mass
      

      nhalf=3
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid cell to mac"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      nstate_test=(SDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1
      nflux_test=nstate_test+nmat*(SDIM+1)
      if (nstate.ne.nstate_test) then
       print *,"nstate invalid"
       stop
      endif
      if (nflux_feed_forward.ne.nflux_test) then
       print *,"nflux_feed_forward invalid"
       stop
      endif

      if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
       print *,"slab_step invalid FORT_FEED_FORWARD_ADVECT"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(mask),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(maskcoef),1,-1,234)

      call checkbound(fablo,fabhi,DIMS(ax),0,normdir,231)
      call checkbound(fablo,fabhi,DIMS(vol),1,-1,234)

      call checkbound(fablo,fabhi,DIMS(xflux),0,normdir,263)

      call checkbound(fablo,fabhi,DIMS(velfab),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(lsfab),1,-1,234)

      call checkbound(fablo,fabhi,DIMS(snew),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(sold),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(lsnew),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(lsold),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(xmacnew),0,normdir,264)

      igeom=(AMREX_SPACEDIM+1)+nmat*num_state_material
      ithermal=(AMREX_SPACEDIM+1)
      ils=(AMREX_SPACEDIM+1)+nmat*num_state_material+ &
       nmat*ngeom_raw+1

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
       print *,"normdir out of range in FORT_FEED_FORWARD_ADVECT"
       stop
      endif 

      if (init_fluxes.eq.1) then
       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0, &
               normdir,2)
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,normdir,1)
        call gridsten_level(xsten_left,i-ii,j-jj,k-kk,level,nhalf)
        call gridsten_level(xsten_right,i,j,k,level,nhalf)

        call CISBOX(xsten_left,nhalf, &
         xlo,dx,i-ii,j-jj,k-kk, &
         bfact,level, &
         vol_left,xcell_cen_left,SDIM)
        call CISBOX(xsten_right,nhalf, &
         xlo,dx,i,j,k, &
         bfact,level, &
         vol_right,xcell_cen_right,SDIM)

        do nc=1,nflux_feed_forward
         if (nc.le.nstate) then
          uLeft(nc)=velfab(D_DECL(i-ii,j-jj,k-kk),nc)
          uRight(nc)=velfab(D_DECL(i,j,k),nc)
         else if ((nc.gt.nstate).and.(nc.le.nflux_feed_forward)) then
          uLeft(nc)=lsfab(D_DECL(i-ii,j-jj,k-kk),nc-nstate)
          uRight(nc)=lsfab(D_DECL(i,j,k),nc-nstate)
         else
          print *,"nc invalid"
          stop
         endif
        enddo ! nc=1..nflux_feed_forward

        call convert_to_conserve(uLeft,conserve_Left,total_mass_left, &
         nstate,nflux_feed_forward,xcell_cen_left,nmat,normdir)
        call convert_to_conserve(uRight,conserve_Right,total_mass_right, &
         nstate,nflux_feed_forward,xcell_cen_right,nmat,normdir)
        if ((total_mass_left.le.zero).or. &
            (total_mass_right.le.zero)) then
         print *,"total_mass_left or total_mass_right invalid"
         stop
        endif
        vol_left=vol(D_DECL(i-ii,j-jj,k-kk))
        vol_right=vol(D_DECL(i,j,k))
        if ((vol_left.lt.zero).or.(vol_right.lt.zero)) then
         print *,"vol_left or vol_right invalid"
         stop
        endif
        if (vol_left+vol_right.le.zero) then
         print *,"vol_left or vol_right invalid"
         stop
        endif
        normal_vel=(vol_left*conserve_Left(normdir+1)+ &
                    vol_right*conserve_Right(normdir+1))/ &
                   (vol_left*total_mass_left+ &
                    vol_right*total_mass_right)
        call shallow_water_LLF(uLeft,uRight,normal_vel, &
         normdir,nstate,nflux_feed_forward,LLF_flux,nmat, &
         xcell_cen_left,xcell_cen_right)

         ! nstate component is normal_vel
        do nc=1,nflux_feed_forward
         xflux(D_DECL(i,j,k),nc)=LLF_flux(nc)
        enddo
   
       enddo
       enddo
       enddo

      else if (init_fluxes.eq.0) then

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        call gridsten_level(xsten,i,j,k,level,nhalf)

        call CISBOX(xsten,nhalf, &
         xlo,dx,i,j,k, &
         bfact,level, &
         vol_cen,xcell_cen,SDIM)

        do nc=1,nflux_feed_forward
         if (nc.le.nstate) then
          primitive(nc)=sold(D_DECL(i,j,k),nc)
         else if ((nc.gt.nstate).and.(nc.le.nflux_feed_forward)) then
          primitive(nc)=lsold(D_DECL(i,j,k),nc-nstate)
         else
          print *,"nc invalid"
          stop
         endif
        enddo ! nc=1..nflux_feed_forward

        call convert_to_conserve(primitive,conserve,total_mass, &
         nstate,nflux_feed_forward,xcell_cen,nmat,normdir)

        if (normdir.eq.0) then
         do nc=1,nflux_feed_forward
          if (nc.eq.nstate) then
           snew(D_DECL(i,j,k),nc)=sold(D_DECL(i,j,k),nc)
          else if (nc.lt.nstate) then
           snew(D_DECL(i,j,k),nc)=conserve(nc)
          else if ((nc.gt.nstate).and.(nc.le.nflux_feed_forward)) then
           lsnew(D_DECL(i,j,k),nc-nstate)=conserve(nc)
          else
           print *,"nc invalid"
           stop
          endif
         enddo ! nc=1..nflux_feed_forward
        else if ((normdir.ge.1).and.(normdir.lt.SDIM)) then
         ! do nothing
        else
         print *,"normdir invalid"
         stop
        endif

        do nc=1,nflux_feed_forward
         fluxdiff(nc)=dt* &
          (ax(D_DECL(i+ii,j+jj,k+kk))*xflux(D_DECL(i+ii,j+jj,k+kk),nc)- &
           ax(D_DECL(i,j,k))*xflux(D_DECL(i,j,k),nc))/vol(D_DECL(i,j,k))
        enddo

        divu=fluxdiff(nstate)

        do veldir=1,SDIM
         snew(D_DECL(i,j,k),veldir)=snew(D_DECL(i,j,k),veldir)- &
           fluxdiff(veldir)
        enddo
        do im=1,nmat
         icomp=ithermal+(im-1)*num_state_material+1 ! density
         snew(D_DECL(i,j,k),icomp)=snew(D_DECL(i,j,k),icomp)-fluxdiff(icomp)
         icomp=icomp+1 ! energy
         snew(D_DECL(i,j,k),icomp)=snew(D_DECL(i,j,k),icomp)-fluxdiff(icomp)
         do irest=1,num_state_material-2
          icomp=ithermal+(im-1)*num_state_material+2+irest
          snew(D_DECL(i,j,k),icomp)=snew(D_DECL(i,j,k),icomp)-fluxdiff(icomp)
         enddo
        enddo ! im=1..nmat

        do im=1,nmat
         vofcomp=igeom+(im-1)*ngeom_raw+1 ! volume fractions
         snew(D_DECL(i,j,k),vofcomp)=snew(D_DECL(i,j,k),vofcomp)- &
          fluxdiff(vofcomp)+divu*conserve(vofcomp)

         do veldir=1,SDIM
          snew(D_DECL(i,j,k),vofcomp+veldir)= &
           snew(D_DECL(i,j,k),vofcomp+veldir)- &
           fluxdiff(vofcomp+veldir)+divu*conserve(vofcomp+veldir)
          if (veldir.eq.normdir+1) then
           snew(D_DECL(i,j,k),vofcomp+veldir)= &
            snew(D_DECL(i,j,k),vofcomp+veldir)+ &
            fluxdiff(vofcomp)
          endif
         enddo ! veldir
         lscomp=ils+(im-1)*(AMREX_SPACEDIM+1)+1
         lsnew(D_DECL(i,j,k),lscomp-ils)= &
          lsnew(D_DECL(i,j,k),lscomp-ils)- &
          fluxdiff(lscomp)+divu*conserve(lscomp)
         do veldir=1,SDIM
          iii=0
          jjj=0
          kkk=0
          if (veldir.eq.1) then
           iii=1
          else if (veldir.eq.2) then
           jjj=1
          else if ((veldir.eq.SDIM).and.(SDIM.eq.3)) then
           kkk=1
          else
           print *,"veldir invalid"
           stop
          endif
           
          lsnew(D_DECL(i,j,k),lscomp-ils+veldir)= &
           lsnew(D_DECL(i,j,k),lscomp-ils+veldir)- &
           fluxdiff(lscomp+veldir)+divu*conserve(lscomp+veldir)
          lsnew(D_DECL(i,j,k),lscomp-ils+veldir)= &
           lsnew(D_DECL(i,j,k),lscomp-ils+veldir)- &
           dt* &
           (sold(D_DECL(i+iii,j+jjj,k+kkk),normdir+1)- &
            sold(D_DECL(i-iii,j-jjj,k-kkk),normdir+1))* &
           conserve(lscomp+normdir+1)/ &
           (xsten(1,veldir)-xsten(-1,veldir))
         enddo ! veldir=1..sdim
        enddo ! im=1..nmat
 
        if (normdir.eq.SDIM-1) then
         do nc=1,nflux_feed_forward
          if (nc.le.nstate) then
           conserve(nc)=snew(D_DECL(i,j,k),nc)
          else if ((nc.gt.nstate).and.(nc.le.nflux_feed_forward)) then
           conserve(nc)=lsnew(D_DECL(i,j,k),nc-nstate)
          else
           print *,"nc invalid"
           stop
          endif
         enddo ! nc=1..nflux_feed_forward
 
         call convert_to_primitive(primitive,conserve, &
           nstate,nflux_feed_forward,xcell_cen,nmat,normdir)

         do nc=1,nflux_feed_forward
          if (nc.eq.nstate) then
           snew(D_DECL(i,j,k),nc)=sold(D_DECL(i,j,k),nc)
          else if (nc.lt.nstate) then
           snew(D_DECL(i,j,k),nc)=primitive(nc)
          else if ((nc.gt.nstate).and.(nc.le.nflux_feed_forward)) then
           lsnew(D_DECL(i,j,k),nc-nstate)=primitive(nc)
          else
           print *,"nc invalid"
           stop
          endif
         enddo ! nc=1..nflux_feed_forward
        else if ((normdir.ge.0).and.(normdir.lt.SDIM-1)) then
         ! do nothing
        else
         print *,"normdir invalid"
         stop
        endif

       enddo
       enddo
       enddo

      else if (init_fluxes.eq.2) then

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0, &
               normdir,3)
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,normdir,2)
        call gridsten_level(xsten_left,i-ii,j-jj,k-kk,level,nhalf)
        call gridsten_level(xsten_right,i,j,k,level,nhalf)
        call CISBOX(xsten_left,nhalf, &
         xlo,dx,i-ii,j-jj,k-kk, &
         bfact,level, &
         vol_left,xcell_cen_left,SDIM)
        call CISBOX(xsten_right,nhalf, &
         xlo,dx,i,j,k, &
         bfact,level, &
         vol_right,xcell_cen_right,SDIM)

        do nc=1,nflux_feed_forward
         if (nc.le.nstate) then
          uLeft(nc)=velfab(D_DECL(i-ii,j-jj,k-kk),nc)
          uRight(nc)=velfab(D_DECL(i,j,k),nc)
         else if ((nc.gt.nstate).and.(nc.le.nflux_feed_forward)) then
          uLeft(nc)=lsfab(D_DECL(i-ii,j-jj,k-kk),nc-nstate)
          uRight(nc)=lsfab(D_DECL(i,j,k),nc-nstate)
         else
          print *,"nc invalid"
          stop
         endif
        enddo ! nc=1..nflux_feed_forward

        call convert_to_conserve(uLeft,conserve_Left,total_mass_left, &
         nstate,nflux_feed_forward,xcell_cen_left,nmat,normdir)
        call convert_to_conserve(uRight,conserve_Right,total_mass_right, &
         nstate,nflux_feed_forward,xcell_cen_right,nmat,normdir)
        if ((total_mass_left.le.zero).or. &
            (total_mass_right.le.zero)) then
         print *,"total_mass_left or total_mass_right invalid"
         stop
        endif
        vol_left=vol(D_DECL(i-ii,j-jj,k-kk))
        vol_right=vol(D_DECL(i,j,k))
        if ((vol_left.lt.zero).or.(vol_right.lt.zero)) then
         print *,"vol_left or vol_right invalid"
         stop
        endif
        if (vol_left+vol_right.le.zero) then
         print *,"vol_left or vol_right invalid"
         stop
        endif
        normal_vel=(vol_left*conserve_Left(normdir+1)+ &
                    vol_right*conserve_Right(normdir+1))/ &
                   (vol_left*total_mass_left+ &
                    vol_right*total_mass_right)
        xmacnew(D_DECL(i,j,k),1)=normal_vel
   
       enddo
       enddo
       enddo

      else
       print *,"init_fluxes invalid"
       stop
      endif

      return
      end subroutine FORT_FEED_FORWARD_ADVECT
