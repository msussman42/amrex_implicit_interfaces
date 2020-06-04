#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "PLIC_F.H"
#include "TECPLOTUTIL_F.H"
#include "MASS_TRANSFER_F.H"
#include "MOF_REDIST_F.H"
#include "SOLIDFLUID_F.H"
#include "AMReX_ArrayLim.H"


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

module tsat_module
implicit none

      INTEGER_T :: nburning
      INTEGER_T :: ntsat
      INTEGER_T :: ncomp_per_burning
      INTEGER_T :: ncomp_per_tsat
      INTEGER_T :: velflag

      INTEGER_T DIMDEC(burnvel)
      INTEGER_T DIMDEC(tsatfab)

      REAL_T, dimension(:,:,:), allocatable :: burnvel
      REAL_T, dimension(:,:,:), allocatable :: tsatfab

contains

end module tsat_module

      subroutine set_dimdec(DIMS(fabdim), &
                      fablo,fabhi,ngrow)
      IMPLICIT NONE

      INTEGER_T DIMDEC(fabdim)
      INTEGER_T fablo(SDIM)
      INTEGER_T fabhi(SDIM)
      INTEGER_T ngrow
      INTEGER_T dir

      dir=1
      ARG_L1(fabdim)=fablo(dir)-ngrow
      ARG_H1(fabdim)=fabhi(dir)+ngrow
      dir=2
      ARG_L2(fabdim)=fablo(dir)-ngrow
      ARG_H2(fabdim)=fabhi(dir)+ngrow
#if (AMREX_SPACEDIM==3)
      ARG_L3(fabdim)=fablo(dir)-ngrow
      ARG_H3(fabdim)=fabhi(dir)+ngrow
      print *,"prototype code only for 2d"
      stop
#elif (AMREX_SPACEDIM==2)
      ! do nothing
#else
      print *,"dimension bust"
      stop
#endif
      return
      end subroutine set_dimdec

      subroutine set_boundary_recon( &
        recon,DIMS(recon), &
        fablo,fabhi, &
        nmat,ngrow)
      USE probcommon_module 
      USE global_utility_module 
      IMPLICIT NONE

      INTEGER_T nmat,ngrow
      INTEGER_T fablo(SDIM)
      INTEGER_T fabhi(SDIM)
      INTEGER_T DIMDEC(recon)
      REAL_T recon(DIMV(recon),nmat*ngeom_recon) ! F,X,order,SL,I x nmat
      INTEGER_T i,j,iofs,jofs,imof,dir,dirtan

      if (ngrow.ge.0) then
       ! do nothing
      else
       print *,"ngrow invalid"
       stop
      endif
      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(recon),ngrow,-1,1221)

      if (SDIM.eq.2) then

        ! top and bottom
       dir=1
       dirtan=2
       do i=fablo(dir),fabhi(dir)
        do imof=1,ngeom_recon*nmat
         do jofs=-1,-ngrow,-1
          j=fablo(dirtan)
          recon(i,j+jofs,imof) = recon(i,j,imof)
         enddo
         do jofs=1,ngrow,1
          j=fabhi(dirtan)
          recon(i,j+jofs,imof) = recon(i,j,imof)
         enddo
        enddo ! imof=1..ngeom_recon*nmat
       enddo ! i=fablo,fabhi

        ! left and right
       dir=2
       dirtan=1
       do j=fablo(dir)-ngrow,fabhi(dir)+ngrow
        do imof=1,ngeom_recon*nmat
         do iofs=-1,-ngrow,-1
          i=fablo(dirtan)
          recon(i+iofs,j,imof) = recon(i,j,imof)
         enddo
         do iofs=1,ngrow,1
          i=fabhi(dirtan)
          recon(i+iofs,j,imof) = recon(i,j,imof)
         enddo
        enddo ! imof=1..ngeom_recon*nmat
       enddo ! j=fablo-ngrow,fabhi+ngrow

      else
       print *,"only sdim==2 supported"
       stop
      endif

      end subroutine set_boundary_recon

      subroutine set_boundary_burning( &
        vel,DIMS(vel), &
        fablo,fabhi, &
        nten,nmat,ncomp_per, &
        nburning,ngrow)
      USE probcommon_module 
      USE global_utility_module 
      IMPLICIT NONE

      INTEGER_T nten,nmat,nburning,ngrow,ncomp_per
      INTEGER_T fablo(SDIM)
      INTEGER_T fabhi(SDIM)
      INTEGER_T DIMDEC(vel)
      REAL_T vel(DIMV(vel),nburning) ! status first nten comp
      INTEGER_T i,j,iofs,jofs,ivel,dir,dirtan

      if (ngrow.ge.0) then
       ! do nothing
      else
       print *,"ngrow invalid"
       stop
      endif
      if (nten.eq.(((nmat-1)*(nmat-1)+nmat-1)/2)) then
       ! do nothing
      else
       print *,"nten invalid"
       stop
      endif
      if (nburning.eq.nten*(ncomp_per+1)) then
       ! do nothing
      else
       print *,"nburning invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(vel),ngrow,-1,1221)

      if (SDIM.eq.2) then

        ! top and bottom
       dir=1
       dirtan=2
       do i=fablo(dir),fabhi(dir)
        do ivel=1,nburning
         do jofs=-1,-ngrow,-1
          j=fablo(dirtan)
          vel(i,j+jofs,ivel) = vel(i,j,ivel)
         enddo
         do jofs=1,ngrow,1
          j=fabhi(dirtan)
          vel(i,j+jofs,ivel) = vel(i,j,ivel)
         enddo
        enddo ! ivel=1..nburning
       enddo ! i=fablo,fabhi

        ! left and right
       dir=2
       dirtan=1
       do j=fablo(dir)-ngrow,fabhi(dir)+ngrow
        do ivel=1,nburning
         do iofs=-1,-ngrow,-1
          i=fablo(dirtan)
          vel(i+iofs,j,ivel) = vel(i,j,ivel)
         enddo
         do iofs=1,ngrow,1
          i=fabhi(dirtan)
          vel(i+iofs,j,ivel) = vel(i,j,ivel)
         enddo
        enddo ! ivel=1..nburning
       enddo ! j=fablo-ngrow,fabhi+ngrow

      else
       print *,"only sdim==2 supported"
       stop
      endif

      end subroutine set_boundary_burning

      subroutine set_boundary_EOS( &
        EOS,DIMS(EOS), &
        fablo,fabhi, &
        nmat,ngrow)
      USE probcommon_module 
      USE global_utility_module 
      IMPLICIT NONE

      INTEGER_T nmat,ngrow
      INTEGER_T fablo(SDIM)
      INTEGER_T fabhi(SDIM)
      INTEGER_T DIMDEC(EOS)
      REAL_T EOS(DIMV(EOS),nmat*2) ! density, temperature x nmat
      INTEGER_T i,j,iofs,jofs,i_eos,dir,dirtan

      if (ngrow.ge.0) then
       ! do nothing
      else
       print *,"ngrow invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(EOS),ngrow,-1,1221)

      if (SDIM.eq.2) then

        ! top and bottom
       dir=1
       dirtan=2
       do i=fablo(dir),fabhi(dir)
        do i_eos=1,2*nmat
         do jofs=-1,-ngrow,-1
          j=fablo(dirtan)
          EOS(i,j+jofs,i_eos) = EOS(i,j,i_eos)
         enddo
         do jofs=1,ngrow,1
          j=fabhi(dirtan)
          EOS(i,j+jofs,i_eos) = EOS(i,j,i_eos)
         enddo
        enddo ! i_eos=1..2*nmat
       enddo ! i=fablo,fabhi

        ! left and right
       dir=2
       dirtan=1
       do j=fablo(dir)-ngrow,fabhi(dir)+ngrow
        do i_eos=1,2*nmat
         do iofs=-1,-ngrow,-1
          i=fablo(dirtan)
          EOS(i+iofs,j,i_eos) = EOS(i,j,i_eos)
         enddo
         do iofs=1,ngrow,1
          i=fabhi(dirtan)
          EOS(i+iofs,j,i_eos) = EOS(i,j,i_eos)
         enddo
        enddo ! i_eos=1..2*nmat
       enddo ! j=fablo-ngrow,fabhi+ngrow

      else
       print *,"only sdim==2 supported"
       stop
      endif

      end subroutine set_boundary_EOS


      subroutine set_boundary_LS( &
        LS,DIMS(LS), &
        fablo,fabhi, &
        nmat,ngrow)
      USE probcommon_module 
      USE global_utility_module 
      IMPLICIT NONE

      INTEGER_T nmat,ngrow
      INTEGER_T fablo(SDIM)
      INTEGER_T fabhi(SDIM)
      INTEGER_T DIMDEC(LS)
      REAL_T LS(DIMV(LS),nmat*(SDIM+1)) ! LS x nmat + slope x nmat
      INTEGER_T i,j,iofs,jofs,i_ls,dir,dirtan

      if (ngrow.ge.0) then
       ! do nothing
      else
       print *,"ngrow invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(LS),ngrow,-1,1221)

      if (SDIM.eq.2) then

        ! top and bottom
       dir=1
       dirtan=2
       do i=fablo(dir),fabhi(dir)
        do i_ls=1,nmat*(SDIM+1)
         do jofs=-1,-ngrow,-1
          j=fablo(dirtan)
          LS(i,j+jofs,i_ls) = LS(i,j,i_ls)
         enddo
         do jofs=1,ngrow,1
          j=fabhi(dirtan)
          LS(i,j+jofs,i_ls) = LS(i,j,i_ls)
         enddo
        enddo ! i_ls=1..nmat*(sdim+1)
       enddo ! i=fablo,fabhi

        ! left and right
       dir=2
       dirtan=1
       do j=fablo(dir)-ngrow,fabhi(dir)+ngrow
        do i_ls=1,nmat*(SDIM+1)
         do iofs=-1,-ngrow,-1
          i=fablo(dirtan)
          LS(i+iofs,j,i_ls) = LS(i,j,i_ls)
         enddo
         do iofs=1,ngrow,1
          i=fabhi(dirtan)
          LS(i+iofs,j,i_ls) = LS(i,j,i_ls)
         enddo
        enddo ! i_ls=1..nmat*(sdim+1)
       enddo ! j=fablo-ngrow,fabhi+ngrow

      else
       print *,"only sdim==2 supported"
       stop
      endif

      end subroutine set_boundary_LS


      subroutine set_boundary_VOF( &
        VOF,DIMS(VOF), &
        fablo,fabhi, &
        nmat,ngrow)
      USE probcommon_module 
      USE global_utility_module 
      IMPLICIT NONE

      INTEGER_T nmat,ngrow
      INTEGER_T fablo(SDIM)
      INTEGER_T fabhi(SDIM)
      INTEGER_T DIMDEC(VOF)
      REAL_T VOF(DIMV(VOF),nmat*ngeom_raw) 
      INTEGER_T i,j,iofs,jofs,i_vof,dir,dirtan

      if (ngrow.ge.0) then
       ! do nothing
      else
       print *,"ngrow invalid"
       stop
      endif
      if (ngeom_raw.eq.SDIM+1) then
       ! do nothing
      else
       print *,"ngeom_raw invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(VOF),ngrow,-1,1221)

      if (SDIM.eq.2) then

        ! top and bottom
       dir=1
       dirtan=2
       do i=fablo(dir),fabhi(dir)
        do i_vof=1,nmat*ngeom_raw
         do jofs=-1,-ngrow,-1
          j=fablo(dirtan)
          VOF(i,j+jofs,i_vof) = VOF(i,j,i_vof)
         enddo
         do jofs=1,ngrow,1
          j=fabhi(dirtan)
          VOF(i,j+jofs,i_vof) = VOF(i,j,i_vof)
         enddo
        enddo ! i_vof=1..nmat*(sdim+1)
       enddo ! i=fablo,fabhi

        ! left and right
       dir=2
       dirtan=1
       do j=fablo(dir)-ngrow,fabhi(dir)+ngrow
        do i_vof=1,nmat*ngeom_raw
         do iofs=-1,-ngrow,-1
          i=fablo(dirtan)
          VOF(i+iofs,j,i_vof) = VOF(i,j,i_vof)
         enddo
         do iofs=1,ngrow,1
          i=fabhi(dirtan)
          VOF(i+iofs,j,i_vof) = VOF(i,j,i_vof)
         enddo
        enddo ! i_vof=1..nmat*(sdim+1)
       enddo ! j=fablo-ngrow,fabhi+ngrow

      else
       print *,"only sdim==2 supported"
       stop
      endif

      end subroutine set_boundary_VOF


      subroutine set_boundary_FSI( &
        FSI_MF,DIMS(FSI_MF), &
        fablo,fabhi, &
        ngrow)
      USE probcommon_module 
      USE probmain_module 
      USE global_utility_module 
      IMPLICIT NONE

      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: fablo(SDIM)
      INTEGER_T, intent(in) :: fabhi(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(FSI_MF)
      REAL_T, intent(inout) :: FSI_MF(DIMV(FSI_MF),nFSI_all) 
      INTEGER_T ::  i,j,iofs,jofs,i_fsi,dir,dirtan

      if (ngrow.ge.0) then
       ! do nothing
      else
       print *,"ngrow invalid"
       stop
      endif

       !nparts x (velocity + LS + temperature + flag+stress)  3D
      if ((nFSI_all.eq.global_nparts*nFSI_sub).and. &
          (nFSI_sub.eq.12).and. &
          (ngrowFSI.eq.3).and. &
          (ngrow.eq.3)) then
       ! do nothing
       else
        print *,"nFSI_all,nFSI_sub, or ngrowFSI invalid"
        stop
       endif
      call checkbound(fablo,fabhi,DIMS(FSI_MF),ngrow,-1,1221)

      if (SDIM.eq.2) then

        ! top and bottom
       dir=1
       dirtan=2
       do i=fablo(dir),fabhi(dir)
        do i_fsi=1,nFSI_all
         do jofs=-1,-ngrow,-1
          j=fablo(dirtan)
          FSI_MF(i,j+jofs,i_fsi) = FSI_MF(i,j,i_fsi)
         enddo
         do jofs=1,ngrow,1
          j=fabhi(dirtan)
          FSI_MF(i,j+jofs,i_fsi) = FSI_MF(i,j,i_fsi)
         enddo
        enddo ! i_fsi=1..nFSI_all
       enddo ! i=fablo,fabhi

        ! left and right
       dir=2
       dirtan=1
       do j=fablo(dir)-ngrow,fabhi(dir)+ngrow
        do i_fsi=1,nFSI_all
         do iofs=-1,-ngrow,-1
          i=fablo(dirtan)
          FSI_MF(i+iofs,j,i_fsi) = FSI_MF(i,j,i_fsi)
         enddo
         do iofs=1,ngrow,1
          i=fabhi(dirtan)
          FSI_MF(i+iofs,j,i_fsi) = FSI_MF(i,j,i_fsi)
         enddo
        enddo ! i_fsi=1..nFSI_all
       enddo ! j=fablo-ngrow,fabhi+ngrow

      else
       print *,"only sdim==2 supported"
       stop
      endif

      end subroutine set_boundary_FSI

      subroutine get_exact_LS(xgrid,time,im,LS)
      USE probcommon_module 
      USE probmain_module 
      USE global_utility_module 
      USE supercooled_exact_sol
      USE variable_temperature_drop
      USE mmat_FVM, only: dist_concentric

      IMPLICIT NONE

      REAL_T xgrid(SDIM)
      REAL_T time
      INTEGER_T im
      REAL_T LS,stefan_time,test_radblob

      if (probtype.eq.3) then
       stefan_time=fort_time_radblob(2)+time
       call solidification_front_radius_driver(stefan_time,test_radblob)
       LS=sqrt((xgrid(1)-xblob)**2+(xgrid(2)-yblob)**2)-test_radblob
       if (im.eq.1) then
        LS=-LS
       else if (im.eq.2) then
        ! do nothing
       else
        print *,"im invalid"
        stop
       endif
      else if (probtype.eq.4) then
       test_radblob=axisymmetric_disk_radblob(2)
       LS=sqrt((xgrid(1)-xblob)**2+(xgrid(2)-yblob)**2)-test_radblob
       if (im.eq.1) then
        LS=-LS
       else if (im.eq.2) then
        ! do nothing
       else
        print *,"im invalid"
        stop
       endif
      else if (probtype.eq.400) then
       call dist_concentric(im,xgrid(1),xgrid(2),LS,probtype)
      else if (probtype.eq.401) then
       call dist_concentric(im,xgrid(1),xgrid(2),LS,probtype)
      else if (probtype.eq.402) then
       call dist_concentric(im,xgrid(1),xgrid(2),LS,probtype)
      else if (probtype.eq.5) then
       ! material 1: left  material 2: right
       LS=0.2d0+time-xblob
       if (im.eq.2) then
        LS=-LS
       else if (im.eq.1) then
        ! do nothing
       else
        print *,"im invalid"
        stop
       endif
      else if (probtype.eq.19) then
       call dist_concentric(im,xgrid(1),xgrid(2),LS,probtype)
      else if (probtype.eq.13) then
       call dist_concentric(im,xgrid(1),xgrid(2),LS,probtype)
      else if (probtype.eq.1) then
       call dist_concentric(im,xgrid(1),xgrid(2),LS,probtype)
      else
       print *,"probtype invalid in get_exact_LS"
       stop
      endif

      end subroutine get_exact_LS

      subroutine get_exact_NRM(xgrid,time,im,NRM)
      USE probcommon_module 
      USE probmain_module 
      USE global_utility_module 
      USE supercooled_exact_sol
      USE variable_temperature_drop
      USE mmat_FVM, only: dist_concentric

      IMPLICIT NONE

      REAL_T xgrid(SDIM)
      REAL_T time
      INTEGER_T im
      REAL_T NRM(SDIM)
      REAL_T mag
      REAL_T stefan_time,test_radblob
      REAL_T h_opt,LSp1,LSm1
      INTEGER_T dir

      if ((im.ge.1).and.(im.le.num_materials)) then
       ! do nothing
      else
       print *,"im invalid get_exact_NRM"
       stop
      endif

      do dir=1,SDIM
       NRM(dir)=zero
      enddo
      if (probtype.eq.3) then
       stefan_time=fort_time_radblob(2)+time
       call solidification_front_radius_driver(stefan_time,test_radblob)
       mag=sqrt((xgrid(1)-xblob)**2+(xgrid(2)-yblob)**2)
       if (mag.gt.zero) then
        NRM(1)=(xgrid(1)-xblob)/mag
        NRM(2)=(xgrid(2)-yblob)/mag
       endif
       if (im.eq.1) then
        do dir=1,SDIM
         NRM(dir)=-NRM(dir)
        enddo
       else if (im.eq.2) then
        ! do nothing
       else
        print *,"im invalid"
        stop
       endif
      else if (probtype.eq.4) then
       ! inside material 1, outside material 2
       test_radblob=axisymmetric_disk_radblob(2)
       mag=sqrt((xgrid(1)-xblob)**2+(xgrid(2)-yblob)**2)
       if (mag.gt.zero) then
        NRM(1)=(xgrid(1)-xblob)/mag
        NRM(2)=(xgrid(2)-yblob)/mag
       endif
       if (im.eq.1) then
        ! material 1 normal points into material 1
        do dir=1,SDIM
         NRM(dir)=-NRM(dir)
        enddo
       else if (im.eq.2) then
        ! do nothing: material 2 normal points into material 2.
       else
        print *,"im invalid"
        stop
       endif

      else if (probtype.eq.5) then 
       ! im=1 left  im=2 right
       NRM(2)=0.0d0
       if (im.eq.1) then
        NRM(1)=-1.0d0
       else if (im.eq.2) then
        NRM(1)=1.0d0
       else
        print *,"im invalid"
        stop
       endif

      else if ((probtype.eq.19).or. &
               (probtype.eq.13).or. &
               (probtype.eq.1).or. &
               (probtype.eq.400).or. &
               (probtype.eq.401).or. &
               (probtype.eq.402)) then
       h_opt=1.0D-6
       dir=1
       call dist_concentric(im,xgrid(1)+h_opt,xgrid(2),LSp1,probtype)
       call dist_concentric(im,xgrid(1)-h_opt,xgrid(2),LSm1,probtype)
       NRM(dir)=(LSp1-LSm1)/(2.0d0*h_opt)
       dir=2
       call dist_concentric(im,xgrid(1),xgrid(2)+h_opt,LSp1,probtype)
       call dist_concentric(im,xgrid(1),xgrid(2)-h_opt,LSm1,probtype)
       NRM(dir)=(LSp1-LSm1)/(2.0d0*h_opt)
       mag=0.0d0
       do dir=1,SDIM
        mag=mag+NRM(dir)**2
       enddo
       mag=sqrt(mag)
       if (mag.gt.0.0d0) then
        do dir=1,SDIM
         NRM(dir)=NRM(dir)/mag
        enddo
       else if (mag.eq.0.0d0) then
        ! do nothing
       else
        print *,"mag invalid"
        stop
       endif
      else
       print *,"probtype invalid in get_exact_NRM"
       stop
      endif

      end subroutine get_exact_NRM

       ! probmain_module has global_nten
       ! probcommon_module has num_materials
      subroutine get_exact_VEL(xgrid,dx,time,iten,VEL)
      USE probcommon_module 
      USE probmain_module 
      USE global_utility_module 
      USE supercooled_exact_sol
      USE variable_temperature_drop
      USE GeneralClass
      IMPLICIT NONE

      REAL_T, intent(in) :: xgrid(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: iten
      REAL_T LS
      REAL_T, intent(out) :: VEL(SDIM)
      REAL_T NRM(SDIM)
      REAL_T xgrid_probe(SDIM)
      REAL_T T_grid,T_probe
      REAL_T h_opt
      REAL_T stefan_time
      REAL_T front_vel,test_front_vel
      REAL_T rgrid,dux
      INTEGER_T diflag
      INTEGER_T dir
      INTEGER_T im
      INTEGER_T local_nten
      REAL_T radial_slope
      REAL_T r1,r2
        ! declared in multimat_FVM.F90
      REAL_T, external   :: exact_temperature

      local_nten=( (num_materials-1)*(num_materials-1)+num_materials-1 )/2
      if (local_nten.eq.global_nten) then
       ! do nothing
      else
       print *,"local_nten or global_nten invalid"
       stop
      endif
      if (probtype.eq.3) then
       if ((iten.eq.1).or.(iten.eq.global_nten+1)) then
        stefan_time=fort_time_radblob(2)+time
        call solidification_front_speed_driver(stefan_time,front_vel)
        LS=sqrt((xgrid(1)-xblob)**2+(xgrid(2)-yblob)**2)
        if (LS.gt.0.0d0) then
         VEL(1)=front_vel*(xgrid(1)-xblob)/LS
         VEL(2)=front_vel*(xgrid(2)-yblob)/LS
        else if (LS.eq.0.0d0) then
         VEL(1)=0.0d0
         VEL(2)=0.0d0
        else
         print *,"LS invalid"
         stop
        endif
       else
        print *,"iten invalid (get exact vel) "
        print *,"iten=",iten
        print *,"probtype=",probtype
        stop
       endif
      else if (probtype.eq.4) then
       if ((iten.eq.1).or.(iten.eq.global_nten+1)) then
        call interp_LS_vel_to_grid(xgrid,2,LS,VEL)
        test_front_vel=0.0d0
        do dir=1,SDIM
         test_front_vel=test_front_vel+VEL(dir)**2
        enddo
        test_front_vel=sqrt(test_front_vel)
        call disk_get_speed(2,front_vel)
        if (1.eq.0) then
         print *,"front_vel,test_front_vel ",front_vel,test_front_vel
        endif
       else
        print *,"iten invalid (get exact vel 2)"
        print *,"iten=",iten
        print *,"probtype=",probtype
        stop
       endif
      else if (probtype.eq.400) then

       if (iten.eq.1) then
        test_front_vel=0.0d0
        VEL(1)=0.0d0
        VEL(2)=0.0d0
       else
        print *,"iten invalid (get exact vel 2)"
        print *,"iten=",iten
        print *,"probtype=",probtype
        stop
       endif

      else if (probtype.eq.401) then
        ! 1=liquid 2=gas 3=ice
        ! 12 13 23 21 31 32
        ! 31 -> ice to liquid
        ! 12 -> liquid to gas
       if ((iten.eq.1).or. &
           (iten.eq.2).or. &
           (iten.eq.3).or. &
           (iten.eq.4).or. &
           (iten.eq.5).or. &
           (iten.eq.6)) then
        test_front_vel=0.0d0
        VEL(1)=0.0d0
        VEL(2)=0.0d0
       else
        print *,"iten invalid (get exact vel 2)"
        print *,"iten=",iten
        print *,"probtype=",probtype
        stop
       endif

      else if (probtype.eq.402) then
        ! 1=liquid 2=gas 3=SUBSTRATE
        ! 12 13 23 21 31 32
        ! 31 -> ice to liquid
        ! 12 -> liquid to gas
       if ((iten.eq.1).or. &
           (iten.eq.2).or. &
           (iten.eq.3).or. &
           (iten.eq.4).or. &
           (iten.eq.5).or. &
           (iten.eq.6)) then
        test_front_vel=0.0d0
        VEL(1)=0.0d0
        VEL(2)=0.0d0
       else
        print *,"iten invalid (get exact vel 2)"
        print *,"iten=",iten
        print *,"probtype=",probtype
        stop
       endif

      else if (probtype.eq.5) then

       if ((iten.eq.1).or.(iten.eq.global_nten+1)) then
        VEL(1)=1.0d0
        VEL(2)=0.0d0
        test_front_vel=0.0d0
        do dir=1,SDIM
         test_front_vel=test_front_vel+VEL(dir)**2
        enddo
        test_front_vel=sqrt(test_front_vel)
        front_vel=1.0d0
        if (1.eq.0) then
         print *,"front_vel,test_front_vel ",front_vel,test_front_vel
        endif
       else
        print *,"iten invalid (get exact vel)"
        print *,"iten=",iten
        print *,"probtype=",probtype
        stop
       endif
      else if (probtype.eq.19) then
       rgrid=sqrt((xgrid(1)-pcenter(1))**2+(xgrid(2)-pcenter(2))**2)
       if (abs(rgrid-rlo).le.abs(rgrid-rhi)) then
        diflag=1  ! inner boundary
       else if (abs(rgrid-rlo).ge.abs(rgrid-rhi)) then
        diflag=2  ! outer boundary
       else
        print *,"rgrid invalid"
        stop
       endif

        ! materials 1,2,3 interfaces 12 13 23
       if ((iten.eq.1).or.(iten.eq.3)) then
        im=2
        ! Np,Mp,upolar,pcenter,rlo,rhi defined in vof_cisl.F90
        call find_polar_cart_inter(Np,Mp,upolar,pcenter,rlo,rhi, &
         xgrid,diflag,dux)
        call get_exact_NRM(xgrid,time,im,NRM)
        do dir=1,SDIM
         VEL(dir)=abs(dux)*NRM(dir)
        enddo
       else if (iten.eq.2) then
        do dir=1,SDIM
         VEL(dir)=zero
        enddo
       else
        print *,"iten invalid"
        stop
       endif
      else if ((probtype.eq.1).or.(probtype.eq.13)) then
       im=2
       if (probtype.eq.1) then
        dux=0.0d0
        r1=radcen-radeps
        r2=radcen+radeps
        if (r2.gt.r1) then
         ! do nothing
        else
         print *,"r2 or r1 invalid"
         stop
        endif
        radial_slope=1.0d0/(r2-r1)

        if (radial_variation.eq.1) then
         ! do nothing
        else if (radial_variation.eq.0) then
         radial_slope=0.0d0
        else
         print *,"radial_variation invalid"
         stop
        endif

        dux=radial_slope

       else if (probtype.eq.13) then
        call get_exact_NRM(xgrid,time,im,NRM)
        h_opt=1.0D-6
        if (h_opt.gt.0.0d0) then
         do dir=1,SDIM
          xgrid_probe(dir)=xgrid(dir)+h_opt*NRM(dir)
         enddo
         T_grid=exact_temperature(xgrid,time,im,probtype,num_materials, &
                fort_heatviscconst)
         T_probe=exact_temperature(xgrid_probe,time,im,probtype, &
                num_materials,fort_heatviscconst)
         dux=abs(T_grid-T_probe)/h_opt
        else
         print *,"h_opt invalid"
         stop
        endif
       else
        print *,"probtype invalid"
        stop
       endif

       call get_exact_NRM(xgrid,time,im,NRM)
       do dir=1,SDIM
        VEL(dir)=abs(dux)*NRM(dir)
       enddo
      else
       print *,"probtype invalid in get_exact_VEL"
       stop
      endif

      end subroutine get_exact_VEL

      subroutine deallocate_FSI()
      USE probcommon_module 
      USE probmain_module 
      USE global_utility_module 
      IMPLICIT NONE

      INTEGER_T ilev

      if (global_nparts.gt.0) then
       do ilev=0,cache_max_level
        deallocate(MG(ilev)%FSI_MF)
        deallocate(MG(ilev)%MASK_NBR_MF)
       enddo
       deallocate(MG)
       deallocate(im_solid_map)
      else if (global_nparts.eq.0) then
       ! do nothing
      else
       print *,"global_nparts invalid"
       stop
      endif

      return
      end subroutine deallocate_FSI


      subroutine ns_header_msg_level(ilev, &
         FSI_operation, &
         FSI_sub_operation, &
         iter)
      USE probcommon_module 
      USE probmain_module 
      USE global_utility_module 
      use MOF_routines_module
      use solidfluid_cpp_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: FSI_operation
      INTEGER_T, intent(in) :: FSI_sub_operation
      INTEGER_T, intent(in) :: ilev
      INTEGER_T, intent(in) :: iter

      INTEGER_T DIMDEC(unitdata)
      INTEGER_T DIMDEC(FSI_MF)

      INTEGER_T local_domlo(SDIM)
      INTEGER_T local_domhi(SDIM)
      INTEGER_T unitlo(SDIM)
      INTEGER_T unithi(SDIM)
      INTEGER_T ngrowFSI_unitfab

      REAL_T problo(SDIM)
      REAL_T probhi(SDIM)
      REAL_T problen(SDIM)

      INTEGER_T velbc(global_nparts*SDIM*SDIM*2)
      INTEGER_T vofbc(SDIM*2)

      REAL_T local_dx(SDIM)
      REAL_T dx_max_level(SDIM)
      REAL_T h_small

      INTEGER_T dir
      INTEGER_T bfact
      INTEGER_T tid,gridno,nthreads,tilenum

      REAL_T start_time,start_dt

      INTEGER_T FSI_refine_factor(num_materials)
      INTEGER_T FSI_bounding_box_ngrow(num_materials)
      INTEGER_T current_step,plot_interval,ioproc
      INTEGER_T partid,im_part,im,ibase,isub
      INTEGER_T klo,khi
      INTEGER_T i,j,k
      INTEGER_T ic,jc,kc
      REAL_T coarse_data

      REAL_T, allocatable, dimension(D_DECL(:,:,:),:) :: unitdata

      if ((ilev.ge.0).and.(ilev.le.cache_max_level)) then
       ! do nothing
      else
       print *,"ilev invalid"
       stop
      endif

      current_step=0
      plot_interval=-1
      ioproc=1

      do im=1,num_materials
       FSI_refine_factor(im)=1
       FSI_bounding_box_ngrow(im)=3
      enddo

      start_time=0.0d0
      start_dt=1.0d0
 
      bfact=1

      problo(1)=problox
      problo(2)=probloy
      probhi(1)=probhix
      probhi(2)=probhiy
      if (SDIM.eq.3) then
       problo(SDIM)=probloz
       probhi(SDIM)=probhiz
      endif

      do dir=1,SDIM
       problen(dir)=probhi(dir)-problo(dir)
       if (problen(dir).gt.0.0d0) then
        ! do nothing
       else
        print *,"problen invalid"
        stop
       endif
      enddo

      do dir=1,SDIM
       dx_max_level(dir)=dxlevel(cache_max_level,dir)
       local_dx(dir)=dxlevel(ilev,dir)
       local_domlo(dir)=domlo_level(ilev,dir)
       local_domhi(dir)=domhi_level(ilev,dir)
      enddo

      call set_dimdec(DIMS(FSI_MF),local_domlo, &
             local_domhi,ngrowFSI)

      if (FSI_operation.eq.0) then !init node locations;generate_new_triangles
       if ((iter.eq.0).and.(FSI_sub_operation.eq.0)) then
        ! do nothing
       else
        print *,"parameters invalid"
        stop
       endif
      else if (FSI_operation.eq.1) then !update node locations
       if ((iter.eq.0).and.(FSI_sub_operation.eq.0)) then
        ! do nothing
       else
        print *,"parameters invalid"
        stop
       endif
      else if (FSI_operation.eq.2) then !make distance in narrow band
       if ((iter.eq.0).and.(FSI_sub_operation.eq.0)) then
        ! do nothing
       else
        print *,"parameters invalid"
        stop
       endif
      else if (FSI_operation.eq.3) then !update the sign
       if ((iter.ge.0).and.(FSI_sub_operation.eq.0)) then
        ! do nothing
       else
        print *,"parameters invalid"
        stop
       endif
      else
       print *,"FSI_operation out of range"
       stop
      endif

      if (iter.eq.0) then
       FSI_touch_flag=0
      else if (iter.gt.0) then
       FSI_touch_flag=0
      else
       print *,"iter invalid"
       stop
      endif

      h_small=dxlevel(cache_max_level,1)

      if (FSI_operation.eq.0) then ! init node locations
       if (ilev.eq.0) then
        elements_generated=0
       else
        elements_generated=1
       endif
      else if (FSI_operation.eq.1) then ! update node locations
       if (ilev.eq.0) then
        elements_generated=0
       else
        elements_generated=1
       endif
      else if ((FSI_operation.ge.2).and.(FSI_operation.le.3)) then
       elements_generated=1
      else
       print *,"FSI_operation invalid"
       stop
      endif
 
      if (nFSI_sub.eq.12) then
       ! do nothing
      else
       print *,"nFSI_sub invalid ns_header_msg_level: ",nFSI_sub
       stop
      endif

      do dir=1,global_nparts*SDIM*2*SDIM
       velbc(dir)=REFLECT_EVEN
      enddo
      do dir=1,2*SDIM
       vofbc(dir)=REFLECT_EVEN
      enddo

      tid=0
      gridno=0
      tilenum=0
      nthreads=1

      if ((FSI_operation.eq.0).or. & ! initialize nodes
          (FSI_operation.eq.1)) then ! update node locations

       if (FSI_sub_operation.eq.0) then
        ! do nothing
       else
        print *,"FSI_sub_operation!=0"
        stop
       endif

       if (elements_generated.eq.0) then ! on the coarsest level
        ngrowFSI_unitfab=0
        do dir=1,SDIM
         unitlo(dir)=0
         unithi(dir)=0
        enddo

        call set_dimdec(DIMS(unitdata),unitlo,unithi,ngrowFSI_unitfab)
        allocate(unitdata(DIMV(unitdata),nFSI_all))

        call FORT_HEADERMSG( &
          tid, &
          tilenum, &
          gridno, &
          nthreads, &
          ilev, &
          cache_max_level, &
          cache_max_level, &
          FSI_operation, & ! 0 or 1 (initialize or update nodes)
          FSI_sub_operation, & ! 0
          unitlo,unithi, &
          unitlo,unithi, &
          bfact, &
          problo, &
          problen,  &
          dx_max_level, &
          problo, &
          probhi, &
          velbc, & 
          vofbc, &
          unitdata, & ! placeholder
          DIMS(unitdata), &
          unitdata, & ! velfab spot
          DIMS(unitdata), &
          unitdata, & ! mnbrfab spot
          DIMS(unitdata), &
          unitdata, & ! mfiner spot
          DIMS(unitdata), &
          nFSI_all, &
          nFSI_sub, &
          ngrowFSI_unitfab, &
          global_nparts, &
          im_solid_map, & ! type: 0..nmat-1
          h_small, &
          start_time,  &
          start_dt, &
          FSI_refine_factor, &
          FSI_bounding_box_ngrow, &
          FSI_touch_flag, &
          CTML_FSI_init, &
          CTML_force_model, &
          iter, &
          current_step, &
          plot_interval, &
          ioproc)

        elements_generated=1
        deallocate(unitdata)

       else if (elements_generated.eq.1) then
        ! do nothing
       else 
        print *,"elements_generated invalid"
        stop
       endif
       elements_generated=1

       CTML_FSI_init=1

      else if ((FSI_operation.eq.2).or. & ! make distance in narrow band
               (FSI_operation.eq.3)) then ! update the sign

       if ((FSI_sub_operation.eq.0).and.(ngrowFSI.eq.3)) then
        ! do nothing
       else
        print *,"FSI_sub_operation!=0 or ngrowFSI!=3"
        stop
       endif

       elements_generated=1;

        !nparts x (velocity + LS + temperature + flag+stress)  3D
       if ((nFSI_all.eq.global_nparts*nFSI_sub).and. &
           (nFSI_sub.eq.12).and. &
           (ngrowFSI.eq.3)) then
        ! do nothing
       else
        print *,"nFSI_all,nFSI_sub, or ngrowFSI invalid"
        stop
       endif
       
       if (FSI_operation.eq.2) then ! make distance in narrow band.

        ! fill coarse patch 
        if (ilev.gt.0) then

         do partid=1,global_nparts

          im_part=im_solid_map(partid)+1

          if ((im_part.ge.1).and. &
              (im_part.le.num_materials)) then

           if (FSI_flag(im_part).eq.7) then 

            if (SDIM.eq.3) then
             klo=local_domlo(SDIM)
             khi=local_domhi(SDIM)
            else if (SDIM.eq.2) then
             klo=0
             khi=0
            else
             print *,"dimension bust"
             stop
            endif

            do i=local_domlo(1),local_domhi(1)
            do j=local_domlo(2),local_domhi(2)
            do k=klo,khi
             ic=i/2
             jc=j/2
             kc=k/2
             ibase=(partid-1)*nFSI_sub
             do isub=1,nFSI_sub
              coarse_data=MG(ilev-1)%FSI_MF(D_DECL(ic,jc,kc),ibase+isub)
              if (((isub.ge.1).and.(isub.le.3)).or. & ! velocity
                  (isub.eq.4).or. & ! LS
                  (isub.eq.5).or. & ! temperature
                  ((isub.ge.7).and.(isub.le.12))) then
               MG(ilev)%FSI_MF(D_DECL(i,j,k),ibase+isub)=coarse_data
              else if (isub.eq.6) then ! flag
               if (ilev.gt.0) then
                MG(ilev)%FSI_MF(D_DECL(i,j,k),ibase+isub)=10.0d0
               else
                print *,"ilev invalid; flag=0 if ilev=0 and t=0"
                stop
               endif
              else
               print *,"isub invalid"
               stop
              endif
             enddo ! isub=1..nFSI_sub
            enddo !k
            enddo !j
            enddo !i
           else
            print *,"FSI_flag invalid"
            stop
           endif

          else 
           print *,"im_part invalid"
           stop
          endif

         enddo ! partid=1..global_nparts
  
        else if (ilev.eq.0) then
         ! do nothing
        else
         print *,"ilev invalid"
         stop
        endif

       else if (FSI_operation.eq.3) then ! update sign
        ! do not fill coarse patch.
       else
        print *,"FSI_operation invalid"
        stop
       endif

       call set_boundary_FSI( &
         MG(ilev)%FSI_MF, &
         DIMS(FSI_MF), &
         local_domlo,local_domhi, &
         ngrowFSI);

        ! MG(ilev)%MASK_NBR_MF 
        ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
        ! (2) =1 interior  =0 otherwise
        ! (3) =1 interior+ngrow-1  =0 otherwise
        ! (4) =1 interior+ngrow    =0 otherwise
 
       call FORT_HEADERMSG( &
          tid, &
          tilenum, &
          gridno, &
          nthreads, &
          ilev, &
          cache_max_level, &
          cache_max_level, &
          FSI_operation, & ! 2 or 3 (make distance or update sign)
          FSI_sub_operation, & ! 0
          local_domlo,local_domhi, &
          local_domlo,local_domhi, &
          bfact, &
          problo, &
          local_dx,  &
          dx_max_level, &
          problo, &
          probhi, &
          velbc, & 
          vofbc, &
          MG(ilev)%FSI_MF, & 
          DIMS(FSI_MF), &
          MG(ilev)%FSI_MF, &  ! velfab spot
          DIMS(FSI_MF), &
          MG(ilev)%MASK_NBR_MF, & 
          DIMS(FSI_MF), &
          MG(ilev)%MASK_NBR_MF, &  ! mfiner spot
          DIMS(FSI_MF), &
          nFSI_all, &
          nFSI_sub, &
          ngrowFSI, &
          global_nparts, &
          im_solid_map, & ! type: 0..nmat-1
          h_small, &
          start_time,  &
          start_dt, &
          FSI_refine_factor, &
          FSI_bounding_box_ngrow, &
          FSI_touch_flag, &
          CTML_FSI_init, &
          CTML_force_model, &
          iter, &
          current_step, &
          plot_interval, &
          ioproc)

       call set_boundary_FSI( &
         MG(ilev)%FSI_MF, &
         DIMS(FSI_MF), &
         local_domlo,local_domhi, &
         ngrowFSI);

      else
       print *,"FSI_operation invalid"
       stop
      endif

      return
      end subroutine ns_header_msg_level

      subroutine convert_lag_to_eul( &
          cache_max_level_in,sdim_in)
      USE probcommon_module 
      USE probmain_module 
      USE global_utility_module 
      use MOF_routines_module
      use solidfluid_cpp_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: cache_max_level_in
      INTEGER_T, intent(in) :: sdim_in

      INTEGER_T im,partid,iter,first_iter,ilev
      INTEGER_T FSI_operation,FSI_sub_operation

      REAL_T start_time,start_dt
      INTEGER_T local_domlo(SDIM)
      INTEGER_T local_domhi(SDIM)
      INTEGER_T klo,khi
      INTEGER_T i,j,k
      INTEGER_T ibase
      INTEGER_T dir
      INTEGER_T interior_flag
      INTEGER_T big_interior_flag
      INTEGER_T tilelo_array(SDIM)
      INTEGER_T tilehi_array(SDIM)
      REAL_T dx_current(SDIM)
      REAL_T dx_max_level(SDIM)
      REAL_T xlo_array(SDIM)
      REAL_T xhi_array(SDIM)
      INTEGER_T gridno_array(100) ! 100 is a placeholder
      INTEGER_T num_tiles_on_thread_proc(100)
      INTEGER_T max_num_tiles_on_thread_proc
      INTEGER_T tile_dim
      INTEGER_T num_grids_on_level
      INTEGER_T local_nthreads

      INTEGER_T DIMDEC(FSI_MF)

      if (sdim_in.eq.SDIM) then
       ! do nothing
      else
       print *,"sdim_in invalid"
       stop
      endif
      if (cache_max_level_in.eq.cache_max_level) then
       ! do nothing
      else
       print *,"cache_max_level_in invalid"
       stop
      endif

      CTML_FSI_init=0
      CTML_FSI_numsolids=0
      CTML_force_model=0

      global_nparts=0
      do im=1,num_materials
       if (FSI_flag(im).eq.7) then
        global_nparts=global_nparts+1
       else if (FSI_flag(im).eq.0) then
        ! do nothing
       else
        print *,"FSI_flag(im) value unexpected"
        stop
       endif
      enddo

      if ((global_nparts.gt.0).and. &
          (global_nparts.le.num_materials)) then
       allocate(im_solid_map(global_nparts))
       partid=0
       do im=1,num_materials
        if (FSI_flag(im).eq.7) then
         partid=partid+1
         im_solid_map(partid)=im-1
        else if (FSI_flag(im).eq.0) then
         ! do nothing
        else
         print *,"FSI_flag(im) value unexpected"
         stop
        endif
       enddo
       if (partid.ne.global_nparts) then
        print *,"partid.ne.global_nparts"
        stop
       endif

        !nparts x (velocity + LS + temperature + flag+stress)  3D
       nFSI_sub=12
       nFSI_all=global_nparts*nFSI_sub
       ngrowFSI=3

       start_time=0.0d0
       start_dt=1.0d0

       iter=0
        !initialize node locations; generate_new_triangles
       FSI_operation=0
       FSI_sub_operation=0
       ilev=0
       call ns_header_msg_level(ilev, &
              FSI_operation, &
              FSI_sub_operation, &
              iter)

       allocate(MG(0:cache_max_level))

       do ilev=0,cache_max_level

        !  new_localMF(FSI_MF,nFSI_all,ngrowFSI,-1);

        do dir=1,SDIM
         local_domlo(dir)=domlo_level(ilev,dir)
         local_domhi(dir)=domhi_level(ilev,dir)
        enddo

        call set_dimdec(DIMS(FSI_MF),local_domlo, &
                local_domhi,ngrowFSI)

        if (SDIM.eq.3) then
         klo=local_domlo(SDIM)-ngrowFSI
         khi=local_domhi(SDIM)+ngrowFSI
        else if (SDIM.eq.2) then
         klo=0
         khi=0
        else
         print *,"dimension bust"
         stop
        endif

        allocate(MG(ilev)%FSI_MF(DIMV(FSI_MF),nFSI_all))
        allocate(MG(ilev)%MASK_NBR_MF(DIMV(FSI_MF),4))

        ! MG(ilev)%MASK_NBR_MF 
        ! (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
        ! (2) =1 interior  =0 otherwise
        ! (3) =1 interior+ngrow-1  =0 otherwise
        ! (4) =1 interior+ngrow    =0 otherwise
        do i=ARG_L1(FSI_MF),ARG_H1(FSI_MF)
        do j=ARG_L2(FSI_MF),ARG_H2(FSI_MF)
        do k=klo,khi
         interior_flag=1
         if ((i.lt.local_domlo(1)).or. &
             (j.lt.local_domlo(2)).or. &
             (i.gt.local_domhi(1)).or. &
             (j.gt.local_domhi(2))) then
          interior_flag=0
         endif
         big_interior_flag=1
         if ((i.lt.local_domlo(1)-ngrowFSI+1).or. &
             (j.lt.local_domlo(2)-ngrowFSI+1).or. &
             (i.gt.local_domhi(1)+ngrowFSI-1).or. &
             (j.gt.local_domhi(2)+ngrowFSI-1)) then
          big_interior_flag=0
         endif
         if (SDIM.eq.3) then
          if ((k.lt.local_domlo(SDIM)).or. &
              (k.gt.local_domhi(SDIM))) then
           interior_flag=0
          endif
          if ((k.lt.local_domlo(SDIM)-ngrowFSI+1).or. &
              (k.gt.local_domhi(SDIM)+ngrowFSI-1)) then
           big_interior_flag=0
          endif
         else if (SDIM.eq.2) then
          ! do nothing
         else
          print *,"dimension bust"
          stop
         endif

         MG(ilev)%MASK_NBR_MF(D_DECL(i,j,k),1)=0.0d0
         if (interior_flag.eq.1) then
          MG(ilev)%MASK_NBR_MF(D_DECL(i,j,k),1)=1.0d0
         endif
         MG(ilev)%MASK_NBR_MF(D_DECL(i,j,k),2)= &
           MG(ilev)%MASK_NBR_MF(D_DECL(i,j,k),1)
         MG(ilev)%MASK_NBR_MF(D_DECL(i,j,k),3)=0.0d0
         if (big_interior_flag.eq.1) then
          MG(ilev)%MASK_NBR_MF(D_DECL(i,j,k),3)=1.0d0
         endif
         MG(ilev)%MASK_NBR_MF(D_DECL(i,j,k),4)=1.0d0
        enddo ! k
        enddo ! j
        enddo ! i

        do partid=1,global_nparts
         ibase=(partid-1)*nFSI_sub

         do i=ARG_L1(FSI_MF),ARG_H1(FSI_MF)
         do j=ARG_L2(FSI_MF),ARG_H2(FSI_MF)
         do k=klo,khi
          do dir=1,3
           MG(ilev)%FSI_MF(D_DECL(i,j,k),ibase+dir)=0.0d0 ! velocity
          enddo
          MG(ilev)%FSI_MF(D_DECL(i,j,k),ibase+3)=-9999.0d0 ! LS
          MG(ilev)%FSI_MF(D_DECL(i,j,k),ibase+4)=0.0d0 ! temp
          MG(ilev)%FSI_MF(D_DECL(i,j,k),ibase+5)=0.0d0 ! mask
          do dir=1,6
           MG(ilev)%FSI_MF(D_DECL(i,j,k),ibase+5+dir)=0.0d0 ! stress
          enddo
         enddo
         enddo
         enddo
        enddo ! partid=1..global_nparts

        do dir=1,SDIM
         tilelo_array(dir)=domlo_level(ilev,dir)
         tilehi_array(dir)=domhi_level(ilev,dir)
         dx_current(dir)=dxlevel(ilev,dir)
         dx_max_level(dir)=dxlevel(cache_max_level,dir)
        enddo
        xlo_array(1)=problox
        xlo_array(2)=probloy
        if (SDIM.eq.3) then
         xlo_array(SDIM)=probloz
        endif
        xhi_array(1)=probhix
        xhi_array(2)=probhiy
        if (SDIM.eq.3) then
         xhi_array(SDIM)=probhiz
        endif
        num_grids_on_level=1
        gridno_array(1)=0
        num_tiles_on_thread_proc(1)=1
        max_num_tiles_on_thread_proc=1
        tile_dim=1
        local_nthreads=1

        ! 1. create lagrangian container data structure within the 
        !    fortran part that recognizes tiles. 
        !    (FILLCONTAINER in SOLIDFLUID.F90)
        ! 2. fill the containers with the Lagrangian information.
        !    (CLSVOF_FILLCONTAINER called from FILLCONTAINER)
        !    i.e. associate to each tile a set of Lagrangian nodes and elements
        !    that are located in or very near the tile.
        call FORT_FILLCONTAINER( &
          ilev, &
          cache_max_level, &
          cache_max_level, &
          start_time, &
          start_dt, &
          tilelo_array, &
          tilehi_array, &
          xlo_array, &
          dx_current, &
          dx_max_level, &
          num_grids_on_level, &
          num_grids_on_level, &
          gridno_array, &
          num_tiles_on_thread_proc, &
          local_nthreads, &
          max_num_tiles_on_thread_proc, &
          tile_dim, &
          num_materials, &
          global_nparts, &
          im_solid_map, & ! type: 0..nmat-1
          xlo_array, &
          xhi_array)

        iter=0 ! touch_flag=0
        FSI_operation=2 ! make distance in narrow band
        FSI_sub_operation=0
        ! 1.FillCoarsePatch
        ! 2.traverse lagrangian elements belonging to each tile and update
        !   cells within "bounding box" of the element.
        call ns_header_msg_level(ilev, &
          FSI_operation, &
          FSI_sub_operation, &
          iter)

        first_iter=1 
        do while ((FSI_touch_flag.eq.1).or.(first_iter.eq.1)) 
         first_iter=0
   
         FSI_operation=3 ! sign update   
         FSI_sub_operation=0
         call ns_header_msg_level(ilev, &
           FSI_operation, &
           FSI_sub_operation, &
           iter)

         iter=iter+1
        enddo

       enddo ! ilev=0..cache_max_level

      else if (global_nparts.eq.0) then
       ! do nothing
      else
       print *,"global_nparts invalid"
       stop
      endif
 
      return
      end subroutine convert_lag_to_eul

      ! input: UOLD
      ! output: UNEW
      subroutine update_interface(UOLD,UNEW,NCELL,state_ncomp, &
            dx,prev_time,dt,nsteps,nten_in,stefan_flag)
      USE probcommon_module 
      USE probmain_module 
      USE global_utility_module 
      use MOF_routines_module
      use tecplotutil_cpp_module
      use mass_transfer_cpp_module
      use mof_redist_cpp_module
      use plic_cpp_module
      use tsat_module

      IMPLICIT NONE

      INTEGER_T nten_in,stefan_flag
      INTEGER_T NCELL,state_ncomp,nsteps
      INTEGER_T fablo(SDIM)
      INTEGER_T fabhi(SDIM)
      REAL_T dx(SDIM)
      REAL_T prev_time,dt
      REAL_T UOLD(-1:NCELL,-1:NCELL,1:state_ncomp)
      REAL_T UNEW(-1:NCELL,-1:NCELL,1:state_ncomp)
      INTEGER_T nten
      INTEGER_T nmat
      INTEGER_T nLS
      INTEGER_T DIMDEC(maskcov)
      INTEGER_T DIMDEC(masknbr)
      INTEGER_T DIMDEC(nodevel)
      INTEGER_T DIMDEC(deltaVOF)
      INTEGER_T DIMDEC(LS)
      INTEGER_T DIMDEC(LSnew)
      INTEGER_T DIMDEC(LS_slopes_FD)
      INTEGER_T DIMDEC(VOFnew)
      INTEGER_T DIMDEC(Snew)
      INTEGER_T DIMDEC(EOS)
      INTEGER_T DIMDEC(recon)
      INTEGER_T DIMDEC(pres)
      INTEGER_T DIMDEC(FD_NRM_ND)
      INTEGER_T DIMDEC(FD_CURV_CELL)
      INTEGER_T DIMDEC(jump_strength)
      INTEGER_T DIMDEC(swept)
      INTEGER_T DIMDEC(dist_touch)
      INTEGER_T DIMDEC(stencil)
      INTEGER_T DIMDEC(facefrac)
      INTEGER_T DIMDEC(facepairX)
      INTEGER_T DIMDEC(facepairY)
      INTEGER_T DIMDEC(facetest)
      INTEGER_T DIMDEC(slopes)
      REAL_T, dimension(:,:), allocatable :: maskcov
      REAL_T, dimension(:,:,:), allocatable :: masknbr
      REAL_T, dimension(:,:,:), allocatable :: nodevel
      REAL_T, dimension(:,:,:), allocatable :: deltaVOF
      REAL_T, dimension(:,:,:), allocatable :: LS
      REAL_T, dimension(:,:,:), allocatable :: LSnew
      REAL_T, dimension(:,:,:), allocatable :: LS_slopes_FD
      REAL_T, dimension(:,:,:), allocatable :: VOFnew
      REAL_T, dimension(:,:,:), allocatable :: Snew
      REAL_T, dimension(:,:,:), allocatable :: EOS
      REAL_T, dimension(:,:,:), allocatable :: recon
      REAL_T, dimension(:,:), allocatable :: pres
        ! (nmat+nten)*(sdim+1), ngrow_distance
      REAL_T, dimension(:,:,:), allocatable :: FD_NRM_ND
        ! 2*(nmat+nten), ngrow_distance
      REAL_T, dimension(:,:,:), allocatable :: FD_CURV_CELL
        ! 2*nten, ngrow_expansion
      REAL_T, dimension(:,:,:), allocatable :: jump_strength 
        ! ncomp=1, ngrow=0
      REAL_T, dimension(:,:), allocatable :: swept
        ! ncomp=nmat, ngrow=0
      REAL_T, dimension(:,:,:), allocatable :: dist_touch
        ! ncomp=9, ngrow=ngrow_distance
      REAL_T, dimension(:,:,:), allocatable :: stencil
        ! ncomp=nface, ngrow=ngrow_distance
      REAL_T, dimension(:,:,:), allocatable :: facefrac
        ! nmat x nmat x 2 ngrow=ngrow_distance
      REAL_T, dimension(:,:,:), allocatable :: facepairX
      REAL_T, dimension(:,:,:), allocatable :: facepairY
        ! ncomp=nmat*sdim, ngrow=ngrow_distance
      REAL_T, dimension(:,:,:), allocatable :: facetest
        ! ncomp=nmat*ngeom_recon, ngrow=ngrow_distance
      REAL_T, dimension(:,:,:), allocatable :: slopes

      REAL_T DVOF(num_materials)
      REAL_T DVOF_local(num_materials)
      REAL_T delta_mass(2*num_materials)
      REAL_T delta_mass_local(2*num_materials)

      INTEGER_T truncate_volume_fractions(num_materials)

      REAL_T blob_array(2)
      INTEGER_T arraysize
      INTEGER_T radius_cutoff(num_materials)
      INTEGER_T update_flag
      INTEGER_T total_calls(num_materials)
      INTEGER_T total_iterations(num_materials)
      INTEGER_T continuous_mof
      INTEGER_T nstar
      INTEGER_T bfact
      INTEGER_T color_count
      REAL_T cur_time
      INTEGER_T scomp
      INTEGER_T dcomp
      INTEGER_T i,j,k
      INTEGER_T ii,jj
      INTEGER_T dir
      INTEGER_T veldir
      INTEGER_T side
      INTEGER_T do_face_decomp
      INTEGER_T max_level,level,finest_level
      INTEGER_T gridno
      INTEGER_T im
      INTEGER_T im_opp
      INTEGER_T im_primary
      INTEGER_T im_primary_side
      INTEGER_T is_phase_change
      INTEGER_T isweep
      INTEGER_T iten
      INTEGER_T nhalf
      REAL_T max_problen
      REAL_T maxLS(num_materials)
      REAL_T minLS(num_materials)
      REAL_T LSerr,LSexact,VELerr,NRMerr,NRM_FD_err,SPEEDerr
      REAL_T NRM_FD(SDIM)
      REAL_T NRM_FD_mag
      INTEGER_T nden,nface,nface_decomp,npair
      INTEGER_T ngrow,ngrow_distance,ngrow_recon
      INTEGER_T n_normal
      INTEGER_T nprocessed
      INTEGER_T nstate
      INTEGER_T num_elements_blobclass
      INTEGER_T rzflag
      INTEGER_T solvability_projection
      INTEGER_T tessellate
      INTEGER_T tid
      INTEGER_T tid_data
      INTEGER_T vofbc(SDIM,2)
      INTEGER_T velbc(SDIM,2,SDIM)
      REAL_T xsten(-3:3,SDIM)
      REAL_T xgrid(SDIM)
      REAL_T xclosest(SDIM)
      REAL_T rgrid
      REAL_T xlo(SDIM)
      REAL_T problo(SDIM)
      REAL_T probhi(SDIM)
      REAL_T problen(SDIM)
      REAL_T VELexact(SDIM)
      REAL_T NRMexact(SDIM)
      REAL_T SPEEDexact
      REAL_T SPEEDapprox
      REAL_T LScen(num_materials)
      REAL_T LSsten(num_materials)

      INTEGER_T microlayer_substrate(num_materials)
      REAL_T microlayer_angle(num_materials)
      REAL_T microlayer_size(num_materials)
      REAL_T macrolayer_size(num_materials)
      REAL_T max_contact_line_size(num_materials)
      REAL_T microlayer_temperature_substrate(num_materials)
      INTEGER_T distribute_from_target(2*nten_in)
      INTEGER_T mass_fraction_id(2*nten_in)
      REAL_T species_evaporation_density(num_materials)

      REAL_T l1_interface(nten_in)
      REAL_T l2_interface(nten_in)
      REAL_T linf_interface(nten_in)
      REAL_T l1_speed(nten_in)
      REAL_T l2_speed(nten_in)
      REAL_T linf_speed(nten_in)
      REAL_T max_speed_approx(nten_in)
      REAL_T max_speed_exact(nten_in)

      REAL_T l1_nrm(nten_in)
      REAL_T l2_nrm(nten_in)
      REAL_T linf_nrm(nten_in)

      REAL_T l1_nrmFD(nten_in)
      REAL_T l2_nrmFD(nten_in)
      REAL_T linf_nrmFD(nten_in)

      INTEGER_T num_error(nten_in)
      INTEGER_T signchange(nten_in)
      INTEGER_T i_inf(nten_in)
      INTEGER_T j_inf(nten_in)

      INTEGER_T i_inf_nrm(nten_in)
      INTEGER_T j_inf_nrm(nten_in)

      INTEGER_T keep_all_interfaces

      INTEGER_T n_root
      character(len=6) :: root_char_array
      INTEGER_T data_dir,SDC_outer_sweeps,slab_step
      INTEGER_T data_id,visual_revolve,visual_option

      INTEGER_T debug_plot_dir,interior_only,diagnostic_output

      diagnostic_output=0
      nhalf=3

      ncomp_per_burning=SDIM
      ncomp_per_tsat=2

      if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"only SDIM==2 supported here"
       stop
      endif 
      if (NCELL.ge.4) then
       bfact=1
       nmat=num_materials
       nten=((nmat-1)*(nmat-1)+(nmat-1))/2
       if (nten.eq.nten_in) then
        ! do nothing
       else
        print *,"nten_in or nten invalid"
        stop
       endif
       nburning=nten*(ncomp_per_burning+1)
       ntsat=nten*(ncomp_per_tsat+1)
       nstar=9
      else
       print *,"NCELL invalid"
       stop
      endif

      nLS=nmat*(SDIM+1)

      n_normal=(nmat+nten)*(SDIM+1)

      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid"
       stop
      endif
      if (ngeom_raw.eq.1+SDIM) then
       ! do nothing
      else
       print *,"ngeom_raw invalid"
       stop
      endif
      if (nten.ne.global_nten) then
       print *,"nten invalid"
       stop
      endif

      do_face_decomp=0
      nface=nmat*SDIM*2*(1+SDIM)
      nface_decomp=0
      npair=nmat*nmat*2

      problo(1)=problox
      problo(2)=probloy
      probhi(1)=probhix
      probhi(2)=probhiy
      if (SDIM.eq.3) then
       problo(SDIM)=probloz
       probhi(SDIM)=probhiz
      endif

      max_problen=0.0d0
      do dir=1,SDIM
       problen(dir)=probhi(dir)-problo(dir)
       if (problen(dir).gt.0.0d0) then
        max_problen=max_problen+problen(dir)**2
       else
        print *,"problen invalid"
        stop
       endif
       do side=1,2
        if (physbc(dir,side).eq.INT_DIR) then
         vofbc(dir,side)=INT_DIR
         do veldir=1,SDIM
          velbc(dir,side,veldir)=INT_DIR
         enddo
        else if (physbc(dir,side).eq.EXT_DIR) then
         vofbc(dir,side)=FOEXTRAP
         do veldir=1,SDIM
          velbc(dir,side,veldir)=FOEXTRAP
         enddo
        else if (physbc(dir,side).eq.REFLECT_EVEN) then
         vofbc(dir,side)=FOEXTRAP
         do veldir=1,SDIM
          velbc(dir,side,veldir)=FOEXTRAP
         enddo
        else
         print *,"physbc invalid"
         stop
        endif
       enddo ! side=1..2
      enddo ! dir=1..sdim

      max_problen=sqrt(max_problen);

      if (max_problen.gt.0.0d0) then
       ! do nothing
      else
       print *,"max_problen invalid"
       stop
      endif

      nstate=SDIM+1+nmat*(2+ngeom_raw)+1

      if (state_ncomp.eq.nmat+nten*SDIM+ngeom_recon*nmat+nmat*(SDIM+1)) then
       ! do nothing
      else
       print *,"state_ncomp invalid"
       stop
      endif
      is_phase_change=0
      do iten=1,2*nten
       if (latent_heat(iten).ne.0.0d0) then
        is_phase_change=1
       else if (latent_heat(iten).eq.0.0d0) then
        ! do nothing
       else
        print *,"latent_heat bust"
        stop
       endif
      enddo

      if (is_phase_change.eq.1) then
       do dir=1,SDIM
        fablo(dir)=0
        fabhi(dir)=NCELL-1
        xlo(dir)=0.0d0
       enddo

       level=0
       finest_level=0
       ngrow_distance=4
       if (ngrow_make_distance.eq.3) then
        ! do nothing
       else
        print *,"ngrow_make_distance invalid"
        stop
       endif
       ngrow=ngrow_distance
       ngrow_expansion=2

       nden=nmat*2

       do im=1,nmat
        truncate_volume_fractions(im)=0
       enddo

       do im=1,nmat
        radius_cutoff(im)=0
        microlayer_substrate(im)=0
        microlayer_angle(im)=0.0d0
        microlayer_size(im)=0.0d0
        macrolayer_size(im)=0.0d0
        max_contact_line_size(im)=0.0d0
        microlayer_temperature_substrate(im)=0.0d0
       enddo
       do iten=1,2*nten
        distribute_from_target(iten)=0
        mass_fraction_id(iten)=0
       enddo
       do im=1,nmat
        species_evaporation_density(im)=1.0d0
       enddo
       num_elements_blobclass= &
          3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
          2*(2*SDIM)+1+ &
          3+1+2*SDIM+1+nmat+nmat*nmat

       arraysize=num_elements_blobclass
       color_count=0

       call set_dimdec(DIMS(maskcov),fablo,fabhi,ngrow_distance)
       call set_dimdec(DIMS(masknbr),fablo,fabhi,ngrow_distance)
       call set_dimdec(DIMS(burnvel),fablo,fabhi,ngrow_make_distance)
       call set_dimdec(DIMS(tsatfab),fablo,fabhi,ngrow_make_distance)
       call set_dimdec(DIMS(nodevel),fablo,fabhi,1)
       call set_dimdec(DIMS(deltaVOF),fablo,fabhi,0)
       call set_dimdec(DIMS(LS),fablo,fabhi,ngrow)
       call set_dimdec(DIMS(LSnew),fablo,fabhi,ngrow)
       call set_dimdec(DIMS(LS_slopes_FD),fablo,fabhi,ngrow)
       call set_dimdec(DIMS(VOFnew),fablo,fabhi,ngrow)
       call set_dimdec(DIMS(Snew),fablo,fabhi,1)
       call set_dimdec(DIMS(EOS),fablo,fabhi,ngrow)
       call set_dimdec(DIMS(recon),fablo,fabhi,ngrow)
       call set_dimdec(DIMS(pres),fablo,fabhi,ngrow)
       call set_dimdec(DIMS(FD_NRM_ND),fablo,fabhi,ngrow_distance)
       call set_dimdec(DIMS(FD_CURV_CELL),fablo,fabhi,ngrow_distance)
       call set_dimdec(DIMS(jump_strength),fablo,fabhi,ngrow_expansion)
       call set_dimdec(DIMS(swept),fablo,fabhi,0)
       call set_dimdec(DIMS(dist_touch),fablo,fabhi,0)
       call set_dimdec(DIMS(stencil),fablo,fabhi,ngrow_distance)
       call set_dimdec(DIMS(facefrac),fablo,fabhi,ngrow_distance)
       call set_dimdec(DIMS(facepairX),fablo,fabhi,ngrow_distance+1)
       call set_dimdec(DIMS(facepairY),fablo,fabhi,ngrow_distance+1)
       call set_dimdec(DIMS(facetest),fablo,fabhi,ngrow_distance)
       call set_dimdec(DIMS(slopes),fablo,fabhi,ngrow_distance)

       allocate(maskcov(DIMV(maskcov)))
       allocate(masknbr(DIMV(masknbr),4))
       allocate(burnvel(DIMV(burnvel),nburning))
       allocate(tsatfab(DIMV(tsatfab),ntsat))
       allocate(nodevel(DIMV(nodevel),2*nten*SDIM))
       allocate(deltaVOF(DIMV(deltaVOF),nmat))
       allocate(LS(DIMV(LS),nmat*(SDIM+1)))
       allocate(LSnew(DIMV(LSnew),nmat*(SDIM+1)))
       allocate(LS_slopes_FD(DIMV(LS_slopes_FD),nmat*SDIM))
       allocate(VOFnew(DIMV(VOFnew),nmat*ngeom_raw))
       allocate(Snew(DIMV(Snew),nstate))
       allocate(EOS(DIMV(EOS),nden))
       allocate(recon(DIMV(recon),nmat*ngeom_recon)) ! F,X,order,SL,I x nmat
       allocate(pres(DIMV(pres)))
       allocate(FD_NRM_ND(DIMV(FD_NRM_ND),n_normal))
       allocate(FD_CURV_CELL(DIMV(FD_CURV_CELL),2*(nmat+nten)))
       allocate(jump_strength(DIMV(jump_strength),2*nten))
       allocate(swept(DIMV(swept)))
       allocate(dist_touch(DIMV(dist_touch),nmat))
       allocate(stencil(DIMV(stencil),nstar))
       allocate(facefrac(DIMV(facefrac),nface))
       allocate(facepairX(DIMV(facepairX),npair))
       allocate(facepairY(DIMV(facepairY),npair))
       allocate(facetest(DIMV(facetest),nmat*SDIM))
       allocate(slopes(DIMV(slopes),nmat*ngeom_recon))

       do im=1,nmat
        DVOF(im)=0.0d0
        DVOF_local(im)=0.0d0
       enddo
       do im=1,2*nmat
        delta_mass(im)=0.0d0
        delta_mass_local(im)=0.0d0
       enddo

       do i=fablo(1)-ngrow_distance,fabhi(1)+ngrow_distance
       do j=fablo(2)-ngrow_distance,fabhi(2)+ngrow_distance
        maskcov(i,j)=1.0d0
        masknbr(i,j,1)=0.0d0
        masknbr(i,j,2)=0.0d0
        masknbr(i,j,3)=0.0d0
        masknbr(i,j,4)=1.0d0
       enddo
       enddo
       do i=fablo(1),fabhi(1)
       do j=fablo(2),fabhi(2)
        masknbr(i,j,1)=1.0d0
        masknbr(i,j,2)=1.0d0
       enddo
       enddo
       do i=fablo(1)-ngrow_distance+1,fabhi(1)+ngrow_distance-1
       do j=fablo(2)-ngrow_distance+1,fabhi(2)+ngrow_distance-1
        masknbr(i,j,3)=1.0d0
       enddo
       enddo
       do i=fablo(1),fabhi(1)
       do j=fablo(2),fabhi(2)
        swept(i,j)=0.0d0
       enddo
       enddo
       do i=fablo(1),fabhi(1)
       do j=fablo(2),fabhi(2)
       do im=1,nmat
        dist_touch(i,j,im)=0.0d0
       enddo
       enddo
       enddo
       do i=fablo(1)-ngrow_distance-1,fabhi(1)+ngrow_distance+1
       do j=fablo(2)-ngrow_distance-1,fabhi(2)+ngrow_distance+1
        do im=1,npair
         facepairX(i,j,im)=0.0d0
         facepairY(i,j,im)=0.0d0
        enddo
       enddo
       enddo
       do i=fablo(1)-ngrow_distance,fabhi(1)+ngrow_distance
       do j=fablo(2)-ngrow_distance,fabhi(2)+ngrow_distance
        do im=1,nstar
         stencil(i,j,im)=0.0d0
        enddo
        do im=1,nface
         facefrac(i,j,im)=0.0d0
        enddo
        do im=1,nmat*SDIM
         facetest(i,j,im)=0.0d0
        enddo
        do im=1,nmat*ngeom_recon
         slopes(i,j,im)=0.0d0
        enddo
       enddo
       enddo
       do i=fablo(1)-ngrow_make_distance,fabhi(1)+ngrow_make_distance
       do j=fablo(2)-ngrow_make_distance,fabhi(2)+ngrow_make_distance
        do im=1,nburning
         burnvel(i,j,im)=0.0d0
        enddo
        do im=1,ntsat
         tsatfab(i,j,im)=0.0d0
        enddo
        do im=1,n_normal
         FD_NRM_ND(i,j,im)=0.0d0
        enddo
        do im=1,2*(nmat+nten)
         FD_CURV_CELL(i,j,im)=0.0d0
        enddo
       enddo
       enddo
       do i=fablo(1)-1,fabhi(1)+1
       do j=fablo(2)-1,fabhi(2)+1
       do im=1,2*nten*SDIM
        nodevel(i,j,im)=0.0d0
       enddo
       enddo
       enddo
       do i=fablo(1)-ngrow_expansion,fabhi(1)+ngrow_expansion
       do j=fablo(2)-ngrow_expansion,fabhi(2)+ngrow_expansion
       do im=1,2*nten
        jump_strength(i,j,im)=0.0d0
       enddo
       enddo
       enddo
       do i=fablo(1)-ngrow,fabhi(1)+ngrow
       do j=fablo(2)-ngrow,fabhi(2)+ngrow
         pres(i,j)=0.0d0
       enddo
       enddo

       do i=fablo(1),fabhi(1)
       do j=fablo(2),fabhi(2)
        do im=1,nmat*(1+SDIM)
         scomp=nmat+nten*SDIM+ngeom_recon*nmat
         LS(i,j,im)=UOLD(i,j,scomp+im)
         LSnew(i,j,im)=UOLD(i,j,scomp+im)
        enddo
        do im=1,nmat
         do dir=1,ngeom_raw
          dcomp=(im-1)*ngeom_raw+dir
          scomp=nmat+nten*SDIM+(im-1)*ngeom_recon+dir
          VOFnew(i,j,dcomp)=UOLD(i,j,scomp)
         enddo
        enddo
        do im=1,nmat
         scomp=(im-1)*2
         EOS(i,j,scomp+1)=fort_denconst(im)
         EOS(i,j,scomp+2)=UOLD(i,j,im)
        enddo
        do im=1,nmat*ngeom_recon
         scomp=nmat+nten*SDIM
         recon(i,j,im)=UOLD(i,j,scomp+im)
        enddo
       enddo
       enddo
       call set_boundary_VOF( &
        VOFnew,DIMS(VOFnew), &
        fablo,fabhi, &
        nmat,ngrow)
       call set_boundary_recon( &
        recon,DIMS(recon), &
        fablo,fabhi, &
        nmat,ngrow)
       call set_boundary_EOS( &
        EOS,DIMS(EOS), &
        fablo,fabhi, &
        nmat,ngrow)
       call set_boundary_LS( &
        LSnew,DIMS(LSnew), &
        fablo,fabhi, &
        nmat,ngrow)
       call set_boundary_LS( &
        LS,DIMS(LS), &
        fablo,fabhi, &
        nmat,ngrow)

       do i=fablo(1)-ngrow,fabhi(1)+ngrow
       do j=fablo(2)-ngrow,fabhi(2)+ngrow
        do im=1,nmat*SDIM
         LS_slopes_FD(i,j,im)=LS(i,j,nmat+im)
        enddo
       enddo
       enddo

       call FORT_FD_NORMAL( &
         level, &
         finest_level, &
         LS, &
         DIMS(LS), &
         LS_slopes_FD, &
         DIMS(LS_slopes_FD), &
         fablo,fabhi, &
         fablo,fabhi, &
         bfact, &
         xlo,dx, &
         nmat)


       if (DEBUG_LS_MOVE_INTERFACE.eq.1) then
        k=0
        do i=fablo(1)-ngrow,fabhi(1)+ngrow
        do j=fablo(2)-ngrow,fabhi(2)+ngrow
         call gridsten_level(xsten,i,j,k,level,nhalf)
         do dir=1,SDIM
          xgrid(dir)=xsten(0,dir)
         enddo
         do im=1,nmat
          call get_exact_LS(xgrid,cur_time,im,LSexact)
          call get_exact_NRM(xgrid,cur_time,im,NRMexact)
          LS(i,j,im)=LSexact
          LSnew(i,j,im)=LSexact
          do dir=1,SDIM
           LS(i,j,nmat+(im-1)*SDIM+dir)=NRMexact(dir)
           LSnew(i,j,nmat+(im-1)*SDIM+dir)=NRMexact(dir)
           LS_slopes_FD(i,j,(im-1)*SDIM+dir)=NRMexact(dir)
          enddo
         enddo
        enddo
        enddo
       else if (DEBUG_LS_MOVE_INTERFACE.eq.0) then
        ! do nothing
       else
        print *,"DEBUG_LS_MOVE_INTERFACE invalid"
        stop
       endif


       do i=fablo(1)-1,fabhi(1)+1
       do j=fablo(2)-1,fabhi(2)+1
        do im=1,nstate
         Snew(i,j,im)=0.0d0
        enddo
        do im=1,nmat
         dcomp=SDIM+1+(im-1)*2
         scomp=(im-1)*2
         Snew(i,j,dcomp+1)=EOS(i,j,scomp+1) 
         Snew(i,j,dcomp+2)=EOS(i,j,scomp+2) 
         dcomp=SDIM+1+nmat*2+(im-1)*ngeom_raw
         scomp=(im-1)*ngeom_recon
         do dir=1,ngeom_raw
          Snew(i,j,dcomp+dir)=recon(i,j,scomp+dir)
         enddo
        enddo ! im=1..nmat
       enddo
       enddo
     
       ngrow_dest=ngrow_distance-1
 
       call FORT_FD_NODE_NORMAL( &
         level, &
         finest_level, &
         LS,DIMS(LS),  & ! ngrow==ngrow_distance
         FD_NRM_ND,DIMS(FD_NRM_ND),  & ! ngrow==ngrow_distance
         fablo,fabhi, &
         fablo,fabhi,bfact, &
         xlo,dx, &
         nmat, &
         nten, &
         n_normal, &
         nrow_dest)

 
       if (1.eq.1) then 
        ! burnvel flag==1 if valid rate of phase change.
        call FORT_RATEMASSCHANGE( &
         stefan_flag, & ! do not update LSnew if stefan_flag==0
         level, &
         finest_level, &
         normal_probe_size, &
         ngrow_distance, &
         nmat, &
         nten, &
         nburning, &
         ntsat, &
         nden, &
         fort_density_floor, &
         fort_density_ceiling, &
         microlayer_substrate, &
         microlayer_angle, &
         microlayer_size, &
         macrolayer_size, &
         max_contact_line_size, &
         latent_heat, &
         use_exact_temperature, &
         reaction_rate, &
         saturation_temp, &
         freezing_model, &
         distribute_from_target, &
         mass_fraction_id, &
         species_evaporation_density, &
         fablo,fabhi, &
         fablo,fabhi,bfact, &
         xlo,dx, &
         prev_time, &
         dt, &
         arraysize, &
         blob_array, &
         num_elements_blobclass, &
         color_count, &
         maskcov,DIMS(maskcov), & ! colorfab 1 grow unused
         maskcov,DIMS(maskcov), & ! ! typefab 1 grow unused
         maskcov,DIMS(maskcov), & ! 1 grow
         burnvel,DIMS(burnvel), & ! ngrow_make_distance
         tsatfab,DIMS(tsatfab), & ! ngrow_make_distance
         LS,DIMS(LS),  & ! ngrow
         LSnew,DIMS(LSnew), &  ! ngrow
         LS_slopes_FD, &
         DIMS(LS_slopes_FD), &  ! ngrow
         EOS,DIMS(EOS), & ! ngrow
         recon,DIMS(recon), & ! ngrow
         pres,DIMS(pres) ) ! ngrow
       endif


       if (DEBUG_LS_MOVE_INTERFACE.eq.1) then
        k=0
        do i=fablo(1)-ngrow,fabhi(1)+ngrow
        do j=fablo(2)-ngrow,fabhi(2)+ngrow
         call gridsten_level(xsten,i,j,k,level,nhalf)
         do dir=1,SDIM
          xgrid(dir)=xsten(0,dir)
         enddo
         do im=1,nmat
          call get_exact_LS(xgrid,cur_time,im,LSexact)
          call get_exact_NRM(xgrid,cur_time,im,NRMexact)
          LS(i,j,im)=LSexact
          LSnew(i,j,im)=LSexact
          do dir=1,SDIM
           LS(i,j,nmat+(im-1)*SDIM+dir)=NRMexact(dir)
           LSnew(i,j,nmat+(im-1)*SDIM+dir)=NRMexact(dir)
           LS_slopes_FD(i,j,(im-1)*SDIM+dir)=NRMexact(dir)
          enddo
         enddo
        enddo
        enddo
       else if (DEBUG_LS_MOVE_INTERFACE.eq.0) then
        ! do nothing
       else
        print *,"DEBUG_LS_MOVE_INTERFACE invalid"
        stop
       endif

       ngrow=normal_probe_size+3
       if (ngrow.eq.4) then
        ! do nothing
       else
        print *,"ngrow invalid"
        stop
       endif

        ! first nten components are the status
       call set_boundary_burning( &
        burnvel, &
        DIMS(burnvel), &
        fablo,fabhi, &
        nten,nmat,ncomp_per_burning, &
        nburning, &
        ngrow_make_distance)

       call set_boundary_burning( &
        tsatfab, &
        DIMS(tsatfab), &
        fablo,fabhi, &
        nten,nmat,ncomp_per_tsat, &
        ntsat, &
        ngrow_make_distance)

       debug_plot_dir=-1
       interior_only=1

       if (diagnostic_output.eq.1) then
        n_root=6
        root_char_array='burnVL'
        data_dir=-1
        SDC_outer_sweeps=0
        slab_step=0
        data_id=0
        visual_revolve=0
        visual_option=-2

        call FORT_TECPLOTFAB_SANITY( &
         root_char_array, &
         n_root, &
         data_dir, &
         bfact, & 
         fablo,fabhi, &
         burnvel, &
         DIMS(burnvel), &
         problo,probhi, &
         dx, &
         SDC_outer_sweeps, &
         slab_step, &
         data_id, &
         nsteps, &
         prev_time, &  ! cur_time will not show on same mesh as prev_time.
         visual_option, &
         visual_revolve, &
         level, &
         finest_level, &
         nburning)

        root_char_array='LVLSET'
        data_id=1

        call FORT_TECPLOTFAB_SANITY( &
         root_char_array, &
         n_root, &
         data_dir, &
         bfact, & 
         fablo,fabhi, &
         LS, &
         DIMS(LS), &
         problo,probhi, &
         dx, &
         SDC_outer_sweeps, &
         slab_step, &
         data_id, &
         nsteps, &
         prev_time, &  ! cur_time will not show on same mesh as prev_time.
         visual_option, &
         visual_revolve, &
         level, &
         finest_level, &
         nLS)

       endif

       !burnvel is cell centered.
       !burnvel flag set from 0 to 2 if
       !foot of characteristic within range.
       velflag=1
       call FORT_EXTEND_BURNING_VEL( &
         velflag, &
         level, &
         finest_level, &
         xlo,dx, &
         nmat, &
         nten, &
         nburning, &
         ngrow, &
         latent_heat, &
         fablo,fabhi, &
         fablo,fabhi,bfact, &
         burnvel,DIMS(burnvel), & ! ngrow_make_distance
         LS,DIMS(LS))

       velflag=0
       call FORT_EXTEND_BURNING_VEL( &
         velflag, &
         level, &
         finest_level, &
         xlo,dx, &
         nmat, &
         nten, &
         ntsat, &
         ngrow, &
         latent_heat, &
         fablo,fabhi, &
         fablo,fabhi,bfact, &
         tsatfab,DIMS(tsatfab), & ! ngrow_make_distance
         LS,DIMS(LS))

        ! first nten components are the status
       call set_boundary_burning( &
        burnvel, &
        DIMS(burnvel), &
        fablo,fabhi, &
        nten,nmat, &
        ncomp_per_burning,nburning, &
        ngrow_make_distance)

       call set_boundary_burning( &
        tsatfab, &
        DIMS(tsatfab), &
        fablo,fabhi, &
        nten,nmat, &
        ncomp_per_tsat,ntsat, &
        ngrow_make_distance)


       do isweep=0,1

        if (isweep.eq.0) then
         call FORT_NODEDISPLACE( &
          nmat, &
          nten, &
          nburning, &
          fablo,fabhi, &
          fablo,fabhi, &
          bfact, &
          velbc, &
          dt, &
          nodevel, &
          DIMS(nodevel), &
          burnvel, &
          DIMS(burnvel), &
          xlo,dx, &
          level,finest_level)
        else if (isweep.eq.1) then
         ! do nothing
        else
         print *,"isweep invalid"
         stop
        endif

        if (isweep.eq.0) then
                ! do nothing
        else if (isweep.eq.1) then
         do im=1,nmat
          DVOF_local(im)=DVOF(im)
         enddo
        else
         print *,"isweep invalid"
         stop
        endif

        tid=0
        tid_data=0
        solvability_projection=0

        if (stefan_flag.eq.1) then
         ! do nothing
        else if (stefan_flag.eq.0) then
         do i=fablo(1)-1,fabhi(1)+1
         do j=fablo(2)-1,fabhi(2)+1
         do im=1,2*nten*SDIM
          nodevel(i,j,im)=0.0d0
         enddo
         enddo
         enddo
        else
         print *,"stefan_flag invalid"
         stop
        endif

        call FORT_CONVERTMATERIAL( &
         tid_data, &
         isweep, &
         solvability_projection, &
         ngrow_expansion, &
         level,finest_level, &
         normal_probe_size, &
         nmat, &
         nten, &
         nden, &
         nstate, &
         ntsat, &
         fort_density_floor, &
         fort_density_ceiling, &
         latent_heat, &
         saturation_temp, &
         freezing_model, &
         mass_fraction_id, &
         species_evaporation_density, &
         distribute_from_target, &
         fablo,fabhi, &
         fablo,fabhi, &
         bfact,  &
         vofbc, &
         xlo,dx, &
         dt, &
         delta_mass_local, &
         DVOF_local, &
         maskcov,DIMS(maskcov), &
         deltaVOF,DIMS(deltaVOF), &
         nodevel,DIMS(nodevel), &
         jump_strength, &
         DIMS(jump_strength), &
         tsatfab,DIMS(tsatfab), &
         LS,DIMS(LS), &
         LSnew,DIMS(LSnew), &
         recon,DIMS(recon), &
         Snew,DIMS(Snew), &
         EOS,DIMS(EOS), &
         swept,DIMS(swept))

        if (isweep.eq.0) then
         do im=1,nmat
          DVOF(im)=DVOF(im)+DVOF_local(im)
         enddo
        else if (isweep.eq.1) then
         do im=1,2*nmat
          delta_mass(im)=delta_mass(im)+delta_mass_local(im)
         enddo
        else
         print *,"isweep invalid"
         stop
        endif

       enddo ! isweep=0..1


       if (DEBUG_LS_MOVE_INTERFACE.eq.1) then
        k=0
        do i=fablo(1)-ngrow,fabhi(1)+ngrow
        do j=fablo(2)-ngrow,fabhi(2)+ngrow
         call gridsten_level(xsten,i,j,k,level,nhalf)
         do dir=1,SDIM
          xgrid(dir)=xsten(0,dir)
         enddo
         do im=1,nmat
          call get_exact_LS(xgrid,cur_time,im,LSexact)
          call get_exact_NRM(xgrid,cur_time,im,NRMexact)
          LS(i,j,im)=LSexact
          LSnew(i,j,im)=LSexact
          do dir=1,SDIM
           LS(i,j,nmat+(im-1)*SDIM+dir)=NRMexact(dir)
           LSnew(i,j,nmat+(im-1)*SDIM+dir)=NRMexact(dir)
           LS_slopes_FD(i,j,(im-1)*SDIM+dir)=NRMexact(dir)
          enddo
         enddo
        enddo
        enddo
       else if (DEBUG_LS_MOVE_INTERFACE.eq.0) then
        ! do nothing
       else
        print *,"DEBUG_LS_MOVE_INTERFACE invalid"
        stop
       endif


       do i=fablo(1),fabhi(1)
       do j=fablo(2),fabhi(2)
        do im=1,nmat*(SDIM+1)
         LS(i,j,im)=LSnew(i,j,im)
        enddo
        do im=1,nmat*ngeom_raw
         scomp=SDIM+1+nmat*2+im
         VOFnew(i,j,im)=Snew(i,j,scomp)
        enddo
       enddo
       enddo
       call set_boundary_VOF( &
        VOFnew,DIMS(VOFnew), &
        fablo,fabhi, &
        nmat,ngrow)
       call set_boundary_LS( &
        LSnew,DIMS(LSnew), &
        fablo,fabhi, &
        nmat,ngrow)
       call set_boundary_LS( &
        LS,DIMS(LS), &
        fablo,fabhi, &
        nmat,ngrow)

       tid=0
       gridno=0
       level=0
       finest_level=0
       max_level=0
       ngrow_recon=1
       update_flag=0
       total_calls=0
       total_iterations=0
       continuous_mof=0

       print *,"first sloperecon"

       call FORT_SLOPERECON( &
        tid, &
        gridno, &
        level, &
        finest_level, &
        max_level, &
        ngrow_recon, &
        fablo,fabhi, &
        fablo,fabhi,bfact, &
        xlo,dx, &
        masknbr,DIMS(masknbr), &
        Snew,DIMS(Snew), &
        VOFnew,DIMS(VOFnew), & ! nmat x ngeom_raw
        LS,DIMS(LS), &  ! nmat
        slopes,DIMS(slopes), & ! nmat x ngeom_recon
        nsteps, &
        prev_time, &
        nmat,nten, &
        latent_heat, &
        update_flag, &
        total_calls, &
        total_iterations, &
        continuous_mof, &
        radius_cutoff)

       call set_boundary_recon( &
        slopes,DIMS(slopes), &
        fablo,fabhi, &
        nmat,ngrow)

       do im=1,nmat
        minLS(im)=max_problen
        maxLS(im)=-max_problen
       enddo

       rzflag=0

       tessellate=0
       cur_time=prev_time+dt

       call FORT_FACEINIT( &
        tid, &
        tessellate, &
        nten, &
        level, &
        finest_level, &
        facefrac, &
        DIMS(facefrac), &
        masknbr, &
        DIMS(masknbr), &
        slopes, &
        DIMS(slopes), &
        fablo,fabhi, &
        fablo,fabhi, &
        bfact, &
        rzflag, &
        xlo,dx, &
        cur_time, &
        ngrow_distance, &
        nmat, &
        nface, &
        nface_decomp)

       dir=0
       call FORT_FACEPROCESS( &
        ngrow_distance, &
        ngrow_distance, &
        tid, &
        dir, &
        tessellate, &
        level, &
        finest_level, &
        facepairX,DIMS(facepairX), &
        facefrac,DIMS(facefrac), &
        slopes,DIMS(slopes), &
        fablo,fabhi, &
        fablo,fabhi,bfact, &
        rzflag, &
        xlo,dx, &
        cur_time, &
        nmat,nface,npair)


       dir=1
       call FORT_FACEPROCESS( &
        ngrow_distance, &
        ngrow_distance, &
        tid, &
        dir, &
        tessellate, &
        level, &
        finest_level, &
        facepairY,DIMS(facepairY), &
        facefrac,DIMS(facefrac), &
        slopes,DIMS(slopes), &
        fablo,fabhi, &
        fablo,fabhi,bfact, &
        rzflag, &
        xlo,dx, &
        cur_time, &
        nmat,nface,npair)

       call FORT_FACEINITTEST(  &
        tid, &
        tessellate, &
        level, &
        finest_level, &
        facefrac, &
        DIMS(facefrac), &
        facetest, &
        DIMS(facetest), &
        masknbr, &
        DIMS(masknbr), &
        slopes, &
        DIMS(slopes), &
        fablo,fabhi, &
        fablo,fabhi,bfact, &
        rzflag, &
        xlo,dx, &
        cur_time, &
        ngrow_distance, &
        nmat, &
        nface)

       call FORT_STENINIT( & 
        level, &
        finest_level, &
        stencil, &
        DIMS(stencil), &
        masknbr, &
        DIMS(masknbr), &
        slopes, &
        DIMS(slopes), &
        fablo,fabhi, &
        fablo,fabhi,bfact, &
        rzflag, &
        xlo,dx, &
        cur_time, &
        ngrow_distance, &
        nmat,nstar)

       nprocessed=0

       keep_all_interfaces=1

       call FORT_LEVELSTRIP(  &
        keep_all_interfaces, &
        nprocessed, &
        minLS, &
        maxLS, &
        max_problen, &
        level, &
        finest_level, &
        truncate_volume_fractions, &
        latent_heat, &
        masknbr, &
        DIMS(masknbr), &
        facepairX, &
        DIMS(facepairX), &
        facepairY, &
        DIMS(facepairY), &
        facepairY, &
        DIMS(facepairY), &
        facefrac, &
        DIMS(facefrac), &
        facetest, &
        DIMS(facetest), &
        stencil, &
        DIMS(stencil), &
        slopes, &
        DIMS(slopes), &
        LS, &
        DIMS(LS), &
        LSnew, &
        DIMS(LSnew), &
        dist_touch, &
        DIMS(dist_touch), &
        dist_touch, &
        DIMS(dist_touch), &
        LSnew, &
        DIMS(LSnew), &
        fablo,fabhi, &
        fablo,fabhi,bfact, &
        vofbc, &
        rzflag, &
        xlo,dx, &
        cur_time, &
        ngrow_distance, &
        nmat,nten,nstar,nface,npair)

       call FORT_CORRECT_UNINIT(  &
        minLS, &
        maxLS, &
        max_problen, &
        level, &
        finest_level, &
        LSnew, &
        DIMS(LSnew), &
        dist_touch, &
        DIMS(dist_touch), &
        fablo,fabhi, &
        fablo,fabhi,bfact, &
        xlo,dx, &
        cur_time, &
        nmat)

       call set_boundary_LS( &
        LSnew,DIMS(LSnew), &
        fablo,fabhi, &
        nmat,ngrow)

       tid=0
       gridno=0
       level=0
       finest_level=0
       max_level=0
       ngrow_recon=1
       update_flag=0
       total_calls=0
       total_iterations=0
       continuous_mof=0

       print *,"second sloperecon"

       call FORT_SLOPERECON( &
        tid, &
        gridno, &
        level, &
        finest_level, &
        max_level, &
        ngrow_recon, &
        fablo,fabhi, &
        fablo,fabhi,bfact, &
        xlo,dx, &
        masknbr,DIMS(masknbr), &
        Snew,DIMS(Snew), &
        VOFnew,DIMS(VOFnew), & ! nmat x ngeom_raw
        LSnew,DIMS(LSnew), &  ! nmat
        slopes,DIMS(slopes), & ! nmat x ngeom_recon
        nsteps, &
        prev_time, &
        nmat,nten, &
        latent_heat, &
        update_flag, &
        total_calls, &
        total_iterations, &
        continuous_mof, &
        radius_cutoff)

       call set_boundary_recon( &
        slopes,DIMS(slopes), &
        fablo,fabhi, &
        nmat,ngrow)

       if (DEBUG_LS_MOVE_INTERFACE.eq.1) then
        k=0
        do i=fablo(1)-ngrow,fabhi(1)+ngrow
        do j=fablo(2)-ngrow,fabhi(2)+ngrow
         call gridsten_level(xsten,i,j,k,level,nhalf)
         do dir=1,SDIM
          xgrid(dir)=xsten(0,dir)
         enddo
         do im=1,nmat
          call get_exact_LS(xgrid,cur_time,im,LSexact)
          call get_exact_NRM(xgrid,cur_time,im,NRMexact)
          LS(i,j,im)=LSexact
          LSnew(i,j,im)=LSexact
          do dir=1,SDIM
           LS(i,j,nmat+(im-1)*SDIM+dir)=NRMexact(dir)
           LSnew(i,j,nmat+(im-1)*SDIM+dir)=NRMexact(dir)
           LS_slopes_FD(i,j,(im-1)*SDIM+dir)=NRMexact(dir)
          enddo
         enddo
        enddo
        enddo
       else if (DEBUG_LS_MOVE_INTERFACE.eq.0) then
        ! do nothing
       else
        print *,"DEBUG_LS_MOVE_INTERFACE invalid"
        stop
       endif

       do i=fablo(1)-1,fabhi(1)+1
       do j=fablo(2)-1,fabhi(2)+1
        do im=1,nmat*(1+SDIM)
         dcomp=nmat+nten*SDIM+ngeom_recon*nmat
         UNEW(i,j,dcomp+im)=LSnew(i,j,im)
        enddo
        do im=1,nmat
         do dir=1,ngeom_raw
          scomp=(im-1)*ngeom_raw+dir
          dcomp=nmat+nten*SDIM+(im-1)*ngeom_recon+dir
          UNEW(i,j,dcomp)=VOFnew(i,j,scomp)
         enddo
         do dir=ngeom_raw+1,ngeom_recon
          scomp=(im-1)*ngeom_recon+dir
          dcomp=nmat+nten*SDIM+(im-1)*ngeom_recon+dir
          UNEW(i,j,dcomp)=slopes(i,j,scomp)
         enddo
        enddo ! im=1..nmat
         ! temperature
        do im=1,nmat
         dcomp=im
         scomp=SDIM+1+(im-1)*2+2
         UNEW(i,j,dcomp)=Snew(i,j,scomp)
        enddo
         ! phase change velocity
        do im=1,nten*SDIM
         dcomp=nmat+im
         scomp=nten+im
         UNEW(i,j,dcomp)=burnvel(i,j,scomp)
        enddo
       enddo
       enddo

       do iten=1,nten
        l1_interface(iten)=0.0d0
        l2_interface(iten)=0.0d0
        linf_interface(iten)=0.0d0
        i_inf(iten)=fablo(1)-1
        j_inf(iten)=fablo(2)-1

        l1_nrm(iten)=0.0d0
        l2_nrm(iten)=0.0d0
        linf_nrm(iten)=0.0d0

        i_inf_nrm(iten)=fablo(1)-1
        j_inf_nrm(iten)=fablo(2)-1

        l1_nrmFD(iten)=0.0d0
        l2_nrmFD(iten)=0.0d0
        linf_nrmFD(iten)=0.0d0

        l1_speed(iten)=0.0d0
        l2_speed(iten)=0.0d0
        linf_speed(iten)=0.0d0
        max_speed_approx(iten)=0.0d0
        max_speed_exact(iten)=0.0d0
        num_error(iten)=0
       enddo

       k=0
       do i=fablo(1),fabhi(1)
       do j=fablo(2),fabhi(2)
        call gridsten_level(xsten,i,j,k,level,nhalf)
        do dir=1,SDIM
         xgrid(dir)=xsten(0,dir)
        enddo
        do im=1,nmat
         LScen(im)=LSnew(i,j,im)
        enddo
        call get_primary_material(LScen,nmat,im_primary)
        if ((im_primary.ge.1).and.(im_primary.le.nmat)) then
         do iten=1,nten
          signchange(iten)=0
         enddo
         do ii=-1,1
         do jj=-1,1
          if (abs(ii)+abs(jj).eq.1) then
           do im=1,nmat
            LSsten(im)=LSnew(i+ii,j+jj,im)
           enddo
           call get_primary_material(LSsten,nmat,im_primary_side)
           if (im_primary_side.ne.im_primary) then
            call get_iten(im_primary,im_primary_side,iten,nmat)
            signchange(iten)=1
           endif
          else if (abs(ii)+abs(jj).eq.2) then
           ! do nothing
          else if (abs(ii)+abs(jj).eq.0) then
           ! do nothing
          else
           print *,"ii or jj invalid"
           stop
          endif 
         enddo 
         enddo 
         do iten=1,nten
          if (signchange(iten).eq.1) then
           call get_exact_LS(xgrid,cur_time,im_primary,LSexact)
           call get_exact_NRM(xgrid,cur_time,im_primary,NRMexact)
           call get_exact_VEL(xgrid,dx,cur_time,iten,VELexact)
           SPEEDexact=0.0d0
           SPEEDapprox=0.0d0
           do dir=1,SDIM
            scomp=nten+(iten-1)*SDIM+dir
            SPEEDexact=SPEEDexact+VELexact(dir)**2
            SPEEDapprox=SPEEDapprox+burnvel(i,j,scomp)**2
           enddo
           SPEEDexact=sqrt(SPEEDexact)
           SPEEDapprox=sqrt(SPEEDapprox)

           rgrid=sqrt((xgrid(1)-xblob)**2+(xgrid(2)-yblob)**2)
           if (1.eq.0) then
            print *,"rgrid,speed_exact,speed_approx ",rgrid, &
             SPEEDexact,SPEEDapprox
           endif

           num_error(iten)=num_error(iten)+1
           LSerr=abs(LSexact-LScen(im_primary))

           VELerr=0.0d0
           do dir=1,SDIM
            scomp=nten+(iten-1)*SDIM+dir
            VELerr=VELerr+(VELexact(dir)-burnvel(i,j,scomp))**2
           enddo
           VELerr=sqrt(VELerr)

           NRMerr=0.0d0
           NRM_FD_err=0.0d0
           do dir=1,SDIM
            NRM_FD(dir)=LS_slopes_FD(i,j,(im_primary-1)*SDIM+dir)
           enddo
           NRM_FD_mag=zero
           do dir=1,SDIM
            NRM_FD_mag=NRM_FD_mag+NRM_FD(dir)**2
           enddo
           NRM_FD_mag=sqrt(NRM_FD_mag)
           if (NRM_FD_mag.gt.zero) then
            do dir=1,SDIM
             NRM_FD(dir)=NRM_FD(dir)/NRM_FD_mag
            enddo
           endif
           do dir=1,SDIM
            scomp=nmat+(im_primary-1)*SDIM+dir
            NRMerr=NRMerr+(NRMexact(dir)-LSnew(i,j,scomp))**2
            NRM_FD_err=NRM_FD_err+(NRMexact(dir)-NRM_FD(dir))**2
           enddo
           NRMerr=sqrt(NRMerr)
           NRM_FD_err=sqrt(NRM_FD_err)

           SPEEDerr=abs(SPEEDexact-SPEEDapprox)

           if (SPEEDexact.gt.max_speed_exact(iten)) then
            max_speed_exact(iten)=SPEEDexact
           endif
           if (SPEEDapprox.gt.max_speed_approx(iten)) then
            max_speed_approx(iten)=SPEEDapprox
           endif
           l1_interface(iten)=l1_interface(iten)+LSerr
           l2_interface(iten)=l2_interface(iten)+LSerr**2
           if (LSerr.gt.linf_interface(iten)) then
            linf_interface(iten)=LSerr
           endif

           l1_nrm(iten)=l1_nrm(iten)+NRMerr
           l2_nrm(iten)=l2_nrm(iten)+NRMerr**2
           if (NRMerr.gt.linf_nrm(iten)) then
            linf_nrm(iten)=NRMerr
            i_inf_nrm(iten)=i
            j_inf_nrm(iten)=j
           endif

           l1_nrmFD(iten)=l1_nrmFD(iten)+NRM_FD_err
           l2_nrmFD(iten)=l2_nrmFD(iten)+NRM_FD_err**2
           if (NRM_FD_err.gt.linf_nrmFD(iten)) then
            linf_nrmFD(iten)=NRM_FD_err
           endif

           l1_speed(iten)=l1_speed(iten)+SPEEDerr
           l2_speed(iten)=l2_speed(iten)+SPEEDerr**2
           if (SPEEDerr.gt.linf_speed(iten)) then
            linf_speed(iten)=SPEEDerr
            i_inf(iten)=i
            j_inf(iten)=j
           endif
          else if (signchange(iten).eq.0) then
           ! do nothing
          else
           print *,"signchange invalid"
           stop
          endif
         enddo ! iten=1..nten
        else
         print *,"im_primary invalid"
         stop
        endif
       enddo
       enddo

       do iten=1,nten
        if (num_error(iten).gt.0) then
         print *,"TIME=",cur_time," iten=",iten," num_error=", &
                 num_error(iten)
         print *,"TIME=",cur_time," iten=",iten," l1_int=", &
                 l1_interface(iten)/num_error(iten)
         print *,"TIME=",cur_time," iten=",iten," l2_int=", &
                 sqrt(l2_interface(iten)/num_error(iten))
         print *,"TIME=",cur_time," iten=",iten," linf_int=", &
                 linf_interface(iten)

         print *,"TIME=",cur_time," iten=",iten," l1_nrm=", &
                 l1_nrm(iten)/num_error(iten)
         print *,"TIME=",cur_time," iten=",iten," l2_nrm=", &
                 sqrt(l2_nrm(iten)/num_error(iten))
         print *,"TIME=",cur_time," iten=",iten," linf_nrm=", &
                 linf_nrm(iten)

         print *,"TIME=",cur_time," iten=",iten," l1_nrmFD=", &
                 l1_nrmFD(iten)/num_error(iten)
         print *,"TIME=",cur_time," iten=",iten," l2_nrmFD=", &
                 sqrt(l2_nrmFD(iten)/num_error(iten))
         print *,"TIME=",cur_time," iten=",iten," linf_nrmFD=", &
                 linf_nrmFD(iten)

         print *,"TIME=",cur_time," iten=",iten," l1_speed=", &
                 l1_speed(iten)/num_error(iten)
         print *,"TIME=",cur_time," iten=",iten," l2_speed=", &
                 sqrt(l2_speed(iten)/num_error(iten))
         print *,"TIME=",cur_time," iten=",iten," linf_speed=", &
                 linf_speed(iten)
         print *,"TIME=",cur_time," iten=",iten," max_speed_approx=", &
                 max_speed_approx(iten)
         print *,"TIME=",cur_time," iten=",iten," max_speed_exact=", &
                 max_speed_exact(iten)

          ! on 256 x 256 grid, cell (44,67) normal has an error because
          ! the interface in cell (45,68) was improperly ignored.
          ! also, the facetest criteria should be reconsidered.
          ! x,y =  0.17382812500000000       0.26367187500000000
          ! x,y (closest)=  0.17578125000000000       0.26562500000000000
         if (i_inf_nrm(iten).ge.fablo(1)) then
          i=i_inf_nrm(iten)
          j=j_inf_nrm(iten)
          print *,"TIME=",cur_time," iten=",iten," i,j (normal inf)=",i,j
          call get_inverse_iten(im,im_opp,iten,nmat)
          call gridsten_level(xsten,i,j,k,level,nhalf)
          do dir=1,SDIM
           xgrid(dir)=xsten(0,dir)
           scomp=nmat+(im-1)*SDIM+dir
           xclosest(dir)=xgrid(dir)-LSnew(i,j,im)*LSnew(i,j,scomp)
          enddo
          print *,"TIME=",cur_time," iten=",iten," x,y (normal inf)=", &
             xgrid(1),xgrid(2)
          print *,"TIME=",cur_time," iten=",iten," x,y (closest)=", &
             xclosest(1),xclosest(2)

          print *,"TIME=",cur_time," iten=",iten," im,im_opp (inf)=", &
            im,im_opp
          do ii=-1,1
          do jj=-1,1
           call gridsten_level(xsten,i+ii,j+jj,k,level,nhalf)
           do dir=1,SDIM
            xgrid(dir)=xsten(0,dir)
           enddo
           call get_exact_VEL(xgrid,dx,cur_time,iten,VELexact)
           print *,"ii,jj,LS(im),LS(im_opp) ",ii,jj, &
              LSnew(i+ii,j+jj,im),LSnew(i+ii,j+jj,im_opp)
           dcomp=nmat+nten*SDIM+(im-1)*ngeom_recon+1
           print *,"ii,jj,vof(im) ",ii,jj, &
              UNEW(i+ii,j+jj,dcomp)
           dcomp=nmat+nten*SDIM+(im_opp-1)*ngeom_recon+1
           print *,"ii,jj,vof(im_opp) ",ii,jj, &
              UNEW(i+ii,j+jj,dcomp)
           scomp=nten+(iten-1)*SDIM
           print *,"ii,jj,u,v ",ii,jj, &
              burnvel(i+ii,j+jj,scomp+1), &
              burnvel(i+ii,j+jj,scomp+2)
           print *,"ii,jj,uexact,vexact ",ii,jj, &
              VELexact(1),VELexact(2)
          enddo ! jj
          enddo ! ii

         endif

         if (i_inf(iten).ge.fablo(1)) then
          i=i_inf(iten)
          j=j_inf(iten)
          k=0
          call gridsten_level(xsten,i,j,k,level,nhalf)
          print *,"TIME=",cur_time," iten=",iten," i,j (inf_spd)=",i,j, &
           " x,y (inf_spd)= ",xsten(0,1),xsten(0,2)

          call get_inverse_iten(im,im_opp,iten,nmat)
          print *,"TIME=",cur_time," iten=",iten," im,im_opp (inf)=", &
            im,im_opp
          do ii=-1,1
          do jj=-1,1
           call gridsten_level(xsten,i+ii,j+jj,k,level,nhalf)
           do dir=1,SDIM
            xgrid(dir)=xsten(0,dir)
           enddo
           call get_exact_VEL(xgrid,dx,cur_time,iten,VELexact)
           print *,"ii,jj,LS(im),LS(im_opp) ",ii,jj, &
              LSnew(i+ii,j+jj,im),LSnew(i+ii,j+jj,im_opp)
           dcomp=nmat+nten*SDIM+(im-1)*ngeom_recon+1
           print *,"ii,jj,vof(im) ",ii,jj, &
              UNEW(i+ii,j+jj,dcomp)
           dcomp=nmat+nten*SDIM+(im_opp-1)*ngeom_recon+1
           print *,"ii,jj,vof(im_opp) ",ii,jj, &
              UNEW(i+ii,j+jj,dcomp)
           scomp=nten+(iten-1)*SDIM
           print *,"ii,jj,u,v ",ii,jj, &
              burnvel(i+ii,j+jj,scomp+1), &
              burnvel(i+ii,j+jj,scomp+2)
           print *,"ii,jj,uexact,vexact ",ii,jj, &
              VELexact(1),VELexact(2)
          enddo ! jj
          enddo ! ii
         endif
        else if (num_error(iten).eq.0) then
         ! do nothing
        else
         print *,"num_error invalid"
         stop
        endif
       enddo ! iten=1..nten
 
       deallocate(maskcov)
       deallocate(masknbr)
       deallocate(burnvel)
       deallocate(tsatfab)
       deallocate(nodevel)
       deallocate(deltaVOF)
       deallocate(LS)
       deallocate(LSnew)
       deallocate(LS_slopes_FD)
       deallocate(VOFnew)
       deallocate(Snew)
       deallocate(EOS)
       deallocate(recon)
       deallocate(pres)
       deallocate(FD_NRM_ND)
       deallocate(FD_CURV_CELL)
       deallocate(jump_strength)
       deallocate(swept)
       deallocate(dist_touch)
       deallocate(stencil)
       deallocate(facefrac)
       deallocate(facepairX)
       deallocate(facepairY)
       deallocate(facetest)
       deallocate(slopes)

      else if (is_phase_change.eq.0) then
       ! do nothing
      else
       print *,"is_phase_change invalid"
       stop
      endif

      return
      end subroutine update_interface
