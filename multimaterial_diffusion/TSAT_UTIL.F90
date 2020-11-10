#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
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
      INTEGER_T :: use_tsatfab
      INTEGER_T :: i_tsat
      INTEGER_T :: j_tsat
      INTEGER_T :: k_tsat
      INTEGER_T :: ireverse_tsat
      INTEGER_T :: iten_tsat
      INTEGER_T :: bfact_tsat ! in init_tsatfab
      INTEGER_T :: level_tsat ! in init_tsatfab
      INTEGER_T :: finest_level_tsat ! in init_tsatfab
      REAL_T :: dx_tsat(SDIM)  ! in init_tsatfab
      REAL_T :: xlo_tsat(SDIM) ! in init_tsatfab
      INTEGER_T :: ngrow_tsat  ! in init_tsatfab
      INTEGER_T :: fablo_tsat(SDIM) ! in init_tsatfab
      INTEGER_T :: fabhi_tsat(SDIM) ! in init_tsatfab

        ! DIMDEC is a macro
        ! burnvelxlo,burnvelxhi,burnvelylo,burnvelyhi in 2d
        ! burnvelxlo,burnvelxhi,burnvelylo,burnvelyhi,
        ! burnvelzlo,burnvelzhi   in 3d
      INTEGER_T DIMDEC(burnvel)
      INTEGER_T DIMDEC(tsatfab)
      INTEGER_T DIMDEC(swept)

      REAL_T, dimension(:,:,:), allocatable :: burnvel
      REAL_T, dimension(:,:,:), allocatable :: tsatfab
        ! ncomp=nmat, ngrow=0
      REAL_T, dimension(:,:,:), allocatable :: swept

contains

      subroutine delete_tsatfab()

      deallocate(tsatfab)
      deallocate(swept)

      return
      end subroutine delete_tsatfab

      subroutine init_tsatfab(NCELL)
      USE probcommon_module 
      use global_utility_module, only: set_dimdec
      IMPLICIT NONE
      INTEGER_T, intent(in) :: NCELL
      INTEGER_T nten,nmat,dir
      INTEGER_T i,j,im

      level_tsat=0
      finest_level_tsat=0
      bfact_tsat=1

      nmat=num_materials
      nten=((nmat-1)*(nmat-1)+(nmat-1))/2

      ncomp_per_tsat=2
      ntsat=nten*(ncomp_per_tsat+1)

      if (nmat.ge.2) then
       ! do nothing
      else
       print *,"require nmat>=2"
       stop
      endif
      if (nten.ge.1) then
       ! do nothing
      else
       print *,"require nten>=1"
       stop
      endif
      if (ntsat.ge.3) then
       ! do nothing
      else
       print *,"ntsat invalid"
       stop
      endif
      if (NCELL.ge.4) then
       ! do nothing
      else
       print *,"need NCELL>=4"
       stop
      endif

      use_tsatfab=0

      do dir=1,SDIM
       fablo_tsat(dir)=0
       fabhi_tsat(dir)=NCELL-1
       xlo_tsat(dir)=0.0d0
      enddo
      dx_tsat(1)=probhix/NCELL
      dx_tsat(2)=probhiy/NCELL
      if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"only 2d supported right now"
       stop
      endif
      ngrow_tsat=ngrow_make_distance
      call set_dimdec(DIMS(tsatfab),fablo_tsat,fabhi_tsat, &
        ngrow_tsat)
      allocate(tsatfab(DIMV(tsatfab),ntsat))
      call set_dimdec(DIMS(swept),fablo_tsat,fabhi_tsat,0)
      allocate(swept(DIMV(swept),nmat))
      do i=0,NCELL-1
      do j=0,NCELL-1
       do im=1,ntsat
        tsatfab(i,j,im)=0.0d0
       enddo
       do im=1,nmat
        swept(i,j,im)=1.0d0
       enddo
      enddo
      enddo

      return
      end subroutine init_tsatfab

      end module tsat_module
