#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_FORT_INTEGER.H>
#include "AMReX_SPACE.H"
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include "AMReX_ArrayLim.H"
#include "CG_F.H"

      module cpp_cg

      contains

!-----------------------------------------------------------------------
!      
!     CGUPDATE: Modify the input arrays as follows:
!     
!     phi = yy + a*pp
!     
!-----------------------------------------------------------------------
      subroutine fort_cgupdate( &
       phi,  &
       DIMS(phi), &
       a, &
       yy,   &
       DIMS(yy), &
       pp,   &
       DIMS(pp), &
       tilelo, tilehi, &
       fablo,fabhi,bfact,bfact_top) &
      bind(c,name='fort_cgupdate')

      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3)
      INTEGER_T :: growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(pp)
      REAL_T, intent(in), target :: pp(DIMV(pp))
      INTEGER_T, intent(in) :: DIMDEC(phi)
      REAL_T, intent(out), target :: phi(DIMV(phi))
      REAL_T, pointer :: phi_ptr(D_DECL(:,:,:))

      INTEGER_T, intent(in) :: DIMDEC(yy)
      REAL_T, intent(in), target :: yy(DIMV(yy))
      REAL_T, intent(in) :: a

      INTEGER_T i,j,k

      phi_ptr=>phi

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif

      call checkbound_array1(fablo,fabhi,phi_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,pp,0,-1)
      call checkbound_array1(fablo,fabhi,yy,0,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       phi(D_DECL(i,j,k)) = yy(D_DECL(i,j,k)) +  &
         a * pp(D_DECL(i,j,k))
      end do
      end do
      end do

      end subroutine fort_cgupdate
!-----------------------------------------------------------------------
!      
!     CGADVCP: Modify the input arrays as follows:
!     
!     pp = rr + b*yy    through range lo:hi
!     
!     pp  <=>
!     rr  <=>
!     b   <=
!     
!-----------------------------------------------------------------------
      subroutine fort_cgadvcp( &
       pp,  &
       DIMS(pp), &
       rr,  &
       DIMS(rr), &
       yy,  &
       DIMS(yy), &
       b, &
       tilelo, tilehi, &
       fablo,fabhi,bfact,bfact_top ) &
      bind(c,name='fort_cgadvcp')

      use global_utility_module
      IMPLICIT NONE
!
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3)
      INTEGER_T :: growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(rr)
      REAL_T, intent(in), target :: rr(DIMV(rr))
      INTEGER_T, intent(in) :: DIMDEC(pp)
      REAL_T, intent(out), target :: pp(DIMV(pp))
      REAL_T, pointer :: pp_ptr(D_DECL(:,:,:))

      INTEGER_T, intent(in) :: DIMDEC(yy)
      REAL_T, intent(in), target :: yy(DIMV(yy))
      REAL_T, intent(in) :: b

      INTEGER_T i,j,k

      pp_ptr=>pp

      call checkbound_array1(fablo,fabhi,pp_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,rr,0,-1)
      call checkbound_array1(fablo,fabhi,yy,0,-1)

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       pp(D_DECL(i,j,k)) = rr(D_DECL(i,j,k)) + b*yy(D_DECL(i,j,k))
      end do
      end do
      end do
!
      end subroutine fort_cgadvcp
      
!-----------------------------------------------------------------------
!      
!     CGAXP: Modify the input as follows:
!     
!     pw = Transpose(pp) . ww    through range lo:hi
!     
!     pw   =>
!     pp  <=
!     ww  <=
!     
!-----------------------------------------------------------------------
      subroutine fort_cgxdoty( &
       ncomp, &
       pw, &
       pp,  &
       DIMS(pp), &
       ww,  &
       DIMS(ww), &
       tilelo, tilehi, &
       fablo,fabhi, &
       bfact,bfact_top) &
      bind(c,name='fort_cgxdoty')

      use global_utility_module
      IMPLICIT NONE
      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3)
      INTEGER_T :: growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(ww)
      REAL_T, intent(in), target :: ww(DIMV(ww),ncomp)
      INTEGER_T, intent(in) :: DIMDEC(pp)
      REAL_T, intent(in), target :: pp(DIMV(pp),ncomp)
      REAL_T, intent(out) :: pw
!
      INTEGER_T i, j, k,veldir
!
      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.ge.1) then
       ! do nothing
      else
       print *,"bfact_top invalid"
       stop
      endif
      if ((ncomp.eq.1).or.(ncomp.eq.AMREX_SPACEDIM)) then
       ! do nothing
      else
       print *,"ncomp invalid"
       stop
      endif
      call checkbound_array(fablo,fabhi,pp,0,-1)
      call checkbound_array(fablo,fabhi,ww,0,-1)

      pw = zero
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       do veldir=1,ncomp
        pw = pw +  &
         pp(D_DECL(i,j,k),veldir)*ww(D_DECL(i,j,k),veldir)
       enddo ! veldir
      end do
      end do
      end do
!
      end subroutine fort_cgxdoty

      end module cpp_cg
