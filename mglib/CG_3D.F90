#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_SPACE.H"
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include "AMReX_ArrayLim.H"
#include "CG_F.H"

      module cpp_cg
      use amrex_fort_module, only : amrex_real

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

      integer, intent(in) :: tilelo(AMREX_SPACEDIM)
      integer, intent(in) :: tilehi(AMREX_SPACEDIM)
      integer, intent(in) :: fablo(AMREX_SPACEDIM)
      integer, intent(in) :: fabhi(AMREX_SPACEDIM)
      integer :: growlo(3)
      integer :: growhi(3)
      integer, intent(in) :: bfact,bfact_top
      integer, intent(in) :: DIMDEC(pp)
      real(amrex_real), intent(in), target :: pp(DIMV(pp))
      integer, intent(in) :: DIMDEC(phi)
      real(amrex_real), intent(out), target :: phi(DIMV(phi))
      real(amrex_real), pointer :: phi_ptr(D_DECL(:,:,:))

      integer, intent(in) :: DIMDEC(yy)
      real(amrex_real), intent(in), target :: yy(DIMV(yy))
      real(amrex_real), intent(in) :: a

      integer i,j,k

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
      integer, intent(in) :: tilelo(AMREX_SPACEDIM)
      integer, intent(in) :: tilehi(AMREX_SPACEDIM)
      integer, intent(in) :: fablo(AMREX_SPACEDIM)
      integer, intent(in) :: fabhi(AMREX_SPACEDIM)
      integer :: growlo(3)
      integer :: growhi(3)
      integer, intent(in) :: bfact,bfact_top
      integer, intent(in) :: DIMDEC(rr)
      real(amrex_real), intent(in), target :: rr(DIMV(rr))
      integer, intent(in) :: DIMDEC(pp)
      real(amrex_real), intent(out), target :: pp(DIMV(pp))
      real(amrex_real), pointer :: pp_ptr(D_DECL(:,:,:))

      integer, intent(in) :: DIMDEC(yy)
      real(amrex_real), intent(in), target :: yy(DIMV(yy))
      real(amrex_real), intent(in) :: b

      integer i,j,k

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
      integer, intent(in) :: ncomp
      integer, intent(in) :: tilelo(AMREX_SPACEDIM)
      integer, intent(in) :: tilehi(AMREX_SPACEDIM)
      integer, intent(in) :: fablo(AMREX_SPACEDIM)
      integer, intent(in) :: fabhi(AMREX_SPACEDIM)
      integer :: growlo(3)
      integer :: growhi(3)
      integer, intent(in) :: bfact,bfact_top
      integer, intent(in) :: DIMDEC(ww)
      real(amrex_real), intent(in), target :: ww(DIMV(ww),ncomp)
      integer, intent(in) :: DIMDEC(pp)
      real(amrex_real), intent(in), target :: pp(DIMV(pp),ncomp)
      real(amrex_real), intent(out) :: pw
!
      integer i, j, k,veldir
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
