!     $Id: CG_3D.F,v 1.5 2000/08/24 16:02:46 car Exp $
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_SPACE.H"
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include "AMReX_ArrayLim.H"
#include "CG_F.H"

!-----------------------------------------------------------------------
!      
!     CGUPDATE: Modify the input arrays as follows:
!     
!     phi = yy + a*pp
!     
!-----------------------------------------------------------------------
      subroutine FORT_CGUPDATE( &
       phi,  &
       DIMS(phi), &
       a, &
       yy,   &
       DIMS(yy), &
       pp,   &
       DIMS(pp), &
       tilelo, tilehi, &
       fablo,fabhi,bfact,bfact_top)
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
      REAL_T, intent(in) :: pp(DIMV(pp))
      INTEGER_T, intent(in) :: DIMDEC(phi)
      REAL_T, intent(out) :: phi(DIMV(phi))
      INTEGER_T, intent(in) :: DIMDEC(yy)
      REAL_T, intent(in) :: yy(DIMV(yy))
      REAL_T, intent(in) :: a

      INTEGER_T i,j,k

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(phi),0,-1,4141)
      call checkbound(fablo,fabhi,DIMS(pp),0,-1,42)
      call checkbound(fablo,fabhi,DIMS(yy),0,-1,42)
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       phi(D_DECL(i,j,k)) = yy(D_DECL(i,j,k)) +  &
         a * pp(D_DECL(i,j,k))
      end do
      end do
      end do

      end subroutine FORT_CGUPDATE
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
      subroutine FORT_CGADVCP( &
       pp,  &
       DIMS(pp), &
       rr,  &
       DIMS(rr), &
       yy,  &
       DIMS(yy), &
       b, &
       tilelo, tilehi, &
       fablo,fabhi,bfact,bfact_top )
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
      REAL_T, intent(in) :: rr(DIMV(rr))
      INTEGER_T, intent(in) :: DIMDEC(pp)
      REAL_T, intent(out) :: pp(DIMV(pp))
      INTEGER_T, intent(in) :: DIMDEC(yy)
      REAL_T, intent(in) :: yy(DIMV(yy))
      REAL_T, intent(in) :: b

      INTEGER_T i,j,k
!
      call checkbound(fablo,fabhi, &
      DIMS(pp) &
      ,0,-1,51)
      call checkbound(fablo,fabhi, &
      DIMS(rr) &
      ,0,-1,52)
      call checkbound(fablo,fabhi, &
      DIMS(yy) &
      ,0,-1,52)

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
      end subroutine FORT_CGADVCP
      
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
      subroutine FORT_CGXDOTY( &
       ncomp, &
       pw, &
       pp,  &
       DIMS(pp), &
       ww,  &
       DIMS(ww), &
       tilelo, tilehi, &
       fablo,fabhi,bfact,bfact_top)
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
      REAL_T, intent(in) :: ww(DIMV(ww),ncomp)
      INTEGER_T, intent(in) :: DIMDEC(pp)
      REAL_T, intent(in) :: pp(DIMV(pp),ncomp)
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
      call checkbound(fablo,fabhi, &
      DIMS(pp) &
      ,0,-1,61)
      call checkbound(fablo,fabhi, &
      DIMS(ww) &
      ,0,-1,62)

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
      end subroutine FORT_CGXDOTY
