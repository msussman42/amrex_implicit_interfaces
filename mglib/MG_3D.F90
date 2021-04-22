!
! $Id: MG_3D.F,v 1.2 1998/12/14 20:23:47 lijewski Exp $
!
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_SPACE.H"
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include "AMReX_ArrayLim.H"
#include "MG_F.H"

!-----------------------------------------------------------------------
      subroutine FORT_AVERAGE ( &
           c,  &
           DIMS(c), &
           f,  &
           DIMS(f), &
           lo, hi, &
           iaverage, &
           bfact_coarse, &
           bfact_fine,bfact_top)
      use global_utility_module

      IMPLICIT NONE
      INTEGER_T iaverage,bfact_coarse,bfact_fine,bfact_top
      INTEGER_T DIMDEC(c)
      INTEGER_T DIMDEC(f)
      INTEGER_T lo(AMREX_SPACEDIM)
      INTEGER_T hi(AMREX_SPACEDIM)
      INTEGER_T growlo(3),growhi(3)
      REAL_T f(DIMV(f))
      REAL_T c(DIMV(c))
!
      INTEGER_T i, i2, i2p1
      INTEGER_T j, j2, j2p1
      INTEGER_T k, k2, k2p1
      REAL_T denom
!
      if (iaverage.eq.1) then
       if (AMREX_SPACEDIM.eq.3) then
        denom=eighth
       else if (AMREX_SPACEDIM.eq.2) then
        denom=fourth
       else
        print *,"dimension bust"
        stop
       endif
      else if (iaverage.eq.0) then
       denom=one
      else
       print *,"iaverage invalid"
       stop
      endif
      if (bfact_coarse.lt.1) then
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.lt.1) then
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      call checkbound(lo,hi, &
       DIMS(c), &
       0,-1,130)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
       k2 = 2*k
       k2p1 = k2 + 1
       do j=growlo(2),growhi(2)
        j2 = 2*j
        j2p1 = j2 + 1
        do i=growlo(1),growhi(1)
         i2 = 2*i
         i2p1 = i2 + 1
         if (AMREX_SPACEDIM.eq.3) then
          c(D_DECL(i,j,k)) =  ( &
            + f(D_DECL(i2p1,j2p1,k2))  &
            + f(D_DECL(i2,j2p1,k2)) &
            + f(D_DECL(i2p1,j2,k2)) &
            + f(D_DECL(i2,j2,k2)) &
            + f(D_DECL(i2p1,j2p1,k2p1)) &
            + f(D_DECL(i2,j2p1,k2p1)) &
            + f(D_DECL(i2p1,j2,k2p1)) &
            + f(D_DECL(i2,j2,k2p1)) )*denom
         else if (AMREX_SPACEDIM.eq.2) then
          c(D_DECL(i,j,k)) =  ( &
            + f(D_DECL(i2p1,j2p1,k2))  &
            + f(D_DECL(i2,j2p1,k2)) &
            + f(D_DECL(i2p1,j2,k2)) &
            + f(D_DECL(i2,j2,k2)) )*denom
         else
          print *,"dimension bust"
          stop
         endif
        end do
       end do
      end do

      end subroutine FORT_AVERAGE

!-----------------------------------------------------------------------
      subroutine FORT_INTERP ( &
       bfact,bfact_f,bfact_top, &
       f,  &
       DIMS(f), &
       c,  &
       DIMS(c), &
       lo, hi )
      use global_utility_module
      IMPLICIT NONE
      INTEGER_T DIMDEC(f)
      INTEGER_T DIMDEC(c)
      INTEGER_T lo(AMREX_SPACEDIM)
      INTEGER_T hi(AMREX_SPACEDIM)
      INTEGER_T growlo(3),growhi(3)
      REAL_T f(DIMV(f))
      REAL_T c(DIMV(c))
      INTEGER_T bfact,bfact_f,bfact_top
!
      INTEGER_T i, i2, i2p1
      INTEGER_T j, j2, j2p1
      INTEGER_T k, k2, k2p1

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_f.lt.1) then
       print *,"bfact_f invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      
!
      call checkbound(lo,hi, &
       DIMS(c), &
       0,-1,140)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
       k2 = 2*k
       k2p1 = k2 + 1
       do j=growlo(2),growhi(2)
        j2 = 2*j
        j2p1 = j2 + 1
        do i=growlo(1),growhi(1)
         i2 = 2*i
         i2p1 = i2 + 1
 
         f(D_DECL(i2p1,j2p1,k2)) = c(D_DECL(i,j,k)) + &
             f(D_DECL(i2p1,j2p1,k2))
         f(D_DECL(i2,j2p1,k2)) = c(D_DECL(i,j,k)) + &
             f(D_DECL(i2,j2p1,k2))
         f(D_DECL(i2p1,j2,k2)) = c(D_DECL(i,j,k)) + &
             f(D_DECL(i2p1,j2,k2))
         f(D_DECL(i2,j2,k2)) = c(D_DECL(i,j,k)) + &
             f(D_DECL(i2,j2,k2))

         if (AMREX_SPACEDIM.eq.3) then
          f(D_DECL(i2p1,j2p1,k2p1)) = c(D_DECL(i,j,k)) + &
             f(D_DECL(i2p1,j2p1,k2p1))
          f(D_DECL(i2,j2p1,k2p1)) = c(D_DECL(i,j,k)) + &
             f(D_DECL(i2,j2p1,k2p1))
          f(D_DECL(i2p1,j2,k2p1)) = c(D_DECL(i,j,k)) + &
             f(D_DECL(i2p1,j2,k2p1))
          f(D_DECL(i2,j2,k2p1)) = c(D_DECL(i,j,k)) + &
             f(D_DECL(i2,j2,k2p1))
         else if (AMREX_SPACEDIM.eq.2) then
          ! do nothing
         else
          print *,"dimension bust"
          stop
         endif
        end do
       end do
      end do

      end subroutine FORT_INTERP


      subroutine FORT_CHECKERBOARD_RB( &
       ncomp, &
       v,  &
       DIMS(v), &
       tilelo, tilehi, &
       fablo,fabhi, &
       bfact,bfact_top)
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
      INTEGER_T, intent(in) :: DIMDEC(v)
      REAL_T, intent(inout) :: v(DIMV(v),ncomp)

      INTEGER_T i, j, k, veldir,ioff

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
      DIMS(v) &
      ,0,-1,61)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       do veldir=1,ncomp
        ioff=MOD(i+j+k+4,2)
        if (ioff.eq.0) then
         v(D_DECL(i,j,k),veldir)=-one
        else if (ioff.eq.1) then
         v(D_DECL(i,j,k),veldir)=one
        else
         print *,"ioff invalid"
         stop
        endif
       enddo ! veldir
      end do
      end do
      end do

      end subroutine FORT_CHECKERBOARD_RB


