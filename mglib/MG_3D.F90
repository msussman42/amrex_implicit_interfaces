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

      module cpp_mg
      use amrex_fort_module, only : amrex_real

      contains
!-----------------------------------------------------------------------
      subroutine fort_average( &
           c,  &
           DIMS(c), &
           f,  &
           DIMS(f), &
           lo, hi, &
           iaverage, &
           bfact_coarse, &
           bfact_fine,bfact_top) &
      bind(c,name='fort_average')

      use global_utility_module

      IMPLICIT NONE
      integer, intent(in) :: iaverage,bfact_coarse,bfact_fine,bfact_top
      integer, intent(in) :: DIMDEC(c)
      integer, intent(in) :: DIMDEC(f)
      integer, intent(in) :: lo(AMREX_SPACEDIM)
      integer, intent(in) :: hi(AMREX_SPACEDIM)
      integer :: growlo(3),growhi(3)
      real(amrex_real), intent(in), target :: f(DIMV(f))
      real(amrex_real), intent(inout), target :: c(DIMV(c))
      real(amrex_real), pointer :: c_ptr(D_DECL(:,:,:))
!
      integer i, i2, i2p1
      integer j, j2, j2p1
      integer k, k2, k2p1
      real(amrex_real) denom

      c_ptr=>c

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
      call checkbound_array1(lo,hi,c_ptr,0,-1)

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

      end subroutine fort_average

!-----------------------------------------------------------------------
      subroutine fort_interp( &
       bfact,bfact_f,bfact_top, &
       f,  &
       DIMS(f), &
       c,  &
       DIMS(c), &
       lo, hi ) &
      bind(c,name='fort_interp')

      use global_utility_module
      IMPLICIT NONE
      integer, intent(in) :: DIMDEC(f)
      integer, intent(in) :: DIMDEC(c)
      integer, intent(in) :: lo(AMREX_SPACEDIM)
      integer, intent(in) :: hi(AMREX_SPACEDIM)
      integer growlo(3),growhi(3)
      real(amrex_real), intent(inout), target :: f(DIMV(f))
      real(amrex_real), pointer :: f_ptr(D_DECL(:,:,:))
      real(amrex_real), intent(in), target :: c(DIMV(c))
      real(amrex_real), pointer :: c_ptr(D_DECL(:,:,:))
      integer, intent(in) :: bfact,bfact_f,bfact_top
!
      integer i, i2, i2p1
      integer j, j2, j2p1
      integer k, k2, k2p1

      f_ptr=>f

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
      
      c_ptr=>c
      call checkbound_array1(lo,hi,c_ptr,0,-1)

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

      end subroutine fort_interp

      end module cpp_mg
