#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_SPACE.H"
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "LO_F.H"



      module cpp_lo

      contains

      subroutine fort_buildmat( &
       level, &
       veldir, &
       nsolve, &
       isweep, &
       solvemask, &  ! ones_mf
       DIMS(solvemask),  &
       a,DIMS(a),  &
       bx,DIMS(bx), &
       by,DIMS(by), &
       bz,DIMS(bz), &
       diagfab, &
       DIMS(work), &
       bxleft,bxright, &
       byleft,byright, &
       bzleft,bzright, &
       icbx,icby,icbz, &
       icdiag, &
       icdiagrb, &
       mask, &
       DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top) &
      bind(c,name='fort_buildmat')

      use global_utility_module
      IMPLICIT NONE
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: veldir
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: isweep
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3)
      INTEGER_T :: growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(solvemask)
      INTEGER_T, intent(in) :: DIMDEC(a)
      INTEGER_T, intent(in) :: DIMDEC(work)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(bx)
      INTEGER_T, intent(in) :: DIMDEC(by)
      INTEGER_T, intent(in) :: DIMDEC(bz)
      REAL_T, intent(inout), target :: bx(DIMV(bx))
      REAL_T, pointer :: bx_ptr(D_DECL(:,:,:))
      REAL_T, intent(inout), target :: by(DIMV(by))
      REAL_T, pointer :: by_ptr(D_DECL(:,:,:))
      REAL_T, intent(inout), target :: bz(DIMV(bz))
      REAL_T, pointer :: bz_ptr(D_DECL(:,:,:))
      REAL_T, intent(inout), target :: solvemask(DIMV(solvemask))
      REAL_T, pointer :: solvemask_ptr(D_DECL(:,:,:))
      REAL_T, intent(inout), target :: a(DIMV(a))
      REAL_T, pointer :: a_ptr(D_DECL(:,:,:))
      REAL_T, intent(inout), target :: diagfab(DIMV(work))
      REAL_T, pointer :: diagfab_ptr(D_DECL(:,:,:))
      REAL_T, intent(inout), target :: bxleft(DIMV(work))
      REAL_T, intent(inout), target :: bxright(DIMV(work))
      REAL_T, intent(inout), target :: byleft(DIMV(work))
      REAL_T, intent(inout), target :: byright(DIMV(work))
      REAL_T, intent(inout), target :: bzleft(DIMV(work))
      REAL_T, intent(inout), target :: bzright(DIMV(work))
      REAL_T, intent(inout), target :: icbx(DIMV(work))
      REAL_T, intent(inout), target :: icby(DIMV(work))
      REAL_T, intent(inout), target :: icbz(DIMV(work))
      REAL_T, intent(inout), target :: icdiag(DIMV(work))
      REAL_T, intent(out), target :: icdiagrb(DIMV(work))
      REAL_T, intent(out), target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))

      INTEGER_T i,j,k,ioff
      REAL_T test_mask
      REAL_T offdiagsum
      REAL_T local_diag
      REAL_T DD

      bx_ptr=>bx
      by_ptr=>by
      bz_ptr=>bz
      solvemask_ptr=>solvemask
      a_ptr=>a
      diagfab_ptr=>diagfab
      mask_ptr=>mask

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
      if (nsolve.ge.1) then
       ! do nothing
      else
       print *,"nsolve invalid"
       stop
      endif
      if ((veldir.ge.0).and.(veldir.lt.nsolve)) then
       ! do nothing
      else
       print *,"veldir invalid"
       stop
      endif
      if (level.ge.0) then
       ! do nothing
      else
       print *,"level invalid"
       stop
      endif

      call checkbound_array1(fablo,fabhi,solvemask_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,a_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,diagfab_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,bx_ptr,0,0)
      call checkbound_array1(fablo,fabhi,by_ptr,0,1)
      call checkbound_array1(fablo,fabhi,bz_ptr,0,AMREX_SPACEDIM-1)


      if (isweep.eq.0) then
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        solvemask(D_DECL(i,j,k))=one
       enddo ! i
       enddo ! j
       enddo ! k

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        ioff=MOD(i+j+k+4,2)
        mask(D_DECL(i,j,k))=1.0-ioff
        diagfab(D_DECL(i,j,k))=one  ! avoid divide by zero
       enddo
       enddo
       enddo

      else if (isweep.eq.1) then
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        bxleft(D_DECL(i,j,k))=bx(D_DECL(i,j,k))
        bxright(D_DECL(i,j,k))=bx(D_DECL(i+1,j,k))
        byleft(D_DECL(i,j,k))=by(D_DECL(i,j,k))
        byright(D_DECL(i,j,k))=by(D_DECL(i,j+1,k))
#if (AMREX_SPACEDIM==3)
        bzleft(D_DECL(i,j,k))=bz(D_DECL(i,j,k))
        bzright(D_DECL(i,j,k))=bz(i,j,k+1)
#endif
        offdiagsum= &
          bxleft(D_DECL(i,j,k))+bxright(D_DECL(i,j,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))+bzright(D_DECL(i,j,k)) &
#endif
          +byleft(D_DECL(i,j,k))+byright(D_DECL(i,j,k))

        local_diag=a(D_DECL(i,j,k))+offdiagsum
        diagfab(D_DECL(i,j,k))=local_diag

        test_mask=solvemask(D_DECL(i,j,k))

        if (offdiagsum.lt.zero) then
         print *,"offdiagsum invalid"
         stop
        else if (offdiagsum.eq.zero) then
         ! do nothing
        else if (offdiagsum.gt.zero) then
         ! do nothing
        else
         print *,"offdiagsum bust"
         stop
        endif

        if (a(D_DECL(i,j,k)).ge.zero) then
         ! do nothing
        else
         print *,"a should be nonneg"
         stop
        endif

        if (test_mask.eq.one) then
         ! do nothing
        else
         print *,"test_mask invalid"
         stop
        endif

        if (local_diag.gt.zero) then ! local_diag=a+offdiagsum
         ! do nothing
        else
         print *,"local_diag invalid"
         stop
        endif

       enddo
       enddo
       enddo
      else if (isweep.eq.2) then
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        test_mask=solvemask(D_DECL(i,j,k))
        local_diag=diagfab(D_DECL(i,j,k))
        if (local_diag.gt.zero) then
         ! do nothing
        else
         print *,"local_diag invalid"
         stop
        endif
        if (test_mask.eq.one) then
         ! do nothing
        else
         print *,"test_mask invalid"
         stop
        endif

        DD=diagfab(D_DECL(i,j,k))
        if (i.gt.tilelo(1)) then
         DD=DD-bxleft(D_DECL(i,j,k))* &
          bxleft(D_DECL(i,j,k))/icdiag(D_DECL(i-1,j,k))
        endif
        if (j.gt.tilelo(2)) then
         DD=DD-byleft(D_DECL(i,j,k))* &
          byleft(D_DECL(i,j,k))/icdiag(D_DECL(i,j-1,k))
        endif
#if (AMREX_SPACEDIM==3)
        if (k.gt.tilelo(AMREX_SPACEDIM)) then
         DD=DD-bzleft(D_DECL(i,j,k))* &
          bzleft(D_DECL(i,j,k))/icdiag(i,j,k-1)
        endif
#endif
        icdiag(D_DECL(i,j,k))=DD
        if (i.gt.tilelo(1)) then
         icbx(D_DECL(i,j,k))=-bxleft(D_DECL(i,j,k))/ &
          icdiag(D_DECL(i-1,j,k)) 
        endif
        if (j.gt.tilelo(2)) then
         icby(D_DECL(i,j,k))=-byleft(D_DECL(i,j,k))/ &
          icdiag(D_DECL(i,j-1,k)) 
        endif
#if (AMREX_SPACEDIM==3)
        if (k.gt.tilelo(AMREX_SPACEDIM)) then
         icbz(D_DECL(i,j,k))=-bzleft(D_DECL(i,j,k))/ &
          icdiag(i,j,k-1) 
        endif
#endif

       end do
       end do
       end do
      else if (isweep.eq.3) then
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        test_mask=solvemask(D_DECL(i,j,k))
        local_diag=diagfab(D_DECL(i,j,k))
        if (test_mask.eq.one) then
         ! do nothing
        else
         print *,"test_mask invalid"
         stop
        endif

        if (local_diag.gt.zero) then ! local_diag=a+offdiagsum
         ! do nothing
        else
         print *,"local_diag invalid"
         stop
        endif

        DD=diagfab(D_DECL(i,j,k))

        if (i.gt.tilelo(1)) then

         local_diag=diagfab(D_DECL(i-1,j,k))
         if (local_diag.gt.zero) then ! local_diag=a+offdiagsum
          ! do nothing
         else
          print *,"local_diag invalid"
          stop
         endif

         DD=DD-bxleft(D_DECL(i,j,k))* &
          bxleft(D_DECL(i,j,k))/local_diag

        endif

        if (j.gt.tilelo(2)) then

         local_diag=diagfab(D_DECL(i,j-1,k))
         if (local_diag.gt.zero) then ! local_diag=a+offdiagsum
          ! do nothing
         else
          print *,"local_diag invalid"
          stop
         endif

         DD=DD-byleft(D_DECL(i,j,k))* &
          byleft(D_DECL(i,j,k))/local_diag

        endif
#if (AMREX_SPACEDIM==3)
        if (k.gt.tilelo(AMREX_SPACEDIM)) then

         local_diag=diagfab(D_DECL(i,j,k-1))
         if (local_diag.gt.zero) then ! local_diag=a+offdiagsum
          ! do nothing
         else
          print *,"local_diag invalid"
          stop
         endif

         DD=DD-bzleft(D_DECL(i,j,k))* &
          bzleft(D_DECL(i,j,k))/local_diag

        endif
#endif
        if (i.lt.tilehi(1)) then

         local_diag=diagfab(D_DECL(i+1,j,k))
         if (local_diag.gt.zero) then ! local_diag=a+offdiagsum
          ! do nothing
         else
          print *,"local_diag invalid"
          stop
         endif

         DD=DD-bxright(D_DECL(i,j,k))* &
          bxright(D_DECL(i,j,k))/local_diag 

        endif
        if (j.lt.tilehi(2)) then

         local_diag=diagfab(D_DECL(i,j+1,k))
         if (local_diag.gt.zero) then ! local_diag=a+offdiagsum
          ! do nothing
         else
          print *,"local_diag invalid"
          stop
         endif

         DD=DD-byright(D_DECL(i,j,k))* &
          byright(D_DECL(i,j,k))/local_diag

        endif
#if (AMREX_SPACEDIM==3)
        if (k.lt.tilehi(AMREX_SPACEDIM)) then

         local_diag=diagfab(D_DECL(i,j,k+1))
         if (local_diag.gt.zero) then ! local_diag=a+offdiagsum
          ! do nothing
         else
          print *,"local_diag invalid"
          stop
         endif

         DD=DD-bzright(D_DECL(i,j,k))* &
          bzright(D_DECL(i,j,k))/local_diag

        endif
#endif

        icdiagrb(D_DECL(i,j,k))=DD
       enddo
       enddo
       enddo
      else
       print *,"isweep invalid"
       stop
      endif

      return
      end subroutine fort_buildmat


      subroutine fort_residl( &
       level, &
       mg_coarsest_level, &
       nsolve, &
       masksing, &
       DIMS(masksing),  &
       res,DIMS(res),  &
       rhs,DIMS(rhs), &
       phi,DIMS(phi), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top) &
      bind(c,name='fort_residl')

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: mg_coarsest_level
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3)
      INTEGER_T :: growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(masksing)
      REAL_T, intent(in), target :: masksing(DIMV(masksing))
      INTEGER_T, intent(in) :: DIMDEC(phi)
      REAL_T, intent(in), target :: phi(DIMV(phi),nsolve)
      INTEGER_T, intent(in) :: DIMDEC(rhs)
      REAL_T, intent(in), target :: rhs(DIMV(rhs),nsolve)
      INTEGER_T, intent(in) :: DIMDEC(res)
      REAL_T, intent(out), target :: res(DIMV(res),nsolve)
      REAL_T, pointer :: res_ptr(D_DECL(:,:,:),:)
!
      INTEGER_T i,j,k,veldir
      REAL_T test_mask

      res_ptr=>res

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      if (nsolve.le.0) then
       print *,"nsolve invalid"
       stop
      endif
      if ((level.ge.0).and.(level.le.mg_coarsest_level)) then
       ! do nothing
      else
       print *,"level or mg_coarsest_level invalid"
       stop
      endif
      call checkbound_array1(fablo,fabhi,masksing,1,-1)
      call checkbound_array(fablo,fabhi,rhs,0,-1)
      call checkbound_array(fablo,fabhi,res_ptr,0,-1)
      call checkbound_array(fablo,fabhi,phi,0,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       test_mask=masksing(D_DECL(i,j,k))
       do veldir=1,nsolve
        if (test_mask.eq.one) then
         res(D_DECL(i,j,k),veldir) = rhs(D_DECL(i,j,k),veldir) - &
          phi(D_DECL(i,j,k),veldir)
        else
         print *,"test_mask invalid in RESIDL"
         stop
        endif
       enddo ! veldir
      enddo
      enddo
      enddo

      return
      end subroutine fort_residl

      subroutine fort_averageec( &
       nsolve, &
       c,  &
       DIMS(c), &
       f,  &
       DIMS(f), &
       lo, hi, &
       cdir,avg, &
       bfact_coarse,bfact_fine,bfact_top) &
      bind(c,name='fort_averageec')

      use global_utility_module
      IMPLICIT NONE
!
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: bfact_coarse,bfact_fine,bfact_top
      INTEGER_T, intent(in) :: avg
      INTEGER_T, intent(in) :: lo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: hi(AMREX_SPACEDIM)
      INTEGER_T growlo(3)
      INTEGER_T growhi(3)
      INTEGER_T, intent(in) :: cdir
      INTEGER_T, intent(in) :: DIMDEC(f)
      REAL_T, intent(in), target :: f(DIMV(f),nsolve)
      INTEGER_T, intent(in) :: DIMDEC(c)
      REAL_T, intent(out), target :: c(DIMV(c),nsolve)
      REAL_T, pointer :: c_ptr(D_DECL(:,:,:),:)
!     
      INTEGER_T i,j,k,veldir
      REAL_T denom

      c_ptr=>c

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
      if (nsolve.le.0) then
       print *,"nsolve invalid"
       stop
      endif
      call checkbound_array(lo,hi,c_ptr,0,cdir)

      if (avg.eq.1) then
       if (AMREX_SPACEDIM.eq.3) then
        denom=0.25d0
       else if (AMREX_SPACEDIM.eq.2) then
        denom=half
       else
        print *,"dimension bust"
        stop
       endif
      else if (avg.eq.0) then
       denom=one
      else
       print *,"avg invalid"
       stop
      endif

      if ( cdir .eq. 0 ) then
       call growntileboxMAC(lo,hi,lo,hi,growlo,growhi,0,cdir,101) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        do veldir=1,nsolve
         if (AMREX_SPACEDIM.eq.3) then
          c(D_DECL(i,j,k),veldir) = denom*( &
           + f(D_DECL(2*i,2*j,2*k),veldir) &
           + f(D_DECL(2*i,2*j+1,2*k),veldir) &
           + f(D_DECL(2*i,2*j,2*k+1),veldir) &
           + f(D_DECL(2*i,2*j+1,2*k+1),veldir))
         else if (AMREX_SPACEDIM.eq.2) then
          c(D_DECL(i,j,k),veldir) = denom*( &
           + f(D_DECL(2*i,2*j,2*k),veldir) &
           + f(D_DECL(2*i,2*j+1,2*k),veldir))
         else
          print *,"dimension bust"
          stop
         endif
        enddo ! veldir
       enddo
       enddo
       enddo
      else if (cdir .eq. 1 ) then
       call growntileboxMAC(lo,hi,lo,hi,growlo,growhi,0,cdir,102) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        do veldir=1,nsolve
         if (AMREX_SPACEDIM.eq.3) then
          c(D_DECL(i,j,k),veldir) = denom*( &
           + f(D_DECL(2*i,2*j,2*k),veldir) &
           + f(D_DECL(2*i+1,2*j,2*k),veldir) &
           + f(D_DECL(2*i,2*j,2*k+1),veldir) &
           + f(D_DECL(2*i+1,2*j,2*k+1),veldir))
         else if (AMREX_SPACEDIM.eq.2) then
          c(D_DECL(i,j,k),veldir) = denom*( &
           + f(D_DECL(2*i,2*j,2*k),veldir) &
           + f(D_DECL(2*i+1,2*j,2*k),veldir))
         else
          print *,"dimension bust"
          stop
         endif
        enddo ! veldir
       enddo
       enddo
       enddo
      else if ((cdir .eq. 2 ).and.(AMREX_SPACEDIM.eq.3)) then
       call growntileboxMAC(lo,hi,lo,hi,growlo,growhi,0,cdir,103) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        do veldir=1,nsolve
         if (AMREX_SPACEDIM.eq.3) then
          c(D_DECL(i,j,k),veldir) = denom*( &
           + f(D_DECL(2*i,2*j,2*k),veldir) &
           + f(D_DECL(2*i+1,2*j,2*k),veldir) &
           + f(D_DECL(2*i,2*j+1,2*k),veldir) &
           + f(D_DECL(2*i+1,2*j+1,2*k),veldir))
         else
          print *,"dimension bust"
          stop
         endif
        enddo ! veldir
       enddo
       enddo
       enddo
      else
       print *,"cdir invalid"
       stop
      endif
!     
      end subroutine fort_averageec
!-----------------------------------------------------------------------

      subroutine fort_averagecc( &
       nsolve, &
       ncomp_expect, &
       c, &
       DIMS(c), &
       f,  &
       DIMS(f), &
       lo, hi, avg, &
       ngrow, &
       bfact_coarse,bfact_fine,bfact_top) &
      bind(c,name='fort_averagecc')


      use global_utility_module 
      IMPLICIT NONE
!
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: ncomp_expect
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: bfact_coarse,bfact_fine,bfact_top
      INTEGER_T, intent(in) :: avg
      INTEGER_T, intent(in) :: DIMDEC(f)
      INTEGER_T, intent(in) :: DIMDEC(c)
      INTEGER_T, intent(in) :: lo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: hi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3)
      INTEGER_T :: growhi(3)
      REAL_T, intent(in), target :: f(DIMV(f),ncomp_expect)
      REAL_T, intent(out), target :: c(DIMV(c),ncomp_expect)
      REAL_T, pointer :: c_ptr(D_DECL(:,:,:),:)
!
      INTEGER_T i,j,k,veldir
      REAL_T one_over_denom
      INTEGER_T sum_mask,max_sum

      c_ptr=>c
      
      if (bfact_coarse.ge.1) then
       ! do nothing
      else
       print *,"bfact_coarse invalid"
       stop
      endif
      if (bfact_fine.ge.1) then
       ! do nothing
      else
       print *,"bfact_fine invalid"
       stop
      endif
      if (bfact_top.ge.1) then
       ! do nothing
      else
       print *,"bfact_top invalid"
       stop
      endif
      if (nsolve.ge.1) then
       ! do nothing
      else
       print *,"nsolve invalid"
       stop
      endif
      if ((avg.eq.0).or. &  ! one_over_denom=1
          (avg.eq.1)) then  ! one_over_denom=1/2^dim
       if (nsolve.eq.ncomp_expect) then
        ! do nothing
       else
        print *,"nsolve.ne.ncomp_expect"
        stop
       endif
       if (ngrow.eq.0) then
        ! do nothing
       else
        print *,"ngrow.ne.0"
        stop
       endif
      else if (avg.eq.2) then ! ones_mf
       if (ncomp_expect.eq.1) then
        ! do nothing
       else
        print *,"ncomp_expect.ne.1"
        stop
       endif
       if (ngrow.eq.1) then
        ! do nothing
       else
        print *,"ngrow.ne.1"
        stop
       endif
      else
       print *,"avg invalid"
       stop
      endif

      call checkbound_array(lo,hi,c_ptr,ngrow,-1)

      if (avg.eq.1) then
       if (AMREX_SPACEDIM.eq.3) then
        one_over_denom=0.125d0
       else if (AMREX_SPACEDIM.eq.2) then
        one_over_denom=0.25d0
       else
        print *,"dimension bust"
        stop
       endif
      else if (avg.eq.0) then
       one_over_denom=one
      else if (avg.eq.2) then
       one_over_denom=one
      else
       print *,"avg invalid"
       stop
      endif

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       do veldir=1,ncomp_expect

        if (AMREX_SPACEDIM.eq.3) then
         c(D_DECL(i,j,k),veldir) =  one_over_denom*( &
          + f(D_DECL(2*i+1,2*j+1,2*k),veldir) &
          + f(D_DECL(2*i,2*j+1,2*k),veldir) &
          + f(D_DECL(2*i+1,2*j,2*k),veldir) &
          + f(D_DECL(2*i,2*j,2*k),veldir) &
          + f(D_DECL(2*i+1,2*j+1,2*k+1),veldir) &
          + f(D_DECL(2*i,2*j+1,2*k+1),veldir) &
          + f(D_DECL(2*i+1,2*j,2*k+1),veldir) &
          + f(D_DECL(2*i,2*j,2*k+1),veldir) )
        else if (AMREX_SPACEDIM.eq.2) then
         c(D_DECL(i,j,k),veldir) =  one_over_denom*( &
          + f(D_DECL(2*i+1,2*j+1,2*k),veldir) &
          + f(D_DECL(2*i,2*j+1,2*k),veldir) &
          + f(D_DECL(2*i+1,2*j,2*k),veldir) &
          + f(D_DECL(2*i,2*j,2*k),veldir))
        else
         print *,"dimension bust"
         stop
        endif
        if (avg.eq.2) then ! ones_mf
         sum_mask=NINT(c(D_DECL(i,j,k),veldir))
         if (AMREX_SPACEDIM.eq.3) then
          max_sum=8
         else if (AMREX_SPACEDIM.eq.2) then
          max_sum=4
         else
          print *,"dimension bust"
          stop
         endif
         if (sum_mask.eq.max_sum) then
          sum_mask=1
         else if ((sum_mask.ge.0).and.(sum_mask.le.max_sum-1)) then
          sum_mask=0
         else
          print *,"sum_mask invalid"
          stop
         endif
         c(D_DECL(i,j,k),veldir)=sum_mask
        else if ((avg.eq.0).or.(avg.eq.1)) then
         ! do nothing
        else
         print *,"avg invalid"
         stop
        endif

       enddo ! veldir=1,ncomp_expect

      enddo
      enddo
      enddo

      end subroutine fort_averagecc



!-----------------------------------------------------------------------


        ! NO TILING 
      subroutine fort_applybc( &
       nsolve, &
       phi,  &
       DIMS(phi), &
       bfab,  &
       DIMS(bfab), &
       mfab,  &
       DIMS(mfab), &
       bcpres, &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top) &
      bind(c,name='fort_applybc')

      use global_utility_module
      IMPLICIT NONE
!
!     mask=0 => fine/fine grid interface....
!
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T :: dir,side
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3)
      INTEGER_T :: growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(phi)
      INTEGER_T, intent(in) :: DIMDEC(bfab)
      INTEGER_T, intent(in) :: DIMDEC(mfab)
      INTEGER_T, intent(in) :: bcpres(AMREX_SPACEDIM,2,nsolve)
      REAL_T, intent(inout), target :: phi(DIMV(phi),nsolve)
      REAL_T, pointer :: phi_ptr(D_DECL(:,:,:),:)

      REAL_T, intent(in), target :: bfab(DIMV(bfab),nsolve)
      REAL_T, intent(in), target :: mfab(DIMV(mfab))
!
      INTEGER_T i,j,k
      INTEGER_T ib,jb,kb
      INTEGER_T iout,jout,kout
      INTEGER_T ireflect,jreflect,kreflect
      INTEGER_T ii,jj,kk
      INTEGER_T bct,veldir,dir_bc
 
      phi_ptr=>phi

      call checkbound_array(fablo,fabhi,phi_ptr,1,-1)
      call checkbound_array(fablo,fabhi,bfab,1,-1)
      call checkbound_array1(fablo,fabhi,mfab,1,-1)

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      if (nsolve.le.0) then
       print *,"nsolve invalid"
       stop
      endif

      do veldir=1,nsolve

       do dir=1,AMREX_SPACEDIM
        ii=0
        jj=0
        kk=0
        if (dir.eq.1) then
         ii=1
        else if (dir.eq.2) then
         jj=1
        else if ((dir.eq.3).and.(AMREX_SPACEDIM.eq.3)) then
         kk=1
        else
         print *,"dir invalid apply bc"
         stop
        endif

        do side=1,2
         bct=bcpres(dir,side,veldir)
         if ((bct.ne.INT_DIR).and. &
             (bct.ne.EXT_DIR).and. &
             (bct.ne.REFLECT_EVEN).and. &
             (bct.ne.FOEXTRAP).and. &
             (bct.ne.REFLECT_ODD)) then
          print *,"bct invalid bct=",bct
          stop
         endif
          ! NO TILING 
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         if (side.eq.1) then
          growhi(dir)=growlo(dir)
         else if (side.eq.2) then
          growlo(dir)=growhi(dir)
         else
          print *,"side invalid"
          stop
         endif
         do dir_bc=1,AMREX_SPACEDIM
          if (dir_bc.ne.dir) then
           growlo(dir_bc)=growlo(dir_bc)-1
           growhi(dir_bc)=growhi(dir_bc)+1
          endif
         enddo

         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          if (side.eq.1) then
           ib=i-ii
           jb=j-jj
           kb=k-kk
          else if (side.eq.2) then
           ib=i+ii
           jb=j+jj
           kb=k+kk
          else
           print *,"side invalid"
           stop
          endif

          if (side.eq.1) then
           iout=i-ii
           jout=j-jj
           kout=k-kk
           ireflect=i
           jreflect=j
           kreflect=k
          else if (side.eq.2) then
           iout=i+ii
           jout=j+jj
           kout=k+kk
           ireflect=i
           jreflect=j
           kreflect=k
          else
           print *,"side invalid"
           stop
          endif
 
           ! mask=0 at fine/fine

          if (bct.eq.EXT_DIR) then
           phi(D_DECL(iout,jout,kout),veldir)= &
             bfab(D_DECL(ib,jb,kb),veldir)
          else if (bct.eq.INT_DIR) then
            ! coarse/fine
           if (mfab(D_DECL(iout,jout,kout)).ne.zero) then
            phi(D_DECL(iout,jout,kout),veldir)= &
              bfab(D_DECL(iout,jout,kout),veldir)
           endif
          else if (bct.eq.REFLECT_ODD) then
           phi(D_DECL(iout,jout,kout),veldir)= &
            -phi(D_DECL(ireflect,jreflect,kreflect),veldir)
          else if (bct.eq.REFLECT_EVEN) then
           phi(D_DECL(iout,jout,kout),veldir)= &
             phi(D_DECL(ireflect,jreflect,kreflect),veldir)
          else if (bct.eq.FOEXTRAP) then
           phi(D_DECL(iout,jout,kout),veldir)=phi(D_DECL(i,j,k),veldir)
          else
           phi(D_DECL(iout,jout,kout),veldir)=zero
          endif
         enddo ! k
         enddo ! j
         enddo ! i
        enddo ! side
       enddo ! dir
      enddo ! veldir

      return
      end subroutine fort_applybc

      end module cpp_lo
