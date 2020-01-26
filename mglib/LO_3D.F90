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



!-----------------------------------------------------------------------


      subroutine FORT_BUILDMAT( &
       level, &
       veldir, &
       nsolve, &
       isweep, &
       offdiag_coeff, &
       check_for_singular, &
       diag_regularization, &
       solvemask, &  ! ones_mf
       DIMS(solvemask),  &
       adual,DIMS(adual),  &
       a,DIMS(a),  &
       bx,DIMS(bx), &
       by,DIMS(by), &
       bz,DIMS(bz), &
       diag_non_sing, &
       DIMS(work), &
       diag_dual, &
       diag_nodual, &
       bxleft,bxright, &
       byleft,byright, &
       bzleft,bzright, &
       icbx,icby,icbz,icdiag,icdiagrb, &
       mask, &
       DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top)
      use global_utility_module
      IMPLICIT NONE
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: veldir
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: isweep
      INTEGER_T, intent(in) :: check_for_singular
      REAL_T, intent(in) :: offdiag_coeff
      REAL_T, intent(in) :: diag_regularization
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3)
      INTEGER_T :: growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(solvemask)
      INTEGER_T, intent(in) :: DIMDEC(adual)
      INTEGER_T, intent(in) :: DIMDEC(a)
      INTEGER_T, intent(in) :: DIMDEC(work)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(bx)
      INTEGER_T, intent(in) :: DIMDEC(by)
      INTEGER_T, intent(in) :: DIMDEC(bz)
      REAL_T, intent(inout) :: bx(DIMV(bx))
      REAL_T, intent(inout) :: by(DIMV(by))
      REAL_T, intent(inout) :: bz(DIMV(bz))
      REAL_T, intent(inout) :: solvemask(DIMV(solvemask))
      REAL_T, intent(inout) :: adual(DIMV(adual))
      REAL_T, intent(inout) :: a(DIMV(a))
      REAL_T, intent(inout) :: diag_non_sing(DIMV(work))
      REAL_T, intent(inout) :: diag_dual(DIMV(work))
      REAL_T, intent(inout) :: diag_nodual(DIMV(work))
      REAL_T, intent(inout) :: bxleft(DIMV(work))
      REAL_T, intent(inout) :: bxright(DIMV(work))
      REAL_T, intent(inout) :: byleft(DIMV(work))
      REAL_T, intent(inout) :: byright(DIMV(work))
      REAL_T, intent(inout) :: bzleft(DIMV(work))
      REAL_T, intent(inout) :: bzright(DIMV(work))
      REAL_T, intent(inout) :: icbx(DIMV(work))
      REAL_T, intent(inout) :: icby(DIMV(work))
      REAL_T, intent(inout) :: icbz(DIMV(work))
      REAL_T, intent(inout) :: icdiag(DIMV(work))
      REAL_T, intent(out) :: icdiagrb(DIMV(work))
      REAL_T, intent(out) :: mask(DIMV(mask))

      INTEGER_T i,j,k,ioff
      REAL_T offdiagsum
      REAL_T local_diag_NONSING
      REAL_T local_diag_dual
      REAL_T local_diag_nodual
      REAL_T test_mask
      REAL_T DD
      REAL_T dtau_factor,dtau
      REAL_T mask_cen
      REAL_T maskXP,maskYP,maskZP
      REAL_T maskXM,maskYM,maskZM

      dtau_factor=0.001

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
      if ((check_for_singular.eq.0).or.(check_for_singular.eq.1)) then
       ! do nothing
      else
       print *,"check_for_singular invalid"
       stop
      endif
      if (offdiag_coeff.gt.zero) then
       ! do nothing
      else
       print *,"offdiag_coeff invalid"
       stop
      endif
      if ((diag_regularization.gt.zero).and. &
          (diag_regularization.le.1.0D-3)) then
       ! do nothing
      else
       print *,"diag_regularization invalid"
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

      call checkbound(fablo,fabhi,DIMS(solvemask),1,-1,81)
      call checkbound(fablo,fabhi,DIMS(adual),0,-1,81)
      call checkbound(fablo,fabhi,DIMS(a),0,-1,81)
      call checkbound(fablo,fabhi,DIMS(work),1,-1,81)
      call checkbound(fablo,fabhi,DIMS(mask),1,-1,81)
      call checkbound(fablo,fabhi,DIMS(bx),0,0,81)
      call checkbound(fablo,fabhi,DIMS(by),0,1,81)
      call checkbound(fablo,fabhi,DIMS(bz),0,AMREX_SPACEDIM-1,81)


      if (isweep.eq.0) then
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        if (level.eq.0) then ! finest level
         if (veldir.eq.0) then
          if (solvemask(D_DECL(i,j,k)).eq.one) then
           ! do nothing
          else
           print *,"solvemask not initialized properly"
           stop
          endif
         else if ((veldir.gt.0).and.(veldir.lt.nsolve)) then
          if ((solvemask(D_DECL(i,j,k)).eq.one).or. &
              (solvemask(D_DECL(i,j,k)).eq.zero)) then
           ! do nothing
          else
           print *,"solvemask not initialized properly"
           stop
          endif
         else
          print *,"veldir invalid"
          stop
         endif
        else if (level.gt.0) then
         if (veldir.eq.0) then
          mask_cen=solvemask(D_DECL(i,j,k))
          maskXP=solvemask(D_DECL(i+1,j,k))
          maskXM=solvemask(D_DECL(i-1,j,k))
          maskYP=solvemask(D_DECL(i,j+1,k))
          maskYM=solvemask(D_DECL(i,j-1,k))
          maskZP=solvemask(D_DECL(i,j,k+1))
          maskZM=solvemask(D_DECL(i,j,k-1))
          if (mask_cen.eq.zero) then
           adual(D_DECL(i,j,k))=zero
           a(D_DECL(i,j,k))=zero
           bx(D_DECL(i,j,k))=zero
           bx(D_DECL(i+1,j,k))=zero
           by(D_DECL(i,j,k))=zero
           by(D_DECL(i,j+1,k))=zero
#if (AMREX_SPACEDIM==3)
           bz(D_DECL(i,j,k))=zero
           bz(D_DECL(i,j,k+1))=zero
#endif
          else if (mask_cen.eq.one) then
           if (maskXP.eq.zero) then
            bx(D_DECL(i+1,j,k))=zero
           else if (maskXP.eq.one) then
            ! do nothing
           else
            print *,"maskXP invalid"
            stop
           endif
           if (maskXM.eq.zero) then
            bx(D_DECL(i,j,k))=zero
           else if (maskXM.eq.one) then
            ! do nothing
           else
            print *,"maskXM invalid"
            stop
           endif
           if (maskYP.eq.zero) then
            by(D_DECL(i,j+1,k))=zero
           else if (maskYP.eq.one) then
            ! do nothing
           else
            print *,"maskYP invalid"
            stop
           endif
           if (maskYM.eq.zero) then
            by(D_DECL(i,j,k))=zero
           else if (maskYM.eq.one) then
            ! do nothing
           else
            print *,"maskYM invalid"
            stop
           endif
#if (AMREX_SPACEDIM==3)
           if (maskZP.eq.zero) then
            bz(D_DECL(i,j,k+1))=zero
           else if (maskZP.eq.one) then
            ! do nothing
           else
            print *,"maskZP invalid"
            stop
           endif
           if (maskZM.eq.zero) then
            bz(D_DECL(i,j,k))=zero
           else if (maskZM.eq.one) then
            ! do nothing
           else
            print *,"maskZM invalid"
            stop
           endif
#endif     
          else
           print *,"mask_cen invalid"
           stop
          endif
  
         else if ((veldir.gt.0).and.(veldir.lt.nsolve)) then
          ! do not change any coefficients
         else
          print *,"veldir invalid"
          stop
         endif

        else
         print *,"level invalid"
         stop
        endif
       enddo ! i
       enddo ! j
       enddo ! k

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        ioff=MOD(i+j+k+4,2)
        mask(D_DECL(i,j,k))=1.0-ioff
        diag_non_sing(D_DECL(i,j,k))=one  ! avoid divide by zero
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

        local_diag_dual=adual(D_DECL(i,j,k))+offdiagsum
        diag_dual(D_DECL(i,j,k))=local_diag_dual

        local_diag_nodual=a(D_DECL(i,j,k))+offdiagsum
        diag_nodual(D_DECL(i,j,k))=local_diag_nodual

        test_mask=solvemask(D_DECL(i,j,k))

        if (offdiagsum.lt.zero) then
         print *,"offdiagsum invalid"
         stop
        else if (offdiagsum.eq.zero) then
         dtau=dtau_factor
        else if (offdiagsum.gt.zero) then
         dtau=offdiagsum*dtau_factor
        else
         print *,"offdiagsum bust"
         stop
        endif
        if (adual(D_DECL(i,j,k)).ge.zero) then
         ! do nothing
        else
         print *,"adual should be nonneg"
         stop
        endif
        if (adual(D_DECL(i,j,k)).gt.dtau) then
         dtau=zero
        endif
        if (a(D_DECL(i,j,k)).ge.zero) then
         ! do nothing
        else
         print *,"a should be nonneg"
         stop
        endif
        if (a(D_DECL(i,j,k)).gt.dtau) then
         dtau=zero
        endif
        if (adual(D_DECL(i,j,k)).ge.a(D_DECL(i,j,k))) then
         ! do nothing
        else
         print *,"must have adual>=a"
         stop
        endif
        diag_non_sing(D_DECL(i,j,k))=diag_dual(D_DECL(i,j,k))+dtau

        if (local_diag_nodual.eq.zero) then ! local_diag_nodual=a+offdiagsum
         if ((test_mask.eq.one).and.(veldir.gt.0)) then
          print *,"solvemask not uniform wrt veldir"
          stop
         else if ((test_mask.eq.zero).or.(veldir.eq.0)) then
          ! do nothing
         else
          print *,"test_mask or veldir invalid"
          stop
         endif
         solvemask(D_DECL(i,j,k))=zero
         if (check_for_singular.eq.1) then
          ! do nothing
         else
          print *,"check_for_singular invalid"
          stop
         endif
        else if (local_diag_nodual.gt.zero) then
         if (test_mask.ne.one) then
          print *,"test_mask invalid"
          stop
         endif
        else
         print *,"local_diag_nodual invalid"
         stop
        endif
        if (adual(D_DECL(i,j,k)).ge.zero) then
         ! do nothing
        else
         print *,"adual invalid: ",adual(D_DECL(i,j,k))
         stop
        endif
        if (a(D_DECL(i,j,k)).ge.zero) then
         ! do nothing
        else
         print *,"a invalid: ",a(D_DECL(i,j,k))
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
        local_diag_nodual=diag_nodual(D_DECL(i,j,k))
        local_diag_dual=diag_dual(D_DECL(i,j,k))
        if (local_diag_dual.ge.local_diag_nodual) then
         ! do nothing
        else
         print *,"local_diag_dual or local_diag_nodual invalid"
         stop
        endif

        if (local_diag_nodual.eq.zero) then

         if (test_mask.eq.zero) then
          ! do nothing
         else
          print *,"test_mask invalid"
          stop
         endif

         if (local_diag_dual.eq.zero) then
          if (check_for_singular.eq.1) then
           local_diag_dual=local_diag_dual+ &
            diag_regularization*offdiag_coeff
          else
           print *,"local_diag_dual or check_for_singular invalid"
           stop
          endif
         else if (local_diag_dual.gt.zero) then
          ! do nothing
         else
          print *,"local_diag_dual invalid"
          stop
         endif
        else if (local_diag_nodual.gt.zero) then
         if (test_mask.eq.one) then
          ! do nothing
         else
          print *,"test_mask invalid"
          stop
         endif
        else
         print *,"local_diag_nodual invalid"
         stop
        endif
        DD=diag_non_sing(D_DECL(i,j,k))
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
        local_diag_dual=diag_dual(D_DECL(i,j,k))
        local_diag_nodual=diag_nodual(D_DECL(i,j,k))
        if (local_diag_dual.ge.local_diag_nodual) then
         ! do nothing
        else
         print *,"local_diag_dual or local_diag_nodual invalid"
         stop
        endif

        if (local_diag_nodual.eq.zero) then

         if (test_mask.eq.zero) then
          ! do nothing
         else
          print *,"test_mask invalid"
          stop
         endif

         if (local_diag_dual.eq.zero) then
          if (check_for_singular.eq.1) then
           local_diag_dual=local_diag_dual+ &
            diag_regularization*offdiag_coeff
          else
           print *,"local_diag_dual or check_for_singular invalid"
           stop
          endif
         else if (local_diag_dual.gt.zero) then
          ! do nothing
         else
          print *,"local_diag_dual invalid"
          stop
         endif
        else if (local_diag_nodual.gt.zero) then
         if (test_mask.eq.one) then
          ! do nothing
         else
          print *,"test_mask invalid"
          stop
         endif
        else
         print *,"local_diag_nodual invalid"
         stop
        endif

        DD=diag_non_sing(D_DECL(i,j,k))

        if (i.gt.tilelo(1)) then

         local_diag_dual=diag_dual(D_DECL(i-1,j,k))
         if (local_diag_dual.eq.zero) then
          if (check_for_singular.eq.1) then
           local_diag_dual=local_diag_dual+ &
             diag_regularization*offdiag_coeff
          else
           print *,"local_diag_dual or check_for_singular invalid"
           stop
          endif
         else if (local_diag_dual.gt.zero) then
          ! do nothing
         else
          print *,"local_diag_dual invalid"
          stop
         endif
         local_diag_NONSING=diag_non_sing(D_DECL(i-1,j,k))

         DD=DD-bxleft(D_DECL(i,j,k))* &
          bxleft(D_DECL(i,j,k))/local_diag_NONSING

        endif

        if (j.gt.tilelo(2)) then

         local_diag_dual=diag_dual(D_DECL(i,j-1,k))
         if (local_diag_dual.eq.zero) then
          if (check_for_singular.eq.1) then
           local_diag_dual=local_diag_dual+ &
             diag_regularization*offdiag_coeff
          else
           print *,"local_diag_dual or check_for_singular invalid"
           stop
          endif
         else if (local_diag_dual.gt.zero) then
          ! do nothing
         else
          print *,"local_diag_dual invalid"
          stop
         endif
         local_diag_NONSING=diag_non_sing(D_DECL(i,j-1,k))

         DD=DD-byleft(D_DECL(i,j,k))* &
          byleft(D_DECL(i,j,k))/local_diag_NONSING

        endif
#if (AMREX_SPACEDIM==3)
        if (k.gt.tilelo(AMREX_SPACEDIM)) then

         local_diag_dual=diag_dual(D_DECL(i,j,k-1))
         if (local_diag_dual.eq.zero) then
          if (check_for_singular.eq.1) then
           local_diag_dual=local_diag_dual+ &
             diag_regularization*offdiag_coeff
          else
           print *,"local_diag_dual or check_for_singular invalid"
           stop
          endif
         else if (local_diag_dual.gt.zero) then
          ! do nothing
         else
          print *,"local_diag_dual invalid"
          stop
         endif
         local_diag_NONSING=diag_non_sing(D_DECL(i,j,k-1))

         DD=DD-bzleft(D_DECL(i,j,k))* &
          bzleft(D_DECL(i,j,k))/local_diag_NONSING

        endif
#endif
        if (i.lt.tilehi(1)) then

         local_diag_dual=diag_dual(D_DECL(i+1,j,k))
         if (local_diag_dual.eq.zero) then
          if (check_for_singular.eq.1) then
           local_diag_dual=local_diag_dual+ &
             diag_regularization*offdiag_coeff
          else
           print *,"local_diag_dual or check_for_singular invalid"
           stop
          endif
         else if (local_diag_dual.gt.zero) then
          ! do nothing
         else
          print *,"local_diag_dual invalid"
          stop
         endif

         local_diag_NONSING=diag_non_sing(D_DECL(i+1,j,k))

         DD=DD-bxright(D_DECL(i,j,k))* &
          bxright(D_DECL(i,j,k))/local_diag_NONSING 

        endif
        if (j.lt.tilehi(2)) then

         local_diag_dual=diag_dual(D_DECL(i,j+1,k))
         if (local_diag_dual.eq.zero) then
          if (check_for_singular.eq.1) then
           local_diag_dual=local_diag_dual+ &
            diag_regularization*offdiag_coeff
          else
           print *,"local_diag_dual or check_for_singular invalid"
           stop
          endif
         else if (local_diag_dual.gt.zero) then
          ! do nothing
         else
          print *,"local_diag_dual invalid"
          stop
         endif
         local_diag_NONSING=diag_non_sing(D_DECL(i,j+1,k))

         DD=DD-byright(D_DECL(i,j,k))* &
          byright(D_DECL(i,j,k))/local_diag_NONSING

        endif
#if (AMREX_SPACEDIM==3)
        if (k.lt.tilehi(AMREX_SPACEDIM)) then

         local_diag_dual=diag_dual(D_DECL(i,j,k+1))
         if (local_diag_dual.eq.zero) then
          if (check_for_singular.eq.1) then
           local_diag_dual=local_diag_dual+ &
             diag_regularization*offdiag_coeff
          else
           print *,"local_diag_dual or check_for_singular invalid"
           stop
          endif
         else if (local_diag_dual.gt.zero) then
          ! do nothing
         else
          print *,"local_diag_dual invalid"
          stop
         endif
         local_diag_NONSING=diag_non_sing(D_DECL(i,j,k+1))

         DD=DD-bzright(D_DECL(i,j,k))* &
          bzright(D_DECL(i,j,k))/local_diag_NONSING

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
      end subroutine FORT_BUILDMAT


      subroutine FORT_RESIDL ( &
       level, &
       mg_coarsest_level, &
       nsolve, &
       masksing, &
       DIMS(masksing),  &
       res,DIMS(res),  &
       rhs,DIMS(rhs), &
       phi,DIMS(phi), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top)
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
      REAL_T, intent(in) :: masksing(DIMV(masksing))
      INTEGER_T, intent(in) :: DIMDEC(phi)
      REAL_T, intent(in) :: phi(DIMV(phi),nsolve)
      INTEGER_T, intent(in) :: DIMDEC(rhs)
      REAL_T, intent(in) :: rhs(DIMV(rhs),nsolve)
      INTEGER_T, intent(in) :: DIMDEC(res)
      REAL_T, intent(out) :: res(DIMV(res),nsolve)
!
      INTEGER_T i,j,k,veldir
      REAL_T test_mask
!
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
      call checkbound(fablo,fabhi,DIMS(masksing),1,-1,81)
      call checkbound(fablo,fabhi,DIMS(rhs),0,-1,81)
      call checkbound(fablo,fabhi,DIMS(res),0,-1,84)
      call checkbound(fablo,fabhi,DIMS(phi),0,-1,85)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       test_mask=masksing(D_DECL(i,j,k))
       do veldir=1,nsolve
        if (test_mask.eq.zero) then
         res(D_DECL(i,j,k),veldir)=zero
        else if (test_mask.eq.one) then
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
      end subroutine FORT_RESIDL

      subroutine FORT_AVERAGEEC ( &
       nsolve, &
       c,  &
       DIMS(c), &
       f,  &
       DIMS(f), &
       lo, hi, &
       cdir,avg, &
       bfact_coarse,bfact_fine,bfact_top)
      use global_utility_module
      IMPLICIT NONE
!
      INTEGER_T nsolve
      INTEGER_T bfact_coarse,bfact_fine,bfact_top
      INTEGER_T avg
      INTEGER_T lo(AMREX_SPACEDIM)
      INTEGER_T hi(AMREX_SPACEDIM)
      INTEGER_T growlo(3)
      INTEGER_T growhi(3)
      INTEGER_T cdir
      INTEGER_T DIMDEC(f)
      REAL_T f(DIMV(f),nsolve)
      INTEGER_T DIMDEC(c)
      REAL_T c(DIMV(c),nsolve)
!     
      INTEGER_T i,j,k,veldir
      REAL_T denom
!     
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
      call checkbound(lo,hi, &
      DIMS(c) &
      ,0,cdir,201)

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
       call growntileboxMAC(lo,hi,lo,hi,growlo,growhi,0,cdir) 
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
       call growntileboxMAC(lo,hi,lo,hi,growlo,growhi,0,cdir) 
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
       call growntileboxMAC(lo,hi,lo,hi,growlo,growhi,0,cdir) 
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
      end subroutine FORT_AVERAGEEC
!-----------------------------------------------------------------------

      subroutine FORT_AVERAGECC ( &
       nsolve, &
       ncomp_expect, &
       c, &
       DIMS(c), &
       f,  &
       DIMS(f), &
       lo, hi, avg, &
       ngrow, &
       bfact_coarse,bfact_fine,bfact_top)
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
      REAL_T, intent(in) :: f(DIMV(f),ncomp_expect)
      REAL_T, intent(out) :: c(DIMV(c),ncomp_expect)
!
      INTEGER_T i,j,k,veldir
      REAL_T denom
      INTEGER_T sum_mask,max_sum
!
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
      if ((avg.eq.0).or.(avg.eq.1)) then
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
      else if (avg.eq.2) then
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

      call checkbound(lo,hi, &
      DIMS(c), &
      ngrow,-1,301)

      if (avg.eq.1) then
       if (AMREX_SPACEDIM.eq.3) then
        denom=0.125d0
       else if (AMREX_SPACEDIM.eq.2) then
        denom=0.25d0
       else
        print *,"dimension bust"
        stop
       endif
      else if (avg.eq.0) then
       denom=one
      else if (avg.eq.2) then
       denom=one
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
         c(D_DECL(i,j,k),veldir) =  denom*( &
          + f(D_DECL(2*i+1,2*j+1,2*k),veldir) &
          + f(D_DECL(2*i,2*j+1,2*k),veldir) &
          + f(D_DECL(2*i+1,2*j,2*k),veldir) &
          + f(D_DECL(2*i,2*j,2*k),veldir) &
          + f(D_DECL(2*i+1,2*j+1,2*k+1),veldir) &
          + f(D_DECL(2*i,2*j+1,2*k+1),veldir) &
          + f(D_DECL(2*i+1,2*j,2*k+1),veldir) &
          + f(D_DECL(2*i,2*j,2*k+1),veldir) )
        else if (AMREX_SPACEDIM.eq.2) then
         c(D_DECL(i,j,k),veldir) =  denom*( &
          + f(D_DECL(2*i+1,2*j+1,2*k),veldir) &
          + f(D_DECL(2*i,2*j+1,2*k),veldir) &
          + f(D_DECL(2*i+1,2*j,2*k),veldir) &
          + f(D_DECL(2*i,2*j,2*k),veldir))
        else
         print *,"dimension bust"
         stop
        endif
        if (avg.eq.2) then
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

      end subroutine FORT_AVERAGECC



!-----------------------------------------------------------------------


        ! NO TILING 
      subroutine FORT_APPLYBC ( &
       nsolve, &
       phi,  &
       DIMS(phi), &
       bfab,  &
       DIMS(bfab), &
       mfab,  &
       DIMS(mfab), &
       bcpres, &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top)
      use global_utility_module
      IMPLICIT NONE
!
!     mask=0 => fine/fine grid interface....
!
      INTEGER_T nsolve
      INTEGER_T dir,side
      INTEGER_T tilelo(AMREX_SPACEDIM)
      INTEGER_T tilehi(AMREX_SPACEDIM)
      INTEGER_T fablo(AMREX_SPACEDIM)
      INTEGER_T fabhi(AMREX_SPACEDIM)
      INTEGER_T growlo(3)
      INTEGER_T growhi(3)
      INTEGER_T bfact,bfact_top
      INTEGER_T DIMDEC(phi)
      INTEGER_T DIMDEC(bfab)
      INTEGER_T DIMDEC(mfab)
      INTEGER_T bcpres(AMREX_SPACEDIM,2,nsolve)
      REAL_T phi(DIMV(phi),nsolve)
      REAL_T bfab(DIMV(bfab),nsolve)
      REAL_T mfab(DIMV(mfab))
!
      INTEGER_T i,j,k
      INTEGER_T ib,jb,kb
      INTEGER_T iout,jout,kout
      INTEGER_T ireflect,jreflect,kreflect
      INTEGER_T ii,jj,kk
      INTEGER_T bct,veldir,dir_bc
  
      call checkbound(fablo,fabhi, &
       DIMS(phi), &
       1,-1,113)
      call checkbound(fablo,fabhi, &
      DIMS(bfab), &
       1,-1,113)
      call checkbound(fablo,fabhi, &
       DIMS(mfab), &
       1,-1,113)

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
      end subroutine FORT_APPLYBC
