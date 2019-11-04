#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "ABec_F.H"


      subroutine FORT_GSRB( &
       isweep, &
       num_sweeps, & 
       offdiag_coeff, &
       check_for_singular, &
       diag_regularization, &
       phi,DIMS(phi), &
       rhs,DIMS(rhs), &
       diagnonsing, &
       DIMS(diagnonsing), &
       diagsing, &
       bxleft,bxright, &
       byleft,byright, &
       bzleft,bzright, &
       icbx,icby,icbz,icdiag,icdiagrb,mask,ax,solnsave, &
       rhssave,redsoln,blacksoln, &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top, &
       smooth_type)
      use global_utility_module

      IMPLICIT NONE
      INTEGER_T isweep
      INTEGER_T num_sweeps
      INTEGER_T check_for_singular
      REAL_T offdiag_coeff
      REAL_T diag_regularization
      INTEGER_T tilelo(AMREX_SPACEDIM), tilehi(AMREX_SPACEDIM)
      INTEGER_T fablo(AMREX_SPACEDIM), fabhi(AMREX_SPACEDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact,bfact_top
      INTEGER_T DIMDEC(phi)
      INTEGER_T DIMDEC(rhs)
      INTEGER_T DIMDEC(diagnonsing)

      REAL_T  phi(DIMV(phi))
      REAL_T  rhs(DIMV(rhs))
      REAL_T  diagnonsing(DIMV(diagnonsing))
      REAL_T  diagsing(DIMV(diagnonsing))
      REAL_T bxleft(DIMV(diagnonsing))
      REAL_T bxright(DIMV(diagnonsing))
      REAL_T byleft(DIMV(diagnonsing))
      REAL_T byright(DIMV(diagnonsing))
      REAL_T bzleft(DIMV(diagnonsing))
      REAL_T bzright(DIMV(diagnonsing))
      REAL_T icbx(DIMV(diagnonsing))
      REAL_T icby(DIMV(diagnonsing))
      REAL_T icbz(DIMV(diagnonsing))
      REAL_T icdiag(DIMV(diagnonsing))
      REAL_T icdiagrb(DIMV(diagnonsing))
      REAL_T mask(DIMV(diagnonsing))
      REAL_T ax(DIMV(diagnonsing))
      REAL_T solnsave(DIMV(diagnonsing))
      REAL_T rhssave(DIMV(diagnonsing))
      REAL_T redsoln(DIMV(diagnonsing))
      REAL_T blacksoln(DIMV(diagnonsing))

      INTEGER_T i,j,k,smooth_type
      REAL_T XX,YY
      REAL_T local_diag

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
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
      if (diag_regularization.gt.zero) then
       ! do nothing
      else
       print *,"diag_regularization invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(phi),1,-1,18)
      call checkbound(fablo,fabhi,DIMS(diagnonsing),1,-1,19)
      call checkbound(fablo,fabhi,DIMS(rhs),0,-1,23)
      if (num_sweeps.le.1) then
       print *,"num_sweeps invalid"
       stop
      endif

      if (isweep.eq.0) then
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        rhssave(D_DECL(i,j,k))=rhs(D_DECL(i,j,k))+ &
         bxleft(D_DECL(i,j,k))*phi(D_DECL(i-1,j,k))+ &
         bxright(D_DECL(i,j,k))*phi(D_DECL(i+1,j,k))+ &
         byleft(D_DECL(i,j,k))*phi(D_DECL(i,j-1,k))+ &
         byright(D_DECL(i,j,k))*phi(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
        +bzleft(D_DECL(i,j,k))*phi(D_DECL(i,j,k-1))+ &
         bzright(i,j,k)*phi(D_DECL(i,j,k+1)) &
#endif
        -diagsing(D_DECL(i,j,k))*phi(D_DECL(i,j,k))
       enddo
       enddo
       enddo

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        solnsave(D_DECL(i,j,k))=0.0
       enddo
       enddo
       enddo
      else if ((isweep.ge.1).and.(isweep.lt.num_sweeps-1)) then

        ! GSRB
       if (smooth_type.eq.0) then

        if (isweep.eq.1) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          ax(D_DECL(i,j,k))=rhssave(D_DECL(i,j,k))+ &
           bxleft(D_DECL(i,j,k))*solnsave(D_DECL(i-1,j,k))+ &
           bxright(D_DECL(i,j,k))*solnsave(D_DECL(i+1,j,k))+ &
           byleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j-1,k))+ &
           byright(D_DECL(i,j,k))*solnsave(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k-1))+ &
           bzright(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k+1))  &
#endif
          -diagsing(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k))
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          redsoln(D_DECL(i,j,k))=zero
          blacksoln(D_DECL(i,j,k))=zero
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          local_diag=diagsing(D_DECL(i,j,k))
          if (local_diag.eq.zero) then
           if (check_for_singular.eq.1) then
            local_diag=local_diag+diag_regularization*offdiag_coeff
           else
            print *,"local_diag or check_for_singular invalid"
            stop
           endif
          else if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid"
           stop
          endif

          local_diag=diagnonsing(D_DECL(i,j,k))

          redsoln(D_DECL(i,j,k))=ax(D_DECL(i,j,k))/local_diag
         enddo
         enddo
         enddo
        else if (isweep.eq.2) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          local_diag=diagsing(D_DECL(i,j,k))
          if (local_diag.eq.zero) then
           if (check_for_singular.eq.1) then
            local_diag=local_diag+diag_regularization*offdiag_coeff
           else
            print *,"local_diag or check_for_singular invalid"
            stop
           endif
          else if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid"
           stop
          endif

          local_diag=diagnonsing(D_DECL(i,j,k))

          blacksoln(D_DECL(i,j,k))=(ax(D_DECL(i,j,k))+ &
           bxleft(D_DECL(i,j,k))*redsoln(D_DECL(i-1,j,k))+ &
           bxright(D_DECL(i,j,k))*redsoln(D_DECL(i+1,j,k))+ &
           byleft(D_DECL(i,j,k))*redsoln(D_DECL(i,j-1,k))+ &
           byright(D_DECL(i,j,k))*redsoln(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))*redsoln(D_DECL(i,j,k-1))+ &
           bzright(D_DECL(i,j,k))*redsoln(D_DECL(i,j,k+1)) &
#endif
          )/local_diag
         enddo
         enddo
         enddo
        else if (isweep.eq.3) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          local_diag=diagsing(D_DECL(i,j,k))
          if (local_diag.eq.zero) then
           if (check_for_singular.eq.1) then
            local_diag=local_diag+diag_regularization*offdiag_coeff
           else
            print *,"local_diag or check_for_singular invalid"
            stop
           endif
          else if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid"
           stop
          endif

          local_diag=diagnonsing(D_DECL(i,j,k))

          redsoln(D_DECL(i,j,k))=(ax(D_DECL(i,j,k))+ &
           bxleft(D_DECL(i,j,k))*blacksoln(D_DECL(i-1,j,k))+ &
           bxright(D_DECL(i,j,k))*blacksoln(D_DECL(i+1,j,k))+ &
           byleft(D_DECL(i,j,k))*blacksoln(D_DECL(i,j-1,k))+ &
           byright(D_DECL(i,j,k))*blacksoln(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))*blacksoln(D_DECL(i,j,k-1))+ &
           bzright(D_DECL(i,j,k))*blacksoln(D_DECL(i,j,k+1)) &
#endif
          )/local_diag
         enddo
         enddo
         enddo
        else if (isweep.eq.4) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          solnsave(D_DECL(i,j,k))=solnsave(D_DECL(i,j,k))+ &
           mask(D_DECL(i,j,k))*redsoln(D_DECL(i,j,k))+ &
           (one-mask(D_DECL(i,j,k)))*blacksoln(D_DECL(i,j,k))
         enddo
         enddo
         enddo
        else
         print *,"isweep invalid"
         stop
        endif

         ! Jacobi
       else if (smooth_type.eq.3) then

        ! AX=B-AX=B-(D-L-U)X=B+(L+U)X-DX
        if (isweep.eq.1) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          ax(D_DECL(i,j,k))=rhssave(D_DECL(i,j,k))+ &
           bxleft(D_DECL(i,j,k))*solnsave(D_DECL(i-1,j,k))+ &
           bxright(D_DECL(i,j,k))*solnsave(D_DECL(i+1,j,k))+ &
           byleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j-1,k))+ &
           byright(D_DECL(i,j,k))*solnsave(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k-1))+ &
           bzright(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k+1))  &
#endif
          -diagsing(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k))
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          redsoln(D_DECL(i,j,k))=0.0
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          ! R=AX/D=(B+(L+U)X)/D-X

          local_diag=diagsing(D_DECL(i,j,k))
          if (local_diag.eq.zero) then
           if (check_for_singular.eq.1) then
            local_diag=local_diag+diag_regularization*offdiag_coeff
           else
            print *,"local_diag or check_for_singular invalid"
            stop
           endif
          else if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid"
           stop
          endif

          local_diag=diagnonsing(D_DECL(i,j,k))

          redsoln(D_DECL(i,j,k))=ax(D_DECL(i,j,k))/local_diag
         enddo
         enddo
         enddo
        else if (isweep.eq.2) then
         ! X=X+(B+(L+U)X)/D-X=(B+(L+U)X)/D  (A=D-L-U)
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          solnsave(D_DECL(i,j,k))=solnsave(D_DECL(i,j,k))+ &
           redsoln(D_DECL(i,j,k))
         enddo
         enddo
         enddo
        else
         print *,"isweep invalid"
         stop
        endif

         ! ILU
       else if (smooth_type.eq.2) then

        if (isweep.eq.1) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          ax(D_DECL(i,j,k))=rhssave(D_DECL(i,j,k))+ &
           bxleft(D_DECL(i,j,k))*solnsave(D_DECL(i-1,j,k))+ &
           bxright(D_DECL(i,j,k))*solnsave(D_DECL(i+1,j,k))+ &
           byleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j-1,k))+ &
           byright(D_DECL(i,j,k))*solnsave(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k-1))+ &
           bzright(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k+1)) &
#endif
          -diagsing(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k))
         enddo
         enddo
         enddo
        else if (isweep.eq.2) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          YY=ax(D_DECL(i,j,k))
          if (i.gt.tilelo(1)) then
           YY=YY-icbx(D_DECL(i,j,k))*redsoln(D_DECL(i-1,j,k))
          endif
          if (j.gt.tilelo(2)) then
           YY=YY-icby(D_DECL(i,j,k))*redsoln(D_DECL(i,j-1,k))
          endif
#if (AMREX_SPACEDIM==3)
          if (k.gt.tilelo(AMREX_SPACEDIM)) then
           YY=YY-icbz(D_DECL(i,j,k))*redsoln(D_DECL(i,j,k-1))
          endif
#endif
          redsoln(D_DECL(i,j,k))=YY
         enddo
         enddo
         enddo
        else if (isweep.eq.3) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          blacksoln(D_DECL(i,j,k))= &
           redsoln(D_DECL(i,j,k))/icdiag(D_DECL(i,j,k))
         enddo 
         enddo 
         enddo 
        else if (isweep.eq.4) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growhi(3),growlo(3),-1
         do j=growhi(2),growlo(2),-1
         do i=growhi(1),growlo(1),-1
          XX=blacksoln(D_DECL(i,j,k))
          if (i.lt.tilehi(1)) then
           XX=XX-icbx(D_DECL(i+1,j,k))*redsoln(D_DECL(i+1,j,k))
          endif
          if (j.lt.tilehi(2)) then
           XX=XX-icby(D_DECL(i,j+1,k))*redsoln(D_DECL(i,j+1,k))
          endif
#if (AMREX_SPACEDIM==3)
          if (k.lt.tilehi(AMREX_SPACEDIM)) then
           XX=XX-icbz(D_DECL(i,j,k+1))*redsoln(D_DECL(i,j,k+1))
          endif
#endif
          redsoln(D_DECL(i,j,k))=XX
         enddo
         enddo
         enddo
        else if (isweep.eq.5) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          solnsave(D_DECL(i,j,k))= &
           solnsave(D_DECL(i,j,k))+redsoln(D_DECL(i,j,k))
         enddo
         enddo
         enddo
        else 
         print *,"isweep invalid"
         stop
        endif

         ! ICRB
       else if (smooth_type.eq.1) then

        if (isweep.eq.1) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          ax(D_DECL(i,j,k))=rhssave(D_DECL(i,j,k))+ &
           bxleft(D_DECL(i,j,k))*solnsave(D_DECL(i-1,j,k))+ &
           bxright(D_DECL(i,j,k))*solnsave(D_DECL(i+1,j,k))+ &
           byleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j-1,k))+ &
           byright(D_DECL(i,j,k))*solnsave(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k-1))+ &
           bzright(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k+1)) &
#endif
          -diagsing(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k))
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          redsoln(D_DECL(i,j,k))=zero
          blacksoln(D_DECL(i,j,k))=zero
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          local_diag=diagsing(D_DECL(i,j,k))
          if (local_diag.eq.zero) then
           if (check_for_singular.eq.1) then
            local_diag=local_diag+diag_regularization*offdiag_coeff
           else
            print *,"local_diag or check_for_singular invalid"
            stop
           endif
          else if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid"
           stop
          endif

          local_diag=diagnonsing(D_DECL(i,j,k))

          redsoln(D_DECL(i,j,k))=ax(D_DECL(i,j,k))/local_diag
         enddo
         enddo
         enddo
        else if (isweep.eq.2) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          blacksoln(D_DECL(i,j,k))=(ax(D_DECL(i,j,k))+ &
           bxleft(D_DECL(i,j,k))*redsoln(D_DECL(i-1,j,k))+ &
           bxright(D_DECL(i,j,k))*redsoln(D_DECL(i+1,j,k))+ &
           byleft(D_DECL(i,j,k))*redsoln(D_DECL(i,j-1,k))+ &
           byright(D_DECL(i,j,k))*redsoln(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))*redsoln(D_DECL(i,j,k-1))+ &
           bzright(D_DECL(i,j,k))*redsoln(D_DECL(i,j,k+1)) &
#endif
          )/icdiagrb(D_DECL(i,j,k))
         enddo
         enddo
         enddo
        else if (isweep.eq.3) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          local_diag=diagsing(D_DECL(i,j,k))
          if (local_diag.eq.zero) then
           if (check_for_singular.eq.1) then
            local_diag=local_diag+diag_regularization*offdiag_coeff
           else
            print *,"local_diag or check_for_singular invalid"
            stop
           endif
          else if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid"
           stop
          endif

          local_diag=diagnonsing(D_DECL(i,j,k))

          redsoln(D_DECL(i,j,k))=(ax(D_DECL(i,j,k))+ &
           bxleft(D_DECL(i,j,k))*blacksoln(D_DECL(i-1,j,k))+ &
           bxright(D_DECL(i,j,k))*blacksoln(D_DECL(i+1,j,k))+ &
           byleft(D_DECL(i,j,k))*blacksoln(D_DECL(i,j-1,k))+ &
           byright(D_DECL(i,j,k))*blacksoln(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
          +bzleft(D_DECL(i,j,k))*blacksoln(D_DECL(i,j,k-1))+ &
           bzright(D_DECL(i,j,k))*blacksoln(D_DECL(i,j,k+1)) &
#endif
          )/local_diag
         enddo
         enddo
         enddo
        else if (isweep.eq.4) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          solnsave(D_DECL(i,j,k))=solnsave(D_DECL(i,j,k))+ &
           mask(D_DECL(i,j,k))*redsoln(D_DECL(i,j,k))+ &
           (one-mask(D_DECL(i,j,k)))*blacksoln(D_DECL(i,j,k))
         enddo
         enddo
         enddo
        else
         print *,"isweep invalid"
         stop
        endif
       else 
        print *,"smooth_type invalid"
        stop
       endif

      else if (isweep.eq.num_sweeps-1) then
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        phi(D_DECL(i,j,k))=phi(D_DECL(i,j,k))+solnsave(D_DECL(i,j,k))
       enddo
       enddo
       enddo
      else
       print *,"isweep invalid"
       stop
      endif

      return
      end subroutine FORT_GSRB

      subroutine FORT_ADOTX( &
       offdiag_coeff, &
       check_for_singular, &
       diag_regularization, &
       y,DIMS(y), &
       x,DIMS(x), &
       diagnonsing, &
       DIMS(diagnonsing), &
       diagsing, &
       bxleft,bxright, &
       byleft,byright, &
       bzleft,bzright, &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top)
      use global_utility_module

      IMPLICIT NONE
      INTEGER_T check_for_singular
      REAL_T offdiag_coeff
      REAL_T diag_regularization
      INTEGER_T tilelo(AMREX_SPACEDIM), tilehi(AMREX_SPACEDIM)
      INTEGER_T fablo(AMREX_SPACEDIM), fabhi(AMREX_SPACEDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact,bfact_top
      INTEGER_T DIMDEC(y)
      INTEGER_T DIMDEC(x)
      INTEGER_T DIMDEC(diagnonsing)

      REAL_T  y(DIMV(y))
      REAL_T  x(DIMV(x))
      REAL_T  diagnonsing(DIMV(diagnonsing))
      REAL_T  diagsing(DIMV(diagnonsing))
      REAL_T bxleft(DIMV(diagnonsing))
      REAL_T bxright(DIMV(diagnonsing))
      REAL_T byleft(DIMV(diagnonsing))
      REAL_T byright(DIMV(diagnonsing))
      REAL_T bzleft(DIMV(diagnonsing))
      REAL_T bzright(DIMV(diagnonsing))

      INTEGER_T i,j,k

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
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
      if (diag_regularization.gt.zero) then
       ! do nothing
      else
       print *,"diag_regularization invalid"
       stop
      endif

      call checkbound(fablo,fabhi, &
      DIMS(y), &
      0,-1,18)
      call checkbound(fablo,fabhi, &
      DIMS(diagnonsing), &
      1,-1,19)
      call checkbound(fablo,fabhi, &
      DIMS(x), &
      1,-1,23)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       y(D_DECL(i,j,k))= &
        diagsing(D_DECL(i,j,k))*x(D_DECL(i,j,k))- &
        bxleft(D_DECL(i,j,k))*x(D_DECL(i-1,j,k))- &
        bxright(D_DECL(i,j,k))*x(D_DECL(i+1,j,k)) &
#if (AMREX_SPACEDIM==3)
        -bzleft(D_DECL(i,j,k))*x(D_DECL(i,j,k-1))- &
        bzright(D_DECL(i,j,k))*x(D_DECL(i,j,k+1)) &
#endif
        -byleft(D_DECL(i,j,k))*x(D_DECL(i,j-1,k))- &
        byright(D_DECL(i,j,k))*x(D_DECL(i,j+1,k)) 
      enddo
      enddo
      enddo

      return
      end subroutine FORT_ADOTX

      subroutine FORT_DIAGSUM( &
       y, &
       DIMS(y), &
       a,  &
       DIMS(a), &
       bX, &
       DIMS(bX), &
       bY, &
       DIMS(bY), &
       bZ, &
       DIMS(bZ), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top)
      use global_utility_module

      IMPLICIT NONE
      INTEGER_T tilelo(AMREX_SPACEDIM), tilehi(AMREX_SPACEDIM)
      INTEGER_T fablo(AMREX_SPACEDIM), fabhi(AMREX_SPACEDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact,bfact_top
      INTEGER_T DIMDEC(y)
      INTEGER_T DIMDEC(a)
      INTEGER_T DIMDEC(bX)
      INTEGER_T DIMDEC(bY)
      INTEGER_T DIMDEC(bZ)

      REAL_T  y(DIMV(y))
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T bZ(DIMV(bZ))

      INTEGER_T i,j,k
      REAL_T bxleft,bxright,byleft,byright,bzleft,bzright
      REAL_T offdiagsum

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      call checkbound(fablo,fabhi, &
      DIMS(y) &
      ,0,-1,18)
      call checkbound(fablo,fabhi, &
      DIMS(a) &
      ,0,-1,19)
      call checkbound(fablo,fabhi, &
      DIMS(bX) &
      ,0,0,20)
      call checkbound(fablo,fabhi, &
      DIMS(bY) &
      ,0,1,21)
      call checkbound(fablo,fabhi, &
      DIMS(bZ) &
      ,0,AMREX_SPACEDIM-1,22)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       bxleft=bX(D_DECL(i,j,k))
       bxright=bX(D_DECL(i+1,j,k))
       byleft=bY(D_DECL(i,j,k))
       byright=bY(D_DECL(i,j+1,k))
       bzleft=bZ(D_DECL(i,j,k))
       bzright=bZ(D_DECL(i,j,k+1))

       offdiagsum=bxleft+bxright &
#if (AMREX_SPACEDIM==3)
        +bzleft+bzright &
#endif
        +byleft+byright
       if (offdiagsum.ne.zero) then
        offdiagsum=one
       endif
       y(D_DECL(i,j,k)) = offdiagsum

      enddo
      enddo
      enddo
      end subroutine FORT_DIAGSUM

