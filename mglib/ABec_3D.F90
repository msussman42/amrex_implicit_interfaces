#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>

#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "AMReX_CONSTANTS.H"
#include "ABec_F.H"

      module cpp_abec
      use amrex_fort_module, only : amrex_real
      
      contains

      subroutine fort_gsrb( &
       level, &
       mg_coarsest_level, &
       isweep, &
       num_sweeps, & 
       masksing, &
       DIMS(masksing), &
       phi,DIMS(phi), &
       rhs,DIMS(rhs), &
       diagfab, &
       DIMS(diagfab), &
       bxleft,bxright, &
       byleft,byright, &
       bzleft,bzright, &
       icbx,icby,icbz,icdiag,icdiagrb, &
       mask,ax,solnsave, &
       rhssave,redsoln,blacksoln, &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top, &
       smooth_type) &
      bind(c,name='fort_gsrb')

      use global_utility_module

      IMPLICIT NONE

      integer, intent(in) :: level
      integer, intent(in) :: mg_coarsest_level
      integer, intent(in) :: isweep
      integer, intent(in) :: num_sweeps
      integer, intent(in) :: tilelo(AMREX_SPACEDIM)
      integer, intent(in) :: tilehi(AMREX_SPACEDIM)
      integer, intent(in) :: fablo(AMREX_SPACEDIM)
      integer, intent(in) :: fabhi(AMREX_SPACEDIM)
      integer :: growlo(3), growhi(3)
      integer, intent(in) :: bfact,bfact_top
      integer, intent(in) :: DIMDEC(masksing)
      integer, intent(in) :: DIMDEC(phi)
      integer, intent(in) :: DIMDEC(rhs)
      integer, intent(in) :: DIMDEC(diagfab)

      real(amrex_real), intent(in), target :: masksing(DIMV(masksing))

      real(amrex_real), intent(inout), target :: phi(DIMV(phi))
      real(amrex_real), pointer :: phi_ptr(D_DECL(:,:,:))

      real(amrex_real), intent(in), target :: rhs(DIMV(rhs))
      real(amrex_real), intent(in), target :: diagfab(DIMV(diagfab))
      real(amrex_real), intent(in), target :: bxleft(DIMV(diagfab))
      real(amrex_real), intent(in), target :: bxright(DIMV(diagfab))
      real(amrex_real), intent(in), target :: byleft(DIMV(diagfab))
      real(amrex_real), intent(in), target :: byright(DIMV(diagfab))
      real(amrex_real), intent(in), target :: bzleft(DIMV(diagfab))
      real(amrex_real), intent(in), target :: bzright(DIMV(diagfab))
      real(amrex_real), intent(in), target :: icbx(DIMV(diagfab))
      real(amrex_real), intent(in), target :: icby(DIMV(diagfab))
      real(amrex_real), intent(in), target :: icbz(DIMV(diagfab))
      real(amrex_real), intent(in), target :: icdiag(DIMV(diagfab))
      real(amrex_real), intent(in), target :: icdiagrb(DIMV(diagfab))
      real(amrex_real), intent(in), target :: mask(DIMV(diagfab))
      real(amrex_real), intent(out), target :: ax(DIMV(diagfab))
      real(amrex_real), intent(inout), target :: solnsave(DIMV(diagfab))
      real(amrex_real), intent(inout), target :: rhssave(DIMV(diagfab))
      real(amrex_real), intent(inout), target :: redsoln(DIMV(diagfab))
      real(amrex_real), intent(inout), target :: blacksoln(DIMV(diagfab))
      integer, intent(in) :: smooth_type

      integer i,j,k
      real(amrex_real) local_diag
      real(amrex_real) test_mask
      real(amrex_real) :: XX,YY

      phi_ptr=>phi

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      if ((level.ge.0).and.(level.le.mg_coarsest_level)) then
       ! do nothing
      else
       print *,"level or mg_coarsest_level invalid"
       stop
      endif

      call checkbound_array1(fablo,fabhi,masksing,1,-1)
      call checkbound_array1(fablo,fabhi,phi_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,diagfab,1,-1)
      call checkbound_array1(fablo,fabhi,rhs,0,-1)

      if (num_sweeps.le.1) then
       print *,"num_sweeps invalid"
       stop
      endif

      if (isweep.eq.0) then
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
         !rhssave=b-Ax^{k}=b-(D-L-U)x^{k}
        if (diagfab(D_DECL(i,j,k)).gt.zero) then
         rhssave(D_DECL(i,j,k))=rhs(D_DECL(i,j,k))+ &
          bxleft(D_DECL(i,j,k))*phi(D_DECL(i-1,j,k))+ &
          bxright(D_DECL(i,j,k))*phi(D_DECL(i+1,j,k))+ &
          byleft(D_DECL(i,j,k))*phi(D_DECL(i,j-1,k))+ &
          byright(D_DECL(i,j,k))*phi(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
         +bzleft(D_DECL(i,j,k))*phi(D_DECL(i,j,k-1))+ &
          bzright(i,j,k)*phi(D_DECL(i,j,k+1)) &
#endif
         -diagfab(D_DECL(i,j,k))*phi(D_DECL(i,j,k))
        else
         print *,"diagfab invalid"
         stop
        endif
       enddo
       enddo
       enddo

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        solnsave(D_DECL(i,j,k))=0.0
       enddo
       enddo
       enddo
      else if ((isweep.ge.1).and.(isweep.lt.num_sweeps-1)) then

        ! GSRB (Gauss Seidel Red Black)
       if (smooth_type.eq.0) then

         !rhssave=b-Ax^{k} in the bulk after isweep.eq.0
         !solnsave=0 in the bulk after isweep.eq.0
        if (isweep.eq.1) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          if (diagfab(D_DECL(i,j,k)).gt.zero) then
           ax(D_DECL(i,j,k))=rhssave(D_DECL(i,j,k))+ &
            bxleft(D_DECL(i,j,k))*solnsave(D_DECL(i-1,j,k))+ &
            bxright(D_DECL(i,j,k))*solnsave(D_DECL(i+1,j,k))+ &
            byleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j-1,k))+ &
            byright(D_DECL(i,j,k))*solnsave(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
           +bzleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k-1))+ &
            bzright(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k+1))  &
#endif
           -diagfab(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k))
          else
           print *,"diagfab invalid"
           stop
          endif
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          redsoln(D_DECL(i,j,k))=zero
          blacksoln(D_DECL(i,j,k))=zero
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          local_diag=diagfab(D_DECL(i,j,k))
          if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid"
           stop
          endif

          redsoln(D_DECL(i,j,k))=ax(D_DECL(i,j,k))/local_diag
         enddo
         enddo
         enddo
        else if (isweep.eq.2) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          local_diag=diagfab(D_DECL(i,j,k))
          if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid"
           stop
          endif

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
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          local_diag=diagfab(D_DECL(i,j,k))

          if (local_diag.gt.zero) then
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
          else
           print *,"local_diag invalid 3"
           stop
          endif
         enddo
         enddo
         enddo
        else if (isweep.eq.4) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
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

         ! Weighted Jacobi
       else if (smooth_type.eq.1) then

        !rhssave=b-Ax^{k} in the bulk after isweep.eq.0
        !solnsave=0 in the bulk after isweep.eq.0
        ! 0=B-AX=B-(D-L-U)X=B+(L+U)X-DX
        ! X^predict=D^{-1}(B+(L+U)X^{k})=
        ! D^{-1}(B+(D-D+L+U)X^{k})=
        ! D^{-1}(B-AX+DX)=D^{-1}R+X^{k}
        ! weighted Jacobi: X^new=(X^predict+X^{k})/2
        ! R = residual = B - AX
        ! Jacobi method: X^{k+1}=X^{k}+D^{-1}(B-A X^{k})=X^{k}+D^{-1}R^{k}
        ! Weighted Jacobi: 
        !  X^{k+1}=X^{k}+D^{-1}R^{k}/2
        if (isweep.eq.1) then
           ! zero grow cells prescribed, so growlo=tilelo,
           ! growhi=tilehi  an instructive diagram is found here:
           ! https://amrex-codes.github.io/amrex/docs_html/GPU ....
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          if (diagfab(D_DECL(i,j,k)).gt.zero) then
             !ax=B-(D-L-U)x^k
           ax(D_DECL(i,j,k))=rhssave(D_DECL(i,j,k))+ &
            bxleft(D_DECL(i,j,k))*solnsave(D_DECL(i-1,j,k))+ &
            bxright(D_DECL(i,j,k))*solnsave(D_DECL(i+1,j,k))+ &
            byleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j-1,k))+ &
            byright(D_DECL(i,j,k))*solnsave(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
           +bzleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k-1))+ &
            bzright(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k+1))  &
#endif
           -diagfab(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k))
          else
           print *,"diagfab invalid"
           stop
          endif
         enddo
         enddo
         enddo
          ! set redsoln=0 in the interior + 1 ghost level
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          redsoln(D_DECL(i,j,k))=0.0
         enddo
         enddo
         enddo
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          local_diag=diagfab(D_DECL(i,j,k))
          if (local_diag.gt.zero) then
           ! do nothing
          else
           print *,"local_diag invalid: ",local_diag
           stop
          endif

          redsoln(D_DECL(i,j,k))=ax(D_DECL(i,j,k))/local_diag
         enddo
         enddo
         enddo
        else if (isweep.eq.2) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          solnsave(D_DECL(i,j,k))=solnsave(D_DECL(i,j,k))+ &
           half*redsoln(D_DECL(i,j,k))
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
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          if (diagfab(D_DECL(i,j,k)).gt.zero) then
           ax(D_DECL(i,j,k))=rhssave(D_DECL(i,j,k))+ &
            bxleft(D_DECL(i,j,k))*solnsave(D_DECL(i-1,j,k))+ &
            bxright(D_DECL(i,j,k))*solnsave(D_DECL(i+1,j,k))+ &
            byleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j-1,k))+ &
            byright(D_DECL(i,j,k))*solnsave(D_DECL(i,j+1,k)) &
#if (AMREX_SPACEDIM==3)
           +bzleft(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k-1))+ &
            bzright(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k+1)) &
#endif
           -diagfab(D_DECL(i,j,k))*solnsave(D_DECL(i,j,k))
          else
           print *,"diagfab invalid"
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
          local_diag=icdiag(D_DECL(i,j,k))
          if (local_diag.ne.zero) then
           blacksoln(D_DECL(i,j,k))= &
            redsoln(D_DECL(i,j,k))/local_diag
          else
           print *,"local_diag invalid 5"
           stop
          endif
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
         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)
          solnsave(D_DECL(i,j,k))= &
           solnsave(D_DECL(i,j,k))+redsoln(D_DECL(i,j,k))
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
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        test_mask=masksing(D_DECL(i,j,k))
        if (test_mask.eq.one) then
         phi(D_DECL(i,j,k))=phi(D_DECL(i,j,k))+solnsave(D_DECL(i,j,k))
        else
         print *,"test_mask invalid in gsrb"
         stop
        endif
       enddo
       enddo
       enddo
      else
       print *,"isweep invalid"
       stop
      endif

      return
      end subroutine fort_gsrb

      subroutine fort_adotx( &
       level, &
       mg_coarsest_level, &
       masksing, &
       DIMS(masksing), &
       y,DIMS(y), &
       x,DIMS(x), &
       diagfab, &
       DIMS(diagfab), &
       bxleft,bxright, &
       byleft,byright, &
       bzleft,bzright, &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top) &
      bind(c,name='fort_adotx')

      use global_utility_module

      IMPLICIT NONE

      integer, intent(in) :: level
      integer, intent(in) :: mg_coarsest_level
      integer, intent(in) :: tilelo(AMREX_SPACEDIM)
      integer, intent(in) :: tilehi(AMREX_SPACEDIM)
      integer, intent(in) :: fablo(AMREX_SPACEDIM)
      integer, intent(in) :: fabhi(AMREX_SPACEDIM)
      integer :: growlo(3), growhi(3)
      integer, intent(in) :: bfact,bfact_top
      integer, intent(in) :: DIMDEC(masksing)
      integer, intent(in) :: DIMDEC(y)
      integer, intent(in) :: DIMDEC(x)
      integer, intent(in) :: DIMDEC(diagfab)

      real(amrex_real), intent(in), target :: masksing(DIMV(masksing))
      real(amrex_real), intent(out), target :: y(DIMV(y))
      real(amrex_real), pointer :: y_ptr(D_DECL(:,:,:))

      real(amrex_real), intent(in), target :: x(DIMV(x))
      real(amrex_real), intent(in), target :: diagfab(DIMV(diagfab))
      real(amrex_real), intent(in), target :: bxleft(DIMV(diagfab))
      real(amrex_real), intent(in), target :: bxright(DIMV(diagfab))
      real(amrex_real), intent(in), target :: byleft(DIMV(diagfab))
      real(amrex_real), intent(in), target :: byright(DIMV(diagfab))
      real(amrex_real), intent(in), target :: bzleft(DIMV(diagfab))
      real(amrex_real), intent(in), target :: bzright(DIMV(diagfab))

      integer i,j,k
      real(amrex_real) test_mask

      y_ptr=>y

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      if ((level.ge.0).and.(level.le.mg_coarsest_level)) then
       ! do nothing
      else
       print *,"level or mg_coarsest_level invalid"
       stop
      endif

      call checkbound_array1(fablo,fabhi,masksing,1,-1)
      call checkbound_array1(fablo,fabhi,y_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,diagfab,1,-1)
      call checkbound_array1(fablo,fabhi,x,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       test_mask=masksing(D_DECL(i,j,k))
       if (test_mask.eq.one) then
        if (diagfab(D_DECL(i,j,k)).gt.zero) then
         y(D_DECL(i,j,k))= &
          diagfab(D_DECL(i,j,k))*x(D_DECL(i,j,k))- &
          bxleft(D_DECL(i,j,k))*x(D_DECL(i-1,j,k))- &
          bxright(D_DECL(i,j,k))*x(D_DECL(i+1,j,k)) &
#if (AMREX_SPACEDIM==3)
          -bzleft(D_DECL(i,j,k))*x(D_DECL(i,j,k-1))- &
          bzright(D_DECL(i,j,k))*x(D_DECL(i,j,k+1)) &
#endif
          -byleft(D_DECL(i,j,k))*x(D_DECL(i,j-1,k))- &
          byright(D_DECL(i,j,k))*x(D_DECL(i,j+1,k)) 
        else
         print *,"diagfab invalid"
         stop
        endif
       else
        print *,"test_mask invalid in adotx"
        stop
       endif
      enddo
      enddo
      enddo

      return
      end subroutine fort_adotx

      subroutine fort_diagsum( &
       y, &
       DIMS(y), &
       bX, &
       DIMS(bX), &
       bY, &
       DIMS(bY), &
       bZ, &
       DIMS(bZ), &
       tilelo,tilehi, &
       fablo,fabhi,bfact,bfact_top) &
      bind(c,name='fort_diagsum')

      use global_utility_module

      IMPLICIT NONE
      integer, intent(in) :: tilelo(AMREX_SPACEDIM)
      integer, intent(in) :: tilehi(AMREX_SPACEDIM)
      integer, intent(in) :: fablo(AMREX_SPACEDIM)
      integer, intent(in) :: fabhi(AMREX_SPACEDIM)
      integer growlo(3), growhi(3)
      integer, intent(in) :: bfact,bfact_top
      integer, intent(in) :: DIMDEC(y)
      integer, intent(in) :: DIMDEC(bX)
      integer, intent(in) :: DIMDEC(bY)
      integer, intent(in) :: DIMDEC(bZ)

      real(amrex_real), intent(out), target ::  y(DIMV(y))
      real(amrex_real), pointer :: y_ptr(D_DECL(:,:,:))

      real(amrex_real), intent(in), target :: bX(DIMV(bX))
      real(amrex_real), intent(in), target :: bY(DIMV(bY))
      real(amrex_real), intent(in), target :: bZ(DIMV(bZ))

      integer i,j,k
      real(amrex_real) bxleft,bxright,byleft,byright,bzleft,bzright
      real(amrex_real) offdiagsum

      y_ptr=>y

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (bfact_top.lt.1) then
       print *,"bfact_top invalid"
       stop
      endif
      call checkbound_array1(fablo,fabhi,y_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,bX,0,0)
      call checkbound_array1(fablo,fabhi,bY,0,1)
      call checkbound_array1(fablo,fabhi,bZ,0,AMREX_SPACEDIM-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

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
      end subroutine fort_diagsum

      end module cpp_abec
