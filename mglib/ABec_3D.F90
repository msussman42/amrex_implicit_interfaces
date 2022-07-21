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

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: mg_coarsest_level
      INTEGER_T, intent(in) :: isweep
      INTEGER_T, intent(in) :: num_sweeps
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(masksing)
      INTEGER_T, intent(in) :: DIMDEC(phi)
      INTEGER_T, intent(in) :: DIMDEC(rhs)
      INTEGER_T, intent(in) :: DIMDEC(diagfab)

      REAL_T, intent(in), target :: masksing(DIMV(masksing))

      REAL_T, intent(inout), target :: phi(DIMV(phi))
      REAL_T, pointer :: phi_ptr(D_DECL(:,:,:))

      REAL_T, intent(in), target :: rhs(DIMV(rhs))
      REAL_T, intent(in), target :: diagfab(DIMV(diagfab))
      REAL_T, intent(in), target :: bxleft(DIMV(diagfab))
      REAL_T, intent(in), target :: bxright(DIMV(diagfab))
      REAL_T, intent(in), target :: byleft(DIMV(diagfab))
      REAL_T, intent(in), target :: byright(DIMV(diagfab))
      REAL_T, intent(in), target :: bzleft(DIMV(diagfab))
      REAL_T, intent(in), target :: bzright(DIMV(diagfab))
      REAL_T, intent(in), target :: icbx(DIMV(diagfab))
      REAL_T, intent(in), target :: icby(DIMV(diagfab))
      REAL_T, intent(in), target :: icbz(DIMV(diagfab))
      REAL_T, intent(in), target :: icdiag(DIMV(diagfab))
      REAL_T, intent(in), target :: icdiagrb(DIMV(diagfab))
      REAL_T, intent(in), target :: mask(DIMV(diagfab))
      REAL_T, intent(out), target :: ax(DIMV(diagfab))
      REAL_T, intent(inout), target :: solnsave(DIMV(diagfab))
      REAL_T, intent(inout), target :: rhssave(DIMV(diagfab))
      REAL_T, intent(inout), target :: redsoln(DIMV(diagfab))
      REAL_T, intent(inout), target :: blacksoln(DIMV(diagfab))
      INTEGER_T, intent(in) :: smooth_type

      INTEGER_T i,j,k
      REAL_T XX,YY
      REAL_T local_diag
      REAL_T test_mask

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
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
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
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

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
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

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

          local_diag=diagfab(D_DECL(i,j,k))

          if (local_diag.gt.zero) then
           redsoln(D_DECL(i,j,k))=ax(D_DECL(i,j,k))/local_diag
          else
           print *,"local_diag invalid 6"
           stop
          endif
         enddo
         enddo
         enddo
        else if (isweep.eq.2) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)
          local_diag=icdiagrb(D_DECL(i,j,k))
          if (local_diag.ne.zero) then
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
          else
           print *,"local_diag invalid 7"
           stop
          endif
         enddo
         enddo
         enddo
        else if (isweep.eq.3) then
         call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

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
           print *,"local_diag invalid 8"
           stop
          endif
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

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: mg_coarsest_level
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(masksing)
      INTEGER_T, intent(in) :: DIMDEC(y)
      INTEGER_T, intent(in) :: DIMDEC(x)
      INTEGER_T, intent(in) :: DIMDEC(diagfab)

      REAL_T, intent(in), target :: masksing(DIMV(masksing))
      REAL_T, intent(out), target :: y(DIMV(y))
      REAL_T, pointer :: y_ptr(D_DECL(:,:,:))

      REAL_T, intent(in), target :: x(DIMV(x))
      REAL_T, intent(in), target :: diagfab(DIMV(diagfab))
      REAL_T, intent(in), target :: bxleft(DIMV(diagfab))
      REAL_T, intent(in), target :: bxright(DIMV(diagfab))
      REAL_T, intent(in), target :: byleft(DIMV(diagfab))
      REAL_T, intent(in), target :: byright(DIMV(diagfab))
      REAL_T, intent(in), target :: bzleft(DIMV(diagfab))
      REAL_T, intent(in), target :: bzright(DIMV(diagfab))

      INTEGER_T i,j,k
      REAL_T test_mask

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
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

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
      INTEGER_T, intent(in) :: tilelo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: tilehi(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fablo(AMREX_SPACEDIM)
      INTEGER_T, intent(in) :: fabhi(AMREX_SPACEDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_top
      INTEGER_T, intent(in) :: DIMDEC(y)
      INTEGER_T, intent(in) :: DIMDEC(bX)
      INTEGER_T, intent(in) :: DIMDEC(bY)
      INTEGER_T, intent(in) :: DIMDEC(bZ)

      REAL_T, intent(out), target ::  y(DIMV(y))
      REAL_T, pointer :: y_ptr(D_DECL(:,:,:))

      REAL_T, intent(in), target :: bX(DIMV(bX))
      REAL_T, intent(in), target :: bY(DIMV(bY))
      REAL_T, intent(in), target :: bZ(DIMV(bZ))

      INTEGER_T i,j,k
      REAL_T bxleft,bxright,byleft,byright,bzleft,bzright
      REAL_T offdiagsum

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
      end subroutine fort_diagsum

      end module cpp_abec
