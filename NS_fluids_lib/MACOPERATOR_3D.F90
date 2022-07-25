#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif


#include "AMReX_FORT_INTEGER.H"
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "EXTRAP_COMP.H"
#include "MACOPERATOR_F.H"

#define DEBUG_THERMAL_WEIGHT 0

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

      module macoperator_module
      use probf90_module

      contains

! maskcov=1 outside domain and on fine-fine ghost cells
! maskcov=1 for interior cells not covered by a finer cell
! fwtx,fwty,fwtz are not averaged down.
! bx,by,bz are equal to bxcoefnoarea; bxcoefnoarea is not averaged down at this
! point.
      subroutine fort_init_mask_sing( &
       level, &
       finest_level, &
       nsolve, &
       project_option, &
       masksolv,DIMS(masksolv), & ! ONES_MF in c++
       maskcov,DIMS(maskcov), &
       alpha,DIMS(alpha), &
       offdiagcheck, &
       DIMS(offdiagcheck), &
       maskdivres, &
       DIMS(maskdivres), &
       maskres,DIMS(maskres), &
       mdot,DIMS(mdot), &
       bx,DIMS(bx), &
       by,DIMS(by), &
       bz,DIMS(bz), &
       fwtx,DIMS(fwtx), &
       fwty,DIMS(fwty), &
       fwtz,DIMS(fwtz), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       bc) &
      bind(c,name='fort_init_mask_sing')
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: nsolve
      INTEGER_T, INTENT(in) :: project_option
      INTEGER_T, INTENT(in) :: DIMDEC(masksolv) ! ONES_MF in c++
      INTEGER_T, INTENT(in) :: DIMDEC(maskcov)
      INTEGER_T, INTENT(in) :: DIMDEC(alpha)
      INTEGER_T, INTENT(in) :: DIMDEC(offdiagcheck)
      INTEGER_T, INTENT(in) :: DIMDEC(maskdivres)
      INTEGER_T, INTENT(in) :: DIMDEC(maskres)
      INTEGER_T, INTENT(in) :: DIMDEC(mdot)
      INTEGER_T, INTENT(in) :: DIMDEC(bx)
      INTEGER_T, INTENT(in) :: DIMDEC(by)
      INTEGER_T, INTENT(in) :: DIMDEC(bz)
      INTEGER_T, INTENT(in) :: DIMDEC(fwtx)
      INTEGER_T, INTENT(in) :: DIMDEC(fwty)
      INTEGER_T, INTENT(in) :: DIMDEC(fwtz)
      INTEGER_T, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T             :: growlo(3), growhi(3)
      INTEGER_T, INTENT(in) :: bc(SDIM,2,nsolve)

       ! ONES_MF in c++
      REAL_T, INTENT(out), target :: masksolv(DIMV(masksolv))
      REAL_T, pointer :: masksolv_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: maskcov(DIMV(maskcov))
      REAL_T, pointer :: maskcov_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: alpha(DIMV(alpha),nsolve)
      REAL_T, pointer :: alpha_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: offdiagcheck(DIMV(offdiagcheck),nsolve)
      REAL_T, pointer :: offdiagcheck_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out), target :: maskdivres(DIMV(maskdivres))
      REAL_T, pointer :: maskdivres_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(out), target :: maskres(DIMV(maskres))
      REAL_T, pointer :: maskres_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: mdot(DIMV(mdot),nsolve)
      REAL_T, pointer :: mdot_ptr(D_DECL(:,:,:),:)
       ! coeff * areafrac * areaface / (dxfrac*dx)
      REAL_T, INTENT(in), target :: bx(DIMV(bx),nsolve) 
      REAL_T, pointer :: bx_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: by(DIMV(by),nsolve)
      REAL_T, pointer :: by_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: bz(DIMV(bz),nsolve)
      REAL_T, pointer :: bz_ptr(D_DECL(:,:,:),:)
       ! coeff * areafrac / dxfrac
      REAL_T, INTENT(in), target :: fwtx(DIMV(fwtx),nsolve)  
      REAL_T, pointer :: fwtx_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: fwty(DIMV(fwty),nsolve)
      REAL_T, pointer :: fwty_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: fwtz(DIMV(fwtz),nsolve)
      REAL_T, pointer :: fwtz_ptr(D_DECL(:,:,:),:)

      INTEGER_T i,j,k
      INTEGER_T inormal
      INTEGER_T ii,jj,kk
      INTEGER_T iface,jface,kface
      INTEGER_T icell,jcell,kcell
      INTEGER_T dir,side
      INTEGER_T veldir
      REAL_T bface
      REAL_T facewtsum,offdiagsum
      REAL_T local_diag
      REAL_T local_diag_check1
      REAL_T local_diag_check2

      if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
       print *,"nsolve invalid24"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid nsgenerate"
       stop
      endif

       ! ONES_MF in c++
      masksolv_ptr=>masksolv
      maskcov_ptr=>maskcov
      call checkbound_array1(fablo,fabhi,masksolv_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      alpha_ptr=>alpha
      call checkbound_array(fablo,fabhi,alpha_ptr,0,-1)
      offdiagcheck_ptr=>offdiagcheck
      call checkbound_array(fablo,fabhi,offdiagcheck_ptr,0,-1)
      maskdivres_ptr=>maskdivres
      call checkbound_array1(fablo,fabhi,maskdivres_ptr,0,-1)
      maskres_ptr=>maskres
      call checkbound_array1(fablo,fabhi,maskres_ptr,0,-1)
      mdot_ptr=>mdot
      call checkbound_array(fablo,fabhi,mdot_ptr,0,-1)
      bx_ptr=>bx
      by_ptr=>by
      bz_ptr=>bz
      call checkbound_array(fablo,fabhi,bx_ptr,0,0)
      call checkbound_array(fablo,fabhi,by_ptr,0,1)
      call checkbound_array(fablo,fabhi,bz_ptr,0,AMREX_SPACEDIM-1)
      fwtx_ptr=>fwtx
      fwty_ptr=>fwty
      fwtz_ptr=>fwtz
      call checkbound_array(fablo,fabhi,fwtx_ptr,0,0)
      call checkbound_array(fablo,fabhi,fwty_ptr,0,1)
      call checkbound_array(fablo,fabhi,fwtz_ptr,0,AMREX_SPACEDIM-1)

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       do veldir=1,nsolve

        facewtsum=fwtx(D_DECL(i,j,k),veldir)+fwtx(D_DECL(i+1,j,k),veldir)+ &
                  fwty(D_DECL(i,j,k),veldir)+fwty(D_DECL(i,j+1,k),veldir)
        if (SDIM.eq.3) then
         facewtsum=facewtsum+ &
          fwtz(D_DECL(i,j,k),veldir)+fwtz(D_DECL(i,j,k+1),veldir)
        endif

        if (maskcov(D_DECL(i,j,k)).eq.zero) then ! covered by finer cell

         offdiagsum=bx(D_DECL(i,j,k),veldir)+bx(D_DECL(i+1,j,k),veldir)+ &
                    by(D_DECL(i,j,k),veldir)+by(D_DECL(i,j+1,k),veldir)
         if (SDIM.eq.3) then
          offdiagsum=offdiagsum+ &
           bz(D_DECL(i,j,k),veldir)+bz(D_DECL(i,j,k+1),veldir)
         endif

         local_diag=alpha(D_DECL(i,j,k),veldir)+offdiagsum
         if (local_diag.gt.zero) then
          masksolv(D_DECL(i,j,k))=one
         else if (local_diag.eq.zero) then
          masksolv(D_DECL(i,j,k))=zero
         else
          print *,"local_diag invalid"
          stop
         endif

        else if (maskcov(D_DECL(i,j,k)).eq.one) then

         offdiagsum=zero

         do dir=1,SDIM

          ii=0
          jj=0
          kk=0
          if (dir.eq.1) then
           ii=1
           inormal=i
          else if (dir.eq.2) then
           jj=1
           inormal=j
          else if ((dir.eq.3).and.(SDIM.eq.3)) then
           kk=1
           inormal=k
          else
           print *,"dir invalid"
           stop
          endif

          do side=1,2

           if (side.eq.1) then
            iface=i
            jface=j
            kface=k
            icell=i-ii
            jcell=j-jj
            kcell=k-kk
           else if (side.eq.2) then
            iface=i+ii
            jface=j+jj
            kface=k+kk
            icell=i+ii
            jcell=j+jj
            kcell=k+kk
           else
            print *,"side invalid"
            stop
           endif

           if (dir.eq.1) then
            bface=bx(D_DECL(iface,jface,kface),veldir)
           else if (dir.eq.2) then
            bface=by(D_DECL(iface,jface,kface),veldir)
           else if ((dir.eq.3).and.(SDIM.eq.3)) then
            bface=bz(D_DECL(iface,jface,kface),veldir)
           else
            print *,"dir invalid"
            stop
           endif
           offdiagsum=offdiagsum+bface

          enddo ! side=1,2

         enddo ! dir=1..sdim

         local_diag=alpha(D_DECL(i,j,k),veldir)+offdiagsum
         if (local_diag.gt.zero) then
          masksolv(D_DECL(i,j,k))=one
         else if (local_diag.eq.zero) then
          masksolv(D_DECL(i,j,k))=zero
         else
          print *,"local_diag invalid"
          stop
         endif

        else
         print *,"maskcov invalid"
         stop
        endif

        if ((facewtsum.ge.zero).and.(offdiagsum.ge.zero)) then

         if ((facewtsum.gt.zero).and.(offdiagsum.gt.zero)) then
          maskdivres(D_DECL(i,j,k))=one
          maskres(D_DECL(i,j,k))=one
         else if ((facewtsum.eq.zero).or.(offdiagsum.eq.zero)) then
          maskdivres(D_DECL(i,j,k))=zero
          if (alpha(D_DECL(i,j,k),veldir).gt.zero) then
           maskres(D_DECL(i,j,k))=one
          else if (alpha(D_DECL(i,j,k),veldir).eq.zero) then
           maskres(D_DECL(i,j,k))=zero
           if (mdot(D_DECL(i,j,k),veldir).eq.zero) then
            ! do nothing
           else
            print *,"mdot invalid in nsgenerate"
            print *,"level,finest_level ",level,finest_level
            print *,"i,j,k,mdot ",i,j,k,mdot(D_DECL(i,j,k),veldir)
            print *,"i,j,k,maskcov ",i,j,k,maskcov(D_DECL(i,j,k))
            stop
           endif
          else
           print *,"alpha invalid"
           print *,"i,j,k= ",i,j,k
           print *,"alpha=",alpha(D_DECL(i,j,k),veldir)
           stop
          endif
         else
          print *,"facewtsum or offdiagsum invalid"
          print *,"i,j,k ",i,j,k
          print *,"facewtsum= ",facewtsum
          print *,"offdiagsum= ",offdiagsum
          print *,"fwtx Left ",fwtx(D_DECL(i,j,k),veldir)
          print *,"fwtx Right ",fwtx(D_DECL(i+1,j,k),veldir)
          print *,"fwty front ",fwty(D_DECL(i,j,k),veldir)
          print *,"fwty back ",fwty(D_DECL(i,j+1,k),veldir)
          print *,"fwtz bottom ",fwtz(D_DECL(i,j,k),veldir)
          print *,"fwtz top ",fwtz(D_DECL(i,j,k+1),veldir)
          stop
         endif

        else
         print *,"facewtsum or offdiagsum invalid"
         stop
        endif

        local_diag=alpha(D_DECL(i,j,k),veldir)+offdiagsum

        local_diag_check1=facewtsum
          ! offdiagcheck is initialized in fort_buildfacewt which
          ! is declared in LEVELSET_3D.F90
        local_diag_check2=offdiagcheck(D_DECL(i,j,k),veldir)

        if ((local_diag_check1.eq.zero).and. &
            (local_diag_check2.eq.zero)) then
         ! do nothing
        else if ((local_diag_check1.gt.zero).and. &
                 (local_diag_check2.gt.zero)) then
         ! do nothing
        else
         print *,"local_diag_check bust"
         print *,"local_diag_check1=",local_diag_check1
         print *,"local_diag_check2=",local_diag_check2
         stop
        endif

          ! sanity check
        if (alpha(D_DECL(i,j,k),veldir).ge.zero) then
          ! do nothing
        else
         print *,"alpha should be nonneg"
         stop
        endif

       enddo ! veldir=1..nsolve
          
      enddo
      enddo
      enddo

      return
      end subroutine fort_init_mask_sing


       subroutine fort_scalarcoeff( &
         nsolve, &
         xlo,dx, &
         offdiagcheck, &
         DIMS(offdiagcheck), &
         cterm,DIMS(cterm), &
         c2,DIMS(c2), &
         DeDT,DIMS(DeDT), &
         lsnew,DIMS(lsnew), &
         den,DIMS(den), &
         mu,DIMS(mu), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         finest_level, &
         visc_coef, &
         angular_velocity, &
         dt, &
         cur_time, &
         project_option, &
         rzflag, &
         solidheat_flag) &
       bind(c,name='fort_scalarcoeff')

       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, INTENT(in) :: nsolve
       INTEGER_T, INTENT(in) :: level
       INTEGER_T, INTENT(in) :: finest_level
       INTEGER_T, INTENT(in) :: solidheat_flag
       REAL_T, INTENT(in) :: xlo(SDIM)
       REAL_T, INTENT(in) :: dx(SDIM)
       INTEGER_T, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
       INTEGER_T, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
       INTEGER_T, INTENT(in) :: bfact
       INTEGER_T :: growlo(3), growhi(3)
       INTEGER_T, INTENT(in) :: DIMDEC(offdiagcheck)
       INTEGER_T, INTENT(in) :: DIMDEC(cterm)
       INTEGER_T, INTENT(in) :: DIMDEC(c2)
       INTEGER_T, INTENT(in) :: DIMDEC(DeDT)
       INTEGER_T, INTENT(in) :: DIMDEC(lsnew)
       INTEGER_T, INTENT(in) :: DIMDEC(den)
       INTEGER_T, INTENT(in) :: DIMDEC(mu)
       REAL_T, INTENT(in) :: visc_coef,angular_velocity
       REAL_T, INTENT(in) :: dt
       REAL_T, INTENT(in) :: cur_time
       INTEGER_T, INTENT(in) :: project_option,rzflag

       REAL_T, INTENT(in),target :: mu(DIMV(mu))
       REAL_T, pointer :: mu_ptr(D_DECL(:,:,:))
       REAL_T, INTENT(in),target :: den(DIMV(den))
       REAL_T, pointer :: den_ptr(D_DECL(:,:,:))
       REAL_T, INTENT(in),target :: offdiagcheck(DIMV(offdiagcheck),nsolve)
       REAL_T, pointer :: offdiagcheck_ptr(D_DECL(:,:,:),:)
       REAL_T, INTENT(out),target :: cterm(DIMV(cterm),nsolve)
       REAL_T, pointer :: cterm_ptr(D_DECL(:,:,:),:)
       REAL_T, INTENT(in),target :: c2(DIMV(c2),2)
       REAL_T, pointer :: c2_ptr(D_DECL(:,:,:),:)
       REAL_T, INTENT(in),target :: DeDT(DIMV(DeDT))
       REAL_T, pointer :: DeDT_ptr(D_DECL(:,:,:))
       REAL_T, INTENT(in),target :: lsnew(DIMV(lsnew),num_materials*(SDIM+1))
       REAL_T, pointer :: lsnew_ptr(D_DECL(:,:,:),:)

       INTEGER_T i,j,k
       INTEGER_T in_rigid
       INTEGER_T in_prescribed
       INTEGER_T veldir
       INTEGER_T dir_local
       REAL_T rigid_mask
       REAL_T prescribed_mask
       REAL_T den_inverse,dedt_inverse
       REAL_T xsten(-1:1,SDIM)
       INTEGER_T nhalf
       REAL_T local_cterm(nsolve)
       INTEGER_T im
       INTEGER_T velcomp
       REAL_T LSTEST
       REAL_T local_diag
       REAL_T xclamped(SDIM)
       REAL_T LS_clamped
       REAL_T vel_clamped(SDIM)
       REAL_T temperature_clamped
       INTEGER_T prescribed_flag
       INTEGER_T is_clamped_cell

       nhalf=1

       mu_ptr=>mu
       call checkbound_array1(fablo,fabhi,mu_ptr,1,-1)
       den_ptr=>den
       call checkbound_array1(fablo,fabhi,den_ptr,1,-1)
       offdiagcheck_ptr=>offdiagcheck
       call checkbound_array(fablo,fabhi,offdiagcheck_ptr,0,-1)
       cterm_ptr=>cterm
       call checkbound_array(fablo,fabhi,cterm_ptr,0,-1)
       c2_ptr=>c2
       call checkbound_array(fablo,fabhi,c2_ptr,0,-1)
       DeDT_ptr=>DeDT
       call checkbound_array1(fablo,fabhi,DeDT_ptr,1,-1)
       lsnew_ptr=>lsnew
       call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)

       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if ((solidheat_flag.lt.0).or. &
           (solidheat_flag.gt.2)) then
        print *,"solidheat_flag invalid"
        stop
       endif
       if (rzflag.ne.levelrz) then
        print *,"rzflag invalid"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid20"
        stop
       endif
       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid scalar coeff"
        stop
       endif

       if (dt.gt.zero) then
        ! do nothing
       else
        print *,"dt invalid"
        stop
       endif
       if (cur_time.ge.zero) then
        ! do nothing
       else
        print *,"cur_time invalid"
        stop
       endif

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        call gridsten_level(xsten,i,j,k,level,nhalf)
        do dir_local=1,SDIM
         xclamped(dir_local)=xsten(0,dir_local)
        enddo
         ! LS>0 if clamped
        call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
          vel_clamped,temperature_clamped,prescribed_flag,dx)

        if (LS_clamped.ge.zero) then
         is_clamped_cell=1
        else if (LS_clamped.lt.zero) then
         is_clamped_cell=0
        else
         print *,"LS_clamped invalid"
         stop
        endif

        rigid_mask=one
        in_rigid=0

        prescribed_mask=one
        in_prescribed=0

        if (is_clamped_cell.eq.1) then
         rigid_mask=1.0D+6
         in_rigid=1
         prescribed_mask=1.0D+6
         in_prescribed=1
        else if (is_clamped_cell.eq.0) then

         do im=1,num_materials
          if (is_rigid(im).eq.1) then
           LSTEST=lsnew(D_DECL(i,j,k),im)
           if (LSTEST.ge.zero) then
            rigid_mask=1.0D+6
            in_rigid=1
           else if (LSTEST.lt.zero) then
            ! do nothing
           else
            print *,"LSTEST invalid"
            stop
           endif
          else if (is_rigid(im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(im) invalid"
           stop
          endif
         enddo ! im=1..num_materials

         do im=1,num_materials
          if (is_prescribed(im).eq.1) then
           if (is_rigid(im).eq.1) then
            LSTEST=lsnew(D_DECL(i,j,k),im)
            if (LSTEST.ge.zero) then
             prescribed_mask=1.0D+6
             in_prescribed=1
            else if (LSTEST.lt.zero) then
             ! do nothing
            else
             print *,"LSTEST invalid"
             stop
            endif
           else
            print *,"is_rigid(im) invalid"
            stop
           endif
          else if (is_prescribed(im).eq.0) then
           ! do nothing
          else
           print *,"is_prescribed(im) invalid"
           stop
          endif
         enddo ! im=1..num_materials

        else
         print *,"is_clamped_cell invalid"
         stop
        endif

        if (project_option_projectionF(project_option).eq.1) then

         if (nsolve.ne.1) then
          print *,"nsolveMM invalid 150"
          stop
         endif

         local_cterm(1)=c2(D_DECL(i,j,k),1) ! 1/(rho c^2 dt^2)

        else if (project_option.eq.SOLVETYPE_PRESEXTRAP) then 

         if (nsolve.ne.1) then
          print *,"nsolveMM invalid 150"
          stop
         endif

         local_diag=offdiagcheck(D_DECL(i,j,k),1)

          !if FACECOMP_FACECUT component (c++) >0 then facecut_extend=0
          !if FACECOMP_FACECUT component (c++) =0 then facecut_extend=1
          !i.e. if both adjoining cells are fluid, then
          !facecut_extend=0, otherwise, if at least
          !one adjoining cell is solid, then 
          !facecut_extend=1.
          !
          !at least one adjoining face with both adjoining
          !cells as fluid cells.
          !i.e. a fluid cell with at least 1 surrounding fluid
          !cell.
         if ((local_diag.ge.zero).and. &
             (local_diag.le.two*SDIM-one)) then
          local_cterm(1)=1.0D+6

          !all adjoining faces have an adjoining solid cell.
          !e.g. fluid cell completely surrounded by solid
          !cells or a solid cell with any kind of neighbor.
         else if (local_diag.eq.two*SDIM) then 
          local_cterm(1)=zero
         else
          print *,"local_diag invalid"
          stop
         endif

        else if (project_option.eq.SOLVETYPE_HEAT) then ! temperature diffusion

         if (dt.gt.zero) then
          ! do nothing
         else
          print *,"dt invalid"
          stop
         endif
          ! The variable "DeDT" stores: 1/(den cv)  
          ! Note that the derivative, De/DT, satifies: De/DT=cv
          ! In the c++ data structure:
          !  DeDT is declared as localMF[CELL_DEDT_MF] and
          !  initialized in: fort_init_physics_vars
          ! fort_init_physics_vars calls "DeDT_material"
          ! "DeDT_material" calls "INTERNAL_material" for default case.
         dedt_inverse=DeDT(D_DECL(i,j,k))
         if (dedt_inverse.gt.zero) then
          ! do nothing
         else
          print *,"dedt_inverse invalid"
          stop
         endif

         local_cterm(1)=one/(dt*dedt_inverse) ! den cv / dt

         if (is_clamped_cell.eq.1) then
          ! do nothing
         else if (is_clamped_cell.eq.0) then

           ! solidheat_flag==0 diffuse in solid
           ! solidheat_flag==1 dirichlet bc at solid-fluid
           ! solidheat_flag==2 insulating bc at solid-fluid
          if (solidheat_flag.eq.0) then
           ! do nothing 
          else if (solidheat_flag.eq.2) then ! face coeff==0.0 at solid/fluid
           ! do nothing 
          else if (solidheat_flag.eq.1) then
           ! T=Tsolid in solid so coeff>>1
           ! rigid_mask>>1 if solid material occupies cell. (rigid_mask==1
           ! otherwise)
           local_cterm(1)=local_cterm(1)*rigid_mask 
          else
           print *,"solidheat_flag invalid"
           stop
          endif

         else
          print *,"is_clamped_cell invalid"
          stop
         endif


         if (DEBUG_THERMAL_WEIGHT.eq.1) then
          if ((j.eq.32).or.(j.eq.96)) then
           if (i.ge.25) then
            print *,"i,j,local_cterm,rigid_mask ", &
                    i,j,local_cterm(1),rigid_mask
           endif
          endif
         endif

        else if ((project_option.ge.SOLVETYPE_VELEXTRAP).and. &
                 (project_option.lt.SOLVETYPE_VELEXTRAP+num_materials)) then 

         if (dt.ge.zero) then
          ! do nothing
         else
          print *,"dt invalid"
          stop
         endif
         if (nsolve.ne.1) then
          print *,"nsolveMM invalid 150"
          stop
         endif
         local_diag=den(D_DECL(i,j,k))
         if ((local_diag.gt.zero).and. &
             (local_diag.le.two*SDIM)) then
          local_cterm(1)=1.0D+6
         else if (local_diag.eq.zero) then
          local_cterm(1)=zero
         else
          print *,"local_diag invalid"
          stop
         endif

        else if (project_option.eq.SOLVETYPE_VISC) then ! viscosity

         if (nsolve.ne.SDIM) then
          print *,"nsolve invalid21"
          stop
         endif

         if (dt.gt.zero) then
          ! do nothing
         else
          print *,"dt invalid"
          stop
         endif
         den_inverse=den(D_DECL(i,j,k)) ! 1/den
         if (den_inverse.gt.zero) then
          ! do nothing
         else
          print *,"den_inverse invalid"
          stop
         endif
         do veldir=0,nsolve-1
          velcomp=veldir+1
          if (in_prescribed.eq.1) then
           local_cterm(velcomp)=one/(den_inverse*dt)
          else if (in_prescribed.eq.0) then
           local_cterm(velcomp)=one/(den_inverse*dt) ! den/dt
           if (levelrz.eq.COORDSYS_CARTESIAN) then
            ! do nothing
           else if (levelrz.eq.COORDSYS_RZ) then
            if (SDIM.ne.2) then
             print *,"dimension bust"
             stop
            endif
            if (veldir.eq.0) then
             if (xsten(0,1).le.zero) then
              print *,"r invalid"
              stop
             endif
            endif
           else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
            if (veldir.eq.0) then
             if (xsten(0,1).le.zero) then
              print *,"r invalid"
              stop
             endif
            endif
           else 
            print *,"levelrz invalid scalarcoeff"
            stop
           endif

          else
           print *,"in_prescribed invalid"
           stop
          endif
          ! vel=velsolid in solid, so cterm>>1 there.
          local_cterm(velcomp)=local_cterm(velcomp)*prescribed_mask
         enddo ! veldir=0..nsolve-1
        else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                 (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then
         if (dt.gt.zero) then
          ! do nothing
         else
          print *,"dt invalid"
          stop
         endif
         den_inverse=den(D_DECL(i,j,k)) ! 1/den
         if (den_inverse.gt.zero) then
          ! do nothing
         else
          print *,"den_inverse invalid"
          stop
         endif
          ! diffuse mass fraction
         local_cterm(1)=one/(den_inverse*dt) ! den/dt
        else
         print *,"project_option invalid scalar coeff"
         stop
        endif 

        do veldir=1,nsolve
         local_diag=offdiagcheck(D_DECL(i,j,k),veldir)
         if (local_diag.ge.zero) then
          cterm(D_DECL(i,j,k),veldir)=local_cterm(veldir)
         else
          print *,"offdiagcheck invalid"
          stop
         endif
        enddo

       enddo
       enddo
       enddo
 
       return
       end subroutine fort_scalarcoeff


       subroutine fort_restore_pres( &
         offdiagcheck, &
         DIMS(offdiagcheck), &
         savepres,DIMS(savepres), &
         newpres,DIMS(newpres), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         finest_level) &
       bind(c,name='fort_restore_pres')

       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, INTENT(in) :: level
       INTEGER_T, INTENT(in) :: finest_level
       INTEGER_T, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
       INTEGER_T, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
       INTEGER_T, INTENT(in) :: bfact
       INTEGER_T :: growlo(3), growhi(3)
       INTEGER_T, INTENT(in) :: DIMDEC(offdiagcheck)
       INTEGER_T, INTENT(in) :: DIMDEC(savepres)
       INTEGER_T, INTENT(in) :: DIMDEC(newpres)

       REAL_T, INTENT(in),target :: offdiagcheck(DIMV(offdiagcheck))
       REAL_T, pointer :: offdiagcheck_ptr(D_DECL(:,:,:))

       REAL_T, INTENT(in),target :: savepres(DIMV(savepres))
       REAL_T, pointer :: savepres_ptr(D_DECL(:,:,:))
       REAL_T, INTENT(out),target :: newpres(DIMV(newpres))
       REAL_T, pointer :: newpres_ptr(D_DECL(:,:,:))

       INTEGER_T i,j,k
       REAL_T local_diag

       offdiagcheck_ptr=>offdiagcheck
       savepres_ptr=>savepres
       newpres_ptr=>newpres
       call checkbound_array1(fablo,fabhi,offdiagcheck_ptr,0,-1)
       call checkbound_array1(fablo,fabhi,savepres_ptr,0,-1)
       call checkbound_array1(fablo,fabhi,newpres_ptr,0,-1)

       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid scalar coeff"
        stop
       endif

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

         local_diag=offdiagcheck(D_DECL(i,j,k))

          !at least one adjoining face with both adjoining
          !cells as fluid cells.
          !i.e. a fluid cell with at least 1 surrounding fluid
          !cell.
         if ((local_diag.ge.zero).and. &
             (local_diag.le.two*SDIM-one)) then
          newpres(D_DECL(i,j,k))=savepres(D_DECL(i,j,k))
         else if (local_diag.eq.two*SDIM) then
          ! do nothing
         else
          print *,"local_diag invalid"
          stop
         endif

       enddo
       enddo
       enddo
 
       return
       end subroutine fort_restore_pres

       subroutine fort_dividedx( &
         nsolve, &
         bx,DIMS(bx), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact,level, &
         xlo,dx,dir) &
       bind(c,name='fort_dividedx')

       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, INTENT(in) :: nsolve
       INTEGER_T, INTENT(in) :: dir,level
       INTEGER_T, INTENT(in) :: DIMDEC(bx)
       INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, INTENT(in) :: bfact
       REAL_T, INTENT(inout),target :: bx(DIMV(bx),nsolve)
       REAL_T, pointer :: bx_ptr(D_DECL(:,:,:),:)
       REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)

       INTEGER_T i,j,k,n
       REAL_T xsten(-1:1,SDIM)
       INTEGER_T nhalf
       REAL_T hx,RR

       nhalf=1

       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid23"
        stop
       endif

       if ((dir.lt.0).or.(dir.ge.SDIM)) then
        print *,"dir invalid dividedx"
        stop
       endif
       bx_ptr=>bx
       call checkbound_array(fablo,fabhi,bx_ptr,0,dir)

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir)
        RR=one
        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
          ! do nothing
        else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
         if (dir.eq.1) then
          RR=xsten(0,1)
         endif
        else
         print *,"levelrz invalid dividedx"
         stop
        endif
        hx=(xsten(1,dir+1)-xsten(-1,dir+1))*RR
        if (hx.le.zero) then
         print *,"hx invalid"
         stop
        endif
        do n=1,nsolve
         bx(D_DECL(i,j,k),n)=bx(D_DECL(i,j,k),n)/hx
        enddo

       enddo
       enddo
       enddo
 
       return
       end subroutine fort_dividedx


       subroutine fort_regularize_bx( &
         nsolve, &
         bx,DIMS(bx), &
         min_interior_coeff, &
         domlo,domhi, &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         xlo,dx, &
         dir) &
       bind(c,name='fort_regularize_bx')
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, INTENT(in) :: nsolve
       INTEGER_T, INTENT(in) :: dir
       INTEGER_T, INTENT(in) :: level
       INTEGER_T, INTENT(in) :: DIMDEC(bx)
       INTEGER_T, INTENT(in) :: domlo(SDIM),domhi(SDIM)
       INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, INTENT(in) :: bfact
       REAL_T, INTENT(in) :: min_interior_coeff

       REAL_T, INTENT(inout),target :: bx(DIMV(bx),nsolve)
       REAL_T, pointer :: bx_ptr(D_DECL(:,:,:),:)

       REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)

       INTEGER_T i,j,k,n
       INTEGER_T inorm
       REAL_T local_bx

       if (min_interior_coeff.gt.zero) then
        ! do nothing
       else
        print *,"min_interior_coeff invalid"
        stop
       endif

       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid23"
        stop
       endif

       if ((dir.lt.0).or.(dir.ge.SDIM)) then
        print *,"dir invalid mult_facewt"
        stop
       endif

       bx_ptr=>bx
       call checkbound_array(fablo,fabhi,bx_ptr,0,dir)

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
       
        if (dir.eq.0) then
         inorm=i
        else if (dir.eq.1) then
         inorm=j
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         inorm=k
        else
         print *,"dir invalid"
         stop
        endif

        do n=1,nsolve
         local_bx=bx(D_DECL(i,j,k),n)
         if (local_bx.eq.zero) then
          if ((inorm.gt.domlo(dir+1)).and. &
              (inorm.lt.domhi(dir+1)+1)) then
           bx(D_DECL(i,j,k),n)=min_interior_coeff
          else if ((inorm.eq.domlo(dir+1)).or. &
                   (inorm.eq.domhi(dir+1)+1)) then
           ! do nothing
          else
           print *,"inorm invalid"
           stop
          endif
         else if (local_bx.gt.zero) then
          ! do nothing
         else
          print *,"local_bx invalid"
          stop
         endif
        enddo ! n=1,nsolve

       enddo
       enddo
       enddo
 
       return
       end subroutine fort_regularize_bx


       subroutine fort_mult_facewt( &
         nsolve, &
         bx,DIMS(bx), &
         facewt,DIMS(facewt), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         xlo,dx, &
         dir) &
       bind(c,name='fort_mult_facewt')
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, INTENT(in) :: nsolve
       INTEGER_T, INTENT(in) :: dir
       INTEGER_T, INTENT(in) :: level
       INTEGER_T, INTENT(in) :: DIMDEC(bx)
       INTEGER_T, INTENT(in) :: DIMDEC(facewt)
       INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, INTENT(in) :: bfact
       REAL_T, INTENT(out),target :: bx(DIMV(bx),nsolve)
       REAL_T, pointer :: bx_ptr(D_DECL(:,:,:),:)
       REAL_T, INTENT(in),target :: facewt(DIMV(facewt),nsolve)
       REAL_T, pointer :: facewt_ptr(D_DECL(:,:,:),:)
       REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)

       INTEGER_T i,j,k,n

       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid23"
        stop
       endif

       if ((dir.lt.0).or.(dir.ge.SDIM)) then
        print *,"dir invalid mult_facewt"
        stop
       endif
       bx_ptr=>bx
       facewt_ptr=>facewt
       call checkbound_array(fablo,fabhi,bx_ptr,0,dir)
       call checkbound_array(fablo,fabhi,facewt_ptr,0,dir)

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        
        do n=1,nsolve
         bx(D_DECL(i,j,k),n)=bx(D_DECL(i,j,k),n)*facewt(D_DECL(i,j,k),n)
        enddo ! n=1,nsolveMM

       enddo
       enddo
       enddo
 
       return
       end subroutine fort_mult_facewt

      subroutine fort_interpmac( &
        bfact,bfact_f, &
        fdata, DIMS(fdata), &
        cdata, DIMS(cdata), &
        lo, hi,  &
        cdiag,DIMS(cdiag) ) &
      bind(c,name='fort_interpmac')

      use global_utility_module
      IMPLICIT NONE
      INTEGER_T, INTENT(in) :: DIMDEC(fdata)
      INTEGER_T, INTENT(in) :: DIMDEC(cdata)
      INTEGER_T, INTENT(in) :: DIMDEC(cdiag)
      INTEGER_T, INTENT(in) :: lo(AMREX_SPACEDIM)
      INTEGER_T, INTENT(in) :: hi(AMREX_SPACEDIM)
      INTEGER_T lof(SDIM),hif(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      REAL_T, INTENT(inout),target :: fdata(DIMV(fdata))
      REAL_T, pointer :: fdata_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: cdata(DIMV(cdata))
      REAL_T, pointer :: cdata_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: cdiag(DIMV(cdiag))
      REAL_T, pointer :: cdiag_ptr(D_DECL(:,:,:))
      INTEGER_T, INTENT(in) :: bfact,bfact_f

      INTEGER_T ic,jc,kc,ifine,jfine,kfine
      INTEGER_T cvalid,dir
      REAL_T fine_value,voltotal,volall
      REAL_T wt(SDIM)
 

      if (bfact.lt.1) then
       print *,"bfact invalid110"
       stop
      endif
      if (bfact_f.lt.1) then
       print *,"bfact_f invalid3 ",bfact_f
       stop
      endif
      if ((bfact.ne.bfact_f).and.(bfact.ne.2*bfact_f)) then
       print *,"bfact or bfact_f invalid"
       stop
      endif
      do dir=1,SDIM
       lof(dir)=2*lo(dir)
       hif(dir)=2*(hi(dir)+1)-1
      enddo

      cdata_ptr=>cdata
      fdata_ptr=>fdata
      cdiag_ptr=>cdiag
      call checkbound_array1(lo,hi,cdata_ptr,0,-1)
      call checkbound_array1(lof,hif,fdata_ptr,0,-1)
      call checkbound_array1(lo,hi,cdiag_ptr,0,-1)

      call growntilebox(lof,hif,lof,hif,growlo,growhi,0) 
 
      do kfine=growlo(3),growhi(3)
      do jfine=growlo(2),growhi(2)
      do ifine=growlo(1),growhi(1)
       fine_value=zero
       voltotal=zero
       call coarse_subelement_stencil(ifine,jfine,kfine,stenlo,stenhi, &
         bfact,bfact_f)
       do ic=stenlo(1),stenhi(1)
        call intersect_weight_interp(ic,ifine,bfact,bfact_f,wt(1))
        if (wt(1).gt.zero) then
         do jc=stenlo(2),stenhi(2)
          call intersect_weight_interp(jc,jfine,bfact,bfact_f,wt(2))
          if (wt(2).gt.zero) then
           do kc=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_interp(kc,kfine,bfact,bfact_f,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
              ! mask_residual>0?
             if (cdiag(D_DECL(ic,jc,kc)).gt.zero) then 
              cvalid=1
             else if (cdiag(D_DECL(ic,jc,kc)).eq.zero) then
              cvalid=0
             else
              print *,"cdiag invalid"
              print *,"cdiag= ",cdiag(D_DECL(ic,jc,kc))
              stop
             endif
             if (cvalid.eq.1) then
              volall=wt(1)
              do dir=2,SDIM
               volall=volall*wt(dir)
              enddo
              fine_value=fine_value+volall*cdata(D_DECL(ic,jc,kc))
              voltotal=voltotal+volall
             endif
            endif
           enddo ! kc
          endif
         enddo ! jc
        endif
       enddo ! ic

       if (voltotal.gt.zero) then
        fine_value=fine_value/voltotal
        fdata(D_DECL(ifine,jfine,kfine))= &
           fdata(D_DECL(ifine,jfine,kfine))+fine_value  
       else if (voltotal.eq.zero) then
        ! do nothing
       else
        print *,"voltotal invalid" 
        stop
       endif

      end do ! ifine
      end do ! jfine
      end do ! kfine

      return
      end subroutine fort_interpmac


! maskcov=1 outside domain and on fine-fine ghost cells
! maskcov=1 for interior cells not covered by a finer cell
! fwtx,fwty,fwtz are not averaged down.
! bx,by,bz are equal to bx_noarea * area/dx and bx_noarea is averaged down.
      subroutine fort_nsgenerate( &
       level, &
       finest_level, &
       nsolve, &
       project_option, &
       alpha,DIMS(alpha), &
       diag_reg, &
       DIMS(diag_reg), &
       bx,DIMS(bx), &
       by,DIMS(by), &
       bz,DIMS(bz), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact) &
      bind(c,name='fort_nsgenerate')
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: nsolve
      INTEGER_T, INTENT(in) :: project_option
      INTEGER_T, INTENT(in) :: DIMDEC(alpha)
      INTEGER_T, INTENT(in) :: DIMDEC(diag_reg)
      INTEGER_T, INTENT(in) :: DIMDEC(bx)
      INTEGER_T, INTENT(in) :: DIMDEC(by)
      INTEGER_T, INTENT(in) :: DIMDEC(bz)
      INTEGER_T, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T             :: growlo(3), growhi(3)

      REAL_T, INTENT(in),target :: alpha(DIMV(alpha),nsolve)
      REAL_T, pointer :: alpha_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out),target :: diag_reg(DIMV(diag_reg),nsolve)
      REAL_T, pointer :: diag_reg_ptr(D_DECL(:,:,:),:)
      ! coeff * areafrac * areaface / (dxfrac*dx)  (if coeff>0.0)
      REAL_T, INTENT(in),target :: bx(DIMV(bx),nsolve) 
      REAL_T, INTENT(in),target :: by(DIMV(by),nsolve)
      REAL_T, INTENT(in),target :: bz(DIMV(bz),nsolve)
      REAL_T, pointer :: bx_ptr(D_DECL(:,:,:),:)
      REAL_T, pointer :: by_ptr(D_DECL(:,:,:),:)
      REAL_T, pointer :: bz_ptr(D_DECL(:,:,:),:)

      INTEGER_T i,j,k
      INTEGER_T veldir
      REAL_T offdiagsum
      REAL_T local_diag

      if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
       print *,"nsolve invalid24"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid nsgenerate"
       stop
      endif

      alpha_ptr=>alpha
      call checkbound_array(fablo,fabhi,alpha_ptr,0,-1)
      diag_reg_ptr=>diag_reg
      call checkbound_array(fablo,fabhi,diag_reg_ptr,0,-1)
      bx_ptr=>bx
      by_ptr=>by
      bz_ptr=>bz
      call checkbound_array(fablo,fabhi,bx_ptr,0,0)
      call checkbound_array(fablo,fabhi,by_ptr,0,1)
      call checkbound_array(fablo,fabhi,bz_ptr,0,AMREX_SPACEDIM-1)

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       do veldir=1,nsolve

        offdiagsum=bx(D_DECL(i,j,k),veldir)+bx(D_DECL(i+1,j,k),veldir)+ &
                   by(D_DECL(i,j,k),veldir)+by(D_DECL(i,j+1,k),veldir)
        if (SDIM.eq.3) then
         offdiagsum=offdiagsum+ &
          bz(D_DECL(i,j,k),veldir)+bz(D_DECL(i,j,k+1),veldir)
        endif

        local_diag=alpha(D_DECL(i,j,k),veldir)+offdiagsum
        if (local_diag.gt.zero) then
         diag_reg(D_DECL(i,j,k),veldir)=local_diag
        else
         print *,"local_diag invalid"
         stop
        endif

        if (alpha(D_DECL(i,j,k),veldir).ge.zero) then
          ! do nothing
        else
         print *,"alpha should be nonneg"
         stop
        endif
        if (offdiagsum.ge.zero) then
          ! do nothing
        else
         print *,"offdiagsum should be nonneg"
         stop
        endif

       enddo ! veldir=1..nsolveMM
          
      enddo
      enddo
      enddo

      return
      end subroutine fort_nsgenerate

      end module macoperator_module
