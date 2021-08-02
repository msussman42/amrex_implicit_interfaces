#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif


#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

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



       subroutine fort_scalarcoeff( &
         nsolve, &
         nmat, &
         xlo,dx, &
         offdiagcheck, &
         DIMS(offdiagcheck), &
         cterm,DIMS(cterm), &
         c2,DIMS(c2), &
         DeDT,DIMS(DeDT), &
         recon,DIMS(recon), &
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

       use probf90_module
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, intent(in) :: nsolve
       INTEGER_T, intent(in) :: level
       INTEGER_T, intent(in) :: finest_level
       INTEGER_T, intent(in) :: nmat
       INTEGER_T, intent(in) :: solidheat_flag
       REAL_T, intent(in) :: xlo(SDIM)
       REAL_T, intent(in) :: dx(SDIM)
       INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
       INTEGER_T, intent(in) :: bfact
       INTEGER_T :: growlo(3), growhi(3)
       INTEGER_T, intent(in) :: DIMDEC(offdiagcheck)
       INTEGER_T, intent(in) :: DIMDEC(cterm)
       INTEGER_T, intent(in) :: DIMDEC(c2)
       INTEGER_T, intent(in) :: DIMDEC(DeDT)
       INTEGER_T, intent(in) :: DIMDEC(recon)
       INTEGER_T, intent(in) :: DIMDEC(lsnew)
       INTEGER_T, intent(in) :: DIMDEC(den)
       INTEGER_T, intent(in) :: DIMDEC(mu)
       REAL_T, intent(in) :: visc_coef,angular_velocity
       REAL_T, intent(in) :: dt
       REAL_T, intent(in) :: cur_time
       INTEGER_T, intent(in) :: project_option,rzflag

       REAL_T, intent(in),target :: mu(DIMV(mu))
       REAL_T, pointer :: mu_ptr(D_DECL(:,:,:))
       REAL_T, intent(in),target :: den(DIMV(den))
       REAL_T, pointer :: den_ptr(D_DECL(:,:,:))
       REAL_T, intent(in),target :: offdiagcheck(DIMV(offdiagcheck),nsolve)
       REAL_T, pointer :: offdiagcheck_ptr(D_DECL(:,:,:),:)
       REAL_T, intent(out),target :: cterm(DIMV(cterm),nsolve)
       REAL_T, pointer :: cterm_ptr(D_DECL(:,:,:),:)
       REAL_T, intent(in),target :: c2(DIMV(c2),2)
       REAL_T, pointer :: c2_ptr(D_DECL(:,:,:),:)
       REAL_T, intent(in),target :: DeDT(DIMV(DeDT))
       REAL_T, pointer :: DeDT_ptr(D_DECL(:,:,:))
       REAL_T, intent(in),target :: recon(DIMV(recon),nmat*ngeom_recon)
       REAL_T, pointer :: recon_ptr(D_DECL(:,:,:),:)
       REAL_T, intent(in),target :: lsnew(DIMV(lsnew),nmat*(SDIM+1))
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
       INTEGER_T is_clamped_cell

       nhalf=1

       mu_ptr=>mu
       call checkbound_array1(fablo,fabhi,mu_ptr,1,-1,33)
       den_ptr=>den
       call checkbound_array1(fablo,fabhi,den_ptr,1,-1,33)
       offdiagcheck_ptr=>offdiagcheck
       call checkbound_array(fablo,fabhi,offdiagcheck_ptr,0,-1,33)
       cterm_ptr=>cterm
       call checkbound_array(fablo,fabhi,cterm_ptr,0,-1,33)
       c2_ptr=>c2
       call checkbound_array(fablo,fabhi,c2_ptr,0,-1,33)
       DeDT_ptr=>DeDT
       call checkbound_array1(fablo,fabhi,DeDT_ptr,1,-1,33)
       lsnew_ptr=>lsnew
       call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1,33)
       recon_ptr=>recon
       call checkbound_array(fablo,fabhi,recon_ptr,1,-1,33)

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
       if (nmat.ne.num_materials) then
        print *,"nmat invalid"
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
                vel_clamped,temperature_clamped)

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

         do im=1,nmat
          if (is_rigid(nmat,im).eq.1) then
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
          else if (is_rigid(nmat,im).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif
         enddo ! im=1..nmat

         do im=1,nmat
          if (is_prescribed(nmat,im).eq.1) then
           if (is_rigid(nmat,im).eq.1) then
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
            print *,"is_rigid(nmat,im) invalid"
            stop
           endif
          else if (is_prescribed(nmat,im).eq.0) then
           ! do nothing
          else
           print *,"is_prescribed(nmat,im) invalid"
           stop
          endif
         enddo ! im=1..nmat

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

        else if (project_option.eq.12) then ! pressure extension

         if (nsolve.ne.1) then
          print *,"nsolveMM invalid 150"
          stop
         endif

         local_diag=offdiagcheck(D_DECL(i,j,k),1)

          !if facecut_index>0 then facecut_extend=0
          !if facecut_index=0 then facecut_extend=1
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

        else if (project_option.eq.2) then  ! temperature diffusion

         if (dt.gt.zero) then
          ! do nothing
         else
          print *,"dt invalid"
          stop
         endif
          ! 1/(den cv)  note: De/DT=cv
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

        else if (project_option.eq.200) then  ! smoothing

         if (dt.gt.zero) then
          ! do nothing
         else
          print *,"dt invalid"
          stop
         endif

         local_cterm(1)=one/dt ! den cv / dt

        else if (project_option.eq.3) then ! viscosity

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
           if (levelrz.eq.0) then
            ! do nothing
           else if (levelrz.eq.1) then
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
           else if (levelrz.eq.3) then
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
        else if ((project_option.ge.100).and. &
                 (project_option.lt.100+num_species_var)) then
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


       subroutine FORT_RESTORE_PRES( &
         offdiagcheck, &
         DIMS(offdiagcheck), &
         savepres,DIMS(savepres), &
         newpres,DIMS(newpres), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         finest_level)
       use probf90_module
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, intent(in) :: level
       INTEGER_T, intent(in) :: finest_level
       INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
       INTEGER_T, intent(in) :: bfact
       INTEGER_T :: growlo(3), growhi(3)
       INTEGER_T, intent(in) :: DIMDEC(offdiagcheck)
       INTEGER_T, intent(in) :: DIMDEC(savepres)
       INTEGER_T, intent(in) :: DIMDEC(newpres)

       REAL_T, intent(in) :: offdiagcheck(DIMV(offdiagcheck))
       REAL_T, intent(in) :: savepres(DIMV(savepres))
       REAL_T, intent(out) :: newpres(DIMV(newpres))

       INTEGER_T i,j,k
       REAL_T local_diag

       call checkbound(fablo,fabhi,DIMS(offdiagcheck),0,-1,33)
       call checkbound(fablo,fabhi,DIMS(savepres),0,-1,33)
       call checkbound(fablo,fabhi,DIMS(newpres),0,-1,33)

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
       end subroutine FORT_RESTORE_PRES

       subroutine FORT_DIVIDEDX ( &
         nsolve, &
         bx,DIMS(bx), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact,level, &
         xlo,dx,dir)
       use probf90_module
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, intent(in) :: nsolve
       INTEGER_T, intent(in) :: dir,level
       INTEGER_T, intent(in) :: DIMDEC(bx)
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, intent(in) :: bfact
       REAL_T, intent(inout) :: bx(DIMV(bx),nsolve)
       REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

       INTEGER_T i,j,k,n
       REAL_T xsten(-1:1,SDIM)
       INTEGER_T nhalf
       REAL_T hx,RR
       INTEGER_T nmat

       nhalf=1

       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid23"
        stop
       endif

       nmat=num_materials

       if ((dir.lt.0).or.(dir.ge.SDIM)) then
        print *,"dir invalid dividedx"
        stop
       endif
       call checkbound(fablo,fabhi,DIMS(bx),0,dir,33)

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi, &
               0,dir,13) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir,20)
        RR=one
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
          ! do nothing
        else if (levelrz.eq.3) then
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
       end subroutine FORT_DIVIDEDX


       subroutine FORT_REGULARIZE_BX ( &
         nsolve, &
         bx,DIMS(bx), &
         min_interior_coeff, &
         domlo,domhi, &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         xlo,dx, &
         dir)
       use probf90_module
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, intent(in) :: nsolve
       INTEGER_T, intent(in) :: dir
       INTEGER_T, intent(in) :: level
       INTEGER_T, intent(in) :: DIMDEC(bx)
       INTEGER_T, intent(in) :: domlo(SDIM),domhi(SDIM)
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, intent(in) :: bfact
       REAL_T, intent(in) :: min_interior_coeff
       REAL_T, intent(inout) :: bx(DIMV(bx),nsolve)
       REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

       INTEGER_T i,j,k,n
       INTEGER_T inorm
       INTEGER_T nmat
       REAL_T local_bx

       nmat=num_materials

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
       call checkbound(fablo,fabhi,DIMS(bx),0,dir,33)

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir,14) 
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
       end subroutine FORT_REGULARIZE_BX


       subroutine FORT_MULT_FACEWT ( &
         nsolve, &
         bx,DIMS(bx), &
         facewt,DIMS(facewt), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         level, &
         xlo,dx, &
         dir)
       use probf90_module
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, intent(in) :: nsolve
       INTEGER_T, intent(in) :: dir
       INTEGER_T, intent(in) :: level
       INTEGER_T, intent(in) :: DIMDEC(bx)
       INTEGER_T, intent(in) :: DIMDEC(facewt)
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, intent(in) :: bfact
       REAL_T, intent(out) :: bx(DIMV(bx),nsolve)
       REAL_T, intent(in) :: facewt(DIMV(facewt),nsolve)
       REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

       INTEGER_T i,j,k,n
       INTEGER_T nmat

       nmat=num_materials

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
       call checkbound(fablo,fabhi,DIMS(bx),0,dir,33)
       call checkbound(fablo,fabhi,DIMS(facewt),0,dir,33)

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir,15) 
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
       end subroutine FORT_MULT_FACEWT

      subroutine FORT_INTERPMAC( &
        bfact,bfact_f, &
        f, DIMS(f), &
        c, DIMS(c), &
        lo, hi,  &
        cdiag,DIMS(cdiag), &
        fdiag,DIMS(fdiag) )
      use global_utility_module
      IMPLICIT NONE
      INTEGER_T DIMDEC(f)
      INTEGER_T DIMDEC(c)
      INTEGER_T DIMDEC(cdiag)
      INTEGER_T DIMDEC(fdiag)
      INTEGER_T lo(AMREX_SPACEDIM)
      INTEGER_T hi(AMREX_SPACEDIM)
      INTEGER_T lof(SDIM),hif(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      REAL_T f(DIMV(f))
      REAL_T c(DIMV(c))
      REAL_T cdiag(DIMV(cdiag))
      REAL_T fdiag(DIMV(fdiag))
      INTEGER_T bfact,bfact_f

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

      call checkbound(lo,hi,DIMS(c),0,-1,140)
      call checkbound(lof,hif,DIMS(f),0,-1,140)
      call checkbound(lo,hi,DIMS(cdiag),0,-1,140)

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
              fine_value=fine_value+volall*c(D_DECL(ic,jc,kc))
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
        f(D_DECL(ifine,jfine,kfine))=f(D_DECL(ifine,jfine,kfine))+fine_value  
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
      end subroutine FORT_INTERPMAC



! maskcov=1 outside domain and on fine-fine ghost cells
! maskcov=1 for interior cells not covered by a finer cell
! fwtx,fwty,fwtz are not averaged down.
! bx,by,bz are equal to bxcoefnoarea; bxcoefnoarea is not averaged down at this
! point.
      subroutine FORT_INIT_MASK_SING( &
       level, &
       finest_level, &
       nsolve, &
       nmat, &
       project_option, &
       ncphys, &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
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
       bc)
      use probf90_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(masksolv) ! ONES_MF in c++
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(alpha)
      INTEGER_T, intent(in) :: DIMDEC(offdiagcheck)
      INTEGER_T, intent(in) :: DIMDEC(maskdivres)
      INTEGER_T, intent(in) :: DIMDEC(maskres)
      INTEGER_T, intent(in) :: DIMDEC(mdot)
      INTEGER_T, intent(in) :: DIMDEC(bx)
      INTEGER_T, intent(in) :: DIMDEC(by)
      INTEGER_T, intent(in) :: DIMDEC(bz)
      INTEGER_T, intent(in) :: DIMDEC(fwtx)
      INTEGER_T, intent(in) :: DIMDEC(fwty)
      INTEGER_T, intent(in) :: DIMDEC(fwtz)
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T             :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bc(SDIM,2,nsolve)

      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: yface(DIMV(yface),ncphys)
      REAL_T, intent(in) :: zface(DIMV(zface),ncphys)
       ! ONES_MF in c++
      REAL_T, intent(out) :: masksolv(DIMV(masksolv))
      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: alpha(DIMV(alpha),nsolve)
      REAL_T, intent(in) :: offdiagcheck(DIMV(offdiagcheck),nsolve)
      REAL_T, intent(out) :: maskdivres(DIMV(maskdivres))
      REAL_T, intent(out) :: maskres(DIMV(maskres))
      REAL_T, intent(in) :: mdot(DIMV(mdot),nsolve)
       ! coeff * areafrac * areaface / (dxfrac*dx)
      REAL_T, intent(in) :: bx(DIMV(bx),nsolve) 
      REAL_T, intent(in) :: by(DIMV(by),nsolve)
      REAL_T, intent(in) :: bz(DIMV(bz),nsolve)
       ! coeff * areafrac / dxfrac
      REAL_T, intent(in) :: fwtx(DIMV(fwtx),nsolve)  
      REAL_T, intent(in) :: fwty(DIMV(fwty),nsolve)
      REAL_T, intent(in) :: fwtz(DIMV(fwtz),nsolve)

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

       ! indexes start at 0
      if (ncphys.lt.8) then
       print *,"ncphys invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
       print *,"nsolve invalid24"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid nsgenerate"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(xface),0,0,244)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,244)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,244)
       ! ONES_MF in c++
      call checkbound(fablo,fabhi,DIMS(masksolv),0,-1,140)
      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,140)
      call checkbound(fablo,fabhi,DIMS(alpha),0,-1,140)
      call checkbound(fablo,fabhi,DIMS(offdiagcheck),0,-1,140)
      call checkbound(fablo,fabhi, &
       DIMS(maskdivres), &
       0,-1,140)
      call checkbound(fablo,fabhi,DIMS(maskres),0,-1,140)
      call checkbound(fablo,fabhi,DIMS(mdot),0,-1,140)
      call checkbound(fablo,fabhi,DIMS(bx),0,0,140)
      call checkbound(fablo,fabhi,DIMS(by),0,1,140)
      call checkbound(fablo,fabhi,DIMS(bz),0,AMREX_SPACEDIM-1,140)
      call checkbound(fablo,fabhi,DIMS(fwtx),0,0,140)
      call checkbound(fablo,fabhi,DIMS(fwty),0,1,140)
      call checkbound(fablo,fabhi,DIMS(fwtz),0,AMREX_SPACEDIM-1,140)

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
      end subroutine FORT_INIT_MASK_SING





! maskcov=1 outside domain and on fine-fine ghost cells
! maskcov=1 for interior cells not covered by a finer cell
! fwtx,fwty,fwtz are not averaged down.
! bx,by,bz are equal to bx_noarea * area/dx and bx_noarea is averaged down.
      subroutine FORT_NSGENERATE( &
       level, &
       finest_level, &
       nsolve, &
       nmat, &
       project_option, &
       ncphys, &
       alpha,DIMS(alpha), &
       diag_reg, &
       DIMS(diag_reg), &
       bx,DIMS(bx), &
       by,DIMS(by), &
       bz,DIMS(bz), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact)
      use probf90_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: DIMDEC(alpha)
      INTEGER_T, intent(in) :: DIMDEC(diag_reg)
      INTEGER_T, intent(in) :: DIMDEC(bx)
      INTEGER_T, intent(in) :: DIMDEC(by)
      INTEGER_T, intent(in) :: DIMDEC(bz)
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T             :: growlo(3), growhi(3)

      REAL_T, intent(in) :: alpha(DIMV(alpha),nsolve)
      REAL_T, intent(out) :: diag_reg(DIMV(diag_reg),nsolve)
      ! coeff * areafrac * areaface / (dxfrac*dx)  (if coeff>0.0)
      REAL_T, intent(in) :: bx(DIMV(bx),nsolve) 
      REAL_T, intent(in) :: by(DIMV(by),nsolve)
      REAL_T, intent(in) :: bz(DIMV(bz),nsolve)

      INTEGER_T i,j,k
      INTEGER_T veldir
      REAL_T offdiagsum
      REAL_T local_diag

       ! indexes start at 0
      if (ncphys.lt.8) then
       print *,"ncphys invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
       print *,"nsolve invalid24"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid nsgenerate"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(alpha),0,-1,140)
      call checkbound(fablo,fabhi,DIMS(diag_reg),0,-1,140)
      call checkbound(fablo,fabhi,DIMS(bx),0,0,140)
      call checkbound(fablo,fabhi,DIMS(by),0,1,140)
      call checkbound(fablo,fabhi,DIMS(bz),0,AMREX_SPACEDIM-1,140)

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
      end subroutine FORT_NSGENERATE

