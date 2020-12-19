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



       subroutine FORT_SCALARCOEFF( &
         num_materials_face, &
         nsolve, &
         nsolveMM, &
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
         solidheat_flag)
       use probf90_module
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, intent(in) :: num_materials_face
       INTEGER_T, intent(in) :: nsolve
       INTEGER_T, intent(in) :: nsolveMM
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

       REAL_T, intent(in) ::  mu(DIMV(mu),nmat+1)
       REAL_T, intent(in) ::  den(DIMV(den),nmat+1)
       REAL_T, intent(in) ::  offdiagcheck(DIMV(offdiagcheck),nsolveMM)
       REAL_T, intent(out) ::  cterm(DIMV(cterm),nsolveMM)
       REAL_T, intent(in) ::  c2(DIMV(c2),2)
       REAL_T, intent(in) ::  DeDT(DIMV(DeDT),nmat+1)
       REAL_T, intent(in) ::  recon(DIMV(recon),nmat*ngeom_recon)
       REAL_T, intent(in) ::  lsnew(DIMV(lsnew),nmat*(SDIM+1))

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
       REAL_T local_cterm(nsolveMM)
       INTEGER_T im_vel,im
       INTEGER_T velcomp
       REAL_T LSTEST
       REAL_T local_diag
       REAL_T xclamped(SDIM)
       REAL_T LS_clamped
       REAL_T vel_clamped(SDIM)
       REAL_T temperature_clamped
       INTEGER_T is_clamped_cell

       nhalf=1

       call checkbound(fablo,fabhi,DIMS(mu),1,-1,33)
       call checkbound(fablo,fabhi,DIMS(den),1,-1,33)
       call checkbound(fablo,fabhi,DIMS(offdiagcheck),0,-1,33)
       call checkbound(fablo,fabhi,DIMS(cterm),0,-1,33)
       call checkbound(fablo,fabhi,DIMS(c2),0,-1,33)
       call checkbound(fablo,fabhi,DIMS(DeDT),1,-1,33)
       call checkbound(fablo,fabhi,DIMS(lsnew),1,-1,33)
       call checkbound(fablo,fabhi,DIMS(recon),1,-1,33)

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
       if (nsolveMM.ne.nsolve*num_materials_face) then
        print *,"nsolveMM invalid 124"
        stop
       endif
       if ((level.lt.0).or.(level.gt.finest_level)) then
        print *,"level invalid scalar coeff"
        stop
       endif
       if ((num_materials_face.ne.1).and.(num_materials_face.ne.nmat)) then
        print *,"num_materials_face invalid"
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

        if ((project_option.eq.0).or. & ! regular project
            (project_option.eq.1).or. & ! initial project
            (project_option.eq.10).or. & ! sync project
            (project_option.eq.13).or. & ! FSI_material_exists 1st project
            (project_option.eq.11)) then ! FSI_material_exists 2nd project

         if (num_materials_vel.ne.1) then
          print *,"num_materials_vel invalid"
          stop
         endif
         if (num_materials_face.ne.1) then
          print *,"num_materials_face invalid"
          stop
         endif

         if (nsolveMM.ne.num_materials_face) then
          print *,"nsolveMM invalid 150"
          stop
         endif

         im_vel=1
         local_cterm(im_vel)=c2(D_DECL(i,j,k),1) ! 1/(rho c^2 dt^2)

        else if (project_option.eq.12) then ! pressure extension

         if (num_materials_vel.ne.1) then
          print *,"num_materials_vel invalid"
          stop
         endif
         if (num_materials_face.ne.1) then
          print *,"num_materials_face invalid"
          stop
         endif

         if (nsolveMM.ne.num_materials_face) then
          print *,"nsolveMM invalid 150"
          stop
         endif

         im_vel=1
         local_diag=offdiagcheck(D_DECL(i,j,k),im_vel)

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
          local_cterm(im_vel)=1.0D+6

          !all adjoining faces have an adjoining solid cell.
          !e.g. fluid cell completely surrounded by solid
          !cells or a solid cell with any kind of neighbor.
         else if (local_diag.eq.two*SDIM) then 
          local_cterm(im_vel)=zero
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
         dedt_inverse=DeDT(D_DECL(i,j,k),1)
         if (dedt_inverse.le.zero) then
          print *,"dedt_inverse invalid"
          stop
         endif
         if (num_materials_face.ne.num_materials_scalar_solve) then
          print *,"num_materials_face invalid"
          stop
         endif

         do im_vel=1,num_materials_scalar_solve

          local_cterm(im_vel)=one/(dt*dedt_inverse) ! den cv / dt

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
            local_cterm(im_vel)=local_cterm(im_vel)*rigid_mask 
           else
            print *,"solidheat_flag invalid"
            stop
           endif

          else
           print *,"is_clamped_cell invalid"
           stop
          endif

         enddo ! im_vel=1,num_materials_scalar_solve

         if (DEBUG_THERMAL_WEIGHT.eq.1) then
          if ((j.eq.32).or.(j.eq.96)) then
           if (i.ge.25) then
            print *,"i,j,local_cterm,rigid_mask ", &
                    i,j,local_cterm(1),rigid_mask
           endif
          endif
         endif

        else if (project_option.eq.3) then ! viscosity

         if (nsolve.ne.SDIM) then
          print *,"nsolve invalid21"
          stop
         endif
         if (num_materials_vel.eq.1) then
          ! do nothing
         else
          print *,"num_materials_vel invalid"
          stop
         endif
         if (num_materials_face.eq.1) then
          ! do nothing
         else
          print *,"num_materials_face invalid"
          stop
         endif

         if (dt.gt.zero) then
          ! do nothing
         else
          print *,"dt invalid"
          stop
         endif
         den_inverse=den(D_DECL(i,j,k),1) ! 1/den
         if (den_inverse.le.zero) then
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
         den_inverse=den(D_DECL(i,j,k),1) ! 1/den
         if (den_inverse.le.zero) then
          print *,"den_inverse invalid"
          stop
         endif
         if (num_materials_face.ne.num_materials_scalar_solve) then
          print *,"num_materials_face invalid"
          stop
         endif
          ! diffuse mass fraction
         do im_vel=1,num_materials_scalar_solve
          local_cterm(im_vel)=one/(den_inverse*dt) ! den/dt
         enddo
        else
         print *,"project_option invalid scalar coeff"
         stop
        endif 

        do veldir=1,nsolveMM
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
       end subroutine FORT_SCALARCOEFF


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
         num_materials_face, &
         nsolve, &
         nsolveMM, &
         bx,DIMS(bx), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact,level, &
         xlo,dx,dir)
       use probf90_module
       use global_utility_module
       IMPLICIT NONE
 
       INTEGER_T, intent(in) :: num_materials_face
       INTEGER_T, intent(in) :: nsolve
       INTEGER_T, intent(in) :: nsolveMM
       INTEGER_T, intent(in) :: dir,level
       INTEGER_T, intent(in) :: DIMDEC(bx)
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, intent(in) :: bfact
       REAL_T, intent(inout) :: bx(DIMV(bx),nsolveMM)
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

       if (num_materials_vel.ne.1) then
        print *,"num_materials_vel.ne.1"
        stop
       endif
       if ((num_materials_face.ne.1).and. &
           (num_materials_face.ne.nmat)) then
        print *,"num_materials_face invalid"
        stop
       endif
       if (nsolveMM.ne.nsolve*num_materials_face) then
        print *,"nsolveMM invalid 595"
        stop
       endif
       if ((dir.lt.0).or.(dir.ge.SDIM)) then
        print *,"dir invalid dividedx"
        stop
       endif
       call checkbound(fablo,fabhi,DIMS(bx),0,dir,33)

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir+1)
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
        do n=1,nsolveMM
         bx(D_DECL(i,j,k),n)=bx(D_DECL(i,j,k),n)/hx
        enddo

       enddo
       enddo
       enddo
 
       return
       end subroutine FORT_DIVIDEDX


       subroutine FORT_REGULARIZE_BX ( &
         num_materials_face, &
         nsolve, &
         nsolveMM, &
         nsolveMM_FACE, &
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
 
       INTEGER_T, intent(in) :: num_materials_face
       INTEGER_T, intent(in) :: nsolve
       INTEGER_T, intent(in) :: nsolveMM
       INTEGER_T, intent(in) :: nsolveMM_FACE
       INTEGER_T :: nsolveMM_FACE_test
       INTEGER_T, intent(in) :: dir
       INTEGER_T, intent(in) :: level
       INTEGER_T, intent(in) :: DIMDEC(bx)
       INTEGER_T, intent(in) :: domlo(SDIM),domhi(SDIM)
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, intent(in) :: bfact
       REAL_T, intent(in) :: min_interior_coeff
       REAL_T, intent(inout) :: bx(DIMV(bx),nsolveMM)
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
       if (num_materials_vel.ne.1) then
        print *,"num_materials_vel.ne.1"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid23"
        stop
       endif
       if (nsolveMM.ne.nsolve*num_materials_face) then
        print *,"nsolveMM invalid 690"
        stop
       endif
       nsolveMM_FACE_test=nsolveMM
       if (num_materials_face.eq.1) then
        ! do nothing
       else if (num_materials_face.eq.nmat) then
        nsolveMM_FACE_test=nsolveMM_FACE_test*2
       else
        print *,"num_materials_face invalid"
        stop
       endif
       if (nsolveMM_FACE_test.ne.nsolveMM_FACE) then
        print *,"nsolveMM_FACE invalid"
        stop
       endif

       if ((dir.lt.0).or.(dir.ge.SDIM)) then
        print *,"dir invalid mult_facewt"
        stop
       endif
       call checkbound(fablo,fabhi,DIMS(bx),0,dir,33)

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

        do n=1,nsolveMM
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
        enddo ! n=1,nsolveMM

       enddo
       enddo
       enddo
 
       return
       end subroutine FORT_REGULARIZE_BX


       subroutine FORT_MULT_FACEWT ( &
         num_materials_face, &
         nsolve, &
         nsolveMM, &
         nsolveMM_FACE, &
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
 
       INTEGER_T, intent(in) :: num_materials_face
       INTEGER_T, intent(in) :: nsolve
       INTEGER_T, intent(in) :: nsolveMM
       INTEGER_T, intent(in) :: nsolveMM_FACE
       INTEGER_T :: nsolveMM_FACE_test
       INTEGER_T, intent(in) :: dir
       INTEGER_T, intent(in) :: level
       INTEGER_T, intent(in) :: DIMDEC(bx)
       INTEGER_T, intent(in) :: DIMDEC(facewt)
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, intent(in) :: bfact
       REAL_T, intent(out) :: bx(DIMV(bx),nsolveMM)
       REAL_T, intent(in) :: facewt(DIMV(facewt),nsolveMM_FACE)
       REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

       INTEGER_T i,j,k,n
       INTEGER_T nmat
       INTEGER_T faceL,faceR

       nmat=num_materials

       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif
       if (num_materials_vel.ne.1) then
        print *,"num_materials_vel.ne.1"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid23"
        stop
       endif
       if (nsolveMM.ne.nsolve*num_materials_face) then
        print *,"nsolveMM invalid 690"
        stop
       endif
       nsolveMM_FACE_test=nsolveMM
       if (num_materials_face.eq.1) then
        ! do nothing
       else if (num_materials_face.eq.nmat) then
        nsolveMM_FACE_test=nsolveMM_FACE_test*2
       else
        print *,"num_materials_face invalid"
        stop
       endif
       if (nsolveMM_FACE_test.ne.nsolveMM_FACE) then
        print *,"nsolveMM_FACE invalid"
        stop
       endif

       if ((dir.lt.0).or.(dir.ge.SDIM)) then
        print *,"dir invalid mult_facewt"
        stop
       endif
       call checkbound(fablo,fabhi,DIMS(bx),0,dir,33)
       call checkbound(fablo,fabhi,DIMS(facewt),0,dir,33)

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        
        do n=1,nsolveMM
         if (num_materials_face.eq.1) then
          faceL=n
          faceR=n
         else if (num_materials_face.eq.nmat) then
          faceL=n
          faceR=n+nsolveMM_FACE/2
         else
          print *,"num_materials_face invalid"
          stop
         endif 
         bx(D_DECL(i,j,k),n)=half*bx(D_DECL(i,j,k),n)* &
          (facewt(D_DECL(i,j,k),faceL)+facewt(D_DECL(i,j,k),faceR))
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
       num_materials_face, &
       level, &
       finest_level, &
       nsolve, &
       nsolveMM, &
       nsolveMM_FACE, &
       nfacefrac, &
       ncellfrac, &
       nmat, &
       project_option, &
       ncphys, &
       cellmm,DIMS(cellmm), &
       xfacemm,DIMS(xfacemm), &
       yfacemm,DIMS(yfacemm), &
       zfacemm,DIMS(zfacemm), &
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

      INTEGER_T, intent(in) :: num_materials_face
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: nsolveMM
      INTEGER_T, intent(in) :: nsolveMM_FACE
      INTEGER_T             :: nsolveMM_FACE_test
      INTEGER_T, intent(in) :: nfacefrac
      INTEGER_T, intent(in) :: ncellfrac
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: DIMDEC(cellmm)
      INTEGER_T, intent(in) :: DIMDEC(xfacemm)
      INTEGER_T, intent(in) :: DIMDEC(yfacemm)
      INTEGER_T, intent(in) :: DIMDEC(zfacemm)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(masksolv)
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
      INTEGER_T, intent(in) :: bc(SDIM,2,nsolveMM)

      REAL_T, intent(in) :: cellmm(DIMV(cellmm),ncellfrac) 
      REAL_T, intent(in) :: xfacemm(DIMV(xfacemm),nfacefrac) 
      REAL_T, intent(in) :: yfacemm(DIMV(yfacemm),nfacefrac) 
      REAL_T, intent(in) :: zfacemm(DIMV(zfacemm),nfacefrac) 
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: yface(DIMV(yface),ncphys)
      REAL_T, intent(in) :: zface(DIMV(zface),ncphys)
      REAL_T, intent(out) :: masksolv(DIMV(masksolv),num_materials_face)
      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: alpha(DIMV(alpha),nsolveMM)
      REAL_T, intent(in) :: offdiagcheck(DIMV(offdiagcheck),nsolveMM)
      REAL_T, intent(out) :: maskdivres(DIMV(maskdivres))
      REAL_T, intent(out) :: maskres(DIMV(maskres))
      REAL_T, intent(in) :: mdot(DIMV(mdot),nsolveMM)
       ! coeff * areafrac * areaface / (dxfrac*dx)
      REAL_T, intent(in) :: bx(DIMV(bx),nsolveMM) 
      REAL_T, intent(in) :: by(DIMV(by),nsolveMM)
      REAL_T, intent(in) :: bz(DIMV(bz),nsolveMM)
       ! coeff * areafrac / dxfrac
      REAL_T, intent(in) :: fwtx(DIMV(fwtx),nsolveMM_FACE)  
      REAL_T, intent(in) :: fwty(DIMV(fwty),nsolveMM_FACE)
      REAL_T, intent(in) :: fwtz(DIMV(fwtz),nsolveMM_FACE)

      INTEGER_T i,j,k
      INTEGER_T inormal
      INTEGER_T ii,jj,kk
      INTEGER_T iface,jface,kface
      INTEGER_T icell,jcell,kcell
      INTEGER_T dir,side
      INTEGER_T veldir
      INTEGER_T im
      INTEGER_T faceL,faceR
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
      if ((num_materials_face.ne.1).and. &
          (num_materials_face.ne.nmat)) then
       print *,"num_materials_face invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif
      if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
       print *,"nsolve invalid24"
       stop
      endif
      if (nsolveMM.ne.nsolve*num_materials_face) then
       print *,"nsolveMM invalid 972"
       stop
      endif
      nsolveMM_FACE_test=nsolveMM
      if (num_materials_face.eq.1) then
       ! do nothing
      else if (num_materials_face.eq.nmat) then
       nsolveMM_FACE_test=nsolveMM_FACE_test*2
      else
       print *,"num_materials_face invalid"
       stop
      endif
      if (nsolveMM_FACE_test.ne.nsolveMM_FACE) then
       print *,"nsolveMM_FACE invalid"
       stop
      endif
       ! (ml,mr,2) frac_pair(ml,mr), dist_pair(ml,mr)
      if (nfacefrac.ne.nmat*nmat*2) then
       print *,"nfacefrac invalid"
       stop
      endif
       !im_inside,im_outside,3+sdim -->
       !area, dist_to_line, dist, line normal.
      if (ncellfrac.ne.nmat*nmat*(3+SDIM)) then
       print *,"ncellfrac invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid nsgenerate"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(cellmm),0,-1,234)
      call checkbound(fablo,fabhi,DIMS(xfacemm),0,0,264)
      call checkbound(fablo,fabhi,DIMS(yfacemm),0,1,264)
      call checkbound(fablo,fabhi,DIMS(zfacemm),0,SDIM-1,264)

      call checkbound(fablo,fabhi,DIMS(xface),0,0,244)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,244)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,244)
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

       do veldir=1,nsolveMM

        if (num_materials_face.eq.1) then
         if (nsolveMM.eq.nsolve) then
          faceL=veldir
          faceR=veldir
          im=1
         else
          print *,"nsolveMM invalid"
          stop
         endif
        else if (num_materials_face.eq.nmat) then
         if ((nsolve.eq.1).and.(nsolveMM.eq.nmat)) then
          faceL=veldir
          faceR=faceL+nsolveMM_FACE/2
          im=veldir
         else
          print *,"nsolve or nsolveMM invalid"
          stop
         endif
        else
         print *,"num_materials_face invalid"
         stop
        endif

        facewtsum=fwtx(D_DECL(i,j,k),faceR)+fwtx(D_DECL(i+1,j,k),faceL)+ &
                  fwty(D_DECL(i,j,k),faceR)+fwty(D_DECL(i,j+1,k),faceL)
        if (SDIM.eq.3) then
         facewtsum=facewtsum+ &
          fwtz(D_DECL(i,j,k),faceR)+fwtz(D_DECL(i,j,k+1),faceL)
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
          masksolv(D_DECL(i,j,k),im)=one
         else if (local_diag.eq.zero) then
          masksolv(D_DECL(i,j,k),im)=zero
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
          masksolv(D_DECL(i,j,k),im)=one
         else if (local_diag.eq.zero) then
          masksolv(D_DECL(i,j,k),im)=zero
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
          print *,"fwtx Left ",fwtx(D_DECL(i,j,k),faceR)
          print *,"fwtx Right ",fwtx(D_DECL(i+1,j,k),faceL)
          print *,"fwty front ",fwty(D_DECL(i,j,k),faceR)
          print *,"fwty back ",fwty(D_DECL(i,j+1,k),faceL)
          print *,"fwtz bottom ",fwtz(D_DECL(i,j,k),faceR)
          print *,"fwtz top ",fwtz(D_DECL(i,j,k+1),faceL)
          stop
         endif

        else
         print *,"facewtsum or offdiagsum invalid"
         stop
        endif

        local_diag=alpha(D_DECL(i,j,k),veldir)+offdiagsum

        local_diag_check1=facewtsum
          ! offdiagcheck is initialized in BUILDFACEWT
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

       enddo ! veldir=1..nsolveMM
          
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
       num_materials_face, &
       level, &
       finest_level, &
       nsolve, &
       nsolveMM, &
       nsolveMM_FACE, &
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

      INTEGER_T, intent(in) :: num_materials_face
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: nsolveMM
      INTEGER_T, intent(in) :: nsolveMM_FACE
      INTEGER_T             :: nsolveMM_FACE_test
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

      REAL_T, intent(in) :: alpha(DIMV(alpha),nsolveMM)
      REAL_T, intent(out) :: diag_reg(DIMV(diag_reg),nsolveMM)
      ! coeff * areafrac * areaface / (dxfrac*dx)  (if coeff>0.0)
      REAL_T, intent(in) :: bx(DIMV(bx),nsolveMM) 
      REAL_T, intent(in) :: by(DIMV(by),nsolveMM)
      REAL_T, intent(in) :: bz(DIMV(bz),nsolveMM)

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
      if ((num_materials_face.ne.1).and. &
          (num_materials_face.ne.nmat)) then
       print *,"num_materials_face invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif
      if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
       print *,"nsolve invalid24"
       stop
      endif
      if (nsolveMM.ne.nsolve*num_materials_face) then
       print *,"nsolveMM invalid 972"
       stop
      endif
      nsolveMM_FACE_test=nsolveMM
      if (num_materials_face.eq.1) then
       ! do nothing
      else if (num_materials_face.eq.nmat) then
       nsolveMM_FACE_test=nsolveMM_FACE_test*2
      else
       print *,"num_materials_face invalid"
       stop
      endif
      if (nsolveMM_FACE_test.ne.nsolveMM_FACE) then
       print *,"nsolveMM_FACE invalid"
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

       do veldir=1,nsolveMM

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

