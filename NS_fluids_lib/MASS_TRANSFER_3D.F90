#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 0

#define DEBUG_TRIPLE 0
#define DEBUG_ACTIVE_CELL 0
#define DEBUG_I 22
#define DEBUG_J 13
#define DEBUG_K 0

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "MASS_TRANSFER_F.H"


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

      module mass_transfer_module
      use probcommon_module

      ! if t1 is a c++ parameter and p1 is a "type" component:
      ! REAL_T, pointer :: p1(DIMV(t1),ncomp)
      ! REAL_T, target :: t1(DIMV(t1),ncomp)
      ! p1=>t1
      type probe_parm_type
       INTEGER_T :: tid
       INTEGER_T :: use_supermesh
       REAL_T, pointer :: Y_TOLERANCE
       INTEGER_T, pointer :: local_freezing_model
       REAL_T, pointer :: LL
       INTEGER_T, pointer :: debugrate
       INTEGER_T, pointer :: i,j,k
       REAL_T, pointer :: xsrc(:)
       REAL_T, pointer :: xdst(:)
       REAL_T, pointer :: xsrc_micro(:)
       REAL_T, pointer :: xdst_micro(:)
       INTEGER_T, pointer :: im_source
       INTEGER_T, pointer :: im_dest
       INTEGER_T, pointer :: tcomp_source
       INTEGER_T, pointer :: Ycomp_source
       INTEGER_T, pointer :: dencomp_source
       REAL_T, pointer :: dxprobe_source
       INTEGER_T, pointer :: tcomp_dest
       INTEGER_T, pointer :: Ycomp_dest
       INTEGER_T, pointer :: dencomp_dest
       REAL_T, pointer :: dxprobe_dest
       REAL_T, pointer, dimension(:) :: LSINT
       INTEGER_T, pointer :: imls_I
       REAL_T, pointer :: dxmaxLS
       INTEGER_T, pointer :: bfact
       INTEGER_T, pointer :: level
       INTEGER_T, pointer :: finest_level
       REAL_T, pointer :: dx(:)
       REAL_T, pointer :: xlo(:)
       REAL_T, pointer :: xI(:)
       INTEGER_T, pointer :: nmat
       INTEGER_T, pointer :: ngrow
       INTEGER_T, pointer :: fablo(:)
       INTEGER_T, pointer :: fabhi(:)
       INTEGER_T :: DIMDEC(EOS)
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: EOS
       INTEGER_T :: DIMDEC(recon)
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: recon
       INTEGER_T :: DIMDEC(LS)
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LS
       INTEGER_T :: DIMDEC(pres)
       REAL_T, pointer, dimension(D_DECL(:,:,:)) :: pres
      end type probe_parm_type

      type TSAT_MASS_FRAC_parm_type
       type(probe_parm_type), pointer :: PROBE_PARMS
       INTEGER_T :: Tanasawa_or_Schrage_or_Kassemi
       REAL_T :: accommodation_coefficient
       REAL_T :: reference_pressure
       REAL_T :: universal_gas_constant_R
       REAL_T :: molar_mass_ambient   
       REAL_T :: molar_mass_vapor
       REAL_T :: TSAT_base
       REAL_T :: YI_min
       REAL_T :: TI_min
       REAL_T :: TI_max
       REAL_T :: D_MASS
       REAL_T :: den_G
       INTEGER_T :: iprobe_vapor
       INTEGER_T, pointer :: material_type_evap(:)
       REAL_T, pointer :: thermal_k(:)
      end type TSAT_MASS_FRAC_parm_type

      contains

      subroutine get_interface_temperature( &
        use_tsatfab, &
        i,j,k, &
        ireverse, &
        iten, &
        ntsat, &
        bfact, &
        level, &
        finest_level, &
        dx,xlo, &
        ngrow_tsat, &
        fablo,fabhi, &
        TSATFAB,DIMS(TSATFAB), &
        T_I, &
        latent_comp, &
        TSAT_array, &
        TSAT_flag_array, &  ! =0 if derived, =1 if T_interface(x,y,z,t) given.
        x_I,time,nmat,nten,caller_id)
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: use_tsatfab
      INTEGER_T, intent(in) :: i,j,k
      INTEGER_T, intent(in) :: ireverse
      INTEGER_T, intent(in) :: iten
      INTEGER_T, intent(in) :: ntsat
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: ngrow_tsat
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(TSATFAB)
      REAL_T, intent(in) :: TSATFAB(DIMV(TSATFAB),ntsat)
      INTEGER_T, intent(in) :: caller_id
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: latent_comp
      REAL_T, intent(out) :: T_I
      REAL_T, intent(in)  :: TSAT_array(2*nten)
      INTEGER_T, intent(in)  :: TSAT_flag_array(2*nten)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: x_I(SDIM)
      REAL_T, intent(in) :: time
      INTEGER_T :: local_flag
      INTEGER_T :: dir
      INTEGER_T :: ncomp_per_tsat
      INTEGER_T :: tsat_comp

#if (STANDALONE==1)
      REAL_T, external   :: exact_temperature
#endif

      call checkbound(fablo,fabhi,DIMS(TSATFAB),ngrow_tsat,-1,122)
      ncomp_per_tsat=2
      tsat_comp=nten+(iten-1)*ncomp_per_tsat+1

      if (ntsat.eq.nten*(1+ncomp_per_tsat)) then
       ! do nothing
      else
       print *,"ntsat invalid"
       stop
      endif
      if (bfact.lt.1) then 
       print *,"bfact invalid114"
       stop
      endif
      if ((iten.ge.1).and.(iten.le.nten)) then
       ! do nothing
      else
       print *,"iten invalid"
       stop
      endif

      if (nmat.eq.num_materials) then
       if (nten.eq.(((nmat-1)*(nmat-1)+nmat-1)/2)) then 
        if ((latent_comp.ge.1).and.(latent_comp.le.2*nten)) then
         if (time.ge.zero) then
          if (use_tsatfab.eq.1) then
           call interpfab_tsat( &
             i,j,k, &
             ireverse, &
             iten, &
             nten, &
             ntsat, &
             bfact, &
             level, &
             finest_level, &
             dx,xlo, &
             x_I, &
             tsat_comp, &
             ngrow_tsat, &
             fablo,fabhi, &
             TSATFAB,DIMS(TSATFAB), &
             T_I)
          else if (use_tsatfab.eq.0) then
           local_flag=TSAT_flag_array(latent_comp)
           if (local_flag.eq.0) then
            T_I=TSAT_array(latent_comp)
           else if ((local_flag.ge.1).and.(local_flag.le.nmat)) then
#if (STANDALONE==0)
            print *,"local_flag <> 0 not supported yet"
            stop
#elif (STANDALONE==1)
            T_I=exact_temperature(x_I,time,local_flag,probtype,nmat, &
                   fort_heatviscconst)
#endif
           else
            print *,"local_flag invalid"
            stop
           endif
          else
           print *,"use_tsatfab invalid"
           stop
          endif

          if (T_I.ge.TEMPERATURE_FLOOR) then
           ! do nothing
          else
           print *,"T_I out of range T_I=",T_I
           print *,"caller_id=",caller_id
           print *,"latent_comp=",latent_comp
           print *,"local_flag=",local_flag
           do dir=1,SDIM
            print *,"dir,x_I ",dir,x_I(dir)
           enddo
           stop
          endif
         else
          print *,"time invalid"
          stop
         endif
        else
         print *,"latent_comp invalid"
         stop
        endif
       else
        print *,"nten invalid"
        stop
       endif
      else
       print *,"nmat invalid"
       stop
      endif

      return
      end subroutine get_interface_temperature

      subroutine adjust_du(du,normdir,rval,map_forward)
      IMPLICIT NONE

      REAL_T du,rval,disc
      INTEGER_T normdir,map_forward


      if ((normdir.lt.0).or.(normdir.ge.SDIM)) then
       print *,"normdir invalid"
       stop
      endif

        ! if map_forward=0,
        ! dr/dt=-rf uf/r
        ! r^2/2=-rf uf t + C
        ! r^2/2=-rf uf t +rf^2/2
        ! r=sqrt(rf^2-2 rf du)  
        ! r=rf-new du
        ! new du=rf-sqrt(rf^2-2 rf du) =2rf du/(rf+sqrt(rf^2-2 rf du))=
        ! du/(1/2+sqrt(1/4-du/(2rf))
      if ((levelrz.eq.1).or. &
          (levelrz.eq.3)) then
       if (normdir.eq.0) then
        if (rval.le.zero) then
         du=zero
        else
         if (map_forward.eq.0) then
          disc=one/four-du/(two*rval)
         else if (map_forward.eq.1) then
          disc=one/four+du/(two*rval)
         else
          print *,"map_forward invalid"
          stop
         endif

         if (disc.gt.zero) then
          du=du/(half+sqrt(disc))
         endif
        endif
       else if ((normdir.eq.1).or.(normdir.eq.SDIM-1)) then
        ! do nothing
       else
        print *,"normdir invalid"
        stop
       endif
      else if (levelrz.eq.0) then
       ! do nothing
      else 
       print *,"levelrz invalid adjust du"
       stop
      endif
 
      return
      end subroutine adjust_du

       ! cc_flag=0  centroid -> center or some other target "xtarget"
       ! cc_flag=1  center -> centroid
       ! tsat_flag=-1 use all cells in the stencil
       ! tsat_flag=0 no tsat
       ! tsat_flag=1 tsat
       ! if normal probe is target: xtarget = xprobe
       ! if centroid to center: xtarget = xcenter
       ! if center to centroid: xtarget = xcentroid
      subroutine center_centroid_interchange( &
       DATA_FLOOR, &
       nsolve, &
       cc_flag, &
       tsat_flag, &
       bfact, &
       level, &
       finest_level, &
       dx,xlo, & 
       xsten,nhalf, &
       T_sten, &
       XC_sten, &
       xI, &  ! closest point on interface to cell center.
       xtarget, &
       VF_sten, &
       LS_sten, &
       TSAT, & ! unused if tsat_flag==-1, 0 (but set to 293.0 for sanity check)
       T_out)
      use global_utility_module
      implicit none

      REAL_T, intent(in) :: DATA_FLOOR
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T :: nc
      INTEGER_T, intent(in) :: cc_flag
      INTEGER_T, intent(in) :: tsat_flag
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: xI(SDIM)
      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
      REAL_T, intent(in) :: TSAT
      REAL_T, intent(inout) :: T_sten(D_DECL(-1:1,-1:1,-1:1),nsolve)
      REAL_T, intent(in) :: XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T, intent(in) :: VF_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T, intent(in) :: LS_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T :: wt_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T, intent(out) :: T_out(nsolve)

      INTEGER_T i,j,k
      INTEGER_T klosten,khisten
      INTEGER_T i1,j1,k1
      INTEGER_T dir
      REAL_T dist
      REAL_T, dimension(:,:), allocatable :: AA
      REAL_T, dimension(:,:), allocatable :: AAcopy
      REAL_T BB(SDIM+1)
      REAL_T BBcopy(SDIM+1)
      REAL_T xtemp(SDIM)
      REAL_T, intent(in) :: xtarget(SDIM)
      REAL_T xbase(SDIM)
      REAL_T xlive(SDIM)
      REAL_T delx(SDIM+1)
      REAL_T GRADTEMP(SDIM+1)
      REAL_T TMIN,TMAX,T_test,T_avg,wt_sum,VF,LS
      INTEGER_T mat_ncomp,matstatus
      REAL_T T_hold
      INTEGER_T own_flag
      REAL_T wt_local

      if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
       print *,"nsolve invalid"
       stop
      endif
      if ((cc_flag.ne.0).and.(cc_flag.ne.1)) then
       print *,"cc_flag invalid"
       stop
      endif
      if ((tsat_flag.ne.-1).and. &
          (tsat_flag.ne.0).and. &
          (tsat_flag.ne.1)) then
       print *,"tsat_flag invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid112"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid center_centroid_interchange"
       stop
      endif
      if (nhalf.ne.3) then
       print *,"nhalf invalid"
       stop
      endif
       ! was TEMPERATURE_FLOOR
      if (TSAT.lt.DATA_FLOOR) then
       print *,"TSAT invalid in center_centroid_interchange"
       print *,"TSAT=",TSAT
       stop
      endif

      if (SDIM.eq.2) then
       klosten=0
       khisten=0
      else if (SDIM.eq.3) then
       klosten=-1
       khisten=1
      else
       print *,"dimension bust"
       stop
      endif

       !   T(x)=T(xbase) + grad T dot (x-xbase)
      do dir=1,SDIM
       if (dx(dir).le.zero) then
        print *,"dx(dir) invalid"
        stop
       endif
       if (tsat_flag.eq.1) then
        xbase(dir)=xI(dir)
       else if (tsat_flag.eq.0) then
        xbase(dir)=xtarget(dir)
       else if (tsat_flag.eq.-1) then
        xbase(dir)=xtarget(dir)
       else
        print *,"tsat_flag invalid"
        stop
       endif
      enddo ! dir=1..sdim

      do nc=1,nsolve

       TMAX=T_sten(D_DECL(0,0,0),nc)
       if (tsat_flag.eq.0) then
        ! do nothing
       else if (tsat_flag.eq.1) then
        TMAX=TSAT
       else if (tsat_flag.eq.-1) then
        ! do nothing
       else
        print *,"tsat_flag invalid"
        stop
       endif
       TMIN=TMAX

       T_avg=zero
       wt_sum=zero
       
       do i=-1,1
       do j=-1,1
       do k=klosten,khisten

        xtemp(1)=xsten(2*i,1)
        xtemp(2)=xsten(2*j,2)
        if (SDIM.eq.3) then
         xtemp(SDIM)=xsten(2*k,SDIM)
        endif
        do dir=1,SDIM
         if (cc_flag.eq.0) then ! centroid-> center or centroid to probe.
          xlive(dir)=XC_sten(D_DECL(i,j,k),dir) ! cell centroid
         else if (cc_flag.eq.1) then ! center-> centroid
          xlive(dir)=xtemp(dir)  ! cell center
         else
          print *,"cc_flag invalid"
          stop
         endif
        enddo ! dir=1..sdim

        if ((tsat_flag.eq.0).or. & ! no tsat
            (tsat_flag.eq.1)) then ! tsat
         VF=VF_sten(D_DECL(i,j,k))
         LS=LS_sten(D_DECL(i,j,k))
        else if (tsat_flag.eq.-1) then ! use all cells in the stencil
         VF=zero
         LS=zero
        else
         print *,"tsat_flag invalid"
         stop
        endif

        own_flag=0
        if (tsat_flag.eq.-1) then
         own_flag=1
        else if (cc_flag.eq.0) then ! centroid->center
         if (abs(VF).le.VOFTOL) then
          own_flag=0
         else if ((VF.ge.VOFTOL).and.(VF.le.one+VOFTOL)) then
          own_flag=1
         else
          print *,"VF invalid"
          stop
         endif
        else if (cc_flag.eq.1) then ! center->centroid
         if ((LS.le.zero).or.(abs(VF).le.VOFTOL)) then
          own_flag=0
         else if ((LS.ge.zero).and.(VF.ge.VOFTOL)) then
          own_flag=1
         else
          print *,"LS or VF invalid"
          stop
         endif
        else
         print *,"cc_flag invalid"
         stop
        endif
 
        wt_local=zero
 
        if (own_flag.eq.0) then
         if (tsat_flag.eq.1) then
          T_sten(D_DECL(i,j,k),nc)=TSAT
         else if (tsat_flag.eq.0) then
          ! do nothing
         else
          print *,"tsat_flag invalid"
          stop
         endif
         wt_local=VOFTOL
        else if (own_flag.eq.1) then
          ! if normal probe is target: xtarget = xprobe
          ! if new supermesh centroid is target: xtarget = x_new_centroid
          ! if centroid to center: xtarget = xcenter
          ! if center to centroid: xtarget = xcentroid
         dist=VOFTOL
         do dir=1,SDIM
          dist=dist+(xlive(dir)-xtarget(dir))**2/(dx(1)**2)
         enddo ! dir
         if (dist.le.zero) then
          print *,"dist invalid"
          stop
         endif
         wt_local=one/dist
         if (tsat_flag.eq.-1) then ! use all cells in the stencil
          ! do nothing
         else if ((tsat_flag.eq.1).or. &
                  (tsat_flag.eq.0)) then
          if (VF.le.zero) then
           print *,"VF invalid"
           stop
          endif
          wt_local=wt_local*(VF**2)
         else
          print *,"tsat_flag invalid"
          stop
         endif

        else
         print *,"own_flag invalid"
         stop
        endif

        if (wt_local.le.zero) then
         print *,"wt_local invalid"
         stop
        endif

        wt_sten(D_DECL(i,j,k))=wt_local

        T_test=T_sten(D_DECL(i,j,k),nc)
        if (T_test.gt.TMAX) then
         TMAX=T_test
        endif
        if (T_test.lt.TMIN) then
         TMIN=T_Test
        endif
        if (abs(T_test).lt.1.0D+50) then
         ! do nothing
        else
         print *,"T_test bust"
         print *,"i,j,k,nc,T_test,TMIN,TMAX ",i,j,k,nc,T_test,TMIN,TMAX
         stop
        endif

        wt_sum=wt_sum+wt_local
        T_avg=T_avg+wt_local*T_test 

       enddo
       enddo
       enddo  ! i,j,k

       if (wt_sum.le.zero) then
        print *,"wt_sum invalid"
        stop
       endif

       do i=-1,1
       do j=-1,1
       do k=klosten,khisten
        wt_sten(D_DECL(i,j,k))=wt_sten(D_DECL(i,j,k))/wt_sum
       enddo
       enddo
       enddo

       if (tsat_flag.eq.1) then
        if (T_avg.lt.DATA_FLOOR) then
         print *,"T_avg invalid"
         print *,"T_avg=",T_avg
         stop
        endif
        T_avg=TSAT
       else if ((tsat_flag.eq.0).or.(tsat_flag.eq.-1)) then
        T_avg=T_avg/wt_sum
       else
        print *,"tsat_flag invalid"
        stop
       endif 

       if ((tsat_flag.eq.0).or.(tsat_flag.eq.-1)) then
        mat_ncomp=SDIM+1
       else if (tsat_flag.eq.1) then
        mat_ncomp=SDIM
       else
        print *,"tsat_flag invalid"
        stop
       endif

       allocate(AA(mat_ncomp,mat_ncomp))
       allocate(AAcopy(mat_ncomp,mat_ncomp))

       do i=1,mat_ncomp
       do j=1,mat_ncomp
        AA(i,j)=zero
       enddo
       enddo
       do i=1,SDIM+1
        BB(i)=zero
       enddo

       do i=1,mat_ncomp
       do j=1,mat_ncomp
        do i1=-1,1
        do j1=-1,1
        do k1=klosten,khisten 
         xtemp(1)=xsten(2*i1,1)
         xtemp(2)=xsten(2*j1,2)
         if (SDIM.eq.3) then
          xtemp(SDIM)=xsten(2*k1,SDIM)
         endif
         do dir=1,SDIM
          if (cc_flag.eq.0) then ! centroid-> center or probe
           xlive(dir)=XC_sten(D_DECL(i1,j1,k1),dir) ! centroid
          else if (cc_flag.eq.1) then ! center-> centroid
           xlive(dir)=xtemp(dir) ! center
          else
           print *,"cc_flag invalid"
           stop
          endif
         enddo ! dir

         do dir=1,SDIM
          delx(dir)=(xlive(dir)-xbase(dir))/dx(dir)
         enddo ! dir
         delx(SDIM+1)=one

         AA(i,j)=AA(i,j)+wt_sten(D_DECL(i1,j1,k1))*delx(i)*delx(j)
         if (j.eq.1) then
          BB(i)=BB(i)+wt_sten(D_DECL(i1,j1,k1))*delx(i)* &
           (T_sten(D_DECL(i1,j1,k1),nc)-T_avg)
         endif
        enddo ! k1
        enddo ! j1
        enddo ! i1
        AAcopy(i,j)=AA(i,j)
        if (j.eq.1) then
         BBcopy(i)=BB(i)
        endif
       enddo   ! j
       enddo   ! i
 
       call matrix_solve(AA,GRADTEMP,BB,matstatus,mat_ncomp)

       if (matstatus.eq.1) then
        T_hold=T_avg

        do dir=1,SDIM
         delx(dir)=(xtarget(dir)-xbase(dir))/dx(dir)
        enddo ! dir
        delx(SDIM+1)=one

        do dir=1,mat_ncomp
         T_hold=T_hold+GRADTEMP(dir)*delx(dir)
        enddo
        if (T_hold.lt.TMIN) then
         T_hold=TMIN
        endif
        if (T_hold.gt.TMAX) then
         T_hold=TMAX
        endif
        if (tsat_flag.eq.1) then
         if (T_hold.lt.zero) then
          print *,"temperature underflow in center_centroid_interchange"
          print *,"T_hold=",T_hold
          stop
         endif
        else if ((tsat_flag.eq.0).or.(tsat_flag.eq.-1)) then
         ! do nothing
        else
         print *,"tsat_flag invalid"
         stop
        endif
       else if (matstatus.eq.0) then
        print *,"WARNING: subroutine center_centroid_interchange"
        print *,"WARNING: linear least squares failed"
        print *,"WARNING: doing piecewise constant interpolation instead"
        print *,"nsolve=",nsolve
        print *,"cc_flag=",cc_flag
        print *,"tsat_flag=",tsat_flag
        print *,"TSAT=",TSAT
        print *,"T_avg=",T_avg
        print *,"wt_sum=",wt_sum
        print *,"mat_ncomp=",mat_ncomp
        print *,"xsten(0,dir) ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
        print *,"xtarget ",xtarget(1),xtarget(2),xtarget(SDIM)
        print *,"dx ",dx(1),dx(2),dx(SDIM)
        print *,"xI ",xI(1),xI(2),xI(SDIM)
        do i=1,mat_ncomp
        do j=1,mat_ncomp
         print *,"i,j,AAcopy,AA ",i,j,AAcopy(i,j),AA(i,j)
        enddo
        enddo
        do i=1,mat_ncomp
         print *,"i,BBcopy,BB ",i,BBcopy(i),BB(i)
        enddo
        
        T_hold=T_avg
       else
        print *,"matstatus has an invalid value"
        stop
       endif
 
       T_out(nc)=T_hold

       deallocate(AA)
       deallocate(AAcopy)

      enddo ! nc

      return
      end subroutine center_centroid_interchange

      subroutine interpfabFWEIGHT( &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo,x, &
       im,nmat,&
       comp, &
       ngrow, &
       lo,hi, &
       data,DIMS(data), &
       recon,DIMS(recon), &
       dest)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: x(SDIM)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T, intent(in) :: comp,ngrow
      INTEGER_T, intent(in) :: im,nmat
       ! datalox,datahix,dataloy,datahiy,dataloz,datahiz
      INTEGER_T, intent(in) :: DIMDEC(data)
      INTEGER_T, intent(in) :: DIMDEC(recon)
       ! datalox:datahix,dataloy:datahiy,dataloz:datahiz
      REAL_T, intent(in) :: data(DIMV(data),comp)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(out) :: dest

      REAL_T :: DATA_FLOOR

      REAL_T :: T_out(1)

      INTEGER_T dir
      INTEGER_T ic,jc,kc
      INTEGER_T i1,j1,k1
      INTEGER_T isten,jsten,ksten
      INTEGER_T k1lo,k1hi
      INTEGER_T nhalf
      INTEGER_T vofcomp
      INTEGER_T cell_index(SDIM)

      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_stencil(-3:3,SDIM)
      REAL_T volcell
      REAL_T cencell(SDIM)
      REAL_T T_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T VF_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      INTEGER_T cc_flag
      INTEGER_T tsat_flag
      INTEGER_T nsolve
      REAL_T Tsat

      DATA_FLOOR=zero

      call checkbound(lo,hi,DIMS(data),ngrow,-1,1221)
      call checkbound(lo,hi,DIMS(recon),ngrow,-1,1222)

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid20"
       stop
      endif
      if (bfact.lt.1) then 
       print *,"bfact invalid113"
       stop
      endif
      if ((comp.lt.1).or.(comp.gt.1000)) then
       print *,"comp out of range"
       stop
      endif

      call containing_cell(bfact,dx,xlo,lo,x,cell_index)

      do dir=1,SDIM
       if (cell_index(dir).lt.lo(dir)-ngrow+1) then
        cell_index(dir)=lo(dir)-ngrow+1
       endif
       if (cell_index(dir).gt.hi(dir)+ngrow-1) then
        cell_index(dir)=hi(dir)+ngrow-1
       endif
      enddo ! dir

      vofcomp=(im-1)*ngeom_recon+1

      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      nhalf=3
      call gridsten_level(xsten,ic,jc,kc,level,nhalf)

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+ic
       jsten=j1+jc
       ksten=k1+kc

       call gridsten_level(xsten_stencil,isten,jsten,ksten,level,nhalf)
       call Box_volumeFAST(bfact,dx,xsten_stencil,nhalf, &
         volcell,cencell,SDIM)

       T_sten(D_DECL(i1,j1,k1))=data(D_DECL(isten,jsten,ksten),comp)
       VF_sten(D_DECL(i1,j1,k1))=recon(D_DECL(isten,jsten,ksten),vofcomp)
       do dir=1,SDIM
        XC_sten(D_DECL(i1,j1,k1),dir)= &
         recon(D_DECL(isten,jsten,ksten),vofcomp+dir)+cencell(dir)
       enddo

      enddo
      enddo
      enddo ! i1,j1,k1

       ! in: interpfabFWEIGHT
      cc_flag=0  ! centroid -> target
      tsat_flag=0 ! do not use TSAT
      nsolve=1
      Tsat=293.0 ! representative value for sanity check
      call center_centroid_interchange( &
       DATA_FLOOR, &
       nsolve, &
       cc_flag, &
       tsat_flag, &
       bfact, &
       level, &
       finest_level, &
       dx,xlo, &
       xsten,nhalf, &
       T_sten, & 
       XC_sten, & 
       x, &   ! xI not used
       x, &   ! xtarget
       VF_sten, & 
       VF_sten, & 
       Tsat, & ! unused if tsat_flag==0
       T_out)

      dest=T_out(1)

      return 
      end subroutine interpfabFWEIGHT

      subroutine interpfabVFRAC_tess( &
       tid, &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo,x, &
       nmat,&
       ngrow, &
       lo,hi, &
       recon,DIMS(recon), & ! fluids tess, solids overlay
       dest)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: x(SDIM)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: DIMDEC(recon)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(out) :: dest(nmat)

      INTEGER_T im
      INTEGER_T dir
      INTEGER_T ic,jc,kc
      INTEGER_T nhalf
      INTEGER_T vofcomp
      INTEGER_T cell_index(SDIM)

      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T xsten(-3:3,SDIM)
      REAL_T volcell
      REAL_T cencell(SDIM)
      INTEGER_T nmax
      INTEGER_T local_tessellate

      call checkbound(lo,hi,DIMS(recon),ngrow,-1,1222)

      if (bfact.lt.1) then 
       print *,"bfact invalid113"
       stop
      endif

      nmax=POLYGON_LIST_MAX 
      if ((nmax.lt.100).or.(nmax.gt.2000)) then
       print *,"nmax invalid"
       stop
      endif

      call containing_cell(bfact,dx,xlo,lo,x,cell_index)

      do dir=1,SDIM
       if (cell_index(dir).lt.lo(dir)-ngrow+1) then
        cell_index(dir)=lo(dir)-ngrow+1
       endif
       if (cell_index(dir).gt.hi(dir)+ngrow-1) then
        cell_index(dir)=hi(dir)+ngrow-1
       endif
      enddo ! dir

      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      nhalf=3
      call gridsten_level(xsten,ic,jc,kc,level,nhalf)
      call Box_volumeFAST(bfact,dx,xsten,nhalf, &
        volcell,cencell,SDIM)

      do im=1,nmat*ngeom_recon
       mofdata(im)=recon(D_DECL(ic,jc,kc),im)
      enddo

      local_tessellate=3
      call multi_get_volume_tessellate( &
        local_tessellate, & ! =3
        bfact, &
        dx, &
        xsten,nhalf, &
        mofdata, &
        geom_xtetlist(1,1,1,tid+1), &
        nmax, &
        nmax, &
        nmat, &
        SDIM, &
        4)

      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       dest(im)=mofdata(vofcomp)
      enddo

      return 
      end subroutine interpfabVFRAC_tess


      subroutine grad_probe_sanity(xI,xprobe,temp_probe,Tsat,LL)
      IMPLICIT NONE

      REAL_T, intent(in) :: LL
      REAL_T, intent(in) :: xI(SDIM)
      REAL_T, intent(in) :: xprobe(SDIM)
      REAL_T, intent(in) :: temp_probe
      REAL_T, intent(in) :: Tsat

      REAL_T :: grad_probe_local
      REAL_T :: mag
      REAL_T :: nrm(SDIM)

      INTEGER_T dir

       ! nrm points from xI to xprobe

      if (LL.ne.zero) then

       mag=zero
       do dir=1,SDIM
        nrm(dir)=xprobe(dir)-xI(dir)
        mag=mag+nrm(dir)**2
       enddo
       mag=sqrt(mag)
       if (mag.le.zero) then
        print *,"mag invalid 1"
        stop
       endif
       do dir=1,SDIM
        nrm(dir)=nrm(dir)/mag
       enddo
       grad_probe_local=(temp_probe-Tsat)/mag

       if (grad_probe_local/LL.ge.zero) then
        ! do nothing (expected case for temperature gradient)
       else if (grad_probe_local/LL.lt.zero) then
        ! do nothing 
       else
        print *,"grad_probe_local bust"
        stop
       endif

      else
       print *,"LL invalid"
       stop
      endif

      return
      end subroutine grad_probe_sanity

      subroutine interpfab( &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       comp, &
       ngrow, &
       lo,hi, &
       data, &
       DIMS(data), &
       dest)
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xtarget(SDIM)
      INTEGER_T, intent(in) :: lo(SDIM)
      INTEGER_T, intent(in) :: hi(SDIM)
      INTEGER_T, intent(in) :: comp
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: DIMDEC(data)
      REAL_T, intent(in) :: data(DIMV(data),comp)
      REAL_T, intent(out) :: dest

      REAL_T :: DATA_FLOOR

      REAL_T :: T_out(1)

      INTEGER_T dir
      INTEGER_T ic,jc,kc
      INTEGER_T i1,j1,k1
      INTEGER_T isten,jsten,ksten
      INTEGER_T k1lo,k1hi
      INTEGER_T nhalf
      INTEGER_T cell_index(SDIM)

      REAL_T xsten(-3:3,SDIM)
      REAL_T T_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      INTEGER_T cc_flag
      INTEGER_T tsat_flag
      INTEGER_T nsolve
      REAL_T Tsat

      DATA_FLOOR=zero

      call checkbound(lo,hi,DIMS(data),ngrow,-1,122)

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      if (bfact.lt.1) then 
       print *,"bfact invalid114"
       stop
      endif
      if ((comp.lt.1).or.(comp.gt.1000)) then
       print *,"comp out of range"
       stop
      endif

      call containing_cell(bfact,dx,xlo,lo,xtarget,cell_index)

      do dir=1,SDIM
       if (cell_index(dir).lt.lo(dir)-ngrow+1) then
        cell_index(dir)=lo(dir)-ngrow+1
       endif
       if (cell_index(dir).gt.hi(dir)+ngrow-1) then
        cell_index(dir)=hi(dir)+ngrow-1
       endif
      enddo ! dir
      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      nhalf=3
      call gridsten_level(xsten,ic,jc,kc,level,nhalf)

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+ic
       jsten=j1+jc
       ksten=k1+kc

       T_sten(D_DECL(i1,j1,k1))=data(D_DECL(isten,jsten,ksten),comp)

      enddo
      enddo
      enddo ! i1,j1,k1

       ! in: interpfab
      cc_flag=1  ! center -> target
      tsat_flag=-1  ! use all cells in the stencil
      nsolve=1
      Tsat=293.0 !unused if tsat_flag==-1 (but set to 293.0 for sanity check)
      call center_centroid_interchange( &
       DATA_FLOOR, &
       nsolve, &
       cc_flag, &
       tsat_flag, &
       bfact, &
       level, &
       finest_level, &
       dx,xlo, &
       xsten,nhalf, &
       T_sten, & 
       XC_sten, & 
       xtarget, & ! xI (not used)
       xtarget, & ! xtarget
       T_sten, &  ! VF_sten (not used)
       T_sten, &  ! LS_Sten (not used)
       Tsat, & ! unused if tsat_flag==-1
       T_out)

      dest=T_out(1)

      return 
      end subroutine interpfab

      subroutine interpfab_tsat( &
       i,j,k, &
       ireverse, &
       iten, &
       nten, &
       ntsat, &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       comp, &
       ngrow, &
       lo,hi, &
       data,DIMS(data), &
       TSAT)
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: i,j,k
      INTEGER_T, intent(in) :: ireverse
      INTEGER_T, intent(in) :: iten
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: ntsat
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xtarget(SDIM)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T, intent(in) :: comp,ngrow
      INTEGER_T, intent(in) :: DIMDEC(data)
      REAL_T, intent(in) :: data(DIMV(data),ntsat)
      REAL_T, intent(out) :: TSAT

      INTEGER_T ncomp_per_tsat
      INTEGER_T k1lo,k1hi
      INTEGER_T cell_index(3)
      INTEGER_T dir
      INTEGER_T nhalf
      REAL_T TSAT_times_weight
      REAL_T TSAT_weight
      INTEGER_T i1,j1,k1
      INTEGER_T isten,jsten,ksten
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T TSAT_FLAG
      REAL_T local_TSAT
      REAL_T local_weight
      REAL_T eps

      call checkbound(lo,hi,DIMS(data),ngrow,-1,122)

      ncomp_per_tsat=2
      if (ntsat.eq.nten*(1+ncomp_per_tsat)) then
       ! do nothing
      else
       print *,"ntsat invalid"
       stop
      endif

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      if (bfact.lt.1) then 
       print *,"bfact invalid114"
       stop
      endif
      if ((comp.gt.nten).and.(comp.le.ntsat)) then
       ! do nothing
      else
       print *,"comp out of range"
       stop
      endif
      if ((iten.ge.1).and.(iten.le.nten)) then
       ! do nothing
      else
       print *,"iten invalid"
       stop
      endif

      cell_index(1)=i
      cell_index(2)=j
      cell_index(3)=k

      do dir=1,SDIM
       if ((cell_index(dir).ge.lo(dir)).and. &
           (cell_index(dir).le.hi(dir))) then
        ! do nothing
       else
        print *,"cell_index out of range"
        stop
       endif
      enddo ! dir

      nhalf=3

      TSAT_times_weight=zero
      TSAT_weight=zero

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+i
       jsten=j1+j
       ksten=k1+k
       call gridsten_level(xsten,isten,jsten,ksten,level,nhalf)

       TSAT_FLAG=NINT(data(D_DECL(isten,jsten,ksten),iten))
       if (ireverse.eq.0) then
        ! do nothing
       else if (ireverse.eq.1) then
        TSAT_FLAG=-TSAT_FLAG
       else
        print *,"ireverse invalid"
        stop
       endif
       if ((TSAT_FLAG.eq.1).or.(TSAT_FLAG.eq.2)) then
        local_TSAT=data(D_DECL(isten,jsten,ksten),comp)
        local_weight=zero
        eps=dx(1)*0.001
        do dir=1,SDIM
         local_weight=local_weight+(xsten(0,dir)-xtarget(dir))**2
        enddo 
        local_weight=local_weight+eps**2
        local_weight=one/local_weight
        TSAT_weight=TSAT_weight+local_weight
        TSAT_times_weight=TSAT_times_weight+local_weight*local_TSAT
       else if ((TSAT_FLAG.eq.-1).or.(TSAT_FLAG.eq.-2)) then
        ! do nothing
       else if (TSAT_FLAG.eq.0) then
        ! do nothing
       else
        print *,"TSAT_FLAG invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i1,j1,k1

      if (TSAT_weight.gt.zero) then
       if (TSAT_times_weight.gt.zero) then
        TSAT=TSAT_times_weight/TSAT_weight
       else
        print *,"TSAT_times_weight invalid"
        stop
       endif
      else 
       print *,"TSAT_weight invalid"
       print *,"VOFTOL=",VOFTOL
       print *,"VOFTOL_REDIST=",VOFTOL_REDIST
       print *,"i,j,k,ireverse,iten,nten,ntsat,bfact ", &
               i,j,k,ireverse,iten,nten,ntsat,bfact
       print *,"level,finest_level,dx ",level,finest_level, &
               dx(1),dx(2),dx(SDIM)
       print *,"xlo,xtarget,comp,ngrow ",xlo(1),xlo(2),xlo(SDIM), &
               xtarget(1),xtarget(2),xtarget(SDIM)
       print *,"TSAT_weight ",TSAT_weight
       stop
      endif

      return 
      end subroutine interpfab_tsat

       ! i,j,k is the cell containing xtarget
      subroutine interpfab_curv( &
       curv_comp, &
       nten, &
       nmat, &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       ngrow, &
       lo,hi, &
       data,DIMS(data), &
       CURV_OUT)
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: curv_comp
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xtarget(SDIM)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: DIMDEC(data)
       ! first nmat+nten components are curvatures
       ! second nmat+nten components are status (0=bad 1=good)
      REAL_T, intent(in) :: data(DIMV(data),2*(nmat+nten))
      REAL_T, intent(out) :: CURV_OUT

      INTEGER_T k1lo,k1hi
      INTEGER_T cell_index(3)
      INTEGER_T sten_index(3)
      INTEGER_T dir
      INTEGER_T nhalf
      REAL_T CURV_times_weight
      REAL_T CURV_weight
      INTEGER_T i1,j1,k1
      INTEGER_T isten,jsten,ksten
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T CURV_FLAG
      REAL_T local_CURV
      REAL_T local_weight
      REAL_T eps
      INTEGER_T nten_test

      call checkbound(lo,hi,DIMS(data),ngrow,-1,122)

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid interpfab_curv nten, nten_test ",nten,nten_test
       stop
      endif

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      if (bfact.lt.1) then 
       print *,"bfact invalid114"
       stop
      endif
      if ((curv_comp.ge.1).and.(curv_comp.le.nmat+nten)) then
       ! do nothing
      else
       print *,"curv_comp out of range"
       stop
      endif

      call containing_cell(bfact,dx,xlo,lo,xtarget,cell_index)

      do dir=1,SDIM
       if (cell_index(dir).lt.lo(dir)-ngrow) then
        cell_index(dir)=lo(dir)-ngrow
       endif
       if (cell_index(dir).gt.hi(dir)+ngrow) then
        cell_index(dir)=hi(dir)+ngrow
       endif
      enddo ! dir

      nhalf=3

      CURV_times_weight=zero
      CURV_weight=zero

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+cell_index(1)
       jsten=j1+cell_index(2)
       ksten=k1+cell_index(3)
       sten_index(1)=isten
       sten_index(2)=jsten
       sten_index(3)=ksten

       CURV_FLAG=1

       do dir=1,SDIM
        if ((sten_index(dir).lt.lo(dir)-ngrow).or. &
            (sten_index(dir).gt.hi(dir)+ngrow)) then
         CURV_FLAG=0
        endif
       enddo ! dir

       if (CURV_FLAG.eq.1) then

        call gridsten_level(xsten,isten,jsten,ksten,level,nhalf)

        CURV_FLAG=NINT(data(D_DECL(isten,jsten,ksten),nmat+nten+curv_comp))

        if (CURV_FLAG.eq.1) then
         local_CURV=data(D_DECL(isten,jsten,ksten),curv_comp)
         local_weight=zero
         eps=dx(1)*0.001
         do dir=1,SDIM
          local_weight=local_weight+(xsten(0,dir)-xtarget(dir))**2
         enddo 
         local_weight=local_weight+eps**2
         local_weight=one/local_weight
         CURV_weight=CURV_weight+local_weight
         CURV_times_weight=CURV_times_weight+local_weight*local_CURV
        else if (CURV_FLAG.eq.0) then
         ! do nothing
        else
         print *,"CURV_FLAG invalid"
         stop
        endif

       else if (CURV_FLAG.eq.0) then
        ! do nothing
       else
        print *,"CURV_FLAG invalid"
        stop
       endif

      enddo
      enddo
      enddo ! i1,j1,k1

      if (CURV_weight.gt.zero) then
       CURV_OUT=CURV_times_weight/CURV_weight
      else if (CURV_weight.eq.zero) then
       CURV_OUT=zero
      else
       print *,"CURV_weight invalid"
       stop
      endif

      return 
      end subroutine interpfab_curv




      subroutine interpfab_piecewise_constant( &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       comp, &
       ngrow, &
       lo,hi, &
       data,DIMS(data), &
       dest)
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xtarget(SDIM)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T, intent(in) :: comp,ngrow
      INTEGER_T, intent(in) :: DIMDEC(data)
      REAL_T, intent(in) :: data(DIMV(data),comp)
      REAL_T, intent(out) :: dest

      INTEGER_T dir
      INTEGER_T ic,jc,kc
      INTEGER_T cell_index(SDIM)

      call checkbound(lo,hi,DIMS(data),ngrow,-1,122)

      if (bfact.lt.1) then 
       print *,"bfact invalid115"
       stop
      endif
      if ((comp.lt.1).or.(comp.gt.1000)) then
       print *,"comp out of range"
       stop
      endif

      call containing_cell(bfact,dx,xlo,lo,xtarget,cell_index)

      do dir=1,SDIM
       if (cell_index(dir).lt.lo(dir)-ngrow) then
        cell_index(dir)=lo(dir)-ngrow
       endif
       if (cell_index(dir).gt.hi(dir)+ngrow) then
        cell_index(dir)=hi(dir)+ngrow
       endif
      enddo ! dir
      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      dest=data(D_DECL(ic,jc,kc),comp)

      return 
      end subroutine interpfab_piecewise_constant

      subroutine interpfabTEMP( &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       x, &  ! normal probe position
       xI, & ! closest point on interface from a cell center.
       Tsat, &
       im,nmat, &
       comp, &
       ngrow, &
       lo,hi, &
       tempfab, &
       DIMS(tempfab), &
       LS,DIMS(LS), &
       recon,DIMS(recon), &
       dest, &
       debugrate)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: debugrate
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: xI(SDIM)
      REAL_T, intent(in) :: Tsat
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T, intent(in) :: im,nmat,comp,ngrow
      INTEGER_T, intent(in) :: DIMDEC(tempfab)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      REAL_T, intent(in) :: tempfab(DIMV(tempfab),comp)
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(1+SDIM))
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(out) :: dest

      REAL_T :: DATA_FLOOR

      REAL_T :: T_out(1)

      INTEGER_T dir
      INTEGER_T ic,jc,kc
      INTEGER_T i1,j1,k1
      INTEGER_T isten,jsten,ksten
      INTEGER_T k1lo,k1hi
      INTEGER_T nhalf
      INTEGER_T vofcomp
      INTEGER_T cell_index(SDIM)
      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_stencil(-3:3,SDIM)
      REAL_T volcell
      REAL_T cencell(SDIM)
      REAL_T T_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T VF_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T LS_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      INTEGER_T cc_flag
      INTEGER_T tsat_flag
      INTEGER_T nsolve

      DATA_FLOOR=zero

      if (bfact.lt.1) then 
       print *,"bfact invalid116"
       stop
      endif
      if (ngrow.lt.1) then
       print *,"ngrow invalid"
       stop
      endif
      if ((comp.lt.1).or.(comp.gt.1000)) then
       print *,"comp out of range"
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid21"
       stop
      endif
      if ((Tsat.ge.zero).and.(Tsat.le.1.0D+99)) then
       ! do nothing
      else
       print *,"Tsat out of range 1"
       print *,"Tsat= ",Tsat
       stop
      endif

      call checkbound(lo,hi,DIMS(tempfab),ngrow,-1,1223)
      call checkbound(lo,hi,DIMS(LS),ngrow,-1,1224)
      call checkbound(lo,hi,DIMS(recon),ngrow,-1,1224)

      if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      call containing_cell(bfact,dx,xlo,lo,x,cell_index)

      do dir=1,SDIM
       if (cell_index(dir).lt.lo(dir)-ngrow+1) then
        cell_index(dir)=lo(dir)-ngrow+1
       endif
       if (cell_index(dir).gt.hi(dir)+ngrow-1) then
        cell_index(dir)=hi(dir)+ngrow-1
       endif
      enddo ! dir

      vofcomp=(im-1)*ngeom_recon+1

      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      nhalf=3
      call gridsten_level(xsten,ic,jc,kc,level,nhalf)

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+ic
       jsten=j1+jc
       ksten=k1+kc

       call gridsten_level(xsten_stencil,isten,jsten,ksten,level,nhalf)
       call Box_volumeFAST(bfact,dx,xsten_stencil,nhalf, &
         volcell,cencell,SDIM)

       T_sten(D_DECL(i1,j1,k1))=tempfab(D_DECL(isten,jsten,ksten),comp)
       VF_sten(D_DECL(i1,j1,k1))=recon(D_DECL(isten,jsten,ksten),vofcomp)
       LS_sten(D_DECL(i1,j1,k1))=LS(D_DECL(isten,jsten,ksten),im)
       do dir=1,SDIM
        XC_sten(D_DECL(i1,j1,k1),dir)= &
         recon(D_DECL(isten,jsten,ksten),vofcomp+dir)+cencell(dir)
       enddo

      enddo
      enddo
      enddo ! i1,j1,k1

       ! in: interpfabTEMP
      cc_flag=0  ! centroid -> target
      tsat_flag=1 ! use TSAT
      nsolve=1
       ! YANG LIUs routine
      call center_centroid_interchange( &
       DATA_FLOOR, &
       nsolve, &
       cc_flag, &
       tsat_flag, &
       bfact, &
       level, &
       finest_level, &
       dx,xlo, &
       xsten,nhalf, &
       T_sten, & 
       XC_sten, & 
       xI, &
       x, &  ! xtarget
       VF_sten, & 
       LS_sten, & 
       Tsat, &
       T_out)

      dest=T_out(1)

      return 
      end subroutine interpfabTEMP

      subroutine interpfab_filament_probe( &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       xI, &
       Tsat, &
       im_target_probe, &
       nmat, &
       comp_probe, &
       ngrow, &
       lo,hi, &
       tempfab, &
       DIMS(tempfab), &
       LS,DIMS(LS), &
       recon,DIMS(recon), &
       dest, &
       dxprobe_target, &
       VOF_pos_probe_counter, &
       use_supermesh)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xtarget(SDIM)
      REAL_T, intent(in) :: xI(SDIM)
      REAL_T, intent(in) :: Tsat
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T, intent(in) :: im_target_probe
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: comp_probe
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: DIMDEC(tempfab)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      REAL_T, intent(in) :: tempfab(DIMV(tempfab),comp_probe)
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(1+SDIM))
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(out) :: dest
      REAL_T, intent(inout) :: dxprobe_target !interpfab_filament_probe
      INTEGER_T, intent(inout) :: VOF_pos_probe_counter
      INTEGER_T, intent(in) :: use_supermesh

      INTEGER_T ic,jc,kc

      INTEGER_T dir
      INTEGER_T nhalf
      INTEGER_T vofcomp
      INTEGER_T cell_index(SDIM)
      REAL_T xsten_stencil(-3:3,SDIM)
      REAL_T volcell
      REAL_T cencell(SDIM)
      REAL_T T_sten
      REAL_T VF_sten
      REAL_T XC_sten(SDIM)
      REAL_T mag
      INTEGER_T material_found_in_cell

      if (bfact.lt.1) then 
       print *,"bfact invalid116"
       stop
      endif
      if ((VOF_pos_probe_counter.eq.0).or. &
          (VOF_pos_probe_counter.eq.1)) then
       ! do nothing
      else
       print *,"VOF_pos_probe_counter invalid"
       stop
      endif

      if (ngrow.lt.1) then
       print *,"ngrow invalid"
       stop
      endif
      if ((comp_probe.lt.1).or.(comp_probe.gt.1000)) then
       print *,"comp_probe out of range"
       stop
      endif
      if ((im_target_probe.lt.1).or. &
          (im_target_probe.gt.nmat)) then
       print *,"im_target_probe invalid21"
       stop
      endif
      if ((Tsat.ge.zero).and.(Tsat.le.1.0D+99)) then
       ! do nothing
      else
       print *,"Tsat out of range 2"
       print *,"Tsat=",Tsat
       stop
      endif

      call checkbound(lo,hi,DIMS(tempfab),ngrow,-1,1223)
      call checkbound(lo,hi,DIMS(LS),ngrow,-1,1224)
      call checkbound(lo,hi,DIMS(recon),ngrow,-1,1224)

        ! cell that contains xtarget
      call containing_cell(bfact,dx,xlo,lo,xtarget,cell_index)

      do dir=1,SDIM
       if (cell_index(dir).lt.lo(dir)-ngrow) then
        cell_index(dir)=lo(dir)-ngrow
       endif
       if (cell_index(dir).gt.hi(dir)+ngrow) then
        cell_index(dir)=hi(dir)+ngrow
       endif
      enddo ! dir

      vofcomp=(im_target_probe-1)*ngeom_recon+1

      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      nhalf=3

      call gridsten_level(xsten_stencil,ic,jc,kc,level,nhalf)
      call Box_volumeFAST(bfact,dx,xsten_stencil,nhalf, &
        volcell,cencell,SDIM)

       ! temperature at the material centroid of cell (isten,jsten,ksten)
      T_sten=tempfab(D_DECL(ic,jc,kc),comp_probe)
      VF_sten=recon(D_DECL(ic,jc,kc),vofcomp)
      mag=zero
      do dir=1,SDIM
       XC_sten(dir)= &
        recon(D_DECL(ic,jc,kc),vofcomp+dir)+cencell(dir)
       mag=mag+(XC_sten(dir)-cencell(dir))**2
      enddo
      mag=sqrt(mag)

      material_found_in_cell=0

      if (mag.ge.zero) then
       if ((VF_sten.ge.-VOFTOL).and. &
           (VF_sten.le.one+VOFTOL)) then
        if (VF_sten.ge.VOFTOL) then

         if (use_supermesh.eq.1) then 

          dxprobe_target=zero
          do dir=1,SDIM
           dxprobe_target=dxprobe_target+(XC_sten(dir)-xI(dir))**2
          enddo
          dxprobe_target=sqrt(dxprobe_target)
          if (dxprobe_target.gt.zero) then
           dest=T_sten
           material_found_in_cell=1
          else if (dxprobe_target.eq.zero) then
           ! do nothing
          else
           print *,"dxprobe_target invalid 1"
           print *,"dxprobe_target= ",dxprobe_target
           stop
          endif

         else if (use_supermesh.eq.0) then ! GFM
          ! do nothing
         else
          print *,"use_supermesh invalid"
          stop
         endif

        else if (VF_sten.le.VOFTOL) then
         ! do nothing
        else
         print *,"VF_sten invalid"
         stop
        endif
       else
        print *,"VF_sten invalid"
        stop
       endif
      else
       print *,"mag invalid MASS_TRANSFER_3D.F90 1953"
       stop
      endif

      if (material_found_in_cell.eq.1) then
       ! do nothing
      else if (material_found_in_cell.eq.0) then
       dest=Tsat
       dxprobe_target=zero
       do dir=1,SDIM
        dxprobe_target=dxprobe_target+(xtarget(dir)-xI(dir))**2
       enddo
       dxprobe_target=sqrt(dxprobe_target)
      else
       print *,"material_found_in_cell invalid"
       stop
      endif
      if (dxprobe_target.gt.zero) then
       ! do nothing
      else
       print *,"dxprobe_target invalid 2"
       print *,"dxprobe_target ",dxprobe_target
       do dir=1,SDIM
        print *,"dir,xtarget ",dir,xtarget(dir)
        print *,"dir,xI ",dir,xI(dir)
       enddo
       stop
      endif
      VOF_pos_probe_counter=VOF_pos_probe_counter+1

      return 
      end subroutine interpfab_filament_probe

      subroutine probe_interpolation( &
       PROBE_PARMS, &
       T_I,Y_I, &
       T_probe,Y_probe, &
       den_I_interp, &
       den_probe, &
       T_I_interp,Y_I_interp, &
       pres_I_interp, &
       vfrac_I, &  ! solids and fluids tessellate
       T_probe_raw, &
       dxprobe_target, &
       interp_valid_flag, &
       at_interface)
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE
 
      type(probe_parm_type), intent(in) :: PROBE_PARMS
      REAL_T, intent(in) :: T_I
      REAL_T, intent(in) :: Y_I
      REAL_T, intent(out) :: T_probe(2)
      REAL_T, intent(out) :: T_probe_raw(2)
      REAL_T, intent(out) :: Y_probe(2)
      REAL_T, intent(out) :: den_I_interp(2)
      REAL_T, intent(out) :: den_probe(2)
      REAL_T, intent(out) :: T_I_interp(2)
      REAL_T, intent(out) :: Y_I_interp(2)
      REAL_T, intent(out) :: pres_I_interp(2)
      REAL_T, intent(out) :: vfrac_I(2)
      REAL_T, intent(out) :: dxprobe_target(2)
      INTEGER_T, intent(out) :: interp_valid_flag(2)
      INTEGER_T, intent(out) :: at_interface
      INTEGER_T :: iprobe
      INTEGER_T LS_pos_probe_counter
      INTEGER_T LS_INT_VERY_CLOSE_counter
      INTEGER_T LS_INT_OWN_counter
      INTEGER_T VOF_pos_probe_counter
      INTEGER_T dummy_VOF_pos_probe_counter
      REAL_T xtarget_probe(SDIM)
      REAL_T xtarget_probe_micro(SDIM)
      INTEGER_T im_target_probe(2) ! source,dest
      INTEGER_T im_primary_probe(2)
      INTEGER_T im_secondary_probe(2)
      INTEGER_T im_target_probe_opp(2)
      INTEGER_T Ycomp_probe(2)
      INTEGER_T tcomp_probe(2)
      INTEGER_T dencomp_probe(2)
      REAL_T dist_probe_sanity
      INTEGER_T imls
      INTEGER_T dir
      INTEGER_T mtype
      INTEGER_T pcomp
      REAL_T LSPROBE(num_materials)
      REAL_T F_tess(num_materials)

      if ((Y_I.ge.zero).and.(Y_I.le.one)) then
       ! do nothing
      else
       print *,"Y_I out of range"
       print *,"Y_I= ",Y_I
       stop
      endif
      if ((T_I.ge.zero).and.(T_I.le.1.0D+99)) then
       ! do nothing
      else
       print *,"T_I out of range in probe_interpolation"
       print *,"T_I= ",T_I
       print *,"Y_I= ",Y_I
       stop
      endif
       ! 0=cannot do least squares interp or supermesh interp.
       ! 1=can do least squares interp
       ! 2=can do supermesh interp.
       ! iprobe=1 source
       ! iprobe=2 dest
      interp_valid_flag(1)=0
      interp_valid_flag(2)=0

      LS_pos_probe_counter=0
      LS_INT_VERY_CLOSE_counter=0
      LS_INT_OWN_counter=0
      VOF_pos_probe_counter=0

       ! tessellating volume fractions.
      call interpfabVFRAC_tess( &
       PROBE_PARMS%tid, &
       PROBE_PARMS%bfact, &
       PROBE_PARMS%level, &
       PROBE_PARMS%finest_level, &
       PROBE_PARMS%dx, &
       PROBE_PARMS%xlo, &
       PROBE_PARMS%xI, &
       PROBE_PARMS%nmat, &
       PROBE_PARMS%ngrow, &
       PROBE_PARMS%fablo, &
       PROBE_PARMS%fabhi, &
       PROBE_PARMS%recon, &
       DIMS(PROBE_PARMS%recon), &
       F_tess)

      do iprobe=1,2 ! iprobe=1 source    iprobe=2 dest

       if (iprobe.eq.1) then ! source
        do dir=1,SDIM
         xtarget_probe(dir)=PROBE_PARMS%xsrc(dir)
         xtarget_probe_micro(dir)=PROBE_PARMS%xsrc_micro(dir)
        enddo
        im_target_probe(iprobe)=PROBE_PARMS%im_source
        im_target_probe_opp(iprobe)=PROBE_PARMS%im_dest
        tcomp_probe(iprobe)=PROBE_PARMS%tcomp_source
        Ycomp_probe(iprobe)=PROBE_PARMS%Ycomp_source
        dencomp_probe(iprobe)=PROBE_PARMS%dencomp_source
        dxprobe_target(iprobe)=PROBE_PARMS%dxprobe_source
       else if (iprobe.eq.2) then  ! dest
        do dir=1,SDIM
         xtarget_probe(dir)=PROBE_PARMS%xdst(dir)
         xtarget_probe_micro(dir)=PROBE_PARMS%xdst_micro(dir)
        enddo
        im_target_probe(iprobe)=PROBE_PARMS%im_dest
        im_target_probe_opp(iprobe)=PROBE_PARMS%im_source
        tcomp_probe(iprobe)=PROBE_PARMS%tcomp_dest
        Ycomp_probe(iprobe)=PROBE_PARMS%Ycomp_dest
        dencomp_probe(iprobe)=PROBE_PARMS%dencomp_dest
        dxprobe_target(iprobe)=PROBE_PARMS%dxprobe_dest
       else
        print *,"iprobe invalid"
        stop
       endif
  
        ! imls_I dominates at the interface. 
       if (PROBE_PARMS%imls_I.eq.im_target_probe(iprobe)) then
        LS_INT_OWN_counter=LS_INT_OWN_counter+1
       else if ((PROBE_PARMS%imls_I.ge.1).and. &
                (PROBE_PARMS%imls_I.le.PROBE_PARMS%nmat)) then
        ! do nothing
       else
        print *,"imls_I invalid"
        stop
       endif

       if (PROBE_PARMS%LSINT(im_target_probe(iprobe)).ge. &
           -PROBE_PARMS%dxmaxLS) then
        LS_INT_VERY_CLOSE_counter=LS_INT_VERY_CLOSE_counter+1
       else if (PROBE_PARMS%LSINT(im_target_probe(iprobe)).le. &
                -PROBE_PARMS%dxmaxLS) then
        ! do nothing
       else
        print *,"LSINT(im_target_probe) invalid"
        stop
       endif

       vfrac_I(iprobe)=F_tess(im_target_probe(iprobe))

       mtype=fort_material_type(im_target_probe(iprobe))
       if ((mtype.ge.0).and. &
           (mtype.le.MAX_NUM_EOS)) then
        call interpfabFWEIGHT( &
         PROBE_PARMS%bfact, &
         PROBE_PARMS%level, &
         PROBE_PARMS%finest_level, &
         PROBE_PARMS%dx, &
         PROBE_PARMS%xlo, &
         PROBE_PARMS%xI, &
         im_target_probe(iprobe), &
         PROBE_PARMS%nmat, &
         dencomp_probe(iprobe), &
         PROBE_PARMS%ngrow, &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &        ! Fortran array box
         DIMS(PROBE_PARMS%EOS), &  ! Fortran array box
         PROBE_PARMS%recon, &
         DIMS(PROBE_PARMS%recon), &
         den_I_interp(iprobe))

        call interpfabFWEIGHT( &
         PROBE_PARMS%bfact, &
         PROBE_PARMS%level, &
         PROBE_PARMS%finest_level, &
         PROBE_PARMS%dx, &
         PROBE_PARMS%xlo, &
         xtarget_probe, &
         im_target_probe(iprobe), &
         PROBE_PARMS%nmat, &
         dencomp_probe(iprobe), &
         PROBE_PARMS%ngrow, &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &        ! Fortran array box
         DIMS(PROBE_PARMS%EOS), &  ! Fortran array box
         PROBE_PARMS%recon, &
         DIMS(PROBE_PARMS%recon), &
         den_probe(iprobe))

       else
        print *,"mtype invalid"
        stop
       endif

       ! centroid -> target (cc_flag==0)
       ! tsat_flag==1
       ! call center_centroid_interchange
       call interpfabTEMP( &
        PROBE_PARMS%bfact, &
        PROBE_PARMS%level, &
        PROBE_PARMS%finest_level, &
        PROBE_PARMS%dx, &
        PROBE_PARMS%xlo, &
        xtarget_probe, &
        PROBE_PARMS%xI, &
        T_I, &
        im_target_probe(iprobe), &
        PROBE_PARMS%nmat, &
        tcomp_probe(iprobe), &
        PROBE_PARMS%ngrow, &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%EOS, &
        DIMS(PROBE_PARMS%EOS), &
        PROBE_PARMS%LS, &
        DIMS(PROBE_PARMS%LS), &
        PROBE_PARMS%recon, &
        DIMS(PROBE_PARMS%recon), &
        T_probe(iprobe), &
        PROBE_PARMS%debugrate)

       if (T_probe(iprobe).lt.zero) then
        print *,"T_probe went negative"
        print *,"T_probe ",T_probe(iprobe)
        stop
       endif

       if (Ycomp_probe(iprobe).ge.1) then

        ! centroid -> target (cc_flag==0)
        ! tsat_flag==1
        ! call center_centroid_interchange
        call interpfabTEMP( &
         PROBE_PARMS%bfact, &
         PROBE_PARMS%level, &
         PROBE_PARMS%finest_level, &
         PROBE_PARMS%dx, &
         PROBE_PARMS%xlo, &
         xtarget_probe, &
         PROBE_PARMS%xI, &
         Y_I, &
         im_target_probe(iprobe), &
         PROBE_PARMS%nmat, &
         Ycomp_probe(iprobe), &
         PROBE_PARMS%ngrow, &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &
         DIMS(PROBE_PARMS%EOS), &
         PROBE_PARMS%LS, &
         DIMS(PROBE_PARMS%LS), &
         PROBE_PARMS%recon, &
         DIMS(PROBE_PARMS%recon), &
         Y_probe(iprobe), &
         PROBE_PARMS%debugrate)

        if ((Y_probe(iprobe).ge.-VOFTOL).and. &
            (Y_probe(iprobe).le.zero)) then
         Y_probe(iprobe)=zero
        else if ((Y_probe(iprobe).ge.zero).and. &
                 (Y_probe(iprobe).le.one)) then
         ! do nothing
        else if ((Y_probe(iprobe).ge.one).and. &
                 (Y_probe(iprobe).le.one+VOFTOL)) then
         Y_probe(iprobe)=one
        else
         print *,"Y_probe out of bounds probe_interpolation"
         print *,"Y_probe ",Y_probe(iprobe)
         stop
        endif

        call interpfabFWEIGHT( &
         PROBE_PARMS%bfact, &
         PROBE_PARMS%level, &
         PROBE_PARMS%finest_level, &
         PROBE_PARMS%dx, &
         PROBE_PARMS%xlo, &
         PROBE_PARMS%xI, &
         im_target_probe(iprobe), &
         PROBE_PARMS%nmat, &
         Ycomp_probe(iprobe), &
         PROBE_PARMS%ngrow, &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &
         DIMS(PROBE_PARMS%EOS), &
         PROBE_PARMS%recon, &
         DIMS(PROBE_PARMS%recon), &
         Y_I_interp(iprobe))
       else if (Ycomp_probe(iprobe).eq.0) then
        Y_probe(iprobe)=one
        Y_I_interp(iprobe)=one
       else
        print *,"Ycomp_probe invalid"
        stop
       endif

       do imls=1,PROBE_PARMS%nmat
         ! center -> target (cc_flag==1)
         ! tsat_flag==-1
         ! call center_centroid_interchange
        call interpfab( &
         PROBE_PARMS%bfact, &
         PROBE_PARMS%level, &
         PROBE_PARMS%finest_level, &
         PROBE_PARMS%dx, &
         PROBE_PARMS%xlo, &
         xtarget_probe, &
         imls, &
         PROBE_PARMS%ngrow, &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%LS, &
         DIMS(PROBE_PARMS%LS), &
         LSPROBE(imls))
       enddo ! imls=1..nmat

       call get_primary_material(LSPROBE, &
        PROBE_PARMS%nmat, &
        im_primary_probe(iprobe))

       if (DEBUG_TRIPLE.eq.1) then
        if ((DEBUG_I.eq.PROBE_PARMS%i).and. &
            (DEBUG_J.eq.PROBE_PARMS%j)) then
         print *,"i,j,im_primary,im_target ", &
          PROBE_PARMS%i,PROBE_PARMS%j, &
          im_primary_probe(iprobe), &
          im_target_probe(iprobe)
        endif
       endif

       if (im_primary_probe(iprobe).eq. &
           im_target_probe(iprobe)) then

        interp_valid_flag(iprobe)=1

        LS_pos_probe_counter=LS_pos_probe_counter+1

        call grad_probe_sanity( &
         PROBE_PARMS%xI, &
         xtarget_probe, &
         T_probe(iprobe), &
         T_I, &
         PROBE_PARMS%LL)

        call grad_probe_sanity( &
         PROBE_PARMS%xI, &
         xtarget_probe, &
         Y_probe(iprobe), &
         Y_I, &
         PROBE_PARMS%LL)

       else if ((im_primary_probe(iprobe).ne. &
                 im_target_probe(iprobe)).and. &
                (im_primary_probe(iprobe).ge.1).and. &
                (im_primary_probe(iprobe).le. &
                PROBE_PARMS%nmat)) then

        ! default value for T_probe
        T_probe(iprobe)=T_I
        ! default value for Y_probe
        Y_probe(iprobe)=Y_I

        call get_secondary_material(LSPROBE, &
         PROBE_PARMS%nmat, &
         im_primary_probe(iprobe), &
         im_secondary_probe(iprobe))

        dist_probe_sanity=two*dxprobe_target(iprobe)

        if ((im_secondary_probe(iprobe).eq. &
             im_target_probe(iprobe)).and. &
            (LSPROBE(im_target_probe(iprobe)).ge. &
             -dist_probe_sanity)) then

         interp_valid_flag(iprobe)=2

         dummy_VOF_pos_probe_counter=VOF_pos_probe_counter

         ! (a) find containing cell for xtarget_probe_micro
         !   xtarget_probe_micro=xI-microscale_probe_size * n * dx
         ! (b) if F(im_target_probe)<TOL in containing cell, then
         !     temp_target_probe=TSAT and dxprobe_target=
         !     ||x_probe-x_I||
         ! (c) if F(im_target_probe)>TOL in containing cell, then
         !     temp_target_probe=T(containing_cell,im_target)
         !     dxprobe_target=||x_centroid-x_I||
         call interpfab_filament_probe( &
          PROBE_PARMS%bfact, &
          PROBE_PARMS%level, &
          PROBE_PARMS%finest_level, &
          PROBE_PARMS%dx, &
          PROBE_PARMS%xlo, &
          xtarget_probe_micro, &
          PROBE_PARMS%xI, &
          T_I, &
          im_target_probe(iprobe), &
          PROBE_PARMS%nmat, &
          tcomp_probe(iprobe), &
          PROBE_PARMS%ngrow, &
          PROBE_PARMS%fablo, &
          PROBE_PARMS%fabhi, &
          PROBE_PARMS%EOS, &
          DIMS(PROBE_PARMS%EOS), &
          PROBE_PARMS%LS, &
          DIMS(PROBE_PARMS%LS), &
          PROBE_PARMS%recon, &
          DIMS(PROBE_PARMS%recon), &
          T_probe(iprobe), &
          dxprobe_target(iprobe), &
          VOF_pos_probe_counter, &
          PROBE_PARMS%use_supermesh)

         if (DEBUG_TRIPLE.eq.1) then
          if ((DEBUG_I.eq.PROBE_PARMS%i).and. &
              (DEBUG_J.eq.PROBE_PARMS%j)) then
           print *,"i,j,VOF_pos_probe_counter,iprobe ", &
            PROBE_PARMS%i,PROBE_PARMS%j,VOF_pos_probe_counter,iprobe
          endif
         endif

         if (Ycomp_probe(iprobe).ge.1) then
          ! find the mass fraction at a probe (centroid)
          ! location.
          call interpfab_filament_probe( &
           PROBE_PARMS%bfact, &
           PROBE_PARMS%level, &
           PROBE_PARMS%finest_level, &
           PROBE_PARMS%dx, &
           PROBE_PARMS%xlo, &
           xtarget_probe_micro, &
           PROBE_PARMS%xI, &
           Y_I, &
           im_target_probe(iprobe), &
           PROBE_PARMS%nmat, &
           Ycomp_probe(iprobe), &
           PROBE_PARMS%ngrow, &
           PROBE_PARMS%fablo, &
           PROBE_PARMS%fabhi, &
           PROBE_PARMS%EOS, &
           DIMS(PROBE_PARMS%EOS), &
           PROBE_PARMS%LS, &
           DIMS(PROBE_PARMS%LS), &
           PROBE_PARMS%recon, &
           DIMS(PROBE_PARMS%recon), &
           Y_probe(iprobe), &
           dxprobe_target(iprobe), &
           dummy_VOF_pos_probe_counter, &
           PROBE_PARMS%use_supermesh)


          if ((Y_probe(iprobe).ge.-VOFTOL).and. &
              (Y_probe(iprobe).le.zero)) then
           Y_probe(iprobe)=zero
          else if ((Y_probe(iprobe).ge.zero).and. &
                   (Y_probe(iprobe).le.one)) then
           ! do nothing
          else if ((Y_probe(iprobe).ge.one).and. &
                   (Y_probe(iprobe).le.one+VOFTOL)) then
           Y_probe(iprobe)=one
          else
           print *,"Y_probe out of bounds probe_interpolation (filament probe)"
           print *,"Y_probe ",Y_probe(iprobe)
           stop
          endif

         else if (Ycomp_probe(iprobe).eq.0) then
          Y_probe(iprobe)=one
         else
          print *,"Ycomp_probe invalid"
          stop
         endif

        else if ((im_secondary_probe(iprobe).ne. &
                  im_target_probe(iprobe)).or. &
                 (LSPROBE(im_target_probe(iprobe)).le. &
                  -dist_probe_sanity)) then
         T_probe(iprobe)=T_I
         Y_probe(iprobe)=Y_I
        else
         print *,"probe parameters bust"
         stop
        endif
       else 
        print *,"im_primary_probe invalid"
        stop
       endif

       ! centroid -> target (cc_flag==0)
       ! tsat_flag==0 do not use TSAT
       ! call center_centroid_interchange
       ! distance and volume fraction weighted linear least
       ! squares.  If matrix system is singular, then
       ! zeroth order least squares is used.
       call interpfabFWEIGHT( &
        PROBE_PARMS%bfact, &
        PROBE_PARMS%level, &
        PROBE_PARMS%finest_level, &
        PROBE_PARMS%dx, &
        PROBE_PARMS%xlo, &
        xtarget_probe, &
        im_target_probe(iprobe), &
        PROBE_PARMS%nmat, &
        tcomp_probe(iprobe), &
        PROBE_PARMS%ngrow, &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%EOS, &
        DIMS(PROBE_PARMS%EOS), &
        PROBE_PARMS%recon, &
        DIMS(PROBE_PARMS%recon), &
        T_probe_raw(iprobe))

       call interpfabFWEIGHT( &
        PROBE_PARMS%bfact, &
        PROBE_PARMS%level, &
        PROBE_PARMS%finest_level, &
        PROBE_PARMS%dx, &
        PROBE_PARMS%xlo, &
        PROBE_PARMS%xI, &
        im_target_probe(iprobe), &
        PROBE_PARMS%nmat, &
        tcomp_probe(iprobe), &
        PROBE_PARMS%ngrow, &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%EOS, &
        DIMS(PROBE_PARMS%EOS), &
        PROBE_PARMS%recon, &
        DIMS(PROBE_PARMS%recon), &
        T_I_interp(iprobe))

       pcomp=1
       call interpfab( &
        PROBE_PARMS%bfact, &
        PROBE_PARMS%level, &
        PROBE_PARMS%finest_level, &
        PROBE_PARMS%dx, &
        PROBE_PARMS%xlo, &
        PROBE_PARMS%xI, &
        pcomp, &
        PROBE_PARMS%ngrow, &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%pres, &
        DIMS(PROBE_PARMS%pres), &
        pres_I_interp(iprobe))

       ! local_freezing_model=0 (sharp interface stefan model)
       ! local_freezing_model=1 (source term model)
       ! local_freezing_model=2 (hydrate model)
       ! local_freezing_model=3 (wildfire)
       ! local_freezing_model=4 (source term model - Tanasawa Model
       !  or Schrage)
       ! local_freezing_model=5 (evaporation/condensation)
       ! local_freezing_model=6 (evaporation/condensation Palmore)
       if ((PROBE_PARMS%local_freezing_model.eq.0).or. & !fully saturated
           (PROBE_PARMS%local_freezing_model.eq.5).or. & !Stefan evap/cond
           (PROBE_PARMS%local_freezing_model.eq.6)) then !Palmore,Desjardins
        ! do nothing
       else if ((PROBE_PARMS%local_freezing_model.eq.1).or. & !source term
                (PROBE_PARMS%local_freezing_model.eq.2).or. & !hydrate
                (PROBE_PARMS%local_freezing_model.eq.4).or. & !Tanasawa,Schrage
                (PROBE_PARMS%local_freezing_model.eq.7)) then !Cavitation

       else if (PROBE_PARMS%local_freezing_model.eq.7) then ! cavitation
        print *,"cavitation model still under construction"
        stop
       else if (PROBE_PARMS%local_freezing_model.ne.2) then
        ! do nothing
       else
        print *,"PROBE_PARMS%local_freezing_model bust"
        stop
       endif  ! hydrate

      enddo ! iprobe=1..2

      at_interface=0

      if (DEBUG_TRIPLE.eq.1) then
       if ((DEBUG_I.eq.PROBE_PARMS%i).and. &
           (DEBUG_J.eq.PROBE_PARMS%j)) then
        print *,"i,j,LS_pos_probe_counter ", &
         PROBE_PARMS%i,PROBE_PARMS%j,LS_pos_probe_counter 
        print *,"i,j,LS_INT_VERY_CLOSE_counter ", &
         PROBE_PARMS%i,PROBE_PARMS%j,LS_INT_VERY_CLOSE_counter
        print *,"i,j,LS_INT_OWN_counter ", &
         PROBE_PARMS%i,PROBE_PARMS%j,LS_INT_OWN_counter
        print *,"i,j,VOF_pos_probe_counter ", &
         PROBE_PARMS%i,PROBE_PARMS%j,VOF_pos_probe_counter
       endif
      endif

      if ((interp_valid_flag(1).ge.1).and. &
          (interp_valid_flag(2).ge.1)) then

       ! LS_pos_probe_counter is incremented when
       ! im_primary_probe(iprobe)==im_target_probe(iprobe)
       if ((LS_pos_probe_counter.eq.1).or. &
           (LS_pos_probe_counter.eq.2)) then
        ! LS_INT_VERY_CLOSE_counter is incremented when
        ! LSINT(im_target_probe(iprobe)).ge.-dxmaxLS
        if (LS_INT_VERY_CLOSE_counter.eq.2) then
         if (LS_INT_OWN_counter.eq.1) then
          if (LS_pos_probe_counter+VOF_pos_probe_counter.eq.2) then
           at_interface=1
          else if (VOF_pos_probe_counter.eq.0) then
           ! do nothing
          else
           print *,"VOF_pos_probe_counter invalid"
           stop
          endif
         else if (LS_INT_OWN_counter.eq.0) then
          ! do nothing
         else
          print *,"LS_INT_OWN_counter invalid"
          stop
         endif 
        else if ((LS_INT_VERY_CLOSE_counter.eq.1).or. &
                 (LS_INT_VERY_CLOSE_counter.eq.0)) then
         ! do nothing
        else
         print *,"LS_INT_VERY_CLOSE_counter invalid"
         stop
        endif
       else if (LS_pos_probe_counter.eq.0) then
        ! do nothing
       else
        print *,"LS_pos_probe_counter invalid"
        stop
       endif

      else if ((interp_valid_flag(1).eq.0).or. &
               (interp_valid_flag(2).eq.0)) then
       ! do nothing
      else
       print *,"interp_valid_flag invalid"
       stop
      endif
 
      return
      end subroutine probe_interpolation

      subroutine mdot_from_Y_probe( &
       TSAT_Y_PARMS, &
       Y_gamma,T_gamma, &
       mdotY_top,mdotY_bot,mdotY)
      use global_utility_module
      IMPLICIT NONE

      type(TSAT_MASS_FRAC_parm_type), intent(in) :: TSAT_Y_PARMS
      REAL_T, intent(in) :: Y_gamma
      REAL_T, intent(in) :: T_gamma
      REAL_T, intent(out) :: mdotY_top,mdotY_bot,mdotY
      REAL_T D_MASS
      REAL_T T_probe(2)
      REAL_T Y_probe(2)
      REAL_T den_I_interp(2)
      REAL_T den_probe(2)
      REAL_T T_I_interp(2)
      REAL_T Y_I_interp(2)
      REAL_T pres_I_interp(2)
      REAL_T vfrac_I(2)
      REAL_T T_probe_raw(2)
      REAL_T dxprobe_target(2)
      INTEGER_T interp_valid_flag(2)
      INTEGER_T at_interface
      REAL_T LL
      INTEGER_T iprobe_vapor
      REAL_T den_G
      REAL_T YI_min
      INTEGER_T Kassemi_flag
      REAL_T Pgamma
      REAL_T density_probe,Tvapor_probe,internal_energy,Pvapor_probe
      INTEGER_T im_probe
      INTEGER_T imattype
      REAL_T massfrac_parm(num_species_var+1)
      INTEGER_T ispec

      YI_min=TSAT_Y_PARMS%YI_min

      if ((Y_gamma.ge.YI_min).and. &
          (YI_min.ge.zero).and. &
          (YI_min.le.one).and. &
          (Y_gamma.le.one)) then

       !iprobe=1 source
       !iprobe=2 dest
       call probe_interpolation( &
        TSAT_Y_PARMS%PROBE_PARMS, &
        T_gamma,Y_gamma, &
        T_probe,Y_probe, &
        den_I_interp, &
        den_probe, &
        T_I_interp,Y_I_interp, &
        pres_I_interp, &
        vfrac_I, &
        T_probe_raw, &
        dxprobe_target, &
        interp_valid_flag, &
        at_interface)

       if (at_interface.eq.1) then

        LL=TSAT_Y_PARMS%PROBE_PARMS%LL
        if (LL.gt.zero) then
         ! do nothing (evaporation)
        else if (LL.lt.zero) then
         ! do nothing (condensation)
        else
         print *,"LL invalid"
         stop
        endif

        iprobe_vapor=TSAT_Y_PARMS%iprobe_vapor

        if ((dxprobe_target(1).gt.zero).and. &
            (dxprobe_target(2).gt.zero)) then

         Kassemi_flag=TSAT_Y_PARMS%Tanasawa_or_Schrage_or_Kassemi

         if (Kassemi_flag.eq.0) then

          D_MASS=TSAT_Y_PARMS%D_MASS
          den_G=TSAT_Y_PARMS%den_G
          if ((D_MASS.gt.zero).and.(den_G.gt.zero)) then

           mdotY_top=Y_gamma-Y_probe(iprobe_vapor)
           mdotY_top=den_G*D_MASS*mdotY_top/ &
                  dxprobe_target(iprobe_vapor)
           mdotY_bot=one-Y_gamma
           if (mdotY_bot.lt.zero) then
            print *,"Y_gamma cannot exceed 1"
            print *,"mdotY_top= ",mdotY_top
            print *,"mdotY_bot= ",mdotY_bot
            print *,"den_G= ",den_G
            print *,"D_mass= ",D_mass
            print *,"Y_gamma= ",Y_gamma
            print *,"iprobe_vapor=",iprobe_vapor
            print *,"Y_probe(iprobe_vapor)=",Y_probe(iprobe_vapor)
            print *,"dxprobe_target(iprobe_vapor)=", &
                   dxprobe_target(iprobe_vapor)
            stop
           else if (mdotY_bot.eq.zero) then
            mdotY=zero
           else if ((mdotY_bot.gt.zero).and. &
                    (mdotY_bot.le.one)) then
            mdotY=mdotY_top/mdotY_bot
           else
            print *,"mdotY_bot invalid"
            stop
           endif
          else
           print *,"D_MASS or den_G invalid"
           stop
          endif

         else if (Kassemi_flag.eq.3) then

          call Pgamma_Clausius_Clapyron(Pgamma, &
           TSAT_Y_PARMS%reference_pressure, &
           T_gamma, &
           TSAT_Y_PARMS%TSAT_base, &
           LL, &
           TSAT_Y_PARMS%universal_gas_constant_R, &
           TSAT_Y_PARMS%molar_mass_vapor)

          density_probe=den_I_interp(iprobe_vapor)
          Tvapor_probe=T_probe(iprobe_vapor)

          if (LL.gt.zero) then 
           im_probe=TSAT_Y_PARMS%PROBE_PARMS%im_dest
          else if (LL.lt.zero) then
           im_probe=TSAT_Y_PARMS%PROBE_PARMS%im_source
          else
           print *,"LL invalid"
           stop
          endif

          call init_massfrac_parm(density_probe,massfrac_parm,im_probe)
          do ispec=1,num_species_var
           massfrac_parm(ispec)=Y_gamma
          enddo

          imattype=TSAT_Y_PARMS%material_type_evap(im_probe)
          call INTERNAL_material(density_probe,massfrac_parm, &
            Tvapor_probe,internal_energy,imattype,im_probe)
          call EOS_material(density_probe,massfrac_parm,internal_energy, &
            Pvapor_probe,imattype,im_probe)

          call MDOT_Kassemi( &
             TSAT_Y_PARMS%accommodation_coefficient, &
             TSAT_Y_PARMS%molar_mass_vapor, &
             TSAT_Y_PARMS%universal_gas_constant_R, &
             Pgamma, &
             Pvapor_probe, &
             T_gamma, &
             Tvapor_probe, &
             mdotY)
          mdotY_top=mdotY
          mdotY_bot=one

         else
          print *,"Kassemi_flag invalid"
          stop
         endif

        else
         print *,"dxprobe_target invalid"
         stop
        endif

       else
        print *,"at_interface invalid in mdot_from_Y_probe ",at_interface
        stop
       endif

      else
       print *,"Y_gamma or YI_min invalid"
       stop
      endif

      end subroutine mdot_from_Y_probe


      subroutine mdot_from_T_probe( &
       TSAT_Y_PARMS, &
       T_gamma,Y_gamma, &
       mdotT)
      IMPLICIT NONE

      type(TSAT_MASS_FRAC_parm_type), intent(in) :: TSAT_Y_PARMS
      REAL_T, intent(in) :: T_gamma
      REAL_T, intent(in) :: Y_gamma
      REAL_T, intent(out) :: mdotT
      REAL_T D_MASS
      REAL_T T_probe(2)
      REAL_T Y_probe(2)
      REAL_T den_I_interp(2)
      REAL_T den_probe(2)
      REAL_T T_I_interp(2)
      REAL_T Y_I_interp(2)
      REAL_T pres_I_interp(2)
      REAL_T vfrac_I(2)
      REAL_T T_probe_raw(2)
      REAL_T dxprobe_target(2)
      INTEGER_T interp_valid_flag(2)
      INTEGER_T at_interface
      REAL_T LL
      INTEGER_T iprobe_vapor
      REAL_T den_G
      REAL_T YI_min
      REAL_T wt(2)
      INTEGER_T iprobe

      YI_min=TSAT_Y_PARMS%YI_min

      if ((Y_gamma.ge.YI_min).and. &
          (YI_min.ge.zero).and. &
          (YI_min.le.one).and. &
          (Y_gamma.le.one)) then

       D_MASS=TSAT_Y_PARMS%D_MASS
       den_G=TSAT_Y_PARMS%den_G
       if ((D_MASS.gt.zero).and.(den_G.gt.zero)) then

        !iprobe=1 source
        !iprobe=2 dest
        call probe_interpolation( &
         TSAT_Y_PARMS%PROBE_PARMS, &
         T_gamma,Y_gamma, &
         T_probe,Y_probe, &
         den_I_interp, &
         den_probe, &
         T_I_interp,Y_I_interp, &
         pres_I_interp, &
         vfrac_I, &
         T_probe_raw, &
         dxprobe_target, &
         interp_valid_flag, &
         at_interface)

        if (at_interface.eq.1) then

         LL=TSAT_Y_PARMS%PROBE_PARMS%LL
         if (LL.gt.zero) then
          ! do nothing (evaporation)
         else if (LL.lt.zero) then
          ! do nothing (condensation)
         else
          print *,"LL invalid"
          stop
         endif
         iprobe_vapor=TSAT_Y_PARMS%iprobe_vapor

         if ((dxprobe_target(1).gt.zero).and. &
             (dxprobe_target(2).gt.zero)) then
        
          do iprobe=1,2
           wt(iprobe)=TSAT_Y_PARMS%thermal_k(iprobe)/ &
                  (abs(LL)*dxprobe_target(iprobe))
           if (wt(iprobe).ge.zero) then
            ! do nothing
           else
            print *,"wt(iprobe) invalid"
            stop
           endif
          enddo ! iprobe=1,2

          if (wt(1)+wt(2).gt.zero) then
           mdotT=wt(1)*(T_probe(1)-T_gamma)+ &
                 wt(2)*(T_probe(2)-T_gamma)
          else
           print *,"expecting wt(1)+wt(2) to be positive"
           stop
          endif
          if (LL.gt.zero) then
           ! do nothing
          else if (LL.lt.zero) then
           mdotT=-mdotT
          else
           print *,"LL invalid"
           stop
          endif 

         else
          print *,"dxprobe_target invalid"
          stop
         endif
        else
         print *,"at_interface invalid in mdot_from_T_probe ",at_interface
         stop
        endif

       else
        print *,"D_MASS or den_G invalid"
        stop
       endif
      else
       print *,"Y_gamma or YI_min invalid"
       stop
      endif

      end subroutine mdot_from_T_probe


      subroutine mdot_diff_func( &
       TSAT_Y_PARMS, &
       Y_gamma, &
       T_gamma, &
       mdot_diff)
      IMPLICIT NONE

      type(TSAT_MASS_FRAC_parm_type), intent(in) :: TSAT_Y_PARMS
      REAL_T, intent(in) :: T_gamma
      REAL_T, intent(in) :: Y_gamma
      REAL_T, intent(out) :: mdot_diff
      REAL_T mdotT
      REAL_T mdotY_top,mdotY_bot,mdotY
      INTEGER_T Kassemi_flag

      Kassemi_flag=TSAT_Y_PARMS%Tanasawa_or_Schrage_or_Kassemi
      if (Kassemi_flag.eq.0) then
       if (Y_gamma.eq.one) then
        print *,"expecting Y_gamma lt one in mdot_diff_func"
        print *,"Y_gamma= ",Y_gamma
        stop
       else if ((Y_gamma.ge.zero).and.(Y_gamma.lt.one)) then
        ! do nothing
       else
        print *,"Y_gamma invalid: ",Y_gamma
        stop
       endif
      else if (Kassemi_flag.eq.3) then
       if (Y_gamma.eq.one) then
        ! do nothing
       else
        print *,"expecting Y_gamma==1 for Kassemi model"
        stop
       endif
      else
       print *,"Kassemi_flag invalid in mdot diff fn"
       stop
      endif

      call mdot_from_T_probe( &
       TSAT_Y_PARMS, &
       T_gamma,Y_gamma, &
       mdotT)

      call mdot_from_Y_probe( &
       TSAT_Y_PARMS, &
       Y_gamma,T_gamma, &
       mdotY_top,mdotY_bot,mdotY)

      if ((mdotT.ge.zero).or. &
          (mdotT.le.zero)) then
       if ((mdotY.ge.zero).or. &
           (mdotY.le.zero)) then
        mdot_diff=mdotT-mdotY
       else
        print *,"mdotY invalid"
        stop
       endif
      else
       print *,"mdotT invalid"
       stop
      endif

      end subroutine mdot_diff_func


      subroutine TSAT_MASS_FRAC_YMIN( &
       TSAT_Y_PARMS, &
       Y_I_MIN)
      use global_utility_module
      IMPLICIT NONE

      type(TSAT_MASS_FRAC_parm_type), intent(in) :: TSAT_Y_PARMS
      REAL_T, intent(out) :: Y_I_MIN
      REAL_T :: X_I_MIN
      REAL_T :: LL,R,TSAT_base,T_I_MAX,WA,WV
      REAL_T YI_min_old

      YI_min_old=TSAT_Y_PARMS%YI_min

      WA=TSAT_Y_PARMS%molar_mass_ambient
      WV=TSAT_Y_PARMS%molar_mass_vapor
      R=TSAT_Y_PARMS%universal_gas_constant_R
      LL=TSAT_Y_PARMS%PROBE_PARMS%LL

      if ((YI_min_old.ge.zero).and.(YI_min_old.le.one)) then

       if ((LL.gt.zero).or.(LL.lt.zero)) then

        TSAT_base=TSAT_Y_PARMS%TSAT_base
        T_I_MAX=TSAT_Y_PARMS%TI_max

        if ((TSAT_base.gt.zero).and. &
            (T_I_MAX.ge.TSAT_base)) then

         if ((WA.gt.zero).and.(WV.gt.zero).and.(R.gt.zero)) then

          if (LL.gt.zero) then
           Y_I_MIN=YI_min_old
          else if (LL.lt.zero) then

             ! X=exp(-(L*WV/R)*(one/Tgamma-one/TSAT))
           call X_from_Tgamma(X_I_MIN,T_I_MAX,TSAT_base,LL,R,WV)

           if ((X_I_MIN.gt.zero).and.(X_I_MIN.le.one)) then

            call massfrac_from_volfrac(X_I_MIN,Y_I_MIN,WA,WV)

            if (Y_I_MIN.lt.YI_min_old) then
             Y_I_MIN=YI_min_old
            endif
           else
            print *,"X_I_MIN invalid"
            stop
           endif
          else
           print *,"LL invalid"
           stop
          endif
          if ((Y_I_MIN.ge.YI_min_old).and.(Y_I_MIN.le.one)) then
           ! do nothing
          else
           print *,"Y_I_MIN invalid"
           stop
          endif
         else
          print *,"WA,WV,or R invalid"
          stop
         endif
        else
         print *,"TSAT_base or T_I_MAX invalid"
         stop
        endif
       else
        print *,"LL invalid"
        stop
       endif

      else
       print *,"YI_min_old invalid"
       stop
      endif
           
      end subroutine TSAT_MASS_FRAC_YMIN



      end module mass_transfer_module

#if (STANDALONE==1)
      module mass_transfer_cpp_module
      contains
#endif

        ! this is for unsplit advection: for phase change.
      subroutine FORT_NODEDISPLACE( &
       nmat, &
       nten, &
       nburning, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       velbc, &
       dt, &
       unode,DIMS(unode), &
       ucell,DIMS(ucell), &
       xlo,dx, &
       level,finest_level)
      use probcommon_module
      use global_utility_module
      use mass_transfer_module

      IMPLICIT NONE
       
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nburning
      INTEGER_T ncomp_per_burning
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: dt

      INTEGER_T, intent(in) :: DIMDEC(unode)
      INTEGER_T, intent(in) :: DIMDEC(ucell)
     
      REAL_T, intent(out) ::  unode(DIMV(unode),2*nten*SDIM) 
      REAL_T, intent(in) ::  ucell(DIMV(ucell),nburning) 
      INTEGER_T, intent(in) :: velbc(SDIM,2,num_materials_vel*SDIM)

      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      REAL_T xstenND(-3:3,SDIM)
     
      INTEGER_T i,j,k,klosten,khisten,i1,j1,k1
      INTEGER_T dir,itencrit

      REAL_T delta,velnd,totalwt
      INTEGER_T scomp,tag_local,nhalf
      INTEGER_T ireverse
      INTEGER_T sign_reverse

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact invalid117"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nten.eq.((nmat-1)*(nmat-1)+nmat-1)/2) then
       ! do nothing
      else
       print *,"nten invalid"
       stop
      endif
      ncomp_per_burning=SDIM
      if (nburning.eq.nten*(ncomp_per_burning+1)) then
       ! do nothing
      else
       print *,"nburning invalid"
       stop
      endif

      if (level.gt.finest_level) then
       print *,"finest_level invalid nodedisplace"
       stop
      else if (level.lt.0) then
       print *,"level invalid nodedisplace"
       stop
      endif

      if (SDIM.eq.3) then
       klosten=0
       khisten=1
      else if (SDIM.eq.2) then
       klosten=0
       khisten=0
      else
       print *,"dimension crash"
       stop
      endif

      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.3) then
       ! do nothing
      else
       print *,"levelrz invalid node displace"
       stop
      endif
       

      call checkbound(fablo,fabhi,DIMS(unode),1,-1,12)
      call checkbound(fablo,fabhi,DIMS(ucell),1,-1,12)

      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif

      call growntileboxNODE(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridstenND_level(xstenND,i,j,k,level,nhalf)
       do itencrit=1,nten
       do ireverse=0,1

        if (ireverse.eq.0) then
         sign_reverse=1
        else if (ireverse.eq.1) then
         sign_reverse=-1
        else
         print *,"ireverse invalid"
         stop
        endif

        do dir=1,ncomp_per_burning
         velnd=zero
         totalwt=zero
         do i1=0,1
         do j1=0,1
         do k1=klosten,khisten
          tag_local=NINT(ucell(D_DECL(i-i1,j-j1,k-k1),itencrit))
          if (tag_local.eq.0) then
           ! do nothing
          else if ((sign_reverse*tag_local.eq.1).or. &
                   (sign_reverse*tag_local.eq.2)) then
           scomp=nten+(itencrit-1)*ncomp_per_burning+dir
           velnd=velnd+ucell(D_DECL(i-i1,j-j1,k-k1),scomp)
           totalwt=totalwt+one
          else if ((-sign_reverse*tag_local.eq.1).or. &
                   (-sign_reverse*tag_local.eq.2)) then
           ! do nothing
          else
           print *,"tag_local invalid3:",tag_local
           stop
          endif
         enddo
         enddo
         enddo
         if (totalwt.ge.one) then
          velnd=velnd/totalwt
         else if (totalwt.eq.zero) then
          ! do nothing
         else
          print *,"totalwt invalid"
          stop
         endif
         delta=dt*velnd
         if ((levelrz.eq.0).or. &
             (levelrz.eq.1)) then
          ! do nothing
         else if (levelrz.eq.3) then
          if (dir.eq.2) then
           delta=delta/xstenND(0,1)
          endif
         else
          print *,"levelrz invalid node displace 2"
          stop
         endif

         if (levelrz.eq.0) then
          ! do nothing
         else if ((levelrz.eq.1).or. &
                  (levelrz.eq.3)) then
          if (dir.eq.1) then
           if (abs(xstenND(0,dir)).le.VOFTOL*dx(dir)) then
            delta=zero
           else
            call adjust_du(delta,dir-1,xstenND(0,dir),0)  ! map_forward=0
           endif
          else if ((dir.eq.2).or.(dir.eq.SDIM)) then
           ! do nothing
          else
           print *,"dir invalid nodedisplace"
           stop
          endif 
         else
          print *,"levelrz invalid node displace 3"
          stop
         endif
         scomp=(itencrit+ireverse*nten-1)*ncomp_per_burning+dir
         unode(D_DECL(i,j,k),scomp)=delta
        enddo ! dir=1..sdim

       enddo ! ireverse=0,1
       enddo ! itencrit=1,nten

      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine FORT_NODEDISPLACE


        ! notes on phase change:
        ! 1. advection (density, temperature, species, etc)
        ! 2. (a) calculate rate of mass transfer
        !    (b) calculates divergence source term in order to
        !        preserve mass.
        ! 3. diffusion (temperature, species) 
        !    either (i) no temperature condition on interface, or
        !           (ii) some kind of Dirichlet temperature condition on
        !                interface that is changing phase.
        ! recon:
        ! vof,ref centroid,order,slope,intercept  x nmat
      subroutine FORT_CONVERTMATERIAL( &
       tid, &
       im_outer, &     ! im_outer and im_opp_outer define an interface
       im_opp_outer, & ! between im_outer and im_opp_outer.
       solvability_projection, &
       ngrow_expansion, &
       level,finest_level, &
       normal_probe_size, &
       nmat, &
       nten, &
       nden, &
       nstate, &
       ntsat, &
       supermesh_flag, &
       latent_heat, &
       saturation_temp, &
       freezing_model, &
       mass_fraction_id, &
       distribute_from_target, &
       constant_density_all_time, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       min_stefan_velocity_for_dt, &
       vofbc, &
       xlo,dx, &
       dt, &
       delta_mass, &
       maskcov,DIMS(maskcov), &
       nodevel,DIMS(nodevel), &
       JUMPFAB,DIMS(JUMPFAB), &
       TgammaFAB,DIMS(TgammaFAB), &
       LSold,DIMS(LSold), &
         ! in FORT_RATEMASSCHANGE
         ! LSnew=LSold - dt USTEFAN dot n, USTEFAN=|U|n  
         ! LSnew=LSold - dt * |U|
       LSnew,DIMS(LSnew), &
       recon,DIMS(recon), &
       snew,DIMS(snew), &
       EOS,DIMS(EOS), &
       swept,DIMS(swept) )
#if (STANDALONE==0)
      use probf90_module
#elif (STANDALONE==1)
      use probcommon_module
#endif
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
#if (STANDALONE==0)
      use hydrateReactor_module
#endif
      use mass_transfer_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: im_outer
      INTEGER_T, intent(in) :: im_opp_outer
      INTEGER_T, intent(in) :: solvability_projection
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: ngrow_expansion
      INTEGER_T, intent(in) :: normal_probe_size
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nden
      INTEGER_T, intent(in) :: nstate
      INTEGER_T, intent(in) :: ntsat
      INTEGER_T, intent(in) :: supermesh_flag
      REAL_T, intent(in) :: latent_heat(2*nten)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: mass_fraction_id(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      INTEGER_T, intent(in) :: constant_density_all_time(nmat)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in),target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: min_stefan_velocity_for_dt
      INTEGER_T, intent(in) :: vofbc(SDIM,2)
      REAL_T, intent(in),target :: xlo(SDIM)
      REAL_T, intent(in),target :: dx(SDIM)
      REAL_T, intent(in) :: dt
      REAL_T, intent(inout) :: delta_mass(2*nmat)
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(nodevel)
      INTEGER_T, intent(in) :: DIMDEC(JUMPFAB)
      INTEGER_T, intent(in) :: DIMDEC(TgammaFAB)
      INTEGER_T, intent(in) :: DIMDEC(LSold)
      INTEGER_T, intent(in) :: DIMDEC(LSnew)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(snew)
      INTEGER_T, intent(in) :: DIMDEC(EOS)
      INTEGER_T, intent(in) :: DIMDEC(swept)

      REAL_T, intent(in) :: maskcov(DIMV(maskcov))

      REAL_T, intent(in) :: nodevel(DIMV(nodevel),2*nten*SDIM)
      REAL_T, intent(out) :: JUMPFAB(DIMV(JUMPFAB),2*nten)
      REAL_T, intent(out) :: TgammaFAB(DIMV(TgammaFAB),ntsat)
      REAL_T, intent(in), target :: LSold(DIMV(LSold),nmat*(1+SDIM))
      REAL_T, intent(out) :: LSnew(DIMV(LSnew),nmat)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(out) :: snew(DIMV(snew),nstate)
      REAL_T, intent(in) :: EOS(DIMV(EOS),nden)
      REAL_T, intent(out) :: swept(DIMV(swept),nmat)

      REAL_T :: denratio_factor

      INTEGER_T i,j,k,dir
      INTEGER_T i1,j1,k1
      INTEGER_T im,im_opp
      INTEGER_T ireverse
      INTEGER_T iten
      INTEGER_T im_source
      INTEGER_T im_dest
      INTEGER_T im_dest_crit
      INTEGER_T im_source_crit
      INTEGER_T im_primary
      INTEGER_T iten_crit
      INTEGER_T iten_outer
      INTEGER_T ireverse_crit
      REAL_T max_velnode
      REAL_T velnode_test
      INTEGER_T nten_test
      INTEGER_T vcompsrc_snew,vcompdst_snew
      REAL_T oldvfrac(nmat)
      REAL_T newvfrac(nmat)
      REAL_T dF,dFdst,dFsrc
      REAL_T den_dF(2)
      REAL_T jump_strength
      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_ofs(-3:3,SDIM)

      REAL_T volgrid
      REAL_T cengrid(SDIM)
      REAL_T new_centroid(nmat,SDIM)
      REAL_T old_centroid(nmat,SDIM)
      REAL_T EBVOFTOL
      REAL_T SWEPTFACTOR
      REAL_T SWEPTFACTOR_GFM
      INTEGER_T SWEPTFACTOR_centroid
      REAL_T LL
      REAL_T Tgamma_default
      REAL_T Ygamma_default
      INTEGER_T Tsat_flag
      REAL_T energy_source
#if (STANDALONE==0)
      INTEGER_T ccomp
      REAL_T amount_used,methaneC_old,methaneC_new
#endif
      REAL_T cvtotal,wttotal,Ftemp
      INTEGER_T im_weight
      INTEGER_T tcomp_wt
      INTEGER_T vofcomp_raw
      INTEGER_T vofcomp_recon
      INTEGER_T vofcomp_raw_dest
      INTEGER_T vofcomp_recon_source
      INTEGER_T local_freezing_model
      INTEGER_T mass_frac_id
      INTEGER_T distribute_from_targ
      INTEGER_T debugrate
      REAL_T F_STEN(nmat)
      REAL_T F_STEN_CENTER(nmat)
      REAL_T unsplit_snew(nmat*ngeom_raw)
      REAL_T unsplit_density(nmat)
      REAL_T unsplit_temperature(nmat)
      REAL_T unsplit_species(nmat*(num_species_var+1))
      REAL_T unsplit_lsnew(nmat)
      REAL_T oldLS_point(nmat*(1+SDIM))
      REAL_T dxmax,dxmaxLS
      INTEGER_T recon_ncomp
      INTEGER_T scomp
 
      INTEGER_T tcomp,dencomp
      REAL_T density_data(nmat)
      REAL_T temperature_data(nmat)
      REAL_T species_data(nmat*(num_species_var+1))
      REAL_T density_mat(nmat)
      REAL_T temperature_mat(nmat)
      REAL_T species_mat(nmat*(num_species_var+1))
      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T mofdata_new(nmat*ngeom_recon)
      REAL_T multi_centroidA(nmat,SDIM)
      REAL_T volmat(nmat)
      REAL_T lsmat(nmat)
      REAL_T lsdata(nmat)
      REAL_T cenmat(SDIM,nmat)
      REAL_T multi_volume(nmat)
      REAL_T multi_area(nmat)
      REAL_T multi_cen(SDIM,nmat)

      REAL_T thermal_k(2)  ! source,dest

      INTEGER_T nmax,nhalf,nhalf0
      INTEGER_T klosten,khisten
      INTEGER_T symmetry_flag,ntetbox
      REAL_T u_xsten_grid(-3:3,SDIM)
      REAL_T dxgrid(SDIM)
      REAL_T u_xsten_updatecell(-3:3,SDIM)
      REAL_T u_xsten_departmap(-1:1,SDIM)
      REAL_T absolute_voltotal
      REAL_T voltotal
      REAL_T multi_volume_total
      INTEGER_T igrid,jgrid,kgrid
      INTEGER_T imac,jmac,kmac
      INTEGER_T inode1,jnode1,knode1,u_imaterial,u_im,itri
      INTEGER_T udir,udir2,inode,id,nlist
      INTEGER_T matrix_status
      REAL_T velnode
      REAL_T xdepartnode(4*(SDIM-1),SDIM)
      REAL_T xtargetnode(4*(SDIM-1),SDIM)
      REAL_T tempdatanode(4*(SDIM-1))
      REAL_T xinttri(SDIM+1,SDIM)
      REAL_T xtri(SDIM+1,SDIM)
      REAL_T xmaptri(SDIM+1,SDIM)
      REAL_T xtargettri(SDIM+1,SDIM)
      REAL_T datatri(SDIM+1)
      REAL_T AA(SDIM+1,SDIM+1)
      REAL_T AS(SDIM,SDIM)
      REAL_T ASINV(SDIM,SDIM)
      REAL_T bb(SDIM+1)
      REAL_T AAINV(SDIM+1,SDIM+1)
      REAL_T bbINV(SDIM+1)
      REAL_T xx(SDIM+1)
      REAL_T coeff(SDIM,SDIM+1)
      REAL_T coeffINV(SDIM,SDIM+1)
      INTEGER_T n,ivert
      REAL_T nn(SDIM)
      REAL_T nnmap(SDIM)
      REAL_T uncaptured_volume
      REAL_T uncaptured_centroid(SDIM)
      INTEGER_T shapeflag
      INTEGER_T tessellate
      REAL_T volcell
      REAL_T volcell_ofs
      REAL_T tempvfrac
      REAL_T temp_new_vfrac
      REAL_T cencell(SDIM)
      REAL_T cencell_ofs(SDIM)
      REAL_T tempcen(SDIM)
      INTEGER_T do_unsplit_advection
      INTEGER_T interface_near(2*nten)
      INTEGER_T local_mask
      INTEGER_T im_primary_new
      INTEGER_T im_primary_old
      INTEGER_T away_from_interface
      REAL_T solid_vof_new,solid_vof_old
      INTEGER_T mtype
      REAL_T vapor_den
      REAL_T condensed_den
      REAL_T local_cv_or_cp
      INTEGER_T speccomp_mod
      INTEGER_T default_comp
      INTEGER_T ncomp_per_tsat
      INTEGER_T ispec
      INTEGER_T iprobe
      INTEGER_T iprobe_vapor
      INTEGER_T iprobe_condensed
      INTEGER_T im_probe
      INTEGER_T im_vapor
      INTEGER_T im_condensed
      INTEGER_T base_index
      INTEGER_T dencomp_probe
      INTEGER_T tcomp_probe
      INTEGER_T mfrac_comp_probe
      INTEGER_T ispec_probe
      REAL_T temp_mix_new(2)
      REAL_T mass_frac_new(2)
      REAL_T density_old(2)
      REAL_T temperature_new(2)
      REAL_T temperature_old(2)
      REAL_T species_new(2)
      REAL_T species_old(2)
      REAL_T delta_mass_local(2) ! iprobe==1: source   iprobe==2: dest
      REAL_T xPOINT_supermesh(SDIM)
      REAL_T xPOINT_GFM(SDIM)
      REAL_T, target :: xstar(SDIM)
      INTEGER_T im_old_crit
      INTEGER_T im_new_crit
      REAL_T DATA_FLOOR
      INTEGER_T combine_flag
      INTEGER_T nsolve_interp
      REAL_T xtarget_interp(SDIM)
      REAL_T old_xI(SDIM)
      REAL_T old_nrm(SDIM)
      INTEGER_T interp_to_new_supermesh

      REAL_T XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T VF_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T LS_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T temperature_sten(D_DECL(-1:1,-1:1,-1:1))
      REAL_T massfrac_sten(D_DECL(-1:1,-1:1,-1:1))

      INTEGER_T continuous_mof_parm
      INTEGER_T cmofsten(D_DECL(-1:1,-1:1,-1:1))
      INTEGER_T use_ls_data
      INTEGER_T mof_verbose
      REAL_T LS_stencil(D_DECL(-1:1,-1:1,-1:1),nmat)
      type(interp_from_grid_parm_type) :: data_in 
      type(interp_from_grid_out_parm_type) :: data_out
      REAL_T, target :: cell_data_interp(nmat*(1+SDIM)) !LS1..LSn,slopes ...
      INTEGER_T order_probe(2)
      REAL_T nslope_probe(SDIM,2)
      REAL_T intercept_probe(2)
      REAL_T nslope_dest(SDIM)
      REAL_T intercept_dest
      REAL_T LS_dest_old,LS_dest_new
      REAL_T mass_frac_limit
      INTEGER_T nhalf_box
      INTEGER_T im_local
      INTEGER_T im_trust
      INTEGER_T im_distrust
      REAL_T LS_MAX_fixed
      REAL_T fixed_vfrac_sum
      REAL_T fixed_centroid_sum(SDIM)
      REAL_T avail_vfrac


      if ((tid.lt.0).or. &
          (tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: CONVERTMATERIAL
      if ((nmax.lt.100).or.(nmax.gt.2000)) then
       print *,"nmax invalid"
       stop
      endif

      nhalf_box=1
      nhalf=3
      nhalf0=1

      recon_ncomp=nmat*ngeom_recon

      ncomp_per_tsat=2

      if ((solvability_projection.ne.0).and. &
          (solvability_projection.ne.1)) then
       print *,"solvability_projection invalid"
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid"
       stop
      endif

      debugrate=0
      if (1.eq.0) then
       debugrate=1
      endif

      if (ngrow_expansion.ne.2) then
       print *,"ngrow_expansion invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in convert_material"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (normal_probe_size.ne.1) then
       print *,"normal_probe_size invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if (min_stefan_velocity_for_dt.ge.zero) then
       ! do nothing
      else
       print *,"min_stefan_velocity_for_dt invalid"
       stop
      endif

      if (im_outer.ne.im_opp_outer) then
       if ((im_outer.ge.1).and. &
           (im_outer.le.nmat)) then
        ! do nothing
       else
        print *,"im_outer invalid"
        stop
       endif
       if ((im_opp_outer.ge.1).and. &
           (im_opp_outer.le.nmat)) then
        ! do nothing
       else
        print *,"im_opp_outer invalid"
        stop
       endif
      else
       print *,"cannot have im_outer==im_opp_outer"
       stop
      endif
      call get_iten(im_outer,im_opp_outer,iten_outer,nmat)

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)

      EBVOFTOL=VOFTOL
      if (dt.lt.one) then
       EBVOFTOL=EBVOFTOL*dt
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid ratemass nten, nten_test ",nten,nten_test
       stop
      endif
       ! For Tsatfab: 1st nten components are the status, then next 2*nten
       ! components go:
       ! T_gamma_1,Y_gamma_1,
       ! T_gamma_2,Y_gamma_2, ....
      if (ntsat.eq.nten*(1+ncomp_per_tsat)) then
       ! do nothing
      else
       print *,"ntsat invalid"
       stop
      endif
      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid in convertmaterial"
       print *,"nden=",nden
       print *,"nmat=",nmat
       print *,"num_state_material=",num_state_material
       stop
      endif
      if (nstate.ne.num_materials_vel*(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif

      do im=1,nmat-1
       do im_opp=im+1,nmat
        do ireverse=0,1
         call get_iten(im,im_opp,iten,nmat)
         local_freezing_model=freezing_model(iten+ireverse*nten)
         distribute_from_targ=distribute_from_target(iten+ireverse*nten)
         LL=latent_heat(iten+ireverse*nten)
         if ((local_freezing_model.eq.0).or. & ! Stefan model
             (local_freezing_model.eq.1).or. & ! source term
             (local_freezing_model.eq.2)) then ! hydrate
          ! do nothing
         else if (local_freezing_model.eq.3) then ! wildfire
          print *,"fix me"
          stop
         else if ((local_freezing_model.eq.4).or. & ! Tanasawa or Schrage
                  (local_freezing_model.eq.5).or. & ! Stefan model evap/cond.
                  (local_freezing_model.eq.6).or. & ! Palmore/Desjardins
                  (local_freezing_model.eq.7)) then ! Cavitation
          mass_frac_id=mass_fraction_id(iten+ireverse*nten)
          if ((mass_frac_id.ge.1).and. &
              (mass_frac_id.le.num_species_var)) then
           ! do nothing
          else
           print *,"mass_frac_id invalid"
           stop
          endif
          if (local_freezing_model.eq.4) then 
           if (LL.ne.zero) then
            if (ireverse.eq.0) then ! evaporation
             if (LL.le.zero) then
              print *,"LL invalid"
              stop
             endif
             if (distribute_from_targ.ne.0) then
              print *,"distribute_from_targ invalid"
              stop
             endif
             im_source=im
             im_dest=im_opp
            else if (ireverse.eq.1) then ! condensation
             if (LL.ge.zero) then
              print *,"LL invalid"
              stop
             endif
             if (distribute_from_targ.ne.1) then
              print *,"distribute_from_targ invalid"
              stop
             endif
             im_source=im_opp
             im_dest=im
            else
             print *,"ireverse invalid"
             stop
            endif
           else if (LL.eq.zero) then
            ! do nothing
           else
            print *,"LL invalid"
            stop
           endif 
          else if (local_freezing_model.eq.5) then  ! stefan evap/cond
           ! do nothing
          else if (local_freezing_model.eq.6) then  ! Palmore/Desjardin
           ! do nothing
          else if (local_freezing_model.eq.7) then  ! Cavitation
           ! do nothing
          else
           print *,"local_freezing_model invalid in convertmaterial"
           stop
          endif

         else
          print *,"local_freezing_model invalid in convertmaterial"
          print *,"local_freezing_model= ",local_freezing_model
          print *,"iten,ireverse,nten ",iten,ireverse,nten
          stop
         endif
        enddo ! ireverse
       enddo ! im_opp
      enddo !im

      call checkbound(fablo,fabhi, &
       DIMS(maskcov),1,-1,1256)

      call checkbound(fablo,fabhi, &
       DIMS(JUMPFAB),ngrow_expansion,-1,1256)
      call checkbound(fablo,fabhi, &
       DIMS(TgammaFAB),ngrow_expansion,-1,1256)
      call checkbound(fablo,fabhi, &
       DIMS(LSold),normal_probe_size+3,-1,1257)
      call checkbound(fablo,fabhi, &
       DIMS(LSnew),1,-1,1258)
      call checkbound(fablo,fabhi, &
       DIMS(recon),1,-1,1259)
      call checkbound(fablo,fabhi, &
       DIMS(snew),1,-1,1261)
      call checkbound(fablo,fabhi, &
       DIMS(EOS),1,-1,1262)
      call checkbound(fablo,fabhi, &
       DIMS(swept),0,-1,1263)
      call checkbound(fablo,fabhi,DIMS(nodevel),1,-1,12)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      if (SDIM.eq.3) then
       klosten=-1
       khisten=1
      else if (SDIM.eq.2) then
       klosten=0
       khisten=0
      else
       print *,"dimension bust"
       stop
      endif

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        do ireverse=0,1 
         if (JUMPFAB(D_DECL(i,j,k),iten_outer+ireverse*nten).eq.zero) then
          ! do nothing
         else
          print *,"JUMPFAB(D_DECL(i,j,k),iten_outer+ireverse*nten) init bad"
          stop
         endif
        enddo

        call gridsten_level(xsten,i,j,k,level,nhalf)
        call gridsten_level(u_xsten_updatecell,i,j,k,level,nhalf)

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volgrid,cengrid,SDIM)

        do_unsplit_advection=0
        do im=1,2*nten
         interface_near(im)=0
        enddo

        do im=1,nmat
         lsmat(im)=LSold(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(lsmat,nmat,im_primary)

        if (is_rigid(nmat,im_primary).eq.0) then

         do im=1,nmat
          vofcomp_recon=(im-1)*ngeom_recon+1
          F_STEN(im)=zero
          F_STEN_CENTER(im)=recon(D_DECL(i,j,k),vofcomp_recon)

          do i1=-1,1
          do j1=-1,1
          do k1=klosten,khisten
           F_STEN(im)=F_STEN(im)+recon(D_DECL(i+i1,j+j1,k+k1),vofcomp_recon)
          enddo
          enddo
          enddo ! i1,j1,k1
         enddo ! im=1..nmat

         do im=1,nmat-1
          do im_opp=im+1,nmat

           call get_iten(im,im_opp,iten,nmat)

           do ireverse=0,1
            if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
             print *,"im or im_opp bust 10"
             stop
            endif

            LL=latent_heat(iten+ireverse*nten)

            if ((is_rigid(nmat,im).eq.1).or. &
                (is_rigid(nmat,im_opp).eq.1)) then
             ! do nothing
            else if (LL.ne.zero) then
             if (ireverse.eq.0) then
              im_source=im
              im_dest=im_opp
             else if (ireverse.eq.1) then
              im_source=im_opp
              im_dest=im
             else
              print *,"ireverse invalid"
              stop
             endif

             if (F_STEN_CENTER(im_source).gt.zero) then

              if ((F_STEN(im_source).gt.zero).and. &
                  (F_STEN(im_dest).gt.zero)) then
               if ((im.eq.im_outer).and.(im_opp.eq.im_opp_outer)) then
                interface_near(iten+ireverse*nten)=1
                do_unsplit_advection=1
               else 
                ! do nothing
               endif
              else if ((F_STEN(im_source).eq.zero).or. &
                       (F_STEN(im_dest).eq.zero)) then
               ! do nothing
              else
               print *,"F_STEN invalid"
               stop
              endif

             else if (F_STEN_CENTER(im_source).eq.zero) then
              ! do nothing
             else
              print *,"F_STEN_CENTER(im_source) invalid"
              stop
             endif

            else if (LL.eq.zero) then
             ! do nothing
            else
             print *,"LL invalid"
             stop
            endif
           enddo ! ireverse
          enddo ! im_opp
         enddo ! im

         if (DEBUG_ACTIVE_CELL.eq.1) then
          if ((i.eq.DEBUG_I).and.(j.eq.DEBUG_J)) then
           print *,"i,j,do_unsplit_advection ",i,j,do_unsplit_advection
          endif
         endif

          ! do the 3x3x3 stencil volume fractions satisfy:
          !  0 < F_source < 1 and
          !  0 < F_dest < 1
         if (do_unsplit_advection.eq.1) then

          ! determine the interface that is changing phase:
          ! find the interface with the largest velocity magnitude.
          ! The velocity, "nodevel", lives at the nodes.
          max_velnode=zero
          im_dest_crit=-1
          im_source_crit=-1
          iten_crit=-1
          ireverse_crit=-1

          do knode1=klosten,khisten,2
          do jnode1=-1,1,2
          do inode1=-1,1,2
           if (inode1.eq.-1) then
            imac=i
           else if (inode1.eq.1) then
            imac=i+1
           else
            print *,"inode1 invalid"
            stop
           endif
           if (jnode1.eq.-1) then
            jmac=j
           else if (jnode1.eq.1) then
            jmac=j+1
           else
            print *,"jnode1 invalid"
            stop
           endif
           if (SDIM.eq.2) then
            kmac=0
           else if (SDIM.eq.3) then
            if (knode1.eq.-1) then
             kmac=k
            else if (knode1.eq.1) then
             kmac=k+1
            else
             print *,"knode1 invalid"
             stop
            endif
           else
            print *,"dimension bust"
            stop
           endif

           do im=1,nmat-1
            do im_opp=im+1,nmat
             call get_iten(im,im_opp,iten,nmat)
             do ireverse=0,1
              if (interface_near(iten+ireverse*nten).eq.1) then
               if (ireverse.eq.0) then
                im_source=im
                im_dest=im_opp
               else if (ireverse.eq.1) then
                im_source=im_opp
                im_dest=im
               else
                print *,"ireverse invalid"
                stop
               endif
               if ((iten.ge.1).and.(iten.le.nten)) then
                scomp=(iten+nten*ireverse-1)*SDIM
                velnode_test=zero
                do udir=1,SDIM
                 velnode_test=velnode_test+ &
                  nodevel(D_DECL(imac,jmac,kmac),scomp+udir)**2
                enddo
               else
                print *,"iten invalid"
                stop
               endif
               if (im_dest_crit.eq.-1) then
                im_dest_crit=im_dest
                im_source_crit=im_source
                iten_crit=iten
                ireverse_crit=ireverse
                max_velnode=velnode_test
               else if ((im_dest_crit.ge.1).and.(im_dest_crit.le.nmat)) then
                if (velnode_test.gt.max_velnode) then
                 max_velnode=velnode_test
                 im_dest_crit=im_dest
                 im_source_crit=im_source
                 iten_crit=iten
                 ireverse_crit=ireverse
                endif
               else
                print *,"im_dest_crit invalid"
                stop
               endif
              else if (interface_near(iten+ireverse*nten).eq.0) then
               ! do nothing
              else
               print *,"interface_near invalid"
               stop
              endif
             enddo ! ireverse
            enddo ! im_opp
           enddo ! im
                  
          enddo
          enddo
          enddo ! inode1,jnode1,knode1

          if (im_dest_crit.eq.im_source_crit) then
           print *,"im_dest_crit or im_source_crit invalid"
           stop
          endif
          if ((im_dest_crit.lt.1).or.(im_dest_crit.gt.nmat)) then
           print *,"im_dest_crit invalid"
           stop
          endif
          if ((im_source_crit.lt.1).or.(im_source_crit.gt.nmat)) then
           print *,"im_source_crit invalid"
           stop
          endif
          if ((iten_crit.lt.1).or.(iten_crit.gt.nten)) then
           print *,"iten_crit invalid"
           stop
          endif
          if ((ireverse_crit.ne.0).and.(ireverse_crit.ne.1)) then
           print *,"ireverse_crit invalid"
           stop
          endif

           ! advection: unsplit_lsnew, unsplit_snew, dF
           !               unsplit_density,unsplit_temperature
           ! unsplit advection: backward tracing of characteristics, but
           !  evaluate volumes, centroids, level set functions, and
           !  temperatures in the target cell.

          symmetry_flag=0 
          call get_ntetbox(ntetbox,symmetry_flag,SDIM)

          absolute_voltotal=zero
          voltotal=zero
          do u_im=1,nmat
           density_mat(u_im)=zero
           temperature_mat(u_im)=zero
           do ispec=1,num_species_var
            species_mat((ispec-1)*nmat+u_im)=zero
           enddo
           volmat(u_im)=zero
           do udir=1,SDIM
            cenmat(udir,u_im)=zero
           enddo
          enddo
          do u_im=1,nmat
           lsmat(u_im)=zero
          enddo

           ! backwards tracing of characteristics:
           ! add up contributions from neighboring cells.
           ! (i,j,k) is the cell to be updated.
          do igrid=-1,1
          do jgrid=-1,1
          do kgrid=klosten,khisten

           call gridsten_level(u_xsten_grid, &
             i+igrid,j+jgrid,k+kgrid,level,nhalf)

           inode=1

            ! index order must be knode1,jnode1,inode1
           do knode1=klosten,khisten,2
           do jnode1=-1,1,2
           do inode1=-1,1,2
            if (inode1.eq.-1) then
             imac=i
            else if (inode1.eq.1) then
             imac=i+1
            else
             print *,"inode1 invalid"
             stop
            endif
            if (jnode1.eq.-1) then
             jmac=j
            else if (jnode1.eq.1) then
             jmac=j+1
            else
             print *,"jnode1 invalid"
             stop
            endif
            if (SDIM.eq.2) then
             kmac=0
            else if (SDIM.eq.3) then
             if (knode1.eq.-1) then
              kmac=k
             else if (knode1.eq.1) then
              kmac=k+1
             else
              print *,"knode1 invalid"
              stop
             endif
            else
             print *,"dimension bust"
             stop
            endif

            scomp=(iten_crit+ireverse_crit*nten-1)*SDIM

             ! node velocity initialized in "nodedisplace"
            do udir=1,SDIM
             velnode=nodevel(D_DECL(imac,jmac,kmac),scomp+udir)

             if (DEBUG_ACTIVE_CELL.eq.1) then
              if ((i.eq.DEBUG_I).and.(j.eq.DEBUG_J)) then
               print *,"iten_crit,ireverse_crit,udir ", &
                       iten_crit,ireverse_crit,udir
               print *,"inode1,jnode1,knode1,velnode ", &
                       inode1,jnode1,knode1,velnode
              endif
             endif

             if (udir.eq.1) then
              xtargetnode(inode,udir)=u_xsten_updatecell(inode1,udir)
             else if (udir.eq.2) then
              xtargetnode(inode,udir)=u_xsten_updatecell(jnode1,udir)
             else if ((udir.eq.3).and.(SDIM.eq.3)) then
              xtargetnode(inode,udir)=u_xsten_updatecell(knode1,udir)
             else
              print *,"udir invalid"
              stop
             endif
             xdepartnode(inode,udir)=xtargetnode(inode,udir)-velnode
            enddo ! udir

            tempdatanode(inode)=one

            inode=inode+1
           enddo
           enddo
           enddo ! inode1,jnode1,knode1

           do id=1,ntetbox
            call extract_tet(xdepartnode,tempdatanode,xtri,datatri, &
              id,symmetry_flag,SDIM)
            call extract_tet(xtargetnode,tempdatanode, &
              xtargettri,datatri,id,symmetry_flag,SDIM)

            call intersect_cube( &
              xtri, &
              u_xsten_grid,nhalf, &
              geom_xtetlist_old(1,1,1,tid+1), &
              geom_xtetlist(1,1,1,tid+1), &
              nmax, &
              nlist,nmax,SDIM)

            if (nlist.gt.0) then

              ! find mapping from departure to target
              ! xt(dir)=coeff(dir,1)*xd(1)+coeff(dir,2)*xd(2)+coeff(dir,3)
              ! and from target back to departure
              ! xd(dir)=coeffINV(dir,1)*xt(1)+coeffINV(dir,2)*xt(2)+
              !         coeffINV(dir,3)
              ! "AS" "short A" -- coefficients of x,y's
              ! phi=n dot (x-x0)+b
              ! phi^{map}=nmap dot (x-x0map)+b
              ! phi^{map}(x)=phi(xmap^{-1}(x))
              ! nmap dot (x-x0map)+b=n dot (AS^{-1}(x-x0map)) +b
              ! nmap^{T}=n^{T}AS^{-1}  or nmap=AS^{-T}n
              ! x0map=A(x0)+gamma

             do udir=1,SDIM
              do itri=1,SDIM+1
               do udir2=1,SDIM
                AA(itri,udir2)=xtri(itri,udir2)           ! xdepart
                AAINV(itri,udir2)=xtargettri(itri,udir2)  ! xtarget
               enddo
               AA(itri,SDIM+1)=one
               AAINV(itri,SDIM+1)=one
               bb(itri)=xtargettri(itri,udir)  ! xtarget
               bbINV(itri)=xtri(itri,udir)     ! xdepart
              enddo  ! itri
              call matrix_solve(AA,xx,bb,matrix_status,SDIM+1)
              if (matrix_status.ne.1) then
               print *,"mapping bust: matrix_solve failed"
               stop
              endif
              do udir2=1,SDIM+1
               coeff(udir,udir2)=xx(udir2)
              enddo
              do udir2=1,SDIM
               AS(udir,udir2)=xx(udir2)
              enddo
              call matrix_solve(AAINV,xx,bbINV,matrix_status,SDIM+1)
              if (matrix_status.eq.0) then
               print *,"mapping bust"
               stop
              endif
              do udir2=1,SDIM+1
               coeffINV(udir,udir2)=xx(udir2)
              enddo
              do udir2=1,SDIM
               ASINV(udir,udir2)=xx(udir2)
              enddo
             enddo  ! udir

              ! find expression for interface reconstruction in the mapped
              ! target cell.  u_xsten_departmap is the mapped location of 
              ! u_xsten_grid.  
              ! Only the u_xsten_departmap(0,udir) information is used:
              ! used for evaluating n dot (x-x0).
              ! interface reconstruction in target cell is
              ! n_mapped dot (x-x0_mapped)  (*)
              ! if (*) mapped back to departure cell, one should
              ! have n dot (x-x0)
              ! given 3 points on the plane in the departure region:
              ! n dot (xi-x0)=0  i=1..3
              ! this implies that nmapped dot (xi_mapped - x0mapped)=0
              ! nmapped dot (A(xi-x0))=0  
              ! nmapped^T A(xi-x0)=0
              ! (A^T nmapped)^T (xi-x0)=0
              ! n=A^T nmapped
              ! nmapped=(A^T)^{-1} n
             do udir=1,SDIM
              dxgrid(udir)=u_xsten_grid(1,udir)-u_xsten_grid(-1,udir)
              if (dxgrid(udir).le.zero) then
               print *,"u_xsten_grid became corrupt" 
               stop
              endif 
             enddo
          
             do udir=1,SDIM
              u_xsten_departmap(0,udir)=coeff(udir,SDIM+1)
              do udir2=1,SDIM
               u_xsten_departmap(0,udir)=u_xsten_departmap(0,udir)+ &
                 coeff(udir,udir2)*u_xsten_grid(0,udir2)
              enddo ! udir2
              u_xsten_departmap(-1,udir)=u_xsten_departmap(0,udir)- &
                half*dxgrid(udir)
              u_xsten_departmap(1,udir)=u_xsten_departmap(0,udir)+ &
                half*dxgrid(udir)
             enddo ! udir

             do u_im=1,nmat
              lsdata(u_im)= &
               LSold(D_DECL(i+igrid,j+jgrid,k+kgrid),u_im) 
             enddo

             do u_im=1,nmat

              dencomp=(u_im-1)*num_state_material+1
              tcomp=dencomp+1
               ! mixture density of gas if u_im corresponds to gas
              density_data(u_im)=EOS(D_DECL(i+igrid,j+jgrid,k+kgrid),dencomp)
              temperature_data(u_im)= &
                EOS(D_DECL(i+igrid,j+jgrid,k+kgrid),tcomp)
              do ispec=1,num_species_var
               species_data((ispec-1)*nmat+u_im)= &
                EOS(D_DECL(i+igrid,j+jgrid,k+kgrid),tcomp+ispec)
              enddo
         
              vofcomp_recon=(u_im-1)*ngeom_recon+1
              mofdata(vofcomp_recon)= &
               recon(D_DECL(i+igrid,j+jgrid,k+kgrid),vofcomp_recon) 
              do udir=1,SDIM 
               mofdata(vofcomp_recon+udir)= &
                recon(D_DECL(i+igrid,j+jgrid,k+kgrid),vofcomp_recon+udir)
              enddo
              mofdata(vofcomp_recon+SDIM+1)= &
               recon(D_DECL(i+igrid,j+jgrid,k+kgrid),vofcomp_recon+SDIM+1) !ord
              mofdata(vofcomp_recon+2*SDIM+2)= &
               recon(D_DECL(i+igrid,j+jgrid,k+kgrid),vofcomp_recon+2*SDIM+2) 
              do udir=1,SDIM
               nn(udir)= &
                recon(D_DECL(i+igrid,j+jgrid,k+kgrid), &
                      vofcomp_recon+SDIM+1+udir) 
              enddo ! udir

              ! ASINV maps target back to departure
              ! e.g. in 1D, if compression => AS<1 => ASINV>1
              ! 1D:AS=(xtarget_left-xtarget_right)/(xdeptleft-xeptright)
              ! ASINV=(xdeptleft-xdeptright)/(xtarget_left-xtarget_right)
              do udir=1,SDIM
               nnmap(udir)=zero
               do udir2=1,SDIM
                nnmap(udir)=nnmap(udir)+ASINV(udir2,udir)*nn(udir2)
               enddo
               mofdata(vofcomp_recon+SDIM+1+udir)=nnmap(udir)
              enddo ! udir
             enddo ! u_im

             do n=1,nlist

               ! xinttri is in the departure region
              do ivert=1,SDIM+1
              do udir=1,SDIM
               xinttri(ivert,udir)=geom_xtetlist(ivert,udir,n,tid+1)
              enddo
              enddo
           
              ! target triangle (xinttri mapped), confined to target cell
              do ivert=1,SDIM+1
               do udir=1,SDIM
                xmaptri(ivert,udir)=coeff(udir,SDIM+1)
                do udir2=1,SDIM
                 xmaptri(ivert,udir)=xmaptri(ivert,udir)+ &
                  coeff(udir,udir2)*xinttri(ivert,udir2)
                enddo
               enddo ! udir
              enddo ! ivert

              call tetrahedron_volume(xmaptri,uncaptured_volume, &
               uncaptured_centroid,SDIM)

              absolute_voltotal=absolute_voltotal+uncaptured_volume

                ! find volumes within u_xsten_updatecell (target)
              shapeflag=1  
              tessellate=0
              call multi_get_volume_grid( &
                tessellate, & ! =0
                bfact,dx, &
                u_xsten_departmap,nhalf0, & ! nhalf0=1
                mofdata, &
                u_xsten_updatecell,nhalf, & ! nhalf=3
                xmaptri, &
                multi_volume,multi_cen,multi_area, &
                geom_xtetlist_uncapt(1,1,1,tid+1),  &
                nmax, &
                nmax, &
                nmat,SDIM,shapeflag,102)

              multi_volume_total=zero
              do u_im=1,nmat
               if (is_rigid(nmat,u_im).eq.0) then 
                multi_volume_total=multi_volume_total+multi_volume(u_im)
               else if (is_rigid(nmat,u_im).eq.1) then 
                ! do nothing
               else
                print *,"is_rigid invalid"
                stop
               endif
               density_mat(u_im)=density_mat(u_im)+ &
                 multi_volume(u_im)*density_data(u_im)
               temperature_mat(u_im)=temperature_mat(u_im)+ &
                 multi_volume(u_im)*temperature_data(u_im)
               do ispec=1,num_species_var
                species_mat((ispec-1)*nmat+u_im)= &
                 species_mat((ispec-1)*nmat+u_im)+ &
                 multi_volume(u_im)*species_data((ispec-1)*nmat+u_im)
               enddo
              enddo ! u_im=1..nmat

              voltotal=voltotal+multi_volume_total

              do u_im=1,nmat
               lsmat(u_im)=lsmat(u_im)+multi_volume_total*lsdata(u_im)
              enddo
              do u_im=1,nmat 
               volmat(u_im)=volmat(u_im)+multi_volume(u_im)
               do udir=1,SDIM
                cenmat(udir,u_im)=cenmat(udir,u_im)+ &
                  multi_volume(u_im)*multi_cen(udir,u_im)
               enddo
              enddo  ! u_im
             enddo ! n - traversing triangles in intersection 
            else if (nlist.eq.0) then
             ! do nothing
            else
             print *,"nlist invalid"
             stop
            endif 
           enddo ! id=1..ntetbox

          enddo
          enddo
          enddo ! igrid,jgrid,kgrid
          
          if (voltotal.le.zero) then
           print *,"voltotal bust"
           stop
          endif

          call Box_volumeFAST(bfact,dx,u_xsten_updatecell,nhalf, &
           volcell,cencell,SDIM)

          do u_imaterial=1,nmat
           unsplit_lsnew(u_imaterial)=lsmat(u_imaterial)/voltotal
          enddo

          iten=iten_crit
          ireverse=ireverse_crit
          im_dest=im_dest_crit
          im_source=im_source_crit 

          if (iten.eq.iten_outer) then
           ! do nothing
          else
           print *,"iten must be equal to iten_outer"
           stop
          endif

          if (interface_near(iten+ireverse*nten).ne.1) then
           print *,"interface_near invalid"
           stop
          endif

          Tgamma_default=saturation_temp(iten+ireverse*nten)
          Ygamma_default=one

          Tsat_flag=NINT(TgammaFAB(D_DECL(i,j,k),iten))
          if (ireverse.eq.0) then
            ! do nothing
          else if (ireverse.eq.1) then
            Tsat_flag=-Tsat_flag
          else
            print *,"ireverse invalid"
            stop
          endif

          LL=latent_heat(iten+ireverse*nten)
          local_freezing_model=freezing_model(iten+ireverse*nten)
          distribute_from_targ=distribute_from_target(iten+ireverse*nten)
          mass_frac_id=mass_fraction_id(iten+ireverse*nten)
          
           ! TgammaFAB has valid values corresponding to the
           ! given value for "ireverse". 
          if ((Tsat_flag.eq.1).or.(Tsat_flag.eq.2)) then
            Tgamma_default=TgammaFAB(D_DECL(i,j,k), &
             nten+(iten-1)*ncomp_per_tsat+1)
            if ((mass_frac_id.ge.1).and. &
                (mass_frac_id.le.num_species_var)) then
             Ygamma_default=TgammaFAB(D_DECL(i,j,k), &
              nten+(iten-1)*ncomp_per_tsat+2)
            else if (mass_frac_id.eq.0) then
             Ygamma_default=one
            else
             print *,"mass_frac_id invalid"
             stop
            endif
           ! TgammaFAB has valid values corresponding to the
           ! given value for 1-ireverse
          else if ((Tsat_flag.eq.-1).or. &
                   (Tsat_flag.eq.-2)) then
            ! do nothing
          else if (Tsat_flag.eq.0) then
            ! do nothing
          else
            print *,"Tsat_flag invalid"
            stop
          endif

          do u_imaterial=1,nmat

           vofcomp_raw=(u_imaterial-1)*ngeom_raw+1

           tempvfrac=volmat(u_imaterial)/voltotal
           if ((tempvfrac.ge.EBVOFTOL).and.(tempvfrac.le.1.1d0)) then
            if (tempvfrac.gt.one) then
             tempvfrac=one
            endif
           else if ((tempvfrac.ge.zero).and.(tempvfrac.le.EBVOFTOL)) then
            tempvfrac=zero
           else
            print *,"tempvfrac bust1 tempvfrac=",tempvfrac
            stop
           endif

           do udir=1,SDIM
            if ((tempvfrac.gt.zero).and.(tempvfrac.le.one)) then
             tempcen(udir)=cenmat(udir,u_imaterial)/volmat(u_imaterial)- &
              cencell(udir)
            else if (tempvfrac.eq.zero) then
             tempcen(udir)=zero
            else
             print *,"tempvfrac bust2 tempvfrac=",tempvfrac
             print *,"udir=",udir
             stop
            endif
           enddo ! udir
           
           unsplit_snew(vofcomp_raw)=tempvfrac
           do udir=1,SDIM
            unsplit_snew(vofcomp_raw+udir)=tempcen(udir)
           enddo

           if ((tempvfrac.gt.zero).and.(tempvfrac.le.one)) then
            unsplit_density(u_imaterial)= &
             density_mat(u_imaterial)/volmat(u_imaterial)
            unsplit_temperature(u_imaterial)= &
             temperature_mat(u_imaterial)/volmat(u_imaterial)
            do ispec=1,num_species_var
             unsplit_species((ispec-1)*nmat+u_imaterial)= &
              species_mat((ispec-1)*nmat+u_imaterial)/volmat(u_imaterial)
            enddo
           else if (tempvfrac.eq.zero) then
            unsplit_density(u_imaterial)=fort_denconst(u_imaterial)
            unsplit_temperature(u_imaterial)=Tgamma_default
            do ispec=1,num_species_var
             unsplit_species((ispec-1)*nmat+u_imaterial)=Ygamma_default
            enddo
           else
            print *,"tempvfrac bust3 tempvfrac=",tempvfrac
            stop
           endif

          enddo  ! u_imaterial=1..nmat

          vcompsrc_snew=num_materials_vel*(SDIM+1)+ &
           nmat*num_state_material+(im_source-1)*ngeom_raw+1
          vcompdst_snew=num_materials_vel*(SDIM+1)+ &
           nmat*num_state_material+(im_dest-1)*ngeom_raw+1

          do u_imaterial=1,nmat
           vofcomp_recon=(u_imaterial-1)*ngeom_recon+1
           oldvfrac(u_imaterial)=recon(D_DECL(i,j,k),vofcomp_recon)
           vofcomp_raw=(u_imaterial-1)*ngeom_raw+1
           newvfrac(u_imaterial)=unsplit_snew(vofcomp_raw)
           do dir=1,SDIM
            old_centroid(u_imaterial,dir)= &
              recon(D_DECL(i,j,k),vofcomp_recon+dir)+cengrid(dir)
            new_centroid(u_imaterial,dir)= &
              unsplit_snew(vofcomp_raw+dir)+cengrid(dir)
           enddo
          enddo ! u_imaterial=1,nmat

          do u_imaterial=1,nmat*(1+SDIM)
           oldLS_point(u_imaterial)=LSold(D_DECL(i,j,k),u_imaterial)
          enddo
          call normalize_LS_normals(nmat,oldLS_point)

          vofcomp_raw_dest=(im_dest-1)*ngeom_raw+1
          vofcomp_recon_source=(im_source-1)*ngeom_recon+1

          do dir=1,SDIM
           if ((newvfrac(im_dest).gt.zero).and. &
               (newvfrac(im_dest).le.one)) then
            new_centroid(im_dest,dir)= &
              unsplit_snew(vofcomp_raw_dest+dir)+cengrid(dir)

            ! all the source material can be converted into
            ! destination material.
           else if ((newvfrac(im_dest).eq.zero).and. &
                    (oldvfrac(im_source).gt.zero)) then
            new_centroid(im_dest,dir)= &
              recon(D_DECL(i,j,k),vofcomp_recon_source+dir)+cengrid(dir)
           else if ((newvfrac(im_dest).eq.zero).and. &
                    (abs(oldvfrac(im_source)).le.VOFTOL)) then
            new_centroid(im_dest,dir)=cengrid(dir)
           else
            print *,"newvfrac(im_dest) or oldvfrac(im_source) invalid"
            stop
           endif
          enddo ! dir=1..sdim

          call get_primary_material(unsplit_lsnew,nmat,im_primary_new)
          call get_primary_material(oldLS_point,nmat,im_primary_old)
          call combine_solid_VOF(newvfrac,nmat,solid_vof_new)
          call combine_solid_VOF(oldvfrac,nmat,solid_vof_old)

          away_from_interface=0
          if ((is_rigid(nmat,im_primary_new).eq.1).or. &
              (is_rigid(nmat,im_primary_old).eq.1)) then
           away_from_interface=1
          endif

          if ((solid_vof_new.ge.half).or. &
              (solid_vof_old.ge.half)) then
           away_from_interface=1
          endif
          if (((abs(unsplit_lsnew(im_source)).gt.dxmaxLS).or. &
               (abs(unsplit_lsnew(im_dest)).gt.dxmaxLS)).and. &
              ((abs(oldLS_point(im_dest)).gt.dxmaxLS).or. &
               (abs(oldLS_point(im_source)).gt.dxmaxLS))) then
           away_from_interface=1
          endif

          if (DEBUG_ACTIVE_CELL.eq.1) then
           if ((i.eq.DEBUG_I).and.(j.eq.DEBUG_J)) then
            print *,"im_source,im_dest,away_from_interface ", &
               im_source,im_dest,away_from_interface
           endif
          endif

          if (away_from_interface.eq.1) then
           newvfrac(im_source)=oldvfrac(im_source)
           newvfrac(im_dest)=oldvfrac(im_dest)
           do udir=1,SDIM 
            new_centroid(im_source,udir)=old_centroid(im_source,udir)
            new_centroid(im_dest,udir)=old_centroid(im_dest,udir)
           enddo
          else if (away_from_interface.eq.0) then

            ! only the interface changing phase should move
           LS_MAX_fixed=-1.0D+20
           do im_local=1,nmat
            lsdata(im_local)=LSold(D_DECL(i,j,k),im_local)
            if ((im_local.eq.im_source).or.(im_local.eq.im_dest)) then
             ! do nothing
            else if ((im_local.ge.1).and.(im_local.le.nmat)) then
             if (is_rigid(nmat,im_local).eq.1) then
              ! do nothing
             else if (is_rigid(nmat,im_local).eq.0) then  
              LS_MAX_fixed=max(LS_MAX_fixed,lsdata(im_local))
             else
              print *,"is_rigid(nmat,im_local) invalid"
              stop
             endif
            else
             print *,"im_local bust"
             stop
            endif
           enddo ! im_local=1..nmat
            ! make sure the level set functions 
            ! LSnew(im_source), LSnew(im_dest),
            ! LS_old(im_fluid_rest) 
            ! are tessellating.
           lsdata(im_source)=half*(unsplit_lsnew(im_source)- &
             max(LS_MAX_fixed,unsplit_lsnew(im_dest)))
           lsdata(im_dest)=half*(unsplit_lsnew(im_dest)- &
             max(LS_MAX_fixed,unsplit_lsnew(im_source)))

           LSnew(D_DECL(i,j,k),im_source)=lsdata(im_source)
           LSnew(D_DECL(i,j,k),im_dest)=lsdata(im_dest)

          else
           print *,"away_from_interface invalid"
           stop
          endif

           ! The new volume fractions should be tessellating, while
           ! at the same time materials not involved in the phase change
           ! should not have a change in volume.
          fixed_vfrac_sum=zero
          do udir=1,SDIM 
           fixed_centroid_sum(udir)=zero
          enddo
          do im_local=1,nmat
           if ((im_local.eq.im_source).or. &
               (im_local.eq.im_dest)) then
            ! do nothing
           else if ((im_local.ge.1).and.(im_local.le.nmat)) then
            newvfrac(im_local)=oldvfrac(im_local)
            do udir=1,SDIM 
             new_centroid(im_local,udir)=old_centroid(im_local,udir)
            enddo
            if (is_rigid(nmat,im_local).eq.0) then
             fixed_vfrac_sum=fixed_vfrac_sum+newvfrac(im_local)
             do udir=1,SDIM 
              fixed_centroid_sum(udir)=fixed_centroid_sum(udir)+ &
                newvfrac(im_local)*new_centroid(im_local,udir)
             enddo
            else if (is_rigid(nmat,im_local).eq.1) then
             ! ignore, solids are embedded
            else
             print *,"is_rigid invalid"
             stop
            endif
           else
            print *,"im_local became corrupt"
            stop
           endif
          enddo ! im_local=1..nmat

          if ((fixed_vfrac_sum.ge.one-VOFTOL).and. &
              (fixed_vfrac_sum.le.one+VOFTOL)) then
           avail_vfrac=zero
          else if ((fixed_vfrac_sum.ge.-VOFTOL).and. &
                   (fixed_vfrac_sum.le.zero)) then
           avail_vfrac=one
          else if ((fixed_vfrac_sum.gt.zero).and. &
                   (fixed_vfrac_sum.le.one-VOFTOL)) then
           avail_vfrac=one-fixed_vfrac_sum
          else
           print *,"fixed_vfrac_sum invalid"
           stop
          endif

          if (DEBUG_ACTIVE_CELL.eq.1) then
           if ((i.eq.DEBUG_I).and.(j.eq.DEBUG_J)) then
            print *,"im_source,im_dest,avail_vfrac,oldvfrac(im_dest) ", &
               im_source,im_dest,avail_vfrac,oldvfrac(im_dest)
           endif
          endif

          if ((avail_vfrac.eq.zero).or. &
              (avail_vfrac.le.oldvfrac(im_dest)+VOFTOL)) then
           dF=zero
           newvfrac(im_dest)=oldvfrac(im_dest)
           newvfrac(im_source)=oldvfrac(im_source)
           do udir=1,SDIM 
            new_centroid(im_source,udir)=old_centroid(im_source,udir)
            new_centroid(im_dest,udir)=old_centroid(im_dest,udir)
           enddo
           LSnew(D_DECL(i,j,k),im_source)=LSold(D_DECL(i,j,k),im_source)
           LSnew(D_DECL(i,j,k),im_dest)=LSold(D_DECL(i,j,k),im_dest)
          else if ((avail_vfrac.gt.zero).and. &
                   (avail_vfrac.le.one)) then
           dFdst=(newvfrac(im_dest)-oldvfrac(im_dest))
           dFsrc=(oldvfrac(im_source)-newvfrac(im_source))

           if (DEBUG_ACTIVE_CELL.eq.1) then
            if ((i.eq.DEBUG_I).and.(j.eq.DEBUG_J)) then
             print *,"im_source,im_dest,dFdst,dFsrc ", &
               im_source,im_dest,dFdst,dFsrc
            endif
           endif

           if ((dFdst.le.zero).and. &
               (dFsrc.le.zero)) then
            dF=zero
            newvfrac(im_dest)=oldvfrac(im_dest)
            newvfrac(im_source)=oldvfrac(im_source)
            do udir=1,SDIM 
             new_centroid(im_source,udir)=old_centroid(im_source,udir)
             new_centroid(im_dest,udir)=old_centroid(im_dest,udir)
            enddo
            LSnew(D_DECL(i,j,k),im_source)=LSold(D_DECL(i,j,k),im_source)
            LSnew(D_DECL(i,j,k),im_dest)=LSold(D_DECL(i,j,k),im_dest)
           else if ((dFdst.gt.zero).or. &
                    (dFsrc.gt.zero)) then

            ! mass fraction equation:
            ! in a given cell with m species.
            ! mass= sum_i=1^m  Y_i overall_mass = 
            !     = sum_i=1^m  density_i F_i V_cell
            ! F_i = volume fraction of material i.
            ! (rho Y_i)_t + div (rho u Y_i) = div rho D_i  grad Y_i
            ! (Y_i)_t + div (u Y_i) = div rho D_i  grad Y_i/rho

             ! if there are 2 materials, dFdst=dFsrc, but if there
             ! are 3 materials or more, the larger value is to be
             ! trusted over the smaller, due to interference from the
             ! extra materials (which hypothetically, should be fixed)
            im_trust=0
            im_distrust=0
            if (dFdst.ge.dFsrc) then
             im_trust=im_dest
             im_distrust=im_source
             dF=dFdst
            else if (dFdst.le.dFsrc) then
             im_trust=im_source
             im_distrust=im_dest
             dF=dFsrc
            else
             print *,"dFdst or dFsrc bust"
             stop
            endif

            if ((LL.gt.zero).or. & !evaporation,boiling,melting,cavitation
                (LL.lt.zero)) then !freezing, condensation
             dF=min(dF,oldvfrac(im_source))
             dF=min(dF,avail_vfrac-oldvfrac(im_dest))
            else
             print *,"LL invalid"
             stop
            endif

            newvfrac(im_dest)=oldvfrac(im_dest)+dF
            newvfrac(im_source)=oldvfrac(im_source)-dF
 
            if (dF.eq.zero) then
             newvfrac(im_dest)=oldvfrac(im_dest)
             newvfrac(im_source)=oldvfrac(im_source)
             do udir=1,SDIM 
              new_centroid(im_source,udir)=old_centroid(im_source,udir)
              new_centroid(im_dest,udir)=old_centroid(im_dest,udir)
             enddo
             LSnew(D_DECL(i,j,k),im_source)=LSold(D_DECL(i,j,k),im_source)
             LSnew(D_DECL(i,j,k),im_dest)=LSold(D_DECL(i,j,k),im_dest)
            else if (dF.gt.zero) then
             if (fixed_vfrac_sum.gt.VOFTOL) then
              do udir=1,SDIM 
               if (newvfrac(im_distrust).gt.VOFTOL) then
                new_centroid(im_distrust,udir)= &
                 (cengrid(udir)- &
                  (new_centroid(im_trust,udir)*newvfrac(im_trust)+ &
                   fixed_centroid_sum(udir)))/newvfrac(im_distrust)
               else if (abs(newvfrac(im_distrust)).le.VOFTOL) then
                new_centroid(im_distrust,udir)=cengrid(udir)
               else
                print *,"newvfrac(im_distrust) invalid"
                stop
               endif
              enddo ! udir=1..sdim
             else if (abs(fixed_vfrac_sum).le.VOFTOL) then
              ! do nothing
             else
              print *,"fixed_vfrac_sum invalid"
              stop
             endif
                
            else
             print *,"dF invalid"
             stop
            endif
           else
            print *,"dF became corrupt1 dF=",dF
            stop
           endif
          else
           print *,"avail_vfrac invalid"
           stop
          endif

          if (DEBUG_ACTIVE_CELL.eq.1) then
           if ((i.eq.DEBUG_I).and.(j.eq.DEBUG_J)) then
            print *,"im_source,im_dest,dF ", &
               im_source,im_dest,dF
            stop
           endif
          endif

           ! 1. the MAC and cell velocity field should be extrapolated
           !    from the old destination material side into the cells/faces
           !    that were swept by the interface. (SWEPTFACTOR>0.0)
           ! 2. the destination temperature and species should be set 
           !    to the interface values at destination centroids which were
           !    swept by the interface. (SWEPTFACTOR_centroid==1)
           ! 3. the source and destination temperature and species should 
           !    be recalculated in non-swept partial cells.  They should 
           !    be interpolated from the old supermesh grid to the new
           !    taking into consideration the interface boundary condition.
          if (dF.gt.EBVOFTOL) then

           do iprobe=1,2

            if (iprobe.eq.1) then ! source
             im_probe=im_source
            else if (iprobe.eq.2) then ! dest
             im_probe=im_dest
            else
             print *,"iprobe invalid"
             stop
            endif

            dencomp_probe=(im_probe-1)*num_state_material+1
            tcomp_probe=dencomp_probe+1
            if ((mass_frac_id.ge.1).and. &
                (mass_frac_id.le.num_species_var)) then
             mfrac_comp_probe=tcomp_probe+mass_frac_id
             ispec_probe=(mass_frac_id-1)*nmat+im_probe
            else if (mass_frac_id.eq.0) then
             mfrac_comp_probe=0
             ispec_probe=0
            else
             print *,"mass_frac_id invalid"
             stop
            endif

            if (newvfrac(im_probe).ge.EBVOFTOL) then
             temperature_new(iprobe)=unsplit_temperature(im_probe)
             if (mfrac_comp_probe.eq.0) then ! no species vars
              species_new(iprobe)=one
             else if (mfrac_comp_probe.gt.0) then
              species_new(iprobe)= &
                unsplit_species((mass_frac_id-1)*nmat+im_probe)
             else
              print *,"mfrac_comp_probe invalid"
              stop
             endif
            else if ((newvfrac(im_probe).le.EBVOFTOL).and. &
                     (newvfrac(im_probe).ge.-EBVOFTOL)) then
             temperature_new(iprobe)=Tgamma_default
             species_new(iprobe)=Ygamma_default
            else
             print *,"newvfrac invalid"
             stop
            endif

            mtype=fort_material_type(im_probe)
            if (mtype.eq.0) then

              ! the density was extrapolated after CISL advection, so 
              ! that a valid value exists even if F_new_unsplit=0.0 
             if (constant_density_all_time(im_probe).eq.1) then
              density_old(iprobe)=fort_denconst(im_probe)
             else if (constant_density_all_time(im_probe).eq.0) then
              density_old(iprobe)=EOS(D_DECL(i,j,k),dencomp_probe)
             else
              print *,"constant_density_all_time(im_probe) invalid"
              stop
             endif

             if (density_old(iprobe).gt.zero) then
              ! do nothing
             else
              print *,"density_old(iprobe) invalid"
              stop
             endif

            else if ((mtype.ge.1).and.(mtype.le.MAX_NUM_EOS)) then
             print *,"only spatially uniform density phase change allowed"
             stop
            else
             print *,"mtype invalid"
             stop
            endif

            if (oldvfrac(im_probe).ge.EBVOFTOL) then
             temperature_old(iprobe)=EOS(D_DECL(i,j,k),tcomp_probe)
             if (mfrac_comp_probe.eq.0) then ! no species var
              species_old(iprobe)=one
             else if (mfrac_comp_probe.gt.0) then
              species_old(iprobe)=EOS(D_DECL(i,j,k),mfrac_comp_probe)
             else
              print *,"mfrac_comp_probe invalid"
              stop
             endif
            else if ((oldvfrac(im_probe).le.EBVOFTOL).and. &
                     (oldvfrac(im_probe).ge.-EBVOFTOL)) then
             temperature_old(iprobe)=Tgamma_default
             species_old(iprobe)=Ygamma_default
            else
             print *,"oldvfrac invalid"
             stop
            endif

           enddo ! iprobe=1,2

           if (LL.gt.zero) then !evaporation
             im_vapor=im_dest
             im_condensed=im_source
             iprobe_vapor=2 ! dest
             iprobe_condensed=1 !source
           else if (LL.lt.zero) then ! condensation
             im_vapor=im_source
             im_condensed=im_dest
             iprobe_vapor=1 ! source
             iprobe_condensed=2 ! dest
           else
             print *,"im_vapor invalid"
             stop
           endif

            ! at CONSTANT DENSITY:
            ! "vapor" is what material transforms TO if system is HEATED.
            ! "condensed" is what material transforms FROM if system is HEATED.
            ! "vapor" is what material transforms FROM if system is COOLED.
            ! "condensed" is what material transforms TO if system is COOLED.
            ! 
           if (local_freezing_model.eq.0) then !Stefan model
             vapor_den=density_old(iprobe_vapor)
             condensed_den=density_old(iprobe_condensed)
           else if (local_freezing_model.eq.1) then !source term
             vapor_den=density_old(iprobe_vapor)
             condensed_den=density_old(iprobe_condensed)
           else if (local_freezing_model.eq.2) then !hydrate
             vapor_den=density_old(iprobe_vapor)
             condensed_den=density_old(iprobe_condensed)
           else if (local_freezing_model.eq.3) then !wildfire
             vapor_den=density_old(iprobe_vapor)
             condensed_den=density_old(iprobe_condensed)
           else if ((local_freezing_model.eq.4).or. & !Tannasawa or Schrage
                    (local_freezing_model.eq.5).or. & !Stefan evap/cond model
                    (local_freezing_model.eq.6).or. & !Palmore/Desjardins
                    (local_freezing_model.eq.7)) then !Cavitation

            if ((mass_frac_id.ge.1).and. &
                (mass_frac_id.le.num_species_var)) then

             vapor_den=density_old(iprobe_vapor)
             condensed_den=density_old(iprobe_condensed)

             if ((vapor_den.gt.zero).and.(condensed_den.gt.zero)) then
              ! do nothing
             else
              print *,"vapor_den or condensed_den invalid"
              stop
             endif

            else
             print *,"mass_frac_id invalid"
             stop
            endif

           else
            print *,"local_freezing_model invalid 1"
            stop
           endif

           temp_mix_new(2)=temperature_new(2)
           temp_mix_new(1)=temperature_new(1)

           if (LL.gt.zero) then ! evaporation or boiling

            temp_new_vfrac=oldvfrac(im_dest)+dF
            if (temp_new_vfrac.gt.EBVOFTOL) then 
             
             if (1.eq.0) then
              print *,"i,j,k,temp_new_vfrac ",i,j,k,temp_new_vfrac
              print *,"denold,oldF,vapor_den,dF ", &
               density_old(2),oldvfrac(im_dest),vapor_den,dF
             endif

             mtype=fort_material_type(im_vapor)
             if (mtype.eq.0) then
              mass_frac_new(2)=Ygamma_default
             else if ((mtype.ge.1).and.(mtype.le.MAX_NUM_EOS)) then
              print *,"only spatially uniform density phase change allowed"
              stop 
             else
              print *,"mtype invalid"
              stop
             endif
            else if (temp_new_vfrac.le.EBVOFTOL) then
             mass_frac_new(2)=Ygamma_default
            else
             print *,"temp_new_vfrac invalid"
             stop
            endif
           
            if (im_source.eq.im_condensed) then 
             ! do nothing
            else
             print *,"expecting im_source==im_condensed"
             stop
            endif
            mass_frac_new(1)=Ygamma_default

           else if (LL.lt.zero) then ! condensation

            if (im_dest.eq.im_condensed) then 
             ! do nothing
            else
             print *,"expecting im_dest==im_condensed"
             stop
            endif

            temp_new_vfrac=oldvfrac(im_source)-dF
            if (temp_new_vfrac.gt.EBVOFTOL) then
             mass_frac_new(1)=Ygamma_default
            else if (temp_new_vfrac.le.EBVOFTOL) then
             mass_frac_new(1)=Ygamma_default
            else
             print *,"temp_new_vfrac invalid"
             stop
            endif
            mass_frac_new(2)=Ygamma_default
           else
            print *,"LL invalid"
            stop
           endif

            ! dF>EBVOFTOL in this section, so now let's see
            ! if the cell center has been swept and also the
            ! destination centroid.
            ! for GFM, if a cell center is not occupied by the
            ! destination material at tn but is occupied at tnp1, then
            ! (theta^{n+1} - theta_{I})/(tnp1 - tswept) = L(theta^{n+1}
            ! tswept is the crossing time.
            ! Suppose 1D and GFM:
            !  at t=tn      F=Fn
            !  at t=tswept  F=1/2
            !  at t=tnp1    F=Fnp1
            !  F(t)=(Fnp1 - Fn)/(tnp1-tn)  * (t-tnp1) + Fnp1
            !  1/2 = dF/(tnp1-tn)   * (tswept-tnp1)  + Fnp1
            !  dF/(tnp1-tn)  * (tswept-tnp1) = Fnp1-1/2
            !  SWEPTFACTOR=(tswept-tnp1)/(tnp1-tn)  = (Fnp1-1/2)/dF
            !  new approach (motivated by supermesh):
            !  1. given an estimate of phase change velocity "USTEFAN",
            !     dt is chosen such that USTEFAN * dt <= dx/4
            !     This assures that a swept cell will NOT be full at tnp1.
            !  2. Strategy is this: 
            !     a) we have the signed distance function at tn.
            !     b) let x^* be either (i) cell center if GFM, or (ii)
            !        cell centroid at tnp1 of the destination material.
            !     c) if phi(tn,x^*)<0 then:
            !         (i) find MOF reconstruction and determine if
            !             phi^reconstruct(tnp1,x^*)>0
            !             if yes, then 
            !             SWEPTFACTOR=1-(phi_np1/(phi_np1-phi_n)     

            ! this is a placeholder for what used to be implemented.
            ! Now, the swept factor is calculated much more precisely.
            !
           SWEPTFACTOR_GFM=one
           if ((oldvfrac(im_dest).lt.half).and. &
               (newvfrac(im_dest).gt.half)) then
            SWEPTFACTOR_GFM=(newvfrac(im_dest)-half)/dF
           else if ((oldvfrac(im_dest).ge.half).or. &
                    (newvfrac(im_dest).le.half)) then
            SWEPTFACTOR_GFM=one
           else
            print *,"oldvfrac or newvfrac invalid"
            stop
           endif

            ! => new_centroid(im_dest,udir) (absolute coord)
           do u_im=1,nmat
            vofcomp_recon=(u_im-1)*ngeom_recon+1

            mofdata(vofcomp_recon)=oldvfrac(u_im)
            mofdata_new(vofcomp_recon)=newvfrac(u_im)

            do udir=1,SDIM 
             mofdata(vofcomp_recon+udir)= &
               old_centroid(u_im,udir)-cengrid(udir)

             mofdata_new(vofcomp_recon+udir)= &
               new_centroid(u_im,udir)-cengrid(udir)
            enddo
            mofdata(vofcomp_recon+SDIM+1)= &
             recon(D_DECL(i,j,k),vofcomp_recon+SDIM+1) !ord

            mofdata_new(vofcomp_recon+SDIM+1)=0  ! placeholder order

            mofdata(vofcomp_recon+2*SDIM+2)= & 
             recon(D_DECL(i,j,k),vofcomp_recon+2*SDIM+2)  !intercept

            mofdata_new(vofcomp_recon+2*SDIM+2)=zero  ! placeholder intercept

            do udir=1,SDIM
             mofdata(vofcomp_recon+SDIM+1+udir)= &
              recon(D_DECL(i,j,k),vofcomp_recon+SDIM+1+udir) !slope

             mofdata_new(vofcomp_recon+SDIM+1+udir)=zero ! placeholder slope
            enddo ! udir
           enddo ! u_im=1..nmat

           mof_verbose=0
           use_ls_data=0
           continuous_mof_parm=0

            ! LS=n dot (x-x0)+intercept
           call multimaterial_MOF( &
            bfact,dx,u_xsten_updatecell,nhalf, &
            mof_verbose, &
            use_ls_data, &
            LS_stencil, &
            geom_xtetlist(1,1,1,tid+1), &
            geom_xtetlist_old(1,1,1,tid+1), &
            nmax, &
            nmax, &
            mofdata_new, &
            multi_centroidA, &
            continuous_mof_parm, &
            cmofsten, &
            nmat,SDIM,202)

            ! xPOINT_supermesh is needed in order to determine
            ! whether to interpolate old temperature and mass fraction
            ! data from the old supermesh to the new supermesh.
           do udir=1,SDIM
            xPOINT_supermesh(udir)=new_centroid(im_dest,udir)
           enddo

            ! now we check if new_centroid(im_dest,dir) is in the
            ! old dest material.
            ! If new_centroid(im_dest) in the old dest material, then
            ! the cell centroid HAS NOT been swept and the temperature
            ! at this point needs to be interpolated from the t=tn
            ! im_dest supermesh.
           SWEPTFACTOR_centroid=0
           if ((newvfrac(im_dest).gt.zero).and. &
               (newvfrac(im_dest).le.one+EBVOFTOL)) then

            tessellate=3
            call multi_get_volumePOINT( &
              tessellate, &
              bfact,dx, &
              u_xsten_updatecell,nhalf, &  ! absolute coordinate system
              mofdata, &
              xPOINT_supermesh, & ! absolute coordinate system
              im_old_crit,nmat,SDIM)

             ! im_old_crit=material at t=tn that occupies the new 
             ! centroid location
            if ((im_old_crit.eq.im_dest).or. &
                (oldvfrac(im_dest).ge.half)) then
             ! do nothing
            else if ((im_old_crit.ge.1).and. &
                     (im_old_crit.le.nmat).and. &
                     (im_old_crit.ne.im_dest).and. &
                     (oldvfrac(im_dest).le.half)) then
             SWEPTFACTOR_centroid=1

             if (1.eq.0) then
              print *,"setting SWEPTFACTOR_centroid=1"
              print *,"im_dest= ",im_dest
              print *,"xPOINT_supermesh= ", &
                xPOINT_supermesh(1),xPOINT_supermesh(2),xPOINT_supermesh(SDIM)
              print *,"xPOINT_supermesh(1)-xsten(0,1)= ", &
                      xPOINT_supermesh(1)-u_xsten_updatecell(0,1)
              print *,"oldvfrac(im_dest)=",oldvfrac(im_dest)
              print *,"newvfrac(im_dest)=",newvfrac(im_dest)
              print *,"im_old_crit=",im_old_crit
              print *,"i,j,k ",i,j,k
              print *,"u_xsten_updatecell xlo ",u_xsten_updatecell(-1,1)
              print *,"u_xsten_updatecell xhi ",u_xsten_updatecell(1,1)
              print *,"u_xsten_updatecell ",u_xsten_updatecell(0,1), &
                      u_xsten_updatecell(0,2), &
                      u_xsten_updatecell(0,SDIM)
              do u_im=1,nmat
               do udir=1,ngeom_recon
                print *,"im,mofcomp,mofdata ",u_im,udir, &
                        mofdata((u_im-1)*ngeom_recon+udir)
               enddo
              enddo
             endif
            else
             print *,"im_old_crit invalid"
             stop
            endif 

           else
            print *,"expecting newvfrac(im_dest)>0 since dF>0"
            stop
           endif

           interp_to_new_supermesh=1

           if (oldLS_point(im_dest).ge.zero) then
            do udir=1,SDIM
             old_nrm(udir)=oldLS_point(nmat+(im_source-1)*SDIM+udir)
             old_xI(udir)=u_xsten_updatecell(0,udir)- &
               oldLS_point(im_source)*old_nrm(udir)
            enddo
           else if (oldLS_point(im_source).ge.zero) then
            do udir=1,SDIM
             old_nrm(udir)=oldLS_point(nmat+(im_dest-1)*SDIM+udir)
             old_xI(udir)=u_xsten_updatecell(0,udir)- &
               oldLS_point(im_dest)*old_nrm(udir)
            enddo
           else if ((oldLS_point(im_dest).le.zero).and. &
                    (oldLS_point(im_source).le.zero)) then
            interp_to_new_supermesh=0
           else
            print *,"oldLS_point(im_dest) or oldLS_point(im_source) invalid"
            stop
           endif

           do iprobe=1,2

            if (iprobe.eq.1) then ! source
             im_probe=im_source
            else if (iprobe.eq.2) then ! dest
             im_probe=im_dest
            else
             print *,"iprobe invalid"
             stop
            endif
            vofcomp_recon=(im_probe-1)*ngeom_recon+1

            dencomp_probe=(im_probe-1)*num_state_material+1
            tcomp_probe=dencomp_probe+1
            if ((mass_frac_id.ge.1).and. &
                (mass_frac_id.le.num_species_var)) then
             mfrac_comp_probe=tcomp_probe+mass_frac_id
            else if (mass_frac_id.eq.0) then
             mfrac_comp_probe=0
            else
             print *,"mass_frac_id invalid"
             stop
            endif
            base_index=num_materials_vel*(SDIM+1)

            if (iprobe.eq.1) then ! source
             delta_mass_local(iprobe)=-density_old(iprobe)*dF
            else if (iprobe.eq.2) then ! dest
             delta_mass_local(iprobe)=density_old(iprobe)*dF
            else
             print *,"iprobe invalid"
             stop
            endif
            if (iprobe.eq.1) then ! source
             if (delta_mass_local(iprobe).le.zero) then
              ! do nothing
             else
              print *,"delta_mass_local invalid"
              stop
             endif
            else if (iprobe.eq.2) then ! dest
             if (delta_mass_local(iprobe).ge.zero) then
              ! do nothing
             else
              print *,"delta_mass_local invalid"
              stop
             endif
            else
             print *,"iprobe invalid"
             stop
            endif

            if ((SWEPTFACTOR_centroid.eq.1).and. &
                (iprobe.eq.2)) then ! destination
             temp_mix_new(iprobe)=Tgamma_default
             mass_frac_new(iprobe)=Ygamma_default
            else if ((SWEPTFACTOR_centroid.eq.0).or. &
                     (iprobe.eq.1)) then ! source
             ! here we interpolate from old supermesh to new, making
             ! sure to take into account Tgamma_default and Ygamma_default
             if (interp_to_new_supermesh.eq.1) then

              do i1=-1,1
              do j1=-1,1
              do k1=klosten,khisten
               call gridsten_level(xsten_ofs,i+i1,j+j1,k+k1,level,nhalf)
               call Box_volumeFAST(bfact,dx,xsten_ofs,nhalf, &
                volcell_ofs,cencell_ofs,SDIM)
               do udir=1,SDIM
                XC_sten(D_DECL(i1,j1,k1),udir)= &
                 recon(D_DECL(i+i1,j+j1,k+k1),vofcomp_recon+udir)+ &
                 cencell_ofs(udir)
               enddo
               VF_sten(D_DECL(i1,j1,k1))= &
                recon(D_DECL(i+i1,j+j1,k+k1),vofcomp_recon)
               LS_sten(D_DECL(i1,j1,k1))= &
                LSold(D_DECL(i+i1,j+j1,k+k1),im_probe)
               temperature_sten(D_DECL(i1,j1,k1))= &
                EOS(D_DECL(i+i1,j+j1,k+k1),tcomp_probe)
               if (mfrac_comp_probe.eq.0) then
                massfrac_sten(D_DECL(i1,j1,k1))=Ygamma_default
               else 
                massfrac_sten(D_DECL(i1,j1,k1))= &
                 EOS(D_DECL(i+i1,j+j1,k+k1),mfrac_comp_probe)
               endif
              enddo
              enddo
              enddo ! i1,j1,k1

              if ((oldvfrac(im_probe).le.VOFTOL).and. &
                  (oldvfrac(im_probe).ge.-VOFTOL)) then
               temp_mix_new(iprobe)=Tgamma_default
               mass_frac_new(iprobe)=Ygamma_default
              else if ((oldvfrac(im_probe).ge.VOFTOL).and. &
                       (oldvfrac(im_probe).le.one+VOFTOL)) then

               DATA_FLOOR=zero
               combine_flag=0
               nsolve_interp=1
               do udir=1,SDIM
                xtarget_interp(udir)=new_centroid(im_probe,udir)
               enddo

               call center_centroid_interchange( &
                DATA_FLOOR, &
                nsolve_interp, &
                combine_flag, & !0=>centroid -> center   1=>center->centroid
                interp_to_new_supermesh, &
                bfact, &
                level, &
                finest_level, &
                dx,xlo, &
                u_xsten_updatecell,nhalf, &
                temperature_sten, &
                XC_sten, &
                old_xI, &
                xtarget_interp, &
                VF_sten, &
                LS_sten, &
                Tgamma_default, &
                temp_mix_new(iprobe))

               if (1.eq.0) then
                print *,"correcting mass fraction"
                print *,"i,j,k ",i,j,k
                print *,"iprobe=",iprobe
                print *,"im_probe=",im_probe
                print *,"oldvfrac(im_probe) ",oldvfrac(im_probe)
                print *,"newvfrac(im_probe) ",newvfrac(im_probe)
               endif

               call center_centroid_interchange( &
                DATA_FLOOR, &
                nsolve_interp, &
                combine_flag, & !0=>centroid -> center   1=>center->centroid
                interp_to_new_supermesh, &
                bfact, &
                level, &
                finest_level, &
                dx,xlo, &
                u_xsten_updatecell,nhalf, &
                massfrac_sten, &
                XC_sten, &
                old_xI, &
                xtarget_interp, &
                VF_sten, &
                LS_sten, &
                Ygamma_default, &
                mass_frac_new(iprobe))

              else
               print *,"oldvfrac(im_probe) invalid"
               stop
              endif

             else if (interp_to_new_supermesh.eq.0) then
              ! do nothing
             else
              print *,"interp_to_new_supermesh invalid"
              stop
             endif

            else
             print *,"SWEPTFACTOR_centroid invalid"
             stop
            endif

            snew(D_DECL(i,j,k),base_index+tcomp_probe)= &
                    temp_mix_new(iprobe)
            if (mfrac_comp_probe.eq.0) then
             ! do nothing
            else if (mfrac_comp_probe.gt.0) then
             snew(D_DECL(i,j,k),base_index+mfrac_comp_probe)=  &
              mass_frac_new(iprobe)
            else
             print *,"mfrac_comp_probe invalid"
             stop
            endif

            if (temp_mix_new(iprobe).ge.TEMPERATURE_FLOOR) then
             ! do nothing
            else
             print *,"temp_mix_new invalid"
             print *,"oldvfrac(im_probe) ",oldvfrac(im_probe)
             print *,"newvfrac(im_probe) ",newvfrac(im_probe)
             print *,"temp_mix_new(im_probe) ",temp_mix_new(im_probe)
            endif

           enddo !iprobe=1,2

            ! volume fractions updated on the 2nd sweep using
            ! deltaVOF which has "dF"
            ! centroids updated here.
           do udir=1,SDIM
            snew(D_DECL(i,j,k),vcompdst_snew+udir)= &
             new_centroid(im_dest,udir)-cengrid(udir)
            snew(D_DECL(i,j,k),vcompsrc_snew+udir)= &
             new_centroid(im_source,udir)-cengrid(udir)
           enddo ! udir
           if (ngeom_raw.eq.SDIM+1) then
            ! do nothing
           else
            print *,"ngeom_raw invalid in convert material"
            print *,"ngeom_raw= ",ngeom_raw
            stop
           endif

           if ((local_freezing_model.lt.0).or. &
               (local_freezing_model.gt.7)) then
            print *,"local_freezing_model invalid 2"
            stop
           endif
           if ((distribute_from_targ.lt.0).or. &
               (distribute_from_targ.gt.1)) then
            print *,"distribute_from_targ invalid"
            stop
           endif

           oldvfrac(im_dest)=recon(D_DECL(i,j,k),(im_dest-1)*ngeom_recon+1)
           oldvfrac(im_source)=recon(D_DECL(i,j,k),(im_source-1)*ngeom_recon+1)
           newvfrac(im_dest)=oldvfrac(im_dest)+dF
           newvfrac(im_source)=oldvfrac(im_source)-dF

           do iprobe=1,2  ! source,dest

            if (iprobe.eq.1) then ! source
             im_probe=im_source
            else if (iprobe.eq.2) then ! dest
             im_probe=im_dest
            else
             print *,"iprobe invalid"
             stop
            endif

            base_index=num_materials_vel*(SDIM+1)
            dencomp_probe=(im_probe-1)*num_state_material+1
            tcomp_probe=dencomp_probe+1

             ! iprobe==1: source
             ! iprobe==2: dest
             ! (rho F)^new - (rho F)^old
            den_dF(iprobe)=delta_mass_local(iprobe)

            if (newvfrac(im_probe).gt.one+VOFTOL) then
             print *,"newvfrac(im_probe) overflow"
             stop
            else if (newvfrac(im_probe).gt.one) then
             newvfrac(im_probe)=one
            else if (newvfrac(im_probe).lt.-VOFTOL) then
             print *,"newvfrac(im_probe) underflow"
             stop
            else if (newvfrac(im_probe).lt.zero) then
             newvfrac(im_probe)=zero
            else if ((newvfrac(im_probe).ge.zero).and. &
                     (newvfrac(im_probe).le.one)) then
             ! do nothing
            else
             print *,"newvfrac(im_probe) invalid"
             stop
            endif

#if (STANDALONE==0)
            thermal_k(iprobe)=get_user_heatviscconst(im_probe)
#elif (STANDALONE==1)
            thermal_k(iprobe)=fort_heatviscconst(im_probe)
#else
            print *,"bust compiling convertmaterial"
            stop
#endif

           enddo ! iprobe=1,2


           jump_strength=zero

           !dF=mdot/rho_source
           !dF_expand_dest=(rho_source/rho_dest - 1)*dF
           !dM=rho_dest * dF_dest + rho_source * dF_source=
           !rho_dest (mdot/rho_source + (rho_source/rho_dest - 1)(dF))+
           !rho_source * (-dF) =
           !((rho_dest/rho_source)mdot+(rho_source-rho_dest)*dF)-rho_source*dF=
           !(rho_dest/rho_source)mdot-rho_dest*dF=0
           if (distribute_from_targ.eq.0) then ! default
            ! distribute div u source to the cells in which F_dest>1/2  

            if (dF.le.EBVOFTOL) then
             denratio_factor=zero
            else if (dF.ge.EBVOFTOL) then
             if ((den_dF(2).gt.zero).and. &
                 (den_dF(1).lt.zero)) then
              denratio_factor=-den_dF(1)/den_dF(2)-one ! den_src/den_dst-1
             else if ((den_dF(2).eq.zero).or. &
                      (den_dF(1).eq.zero)) then
              ! do nothing
             else
              print *,"den_dF invalid 1"
              print *,"dF= ",dF
              print *,"den_dF(1) = ",den_dF(1)
              print *,"den_dF(2) = ",den_dF(2)
              print *,"EBVOFTOL = ",EBVOFTOL
              print *,"dt = ",dt
              stop
             endif
            else
             print *,"dF is corrupt"
             stop
            endif


           !dF=mdot/rho_dest
           !dF_expand_source=(1-rho_dest/rho_source)*dF
           !dM=rho_dest * dF_dest + rho_source * dF_source=
           !rho_dest(mdot/rho_dest)+ 
           !rho_source (-mdot/rho_dest+(1-rho_dest/rho_source)(dF))=
           !mdot-(rho_source/rho_dest)mdot+(rho_source-rho_dest)dF=
           !mdot-(rho_source/rho_dest)mdot+(rho_source/rho_dest)mdot-mdot=0
           else if (distribute_from_targ.eq.1) then
            ! distribute div u source to the cells in which F_dest<1/2  

            if (dF.le.EBVOFTOL) then
             denratio_factor=zero
            else if (dF.ge.EBVOFTOL) then
             if ((den_dF(2).gt.zero).and. &
                 (den_dF(1).lt.zero)) then
              denratio_factor=one+den_dF(2)/den_dF(1) ! 1-den_dst/den_src
             else if ((den_dF(2).eq.zero).or. &
                      (den_dF(1).eq.zero)) then
              ! do nothing
             else
              print *,"den_dF invalid 2"
              print *,"dF= ",dF
              print *,"den_dF(1) = ",den_dF(1)
              print *,"den_dF(2) = ",den_dF(2)
              print *,"EBVOFTOL = ",EBVOFTOL
              print *,"dt = ",dt
              stop
             endif
            else
             print *,"dF is corrupt"
             stop
            endif

           else 
            print *,"distribute_from_targ invalid"
            stop
           endif

           if (abs(denratio_factor).le.VOFTOL) then
            denratio_factor=zero
           endif

           jump_strength=denratio_factor/dt 

            !for distribute_from_targ==0:
            !initially: jump_strength=(den_src/den_dst - 1)/dt
            !ultimately jump_strength has units of cm^{3}/s^{2}
            !source term is jump_strength
            !dF has no units.  dF=F^new_dest - F^old_dest
            !note: vol div u/dt = cm^3 (cm/s) (1/cm) (1/s)=cm^3 / s^2
            !for boiling: jump_strength>0
           jump_strength=jump_strength*dF*volgrid/dt

           JUMPFAB(D_DECL(i,j,k),iten+ireverse*nten)=jump_strength
         
           if (dF.le.EBVOFTOL) then
            print *,"expecting dF>EBVOFTOL"
            stop
           else if (dF.ge.zero) then

             ! find mass weighted average of cv

            cvtotal=zero
            wttotal=zero
            do im_weight=1,nmat
             vofcomp_recon=(im_weight-1)*ngeom_recon+1
             Ftemp=recon(D_DECL(i,j,k),vofcomp_recon)*fort_denconst(im_weight)
#if (STANDALONE==0)
             local_cv_or_cp=get_user_stiffCP(im_weight)
#elif (STANDALONE==1)
             local_cv_or_cp=fort_stiffCP(im_weight)
#else
             print *,"bust compiling convertmaterial"
             stop
#endif
             cvtotal=cvtotal+Ftemp*local_cv_or_cp
             wttotal=wttotal+Ftemp
            enddo
            if (wttotal.le.zero) then
             print *,"wttotal invalid"
             stop
            endif
            cvtotal=cvtotal/wttotal
            if (cvtotal.le.zero) then
             print *,"cvtotal invalid"
             stop
            endif

            ! divide and conquer temperature equation (TSAT Dirichlet BC)
            if ((local_freezing_model.eq.0).or. &
                (local_freezing_model.eq.5).or. & ! Stefan evap/cond.
                (local_freezing_model.eq.6)) then ! Palmore/Desjardins 

              !F dt=Fn (tnp1-t) + Fnp1 (t-tn)
              !dt/2 = Fn tnp1 - Fnp1 tn + t(Fnp1-Fn)
              !t=(dt/2 - Fn tnp1 + Fnp1 tn)/(Fnp1-Fn)=
              !  (dt/2 - Fn (tn+dt) + Fnp1 tn)/dF=
              !  (dt/2 + dF tn - Fn dt)/dF
              !t-tn=dt(1/2 - Fn)/dF
              !tnp1-t=-dt(1/2-Fn)/dF+dt=
              !dt(dF-1/2+Fn)/dF=dt(Fnp1-1/2)/dF

             vofcomp_recon=(im_dest-1)*ngeom_recon+1
             do udir=1,SDIM
              xPOINT_supermesh(udir)=mofdata_new(vofcomp_recon+udir)+ &
                     cengrid(udir)
             enddo
             do udir=1,SDIM
              xPOINT_GFM(udir)=u_xsten_updatecell(0,udir)
             enddo

             if (supermesh_flag.eq.1) then
              do udir=1,SDIM
               xstar(udir)=xPOINT_supermesh(udir)
              enddo
             else if (supermesh_flag.eq.0) then
              do udir=1,SDIM
               xstar(udir)=xPOINT_GFM(udir)
              enddo
             else
              print *,"supermesh_flag invalid"
              stop
             endif

             data_out%data_interp=>cell_data_interp

             data_in%level=level
             data_in%finest_level=finest_level
             data_in%bfact=bfact
             data_in%nmat=nmat
             data_in%im_PLS=0 !0=> do not weight using LS
             data_in%dx=>dx
             data_in%xlo=>xlo
             data_in%fablo=>fablo
             data_in%fabhi=>fabhi
             data_in%ngrowfab=normal_probe_size+3

             data_in%state=>LSold
             data_in%LS=>LSold

             data_in%ncomp=nmat*(1+SDIM)
             data_in%scomp=1

             data_in%xtarget=>xstar
             data_in%interp_foot_flag=0 !=1 if interp xfoot from xdisp data

             call interp_from_grid_util(data_in,data_out)

             LS_dest_old=cell_data_interp(im_dest)

             tessellate=3
             call multi_get_volumePOINT( &
              tessellate, &
              bfact,dx, &
              u_xsten_updatecell,nhalf, &  ! absolute coordinate system
              mofdata, &
              xstar, & ! absolute coordinate system
              im_old_crit,nmat,SDIM)

             tessellate=3
             call multi_get_volumePOINT( &
               tessellate, &
               bfact,dx, &
               u_xsten_updatecell,nhalf, &  ! absolute coordinate system
               mofdata_new, &
               xstar, & ! absolute coordinate system
               im_new_crit,nmat,SDIM)

             if ((newvfrac(im_dest).gt.zero).and. &
                 (newvfrac(im_dest).le.one+EBVOFTOL)) then

               ! im_new_crit "owns" xstar
              if ((im_new_crit.eq.im_dest).or. &
                  (newvfrac(im_dest).ge.half)) then
               ! determine slope and intercept of the interface separating
               ! the source material from the destination.

               do iprobe=1,2

                if (iprobe.eq.1) then ! source
                 im_probe=im_source
                else if (iprobe.eq.2) then ! dest
                 im_probe=im_dest
                else
                 print *,"iprobe invalid"
                 stop
                endif
                vofcomp_recon=(im_probe-1)*ngeom_recon+1

                order_probe(iprobe)=NINT(mofdata_new(vofcomp_recon+SDIM+1))
                do udir=1,SDIM 
                 nslope_probe(udir,iprobe)= &
                   recon(D_DECL(i,j,k),vofcomp_recon+SDIM+1+udir) !slope
                enddo
                intercept_probe(iprobe)= &
                  mofdata_new(vofcomp_recon+2*SDIM+2)

               enddo ! iprobe=1,2

               if (order_probe(1).eq.0) then
                do udir=1,SDIM 
                 nslope_dest(udir)=nslope_probe(udir,2)
                enddo
                intercept_dest=intercept_probe(2)
               else if (order_probe(2).eq.0) then
                do udir=1,SDIM 
                 nslope_dest(udir)=-nslope_probe(udir,1)
                enddo
                intercept_dest=-intercept_probe(1)
               else if (order_probe(1).lt.order_probe(2)) then
                do udir=1,SDIM 
                 nslope_dest(udir)=-nslope_probe(udir,1)
                enddo
                intercept_dest=-intercept_probe(1)
               else if (order_probe(2).le.order_probe(1)) then
                do udir=1,SDIM 
                 nslope_dest(udir)=nslope_probe(udir,2)
                enddo
                intercept_dest=intercept_probe(2)
               else
                print *,"order_probe bust"
                stop
               endif

               LS_dest_new=intercept_dest 
               do udir=1,SDIM 
                LS_dest_new=LS_dest_new+nslope_dest(udir)* &
                       (xstar(udir)-u_xsten_updatecell(0,udir))
               enddo

              else if ((im_new_crit.ge.1).and. &
                       (im_new_crit.le.nmat).and. &
                       (im_new_crit.ne.im_dest).and. &
                       (newvfrac(im_dest).le.half)) then
               ! nothing is swept
              else
               print *,"im_new_crit invalid"
               stop
              endif

             else
              print *,"expecting newvfrac(im_dest)>0 since dF>0"
              stop
             endif

             if ((dF.gt.zero).and.(dF.le.one+VOFTOL)) then
              ! do nothing
             else
              print *,"dF invalid"
              stop
             endif

             if (((im_new_crit.eq.im_dest).or. &
                  (newvfrac(im_dest).ge.half)).and. &
                 ((im_old_crit.eq.im_source).or. &
                  (oldvfrac(im_source).ge.half))) then

                 ! if order_probe(1) or (2) == 0 => no slope found
                 ! in the reconstruction of im_source or im_dest
                 ! materials.
              if (newvfrac(im_dest).ge.one-VOFTOL) then
               SWEPTFACTOR=one
               !(1) order in slope recon of im_source
               !(2) order in slope recon of im_dest
              else if (oldvfrac(im_source).ge.one-VOFTOL) then
               if (supermesh_flag.eq.0) then
                SWEPTFACTOR=zero
               else if (supermesh_flag.eq.1) then

                if (((order_probe(2).ge.1).and. &
                     (order_probe(2).le.nmat)).or. &
                    ((order_probe(1).ge.1).and. &
                     (order_probe(1).le.nmat))) then

                 if ((LS_dest_old.eq.zero).and.(LS_dest_new.gt.zero)) then
                  SWEPTFACTOR=one
                 else if ((LS_dest_old.lt.zero).and.(LS_dest_new.eq.zero)) then
                  SWEPTFACTOR=LSTOL
                 else if ((LS_dest_old.ge.zero).or. &
                          (LS_dest_new.le.zero)) then
                  SWEPTFACTOR=one
                 else if (LS_dest_new-LS_dest_old.gt.zero) then
                  SWEPTFACTOR=-LS_dest_old/ &
                       (LS_dest_new-LS_dest_old)
                 else
                  print *,"LS_dest_new or LS_dest_old invalid"
                  stop
                 endif

                else
                 print *,"im_dest material disappeared at tnp1"
                 stop
                endif

               else
                print *,"supermesh_flag invalid"
                stop
               endif

               if (SWEPTFACTOR.le.LSTOL) then
                SWEPTFACTOR=LSTOL
               endif

              else if ((order_probe(1).eq.0).and. &
                       (order_probe(2).eq.0)) then 
               SWEPTFACTOR=one ! default
              else if ((order_probe(1).gt.0).or. &
                       (order_probe(2).gt.0)) then

               if ((LS_dest_old.eq.zero).and.(LS_dest_new.gt.zero)) then
                SWEPTFACTOR=one
               else if ((LS_dest_old.lt.zero).and.(LS_dest_new.eq.zero)) then
                SWEPTFACTOR=LSTOL
               else if ((LS_dest_old.ge.zero).or. &
                        (LS_dest_new.le.zero)) then
                SWEPTFACTOR=one
               else if (LS_dest_new-LS_dest_old.gt.zero) then
                SWEPTFACTOR=-LS_dest_old/ &
                       (LS_dest_new-LS_dest_old)
               else
                print *,"LS_dest_new or LS_dest_old invalid"
                stop
               endif
               if (SWEPTFACTOR.le.LSTOL) then
                SWEPTFACTOR=LSTOL
               endif

              else
               print *,"order_probe bust"
               stop
              endif

              if ((SWEPTFACTOR.ge.LSTOL).and.(SWEPTFACTOR.le.one)) then
               ! do nothing
              else
               print *,"SWEPTFACTOR invalid: ",SWEPTFACTOR
               print *,"dF=",dF
               print *,"im_dest=",im_dest
               print *,"newvfrac(im_dest) ",newvfrac(im_dest)
               print *,"oldvfrac(im_dest) ",oldvfrac(im_dest)
               print *,"LSTOL ",LSTOL
               stop
              endif
              swept(D_DECL(i,j,k),im_dest)=SWEPTFACTOR

             else if (((im_new_crit.ne.im_dest).and. &
                       (newvfrac(im_dest).le.half)).or. &
                      ((im_old_crit.ne.im_source).and. &
                       (oldvfrac(im_source).le.half))) then
              ! do nothing
             else
              print *,"im_new_crit, im_old_crit"
              stop
             endif

! single (continuum method) temperature equation for both phases.
! source term at the interface.
! latent_heat<0 condensation or solidification
! latent_heat>0 boiling or melting
! units of specific heat: J/(kg K)
! units of latent heat: J/kg
            else if (local_freezing_model.eq.1) then ! source term

             if (dF.gt.zero) then
              energy_source=-LL*dF
              do im_weight=1,nmat
               tcomp_wt=num_materials_vel*(SDIM+1)+ &
                (im_weight-1)*num_state_material+2
               snew(D_DECL(i,j,k),tcomp_wt)= &
                snew(D_DECL(i,j,k),tcomp_wt)+energy_source/cvtotal
              enddo
             else if (dF.eq.zero) then
              ! do nothing
             else
              print *,"dF invalid"
              stop
             endif

! "single (continuum method) temperature equation for both phases.
! source term at the interface.  Hydrates.
! rho c T^new - rho c T^old = rho (dt A LL/V) dS/dt = LL * dF
! c T^new - c T^old = (dt A LL/V) dS/dt = LL * dF 
            else if (local_freezing_model.eq.2) then ! hydrate

#if (STANDALONE==0)

             if (distribute_from_targ.ne.0) then
              print *,"distribute_from_targ invalid"
              stop
             endif
             if (dF.gt.zero) then
              if (num_species_var.ne.1) then
               print *,"num_species_var invalid"
               stop
              endif
   
              call Hydrate_energy_source_term(dF,dt, &
               thermal_k(1), &  ! source
               energy_source,LL)
              call Methane_usage(dF,dt, &
               fort_speciesviscconst(im_dest),amount_used)

              ccomp=num_materials_vel*(SDIM+1)+ &
               (im_dest-1)*num_state_material+3
              methaneC_old=snew(D_DECL(i,j,k),ccomp)*oldvfrac(im_dest)
              if (methaneC_old.ge.amount_used) then
               methaneC_old=methaneC_old-amount_used
              else
               methaneC_old=zero
              endif 
              if (newvfrac(im_dest).eq.zero) then
               print *,"newvfrac(im_dest) invalid"
               stop
              endif
              methaneC_new=methaneC_old/newvfrac(im_dest)
              snew(D_DECL(i,j,k),ccomp)=methaneC_new
     
              do im_weight=1,nmat
               tcomp_wt=num_materials_vel*(SDIM+1)+ &
                (im_weight-1)*num_state_material+2
               snew(D_DECL(i,j,k),tcomp_wt)= &
                snew(D_DECL(i,j,k),tcomp_wt)+energy_source/cvtotal
              enddo

             else if (dF.eq.zero) then
              ! do nothing
             else
              print *,"dF invalid"
              stop
             endif
#elif (STANDALONE==1)
             print *,"local_freezing_model cannot be 2 (convertmaterial)"
             stop
#else
             print *,"bust compiling convertmaterial"
             stop
#endif

            else if (local_freezing_model.eq.4) then ! Tanasawa or Schrage
              ! if LL>0 => evaporation => delete energy 
              ! if LL<0 => condensation => add energy 
              ! latent_heat: erg/g
              ! cv: erg/(g Kelvin)
              ! 
             if (dF.gt.zero) then
              energy_source=-LL*dF
              do im_weight=1,nmat
               tcomp_wt=num_materials_vel*(SDIM+1)+ &
                (im_weight-1)*num_state_material+2
               snew(D_DECL(i,j,k),tcomp_wt)= &
                snew(D_DECL(i,j,k),tcomp_wt)+energy_source/cvtotal
              enddo
             else if (dF.eq.zero) then
              ! do nothing
             else
              print *,"dF invalid"
              stop
             endif
            else if (local_freezing_model.eq.7) then ! Cavitation
             print *,"FIX ME"
             stop
            else
             print *,"local_freezing_model invalid in convertmaterial(2)"
             print *,"local_freezing_model= ",local_freezing_model
             print *,"iten,ireverse,nten ",iten,ireverse,nten
             stop
            endif

            do iprobe=1,2

             if (iprobe.eq.1) then ! source
              im_probe=im_source
             else if (iprobe.eq.2) then ! dest
              im_probe=im_dest
             else
              print *,"iprobe invalid"
              stop
             endif
             base_index=num_materials_vel*(SDIM+1)
             snew(D_DECL(i,j,k), &
                 base_index+nmat*num_state_material+ &
                 (im_probe-1)*ngeom_raw+1)=newvfrac(im_probe)

            enddo ! iprobe=1,2

            base_index=num_materials_vel*(SDIM+1)+ &
             nmat*num_state_material

            do u_im=1,nmat*ngeom_recon
             mofdata(u_im)=zero
            enddo

            do u_im=1,nmat
             vofcomp_recon=(u_im-1)*ngeom_recon+1
             vofcomp_raw=(u_im-1)*ngeom_raw+1
             do dir=0,SDIM
              mofdata(vofcomp_recon+dir)= &
                snew(D_DECL(i,j,k),base_index+vofcomp_raw+dir)
             enddo
            enddo ! u_im=1..nmat

            ! sum of F_fluid=1
            ! sum of F_rigid<=1
            tessellate=0
            call make_vfrac_sum_ok_base( &
              cmofsten, &
              u_xsten_updatecell,nhalf,nhalf_box, &
              bfact,dx, &
              tessellate,mofdata,nmat,SDIM,106)

            do u_im=1,nmat
             vofcomp_recon=(u_im-1)*ngeom_recon+1
             vofcomp_raw=(u_im-1)*ngeom_raw+1
             do dir=0,SDIM
              snew(D_DECL(i,j,k),base_index+vofcomp_raw+dir)= &
                 mofdata(vofcomp_recon+dir)
             enddo
            enddo ! u_im=1..nmat

            delta_mass(im_source)=delta_mass(im_source)+ &
             volgrid*(newvfrac(im_source)-oldvfrac(im_source))
            delta_mass(im_dest+nmat)=delta_mass(im_dest+nmat)+ &
             volgrid*(newvfrac(im_dest)-oldvfrac(im_dest))

           else
            print *,"dF bust"
            stop
           endif

          else if ((dF.ge.zero).and.(dF.le.EBVOFTOL)) then
           ! do nothing
          else
           print *,"dF became corrupt2 dF=",dF
           print *,"dF invalid"
           stop
          endif 

         else if (do_unsplit_advection.eq.0) then
          ! do nothing
         else
          print *,"do_unsplit_advection invalid"
          stop
         endif

         ! set vapor mass fraction to Ygamma where possible 
         ! in the corresponding liquid material.
         do im=1,nmat-1
          do im_opp=im+1,nmat

           call get_iten(im,im_opp,iten,nmat)

           do ireverse=0,1
            if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
             print *,"im or im_opp bust 10"
             stop
            endif

            LL=latent_heat(iten+ireverse*nten)

            if ((is_rigid(nmat,im).eq.1).or. &
                (is_rigid(nmat,im_opp).eq.1)) then
             ! do nothing
            else if ((is_rigid(nmat,im).eq.0).and. &
                     (is_rigid(nmat,im_opp).eq.0)) then

             if (LL.ne.zero) then
              if (ireverse.eq.0) then
               im_source=im
               im_dest=im_opp
              else if (ireverse.eq.1) then
               im_source=im_opp
               im_dest=im
              else
               print *,"ireverse invalid"
               stop
              endif
              local_freezing_model=freezing_model(iten+ireverse*nten)
              mass_frac_id=mass_fraction_id(iten+ireverse*nten)

              if (local_freezing_model.eq.0) then ! standard Stefan model
               ! do nothing
              else if (local_freezing_model.eq.1) then ! source term
               ! do nothing
              else if (local_freezing_model.eq.2) then !hydrate
               ! do nothing
              else if (local_freezing_model.eq.3) then !wildfire
               ! do nothing
              else if ((local_freezing_model.eq.4).or. & ! Tannasawa or Schrage
                       (local_freezing_model.eq.5).or. & ! Stefan model
                       (local_freezing_model.eq.6).or. & ! Palmore/Desjardins
                       (local_freezing_model.eq.7)) then ! Cavitation

               if ((mass_frac_id.ge.1).and. &
                   (mass_frac_id.le.num_species_var)) then

                if (LL.gt.zero) then ! evaporation
                 im_condensed=im_source
                else if (LL.lt.zero) then ! condensation
                 im_condensed=im_dest
                else
                 print *,"LL invalid"
                 stop
                endif
                       
                speccomp_mod=num_materials_vel*(SDIM+1)+ &
                 (im_condensed-1)*num_state_material+num_state_base+ &
                 mass_frac_id
          
                default_comp=(mass_frac_id-1)*nmat+im_condensed 
                Ygamma_default=fort_speciesconst(default_comp)

                Tsat_flag=NINT(TgammaFAB(D_DECL(i,j,k),iten))
                if (ireverse.eq.0) then
                 ! do nothing
                else if (ireverse.eq.1) then
                 Tsat_flag=-Tsat_flag
                else
                 print *,"ireverse invalid"
                 stop
                endif
            
                if ((Tsat_flag.eq.1).or.(Tsat_flag.eq.2)) then
                 Ygamma_default=TgammaFAB(D_DECL(i,j,k), &
                   nten+(iten-1)*ncomp_per_tsat+2)
                else if ((Tsat_flag.eq.-1).or. &
                         (Tsat_flag.eq.-2)) then
                 ! do nothing
                else if (Tsat_flag.eq.0) then
                 ! do nothing
                else
                 print *,"Tsat_flag invalid"
                 stop
                endif

                snew(D_DECL(i,j,k),speccomp_mod)=Ygamma_default

                 ! for compressible flows, due to round off error, the
                 ! mass fraction might go out of bounds [0,1] ?
                do iprobe=1,2

                 if (iprobe.eq.1) then ! source
                  im_probe=im_source
                 else if (iprobe.eq.2) then ! dest
                  im_probe=im_dest
                 else
                  print *,"iprobe invalid"
                  stop
                 endif
                 speccomp_mod=num_materials_vel*(SDIM+1)+ &
                  (im_probe-1)*num_state_material+num_state_base+ &
                  mass_frac_id
                 mass_frac_limit=snew(D_DECL(i,j,k),speccomp_mod)
                 if ((mass_frac_limit.ge.-VOFTOL).and. &
                     (mass_frac_limit.le.zero)) then
                  mass_frac_limit=zero
                 else if ((mass_frac_limit.ge.zero).and. &
                          (mass_frac_limit.le.one)) then
                  ! do nothing
                 else if ((mass_frac_limit.ge.one).and. &
                          (mass_frac_limit.le.one+VOFTOL)) then
                  mass_frac_limit=one
                 else
                  print *,"mass_frac_limit invalid: ",mass_frac_limit
                  print *,"is advection_order=1? "
                  stop
                 endif
                 snew(D_DECL(i,j,k),speccomp_mod)=mass_frac_limit
                enddo ! iprobe=1,2

               else
                print *,"mass_frac_id invalid"
                stop
               endif

              else
               print *,"local_freezing_model invalid 1"
               stop
              endif

             else if (LL.eq.zero) then
              ! do nothing
             else
              print *,"LL invalid"
              stop
             endif
            else
             print *,"is_rigid invalid"
             stop
            endif
           enddo ! ireverse=0,1
          enddo ! im_opp
         enddo ! im
        else if (is_rigid(nmat,im_primary).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(nmat,im_primary) invalid"
         stop
        endif

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine FORT_CONVERTMATERIAL


        ! ngrow corresponds to normal_probe_size+3
      subroutine FORT_EXTEND_BURNING_VEL( &
        velflag, &
        level, &
        finest_level, &
        xlo,dx, &
        nmat, &
        nten, &
        nburning, &
        ngrow, &
        latent_heat, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        vel,DIMS(vel), &
        LS,DIMS(LS))  ! old level set function before phase change
      use probcommon_module
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: velflag
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level 
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nburning
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact

      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      ! first nten components are the status.
      REAL_T, intent(inout) :: vel(DIMV(vel),nburning)
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(SDIM+1))

      INTEGER_T im,im_opp
      INTEGER_T iten,ireverse,sign_local
      INTEGER_T im_dest,im_source
      INTEGER_T i,j,k,dir
      INTEGER_T i_sp,j_sp,k_sp
      INTEGER_T sten_lo(3)
      INTEGER_T sten_hi(3)
      INTEGER_T indexcp(SDIM)
      INTEGER_T tag_local
      REAL_T rtag_local

      REAL_T xijk(-3:3,SDIM)
      REAL_T xcp(SDIM)
      REAL_T xsp(-3:3,SDIM)
      REAL_T nrm(SDIM)
      REAL_T weight,total_weight
      REAL_T vel_sum(SDIM)
      REAL_T dxmax,dxmaxLS,eps,extensionwidth,LL
      REAL_T ls_local(nmat)
      INTEGER_T im_local
      INTEGER_T nten_test
      INTEGER_T scomp,nhalf
      INTEGER_T ncomp_per

      nhalf=3

      if (velflag.eq.1) then
       ncomp_per=SDIM
      else if (velflag.eq.0) then
       ncomp_per=2 ! interface temperature, mass fraction
      else
       print *,"velflag invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid118"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in extend_burning_vel"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (ngrow.ne.4) then
       print *,"expecting ngrow==4 in FORT_EXTEND_BURNING_VEL"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"expecting ngrow_make_distance==3 in FORT_EXTEND_BURNING_VEL"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid velextend nten, nten_test ",nten,nten_test
       stop
      endif
      if (nburning.eq.nten*(ncomp_per+1)) then
       ! do nothing
      else
       print *,"nburning invalid"
       stop
      endif

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)
 
       ! Guard against division zero in the weight calculation
      eps=dxmaxLS*1.0E-4

      extensionwidth=dxmaxLS*ngrow_make_distance 

      call checkbound(fablo,fabhi,DIMS(vel),ngrow_make_distance,-1,1250)
      call checkbound(fablo,fabhi,DIMS(LS),ngrow,-1,1250)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do im=1,nmat-1
       do im_opp=im+1,nmat
        do ireverse=0,1
         if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
          print *,"im or im_opp bust"
          stop
         endif
         if (ireverse.eq.0) then
          im_source=im
          im_dest=im_opp
          sign_local=1
         else if (ireverse.eq.1) then
          im_source=im_opp
          im_dest=im
          sign_local=-1
         else
          print *,"ireverse invalid"
          stop
         endif
  
         call get_iten(im,im_opp,iten,nmat)
         LL=latent_heat(iten+ireverse*nten)

         if (LL.ne.zero) then
  
          ! iterate over domain and check tag
          do i=growlo(1),growhi(1)
          do j=growlo(2),growhi(2)
          do k=growlo(3),growhi(3)
           do im_local=1,nmat
            ls_local(im_local)=LS(D_DECL(i,j,k),im_local)
           enddo
           call get_primary_material(ls_local,nmat,im_local)
           if (is_rigid(nmat,im_local).eq.0) then

            if ((im_dest.ge.1).and.(im_dest.le.nmat)) then
             rtag_local=vel(D_DECL(i,j,k),iten)
             tag_local=NINT(rtag_local) 
              ! tag_local==0 if phase change rate not yet init.
             if ((tag_local.eq.0).and.(rtag_local.eq.zero)) then
              ! if close enough to both materials
              if (max(abs(ls_local(im_source)),abs(ls_local(im_dest))).le. &
                  extensionwidth) then
               ! Project a point in normal direction of 
               ! the level set function with least magnitude
               call gridsten_level(xijk,i,j,k,level,nhalf)

               if ((ls_local(im_source).ge.zero).and. &
                   (ls_local(im_dest).le.zero)) then
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),nmat+(im_dest-1)*SDIM+dir)
                 xcp(dir)=xijk(0,dir)-ls_local(im_dest)*nrm(dir)
                enddo ! dir
               else if ((ls_local(im_dest).ge.zero).and. &
                        (ls_local(im_source).le.zero)) then
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),nmat+(im_source-1)*SDIM+dir)
                 xcp(dir)=xijk(0,dir)-ls_local(im_source)*nrm(dir)
                enddo ! dir
               else if (abs(ls_local(im_source)).le. &
                        abs(ls_local(im_dest))) then

                ! closest point to "dest" material is farther
                ! from the "source" material closest point, but
                ! most importantly the "dest" material closest 
                ! point is at the triple point.
                ! ALSO: this procedure is more stable since if 
                ! ls_source perturbed slightly from positive to negative,
                ! then the algorithm to find the closest point will
                ! not change.
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),nmat+(im_dest-1)*SDIM+dir)
                 xcp(dir)=xijk(0,dir)-ls_local(im_dest)*nrm(dir)
                enddo ! dir
               else if (abs(ls_local(im_dest)).le. &
                        abs(ls_local(im_source))) then

                ! closest point to "source" material is farther
                ! from the "dest" material closest point, but
                ! most importantly the "source" material closest 
                ! point is at the triple point.
                ! ALSO: this procedure is more stable since if 
                ! ls_dest perturbed slightly from positive to negative,
                ! then the algorithm to find the closest point will
                ! not change.
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),nmat+(im_source-1)*SDIM+dir)
                 xcp(dir)=xijk(0,dir)-ls_local(im_source)*nrm(dir)
                enddo ! dir
               else
                print *,"ls_local bust"
                stop
               endif

               sten_lo(3)=0
               sten_hi(3)=0
               call containing_cell(bfact,dx,xlo,fablo,xcp,indexcp)
               do dir=1,SDIM
                sten_lo(dir)=indexcp(dir)-1
                if (sten_lo(dir).lt.fablo(dir)-ngrow_make_distance) then
                 sten_lo(dir)=fablo(dir)-ngrow_make_distance
                endif
                sten_hi(dir)=indexcp(dir)+1
                if (sten_hi(dir).gt.fabhi(dir)+ngrow_make_distance) then
                 sten_hi(dir)=fabhi(dir)+ngrow_make_distance
                endif
               enddo ! dir=1..sdim

               do dir=1,ncomp_per
                vel_sum(dir)=zero
               enddo

               total_weight = zero

               ! gather support point/weights
               do i_sp = sten_lo(1),sten_hi(1)
               do j_sp = sten_lo(2),sten_hi(2)
               do k_sp = sten_lo(3),sten_hi(3)
                rtag_local=vel(D_DECL(i_sp,j_sp,k_sp),iten) 
                tag_local=NINT(rtag_local) 

                if ((tag_local.eq.sign_local).and. &
                    (rtag_local.eq.sign_local)) then
                 call gridsten_level(xsp,i_sp,j_sp,k_sp,level,nhalf)

                 weight=zero
                 do dir=1,SDIM
                  weight = weight + (xsp(0,dir)-xcp(dir))**two
                 enddo
                 weight=sqrt(weight)
                 weight =one/((weight+eps)**four)
                 total_weight = total_weight+weight
                 do dir=1,ncomp_per
                  scomp=nten+(iten-1)*ncomp_per+dir
                  vel_sum(dir) = vel_sum(dir)+ &
                   weight*vel(D_DECL(i_sp,j_sp,k_sp),scomp)
                 enddo
                else if ((tag_local.eq.0).and.(rtag_local.eq.zero)) then
                 ! do nothing
                else if ((abs(tag_local).eq.2).and. &
                         (abs(rtag_local).eq.two)) then
                 ! do nothing (this is a newly extended rate)
                else if ((tag_local.eq.-sign_local).and. &
                         (rtag_local.eq.-sign_local)) then
                 ! do nothing (this is opposite sign)
                else
                 print *,"tag_local or rtag_local invalid1:", &
                    tag_local,rtag_local
                 stop
                endif
               enddo
               enddo
               enddo ! i_sp,j_sp,k_sp

               if (total_weight.gt.zero) then
                do dir=1,ncomp_per
                 scomp=nten+(iten-1)*ncomp_per+dir
                 vel(D_DECL(i,j,k),scomp)=vel_sum(dir)/total_weight
                enddo
                vel(D_DECL(i,j,k),iten)=two*sign_local
               else if (total_weight.eq.zero) then
                ! do nothing
               else
                print *,"total_weight invalid"
                stop
               endif
              endif ! min(|ls_w|,|ls_i|)<extension width
             else if ((tag_local.eq.1).and.(rtag_local.eq.one)) then
              ! do nothing
             else if ((tag_local.eq.-1).and.(rtag_local.eq.-one)) then
              ! do nothing
             else
              print *,"tag_local or rtag_local invalid2:", &
                 tag_local,rtag_local
              print *,"i,j,k ",i,j,k
              print *,"LL= ",LL
              print *,"iten=",iten
              print *,"im,im_opp,ireverse,im_source,im_dest ", &
               im,im_opp,ireverse,im_source,im_dest
              print *,"raw tag data ",vel(D_DECL(i,j,k),iten)
              print *,"level,finest_level ",level,finest_level
              print *,"nmat,nten,nburning,ngrow ",nmat,nten,nburning,ngrow
              stop
             endif
            else
             print *,"im_dest invalid"
             stop
            endif

           else if (is_rigid(nmat,im_local).eq.1) then
            ! do nothing
           else
            print *,"is_rigid(nmat,im_local) invalid"
            stop
           endif
          enddo ! k
          enddo ! j
          enddo ! i
         else if (LL.eq.zero) then
          ! do nothing
         else
          print *,"LL bust"
          stop
         endif 
        enddo ! ireverse=0..1
       enddo ! im_opp=im+1..nmat
      enddo ! im=1..nmat-1

      return 
      end subroutine FORT_EXTEND_BURNING_VEL

 
      ! vof,ref centroid,order,slope,intercept  x nmat
      ! LS_slopes_FD comes from FORT_FD_NORMAL (MOF_REDIST)
      ! FORT_FD_NORMAL calls find_cut_geom_slope_CLSVOF:
      ! finds grad phi/|grad phi| where grad=(d/dx,d/dy,d/dz) or
      ! grad=(d/dr,d/dz) or
      ! grad=(d/dr,d/dtheta,d/dz)
      ! for evaporation the following equations are needed:
      ! T_interface=f_{1}(Y_interface)  (Clasius Clapyron condition)
      ! mdot_T = mdot_Y
      ! T_interface=f_{2}(Y_interface)
      ! f_{1} is an increasing function
      ! f_{2} is a decreasing function.
      ! Look for intersection of f_{1} and f_{2}
      ! g(Y)=f_{1}(y)-f_{2}(y)
      ! Y0 given
      ! Y_{n+1} = Y_{n} - g(Y_n)/( (g(Y_{n})-g(Y_{n-1}))/(Y_{n}-Y_{n-1}))
      ! Palmore and Desjardins
      ! Secant method will be implemented.
      subroutine FORT_RATEMASSCHANGE( &
       tid, &
       nucleation_flag, &
       stefan_flag, &  ! do not update LSnew if stefan_flag==0
       level, &
       finest_level, &
       normal_probe_size, &
       ngrow_distance, &
       nstate, &
       nmat, &
       nten, &
       nburning, &
       ntsat, &
       nden, &
       custom_nucleation_model, &
       do_the_nucleate, &
       nucleate_pos, &
       nucleate_pos_size, &
       nucleation_temp, &
       nucleation_pressure, &
       nucleation_pmg, &
       nucleation_mach, &
       cavitation_pressure, &
       cavitation_vapor_density, &
       cavitation_tension, &
       microlayer_substrate, &
       microlayer_angle, &
       microlayer_size, &
       macrolayer_size, &
       max_contact_line_size, &
       R_Palmore_Desjardins, &
       latent_heat, &
       use_exact_temperature, &
       reaction_rate, &
       hardwire_Y_gamma, &
       hardwire_T_gamma, &
       accommodation_coefficient, &
       reference_pressure, &
       saturation_temp, &
       saturation_temp_curv, &
       saturation_temp_vel, &
       saturation_temp_min, &
       saturation_temp_max, &
       freezing_model, &
       Tanasawa_or_Schrage_or_Kassemi, &
       distribute_from_target, &
       mass_fraction_id, &
       constant_density_all_time, & ! 1..nmat
       material_type_evap, &
       molar_mass, &
       species_molar_mass, &
       use_supermesh, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx, &
       prev_time, &
       dt, &
       arraysize, &
       blob_array, &
       num_elements_blobclass, &
       color_count, &
       colorfab,DIMS(colorfab), &
       typefab,DIMS(typefab), &
       maskcov,DIMS(maskcov), &
       burnvel,DIMS(burnvel), &
       Tsatfab,DIMS(Tsatfab), &
       LS,DIMS(LS),  &
       LSnew,DIMS(LSnew), & 
       Snew,DIMS(Snew), & 
       LS_slopes_FD, &
       DIMS(LS_slopes_FD), & 
       EOS,DIMS(EOS), &
       recon,DIMS(recon), &
       pres,DIMS(pres), &
       pres_eos,DIMS(pres_eos), &
       curvfab,DIMS(curvfab) )
#if (STANDALONE==0)
      use probf90_module
#elif (STANDALONE==1)
      use probcommon_module
#endif
      use global_utility_module
      use MOF_routines_module
      use mass_transfer_module
#if (STANDALONE==0)
      use hydrateReactor_module
#endif

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: nucleation_flag
      INTEGER_T, intent(in) :: stefan_flag
      INTEGER_T, target, intent(in) :: level
      INTEGER_T, target, intent(in) :: finest_level
      INTEGER_T, intent(in) :: normal_probe_size
      REAL_T :: microscale_probe_size
      INTEGER_T, intent(in) :: ngrow_distance
      INTEGER_T, intent(in) :: nstate
      INTEGER_T, target, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nburning
      INTEGER_T, intent(in) :: ntsat
      INTEGER_T, intent(in) :: nden
      INTEGER_T, intent(in) :: custom_nucleation_model
      INTEGER_T, intent(in) :: do_the_nucleate
      INTEGER_T, intent(in) :: nucleate_pos_size
      REAL_T, target, intent(in) :: nucleate_pos(nucleate_pos_size)
      REAL_T, target, intent(in) :: nucleation_temp(2*nten)
      REAL_T, target, intent(in) :: nucleation_pressure(2*nten)
      REAL_T, target, intent(in) :: nucleation_pmg(2*nten)
      REAL_T, target, intent(in) :: nucleation_mach(2*nten)
      REAL_T, target, intent(in) :: cavitation_pressure(nmat)
      REAL_T, target, intent(in) :: cavitation_vapor_density(nmat)
      REAL_T, target, intent(in) :: cavitation_tension(nmat)
      INTEGER_T, intent(in) ::  microlayer_substrate(nmat)
      REAL_T, intent(in) :: microlayer_angle(nmat)
      REAL_T, intent(in) :: microlayer_size(nmat)
      REAL_T, intent(in) :: macrolayer_size(nmat)
      REAL_T, intent(in) :: max_contact_line_size(nmat)
      REAL_T, intent(in) :: R_Palmore_Desjardins
      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: use_exact_temperature(2*nten)
      REAL_T, intent(in) :: reaction_rate(2*nten)
      REAL_T :: K_f(0:1)
      REAL_T, intent(in) :: hardwire_Y_gamma(2*nten)
      REAL_T, intent(in) :: hardwire_T_gamma(2*nten)
      REAL_T, intent(in) :: accommodation_coefficient(2*nten)
      REAL_T, intent(in) :: reference_pressure(2*nten)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      REAL_T, intent(in) :: saturation_temp_curv(2*nten)
      REAL_T, intent(in) :: saturation_temp_vel(2*nten)
      REAL_T, intent(in) :: saturation_temp_min(2*nten)
      REAL_T, intent(in) :: saturation_temp_max(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: Tanasawa_or_Schrage_or_Kassemi(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      INTEGER_T, intent(in) :: mass_fraction_id(2*nten)
      INTEGER_T, intent(in), target :: material_type_evap(nmat)
      REAL_T, intent(in) :: molar_mass(nmat)
      REAL_T, intent(in) :: species_molar_mass(num_species_var+1)
      INTEGER_T, intent(in) :: constant_density_all_time(nmat)
      INTEGER_T, intent(in) :: use_supermesh
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, target, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, target, intent(in) :: bfact
      REAL_T, target, intent(in) :: xlo(SDIM)
      REAL_T, target, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: prev_time
      REAL_T :: cur_time
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: arraysize
      REAL_T, intent(in) :: blob_array(arraysize)
      INTEGER_T, intent(in) :: num_elements_blobclass
      INTEGER_T, intent(in) :: color_count
      INTEGER_T, intent(in) :: DIMDEC(colorfab)
      INTEGER_T, intent(in) :: DIMDEC(typefab)
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(burnvel)
      INTEGER_T, intent(in) :: DIMDEC(Tsatfab)
      INTEGER_T, intent(in) :: DIMDEC(LS) ! declare the x,y,z dimensions of LS
      INTEGER_T, intent(in) :: DIMDEC(LSnew)
      INTEGER_T, intent(in) :: DIMDEC(Snew)
      INTEGER_T, intent(in) :: DIMDEC(LS_slopes_FD)
      INTEGER_T, intent(in) :: DIMDEC(EOS)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(pres)
      INTEGER_T, intent(in) :: DIMDEC(pres_eos)
      INTEGER_T, intent(in) :: DIMDEC(curvfab)

      REAL_T, intent(in) :: typefab(DIMV(typefab))
      REAL_T, intent(in) :: colorfab(DIMV(colorfab))

      REAL_T, intent(in) :: maskcov(DIMV(maskcov)) 

        ! destination vel: first nten components are the status.
      REAL_T, intent(out) :: burnvel(DIMV(burnvel),nburning)
        ! TSAT : first nten components are the status.
      REAL_T, intent(out) :: Tsatfab(DIMV(Tsatfab),ntsat)
        ! LS1,LS2,..,LSn,normal1,normal2,...normal_n 
        ! normal points from negative to positive
        !DIMV(LS)=x,y,z  nmat=num. materials
      REAL_T, target, intent(in) :: LS(DIMV(LS),nmat*(SDIM+1)) 
      REAL_T, target, intent(inout) :: LSnew(DIMV(LSnew),nmat*(SDIM+1))
      REAL_T, target, intent(inout) :: Snew(DIMV(Snew),nstate)
      REAL_T, intent(in) :: LS_slopes_FD(DIMV(LS_slopes_FD),nmat*SDIM)
      REAL_T, target, intent(in) :: EOS(DIMV(EOS),nden)
       ! F,X,order,SL,I x nmat
      REAL_T, target, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon) 
      REAL_T, target, intent(in) :: pres(DIMV(pres)) 
      REAL_T, target, intent(in) :: pres_eos(DIMV(pres_eos)) 
      REAL_T, intent(in) :: curvfab(DIMV(curvfab),2*(nmat+nten)) 

      INTEGER_T, target :: i,j,k
      INTEGER_T dir,dir2
      INTEGER_T im,im_opp,ireverse,iten
      INTEGER_T imls
      INTEGER_T im_ambient
      INTEGER_T im_primary
      INTEGER_T, target :: imls_I
      INTEGER_T im_substrate_source
      INTEGER_T im_substrate_dest
      INTEGER_T, target :: im_source,im_dest
      INTEGER_T, target :: ngrow
      INTEGER_T nten_test
      INTEGER_T, target :: tcomp_source
      INTEGER_T, target :: tcomp_dest
      INTEGER_T, target :: Ycomp_source
      INTEGER_T, target :: Ycomp_dest
      REAL_T velmag_sum,local_velmag
      INTEGER_T burnflag
      REAL_T dxmin,dxmax
      REAL_T, target :: dxmaxLS
      REAL_T xsten(-3:3,SDIM)
      REAL_T, target :: xI(SDIM)
      REAL_T, target :: xsrc(SDIM)
      REAL_T, target :: xdst(SDIM)
      REAL_T, target :: xsrc_micro(SDIM)
      REAL_T, target :: xdst_micro(SDIM)
      REAL_T nrmCP(SDIM)  ! closest point normal
      REAL_T nrmFD(SDIM)  ! finite difference normal
      REAL_T nrmPROBE(SDIM)  ! must choose nrmCP is microstructure.
      REAL_T theta_nrmPROBE(SDIM)
      REAL_T, target :: LSINT(nmat*(SDIM+1))
      REAL_T LShere(nmat)
      REAL_T T_probe(2) ! iprobe=1 source; iprobe=2 dest.
      REAL_T T_probe_raw(2) ! iprobe=1 source; iprobe=2 dest.
      REAL_T Y_probe(2)
      REAL_T den_I_interp(2)
      REAL_T den_probe(2)
      REAL_T den_I_interp_SAT(2)
      REAL_T T_I_interp(2)
      REAL_T Y_I_interp(2)
      REAL_T pres_I_interp(2)
      REAL_T vfrac_I(2)
      REAL_T local_hardwire_T(0:1)
      REAL_T local_hardwire_Y(0:1)
      INTEGER_T hardwire_flag(0:1)  ! =0 do not hardwire  =1 hardwire
      REAL_T local_Tsat(0:1)
      REAL_T delta_Tsat
      REAL_T local_Tsat_base(0:1)
      REAL_T vel_phasechange(0:1)
      REAL_T, target :: LL(0:1)
      INTEGER_T valid_phase_change(0:1)
      REAL_T, target :: dxprobe_source
      REAL_T, target :: dxprobe_dest
      REAL_T dxprobe_target(2)
      REAL_T, target :: thermal_k(2) ! source,dest
      REAL_T LS_pos
      REAL_T C_w0
      INTEGER_T, target :: local_freezing_model
      INTEGER_T local_Tanasawa_or_Schrage_or_Kassemi
      INTEGER_T distribute_from_targ
      INTEGER_T at_interface
      INTEGER_T vofcomp_source,vofcomp_dest
      REAL_T Fsource,Fdest
      REAL_T LSSIGN,SIGNVEL
      INTEGER_T found_path
      INTEGER_T, target :: debugrate
      INTEGER_T nhalf
      REAL_T RR,mag
      INTEGER_T for_estdt
      INTEGER_T local_mask
      INTEGER_T microlayer_substrate_source
      INTEGER_T microlayer_substrate_dest
      REAL_T gradphi_substrate(SDIM)
      REAL_T newphi_substrate
      INTEGER_T, target :: dencomp_source,dencomp_dest
      INTEGER_T ispec
      REAL_T vapor_den
      REAL_T source_perim_factor,dest_perim_factor
      REAL_T contact_line_perim
      INTEGER_T icolor,base_type,ic,im1,im2
      REAL_T normal_probe_factor
      INTEGER_T iprobe_vapor
      REAL_T, target :: Y_TOLERANCE
       ! iten=1..nten  ireverse=0..1
      REAL_T temp_target_probe_history(2*nten,2)
      REAL_T dxprobe_target_history(2*nten,2)

      INTEGER_T use_tsatfab
      INTEGER_T ncomp_per_burning
      INTEGER_T ncomp_per_tsat

      REAL_T CURV_OUT_I

      REAL_T VEL_predict,VEL_correct
      REAL_T X_predict
      REAL_T Y_predict
      REAL_T TSAT_predict,TSAT_correct
      REAL_T TSAT_ERR,TSAT_INIT_ERR
      INTEGER_T TSAT_iter,TSAT_converge
      REAL_T :: Y_interface_min
      REAL_T :: TI_min
      REAL_T :: TI_max
      REAL_T FicksLawD(2)  ! iprobe=1 source iprobe=2 dest 
      REAL_T molar_mass_ambient
      REAL_T molar_mass_vapor
      INTEGER_T interp_valid_flag(2) ! iprobe=1 source iprobe=2 dest
      INTEGER_T interp_valid_flag_initial(2) ! iprobe=1 source iprobe=2 dest
      INTEGER_T interface_resolved
      type(probe_parm_type), target :: PROBE_PARMS
      type(TSAT_MASS_FRAC_parm_type) :: TSAT_Y_PARMS
      type(nucleation_parm_type_input) :: create_in
      type(nucleation_parm_type_inout) :: create_inout
      INTEGER_T iprobe,im_probe,microlayer_substrate_probe
      REAL_T x_gamma_a,x_gamma_b,x_gamma_c
      REAL_T Y_gamma_a,Y_gamma_b,Y_gamma_c
      REAL_T T_gamma_a,T_gamma_b,T_gamma_c
      REAL_T T_gamma_a_init,T_gamma_b_init
      REAL_T Y_gamma_a_init,Y_gamma_b_init
      REAL_T mdot_diff_a,mdot_diff_b,mdot_diff_c
      INTEGER_T fully_saturated
      REAL_T mdotT_debug
      REAL_T mdotY_top_debug,mdotY_bot_debug,mdotY_debug

#if (STANDALONE==1)
      REAL_T DTsrc,DTdst,velsrc,veldst,velsum
#endif

      nhalf=3

      if (prev_time.ge.zero) then
       cur_time=prev_time+dt
      else 
       print *,"prev_time cannot be negative, prev_time=",prev_time
       stop
      endif

      if (R_Palmore_Desjardins.eq.fort_R_Palmore_Desjardins) then
       ! do nothing
      else
       print *,"R_Palmore_Desjardins.ne.fort_R_Palmore_Desjardins"
       stop
      endif

      Y_TOLERANCE=0.01D0

      if ((stefan_flag.eq.0).or.(stefan_flag.eq.1)) then
       ! do nothing
      else
       print *,"stefan_flag invalid"
       stop
      endif

      if ((use_supermesh.eq.0).or. &
          (use_supermesh.eq.1)) then
       ! do nothing
      else
       print *,"use_supermesh invalid"
       stop
      endif

      if ((nucleation_flag.eq.0).or.(nucleation_flag.eq.1)) then
       ! do nothing
      else
       print *,"nucleation_flag invalid"
       stop
      endif

      debugrate=0
      if ((level.eq.finest_level).and.(1.eq.0)) then
       debugrate=1
      endif

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in ratemasschange"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (normal_probe_size.ne.1) then
       print *,"normal_probe_size invalid"
       stop
      endif

      microscale_probe_size=1.0D-2

      if (nucleate_pos_size.lt.4) then
       print *,"nucleate_pos_size invalid: ",nucleate_pos_size
       stop
      endif
      if (n_sites.gt.0) then
       if (nucleate_pos_size.ne.n_sites*4) then
        print *,"nucleate_pos_size invalid"
        stop
       endif
      endif
      if ((custom_nucleation_model.ne.0).and. &
          (custom_nucleation_model.ne.1)) then
       print *,"custom_nucleation_model invalid"
       stop
      endif

      ngrow=normal_probe_size+3

      if (ngrow_distance.ne.4) then
       print *,"expecting ngrow_distance==4 in RATEMASSCHANGE"
       stop
      endif

      if (ngrow.ne.ngrow_distance) then
       print *,"ngrow or ngrow_distance invalid in RATEMASSCHANGE"
       print *,"ngrow=",ngrow
       print *,"ngrow_distance=",ngrow_distance
       stop
      endif

      if (ngrow_make_distance.ne.3) then
       print *,"expecting ngrow_make_distance==3 in FORT_RATEMASSCHANGE"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid ratemass nten, nten_test ",nten,nten_test
       stop
      endif
      ncomp_per_burning=SDIM
      if (nburning.eq.nten*(ncomp_per_burning+1)) then
       ! do nothing
      else
       print *,"nburning invalid"
       stop
      endif
      ncomp_per_tsat=2 ! interface temperature, mass fraction
      if (ntsat.eq.nten*(ncomp_per_tsat+1)) then
       ! do nothing
      else
       print *,"ntsat invalid"
       stop
      endif
      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid in rate mass change"
       print *,"nden=",nden
       print *,"nmat=",nmat
       print *,"num_state_material=",num_state_material
       stop
      endif
      if (dt.le.zero) then
       print *,"dt invalid"
       stop
      endif
      
      if (nstate.eq.num_materials_vel*(SDIM+1)+nmat* &
          (num_state_material+ngeom_raw)+1) then
       ! do nothing
      else 
       print *,"nstate invalid"
       stop
      endif

      if (num_materials_scalar_solve.eq.1) then ! GFM
       normal_probe_factor=half
      else if (num_materials_scalar_solve.eq.nmat) then ! FVM multimat
       normal_probe_factor=half
      else
       print *,"num_materials_scalar_solve invalid"
       stop
      endif
      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,1251)
      call checkbound(fablo,fabhi,DIMS(LSnew),1,-1,1253)
      call checkbound(fablo,fabhi,DIMS(Snew),1,-1,1253)
      call checkbound(fablo,fabhi,DIMS(EOS),ngrow,-1,1254)
      call checkbound(fablo,fabhi,DIMS(pres),ngrow,-1,1255)
      call checkbound(fablo,fabhi,DIMS(pres_eos),1,-1,1255)


      if (nucleation_flag.eq.0) then

       call checkbound(fablo,fabhi,DIMS(typefab),1,-1,6625)
       call checkbound(fablo,fabhi,DIMS(colorfab),1,-1,6626)

       ! nmat x (sdim+1) components
       call checkbound(fablo,fabhi, &
        DIMS(burnvel), &
        ngrow_make_distance,-1,1250)
       call checkbound(fablo,fabhi, &
        DIMS(Tsatfab), &
        ngrow_make_distance,-1,1250)

       call checkbound(fablo,fabhi, &
        DIMS(curvfab), &
        ngrow_make_distance,-1,1250)

       call checkbound(fablo,fabhi,DIMS(recon),ngrow,-1,1251)
       call checkbound(fablo,fabhi,DIMS(LS),ngrow,-1,1252)
       call checkbound(fablo,fabhi,DIMS(LS_slopes_FD),1,-1,1253)

      else if (nucleation_flag.eq.1) then
       ! do nothing
      else
       print *,"nucleation_flag invalid"
       stop
      endif

      !blob_matrix,blob_RHS,blob_velocity,
      !blob_integral_momentum,blob_energy,
      !blob_mass_for_velocity (3 comp)
      !blob_volume, 
      !blob_center_integral,blob_center_actual
      !blob_perim, blob_perim_mat, blob_triple_perim, 
      !blob_cell_count
      !blob_cellvol_count
      !blob_mass
      if (num_elements_blobclass.ne. &
          3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
          2*(2*SDIM)+1+ &
          3+ & ! blob_mass_for_velocity
          1+ & ! volume
          2*SDIM+ & ! centroid integral, actual
          1+ & ! perim 
          nmat+ & ! perim_mat
          nmat*nmat+ & ! blob_triple_perim 
          1+1+1) then
       print *,"num_elements_blobclass invalid rate mass change:", &
         num_elements_blobclass
       print *,"blob_cell_count readded Febrary 11, 2020"
       print *,"blob_cellvol_count added December 6, 2020"
       print *,"blob_mass added January 23, 2021"
       stop
      endif

      if (STANDALONE.eq.0) then
       if (arraysize.ne.num_elements_blobclass*color_count) then
        print *,"arraysize invalid rate mass change (get stat==1)"
        print *,"arraysize=",arraysize
        print *,"num_elements_blobclass=",num_elements_blobclass
        print *,"color_count=",color_count
        stop
       endif
      else if (STANDALONE.eq.1) then
       ! check nothing
      else
       print *,"STANDALONE invalid"
       stop
      endif

         ! SANDIPAN HOOK HERE
         ! pseudo code:
         ! if typefab(D_DECL(i,j,k))=im_vapor then
         !  color = colorfab(D_DECL(i,j,k))
         !  vapor bubble statistics are in 
         !   blob_array((color-1)*num_elements_blobclass + l)
         !  l=1..num_elements_blobclass
         ! for example l=3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+
         !               2*(2*SDIM)+1+
         !               3+1  corresponds to blob_volume.
         ! MDOT=(1-den_vapor/den_liquid)*(F^Vapor_new - F^Vapor_old)*
         !  volume/dt^2
         ! =(1-den_vapor/den_liquid)*(Vvapor_new - Vvapor_old)/dt^2
         ! MDOT should have the same units as volume * div u/dt
         ! units of volume * div u/dt = m^3 (1/m)(m/s)(1/s)=1/s^2
         ! for standard phase change:
         ! units: cm^3/s^2
         ! jump_strength=(denratio_factor/dt)*dF*volgrid/dt
         ! if evaporation then expansion_term should be positive.
         ! FOR CODY: if pressure falls below some cavitation pressure, then
         ! material is cavitated. (velocity of phase change is dx/dt?)


      call get_dxmin(dx,bfact,dxmin)
      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do iten=1,2*nten
       temp_target_probe_history(iten,1)=zero
       dxprobe_target_history(iten,1)=zero
       temp_target_probe_history(iten,2)=zero
       dxprobe_target_history(iten,2)=zero
      enddo
        ! copy_dimdec(dest,source), in: GLOBALUTIL.F90
      call copy_dimdec( &
        DIMS(create_in%EOS), &
        DIMS(EOS))
      call copy_dimdec( &
        DIMS(create_in%pres), &
        DIMS(pres))
      call copy_dimdec( &
        DIMS(create_in%pres_eos), &
        DIMS(pres_eos))
      call copy_dimdec( &
        DIMS(create_in%Snew), &
        DIMS(Snew))
      call copy_dimdec( &
        DIMS(create_in%LSnew), &
        DIMS(Snew))
      create_in%EOS=>EOS
      create_in%pres=>pres
      create_in%pres_eos=>pres_eos
      create_inout%Snew=>Snew
      create_inout%LSnew=>LSnew

      create_in%tid=tid
      create_in%dxmaxLS=dxmaxLS
      create_in%bfact=bfact
      create_in%level=level
      create_in%finest_level=finest_level
      create_in%dx=>dx
      create_in%xlo=>xlo
      create_in%nmat=nmat
      create_in%nten=nten
      create_in%nstate=nstate
      create_in%fablo=>fablo
      create_in%fabhi=>fabhi
      create_in%custom_nucleation_model=custom_nucleation_model
      create_in%do_the_nucleate=do_the_nucleate
      create_in%nucleate_pos_size=nucleate_pos_size
      create_in%nucleate_pos=>nucleate_pos
      create_in%nucleation_temp=>nucleation_temp
      create_in%nucleation_pressure=>nucleation_pressure
      create_in%nucleation_pmg=>nucleation_pmg
      create_in%nucleation_mach=>nucleation_mach
      create_in%cavitation_pressure=>cavitation_pressure
      create_in%cavitation_vapor_density=>cavitation_vapor_density
      create_in%cavitation_tension=>cavitation_tension
      create_in%prev_time=prev_time
      create_in%cur_time=cur_time
      create_in%dt=dt

        ! copy_dimdec(dest,source), in: GLOBALUTIL.F90
      call copy_dimdec( &
        DIMS(PROBE_PARMS%EOS), &
        DIMS(EOS))
      call copy_dimdec( &
        DIMS(PROBE_PARMS%recon), &
        DIMS(recon))
      call copy_dimdec( &
        DIMS(PROBE_PARMS%LS), &
        DIMS(LS))
      call copy_dimdec( &
        DIMS(PROBE_PARMS%pres), &
        DIMS(pres))

      PROBE_PARMS%tid=tid
      PROBE_PARMS%use_supermesh=use_supermesh

      PROBE_PARMS%Y_TOLERANCE=>Y_TOLERANCE

      PROBE_PARMS%debugrate=>debugrate
      PROBE_PARMS%EOS=>EOS 
      PROBE_PARMS%LS=>LS  ! PROBE_PARMS%LS is pointer, LS is target
      PROBE_PARMS%recon=>recon
      PROBE_PARMS%pres=>pres
      PROBE_PARMS%dxmaxLS=>dxmaxLS
      PROBE_PARMS%bfact=>bfact
      PROBE_PARMS%level=>level
      PROBE_PARMS%finest_level=>finest_level
      PROBE_PARMS%dx=>dx
      PROBE_PARMS%xlo=>xlo
      PROBE_PARMS%nmat=>nmat
      PROBE_PARMS%ngrow=>ngrow
      PROBE_PARMS%fablo=>fablo
      PROBE_PARMS%fabhi=>fabhi

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       create_in%i=i
       create_in%j=j
       create_in%k=k

       call gridsten_level(xsten,i,j,k,level,nhalf)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        if (nucleation_flag.eq.0) then

         ! LEVELSET FUNCTION AT CELL CENTERS.
         do im=1,nmat
          LShere(im)=LS(D_DECL(i,j,k),im)
         enddo
         call get_primary_material(LShere,nmat,im_primary)

         if (is_rigid(nmat,im_primary).eq.0) then

          do im=1,nmat-1
           do im_opp=im+1,nmat

            if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
             print *,"im or im_opp bust 9"
             stop
            endif
            call get_iten(im,im_opp,iten,nmat)

            do ireverse=0,1

             temp_target_probe_history(iten+ireverse*nten,1)=zero
             dxprobe_target_history(iten+ireverse*nten,1)=zero
             temp_target_probe_history(iten+ireverse*nten,2)=zero
             dxprobe_target_history(iten+ireverse*nten,2)=zero

             valid_phase_change(ireverse)=0
             LL(ireverse)=latent_heat(iten+ireverse*nten)
             K_f(ireverse)=reaction_rate(iten+ireverse*nten)

             local_freezing_model=freezing_model(iten+ireverse*nten)
             local_Tanasawa_or_Schrage_or_Kassemi= &
               Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*nten)

             ispec=mass_fraction_id(iten+ireverse*nten)
             if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
              molar_mass_vapor=species_molar_mass(ispec)
             else if (ispec.eq.0) then
              molar_mass_vapor=zero
             else
              print *,"ispec invalid"
              stop
             endif

             if (ireverse.eq.0) then
              im_source=im
              im_dest=im_opp
             else if (ireverse.eq.1) then
              im_source=im_opp
              im_dest=im
             else
              print *,"ireverse invalid"
              stop
             endif

             dencomp_source=(im_source-1)*num_state_material+1
             dencomp_dest=(im_dest-1)*num_state_material+1

             if ((ispec.ge.0).and.(ispec.le.num_species_var)) then
              ! do nothing
             else
              print *,"ispec invalid"
              stop
             endif

             vapor_den=fort_denconst(im_dest) ! default

             if (LL(ireverse).eq.zero) then
              ! do nothing
             else if (LL(ireverse).gt.zero) then
              if (constant_density_all_time(im_dest).eq.1) then
               vapor_den=fort_denconst(im_dest) 
              else if (constant_density_all_time(im_dest).eq.0) then
               vapor_den=EOS(D_DECL(i,j,k),dencomp_dest)
              else
               print *,"constant_density_all_time(im_dest) invalid"
               stop
              endif
             else if (LL(ireverse).lt.zero) then
              if (constant_density_all_time(im_source).eq.1) then
               vapor_den=fort_denconst(im_source) 
              else if (constant_density_all_time(im_source).eq.0) then
               vapor_den=EOS(D_DECL(i,j,k),dencomp_source)
              else
               print *,"constant_density_all_time(im_source) invalid"
               stop
              endif
             else
              print *,"LL invalid"
              stop
             endif

             if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
              ! do nothing
             else if (ispec.eq.0) then
              if ((local_freezing_model.eq.2).or. & !hydrate
                  (local_freezing_model.eq.4).or. & !Tanasawa or Schrage
                  (local_freezing_model.eq.5).or. & !stefan evap/cond
                  (local_freezing_model.eq.6)) then !Palmore/Desjardins
               print *,"ispec invalid"
               stop
              endif
             else
              print *,"ispec invalid"
              stop
             endif

             distribute_from_targ=distribute_from_target(iten+ireverse*nten)

             if ((distribute_from_targ.ne.0).and. &
                 (distribute_from_targ.ne.1)) then
              print *,"distribute_from_targ invalid"
              stop
             endif

             if (LL(ireverse).ne.zero) then
            
              if (im_dest.lt.im_source) then
               LSSIGN=one
              else if (im_dest.gt.im_source) then
               LSSIGN=-one
              else
               print *,"im_dest<>im_source not satisfied"
               stop
              endif

              if ((is_rigid(nmat,im).eq.1).or. &
                  (is_rigid(nmat,im_opp).eq.1)) then

               ! do nothing

              else if (LL(ireverse).ne.zero) then

               local_hardwire_T(ireverse)=hardwire_T_gamma(iten+ireverse*nten)
               local_hardwire_Y(ireverse)=hardwire_Y_gamma(iten+ireverse*nten)
               if ((local_hardwire_T(ireverse).gt.zero).and. &
                   (local_hardwire_Y(ireverse).gt.zero)) then
                hardwire_flag(ireverse)=1
               else if ((local_hardwire_T(ireverse).eq.zero).and. &
                        (local_hardwire_Y(ireverse).eq.zero)) then
                hardwire_flag(ireverse)=0
               else
                print *,"local_hardwire vars invalid"
                stop
               endif

               local_Tsat(ireverse)=saturation_temp(iten+ireverse*nten)
               local_Tsat_base(ireverse)=saturation_temp(iten+ireverse*nten)

               found_path=0

                ! FOR YANG:
                ! BOTH LEVELSET FUNCTIONS WITHIN 2 dx of interface and
                ! at least one of them is positive.
                ! note: is_rigid(nmat,im_primary).eq.0 so we are not in 
                ! a solid material.
               if ((abs(LShere(im_source)).le.two*dxmaxLS).and. &
                   (abs(LShere(im_dest)).le.two*dxmaxLS)) then

                 ! The phase change velocity is defined at the cell
                 ! centers, therefore, we find the jump in heat
                 ! flux at the closest point on the interface from the
                 ! cell center (xsten(0,dir)).
                if ((LShere(im_dest).ge.zero).or. &
                    (LShere(im_source).ge.zero)) then

                 if (LShere(im_dest).ge.zero) then
                  LS_pos=LShere(im_source)
                  do dir=1,SDIM
                     ! xCP=x-phi grad phi   grad phi=(x-xCP)/phi
                   nrmCP(dir)=LS(D_DECL(i,j,k),nmat+(im_source-1)*SDIM+dir)
                     ! Least squares slope: see Sussman and Puckett (2000)
                   nrmFD(dir)= &
                    LS_slopes_FD(D_DECL(i,j,k),(im_source-1)*SDIM+dir)
                   xI(dir)=xsten(0,dir)-LS_pos*nrmCP(dir)

                   nrmPROBE(dir)=nrmCP(dir)

                    ! e.g.
                    ! normal_probe_factor=1/2
                    ! normal_probe_size=1
                   xdst(dir)=xI(dir)- &
                    normal_probe_factor*normal_probe_size*dxmin*nrmPROBE(dir) 
                   xsrc(dir)=xI(dir)+ &
                    normal_probe_factor*normal_probe_size*dxmin*nrmPROBE(dir)
                   xdst_micro(dir)=xI(dir)- &
                    microscale_probe_size*dxmin*nrmPROBE(dir) 
                   xsrc_micro(dir)=xI(dir)+ &
                    microscale_probe_size*dxmin*nrmPROBE(dir)
                  enddo ! dir=1..sdim
                 else if (LShere(im_source).ge.zero) then
                  LS_pos=LShere(im_dest)
                  do dir=1,SDIM
                   nrmCP(dir)=LS(D_DECL(i,j,k),nmat+(im_dest-1)*SDIM+dir)
                   nrmFD(dir)= &
                    LS_slopes_FD(D_DECL(i,j,k),(im_dest-1)*SDIM+dir)

                   nrmPROBE(dir)=nrmCP(dir)

                   xI(dir)=xsten(0,dir)-LS_pos*nrmCP(dir)
                   xdst(dir)=xI(dir)+ &
                      normal_probe_factor*normal_probe_size*dxmin*nrmPROBE(dir) 
                   xsrc(dir)=xI(dir)- &
                      normal_probe_factor*normal_probe_size*dxmin*nrmPROBE(dir)
                   xdst_micro(dir)=xI(dir)+ &
                      microscale_probe_size*dxmin*nrmPROBE(dir) 
                   xsrc_micro(dir)=xI(dir)- &
                      microscale_probe_size*dxmin*nrmPROBE(dir)
                  enddo ! dir=1..sdim
                 else
                  print *,"LShere bust"
                  print *,"LShere(im_dest) ",LShere(im_dest)
                  print *,"LShere(im_source) ",LShere(im_source)
                  stop
                 endif

                 found_path=1

                else if ((LShere(im_dest).lt.zero).and. &
                         (LShere(im_source).lt.zero)) then

                 found_path=0

                else
                 print *,"LShere bust"
                 print *,"LShere(im_dest) ",LShere(im_dest)
                 print *,"LShere(im_source) ",LShere(im_source)
                 stop
                endif
            
                if (found_path.eq.1) then  

                 vofcomp_source=(im_source-1)*ngeom_recon+1
                 vofcomp_dest=(im_dest-1)*ngeom_recon+1
                 Fsource=recon(D_DECL(i,j,k),vofcomp_source)
                 Fdest=recon(D_DECL(i,j,k),vofcomp_dest)

                 C_w0=fort_denconst(1)  ! density of water
                 pres_I_interp(1)=2.0D+19 ! source
                 pres_I_interp(2)=2.0D+19 ! dest
                 Y_I_interp(1)=zero
                 Y_I_interp(2)=zero ! destination, C_methane_in_hyd

                 tcomp_source=(im_source-1)*num_state_material+2
                 tcomp_dest=(im_dest-1)*num_state_material+2

                 im_ambient=0
                 Ycomp_source=0
                 Ycomp_dest=0

                 if ((local_freezing_model.eq.2).or. & !hydrate
                     (local_freezing_model.eq.4).or. & !Tanasawa or Schrage
                     (local_freezing_model.eq.5).or. & !Stefan evap/cond
                     (local_freezing_model.eq.6).or. & !Palmore/Desjardins
                     (local_freezing_model.eq.7)) then !Cavitation
                  if (LL(ireverse).gt.zero) then ! evaporation
                   im_ambient=im_dest
                  else if (LL(ireverse).lt.zero) then ! condensation
                   im_ambient=im_source
                  else
                   print *,"LL invalid"
                   stop
                  endif
                    
                  if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
                   Ycomp_source=(im_source-1)*num_state_material+2+ispec
                   Ycomp_dest=(im_dest-1)*num_state_material+2+ispec
                   FicksLawD(1)= &
                    fort_speciesviscconst((ispec-1)*nmat+im_source)
                   FicksLawD(2)= &
                    fort_speciesviscconst((ispec-1)*nmat+im_dest)
                  else
                   print *,"ispec invalid"
                   stop
                  endif
                 else if ((local_freezing_model.eq.0).or. & ! Stefan model
                          (local_freezing_model.eq.1).or. & ! source term
                          (local_freezing_model.eq.3)) then ! wild fire

                  if (ispec.eq.0) then
                   Ycomp_source=0
                   Ycomp_dest=0
                  else if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
                   print *,"expecting ispec ==0 "
                   stop
                  endif
                 else
                  print *,"local_freezing_model invalid 3"
                  stop
                 endif

                 call interpfab_curv( &
                  nmat+iten, &
                  nten, &
                  nmat, &
                  bfact, &
                  level, &
                  finest_level, &
                  dx,xlo, &
                  xI, &
                  ngrow_make_distance, &
                  fablo,fabhi, &
                  curvfab,DIMS(curvfab), &
                  CURV_OUT_I)

                 use_tsatfab=0

                  ! use_exact_temperature==0 if regular Stefan problem.
                  ! subroutine get_interface_temperature defined here in
                  ! MASS_TRANSFER_3D.F90
                 if (hardwire_flag(ireverse).eq.0) then
                  call get_interface_temperature( &
                   use_tsatfab, &
                   i,j,k, &
                   ireverse, &
                   iten, &
                   ntsat, &
                   bfact, &
                   level, &
                   finest_level, &
                   dx,xlo, &
                   ngrow_make_distance, & ! ngrow_tsat
                   fablo,fabhi, &
                   Tsatfab,DIMS(Tsatfab), & ! not used since use_tsatfab==0
                   local_Tsat(ireverse), & !user def. interface temperature 
                   iten+ireverse*nten, &
                   saturation_temp, &
                   use_exact_temperature, &
                   xI,cur_time,nmat,nten,7)
                 else if (hardwire_flag(ireverse).eq.1) then
                  local_Tsat(ireverse)=local_hardwire_T(ireverse)
                 else
                  print *,"hardwire_flag(ireverse) invalid"
                  stop
                 endif

                 if (local_Tsat(ireverse).ge.zero) then
                  ! do nothing
                 else
                  print *,"local_Tsat(ireverse) bad:",local_Tsat(ireverse)
                  print *,"ireverse=",ireverse
                  stop
                 endif

                 dxprobe_source=zero
                 do dir=1,SDIM
                  dxprobe_source=dxprobe_source+(xdst(dir)-xsrc(dir))**2
                 enddo
                 dxprobe_source=half*sqrt(dxprobe_source)
                 dxprobe_dest=dxprobe_source

                 RR=one
                 call prepare_normal(nrmPROBE,RR,mag)

                 if (levelrz.eq.0) then
                  ! do nothing
                 else if (levelrz.eq.1) then
                  if (SDIM.ne.2) then
                   print *,"dimension bust"
                   stop
                  endif
                 else if (levelrz.eq.3) then
                  if (mag.gt.zero) then
                   do dir=1,SDIM
                    theta_nrmPROBE(dir)=nrmPROBE(dir)
                   enddo
                   RR=xsten(0,1)
                   call prepare_normal(theta_nrmPROBE,RR,mag)
                   if (mag.gt.zero) then
                    ! mag=theta_nrmPROBE dot nrmPROBE
                    mag=zero
                    do dir=1,SDIM
                     mag=mag+theta_nrmPROBE(dir)*nrmPROBE(dir)
                    enddo
                    if (abs(mag).gt.one+VOFTOL) then
                     print *,"dot product bust"
                     stop
                    endif
                    if (abs(mag).gt.zero) then
                     dxprobe_source=dxprobe_source*abs(mag)
                     dxprobe_dest=dxprobe_dest*abs(mag)
                    endif
                   else if (mag.eq.zero) then
                    ! do nothing
                   else
                    print *,"mag bust"
                    stop
                   endif 
                  else if (mag.eq.zero) then
                   ! do nothing
                  else
                   print *,"mag bust"
                   stop
                  endif 
                 else
                  print *,"levelrz invalid rate mass change"
                  stop
                 endif

                 do imls=1,nmat*(SDIM+1)
                  call interpfab( &
                   bfact, &
                   level, &
                   finest_level, &
                   dx, &
                   xlo,xI, &
                   imls, &
                   ngrow, &
                   fablo,fabhi, &
                   LS,DIMS(LS), &
                   LSINT(imls))
                 enddo ! imls=1..nmat*(SDIM+1)

                 call get_primary_material(LSINT,nmat,imls_I)

                 TI_min=saturation_temp_min(iten+ireverse*nten)
                 TI_max=saturation_temp_max(iten+ireverse*nten)

                 if ((local_freezing_model.eq.4).or. & !Tanasawa or Schrage
                     (local_freezing_model.eq.5).or. & !Stefan evap/cond
                     (local_freezing_model.eq.6)) then !Palmore/Desjardins

                  if (LL(ireverse).gt.zero) then ! evaporation
                   TI_max=local_Tsat(ireverse)
                  else if (LL(ireverse).lt.zero) then ! condensation
                   TI_min=local_Tsat(ireverse)
                  else
                   print *,"LL(ireverse) invalid"
                   stop
                  endif

                 else if ((local_freezing_model.eq.0).or. & ! Stefan model
                          (local_freezing_model.eq.1).or. & ! source term
                          (local_freezing_model.eq.2).or. & ! hydrate
                          (local_freezing_model.eq.3).or. & ! wildfire
                          (local_freezing_model.eq.7)) then ! cavitation
                  ! do nothing
                 else
                  print *,"local_freezing_model invalid"
                  stop
                 endif
                          
                 if ((TI_max.ge.TI_min).and.(TI_min.ge.zero)) then
                  ! do nothing
                 else
                  print *,"TI_max or TI_min invalid"
                  stop
                 endif

                 PROBE_PARMS%dxprobe_source=>dxprobe_source
                 PROBE_PARMS%dxprobe_dest=>dxprobe_dest
                 PROBE_PARMS%local_freezing_model=>local_freezing_model
                 PROBE_PARMS%LL=>LL(ireverse)
                 PROBE_PARMS%i=>i
                 PROBE_PARMS%j=>j
                 PROBE_PARMS%k=>k
                 PROBE_PARMS%xsrc=>xsrc 
                 PROBE_PARMS%xdst=>xdst
                 PROBE_PARMS%xsrc_micro=>xsrc_micro 
                 PROBE_PARMS%xdst_micro=>xdst_micro
                 PROBE_PARMS%LSINT=>LSINT
                 PROBE_PARMS%imls_I=>imls_I
                 PROBE_PARMS%im_source=>im_source
                 PROBE_PARMS%im_dest=>im_dest
                 PROBE_PARMS%tcomp_source=>tcomp_source
                 PROBE_PARMS%Ycomp_source=>Ycomp_source
                 PROBE_PARMS%dencomp_source=>dencomp_source
                 PROBE_PARMS%tcomp_dest=>tcomp_dest
                 PROBE_PARMS%Ycomp_dest=>Ycomp_dest
                 PROBE_PARMS%dencomp_dest=>dencomp_dest
                 PROBE_PARMS%xI=>xI

#if (STANDALONE==0)
                 thermal_k(1)=get_user_heatviscconst(im_source)
                 thermal_k(2)=get_user_heatviscconst(im_dest)
#elif (STANDALONE==1)
                 thermal_k(1)=fort_heatviscconst(im_source)
                 thermal_k(2)=fort_heatviscconst(im_dest)
#else
                 print *,"bust compiling ratemasschange"
                 stop
#endif

                 if (LL(ireverse).gt.zero) then ! evaporation
                  iprobe_vapor=2  ! destination
                  molar_mass_ambient=molar_mass(im_dest)
                 else if (LL(ireverse).lt.zero) then ! condensation
                  iprobe_vapor=1  ! source
                  molar_mass_ambient=molar_mass(im_source)
                 else
                  print *,"LL invalid"
                  stop
                 endif


                  ! for T_INTERFACE= TSAT - eps1 K - eps2 V
                  ! 1. T_I^(0)=TSAT - eps1 K - eps2 ( 0 )
                  ! 2. n=0, V^{0}=0
                  ! 2. while not converged
                  ! 3.  V^{n+1}=-[k grad T dot n]/L=
                  !     [(Tprobe^1(T_I^{n}) - T_I^{n})/dxprobe^1 -
                  !     (Tprobe^2(T_I^{n}) - T_I^{n})/dxprobe^2 ]/L
                  !     n points from material 1 to material 2.
                  ! 4.  T_I^{n+1}=T_I^{n} - eps2 (V^{n+1}-V^{n})
                  ! 5.  convergence when |T_{I}^{n+1}-T_{I}^{n}|<tol
                  !
                 if (hardwire_flag(ireverse).eq.0) then

                  delta_Tsat= &
                    saturation_temp_curv(iten+ireverse*nten)*CURV_OUT_I

                   ! T_Gamma = Tsat - eps1 * kappa
                  local_Tsat(ireverse)=local_Tsat(ireverse)-delta_Tsat

                  if (local_Tsat(ireverse).lt.TI_min) then
                   local_Tsat(ireverse)=TI_min
                  endif
                  if (local_Tsat(ireverse).gt.TI_max) then
                   local_Tsat(ireverse)=TI_max
                  endif

                  Y_predict=one

                 else if (hardwire_flag(ireverse).eq.1) then

                  local_Tsat(ireverse)=local_hardwire_T(ireverse)
                  Y_predict=local_hardwire_Y(ireverse)

                 else
                  print *,"hardwire_flag(ireverse) invalid"
                  stop
                 endif
               
                 !iprobe=1 source
                 !iprobe=2 dest
                 call probe_interpolation( &
                  PROBE_PARMS, &
                  local_Tsat(ireverse), &
                  Y_predict, &
                  T_probe,Y_probe, &
                  den_I_interp, &
                  den_probe, &
                  T_I_interp,Y_I_interp, &
                  pres_I_interp, &
                  vfrac_I, &
                  T_probe_raw, &
                  dxprobe_target, &
                  interp_valid_flag_initial, &
                  at_interface)

                 interface_resolved=1

                 if ((interp_valid_flag_initial(1).eq.1).and. &
                     (interp_valid_flag_initial(2).eq.1)) then
                  ! do nothing (valid probe on both sides of Gamma)
                 else if ((interp_valid_flag_initial(1).eq.0).or. &
                          (interp_valid_flag_initial(1).eq.2).or. &
                          (interp_valid_flag_initial(2).eq.0).or. &
                          (interp_valid_flag_initial(2).eq.2)) then
                  interface_resolved=0
                 else
                  print *,"interp_valid_flag_initial invalid"
                  stop
                 endif

                 TSAT_predict=local_Tsat(ireverse)

                 TSAT_correct=TSAT_predict

                 VEL_predict=zero
                 VEL_correct=zero

                 TSAT_iter=0
                 TSAT_converge=0


!    local_freezing_model=4  Tanasawa or Schrage
!    local_freezing_model=5  fully saturated evaporation?
!    local_freezing_model=6  Palmore/Desjardins
!    local_freezing_model=7  Cavitation (a seed must exist)

                 do while (TSAT_converge.eq.0) 

                  !iprobe=1 source
                  !iprobe=2 dest
                  call probe_interpolation( &
                   PROBE_PARMS, &
                   TSAT_predict, &
                   Y_predict, &
                   T_probe,Y_probe, &
                   den_I_interp, &
                   den_probe, &
                   T_I_interp,Y_I_interp, &
                   pres_I_interp, &
                   vfrac_I, &
                   T_probe_raw, &
                   dxprobe_target, &
                   interp_valid_flag, &
                   at_interface)

                  if (TSAT_iter.eq.0) then
                   fully_saturated=0

                   T_gamma_a=TSAT_predict
                   T_gamma_b=TSAT_predict
                   Y_gamma_a=Y_predict
                   Y_gamma_b=Y_predict

                   if (at_interface.eq.1) then
                    if (hardwire_flag(ireverse).eq.0) then
                     if (local_Tanasawa_or_Schrage_or_Kassemi.eq.3) then
                      fully_saturated=2
                      Y_predict=one
                      X_predict=one
                      T_gamma_a=TSAT_predict
                      T_gamma_b=TSAT_predict
                      Y_gamma_a=Y_predict
                      Y_gamma_b=Y_predict
                     else if ((Y_probe(iprobe_vapor).ge.one-Y_TOLERANCE).and. &
                              (Y_probe(iprobe_vapor).le.one).and. &
                              (Y_predict.eq.one)) then
                      fully_saturated=1
                     else if ((Y_probe(iprobe_vapor).le.one-Y_TOLERANCE).and. &
                              (Y_probe(iprobe_vapor).ge.zero).and. &
                              (Y_predict.eq.one)) then
                       ! declared in PROBCOMMON.F90
                      Y_predict=one-EVAP_BISECTION_TOL
                      call volfrac_from_massfrac(X_predict,Y_predict, &
                        molar_mass_ambient,molar_mass_vapor) ! WA,WV
                      call Tgamma_from_TSAT_and_X(TSAT_predict, &
                       local_Tsat(ireverse), &
                       X_predict,LL(ireverse),R_Palmore_Desjardins, &
                       molar_mass_vapor,TI_min,TI_max)

                      T_gamma_a=TSAT_predict
                      T_gamma_b=TSAT_predict
                      Y_gamma_a=Y_predict
                      Y_gamma_b=Y_predict
                     else
                      print *,"mass fraction bust"
                      stop
                     endif

                     if (fully_saturated.eq.1) then
                      ! do nothing
                     else if (fully_saturated.eq.2) then
                      if (LL(ireverse).gt.zero) then ! evaporation
                       T_gamma_a=TI_min
                       X_gamma_a=one
                       Y_gamma_a=one
                      else if (LL(ireverse).lt.zero) then ! condensation
                       T_gamma_b=TI_max
                       X_gamma_b=one
                       Y_gamma_b=one
                      else
                       print *,"LL(ireverse) invalid"
                       stop
                      endif
                     else if (fully_saturated.eq.0) then
                      if (LL(ireverse).gt.zero) then ! evaporation
                       T_gamma_a=TI_min
                       call X_from_Tgamma(X_gamma_a,T_gamma_a, &
                         local_Tsat(ireverse), &
                         LL(ireverse),R_Palmore_Desjardins, &
                         molar_mass_vapor) ! WV
                       call massfrac_from_volfrac(X_gamma_a,Y_gamma_a, &
                        molar_mass_ambient,molar_mass_vapor) ! WA,WV
                      else if (LL(ireverse).lt.zero) then ! condensation
                       T_gamma_b=TI_max
                       call X_from_Tgamma(X_gamma_b,T_gamma_b, &
                         local_Tsat(ireverse), &
                         LL(ireverse),R_Palmore_Desjardins, &
                         molar_mass_vapor) ! WV
                       call massfrac_from_volfrac(X_gamma_b,Y_gamma_b, &
                        molar_mass_ambient,molar_mass_vapor) ! WA,WV
                      else
                       print *,"LL(ireverse) invalid"
                       stop
                      endif
                     else
                      print *,"fully_saturated invalid"
                      stop
                     endif

                    else if (hardwire_flag(ireverse).eq.1) then

                     if (Y_predict.eq.one) then
                      fully_saturated=1
                     else if ((Y_predict.ge.zero).and. &
                              (Y_predict.lt.one)) then
                      ! do nothing
                     else
                      print *,"Y_predict invalid"
                      stop
                     endif

                    else
                     print *,"hardwire_flag(ireverse) invalid"
                     stop
                    endif
                   else if (at_interface.eq.0) then
                    ! do nothing
                   else
                    print *,"at_interface invalid FORT_RATEMASSCHANGE"
                    print *,"at_interface=",at_interface
                    stop
                   endif

                   T_gamma_a_init=T_gamma_a
                   T_gamma_b_init=T_gamma_b
                   Y_gamma_a_init=Y_gamma_a
                   Y_gamma_b_init=Y_gamma_b
                  else if (TSAT_iter.gt.0) then
                   ! do nothing
                  else
                   print *,"TSAT_iter invalid"
                   stop
                  endif

                  den_I_interp_SAT(1)=den_I_interp(1) ! iprobe=1 source
                  den_I_interp_SAT(2)=den_I_interp(2) ! iprobe=2 dest

                  !iprobe=1 source
                  !iprobe=2 dest

                  if (at_interface.eq.1) then
                      
                   microlayer_substrate_source=0
                   microlayer_substrate_dest=0

                   im_substrate_source=microlayer_substrate(im_source)
                   im_substrate_dest=microlayer_substrate(im_dest)

                   if (microlayer_size(im_source).gt.zero) then
                    if (im_substrate_source.gt.0) then
                     if (is_rigid(nmat,im_substrate_source).ne.1) then
                      print *,"is_rigid(nmat,im_substrate_source).ne.1"
                      stop
                     endif
                     do dir=1,SDIM
                      dir2=nmat+(im_substrate_source-1)*SDIM+dir
                      gradphi_substrate(dir)=LSINT(dir2)
                     enddo
                     call get_physical_dist(xI,LSINT(im_substrate_source), &
                      gradphi_substrate,newphi_substrate);
                     if (newphi_substrate.ge.-macrolayer_size(im_source)) then
                      microlayer_substrate_source=im_substrate_source
                     endif
                    else if (im_substrate_source.eq.0) then
                     print *,"microlayer material must have corresponding"
                     print *,"substrate. im_source,im_dest=",im_source,im_dest
                     stop
                    else
                     print *,"im_substrate_source invalid"
                     stop
                    endif
                   else if (microlayer_size(im_source).eq.zero) then
                    ! do nothing
                   else
                    print *,"microlayer_size(im_source) invalid"
                    stop
                   endif

                   if (microlayer_size(im_dest).gt.zero) then
                    if (im_substrate_dest.gt.0) then
                     if (is_rigid(nmat,im_substrate_dest).ne.1) then
                      print *,"is_rigid(nmat,im_substrate_dest).ne.1"
                      stop
                     endif
                     do dir=1,SDIM
                      dir2=nmat+(im_substrate_dest-1)*SDIM+dir
                      gradphi_substrate(dir)=LSINT(dir2)
                     enddo
                     call get_physical_dist(xI,LSINT(im_substrate_dest), &
                      gradphi_substrate,newphi_substrate);
                     if (newphi_substrate.ge.-macrolayer_size(im_dest)) then
                      microlayer_substrate_dest=im_substrate_dest
                     endif
                    else if (im_substrate_dest.eq.0) then
                     print *,"microlayer material must have corresponding"
                     print *,"substrate. im_source,im_dest=",im_source,im_dest
                     stop
                    else
                     print *,"im_substrate_dest invalid"
                     stop
                    endif
                   else if (microlayer_size(im_dest).eq.zero) then
                    ! do nothing
                   else
                    print *,"microlayer_size(im_dest) invalid"
                    stop
                   endif

                   if ((local_freezing_model.eq.5).or. & !stefan evap/cond
                       (local_freezing_model.eq.6)) then !Palmore/Desjardins
                    if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
                     if (vapor_den.gt.zero) then
                      if (LL(ireverse).gt.zero) then ! evaporation
                       den_I_interp_SAT(2)=vapor_den ! dest
                      else if (LL(ireverse).lt.zero) then ! condensation
                       den_I_interp_SAT(1)=vapor_den ! source
                      else
                       print *,"LL invalid"
                       stop
                      endif
                     else
                      print *,"vapor_den invalid"
                      stop
                     endif  
                    else
                     print *,"ispec invalid"
                     stop
                    endif
                   else if (local_freezing_model.eq.4) then ! Tanasawa/Schrage
                    ! do nothing
                   else if (local_freezing_model.eq.7) then ! cavitation
                    ! do nothing
                   else if (local_freezing_model.eq.3) then ! wild fire
                    ! do nothing
                   else if (local_freezing_model.eq.2) then ! hydrate
                    ! do nothing
                   else if (local_freezing_model.eq.1) then ! source term
                    ! do nothing
                   else if (local_freezing_model.eq.0) then ! stefan cond
                    ! do nothing
                   else
                    print *,"local_freezing_model invalid 14"
                    stop
                   endif 

                   source_perim_factor=one
                   dest_perim_factor=one

                   do iprobe=1,2

                    if (iprobe.eq.1) then
                     im_probe=im_source
                     microlayer_substrate_probe=microlayer_substrate_source
                    else if (iprobe.eq.2) then
                     im_probe=im_dest
                     microlayer_substrate_probe=microlayer_substrate_dest
                    else
                     print *,"iprobe invalid"
                     stop
                    endif
 
                    if ((local_freezing_model.eq.0).and. &
                        (max_contact_line_size(im_probe).gt.zero).and. &
                        (microlayer_size(im_probe).gt.zero).and. &
                        (macrolayer_size(im_probe).gt.zero).and. &
                        (microlayer_substrate_probe.ge.1).and. &
                        (microlayer_substrate_probe.le.nmat)) then
                     icolor=NINT(colorfab(D_DECL(i,j,k)))
                     if ((icolor.gt.color_count).or.(icolor.le.0)) then
                      print *,"icolor invalid in RATEMASSCHANGE icolor=",icolor
                      print *,"i,j,k ",i,j,k
                      stop
                     endif
                     base_type=NINT(typefab(D_DECL(i,j,k)))
                     if ((base_type.lt.1).or.(base_type.gt.nmat)) then
                      print *,"base_type invalid"
                      stop
                     endif
                     ! blob_matrix,blob_RHS,blob_velocity,
                     ! blob_integral_momentum,blob_energy,
                     ! blob_mass_for_velocity, (3 comp)
                     ! volume, 
                     ! centroid_integral, centroid_actual, 
                     ! perim, perim_mat, 
                     ! blob_triple_perim,
                     ! blob_cell_count
                     ! blob_cellvol_count
                     ! blob_mass

                      ! ic+1 is blob_triple_perim index
                     ic=(icolor-1)*num_elements_blobclass+ &
                      3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
                      2*(2*SDIM)+1+ &
                      3+ &  ! blob_mass_for_velocity
                      1+ &  ! volume
                      2*SDIM+ & ! centroid integral, centroid actual
                      1+ & ! perim 
                      nmat ! perim_mat

                     im2=microlayer_substrate_probe
                     if ((im2.ge.1).and.(im2.le.nmat)) then
                      if (base_type.eq.im_source) then
                       im1=im_dest
                      else if (base_type.eq.im_dest) then
                       im1=im_source
                      else
                       print *,"base_type invalid"
                       stop
                      endif
                      ic=ic+(im1-1)*nmat+im2
                      contact_line_perim=blob_array(ic)
                     else 
                      print *,"im2 invalid"
                      stop
                     endif

                     if ((contact_line_perim.ge.zero).and. &
                         (contact_line_perim.le. &
                          max_contact_line_size(im_probe))) then
                      ! do nothing
                     else if (contact_line_perim.gt. &
                              max_contact_line_size(im_probe)) then
                      source_perim_factor= &
                       max_contact_line_size(im_probe)/contact_line_perim
                     else
                      print *,"contact_line_perim invalid"
                      stop
                     endif
                    else if ((local_freezing_model.gt.0).or. &
                             (max_contact_line_size(im_probe).eq.zero).or. &
                             (microlayer_size(im_probe).eq.zero).or. &
                             (macrolayer_size(im_probe).eq.zero).or. &
                             (microlayer_substrate_probe.eq.0)) then
                     ! do nothing
                    else
                     print *,"microlayer parameters invalid"
                     stop
                    endif
       
                   enddo ! iprobe=1,2


                    ! V dt * L * L = volume of material change of phase
                    ! rho_src * V dt L^2 = rho_dst * (V+Vexpand) * dt L^2
                    ! S=V+Vexpand=rho_src V/rho_dst=[k grad T dot n]/(L rho_dst)
                    ! Vexpand=(rho_src/rho_dst-1)V
                    ! note: rho_src V = MDOT

                    ! this call to get_vel_phasechange is not for the purpose
                    ! of estimating the timestep dt.
                   for_estdt=0

#if (STANDALONE==0)
                    ! if local_freezing_model==0 stefan problem
                    !  5, some kind of evaporation model,
                    !  or 1, then
                    !  DTsrc=(Tsrc-TSAT_predict)
                    !  DTdst=(Tdst-TSAT_predict)
                    !  velsrc=ksrc*DTsrc/(LL * dxprobe_src)
                    !  veldst=kdst*DTdst/(LL * dxprobe_dest)
                    !  in: PROB.F90
                   call get_vel_phasechange( &
                     for_estdt, &
                     xI, &
                     ispec, &
                     molar_mass, &
                     species_molar_mass, &
                     local_freezing_model, &
                     local_Tanasawa_or_Schrage_or_Kassemi, &
                     distribute_from_targ, &
                     VEL_correct, & ! vel
                     den_I_interp_SAT(1), & ! source 
                     den_I_interp_SAT(2), & ! dest
                     den_probe(1), & ! source 
                     den_probe(2), & ! dest
                     thermal_k(1), &
                     thermal_k(2), & ! ksrc,kdst
                     T_probe(1), & ! source
                     T_probe(2), & ! dest
                     TSAT_predict, &
                     T_I_interp(1), & !source
                     T_I_interp(2), & !dest
                     LL(ireverse), &
                     source_perim_factor, &
                     dest_perim_factor, &
                     microlayer_substrate_source, &
                     microlayer_angle(im_source), &
                     microlayer_size(im_source), &
                     macrolayer_size(im_source), &
                     microlayer_substrate_dest, &
                     microlayer_angle(im_dest), &
                     microlayer_size(im_dest), &
                     macrolayer_size(im_dest), &
                     dxprobe_target(1), & ! source
                     dxprobe_target(2), & ! dest
                     im_source,im_dest, &
                     prev_time,dt, &
                     fort_alpha(iten+ireverse*nten), &
                     fort_beta(iten+ireverse*nten), &
                     fort_expansion_factor(iten+ireverse*nten), &
                     K_f(ireverse), &
                     Y_I_interp(2), & ! Cmethane_in_hydrate (dest)
                     C_w0, &
                     pres_I_interp(1), & ! PHYDWATER
                     Fsource,Fdest)
#elif (STANDALONE==1)
                   if (local_freezing_model.eq.0) then
                     DTsrc=T_probe(1)-TSAT_predict
                     DTdst=T_probe(2)-TSAT_predict
                     velsrc=thermal_k(1)*DTsrc/(LL(ireverse)*dxprobe_target(1))
                     veldst=thermal_k(2)*DTdst/(LL(ireverse)*dxprobe_target(2))
                   
                     velsum=velsrc+veldst
                     if (velsum.gt.zero) then
                      ! do nothing
                     else if (velsum.le.zero) then
                      velsum=zero
                     else
                      print *,"velsum invalid in stand_alone code"
                      print *,"velsum= ",velsum
                      stop
                     endif 
                     VEL_correct=velsum
                   else
                     print *,"local_freezing_model invalid 4"
                     stop
                   endif
#else
                   print *,"bust compiling ratemasschange"
                   stop
#endif
                  else if (at_interface.eq.0) then
                   ! do nothing
                  else
                   print *,"at_interface inalid"
                   stop
                  endif

                  if (VEL_correct.lt.zero) then
                   VEL_correct=zero
                  endif

                  TSAT_correct=TSAT_predict

                  if (at_interface.eq.1) then

                   if (local_freezing_model.eq.6) then ! Palmore/Desjardins

                     ! type(TSAT_MASS_FRAC_parm_type)
                    Y_interface_min=zero
                    TSAT_Y_PARMS%PROBE_PARMS=>PROBE_PARMS
                    TSAT_Y_PARMS%Tanasawa_or_Schrage_or_Kassemi= &
                           local_Tanasawa_or_Schrage_or_Kassemi
                    TSAT_Y_PARMS%accommodation_coefficient= &
                          accommodation_coefficient(iten+ireverse*nten) 
                    TSAT_Y_PARMS%reference_pressure= &
                          reference_pressure(iten+ireverse*nten) 
                    TSAT_Y_PARMS%YI_min=Y_interface_min
                    TSAT_Y_PARMS%TI_min=TI_min
                    TSAT_Y_PARMS%TI_max=TI_max
                    TSAT_Y_PARMS%universal_gas_constant_R= &
                            R_Palmore_Desjardins
                    TSAT_Y_PARMS%molar_mass_ambient=molar_mass_ambient
                    TSAT_Y_PARMS%molar_mass_vapor=molar_mass_vapor
                    TSAT_Y_PARMS%TSAT_base=local_Tsat(ireverse)
                    TSAT_Y_PARMS%D_MASS=FicksLawD(iprobe_vapor)
                    TSAT_Y_PARMS%den_G=den_I_interp_SAT(iprobe_vapor)
                    TSAT_Y_PARMS%thermal_k=>thermal_k
                    TSAT_Y_PARMS%iprobe_vapor=iprobe_vapor
                    TSAT_Y_PARMS%material_type_evap=>material_type_evap

                    if (hardwire_flag(ireverse).eq.0) then

                         ! Kassemi
                     if (fully_saturated.eq.2) then

                      T_gamma_c=half*(T_gamma_a+T_gamma_b)
                      X_gamma_c=one
                      Y_gamma_c=one

                      call mdot_diff_func(TSAT_Y_PARMS, &
                         Y_gamma_a,T_gamma_a,mdot_diff_a)
                      call mdot_diff_func(TSAT_Y_PARMS, &
                         Y_gamma_b,T_gamma_b,mdot_diff_b)
                      call mdot_diff_func(TSAT_Y_PARMS, &
                         Y_gamma_c,T_gamma_c,mdot_diff_c)

                      if (mdot_diff_a*mdot_diff_b.le.zero) then
                       if (mdot_diff_a*mdot_diff_c.gt.zero) then
                        T_gamma_a=T_gamma_c
                       else if (mdot_diff_a*mdot_diff_c.le.zero) then
                        T_gamma_b=T_gamma_c
                       else
                        print *,"mdot_diff_a or mdot_diff_c invalid"
                        stop
                       endif
                      else
                       print *,"bracketing interval lost Kassemi model"
                       print *,"T_gamma init: a,b ", &
                          T_gamma_a_init,T_gamma_b_init
                       print *,"Y_gamma init: a,b ", &
                          Y_gamma_a_init,Y_gamma_b_init
                       print *,"T_gamma a,b,c ", &
                          T_gamma_a,T_gamma_b,T_gamma_c
                       print *,"Y_gamma a,b,c ", &
                          Y_gamma_a,Y_gamma_b,Y_gamma_c
                       print *,"mdot_diff a,b,c ", &
                          mdot_diff_a,mdot_diff_b,mdot_diff_c
                       print *,"TSAT_iter=",TSAT_iter
                       print *,"Y_interface_min ",Y_interface_min
                       print *,"TI_min ",TI_min
                       print *,"TI_max ",TI_max
                       print *,"molar_mass_ambient ",molar_mass_ambient
                       print *,"molar_mass_vapor ",molar_mass_vapor
                       print *,"iprobe_vapor = ",iprobe_vapor
                       print *,"FicksLawD(iprobe_vapor) ", &
                               FicksLawD(iprobe_vapor)
                       print *,"den_I_interp_sat ", &
                               den_I_interp_SAT(iprobe_vapor)
                       print *," LL(ireverse) ",LL(ireverse)
                       print *," local_Tsat(ireverse) ",local_Tsat(ireverse)
                       print *,"im_source,im_dest ",im_source,im_dest
                       stop
                      endif

                      TSAT_correct=T_gamma_c
                      Y_predict=Y_gamma_c

                       ! Palmore and Desjardins
                     else if ((molar_mass_ambient.gt.zero).and. &
                              (molar_mass_vapor.gt.zero).and. &
                              (R_Palmore_Desjardins.gt.zero).and. &
                              (TSAT_Y_PARMS%den_G.gt.zero).and. &
                              ((fully_saturated.eq.0).or. &
                               (fully_saturated.eq.1))) then

                      if (fully_saturated.eq.1) then
                       Y_interface_min=one
                       Y_predict=one
                       TSAT_correct=local_Tsat(ireverse)

                      else if ((Y_probe(iprobe_vapor).ge. &
                                one-Y_TOLERANCE).and. &
                               (Y_probe(iprobe_vapor).le.one)) then
                       Y_interface_min=one
                       Y_predict=one
                       TSAT_correct=local_Tsat(ireverse)
                      else if (TSAT_Y_PARMS%D_MASS.eq.zero) then
                       Y_interface_min=one
                       Y_predict=one
                       TSAT_correct=local_Tsat(ireverse)
                      else if ((Y_probe(iprobe_vapor).le.one-Y_TOLERANCE).and. &
                               (Y_probe(iprobe_vapor).ge.zero).and. &
                               (TSAT_Y_PARMS%D_MASS.gt.zero)) then

                       if (fully_saturated.eq.0) then
                        T_gamma_c=half*(T_gamma_a+T_gamma_b)
                        call X_from_Tgamma(X_gamma_c,T_gamma_c, &
                         local_Tsat(ireverse), &
                         LL(ireverse),R_Palmore_Desjardins, &
                         molar_mass_vapor) ! WV
                        call massfrac_from_volfrac(X_gamma_c,Y_gamma_c, &
                         molar_mass_ambient,molar_mass_vapor) ! WA,WV

                        call mdot_diff_func(TSAT_Y_PARMS, &
                           Y_gamma_a,T_gamma_a,mdot_diff_a)
                        call mdot_diff_func(TSAT_Y_PARMS, &
                           Y_gamma_b,T_gamma_b,mdot_diff_b)
                        call mdot_diff_func(TSAT_Y_PARMS, &
                           Y_gamma_c,T_gamma_c,mdot_diff_c)

                        if (mdot_diff_a*mdot_diff_b.le.zero) then
                         if (mdot_diff_a*mdot_diff_c.gt.zero) then
                          T_gamma_a=T_gamma_c
                          Y_gamma_a=Y_gamma_c
                         else if (mdot_diff_a*mdot_diff_c.le.zero) then
                          T_gamma_b=T_gamma_c
                          Y_gamma_b=Y_gamma_c 
                         else
                          print *,"mdot_diff_a or mdot_diff_c invalid"
                          stop
                         endif
                        else
                         print *,"bracketing interval lost"
                         print *,"T_gamma init: a,b ", &
                            T_gamma_a_init,T_gamma_b_init
                         print *,"Y_gamma init: a,b ", &
                            Y_gamma_a_init,Y_gamma_b_init
                         print *,"T_gamma a,b,c ", &
                            T_gamma_a,T_gamma_b,T_gamma_c
                         print *,"Y_gamma a,b,c ", &
                            Y_gamma_a,Y_gamma_b,Y_gamma_c
                         print *,"mdot_diff a,b,c ", &
                            mdot_diff_a,mdot_diff_b,mdot_diff_c
                         print *,"TSAT_iter=",TSAT_iter
                         print *,"Y_interface_min ",Y_interface_min
                         print *,"TI_min ",TI_min
                         print *,"TI_max ",TI_max
                         print *,"molar_mass_ambient ",molar_mass_ambient
                         print *,"molar_mass_vapor ",molar_mass_vapor
                         print *,"iprobe_vapor = ",iprobe_vapor
                         print *,"FicksLawD(iprobe_vapor) ", &
                                 FicksLawD(iprobe_vapor)
                         print *,"den_I_interp_sat ", &
                                 den_I_interp_SAT(iprobe_vapor)
                         print *," LL(ireverse) ",LL(ireverse)
                         print *," local_Tsat(ireverse) ",local_Tsat(ireverse)
                         print *,"im_source,im_dest ",im_source,im_dest
                         stop
                        endif

                        TSAT_correct=T_gamma_c
                        Y_predict=Y_gamma_c

                       else 
                        print *,"expecting fully_saturated==0"
                        stop
                       endif

                      else
                       print *,"Y_probe invalid"
                       stop
                      endif

                     else
                      print *,"molar masses, den_G, fully_saturated, or R bad"
                      print *,"molar_mass_ambient=",molar_mass_ambient
                      print *,"molar_mass_vapor=",molar_mass_vapor
                      print *,"R_Palmore_Desjardins=",R_Palmore_Desjardins
                      print *,"TSAT_Y_PARMS%den_G=",TSAT_Y_PARMS%den_G
                      print *,"fully_saturated=",fully_saturated
                      print *,"local_freezing_model=",local_freezing_model
                      stop
                     endif

                    else if (hardwire_flag(ireverse).eq.1) then
                     ! do nothing
                    else
                     print *,"hardwire_flag(ireverse) invalid"
                     stop
                    endif
 
                    call mdot_from_T_probe( &
                     TSAT_Y_PARMS, &
                     TSAT_correct,Y_predict, &
                     mdotT_debug)

                    call mdot_from_Y_probe( &
                     TSAT_Y_PARMS, &
                     Y_predict,TSAT_correct, &
                     mdotY_top_debug,mdotY_bot_debug,mdotY_debug)

                   else if (local_freezing_model.ge.0) then

                    mdotT_debug=zero
                    mdotY_top_debug=zero
                    mdotY_bot_debug=zero
                    mdotY_debug=zero

                   else
                    print *,"local_freezing_model invalid 7"
                    stop
                   endif

                  else if (at_interface.eq.0) then

                   mdotT_debug=zero
                   mdotY_top_debug=zero
                   mdotY_bot_debug=zero
                   mdotY_debug=zero

                  else
                   print *,"at_interface invalid FORT_RATEMASSCHANGE (2) "
                   print *,"at_interface=",at_interface
                   stop
                  endif

                  if (hardwire_flag(ireverse).eq.0) then

                   if ((interface_resolved.eq.1).or. &
                       (TSAT_iter.eq.0)) then
                    delta_Tsat= &
                     saturation_temp_vel(iten+ireverse*nten)* &
                     (VEL_correct-VEL_predict)
                   else if ((interface_resolved.eq.0).and. &
                            (TSAT_iter.gt.0)) then
                    delta_Tsat=zero
                   else
                    print *,"interface_resolved or TSAT_iter invalid"
                    stop
                   endif

                    ! if interface_resolved==1 or first iteration,
                    !  Tgamma^{k+1}=Tgamma^{k}- 
                    !    eps2 * (V^{k+1}(Tgamma^{k})-
                    !            V^{k}(Tgamma^{k-1}) ) 
                    !  V^{0}=0
                   TSAT_correct=TSAT_correct-delta_Tsat

                   if (TSAT_correct.lt.TI_min) then
                    TSAT_correct=TI_min
                   endif
                   if (TSAT_correct.gt.TI_max) then
                    TSAT_correct=TI_max
                   endif

                  else if (hardwire_flag(ireverse).eq.1) then
                   TSAT_correct=local_hardwire_T(ireverse)
                  else
                   print *,"hardwire_flag(ireverse) invalid"
                   stop
                  endif
             
                  VEL_predict=VEL_correct
  
                  TSAT_ERR=abs(TSAT_correct-TSAT_predict)

                  TSAT_predict=TSAT_correct

                  if (TSAT_iter.eq.0) then
                    TSAT_INIT_ERR=TSAT_ERR
                  endif

                  if (1.eq.0) then
                   if (j.eq.8) then
                    print *,"PD DEBUG prev_time,dt ",prev_time,dt
                    print *,"i,j,k,x ",i,j,k,xsten(0,1),xsten(0,2), &
                      xsten(0,SDIM)
                    print *,"TSAT_iter,TSAT_correct,TSAT_predict ", &
                     TSAT_iter,TSAT_correct,TSAT_predict
                    print *,"Y_predict ",Y_predict
                    print *,"VEL_correct ",VEL_correct
                    print *,"mdotT_debug ",mdotT_debug
                    print *,"mdotY_debug ",mdotY_debug
                    print *,"mdotY_top_debug ",mdotY_top_debug
                    print *,"mdotY_bot_debug ",mdotY_bot_debug
                   endif
                  endif

                  TSAT_converge=0
                  if (TSAT_ERR.eq.zero) then
                    TSAT_converge=1
                  endif
                   ! TSAT_iter starts at 0
                  TSAT_iter=TSAT_iter+1
                  if (TSAT_iter.gt.EVAPORATION_iter_max) then
                    TSAT_converge=1
                  endif
                  if (TSAT_iter.gt.1) then
                    if (TSAT_err.lt.EVAPORATION_TOL*TSAT_INIT_ERR) then
                     TSAT_converge=1
                    endif
                  endif
                  if (at_interface.eq.0) then
                    TSAT_converge=1
                  endif

                  if (at_interface.eq.1) then
                   if (DEBUG_EVAPORATION.eq.1) then
                    print *,"DEBUG_EVAPORATION STATEMENT 3"
                    print *,"i,j,k,TSAT_iter,TSAT_ERR ", &
                     i,j,k,TSAT_iter,TSAT_ERR
                    print *,"TSAT_correct ",TSAT_correct
                    print *,"Y_predict ",Y_predict
                   endif
                  endif

                 enddo ! do while (TSAT_converge.eq.0)

                 if (at_interface.eq.1) then

                  local_Tsat(ireverse)=TSAT_correct
                  vel_phasechange(ireverse)=VEL_correct

                    ! source
                  temp_target_probe_history(iten+ireverse*nten,1)=T_probe(1)
                  dxprobe_target_history(iten+ireverse*nten,1)=dxprobe_target(1)
                    ! dest
                  temp_target_probe_history(iten+ireverse*nten,2)=T_probe(2)
                  dxprobe_target_history(iten+ireverse*nten,2)=dxprobe_target(2)

                  if (debugrate.eq.1) then
                   print *,"i,j,k,ireverse,vel_phasechange ", &
                    i,j,k,ireverse,vel_phasechange(ireverse)
                  endif
                  if (1.eq.0) then
                   if (local_freezing_model.eq.4) then
                    print *,"Tanasawa or Schrage"
                    print *,"i,j,k,ireverse,vel_phasechange ", &
                     i,j,k,ireverse,vel_phasechange(ireverse)
                    print *,"im_source,im_dest ",im_source,im_dest
                    print *,"local_Tsat(ireverse) ",local_Tsat(ireverse)
                    print *,"T_probe(1),T_I_interp(1) ",T_probe(1),T_I_interp(1)
                    print *,"T_probe(2),T_I_interp(2) ",T_probe(2),T_I_interp(2)
                   endif
                  endif

                  ! Li-Shi Luo - invite him to talk with FSU?
                  ! LATTICE BOLTZMANN TREATMENT:
                  ! max velocity cannot exceed dxmin/(2 dt)
                  if (dt.gt.zero) then
                   if (dxmin.gt.zero) then
                    if (vel_phasechange(ireverse).ge.dxmin/(two*dt)) then
                     if (stefan_flag.eq.1) then
                      vel_phasechange(ireverse)=dxmin/(two*dt)
                     endif
                    else if (vel_phasechange(ireverse).ge.zero) then
                     ! do nothing
                    else
                     print *,"vel_phasechange bust"
                     stop
                    endif
                   else
                    print *,"dxmin invalid"
                    stop
                   endif
                  else
                   print *,"dt invalid"
                   stop
                  endif

                  valid_phase_change(ireverse)=1
              
                  if (debugrate.eq.1) then
                   print *,"i,j,k,im_source,im_dest ",i,j,k, &
                    im_source,im_dest 
                   print *,"dt,vel_phasechange(ireverse) ", &
                    dt,vel_phasechange(ireverse)
                   print *,"LL,dxmin ",LL(ireverse),dxmin
                   print *,"dxprobe_target(1)=",dxprobe_target(1)
                   print *,"dxprobe_target(2)=",dxprobe_target(2)
                   print *,"thermal_k ",thermal_k(1),thermal_k(2)
                   print *,"local_Tsat(ireverse) ",local_Tsat(ireverse)
                   print *,"T_Probe(1),T_probe(2) ",T_Probe(1),T_probe(2)
                   print *,"den_I_interp(1) ",den_I_interp(1)
                   print *,"den_I_interp_SAT(1) ",den_I_interp_SAT(1)
                   print *,"den_I_interp_SAT(2) ",den_I_interp_SAT(2)
                   print *,"LSINTsrc,LSINTdst ",LSINT(im_source),LSINT(im_dest)
                   print *,"nrmCP ",nrmCP(1),nrmCP(2),nrmCP(SDIM)
                   print *,"nrmFD ",nrmFD(1),nrmFD(2),nrmFD(SDIM)
                   print *,"nrmPROBE ",nrmPROBE(1),nrmPROBE(2),nrmPROBE(SDIM)
                   print *,"im_dest= ",im_dest
                  endif
    
                 else if (at_interface.eq.0) then
                  ! do nothing
                 else
                  print *,"at_interface invalid in FORT_RATEMASSCHANGE (3)"
                  print *,"at_interface=",at_interface
                  stop
                 endif

                else if (found_path.eq.0) then
                 ! do nothing
                else
                 print *,"found_path invalid: ",found_path
                 stop
                endif

               else if ((abs(LShere(im_source)).gt.two*dxmaxLS).or. &
                        (abs(LShere(im_dest)).gt.two*dxmaxLS)) then
                ! do nothing
               else
                print *,"LShere bust:"
                print *,"LShere(im_source) ",LShere(im_source)
                print *,"LShere(im_dest) ",LShere(im_dest)
                stop
               endif

              else
               print *,"LL(ireverse) bust: ",LL(ireverse)
               stop
              endif

             else if (LL(ireverse).eq.zero) then
              ! do nothing
             else
              print *,"LL bust"
              stop
             endif

            enddo ! ireverse=0,...,1
             
            ireverse=-1 
            if ((valid_phase_change(0).eq.0).and. &
                (valid_phase_change(1).eq.0)) then
             ! do nothing
            else if ((valid_phase_change(0).eq.1).and. &
                     (valid_phase_change(1).eq.0)) then
             ireverse=0
            else if ((valid_phase_change(0).eq.0).and. &
                     (valid_phase_change(1).eq.1)) then
             ireverse=1
            else if ((valid_phase_change(0).eq.1).and. &
                     (valid_phase_change(1).eq.1)) then
             if (vel_phasechange(0).gt.vel_phasechange(1)) then
              ireverse=0
             else if (vel_phasechange(1).gt.vel_phasechange(0)) then
              ireverse=1
             else if (vel_phasechange(1).eq.vel_phasechange(0)) then
              ireverse=-1
             else
              print *,"vel_phasechange bust"
              stop
             endif
            else
             print *,"valid_phase_change bust"
             stop
            endif

            if ((ireverse.eq.0).or.(ireverse.eq.1)) then

             if (ireverse.eq.0) then
              im_source=im
              im_dest=im_opp
             else if (ireverse.eq.1) then
              im_source=im_opp
              im_dest=im
             else
              print *,"ireverse invalid"
              stop
             endif

             if (im_dest.lt.im_source) then
              LSSIGN=one
             else if (im_dest.gt.im_source) then
              LSSIGN=-one
             else
              print *,"im_dest<>im_source not satisfied"
              stop
             endif

             if (stefan_flag.eq.1) then
              LSnew(D_DECL(i,j,k),im_source)= &
               LSnew(D_DECL(i,j,k),im_source)-dt*vel_phasechange(ireverse)
              LSnew(D_DECL(i,j,k),im_dest)= &
               LSnew(D_DECL(i,j,k),im_dest)+dt*vel_phasechange(ireverse)
             else if (stefan_flag.eq.0) then
              ! do nothing
             else
              print *,"stefan_flag invalid in rate mass change"
              stop
             endif

             do dir=1,SDIM
              if (LShere(im_dest).ge.zero) then
               SIGNVEL=one
               nrmFD(dir)= &
                  LS_slopes_FD(D_DECL(i,j,k),(im_source-1)*SDIM+dir)
               nrmCP(dir)= &
                  LS(D_DECL(i,j,k),nmat+(im_source-1)*SDIM+dir)
              else if (LShere(im_source).ge.zero) then
               SIGNVEL=-one
               nrmFD(dir)= &
                  LS_slopes_FD(D_DECL(i,j,k),(im_dest-1)*SDIM+dir)
               nrmCP(dir)= &
                  LS(D_DECL(i,j,k),nmat+(im_dest-1)*SDIM+dir)
              else
               print *,"LShere bust"
               stop
              endif
              nrmPROBE(dir)=nrmCP(dir)
             enddo ! dir=1..sdim

              ! units of k (thermal conductivity): watts/(m kelvin)
              ! watts=kg m^2/s^3
              ! units of k=kg m/(s^3 kelvin)
              ! units of L=Joule/kg
              ! joule=kg m^2/s^2
              ! units of L=m^2/s^2
              ! units of k/L=kg (1/m)(1/s)(1/kelvin)
              ! units of (1/rho_source)[k grad T]/L =
              ! (m^3/kg)(kelvin/m)kg(1/m)(1/s)(1/kelvin)=
              ! m/s 
             if ((im_dest.ge.1).and.(im_dest.le.nmat)) then
              if (ireverse.eq.0) then
               burnvel(D_DECL(i,j,k),iten)=one
               Tsatfab(D_DECL(i,j,k),iten)=one
              else if (ireverse.eq.1) then
               burnvel(D_DECL(i,j,k),iten)=-one
               Tsatfab(D_DECL(i,j,k),iten)=-one
              else
               print *,"ireverse invalid"
               stop
              endif

              if (local_Tsat(ireverse).gt.zero) then
               Tsatfab(D_DECL(i,j,k),nten+ncomp_per_tsat*(iten-1)+1)= &
                local_Tsat(ireverse)
               if ((Y_predict.ge.zero).and.(Y_predict.le.one)) then
                Tsatfab(D_DECL(i,j,k),nten+ncomp_per_tsat*(iten-1)+2)= &
                  Y_predict  ! default mass fraction=1 (saturated)
               else
                print *,"Y_predict invalid"
                stop
               endif
              else
               print *,"local_Tsat(ireverse) should be positive"
               stop
              endif

              if ((vel_phasechange(ireverse).ge.zero).or. &
                  (vel_phasechange(ireverse).le.zero)) then
               do dir=1,ncomp_per_burning
                burnvel(D_DECL(i,j,k),nten+(iten-1)*ncomp_per_burning+dir)= &
                 SIGNVEL*nrmPROBE(dir)*vel_phasechange(ireverse)
               enddo
              else
               print *,"vel_phasechange(ireverse) cannot be NaN"
               stop
              endif
             else
              print *,"im_dest bust"
              stop
             endif

            else if (ireverse.eq.-1) then
             ! do nothing
            else
             print *,"ireverse invalid"
             stop
            endif

           enddo ! im_opp=im+1,nmat
          enddo ! im=1,nmat-1

          velmag_sum=zero
          do im=1,nmat-1
           do im_opp=im+1,nmat
            call get_iten(im,im_opp,iten,nmat)
            burnflag=NINT(burnvel(D_DECL(i,j,k),iten))
            if (burnflag.eq.0) then
             ! do nothing
            else if ((burnflag.eq.1).or. &
                     (burnflag.eq.-1)) then
             local_velmag=zero
             do dir=1,ncomp_per_burning
              local_velmag=local_velmag+ &
               burnvel(D_DECL(i,j,k),nten+(iten-1)*ncomp_per_burning+dir)**2
             enddo
             local_velmag=sqrt(local_velmag)
             velmag_sum=velmag_sum+local_velmag
            else
             print *,"burnflag invalid"
             stop
            endif 
           enddo !im_opp=im+1..nmat
          enddo !im=1..nmat-1

           ! factor of 2 in order to guarantee that characteristics do not
           ! collide.
           ! LATTICE BOLTZMANN TREATMENT:
           ! max velocity cannot exceed dxmin/(2 dt)
           ! (see above, search for "Li-Shi Luo")
          if ((two*velmag_sum*dt.le.(one+VOFTOL)*dxmin).and. &
              (velmag_sum.ge.zero)) then
           ! do nothing
          else if (two*velmag_sum*dt.gt.dxmin) then
           if (stefan_flag.eq.1) then
            print *,"phase change velocity exceeds cfl limits"
            print *,"velmag_sum ",velmag_sum
            print *,"dxmin ",dxmin
            print *,"dt ",dt
            print *,"velmag_sum x dt ",velmag_sum*dt
            print *,"i,j,k ",i,j,k
            print *,"im_primary ",im_primary
            do imls=1,nmat
             print *,"imls,LShere ",imls,LShere(imls)
            enddo
            print *,"xsten(0,1-3) ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
            do im=1,nmat-1
             do im_opp=im+1,nmat
              call get_iten(im,im_opp,iten,nmat)
              do ireverse=0,1
               print *,"im,im_opp,ireverse,latent_heat ",im,im_opp,ireverse, &
                       latent_heat(iten+ireverse*nten)
               print *,"im,im_opp,ireverse,saturation_temp ", &
                       im,im_opp,ireverse, &
                       saturation_temp(iten+ireverse*nten)
              enddo
             enddo
            enddo
            do im=1,nmat
             print *,"im,T ",im,EOS(D_DECL(i,j,k),(im-1)*num_state_material+2)
            enddo 
            do iten=1,nten
             do ireverse=0,1
              print *,"iten,ireverse,temp_probe(1) ", &
                iten,ireverse,temp_target_probe_history(iten+ireverse*nten,1)
              print *,"iten,ireverse,temp_probe(2) ", &
                iten,ireverse,temp_target_probe_history(iten+ireverse*nten,2)
              print *,"iten,ireverse,dxprobe(1) ", &
                iten,ireverse,dxprobe_target_history(iten+ireverse*nten,1)
              print *,"iten,ireverse,dxprobe(2) ", &
                iten,ireverse,dxprobe_target_history(iten+ireverse*nten,2)
             enddo
            enddo
            stop
           else if (stefan_flag.eq.0) then
            ! do nothing
           else
            print *,"stefan_flag invalid"
            stop
           endif
          else
           print *,"velmag_sum invalid"
           stop
          endif

         else if (is_rigid(nmat,im_primary).eq.1) then

          ! do nothing

         else 
          print *,"is_rigid invalid"
          stop
         endif

        else if (nucleation_flag.eq.1) then
         ! see RatePhaseChange in PROB.F90
         ! LEVELSET FUNCTION AT CELL CENTERS.
         do im=1,nmat
          LShere(im)=LSnew(D_DECL(i,j,k),im)
         enddo
         call get_primary_material(LShere,nmat,im_primary)
         if (is_rigid(nmat,im_primary).eq.0) then
          do im=1,nmat-1
           do im_opp=im+1,nmat
            if ((im.gt.nmat).or.(im_opp.gt.nmat)) then
             print *,"im or im_opp bust 9"
             stop
            endif
            call get_iten(im,im_opp,iten,nmat)
            do ireverse=0,1
             LL(ireverse)=latent_heat(iten+ireverse*nten)
             local_freezing_model=freezing_model(iten+ireverse*nten)
             local_Tsat(ireverse)=saturation_temp(iten+ireverse*nten)

             if (LL(ireverse).ne.zero) then
            
              if (ireverse.eq.0) then
               im_source=im
               im_dest=im_opp
              else if (ireverse.eq.1) then
               im_source=im_opp
               im_dest=im
              else
               print *,"ireverse invalid"
               stop
              endif

              if ((is_rigid(nmat,im).eq.1).or. &
                  (is_rigid(nmat,im_opp).eq.1)) then

               ! do nothing

              else if ((is_rigid(nmat,im).eq.0).and. &
                       (is_rigid(nmat,im_opp).eq.0)) then
               if (im_primary.eq.im_source) then
                create_in%LL=LL(ireverse)
                create_in%local_freezing_model=local_freezing_model
                create_in%local_TSAT=local_Tsat(ireverse)
                create_in%im_source=im_source
                create_in%im_dest=im_dest
#if (STANDALONE==0)
                call get_vel_phasechange_NUCLEATE( &
                 create_in,create_inout)
#elif (STANDALONE==1)
                print *,"should not call get_vel_phasechange_NUCLEATE for"
                print *,"stand alone version"
                stop
#else
                print *,"bust compiling ratemasschange"
                stop
#endif
               else if ((im_primary.ge.1).and.(im_primary.le.nmat)) then
                ! do nothing
               else
                print *,"im_primary invalid"
                stop
               endif
              else
               print *,"is_rigid(nmat,im or im_opp) invalid"
               stop
              endif
             else if (LL(ireverse).eq.zero) then
              ! do nothing
             else
              print *,"LL invalid"
              stop
             endif
            enddo ! ireverse=0,1
           enddo ! im_opp
          enddo ! im

         else if (is_rigid(nmat,im_primary).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(nmat,im_primary) invalid"
          stop
         endif

        else
         print *,"nucleation_flag invalid"
         stop
        endif

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine FORT_RATEMASSCHANGE

#if (STANDALONE==1)
      end module mass_transfer_cpp_module
#endif

#undef STANDALONE


