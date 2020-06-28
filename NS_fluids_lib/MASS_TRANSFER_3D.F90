#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 0

#define DEBUG_TRIPLE 0
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
       REAL_T, pointer, dimension(:) :: density_floor_expansion
       REAL_T, pointer, dimension(:) :: density_ceiling_expansion
      end type probe_parm_type

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
        if (abs(T_test).ge.1.0D+50) then
         print *,"T_test bust"
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
       print *,"Tsat out of range"
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
       VOF_pos_probe_counter)
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
       print *,"Tsat out of range"
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
       print *,"mag invalid"
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
       den_I_interp,T_I_interp,Y_I_interp, &
       pres_I_interp, &
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
      REAL_T, intent(out) :: T_I_interp(2)
      REAL_T, intent(out) :: Y_I_interp(2)
      REAL_T, intent(out) :: pres_I_interp(2)
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

      do iprobe=1,2

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

       mtype=fort_material_type(im_target_probe(iprobe))
       if ((mtype.ge.0).and. &
           (mtype.le.fort_max_num_eos)) then
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
         PROBE_PARMS%EOS, &
         DIMS(PROBE_PARMS%EOS), &
         PROBE_PARMS%recon, &
         DIMS(PROBE_PARMS%recon), &
         den_I_interp(iprobe))

        if (den_I_interp(iprobe).lt. &
            PROBE_PARMS%density_floor_expansion(im_target_probe(iprobe))) then
         den_I_interp(iprobe)= &
          PROBE_PARMS%density_floor_expansion(im_target_probe(iprobe))
        endif
        if (den_I_interp(iprobe).gt. &
            PROBE_PARMS%density_ceiling_expansion(im_target_probe(iprobe))) then
         den_I_interp(iprobe)= &
          PROBE_PARMS%density_ceiling_expansion(im_target_probe(iprobe))
        endif
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

        if ((Y_probe(iprobe).lt.zero).or. &
            (Y_probe(iprobe).gt.one)) then
         print *,"Y_probe out of bounds"
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
          VOF_pos_probe_counter)

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
           dummy_VOF_pos_probe_counter)

          if ((Y_probe(iprobe).ge.zero).and. &
              (Y_probe(iprobe).le.one)) then
           ! do nothing
          else
           print *,"Y_probe out of bounds"
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
       else if ((PROBE_PARMS%local_freezing_model.eq.1).or. &
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
       isweep, &
       solvability_projection, &
       ngrow_expansion, &
       level,finest_level, &
       normal_probe_size, &
       nmat, &
       nten, &
       nden, &
       nstate, &
       ntsat, &
       density_floor_expansion, &
       density_ceiling_expansion, &
       latent_heat, &
       saturation_temp, &
       freezing_model, &
       mass_fraction_id, &
       species_evaporation_density, &
       distribute_from_target, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       vofbc, &
       xlo,dx, &
       dt, &
       delta_mass, &
       DVOF, &
       maskcov,DIMS(maskcov), &
       deltaVOF,DIMS(deltaVOF), &
       nodevel,DIMS(nodevel), &
       JUMPFAB,DIMS(JUMPFAB), &
       TSATFAB,DIMS(TSATFAB), &
       LSold,DIMS(LSold), &
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

      IMPLICIT NONE

      INTEGER_T, intent(in) :: isweep,solvability_projection,tid
      INTEGER_T, intent(in) :: level,finest_level,ngrow_expansion
      INTEGER_T, intent(in) :: normal_probe_size
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nden
      INTEGER_T, intent(in) :: nstate
      INTEGER_T, intent(in) :: ntsat
      REAL_T, intent(in) :: density_floor_expansion(nmat)
      REAL_T, intent(in) :: density_ceiling_expansion(nmat)
      REAL_T, intent(in) :: latent_heat(2*nten)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: mass_fraction_id(2*nten)
      REAL_T, intent(in) :: species_evaporation_density(num_species_var+1)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: vofbc(SDIM,2)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: dt
      REAL_T, intent(inout) :: delta_mass(2*nmat)
      REAL_T, intent(inout) :: DVOF(nmat)
      REAL_T :: DVOF_FACT(nmat)
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(deltaVOF)
      INTEGER_T, intent(in) :: DIMDEC(nodevel)
      INTEGER_T, intent(in) :: DIMDEC(JUMPFAB)
      INTEGER_T, intent(in) :: DIMDEC(TSATFAB)
      INTEGER_T, intent(in) :: DIMDEC(LSold)
      INTEGER_T, intent(in) :: DIMDEC(LSnew)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(snew)
      INTEGER_T, intent(in) :: DIMDEC(EOS)
      INTEGER_T, intent(in) :: DIMDEC(swept)

      REAL_T, intent(in) :: maskcov(DIMV(maskcov))

      REAL_T, intent(inout) :: deltaVOF(DIMV(deltaVOF),nmat)

      REAL_T, intent(in) :: nodevel(DIMV(nodevel),2*nten*SDIM)
      REAL_T, intent(out) :: JUMPFAB(DIMV(JUMPFAB),2*nten)
      REAL_T, intent(out) :: TSATFAB(DIMV(TSATFAB),2*nten)
      REAL_T, intent(in) :: LSold(DIMV(LSold),nmat*(1+SDIM))
      REAL_T, intent(out) :: LSnew(DIMV(LSnew),nmat)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(out) :: snew(DIMV(snew),nstate)
      REAL_T, intent(in) :: EOS(DIMV(EOS),nden)
      REAL_T, intent(out) :: swept(DIMV(swept))

      INTEGER_T i,j,k,dir
      INTEGER_T i1,j1,k1
      INTEGER_T im,im_opp,ireverse,iten
      INTEGER_T im_source
      INTEGER_T im_dest
      INTEGER_T im_dest_crit
      INTEGER_T im_source_crit
      INTEGER_T iten_crit
      INTEGER_T ireverse_crit
      REAL_T max_velnode
      REAL_T velnode_test
      INTEGER_T nten_test
      INTEGER_T vcompsrc_snew,vcompdst_snew
      INTEGER_T dcompsrc,dcompdst
      INTEGER_T dcompdst_snew,tcompdst_snew
      REAL_T densrc,dendst
      REAL_T densrc_restrict,dendst_restrict
      REAL_T denratio_factor
      REAL_T oldvfrac(nmat)
      REAL_T newvfrac(nmat)
      REAL_T dF,dFdst,dFsrc
      REAL_T jump_strength
      REAL_T xsten(-3:3,SDIM)

      REAL_T volgrid
      REAL_T cengrid(SDIM)
      REAL_T new_centroid(nmat,SDIM)
      REAL_T EBVOFTOL
      REAL_T SWEPTFACTOR
      REAL_T LL
      REAL_T Tsat_default
      INTEGER_T Tsat_flag
      REAL_T ksource,kdest
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
      INTEGER_T local_freezing_model
      INTEGER_T mass_frac_id
      INTEGER_T distribute_from_targ
      INTEGER_T debugrate
      REAL_T F_STEN(nmat)
      REAL_T unsplit_snew(nmat*ngeom_raw)
      REAL_T unsplit_density(nmat)
      REAL_T unsplit_temperature(nmat)
      REAL_T unsplit_lsnew(nmat)
      REAL_T oldLS_point(nmat)
      REAL_T dxmax,dxmaxLS
      INTEGER_T recon_ncomp
      INTEGER_T scomp
 
      INTEGER_T tcomp,dencomp
      REAL_T density_data(nmat)
      REAL_T temperature_data(nmat)
      REAL_T density_mat(nmat)
      REAL_T temperature_mat(nmat)
      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T volmat(nmat)
      REAL_T lsmat(nmat)
      REAL_T lsdata(nmat)
      REAL_T cenmat(SDIM,nmat)
      REAL_T multi_volume(nmat)
      REAL_T multi_area(nmat)
      REAL_T multi_cen(SDIM,nmat)

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
      REAL_T volcell,tempvfrac
      REAL_T cencell(SDIM)
      REAL_T tempcen(SDIM)
      INTEGER_T do_unsplit_advection
      INTEGER_T interface_near(2*nten)
      INTEGER_T local_mask
      INTEGER_T im_primary_new
      INTEGER_T im_primary_old
      INTEGER_T away_from_interface
      REAL_T solid_vof_new,solid_vof_old
      INTEGER_T mtype
      REAL_T Yfrac_vapor_to_gas
      REAL_T vfrac_vapor_to_gas
      REAL_T new_vfrac_vapor
      REAL_T Yfrac_old
      REAL_T density_new,density_old
      REAL_T density_dest,density_source
      REAL_T evap_den
      REAL_T local_cv_or_cp
      INTEGER_T speccompsrc,speccompdst
      INTEGER_T ncomp_per_tsat

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

      do im=1,nmat
       if ((density_floor_expansion(im).gt.zero).and. &
           (density_floor_expansion(im).le.fort_denconst(im))) then
        ! do nothing
       else
        print *,"density_floor_expansion invalid"
        stop
       endif
       if ((density_ceiling_expansion(im).gt.zero).and. &
           (density_ceiling_expansion(im).ge.fort_denconst(im))) then
        ! do nothing
       else
        print *,"density_ceiling_expansion invalid"
        stop
       endif
      enddo ! im=1..nmat

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

       ! DVOF is accumulated material to change phase.
       ! for heat pipe: we must have sum_m (sign L_m) DVOF_m = 0
      do im=1,nmat
       if (DVOF(im).lt.zero) then
        print *,"DVOF invalid"
        stop
       endif
       DVOF_FACT(im)=one
      enddo ! im

      do im=1,nmat-1
       do im_opp=im+1,nmat
        do ireverse=0,1
         call get_iten(im,im_opp,iten,nmat)
         local_freezing_model=freezing_model(iten+ireverse*nten)
         distribute_from_targ=distribute_from_target(iten+ireverse*nten)
         LL=latent_heat(iten+ireverse*nten)
         if ((local_freezing_model.eq.0).or. &
             (local_freezing_model.eq.1).or. &
             (local_freezing_model.eq.2)) then
          ! do nothing
         else if ((local_freezing_model.eq.4).or. & ! Tanasawa or Schrage
                  (local_freezing_model.eq.5).or. & ! Stefan model evap/cond.
                  (local_freezing_model.eq.6).or. & ! Palmore/Desjardins
                  (local_freezing_model.eq.7)) then ! Cavitation
          mass_frac_id=mass_fraction_id(iten+ireverse*nten)
          if ((mass_frac_id.ge.1).and.(mass_frac_id.le.num_species_var)) then
           evap_den=species_evaporation_density(mass_frac_id)
           if (evap_den.gt.zero) then
            ! do nothing
           else
            print *,"evap_den invalid"
            stop
           endif
          else
           print *,"mass_frac_id invalid"
           stop
          endif
           ! require V_evaporate = V_condense (Tanasawa model or Schrage)
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

            if (isweep.eq.0) then
             ! do nothing
            else if (isweep.eq.1) then

             if (solvability_projection.eq.0) then
              ! do nothing
             else if (solvability_projection.eq.1) then
              if (DVOF(im_dest).lt.DVOF(im_source)) then
               DVOF_FACT(im_source)=DVOF(im_dest)/DVOF(im_source)
              else if (DVOF(im_source).lt.DVOF(im_dest)) then
               DVOF_FACT(im_dest)=DVOF(im_source)/DVOF(im_dest)
              endif
             else
              print *,"solvability_projection invalid"
              stop
             endif

            else
             print *,"isweep invalid"
             stop
            endif

           endif ! LL <> 0
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
       DIMS(deltaVOF),0,-1,1256)

      call checkbound(fablo,fabhi, &
       DIMS(JUMPFAB),ngrow_expansion,-1,1256)
      call checkbound(fablo,fabhi, &
       DIMS(TSATFAB),ngrow_expansion,-1,1256)
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

        do im=1,2*nten
         JUMPFAB(D_DECL(i,j,k),im)=zero
        enddo

        call gridsten_level(xsten,i,j,k,level,nhalf)
        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volgrid,cengrid,SDIM)

        do_unsplit_advection=0
        do im=1,2*nten
         interface_near(im)=0
        enddo

        do im=1,nmat
         lsmat(im)=LSold(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(lsmat,nmat,im)

        if (is_rigid(nmat,im).eq.0) then

         do im=1,nmat
          F_STEN(im)=zero

          vofcomp_recon=(im-1)*ngeom_recon+1

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

             if ((F_STEN(im_source).gt.zero).and. &
                 (F_STEN(im_dest).gt.zero)) then
              interface_near(iten+ireverse*nten)=1
              do_unsplit_advection=1
             else if ((F_STEN(im_source).eq.zero).or. &
                      (F_STEN(im_dest).eq.zero)) then
              ! do nothing
             else
              print *,"F_STEN invalid"
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

           ! sweep==0: advection: unsplit_lsnew, unsplit_snew, dF, deltaVOF
           !                      unsplit_density,unsplit_temperature
           ! unsplit advection: backward tracing of characteristics, but
           !  evaluate volumes, centroids, level set functions, and
           !  temperatures in the target cell.

          if (isweep.eq.0) then

           symmetry_flag=0 
           call get_ntetbox(ntetbox,symmetry_flag,SDIM)

           call gridsten_level(u_xsten_updatecell,i,j,k,level,nhalf)

           absolute_voltotal=zero
           voltotal=zero
           do u_im=1,nmat
            density_mat(u_im)=zero
            temperature_mat(u_im)=zero
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
               density_data(u_im)=EOS(D_DECL(i+igrid,j+jgrid,k+kgrid),dencomp)
               temperature_data(u_im)= &
                 EOS(D_DECL(i+igrid,j+jgrid,k+kgrid),tcomp)
          
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
                 tessellate, &
                 bfact,dx, &
                 u_xsten_departmap,nhalf0, &
                 mofdata, &
                 u_xsten_updatecell,nhalf, &
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

           if (interface_near(iten+ireverse*nten).ne.1) then
            print *,"interface_near invalid"
            stop
           endif

           Tsat_default=saturation_temp(iten+ireverse*nten)
           Tsat_flag=NINT(TSATFAB(D_DECL(i,j,k),iten))
           if (ireverse.eq.0) then
             ! do nothing
           else if (ireverse.eq.1) then
             Tsat_flag=-Tsat_flag
           else
             print *,"ireverse invalid"
             stop
           endif
            
           if ((Tsat_flag.eq.1).or.(Tsat_flag.eq.2)) then
             Tsat_default=TSATFAB(D_DECL(i,j,k), &
              nten+(iten-1)*ncomp_per_tsat+1)
           else if ((Tsat_flag.eq.-1).or. &
                    (Tsat_flag.eq.-2)) then
             ! do nothing
           else if (Tsat_flag.eq.0) then
             ! do nothing
           else
             print *,"Tsat_flag invalid"
             stop
           endif

           LL=latent_heat(iten+ireverse*nten)
           local_freezing_model=freezing_model(iten+ireverse*nten)
           distribute_from_targ=distribute_from_target(iten+ireverse*nten)
           mass_frac_id=mass_fraction_id(iten+ireverse*nten)

           do u_imaterial=1,nmat

            vofcomp_raw=(u_imaterial-1)*ngeom_raw+1

            tempvfrac=volmat(u_imaterial)/voltotal
            if ((tempvfrac.ge.EBVOFTOL).and.(tempvfrac.le.1.1)) then
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
            else if (tempvfrac.eq.zero) then
             unsplit_density(u_imaterial)=fort_denconst(u_imaterial)
             unsplit_temperature(u_imaterial)=Tsat_default
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
             new_centroid(u_imaterial,dir)= &
               unsplit_snew(vofcomp_raw+dir)+cengrid(dir)
            enddo
            oldLS_point(u_imaterial)=LSold(D_DECL(i,j,k),u_imaterial)
           enddo ! u_imaterial=1,nmat

           vofcomp_raw=(im_dest-1)*ngeom_raw+1
           vofcomp_recon=(im_source-1)*ngeom_recon+1

           do dir=1,SDIM
            if ((newvfrac(im_dest).gt.zero).and. &
                (newvfrac(im_dest).le.one)) then
             new_centroid(im_dest,dir)= &
               unsplit_snew(vofcomp_raw+dir)+cengrid(dir)

             ! all the source material can be converted into
             ! destination material.
            else if ((newvfrac(im_dest).eq.zero).and. &
                     (oldvfrac(im_source).gt.zero)) then
             new_centroid(im_dest,dir)= &
               recon(D_DECL(i,j,k),vofcomp_recon+dir)+cengrid(dir)
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

           if (away_from_interface.eq.1) then
            newvfrac(im_source)=oldvfrac(im_source)
            newvfrac(im_dest)=oldvfrac(im_dest)
           else if (away_from_interface.eq.0) then
            LSnew(D_DECL(i,j,k),im_source)=unsplit_lsnew(im_source)
            LSnew(D_DECL(i,j,k),im_dest)=unsplit_lsnew(im_dest)
           else
            print *,"away_from_interface invalid"
            stop
           endif

            ! old mass and volume fraction of vapor to gas
           Yfrac_vapor_to_gas=one
           vfrac_vapor_to_gas=one

           mtype=fort_material_type(im_dest)
           if (mtype.eq.0) then
            density_dest=fort_denconst(im_dest)
           else if ((mtype.ge.1).and.(mtype.le.fort_max_num_eos)) then
            if (newvfrac(im_dest).gt.EBVOFTOL) then
             density_dest=unsplit_density(im_dest) 
            else if (newvfrac(im_dest).le.EBVOFTOL) then
             density_dest=fort_denconst(im_dest)
            else
             print *,"newvfrac(im_dest) invalid"
             stop
            endif
           else
            print *,"mtype invalid"
            stop
           endif

           mtype=fort_material_type(im_source)
           if (mtype.eq.0) then
            density_source=fort_denconst(im_source)
           else if ((mtype.ge.1).and.(mtype.le.fort_max_num_eos)) then
            if (newvfrac(im_source).gt.EBVOFTOL) then
             density_source=unsplit_density(im_source) 
            else if (newvfrac(im_source).le.EBVOFTOL) then
             if (oldvfrac(im_source).gt.EBVOFTOL) then
              dencomp=(im_source-1)*num_state_material+1
              density_source=EOS(D_DECL(i,j,k),dencomp) 
             else if (oldvfrac(im_source).le.EBVOFTOL) then
              density_source=fort_denconst(im_source)
             else
              print *,"oldvfrac(im_source) invalid"
              stop
             endif
            else
             print *,"newvfrac(im_source) invalid"
             stop
            endif
           else
            print *,"mtype invalid"
            stop
           endif
            
           dFdst=(newvfrac(im_dest)-oldvfrac(im_dest))
           dFsrc=(oldvfrac(im_source)-newvfrac(im_source))

           if (local_freezing_model.eq.0) then ! standard Stefan model
            ! do nothing
           else if ((local_freezing_model.eq.4).or. & ! Tannasawa or Schrage
                    (local_freezing_model.eq.5).or. & ! Stefan evap/cond model
                    (local_freezing_model.eq.6).or. & ! Palmore/Desjardins
                    (local_freezing_model.eq.7)) then ! Cavitation

            if ((mass_frac_id.ge.1).and. &
                (mass_frac_id.le.num_species_var)) then

             evap_den=species_evaporation_density(mass_frac_id)

             if (evap_den.gt.zero) then

              if (LL.gt.zero) then ! evaporation

               if (oldvfrac(im_dest).gt.EBVOFTOL) then ! valid Y,F vapor
                speccompdst=(im_dest-1)*num_state_material+ &
                 num_state_base+mass_frac_id
                Yfrac_old=EOS(D_DECL(i,j,k),speccompdst)
                if ((Yfrac_old.ge.zero).and.(Yfrac_old.le.one)) then
                 density_old=density_dest !denconst if material_id==0
                 call make_mixture_density(Yfrac_old, &
                   density_old,evap_den)
                 Yfrac_vapor_to_gas=Yfrac_old 
                 vfrac_vapor_to_gas=density_old*Yfrac_old/evap_den
                else
                 print *,"Yfrac_old invalid"
                 stop
                endif
               else if (oldvfrac(im_dest).le.EBVOFTOL) then ! default Y,F vapor
                Yfrac_vapor_to_gas=zero
                vfrac_vapor_to_gas=zero
               else
                print *,"oldvfrac(im_dest) invalid"
                stop
               endif

              else if (LL.lt.zero) then ! condensation

                  ! mass fraction equation:
                  ! in a given cell with m species.
                  ! mass= sum_i=1^m  Y_i overall_mass = 
                  !     = sum_i=1^m  density_i F_i V_cell
                  ! F_i = volume fraction of material i.
                  ! (rho Y_i)_t + div (rho u Y_i) = div rho D_i  grad Y_i
                  ! (Y_i)_t + div (u Y_i) = div rho D_i  grad Y_i/rho
               if (oldvfrac(im_source).gt.EBVOFTOL) then
                speccompsrc=(im_source-1)*num_state_material+ &
                 num_state_base+mass_frac_id
                Yfrac_old=EOS(D_DECL(i,j,k),speccompsrc)
                if ((Yfrac_old.ge.zero).and.(Yfrac_old.le.one)) then
                 density_old=density_source
                 call make_mixture_density(Yfrac_old, &
                   density_old,evap_den)
                 Yfrac_vapor_to_gas=Yfrac_old
                 vfrac_vapor_to_gas=density_old*Yfrac_old/evap_den
                else
                 print *,"Yfrac_old invalid"
                 stop
                endif
               else if (oldvfrac(im_source).le.EBVOFTOL) then
                Yfrac_vapor_to_gas=zero
                vfrac_vapor_to_gas=zero
               else
                print *,"oldvfrac(im_source) invalid"
                stop
               endif

              else
               print *,"LL invalid"
               stop
              endif

             else
              print *,"evap_den invalid"
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
             
           if (1.eq.0) then
            print *,"i,j,k,dFdst,dFsrc,olddst,oldsrc ",i,j,k,dFdst,dFsrc, &
             oldvfrac(im_dest),oldvfrac(im_source)
           endif

           dF=max(dFdst,dFsrc)
           if (LL.gt.zero) then ! evaporation, boiling, melting
            dF=min(dF,oldvfrac(im_source))
           else if (LL.lt.zero) then ! freezing, condensation
            dF=min(dF,vfrac_vapor_to_gas*oldvfrac(im_source))
           else
            print *,"LL invalid"
            stop
           endif

           newvfrac(im_dest)=oldvfrac(im_dest)
           newvfrac(im_source)=oldvfrac(im_source)

           if (dF.lt.zero) then
            dF=zero
           endif
           newvfrac(im_dest)=oldvfrac(im_dest)+dF
           if (newvfrac(im_dest).gt.one) then
            dF=one-oldvfrac(im_dest)
            newvfrac(im_dest)=one
           endif
           newvfrac(im_source)=oldvfrac(im_source)-dF
           if (newvfrac(im_source).lt.zero) then
            dF=oldvfrac(im_source)
            newvfrac(im_source)=zero
            newvfrac(im_dest)=oldvfrac(im_dest)+dF
            if (newvfrac(im_dest).gt.one) then
             newvfrac(im_dest)=one
            endif
           endif
           if (dF.le.zero) then
            dF=zero
           endif
           if (dF.ge.one) then
            dF=one
           endif

           dF=newvfrac(im_dest)-oldvfrac(im_dest)
           DVOF(im_dest)=DVOF(im_dest)+dF
           deltaVOF(D_DECL(i,j,k),im_dest)=dF

           if (abs(dF).gt.EBVOFTOL) then

            dencomp=num_materials_vel*(SDIM+1)+ &
                    (im_dest-1)*num_state_material+1
            tcomp=dencomp+1

            if (newvfrac(im_dest).gt.EBVOFTOL) then

             snew(D_DECL(i,j,k),tcomp)=unsplit_temperature(im_dest) 

              ! evaporation: im_source=liquid im_dest=surrounding gas
              ! condensation: im_source=surrounding gas im_dest=liquid
             if ((local_freezing_model.eq.4).or. & !Tanasawa or Schrage
                 (local_freezing_model.eq.5).or. & !Stefan model evap/cond.
                 (local_freezing_model.eq.6).or. & !Palmore/Desjardins
                 (local_freezing_model.eq.7)) then !Cavitation
              speccompdst=num_materials_vel*(SDIM+1)+ &
                  (im_dest-1)*num_state_material+2+mass_frac_id
            
               ! evaporation 
               ! assume all liquid becomes pure vapor: dF (expansion comes 
               ! later)
               ! total new gases: oldvfrac(im_dest)+dF
               ! total new vapor: oldvfrac(im_dest)*vfrac_vapor_to_gas+dF
               ! new vapor/(total new gases) = new_vfrac_vapor
              if (LL.gt.zero) then
               new_vfrac_vapor=(oldvfrac(im_dest)*vfrac_vapor_to_gas+dF)/ &
                (oldvfrac(im_dest)+dF)
               if ((new_vfrac_vapor.ge.zero).and. &
                   (new_vfrac_vapor.le.one)) then
                density_new=new_vfrac_vapor*evap_den+ &
                 (one-new_vfrac_vapor)*density_dest
                snew(D_DECL(i,j,k),speccompdst)= &
                 new_vfrac_vapor*evap_den/density_new
               else
                print *,"new_vfrac_vapor invalid"
                stop
               endif

                ! condensation: vapor fraction=1 in the liquid.
              else if (LL.lt.zero) then 
               snew(D_DECL(i,j,k),speccompdst)=one
              else
               print *,"LL invalid"
               stop
              endif
             endif

             mtype=fort_material_type(im_dest)
             if (mtype.eq.0) then
              ! density modified in FORT_DENCOR
             else if ((mtype.ge.1).and.(mtype.le.fort_max_num_eos)) then
              snew(D_DECL(i,j,k),dencomp)=density_dest 
             else
              print *,"mtype invalid"
              stop
             endif

             if ((unsplit_temperature(im_dest).lt.TEMPERATURE_FLOOR).or. &
                 (unsplit_density(im_dest).le.zero)) then
              print *,"unsplit_temperature(im_dest) invalid   or"
              print *,"unsplit_density(im_dest) invalid"
              print *,"oldvfrac(im_dest) ",oldvfrac(im_dest)
              print *,"newvfrac(im_dest) ",newvfrac(im_dest)
              print *,"unsplit_density(im_dest)=", &
               unsplit_density(im_dest)
              print *,"unsplit_temperature(im_dest)=", &
               unsplit_temperature(im_dest)
              print *,"oldvfrac(im_source) ",oldvfrac(im_source)
              print *,"newvfrac(im_source) ",newvfrac(im_source)
              print *,"unsplit_density(im_source)=", &
               unsplit_density(im_source)
              print *,"unsplit_temperature(im_source)=", &
               unsplit_temperature(im_source)
              stop
             endif

            else if (newvfrac(im_dest).le.EBVOFTOL) then
             snew(D_DECL(i,j,k),tcomp)=Tsat_default
             if ((local_freezing_model.eq.4).or. & ! Tanasawa or Schrage
                 (local_freezing_model.eq.5).or. & ! Stefan Evap/Cond.
                 (local_freezing_model.eq.6).or. & ! Palmore/Desjardins
                 (local_freezing_model.eq.7)) then ! Cavitation
              speccompdst=num_materials_vel*(SDIM+1)+ &
                   (im_dest-1)*num_state_material+2+mass_frac_id

              if (LL.gt.zero) then ! evaporation
               snew(D_DECL(i,j,k),speccompdst)=zero
              else if (LL.lt.zero) then ! condensation
               snew(D_DECL(i,j,k),speccompdst)=one
              else
               print *,"LL invalid"
               stop
              endif

             endif
            else
             print *,"newvfrac(im_dest) invalid"
             stop
            endif

            dencomp=num_materials_vel*(SDIM+1)+ &
                    (im_source-1)*num_state_material+1
            tcomp=dencomp+1

            if (newvfrac(im_source).gt.EBVOFTOL) then

             mtype=fort_material_type(im_source)
             if (mtype.eq.0) then
              ! source density modified in FORT_DENCOR
             else if ((mtype.ge.1).and.(mtype.le.fort_max_num_eos)) then
              snew(D_DECL(i,j,k),dencomp)=density_source
             else
              print *,"mtype invalid"
              stop
             endif

             snew(D_DECL(i,j,k),tcomp)=unsplit_temperature(im_source) 

             if ((local_freezing_model.eq.4).or. & ! Tanasawa or Schrage
                 (local_freezing_model.eq.5).or. & ! Stefan evap/cond.
                 (local_freezing_model.eq.6).or. & ! Palmore/Desjardins
                 (local_freezing_model.eq.7)) then ! Cavitation
              speccompsrc=num_materials_vel*(SDIM+1)+ &
                   (im_source-1)*num_state_material+2+mass_frac_id
            
               ! condensation
              if (LL.lt.zero) then
               new_vfrac_vapor=(oldvfrac(im_source)*vfrac_vapor_to_gas-dF)/ &
                (oldvfrac(im_source)-dF)
               if ((new_vfrac_vapor.ge.zero).and. &
                   (new_vfrac_vapor.le.one)) then
                density_new=new_vfrac_vapor*evap_den+ &
                 (one-new_vfrac_vapor)*density_source
                snew(D_DECL(i,j,k),speccompsrc)= &
                 new_vfrac_vapor*evap_den/density_new
               else
                print *,"new_vfrac_vapor invalid"
                stop
               endif
              else if (LL.gt.zero) then ! evaporation
               snew(D_DECL(i,j,k),speccompsrc)=one
              else
               print *,"LL invalid"
               stop
              endif
             endif

             if ((unsplit_temperature(im_source).lt.TEMPERATURE_FLOOR).or. &
                 (unsplit_density(im_source).le.zero)) then
              print *,"unsplit_temperature(im_source) invalid    or "
              print *,"unsplit_density(im_source) invalid "
              stop
             endif
            else if (newvfrac(im_source).le.EBVOFTOL) then
             snew(D_DECL(i,j,k),tcomp)=Tsat_default
             if ((local_freezing_model.eq.4).or. & ! Tanasawa or Schrage
                 (local_freezing_model.eq.5).or. & ! Stefan evap/cond.
                 (local_freezing_model.eq.6).or. & ! Palmore/Desjardins
                 (local_freezing_model.eq.7)) then ! Cavitation
              speccompsrc=num_materials_vel*(SDIM+1)+ &
                       (im_source-1)*num_state_material+2+mass_frac_id
              if (LL.gt.zero) then ! evaporation
               snew(D_DECL(i,j,k),speccompsrc)=one
              else if (LL.lt.zero) then ! condensation
               snew(D_DECL(i,j,k),speccompsrc)=zero
              else
               print *,"LL invalid"
               stop
              endif

             endif
            else
             print *,"newvfrac(im_source) invalid"
             stop
            endif

            do dir=1,SDIM
             snew(D_DECL(i,j,k),vcompdst_snew+dir)= &
              new_centroid(im_dest,dir)-cengrid(dir)
             snew(D_DECL(i,j,k),vcompsrc_snew+dir)= &
              new_centroid(im_source,dir)-cengrid(dir)
            enddo ! dir
            if (ngeom_raw.eq.SDIM+1) then
             ! do nothing
            else
             print *,"ngeom_raw invalid in convert material"
             print *,"ngeom_raw= ",ngeom_raw
             stop
            endif

           endif ! |dF|>EBVOFTOL

          else if (isweep.eq.1) then

           iten=iten_crit
           ireverse=ireverse_crit
           im_dest=im_dest_crit
           im_source=im_source_crit 

           if (interface_near(iten+ireverse*nten).ne.1) then
            print *,"interface_near invalid"
            stop
           endif

           Tsat_default=saturation_temp(iten+ireverse*nten)
           Tsat_flag=NINT(TSATFAB(D_DECL(i,j,k),iten))
           if (ireverse.eq.0) then
             ! do nothing
           else if (ireverse.eq.1) then
             Tsat_flag=-Tsat_flag
           else
             print *,"ireverse invalid"
             stop
           endif
            
           if ((Tsat_flag.eq.1).or.(Tsat_flag.eq.2)) then
             Tsat_default=TSATFAB(D_DECL(i,j,k), &
              nten+(iten-1)*ncomp_per_tsat+1)
           else if ((Tsat_flag.eq.-1).or. &
                    (Tsat_flag.eq.-2)) then
             ! do nothing
           else if (Tsat_flag.eq.0) then
             ! do nothing
           else
             print *,"Tsat_flag invalid"
             stop
           endif

           LL=latent_heat(iten+ireverse*nten)
           local_freezing_model=freezing_model(iten+ireverse*nten)
           distribute_from_targ=distribute_from_target(iten+ireverse*nten)
           mass_frac_id=mass_fraction_id(iten+ireverse*nten)

           if ((local_freezing_model.lt.0).or.(local_freezing_model.gt.7)) then
            print *,"local_freezing_model invalid 2"
            stop
           endif
           if ((distribute_from_targ.lt.0).or.(distribute_from_targ.gt.1)) then
            print *,"distribute_from_targ invalid"
            stop
           endif

           vcompsrc_snew=num_materials_vel*(SDIM+1)+ &
            nmat*num_state_material+(im_source-1)*ngeom_raw+1
           vcompdst_snew=num_materials_vel*(SDIM+1)+ &
            nmat*num_state_material+(im_dest-1)*ngeom_raw+1

           dcompsrc=(im_source-1)*num_state_material+1
           dcompdst=(im_dest-1)*num_state_material+1

           dcompdst_snew=num_materials_vel*(SDIM+1)+dcompdst
           tcompdst_snew=dcompdst_snew+1
    
           do u_imaterial=1,nmat
            oldLS_point(u_imaterial)=LSold(D_DECL(i,j,k),u_imaterial)
            vofcomp_recon=(u_imaterial-1)*ngeom_recon+1
            oldvfrac(u_imaterial)=recon(D_DECL(i,j,k),vofcomp_recon)
           enddo

            ! deltaVOF init when isweep==0
           dF=deltaVOF(D_DECL(i,j,k),im_dest)
           dF=DVOF_FACT(im_dest)*dF
           if ((DVOF_FACT(im_dest).lt.zero).or. &
               (DVOF_FACT(im_dest).gt.one)) then
            print *,"DVOF_FACT invalid"
            stop
           endif
           newvfrac(im_dest)=oldvfrac(im_dest)+dF
           newvfrac(im_source)=oldvfrac(im_source)-dF

           if (newvfrac(im_dest).gt.one+VOFTOL) then
            print *,"newvfrac(im_dest) overflow"
            stop
           else if (newvfrac(im_dest).gt.one) then
            newvfrac(im_dest)=one
           endif
           if (newvfrac(im_source).lt.-VOFTOL) then
            print *,"newvfrac(im_source) underflow"
            stop
           else if (newvfrac(im_source).lt.zero) then
            newvfrac(im_source)=zero
           endif

#if (STANDALONE==0)
           ksource=get_user_heatviscconst(im_source)
           kdest=get_user_heatviscconst(im_dest)
#elif (STANDALONE==1)
           ksource=fort_heatviscconst(im_source)
           kdest=fort_heatviscconst(im_dest)
#else
           print *,"bust compiling convertmaterial"
           stop
#endif

           densrc=EOS(D_DECL(i,j,k),dcompsrc)
           dendst=EOS(D_DECL(i,j,k),dcompdst)

           densrc_restrict=densrc
           dendst_restrict=dendst

           if ((densrc.le.zero).or.(dendst.le.zero)) then
            print *,"densrc and dendst must be positive"
            stop
           endif

           mtype=fort_material_type(im_source)
           if (mtype.eq.0) then
            ! do nothing
           else if ((mtype.ge.1).and.(mtype.le.fort_max_num_eos)) then
            if (densrc_restrict.lt.density_floor_expansion(im_source)) then
             densrc_restrict=density_floor_expansion(im_source)
            endif
            if (densrc_restrict.gt.density_ceiling_expansion(im_source)) then
             densrc_restrict=density_ceiling_expansion(im_source)
            endif
           else
            print *,"mtype invalid"
            stop
           endif

           mtype=fort_material_type(im_dest)
           if (mtype.eq.0) then
            ! do nothing
           else if ((mtype.ge.1).and.(mtype.le.fort_max_num_eos)) then
            if (dendst_restrict.lt.density_floor_expansion(im_dest)) then
             dendst_restrict=density_floor_expansion(im_dest)
            endif
            if (dendst_restrict.gt.density_ceiling_expansion(im_dest)) then
             dendst_restrict=density_ceiling_expansion(im_dest)
            endif
           else
            print *,"mtype invalid"
            stop
           endif

           if ((local_freezing_model.eq.4).or. & ! Tanasawa or Schrage
               (local_freezing_model.eq.5).or. & ! Stefan evap/cond.
               (local_freezing_model.eq.6).or. & ! Palmore/Desjardins
               (local_freezing_model.eq.7)) then ! cavitation
            if ((mass_frac_id.ge.1).and. &
                (mass_frac_id.le.num_species_var)) then
             evap_den=species_evaporation_density(mass_frac_id)
             if (evap_den.gt.zero) then
              if (LL.gt.zero) then ! evaporation
               dendst_restrict=evap_den
              else if (LL.lt.zero) then ! condensation
               densrc_restrict=evap_den
              else
               print *,"LL invalid"
               stop
              endif
             else
              print *,"evap_den invalid"
              stop
             endif
            else
             print *,"mass_frac_id invalid"
             stop
            endif
           endif ! local_freezing_model==4,5,6,7

           if ((densrc_restrict.le.zero).or.(dendst_restrict.le.zero)) then
            print *,"densrc_restrict and dendst_restrict must be positive"
            stop
           endif

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
            denratio_factor=densrc_restrict/dendst_restrict-one

           !dF=mdot/rho_dest
           !dF_expand_source=(1-rho_dest/rho_source)*dF
           !dM=rho_dest * dF_dest + rho_source * dF_source=
           !rho_dest(mdot/rho_dest)+ 
           !rho_source (-mdot/rho_dest+(1-rho_dest/rho_source)(dF))=
           !mdot-(rho_source/rho_dest)mdot+(rho_source-rho_dest)dF=
           !mdot-(rho_source/rho_dest)mdot+(rho_source/rho_dest)mdot-mdot=0
           else if (distribute_from_targ.eq.1) then
            ! distribute div u source to the cells in which F_dest<1/2  

            denratio_factor=one-dendst_restrict/densrc_restrict

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
            ! do nothing
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

             if ((oldvfrac(im_dest).lt.half).and. &
                 (newvfrac(im_dest).gt.half)) then
              if (dF.le.zero) then
               print *,"dF invalid"
               stop
              endif
              SWEPTFACTOR=(newvfrac(im_dest)-half)/dF
              if (SWEPTFACTOR.le.LSTOL) then
               SWEPTFACTOR=LSTOL
              endif
              if ((SWEPTFACTOR.ge.LSTOL).and.(SWEPTFACTOR.le.one)) then
               ! do nothing
              else
               print *,"SWEPTFACTOR invalid"
               stop
              endif
              swept(D_DECL(i,j,k))=SWEPTFACTOR
             else if ((oldvfrac(im_dest).ge.half).or. &
                      (newvfrac(im_dest).le.half)) then
              ! do nothing
             else
              print *,"oldvfrac(im_dest) or newvfrac(im_dest) invalid"
              stop
             endif

! single (continuum method) temperature equation for both phases.
! source term at the interface.
! latent_heat<0 condensation or solidification
! latent_heat>0 boiling or melting
! units of specific heat: J/(kg K)
! units of latent heat: J/kg
            else if (local_freezing_model.eq.1) then

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
            else if (local_freezing_model.eq.2) then

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
   
              call Hydrate_energy_source_term(dF,dt,ksource, &
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

            snew(D_DECL(i,j,k),vcompsrc_snew)=newvfrac(im_source)
            snew(D_DECL(i,j,k),vcompdst_snew)=newvfrac(im_dest)

            delta_mass(im_source)=delta_mass(im_source)+ &
             volgrid*(newvfrac(im_source)-oldvfrac(im_source))
            delta_mass(im_dest+nmat)=delta_mass(im_dest+nmat)+ &
             volgrid*(newvfrac(im_dest)-oldvfrac(im_dest))

           else
            print *,"dF bust"
            stop
           endif
       
          else
           print *,"isweep invalid"
           stop
          endif

         else if (do_unsplit_advection.eq.0) then
          ! do nothing
         else
          print *,"do_unsplit_advection invalid"
          stop
         endif

        else if (is_rigid(nmat,im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(nmat,im) invalid"
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
       density_floor_expansion, &
       density_ceiling_expansion, &
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
       latent_heat, &
       use_exact_temperature, &
       reaction_rate, &
       saturation_temp, &
       saturation_temp_curv, &
       saturation_temp_vel, &
       freezing_model, &
       Tanasawa_or_Schrage, &
       distribute_from_target, &
       mass_fraction_id, &
       species_evaporation_density, &
       molar_mass, &
       species_molar_mass, &
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


      INTEGER_T, intent(in) :: stefan_flag
      INTEGER_T, target, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: normal_probe_size
      REAL_T :: microscale_probe_size
      INTEGER_T, intent(in) :: ngrow_distance
      INTEGER_T, intent(in) :: nstate
      INTEGER_T, target, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nburning
      INTEGER_T, intent(in) :: ntsat
      INTEGER_T, intent(in) :: nden
      REAL_T, target, intent(in) :: density_floor_expansion(nmat)
      REAL_T, target, intent(in) :: density_ceiling_expansion(nmat)
      INTEGER_T, intent(in) :: custom_nucleation_model
      INTEGER_T, intent(in) :: do_the_nucleate
      INTEGER_T, intent(in) :: nucleate_pos_size
      REAL_T, intent(in) :: nucleate_pos(nucleate_pos_size)
      REAL_T, intent(in) :: nucleation_temp(2*nten)
      REAL_T, intent(in) :: nucleation_pressure(2*nten)
      REAL_T, intent(in) :: nucleation_pmg(2*nten)
      REAL_T, intent(in) :: nucleation_mach(2*nten)
      REAL_T, intent(in) :: cavitation_pressure(nmat)
      REAL_T, intent(in) :: cavitation_vapor_density(nmat)
      REAL_T, intent(in) :: cavitation_tension(nmat)
      INTEGER_T, intent(in) ::  microlayer_substrate(nmat)
      REAL_T, intent(in) :: microlayer_angle(nmat)
      REAL_T, intent(in) :: microlayer_size(nmat)
      REAL_T, intent(in) :: macrolayer_size(nmat)
      REAL_T, intent(in) :: max_contact_line_size(nmat)
      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: use_exact_temperature(2*nten)
      REAL_T, intent(in) :: reaction_rate(2*nten)
      REAL_T :: K_f(0:1)
      REAL_T, intent(in) :: saturation_temp(2*nten)
      REAL_T, intent(in) :: saturation_temp_curv(2*nten)
      REAL_T, intent(in) :: saturation_temp_vel(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: Tanasawa_or_Schrage(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      INTEGER_T, intent(in) :: mass_fraction_id(2*nten)
      REAL_T, intent(in) :: molar_mass(nmat)
      REAL_T, intent(in) :: species_molar_mass(num_species_var+1)
      REAL_T, intent(in) :: species_evaporation_density(num_species_var+1)
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
      REAL_T, intent(inout) :: LSnew(DIMV(LSnew),nmat*(SDIM+1))
      REAL_T, intent(inout) :: Snew(DIMV(Snew),nstate)
      REAL_T, intent(in) :: LS_slopes_FD(DIMV(LS_slopes_FD),nmat*SDIM)
      REAL_T, target, intent(in) :: EOS(DIMV(EOS),nden)
       ! F,X,order,SL,I x nmat
      REAL_T, target, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon) 
      REAL_T, target, intent(in) :: pres(DIMV(pres)) 
      REAL_T, intent(in) :: pres_eos(DIMV(pres_eos)) 
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
      REAL_T T_I_interp(2)
      REAL_T Y_I_interp(2)
      REAL_T pres_I_interp(2)
      REAL_T local_Tsat(0:1)
      REAL_T local_Tsat_base(0:1)
      REAL_T vel_phasechange(0:1)
      REAL_T, target :: LL(0:1)
      INTEGER_T valid_phase_change(0:1)
      REAL_T, target :: dxprobe_source
      REAL_T, target :: dxprobe_dest
      REAL_T dxprobe_target(2)
      REAL_T ksource,kdest
      REAL_T LS_pos
      REAL_T C_w0
      INTEGER_T, target :: local_freezing_model
      INTEGER_T local_Tanasawa_or_Schrage
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
      REAL_T evap_den
      REAL_T source_perim_factor,dest_perim_factor
      REAL_T contact_line_perim
      INTEGER_T icolor,base_type,ic,im1,im2
      REAL_T normal_probe_factor
      INTEGER_T iprobe
      REAL_T Y_TOLERANCE
       ! iten=1..nten  ireverse=0..1
      REAL_T temp_target_probe_history(2*nten,2)
      REAL_T dxprobe_target_history(2*nten,2)

      INTEGER_T use_tsatfab
      INTEGER_T ncomp_per_burning
      INTEGER_T ncomp_per_tsat

      REAL_T CURV_OUT_I

      REAL_T VEL_predict,VEL_correct
      REAL_T Y_predict
      REAL_T Y_correct
      REAL_T X_predict
      REAL_T X_correct
      REAL_T TSAT_predict,TSAT_correct
      REAL_T TSAT_ERR,TSAT_INIT_ERR
      INTEGER_T TSAT_iter,TSAT_converge,TSAT_iter_max
      INTEGER_T YMIN_converge
      REAL_T Y_interface_min
      REAL_T X_interface_min
      INTEGER_T YMIN_iter
      INTEGER_T YMIN_iter_max
      REAL_T denom
      REAL_T FicksLawD(2)  ! iprobe=1 source iprobe=2 dest 
      REAL_T Tprobe_avg 
      REAL_T molar_mass_ambient
      REAL_T molar_mass_vapor
      REAL_T T_interface_min
      REAL_T YMIN_ERR
      REAL_T YMIN_INIT_ERR
      REAL_T GRAD_Y_dot_n
      INTEGER_T interp_valid_flag(2) ! iprobe=1 source iprobe=2 dest
      type(probe_parm_type) :: PROBE_PARMS

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

      Y_TOLERANCE=0.01

      if ((stefan_flag.eq.0).or.(stefan_flag.eq.1)) then
       ! do nothing
      else
       print *,"stefan_flag invalid"
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

      do im=1,nmat
       if ((density_floor_expansion(im).le.zero).or. &
           (density_floor_expansion(im).gt.fort_denconst(im))) then
        print *,"density_floor_expansion invalid"
        stop
       endif
       if ((density_ceiling_expansion(im).le.zero).or. &
           (density_ceiling_expansion(im).lt.fort_denconst(im))) then
        print *,"density_ceiling_expansion invalid"
        stop
       endif
      enddo ! im=1..nmat

      if (num_materials_scalar_solve.eq.1) then ! GFM
       normal_probe_factor=half
      else if (num_materials_scalar_solve.eq.nmat) then ! FVM multimat
       normal_probe_factor=half
      else
       print *,"num_materials_scalar_solve invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(typefab),1,-1,6625)
      call checkbound(fablo,fabhi,DIMS(colorfab),1,-1,6626)

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,1251)

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
      call checkbound(fablo,fabhi,DIMS(LSnew),1,-1,1253)
      call checkbound(fablo,fabhi,DIMS(Snew),1,-1,1253)
      call checkbound(fablo,fabhi,DIMS(LS_slopes_FD),1,-1,1253)
      call checkbound(fablo,fabhi,DIMS(EOS),ngrow,-1,1254)
      call checkbound(fablo,fabhi,DIMS(pres),ngrow,-1,1255)
      call checkbound(fablo,fabhi,DIMS(pres_eos),1,-1,1255)

      !blob_matrix,blob_RHS,blob_velocity,
      !blob_integral_momentum,blob_energy,
      !blob_mass_for_velocity (3 comp)
      !blob_volume, 
      !blob_center_integral,blob_center_actual
      !blob_perim, blob_perim_mat, blob_triple_perim, 
      if (num_elements_blobclass.ne. &
          3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
          2*(2*SDIM)+1+ &
          3+1+2*SDIM+1+nmat+nmat*nmat) then
       print *,"num_elements_blobclass invalid rate mass change:", &
         num_elements_blobclass
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

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

         ! LEVELSET FUNCTION AT CELL CENTERS YANG.
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
            local_Tanasawa_or_Schrage=Tanasawa_or_Schrage(iten+ireverse*nten)

            ispec=mass_fraction_id(iten+ireverse*nten)

            evap_den=one

            if ((ispec.ge.0).and.(ispec.le.num_species_var)) then
             ! do nothing
            else
             print *,"ispec invalid"
             stop
            endif

            if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
             evap_den=species_evaporation_density(ispec)
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

             if ((is_rigid(nmat,im).eq.1).or. &
                 (is_rigid(nmat,im_opp).eq.1)) then

              ! do nothing

             else if (LL(ireverse).ne.zero) then

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

                dencomp_source=(im_source-1)*num_state_material+1
                dencomp_dest=(im_dest-1)*num_state_material+1

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
                  Tsatfab,DIMS(Tsatfab), &
                  local_Tsat(ireverse), &
                  iten+ireverse*nten, &
                  saturation_temp, &
                  use_exact_temperature, &
                  xI,cur_time,nmat,nten,7)

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

                local_Tsat(ireverse)=local_Tsat(ireverse)- &
                  saturation_temp_curv(iten+ireverse*nten)* &
                  CURV_OUT_I

                Y_predict=one
                Y_correct=Y_predict
               
                X_predict=one
                X_correct=X_predict

                TSAT_predict=local_Tsat(ireverse)
                TSAT_correct=TSAT_predict

                VEL_predict=zero
                VEL_correct=zero

                TSAT_iter=0
                TSAT_iter_max=5
                YMIN_iter_max=5
                TSAT_converge=0

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
                PROBE_PARMS%dxprobe_source=>dxprobe_source
                PROBE_PARMS%dxprobe_dest=>dxprobe_dest
                PROBE_PARMS%local_freezing_model=>local_freezing_model
                PROBE_PARMS%LL=>LL(ireverse)
                PROBE_PARMS%debugrate=>debugrate
                PROBE_PARMS%i=>i
                PROBE_PARMS%j=>j
                PROBE_PARMS%k=>k
                PROBE_PARMS%EOS=>EOS 
                PROBE_PARMS%LS=>LS  ! PROBE_PARMS%LS is pointer, LS is target
                PROBE_PARMS%recon=>recon
                PROBE_PARMS%pres=>pres
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
                PROBE_PARMS%dxmaxLS=>dxmaxLS
                PROBE_PARMS%bfact=>bfact
                PROBE_PARMS%level=>level
                PROBE_PARMS%finest_level=>finest_level
                PROBE_PARMS%dx=>dx
                PROBE_PARMS%xlo=>xlo
                PROBE_PARMS%xI=>xI
                PROBE_PARMS%nmat=>nmat
                PROBE_PARMS%ngrow=>ngrow
                PROBE_PARMS%fablo=>fablo
                PROBE_PARMS%fabhi=>fabhi
                PROBE_PARMS%density_floor_expansion=>density_floor_expansion
                PROBE_PARMS%density_ceiling_expansion=>density_ceiling_expansion

!FIX ME
! 1. Y BC in diffusion solver
! 2. div ( rho D grad Y )/rho
!    local_freezing_model=4  Tanasawa or Schrage
!    local_freezing_model=5  fully saturated evaporation?
!    local_freezing_model=6  partially saturated evaporation?
!    local_freezing_model=7  Cavitation (a seed must exist)

                do while (TSAT_converge.eq.0) 

                 !iprobe=1 source
                 !iprobe=2 dest
                 call probe_interpolation( &
                  PROBE_PARMS, &
                  TSAT_predict,Y_predict, &
                  T_probe,Y_probe, &
                  den_I_interp,T_I_interp,Y_I_interp, &
                  pres_I_interp, &
                  T_probe_raw, &
                  dxprobe_target, &
                  interp_valid_flag, &
                  at_interface)

                 !iprobe=1 source
                 !iprobe=2 dest

                 if (at_interface.eq.1) then
                     
#if (STANDALONE==0)
                  ksource=get_user_heatviscconst(im_source)
                  kdest=get_user_heatviscconst(im_dest)
#elif (STANDALONE==1)
                  ksource=fort_heatviscconst(im_source)
                  kdest=fort_heatviscconst(im_dest)
#else
                  print *,"bust compiling ratemasschange"
                  stop
#endif

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
                    if (evap_den.gt.zero) then
                     if (LL(ireverse).gt.zero) then ! evaporation
                      den_I_interp(2)=evap_den ! dest
                     else if (LL(ireverse).lt.zero) then ! condensation
                      den_I_interp(1)=evap_den ! source
                     else
                      print *,"LL invalid"
                      stop
                     endif
                    else
                     print *,"evap_den invalid"
                     stop
                    endif  
                   else
                    print *,"ispec invalid"
                    stop
                   endif
                  endif ! local_freezing_model==5,6

                  source_perim_factor=one
                  dest_perim_factor=one

                  if ((max_contact_line_size(im_source).gt.zero).and. &
                      (microlayer_size(im_source).gt.zero).and. &
                      (macrolayer_size(im_source).gt.zero).and. &
                      (microlayer_substrate_source.ge.1).and. &
                      (microlayer_substrate_source.le.nmat)) then
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
                   ! perim, perim_mat
                   ic=(icolor-1)*num_elements_blobclass+ &
                    3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
                    2*(2*SDIM)+1+ &
                    3+1+2*SDIM+1+nmat

                   im2=microlayer_substrate_source
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
                        max_contact_line_size(im_source))) then
                    ! do nothing
                   else if (contact_line_perim.gt. &
                            max_contact_line_size(im_source)) then
                    source_perim_factor= &
                      max_contact_line_size(im_source)/contact_line_perim
                   else
                    print *,"contact_line_perim invalid"
                    stop
                   endif
                  else if ((max_contact_line_size(im_source).eq.zero).or. &
                           (microlayer_size(im_source).eq.zero).or. &
                           (macrolayer_size(im_source).eq.zero).or. &
                           (microlayer_substrate_source.eq.0)) then
                   ! do nothing
                  else
                   print *,"microlayer parameters invalid"
                   stop
                  endif

                  if ((max_contact_line_size(im_dest).gt.zero).and. &
                      (microlayer_size(im_dest).gt.zero).and. &
                      (macrolayer_size(im_dest).gt.zero).and. &
                      (microlayer_substrate_dest.ge.1).and. &
                      (microlayer_substrate_dest.le.nmat)) then
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
                   ! perim, perim_mat
                   ic=(icolor-1)*num_elements_blobclass+ &
                    3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
                    2*(2*SDIM)+1+ &
                    3+1+2*SDIM+1+nmat

                   im2=microlayer_substrate_dest
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
                        max_contact_line_size(im_dest))) then
                     ! do nothing
                   else if (contact_line_perim.gt. &
                            max_contact_line_size(im_dest)) then
                    dest_perim_factor= &
                      max_contact_line_size(im_dest)/contact_line_perim
                   else
                    print *,"contact_line_perim invalid"
                    stop
                   endif
                  else if ((max_contact_line_size(im_dest).eq.zero).or. &
                           (microlayer_size(im_dest).eq.zero).or. &
                           (macrolayer_size(im_dest).eq.zero).or. &
                           (microlayer_substrate_dest.eq.0)) then
                   ! do nothing
                  else
                   print *,"microlayer parameters invalid"
                   stop
                  endif

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
                    local_freezing_model, &
                    local_Tanasawa_or_Schrage, &
                    evap_den, &
                    distribute_from_targ, &
                    VEL_correct, & ! vel
                    den_I_interp(1), & ! source 
                    den_I_interp(2), & ! dest
                    ksource,kdest, & ! ksrc,kdst
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
                    velsrc=ksource*DTsrc/(LL(ireverse)*dxprobe_target(1))
                    veldst=kdest*DTdst/(LL(ireverse)*dxprobe_target(2))
                  
                    velsum=velsrc+veldst
                    if (velsum.gt.zero) then
                     ! do nothing
                    else if (velsum.le.zero) then
                     velsum=zero
                    else
                     print *,"velsum invalid"
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
                   molar_mass_vapor=species_molar_mass(ispec)
                   if (LL(ireverse).gt.zero) then ! evaporation
                    iprobe=2  ! destination
                    molar_mass_ambient=molar_mass(im_dest)
                   else if (LL(ireverse).lt.zero) then ! condensation
                    iprobe=1  ! source
                    molar_mass_ambient=molar_mass(im_source)
                   else
                    print *,"LL invalid"
                    stop
                   endif

                   if ((molar_mass_ambient.gt.zero).and. &
                       (molar_mass_vapor.gt.zero).and. &
                       (R_Palmore_Desjardins.gt.zero)) then

                    if ((Y_probe(iprobe).ge.one-Y_TOLERANCE).and. &
                        (Y_probe(iprobe).le.one)) then
                     ! do nothing
                    else if ((Y_probe(iprobe).le.one-Y_TOLERANCE).and. &
                             (Y_probe(iprobe).ge.zero)) then
                      !Y_probe<=Y_interface<=1
                     if (TSAT_iter.eq.0) then
                      YMIN_converge=0
                      Y_interface_min=zero
                      YMIN_iter=0
                      do while (YMIN_converge.eq.0)

                       call probe_interpolation( &
                        PROBE_PARMS, &
                        TSAT_predict,Y_interface_min, &
                        T_probe,Y_probe, &
                        den_I_interp,T_I_interp,Y_I_interp, &
                        pres_I_interp, &
                        T_probe_raw, &
                        dxprobe_target, &
                        interp_valid_flag, &
                        at_interface)

                       YMIN_ERR=abs(Y_probe(iprobe)-Y_interface_min)
                       Y_interface_min=Y_probe(iprobe)
                       if (YMIN_iter.eq.0) then
                        YMIN_INIT_ERR=YMIN_ERR
                       endif
                       YMIN_converge=0
                       if (YMIN_ERR.eq.zero) then
                        YMIN_converge=1
                       endif
                       YMIN_iter=YMIN_iter+1
                       if (YMIN_iter.gt.YMIN_iter_max) then
                        YMIN_converge=1
                       endif
                       if ((Y_probe(iprobe).ge.one-Y_TOLERANCE).and. &
                           (Y_probe(iprobe).le.one)) then
                        YMIN_converge=1
                       else if ((Y_probe(iprobe).ge.zero).and. &
                                (Y_probe(iprobe).le.one)) then
                        ! do nothing
                       else
                        print *,"Y_probe invalid"
                        stop
                       endif
                       if (YMIN_iter.gt.1) then
                        if (YMIN_err.lt.(0.001d0)*YMIN_INIT_ERR) then
                         YMIN_converge=1
                        endif
                       endif
                      enddo ! do while (YMIN_converge.eq.0)

                      if ((Y_probe(iprobe).ge.one-Y_TOLERANCE).and. &
                          (Y_probe(iprobe).le.one)) then
                       ! do nothing
                      else if ((Y_probe(iprobe).le. &
                                one-Y_TOLERANCE).and. &
                               (Y_probe(iprobe).ge.zero)) then
                       X_interface_min=molar_mass_ambient*Y_interface_min/ &
                        ((one-Y_interface_min)*molar_mass_vapor+ &
                         Y_interface_min*molar_mass_ambient)
                       if ((X_interface_min.ge.one-Y_TOLERANCE).and. &
                           (X_interface_min.le.one)) then
                        ! do nothing
                       else if ((X_interface_min.le. &
                                 one-Y_TOLERANCE).and. &
                                (X_interface_min.gt.zero)) then
                        TSAT_correct=one/ &
                         (one/local_Tsat_base(ireverse)- &
                         R_Palmore_Desjardins*log(X_interface_min)/ &
                         (abs(LL(ireverse))*molar_mass_vapor)) 
                        T_interface_min=TSAT_correct
                        X_correct=X_interface_min
                        Y_correct=Y_interface_min
                       else if (X_interface_min.eq.zero) then
                        TSAT_correct=zero
                        T_interface_min=zero
                        X_correct=zero
                        Y_correct=zero
                       else
                        print *,"X_interface_min invalid"
                        stop
                       endif
                      else
                       print *,"Y_probe invalid"
                       stop
                      endif
                     else if (TSAT_iter.ge.1) then
                      denom=one/dxprobe_target(1)+one/dxprobe_target(2)
                      if (denom.gt.zero) then
                       Tprobe_avg=T_probe(1)/dxprobe_target(1)+ &
                           T_probe(2)/dxprobe_target(2)
                       if (Tprobe_avg.ge.zero) then
                        Tprobe_avg=Tprobe_avg/denom
                        if (den_I_interp(iprobe).gt.zero) then
                          ! LL>0 melting   LL<0 freezing
                          ! do nothing
                        else
                         print *,"den_I_interp(iprobe) invalid"
                         stop
                        endif
                        if (TSAT_correct.gt.zero) then
                         X_correct=exp(-abs(LL(ireverse))*molar_mass_vapor/ &
                          R_Palmore_Desjardins)*(one/TSAT_correct- &
                           one/local_Tsat_base(ireverse))
                        else if (TSAT_correct.eq.zero) then
                         X_correct=zero
                         Y_correct=zero
                        else
                         print *,"TSAT_correct invalid"
                         stop
                        endif
                        if ((X_correct.ge.zero).and.(X_correct.le.one)) then
                         Y_correct=X_correct*molar_mass_vapor/ &
                          ((one-X_correct)*molar_mass_ambient+ &
                           X_correct*molar_mass_vapor)
                         if (Y_correct.le.Y_interface_min) then
                          Y_correct=Y_interface_min
                          X_correct=X_interface_min
                          TSAT_correct=T_interface_min
                         else if (Y_correct.ge.one-Y_TOLERANCE) then
                          Y_correct=one
                          X_correct=one
                          TSAT_correct=local_Tsat(ireverse)
                         else if ((Y_correct.gt.Y_interface_min).and. &
                                  (Y_correct.lt.one-Y_TOLERANCE)) then
                          if (dxprobe_target(iprobe).gt.zero) then
                           GRAD_Y_dot_n= &
                             (Y_correct-Y_probe(iprobe))/ &
                             dxprobe_target(iprobe) 
                           if (GRAD_Y_dot_n.ge.zero) then
                            TSAT_correct=Tprobe_avg- &
                              (one/denom)* &
                              FicksLawD(iprobe)* &
                              den_I_interp(iprobe)* &
                              LL(ireverse)* &
                              GRAD_Y_dot_n/(one-Y_correct)
                           else
                            print *,"GRAD_Y_dot_n invalid"
                            stop
                           endif 
                          else
                           print *,"dxprobe_target(iprobe) invalid"
                           stop
                          endif
                         else 
                          print *,"Y_correct invalid"
                          stop
                         endif
                        else
                         print *,"X_correct invalid"
                         stop
                        endif 
                       else
                        print *,"Tprobe_avg invalid"
                        stop
                       endif
                      else
                       print *,"denom invalid"
                       stop
                      endif
                     else
                      print *,"TSAT_iter invalid"
                      stop
                     endif
                    else
                     print *,"Y_probe invalid"
                     stop
                    endif

                   else
                    print *,"molar masses, or R invalid"
                    stop
                   endif
 
                  else if (local_freezing_model.ge.0) then
                   ! do nothing
                  else
                   print *,"local_freezing_model invalid 7"
                   stop
                  endif

                 else if (at_interface.eq.0) then

                 else
                  print *,"at_interface invalid"
                  stop
                 endif

                 TSAT_correct=TSAT_correct- &
                    saturation_temp_vel(iten+ireverse*nten)* &
                    (VEL_correct-VEL_predict)
            
                 VEL_predict=VEL_correct
 
                 TSAT_ERR=abs(TSAT_correct-TSAT_predict)

                 TSAT_predict=TSAT_correct

                 if (TSAT_iter.eq.0) then
                   TSAT_INIT_ERR=TSAT_ERR
                 endif

                 TSAT_converge=0
                 if (TSAT_ERR.eq.zero) then
                   TSAT_converge=1
                 endif
                  ! TSAT_iter starts at 0
                 TSAT_iter=TSAT_iter+1
                 if (TSAT_iter.gt.TSAT_iter_max) then
                   TSAT_converge=1
                 endif
                 if (TSAT_iter.gt.1) then
                   if (TSAT_err.lt.(0.001d0)*TSAT_INIT_ERR) then
                    TSAT_converge=1
                   endif
                 endif
                 if (at_interface.eq.0) then
                   TSAT_converge=1
                 endif

                 if (at_interface.eq.1) then
                  if (1.eq.0) then
                   print *,"i,j,k,TSAT_iter,TSAT_iter_max,TSAT_ERR ", &
                    i,j,k,TSAT_iter,TSAT_iter_max,TSAT_ERR
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
                  print *,"ksource,kdest,local_Tsat(ireverse) ", &
                         ksource,kdest,local_Tsat(ireverse)
                  print *,"T_Probe(1),T_probe(2) ",T_Probe(1),T_probe(2)
                  print *,"den_I_interp(1) ",den_I_interp(1)
                  print *,"LSINTsrc,LSINTdst ",LSINT(im_source),LSINT(im_dest)
                  print *,"nrmCP ",nrmCP(1),nrmCP(2),nrmCP(SDIM)
                  print *,"nrmFD ",nrmFD(1),nrmFD(2),nrmFD(SDIM)
                  print *,"nrmPROBE ",nrmPROBE(1),nrmPROBE(2),nrmPROBE(SDIM)
                  print *,"im_dest= ",im_dest
                 endif
   
                else if (at_interface.eq.0) then
                 ! do nothing
                else
                 print *,"at_interface invalid"
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
              Tsatfab(D_DECL(i,j,k),nten+ncomp_per_tsat*(iten-1)+2)= &
               one  ! default mass fraction=1 (saturated)
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

