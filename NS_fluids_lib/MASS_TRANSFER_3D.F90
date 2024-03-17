#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

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

#include "EXTRAP_COMP.H"
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
      use amrex_fort_module, only : amrex_real
      use probcommon_module

      ! if t1 is a c++ parameter and p1 is a "type" component:
      ! real(amrex_real), pointer :: p1(DIMV(t1),ncomp)
      ! real(amrex_real), target :: t1(DIMV(t1),ncomp)
      ! p1=>t1
      type probe_parm_type
       integer :: tid
       integer :: probe_constrain
       real(amrex_real), pointer :: Y_TOLERANCE
       integer, pointer :: local_freezing_model
       real(amrex_real), pointer :: LL
       integer, pointer :: debugrate
       integer, pointer :: i,j,k
       real(amrex_real), pointer :: xsrc(:)
       real(amrex_real), pointer :: xdst(:)
       real(amrex_real), pointer :: xsrc_micro(:)
       real(amrex_real), pointer :: xdst_micro(:)
       integer, pointer :: im_source
       integer, pointer :: im_dest
       integer, pointer :: tcomp_source
       integer, pointer :: Ycomp_source
       integer, pointer :: dencomp_source
       real(amrex_real), pointer :: dxprobe_source
       integer, pointer :: tcomp_dest
       integer, pointer :: Ycomp_dest
       integer, pointer :: dencomp_dest
       real(amrex_real), pointer :: dxprobe_dest
       real(amrex_real), pointer, dimension(:) :: LSINT
       real(amrex_real), pointer :: dxmaxLS
       integer, pointer :: bfact
       integer, pointer :: level
       integer, pointer :: finest_level
       real(amrex_real), pointer :: dx(:)
       real(amrex_real), pointer :: xlo(:)
       real(amrex_real), pointer :: xI(:)
       integer, pointer :: fablo(:)
       integer, pointer :: fabhi(:)
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: EOS
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: recon
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: LS
       real(amrex_real), pointer, dimension(D_DECL(:,:,:)) :: pres
      end type probe_parm_type

      type probe_out_type
      real(amrex_real) :: T_probe(2)
      real(amrex_real) :: T_probe_raw(2)
      real(amrex_real) :: Y_probe(2)
      real(amrex_real) :: den_I_interp(2)
      real(amrex_real) :: den_probe(2)
      real(amrex_real) :: T_I_interp(2)
      real(amrex_real) :: Y_I_interp(2)
      real(amrex_real) :: pres_I_interp(2)
      real(amrex_real) :: vfrac_I(2) !solids and fluids tessellate
      real(amrex_real) :: dxprobe_target(2)
      integer :: interp_valid_flag(2)
      end type probe_out_type

      type TSAT_MASS_FRAC_parm_type
       type(probe_parm_type), pointer :: PROBE_PARMS
       integer :: Tanasawa_or_Schrage_or_Kassemi
       real(amrex_real) :: accommodation_coefficient
       real(amrex_real) :: reference_pressure
       real(amrex_real) :: universal_gas_constant_R
       real(amrex_real) :: molar_mass_ambient   
       real(amrex_real) :: molar_mass_vapor
       real(amrex_real) :: Clausius_Clapyron_Tsat
       real(amrex_real) :: TSAT_base
       real(amrex_real) :: YI_min
       real(amrex_real) :: TI_min
       real(amrex_real) :: TI_max
       real(amrex_real) :: D_MASS
       real(amrex_real) :: den_G
       integer :: iprobe_vapor
       integer, pointer :: material_type_evap(:)
       real(amrex_real), pointer :: thermal_k(:)
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
        fablo,fabhi, &
        TSATFAB, &
        T_I, & !intent(out)
        latent_comp, &
        TSAT_array, & !intent(in)
        TSAT_flag_array, &  ! =0 if derived, =1 if T_interface(x,y,z,t) given.
        x_I,time)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: use_tsatfab
      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: ireverse
      integer, INTENT(in) :: iten
      integer, INTENT(in) :: ntsat
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      real(amrex_real), pointer, INTENT(in) :: TSATFAB(D_DECL(:,:,:),:)
      integer, INTENT(in) :: latent_comp
      real(amrex_real), INTENT(out) :: T_I
      real(amrex_real), INTENT(in)  :: TSAT_array(2*num_interfaces)
      integer, INTENT(in)  :: TSAT_flag_array(2*num_interfaces)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: x_I(SDIM)
      real(amrex_real), INTENT(in) :: time
      integer :: local_flag
      integer :: dir
      integer :: ncomp_per_tsat
      integer :: tsat_comp

      ncomp_per_tsat=EXTRAP_PER_TSAT
      tsat_comp=num_interfaces+(iten-1)*ncomp_per_tsat+1

      if (ntsat.eq.EXTRAP_NCOMP_TSAT) then
       ! do nothing
      else
       print *,"ntsat invalid"
       stop
      endif
      if (bfact.lt.1) then 
       print *,"bfact invalid114"
       stop
      endif
      if ((iten.ge.1).and.(iten.le.num_interfaces)) then
       ! do nothing
      else
       print *,"iten invalid"
       stop
      endif

      if ((latent_comp.ge.1).and.(latent_comp.le.2*num_interfaces)) then
       if (time.ge.zero) then
        if (use_tsatfab.eq.1) then
         call interpfab_tsat( &
           i,j,k, &
           ireverse, &
           iten, &
           ntsat, &
           bfact, &
           level, &
           finest_level, &
           dx,xlo, &
           x_I, &
           tsat_comp, &
           fablo,fabhi, &
           TSATFAB, &
           T_I)
        else if (use_tsatfab.eq.0) then
         local_flag=TSAT_flag_array(latent_comp)
         if (local_flag.eq.0) then
          T_I=TSAT_array(latent_comp)
         else if ((local_flag.ge.1).and.(local_flag.le.num_materials)) then
          print *,"local_flag <> 0 not supported yet"
          stop
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
         print *,"put breakpoint here to see the caller"
         print *,"T_I out of range T_I=",T_I
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

      return
      end subroutine get_interface_temperature

      subroutine adjust_du(du,normdir,rval,map_forward)
      IMPLICIT NONE

      real(amrex_real) du,rval,disc
      integer normdir,map_forward


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
      if ((levelrz.eq.COORDSYS_RZ).or. &
          (levelrz.eq.COORDSYS_CYLINDRICAL)) then
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
      else if (levelrz.eq.COORDSYS_CARTESIAN) then
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
       im_critical, &
       im_primary_sten, &
       VF_sten, &
       LS_sten, &
       TSAT, & !unused if tsat_flag==-1, 0 (but set to 293.0 for sanity check)
       T_out)
      use global_utility_module
      implicit none

      real(amrex_real), INTENT(in) :: DATA_FLOOR
      integer, INTENT(in) :: nsolve
      integer :: nc
      integer, INTENT(in) :: cc_flag
      integer, INTENT(in) :: tsat_flag
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: xI(SDIM)
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: TSAT
      real(amrex_real), INTENT(inout) :: T_sten(D_DECL(-1:1,-1:1,-1:1),nsolve)
      real(amrex_real), INTENT(in) :: XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      integer, INTENT(in) :: im_critical
      integer :: im_primary
      integer, INTENT(in) :: im_primary_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real), INTENT(in) :: VF_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real), INTENT(in) :: LS_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) :: wt_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real), INTENT(out) :: T_out(nsolve)

      integer i,j,k
      integer klosten,khisten
      integer i1,j1,k1
      integer dir
      real(amrex_real) dist
      real(amrex_real), dimension(:,:), allocatable :: AA
      real(amrex_real), dimension(:,:), allocatable :: AAcopy
      real(amrex_real) BB(SDIM+1)
      real(amrex_real) BBcopy(SDIM+1)
      real(amrex_real) xtemp(SDIM)
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      real(amrex_real) xbase(SDIM)
      real(amrex_real) xlive(SDIM)
      real(amrex_real) delx(SDIM+1)
      real(amrex_real) GRADTEMP(SDIM+1)

      real(amrex_real) TMIN,TMAX
      integer TBOUNDS_INIT

      real(amrex_real) T_test,T_avg,wt_sum,VF,LS
      integer mat_ncomp,matstatus
      real(amrex_real) T_hold
      integer own_flag
      real(amrex_real) wt_local

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
      if (TSAT.ge.DATA_FLOOR) then
       ! do nothing
      else
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
       if (dx(dir).gt.zero) then 
        ! do nothing
       else
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

       TMAX=-1.0D+20
       TMIN=1.0D+20
       TBOUNDS_INIT=0

       if (tsat_flag.eq.0) then
        ! do nothing
       else if (tsat_flag.eq.1) then
        TMAX=TSAT
        TMIN=TSAT
        TBOUNDS_INIT=1
       else if (tsat_flag.eq.-1) then
        ! do nothing
       else
        print *,"tsat_flag invalid"
        stop
       endif

       T_avg=zero
       wt_sum=zero
       
       do k=klosten,khisten
       do j=-1,1
       do i=-1,1

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
         im_primary=im_primary_sten(D_DECL(i,j,k))

         if ((im_primary.ge.1).and.(im_primary.le.num_materials)) then
          ! do nothing
         else if (im_primary.eq.-1) then
          ! do nothing
         else
          print *,"im_primary invalid"
          stop
         endif

         VF=VF_sten(D_DECL(i,j,k))
         LS=LS_sten(D_DECL(i,j,k))
        else if (tsat_flag.eq.-1) then ! use all cells in the stencil
         im_primary=-1
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
         if ((VF.ge.-EPS1).and. &
             (VF.le.VOFTOL)) then
          own_flag=0
         else if ((VF.ge.VOFTOL).and.(VF.le.one+EPS1)) then
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
          print *,"tsat_flag invalid:",tsat_flag
          stop
         endif
         wt_local=EPS_8_4
        else if (own_flag.eq.1) then
          ! if normal probe is target: xtarget = xprobe
          ! if new supermesh centroid is target: xtarget = x_new_centroid
          ! if centroid to center: xtarget = xcenter
          ! if center to centroid: xtarget = xcentroid
         dist=EPS_8_4
         do dir=1,SDIM
          dist=dist+(xlive(dir)-xtarget(dir))**2/(dx(1)**2)
         enddo ! dir
         if (dist.gt.zero) then
          ! do nothing
         else
          print *,"dist invalid"
          stop
         endif
         wt_local=one/dist
         if (tsat_flag.eq.-1) then ! use all cells in the stencil
          ! do nothing
         else if ((tsat_flag.eq.1).or. &
                  (tsat_flag.eq.0)) then
          if (VF.gt.zero) then
           ! do nothing
          else
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

        if (wt_local.gt.zero) then
         ! do nothing
        else
         print *,"wt_local invalid"
         stop
        endif

        wt_sten(D_DECL(i,j,k))=wt_local

        T_test=T_sten(D_DECL(i,j,k),nc)

        if ((im_primary.eq.im_critical).or. &
            (im_primary.eq.-1)) then

         if (TBOUNDS_INIT.eq.0) then
          TBOUNDS_INIT=1
          TMAX=T_test
          TMIN=T_test
         else if (TBOUNDS_INIT.eq.1) then
          if (T_test.gt.TMAX) then
           TMAX=T_test
          endif
          if (T_test.lt.TMIN) then
           TMIN=T_Test
          endif
         else 
          print *,"TBOUNDS_INIT invalid"
          stop
         endif

        else if ((im_primary.ne.im_critical).and. &
                 (im_primary.ge.1).and. &
                 (im_primary.le.num_materials)) then

         ! do nothing

        else
         print *,"im_primary invalid"
         stop
        endif

        if (abs(T_test).lt.1.0D+30) then
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

       if (wt_sum.gt.zero) then
        ! do nothing
       else
        print *,"wt_sum invalid"
        stop
       endif

       do k=klosten,khisten
       do j=-1,1
       do i=-1,1
        wt_sten(D_DECL(i,j,k))=wt_sten(D_DECL(i,j,k))/wt_sum
       enddo
       enddo
       enddo

       if (tsat_flag.eq.1) then
        T_avg=T_avg/wt_sum
        if (T_avg.ge.DATA_FLOOR) then
         ! do nothing
        else
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

       if (TBOUNDS_INIT.eq.0) then
        TMIN=T_avg
        TMAX=T_avg
       else if (TBOUNDS_INIT.eq.1) then
        ! do nothing
       else
        print *,"TBOUNDS_INIT invalid"
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
         if (T_hold.ge.zero) then
          ! do nothing
         else
          print *,"temperature underflow in center_centroid_interchange"
          print *,"xI ",xI(1),xI(2),xI(SDIM)
          print *,"xtarget ",xtarget(1),xtarget(2),xtarget(SDIM)
          print *,"nsolve=",nsolve
          print *,"cc_flag=",cc_flag
          print *,"tsat_flag=",tsat_flag
          print *,"bfact=",bfact
          print *,"level=",level
          print *,"finest_level=",finest_level
          print *,"TSAT = ",TSAT
          print *,"TMIN = ",TMIN
          print *,"TMAX = ",TMAX
          print *,"T_avg = ",T_avg
          print *,"DATA_FLOOR= ",DATA_FLOOR
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

      enddo ! do nc=1,nsolve

      return
      end subroutine center_centroid_interchange


      subroutine interpfabFWEIGHT( &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo,x, &
       im, &
       comp, &
       lo,hi, &
       data, &
       recon, &
       dest)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: x(SDIM)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
      integer, INTENT(in) :: comp
      integer, INTENT(in) :: im
      integer :: im_local
       ! datalox:datahix,dataloy:datahiy,dataloz:datahiz
      real(amrex_real), pointer, INTENT(in) :: data(D_DECL(:,:,:),:)
      real(amrex_real), pointer, INTENT(in) :: recon(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out) :: dest

      real(amrex_real) :: DATA_FLOOR

      real(amrex_real) :: T_out(1)

      integer dir
      integer ic,jc,kc
      integer i1,j1,k1
      integer isten,jsten,ksten
      integer k1lo,k1hi
      integer, parameter :: nhalf=3
      integer vofcomp
      integer cell_index(SDIM)

      real(amrex_real) :: local_VOF(num_materials)

      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_stencil(-nhalf:nhalf,SDIM)
      real(amrex_real) volcell
      real(amrex_real) cencell(SDIM)
      integer im_primary_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) T_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) VF_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      integer cc_flag
      integer tsat_flag
      integer nsolve
      real(amrex_real) Tsat
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:),:)
      real(amrex_real) local_data_out

      DATA_FLOOR=zero

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

      if ((im.lt.1).or.(im.gt.num_materials)) then
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

       ! declared in GLOBALUTIL.F90
      call containing_cell(bfact,dx,xlo,lo,x,cell_index)

      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

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

       local_data_fab=>data
       call safe_data(isten,jsten,ksten,comp,local_data_fab,local_data_out)
       T_sten(D_DECL(i1,j1,k1))=local_data_out

       do im_local=1,num_materials
        vofcomp=(im_local-1)*ngeom_recon+1
        local_data_fab=>recon
        call safe_data(isten,jsten,ksten,vofcomp,local_data_fab,local_data_out)
        local_VOF(im_local)=local_data_out
       enddo !im_local=1..num_materials

       call get_primary_material_VFRAC(local_VOF, &
         im_primary_sten(D_DECL(i1,j1,k1)))

       VF_sten(D_DECL(i1,j1,k1))=local_VOF(im)

       vofcomp=(im-1)*ngeom_recon+1
       local_data_fab=>recon
       do dir=1,SDIM
        call safe_data(isten,jsten,ksten,vofcomp+dir, &
          local_data_fab,local_data_out)
        XC_sten(D_DECL(i1,j1,k1),dir)=local_data_out+cencell(dir)
       enddo

      enddo
      enddo
      enddo ! i1,j1,k1

       ! in: interpfabFWEIGHT
      cc_flag=0  ! centroid -> target
      tsat_flag=0 ! do not use TSAT
      nsolve=1
      Tsat=room_temperature ! 293.0d0 if double prec.
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
       im, &
       im_primary_sten, & 
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
       lo,hi, &
       recon, & ! fluids tess, solids overlay
       dest)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: x(SDIM)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
      real(amrex_real), pointer, INTENT(in) :: recon(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out) :: dest(num_materials)

      integer im
      integer ic,jc,kc
      integer, parameter :: nhalf=3
      integer vofcomp
      integer cell_index(SDIM)

      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) volcell
      real(amrex_real) cencell(SDIM)
      integer nmax
      integer local_tessellate
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:),:)
      real(amrex_real) local_data_out

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

      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      call gridsten_level(xsten,ic,jc,kc,level,nhalf)
      call Box_volumeFAST(bfact,dx,xsten,nhalf, &
        volcell,cencell,SDIM)

      local_data_fab=>recon
      do im=1,num_materials*ngeom_recon
       call safe_data(ic,jc,kc,im,local_data_fab,local_data_out)
       mofdata(im)=local_data_out
      enddo

      local_tessellate=3
       !EPS2
      call multi_get_volume_tessellate( &
        local_tessellate, & ! =3
        bfact, &
        dx, &
        xsten,nhalf, &
        mofdata, &
        geom_xtetlist(1,1,1,tid+1), &
        nmax, &
        nmax, &
        SDIM)

      do im=1,num_materials
       vofcomp=(im-1)*ngeom_recon+1
       dest(im)=mofdata(vofcomp)
      enddo

      return 
      end subroutine interpfabVFRAC_tess


      subroutine grad_probe_sanity(xI,xprobe,temp_probe,Tsat,LL)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: LL
      real(amrex_real), INTENT(in) :: xI(SDIM)
      real(amrex_real), INTENT(in) :: xprobe(SDIM)
      real(amrex_real), INTENT(in) :: temp_probe
      real(amrex_real), INTENT(in) :: Tsat

      real(amrex_real) :: grad_probe_local
      real(amrex_real) :: mag
      real(amrex_real) :: nrm(SDIM)

      integer dir

       ! nrm points from xI to xprobe

      if (LL.ne.zero) then

       mag=zero
       do dir=1,SDIM
        nrm(dir)=xprobe(dir)-xI(dir)
        mag=mag+nrm(dir)**2
       enddo
       mag=sqrt(mag)
       if (mag.gt.zero) then
        ! do nothing
       else
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
       print *,"LL invalid: ",LL
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
       lo,hi, &
       data, &
       dest)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      integer, INTENT(in) :: lo(SDIM)
      integer, INTENT(in) :: hi(SDIM)
      integer, INTENT(in) :: comp
      real(amrex_real), pointer, INTENT(in) :: data(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out) :: dest
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:),:)
      real(amrex_real) local_data_out

      real(amrex_real) :: DATA_FLOOR

      real(amrex_real) :: T_out(1)

      integer ic,jc,kc
      integer i1,j1,k1
      integer isten,jsten,ksten
      integer k1lo,k1hi
      integer, parameter :: nhalf=3
      integer cell_index(SDIM)

      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer, PARAMETER :: im_critical=-1
      integer im_primary_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) T_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      integer cc_flag
      integer tsat_flag
      integer nsolve
      real(amrex_real) Tsat

      DATA_FLOOR=zero

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
      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      call gridsten_level(xsten,ic,jc,kc,level,nhalf)

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+ic
       jsten=j1+jc
       ksten=k1+kc
       local_data_fab=>data
       call safe_data(isten,jsten,ksten,comp,local_data_fab,local_data_out)
       T_sten(D_DECL(i1,j1,k1))=local_data_out

      enddo
      enddo
      enddo ! i1,j1,k1

       ! in: interpfab
      cc_flag=1  ! center -> target
      tsat_flag=-1  ! use all cells in the stencil
      nsolve=1
      !unused if tsat_flag==-1 
      !(but set to room_temperature=293.0d0 for sanity check)
      Tsat=room_temperature 
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
       im_critical, &
       im_primary_sten, &  ! (not used)
       T_sten, &  ! VF_sten (not used)
       T_sten, &  ! LS_Sten (not used)
       Tsat, & ! unused if tsat_flag==-1
       T_out)

      dest=T_out(1)

      return 
      end subroutine interpfab

      subroutine single_interpfab( &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       lo,hi, &
       data, &
       dest)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      integer, INTENT(in) :: lo(SDIM)
      integer, INTENT(in) :: hi(SDIM)
      real(amrex_real), pointer, INTENT(in) :: data(D_DECL(:,:,:))
      real(amrex_real), INTENT(out) :: dest
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:))
      real(amrex_real) local_data_out

      real(amrex_real) :: DATA_FLOOR

      real(amrex_real) :: T_out(1)

      integer ic,jc,kc
      integer i1,j1,k1
      integer isten,jsten,ksten
      integer k1lo,k1hi
      integer, parameter :: nhalf=3
      integer cell_index(SDIM)

      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer, PARAMETER :: im_critical=-1
      integer im_primary_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) T_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      integer cc_flag
      integer tsat_flag
      integer nsolve
      real(amrex_real) Tsat

      DATA_FLOOR=zero

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

      call containing_cell(bfact,dx,xlo,lo,xtarget,cell_index)
      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      call gridsten_level(xsten,ic,jc,kc,level,nhalf)

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+ic
       jsten=j1+jc
       ksten=k1+kc

       local_data_fab=>data
       call safe_data_single(isten,jsten,ksten,local_data_fab,local_data_out)
       T_sten(D_DECL(i1,j1,k1))=local_data_out

      enddo
      enddo
      enddo ! i1,j1,k1

       ! in: interpfab
      cc_flag=1  ! center -> target
      tsat_flag=-1  ! use all cells in the stencil
      nsolve=1
      !unused if tsat_flag==-1 
      !(but set to room_temperature=293.0d0 for sanity check)
      Tsat=room_temperature 
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
       im_critical, &  ! (not used)
       im_primary_sten, &  ! (not used)
       T_sten, &  ! VF_sten (not used)
       T_sten, &  ! LS_Sten (not used)
       Tsat, & ! unused if tsat_flag==-1
       T_out)

      dest=T_out(1)

      return 
      end subroutine single_interpfab



      subroutine interpfab_tsat( &
       i,j,k, &
       ireverse, &
       iten, &
       ntsat, &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       comp, &
       lo,hi, &
       data, &
       TSAT)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: ireverse
      integer, INTENT(in) :: iten
      integer, INTENT(in) :: ntsat
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
      integer, INTENT(in) :: comp
      real(amrex_real), pointer, INTENT(in) :: data(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out) :: TSAT

      integer ncomp_per_tsat
      integer k1lo,k1hi
      integer cell_index(3)
      integer dir
      integer, parameter :: nhalf=3
      real(amrex_real) TSAT_times_weight
      real(amrex_real) TSAT_weight
      integer i1,j1,k1
      integer isten,jsten,ksten
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer TSAT_FLAG
      real(amrex_real) local_TSAT
      real(amrex_real) local_weight
      real(amrex_real) eps
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:),:)
      real(amrex_real) local_data_out

      ncomp_per_tsat=EXTRAP_PER_TSAT
      if (ntsat.eq.EXTRAP_NCOMP_TSAT) then
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
      if ((comp.gt.num_interfaces).and.(comp.le.ntsat)) then
       ! do nothing
      else
       print *,"comp out of range"
       stop
      endif
      if ((iten.ge.1).and.(iten.le.num_interfaces)) then
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

      TSAT_times_weight=zero
      TSAT_weight=zero

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+i
       jsten=j1+j
       ksten=k1+k
       call gridsten_level(xsten,isten,jsten,ksten,level,nhalf)

       local_data_fab=>data
       call safe_data(isten,jsten,ksten,iten,local_data_fab,local_data_out)
       TSAT_FLAG=NINT(local_data_out)
       if (ireverse.eq.0) then
        ! do nothing
       else if (ireverse.eq.1) then
        TSAT_FLAG=-TSAT_FLAG
       else
        print *,"ireverse invalid"
        stop
       endif
       if ((TSAT_FLAG.eq.1).or.(TSAT_FLAG.eq.2)) then
        call safe_data(isten,jsten,ksten,comp,local_data_fab,local_data_out)
        local_TSAT=local_data_out
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
       print *,"i,j,k,ireverse,iten,num_interfaces,ntsat,bfact ", &
               i,j,k,ireverse,iten,num_interfaces,ntsat,bfact
       print *,"level,finest_level,dx ",level,finest_level, &
               dx(1),dx(2),dx(SDIM)
       print *,"xlo,xtarget,comp ",xlo(1),xlo(2),xlo(SDIM), &
               xtarget(1),xtarget(2),xtarget(SDIM)
       print *,"TSAT_weight ",TSAT_weight
       stop
      endif

      return 
      end subroutine interpfab_tsat

       ! i,j,k is the cell containing xtarget
      subroutine interpfab_curv( &
       interp_status, &  !interp_status=1 good, =0 not valid.
       curv_comp, &
       bfact, &
       level, &
       finest_level, &
       dx, &
       xlo, &
       xtarget, &
       lo,hi, &
       data_in, &
       CURV_OUT)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(out) :: interp_status
      integer, INTENT(in) :: curv_comp
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
       ! first num_materials+num_interfaces components are curvatures
       ! second num_materials+num_interfaces components are status (0=bad 1=good)
      real(amrex_real), pointer, INTENT(in) :: data_in(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out) :: CURV_OUT

      integer k1lo,k1hi
      integer cell_index(3)
      integer dir
      integer, parameter :: nhalf=3
      integer CURV_FLAG
      real(amrex_real) CURV_times_weight
      real(amrex_real) CURV_weight
      integer i1,j1,k1
      integer isten,jsten,ksten
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) local_CURV
      real(amrex_real) local_weight
      real(amrex_real) eps
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:),:)
      real(amrex_real) local_data_out

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
      if ((curv_comp.ge.1).and.(curv_comp.le.num_materials+num_interfaces)) then
       ! do nothing
      else
       print *,"curv_comp out of range"
       stop
      endif

      call containing_cell(bfact,dx,xlo,lo,xtarget,cell_index)

      CURV_times_weight=zero
      CURV_weight=zero

      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi

       isten=i1+cell_index(1)
       jsten=j1+cell_index(2)
       ksten=k1+cell_index(3)

       call gridsten_level(xsten,isten,jsten,ksten,level,nhalf)

       local_data_fab=>data_in
       call safe_data(isten,jsten,ksten,num_materials+num_interfaces+curv_comp, &
           local_data_fab,local_data_out)
       CURV_FLAG=NINT(local_data_out)

       if (CURV_FLAG.eq.1) then
         call safe_data(isten,jsten,ksten,curv_comp, &
           local_data_fab,local_data_out)
         local_CURV=local_data_out
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

      enddo
      enddo
      enddo ! i1,j1,k1

      if (CURV_weight.gt.zero) then
       CURV_OUT=CURV_times_weight/CURV_weight
       interp_status=1
      else if (CURV_weight.eq.zero) then
       CURV_OUT=zero
       interp_status=0
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
       lo,hi, &
       data, &
       dest)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
      integer, INTENT(in) :: comp
      real(amrex_real), pointer, INTENT(in) :: data(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out) :: dest

      integer ic,jc,kc
      integer cell_index(SDIM)
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:),:)
      real(amrex_real) local_data_out

      if (bfact.lt.1) then 
       print *,"bfact invalid115"
       stop
      endif
      if ((comp.lt.1).or.(comp.gt.1000)) then
       print *,"comp out of range"
       stop
      endif

      call containing_cell(bfact,dx,xlo,lo,xtarget,cell_index)
      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

      local_data_fab=>data
      call safe_data(ic,jc,kc,comp,local_data_fab,local_data_out)
      dest=local_data_out

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
       im, &
       comp, &
       lo,hi, &
       tempfab, &
       LS, &
       recon, &
       dest, &
       debugrate)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      integer, INTENT(in) :: debugrate
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: xI(SDIM)
      real(amrex_real), INTENT(in) :: Tsat
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
      integer, INTENT(in) :: comp
      integer, INTENT(in) :: im
      integer :: im_local
      real(amrex_real), pointer, INTENT(in) :: tempfab(D_DECL(:,:,:),:)
      real(amrex_real), pointer, INTENT(in) :: LS(D_DECL(:,:,:),:)
      real(amrex_real), pointer, INTENT(in) :: recon(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(out) :: dest
      real(amrex_real), pointer :: local_data_fab(D_DECL(:,:,:),:)
      real(amrex_real) local_data_out

      real(amrex_real) :: DATA_FLOOR

      real(amrex_real) :: T_out(1)

      integer dir
      integer ic,jc,kc
      integer i1,j1,k1
      integer isten,jsten,ksten
      integer k1lo,k1hi
      integer, parameter :: nhalf=3
      integer vofcomp
      integer cell_index(SDIM)

      real(amrex_real) :: local_VOF(num_materials)

      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_stencil(-nhalf:nhalf,SDIM)
      real(amrex_real) volcell
      real(amrex_real) cencell(SDIM)
      integer im_primary_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) T_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) VF_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) LS_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      integer cc_flag
      integer tsat_flag
      integer nsolve

      DATA_FLOOR=zero

      if (bfact.lt.1) then 
       print *,"bfact invalid116"
       stop
      endif
      if ((comp.lt.1).or.(comp.gt.1000)) then
       print *,"comp out of range"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid21"
       stop
      endif
      if ((Tsat.ge.zero).and.(Tsat.le.1.0D+30)) then
       ! do nothing
      else
       print *,"Tsat out of range 1"
       print *,"Tsat= ",Tsat
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

      call containing_cell(bfact,dx,xlo,lo,x,cell_index)

      ic=cell_index(1)
      jc=cell_index(2)
      kc=cell_index(SDIM)

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

       local_data_fab=>tempfab
       call safe_data(isten,jsten,ksten,comp,local_data_fab,local_data_out)
       T_sten(D_DECL(i1,j1,k1))=local_data_out

       do im_local=1,num_materials
        vofcomp=(im_local-1)*ngeom_recon+1
        local_data_fab=>recon
        call safe_data(isten,jsten,ksten,vofcomp,local_data_fab,local_data_out)
        local_VOF(im_local)=local_data_out
       enddo !im_local=1..num_materials

       call get_primary_material_VFRAC(local_VOF, &
         im_primary_sten(D_DECL(i1,j1,k1)))

       VF_sten(D_DECL(i1,j1,k1))=local_VOF(im)

       vofcomp=(im-1)*ngeom_recon+1
       local_data_fab=>recon
       do dir=1,SDIM
        call safe_data(isten,jsten,ksten,vofcomp+dir, &
          local_data_fab,local_data_out)
        XC_sten(D_DECL(i1,j1,k1),dir)=local_data_out+cencell(dir)
       enddo

       local_data_fab=>LS
       call safe_data(isten,jsten,ksten,im,local_data_fab,local_data_out)
       LS_sten(D_DECL(i1,j1,k1))=local_data_out

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
       im, &
       im_primary_sten, & 
       VF_sten, & 
       LS_sten, & 
       Tsat, &
       T_out)

      dest=T_out(1)

      return 
      end subroutine interpfabTEMP

       ! called from: mdot_from_Y_probe
       !              mdot_from_T_probe
       !              fort_ratemasschange
      subroutine probe_interpolation( &
       PROBE_PARMS, & !intent(in)
       T_I,Y_I, & !intent(in)
       POUT) !intent(out)
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE
 
      type(probe_parm_type), INTENT(in) :: PROBE_PARMS
      real(amrex_real), INTENT(in) :: T_I
      real(amrex_real), INTENT(in) :: Y_I
      type(probe_out_type), INTENT(out) :: POUT
      real(amrex_real) :: T_probe_no_constrain
      real(amrex_real) :: Y_probe_no_constrain
      integer :: iprobe
      real(amrex_real) xtarget_probe(SDIM)
      integer im_target_probe(2) ! source,dest
      integer im_primary_probe(2)
      integer im_target_probe_opp(2)
      integer Ycomp_probe(2)
      integer tcomp_probe(2)
      integer dencomp_probe(2)
      integer imls
      integer dir
      integer mtype
      real(amrex_real) LSPROBE(num_materials)
      real(amrex_real) F_tess(num_materials)

      if ((Y_I.ge.zero).and.(Y_I.le.one)) then
       ! do nothing
      else
       print *,"Y_I out of range"
       print *,"Y_I= ",Y_I
       stop
      endif
      if ((T_I.ge.zero).and.(T_I.le.1.0D+30)) then
       ! do nothing
      else
       print *,"T_I out of range in probe_interpolation"
       print *,"T_I= ",T_I
       print *,"Y_I= ",Y_I
       stop
      endif
       ! 0="im_primary_probe(iprobe).ne.im_target_probe(iprobe)"
       ! 1="im_primary_probe(iprobe).eq.im_target_probe(iprobe)"
       ! iprobe=1 source
       ! iprobe=2 dest
      POUT%interp_valid_flag(1)=0
      POUT%interp_valid_flag(2)=0

       ! tessellating volume fractions at xI.
      call interpfabVFRAC_tess( &
       PROBE_PARMS%tid, &
       PROBE_PARMS%bfact, &
       PROBE_PARMS%level, &
       PROBE_PARMS%finest_level, &
       PROBE_PARMS%dx, &
       PROBE_PARMS%xlo, &
       PROBE_PARMS%xI, &
       PROBE_PARMS%fablo, &
       PROBE_PARMS%fabhi, &
       PROBE_PARMS%recon, &
       F_tess)

      do iprobe=1,2 ! iprobe=1 source    iprobe=2 dest

       if (iprobe.eq.1) then ! source
        do dir=1,SDIM
         xtarget_probe(dir)=PROBE_PARMS%xsrc(dir)
        enddo
        im_target_probe(iprobe)=PROBE_PARMS%im_source
        im_target_probe_opp(iprobe)=PROBE_PARMS%im_dest
        tcomp_probe(iprobe)=PROBE_PARMS%tcomp_source
        Ycomp_probe(iprobe)=PROBE_PARMS%Ycomp_source
        dencomp_probe(iprobe)=PROBE_PARMS%dencomp_source
        POUT%dxprobe_target(iprobe)=PROBE_PARMS%dxprobe_source
       else if (iprobe.eq.2) then  ! dest
        do dir=1,SDIM
         xtarget_probe(dir)=PROBE_PARMS%xdst(dir)
        enddo
        im_target_probe(iprobe)=PROBE_PARMS%im_dest
        im_target_probe_opp(iprobe)=PROBE_PARMS%im_source
        tcomp_probe(iprobe)=PROBE_PARMS%tcomp_dest
        Ycomp_probe(iprobe)=PROBE_PARMS%Ycomp_dest
        dencomp_probe(iprobe)=PROBE_PARMS%dencomp_dest
        POUT%dxprobe_target(iprobe)=PROBE_PARMS%dxprobe_dest
       else
        print *,"iprobe invalid"
        stop
       endif
  
       POUT%vfrac_I(iprobe)=F_tess(im_target_probe(iprobe))

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
         dencomp_probe(iprobe), &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &        ! Fortran array box
         PROBE_PARMS%recon, &
         POUT%den_I_interp(iprobe))

        call interpfabFWEIGHT( &
         PROBE_PARMS%bfact, &
         PROBE_PARMS%level, &
         PROBE_PARMS%finest_level, &
         PROBE_PARMS%dx, &
         PROBE_PARMS%xlo, &
         xtarget_probe, &
         im_target_probe(iprobe), &
         dencomp_probe(iprobe), &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &        ! Fortran array box
         PROBE_PARMS%recon, &
         POUT%den_probe(iprobe))

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
        tcomp_probe(iprobe), &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%EOS, &
        PROBE_PARMS%LS, &
        PROBE_PARMS%recon, &
        POUT%T_probe(iprobe), & !T_probe constrained by T_I
        PROBE_PARMS%debugrate)

       if (POUT%T_probe(iprobe).ge.zero) then
        ! do nothing
       else
        print *,"T_probe went negative"
        print *,"T_probe ",POUT%T_probe(iprobe)
        stop
       endif

        ! the least squares interpolant is limited by the stencil values.
        ! interpfabFWEIGHT calls center_centroid_interchange with
        ! cc_flag=0 (centroid->target), tsat_flag=0,nsolve=1,Tsat=293
        ! (placeholder).
       call interpfabFWEIGHT( &
        PROBE_PARMS%bfact, &
        PROBE_PARMS%level, &
        PROBE_PARMS%finest_level, &
        PROBE_PARMS%dx, &
        PROBE_PARMS%xlo, &
        xtarget_probe, &
        im_target_probe(iprobe), &
        tcomp_probe(iprobe), &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%EOS, &   
        PROBE_PARMS%recon, &
        T_probe_no_constrain)

       if (T_probe_no_constrain.ge.zero) then
        ! do nothing
       else
        print *,"T_probe_no_constrain went negative"
        print *,"T_probe_no_constrain ",T_probe_no_constrain
        stop
       endif

       if (1.eq.0) then
        print *,"probe_constrain,iprobe,T_probe,T_probe_new ", &
         PROBE_PARMS%probe_constrain, &
         iprobe,POUT%T_probe(iprobe),T_probe_no_constrain
       endif

       if (PROBE_PARMS%probe_constrain.eq.1) then
        ! do nothing
       else if (PROBE_PARMS%probe_constrain.eq.0) then
        POUT%T_probe(iprobe)=T_probe_no_constrain
       else 
        print *,"PROBE_PARMS%probe_constrain invalid"
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
         Ycomp_probe(iprobe), &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &
         PROBE_PARMS%LS, &
         PROBE_PARMS%recon, &
         POUT%Y_probe(iprobe), &
         PROBE_PARMS%debugrate)

        if ((POUT%Y_probe(iprobe).ge.-EPS_8_4).and. &
            (POUT%Y_probe(iprobe).le.zero)) then
         POUT%Y_probe(iprobe)=zero
        else if ((POUT%Y_probe(iprobe).ge.zero).and. &
                 (POUT%Y_probe(iprobe).le.one)) then
         ! do nothing
        else if ((POUT%Y_probe(iprobe).ge.one).and. &
                 (POUT%Y_probe(iprobe).le.one+EPS_8_4)) then
         POUT%Y_probe(iprobe)=one
        else
         print *,"Y_probe out of bounds probe_interpolation"
         print *,"Y_probe ",POUT%Y_probe(iprobe)
         stop
        endif

        call interpfabFWEIGHT( &
         PROBE_PARMS%bfact, &
         PROBE_PARMS%level, &
         PROBE_PARMS%finest_level, &
         PROBE_PARMS%dx, &
         PROBE_PARMS%xlo, &
         xtarget_probe, &
         im_target_probe(iprobe), &
         Ycomp_probe(iprobe), &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &        ! Fortran array box
         PROBE_PARMS%recon, &
         Y_probe_no_constrain)

        if ((Y_probe_no_constrain.ge.-EPS_8_4).and. &
            (Y_probe_no_constrain.le.zero)) then
         Y_probe_no_constrain=zero
        else if ((Y_probe_no_constrain.ge.zero).and. &
                 (Y_probe_no_constrain.le.one)) then
         ! do nothing
        else if ((Y_probe_no_constrain.ge.one).and. &
                 (Y_probe_no_constrain.le.one+EPS_8_4)) then
         Y_probe_no_constrain=one
        else
         print *,"Y_probe_no_constrain out of bounds probe_interpolation"
         print *,"Y_probe_no_constrain ",Y_probe_no_constrain
         stop
        endif

        if (PROBE_PARMS%probe_constrain.eq.1) then
         ! do nothing
        else if (PROBE_PARMS%probe_constrain.eq.0) then
         POUT%Y_probe(iprobe)=Y_probe_no_constrain
        else 
         print *,"PROBE_PARMS%probe_constrain invalid"
         stop
        endif

         ! cc_flag=0 centroid->target
         ! tsat_flag=0 do not use TSAT
        call interpfabFWEIGHT( &
         PROBE_PARMS%bfact, &
         PROBE_PARMS%level, &
         PROBE_PARMS%finest_level, &
         PROBE_PARMS%dx, &
         PROBE_PARMS%xlo, &
         PROBE_PARMS%xI, &
         im_target_probe(iprobe), &
         Ycomp_probe(iprobe), &
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%EOS, &
         PROBE_PARMS%recon, &
         POUT%Y_I_interp(iprobe))
       else if (Ycomp_probe(iprobe).eq.0) then
        POUT%Y_probe(iprobe)=one
        POUT%Y_I_interp(iprobe)=one
       else
        print *,"Ycomp_probe invalid"
        stop
       endif

       do imls=1,num_materials
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
         PROBE_PARMS%fablo, &
         PROBE_PARMS%fabhi, &
         PROBE_PARMS%LS, &
         LSPROBE(imls))
       enddo ! imls=1..num_materials

       call get_primary_material(LSPROBE,im_primary_probe(iprobe))

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

        POUT%interp_valid_flag(iprobe)=1

       else if (im_primary_probe(iprobe).ne. &
                im_target_probe(iprobe)) then

        POUT%interp_valid_flag(iprobe)=0

       else
        print *,"im_primary_probe or im_target_probe invalid"
        stop
       endif

        ! check for NaN
       call grad_probe_sanity( &
        PROBE_PARMS%xI, &
        xtarget_probe, &
        POUT%T_probe(iprobe), &
        T_I, &
        PROBE_PARMS%LL)

        ! check for NaN
       call grad_probe_sanity( &
        PROBE_PARMS%xI, &
        xtarget_probe, &
        POUT%Y_probe(iprobe), &
        Y_I, &
        PROBE_PARMS%LL)

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
        tcomp_probe(iprobe), &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%EOS, &
        PROBE_PARMS%recon, &
        POUT%T_probe_raw(iprobe)) ! T_probe_raw same as T_probe_no_constrain

       if (POUT%T_probe_raw(iprobe).ge.zero) then
        ! do nothing
       else
        print *,"T_probe_raw went negative"
        print *,"T_probe_raw ",POUT%T_probe_raw(iprobe)
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
        tcomp_probe(iprobe), &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%EOS, &
        PROBE_PARMS%recon, &
        POUT%T_I_interp(iprobe))

       call single_interpfab( &
        PROBE_PARMS%bfact, &
        PROBE_PARMS%level, &
        PROBE_PARMS%finest_level, &
        PROBE_PARMS%dx, &
        PROBE_PARMS%xlo, &
        PROBE_PARMS%xI, &
        PROBE_PARMS%fablo, &
        PROBE_PARMS%fabhi, &
        PROBE_PARMS%pres, &
        POUT%pres_I_interp(iprobe))

       ! local_freezing_model=0 (sharp interface stefan model)
       ! local_freezing_model=1 (source term model)
       ! local_freezing_model=2 (hydrate model)
       ! local_freezing_model=3 (wildfire)
       ! local_freezing_model=5 (evaporation/condensation)
       ! local_freezing_model=6 (evaporation/condensation Palmore)
       ! local_freezing_model=7 (cavitation (under construction))
       if &
        (is_valid_freezing_modelF(PROBE_PARMS%local_freezing_model).eq.1) then
        ! do nothing
       else
        print *,"PROBE_PARMS%local_freezing_model invalid"
        stop
       endif

      enddo ! iprobe=1..2

      return
      end subroutine probe_interpolation


      subroutine apply_TI_limiter( &
       TI_min,TI_max, &
       PROBE_PARMS, &
       T_I,Y_I, &
       POUT)
      use global_utility_module
      use MOF_routines_module

      IMPLICIT NONE

      real(amrex_real), INTENT(inout) :: TI_min,TI_max 
      type(probe_parm_type), INTENT(in) :: PROBE_PARMS
      real(amrex_real), INTENT(inout) :: T_I
      real(amrex_real), INTENT(inout) :: Y_I
      type(probe_out_type), INTENT(in) :: POUT
      real(amrex_real) :: local_TI_min,local_TI_max 
      integer :: iprobe

      if ((Y_I.ge.zero).and.(Y_I.le.one)) then
       ! do nothing
      else
       print *,"Y_I out of range"
       print *,"Y_I= ",Y_I
       stop
      endif
      if ((T_I.ge.zero).and.(T_I.le.1.0D+30)) then
       ! do nothing
      else
       print *,"T_I out of range in apply_TI_limiter"
       print *,"T_I= ",T_I
       print *,"Y_I= ",Y_I
       stop
      endif

      local_TI_min=POUT%T_probe(1)
      local_TI_max=POUT%T_probe(1)

      if (POUT%T_probe(2).lt.local_TI_min) then
       local_TI_min=POUT%T_probe(2)
      endif
      if (POUT%T_probe(2).gt.local_TI_max) then
       local_TI_max=POUT%T_probe(2)
      endif

      do iprobe=1,2 ! iprobe=1 source    iprobe=2 dest

       if (POUT%T_I_interp(iprobe).lt.local_TI_min) then
        local_TI_min=POUT%T_I_interp(iprobe)
       endif
       if (POUT%T_I_interp(iprobe).gt.local_TI_max) then
        local_TI_max=POUT%T_I_interp(iprobe)
       endif

      enddo ! iprobe=1,2

      if (local_TI_min.gt.TI_min) then
       TI_min=local_TI_min
      endif
      if (local_TI_max.lt.TI_max) then
       TI_max=local_TI_max
      endif

      if (T_I.lt.TI_min) then
       T_I=TI_min
      endif
      if (T_I.gt.TI_max) then
       T_I=TI_max
      endif
      if (TI_min.le.TI_max) then
       ! do nothing
      else
       print *,"TI_min>TI_max"
       stop
      endif

      return
      end subroutine apply_TI_limiter

      subroutine mdot_from_Y_probe( &
       prescribed_mdot, &
       probe_ok, &
       TSAT_Y_PARMS, &
       POUT, &
       Y_gamma,T_gamma, &
       mdotY_top, &
       mdotY_bot, &
       mdotY, &
       Y_PROBE_VAPOR)
      use global_utility_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: prescribed_mdot
      integer, INTENT(in) :: probe_ok
      type(TSAT_MASS_FRAC_parm_type), INTENT(in) :: TSAT_Y_PARMS
      type(probe_out_type), INTENT(inout) :: POUT
      real(amrex_real), INTENT(inout) :: Y_gamma
      real(amrex_real), INTENT(inout) :: T_gamma
      real(amrex_real), INTENT(out) :: mdotY_top,mdotY_bot,mdotY
      real(amrex_real), INTENT(out) :: Y_PROBE_VAPOR
      real(amrex_real) D_MASS
      real(amrex_real) LL
      integer iprobe_vapor
      real(amrex_real) den_G
      real(amrex_real) YI_min
      integer Kassemi_flag
      real(amrex_real) Pgamma
      real(amrex_real) density_probe
      real(amrex_real) Tvapor_probe
      real(amrex_real) Pvapor_probe
      real(amrex_real) internal_energy
      integer im_probe
      integer imattype
      real(amrex_real) massfrac_parm(num_species_var+1)
      integer ispec

      Y_PROBE_VAPOR=one

      YI_min=TSAT_Y_PARMS%YI_min

      if ((Y_gamma.ge.YI_min).and. &
          (YI_min.ge.zero).and. &
          (YI_min.le.one).and. &
          (Y_gamma.le.one)) then

       !iprobe=1 source
       !iprobe=2 dest
       if (probe_ok.eq.1) then ! probes do not depend on TI
        ! do nothing
       else if (probe_ok.eq.0) then
        call probe_interpolation( &
         TSAT_Y_PARMS%PROBE_PARMS, &
         T_gamma,Y_gamma, &
         POUT)
       else
        print *,"probe_ok invalid"
        stop
       endif

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

       if ((POUT%dxprobe_target(1).gt.zero).and. &
           (POUT%dxprobe_target(2).gt.zero)) then

        Kassemi_flag=TSAT_Y_PARMS%Tanasawa_or_Schrage_or_Kassemi

        if (Kassemi_flag.eq.0) then ! Palmore/Desjardins

         D_MASS=TSAT_Y_PARMS%D_MASS
         den_G=TSAT_Y_PARMS%den_G
         if ((D_MASS.gt.zero).and.(den_G.gt.zero)) then

          Y_PROBE_VAPOR=POUT%Y_probe(iprobe_vapor)

          if (prescribed_mdot.eq.zero) then
           ! do nothing
          else if (prescribed_mdot.gt.zero) then
           !mdot=(Y-Yprobe)*rho*D/((1-Y)*dx)
           !(1-Y)dx mdot=(Y-Yprobe)*rho*D
           !Yprobe*rho*D+dx mdot=Y(rho*D+dx*mdot)
           !Y=(dx mdot + Yprobe*rho*D)/(dx mdot+rho*D)
           Y_gamma= &
             (POUT%dxprobe_target(iprobe_vapor)*prescribed_mdot+ &
             Y_PROBE_VAPOR*den_G*D_MASS)/ &
             (POUT%dxprobe_target(iprobe_vapor)*prescribed_mdot+ &
             den_G*D_MASS)
           if ((Y_gamma.ge.zero).and.(Y_gamma.le.one)) then
            !do nothing
           else
            print *,"Y_gamma invalid: ",Y_gamma
            stop
           endif
          else
           print *,"prescribed_mdot invalid: ",prescribed_mdot
           stop
          endif

          mdotY_top=Y_gamma-Y_PROBE_VAPOR
          mdotY_top=den_G*D_MASS*mdotY_top/ &
                 POUT%dxprobe_target(iprobe_vapor)
          mdotY_bot=one-Y_gamma
          if (mdotY_bot.lt.zero) then
           print *,"Y_gamma cannot exceed 1; try increasing TSAT"
           print *,"mdotY_top= ",mdotY_top
           print *,"mdotY_bot= ",mdotY_bot
           print *,"den_G= ",den_G
           print *,"D_mass= ",D_mass
           print *,"Y_gamma= ",Y_gamma
           print *,"iprobe_vapor=",iprobe_vapor
           print *,"Y_probe(iprobe_vapor)=",POUT%Y_probe(iprobe_vapor)
           print *,"dxprobe_target(iprobe_vapor)=", &
                  POUT%dxprobe_target(iprobe_vapor)
           stop
          else if (mdotY_bot.eq.zero) then
           mdotY=zero
          else if ((mdotY_bot.gt.zero).and. &
                   (mdotY_bot.le.one)) then
           mdotY=mdotY_top/mdotY_bot
          else
           print *,"mdotY_bot invalid: ",mdotY_bot
           stop
          endif
         else
          print *,"D_MASS or den_G invalid"
          print *,"Kassemi_flag=",Kassemi_flag
          print *,"D_MASS=",D_MASS
          print *,"den_G=",den_G
          stop
         endif

        else if (Kassemi_flag.eq.3) then ! Kassemi model

         if (prescribed_mdot.eq.zero) then

          ! Pgamma_Clausius_Clapyron is declared in: GLOBALUTIL.F90
          call Pgamma_Clausius_Clapyron( &
           Pgamma, & !intent(out)
           TSAT_Y_PARMS%reference_pressure, & !intent(in)
           T_gamma, & !intent(in)
           TSAT_Y_PARMS%Clausius_Clapyron_Tsat, & !intent(in)
           LL, & !intent(in)
           TSAT_Y_PARMS%universal_gas_constant_R, & !intent(in)
           TSAT_Y_PARMS%molar_mass_vapor) !intent(in)

          density_probe=POUT%den_I_interp(iprobe_vapor)
          Tvapor_probe=POUT%T_probe(iprobe_vapor)

          if (LL.gt.zero) then  ! evaporation (destination=vapor)
           im_probe=TSAT_Y_PARMS%PROBE_PARMS%im_dest
          else if (LL.lt.zero) then ! condensation (source=vapor)
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
           T_gamma,internal_energy,imattype,im_probe)
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

         else if (prescribed_mdot.gt.zero) then
       
          Y_gamma=one
          mdotY=prescribed_mdot

         else
          print *,"prescribed_mdot invalid: ",prescribed_mdot
          stop
         endif

         mdotY_top=mdotY
         mdotY_bot=one

        else
         print *,"Kassemi_flag invalid in mdot_from_Y_probe"
         stop
        endif

       else
        print *,"dxprobe_target invalid"
        stop
       endif

      else
       print *,"Y_gamma or YI_min invalid"
       stop
      endif

      end subroutine mdot_from_Y_probe


      subroutine mdot_from_T_probe( &
       prescribed_mdot, &
       probe_ok, &
       TSAT_Y_PARMS, &
       POUT, &
       T_gamma,Y_gamma, &
       mdotT, &
       TEMP_PROBE_source, &
       TEMP_PROBE_dest)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: prescribed_mdot
      integer, INTENT(in) :: probe_ok
      type(TSAT_MASS_FRAC_parm_type), INTENT(in) :: TSAT_Y_PARMS
      type(probe_out_type), INTENT(inout) :: POUT
      real(amrex_real), INTENT(inout) :: T_gamma
      real(amrex_real), INTENT(inout) :: Y_gamma
      real(amrex_real), INTENT(out) :: mdotT
      real(amrex_real), INTENT(out) :: TEMP_PROBE_source
      real(amrex_real), INTENT(out) :: TEMP_PROBE_dest
      real(amrex_real) D_MASS
      integer Kassemi_flag
      real(amrex_real) LL
      integer iprobe_vapor
      real(amrex_real) den_G
      real(amrex_real) YI_min
      real(amrex_real) wt(2)
      real(amrex_real) probediff(2)
      integer iprobe

      YI_min=TSAT_Y_PARMS%YI_min

      if ((Y_gamma.ge.YI_min).and. &
          (YI_min.ge.zero).and. &
          (YI_min.le.one).and. &
          (Y_gamma.le.one)) then

       Kassemi_flag=TSAT_Y_PARMS%Tanasawa_or_Schrage_or_Kassemi
       D_MASS=TSAT_Y_PARMS%D_MASS
       den_G=TSAT_Y_PARMS%den_G
       if (Kassemi_flag.eq.0) then
        if (D_MASS.gt.zero) then
         ! do nothing
        else 
         print *,"D_MASS invalid in mdot_from_T_probe"
         stop
        endif
       else if (Kassemi_flag.eq.3) then
        if (D_MASS.eq.zero) then
         ! do nothing
        else
         print *,"D_MASS invalid in mdot_from_T_probe"
         stop
        endif
       else
        print *,"Kassemi_flag invalid in mdot_from_T_probe"
        stop
       endif

       if ((D_MASS.ge.zero).and.(den_G.gt.zero)) then

        if (probe_ok.eq.1) then ! probes do not depend on TI
         ! do nothing
        else if (probe_ok.eq.0) then
         !iprobe=1 source
         !iprobe=2 dest
         call probe_interpolation( &
          TSAT_Y_PARMS%PROBE_PARMS, &
          T_gamma,Y_gamma, &
          POUT)
        else
         print *,"probe_ok invalid"
         stop
        endif

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

        if ((POUT%dxprobe_target(1).gt.zero).and. &
            (POUT%dxprobe_target(2).gt.zero)) then
        
         do iprobe=1,2
          wt(iprobe)=TSAT_Y_PARMS%thermal_k(iprobe)/ &
                 (abs(LL)*POUT%dxprobe_target(iprobe))
          if (wt(iprobe).ge.zero) then
           ! do nothing
          else
           print *,"wt(iprobe) invalid"
           stop
          endif
         enddo ! iprobe=1,2

         TEMP_PROBE_source=POUT%T_probe(1)
         TEMP_PROBE_dest=POUT%T_probe(2)

         if (wt(1)+wt(2).gt.zero) then

          if (prescribed_mdot.eq.zero) then
           ! do nothing
          else if (prescribed_mdot.gt.zero) then
           !mdot=w1 (TS-T)+w2 (TD-T)=w1 TS + w2 TD-(w1+w2)T
           !T=(w1 TS + w2 TD - mdot)/(w1+w2)
           T_gamma=(wt(1)*TEMP_PROBE_source+wt(2)*TEMP_PROBE_dest- &
                    prescribed_mdot)/(wt(1)+wt(2))

           if (T_gamma.gt.zero) then
            !do nothing
           else
            print *,"T_gamma invalid: ",T_gamma
            print *,"prescribed_mdot: ",prescribed_mdot
            stop
           endif

          else
           print *,"prescribed_mdot invalid: ",prescribed_mdot
           stop
          endif

          probediff(1)=TEMP_PROBE_source-T_gamma
          probediff(2)=TEMP_PROBE_dest-T_gamma
          mdotT=wt(1)*probediff(1)+wt(2)*probediff(2)
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
        print *,"D_MASS or den_G invalid"
        stop
       endif
      else
       print *,"Y_gamma or YI_min invalid"
       stop
      endif

      end subroutine mdot_from_T_probe


      subroutine mdot_diff_func( &
       prescribed_mdot, &
       probe_ok, &
       TSAT_Y_PARMS, &
       POUT, &
       Y_gamma, &
       T_gamma, &
       mdot_diff)
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: prescribed_mdot
      integer, INTENT(in) :: probe_ok
      type(TSAT_MASS_FRAC_parm_type), INTENT(in) :: TSAT_Y_PARMS
      type(probe_out_type), INTENT(inout) :: POUT
      real(amrex_real), INTENT(inout) :: T_gamma
      real(amrex_real), INTENT(inout) :: Y_gamma
      real(amrex_real), INTENT(out) :: mdot_diff
      real(amrex_real) mdotT
      real(amrex_real) mdotY_top,mdotY_bot,mdotY
      integer Kassemi_flag
      real(amrex_real) TEMP_PROBE_source
      real(amrex_real) TEMP_PROBE_dest
      real(amrex_real) Y_PROBE_VAPOR

      Kassemi_flag=TSAT_Y_PARMS%Tanasawa_or_Schrage_or_Kassemi
      if (Kassemi_flag.eq.0) then
       if (Y_gamma.eq.one) then
        print *,"expecting Y_gamma lt one in mdot_diff_func"
        print *,"Y_gamma= ",Y_gamma
        print *,"probe_ok= ",probe_ok
        print *,"Clausius_Clapyron_Tsat= ",TSAT_Y_PARMS%Clausius_Clapyron_Tsat
        print *,"TSAT_base= ",TSAT_Y_PARMS%TSAT_base
        print *,"YI_min= ",TSAT_Y_PARMS%YI_min
        print *,"den_G= ",TSAT_Y_PARMS%den_G
        print *,"molar_mass_ambient= ",TSAT_Y_PARMS%molar_mass_ambient
        print *,"molar_mass_vapor= ",TSAT_Y_PARMS%molar_mass_vapor
        print *,"POUT%Y_probe(1)= ",POUT%Y_probe(1)
        print *,"POUT%Y_probe(2)= ",POUT%Y_probe(2)
        print *,"TSAT_Y_PARMS%PROBE_PARMS%i= ", &
                TSAT_Y_PARMS%PROBE_PARMS%i
        print *,"TSAT_Y_PARMS%PROBE_PARMS%j= ", &
                TSAT_Y_PARMS%PROBE_PARMS%j
        print *,"TSAT_Y_PARMS%PROBE_PARMS%k= ", &
                TSAT_Y_PARMS%PROBE_PARMS%k
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
       prescribed_mdot, &
       probe_ok, &
       TSAT_Y_PARMS, &
       POUT, &
       T_gamma,Y_gamma, &
       mdotT, &
       TEMP_PROBE_source, &
       TEMP_PROBE_dest)

      call mdot_from_Y_probe( &
       prescribed_mdot, &
       probe_ok, &
       TSAT_Y_PARMS, &
       POUT, &
       Y_gamma,T_gamma, &
       mdotY_top, &
       mdotY_bot, &
       mdotY, &
       Y_PROBE_VAPOR)

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

      type(TSAT_MASS_FRAC_parm_type), INTENT(in) :: TSAT_Y_PARMS
      real(amrex_real), INTENT(out) :: Y_I_MIN
      real(amrex_real) :: X_I_MIN
      real(amrex_real) :: LL,R
      real(amrex_real) :: Clausius_Clapyron_Tsat
      real(amrex_real) :: TSAT_base
      real(amrex_real) :: T_I_MAX,WA,WV
      real(amrex_real) YI_min_old

      YI_min_old=TSAT_Y_PARMS%YI_min

      WA=TSAT_Y_PARMS%molar_mass_ambient
      WV=TSAT_Y_PARMS%molar_mass_vapor
      R=TSAT_Y_PARMS%universal_gas_constant_R
      LL=TSAT_Y_PARMS%PROBE_PARMS%LL

      if ((YI_min_old.ge.zero).and.(YI_min_old.le.one)) then

       if ((LL.gt.zero).or.(LL.lt.zero)) then

        TSAT_base=TSAT_Y_PARMS%TSAT_base
        Clausius_Clapyron_Tsat=TSAT_Y_PARMS%Clausius_Clapyron_Tsat
        T_I_MAX=TSAT_Y_PARMS%TI_max

        if ((TSAT_base.gt.zero).and. &
            (T_I_MAX.ge.TSAT_base)) then

         if ((WA.gt.zero).and.(WV.gt.zero).and.(R.gt.zero)) then

          if (LL.gt.zero) then
           Y_I_MIN=YI_min_old
          else if (LL.lt.zero) then

             ! X=exp(-(L*WV/R)*(one/Tgamma-one/TSAT))
           call X_from_Tgamma(X_I_MIN,T_I_MAX,Clausius_Clapyron_Tsat,LL,R,WV)

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

      module mass_transfer_cpp_module
      use amrex_fort_module, only : amrex_real
      implicit none

      integer, PARAMETER :: MAX_TI_YI_logfile=100
      integer, PARAMETER :: MAX_TI_YI_trials=20

      contains

        ! this is for unsplit advection: for phase change.
        ! Called from NavierStokes.cpp: 
        !   NavierStokes::level_phase_change_convert
      subroutine fort_nodedisplace( &
       nburning, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       velbc, &
       dt, &
       unode,DIMS(unode), &
       ucell,DIMS(ucell), &
       oldLS,DIMS(oldLS), &
       xlo,dx, &
       level,finest_level) &
      bind(c,name='fort_nodedisplace')

      use probcommon_module
      use global_utility_module
      use mass_transfer_module

      IMPLICIT NONE
       
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: nburning
      integer ncomp_per_burning
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: dt

      integer, INTENT(in) :: DIMDEC(unode)
      integer, INTENT(in) :: DIMDEC(ucell)
      integer, INTENT(in) :: DIMDEC(oldLS)
     
      real(amrex_real), target, INTENT(out) :: &
              unode(DIMV(unode),2*num_interfaces*SDIM) 
      real(amrex_real), pointer :: unode_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(in) ::  ucell(DIMV(ucell),nburning) 
      real(amrex_real), pointer :: ucell_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: &
              oldLS(DIMV(oldLS),num_materials*(1+SDIM))
      real(amrex_real), pointer :: oldLS_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: velbc(SDIM,2,SDIM)

      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      integer, parameter :: nhalf=3

      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
     
      integer i,j,k,klosten,khisten,i1,j1,k1
      integer dir,itencrit

      real(amrex_real) delta,velnd,totalwt
      integer scomp,tag_local
      integer ireverse
      integer sign_reverse
      integer rigid_in_stencil
      integer im
      integer im_primary
      real(amrex_real) LS_local(num_materials)

      unode_ptr=>unode
      ucell_ptr=>ucell
      oldLS_ptr=>oldLS

      if (bfact.lt.1) then
       print *,"bfact invalid117"
       stop
      endif

      ncomp_per_burning=EXTRAP_PER_BURNING
      if (ncomp_per_burning.eq.AMREX_SPACEDIM) then
       ! do nothing
      else
       print *,"expecting ncomp_per_burning==sdim"
       stop
      endif

      if (nburning.eq.EXTRAP_NCOMP_BURNING) then
       ! do nothing
      else
       print *,"nburning invalid"
       stop
      endif

      if (nburning.eq.(1+AMREX_SPACEDIM)*num_interfaces) then
       ! do nothing
      else
       print *,"expecting nburning=(1+sdim)*num_interfaces"
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

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"levelrz invalid node displace"
       stop
      endif
       
      if (ngrow_distance.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance invalid"
       stop
      endif

      call checkbound_array(fablo,fabhi,unode_ptr,1,-1)
      call checkbound_array(fablo,fabhi,ucell_ptr,1,-1)
      call checkbound_array(fablo,fabhi,oldLS_ptr,ngrow_distance,-1)

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      call growntileboxNODE(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridstenND_level(xstenND,i,j,k,level,nhalf)

       rigid_in_stencil=0

       do k1=klosten,khisten
       do j1=0,1
       do i1=0,1
        do im=1,num_materials
         LS_local(im)=oldLS(D_DECL(i-i1,j-j1,k-k1),im)
        enddo
        call get_primary_material(LS_local,im_primary)
        if (is_rigid(im_primary).eq.0) then
         ! do nothing
        else if (is_rigid(im_primary).eq.1) then
         rigid_in_stencil=1
        else
         print *,"is_rigid invalid MASS_TRANSFER_3D.F90"
         stop
        endif
       enddo
       enddo
       enddo

       do itencrit=1,num_interfaces
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
         do k1=klosten,khisten
         do j1=0,1
         do i1=0,1

          tag_local=NINT(ucell(D_DECL(i-i1,j-j1,k-k1),itencrit))
          if (tag_local.eq.0) then
           ! do nothing
          else if ((sign_reverse*tag_local.eq.1).or. &
                   (sign_reverse*tag_local.eq.2)) then
           scomp=num_interfaces+(itencrit-1)*ncomp_per_burning+dir
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
         if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
             (levelrz.eq.COORDSYS_RZ)) then
          ! do nothing
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
          if (dir.eq.2) then
           delta=delta/xstenND(0,1)
          endif
         else
          print *,"levelrz invalid node displace 2"
          stop
         endif

         if (levelrz.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if ((levelrz.eq.COORDSYS_RZ).or. &
                  (levelrz.eq.COORDSYS_CYLINDRICAL)) then
          if (dir.eq.1) then
           if (abs(xstenND(0,dir)).le.EPS_8_4*dx(dir)) then
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

         if (rigid_in_stencil.eq.1) then
          delta=zero
         else if (rigid_in_stencil.eq.0) then
          ! do nothing
         else
          print *,"rigid_in_stencil invalid"
          stop
         endif

         scomp=(itencrit+ireverse*num_interfaces-1)*ncomp_per_burning+dir
         unode(D_DECL(i,j,k),scomp)=delta
        enddo ! dir=1..sdim

       enddo ! ireverse=0,1
       enddo ! itencrit=1,num_interfaces

      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine fort_nodedisplace

      subroutine fort_apply_reaction( &
       tid, &
       level,finest_level, &
       nstate, &
       speciesreactionrate, &
       rigid_fraction_id, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       initialize_flag, &
       dt, &
       time, &
       maskcov,DIMS(maskcov), &
       LSnew,DIMS(LSnew), &
       snew,DIMS(snew)) &
      bind(c,name='fort_apply_reaction')

      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nstate
      real(amrex_real), INTENT(in) :: &
              speciesreactionrate(num_species_var*num_materials)
      integer, INTENT(in) :: rigid_fraction_id(num_materials)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in),target :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in),target :: xlo(SDIM)
      real(amrex_real), INTENT(in),target :: dx(SDIM)
      integer, INTENT(in) :: initialize_flag
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(LSnew)
      integer, INTENT(in) :: DIMDEC(snew)
      real(amrex_real), target, INTENT(in) :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), target, INTENT(out) :: LSnew(DIMV(LSnew),num_materials)
      real(amrex_real), pointer :: LSnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(out) :: snew(DIMV(snew),nstate)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)
      integer i,j,k
      integer im,im_inner,im_melt,im_opp,im_primary,ispec
      integer ireverse,iten
      real(amrex_real) LL
      real(amrex_real) LSMELT
      integer local_mask
      real(amrex_real) local_LS(num_materials)
      real(amrex_real) local_VOF(num_materials)
      real(amrex_real) local_DEN(num_materials)
      real(amrex_real) local_MASS(num_materials)
      real(amrex_real) species_vfrac_sum
      real(amrex_real) species_mass_sum
      real(amrex_real) species_avg
      real(amrex_real) speciesconst_avg
      real(amrex_real) species_base
      real(amrex_real) local_rate
      integer vofcomp
      integer spec_comp
      integer dencomp
      real(amrex_real) spec_old,spec_new
      real(amrex_real) local_cutoff
      real(amrex_real) species_scale
      real(amrex_real), PARAMETER :: species_max=1.0d0
      real(amrex_real), PARAMETER :: MUSHY_THICK=2.0d0

      LSnew_ptr=>LSnew
      snew_ptr=>snew
      maskcov_ptr=>maskcov

      if ((tid.lt.0).or. &
          (tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in reactionrate: ",level,finest_level
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif

       !called from NavierStokes::veldiffuseALL
      if (initialize_flag.eq.0) then

       if (dt.gt.zero) then
        ! do nothing
       else
        print *,"dt invalid (fort_apply_reaction, init_flag==0): ",dt
        stop
       endif
       if (time.gt.zero) then
        ! do nothing
       else
        print *,"time invalid (fort_apply_reaction, init_flag==0)"
        stop
       endif

       !called from NavierStokes::initData()
      else if (initialize_flag.eq.1) then

       if (dt.eq.zero) then
        ! do nothing
       else
        print *,"dt invalid (fort_apply_reaction, init_flag==1) ",dt
        stop
       endif
       if (time.eq.zero) then
        ! do nothing
       else
        print *,"time invalid (fort_apply_reaction, init_flag==1) ",time
        stop
       endif

      else
       print *,"initialize_flag invalid"
       stop
      endif

      if (species_max.gt.zero) then
       !do nothing
      else
       print *,"species_max invalid"
       stop
      endif

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LSnew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       if (initialize_flag.eq.1) then
        local_mask=1
       else if (initialize_flag.eq.0) then
        local_mask=NINT(maskcov(D_DECL(i,j,k)))
       else
        print *,"initialize_flag invalid"
        stop
       endif

       if (local_mask.eq.1) then

        species_vfrac_sum=zero
        species_mass_sum=zero

        do im=1,num_materials

         vofcomp=STATECOMP_MOF+(im-1)*ngeom_raw+1
         dencomp=STATECOMP_STATES+(im-1)*num_state_material+1+ENUM_DENVAR
         local_LS(im)=LSnew(D_DECL(i,j,k),im)
         local_VOF(im)=snew(D_DECL(i,j,k),vofcomp)
         local_DEN(im)=snew(D_DECL(i,j,k),dencomp)
         local_MASS(im)=local_VOF(im)*local_DEN(im)

         if (local_DEN(im).gt.zero) then
          ! do nothing
         else
          print *,"local_DEN(im) invalid: ",local_DEN(im)
          stop
         endif

         if ((local_VOF(im).ge.-EPS1).and. &
             (local_VOF(im).le.VOFTOL)) then
          local_VOF(im)=zero
         else if ((local_VOF(im).ge.one-VOFTOL).and. &
                  (local_VOF(im).le.one+EPS1)) then
          local_VOF(im)=one
         else if ((local_VOF(im).ge.zero).and. &
                  (local_VOF(im).le.one)) then
          ! do nothing
         else
          print *,"local_VOF invalid:",im,local_VOF(im)
          stop
         endif
         if (is_rigid(im).eq.0) then
          species_vfrac_sum=species_vfrac_sum+local_VOF(im)
          species_mass_sum=species_mass_sum+local_MASS(im)
         else if (is_rigid(im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif

        enddo !im=1..num_materials

        call get_primary_material(local_LS,im_primary)

        if (abs(species_vfrac_sum-one).le.EPS3) then
         ! do nothing
        else
         print *,"species_vfrac_sum invalid: ",species_vfrac_sum
         stop
        endif

        if (species_mass_sum.gt.zero) then
         ! do nothing
        else
         print *,"species_mass_sum invalid: ",species_mass_sum
         stop
        endif

        do ispec=1,num_species_var

         species_avg=zero
         speciesconst_avg=zero
         do im=1,num_materials
          local_rate=speciesreactionrate((ispec-1)*num_materials+im)

          spec_comp=STATECOMP_STATES+(im-1)*num_state_material+ &
               ENUM_SPECIESVAR+ispec

          species_scale=fort_speciesconst(im+num_materials*(ispec-1))
          speciesconst_avg=speciesconst_avg+species_scale

          if (species_scale.eq.zero) then
           species_scale=one
          else if ((species_scale.gt.zero).and. &
                   (species_scale.le.species_max)) then
           !do nothing
          else
           print *,"species_scale invalid: ",species_scale
           stop
          endif

          if (local_rate.ge.zero) then
           ! Y'=r(species_max-Y)
           ! Ynew=Yold+dt * r * (species_max-Ynew)
           ! Ynew=(Yold+species_max*r*dt)/(1+r*dt)
           spec_old=snew(D_DECL(i,j,k),spec_comp)
           if (abs(spec_old).le.EPS3*species_scale) then
            spec_old=zero
           else if (abs(spec_old-species_max).le.EPS3) then
            spec_old=species_max
           else if ((spec_old.ge.zero).and. &
                    (spec_old.le.species_max)) then
            ! do nothing
           else
            print *,"spec_old invalid: ",spec_old
            stop
           endif
           spec_new=(spec_old+species_max*local_rate*dt)/(one+local_rate*dt)

           if (abs(spec_new).le.EPS3*species_scale) then
            spec_new=zero
           else if (abs(spec_new-species_max).le.EPS3) then
            spec_new=species_max
           else if ((spec_new.ge.zero).and. &
                    (spec_new.le.species_max)) then
            ! do nothing
           else
            print *,"spec_new invalid: ",spec_new
            stop
           endif
           snew(D_DECL(i,j,k),spec_comp)=spec_new
          else
           print *,"local_rate invalid: ",local_rate
           stop
          endif

          if (is_rigid(im).eq.0) then
           species_avg=species_avg+spec_new*local_MASS(im)
          else if (is_rigid(im).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif
         enddo ! im=1,num_materials

         speciesconst_avg=speciesconst_avg/num_materials
         species_scale=speciesconst_avg
       
         if (species_scale.eq.zero) then
          species_scale=one
         else if ((species_scale.gt.zero).and. &
                  (species_scale.le.species_max)) then
          !do nothing
         else
          print *,"species_scale invalid (avg): ",species_scale
          stop
         endif

         if ((species_vfrac_sum.gt.zero).and. &
             (species_mass_sum.gt.zero)) then
          species_avg=species_avg/species_mass_sum
         else
          print *,"species_vfrac_sum or species_mass_sum invalid"
          print *,"species_vfrac_sum: ",species_vfrac_sum
          print *,"species_mass_sum: ",species_mass_sum
          stop
         endif
         if ((species_avg.ge.zero).and. &
             (species_avg.le.EPS3*species_scale)) then
          species_avg=zero
         else if (abs(species_avg-species_max).le.EPS3) then
          species_avg=species_max
         else if ((species_avg.gt.zero).and. &
                  (species_avg.lt.species_max))  then
          ! do nothing
         else
          print *,"species_avg invalid: ",species_avg
          stop
         endif

         do im=1,num_materials
          spec_comp=STATECOMP_STATES+(im-1)*num_state_material+ &
                ENUM_SPECIESVAR+ispec
          if (is_rigid(im).eq.0) then
           if (local_VOF(im).eq.zero) then
            snew(D_DECL(i,j,k),spec_comp)=species_avg
           else if ((local_VOF(im).gt.zero).and. &
                    (local_VOF(im).le.one)) then
            ! do nothing
           else
            print *,"local_VOF invalid"
            print *,"im=",im
            print *,"local_VOF(im)=",local_VOF(im)
            stop
           endif
          else if (is_rigid(im).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif

         enddo ! im=1,num_materials

        enddo ! ispec=1,num_species_var

        do im=1,num_materials
         if (is_ice(im).eq.1) then
          ispec=rigid_fraction_id(im)
          if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
           !do nothing 
          else
           print *,"ispec invalid"
           stop
          endif

          species_base=fort_speciesconst((ispec-1)*num_materials+im)

          if ((species_base.gt.zero).and. &
              (species_base.le.one)) then

           spec_comp=STATECOMP_STATES+(im-1)*num_state_material+ &
               ENUM_SPECIESVAR+ispec

           spec_old=snew(D_DECL(i,j,k),spec_comp)
           spec_new=spec_old

           if ((spec_old.ge.zero).and. &
               (spec_old.lt.species_base)) then
            spec_old=species_base
            spec_new=spec_old
           else if ((spec_old.ge.species_base).and. &
                    (spec_old.le.species_max)) then

            im_melt=0

            do im_opp=1,num_materials
             if (im_opp.ne.im) then
              do ireverse=0,1
               call get_iten(im,im_opp,iten)
               LL=get_user_latent_heat(iten+ireverse*num_interfaces, &
                                       room_temperature,1)
               if (LL.eq.zero) then
                ! do nothing
               else if (LL.ne.zero) then
                im_melt=im_opp
               else
                print *,"LL invalid"
                stop
               endif
              enddo !ireverse=0,1
             else if (im_opp.eq.im) then
              ! do nothing
             else
              print *,"im_opp invalid"
              stop
             endif
            enddo !im_opp=1,num_materials

            if ((im_melt.ge.1).and.(im_melt.le.num_materials)) then

             if (is_rigid(im_primary).eq.1) then

              spec_new=one

              if (1.eq.0) then
               print *,"i,j,k,im,im_melt,im_primary,species_base ", &
                  i,j,k,im,im_melt,im_primary,species_base
              endif

             else if (is_rigid(im_primary).eq.0) then

              if (local_LS(im).lt.zero) then
               spec_new=species_base
              else if (local_LS(im).ge.zero) then
               if (local_LS(im).ge.MUSHY_THICK*dx(1)) then
                spec_new=species_max
               else if ((local_LS(im).le.MUSHY_THICK*dx(1)).and. &
                        (local_LS(im).ge.zero)) then
                if (local_LS(im_melt).le.-MUSHY_THICK*dx(1)) then
                 spec_new=species_max
                else if (local_LS(im_melt).ge.zero) then
                 spec_new=species_base
                else if ((local_LS(im_melt).ge.-MUSHY_THICK*dx(1)).and. &
                         (local_LS(im_melt).le.zero)) then
                 local_cutoff=half*MUSHY_THICK*dx(1)
                 LSMELT=abs(local_LS(im_melt))-local_cutoff
                 spec_new=species_base+ &
                    (one-species_base)*hs(LSMELT,local_cutoff)
                else
                 print *,"local_LS(im_melt) invalid: ",local_LS(im_melt)
                 stop
                endif
               else
                print *,"local_LS(im) invalid(1): ",local_LS(im)
                stop
               endif
              else
               print *,"local_LS(im) invalid(2): ",local_LS(im)
               stop
              endif

             else
              print *,"is_rigid(im_primary) invalid"
              stop
             endif

             if ((spec_new.ge.(one-EPS_8_4)*species_base).and. &
                 (spec_new.le.species_max)) then
              ! do nothing
             else
              print *,"spec_new invalid: ",spec_new
              stop
             endif

             spec_new=max(spec_new,spec_old)

            else
             print *,"im_melt invalid: ",im_melt
             stop
            endif

            do im_inner=1,num_materials
             spec_comp=STATECOMP_STATES+(im_inner-1)*num_state_material+ &
               ENUM_SPECIESVAR+ispec
             snew(D_DECL(i,j,k),spec_comp)=spec_new
            enddo !do im_inner=1..num_materials

           else
            print *,"spec_old invalid(1): ",spec_old
            print *,"species_base: ",species_base
            print *,"species_max: ",species_max
            stop
           endif

          else
           print *,"species_base invalid: ",species_base
           stop
          endif

         else if (is_ice(im).eq.0) then
          ! do nothing
         else
          print *,"is_ice invalid"
          stop
         endif
        enddo ! im=1,num_materials

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
      end subroutine fort_apply_reaction

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
        ! vof,ref centroid,order,slope,intercept  x num_materials
      subroutine fort_convertmaterial( &
       tid, &
       im_outer, &     ! im_outer and im_opp_outer define an interface
       im_opp_outer, & ! between im_outer and im_opp_outer.
       level,finest_level, &
       nden, &
       nstate, &
       ntsat, &
       saturation_temp, &
       freezing_model, &
       Tanasawa_or_Schrage_or_Kassemi, &
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
       conductstate,DIMS(conductstate), & ! num_materials components
       nodevel,DIMS(nodevel), &
       JUMPFAB,DIMS(JUMPFAB), &
       TgammaFAB,DIMS(TgammaFAB), &
       LSold,DIMS(LSold), &
         ! in fort_ratemasschange
         ! LSnew=LSold - dt USTEFAN dot n, USTEFAN=|U|n  
         ! LSnew=LSold - dt * |U|
       LSnew,DIMS(LSnew), &
       recon,DIMS(recon), &
       snew,DIMS(snew), &
       EOS,DIMS(EOS), &
       swept,DIMS(swept) ) &
      bind(c,name='fort_convertmaterial')

      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use hydrateReactor_module
      use mass_transfer_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: im_outer
      integer, INTENT(in) :: im_opp_outer
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nden
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: ntsat
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      integer, INTENT(in) :: Tanasawa_or_Schrage_or_Kassemi(2*num_interfaces)
      integer, INTENT(in) :: mass_fraction_id(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      integer, INTENT(in) :: constant_density_all_time(num_materials)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in),target :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: min_stefan_velocity_for_dt
      integer, INTENT(in) :: vofbc(SDIM,2)
      real(amrex_real), INTENT(in),target :: xlo(SDIM)
      real(amrex_real), INTENT(in),target :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(inout) :: delta_mass(2*num_materials)
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(conductstate)
      integer, INTENT(in) :: DIMDEC(nodevel)
      integer, INTENT(in) :: DIMDEC(JUMPFAB)
      integer, INTENT(in) :: DIMDEC(TgammaFAB)
      integer, INTENT(in) :: DIMDEC(LSold)
      integer, INTENT(in) :: DIMDEC(LSnew)
      integer, INTENT(in) :: DIMDEC(recon)
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(EOS)
      integer, INTENT(in) :: DIMDEC(swept)

      real(amrex_real), target, INTENT(in) :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))

      real(amrex_real), target, INTENT(in) :: &
          conductstate(DIMV(conductstate),num_materials)
      real(amrex_real), pointer :: conductstate_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(in) :: nodevel(DIMV(nodevel),2*num_interfaces*SDIM)
      real(amrex_real), pointer :: nodevel_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(out) :: JUMPFAB(DIMV(JUMPFAB),2*num_interfaces)
      real(amrex_real), target, INTENT(out) :: TgammaFAB(DIMV(TgammaFAB),ntsat)
      real(amrex_real), pointer :: JUMPFAB_ptr(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: TgammaFAB_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: &
            LSold(DIMV(LSold),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LSold_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(out) :: &
            LSnew(DIMV(LSnew),num_materials)
      real(amrex_real), pointer :: LSnew_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(in) :: recon(DIMV(recon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: recon_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(out) :: snew(DIMV(snew),nstate)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(in) :: EOS(DIMV(EOS),nden)
      real(amrex_real), pointer :: EOS_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(out) :: swept(DIMV(swept),num_materials)
      real(amrex_real), pointer :: swept_ptr(D_DECL(:,:,:),:)

      real(amrex_real) :: denratio_factor

      integer i,j,k,dir
      integer i1,j1,k1
      integer im,im_opp
      integer ireverse
      integer iten
      integer im_source
      integer im_dest
      integer im_dest_crit
      integer im_source_crit
      integer im_primary
      integer iten_crit
      integer iten_outer
      integer ireverse_crit
      real(amrex_real) max_velnode
      real(amrex_real) velnode_test
      integer vcompsrc_snew,vcompdst_snew

      real(amrex_real) local_VOF(num_materials)

      real(amrex_real) oldvfrac(num_materials)
      real(amrex_real) newvfrac(num_materials)
      real(amrex_real) vof_super(num_materials)
      real(amrex_real) dF,dFdst,dFsrc
      real(amrex_real) den_dF(2)
      real(amrex_real) jump_strength

      real(amrex_real) volgrid
      real(amrex_real) cengrid(SDIM)
      real(amrex_real) new_centroid(num_materials,SDIM)
      real(amrex_real) old_centroid(num_materials,SDIM)
      real(amrex_real) EBVOFTOL
      real(amrex_real) SWEPTFACTOR
      real(amrex_real) SWEPTFACTOR_GFM
      integer SWEPTFACTOR_centroid
      real(amrex_real) LL
      real(amrex_real) Tgamma_default
      real(amrex_real) Ygamma_default
      integer Tsat_flag
      real(amrex_real) energy_source
      integer ccomp
      real(amrex_real) amount_used,methaneC_old,methaneC_new
      real(amrex_real) cvtotal,wttotal,Ftemp
      integer im_weight
      integer tcomp_wt
      integer vofcomp_raw
      integer vofcomp_recon
      integer vofcomp_raw_dest
      integer vofcomp_recon_source
      integer local_freezing_model
      integer mass_frac_id
      integer distribute_from_targ
      integer debugrate
      real(amrex_real) F_STEN(num_materials)
      real(amrex_real) F_STEN_CENTER(num_materials)
      real(amrex_real) unsplit_snew(num_materials*ngeom_raw)
      real(amrex_real) unsplit_lsnew(num_materials)
      real(amrex_real) oldLS_point(num_materials*(1+SDIM))
      real(amrex_real) dxmax,dxmaxLS
      integer recon_ncomp
      integer scomp
 
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdata_new(num_materials*ngeom_recon)
      real(amrex_real) multi_centroidA(num_materials,SDIM)
      real(amrex_real) volmat(num_materials)
      real(amrex_real) lsmat(num_materials)
      real(amrex_real) lsdata(num_materials)
      real(amrex_real) centroid_mat(SDIM,num_materials)
      real(amrex_real) multi_volume(num_materials)
      real(amrex_real) multi_area(num_materials)
      real(amrex_real) multi_cen(SDIM,num_materials)

      real(amrex_real) thermal_k_model_predict(2)  ! source,dest
      real(amrex_real) thermal_k_physical_base(2)  ! source,dest

      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_ofs(-nhalf:nhalf,SDIM)

      integer nmax,nhalf0
      integer klosten,khisten
      integer symmetry_flag,ntetbox
      real(amrex_real) u_xsten_grid(-nhalf:nhalf,SDIM)
      real(amrex_real) dxgrid(SDIM)
      real(amrex_real) u_xsten_updatecell(-nhalf:nhalf,SDIM)
      real(amrex_real) u_xsten_departmap(-1:1,SDIM)
      real(amrex_real) absolute_voltotal
      real(amrex_real) voltotal
      real(amrex_real) multi_volume_total
      integer igrid,jgrid,kgrid
      integer imac,jmac,kmac
      integer inode1,jnode1,knode1,u_imaterial,u_im,itri
      integer udir,udir2,inode,id,nlist
      integer matrix_status
      real(amrex_real) velnode
      real(amrex_real) xdepartnode(4*(SDIM-1),SDIM)
      real(amrex_real) xtargetnode(4*(SDIM-1),SDIM)
      real(amrex_real) tempdatanode(4*(SDIM-1))
      real(amrex_real) xinttri(SDIM+1,SDIM)
      real(amrex_real) xtri(SDIM+1,SDIM)
      real(amrex_real) xmaptri(SDIM+1,SDIM)
      real(amrex_real) xtargettri(SDIM+1,SDIM)
      real(amrex_real) datatri(SDIM+1)
      real(amrex_real) AA(SDIM+1,SDIM+1)
      real(amrex_real) AS(SDIM,SDIM)
      real(amrex_real) ASINV(SDIM,SDIM)
      real(amrex_real) bb(SDIM+1)
      real(amrex_real) AAINV(SDIM+1,SDIM+1)
      real(amrex_real) bbINV(SDIM+1)
      real(amrex_real) xx(SDIM+1)
      real(amrex_real) coeff(SDIM,SDIM+1)
      real(amrex_real) coeffINV(SDIM,SDIM+1)
      integer n,ivert
      real(amrex_real) nn(SDIM)
      real(amrex_real) nnmap(SDIM)
      real(amrex_real) uncaptured_volume
      real(amrex_real) uncaptured_centroid(SDIM)
      integer shapeflag
      integer tessellate
      real(amrex_real) volcell
      real(amrex_real) volcell_ofs
      real(amrex_real) tempvfrac
      real(amrex_real) cencell(SDIM)
      real(amrex_real) cencell_ofs(SDIM)
      real(amrex_real) tempcen(SDIM)
      integer do_unsplit_advection
      integer interface_near(2*num_interfaces)
      integer local_mask
      integer im_primary_new
      integer im_primary_old
      integer im_primary_local
      integer away_from_interface
      real(amrex_real) solid_vof_new,solid_vof_old
      integer mtype
      real(amrex_real) vapor_den
      real(amrex_real) condensed_den
      real(amrex_real) local_cv_or_cp
      integer speccomp_mod
      integer default_comp
      integer ncomp_per_tsat
      integer iprobe
      integer iprobe_vapor
      integer iprobe_condensed
      integer im_probe
      integer im_vapor
      integer im_condensed
      integer dencomp_probe
      integer tcomp_probe
      integer mfrac_comp_probe
      integer ispec_probe
      real(amrex_real) temp_mix_new(2)
      real(amrex_real) mass_frac_new(2)
      real(amrex_real) density_old(2)
      real(amrex_real) temperature_old(2)
      real(amrex_real) species_old(2)
      real(amrex_real) delta_mass_local(2) ! iprobe==1: source   iprobe==2: dest
      real(amrex_real) :: xPOINT_supermesh(SDIM)
      real(amrex_real) :: xPOINT_GFM(SDIM)
      integer im_old_crit
      integer im_new_crit
      real(amrex_real) DATA_FLOOR
      integer combine_flag
      integer nsolve_interp
      real(amrex_real) xtarget_interp(SDIM)
      real(amrex_real) old_xI(SDIM)
      real(amrex_real) old_nrm(SDIM)
      integer interp_to_new_supermesh

      integer im_primary_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) XC_sten(D_DECL(-1:1,-1:1,-1:1),SDIM)
      real(amrex_real) VF_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) LS_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) temperature_sten(D_DECL(-1:1,-1:1,-1:1))
      real(amrex_real) massfrac_sten(D_DECL(-1:1,-1:1,-1:1))

      integer, parameter :: continuous_mof_parm=STANDARD_MOF
      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))

      integer :: grid_index(SDIM)
      integer, parameter :: grid_level=-1

      integer, parameter :: use_ls_data=0
      integer, parameter :: mof_verbose=0
      real(amrex_real) LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)
      integer order_probe(2)
      real(amrex_real) nslope_probe(SDIM,2)
      real(amrex_real) intercept_probe(2)
      real(amrex_real) nslope_dest(SDIM)
      real(amrex_real) intercept_dest
      real(amrex_real) LS_dest_old,LS_dest_new
      real(amrex_real) mass_frac_limit
      integer vofcomp_local
      integer im_local
      integer im_trust
      integer im_distrust
      real(amrex_real) fixed_vfrac_sum
      real(amrex_real) fixed_centroid_sum(SDIM)
      real(amrex_real) avail_vfrac

      JUMPFAB_ptr=>JUMPFAB
      TgammaFAB_ptr=>TgammaFAB
      LSnew_ptr=>LSnew
      snew_ptr=>snew
      swept_ptr=>swept
      maskcov_ptr=>maskcov
      LSold_ptr=>LSold
      recon_ptr=>recon
      EOS_ptr=>EOS
      nodevel_ptr=>nodevel
      conductstate_ptr=>conductstate

      if ((tid.lt.0).or. &
          (tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: fort_convertmaterial
      if ((nmax.lt.100).or.(nmax.gt.2000)) then
       print *,"nmax invalid"
       stop
      endif

      nhalf0=1

      recon_ncomp=num_materials*ngeom_recon

      ncomp_per_tsat=EXTRAP_PER_TSAT

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

      if (ngrow_distance.ne.4) then
       print *,"ngrow_distance invalid"
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
           (im_outer.le.num_materials)) then
        ! do nothing
       else
        print *,"im_outer invalid"
        stop
       endif
       if ((im_opp_outer.ge.1).and. &
           (im_opp_outer.le.num_materials)) then
        ! do nothing
       else
        print *,"im_opp_outer invalid"
        stop
       endif
      else
       print *,"cannot have im_outer==im_opp_outer"
       stop
      endif
      call get_iten(im_outer,im_opp_outer,iten_outer)

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)

      EBVOFTOL=VOFTOL
      if (dt.lt.one) then
       EBVOFTOL=EBVOFTOL*dt
      endif

       ! For Tsatfab: 1st num_interfaces 
       ! components are the status, then next 2*num_interfaces
       ! components go:
       ! T_gamma_1,Y_gamma_1,
       ! T_gamma_2,Y_gamma_2, ....
      if (ntsat.eq.EXTRAP_NCOMP_TSAT) then
       ! do nothing
      else
       print *,"ntsat invalid"
       stop
      endif
      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid in fort_convertmaterial"
       print *,"nden=",nden
       print *,"num_materials=",num_materials
       print *,"num_state_material=",num_state_material
       stop
      endif
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      do im=1,num_materials-1
       do im_opp=im+1,num_materials
        do ireverse=0,1
         call get_iten(im,im_opp,iten)
         local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
         distribute_from_targ= &
            distribute_from_target(iten+ireverse*num_interfaces)
         LL=get_user_latent_heat(iten+ireverse*num_interfaces,room_temperature,1)
         mass_frac_id=0

         if (is_multi_component_evapF(local_freezing_model, &
              Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces), &
              LL).eq.0) then 
          ! do nothing
         else if (is_multi_component_evapF(local_freezing_model, &
               Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces), &
               LL).eq.1) then 

          mass_frac_id=mass_fraction_id(iten+ireverse*num_interfaces)
          if ((mass_frac_id.ge.1).and. &
              (mass_frac_id.le.num_species_var)) then
           ! do nothing
          else
           print *,"mass_frac_id invalid"
           stop
          endif
         else
          print *,"is_multi_component_evapF invalid"
          stop
         endif

         if (is_valid_freezing_modelF(local_freezing_model).eq.1) then
          ! do nothing
         else
          print *,"local_freezing_model invalid in fort_convertmaterial"
          print *,"local_freezing_model= ",local_freezing_model
          print *,"iten,ireverse,num_interfaces ",iten,ireverse,num_interfaces
          stop
         endif
        enddo ! ireverse
       enddo ! im_opp
      enddo !im

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)

      call checkbound_array(fablo,fabhi,JUMPFAB_ptr,ngrow_distance,-1)
      call checkbound_array(fablo,fabhi,TgammaFAB_ptr,ngrow_distance,-1)
      call checkbound_array(fablo,fabhi,LSold_ptr,ngrow_distance,-1)
      call checkbound_array(fablo,fabhi,LSnew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,recon_ptr,1,-1)
      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,EOS_ptr,1,-1)
      call checkbound_array(fablo,fabhi,swept_ptr,0,-1)
      call checkbound_array(fablo,fabhi,nodevel_ptr,1,-1)
      call checkbound_array(fablo,fabhi,conductstate_ptr,1,-1)

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

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       grid_index(1)=i
       grid_index(2)=j
       if (SDIM.eq.3) then
        grid_index(SDIM)=k
       endif

       if (local_mask.eq.1) then

        do ireverse=0,1 
         if (JUMPFAB(D_DECL(i,j,k), &
             iten_outer+ireverse*num_interfaces).eq.zero) then
          ! do nothing
         else
          print *, &
           "JUMPFAB(D_DECL(i,j,k),iten_outer+ireverse*num_interfaces) bad"
          stop
         endif
        enddo

        call gridsten_level(xsten,i,j,k,level,nhalf)
        call gridsten_level(u_xsten_updatecell,i,j,k,level,nhalf)

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volgrid,cengrid,SDIM)

        do_unsplit_advection=0
        do im=1,2*num_interfaces
         interface_near(im)=0
        enddo

        do im=1,num_materials
         lsmat(im)=LSold(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(lsmat,im_primary)

        if (is_rigid(im_primary).eq.0) then

         do im=1,num_materials
          vofcomp_recon=(im-1)*ngeom_recon+1
          F_STEN(im)=zero
          F_STEN_CENTER(im)=recon(D_DECL(i,j,k),vofcomp_recon)

          do k1=klosten,khisten
          do j1=-1,1
          do i1=-1,1
           F_STEN(im)=F_STEN(im)+recon(D_DECL(i+i1,j+j1,k+k1),vofcomp_recon)
          enddo
          enddo
          enddo ! i1,j1,k1
         enddo ! im=1..num_materials

         do im=1,num_materials-1
          do im_opp=im+1,num_materials

           call get_iten(im,im_opp,iten)

           do ireverse=0,1
            if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
             print *,"im or im_opp bust 10"
             stop
            endif

            LL=get_user_latent_heat(iten+ireverse*num_interfaces,room_temperature,1)

            if ((is_rigid(im).eq.1).or. &
                (is_rigid(im_opp).eq.1)) then
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
                interface_near(iten+ireverse*num_interfaces)=1
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

           do im=1,num_materials-1
            do im_opp=im+1,num_materials
             call get_iten(im,im_opp,iten)
             do ireverse=0,1
               !interface_near==0 if im<>im_outer or im_opp<>im_opp_outer
              if (interface_near(iten+ireverse*num_interfaces).eq.1) then
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
               if ((iten.ge.1).and.(iten.le.num_interfaces)) then
                scomp=(iten+num_interfaces*ireverse-1)*SDIM
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
               else if ((im_dest_crit.ge.1).and. &
                        (im_dest_crit.le.num_materials)) then
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
              else if (interface_near(iten+ireverse*num_interfaces).eq.0) then
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
          if ((im_dest_crit.lt.1).or.(im_dest_crit.gt.num_materials)) then
           print *,"im_dest_crit invalid"
           stop
          endif
          if ((im_source_crit.lt.1).or.(im_source_crit.gt.num_materials)) then
           print *,"im_source_crit invalid"
           stop
          endif
          if ((iten_crit.lt.1).or.(iten_crit.gt.num_interfaces)) then
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

          if (max_velnode.gt.zero) then

           symmetry_flag=0 
           call get_ntetbox(ntetbox,symmetry_flag,SDIM)

           absolute_voltotal=zero
           voltotal=zero
           do u_im=1,num_materials
            volmat(u_im)=zero
            do udir=1,SDIM
             centroid_mat(udir,u_im)=zero
            enddo
            lsmat(u_im)=zero
           enddo

            ! backwards tracing of characteristics:
            ! add up contributions from neighboring cells.
            ! (i,j,k) is the cell to be updated.
           do kgrid=klosten,khisten
           do jgrid=-1,1
           do igrid=-1,1

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

             scomp=(iten_crit+ireverse_crit*num_interfaces-1)*SDIM

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

              do u_im=1,num_materials

               lsdata(u_im)= &
                LSold(D_DECL(i+igrid,j+jgrid,k+kgrid),u_im) 

               vofcomp_recon=(u_im-1)*ngeom_recon+1
               mofdata(vofcomp_recon)= &
                recon(D_DECL(i+igrid,j+jgrid,k+kgrid),vofcomp_recon) 
               do udir=1,SDIM 
                mofdata(vofcomp_recon+udir)= &
                 recon(D_DECL(i+igrid,j+jgrid,k+kgrid),vofcomp_recon+udir)
               enddo
               mofdata(vofcomp_recon+SDIM+1)= &
                recon(D_DECL(i+igrid,j+jgrid,k+kgrid),vofcomp_recon+SDIM+1)!ord
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
              enddo ! u_im=1..num_materials

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
                 EPS_8_4, &
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
                 SDIM,shapeflag)

               multi_volume_total=zero
               do u_im=1,num_materials
                if (is_rigid(u_im).eq.0) then 
                 multi_volume_total=multi_volume_total+multi_volume(u_im)
                else if (is_rigid(u_im).eq.1) then 
                 ! do nothing
                else
                 print *,"is_rigid invalid MASS_TRANSFER_3D.F90"
                 stop
                endif

                volmat(u_im)=volmat(u_im)+multi_volume(u_im)
                do udir=1,SDIM
                 centroid_mat(udir,u_im)=centroid_mat(udir,u_im)+ &
                   multi_volume(u_im)*multi_cen(udir,u_im)
                enddo

               enddo ! u_im=1..num_materials

               voltotal=voltotal+multi_volume_total

               do u_im=1,num_materials
                lsmat(u_im)=lsmat(u_im)+multi_volume_total*lsdata(u_im)
               enddo
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

           do u_imaterial=1,num_materials
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

           if (interface_near(iten+ireverse*num_interfaces).ne.1) then
            print *,"interface_near invalid"
            stop
           endif

           Tgamma_default=saturation_temp(iten+ireverse*num_interfaces)
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

           LL=get_user_latent_heat(iten+ireverse*num_interfaces,room_temperature,1)
           local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
           distribute_from_targ= &
             distribute_from_target(iten+ireverse*num_interfaces)
           mass_frac_id=mass_fraction_id(iten+ireverse*num_interfaces)
           
            ! TgammaFAB has valid values corresponding to the
            ! given value for "ireverse". 
           if ((Tsat_flag.eq.1).or.(Tsat_flag.eq.2)) then
             Tgamma_default=TgammaFAB(D_DECL(i,j,k), &
              num_interfaces+(iten-1)*ncomp_per_tsat+1)
             if ((mass_frac_id.ge.1).and. &
                 (mass_frac_id.le.num_species_var)) then
              Ygamma_default=TgammaFAB(D_DECL(i,j,k), &
               num_interfaces+(iten-1)*ncomp_per_tsat+2)
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

           do u_imaterial=1,num_materials

            vofcomp_raw=(u_imaterial-1)*ngeom_raw+1

            tempvfrac=volmat(u_imaterial)/voltotal
            if ((tempvfrac.ge.EBVOFTOL).and.(tempvfrac.le.one+EPS1)) then
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
              tempcen(udir)= &
               centroid_mat(udir,u_imaterial)/volmat(u_imaterial)-cencell(udir)
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

           enddo  ! u_imaterial=1..num_materials

           vcompsrc_snew=STATECOMP_MOF+(im_source-1)*ngeom_raw+1
           vcompdst_snew=STATECOMP_MOF+(im_dest-1)*ngeom_raw+1

           do u_imaterial=1,num_materials
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
           enddo ! u_imaterial=1,num_materials

           do u_imaterial=1,num_materials*(1+SDIM)
            oldLS_point(u_imaterial)=LSold(D_DECL(i,j,k),u_imaterial)
           enddo
           call normalize_LS_normals(oldLS_point)

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

            ! declared in MOF.F90: checks both is_rigid and non is_rigid
            ! materials.
           call get_primary_material(unsplit_lsnew,im_primary_new)
           call get_primary_material(oldLS_point,im_primary_old)
           call combine_solid_VOF(newvfrac,solid_vof_new,im_primary_local)
           call combine_solid_VOF(oldvfrac,solid_vof_old,im_primary_local)

           away_from_interface=0
           
           if ((is_rigid(im_primary_new).eq.1).or. &
               (is_rigid(im_primary_old).eq.1)) then
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

            if ((im_primary_old.ne.im_source).and. &
                (im_primary_old.ne.im_dest)) then

             if (im_primary_old.ne.im_primary_new) then
              ! do not touch LSnew if im_primary_old<>
              ! im_source,im_dest,im_primary_new
             else if (im_primary_old.eq.im_primary_new) then
              LSnew(D_DECL(i,j,k),im_source)=unsplit_lsnew(im_source)
              LSnew(D_DECL(i,j,k),im_dest)=unsplit_lsnew(im_dest)
             else
              print *,"im_primary_new or im_primary_old invalid"
              stop
             endif

            else if ((im_primary_old.eq.im_source).or. &
                     (im_primary_old.eq.im_dest)) then

             if ((im_primary_new.ne.im_source).and. &
                 (im_primary_new.ne.im_dest)) then
              ! do nothing
             else if ((im_primary_new.eq.im_source).or. &
                      (im_primary_new.eq.im_dest)) then
              if (unsplit_lsnew(im_dest).lt.oldLS_point(im_dest)) then
               ! do nothing
              else if (unsplit_lsnew(im_source).gt.oldLS_point(im_source)) then
               ! do nothing
              else if ((unsplit_lsnew(im_dest).ge. &
                        oldLS_point(im_dest)).and. &
                       (unsplit_lsnew(im_source).le. &
                        oldLS_point(im_source))) then
               LSnew(D_DECL(i,j,k),im_source)=unsplit_lsnew(im_source)
               LSnew(D_DECL(i,j,k),im_dest)=unsplit_lsnew(im_dest)
              else
               print *,"unsplit_lsnew or lsdata invalid"
               stop
              endif
             else
              print *,"im_primary_new invalid"
              stop
             endif
            else
             print *,"im_primary_old invalid"
             stop
            endif

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
           do im_local=1,num_materials
            if ((im_local.eq.im_source).or. &
                (im_local.eq.im_dest)) then
             ! do nothing
            else if ((im_local.ge.1).and.(im_local.le.num_materials)) then
             newvfrac(im_local)=oldvfrac(im_local)
             do udir=1,SDIM 
              new_centroid(im_local,udir)=old_centroid(im_local,udir)
             enddo
             if (is_rigid(im_local).eq.0) then
              fixed_vfrac_sum=fixed_vfrac_sum+newvfrac(im_local)
              do udir=1,SDIM 
               fixed_centroid_sum(udir)=fixed_centroid_sum(udir)+ &
                 newvfrac(im_local)*new_centroid(im_local,udir)
              enddo
             else if (is_rigid(im_local).eq.1) then
              ! ignore, solids are embedded
             else
              print *,"is_rigid invalid MASS_TRANSFER_3D.F90"
              stop
             endif
            else
             print *,"im_local became corrupt"
             stop
            endif
           enddo ! im_local=1..num_materials

           if ((fixed_vfrac_sum.ge.one-VOFTOL).and. &
               (fixed_vfrac_sum.le.one+EPS1)) then
            avail_vfrac=zero
           else if ((fixed_vfrac_sum.ge.-EPS1).and. &
                    (fixed_vfrac_sum.le.zero)) then
            avail_vfrac=one
           else if ((fixed_vfrac_sum.gt.zero).and. &
                    (fixed_vfrac_sum.le.one-VOFTOL)) then
            avail_vfrac=one-fixed_vfrac_sum
           else
            print *,"fixed_vfrac_sum invalid:",fixed_vfrac_sum
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
                else if ((newvfrac(im_distrust).ge.-EPS1).and. &
                         (newvfrac(im_distrust).le.VOFTOL)) then
                 new_centroid(im_distrust,udir)=cengrid(udir)
                else
                 print *,"newvfrac(im_distrust) invalid"
                 stop
                endif
               enddo ! udir=1..sdim
              else if ((fixed_vfrac_sum.ge.-EPS1).and. &
                       (fixed_vfrac_sum.le.VOFTOL)) then
               ! do nothing
              else
               print *,"fixed_vfrac_sum invalid: ",fixed_vfrac_sum
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

             dencomp_probe=(im_probe-1)*num_state_material+1+ENUM_DENVAR
             tcomp_probe=dencomp_probe+1
             if ((mass_frac_id.ge.1).and. &
                 (mass_frac_id.le.num_species_var)) then
              mfrac_comp_probe=tcomp_probe+mass_frac_id
              ispec_probe=(mass_frac_id-1)*num_materials+im_probe
             else if (mass_frac_id.eq.0) then
              mfrac_comp_probe=0
              ispec_probe=0
             else
              print *,"mass_frac_id invalid"
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
               print *,"density_old(iprobe) invalid: ",density_old(iprobe)
               stop
              endif

             else if ((mtype.ge.1).and.(mtype.le.MAX_NUM_EOS)) then

              density_old(iprobe)=fort_denconst(im_probe)

             else
              print *,"mtype invalid: ",mtype
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
                      (oldvfrac(im_probe).ge.-EPS1)) then
              temperature_old(iprobe)=Tgamma_default
              species_old(iprobe)=Ygamma_default
             else
              print *,"oldvfrac invalid: ",oldvfrac(im_probe)
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
            vapor_den=density_old(iprobe_vapor)
            condensed_den=density_old(iprobe_condensed)
            if ((vapor_den.gt.zero).and.(condensed_den.gt.zero)) then
             ! do nothing
            else
             print *,"vapor_den or condensed_den invalid"
             stop
            endif
            if (is_multi_component_evapF(local_freezing_model, &
              Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces), &
              LL).eq.0) then 
             ! do nothing
            else if (is_multi_component_evapF(local_freezing_model, &
              Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces), &
              LL).eq.1) then 
             if ((mass_frac_id.ge.1).and. &
                 (mass_frac_id.le.num_species_var)) then
              ! do nothing
             else
              print *,"mass_frac_id invalid"
              stop
             endif

            else
             print *,"is_multi_component_evapF invalid 1"
             stop
            endif

            mass_frac_new(1)=Ygamma_default ! source
            mass_frac_new(2)=Ygamma_default ! dest
            temp_mix_new(1)=Tgamma_default
            temp_mix_new(2)=Tgamma_default

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
             !  dF/(tnp1-tn)  * (tswept-tnp1) = 1/2-Fnp1
             !  SWEPTFACTOR=(tnp1-tswept)/(tnp1-tn)  = 
             !              (Fnp1-1/2)/(Fnp1-Fn)
             !  new approach (motivated by supermesh):
             !  1. given an estimate of phase change velocity "USTEFAN",
             !     dt is chosen such that USTEFAN * dt <= dx/4
             !     This assures that a swept cell will NOT be full at tnp1.
             !  2. Strategy is this: 
             !     a) we have the signed distance function at tn.
             !     b) let x^* be either (i) cell center if GFM, or (ii)
             !        cell centroid at tnp1 of the destination material.
             !     c) if phi(tn,x^*)<0 then:
             !      find MOF reconstruction and determine if
             !      phi^reconstruct(tnp1,x^*)>0
             !      if yes, then:
             !       phi(t)=(phi_np1-phi_n)/(tnp1-tn) * (t-tnp1)+phi_np1
             !       0=(phi_np1-phi_n)/(tnp1-tn)  * (tswept-tnp1)+phi_np1 
             !       SWEPTFACTOR=(tnp1-tswept)/(tnp1-tn)=
             !          phi_np1/(phi_np1-phi_n)

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
            do u_im=1,num_materials
             vofcomp_recon=(u_im-1)*ngeom_recon+1

             mofdata(vofcomp_recon)=oldvfrac(u_im)
             mofdata_new(vofcomp_recon)=newvfrac(u_im)
             vof_super(u_im)=newvfrac(u_im)

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
            enddo ! u_im=1..num_materials

             ! LS=n dot (x-x0)+intercept
            call multimaterial_MOF( &
             bfact,dx,u_xsten_updatecell,nhalf, &
             mof_verbose, & !=0
             use_ls_data, & !=0
             LS_stencil, &
             geom_xtetlist(1,1,1,tid+1), &
             geom_xtetlist_old(1,1,1,tid+1), &
             nmax, &
             nmax, &
             mofdata_new, & !intent(inout)
             vof_super, &
             multi_centroidA, &
             continuous_mof_parm, & !=STANDARD_MOF
             cmofsten, &
             grid_index, &
             grid_level, &
             SDIM)

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
             ! If the new centroid is NOT in the old dest material, then
             ! we initialize the temperature to be the interface temperature.
            SWEPTFACTOR_centroid=0
            if ((newvfrac(im_dest).gt.zero).and. &
                (newvfrac(im_dest).le.one+EPS1)) then

             tessellate=3
             call multi_get_volumePOINT( &
               tessellate, &
               bfact,dx, &
               u_xsten_updatecell,nhalf, &  ! absolute coordinate system
               mofdata, &
               xPOINT_supermesh, & !new centroid: absolute coordinate system
               im_old_crit, &
               SDIM)

              ! im_old_crit=material at t=tn that occupies the new 
              ! centroid location
             if ((im_old_crit.eq.im_dest).or. &
                 (oldvfrac(im_dest).ge.half)) then
              ! do nothing
             else if ((im_old_crit.ge.1).and. &
                      (im_old_crit.le.num_materials).and. &
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
               do u_im=1,num_materials
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

            !dF.gt.EBVOFTOL
            !oldLS_point(u_imaterial)=LSold(D_DECL(i,j,k),u_imaterial)
            if (oldLS_point(im_dest).ge.zero) then
             do udir=1,SDIM
              old_nrm(udir)=oldLS_point(num_materials+(im_source-1)*SDIM+udir)
              old_xI(udir)=u_xsten_updatecell(0,udir)- &
                oldLS_point(im_source)*old_nrm(udir)
             enddo
            else if (oldLS_point(im_source).ge.zero) then
             do udir=1,SDIM
              old_nrm(udir)=oldLS_point(num_materials+(im_dest-1)*SDIM+udir)
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

             dencomp_probe=(im_probe-1)*num_state_material+ENUM_DENVAR+1
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

               do k1=klosten,khisten
               do j1=-1,1
               do i1=-1,1
                call gridsten_level(xsten_ofs,i+i1,j+j1,k+k1,level,nhalf)
                call Box_volumeFAST(bfact,dx,xsten_ofs,nhalf, &
                 volcell_ofs,cencell_ofs,SDIM)

                do im_local=1,num_materials
                 vofcomp_local=(im_local-1)*ngeom_recon+1
                 local_VOF(im_local)= &
                   recon(D_DECL(i+i1,j+j1,k+k1),vofcomp_local)
                enddo !im_local=1..num_materials

                call get_primary_material_VFRAC(local_VOF, &
                  im_primary_sten(D_DECL(i1,j1,k1)))

                do udir=1,SDIM
                 XC_sten(D_DECL(i1,j1,k1),udir)= &
                  recon(D_DECL(i+i1,j+j1,k+k1),vofcomp_recon+udir)+ &
                  cencell_ofs(udir)
                enddo
                VF_sten(D_DECL(i1,j1,k1))=local_VOF(im_probe)

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
                   (oldvfrac(im_probe).ge.-EPS1)) then
                temp_mix_new(iprobe)=Tgamma_default
                mass_frac_new(iprobe)=Ygamma_default
               else if ((oldvfrac(im_probe).ge.VOFTOL).and. &
                        (oldvfrac(im_probe).le.one+EPS1)) then

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
                 interp_to_new_supermesh, &!interp_to_new_supermesh(tsatflag)=1
                 bfact, &
                 level, &
                 finest_level, &
                 dx,xlo, &
                 u_xsten_updatecell,nhalf, &
                 temperature_sten, &
                 XC_sten, &
                 old_xI, &
                 xtarget_interp, &
                 im_probe, & 
                 im_primary_sten, & 
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
                 im_probe, & 
                 im_primary_sten, & 
                 VF_sten, &
                 LS_sten, &
                 Ygamma_default, &
                 mass_frac_new(iprobe))

               else
                print *,"oldvfrac(im_probe) invalid: ",oldvfrac(im_probe)
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

             snew(D_DECL(i,j,k),STATECOMP_STATES+tcomp_probe)= &
                     temp_mix_new(iprobe)
             if (mfrac_comp_probe.eq.0) then
              ! do nothing
             else if (mfrac_comp_probe.gt.0) then
              snew(D_DECL(i,j,k),STATECOMP_STATES+mfrac_comp_probe)=  &
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

            if (is_valid_freezing_modelF(local_freezing_model).eq.1) then
             ! do nothing
            else
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

             dencomp_probe=(im_probe-1)*num_state_material+ENUM_DENVAR+1
             tcomp_probe=dencomp_probe+1

              ! iprobe==1: source
              ! iprobe==2: dest
              ! (rho F)^new - (rho F)^old
             den_dF(iprobe)=delta_mass_local(iprobe)

             if (newvfrac(im_probe).gt.one+EPS1) then
              print *,"newvfrac(im_probe) overflow"
              stop
             else if (newvfrac(im_probe).gt.one) then
              newvfrac(im_probe)=one
             else if (newvfrac(im_probe).lt.-EPS1) then
              print *,"newvfrac(im_probe) underflow"
              stop
             else if (newvfrac(im_probe).lt.zero) then
              newvfrac(im_probe)=zero
             else if ((newvfrac(im_probe).ge.zero).and. &
                      (newvfrac(im_probe).le.one)) then
              ! do nothing
             else
              print *,"newvfrac(im_probe) invalid: ",newvfrac(im_probe)
              stop
             endif

             thermal_k_model_predict(iprobe)= &
               conductstate(D_DECL(i,j,k),im_probe)* &
               fort_heatflux_factor(im_probe)
             thermal_k_physical_base(iprobe)= &
               get_user_heatviscconst(im_probe)

             if (fort_heatflux_factor(im_probe).ge.zero) then
              ! do nothing
             else
              print *,"fort_heatflux_factor(im_probe) invalid"
              stop
             endif

             if (thermal_k_model_predict(iprobe).ge.zero) then
              ! do nothing
             else
              print *,"thermal_k_model_predict(iprobe) invalid"
              stop
             endif
             if (thermal_k_physical_base(iprobe).ge.zero) then
              ! do nothing
             else
              print *,"thermal_k_physical_base(iprobe) invalid"
              stop
             endif

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

            if (abs(denratio_factor).le.EPS_8_4) then
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

            JUMPFAB(D_DECL(i,j,k),iten+ireverse*num_interfaces)=jump_strength
          
            if (dF.le.EBVOFTOL) then
             print *,"expecting dF>EBVOFTOL"
             stop
            else if (dF.ge.zero) then

              ! find mass weighted average of cv

             cvtotal=zero
             wttotal=zero
             do im_weight=1,num_materials
              vofcomp_recon=(im_weight-1)*ngeom_recon+1
              Ftemp=recon(D_DECL(i,j,k),vofcomp_recon)*fort_denconst(im_weight)
              local_cv_or_cp=get_user_stiffCP(im_weight)
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
             if (is_GFM_freezing_modelF(local_freezing_model).eq.1) then

               !F dt=Fn (tnp1-t) + Fnp1 (t-tn)
               !dt/2 = Fn tnp1 - Fnp1 tn + tswept(Fnp1-Fn)
               !tswept=(dt/2 - Fn tnp1 + Fnp1 tn)/(Fnp1-Fn)=
               !  (dt/2 - Fn (tn+dt) + Fnp1 tn)/dF=
               !  (dt/2 + dF tn - Fn dt)/dF
               !tswept-tn=dt(1/2 - Fn)/dF
               !tnp1-tswept=-dt(1/2-Fn)/dF+dt=
               !  dt(dF-1/2+Fn)/dF=dt(Fnp1-1/2)/dF=
               !  (tnp1-tn)*(Fnp1-1/2)/(Fnp1-Fn)
               !sweptfactor=(tnp1-tswept)/(tnp1-tn)=
               !  (Fnp1-1/2)/(Fnp1-Fn)

              do udir=1,SDIM
               xPOINT_GFM(udir)=u_xsten_updatecell(0,udir)
              enddo

              LS_dest_old=LSold(D_DECL(i,j,k),im_dest)

              tessellate=3
              call multi_get_volumePOINT( &
               tessellate, &
               bfact,dx, &
               u_xsten_updatecell,nhalf, &  ! absolute coordinate system
               mofdata, &
               xPOINT_GFM, & ! absolute coordinate system
               im_old_crit,SDIM)

              tessellate=3
              call multi_get_volumePOINT( &
                tessellate, &
                bfact,dx, &
                u_xsten_updatecell,nhalf, &  ! absolute coordinate system
                mofdata_new, &
                xPOINT_GFM, & ! absolute coordinate system
                im_new_crit,SDIM)

              if ((newvfrac(im_dest).gt.zero).and. &
                  (newvfrac(im_dest).le.one+EPS1)) then

                ! im_new_crit "owns" xPOINT_GFM
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

                if (order_probe(1).eq.0) then ! source order = 0
                 do udir=1,SDIM 
                  nslope_dest(udir)=nslope_probe(udir,2)
                 enddo
                 intercept_dest=intercept_probe(2)
                else if (order_probe(2).eq.0) then  ! dest order = 0
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
                     (xPOINT_GFM(udir)-u_xsten_updatecell(0,udir))
                enddo

               else if ((im_new_crit.ge.1).and. &
                        (im_new_crit.le.num_materials).and. &
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

              if ((dF.gt.zero).and. &
                  (dF.le.one+EPS1)) then
               ! do nothing
              else
               print *,"dF invalid: ",dF
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
                SWEPTFACTOR=EPS1
               else if ((order_probe(1).eq.0).and. &
                        (order_probe(2).eq.0)) then 
                SWEPTFACTOR=one ! default
               else if ((order_probe(1).gt.0).or. &
                        (order_probe(2).gt.0)) then

                if ((LS_dest_old.eq.zero).and.(LS_dest_new.gt.zero)) then
                 SWEPTFACTOR=one
                else if ((LS_dest_old.lt.zero).and.(LS_dest_new.eq.zero)) then
                 SWEPTFACTOR=EPS1
                else if ((LS_dest_old.ge.zero).or. &
                         (LS_dest_new.le.zero)) then
                 SWEPTFACTOR=one
                else if ((LS_dest_new-LS_dest_old.gt.zero).and. &
                         (LS_dest_new.ge.zero).and. &
                         (LS_dest_old.le.zero)) then
                 SWEPTFACTOR=LS_dest_new/ &
                        (LS_dest_new-LS_dest_old)
                else
                 print *,"LS_dest_new or LS_dest_old invalid"
                 print *,"LS_dest_new ",LS_dest_new
                 print *,"LS_dest_old ",LS_dest_old
                 stop
                endif
                if (SWEPTFACTOR.le.EPS1) then
                 SWEPTFACTOR=EPS1
                endif

               else
                print *,"order_probe bust"
                stop
               endif

               if ((SWEPTFACTOR.ge.EPS2).and. &
                   (SWEPTFACTOR.le.one)) then
                ! do nothing
               else
                print *,"SWEPTFACTOR invalid: ",SWEPTFACTOR
                print *,"dF=",dF
                print *,"im_dest=",im_dest
                print *,"newvfrac(im_dest) ",newvfrac(im_dest)
                print *,"oldvfrac(im_dest) ",oldvfrac(im_dest)
                print *,"EPS2 ",EPS2
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
               do im_weight=1,num_materials
                tcomp_wt=STATECOMP_STATES+(im_weight-1)*num_state_material+ &
                        ENUM_TEMPERATUREVAR+1
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

             else if &
               (is_hydrate_freezing_modelF(local_freezing_model).eq.1) then 

              if (distribute_from_targ.ne.0) then
               print *,"distribute_from_targ invalid"
               stop
              endif
              if (dF.gt.zero) then
               if (num_species_var.ne.1) then
                print *,"num_species_var invalid"
                stop
               endif
   
               ! Hydrate_energy_soure_term is declared in HYDRATE_REACTOR.F90 
               call Hydrate_energy_source_term( &
                dF,dt, &
                thermal_k_physical_base(1), & 
                thermal_k_model_predict(1), & 
                energy_source,LL)
               call Methane_usage(dF,dt, &
                fort_speciesviscconst(im_dest),amount_used)

               ccomp=STATECOMP_STATES+(im_dest-1)*num_state_material+ &
                       ENUM_SPECIESVAR+1
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
      
               do im_weight=1,num_materials
                tcomp_wt=STATECOMP_STATES+(im_weight-1)*num_state_material+ &
                    ENUM_TEMPERATUREVAR+1
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
              print *,"FIX ME: cavitation case not considered here"
              stop
             else
              print *,"local_freezing_model invalid in fort_convertmaterial(2)"
              print *,"local_freezing_model= ",local_freezing_model
              print *,"iten,ireverse,num_interfaces ", &
                 iten,ireverse,num_interfaces
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
              snew(D_DECL(i,j,k),STATECOMP_MOF+(im_probe-1)*ngeom_raw+1)= &
                      newvfrac(im_probe)

             enddo ! iprobe=1,2

             do u_im=1,num_materials*ngeom_recon
              mofdata(u_im)=zero
             enddo

             do u_im=1,num_materials
              vofcomp_recon=(u_im-1)*ngeom_recon+1
              vofcomp_raw=(u_im-1)*ngeom_raw+1
              do dir=0,SDIM
               mofdata(vofcomp_recon+dir)= &
                 snew(D_DECL(i,j,k),STATECOMP_MOF+vofcomp_raw+dir)
              enddo
             enddo ! u_im=1..num_materials

             ! sum of F_fluid=1
             ! sum of F_rigid<=1
             tessellate=0
             call make_vfrac_sum_ok_base( &
               cmofsten, &
               u_xsten_updatecell,nhalf, &
               continuous_mof_parm, &
               bfact,dx, &
               tessellate,mofdata,SDIM)

             do u_im=1,num_materials
              vofcomp_recon=(u_im-1)*ngeom_recon+1
              vofcomp_raw=(u_im-1)*ngeom_raw+1
              do dir=0,SDIM
               snew(D_DECL(i,j,k),STATECOMP_MOF+vofcomp_raw+dir)= &
                  mofdata(vofcomp_recon+dir)
              enddo
             enddo ! u_im=1..num_materials

             delta_mass(im_source)=delta_mass(im_source)+ &
              volgrid*(newvfrac(im_source)-oldvfrac(im_source))
             delta_mass(im_dest+num_materials)= &
              delta_mass(im_dest+num_materials)+ &
              volgrid*(newvfrac(im_dest)-oldvfrac(im_dest))

            else
             print *,"dF bust"
             stop
            endif

           else if ((dF.ge.zero).and. &
                    (dF.le.EBVOFTOL)) then
            ! do nothing
           else
            print *,"dF became corrupt2 dF=",dF
            print *,"dF invalid"
            stop
           endif 

          else if (max_velnode.eq.zero) then
           ! do nothing
          else
           print *,"max_velnode invalid"
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
         do im=1,num_materials-1
          do im_opp=im+1,num_materials

           call get_iten(im,im_opp,iten)

           do ireverse=0,1
            if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
             print *,"im or im_opp bust 10"
             stop
            endif

            LL=get_user_latent_heat(iten+ireverse*num_interfaces,room_temperature,1)

            if ((is_rigid(im).eq.1).or. &
                (is_rigid(im_opp).eq.1)) then
             ! do nothing
            else if ((is_rigid(im).eq.0).and. &
                     (is_rigid(im_opp).eq.0)) then

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
              local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
              mass_frac_id=mass_fraction_id(iten+ireverse*num_interfaces)

              if (LL.gt.zero) then ! evaporation
               im_condensed=im_source
              else if (LL.lt.zero) then ! condensation
               im_condensed=im_dest
              else
               print *,"LL invalid"
               stop
              endif

              if (is_GFM_freezing_modelF(local_freezing_model).eq.1) then

               Tsat_flag=NINT(TgammaFAB(D_DECL(i,j,k),iten))
               if (ireverse.eq.0) then
                ! do nothing
               else if (ireverse.eq.1) then
                Tsat_flag=-Tsat_flag
               else
                print *,"ireverse invalid"
                stop
               endif

              else if (is_GFM_freezing_modelF(local_freezing_model).eq.0) then

               Tsat_flag=0

              else
               print *,"is_GFM_freezing_modelF invalid"
               stop
              endif

              if (is_multi_component_evapF(local_freezing_model, &
                   Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces), &
                   LL).eq.0) then 
               ! do nothing
              else if (is_multi_component_evapF(local_freezing_model, &
                        Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces), &
                        LL).eq.1) then 

               if ((mass_frac_id.ge.1).and. &
                   (mass_frac_id.le.num_species_var)) then

                speccomp_mod=STATECOMP_STATES+ &
                 (im_condensed-1)*num_state_material+num_state_base+ &
                 mass_frac_id
          
                default_comp=(mass_frac_id-1)*num_materials+im_condensed 
                Ygamma_default=fort_speciesconst(default_comp)

                if ((Tsat_flag.eq.1).or.(Tsat_flag.eq.2)) then
                 Ygamma_default=TgammaFAB(D_DECL(i,j,k), &
                   num_interfaces+(iten-1)*ncomp_per_tsat+2)
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
                 speccomp_mod=STATECOMP_STATES+ &
                  (im_probe-1)*num_state_material+num_state_base+ &
                  mass_frac_id
                 mass_frac_limit=snew(D_DECL(i,j,k),speccomp_mod)
                 if ((mass_frac_limit.ge.-EPS_8_4).and. &
                     (mass_frac_limit.le.zero)) then
                  mass_frac_limit=zero
                 else if ((mass_frac_limit.ge.zero).and. &
                          (mass_frac_limit.le.one)) then
                  ! do nothing
                 else if ((mass_frac_limit.ge.one).and. &
                          (mass_frac_limit.le.one+EPS_8_4)) then
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
             print *,"is_rigid invalid MASS_TRANSFER_3D.F90"
             stop
            endif
           enddo ! ireverse=0,1
          enddo ! im_opp
         enddo ! im
        else if (is_rigid(im_primary).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(im_primary) invalid"
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
      end subroutine fort_convertmaterial

      subroutine fort_extend_burning_vel( &
        velflag, &
        level, &
        finest_level, &
        xlo,dx, &
        nburning, &
        ngrow, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        vel,DIMS(vel), &
        LS,DIMS(LS)) &  ! old level set function before phase change
      bind(c,name='fort_extend_burning_vel')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: velflag
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level 
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: nburning
      integer, INTENT(in) :: ngrow
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact

      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(LS)

      ! first num_interfaces components are the status.
      real(amrex_real), target, INTENT(inout) :: vel(DIMV(vel),nburning)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(in) :: LS(DIMV(LS),num_materials*(SDIM+1))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)

      integer im,im_opp
      integer iten,ireverse,sign_local
      integer im_dest,im_source
      integer i,j,k,dir
      integer i_sp,j_sp,k_sp
      integer sten_lo(3)
      integer sten_hi(3)
      integer indexcp(SDIM)
      integer tag_local
      real(amrex_real) rtag_local

      integer, parameter :: nhalf=3

      real(amrex_real) xijk(-nhalf:nhalf,SDIM)
      real(amrex_real) xcp(SDIM)
      real(amrex_real) xsp(-nhalf:nhalf,SDIM)
      real(amrex_real) nrm(SDIM)
      real(amrex_real) weight,total_weight
      real(amrex_real) vel_sum(SDIM)
      real(amrex_real) vel_local
      real(amrex_real) dxmax,dxmaxLS,eps,extensionwidth,LL
      real(amrex_real) ls_local(num_materials)
      integer im_local
      integer scomp
      integer ncomp_per

      vel_ptr=>vel
      LS_ptr=>LS

      if (velflag.eq.1) then
       ncomp_per=EXTRAP_PER_BURNING
      else if (velflag.eq.0) then
       ncomp_per=EXTRAP_PER_TSAT ! interface temperature, mass fraction
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
       print *,"expecting ngrow==4 in fort_extend_burning_vel"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"expecting ngrow_make_distance==3 in fort_extend_burning_vel"
       stop
      endif
      if (nburning.eq.num_interfaces*(ncomp_per+1)) then
       ! do nothing
      else
       print *,"nburning invalid"
       stop
      endif

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)
 
       ! Guard against division zero in the weight calculation
      eps=dxmaxLS*EPS4

      extensionwidth=dxmaxLS*ngrow_make_distance 

      call checkbound_array(fablo,fabhi,vel_ptr,ngrow_distance,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,ngrow,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do im=1,num_materials-1
       do im_opp=im+1,num_materials
        do ireverse=0,1
         if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
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
  
         call get_iten(im,im_opp,iten)
         LL=get_user_latent_heat(iten+ireverse*num_interfaces,room_temperature,1)

         if (LL.ne.zero) then
  
          ! iterate over domain and check tag
          do k=growlo(3),growhi(3)
          do j=growlo(2),growhi(2)
          do i=growlo(1),growhi(1)
           do im_local=1,num_materials
            ls_local(im_local)=LS(D_DECL(i,j,k),im_local)
           enddo
           call get_primary_material(ls_local,im_local)
           if (is_rigid(im_local).eq.0) then

            if ((im_dest.ge.1).and.(im_dest.le.num_materials)) then
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
                 nrm(dir)=LS(D_DECL(i,j,k),num_materials+(im_dest-1)*SDIM+dir)
                 xcp(dir)=xijk(0,dir)-ls_local(im_dest)*nrm(dir)
                enddo ! dir
               else if ((ls_local(im_dest).ge.zero).and. &
                        (ls_local(im_source).le.zero)) then
                do dir=1,SDIM
                 nrm(dir)=LS(D_DECL(i,j,k),num_materials+(im_source-1)*SDIM+dir)
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
                 nrm(dir)=LS(D_DECL(i,j,k),num_materials+(im_dest-1)*SDIM+dir)
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
                 nrm(dir)=LS(D_DECL(i,j,k),num_materials+(im_source-1)*SDIM+dir)
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
                sten_hi(dir)=indexcp(dir)+1
               enddo ! dir=1..sdim

               do dir=1,ncomp_per
                vel_sum(dir)=zero
               enddo

               total_weight = zero

               ! gather support point/weights
               do k_sp = sten_lo(3),sten_hi(3)
               do j_sp = sten_lo(2),sten_hi(2)
               do i_sp = sten_lo(1),sten_hi(1)
                call safe_data(i_sp,j_sp,k_sp,iten,vel_ptr,rtag_local)
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
                  scomp=num_interfaces+(iten-1)*ncomp_per+dir
                  call safe_data(i_sp,j_sp,k_sp,scomp,vel_ptr,vel_local)
                  vel_sum(dir) = vel_sum(dir)+weight*vel_local 
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
                 scomp=num_interfaces+(iten-1)*ncomp_per+dir
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
              print *,"num_materials,num_interfaces,nburning,ngrow ",num_materials,num_interfaces,nburning,ngrow
              stop
             endif
            else
             print *,"im_dest invalid"
             stop
            endif

           else if (is_rigid(im_local).eq.1) then
            ! do nothing
           else
            print *,"is_rigid(im_local) invalid"
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
       enddo ! im_opp=im+1..num_materials
      enddo ! im=1..num_materials-1

      return 
      end subroutine fort_extend_burning_vel

        ! ngrow corresponds to ngrow_distance
      subroutine fort_extend_drag( &
        level, &
        finest_level, &
        xlo,dx, &
        ncomp, &
        ngrow, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        drag,DIMS(drag), &
        LS,DIMS(LS)) &  
      bind(c,name='fort_extend_drag')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level 
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: ngrow
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact

      integer, INTENT(in) :: DIMDEC(drag)
      integer, INTENT(in) :: DIMDEC(LS)

      real(amrex_real), target, INTENT(inout) :: drag(DIMV(drag),ncomp)
      real(amrex_real), pointer :: drag_ptr(D_DECL(:,:,:),:)

      real(amrex_real), target, INTENT(in) :: LS(DIMV(LS),num_materials*(SDIM+1))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)

      integer im
      integer im_local
      integer i,j,k,dir
      integer drag_comp
      integer drag_im
      integer drag_type
      integer i_sp,j_sp,k_sp
      integer sten_lo(3)
      integer sten_hi(3)
      integer indexcp(SDIM)
      integer tag_local
      real(amrex_real) rtag_local

      integer, parameter :: nhalf=3
      real(amrex_real) xijk(-nhalf:nhalf,SDIM)
      real(amrex_real) xcp(SDIM)
      real(amrex_real) xsp(-nhalf:nhalf,SDIM)
      real(amrex_real) nrm(SDIM)
      real(amrex_real) weight,total_weight
      real(amrex_real) drag_sum(ncomp)
      real(amrex_real) drag_local(ncomp)
      real(amrex_real) dxmax,dxmaxLS,eps,extensionwidth
      real(amrex_real) ls_local(num_materials)

      drag_ptr=>drag
      LS_ptr=>LS

      if (bfact.lt.1) then
       print *,"bfact invalid118"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in extend_drag"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (ngrow.ne.4) then
       print *,"expecting ngrow==4 in fort_extend_drag"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"expecting ngrow_make_distance==3 in fort_extend_drag"
       stop
      endif
      if (ncomp.eq.N_DRAG) then
       ! do nothing
      else
       print *,"ncomp invalid"
       stop
      endif

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)
 
       ! Guard against division zero in the weight calculation
      eps=dxmaxLS*1.0E-4

      extensionwidth=dxmaxLS*ngrow_make_distance 

      call checkbound_array(fablo,fabhi,drag_ptr,ngrow_make_distance,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,ngrow,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       ! im is the material on which a force is applied.
      do im=1,num_materials
  
       ! iterate over domain and check tag
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        do im_local=1,num_materials
         ls_local(im_local)=LS(D_DECL(i,j,k),im_local)
        enddo
        rtag_local=drag(D_DECL(i,j,k),DRAGCOMP_FLAG+im)
        tag_local=NINT(rtag_local) 
         ! tag_local==0 if drag (flux) not yet init.
        if ((tag_local.eq.0).and.(rtag_local.eq.zero)) then
         ! if close enough to force recipient interface.
         if (abs(ls_local(im)).le.extensionwidth) then
          ! Project a point in normal direction of 
          ! the recipient level set function interface.
          call gridsten_level(xijk,i,j,k,level,nhalf)
          do dir=1,SDIM
           nrm(dir)=LS(D_DECL(i,j,k),num_materials+(im-1)*SDIM+dir)
           xcp(dir)=xijk(0,dir)-ls_local(im)*nrm(dir)
          enddo ! dir

          sten_lo(3)=0
          sten_hi(3)=0
          call containing_cell(bfact,dx,xlo,fablo,xcp,indexcp)
          do dir=1,SDIM
           sten_lo(dir)=indexcp(dir)-1
           sten_hi(dir)=indexcp(dir)+1
          enddo ! dir=1..sdim

          do drag_comp=1,ncomp
           drag_sum(drag_comp)=zero
          enddo

          total_weight = zero

          ! gather support point/weights
          do k_sp = sten_lo(3),sten_hi(3)
          do j_sp = sten_lo(2),sten_hi(2)
          do i_sp = sten_lo(1),sten_hi(1)
           do drag_comp=1,ncomp
            call safe_data(i_sp,j_sp,k_sp,drag_comp,drag_ptr, &
              drag_local(drag_comp))
           enddo
           rtag_local=drag_local(DRAGCOMP_FLAG+im) 
           tag_local=NINT(rtag_local) 

           if ((tag_local.eq.1).and. &
               (rtag_local.eq.one)) then
            call gridsten_level(xsp,i_sp,j_sp,k_sp,level,nhalf)

            weight=zero
            do dir=1,SDIM
             weight = weight + (xsp(0,dir)-xcp(dir))**two
            enddo
            weight=sqrt(weight)
            weight =one/((weight+eps)**four)
            total_weight = total_weight+weight
            do drag_comp=0,ncomp-1
             drag_type=fort_drag_type(drag_comp,drag_im)
             if (drag_im+1.eq.im) then
              if ((drag_type.ge.0).and. &
                  (drag_type.lt.DRAG_TYPE_NEXT).and. &
                  (drag_type.ne.DRAG_TYPE_FLAG)) then 
               drag_sum(drag_comp+1) = drag_sum(drag_comp+1)+ &
                 weight*drag_local(drag_comp+1)
              else if (drag_type.eq.DRAG_TYPE_FLAG) then
               if (drag_comp+1.eq.DRAGCOMP_FLAG+im) then
                ! do nothing
               else
                print *,"drag_comp invalid MASS_TRANSFER (1), drag_comp=", &
                   drag_comp
                stop
               endif
              else
               print *,"drag_type invalid"
               stop
              endif
             else if ((drag_im.ge.0).and.(drag_im.lt.num_materials)) then
              ! do nothing
             else
              print *,"drag_im invalid"
              stop
             endif
            enddo ! do drag_comp=0,ncomp-1
           else if ((tag_local.eq.0).and.(rtag_local.eq.zero)) then
            ! do nothing
           else if ((tag_local.eq.2).and.(rtag_local.eq.two)) then
            ! do nothing (this is a newly extended drag)
           else
            print *,"tag_local or rtag_local invalid1:", &
               tag_local,rtag_local
            stop
           endif
          enddo
          enddo
          enddo ! i_sp,j_sp,k_sp

          if (total_weight.gt.zero) then
           do drag_comp=0,ncomp-1
            drag_type=fort_drag_type(drag_comp,drag_im)
            if (drag_im+1.eq.im) then
             if ((drag_type.ge.0).and. &
                 (drag_type.lt.DRAG_TYPE_NEXT).and. &
                 (drag_type.ne.DRAG_TYPE_FLAG)) then 
              drag(D_DECL(i,j,k),drag_comp+1)=drag_sum(drag_comp+1)/ &
                      total_weight
             else if (drag_type.eq.DRAG_TYPE_FLAG) then
              if (drag_comp+1.eq.DRAGCOMP_FLAG+im) then
               ! do nothing
              else
               print *,"drag_comp invalid MASS_TRANSFER (2), drag_comp=", &
                   drag_comp
               stop
              endif
             else
              print *,"drag_type invalid"
              stop
             endif
            else if ((drag_im.ge.0).and.(drag_im.lt.num_materials)) then
             ! do nothing
            else
             print *,"drag_im invalid"
             stop
            endif
           enddo ! do drag_comp=0,ncomp-1
           drag(D_DECL(i,j,k),DRAGCOMP_FLAG+im)=two
          else if (total_weight.eq.zero) then
           ! do nothing
          else
           print *,"total_weight invalid"
           stop
          endif
         else if (abs(ls_local(im)).gt.extensionwidth) then
          ! do nothing
         else
          print *,"ls_local(im) invalid"
          stop
         endif
        else if ((tag_local.eq.1).and.(rtag_local.eq.one)) then
         ! do nothing
        else
         print *,"tag_local or rtag_local invalid2:", &
           tag_local,rtag_local
         print *,"i,j,k ",i,j,k
         print *,"level,finest_level ",level,finest_level
         print *,"num_materials,num_interfaces,ngrow ",num_materials,num_interfaces,ngrow
         stop
        endif

       enddo ! k
       enddo ! j
       enddo ! i
      enddo !im=1,num_materials

      return 
      end subroutine fort_extend_drag

      subroutine add_to_TI_YI( &
        T_history,Y_history,VEL_history, &
        mdot_diff_history, &
        TI_YI_ptr, &
        TI_YI_counter, &
        TI_YI_best_guess_index)
      IMPLICIT NONE
      real(amrex_real), INTENT(in), pointer :: TI_YI_ptr(:,:)
      integer, INTENT(inout) :: TI_YI_counter
      integer, INTENT(inout) :: TI_YI_best_guess_index
      real(amrex_real), INTENT(in) :: T_history
      real(amrex_real), INTENT(in) :: Y_history
      real(amrex_real), INTENT(in) :: VEL_history
      real(amrex_real), INTENT(in) :: mdot_diff_history
     
      if (TI_YI_counter.ge.MAX_TI_YI_logfile) then
       print *,"TI_YI_counter too large"
       stop
      else if ((TI_YI_counter.ge.0).and. &
               (TI_YI_counter.lt.MAX_TI_YI_logfile)) then
       if ((TI_YI_best_guess_index.ge.0).and. &
           (TI_YI_best_guess_index.le.TI_YI_counter)) then
        TI_YI_counter=TI_YI_counter+1
        TI_YI_ptr(TI_YI_counter,1)=T_history
        TI_YI_ptr(TI_YI_counter,2)=Y_history
        TI_YI_ptr(TI_YI_counter,3)=VEL_history
        TI_YI_ptr(TI_YI_counter,4)=mdot_diff_history
        if (TI_YI_best_guess_index.eq.0) then
         TI_YI_best_guess_index=TI_YI_counter
        else if ((TI_YI_best_guess_index.ge.1).and. &
                 (TI_YI_best_guess_index.le.TI_YI_counter-1)) then
         if (abs(mdot_diff_history).lt. &
             abs(TI_YI_ptr(TI_YI_best_guess_index,4))) then
          TI_YI_best_guess_index=TI_YI_counter
         else if (abs(mdot_diff_history).ge. &
                  abs(TI_YI_ptr(TI_YI_best_guess_index,4))) then
          ! do nothing
         else
          print *,"mdot_diff_history or TI_YI_ptr is NaN"
          stop
         endif
        else
         print *,"TI_YI_best_guess_index invalid"
         stop
        endif
       else
        print *,"TI_YI_best_guess_index invalid"
        stop
       endif
      else
       print *,"TI_YI_counter invalid"
       stop
      endif

      end subroutine add_to_TI_YI

      subroutine advance_TY_gamma( &
        prescribed_mdot, &
        TSAT_Y_PARMS, &
        POUT, &
        fully_saturated, &
        trial_and_error, &
        probe_ok, &
        TSAT_iter, &
        TSAT_predict, &
        TSAT_correct, &
        Y_predict, &
        VEL_correct, &
        Clausius_Clapyron_Tsat, & !intent(in)
        TSAT_base, & !intent(in)
        LL,RR,WA,WV, &
        TI_YI_ptr,TI_YI_counter,TI_YI_best_guess_index, &
        TI_min,TI_max, &
        T_gamma_a,T_gamma_b,T_gamma_c, &
        Y_gamma_a,Y_gamma_b,Y_gamma_c, &
        X_gamma_a,X_gamma_b,X_gamma_c)
      use global_utility_module
      use mass_transfer_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: prescribed_mdot
      real(amrex_real), INTENT(inout) :: X_gamma_a,X_gamma_b,X_gamma_c
      real(amrex_real), INTENT(inout) :: Y_gamma_a,Y_gamma_b,Y_gamma_c
      real(amrex_real), INTENT(inout) :: T_gamma_a,T_gamma_b,T_gamma_c
      real(amrex_real), INTENT(in) :: TI_min,TI_max
      real(amrex_real), INTENT(in), pointer :: TI_YI_ptr(:,:)
      integer, INTENT(inout) :: TI_YI_counter
      integer, INTENT(inout) :: TI_YI_best_guess_index
      real(amrex_real), INTENT(in) :: LL,RR,WA,WV
      real(amrex_real), INTENT(in) :: Clausius_Clapyron_Tsat
      real(amrex_real), INTENT(in) :: TSAT_base
      real(amrex_real), INTENT(inout) :: TSAT_predict
      real(amrex_real), INTENT(inout) :: TSAT_correct
      real(amrex_real), INTENT(inout) :: Y_predict
      real(amrex_real), INTENT(in) :: VEL_correct
      integer, INTENT(in) :: TSAT_iter
      integer, INTENT(in) :: probe_ok
      integer, INTENT(in) :: trial_and_error
      integer, INTENT(in) :: fully_saturated
      type(TSAT_MASS_FRAC_parm_type), INTENT(in) :: TSAT_Y_PARMS
      type(probe_out_type), INTENT(inout) :: POUT
      real(amrex_real) X_history
      real(amrex_real) T_history,Y_history,mdot_diff_history
      real(amrex_real) mdot_diff_a
      real(amrex_real) mdot_diff_b
      real(amrex_real) mdot_diff_c

      if (trial_and_error.eq.0) then
       T_gamma_c=half*(T_gamma_a+T_gamma_b)
      else if (trial_and_error.eq.1) then

       if (prescribed_mdot.eq.zero) then
        !do nothing
       else
        print *,"expecting prescribed_mdot.eq.zero"
        stop
       endif

       if (TSAT_iter.eq.0) then
        T_gamma_c=TSAT_predict
       else if ((TSAT_iter.gt.0).and. &
                (TSAT_iter.le.MAX_TI_YI_trials+1)) then
        if (TI_min.lt.TI_max) then
         T_gamma_c=TI_min+((TI_max-TI_min)/MAX_TI_YI_trials)* &
                 (TSAT_iter-1)
        else if (TI_min.eq.TI_max) then
         T_gamma_c=TI_min
        else
         print *,"TI_min or TI_max invalid"
         stop
        endif
       else
        print *,"TSAT_iter invalid"
        stop
       endif
      else
       print *,"trial_and_error invalid"
       stop
      endif

      T_history=TSAT_predict

      if (fully_saturated.eq.2) then
       X_gamma_c=one
       Y_gamma_c=one
       Y_history=one
      else if (fully_saturated.eq.0) then
       call X_from_Tgamma(X_gamma_c,T_gamma_c, &
         Clausius_Clapyron_Tsat, &
         LL,RR,WV) 
       call massfrac_from_volfrac(X_gamma_c,Y_gamma_c,WA,WV)

       call X_from_Tgamma(X_history,T_history, &
         Clausius_Clapyron_Tsat, &
         LL,RR,WV) 
       call massfrac_from_volfrac(X_history,Y_history,WA,WV)
      else
       print *,"fully_saturated invalid"
       stop
      endif

      call mdot_diff_func( &
       prescribed_mdot, &
       probe_ok,TSAT_Y_PARMS, &
       POUT, &
       Y_gamma_a,T_gamma_a,mdot_diff_a)

      call mdot_diff_func( &
       prescribed_mdot, &
       probe_ok,TSAT_Y_PARMS, &
       POUT, &
       Y_gamma_b,T_gamma_b,mdot_diff_b)

      call mdot_diff_func( &
       prescribed_mdot, &
       probe_ok,TSAT_Y_PARMS, &
       POUT, &
       Y_gamma_c,T_gamma_c,mdot_diff_c)

      call mdot_diff_func( &
       prescribed_mdot, &
       probe_ok,TSAT_Y_PARMS, &
       POUT, &
       Y_history,T_history, &
       mdot_diff_history)

      call add_to_TI_YI( &
       T_history,Y_history,VEL_correct, &
       mdot_diff_history, &
       TI_YI_ptr, &
       TI_YI_counter, &
       TI_YI_best_guess_index)

      if (trial_and_error.eq.0) then

       if (prescribed_mdot.eq.zero) then

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
         print *,"T_gamma a,b,c ", &
          T_gamma_a,T_gamma_b,T_gamma_c
         print *,"Y_gamma a,b,c ", &
          Y_gamma_a,Y_gamma_b,Y_gamma_c
         print *,"mdot_diff a,b,c ", &
          mdot_diff_a,mdot_diff_b,mdot_diff_c
         print *,"TSAT_iter=",TSAT_iter
         print *,"TI_min ",TI_min
         print *,"TI_max ",TI_max
         print *,"WA ",WA
         print *,"WV ",WV
         print *," LL ",LL
         print *," Clausius_Clapyron_Tsat ",Clausius_Clapyron_Tsat
         print *," TSAT_base ",TSAT_base
         stop
        endif

       else if (prescribed_mdot.gt.zero) then
        T_gamma_a=T_gamma_c
        T_gamma_b=T_gamma_c
        Y_gamma_a=Y_gamma_c
        Y_gamma_b=Y_gamma_c
       else
        print *,"prescribed_mdot invalid: ",prescribed_mdot
        stop
       endif

      else if (trial_and_error.eq.1) then
       T_gamma_a=T_gamma_c
       T_gamma_b=T_gamma_c
       Y_gamma_a=Y_gamma_c
       Y_gamma_b=Y_gamma_c
      else
       print *,"trial_and_error invalid"
       stop
      endif

      TSAT_correct=T_gamma_c
      Y_predict=Y_gamma_c

      end subroutine advance_TY_gamma
 
      ! vof,ref centroid,order,slope,intercept  x num_materials
      ! fort_fd_normal calls find_cut_geom_slope_CLSVOF:
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
      subroutine fort_ratemasschange( &
       tid, &
       nucleation_flag, &
       level, &
       finest_level, &
       ngrow_distance_in, &
       nstate, &
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
       prescribed_mdot, &
       observe_initial_mdot, &
       Tanasawa_or_Schrage_or_Kassemi, &
       interface_mass_transfer_model, &
       distribute_from_target, &
       mass_fraction_id, &
       constant_density_all_time, & ! 1..num_materials
       material_type_evap, &
       molar_mass, &
       species_molar_mass, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx, &
       prev_time, &
       dt, &
       arraysize, &
       blob_array, &
       color_count, &
       colorfab,DIMS(colorfab), &
       typefab,DIMS(typefab), &
       maskcov,DIMS(maskcov), &
       conductstate,DIMS(conductstate), &
       burnvel,DIMS(burnvel), &
       Tsatfab,DIMS(Tsatfab), &
       LS,DIMS(LS),  & !if nucleation_flag==0: localMF[HOLD_LS_DATA_MF]
       LSnew,DIMS(LSnew), & ! get_new_data(LS_Type,slab_step+1);
       Snew,DIMS(Snew), & 
       EOS,DIMS(EOS), &
       recon,DIMS(recon), &
       pres,DIMS(pres), &
       pres_eos,DIMS(pres_eos), &
       curvfab,DIMS(curvfab) ) &
      bind(c,name='fort_ratemasschange')

      use probf90_module
      use global_utility_module
      use MOF_routines_module
      use mass_transfer_module
      use hydrateReactor_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: nucleation_flag
      integer, target, INTENT(in) :: level
      integer, target, INTENT(in) :: finest_level
      integer :: local_probe_constrain
      integer :: user_override_TI_YI

      real(amrex_real), target :: TI_YI(MAX_TI_YI_logfile+1,4) !T,Y,VEL,mdot_diff
      real(amrex_real), pointer :: TI_YI_ptr(:,:)
      integer :: TI_YI_counter
      integer :: TI_YI_best_guess_index
      integer :: trial_and_error
      integer :: TI_YI_loop

      integer :: probe_ok
      real(amrex_real) :: microscale_probe_size
      integer, INTENT(in) :: ngrow_distance_in
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: nburning
      integer, INTENT(in) :: ntsat
      integer, INTENT(in) :: nden
      integer, INTENT(in) :: custom_nucleation_model
      integer, INTENT(in) :: do_the_nucleate
      integer, INTENT(in) :: nucleate_pos_size
      real(amrex_real), target, INTENT(in) :: nucleate_pos(nucleate_pos_size)
      real(amrex_real), target, INTENT(in) :: nucleation_temp(2*num_interfaces)
      real(amrex_real), target, INTENT(in) :: nucleation_pressure(2*num_interfaces)
      real(amrex_real), target, INTENT(in) :: nucleation_pmg(2*num_interfaces)
      real(amrex_real), target, INTENT(in) :: nucleation_mach(2*num_interfaces)
      real(amrex_real), target, INTENT(in) :: cavitation_pressure(num_materials)
      real(amrex_real), target, INTENT(in) :: cavitation_vapor_density(num_materials)
      real(amrex_real), target, INTENT(in) :: cavitation_tension(num_materials)
      integer, INTENT(in) ::  microlayer_substrate(num_materials)
      real(amrex_real), INTENT(in) :: microlayer_angle(num_materials)
      real(amrex_real), INTENT(in) :: microlayer_size(num_materials)
      real(amrex_real), INTENT(in) :: macrolayer_size(num_materials)
      real(amrex_real), INTENT(in) :: max_contact_line_size(num_materials)
      real(amrex_real), INTENT(in) :: R_Palmore_Desjardins
      integer, INTENT(in) :: use_exact_temperature(2*num_interfaces)
      real(amrex_real), INTENT(in) :: reaction_rate(2*num_interfaces)
      real(amrex_real) :: K_f(0:1)
      real(amrex_real), INTENT(in) :: hardwire_Y_gamma(2*num_interfaces)
      real(amrex_real), INTENT(in) :: hardwire_T_gamma(2*num_interfaces)
      real(amrex_real), INTENT(in) :: &
        accommodation_coefficient(2*num_interfaces)
      real(amrex_real), INTENT(in) :: reference_pressure(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp_curv(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp_vel(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp_min(2*num_interfaces)
      real(amrex_real), INTENT(in) :: saturation_temp_max(2*num_interfaces)
      integer, INTENT(in) :: freezing_model(2*num_interfaces)
      real(amrex_real), INTENT(in) :: prescribed_mdot(2*num_interfaces)
      integer, INTENT(in) :: observe_initial_mdot
      integer, INTENT(in) :: Tanasawa_or_Schrage_or_Kassemi(2*num_interfaces)
      integer, INTENT(in) :: interface_mass_transfer_model(2*num_interfaces)
      integer, INTENT(in) :: distribute_from_target(2*num_interfaces)
      integer, INTENT(in) :: mass_fraction_id(2*num_interfaces)
      integer, INTENT(in), target :: material_type_evap(num_materials)
      real(amrex_real), INTENT(in) :: molar_mass(num_materials)
      real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var+1)
      integer, INTENT(in) :: constant_density_all_time(num_materials)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, target, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, target, INTENT(in) :: bfact
      real(amrex_real), target, INTENT(in) :: xlo(SDIM)
      real(amrex_real), target, INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: prev_time
      real(amrex_real) :: cur_time
      real(amrex_real), INTENT(in) :: dt
      integer, INTENT(in) :: arraysize
      real(amrex_real), INTENT(in) :: blob_array(arraysize)
      integer, INTENT(in) :: color_count
      integer, INTENT(in) :: DIMDEC(colorfab)
      integer, INTENT(in) :: DIMDEC(typefab)
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(conductstate)
      integer, INTENT(in) :: DIMDEC(burnvel)
      integer, INTENT(in) :: DIMDEC(Tsatfab)
      integer, INTENT(in) :: DIMDEC(LS) ! declare the x,y,z dimensions of LS
      integer, INTENT(in) :: DIMDEC(LSnew)
      integer, INTENT(in) :: DIMDEC(Snew)
      integer, INTENT(in) :: DIMDEC(EOS)
      integer, INTENT(in) :: DIMDEC(recon)
      integer, INTENT(in) :: DIMDEC(pres)
      integer, INTENT(in) :: DIMDEC(pres_eos)
      integer, INTENT(in) :: DIMDEC(curvfab)

      real(amrex_real), INTENT(in), target :: typefab(DIMV(typefab))
      real(amrex_real), pointer :: typefab_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: colorfab(DIMV(colorfab))
      real(amrex_real), pointer :: colorfab_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov)) 
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in), target :: &
        conductstate(DIMV(conductstate),num_materials)
      real(amrex_real), pointer :: conductstate_ptr(D_DECL(:,:,:),:)

        ! destination vel: first num_interfaces components are the status.
      real(amrex_real), INTENT(out), target :: burnvel(DIMV(burnvel),nburning)
      real(amrex_real), pointer :: burnvel_ptr(D_DECL(:,:,:),:)

        ! TSAT : first num_interfaces components are the status.
      real(amrex_real), INTENT(out), target :: Tsatfab(DIMV(Tsatfab),ntsat)
      real(amrex_real), pointer :: Tsatfab_ptr(D_DECL(:,:,:),:)

        ! LS1,LS2,..,LSn,normal1,normal2,...normal_n 
        ! normal points from negative to positive
        !DIMV(LS)=x,y,z  
      real(amrex_real), target, INTENT(in) :: LS(DIMV(LS),num_materials*(SDIM+1)) 
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(inout) :: LSnew(DIMV(LSnew),num_materials*(SDIM+1))
      real(amrex_real), pointer :: LSnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(inout) :: Snew(DIMV(Snew),nstate)
      real(amrex_real), pointer :: Snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: EOS(DIMV(EOS),nden)
      real(amrex_real), pointer :: EOS_ptr(D_DECL(:,:,:),:)
       ! F,X,order,SL,I x num_materials
      real(amrex_real), target, INTENT(in) :: &
              recon(DIMV(recon),num_materials*ngeom_recon) 
      real(amrex_real), pointer :: recon_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: pres(DIMV(pres)) 
      real(amrex_real), pointer :: pres_ptr(D_DECL(:,:,:))
      real(amrex_real), target, INTENT(in) :: pres_eos(DIMV(pres_eos)) 
      real(amrex_real), pointer :: pres_eos_ptr(D_DECL(:,:,:))
      real(amrex_real), target, INTENT(in) :: &
          curvfab(DIMV(curvfab),2*(num_materials+num_interfaces)) 
      real(amrex_real), pointer :: curvfab_ptr(D_DECL(:,:,:),:)

      integer, parameter :: nhalf=3

      integer, target :: i,j,k
      integer dir,dir2
      integer im,im_opp,ireverse,iten
      integer imls
      integer im_ambient
      integer im_primary
      integer :: imls_I
      integer :: imls_I2
      integer im_substrate_source
      integer im_substrate_dest
      integer, target :: im_source,im_dest
      integer, target :: tcomp_source
      integer, target :: tcomp_dest
      integer, target :: Ycomp_source
      integer, target :: Ycomp_dest
      real(amrex_real) velmag_sum,local_velmag
      integer burnflag
      real(amrex_real) dxmin,dxmax
      real(amrex_real), target :: dxmaxLS
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), target :: xI(SDIM)
      real(amrex_real), target :: xsrc(SDIM)
      real(amrex_real), target :: xdst(SDIM)
      real(amrex_real), target :: xsrc_micro(SDIM)
      real(amrex_real), target :: xdst_micro(SDIM)
      real(amrex_real) nrmCP(SDIM)  ! closest point normal
      real(amrex_real) theta_nrmCP(SDIM)
      real(amrex_real), target :: LSINT(num_materials*(SDIM+1))
      real(amrex_real) LShere(num_materials)
      type(probe_out_type) :: POUT
      real(amrex_real) den_I_interp_SAT(2)
      real(amrex_real) local_hardwire_T(0:1)
      real(amrex_real) local_hardwire_Y(0:1)
      integer hardwire_flag(0:1)  ! =0 do not hardwire  =1 hardwire
      real(amrex_real) local_Tsat(0:1)
      real(amrex_real) Clausius_Clapyron_Tsat(0:1)
      real(amrex_real) delta_Tsat
      real(amrex_real) local_Tsat_base(0:1)
      real(amrex_real) vel_phasechange(0:1)
      real(amrex_real), target :: LL(0:1)
      real(amrex_real) :: latent_heat_temperature
      real(amrex_real), target :: dxprobe_source
      real(amrex_real), target :: dxprobe_dest
      real(amrex_real), target :: thermal_k_model_predict(2) ! source,dest
      real(amrex_real), target :: thermal_k_model_correct(2) ! source,dest
      real(amrex_real), target :: thermal_k_physical_base(2) ! source,dest
      real(amrex_real) LS_pos
      real(amrex_real) C_w0
      integer, target :: local_freezing_model
      integer local_Tanasawa_or_Schrage_or_Kassemi
      integer distribute_from_targ
      integer vofcomp_source,vofcomp_dest
      real(amrex_real) Fsource,Fdest ! used for the hydrate model
      real(amrex_real) LSSIGN,SIGNVEL
      integer found_path
      integer, target :: debugrate
      real(amrex_real) RR,mag
      integer for_estdt
      integer local_mask
      integer microlayer_substrate_source
      integer microlayer_substrate_dest
      real(amrex_real) gradphi_substrate(SDIM)
      real(amrex_real) newphi_substrate
      integer, target :: dencomp_source,dencomp_dest
      integer ispec
      real(amrex_real) vapor_den
      real(amrex_real) source_perim_factor,dest_perim_factor
      real(amrex_real) contact_line_perim
      integer icolor,base_type,ic,im1,im2
      real(amrex_real) normal_probe_factor
      integer iprobe_vapor
      integer im_vapor
      real(amrex_real), target :: Y_TOLERANCE
       ! iten=1..num_interfaces  ireverse=0..1
      real(amrex_real) temp_target_probe_history(2*num_interfaces,2)
      real(amrex_real) dxprobe_target_history(2*num_interfaces,2)

      integer use_tsatfab
      integer ncomp_per_burning
      integer ncomp_per_tsat

      real(amrex_real) CURV_OUT_I

      real(amrex_real) VEL_predict,VEL_correct
      real(amrex_real) X_predict
      real(amrex_real) Y_predict
      real(amrex_real) Y_predict_hold(0:1)
      real(amrex_real) TSAT_predict,TSAT_correct
      real(amrex_real) TSAT_ERR,TSAT_INIT_ERR
      integer TSAT_iter,TSAT_converge
      real(amrex_real) :: Y_interface_min
      real(amrex_real) :: TI_min
      real(amrex_real) :: TI_max
      real(amrex_real) FicksLawD(2)  ! iprobe=1 source iprobe=2 dest 
      real(amrex_real) molar_mass_ambient
      real(amrex_real) molar_mass_vapor
      integer interp_valid_flag_initial(2) ! iprobe=1 source iprobe=2 dest
      integer interface_resolved
      integer interp_status
      type(probe_parm_type), target :: PROBE_PARMS
      type(TSAT_MASS_FRAC_parm_type) :: TSAT_Y_PARMS
      type(nucleation_parm_type_input) :: create_in
      type(nucleation_parm_type_inout) :: create_inout
      integer iprobe
      integer im_probe
      integer microlayer_substrate_probe
      real(amrex_real) x_gamma_a,x_gamma_b,x_gamma_c
      real(amrex_real) Y_gamma_a,Y_gamma_b,Y_gamma_c
      real(amrex_real) T_gamma_a,T_gamma_b,T_gamma_c
      real(amrex_real) T_gamma_a_init,T_gamma_b_init
      real(amrex_real) Y_gamma_a_init,Y_gamma_b_init
      integer fully_saturated
      real(amrex_real) mdotT_debug
      real(amrex_real) mdotY_top_debug,mdotY_bot_debug,mdotY_debug
      real(amrex_real) TEMP_PROBE_source
      real(amrex_real) TEMP_PROBE_dest
      real(amrex_real) Y_PROBE_VAPOR

      integer debug_limiter

      TI_YI_ptr=>TI_YI

      maskcov_ptr=>maskcov
      conductstate_ptr=>conductstate
      burnvel_ptr=>burnvel
      Tsatfab_ptr=>Tsatfab
      LS_ptr=>LS
      LSnew_ptr=>LSnew
      Snew_ptr=>Snew
      curvfab_ptr=>curvfab
      EOS_ptr=>EOS
      pres_ptr=>pres
      pres_eos_ptr=>pres_eos

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

      if (nucleation_flag.eq.0) then
       ! do nothing
      else if (nucleation_flag.eq.1) then
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

      microscale_probe_size=EPS2

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

      if (ngrow_distance.ne.4) then
       print *,"expecting ngrow_distance==4 in fort_ratemasschange"
       stop
      endif
      if (ngrow_distance_in.ne.4) then
       print *,"expecting ngrow_distance_in==4 in fort_ratemasschange"
       stop
      endif

      if (ngrow_make_distance.ne.3) then
       print *,"expecting ngrow_make_distance==3 in fort_ratemasschange"
       stop
      endif

      ncomp_per_burning=EXTRAP_PER_BURNING
      if (nburning.eq.EXTRAP_NCOMP_BURNING) then
       ! do nothing
      else
       print *,"nburning invalid"
       stop
      endif
      ncomp_per_tsat=EXTRAP_PER_TSAT ! interface temperature, mass fraction
      if (ntsat.eq.EXTRAP_NCOMP_TSAT) then
       ! do nothing
      else
       print *,"ntsat invalid"
       stop
      endif
      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid in rate mass change"
       print *,"nden=",nden
       print *,"num_materials=",num_materials
       print *,"num_state_material=",num_state_material
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif
      
      if (nstate.eq.STATE_NCOMP) then
       ! do nothing
      else 
       print *,"nstate invalid"
       stop
      endif

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LSnew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,Snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,EOS_ptr,ngrow_distance,-1)
      call checkbound_array1(fablo,fabhi,pres_ptr,ngrow_distance,-1)
      call checkbound_array1(fablo,fabhi,pres_eos_ptr,1,-1)


      if (nucleation_flag.eq.0) then

       typefab_ptr=>typefab
       call checkbound_array1(fablo,fabhi,typefab_ptr,1,-1)
       colorfab_ptr=>colorfab
       call checkbound_array1(fablo,fabhi,colorfab_ptr,1,-1)

       call checkbound_array(fablo,fabhi,conductstate_ptr,1,-1)

       ! num_materials x (sdim+1) components
       call checkbound_array(fablo,fabhi, &
        burnvel_ptr, &
        ngrow_distance,-1)
       call checkbound_array(fablo,fabhi, &
        Tsatfab_ptr, &
        ngrow_distance,-1)

       call checkbound_array(fablo,fabhi, &
        curvfab_ptr, &
        ngrow_make_distance,-1)

       recon_ptr=>recon
       call checkbound_array(fablo,fabhi,recon_ptr,ngrow_distance,-1)
       call checkbound_array(fablo,fabhi,LS_ptr,ngrow_distance,-1)

      else if (nucleation_flag.eq.1) then
       ! do nothing
      else
       print *,"nucleation_flag invalid"
       stop
      endif

      if (arraysize.ne.num_elements_blobclass*color_count) then
       print *,"arraysize invalid rate mass change (get stat==1)"
       print *,"arraysize=",arraysize
       print *,"num_elements_blobclass=",num_elements_blobclass
       print *,"color_count=",color_count
       stop
      endif

       ! SANDIPAN HOOK HERE
       ! pseudo code:
       ! if typefab(D_DECL(i,j,k))=im_vapor then
       !  color = colorfab(D_DECL(i,j,k))
       !  vapor bubble statistics are in 
       !   blob_array((color-1)*num_elements_blobclass + l)
       !  l=1..num_elements_blobclass
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

      do iten=1,2*num_interfaces
       temp_target_probe_history(iten,1)=zero
       dxprobe_target_history(iten,1)=zero
       temp_target_probe_history(iten,2)=zero
       dxprobe_target_history(iten,2)=zero
      enddo
      create_in%EOS=>EOS
      create_in%Snew=>Snew
      create_in%LSnew=>LSnew
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

      PROBE_PARMS%tid=tid

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
      PROBE_PARMS%fablo=>fablo
      PROBE_PARMS%fabhi=>fabhi

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       create_in%i=i
       create_in%j=j
       create_in%k=k

       call gridsten_level(xsten,i,j,k,level,nhalf)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        if (nucleation_flag.eq.0) then

         ! LEVELSET FUNCTION AT CELL CENTERS.
         do im=1,num_materials
          LShere(im)=LS(D_DECL(i,j,k),im)
         enddo
         call get_primary_material(LShere,im_primary)

         if (is_rigid(im_primary).eq.0) then

          do im=1,num_materials-1
           do im_opp=im+1,num_materials

            vel_phasechange(0)=zero
            vel_phasechange(1)=zero

            if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
             print *,"im or im_opp bust 9"
             stop
            endif
            call get_iten(im,im_opp,iten)

            do ireverse=0,1

             temp_target_probe_history(iten+ireverse*num_interfaces,1)=zero
             dxprobe_target_history(iten+ireverse*num_interfaces,1)=zero
             temp_target_probe_history(iten+ireverse*num_interfaces,2)=zero
             dxprobe_target_history(iten+ireverse*num_interfaces,2)=zero

             latent_heat_temperature= &
              EOS(D_DECL(i,j,k), &
                  (im_primary-1)*num_state_material+ENUM_TEMPERATUREVAR+1)

             LL(ireverse)=get_user_latent_heat( &
               iten+ireverse*num_interfaces,latent_heat_temperature,0)

             K_f(ireverse)=reaction_rate(iten+ireverse*num_interfaces)

             local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
             local_Tanasawa_or_Schrage_or_Kassemi= &
               Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces)

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

             if (LL(ireverse).eq.zero) then
              im_vapor=im_dest ! placeholder
             else if (LL(ireverse).gt.zero) then
              im_vapor=im_dest
             else if (LL(ireverse).lt.zero) then
              im_vapor=im_source
             else
              print *,"LL invalid"
              stop
             endif

             ispec=mass_fraction_id(iten+ireverse*num_interfaces)
             if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
              molar_mass_vapor=species_molar_mass(ispec)
             else if (ispec.eq.0) then
              molar_mass_vapor=molar_mass(im_vapor)
             else
              print *,"ispec invalid"
              stop
             endif

             dencomp_source=(im_source-1)*num_state_material+1+ENUM_DENVAR
             dencomp_dest=(im_dest-1)*num_state_material+1+ENUM_DENVAR

             if ((ispec.ge.0).and.(ispec.le.num_species_var)) then
              ! do nothing
             else
              print *,"ispec invalid"
              stop
             endif

             vapor_den=fort_denconst(im_vapor) ! default

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

              if (is_multi_component_evapF(local_freezing_model, &
                   Tanasawa_or_Schrage_or_Kassemi( &
                   iten+ireverse*num_interfaces), &
                   LL(ireverse)).eq.0) then
               ! do nothing
              else if (is_multi_component_evapF(local_freezing_model, &
               Tanasawa_or_Schrage_or_Kassemi(iten+ireverse*num_interfaces), &
               LL(ireverse)).eq.1) then
               print *,"ispec invalid"
               stop
              endif
             else
              print *,"ispec invalid"
              stop
             endif

             distribute_from_targ= &
               distribute_from_target(iten+ireverse*num_interfaces)

             if ((distribute_from_targ.ne.0).and. &
                 (distribute_from_targ.ne.1)) then
              print *,"distribute_from_targ invalid"
              stop
             endif

              !Kassemi or Palmore and Desjardins
             if (local_freezing_model.eq.6) then 
              local_probe_constrain=0
              normal_probe_factor=one
             else if (is_valid_freezing_modelF(local_freezing_model).eq.1) then
              local_probe_constrain=1
              normal_probe_factor=half
             else
              print *,"local_freezing_model invalid"
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

              if ((is_rigid(im).eq.1).or. &
                  (is_rigid(im_opp).eq.1)) then

               ! do nothing

              else if (LL(ireverse).ne.zero) then

               local_hardwire_T(ireverse)= &
                  hardwire_T_gamma(iten+ireverse*num_interfaces)
               local_hardwire_Y(ireverse)= &
                  hardwire_Y_gamma(iten+ireverse*num_interfaces)
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

               local_Tsat(ireverse)= &
                 saturation_temp(iten+ireverse*num_interfaces)
               local_Tsat_base(ireverse)= &
                 saturation_temp(iten+ireverse*num_interfaces)
               Clausius_Clapyron_Tsat(ireverse)= &
                 saturation_temp(iten+ireverse*num_interfaces)

               if (1.eq.0) then
                print *,"(1)local_Tsat: iten,ireverse,T ", &
                   iten,ireverse,local_Tsat(ireverse)
               endif

               debug_limiter=0
               if (1.eq.0) then
                if ((level.eq.finest_level).and. &
                    (im_source.eq.1).and.(im_dest.eq.2).and. &
                    (i.eq.1).and.(j.eq.432)) then
                 debug_limiter=1
                endif
               endif

               found_path=0

                ! FOR YANG:
                ! BOTH LEVELSET FUNCTIONS WITHIN 2 dx of interface and
                ! at least one of them is positive.
                ! note: is_rigid(im_primary).eq.0 so we are not in 
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
                   nrmCP(dir)= &
                     LS(D_DECL(i,j,k),num_materials+(im_source-1)*SDIM+dir)
                     ! Least squares slope: see Sussman and Puckett (2000)
                   xI(dir)=xsten(0,dir)-LS_pos*nrmCP(dir)

                    ! normal_probe_factor=1/2 or 1
                   xdst(dir)=xI(dir)- &
                    normal_probe_factor*dxmin*nrmCP(dir) 
                   xsrc(dir)=xI(dir)+ &
                    normal_probe_factor*dxmin*nrmCP(dir)
                   xdst_micro(dir)=xI(dir)- &
                    microscale_probe_size*dxmin*nrmCP(dir) 
                   xsrc_micro(dir)=xI(dir)+ &
                    microscale_probe_size*dxmin*nrmCP(dir)
                  enddo ! dir=1..sdim
                 else if (LShere(im_source).ge.zero) then
                  LS_pos=LShere(im_dest)
                  do dir=1,SDIM
                   nrmCP(dir)= &
                     LS(D_DECL(i,j,k),num_materials+(im_dest-1)*SDIM+dir)

                   xI(dir)=xsten(0,dir)-LS_pos*nrmCP(dir)
                   xdst(dir)=xI(dir)+ &
                      normal_probe_factor*dxmin*nrmCP(dir) 
                   xsrc(dir)=xI(dir)- &
                      normal_probe_factor*dxmin*nrmCP(dir)
                   xdst_micro(dir)=xI(dir)+ &
                      microscale_probe_size*dxmin*nrmCP(dir) 
                   xsrc_micro(dir)=xI(dir)- &
                      microscale_probe_size*dxmin*nrmCP(dir)
                  enddo ! dir=1..sdim
                 else
                  print *,"LShere bust1"
                  print *,"LShere(im_dest) ",LShere(im_dest)
                  print *,"LShere(im_source) ",LShere(im_source)
                  stop
                 endif

                 found_path=1

                else if ((LShere(im_dest).lt.zero).and. &
                         (LShere(im_source).lt.zero)) then

                 found_path=0

                else
                 print *,"LShere bust2"
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
                 POUT%pres_I_interp(1)=2.0D+19 ! source
                 POUT%pres_I_interp(2)=2.0D+19 ! dest
                 POUT%Y_I_interp(1)=zero
                 POUT%Y_I_interp(2)=zero ! destination, C_methane_in_hyd

                 tcomp_source=(im_source-1)*num_state_material+ &
                   ENUM_TEMPERATUREVAR+1
                 tcomp_dest=(im_dest-1)*num_state_material+ &
                   ENUM_TEMPERATUREVAR+1

                 im_ambient=0
                 Ycomp_source=0
                 Ycomp_dest=0

                 if (LL(ireverse).gt.zero) then ! evaporation
                  im_ambient=im_dest
                 else if (LL(ireverse).lt.zero) then ! condensation
                  im_ambient=im_source
                 else
                  print *,"LL invalid"
                  stop
                 endif
                    
                 if (is_multi_component_evapF(local_freezing_model, &
                   Tanasawa_or_Schrage_or_Kassemi( &
                   iten+ireverse*num_interfaces),LL(ireverse)).eq.1) then
                  if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
                   Ycomp_source=(im_source-1)*num_state_material+2+ispec
                   Ycomp_dest=(im_dest-1)*num_state_material+2+ispec
                   FicksLawD(1)= &
                    fort_speciesviscconst((ispec-1)*num_materials+im_source)
                   FicksLawD(2)= &
                    fort_speciesviscconst((ispec-1)*num_materials+im_dest)
                  else
                   print *,"ispec invalid"
                   stop
                  endif
                 else if (is_multi_component_evapF(local_freezing_model, &
                  Tanasawa_or_Schrage_or_Kassemi( &
                  iten+ireverse*num_interfaces),LL(ireverse)).eq.0) then

                  if (ispec.eq.0) then
                   Ycomp_source=0
                   Ycomp_dest=0
                   FicksLawD(1)=zero
                   FicksLawD(2)=zero
                  else if ((ispec.ge.1).and.(ispec.le.num_species_var)) then
                   print *,"expecting ispec ==0 "
                   stop
                  endif
                 else
                  print *,"local_freezing_model invalid 3"
                  stop
                 endif

                  ! interp_status=1 success
                  ! interp_status=0 no proper stencil available
                 call interpfab_curv( &
                  interp_status, & !interp_status=0 if CURV_weight==0.0
                  num_materials+iten, &
                  bfact, &
                  level, &
                  finest_level, &
                  dx,xlo, &
                  xI, &
                  fablo,fabhi, &
                  curvfab_ptr, &
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
                   fablo,fabhi, &
                   Tsatfab_ptr, & ! not used since use_tsatfab==0
                   local_Tsat(ireverse), & !intent(out)
                   iten+ireverse*num_interfaces, &
                   saturation_temp, &
                   use_exact_temperature, &
                   xI,cur_time)
                 else if (hardwire_flag(ireverse).eq.1) then
                  local_Tsat(ireverse)=local_hardwire_T(ireverse)
                 else
                  print *,"hardwire_flag(ireverse) invalid"
                  stop
                 endif

                 if (1.eq.0) then
                  print *,"(2)local_Tsat: iten,ireverse,T ", &
                   iten,ireverse,local_Tsat(ireverse)
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
                 call prepare_normal(nrmCP,RR,mag,SDIM)

                 if (levelrz.eq.COORDSYS_CARTESIAN) then
                  ! do nothing
                 else if (levelrz.eq.COORDSYS_RZ) then
                  if (SDIM.ne.2) then
                   print *,"dimension bust"
                   stop
                  endif
                 else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
                  if (mag.gt.zero) then
                   do dir=1,SDIM
                    theta_nrmCP(dir)=nrmCP(dir)
                   enddo
                   RR=xsten(0,1)
                   call prepare_normal(theta_nrmCP,RR,mag,SDIM)
                   if (mag.gt.zero) then
                    ! mag=theta_nrmCP dot nrmCP
                    mag=zero
                    do dir=1,SDIM
                     mag=mag+theta_nrmCP(dir)*nrmCP(dir)
                    enddo
                    if (abs(mag).gt.one+EPS_8_4) then
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

                 do imls=1,num_materials*(SDIM+1)
                  call interpfab( &
                   bfact, &
                   level, &
                   finest_level, &
                   dx, &
                   xlo,xI, &
                   imls, &
                   fablo,fabhi, &
                   LS_ptr, &
                   LSINT(imls))
                 enddo ! imls=1..num_materials*(SDIM+1)

                  ! checks rigid materials and fluids (all materials)
                 call get_primary_material(LSINT,imls_I)
                 if (found_path.eq.1) then
                  if ((imls_I.eq.im_dest).or. &
                      (imls_I.eq.im_source)) then
                   call get_secondary_material(LSINT,imls_I,imls_I2)
                   if (imls_I2.ne.imls_I) then
                    if ((imls_I2.eq.im_dest).or. &
                        (imls_I2.eq.im_source)) then
                     ! do nothing
                    else 
                     found_path=0
                    endif
                   else
                    print *,"expecting imls_I2.ne.imls_I"
                    stop
                   endif
                  else
                   found_path=0
                  endif
                 else
                  print *,"expecting found_path.eq.1"
                  stop
                 endif

                 if (found_path.eq.1) then

                  TI_min=saturation_temp_min(iten+ireverse*num_interfaces)
                  TI_max=saturation_temp_max(iten+ireverse*num_interfaces)

                  if ((local_freezing_model.eq.5).or. & !Stefan evap/cond
                      (local_freezing_model.eq.6)) then !Palmore/Desjardins

                   if (LL(ireverse).gt.zero) then ! evaporation
                    TI_max=local_Tsat(ireverse)
                   else if (LL(ireverse).lt.zero) then ! condensation
                    TI_min=local_Tsat(ireverse)
                   else
                    print *,"LL(ireverse) invalid"
                    stop
                   endif

                  else if &
                    (is_valid_freezing_modelF(local_freezing_model).eq.1) then
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

                  PROBE_PARMS%probe_constrain=local_probe_constrain

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
                  PROBE_PARMS%im_source=>im_source
                  PROBE_PARMS%im_dest=>im_dest
                  PROBE_PARMS%tcomp_source=>tcomp_source
                  PROBE_PARMS%Ycomp_source=>Ycomp_source
                  PROBE_PARMS%dencomp_source=>dencomp_source
                  PROBE_PARMS%tcomp_dest=>tcomp_dest
                  PROBE_PARMS%Ycomp_dest=>Ycomp_dest
                  PROBE_PARMS%dencomp_dest=>dencomp_dest
                  PROBE_PARMS%xI=>xI

                  do iprobe=1,2

                   if (iprobe.eq.1) then
                    im_probe=im_source
                   else if (iprobe.eq.2) then
                    im_probe=im_dest
                   else
                    print *,"iprobe invalid"
                    stop
                   endif

                   thermal_k_model_predict(iprobe)= &
                    conductstate(D_DECL(i,j,k),im_probe)* &
                    fort_heatflux_factor(im_probe)
                   thermal_k_model_correct(iprobe)= &
                       thermal_k_model_predict(iprobe)
                   thermal_k_physical_base(iprobe)= &
                       get_user_heatviscconst(im_probe)

                   if (fort_heatflux_factor(im_probe).ge.zero) then
                    ! do nothing
                   else
                    print *,"fort_heatflux_factor(im_probe) invalid"
                    stop
                   endif
                   if (thermal_k_model_predict(iprobe).ge.zero) then
                    ! do nothing
                   else
                    print *,"thermal_k_model_predict(iprobe) invalid"
                    stop
                   endif
                   if (thermal_k_physical_base(iprobe).ge.zero) then
                    ! do nothing
                   else
                    print *,"thermal_k_physical_base(iprobe) invalid"
                    stop
                   endif

                  enddo ! iprobe=1,2


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

                   if (interp_status.eq.1) then
                    delta_Tsat= &
                     saturation_temp_curv(iten+ireverse*num_interfaces)* &
                     CURV_OUT_I
                   else if (interp_status.eq.0) then
                    delta_Tsat=zero
                   else
                    print *,"interp_status invalid"
                    stop
                   endif

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

                  if (1.eq.0) then
                   print *,"(3)local_Tsat: iten,ireverse,T ", &
                    iten,ireverse,local_Tsat(ireverse)
                  endif
              
                   ! initially probe_ok==0 since the probe values have
                   ! never been calculated at this point.
                   ! if local_probe_constrain==0, then the probes do not
                   ! have to be calculated again. 
                  probe_ok=0
  
                  !iprobe=1 source
                  !iprobe=2 dest
                  if (probe_ok.eq.1) then
                   ! do nothing
                  else if (probe_ok.eq.0) then
                   ! if local_probe_constrain==1, then
                   !  probe values depend on interface values.
                   ! if local_probe_constrain==0, then
                   !  probe values are insensitive to the interface values.
                   call probe_interpolation( &
                    PROBE_PARMS, & !intent(in)
                    local_Tsat(ireverse), & !intent(in)
                    Y_predict, & !intent(in)
                    POUT) !intent(out)
                  else
                   print *,"probe_ok invalid"
                   stop
                  endif

                  if (local_probe_constrain.eq.0) then
                   probe_ok=1 ! do not have to recompute the probes
                  else if (local_probe_constrain.eq.1) then
                   probe_ok=0 ! probes depend on TI, so they must be recalc.
                  else
                   print *,"local_probe_constrain invalid"
                   stop
                  endif

                  interp_valid_flag_initial(1)=POUT%interp_valid_flag(1)
                  interp_valid_flag_initial(2)=POUT%interp_valid_flag(2)

                  interface_resolved=1

                  if ((interp_valid_flag_initial(1).eq.1).and. &
                      (interp_valid_flag_initial(2).eq.1)) then
                   ! do nothing
                  else if ((interp_valid_flag_initial(1).eq.0).or. &
                           (interp_valid_flag_initial(2).eq.0)) then
                   interface_resolved=0
                  else
                   print *,"interp_valid_flag_initial invalid"
                   stop
                  endif

                  user_override_TI_YI=0

                  if (hardwire_flag(ireverse).eq.0) then

                   if (is_in_probtype_list().eq.1) then
                    ! do not call "mdot_diff_func" (below) if 
                    ! user_override_TI_YI=1
                    call SUB_INTERFACE_TEMPERATURE( &
                     interface_mass_transfer_model( &
                       iten+ireverse*num_interfaces), &
                     local_probe_constrain, &
                     ireverse, &
                     iten, &
                     xI, &
                     cur_time, &
                     prev_time, &
                     dt, &
                     local_Tsat(ireverse), & ! INTENT(inout)
                     Y_predict, &  ! INTENT(inout)
                     user_override_TI_YI, & ! INTENT(inout)
                     molar_mass, & ! index: 1..num_materials
                     species_molar_mass, & ! index: 1..num_species_var
                     thermal_k_model_predict(1), &
                     thermal_k_model_predict(2), & ! ksrc,kdst
                     thermal_k_physical_base(1), &
                     thermal_k_physical_base(2), & ! ksrc,kdst
                     POUT%T_probe(1), & ! source
                     POUT%T_probe(2), & ! dest
                     LL(ireverse), &
                     POUT%dxprobe_target(1), & ! source
                     POUT%dxprobe_target(2), & ! dest
                     num_materials, &
                     num_species_var)
                   endif

                  else if (hardwire_flag(ireverse).eq.1) then
                   ! do nothing
                  else
                   print *,"hardwire_flag invalid"
                   stop
                  endif

                  trial_and_error=0

                  if (local_probe_constrain.eq.0) then
                   if ((interface_resolved.eq.1).and. & ! valid probe
                       (interp_status.eq.1)) then ! valid curvature
                    ! do nothing
                   else if ((interface_resolved.eq.0).or. &
                            (interp_status.eq.0)) then
                    if (hardwire_flag(ireverse).eq.0) then
                     if (user_override_TI_YI.eq.0) then

                      if (1.eq.0) then
                       print *,"(3.8)local_Tsat: iten,ireverse,T ", &
                        iten,ireverse,local_Tsat(ireverse)
                      endif

                       ! TI must be in [TI_min,TI_max]
                      call apply_TI_limiter( &
                       TI_min,TI_max, & !intent(inout)
                       PROBE_PARMS, & !intent(in)
                       local_Tsat(ireverse), & !intent(inout)
                       Y_predict, & !intent(inout)
                       POUT) !intent(in)

                      if (1.eq.0) then
                       print *,"(3.9)local_Tsat: iten,ireverse,T ", &
                        iten,ireverse,local_Tsat(ireverse)
                      endif

                      if (prescribed_mdot(iten+ireverse*num_interfaces).eq. &
                          zero) then
                       trial_and_error=1
                      endif

                     else if (user_override_TI_YI.eq.1) then
                      ! do nothing
                     else
                      print *,"user_override_TI_YI invalid"
                      stop
                     endif

                    else if (hardwire_flag(ireverse).eq.1) then
                     ! do nothing
                    else
                     print *,"hardwire_flag(ireverse) invalid"
                     stop
                    endif
                   else
                    print *,"interface_resolved or interp_status invalid"
                    stop
                   endif
  
                  else if (local_probe_constrain.eq.1) then
                   ! do nothing; interface temperature does not depend on
                   ! probe values.
                  else
                   print *,"local_probe_constrain invalid"
                   stop
                  endif

                  if (1.eq.0) then
                   print *,"(4)local_Tsat: iten,ireverse,T ", &
                    iten,ireverse,local_Tsat(ireverse)
                  endif

                  TSAT_predict=local_Tsat(ireverse)

                  TSAT_correct=TSAT_predict

                  if (debug_limiter.eq.1) then
                   print *,"debug_limiter=1 (1) "
                   print *,"i,j,k,bfact,level,finest_level,im_source,im_dest ",&
                     i,j,k,bfact, &
                     level,finest_level,im_source,im_dest
                   print *,"ireverse, local_probe_constrain ",ireverse, &
                    local_probe_constrain
                   print *,"interp_status,interface_resolved ",interp_status, &
                    interface_resolved
                   print *,"TI_min,TI_max,local_Tsat(ireverse) ", &
                    TI_min,TI_max,local_Tsat(ireverse) 
                   print *,"hardwire_flag(ireverse) ",hardwire_flag(ireverse)
                   print *,"user_override_TI_YI ",user_override_TI_YI
                   print *,"T_I_interp(1) ",POUT%T_I_interp(1)
                   print *,"T_I_interp(2) ",POUT%T_I_interp(2)
                   print *,"T_probe(1) ",POUT%T_probe(1)
                   print *,"T_probe(2) ",POUT%T_probe(2)
                   print *,"POUT%interp_valid_flag(1) ", &
                     POUT%interp_valid_flag(1)
                   print *,"POUT%interp_valid_flag(2) ", &
                     POUT%interp_valid_flag(2)
                  endif

                  VEL_predict=zero
                  VEL_correct=zero

                  TSAT_iter=0
                  TSAT_converge=0

                  TI_YI_counter=0
                  TI_YI_best_guess_index=0

!    local_freezing_model=5  fully saturated evaporation?
!    local_freezing_model=6  Palmore/Desjardins
!    local_freezing_model=7  Cavitation (a seed must exist)

                  do while (TSAT_converge.eq.0) 

                   TSAT_Y_PARMS%Tanasawa_or_Schrage_or_Kassemi= &
                         local_Tanasawa_or_Schrage_or_Kassemi
                   TSAT_Y_PARMS%accommodation_coefficient= &
                         accommodation_coefficient(iten+ireverse*num_interfaces) 

                   TSAT_Y_PARMS%reference_pressure= &
                         reference_pressure(iten+ireverse*num_interfaces) 
                   if (TSAT_Y_PARMS%reference_pressure.gt.zero) then
                    ! do nothing
                   else
                    print *,"expecting TSAT_Y_PARMS%reference_pressure>0"
                    print *,"set ns.reference_pressure accordingly"
                    print *,"iten= ",iten
                    print *,"ireverse= ",ireverse
                    print *,"TSAT_Y_PARMS%reference_pressure=", &
                      TSAT_Y_PARMS%reference_pressure
                    stop
                   endif

                   TSAT_Y_PARMS%universal_gas_constant_R= &
                         R_Palmore_Desjardins
                   TSAT_Y_PARMS%molar_mass_ambient=molar_mass_ambient
                   TSAT_Y_PARMS%molar_mass_vapor=molar_mass_vapor
                   TSAT_Y_PARMS%iprobe_vapor=iprobe_vapor
                   TSAT_Y_PARMS%material_type_evap=>material_type_evap

                   !iprobe=1 source
                   !iprobe=2 dest
                   if (probe_ok.eq.1) then
                    ! do nothing (probes do not depend on TI)
                   else if (probe_ok.eq.0) then ! probes depend on TI
                    call probe_interpolation( &
                     PROBE_PARMS, &
                     TSAT_predict, &
                     Y_predict, &
                     POUT)
                   else
                    print *,"probe_ok invalid"
                    stop
                   endif

                   if (TSAT_iter.eq.0) then
                    fully_saturated=0

                    T_gamma_a=TSAT_predict
                    T_gamma_b=TSAT_predict
                    Y_gamma_a=Y_predict
                    Y_gamma_b=Y_predict

                    if (hardwire_flag(ireverse).eq.0) then

                     if (user_override_TI_YI.eq.0) then

                       ! 3=Kassemi model
                      if (local_Tanasawa_or_Schrage_or_Kassemi.eq.3) then
                       fully_saturated=2
                       Y_predict=one
                       X_predict=one
                       Y_gamma_a=Y_predict
                       Y_gamma_b=Y_predict

                       T_gamma_a=TSAT_predict
                       T_gamma_b=TSAT_predict
                      else if ((POUT%Y_probe(iprobe_vapor).ge. &
                                one-Y_TOLERANCE).and. &
                               (POUT%Y_probe(iprobe_vapor).le.one).and. &
                               (Y_predict.eq.one)) then
                       fully_saturated=1
                      else if ((POUT%Y_probe(iprobe_vapor).le. &
                                one-Y_TOLERANCE).and. &
                               (POUT%Y_probe(iprobe_vapor).ge.zero).and. &
                               (Y_predict.eq.one)) then
                        ! declared in PROBCOMMON.F90
                       Y_predict=one-EVAP_BISECTION_TOL
                       call volfrac_from_massfrac(X_predict,Y_predict, &
                        molar_mass_ambient,molar_mass_vapor) ! WA,WV
                       call Tgamma_from_TSAT_and_X( &
                        TSAT_predict, & !intent(out)
                        Clausius_Clapyron_Tsat(ireverse), & !intent(in)
                        X_predict,LL(ireverse), & !intent(in)
                        R_Palmore_Desjardins, & !intent(in)
                        molar_mass_vapor,TI_min,TI_max) !intent(in)

                       T_gamma_a=TSAT_predict
                       T_gamma_b=TSAT_predict
                       Y_gamma_a=Y_predict
                       Y_gamma_b=Y_predict
                      else
                       print *,"mass fraction bust"
                       stop
                      endif

                      if (fully_saturated.eq.1) then ! Tgamma=Tboil
                       ! do nothing
                      else if (fully_saturated.eq.2) then ! Kassemi
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
                         Clausius_Clapyron_Tsat(ireverse), &
                         LL(ireverse),R_Palmore_Desjardins, &
                         molar_mass_vapor) ! WV
                        call massfrac_from_volfrac(X_gamma_a,Y_gamma_a, &
                         molar_mass_ambient,molar_mass_vapor) ! WA,WV
                       else if (LL(ireverse).lt.zero) then ! condensation
                        T_gamma_b=TI_max
                        call X_from_Tgamma(X_gamma_b,T_gamma_b, &
                         Clausius_Clapyron_Tsat(ireverse), &
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

                     else if (user_override_TI_YI.eq.1) then

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
                      print *,"user_override_TI_YI invalid"
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

                   den_I_interp_SAT(1)=POUT%den_I_interp(1) ! iprobe=1 source
                   den_I_interp_SAT(2)=POUT%den_I_interp(2) ! iprobe=2 dest

                   !iprobe=1 source
                   !iprobe=2 dest

                   microlayer_substrate_source=0
                   microlayer_substrate_dest=0

                   im_substrate_source=microlayer_substrate(im_source)
                   im_substrate_dest=microlayer_substrate(im_dest)

                   if (microlayer_size(im_source).gt.zero) then
                    if (im_substrate_source.gt.0) then
                     if (is_rigid(im_substrate_source).ne.1) then
                      print *,"is_rigid(im_substrate_source).ne.1"
                      stop
                     endif
                     do dir=1,SDIM
                      dir2=num_materials+(im_substrate_source-1)*SDIM+dir
                      gradphi_substrate(dir)=LSINT(dir2)
                     enddo
                     call get_physical_dist(xI,LSINT(im_substrate_source), &
                      gradphi_substrate,newphi_substrate)
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
                     if (is_rigid(im_substrate_dest).ne.1) then
                      print *,"is_rigid(im_substrate_dest).ne.1"
                      stop
                     endif
                     do dir=1,SDIM
                      dir2=num_materials+(im_substrate_dest-1)*SDIM+dir
                      gradphi_substrate(dir)=LSINT(dir2)
                     enddo
                     call get_physical_dist(xI,LSINT(im_substrate_dest), &
                      gradphi_substrate,newphi_substrate)
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
                   else if &
                    (is_valid_freezing_modelF(local_freezing_model).eq.1) then
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
                        (microlayer_substrate_probe.le.num_materials)) then
                     icolor=NINT(colorfab(D_DECL(i,j,k)))
                     if ((icolor.gt.color_count).or.(icolor.le.0)) then
                      print *,"icolor invalid in fort_ratemasschange icolor=",&
                        icolor
                      print *,"i,j,k ",i,j,k
                      stop
                     endif
                     base_type=NINT(typefab(D_DECL(i,j,k)))
                     if ((base_type.lt.1).or.(base_type.gt.num_materials)) then
                      print *,"base_type invalid"
                      stop
                     endif

                      ! ic+1 is blob_triple_perim index
                     ic=(icolor-1)*num_elements_blobclass+BLB_TRIPLE_PERIM

                     im2=microlayer_substrate_probe
                     if ((im2.ge.1).and.(im2.le.num_materials)) then
                      if (base_type.eq.im_source) then
                       im1=im_dest
                      else if (base_type.eq.im_dest) then
                       im1=im_source
                      else
                       print *,"base_type invalid"
                       stop
                      endif
                      ic=ic+(im1-1)*num_materials+im2
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
                    else if &
                      (is_valid_freezing_modelF(local_freezing_model).eq.1) &
                      then
                     if ((max_contact_line_size(im_probe).eq.zero).or. &
                         (microlayer_size(im_probe).eq.zero).or. &
                         (macrolayer_size(im_probe).eq.zero).or. &
                         (microlayer_substrate_probe.eq.0)) then
                      ! do nothing
                     else
                      print *,"microlayer parameters invalid"
                      stop
                     endif
                    else
                     print *,"is_valid_freezing_modelF invalid"
                     stop
                    endif
        
                   enddo ! iprobe=1,2


                    ! V dt * L * L = volume of material change of phase
                    ! rho_src * V dt L^2 = rho_dst * (V+Vexpand) * dt L^2
                    ! S=V+Vexpand=
                    !  rho_src V/rho_dst=[k grad T dot n]/(L rho_dst)
                    ! Vexpand=(rho_src/rho_dst-1)V
                    ! note: rho_src V = MDOT

                    ! this call to get_vel_phasechange is not for the purpose
                    ! of estimating the timestep dt.
                   for_estdt=0

                   if (is_in_probtype_list().eq.1) then
                    call SUB_K_EFFECTIVE( &
                     interface_mass_transfer_model( &
                       iten+ireverse*num_interfaces), &
                     ireverse, &
                     iten, &
                     molar_mass, & ! index: 1..num_materials
                     species_molar_mass, & ! index: 1..num_species_var
                     thermal_k_model_predict, &
                     thermal_k_model_correct, &
                     thermal_k_physical_base, &
                     POUT%T_probe(1), & ! source
                     POUT%T_probe(2), & ! dest
                     POUT%dxprobe_target(1), & ! source
                     POUT%dxprobe_target(2), & ! dest
                     LL(ireverse), &
                     num_materials, &
                     num_species_var)
                   endif

                    ! if local_freezing_model==0 stefan problem
                    !  5, some kind of evaporation model,
                    !  or 1, then
                    !  DTsrc=(Tsrc-TSAT_predict)
                    !  DTdst=(Tdst-TSAT_predict)
                    !  velsrc=ksrc*DTsrc/(LL * dxprobe_src)
                    !  veldst=kdst*DTdst/(LL * dxprobe_dest)
                    !  in: PROB.F90
                   call get_vel_phasechange( &
                     interface_mass_transfer_model( &
                       iten+ireverse*num_interfaces), &
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
                     POUT%den_probe(1), & ! source 
                     POUT%den_probe(2), & ! dest
                     thermal_k_model_correct(1), &
                     thermal_k_model_correct(2), & ! ksrc,kdst
                     thermal_k_physical_base(1), &
                     thermal_k_physical_base(2), & ! ksrc,kdst
                     POUT%T_probe(1), & ! source
                     POUT%T_probe(2), & ! dest
                     TSAT_predict, &
                     POUT%T_I_interp(1), & !source
                     POUT%T_I_interp(2), & !dest
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
                     POUT%dxprobe_target(1), & ! source
                     POUT%dxprobe_target(2), & ! dest
                     im_source,im_dest, &
                     prev_time,dt, &
                     fort_alpha(iten+ireverse*num_interfaces), &
                     fort_beta(iten+ireverse*num_interfaces), &
                     fort_expansion_factor(iten+ireverse*num_interfaces), &
                     K_f(ireverse), &
                     POUT%Y_I_interp(2), & ! Cmethane_in_hydrate (dest)
                     C_w0, &
                     POUT%pres_I_interp(1), & ! PHYDWATER
                     Fsource,Fdest)

                   if (VEL_correct.ge.zero) then
                    ! do nothing
                   else if (VEL_correct.lt.zero) then
                    VEL_correct=zero
                   else
                    print *,"VEL_correct bust"
                    stop
                   endif

                   TSAT_correct=TSAT_predict

                   if (user_override_TI_YI.eq.0) then

                    if (local_freezing_model.eq.6) then ! Palmore/Desjardins

                     ! type(TSAT_MASS_FRAC_parm_type)
                     Y_interface_min=zero
                     TSAT_Y_PARMS%PROBE_PARMS=>PROBE_PARMS
                     TSAT_Y_PARMS%YI_min=Y_interface_min
                     TSAT_Y_PARMS%TI_min=TI_min
                     TSAT_Y_PARMS%TI_max=TI_max
                     TSAT_Y_PARMS%Clausius_Clapyron_Tsat= &
                         Clausius_Clapyron_Tsat(ireverse)
                     TSAT_Y_PARMS%TSAT_base=local_Tsat(ireverse)
                     TSAT_Y_PARMS%D_MASS=FicksLawD(iprobe_vapor)
                     TSAT_Y_PARMS%den_G=den_I_interp_SAT(iprobe_vapor)
                     TSAT_Y_PARMS%thermal_k=>thermal_k_model_correct

                     if (hardwire_flag(ireverse).eq.0) then

                       ! Kassemi
                      if (fully_saturated.eq.2) then

                       call advance_TY_gamma( &
                        prescribed_mdot(iten+ireverse*num_interfaces), &
                        TSAT_Y_PARMS, &
                        POUT, &
                        fully_saturated, &
                        trial_and_error, &
                        probe_ok, &
                        TSAT_iter, &
                        TSAT_predict, &
                        TSAT_correct, &
                        Y_predict, &
                        VEL_correct, &
                        Clausius_Clapyron_Tsat(ireverse), & !intent(in)
                        local_Tsat(ireverse), & !intent(in)
                        LL(ireverse),R_Palmore_Desjardins, &
                        molar_mass_ambient,molar_mass_vapor, &
                        TI_YI_ptr,TI_YI_counter,TI_YI_best_guess_index, &
                        TI_min,TI_max, &
                        T_gamma_a,T_gamma_b,T_gamma_c, &
                        Y_gamma_a,Y_gamma_b,Y_gamma_c, &
                        X_gamma_a,X_gamma_b,X_gamma_c)

                       ! Palmore and Desjardins, or
                       ! Kassemi
                      else if ((molar_mass_ambient.gt.zero).and. &
                               (molar_mass_vapor.gt.zero).and. &
                               (R_Palmore_Desjardins.gt.zero).and. &
                               (TSAT_Y_PARMS%den_G.gt.zero).and. &
                               ((fully_saturated.eq.0).or. &
                                (fully_saturated.eq.1))) then

                       ! Palmore and Desjardins, Y=X=1
                       if (fully_saturated.eq.1) then 

                        trial_and_error=0
                        Y_interface_min=one
                        Y_predict=one
                        TSAT_correct=local_Tsat(ireverse)

                       else if ((fully_saturated.eq.0).and. &
                                (POUT%Y_probe(iprobe_vapor).ge. &
                                 one-Y_TOLERANCE).and. &
                                (POUT%Y_probe(iprobe_vapor).le.one)) then

                        trial_and_error=0
                        Y_interface_min=one
                        Y_predict=one
                        TSAT_correct=local_Tsat(ireverse)

                       else if ((fully_saturated.eq.0).and. &
                                (TSAT_Y_PARMS%D_MASS.eq.zero)) then

                        trial_and_error=0
                        Y_interface_min=one
                        Y_predict=one
                        TSAT_correct=local_Tsat(ireverse)

                       else if ((fully_saturated.eq.0).and. &
                                (POUT%Y_probe(iprobe_vapor).le. &
                                 one-Y_TOLERANCE).and. &
                                (POUT%Y_probe(iprobe_vapor).ge.zero).and. &
                                (TSAT_Y_PARMS%D_MASS.gt.zero)) then

                        if (fully_saturated.eq.0) then

                         call advance_TY_gamma( &
                          prescribed_mdot(iten+ireverse*num_interfaces), &
                          TSAT_Y_PARMS, &
                          POUT, &
                          fully_saturated, &
                          trial_and_error, &
                          probe_ok, &
                          TSAT_iter, &
                          TSAT_predict, &
                          TSAT_correct, &
                          Y_predict, &
                          VEL_correct, &
                          Clausius_Clapyron_Tsat(ireverse), & !intent(in)
                          local_Tsat(ireverse), & !intent(in)
                          LL(ireverse),R_Palmore_Desjardins, &
                          molar_mass_ambient,molar_mass_vapor, &
                          TI_YI_ptr,TI_YI_counter,TI_YI_best_guess_index, &
                          TI_min,TI_max, &
                          T_gamma_a,T_gamma_b,T_gamma_c, &
                          Y_gamma_a,Y_gamma_b,Y_gamma_c, &
                          X_gamma_a,X_gamma_b,X_gamma_c)

                        else 
                         print *,"expecting fully_saturated==0"
                         stop
                        endif

                       else
                        print *,"fully_saturated or Y_probe invalid"
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
                      trial_and_error=0
                     else
                      print *,"hardwire_flag(ireverse) invalid"
                      stop
                     endif
  
                     call mdot_from_T_probe( &
                      prescribed_mdot(iten+ireverse*num_interfaces), &
                      probe_ok, &
                      TSAT_Y_PARMS, &
                      POUT, &
                      TSAT_correct,Y_predict, &
                      mdotT_debug, &
                      TEMP_PROBE_source, &
                      TEMP_PROBE_dest)

                     call mdot_from_Y_probe( &
                      prescribed_mdot(iten+ireverse*num_interfaces), &
                      probe_ok, &
                      TSAT_Y_PARMS, &
                      POUT, &
                      Y_predict,TSAT_correct, &
                      mdotY_top_debug, &
                      mdotY_bot_debug, &
                      mdotY_debug, &
                      Y_PROBE_VAPOR)

                    else if &
                      (is_valid_freezing_modelF(local_freezing_model).eq.1) &
                     then

                     trial_and_error=0
                     mdotT_debug=zero
                     mdotY_top_debug=zero
                     mdotY_bot_debug=zero
                     mdotY_debug=zero

                    else
                     print *,"local_freezing_model invalid 7"
                     stop
                    endif

                   else if (user_override_TI_YI.eq.1) then
                    trial_and_error=0
                    mdotT_debug=zero
                    mdotY_top_debug=zero
                    mdotY_bot_debug=zero
                    mdotY_debug=zero
                   else
                    print *,"(user_override_TI_YI invalid)"
                    stop
                   endif

                   if (hardwire_flag(ireverse).eq.0) then

                    if (user_override_TI_YI.eq.0) then

                     if ((interface_resolved.eq.1).and. &
                         (interp_status.eq.1)) then
                      delta_Tsat= &
                       saturation_temp_vel(iten+ireverse*num_interfaces)* &
                       (VEL_correct-VEL_predict)
                     else if ((interface_resolved.eq.0).or. &
                              (interp_status.eq.0)) then
                      delta_Tsat=zero
                     else
                      print *,"interface_resolved or interp_status invalid"
                      stop
                     endif

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

                    else if (user_override_TI_YI.eq.1) then
                     delta_Tsat=zero
                    else
                     print *,"user_override_TI_YI invalid"
                     stop
                    endif

                   else if (hardwire_flag(ireverse).eq.1) then
                    delta_Tsat=zero
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
                   if (1.eq.0) then
                    print *,"TSAT:it,Tnp1,Tn,mdot,Tsrc,Tdst ", &
                      TSAT_iter,TSAT_correct,TSAT_predict,mdotT_debug, &
                      TEMP_PROBE_source,TEMP_PROBE_dest
                   endif

                   TSAT_converge=0
                    ! TSAT_iter starts at 0
                   TSAT_iter=TSAT_iter+1

                   if (trial_and_error.eq.0) then 

                    if (TSAT_ERR.eq.zero) then
                     TSAT_converge=1
                    endif
                    if (TSAT_iter.gt.EVAPORATION_iter_max) then
                     TSAT_converge=1
                    endif
                    if (TSAT_iter.gt.1) then
                     if (TSAT_ERR.lt.EVAPORATION_TOL*TSAT_INIT_ERR) then
                      TSAT_converge=1
                     endif
                    endif

                   else if (trial_and_error.eq.1) then 

                    if ((TSAT_iter.ge.1).and. &
                        (TSAT_iter.le.MAX_TI_YI_trials+1)) then
                     ! do nothing
                    else if (TSAT_iter.eq.MAX_TI_YI_trials+2) then
                     TSAT_converge=1
                     TSAT_correct=TI_YI(TI_YI_best_guess_index,1)
                     Y_predict=TI_YI(TI_YI_best_guess_index,2)
                     VEL_correct=TI_YI(TI_YI_best_guess_index,3)
                    else 
                     print *,"TSAT_iter invalid"
                     stop
                    endif

                   else
                    print *,"trial_and_error invalid"
                    stop
                   endif

                   if (DEBUG_EVAPORATION.eq.1) then
                    print *,"DEBUG_EVAPORATION STATEMENT 3"
                    print *,"i,j,k,TSAT_iter,TSAT_ERR ", &
                     i,j,k,TSAT_iter,TSAT_ERR
                    print *,"TSAT_correct ",TSAT_correct
                    print *,"Y_predict ",Y_predict
                   endif

                  enddo ! do while (TSAT_converge.eq.0)

                  if (observe_initial_mdot.eq.1) then
                   if (user_override_TI_YI.eq.0) then
                    if (local_freezing_model.eq.6) then ! Palmore/Desjardins
                     if (hardwire_flag(ireverse).eq.0) then
                      print *,"i,j,k,iten,ireverse ", &
                        i,j,k,iten,ireverse
                      print *,"mdotT_debug ",mdotT_debug
                      print *,"mdotY_debug ",mdotY_debug
                      print *,"TEMP_PROBE_source ",TEMP_PROBE_source
                      print *,"TEMP_PROBE_dest ",TEMP_PROBE_dest
                      print *,"Clausius_Clapyron_Tsat ", &
                        Clausius_Clapyron_Tsat(ireverse)
                      print *,"Y_PROBE_VAPOR ",Y_PROBE_VAPOR
                      print *,"Y_predict (Y_interface) ",Y_predict
                      print *,"TSAT_correct (T_interface) ",TSAT_correct
                      print *,"n dot rhoG(U_V-U_I)=n dot rhoL(U_L-U_I)=mdot"
                      print *,"VEL_correct (U_I-U_L or U_I-U_V) ",VEL_correct
                      print *,"dxprobe_target(1) (source): ", &
                        POUT%dxprobe_target(1)
                      print *,"dxprobe_target(2) (dest): ", &
                        POUT%dxprobe_target(2)
                     endif
                    endif
                   endif
                  endif

                     
                  if (TSAT_iter.eq.1) then
                   ! check nothing
                  else if (TSAT_iter.gt.1) then

                   if (debug_limiter.eq.1) then
                    print *,"debug_limiter=1 (2) "
                    print *,"TSAT_correct=",TSAT_correct
                    print *,"TI_YI(TI_YI_best_guess_index,1)=", &
                      TI_YI(TI_YI_best_guess_index,1)
                    print *,"TI_YI_best_guess_index=", &
                      TI_YI_best_guess_index
                    do TI_YI_loop=1,TI_YI_counter
                     print *,"TI_YI idx,T,Y,VEL,mdotdiff ", &
                      TI_YI_loop, &
                      TI_YI(TI_YI_loop,1), & 
                      TI_YI(TI_YI_loop,2), & 
                      TI_YI(TI_YI_loop,3), & 
                      TI_YI(TI_YI_loop,4)
                    enddo
                   endif

                   if (TSAT_iter.eq.TI_YI_counter) then
                    if (abs(TSAT_correct-TI_YI(TI_YI_best_guess_index,1)).le. &
                        4.0d0*TSAT_ERR) then
                     ! do nothing
                    else
                     print *,"best guess too far from bisection method guess"
                     print *,"TSAT_correct=",TSAT_correct
                     print *,"TI_YI(TI_YI_best_guess_index,1)=", &
                       TI_YI(TI_YI_best_guess_index,1)
                     print *,"TI_YI_best_guess_index=", &
                       TI_YI_best_guess_index
                     do TI_YI_loop=1,TI_YI_counter
                      print *,"TI_YI idx,T,Y,VEL,mdotdiff ", &
                       TI_YI_loop, &
                       TI_YI(TI_YI_loop,1), & 
                       TI_YI(TI_YI_loop,2), & 
                       TI_YI(TI_YI_loop,3), & 
                       TI_YI(TI_YI_loop,4)
                     enddo

                     stop
                    endif
                   else
                    print *,"TSAT_iter.eq.TI_YI_counter failed"
                    print *,"TSAT_iter=",TSAT_iter
                    print *,"TI_YI_counter=",TI_YI_counter
                    stop
                   endif
                  else
                   print *,"TSAT_iter invalid"
                   stop
                  endif

                  Y_predict_hold(ireverse)=Y_predict
                  local_Tsat(ireverse)=TSAT_correct
                  vel_phasechange(ireverse)=VEL_correct

                    ! source
                  temp_target_probe_history(iten+ireverse*num_interfaces,1)= &
                    POUT%T_probe(1)
                  dxprobe_target_history(iten+ireverse*num_interfaces,1)= &
                    POUT%dxprobe_target(1)
                    ! dest
                  temp_target_probe_history(iten+ireverse*num_interfaces,2)= &
                    POUT%T_probe(2)
                  dxprobe_target_history(iten+ireverse*num_interfaces,2)= &
                    POUT%dxprobe_target(2)

                  if (debugrate.eq.1) then
                   print *,"i,j,k,ireverse,vel_phasechange ", &
                    i,j,k,ireverse,vel_phasechange(ireverse)
                  endif

                  ! max velocity cannot exceed dxmin/(2 dt)
                  ! i.e. umax * dt <= dxmin/2
                  if (dt.gt.zero) then
                   if (dxmin.gt.zero) then
                    if (vel_phasechange(ireverse).ge.dxmin/(two*dt)) then
                     vel_phasechange(ireverse)=dxmin/(two*dt)
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

                  if (debugrate.eq.1) then
                   print *,"i,j,k,im_source,im_dest ",i,j,k, &
                    im_source,im_dest 
                   print *,"dt,vel_phasechange(ireverse) ", &
                    dt,vel_phasechange(ireverse)
                   print *,"LL,dxmin ",LL(ireverse),dxmin
                   print *,"dxprobe_target(1)=",POUT%dxprobe_target(1)
                   print *,"dxprobe_target(2)=",POUT%dxprobe_target(2)
                   print *,"thermal_k_model_correct ", &
                    thermal_k_model_correct(1),thermal_k_model_correct(2)
                   print *,"local_Tsat(ireverse) ",local_Tsat(ireverse)
                   print *,"T_Probe(1),T_probe(2) ",POUT%T_Probe(1), &
                    POUT%T_probe(2)
                   print *,"den_I_interp(1) ",POUT%den_I_interp(1)
                   print *,"den_I_interp_SAT(1) ",den_I_interp_SAT(1)
                   print *,"den_I_interp_SAT(2) ",den_I_interp_SAT(2)
                   print *,"LSINTsrc,LSINTdst ",LSINT(im_source),LSINT(im_dest)
                   print *,"nrmCP ",nrmCP(1),nrmCP(2),nrmCP(SDIM)
                   print *,"im_dest= ",im_dest
                  endif
     
                 else if (found_path.eq.0) then
                  ! do nothing
                 else
                  print *,"found_path invalid: ",found_path
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
                print *,"LShere bust3:"
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

            if ((ireverse.eq.0).or. &
                (ireverse.eq.1)) then

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

             LSnew(D_DECL(i,j,k),im_source)= &
               LSnew(D_DECL(i,j,k),im_source)-dt*vel_phasechange(ireverse)
             LSnew(D_DECL(i,j,k),im_dest)= &
               LSnew(D_DECL(i,j,k),im_dest)+dt*vel_phasechange(ireverse)

             do dir=1,SDIM
              if (LShere(im_dest).ge.zero) then
               SIGNVEL=one
               nrmCP(dir)= &
                  LS(D_DECL(i,j,k),num_materials+(im_source-1)*SDIM+dir)
              else if (LShere(im_source).ge.zero) then
               SIGNVEL=-one
               nrmCP(dir)= &
                  LS(D_DECL(i,j,k),num_materials+(im_dest-1)*SDIM+dir)
              else
               print *,"LShere bust4"
               print *,"LShere(im_dest) ",LShere(im_dest)
               print *,"LShere(im_source) ",LShere(im_source)
               stop
              endif
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
             if ((im_dest.ge.1).and.(im_dest.le.num_materials)) then
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
               Tsatfab(D_DECL(i,j,k), &
                 num_interfaces+ncomp_per_tsat*(iten-1)+1)= &
                local_Tsat(ireverse)
               if ((Y_predict_hold(ireverse).ge.zero).and. &
                   (Y_predict_hold(ireverse).le.one)) then
                Tsatfab(D_DECL(i,j,k), &
                  num_interfaces+ncomp_per_tsat*(iten-1)+2)= &
                 Y_predict_hold(ireverse) !default mass fraction=1 (saturated)
               else
                print *,"Y_predict_hold(ireverse) invalid"
                stop
               endif
              else
               print *,"local_Tsat(ireverse) should be positive"
               stop
              endif

              if ((vel_phasechange(ireverse).ge.zero).or. &
                  (vel_phasechange(ireverse).le.zero)) then

               if (1.eq.0) then
                print *,"i,j,k,prev_time,dt,iten,ireverse,vel_phasechange ", &
                 i,j,k,prev_time,dt,iten,ireverse,vel_phasechange(ireverse)
               endif

               do dir=1,ncomp_per_burning

                burnvel(D_DECL(i,j,k), &
                  num_interfaces+(iten-1)*ncomp_per_burning+dir)= &
                 SIGNVEL*nrmCP(dir)*vel_phasechange(ireverse)

                if (1.eq.0) then
                 print *,"i,j,k,prev_time,dt,iten,ireverse,dir,burn ", &
                  i,j,k,prev_time,dt,iten,ireverse,dir, &
                  burnvel(D_DECL(i,j,k), &
                    num_interfaces+(iten-1)*ncomp_per_burning+dir)
                endif

               enddo !dir=1,ncomp_per_burning (1..sdim)

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

           enddo ! im_opp=im+1,num_materials
          enddo ! im=1,num_materials-1

          velmag_sum=zero
          do im=1,num_materials-1
           do im_opp=im+1,num_materials
            call get_iten(im,im_opp,iten)
            burnflag=NINT(burnvel(D_DECL(i,j,k),iten))
            if (burnflag.eq.0) then
             ! do nothing
            else if ((burnflag.eq.1).or. &
                     (burnflag.eq.-1)) then
             local_velmag=zero
             do dir=1,ncomp_per_burning
              local_velmag=local_velmag+ &
               burnvel(D_DECL(i,j,k), &
                 num_interfaces+(iten-1)*ncomp_per_burning+dir)**2
             enddo
             local_velmag=sqrt(local_velmag)
             velmag_sum=velmag_sum+local_velmag
            else
             print *,"burnflag invalid"
             stop
            endif 
           enddo !im_opp=im+1..num_materials
          enddo !im=1..num_materials-1

           ! factor of 2 in order to guarantee that characteristics do not
           ! collide.
           ! LATTICE BOLTZMANN TREATMENT:
           ! max velocity cannot exceed dxmin/(2 dt)
           ! (see above, search for "Li-Shi Luo")
          if ((two*velmag_sum*dt.le.(one+EPS_8_4)*dxmin).and. &
              (velmag_sum.ge.zero)) then
           ! do nothing
          else if (two*velmag_sum*dt.gt.dxmin) then
           print *,"phase change velocity exceeds cfl limits"
           print *,"velmag_sum ",velmag_sum
           print *,"dxmin ",dxmin
           print *,"dt ",dt
           print *,"velmag_sum x dt ",velmag_sum*dt
           print *,"i,j,k ",i,j,k
           print *,"im_primary ",im_primary
           do imls=1,num_materials
            print *,"imls,LShere ",imls,LShere(imls)
           enddo
           print *,"xsten(0,1-3) ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
           do im=1,num_materials-1
            do im_opp=im+1,num_materials
             call get_iten(im,im_opp,iten)
             latent_heat_temperature= &
              EOS(D_DECL(i,j,k), &
                  (im_primary-1)*num_state_material+ENUM_TEMPERATUREVAR+1)
             do ireverse=0,1
              print *,"im,im_opp,ireverse,latent_heat ",im,im_opp,ireverse, &
                get_user_latent_heat(iten+ireverse*num_interfaces, &
                 latent_heat_temperature,0)
              print *,"im,im_opp,ireverse,saturation_temp ", &
                im,im_opp,ireverse, &
                saturation_temp(iten+ireverse*num_interfaces)
             enddo
            enddo
           enddo
           do im=1,num_materials
            print *,"im,T ",im,EOS(D_DECL(i,j,k), &
                    (im-1)*num_state_material+ENUM_TEMPERATUREVAR+1)
           enddo 
           do iten=1,num_interfaces
            do ireverse=0,1
             print *,"iten,ireverse,temp_probe(1) ", &
               iten,ireverse, &
               temp_target_probe_history(iten+ireverse*num_interfaces,1)
             print *,"iten,ireverse,temp_probe(2) ", &
               iten,ireverse, &
               temp_target_probe_history(iten+ireverse*num_interfaces,2)
             print *,"iten,ireverse,dxprobe(1) ", &
               iten,ireverse, &
               dxprobe_target_history(iten+ireverse*num_interfaces,1)
             print *,"iten,ireverse,dxprobe(2) ", &
               iten,ireverse, &
               dxprobe_target_history(iten+ireverse*num_interfaces,2)
            enddo
           enddo
           stop
          else
           print *,"velmag_sum invalid"
           stop
          endif

         else if (is_rigid(im_primary).eq.1) then

          ! do nothing

         else 
          print *,"is_rigid invalid MASS_TRANSFER_3D.F90"
          stop
         endif

        else if (nucleation_flag.eq.1) then
         ! see RatePhaseChange in PROB.F90
         ! LEVELSET FUNCTION AT CELL CENTERS.
         do im=1,num_materials
          LShere(im)=LSnew(D_DECL(i,j,k),im)
         enddo
         call get_primary_material(LShere,im_primary)
         if (is_rigid(im_primary).eq.0) then
          do im=1,num_materials-1
           do im_opp=im+1,num_materials
            if ((im.gt.num_materials).or.(im_opp.gt.num_materials)) then
             print *,"im or im_opp bust 9"
             stop
            endif
            call get_iten(im,im_opp,iten)
            latent_heat_temperature= &
             EOS(D_DECL(i,j,k), &
                 (im_primary-1)*num_state_material+ENUM_TEMPERATUREVAR+1)
            do ireverse=0,1
             LL(ireverse)=get_user_latent_heat(iten+ireverse*num_interfaces, &
                     latent_heat_temperature,0)
             local_freezing_model=freezing_model(iten+ireverse*num_interfaces)
             local_Tsat(ireverse)=saturation_temp(iten+ireverse*num_interfaces)

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

              if ((is_rigid(im).eq.1).or. &
                  (is_rigid(im_opp).eq.1)) then

               ! do nothing

              else if ((is_rigid(im).eq.0).and. &
                       (is_rigid(im_opp).eq.0)) then
               if (im_primary.eq.im_source) then
                create_in%LL=LL(ireverse)
                create_in%local_freezing_model=local_freezing_model
                create_in%local_TSAT=local_Tsat(ireverse)
                create_in%im_source=im_source
                create_in%im_dest=im_dest
                 ! get_vel_phasechange_NUCLEATE is declared in: PROB.F90
                call get_vel_phasechange_NUCLEATE( &
                 create_in,create_inout)
               else if ((im_primary.ge.1).and. &
                        (im_primary.le.num_materials)) then
                ! do nothing
               else
                print *,"im_primary invalid"
                stop
               endif
              else
               print *,"is_rigid(im or im_opp) invalid"
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

         else if (is_rigid(im_primary).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(im_primary) invalid"
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
      end subroutine fort_ratemasschange

      end module mass_transfer_cpp_module

