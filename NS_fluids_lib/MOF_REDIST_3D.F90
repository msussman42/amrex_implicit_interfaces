#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 0

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "EXTRAP_COMP.H"

#include "MOF_REDIST_F.H"

#define nsum 64
#define nsum2 32

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

        module mof_redist_module
        use amrex_fort_module, only : amrex_real
        use probcommon_module

        contains

      subroutine update_closest( &
        xsten_accept,xsten_donate,nhalf, &
        dx,xlo,bfact,level,fablo, &
        mofdata, &
        LSslope_center, &
        imslope_center, &
        nstar, &
        idon,jdon,kdon, & ! donate index
        i1,j1,k1, &  ! accept index: idon+i1,jdon+j1,kdon+k1
        newLS, &
        touch_hold, &
        minLS, &
        maxLS, &
        donateflag, & !1..num_materials+1+nstar
        time)
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: nhalf
      integer, INTENT(in) :: bfact,nstar
      integer, INTENT(in) :: idon,jdon,kdon
      integer, INTENT(in) :: i1,j1,k1
      integer, INTENT(in) :: fablo(SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: xsten_accept(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: xsten_donate(-nhalf:nhalf,SDIM)
      real(amrex_real) :: xsten_vert(-nhalf:nhalf,SDIM)
      real(amrex_real) :: xdonate_vert(SDIM)
      real(amrex_real) :: xdonate_point(SDIM)
      real(amrex_real) :: xaccept_point(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: mofdata(num_materials*ngeom_recon)
      real(amrex_real), INTENT(inout) :: newLS(num_materials*(1+SDIM))
      integer, INTENT(inout) :: touch_hold(num_materials)
      real(amrex_real), INTENT(inout) :: minLS(num_materials)
      real(amrex_real), INTENT(inout) :: maxLS(num_materials)
      integer, INTENT(in) :: donateflag(num_materials+1+nstar)
      integer :: center_stencil
      integer :: im0_center
      integer, INTENT(in) :: imslope_center
      integer :: imslope

      integer nstar_test
      integer istar_array(3)
      integer istar,i2,j2,k2,donateIND
      integer klosten,khisten
      integer dir
      real(amrex_real), INTENT(in) :: LSslope_center(SDIM)
      real(amrex_real) LSslope(SDIM)
      integer n_im
      integer im_test(6)
      integer i_DEB_DIST
      integer j_DEB_DIST
      integer k_DEB_DIST

      if (1.eq.0) then    
       i_DEB_DIST=45
       j_DEB_DIST=68
       k_DEB_DIST=0
      else
       i_DEB_DIST=-9999
       j_DEB_DIST=-9999
       k_DEB_DIST=-9999
      endif

      if (nhalf.lt.3) then
       print *,"nhalf invalid update closest"
       stop
      endif
      do dir=1,SDIM
       LSslope(dir)=zero
      enddo
      imslope=0

      if (bfact.lt.1) then
       print *,"bfact140 invalid"
       stop
      endif
      nstar_test=9
      if (SDIM.eq.3) then
       nstar_test=nstar_test*3
      endif
      if (nstar_test.ne.nstar) then
       print *,"nstar invalid updated nstar nstar_test ",nstar,nstar_test
       stop
      endif

      center_stencil=0
      if ((i1.eq.0).and.(j1.eq.0).and.(k1.eq.0)) then

       center_stencil=1
        ! get the slope of the fluid material whose interface is closest to
        ! the center of the cell.  Slope comes from mofdata and points
        ! towards material "imslope"
        ! get_primary_slope is declared in: MOF.F90
       imslope=imslope_center
       do dir=1,SDIM
        LSslope(dir)=LSslope_center(dir)
       enddo

      else if ((i1.ne.0).or.(j1.ne.0).or.(k1.ne.0)) then
       ! do nothing
      else
       print *,"i1,j1, or k1 bust"
       stop
      endif

      do dir=1,SDIM
       xaccept_point(dir)=xsten_accept(0,dir)
       xdonate_point(dir)=xsten_donate(0,dir)
      enddo

      do dir=1,3
       istar_array(dir)=0
      enddo
      call put_istar(istar,istar_array) 
      im0_center=donateflag(num_materials+1+istar)

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

      call gridsten_level(xsten_vert,idon,jdon,kdon,level,nhalf)
      do k2=klosten,khisten
      do j2=-1,1
      do i2=-1,1
        dir=1
        xdonate_vert(dir)=xsten_vert(i2,dir)
        dir=2
        xdonate_vert(dir)=xsten_vert(j2,dir)
        if (SDIM.eq.3) then
         dir=SDIM
         xdonate_vert(dir)=xsten_vert(k2,dir)
        endif

        istar_array(1)=i2
        istar_array(2)=j2
        istar_array(3)=k2
        call put_istar(istar,istar_array) 
        donateIND=donateflag(num_materials+1+istar)

        if (donateIND.eq.0) then
         ! do nothing (corner is on FAB boundary or corner is
         ! occupied by flotsam that should be ignored)
        else if ((donateIND.ge.1).and. &
                 (donateIND.le.num_materials)) then

         ! this routine will not update solid components.
         ! LSslope only used if xaccept=xdonate=xcell
         ! slope points into the imsource (phi(imsource)>0) region.
         ! phi= n dot (x-xcell) + intercept
         n_im=1
         im_test(1)=donateIND
         if (center_stencil.eq.1) then
          n_im=n_im+1
          im_test(n_im)=im0_center
         endif
          ! compare_distance is declared in: MOF.F90
         call compare_distance( &
          bfact,dx,xsten_donate,nhalf, &
          xaccept_point, &
          xdonate_vert, &
          newLS, & !intent(inout)
          touch_hold, &
          minLS, & !intent(inout)
          maxLS, & !intent(inout) 
          im_test,n_im, &
          LSslope, &  ! slope if xaccept=xdonate=xcell intent(in)
          imslope, &
          im0_center, &
          SDIM, &
          center_stencil, &
          donateflag)

        else
         print *,"donateIND invalid donateIND=",donateIND
         print *,"istar= ",istar
         print *,"idon,jdon,kdon ",idon,jdon,kdon
         print *,"i1,j1,k1 ",i1,j1,k1
         print *,"i2,j2,k2 ",i2,j2,k2
         print *,"center_stencil=",center_stencil
         stop
        endif

      enddo
      enddo
      enddo ! i2,j2,k2

       ! only uses donateflag(1..num_materials+1)
       ! multi_get_distance is declared in: MOF.F90
      call multi_get_distance( &
        bfact,dx,xsten_donate,nhalf,xaccept_point, &
        mofdata, &
        newLS, &
        touch_hold, &
        minLS, &
        maxLS, &
        SDIM, &
        center_stencil, &
        donateflag)

      return
      end subroutine update_closest

      end module mof_redist_module

      module mof_redist_cpp_module
      contains

       ! prior to calling this routine, copy LS_new normal information
       ! to LS_NRM_FD.
       ! called from: NavierStokes::build_NRM_FD_MF (NavierStokes.cpp)
       ! The output from this routine is used by the GNBC algorithm
      subroutine fort_fd_normal( &
       level, &
       finest_level, &
       LS_new, &
       DIMS(LS_new), &
       LS_NRM_FD, &
       DIMS(LS_NRM_FD), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx) &
      bind(c,name='fort_fd_normal')

      use global_utility_module
      use probcommon_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: DIMDEC(LS_new)
      integer, INTENT(in) :: DIMDEC(LS_NRM_FD)

      real(amrex_real), INTENT(in), target :: &
           LS_new(DIMV(LS_new),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_new_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(out), target :: &
           LS_NRM_FD(DIMV(LS_NRM_FD),num_materials*SDIM)
      real(amrex_real), pointer :: LS_NRM_FD_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      integer nhalf
      real(amrex_real) xsten(-3:3,SDIM)
      real(amrex_real) :: centroid_absolute(SDIM)

      integer dir
      integer im

      integer, parameter :: num_particles=0
      real(amrex_real) :: particle_list(1,SDIM+1)

      real(amrex_real) ls_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)
      real(amrex_real) lsnormal(num_materials,SDIM)
      integer lsnormal_valid(num_materials)
      real(amrex_real) ls_intercept(num_materials)
      real(amrex_real) dxmaxLS
      integer k1lo,k1hi
      integer i,j,k
      integer i1,j1,k1
      integer dcomp
      real(amrex_real) local_LS(num_materials)
      integer im_primary,im_secondary,triple_point_flag
      integer, PARAMETER :: continuous_mof=STANDARD_MOF

      nhalf=3 

      LS_NRM_FD_ptr=>LS_NRM_FD
      LS_new_ptr=>LS_new

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid141"
       stop
      endif
      if ((level.le.finest_level).and.(level.ge.0)) then
       ! do nothing
      else
       print *,"level invalid in fort_fd_normal"
       stop
      endif

      if (ngrow_make_distance.eq.3) then
       ! do nothing
      else
       print *,"ngrow_make_distance invalid fort_fd_normal"
       stop
      endif

      call checkbound_array(fablo,fabhi,LS_new_ptr, &
              ngrow_make_distance+1,-1)
      call checkbound_array(fablo,fabhi,LS_NRM_FD_ptr, &
              ngrow_make_distance+1,-1)
      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid fort_fd_normal"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.eq.SDIM+1) then
       ! do nothing
      else
       print *,"ngeom_raw invalid fort_fd_normal"
       print *,"ngeom_raw=",ngeom_raw
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
       print *,"levelrz invalid in fort_fd_normal"
       stop
      endif
      k1lo=0
      k1hi=0
      if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif
      call get_dxmaxLS(dx,bfact,dxmaxLS)

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       do k1=k1lo,k1hi
       do j1=-1,1
       do i1=-1,1
       do im=1,num_materials
        ls_stencil(D_DECL(i1,j1,k1),im)= &
                LS_new(D_DECL(i+i1,j+j1,k+k1),im)
       enddo
       enddo
       enddo
       enddo

       do im=1,num_materials
        local_LS(im)=ls_stencil(D_DECL(0,0,0),im)
       enddo
       call get_primary_material(local_LS,im_primary)
       call get_secondary_material(local_LS,im_primary,im_secondary)
       triple_point_flag=0
       do im=1,num_materials
        if ((im.ne.im_primary).and.(im.ne.im_secondary)) then
         if (abs(local_LS(im)).le.dxmaxLS) then
          triple_point_flag=1
         else if (abs(local_LS(im)).ge.dxmaxLS) then
          ! do nothing
         else
          print *,"local_LS bust"
          stop
         endif
        else if ((im.ge.1).and.(im.le.num_materials)) then
         ! do nothing
        else
         print *,"im bust"
         stop
        endif
       enddo ! im=1..num_materials
      
       if (triple_point_flag.eq.0) then 

        do dir=1,SDIM
         centroid_absolute(dir)=xsten(0,dir)
        enddo

        do im=1,num_materials
         if (is_rigid(im).eq.0) then
          if (abs(local_LS(im)).le.two*dxmaxLS) then
           if ((im.eq.im_primary).or.(im.eq.im_secondary)) then
            call find_cut_geom_slope_CLSVOF( &
             continuous_mof, & !STANDARD_MOF
             ls_stencil, & ! (-1,1)^3,num_materials
             particle_list, &
             num_particles, &
             lsnormal, &  ! (num_materials,sdim)
             lsnormal_valid, &  ! num_materials
             ls_intercept, & ! num_materials
             bfact,dx, &
             xsten,nhalf, &
             centroid_absolute, &
             im, &
             dxmaxLS, &
             SDIM)

            if (lsnormal_valid(im).eq.1) then
             do dir=1,SDIM
              dcomp=(im-1)*SDIM+dir
              LS_NRM_FD(D_DECL(i,j,k),dcomp)=lsnormal(im,dir)
             enddo
            else if (lsnormal_valid(im).eq.0) then
             ! do nothing
            else
             print *,"lsnormal_valid invalid"
             stop
            endif
           else if ((im.ge.1).and.(im.le.num_materials)) then
            ! do nothing
           else
            print *,"im invalid 110"
            stop
           endif

          else if (abs(local_LS(im)).ge.two*dxmaxLS) then
           ! do nothing
          else
           print *,"local_LS invalid"
           stop
          endif
         else if (is_rigid(im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid MOF_REDIST_3D.F90"
          stop
         endif
        enddo ! im=1..num_materials
       else if (triple_point_flag.eq.1) then
        ! do nothing
       else
        print *,"triple_point_flag invalid"
        stop
       endif

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine fort_fd_normal

      subroutine fort_fd_node_normal( &
       level, &
       finest_level, &
       LS_new, &
       DIMS(LS_new), &
       FD_NRM_ND_fab, &
       DIMS(FD_NRM_ND_fab), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       n_normal, &
       ngrow_make_distance_in) &
      bind(c,name='fort_fd_node_normal')

      use global_utility_module
      use probcommon_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: n_normal
      integer, INTENT(in) :: ngrow_make_distance_in
      integer, INTENT(in) :: DIMDEC(LS_new)
      integer, INTENT(in) :: DIMDEC(FD_NRM_ND_fab)
      real(amrex_real), INTENT(in), target :: &
              LS_new(DIMV(LS_new),num_materials*(1+SDIM))
      real(amrex_real), INTENT(out), target :: &
              FD_NRM_ND_fab(DIMV(FD_NRM_ND_fab),n_normal)
      real(amrex_real), pointer :: FD_NRM_ND_fab_ptr(D_DECL(:,:,:),:)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      real(amrex_real) xsten_nd(-3:3,SDIM)
      integer n_normal_test
      integer nhalf
      integer k1hi
      integer i,j,k
      integer i1,j1,k1
      integer im,im1,im2
      integer iten
      integer ibase
      integer dir
      real(amrex_real) local_normal(SDIM)
      real(amrex_real) local_LS
      real(amrex_real) local_mag
      real(amrex_real) sign_nm
      real(amrex_real) xplus,xminus,RR

      nhalf=3 

      FD_NRM_ND_fab_ptr=>FD_NRM_ND_fab

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid141"
       stop
      endif
      if ((level.le.finest_level).and.(level.ge.0)) then
       ! do nothing
      else
       print *,"level invalid in fort_fd_node_normal"
       stop
      endif

      if (ngrow_make_distance.ne.3) then
       print *,"ngrow_make_distance.ne.3 fort_fd_node_normal"
       stop
      endif
      if (ngrow_make_distance_in.ne.3) then
       print *,"ngrow_make_distance_in.ne.3 fort_fd_node_normal"
       stop
      endif
      if (ngrow_distance.ne.4) then
       print *,"ngrow_distance.ne.4 fort_fd_node_normal"
       stop
      endif
      n_normal_test=(SDIM+1)*(num_interfaces+num_materials)
      if (n_normal_test.ne.n_normal) then
       print *,"n_normal invalid fd_node_normal n_normal ",n_normal
       stop
      endif
 
      call checkbound_array(fablo,fabhi,LS_new,ngrow_distance,-1)
      call checkbound_array(fablo,fabhi,FD_NRM_ND_fab_ptr, &
              ngrow_distance,-1)

      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid fort_fd_node_normal"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.eq.SDIM+1) then
       ! do nothing
      else
       print *,"ngeom_raw invalid fort_fd_node_normal"
       print *,"ngeom_raw=",ngeom_raw
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
       print *,"levelrz invalid in fort_fd_node_normal"
       stop
      endif

      k1hi=0
      if (SDIM.eq.3) then
       k1hi=1
      endif

      call growntileboxNODE(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_make_distance)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridstenND_level(xsten_nd,i,j,k,level,nhalf)

       do im=1,num_materials+num_interfaces
        do dir=1,SDIM
         local_normal(dir)=zero
        enddo
        local_mag=zero

        do dir=1,SDIM
         do k1=0,k1hi
         do j1=0,1
         do i1=0,1
          if ((im.ge.1).and.(im.le.num_materials)) then
           local_LS=LS_new(D_DECL(i+i1-1,j+j1-1,k+k1-1),im) 
          else if ((im.ge.num_materials+1).and. &
                   (im.le.num_materials+num_interfaces)) then
           iten=im-num_materials
           call get_inverse_iten(im1,im2,iten)
           local_LS=half*(LS_new(D_DECL(i+i1-1,j+j1-1,k+k1-1),im1)- &
                LS_new(D_DECL(i+i1-1,j+j1-1,k+k1-1),im2))
          else
           print *,"im invalid 111"
           stop
          endif
          sign_nm=one
          if (((dir.eq.1).and.(i1.eq.0)).or. &
              ((dir.eq.2).and.(j1.eq.0)).or. &
              ((dir.eq.3).and.(SDIM.eq.3).and.(k1.eq.0))) then
           sign_nm=-one
          endif
          local_normal(dir)=local_normal(dir)+sign_nm*local_LS
         enddo !  k1
         enddo !  j1
         enddo !  i1
         xplus=xsten_nd(1,dir)
         xminus=xsten_nd(-1,dir)
         if (xplus.gt.xminus) then
          local_normal(dir)=local_normal(dir)/(xplus-xminus)
         else
          print *,"xplus or xminus invalid"
          stop
         endif
         if (levelrz.eq.COORDSYS_CARTESIAN) then
          RR=one
         else if (levelrz.eq.COORDSYS_RZ) then
          RR=one
          if ((dir.eq.1).and.(xminus.lt.zero)) then
           local_normal(dir)=zero
          endif
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
          RR=one
          if ((dir.eq.1).and.(xminus.lt.zero)) then
           local_normal(dir)=zero
          endif
          if (dir.eq.2) then ! theta direction
           RR=xsten_nd(0,1)
           if (RR.gt.zero) then
            ! do nothing
           else
            print *,"RR invalid"
            stop
           endif
          endif 
         else
          print *,"levelrz invalid"
          stop
         endif
         local_normal(dir)=local_normal(dir)/RR
         local_mag=local_mag+local_normal(dir)**2
        enddo ! dir=1..sdim 
        if (local_mag.eq.zero) then
         ! do nothing
        else if (local_mag.gt.zero) then
         local_mag=sqrt(local_mag)
         do dir=1,SDIM
          local_normal(dir)=local_normal(dir)/local_mag
         enddo
        else
         print *,"local_mag invalid"
         stop
        endif
        ibase=(im-1)*(SDIM+1)
        do dir=1,SDIM
         FD_NRM_ND_fab(D_DECL(i,j,k),ibase+dir)=local_normal(dir)
        enddo 
        FD_NRM_ND_fab(D_DECL(i,j,k),ibase+SDIM+1)=local_mag
       enddo ! im=1..num_materials+num_interfaces

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine fort_fd_node_normal

      subroutine simple_htfunc_sum( &
        stenlo,stenhi,vofsten, &
        localsum,sign_change_dir, &
        num_sign_changes_plus, &
        num_sign_changes_minus)

      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: stenlo(3)
      integer, INTENT(in) :: stenhi(3)
      real(amrex_real), INTENT(in), &
              dimension(-ngrow_distance:ngrow_distance, &
                -ngrow_distance:ngrow_distance, &
                -ngrow_distance:ngrow_distance) :: vofsten
      real(amrex_real), INTENT(out) :: localsum
      integer, INTENT(in) :: sign_change_dir
      integer, INTENT(out) :: num_sign_changes_plus
      integer, INTENT(out) :: num_sign_changes_minus
      integer :: i1,j1,k1
      integer :: dir
      integer :: ibase(3),itop(3)
      real(amrex_real) LS_base,LS_top

      num_sign_changes_plus=0
      num_sign_changes_minus=0
      localsum=zero
      do k1=stenlo(3),stenhi(3)
      do j1=stenlo(2),stenhi(2)
      do i1=stenlo(1),stenhi(1)
       ibase(1)=i1
       ibase(2)=j1
       ibase(3)=k1
       localsum=localsum+vofsten(i1,j1,k1)
       if (sign_change_dir.eq.0) then
        ! do nothing
       else if ((sign_change_dir.ge.1).and. &
                (sign_change_dir.le.SDIM)) then 
        do dir=1,3
         itop(dir)=ibase(dir)
        enddo
        itop(sign_change_dir)=itop(sign_change_dir)+1
        if (itop(sign_change_dir).le.ngrow_distance) then
         LS_base=vofsten(i1,j1,k1)-half
         LS_top=vofsten(itop(1),itop(2),itop(3))-half
         if ((LS_base.lt.zero).and. &
             (LS_top.ge.zero).and. &
             (LS_base.lt.LS_top)) then
          num_sign_changes_plus=num_sign_changes_plus+1
         else if ((LS_base.ge.zero).and. &
                  (LS_top.lt.zero).and. &
                  (LS_base.gt.LS_top)) then
          num_sign_changes_minus=num_sign_changes_minus+1
         else if ((LS_base.eq.zero).and. &
                  (LS_top.eq.zero)) then
          num_sign_changes_plus=num_sign_changes_plus+1
          num_sign_changes_minus=num_sign_changes_minus+1
         else if ((LS_base.ge.zero).and. &
                  (LS_top.ge.zero)) then
          ! do nothing
         else if ((LS_base.lt.zero).and. &
                  (LS_top.lt.zero)) then
          ! do nothing
         else
          print *,"LS_base or LS_top invalid"
          stop
         endif
        else if (itop(sign_change_dir).eq.ngrow_distance+1) then
         ! do nothing
        else
         print *,"itop(sign_change_dir) invalid"
         stop
        endif
       else
        print *,"sign_change_dir invalid"
        stop
       endif
      enddo !k1
      enddo !j1
      enddo !i1

      end subroutine simple_htfunc_sum


       ! since the FillBoundary command is issued after this command,
       ! "localMF[FD_CURV_CELL_MF]->FillBoundary(geom.periodicity());",
       ! we zero out "CURV_CELL" in the ghost cells.  i.e. the status
       ! should be 0 ("invalid") at coarse/fine borders and domain
       ! borders.  At fine-fine borders, the status should be corrected
       ! after the FillBoundary command.
       ! "fort_node_to_cell" is called after "fort_fd_node_normal"
      subroutine fort_node_to_cell( &
       tid_current, &
       level, &
       finest_level, &
       height_function_flag, &  ! 1=> use height function 0 => use FD
       F_new, &  !F_new(i,j,k,im)  im=1..num_materials*ngeom_recon
       DIMS(F_new), &
       LS_new, &
       DIMS(LS_new), &
       FD_NRM_ND_fab, &
       DIMS(FD_NRM_ND_fab), &
       CURV_CELL, &
       DIMS(CURV_CELL), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       n_normal, &
       ngrow_make_distance_in) &
      bind(c,name='fort_node_to_cell')

      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_current
      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: height_function_flag
      integer, INTENT(in) :: n_normal
      integer, INTENT(in) :: ngrow_make_distance_in
      integer, INTENT(in) :: DIMDEC(LS_new)
      integer, INTENT(in) :: DIMDEC(F_new)
      integer, INTENT(in) :: DIMDEC(FD_NRM_ND_fab)
      integer, INTENT(in) :: DIMDEC(CURV_CELL)
      real(amrex_real), INTENT(in), target :: &
        F_new(DIMV(F_new),num_materials*ngeom_recon)
      real(amrex_real), pointer :: F_new_ptr(D_DECL(:,:,:),:)

      real(amrex_real), allocatable, target :: F_tess(D_DECL(:,:,:),:)
      real(amrex_real), pointer :: F_tess_ptr(D_DECL(:,:,:),:)
      integer DIMDEC(F_tess)

      real(amrex_real), INTENT(in), target :: &
         LS_new(DIMV(LS_new),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_new_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              FD_NRM_ND_fab(DIMV(FD_NRM_ND_fab),n_normal)
      real(amrex_real), pointer :: FD_NRM_ND_fab_ptr(D_DECL(:,:,:),:)
       ! first num_materials+num_interfaces components are curvature
       ! second num_materials+num_interfaces components are 
       !   status (0=bad 1=good)
      real(amrex_real), INTENT(out), target ::  &
              CURV_CELL(DIMV(CURV_CELL),2*(num_materials+num_interfaces))
      real(amrex_real), pointer :: CURV_CELL_ptr(D_DECL(:,:,:),:)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      real(amrex_real) xsten(-3:3,SDIM)
      real(amrex_real) xsten_grow(-(2*ngrow_distance+1):(2*ngrow_distance+1),SDIM)
      real(amrex_real) xcenter(SDIM)

      real(amrex_real) dx_col(SDIM)
      real(amrex_real) x_col(SDIM)
      real(amrex_real) x_col_avg(SDIM)

      integer nhalf_grow
      integer n_normal_test
      integer nhalf
      integer k1hi
      integer i,j,k
      integer i1,j1,k1
      integer ii,jj,kk
      integer im
      integer im_local
      integer im1,im2
      integer im_primary,im_secondary
      integer iten
      integer ibase
      integer dir
      integer dirloc
      integer dir2
      integer side
      real(amrex_real) local_curv(SDIM)
      real(amrex_real) local_normal(SDIM)
      real(amrex_real) local_mag
      real(amrex_real) sign_nm
      real(amrex_real) xplus,xminus,xmiddle,RR
      real(amrex_real) denom_factor
      real(amrex_real) total_curv
      integer local_status
      integer crossing_status
      real(amrex_real) local_LS(num_materials)
      real(amrex_real), dimension(-ngrow_distance:ngrow_distance, &
                -ngrow_distance:ngrow_distance, &
                -ngrow_distance:ngrow_distance) :: vofsten
      real(amrex_real), dimension(-ngrow_distance:ngrow_distance, &
                -ngrow_distance:ngrow_distance, &
                -ngrow_distance:ngrow_distance) :: lssten
      integer ngrow_null
      real(amrex_real) curv_LS
      real(amrex_real) curv_VOF
      real(amrex_real) curv_choice
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      integer tessellate
      integer nmax
      integer LSstenlo(3),LSstenhi(3)
      integer HTstenlo(3),HTstenhi(3)
      real(amrex_real) vof_local(num_materials)
      real(amrex_real) slopesum(SDIM,-1:1)
      real(amrex_real) localsum
      integer sign_change
      integer num_sign_changes
      integer num_sign_changes_plus
      integer num_sign_changes_minus
      integer total_num_sign_changes_plus
      integer total_num_sign_changes_minus
      integer normal_dir
      real(amrex_real) normal_max
      real(amrex_real) normal_test(SDIM)
      integer sign_change_dir
      integer vofcomp
      real(amrex_real) htfunc_LS(-1:1,-1:1)
      real(amrex_real) htfunc_VOF(-1:1,-1:1)
      real(amrex_real) col_ht_LS
      real(amrex_real) col_ht_VOF
      integer icol,jcol,kcol
      integer ivert,itan,jtan
      integer iwidth
      integer iwidthnew
      integer jwidth
      real(amrex_real) n1d
      real(amrex_real) ls_column(-ngrow_distance:ngrow_distance)
      real(amrex_real) vof_column(-ngrow_distance:ngrow_distance)

      if ((tid_current.lt.0).or.(tid_current.ge.geom_nthreads)) then
       print *,"tid_current invalid"
       stop
      endif

      ngrow_null=0
      nhalf=3 
      nhalf_grow=2*ngrow_distance+1
 
      nmax=POLYGON_LIST_MAX ! in: fort_node_to_cell

      CURV_CELL_ptr=>CURV_CELL

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid141"
       stop
      endif
      if ((level.le.finest_level).and.(level.ge.0)) then
       ! do nothing
      else
       print *,"level invalid in fort_node_to_cell"
       stop
      endif

      if (ngrow_make_distance.ne.3) then
       print *,"ngrow_make_distance.ne.3 fort_node_to_cell"
       stop
      endif
      if (ngrow_make_distance_in.ne.3) then
       print *,"ngrow_make_distance_in.ne.3 fort_node_to_cell"
       stop
      endif
      if (ngrow_distance.ne.4) then
       print *,"ngrow_distance invalid"
       stop
      endif

      n_normal_test=(SDIM+1)*(num_interfaces+num_materials)
      if (n_normal_test.ne.n_normal) then
       print *,"n_normal invalid node_to_cell n_normal ",n_normal
       stop
      endif

      F_new_ptr=>F_new 
      call checkbound_array(fablo,fabhi,F_new_ptr,ngrow_distance,-1)
      LS_new_ptr=>LS_new 
      call checkbound_array(fablo,fabhi,LS_new_ptr,ngrow_distance,-1)
      FD_NRM_ND_fab_ptr=>FD_NRM_ND_fab
      call checkbound_array(fablo,fabhi,FD_NRM_ND_fab_ptr, &
              ngrow_distance,-1)
      call checkbound_array(fablo,fabhi,CURV_CELL_ptr, &
              ngrow_make_distance,-1)

      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid fort_node_to_cell"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.eq.SDIM+1) then
       ! do nothing
      else
       print *,"ngeom_raw invalid fort_node_to_cell"
       print *,"ngeom_raw=",ngeom_raw
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
       print *,"levelrz invalid in fort_node_to_cell"
       stop
      endif

      k1hi=0
      if (SDIM.eq.2) then
       k1hi=0
       denom_factor=two
      else if (SDIM.eq.3) then
       k1hi=1
       denom_factor=four
      else
       print *,"dimension bust"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_distance)

       ! input: growlo,growhi
       ! output: DIMS(F_tess)
      call box_to_dim( &
        DIMS(F_tess), &
        growlo,growhi)
      allocate(F_tess(DIMV(F_tess),num_materials*ngeom_recon))
      
      F_tess_ptr=>F_tess
      call checkbound_array(tilelo,tilehi,F_tess_ptr,ngrow_distance,-1)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       do im=1,num_materials*ngeom_recon
        mofdata(im)=F_new(D_DECL(i,j,k),im)
       enddo

       ! before (mofdata): fluids tessellate, solids are embedded.
       ! after  (mofdata): fluids and solids tessellate
       ! if tessellate==3:
       !  if solid material(s) dominate the cell, then F_solid_raster=1
       !  and F_fluid=0.
       !  if fluid material(s) dominate the cell, then F_solid=0,
       !  sum F_fluid=1
       tessellate=3

        !EPS2
       call multi_get_volume_tessellate( &
         tessellate, & ! =1 or 3
         bfact, &
         dx, &
         xsten,nhalf, &
         mofdata, &
         geom_xtetlist(1,1,1,tid_current+1), &
         nmax, &
         nmax, &
         SDIM)

       do im=1,num_materials*ngeom_recon
        F_tess(D_DECL(i,j,k),im)=mofdata(im)
       enddo

      enddo
      enddo
      enddo

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_make_distance)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       ! first num_materials+num_interfaces components are curvature
       ! second num_materials+num_interfaces components are 
       !   status (0=bad 1=good)
       do im=1,2*(num_materials+num_interfaces)
        CURV_CELL(D_DECL(i,j,k),im)=zero
       enddo
      enddo
      enddo
      enddo

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_null)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       ! first num_materials+num_interfaces components are curvature
       ! second num_materials+num_interfaces components are 
       !   status (0=bad 1=good)

       call gridsten_level(xsten,i,j,k,level,nhalf)
       call gridsten_level(xsten_grow,i,j,k,level,nhalf_grow)
       do dir=1,SDIM
        xcenter(dir)=xsten(0,dir)
       enddo

        ! fort_ratemasschange will not consider cells in which
        ! (is_rigid(im_primary).eq.0) 
        !
       do im=1,num_materials+num_interfaces

        do im_local=1,num_materials
         local_LS(im_local)=LS_new(D_DECL(i,j,k),im_local)
        enddo
        call get_primary_material(local_LS,im_primary)
        call get_secondary_material(local_LS,im_primary,im_secondary)

        local_status=1

        if (is_rigid(im_primary).eq.1) then

         local_status=0

        else if (is_rigid(im_primary).eq.0) then

         if ((im.ge.1).and.(im.le.num_materials)) then

          if (is_rigid(im).eq.1) then

           local_status=0

          else if (is_rigid(im).eq.0) then

           if ((im.ne.im_primary).and. &
               (im.ne.im_secondary)) then
            local_status=0
           endif

          else
           print *,"is_rigid(im) invalid"
           stop
          endif

         else if ((im.ge.num_materials+1).and. &
                  (im.le.num_materials+num_interfaces)) then

          iten=im-num_materials
          call get_inverse_iten(im1,im2,iten)
          if (is_rigid(im1).eq.1) then
           local_status=0
          else if (is_rigid(im1).eq.0) then
           if (is_rigid(im2).eq.1) then
            local_status=0
           else if (is_rigid(im2).eq.0) then
           
            if ((im1.ne.im_primary).and. &
                (im1.ne.im_secondary)) then
             local_status=0
            endif
            if ((im2.ne.im_primary).and. &
                (im2.ne.im_secondary)) then
             local_status=0
            endif

           else
            print *,"is_rigid(im2) invalid"
            stop
           endif

          else
           print *,"is_rigid(im1) invalid"
           stop
          endif

         else
          print *,"im invalid 112"
          stop
         endif

         if (local_status.eq.0) then

          ! do nothing

         else if (local_status.eq.1) then

          do dir=1,SDIM
           local_curv(dir)=zero
          enddo

          ibase=(im-1)*(SDIM+1)

          do dir=1,SDIM
           do k1=0,k1hi
           do j1=0,1
           do i1=0,1
            local_normal(dir)=FD_NRM_ND_fab(D_DECL(i+i1,j+j1,k+k1),ibase+dir)
            local_mag=FD_NRM_ND_fab(D_DECL(i+i1,j+j1,k+k1),ibase+SDIM+1)
            if (local_mag.gt.zero) then
             ! do nothing
            else if (local_mag.eq.zero) then
             local_status=0
            else
             print *,"local_mag invalid"
             stop
            endif

            sign_nm=one
            if (((dir.eq.1).and.(i1.eq.0)).or. &
                ((dir.eq.2).and.(j1.eq.0)).or. &
                ((dir.eq.3).and.(SDIM.eq.3).and.(k1.eq.0))) then
             sign_nm=-one
            endif
            RR=one
            if (levelrz.eq.COORDSYS_CARTESIAN) then
             ! do nothing
            else if (levelrz.eq.COORDSYS_RZ) then
             if (dir.eq.1) then
              RR=xsten(2*i1-1,dir)
              if (RR.lt.zero) then
               RR=zero
              endif
             endif
            else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
             if (dir.eq.1) then
              RR=xsten(2*i1-1,dir)
              if (RR.lt.zero) then
               RR=zero
              endif
             endif
            else
             print *,"dir invalid"
             stop
            endif

            local_curv(dir)=local_curv(dir)+sign_nm*RR*local_normal(dir)
           enddo !  k1
           enddo !  j1
           enddo !  i1
           xplus=xsten(1,dir)
           xminus=xsten(-1,dir)
           xmiddle=xsten(0,dir)
           if ((xplus.gt.xminus).and. &
               (xplus.gt.xmiddle).and. &
               (xmiddle.gt.xminus)) then
            local_curv(dir)=local_curv(dir)/(denom_factor*(xplus-xminus))
           else
            print *,"xplus or xminus invalid"
            stop
           endif
           RR=one
           if (levelrz.eq.COORDSYS_CARTESIAN) then
            RR=one
           else if (levelrz.eq.COORDSYS_RZ) then
            RR=one
            if (dir.eq.1) then
             if (xmiddle.le.zero) then
              local_curv(dir)=zero
             else if (xmiddle.gt.zero) then
              RR=xmiddle
             else
              print *,"xmiddle invalid"
              stop
             endif
            endif
           else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
            RR=one
            if (dir.eq.1) then
             if (xmiddle.le.zero) then
              local_curv(dir)=zero
             else if (xmiddle.gt.zero) then
              RR=xmiddle
             else
              print *,"xmiddle invalid"
              stop
             endif
            endif
            if (dir.eq.2) then ! theta direction
             if (xsten(0,1).le.zero) then
              local_curv(dir)=zero
             else if (xsten(0,1).gt.zero) then
              RR=xsten(0,1)
             else
              print *,"xsten(0,1) invalid"
              stop
             endif
            endif 
           else
            print *,"levelrz invalid"
            stop
           endif
           local_curv(dir)=local_curv(dir)/RR
          enddo ! dir=1..sdim 

          total_curv=zero
          do dir=1,SDIM
           total_curv=total_curv+local_curv(dir)
          enddo
          curv_LS=total_curv
          curv_VOF=total_curv

           ! for 1<=im<=num_materials, use F_m,  L_m=F_m-1/2, 
           !   L_m should change sign in the "cross" stencil.
           !   L_m should be primary or secondary in the center,
           !   and center should not be dominated by an "is_rigid" material.
           ! for num_materials+1<=im<=num_materials+num_interfaces:
           !  iten=im-num_materials
           !  iten => m1,m2
           !  use 
           !(F_m1-F_m2+1)/2,L_m1=(1/2)((F_m1-1/2)-(F_m2-1/2))=(F_m1-F_m2)/2
           !   L_m1 should change sign in the "cross" stencil.
           !   both L_m1 and L_m2 should be primary or secondary in the center,
           ! If a height function curvature can be succesfully computed
           ! (assuming hypothetically that X-Y or X-Y-Z coordinates are used)
           ! then no need to limit the interface temperature for evaporation or
           ! condensation, otherwise the interface temperature (for evap or 
           ! or condensation) must be limited. 
   
          LSstenlo(3)=0
          LSstenhi(3)=0
          do dir=1,SDIM
           LSstenlo(dir)=-ngrow_distance
           LSstenhi(dir)=ngrow_distance
          enddo

          ! i1,j1,k1=-ngrow_distance ... ngrow_distance
          do k1=LSstenlo(3),LSstenhi(3)
          do j1=LSstenlo(2),LSstenhi(2)
          do i1=LSstenlo(1),LSstenhi(1)
           do im_local=1,num_materials
            vofcomp=(im_local-1)*ngeom_recon+1
             ! we do not look at the tessellating volume fractions since
             ! we need the curvature near embedded walls.
            call safe_data(i+i1,j+j1,k+k1,vofcomp, &
              F_new_ptr, &
              vof_local(im_local))
           enddo
           if ((im.ge.1).and.(im.le.num_materials)) then
            vofsten(i1,j1,k1)=vof_local(im)
            lssten(i1,j1,k1)=vof_local(im)-half
           else if ((im.ge.num_materials+1).and. &
                    (im.le.num_materials+num_interfaces)) then
            iten=im-num_materials
            call get_inverse_iten(im1,im2,iten)
            vofsten(i1,j1,k1)= &
               half*(vof_local(im1)-vof_local(im2)+one) 
            lssten(i1,j1,k1)= &
               half*(vof_local(im1)-vof_local(im2))
           else
            print *,"im invalid"
            stop
           endif 
          enddo !k1
          enddo !j1
          enddo !i1

          ! check for sign change in cross stencil
          sign_change=0
          normal_max=zero
          normal_dir=0
          do dir=1,SDIM
           ii=0
           jj=0
           kk=0
           if (dir.eq.1) then
            ii=1
           else if (dir.eq.2) then
            jj=1
           else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
            kk=1
           else
            print *,"dir invalid"
            stop
           endif
           do side=-1,1,2
            if (lssten(0,0,0)* &
                lssten(ii*side,jj*side,kk*side).le.zero) then
                sign_change=sign_change+1
            else if (lssten(0,0,0)* &
                     lssten(ii*side,jj*side,kk*side).gt.zero) then
             !do nothing
            else
             print *,"lssten bust"
             stop
            endif
            HTstenlo(3)=0
            HTstenhi(3)=0
            do dirloc=1,SDIM
             HTstenlo(dirloc)=-1
             HTstenhi(dirloc)=1
            enddo
            HTstenlo(dir)=side
            HTstenhi(dir)=side
            sign_change_dir=0
            call simple_htfunc_sum(HTstenlo,HTstenhi,vofsten, &
              slopesum(dir,side),sign_change_dir, &
              num_sign_changes_plus, &
              num_sign_changes_minus)
           enddo ! side=-1,1,2

           normal_test(dir)=abs(slopesum(dir,1)-slopesum(dir,-1))

           if (dir.eq.1) then
            if (levelrz.eq.COORDSYS_CARTESIAN) then
             ! do nothing
            else if ((levelrz.eq.COORDSYS_RZ).or.(levelrz.eq.COORDSYS_CYLINDRICAL)) then

             if (xsten_grow(-2*ngrow_distance,1).lt.zero) then
              normal_test(dir)=zero
             else if (xsten_grow(-2*ngrow_distance,1).ge.zero) then
              ! do nothing
             else
              print *,"xsten_grow is NaN"
              stop
             endif

            else
             print *,"levelrz invalid"
             stop
            endif
           else if ((dir.eq.2).or.(dir.eq.SDIM)) then
            ! do nothing
           else
            print *,"dir invalid"
            stop
           endif

           if (normal_test(dir).gt.normal_max) then
            normal_max=normal_test(dir)
            normal_dir=dir
           else if (normal_test(dir).le.normal_max) then
            ! do nothing
           else
            print *,"normal_test invalid"
            stop
           endif
          enddo ! dir=1..sdim
          if (sign_change.eq.0) then
           local_status=0
          endif
          if (normal_max.eq.zero) then
           local_status=0
          else if (normal_max.gt.zero) then
           if ((normal_dir.ge.1).and.(normal_dir.le.SDIM)) then
            LSstenlo(3)=0
            LSstenhi(3)=0
            do dir=1,SDIM
             LSstenlo(dir)=-1
             LSstenhi(dir)=1
            enddo
            LSstenlo(normal_dir)=0
            LSstenhi(normal_dir)=0

            total_num_sign_changes_plus=0
            total_num_sign_changes_minus=0

            do k1=LSstenlo(3),LSstenhi(3)
            do j1=LSstenlo(2),LSstenhi(2)
            do i1=LSstenlo(1),LSstenhi(1)
             HTstenlo(1)=i1
             HTstenhi(1)=i1
             HTstenlo(2)=j1
             HTstenhi(2)=j1
             HTstenlo(3)=k1
             HTstenhi(3)=k1
             HTstenlo(normal_dir)=-ngrow_distance
             HTstenhi(normal_dir)=ngrow_distance
           
             sign_change_dir=normal_dir
             call simple_htfunc_sum(HTstenlo,HTstenhi,vofsten, &
               localsum,sign_change_dir, &
               num_sign_changes_plus, &
               num_sign_changes_minus)
             total_num_sign_changes_plus= &
               total_num_sign_changes_plus+num_sign_changes_plus
             total_num_sign_changes_minus= &
               total_num_sign_changes_minus+num_sign_changes_minus

             num_sign_changes= &
                 num_sign_changes_plus+num_sign_changes_minus
             if (num_sign_changes.eq.1) then
              if ((height_function_flag.eq.0).or. &
                  (local_status.eq.0)) then
               ! do nothing (curv_LS, curv_VOF = total_curv above)
              else if ((height_function_flag.eq.1).and. &
                       (local_status.eq.1)) then

               do ivert=-ngrow_distance,ngrow_distance
                icol=i1  
                jcol=j1  
                kcol=k1  

                if (normal_dir.eq.1) then
                 icol=ivert
                 iwidth=jcol
                 jwidth=kcol
                 itan=2
                 jtan=3
                else if (normal_dir.eq.2) then
                 jcol=ivert
                 iwidth=icol
                 jwidth=kcol
                 itan=1
                 jtan=3
                else if ((normal_dir.eq.3).and.(SDIM.eq.3)) then
                 kcol=ivert
                 iwidth=icol
                 jwidth=jcol
                 itan=1
                 jtan=2
                else
                 print *,"normal_dir invalid"
                 stop
                endif
                ls_column(ivert)=lssten(icol,jcol,kcol)
                vof_column(ivert)=vofsten(icol,jcol,kcol)
               enddo !ivert=-ngrow_distance,ngrow_distance

               if (num_sign_changes_plus.eq.1) then
                n1d=one
               else if (num_sign_changes_minus.eq.1) then
                n1d=-one
               else
                print *,"num_sign_changes bust"
                stop
               endif
         
               iwidthnew=iwidth
          
               if (levelrz.eq.COORDSYS_CARTESIAN) then
                ! do nothing
               else if (levelrz.eq.COORDSYS_RZ) then
                if (SDIM.ne.2) then
                 print *,"dimension bust"
                 stop
                endif
                if (itan.eq.1) then ! vertical columns
                 if (iwidth.eq.-1) then
                  if (xsten_grow(2*iwidth,1).le.zero) then
                   iwidthnew=0
                  endif
                 endif
                endif
               else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
                if (xsten_grow(-2,1).le.zero) then
                 print *,"xsten_grow cannot be negative for levelrz==COORDSYS_CYLINDRICAL"
                 stop
                endif
               else
                print *,"levelrz invalid node_to_cell"
                stop
               endif
         
               do dir2=1,SDIM
                x_col(dir2)=xsten_grow(0,dir2)
                x_col_avg(dir2)=half*(xsten_grow(1,dir2)+xsten_grow(-1,dir2))
                dx_col(dir2)=xsten_grow(1,dir2)-xsten_grow(-1,dir2)
               enddo
               x_col(itan)=xsten_grow(2*iwidthnew,itan)
               dx_col(itan)=xsten_grow(2*iwidthnew+1,itan)- &
                            xsten_grow(2*iwidthnew-1,itan)
               x_col_avg(itan)=half*(xsten_grow(2*iwidthnew+1,itan)+ &
                                     xsten_grow(2*iwidthnew-1,itan))

               if ((SDIM.eq.3).or. &
                   ((SDIM.eq.2).and.(jtan.ge.1).and.(jtan.le.SDIM))) then
                x_col(jtan)=xsten_grow(2*jwidth,jtan)
                dx_col(jtan)=xsten_grow(2*jwidth+1,jtan)- &
                             xsten_grow(2*jwidth-1,jtan)
                x_col_avg(jtan)=half*(xsten_grow(2*jwidth+1,jtan)+ &
                                      xsten_grow(2*jwidth-1,jtan))
               else if ((SDIM.eq.2).and.(jtan.eq.3)) then
                !do nothing
               else
                print *,"sdim or jtan invalid"
                stop
               endif

               call get_col_ht_LS( &
                height_function_flag, &
                crossing_status, &  !crossing_status=1=>success
                bfact, &
                dx, &
                xsten_grow, &
                dx_col, &
                x_col, &
                x_col_avg, &
                ls_column, & 
                vof_column, & 
                col_ht_LS, &
                col_ht_VOF, &
                normal_dir, &
                n1d, & !n1d==1.0d0=>im on top,n1d==-1.0d0 => im on bot
                SDIM)

               if (crossing_status.eq.0) then
                local_status=0
               else if (crossing_status.eq.1) then
                htfunc_LS(iwidth,jwidth)=col_ht_LS
                htfunc_VOF(iwidth,jwidth)=col_ht_VOF
               else
                print *,"crossing_status invalid"
                stop
               endif
              else
               print *,"height_function_flag or local_status invalid"
               stop
              endif

             else if ((num_sign_changes.eq.0).or. &
                      (num_sign_changes.gt.1)) then
              local_status=0
             else
              print *,"num_sign_changes invalid"
              stop
             endif

            enddo !k1
            enddo !j1
            enddo !i1

            if ((total_num_sign_changes_plus.gt.0).and. &
                (total_num_sign_changes_minus.gt.0)) then
             local_status=0
            else if ((total_num_sign_changes_plus.eq.0).and. &
                     (total_num_sign_changes_minus.eq.0)) then
             local_status=0
            else if ((total_num_sign_changes_plus.eq.0).and. &
                     (total_num_sign_changes_minus.gt.0)) then
             ! do nothing
            else if ((total_num_sign_changes_minus.eq.0).and. &
                     (total_num_sign_changes_plus.gt.0)) then
             ! do nothing
            else
             print *,"total_num_sign_changes invalid"
             stop
            endif

           else
            print *,"normal_dir invalid"
            stop
           endif
          else
           print *,"normal_max invalid"
           stop
          endif

          if ((height_function_flag.eq.0).or. &
              (local_status.eq.0)) then
           curv_choice=curv_LS
          else if ((height_function_flag.eq.1).and. &
                   (local_status.eq.1)) then
           call analyze_heights( &
            htfunc_LS, &
            htfunc_VOF, &
            xsten_grow, &
            nhalf_grow, &
            itan,jtan, &
            curv_LS, &
            curv_VOF, &
            curv_choice, &
            normal_dir, &
            xcenter, &
            n1d, &
            local_status, &
            height_function_flag)
          else
           print *,"height_function_flag or local_status invalid"
           stop
          endif

         else
          print *,"local_status invalid"
          stop
         endif

        else
         print *,"is_rigid(im_primary) invalid"
         stop
        endif
 
        if (local_status.eq.0) then
         curv_choice=zero
        else if (local_status.eq.1) then
         ! do nothing
        else
         print *,"local_status invalid"
         stop
        endif
 
        CURV_CELL(D_DECL(i,j,k),im)=curv_choice
        CURV_CELL(D_DECL(i,j,k),num_materials+num_interfaces+im)=local_status

       enddo ! im=1..num_materials+num_interfaces

      enddo
      enddo
      enddo  !i,j,k 

      deallocate(F_tess)

      return
      end subroutine fort_node_to_cell

        ! vofrecon=vof,ref centroid,order,slope,intercept
        ! newfab has num_materials*(sdim+1) components
        !
      subroutine fort_levelstrip( &
         nprocessed, &
         minLS, &
         maxLS, &
         max_problen, &
         level, &
         finest_level, &
         maskfab,DIMS(maskfab), &
         stenfab,DIMS(stenfab), &
         vofrecon,DIMS(vofrecon), &
         newfab,DIMS(newfab), &
         touchfab,DIMS(touchfab), &
         crsetouch,DIMS(crsetouch), & !traverse coarsest to finest level.
         crsedist,DIMS(crsedist), & !traverse coarsest to finest level.
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         bc, &
         rz_flag, &
         xlo,dx, &
         time, &
         ngrow_distance_in, &
         nstar) &
      bind(c,name='fort_levelstrip')

      use global_utility_module
      use probcommon_module
      use MOF_routines_module
      use mof_redist_module

      IMPLICIT NONE

      integer, INTENT(inout) :: nprocessed
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nstar
      integer, INTENT(in) :: ngrow_distance_in

      integer, PARAMETER :: ngrow_make_distance_accept=3

      real(amrex_real), INTENT(inout) :: minLS(num_materials)
      real(amrex_real), INTENT(inout) :: maxLS(num_materials)
      real(amrex_real), INTENT(in) :: max_problen
      integer, INTENT(in) :: DIMDEC(maskfab)
      integer, INTENT(in) :: DIMDEC(stenfab)
      integer, INTENT(in) :: DIMDEC(vofrecon)
      integer, INTENT(in) :: DIMDEC(newfab)
      integer, INTENT(in) :: DIMDEC(touchfab)
      integer, INTENT(in) :: DIMDEC(crsetouch)
      integer, INTENT(in) :: DIMDEC(crsedist)

      real(amrex_real), INTENT(in), target :: maskfab(DIMV(maskfab),4)
      real(amrex_real), pointer :: maskfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: stenfab(DIMV(stenfab),nstar)
      real(amrex_real), pointer :: stenfab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: &
       vofrecon(DIMV(vofrecon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: vofrecon_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
        newfab(DIMV(newfab),num_materials*(1+SDIM))
      real(amrex_real), pointer :: newfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              touchfab(DIMV(touchfab),num_materials)
      real(amrex_real), pointer :: touchfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              crsetouch(DIMV(crsetouch),num_materials)
      real(amrex_real), pointer :: crsetouch_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
        crsedist(DIMV(crsedist),num_materials*(1+SDIM))
      real(amrex_real), pointer :: crsedist_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: rz_flag
      integer, INTENT(in) :: bc(SDIM,2)
      real(amrex_real), INTENT(in) :: time

      integer i,j,k,i1,j1,k1
      real(amrex_real) vcenter(num_materials)
      real(amrex_real) VFRAC_TEMP
      integer nhalf
      real(amrex_real) xsten_accept(-3:3,SDIM)
      real(amrex_real) xsten_donate(-3:3,SDIM)
      integer dir,dir2

      integer isten,jsten,ksten
      integer iside,jside,kside
      integer klosten,khisten
      integer fluid_materials_in_cell_stencil
      integer im
      integer vofcomp
      real(amrex_real) newfab_hold(num_materials*(1+SDIM))
      integer touch_hold(num_materials)
      integer ilocut(3),ihicut(3)
      integer icur(3)
      integer istar_array(3)
      integer istar_array_offset(3)
      integer sorted_list(num_materials)
      integer nstar_test
      integer FSI_exclude
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      integer istar
      integer im_crit
      integer stencil_test(num_materials)
      integer cell_test(num_materials)
      integer full_neighbor(num_materials)
      integer i3,j3,k3
      integer i4,j4,k4
      integer i4low(3)
      integer i4high(3)
      integer i4_array(3)
      integer im_corner
      integer im_test_stencil
      integer im_test_center
      real(amrex_real) FSUM(num_materials)
      integer on_border
       ! 1..num_materials,fluid materials in cell,nstar
       ! donateflag(num_materials+2 ... num_materials+1+nstar)=
       ! fluid material id that owns
       ! the respective star stencil position.
       ! donateflag(1..num_materials)=1 if the respective material is a fluid
       ! which has a "non-flotsam" presence in the cell.
       ! If (cell_test(im)==1) and 
       !    (the local star point is owned by material im) then
       !  donateflag(num_materials+1+istar)=im 
      integer donateflag(num_materials+1+nstar)
      integer crse_dist_valid
      integer ctouch
      real(amrex_real) init_dist_from_crse
      integer rigid_in_stencil
      real(amrex_real) VFRAC_STENCIL_CUTOFF
      real(amrex_real) VFRAC_INTERP,theta_nbr,theta_cen
      real(amrex_real) xnbr(SDIM)
      real(amrex_real) xmid(SDIM)
      integer i_DEB_DIST
      integer j_DEB_DIST
      integer k_DEB_DIST

      integer keep_flotsam
      integer legitimate_material
      real(amrex_real) LSslope_center(SDIM)
      integer imslope_center
      real(amrex_real) bypass_cutoff
      real(amrex_real) crude_dist,dotprod,crude_normal
      integer bypass_update_closest

      newfab_ptr=>newfab
      touchfab_ptr=>touchfab

      if (1.eq.0) then    
       i_DEB_DIST=45 
       j_DEB_DIST=68
       k_DEB_DIST=0
      else
       i_DEB_DIST=-9999
       j_DEB_DIST=-9999
       k_DEB_DIST=-9999
      endif

      nhalf=3 
     
      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid142"
       stop
      endif

      if (max_problen.gt.zero) then
       !do nothing
      else
       print *,"max_problen invalid: ",max_problen
       stop
      endif

      nstar_test=9
      if (SDIM.eq.3) then
       nstar_test=nstar_test*3
      endif
      if (nstar_test.ne.nstar) then
       print *,"nstar invalid levelstrip nstar nstar_test ",nstar,nstar_test
       stop
      endif

      if ((level.le.finest_level).and.(level.ge.0)) then
       ! do nothing
      else
       print *,"level invalid in levelstrip"
       stop
      endif

      if (ngrow_distance.ne.4) then
       print *,"ngrow_distance<>4 error in levelstrip"
       stop
      endif
      if (ngrow_distance_in.ne.4) then
       print *,"ngrow_distance_in<>4 error in levelstrip"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"ngrow_make_distance<>3 error in levelstrip"
       stop
      endif
       ! 7x7x7 stencil guarantees linearity preserving property.
      if (ngrow_make_distance_accept.eq.3) then
       ! do nothing
      else
       print *,"ngrow_make_distance_accept<>3 in fort_levelstrip"
       stop
      endif

      maskfab_ptr=>maskfab
      call checkbound_array(fablo,fabhi,maskfab_ptr,ngrow_distance,-1)
      stenfab_ptr=>stenfab
      call checkbound_array(fablo,fabhi,stenfab_ptr,ngrow_distance,-1)
      vofrecon_ptr=>vofrecon
      call checkbound_array(fablo,fabhi,vofrecon_ptr,ngrow_distance,-1)
      call checkbound_array(fablo,fabhi,newfab_ptr,1,-1)
      call checkbound_array(fablo,fabhi,touchfab_ptr,0,-1)
      crsetouch_ptr=>crsetouch
      crsedist_ptr=>crsedist
      call checkbound_array(fablo,fabhi,crsetouch_ptr,0,-1)
      call checkbound_array(fablo,fabhi,crsedist_ptr,0,-1)
      
      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid level strip"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.eq.SDIM+1) then
       ! do nothing
      else
       print *,"ngeom_raw invalid level strip"
       print *,"ngeom_raw=",ngeom_raw
       stop
      endif

      if (rz_flag.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (rz_flag.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"rz_flag invalid in levelstrip"
       stop
      endif
     
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

       ! the maximum distance error that is created by using a number
       ! less than 1/2 here is VOFTOL^{1/d} where d=2 or 3.
      VFRAC_STENCIL_CUTOFF=half-VOFTOL

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0)

       ! initialize LS
       ! if is_rigid==0, then use coarse data or init to + or - bigdist
       ! if is_rigid==1, then use existing solid dist funct.
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridsten_level(xsten_donate,i,j,k,level,nhalf)

        ! im_crit fluid material occupies the center of the cell.
       do dir=1,3
        istar_array(dir)=0
       enddo
       call put_istar(istar,istar_array) 
       im_crit=NINT(stenfab(D_DECL(i,j,k),istar))

       if ((i.eq.i_DEB_DIST).and. &
           (j.eq.j_DEB_DIST).and. &
           (k.eq.k_DEB_DIST)) then
        print *,"DEB_DIST: i,j,k=",i,j,k
        print *,"DEB_DIST: x,y,z=",xsten_donate(0,1),xsten_donate(0,2), &
         xsten_donate(0,SDIM)
        print *,"DEB_DIST: im_crit=",im_crit
       endif 

       do im=1,num_materials

        if (is_rigid(im).eq.0) then

         crse_dist_valid=1
         if (level.eq.0) then
          crse_dist_valid=0
         else if ((level.gt.0).and.(level.le.finest_level)) then
          init_dist_from_crse=crsedist(D_DECL(i,j,k),im)
          if ((im.eq.im_crit).and. &
              (init_dist_from_crse.lt.zero)) then
           crse_dist_valid=0
          else if ((im.ne.im_crit).and. &
                   (init_dist_from_crse.ge.zero)) then
           crse_dist_valid=0
          else if ((im.eq.im_crit).and. &
                   (init_dist_from_crse.ge.zero)) then
           ! do nothing
          else if ((im.ne.im_crit).and. &
                   (init_dist_from_crse.le.zero)) then
           ! do nothing
          else
           print *,"im,im_crit, or init_dist_from_crse invalid: ", &
                 im,im_crit,init_dist_from_crse
           stop
          endif 
          ctouch=NINT(crsetouch(D_DECL(i,j,k),im))
          if (ctouch.eq.0) then
           crse_dist_valid=0
          else if ((ctouch.eq.1).or.(ctouch.eq.2)) then
           ! do nothing
          else
           print *,"ctouch invalid: ",ctouch
           stop
          endif
         else
          print *,"level invalid: ",level
          stop
         endif
        
         if (crse_dist_valid.eq.1) then
          do dir=1,SDIM
           dir2=num_materials+(im-1)*SDIM+dir
           newfab(D_DECL(i,j,k),dir2)=crsedist(D_DECL(i,j,k),dir2)
          enddo
          init_dist_from_crse=crsedist(D_DECL(i,j,k),im)
           ! touchfab=1 => modify with candidate distance even if 
           ! new distance will be larger.
          touchfab(D_DECL(i,j,k),im)=one
          if (init_dist_from_crse.lt.minLS(im)) then
           minLS(im)=init_dist_from_crse
          endif
          if (init_dist_from_crse.gt.maxLS(im)) then
           maxLS(im)=init_dist_from_crse
          endif
         else if (crse_dist_valid.eq.0) then
          if (im.eq.im_crit) then
           init_dist_from_crse=max_problen
          else
           init_dist_from_crse=-max_problen
          endif
         else
          print *,"crse_dist_valid invalid: ",crse_dist_valid
          stop
         endif
 
         newfab(D_DECL(i,j,k),im)=init_dist_from_crse

         if ((i.eq.i_DEB_DIST).and. &
             (j.eq.j_DEB_DIST).and. &
             (k.eq.k_DEB_DIST)) then
          print *,"DEB_DIST: im,init_dist_from_crse=", &
                  im,init_dist_from_crse
         endif 

        else if (is_rigid(im).eq.1) then

         init_dist_from_crse=newfab(D_DECL(i,j,k),im)
         touchfab(D_DECL(i,j,k),im)=two
         if (init_dist_from_crse.lt.minLS(im)) then
          minLS(im)=init_dist_from_crse
         endif
         if (init_dist_from_crse.gt.maxLS(im)) then
          maxLS(im)=init_dist_from_crse
         endif

        else
         print *,"is_rigid invalid MOF_REDIST_3D.F90"
         stop
        endif

       enddo ! im=1..num_materials

      enddo
      enddo
      enddo  ! initialize + or -

       ! ngrow_make_distance=3
      call growntilebox_TILE(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_make_distance)
 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       on_border=1-NINT(maskfab(D_DECL(i,j,k),3))

       nprocessed=nprocessed+1

       call gridsten_level(xsten_donate,i,j,k,level,nhalf)

       ! im=1..num_materials: donateflag(im)=1 => find closest distance to the 
       !             im interface.
       ! donateflag(num_materials+1)=number fluid materials in the cell.
       ! donateflag(num_materials+2 ... num_materials+1+nstar)=
       ! fluid material id that owns
       ! the respective star stencil position.
       ! If (cell_test(im)==1) and 
       !    (the LOCAL star point is owned by material im) then
       !  donateflag(num_materials+1+istar)=im 

       do im=1,num_materials+1+nstar
        donateflag(im)=0
       enddo

       fluid_materials_in_cell_stencil=0

       do im=1,num_materials*ngeom_recon
        mofdata(im)=vofrecon(D_DECL(i,j,k),im)
       enddo

       do im=1,num_materials
        vofcomp=(im-1)*ngeom_recon+1
        vcenter(im)=mofdata(vofcomp)
       enddo

       FSI_exclude=1
       call sort_volume_fraction(vcenter,FSI_exclude,sorted_list)
       im_crit=sorted_list(1)
       if ((im_crit.ge.1).and.(im_crit.le.num_materials)) then
        ! do nothing
       else
        print *,"im_crit invalid: ",im_crit
        stop
       endif
       if (is_rigid(im_crit).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid MOF_REDIST_3D.F90"
        stop
       endif

        ! in Sussman and Puckett, 
        ! a reconstructed interface segment was valid if
        ! LS_{ij}(LS_{ij}+LS_{i+i',j+j'}) <= 0 for
        ! some |i'|<=1 |j'|<=1
        ! in the present algorithm, a reconstructed segment is valid if
        ! (A) F>EPS3 and
        ! (B)        
       do im=1,num_materials
        cell_test(im)=0 !F>EPS3?
        full_neighbor(im)=0 !neighbor F(im)>1-EPS3 ?
        stencil_test(im)=0 ! F(im)>1/2-eps on cell bdry?
       enddo

         ! initialize: cell_test
       do im=1,num_materials
        if (is_rigid(im).eq.0) then
         if (vcenter(im).gt.EPS3) then
          cell_test(im)=1
         endif
        else if (is_rigid(im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid MOF_REDIST_3D.F90"
         stop
        endif

        if ((i.eq.i_DEB_DIST).and. &
            (j.eq.j_DEB_DIST).and. &
            (k.eq.k_DEB_DIST)) then
         print *,"DEB_DIST: im,cell_test=",im,cell_test(im)
        endif 

       enddo  ! im=1..num_materials

       do dir=1,3
        istar_array(dir)=0
       enddo
       call put_istar(istar,istar_array)
       im_test_center=NINT(stenfab(D_DECL(i,j,k),istar))
       if ((im_test_center.ge.1).and. &
           (im_test_center.le.num_materials)) then
        ! do nothing
       else
        print *,"im_test_center out of range (0): ",im_test_center
        stop
       endif
        ! im_test_center must be a fluid material.
       if (is_rigid(im_test_center).eq.0) then 
        ! do nothing
       else
        print *,"is_rigid(im_test_center).ne.0 (0)"
        stop
       endif
        ! 1..num_materials,fluid materials in cell, nstar
       donateflag(num_materials+1+istar)=im_test_center

       if ((i.eq.i_DEB_DIST).and. &
           (j.eq.j_DEB_DIST).and. &
           (k.eq.k_DEB_DIST)) then
        print *,"DEB_DIST: im_test=",im_test_center
       endif 

       rigid_in_stencil=0

       if (on_border.eq.0) then

         ! stencil_test and rigid_in_stencil
         ! stencil_test(im)=1 if the interpolated volume fraction
         ! on one of the 27 cell boundary points exceeds a threshold.
        do k1=klosten+2,khisten+2
        do j1=1,3
        do i1=1,3
         isten=i+i1-2
         jsten=j+j1-2
         ksten=k+k1-2

         icur(1)=i1-2 
         icur(2)=j1-2 
         icur(3)=k1-2 
         theta_nbr=zero
         theta_cen=zero
         do dir=1,SDIM
          xnbr(dir)=xsten_donate(2*icur(dir),dir)
          xmid(dir)=xsten_donate(icur(dir),dir)
          theta_nbr=theta_nbr+(xnbr(dir)-xmid(dir))**2
          theta_cen=theta_cen+(xmid(dir)-xsten_donate(0,dir))**2
         enddo
         theta_nbr=sqrt(theta_nbr)
         theta_cen=sqrt(theta_cen)

         do im=1,num_materials
          vofcomp=(im-1)*ngeom_recon+1
          VFRAC_TEMP=vofrecon(D_DECL(isten,jsten,ksten),vofcomp)
          if (is_rigid(im).eq.0) then
           if ((icur(1).eq.0).and.(icur(2).eq.0).and.(icur(3).eq.0)) then
            VFRAC_INTERP=vcenter(im)
           else if ((theta_nbr.gt.zero).and.(theta_cen.gt.zero)) then
            VFRAC_INTERP=(theta_nbr*vcenter(im)+theta_cen*VFRAC_TEMP)/ &
                         (theta_nbr+theta_cen)
           else
            print *,"theta_nbr or theta_cen invalid"
            stop
           endif
           if (VFRAC_INTERP.ge.VFRAC_STENCIL_CUTOFF) then
            stencil_test(im)=1
           endif
          else if (is_rigid(im).eq.1) then
           if (VFRAC_TEMP.ge.EPS3) then
            rigid_in_stencil=1
           endif
          else
           print *,"is_rigid(im) invalid"
           stop
          endif
         enddo  ! im=1..num_materials

        enddo
        enddo
        enddo  ! i1,j1,k1 (test VFRAC_INTERP>=VFRAC_STENCIL_CUTOFF)

        if ((i.eq.i_DEB_DIST).and. &
            (j.eq.j_DEB_DIST).and. &
            (k.eq.k_DEB_DIST)) then
         do im=1,num_materials
          print *,"DEB_DIST: im,stencil_test=",im,stencil_test(im)
         enddo
        endif

         ! investigate all points coinciding at the intersection of the
         ! line connecting cells (i,j,k) and (i+i3,j+j3,k+k3) with the
         ! boundary of cell (i,j,k). 
         ! 
        do k3=klosten,khisten
        do j3=-1,1
        do i3=-1,1
         istar_array(1)=i3
         istar_array(2)=j3
         istar_array(3)=k3
         do dir=1,3
          if (istar_array(dir).eq.0) then
           i4low(dir)=0
           i4high(dir)=0
          else if (istar_array(dir).eq.-1) then
           i4low(dir)=-1
           i4high(dir)=0
          else if (istar_array(dir).eq.1) then
           i4low(dir)=0
           i4high(dir)=1
          else
           print *,"istar_array bust"
           stop
          endif
         enddo ! dir=1..3
         im_corner=0
         do im=1,num_materials
          FSUM(im)=zero
         enddo
         ! (i4,j4,k4)=(0,0,0) if (i3,j3,k3)=(0,0,0)
         do k4=i4low(3),i4high(3)
         do j4=i4low(2),i4high(2)
         do i4=i4low(1),i4high(1)
          isten=i4+i
          jsten=j4+j
          ksten=k4+k
          
          i4_array(1)=i4
          i4_array(2)=j4
          i4_array(3)=k4
          do dir=1,3
           if (i4_array(dir).eq.0) then
            istar_array_offset(dir)=istar_array(dir)
           else if (i4_array(dir).eq.1) then
            istar_array_offset(dir)=-1
           else if (i4_array(dir).eq.-1) then
            istar_array_offset(dir)=1
           else
            print *,"i4_array invalid"
            stop
           endif
          enddo ! dir

           ! fort_steninit called with 
           ! tessellate==0 prior to fort_levelstrip.
           ! im_test_stencil is the material that occupies the
           ! node in question.
          call put_istar(istar,istar_array_offset)
          im_test_stencil=NINT(stenfab(D_DECL(isten,jsten,ksten),istar))
          if ((im_test_stencil.ge.1).and. &
              (im_test_stencil.le.num_materials)) then
           ! do nothing
          else
           print *,"im_test_stencil out of range 1: ",im_test_stencil
           stop
          endif
          if (is_rigid(im_test_stencil).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(im_test_stencil).ne.0 (1)"
           stop
          endif

          if ((i.eq.i_DEB_DIST).and. &
              (j.eq.j_DEB_DIST).and. &
              (k.eq.k_DEB_DIST)) then
           print *,"DEB_DIST: i3,j3,k3 ",i3,j3,k3
           print *,"DEB_DIST: i4,j4,k4 ",i4,j4,k4
           print *,"DEB_DIST: im_test_stencil ",im_test_stencil
          endif

          if (im_corner.eq.0) then
           im_corner=im_test_stencil
          else if ((im_corner.ge.1).and.(im_corner.le.num_materials)) then
           if (im_corner.eq.im_test_stencil) then
            ! do nothing
           else
            im_corner=-1
           endif
          else if (im_corner.eq.-1) then
           ! do nothing
          else
           print *,"im_corner invalid: ",im_corner
           stop
          endif
          do im=1,num_materials
           vofcomp=(im-1)*ngeom_recon+1
           VFRAC_TEMP=vofrecon(D_DECL(isten,jsten,ksten),vofcomp)
           FSUM(im)=FSUM(im)+VFRAC_TEMP
          enddo
         enddo ! k4
         enddo ! j4
         enddo ! i4

            ! if im_corner=-1, then there is a jump in the material
            ! type at the node, so the material type is assigned
            ! to be the material with the dominant "nodal volume fraction"
         if (im_corner.eq.-1) then
          FSI_exclude=1
          call sort_volume_fraction(FSUM,FSI_exclude,sorted_list)
          im_corner=sorted_list(1)
          if (is_rigid(im_corner).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(im_corner).ne.0"
           stop
          endif
         else if (im_corner.eq.0) then
          print *,"im_corner invalid"
          stop
         else if ((im_corner.ge.1).and.(im_corner.le.num_materials)) then
          ! do nothing
         else
          print *,"im_corner invalid"
          stop
         endif
         if ((im_corner.ge.1).and. &
             (im_corner.le.num_materials)) then

           ! a tougher test for flotsam near contact lines.
          if ((stencil_test(im_corner).eq.1).or. &
              (rigid_in_stencil.eq.0)) then
           call put_istar(istar,istar_array) 
           donateflag(num_materials+1+istar)=im_corner
          else if ((stencil_test(im_corner).eq.0).and. &
                   (rigid_in_stencil.eq.1)) then
           ! do nothing
          else
           print *,"stencil_test or rigid_in_stencil invalid"
           stop
          endif

         else
          print *,"im_corner invalid"
          stop
         endif

         if ((i.eq.i_DEB_DIST).and. &
             (j.eq.j_DEB_DIST).and. &
             (k.eq.k_DEB_DIST)) then
          call put_istar(istar,istar_array) 
          print *,"DEB_DIST:ijk 3,istar,donateflag(num_materials+1+istar) ", &
           i3,j3,k3,istar,donateflag(num_materials+1+istar)
         endif

        enddo
        enddo
        enddo  ! i3,j3,k3

         ! full_neighbor
        do k3=klosten,khisten
        do j3=-1,1
        do i3=-1,1
         if ((i3.eq.0).and.(j3.eq.0).and.(k3.eq.0)) then
          ! do nothing
         else if ((abs(i3).eq.1).or.(abs(j3).eq.1).or.(abs(k3).eq.1)) then
          iside=i+i3
          jside=j+j3
          kside=k+k3
          istar_array(1)=i3
          istar_array(2)=j3
          istar_array(3)=k3

          do im=1,num_materials
           if (is_rigid(im).eq.0) then
            vofcomp=(im-1)*ngeom_recon+1
            VFRAC_TEMP=vofrecon(D_DECL(iside,jside,kside),vofcomp)
            if (VFRAC_TEMP.ge.one-EPS3) then
             full_neighbor(im)=1
             call put_istar(istar,istar_array) 
             donateflag(num_materials+1+istar)=im
            endif
           else if (is_rigid(im).eq.1) then
            ! do nothing
           else
            print *,"is_rigid invalid MOF_REDIST_3D.F90"
            stop
           endif
          enddo ! im=1..num_materials 

         else
          print *,"i3,j3,k3 bust"
          stop
         endif
        enddo
        enddo
        enddo ! i3,j3,k3

       else if (on_border.eq.1) then

        ! do nothing 

       else
        print *,"on_border invalid: ",on_border
        stop
       endif

       do im=1,num_materials

        if (is_rigid(im).eq.0) then

         keep_flotsam=0
         if (cell_test(im).eq.1) then !F_{im}>EPS3?
          keep_flotsam=1
         else if (cell_test(im).eq.0) then !F_{im}<EPS3?
          keep_flotsam=0
         else
          print *,"cell_test invalid: ",im,cell_test(im)
          stop
         endif

         legitimate_material=0
         if ((vcenter(im).ge.half).or. &
             (im.eq.im_crit).or. & !im_crit=argmax_{im} F_{im}
             (full_neighbor(im).eq.1).or. &
             (keep_flotsam.eq.1)) then!keep_flotsam=1 if Fm>EPS3 and no trunc.
          legitimate_material=1
         else if ((vcenter(im).le.half).and. &
                  (im.ne.im_crit).and. &
                  (full_neighbor(im).eq.0).and. &
                  (keep_flotsam.eq.0)) then
          legitimate_material=0
         else
          print *,"legitimate check failed"
          stop
         endif

         if (legitimate_material.eq.1) then

          do k3=klosten,khisten
          do j3=-1,1
          do i3=-1,1
           istar_array(1)=i3
           istar_array(2)=j3
           istar_array(3)=k3
           call put_istar(istar,istar_array)
           if ((istar.ge.1).and.(istar.le.nstar)) then
            im_test_stencil=NINT(stenfab(D_DECL(i,j,k),istar))
            if (is_rigid(im_test_stencil).eq.0) then
             if (im_test_stencil.eq.im) then
               ! 1<=istar<=nstar
              donateflag(num_materials+1+istar)=im
             endif
            else
             print *,"is_rigid(im_test_stencil).ne.0 (1)"
             stop
            endif
           else
            print *,"istar invalid: ",istar
            stop
           endif
          enddo ! k3
          enddo ! j3
          enddo ! i3

         else if (legitimate_material.eq.0) then
          ! do nothing
         else
          print *,"legitimate_material invalid 2902: ",legitimate_material
          stop
         endif

         if (legitimate_material.eq.1) then
          fluid_materials_in_cell_stencil= &
           fluid_materials_in_cell_stencil+1
          donateflag(im)=1
         else if (legitimate_material.eq.0) then
          ! do nothing
         else
          print *,"legitimate_material invalid 2913: ",legitimate_material
          stop
         endif 

        else if (is_rigid(im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid MOF_REDIST_3D.F90"
         stop
        endif
       enddo ! im=1..num_materials

       donateflag(num_materials+1)=fluid_materials_in_cell_stencil

       if ((i.eq.i_DEB_DIST).and. &
           (j.eq.j_DEB_DIST).and. &
           (k.eq.k_DEB_DIST)) then
        print *,"DEB_DIST: fluid_materials_in_cell_stencil=", &
         fluid_materials_in_cell_stencil
       endif

       if ((fluid_materials_in_cell_stencil.lt.1).or. &
           (fluid_materials_in_cell_stencil.gt.num_materials)) then
        print *,"fluid_materials_in_cell_stencil invalid:", &
          fluid_materials_in_cell_stencil
        stop
       else if (fluid_materials_in_cell_stencil.eq.1) then
         ! do nothing
       else if ((fluid_materials_in_cell_stencil.ge.2).and. &
                (fluid_materials_in_cell_stencil.le.num_materials)) then

        ilocut(3)=0
        ihicut(3)=0
        icur(1)=i
        icur(2)=j
        icur(3)=k

         ! The main loop stencil,
         ! traversing valid interface segments,
         ! is derived from "growntilebox_TILE(ngrow_make_distance=3)
         ! This inner loop cannot extend beyond the "owned"
         ! tile otherwise there will be competition between
         ! threads for modifying the same piece of data.
         ! i.e. if tile A has the interface and tile B distance
         ! needs to be modified, then the tile B distance should
         ! only be touched by "thread B"
         ! 
        do dir=1,SDIM
         ilocut(dir)=-ngrow_make_distance_accept
         ihicut(dir)=ngrow_make_distance_accept
         if (icur(dir)+ilocut(dir).lt.tilelo(dir)) then
          ilocut(dir)=tilelo(dir)-icur(dir)
         endif
         if (icur(dir)+ihicut(dir).gt.tilehi(dir)) then
          ihicut(dir)=tilehi(dir)-icur(dir)
         endif
        enddo ! dir=1..sdim

         !LSslope_center is normalized
        call get_primary_slope( &
         bfact,dx,xsten_donate,nhalf, &
         mofdata, &
         LSslope_center,imslope_center,SDIM) 

        bypass_cutoff=sqrt(three)/two

        do k1=ilocut(3),ihicut(3)
        do j1=ilocut(2),ihicut(2)
        do i1=ilocut(1),ihicut(1)

         call gridsten_level(xsten_accept,i+i1,j+j1,k+k1,level,nhalf)

         bypass_update_closest=1

         if (1.eq.0) then
          bypass_update_closest=0
         endif

         if (levelrz.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if (levelrz.eq.COORDSYS_RZ) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
          bypass_update_closest=0
         else
          print *,"levelrz invalid in fort_levelstrip: ",levelrz
          stop
         endif

         if ((abs(i1).le.1).and.(abs(j1).le.1).and.(abs(k1).le.1)) then
          bypass_update_closest=0
         endif

         if (fluid_materials_in_cell_stencil.eq.2) then
          ! do nothing
         else if ((fluid_materials_in_cell_stencil.ge.3).and. &
                  (fluid_materials_in_cell_stencil.le.num_materials)) then
          bypass_update_closest=0
         else
          print *,"fluid_materials_in_cell_stencil invalid: ", &
            fluid_materials_in_cell_stencil
          stop
         endif

         if (imslope_center.eq.0) then
          bypass_update_closest=0
         else if ((imslope_center.ge.1).and. &
                  (imslope_center.le.num_materials)) then
          ! do nothing
         else
          print *,"imslope_center invalid: ",imslope_center
          stop
         endif

         if (bypass_update_closest.eq.1) then

          crude_dist=zero
          dotprod=zero
          do dir=1,SDIM
           crude_normal=xsten_donate(0,dir)-xsten_accept(0,dir)
           crude_dist=crude_dist+crude_normal**2 
           dotprod=dotprod+crude_normal*LSslope_center(dir)
          enddo
          crude_dist=sqrt(crude_dist)
          if (crude_dist.gt.zero) then
           dotprod=abs(dotprod)/crude_dist
          else
           print *,"crude_dist invalid: ",crude_dist
           stop
          endif
           ! Sean Mauch type of algorithm
          if ((dotprod.ge.zero).and.(dotprod.le.bypass_cutoff)) then
           ! do nothing
          else if ((dotprod.ge.bypass_cutoff).and. &
                   (dotprod.le.one+EPS2)) then
           bypass_update_closest=0
          else
           print *,"dotprod invalid: ",dotprod
           stop
          endif
         else if (bypass_update_closest.eq.0) then
          !do nothing
         else
          print *,"bypass_update_closest invalid: ",bypass_update_closest
          stop
         endif

         if (bypass_update_closest.eq.0) then
         
          do im=1,num_materials*(1+SDIM)
           newfab_hold(im)=newfab(D_DECL(i+i1,j+j1,k+k1),im)
          enddo  ! im
          do im=1,num_materials
           touch_hold(im)=NINT(touchfab(D_DECL(i+i1,j+j1,k+k1),im))
          enddo

          call update_closest( &
           xsten_accept, &
           xsten_donate, &
           nhalf, &
           dx,xlo,bfact, &
           level, &
           fablo, &
           mofdata, & !intent(in)
           LSslope_center, &
           imslope_center, &
           nstar, &
           i,j,k, &  ! donate index
           i1,j1,k1, & ! accept index: i+i1,j+j1,k+k1
           newfab_hold, & !intent(inout)
           touch_hold, &  !intent(inout)
           minLS, & !intent(inout)
           maxLS, & !intent(inout)
           donateflag, & !intent(in)
           time)

          do im=1,num_materials*(1+SDIM)
           newfab(D_DECL(i+i1,j+j1,k+k1),im)=newfab_hold(im)
          enddo  ! im
          do im=1,num_materials
           touchfab(D_DECL(i+i1,j+j1,k+k1),im)=touch_hold(im)
          enddo

         else if (bypass_update_closest.eq.1) then
          !do nothing
         else
          print *,"bypass_update_closest invalid: ",bypass_update_closest
          stop
         endif

        enddo
        enddo
        enddo ! i1,j1,k1

       else
        print *,"parameter bust in fort_levelstrip"
        print *,"fluid_materials_in_cell_stencil=", &
          fluid_materials_in_cell_stencil
        stop
       endif

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine fort_levelstrip

      subroutine fort_correct_uninit( &
         minLS, &
         maxLS, &
         max_problen, &
         level, &
         finest_level, &
         newfab,DIMS(newfab), &
         touchfab,DIMS(touchfab), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         xlo,dx, &
         time) &
      bind(c,name='fort_correct_uninit')

      use global_utility_module
      use probcommon_module
      use mof_redist_module

      IMPLICIT NONE

      integer, INTENT(in) :: level,finest_level
      real(amrex_real), INTENT(in) :: minLS(num_materials)
      real(amrex_real), INTENT(in) :: maxLS(num_materials)
      real(amrex_real), INTENT(in) :: max_problen
      integer, INTENT(in) :: DIMDEC(newfab)
      integer, INTENT(in) :: DIMDEC(touchfab)

      real(amrex_real), INTENT(inout), target :: &
        newfab(DIMV(newfab),num_materials*(1+SDIM))
      real(amrex_real), pointer :: newfab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: &
              touchfab(DIMV(touchfab),num_materials)

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      real(amrex_real), INTENT(in) :: time

      integer i,j,k
      integer im
      integer ctouch
      real(amrex_real) init_dist
     
      newfab_ptr=>newfab

      if (bfact.lt.1) then
       print *,"bfact invalid143"
       stop
      endif
      if (max_problen.le.zero) then
       print *,"max_problen invalid"
       stop
      endif

      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid in levelstrip"
       stop
      endif

      call checkbound_array(fablo,fabhi,newfab_ptr,1,-1)
      call checkbound_array(fablo,fabhi,touchfab,0,-1)
      
      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0)
 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       do im=1,num_materials

         ! touchfab=1 => modify with candidate distance even if 
         ! new distance will be larger.
        ctouch=NINT(touchfab(D_DECL(i,j,k),im))

        if (is_rigid(im).eq.0) then

         if (ctouch.eq.0) then
          init_dist=newfab(D_DECL(i,j,k),im)
          if (init_dist.lt.zero) then
           if (minLS(im).ge.zero) then
            init_dist=-max_problen
           else
            init_dist=minLS(im)
           endif
          else if (init_dist.ge.zero) then
           if (maxLS(im).lt.zero) then
            init_dist=max_problen
           else
            init_dist=maxLS(im)
           endif
          else
           print *,"init_dist bust: ",init_dist
           stop
          endif
          newfab(D_DECL(i,j,k),im)=init_dist
         else if ((ctouch.eq.1).or.(ctouch.eq.2)) then
          ! do nothing
         else
          print *,"ctouch invalid: ",ctouch
          stop
         endif
 
        else if (is_rigid(im).eq.1) then

         if (ctouch.ne.2) then
          print *,"ctouch invalid: ",ctouch
          stop
         endif

        else
         print *,"is_rigid invalid MOF_REDIST_3D.F90"
         stop
        endif

       enddo ! im=1..num_materials

      enddo
      enddo
      enddo  ! replace uninit with minLS or maxLS

      return
      end subroutine fort_correct_uninit


      subroutine fort_steninit( &
       level, &
       finest_level, &
       stenfab,DIMS(stenfab), &
       maskfab,DIMS(maskfab), &
       vofrecon,DIMS(vofrecon), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       rz_flag, &
       xlo,dx, &
       time, &
       ngrow_distance_in, &
       nstar) &
      bind(c,name='fort_steninit')

      use global_utility_module
      use probcommon_module
      use MOF_routines_module
      use mof_redist_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nstar
      integer, INTENT(in) :: ngrow_distance_in
      integer, INTENT(in) :: DIMDEC(stenfab)
      integer, INTENT(in) :: DIMDEC(maskfab)
      integer, INTENT(in) :: DIMDEC(vofrecon)

      real(amrex_real), INTENT(out), target :: stenfab(DIMV(stenfab),nstar)
      real(amrex_real), pointer :: stenfab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: maskfab(DIMV(maskfab),2)
      real(amrex_real), pointer :: maskfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
        vofrecon(DIMV(vofrecon),num_materials*ngeom_recon)
      real(amrex_real), pointer :: vofrecon_ptr(D_DECL(:,:,:),:)

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: rz_flag
      real(amrex_real), INTENT(in) :: time

      integer i,j,k
      real(amrex_real) vcenter(num_materials)
      real(amrex_real) xsten(-3:3,SDIM)
      integer nhalf

      integer klosten,khisten
      integer im
      integer vofcomp
      integer istar_array(3)
      integer nstar_test
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      integer istar
      integer im_crit
      integer i3,j3,k3
      real(amrex_real) xcorner(SDIM)
      integer im_test
      integer dir
      integer mask1,mask2
      integer tessellate
 
      nhalf=3      
   
      stenfab_ptr=>stenfab

      tessellate=0
 
      if (bfact.lt.1) then
       print *,"bfact invalid144"
       stop
      endif

      nstar_test=9
      if (SDIM.eq.3) then
       nstar_test=nstar_test*3
      endif
      if (nstar_test.ne.nstar) then
       print *,"nstar invalid steninit nstar nstar_test ",nstar,nstar_test
       stop
      endif

      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid in steninit"
       stop
      endif

      if (ngrow_distance.ne.4) then
       print *,"ngrow_distance<>4 error in steninit"
       stop
      endif
      if (ngrow_distance_in.ne.4) then
       print *,"ngrow_distance_in<>4 error in steninit"
       stop
      endif

      call checkbound_array(fablo,fabhi,stenfab_ptr,ngrow_distance,-1)
      maskfab_ptr=>maskfab
      call checkbound_array(fablo,fabhi,maskfab_ptr,ngrow_distance,-1)
      vofrecon_ptr=>vofrecon
      call checkbound_array(fablo,fabhi,vofrecon_ptr,ngrow_distance,-1)
      
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid steninit"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid steninit"
       print *,"ngeom_raw=",ngeom_raw
       stop
      endif

      if (rz_flag.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (rz_flag.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.COORDSYS_CYLINDRICAL) then
       ! do nothing
      else
       print *,"rz_flag invalid in steninit"
       stop
      endif

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

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_distance) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       ! mask1=1 at interior cells or fine/fine ghost cells
       ! mask1=0 at coarse/fine ghost cells or outside domain.
       ! mask2=1 at interior cells
       mask1=NINT(maskfab(D_DECL(i,j,k),1))
       mask2=NINT(maskfab(D_DECL(i,j,k),2))

       if ((mask2.eq.1).or.(mask1.eq.0)) then

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do im=1,num_materials*ngeom_recon
         mofdata(im)=vofrecon(D_DECL(i,j,k),im)
        enddo

         ! vcenter = volume fraction 
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         vcenter(im)=mofdata(vofcomp)
        enddo ! im

         ! uses VOFTOL
        call check_full_cell_vfrac( &
          vcenter, &
          tessellate, &  ! =0
          im_crit)

        if ((im_crit.ge.1).and.(im_crit.le.num_materials)) then

         do istar=1,nstar
          stenfab(D_DECL(i,j,k),istar)=im_crit
         enddo

        else if (im_crit.eq.0) then

         do k3=klosten,khisten
         do j3=-1,1
         do i3=-1,1

          dir=1
          xcorner(dir)=xsten(i3,dir)
          dir=2
          xcorner(dir)=xsten(j3,dir)
          if (SDIM.eq.3) then
           dir=SDIM
           xcorner(dir)=xsten(k3,dir)
          endif

          istar_array(1)=i3
          istar_array(2)=j3
          istar_array(3)=k3

          call multi_get_volumePOINT( &
           tessellate, &  ! =0
           bfact,dx,xsten,nhalf, &
           mofdata,xcorner, &
           im_test,SDIM)
          if ((im_test.ge.1).and.(im_test.le.num_materials)) then
           call put_istar(istar,istar_array) 
           stenfab(D_DECL(i,j,k),istar)=im_test
          else
           print *,"im_test invalid: ",im_test
           stop
          endif
         enddo
         enddo
         enddo  ! i3,j3,k3
 
        else
         print *,"im_crit invalid"
         stop
        endif

       else if ((mask2.eq.0).and.(mask1.eq.1)) then
        ! do nothing
       else
        print *,"mask invalid"
        stop
       endif

      enddo
      enddo
      enddo  !i,j,k 


      return
      end subroutine fort_steninit


       ! fort_faceinit is called from NavierStokes.cpp,
       !  NavierStokes::makeFaceFrac
       ! makeFaceFrac is called from NavierStokes3.cpp,
       !  NavierStokes::ColorSum
       ! for finding areas internal to a cell, perturb each internal 
       ! interface, find areas and volumes, then check for the difference 
       ! in volumes divided by eps times the area.
       ! 
      subroutine fort_faceinit( &
         tid, &
         tessellate, &  ! =0,1, or 3
         level, &
         finest_level, &
         facefab,DIMS(facefab), &
         maskfab,DIMS(maskfab), &
         vofrecon,DIMS(vofrecon), &
         tilelo,tilehi, &
         fablo,fabhi,bfact, &
         rz_flag, &
         xlo,dx, &
         time, &
         ngrow, &
         nface) &
      bind(c,name='fort_faceinit')

      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module
      use mof_redist_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: tessellate ! 0,1, or 3
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nface
      integer, INTENT(in) :: ngrow
      integer, INTENT(in) :: DIMDEC(facefab)
      integer, INTENT(in) :: DIMDEC(maskfab)
      integer, INTENT(in) :: DIMDEC(vofrecon)

      real(amrex_real), INTENT(out), target :: facefab(DIMV(facefab),nface)
      real(amrex_real), pointer :: facefab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: maskfab(DIMV(maskfab),2)
      real(amrex_real), INTENT(in), target :: &
          vofrecon(DIMV(vofrecon),num_materials*ngeom_recon)

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: rz_flag
      real(amrex_real), INTENT(in) :: time

      integer i,j,k
      integer im
      integer im_crit
      integer iface
      integer dir,side

      integer nface_test
      integer vofcomp
      real(amrex_real) vcenter(num_materials)
      real(amrex_real) vcenter_thin(num_materials)
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdatavalid(num_materials*ngeom_recon)
      real(amrex_real) mofdataproject(num_materials*ngeom_recon)
      real(amrex_real) localface(num_materials,SDIM,2)
      real(amrex_real) totalface(SDIM,2)

      real(amrex_real) xsten(-3:3,SDIM)
      real(amrex_real) xsten_thin(-1:1,SDIM)
      real(amrex_real) dxthin
      integer nhalf,nhalf_thin,isten
      real(amrex_real) dummy_tri(SDIM+1,SDIM)
      integer nmax,ivert
      integer dir2
      integer shapeflag
      real(amrex_real) multi_volume(num_materials)
      real(amrex_real) multi_cen(SDIM,num_materials)
      real(amrex_real) multi_area(num_materials)
      real(amrex_real) total_vol
      integer mask1,mask2
      integer normalize_tessellate
      integer local_tessellate
      integer, parameter :: continuous_mof=STANDARD_MOF
      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))

      facefab_ptr=>facefab

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if ((tessellate.ne.0).and. &
          (tessellate.ne.1).and. &
          (tessellate.ne.3)) then
       print *,"tessellate invalid45"
       stop
      endif

      nhalf=3 
      nhalf_thin=1 
      nmax=POLYGON_LIST_MAX ! in: FACEINIT

      do ivert=1,SDIM+1
      do dir2=1,SDIM
       dummy_tri(ivert,dir2)=zero
      enddo
      enddo
      if (bfact.lt.1) then
       print *,"bfact invalid145"
       stop
      endif

       ! area for all faces of a cell.
       ! (num_materials,sdim,2)
      nface_test=num_materials*SDIM*2
      if (nface_test.ne.nface) then
       print *,"nface invalid faceinit nface nface_test ",nface,nface_test
       stop
      endif

      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid in faceinit"
       stop
      endif

      if (ngrow.lt.1) then
       print *,"ngrow<1 error in faceinit"
       stop
      endif

      call checkbound_array(fablo,fabhi,facefab_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,maskfab,ngrow,-1)
      call checkbound_array(fablo,fabhi,vofrecon,ngrow,-1)
      
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid faceinit"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid faceinit"
       print *,"ngeom_raw=",ngeom_raw
       stop
      endif

      if (rz_flag.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (rz_flag.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.COORDSYS_CYLINDRICAL) then 
       ! do nothing
      else
       print *,"rz_flag invalid in faceinit"
       stop
      endif
 
      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

        ! mask1=1 at interior cells or fine/fine ghost cells
        ! mask2=1 at interior cells
       mask1=NINT(maskfab(D_DECL(i,j,k),1))
       mask2=NINT(maskfab(D_DECL(i,j,k),2))

       if ((mask2.eq.1).or.(mask1.eq.0)) then

         ! xsten(0,dir)  dir=1,2  is center of cell
         ! xsten(1,dir)  is right coordinate in dir direction
         ! xsten(-1,dir) is left coordinate in dir direction.
        call gridsten_level(xsten,i,j,k,level,nhalf)

        do im=1,num_materials*ngeom_recon
         mofdata(im)=vofrecon(D_DECL(i,j,k),im)
        enddo

        normalize_tessellate=0
        call make_vfrac_sum_ok_copy( &
         cmofsten, &
         xsten,nhalf, &
         continuous_mof, &
         bfact,dx, &
         normalize_tessellate, &  ! =0
         mofdata,mofdatavalid, &
         SDIM)

        if (tessellate.eq.3) then

         local_tessellate=2

          !EPS2
         call multi_get_volume_tessellate( &
          tessellate, &  ! =3
          bfact,dx, &
          xsten,nhalf, &
          mofdatavalid, &
          geom_xtetlist_uncapt(1,1,1,tid+1), &
          nmax, &
          nmax, &
          SDIM) 

        else if ((tessellate.eq.0).or. &
                 (tessellate.eq.1)) then
         local_tessellate=tessellate
        else
         print *,"tessellate invalid"
         stop
        endif

          ! vcenter = volume fraction 
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         vcenter(im)=mofdatavalid(vofcomp)
        enddo ! im

        call check_full_cell_vfrac( &
          vcenter, &
          tessellate, & ! 0,1, or 3
          im_crit)

        do dir=1,SDIM
         do side=1,2
          totalface(dir,side)=zero
         enddo
        enddo

        iface=0
        do im=1,num_materials
         do dir=1,SDIM
          do side=1,2
           iface=iface+1
           localface(im,dir,side)=zero
          enddo ! side
         enddo ! dir
        enddo ! im

        if (iface.ne.nface) then
         print *,"iface invalid"
         stop
        endif

        do dir=1,SDIM
        do side=1,2
         call gridarea(xsten,nhalf,levelrz,dir-1,side-1,totalface(dir,side))
        enddo
        enddo ! dir,side

        if ((im_crit.ge.1).and.(im_crit.le.num_materials)) then

         do dir=1,SDIM
         do side=1,2
          localface(im_crit,dir,side)=one ! area fraction
         enddo
         enddo  ! dir,side

        else if (im_crit.eq.0) then

         shapeflag=0

         do dir=1,SDIM
         do side=1,2

          do isten=-1,1
           do dir2=1,SDIM
            xsten_thin(isten,dir2)=xsten(isten,dir2)
           enddo
          enddo ! isten

          dxthin=EPS_3_2*(xsten(1,dir)-xsten(-1,dir))
          if (side.eq.1) then
           xsten_thin(1,dir)=xsten(-1,dir)+dxthin
          else if (side.eq.2) then
           xsten_thin(-1,dir)=xsten(1,dir)-dxthin
          else
           print *,"side invalid"
           stop
          endif
          xsten_thin(0,dir)=half*(xsten_thin(-1,dir)+xsten_thin(1,dir))
           ! find volumes and areas (not scaled) of materials in
           ! xsten_thin box: xsten_thin(0,dir) = center of thin box
           ! xsten_thin(1,dir) right side in dir direction
           ! xsten_thin(-1,dir) left side
           ! multi_cen(sdim,num_materials) is "absolute" 
          call project_slopes_to_face( &
           bfact,dx,xsten,nhalf, &
           mofdatavalid,mofdataproject, &
           SDIM,dir,side)

           ! base case (also area fractions)
           ! multi_cen in absolute coordinate system (not relative to cell
           ! centroid)
           ! in: fort_faceinit
          call multi_get_volume_grid( &
            EPS2, &
            local_tessellate, &  ! 0,1, or 2
            bfact,dx,xsten,nhalf, &
            mofdataproject, &
            xsten_thin,nhalf_thin, &
            dummy_tri, &
            multi_volume, &
            multi_cen, &
            multi_area, &
            geom_xtetlist_uncapt(1,1,1,tid+1), &
            nmax, &
            nmax, &
            SDIM, &
            shapeflag) 

          total_vol=zero
          do im=1,num_materials
           if ((tessellate.eq.1).or. &
               (tessellate.eq.3)) then
            total_vol=total_vol+multi_volume(im)
           else if (tessellate.eq.0) then
            if (is_rigid(im).eq.0) then
             total_vol=total_vol+multi_volume(im)
            else if (is_rigid(im).eq.1) then
             ! do nothing
            else
             print *,"is_rigid invalid MOF_REDIST_3D.F90"
             stop
            endif
           else
            print *,"tessellate invalid46"
            stop
           endif
          enddo ! im=1..num_materials

          do im=1,num_materials
           if (total_vol.gt.zero) then
            vcenter_thin(im)=multi_volume(im)/total_vol
           else if (total_vol.eq.zero) then
            vcenter_thin(im)=zero
           else
            print *,"total_vol invalid: ",total_vol
            stop
           endif 
           localface(im,dir,side)=vcenter_thin(im)
          enddo ! im=1..num_materials

         enddo
         enddo ! dir,side
    
        else
         print *,"im_crit out of range: ",im_crit
         stop
        endif

        do im=1,num_materials
         do dir=1,SDIM
          do side=1,2
           localface(im,dir,side)= &
            localface(im,dir,side)*totalface(dir,side)
          enddo
         enddo
        enddo 

        iface=0
        do im=1,num_materials
         do dir=1,SDIM
          do side=1,2
           iface=iface+1
           facefab(D_DECL(i,j,k),iface)=localface(im,dir,side)
          enddo ! side
         enddo ! dir
        enddo ! im

        if (iface.ne.nface) then
         print *,"iface invalid iface,nface: ",iface,nface
         stop
        endif

       else if ((mask2.eq.0).and.(mask1.eq.1)) then
        ! do nothing
       else
        print *,"mask invalid mask1,mask2: ",mask1,mask2
        stop
       endif

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine fort_faceinit


       ! 1. NavierStokes::makeFaceFrac
       !      -> fort_faceinit
       ! 2. NavierStokes::ProcessFaceFrac
       !      -> fort_faceprocess
       ! called from: NavierStokes::ProcessFaceFrac (NavierStokes.cpp)
       ! facefab is initialized in fort_faceinit
      subroutine fort_faceprocess( &
       ngrow_source, &
       ngrow_dest, &
       tid, &
       dir, &
       tessellate, & ! 0,1, or 3
       level, &
       finest_level, &
       dstfab,DIMS(dstfab), &
       facefab,DIMS(facefab), &
       vofrecon,DIMS(vofrecon), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       rz_flag, &
       xlo,dx, &
       time, &
       nface_src,nface_dst) &
      bind(c,name='fort_faceprocess')

      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: ngrow_source
      integer, INTENT(in) :: ngrow_dest
      integer, INTENT(in) :: tid
      integer, INTENT(in) :: dir
      integer, INTENT(in) :: tessellate ! 0,1, or 3
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nface_src,nface_dst
      integer, INTENT(in) :: DIMDEC(dstfab)
      integer, INTENT(in) :: DIMDEC(facefab)
      integer, INTENT(in) :: DIMDEC(vofrecon)

      real(amrex_real), INTENT(out), target :: dstfab(DIMV(dstfab),nface_dst)
      real(amrex_real), pointer :: dstfab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: facefab(DIMV(facefab),nface_src)
      real(amrex_real), INTENT(in), target :: &
        vofrecon(DIMV(vofrecon),num_materials*ngeom_recon)

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: rz_flag
      real(amrex_real), INTENT(in) :: time

      integer i,j,k
      integer iface_left,iface_right
      integer inormal
      integer left_face_ok,right_face_ok
      integer side
      integer ii,jj,kk
      integer im,im_opp

      integer nface_src_test,nface_dst_test,vofcomp
      real(amrex_real) mofdata_left(num_materials*ngeom_recon)
      real(amrex_real) mofdata_right(num_materials*ngeom_recon)
      real(amrex_real) xsten_left(-3:3,SDIM)
      real(amrex_real) xsten_right(-3:3,SDIM)
      integer nhalf
      integer dir2
      real(amrex_real) multi_cen_left(SDIM,num_materials)
      real(amrex_real) multi_cen_right(SDIM,num_materials)
      real(amrex_real) cencell_left(SDIM)
      real(amrex_real) cencell_right(SDIM)
      real(amrex_real) volcell_left
      real(amrex_real) volcell_right
      real(amrex_real) leftface(num_materials)
      real(amrex_real) rightface(num_materials)
      real(amrex_real) frac_left(num_materials)
      real(amrex_real) frac_right(num_materials)
      real(amrex_real) left_total
      real(amrex_real) right_total
      real(amrex_real) vol_total
      integer ml,mr
      real(amrex_real) frac_pair(num_materials,num_materials)  !(m_left,m_right)
      real(amrex_real) dist_pair(num_materials,num_materials)
      real(amrex_real) dpair
      real(amrex_real) delta,RR
      real(amrex_real) xstenMAC(-1:1,SDIM)
      integer nhalfMAC
      integer at_RZ_face
      real(amrex_real) L_face
      integer nmax
      integer local_tessellate

      dstfab_ptr=>dstfab

      L_face=dx(dir+1)

      nmax=POLYGON_LIST_MAX  ! in: FACE_PROCESS

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif
 
      nhalf=3 
      nhalfMAC=1 

      if (bfact.lt.1) then
       print *,"bfact invalid87"
       stop
      endif
      if ((tessellate.ne.0).and. &
          (tessellate.ne.1).and. &
          (tessellate.ne.3)) then  ! raster
       print *,"tessellate invalid51"
       stop
      endif

      nface_src_test=num_materials*SDIM*2
      if (nface_src_test.ne.nface_src) then
       print *,"nface_src invalid nface_src nface_src_test ", &
         nface_src,nface_src_test
       stop
      endif
      nface_dst_test=num_materials*num_materials*2
      if (nface_dst_test.ne.nface_dst) then
       print *,"nface_dst invalid nface_dst nface_dst_test ", &
         nface_dst,nface_dst_test
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in faceprocess"
       stop
      endif

      call checkbound_array(fablo,fabhi,dstfab_ptr,ngrow_dest,dir)
      call checkbound_array(fablo,fabhi,facefab,ngrow_source,-1)
      call checkbound_array(fablo,fabhi,vofrecon,ngrow_source,-1)
      
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid faceprocess"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid faceprocess"
       print *,"ngeom_raw=",ngeom_raw
       stop
      endif

      if (rz_flag.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (rz_flag.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.COORDSYS_CYLINDRICAL) then 
       ! do nothing
      else
       print *,"rz_flag invalid in faceinit"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dir.eq.0) then
       ii=1
      else if (dir.eq.1) then
       jj=1
      else if ((dir.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dir invalid faceprocess"
       stop
      endif
 
      call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_dest,dir) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalfMAC,dir)

       if (dir.eq.0) then
        inormal=i
       else if (dir.eq.1) then
        inormal=j
       else if ((dir.eq.2).and.(SDIM.eq.3)) then
        inormal=k
       else
        print *,"dir invalid"
        stop
       endif

       if ((inormal-1.ge.fablo(dir+1)-ngrow_source).and. &
           (inormal-1.le.fabhi(dir+1)+ngrow_source)) then
        left_face_ok=1
       else if (inormal-1.eq.fablo(dir+1)-ngrow_source-1) then
        left_face_ok=0
       else
        print *,"inormal invalid"
        stop
       endif

       if ((inormal.ge.fablo(dir+1)-ngrow_source).and. &
           (inormal.le.fabhi(dir+1)+ngrow_source)) then
        right_face_ok=1
       else if (inormal.eq.fabhi(dir+1)+ngrow_source+1) then
        right_face_ok=0
       else
        print *,"inormal invalid"
        stop
       endif

       at_RZ_face=0
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (dir.eq.0) then
         if (abs(xstenMAC(0,1)).le.EPS2*dx(1)) then
          at_RZ_face=1
         endif
         if (inormal.eq.i) then
          if (i.le.0) then
           at_RZ_face=1
          endif
         else
          print *,"inormal invalid"
          stop
         endif
        endif
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        if (dir.eq.0) then
         if (abs(xstenMAC(0,1)).le.EPS2*dx(1)) then
          at_RZ_face=1
         endif
         if (inormal.eq.i) then
          if (xstenMAC(0,1).lt.zero) then
           at_RZ_face=1
          endif
         else
          print *,"inormal invalid"
          stop
         endif
        endif
       else
        print *,"levelrz invalid"
        stop
       endif

       if (at_RZ_face.eq.0) then

        left_total=zero
        right_total=zero
         ! (im,dir,side)
        do im=1,num_materials
         side=2
         call abs_array_index3(im,dir+1,side, &
          num_materials,SDIM,2,iface_left)
         side=1
         call abs_array_index3(im,dir+1,side, &
          num_materials,SDIM,2,iface_right)

         if (left_face_ok.eq.1) then

          leftface(im)=facefab(D_DECL(i-ii,j-jj,k-kk),iface_left)
          if (right_face_ok.eq.1) then
           ! do nothing
          else if (right_face_ok.eq.0) then
           rightface(im)=leftface(im)
          else
           print *,"right_face_ok invalid"
           stop
          endif
           
         else if (left_face_ok.eq.0) then
          ! do nothing
         else
          print *,"left_face_ok invalid"
          stop
         endif

         if (right_face_ok.eq.1) then

          rightface(im)=facefab(D_DECL(i,j,k),iface_right)
          if (left_face_ok.eq.1) then
           ! do nothing
          else if (left_face_ok.eq.0) then
           leftface(im)=rightface(im)
          else
           print *,"left_face_ok invalid"
           stop
          endif

         else if (right_face_ok.eq.0) then
          ! do nothing
         else
          print *,"right_face_ok invalid"
          stop
         endif

         frac_left(im)=leftface(im)
         frac_right(im)=rightface(im)

         if ((tessellate.eq.1).or. &
             (tessellate.eq.3).or. &
             (is_rigid(im).eq.0)) then
          left_total=left_total+frac_left(im)
          right_total=right_total+frac_right(im)
         else if ((tessellate.eq.0).and.(is_rigid(im).eq.1)) then
          ! do nothing
         else 
          print *,"tessellate or is_rigid invalid MOF_REDIST_3D.F90"
          stop
         endif 

        enddo ! im=1..num_materials

        do im=1,num_materials
         if (left_total.gt.zero) then
          frac_left(im)=frac_left(im)/left_total
         else if (left_total.eq.zero) then
          frac_left(im)=zero
         else
          print *,"left_total invalid: ",left_total
          stop
         endif
         if (right_total.gt.zero) then
          frac_right(im)=frac_right(im)/right_total
         else if (right_total.eq.zero) then
          frac_right(im)=zero
         else
          print *,"right_total invalid: ",right_total
          stop
         endif
        enddo ! im=1,num_materials

        call CISBOX(xsten_left,nhalf, &
         xlo,dx,i-ii,j-jj,k-kk, &
         bfact,level, &
         volcell_left,cencell_left,SDIM) 

        call CISBOX(xsten_right,nhalf, &
         xlo,dx,i,j,k, &
         bfact,level, &
         volcell_right,cencell_right,SDIM) 

        do im = 1,num_materials

         if ((tessellate.eq.1).or. &
             (tessellate.eq.3).or. &
             (is_rigid(im).eq.0)) then

          if ((frac_left(im).gt.one+EPS3).or. &
              (frac_left(im).lt.zero).or. &
              (frac_right(im).gt.one+EPS3).or. &
              (frac_right(im).lt.zero)) then
           print *,"frac_left or frac_right out of range"
           print *,"im,frac_left ",im,frac_left(im)
           print *,"im,frac_right ",im,frac_right(im)
           stop
          endif
          if (frac_left(im).le.EPS_3_2) then
           frac_left(im) = zero
          else if (frac_left(im).ge.one-EPS_3_2) then
           frac_left(im) = one
          endif
          if (frac_right(im).lt.EPS_3_2) then
           frac_right(im) = zero
          else if (frac_right(im).gt.one-EPS_3_2)then
           frac_right(im) = one
          endif

         else if ((tessellate.eq.0).and.(is_rigid(im).eq.1)) then
          ! do nothing
         else
          print *,"tessellate or is_rigid invalid MOF_REDIST_3D.F90"
          stop
         endif 

        enddo ! im=1..num_materials

        if ((left_face_ok.eq.1).and. &  ! cell to left of face is in FAB
            (right_face_ok.eq.1)) then  ! cell to right of face is in FAB

         do im=1,num_materials*ngeom_recon
          mofdata_left(im)=vofrecon(D_DECL(i-ii,j-jj,k-kk),im)
          mofdata_right(im)=vofrecon(D_DECL(i,j,k),im)
         enddo

         do im=1,num_materials
          vofcomp=(im-1)*ngeom_recon+1
          if ((tessellate.eq.1).or. &
              (tessellate.eq.3).or. &
              (is_rigid(im).eq.0)) then
           ! do nothing
          else if ((tessellate.eq.0).and. &
                   (is_rigid(im).eq.1)) then
           do dir2=1,ngeom_recon
            mofdata_left(vofcomp+dir2-1)=zero 
            mofdata_right(vofcomp+dir2-1)=zero 
           enddo
          else
           print *,"tessellate or is_rigid invalid MOF_REDIST_3D.F90"
           stop
          endif

          do dir2=1,SDIM
           multi_cen_left(dir2,im)=mofdata_left(vofcomp+dir2)
           multi_cen_right(dir2,im)=mofdata_right(vofcomp+dir2)
          enddo

         enddo ! im=1..num_materials

         if (tessellate.eq.0) then
          local_tessellate=1 ! is_rigid data has been zeroed out.
         else if ((tessellate.eq.1).or. &
                  (tessellate.eq.3)) then
          local_tessellate=tessellate
         else
          print *,"tessellate invalid"
          stop
         endif

          ! if local_tessellate==1 or 3:
          !  a) call multi_get_volume_tessellate
          !  b) local_tessellate_in=2
          !  c) is_rigid_local=0 for all materials
          ! else if local_tessellate==0
          !  a) local_tessellate_in=0
          ! 
         call multi_get_area_pairs( &
           local_tessellate, & ! =1 or 3
           bfact, &
           dx, &
           xsten_right, &
           xsten_left, &
           nhalf, &
           mofdata_right, &
           mofdata_left, &
           dir+1, &
           frac_pair, & ! left,right
           SDIM, &
           geom_xtetlist(1,1,1,tid+1), &
           nmax, &
           geom_xtetlist_old(1,1,1,tid+1), &
           nmax, &
           nmax)

        else if ((left_face_ok.eq.0).or. &
                 (right_face_ok.eq.0)) then

         do im=1,num_materials
          vofcomp=(im-1)*ngeom_recon+1
          do dir2=1,SDIM
           multi_cen_left(dir2,im)=zero
           multi_cen_right(dir2,im)=zero
          enddo
         enddo ! im

         do ml = 1, num_materials
          do mr = 1, num_materials
           frac_pair(ml,mr)=zero
          enddo
         enddo
         do ml=1,num_materials
          frac_pair(ml,ml)=frac_left(ml)
         enddo

        else
         print *,"left or right face ok invalid"
         stop
        endif

        vol_total=zero
        do im=1,num_materials
        do im_opp=1,num_materials
         vol_total=vol_total+frac_pair(im,im_opp)
        enddo
        enddo

        do im=1,num_materials
        do im_opp=1,num_materials
         if (vol_total.gt.zero) then
          frac_pair(im,im_opp)=frac_pair(im,im_opp)/vol_total
         else if (vol_total.eq.zero) then
          frac_pair(im,im_opp)=zero
         else
          print *,"vol_total invalid: ",vol_total
          stop
         endif
        enddo !im_opp=1,num_materials
        enddo !im=1,num_materials

        delta=xsten_right(0,dir+1)-xsten_left(0,dir+1)
        if (delta.gt.zero) then
         ! do nothing
        else
         print *,"delta invalid faceprocess: ",delta
         stop
        endif 

        do ml = 1, num_materials

         if ((tessellate.eq.1).or. &
             (tessellate.eq.3).or. &
             (is_rigid(ml).eq.0)) then

          do mr = 1, num_materials

           if ((tessellate.eq.1).or. &
               (tessellate.eq.3).or. &
               (is_rigid(mr).eq.0)) then

            if ((frac_pair(ml,mr).lt.-EPS1).or. &
                (frac_pair(ml,mr).gt.one+EPS1)) then
             print *,"frac_pair invalid: ml,mr,frac_pair ", &
               ml,mr,frac_pair(ml,mr)
             stop
            endif
  
            if (frac_pair(ml,mr).lt.EPS_3_2) then
             frac_pair(ml,mr)=zero
             dist_pair(ml,mr)=delta
            else
             dpair=zero
             do dir2=1,SDIM
              dpair=dpair+ &
               (xsten_right(0,dir2)+multi_cen_right(dir2,mr)- &
                xsten_left(0,dir2)-multi_cen_left(dir2,ml))**2
             enddo
             dpair=sqrt(dpair)
             dist_pair(ml,mr)=dpair
            endif
    
            if (dir.eq.1) then ! theta direction
   
             if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
                 (levelrz.eq.COORDSYS_RZ)) then !XY or RZ
              RR=one 
             else if (levelrz.eq.COORDSYS_CYLINDRICAL) then ! R-theta
              RR=half*(xsten_right(0,1)+xsten_left(0,1))
             else
              print *,"levelrz invalid"
              stop
             endif
             dist_pair(ml,mr)=dist_pair(ml,mr)*RR

            else if ((dir.eq.0).or.(dir.eq.SDIM-1)) then
             ! do nothing
            else
             print *,"dir invalid face process"
             stop
            endif ! dir==1 ?

           else if ((tessellate.eq.0).and.(is_rigid(mr).eq.1)) then
            ! do nothing
           else
            print *,"tessellate or is_rigid invalid MOF_REDIST_3D.F90"
            stop
           endif 

          enddo ! mr=1..num_materials

         else if ((tessellate.eq.0).and.(is_rigid(ml).eq.1)) then
          ! do nothing
         else
          print *,"tessellate or is_rigid invalid MOF_REDIST_3D.F90"
          stop
         endif 

        enddo ! ml=1..num_materials

       else if (at_RZ_face.eq.1) then

        do ml = 1, num_materials
         do mr = 1, num_materials
          frac_pair(ml,mr)=zero
          dist_pair(ml,mr)=zero
         enddo
        enddo

       else
        print *,"at_RZ_face invalid"
        stop
       endif

       iface_left=0
       do ml = 1, num_materials
        do mr = 1, num_materials
          !(ml,mr,2) 
         iface_left=iface_left+1
         dstfab(D_DECL(i,j,k),iface_left)=frac_pair(ml,mr)
         iface_left=iface_left+1
         dstfab(D_DECL(i,j,k),iface_left)=dist_pair(ml,mr)
        enddo ! mr
       enddo ! ml
       if (iface_left.ne.nface_dst) then
        print *,"iface_left invalid"
        stop
       endif

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine fort_faceprocess

      end module mof_redist_cpp_module

#undef STANDALONE

