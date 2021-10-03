#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 1

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "MOF_REDIST_F.H"

#define nsum 64
#define nsum2 32
#define VISCINVTOL (1.0D-8)
#define CURVWT (1.0D-3)

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

        module mof_redist_module
        use probcommon_module

        contains

      subroutine update_closest( &
        xsten_accept,xsten_donate,nhalf, &
        dx,xlo,bfact,level,fablo, &
        mofdata,nmat,nstar, &
        idon,jdon,kdon, & ! donate index
        i1,j1,k1, &  ! accept index: idon+i1,jdon+j1,kdon+k1
        newLS, &
        touch_hold, &
        minLS, &
        maxLS, &
        donateflag, &
        time)
      use global_utility_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: nhalf
      INTEGER_T, intent(in) :: bfact,nmat,nstar
      INTEGER_T, intent(in) :: idon,jdon,kdon
      INTEGER_T, intent(in) :: i1,j1,k1
      INTEGER_T, intent(in) :: fablo(SDIM)
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: xsten_accept(-nhalf:nhalf,SDIM)
      REAL_T, intent(in) :: xsten_donate(-nhalf:nhalf,SDIM)
      REAL_T :: xsten_vert(-nhalf:nhalf,SDIM)
      REAL_T :: xdonate_vert(SDIM)
      REAL_T :: xdonate_point(SDIM)
      REAL_T :: xaccept_point(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: mofdata(nmat*ngeom_recon)
      REAL_T, intent(inout) :: newLS(nmat*(1+SDIM))
      INTEGER_T, intent(inout) :: touch_hold(nmat)
      REAL_T, intent(inout) :: minLS(nmat)
      REAL_T, intent(inout) :: maxLS(nmat)
      INTEGER_T, intent(in) :: donateflag(nmat+1+nstar)
      INTEGER_T :: center_stencil
      INTEGER_T :: im0_center
      INTEGER_T :: imslope

      INTEGER_T nstar_test
      INTEGER_T istar_array(3)
      INTEGER_T istar,i2,j2,k2,donateIND
      INTEGER_T klosten,khisten
      INTEGER_T dir
      REAL_T LSslope(SDIM)
      INTEGER_T n_im
      INTEGER_T im_test(6)
      INTEGER_T i_DEB_DIST
      INTEGER_T j_DEB_DIST
      INTEGER_T k_DEB_DIST

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
       call get_primary_slope( &
        bfact,dx,xsten_donate,nhalf, &
        mofdata, &
        LSslope,imslope,nmat,SDIM) 
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
      im0_center=donateflag(nmat+1+istar)

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
      do i2=-1,1
      do j2=-1,1
      do k2=klosten,khisten
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
        donateIND=donateflag(nmat+1+istar)

        if (donateIND.eq.0) then
         ! do nothing (corner is on FAB boundary or corner is
         ! occupied by flotsam that should be ignored)
        else if ((donateIND.ge.1).and.(donateIND.le.nmat)) then

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
         call compare_distance( &
          bfact,dx,xsten_donate,nhalf, &
          nmat, &
          xaccept_point, &
          xdonate_vert, &
          newLS, &
          touch_hold, &
          minLS, &
          maxLS, &
          im_test,n_im,LSslope, &
          imslope, &
          im0_center, &
          SDIM, &
          center_stencil, &
          donateflag,6)

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

       ! only uses donateflag(1..nmat+1)
      call multi_get_distance( &
        bfact,dx,xsten_donate,nhalf,xaccept_point, &
        mofdata, &
        newLS, &
        touch_hold, &
        minLS, &
        maxLS, &
        nmat,SDIM, &
        center_stencil, &
        donateflag)

      return
      end subroutine update_closest

      end module mof_redist_module

#if (STANDALONE==1)
      module mof_redist_cpp_module
      contains
#endif

       ! prior to calling this routine, copy LS_new normal information
       ! to LS_NRM_FD.
       ! called from: NavierStokes::build_NRM_FD_MF (NavierStokes.cpp)
      subroutine FORT_FD_NORMAL( &
       level, &
       finest_level, &
       LS_new, &
       DIMS(LS_new), &
       LS_NRM_FD, &
       DIMS(LS_NRM_FD), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       nmat)
      use global_utility_module
      use probcommon_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: DIMDEC(LS_new)
      INTEGER_T, intent(in) :: DIMDEC(LS_NRM_FD)
      REAL_T, intent(in) :: LS_new(DIMV(LS_new),nmat*(1+SDIM))
      REAL_T, intent(out) :: LS_NRM_FD(DIMV(LS_NRM_FD),nmat*SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T nhalf
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T dir
      INTEGER_T im
      REAL_T ls_stencil(D_DECL(-1:1,-1:1,-1:1),nmat)
      REAL_T lsnormal(nmat,SDIM)
      INTEGER_T lsnormal_valid(nmat)
      REAL_T ls_intercept(nmat)
      REAL_T dxmaxLS
      INTEGER_T k1lo,k1hi
      INTEGER_T i,j,k
      INTEGER_T i1,j1,k1
      INTEGER_T dcomp
      REAL_T local_LS(nmat)
      INTEGER_T im_primary,im_secondary,triple_point_flag

      nhalf=3 

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid141"
       stop
      endif
      if ((level.le.finest_level).and.(level.ge.0)) then
       ! do nothing
      else
       print *,"level invalid in FORT_FD_NORMAL"
       stop
      endif
      call checkbound(fablo,fabhi,DIMS(LS_new),1,-1,2871)
      call checkbound(fablo,fabhi,DIMS(LS_NRM_FD),0,-1,2872)
      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid FORT_FD_NORMAL"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.eq.SDIM+1) then
       ! do nothing
      else
       print *,"ngeom_raw invalid FORT_FD_NORMAL"
       print *,"ngeom_raw=",ngeom_raw
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
       print *,"levelrz invalid in FORT_FD_NORMAL"
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

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       do i1=-1,1
       do j1=-1,1
       do k1=k1lo,k1hi
       do im=1,nmat
        ls_stencil(D_DECL(i1,j1,k1),im)= &
                LS_new(D_DECL(i+i1,j+j1,k+k1),im)
       enddo
       enddo
       enddo
       enddo

       do im=1,nmat
        local_LS(im)=ls_stencil(D_DECL(0,0,0),im)
       enddo
       call get_primary_material(local_LS,nmat,im_primary)
       call get_secondary_material(local_LS,nmat,im_primary,im_secondary)
       triple_point_flag=0
       do im=1,nmat
        if ((im.ne.im_primary).and.(im.ne.im_secondary)) then
         if (abs(local_LS(im)).le.dxmaxLS) then
          triple_point_flag=1
         else if (abs(local_LS(im)).ge.dxmaxLS) then
          ! do nothing
         else
          print *,"local_LS bust"
          stop
         endif
        else if ((im.ge.1).and.(im.le.nmat)) then
         ! do nothing
        else
         print *,"im bust"
         stop
        endif
       enddo ! im=1..nmat
      
       if (triple_point_flag.eq.0) then 
        do im=1,nmat
         if (is_rigid(nmat,im).eq.0) then
          if (abs(local_LS(im)).le.two*dxmaxLS) then
           if ((im.eq.im_primary).or.(im.eq.im_secondary)) then
            call find_cut_geom_slope_CLSVOF( &
             ls_stencil, & ! (-1,1)^3,nmat
             lsnormal, &  ! (nmat,sdim)
             lsnormal_valid, &  ! nmat
             ls_intercept, & ! nmat
             bfact,dx, &
             xsten,nhalf, &
             im, &
             dxmaxLS, &
             nmat,SDIM)

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
           else if ((im.ge.1).and.(im.le.nmat)) then
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
         else if (is_rigid(nmat,im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif
        enddo ! im=1..nmat
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
      end subroutine FORT_FD_NORMAL

      subroutine FORT_FD_NODE_NORMAL( &
       level, &
       finest_level, &
       LS_new, &
       DIMS(LS_new), &
       LS_NRM_FD, &
       DIMS(LS_NRM_FD), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       nmat, &
       nten, &
       n_normal, &
       ngrow_dest)
      use global_utility_module
      use probcommon_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: n_normal
      INTEGER_T, intent(in) :: ngrow_dest
      INTEGER_T, intent(in) :: DIMDEC(LS_new)
      INTEGER_T, intent(in) :: DIMDEC(LS_NRM_FD)
      REAL_T, intent(in) :: LS_new(DIMV(LS_new),nmat*(1+SDIM))
      REAL_T, intent(out) :: LS_NRM_FD(DIMV(LS_NRM_FD),n_normal)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      REAL_T xsten_nd(-3:3,SDIM)
      INTEGER_T nten_test
      INTEGER_T n_normal_test
      INTEGER_T nhalf
      INTEGER_T k1hi
      INTEGER_T i,j,k
      INTEGER_T i1,j1,k1
      INTEGER_T im,im1,im2
      INTEGER_T iten
      INTEGER_T ibase
      INTEGER_T dir
      REAL_T local_normal(SDIM)
      REAL_T local_LS
      REAL_T local_mag
      REAL_T sign_nm
      REAL_T xplus,xminus,RR


      nhalf=3 

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid141"
       stop
      endif
      if ((level.le.finest_level).and.(level.ge.0)) then
       ! do nothing
      else
       print *,"level invalid in FORT_FD_NODE_NORMAL"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid fd_node_normal nten nten_test ",nten,nten_test
       stop
      endif
      n_normal_test=(SDIM+1)*(nten+nmat)
      if (n_normal_test.ne.n_normal) then
       print *,"n_normal invalid fd_node_normal n_normal ",n_normal
       stop
      endif
 
      if (ngrow_dest.ge.0) then
       ! do nothing
      else
       print *,"ngrow_dest invalid"
       stop
      endif
      
      call checkbound(fablo,fabhi,DIMS(LS_new),ngrow_dest+1,-1,2873)
      call checkbound(fablo,fabhi,DIMS(LS_NRM_FD),ngrow_dest+1,-1,2874)
      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid FORT_FD_NODE_NORMAL"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.eq.SDIM+1) then
       ! do nothing
      else
       print *,"ngeom_raw invalid FORT_FD_NODE_NORMAL"
       print *,"ngeom_raw=",ngeom_raw
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
       print *,"levelrz invalid in FORT_FD_NODE_NORMAL"
       stop
      endif

      k1hi=0
      if (SDIM.eq.3) then
       k1hi=1
      endif

      call growntileboxNODE(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_dest)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridstenND_level(xsten_nd,i,j,k,level,nhalf)

       do im=1,nmat+nten
        do dir=1,SDIM
         local_normal(dir)=zero
        enddo
        local_mag=zero

        do dir=1,SDIM
         do i1=0,1
         do j1=0,1
         do k1=0,k1hi
          if ((im.ge.1).and.(im.le.nmat)) then
           local_LS=LS_new(D_DECL(i+i1-1,j+j1-1,k+k1-1),im) 
          else if ((im.ge.nmat+1).and.(im.le.nmat+nten)) then
           iten=im-nmat
           call get_inverse_iten(im1,im2,iten,nmat)
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
         if (levelrz.eq.0) then
          RR=one
         else if (levelrz.eq.1) then
          RR=one
          if ((dir.eq.1).and.(xminus.lt.zero)) then
           local_normal(dir)=zero
          endif
         else if (levelrz.eq.3) then
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
         LS_NRM_FD(D_DECL(i,j,k),ibase+dir)=local_normal(dir)
        enddo 
        LS_NRM_FD(D_DECL(i,j,k),ibase+SDIM+1)=local_mag
       enddo ! im=1..nmat+nten

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine FORT_FD_NODE_NORMAL

      subroutine FORT_NODE_TO_CELL( &
       level, &
       finest_level, &
       height_function_flag, &  ! 1=> use height function 0 => use FD
       F_new, &  !F_new(i,j,k,im)  im=1..nmat
       DIMS(F_new), &
       LS_new, &
       DIMS(LS_new), &
       LS_NRM_FD, &
       DIMS(LS_NRM_FD), &
       CURV_CELL, &
       DIMS(CURV_CELL), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       nmat, &
       nten, &
       n_normal, &
       ngrow_nrm)
      use global_utility_module
      use probcommon_module
      use MOF_routines_module
      use height_method_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: height_function_flag
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: n_normal
      INTEGER_T, intent(in) :: ngrow_nrm
      INTEGER_T, intent(in) :: DIMDEC(LS_new)
      INTEGER_T, intent(in) :: DIMDEC(F_new)
      INTEGER_T, intent(in) :: DIMDEC(LS_NRM_FD)
      INTEGER_T, intent(in) :: DIMDEC(CURV_CELL)
      REAL_T, intent(in) :: F_new(DIMV(F_new),nmat)
      REAL_T, intent(in) :: LS_new(DIMV(LS_new),nmat*(1+SDIM))
      REAL_T, intent(in) :: LS_NRM_FD(DIMV(LS_NRM_FD),n_normal)
       ! first nmat+nten components are curvature
       ! second nmat+nten components are status (0=bad 1=good)
      REAL_T, intent(out) :: CURV_CELL(DIMV(CURV_CELL),2*(nmat+nten))
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nten_test
      INTEGER_T n_normal_test
      INTEGER_T nhalf
      INTEGER_T k1hi
      INTEGER_T i,j,k
      INTEGER_T i1,j1,k1
      INTEGER_T im
      INTEGER_T im_local
      INTEGER_T im1,im2
      INTEGER_T im_primary,im_secondary
      INTEGER_T iten
      INTEGER_T ibase
      INTEGER_T dir
      REAL_T local_curv(SDIM)
      REAL_T local_normal(SDIM)
      REAL_T local_mag
      REAL_T sign_nm
      REAL_T xplus,xminus,xmiddle,RR
      REAL_T denom_factor
      REAL_T total_curv
      INTEGER_T local_status
      REAL_T local_LS(nmat)
      REAL_T kappa(nmat+1)  !nmat+nten
      REAL_T F_local
      REAL_T LS_local
      INTEGER_T curv_valid
      REAL_T, dimension(nmat,-3:3,-3:3) :: vf_curv
      REAL_T, dimension(nmat,-3:3,-3:3) :: ls_curv

      nhalf=3 

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid141"
       stop
      endif
      if ((level.le.finest_level).and.(level.ge.0)) then
       ! do nothing
      else
       print *,"level invalid in FORT_NODE_TO_CELL"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid node_to_cell nten nten_test ",nten,nten_test
       stop
      endif
      n_normal_test=(SDIM+1)*(nten+nmat)
      if (n_normal_test.ne.n_normal) then
       print *,"n_normal invalid node_to_cell n_normal ",n_normal
       stop
      endif
 
      if (ngrow_nrm.ge.0) then
       ! do nothing
      else
       print *,"ngrow_nrm invalid"
       stop
      endif
      
      call checkbound(fablo,fabhi,DIMS(F_new),ngrow_nrm+1,-1,2875)
      call checkbound(fablo,fabhi,DIMS(LS_new),ngrow_nrm+1,-1,2875)
      call checkbound(fablo,fabhi,DIMS(LS_NRM_FD),ngrow_nrm+1,-1,2876)
      call checkbound(fablo,fabhi,DIMS(CURV_CELL),ngrow_nrm,-1,2877)
      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.eq.2*SDIM+3) then
       ! do nothing
      else
       print *,"ngeom_recon invalid FORT_NODE_TO_CELL"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.eq.SDIM+1) then
       ! do nothing
      else
       print *,"ngeom_raw invalid FORT_NODE_TO_CELL"
       print *,"ngeom_raw=",ngeom_raw
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
       print *,"levelrz invalid in FORT_NODE_TO_CELL"
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
        growlo,growhi,ngrow_nrm)


      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       if (height_function_flag.eq.0) then

        do im=1,nmat+nten

         do im_local=1,nmat
          local_LS(im_local)=LS_new(D_DECL(i,j,k),im_local)
         enddo
         call get_primary_material(local_LS,nmat,im_primary)
         call get_secondary_material(local_LS,nmat,im_primary,im_secondary)

         local_status=1
         if ((im.ge.1).and.(im.le.nmat)) then
          if ((im.ne.im_primary).and. &
              (im.ne.im_secondary)) then
           local_status=0
          endif
          if (is_rigid(nmat,im).eq.1) then
           ! do nothing
          else if (is_rigid(nmat,im).eq.0) then
           if (is_rigid(nmat,im_primary).eq.0) then
            ! do nothing
           else if (is_rigid(nmat,im_primary).eq.1) then
            local_status=0
           else
            print *,"is_rigid(nmat,im_primary) invalid"
            stop
           endif
          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif
         else if ((im.ge.nmat+1).and.(im.le.nmat+nten)) then
          iten=im-nmat
          call get_inverse_iten(im1,im2,iten,nmat)
          if ((im1.ne.im_primary).and. &
              (im1.ne.im_secondary)) then
           local_status=0
          endif
          if ((im2.ne.im_primary).and. &
              (im2.ne.im_secondary)) then
           local_status=0
          endif
         else
          print *,"im invalid 112"
          stop
         endif
        
         do dir=1,SDIM
          local_curv(dir)=zero
         enddo

         ibase=(im-1)*(SDIM+1)

         do dir=1,SDIM
          do i1=0,1
          do j1=0,1
          do k1=0,k1hi
           local_normal(dir)=LS_NRM_FD(D_DECL(i+i1,j+j1,k+k1),ibase+dir)
           local_mag=LS_NRM_FD(D_DECL(i+i1,j+j1,k+k1),ibase+SDIM+1)
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
           if (levelrz.eq.0) then
            ! do nothing
           else if (levelrz.eq.1) then
            if (dir.eq.1) then
             RR=xsten(2*i1-1,dir)
             if (RR.lt.zero) then
              RR=zero
             endif
            endif
           else if (levelrz.eq.3) then
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
          if (levelrz.eq.0) then
           RR=one
          else if (levelrz.eq.1) then
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
          else if (levelrz.eq.3) then
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
         CURV_CELL(D_DECL(i,j,k),im)=total_curv
         CURV_CELL(D_DECL(i,j,k),nmat+nten+im)=local_status
        enddo ! im=1..nmat+nten

       else if (height_function_flag.eq.1) then

        if (SDIM.eq.2) then
         ! do nothing
        else
         print *,"this is 2d code only for now"
         stop
        endif
        k1=0

        if (nmat.eq.2) then

         do im=1,nmat+1
          kappa(im)=0.0d0
         enddo
         curv_valid=1

         do i1=-3,3
         do j1=-3,3
          if ((i+i1.ge.fablo(1)).and. &
              (i+i1.le.fabhi(1)).and. &
              (j+j1.le.fabhi(2)).and. &
              (j+j1.ge.fablo(2))) then
           F_local=F_new(D_DECL(i+i1,j+j1,k+k1),1)
           LS_local=F_local-0.5d0
          else
#if (STANDALONE==1)
           curv_valid=0
           F_local=1.0d0
           LS_local=F_local-0.5d0
#else
           print *,"only use height function in prototype code"
           stop
#endif

          endif
          vf_curv(1,i1,j1)=F_local
          ls_curv(1,i1,j1)=F_local-0.5d0
          vf_curv(2,i1,j1)=one-F_local
          ls_curv(2,i1,j1)=0.5d0-F_local
         enddo
         enddo

         if (curv_valid.eq.1) then

          curv_valid=0

          if (ls_curv(1,0,0).ge.0.0d0) then

           if ((ls_curv(1,1,0).lt.0.0d0).and. &
               (ls_curv(1,-1,0).lt.0.0d0)) then
            curv_valid=-1
           else if ((ls_curv(1,1,0).ge.0.0d0).and. &
                    (ls_curv(1,-1,0).ge.0.0d0)) then
            ! do nothing
           else if ((ls_curv(1,1,0).lt.0.0d0).and. &
                    (ls_curv(1,-1,0).ge.0.0d0)) then
            curv_valid=1
           else if ((ls_curv(1,1,0).ge.0.0d0).and. &
                    (ls_curv(1,-1,0).lt.0.0d0)) then
            curv_valid=1
           else
            print *,"ls_curv NaN"
            stop
           endif

           if (curv_valid.ge.0) then

            if ((ls_curv(1,0,1).lt.0.0d0).and. &
                (ls_curv(1,0,-1).lt.0.0d0)) then
             curv_valid=-1
            else if ((ls_curv(1,0,1).ge.0.0d0).and. &
                     (ls_curv(1,0,-1).ge.0.0d0)) then
             ! do nothing
            else if ((ls_curv(1,0,1).lt.0.0d0).and. &
                     (ls_curv(1,0,-1).ge.0.0d0)) then
             curv_valid=1
            else if ((ls_curv(1,0,1).ge.0.0d0).and. &
                     (ls_curv(1,0,-1).lt.0.0d0)) then
             curv_valid=1
            else
             print *,"ls_curv NaN"
             stop
            endif

           else if (curv_valid.eq.-1) then
            ! do nothing
           else
            print *,"curv_valid invalid"
            stop
           endif

          else if (ls_curv(1,0,0).lt.0.0d0) then

           if ((ls_curv(1,1,0).ge.0.0d0).and. &
               (ls_curv(1,-1,0).ge.0.0d0)) then
            curv_valid=-1
           else if ((ls_curv(1,1,0).lt.0.0d0).and. &
                    (ls_curv(1,-1,0).lt.0.0d0)) then
            ! do nothing
           else if ((ls_curv(1,1,0).lt.0.0d0).and. &
                    (ls_curv(1,-1,0).ge.0.0d0)) then
            curv_valid=1
           else if ((ls_curv(1,1,0).ge.0.0d0).and. &
                    (ls_curv(1,-1,0).lt.0.0d0)) then
            curv_valid=1
           else
            print *,"ls_curv NaN"
            stop
           endif

           if (curv_valid.ge.0) then

            if ((ls_curv(1,0,1).ge.0.0d0).and. &
                (ls_curv(1,0,-1).ge.0.0d0)) then
             curv_valid=-1
            else if ((ls_curv(1,0,1).lt.0.0d0).and. &
                     (ls_curv(1,0,-1).lt.0.0d0)) then
             ! do nothing
            else if ((ls_curv(1,0,1).lt.0.0d0).and. &
                     (ls_curv(1,0,-1).ge.0.0d0)) then
             curv_valid=1
            else if ((ls_curv(1,0,1).ge.0.0d0).and. &
                     (ls_curv(1,0,-1).lt.0.0d0)) then
             curv_valid=1
            else
             print *,"ls_curv NaN"
             stop
            endif

           else if (curv_valid.eq.-1) then
            ! do nothing
           else
            print *,"curv_valid invalid"
            stop
           endif

          else
           print *,"ls_curv(1,0,0) is NaN"
           stop
          endif

          if (curv_valid.eq.1) then

           call get_curvature_heightf(nmat,ls_curv,vf_curv,dx(1), &
              kappa,curv_valid)

           if (curv_valid.eq.1) then
            kappa(nmat+1)=kappa(1)
           endif

          else if ((curv_valid.eq.-1).or.(curv_valid.eq.0)) then
           ! do nothing
          else
           print *,"curv_valid invalid"
           stop
          endif

         else if (curv_valid.eq.0) then
          ! do nothing
         else
          print *,"curv_valid invalid"
          stop
         endif

         if (curv_valid.eq.-1) then
          curv_valid=0
         endif

         do im=1,nmat+1
          CURV_CELL(D_DECL(i,j,k),im)=kappa(im)
         enddo
          ! status=1 good status=0 bad
         do im=1,nmat+1
          CURV_CELL(D_DECL(i,j,k),nmat+nten+im)=curv_valid
         enddo

        else
         print *,"only two materials in 2d supported"
         stop
        endif

       else
        print *,"height_function_flag invalid"
        stop
       endif

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine FORT_NODE_TO_CELL



        ! vofrecon=vof,ref centroid,order,slope,intercept
        ! newfab has nmat*(sdim+1) components
        !
      subroutine FORT_LEVELSTRIP( &
         keep_all_interfaces, &
         nprocessed, &
         minLS, &
         maxLS, &
         max_problen, &
         level, &
         finest_level, &
         truncate_volume_fractions, &
         latent_heat, &
         maskfab,DIMS(maskfab), &
         facepairX,DIMS(facepairX), &
         facepairY,DIMS(facepairY), &
         facepairZ,DIMS(facepairZ), &
         facefab,DIMS(facefab), &
         facetest,DIMS(facetest), &
         stenfab,DIMS(stenfab), &
         vofrecon,DIMS(vofrecon), &
         origdist,DIMS(origdist), &
         newfab,DIMS(newfab), &
         touchfab,DIMS(touchfab), &
         crsetouch,DIMS(crsetouch), &
         crsedist,DIMS(crsedist), &
         tilelo,tilehi, &
         fablo,fabhi, &
         bfact, &
         bc, &
         rz_flag, &
         xlo,dx, &
         time, &
         ngrow_distance, &
         nmat,nten, &
         nstar, &
         nface, &
         nface_dst)

      use global_utility_module
      use probcommon_module
      use MOF_routines_module
      use mof_redist_module

      IMPLICIT NONE

      INTEGER_T, intent(inout) :: nprocessed
      INTEGER_T, intent(in) :: keep_all_interfaces 
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nstar
      INTEGER_T, intent(in) :: nface
      INTEGER_T, intent(in) :: nface_dst
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: ngrow_distance
      REAL_T, intent(inout) :: minLS(nmat)
      REAL_T, intent(inout) :: maxLS(nmat)
      REAL_T, intent(in) :: max_problen
      INTEGER_T, intent(in) :: truncate_volume_fractions(nmat)
      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: DIMDEC(maskfab)
      INTEGER_T, intent(in) :: DIMDEC(facepairX)
      INTEGER_T, intent(in) :: DIMDEC(facepairY)
      INTEGER_T, intent(in) :: DIMDEC(facepairZ)
      INTEGER_T, intent(in) :: DIMDEC(facefab)
      INTEGER_T, intent(in) :: DIMDEC(facetest)
      INTEGER_T, intent(in) :: DIMDEC(stenfab)
      INTEGER_T, intent(in) :: DIMDEC(vofrecon)
      INTEGER_T, intent(in) :: DIMDEC(origdist)
      INTEGER_T, intent(in) :: DIMDEC(newfab)
      INTEGER_T, intent(in) :: DIMDEC(touchfab)
      INTEGER_T, intent(in) :: DIMDEC(crsetouch)
      INTEGER_T, intent(in) :: DIMDEC(crsedist)

      REAL_T, intent(in) :: maskfab(DIMV(maskfab),4)
      REAL_T, intent(in) :: facepairX(DIMV(facepairX),nface_dst)
      REAL_T, intent(in) :: facepairY(DIMV(facepairY),nface_dst)
      REAL_T, intent(in) :: facepairZ(DIMV(facepairZ),nface_dst)
      REAL_T, intent(in) :: facefab(DIMV(facefab),nface)
      REAL_T, intent(in) :: facetest(DIMV(facetest),nmat*SDIM)
      REAL_T, intent(in) :: stenfab(DIMV(stenfab),nstar)

      REAL_T, intent(in) :: vofrecon(DIMV(vofrecon),nmat*ngeom_recon)
      REAL_T, intent(in) :: origdist(DIMV(origdist),nmat*(1+SDIM))
      REAL_T, intent(inout) :: newfab(DIMV(newfab),nmat*(1+SDIM))
      REAL_T, intent(inout) :: touchfab(DIMV(touchfab),nmat)
      REAL_T, intent(in) :: crsetouch(DIMV(crsetouch),nmat)
      REAL_T, intent(in) :: crsedist(DIMV(crsedist),nmat*(1+SDIM))

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: rz_flag
      INTEGER_T, intent(in) :: bc(SDIM,2)
      REAL_T, intent(in) :: time

      INTEGER_T i,j,k,i1,j1,k1
      REAL_T vcenter(nmat)
      REAL_T VFRAC_TEMP
      INTEGER_T nhalf
      REAL_T xsten_accept(-3:3,SDIM)
      REAL_T xsten_donate(-3:3,SDIM)
      INTEGER_T dir,dir2,side
      INTEGER_T local_facetest

      INTEGER_T isten,jsten,ksten
      INTEGER_T iside,jside,kside
      INTEGER_T klosten,khisten
      INTEGER_T fluid_materials_in_cell_stencil
      INTEGER_T im
      INTEGER_T vofcomp
      REAL_T newfab_hold(nmat*(1+SDIM))
      INTEGER_T touch_hold(nmat)
      INTEGER_T ilocut(3),ihicut(3)
      INTEGER_T icur(3)
      INTEGER_T istar_array(3)
      INTEGER_T istar_array_offset(3)
      INTEGER_T sorted_list(nmat)
      INTEGER_T nten_test,nstar_test,nface_test
      INTEGER_T FSI_exclude
      REAL_T mofdata(nmat*ngeom_recon)
      INTEGER_T istar
      INTEGER_T im_crit
      INTEGER_T stencil_test(nmat)
      INTEGER_T cell_test(nmat)
      INTEGER_T face_test(nmat)
      INTEGER_T stringent_test_passed(nmat)
      INTEGER_T full_neighbor(nmat)
      INTEGER_T ii,jj,kk
      INTEGER_T i3,j3,k3
      INTEGER_T i4,j4,k4
      INTEGER_T i4low(3)
      INTEGER_T i4high(3)
      INTEGER_T i4_array(3)
      INTEGER_T im_corner
      INTEGER_T im_test_stencil
      INTEGER_T im_test_center
      REAL_T FSUM(nmat)
      INTEGER_T on_border
       ! 1..nmat,fluid materials in cell,nstar
       ! donateflag(nmat+2 ... nmat+1+nstar)=fluid material id that owns
       ! the respective star stencil position.
       ! donateflag(1..nmat)=1 if the respective material is a fluid
       ! which has a "non-flotsam" presence in the cell.
       ! If (cell_test(im)==1) and 
       !    ((keep_all_interfaces==1)or
       !     (truncate_volume_fractions(im)==0)) and
       !     (the local star point is owned by material im) then
       !  donateflag(nmat+1+istar)=im 
      INTEGER_T donateflag(nmat+1+nstar)
      INTEGER_T crse_dist_valid
      INTEGER_T ctouch
      REAL_T init_dist_from_crse
      INTEGER_T rigid_in_stencil
      REAL_T VFRAC_STENCIL_CUTOFF
      REAL_T VFRAC_INTERP,theta_nbr,theta_cen
      REAL_T xnbr(SDIM)
      REAL_T xmid(SDIM)
      INTEGER_T i_DEB_DIST
      INTEGER_T j_DEB_DIST
      INTEGER_T k_DEB_DIST

      INTEGER_T versionA_test
      INTEGER_T versionB_test
      INTEGER_T version_of_choice
      INTEGER_T keep_flotsam
      INTEGER_T legitimate_material
      INTEGER_T height_check(nmat)
      INTEGER_T boundary_face_count(nmat)
      INTEGER_T center_face_count(nmat)
      INTEGER_T iface,jface,kface
      INTEGER_T f_index(3)
      INTEGER_T ml,mr,ifacepair
      REAL_T frac_pair(nmat,nmat)  !(m_left,m_right)

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

      if ((keep_all_interfaces.eq.0).or. &
          (keep_all_interfaces.eq.1)) then
       ! do nothing
      else
       print *,"keep_all_interfaces invalid"
       stop
      endif

      if (max_problen.le.zero) then
       print *,"max_problen invalid"
       stop
      endif

       ! areas and centroids on face
      nface_test=nmat*SDIM*2*(1+SDIM)
      if (nface_test.ne.nface) then
       print *,"nface invalid levelstrip nface nface_test ",nface,nface_test
       stop
      endif
      if (nface_dst.eq.2*nmat*nmat) then
       ! do nothing
      else
       print *,"nface_dst invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid levelstrip nten nten_test ",nten,nten_test
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
      if (ngrow_make_distance.ne.3) then
       print *,"ngrow_make_distance<>3 error in levelstrip"
       stop
      endif

      do im=1,nmat
       if (truncate_volume_fractions(im).eq.0) then
        ! do nothing
       else if (truncate_volume_fractions(im).eq.1) then
        ! do nothing
       else
        print *,"truncate_volume_fractions invalid"
        stop
       endif
      enddo !im=1..nmat

      call checkbound(fablo,fabhi,DIMS(maskfab),ngrow_distance,-1,2878)
      call checkbound(fablo,fabhi,DIMS(facepairX),ngrow_distance,0,1871)
      call checkbound(fablo,fabhi,DIMS(facepairY),ngrow_distance,1,1871)
      call checkbound(fablo,fabhi,DIMS(facepairZ),ngrow_distance,SDIM-1,1871)
      call checkbound(fablo,fabhi,DIMS(facefab),ngrow_distance,-1,1871)
      call checkbound(fablo,fabhi,DIMS(facetest),ngrow_distance,-1,1872)
      call checkbound(fablo,fabhi,DIMS(stenfab),ngrow_distance,-1,1873)
      call checkbound(fablo,fabhi,DIMS(vofrecon),ngrow_distance,-1,1874)
      call checkbound(fablo,fabhi,DIMS(origdist),ngrow_distance,-1,1875)
      call checkbound(fablo,fabhi,DIMS(newfab),1,-1,1876)
      call checkbound(fablo,fabhi,DIMS(touchfab),0,-1,1876)
      call checkbound(fablo,fabhi,DIMS(crsetouch),0,-1,1876)
      call checkbound(fablo,fabhi,DIMS(crsedist),0,-1,1876)
      
      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif
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

      if (rz_flag.eq.0) then
       ! do nothing
      else if (rz_flag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.3) then
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
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

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

       do im=1,nmat

        if (is_rigid(nmat,im).eq.0) then

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
           print *,"im,im_crit, or init_dist_from_crse invalid"
           stop
          endif 
          ctouch=NINT(crsetouch(D_DECL(i,j,k),im))
          if (ctouch.eq.0) then
           crse_dist_valid=0
          else if ((ctouch.eq.1).or.(ctouch.eq.2)) then
           ! do nothing
          else
           print *,"ctouch invalid"
           stop
          endif
         else
          print *,"level invalid"
          stop
         endif
        
         if (crse_dist_valid.eq.1) then
          do dir=1,SDIM
           dir2=nmat+(im-1)*SDIM+dir
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
          print *,"crse_dist_valid invalid"
          stop
         endif
 
         newfab(D_DECL(i,j,k),im)=init_dist_from_crse

         if ((i.eq.i_DEB_DIST).and. &
             (j.eq.j_DEB_DIST).and. &
             (k.eq.k_DEB_DIST)) then
          print *,"DEB_DIST: im,init_dist_from_crse=", &
                  im,init_dist_from_crse
         endif 

        else if (is_rigid(nmat,im).eq.1) then

         init_dist_from_crse=newfab(D_DECL(i,j,k),im)
         touchfab(D_DECL(i,j,k),im)=two
         if (init_dist_from_crse.lt.minLS(im)) then
          minLS(im)=init_dist_from_crse
         endif
         if (init_dist_from_crse.gt.maxLS(im)) then
          maxLS(im)=init_dist_from_crse
         endif

        else
         print *,"is_rigid invalid"
         stop
        endif

       enddo ! im=1..nmat

      enddo
      enddo
      enddo  ! initialize + or -

      call growntilebox_TILE(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow_make_distance)
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       on_border=1-NINT(maskfab(D_DECL(i,j,k),3))

       nprocessed=nprocessed+1

       call gridsten_level(xsten_donate,i,j,k,level,nhalf)

       ! im=1..nmat: donateflag(im)=1 => find closest distance to the 
       !             im interface.
       ! donateflag(nmat+1)=number fluid materials in the cell.
       ! donateflag(nmat+2 ... nmat+1+nstar)=fluid material id that owns
       ! the respective star stencil position.
       ! If (cell_test(im)==1) and 
       !    ((keep_all_interfaces==1)or
       !     (truncate_volume_fractions(im)==0)) and
       !     (the LOCAL star point is owned by material im) then
       !  donateflag(nmat+1+istar)=im 

       do im=1,nmat+1+nstar
        donateflag(im)=0
       enddo

       fluid_materials_in_cell_stencil=0

       do im=1,nmat*ngeom_recon
        mofdata(im)=vofrecon(D_DECL(i,j,k),im)
       enddo

       do im=1,nmat
        vofcomp=(im-1)*ngeom_recon+1
        vcenter(im)=mofdata(vofcomp)
       enddo

       FSI_exclude=1
       call sort_volume_fraction(vcenter,FSI_exclude,sorted_list,nmat)
       im_crit=sorted_list(1)
       if ((im_crit.ge.1).and.(im_crit.le.nmat)) then
        ! do nothing
       else
        print *,"im_crit invalid"
        stop
       endif
       if (is_rigid(nmat,im_crit).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif

        ! in Sussman and Puckett, 
        ! a reconstructed interface segment was valid if
        ! LS_{ij}(LS_{ij}+LS_{i+i',j+j'}) <= 0 for
        ! some |i'|<=1 |j'|<=1
        ! in the present algorithm, a reconstructed segment is valid if
        ! (A) F>VOFTOL_REDIST and
        ! (B)        
       do im=1,nmat
        cell_test(im)=0 !F>VOFTOL_REDIST?
        face_test(im)=0 !face areafrac between cells consistent?
        full_neighbor(im)=0 !neighbor F(im)>1-facetol ?
        stencil_test(im)=0 ! F(im)>1/2-eps on cell bdry?
        stringent_test_passed(im)=0 ! stenfab consistent between cells?
       enddo

         ! initialize: cell_test
       do im=1,nmat
        if (is_rigid(nmat,im).eq.0) then
         if (vcenter(im).gt.VOFTOL_REDIST) then
          cell_test(im)=1
         endif
        else if (is_rigid(nmat,im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid"
         stop
        endif

        if ((i.eq.i_DEB_DIST).and. &
            (j.eq.j_DEB_DIST).and. &
            (k.eq.k_DEB_DIST)) then
         print *,"DEB_DIST: im,cell_test=",im,cell_test(im)
        endif 

       enddo  ! im=1..nmat

       do dir=1,3
        istar_array(dir)=0
       enddo
       call put_istar(istar,istar_array)
       im_test_center=NINT(stenfab(D_DECL(i,j,k),istar))
       if ((im_test_center.ge.1).and.(im_test_center.le.nmat)) then
        ! do nothing
       else
        print *,"im_test_center out of range (0)"
        stop
       endif
        ! im_test_center must be a fluid material.
       if (is_rigid(nmat,im_test_center).eq.0) then 
        ! do nothing
       else
        print *,"is_rigid(nmat,im_test_center).ne.0 (0)"
        stop
       endif
        ! 1..nmat,fluid materials in cell, nstar
       donateflag(nmat+1+istar)=im_test_center

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
        do i1=1,3
        do j1=1,3
        do k1=klosten+2,khisten+2
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

         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          VFRAC_TEMP=vofrecon(D_DECL(isten,jsten,ksten),vofcomp)
          if (is_rigid(nmat,im).eq.0) then
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
          else if (is_rigid(nmat,im).eq.1) then
           if (VFRAC_TEMP.ge.VOFTOL_REDIST) then
            rigid_in_stencil=1
           endif
          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif
         enddo  ! im=1..nmat

        enddo
        enddo
        enddo  ! i1,j1,k1 (test VFRAC_INTERP>=VFRAC_STENCIL_CUTOFF)

        if ((i.eq.i_DEB_DIST).and. &
            (j.eq.j_DEB_DIST).and. &
            (k.eq.k_DEB_DIST)) then
         do im=1,nmat
          print *,"DEB_DIST: im,stencil_test=",im,stencil_test(im)
         enddo
        endif

         ! face_test
        do im=1,nmat
         if (is_rigid(nmat,im).eq.0) then

          if (cell_test(im).eq.1) then

           if ((stencil_test(im).eq.1).or. &
               (rigid_in_stencil.eq.0)) then
 
            face_test(im)=1
            do dir=1,SDIM
             ii=0
             jj=0
             kk=0
             if (dir.eq.1) then
              ii=1
             else if (dir.eq.2) then
              jj=1
             else if ((dir.eq.3).and.(SDIM.eq.3)) then
              kk=1
             else
              print *,"dir invalid levelstrip"
              stop
             endif

             do side=1,2

              if (side.eq.1) then
               iside=i
               jside=j
               kside=k
              else if (side.eq.2) then
               iside=i+ii
               jside=j+jj
               kside=k+kk
              else
               print *,"side invalid"
               stop
              endif

              local_facetest= &
               NINT(facetest(D_DECL(iside,jside,kside),(dir-1)*nmat+im))
              if (local_facetest.eq.0) then
               face_test(im)=0
              else if (local_facetest.eq.1) then
               ! do nothing
              else
               print *,"local_facetest invalid"
               stop
              endif

             enddo ! side=1,2

            enddo ! dir=1..sdim

           else if ((stencil_test(im).eq.0).and. &
                    (rigid_in_stencil.eq.1)) then
            ! do nothing
           else
            print *,"stencil_test or rigid_in_stencil invalid"
            stop
           endif

          else if (cell_test(im).eq.0) then
           ! do nothing
          else
           print *,"cell_test invalid"
           stop
          endif
         else if (is_rigid(nmat,im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif
        enddo ! im=1..nmat


         ! stringent_test_passed (on_border==0):
         ! investigate all points coinciding at the intersection of the
         ! line connecting cells (i,j,k) and (i+i3,j+j3,k+k3) with the
         ! boundary of cell (i,j,k). 
         ! 
        do i3=-1,1
        do j3=-1,1
        do k3=klosten,khisten
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
         do im=1,nmat
          FSUM(im)=zero
         enddo
         ! (i4,j4,k4)=(0,0,0) if (i3,j3,k3)=(0,0,0)
         do i4=i4low(1),i4high(1)
         do j4=i4low(2),i4high(2)
         do k4=i4low(3),i4high(3)
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

           ! STENINIT called with tessellate==0 prior to LEVELSTRIP.
           ! im_test_stencil is the material that occupies the
           ! node in question.
          call put_istar(istar,istar_array_offset)
          im_test_stencil=NINT(stenfab(D_DECL(isten,jsten,ksten),istar))
          if ((im_test_stencil.ge.1).and.(im_test_stencil.le.nmat)) then
           ! do nothing
          else
           print *,"im_test_stencil out of range 1"
           stop
          endif
          if (is_rigid(nmat,im_test_stencil).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(nmat,im_test_stencil).ne.0 (1)"
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
          else if ((im_corner.ge.1).and.(im_corner.le.nmat)) then
           if (im_corner.eq.im_test_stencil) then
            ! do nothing
           else
            im_corner=-1
           endif
          else if (im_corner.eq.-1) then
           ! do nothing
          else
           print *,"im_corner invalid"
           stop
          endif
          do im=1,nmat
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
          call sort_volume_fraction(FSUM,FSI_exclude,sorted_list,nmat)
          im_corner=sorted_list(1)
          if (is_rigid(nmat,im_corner).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(nmat,im_corner).ne.0"
           stop
          endif
         else if (im_corner.eq.0) then
          print *,"im_corner invalid"
          stop
         else if ((im_corner.ge.1).and.(im_corner.le.nmat)) then
          ! do nothing
         else
          print *,"im_corner invalid"
          stop
         endif
         if ((im_corner.ge.1).and. &
             (im_corner.le.nmat)) then

           ! a tougher test for flotsam near contact lines.
          if ((stencil_test(im_corner).eq.1).or. &
              (rigid_in_stencil.eq.0)) then
           call put_istar(istar,istar_array) 
           donateflag(nmat+1+istar)=im_corner
           stringent_test_passed(im_corner)=1
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
          print *,"DEB_DIST: i3,j3,k3,istar,donateflag(nmat+1+istar) ", &
           i3,j3,k3,istar,donateflag(nmat+1+istar)
         endif

        enddo
        enddo
        enddo  ! i3,j3,k3

         ! dir=1..sdim, side=1..2
         ! on_border==0
        do im=1,nmat
         height_check(im)=0
        enddo

        do dir=1,SDIM

         ii=0
         jj=0
         kk=0
         if (dir.eq.1) then
          ii=1
         else if (dir.eq.2) then
          jj=1
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          kk=1
         else
          print *,"dir invalid levelstrip"
          stop
         endif

         do side=1,2

          do im=1,nmat
           boundary_face_count(im)=0
           center_face_count(im)=0
          enddo

          if (side.eq.1) then
           iface=i
           jface=j
           kface=k
          else if (side.eq.2) then
           iface=i+ii
           jface=j+jj
           kface=k+kk
          else
           print *,"side invalid"
           stop
          endif
          do i3=-1,1
          do j3=-1,1
          do k3=klosten,khisten
      
           f_index(1)=i3 
           f_index(2)=j3 
           f_index(3)=k3 

           i4=i3+iface
           j4=j3+jface
           k4=k3+kface

           ifacepair=1
           do ml = 1, nmat
           do mr = 1, nmat
            if (dir.eq.1) then
             frac_pair(ml,mr)=facepairX(D_DECL(i4,j4,k4),ifacepair)
            else if (dir.eq.2) then
             frac_pair(ml,mr)=facepairY(D_DECL(i4,j4,k4),ifacepair)
            else if ((dir.eq.3).and.(SDIM.eq.3)) then
             frac_pair(ml,mr)=facepairZ(D_DECL(i4,j4,k4),ifacepair)
            else
             print *,"dir invalid"
             stop
            endif
            ifacepair=ifacepair+2
           enddo
           enddo
           if (ifacepair.eq.nface_dst+1) then
            ! do nothing
           else
            print *,"ifacepair invalid"
            stop
           endif 

           do im=1,nmat
            if (is_rigid(nmat,im).eq.0) then
             if ((frac_pair(im,im).ge.VOFTOL_REDIST).and. &
                 (frac_pair(im,im).le.one+VOFTOL_REDIST)) then 
              if ((i3.eq.0).and. &
                  (j3.eq.0).and. &
                  (k3.eq.0).and. &
                  (f_index(dir).eq.0)) then
               center_face_count(im)=1
              else if (f_index(dir).eq.0) then
               ! do nothing 
              else if ((f_index(dir).eq.1).or. &
                       (f_index(dir).eq.-1)) then
               boundary_face_count(im)=1
              else
               print *,"f_index invalid"
               stop
              endif
             else if ((frac_pair(im,im).ge.zero).and. &
                      (frac_pair(im,im).le.VOFTOL_REDIST)) then
              ! do nothing
             else
              print *,"frac_pair invalid"
              stop
             endif
            else if (is_rigid(nmat,im).eq.1) then
             ! do nothing
            else
             print *,"is_rigid(nmat,im) invalid"
             stop
            endif
           enddo ! im=1..nmat
          enddo
          enddo
          enddo ! i3,j3,k3
          do im=1,nmat
           if ((center_face_count(im).eq.1).and. &
               (boundary_face_count(im).eq.1)) then
            height_check(im)=1
           else if ((center_face_count(im).eq.0).or. &
                    (boundary_face_count(im).eq.0)) then
            ! do nothing
           else
            print *,"center_face_count or boundary_face_count bad"
            stop
           endif
          enddo ! im=1..nmat
         enddo ! side=1,2     
        enddo ! dir=1..sdim
        
         ! full_neighbor
        do i3=-1,1
        do j3=-1,1
        do k3=klosten,khisten
         if ((i3.eq.0).and.(j3.eq.0).and.(k3.eq.0)) then
          ! do nothing
         else if ((abs(i3).eq.1).or.(abs(j3).eq.1).or.(abs(k3).eq.1)) then
          iside=i+i3
          jside=j+j3
          kside=k+k3
          istar_array(1)=i3
          istar_array(2)=j3
          istar_array(3)=k3

          do im=1,nmat
           if (is_rigid(nmat,im).eq.0) then
            vofcomp=(im-1)*ngeom_recon+1
            VFRAC_TEMP=vofrecon(D_DECL(iside,jside,kside),vofcomp)
            if (VFRAC_TEMP.ge.one-VOFTOL) then
             full_neighbor(im)=1
             call put_istar(istar,istar_array) 
             donateflag(nmat+1+istar)=im
            endif
           else if (is_rigid(nmat,im).eq.1) then
            ! do nothing
           else
            print *,"is_rigid invalid"
            stop
           endif
          enddo ! im=1..nmat 

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
        print *,"on_border invalid"
        stop
       endif

       do im=1,nmat

        if ((i.eq.i_DEB_DIST).and. &
            (j.eq.j_DEB_DIST).and. &
            (k.eq.k_DEB_DIST)) then
         print *,"DEB_DIST: im,stringent_test_passed,face_test ", &
          im,stringent_test_passed(im),face_test(im)
        endif

         ! face_test=0 if cell_test==0
        if (is_rigid(nmat,im).eq.0) then

         if ((stringent_test_passed(im).eq.1).or. &
             (face_test(im).eq.1)) then
          versionA_test=1
         else
          versionA_test=0
         endif

         if ((height_check(im).eq.1).and. &
             (cell_test(im).eq.1)) then
          versionB_test=1
         else
          versionB_test=0
         endif

         if (1.eq.1) then
          version_of_choice=versionB_test
         else
          version_of_choice=versionA_test
         endif

         keep_flotsam=0
         if (cell_test(im).eq.1) then
          if ((keep_all_interfaces.eq.1).or. &
              (truncate_volume_fractions(im).eq.0)) then
           keep_flotsam=1
          else if ((keep_all_interfaces.eq.0).and. &
                   (truncate_volume_fractions(im).eq.1)) then
           keep_flotsam=0
          else
           print *,"keep_all_interfaces, truncate_volume_fraction, err"
           stop
          endif
         else if (cell_test(im).eq.0) then
          keep_flotsam=0
         else
          print *,"cell_test invalid"
          stop
         endif

         legitimate_material=0
         if ((vcenter(im).ge.half).or. &
             (im.eq.im_crit).or. &
             (version_of_choice.eq.1).or. &
             (full_neighbor(im).eq.1).or. &
             (keep_flotsam.eq.1)) then
          legitimate_material=1
         else if ((vcenter(im).le.half).and. &
                  (im.ne.im_crit).and. &
                  (version_of_choice.eq.0).and. &
                  (full_neighbor(im).eq.0).and. &
                  (keep_flotsam.eq.0)) then
          legitimate_material=0
         else
          print *,"legitimate check failed"
          stop
         endif

         if (legitimate_material.eq.1) then

          do i3=-1,1
          do j3=-1,1
          do k3=klosten,khisten
           istar_array(1)=i3
           istar_array(2)=j3
           istar_array(3)=k3
           call put_istar(istar,istar_array)
           if ((istar.ge.1).and.(istar.le.nstar)) then
            im_test_stencil=NINT(stenfab(D_DECL(i,j,k),istar))
            if (is_rigid(nmat,im_test_stencil).eq.0) then
             if (im_test_stencil.eq.im) then
               ! 1<=istar<=nstar
              donateflag(nmat+1+istar)=im
             endif
            else
             print *,"is_rigid(nmat,im_test_stencil).ne.0 (1)"
             stop
            endif
           else
            print *,"istar invalid"
            stop
           endif
          enddo ! k3
          enddo ! j3
          enddo ! i3

         else if (legitimate_material.eq.0) then
          ! do nothing
         else
          print *,"legitimate_material invalid"
          stop
         endif

         if (legitimate_material.eq.1) then
          fluid_materials_in_cell_stencil= &
           fluid_materials_in_cell_stencil+1
          donateflag(im)=1
         else if (legitimate_material.eq.0) then
          ! do nothing
         else
          print *,"legitimate_material invalid"
          stop
         endif 

        else if (is_rigid(nmat,im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid"
         stop
        endif
       enddo ! im=1..nmat

       donateflag(nmat+1)=fluid_materials_in_cell_stencil

       if ((i.eq.i_DEB_DIST).and. &
           (j.eq.j_DEB_DIST).and. &
           (k.eq.k_DEB_DIST)) then
        print *,"DEB_DIST: fluid_materials_in_cell_stencil=", &
         fluid_materials_in_cell_stencil
       endif

       if ((fluid_materials_in_cell_stencil.lt.1).or. &
           (fluid_materials_in_cell_stencil.gt.nmat)) then
        print *,"fluid_materials_in_cell_stencil invalid:", &
          fluid_materials_in_cell_stencil
        stop
       else if (fluid_materials_in_cell_stencil.eq.1) then
         ! do nothing
       else if ((fluid_materials_in_cell_stencil.ge.2).and. &
                (fluid_materials_in_cell_stencil.le.nmat)) then

        ilocut(3)=0
        ihicut(3)=0
        icur(1)=i
        icur(2)=j
        icur(3)=k

         ! we do not have to check outside the tile since tiles
         ! outside the present tile might traverse interface segments which
         ! are in the present tile.
        do dir=1,SDIM
         ilocut(dir)=-ngrow_make_distance
         ihicut(dir)=ngrow_make_distance
         if (icur(dir)+ilocut(dir).lt.tilelo(dir)) then
          ilocut(dir)=tilelo(dir)-icur(dir)
         endif
         if (icur(dir)+ihicut(dir).gt.tilehi(dir)) then
          ihicut(dir)=tilehi(dir)-icur(dir)
         endif
        enddo ! dir

        do i1=ilocut(1),ihicut(1)
        do j1=ilocut(2),ihicut(2)
        do k1=ilocut(3),ihicut(3)

         call gridsten_level(xsten_accept,i+i1,j+j1,k+k1,level,nhalf)

         do im=1,nmat*(1+SDIM)
          newfab_hold(im)=newfab(D_DECL(i+i1,j+j1,k+k1),im)
         enddo  ! im
         do im=1,nmat
          touch_hold(im)=NINT(touchfab(D_DECL(i+i1,j+j1,k+k1),im))
         enddo

         call update_closest( &
          xsten_accept,xsten_donate,nhalf, &
          dx,xlo,bfact,level,fablo, &
          mofdata,nmat,nstar, &
          i,j,k, &  ! donate index
          i1,j1,k1, & ! accept index: i+i1,j+j1,k+k1
          newfab_hold, &
          touch_hold, &
          minLS, &
          maxLS, &
          donateflag, &
          time)

         do im=1,nmat*(1+SDIM)
          newfab(D_DECL(i+i1,j+j1,k+k1),im)=newfab_hold(im)
         enddo  ! im
         do im=1,nmat
          touchfab(D_DECL(i+i1,j+j1,k+k1),im)=touch_hold(im)
         enddo

        enddo
        enddo
        enddo ! i1,j1,k1

       else
        print *,"parameter bust in FORT_LEVELSTRIP"
        print *,"fluid_materials_in_cell_stencil=", &
          fluid_materials_in_cell_stencil
        stop
       endif

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine FORT_LEVELSTRIP


      subroutine FORT_CORRECT_UNINIT( &
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
         time, &
         nmat)

      use global_utility_module
      use probcommon_module
      use mof_redist_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: minLS(nmat)
      REAL_T, intent(in) :: maxLS(nmat)
      REAL_T, intent(in) :: max_problen
      INTEGER_T, intent(in) :: DIMDEC(newfab)
      INTEGER_T, intent(in) :: DIMDEC(touchfab)

      REAL_T, intent(inout) :: newfab(DIMV(newfab),nmat*(1+SDIM))
      REAL_T, intent(in) :: touchfab(DIMV(touchfab),nmat)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      REAL_T, intent(in) :: time

      INTEGER_T i,j,k
      INTEGER_T im
      INTEGER_T ctouch
      REAL_T init_dist
     
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

      call checkbound(fablo,fabhi,DIMS(newfab),1,-1,2876)
      call checkbound(fablo,fabhi,DIMS(touchfab),0,-1,2876)
      
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0)
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       do im=1,nmat

        ctouch=NINT(touchfab(D_DECL(i,j,k),im))

        if (is_rigid(nmat,im).eq.0) then

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
           print *,"init_dist bust"
           stop
          endif
          newfab(D_DECL(i,j,k),im)=init_dist
         else if ((ctouch.eq.1).or.(ctouch.eq.2)) then
          ! do nothing
         else
          print *,"ctouch invalid"
          stop
         endif
 
        else if (is_rigid(nmat,im).eq.1) then

         if (ctouch.ne.2) then
          print *,"ctouch invalid"
          stop
         endif

        else
         print *,"is_rigid invalid"
         stop
        endif

       enddo ! im=1..nmat

      enddo
      enddo
      enddo  ! replace uninit with minLS or maxLS

      return
      end subroutine FORT_CORRECT_UNINIT


      subroutine FORT_STENINIT( &
       level, &
       finest_level, &
       stenfab,DIMS(stenfab), &
       maskfab,DIMS(maskfab), &
       vofrecon,DIMS(vofrecon), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       rz_flag, &
       xlo,dx, &
       time,ngrow_distance, &
       nmat,nstar)

      use global_utility_module
      use probcommon_module
      use MOF_routines_module
      use mof_redist_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nstar
      INTEGER_T, intent(in) :: nmat,ngrow_distance
      INTEGER_T, intent(in) :: DIMDEC(stenfab)
      INTEGER_T, intent(in) :: DIMDEC(maskfab)
      INTEGER_T, intent(in) :: DIMDEC(vofrecon)

      REAL_T, intent(out) :: stenfab(DIMV(stenfab),nstar)
      REAL_T, intent(in) :: maskfab(DIMV(maskfab),2)
      REAL_T, intent(in) :: vofrecon(DIMV(vofrecon),nmat*ngeom_recon)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: rz_flag
      REAL_T, intent(in) :: time

      INTEGER_T i,j,k
      REAL_T vcenter(nmat)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      INTEGER_T klosten,khisten
      INTEGER_T im
      INTEGER_T vofcomp
      INTEGER_T istar_array(3)
      INTEGER_T nstar_test
      REAL_T mofdata(nmat*ngeom_recon)
      INTEGER_T istar
      INTEGER_T im_crit
      INTEGER_T i3,j3,k3
      REAL_T xcorner(SDIM)
      INTEGER_T im_test
      INTEGER_T dir
      INTEGER_T mask1,mask2
      INTEGER_T tessellate
 
      nhalf=3      
    
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

      call checkbound(fablo,fabhi,DIMS(stenfab),ngrow_distance,-1,2877)
      call checkbound(fablo,fabhi,DIMS(maskfab),ngrow_distance,-1,2878)
      call checkbound(fablo,fabhi,DIMS(vofrecon),ngrow_distance,-1,2879)
      
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
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

      if (rz_flag.eq.0) then
       ! do nothing
      else if (rz_flag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.3) then
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
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       ! mask1=1 at interior cells or fine/fine ghost cells
       ! mask1=0 at coarse/fine ghost cells or outside domain.
       ! mask2=1 at interior cells
       mask1=NINT(maskfab(D_DECL(i,j,k),1))
       mask2=NINT(maskfab(D_DECL(i,j,k),2))

       if ((mask2.eq.1).or.(mask1.eq.0)) then

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do im=1,nmat*ngeom_recon
         mofdata(im)=vofrecon(D_DECL(i,j,k),im)
        enddo

         ! vcenter = volume fraction 
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         vcenter(im)=mofdata(vofcomp)
        enddo ! im

         ! uses VOFTOL
        call check_full_cell_vfrac( &
                vcenter, &
                tessellate, &  ! =0
                nmat,im_crit)

        if ((im_crit.ge.1).and.(im_crit.le.nmat)) then

         do istar=1,nstar
          stenfab(D_DECL(i,j,k),istar)=im_crit
         enddo

        else if (im_crit.eq.0) then

         do i3=-1,1
         do j3=-1,1
         do k3=klosten,khisten

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
           im_test,nmat,SDIM)
          if ((im_test.ge.1).and.(im_test.le.nmat)) then
           call put_istar(istar,istar_array) 
           stenfab(D_DECL(i,j,k),istar)=im_test
          else
           print *,"im_test invalid"
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
      end subroutine FORT_STENINIT


       ! for finding areas internal to a cell, perturb each internal 
       ! interface, find areas and volumes, then check for the difference 
       ! in volumes divided by eps times the area.
       ! 
      subroutine FORT_FACEINIT( &
         tid, &
         tessellate, &  ! =0,1, or 3
         nten, &
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
         nmat, &
         nface, &
         nface_decomp)

      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module
      use mof_redist_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: tessellate ! 0,1, or 3
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nface
      INTEGER_T, intent(in) :: nface_decomp
      INTEGER_T :: nface_decomp_test
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: DIMDEC(facefab)
      INTEGER_T, intent(in) :: DIMDEC(maskfab)
      INTEGER_T, intent(in) :: DIMDEC(vofrecon)

      REAL_T, intent(out) :: facefab(DIMV(facefab),nface+nface_decomp)
      REAL_T, intent(in) :: maskfab(DIMV(maskfab),2)
      REAL_T, intent(in) :: vofrecon(DIMV(vofrecon),nmat*ngeom_recon)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: rz_flag
      REAL_T, intent(in) :: time

      INTEGER_T i,j,k
      INTEGER_T im,im_opp
      INTEGER_T im1,im2
      INTEGER_T im_crit
      INTEGER_T im_crit_thin
      INTEGER_T iface
      INTEGER_T dir,side

      INTEGER_T nface_test
      INTEGER_T vofcomp
      REAL_T vcenter(nmat)
      REAL_T vcenter_thin(nmat)
      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T mofdatavalid(nmat*ngeom_recon)
      REAL_T mofdataproject(nmat*ngeom_recon)
      REAL_T localface(nmat,SDIM,2,SDIM+1)
      REAL_T localface_line(nmat,nmat,SDIM+1,SDIM,2)
      REAL_T totalface_line(nmat,SDIM,2)
      REAL_T totalface(SDIM,2,SDIM+1)

      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten_thin(-1:1,SDIM)
      REAL_T dxthin
      INTEGER_T nhalf,nhalf_thin,isten
      REAL_T dummy_tri(SDIM+1,SDIM)
      INTEGER_T nmax,ivert
      INTEGER_T dir2
      INTEGER_T shapeflag
      REAL_T multi_volume(nmat)
      REAL_T multi_cen(SDIM,nmat)
      REAL_T multi_area(nmat)
      REAL_T multi_area_line(nmat)
      REAL_T multi_volume_offset(nmat)
      REAL_T multi_cen_offset(SDIM,nmat)
      REAL_T total_vol
      REAL_T total_line
      INTEGER_T mask1,mask2
      REAL_T areacen(SDIM)
      REAL_T uncaptured_volume_fraction
      REAL_T vfrac_fluid_sum
      REAL_T vfrac_solid_sum
      INTEGER_T loop_counter
      INTEGER_T num_processed_fluid
      INTEGER_T num_processed_solid
      INTEGER_T nmat_fluid
      INTEGER_T nmat_rigid
      INTEGER_T testflag
      INTEGER_T iten
      INTEGER_T nten_test
      INTEGER_T is_processed(nten)
      REAL_T slope(SDIM)
      REAL_T intercept
      REAL_T local_linefrac(nmat)
      REAL_T F1,F2,FSTRIP
      REAL_T X1(SDIM)
      REAL_T X2(SDIM)
      REAL_T XSTRIP(SDIM)
      INTEGER_T normalize_tessellate
      INTEGER_T local_tessellate
      INTEGER_T is_rigid_local(nmat)
      INTEGER_T nhalf_box
      INTEGER_T cmofsten(D_DECL(-1:1,-1:1,-1:1))

      nhalf_box=1
  
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

      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid FACEINIT nten nten_test ",nten,nten_test
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

       ! area+centroid for all faces of a cell.
       ! (nmat,sdim,2,sdim+1)
      nface_test=nmat*SDIM*2*(1+SDIM)
      if (nface_test.ne.nface) then
       print *,"nface invalid faceinit nface nface_test ",nface,nface_test
       stop
      endif

       ! areas and centroids after projecting surface to cell faces.
       ! inside,outside,area+centroid,dir,side 
       ! (nmat,nmat,sdim+1,sdim,2)
      nface_decomp_test=nmat*nmat*(SDIM+1)*SDIM*2

      if ((nface_decomp.ne.0).and. &
          (nface_decomp.ne.nface_decomp_test)) then
       print *,"nface_decomp invalid"
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

      call checkbound(fablo,fabhi,DIMS(facefab),ngrow,-1,2883)
      call checkbound(fablo,fabhi,DIMS(maskfab),ngrow,-1,2884)
      call checkbound(fablo,fabhi,DIMS(vofrecon),ngrow,-1,2885)
      
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
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

      if (rz_flag.eq.0) then
       ! do nothing
      else if (rz_flag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.3) then 
       ! do nothing
      else
       print *,"rz_flag invalid in faceinit"
       stop
      endif
 
      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

        ! mask1=1 at interior cells or fine/fine ghost cells
        ! mask2=1 at interior cells
       mask1=NINT(maskfab(D_DECL(i,j,k),1))
       mask2=NINT(maskfab(D_DECL(i,j,k),2))

       if ((mask2.eq.1).or.(mask1.eq.0)) then

         ! xsten(0,dir)  dir=1,2  is center of cell
         ! xsten(1,dir)  is right coordinate in dir direction
         ! xsten(-1,dir) is left coordinate in dir direction.
        call gridsten_level(xsten,i,j,k,level,nhalf)

        do im=1,nmat*ngeom_recon
         mofdata(im)=vofrecon(D_DECL(i,j,k),im)
        enddo

        normalize_tessellate=0
        call make_vfrac_sum_ok_copy( &
         cmofsten, &
         xsten,nhalf,nhalf_box, &
         bfact,dx, &
         normalize_tessellate, &  ! =0
         mofdata,mofdatavalid, &
         nmat,SDIM,3000)

        if (tessellate.eq.3) then

         local_tessellate=2

         call multi_get_volume_tessellate( &
          tessellate, &  ! =3
          bfact,dx, &
          xsten,nhalf, &
          mofdatavalid, &
          geom_xtetlist_uncapt(1,1,1,tid+1), &
          nmax, &
          nmax, &
          nmat, &
          SDIM, &
          3)  ! caller_id=3

        else if ((tessellate.eq.0).or. &
                 (tessellate.eq.1)) then
         local_tessellate=tessellate
        else
         print *,"tessellate invalid"
         stop
        endif

          ! vcenter = volume fraction 
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         vcenter(im)=mofdatavalid(vofcomp)
        enddo ! im

        call check_full_cell_vfrac( &
                vcenter, &
                tessellate, & ! 0,1, or 3
                nmat,im_crit)

        do dir=1,SDIM
         do side=1,2
          do dir2=1,SDIM+1
           totalface(dir,side,dir2)=zero
          enddo
         enddo
        enddo

        iface=0
        do im=1,nmat
         do dir=1,SDIM
          do side=1,2
           do dir2=1,SDIM+1
            iface=iface+1
            localface(im,dir,side,dir2)=zero
           enddo ! dir2
          enddo ! side
         enddo ! dir
        enddo ! im

        if (iface.ne.nface) then
         print *,"iface invalid"
         stop
        endif

        if (nface_decomp.eq.0) then
         ! do nothing
        else if (nface_decomp.eq.nface_decomp_test) then

         do im1=1,nmat
          do dir=1,SDIM
           do side=1,2
            totalface_line(im1,dir,side)=zero
           enddo
          enddo
         enddo

         iface=0
         do im1=1,nmat
          do im2=1,nmat
           do dir2=1,SDIM+1
            do dir=1,SDIM
             do side=1,2
              iface=iface+1
              localface_line(im1,im2,dir2,dir,side)=zero
             enddo ! side
            enddo ! dir
           enddo ! dir2
          enddo ! im2
         enddo ! im1
         if (nface_decomp.ne.iface) then
          print *,"nface_decomp invalid"
          stop
         endif
        else
         print *,"nface_decomp invalid"
         stop
        endif 

        do dir=1,SDIM
        do side=1,2
          ! areacen in absolute coordinate system.
         call gridarea(xsten,nhalf,levelrz,dir-1,side-1, &
          totalface(dir,side,1),areacen)
         do dir2=1,SDIM
          totalface(dir,side,1+dir2)=areacen(dir2)
         enddo
        enddo
        enddo ! dir,side

        if ((im_crit.ge.1).and.(im_crit.le.nmat)) then

         do dir=1,SDIM
         do side=1,2
          localface(im_crit,dir,side,1)=one ! area fraction
           ! centroid of face (absolute coordinate system)
          do dir2=1,SDIM
           localface(im_crit,dir,side,1+dir2)=totalface(dir,side,1+dir2)
          enddo
         enddo
         enddo  ! dir,side

         if (nface_decomp.gt.0) then
          ! do nothing, there are no interfaces along faces.
         else if (nface_decomp.eq.0) then
          ! do nothing
         else
          print *,"nface_decomp invalid"
          stop
         endif

        else if (im_crit.eq.0) then

         shapeflag=0

         do dir=1,SDIM
         do side=1,2

          do isten=-1,1
           do dir2=1,SDIM
            xsten_thin(isten,dir2)=xsten(isten,dir2)
           enddo
          enddo ! isten

          dxthin=FACETOL_DVOL*(xsten(1,dir)-xsten(-1,dir))
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
           ! multi_cen(sdim,nmat) is "absolute" 
          call project_slopes_to_face( &
           bfact,dx,xsten,nhalf, &
           mofdatavalid,mofdataproject, &
           nmat,SDIM,dir,side)

           ! base case (also area fractions)
           ! multi_cen in absolute coordinate system (not relative to cell
           ! centroid)
           ! in: FORT_FACEINIT
          call multi_get_volume_grid( &
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
            nmat,SDIM, &
            shapeflag,3) 

          total_vol=zero
          do im=1,nmat
           if ((tessellate.eq.1).or. &
               (tessellate.eq.3)) then
            total_vol=total_vol+multi_volume(im)
           else if (tessellate.eq.0) then
            if (is_rigid(nmat,im).eq.0) then
             total_vol=total_vol+multi_volume(im)
            else if (is_rigid(nmat,im).eq.1) then
             ! do nothing
            else
             print *,"is_rigid invalid"
             stop
            endif
           else
            print *,"tessellate invalid46"
            stop
           endif
          enddo ! im=1..nmat

          if (total_vol.gt.zero) then
           do im=1,nmat
            vcenter_thin(im)=multi_volume(im)/total_vol
            localface(im,dir,side,1)=vcenter_thin(im)
             ! absolute coordinate system.
            do dir2=1,SDIM
             localface(im,dir,side,1+dir2)=multi_cen(dir2,im)
            enddo
             ! centroid in normal direction to face=centroid of face.
             ! absolute coordinate system.
            localface(im,dir,side,1+dir)=totalface(dir,side,1+dir)
           enddo ! im=1..nmat
          else
           print *,"total_vol invalid"
           stop
          endif 

          if (nface_decomp.gt.0) then

           if (tessellate.eq.0) then
            ! do nothing
           else if ((tessellate.eq.1).or. &
                    (tessellate.eq.3)) then
            print *,"expecting nface_decomp=0 if tessellate=1,3"
            stop
           else
            print *,"tessellate invalid"
            stop
           endif

           call check_full_cell_vfrac( &
             vcenter_thin, &
             tessellate, & ! 0,1, or 3
             nmat,im_crit_thin)

           if ((im_crit_thin.ge.1).and.(im_crit_thin.le.nmat)) then
            ! do nothing, there are no internal interfaces
           else if (im_crit_thin.eq.0) then

            do im=1,nmat
             is_rigid_local(im)=is_rigid(nmat,im)
             if (local_tessellate.eq.2) then
              is_rigid_local(im)=0
             else if (local_tessellate.eq.0) then 
              ! do nothing
             else if (local_tessellate.eq.1) then 
              ! do nothing
             else
              print *,"local_tessellate invalid4"
              stop
             endif
            enddo ! im=1..nmat

            do iten=1,nten
             is_processed(iten)=0
            enddo

            uncaptured_volume_fraction=one

            nmat_rigid=0
            nmat_fluid=0
            vfrac_fluid_sum=zero
            vfrac_solid_sum=zero

            do im=1,nmat
             if (is_rigid_local(im).eq.0) then
              nmat_fluid=nmat_fluid+1
              vfrac_fluid_sum=vfrac_fluid_sum+vcenter_thin(im)
             else if (is_rigid_local(im).eq.1) then
              nmat_rigid=nmat_rigid+1
              vfrac_solid_sum=vfrac_solid_sum+vcenter_thin(im)
             else
              print *,"is_rigid_local invalid"
              stop
             endif
            enddo ! im=1..nmat

            if ((tessellate.eq.0).or. &
                (tessellate.eq.3)) then
             if (abs(one-vfrac_fluid_sum).le.VOFTOL) then
              ! do nothing
             else
              print *,"vfrac_fluid_sum invalid"
              stop
             endif
            else if (tessellate.eq.1) then
             ! do nothing
            else
             print *,"tessellate invalid"
             stop
            endif

            if ((vfrac_solid_sum.le.one+VOFTOL).and. &
                (vfrac_solid_sum.ge.zero)) then
             ! do nothing
            else
             print *,"vfrac_solid_sum invalid"
             stop
            endif

            if (nmat_fluid+nmat_rigid.ne.nmat) then
             print *,"nmat_fluid or nmat_rigid invalid"
             stop
            endif
            num_processed_solid=0
            num_processed_fluid=0

            if ((tessellate.eq.1).or. &
                (tessellate.eq.3)) then

             loop_counter=0
             do while ((loop_counter.lt.nmat_rigid).and. &
                       (num_processed_solid.lt.nmat_rigid).and. &
                       (uncaptured_volume_fraction.gt. &
                        one-vfrac_solid_sum)) 

               ! F,CEN,ORDER,SLOPE,INTERCEPT
              do im=1,nmat
               if (is_rigid_local(im).eq.1) then
                vofcomp=(im-1)*ngeom_recon+1
                testflag=NINT(mofdataproject(vofcomp+SDIM+1))

                if (testflag.eq.1) then

                 num_processed_solid=num_processed_solid+1

                 uncaptured_volume_fraction= &
                  uncaptured_volume_fraction-vcenter_thin(im)
                 if (uncaptured_volume_fraction.lt.FACETOL_DVOL) then
                  uncaptured_volume_fraction=zero
                 endif

                  ! vcenter_thin(im)=multi_volume(im)/total_vol
                 if ((vcenter_thin(im).gt.FACETOL_DVOL).and. &
                     (uncaptured_volume_fraction.gt.zero)) then

                  do dir2=1,SDIM
                   slope(dir2)=mofdataproject(vofcomp+SDIM+1+dir2)
                  enddo
                  intercept=mofdataproject(vofcomp+2*SDIM+2)

                  ! dist=intercept+n dot (x-x0)
                  ! n points into im
                  ! perturb interface into the other materials
                  ! for shifting interface in direction of normal:
                  ! x0'=x0-eps n
                  ! dist=intercept+n dot (x-x0')=intercept+eps+n dot(x-x0)
                  ! for shifting in direction of n':
                  ! dist=intercept+n dot (x-x0)+eps n dot n'
                  if (slope(dir).eq.zero) then

                   mofdataproject(vofcomp+2*SDIM+2)= &
                    intercept+half*FACETOL_DVOL*dx(dir)

                    ! rigid material
                    ! in: FORT_FACEINIT
                   call multi_get_volume_grid_simple( &
                    local_tessellate, &  !=0,1, or 2
                    bfact,dx,xsten,nhalf, &
                    mofdataproject, &
                    xsten_thin,nhalf_thin, &
                    multi_volume_offset, &
                    multi_cen_offset, &
                    geom_xtetlist_uncapt(1,1,1,tid+1), &
                    nmax, &
                    nmax, &
                    nmat,SDIM,3)

                   mofdataproject(vofcomp+2*SDIM+2)=intercept

                   if (multi_volume_offset(im).gt.multi_volume(im)) then

                    if (multi_area(im).gt.zero) then

                     if (dxthin.gt.zero) then
                      multi_area_line(im)=multi_area(im)/dxthin
                     else
                      print *,"dxthin invalid"
                      stop
                     endif

                     do im_opp=1,nmat
                      local_linefrac(im_opp)=zero
                     enddo
                     total_line=zero
                     do im_opp=1,nmat
                      if (im_opp.ne.im) then
                       call get_iten(im,im_opp,iten,nmat)
                       if (is_processed(iten).eq.0) then
                        is_processed(iten)=1
                        F1=multi_volume(im_opp)
                        F2=multi_volume_offset(im_opp)
                        do dir2=1,SDIM
                         X1(dir2)=multi_cen(dir2,im)
                         X2(dir2)=multi_cen_offset(dir2,im)
                        enddo
                        if (F2.ne.F1) then
                         FSTRIP=F2-F1
                         do dir2=1,SDIM
                          XSTRIP(dir2)=(F2*X2(dir2)-F1*X1(dir2))/FSTRIP
                         enddo
                        else if (F2.eq.F1) then
                         do dir2=1,SDIM
                          XSTRIP(dir2)=zero
                         enddo
                        else
                         print *,"F2 and F1 are invalid"
                         stop
                        endif
                        if (side.eq.1) then
                         XSTRIP(dir)=xsten(-1,dir)
                        else if (side.eq.2) then
                         XSTRIP(dir)=xsten(1,dir)
                        else
                         print *,"side invalid"
                         stop
                        endif
                        do dir2=1,SDIM
                         localface_line(im,im_opp,1+dir2,dir,side)=XSTRIP(dir2)
                        enddo
   
                        local_linefrac(im_opp)=abs(F1-F2)
                        total_line=total_line+local_linefrac(im_opp)
                       else if (is_processed(iten).eq.1) then
                        ! do nothing
                       else
                        print *,"is_processed invalid"
                        stop
                       endif
                      else if (im_opp.eq.im) then
                       ! do nothing
                      else
                       print *,"im_opp or im bust"
                       stop
                      endif
                     enddo !im_opp=1..nmat
   
                     if (total_line.gt.zero) then
                      do im_opp=1,nmat
                       if (im_opp.ne.im) then
                        local_linefrac(im_opp)=local_linefrac(im_opp)/total_line
                        localface_line(im,im_opp,1,dir,side)= &
                         local_linefrac(im_opp)
                       else if (im_opp.eq.im) then
                        ! do nothing
                       else
                        print *,"im_opp or im bust"
                        stop
                       endif
                      enddo ! im_opp

                      totalface_line(im,dir,side)=multi_area_line(im)
                     else
                      print *,"total_line bust: total_face=",total_line
                      stop
                     endif
                    else
                     print *,"im boundary disappeared 1"
                     print *,"loop_counter=",loop_counter
                     print *,"num_processed_solid=",num_processed_solid
                     print *,"uncaptured_volume_fraction ", &
                       uncaptured_volume_fraction
                     print *,"vfrac_solid_sum ",vfrac_solid_sum
                     print *,"im=",im
                     print *,"multi_volume(im)=",multi_volume(im)
                     print *,"multi_volume_offset(im)=",multi_volume_offset(im)
                     print *,"multi_area(im)=",multi_area(im)
                     stop
                    endif

                   else
                    print *,"im region should grow 1"
                    stop
                   endif

                  else if ((abs(slope(dir)-one).le.VOFTOL).or. &
                           (abs(slope(dir)+one).le.VOFTOL)) then
                   ! do nothing
                  else
                   print *,"slope invalid"
                   stop
                  endif

                 else if ((uncaptured_volume_fraction.eq.zero).or. &
                          (vcenter_thin(im).le.FACETOL_DVOL)) then
                  ! do nothing
                 else
                  print *,"uncaptured_volume_fraction or vcenter_thin bad"
                  stop 
                 endif 

                else if (testflag.eq.0) then
                 ! do nothing
                else
                 print *,"testflag invalid"
                 stop
                endif  
               else if (is_rigid_local(im).eq.0) then
                ! do nothing
               else
                print *,"is_rigid_local invalid"
                stop
               endif
              enddo ! im=1..nmat
              loop_counter=loop_counter+1
             enddo  ! while 
                    ! loop_counter<nmat_rigid and
                    ! num_processed_solid<nmat_rigid and
                    ! uncaptured_volume_fraction>1-vfrac_solid_sum

            else if (tessellate.eq.0) then
             ! do nothing
            else
             print *,"tessellate invalid47"
             stop
            endif


            loop_counter=0
            do while ((loop_counter.lt.nmat_fluid).and. &
                      (num_processed_fluid.lt.nmat_fluid).and. &
                      (uncaptured_volume_fraction.gt.zero))

              ! F,CEN,ORDER,SLOPE,INTERCEPT
             do im=1,nmat
              if (is_rigid_local(im).eq.0) then
               vofcomp=(im-1)*ngeom_recon+1
               testflag=NINT(mofdataproject(vofcomp+SDIM+1))

               if (testflag.eq.num_processed_fluid+1) then

                num_processed_fluid=num_processed_fluid+1

                uncaptured_volume_fraction= &
                  uncaptured_volume_fraction-vcenter_thin(im)
                if (uncaptured_volume_fraction.lt.FACETOL_DVOL) then
                 uncaptured_volume_fraction=zero
                endif

                  ! vcenter_thin(im)=multi_volume(im)/total_vol
                if ((vcenter_thin(im).gt.FACETOL_DVOL).and. &
                    (uncaptured_volume_fraction.gt.zero)) then

                 do dir2=1,SDIM
                  slope(dir2)=mofdataproject(vofcomp+SDIM+1+dir2)
                 enddo
                 intercept=mofdataproject(vofcomp+2*SDIM+2)

                 ! dist=intercept+n dot (x-x0)
                 ! n points into im
                 ! perturb interface into the other materials
                 ! for shifting interface in direction of normal:
                 ! x0'=x0-eps n
                 ! dist=intercept+n dot (x-x0')=intercept+eps+n dot(x-x0)
                 ! for shifting in direction of n':
                 ! dist=intercept+n dot (x-x0)+eps n dot n'
                 if (slope(dir).eq.zero) then

                  mofdataproject(vofcomp+2*SDIM+2)= &
                    intercept+half*FACETOL_DVOL*dx(dir)

                   ! fluid material
                   ! in: FORT_FACEINIT
                  call multi_get_volume_grid_simple( &
                   local_tessellate, &  ! 0,1, or 2
                   bfact,dx,xsten,nhalf, &
                   mofdataproject, &
                   xsten_thin,nhalf_thin, &
                   multi_volume_offset, &
                   multi_cen_offset, &
                   geom_xtetlist_uncapt(1,1,1,tid+1), &
                   nmax, &
                   nmax, &
                   nmat,SDIM,3)

                  mofdataproject(vofcomp+2*SDIM+2)=intercept

                   ! vcenter_thin(im)=multi_volume(im)/total_vol>FACETOL_DVOL
                  if (multi_volume_offset(im).gt.multi_volume(im)) then

                   if (multi_area(im).gt.zero) then

                    if (dxthin.gt.zero) then
                     multi_area_line(im)=multi_area(im)/dxthin
                    else
                     print *,"dxthin invalid"
                     stop
                    endif

                    do im_opp=1,nmat
                     local_linefrac(im_opp)=zero
                    enddo
                    total_line=zero

                    do im_opp=1,nmat

                     if ((is_rigid_local(im_opp).eq.0).or. &
                         (local_tessellate.eq.1)) then

                      if (im_opp.ne.im) then
                       call get_iten(im,im_opp,iten,nmat)
                       if (is_processed(iten).eq.0) then
                        is_processed(iten)=1
                        F1=multi_volume(im_opp)
                        F2=multi_volume_offset(im_opp)
                        do dir2=1,SDIM
                         X1(dir2)=multi_cen(dir2,im)
                         X2(dir2)=multi_cen_offset(dir2,im)
                        enddo
                        if (F2.ne.F1) then
                         FSTRIP=F2-F1
                         do dir2=1,SDIM
                          XSTRIP(dir2)=(F2*X2(dir2)-F1*X1(dir2))/FSTRIP
                         enddo
                        else if (F2.eq.F1) then
                         do dir2=1,SDIM
                          XSTRIP(dir2)=zero
                         enddo
                        else
                         print *,"F2 and F1 are invalid"
                         stop
                        endif
                        if (side.eq.1) then
                         XSTRIP(dir)=xsten(-1,dir)
                        else if (side.eq.2) then
                         XSTRIP(dir)=xsten(1,dir)
                        else
                         print *,"side invalid"
                         stop
                        endif
                        do dir2=1,SDIM
                         localface_line(im,im_opp,1+dir2,dir,side)=XSTRIP(dir2)
                        enddo

                        local_linefrac(im_opp)=abs(F1-F2)
                        total_line=total_line+local_linefrac(im_opp)
                       else if (is_processed(iten).eq.1) then
                        ! do nothing
                       else
                        print *,"is_processed invalid"
                        stop
                       endif
                      else if (im_opp.eq.im) then
                       ! do nothing
                      else
                       print *,"im_opp or im bust"
                       stop
                      endif
                     else if ((is_rigid_local(im_opp).eq.1).and. &
                              (local_tessellate.eq.0)) then
                      ! do nothing
                     else
                      print *,"is_rigid_local or local_tessellate invalid48"
                      stop
                     endif
                    enddo !im_opp

                    if (total_line.gt.zero) then
                     do im_opp=1,nmat
                      if (im_opp.ne.im) then
                       local_linefrac(im_opp)=local_linefrac(im_opp)/total_line
                       localface_line(im,im_opp,1,dir,side)= &
                         local_linefrac(im_opp)
                      else if (im_opp.eq.im) then
                       ! do nothing
                      else
                       print *,"im_opp or im bust"
                       stop
                      endif
                     enddo ! im_opp

                     totalface_line(im,dir,side)=multi_area_line(im)
                    else
                     print *,"total_line bust: total_face=",total_line
                     stop
                    endif
                   else
                    print *,"im boundary disappeared 2"
                    print *,"loop_counter=",loop_counter
                    print *,"num_processed_solid=",num_processed_solid
                    print *,"num_processed_fluid=",num_processed_fluid
                    print *,"uncaptured_volume_fraction ", &
                       uncaptured_volume_fraction
                    print *,"vfrac_solid_sum ",vfrac_solid_sum
                    print *,"im=",im
                    print *,"multi_volume(im)=",multi_volume(im)
                    print *,"multi_volume_offset(im)=",multi_volume_offset(im)
                    print *,"multi_area(im)=",multi_area(im)
                    stop
                   endif

                  else
                   print *,"im region should grow 2"
                   stop
                  endif

                 else if ((abs(slope(dir)-one).le.VOFTOL).or. &
                          (abs(slope(dir)+one).le.VOFTOL)) then
                  ! do nothing
                 else
                  print *,"slope invalid"
                  stop
                 endif

                else if ((uncaptured_volume_fraction.eq.zero).or. &
                         (vcenter_thin(im).le.FACETOL_DVOL)) then
                 ! do nothing
                else
                 print *,"uncaptured_volume or vcenter_thin invalid"
                 stop 
                endif 

               endif  ! testflag==num_processed_fluid+1
           
              else if (is_rigid_local(im).eq.1) then
               ! do nothing
              else
               print *,"is_rigid_local invalid"
               stop
              endif

             enddo ! im=1..nmat

             loop_counter=loop_counter+1
            enddo  ! while 
                   ! loop_counter<nmat_fluid and
                   ! num_processed_fluid<nmat_fluid and
                   ! uncaptured_volume_fraction>0

           else
            print *,"im_crit_thin out of range"
            stop
           endif

          else if (nface_decomp.eq.0) then
           ! do nothing
          else
           print *,"nface_decomp invalid"
           stop
          endif

         enddo
         enddo ! dir,side
    
        else
         print *,"im_crit out of range"
         stop
        endif

        dir2=1
        do im=1,nmat
         do dir=1,SDIM
          do side=1,2
           localface(im,dir,side,dir2)= &
            localface(im,dir,side,dir2)*totalface(dir,side,1)
          enddo
         enddo
        enddo 

         ! centroid in absolute coordinate system.
         ! centroid is 2,...,sdim+1 components
        iface=0
        do im=1,nmat
         do dir=1,SDIM
          do side=1,2
           do dir2=1,SDIM+1
            iface=iface+1
            facefab(D_DECL(i,j,k),iface)=localface(im,dir,side,dir2)
           enddo ! dir2
          enddo ! side
         enddo ! dir
        enddo ! im

        if (iface.ne.nface) then
         print *,"iface invalid"
         stop
        endif

        iface=0

        if (nface_decomp.eq.1) then

         dir2=1
         do im1=1,nmat
          do im2=1,nmat
           do dir=1,SDIM
            do side=1,2
             localface_line(im1,im2,dir2,dir,side)= &
              localface_line(im1,im2,dir2,dir,side)* &
              totalface_line(im1,dir,side)
            enddo ! side
           enddo ! dir
          enddo ! im2
         enddo ! im1

         do im1=1,nmat
          do im2=1,nmat
           do dir2=1,SDIM+1
            do dir=1,SDIM
             do side=1,2
              iface=iface+1
              facefab(D_DECL(i,j,k),nface+iface)= &
                localface_line(im1,im2,dir2,dir,side)
             enddo ! side
            enddo ! dir
           enddo ! dir2
          enddo ! im2
         enddo ! im1
 
        else if (nface_decomp.eq.0) then
         ! do nothing
        else
         print *,"nface_decomp invalid"
         stop
        endif

        if (iface.ne.nface_decomp) then
         print *,"iface invalid"
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
      end subroutine FORT_FACEINIT


      subroutine FORT_FACEINITTEST( &
       tid, &
       tessellate, &
       level, &
       finest_level, &
       facefab,DIMS(facefab), &
       facetest,DIMS(facetest), &
       maskfab,DIMS(maskfab), &
       vofrecon,DIMS(vofrecon), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       rz_flag, &
       xlo,dx, &
       time, &
       ngrow, &
       nmat, &
       nface)

      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module
      use mof_redist_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: tessellate
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nface
      INTEGER_T, intent(in) :: nmat,ngrow
      INTEGER_T, intent(in) :: DIMDEC(facefab)
      INTEGER_T, intent(in) :: DIMDEC(facetest)
      INTEGER_T, intent(in) :: DIMDEC(maskfab)
      INTEGER_T, intent(in) :: DIMDEC(vofrecon)

      REAL_T, intent(in) :: facefab(DIMV(facefab),nface)
      REAL_T, intent(out) :: facetest(DIMV(facetest),nmat*SDIM)
      REAL_T, intent(in) :: maskfab(DIMV(maskfab),2)
      REAL_T, intent(in) :: vofrecon(DIMV(vofrecon),nmat*ngeom_recon)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: rz_flag
      REAL_T, intent(in) :: time

      INTEGER_T i,j,k
      INTEGER_T icell,jcell,kcell
      INTEGER_T ii,jj,kk
      INTEGER_T im
      INTEGER_T dir
      INTEGER_T side
      INTEGER_T side_cell
      INTEGER_T iface
      REAL_T total_face
      REAL_T facefrac(nmat)
      REAL_T faceleft(nmat)
      REAL_T faceright(nmat)
      REAL_T local_face_test(nmat)
      INTEGER_T nface_test
      REAL_T xstenMAC(-1:1,SDIM)
      INTEGER_T nhalf
      INTEGER_T at_RZ_face
 
      nhalf=1
 
      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if ((tessellate.ne.0).and.(tessellate.ne.1)) then
       print *,"tessellate invalid49"
       stop
      endif
 
      if (bfact.lt.1) then
       print *,"bfact invalid146"
       stop
      endif

       ! (nmat,sdim,2,sdim+1)
      nface_test=nmat*SDIM*2*(1+SDIM)
      if (nface_test.ne.nface) then
       print *,"nface invalid faceinittest nface nface_test ",nface,nface_test
       stop
      endif

      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid in faceinittest"
       stop
      endif

      if (ngrow.lt.1) then
       print *,"ngrow<1 error in faceinittest"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(facetest),ngrow,-1,2886)
      call checkbound(fablo,fabhi,DIMS(facefab),ngrow,-1,2887)
      call checkbound(fablo,fabhi,DIMS(maskfab),ngrow,-1,2888)
      call checkbound(fablo,fabhi,DIMS(vofrecon),ngrow,-1,2889)
      
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid faceinittest"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid faceinittest"
       print *,"ngeom_raw=",ngeom_raw
       stop
      endif

      if (rz_flag.eq.0) then
       ! do nothing
      else if (rz_flag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.3) then 
       ! do nothing
      else
       print *,"rz_flag invalid in faceinittest"
       stop
      endif

       ! facetest is initialized to zero in NavierStokes::makeFaceTest.
      do dir=1,SDIM
 
       if ((level.ge.0).and.(level.le.finest_level)) then

        ii=0
        jj=0
        kk=0
        if (dir.eq.1) then
         ii=1
        else if (dir.eq.2) then
         jj=1
        else if ((dir.eq.3).and.(SDIM.eq.3)) then
         kk=1
        else
         print *,"dir invalid in faceinittest"
         stop
        endif

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
         growlo,growhi,ngrow-1,dir-1) 

        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir)

         at_RZ_face=0
         if (levelrz.eq.0) then
          ! do nothing
         else if (levelrz.eq.1) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          if (dir.eq.1) then
           if (abs(xstenMAC(0,1)).le.VOFTOL*dx(1)) then
            at_RZ_face=1
           endif
          endif
         else if (levelrz.eq.3) then
          if (dir.eq.1) then
           if (abs(xstenMAC(0,1)).le.VOFTOL*dx(1)) then
            at_RZ_face=1
           endif
          endif
         else
          print *,"levelrz invalid"
          stop
         endif

         if (at_RZ_face.eq.1) then
          do im=1,nmat
           local_face_test(im)=1
          enddo
         else if (at_RZ_face.eq.0) then

          do side=1,2

           if (side.eq.1) then
            icell=i-ii
            jcell=j-jj
            kcell=k-kk
            side_cell=2
           else if (side.eq.2) then
            icell=i
            jcell=j
            kcell=k
            side_cell=1
           else
            print *,"side invalid"
            stop
           endif

           total_face=zero
           do im=1,nmat
            ! (im,dir,side,dir2)
            iface=(im-1)*SDIM*2*(SDIM+1)+(dir-1)*2*(SDIM+1)+ &
                  (side_cell-1)*(SDIM+1)+1
            facefrac(im)=facefab(D_DECL(icell,jcell,kcell),iface)
            if (tessellate.eq.1) then
             total_face=total_face+facefrac(im)
            else if (tessellate.eq.0) then
             if (is_rigid(nmat,im).eq.0) then
              total_face=total_face+facefrac(im)
             else if (is_rigid(nmat,im).eq.1) then
              ! do nothing
             else
              print *,"is_rigid invalid"
              stop
             endif
            else
             print *,"tessellate invalid50"
             stop
            endif
           enddo !im=1..nmat 

           if (total_face.le.zero) then
            print *,"total_face invalid"
            stop
           endif
           do im=1,nmat
            if (side.eq.1) then
             faceleft(im)=facefrac(im)/total_face
            else if (side.eq.2) then
             faceright(im)=facefrac(im)/total_face
            else
             print *,"side invalid"
             stop
            endif
           enddo ! im
          enddo ! side
        
          do im=1,nmat
           if (abs(faceleft(im)-faceright(im)).le.FACETOL_REDIST) then
            local_face_test(im)=1
           else if (abs(faceleft(im)-faceright(im)).ge.FACETOL_REDIST) then
            local_face_test(im)=0
           else
            print *,"faceleft or faceright bust"
            stop
           endif
          enddo ! im=1..nmat
         else 
          print *,"at_RZ_face invalid"
          stop
         endif

         do im=1,nmat
          facetest(D_DECL(i,j,k),nmat*(dir-1)+im)=local_face_test(im)
         enddo ! im

        enddo ! k
        enddo ! j
        enddo ! i

       else
        print *,"level invalid faceinittest 2"
        stop
       endif

      enddo ! dir

      return
      end subroutine FORT_FACEINITTEST

       ! 1. NavierStokes::makeFaceFrac
       !      -> FORT_FACEINIT
       ! 2. NavierStokes::ProcessFaceFrac
       !      -> FORT_FACEPROCESS
       ! called from: NavierStokes::ProcessFaceFrac (NavierStokes.cpp)
       ! facefab is initialized in FORT_FACEINIT
       ! centroids in facefab are in an absolute coordinate system.
      subroutine FORT_FACEPROCESS( &
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
       nmat,nface_src,nface_dst)

      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: ngrow_source
      INTEGER_T, intent(in) :: ngrow_dest
      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: dir
      INTEGER_T, intent(in) :: tessellate ! 0,1, or 3
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nface_src,nface_dst
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: DIMDEC(dstfab)
      INTEGER_T, intent(in) :: DIMDEC(facefab)
      INTEGER_T, intent(in) :: DIMDEC(vofrecon)

      REAL_T, intent(out) :: dstfab(DIMV(dstfab),nface_dst)
      REAL_T, intent(in) :: facefab(DIMV(facefab),nface_src)
      REAL_T, intent(in) :: vofrecon(DIMV(vofrecon),nmat*ngeom_recon)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: rz_flag
      REAL_T, intent(in) :: time

      INTEGER_T i,j,k
      INTEGER_T iface_left,iface_right
      INTEGER_T inormal
      INTEGER_T left_face_ok,right_face_ok
      INTEGER_T side
      INTEGER_T ii,jj,kk
      INTEGER_T im,im_opp

      INTEGER_T nface_src_test,nface_dst_test,vofcomp
      REAL_T mofdata_left(nmat*ngeom_recon)
      REAL_T mofdata_right(nmat*ngeom_recon)
      REAL_T xsten_left(-3:3,SDIM)
      REAL_T xsten_right(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T dir2
      REAL_T multi_cen_left(SDIM,nmat)
      REAL_T multi_cen_right(SDIM,nmat)
      REAL_T cencell_left(SDIM)
      REAL_T cencell_right(SDIM)
      REAL_T volcell_left
      REAL_T volcell_right
      REAL_T leftface(nmat,SDIM+1)
      REAL_T rightface(nmat,SDIM+1)
      REAL_T frac_left(nmat)
      REAL_T frac_right(nmat)
      REAL_T left_total
      REAL_T right_total
      REAL_T vol_total
      REAL_T x_left(SDIM,nmat)
      REAL_T x_right(SDIM,nmat)
      INTEGER_T ml,mr
      REAL_T frac_pair(nmat,nmat)  !(m_left,m_right)
      REAL_T dist_pair(nmat,nmat)
      REAL_T x_pair(nmat,nmat,SDIM)  ! (m_left,m_right)
      REAL_T dpair
      REAL_T delta,RR
      REAL_T xstenMAC(-1:1,SDIM)
      INTEGER_T nhalfMAC
      INTEGER_T at_RZ_face
      REAL_T L_face
      INTEGER_T nmax
      INTEGER_T local_tessellate
      INTEGER_T caller_id
 
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

      nface_src_test=nmat*SDIM*2*(1+SDIM)
      if (nface_src_test.ne.nface_src) then
       print *,"nface_src invalid nface_src nface_src_test ", &
         nface_src,nface_src_test
       stop
      endif
      nface_dst_test=nmat*nmat*2
      if (nface_dst_test.ne.nface_dst) then
       print *,"nface_dst invalid nface_dst nface_dst_test ", &
         nface_dst,nface_dst_test
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in faceprocess"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(dstfab),ngrow_dest,dir,2880)
      call checkbound(fablo,fabhi,DIMS(facefab),ngrow_source,-1,2881)
      call checkbound(fablo,fabhi,DIMS(vofrecon),ngrow_source,-1,2882)
      
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
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

      if (rz_flag.eq.0) then
       ! do nothing
      else if (rz_flag.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
      else if (rz_flag.eq.3) then 
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

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalfMAC,dir+1)

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
       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (dir.eq.0) then
         if (abs(xstenMAC(0,1)).le.VOFTOL*dx(1)) then
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
       else if (levelrz.eq.3) then
        if (dir.eq.0) then
         if (abs(xstenMAC(0,1)).le.VOFTOL*dx(1)) then
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
         ! (im,dir,side,dir2)
        do im=1,nmat
         do dir2=1,SDIM+1
          side=2
          call abs_array_index4(im,dir+1,side,dir2, &
           nmat,SDIM,2,SDIM+1,iface_left)
          side=1
          call abs_array_index4(im,dir+1,side,dir2, &
           nmat,SDIM,2,SDIM+1,iface_right)

          if (left_face_ok.eq.1) then

           leftface(im,dir2)=facefab(D_DECL(i-ii,j-jj,k-kk),iface_left)
           if (right_face_ok.eq.1) then
            ! do nothing
           else if (right_face_ok.eq.0) then
            rightface(im,dir2)=leftface(im,dir2)
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

           rightface(im,dir2)=facefab(D_DECL(i,j,k),iface_right)
           if (left_face_ok.eq.1) then
            ! do nothing
           else if (left_face_ok.eq.0) then
            leftface(im,dir2)=rightface(im,dir2)
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

         enddo ! dir2

         frac_left(im)=leftface(im,1)
         frac_right(im)=rightface(im,1)

         if ((tessellate.eq.1).or. &
             (tessellate.eq.3).or. &
             (is_rigid(nmat,im).eq.0)) then
          left_total=left_total+frac_left(im)
          right_total=right_total+frac_right(im)
         else if ((tessellate.eq.0).and.(is_rigid(nmat,im).eq.1)) then
          ! do nothing
         else 
          print *,"tessellate or is_rigid invalid"
          stop
         endif 

          ! absolute coordinate system.
         do dir2=1,SDIM
          x_left(dir2,im)=leftface(im,dir2+1)
          x_right(dir2,im)=rightface(im,dir2+1)
         enddo
        enddo ! im=1..nmat

        if ((left_total.le.zero).or.(right_total.le.zero)) then
         print *,"left_total or right_total invalid"
         stop
        endif
        do im=1,nmat
         frac_left(im)=frac_left(im)/left_total
         frac_right(im)=frac_right(im)/right_total
        enddo ! im

        call CISBOX(xsten_left,nhalf, &
         xlo,dx,i-ii,j-jj,k-kk, &
         bfact,level, &
         volcell_left,cencell_left,SDIM) 

        call CISBOX(xsten_right,nhalf, &
         xlo,dx,i,j,k, &
         bfact,level, &
         volcell_right,cencell_right,SDIM) 

        do im = 1,nmat

         if ((tessellate.eq.1).or. &
             (tessellate.eq.3).or. &
             (is_rigid(nmat,im).eq.0)) then

          if ((frac_left(im).gt.one+FACETOL_SANITY).or. &
              (frac_left(im).lt.zero).or. &
              (frac_right(im).gt.one+FACETOL_SANITY).or. &
              (frac_right(im).lt.zero)) then
           print *,"frac_left or frac_right out of range"
           stop
          endif
          if (frac_left(im).le.FACETOL_DVOL) then
           frac_left(im) = zero
          else if (frac_left(im).ge.one-FACETOL_DVOL) then
           frac_left(im) = one
          endif
          if (frac_right(im).lt.FACETOL_DVOL) then
           frac_right(im) = zero
          else if (frac_right(im).gt.one-FACETOL_DVOL)then
           frac_right(im) = one
          endif

         else if ((tessellate.eq.0).and.(is_rigid(nmat,im).eq.1)) then
          ! do nothing
         else
          print *,"tessellate or is_rigid invalid"
          stop
         endif 

        enddo ! im=1..nmat

        if ((left_face_ok.eq.1).and. &  ! cell to left of face is in FAB
            (right_face_ok.eq.1)) then  ! cell to right of face is in FAB

         do im=1,nmat*ngeom_recon
          mofdata_left(im)=vofrecon(D_DECL(i-ii,j-jj,k-kk),im)
          mofdata_right(im)=vofrecon(D_DECL(i,j,k),im)
         enddo

         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          if ((tessellate.eq.1).or. &
              (tessellate.eq.3).or. &
              (is_rigid(nmat,im).eq.0)) then
           ! do nothing
          else if ((tessellate.eq.0).and. &
                   (is_rigid(nmat,im).eq.1)) then
           do dir2=1,ngeom_recon
            mofdata_left(vofcomp+dir2-1)=zero 
            mofdata_right(vofcomp+dir2-1)=zero 
           enddo
          else
           print *,"tessellate or is_rigid invalid"
           stop
          endif

          do dir2=1,SDIM
           multi_cen_left(dir2,im)=mofdata_left(vofcomp+dir2)
           multi_cen_right(dir2,im)=mofdata_right(vofcomp+dir2)
          enddo

         enddo ! im=1..nmat

           ! x_pair in absolute coordinate system.
         caller_id=12
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
           nmat, &
           dir+1, &
           frac_pair, & ! left,right
           x_pair, & ! left,right (unused for now)
           SDIM, &
           geom_xtetlist(1,1,1,tid+1), &
           nmax, &
           geom_xtetlist_old(1,1,1,tid+1), &
           nmax, &
           nmax, &
           caller_id)

        else if ((left_face_ok.eq.0).or. &
                 (right_face_ok.eq.0)) then

         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          do dir2=1,SDIM
           multi_cen_left(dir2,im)=zero
           multi_cen_right(dir2,im)=zero
          enddo
         enddo ! im

         do ml = 1, nmat
          do mr = 1, nmat
           frac_pair(ml,mr)=zero
          enddo
         enddo
         do ml=1,nmat
          frac_pair(ml,ml)=frac_left(ml)
         enddo

        else
         print *,"left or right face ok invalid"
         stop
        endif

        vol_total=zero
        do im=1,nmat
        do im_opp=1,nmat
         vol_total=vol_total+frac_pair(im,im_opp)
        enddo
        enddo

        if (vol_total.gt.zero) then

         do im=1,nmat
         do im_opp=1,nmat
          frac_pair(im,im_opp)=frac_pair(im,im_opp)/vol_total
         enddo
         enddo

        else
         print *,"vol_total invalid"
         stop
        endif

        delta=xsten_right(0,dir+1)-xsten_left(0,dir+1)
        if (delta.le.zero) then
         print *,"delta invalid faceprocess"
         stop
        endif 

        do ml = 1, nmat

         if ((tessellate.eq.1).or. &
             (tessellate.eq.3).or. &
             (is_rigid(nmat,ml).eq.0)) then

          do mr = 1, nmat

           if ((tessellate.eq.1).or. &
               (tessellate.eq.3).or. &
               (is_rigid(nmat,mr).eq.0)) then

            if ((frac_pair(ml,mr).lt.-FACETOL_SANITY).or. &
                (frac_pair(ml,mr).gt.one+FACETOL_SANITY)) then
             print *,"frac_pair invalid"
             stop
            endif
  
            if (frac_pair(ml,mr).lt.FACETOL_DVOL) then
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
   
             if ((levelrz.eq.0).or. &
                 (levelrz.eq.1)) then !XY or RZ
              RR=one 
             else if (levelrz.eq.3) then ! R-theta
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

           else if ((tessellate.eq.0).and.(is_rigid(nmat,mr).eq.1)) then
            ! do nothing
           else
            print *,"tessellate or is_rigid invalid"
            stop
           endif 

          enddo ! mr=1..nmat

         else if ((tessellate.eq.0).and.(is_rigid(nmat,ml).eq.1)) then
          ! do nothing
         else
          print *,"tessellate or is_rigid invalid"
          stop
         endif 

        enddo ! ml=1..nmat

       else if (at_RZ_face.eq.1) then

        do ml = 1, nmat
         do mr = 1, nmat
          frac_pair(ml,mr)=zero
          dist_pair(ml,mr)=zero
         enddo
        enddo

       else
        print *,"at_RZ_face invalid"
        stop
       endif

       iface_left=0
       do ml = 1, nmat
        do mr = 1, nmat
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
      end subroutine FORT_FACEPROCESS



#if (STANDALONE==1)
      end module mof_redist_cpp_module
#endif

#undef STANDALONE

