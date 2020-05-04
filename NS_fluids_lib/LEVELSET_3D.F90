#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "LEVEL_F.H"

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

        module levelset_module
        use probf90_module

        contains


        INTEGER_T function is_default(contactangle)
        IMPLICIT NONE

        REAL_T contactangle

        if ((contactangle.lt.zero).or.(contactangle.gt.Pi)) then
         print *,"contactangle invalid"
         stop
        else if (abs(contactangle-half*Pi).le.0.1) then
         is_default=1
        else 
         is_default=0
        endif

        return
        end function is_default
 

      subroutine initpforce( &
        bfact,dx,xsten0, &
        RD,RDx,RD_HEIGHT, &
        time, &
        lsdata, &
        dir,nmat, &
        im,im_opp, &
        pforce)
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T bfact
      INTEGER_T RD,RDx,RD_HEIGHT
      REAL_T dx(SDIM)
      REAL_T xsten0(-RDx:RDx,SDIM)
      REAL_T time
      INTEGER_T nmat,im,im_opp
      INTEGER_T dir
      REAL_T lsdata( &
       D_DECL(-RD:RD,-RD:RD,-RD:RD),nmat)

      REAL_T pforce

      INTEGER_T is_pforce_probtype
      REAL_T columnLS(-RD:RD)
      INTEGER_T icolumn
      REAL_T n1d,col_ht,xforce
      INTEGER_T crossing_status

      if ((dir.lt.0).or.(dir.ge.SDIM)) then
       print *,"dir invalid in initpforce, dir=",dir
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid84"
       stop
      endif

       ! Sep 16, 2018: returns 1 if probtype.eq.90 and 2D
       !               returns 0 otherwise
      call pforce_probtype_flag_init(is_pforce_probtype)

      if (is_pforce_probtype.eq.1) then

       if (RD.lt.4) then
        print *,"RD invalid in initpforce RD=",RD
        stop
       endif
       if (RD_HEIGHT.lt.3) then
        print *,"RD_HEIGHT invalid ",RD_HEIGHT
        stop
       endif
       if (RD_HEIGHT.lt.RD-1) then
        print *,"RD_HEIGHT<RD-1 RD_HEIGHT= ",RD_HEIGHT
        stop
       endif
       if (RDx.ne.2*RD+1) then
        print *,"RDx invalid"
        stop
       endif

       if (nmat.ne.2) then
        print *,"this routine not designed for nmat<>2 yet"
        print *,"if nmat>2, then extrapolate_LS must be called"
        stop
       endif
       if (nmat.ne.num_materials) then
        print *,"nmat invalid"
        stop
       endif
       if (im.ge.im_opp) then
        print *,"im and im_opp invalid"
        stop
       endif
       if ((im.lt.1).or.(im.gt.nmat)) then
        print *,"im invalid30"
        stop
       endif
       if ((im_opp.lt.1).or.(im_opp.gt.nmat)) then
        print *,"im_opp invalid in init pforce im_opp=",im_opp
        stop
       endif

       do icolumn=-RD_HEIGHT,RD_HEIGHT

        if (dir.eq.0) then
         columnLS(icolumn)=lsdata(D_DECL(icolumn,0,0),im)- &
                           lsdata(D_DECL(icolumn,0,0),im_opp)
        else if (dir.eq.1) then
         columnLS(icolumn)=lsdata(D_DECL(0,icolumn,0),im)- &
                           lsdata(D_DECL(0,icolumn,0),im_opp)
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         columnLS(icolumn)=lsdata(D_DECL(0,0,icolumn),im)- &
                           lsdata(D_DECL(0,0,icolumn),im_opp)
        else
         print *,"dir invalid initpforce"
         stop
        endif

       enddo !icolumn

       if (columnLS(1).ge.columnLS(-1)) then
        n1d=one ! "im" material on top
       else
        n1d=-one  ! "im" material on bottom
       endif

       call get_col_ht_LS( &
        crossing_status, &
        bfact,dx,xsten0, &
        RD,RDx,RD_HEIGHT, &
        columnLS,col_ht, &
        dir+1,n1d, &
        nmat, &
        SDIM)

       if (crossing_status.eq.1) then
        ! do nothing
       else if (crossing_status.eq.0) then
        ! do nothing
       else
        print *,"crossing_status invalid"
        stop
       endif

       if (dir.eq.0) then
        xforce=col_ht
       else if ((dir.gt.0).and.(dir.lt.SDIM)) then
        xforce=xsten0(0,1)
       else
        print *,"dir invalid initpforce 2"
        stop
       endif

       call pressure_force(xforce,time,pforce)

      else if (is_pforce_probtype.eq.0) then

       pforce=zero

      else
       print *,"is_pforce_problem invalid" 
       stop
      endif

      return
      end subroutine initpforce

       ! called from FORT_CURVSTRIP
      subroutine initheightLS( &
        conservative_tension_force, &
        icenter,jcenter,kcenter, &
        level, &
        finest_level, &
        bfact,dx, &
        xcenter, &
        nrmcenter, & ! sdim x nmat components
        dircrit, &
        side, &
        signside, &
        time, &
        xsten, &
        velsten, &
        mgoni_temp, &
        lssten, &
        nrmsten, &
        vol_sten, &
        area_sten, &
        curvHT, &
        curvFD, &
        mgoni_force, &
        ZEYU_thet_d, &
        ZEYU_u_cl, &
        im3, &
        nmat, &
        visc_coef,nten, &
        im,im_opp,iten, &
        RD,RDx,RD_HEIGHT)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: conservative_tension_force
      INTEGER_T, intent(in) :: icenter,jcenter,kcenter
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xcenter(SDIM)
      INTEGER_T, intent(in) :: dircrit
      INTEGER_T, intent(in) :: side
      INTEGER_T, intent(in) :: signside
      INTEGER_T, intent(in) :: nten
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: im,im_opp
      INTEGER_T, intent(out) :: im3
      INTEGER_T :: imhold
      INTEGER_T :: im3_present_node
      INTEGER_T, intent(in) :: iten
      INTEGER_T :: iten_test
      INTEGER_T, intent(in) :: RD,RDx,RD_HEIGHT
      REAL_T, intent(in) :: visc_coef
      REAL_T user_tension(nten)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: vol_sten
      REAL_T, intent(in) :: area_sten(SDIM,2)
      REAL_T, intent(out) :: curvHT
      REAL_T, intent(out) :: curvFD
      REAL_T, intent(out) :: mgoni_force(SDIM)
 
      INTEGER_T dir2

      REAL_T columnLS(-RD:RD)

      REAL_T lsdata( &
       D_DECL(-RD:RD,-RD:RD,-RD:RD))  

      REAL_T htfunc(-1:1,-1:1)

      REAL_T, intent(in) :: xsten(-RDx:RDx,SDIM)
      REAL_T xsten_curv(-2:2,SDIM)

      REAL_T, intent(in) :: velsten( &
       D_DECL(-1:1,-1:1,-1:1),SDIM)

      REAL_T, intent(in) :: mgoni_temp( &
       D_DECL(-1:1,-1:1,-1:1),nmat)

      REAL_T, intent(in) :: lssten( &
       D_DECL(-RD:RD,-RD:RD,-RD:RD),nmat)

      REAL_T, intent(in) :: nrmsten( &
       D_DECL(-1:1,-1:1,-1:1),SDIM*nmat)

      REAL_T, intent(in) :: nrmcenter(SDIM*nmat)
      REAL_T nrmtest(SDIM*nmat)

      REAL_T ngrid(SDIM)

      INTEGER_T itanlo,itanhi,jtanlo,jtanhi
      INTEGER_T itan,jtan
      INTEGER_T iwidth,jwidth,kheight
      INTEGER_T iofs,jofs,kofs
      INTEGER_T iwidthnew
      INTEGER_T icell,jcell,kcell
      INTEGER_T inode,jnode,knode
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      INTEGER_T i2,j2,k2
      INTEGER_T node_index(3)
      INTEGER_T cell_index(3)

      REAL_T hxR,hxL,hx,hxx
      REAL_T hyR,hyL,hy,hyy
      REAL_T hxy
      REAL_T arclen,arclenx,arcleny,arclenr,g,gx,gy,gr

      REAL_T col_ht

      REAL_T n1d

      INTEGER_T iten_13,iten_23
      REAL_T gamma1,gamma2
      REAL_T cos_angle,sin_angle
      INTEGER_T nten_test

      INTEGER_T klo_sten_short,khi_sten_short
      INTEGER_T klo_sten_ht,khi_sten_ht

      REAL_T dotprod,udotn,totaludotn
      REAL_T totalwt,wt,wtnode
      REAL_T liquid_viscosity
      REAL_T nproject(SDIM)
      INTEGER_T use_DCA
      REAL_T LSTEST(nmat)
      REAL_T LSTEST_EXTEND
      REAL_T LSMAX
      REAL_T mag,mag1,mag2,mag3
 

      REAL_T nperp(SDIM) 
      REAL_T nghost(SDIM) 
      REAL_T nmain(SDIM) 
      REAL_T nopp(SDIM) 
      REAL_T nfluid(SDIM) 
      REAL_T nfluid_cen(SDIM) 
      REAL_T nfluid_def1(SDIM) 
      REAL_T nfluid_def2(SDIM) 
      REAL_T nfluid_def3(SDIM) 
      REAL_T nfluid_save(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T nmain_save(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T nopp_save(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T ngrid_save(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T nsolid(SDIM) 
      REAL_T nsolid_save(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T ncurv1_save(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T ncurv2_save(D_DECL(-1:1,-1:1,-1:1),SDIM)
      REAL_T n1
      REAL_T n2
      REAL_T n_node1(SDIM)
      REAL_T n_node2(SDIM)
      REAL_T n_node1LS(SDIM)
      REAL_T n_node2LS(SDIM)
      REAL_T LSmain,LSopp
      REAL_T LS1_save(D_DECL(-1:1,-1:1,-1:1))
      REAL_T LS2_save(D_DECL(-1:1,-1:1,-1:1))
      INTEGER_T im_sort
      REAL_T dxmax,dxmaxLS
      REAL_T dxmin
      REAL_T delta_mgoni
      REAL_T volpos,facearea
      REAL_T cenpos(SDIM)
      REAL_T areacentroid(SDIM)
      REAL_T temperature_cen(nmat)
      REAL_T gradT(SDIM)
      REAL_T RR
      REAL_T dnrm(SDIM)
      REAL_T dxsten(SDIM)
      INTEGER_T im_primary_sten(D_DECL(-RD:RD,-RD:RD,-RD:RD))

      INTEGER_T crossing_status

      INTEGER_T cell_lo(3),cell_hi(3)

      REAL_T LS_CENTER(nmat)
      REAL_T LS_OPP(nmat)
      REAL_T LS_CENTER_EXTEND
      REAL_T LS_OPP_EXTEND

      REAL_T pm_val

      INTEGER_T ZEYU_imodel
      INTEGER_T ZEYU_ifgnbc
      REAL_T ZEYU_lambda
      REAL_T ZEYU_l_macro
      REAL_T ZEYU_l_micro
      REAL_T ZEYU_dgrid
      REAL_T ZEYU_d_closest
      REAL_T ZEYU_thet_s
      REAL_T ZEYU_thet_d_apparent
      REAL_T, intent(out) :: ZEYU_thet_d
      REAL_T, intent(out) :: ZEYU_u_cl
      REAL_T ZEYU_u_slip
      REAL_T ZEYU_mu_l
      REAL_T ZEYU_mu_g
      REAL_T ZEYU_sigma
      REAL_T angle_im
      REAL_T dist_to_CL
      INTEGER_T im_liquid,im_vapor

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)

      dxmin=dx(1)
      if (dx(2).lt.dxmin) then
       dxmin=dx(2)
      endif
      if (dx(SDIM).lt.dxmin) then
       dxmin=dx(SDIM)
      endif
      if (dxmin.gt.zero) then
       ! do nothing
      else
       print *,"dxmin invalid"
       stop
      endif

       ! -1 if use static angle
      call get_use_DCA(use_DCA)
      ZEYU_thet_d=zero
      totaludotn=zero

      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid initheightLS nten nten test", &
         nten,nten_test
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid initheightLS"
       stop
      endif
      if ((side.ne.1).and.(side.ne.-1)) then
       print *,"side invalid"
       stop
      endif
      if ((signside.ne.1).and.(signside.ne.-1)) then
       print *,"signside invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid85"
       stop
      endif
      if (RD.ne.4) then
       print *,"expecting RD=4=ngrow_distance in initheightLS RD=",RD
       stop
      endif
      if (RDx.ne.2*RD+1) then
       print *,"RDx invalid"
       stop
      endif
      if ((RD_HEIGHT.ne.RD).and.(RD_HEIGHT.ne.RD-1)) then
       print *,"RD_HEIGHT not RD or RD-1. RD_HEIGHT= ",RD_HEIGHT
       stop
      endif
      if ((conservative_tension_force.ne.0).and. &
          (conservative_tension_force.ne.1)) then
       print *,"conservative_tension_force invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (im.ge.im_opp) then
       print *,"im and im_opp invalid"
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid31"
       stop
      endif
      if ((im_opp.lt.1).or.(im_opp.gt.nmat)) then
       print *,"im_opp invalid init height mof im_opp=",im_opp
       stop
      endif
      if ((iten.lt.1).or.(iten.gt.nten)) then
       print *,"iten invalid"
       stop
      endif
      call get_iten(im,im_opp,iten_test,nmat)
      if (iten.ne.iten_test) then
       print *,"iten and iten_test differ"
       stop
      endif
      if (vol_sten.le.zero) then
       print *,"vol_sten invalid"
       stop
      endif

      if ((dircrit.lt.1).or.(dircrit.gt.SDIM)) then
       print *,"dircrit invalid initheightLS dircrit=",dircrit
       stop
      endif

      if (SDIM.eq.3) then
       klo_sten_short=-1
       khi_sten_short=1
       klo_sten_ht=-RD_HEIGHT
       khi_sten_ht=RD_HEIGHT
      else if (SDIM.eq.2) then
       klo_sten_short=0
       khi_sten_short=0
       klo_sten_ht=0
       khi_sten_ht=0
      else
       print *,"dimension bust"
       stop
      endif

      ii=0
      jj=0
      kk=0
      if (dircrit.eq.1) then
       ii=side
      else if (dircrit.eq.2) then
       jj=side
      else if ((dircrit.eq.3).and.(SDIM.eq.3)) then
       kk=side
      else
       print *,"dircrit invalid"
       stop
      endif

      do im_sort=1,nmat
       LS_CENTER(im_sort)=lssten(D_DECL(0,0,0),im_sort)
       LS_OPP(im_sort)=lssten(D_DECL(ii,jj,kk),im_sort)
      enddo

      if (LS_CENTER(im).ge.zero) then
       ! do nothing
      else if (LS_CENTER(im_opp).ge.zero) then
       ! do nothing
      else
       print *,"LS_CENTER is corrupt"
       do im_sort=1,nmat
        print *,"im_sort,LS_CENTER,LS_OPP ",im_sort,LS_CENTER(im_sort), &
           LS_OPP(im_sort)
       enddo
       print *,"level,finest_level ",level,finest_level
       print *,"bfact ",bfact
       print *,"dircrit= ",dircrit
       print *,"side,signside ",side,signside
       print *,"ii,jj,kk ",ii,jj,kk
       print *,"time= ",time
       print *,"im,im_opp,iten ",im,im_opp,iten
       print *,"RD,RDx,RD_HEIGHT ",RD,RDx,RD_HEIGHT
       stop
      endif 

      call get_LS_extend(LS_CENTER,nmat,iten,LS_CENTER_EXTEND)
      call get_LS_extend(LS_OPP,nmat,iten,LS_OPP_EXTEND)

      if (LS_CENTER_EXTEND*LS_OPP_EXTEND.gt.zero) then
       print *,"level set does not change sign"
       print *,"nmat,iten= ",nmat,iten
       print *,"LS_CENTER_EXTEND ",LS_CENTER_EXTEND
       print *,"LS_OPP_EXTEND ",LS_OPP_EXTEND
       print *,"icenter,jcenter,kcenter ",icenter,jcenter,kcenter
       print *,"level,finest_level ",level,finest_level
       print *,"LS_CENTER(im),LS_CENTER(im_opp) ", &
        LS_CENTER(im),LS_CENTER(im_opp)
       print *,"LS_OPP(im),LS_OPP(im_opp) ", &
        LS_OPP(im),LS_OPP(im_opp)
       stop
      endif

      if (levelrz.eq.0) then
       ! do nothing
      else if ((levelrz.eq.1).or. &
               (levelrz.eq.3)) then
       if (xcenter(1).le.zero) then
        print *,"xcenter(1) must be positive for RZ or RT"
        stop
       endif
      else
       print *,"initheightLS: levelrz invalid"
       stop
      endif

      do im_sort=1,nmat
       temperature_cen(im_sort)=mgoni_temp(D_DECL(0,0,0),im_sort)
      enddo

      call get_user_tension(xcenter,time, &
       fort_tension,user_tension, &
       temperature_cen, &
       nmat,nten,2)

      do dir2=1,SDIM
       mgoni_force(dir2)=zero
      enddo

      do i=-1,1
      do j=-1,1
      do k=klo_sten_short,khi_sten_short
       do imhold=1,nmat
        LSTEST(imhold)=lssten(D_DECL(i,j,k),imhold)
       enddo
       call get_LS_extend(LSTEST,nmat,iten,LS1_save(D_DECL(i,j,k)))
      enddo
      enddo
      enddo ! i,j,k (-1 ... 1)

       ! centroid in absolute coordinate system
       ! returns a volume fraction
      call getvolume( &
       bfact,dx,xsten,RDx, &
       LS1_save,volpos,facearea, &
       cenpos,areacentroid,VOFTOL,SDIM)

      if (facearea.ge.zero) then
       delta_mgoni=facearea/vol_sten
      else
       print *,"facearea invalid"
       stop
      endif

      call get_LSNRM_extend(LS_CENTER,nrmcenter,nmat,iten,nfluid)
      RR=one
      call prepare_normal(nfluid,RR,mag)
      if (mag.le.zero) then
       print *,"nfluid mag became corrupt"
       stop
      endif

       ! Marangoni force:
       ! (I-nn^T)(grad sigma) delta=
       ! (grad sigma - (grad sigma dot n)n ) delta
      if (fort_tension_slope(iten).ne.zero) then

       dotprod=zero

       do dir2=1,SDIM
        iofs=0
        jofs=0
        kofs=0
        RR=one
        if (dir2.eq.1) then
         iofs=1
         RR=one
        else if (dir2.eq.2) then 
         jofs=1
         if ((levelrz.eq.0).or. &
             (levelrz.eq.1)) then
          RR=one
         else if (levelrz.eq.3) then
          RR=xcenter(1)
          if (RR.le.zero) then
           print *,"RR invalid"
           stop
          endif
         else
          print *,"gradT: levelrz invalid"
          stop
         endif
        else if ((dir2.eq.3).and.(SDIM.eq.3)) then
         kofs=1
         RR=one
        else
         print *,"dir2 invalid"
         stop
        endif 

        gradT(dir2)=( &
         mgoni_temp(D_DECL(iofs,jofs,kofs),im)+ &
         mgoni_temp(D_DECL(iofs,jofs,kofs),im_opp)- &
         mgoni_temp(D_DECL(-iofs,-jofs,-kofs),im)- &
         mgoni_temp(D_DECL(-iofs,-jofs,-kofs),im_opp))/ &
         (two*RR*(xsten(2,dir2)-xsten(-2,dir2)))

        dotprod=dotprod+ &
         gradT(dir2)*fort_tension_slope(iten)*nfluid(dir2)
       enddo ! dir2

        ! tension=sigma_0 + slope*(T-T0)
       do dir2=1,SDIM
        mgoni_force(dir2)= &
         (fort_tension_slope(iten)*gradT(dir2)- &
          nfluid(dir2)*dotprod)*delta_mgoni
       enddo
      else if (fort_tension_slope(iten).eq.zero) then
       ! do nothing
      else
       print *,"fort_tension_slopes bust"
       stop
      endif 

      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       if (xcenter(1).le.zero) then
        print *,"xcenter invalid"
        stop
       endif
      else if (levelrz.eq.3) then
       if (xcenter(1).le.zero) then
        print *,"xcenter invalid"
        stop
       endif
      else 
       print *,"levelrz invalid init height ls 2"
       stop
      endif

      im3=0
      LSMAX=-1.0D+10
      do i=-RD_HEIGHT,RD_HEIGHT
      do j=-RD_HEIGHT,RD_HEIGHT
      do k=klo_sten_ht,khi_sten_ht

       do imhold=1,nmat
        LSTEST(imhold)=lssten(D_DECL(i,j,k),imhold)
       enddo
       call get_LS_extend(LSTEST,nmat,iten,lsdata(D_DECL(i,j,k)))
       call get_primary_material(LSTEST,nmat,imhold)
       im_primary_sten(D_DECL(i,j,k))=imhold

       if ((abs(i).le.1).and.(abs(j).le.1).and.(abs(k).le.1)) then
  
        if ((imhold.ge.1).and.(imhold.le.nmat)) then
         if ((imhold.ne.im).and. &
             (imhold.ne.im_opp)) then
          if (im3.eq.0) then
           im3=imhold
           LSMAX=LSTEST(imhold)
          else if ((im3.ge.1).and.(im3.le.nmat)) then
           if (LSTEST(imhold).gt.LSMAX) then
            im3=imhold
            LSMAX=LSTEST(imhold)
           endif
          else
           print *,"im3 invalid"
           stop
          endif
         else if ((imhold.eq.im).or.(imhold.eq.im_opp)) then
          ! do nothing
         else
          print *,"imhold invalid"
          stop
         endif
        else
         print *,"imhold invalid"
         stop
        endif
       endif ! abs(i)<=1 abs(j)<=1 abs(k)<=1
      enddo
      enddo
      enddo ! i,j,k (-RD_HEIGHT .. RD_HEIGHT)
     
      cos_angle=zero
      sin_angle=zero
      iten_13=0
      iten_23=0
      if ((im3.ge.1).and.(im3.le.nmat)) then
        ! sigma_{i,j}cos(theta_{i,k})=sigma_{j,k}-sigma_{i,k}
        ! theta_{ik}=0 => material i wets material k.
        ! im is material "i"  ("fluid" material)
        ! im_opp is material "j"
       call get_CL_iten(im,im_opp,im3,iten_13,iten_23, &
        user_tension,nten,cos_angle,sin_angle)

       if ((sin_angle.ge.zero).and.(cos_angle.ge.zero)) then
        angle_im=asin(sin_angle)
        ! sin_angle=sin(a)  cos_angle=cos(a)
        ! a=pi-asin(sin_angle)
        ! sin(pi-asin(sin_angle))=sin(pi)cos(-asin(sin_angle)+
        !  cos(pi)sin(-asin(sin_angle))=-sin(-asin(sin_angle))=sin_angle
       else if ((sin_angle.ge.zero).and.(cos_angle.le.zero)) then
        angle_im=Pi-asin(sin_angle)
       else
        print *,"sin_angle or cos_angle invalid"
        stop
       endif
      else if (im3.eq.0) then
       ! do nothing
      else
       print *,"im3 invalid"
       stop
      endif

      if (nmat.eq.1) then
       print *,"nmat==1 not supported"
       stop
      else if (nmat.eq.2) then
       if (im3.ne.0) then
        print *,"im3 invalid, nmat=",nmat
        stop
       endif
      else if (nmat.gt.2) then
       ! do nothing
      else
       print *,"nmat invalid"
       stop
      endif

      ! first: standard height function technique

! normal points towards "im"
! n=grad LS/|grad LS|

      imhold=im_primary_sten(D_DECL(0,0,0))
      do dir2=1,SDIM
       if (imhold.eq.im) then
        nfluid(dir2)=-nrmcenter(SDIM*(im_opp-1)+dir2)
       else if (imhold.eq.im_opp) then
        nfluid(dir2)=nrmcenter(SDIM*(im-1)+dir2)
       else
        print *,"either im or im_opp material must be at center"
        stop
       endif
      enddo ! dir2
      RR=one
      call prepare_normal(nfluid,RR,mag)
      if (mag.le.zero) then
       print *,"nfluid mag became corrupt"
       stop
      endif
 
       ! signside points towards im 
      if (nfluid(dircrit)*signside.le.zero) then
       print *,"nfluid or signside has wrong sign"
       stop
      endif

      n1d=signside

      itanlo=-1
      itanhi=1
      jtanlo=-1
      jtanhi=1

      if (dircrit.eq.1) then
       itan=2
       jtan=3
       if (SDIM.eq.2) then
        jtanlo=0 
        jtanhi=0 
       endif
      else if (dircrit.eq.2) then
       itan=1
       jtan=3
       if (SDIM.eq.2) then
        jtanlo=0 
        jtanhi=0 
       endif
      else if ((dircrit.eq.3).and.(SDIM.eq.3)) then
       itan=1
       jtan=2
      else
       print *,"dircrit invalid"
       stop
      endif

      do iwidth=itanlo,itanhi
      do jwidth=jtanlo,jtanhi

       iwidthnew=iwidth
        
       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        if (SDIM.ne.2) then
         print *,"dimension bust"
         stop
        endif
        if (itan.eq.1) then
         if (iwidth.eq.-1) then
          if (xsten(2*iwidth,1).le.zero) then
           iwidthnew=0
          endif
         endif
        endif
       else if (levelrz.eq.3) then
        if (xsten(-2,1).le.zero) then
         print *,"xsten cannot be negative for levelrz==3"
         stop
        endif
       else
        print *,"levelrz invalid initheight ls 4"
        stop
       endif
        
       do kheight=-RD_HEIGHT,RD_HEIGHT

        if (dircrit.eq.1) then
         icell=kheight
         jcell=iwidthnew
         kcell=jwidth
        else if (dircrit.eq.2) then
         icell=iwidthnew
         jcell=kheight
         kcell=jwidth
        else if ((dircrit.eq.3).and.(SDIM.eq.3)) then
         icell=iwidthnew
         jcell=jwidth
         kcell=kheight
        else
         print *,"dircrit invalid"
         stop
        endif

        columnLS(kheight)=lsdata(D_DECL(icell,jcell,kcell))
          
       enddo ! kheight

       call get_col_ht_LS( &
        crossing_status, &
        bfact,dx,xsten, &
        RD,RDx,RD_HEIGHT, &
        columnLS,col_ht, &
        dircrit,n1d, &
        nmat, &
        SDIM)

       if (crossing_status.eq.1) then
        ! do nothing
       else if (crossing_status.eq.0) then
        if (1.eq.0) then
         print *,"no crossing found iwidth,jwidth= ",iwidth,jwidth
        endif 
       else
        print *,"crossing_status invalid"
        stop
       endif

       htfunc(iwidth,jwidth)=col_ht
      enddo ! jwidth
      enddo ! iwidth

       ! dark material on the bottom
       ! phi=h(x,y)-z  z=h(x,y)   (dircrit=2)
       ! n=grad phi/|grad phi| = (h_x,h_y,-1)/sqrt(h_x^2+h_y^2+1)
       ! div n = (n_x)_x + (n_y)_y

      hxR=(htfunc(1,0)-htfunc(0,0))/(xsten(2,itan)-xsten(0,itan))
      hxL=(htfunc(0,0)-htfunc(-1,0))/(xsten(0,itan)-xsten(-2,itan))
      hx=(htfunc(1,0)-htfunc(-1,0))/(xsten(2,itan)-xsten(-2,itan))
      hxx=(hxR-hxL)/(xsten(1,itan)-xsten(-1,itan))
      hyR=zero
      hyL=zero
      hy=zero
      hxy=zero
      hyy=zero

      if (SDIM.eq.3) then
       hyR=(htfunc(0,1)-htfunc(0,0))/(xsten(2,jtan)-xsten(0,jtan))
       hyL=(htfunc(0,0)-htfunc(0,-1))/(xsten(0,jtan)-xsten(-2,jtan))
       hy=(htfunc(0,1)-htfunc(0,-1))/(xsten(2,jtan)-xsten(-2,jtan))
       hxy=(htfunc(1,1)-htfunc(-1,1)-htfunc(1,-1)+htfunc(-1,-1))/ &
           ( (xsten(2,itan)-xsten(-2,itan))* &
             (xsten(2,jtan)-xsten(-2,jtan)) )
       hyy=(hyR-hyL)/(xsten(1,jtan)-xsten(-1,jtan))
      endif
      arclen=one+hx**2+hy**2
      arclenx=two*hx*hxx+two*hy*hxy
      arcleny=two*hx*hxy+two*hy*hyy
      g=one/sqrt(arclen)
      gx=-half*(arclen**(-1.5))*arclenx
      gy=-half*(arclen**(-1.5))*arcleny

      if (levelrz.eq.0) then
        ! phi=h(x,y)-z  z=h(x,y)   (dircrit=2)
        !arclen=1+hx^2+hy^2
        !g(x,y)=arclen^{-1/2}  
        !n=grad phi/|grad phi|=(hx g,hy g,-g)
        !div n=hxx g + hx gx +hyy g + hy gy - gz
        !gx=(-1/2)arclen^{-3/2}arclenx
       curvHT=hxx*g+hx*gx+hyy*g+hy*gy
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       if (dircrit.eq.1) then
         ! phi=h(z)-r  
         ! arclen=1+hz^2 
         ! g(z)=arclen^(-1/2)
         ! n=(-g,hz g) 
         ! div n =-(r g)_r/r + (hz g)_z=-g/r+hzz g + gz hz
        RR=htfunc(0,0)
        if (RR.le.zero) then
         print *,"RR invalid"
         stop
        endif
        curvHT=-g/RR+hxx*g+gx*hx
       else if (dircrit.eq.2) then
         ! phi=h(r)-z
         ! arclen=1+hr^2
         ! g(r)=arclen^(-1/2)
         ! n=(hr g,-g)
         ! div n = (r hr g)_r/r = hr g/r+hrr g + hr gr
        RR=xcenter(1) 
        curvHT=hx*g/RR+hxx*g+hx*gx
       else
        print *,"dircrit invalid"
        stop
       endif
      else if (levelrz.eq.3) then
       if (dircrit.eq.1) then
         ! phi=h(y,z)-r
         ! grad phi=(-1,hy/r,hz)
         ! arclen=1+(hy/r)^2+hz^2
         ! g(r,y,z)=arclen^(-1/2)
         ! n=(-g,hy g/r,hz g)
         ! div n=(-g r)_r/r+(hy g/r)_y/r+(hz g)_z=
         ! -g_r-g/r+(hyy g + gy hy)/r^2+hzz g + hz gz
         ! g_r=-1/2 arclen^(-3/2) arclen_r
         ! arclen_r=hy^2 (-2/r^3)
         ! g_y=-1/2 arclen^(-3/2) arclen_y
         ! arclen_y=2hy hyy/r^2 +2hz hzy
         ! g_z=-1/2 arclen^(-3/2) arclen_z
         ! arclen_z=2hy hyz/r^2 +2hz hzz
        RR=htfunc(0,0)
        arclen=one+(hx/RR)**2+hy**2
        g=one/sqrt(arclen)
        arclenr=-two*(hx**2)/(RR**3)
        gr=-half*(arclen**(-1.5))*arclenr
        arclenx=two*hx*hxx/(RR**2)+two*hy*hxy
        arcleny=two*hx*hxy/(RR**2)+two*hy*hyy
        gx=-half*(arclen**(-1.5))*arclenx
        gy=-half*(arclen**(-1.5))*arcleny
        curvHT=-gr-g/RR+(hxx*g+hx*gx)/(RR**2)+hyy*g+hy*gy
       else if (dircrit.eq.2) then 
         ! phi=h(r,z)-y
         ! grad phi=(hr,-1/r,hz)
         ! arclen=hr^2+(1/r)^2+hz^2
         ! g(r,z)=arclen^(-1/2)
         ! n=(g hr,-g/r,hz g)
         ! div n=(g hr r)_r/r-(g/r)_y/r+(hz g)_z=
         !  g_r hr + g hrr + g hr/r +hzz g + hz gz=
         ! g_r=-1/2 arclen^(-3/2) arclen_r
         ! arclen_r=2 hr hrr -2/r^3 + 2hz hzr
         ! g_z=-1/2 arclen^(-3/2) arclen_z
         ! arclen_z=2 hr hrz +2hz hzz
        RR=xcenter(1) 
        arclen=hx**2+hy**2+(one/RR)**2
        g=one/sqrt(arclen)
        arclenx=two*hx*hxx-two/(RR**3)+two*hy*hxy
        arcleny=two*hx*hxy+two*hy*hyy
        gx=-half*(arclen**(-1.5))*arclenx
        gy=-half*(arclen**(-1.5))*arcleny
        curvHT=gx*hx+g*hxx+g*hx/RR+hyy*g+hy*gy
       else if ((dircrit.eq.3).and.(SDIM.eq.3)) then
         ! phi=h(r,y)-z
         ! grad phi=(hr,hy/r,-1)
         ! arclen=hr^2+(hy/r)^2+1
         ! g(r,y)=arclen^(-1/2)
         ! n=(g hr,g hy/r,-g)
         ! div n=(g hr r)_r/r+(g hy/r)_y/r=
         !  g_r h_r + g h_rr + g hr/r +(gy hy+g hyy)/r^2
         ! g_r=-1/2 arclen^(-3/2) arclen_r
         ! arclen_r=2 hr hrr + 2(hy/r)(hry r-hy)/r^2 =
         !          2 hr hrr + 2 hy hry/r^2 - 2 hy^2/r^3
         ! g_y=-1/2 arclen^(-3/2) arclen_y
         ! arclen_y=2 hr hry +2hy hyy/r^2
        RR=xcenter(1) 
        arclen=hx**2+(hy/RR)**2+one
        g=one/sqrt(arclen)
        arclenx=two*hx*hxx+two*hy*hxy/(RR**2)- &
                two*(hy**2)/(RR**3)
        arcleny=two*hx*hxy+two*hy*hyy/(RR**2)
        gx=-half*(arclen**(-1.5))*arclenx
        gy=-half*(arclen**(-1.5))*arcleny
        curvHT=gx*hx+g*hxx+g*hx/RR+(hyy*g+hy*gy)/(RR**2)
       else
        print *,"dircrit invalid"
        stop
       endif
      else
       print *,"levelrz invalid init height ls 5"
       stop
      endif


      if (n1d.eq.one) then
       curvHT=-curvHT
      else if (n1d.eq.-one) then
       ! do nothing
      else 
       print *,"n1d invalid"
       stop
      endif


      ! above: use height function 
      ! below: use finite difference 

      if ((im3.lt.0).or.(im3.gt.nmat).or. &
          (im3.eq.im).or.(im3.eq.im_opp)) then
       print *,"im3 invalid"
       stop
      endif

!    \
!     \
!   im \  im_opp
!  -----------
!    im3
! sigma_{im,im_opp}cos(theta_{im,im3})=sigma_{im_opp,im3}-sigma_{im,im3}
! theta_{im,im3}=0 => material im wets material im3 (angle between im and im3)
! iten_13 corresponds to "im,im3"
! iten_23 corresponds to "im_opp,im3"
! if nmat=4, 12 13 14 23 24 34

! nsolid points into the solid material (im3)
! nfluid points into im (outward from im_opp)
!  (nfluid derived from distance function)

! nfluid_cen is a normal to the im,im_opp interface
! and points towards the im material.

      do dir2=1,SDIM

       nfluid_def1(dir2)=nrmcenter(SDIM*(im-1)+dir2)
       nfluid_def2(dir2)=nrmcenter(SDIM*(im_opp-1)+dir2)
        
       if (imhold.eq.im) then
        nfluid_cen(dir2)=-nrmcenter(SDIM*(im_opp-1)+dir2)
       else if (imhold.eq.im_opp) then
        nfluid_cen(dir2)=nrmcenter(SDIM*(im-1)+dir2)
       else 
        print *,"either im or im_opp material must be at center"
        stop
       endif

       if (im3.eq.0) then
        nfluid_def3(dir2)=nfluid_cen(dir2)
       else if ((im3.ge.1).and.(im3.le.nmat)) then
        nfluid_def3(dir2)=nrmcenter(SDIM*(im3-1)+dir2)
       else
        print *,"im3 invalid"
        stop
       endif
     
      enddo ! dir2=1..sdim

      RR=one
      call prepare_normal(nfluid_cen,RR,mag)
      call prepare_normal(nfluid_def1,RR,mag1)
      call prepare_normal(nfluid_def2,RR,mag2)
      call prepare_normal(nfluid_def3,RR,mag3)

      if ((mag1.le.zero).or. &
          (mag2.le.zero).or. &
          (mag3.le.zero).or. &
          (mag.le.zero)) then
       print *,"err:nfluid_def1, nfluid_def2, nfluid_def3, or nfluid_cen"
       print *,"mag1,mag2,mag3,mag=",mag1,mag2,mag3,mag
       stop
      endif

      do i=-1,1
      do j=-1,1
      do k=klo_sten_short,khi_sten_short 

       do im_sort=1,nmat
        LSTEST(im_sort)=lssten(D_DECL(i,j,k),im_sort) 
       enddo

       ! nfluid is a normal to the im,im_opp interface
       ! and points towards the im material.

       do imhold=1,nmat
        do dir2=1,SDIM
         nrmtest(SDIM*(imhold-1)+dir2)= &
           nrmsten(D_DECL(i,j,k),SDIM*(imhold-1)+dir2)
        enddo
       enddo ! imhold

       call get_LSNRM_extend(LSTEST,nrmtest,nmat,iten,nfluid)

       RR=one
       call prepare_normal(nfluid,RR,mag)
       if (mag.eq.zero) then
        do dir2=1,SDIM
         nfluid(dir2)=nfluid_cen(dir2)
        enddo
       else if (mag.gt.zero) then
        ! do nothing
       else
        print *,"mag invalid"
        stop
       endif  

       do dir2=1,SDIM
        nfluid_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
       enddo

       if ((im3.ge.1).and.(im3.le.nmat)) then

        if (is_rigid(nmat,im3).eq.1) then

         if ((im3.eq.im).or.(im3.eq.im_opp)) then
          print *,"im3 invalid" 
          stop
         endif

         ! nsolid points into the solid
         do dir2=1,SDIM
          nsolid(dir2)=nrmsten(D_DECL(i,j,k),SDIM*(im3-1)+dir2)
         enddo
         RR=one
         call prepare_normal(nsolid,RR,mag)
         if (mag.le.zero) then
          do dir2=1,SDIM
           nsolid(dir2)=nfluid_def3(dir2)
          enddo
         endif
         do dir2=1,SDIM
          nsolid_save(D_DECL(i,j,k),dir2)=nsolid(dir2)
         enddo

        else if (is_rigid(nmat,im3).eq.0) then
         do dir2=1,SDIM
          nsolid_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
         enddo
        else 
         print *,"is_rigid invalid"
         stop
        endif

       else if (im3.eq.0) then

        do dir2=1,SDIM
         nsolid_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
        enddo

       else 
        print *,"im3 invalid"
        stop
       endif

       do dir2=1,SDIM
        nfluid(dir2)=nrmsten(D_DECL(i,j,k),SDIM*(im-1)+dir2)
       enddo
       RR=one
       call prepare_normal(nfluid,RR,mag)
       if (mag.gt.zero) then
        ! do nothing
       else if (mag.eq.zero) then
        do dir2=1,SDIM
         nfluid(dir2)=nfluid_def1(dir2)
        enddo
       else
        print *,"mag invalid"
        stop
       endif 
       do dir2=1,SDIM
        nmain_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
       enddo

       do dir2=1,SDIM
        nfluid(dir2)=nrmsten(D_DECL(i,j,k),SDIM*(im_opp-1)+dir2)
       enddo
       RR=one
       call prepare_normal(nfluid,RR,mag)
       if (mag.gt.zero) then
        ! do nothing
       else if (mag.eq.zero) then
        do dir2=1,SDIM
         nfluid(dir2)=nfluid_def2(dir2)
        enddo
       else
        print *,"mag invalid"
        stop
       endif 
       do dir2=1,SDIM
        nopp_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
       enddo

      enddo
      enddo
      enddo ! i,j,k=-1,..,1

      if ((im3.ge.1).and.(im3.le.nmat)) then

       if (is_rigid(nmat,im3).eq.1) then
        
        if ((im3.eq.im).or.(im3.eq.im_opp)) then
         print *,"im3 invalid" 
         stop
        endif
        if (user_tension(iten).eq.zero) then  
         ! do nothing
        else if (user_tension(iten).gt.zero) then

         ! implement dynamic contact angle algorithm here.
         ! first project nfluid onto the solid (im3) material

         if ((use_DCA.eq.-1).or.(use_DCA.ge.0)) then

          dotprod=zero
          do dir2=1,SDIM
           nfluid(dir2)=nfluid_save(D_DECL(0,0,0),dir2)
           nsolid(dir2)=nsolid_save(D_DECL(0,0,0),dir2) 
           dotprod=dotprod+nfluid(dir2)*nsolid(dir2) 
          enddo
           ! nproject=(I-nsolid nsolid^T)nfluid
          do dir2=1,SDIM
           nproject(dir2)=nfluid(dir2)-dotprod*nsolid(dir2)
          enddo
          RR=one
          call prepare_normal(nproject,RR,mag)

          if (mag.gt.zero) then
           ! find u dot nproject
           totaludotn=zero
           totalwt=zero

           do i2=-1,1
           do j2=-1,1
           do k2=klo_sten_short,khi_sten_short

            ! positive in the solid/ice
            LSTEST_EXTEND=lssten(D_DECL(i2,j2,k2),im3)
  
            if (LSTEST_EXTEND.lt.zero) then
             wt=one
            else 
             wt=CURVWT
            endif
            totalwt=totalwt+wt

            udotn=zero
            do dir2=1,SDIM
             udotn=udotn+velsten(D_DECL(i2,j2,k2),dir2)*nproject(dir2)
            enddo
            totaludotn=totaludotn+wt*udotn
           enddo
           enddo
           enddo  ! i2,j2,k2=-1...1

           if (totalwt.le.zero) then
            print *,"totalwt invalid"
            stop
           endif
           totaludotn=totaludotn/totalwt  ! contact line velocity
           liquid_viscosity= &
             visc_coef* &
             get_user_viscconst(1,fort_denconst(1),temperature_cen(1)) 

           if (1.eq.0) then
            do dir2=1,SDIM
             print *,"dir,x ",dir2,xcenter(dir2)
             print *,"dir,normal pointing into solid ", &
               dir2,nsolid(dir2)
             print *,"dir,CL normal pointing into liquid ", &
               dir2,nproject(dir2)
            enddo
            print *,"CL velocity ",totaludotn
            print *," cos static angle ",cos_angle
           endif 

            ! ZEYU_u_cl is positive if the contact line is advancing into
            ! the gas.
            ! nproject points towards the im material
           if (fort_denconst(im).ge. &
               fort_denconst(im_opp)) then
            im_liquid=im
            im_vapor=im_opp
            ZEYU_thet_s=angle_im  ! thet_s in the liquid.
            ZEYU_u_cl=-totaludotn
           else if (fort_denconst(im_opp).ge. &
                    fort_denconst(im)) then
            im_liquid=im_opp
            im_vapor=im
            ZEYU_thet_s=Pi-angle_im
            ZEYU_u_cl=totaludotn
           else
            print *,"fort_denconst bust"
            stop
           endif
           ZEYU_mu_l=fort_viscconst(im_liquid)
           ZEYU_mu_g=fort_viscconst(im_vapor)
           ZEYU_sigma=user_tension(iten)
           ZEYU_thet_d_apparent=ZEYU_thet_s
           ZEYU_thet_d=ZEYU_thet_d_apparent
           ZEYU_u_slip=zero ! unused for standard CL models.

           ZEYU_ifgnbc=0
           ZEYU_lambda=8.0D-7
           ZEYU_l_macro=dxmin
           ZEYU_l_micro=1.0D-9
           ZEYU_dgrid=dxmin
           dist_to_CL=zero
           ZEYU_d_closest=abs(dist_to_CL)

           if (use_DCA.eq.-1) then

            ZEYU_thet_d=acos(cos_angle)

            ! modify cos_angle (initialized above as static angle)
            ! use_DCA=-1 static angle
            ! use_DCA=0 static angle
            ! use_DCA=1 Jiang
            ! use_DCA=2 Kistler
           else if ((use_DCA.eq.0).or.(use_DCA.eq.1).or.(use_DCA.eq.2)) then

            call DCA_select_model(nproject,totaludotn,cos_angle, &
             liquid_viscosity,user_tension(iten),cos_angle,use_DCA)

            ZEYU_thet_d=acos(cos_angle)

           else if ((use_DCA.ge.101).and. & ! fort_ZEYU_DCA_SELECT>=1
                    (use_DCA.le.108)) then
            if (use_DCA.eq.101) then
             ! do nothing (GNBC model)

             ! cases 2,...,8 for Zeyu's code.
             ! case 2 Jiang 1970 ...
             ! case 8 model=Cox 1986
            else if ((use_DCA.ge.102).and.(use_DCA.le.108)) then
             ZEYU_imodel=use_DCA-100

             call dynamic_contact_angle(ZEYU_mu_l, ZEYU_mu_g, ZEYU_sigma, &
               ZEYU_thet_s, &
               ZEYU_imodel, ZEYU_ifgnbc, ZEYU_lambda, &
               ZEYU_l_macro, ZEYU_l_micro, &
               ZEYU_dgrid, ZEYU_d_closest, ZEYU_thet_d_apparent, &
               ZEYU_u_cl, ZEYU_u_slip, ZEYU_thet_d)
             if (im.eq.im_liquid) then
              cos_angle=cos(ZEYU_thet_d)
             else if (im.eq.im_vapor) then
              cos_angle=Pi-cos(ZEYU_thet_d)
             else
              print *,"im invalid"
              stop
             endif
                     
            else
             print *,"use_DCA bust"
             stop
            endif
           else
            print *,"use_DCA invalid"
            stop
           endif

           if (cos_angle.gt.one) then 
            cos_angle=one
           else if (cos_angle.lt.-one) then
            cos_angle=-one
           endif

           if (1.eq.0) then
            print *," cos dynamic angle ",cos_angle
           endif
          else if (mag.eq.zero) then
           ! do nothing (nproject has mag=0)
          else
           print *,"mag cannot be negative"
           stop
          endif 

         else
          print *,"use_DCA invalid"
          stop
         endif 
 
        else
         print *,"user_tension(iten) cannot be negative"
         stop
        endif

       else if (is_rigid(nmat,im3).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid"
        stop
       endif

      else if (im3.eq.0) then
       ! do nothing
      else
       print *,"im3 invalid"
       stop
      endif

      if (user_tension(iten).eq.zero) then
       gamma1=zero
       gamma2=zero
      else if (user_tension(iten).gt.zero) then

       if ((im3.ge.1).and.(im3.le.nmat)) then

        if ((im3.eq.im).or.(im3.eq.im_opp)) then
         print *,"im3 invalid" 
         stop
        endif

        if (is_rigid(nmat,im3).eq.1) then

         ! cos(theta_1)=(sigma_23-sigma_13)/sigma_12
         ! cos(theta_2)=(-sigma_23+sigma_13)/sigma_12
         if (use_DCA.eq.101) then ! GNBC
          gamma1=half*user_tension(iten)
          gamma2=half*user_tension(iten)
         else if (use_DCA.ge.-1) then
          gamma1=half*(one-cos_angle)
          gamma2=half*(one+cos_angle)
         else
          print *,"use_DCA invalid"
          stop
         endif 

        else if (is_rigid(nmat,im3).eq.0) then

         gamma1=half*(user_tension(iten)-user_tension(iten_23)+ &
           user_tension(iten_13))/user_tension(iten)
         gamma2=half*(user_tension(iten)+user_tension(iten_23)- &
           user_tension(iten_13))/user_tension(iten)

        else
         print *,"is_rigid invalid"
         stop
        endif

       else if (im3.eq.0) then

        gamma1=half
        gamma2=half

       else
        print *,"im3 invalid"
        stop
       endif
  
      else
       print *,"user_tension coeff invalid"
       stop
      endif

      do i=-1,1
      do j=-1,1
      do k=klo_sten_short,khi_sten_short

       do imhold=1,nmat
        LSTEST(imhold)=lssten(D_DECL(i,j,k),imhold)
       enddo

       ! FTEN=gamma1 K1 grad H1 + gamma2 K2 grad H2=
       !      (gamma1 K1 - gamma2 K2) grad H1=
       !      div (gamma1 n1 - gamma2 n2) grad H1 
       if ((gamma1.eq.zero).and. &
           (gamma2.eq.zero)) then

        LSmain=LSTEST(im)
        LSopp=LSTEST(im_opp)

        do dir2=1,SDIM
         nmain(dir2)=nmain_save(D_DECL(i,j,k),dir2)
         nopp(dir2)=nopp_save(D_DECL(i,j,k),dir2)
        enddo

       else if (im3.eq.0) then

        LSmain=LSTEST(im)
        LSopp=LSTEST(im_opp)

        do dir2=1,SDIM
         nmain(dir2)=nmain_save(D_DECL(i,j,k),dir2)
         nopp(dir2)=nopp_save(D_DECL(i,j,k),dir2)
        enddo

       else if (is_rigid(nmat,im3).eq.0) then

        LSmain=LSTEST(im)
        LSopp=LSTEST(im_opp)
        do dir2=1,SDIM
         nmain(dir2)=nmain_save(D_DECL(i,j,k),dir2)
         nopp(dir2)=nopp_save(D_DECL(i,j,k),dir2)
        enddo

       else if (is_rigid(nmat,im3).eq.1) then

        call get_LS_extend(LSTEST,nmat,iten,LSmain)
        LSopp=LSmain
        do dir2=1,SDIM
          ! nmain points from material im_opp into material im.
         nmain(dir2)=nfluid_save(D_DECL(i,j,k),dir2)
          ! nsolid_save points into the solid
          ! nopp points aways from the solid.
         nopp(dir2)=-nsolid_save(D_DECL(i,j,k),dir2)
        enddo
        if (use_DCA.eq.101) then
         do dir2=1,SDIM
          nghost(dir2)=nmain(dir2)
          nopp(dir2)=nmain(dir2)
         enddo
        else if (use_DCA.ge.-1) then
          ! nghost points from material im_opp into material im.
         call ghostnormal(nmain,nopp,cos_angle,nghost,nperp)
         do dir2=1,SDIM
          nopp(dir2)=nghost(dir2)
         enddo
        else
         print *,"use_DCA invalid"
         stop
        endif 

       else
        print *,"is_rigid(nmat,im3) invalid"
        stop
       endif
 
       do dir2=1,SDIM
        ncurv1_save(D_DECL(i,j,k),dir2)=nmain(dir2)
        ncurv2_save(D_DECL(i,j,k),dir2)=nopp(dir2)
       enddo  ! dir2
       LS1_save(D_DECL(i,j,k))=LSmain
       LS2_save(D_DECL(i,j,k))=LSopp

      enddo
      enddo
      enddo ! i,j,k (-1..1)

      do dir2=1,SDIM
       dnrm(dir2)=zero
      enddo

      do i=-1,1
      do j=-1,1
      do k=klo_sten_short,khi_sten_short

       imhold=im_primary_sten(D_DECL(i,j,k))

       do dir2=1,SDIM 
        n1=ncurv1_save(D_DECL(i,j,k),dir2)
        n2=ncurv2_save(D_DECL(i,j,k),dir2)
        if (im3.eq.0) then
         ngrid(dir2)=gamma1*n1-gamma2*n2
        else if (is_rigid(nmat,im3).eq.0) then
         ngrid(dir2)=gamma1*n1-gamma2*n2
        else if (is_rigid(nmat,im3).eq.1) then
         if (imhold.eq.im3) then
          ngrid(dir2)=n2
         else if ((imhold.ne.im3).and. &
                  (imhold.ge.1).and. &
                  (imhold.le.nmat)) then 
          ngrid(dir2)=n1
         else
          print *,"imhold invalid"
          stop
         endif
        else
         print *,"is_rigid(nmat,im3) invalid"
         stop
        endif
        ngrid_save(D_DECL(i,j,k),dir2)=ngrid(dir2)
       enddo ! dir2=1..sdim

      enddo
      enddo
      enddo ! i,j,k=-1 ... 1

      do dir2=1,SDIM
       xsten_curv(-2,dir2)=xsten(-2,dir2)
       xsten_curv(2,dir2)=xsten(2,dir2)
       xsten_curv(-1,dir2)=(xsten(0,dir2)+xsten(-2,dir2))/two
       xsten_curv(1,dir2)=(xsten(0,dir2)+xsten(2,dir2))/two
       xsten_curv(0,dir2)=(xsten_curv(1,dir2)+xsten_curv(-1,dir2))/two
      enddo ! dir2

      totalwt=zero

      do inode=-1,1,2
      do jnode=-1,1,2
      do knode=klo_sten_short,khi_sten_short,2

       node_index(1)=inode 
       node_index(2)=jnode 
       node_index(3)=knode 

       do dir2=1,SDIM 
        n_node1(dir2)=zero
        n_node2(dir2)=zero
        n_node1LS(dir2)=zero
        n_node2LS(dir2)=zero
       enddo

       im3_present_node=0

       cell_lo(3)=0
       cell_hi(3)=0

       do dir2=1,SDIM
        if (node_index(dir2).eq.-1) then
         cell_lo(dir2)=-1
         cell_hi(dir2)=0
        else if (node_index(dir2).eq.1) then
         cell_lo(dir2)=0
         cell_hi(dir2)=1
        else
         print *,"node_index invalid"
         stop
        endif
       enddo ! dir2=1..sdim

       wtnode=zero

       do i=cell_lo(1),cell_hi(1)
       do j=cell_lo(2),cell_hi(2)
       do k=cell_lo(3),cell_hi(3)

        wtnode=wtnode+one

        cell_index(1)=i 
        cell_index(2)=j
        cell_index(3)=k

        do dir2=1,SDIM

         if (cell_index(dir2).eq.cell_lo(dir2)) then
          pm_val=-one
         else if (cell_index(dir2).eq.cell_hi(dir2)) then
          pm_val=one
         else
          print *,"cell_lo or cell_hi invalid"
          stop
         endif

         n_node1LS(dir2)=n_node1LS(dir2)+pm_val*LS1_save(D_DECL(i,j,k))
         n_node2LS(dir2)=n_node2LS(dir2)+pm_val*LS2_save(D_DECL(i,j,k))

         n_node1(dir2)=n_node1(dir2)+ncurv1_save(D_DECL(i,j,k),dir2)
         n_node2(dir2)=n_node2(dir2)+ncurv2_save(D_DECL(i,j,k),dir2)

        enddo ! dir2

        imhold=im_primary_sten(D_DECL(i,j,k))
        if (imhold.eq.im3) then
         im3_present_node=1
        else if ((imhold.ne.im3).and. &
                 (imhold.ge.1).and. &
                 (imhold.le.nmat)) then 
         ! do nothing
        else
         print *,"imhold invalid"
         stop
        endif

       enddo
       enddo
       enddo ! i,j,k=cell_lo ... cell_hi

       if (wtnode.le.zero) then
        print *,"wtnode invalid"
        stop
       endif

       do dir2=1,SDIM

        dxsten(dir2)=xsten(2*cell_hi(dir2),dir2)- &
                     xsten(2*cell_lo(dir2),dir2)
        if (dxsten(dir2).le.zero) then
         print *,"dxsten invalid"
         stop
        endif
        RR=one
        if (dir2.eq.1) then
         ! do nothing
        else if (dir2.eq.2) then ! theta direction in cylindrical coord.
         if (levelrz.eq.3) then 
          RR=abs(xsten_curv(node_index(1),1))
         endif
        else if ((dir2.eq.3).and.(SDIM.eq.3)) then
         ! do nothing
        else
         print *,"dir2 invalid"
         stop
        endif
        gx=two/(RR*wtnode*dxsten(dir2))

        n_node1LS(dir2)=n_node1LS(dir2)*gx
        n_node2LS(dir2)=n_node2LS(dir2)*gx
        n_node1(dir2)=n_node1(dir2)/wtnode
        n_node2(dir2)=n_node2(dir2)/wtnode

       enddo ! dir2

       RR=one
       call prepare_normal(n_node1LS,RR,mag)
       if (mag.eq.zero) then
        do dir2=1,SDIM
         n_node1LS(dir2)=n_node1(dir2)
        enddo
       else if (mag.gt.zero) then
        ! do nothing
       else
        print *,"mag invalid"
        stop
       endif
       call prepare_normal(n_node2LS,RR,mag)
       if (mag.eq.zero) then
        do dir2=1,SDIM
         n_node2LS(dir2)=n_node2(dir2)
        enddo
       else if (mag.gt.zero) then
        ! do nothing
       else
        print *,"mag invalid"
        stop
       endif
       
       do dir2=1,SDIM
        n_node1(dir2)=n_node1LS(dir2)
       
        if (im3.eq.0) then
         n_node2(dir2)=n_node2LS(dir2)
        else if (is_rigid(nmat,im3).eq.0) then
         n_node2(dir2)=n_node2LS(dir2)
        else if (is_rigid(nmat,im3).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(nmat,im3) invalid"
         stop
        endif
       enddo ! dir2

       RR=one
       call prepare_normal(n_node1,RR,mag)
       call prepare_normal(n_node2,RR,mag)

       do dir2=1,SDIM
        n1=n_node1(dir2)
        n2=n_node2(dir2)

        if (im3.eq.0) then
         ngrid(dir2)=gamma1*n1-gamma2*n2
        else if (is_rigid(nmat,im3).eq.0) then
         ngrid(dir2)=gamma1*n1-gamma2*n2
        else if (is_rigid(nmat,im3).eq.1) then
         if (im3_present_node.eq.1) then
          ngrid(dir2)=n2
         else if (im3_present_node.eq.0) then 
          ngrid(dir2)=n1
         else
          print *,"im3_present_node invalid"
          stop
         endif
        else
         print *,"is_rigid(nmat,im3) invalid"
         stop
        endif
       enddo ! dir2=1..sdim

       do dir2=1,SDIM
        if (dir2.eq.1) then
         if (levelrz.eq.0) then
          RR=one
         else if ((levelrz.eq.1).or. &
                  (levelrz.eq.3)) then
          RR=abs(xsten_curv(node_index(1),1))
         else
          print *,"levelrz invalid initheightLS: RR"
          stop
         endif
        else if ((dir2.eq.2).or.(dir2.eq.SDIM)) then
         RR=one
        else
         print *,"dir2 invalid"
         stop
        endif

        dnrm(dir2)=dnrm(dir2)+node_index(dir2)*RR*ngrid(dir2)
       enddo ! dir2

       totalwt=totalwt+one

      enddo
      enddo
      enddo  ! inode,jnode,knode=-1,1,2

      if (totalwt.le.zero) then
       print *,"totalwt invalid in initheightLS"
       stop
      endif

       ! dxsten=(xsten(2)+xsten(0))/2-(xsten(-2)+xsten(0))/2=
       !   (xsten(2)-xsten(-2))/2
      do dir2=1,SDIM
        dxsten(dir2)=xsten_curv(1,dir2)-xsten_curv(-1,dir2)
        if (dxsten(dir2).le.zero) then
         print *,"dxsten invalid"
         stop
        endif 
      enddo ! dir2

      do dir2=1,SDIM

        if (dir2.eq.1) then
         if (levelrz.eq.0) then
          RR=one
         else if ((levelrz.eq.1).or. &
                  (levelrz.eq.3)) then
          RR=abs(xsten_curv(0,1))
         else
          print *,"levelrz invalid initheightLS: RR 3"
          stop
         endif
        else if (dir2.eq.2) then
         if (levelrz.eq.0) then
          RR=one
         else if (levelrz.eq.1) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          RR=one
         else if (levelrz.eq.3) then
          RR=abs(xsten_curv(0,1))
         else
          print *,"levelrz invalid initheightLS: RR 4"
          stop
         endif
        else if ((dir2.eq.3).and.(SDIM.eq.3)) then
         RR=one
        else
         print *,"dir2 invalid"
         stop
        endif

        dnrm(dir2)=two*dnrm(dir2)/(totalwt*RR*dxsten(dir2))

      enddo ! dir2=1..sdim

      curvFD=zero
      do dir2=1,SDIM
       curvFD=curvFD+dnrm(dir2)
      enddo

      if (conservative_tension_force.eq.1) then
       if (user_tension(iten).gt.zero) then
        do dir2=1,SDIM
         mgoni_force(dir2)=mgoni_force(dir2)- &
          user_tension(iten)*delta_mgoni*nfluid(dir2)*curvHT 
        enddo
       else if (user_tension(iten).eq.zero) then
        ! do nothing
       else
        print *,"user_tension(iten) invalid"
        stop
       endif
      else if (conservative_tension_force.eq.0) then
       ! do nothing
      else
       print *,"conservative_tension_force invalid"
       stop
      endif

      if (1.eq.0) then
       print *,"xcenter ",xcenter(1),xcenter(2),xcenter(SDIM)
       print *,"dircrit,side,signside ",dircrit,side,signside
       print *,"im3,curvFD,curvHT ",im3,curvFD,curvHT
      endif

      return
      end subroutine initheightLS


      end module levelset_module


      subroutine FORT_L1_DISTANCE( &
       level, &
       finest_level, &
       maskfab,DIMS(maskfab), &
       lsfab,DIMS(lsfab), &
       l1lsfab,DIMS(l1lsfab), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       ngrow_distance, &
       nmat)

      use global_utility_module
      use probcommon_module
      use MOF_routines_module
      use levelset_module

      IMPLICIT NONE

      INTEGER_T level,finest_level
      INTEGER_T nmat,ngrow_distance
      INTEGER_T DIMDEC(maskfab)
      INTEGER_T DIMDEC(lsfab)
      INTEGER_T DIMDEC(l1lsfab)

      REAL_T maskfab(DIMV(maskfab),4)
      REAL_T lsfab(DIMV(lsfab),nmat)
      REAL_T l1lsfab(DIMV(lsfab),nmat)

      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      REAL_T xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k,i1,j1,k1
      INTEGER_T im,im_primary,im_opp
      INTEGER_T klosten,khisten
      REAL_T LS_center(nmat)
      REAL_T LS_opp(nmat)
      INTEGER_T radmin(nmat)
      INTEGER_T radcurrent

      if (bfact.lt.1) then
       print *,"bfact invalid86"
       stop
      endif

      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid in levelstrip"
       stop
      endif

      if (ngrow_distance.ne.4) then
       print *,"ngrow_distance<>4 error in L1_DISTANCE"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskfab),ngrow_distance,-1,2870)
      call checkbound(fablo,fabhi,DIMS(lsfab),ngrow_distance,-1,2871)
      call checkbound(fablo,fabhi,DIMS(l1lsfab),ngrow_distance,-1,2872)
      
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      if (SDIM.eq.3) then
       klosten=-ngrow_distance
       khisten=ngrow_distance
      else if (SDIM.eq.2) then
       klosten=0
       khisten=0
      else
       print *,"dimension bust"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0)
 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       do im=1,nmat
        LS_center(im)=lsfab(D_DECL(i,j,k),im)
        radmin(im)=ngrow_distance+1
       enddo
       call get_primary_material(LS_center,nmat,im_primary)
       do i1=-ngrow_distance,ngrow_distance
       do j1=-ngrow_distance,ngrow_distance
       do k1=klosten,khisten

        radcurrent=abs(i1)
        if (radcurrent.lt.abs(j1)) then
         radcurrent=abs(j1)
        endif
        if (radcurrent.lt.abs(k1)) then
         radcurrent=abs(k1)
        endif

        if (radcurrent.eq.0) then

         ! do nothing

        else if ((radcurrent.ge.1).and. &
                 (radcurrent.le.ngrow_distance)) then 

         do im=1,nmat
          LS_opp(im)=lsfab(D_DECL(i+i1,j+j1,k+k1),im)
         enddo
         call get_primary_material(LS_opp,nmat,im_opp)
         do im=1,nmat
          if (((im.eq.im_primary).and. &
               (im.ne.im_opp)).or. &
              ((im.eq.im_opp).and. &
               (im.ne.im_primary))) then
           if (radmin(im).eq.ngrow_distance+1) then
            radmin(im)=radcurrent 
           else if (radmin(im).gt.radcurrent) then
            radmin(im)=radcurrent
           else if ((radmin(im).ge.1).and. &
                    (radmin(im).le.radcurrent)) then
            ! do nothing
           else
            print *,"radmin invalid"
            stop
           endif
          endif
         enddo ! im=1..nmat

        else
         print *,"radcurrent invalid"
         stop
        endif
       enddo
       enddo
       enddo ! i1,j1,k1

       do im=1,nmat
        if ((radmin(im).ge.1).and. &
            (radmin(im).le.ngrow_distance+1)) then
         if (im.eq.im_primary) then
          ! do nothing
         else if (im.ne.im_primary) then
          radmin(im)=-radmin(im)
         endif
        else
         print *,"radmin(im) invalid"
         stop
        endif
        l1lsfab(D_DECL(i,j,k),im)=radmin(im)
       enddo ! im=1..nmat

      enddo
      enddo
      enddo

      return
      end subroutine FORT_L1_DISTANCE

       ! facefab is initialized in FORT_FACEINIT
       ! centroids in facefab are in an absolute coordinate system.
      subroutine FORT_FACEPROCESS( &
       tid, &
       dir, &
       tessellate, &
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
      use global_distance_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module
      use levelset_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: dir
      INTEGER_T, intent(in) :: tessellate
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

      INTEGER_T i,j,k,iface,side
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
      if ((tessellate.ne.0).and.(tessellate.ne.1)) then
       print *,"tessellate invalid"
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

      call checkbound(fablo,fabhi,DIMS(dstfab),0,dir,2880)
      call checkbound(fablo,fabhi,DIMS(facefab),1,-1,2881)
      call checkbound(fablo,fabhi,DIMS(vofrecon),1,-1,2882)
      
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
        growlo,growhi,0,dir) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalfMAC,dir+1)

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
        endif
       else if (levelrz.eq.3) then
        if (dir.eq.0) then
         if (abs(xstenMAC(0,1)).le.VOFTOL*dx(1)) then
          at_RZ_face=1
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
           nmat,SDIM,2,SDIM+1,iface)
          leftface(im,dir2)=facefab(D_DECL(i-ii,j-jj,k-kk),iface)
          side=1
          call abs_array_index4(im,dir+1,side,dir2, &
           nmat,SDIM,2,SDIM+1,iface)
          rightface(im,dir2)=facefab(D_DECL(i,j,k),iface)
         enddo ! dir2
         frac_left(im)=leftface(im,1)
         frac_right(im)=rightface(im,1)

         if ((tessellate.eq.1).or.(is_rigid(nmat,im).eq.0)) then
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

        do im=1,nmat*ngeom_recon
         mofdata_left(im)=vofrecon(D_DECL(i-ii,j-jj,k-kk),im)
         mofdata_right(im)=vofrecon(D_DECL(i,j,k),im)
        enddo
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         do dir2=1,SDIM
          multi_cen_left(dir2,im)=mofdata_left(vofcomp+dir2)
          multi_cen_right(dir2,im)=mofdata_right(vofcomp+dir2)
         enddo
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

         if ((tessellate.eq.1).or.(is_rigid(nmat,im).eq.0)) then

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

          ! x_pair in absolute coordinate system.
        caller_id=12
        call multi_get_area_pairs( &
          bfact,dx, &
          xsten_right, &
          xsten_left, &
          nhalf, &
          mofdata_right, &
          mofdata_left, &
          nmat, &
          dir+1, &
          frac_pair, & ! left,right
          x_pair, & ! left,right
          SDIM, &
          geom_xtetlist(1,1,1,tid+1), &
          nmax, &
          geom_xtetlist_old(1,1,1,tid+1), &
          nmax, &
          nmax, &
          caller_id)

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

         if ((tessellate.eq.1).or.(is_rigid(nmat,ml).eq.0)) then

          do mr = 1, nmat

           if ((tessellate.eq.1).or.(is_rigid(nmat,mr).eq.0)) then

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

          enddo ! mr

         else if ((tessellate.eq.0).and.(is_rigid(nmat,ml).eq.1)) then
          ! do nothing
         else
          print *,"tessellate or is_rigid invalid"
          stop
         endif 

        enddo ! ml

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

       iface=0
       do ml = 1, nmat
        do mr = 1, nmat
          !(ml,mr,2) 
         iface=iface+1
         dstfab(D_DECL(i,j,k),iface)=frac_pair(ml,mr)
         iface=iface+1
         dstfab(D_DECL(i,j,k),iface)=dist_pair(ml,mr)
        enddo ! mr
       enddo ! ml
       if (iface.ne.nface_dst) then
        print *,"iface invalid"
        stop
       endif

      enddo
      enddo
      enddo  !i,j,k 

      return
      end subroutine FORT_FACEPROCESS


      subroutine FORT_DOTMASK_BUILD( &
       num_materials_face, &
       level, &
       finest_level, &
       dotmask,DIMS(dotmask), &
       maskfab,DIMS(maskfab), &
       vofrecon,DIMS(vofrecon), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       time, &
       nmat)

      use global_utility_module
      use probcommon_module
      use levelset_module

      IMPLICIT NONE

      INTEGER_T num_materials_face
      INTEGER_T level,finest_level
      INTEGER_T nmat
      INTEGER_T DIMDEC(dotmask)
      INTEGER_T DIMDEC(maskfab)
      INTEGER_T DIMDEC(vofrecon)

      REAL_T dotmask(DIMV(dotmask),num_materials_face)
      REAL_T maskfab(DIMV(maskfab))
      REAL_T vofrecon(DIMV(vofrecon),nmat*ngeom_recon)

      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      REAL_T xlo(SDIM),dx(SDIM)
      REAL_T time

      INTEGER_T i,j,k
      INTEGER_T im
      INTEGER_T imaskcov
      INTEGER_T vofcomp
      REAL_T voftest
 
      if (bfact.lt.1) then
       print *,"bfact invalid88"
       stop
      endif

      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid in dotmask_build"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(dotmask),0,-1,2890)
      call checkbound(fablo,fabhi,DIMS(maskfab),1,-1,2891)
      call checkbound(fablo,fabhi,DIMS(vofrecon),1,-1,2892)
      
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid dotmask_build"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid dotmask_build"
       print *,"ngeom_raw=",ngeom_raw
       stop
      endif

      if (num_materials_face.ne.nmat) then
       print *,"num_materials_face invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       ! imaskcov=tag if not covered by level+1 or outside the domain.
       imaskcov=NINT(maskfab(D_DECL(i,j,k)))

       if (imaskcov.eq.1) then
       
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         voftest=vofrecon(D_DECL(i,j,k),vofcomp)
         if (voftest.le.VOFTOL) then
          dotmask(D_DECL(i,j,k),im)=zero
         else
          dotmask(D_DECL(i,j,k),im)=one
         endif
        enddo ! im

       else if (imaskcov.eq.0) then

        do im=1,nmat
         dotmask(D_DECL(i,j,k),im)=zero
        enddo

       else
        print *,"imaskcov invalid"
        stop
       endif

      enddo ! k
      enddo ! j
      enddo ! i

      return
      end subroutine FORT_DOTMASK_BUILD

       ! for finding areas internal to a cell, perturb each internal 
       ! interface, find areas and volumes, then check for the difference 
       ! in volumes divided by eps times the area.
       ! 
      subroutine FORT_CELLFACEINIT( &
         tid, &
         tessellate, &
         nten, &
         level, &
         finest_level, &
         cface,DIMS(cface), &
         maskfab,DIMS(maskfab), &
         vofrecon,DIMS(vofrecon), &
         tilelo,tilehi, &
         fablo,fabhi,bfact, &
         rz_flag, &
         xlo,dx, &
         time,ngrow, &
         nmat,ncellfrac)

      use global_utility_module
      use global_distance_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module
      use levelset_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: tessellate
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: ncellfrac
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: DIMDEC(cface)
      INTEGER_T, intent(in) :: DIMDEC(maskfab)
      INTEGER_T, intent(in) :: DIMDEC(vofrecon)

      REAL_T, intent(out) :: cface(DIMV(cface),ncellfrac)
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
      INTEGER_T im
      INTEGER_T im_local
      INTEGER_T im_crit
      INTEGER_T im_opp
      INTEGER_T iface

      INTEGER_T ncellfrac_test
      INTEGER_T vofcomp
      INTEGER_T vofcomp2
      REAL_T vcenter(nmat)
      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T mofdatavalid(nmat*ngeom_recon)
      REAL_T local_facearea(nmat,nmat)
      REAL_T local_facearea_dimensional(nmat,nmat)
      REAL_T local_dist_to_line(nmat,nmat)
      REAL_T local_dist(nmat,nmat)
      REAL_T local_normal(nmat,nmat,SDIM)
      REAL_T local_facefrac(nmat)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      REAL_T dummy_tri(SDIM+1,SDIM)
      INTEGER_T nmax,ivert
      INTEGER_T dir2
      INTEGER_T shapeflag
      REAL_T multi_volume(nmat)
      REAL_T multi_volume_offset(nmat)
      REAL_T multi_cen(SDIM,nmat)
      REAL_T multi_cen_offset(SDIM,nmat)
      REAL_T multi_area(nmat)
      REAL_T total_facearea
      REAL_T total_facearea_mat(nmat)
      REAL_T uncaptured_volume_fraction
      REAL_T vfrac_fluid_sum
      REAL_T vfrac_solid_sum
      INTEGER_T loop_counter
      INTEGER_T num_processed_fluid
      INTEGER_T num_processed_solid
      INTEGER_T nmat_fluid
      INTEGER_T nmat_rigid
      INTEGER_T testflag
      REAL_T intercept
      REAL_T volcell
      REAL_T cencell(SDIM)
      REAL_T dist_tol,dxmax
      REAL_T dpair
      REAL_T areafrac
      INTEGER_T mask1,mask2
      INTEGER_T iten
      INTEGER_T nten_test
      INTEGER_T is_processed(nten)
      INTEGER_T caller_id
 
      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if ((tessellate.ne.0).and.(tessellate.ne.1)) then
       print *,"tessellate invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid CELLFACEINIT nten nten_test ",nten,nten_test
       stop
      endif

      nhalf=3
 
      nmax=POLYGON_LIST_MAX ! in: CELLFACEINIT

      do ivert=1,SDIM+1
      do dir2=1,SDIM
       dummy_tri(ivert,dir2)=zero
      enddo
      enddo
      if (bfact.lt.1) then
       print *,"bfact invalid89"
       stop
      endif

      ! (nmat,nmat,3+sdim)
      ! im_inside,im_outside,3+sdim --> area, dist_to_line, dist, line normal.
      ncellfrac_test=nmat*nmat*(3+SDIM)
      if (ncellfrac_test.ne.ncellfrac) then
       print *,"ncellfrac invalid cellfaceinit: ", &
        ncellfrac,ncellfrac_test
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in cellfaceinit"
       stop
      endif

      if (ngrow.lt.0) then
       print *,"ngrow<1 error in cellfaceinit"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(cface),ngrow,-1,2893)
      call checkbound(fablo,fabhi,DIMS(maskfab),ngrow,-1,2894)
      call checkbound(fablo,fabhi,DIMS(vofrecon),ngrow,-1,2895)
      
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid cellfaceinit"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid cellfaceinit"
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
       print *,"rz_flag invalid in cellfaceinit"
       stop
      endif

      call get_dxmax(dx,bfact,dxmax)
      dist_tol=FACETOL_DVOL*dxmax

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,ngrow) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       mask1=NINT(maskfab(D_DECL(i,j,k),1))
       mask2=NINT(maskfab(D_DECL(i,j,k),2))

       if ((mask2.eq.1).or.(mask1.eq.0)) then

        ! xsten(0,dir)  dir=1..sdim  is center of cell
        ! xsten(1,dir)  is right coordinate in dir direction
        ! xsten(-1,dir) is left coordinate in dir direction.
        call CISBOX(xsten,nhalf, &
         xlo,dx,i,j,k, &
         bfact,level, &
         volcell,cencell,SDIM)   
        
        do im=1,nmat*ngeom_recon
         mofdata(im)=vofrecon(D_DECL(i,j,k),im)
        enddo

         ! vcenter = volume fraction 
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         vcenter(im)=mofdata(vofcomp)
        enddo ! im

        call check_full_cell_vfrac(vcenter,tessellate,nmat,im_crit)

        iface=0
        do im=1,nmat
         do im_opp=1,nmat
          iface=iface+1
          local_facearea(im,im_opp)=zero
          local_facearea_dimensional(im,im_opp)=zero
          local_dist_to_line(im,im_opp)=zero
          local_dist(im,im_opp)=zero
          do dir2=1,SDIM
           local_normal(im,im_opp,dir2)=zero
          enddo
         enddo ! im_opp
        enddo ! im

        if (iface.ne.nmat*nmat) then
         print *,"iface invalid"
         stop
        endif

        if ((im_crit.ge.1).and.(im_crit.le.nmat)) then
         ! do nothing, there are no internal faces
        else if (im_crit.eq.0) then

         ! sum F_fluid=1  sum F_solid<=1
         call make_vfrac_sum_ok_copy(tessellate,mofdata,mofdatavalid, &
           nmat,SDIM,3)

         shapeflag=0

          ! base case
          ! in: FORT_CELLFACEINIT
         call multi_get_volume_grid( &
          tessellate, &
          bfact,dx, &
          xsten,nhalf, &
          mofdatavalid, &
          xsten,nhalf, &
          dummy_tri, &
          multi_volume, &
          multi_cen, &
          multi_area, &
          geom_xtetlist_uncapt(1,1,1,tid+1), &
          nmax, &
          nmax, &
          nmat,SDIM, &
          shapeflag,3)

         do iten=1,nten
          is_processed(iten)=0
         enddo

         do im=1,nmat
          total_facearea_mat(im)=zero
         enddo

         uncaptured_volume_fraction=one

         nmat_rigid=0
         nmat_fluid=0
         vfrac_fluid_sum=zero
         vfrac_solid_sum=zero

         do im=1,nmat
          if (is_rigid(nmat,im).eq.0) then
           nmat_fluid=nmat_fluid+1
           vfrac_fluid_sum=vfrac_fluid_sum+vcenter(im)
          else if (is_rigid(nmat,im).eq.1) then
           nmat_rigid=nmat_rigid+1
           vfrac_solid_sum=vfrac_solid_sum+vcenter(im)
          else
           print *,"is_rigid invalid"
           stop
          endif
         enddo ! im=1..nmat

         if (abs(one-vfrac_fluid_sum).gt.VOFTOL) then
          print *,"vfrac_fluid_sum invalid"
          stop
         endif
         if ((vfrac_solid_sum.gt.one+VOFTOL).or. &
             (vfrac_solid_sum.lt.zero)) then
          print *,"vfrac_solid_sum invalid"
          stop
         endif

         if (nmat_fluid+nmat_rigid.ne.nmat) then
          print *,"nmat_fluid or nmat_rigid invalid"
          stop
         endif
         num_processed_solid=0
         num_processed_fluid=0

         if ((tessellate.eq.1).and.(vfrac_solid_sum.gt.zero)) then

          loop_counter=0
          do while ((loop_counter.lt.nmat_rigid).and. &
                    (num_processed_solid.lt.nmat_rigid).and. &
                    (uncaptured_volume_fraction.gt. &
                     one-vfrac_solid_sum))

            ! F,CEN,ORDER,SLOPE,INTERCEPT
           do im=1,nmat
            if (is_rigid(nmat,im).eq.1) then
             vofcomp=(im-1)*ngeom_recon+1
             testflag=NINT(mofdatavalid(vofcomp+SDIM+1))

             if (testflag.eq.1) then

              num_processed_solid=num_processed_solid+1

              intercept=mofdatavalid(vofcomp+2*SDIM+2)
  
              uncaptured_volume_fraction= &
                uncaptured_volume_fraction-mofdatavalid(vofcomp)
              if (uncaptured_volume_fraction.lt.FACETOL_DVOL) then
               uncaptured_volume_fraction=zero
              endif
              if (uncaptured_volume_fraction.gt.zero) then ! valid interface.

               ! dist=intercept+slopes dot (x-x0)
               ! perturb interface into the other materials
               mofdatavalid(vofcomp+2*SDIM+2)=intercept+half*FACETOL_DVOL*dx(1)
       
               caller_id=3
 
               if (1.eq.0) then
                if ((SDIM.eq.3).and. &
                    (i.eq.12).and.(j.eq.56).and.(k.eq.9)) then
                 caller_id=101
                endif
               endif
 
                ! solid case
                ! in: FORT_CELLFACEINIT
                ! no need to compute multi_area here.
                ! also, target volume is a cube, not a tet.
               call multi_get_volume_grid_simple( &
                tessellate, &
                bfact,dx,xsten,nhalf, &
                mofdatavalid, &
                xsten,nhalf, &
                multi_volume_offset, &
                multi_cen_offset, &
                geom_xtetlist_uncapt(1,1,1,tid+1), &
                nmax, &
                nmax, &
                nmat, &
                SDIM,caller_id)

               mofdatavalid(vofcomp+2*SDIM+2)=intercept

               if (multi_volume_offset(im).gt.multi_volume(im)) then

                if (multi_area(im).gt.zero) then

                 do im_opp=1,nmat
                  local_facefrac(im_opp)=zero
                 enddo
                 total_facearea=zero

                 do im_opp=1,nmat
                  if (im_opp.ne.im) then
                   call get_iten(im,im_opp,iten,nmat)
                   if (is_processed(iten).eq.0) then 
                    local_facefrac(im_opp)= &
                     abs(multi_volume(im_opp)-multi_volume_offset(im_opp))
                    total_facearea=total_facearea+local_facefrac(im_opp)
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

                 if (total_facearea.gt.zero) then
                  do im_opp=1,nmat
                   if (im_opp.ne.im) then
                    call get_iten(im,im_opp,iten,nmat)
                    if (is_processed(iten).eq.0) then 
                     is_processed(iten)=1
                     local_facefrac(im_opp)= &
                       local_facefrac(im_opp)/total_facearea
                     local_facearea(im,im_opp)=local_facefrac(im_opp)

                     if (local_facefrac(im_opp).gt.zero) then
                      ! vfrac,cen,order,slope,intercept
                      vofcomp=(im-1)*ngeom_recon+1
                      do dir2=1,SDIM
                       local_normal(im,im_opp,dir2)= &
                        -mofdata(vofcomp+SDIM+1+dir2)
                       local_normal(im_opp,im,dir2)= &
                        mofdata(vofcomp+SDIM+1+dir2)
                      enddo
                      call dist_centroid_line( &
                       xsten,nhalf, &
                       im,im, & ! distance from line(im) to point(im)
                       mofdata, &
                       nmat, &
                       cencell, &  ! cell centroid. 
                       dist_tol, &
                       SDIM, &
                       local_dist_to_line(im,im_opp))
 
                      call dist_centroid_line( &
                       xsten,nhalf, &
                       im,im_opp, & ! distance from line(im) to point(im_opp)
                       mofdata, &
                       nmat, &
                       cencell, &  ! cell centroid. 
                       dist_tol, &
                       SDIM, &
                       local_dist_to_line(im_opp,im))
 
                      dpair=zero
                      vofcomp=(im-1)*ngeom_recon+1
                      vofcomp2=(im_opp-1)*ngeom_recon+1
                      do dir2=1,SDIM
                       dpair=dpair+ &
                        (mofdatavalid(vofcomp+dir2)- &
                         mofdatavalid(vofcomp2+dir2))**2
                      enddo
                      dpair=sqrt(dpair)
                      local_dist(im,im_opp)=dpair
                      local_dist(im_opp,im)=dpair
                     endif ! local_facefrac(im_opp)>0.0
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
                  enddo ! im_opp=1..nmat

                  total_facearea_mat(im)=multi_area(im)
                 else
                  print *,"opposite material bad:total_facearea=",total_facearea
                  stop
                 endif
                else
                 print *,"im boundary disappeared 3"
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
                print *,"im region should grow 3"
                print *,"im=",im
                print *,"multi_volume(im)=",multi_volume(im)
                print *,"multi_volume_offset(im)=",multi_volume_offset(im)
                print *,"multi_area(im)=",multi_area(im)
                print *,"intercept=",intercept
                print *,"mofdatavalid(vofcomp)=",mofdatavalid(vofcomp)
                print *,"tessellate=",tessellate
                print *,"nmat_rigid=",nmat_rigid
                print *,"nmat_fluid=",nmat_fluid
                print *,"vfrac_fluid_sum ",vfrac_fluid_sum
                print *,"vfrac_solid_sum ",vfrac_solid_sum
                print *,"loop_counter ",loop_counter
                print *,"num_processed_solid ",num_processed_solid
                print *,"num_processed_fluid ",num_processed_fluid
                print *,"uncaptured_volume_fraction ", &
                   uncaptured_volume_fraction
                print *,"FACETOL_DVOL ",FACETOL_DVOL
                print *,"bfact,level ",bfact,level
                print *,"dx(1) ",dx(1)
                print *,"i,j,k,xlo,volcell,xsten(cen) ", &
                   i,j,k,xlo(1),xlo(2),xlo(SDIM),volcell, &
                   xsten(0,1),xsten(1,1),xsten(2,1)
                do im_local=1,nmat
                 print *,"in loop: im_local=",im_local
                 vofcomp=(im_local-1)*ngeom_recon+1
                 print *,"mofdatavalid(vofcomp)=",mofdatavalid(vofcomp)
                 print *,"mofdatavalid(vofcomp+SDIM+1)=", &
                    mofdatavalid(vofcomp+SDIM+1)
                 print *,"mofdatavalid(vofcomp+2*SDIM+2)=", &
                    mofdatavalid(vofcomp+2*SDIM+2)
                 do dir2=1,SDIM
                  print *,"mofdatavalid(vofcomp+SDIM+1+dir2)=", &
                    mofdatavalid(vofcomp+SDIM+1+dir2)
                 enddo
                enddo  ! im=1..nmat
                stop
               endif
  
              endif ! uncaptured_volume_fraction>0

             else if (testflag.eq.0) then
              ! do nothing
             else
              print *,"testflag invalid"
              stop
             endif 

            else if (is_rigid(nmat,im).eq.0) then
             ! do nothing
            else
             print *,"is_rigid invalid"
             stop
            endif

           enddo ! im=1..nmat
           loop_counter=loop_counter+1
          enddo  ! while 
                 ! loop_counter<nmat_fluid and
                 ! num_processed_fluid<nmat_fluid and
                 ! uncaptured_volume_fraction>0

         else if ((tessellate.eq.0).or.(vfrac_solid_sum.eq.zero)) then
          ! do nothing
         else
          print *,"tessellate or vfrac_solid_sum invalid"
          stop
         endif

         loop_counter=0
         do while ((loop_counter.lt.nmat_fluid).and. &
                   (num_processed_fluid.lt.nmat_fluid).and. &
                   (uncaptured_volume_fraction.gt.zero))

           ! F,CEN,ORDER,SLOPE,INTERCEPT
          do im=1,nmat
           if (is_rigid(nmat,im).eq.0) then
            vofcomp=(im-1)*ngeom_recon+1
            testflag=NINT(mofdatavalid(vofcomp+SDIM+1))

            if (testflag.eq.num_processed_fluid+1) then

             num_processed_fluid=num_processed_fluid+1

             intercept=mofdatavalid(vofcomp+2*SDIM+2)
 
             uncaptured_volume_fraction= &
               uncaptured_volume_fraction-mofdatavalid(vofcomp)
             if (uncaptured_volume_fraction.lt.FACETOL_DVOL) then
              uncaptured_volume_fraction=zero
             endif
             if (uncaptured_volume_fraction.gt.zero) then ! valid interface.

              ! dist=intercept+slopes dot (x-x0)
              ! perturb interface into the other materials
              mofdatavalid(vofcomp+2*SDIM+2)=intercept+half*FACETOL_DVOL*dx(1)

               ! fluid case
               ! in: FORT_CELLFACEINIT
              call multi_get_volume_grid_simple( &
               tessellate, &
               bfact,dx,xsten,nhalf, &
               mofdatavalid, &
               xsten,nhalf, &
               multi_volume_offset, &
               multi_cen_offset, &
               geom_xtetlist_uncapt(1,1,1,tid+1), &
               nmax, &
               nmax, &
               nmat, &
               SDIM,3)

              mofdatavalid(vofcomp+2*SDIM+2)=intercept

              if (multi_volume(im).gt.zero) then

               if (multi_volume_offset(im).gt.multi_volume(im)) then

                if (multi_area(im).gt.zero) then

                 do im_opp=1,nmat
                  local_facefrac(im_opp)=zero
                 enddo
                 total_facearea=zero

                 do im_opp=1,nmat

                  if ((is_rigid(nmat,im_opp).eq.0).or. &
                      (tessellate.eq.1)) then

                   if (im_opp.ne.im) then
                    call get_iten(im,im_opp,iten,nmat)
                    if (is_processed(iten).eq.0) then 
                     local_facefrac(im_opp)= &
                      abs(multi_volume(im_opp)-multi_volume_offset(im_opp))
                     total_facearea=total_facearea+local_facefrac(im_opp)
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
                  else if ((is_rigid(nmat,im_opp).eq.1).and. &
                           (tessellate.eq.0)) then
                   ! do nothing
                  else
                   print *,"is_rigid or tessellate invalid"
                   stop
                  endif

                 enddo !im_opp=1..nmat

                 if (total_facearea.gt.zero) then
                  do im_opp=1,nmat
                   if ((is_rigid(nmat,im_opp).eq.0).or. &
                       (tessellate.eq.1)) then

                    if (im_opp.ne.im) then
                     call get_iten(im,im_opp,iten,nmat)
                     if (is_processed(iten).eq.0) then 
                      is_processed(iten)=1
                      local_facefrac(im_opp)= &
                        local_facefrac(im_opp)/total_facearea
                      local_facearea(im,im_opp)=local_facefrac(im_opp)

                      if (local_facefrac(im_opp).gt.zero) then
                       ! vfrac,cen,order,slope,intercept
                       vofcomp=(im-1)*ngeom_recon+1
                       do dir2=1,SDIM
                        local_normal(im,im_opp,dir2)= &
                         -mofdata(vofcomp+SDIM+1+dir2)
                        local_normal(im_opp,im,dir2)= &
                         mofdata(vofcomp+SDIM+1+dir2)
                       enddo
                       call dist_centroid_line( &
                        xsten,nhalf, &
                        im,im, & ! distance from line(im) to point(im)
                        mofdata, &
                        nmat, &
                        cencell, &  ! cell centroid. 
                        dist_tol, &
                        SDIM, &
                        local_dist_to_line(im,im_opp))
  
                       call dist_centroid_line( &
                        xsten,nhalf, &
                        im,im_opp, & ! distance from line(im) to point(im_opp)
                        mofdata, &
                        nmat, &
                        cencell, &  ! cell centroid. 
                        dist_tol, &
                        SDIM, &
                        local_dist_to_line(im_opp,im))
  
                       dpair=zero
                       vofcomp=(im-1)*ngeom_recon+1
                       vofcomp2=(im_opp-1)*ngeom_recon+1
                       do dir2=1,SDIM
                        dpair=dpair+ &
                         (mofdatavalid(vofcomp+dir2)- &
                          mofdatavalid(vofcomp2+dir2))**2
                       enddo
                       dpair=sqrt(dpair)
                       local_dist(im,im_opp)=dpair
                       local_dist(im_opp,im)=dpair
                      endif ! local_facefrac(im_opp)>0.0
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
                   else if ((is_rigid(nmat,im_opp).eq.1).and. &
                            (tessellate.eq.0)) then
                    ! do nothing
                   else
                    print *,"is_rigid or tessellate invalid"
                    stop
                   endif
                  enddo ! im_opp

                  total_facearea_mat(im)=multi_area(im)
                 else
                  print *,"opposite material bad:total_facearea=",total_facearea
                  stop
                 endif
                else
                 print *,"im boundary disappeared 4"
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
                 do vofcomp=1,nmat*ngeom_recon
                  print *,"vofcomp,mofdatavalid(vofcomp) ",vofcomp, &
                   mofdatavalid(vofcomp)
                 enddo
                 stop
                endif

               else if (multi_volume(im).eq.zero) then
                ! do nothing, material is completely hidden by solid(s)
               else
                print *,"im region should grow 4"
                print *,"im= ",im
                print *,"tessellate=",tessellate
                print *,"uncaptured_volume_fraction=", &
                 uncaptured_volume_fraction
                print *,"multi_volume_offset(im)=", &
                  multi_volume_offset(im)
                print *,"multi_volume(im)=", &
                  multi_volume(im)
                print *,"num_processed_fluid=",num_processed_fluid
                print *,"loop_counter=",loop_counter
                do vofcomp=1,nmat*ngeom_recon
                 print *,"vofcomp,mofdatavalid(vofcomp) ",vofcomp, &
                   mofdatavalid(vofcomp)
                enddo
                stop
               endif

              else if (multi_volume(im).eq.zero) then
               if ((vcenter(im).gt.FACETOL_DVOL).and. &
                   (vfrac_solid_sum.eq.zero)) then
                print *,"multi_volume(im) should not be zero"
                stop
               endif
              else
               print *,"multi_volume(im) invalid"
               stop
              endif
 
             endif ! uncaptured_volume_fraction>0

            endif  ! testflag==num_processed_fluid+1

           else if (is_rigid(nmat,im).eq.1) then
            ! do nothing
           else
            print *,"is_rigid invalid"
            stop
           endif

          enddo ! im=1..nmat
          loop_counter=loop_counter+1
         enddo  ! while 
                ! loop_counter<nmat_fluid and
                ! num_processed_fluid<nmat_fluid and
                ! uncaptured_volume_fraction>0

        else
         print *,"im_crit out of range"
         stop
        endif

        do im=1,nmat
         do im_opp=1,nmat
          if (im.ne.im_opp) then
           areafrac=local_facearea(im,im_opp)
           if (areafrac.lt.zero) then
            print *,"areafrac invalid1 areafrac=",areafrac
            print *,"im,im_opp ",im,im_opp 
            stop
           else if (areafrac.eq.zero) then
            ! do nothing
           else if ((areafrac.gt.zero).and. &
                    (areafrac.le.one)) then
            if (total_facearea_mat(im).le.zero) then
             print *,"total_facearea_mat(im) invalid"
             stop
            endif
            local_facearea_dimensional(im,im_opp)= &
                areafrac*total_facearea_mat(im)
            local_facearea_dimensional(im_opp,im)= &
                local_facearea_dimensional(im,im_opp)
           else
            print *,"areafrac invalid2 areafrac=",areafrac
            print *,"im,im_opp ",im,im_opp 
            print *,"is_processed im,im_opp ",is_processed(im), &
              is_processed(im_opp)
            print *,"total_facearea_mat(im) ",total_facearea_mat(im)  
            print *,"total_facearea_mat(im_opp) ",total_facearea_mat(im_opp) 
            print *,"multi_area(im) ",multi_area(im) 
            print *,"multi_area(im_opp) ",multi_area(im_opp) 
            print *,"uncaptured_volume_fraction ",uncaptured_volume_fraction
            stop
           endif
          else if (im.eq.im_opp) then
           ! do nothing
          else
           print *,"im or im_opp invalid"
           stop
          endif
         enddo ! im_opp
        enddo ! im

        iface=0
        do im=1,nmat
         do im_opp=1,nmat
          iface=iface+1
          cface(D_DECL(i,j,k),iface)=local_facearea_dimensional(im,im_opp)
          iface=iface+1
          cface(D_DECL(i,j,k),iface)=local_dist_to_line(im,im_opp)
          iface=iface+1
          cface(D_DECL(i,j,k),iface)=local_dist(im,im_opp)
          do dir2=1,SDIM
           iface=iface+1
           cface(D_DECL(i,j,k),iface)=local_normal(im,im_opp,dir2)
          enddo
         enddo ! im_opp
        enddo ! im
           
        if (iface.ne.ncellfrac) then
         print *,"ncellfrac invalid"
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
      end subroutine FORT_CELLFACEINIT


      subroutine FORT_CURVSTRIP( &
       post_restart_flag, &
       conservative_tension_force, &
       level, &
       finest_level, &
       curv_min, &
       curv_max, &
       nhistory, &
       history_dat, &
       DIMS(history_dat), &
       maskcov,DIMS(maskcov), &
       vol,DIMS(vol), &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz), &
       masknbr,DIMS(masknbr), &
       LSPC,DIMS(LSPC), &
       LSHO,DIMS(LSHO), &
       curvfab,DIMS(curvfab), &
       velfab,DIMS(velfab), &
       denfab,DIMS(denfab), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       bfact_grid, &
       rz_flag, &
       xlo,dx, &
       time, &
       visc_coef, &
       nmat,nten, & 
       num_curv, &
       ngrow_distance, &
       RD)

      use global_utility_module
      use global_distance_module
      use probcommon_module
      use MOF_routines_module
      use levelset_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: nhistory
      INTEGER_T, intent(in) :: post_restart_flag
      INTEGER_T, intent(in) :: conservative_tension_force
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: ngrow_distance
      INTEGER_T, intent(in) :: RD
      INTEGER_T, intent(in) :: num_curv
      INTEGER_T icurv
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: visc_coef

      INTEGER_T, intent(in) :: DIMDEC(history_dat)
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(masknbr)
      INTEGER_T, intent(in) :: DIMDEC(LSPC)
      INTEGER_T, intent(in) :: DIMDEC(LSHO)
      INTEGER_T, intent(in) :: DIMDEC(curvfab)
      INTEGER_T, intent(in) :: DIMDEC(velfab)
      INTEGER_T, intent(in) :: DIMDEC(denfab)

      INTEGER_T, intent(in) :: DIMDEC(vol)
      INTEGER_T, intent(in) :: DIMDEC(areax)
      INTEGER_T, intent(in) :: DIMDEC(areay)
      INTEGER_T, intent(in) :: DIMDEC(areaz)

      REAL_T, intent(out) :: curv_min
      REAL_T, intent(out) :: curv_max

      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: vol(DIMV(vol))
      REAL_T, intent(in) :: areax(DIMV(areax))
      REAL_T, intent(in) :: areay(DIMV(areay))
      REAL_T, intent(in) :: areaz(DIMV(areaz))

      REAL_T, intent(out) :: history_dat(DIMV(history_dat),nhistory)
      REAL_T, intent(in) :: masknbr(DIMV(masknbr),4)
      REAL_T, intent(in) :: LSPC(DIMV(LSPC),nmat*(1+SDIM))
      REAL_T, intent(in) :: LSHO(DIMV(LSHO),nmat*(1+SDIM))
      REAL_T, intent(out) :: curvfab(DIMV(curvfab),num_curv)
      REAL_T, intent(in) :: velfab(DIMV(velfab),SDIM)
      REAL_T, intent(in) :: denfab(DIMV(denfab),nmat*num_state_material)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T istenlo(3),istenhi(3)
      INTEGER_T LSstenlo(3),LSstenhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: bfact_grid
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: rz_flag
      REAL_T, intent(in) :: time

      INTEGER_T i,j,k
      INTEGER_T iside,jside,kside
      INTEGER_T iface,jface,kface
      INTEGER_T i1,j1,k1
      INTEGER_T ii,jj,kk
      REAL_T LS(nmat)
      REAL_T LS_fixed(nmat)
      REAL_T LSSIDE(nmat)
      REAL_T LSSIDE_fixed(nmat)
      REAL_T xcenter(SDIM)
      INTEGER_T dirloc,dircrossing,dirstar
      INTEGER_T sidestar
      INTEGER_T at_RZ_axis

      INTEGER_T im,im_opp
      INTEGER_T im_opp_test
      INTEGER_T im_curv
      INTEGER_T im_majority
      INTEGER_T im_main,im_main_opp
      INTEGER_T nten_test
      INTEGER_T iten
      INTEGER_T inormal
      REAL_T nrmPROBE(SDIM*nmat)
      REAL_T nrmFD(SDIM*nmat)
      REAL_T nrm_local(SDIM*nmat)
      REAL_T nrm_mat(SDIM)
      REAL_T nrm_test(SDIM)
      REAL_T nrm_center(SDIM)
      REAL_T curv_cellHT
      REAL_T curv_cellFD
      REAL_T pforce_cell
      REAL_T mag,RR
      INTEGER_T itemperature
      INTEGER_T donate_flag
      INTEGER_T signcrossing
      INTEGER_T signside
      INTEGER_T critsign
      REAL_T dxcrossing
      REAL_T dxside
      REAL_T mgoni_force(SDIM)

      INTEGER_T im3

      REAL_T LSCEN_hold(nmat)
      REAL_T LSCEN_hold_fixed(nmat)

      REAL_T LCEN,LSIDE
      REAL_T LSLEFT_EXTEND,LSRIGHT_EXTEND
      REAL_T LSLEFT_fixed(nmat)
      REAL_T LSRIGHT_fixed(nmat)
      REAL_T XLEFT,XRIGHT,XCEN

      REAL_T LS_STAR_FIXED(-1:1,SDIM,nmat)
      INTEGER_T im_star_majority(-1:1,SDIM)

      REAL_T velsten( &
       D_DECL(-1:1,-1:1,-1:1),SDIM)

      REAL_T nrmsten( &
       D_DECL(-1:1,-1:1,-1:1),SDIM*nmat)

      REAL_T mgoni_temp( &
       D_DECL(-1:1,-1:1,-1:1),nmat)

      REAL_T lssten( &
        D_DECL(-RD:RD,-RD:RD,-RD:RD), &
        nmat)

      REAL_T, dimension(:,:), allocatable :: xsten0
      REAL_T, dimension(:,:), allocatable :: xsten_curv

      REAL_T x1dcen,x1dside,x1dcross
      REAL_T vol_sten
      REAL_T area_sten(SDIM,2)
      INTEGER_T side_index

      INTEGER_T nhalf
      INTEGER_T RDx
      INTEGER_T RD_HEIGHT
      INTEGER_T mask1,mask2
      INTEGER_T local_mask
      INTEGER_T ihist
      REAL_T ZEYU_thet_d,ZEYU_u_cl
    
      if (ngrow_distance.ne.4) then
       print *,"ngrow_distance invalid in curvstrip"
       stop
      endif 
      if (RD.ne.4) then
       print *,"RD invalid in curvstrip"
       stop
      endif
      if (RD.ne.ngrow_distance) then
       print *,"expecting RD==ngrow_distance"
       stop
      endif
      RDx=2*RD+1
 
      nhalf=3
      allocate(xsten0(-nhalf:nhalf,SDIM))
      allocate(xsten_curv(-RDx:RDx,SDIM))
 
      if (bfact.lt.1) then
       print *,"bfact invalid90"
       stop
      endif
      if (bfact_grid.lt.4) then
       print *,"bfact_grid invalid in curvstrip"
       stop
      endif

      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid curvstrip nten nten_test ",nten,nten_test
       stop
      endif
      if (nhistory.eq.nten*2) then
       ! do nothing
      else
       print *,"nhistory invalid"
       stop
      endif

      if ((conservative_tension_force.ne.0).and. &
          (conservative_tension_force.ne.1)) then
       print *,"conservative_tension_force invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in curvstrip"
       stop
      endif
  
      if (post_restart_flag.eq.0) then 
       do dirloc=1,SDIM
        if ((fablo(dirloc)/bfact_grid)*bfact_grid.ne.fablo(dirloc)) then
         print *,"fablo mod bfact_grid not 0 in CURVSTRIP"
         stop
        endif
        if (((fabhi(dirloc)+1)/bfact_grid)*bfact_grid.ne.fabhi(dirloc)+1) then
         print *,"fabhi+1 mod bfact_grid not 0 in CURVSTRIP"
         stop
        endif
       enddo ! dirloc=1..sdim
      else if (post_restart_flag.eq.1) then
       ! do nothing
      else
       print *,"post_restart_flag invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,2896)
      call checkbound(fablo,fabhi,DIMS(LSPC),RD,-1,2897)
      call checkbound(fablo,fabhi,DIMS(LSHO),RD,-1,2898)
      call checkbound(fablo,fabhi,DIMS(masknbr),1,-1,2899)
      call checkbound(fablo,fabhi,DIMS(curvfab),1,-1,2900)
      call checkbound(fablo,fabhi,DIMS(velfab),2,-1,2901)
      call checkbound(fablo,fabhi,DIMS(denfab),2,-1,2902)

      call checkbound(fablo,fabhi,DIMS(vol),1,-1,6611)
      call checkbound(fablo,fabhi,DIMS(areax),1,0,6612)
      call checkbound(fablo,fabhi,DIMS(areay),1,1,6613)
      call checkbound(fablo,fabhi,DIMS(areaz),1,SDIM-1,6614)

      call checkbound(fablo,fabhi,DIMS(history_dat),1,-1,2902)

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid curv strip"
       print *,"ngeom_recon=",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid curv strip"
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
       print *,"rz_flag invalid in curvstrip"
       stop
      endif

      ! height function curvature
      ! finite difference curvature
      ! pforce
      ! marangoni force
      ! dir/side flag
      ! im3
      ! x nten

      if (num_curv.ne.nten*(5+SDIM)) then
       print *,"num_curv invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,1) 

       ! curvature stencil is 9x9x9
       ! loop through all cells
       !
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_mask=NINT(maskcov(D_DECL(i,j,k)))

        ! mask=tag if not covered by level+1 or outside the domain.
       if (local_mask.eq.1) then

        ! mask1=1 at interior cells or fine/fine ghost cells
        ! mask2=1 at interior cells
        mask1=NINT(masknbr(D_DECL(i,j,k),1))
        mask2=NINT(masknbr(D_DECL(i,j,k),2))

        if ((mask2.eq.1).or.(mask1.eq.0)) then

         call gridsten_level(xsten0,i,j,k,level,nhalf)

         ! center of cell
         do dirloc=1,SDIM
          xcenter(dirloc)=xsten0(0,dirloc)
         enddo

         vol_sten=vol(D_DECL(i,j,k))

         do icurv=1,num_curv
          curvfab(D_DECL(i,j,k),icurv)=zero
         enddo 
         do im=1,nmat
          LS(im)=LSHO(D_DECL(i,j,k),im)
         enddo

         call FIX_LS_tessellate(LS,LS_fixed,nmat)

         call get_primary_material(LS_fixed,nmat,im_majority)

         if (is_rigid(nmat,im_majority).eq.1) then

          ! do nothing, all interface forces are 0
  
         else if ((xsten0(0,1).le.VOFTOL*dx(1)).and. &
                  (levelrz.eq.1)) then
         
          ! do nothing, all interface forces are 0
         
         else if (is_rigid(nmat,im_majority).eq.0) then

          if (vol_sten.le.zero) then
           print *,"vol_sten invalid: cell volume should be positive"
           stop
          endif

          do dirstar=1,SDIM

           ii=0
           jj=0
           kk=0
           if (dirstar.eq.1) then
            ii=1
           else if (dirstar.eq.2) then
            jj=1
           else if ((dirstar.eq.3).and.(SDIM.eq.3)) then
            kk=1
           else
            print *,"dirstar invalid"
            stop
           endif

           do sidestar=-1,1,2

            iside=i+sidestar*ii
            jside=j+sidestar*jj
            kside=k+sidestar*kk

            do im=1,nmat
             LSSIDE(im)=LSHO(D_DECL(iside,jside,kside),im)
            enddo
            call FIX_LS_tessellate(LSSIDE,LSSIDE_fixed,nmat)
            call get_primary_material(LSSIDE_fixed,nmat,im_opp)
            do im=1,nmat
             LS_STAR_FIXED(sidestar,dirstar,im)=LSSIDE_fixed(im)
            enddo
            im_star_majority(sidestar,dirstar)=im_opp

           enddo ! sidestar=-1,1,2

          enddo ! dirstar=1..sdim

           ! loop through all possible interfaces involving im_majority
           ! and initialize curvfab
          do im_opp=1,nmat

           donate_flag=0
  
           if (im_opp.eq.im_majority) then
            ! do nothing
           else if (is_rigid(nmat,im_opp).eq.1) then
            ! do nothing
           else if (is_rigid(nmat,im_opp).eq.0) then

            if (im_majority.lt.im_opp) then
             im_main=im_majority
             im_main_opp=im_opp
            else if (im_majority.gt.im_opp) then
             im_main=im_opp
             im_main_opp=im_majority
            else
             print *,"im_majority bust"
             stop
            endif
            call get_iten(im_main,im_main_opp,iten,nmat)

            do dirloc=1,SDIM
       
             do im=1,nmat 
              LSLEFT_fixed(im)=LS_STAR_FIXED(-1,dirloc,im)
              LSRIGHT_fixed(im)=LS_STAR_FIXED(1,dirloc,im)
             enddo

             call get_LS_extend(LSLEFT_fixed,nmat,iten,LSLEFT_EXTEND)
             call get_LS_extend(LSRIGHT_fixed,nmat,iten,LSRIGHT_EXTEND)

             XRIGHT=xsten0(2,dirloc)
             XLEFT=xsten0(-2,dirloc)
             XCEN=xsten0(0,dirloc)
             if ((XRIGHT-XCEN.gt.zero).and. &
                 (XCEN-XLEFT.gt.zero)) then
              ! do nothing
             else
              print *,"position bust"
              stop
             endif
   
             nrm_center(dirloc)= &
              (LSRIGHT_EXTEND-LSLEFT_EXTEND)/(XRIGHT-XLEFT)
            enddo !dirloc=1..sdim

            ! if R-Theta, then N(2) -> N(2)/RR + renormalize.
            RR=xcenter(1)
            call prepare_normal(nrm_center,RR,mag)

            dircrossing=0
            ! 1D normal pointing towards the center cell.
            signcrossing=0
            ! abs(x1dcross-x1dcen)
            dxcrossing=zero
            im3=0

            do dirstar=1,SDIM

             x1dcen=xsten0(0,dirstar)

             ii=0
             jj=0
             kk=0
             if (dirstar.eq.1) then
              ii=1
             else if (dirstar.eq.2) then
              jj=1
             else if ((dirstar.eq.3).and.(SDIM.eq.3)) then
              kk=1
             else
              print *,"dirstar invalid"
              stop
             endif

             do sidestar=-1,1,2

              if (sidestar.eq.-1) then
               iface=i
               jface=j
               kface=k
               side_index=1
              else if (sidestar.eq.1) then
               iface=i+ii
               jface=j+jj
               kface=k+kk
               side_index=2
              else
               print *,"sidestar invalid"
               stop
              endif

              if (dirstar.eq.1) then
               area_sten(dirstar,side_index)=areax(D_DECL(iface,jface,kface))
              else if (dirstar.eq.2) then
               area_sten(dirstar,side_index)=areay(D_DECL(iface,jface,kface))
              else if ((dirstar.eq.3).and.(SDIM.eq.3)) then
               area_sten(dirstar,side_index)=areaz(D_DECL(iface,jface,kface))
              else
               print *,"dirstar invalid curvstrip"
               stop
              endif
            
              x1dside=xsten0(2*sidestar,dirstar)
  
              do im=1,nmat
               LSSIDE_fixed(im)=LS_STAR_FIXED(sidestar,dirstar,im)
              enddo
              im_opp_test=im_star_majority(sidestar,dirstar)

              if (im_opp_test.eq.im_opp) then

               at_RZ_axis=0
               if ((levelrz.eq.1).and. &
                   (xsten0(-1,1).le.VOFTOL*dx(1)).and. &
                   (dirstar.eq.1).and. &
                   (sidestar.eq.-1)) then
                at_RZ_axis=1
               endif

               LCEN=LS_fixed(im_majority)
               LSIDE=LSSIDE_fixed(im_opp)

               if ((LCEN*LSIDE.ge.zero).and. &
                   (abs(LCEN)+abs(LSIDE).gt.zero).and. &
                   (at_RZ_axis.eq.0)) then

                LCEN=-LS_fixed(im_opp)
                LSIDE=-LSSIDE_fixed(im_majority)

                if ((LCEN*LSIDE.ge.zero).and. &
                    (abs(LCEN)+abs(LSIDE).gt.zero)) then

                 call get_crossing(x1dcross,x1dcen,x1dside,LCEN,LSIDE)

                 dxside=abs(x1dcross-x1dcen)

                  ! signcrossing points to im_majority material
                 if (donate_flag.eq.0) then
                  donate_flag=1
                  dircrossing=dirstar
                  signcrossing=-sidestar
                  dxcrossing=dxside
                 else if (donate_flag.eq.1) then
                  if (dxside.lt.(one-VOFTOL)*dxcrossing) then
                   dircrossing=dirstar
                   signcrossing=-sidestar
                   dxcrossing=dxside
                  else if ((dxside.le.(one+VOFTOL)*dxcrossing).and. &
                           (abs(nrm_center(dirstar)).gt. &
                            abs(nrm_center(dircrossing)))) then 
                   dircrossing=dirstar
                   signcrossing=-sidestar
                   dxcrossing=dxside
                  endif
                 else
                  print *,"donate_flag invalid"
                  stop
                 endif

                else if (LSIDE*LCEN.le.zero) then
                 ! do nothing
                else
                 print *,"LSIDE or LCEN bust"
                 stop
                endif

               else if ((LSIDE*LCEN.le.zero).or. &
                        (at_RZ_axis.eq.1)) then
                ! do nothing
               else
                print *,"LSIDE, LCEN, or at_RZ_axis bust"
                stop
               endif 

              endif ! im_opp_test==im_opp?

             enddo ! sidestar=-1,1,2
            enddo ! dirstar=1..sdim

            if (donate_flag.eq.1) then

              ! sidestar points away from im_majority and towards the
              ! opposite material.
              ! signcrossing points towards im_majority.
             sidestar=-signcrossing

             if (im_majority.lt.im_opp) then
               ! signside points to im_majority=im_main
              signside=signcrossing
              if ((im_majority.eq.im_main).and. &
                  (im_opp.eq.im_main_opp)) then
               ! do nothing
              else
               print *,"im_majority or im_opp invalid"
               stop
              endif
             else if (im_majority.gt.im_opp) then
               ! signcrossing points towards im_majority.
               ! signside points away from im_majority
               ! signside points towards im_opp=im_main
              signside=-signcrossing
              if ((im_majority.eq.im_main_opp).and. &
                  (im_opp.eq.im_main)) then
               ! do nothing
              else
               print *,"im_majority or im_opp invalid"
               stop
              endif
             else
              print *,"im_majority bust"
              stop
             endif

             dxside=dxcrossing
             dirstar=dircrossing
             if ((dirstar.lt.1).or.(dirstar.gt.SDIM)) then
              print *,"dirstar invalid"
              stop
             endif
             if ((signside.ne.1).and.(signside.ne.-1)) then
              print *,"signside invalid"
              stop
             endif

             istenlo(3)=0
             istenhi(3)=0
             do dirloc=1,SDIM
              istenlo(dirloc)=-1
              istenhi(dirloc)=1
             enddo

             if (mask2.eq.0) then ! mask2==0 => not interior cell 
              RD_HEIGHT=ngrow_distance-1
             else if (mask2.eq.1) then ! mask2==1 => interior cell
              RD_HEIGHT=RD
             else
              print *,"mask2 invalid"
              stop
             endif

             LSstenlo(3)=0
             LSstenhi(3)=0
             do dirloc=1,SDIM
              LSstenlo(dirloc)=-RD_HEIGHT
              LSstenhi(dirloc)=RD_HEIGHT
             enddo
    
             call gridsten_level(xsten_curv,i,j,k,level,RDx)

             ! get normals at the cell center.
             do inormal=1,SDIM*nmat
              nrmPROBE(inormal)=LSHO(D_DECL(i,j,k),nmat+inormal)
             enddo ! inormal

              ! get normals at the cell center using finite differences. 
             do dirstar=1,SDIM

              XRIGHT=xsten_curv(2,dirstar)
              XLEFT=xsten_curv(-2,dirstar)
              XCEN=xsten_curv(0,dirstar)
              if ((XRIGHT-XCEN.gt.zero).and. &
                  (XCEN-XLEFT.gt.zero)) then
               ! do nothing
              else
               print *,"position bust"
               stop
              endif

              do im_curv=1,nmat
               LSRIGHT_EXTEND=LS_STAR_FIXED(1,dirstar,im_curv)
               LSLEFT_EXTEND=LS_STAR_FIXED(-1,dirstar,im_curv)
               LCEN=LS_fixed(im_curv)

               inormal=(im_curv-1)*SDIM+dirstar

               if (dirstar.eq.dircrossing) then
                if (sidestar.eq.1) then
                 nrmFD(inormal)=(LSRIGHT_EXTEND-LCEN)/(XRIGHT-XCEN)
                else if (sidestar.eq.-1) then
                 nrmFD(inormal)=(LCEN-LSLEFT_EXTEND)/(XCEN-XLEFT)
                else
                 print *,"sidestar invalid"
                 stop
                endif
               else if (dirstar.ne.dircrossing) then
                nrmFD(inormal)= &
                  (LSRIGHT_EXTEND-LSLEFT_EXTEND)/(XRIGHT-XLEFT)
               else
                print *,"dirstar invalid"
                stop
               endif
              enddo ! im_curv=1..nmat

             enddo ! dirstar=1..sdim

             do im_curv=1,nmat
              do dirloc=1,SDIM
               inormal=(im_curv-1)*SDIM+dirloc
               nrm_mat(dirloc)=nrmPROBE(inormal)
               nrm_test(dirloc)=nrmFD(inormal)
              enddo ! dirloc=1..sdim
              RR=one
              call prepare_normal(nrm_test,RR,mag)

              ! fix probe normal if it is inconsistent.
              ! iten: im_main,im_main_opp; normal points towards im_main
              if ((im_curv.eq.im_main).or. &
                  (im_curv.eq.im_main_opp)) then
  
               ! nrm_test points towards the im_curv material 
               ! signside points towards im_main. 
               if (im_curv.eq.im_main) then
                critsign=signside
               else if (im_curv.eq.im_main_opp) then
                critsign=-signside
               else
                print *,"im_curv invalid"
                stop
               endif
               if (nrm_test(dircrossing)*critsign.le.zero) then
                print *,"critsign and nrm_test mismatch"
                print *,"dircrossing = ",dircrossing
                print *,"sidestar= ",sidestar
                print *,"signcrossing= ",signcrossing
                print *,"dxcrossing= ",dxcrossing
                print *,"critsign=",critsign
                print *,"im_curv=",im_curv
                print *,"im_majority= ",im_majority
                print *,"im_main=",im_main
                print *,"im_main_opp=",im_main_opp
                print *,"iten= ",iten
                print *,"im_opp=",im_opp
                print *,"nrm_test points towards im_curv= ",im_curv
                do dirloc=1,SDIM
                 print *,"dirloc,nrm_test ",dirloc,nrm_test(dirloc)
                enddo
                print *,"nrm_center points towards im_main= ",im_main
                do dirloc=1,SDIM
                 print *,"dirloc,nrm_center ",dirloc,nrm_center(dirloc)
                enddo
                print *,"nrmPROBE points towards im_curv= ",im_curv
                do dirloc=1,SDIM
                 inormal=(im_curv-1)*SDIM+dirloc
                 print *,"dirloc,nrmPROBE ",dirloc,nrmPROBE(inormal)
                enddo
                print *,"nrmFD points towards im_curv= ",im_curv
                do dirloc=1,SDIM
                 inormal=(im_curv-1)*SDIM+dirloc
                 print *,"dirloc,nrmFD ",dirloc,nrmFD(inormal)
                enddo
                stop
               endif

                ! nrm_mat: from PROBE normal
                ! nrm_test: from FD normal
               if ((nrm_mat(dircrossing)*nrm_test(dircrossing).le.zero).or. &
                   (abs(nrm_mat(dircrossing)).le. &
                    half*abs(nrm_test(dircrossing)))) then
                do dirloc=1,SDIM
                 nrm_mat(dirloc)=nrm_test(dirloc)
                enddo
               endif
  
              endif ! im_curv=im_main or im_main_opp

               ! if R-Theta, then N(2) -> N(2)/RR + renormalize.
              RR=xcenter(1)
              call prepare_normal(nrm_mat,RR,mag)
              do dirloc=1,SDIM
               inormal=(im_curv-1)*SDIM+dirloc
               nrmPROBE(inormal)=nrm_mat(dirloc)
              enddo

             enddo ! im_curv=1..nmat

             if (1.eq.0) then
              print *,"xcenter ",xcenter(1),xcenter(2),xcenter(SDIM)
              print *,"dircrossing ",dircrossing
              print *,"im_majority,im_opp,im_main,im_main_opp ", &
               im_majority,im_opp,im_main,im_main_opp
             endif

             ! i1,j1,k1=-RD_HEIGHT..RD_HEIGHT
             do i1=LSstenlo(1),LSstenhi(1)
             do j1=LSstenlo(2),LSstenhi(2)
             do k1=LSstenlo(3),LSstenhi(3)

              if ((abs(i1).le.1).and.(abs(j1).le.1).and.(abs(k1).le.1)) then
               do inormal=1,SDIM*nmat
                nrm_local(inormal)=LSHO(D_DECL(i+i1,j+j1,k+k1),nmat+inormal)
               enddo
              else if ((abs(i1).le.RD_HEIGHT).and. &
                       (abs(j1).le.RD_HEIGHT).and. &
                       (abs(k1).le.RD_HEIGHT)) then
               ! do nothing
              else
               print *,"i1,j1, or k1 invalid"
               stop
              endif

              do im_curv=1,nmat
               LSCEN_hold(im_curv)=LSHO(D_DECL(i+i1,j+j1,k+k1),im_curv)
              enddo
              call FIX_LS_tessellate(LSCEN_hold,LSCEN_hold_fixed,nmat)
 
              do im_curv=1,nmat

               lssten(D_DECL(i1,j1,k1),im_curv)=LSCEN_hold_fixed(im_curv)

               if ((abs(i1).le.1).and.(abs(j1).le.1).and.(abs(k1).le.1)) then
                do dirloc=1,SDIM
                 inormal=(im_curv-1)*SDIM+dirloc
                 nrm_mat(dirloc)=nrm_local(inormal)
                enddo
                RR=one
                if (levelrz.eq.0) then
                 ! do nothing
                else if (levelrz.eq.1) then
                 if (SDIM.ne.2) then
                  print *,"levelrz invalid"
                  stop
                 endif
                else if (levelrz.eq.3) then
                 RR=xsten_curv(2*i1,1)
                else
                 print *,"transformed normal: levelrz invalid"
                 stop
                endif
                call prepare_normal(nrm_mat,RR,mag)
                do dirloc=1,SDIM
                 inormal=(im_curv-1)*SDIM+dirloc
                 nrmsten(D_DECL(i1,j1,k1),inormal)=nrm_mat(dirloc)
                enddo
               else if ((abs(i1).le.RD_HEIGHT).and. &
                        (abs(j1).le.RD_HEIGHT).and. &
                        (abs(k1).le.RD_HEIGHT)) then
                ! do nothing
               else
                print *,"i1,j1, or k1 invalid"
                stop
               endif

              enddo ! im_curv=1..nmat

             enddo
             enddo
             enddo ! i1,j1,k1 (init nrmsten and lssten)
 
             ! i1,j1,k1=-1..1
             do i1=istenlo(1),istenhi(1)
             do j1=istenlo(2),istenhi(2)
             do k1=istenlo(3),istenhi(3)

              do dirloc=1,SDIM
               velsten(D_DECL(i1,j1,k1),dirloc)= &
                velfab(D_DECL(i+i1,j+j1,k+k1),dirloc)
              enddo

              do im_curv=1,nmat
               itemperature=(im_curv-1)*num_state_material+2
               mgoni_temp(D_DECL(i1,j1,k1),im_curv)= &
                denfab(D_DECL(i+i1,j+j1,k+k1),itemperature)
              enddo
  
             enddo
             enddo
             enddo  ! i1,j1,k1 (init velsten)

             ! tension used to find contact angle (scaling not 
             ! necessary)
             call initheightLS( &
              conservative_tension_force, &
              i,j,k, &
              level, &
              finest_level, &
              bfact,dx, &
              xcenter, &
              nrmPROBE, & ! nmat x sdim components
              dircrossing, &
              sidestar, &
              signside, &
              time, &
              xsten_curv, &
              velsten, &
              mgoni_temp, &
              lssten, &
              nrmsten, &
              vol_sten, &
              area_sten, &
              curv_cellHT, &
              curv_cellFD, &
              mgoni_force, & !(I-nn^T)(grad sigma) delta
              ZEYU_thet_d, &
              ZEYU_u_cl, &
              im3, &
              nmat, &
              visc_coef,nten, &
              im_main,im_main_opp,iten, &
              RD,RDx,RD_HEIGHT)

             if (1.eq.0) then
              print *,"i,j,k,dircrossing,sidestar,nrm ", &
               i,j,k,dircrossing,sidestar, &
               nrm_center(1),nrm_center(2),nrm_center(SDIM)
              print *,"RD,RDx,RD_HEIGHT ",RD,RDx,RD_HEIGHT
              print *,"curv_cellHT,curv_cellFD ",curv_cellHT,curv_cellFD
             endif
          
             call initpforce( &
              bfact,dx,xsten0, &
              RD,RDx,RD_HEIGHT, &
              time, &
              lssten, &
              dircrossing-1,nmat, &
              im_main,im_main_opp, &
              pforce_cell)

             ihist=(iten-1)*2
             history_dat(D_DECL(i,j,k),ihist+1)=ZEYU_thet_d
             history_dat(D_DECL(i,j,k),ihist+2)=ZEYU_u_cl

             icurv=(iten-1)*(5+SDIM)
             curvfab(D_DECL(i,j,k),icurv+1)=curv_cellHT
             curvfab(D_DECL(i,j,k),icurv+2)=curv_cellFD
             curvfab(D_DECL(i,j,k),icurv+3)=pforce_cell
             do dirloc=1,SDIM
              curvfab(D_DECL(i,j,k),icurv+3+dirloc)=mgoni_force(dirloc)
             enddo
              ! dir=1..sdim
              ! side=-1 or 1
             curvfab(D_DECL(i,j,k),icurv+4+SDIM)=dircrossing*sidestar
             curvfab(D_DECL(i,j,k),icurv+5+SDIM)=im3

             if (curv_min.gt.curv_cellHT) then
              curv_min=curv_cellHT
             endif
             if (curv_max.lt.curv_cellHT) then
              curv_max=curv_cellHT
             endif

            else if (donate_flag.eq.0) then
             ! do nothing
            else
             print *,"donate_flag invalid"
             stop
            endif

           else
            print *,"im_opp invalid"
            stop
           endif

          enddo ! im_opp

         else
          print *,"im_majority invalid"
          stop
         endif
  
        else if ((mask2.eq.0).and.(mask1.eq.1)) then
         ! do nothing
        else
         print *,"mask2 or mask1 invalid"
         stop
        endif

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo
      enddo
      enddo  ! i,j,k


      deallocate(xsten0)
      deallocate(xsten_curv)

      return
      end subroutine FORT_CURVSTRIP

      subroutine FORT_GETTYPEFAB( &
       LS,DIMS(LS), &
       typefab,DIMS(typefab), &
       xlo,dx, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       type_flag,nmat)
      use probf90_module
      use global_utility_module
 

      IMPLICIT NONE

      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(out) :: type_flag(nmat)

      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact

      INTEGER_T, intent(in) ::  DIMDEC(LS)
      INTEGER_T, intent(in) ::  DIMDEC(typefab)

      REAL_T, intent(in) :: LS(DIMV(LS),nmat)
      REAL_T, intent(out) :: typefab(DIMV(typefab))
      INTEGER_T i,j,k,im,base_type


      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid91"
       stop
      endif

      call checkbound(fablo,fabhi, &
       DIMS(LS), &
       1,-1,4001)
      call checkbound(fablo,fabhi, &
       DIMS(typefab), &
       1,-1,4001)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       base_type=1
       do im=2,nmat
        if (LS(D_DECL(i,j,k),im).gt. &
            LS(D_DECL(i,j,k),base_type)) then
         base_type=im
        endif
       enddo

       do im=1,nmat
        if (is_rigid(nmat,im).eq.1) then
         if (LS(D_DECL(i,j,k),im).ge.zero) then
          base_type=im
         endif
        else if (is_rigid(nmat,im).eq.0) then
         ! do nothing
        else
         print *,"is_rigid invalid"
         stop
        endif
       enddo ! im=1..nmat

       typefab(D_DECL(i,j,k))=base_type
       if ((base_type.gt.nmat).or.(base_type.lt.1)) then
        print *,"base_type invalid"
        stop
       else
        type_flag(base_type)=1
       endif
      enddo
      enddo
      enddo

      return
      end subroutine FORT_GETTYPEFAB

      subroutine FORT_GETCOLORSUM( &
       sweep_num, &
       dx,xlo, &
       nmat, &
       LS,DIMS(LS), &
       VEL,DIMS(VEL), &
       DEN,DIMS(DEN), &
       VOF,DIMS(VOF), &
       facefab,DIMS(facefab), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz), &
       cellfab,DIMS(cellfab), &
       typefab,DIMS(typefab), &
       color,DIMS(color), &
       mask,DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level, &
       finest_level, &
       rzflag, &
       num_colors, &
       cum_blobdata, &
       level_blobdata, &
       level_blobtypedata, &
       arraysize,num_elements_blobclass, &
       levelbc, &
       velbc, &
       nface,nface_dst,ncellfrac)
      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: sweep_num
      INTEGER_T, intent(in) :: nface,nface_dst,ncellfrac
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo(SDIM)
      INTEGER_T, intent(in) :: levelbc(SDIM,2)
      INTEGER_T, intent(in) :: velbc(SDIM,2,SDIM)
      INTEGER_T, intent(in) :: num_elements_blobclass

      INTEGER_T :: i,j,k
      INTEGER_T :: ii,jj,kk
      INTEGER_T :: iface,jface,kface
      INTEGER_T :: face_index
 
      INTEGER_T, intent(in) :: rzflag,num_colors,arraysize
      INTEGER_T, intent(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(VEL)
      INTEGER_T, intent(in) :: DIMDEC(DEN)
      INTEGER_T, intent(in) :: DIMDEC(VOF)
      INTEGER_T, intent(in) :: DIMDEC(facefab)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(areax)
      INTEGER_T, intent(in) :: DIMDEC(areay)
      INTEGER_T, intent(in) :: DIMDEC(areaz)
      INTEGER_T, intent(in) :: DIMDEC(cellfab)
      INTEGER_T, intent(in) :: DIMDEC(typefab)
      INTEGER_T, intent(in) :: DIMDEC(color)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      REAL_T, intent(inout) :: level_blobdata(arraysize)
      REAL_T, intent(in) :: cum_blobdata(arraysize)
      INTEGER_T, intent(inout) :: level_blobtypedata(num_colors)

      REAL_T, intent(in) :: typefab(DIMV(typefab))
      REAL_T, intent(in) :: LS(DIMV(LS),nmat*(1+SDIM))
      REAL_T, intent(in) :: VEL(DIMV(VEL),SDIM)
      REAL_T, intent(in) :: DEN(DIMV(DEN),nmat*num_state_material)
      REAL_T, intent(in) :: VOF(DIMV(VOF),nmat*ngeom_recon)
      REAL_T, intent(in) :: facefab(DIMV(facefab),nface)
      REAL_T, intent(in) :: xface(DIMV(xface),nface_dst)
      REAL_T, intent(in) :: yface(DIMV(yface),nface_dst)
      REAL_T, intent(in) :: zface(DIMV(zface),nface_dst)
      REAL_T, intent(in) :: areax(DIMV(areax))
      REAL_T, intent(in) :: areay(DIMV(areay))
      REAL_T, intent(in) :: areaz(DIMV(areaz))
      REAL_T, intent(in) :: cellfab(DIMV(cellfab),ncellfrac)
      REAL_T, intent(in) :: color(DIMV(color))
      REAL_T, intent(in) :: mask(DIMV(mask))

      INTEGER_T dir,side
      INTEGER_T dir2
      INTEGER_T vofcomp
      INTEGER_T dencomp
      INTEGER_T base_type
      INTEGER_T opposite_color(nmat)
      INTEGER_T typeside
      INTEGER_T colorside
      INTEGER_T ic
      INTEGER_T ic_center
      INTEGER_T ic_base
      INTEGER_T icolor
      REAL_T local_facearea(nmat,nmat)
      REAL_T local_dist_to_line(nmat,nmat)
      REAL_T local_normal(nmat,nmat,SDIM)
      REAL_T local_dist(nmat,nmat)
      REAL_T vfrac
      REAL_T vol
      REAL_T mass,den_mat
      REAL_T xsten(-3:3,SDIM)
      REAL_T dx_sten(SDIM)
      INTEGER_T nhalf
      REAL_T cencell(SDIM)
      INTEGER_T i1,j1,k1
      INTEGER_T k1lo,k1hi
      INTEGER_T im,im_opp
      INTEGER_T im1,im2
      INTEGER_T ml,mr
      INTEGER_T local_mask
      REAL_T frac_pair(nmat,nmat)
      REAL_T dist_pair(nmat,nmat)
      REAL_T areaface
      INTEGER_T LS_change_sign(nmat)
      INTEGER_T LS_plus(nmat)
      INTEGER_T LS_minus(nmat)
      REAL_T LScen(nmat)
      REAL_T LSside(nmat)
      REAL_T RR,dperim
      REAL_T fluid_velocity(SDIM)
      REAL_T solid_velocity(SDIM)
      REAL_T solid_fraction
      REAL_T dotprod
      INTEGER_T local_solid
      INTEGER_T im_side_majority
      REAL_T blob_center_actual(3)
      REAL_T blob_x0(3)
      REAL_T phi_row(3)
      REAL_T phi_col(3)
      INTEGER_T irow,icol,veltype
      REAL_T DXMAXLS,cutoff
      REAL_T im_interior_wt(3)

      nhalf=3

      if (bfact.lt.1) then
       print *,"bfact invalid92"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nface.ne.nmat*SDIM*2*(1+SDIM)) then
       print *,"nface invalid"
       stop
      endif
      if (nface_dst.ne.nmat*nmat*2) then
       print *,"nface_dst invalid"
       stop
      endif
      if (ncellfrac.ne.nmat*nmat*(3+SDIM)) then
       print *,"ncellfrac invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid get color sum"
       stop
      endif

       ! in: FORT_GETCOLORSUM
      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=two*DXMAXLS

      call checkbound(fablo,fabhi,DIMS(LS),1,-1,6615)
      call checkbound(fablo,fabhi,DIMS(VEL),1,-1,6615)
      call checkbound(fablo,fabhi,DIMS(DEN),1,-1,6615)
      call checkbound(fablo,fabhi,DIMS(VOF),1,-1,6616)
      call checkbound(fablo,fabhi,DIMS(facefab),1,-1,6617)
      call checkbound(fablo,fabhi,DIMS(xface),0,0,6618)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,6619)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,6620)
      call checkbound(fablo,fabhi,DIMS(areax),0,0,6621)
      call checkbound(fablo,fabhi,DIMS(areay),0,1,6622)
      call checkbound(fablo,fabhi,DIMS(areaz),0,SDIM-1,6623)
      call checkbound(fablo,fabhi,DIMS(cellfab),0,-1,6624)
      call checkbound(fablo,fabhi,DIMS(typefab),1,-1,6625)
      call checkbound(fablo,fabhi,DIMS(color),1,-1,6626)
      call checkbound(fablo,fabhi,DIMS(mask),1,-1,6627)
  
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
       print *,"num_elements_blobclass invalid"
       stop
      endif
      if (arraysize.ne.num_elements_blobclass*num_colors) then
       print *,"arraysize invalid"
       stop
      endif

      if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else
       print *,"dimension bust"
       stop
      endif

      do dir=1,SDIM
       if (fabhi(dir)-fablo(dir).le.0) then
        print *,"fablo,fabhi violates blocking factor"
        stop
       endif
      enddo

! color>0 if mask=1 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       local_mask=NINT(mask(D_DECL(i,j,k)))
       if (local_mask.eq.1) then
        icolor=NINT(color(D_DECL(i,j,k)))
        if ((icolor.gt.num_colors).or.(icolor.le.0)) then
         print *,"icolor invalid in GETCOLORSUM icolor=",icolor
         print *,"i,j,k ",i,j,k
         stop
        endif
        base_type=NINT(typefab(D_DECL(i,j,k)))
        if ((base_type.lt.1).or.(base_type.gt.nmat)) then
         print *,"base_type invalid"
         stop
        endif
        if (level_blobtypedata(icolor).eq.0) then
         level_blobtypedata(icolor)=base_type
        else
         if (level_blobtypedata(icolor).ne.base_type) then
          print *,"type problems in GETCOLORSUM"
          print *,"level,finest_level ",level,finest_level
          print *,"num_colors= ",num_colors
          print *,"current blobtype ",level_blobtypedata(icolor)
          print *,"base_type ",base_type
          print *,"icolor=",icolor
          print *,"i,j,k ",i,j,k
          print *,"growlo,growhi ", &
           growlo(1),growlo(2),growlo(3), &
           growhi(1),growhi(2),growhi(3)
          stop
         endif
        endif

        do im=1,nmat
         opposite_color(im)=0
        enddo
        opposite_color(base_type)=icolor

        do im=1,nmat
         LS_change_sign(im)=0
         LS_plus(im)=0
         LS_minus(im)=0
        enddo

        do im=1,nmat
         LScen(im)=LS(D_DECL(i,j,k),im)
        enddo

        do dir=1,SDIM
         fluid_velocity(dir)=VEL(D_DECL(i,j,k),dir)
        enddo

        do dir=1,SDIM
         solid_velocity(dir)=zero
        enddo
        solid_fraction=zero

        do i1=-1,1
        do j1=-1,1
        do k1=k1lo,k1hi

         local_solid=0

         do im=1,nmat
          LSside(im)=LS(D_DECL(i+i1,j+j1,k+k1),im)
          if (LScen(im)+LSside(im).ge.zero) then
           LS_plus(im)=1
          else if (LScen(im)+LSside(im).lt.zero) then
           LS_minus(im)=1
          else
           print *,"LScen or LSside invalid"
           stop
          endif
         enddo ! im=1..nmat

         call get_primary_material(LSside,nmat,im_side_majority)
         if (is_rigid(nmat,im_side_majority).eq.1) then
          local_solid=1
         else if (is_rigid(nmat,im_side_majority).eq.0) then
          ! do nothing
         else
          print *,"is_rigid(nmat,im_side_majority) invalid"
          stop
         endif

         typeside=NINT(typefab(D_DECL(i+i1,j+j1,k+k1)))
         colorside=NINT(color(D_DECL(i+i1,j+j1,k+k1)))

          ! we do not consider neighbor cells outside the physical domain
          ! for incrementing volume, surface area, or contact line perimeter.
         dir=1
         side=1
         if (i+i1.lt.fablo(dir)) then
          if (levelbc(dir,side).ne.INT_DIR) then
           typeside=0
          endif
          if (velbc(dir,side,dir).eq.EXT_DIR) then
           local_solid=1
          endif
         endif
         side=2
         if (i+i1.gt.fabhi(dir)) then
          if (levelbc(dir,side).ne.INT_DIR) then
           typeside=0
          endif
          if (velbc(dir,side,dir).eq.EXT_DIR) then
           local_solid=1
          endif
         endif

         dir=2
         side=1
         if (j+j1.lt.fablo(dir)) then
          if (levelbc(dir,side).ne.INT_DIR) then
           typeside=0
          endif
          if (velbc(dir,side,dir).eq.EXT_DIR) then
           local_solid=1
          endif
         endif
         side=2
         if (j+j1.gt.fabhi(dir)) then
          if (levelbc(dir,side).ne.INT_DIR) then
           typeside=0
          endif
          if (velbc(dir,side,dir).eq.EXT_DIR) then
           local_solid=1
          endif
         endif

         if (SDIM.eq.3) then
          dir=SDIM
          side=1
          if (k+k1.lt.fablo(dir)) then
           if (levelbc(dir,side).ne.INT_DIR) then
            typeside=0
           endif
           if (velbc(dir,side,dir).eq.EXT_DIR) then
            local_solid=1
           endif
          endif
          side=2
          if (k+k1.gt.fabhi(dir)) then
           if (levelbc(dir,side).ne.INT_DIR) then
            typeside=0
           endif
           if (velbc(dir,side,dir).eq.EXT_DIR) then
            local_solid=1
           endif
          endif
         else if (SDIM.eq.2) then
          ! do nothing
         else
          print *,"dimension bust"
          stop
         endif

         if (typeside.eq.0) then
          ! do nothing
         else if (typeside.eq.base_type) then
           ! do nothing
         else if ((typeside.ge.1).and.(typeside.le.nmat)) then
          if (colorside.eq.0) then
           ! do nothing
          else if ((colorside.gt.0).and. &
                   (colorside.le.num_colors)) then
           opposite_color(typeside)=colorside
          else
           print *,"colorside invalid in getcolorsum"
           print *,"colorside= ",colorside
           print *,"level,finest_level ",level,finest_level
           stop
          endif
         else
          print *,"typeside invalid"
          stop
         endif

         if (local_solid.eq.1) then
          solid_fraction=solid_fraction+one
          do dir=1,SDIM
           solid_velocity(dir)=solid_velocity(dir)+ &
            VEL(D_DECL(i+i1,j+j1,k+k1),dir)
          enddo
         else if (local_solid.eq.0) then
          ! do nothing
         else
          print *,"local_solid invalid"
          stop
         endif

        enddo
        enddo
        enddo ! i1,j1,k1

        if (solid_fraction.eq.zero) then
         ! do nothing
        else if (solid_fraction.gt.zero) then
         do dir=1,SDIM
          solid_velocity(dir)=solid_velocity(dir)/solid_fraction
         enddo
         solid_fraction=one
        else
         print *,"solid_fraction invalid"
         stop
        endif

        call gridsten_level(xsten,i,j,k,level,nhalf)
        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
         vol,cencell,SDIM)
        mass=zero
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         vfrac=VOF(D_DECL(i,j,k),vofcomp)
         if (vfrac.ge.-VOFTOL) then
          if (vfrac.le.one+VOFTOL) then
           if (vfrac.lt.zero) then
            vfrac=zero
           endif
           if (vfrac.gt.one) then
            vfrac=one
           endif
           dencomp=(im-1)*num_state_material+1
           den_mat=DEN(D_DECL(i,j,k),dencomp)
           if (den_mat.ge.(one-VOFTOL)*fort_density_floor(im)) then
            if (den_mat.le.(one+VOFTOL)*fort_density_ceiling(im)) then
             mass=mass+den_mat*vfrac
            else
             print *,"den_mat overflow"
             print *,"den_mat= ",den_mat
             print *,"fort_density_ceiling(im)=",fort_density_ceiling(im)
             stop
            endif
           else
            print *,"den_mat underflow"
            stop
           endif
          else
           print *,"vfrac overflow"
           stop
          endif
         else
          print *,"vfrac underflow"
          stop
         endif
        enddo ! im=1..nmat

        mass=mass*vol
        if ((mass.le.zero).or.(mass.gt.1.0D+20)) then
         print *,"mass: floating point bust"
         stop
        endif

        if ((base_type.ge.1).and. &
            (base_type.le.nmat)) then

         do im=1,nmat
          if ((LS_plus(im).eq.1).and.(LS_minus(im).eq.1)) then
           LS_change_sign(im)=1
          else if ((LS_plus(im).eq.1).and.(LS_minus(im).eq.0)) then
           LS_change_sign(im)=0
          else if ((LS_plus(im).eq.0).and.(LS_minus(im).eq.1)) then
           LS_change_sign(im)=0
          else
           print *,"LS_plus or LS_minus invalid"
           stop
          endif
         enddo ! im=1..nmat

         do dir=1,SDIM
          dx_sten(dir)=xsten(1,dir)-xsten(-1,dir)
          if (dx_sten(dir).le.zero) then
           print *,"dx_sten invalid"
           stop
          endif
         enddo
         RR=xsten(0,1)
         if (levelrz.eq.0) then
          ! do nothing
         else if (levelrz.eq.1) then
          ! do nothing
         else if (levelrz.eq.3) then
          dx_sten(2)=dx_sten(2)*RR
          if (RR.le.zero) then
           print *,"RR invalid"
           stop
          endif
         else
          print *,"levelrz invalid"
          stop
         endif

         if (SDIM.eq.2) then
          if (levelrz.eq.0) then
           dperim=one
          else if (levelrz.eq.1) then
           dperim=two*Pi*RR
          else if (levelrz.eq.3) then
           dperim=one
          else
           print *,"levelrz invalid"
           stop
          endif
         else if (SDIM.eq.3) then
          dperim=sqrt(dx_sten(1)**2+dx_sten(2)**2+dx_sten(SDIM)**2)
         else
          print *,"dimension bust"
          stop
         endif
 
         ! init temporary cell data: local_facearea, loca_dist_to_line,
         !   local_dist, and local_normal
         face_index=0
         do im=1,nmat
          do im_opp=1,nmat
           face_index=face_index+1
            ! the following quantities were found in FORT_CELLFACEINIT
            ! face area between im and im_opp
           local_facearea(im,im_opp)=cellfab(D_DECL(i,j,k),face_index)
           face_index=face_index+1
            ! distance from centroid(im) to the im/im_opp line
            ! (im_opp,im) is distance from centroid(im_opp) to the
            ! im/im_opp line.
           local_dist_to_line(im,im_opp)=cellfab(D_DECL(i,j,k),face_index)
           face_index=face_index+1
            ! distance between centroid(im) and centroid(im_opp)
           local_dist(im,im_opp)=cellfab(D_DECL(i,j,k),face_index)
           do dir2=1,SDIM
            face_index=face_index+1
            local_normal(im,im_opp,dir2)=cellfab(D_DECL(i,j,k),face_index)
           enddo ! dir2
          enddo ! im_opp
         enddo ! im
           
         if (face_index.ne.ncellfrac) then
          print *,"ncellfrac invalid"
          stop
         endif

         do im=1,nmat
          if ((opposite_color(im).ge.1).and. &
              (opposite_color(im).le.num_colors)) then

           ic_center=(opposite_color(im)-1)*num_elements_blobclass+ &
            3*(2*SDIM)*(2*SDIM)+3*(2*SDIM)+3*(2*SDIM)+ &
            2*(2*SDIM)+1+ &
            3+1+SDIM+1

           do dir=1,3
            blob_center_actual(dir)=zero
            blob_x0(dir)=zero
           enddo

           do dir=1,SDIM
            blob_center_actual(dir)=cum_blobdata(ic_center)
            ic_center=ic_center+1
            blob_x0(dir)=xsten(0,dir)
           enddo ! dir=1..sdim

            ! index for first component associated with material im
           ic_base=(opposite_color(im)-1)*num_elements_blobclass
           ic=ic_base+1

           if (LScen(im).ge.cutoff) then ! cutoff=2 * DXMAXLS
            im_interior_wt(1)=one
           else if (LScen(im).lt.cutoff) then
            im_interior_wt(1)=1.0E-3
           else
            print *,"LScen(im) invalid"
            stop
           endif

           im_interior_wt(2)=solid_fraction

           if (LScen(im).ge.zero) then
            im_interior_wt(3)=one
           else if (LScen(im).lt.zero) then
            im_interior_wt(3)=1.0E-3
           else
            print *,"LScen(im) invalid"
            stop
           endif

            ! blob_matrix
           do veltype=1,3
            do irow=1,2*SDIM
             do icol=1,2*SDIM
              call init_basis(blob_center_actual, &
               blob_x0,irow,phi_row)
              call init_basis(blob_center_actual, &
               blob_x0,icol,phi_col)
              dotprod=zero
              do dir=1,SDIM
               dotprod=dotprod+phi_row(dir)*phi_col(dir)
              enddo
              dotprod=dotprod*im_interior_wt(veltype)
              level_blobdata(ic)=level_blobdata(ic)+mass*dotprod
              ic=ic+1
             enddo ! icol
            enddo ! irow
           enddo ! veltype=1..3 (blob_matrix)

            ! blob_RHS
           do veltype=1,3
            do irow=1,2*SDIM
             call init_basis(blob_center_actual, &
              blob_x0,irow,phi_row)
             dotprod=zero
             do dir=1,SDIM
              if (veltype.eq.1) then
               dotprod=dotprod+phi_row(dir)*fluid_velocity(dir)
              else if (veltype.eq.2) then
               dotprod=dotprod+phi_row(dir)*solid_velocity(dir)
              else if (veltype.eq.3) then
               dotprod=dotprod+phi_row(dir)*fluid_velocity(dir)
              else
               print *,"veltype invalid"
               stop
              endif
             enddo ! dir=1..sdim

             dotprod=dotprod*im_interior_wt(veltype)

             level_blobdata(ic)=level_blobdata(ic)+mass*dotprod
             ic=ic+1
            enddo ! irow=1..2*sdim
           enddo ! veltype=1..3

           ic=ic+3*(2*SDIM) ! skip over blob_velocity 

            ! blob_integral_momentum  2 * (2 * sdim) components
            ! first group: momentum
            ! second group: mass for momentum
           veltype=3
           do irow=1,2*SDIM
            call init_basis(blob_center_actual, &
             blob_x0,irow,phi_row)
            dotprod=zero
            do dir=1,SDIM
             dotprod=dotprod+phi_row(dir)*fluid_velocity(dir)
            enddo
            dotprod=dotprod*im_interior_wt(veltype)
            level_blobdata(ic)=level_blobdata(ic)+mass*dotprod
            ic=ic+1
           enddo ! irow=1..2 * sdim

           do irow=1,2*SDIM
            call init_basis(blob_center_actual, &
             blob_x0,irow,phi_row)
            dotprod=zero
            do dir=1,SDIM
             dotprod=dotprod+phi_row(dir)**2
            enddo
            dotprod=dotprod*im_interior_wt(veltype)
            level_blobdata(ic)=level_blobdata(ic)+mass*dotprod
            ic=ic+1
           enddo ! irow=1..2 * sdim

            ! blob_energy
           dotprod=zero
           do dir=1,SDIM
            dotprod=dotprod+fluid_velocity(dir)**2
           enddo
           dotprod=dotprod*im_interior_wt(veltype)
           level_blobdata(ic)=level_blobdata(ic)+half*mass*dotprod
           ic=ic+1

            ! blob_mass_for_velocity
           do veltype=1,3
            dotprod=im_interior_wt(veltype)
            level_blobdata(ic)=level_blobdata(ic)+mass*dotprod
            ic=ic+1
           enddo ! veltype=1..3

           if (ic.ne.ic_base+3*(2*SDIM)*(2*SDIM)+ &
               3*(2*SDIM)+3*(2*SDIM)+2*(2*SDIM)+1+3+1) then
            print *,"ic invalid ic=",ic
            print *,"expecting ic=", &
             ic_base+3*(2*SDIM)*(2*SDIM)+ &
               3*(2*SDIM)+3*(2*SDIM)+2*(2*SDIM)+1+3+1
            stop
           endif

            ! blob_volume
           vofcomp=(im-1)*ngeom_recon+1
           vfrac=VOF(D_DECL(i,j,k),vofcomp)
           level_blobdata(ic)=level_blobdata(ic)+vol*vfrac

            ! centroid integral
           ic=ic+1
           do dir=1,SDIM
            level_blobdata(ic)=level_blobdata(ic)+ &
             (cencell(dir)+VOF(D_DECL(i,j,k),vofcomp+dir))*vfrac*vol
            ic=ic+1
           enddo

           ! centroid actual
           ic=ic+SDIM
           
            ! perimeter (internal faces)
           do im_opp=1,nmat
            if (im_opp.ne.im) then
             level_blobdata(ic)=level_blobdata(ic)+local_facearea(im,im_opp) 
            endif
           enddo

            ! perimeter decomposed by material.
           ic=ic+1
           do im_opp=1,nmat
            if (im_opp.ne.im) then
             level_blobdata(ic)=level_blobdata(ic)+local_facearea(im,im_opp) 
            endif
            ic=ic+1
           enddo ! im_opp

            ! contact line perimeter
           do im1=1,nmat
           do im2=1,nmat
            if ((im1.eq.im).or.(im2.eq.im)) then
             ! do nothing
            else if (im1.eq.im2) then
             ! do nothing
            else if ((im1.ne.im).and. &
                     (im2.ne.im).and. &
                     (im1.ne.im2)) then
             if ((LS_change_sign(im).eq.1).and. &
                 (LS_change_sign(im1).eq.1).and. &
                 (LS_change_sign(im2).eq.1)) then
              level_blobdata(ic)=level_blobdata(ic)+dperim
             else if ((LS_change_sign(im).eq.0).or. &
                      (LS_change_sign(im1).eq.0).or. &
                      (LS_change_sign(im2).eq.0)) then
              ! do nothing
             else
              print *,"LS_change_sign invalid"
              stop
             endif 
            else
             print *,"im1 or im2 invalid"
             stop
            endif
            ic=ic+1
           enddo ! im2=1..nmat
           enddo ! im1=1..nmat

           if (ic.ne.opposite_color(im)*num_elements_blobclass+1) then
            print *,"ic invalid in getcolorsum"
            print *,"ic=",ic
            print *,"im=",im
            print *,"opposite_color(im)=",opposite_color(im)
            print *,"num_elements_blobclass=",num_elements_blobclass
            stop
           endif

            ! perimeter (cell faces)
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
             print *,"dir invalid getcolorsum"
             stop
            endif

            do side=1,2
             iface=i
             jface=j
             kface=k
             if (side.eq.1) then
              ! do nothing
             else if (side.eq.2) then
              iface=i+ii
              jface=j+jj
              kface=k+kk
             else
              print *,"side invalid"
              stop
             endif
             if (dir.eq.1) then
              areaface=areax(D_DECL(iface,jface,kface))
              face_index=0
              do ml = 1, nmat
              do mr = 1, nmat
               !(ml,mr,2) 
               face_index=face_index+1
               frac_pair(ml,mr)=xface(D_DECL(iface,jface,kface),face_index)
               face_index=face_index+1
               dist_pair(ml,mr)=xface(D_DECL(iface,jface,kface),face_index)
              enddo ! mr
              enddo ! ml
              if (face_index.ne.nface_dst) then
               print *,"face_index invalid"
               stop
              endif
             else if (dir.eq.2) then
              areaface=areay(D_DECL(iface,jface,kface))
              face_index=0
              do ml = 1, nmat
              do mr = 1, nmat
               !(ml,mr,2) 
               face_index=face_index+1
               frac_pair(ml,mr)=yface(D_DECL(iface,jface,kface),face_index)
               face_index=face_index+1
               dist_pair(ml,mr)=yface(D_DECL(iface,jface,kface),face_index)
              enddo ! mr
              enddo ! ml
              if (face_index.ne.nface_dst) then
               print *,"face_index invalid"
               stop
              endif
             else if ((dir.eq.3).and.(SDIM.eq.3)) then
              areaface=areaz(D_DECL(iface,jface,kface))
              face_index=0
              do ml = 1, nmat
              do mr = 1, nmat
               !(ml,mr,2) 
               face_index=face_index+1
               frac_pair(ml,mr)=zface(D_DECL(iface,jface,kface),face_index)
               face_index=face_index+1
               dist_pair(ml,mr)=zface(D_DECL(iface,jface,kface),face_index)
              enddo ! mr
              enddo ! ml
              if (face_index.ne.nface_dst) then
               print *,"face_index invalid"
               stop
              endif
             else 
              print *,"dir invalid getcolorsum"
              stop
             endif

              ! F,CEN_INTEGRATE,CEN_ACTUAL,AREA,AREA(im_opp)
             ic=(opposite_color(im)-1)*num_elements_blobclass+ &
              3*(2*SDIM)*(2*SDIM)+ &
              3*(2*SDIM)+3*(2*SDIM)+2*(2*SDIM)+1+3+1+2*SDIM+1
             do im_opp=1,nmat
              if (im_opp.ne.im) then
               if (side.eq.1) then
                ml=im_opp
                mr=im
               else if (side.eq.2) then
                mr=im_opp
                ml=im
               else
                print *,"side invalid"
                stop
               endif
                ! overall perimeter
               level_blobdata(ic)=level_blobdata(ic)+frac_pair(ml,mr)*areaface 
              endif
             enddo !im_opp=1..nmat

              ! perimeter decomposed (nmat components)
             ic=ic+1
             do im_opp=1,nmat
              if (im_opp.ne.im) then
               if (side.eq.1) then
                ml=im_opp
                mr=im
               else if (side.eq.2) then
                mr=im_opp
                ml=im
               else
                print *,"side invalid"
                stop
               endif
               level_blobdata(ic)=level_blobdata(ic)+frac_pair(ml,mr)*areaface 
              endif
              ic=ic+1
             enddo ! im_opp=1..nmat
             ic=ic+nmat*nmat
 
             if (ic.ne.opposite_color(im)*num_elements_blobclass+1) then
              print *,"ic invalid in getcolorsum 2"
              print *,"ic=",ic
              print *,"im=",im
              print *,"opposite_color(im)=",opposite_color(im)
              print *,"num_elements_blobclass=",num_elements_blobclass
              stop
             endif

            enddo ! side
           enddo ! dir

          else if (opposite_color(im).eq.0) then
           ! do nothing
          else
           print *,"opposite_color invalid"
           stop
          endif

         enddo ! im=1..nmat

        else
         print *,"base_type invalid"
         stop
        endif
       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo
      enddo
      enddo

      return
      end subroutine FORT_GETCOLORSUM


      subroutine FORT_LEVELRECOLOR( &
        color, &
        DIMS(color), &
        xlo,dx, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        domaincolormap, &
        max_colors_level, &
        level,base_level,arrsize)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
 
      INTEGER_T, intent(in) ::  DIMDEC(color)
      INTEGER_T, intent(in) :: max_colors_level
      INTEGER_T, intent(in) :: level,base_level,arrsize
      INTEGER_T, intent(in) :: domaincolormap(arrsize)
      REAL_T, intent(inout) :: color(DIMV(color))
      INTEGER_T  i, j, k, m, icolor

      call checkbound(fablo,fabhi, &
       DIMS(color), &
       1,-1,4000)

      if (bfact.lt.1) then
       print *,"bfact invalid93"
       stop
      endif
      if (arrsize.ne.2*max_colors_level) then
       print *,"arrsize invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
 
       icolor=NINT(color(D_DECL(i,j,k)))
       if (icolor.ne.0) then
        m=icolor
        if (level.gt.base_level) then
         m=m+max_colors_level
        endif
        color(D_DECL(i,j,k))=domaincolormap(m)
       endif
      enddo
      enddo
      enddo

      return
      end


      subroutine FORT_COLORFILL( &
       mask,DIMS(mask), &
       typefab,DIMS(typefab), &
       color,DIMS(color), &
       ijk,DIMS(ijk), &
       lo,hi, &
       ipass,number_grids,color_per_grid, &
       gridno,max_colors_grid,typedim)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: typedim
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(typefab)
      INTEGER_T, intent(in) :: DIMDEC(color)
      INTEGER_T, intent(in) :: DIMDEC(ijk)
      INTEGER_T, intent(in) :: number_grids,ipass
      INTEGER_T, intent(inout) :: color_per_grid(number_grids)
      INTEGER_T, intent(in) :: gridno
      INTEGER_T, intent(in) :: max_colors_grid

      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: typefab(DIMV(typefab))
      REAL_T, intent(inout) :: color(DIMV(color))
      INTEGER_T, intent(inout) :: ijk(DIMV(ijk),SDIM)
      INTEGER_T i,j,k,icolor,i1,j1,k1
      INTEGER_T istack,ii,jj,kk
      INTEGER_T iprime,jprime,kprime,base_type,test_type
      INTEGER_T ilocal,jlocal,klocal
      INTEGER_T k1lo,k1hi

      call checkbound(lo,hi,DIMS(mask),0,-1,4001)
      call checkbound(lo,hi,DIMS(typefab),0,-1,4001)
      call checkbound(lo,hi,DIMS(color),0,-1,4001)
      call checkbound(lo,hi,DIMS(ijk),0,-1,4001)

      if (gridno.ge.number_grids) then
       print *,"number_grids invalid"
       stop
      endif

      if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else
       print *,"dimension bust"
       stop
      endif

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       if (ipass.eq.1) then
        icolor=NINT(color(D_DECL(i,j,k)))
        if (icolor.gt.0) then
         color(D_DECL(i,j,k))=max_colors_grid*gridno+icolor
        endif
       else if (ipass.eq.0) then
        if (color(D_DECL(i,j,k)).lt.zero) then
         print *,"color has invalid value"
         stop
        endif
        if ((mask(D_DECL(i,j,k)).eq.one).and. &
            (color(D_DECL(i,j,k)).eq.zero)) then
         color_per_grid(gridno+1)=color_per_grid(gridno+1)+1
         color(D_DECL(i,j,k))=color_per_grid(gridno+1)
         base_type=NINT(typefab(D_DECL(i,j,k)))
         istack=0

! recursive stuff begins here ...

         do i1=-1,1
         do j1=-1,1
         do k1=k1lo,k1hi
          if ((i+i1.ge.growlo(1)).and.(i+i1.le.growhi(1)).and. &
              (j+j1.ge.growlo(2)).and.(j+j1.le.growhi(2)).and. &
              (k+k1.ge.growlo(3)).and.(k+k1.le.growhi(3))) then

            ilocal=i+i1
            jlocal=j+j1
            klocal=k+k1

            if (color(D_DECL(ilocal,jlocal,klocal)).eq.zero) then
             color(D_DECL(ilocal,jlocal,klocal))=-one
             if (istack.eq.0) then
              ii=lo(1)
              jj=lo(2)
              kk=lo(SDIM)
              ijk(D_DECL(ii,jj,kk),1)=ilocal
              ijk(D_DECL(ii,jj,kk),2)=jlocal
              if (SDIM.eq.3) then
               ijk(D_DECL(ii,jj,kk),SDIM)=klocal
              endif
              istack=1
             else if (istack.eq.1) then
              ii=ii+1
              if (ii.gt.hi(1)) then
               ii=lo(1)
               jj=jj+1
               if ((jj.gt.hi(2)).and.(SDIM.eq.3)) then
                jj=lo(2)
                kk=kk+1
               endif
              endif
              ijk(D_DECL(ii,jj,kk),1)=ilocal
              ijk(D_DECL(ii,jj,kk),2)=jlocal
              if (SDIM.eq.3) then
               ijk(D_DECL(ii,jj,kk),SDIM)=klocal
              endif
             endif
            endif  ! color=0

          endif ! push neighbors
         enddo
         enddo
         enddo ! i1,j1,k1

         do while (istack.eq.1)

          iprime=ijk(D_DECL(ii,jj,kk),1)
          jprime=ijk(D_DECL(ii,jj,kk),2)
          kprime=0
          if (SDIM.eq.3) then
           kprime=ijk(D_DECL(ii,jj,kk),SDIM)
          endif
          if (color(D_DECL(iprime,jprime,kprime)).lt.zero) then
            color(D_DECL(iprime,jprime,kprime))=zero
          endif
          ii=ii-1
          if (ii.lt.lo(1)) then
            ii=hi(1)
            jj=jj-1
            if (jj.lt.lo(2)) then
             if (SDIM.eq.2) then
              istack=0
             else if (SDIM.eq.3) then 
              jj=hi(2)
              kk=kk-1
              if (kk.lt.lo(SDIM)) then
               istack=0
              endif
             else
              print *,"dimension bust"
              stop
             endif
            endif
          endif

          test_type=NINT(typefab(D_DECL(iprime,jprime,kprime)))
          if ((test_type.eq.base_type).and. &
              (mask(D_DECL(iprime,jprime,kprime)).eq.one).and. &
              (color(D_DECL(iprime,jprime,kprime)).eq.zero)) then
           color(D_DECL(iprime,jprime,kprime))=color(D_DECL(i,j,k))
   
           do i1=-1,1
           do j1=-1,1
           do k1=k1lo,k1hi
            if ((iprime+i1.ge.growlo(1)).and.(iprime+i1.le.growhi(1)).and. &
                (jprime+j1.ge.growlo(2)).and.(jprime+j1.le.growhi(2)).and. &
                (kprime+k1.ge.growlo(3)).and.(kprime+k1.le.growhi(3))) then

              ilocal=iprime+i1
              jlocal=jprime+j1
              klocal=kprime+k1

              if (color(D_DECL(ilocal,jlocal,klocal)).eq.zero) then
               color(D_DECL(ilocal,jlocal,klocal))=-one

               if (istack.eq.0) then
                ii=lo(1)
                jj=lo(2)
                if (SDIM.eq.2) then
                 kk=0
                else if (SDIM.eq.3) then
                 kk=lo(SDIM)
                else
                 print *,"dimension bust"
                 stop
                endif
                ijk(D_DECL(ii,jj,kk),1)=ilocal
                ijk(D_DECL(ii,jj,kk),2)=jlocal
                if (SDIM.eq.3) then
                 ijk(D_DECL(ii,jj,kk),SDIM)=klocal
                endif
                istack=1
               else if (istack.eq.1) then
                ii=ii+1
                if (ii.gt.hi(1)) then
                 ii=lo(1)
                 jj=jj+1
                 if ((jj.gt.hi(2)).and.(SDIM.eq.3)) then
                  jj=lo(2)
                  kk=kk+1
                 endif
                endif
                ijk(D_DECL(ii,jj,kk),1)=ilocal
                ijk(D_DECL(ii,jj,kk),2)=jlocal
                if (SDIM.eq.3) then
                 ijk(D_DECL(ii,jj,kk),SDIM)=klocal
                endif
               endif

              endif ! color=0

            endif  ! test neighbor

           enddo
           enddo
           enddo ! i1,j1,k1

          endif ! push neighbors
 
         enddo ! end while

! recursive stuff ends here

        endif
       else 
        print *,"ipass invalid"
        stop
       endif
      enddo
      enddo
      enddo

      return
      end subroutine FORT_COLORFILL


      subroutine FORT_GRIDRECOLOR( &
       mask, &
       DIMS(mask), &
       color, &
       DIMS(color), &
       xlo,dx, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       levelcolormap,max_colors_grid,number_grids, &
       arrsize)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(color)
      INTEGER_T, intent(in) :: max_colors_grid,number_grids,arrsize
      INTEGER_T, intent(in) :: levelcolormap(arrsize)
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(inout) :: color(DIMV(color))
      INTEGER_T i,j,k,icolor,testsize

      call checkbound(fablo,fabhi, &
       DIMS(mask), &
       1,-1,4000)
      call checkbound(fablo,fabhi, &
       DIMS(color), &
       1,-1,4000)

      if (bfact.lt.1) then
       print *,"bfact invalid94"
       stop
      endif

      testsize=number_grids*max_colors_grid
      if (arrsize.ne.testsize) then
       print *,"arrsize invalid"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       if (mask(D_DECL(i,j,k)).ne.zero) then
        icolor=NINT(color(D_DECL(i,j,k)))
        if ((icolor.le.0).or.(icolor.gt.arrsize)) then
         print *,"icolor invalid in GRIDRECOLOR icolor=",icolor
         print *,"i,j,k ",i,j,k
         print *,"arrsize= ",arrsize
         stop
        endif
        if (levelcolormap(icolor).le.0) then
         print *,"levelcolormap invalid"
         stop
        endif
        color(D_DECL(i,j,k))=levelcolormap(icolor)
       endif
      enddo
      enddo
      enddo
           
      return
      end subroutine FORT_GRIDRECOLOR


! components are: color1,type1,color2,type2,color3,type3
      subroutine FORT_LEVELCOLORINIT( &
        mask,DIMS(mask), &
        color,DIMS(color), &
        xlo,dx, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        level_color, &
        max_colors_level, &
        arrsize,check_corners)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: check_corners
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(color)
      INTEGER_T, intent(in) :: max_colors_level,arrsize
      INTEGER_T, intent(out) :: level_color(arrsize,arrsize)
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: color(DIMV(color),6)
      INTEGER_T i,j,k,icolor,jcolor,testsize
      INTEGER_T k1lo,k1hi
      INTEGER_T ii,jj,kk,base_type,near_type
      INTEGER_T nbase,nbase2
      REAL_T mask2

      call checkbound(fablo,fabhi,DIMS(mask),1,-1,4000)
      call checkbound(fablo,fabhi,DIMS(color),1,-1,4000)

      if (bfact.lt.1) then
       print *,"bfact invalid95"
       stop
      endif

      if ((check_corners.ne.0).and.(check_corners.ne.1)) then
       print *,"check_corners invalid"
       stop
      endif

! max_colors_level=max(colormax[level],colormax[level+1])
! add max_colors_level to level+1 colors (mask=0) in order to differentiate
! from level colors.

      testsize=2*max_colors_level
      if (arrsize.ne.testsize) then
       print *,"arrsize invalid LEVELCOLORINIT"
       stop
      endif

      if (SDIM.eq.3) then
       k1lo=-1
       k1hi=1
      else if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
      else
       print *,"dimension bust"
       stop
      endif

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       do nbase=1,3
        icolor=NINT(color(D_DECL(i,j,k),2*nbase-1))
        base_type=NINT(color(D_DECL(i,j,k),2*nbase))
        if ((icolor.le.0).and.(mask(D_DECL(i,j,k)).ne.zero)) then
         print *,"icolor invalid in LEVELCOLORINIT icolor=",icolor
         print *,"i,j,k ",i,j,k
         stop
        endif
        if (icolor.ne.0) then
         if (mask(D_DECL(i,j,k)).eq.zero) then
          icolor=icolor+max_colors_level
         else if (mask(D_DECL(i,j,k)).ne.one) then
          print *,"mask invalid"
          stop
         endif
         do ii=-1,1
         do jj=-1,1
         do kk=k1lo,k1hi

          if ((check_corners.eq.1).or. &
              (abs(ii)+abs(jj)+abs(kk).le.1)) then

           mask2=mask(D_DECL(i+ii,j+jj,k+kk)) 
           do nbase2=1,3
            jcolor=NINT(color(D_DECL(i+ii,j+jj,k+kk),2*nbase2-1))
            near_type=NINT(color(D_DECL(i+ii,j+jj,k+kk),2*nbase2))
            if ((jcolor.eq.0).and.(mask2.ne.zero)) then
             print *,"jcolor invalid in LEVELCOLORINIT jcolor=",jcolor
             stop
            endif
            if (jcolor.ne.0) then
             if (near_type.eq.base_type) then
              if (mask2.eq.zero) then
               jcolor=jcolor+max_colors_level
              else if (mask2.ne.one) then
               print *,"mask2 invalid"
               stop
              endif
              if ((icolor.gt.arrsize).or.(jcolor.gt.arrsize)) then
               print *,"icolor or jcolor invalid"
               stop
              endif
              if ((mask2.eq.one).or.(mask(D_DECL(i,j,k)).eq.one)) then
               level_color(jcolor,icolor)=1
               level_color(icolor,jcolor)=1
              endif
             endif  !near_type=base_type
            endif ! jcolor.ne.0
           enddo ! nbase2
 
          endif ! check_corners

         enddo
         enddo
         enddo
        endif ! icolor.ne.0
       enddo ! nbase
      enddo
      enddo
      enddo
           
      return
      end subroutine FORT_LEVELCOLORINIT

      subroutine FORT_AVGDOWNCOLOR( &
        problo,dxf, &
        bfact_f,bfact, &
        xlo_fine,dx, &
        crse,DIMS(crse), &
        fine,DIMS(fine), &
        typef,DIMS(typef), &
        typec,DIMS(typec), &
        clo,chi)
      use global_utility_module
      IMPLICIT NONE

      REAL_T problo(SDIM)
      REAL_T dxf(SDIM)
      INTEGER_T bfact_f,bfact
      REAL_T xlo_fine(SDIM)
      REAL_T dx(SDIM)
      INTEGER_T  clo(SDIM),chi(SDIM)
      INTEGER_T  growlo(3),growhi(3)
      INTEGER_T  stenlo(3),stenhi(3)
      INTEGER_T  DIMDEC(crse)
      INTEGER_T  DIMDEC(fine)
      INTEGER_T  DIMDEC(typef)
      INTEGER_T  DIMDEC(typec)
      REAL_T crse(DIMV(crse))
      REAL_T fine(DIMV(fine))
      REAL_T typef(DIMV(typef))
      REAL_T typec(DIMV(typec))
      INTEGER_T icolor
      INTEGER_T  i, j, k, ic, jc, kc
      INTEGER_T  coarse_type,fine_type
      REAL_T wt(SDIM)


      if (bfact_f.lt.1) then
       print *,"bfact_f invalid1 ",bfact_f
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid96"
       stop
      endif
      if ((bfact.ne.bfact_f).and. &
          (bfact.ne.2*bfact_f)) then
       print *,"bfact invalid97"
       stop
      endif


      call checkbound(clo,chi,DIMS(crse),0,-1,410)
      call checkbound(clo,chi,DIMS(typec),0,-1,410)
 
      call growntilebox(clo,chi,clo,chi,growlo,growhi,0) 
      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)
       crse(D_DECL(ic,jc,kc))=zero
       coarse_type=NINT(typec(D_DECL(ic,jc,kc)))
       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi,bfact,bfact_f)
       do i=stenlo(1),stenhi(1)
        call intersect_weight_avg(ic,i,bfact,bfact_f,wt(1))
        if (wt(1).gt.zero) then
         do j=stenlo(2),stenhi(2)
          call intersect_weight_avg(jc,j,bfact,bfact_f,wt(2))
          if (wt(2).gt.zero) then
           do k=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_avg(kc,k,bfact,bfact_f,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             fine_type=NINT(typef(D_DECL(i,j,k)))
             if (coarse_type.eq.fine_type) then
              icolor=NINT(fine(D_DECL(i,j,k)))
              crse(D_DECL(ic,jc,kc))=icolor
             endif
            endif
           enddo ! k
          endif
         enddo ! j
        endif
       enddo ! i
      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine FORT_AVGDOWNCOLOR


! components are: color1,type1,color2,type2,color3,type3
      subroutine FORT_COPYFINECOARSECOLOR( &
        problo,dxf,bfact_f,bfact,xlo_fine,dx, &
        crse,DIMS(crse), &
        fine,DIMS(fine), &
        typef,DIMS(typef), &
        maskf,DIMS(maskf), &
        clo,chi,flo,fhi)
      use global_utility_module
      use probf90_module
      IMPLICIT NONE

      REAL_T problo(SDIM)
      REAL_T dxf(SDIM)
      INTEGER_T bfact_f,bfact
      REAL_T xlo_fine(SDIM)
      REAL_T dx(SDIM)
      INTEGER_T  flo(SDIM),fhi(SDIM)
      INTEGER_T  clo(SDIM),chi(SDIM)
      INTEGER_T  growlo(3),growhi(3)
      INTEGER_T  stenlo(3),stenhi(3)
      INTEGER_T  DIMDEC(crse)
      INTEGER_T  DIMDEC(fine)
      INTEGER_T  DIMDEC(typef)
      INTEGER_T  DIMDEC(maskf)
      REAL_T crse(DIMV(crse),6)
      REAL_T fine(DIMV(fine))
      REAL_T typef(DIMV(typef))
      REAL_T maskf(DIMV(maskf))
      INTEGER_T i, j, k, ic, jc, kc,n
      INTEGER_T fine_type,fine_color
      INTEGER_T icrse,jcrse,alreadyhit,crse_color,crse_type
      INTEGER_T nmat
      REAL_T masktest
      REAL_T wt(SDIM)


      if (bfact_f.lt.1) then
       print *,"bfact_f invalid2  bfact_f=",bfact_f
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid98"
       stop
      endif
      if ((bfact.ne.bfact_f).and. &
          (bfact.ne.2*bfact_f)) then
       print *,"bfact invalid99"
       stop
      endif

      nmat=num_materials

      call checkbound(clo,chi,DIMS(crse),0,-1,410)
      call checkbound(flo,fhi,DIMS(fine),0,-1,410)
      call checkbound(flo,fhi,DIMS(typef),0,-1,410)
      call checkbound(flo,fhi,DIMS(maskf),0,-1,410)
 
      call growntilebox(clo,chi,clo,chi,growlo,growhi,0) 
      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)
       do n=1,6
        crse(D_DECL(ic,jc,kc),n)=zero
       enddo
       icrse=0

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi,bfact,bfact_f)
       do i=stenlo(1),stenhi(1)
        call intersect_weight_avg(ic,i,bfact,bfact_f,wt(1))
        if (wt(1).gt.zero) then
         do j=stenlo(2),stenhi(2)
          call intersect_weight_avg(jc,j,bfact,bfact_f,wt(2))
          if (wt(2).gt.zero) then
           do k=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_avg(kc,k,bfact,bfact_f,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then

             masktest=maskf(D_DECL(i,j,k))
             fine_type=NINT(typef(D_DECL(i,j,k)))
             fine_color=NINT(fine(D_DECL(i,j,k)))

             if (masktest.eq.one) then

              if (fine_color.le.0) then
               print *,"fine_color invalid"
               stop
              endif
              if ((fine_type.le.0).or.(fine_type.gt.nmat)) then
               print *,"fine_type invalid"
               stop
              else
               alreadyhit=0
               do jcrse=1,icrse
                crse_color=NINT(crse(D_DECL(ic,jc,kc),2*jcrse-1))
                crse_type=NINT(crse(D_DECL(ic,jc,kc),2*jcrse))
                if ((crse_type.eq.fine_type).and. &
                    (crse_color.eq.fine_color)) then
                 alreadyhit=1
                endif
               enddo ! jcrse
               if (alreadyhit.eq.0) then
                icrse=icrse+1
                if (icrse.gt.3) then
                 print *,"WARNING type overflow"
                 icrse=3
                else
                 crse(D_DECL(ic,jc,kc),2*icrse-1)=fine_color
                 crse(D_DECL(ic,jc,kc),2*icrse)=fine_type
                endif
               endif
              endif ! nmat>=fine_type>=1

             else if (masktest.ne.zero) then
              print *,"masktest invalid"
              stop
             endif
            endif
           enddo !k
          endif
         enddo !j
        endif
       enddo ! i
      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine FORT_COPYFINECOARSECOLOR



! faceden=1/rho
! facecut=A
! icefacecut=1
! mside=mass
! slope: vof,ref centroid,order,slope,intercept  x nmat
! vof: piecewise constant interp at coarse/fine borders
! mask=1 at interior cells, physical border cells, and fine-fine border cells.
! mask=0 at coarse-fine border cells.
! mask has 3 ghost cells.

      subroutine FORT_INIT_PHYSICS_VARS( &
       tid, &
       FD_curv_interp, &
       curv_min, &
       curv_max, &
       isweep, &
       nrefine_vof, &
       den_interface, &
       visc_interface, &
       heatvisc_interface, &
       speciesvisc_interface, &
       diffusionface_flag, &
       temperatureface_flag, &
       curv_index, &
       pforce_index, &
       faceden_index, &
       facecut_index, &
       icefacecut_index, &
       icemask_index, &
       facevisc_index, &
       faceheat_index, &
       facevel_index, &
       facespecies_index, &
       massface_index, &
       vofface_index, &
       ncphys, &
       latent_heat, &
       freezing_model, &
       distribute_from_target, &
       solidheat_flag, &
       microlayer_size, & ! 1..nmat
       microlayer_substrate, & ! 1..nmat
       microlayer_temperature_substrate, & ! 1..nmat
       spec_material_id, & ! 1..num_species_var+1
       mass_fraction_id, &
       species_evaporation_density, &
       cavitation_vapor_density, &
       cavitation_species, &
       override_density, &
       time, &
       project_option, &
       problo,probhi, &
       visc_coef, &
       nten, &
       xlo,dx, &
       maskcov,DIMS(maskcov), &
       masknbr,DIMS(masknbr), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       curv,DIMS(curv), &
       slope,DIMS(slope), &
       denstate,DIMS(denstate), &
       viscstate,DIMS(viscstate), &
       sol,DIMS(sol), &
        ! voltotal/DeDT_total 1/(rho cv)
       cenDeDT, &
       DIMS(cenDeDT), & 
        ! voltotal/mass_total (1/rho)
       cenden, &
       DIMS(cenden), &   
       cenvof,DIMS(cenvof), &   
       cenvisc,DIMS(cenvisc), &
       vol,DIMS(vol), &
       levelPC,DIMS(levelPC), &
       vofC,DIMS(vofC), &
       vofF,DIMS(vofF), &
       massF,DIMS(massF), &
       modvisc,DIMS(modvisc), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       presbc_arr, &
       velbc, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       num_curv, &
       prescribed_solid_scale, &
       level, &
       finest_level)
      use global_utility_module
      use probf90_module
      use levelset_module
      use godunov_module
      use geometry_intersect_module
      use MOF_routines_module
      use CISL_SANITY_MODULE

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)
      INTEGER_T, intent(in) :: nten
      REAL_T :: curv_min
      REAL_T :: curv_max
      INTEGER_T, intent(in) :: isweep
      INTEGER_T, intent(in) :: nrefine_vof
      REAL_T, intent(in) :: prescribed_solid_scale(nmat)
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: diffusionface_flag
      INTEGER_T, intent(in) :: temperatureface_flag
      INTEGER_T, intent(in) :: curv_index
      INTEGER_T, intent(in) :: pforce_index
      INTEGER_T, intent(in) :: faceden_index
      INTEGER_T, intent(in) :: facecut_index
      INTEGER_T, intent(in) :: icefacecut_index
      INTEGER_T, intent(in) :: icemask_index
      INTEGER_T, intent(in) :: facevisc_index
      INTEGER_T, intent(in) :: faceheat_index
      INTEGER_T, intent(in) :: facevel_index
      INTEGER_T, intent(in) :: facespecies_index
      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys

      INTEGER_T, intent(in) :: solidheat_flag
      REAL_T, intent(in) :: microlayer_size(nmat)
      INTEGER_T, intent(in) :: microlayer_substrate(nmat)
      REAL_T, intent(in) :: microlayer_temperature_substrate(nmat)

      INTEGER_T, intent(in) :: num_curv
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: project_option
      REAL_T, intent(in) :: problo(SDIM),probhi(SDIM)
      REAL_T, intent(in) :: visc_coef
      REAL_T, intent(in) :: latent_heat(2*nten)
      INTEGER_T, intent(in) :: freezing_model(2*nten)
      INTEGER_T, intent(in) :: distribute_from_target(2*nten)
      INTEGER_T :: veldir
      INTEGER_T, intent(in) :: override_density(nmat)
      INTEGER_T, intent(in) :: spec_material_id(num_species_var+1)
      INTEGER_T, intent(in) :: mass_fraction_id(2*nten)
      REAL_T, intent(in) :: species_evaporation_density(num_species_var+1)
      REAL_T, intent(in) :: cavitation_vapor_density(nmat)
      INTEGER_T, intent(in) :: cavitation_species(nmat)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: presbc_arr(SDIM,2)
      INTEGER_T, intent(in) :: velbc(SDIM,2,num_materials_vel*SDIM)
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(masknbr)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(curv)
      INTEGER_T, intent(in) :: DIMDEC(slope)
      INTEGER_T, intent(in) :: DIMDEC(denstate)
      INTEGER_T, intent(in) :: DIMDEC(viscstate)
      INTEGER_T, intent(in) :: DIMDEC(sol)
      INTEGER_T, intent(in) :: DIMDEC(cenDeDT)
      INTEGER_T, intent(in) :: DIMDEC(cenden)
      INTEGER_T, intent(in) :: DIMDEC(cenvof)
      INTEGER_T, intent(in) :: DIMDEC(vol)
      INTEGER_T, intent(in) :: DIMDEC(levelPC)
      INTEGER_T, intent(in) :: DIMDEC(cenvisc)
      INTEGER_T, intent(in) :: DIMDEC(vofC)
      INTEGER_T, intent(in) :: DIMDEC(vofF)
      INTEGER_T, intent(in) :: DIMDEC(massF)
      INTEGER_T, intent(in) :: DIMDEC(modvisc)
      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: masknbr(DIMV(masknbr),4)
      REAL_T, intent(out) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(out) :: yface(DIMV(yface),ncphys)
      REAL_T, intent(out) :: zface(DIMV(zface),ncphys)
      REAL_T, intent(in) :: curv(DIMV(curv),num_curv) 
      REAL_T, intent(in) :: slope(DIMV(slope),nmat*ngeom_recon) 
      REAL_T, intent(in) :: denstate(DIMV(denstate),nmat*num_state_material) 
      REAL_T, intent(in) :: viscstate(DIMV(viscstate),nmat) 
      REAL_T, intent(in) :: sol(DIMV(sol),nparts_def*SDIM) 
      REAL_T, intent(out) :: cenDeDT(DIMV(cenDeDT),nmat+1)
      REAL_T, intent(out) :: cenden(DIMV(cenden),nmat+1)
      REAL_T, intent(out) :: cenvof(DIMV(cenvof),nmat)  
      REAL_T, intent(out) :: cenvisc(DIMV(cenvisc),nmat+1)
      REAL_T, intent(in) :: vol(DIMV(vol))
      REAL_T, intent(in) :: levelPC(DIMV(levelPC),nmat)
      REAL_T, intent(in) :: vofC(DIMV(vofC),nmat)
      REAL_T, intent(in) :: vofF(DIMV(vofF),nrefine_vof)
      REAL_T, intent(in) :: massF(DIMV(massF),nrefine_vof)
      REAL_T, intent(in) :: modvisc(DIMV(modvisc),nmat)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      REAL_T, intent(in) :: den_interface(nten)
      REAL_T, intent(in) :: visc_interface(nten)
      REAL_T, intent(in) :: heatvisc_interface(nten)
      REAL_T, intent(in) :: speciesvisc_interface(nten*num_species_var)

      INTEGER_T im1,jm1,km1
      INTEGER_T i,j,k,ii,jj,kk
      INTEGER_T dir2,inorm,presbclo,presbchi
       ! mask_boundary=1 at left neumann boundary
       ! mask_boundary=2 at right neumann boundary
       ! mask_boundary=0 otherwise
       ! mask_boundary<>0 means the face coefficient (1/rho) should be
       ! equal to 0.0
      INTEGER_T mask_boundary
      INTEGER_T icell,jcell,kcell
      INTEGER_T iside

      INTEGER_T im
      INTEGER_T nmax
      INTEGER_T vofcomp
      REAL_T visc_total,heat_total
      REAL_T spec_total(num_species_var+1) ! +1 to avoid 0 size 
      REAL_T DeDT,DeDT_total

      REAL_T volmat(nmat)
      REAL_T voltotal
      REAL_T voldepart
      REAL_T mass_total
      REAL_T delta_mass
      REAL_T den
      REAL_T density_for_mass_fraction_diffusion
      REAL_T TEMPERATURE

      INTEGER_T LS_consistent
      INTEGER_T LS_consistent_tension

      REAL_T LSplus(nmat)
      REAL_T LSminus(nmat)
      INTEGER_T im_main,im_main_opp,im_opp
      INTEGER_T ireverse
      INTEGER_T iten
      INTEGER_T zeroradius_flag
      INTEGER_T iten_tension,iten_micro
      INTEGER_T im_tension,im_opp_tension
      INTEGER_T im_fluid_micro,im_solid_micro

      REAL_T gradh,gradh_tension,sign_test
      INTEGER_T nten_test
      REAL_T one_over_mu
      REAL_T localvisc(nmat)
      INTEGER_T null_viscosity
       !0=>interior wall 1..nmat=>embedded wall 
       !nmat+1=>domain wall 
      INTEGER_T wall_flag_face 
      REAL_T wallVOF_face
      INTEGER_T imattype
      INTEGER_T dencomp,tempcomp
      REAL_T facevisc_local,faceheat_local
      REAL_T facespecies_local(num_species_var+1)
      REAL_T theta,visc1,visc2,heat1,heat2
      REAL_T spec1(num_species_var+1)
      REAL_T spec2(num_species_var+1)
      REAL_T spec_test
      REAL_T localvisc_plus(nmat)
      REAL_T localvisc_minus(nmat)
      INTEGER_T implus_majority,imminus_majority

      REAL_T local_face(ncphys)
      REAL_T local_volumes(2,nmat)
      REAL_T local_masses(2,nmat)
      INTEGER_T igridlo(3),igridhi(3)

      INTEGER_T dir_refine
      REAL_T, dimension(:,:), allocatable :: comparemassface
      REAL_T, dimension(:,:), allocatable :: comparedenface
      INTEGER_T noslip_wall,velbclo,velbchi
      INTEGER_T imspec,is_zero_visc,is_zero_heat
      INTEGER_T is_zero_spec(num_species_var+1)
      REAL_T solid_velocity

      REAL_T LSIDE(2)
      REAL_T LSIDE_tension(2)
      REAL_T LSIDE_MAT(nmat)
      REAL_T LSIDE_tension_MAT(nmat)

      REAL_T DXMAXLS
      REAL_T FFACE(nmat)
      INTEGER_T irefine
      INTEGER_T solid_present_flag
      REAL_T wtL,wtR,wtsum
      REAL_T xstenMAC(-3:3,SDIM)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T for_hydrostatic_pres

      INTEGER_T icurv,icurv_ofs
      INTEGER_T curv_interp_flag
      INTEGER_T dirL,sideL,im3L,orientL
      INTEGER_T dirR,sideR,im3R,orientR
      REAL_T curvL(5+SDIM)
      REAL_T curvR(5+SDIM)
      REAL_T mu
      INTEGER_T local_maskL,local_maskR,local_masknbr
      INTEGER_T covered_face,coarse_fine_face
      INTEGER_T prescribed_flag
      INTEGER_T FD_curv_interp
      REAL_T mofdata(nmat*ngeom_recon)
      INTEGER_T micro_table(nmat,nmat)
      REAL_T predict_face_afrac_solid
      REAL_T predict_face_afrac_prescribed
      INTEGER_T is_solid_face
      INTEGER_T is_prescribed_face
      INTEGER_T im_solid
      INTEGER_T im_prescribed
      INTEGER_T im_solid_valid
      INTEGER_T im_prescribed_valid
      INTEGER_T partid_solid
      INTEGER_T partid_prescribed
      REAL_T total_density,evap_den,local_mfrac
      INTEGER_T ispecies,im_species,speccomp
      INTEGER_T im_prescribed_primary

! INIT_PHYSICS_VARS code starts here:

      for_hydrostatic_pres=0

      nhalf=3

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if (nrefine_vof.ne.2*nmat*SDIM) then
       print *,"nrefine_vof invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level or finest_level invalid init physics vars"
       stop
      endif
      if ((FD_curv_interp.ne.0).and. &
          (FD_curv_interp.ne.1)) then
       print *,"FD_curv_interp invalid"
       stop
      endif

      ! height function curvature
      ! finite difference curvature
      ! pforce
      ! marangoni force
      ! dir * side dir=1..sdim side=-1 or 1
      ! im3
      ! x nten
      if (num_curv.ne.nten*(5+SDIM)) then
       print *,"num_curv invalid"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_INIT_PHYSICS_VARS"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_INIT_PHYSICS_VARS"
       stop
      endif

      dir_refine=SDIM

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if ((diffusionface_flag.ne.0).and. &
          (diffusionface_flag.ne.1)) then
       print *,"diffusionface_flag invalid"
       stop
      endif
      if ((temperatureface_flag.ne.0).and. &
          (temperatureface_flag.ne.1)) then
       print *,"temperatureface_flag invalid"
       stop
      endif

      if ((solidheat_flag.lt.0).or. &
          (solidheat_flag.gt.2)) then
       print *,"solidheat_flag invalid"
       stop
      endif

      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid phys vars nten nten_test ",nten,nten_test
       stop
      endif

      if ((project_option.eq.0).or. &
          (project_option.eq.1).or. &
          (project_option.eq.10)) then
       ! do nothing
      else
       print *,"project_option invalid PHYSVARS"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid init phys vars "
       print *,"ngeom_recon= ",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid init phys vars "
       print *,"ngeom_raw= ",ngeom_raw
       stop
      endif
      if (curv_index.ne.0) then
       print *,"curv_index invalid"
       stop
      endif
      if ((pforce_index.ne.1).or. &
          (faceden_index.ne.2).or. &
          (facecut_index.ne.3).or. &
          (icefacecut_index.ne.4).or. &
          (icemask_index.ne.5).or. &
          (facevisc_index.ne.6).or. &
          (faceheat_index.ne.7).or. &
          (facevel_index.ne.8).or. &
          (facespecies_index.ne.9).or. &
          (massface_index.ne.facespecies_index+num_species_var).or. &
          (vofface_index.ne.massface_index+2*nmat)) then
       print *,"face_index bust 4"
       stop
      endif
      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: INIT_PHYSICS_VARS

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,213)
      call checkbound(fablo,fabhi,DIMS(masknbr),1,-1,213)

      call checkbound(fablo,fabhi,DIMS(xface),0,0,213)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,214)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,215)

      call checkbound(fablo,fabhi, &
       DIMS(cenDeDT), &
       1,-1,216)
      call checkbound(fablo,fabhi,DIMS(cenden),1,-1,217)
      call checkbound(fablo,fabhi,DIMS(cenvof),1,-1,217)
      call checkbound(fablo,fabhi, &
       DIMS(cenvisc), &
       1,-1,218)

      call checkbound(fablo,fabhi,DIMS(slope),1,-1,219)
      call checkbound(fablo,fabhi,DIMS(curv),1,-1,221)
      call checkbound(fablo,fabhi, &
       DIMS(denstate), &
       1,-1,223)
      call checkbound(fablo,fabhi, &
       DIMS(viscstate), &
       1,-1,224)
      call checkbound(fablo,fabhi,DIMS(sol),1,-1,225)
      call checkbound(fablo,fabhi,DIMS(vol),1,-1,227)
      call checkbound(fablo,fabhi, &
       DIMS(levelPC), &
       2,-1,229) 
      call checkbound(fablo,fabhi,DIMS(vofC),1,-1,227)
      call checkbound(fablo,fabhi,DIMS(vofF),1,-1,227)
      call checkbound(fablo,fabhi,DIMS(massF),1,-1,227)
      call checkbound(fablo,fabhi,DIMS(modvisc),1,-1,227)

      do im=1,nten

       if ((freezing_model(im).lt.0).or. &
           (freezing_model(im).gt.5)) then
        print *,"freezing_model invalid"
        stop
       endif
       if ((freezing_model(im+nten).lt.0).or. &
           (freezing_model(im+nten).gt.5)) then
        print *,"freezing_model invalid"
        stop
       endif

       if ((distribute_from_target(im).lt.0).or. &
           (distribute_from_target(im).gt.1)) then
        print *,"distribute_from_target invalid"
        stop
       endif
       if ((distribute_from_target(im+nten).lt.0).or. &
           (distribute_from_target(im+nten).gt.1)) then
        print *,"distribute_from_target invalid"
        stop
       endif

       if ((den_interface(im).lt.zero).or. &
           (visc_interface(im).lt.zero).or. &
           (heatvisc_interface(im).lt.zero)) then
        print *,"den,visc, or heat interface coeff wrong"
        stop
       endif
       do imspec=1,num_species_var
        if (speciesvisc_interface(nten*(imspec-1)+im).lt.zero) then
         print *,"species interface coeff wrong"
         stop
        endif
       enddo

      enddo ! im=1..nten

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).gt.0).and. &
                (fort_material_type(im).le.fort_max_num_eos)) then
        ! do nothing
       else 
        print *,"material_type invalid"
        stop
       endif

       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"density must be positive init_physics_vars"
        print *,"im,denconst ",im,fort_denconst(im)
        stop
       endif
       if (fort_tempconst(im).gt.zero) then
        ! do nothing
       else
        print *,"INIT_PHYSICS_VAR:temperature must be positive"
        print *,"im,fort_tempconst : ",im,fort_tempconst(im)
        stop
       endif
       if (fort_energyconst(im).gt.zero) then
        ! do nothing
       else
        print *,"energy must be positive in FORT_INIT_PHYSICS_VARS"
        print *,"im= ",im
        print *,"fort_energyconst(im)= ",fort_energyconst(im)
        stop
       endif
       mu=get_user_viscconst(im,fort_denconst(im),fort_tempconst(im))
       if (mu.ge.zero) then
        ! do nothing
       else
        print *,"viscosity cannot be negative"
        stop
       endif

       if ((override_density(im).ne.0).and. &
           (override_density(im).ne.1).and. &
           (override_density(im).ne.2)) then
        print *,"override_density invalid"
        stop
       endif

      enddo ! im=1..nmat  (checking parameters)

       ! create look-up table for thin liquid layer model 4
      do im_fluid_micro=1,nmat
      do im_solid_micro=1,nmat
       micro_table(im_fluid_micro,im_solid_micro)=0
      enddo 
      enddo 

      if (solidheat_flag.eq.0) then ! diffuse in solid
       do im=1,nmat-1
        if (is_rigid(nmat,im).eq.0) then
         do im_opp=im+1,nmat
          if (is_rigid(nmat,im_opp).eq.0) then
           call get_iten(im,im_opp,iten,nmat)
           do ireverse=0,1
            if (latent_heat(iten+ireverse*nten).ne.zero) then
             if (freezing_model(iten+ireverse*nten).eq.0) then
              im_solid_micro=microlayer_substrate(im)
              if ((im_solid_micro.ge.1).and. &
                  (im_solid_micro.le.nmat)) then
               if (is_rigid(nmat,im_solid_micro).ne.1) then
                print *,"is_rigid(nmat,im_solid_micro) invalid"
                stop
               endif 
               micro_table(im_opp,im_solid_micro)=iten
              else if (im_solid_micro.eq.0) then
               ! do nothing
              else
               print *,"im_solid_micro invalid"
               stop
              endif 
              im_solid_micro=microlayer_substrate(im_opp)
              if ((im_solid_micro.ge.1).and. &
                  (im_solid_micro.le.nmat)) then
               if (is_rigid(nmat,im_solid_micro).ne.1) then
                print *,"is_rigid(nmat,im_solid_micro) invalid"
                stop
               endif 
               micro_table(im,im_solid_micro)=iten
              else if (im_solid_micro.eq.0) then
               ! do nothing
              else
               print *,"im_solid_micro invalid"
               stop
              endif 
             endif
            endif
           enddo ! ireverse=0..1
          else if (is_rigid(nmat,im_opp).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid"
           stop
          endif
         enddo !im_opp=im+1,nmat
        else if (is_rigid(nmat,im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(nmat,im) invalid"
         stop
        endif
       enddo ! do im=1,nmat-1
      else if ((solidheat_flag.eq.1).or. & ! dirichlet
               (solidheat_flag.eq.2)) then ! neumann
       ! do nothing
      else
       print *,"solidheat_flag invalid"
       stop
      endif

       ! face centered quantities

      if (isweep.eq.0) then

       do veldir=0,SDIM-1

        ii=0
        jj=0
        kk=0
        if (veldir.eq.0) then
         ii=1
        else if (veldir.eq.1) then
         jj=1
        else if ((veldir.eq.2).and.(SDIM.eq.3)) then
         kk=1
        else
         print *,"veldir invalid"
         stop
        endif

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi, &
         igridlo,igridhi,0,veldir) 

          ! first init xface,yface,zface (vofface_index+1,...)
          ! then follow with the rest...
        do i=igridlo(1),igridhi(1)
        do j=igridlo(2),igridhi(2)
        do k=igridlo(3),igridhi(3)

         local_maskL=NINT(maskcov(D_DECL(i-ii,j-jj,k-kk)))
         local_maskR=NINT(maskcov(D_DECL(i,j,k)))
         if ((local_maskL.eq.0).and.(local_maskR.eq.0)) then
          covered_face=2
         else if ((local_maskL.eq.1).and.(local_maskR.eq.1)) then
          covered_face=0
         else if ((local_maskL.eq.1).and.(local_maskR.eq.0)) then
          covered_face=1
         else if ((local_maskL.eq.0).and.(local_maskR.eq.1)) then
          covered_face=1
         else
          print *,"local_maskL or local_maskR invalid"
          stop
         endif

          ! in: FORT_INIT_PHYSICS_VARS 
         local_face(icemask_index+1)=one
         local_face(curv_index+1)=zero
         local_face(pforce_index+1)=zero
         local_face(facevel_index+1)=zero

         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,veldir+1)

         wtL=(xstenMAC(0,veldir+1)-xstenMAC(-1,veldir+1))
         wtR=(xstenMAC(1,veldir+1)-xstenMAC(0,veldir+1))
         wtsum=wtL+wtR
         if ((wtL.le.zero).or.(wtR.le.zero).or.(wtsum.le.zero)) then
          print *,"wtL, wtR, or wtsum invalid" 
          stop
         endif
         wtL=wtL/wtsum
         wtR=wtR/wtsum

         do iside=0,1
          if (iside.eq.0) then
           icell=i-ii
           jcell=j-jj
           kcell=k-kk
          else if (iside.eq.1) then
           icell=i
           jcell=j
           kcell=k
          else
           print *,"iside invalid"
           stop
          endif

          do im=1,nmat
           irefine=veldir*2*nmat+(1-iside)*nmat+im
           local_volumes(iside+1,im)=vofF(D_DECL(icell,jcell,kcell),irefine)
           local_masses(iside+1,im)=massF(D_DECL(icell,jcell,kcell),irefine)
          enddo ! im=1..nmat
         enddo ! iside=0..1

         do iside=0,1
         do im=1,nmat
          local_face(vofface_index+2*(im-1)+iside+1)= &
            local_volumes(iside+1,im)
          local_face(massface_index+2*(im-1)+iside+1)= &
            local_masses(iside+1,im)
         enddo ! im=1..nmat
         enddo ! iside=0..1

         if (veldir.eq.0) then
          inorm=i
         else if (veldir.eq.1) then
          inorm=j
         else if ((veldir.eq.2).and.(SDIM.eq.3)) then
          inorm=k
         else
          print *,"veldir invalid physics_var"
          stop
         endif

         if (levelrz.eq.0) then
          im1=i-ii
         else if (levelrz.eq.1) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          if (xstenMAC(-1,1).le.VOFTOL*dx(1)) then
           im1=i
          else
           im1=i-ii
          endif
         else if (levelrz.eq.3) then
          im1=i-ii
         else
          print *,"levelrz invalid init physics var"
          stop
         endif

         jm1=j-jj

         if (SDIM.eq.3) then
          km1=k-kk
         else if (SDIM.eq.2) then
          km1=0
         else
          print *,"dimension bust"
          stop
         endif

         do im=1,nmat
          LSplus(im)=levelPC(D_DECL(i,j,k),im)
          LSminus(im)=levelPC(D_DECL(im1,jm1,km1),im)
         enddo
          ! checks rigid and non-rigid materials.
         call get_primary_material(LSplus,nmat,implus_majority)
         call get_primary_material(LSminus,nmat,imminus_majority)

         if ((implus_majority.lt.1).or.(implus_majority.gt.nmat).or. &
             (imminus_majority.lt.1).or.(imminus_majority.gt.nmat)) then
          print *,"implus_majority or imminus_majority invalid"
          stop
         endif

         if ((is_rigid(nmat,implus_majority).eq.1).or. &
             (is_rigid(nmat,imminus_majority).eq.1)) then
          predict_face_afrac_solid=zero
         else if ((is_rigid(nmat,implus_majority).eq.0).and. &
                  (is_rigid(nmat,imminus_majority).eq.0)) then
          predict_face_afrac_solid=one
         else
          print *,"implus_majority or imminus_majority invalid"
          stop
         endif

         if ((is_prescribed(nmat,implus_majority).eq.1).or. &
             (is_prescribed(nmat,imminus_majority).eq.1)) then
          predict_face_afrac_prescribed=zero
         else if ((is_prescribed(nmat,implus_majority).eq.0).and. &
                  (is_prescribed(nmat,imminus_majority).eq.0)) then
          predict_face_afrac_prescribed=one
         else
          print *,"implus_majority or imminus_majority invalid"
          stop
         endif

         call fixed_face( &
          nmat, &
          predict_face_afrac_solid, &
          predict_face_afrac_prescribed, &
          LSminus,LSplus, &
          is_solid_face, &
          is_prescribed_face, &
          im_solid, &
          im_prescribed, &
          im_solid_valid, &
          im_prescribed_valid, &
          partid_solid, &
          partid_prescribed) 

         if (is_prescribed_face.eq.1) then
          if (is_solid_face.eq.1) then
           if (im_prescribed_valid.eq.1) then
            if (im_solid_map(partid_prescribed+1)+1.eq.im_prescribed) then
             dir2=partid_prescribed*SDIM+veldir+1
             solid_velocity=wtR*sol(D_DECL(i,j,k),dir2)+ &
                            wtL*sol(D_DECL(im1,jm1,km1),dir2)
            else
             print *,"im_solid_map(partid_prescribed+1)+1.ne.im_prescribed"
             stop
            endif
           else 
            print *,"im_prescribed_valid bust"
            stop
           endif 
          else
           print *,"is_solid_face invalid 2 ",is_solid_face
           stop
          endif
         else if (is_prescribed_face.eq.0) then
          solid_velocity=zero
         else
          print *,"is_prescribed_face invalid"
          stop
         endif

         presbclo=presbc_arr(veldir+1,1)
         presbchi=presbc_arr(veldir+1,2)
         velbclo=velbc(veldir+1,1,veldir+1)
         velbchi=velbc(veldir+1,2,veldir+1)

          ! mask_boundary=1 at left neumann boundary
          ! mask_boundary=2 at right neumann boundary
          ! mask_boundary=0 otherwise
          ! mask_boundary<>0 means the face coefficient (1/rho) should be
          ! equal to 0.0
         mask_boundary=0

         noslip_wall=0
         coarse_fine_face=0

         if (inorm.eq.fablo(veldir+1)) then

          local_masknbr=NINT(masknbr(D_DECL(i-ii,j-jj,k-kk),1))
          if ((local_masknbr.ne.0).and.(local_masknbr.ne.1)) then
           print *,"local_masknbr invalid"
           stop
          endif
          if ((presbclo.eq.INT_DIR).and.(local_masknbr.eq.0)) then
           coarse_fine_face=1
          endif

          if ((presbclo.eq.INT_DIR).or.(presbclo.eq.EXT_DIR)) then
           ! do nothing
          else if ((presbclo.eq.REFLECT_EVEN).or. &
                   (presbclo.eq.FOEXTRAP)) then
           mask_boundary=1
          else
           print *,"presbclo invalid"
           stop
          endif

          if (velbclo.eq.REFLECT_ODD) then
           local_face(facevel_index+1)=zero
          else if (velbclo.eq.EXT_DIR) then
           iside=1
           call velbc_override(time,veldir+1,iside,veldir+1, &
            local_face(facevel_index+1), &
            xstenMAC,nhalf,dx,bfact)
           noslip_wall=1
          else if ((velbclo.eq.INT_DIR).or. &
                   (velbclo.eq.REFLECT_EVEN).or. &
                   (velbclo.eq.FOEXTRAP)) then
           ! do nothing
          else
           print *,"velbclo invalid"
           stop
          endif

         else if (inorm.eq.fabhi(veldir+1)+1) then

          local_masknbr=NINT(masknbr(D_DECL(i,j,k),1))
          if ((local_masknbr.ne.0).and.(local_masknbr.ne.1)) then
           print *,"local_masknbr invalid"
           stop
          endif
          if ((presbchi.eq.INT_DIR).and.(local_masknbr.eq.0)) then
           coarse_fine_face=1
          endif

          if ((presbchi.eq.INT_DIR).or.(presbchi.eq.EXT_DIR)) then
           ! do nothing
          else if ((presbchi.eq.REFLECT_EVEN).or. &
                   (presbchi.eq.FOEXTRAP)) then
           mask_boundary=2
          else
           print *,"presbclo invalid"
           stop
          endif

          if (velbchi.eq.REFLECT_ODD) then
           local_face(facevel_index+1)=zero
          else if (velbchi.eq.EXT_DIR) then
           iside=2
           call velbc_override(time,veldir+1,iside,veldir+1, &
            local_face(facevel_index+1), &
            xstenMAC,nhalf,dx,bfact)
           noslip_wall=1
          else if ((velbchi.eq.INT_DIR).or. &
                   (velbchi.eq.REFLECT_EVEN).or. &
                   (velbchi.eq.FOEXTRAP)) then
           ! do nothing
          else
           print *,"velbchi invalid"
           stop
          endif

         else if ((inorm.gt.fablo(veldir+1)).and. &
                  (inorm.lt.fabhi(veldir+1)+1)) then
          ! do nothing
         else
          print *,"inorm invalid"
          stop
         endif

         if (is_prescribed_face.eq.1) then
          if (is_solid_face.eq.1) then
           local_face(facevel_index+1)=solid_velocity
          else
           print *,"is_solid_face invalid 3 ",is_solid_face
           stop
          endif
         else if (is_prescribed_face.eq.0) then
          ! do nothing
         else
          print *,"is_prescribed_face invalid"
          stop
         endif

          ! does not consider materials with FSI_flag=1,2,4
          ! gradh>0 => im_main material dominates right and im_main_opp left.
          ! gradh<0 => im_main material dominates left and im_main_opp right.
          ! im_main<im_main_opp
         if (is_solid_face.eq.1) then
          gradh=zero
          gradh_tension=zero
         else if (is_solid_face.eq.0) then
          call fluid_interface(LSminus,LSplus,gradh,im_main_opp,im_main,nmat)
          call get_dxmaxLS(dx,bfact,DXMAXLS)

          if (covered_face.eq.2) then
           gradh_tension=zero
          else if ((covered_face.eq.0).or.(covered_face.eq.1)) then
           call fluid_interface_tension(LSminus,LSplus,gradh_tension, &
            im_opp_tension,im_tension,nmat,nten)
          else
           print *,"covered_face invalid"
           stop
          endif
         else
          print *,"is_solid_face invalid 4 ",is_solid_face
          stop
         endif

         facevisc_local=zero
         faceheat_local=zero
         do imspec=1,num_species_var
          facespecies_local(imspec)=zero
         enddo
         do im=1,nmat
          localvisc_plus(im)=modvisc(D_DECL(i,j,k),im)
          localvisc_minus(im)=modvisc(D_DECL(im1,jm1,km1),im)
         enddo

         if (gradh.ne.zero) then
          if (im_main.ge.im_main_opp) then
           print *,"fluid_interface bust"
           stop
          endif
          call get_iten(im_main,im_main_opp,iten,nmat)

          do im=1,nmat
           LSIDE_MAT(im)=levelPC(D_DECL(im1,jm1,km1),im)
          enddo
          call get_LS_extend(LSIDE_MAT,nmat,iten,LSIDE(1))

          do im=1,nmat
           LSIDE_MAT(im)=levelPC(D_DECL(i,j,k),im)
          enddo
          call get_LS_extend(LSIDE_MAT,nmat,iten,LSIDE(2))

          if (LSIDE(1)*LSIDE(2).le.zero) then
           LS_consistent=1
          else
           LS_consistent=0
          endif
         else if (gradh.eq.zero) then
          LS_consistent=0
          LSIDE(1)=zero
          LSIDE(2)=zero
         else
          print *,"gradh bust"
          stop
         endif

         if (gradh_tension.ne.zero) then
          if (im_tension.ge.im_opp_tension) then
           print *,"fluid_interface_tension bust"
           stop
          endif
          call get_iten(im_tension,im_opp_tension,iten_tension,nmat)

          do im=1,nmat
           LSIDE_tension_MAT(im)=levelPC(D_DECL(im1,jm1,km1),im)
          enddo
          call get_LS_extend(LSIDE_tension_MAT,nmat,iten_tension, &
           LSIDE_tension(1))

          do im=1,nmat
           LSIDE_tension_MAT(im)=levelPC(D_DECL(i,j,k),im)
          enddo
          call get_LS_extend(LSIDE_tension_MAT,nmat,iten_tension, &
           LSIDE_tension(2))

          sign_test=gradh_tension*(LSIDE_tension(2)-LSIDE_tension(1))
          if (sign_test.ge.zero) then
           LS_consistent_tension=1
          else
           LS_consistent_tension=0
          endif
         else if (gradh_tension.eq.zero) then
          LS_consistent_tension=0
          LSIDE_tension(1)=zero
          LSIDE_tension(2)=zero
         else
          print *,"gradh_tension bust"
          stop
         endif

         solid_present_flag=0
         if ((is_rigid(nmat,implus_majority).eq.1).and. &
             (is_rigid(nmat,imminus_majority).eq.1)) then
          solid_present_flag=1
         else if ((is_rigid(nmat,implus_majority).eq.1).and. &
                  (is_rigid(nmat,imminus_majority).eq.0)) then
          solid_present_flag=2
         else if ((is_rigid(nmat,implus_majority).eq.0).and. &
                  (is_rigid(nmat,imminus_majority).eq.1)) then
          solid_present_flag=3
         else if ((is_rigid(nmat,implus_majority).eq.0).and. &
                  (is_rigid(nmat,imminus_majority).eq.0)) then
          solid_present_flag=0
         else
          print *,"implus_majority or imminus_majority invalid"
          stop
         endif

          ! both adjoining cells are solid cells.
         if (solid_present_flag.eq.1) then

          facevisc_local= &
           half*(localvisc_plus(implus_majority)+ &
                 localvisc_minus(imminus_majority))

          facevisc_local=zero  ! dirichlet cond for velocity at the solid.

          faceheat_local= &
           half*(get_user_heatviscconst(implus_majority)+ &
                 get_user_heatviscconst(imminus_majority))

          do imspec=1,num_species_var
           facespecies_local(imspec)= &
             half*(fort_speciesviscconst((imspec-1)*nmat+implus_majority)+ &
                   fort_speciesviscconst((imspec-1)*nmat+imminus_majority))
          enddo

           ! both adjoining cells are solid cells.
          if (solidheat_flag.eq.0) then ! diffuse in solid
           ! do nothing (faceheat_local already defined)
          else if (solidheat_flag.eq.1) then ! dirichlet
           faceheat_local=zero
          else if (solidheat_flag.eq.2) then ! neumann
           faceheat_local=zero
          else
           print *,"solidheat_flag invalid"
           stop
          endif

          ! one adjoining cell is a solid cell.
         else if ((solid_present_flag.eq.2).or. &
                  (solid_present_flag.eq.3)) then
         
          facevisc_local=wtR*localvisc_plus(implus_majority)+ &
                         wtL*localvisc_minus(imminus_majority)
          faceheat_local=wtR*get_user_heatviscconst(implus_majority)+ &
                         wtL*get_user_heatviscconst(imminus_majority)
          do imspec=1,num_species_var
           facespecies_local(imspec)= &
            wtR*fort_speciesviscconst((imspec-1)*nmat+implus_majority)+ &
            wtL*fort_speciesviscconst((imspec-1)*nmat+imminus_majority)
          enddo

          iten_micro=0
          if (solidheat_flag.eq.0) then ! diffuse in solid
           im_solid_micro=0
           im_fluid_micro=0
           if (is_rigid(nmat,implus_majority).eq.1) then
            im_solid_micro=implus_majority
           else if (is_rigid(nmat,implus_majority).eq.0) then
            im_fluid_micro=implus_majority
           else
            print *,"is_rigid(nmat,implus_majority) invalid"
            stop
           endif
           if (is_rigid(nmat,imminus_majority).eq.1) then
            im_solid_micro=imminus_majority
           else if (is_rigid(nmat,imminus_majority).eq.0) then
            im_fluid_micro=imminus_majority
           else
            print *,"is_rigid(nmat,imminus_majority) invalid"
            stop
           endif
           if ((im_fluid_micro.ge.1).and. &
               (im_fluid_micro.le.nmat).and. &
               (im_solid_micro.ge.1).and. &
               (im_solid_micro.le.nmat)) then
            iten_micro=micro_table(im_fluid_micro,im_solid_micro)
           else if ((im_fluid_micro.eq.0).or. &
                    (im_solid_micro.eq.0)) then
            ! do nothing
           else
            print *,"im_fluid_micro or im_solid_micro invalid"
            stop
           endif
          else if ((solidheat_flag.eq.1).or. & ! dirichlet
                   (solidheat_flag.eq.2)) then ! neumann
           ! do nothing
          else
           print *,"solidheat_flag invalid"
           stop
          endif
 
          if (solidheat_flag.eq.0) then
           ! temperature diffusion in the solid
           if (iten_micro.eq.0) then
            ! do nothing
           else if ((iten_micro.ge.1).and.(iten_micro.le.nten)) then
            faceheat_local=zero ! internal TSAT dirichlet bc 
           else
            print *,"iten_micro invalid"
            stop
           endif  
          else if (solidheat_flag.eq.1) then
           ! dirichlet temperature bc
           call microcell_heat_model(faceheat_local,dx,veldir)
          else if (solidheat_flag.eq.2) then
           faceheat_local=zero  ! neumann temperature bc
          else
           print *,"solidheat_flag invalid"
           stop
          endif

           ! check if face coefficient needs to be replaced
           ! prescribed interface coefficient.
          if (implus_majority.ne.imminus_majority) then
           call get_iten(imminus_majority,implus_majority,iten,nmat)

           if (heatvisc_interface(iten).eq.zero) then
            ! do nothing
           else if (heatvisc_interface(iten).gt.zero) then

            if (latent_heat(iten).ne.zero) then
             if ((freezing_model(iten).eq.0).or. &
                 (freezing_model(iten).eq.5)) then
              print *,"heatvisc_interface invalid"
              stop
             endif 
            endif 
            if (latent_heat(iten+nten).ne.zero) then
             if ((freezing_model(iten+nten).eq.0).or. &
                 (freezing_model(iten+nten).eq.5)) then
              print *,"heatvisc_interface invalid"
              stop
             endif 
            endif 

            faceheat_local=heatvisc_interface(iten)
           else
            print *,"heatvisc_interface invalid"
            stop
           endif
          endif ! implus <> imminus

         else if (solid_present_flag.eq.0) then

          if (gradh.eq.zero) then

           facevisc_local=wtR*localvisc_plus(implus_majority)+ &
                          wtL*localvisc_minus(imminus_majority)
           faceheat_local=wtR*get_user_heatviscconst(implus_majority)+ &
                          wtL*get_user_heatviscconst(imminus_majority)
           do imspec=1,num_species_var
            facespecies_local(imspec)= &
             wtL*fort_speciesviscconst((imspec-1)*nmat+imminus_majority)+ &
             wtR*fort_speciesviscconst((imspec-1)*nmat+implus_majority)
           enddo

          else if (gradh.ne.zero) then

           if ((im_main.gt.nmat).or.(im_main_opp.gt.nmat)) then
            print *,"im_main or im_main_opp bust 3"
            stop
           endif
           if (im_main.ge.im_main_opp) then
            print *,"im_main or im_main_opp invalid"
            stop
           endif
           call get_iten(im_main,im_main_opp,iten,nmat)

           if (LSIDE(1).ge.LSIDE(2)) then
            visc1=localvisc_minus(im_main)
            visc2=localvisc_plus(im_main_opp)
           else if (LSIDE(1).le.LSIDE(2)) then
            visc1=localvisc_plus(im_main)
            visc2=localvisc_minus(im_main_opp)
           else
            print *,"LSIDE bust"
            stop
           endif

           heat1=get_user_heatviscconst(im_main) 
           heat2=get_user_heatviscconst(im_main_opp) 
           if ((heat1.lt.zero).or.(heat2.lt.zero)) then
            print *,"heat1 or heat2 invalid"
            stop
           endif

           do imspec=1,num_species_var
            spec1(imspec)=fort_speciesviscconst((imspec-1)*nmat+im_main)
            spec2(imspec)=fort_speciesviscconst((imspec-1)*nmat+im_main_opp)
           enddo
  
             ! 1/s = (wL/sL + wR/sR)/(wL+wR)=(wL sR + wR sL)/(sL sR)*1/(wL+wR)
           if ((visc1.le.zero).or.(visc2.le.zero)) then
            facevisc_local=zero
           else if ((LSIDE(1).eq.zero).and.(LSIDE(2).eq.zero)) then
            facevisc_local=two*visc1*visc2/(visc1+visc2)
           else if ((LSIDE(1).ge.zero).and.(LSIDE(2).ge.zero))  then
            facevisc_local=visc1
           else if ((LSIDE(1).le.zero).and.(LSIDE(2).le.zero)) then
            facevisc_local=visc2
           else if (LSIDE(2).gt.LSIDE(1)) then
            theta=LSIDE(2)/(LSIDE(2)-LSIDE(1))
            facevisc_local=theta/visc1+(one-theta)/visc2
            facevisc_local=one/facevisc_local
           else if (LSIDE(1).gt.LSIDE(2)) then
            theta=LSIDE(1)/(LSIDE(1)-LSIDE(2))
            facevisc_local=theta/visc1+(one-theta)/visc2
            facevisc_local=one/facevisc_local
           else
            print *,"LSIDE bust"
            stop
           endif

           if (visc_interface(iten).eq.zero) then
            ! do nothing
           else if (visc_interface(iten).gt.zero) then
            facevisc_local=visc_interface(iten)
           else
            print *,"visc_interface invalid"
            stop
           endif

            ! STEFANSOLVER will set the thermal face coefficient
            ! to zero where appropriate.
           if ((heat1.eq.zero).or.(heat2.eq.zero)) then
            faceheat_local=zero
           else if ((LSIDE(1).eq.zero).and.(LSIDE(2).eq.zero)) then
            faceheat_local=two*heat1*heat2/(heat1+heat2)
           else if ((LSIDE(1).ge.zero).and.(LSIDE(2).ge.zero))  then
            faceheat_local=heat1
           else if ((LSIDE(1).le.zero).and.(LSIDE(2).le.zero)) then
            faceheat_local=heat2
           else if (LSIDE(2).gt.LSIDE(1)) then
            theta=LSIDE(2)/(LSIDE(2)-LSIDE(1))
            faceheat_local=theta/heat1+(one-theta)/heat2
            faceheat_local=one/faceheat_local
           else if (LSIDE(1).gt.LSIDE(2)) then
            theta=LSIDE(1)/(LSIDE(1)-LSIDE(2))
            faceheat_local=theta/heat1+(one-theta)/heat2
            faceheat_local=one/faceheat_local
           else
            print *,"LSIDE bust"
            stop
           endif

           if (heatvisc_interface(iten).eq.zero) then
            ! do nothing
           else if (heatvisc_interface(iten).gt.zero) then

            if (latent_heat(iten).ne.zero) then
             if ((freezing_model(iten).eq.0).or. &
                 (freezing_model(iten).eq.5)) then
              print *,"heatvisc_interface invalid"
              stop
             endif 
            endif 
            if (latent_heat(iten+nten).ne.zero) then
             if ((freezing_model(iten+nten).eq.0).or. &
                 (freezing_model(iten+nten).eq.5)) then
              print *,"heatvisc_interface invalid"
              stop
             endif 
            endif 

            faceheat_local=heatvisc_interface(iten)
           else
            print *,"heatvisc_interface invalid"
            stop
           endif

           do imspec=1,num_species_var
            if ((spec1(imspec).le.zero).or.(spec2(imspec).le.zero)) then
             facespecies_local(imspec)=zero
            else if ((LSIDE(1).eq.zero).and.(LSIDE(2).eq.zero)) then
             facespecies_local(imspec)= &
               two*spec1(imspec)*spec2(imspec)/(spec1(imspec)+spec2(imspec))
            else if ((LSIDE(1).ge.zero).and.(LSIDE(2).ge.zero))  then
             facespecies_local(imspec)=spec1(imspec)
            else if ((LSIDE(1).le.zero).and.(LSIDE(2).le.zero)) then
             facespecies_local(imspec)=spec2(imspec)
            else if (LSIDE(2).gt.LSIDE(1)) then
             theta=LSIDE(2)/(LSIDE(2)-LSIDE(1))
             facespecies_local(imspec)= &
               theta/spec1(imspec)+(one-theta)/spec2(imspec)
             facespecies_local(imspec)=one/facespecies_local(imspec)
            else if (LSIDE(1).gt.LSIDE(2)) then
             theta=LSIDE(1)/(LSIDE(1)-LSIDE(2))
             facespecies_local(imspec)= &
               theta/spec1(imspec)+(one-theta)/spec2(imspec)
             facespecies_local(imspec)=one/facespecies_local(imspec)
            else
             print *,"LSIDE bust"
             stop
            endif

            spec_test=speciesvisc_interface((imspec-1)*nten+iten)
            if (spec_test.eq.zero) then
             ! do nothing
            else if (spec_test.gt.zero) then
             facespecies_local(imspec)=spec_test
            else
             print *,"spec_test invalid"
             stop
            endif

           enddo !imspec=1..num_species_var

          else
           print *,"gradh bust"
           stop
          endif

          if ((diffusionface_flag.eq.0).and. &
              (temperatureface_flag.eq.0)) then
           ! do nothing (use LS)
          else if ((diffusionface_flag.eq.1).or. &
                   (temperatureface_flag.eq.1)) then ! use VOF

           ! mu_face=sum F_i/(sum F_i/mu_i)
           voltotal=zero
           visc_total=zero
           heat_total=zero
           do im=1,nmat
            FFACE(im)=zero
           enddo
           do imspec=1,num_species_var
            spec_total(imspec)=zero
           enddo
           is_zero_visc=0
           is_zero_heat=0
           do imspec=1,num_species_var
            is_zero_spec(imspec)=0
           enddo
           do iside=0,1
           do im=1,nmat
            voldepart=local_face(vofface_index+2*(im-1)+iside+1)
            voltotal=voltotal+voldepart
            FFACE(im)=FFACE(im)+voldepart
            heat1=get_user_heatviscconst(im)
            do imspec=1,num_species_var
             spec1(imspec)=fort_speciesviscconst((imspec-1)*nmat+im)
            enddo
            if (iside.eq.0) then
             visc1=localvisc_minus(im)
            else if (iside.eq.1) then
             visc1=localvisc_plus(im)
            else
             print *,"iside invalid"
             stop
            endif
            if (voldepart.gt.zero) then
             if (visc1.eq.zero) then
              is_zero_visc=1
             else if (visc1.gt.zero) then
              visc_total=visc_total+voldepart/visc1
             else
              print *,"visc1 invalid"
              stop
             endif
             if (heat1.eq.zero) then
              is_zero_heat=1
             else if (heat1.gt.zero) then
              heat_total=heat_total+voldepart/heat1
             else
              print *,"heat1 invalid"
              stop
             endif
             do imspec=1,num_species_var
              if (spec1(imspec).eq.zero) then
               is_zero_spec(imspec)=1
              else if (spec1(imspec).gt.zero) then
               spec_total(imspec)=spec_total(imspec)+voldepart/spec1(imspec)
              else
               print *,"spec1 invalid"
               stop
              endif
             enddo ! imspec
            else if (voldepart.eq.zero) then
             ! do nothing
            else 
             print *,"voldepart invalid"
             stop
            endif   
           enddo ! im=1..nmat
           enddo ! iside=0..1

           if (voltotal.gt.zero) then
            do im=1,nmat
             FFACE(im)=FFACE(im)/voltotal
            enddo
           else
            print *,"voltotal invalid voltotal= ",voltotal
            stop
           endif

           if (diffusionface_flag.eq.1) then ! use VFRAC
            if (is_zero_visc.eq.0) then
             facevisc_local=voltotal/visc_total
            else if (is_zero_visc.eq.1) then
             facevisc_local=zero
            else
             print *,"is_zero_visc invalid"
             stop
            endif
           else if (diffusionface_flag.eq.0) then ! use LS
            ! do nothing
           else
            print *,"diffusionface_flag invalid"
            stop
           endif

           if (temperatureface_flag.eq.0) then !GFM (use LS)

            ! do nothing

           else if (temperatureface_flag.eq.1) then ! FVM (use VFRAC)

            if (is_zero_heat.eq.0) then
             faceheat_local=voltotal/heat_total
            else if (is_zero_heat.eq.1) then
             faceheat_local=zero
            else
             print *,"is_zero_heat invalid"
             stop
            endif

            do imspec=1,num_species_var
             if (is_zero_spec(imspec).eq.0) then
              facespecies_local(imspec)=voltotal/spec_total(imspec)
             else if (is_zero_spec(imspec).eq.1) then
              facespecies_local(imspec)=zero
             else
              print *,"is_zero_spec invalid"
              stop
             endif
            enddo ! imspec

           else
            print *,"temperatureface_flag invalid"
            stop
           endif

           do im=1,nmat
            do im_opp=im+1,nmat 
 
             if ((FFACE(im).gt.VOFTOL).and. &
                 (FFACE(im_opp).gt.VOFTOL)) then
              call get_iten(im,im_opp,iten,nmat)

              if (diffusionface_flag.eq.1) then
               if (visc_interface(iten).eq.zero) then
                ! do nothing
               else if (visc_interface(iten).gt.zero) then
                facevisc_local=visc_interface(iten)
               else
                print *,"visc_interface invalid"
                stop
               endif
              else if (diffusionface_flag.eq.0) then
               ! do nothing
              else
               print *,"diffusionface_flag invalid"
               stop
              endif

              if (temperatureface_flag.eq.0) then ! use LS
               ! do nothing
              else if (temperatureface_flag.eq.1) then ! use VFRAC

               if (heatvisc_interface(iten).eq.zero) then
                ! do nothing
               else if (heatvisc_interface(iten).gt.zero) then
                faceheat_local=heatvisc_interface(iten)
               else
                print *,"heatvisc_interface invalid"
                stop
               endif
               do imspec=1,num_species_var
                spec_test=speciesvisc_interface((imspec-1)*nten+iten)
                if (spec_test.eq.zero) then
                 ! do nothing
                else if (spec_test.gt.zero) then
                 facespecies_local(imspec)=spec_test
                else
                 print *,"spec_test invalid"
                 stop
                endif
               enddo !imspec
              else
               print *,"temperatureface_flag invalid"
               stop
              endif
             else if ((FFACE(im).gt.-VOFTOL).and. &
                      (FFACE(im_opp).gt.-VOFTOL)) then
              ! do nothing
             else
              print *,"FFACE invalid"
              stop
             endif

            enddo ! im_opp=1..nmat
           enddo ! im=1..nmat

          else
           print *,"diffusionface_flag or temperatureface_flag invalid"
           stop
          endif
  
         else
          print *,"solid_present_flag bust"
          stop
         endif

          ! initialize wall_flag_face

         wall_flag_face=0
         wallVOF_face=zero

         do iside=0,1

          icell=i+ii*(iside-1)  ! i-ii,i
          jcell=j+jj*(iside-1)  ! j-jj,j

          if (SDIM.eq.3) then
           kcell=k+kk*(iside-1)  ! k-kk,k
          else if (SDIM.eq.2) then
           kcell=0
          else
           print *,"dimension bust"
           stop
          endif

          call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)

          if (levelrz.eq.0) then
           ! do nothing
          else if (levelrz.eq.1) then
           if (SDIM.eq.2) then
            if (veldir.eq.0) then
             if (icell.eq.-1) then 
              if (xsten(0,1).lt.zero) then ! reflecting BC
               icell=0
              endif
             else if (icell.ge.0) then
              if (xsten(0,1).gt.zero) then
               ! do nothing
              else
               print *,"xsten bust"
               stop
              endif
             else
              print *,"icell invalid"
              stop
             endif
            else if (veldir.eq.1) then
             if (xsten(0,1).gt.zero) then
              ! do nothing
             else
              print *,"xsten bust"
              stop
             endif
            else 
             print *,"veldir invalid"
             stop
            endif
           else
            print *,"dimension bust"
            stop
           endif
          else if (levelrz.eq.3) then
           if (veldir.eq.0) then
            if (icell.eq.-1) then 
             if (xsten(0,1).lt.zero) then ! reflecting BC
              icell=0
             endif
            else if (icell.ge.0) then
             if (xsten(0,1).gt.zero) then
              ! do nothing
             else
              print *,"xsten bust"
              stop
             endif
            else
             print *,"icell invalid"
             stop
            endif
           else if ((veldir.eq.1).or.(veldir.eq.SDIM-1)) then
            if (xsten(0,1).gt.zero) then
             ! do nothing
            else
             print *,"xsten bust"
             stop
            endif
           else 
            print *,"veldir invalid"
            stop
           endif
          else
           print *,"levelrz invalid init physics vars 2"
           stop
          endif

          if (levelrz.eq.0) then
           ! do nothing
          else if ((levelrz.eq.1).or. &
                   (levelrz.eq.3)) then
           if (veldir.eq.0) then
            if (iside.eq.0) then
             if (abs(xstenMAC(0,1)).le.VOFTOL*dx(1)) then
              wall_flag_face=nmat+1 !Neumann BC, RZ axis of symmetry
             else if (xstenMAC(0,1).gt.VOFTOL*dx(1)) then
              ! do nothing
             else
              print *,"xstenMAC invalid"
              stop
             endif
            else if (iside.eq.1) then
             if (xstenMAC(0,1).gt.-VOFTOL*dx(1)) then
              ! do nothing
             else
              print *,"xstenMAC invalid"
              stop
             endif
            else
             print *,"iside invalid"
             stop
            endif
           else if ((veldir.eq.1).or.(veldir.eq.SDIM-1)) then
            ! do nothing
           else
            print *,"veldir invalid"
            stop
           endif
          else
           print *,"levelrz invalid init physics vars 3"
           stop
          endif

          do im=1,nmat
           volmat(im)=vofC(D_DECL(icell,jcell,kcell),im)
          enddo
           ! voltotal=sum F_prescribed
           ! im_prescribed_primary=argmax_im F_prescribed
          call combine_prescribed_VOF(volmat,nmat,voltotal, &
                 im_prescribed_primary)
          if (voltotal.ge.one-0.01D0) then
           if ((im_prescribed_primary.ge.1).and. &
               (im_prescribed_primary.le.nmat)) then
            if (wall_flag_face.eq.0) then
             wall_flag_face=im_prescribed_primary !Neumann interior wall
             wallVOF_face=volmat(im_prescribed_primary)
            else if ((wall_flag_face.ge.1).and. &
                     (wall_flag_face.le.nmat)) then
             if (volmat(im_prescribed_primary).gt.wallVOF_face) then
              wall_flag_face=im_prescribed_primary
              wallVOF_face=volmat(im_prescribed_primary)
             endif
            else if (wall_flag_face.eq.nmat+1) then
             ! do nothing
            else
             print *,"wall_flag_face invalid"
             stop
            endif
           else
            print *,"im_prescribed_primary invalid"
            stop
           endif
          endif

          if ((wall_flag_face.ge.0).and. &
              (wall_flag_face.le.nmat)) then

           prescribed_flag=prescribed_exists(nmat)

           if (prescribed_flag.eq.1) then 

             ! at least one of the face's adjoining cells is dominated
             ! by a prescribed material.
            if (is_prescribed_face.eq.1) then 
             if ((im_prescribed.ge.1).and. &
                 (im_prescribed.le.nmat)) then
              if (wall_flag_face.eq.0) then
               wall_flag_face=im_prescribed
               wallVOF_face=volmat(im_prescribed)
              else if ((wall_flag_face.ge.1).and. &
                       (wall_flag_face.le.nmat)) then
               if (volmat(im_prescribed).gt.wallVOF_face) then
                wall_flag_face=im_prescribed
                wallVOF_face=volmat(im_prescribed)
               endif
              else
               print *,"wall_flag_face invalid"
               stop
              endif
             else
              print *,"im_prescribed invalid"
              stop
             endif
            else if (is_prescribed_face.eq.0) then
             ! do nothing
            else
             print *,"is_prescribed_face invalid"
             stop
            endif

           else if (prescribed_flag.eq.0) then
            ! do nothing
           else
            print *,"prescribed_flag invalid 70"
            stop
           endif

          else if (wall_flag_face.eq.nmat+1) then
           ! do nothing
          else
           print *,"wall_flag_face invalid"
           stop
          endif

         enddo ! iside=0..1   (initializing wall_flag_side)

         voltotal=zero
         mass_total=zero
         do im=1,nmat
          FFACE(im)=zero
         enddo

         do iside=0,1
         do im=1,nmat

          voldepart=local_face(vofface_index+2*(im-1)+iside+1)
          voltotal=voltotal+voldepart
          FFACE(im)=FFACE(im)+voldepart
          mass_total=mass_total+  &
           local_face(massface_index+2*(im-1)+iside+1)
          
         enddo
         enddo ! iside=0..1

         if ((mass_total.le.zero).or. &
             (voltotal.le.zero)) then
          print *,"mass_total or voltotal invalid"
          print *,"i,j,k,veldir ",i,j,k,veldir
          print *,"vol,mass,iside ",voltotal,mass_total,iside
          do dir2=1,SDIM
           print *,"dir2,fablo,fabhi ",dir2,fablo(dir2),fabhi(dir2)
          enddo
          stop
         endif

         if (voltotal.gt.zero) then
          do im=1,nmat
           FFACE(im)=FFACE(im)/voltotal
          enddo
         else 
          print *,"voltotal invalid (faceden) voltotal= ",voltotal
          stop
         endif

         density_for_mass_fraction_diffusion=mass_total/voltotal

         local_face(faceden_index+1)=one/density_for_mass_fraction_diffusion

         do im=1,nmat
          do im_opp=im+1,nmat
           if ((FFACE(im).gt.VOFTOL).and. &
               (FFACE(im_opp).gt.VOFTOL)) then
            call get_iten(im,im_opp,iten,nmat)
            if (den_interface(iten).eq.zero) then
             ! do nothing
            else if (den_interface(iten).gt.zero) then
             local_face(faceden_index+1)=one/den_interface(iten)
            else
             print *,"den_interface invalid"
             stop
            endif
           else if ((FFACE(im).gt.-VOFTOL).and. &
                    (FFACE(im_opp).gt.-VOFTOL)) then
            ! do nothing
           else
            print *,"FFACE invalid"
            stop
           endif
          enddo ! im_opp=im+1..nmat
         enddo ! im=1..nmat

         local_face(facevisc_index+1)=facevisc_local
         local_face(faceheat_index+1)=faceheat_local
         do imspec=1,num_species_var
          local_face(facespecies_index+imspec)= &
            density_for_mass_fraction_diffusion*facespecies_local(imspec)
         enddo


          ! mask_boundary=1 at left neumann boundary
          ! mask_boundary=2 at right neumann boundary
          ! mask_boundary=0 otherwise
         if (mask_boundary.eq.1) then
          wall_flag_face=nmat+1
         else if (mask_boundary.eq.2) then
          wall_flag_face=nmat+1
         else if (mask_boundary.eq.0) then
          ! do nothing
         else 
          print *,"mask_boundary invalid"
          stop
         endif  

         local_face(icefacecut_index+1)=one

! local_face(facecut_index+1)=0.0 if presbc=REFLECT_EVEN,LO_EXTRAP
! local_face(facecut_index+1)=prescribed_solid_scale if face has adjoining
!  prescribed solid.
         if (wall_flag_face.eq.0) then
          local_face(facecut_index+1)=one
         else if ((wall_flag_face.ge.1).and. &
                  (wall_flag_face.le.nmat)) then
          local_face(facecut_index+1)=prescribed_solid_scale(wall_flag_face)
         else if (wall_flag_face.eq.nmat+1) then
          local_face(facecut_index+1)=zero
         else
          print *,"wall_flag_face invalid"
          stop
         endif

         ! compute curvature at the faces
         ! also compute model based pressure forcing at the faces

         zeroradius_flag=0

         if (levelrz.eq.0) then
          ! do nothing
         else if (levelrz.eq.1) then
          if (SDIM.eq.2) then
           if (veldir.eq.0) then
            if (abs(xstenMAC(0,1)).le.VOFTOL*dx(1)) then ! at r=0 face
             zeroradius_flag=1
            else if (xstenMAC(0,1).ge.-VOFTOL*dx(1)) then
             ! do nothing
            else
             print *,"xstenMAC invalid"
             stop
            endif
           else if (veldir.eq.1) then
            ! do nothing
           else
            print *,"veldir invalid"
            stop
           endif
          else
           print *,"dimension bust"
           stop
          endif
         else if (levelrz.eq.3) then
          ! do nothing
         else
          print *,"levelrz invalid init physics vars 4"
          stop
         endif

         local_face(curv_index+1)=zero
         local_face(pforce_index+1)=zero

         if (local_face(facecut_index+1).lt.zero) then
          print *,"local_face(facecut_index+1).lt.zero"
          stop
         else if ((local_face(facecut_index+1).ge.zero).and. &
                  (local_face(facecut_index+1).le.half)) then
          ! do nothing
         else if (zeroradius_flag.eq.1) then
          ! do nothing
         else if ((zeroradius_flag.eq.0).and. &
                  (local_face(facecut_index+1).ge.half).and. &
                  (local_face(facecut_index+1).le.one)) then

          if (gradh_tension.ne.zero) then

           if ((im_tension.gt.nmat).or.(im_opp_tension.gt.nmat).or. &
               (im_tension.lt.1).or.(im_opp_tension.lt.1)) then
            print *,"im_tension or im_opp_tension bust 4"
            stop
           endif
           call get_iten(im_tension,im_opp_tension,iten_tension,nmat)
  
! vof,ref centroid,order,slope,intercept  x num_materials

           if (LS_consistent_tension.eq.1) then

            if (is_solid_face.eq.1) then
             local_face(curv_index+1)=zero
             local_face(pforce_index+1)=zero
            else if (is_solid_face.eq.0) then
             icurv=(iten_tension-1)*(SDIM+5)

             if ((level.ge.0).and.(level.le.finest_level)) then

              if ((project_option.eq.0).or. &
                  (project_option.eq.1)) then

                ! curv_cellHT
                ! curv_cellFD
                ! pforce_cell
                ! mgoni_force(1..sdim)
                ! dir x side  dir=1..sdim  side=-1 or 1
                ! im3
               do icurv_ofs=1,SDIM+5
                curvL(icurv_ofs)=curv(D_DECL(im1,jm1,km1),icurv+icurv_ofs)
                curvR(icurv_ofs)=curv(D_DECL(i,j,k),icurv+icurv_ofs)
               enddo
               dirL=NINT(curvL(4+SDIM))
               sideL=1
               if (dirL.lt.0) then
                dirL=-dirL
                sideL=-sideL
               endif
               dirR=NINT(curvR(4+SDIM))
               sideR=1
               if (dirR.lt.0) then
                dirR=-dirR
                sideR=-sideR
               endif
               im3L=NINT(curvL(5+SDIM))
               im3R=NINT(curvR(5+SDIM))

                ! dirL or dirR = sdim+1 if curvature record averaged down.
               if ((dirL.eq.SDIM+1).or.(dirL.eq.0)) then
                orientL=0
               else if ((dirL.eq.veldir+1).and.(sideL.eq.1)) then
                orientL=1
               else if ((dirL.ge.1).and.(dirL.le.SDIM)) then
                orientL=0
               else
                print *,"dirL invalid"
                stop
               endif

               if ((dirR.eq.SDIM+1).or.(dirR.eq.0)) then
                orientR=0
               else if ((dirR.eq.veldir+1).and.(sideR.eq.-1)) then
                orientR=1
               else if ((dirR.ge.1).and.(dirR.le.SDIM)) then
                orientR=0
               else
                print *,"dirR invalid"
                stop
               endif
               
                ! note: gradh_tension=0.0 if covered_face==2.

               if ((dirL.eq.0).or. &
                   (dirR.eq.0).or. &
                   (dirL.eq.SDIM+1).or. &
                   (dirR.eq.SDIM+1)) then
                if ((coarse_fine_face.eq.1).or. & ! INT_DIR, coarse-fine
                    (covered_face.eq.1)) then     ! neighbor is covered.
                 im3L=0
                 im3R=0 
                 if ((dirL.eq.0).or. &
                     (dirR.eq.0)) then
                  do icurv_ofs=1,SDIM+5
                   curvL(icurv_ofs)=zero
                   curvR(icurv_ofs)=zero
                  enddo
                 endif
                else
                 print *,"dirL or dirR invalid"
                 stop
                endif
               endif

               if ((im3L.lt.0).or.(im3L.gt.nmat).or. &
                   (im3R.lt.0).or.(im3R.gt.nmat).or. &
                   (dirL.lt.0).or.(dirL.gt.SDIM+1).or. &
                   (dirR.lt.0).or.(dirR.gt.SDIM+1)) then
                print *,"curvature parameters invalid"
                print *,"im3L,im3R ",im3L,im3R
                print *,"dirL,dirR ",dirL,dirR
                print *,"im_tension ",im_tension
                print *,"im_opp_tension ",im_opp_tension
                print *,"iten_tension ",iten_tension
                print *,"i,j,k,veldir ",i,j,k,veldir
                do im=1,nmat
                 print *,"im,LS_LEFT,LS_RIGHT ",im, &
                  levelPC(D_DECL(im1,jm1,km1),im),levelPC(D_DECL(i,j,k),im)
                enddo
                stop
               endif

                ! inputs.curvature_converge, continuous_mof=2
                ! March 10, 2018: 1.99, 2.03 RZ 24x48 HT
                ! March 10, 2018: 1.00, 1.01 XY 24x48 HT
                ! March 10, 2018: 1.93, 2.07 XYZ 32x32x32 HT
               if (((im3L.ge.1).and.(im3L.le.nmat)).or. &
                   ((im3R.ge.1).and.(im3R.le.nmat))) then
                if (FD_curv_interp.eq.1) then
                 curv_interp_flag=0 ! interpolate curvFD
                else if (FD_curv_interp.eq.0) then
                 curv_interp_flag=5 ! closest curvFD
                else
                 print *,"FD_curv_interp invalid"
                 stop
                endif 
               else if ((im3L.eq.0).and.(im3R.eq.0)) then
                if ((orientL.eq.1).and.(orientR.eq.1)) then
                 curv_interp_flag=1 ! closest curvHT
                else if ((orientL.eq.1).and.(orientR.eq.0)) then
                 curv_interp_flag=2 ! left curvHT
                else if ((orientL.eq.0).and.(orientR.eq.1)) then
                 curv_interp_flag=3 ! right curvHT
                else if ((orientL.eq.0).and.(orientR.eq.0)) then
                 curv_interp_flag=4 ! interpolate curvHT
                else
                 print *,"orientL or orientR invalid"
                 stop
                endif
               else
                print *,"im3L or im3R invalid"
                stop
               endif

               wtL=abs(LSIDE_tension(2))
               wtR=abs(LSIDE_tension(1))
               wtsum=wtL+wtR
               if (wtsum.eq.zero) then
                wtL=half
                wtR=half
               else if (wtsum.gt.zero) then 
                wtL=wtL/wtsum
                wtR=wtR/wtsum
               else
                print *,"wtsum invalid"
                stop
               endif 
                
               if (curv_interp_flag.eq.0) then
                local_face(curv_index+1)=wtL*curvL(2)+wtR*curvR(2)
               else if (curv_interp_flag.eq.1) then   
                if (wtL.gt.wtR) then
                 local_face(curv_index+1)=curvL(1)
                else if (wtR.gt.wtL) then
                 local_face(curv_index+1)=curvR(1)
                else 
                 local_face(curv_index+1)=wtL*curvL(1)+wtR*curvR(1)
                endif
               else if (curv_interp_flag.eq.2) then
                local_face(curv_index+1)=curvL(1)
               else if (curv_interp_flag.eq.3) then
                local_face(curv_index+1)=curvR(1)
               else if (curv_interp_flag.eq.4) then
                local_face(curv_index+1)=wtL*curvL(1)+wtR*curvR(1)
               else if (curv_interp_flag.eq.5) then
                if (wtL.gt.wtR) then
                 local_face(curv_index+1)=curvL(2)
                else if (wtR.gt.wtL) then
                 local_face(curv_index+1)=curvR(2)
                else 
                 local_face(curv_index+1)=wtL*curvL(2)+wtR*curvR(2)
                endif
               else
                print *,"curv_interp_flag invalid"
                stop
               endif

               if (curv_min.gt.local_face(curv_index+1)) then
                curv_min=local_face(curv_index+1)
               endif
               if (curv_max.lt.local_face(curv_index+1)) then
                curv_max=local_face(curv_index+1)
               endif

               if (wtL.gt.wtR) then
                local_face(pforce_index+1)=curvL(3)
               else if (wtR.gt.wtL) then
                local_face(pforce_index+1)=curvR(3)
               else 
                local_face(pforce_index+1)=wtL*curvL(3)+wtR*curvR(3)
               endif

              else if (project_option.eq.10) then
               ! do nothing
              else if ((project_option.eq.11).or. & ! FSI_material_exists 2nd 
                       (project_option.eq.13)) then ! FSI_material_exists 1st
               print *,"FORT_INIT_PHYSICS_VARS should not be called here"
               stop
              else
               print *,"project_option invalid"
               stop
              endif

             else
              print *,"level invalid init physics vars 2"
              stop
             endif

            else
             print *,"is_solid_face invalid 5 ",is_solid_face
             stop
            endif

           else if (LS_consistent_tension.eq.0) then
            print *,"WARNING inconsistent signs for levelset function"
            print *,"i,j,k,veldir ",i,j,k,veldir
            print *,"level,finest_level ",level,finest_level
            print *,"im_tension,im_opp_tension ",im_tension,im_opp_tension
            print *,"gradh_tension=",gradh_tension
            print *,"LSIDE_tension(1) ",LSIDE_tension(1) 
            print *,"LSIDE_tension(2) ",LSIDE_tension(2) 
            do im=1,nmat
             print *,"im,lsminus,lsplus ",im,LSminus(im),LSplus(im)
            enddo
            print *,"level,finest_level ",level,finest_level
           else
            print *,"LS_consistent_tension invalid"
            stop
           endif

          else if (gradh_tension.eq.zero) then
           ! do nothing
          else
           print *,"gradh_tension bust"
           stop
          endif

         else
          print *,"zeroradius_flag or local_face(facecut_index+1) invalid"
          stop
         endif 

         do im=1,ncphys
          if (veldir.eq.0) then
           xface(D_DECL(i,j,k),im)=local_face(im)
          else if (veldir.eq.1) then
           yface(D_DECL(i,j,k),im)=local_face(im)
          else if ((veldir.eq.2).and.(SDIM.eq.3)) then
           zface(D_DECL(i,j,k),im)=local_face(im)
          else
           print *,"veldir invalid"
           stop
          endif
         enddo  ! im=1,ncphys

        enddo
        enddo
        enddo  ! i,j,k
       enddo ! veldir

      else if (isweep.eq.1) then


        ! cenden, cenvof, cenDeDT, cenvisc, initialized in this loop.

       call growntilebox(tilelo,tilehi,fablo,fabhi,igridlo,igridhi,1) 

       do i=igridlo(1),igridhi(1)
       do j=igridlo(2),igridhi(2)
       do k=igridlo(3),igridhi(3)

        call gridsten_level(xsten,i,j,k,level,nhalf)
        do im=1,nmat*ngeom_recon
         mofdata(im)=slope(D_DECL(i,j,k),im)
        enddo
         ! before (mofdata): fluids tessellate
         ! after  (mofdata): fluids and solids tessellate
        call multi_get_volume_tessellate( &
         bfact,dx,xsten,nhalf, &
         mofdata, &
         geom_xtetlist(1,1,1,tid+1), &
         nmax, &
         nmax, &
         nmat,SDIM,3)

        voltotal=zero
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         volmat(im)=mofdata(vofcomp)
         voltotal=voltotal+volmat(im)
        enddo ! im=1..nmat

        if (abs(voltotal-one).gt.LSTOL) then
         print *,"voltotal invalid"
         stop
        endif

        DeDT_total=zero

        null_viscosity=0
        visc_total=zero

        mass_total=zero

        do im=1,nmat

         dencomp=(im-1)*num_state_material+1
         tempcomp=dencomp+1

         delta_mass=denstate(D_DECL(i,j,k),dencomp)*volmat(im)
         voldepart=volmat(im)
          ! if is_rigid(im), density=fort_denconst(im)
          ! if incompressible,
          !   if override_density=0,2 then density=fort_denconst(im)
          !   if override_density=1 then density=mass_depart/vol_depart
          ! if compressible,
          !   density=massdepart/voltarget
         call derive_density( &
          voldepart,voldepart,voltotal, &
          override_density,delta_mass, &
          im,nmat,den, &
          for_hydrostatic_pres) 

         if (den.le.zero) then
          print *,"density must be positive init_phyiscs_vars2"
          print *,"i,j,k,im,den ",i,j,k,im,den
          stop
         endif

         total_density=den

         ispecies=cavitation_species(im)
         if (ispecies.eq.0) then
          ! do nothing
         else if ((ispecies.ge.1).and.(ispecies.le.num_species_var)) then
          speccomp=(im-1)*num_state_material+num_state_base+ispecies
          local_mfrac=denstate(D_DECL(i,j,k),speccomp)
          if ((local_mfrac.ge.zero).and.(local_mfrac.le.one)) then
           evap_den=cavitation_vapor_density(im)
           if (evap_den.gt.zero) then
            total_density=den*evap_den/ &
             (local_mfrac*den+(one-local_mfrac)*evap_den)
            den=total_density
           else
            print *,"evap_den invalid"
            stop
           endif
          else
           print *,"local_mfrac invalid 1"
           print *,"local_mfrac=",local_mfrac
           stop
          endif
         else
          print *,"ispecies invalid"
          stop
         endif
 
         do ispecies=1,num_species_var
          im_species=spec_material_id(ispecies)
          if (im_species.eq.0) then
           ! do nothing
          else if (im_species.eq.im) then
           ! do nothing
          else if ((im_species.ge.1).and. &
                   (im_species.le.nmat).and. &
                   (im_species.ne.im)) then
           speccomp=(im-1)*num_state_material+num_state_base+ispecies
           local_mfrac=denstate(D_DECL(i,j,k),speccomp)
           if ((local_mfrac.ge.zero).and.(local_mfrac.le.one)) then
            evap_den=species_evaporation_density(ispecies)
            if (evap_den.gt.zero) then
             total_density=den*evap_den/ &
              (local_mfrac*den+(one-local_mfrac)*evap_den)
            else
             print *,"evap_den invalid"
             stop
            endif
           else
            print *,"local_mfrac invalid 2"
            print *,"local_mfrac=",local_mfrac
            stop
           endif
          else
           print *,"im_species (spec_material_id(ispecies)) invalid1"
           print *,"im_species= ",im_species
           print *,"ispecies=",ispecies
           stop
          endif
         enddo ! ispecies=1,num_species_var

         den=total_density

         cenden(D_DECL(i,j,k),im+1)=one/den 

         localvisc(im)=modvisc(D_DECL(i,j,k),im)
         if (localvisc(im).lt.zero) then
          print *,"viscstate gone negative"
          stop
         else if (localvisc(im).eq.zero) then
          if (volmat(im).gt.zero) then
           null_viscosity=1
          endif
         else if (localvisc(im).gt.zero) then
          ! do nothing
         else
          print *,"localvisc invalid"
          stop
         endif
         cenvisc(D_DECL(i,j,k),im+1)=localvisc(im)

         one_over_mu=one/(localvisc(im)+VISCINVTOL)
         visc_total=visc_total+one_over_mu*volmat(im)

         imattype=fort_material_type(im)
         TEMPERATURE=denstate(D_DECL(i,j,k),tempcomp)
         if (TEMPERATURE.le.zero) then
          print *,"PHYSICS_VAR: temperature must be positive"
          print *,"num_materials ",num_materials
          print *,"num_state_material ",num_state_material
          print *,"fort_initial_temperature(im) ",fort_initial_temperature(im)
          print *,"fort_tempconst(im) ",fort_tempconst(im)
          print *,"i,j,k,im ",i,j,k,im
          print *,"den,TEMPERATURE ",den,TEMPERATURE
          stop
         endif
          ! DeDT = cv
         call DeDT_material(den,TEMPERATURE,DeDT,imattype,im)
         if (DeDT.le.zero) then
          print *,"DeDT must be positive"
          stop
         endif
         cenDeDT(D_DECL(i,j,k),im+1)=one/(den*DeDT)

         delta_mass=den*volmat(im)
         mass_total=mass_total+delta_mass
         DeDT_total=DeDT_total+DeDT*delta_mass
        enddo ! im=1..nmat

        if (mass_total.le.zero) then
         print *,"mass_total invalid"
         stop
        endif

        if (DeDT_total.le.zero) then
         print *,"DeDT_total must be positive"
         stop
        endif

        do im=1,nmat
         cenvof(D_DECL(i,j,k),im)=volmat(im)/voltotal
        enddo

        cenden(D_DECL(i,j,k),1)=voltotal/mass_total 
        cenDeDT(D_DECL(i,j,k),1)=voltotal/DeDT_total

        if (null_viscosity.eq.1) then
         cenvisc(D_DECL(i,j,k),1)=zero
        else if (null_viscosity.eq.0) then
         cenvisc(D_DECL(i,j,k),1)= &
          one/((visc_total/voltotal)+VISCINVTOL)  ! mu
        else
         print *,"null_viscosity invalid"
         stop
        endif

       enddo
       enddo
       enddo  ! i,j,k (cell centered quantities)

      else
       print *,"isweep invalid"
       stop
      endif

      if (DO_SANITY_CHECK.eq.1) then
       call mass_face_weight()
       allocate(comparemassface(fablo(1)-1:fabhi(1)+1,5))
       allocate(comparedenface(fablo(1)-1:fabhi(1)+1,5))
       j=1
       k=0
       do i=fablo(1),fabhi(1)+1
        comparemassface(i,1)=xface(D_DECL(i,j,k),massface_index+1)/dx(2)
        comparemassface(i,2)=xface(D_DECL(i,j,k),massface_index+2)/dx(2)
        comparedenface(i,1)=one/xface(D_DECL(i,j,k),faceden_index+1)
       enddo
       call compare_sanity(comparemassface,1,2,3)
       call compare_sanity(comparedenface,1,1,4)
       deallocate(comparemassface)
       deallocate(comparedenface)
      endif
 
      return
      end subroutine FORT_INIT_PHYSICS_VARS


      subroutine FORT_BUILD_SEMIREFINEVOF( &
       tid, &
       tessellate, &
       ngrow_refine, &
       nrefine_vof, &
       nrefine_cen, &
       nten, &
       spec_material_id, &
       mass_fraction_id, &
       species_evaporation_density, &
       cavitation_vapor_density, &
       cavitation_species, &
       override_density, &
       xlo,dx, &
       slope,DIMS(slope), &
       denstate,DIMS(denstate), &
       vofF,DIMS(vofF), &
       cenF,DIMS(cenF), &
       massF,DIMS(massF), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       nmat, &
       level,finest_level)
      use global_utility_module
      use probf90_module
      use geometry_intersect_module
      use MOF_routines_module
      use godunov_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tessellate
      INTEGER_T, intent(in) :: ngrow_refine
      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: nrefine_vof
      INTEGER_T, intent(in) :: nrefine_cen
      INTEGER_T, intent(in) :: nten
      INTEGER_T :: nten_test
      INTEGER_T, intent(in) :: spec_material_id(num_species_var+1)
      INTEGER_T, intent(in) :: mass_fraction_id(2*nten)
      REAL_T, intent(in) :: species_evaporation_density(num_species_var+1)
      REAL_T, intent(in) :: cavitation_vapor_density(nmat)
      INTEGER_T, intent(in) :: cavitation_species(nmat)
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nmat
      INTEGER_T :: veldir
      INTEGER_T, intent(in) :: override_density(nmat)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(slope)
      INTEGER_T, intent(in) :: DIMDEC(denstate)
      INTEGER_T, intent(in) :: DIMDEC(vofF)
      INTEGER_T, intent(in) :: DIMDEC(cenF)
      INTEGER_T, intent(in) :: DIMDEC(massF)
      REAL_T, intent(in) :: slope(DIMV(slope),nmat*ngeom_recon) 
      REAL_T, intent(in) :: denstate(DIMV(denstate),nmat*num_state_material) 
      REAL_T, intent(out) :: vofF(DIMV(vofF),nrefine_vof)
      REAL_T, intent(out) :: cenF(DIMV(cenF),nrefine_cen)
      REAL_T, intent(out) :: massF(DIMV(massF),nrefine_vof)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k
      INTEGER_T dir2
      INTEGER_T iside

      INTEGER_T im,nmax
      REAL_T mofdata(nmat*ngeom_recon)

      REAL_T multi_volume(nmat)
      REAL_T multi_cen(SDIM,nmat)
      REAL_T voldepart
      REAL_T voltotal_fluid
      REAL_T mass_total_fluid
      REAL_T voltotal_solid
      REAL_T mass_total_solid
      REAL_T voltotal
      REAL_T mass_total
      REAL_T delta_mass
      REAL_T den

      INTEGER_T dencomp

      INTEGER_T igridlo(3),igridhi(3)

      REAL_T volrecon
      REAL_T cenrecon(SDIM)

      REAL_T voldonate
      REAL_T cendonate(SDIM)

      INTEGER_T check_donate
      INTEGER_T irefine
      INTEGER_T irefinecen
      REAL_T xsten_recon(-1:1,SDIM)
      REAL_T xsten_donate(-1:1,SDIM)
      INTEGER_T for_hydrostatic_pres
      REAL_T mu
      REAL_T total_density,evap_den,local_mfrac
      INTEGER_T ispecies,im_species,speccomp

      for_hydrostatic_pres=0

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif
      if ((tessellate.ne.0).and.(tessellate.ne.1)) then
       print *,"tessellate invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level or finest_level invalid build semi refine vof"
       stop
      endif
      if (ngrow_refine.lt.1) then
       print *,"ngrow_refine invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid semirefine_vof "
       print *,"ngeom_recon= ",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid semirefine_vof "
       print *,"ngeom_raw= ",ngeom_raw
       stop
      endif
      if (nrefine_vof.ne.2*nmat*SDIM) then
       print *,"nrefine_vof invalid"
       stop
      endif
      if (nrefine_cen.ne.2*nmat*SDIM*SDIM) then
       print *,"nrefine_cen invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid semi refine vof nten nten_test ",nten,nten_test
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: BUILD_SEMIREFINEVOF

      call checkbound(fablo,fabhi,DIMS(slope),ngrow_refine,-1,219)
      call checkbound(fablo,fabhi, &
       DIMS(denstate), &
       ngrow_refine,-1,223)
      call checkbound(fablo,fabhi,DIMS(vofF),ngrow_refine,-1,227)
      call checkbound(fablo,fabhi,DIMS(cenF),ngrow_refine,-1,227)
      call checkbound(fablo,fabhi,DIMS(massF),ngrow_refine,-1,227)

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        if (is_rigid(nmat,im).ne.1) then
         print *,"is_rigid(nmat,im).ne.1"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (fort_material_type(im).le.fort_max_num_eos)) then
        ! do nothing
       else 
        print *,"material_type invalid"
        stop
       endif

       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"density must be positive semirefinevof"
        print *,"im,denconst ",im,fort_denconst(im)
        stop
       endif
       if (fort_tempconst(im).gt.zero) then
        ! do nothing
       else
        print *,"semirefinevof:temperature must be positive"
        print *,"im,fort_tempconst : ",im,fort_tempconst(im)
        stop
       endif
       if (fort_energyconst(im).gt.zero) then
        ! do nothing
       else
        print *,"energy must be positive in FORT_BUILD_SEMIREFINEVOF"
        print *,"im= ",im
        print *,"fort_energyconst(im)= ",fort_energyconst(im)
        stop
       endif
       mu=get_user_viscconst(im,fort_denconst(im),fort_tempconst(im))
       if (mu.ge.zero) then
        ! do nothing
       else
        print *,"viscosity cannot be negative"
        stop
       endif

       if ((override_density(im).ne.0).and. &
           (override_density(im).ne.1).and. &
           (override_density(im).ne.2)) then
        print *,"override_density invalid"
        stop
       endif

      enddo ! im  (checking parameters)

      call growntilebox(tilelo,tilehi,fablo,fabhi,igridlo,igridhi, &
        ngrow_refine) 

      do veldir=0,SDIM-1

        do i=igridlo(1),igridhi(1)
        do j=igridlo(2),igridhi(2)
        do k=igridlo(3),igridhi(3)

         call CISBOX(xsten_recon,1, &
          xlo,dx,i,j,k, &
          bfact,level, &
          volrecon,cenrecon,SDIM)

         do iside=-1,1,2
 
          check_donate=1
         
          if (levelrz.eq.0) then
           ! do nothing
          else if ((levelrz.eq.1).or. &
                   (levelrz.eq.3)) then
           if (xsten_recon(0,1).lt.VOFTOL*dx(1)) then
            check_donate=0
           endif
          else
           print *,"levelrz invalid build semi refine vof"
           stop
          endif
 
          if (check_donate.eq.0) then
           do im=1,nmat

            if (iside.eq.-1) then
             irefine=veldir*2*nmat+im
            else if (iside.eq.1) then
             irefine=veldir*2*nmat+nmat+im
            else
             print *,"iside invalid"
             stop
            endif

            vofF(D_DECL(i,j,k),irefine)=zero
            massF(D_DECL(i,j,k),irefine)=zero
            do dir2=1,SDIM

             if (iside.eq.-1) then
              irefinecen=veldir*2*nmat*SDIM+ &
               (im-1)*SDIM+dir2
             else if (iside.eq.1) then
              irefinecen=veldir*2*nmat*SDIM+ &
               nmat*SDIM+(im-1)*SDIM+dir2
             else
              print *,"iside invalid"
              stop
             endif
 
             cenF(D_DECL(i,j,k),irefinecen)=zero
            enddo ! dir2 
           enddo ! im
          else if (check_donate.eq.1) then

           call CISBOXHALF(xsten_donate,1, &
            xlo,dx,i,j,k,iside,veldir+1, &
            bfact,level, & 
            voldonate,cendonate,SDIM)

           do dir2=1,nmat*ngeom_recon
            mofdata(dir2)=slope(D_DECL(i,j,k),dir2)
           enddo 
           ! multi_cen is "absolute" (not relative to cell centroid)
           call multi_get_volume_grid_simple( &
             tessellate, &
             bfact,dx,xsten_recon,1, &
             mofdata, &
             xsten_donate,1, &
             multi_volume,multi_cen, &
             geom_xtetlist(1,1,1,tid+1), &
             nmax, &
             nmax, &
             nmat,SDIM,2)
       
           mass_total_fluid=zero
           voltotal_fluid=zero 
           mass_total_solid=zero
           voltotal_solid=zero 
           do im=1,nmat
            dencomp=(im-1)*num_state_material+1
            den=denstate(D_DECL(i,j,k),dencomp)
            delta_mass=den*multi_volume(im)
            voldepart=multi_volume(im)
            ! if is_rigid(im), density=fort_denconst(im)
            ! if incompressible,
            !   if override_density=0,2 then density=fort_denconst(im)
            !   if override_density=1 then density=mass_depart/vol_depart
            ! if compressible,
            !   density=massdepart/voltarget
            call derive_density( &
             voldepart,voldepart,voldonate, &
             override_density,delta_mass, &
             im,nmat,den,for_hydrostatic_pres)
 
            total_density=den

            ispecies=cavitation_species(im)
            if (ispecies.eq.0) then
             ! do nothing
            else if ((ispecies.ge.1).and.(ispecies.le.num_species_var)) then
             speccomp=(im-1)*num_state_material+num_state_base+ispecies
             local_mfrac=denstate(D_DECL(i,j,k),speccomp)
             if ((local_mfrac.ge.zero).and.(local_mfrac.le.one)) then
              evap_den=cavitation_vapor_density(im)
              if (evap_den.gt.zero) then
               total_density=den*evap_den/ &
                (local_mfrac*den+(one-local_mfrac)*evap_den)
               den=total_density
              else
               print *,"evap_den invalid"
               stop
              endif
             else
              print *,"local_mfrac invalid 3 FORT_BUILD_SEMIREFINEVOF"
              print *,"local_mfrac=",local_mfrac
              print *,"ngrow_refine=",ngrow_refine
              print *,"im=",im
              print *,"speccomp=",speccomp
              print *,"ispecies=",ispecies
              print *,"den=",den
              print *,"i,j,k ",i,j,k
              print *,"level,finest_level ",level,finest_level
              stop
             endif
            else
             print *,"ispecies invalid"
             stop
            endif

            do ispecies=1,num_species_var
             im_species=spec_material_id(ispecies)
             if (im_species.eq.0) then
              ! do nothing
             else if (im_species.eq.im) then
              ! do nothing
             else if ((im_species.ge.1).and. &
                      (im_species.le.nmat).and. &
                      (im_species.ne.im)) then
              speccomp=(im-1)*num_state_material+num_state_base+ispecies
              local_mfrac=denstate(D_DECL(i,j,k),speccomp)
              if ((local_mfrac.ge.zero).and.(local_mfrac.le.one)) then
               evap_den=species_evaporation_density(ispecies)
               if (evap_den.gt.zero) then
                total_density=den*evap_den/ &
                 (local_mfrac*den+(one-local_mfrac)*evap_den)
               else
                print *,"evap_den invalid"
                stop
               endif
              else
               print *,"local_mfrac invalid 4"
               print *,"local_mfrac=",local_mfrac
               stop
              endif
             else
              print *,"im_species (spec_material_id(ispecies)) invalid2"
              print *,"im_species= ",im_species
              print *,"ispecies=",ispecies
              stop
             endif
            enddo ! ispecies=1,num_species_var

            den=total_density

            if (is_rigid(nmat,im).eq.0) then
             voltotal_fluid=voltotal_fluid+multi_volume(im)
             mass_total_fluid=mass_total_fluid+den*multi_volume(im)
            else if (is_rigid(nmat,im).eq.1) then
             voltotal_solid=voltotal_solid+multi_volume(im)
             mass_total_solid=mass_total_solid+den*multi_volume(im)
            else
             print *,"is_rigid invalid"
             stop
            endif
  
            if (iside.eq.-1) then 
             irefine=veldir*2*nmat+im
            else if (iside.eq.1) then
             irefine=veldir*2*nmat+nmat+im
            else
             print *,"iside invalid"
             stop
            endif

            vofF(D_DECL(i,j,k),irefine)=multi_volume(im)
            massF(D_DECL(i,j,k),irefine)=den*multi_volume(im)
            do dir2=1,SDIM
             if (iside.eq.-1) then
              irefinecen=veldir*2*nmat*SDIM+ &
               (im-1)*SDIM+dir2
             else if (iside.eq.1) then
              irefinecen=veldir*2*nmat*SDIM+ &
               nmat*SDIM+(im-1)*SDIM+dir2
             else
              print *,"iside invalid"
              stop
             endif
             cenF(D_DECL(i,j,k),irefinecen)=multi_cen(dir2,im)
            enddo ! dir2 
           enddo ! im=1,nmat

           if (tessellate.eq.0) then
            voltotal=voltotal_fluid
            mass_total=mass_total_fluid
           else if (tessellate.eq.1) then
            voltotal=voltotal_fluid+voltotal_solid
            mass_total=mass_total_fluid+mass_total_solid
           else
            print *,"tessellate invalid"
            stop
           endif
 
           if ((voltotal.le.zero).or.(mass_total.le.zero)) then
            print *,"voltotal or mass_total invalid"
            print *,"voltotal, mass_total, nmat ",voltotal,mass_total,nmat
            print *,"voldonate,volrecon ",voldonate,volrecon
            print *,"veldir ",veldir
            print *,"fablo ",fablo(1),fablo(2),fablo(SDIM)
            print *,"fabhi ",fabhi(1),fabhi(2),fabhi(SDIM)
            print *,"i,j,k ",i,j,k
            stop
           endif
           if ((voltotal_solid.lt.zero).or.(mass_total_solid.lt.zero)) then
            print *,"voltotal_solid or mass_total_solid invalid"
            stop
           endif
           if ((voltotal_fluid.lt.zero).or.(mass_total_fluid.lt.zero)) then
            print *,"voltotal_fluid or mass_total_fluid invalid"
            stop
           endif

          else
           print *,"check_donate invalid"
           stop
          endif

         enddo ! iside

        enddo
        enddo
        enddo ! i,j,k 

      enddo ! veldir

      return
      end subroutine FORT_BUILD_SEMIREFINEVOF

      ! unsplit_flag=
      ! 0=directionally split
      ! 1=unsplit everywhere
      ! 2=unsplit in incompressible zones
      ! 3=unsplit in fluid cells that neighbor a prescribed solid cell
      subroutine FORT_BUILD_MASK_UNSPLIT( &
       unsplit_flag, &
       make_interface_incomp, &
       xlo,dx, &
       maskunsplit, &
       DIMS(maskunsplit), &
       vofls0,DIMS(vofls0), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       nmat, &
       level,finest_level)
      use global_utility_module
      use MOF_routines_module
      use probf90_module

      IMPLICIT NONE

      INTEGER_T unsplit_flag
      INTEGER_T make_interface_incomp
      INTEGER_T level,finest_level
      INTEGER_T nmat
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T bfact
      INTEGER_T DIMDEC(maskunsplit)
      INTEGER_T DIMDEC(vofls0)
      REAL_T vofls0(DIMV(vofls0),2*nmat) 
      REAL_T maskunsplit(DIMV(maskunsplit)) 
      REAL_T xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k
      INTEGER_T im
      INTEGER_T im_opp
      INTEGER_T im_primary

      INTEGER_T growlo(3),growhi(3)
      INTEGER_T all_incomp,local_incomp
      INTEGER_T local_mask
      REAL_T cutoff,DXMAXLS
      REAL_T localLS(nmat)

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact too small"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level or finest_level invalid build mask_unsplit"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid build mask_unsplit "
       print *,"ngeom_recon= ",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid build mask_unsplit "
       print *,"ngeom_raw= ",ngeom_raw
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskunsplit),1,-1,219)
      call checkbound(fablo,fabhi,DIMS(vofls0),1,-1,219)

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        if (is_rigid(nmat,im).ne.1) then
         print *,"is_rigid(nmat,im).ne.1"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (fort_material_type(im).le.fort_max_num_eos)) then
        ! do nothing
       else 
        print *,"material_type invalid"
        stop
       endif

      enddo ! im=1..nmat  (checking parameters)

      all_incomp=1

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).ge.1).and. &
                (fort_material_type(im).le.fort_max_num_eos)) then
        all_incomp=0
       else
        print *,"fort_material_type invalid"
        stop
       endif

      enddo  ! im=1..nmat

      if ((make_interface_incomp.ne.0).and. &
          (make_interface_incomp.ne.1).and. & ! Drho/Dt=0 in incomp zone
          (make_interface_incomp.ne.2)) then ! rho_t + div(rho u)=0 " "
       print *,"make_interface_incomp invalid"
       stop
      endif

      if (all_incomp.eq.0) then
       if (unsplit_flag.eq.1) then ! unsplit everywhere
        print *,"unsplit is not conservative, cannot use in shock regions"
        stop
       else if (unsplit_flag.eq.0) then ! never unsplit
        ! do nothing
       else if (unsplit_flag.eq.2) then ! unsplit in incomp zones
        ! do nothing
        
        ! unsplit in fluid cells that neighbor a prescribed solid cell
       else if (unsplit_flag.eq.3) then  
        if (make_interface_incomp.eq.0) then
         print *,"unsplit is not conservative, cannot use in shock regions"
         stop
        else if ((make_interface_incomp.eq.1).or. &
                 (make_interface_incomp.eq.2)) then
         ! do nothing
        else
         print *,"make_interface_incomp invalid"
         stop
        endif
       else
        print *,"unsplit_flag invalid"
        stop
       endif
      else if (all_incomp.eq.1) then
       if ((unsplit_flag.eq.0).or. &
           (unsplit_flag.eq.1).or. &
           (unsplit_flag.eq.2).or. &
           (unsplit_flag.eq.3)) then
        ! do nothing
       else
        print *,"unsplit_flag invalid"
        stop
       endif
      else
       print *,"all_incomp invalid"
       stop
      endif

       ! the band thickness in FORT_ADVECTIVE_PRESSURE is 2 * DXMAXLS
      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=DXMAXLS
         
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       do im=1,nmat
        localLS(im)=vofls0(D_DECL(i,j,k),nmat+im)
       enddo
       call get_primary_material(localLS,nmat,im_primary)

       local_incomp=0
       if (all_incomp.eq.1) then
        local_incomp=1
       else if (all_incomp.eq.0) then
        if ((make_interface_incomp.eq.1).or. & !Drho/Dt=0 incomp zone
            (make_interface_incomp.eq.2)) then !rho_t+div(rho u) incomp zone
         do im=1,nmat
          if (abs(localLS(im)).le.cutoff) then ! cutoff=DXMAXLS
           do im_opp=1,nmat
            if (im.ne.im_opp) then
             if (abs(localLS(im_opp)).le.cutoff) then
              local_incomp=1
             else if (abs(localLS(im_opp)).gt.cutoff) then
              ! do nothing
             else
              print *,"localLS invalid"
              stop
             endif
            else if (im.eq.im_opp) then
             ! do nothing
            else
             print *,"im or im_opp invalid"
             stop
            endif
           enddo ! im_opp=1..nmat 
          else if (abs(localLS(im)).gt.cutoff) then ! cutoff=DXMAXLS
           ! do nothing
          else
           print *,"localLS invalid"
           stop
          endif
         enddo ! im=1..nmat
        else if (make_interface_incomp.eq.0) then
         ! do nothing
        else
         print *,"make_interface_incomp invalid"
         stop
        endif
       else
        print *,"all_incomp invalid"
        stop
       endif 

       local_mask=0

       if (unsplit_flag.eq.0) then ! split
        local_mask=0
       else if (unsplit_flag.eq.1) then ! unsplit everywhere
        local_mask=1
       else if (unsplit_flag.eq.2) then ! unsplit in incomp. zones
        if (local_incomp.eq.1) then
         local_mask=1
        else if (local_incomp.eq.0) then
         if ((fort_material_type(im_primary).eq.0).or. &
             (fort_material_type(im_primary).eq.999)) then
          local_mask=1
         else if ((fort_material_type(im_primary).gt.0).and. &
                  (fort_material_type(im_primary).le.fort_max_num_eos)) then
          ! do nothing 
         else 
          print *,"material_type invalid"
          stop
         endif
        else
         print *,"local_incomp invalid"
         stop
        endif
        ! unsplit in fluid cells that neighbor a prescribed solid cell
       else if (unsplit_flag.eq.3) then 
        do im=1,nmat
         if (is_rigid(nmat,im).eq.1) then
          if (is_prescribed(nmat,im).eq.1) then
           if (abs(localLS(im)).le.cutoff) then
            local_mask=1
           else if (abs(localLS(im)).ge.cutoff) then
            ! do nothing
           else
            print *,"localLS invalid"
            stop
           endif
          else if (is_prescribed(nmat,im).eq.0) then
           ! do nothing
          else
           print *,"is_prescribed(nmat,im) invalid"
           stop
          endif
         else if (is_rigid(nmat,im).eq.0) then
          ! do nothing
         else
          print *,"is_rigid(nmat,im) invalid"
          stop
         endif
        enddo ! im=1..nmat
       else
        print *,"unsplit_flag invalid"
        stop
       endif 

       maskunsplit(D_DECL(i,j,k))=local_incomp

      enddo
      enddo
      enddo ! i,j,k 

      return
      end subroutine FORT_BUILD_MASK_UNSPLIT


      subroutine FORT_BUILD_MODVISC( &
       ngrow_visc, &
       time, &
       problo,probhi, &
       visc_coef, &
       nten, &
       xlo,dx, &
       slope,DIMS(slope), &
       denstate,DIMS(denstate), &
       viscstate,DIMS(viscstate), &
       levelPC,DIMS(levelPC), &
       modvisc,DIMS(modvisc), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       nmat, &
       level,finest_level)
      use global_utility_module
      use probf90_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: ngrow_visc
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: nten
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: problo(SDIM),probhi(SDIM)
      REAL_T, intent(in) :: visc_coef
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(slope)
      INTEGER_T, intent(in) :: DIMDEC(denstate)
      INTEGER_T, intent(in) :: DIMDEC(viscstate)
      INTEGER_T, intent(in) :: DIMDEC(levelPC)
      INTEGER_T, intent(in) :: DIMDEC(modvisc)
      REAL_T, intent(in) :: slope(DIMV(slope),nmat*ngeom_recon) 
      REAL_T, intent(in) :: denstate(DIMV(denstate),nmat*num_state_material) 
      REAL_T, intent(in) :: viscstate(DIMV(viscstate),nmat) 
      REAL_T, intent(in) :: levelPC(DIMV(levelPC),nmat)
      REAL_T, intent(out) :: modvisc(DIMV(modvisc),nmat)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k

      INTEGER_T im

      INTEGER_T nten_test
      REAL_T localvisc(nmat)
      INTEGER_T check_accept

      INTEGER_T igridlo(3),igridhi(3)

      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T mu

      nhalf=1

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level or finest_level invalid build modvisc"
       stop
      endif

 
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid build modvisc nten nten_test ",nten,nten_test
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid build modvisc "
       print *,"ngeom_recon= ",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid build modvisc "
       print *,"ngeom_raw= ",ngeom_raw
       stop
      endif
      if (ngrow_visc.lt.1) then
       print *,"ngrow_visc out of range"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(slope),ngrow_visc,-1,219)
      call checkbound(fablo,fabhi, &
       DIMS(denstate), &
       ngrow_visc,-1,223)
      call checkbound(fablo,fabhi, &
       DIMS(viscstate), &
       ngrow_visc,-1,224)
      call checkbound(fablo,fabhi, &
       DIMS(levelPC), &
       ngrow_visc+1,-1,229) 
      call checkbound(fablo,fabhi,DIMS(modvisc),ngrow_visc,-1,227)

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).gt.0).and. &
                (fort_material_type(im).le.fort_max_num_eos)) then
        ! do nothing
       else 
        print *,"material_type invalid"
        stop
       endif

       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"density must be positive build modvisc"
        print *,"im,denconst ",im,fort_denconst(im)
        stop
       endif
       if (fort_tempconst(im).gt.zero) then
        ! do nothing
       else
        print *,"build modvisc:temperature must be positive"
        print *,"im,fort_tempconst : ",im,fort_tempconst(im)
        stop
       endif
       if (fort_energyconst(im).gt.zero) then
        ! do nothing
       else
        print *,"energy must be positive in FORT_BUILD_MODVISC"
        print *,"im= ",im
        print *,"fort_energyconst(im)= ",fort_energyconst(im)
        stop
       endif
       mu=get_user_viscconst(im,fort_denconst(im),fort_tempconst(im))
       if (mu.ge.zero) then
        ! do nothing
       else
        print *,"viscosity cannot be negative"
        stop
       endif

      enddo ! im  (checking parameters)

      call growntilebox(tilelo,tilehi,fablo,fabhi,igridlo,igridhi, &
       ngrow_visc) 

      do i=igridlo(1),igridhi(1)
      do j=igridlo(2),igridhi(2)
      do k=igridlo(3),igridhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       check_accept=1
       if (levelrz.eq.0) then
        ! do nothing
       else if ((levelrz.eq.1).or. &
                (levelrz.eq.3)) then
        if (xsten(0,1).lt.VOFTOL*dx(1)) then
         check_accept=0
        endif
       else
        print *,"levelrz invalid build mod visc"
        stop
       endif

       if (check_accept.eq.0) then

        do im=1,nmat
         modvisc(D_DECL(i,j,k),im)=zero
        enddo

       else if (check_accept.eq.1) then

        do im=1,nmat
         localvisc(im)=viscstate(D_DECL(i,j,k),im)
         if (localvisc(im).lt.zero) then
          print *,"viscstate gone negative"
          stop
         endif
         modvisc(D_DECL(i,j,k),im)=localvisc(im)
        enddo

       else
        print *,"check_accept invalid"
        stop
       endif

      enddo
      enddo
      enddo  ! i,j,k (modvisc)

 
      return
      end subroutine FORT_BUILD_MODVISC

       ! operation_flag=0 (right hand side for solver)
       ! operation_flag=1 (divergence)
       ! operation_flag=2 (mac -> cell velocity in solver or MAC_TO_CELL)
       ! operation_flag=3 (cell pressure gradient, cell energy)
       ! operation_flag=4 (gravity and surface tension force at cell)
       ! in the face_gradients routine, operation_flag=5 (interp grad U^T)
       ! operation_flag=6 (advection)
      subroutine FORT_MAC_TO_CELL( &
       nsolveMM_FACE, &
       num_materials_face, &
       ns_time_order, &
       divu_outer_sweeps, &
       num_divu_outer_sweeps, &
       operation_flag, &
       energyflag, &
       temperature_primitive_variable, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       prescribed_solid_scale, &
       added_weight, &
       nten, &
       level, &
       finest_level, &
       face_flag, &
       make_interface_incomp, &
       local_solvability_projection, &
       project_option, &
       enable_spectral, &
       fluxvel_index, &
       fluxden_index, &
       facevel_index, &
       facecut_index, &
       icefacecut_index, &
       curv_index, &
       conservative_tension_force, &
       conservative_div_uu, &
       pforce_index, &
       faceden_index, &
       icemask_index, &
       massface_index, &
       vofface_index, &
       ncphys, &
       velbc_in, &
       presbc_in, &
       cur_time, &
       slab_step, &
       dt, &
       xlo,dx, &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       xp,DIMS(xp), &
       yp,DIMS(yp), &
       zp,DIMS(zp), &
       xvel,DIMS(xvel), &
       yvel,DIMS(yvel), &
       zvel,DIMS(zvel), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       ax,DIMS(ax), &
       ay,DIMS(ay), &
       az,DIMS(az), &
       vol,DIMS(vol), &
       rhs,DIMS(rhs), & ! destination for div(up) if SEM.
       veldest,DIMS(veldest), &
       dendest,DIMS(dendest), &
       mask,DIMS(mask), & ! 1=fine/fine  0=coarse/fine
       maskcoef,DIMS(maskcoef), & ! 1=not covered  0=covered
       maskSEM,DIMS(maskSEM), &
       levelPC,DIMS(levelPC), &
       sol,DIMS(sol), &
       cterm,DIMS(cterm), &
       pold,DIMS(pold), &
       denold,DIMS(denold), &
       ustar,DIMS(ustar), &
       recon,DIMS(recon), &
       mdotcell,DIMS(mdotcell), & ! holds velocity if operation_flag==6
       maskdivres,DIMS(maskdivres), &
       maskres,DIMS(maskres), &
       SDC_outer_sweeps, &
       homflag, &
       use_VOF_weight, &
       nsolve, &
       ncomp_denold, &
       ncomp_veldest, &
       ncomp_dendest, &
       SEM_advection_algorithm)
       use probf90_module
       use global_utility_module
       use MOF_routines_module
       use CISL_SANITY_MODULE
       IMPLICIT NONE

      INTEGER_T, intent(in) :: ncomp_denold
      INTEGER_T, intent(in) :: ncomp_veldest
      INTEGER_T, intent(in) :: ncomp_dendest
      INTEGER_T, intent(in) :: num_materials_face
      INTEGER_T, intent(in) :: nsolveMM_FACE
      INTEGER_T, intent(in) :: ns_time_order
      INTEGER_T, intent(in) :: divu_outer_sweeps
      INTEGER_T, intent(in) :: num_divu_outer_sweeps
      INTEGER_T, intent(in) :: SEM_advection_algorithm
      INTEGER_T :: high_order_time_advection
      INTEGER_T, intent(in) :: operation_flag
      INTEGER_T, intent(in) :: slab_step
      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: SDC_outer_sweeps 
      INTEGER_T, intent(in) :: face_flag 
      INTEGER_T, intent(in) :: make_interface_incomp
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: energyflag 
      INTEGER_T, intent(in) :: temperature_primitive_variable(nmat)
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)
      REAL_T, intent(in) :: prescribed_solid_scale(nmat)
      REAL_T, intent(in) :: added_weight(nmat)
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: homflag
      INTEGER_T, intent(in) :: use_VOF_weight
      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: local_solvability_projection
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: fluxvel_index
      INTEGER_T, intent(in) :: fluxden_index
      INTEGER_T, intent(in) :: facevel_index
      INTEGER_T, intent(in) :: facecut_index
      INTEGER_T, intent(in) :: icefacecut_index
      INTEGER_T, intent(in) :: curv_index
      INTEGER_T, intent(in) :: conservative_tension_force
      INTEGER_T, intent(in) :: conservative_div_uu
      INTEGER_T, intent(in) :: pforce_index
      INTEGER_T, intent(in) :: faceden_index
      INTEGER_T, intent(in) :: icemask_index
      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: velbc_in(SDIM,2,SDIM*num_materials_vel)
      INTEGER_T, intent(in) :: presbc_in(SDIM,2,num_materials_face)
      REAL_T, intent(in) :: cur_time,dt
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(xp)
      INTEGER_T, intent(in) :: DIMDEC(yp)
      INTEGER_T, intent(in) :: DIMDEC(zp)
      INTEGER_T, intent(in) :: DIMDEC(xvel)
      INTEGER_T, intent(in) :: DIMDEC(yvel)
      INTEGER_T, intent(in) :: DIMDEC(zvel)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(ax)
      INTEGER_T, intent(in) :: DIMDEC(ay)
      INTEGER_T, intent(in) :: DIMDEC(az)
      INTEGER_T, intent(in) :: DIMDEC(vol)
      INTEGER_T, intent(in) :: DIMDEC(rhs)
      INTEGER_T, intent(in) :: DIMDEC(veldest)
      INTEGER_T, intent(in) :: DIMDEC(dendest)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(maskcoef)
      INTEGER_T, intent(in) :: DIMDEC(maskSEM)
      INTEGER_T, intent(in) :: DIMDEC(levelPC)
      INTEGER_T, intent(in) :: DIMDEC(sol)
      INTEGER_T, intent(in) :: DIMDEC(cterm)
      INTEGER_T, intent(in) :: DIMDEC(pold)
      INTEGER_T, intent(in) :: DIMDEC(denold)
      INTEGER_T, intent(in) :: DIMDEC(ustar)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(mdotcell)
      INTEGER_T, intent(in) :: DIMDEC(maskdivres)
      INTEGER_T, intent(in) :: DIMDEC(maskres)

      REAL_T, intent(in) ::  xp(DIMV(xp),2+nsolveMM_FACE)
      REAL_T, intent(in) ::  yp(DIMV(yp),2+nsolveMM_FACE)
      REAL_T, intent(in) ::  zp(DIMV(zp),2+nsolveMM_FACE)

      REAL_T, intent(in) ::  xvel(DIMV(xvel),nsolveMM_FACE)
      REAL_T, intent(in) ::  yvel(DIMV(yvel),nsolveMM_FACE)
      REAL_T, intent(in) ::  zvel(DIMV(zvel),nsolveMM_FACE)

      REAL_T, intent(in) ::  xface(DIMV(xface),ncphys)
      REAL_T, intent(in) ::  yface(DIMV(yface),ncphys)
      REAL_T, intent(in) ::  zface(DIMV(zface),ncphys)

      REAL_T, intent(in) ::  ax(DIMV(ax))
      REAL_T, intent(in) ::  ay(DIMV(ay))
      REAL_T, intent(in) ::  az(DIMV(az))

      REAL_T, intent(in) :: vol(DIMV(vol))
      REAL_T, intent(inout) :: rhs(DIMV(rhs),nsolve*num_materials_face)
      REAL_T, intent(inout) :: veldest(DIMV(veldest),ncomp_veldest)
      REAL_T, intent(inout) :: dendest(DIMV(dendest),ncomp_dendest)
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: maskcoef(DIMV(maskcoef))
      REAL_T, intent(in) :: maskSEM(DIMV(maskSEM))
      REAL_T, intent(in) :: levelPC(DIMV(levelPC),nmat*(SDIM+1))
      REAL_T, intent(in) :: sol(DIMV(sol),SDIM*nparts_def)
      REAL_T, intent(inout) :: cterm(DIMV(cterm),nsolve*num_materials_face)
      REAL_T, intent(in) :: pold(DIMV(pold),nsolve*num_materials_face)
      REAL_T, intent(in) :: denold(DIMV(denold),ncomp_denold)
      REAL_T, intent(inout) :: ustar(DIMV(ustar),SDIM*num_materials_vel) 
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T, intent(in) :: mdotcell(DIMV(mdotcell),nsolve*num_materials_face)
      REAL_T, intent(in) :: maskdivres(DIMV(maskdivres))
      REAL_T, intent(in) :: maskres(DIMV(maskres))

      REAL_T DXMAXLS,cutoff
      INTEGER_T all_incomp
      INTEGER_T local_primitive,local_incomp
      REAL_T Eforce,Eforce_primitive,cell_pressure
      REAL_T KE_diff

      INTEGER_T cell_is_ice
      INTEGER_T cell_is_FSI_rigid

      INTEGER_T i,j,k
      INTEGER_T dir,dir2,side
      INTEGER_T veldir
      INTEGER_T veldir_left
      INTEGER_T veldir_right
      INTEGER_T im,im_opp
      INTEGER_T vofcomp
      INTEGER_T sidecomp,ibase
      INTEGER_T ii,jj,kk
      INTEGER_T iface,jface,kface
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T nten_test
      INTEGER_T imattype
      REAL_T AXL,AXR
      REAL_T AYL,AYR
      REAL_T AZL,AZR
      REAL_T pgrad
      REAL_T VOLTERM,hx,RR
      REAL_T dencell,dencellgrav,rho
      REAL_T TEMPERATURE,internal_e
      REAL_T NEW_TEMPERATURE
      REAL_T CC,CC_DUAL,MSKDV,MSKRES,MDOT,divu,dp
      REAL_T local_rhs
      REAL_T local_POLD
      REAL_T local_POLD_DUAL

      REAL_T DIAG_REGULARIZE
      REAL_T uface(2,num_materials_face)
      REAL_T aface(2)
      REAL_T pfacegrav(2)
      REAL_T pfacetenleft(2)
      REAL_T pfacetenright(2)
      REAL_T pfaceten(2)
      REAL_T pres_face(2)

       ! 0=no gp or div(up)
       ! 1=no gp but div(up)
       ! 2=no div(up) but gp
       ! 3=both
      INTEGER_T use_face_pres_cen
       ! REALs since they are stored with other vars in "xp","yp","zp"
      REAL_T use_face_pres(2)  ! faces that are on either side of a cell.
      REAL_T use_face_pres_combine
      REAL_T coarse_fine_face(2)
      REAL_T ASIDE(2,ncphys)
      REAL_T mass_side(2)
      REAL_T masscell
      REAL_T fluid_velocity
      INTEGER_T local_maskSEM
      INTEGER_T stripstat
      INTEGER_T elemlo(3),elemhi(3)
      INTEGER_T ielem,jelem,kelem
      INTEGER_T scomp,scomp_bc,dcomp,ncomp ! in: mac_to_cell
      INTEGER_T ncomp_xvel
      INTEGER_T ncomp_cterm
      REAL_T vfrac(nmat)
      REAL_T LStest(nmat)
      INTEGER_T im_vel
      INTEGER_T im_vel_left
      INTEGER_T im_vel_right
      INTEGER_T velcomp
      INTEGER_T nsolveMM_FACE_test
      INTEGER_T partid,nparts_temp,im_solid
      INTEGER_T cell_velocity_override

      REAL_T local_div_val

      REAL_T, dimension(:,:), allocatable :: comparepface
      REAL_T, dimension(:,:), allocatable :: comparevelface
      REAL_T, dimension(:,:), allocatable :: comparestate

      nhalf=3
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid MAC_TO_CELL"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid MAC_TO_CELL"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((ns_time_order.ge.1).and.(ns_time_order.le.32)) then
       ! do nothing
      else
       print *,"ns_time_order invalid"
       stop
      endif
      if (num_divu_outer_sweeps.lt.1) then
       print *,"num_divu_outer_sweeps invalid"
       stop
      endif
      if ((divu_outer_sweeps.lt.0).or. &
          (divu_outer_sweeps.ge.num_divu_outer_sweeps)) then
       print *,"divu_outer_sweeps invalid"
       stop
      endif

      if ((SEM_advection_algorithm.eq.0).or. &
          (SEM_advection_algorithm.eq.1)) then
       ! do nothing
      else
       print *,"SEM_advection_algorithm invalid"
       stop
      endif
       
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid mac to cell"
       stop
      endif

      if ((SDC_outer_sweeps.ge.0).and. &
          (SDC_outer_sweeps.lt.ns_time_order)) then
       ! do nothing
      else
       print *,"SDC_outer_sweeps invalid in mac to cell"
       print *,"SDC_outer_sweeps= ",SDC_outer_sweeps
       stop
      endif

      if ((slab_step.lt.-1).or. &
          (slab_step.gt.bfact_time_order)) then
       print *,"slab_step invalid mac to cell"
       stop
      endif

      if (1.eq.0) then
       print *,"in: mac_to_cell: operation_flag=",operation_flag
      endif

      if ((use_VOF_weight.eq.0).or.(use_VOF_weight.eq.1)) then
       ! do nothing
      else
       print *,"use_VOF_weight invalid"
       stop
      endif 

      if ((enable_spectral.ge.0).and. &
          (enable_spectral.le.3)) then
       ! do nothing
      else
       print *,"enable_spectral invalid mac_to_cell"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten.ne.nten_test) then
       print *,"nten invalid mac_to_cell nten nten_test ",nten,nten_test
       stop
      endif
      if ((local_solvability_projection.ne.0).and. &
          (local_solvability_projection.ne.1)) then
       print *,"local_solvability_projection invalid"
       stop
      endif
      if ((make_interface_incomp.ne.0).and. &
          (make_interface_incomp.ne.1).and. &
          (make_interface_incomp.ne.2)) then
       print *,"make_interface_incomp invalid"
       stop
      endif

      do im=1,nmat
       imattype=fort_material_type(im)
       if (imattype.eq.999) then
        ! do nothing
       else if ((imattype.ge.0).and. &
                (imattype.le.fort_max_num_eos)) then
        ! do nothing
       else
        print *,"imattype invalid"
        stop
       endif
       if ((prescribed_solid_scale(im).ge.zero).and. &
           (prescribed_solid_scale(im).lt.half)) then
        ! do nothing
       else
        print *,"prescribed_solid_scale invalid"
        stop
       endif
       if (added_weight(im).gt.zero) then
        ! do nothing
       else
        print *,"added_weight invalid"
        stop
       endif
       if ((temperature_primitive_variable(im).ne.0).and. &
           (temperature_primitive_variable(im).ne.1)) then
        print *,"temperature_primitive_variable invalid"
        stop
       endif
      enddo ! im=1..nmat
 
      ! indexes start at 0
      if ((curv_index.ne.0).or. &
          (pforce_index.ne.1).or. &
          (facecut_index.ne.3).or. &
          (icefacecut_index.ne.4).or. &
          (icemask_index.ne.5).or. &
          (facevel_index.ne.8).or. &
          (faceden_index.ne.2).or. &
          (vofface_index.ne.massface_index+2*nmat).or. &
          (fluxvel_index.ne.0).or. &
          (fluxden_index.ne.SDIM)) then
       print *,"face_index bust 2"
       stop
      endif

      if ((conservative_tension_force.ne.0).and. &
          (conservative_tension_force.ne.1)) then
       print *,"conservative_tension_force invalid"
       stop
      endif
      if ((conservative_div_uu.eq.0).or. &
          (conservative_div_uu.eq.1).or. &
          (conservative_div_uu.eq.2)) then
       ! do nothing
      else
       print *,"conservative_div_uu invalid"
       stop
      endif
   
      ! mac -> cell in solver (apply_cell_pressure_gradient) or VELMAC_TO_CELL
      if (operation_flag.eq.2) then 

       if (ncomp_veldest.ge. &
           num_materials_vel*SDIM+num_state_material*nmat) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.ge.num_state_material*nmat) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif
       if (nsolve.ne.1) then
        print *,"nsolve invalid 2"
        stop
       endif
       if (num_materials_face.ne.1) then
        print *,"num_materials_face invalid"
        stop
       endif
       nsolveMM_FACE_test=nsolve*num_materials_face

       if ((ncomp_denold.eq.nsolve*num_materials_face).or. &
           (ncomp_denold.eq.1)) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

      else if (operation_flag.eq.3) then ! cell grad p, cell energy
 
       if (ncomp_veldest.ge. &
           num_materials_vel*SDIM+num_state_material*nmat) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.ge.num_state_material*nmat) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif
       if (nsolve.ne.1) then
        print *,"nsolve invalid 2"
        stop
       endif
       if (num_materials_face.ne.1) then
        print *,"num_materials_face invalid"
        stop
       endif
       nsolveMM_FACE_test=nsolve*num_materials_face

       if (ncomp_denold.eq.nsolveMM_FACE_test) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

      else if (operation_flag.eq.0) then ! rhs of solver

       if (ncomp_veldest.eq.nsolve) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.eq.nsolve) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif
       if ((nsolve.ne.1).and. &
           (nsolve.ne.SDIM)) then
        print *,"nsolve invalid 2"
        stop
       endif
       nsolveMM_FACE_test=nsolve*num_materials_face
       if (num_materials_face.eq.1) then
        ! do nothing
       else if (num_materials_face.eq.nmat) then
        nsolveMM_FACE_test=nsolveMM_FACE_test*2
       else
        print *,"num_materials_face invalid"
        stop
       endif

       if (ncomp_denold.eq.nsolve*num_materials_face) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

      else if (operation_flag.eq.1) then ! divergence

       if (ncomp_veldest.eq.num_materials_vel) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.eq.num_materials_vel) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncphys.eq.vofface_index+2*nmat) then
        ! do nothing
       else
        print *,"ncphys invalid"
        stop
       endif
       nsolveMM_FACE_test=nsolve*num_materials_face
       if (num_materials_face.eq.1) then
        ! do nothing
       else
        print *,"num_materials_face invalid"
        stop
       endif
       if (nsolve.eq.1) then
        ! do nothing
       else
        print *,"nsolve invalid"
        stop
       endif
       if (ncomp_denold.eq.1) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

      else if (operation_flag.eq.4) then ! gravity and surface tension

       if (ncomp_veldest.eq.SDIM) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.eq.SDIM) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif
       if (nsolve.eq.1) then
        ! do nothing
       else
        print *,"nsolve invalid"
        stop
       endif
       if (num_materials_face.eq.1) then
        ! do nothing
       else
        print *,"num_materials_face invalid"
        stop
       endif
       nsolveMM_FACE_test=num_materials_face

       if (ncomp_denold.eq.1) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

      else if (operation_flag.eq.5) then

       print *,"grad U^MAC -> grad U^CELL called from FACE_GRADIENTS"
       stop

      else if (operation_flag.eq.6) then ! advection

       if (ncomp_veldest.ge. &
           num_materials_vel*SDIM+num_state_material*nmat) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.ge.num_state_material*nmat) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif

       if ((nsolve.ne.SDIM+num_state_base).or. &
           (ncphys.ne.SDIM+num_state_base)) then
        print *,"nsolve or ncphys invalid"
        stop
       endif
       if (num_materials_face.ne.1) then
        print *,"num_materials_face invalid"
        stop
       endif
       nsolveMM_FACE_test=num_materials_face

       if (ncomp_denold.eq.nmat*num_state_material) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

      else
       print *,"operation_flag invalid6"
       stop
      endif

      if (nsolveMM_FACE_test.ne.nsolveMM_FACE) then
       print *,"nsolveMM_FACE invalid"
       stop
      endif

      if (operation_flag.eq.0) then  ! rhs for solver

       if (energyflag.eq.0) then
        ! do nothing
       else
        print *,"energyflag invalid"
        stop
       endif
       if ((homflag.ge.0).and.(homflag.le.4)) then
        ! do nothing
       else
        print *,"homflag invalid in mac to cell homflag=",homflag
        stop
       endif

      else if (operation_flag.eq.1) then ! divergence

       if (energyflag.eq.0) then
        ! do nothing
       else
        print *,"energyflag invalid"
        stop
       endif
       if (homflag.eq.0) then
        ! do nothing
       else
        print *,"homflag invalid"
        stop
       endif
       if (nsolve.eq.1) then
        ! do nothing
       else
        print *,"nsolve invalid 3"
        stop
       endif

       ! umac->ucell in solver or VELMAC_TO_CELL
      else if (operation_flag.eq.2) then 

       if (homflag.ne.0) then
        print *,"homflag invalid"
        stop
       endif
       if ((energyflag.ne.0).and. &
           (energyflag.ne.1)) then
        print *,"energyflag invalid"
        stop
       endif
       if (nsolve.ne.1) then
        print *,"nsolve invalid5"
        stop
       endif

      else if (operation_flag.eq.3) then ! gradp^cell, div(up)
       if (homflag.ne.0) then
        print *,"homflag invalid"
        stop
       endif
       if ((energyflag.ne.0).and. &
           (energyflag.ne.1).and. &
           (energyflag.ne.2)) then
        print *,"energyflag invalid"
        stop
       endif
       if (nsolve.ne.1) then
        print *,"nsolve invalid6"
        stop
       endif

      else if (operation_flag.eq.4) then ! gravity and surface tension at cell

       if (energyflag.ne.0) then
        print *,"energyflag invalid"
        stop
       endif
       if (homflag.ne.0) then
        print *,"homflag invalid"
        stop
       endif

      else if (operation_flag.eq.6) then ! advection

        ! "source_term" 
       if ((homflag.ne.0).and.(homflag.ne.1)) then
        print *,"homflag invalid in mac to cell homflag=",homflag
        stop
       endif
        ! "advect_iter"
       if ((energyflag.ne.0).and.(energyflag.ne.1)) then
        print *,"energyflag invalid"
        stop
       endif
      else
       print *,"operation_flag invalid7"
       stop
      endif
 
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((project_option.eq.0).or. &
          (project_option.eq.1).or. &
          (project_option.eq.10).or. &
          (project_option.eq.11).or. & ! FSI_material_exists 2nd project
          (project_option.eq.13).or. & ! FSI_material_exists 1st project
          (project_option.eq.12).or. & ! pressure extension
          (project_option.eq.2).or. &
          (project_option.ge.3).or. &
          ((project_option.ge.100).and. &
           (project_option.lt.100+num_species_var))) then
       ! do nothing
      else
       print *,"project_option invalid"
       stop
      endif
      do im=1,nmat
       if (fort_denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif
      enddo

      call checkbound(fablo,fabhi,DIMS(xp),0,0,33)
      call checkbound(fablo,fabhi,DIMS(yp),0,1,33)
      call checkbound(fablo,fabhi,DIMS(zp),0,SDIM-1,33)

      call checkbound(fablo,fabhi,DIMS(xvel),0,0,33)
      call checkbound(fablo,fabhi,DIMS(yvel),0,1,33)
      call checkbound(fablo,fabhi,DIMS(zvel),0,SDIM-1,33)

      call checkbound(fablo,fabhi,DIMS(xface),0,0,33)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,33)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,33)

      call checkbound(fablo,fabhi,DIMS(ax),0,0,33)
      call checkbound(fablo,fabhi,DIMS(ay),0,1,33)
      call checkbound(fablo,fabhi,DIMS(az),0,SDIM-1,33)

      call checkbound(fablo,fabhi,DIMS(vol),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(rhs),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(veldest),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(dendest),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)
      call checkbound(fablo,fabhi,DIMS(mask),1,-1,133)
      call checkbound(fablo,fabhi,DIMS(maskcoef),1,-1,134)
      call checkbound(fablo,fabhi,DIMS(levelPC),1,-1,135)
      call checkbound(fablo,fabhi,DIMS(sol),1,-1,136)
      call checkbound(fablo,fabhi,DIMS(cterm),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(pold),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(denold),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(ustar),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(recon),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(mdotcell),0,-1,33)
      call checkbound(fablo,fabhi, &
       DIMS(maskdivres), &
       0,-1,137)
      call checkbound(fablo,fabhi, &
       DIMS(maskres), &
       0,-1,138)

      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=DXMAXLS

      all_incomp=1

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).ge.1).and. &
                (fort_material_type(im).le.fort_max_num_eos)) then
        all_incomp=0
       else
        print *,"fort_material_type invalid"
        stop
       endif

      enddo  ! im=1..nmat

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       AXL=ax(D_DECL(i,j,k))
       AXR=ax(D_DECL(i+1,j,k))
       AYL=ay(D_DECL(i,j,k))
       AYR=ay(D_DECL(i,j+1,k))
       AZL=az(D_DECL(i,j,k))
       AZR=az(D_DECL(i,j,k+1))
       VOLTERM=vol(D_DECL(i,j,k))
       if (VOLTERM.gt.zero) then
        ! do nothing
       else
        print *,"VOLTERM invalid"
        stop
       endif

       if (operation_flag.eq.1) then ! DIV

        im_vel=1
        im_vel_left=im_vel
        im_vel_right=im_vel

        divu= &
         AXR*xvel(D_DECL(i+1,j,k),im_vel_left)-  &
         AXL*xvel(D_DECL(i,j,k),im_vel_right)+ &
         AYR*yvel(D_DECL(i,j+1,k),im_vel_left)-  &
         AYL*yvel(D_DECL(i,j,k),im_vel_right)
        if (SDIM.eq.3) then
         divu=divu+ &
          AZR*zvel(D_DECL(i,j,k+1),im_vel_left)-  &
          AZL*zvel(D_DECL(i,j,k),im_vel_right)
        endif
        divu=divu/VOLTERM
        rhs(D_DECL(i,j,k),im_vel)=divu

       else if (operation_flag.eq.0) then ! RHS

        ! (cterm)*p-vol grad dot grad p/rho=-vol div u/dt + mdot +
        !    cterm * p^adv 
        ! cterm=vol/(rho c^2 dt*dt)

        if (maskcoef(D_DECL(i,j,k)).eq.one) then ! not covered

         do veldir=1,nsolve*num_materials_face
    
          CC=cterm(D_DECL(i,j,k),veldir)
          CC_DUAL=veldest(D_DECL(i,j,k),veldir)
          if (CC_DUAL.eq.CC) then
           ! do nothing
          else
           print *,"CC_DUAL invalid"
           stop
          endif
          MSKDV=maskdivres(D_DECL(i,j,k))
          MSKRES=maskres(D_DECL(i,j,k))
          MDOT=mdotcell(D_DECL(i,j,k),veldir)

          DIAG_REGULARIZE=denold(D_DECL(i,j,k),veldir)

          if (DIAG_REGULARIZE.gt.zero) then
           ! check nothing
          else
           print *,"DIAG_REGULARIZE invalid"
           stop
          endif

          if ((project_option.eq.0).or. &
              (project_option.eq.1).or. &
              (project_option.eq.10).or. &
              (project_option.eq.11).or. & ! FSI_material_exists 2nd project
              (project_option.eq.13).or. & ! FSI_material_exists 1st project
              (project_option.eq.12)) then ! pressure extension

           if (MDOT.eq.zero) then
            ! check nothing
           else if (MDOT.ne.zero) then
            if (MSKRES.eq.zero) then
             print *,"cannot have MDOT<>0 and MSKRES==0"
             stop
            else if (MSKRES.ne.zero) then
             ! do nothing
            else
             print *,"MSKRES bust"
             stop
            endif
           else
            print *,"MDOT bust"
            stop
           endif
           if (local_solvability_projection.eq.1) then
            if (CC.eq.zero) then
             ! do nothing
            else
             print *,"CC invalid"
             stop
            endif
            if (CC_DUAL.eq.zero) then
             ! do nothing
            else
             print *,"CC_DUAL invalid"
             stop
            endif
           else if (local_solvability_projection.eq.0) then
            if ((CC.ge.zero).and.(CC_DUAL.ge.zero)) then
             ! do nothing
            else
             print *,"CC or CC_DUAL invalid"
             stop
            endif
           else
            print *,"local_solvability_projection invalid"
            stop
           endif
          else if ((project_option.eq.2).or. & ! thermal diffusion
                   (project_option.eq.3).or. & ! viscosity
                   ((project_option.ge.100).and. &
                    (project_option.lt.100+num_species_var))) then
           ! do nothing
          else
           print *,"project_option invalid"
           stop
          endif 

          local_POLD=pold(D_DECL(i,j,k),veldir)
          local_POLD_DUAL=dendest(D_DECL(i,j,k),veldir)
          if ((homflag.eq.0).or.(homflag.eq.1)) then
           if (local_POLD.eq.local_POLD_DUAL) then
            rhs(D_DECL(i,j,k),veldir)=local_POLD*CC
           else
            print *,"local_POLD invalid"
            stop
           endif
          else if (homflag.eq.2) then
           if (local_POLD.eq.local_POLD_DUAL) then
            rhs(D_DECL(i,j,k),veldir)=-local_POLD*CC
           else
            print *,"local_POLD invalid"
            stop
           endif
          else if (homflag.eq.3) then
           if ((local_POLD.eq.zero).and.(local_POLD_DUAL.eq.zero)) then
            rhs(D_DECL(i,j,k),veldir)=zero
           else
            print *,"local_POLD or local_POLD_DUAL invalid"
            stop
           endif
          else if (homflag.eq.4) then
           if (local_POLD.eq.local_POLD_DUAL) then
            rhs(D_DECL(i,j,k),veldir)=zero
           else
            print *,"local_POLD invalid"
            stop
           endif
          else
           print *,"homflag invalid"
           stop
          endif

          if (MSKDV.eq.zero) then
           divu=zero
          else if (MSKDV.gt.zero) then

           if (num_materials_face.eq.1) then
            veldir_left=veldir
            veldir_right=veldir
           else if (num_materials_face.eq.nmat) then  
            if (nsolveMM_FACE/2.ne.nsolve*num_materials_face) then
             print *,"nsolveMM_FACE invalid"
             stop
            endif
            veldir_left=veldir
            veldir_right=veldir+nsolveMM_FACE/2
           else
            print *,"num_materials_face invalid"
            stop
           endif

            ! AXR,AXL,AYR,AYL,AZR,AZL are face areas.
           divu= &
            AXR*xvel(D_DECL(i+1,j,k),veldir_left)-  &
            AXL*xvel(D_DECL(i,j,k),veldir_right)+ &
            AYR*yvel(D_DECL(i,j+1,k),veldir_left)-  &
            AYL*yvel(D_DECL(i,j,k),veldir_right)
           if (SDIM.eq.3) then
            divu=divu+ &
             AZR*zvel(D_DECL(i,j,k+1),veldir_left)-  &
             AZL*zvel(D_DECL(i,j,k),veldir_right)
           endif
          else
           print *,"maskdivres invalid" 
           stop
          endif

           ! divu=-dt VOLTERM * div(k grad T)  project_option==2
           ! divu=-dt VOLTERM * visc_coef div(2 mu D) project_option==3
           ! use_dt=1 dir=-1
           ! use_HO=0
           ! constant_viscosity=1
          local_div_val=divu/VOLTERM
          call SEM_VISC_SANITY(110,dt,xsten,nhalf,local_div_val, &
            -1,veldir,1,0,project_option,bfact,enable_spectral,1)

          call SEM_VISC_SANITY_CC(1,dt,CC,MSKDV,MSKRES,MDOT, &
           VOLTERM,project_option,xsten,nhalf,veldir)

          local_rhs=rhs(D_DECL(i,j,k),veldir)

          if (homflag.eq.0) then
           rhs(D_DECL(i,j,k),veldir)=local_rhs-divu/dt+MDOT
          else if (homflag.eq.1) then
           rhs(D_DECL(i,j,k),veldir)=local_rhs+divu/dt
          else if (homflag.eq.2) then
           rhs(D_DECL(i,j,k),veldir)=local_rhs-divu/dt+MDOT
          else if (homflag.eq.3) then
           rhs(D_DECL(i,j,k),veldir)=-divu/dt+MDOT
           if (level.eq.finest_level) then
            if (divu.eq.zero) then
             ! do nothing
            else 
             print *,"divu invalid"
             stop
            endif
           else if ((level.ge.0).and.(level.lt.finest_level)) then
            ! do nothing
           else
            print *,"level invalid"
            stop
           endif
          else if (homflag.eq.4) then
           if (local_rhs.eq.zero) then
            rhs(D_DECL(i,j,k),veldir)=divu/VOLTERM
           else
            print *,"local_rhs invalid"
            stop
           endif
          else
           print *,"homflag invalid"
           stop
          endif

          if (MSKRES.eq.zero) then
           rhs(D_DECL(i,j,k),veldir)=zero
          else if (MSKRES.gt.zero) then
           ! do nothing
          else if (MSKRES.lt.zero) then
           print *,"maskres invalid"
           stop
          else
           print *,"MSKRES bust"
           stop
          endif 

         enddo ! veldir=1..nsolve*num_materials_face

        else if (maskcoef(D_DECL(i,j,k)).eq.zero) then
         ! do nothing (covered)
        else 
         print *,"mask invalid"
         stop
        endif

       ! mac -> cell in solver (apply_cell_pressure_gradient) or VELMAC_TO_CELL
       else if (operation_flag.eq.2) then

        do dir=0,SDIM-1 
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
          print *,"dir out of range in to_cell routine"
          stop
         endif

         ! side=1 left half of cell, side=2 right half of cell
         do side=1,2

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
  
          if (dir.eq.0) then
           do im=1,ncphys
            ASIDE(side,im)=xface(D_DECL(iface,jface,kface),im)
           enddo
          else if (dir.eq.1) then
           do im=1,ncphys
            ASIDE(side,im)=yface(D_DECL(iface,jface,kface),im)
           enddo
          else if ((dir.eq.2).and.(SDIM.eq.3)) then
           do im=1,ncphys
            ASIDE(side,im)=zface(D_DECL(iface,jface,kface),im)
           enddo
          else
           print *,"dir invalid mac to cell 3"
           stop
          endif

         enddo ! side

         ! side=1 left half of cell, side=2 right half of cell
         do side=1,2

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
 
          im_vel=1 

          if (dir.eq.0) then
           uface(side,im_vel)=xvel(D_DECL(iface,jface,kface),im_vel)
           if (SDIM.eq.2) then
            if (levelrz.eq.0) then
             ! do nothing
            else if (levelrz.eq.1) then
             if (iface.eq.0) then
              if (xsten(-1,1).lt.zero) then
               uface(side,im_vel)=zero
              endif
             else if (iface.gt.0) then
              if (xsten(0,1).gt.zero) then
               ! do nothing
              else
               print *,"xsten invalid"
               stop
              endif 
             else
              print *,"iface invalid"
              stop
             endif
            else if (levelrz.eq.3) then
             if (xsten(0,1).gt.zero) then
              ! do nothing
             else
              print *,"xsten invalid"
              stop
             endif 
            else
             print *,"levelrz invalid"
             stop
            endif
           else if (SDIM.eq.3) then
            ! do nothing
           else
            print *,"dimension bust"
            stop
           endif
          else if (dir.eq.1) then
           uface(side,im_vel)=yvel(D_DECL(iface,jface,kface),im_vel)
          else if ((dir.eq.2).and.(SDIM.eq.3)) then
           uface(side,im_vel)=zvel(D_DECL(iface,jface,kface),im_vel)
          else
           print *,"dir invalid mac to cell 3"
           stop
          endif

          mass_side(side)=zero
          do im=1,nmat 
           if (side.eq.1) then  ! left half of cell
            sidecomp=massface_index+2*(im-1)+2
           else if (side.eq.2) then ! right half of cell
            sidecomp=massface_index+2*(im-1)+1
           else
            print *,"side invalid"
            stop
           endif
           if (added_weight(im).gt.zero) then
            mass_side(side)=mass_side(side)+ &
             ASIDE(side,sidecomp)*added_weight(im) 
           else
            print *,"added_weight invalid"
            stop
           endif
          enddo ! im=1..nmat

          if (use_VOF_weight.eq.1) then
           ! do nothing
          else if (use_VOF_weight.eq.0) then
           mass_side(side)=one
          else
           print *,"use_VOF_weight invalid"
           stop
          endif

         enddo ! side=1..2

         masscell=mass_side(1)+mass_side(2)

         if ((mass_side(1).le.zero).or. &
             (mass_side(2).le.zero).or. &
             (masscell.le.zero)) then
          print *,"mass invalid"
          stop
         endif

         partid=-1
         nparts_temp=0
         im_solid=0
         do im=1,nmat
          LStest(im)=levelPC(D_DECL(i,j,k),im)
          if (is_lag_part(nmat,im).eq.1) then

           if (is_rigid(nmat,im).eq.1) then
            if (is_prescribed(nmat,im).eq.1) then
             if (im_solid.eq.0) then
              im_solid=im
              partid=nparts_temp
             else if ((im_solid.ge.1).and.(im_solid.le.nmat)) then
              if (LStest(im).gt.LStest(im_solid)) then
               im_solid=im
               partid=nparts_temp
              endif
             else
              print *,"im_solid invalid 5"
              stop
             endif
            else if (is_prescribed(nmat,im).eq.0) then
             ! do nothing
            else
             print *,"is_prescribed(nmat,im) invalid"
             stop
            endif
           else if (is_rigid(nmat,im).eq.0) then
            ! do nothing
           else
            print *,"is_rigid invalid"
            stop
           endif
           nparts_temp=nparts_temp+1

          else if (is_lag_part(nmat,im).eq.0) then
           if (is_rigid(nmat,im).eq.0) then
            ! do nothing
           else
            print *,"is_rigid invalid"
            stop
           endif
          else
           print *,"is_lag_part invalid"
           stop
          endif

         enddo ! im=1..nmat

         if (nparts_temp.ne.nparts) then
          print *,"nparts_temp invalid"
          stop
         endif
         if ((partid.ge.0).and.(partid.lt.nparts)) then
          if (im_solid_map(partid+1)+1.ne.im_solid) then
           print *,"im_solid invalid 6"
           stop
          endif
          if (LStest(im_solid).lt.zero) then
           partid=-1
          endif
         else if (partid.eq.-1) then
          ! do nothing
         else
          print *,"partid invalid"
          stop
         endif
 
         im_vel=1
         velcomp=dir+1

         cell_velocity_override=0
         if ((partid.ge.0).and.(partid.lt.nparts)) then
          if (prescribed_solid_scale(im_solid).eq.zero) then
           cell_velocity_override=1
          else if (prescribed_solid_scale(im_solid).gt.zero) then
           ! do nothing
          else 
           print *,"prescribed_solid_scale(im_solid) invalid"
           stop
          endif
         else if (partid.eq.-1) then
          ! do nothing
         else
          print *,"partid invalid"
          stop
         endif

         if (cell_velocity_override.eq.1) then
          veldest(D_DECL(i,j,k),velcomp)=sol(D_DECL(i,j,k),partid*SDIM+dir+1)
         else if (cell_velocity_override.eq.0) then
          fluid_velocity=(mass_side(1)*uface(1,im_vel)+ &
                          mass_side(2)*uface(2,im_vel))/masscell
          veldest(D_DECL(i,j,k),velcomp)=fluid_velocity
         else
          print *,"cell_velocity_override invalid"
          stop
         endif

         if (1.eq.0) then
          print *,"AFTER UPDATING CELL FROM FACE (velocity is scaled)"
          print *,"i,j,k,dir,velcomp ",i,j,k,dir,velcomp
          print *,"snew= ",veldest(D_DECL(i,j,k),velcomp)
         endif

        enddo  ! dir=0..sdim-1

        ! in: FORT_MAC_TO_CELL
       else if (operation_flag.eq.3) then ! (grad p)_CELL, div(up)

        use_face_pres_cen=3 ! both gp and div(up)

         ! note, in FORT_BUILD_CONSERVE, if 
         ! temperature_primitive_variable==1,
         ! then (rho T) is advected instead of (rho cv T + rho u dot u/2)
        do im=1,nmat
         imattype=fort_material_type(im)
         vofcomp=(im-1)*ngeom_recon+1
         vfrac(im)=recon(D_DECL(i,j,k),vofcomp)
         LStest(im)=levelPC(D_DECL(i,j,k),im)
         if (vfrac(im).ge.VOFTOL) then
          if (is_rigid(nmat,im).eq.1) then
           use_face_pres_cen=2  ! no div(up)
          else if (is_rigid(nmat,im).eq.0) then
           if (is_ice(nmat,im).eq.1) then
            use_face_pres_cen=2 ! no div(up)
           else if (is_FSI_rigid(nmat,im).eq.1) then
            use_face_pres_cen=2 ! no div(up)
           else if (imattype.eq.0) then
            use_face_pres_cen=2 ! no div(up)
           else if ((imattype.ge.1).and. &
                    (imattype.le.fort_max_num_eos)) then
            ! do nothing
           else
            print *,"imattype or FSI_flag invalid"
            stop
           endif
          else
           print *,"is_rigid(nmat,im) invalid"
           stop
          endif 
         endif ! vfrac(im).ge.voftol
        enddo ! im=1..nmat

        local_primitive=0
        local_incomp=0
        if (all_incomp.eq.1) then
         local_incomp=1
         local_primitive=1
        else if (all_incomp.eq.0) then
         if ((make_interface_incomp.eq.1).or. &
             (make_interface_incomp.eq.2)) then
          do im=1,nmat
           if (abs(LStest(im)).le.cutoff) then ! cutoff=DXMAXLS
            do im_opp=1,nmat
             if (im.ne.im_opp) then
              if (abs(LStest(im_opp)).le.cutoff) then ! cutoff=DXMAXLS
               local_incomp=1
              else if (abs(LStest(im_opp)).gt.cutoff) then
               ! do nothing
              else
               print *,"localLS invalid"
               stop
              endif
             else if (im.eq.im_opp) then
              ! do nothing
             else
              print *,"im or im_opp invalid"
              stop
             endif
            enddo ! im_opp=1..nmat 
           else if (abs(LStest(im)).gt.cutoff) then ! cutoff=DXMAXLS
            ! do nothing
           else
            print *,"LStest invalid"
            stop
           endif
          enddo ! im=1..nmat
         else if (make_interface_incomp.eq.0) then
          ! do nothing
         else
          print *,"make_interface_incomp invalid"
          stop
         endif
        else
         print *,"all_incomp invalid"
         stop
        endif 
        do im=1,nmat
         if (LStest(im).ge.zero) then
          if (temperature_primitive_variable(im).eq.1) then
           local_primitive=1
          else if (temperature_primitive_variable(im).eq.0) then
           ! do nothing
          else
           print *,"temperature_primitive_variable(im) invalid"
           stop
          endif
         else if (LStest(im).lt.zero) then
          ! do nothing
         else
          print *,"LStest(im) invalid"
          stop
         endif
        enddo ! im=1..nmat
     

        im_vel=1
        Eforce=zero
        Eforce_primitive=zero
        cell_pressure=pold(D_DECL(i,j,k),im_vel)
        if (cell_pressure.lt.zero) then
         cell_pressure=zero
        endif

        do dir=0,SDIM-1 
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
          print *,"dir out of range in EDGEPRESSURE"
          stop
         endif

         hx=xsten(1,dir+1)-xsten(-1,dir+1)
         RR=one
         if ((levelrz.eq.0).or.(levelrz.eq.1)) then
          RR=one
         else if (levelrz.eq.3) then
          if (dir.eq.1) then ! theta direction
           RR=xsten(0,1)
          else if ((dir.eq.0).or.(dir.eq.SDIM-1)) then
           RR=one
          else
           print *,"dir invalid mac to cell"
           stop
          endif
         else
          print *,"levelrz invalid edge pressure 2"
          stop
         endif
         hx=hx*RR
         if (hx.le.zero) then
          print *,"hx invalid"
          stop
         endif

          ! side=1 left half of cell, side=2 right half of cell
         do side=1,2

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

          im_vel=1 

          if (dir.eq.0) then
           aface(side)=ax(D_DECL(iface,jface,kface))
           do im=1,ncphys
            ASIDE(side,im)=xface(D_DECL(iface,jface,kface),im)
           enddo
           uface(side,im_vel)=xvel(D_DECL(iface,jface,kface),im_vel)
           pres_face(side)=xp(D_DECL(iface,jface,kface),3)
           use_face_pres(side)=xp(D_DECL(iface,jface,kface),1)
           coarse_fine_face(side)=xp(D_DECL(iface,jface,kface),2)
          else if (dir.eq.1) then
           aface(side)=ay(D_DECL(iface,jface,kface))
           do im=1,ncphys
            ASIDE(side,im)=yface(D_DECL(iface,jface,kface),im)
           enddo
           uface(side,im_vel)=yvel(D_DECL(iface,jface,kface),im_vel)
           pres_face(side)=yp(D_DECL(iface,jface,kface),3)
           use_face_pres(side)=yp(D_DECL(iface,jface,kface),1)
           coarse_fine_face(side)=yp(D_DECL(iface,jface,kface),2)
          else if ((dir.eq.2).and.(SDIM.eq.3)) then
           aface(side)=az(D_DECL(iface,jface,kface))
           do im=1,ncphys
            ASIDE(side,im)=zface(D_DECL(iface,jface,kface),im)
           enddo
           uface(side,im_vel)=zvel(D_DECL(iface,jface,kface),im_vel)
           pres_face(side)=zp(D_DECL(iface,jface,kface),3)
           use_face_pres(side)=zp(D_DECL(iface,jface,kface),1)
           coarse_fine_face(side)=zp(D_DECL(iface,jface,kface),2)
          else
           print *,"dir invalid mac to cell 3"
           stop
          endif

          mass_side(side)=zero
          do im=1,nmat 
           if (side.eq.1) then  ! left half of cell
            sidecomp=massface_index+2*(im-1)+2
           else if (side.eq.2) then ! right half of cell
            sidecomp=massface_index+2*(im-1)+1
           else
            print *,"side invalid"
            stop
           endif
           mass_side(side)=mass_side(side)+ASIDE(side,sidecomp) 
          enddo ! im

         enddo ! side=1..2

         masscell=mass_side(1)+mass_side(2)

         if ((mass_side(1).le.zero).or. &
             (mass_side(2).le.zero).or. &
             (masscell.le.zero)) then
          print *,"mass invalid"
          stop
         endif

         dencell=masscell/VOLTERM
         if (dencell.le.zero) then
          print *,"dencell invalid"
          stop
         endif

         if (1.eq.0) then
          if (dir.eq.SDIM-1) then
           print *,"grad p^cell  i,j,k,dir,dencell ",i,j,k,dir,dencell
          endif
         endif

         cell_is_ice=0
         cell_is_FSI_rigid=0

         do im=1,nmat
          if (is_rigid(nmat,im).eq.1) then
           ! do nothing
          else if ((FSI_flag(im).eq.0).or. &
                   (FSI_flag(im).eq.7)) then
           ! do nothing
          else if (is_ice(nmat,im).eq.1) then
           if (LStest(im).ge.zero) then
            cell_is_ice=1
           endif
          else if (is_FSI_rigid(nmat,im).eq.1) then
           if (LStest(im).ge.zero) then
            cell_is_FSI_rigid=1
           endif
          else
           print *,"FSI_flag invalid"
           stop
          endif
         enddo ! im=1..nmat

         if (use_face_pres(1).eq.use_face_pres(2)) then
          use_face_pres_combine=use_face_pres(1)
         else if ((use_face_pres(1).eq.three).or. &
                  (use_face_pres(2).eq.three)) then
          use_face_pres_combine=min(use_face_pres(1),use_face_pres(2))
         else if ((use_face_pres(1).eq.zero).or. &
                  (use_face_pres(2).eq.zero)) then
          use_face_pres_combine=zero
         else if ((use_face_pres(1).eq.one).and. &
                  (use_face_pres(2).eq.two)) then
          use_face_pres_combine=zero
         else if ((use_face_pres(1).eq.two).and. &
                  (use_face_pres(2).eq.one)) then
          use_face_pres_combine=zero
         else
          print *,"use_face_pres invalid"
          stop
         endif
        
         if (use_face_pres_cen.eq.3) then
          use_face_pres_cen=NINT(use_face_pres_combine)
         else if (use_face_pres_cen.eq.0) then
          ! do nothing
         else if (use_face_pres_cen.eq.2) then
          if (use_face_pres_combine.eq.zero) then
           use_face_pres_cen=0
          else if (use_face_pres_combine.eq.one) then
           use_face_pres_cen=0
          else if (use_face_pres_combine.eq.two) then
           use_face_pres_cen=2
          else if (use_face_pres_combine.eq.three) then
           use_face_pres_cen=2
          else
           print *,"use_face_pres_combine invalid"
           stop
          endif
         else if (use_face_pres_cen.eq.1) then
          if (use_face_pres_combine.eq.zero) then
           use_face_pres_cen=0
          else if (use_face_pres_combine.eq.one) then
           use_face_pres_cen=1
          else if (use_face_pres_combine.eq.two) then
           use_face_pres_cen=0
          else if (use_face_pres_combine.eq.three) then
           use_face_pres_cen=1
          else
           print *,"use_face_pres_combine invalid"
           stop
          endif

         else
          print *,"use_face_pres_cen invalid"
          print *,"use_face_pres_cen=",use_face_pres_cen
          stop
         endif

         if ((cell_is_ice.eq.1).or. &
             (cell_is_FSI_rigid.eq.1)) then
          use_face_pres_cen=0
         endif

         ! do not overwrite veldest with (un-grad p)^Cell
         if ((use_face_pres_cen.eq.0).or. &
             (use_face_pres_cen.eq.1)) then 

          if (energyflag.eq.2) then  ! space-time, set grad p^cell = 0
           im_vel=1
           velcomp=dir+1
           ! this holds grad p^CELL when energyflag==2
           ustar(D_DECL(i,j,k),velcomp)=zero
          endif ! energyflag.eq.2

         ! overwrite veldest with (un-grad p)^Cell if use_face_pres_cen==2,3
         else if ((use_face_pres_cen.eq.2).or. &
                  (use_face_pres_cen.eq.3)) then 

          ! density=masscell/vol
          if ((energyflag.eq.0).or. &
              (energyflag.eq.1)) then

           im_vel=1
           velcomp=dir+1

           dp=dt*(pres_face(2)-pres_face(1))/(dencell*hx)
           veldest(D_DECL(i,j,k),velcomp)=ustar(D_DECL(i,j,k),velcomp)-dp

           if (1.eq.0) then
            if ((i.eq.0).and.(dir.eq.0)) then
             print *,"i,j,k ",i,j,k
             print *,"ustar: ",ustar(D_DECL(i,j,k),velcomp)
             print *,"veldest: ",veldest(D_DECL(i,j,k),velcomp)
             print *,"dp ",dp
            endif
           endif

           ! this is for space-time
          else if (energyflag.eq.2) then

           im_vel=1
           velcomp=dir+1
           ustar(D_DECL(i,j,k),velcomp)=(pres_face(2)-pres_face(1))/hx

          else
           print *,"energyflag invalid"
           stop
          endif

         else
          print *,"use_face_pres_cen bust"
          stop
         endif

         if (1.eq.0) then
          print *,"APPLYING GRADP(velocity is scaled)"
          print *,"i,j,k,dir ",i,j,k,dir
          print *,"veldest= ",veldest(D_DECL(i,j,k),dir+1)
          print *,"ustar= ",ustar(D_DECL(i,j,k),dir+1)
          print *,"level,finest_level ",level,finest_level
         endif

         if (energyflag.eq.0) then
          ! do nothing
         else if (energyflag.eq.1) then
          im_vel=1
          Eforce=Eforce- &
           dt*(aface(2)*uface(2,im_vel)*pres_face(2)- &
               aface(1)*uface(1,im_vel)*pres_face(1))/ &
              (dencell*VOLTERM)
          if (local_incomp.eq.0) then
           Eforce_primitive=Eforce_primitive- &
            dt*(aface(2)*uface(2,im_vel)- &
                aface(1)*uface(1,im_vel))*cell_pressure/ &
               (dencell*VOLTERM)
          else if (local_incomp.eq.1) then
           ! do nothing
          else
           print *,"local_incomp invalid"
           stop
          endif
         else if (energyflag.eq.2) then
          im_vel=1
          Eforce=Eforce+ &
           (aface(2)*uface(2,im_vel)*pres_face(2)- &
            aface(1)*uface(1,im_vel)*pres_face(1))/VOLTERM
          if (local_incomp.eq.0) then
           Eforce_primitive=Eforce_primitive+ &
            (aface(2)*uface(2,im_vel)- &
             aface(1)*uface(1,im_vel))*cell_pressure/VOLTERM
          else if (local_incomp.eq.1) then
           ! do nothing
          else
           print *,"local_incomp invalid"
           stop
          endif
         else
          print *,"energyflag invalid"
          stop
         endif

        enddo ! dir=0..sdim-1 (operation_flag.eq.3)  (grad p)_CELL, div(up)

        im_vel=1
        rhs(D_DECL(i,j,k),im_vel)=zero

         ! update the total energy in partial cells (regular project)
        if ((project_option.eq.0).or. &
            (project_option.eq.13)) then ! FSI_material_exists 1st project

         im_vel=1
 
         if ((use_face_pres_cen.eq.0).or. &
             (use_face_pres_cen.eq.2)) then

          Eforce=zero
          Eforce_primitive=zero
          rhs(D_DECL(i,j,k),im_vel)=zero

         else if ((use_face_pres_cen.eq.1).or. &
                  (use_face_pres_cen.eq.3)) then          

          if (local_primitive.eq.0) then
           rhs(D_DECL(i,j,k),im_vel)=Eforce
          else if (local_primitive.eq.1) then
           rhs(D_DECL(i,j,k),im_vel)=Eforce_primitive
          else
           print *,"local_primitve invalid"
           stop
          endif

          if (energyflag.eq.1) then ! update the temperature

           do im=1,nmat

            KE_diff=zero
            do dir2=1,SDIM
             velcomp=dir2
             KE_diff=KE_diff+ &
              half*ustar(D_DECL(i,j,k),velcomp)**2- &
              half*veldest(D_DECL(i,j,k),velcomp)**2
            enddo ! dir2

            if (vfrac(im).gt.VOFTOL) then

             imattype=fort_material_type(im)

             ibase=(im-1)*num_state_material
            
             rho=dendest(D_DECL(i,j,k),ibase+1)
             TEMPERATURE=dendest(D_DECL(i,j,k),ibase+2)

             if (rho.le.zero) then
              print *,"density underflow"
              print *,"i,j,k,im,rho ",i,j,k,im,rho
              stop
             endif
             if (TEMPERATURE.le.zero) then
              print *,"temperature underflow"
              stop
             endif

             call INTERNAL_material(rho,TEMPERATURE,internal_e, &
               imattype,im)

             if (temperature_primitive_variable(im).eq.0) then
              internal_e=internal_e+KE_diff+Eforce
             else if (temperature_primitive_variable(im).eq.1) then
              internal_e=internal_e+Eforce_primitive
             else
              print *,"temp prim var invalid"
              stop
             endif

             if (internal_e.le.zero) then
              NEW_TEMPERATURE=TEMPERATURE
             else
              call TEMPERATURE_material(rho,NEW_TEMPERATURE, &
               internal_e,imattype,im)
             endif

             if (NEW_TEMPERATURE.gt.zero) then
              dendest(D_DECL(i,j,k),ibase+2)=NEW_TEMPERATURE
             else
              print *,"NEW_TEMPERATURE must be positive"
              stop
             endif
            endif ! vfrac(im)>voftol
           enddo ! im=1..nmat

          else if (energyflag.eq.2) then
           ! do nothing
          else if (energyflag.eq.0) then
           ! do nothing
          else
           print *,"energyflag invalid" 
           stop
          endif

         else 
          print *,"use_face_pres_cen invalid(2)"
          print *,"use_face_pres_cen=",use_face_pres_cen
          stop
         endif

        else if (project_option.eq.1) then
         ! do nothing if initial project
        else if (project_option.eq.10) then
         ! do nothing if sync project
        else if (project_option.eq.11) then !FSI_material_exists 2nd project
         ! do nothing if rigid body project
        else
         print *,"project_option invalid edge pressure 5"
         stop
        endif

       else if (operation_flag.eq.4) then ! (grad ppot)_CELL and surface ten.

        if (nsolve.ne.1) then
         print *,"nsolve invalid7"
         stop
        endif

        do dir=0,SDIM-1 
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
          print *,"dir out of range in mac to cell"
          stop
         endif

         hx=xsten(1,dir+1)-xsten(-1,dir+1)
         RR=one
         if ((levelrz.eq.0).or.(levelrz.eq.1)) then
          RR=one
         else if (levelrz.eq.3) then
          if (dir.eq.1) then ! theta direction
           RR=xsten(0,1)
          else if ((dir.eq.0).or.(dir.eq.SDIM-1)) then
           RR=one
          else
           print *,"dir invalid mac to cell"
           stop
          endif
         else
          print *,"levelrz invalid mac to cell 3"
          stop
         endif
         hx=hx*RR
         if (hx.le.zero) then
          print *,"hx invalid"
          stop
         endif

         dencellgrav=denold(D_DECL(i,j,k),1) ! hydrostatic density

         if (dencellgrav.le.zero) then
          print *,"hydrostatic density must be positive"
          stop
         endif

         ! side=1 left half of cell, side=2 right half of cell
         do side=1,2

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
 
          if (dir.eq.0) then 
           do im=1,ncphys
            ASIDE(side,im)=xface(D_DECL(iface,jface,kface),im)
           enddo
           pfacegrav(side)=xp(D_DECL(iface,jface,kface),1)
           pfacetenleft(side)=xp(D_DECL(iface,jface,kface),2)
           pfacetenright(side)=xp(D_DECL(iface,jface,kface),3)
          else if (dir.eq.1) then
           do im=1,ncphys
            ASIDE(side,im)=yface(D_DECL(iface,jface,kface),im)
           enddo
           pfacegrav(side)=yp(D_DECL(iface,jface,kface),1)
           pfacetenleft(side)=yp(D_DECL(iface,jface,kface),2)
           pfacetenright(side)=yp(D_DECL(iface,jface,kface),3)
          else if ((dir.eq.2).and.(SDIM.eq.3)) then
           do im=1,ncphys
            ASIDE(side,im)=zface(D_DECL(iface,jface,kface),im)
           enddo
           pfacegrav(side)=zp(D_DECL(iface,jface,kface),1)
           pfacetenleft(side)=zp(D_DECL(iface,jface,kface),2)
           pfacetenright(side)=zp(D_DECL(iface,jface,kface),3)
          else
           print *,"dir invalid mac to cell 5"
           stop
          endif

          if (side.eq.1) then
           pfaceten(side)=pfacetenright(side) ! surf ten, rt side of lt face
          else if (side.eq.2) then
           pfaceten(side)=pfacetenleft(side)  ! surf ten, lt side of rt face
          else
           print *,"side invalid"
           stop
          endif

          mass_side(side)=zero
          do im=1,nmat 
           if (side.eq.1) then  ! left half of cell
            sidecomp=massface_index+2*(im-1)+2
           else if (side.eq.2) then ! right half of cell
            sidecomp=massface_index+2*(im-1)+1
           else
            print *,"side invalid"
            stop
           endif
           mass_side(side)=mass_side(side)+ASIDE(side,sidecomp) 
          enddo ! im

         enddo ! side

         masscell=mass_side(1)+mass_side(2)

         if ((mass_side(1).le.zero).or.(mass_side(2).le.zero).or. &
             (masscell.le.zero)) then
          print *,"mass invalid"
          stop
         endif

         dencell=masscell/VOLTERM
         if (dencell.le.zero) then
          print *,"dencell invalid"
          stop
         endif

         if (1.eq.0) then
          if (dir.eq.SDIM-1) then
           print *,"grad pot^cell  i,j,k,dir,dencell,dencellgrav ", &
            i,j,k,dir,dencell,dencellgrav
          endif
         endif

          ! p=-|g| dt z
          ! force=dp/dz=-|g| dt
         pgrad=(pfacegrav(2)-pfacegrav(1))/(dencellgrav*hx)

         dp=dt*(pfaceten(2)-pfaceten(1))/(dencell*hx)

         pgrad=pgrad-dp

         veldest(D_DECL(i,j,k),dir+1)=pgrad

         if (1.eq.0) then
          if ((i.eq.0).or.(i.eq.1)) then
           print *,"i,j,k,dir ",i,j,k,dir
           print *,"dt,pfacegrav(1),pfacegrav(2),dencellgrav,hx ", &
            dt,pfacegrav(1),pfacegrav(2),dencellgrav,hx
           print *,"veldest: ",veldest(D_DECL(i,j,k),dir+1)
           print *,"dp ",dp
          endif
         endif

        enddo ! dir=0..sdim-1
      
       else if (operation_flag.eq.6) then ! advection
        ! low order approximation: CISL or sem_mac_to_cell
        ! high order approximation: sem_mac_to_cell
       else
        print *,"operation_flag invalid8"
        stop
       endif

      enddo
      enddo
      enddo

       ! for advection:
       !  1. low order fluxes are calculated in FORT_CELL_TO_MAC
       !     and high order fluxes are calculated in SEM_CELL_TO_MAC
       !  2. divergence term and SDC correction term are both
       !     calculated in SEM_MAC_TO_CELL (PROB.F90) regardless of the
       !     order.
      high_order_time_advection=0
      if (operation_flag.eq.6) then ! advection
       if ((ns_time_order.ge.2).and. &
           (ns_time_order.le.32)) then
        high_order_time_advection=1
       else if (ns_time_order.eq.1) then
        high_order_time_advection=0
       else
        print *,"ns_time_order invalid"
        stop
       endif
      endif
      if (ns_time_order.eq.1) then
       if ((enable_spectral.eq.0).or.(enable_spectral.eq.2)) then
        ! do nothing
       else
        print *,"enable_spectral and ns_time_order are not consistent"
        stop
       endif
      else if ((ns_time_order.ge.2).and.(ns_time_order.le.32)) then
       if ((enable_spectral.eq.0).or. &  ! disable high order time for visc
           (enable_spectral.eq.1).or. &
           (enable_spectral.eq.3)) then
        ! do nothing
       else
        print *,"enable_spectral and ns_time_order are not consistent"
        stop
       endif
      else
       print *,"ns_time_order invalid"
       stop
      endif

       ! in: FORT_MAC_TO_CELL
      if ((enable_spectral.eq.1).or. &  ! SEM space and time
          (enable_spectral.eq.2).or. &  ! SEM space
          (high_order_time_advection.eq.1)) then

       if ((bfact.ge.2).or. &
           (high_order_time_advection.eq.1)) then

        call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

         local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))

         if ((local_maskSEM.ge.1).and. &
             (local_maskSEM.le.nmat)) then

          call strip_status(i,j,k,bfact,stripstat)

          if (stripstat.eq.1) then

           do dir=0,SDIM-1

            call elementbox(i,j,k,bfact,dir,elemlo,elemhi)
            do ielem=elemlo(1),elemhi(1)
            do jelem=elemlo(2),elemhi(2)
            do kelem=elemlo(3),elemhi(3)
          
             if (operation_flag.eq.0) then ! RHS for solver
              scomp=1
              scomp_bc=dir+1
              dcomp=1
              ncomp=nsolve
              ncomp_xvel=nsolveMM_FACE
              ncomp_cterm=nsolve*num_materials_face
             else if (operation_flag.eq.1) then ! divergence
              scomp=1
              scomp_bc=dir+1
              dcomp=1
              ncomp=1
              ncomp_xvel=nsolveMM_FACE
              ncomp_cterm=1
             ! MAC->CELL in solver or VELMAC_to_CELL
             else if (operation_flag.eq.2) then 
              scomp=1
              scomp_bc=dir+1
              dcomp=dir+1
              ncomp=1
              ncomp_xvel=nsolveMM_FACE
              ncomp_cterm=1
             else if (operation_flag.eq.3) then ! (grad p)^CELL, div(up)
              scomp=1
              scomp_bc=dir+1
              dcomp=dir+1
              ncomp=1
              ncomp_xvel=nsolveMM_FACE
              ncomp_cterm=1
             else if (operation_flag.eq.4) then ! (grad pot)^CELL
              scomp=1
              scomp_bc=dir+1
              dcomp=dir+1
              ncomp=1
              ncomp_xvel=nsolveMM_FACE
              ncomp_cterm=1
             else if (operation_flag.eq.6) then ! advection
              scomp=1
              scomp_bc=1
              dcomp=1
              ncomp=ncphys
              ncomp_xvel=nsolveMM_FACE
              ncomp_cterm=SDIM+num_state_base
             else
              print *,"operation_flag invalid9"
              stop
             endif
              
             if (dir.eq.0) then 

              call SEM_MAC_TO_CELL( &
               ncomp_denold, &
               ncomp_veldest, &
               ncomp_dendest, &
               conservative_div_uu, &
               nsolveMM_FACE, &
               num_materials_face, &
               ns_time_order, &
               divu_outer_sweeps, &
               num_divu_outer_sweeps, &
               SDC_outer_sweeps, &
               SEM_advection_algorithm, &
               level, &
               finest_level, &
               nmat, &
               operation_flag, & 
               project_option, &
               energyflag, &
               temperature_primitive_variable, &
               face_flag, &
               homflag, &
               local_maskSEM, &
               cur_time, &
               slab_step, &
               dt, &
               ielem,jelem,kelem, &
               tilelo,tilehi, &
               fablo,fabhi, &
               xlo,dx,dir+1,bfact, &
               velbc_in, &  
               presbc_in, &
               scomp,scomp_bc, &
               dcomp, &
               ncomp, &
               ncomp_xvel, &
               ncomp_cterm, &
               vol,DIMS(vol), &
               xface,DIMS(xface), &
               xp,DIMS(xp), &
               xvel,DIMS(xvel), &
               maskcoef,DIMS(maskcoef), & ! 1=not covered, 0=covered
               cterm,DIMS(cterm), &
               mdotcell,DIMS(mdotcell), &
               pold,DIMS(pold), &
               denold,DIMS(denold), &
               ustar,DIMS(ustar), &
               veldest,DIMS(veldest), &
               dendest,DIMS(dendest), &
               rhs,DIMS(rhs) )  ! divdest

             else if (dir.eq.1) then

              call SEM_MAC_TO_CELL( &
               ncomp_denold, &
               ncomp_veldest, &
               ncomp_dendest, &
               conservative_div_uu, &
               nsolveMM_FACE, &
               num_materials_face, &
               ns_time_order, &
               divu_outer_sweeps, &
               num_divu_outer_sweeps, &
               SDC_outer_sweeps, &
               SEM_advection_algorithm, &
               level, &
               finest_level, &
               nmat, &
               operation_flag, & 
               project_option, &
               energyflag, &
               temperature_primitive_variable, &
               face_flag, &
               homflag, &
               local_maskSEM, &
               cur_time, &
               slab_step, &
               dt, &
               ielem,jelem,kelem, &
               tilelo,tilehi, &
               fablo,fabhi, &
               xlo,dx,dir+1,bfact, &
               velbc_in, &  
               presbc_in, &
               scomp,scomp_bc, &
               dcomp, &
               ncomp, &
               ncomp_xvel, &
               ncomp_cterm, &
               vol,DIMS(vol), &
               yface,DIMS(yface), &
               yp,DIMS(yp), &
               yvel,DIMS(yvel), &
               maskcoef,DIMS(maskcoef), & ! 1=not covered, 0=covered
               cterm,DIMS(cterm), &
               mdotcell,DIMS(mdotcell), &
               pold,DIMS(pold), &
               denold,DIMS(denold), &
               ustar,DIMS(ustar), &
               veldest,DIMS(veldest), &
               dendest,DIMS(dendest), &
               rhs,DIMS(rhs) )  ! divdest

             else if ((dir.eq.2).and.(SDIM.eq.3)) then

              call SEM_MAC_TO_CELL( &
               ncomp_denold, &
               ncomp_veldest, &
               ncomp_dendest, &
               conservative_div_uu, &
               nsolveMM_FACE, &
               num_materials_face, &
               ns_time_order, &
               divu_outer_sweeps, &
               num_divu_outer_sweeps, &
               SDC_outer_sweeps, &
               SEM_advection_algorithm, &
               level, &
               finest_level, &
               nmat, &
               operation_flag, & 
               project_option, &
               energyflag, &
               temperature_primitive_variable, &
               face_flag, &
               homflag, &
               local_maskSEM, &
               cur_time, &
               slab_step, &
               dt, &
               ielem,jelem,kelem, &
               tilelo,tilehi, &
               fablo,fabhi, &
               xlo,dx,dir+1,bfact, &
               velbc_in, &  
               presbc_in, &
               scomp,scomp_bc, &
               dcomp, &
               ncomp, &
               ncomp_xvel, &
               ncomp_cterm, &
               vol,DIMS(vol), &
               zface,DIMS(zface), &
               zp,DIMS(zp), &
               zvel,DIMS(zvel), &
               maskcoef,DIMS(maskcoef), & ! 1=not covered, 0=covered
               cterm,DIMS(cterm), &
               mdotcell,DIMS(mdotcell), &
               pold,DIMS(pold), &
               denold,DIMS(denold), &
               ustar,DIMS(ustar), &
               veldest,DIMS(veldest), &
               dendest,DIMS(dendest), &
               rhs,DIMS(rhs) ) ! divdest

             else
              print *,"dir invalid mac_to_cell2 "
              stop
             endif

            enddo
            enddo
            enddo !ielem,jelem,kelem
  
           enddo ! dir=0..sdim-1

          else if (stripstat.eq.0) then
           ! do nothing
          else
           print *,"stripstat invalid"
           stop
          endif

         else if (local_maskSEM.eq.0) then

          ! do nothing

         else
          print *,"local_maskSEM invalid"
          stop
         endif 

        enddo
        enddo
        enddo ! i,j,k

       else if ((bfact.eq.1).and. &
                (high_order_time_advection.eq.0)) then
        ! do nothing
       else
        print *,"bfact or high_order_time_advection invalid"
        stop
       endif

      else if (((enable_spectral.eq.0).or. &
                (enable_spectral.eq.3)).and. &
               (high_order_time_advection.eq.0)) then
       ! do nothing
      else
       print *,"enable_spectral or high_order_time_advection invalid"
       stop
      endif

      if (DO_SANITY_CHECK.eq.1) then
       if ((project_option.eq.0).and. &
           (operation_flag.eq.3)) then
        call CISL_projection(dt,divu_outer_sweeps+1)
        allocate(comparepface(fablo(1)-1:fabhi(1)+1,5))
        allocate(comparevelface(fablo(1)-1:fabhi(1)+1,5))
        allocate(comparestate(fablo(1)-1:fabhi(1)+1,5))
        j=1
        k=0
        do i=fablo(1),fabhi(1)
         comparepface(i,1)=xp(D_DECL(i,j,k),2)
        enddo
        i=fabhi(1)+1
        comparepface(i,1)=xp(D_DECL(i,j,k),1)
        do i=fablo(1),fabhi(1)+1
         comparevelface(i,1)=xvel(D_DECL(i,j,k),1)
        enddo
        do i=fablo(1),fabhi(1)
         comparestate(i,1)=veldest(D_DECL(i,j,k),1)
         ibase=0
         comparestate(i,2)=dendest(D_DECL(i,j,k),ibase+1)
         comparestate(i,3)=dendest(D_DECL(i,j,k),ibase+2)
        enddo
        call compare_sanity(comparepface,1,1,5)
        call compare_sanity(comparevelface,1,1,6)
        call compare_sanity(comparestate,1,3,7)
        deallocate(comparepface)
        deallocate(comparevelface)
        deallocate(comparestate)
       endif ! project_option=0 dir=0
      endif  ! do_sanity_check=true
 
      return
      end subroutine FORT_MAC_TO_CELL

!if temperature_primitive_var==0,
! add beta * (1/cv) * (u dot u/2) to temp
      subroutine FORT_INC_TEMP( &
       beta, &
       temperature_primitive_variable, &
       nmat, &
       level, &
       finest_level, &
       ncomp_state, &
       tilelo,tilehi, &
       fablo,fabhi, &
       state,DIMS(state), &
       maskcoef,DIMS(maskcoef)) ! 1=not covered  0=covered
       use probf90_module
       use global_utility_module
       IMPLICIT NONE

      REAL_T beta
      INTEGER_T nmat
      INTEGER_T ncomp_state
      INTEGER_T temperature_primitive_variable(nmat)
      INTEGER_T level,finest_level
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T DIMDEC(state)
      INTEGER_T DIMDEC(maskcoef)
      REAL_T state(DIMV(state),ncomp_state)
      REAL_T  maskcoef(DIMV(maskcoef))

      INTEGER_T i,j,k
      INTEGER_T dir
      INTEGER_T im
      INTEGER_T imattype
      INTEGER_T local_mask
      INTEGER_T vofcomp,dencomp
      REAL_T vof,KE,rho,TEMPERATURE,internal_e

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid INC_TEMP"
       stop
      endif

      do im=1,nmat
       imattype=fort_material_type(im)
       if (imattype.eq.999) then
        ! do nothing
       else if ((imattype.ge.0).and. &
                (imattype.le.fort_max_num_eos)) then
        ! do nothing
       else
        print *,"imattype invalid"
        stop
       endif
       if ((temperature_primitive_variable(im).ne.0).and. &
           (temperature_primitive_variable(im).ne.1)) then
        print *,"temperature_primitive_variable invalid"
        stop
       endif

       if ((fort_material_type(im).eq.0).or. &
           (is_rigid(nmat,im).eq.1).or. &
           (fort_material_type(im).eq.999)) then
        if (temperature_primitive_variable(im).ne.1) then
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (is_rigid(nmat,im).eq.0).and. &
                (fort_material_type(im).ne.999)) then
        if ((temperature_primitive_variable(im).eq.0).or. &
            (temperature_primitive_variable(im).eq.1)) then
         ! do nothing
        else
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else
        print *,"fort_material_type(im) or is_rigid invalid"
        stop
       endif

      enddo ! im=1..nmat

      if ((beta.ne.one).and.(beta.ne.-one)) then
       print *,"beta invalid"
       stop
      endif 
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      do im=1,nmat
       if (fort_denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif
      enddo

      if (ncomp_state.ne.(SDIM+1)*num_materials_vel+ &
          nmat*num_state_material+nmat*ngeom_raw+1) then
       print *,"ncomp_state invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(state),0,-1,33)
      call checkbound(fablo,fabhi,DIMS(maskcoef),0,-1,134)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_mask=NINT(maskcoef(D_DECL(i,j,k)))
       if (local_mask.eq.1) then ! not covered

        do im=1,nmat
         if (temperature_primitive_variable(im).eq.1) then
          ! do nothing
         else if (temperature_primitive_variable(im).eq.0) then
          vofcomp=(SDIM+1)*num_materials_vel+nmat*num_state_material+ &
            (im-1)*ngeom_raw+1
          vof=state(D_DECL(i,j,k),vofcomp)
          if ((vof.ge.-VOFTOL).and.(vof.le.one+VOFTOL)) then

           if (vof.lt.VOFTOL) then
            ! do nothing
           else if (vof.ge.VOFTOL) then
            KE=zero
            do dir=1,SDIM
             KE=KE+half*(state(D_DECL(i,j,k),dir)**2)
            enddo
            dencomp=(SDIM+1)*num_materials_vel+ &
             (im-1)*num_state_material+1
            rho=state(D_DECL(i,j,k),dencomp)
            TEMPERATURE=state(D_DECL(i,j,k),dencomp+1)
            if (rho.le.zero) then
             print *,"density underflow"
             print *,"i,j,k,im,rho ",i,j,k,im,rho
             stop
            endif
             ! the first part is always to add the old kinetic energy.
            if (TEMPERATURE.le.zero) then
             print *,"temperature underflow"
             stop
            endif
            imattype=fort_material_type(im)
            if ((imattype.gt.0).and. &
                (imattype.le.fort_max_num_eos)) then
             call INTERNAL_material(rho,TEMPERATURE,internal_e, &
               imattype,im)
             internal_e=internal_e+beta*KE
             if (internal_e.le.zero) then
              print *,"internal_e.le.zero in INC_TEMP"
              stop
             endif
             call TEMPERATURE_material(rho,TEMPERATURE, &
               internal_e,imattype,im)
             state(D_DECL(i,j,k),dencomp+1)=TEMPERATURE
            else
             print *,"imattype invalid"
             stop
            endif
           else
            print *,"vof invalid"
            stop
           endif

          else
           print *,"vof out of range"
           stop
          endif

         else
          print *,"temperature_primitive_variable invalid"
          stop
         endif
        enddo ! im=1..nmat

       else if (local_mask.eq.0) then ! covered
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif
      enddo !k
      enddo !j
      enddo !i
 
      return
      end subroutine FORT_INC_TEMP

! operation_flag=0  pressure gradient on MAC grid
! operation_flag=1  interpolate pressure from cell to MAC grid.
! operation_flag=2  potential gradient on MAC grid, 
!                   interpolate potential from cell to MAC grid
!                   surface tension on MAC grid
!                   left and right surface tension on MAC grid.
! operation_flag=3  unew^MAC=unew^CELL->MAC
! operation_flag=4  unew^MAC=uSOLID^MAC or uFLUID^MAC
! operation_flag=5  unew^MAC=unew^MAC+beta diffuse_ref^CELL->MAC
! (operation_flag=6 reserved for rate of strain tensor)
! operation_flag=7  advection.
! operation_flag=8  reserved for coupling terms in CROSSTERM
! operation_flag=9  density cell -> MAC
! operation_flag=10 hybrid of 3 and 4
! operation_flag=11 
!   (i) unew^{f} in incompressible non-solid regions
!   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions or
!        compressible regions.
!   (iii) usolid in solid regions

      subroutine FORT_CELL_TO_MAC( &
       ncomp_xp, &  !local_MF[AMRSYNC_PRES_MF]->nComp() if operation_flag==0
       ncomp_xgp, &
       simple_AMR_BC_flag, &
       nsolveMM_FACE, &
       num_materials_face, &
       tileloop, &
       dir, &
       operation_flag, & 
       energyflag, & 
       beta, &
       visc_coef, &
       face_flag, & 
       temperature_primitive_variable, &
       enable_spectral, &
       fluxvel_index, &  
       fluxden_index, &  
       facevel_index, &  
       facecut_index, &
       icefacecut_index, &
       curv_index, &
       conservative_tension_force, &
       conservative_div_uu, &
       pforce_index, &
       faceden_index, &  
       icemask_index, &
       massface_index, &
       vofface_index, &
       ncphys, &  ! nflux for advection
       make_interface_incomp, &
       override_density, &
       solvability_projection, &
       presbc_in, &  ! denbc for advection
       velbc_in, &
       slab_step, &
       dt, &
       time, &
       xlo,dx, &
       spectral_loop, &
       ncfluxreg, &
       semflux,DIMS(semflux), &
       mask,DIMS(mask), & ! 1=fine/fine  0=coarse/fine
       maskcoef,DIMS(maskcoef), & ! 1=not cov. or outside domain  0=covered
       maskSEM,DIMS(maskSEM), &
       levelPC,DIMS(levelPC), &
       sol,DIMS(sol), &
       xcut,DIMS(xcut), &   ! coeff*areafrac, local_xface in increment_face_vel
       xface,DIMS(xface), &  ! xflux for advection
       xfacemm,DIMS(xfacemm), &  
       xcellmm,DIMS(xcellmm), &  
       recon,DIMS(recon), &  
       xgp,DIMS(xgp), & ! holds Umac_old if operation_flag==5 or 11
       xp,DIMS(xp), & ! holds AMRSYNC_PRES if operation_flag==0
       xvel,DIMS(xvel), &
       vel,DIMS(vel), &
       pres,DIMS(pres), & ! holds U_old(dir) if operation_flag==11
       den,DIMS(den), &
       mgoni,DIMS(mgoni), &
       colorfab,DIMS(colorfab), &
       typefab,DIMS(typefab), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact,bfact_c,bfact_f, &
       level,finest_level, &
       rz_flag, &
       domlo,domhi, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       prescribed_solid_scale, &
       added_weight, &
       blob_array, &
       blob_array_size, &
       num_elements_blobclass, &
       num_colors, &
       nten, &
       nfacefrac, &
       ncellfrac, &
       project_option, &
       SEM_upwind, &
       SEM_advection_algorithm)
      use global_utility_module
      use MOF_routines_module
      use probf90_module
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, intent(in) :: ncomp_xp
      INTEGER_T, intent(in) :: ncomp_xgp
      INTEGER_T, intent(in) :: simple_AMR_BC_flag
      INTEGER_T, intent(in) :: num_materials_face
      INTEGER_T, intent(in) :: nsolveMM_FACE
      INTEGER_T, intent(in) :: tileloop
      INTEGER_T, intent(in) :: spectral_loop
      INTEGER_T, intent(in) :: ncfluxreg
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)

      INTEGER_T, intent(in) :: blob_array_size
      INTEGER_T, intent(in) :: num_elements_blobclass
      INTEGER_T, intent(in) :: num_colors
      REAL_T, intent(in) :: blob_array(blob_array_size)
        
      REAL_T, intent(in) :: prescribed_solid_scale(nmat)
      REAL_T, intent(in) :: added_weight(nmat)
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: nfacefrac
      INTEGER_T, intent(in) :: ncellfrac
      INTEGER_T, intent(in) :: slab_step
      INTEGER_T, intent(in) :: face_flag 
      INTEGER_T, intent(in) :: temperature_primitive_variable(nmat)
      INTEGER_T, intent(in) :: operation_flag
      INTEGER_T, intent(in) :: energyflag
      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: fluxvel_index 
      INTEGER_T, intent(in) :: fluxden_index 
      INTEGER_T, intent(in) :: facevel_index 
      INTEGER_T, intent(in) :: facecut_index
      INTEGER_T, intent(in) :: icefacecut_index
      INTEGER_T, intent(in) :: curv_index
      INTEGER_T, intent(in) :: conservative_tension_force
      INTEGER_T, intent(in) :: conservative_div_uu
      INTEGER_T, intent(in) :: pforce_index
      INTEGER_T, intent(in) :: faceden_index 
      INTEGER_T, intent(in) :: icemask_index
      INTEGER_T, intent(in) :: massface_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys  ! nflux for advection
      INTEGER_T, intent(in) :: make_interface_incomp
      INTEGER_T, intent(in) :: override_density(nmat)
      INTEGER_T, intent(in) :: solvability_projection
      REAL_T, intent(in) :: dt,time,beta,visc_coef
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(semflux)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(maskcoef)
      INTEGER_T, intent(in) :: DIMDEC(maskSEM)
      INTEGER_T, intent(in) :: DIMDEC(xcut)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(xfacemm)
      INTEGER_T, intent(in) :: DIMDEC(xcellmm)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(xgp)
      INTEGER_T, intent(in) :: DIMDEC(xp)
      INTEGER_T, intent(in) :: DIMDEC(xvel)
      INTEGER_T, intent(in) :: DIMDEC(pres)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(den)
      INTEGER_T, intent(in) :: DIMDEC(mgoni)
      INTEGER_T, intent(in) :: DIMDEC(levelPC)
      INTEGER_T, intent(in) :: DIMDEC(sol)
      INTEGER_T, intent(in) :: DIMDEC(colorfab)
      INTEGER_T, intent(in) :: DIMDEC(typefab)

       ! denbc for advect
      INTEGER_T, intent(in) :: presbc_in(SDIM,2,nmat*num_state_material) 
      INTEGER_T, intent(in) :: velbc_in(SDIM,2,num_materials_vel*SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact,bfact_c,bfact_f
      INTEGER_T, intent(in) :: rz_flag
      INTEGER_T, intent(in) :: domlo(SDIM),domhi(SDIM)
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: SEM_upwind
      INTEGER_T, intent(in) :: SEM_advection_algorithm

      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(in) :: maskcoef(DIMV(maskcoef))

      REAL_T, intent(in) :: maskSEM(DIMV(maskSEM))
      REAL_T, intent(in) :: levelPC(DIMV(levelPC),nmat*(1+SDIM))
      REAL_T, intent(in) :: sol(DIMV(sol),nparts_def*SDIM)
      REAL_T, intent(inout) :: semflux(DIMV(semflux),ncfluxreg)
      REAL_T, intent(inout) :: xcut(DIMV(xcut),1)
      REAL_T, intent(inout) :: xface(DIMV(xface),ncphys) ! xflux for advection
      REAL_T, intent(inout) :: xfacemm(DIMV(xfacemm),nfacefrac)
      REAL_T, intent(inout) :: xcellmm(DIMV(xfacemm),ncellfrac)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon)
       !holds Umac_old if operation_flag==5 or 11
      REAL_T, intent(inout) :: xgp(DIMV(xgp),ncomp_xgp) 

        ! for regular edge pressure operation:
        ! 1st component reserved for cell velocity override indicator.
        ! 2nd component reserved for coarse/fine indicator.
        ! for gravity/surface tension:
        ! 1st component: gravity edge pressure
        ! 2nd and 3rd components: surface tension edge pressures
      REAL_T, intent(inout) :: xp(DIMV(xp),ncomp_xp)
       ! xvel is destination for: density CELL->MAC (xvel=1/rho)
      REAL_T, intent(inout) :: xvel(DIMV(xvel),nsolveMM_FACE)
      REAL_T, intent(in) :: vel(DIMV(vel),num_materials_face*SDIM)
       ! holds U_old if operation_flag==11
      REAL_T, intent(in) :: pres(DIMV(pres),num_materials_face)
       ! den is the source for: density CELL->MAC
      REAL_T, intent(in) :: den(DIMV(den),nmat*num_state_material)
      REAL_T, intent(in) :: mgoni(DIMV(mgoni),nmat*num_state_material)
      REAL_T, intent(in) :: typefab(DIMV(typefab))
      REAL_T, intent(in) :: colorfab(DIMV(colorfab))
  
      INTEGER_T i,j,k,ii,jj,kk
      REAL_T pplus(num_materials_face)
      REAL_T pminus(num_materials_face)
      REAL_T pgrad(nsolveMM_FACE)
      REAL_T pgrad_gravity
      REAL_T pgrad_tension
      REAL_T gradh
      REAL_T dplus,dminus
       ! operation_flag==1: (1)use_face_pres,(2) grid flag, 2+im_vel,
       !   im_vel=1..nsolveMM_FACE
       ! operation_flag==2: (1)potential, (2-3) surface tension
      REAL_T plocal(2+nsolveMM_FACE)
      INTEGER_T im1,jm1,km1
      INTEGER_T im,im_opp,im_heat,tcomp,iten
      INTEGER_T dir,dir2,side
      INTEGER_T velcomp,iboundary
      REAL_T cutedge,RR
      REAL_T xstenMAC(-3:3,SDIM)
      REAL_T xmac(SDIM)
      INTEGER_T nhalf,nten_test
      REAL_T DXMAXLS
      REAL_T local_vel(nsolveMM_FACE)
      REAL_T local_vel_old(nsolveMM_FACE)
      REAL_T uedge(nsolveMM_FACE)
      REAL_T uedge_rigid
      REAL_T local_face(ncphys)
      INTEGER_T idx
      INTEGER_T is_solid_face
      INTEGER_T is_prescribed_face
      INTEGER_T at_RZ_face
      INTEGER_T ic,jc,kc
      REAL_T fluid_volface
      REAL_T not_prescribed_volface
      REAL_T volface
      REAL_T local_volume
      REAL_T massface,denface
      REAL_T mass(2)
      REAL_T vol_local(2)
      REAL_T den_local(2)
      REAL_T velsum
      REAL_T mass_sum
      REAL_T DMface
      REAL_T velmaterial
      REAL_T velmaterialMAC
      REAL_T solid_velocity
      INTEGER_T mask_covered(2)
      INTEGER_T mask_coarsefine(2)
      INTEGER_T at_reflect_wall,at_wall,at_ext_wall
      INTEGER_T use_face_pres
      INTEGER_T at_coarse_fine_wallF
      INTEGER_T at_coarse_fine_wallC
      REAL_T user_tension(nten)
      REAL_T tension_scaled
      REAL_T pforce_scaled
      REAL_T LSleft(nmat)
      REAL_T LSright(nmat)
      REAL_T LSupwind(nmat)
      REAL_T localLS(nmat)
      REAL_T mgoni_temp(nmat)
      INTEGER_T local_maskSEM
      INTEGER_T maskcov
      REAL_T hx
      INTEGER_T scomp,scomp_bc,dcomp
      INTEGER_T ncomp_dest,ncomp_source
      INTEGER_T update_right_flux
      INTEGER_T nc ! in: cell_to_mac
      INTEGER_T ibase,idonate,jdonate,kdonate
      INTEGER_T stripstat
      INTEGER_T elemlo(3),elemhi(3)
      INTEGER_T ielem,jelem,kelem
      REAL_T AFACE
      REAL_T AFACE_ICE
      REAL_T denlocal,templocal
      INTEGER_T imattype
      INTEGER_T im_vel
      INTEGER_T nsolveMM_FACE_test
      REAL_T test_velocity_FACE
      INTEGER_T ok_to_HO_interp
      INTEGER_T im_solid
      INTEGER_T im_prescribed
      INTEGER_T im_solid_valid
      INTEGER_T im_prescribed_valid
      INTEGER_T partid_solid
      INTEGER_T partid_prescribed
      INTEGER_T partid_check
      REAL_T wtleft,wtright
      INTEGER_T all_incomp,local_incomp
      REAL_T cutoff
      REAL_T local_tension_force
      INTEGER_T typeleft,typeright,typeface
      INTEGER_T colorleft,colorright,colorface
      INTEGER_T face_velocity_override
      INTEGER_T FSI_prescribed_flag
      REAL_T test_local_icefacecut
      REAL_T test_current_icefacecut
      REAL_T test_current_icemask
      INTEGER_T DEBUG_PRESCRIBED
      REAL_T DEBUG_PRESCRIBED_VEL_TOT
      REAL_T DEBUG_PRESCRIBED_VEL_DEN

      DEBUG_PRESCRIBED=0
      DEBUG_PRESCRIBED_VEL_TOT=zero
      DEBUG_PRESCRIBED_VEL_DEN=zero

      nhalf=3
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (ncomp_xp.lt.1) then
       print *,"ncomp_xp invalid"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_CELL_TO_MAC"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_CELL_TO_MAC"
       stop
      endif
      if ((simple_AMR_BC_flag.eq.0).or.(simple_AMR_BC_flag.eq.1)) then
       ! do nothing
      else
       print *,"simple_AMR_BC_flag invalid"
       stop
      endif

      nsolveMM_FACE_test=num_materials_face
      if (num_materials_face.eq.1) then
       ! do nothing
      else if (num_materials_face.eq.nmat) then
       nsolveMM_FACE_test=2*nsolveMM_FACE_test
      else
       print *,"num_materials_face invalid"
       stop
      endif
      if (nsolveMM_FACE_test.ne.nsolveMM_FACE) then
       print *,"nsolveMM_FACE invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((bfact_c.ne.bfact).and.(bfact_c.ne.2*bfact)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_f.ne.bfact).and.(bfact.ne.2*bfact_f)) then
       print *,"bfact_f invalid"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid cell to mac"
       stop
      endif
      if ((SEM_upwind.ne.0).and. &
          (SEM_upwind.ne.1)) then
       print *,"SEM_upwind invalid cell to mac"
       stop
      endif
      if ((SEM_advection_algorithm.eq.0).or. &
          (SEM_advection_algorithm.eq.1)) then
       ! do nothing
      else
       print *,"SEM_advection_algorithm invalid"
       stop
      endif
 
      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.3)) then
       print *,"enable_spectral invalid cell to mac"
       stop
      endif
      if ((spectral_loop.ne.0).and. &
          (spectral_loop.ne.1)) then
       print *,"spectral_loop invalid"
       stop
      endif
      if ((ncfluxreg/SDIM)*SDIM.ne.ncfluxreg) then
       print *,"ncfluxreg invalid11 ",ncfluxreg
       stop
      endif
      if (ncfluxreg.lt.SDIM) then
       print *,"ncfluxreg invalid12 ",ncfluxreg
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
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten.ne.nten_test) then
       print *,"nten invalid edge grad nten nten_test ",nten,nten_test
       stop
      endif
      if ((solvability_projection.ne.0).and. &
          (solvability_projection.ne.1)) then
       print *,"solvability_projection invalid"
       stop
      endif

      if ((make_interface_incomp.ne.0).and. &
          (make_interface_incomp.ne.1).and. &
          (make_interface_incomp.ne.2)) then
       print *,"make_interface_incomp invalid"
       stop
      endif

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
       print *,"num_elements_blobclass invalid"
       stop
      endif

      if ((blob_array_size.ne.1).and. &
          (blob_array_size.ne.num_colors*num_elements_blobclass)) then
       print *,"blob_array_size invalid"
       stop
      endif

      if (num_colors.ge.0) then
       ! do nothing
      else
       print *,"num_colors invalid"
       stop
      endif

      if ((operation_flag.ge.0).and. &
          (operation_flag.le.11)) then
       ! do nothing
      else
       print *,"operation_flag invalid in FORT_CELL_TO_MAC 10"
       stop
      endif

      if (operation_flag.eq.7) then ! advection

       if (ncomp_xp.ne.SDIM+num_state_base) then
        print *,"ncomp_xp invalid"
        stop
       endif
       if (ncomp_xgp.ne.SDIM+num_state_base) then
        print *,"ncomp_xp invalid"
        stop
       endif
       
       if (ncphys.ne.SDIM+num_state_base) then
        print *,"ncphys invalid"
        stop
       endif
       if (num_materials_face.ne.1) then
        print *,"num_materials_face invalid"
        stop
       endif

      else if (operation_flag.eq.1) then

       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif
       if (ncomp_xp.ne.2+nsolveMM_FACE) then
        print *,"ncomp_xp invalid"
        stop
       endif

      else if ((operation_flag.ge.0).and. &
               (operation_flag.le.2)) then
       
       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif

      else if ((operation_flag.ge.3).and. &
               (operation_flag.le.5)) then

       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif

      else if (operation_flag.eq.6) then

       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif

      else if (operation_flag.eq.10) then

       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif

      else if (operation_flag.eq.11) then

       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif

      else if (operation_flag.eq.9) then

       if (ncphys.ne.vofface_index+2*nmat) then
        print *,"ncphys invalid"
        stop
       endif

      else 
       print *,"operation_flag invalid11"
       stop
      endif

      ! indexes start at 0
      if ((curv_index.ne.0).or. &
          (pforce_index.ne.1).or. &
          (facecut_index.ne.3).or. &
          (icefacecut_index.ne.4).or. &
          (icemask_index.ne.5).or. &
          (facevel_index.ne.8).or. &
          (faceden_index.ne.2).or. &
          (vofface_index.ne.massface_index+2*nmat).or. &
          (fluxvel_index.ne.0).or. &
          (fluxden_index.ne.SDIM)) then
       print *,"face_index bust 2"
       stop
      endif

      if ((conservative_tension_force.ne.0).and. &
          (conservative_tension_force.ne.1)) then
       print *,"conservative_tension_force invalid"
       stop
      endif
      if ((conservative_div_uu.eq.0).or. &
          (conservative_div_uu.eq.1).or. &
          (conservative_div_uu.eq.2)) then
       ! do nothing
      else
       print *,"conservative_div_uu invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (operation_flag.eq.1) then ! p^CELL->MAC
       if (num_materials_face.ne.1) then
        print *,"num_materials_face invalid"
        stop
       endif
       if (energyflag.ne.0) then
        print *,"energyflag invalid"
        stop
       endif
      else if (operation_flag.eq.2) then !(grd ppot/den)_MAC,ppot^CELL->MAC,ten.
       if (energyflag.ne.0) then
        print *,"energyflag invalid"
        stop
       endif
       if ((ncomp_xp.ne.3).or.(ncomp_xgp.ne.1)) then
        print *,"ncomp_xp or ncomp_xgp invalid"
        stop
       endif
      else if ((operation_flag.eq.3).or. & !unew^CELL->MAC
               (operation_flag.eq.4).or. & !unew^MAC=uSOLID^MAC / uFLUID^MAC
               (operation_flag.eq.5).or. & !unew^MAC+beta FVISC^CELL->MAC
               (operation_flag.eq.10).or. & !unew^{CELL,MAC} -> MAC
               (operation_flag.eq.11)) then !unew^{CELL diff,MAC} -> MAC
       if (energyflag.ne.0) then
        print *,"energyflag invalid"
        stop
       endif
       if (num_materials_face.ne.1) then
        print *,"num_materials_face invalid"
        stop
       endif
       if ((project_option.eq.0).or. &
           (project_option.eq.1).or. &
           (project_option.eq.3).or. & 
           (project_option.eq.10).or. &
           (project_option.eq.13).or. & !FSI_material_exists 1st project
           (project_option.eq.11)) then !FSI_material_exists 2nd project
        ! do nothing
       else
        print *,"project_option invalid"
        stop
       endif
      else if (operation_flag.eq.9) then ! density: CELL->MAC
       if (energyflag.ne.0) then
        print *,"energyflag invalid"
        stop
       endif
      else if (operation_flag.eq.0) then ! (grad p)_MAC
       if ((energyflag.ne.0).and.(energyflag.ne.2)) then
        print *,"energyflag invalid"
        stop
       endif
      else if (operation_flag.eq.7) then ! advection

       if (energyflag.ne.0) then
        print *,"energyflag invalid"
        stop
       endif

      else if (operation_flag.eq.6) then

       print *,"fort_face_gradients calls sem_cell_to_mac with op=6"
       stop

      else
       print *,"operation_flag invalid12"
       stop
      endif

      if ((project_option.eq.0).or. &
          (project_option.eq.1).or. &
          (project_option.eq.10).or. &
          (project_option.eq.11).or. & !FSI_material_exists 2nd project
          (project_option.eq.13).or. & !FSI_material_exists 1st project
          (project_option.eq.12).or. & !pressure extrapolation
          (project_option.eq.2).or. &  ! thermal diffusion
          (project_option.ge.3).or. &  ! viscosity
          ((project_option.ge.100).and. &
           (project_option.lt.100+num_species_var))) then
       ! do nothing
      else
       print *,"project_option invalid"
       stop
      endif

      do im=1,nmat

       if ((prescribed_solid_scale(im).ge.zero).and. &
           (prescribed_solid_scale(im).lt.half)) then
        ! do nothing
       else
        print *,"prescribed_solid_scale invalid"
        stop
       endif

       if (added_weight(im).gt.zero) then
        ! do nothing
       else
        print *,"added_weight invalid"
        stop
       endif

       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"denconst invalid"
        stop
       endif
       if ((override_density(im).eq.0).or. & ! Dp/Drho=-rho c^2 div u
           (override_density(im).eq.1).or. & ! rho=rho(T,z)
           (override_density(im).eq.2)) then ! Boussinesq
        ! do nothing
       else
        print *,"override_density invalid"
        stop
       endif

       if ((fort_material_type(im).eq.0).or. &
           (is_rigid(nmat,im).eq.1).or. &
           (fort_material_type(im).eq.999)) then
        if (temperature_primitive_variable(im).ne.1) then
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (is_rigid(nmat,im).eq.0).and. &
                (fort_material_type(im).ne.999)) then
        if ((temperature_primitive_variable(im).eq.0).or. &
            (temperature_primitive_variable(im).eq.1)) then
         ! do nothing
        else
         print *,"temperature_primitive_variable(im) invalid"
         stop
        endif
       else
        print *,"fort_material_type(im) or is_rigid invalid"
        stop
       endif

      enddo ! im=1..nmat

      if ((slab_step.lt.-1).or.(slab_step.gt.bfact_time_order)) then
       print *,"slab_step invalid cell to mac"
       stop
      endif

      if ((tileloop.eq.0).and.(spectral_loop.eq.0)) then
       call checkbound(fablo,fabhi,DIMS(xcut),0,dir,231)
       call checkbound(fablo,fabhi,DIMS(xface),0,dir,263)
       call checkbound(fablo,fabhi,DIMS(xfacemm),0,dir,264)
       call checkbound(fablo,fabhi,DIMS(xgp),0,dir,2330)
       call checkbound(fablo,fabhi,DIMS(xp),0,dir,2331)
       call checkbound(fablo,fabhi,DIMS(xvel),0,dir,2332)
       if (dir.eq.0) then
        call checkbound(fablo,fabhi,DIMS(semflux),1,-1,231)
        call checkbound(fablo,fabhi,DIMS(vel),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(pres),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(den),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(sol),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(mgoni),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(typefab),1,-1,6625)
        call checkbound(fablo,fabhi,DIMS(colorfab),1,-1,6626)
        call checkbound(fablo,fabhi,DIMS(levelPC),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(xcellmm),0,-1,234)
        call checkbound(fablo,fabhi,DIMS(recon),0,-1,234)
        call checkbound(fablo,fabhi,DIMS(mask),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(maskcoef),1,-1,234)
        call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,1264)
       endif
      endif

      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=DXMAXLS

      all_incomp=1

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).ge.1).and. &
                (fort_material_type(im).le.fort_max_num_eos)) then
        all_incomp=0
       else
        print *,"fort_material_type invalid"
        stop
       endif

      enddo  ! im=1..nmat

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
       print *,"dir out of range in EDGEGRADP"
       stop
      endif 

       ! first low order gradients.
      if (tileloop.eq.0) then

       if (spectral_loop.eq.0) then

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)
        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir+1)

         if (levelrz.eq.0) then
          RR=one
         else if (levelrz.eq.1) then
          RR=one
         else if (levelrz.eq.3) then
          if (dir.eq.1) then
           RR=xstenMAC(0,1)
          else
           RR=one
          endif
         else
          print *,"levelrz invalid edgegradp"
          stop
         endif 

         hx=xstenMAC(1,dir+1)-xstenMAC(-1,dir+1)

         wtleft=(xstenMAC(1,dir+1)-xstenMAC(0,dir+1))/hx
         wtright=(xstenMAC(0,dir+1)-xstenMAC(-1,dir+1))/hx
         if ((wtleft.le.zero).or.(wtleft.ge.one)) then
          print *,"wtleft invalid"
          stop
         endif
         if ((wtright.le.zero).or.(wtright.ge.one)) then
          print *,"wtright invalid"
          stop
         endif
         if (abs(one-wtleft-wtright).ge.VOFTOL) then
          print *,"wtleft or wtright invalid"
          stop
         endif

         hx=hx*RR

         if (hx.le.zero) then
          print *,"hx invalid"
          stop
         endif

         im1=i-ii
         jm1=j-jj
         km1=k-kk

         if (dir.eq.0) then
          idx=i
         else if (dir.eq.1) then
          idx=j
         else if ((dir.eq.2).and.(SDIM.eq.3)) then
          idx=k
         else
          print *,"dir invalid tfrmac"
          stop
         endif

         if (operation_flag.eq.1) then ! p^CELL->MAC

          do im=1,ncphys
           local_face(im)=xface(D_DECL(i,j,k),im)
          enddo
           ! newly projected face velocity might be overwritten
           ! with the solid velocity.
          do im_vel=1,nsolveMM_FACE
           local_vel(im_vel)=xvel(D_DECL(i,j,k),im_vel)
          enddo

         else if (operation_flag.eq.2) then !(grd ppot)_MAC,ppot^CELL->MAC,ten.

          do im=1,ncphys
           local_face(im)=xface(D_DECL(i,j,k),im)
          enddo

         else if (operation_flag.eq.9) then ! density: CELL->MAC (1/den)

          do im=1,ncphys
           local_face(im)=xface(D_DECL(i,j,k),im)
          enddo
          local_vel(1)=xvel(D_DECL(i,j,k),1) ! low order 1/rho

         else if ((operation_flag.eq.3).or. & !unew^CELL->MAC
                  (operation_flag.eq.4).or. & !unew^MAC=uSOLID^MAC / uFLUID^MAC
                  (operation_flag.eq.5).or. & !unew^MAC+beta FVISC^CELL->MAC
                  (operation_flag.eq.10).or. & !unew^{CELL,MAC} -> MAC
                  (operation_flag.eq.11)) then !unew^{CELL DIFF,MAC} -> MAC

          do im=1,ncphys
           local_face(im)=xface(D_DECL(i,j,k),im)
          enddo
          do im_vel=1,nsolveMM_FACE
           local_vel(im_vel)=xvel(D_DECL(i,j,k),im_vel)
           if ((operation_flag.eq.5).or. &
               (operation_flag.eq.11)) then
            local_vel_old(im_vel)=xgp(D_DECL(i,j,k),im_vel)
           else if ((operation_flag.eq.3).or. &
                    (operation_flag.eq.4).or. &
                    (operation_flag.eq.10)) then
            local_vel_old(im_vel)=zero
           else
            print *,"operation_flag invalid13"
            stop
           endif
          enddo ! im_vel=1,nsolveMM_FACE

         else if (operation_flag.eq.0) then ! (grad p)_MAC

          ! do nothing

         else if (operation_flag.eq.7) then ! advection

          if (nsolveMM_FACE.eq.1) then

           im_vel=1
           local_vel(im_vel)=xvel(D_DECL(i,j,k),im_vel)

          else
           print *,"nsolveMM_FACE invalid"
           stop
          endif

      
         else if (operation_flag.eq.6) then

          print *,"fort_face_gradients calls sem_cell_to_mac with op=6"
          stop
 
         else
          print *,"operation_flag invalid14"
          stop
         endif

          ! set local_incomp, LSleft, LSright, localLS, xmac
         local_incomp=-1
         if ((operation_flag.eq.2).or. &
             (operation_flag.eq.3).or. &
             (operation_flag.eq.4).or. &
             (operation_flag.eq.5).or. &
             (operation_flag.eq.10).or. & ! u^MAC,CELL -> u^MAC
             (operation_flag.eq.11).or. & ! u^MAC,CELL DIFF -> u^MAC
             (operation_flag.eq.1)) then  ! p^CELL->MAC

          ! levelPC() has piecewise constant BC at coarse/fine borders.
          do im=1,nmat
           LSleft(im)=levelPC(D_DECL(im1,jm1,km1),im)
           LSright(im)=levelPC(D_DECL(i,j,k),im)
           localLS(im)=half*(LSright(im)+LSleft(im))
          enddo

          do dir2=1,SDIM
           xmac(dir2)=xstenMAC(0,dir2)
          enddo

          local_incomp=0
          if (all_incomp.eq.1) then
           local_incomp=1
          else if (all_incomp.eq.0) then
           if ((make_interface_incomp.eq.1).or. &
               (make_interface_incomp.eq.2)) then
            do im=1,nmat
             if (abs(localLS(im)).le.cutoff) then ! cutoff=DXMAXLS
              do im_opp=1,nmat
               if (im.ne.im_opp) then
                if (abs(localLS(im_opp)).le.cutoff) then ! cutoff=DXMAXLS
                 local_incomp=1
                else if (abs(localLS(im_opp)).gt.cutoff) then
                 ! do nothing
                else
                 print *,"localLS invalid"
                 stop
                endif
               else if (im.eq.im_opp) then
                ! do nothing
               else
                print *,"im or im_opp invalid"
                stop
               endif
              enddo ! im_opp=1..nmat 
             else if (abs(localLS(im)).gt.cutoff) then ! cutoff=DXMAXLS
              ! do nothing
             else
              print *,"localLS invalid"
              stop
             endif
            enddo ! im=1..nmat
           else if (make_interface_incomp.eq.0) then
            ! do nothing
           else
            print *,"make_interface_incomp invalid"
            stop
           endif
          else
           print *,"all_incomp invalid"
           stop
          endif 

         else if ((operation_flag.eq.0).or. &
                  ((operation_flag.ge.6).and. &
                   (operation_flag.le.9))) then
          ! do nothing
         else 
          print *,"operation_flag invalid15"
          stop
         endif

          ! advection
          ! The low order fluxes here might be needed by the high order
          ! divergence operator if there is a spectral element next
          ! to a low order element.
         if (operation_flag.eq.7) then 

          if (im_vel.ne.1) then
           print *,"im_vel invalid"
           stop
          endif

          do nc=1,ncphys
           local_face(nc)=zero
          enddo

          at_RZ_face=0
          if (levelrz.eq.0) then
           ! do nothing
          else if (levelrz.eq.1) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if ((dir.eq.0).and. &
               (xstenMAC(0,1).le.VOFTOL*dx(1))) then
            at_RZ_face=1
           endif
          else if (levelrz.eq.3) then
           if ((dir.eq.0).and. &
               (xstenMAC(0,1).le.VOFTOL*dx(1))) then
            at_RZ_face=1
           endif
          else
           print *,"levelrz invalid tfrmac"
           stop
          endif 

          if (at_RZ_face.eq.1) then

           test_velocity_FACE=zero
           idonate=i
           jdonate=j
           kdonate=k
           do im=1,nmat
            LSupwind(im)=levelPC(D_DECL(idonate,jdonate,kdonate),im)
           enddo
           call get_primary_material(LSupwind,nmat,im)
           if ((im.ge.1).and.(im.le.nmat)) then
            ibase=num_state_material*(im-1) 
            denlocal=den(D_DECL(idonate,jdonate,kdonate),ibase+1) 
            templocal=den(D_DECL(idonate,jdonate,kdonate),ibase+2) 
            nc=SDIM+1  ! density
            local_face(nc)=denlocal  ! NONCONSERVATIVE
            nc=SDIM+2  ! temperature
            local_face(nc)=templocal  ! NONCONSERVATIVE
           else
            print *,"im invalid"
            stop
           endif

          else if (at_RZ_face.eq.0) then

           test_velocity_FACE=local_vel(im_vel)

           if (test_velocity_FACE.ge.zero) then
            idonate=im1
            jdonate=jm1
            kdonate=km1
           else
            idonate=i
            jdonate=j
            kdonate=k
           endif
           do im=1,nmat
            LSupwind(im)=levelPC(D_DECL(idonate,jdonate,kdonate),im)
           enddo
           call get_primary_material(LSupwind,nmat,im)
           if ((im.ge.1).and.(im.le.nmat)) then
            ibase=num_state_material*(im-1) 
            do nc=1,ncphys
             denlocal=den(D_DECL(idonate,jdonate,kdonate),ibase+1) 
             templocal=den(D_DECL(idonate,jdonate,kdonate),ibase+2) 

             if ((nc.ge.1).and.(nc.le.SDIM)) then
              velcomp=nc
              local_face(nc)= &
                vel(D_DECL(idonate,jdonate,kdonate),velcomp) !NONCONSERVATIVE
             else if (nc.eq.SDIM+1) then
              local_face(nc)=denlocal   !NONCONSERVATIVE
             else if (nc.eq.SDIM+2) then
              local_face(nc)=templocal  !NONCONSERVATIVE
             else
              print *,"nc invalid"
              stop
             endif

             if ((nc.ge.1).and.(nc.le.SDIM)) then
              local_face(nc)=local_face(nc)*test_velocity_FACE
             else if (nc.eq.SDIM+1) then ! density
              ! do nothing
             else if (nc.eq.SDIM+2) then ! temperature
              ! do nothing
             else
              print *,"nc invalid"
              stop
             endif

            enddo ! nc=1..ncphys
           else
            print *,"im invalid32"
            stop
           endif
          else
           print *,"at_RZ_face invalid"
           stop
          endif

         else if (operation_flag.eq.9) then ! density: CELL->MAC (1/rho)

          if (local_vel(1).le.zero) then
           print *,"low order 1/rho invalid"
           stop
          endif
 
         else if ((operation_flag.eq.3).or. & !u^MAC=u^CELL->MAC
                  (operation_flag.eq.4).or. & !u^MAC=uSOLID^MAC or uFLUID^MAC
                  (operation_flag.eq.5).or. & !u^MAC=u^MAC+beta diff^CELL->MAC
                  (operation_flag.eq.10).or. & !u^MAC=u^{CELL,MAC}->MAC
                  (operation_flag.eq.11)) then !u^MAC=u^{CELL DIFF,MAC}->MAC

          im_vel=1
          face_velocity_override=0

          at_RZ_face=0
          if (levelrz.eq.0) then
           ! do nothing
          else if (levelrz.eq.1) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if ((dir.eq.0).and. &
               (xstenMAC(0,1).le.VOFTOL*dx(1))) then
            at_RZ_face=1
           endif
          else if (levelrz.eq.3) then
           if ((dir.eq.0).and. &
               (xstenMAC(0,1).le.VOFTOL*dx(1))) then
            at_RZ_face=1
           endif
          else
           print *,"levelrz invalid tfrmac"
           stop
          endif 
          if (num_materials_face.ne.1) then
           print *,"num_materials_face invalid"
           stop
          endif

          uedge(im_vel)=zero

          if (at_RZ_face.eq.1) then
           face_velocity_override=1
           uedge(im_vel)=zero
          else if (at_RZ_face.eq.0) then

           fluid_volface=zero
           not_prescribed_volface=zero
           volface=zero

           do side=1,2
            do im=1,nmat
             local_volume=local_face(vofface_index+2*(im-1)+side)
             volface=volface+local_volume
             if (is_prescribed(nmat,im).eq.0) then
              not_prescribed_volface=not_prescribed_volface+local_volume
             else if (is_prescribed(nmat,im).eq.1) then
              ! do nothing
             else
              print *,"is_prescribed(nmat,im) invalid"
              stop
             endif
             if (is_rigid(nmat,im).eq.0) then
              if (is_prescribed(nmat,im).eq.0) then
               fluid_volface=fluid_volface+local_volume
              else
               print *,"is_prescribed(nmat,im) invalid"
               stop
              endif
             else if (is_rigid(nmat,im).eq.1) then
              ! do nothing
             else
              print *,"is_rigid invalid"
              stop
             endif
            enddo ! im=1..nmat
           enddo ! side=1,2

           if (volface.le.zero) then
            print *,"not_prescribed_volface ",not_prescribed_volface
            print *,"fluid_volface ",fluid_volface
            print *,"volface ",volface
            print *,"volface bust tfrmac"
            stop
           endif

           fluid_volface=fluid_volface/volface
           not_prescribed_volface=not_prescribed_volface/volface

           if ((local_face(facecut_index+1).ge.zero).and. &
               (local_face(facecut_index+1).le.half)) then
            fluid_volface=zero
            not_prescribed_volface=zero
           else if ((local_face(facecut_index+1).ge.half).and. &
                    (local_face(facecut_index+1).le.one)) then
            fluid_volface= &
             min(fluid_volface,local_face(facecut_index+1))
            not_prescribed_volface= &
             min(not_prescribed_volface,local_face(facecut_index+1))
           else
            print *,"local_face(facecut_index+1) invalid"
            stop
           endif
         
           ! is_solid_face==1 if:
           !   0.0<=fluid_volface<=VOFTOL_AREAFRAC  or
           !   max(LSleft(im_solid),LSright(im_solid))>=0.0
           ! is_prescribed_face==1 if:
           !   0.0<=not_prescribed_prescribed<=VOFTOL_AREAFRAC  or
           !   max(LSleft(im_prescribed),LSright(im_prescribed))>=0.0
           call fixed_face( &
            nmat, &
            fluid_volface, &
            not_prescribed_volface, &
            LSleft,LSright, &
            is_solid_face, &
            is_prescribed_face, &
            im_solid, &
            im_prescribed, &
            im_solid_valid, &
            im_prescribed_valid, &
            partid_solid, &
            partid_prescribed)
          
           if (is_prescribed_face.eq.1) then
            if (is_solid_face.eq.1) then 
             if (im_prescribed_valid.eq.1) then
              if (im_prescribed.ne.im_solid_map(partid_prescribed+1)+1) then
               print *,"im_solid_map invalid"
               stop
              endif
              if (prescribed_solid_scale(im_prescribed).eq.zero) then
               face_velocity_override=1
              else if ((prescribed_solid_scale(im_prescribed).gt.zero).and. &
                       (prescribed_solid_scale(im_prescribed).le.half)) then
               ! do nothing
              else
               print *,"prescribed_solid_scale(im_prescribed) invalid"
               stop
              endif
             else if (im_prescribed_valid.eq.0) then
              face_velocity_override=1
             else
              print *,"im_prescribed_valid invalid"
              stop
             endif
            else
             print *,"is_solid_face invalid 6 ",is_solid_face
             stop
            endif
           else if (is_prescribed_face.eq.0) then
            ! do nothing
           else
            print *,"is_prescribed_face invalid"
            stop
           endif
 
           if ((is_prescribed_face.eq.1).and. &
               (face_velocity_override.eq.1)) then

            if (is_solid_face.eq.1) then 

             if (im_prescribed_valid.eq.1) then
              if (im_prescribed.ne.im_solid_map(partid_prescribed+1)+1) then
               print *,"im_solid_map invalid"
               stop
              endif
              velcomp=partid_prescribed*SDIM+dir+1 
              uedge(im_vel)=wtright*sol(D_DECL(i,j,k),velcomp)+ &
                            wtleft*sol(D_DECL(im1,jm1,km1),velcomp)
              DEBUG_PRESCRIBED_VEL_TOT=DEBUG_PRESCRIBED_VEL_TOT+uedge(im_vel)
              DEBUG_PRESCRIBED_VEL_DEN=DEBUG_PRESCRIBED_VEL_DEN+one
             else if (im_prescribed_valid.eq.0) then
              uedge(im_vel)=zero
             else
              print *,"im_prescribed_valid invalid"
              stop
             endif
            else
             print *,"is_solid_face invalid 6 ",is_solid_face
             stop
            endif

           else if ((is_prescribed_face.eq.0).or. &
                    (face_velocity_override.eq.0)) then

            test_current_icefacecut= &
                 xface(D_DECL(i,j,k),icefacecut_index+1)
          
            if ((test_current_icefacecut.ge.zero).and. &
                (test_current_icefacecut.le.one)) then

             test_current_icemask=xface(D_DECL(i,j,k),icemask_index+1)

             if ((test_current_icemask.eq.zero).or. &
                 (test_current_icemask.eq.one)) then

              velsum=zero
              mass_sum=zero

              do side=1,2

               partid_check=0
               do im=1,nmat

                DMface=local_face(massface_index+2*(im-1)+side)
                if (added_weight(im).gt.zero) then
                 DMface=DMface*added_weight(im)
                else
                 print *,"added_weight invalid"
                 stop
                endif

                if (side.eq.1) then
                 ic=im1
                 jc=jm1
                 kc=km1
                else if (side.eq.2) then
                 ic=i
                 jc=j
                 kc=k
                else
                 print *,"side invalid"
                 stop
                endif

                if (operation_flag.eq.3) then ! cell -> MAC
                 velcomp=dir+1
                 velmaterial=vel(D_DECL(ic,jc,kc),velcomp)
                else if (operation_flag.eq.4) then ! mac -> MAC
                 velcomp=1
                 velmaterial=local_vel(1)
                else if (operation_flag.eq.5) then ! MAC+=(cell->MAC)
                 velcomp=dir+1
                  ! local_vel=xvel
                  ! local_vel_old=xgp (a copy of xvel)
                 velmaterialMAC=local_vel_old(1)
                 velmaterial=velmaterialMAC+beta*vel(D_DECL(ic,jc,kc),velcomp)
                 if ((beta.ne.-one).and.(beta.ne.one)) then
                  print *,"beta invalid"
                  stop
                 endif
                else if (operation_flag.eq.10) then ! cell,MAC -> MAC
                 if (face_flag.ne.1) then
                  print *,"face_flag invalid"
                  stop
                 endif
                 if (local_incomp.eq.1) then
                  velcomp=1
                  velmaterial=local_vel(1)
                 else if (local_incomp.eq.0) then
                  velcomp=dir+1
                  velmaterial=vel(D_DECL(ic,jc,kc),velcomp)
                 else
                  print *,"local_incomp invalid"
                  stop
                 endif
                else if (operation_flag.eq.11) then ! cell diff,MAC -> MAC
                 if (face_flag.ne.1) then
                  print *,"face_flag invalid"
                  stop
                 endif
                 if (local_incomp.eq.1) then
                  velcomp=1
                  velmaterial=local_vel(1)
                 else if (local_incomp.eq.0) then
                  velcomp=dir+1
                  velmaterial=local_vel_old(1)+vel(D_DECL(ic,jc,kc),velcomp)
                 else
                  print *,"local_incomp invalid"
                  stop
                 endif
                else
                 print *,"operation_flag invalid16"
                 stop
                endif
                if (is_lag_part(nmat,im).eq.1) then
                 partid_check=partid_check+1
                else if (is_lag_part(nmat,im).eq.0) then
                 ! do nothing
                else
                 print *,"is_lag_part(nmat,im) invalid"
                 stop
                endif
                mass_sum=mass_sum+DMface
                velsum=velsum+DMface*velmaterial

               enddo ! im=1..nmat

               if (partid_check.ne.nparts) then
                print *,"partid_check invalid"
                stop
               endif

              enddo ! side=1,2

              if (mass_sum.gt.zero) then

               uedge(im_vel)=velsum/mass_sum

              else
               print *,"mass_sum invalid tfrmac"
               print *,"operation_flag=",operation_flag
               print *,"i,j,k,dir ",i,j,k,dir
               stop
              endif

              if (num_colors.eq.0) then
               if ((project_option.eq.11).or. & !FSI_material_exists 2nd project
                   (project_option.eq.13)) then !FSI_material_exists 1st project
                print *,"project_option invalid"
                stop
               endif
              else if (num_colors.gt.0) then
               ! do nothing
              else
               print *,"num_colors invalid"
               stop
              endif

              if (num_colors.gt.0) then 

                ! type init in FORT_GETTYPEFAB
                typeleft=NINT(typefab(D_DECL(im1,jm1,km1)))
                typeright=NINT(typefab(D_DECL(i,j,k)))
                colorleft=NINT(colorfab(D_DECL(im1,jm1,km1)))
                colorright=NINT(colorfab(D_DECL(i,j,k)))
                if ((typeleft.ge.1).and.(typeleft.le.nmat).and. &
                    (typeright.ge.1).and.(typeright.le.nmat)) then

                 if ((is_ice(nmat,typeleft).eq.1).or. &
                     (is_ice(nmat,typeright).eq.1).or. &
                     (is_FSI_rigid(nmat,typeleft).eq.1).or. &
                     (is_FSI_rigid(nmat,typeright).eq.1)) then

                  if (((is_ice(nmat,typeleft).eq.1).and. &
                       (is_ice(nmat,typeright).eq.1)).or. &
                      ((is_FSI_rigid(nmat,typeleft).eq.1).and. &
                       (is_FSI_rigid(nmat,typeright).eq.1))) then
                   if (LSleft(typeleft).ge.LSright(typeright)) then  
                    typeface=typeleft
                    colorface=colorleft
                   else if (LSleft(typeleft).lt.LSright(typeright)) then
                    typeface=typeright
                    colorface=colorright
                   else
                    print *,"LSleft or LSright bust"
                    stop
                   endif
                  else if ((is_ice(nmat,typeleft).eq.1).or. &
                           (is_FSI_rigid(nmat,typeleft).eq.1)) then
                   typeface=typeleft
                   colorface=colorleft
                  else if ((is_ice(nmat,typeright).eq.1).or. &
                           (is_FSI_rigid(nmat,typeright).eq.1)) then
                   typeface=typeright
                   colorface=colorright
                  else
                   print *,"typeleft or typeright bust"
                   stop
                  endif

                  if ((is_ice(nmat,typeface).eq.1).or. &
                      (is_FSI_rigid(nmat,typeface).eq.1)) then
                   if ((colorface.ge.1).and.(colorface.le.num_colors)) then
                     ! in: GLOBALUTIL.F90
                    call get_rigid_velocity( &
                     FSI_prescribed_flag, &
                     colorface,dir+1,uedge_rigid, &
                     xmac,blob_array, &
                     blob_array_size,num_colors,num_elements_blobclass) 

                    uedge(im_vel)=test_current_icemask*uedge(im_vel)+ &
                            (one-test_current_icemask)*uedge_rigid

                    if (project_option.eq.13) then

                     if (operation_flag.eq.4) then
                      test_local_icefacecut=xcut(D_DECL(i,j,k),1)
                      if ((test_local_icefacecut.ge.zero).and. &
                          (test_local_icefacecut.le.one)) then
                       if (FSI_prescribed_flag.eq.0) then
                        if (test_local_icefacecut.eq.one) then
                         ! do nothing
                        else
                         print *,"test_local_icefacecut invalid"
                         stop
                        endif
                       else if (FSI_prescribed_flag.eq.1) then
                        xcut(D_DECL(i,j,k),1)=test_current_icefacecut
                       else
                        print *,"FSI_prescribed_flag invalid"
                        stop
                       endif
                      else
                       print *,"test_local_icefacecut invalid"
                       stop
                      endif
                     else
                      print *,"operation_flag invalid"
                      stop
                     endif

                    endif ! project_option==13

                   else if (colorface.eq.0) then
                    ! do nothing
                   else
                    print *,"colorface invalid"
                    stop
                   endif 
                  else
                   print *,"is_ice or is_FSI_rigid bust"
                   stop
                  endif

                 else if ((is_ice(nmat,typeleft).eq.0).and. &
                          (is_ice(nmat,typeright).eq.0).and. &
                          (is_FSI_rigid(nmat,typeleft).eq.0).and. &
                          (is_FSI_rigid(nmat,typeright).eq.0)) then
                  ! do nothing
                 else
                  print *,"is_ice or is_FSI_rigid invalid"
                  stop
                 endif

                else
                 print *,"typeleft or typeright invalid"
                 stop
                endif

              else if (num_colors.eq.0) then
               ! do nothing
              else
               print *,"num_colors invalid"
               stop
              endif 

             else
              print *,"test_current_icemask bad: ",test_current_icemask
              stop
             endif

            else
             print *,"test_current_icefacecut bad: ",test_current_icefacecut
             stop
            endif

           else
            print *,"is_prescribed_face or (7) ",is_prescribed_face
            print *,"face_velocity_override invalid ",face_velocity_override
            stop
           endif  

           iboundary=0
           side=0
           if (idx.eq.fablo(dir+1)) then
            iboundary=1
            side=1
           else if (idx.eq.fabhi(dir+1)+1) then
            iboundary=1
            side=2
           else if ((idx.gt.fablo(dir+1)).and. &
                    (idx.lt.fabhi(dir+1)+1)) then
            ! do nothing
           else
            print *,"idx invalid"
            stop 
           endif
           if ((iboundary.eq.1).and. &
               ((side.eq.1).or.(side.eq.2))) then
            if (velbc_in(dir+1,side,dir+1).eq.REFLECT_ODD) then
             face_velocity_override=1
             uedge(im_vel)=zero
            else if (velbc_in(dir+1,side,dir+1).eq.EXT_DIR) then
             face_velocity_override=1
             call velbc_override(time,dir+1,side,dir+1, &
              uedge(im_vel), &
              xstenMAC,nhalf,dx,bfact)
            endif
           else if ((iboundary.eq.0).and.(side.eq.0)) then
            ! do nothing
           else 
            print *,"iboundary or side invalid"
            stop
           endif

          else
           print *,"at_RZ_face invalid"
           stop
          endif 

          ! local_vel initialized above with the current MAC velocity
          ! contents.
         else if (operation_flag.eq.1) then ! p^CELL->MAC

          im_vel=1
          pplus(im_vel)=pres(D_DECL(i,j,k),im_vel)
          pminus(im_vel)=pres(D_DECL(im1,jm1,km1),im_vel)

           ! mask=1 if not covered
           ! mask=0 if covered
          mask_covered(1)=NINT(maskcoef(D_DECL(im1,jm1,km1)))
          mask_covered(2)=NINT(maskcoef(D_DECL(i,j,k)))

            ! mask=0 coarse/fine
            ! mask=1 fine/fine
          mask_coarsefine(1)=NINT(mask(D_DECL(im1,jm1,km1)))
          mask_coarsefine(2)=NINT(mask(D_DECL(i,j,k)))

          at_reflect_wall=0
          at_wall=0
          at_ext_wall=0
          at_coarse_fine_wallF=0
          at_coarse_fine_wallC=0

          use_face_pres=three ! use both gp and div(up)
          solid_velocity=local_face(facevel_index+1)

          face_velocity_override=0

          if (face_flag.eq.1) then

           if (local_incomp.eq.1) then
            use_face_pres=one  ! div(up) ok, but not gp
           else if (local_incomp.eq.0) then
            use_face_pres=one  ! div(up) ok, but not gp
           else
            print *,"local_incomp invalid"
            stop
           endif
 
          else if (face_flag.eq.0) then
           use_face_pres=three  ! div(up) and gp ok.
          else
           print *,"face_flag invalid"
           stop
          endif

          do side=1,2

           if (((side.eq.1).and.(idx.eq.fablo(dir+1))).or. &
               ((side.eq.2).and.(idx.eq.fabhi(dir+1)+1))) then

            if ((presbc_in(dir+1,side,1).eq.REFLECT_EVEN).or. &
                (presbc_in(dir+1,side,1).eq.FOEXTRAP)) then
             at_wall=1
            else if ((presbc_in(dir+1,side,1).eq.INT_DIR).or. &
                     (presbc_in(dir+1,side,1).eq.EXT_DIR)) then
             ! do nothing
            else
             print *,"presbc_in invalid"
             stop
            endif

            if (velbc_in(dir+1,side,dir+1).eq.REFLECT_ODD) then
             at_reflect_wall=side
             face_velocity_override=1
             solid_velocity=zero
            else if (velbc_in(dir+1,side,dir+1).eq.EXT_DIR) then
             face_velocity_override=1
             use_face_pres=one ! do not use gp 
             at_ext_wall=side
            else if (velbc_in(dir+1,side,dir+1).eq.INT_DIR) then
             if (mask_coarsefine(side).eq.0) then  
              at_coarse_fine_wallF=side
              if ((project_option.eq.1).or. &
                  (COARSE_FINE_VELAVG.eq.1)) then
               use_face_pres=one ! do not use gp 
              else if (((project_option.eq.0).or. &
                        (project_option.eq.10).or. &
                        (project_option.eq.13).or. &  !FSI_material_exists 1st
                        (project_option.eq.11)).and. &!FSI_material_exists 2nd
                       (COARSE_FINE_VELAVG.eq.0)) then
               ! do nothing
              else
               print *,"project_option invalid edge pressure 2"
               stop
              endif
             else if (mask_coarsefine(side).eq.1) then
              ! do nothing
             else
              print *,"mask_coarsefine invalid"
              stop
             endif
            endif ! int_dir case
           endif ! idx=fablo or idx=fabhi+1

          enddo ! side=1..2

           ! sanity check
          if (levelrz.eq.0) then
           ! do nothing
          else if (levelrz.eq.1) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if ((xstenMAC(0,1).le.VOFTOL*dx(1)).and. &
               (dir.eq.0)) then
            if (at_reflect_wall.ne.1) then
             print *,"at_reflect_wall fails sanity check"
             stop
            endif
           endif
          else if (levelrz.eq.3) then
           ! do nothing
          else
           print *,"levelrz invalid grad potential 2"
           stop
          endif

          if ((local_face(facecut_index+1).ge.zero).and. &
              (local_face(facecut_index+1).le.half)) then
           AFACE=zero
          else if ((local_face(facecut_index+1).ge.half).and. &
                   (local_face(facecut_index+1).le.one)) then
           AFACE=local_face(facecut_index+1)*local_face(icefacecut_index+1)
          else
           print *,"local_face(facecut_index+1) invalid"
           stop
          endif

          if ((AFACE.ge.zero).and. &
              (AFACE.le.half)) then
           if (at_reflect_wall.eq.0) then
            use_face_pres=one ! do not use gp
           else if ((at_reflect_wall.eq.1).or. &
                    (at_reflect_wall.eq.2)) then
            ! do nothing
           else
            print *,"at_reflect_wall invalid"
            stop
           endif
          else if ((AFACE.ge.half).and.(AFACE.le.one)) then
           ! do nothing
          else
           print *,"AFACE invalid"
           stop 
          endif

          do im=1,nmat
           if ((is_ice(nmat,im).eq.1).or. &  
               (is_FSI_rigid(nmat,im).eq.1).or. &
               (CTML_FSI_mat(nmat,im).eq.1)) then  
            if ((levelPC(D_DECL(i,j,k),im).ge.zero).or. &
                (levelPC(D_DECL(im1,jm1,km1),im).ge.zero)) then
             use_face_pres=zero ! do not use gp or div(up)
            endif
           else if ((FSI_flag(im).eq.0).or. &
                    (FSI_flag(im).eq.7)) then ! regular fluid material
            ! do nothing
           else if (is_prescribed(nmat,im).eq.1) then
            ! do nothing (use_face_pres modified above for this case)
           else
            print *,"FSI_flag invalid"
            stop
           endif
          enddo ! im=1..nmat

          if (local_face(icemask_index+1).eq.zero) then
           use_face_pres=zero ! do not use gp or div(up)
          else if (local_face(icemask_index+1).eq.one) then
           ! do nothing
          else
           print *,"icemask invalid in FORT_CELL_TO_MAC"
           print *,"This is the p^CELL->MAC operation"
           print *,"operation_flag (=1) = ",operation_flag
           print *,"icemask_index= ",icemask_index
           print *,"icemask= ",local_face(icemask_index+1)
           stop
          endif

          if ((mask_covered(1).eq.0).and. &
              (mask_covered(2).eq.1)) then
           at_coarse_fine_wallC=1
          else if ((mask_covered(1).eq.1).and. &
                   (mask_covered(2).eq.0)) then
           at_coarse_fine_wallC=2
          else if ((mask_covered(1).eq.0).and. &
                   (mask_covered(2).eq.0)) then
           ! do nothing
          else if ((mask_covered(1).eq.1).and. &
                   (mask_covered(2).eq.1)) then
           ! do nothing
          else
           print *,"mask_covered invalid"
           stop
          endif

          if ((project_option.eq.1).or.(COARSE_FINE_VELAVG.eq.1)) then
            ! at least 1 side is covered
           if ((mask_covered(1).eq.0).or. &
               (mask_covered(2).eq.0)) then
            use_face_pres=one ! do not use gp 
           else if ((mask_covered(1).eq.1).and. &
                    (mask_covered(2).eq.1)) then
            ! do nothing
           else
            print *,"mask_covered invalid"
            stop
           endif
          else if (((project_option.eq.0).or. &
                    (project_option.eq.10).or. &
                    (project_option.eq.13).or. &  !FSI_material_exists 1st
                    (project_option.eq.11)).and. &!FSI_material_exists 2nd
                   (COARSE_FINE_VELAVG.eq.0)) then
           ! both sides are covered
           if ((mask_covered(1).eq.0).and. &
               (mask_covered(2).eq.0)) then
            use_face_pres=zero  ! do not use gp or div(up)
           else if ((mask_covered(1).eq.1).or. &
                    (mask_covered(2).eq.1)) then
            ! do nothing
           else
            print *,"mask_covered invalid"
            stop
           endif
          else
           print *,"project_option invalid FORT_CELL_TO_MAC 4"
           stop
          endif

          ! side=1 is right half of the cell that is to the left of the face.
          ! side=2 is left half of the cell that is to the right of the face.
          massface=zero
          volface=zero
          denface=zero

          do side=1,2
           mass(side)=zero
           vol_local(side)=zero
           do im=1,nmat
            mass(side)=mass(side)+ &
             local_face(massface_index+2*(im-1)+side)
            vol_local(side)=vol_local(side)+ &
             local_face(vofface_index+2*(im-1)+side)
           enddo 
           if (mass(side).lt.zero) then
            print *,"mass(side) invalid"
            stop
           endif
           if (vol_local(side).lt.zero) then
            print *,"vol_local(side) invalid"
            stop
           else if (vol_local(side).eq.zero) then
            den_local(side)=zero
           else
            den_local(side)=mass(side)/vol_local(side)
           endif
           massface=massface+mass(side)
           volface=volface+vol_local(side)
           denface=denface+den_local(side)
          enddo  ! side=1,2

           ! 1=use_face_pres  2=coarse fine flag 2=face pressure
          do im=1,2+nsolveMM_FACE
           plocal(im)=zero
          enddo

          fluid_volface=one
          not_prescribed_volface=one
          call fixed_face( &
           nmat, &
           fluid_volface, &
           not_prescribed_volface, &
           LSleft,LSright, &
           is_solid_face, &
           is_prescribed_face, &
           im_solid, &
           im_prescribed, &
           im_solid_valid, &
           im_prescribed_valid, &
           partid_solid, &
           partid_prescribed) 

          if (is_prescribed_face.eq.1) then

           if (is_solid_face.eq.1) then
            if (im_prescribed_valid.eq.1) then
             if ((im_prescribed.ge.1).and. &
                 (im_prescribed.le.nmat)) then
              if (im_solid_map(partid_prescribed+1)+1.eq.im_prescribed) then
               use_face_pres=zero ! do not use gp or div(up)
               if (prescribed_solid_scale(im_prescribed).eq.zero) then
                face_velocity_override=1
               else if (prescribed_solid_scale(im_prescribed).gt.zero) then
                ! do nothing
               else
                print *,"prescribed_solid_scale(im_prescribed) invalid"
                stop
               endif
              else
               print *,"im_solid_map(partid_prescribed+1)+1 invalid"
               stop
              endif
             else
              print *,"im_prescribed invalid"
              stop
             endif
            else
             print *,"im_prescribed_valid corrupt"
             stop
            endif
           else
            print *,"is_solid_face invalid 8 ",is_solid_face
            stop
           endif

          else if (is_prescribed_face.eq.0) then
           ! do nothing
          else
           print *,"is_prescribed_face invalid"
           stop
          endif

          if (at_reflect_wall.eq.0) then
           if (denface.le.zero) then
            print *,"denface must be positive"
            stop
           endif
           im=im_vel
            ! face pressure
           plocal(2+im_vel)= &
            (den_local(1)*pplus(im)+den_local(2)*pminus(im))/denface
          else if (at_reflect_wall.eq.1) then ! left wall

           im=im_vel
            ! face pressure
           plocal(2+im_vel)=pplus(im)

          else if (at_reflect_wall.eq.2) then ! right wall

           im=im_vel
            ! face pressure
           plocal(2+im_vel)=pminus(im)

          else
           print *,"at reflect wall invalid"
           stop
          endif
          plocal(1)=use_face_pres
          plocal(2)=at_coarse_fine_wallF+at_coarse_fine_wallC*10

          if (face_velocity_override.eq.1) then
           local_vel(im_vel)=solid_velocity
          else if (face_velocity_override.eq.0) then
           if (at_wall.eq.0) then
            ! do nothing
           else
            print *,"at_wall or face_velocity_override invalid"
            stop
           endif
          else
           print *,"face_velocity_override invalid"
           stop
          endif

         else if (operation_flag.eq.2) then !(grd pot)_MAC,pot^CELL->MAC,ten.

           ! hydrostatic pressure
          pplus(1)=pres(D_DECL(i,j,k),1)
          pminus(1)=pres(D_DECL(im1,jm1,km1),1)

           ! hydrostatic density
          dplus=den(D_DECL(i,j,k),1)
          dminus=den(D_DECL(im1,jm1,km1),1)
          if ((dplus.le.zero).or. &
              (dminus.le.zero)) then
           print *,"hydrostatic density must be positive"
           stop
          endif

          at_reflect_wall=0
          at_wall=0

          do side=1,2

           if (((side.eq.1).and.(idx.eq.fablo(dir+1))).or. &
               ((side.eq.2).and.(idx.eq.fabhi(dir+1)+1))) then

            if ((presbc_in(dir+1,side,1).eq.REFLECT_EVEN).or. &
                (presbc_in(dir+1,side,1).eq.FOEXTRAP)) then
             at_wall=1
            else if ((presbc_in(dir+1,side,1).eq.INT_DIR).or. &
                     (presbc_in(dir+1,side,1).eq.EXT_DIR)) then
             ! do nothing
            else
             print *,"presbc_in invalid"
             stop
            endif

            if (velbc_in(dir+1,side,dir+1).eq.REFLECT_ODD) then
             at_reflect_wall=side
            else if (velbc_in(dir+1,side,dir+1).eq.EXT_DIR) then
             ! do nothing
            else if (velbc_in(dir+1,side,dir+1).eq.INT_DIR) then
             ! do nothing
            endif ! int_dir case
           endif ! idx=fablo or idx=fabhi+1

          enddo ! side

           ! sanity check
          if (levelrz.eq.0) then
           ! do nothing
          else if (levelrz.eq.1) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if ((xstenMAC(0,1).le.VOFTOL*dx(1)).and. &
               (dir.eq.0)) then
            if (at_reflect_wall.ne.1) then
             print *,"at_reflect_wall fails sanity check"
             stop
            endif
           endif
          else if (levelrz.eq.3) then
           ! do nothing
          else
           print *,"levelrz invalid grad potential 2"
           stop
          endif

          AFACE=local_face(facecut_index+1)
          if ((AFACE.ge.zero).and.(AFACE.le.half)) then
           AFACE=zero
          else if ((AFACE.ge.half).and.(AFACE.le.one)) then
           ! do nothing
          else
           print *,"AFACE invalid"
           stop 
          endif

          AFACE_ICE=local_face(icefacecut_index+1)
          if ((AFACE_ICE.ge.zero).and.(AFACE_ICE.le.one)) then
           ! do nothing
          else
           print *,"AFACE_ICE invalid"
           stop
          endif

          ! side=1 is right half of the cell that is to the left of the face.
          ! side=2 is left half of the cell that is to the right of the face.
          massface=zero
          volface=zero
          denface=zero

          do side=1,2
           mass(side)=zero
           vol_local(side)=zero
           do im=1,nmat
            mass(side)=mass(side)+ &
             local_face(massface_index+2*(im-1)+side)
            vol_local(side)=vol_local(side)+ &
             local_face(vofface_index+2*(im-1)+side)
           enddo 
           if (mass(side).lt.zero) then
            print *,"mass(side) invalid"
            stop
           endif
           if (vol_local(side).lt.zero) then
            print *,"vol_local(side) invalid"
            stop
           else if (vol_local(side).eq.zero) then
            den_local(side)=zero
           else
            den_local(side)=mass(side)/vol_local(side)
           endif
           massface=massface+mass(side)
           volface=volface+vol_local(side)
           denface=denface+den_local(side)
          enddo  ! side=1,2

          pgrad_gravity=zero
          gradh=zero
          pgrad_tension=zero
           ! (1) hydrostatic pressure on the face.
           ! (2) surface tension left of face
           ! (3) surface tension right of face
          do im=1,3
           plocal(im)=zero
          enddo

          call fixed_face( &
           nmat, &
           AFACE, &
           AFACE, &
           LSleft,LSright, &
           is_solid_face, &
           is_prescribed_face, &
           im_solid, &
           im_prescribed, &
           im_solid_valid, &
           im_prescribed_valid, &
           partid_solid, &
           partid_prescribed) 

          if (solvability_projection.eq.1) then
           ! do nothing
          else if (solvability_projection.eq.0) then
           ! do nothing
          else
           print *,"solvability_projection invalid"
           stop
          endif
 
           ! at_wall==1 if FOEXTRAP or REFLECT_EVEN BC for pressure.
           ! gradh represents (H_{i}-H_{i-1})
          if (at_wall.eq.1) then
           ! do nothing, gradh=0 on a wall
          else if (at_reflect_wall.eq.1) then
           ! do nothing, gradh=0 at a reflecting wall 
          else if (at_reflect_wall.eq.2) then
           ! do nothing, gradh=0 at a reflecting wall 
          else if ((at_reflect_wall.eq.0).and. &
                   (at_wall.eq.0)) then

           ! gradh=0 if FSI_flag(im) or FSI_flag(im_opp) = 1,2,4
           if (is_solid_face.eq.1) then
            gradh=zero
           else if (is_solid_face.eq.0) then
            call fluid_interface_tension(LSleft,LSright,gradh, &
             im_opp,im,nmat,nten)
           else
            print *,"is_solid_face invalid 9 ",is_solid_face
            stop
           endif

           if (gradh.ne.zero) then

            do im_heat=1,nmat
             tcomp=(im_heat-1)*num_state_material+2
             mgoni_temp(im_heat)=half*(mgoni(D_DECL(i,j,k),tcomp)+ &
              mgoni(D_DECL(im1,jm1,km1),tcomp))
            enddo ! im_heat

            call get_user_tension(xmac,time, &
             fort_tension,user_tension, &
             mgoni_temp,nmat,nten,3)

            call get_iten(im,im_opp,iten,nmat)
            call get_scaled_tension(user_tension(iten),tension_scaled)
            call get_scaled_pforce(pforce_scaled)

             ! in: FORT_CELL_TO_MAC, operation_flag=2,
             !     surface tension on MAC grid ...
            if (conservative_tension_force.eq.0) then
             local_tension_force=tension_scaled*local_face(curv_index+1)
            else if (conservative_tension_force.eq.1) then
             local_tension_force=zero
            else
             print *,"conservative_tension_force invalid"
             stop
            endif

            ! side=1 is right half of cell that is to left of face.
            ! side=2 is left half of cell that is to right of face.
            ! the pressure gradient force for cell i will be
            ! (plocal(i+1,2)-plocal(i,3))/dencell
            side=1
            plocal(side+1)= &
             (local_tension_force+ &
              pforce_scaled*local_face(pforce_index+1))*gradh* &
             den_local(side)/denface
            side=2
            plocal(side+1)=- &
             (local_tension_force+ &
              pforce_scaled*local_face(pforce_index+1))*gradh* &
             den_local(side)/denface

             ! pgrad_tension is added to pgrad_gravity at the very end.
            pgrad_tension=-(local_tension_force+ &
                    pforce_scaled*local_face(pforce_index+1))*gradh
            pgrad_tension=dt*pgrad_tension/hx
            if ((local_face(facecut_index+1).ge.zero).and. &
                (local_face(facecut_index+1).le.half)) then
             pgrad_tension=zero
            else if ((local_face(facecut_index+1).ge.half).and. &
                     (local_face(facecut_index+1).le.one)) then
             pgrad_tension=pgrad_tension* &
              local_face(facecut_index+1)* &
              local_face(faceden_index+1)
            else
             print *,"local_face(facecut_index+1) invalid"
             stop
            endif

           else if (gradh.eq.zero) then
             ! do nothing
           else
            print *,"gradh bust"
            stop
           endif 
          else
           print *,"at_reflect_wall or at_wall invalid"
           stop
          endif

           ! p=dt( -|g| z + (1/2)Omega^2 r^2 )
           ! force=grad p=dt( -|g| z^hat + Omega^2 r r^hat )
          cutedge=half*(dplus+dminus)
          if (cutedge.le.zero) then
           print *,"cutedge invalid"
           stop
          endif

           ! plocal(1) is the hydrostatic pressure on the MAC grid.        
          if (at_reflect_wall.eq.0) then
           plocal(1)=(pplus(1)*dminus+pminus(1)*dplus)/(two*cutedge)
          else if (at_reflect_wall.eq.1) then ! left wall
             ! limit as dminus->infinity
           plocal(1)=pplus(1)
          else if (at_reflect_wall.eq.2) then ! right wall
             ! limit as dplus->infinity
           plocal(1)=pminus(1)
          else
           print *,"at reflect wall invalid"
           stop
          endif

           ! hydrostatic pressure gradient on the MAC grid
          pgrad_gravity=(pplus(1)-pminus(1))/(hx*cutedge)

          ! -dt k (grad p)_MAC (energyflag=0)
          ! (grad p)_MAC (energyflag=2)
         else if (operation_flag.eq.0) then 
        
           ! xcut=(*localMF[FACE_WEIGHT_MF+dir])[mfi] 
          if ((project_option.eq.0).or. &
              (project_option.eq.1).or. &
              (project_option.eq.10).or. &
              (project_option.eq.11).or. & !FSI_material_exists 2nd project
              (project_option.eq.13).or. & !FSI_material_exists 1st project
              (project_option.eq.12).or. & !pressure extrapolation
              (project_option.eq.2).or. &  ! thermal diffusion
              (project_option.ge.3).or. &  ! viscosity
              ((project_option.ge.100).and. &
               (project_option.lt.100+num_species_var))) then
           cutedge=xcut(D_DECL(i,j,k),1)  ! e.g. A/rho
          else
           print *,"project_option invalid"
           stop
          endif

          do im_vel=1,nsolveMM_FACE
 
           if (num_materials_face.eq.1) then
            im=im_vel
           else if (num_materials_face.eq.nmat) then
            im=im_vel
            if (im_vel.gt.nmat) then
             im=im_vel-nmat
            endif
           else
            print *,"num_materials_face invalid"
            stop
           endif
           pplus(im)=pres(D_DECL(i,j,k),im)
           pminus(im)=pres(D_DECL(im1,jm1,km1),im)

            ! regular solver or SDC viscosity or thermal flux.
           if (energyflag.eq.0) then  
            pgrad(im_vel)=-dt*cutedge*(pplus(im)-pminus(im))/hx

            ! for SDC pressure gradient
           else if (energyflag.eq.2) then  
            pgrad(im_vel)=(pplus(im)-pminus(im))/hx

           else
            print *,"energyflag invalid"
            stop
           endif

          enddo ! im_vel=1..nsolveMM_FACE

         else
          print *,"operation_flag invalid17"
          stop
         endif

         if (operation_flag.eq.7) then ! advection

          if (ncphys.ne.SDIM+num_state_base) then
           print *,"ncphys invalid"
           stop
          endif

          do nc=1,ncphys
           xface(D_DECL(i,j,k),nc)=local_face(nc)
          enddo ! nc

         else if (operation_flag.eq.9) then ! density: CELL->MAC

          ! do nothing (low order already done in INIT_PHYSICS_VARS)

         else if ((operation_flag.eq.3).or. & !u^CELL->MAC
                  (operation_flag.eq.4).or. & !u^MAC=uSOLID^MAC or uFLUID^MAC
                  (operation_flag.eq.5).or. & !u^MAC+beta diff_ref^CELL->MAC
                  (operation_flag.eq.10).or. & ! u^CELL,MAC -> MAC
                  (operation_flag.eq.11)) then ! u^CELL DIFF,MAC -> MAC

          im_vel=1
          xvel(D_DECL(i,j,k),im_vel)=uedge(im_vel)
  
         else if (operation_flag.eq.0) then ! (grad p)_MAC

          do im_vel=1,nsolveMM_FACE
           xgp(D_DECL(i,j,k),im_vel)=pgrad(im_vel)
          enddo

         else if (operation_flag.eq.2) then ! potential grad+surface tension
       
          pgrad_gravity=pgrad_gravity+pgrad_tension

          xgp(D_DECL(i,j,k),1)=pgrad_gravity
          do im=1,3
           xp(D_DECL(i,j,k),im)=plocal(im)
          enddo

         else if (operation_flag.eq.1) then !p^CELL->MAC

          do im=1,2+nsolveMM_FACE
           xp(D_DECL(i,j,k),im)=plocal(im)
          enddo
          im_vel=1
          xvel(D_DECL(i,j,k),im_vel)=local_vel(im_vel)

         else
          print *,"operation_flag invalid18"
          stop
         endif

        enddo
        enddo
        enddo ! i,j,k (MAC grid, zero ghost cells)

        if (DEBUG_PRESCRIBED.eq.1) then
         if (DEBUG_PRESCRIBED_VEL_DEN.gt.zero) then
          call FLUSH(6) ! 6=screen
          print *,"DEBUG_PRESCRIBED ",DEBUG_PRESCRIBED
          print *,"DEBUG_PRESCRIBED_VEL_DEN ",DEBUG_PRESCRIBED_VEL_DEN
          print *,"DEBUG_PRESCRIBED_VEL_TOT ",DEBUG_PRESCRIBED_VEL_TOT
          print *,"tileloop,spectral_loop ",tileloop,spectral_loop
          print *,"dir ",dir
          print *,"project_option=",project_option
          print *,"operation_flag ",operation_flag
          call FLUSH(6) ! 6=screen
         endif
        endif

       else if (spectral_loop.eq.1) then
        ! do nothing
       else
        print *,"spectral_loop invalid"
        stop
       endif

       ! second: high order gradients
      else if (tileloop.eq.1) then

        ! operation_flag=4  unew^MAC=uSOLID^MAC or uFLUID^MAC
       if (operation_flag.ne.4) then

        if ((enable_spectral.eq.1).or. &
            (enable_spectral.eq.2)) then

         if (bfact.ge.2) then

          if (1.eq.0) then
           print *,"in cell to mac"
           print *,"time,dt ",time,dt
           print *,"dir:",dir
           print *,"operation_flag= ",operation_flag
          endif

          if ((dir.lt.0).or.(dir.ge.SDIM)) then
           print *,"dir invalid cell to mac 2"
           stop
          endif

          call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
          do i=growlo(1),growhi(1)
          do j=growlo(2),growhi(2)
          do k=growlo(3),growhi(3)

           call strip_status(i,j,k,bfact,stripstat)

           if (stripstat.eq.1) then

            local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))
            maskcov=NINT(maskcoef(D_DECL(i,j,k)))

            ok_to_HO_interp=0

             ! local_maskSEM==0 for rigid materials or ice
            if ((local_maskSEM.ge.1).and. &
                (local_maskSEM.le.nmat).and. &
                (maskcov.eq.1)) then

             ok_to_HO_interp=1

              ! operation_flag=9  density cell -> MAC
             if (operation_flag.eq.9) then

              if (override_density(local_maskSEM).eq.0) then
               ! do nothing

               ! rho=rho(T,z)
              else if (override_density(local_maskSEM).eq.1) then
               ok_to_HO_interp=0

               ! temperature dependent buoyancy source term.
              else if (override_density(local_maskSEM).eq.2) then
               ok_to_HO_interp=0
              else
               print *,"override_density invalid"
               stop
              endif
              imattype=fort_material_type(local_maskSEM)
              if (imattype.eq.0) then ! incompressible
               ok_to_HO_interp=0
              else if (imattype.eq.999) then ! rigid material
               ok_to_HO_interp=0
              else if ((imattype.gt.0).and. &
                       (imattype.le.fort_max_num_eos)) then
               ! do nothing
              else
               print *,"imattype invalid"
               stop
              endif
             else if ((operation_flag.ge.0).and. &
                      (operation_flag.le.8)) then
              ! do nothing
             else if ((operation_flag.eq.10).or. &
                      (operation_flag.eq.11)) then
              ! do nothing
             else
              print *,"operation_flag invalid19"
              stop
             endif

             if (ok_to_HO_interp.eq.1) then

              call elementbox(i,j,k,bfact,dir,elemlo,elemhi)
              do ielem=elemlo(1),elemhi(1)
              do jelem=elemlo(2),elemhi(2)
              do kelem=elemlo(3),elemhi(3)

               update_right_flux=0

               if (operation_flag.eq.0) then  ! pressure gradient
                if (num_materials_face.eq.1) then
                 scomp=1
                 dcomp=1
                else if (num_materials_face.eq.nmat) then
                 scomp=local_maskSEM
                 dcomp=local_maskSEM
                 update_right_flux=1
                else
                 print *,"num_materials_face invalid"
                 stop
                endif
                ncomp_dest=1
                ncomp_source=1
                scomp_bc=1
               else if (operation_flag.eq.1) then ! pressure CELL->MAC

                scomp=1
                dcomp=3

                ncomp_dest=1
                ncomp_source=1
                scomp_bc=1

               !potential CELL->MAC,grad ppot/den 
               else if (operation_flag.eq.2) then 
                scomp=1
                dcomp=1
                ncomp_dest=1
                ncomp_source=1
                scomp_bc=1

               else if (operation_flag.eq.9) then ! den CELL->MAC

                scomp=(local_maskSEM-1)*num_state_material+1
                dcomp=1
                scomp_bc=1
                ncomp_dest=1
                ncomp_source=nmat*num_state_material

               else if (operation_flag.eq.3) then ! vel CELL->MAC

                scomp=dir+1
                dcomp=1
                ncomp_dest=1
                ncomp_source=1
                scomp_bc=dir+1
               else if (operation_flag.eq.5) then ! u^MAC=u^MAC+beta F^CELL->MAC
                scomp=dir+1
                dcomp=1
                ncomp_dest=1
                ncomp_source=1
                scomp_bc=dir+1
               else if (operation_flag.eq.10) then ! vel CELL,MAC -> MAC
                scomp=dir+1
                dcomp=1
                ncomp_dest=1
                ncomp_source=1
                scomp_bc=dir+1
               else if (operation_flag.eq.11) then ! vel CELL DIFF,MAC -> MAC
                scomp=dir+1
                dcomp=1
                ncomp_dest=1
                ncomp_source=1
                scomp_bc=dir+1
               else if (operation_flag.eq.7) then ! advection
                scomp=1
                dcomp=1
                ncomp_dest=ncphys
                ncomp_source=SDIM*num_materials_face
                scomp_bc=1
               else
                print *,"operation_flag invalid20"
                stop
               endif
    
               call SEM_CELL_TO_MAC( &
                conservative_div_uu, &
                ncomp_xp, &
                simple_AMR_BC_flag, &
                nsolveMM_FACE, &
                num_materials_face, &
                level, &
                finest_level, &
                nmat, &
                operation_flag, & 
                energyflag, &
                temperature_primitive_variable, &
                project_option, &
                SEM_upwind, &
                SEM_advection_algorithm, &
                beta, &
                visc_coef, &
                time, &
                dt, &
                ielem,jelem,kelem, &
                tilelo,tilehi, &
                fablo,fabhi, &
                xlo,dx,dir+1, &
                bfact,bfact_c,bfact_f, &
                presbc_in, &
                velbc_in, &
                scomp, &
                scomp_bc, &
                dcomp, &
                update_right_flux, &
                ncomp_dest, &
                ncomp_source, &
                ncomp_xgp, &
                ncphys, &
                spectral_loop, &
                ncfluxreg, &
                semflux,DIMS(semflux), &
                mask,DIMS(mask), & !mask=1.0 at interior fine bc ghost cells
                maskcoef,DIMS(maskcoef), & ! 1=not cov. or outside domain  
                vel,DIMS(vel), &
                pres,DIMS(pres), &
                den,DIMS(den), &
                xface,DIMS(xface), &
                xgp,DIMS(xgp), & ! holds Umac_old if operation_flag==5 or 11.
                xcut,DIMS(xcut), &   ! coeff*areafrac
                xp,DIMS(xp), &
                xvel,DIMS(xvel), &
                maskSEM,DIMS(maskSEM))

              enddo 
              enddo 
              enddo  ! ielem,jelem,kelem

             else if (ok_to_HO_interp.eq.0) then
              ! do nothing
             else
              print *,"ok_to_HO_interp invalid"
              stop
             endif

            else if ((local_maskSEM.eq.0).or. &
                     (maskcov.eq.0)) then

             ! do nothing
   
            else
             print *,"local_maskSEM or maskcov invalid"
             stop
            endif 
           
           else if (stripstat.eq.0) then
            ! do nothing
           else
            print *,"stripstat invalid"
            stop
           endif

          enddo
          enddo
          enddo ! i,j,k

         else if (bfact.eq.1) then
          ! do nothing
         else
          print *,"bfact invalid100"
          stop
         endif

        else if ((enable_spectral.eq.0).or. &
                 (enable_spectral.eq.3)) then
         ! do nothing
        else
         print *,"enable_spectral invalid"
         stop
        endif

       else if (operation_flag.eq.4) then
        ! do nothing
       else
        print *,"operation_flag invalid"
        stop
       endif

      else
       print *,"tileloop invalid"
       stop
      endif

      return
      end subroutine FORT_CELL_TO_MAC


      subroutine FORT_GFMUPDATE( &
       tid, &
       dirmac, &
       gfm_error, &
       velbc_in, &
       slab_step, &
       dt, &
       time, &
       xlo,dx, &
       masknbr,DIMS(masknbr), & ! 1=fine/fine  0=coarse/fine
       maskcov,DIMS(maskcov), & ! 1=not cov. or outside domain  0=covered
       levelPC,DIMS(levelPC), &
       sol,DIMS(sol), &
       recon,DIMS(recon), &  
       xvel,DIMS(xvel), &
       xsrc,DIMS(xsrc), &
       ysrc,DIMS(ysrc), &
       zsrc,DIMS(zsrc), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact,bfact_c,bfact_f, &
       level,finest_level, &
       rz_flag, &
       domlo,domhi, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       prescribed_solid_scale, &
       prescribed_solid_method, & !0=stair 1=2nd order 2=mod stair 3=mod normal
       project_option)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probf90_module
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T tid
      INTEGER_T dirmac
      INTEGER_T nmat
      INTEGER_T nparts
      INTEGER_T nparts_def
      INTEGER_T im_solid_map(nparts_def)

      REAL_T gfm_error
      REAL_T prescribed_solid_scale(nmat)
      INTEGER_T prescribed_solid_method
      INTEGER_T slab_step
      INTEGER_T level
      INTEGER_T finest_level
      REAL_T dt,time
      REAL_T xlo(SDIM),dx(SDIM)
      INTEGER_T DIMDEC(masknbr)
      INTEGER_T DIMDEC(maskcov)
      INTEGER_T DIMDEC(recon)
      INTEGER_T DIMDEC(xvel)
      INTEGER_T DIMDEC(xsrc)
      INTEGER_T DIMDEC(ysrc)
      INTEGER_T DIMDEC(zsrc)
      INTEGER_T DIMDEC(levelPC)
      INTEGER_T DIMDEC(sol)

      INTEGER_T velbc_in(SDIM,2,num_materials_vel*SDIM)
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact,bfact_c,bfact_f
      INTEGER_T rz_flag
      INTEGER_T domlo(SDIM),domhi(SDIM)
      INTEGER_T project_option

      REAL_T masknbr(DIMV(masknbr))
      REAL_T maskcov(DIMV(maskcov))

      REAL_T levelPC(DIMV(levelPC),nmat*(1+SDIM))
      REAL_T sol(DIMV(sol),nparts_def*SDIM)
      REAL_T recon(DIMV(recon),nmat*ngeom_recon)

      REAL_T xvel(DIMV(xvel))
      REAL_T xsrc(DIMV(xsrc))
      REAL_T ysrc(DIMV(ysrc))
      REAL_T zsrc(DIMV(zsrc))
 
      INTEGER_T im
      INTEGER_T im_left,im_right
      INTEGER_T im_prescribed_left,im_prescribed_right
      INTEGER_T im_prescribed_main
      INTEGER_T im_prescribed_both
      INTEGER_T im_fluid
      INTEGER_T im_adjoin
      INTEGER_T partid 
      INTEGER_T atwall 
      INTEGER_T idx
      INTEGER_T i,j,k
      INTEGER_T imac,jmac,kmac
      INTEGER_T ifluid,jfluid,kfluid
      INTEGER_T isolid,jsolid,ksolid
      INTEGER_T i_adjoin,j_adjoin,k_adjoin
      INTEGER_T im1,jm1,km1
      INTEGER_T ii,jj,kk
      INTEGER_T iii,jjj,kkk
      INTEGER_T nhalf
      REAL_T xstenMAC(-3:3,SDIM)
      REAL_T xsten(-3:3,SDIM)
      REAL_T LSleft(nmat)
      REAL_T LSright(nmat)
      REAL_T LSadjoin(nmat)
      INTEGER_T at_RZ_face
      INTEGER_T localbc
      INTEGER_T velcomp
      INTEGER_T nrmcomp
      INTEGER_T side
      INTEGER_T local_dir
      INTEGER_T fluid_stencil_ok
      INTEGER_T nvalid
      INTEGER_T maskleft,maskright
      REAL_T velsum
      REAL_T veltest
      REAL_T local_solid_vel(SDIM)
      REAL_T local_fluid_vel(SDIM)
      REAL_T u_ghost
      REAL_T un_fluid,un_solid
      REAL_T un_fluid_simple,un_solid_simple
      REAL_T nsolid(SDIM)
      REAL_T local_err
      REAL_T mag
      REAL_T RR
      REAL_T default_normal
      INTEGER_T normal_ok,charfn,charfn_left,charfn_right
      REAL_T simple_normal(SDIM)

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      nhalf=3
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid GFMUPDATE"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid GFMUPDATE"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((bfact_c.ne.bfact).and.(bfact_c.ne.2*bfact)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((bfact_f.ne.bfact).and.(bfact.ne.2*bfact_f)) then
       print *,"bfact_f invalid"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid GFMUPDATE"
       stop
      endif
      if ((dirmac.ge.0).and.(dirmac.le.SDIM-1)) then
       ! do nothing
      else
       print *,"dirmac invalid"
       stop
      endif
      if (gfm_error.ge.zero) then
       ! do nothing
      else
       print *,"gfm_error invalid"
       stop
      endif
      if ((project_option.eq.0).or. &
          (project_option.eq.1).or. &
          (project_option.eq.10).or. &
          (project_option.eq.13).or. & !FSI_material_exists 1st project
          (project_option.eq.11)) then !FSI_material_exists 2nd project
       ! do nothing
      else
       print *,"project_option invalid"
       stop
      endif

      if ((prescribed_solid_method.ge.0).and. &
          (prescribed_solid_method.le.3)) then 
       ! do nothing
      else
       print *,"prescribed_solid_method invalid"
       stop
      endif

      do im=1,nmat

       if ((prescribed_solid_scale(im).ge.zero).and. &
           (prescribed_solid_scale(im).lt.half)) then
        ! do nothing
       else
        print *,"prescribed_solid_scale invalid"
        stop
       endif

       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"denconst invalid"
        stop
       endif

      enddo ! im=1..nmat

      if ((slab_step.ge.0).and.(slab_step.lt.bfact_time_order)) then
       ! do nothing
      else
       print *,"slab_step invalid gfmupdate"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(xvel),0,dirmac,231)
      call checkbound(fablo,fabhi,DIMS(masknbr),1,-1,231)
      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,231)
      call checkbound(fablo,fabhi,DIMS(recon),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(xsrc),1,0,234)
      call checkbound(fablo,fabhi,DIMS(ysrc),1,1,234)
      call checkbound(fablo,fabhi,DIMS(zsrc),1,SDIM-1,234)
      call checkbound(fablo,fabhi,DIMS(sol),1,-1,234)
      call checkbound(fablo,fabhi,DIMS(levelPC),2,-1,234)

      do im=1,nmat

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).ge.1).and. &
                (fort_material_type(im).le.fort_max_num_eos)) then
        ! do nothing
       else
        print *,"fort_material_type invalid"
        stop
       endif

      enddo  ! im=1..nmat

      ii=0
      jj=0
      kk=0
      if (dirmac.eq.0) then
       ii=1
      else if (dirmac.eq.1) then
       jj=1
      else if ((dirmac.eq.2).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"dirmac out of range in GFMUPDATE"
       stop
      endif 

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dirmac)
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       im1=i-ii
       jm1=j-jj
       km1=k-kk

       maskleft=NINT(maskcov(D_DECL(im1,jm1,km1)))
       maskright=NINT(maskcov(D_DECL(i,j,k)))

       if ((maskleft.eq.1).and.(maskright.eq.1)) then

        call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dirmac+1)

        if (dirmac.eq.0) then
         idx=i
        else if (dirmac.eq.1) then
         idx=j
        else if ((dirmac.eq.2).and.(SDIM.eq.3)) then
         idx=k
        else
         print *,"dirmac invalid GFMUPDATE"
         stop
        endif

        ! levelPC() has piecewise constant BC at coarse/fine borders.
        do im=1,nmat
         LSleft(im)=levelPC(D_DECL(im1,jm1,km1),im)
         LSright(im)=levelPC(D_DECL(i,j,k),im)
        enddo

        at_RZ_face=0
        if (levelrz.eq.0) then
         ! do nothing
        else if (levelrz.eq.1) then
         if (SDIM.eq.2) then
          if (dirmac.eq.0) then
           if (xstenMAC(0,1).gt.VOFTOL*dx(1)) then
            ! do nothing
           else if (abs(xstenMAC(0,1)).le.VOFTOL*dx(1)) then 
            at_RZ_face=1
           else
            print *,"xstenMAC bust"
            stop
           endif
          else if (dirmac.eq.1) then
           ! do nothing
          else
           print *,"dirmac invalid"
           stop
          endif
         else
          print *,"dimension bust"
          stop
         endif
        else if (levelrz.eq.3) then
         if ((dirmac.eq.0).and. &
             (xstenMAC(0,1).le.VOFTOL*dx(1))) then
          at_RZ_face=1
         endif
        else
         print *,"levelrz invalid GFMUPDATE"
         stop
        endif 

        if (at_RZ_face.eq.1) then
         ! do nothing
        else if (at_RZ_face.eq.0) then

         call get_primary_material(LSleft,nmat,im_left)
         call get_primary_material(LSright,nmat,im_right)

         atwall=0
         do side=1,2
          localbc=velbc_in(dirmac+1,side,dirmac+1)
          if ( ((side.eq.1).and.(idx.eq.fablo(dirmac+1))).or. &
               ((side.eq.2).and.(idx.eq.fabhi(dirmac+1)+1)) ) then
           if ((localbc.eq.EXT_DIR).or. &
               (localbc.eq.REFLECT_ODD).or. &
               (localbc.eq.REFLECT_EVEN).or. &
               (localbc.eq.FOEXTRAP)) then
            atwall=1
           else if (localbc.eq.INT_DIR) then
            ! do nothing
           else
            print *,"localbc invalid"
            stop
           endif
          else if ((side.eq.1).and. &
                   (idx.gt.fablo(dirmac+1)).and. &
                   (idx.le.fabhi(dirmac+1)+1)) then
           ! do nothing
          else if ((side.eq.2).and. &
                   (idx.ge.fablo(dirmac+1)).and. &
                   (idx.lt.fabhi(dirmac+1)+1)) then
           ! do nothing
          else
           print *,"idx or side invalid"
           stop
          endif
         enddo ! side=1,2

         if (atwall.eq.1) then
          ! do nothing
         else if (atwall.eq.0) then

          im_prescribed_left=0
          im_prescribed_right=0

          if (is_rigid(nmat,im_left).eq.1) then
            if (is_prescribed(nmat,im_left).eq.1) then
             if ((prescribed_solid_scale(im_left).gt.zero).and. &
                 (prescribed_solid_scale(im_left).le.half)) then
              im_prescribed_left=im_left
             else if (prescribed_solid_scale(im_left).eq.zero) then
              ! do nothing
             else
              print *,"prescribed_solid_scale(im_left) invalid"
              stop
             endif
            else if (is_prescribed(nmat,im_left).eq.0) then
             ! do nothing
            else
             print *,"is_prescribed(nmat,im_left) invalid"
             stop
            endif
          else if (is_rigid(nmat,im_left).eq.0) then
            ! do nothing
          else
            print *,"is_rigid(nmat,im_left) invalid"
            stop
          endif

          if (is_rigid(nmat,im_right).eq.1) then
            if (is_prescribed(nmat,im_right).eq.1) then
             if ((prescribed_solid_scale(im_right).gt.zero).and. &
                 (prescribed_solid_scale(im_right).le.half)) then
              im_prescribed_right=im_right
             else if (prescribed_solid_scale(im_right).eq.zero) then
              ! do nothing
             else
              print *,"prescribed_solid_scale(im_right) invalid"
              stop
             endif
            else if (is_prescribed(nmat,im_right).eq.0) then
             ! do nothing
            else
             print *,"is_prescribed(nmat,im_right) invalid"
             stop
            endif
          else if (is_rigid(nmat,im_right).eq.0) then
            ! do nothing
          else
            print *,"is_rigid(nmat,im_right) invalid"
            stop
          endif

          im_prescribed_main=0
          im_fluid=0
          im_prescribed_both=0
          default_normal=zero

          if ((im_prescribed_left.eq.0).and. &
              (im_prescribed_right.eq.0)) then
           ! do nothing
          else if ((im_prescribed_left.eq.im_prescribed_right).and. &
                   (im_prescribed_left.ge.1).and. &
                   (im_prescribed_left.le.nmat)) then
           im_prescribed_both=im_prescribed_left
          else if ((im_prescribed_left.ge.1).and. &
                   (im_prescribed_left.le.nmat).and. &
                   (im_prescribed_right.eq.0)) then
           im_prescribed_main=im_prescribed_left
           default_normal=-one
           isolid=im1
           jsolid=jm1
           ksolid=km1
           im_fluid=im_right
           ifluid=i
           jfluid=j
           kfluid=k
          else if ((im_prescribed_right.ge.1).and. &
                   (im_prescribed_right.le.nmat).and. &
                   (im_prescribed_left.eq.0)) then
           im_prescribed_main=im_prescribed_right
           default_normal=one
           isolid=i
           jsolid=j
           ksolid=k
           im_fluid=im_left
           ifluid=im1
           jfluid=jm1
           kfluid=km1
          else
           print *,"im_prescribed invalid"
           stop
          endif

          if ((im_prescribed_main.ge.1).and. &
              (im_prescribed_main.le.nmat).and. &
              (im_fluid.ge.1).and. &
              (im_fluid.le.nmat)) then

           if (is_rigid(nmat,im_fluid).eq.0) then

            partid=0
            do im=1,nmat
             if (is_lag_part(nmat,im).eq.1) then
              if (im.le.im_prescribed_main) then
               partid=partid+1
              endif
             else if (is_lag_part(nmat,im).eq.0) then
              ! do nothing
             else
              print *,"is_lag_part invalid"
              stop
             endif
            enddo ! im=1..nmat
            if (im_prescribed_main.ne. &
                im_solid_map(partid)+1) then
             print *,"im_solid_map invalid"
             stop
            endif

            call gridsten_level(xsten,isolid,jsolid,ksolid,level,nhalf)

            normal_ok=1
            do local_dir=1,SDIM
             iii=0
             jjj=0
             kkk=0
             if (local_dir.eq.1) then
              iii=1
             else if (local_dir.eq.2) then
              jjj=1
             else if ((local_dir.eq.3).and.(SDIM.eq.3)) then
              kkk=1
             else
              print *,"local_dir invalid"
              stop
             endif
             
             do side=-1,1,2
              i_adjoin=ifluid+side*iii
              j_adjoin=jfluid+side*jjj
              k_adjoin=kfluid+side*kkk

              do im=1,nmat
               LSadjoin(im)=levelPC(D_DECL(i_adjoin,j_adjoin,k_adjoin),im)
              enddo
              call get_primary_material(LSadjoin,nmat,im_adjoin)

              if (im_adjoin.eq.im_prescribed_main) then
               charfn=1
              else if ((im_adjoin.ge.1).and. &
                       (im_adjoin.le.nmat).and. &
                       (im_adjoin.ne.im_prescribed_main)) then
               charfn=0
              else
               print *,"im_adjoin invalid"
               stop
              endif

              if (side.eq.-1) then
               charfn_left=charfn
              else if (side.eq.1) then
               charfn_right=charfn
              else
               print *,"side invalid"
               stop
              endif
             enddo ! side=-1,1,2
        
             if ((charfn_left.eq.1).and.(charfn_right.eq.1)) then
              normal_ok=0
              simple_normal(local_dir)=zero
             else if ((charfn_left.eq.1).and.(charfn_right.eq.0)) then
              simple_normal(local_dir)=-one
             else if ((charfn_left.eq.0).and.(charfn_right.eq.1)) then
              simple_normal(local_dir)=one
             else if ((charfn_left.eq.0).and.(charfn_right.eq.0)) then
              simple_normal(local_dir)=zero
             else
              print *,"charfn_left or charfn_right invalid"
              stop
             endif
            enddo ! local_dir=1..sdim

            if (simple_normal(dirmac+1).ne.default_normal) then
             normal_ok=0
            endif

            if (normal_ok.eq.0) then
             do local_dir=1,SDIM
              simple_normal(local_dir)=zero
             enddo
             simple_normal(dirmac+1)=default_normal
            else if (normal_ok.eq.1) then
             mag=zero
             do local_dir=1,SDIM
              mag=mag+simple_normal(local_dir)**2
             enddo 
             if (mag.ge.one) then
              mag=sqrt(mag)
              do local_dir=1,SDIM
               simple_normal(local_dir)=simple_normal(local_dir)/mag
              enddo
             else
              print *,"mag invalid"
              stop
             endif
            else
             print *,"normal_ok invalid"
             stop
            endif

            do local_dir=1,SDIM
             velcomp=(partid-1)*SDIM+local_dir
             nrmcomp=nmat+(im_prescribed_main-1)*SDIM+local_dir 
             if ((im_prescribed_main.ge.1).and. &
                 (im_prescribed_main.le.nmat)) then
              local_solid_vel(local_dir)= &
                 half*(sol(D_DECL(isolid,jsolid,ksolid),velcomp)+ &
                       sol(D_DECL(ifluid,jfluid,kfluid),velcomp))
              nsolid(local_dir)=levelPC(D_DECL(isolid,jsolid,ksolid),nrmcomp)
             else 
              print *,"im_prescribed_main invalid"
              stop
             endif
            enddo ! local_dir=1..sdim

            fluid_stencil_ok=1

            RR=one
            if (levelrz.eq.0) then
             ! do nothing
            else if (levelrz.eq.1) then
             RR=xsten(0,1)
             if (RR.gt.zero) then
              ! do nothing
             else
              print *,"RR invalid"
              stop
             endif
            else if (levelrz.eq.3) then
             ! if R-Theta, then N(2) -> N(2)/RR + renormalize.
             RR=xsten(0,1)
             if (RR.gt.zero) then
              ! do nothing
             else
              print *,"RR invalid"
              stop
             endif
            else
             print *,"levelrz invalid"
             stop
            endif

            call prepare_normal(nsolid,RR,mag)
            if (mag.eq.zero) then
             fluid_stencil_ok=0
            else if (mag.gt.zero) then
             ! do nothing
            else
             print *,"mag invalid"
             stop
            endif

            if (fluid_stencil_ok.eq.1) then

             un_solid=zero
             un_solid_simple=zero
             do local_dir=1,SDIM
              un_solid=un_solid+local_solid_vel(local_dir)*nsolid(local_dir)
              un_solid_simple=un_solid_simple+ &
                      local_solid_vel(local_dir)*simple_normal(local_dir)
             enddo

             do local_dir=1,SDIM

              nvalid=0
              velsum=zero

              iii=0
              jjj=0
              kkk=0
              if (local_dir.eq.1) then
               iii=1
              else if (local_dir.eq.2) then
               jjj=1
              else if ((local_dir.eq.3).and.(SDIM.eq.3)) then
               kkk=1
              else
               print *,"local_dir invalid"
               stop
              endif
             
              do side=-1,1,2
               i_adjoin=ifluid+side*iii
               j_adjoin=jfluid+side*jjj
               k_adjoin=kfluid+side*kkk

               do im=1,nmat
                LSadjoin(im)=levelPC(D_DECL(i_adjoin,j_adjoin,k_adjoin),im)
               enddo
               call get_primary_material(LSadjoin,nmat,im_adjoin)

               if (is_rigid(nmat,im_adjoin).eq.0) then
                imac=ifluid+(side+1)*iii/2
                jmac=jfluid+(side+1)*jjj/2
                kmac=kfluid+(side+1)*kkk/2

                if (local_dir.eq.1) then
                 veltest=xsrc(D_DECL(imac,jmac,kmac))
                else if (local_dir.eq.2) then
                 veltest=ysrc(D_DECL(imac,jmac,kmac))
                else if ((local_dir.eq.3).and.(SDIM.eq.3)) then
                 veltest=zsrc(D_DECL(imac,jmac,kmac))
                else
                 print *,"local_dir invalid"
                 stop
                endif

                nvalid=nvalid+1
                velsum=velsum+veltest
               else if (is_rigid(nmat,im_adjoin).eq.1) then
                ! do nothing
               else
                print *,"is_rigid(nmat,im_adjoin) invalid"
                stop
               endif
              enddo ! side=-1,1,2

              if (nvalid.eq.0) then
               fluid_stencil_ok=0
              else if ((nvalid.eq.1).or.(nvalid.eq.2)) then
               local_fluid_vel(local_dir)=velsum/nvalid
              else
               print *,"nvalid invalid"
               stop
              endif
             enddo ! local_dir=1..sdim

             if (fluid_stencil_ok.eq.1) then

               un_fluid=zero
               un_fluid_simple=zero
               do local_dir=1,SDIM
                un_fluid=un_fluid+local_fluid_vel(local_dir)*nsolid(local_dir)
                un_fluid_simple=un_fluid_simple+ &
                        local_fluid_vel(local_dir)*simple_normal(local_dir)
               enddo

               if (normal_ok.eq.1) then
                if (prescribed_solid_method.eq.0) then ! stair-step
                 u_ghost=local_solid_vel(dirmac+1)
                else if (prescribed_solid_method.eq.1) then ! 2nd order
                 u_ghost=local_fluid_vel(dirmac+1)+ &
                  (un_solid-un_fluid)*nsolid(dirmac+1)
                else if (prescribed_solid_method.eq.2) then ! mod stair-step
                 u_ghost=local_solid_vel(dirmac+1)
                else if (prescribed_solid_method.eq.3) then ! mod normal
                 u_ghost=local_fluid_vel(dirmac+1)+ &
                  (un_solid_simple-un_fluid_simple)*simple_normal(dirmac+1)
                else
                 print *,"prescribed_solid_method invalid"
                 stop
                endif
               else if (normal_ok.eq.0) then
                u_ghost=local_solid_vel(dirmac+1)
               else
                print *,"normal_ok invalid"
                stop
               endif

               local_err=abs(u_ghost-xvel(D_DECL(i,j,k)))
               if (local_err.gt.gfm_error) then
                gfm_error=local_err
               endif
               xvel(D_DECL(i,j,k))=u_ghost

             else if (fluid_stencil_ok.eq.0) then
              u_ghost=local_solid_vel(dirmac+1)
              xvel(D_DECL(i,j,k))=u_ghost
             else
              print *,"fluid_stencil_ok invalid"
              stop
             endif

            else if (fluid_stencil_ok.eq.0) then
             u_ghost=local_solid_vel(dirmac+1)
             xvel(D_DECL(i,j,k))=u_ghost
            else
             print *,"fluid_stencil_ok invalid"
             stop
            endif

           else if (is_rigid(nmat,im_fluid).eq.1) then
            ! do nothing
           else
            print *,"is_rigid(nmat,im_fluid) invalid"
            stop
           endif

          else if (im_prescribed_main.eq.0) then

           if (im_fluid.eq.0) then

            if ((im_prescribed_both.ge.1).and. &
                (im_prescribed_both.le.nmat)) then

             partid=0
             do im=1,nmat
              if (is_lag_part(nmat,im).eq.1) then
               if (im.le.im_prescribed_both) then
                partid=partid+1
               endif
              else if (is_lag_part(nmat,im).eq.0) then
               ! do nothing
              else
               print *,"is_lag_part invalid"
               stop
              endif
             enddo ! im=1..nmat
             if (im_prescribed_both.ne. &
                 im_solid_map(partid)+1) then
              print *,"im_solid_map invalid"
              stop
             endif

             do local_dir=1,SDIM
              velcomp=(partid-1)*SDIM+local_dir
              if ((im_prescribed_both.ge.1).and. &
                  (im_prescribed_both.le.nmat)) then
               local_solid_vel(local_dir)= &
                 half*(sol(D_DECL(i,j,k),velcomp)+ &
                       sol(D_DECL(im1,jm1,km1),velcomp))
              else 
               print *,"im_prescribed_both invalid"
               stop
              endif
             enddo ! local_dir=1..sdim

             u_ghost=local_solid_vel(dirmac+1)
             xvel(D_DECL(i,j,k))=u_ghost

            else if (im_prescribed_both.eq.0) then
             ! do nothing
            else
             print *,"im_prescribed_both invalid"
             stop
            endif

           else
            print *,"im_fluid invalid in GFMUPDATE, im_fluid=",im_fluid
            stop
           endif
          else
           print *,"im_prescribed_main or im_fluid invalid"
           stop
          endif

         else
          print *,"atwall invalid"
          stop
         endif

        else
         print *,"at_RZ_face invalid"
         stop
        endif

       else if ((maskleft.eq.0).or.(maskright.eq.0)) then
        ! do nothing
       else
        print *,"maskleft or maskright invalid"
        stop
       endif

      enddo
      enddo
      enddo

      return
      end subroutine FORT_GFMUPDATE


      ! called from: NavierStokes::allocate_FACE_WEIGHT (NavierStokes3.cpp)
      !  which is called from:
      !   NavierStokes::update_SEM_forcesALL
      !   NavierStokes::multiphase_project
      !   NavierStokes::diffusion_heatingALL 
      ! mask=1 at fine-fine boundaries
      subroutine FORT_BUILDFACEWT( &
       facewt_iter, &
       num_materials_face, &
       level, &
       finest_level, &
       nsolve, &
       nsolveMM, &
       nsolveMM_FACE, &
       nfacefrac, &
       ncellfrac, &
       local_face_index, &
       facecut_index, &
       icefacecut_index, &
       ncphys, &
       nmat, &
       xlo,dx, &
       offdiagcheck, &
       DIMS(offdiagcheck), &
       recon,DIMS(recon), &
       cenden,DIMS(cenden), &
       cenvisc,DIMS(cenvisc), &
       cellmm,DIMS(cellmm), &
       xfacemm,DIMS(xfacemm), &
       yfacemm,DIMS(yfacemm), &
       zfacemm,DIMS(zfacemm), &
       xfwt,DIMS(xfwt), &
       yfwt,DIMS(yfwt), &
       zfwt,DIMS(zfwt), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       alt_xface, &
       DIMS(alt_xface), &
       alt_yface, &
       DIMS(alt_yface), &
       alt_zface, &
       DIMS(alt_zface), &
       mask,DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       min_face_wt, & ! static variable
       max_face_wt, & ! static variable
       singular_possible, &
       solvability_projection, &
       solvability_level_flag, &
       presbc_arr, &
       visc_coef, &
       constant_viscosity, &
       prescribed_solid_scale, &
       project_option)
      use global_utility_module
      use probcommon_module
      use probf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: facewt_iter
      INTEGER_T, intent(in) :: num_materials_face
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: nsolveMM
      INTEGER_T, intent(in) :: nsolveMM_FACE
      INTEGER_T :: nsolveMM_FACE_test
      INTEGER_T, intent(in) :: nfacefrac
      INTEGER_T, intent(in) :: ncellfrac
      REAL_T, intent(in) :: prescribed_solid_scale(nmat)
      INTEGER_T, intent(in) :: local_face_index
      INTEGER_T, intent(in) :: facecut_index
      INTEGER_T, intent(in) :: icefacecut_index
      INTEGER_T, intent(in) :: ncphys
      REAL_T, intent(in) :: visc_coef
      INTEGER_T, intent(in) :: constant_viscosity
      INTEGER_T, intent(in) :: project_option
      REAL_T, intent(inout) :: min_face_wt(4)
      REAL_T, intent(inout) :: max_face_wt(4)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: DIMDEC(offdiagcheck)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(cenden)
      INTEGER_T, intent(in) :: DIMDEC(cenvisc)
      INTEGER_T, intent(in) :: DIMDEC(cellmm)
      INTEGER_T, intent(in) :: DIMDEC(xfacemm)
      INTEGER_T, intent(in) :: DIMDEC(yfacemm)
      INTEGER_T, intent(in) :: DIMDEC(zfacemm)
      INTEGER_T, intent(in) :: DIMDEC(xfwt)
      INTEGER_T, intent(in) :: DIMDEC(yfwt)
      INTEGER_T, intent(in) :: DIMDEC(zfwt)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(alt_xface)
      INTEGER_T, intent(in) :: DIMDEC(alt_yface)
      INTEGER_T, intent(in) :: DIMDEC(alt_zface)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      REAL_T, intent(inout) :: offdiagcheck(DIMV(offdiagcheck),nsolveMM) 
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon) 
      REAL_T, intent(in) :: cenden(DIMV(cenden),nmat+1) 
      REAL_T, intent(in) :: cenvisc(DIMV(cenvisc),nmat+1) 
      REAL_T, intent(in) :: cellmm(DIMV(cellmm),ncellfrac) 
      REAL_T, intent(in) :: xfacemm(DIMV(xfacemm),nfacefrac) 
      REAL_T, intent(in) :: yfacemm(DIMV(yfacemm),nfacefrac) 
      REAL_T, intent(in) :: zfacemm(DIMV(zfacemm),nfacefrac) 
      REAL_T, intent(out) :: xfwt(DIMV(xfwt),nsolveMM_FACE)
      REAL_T, intent(out) :: yfwt(DIMV(yfwt),nsolveMM_FACE)
      REAL_T, intent(out) :: zfwt(DIMV(zfwt),nsolveMM_FACE)
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: yface(DIMV(yface),ncphys)
      REAL_T, intent(in) :: zface(DIMV(zface),ncphys)
      REAL_T, intent(in) :: alt_xface(DIMV(alt_xface))
      REAL_T, intent(in) :: alt_yface(DIMV(alt_yface))
      REAL_T, intent(in) :: alt_zface(DIMV(alt_zface))
      REAL_T, intent(in) :: mask(DIMV(mask))
      INTEGER_T, intent(in) :: singular_possible
      INTEGER_T, intent(in) :: solvability_projection
      INTEGER_T, intent(inout) :: solvability_level_flag
      INTEGER_T, intent(in) :: presbc_arr(SDIM,2,nsolveMM)
  
      INTEGER_T i,j,k
      INTEGER_T iface,jface,kface
      INTEGER_T ii,jj,kk
      INTEGER_T inorm
      REAL_T dd,dd_group
      REAL_T cc,cc_group
      REAL_T cc_ice
      REAL_T alt_cc_ice
      REAL_T local_cc_ice
      INTEGER_T side
      INTEGER_T dir
      INTEGER_T veldir
      INTEGER_T velcompL
      INTEGER_T velcompR
      INTEGER_T velcomp
      REAL_T local_wt(nsolve)
      REAL_T local_fwt
      INTEGER_T local_presbc
      INTEGER_T im_vel
      INTEGER_T im
      REAL_T local_mask
      INTEGER_T caller_id

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level or finest_level invalid build face wt"
       stop
      endif

       ! indexes start at 0
      if ((facecut_index.ne.3).or. &
          (icefacecut_index.ne.4).or. &
          (local_face_index+1.gt.ncphys)) then
       print *,"face_index bust 6"
       stop
      endif
      if (ncphys.lt.8) then
       print *,"ncphys invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid101"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      do im=1,nmat
       if ((prescribed_solid_scale(im).ge.zero).and. &
           (prescribed_solid_scale(im).lt.one)) then
        ! do nothing
       else
        print *,"prescribed_solid_scale(im) invalid"
        stop
       endif
      enddo ! im=1..nmat
      if ((singular_possible.ne.0).and. &
          (singular_possible.ne.1)) then
       print *,"singular_possible invalid"
       stop
      endif
      if ((solvability_projection.ne.0).and. &
          (solvability_projection.ne.1)) then
       print *,"solvability projection invalid"
       stop
      endif
      if ((solvability_level_flag.ne.0).and. &
          (solvability_level_flag.ne.1)) then
       print *,"solvability level flag invalid"
       stop
      endif
      if ((constant_viscosity.ne.0).and. &
          (constant_viscosity.ne.1)) then
       print *,"constant_viscosity invalid"
       stop
      endif
      if ((num_materials_face.ne.1).and. &
          (num_materials_face.ne.nmat)) then
       print *,"num_materials_face invalid"
       stop
      endif

      if (visc_coef.ge.zero) then
       ! do nothing
      else
       print *,"visc_coef invalid"
       stop
      endif
      if ((nsolve.ne.1).and.(nsolve.ne.AMREX_SPACEDIM)) then
       print *,"nsolve invalid8"
       stop
      endif 
      if (nsolveMM.ne.nsolve*num_materials_face) then
       print *,"nsolveMM invalid 13441"
       stop
      endif
      nsolveMM_FACE_test=nsolveMM
      if (num_materials_face.eq.1) then
       ! do nothing
      else if (num_materials_face.eq.nmat) then
       nsolveMM_FACE_test=2*nsolveMM_FACE_test
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

      call checkbound(fablo,fabhi,DIMS(offdiagcheck),0,-1,241)
      call checkbound(fablo,fabhi,DIMS(recon),1,-1,241)
      call checkbound(fablo,fabhi,DIMS(cenden),1,-1,241)
      call checkbound(fablo,fabhi,DIMS(cenvisc),1,-1,241)
      call checkbound(fablo,fabhi,DIMS(cellmm),0,-1,234)
      call checkbound(fablo,fabhi,DIMS(xfacemm),0,0,264)
      call checkbound(fablo,fabhi,DIMS(yfacemm),0,1,264)
      call checkbound(fablo,fabhi,DIMS(zfacemm),0,SDIM-1,264)
      call checkbound(fablo,fabhi,DIMS(xfwt),0,0,242)
      call checkbound(fablo,fabhi,DIMS(yfwt),0,1,242)
      call checkbound(fablo,fabhi,DIMS(zfwt),0,SDIM-1,242)
      call checkbound(fablo,fabhi,DIMS(xface),0,0,244)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,244)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,244)
      call checkbound(fablo,fabhi,DIMS(alt_xface),0,0,244)
      call checkbound(fablo,fabhi,DIMS(alt_yface),0,1,244)
      call checkbound(fablo,fabhi,DIMS(alt_zface),0,SDIM-1,244)
      call checkbound(fablo,fabhi,DIMS(mask),1,-1,246)

      if (facewt_iter.eq.0) then

       do dir=0,SDIM-1

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
         print *,"dir out of range in BUILDFACEWT"
         stop
        endif 
  
        call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)

        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)
  
          ! projection: dedge is 1/rho  (faceden_index) 
          ! viscosity: dedge is facevisc_index
          ! temperature: dedge is faceheat_index
          ! species: dedge is facespecies_index

          if (dir.eq.0) then
           inorm=i
          else if (dir.eq.1) then
           inorm=j
          else if ((dir.eq.2).and.(SDIM.eq.3)) then
           inorm=k
          else
           print *,"dir invalid buildfacewt"
           stop
          endif

          side=0
          local_presbc=INT_DIR
          local_mask=zero

          if (inorm.eq.fablo(dir+1)) then
           side=1
           local_mask=mask(D_DECL(i-ii,j-jj,k-kk))
          else if (inorm.eq.fabhi(dir+1)+1) then
           side=2
           local_mask=mask(D_DECL(i,j,k))
          else if ((inorm.gt.fablo(dir+1)).and. &
                   (inorm.lt.fabhi(dir+1)+1)) then
           ! do nothing
          else
           print *,"inorm invalid"
           print *,"dir,inorm,fablo,fabhi ", &
               dir,inorm,fablo(dir+1),fabhi(dir+1)
           stop
          endif

          if ((project_option.eq.0).or. &  !regular project
              (project_option.eq.1).or. &  !initial project
              (project_option.eq.10).or. & !sync project
              (project_option.eq.13).or. & !FSI_material_exists 1st project
              (project_option.eq.11)) then !FSI_material_exists 2nd project

           if (singular_possible.eq.1) then

            if (solvability_projection.eq.0) then
             solvability_level_flag=0
            else if (solvability_projection.eq.1) then
             if (side.eq.0) then
              ! do nothing
             else if ((side.eq.1).or.(side.eq.2)) then
              local_presbc=presbc_arr(dir+1,side,1)
              if (local_presbc.eq.INT_DIR) then
               if (local_mask.eq.zero) then 
                solvability_level_flag=0  ! coarse/fine BC
               else if (local_mask.eq.one) then
                ! do nothing (periodic BC)
               else
                print *,"local_mask invalid"
                stop
               endif
              else if (local_presbc.eq.EXT_DIR) then
               print *,"cannot have outflow bc with solvability constraint"
               stop
              else if ((local_presbc.eq.REFLECT_EVEN).or. &
                       (local_presbc.eq.FOEXTRAP)) then
               ! do nothing
              else
               print *,"local_presbc invalid"
               stop
              endif
             else
              print *,"side invalid"
              stop
             endif
            else
             print *,"solvability projection invalid"
             stop
            endif

           else
            print *,"singular_possible invalid"
            stop
           endif

          else if (project_option.eq.12) then ! pressure extrapolation
           if (singular_possible.eq.1) then
            solvability_level_flag=0
           else
            print *,"singular_possible invalid"
            stop
           endif
          else if (project_option.eq.2) then ! temperature
           if (singular_possible.eq.0) then
            solvability_level_flag=0
           else
            print *,"singular_possible invalid"
            stop
           endif
          else if ((project_option.ge.100).and. &
                   (project_option.lt.100+num_species_var)) then 
           if (singular_possible.eq.0) then
            solvability_level_flag=0
           else
            print *,"singular_possible invalid"
            stop
           endif
          else if (project_option.eq.3) then ! viscosity
           if (singular_possible.eq.0) then
            solvability_level_flag=0
           else
            print *,"singular_possible invalid"
            stop
           endif
          else
           print *,"project_option invalid"
           stop
          endif

          do im_vel=1,num_materials_face
           do veldir=1,nsolve

            if (nsolve.eq.1) then
             velcomp=im_vel
            else if (nsolve.eq.SDIM) then
             velcomp=veldir
            else
             print *,"nsolve invalid"
             stop
            endif

            if (side.eq.0) then
             ! do nothing
            else if ((side.eq.1).or.(side.eq.2)) then
             local_presbc=presbc_arr(dir+1,side,velcomp)
            else
             print *,"side invalid, side= ",side
             stop
            endif

            if (dir.eq.0) then
             cc=xface(D_DECL(i,j,k),facecut_index+1)
             cc_ice=xface(D_DECL(i,j,k),icefacecut_index+1)
             alt_cc_ice=alt_xface(D_DECL(i,j,k))
            else if (dir.eq.1) then
             cc=yface(D_DECL(i,j,k),facecut_index+1)
             cc_ice=yface(D_DECL(i,j,k),icefacecut_index+1)
             alt_cc_ice=alt_yface(D_DECL(i,j,k))
            else if ((dir.eq.2).and.(SDIM.eq.3)) then
             cc=zface(D_DECL(i,j,k),facecut_index+1)
             cc_ice=zface(D_DECL(i,j,k),icefacecut_index+1)
             alt_cc_ice=alt_zface(D_DECL(i,j,k))
            else
             print *,"dir invalid buildfacewt"
             stop
            endif
          
            if (dir.eq.0) then
             dd=xface(D_DECL(i,j,k),local_face_index+1)
            else if (dir.eq.1) then
             dd=yface(D_DECL(i,j,k),local_face_index+1)
            else if ((dir.eq.2).and.(SDIM.eq.3)) then
             dd=zface(D_DECL(i,j,k),local_face_index+1)
            else
             print *,"dir invalid buildfacewt"
             stop
            endif

            if (project_option.eq.13) then !FSI_material_exists 1st project
             local_cc_ice=alt_cc_ice
            else
             local_cc_ice=cc_ice
            endif

             ! in: PROB.F90
            caller_id=1
            call eval_face_coeff( &
             caller_id, &
             level,finest_level, &
             cc,local_cc_ice,cc_group, &
             dd,dd_group, &
             visc_coef, &
             nsolve,nsolveMM,im_vel,dir,veldir,project_option, &
             constant_viscosity,side,local_presbc,local_wt)

            if (dd_group.lt.min_face_wt(1)) then
             min_face_wt(1)=dd_group
            endif
            if (cc_group.lt.min_face_wt(2)) then
             min_face_wt(2)=cc_group
            endif
            if (dd_group.gt.max_face_wt(1)) then
             max_face_wt(1)=dd_group
            endif
            if (cc_group.gt.max_face_wt(2)) then
             max_face_wt(2)=cc_group
            endif

            if ((dd_group.ge.zero).and.(cc_group.ge.zero)) then
             ! do nothing
            else
             print *,"cannot have negative wts"
             stop
            endif

            if (local_wt(veldir).lt.min_face_wt(4)) then
             min_face_wt(4)=local_wt(veldir)
            endif
            if (local_wt(veldir).gt.max_face_wt(4)) then
             max_face_wt(4)=local_wt(veldir)
            endif

            if (num_materials_face.eq.1) then
             velcompL=veldir
             velcompR=veldir
            else if (num_materials_face.eq.nmat) then
             if (nsolve.eq.1) then
              velcompL=im_vel
              velcompR=velcompL+nsolveMM_FACE/2 
             else if (nsolve.eq.SDIM) then
              print *,"nsolve invalid nsolve=",nsolve
              stop
             else
              print *,"nsolve invalid nsolve=",nsolve
              stop
             endif
            else
             print *,"num_materials_face invalid"
             stop
            endif

            if (dir.eq.0) then
             xfwt(D_DECL(i,j,k),velcompL)=local_wt(veldir)
             xfwt(D_DECL(i,j,k),velcompR)=local_wt(veldir)
            else if (dir.eq.1) then
             yfwt(D_DECL(i,j,k),velcompL)=local_wt(veldir)
             yfwt(D_DECL(i,j,k),velcompR)=local_wt(veldir)
            else if ((dir.eq.2).and.(SDIM.eq.3)) then
             zfwt(D_DECL(i,j,k),velcompL)=local_wt(veldir)
             zfwt(D_DECL(i,j,k),velcompR)=local_wt(veldir)
            else
             print *,"dir invalid buildfacewt 2"
             stop
            endif
           enddo ! veldir=1..nsolve

          enddo ! im_vel=1..num_materials_face

        enddo ! k
        enddo ! j
        enddo ! i
       enddo ! dir=0,..,sdim-1

      else if (facewt_iter.eq.1) then

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        do im_vel=1,num_materials_face
         do veldir=1,nsolve

          if (num_materials_face.eq.1) then
           velcompL=veldir
           velcompR=veldir
          else if (num_materials_face.eq.nmat) then
           if (nsolve.eq.1) then
            velcompL=im_vel
            velcompR=velcompL+nsolveMM_FACE/2 
           else if (nsolve.eq.SDIM) then
            print *,"nsolve invalid"
            stop
           else
            print *,"nsolve invalid"
            stop
           endif
          else
           print *,"num_materials_face invalid"
           stop
          endif

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
            print *,"dir invalid"
            stop
           endif
           do side=1,2
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
            if (side.eq.1) then
             velcomp=velcompR
            else if (side.eq.2) then
             velcomp=velcompL
            else
             print *,"side invalid"
             stop
            endif
            if (dir.eq.1) then
             local_fwt=xfwt(D_DECL(iface,jface,kface),velcomp)
            else if (dir.eq.2) then
             local_fwt=yfwt(D_DECL(iface,jface,kface),velcomp)
            else if ((dir.eq.3).and.(SDIM.eq.3)) then
             local_fwt=zfwt(D_DECL(iface,jface,kface),velcomp)
            else
             print *,"dir invalid"
             stop
            endif
            if (nsolve.eq.1) then
             velcomp=im_vel
            else if (nsolve.eq.SDIM) then
             velcomp=veldir
            else
             print *,"nsolve invalid"
             stop
            endif
            offdiagcheck(D_DECL(i,j,k),velcomp)= &
               offdiagcheck(D_DECL(i,j,k),velcomp)+local_fwt
           enddo ! side=1..2
          enddo ! dir=1..sdim
         enddo ! veldir=1,nsolve
        enddo ! im_vel=1..num_materials_face

       enddo
       enddo
       enddo
      else
       print *,"facewt_iter invalid"
       stop
      endif

      return
      end subroutine FORT_BUILDFACEWT



       ! solid: velx,vely,velz,dist  (dist<0 in solid)
       ! called from: NavierStokes::prescribe_solid_geometry
       !   (declared in NavierStokes2.cpp)
      subroutine FORT_RENORMALIZE_PRESCRIBE( &
       tid, &
       level,finest_level, &
       solid_time, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       vofnew,DIMS(vofnew), &
       solid,DIMS(solid), &
       maskcov,DIMS(maskcov), &
       LS,DIMS(LS), & ! getStateDist(time)
       mofdata,DIMS(mofdata), &
       den,DIMS(den), &
       vel,DIMS(vel), &
       velnew,DIMS(velnew), &
       dennew,DIMS(dennew), &
       lsnew,DIMS(lsnew), &
       xlo,dx, &
       time, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       renormalize_only, &
       solidheat_flag)
      use global_utility_module
      use global_distance_module
      use probf90_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: solidheat_flag

      INTEGER_T, intent(in) :: renormalize_only
      INTEGER_T, intent(in) :: level,finest_level
      REAL_T, intent(in) :: solid_time

      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: nparts_def
      INTEGER_T, intent(in) :: im_solid_map(nparts_def)
      INTEGER_T, intent(in) :: bfact

      INTEGER_T, intent(in) :: DIMDEC(vofnew)
      INTEGER_T, intent(in) :: DIMDEC(solid)
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(LS)
      INTEGER_T, intent(in) :: DIMDEC(mofdata)
      INTEGER_T, intent(in) :: DIMDEC(den)
      INTEGER_T, intent(in) :: DIMDEC(vel)
      INTEGER_T, intent(in) :: DIMDEC(velnew)
      INTEGER_T, intent(in) :: DIMDEC(dennew)
      INTEGER_T, intent(in) :: DIMDEC(lsnew)
      REAL_T, intent(inout) ::  vofnew(DIMV(vofnew),nmat*ngeom_raw)
      REAL_T, intent(in) ::  vel(DIMV(vel), &
       num_materials_vel*(SDIM+1))
      REAL_T, intent(out) ::  velnew(DIMV(velnew), &
       num_materials_vel*(SDIM+1))
      REAL_T, intent(in) ::  LS(DIMV(LS),nmat*(1+SDIM))
      REAL_T, intent(in) ::  mofdata(DIMV(mofdata),nmat*ngeom_raw)
      REAL_T, intent(in) ::  den(DIMV(den),nmat*num_state_material)
      REAL_T, intent(inout) ::  dennew(DIMV(dennew),nmat*num_state_material)
      REAL_T, intent(inout) ::  lsnew(DIMV(lsnew),nmat*(1+SDIM))
      REAL_T, intent(in) ::  solid(DIMV(solid),SDIM*nparts_def)
      REAL_T, intent(in) ::  maskcov(DIMV(maskcov))
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)

      INTEGER_T i,j,k,dir
      INTEGER_T im,im_opp
      INTEGER_T im_primary_stencil
      INTEGER_T im_primary_sub_stencil
      INTEGER_T im_solid_max
      INTEGER_T vofcomp,vofcompraw
      INTEGER_T i1,j1,k1
      INTEGER_T i2,j2,k2
      REAL_T centroid(SDIM)
      REAL_T volcell
      REAL_T cencell(SDIM)
      REAL_T censolid_new(nmat,SDIM)

      REAL_T xsten(-7:7,SDIM)
      INTEGER_T nhalf
      REAL_T mofnew(nmat*ngeom_recon)
      INTEGER_T istenlo(3),istenhi(3)
      INTEGER_T LSstenlo(3),LSstenhi(3)
      REAL_T LS_solid_new(nmat)
      INTEGER_T local_maskcov
      REAL_T vfrac_solid_new(nmat)
      REAL_T F_stencil
      REAL_T F_stencil_sum
      INTEGER_T statecomp
      INTEGER_T statecomp_solid
      INTEGER_T istate,ispecies
      INTEGER_T dencomp,tempcomp,speccomp
      REAL_T den_hold(nmat*num_state_material)
      REAL_T state_stencil(num_state_material)
      REAL_T state_stencil_sum(num_state_material)
      REAL_T nslope_solid(SDIM)
      INTEGER_T nmax
      REAL_T LS_extend(D_DECL(-1:1,-1:1,-1:1),nmat)
      REAL_T LS_temp(D_DECL(-1:1,-1:1,-1:1))
      REAL_T LSfacearea
      REAL_T LScentroid(SDIM)
      REAL_T LSareacentroid(SDIM)
      INTEGER_T nrefine_geom
      REAL_T dxmaxLS
      REAL_T ls_hold(nmat*(1+SDIM))
      REAL_T max_solid_LS
      REAL_T sum_vfrac_solid_new
      REAL_T LS_predict(nmat)
      REAL_T LS_virtual(nmat)
      REAL_T LS_virtual_new(nmat)
      REAL_T LS_virtual_max
      REAL_T near_fluid
      REAL_T WT
      INTEGER_T nmat_fluid,nmat_solid,nmat_lag
      INTEGER_T at_center
      INTEGER_T ibase
      INTEGER_T partid
      INTEGER_T partid_max
      REAL_T xfluid,xghost
      INTEGER_T tessellate

      tessellate=0

      nmax=POLYGON_LIST_MAX  ! in: RENORMALIZE_PRESCRIBE
      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      nhalf=7

      if (bfact.lt.1) then
       print *,"bfact invalid102"
       stop
      endif 

      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if (solidheat_flag.eq.0) then 
       !do nothing (heat conduction in solid)
      else if (solidheat_flag.eq.1) then
       !do nothing (dirichlet at solid interface)
      else if (solidheat_flag.eq.2) then
       !do nothing (neumann at solid interface)
      else
       print *,"solidheat_flag invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat bust"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if ((renormalize_only.ne.0).and.(renormalize_only.ne.1)) then
       print *,"renormalize_only invalid"
       stop
      endif
      if (level.lt.0) then
       print *,"level invalid renormalize prescribe 1"
       stop
      endif
      if (level.gt.finest_level) then
       print *,"level invalid renormalize prescribe 2"
       stop
      endif

      if ((time.ge.zero).and.(time.le.1.0D+20)) then
       ! do nothing
      else if (time.ge.1.0D+20) then
       print *,"WARNING time.ge.1.0D+20 in renormalize"
      else if (time.lt.zero) then
       print *,"time invalid in renormalize"
       stop
      else
       print *,"time bust in renormalize"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid RENORMALIZE_PRESCRIBE"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid RENORMALIZE_PRESCRIBE"
       stop
      endif

      nmat_fluid=0
      nmat_solid=0
      nmat_lag=0

      do im=1,nmat

       if (num_state_material.ne. &
           num_state_base+num_species_var) then
        print *,"num_state_material invalid"
        stop
       endif

       if (is_lag_part(nmat,im).eq.1) then
        nmat_lag=nmat_lag+1
        if (is_rigid(nmat,im).eq.1) then
         nmat_solid=nmat_solid+1
        else if (is_rigid(nmat,im).eq.0) then
         nmat_fluid=nmat_fluid+1
        else
         print *,"is_rigid(nmat,im) invalid"
         stop
        endif
       else if (is_lag_part(nmat,im).eq.0) then
        if (is_rigid(nmat,im).eq.0) then
         nmat_fluid=nmat_fluid+1
        else
         print *,"is_rigid(nmat,im) invalid"
         stop
        endif
       else
        print *,"is_lag_part(nmat,im) invalid"
        stop
       endif

      enddo ! im=1..nmat

      if (nmat_lag.ne.nparts) then
       print *,"nmat_lag invalid"
       stop
      endif

      if (nmat_fluid+nmat_solid.ne.nmat) then
       print *,"nmat_fluid and/or nmat_solid invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(vofnew),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(velnew),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(solid),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(maskcov),0,-1,26)
      call checkbound(fablo,fabhi,DIMS(LS),3,-1,26)
      call checkbound(fablo,fabhi,DIMS(mofdata),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(vel),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(den),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(dennew),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(lsnew),1,-1,26)

      istenlo(3)=0
      istenhi(3)=0
      do dir=1,SDIM
       istenlo(dir)=-1
       istenhi(dir)=1
      enddo

      LSstenlo(3)=0
      LSstenhi(3)=0
      do dir=1,SDIM
       LSstenlo(dir)=-2
       LSstenhi(dir)=2
      enddo

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_maskcov=NINT(maskcov(D_DECL(i,j,k)))

       call gridsten_level(xsten,i,j,k,level,nhalf)

       call Box_volumeFAST(bfact,dx,xsten,nhalf,volcell,cencell,SDIM)

       if (local_maskcov.eq.1) then

        ! --------------------------------------------------------- 
        ! first: fluid state variable extrapolation into empty cells.
        ! ----------------------------------------------------------
        do im=1,nmat

         vofcompraw=(im-1)*ngeom_raw+1
         F_stencil=mofdata(D_DECL(i,j,k),vofcompraw)

         if (is_rigid(nmat,im).eq.1) then
          ! do nothing
         else if (is_rigid(nmat,im).eq.0) then

          if (F_stencil.gt.VOFTOL) then
           ! do nothing
          else if (abs(F_stencil).le.VOFTOL) then
           ! extrapolate into the empty cell.

           F_stencil_sum=zero
           do istate=1,num_state_material
            state_stencil_sum(istate)=zero
           enddo
  
           do i1=istenlo(1),istenhi(1)
           do j1=istenlo(2),istenhi(2)
           do k1=istenlo(3),istenhi(3)

            F_stencil=mofdata(D_DECL(i+i1,j+j1,k+k1),vofcompraw)

             ! in: subroutine FORT_RENORMALIZE_PRESCRIBE

            istate=1
            do while (istate.le.num_state_material)

             if (istate.eq.1) then
              dencomp=(im-1)*num_state_material+istate
               ! incompressible
              if (fort_material_type(im).eq.0) then 
               state_stencil(istate)=fort_denconst(im)
               ! compressible
              else if ((fort_material_type(im).gt.0).and. &
                       (fort_material_type(im).le.fort_max_num_eos)) then 
               state_stencil(istate)=den(D_DECL(i+i1,j+j1,k+k1),dencomp)
              else
               print *,"fort_material_type invalid"
               stop
              endif
              istate=istate+1
             else if (istate.eq.2) then
              tempcomp=(im-1)*num_state_material+istate
              state_stencil(istate)=den(D_DECL(i+i1,j+j1,k+k1),tempcomp)
              istate=istate+1
             else if ((istate.eq.num_state_base+1).and. &
                      (num_species_var.gt.0)) then 
              do ispecies=1,num_species_var
               speccomp=(im-1)*num_state_material+num_state_base+ispecies
               state_stencil(istate)=den(D_DECL(i+i1,j+j1,k+k1),speccomp)
               istate=istate+1
              enddo
             else 
              print *,"istate invalid"
              stop
             endif

            enddo ! do while (istate.le.num_state_material)

            F_stencil_sum=F_stencil_sum+F_stencil
            do istate=1,num_state_material
             state_stencil_sum(istate)= &
              state_stencil_sum(istate)+ &
              F_stencil*state_stencil(istate)
            enddo ! istate

           enddo !i1,j1,k1=-1..1 (init: state_stencil_sum, F_stencil_sum)
           enddo
           enddo

           if (F_stencil_sum.gt.VOFTOL) then
            do istate=1,num_state_material
             statecomp=(im-1)*num_state_material+istate
             dennew(D_DECL(i,j,k),statecomp)= &
               state_stencil_sum(istate)/F_stencil_sum
            enddo ! istate
           else if ((F_stencil_sum.ge.zero).and. &
                    (F_stencil_sum.le.VOFTOL)) then
            ! do nothing
           else
            print *,"F_stencil_sum invalid:",F_stencil_sum
            stop
           endif

          else
           print *,"F_stencil must be >= 0 F_stencil= ",F_stencil
           stop
          endif

         else
          print *,"is_rigid invalid"
          stop
         endif

        enddo ! im=1..nmat (extrapolation loop)

        do dir=1,nmat*ngeom_recon
         mofnew(dir)=zero
        enddo

        do im=1,nmat

         vofcomp=(im-1)*ngeom_recon+1
         vofcompraw=(im-1)*ngeom_raw+1

         mofnew(vofcomp)=vofnew(D_DECL(i,j,k),vofcompraw)
         do dir=1,SDIM
          mofnew(vofcomp+dir)=vofnew(D_DECL(i,j,k),vofcompraw+dir)
         enddo

         do istate=1,num_state_material
          statecomp=(im-1)*num_state_material+istate
          den_hold(statecomp)=dennew(D_DECL(i,j,k),statecomp)
         enddo

        enddo  ! im=1..nmat
       
         ! 1. prescribe solid materials. (F,X,LS,velocity,temperature)
         ! 2. extend the fluid level set functions into the solids.
         !   (F,X,LS fluid)
        if (renormalize_only.eq.0) then

         do im=1,nmat*(1+SDIM)
          ls_hold(im)=lsnew(D_DECL(i,j,k),im)
         enddo

         sum_vfrac_solid_new=zero
         max_solid_LS=-99999.0
         im_solid_max=0
         partid_max=0

         if ((nparts.lt.0).or.(nparts.gt.nmat)) then
          print *,"nparts invalid FORT_RENORMALIZE_PRESCRIBE"
          stop
         endif

         do partid=1,nparts

          im=im_solid_map(partid)+1
          if ((im.lt.1).or.(im.gt.nmat)) then
           print *,"im invalid33"
           stop
          endif

          if (is_lag_part(nmat,im).eq.0) then
           print *,"is_lag_part(nmat,im).eq.0"
           stop
          else if (is_lag_part(nmat,im).eq.1) then
           if (is_rigid(nmat,im).eq.0) then
            ! do nothing
           else if (is_rigid(nmat,im).eq.1) then

            ! positive in the rigid body
            call materialdistsolid( &
             xsten(0,1),xsten(0,2),xsten(0,SDIM), &
             LS_solid_new(im),time,im)

            if ((FSI_flag(im).eq.2).or. & ! prescribed solid CAD
                (FSI_flag(im).eq.4)) then ! CTML FSI
             LS_solid_new(im)=LS(D_DECL(i,j,k),im)
            else if (FSI_flag(im).eq.1) then ! prescribed solid EUL
             ! do nothing
            else
             print *,"FSI_flag invalid"
             stop
            endif
 
            if (LS_solid_new(im).gt.max_solid_LS) then
             max_solid_LS=LS_solid_new(im)
             im_solid_max=im
             partid_max=partid
            endif

            nrefine_geom=1

            ! centroid in absolute coordinate system
            ! returns a volume fraction
            call find_LS_stencil_volume( &
             bfact, &
             dx, &
             xsten,nhalf, &
             nrefine_geom, &
             time, &
             vfrac_solid_new(im), &
             centroid, &
             im)

            if ((FSI_flag(im).eq.2).or. & ! prescribed solid CAD
                (FSI_flag(im).eq.4)) then ! CTML FSI
             vofcompraw=(im-1)*ngeom_raw+1
             vfrac_solid_new(im)=vofnew(D_DECL(i,j,k),vofcompraw)
             do dir=1,SDIM
              centroid(dir)=vofnew(D_DECL(i,j,k),vofcompraw+dir)+cencell(dir)
             enddo
            else if (FSI_flag(im).eq.1) then ! prescribed solid EUL
             ! do nothing
            else
             print *,"FSI_flag invalid"
             stop
            endif

            do dir=1,SDIM
             censolid_new(im,dir)=centroid(dir)-cencell(dir)
            enddo

            ! slope points into the solid
            ! this slope ignores 1/R term for dphi/dtheta
            call find_LS_stencil_slope( &
             bfact, &
             dx, &
             xsten,nhalf, &
             nslope_solid, &
             time,im)

            if ((FSI_flag(im).eq.2).or. & ! prescribed solid CAD
                (FSI_flag(im).eq.4)) then ! CTML FSI
             do dir=1,SDIM
              nslope_solid(dir)=LS(D_DECL(i,j,k),nmat+SDIM*(im-1)+dir)
             enddo
            else if (FSI_flag(im).eq.1) then ! prescribed solid EUL
             ! do nothing
            else
             print *,"FSI_flag invalid"
             stop
            endif

            if (vfrac_solid_new(im).le.VOFTOL) then
             vfrac_solid_new(im)=zero
            else if (vfrac_solid_new(im).ge.one-VOFTOL) then
             vfrac_solid_new(im)=one
            else if ((vfrac_solid_new(im).ge.VOFTOL).and. &
                     (vfrac_solid_new(im).le.one-VOFTOL)) then
             ! do nothing
            else
             print *,"vfrac_solid_new bust"
             stop
            endif  
            
            sum_vfrac_solid_new=sum_vfrac_solid_new+vfrac_solid_new(im)

            ! level set for rigid solid
            ls_hold(im)=LS_solid_new(im)
            do dir=1,SDIM
             ls_hold(nmat+SDIM*(im-1)+dir)=nslope_solid(dir)
            enddo

            if ((vfrac_solid_new(im).le.VOFTOL).and. &
                (LS_solid_new(im).ge.zero)) then
             print *,"cannot have F<eps and LS>0"
             print *,"i,j,k,im ",i,j,k,im
             print *,"vfrac_solid_new(im)=",vfrac_solid_new(im)
             print *,"LS_solid_new(im)=",LS_solid_new(im)
             stop
            endif
            if ((vfrac_solid_new(im).ge.one-VOFTOL).and. &
                (LS_solid_new(im).le.zero)) then
             print *,"cannot have F>1-eps and LS<0"
             print *,"i,j,k,im ",i,j,k,im
             print *,"vfrac_solid_new(im)=",vfrac_solid_new(im)
             print *,"LS_solid_new(im)=",LS_solid_new(im)
             stop
            endif

            ! temperature in rigid solid.
            if ((LS_solid_new(im).ge.zero).or. &
                (vfrac_solid_new(im).ge.half)) then
             istate=2
             statecomp=(im-1)*num_state_material+istate
             if (solidheat_flag.eq.0) then
              ! do nothing (heat conduction in solid)
             else if (solidheat_flag.eq.2) then ! neumann at solid/fluid
              if ((FSI_flag(im).eq.2).or. & ! prescribed solid CAD
                  (FSI_flag(im).eq.4)) then ! CTML FSI
               ! den_hold(statecomp) already has the solid temperature
              else if (FSI_flag(im).eq.1) then ! prescribed solid EUL
               call tempsolid(xsten(0,1),xsten(0,2),xsten(0,SDIM), &
                den_hold(statecomp),time,im)
              else
               print *,"FSI_flag invalid"
               stop
              endif
             else if (solidheat_flag.eq.1) then ! dirichlet at solid/fluid
              if ((FSI_flag(im).eq.2).or. & ! prescribed solid CAD
                  (FSI_flag(im).eq.4)) then ! CTML FSI
               ! den_hold(statecomp) already has the solid temperature
              else if (FSI_flag(im).eq.1) then ! prescribed solid EUL
               call tempsolid(xsten(0,1),xsten(0,2),xsten(0,SDIM), &
                den_hold(statecomp),time,im)
              else
               print *,"FSI_flag invalid"
               stop
              endif
             else
              print *,"solidheat_flag out of range"
              stop
             endif
            else if ((LS_solid_new(im).le.zero).and. &
                     (vfrac_solid_new(im).le.half)) then
             ! do nothing
            else
             print *,"LS_solid_new or vfrac_solid_new invalid"
             stop
            endif

           else
            print *,"is_rigid invalid"
            stop
           endif
          else 
           print *,"is_lag_part(nmat,im) invalid"
           stop
          endif

         enddo ! partid=1..nparts

           ! velocity and temperature in rigid solid.
         if ((sum_vfrac_solid_new.ge.half).or. &
             (max_solid_LS.ge.zero)) then 

          if ((im_solid_max.lt.1).or. &
              (im_solid_max.gt.nmat).or. &
              (partid_max.lt.1).or. &
              (partid_max.gt.nparts).or. &
              (is_rigid(nmat,im_solid_max).ne.1)) then
           print *,"im_solid_max or partid_max became corrupt"
           stop
          endif

           ! solid velocity
          ibase=(partid_max-1)*SDIM

           ! velocity
          if (is_prescribed(nmat,im_solid_max).eq.1) then
           do dir=1,SDIM
            velnew(D_DECL(i,j,k),dir)=solid(D_DECL(i,j,k),ibase+dir)
           enddo
          else if (is_prescribed(nmat,im_solid_max).eq.0) then
           ! do nothing
          else
           print *,"is_prescribed(nmat,im_solid_max) invalid"
           stop
          endif

           ! initialize the fluid temperature with the solid temperature
           ! in the solid regions.
          do im=1,nmat
           if (is_rigid(nmat,im).eq.1) then
            ! do nothing
           else if (is_rigid(nmat,im).eq.0) then
            istate=2
            statecomp=(im-1)*num_state_material+istate
            statecomp_solid=(im_solid_max-1)*num_state_material+istate
            den_hold(statecomp)=den_hold(statecomp_solid)
           else
            print *,"is_rigid invalid"
            stop
           endif
          enddo ! im=1..nmat

         else if ((sum_vfrac_solid_new.le.half).and. &
                  (max_solid_LS.le.zero)) then
          ! do nothing
         else
          print *,"sum_vfrac_solid_new or max_solid_LS invalid"
          stop
         endif

         do im=1,nmat
          do istate=1,num_state_material
           statecomp=(im-1)*num_state_material+istate
           dennew(D_DECL(i,j,k),statecomp)=den_hold(statecomp)
          enddo
         enddo

          ! solid volume fractions and centroids
         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          if (is_rigid(nmat,im).eq.0) then
           ! do nothing
          else if (is_rigid(nmat,im).eq.1) then
           mofnew(vofcomp)=vfrac_solid_new(im)
           do dir=1,SDIM
            mofnew(vofcomp+dir)=censolid_new(im,dir)
           enddo
           do istate=SDIM+2,ngeom_recon
            mofnew(vofcomp+istate-1)=zero
           enddo
          else
           print *,"is_rigid invalid"
           stop
          endif
         enddo ! im=1..nmat

          ! extend fluid LS,F,X into the solid.
         if ((im_solid_max.ge.1).and.(im_solid_max.le.nmat)) then

          if ((LS_solid_new(im_solid_max).ge.zero).or. &
              (sum_vfrac_solid_new.ge.half)) then

            ! default radius: 1 cell
            ! this inner loop is needed since the volume fraction
            ! at (i,j,k) depends on the levelset function values
            ! in the (i+i1,j+j1,k+k1) stencil.
           do i1=istenlo(1),istenhi(1)
           do j1=istenlo(2),istenhi(2)
           do k1=istenlo(3),istenhi(3)

            near_fluid=zero

            if ((i1.eq.0).and.(j1.eq.0).and.(k1.eq.0)) then
             at_center=1
            else
             at_center=0
            endif

            do im=1,nmat
             LS_predict(im)=LS(D_DECL(i+i1,j+j1,k+k1),im)
            enddo
            call get_primary_material(LS_predict,nmat,im_primary_stencil)

             !fluid stencil cell
            if ((is_rigid(nmat,im_primary_stencil).eq.0).and. & 
                (at_center.eq.0)) then

             do im=1,nmat
              LS_virtual_new(im)=LS_predict(im)
             enddo

             !solid stencil cell
            else if ((is_rigid(nmat,im_primary_stencil).eq.1).or. & 
                     (at_center.eq.1)) then

             do im=1,nmat
              LS_virtual_new(im)=zero
             enddo

              ! default radius: 2 cells.
              ! In order to find the level set function values at the
              ! (i+i1,j+j1,k+k1) stencil points in the solid, 
              ! extrapolation from the
              ! fluid region must be done.
             do i2=LSstenlo(1),LSstenhi(1)
             do j2=LSstenlo(2),LSstenhi(2)
             do k2=LSstenlo(3),LSstenhi(3)
     
              do im=1,nmat
               LS_virtual(im)=LS(D_DECL(i+i1+i2,j+j1+j2,k+k1+k2),im)
              enddo
               ! the fluid cells closest to the target cell have the
               ! most weight.
              call get_primary_material(LS_virtual,nmat,im_primary_sub_stencil)
              if (is_rigid(nmat,im_primary_sub_stencil).eq.0) then
                ! WT=1/((dx/10)^2 + ||xfluid-xghost||^2)
               WT=0.01*(dxmaxLS**2)
                ! nhalf==7
                ! xfluid: fluid cell
                ! xghost: solid cell
               do dir=1,SDIM
                if (dir.eq.1) then
                 xfluid=xsten(2*(i1+i2),dir)
                 xghost=xsten(2*i1,dir)
                else if (dir.eq.2) then 
                 xfluid=xsten(2*(j1+j2),dir)
                 xghost=xsten(2*j1,dir)
                else if ((dir.eq.3).and.(SDIM.eq.3)) then
                 xfluid=xsten(2*(k1+k2),dir)
                 xghost=xsten(2*k1,dir)
                else
                 print *,"dir invalid"
                 stop
                endif

                WT=WT+(xfluid-xghost)**2
               enddo ! dir=1..sdim

               if (WT.gt.zero) then
                ! do nothing
               else
                print *,"WT invalid"
                stop
               endif
               WT=one/WT
 
               near_fluid=near_fluid+WT
               do im=1,nmat
                LS_virtual_new(im)=LS_virtual_new(im)+WT*LS_virtual(im)
               enddo
              else if (is_rigid(nmat,im_primary_sub_stencil).eq.1) then
               ! do nothing
              else
               print *,"is_rigid(nmat,im_primary_sub_stencil) invalid"
               stop 
              endif
          
             enddo 
             enddo 
             enddo  !i2,j2,k2=-2..2
           
             if (near_fluid.gt.zero) then
              do im=1,nmat
               LS_virtual_new(im)=LS_virtual_new(im)/near_fluid
              enddo
             else if (near_fluid.eq.zero) then
              ! LS_predict=LS(i+i1,j+j1,k+k1)
              do im=1,nmat
               LS_virtual_new(im)=LS_predict(im)
              enddo
             else
              print *,"near_fluid invalid  near_fluid=",near_fluid
              stop
             endif
            else
             print *,"is_rigid(nmat,im_primary_stencil) or at_center invalid"
             stop
            endif

            do im=1,nmat
             LS_extend(D_DECL(i1,j1,k1),im)=LS_virtual_new(im)
            enddo
                
           enddo
           enddo
           enddo ! i1,j1,k1=-1..1

            ! default radius: 1 cell
            ! make the extension fluid level set tessellating
           do i1=istenlo(1),istenhi(1)
           do j1=istenlo(2),istenhi(2)
           do k1=istenlo(3),istenhi(3)

            do im=1,nmat
             LS_virtual(im)=LS_extend(D_DECL(i1,j1,k1),im)
             LS_virtual_new(im)=LS_virtual(im)
            enddo
 
            do im=1,nmat
             if (is_rigid(nmat,im).eq.1) then
              ! do nothing
             else if (is_rigid(nmat,im).eq.0) then

              LS_virtual_max=-99999.0
              im_primary_stencil=0
              do im_opp=1,nmat
               if (is_rigid(nmat,im_opp).eq.1) then
                ! do nothing
               else if (is_rigid(nmat,im_opp).eq.0) then
                if (im.ne.im_opp) then
                 if (im_primary_stencil.eq.0) then
                  im_primary_stencil=im_opp
                  LS_virtual_max=LS_virtual(im_opp)
                 else if ((im_primary_stencil.ge.1).and. &
                          (im_primary_stencil.le.nmat)) then
                  if (LS_virtual(im_opp).gt.LS_virtual_max) then
                   im_primary_stencil=im_opp
                   LS_virtual_max=LS_virtual(im_opp)
                  endif
                 else
                  print *,"im_primary_stencil invalid"
                  stop
                 endif
                else if (im.eq.im_opp) then
                 ! do nothing
                else
                 print *,"im_opp invalid"
                 stop
                endif
               else
                print *,"is_rigid invalid"
                stop
               endif
              enddo ! im_opp=1..nmat
 
              if (im_primary_stencil.eq.0) then
               if (nmat_fluid.ne.1) then
                print *,"nmat_fluid invalid"
                stop
               endif
               if (LS_virtual(im).le.zero) then
                print *,"LS_virtual(im) invalid"
                stop
               endif
              else if ((im_primary_stencil.ge.1).and. &
                       (im_primary_stencil.le.nmat)) then
               LS_virtual_new(im)=half*(LS_virtual(im)-LS_virtual_max)
              else
               print *,"im_primary_stencil invalid"
               stop
              endif

             else
              print *,"is_rigid invalid"
              stop
             endif

            enddo ! im=1..nmat

            do im=1,nmat
             LS_extend(D_DECL(i1,j1,k1),im)=LS_virtual_new(im)
            enddo

           enddo
           enddo
           enddo ! i1,j1,k1=-1..1
               
           do im=1,nmat

            if (is_rigid(nmat,im).eq.1) then
             ! do nothing
            else if (is_rigid(nmat,im).eq.0) then
             do i1=istenlo(1),istenhi(1)
             do j1=istenlo(2),istenhi(2)
             do k1=istenlo(3),istenhi(3)
              LS_temp(D_DECL(i1,j1,k1))=LS_extend(D_DECL(i1,j1,k1),im)
             enddo 
             enddo 
             enddo 

             vofcomp=(im-1)*ngeom_recon+1
              
             call getvolume(bfact,dx,xsten,nhalf, &
              LS_temp,mofnew(vofcomp),LSfacearea, &
              LScentroid,LSareacentroid,VOFTOL,SDIM)

             do dir=1,SDIM
              mofnew(vofcomp+dir)=LScentroid(dir)-cencell(dir)
             enddo

             ls_hold(im)=LS_temp(D_DECL(0,0,0))
 
             if (mofnew(vofcomp).lt.zero) then
              print *,"mofnew(vofcomp) invalid"
              stop
             else if (mofnew(vofcomp).le.VOFTOL) then
              if (ls_hold(im).ge.zero) then
               mofnew(vofcomp)=VOFTOL_SLOPES
              endif
             else if ((mofnew(vofcomp).gt.zero).and. &
                      (mofnew(vofcomp).le.one+VOFTOL)) then
              ! do nothing
             else
              print *,"mofnew(vofcomp) invalid"
              stop
             endif
            else
             print *,"is_rigid invalid"
             stop
            endif

           enddo ! im=1..nmat

          else if ((LS_solid_new(im_solid_max).le.zero).and. &
                   (sum_vfrac_solid_new.le.half)) then
           ! do nothing
          else
           print *,"LS_solid_new or sum_vfrac_Solid_new invalid"
           stop
          endif

         else if (im_solid_max.eq.0) then
          ! do nothing
         else
          print *,"im_solid_max invalid"
          stop
         endif

          ! sum of F_fluid=1
          ! sum of F_rigid<=1
         call make_vfrac_sum_ok_base(tessellate,mofnew,nmat,SDIM,12)

         do im=1,nmat*(1+SDIM)
          lsnew(D_DECL(i,j,k),im)=ls_hold(im)
         enddo
 
         ! above: prescribe new solid location
         ! below: just normalize the volume fractions.
        else if (renormalize_only.eq.1) then

         if (1.eq.0) then
          print *,"i,j,k ",i,j,k
         endif

         call make_vfrac_sum_ok_base(tessellate,mofnew,nmat,SDIM,13)
        else
         print *,"renormalize only invalid"
         stop
        endif

        if (ngeom_raw.ne.SDIM+1) then
         print *,"ngeom_raw invalid"
         stop
        endif

        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         vofcompraw=(im-1)*ngeom_raw+1
         vofnew(D_DECL(i,j,k),vofcompraw)=mofnew(vofcomp)
         do dir=1,SDIM
          vofnew(D_DECL(i,j,k),vofcompraw+dir)=mofnew(vofcomp+dir)
         enddo
        enddo ! im

       else if (local_maskcov.eq.0) then
        ! do nothing
       else
        print *,"local_maskcov invalid"
        stop
       endif

      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine FORT_RENORMALIZE_PRESCRIBE



      subroutine FORT_PURGEFLOTSAM( &
       delta_mass, &
       truncate_volume_fractions, &
       truncate_thickness, &
       level,finest_level, &
       time, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       maskcov,DIMS(maskcov), &
       vofnew,DIMS(vofnew), &
       LS,DIMS(LS), &
       xlo,dx,nmat)
      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE


      INTEGER_T level,finest_level

      REAL_T time
      REAL_T xlo(SDIM)
      REAL_T dx(SDIM)
      INTEGER_T nmat
      REAL_T delta_mass(nmat)
      INTEGER_T DIMDEC(maskcov)
      INTEGER_T DIMDEC(vofnew)
      INTEGER_T DIMDEC(LS)
      INTEGER_T truncate_volume_fractions(nmat)
      REAL_T  maskcov(DIMV(maskcov))
      REAL_T  vofnew(DIMV(vofnew),nmat*ngeom_raw)
      REAL_T  LS(DIMV(LS),nmat)
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      INTEGER_T i,j,k,dir
      INTEGER_T im,im2,imcrit
      INTEGER_T vofcomprecon,vofcompraw
      REAL_T truncate_thickness

      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T volmat(nmat)
      REAL_T lspoint(nmat)
      INTEGER_T sorted_list(nmat)
      REAL_T dxmax,dxmaxLS,LSbandsize,restore_sum
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      REAL_T volgrid
      REAL_T cengrid(SDIM)
      INTEGER_T mask_test
      INTEGER_T FSI_exclude
      INTEGER_T tessellate

      tessellate=0

      nhalf=3
 
      if (bfact.lt.1) then
       print *,"bfact invalid103"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat bust"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (level.lt.0) then
       print *,"level invalid purge flotsam 1"
       stop
      endif
      if (level.gt.finest_level) then
       print *,"finest_level invalid purge flotsam 2"
       stop
      endif

      if (truncate_thickness.lt.one) then
       print *,"truncate_thickness too small"
       stop
      endif

      if (num_state_material.ne. &
          num_state_base+num_species_var) then
       print *,"num_state_material invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(vofnew),1,-1,26)
      call checkbound(fablo,fabhi,DIMS(LS),1,-1,26)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      call get_dxmax(dx,bfact,dxmax)
      call get_dxmaxLS(dx,bfact,dxmaxLS)
      LSbandsize=truncate_thickness*dxmaxLS

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       mask_test=NINT(maskcov(D_DECL(i,j,k)))

       if (mask_test.eq.1) then

        call gridsten_level(xsten,i,j,k,level,nhalf)
        call Box_volumeFAST(bfact,dx,xsten,nhalf,volgrid,cengrid,SDIM)
        if (volgrid.le.zero) then
         print *,"volgrid invalid"
         stop
        endif
        do dir=1,nmat*ngeom_recon
         mofdata(dir)=zero 
        enddo
        do im=1,nmat
         vofcomprecon=(im-1)*ngeom_recon+1
         vofcompraw=(im-1)*ngeom_raw+1
         volmat(im)=vofnew(D_DECL(i,j,k),vofcompraw)
         lspoint(im)=LS(D_DECL(i,j,k),im)
         mofdata(vofcomprecon)=volmat(im)
         do dir=1,SDIM
          mofdata(vofcomprecon+dir)=vofnew(D_DECL(i,j,k),vofcompraw+dir)
         enddo
        enddo ! im
        FSI_exclude=1
        call sort_volume_fraction(volmat,FSI_exclude,sorted_list,nmat)
        imcrit=sorted_list(1)
        if (is_rigid(nmat,imcrit).ne.0) then
         print *,"is_rigid(nmat,imcrit) invalid"
         stop
        endif
        do im=1,nmat
         if (is_rigid(nmat,im).eq.0) then

          if (lspoint(im).gt.LSbandsize) then
           vofcomprecon=(im-1)*ngeom_recon+1
           mofdata(vofcomprecon)=one
           do dir=1,SDIM
            mofdata(vofcomprecon+dir)=zero
           enddo

           restore_sum=zero

           do im2=1,nmat
            if (is_rigid(nmat,im2).eq.0) then

             if (im2.ne.im) then
              vofcomprecon=(im2-1)*ngeom_recon+1

              if (truncate_volume_fractions(im2).eq.1) then

               mofdata(vofcomprecon)=zero
               do dir=1,SDIM
                mofdata(vofcomprecon+dir)=zero
               enddo

              else if (truncate_volume_fractions(im2).eq.0) then
               restore_sum=restore_sum+mofdata(vofcomprecon)
              else
               print *,"truncate_volume_fractions invalid"
               stop
              endif

             endif ! im2<>im

            else if (is_rigid(nmat,im2).eq.1) then
             ! do nothing
            else
             print *,"is_rigid(nmat,im2) invalid"
             stop
            endif
           enddo ! im2

           vofcomprecon=(im-1)*ngeom_recon+1
           mofdata(vofcomprecon)=one-restore_sum
          
          else if (lspoint(im).lt.-LSbandsize) then
           if (im.ne.imcrit) then

            if (truncate_volume_fractions(im).eq.1) then
             vofcomprecon=(im-1)*ngeom_recon+1

             mofdata(vofcomprecon)=zero
             do dir=1,SDIM
              mofdata(vofcomprecon+dir)=zero
             enddo
            else if (truncate_volume_fractions(im).eq.0) then
             ! do nothing
            else
             print *,"truncate_volume_fractions invalid"
             stop
            endif
           else if (im.eq.imcrit) then
            ! do nothing
           else
            print *,"im invalid34"
            stop
           endif
          else if ((lspoint(im).ge.-LSbandsize).and. &
                   (lspoint(im).le.LSbandsize)) then
           ! do nothing
          else
           print *,"lspoint bust"
           stop
          endif 

         else if (is_rigid(nmat,im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(nmat,im) invalid"
          stop
         endif
        enddo ! im=1...nmat
         ! sum F_fluid=1  sum F_solid <=1
        call make_vfrac_sum_ok_base(tessellate,mofdata,nmat,SDIM,14)
        do im=1,nmat
         vofcompraw=(im-1)*ngeom_raw+1
         vofcomprecon=(im-1)*ngeom_recon+1
         vofnew(D_DECL(i,j,k),vofcompraw)=mofdata(vofcomprecon)
         do dir=1,SDIM
          vofnew(D_DECL(i,j,k),vofcompraw+dir)=mofdata(vofcomprecon+dir)
         enddo

         delta_mass(im)=delta_mass(im)+ &
           volgrid*(mofdata(vofcomprecon)-volmat(im))
        enddo ! im=1..nmat

       else if (mask_test.eq.0) then
        ! do nothing
       else
        print *,"mask_test invalid"
        stop
       endif

      enddo
      enddo
      enddo

      return
      end subroutine FORT_PURGEFLOTSAM



      subroutine FORT_INITRECALESCE( &
       recalesce_material_in, &
       recalesce_state_old_in, &
       recalesce_num_state_in,nmat)
      use probcommon_module
      use probf90_module

      IMPLICIT NONE


      INTEGER_T recalesce_num_state_in,nmat
      INTEGER_T i
      INTEGER_T recalesce_material_in(nmat)
      REAL_T recalesce_state_old_in(recalesce_num_state*nmat)

      if (nmat.ne.num_materials) then
       print *,"nmat bust"
       stop
      endif
      if (nmat.gt.100) then
       print *,"too many materials"
       stop
      endif
      if (recalesce_num_state.ne.recalesce_num_state_in) then
       print *,"recalesce_num_state_in invalid"
       stop
      endif

      do i=1,nmat
       recalesce_material(i)=recalesce_material_in(i)
      enddo
      do i=1,recalesce_num_state*nmat
       recalesce_state_old(i)=recalesce_state_old_in(i)
      enddo

      return
      end subroutine FORT_INITRECALESCE



