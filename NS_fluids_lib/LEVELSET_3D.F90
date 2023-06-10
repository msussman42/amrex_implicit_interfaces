! in .vmrc file in home dir: set tabstop=1, set shiftwidth=1, set expandtab
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

! N_EXTRA_REAL.H is in the amrlib directory.
#include "N_EXTRA_REAL.H"
#include "EXTRAP_COMP.H"
#include "LEVEL_F.H"

#define nsum 64
#define nsum2 32
#define CURVWT (1.0D-3)

#define DEBUG_THERMAL_COEFF 0
#define DEBUG_CURVATURE 0

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

        type cell_CP_parm_type
         REAL_T dxmaxLS
         INTEGER_T i,j,k
         INTEGER_T bfact
         INTEGER_T level
         INTEGER_T finest_level
         INTEGER_T, pointer :: fablo(:)
         INTEGER_T, pointer :: fabhi(:)
         REAL_T, pointer :: dx(:)
         REAL_T :: time
         INTEGER_T :: im_solid_max
         INTEGER_T :: least_sqr_radius
         INTEGER_T :: least_sqrZ
         REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LS

        end type cell_CP_parm_type

        contains

        subroutine cell_xCP(cell_CP_parm,xCP,xSOLID_BULK)
        use global_utility_module
        use global_distance_module
        use probf90_module
        use geometry_intersect_module
        use MOF_routines_module
        IMPLICIT NONE
        type(cell_CP_parm_type), INTENT(in) :: cell_CP_parm
        REAL_T, INTENT(out) :: xCP(SDIM)
        REAL_T, INTENT(out) :: xSOLID_BULK(SDIM)
        INTEGER_T, parameter :: nhalf=3
        REAL_T :: xsten(-nhalf:nhalf,SDIM)
        INTEGER_T :: dir
        REAL_T :: nslope_cell(SDIM)
        REAL_T :: LS_cell
        REAL_T :: mag

        ASSOCIATE(CP=>cell_CP_parm)

        if ((CP%im_solid_max.ge.1).and. &
            (CP%im_solid_max.le.num_materials)) then

         if (is_rigid(CP%im_solid_max).eq.1) then

          call gridsten_level(xsten,CP%i,CP%j,CP%k,CP%level,nhalf)

           ! positive in the rigid body
          call materialdistsolid( &
            xsten(0,1),xsten(0,2),xsten(0,SDIM), &
            LS_cell,CP%time,CP%im_solid_max)

           ! xCP=x-phi grad phi   grad phi=(x-xCP)/phi
           ! slope points into the solid
           ! this slope ignores 1/R term for dphi/dtheta
          call find_LS_stencil_slope( &
            CP%bfact, &
            CP%dx, &
            xsten,nhalf, &
            nslope_cell, &
            CP%time, &
            CP%im_solid_max)

          if ((FSI_flag(CP%im_solid_max).eq.FSI_PRESCRIBED_NODES).or. & 
              (FSI_flag(CP%im_solid_max).eq.FSI_SHOELE_PRESVEL).or. & 
              (FSI_flag(CP%im_solid_max).eq.FSI_SHOELE_VELVEL)) then 
           LS_cell=CP%LS(D_DECL(CP%i,CP%j,CP%k),CP%im_solid_max)
           do dir=1,SDIM
            nslope_cell(dir)= &
              CP%LS(D_DECL(CP%i,CP%j,CP%k), &
              num_materials+SDIM*(CP%im_solid_max-1)+dir)
           enddo
          else if (FSI_flag(CP%im_solid_max).eq.FSI_PRESCRIBED_PROBF90) then 
           ! do nothing
          else
           print *,"FSI_flag invalid"
           stop
          endif

          mag=zero
          do dir=1,SDIM
           xSOLID_BULK(dir)=xsten(0,dir)
           xCP(dir)=xSOLID_BULK(dir)-LS_cell*nslope_cell(dir)
           mag=mag+nslope_cell(dir)**2
          enddo
          mag=sqrt(mag)
          if (abs(mag-one).lt.VOFTOL) then
           ! do nothing
          else
           print *,"nslope_cell (mag) invalid"
           stop
          endif

         else 
          print *,"is_rigid(CP%im_solid_max) invalid"
          stop
         endif
        else
         print *,"CP%im_solid_max invalid"
         stop
        endif

        END ASSOCIATE

        return 
        end subroutine cell_xCP


        subroutine interp_fluid_LS( &
         cell_CP_parm, &
         xCP, &
         xSOLID_BULK, &
         cell_index, &
         LS_interp, &
         im_solid, &
         im_fluid_critical)
        use global_utility_module
        use global_distance_module
        use probf90_module
        IMPLICIT NONE
        type(cell_CP_parm_type), INTENT(in) :: cell_CP_parm
        REAL_T, INTENT(in) :: xCP(SDIM)
        REAL_T, INTENT(in) :: xSOLID_BULK(SDIM)
        INTEGER_T, INTENT(in) :: cell_index(SDIM)
        REAL_T, INTENT(out) :: LS_interp(num_materials)
        INTEGER_T, INTENT(in) :: im_solid 
        INTEGER_T, INTENT(out) :: im_fluid_critical
        INTEGER_T, parameter :: nhalf=1
        REAL_T :: xsten(-nhalf:nhalf,SDIM)
        INTEGER_T :: dir
        INTEGER_T :: local_index(SDIM)
        INTEGER_T :: i2,j2,k2
        INTEGER_T :: i_safe,j_safe,k_safe
        INTEGER_T :: isten,jsten,ksten
        REAL_T :: LS_virtual(num_materials)
        REAL_T :: LS_center_stencil(num_materials)
        INTEGER_T im_local
        INTEGER_T im_primary_sub_stencil
        REAL_T shortest_dist_to_fluid
        REAL_T dist_stencil_to_bulk
        REAL_T :: LS_interp_low_order(num_materials)
        REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: local_data_fab
        REAL_T local_data_out

        INTEGER_T LSstenlo(3)
        INTEGER_T LSstenhi(3)

        ASSOCIATE(CP=>cell_CP_parm)

        LSstenlo(3)=0
        LSstenhi(3)=0

         ! cell_index is containing cell for xCP
        do dir=1,SDIM
         local_index(dir)=cell_index(dir)
         LSstenlo(dir)=-1
         LSstenhi(dir)=1
        enddo ! dir=1..sdim

        shortest_dist_to_fluid=-one
        im_fluid_critical=0
        local_data_fab=>CP%LS

         ! stencil radius is 1.
        do i2=LSstenlo(1),LSstenhi(1)
        do j2=LSstenlo(2),LSstenhi(2)
        do k2=LSstenlo(3),LSstenhi(3)

          ! local_index is the containing cell for xCP
         isten=local_index(1)+i2
         jsten=local_index(2)+j2
         ksten=local_index(SDIM)+k2

         call safe_data_index(isten,jsten,ksten,i_safe,j_safe,k_safe, &
                 local_data_fab)

         call gridsten_level(xsten, &
          i_safe,j_safe,k_safe, &
          CP%level,nhalf)

         dist_stencil_to_bulk=zero

          ! xCP=xSOLID_BULK(dir)-LS_cell*nslope_cell(dir)
          ! xSOLID_BULK usually in the solid, but it might be
          ! in the fluid, at most 1 cell away from a solid cell.
          ! NOTE: the output from this routine is ignored if xSOLID_BULK
          ! in a fluid cell.
         do dir=1,SDIM
          dist_stencil_to_bulk=dist_stencil_to_bulk+ &
                  (xsten(0,dir)-xSOLID_BULK(dir))**2
         enddo
         dist_stencil_to_bulk=sqrt(dist_stencil_to_bulk)

         do im_local=1,num_materials
          call safe_data(isten,jsten,ksten,im_local, &
           local_data_fab,local_data_out)
          LS_virtual(im_local)=local_data_out
          if ((i2.eq.0).and.(j2.eq.0).and.(k2.eq.0)) then
           LS_center_stencil(im_local)=local_data_out
          endif
         enddo !im_local=1,num_materials

         ! the fluid cells closest to the substrate, but not
         ! in the substrate, have the most weight.
         call get_primary_material(LS_virtual,im_primary_sub_stencil)

         if (is_rigid(im_primary_sub_stencil).eq.0) then

          if (shortest_dist_to_fluid.eq.-one) then
           im_fluid_critical=im_primary_sub_stencil
           shortest_dist_to_fluid=dist_stencil_to_bulk
           do im_local=1,num_materials
            LS_interp_low_order(im_local)=LS_virtual(im_local)
           enddo
          else if (shortest_dist_to_fluid.ge.zero) then
           if (dist_stencil_to_bulk.lt. &
               shortest_dist_to_fluid) then
            im_fluid_critical=im_primary_sub_stencil
            shortest_dist_to_fluid=dist_stencil_to_bulk
            do im_local=1,num_materials
             LS_interp_low_order(im_local)=LS_virtual(im_local)
            enddo
           else if (dist_stencil_to_bulk.ge. &
                    shortest_dist_to_fluid) then
            ! do nothing
           else
            print *,"dist_stencil_bulk invalid"
            stop
           endif
          else
           print *,"shortest_dist_to_fluid invalid"
           stop
          endif

         else if (is_rigid(im_primary_sub_stencil).eq.1) then

          ! do nothing

         else
          print *,"is_rigid(im_primary_sub_stencil) invalid"
          stop 
         endif
         
        enddo 
        enddo 
        enddo  !i2,j2,k2=LSstenlo ... LSstenhi
       
        if (shortest_dist_to_fluid.ge.zero) then

         do im_local=1,num_materials
          LS_interp(im_local)=LS_interp_low_order(im_local)
         enddo
          !no fluid cells in stencil
        else if (shortest_dist_to_fluid.eq.-one) then 
          !no fluid cells in stencil
         do im_local=1,num_materials
          LS_interp(im_local)=LS_center_stencil(im_local)
         enddo

        else
         print *,"shortest_dist_to_fluid invalid"
         stop
        endif


        END ASSOCIATE

        return 
        end subroutine interp_fluid_LS



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
 

       ! called from fort_curvstrip
      subroutine initheightLS( &
        static_flag, & !intent(in)
        vof_height_function, & !intent(in)
        icenter,jcenter,kcenter, & !intent(in)
        level, & ! intent(in)
        finest_level, & !intent(in)
        bfact, & !intent(in)
        dx, & !intent(in)
        xcenter, & !intent(in)
        nrmcenter, & !intent(in)  sdim x num_materials components
        dircrit, & !intent(in)  1..sdim
        side, & !intent(in)  -1 or 1
        signside, & !intent(in)
        time, & !intent(in)
        xsten, & !intent(in)
        velsten, & !intent(in)
        mgoni_temp, & !intent(in)
        lssten, & !intent(in)
        vofsten, & !intent(in)
        nrmsten, & !intent(in)
        vol_sten, & !intent(in)
        area_sten, & !intent(in)
        curvHT_choice, & !intent(out)
        curvFD, & !intent(out)
        mgoni_force, & !intent(out)
        ZEYU_thet_d, & !intent(out)
        ZEYU_u_cl, & !intent(out)
        im3, & !intent(out)
        visc_coef, & !intent(in)
        unscaled_min_curvature_radius, & !intent(in)
        im, & !intent(in)
        im_opp, & !intent(in)
        iten) !intent(in)
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: static_flag
      INTEGER_T, INTENT(in) :: vof_height_function
      INTEGER_T, INTENT(in) :: icenter,jcenter,kcenter
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: bfact
      REAL_T, INTENT(in) :: xcenter(SDIM)
      INTEGER_T, INTENT(in) :: dircrit ! 1..SDIM
      INTEGER_T, INTENT(in) :: side ! 1 or -1
      INTEGER_T, INTENT(in) :: signside
      REAL_T, INTENT(in) :: time
      INTEGER_T, INTENT(in) :: im,im_opp
      INTEGER_T, INTENT(out) :: im3
      INTEGER_T :: imhold
      INTEGER_T :: im3_present_node
      INTEGER_T, INTENT(in) :: iten
      INTEGER_T :: iten_test
      REAL_T, INTENT(in) :: visc_coef
      REAL_T, INTENT(in) :: unscaled_min_curvature_radius
      REAL_T user_tension(num_interfaces)
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: vol_sten
      REAL_T, INTENT(in) :: area_sten(SDIM,2)
      REAL_T :: curvHT_LS
      REAL_T :: curvHT_VOF
      REAL_T, INTENT(out) :: curvHT_choice
      REAL_T, INTENT(out) :: curvFD
      REAL_T, INTENT(out) :: mgoni_force(SDIM)
 
      INTEGER_T dir2

      REAL_T columnLS(-ngrow_distance:ngrow_distance)
      REAL_T columnVOF(-ngrow_distance:ngrow_distance)

      REAL_T lsdata( &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance)
      REAL_T vofdata( &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance)

      REAL_T htfunc_LS(-1:1,-1:1)
      REAL_T htfunc_VOF(-1:1,-1:1)

      REAL_T, INTENT(in) :: xsten( &
       -(2*ngrow_distance+1):(2*ngrow_distance+1), &
       SDIM)
              
      REAL_T xsten_curv(-2:2,SDIM)

      REAL_T, INTENT(in) :: velsten( &
       -1:1,-1:1,-1:1,SDIM)

      REAL_T, INTENT(in) :: mgoni_temp( &
       -1:1,-1:1,-1:1,num_materials)

      REAL_T :: mgoni_tension(-1:1,-1:1,-1:1)
      REAL_T :: local_tension(num_interfaces)

      REAL_T, INTENT(in) :: lssten( &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        num_materials)

      REAL_T, INTENT(in) :: vofsten( &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        num_materials)

      REAL_T, INTENT(in) :: nrmsten( &
       -1:1,-1:1,-1:1,SDIM*num_materials)

      REAL_T, INTENT(in) :: nrmcenter(SDIM*num_materials)
      REAL_T nrmtest(SDIM*num_materials)

      REAL_T ngrid(SDIM)

      INTEGER_T itanlo,itanhi,jtanlo,jtanhi
      INTEGER_T itan,jtan
      INTEGER_T iwidth,jwidth,kheight
      INTEGER_T lmin,lmax
      REAL_T xtop,xbottom
      INTEGER_T iofs,jofs,kofs
      INTEGER_T iwidthnew
      INTEGER_T icell,jcell,kcell
      INTEGER_T inode,jnode,knode
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      INTEGER_T i2,j2,k2
      INTEGER_T node_index(3)
      INTEGER_T cell_index(3)

      REAL_T dx_col(SDIM)
      REAL_T x_col(SDIM)
      REAL_T x_col_avg(SDIM)
     
      REAL_T col_ht_LS
      REAL_T col_ht_VOF

      REAL_T n1d

      INTEGER_T iten_13,iten_23
      REAL_T gamma1,gamma2
      REAL_T cos_angle,sin_angle

      INTEGER_T klo_sten_short,khi_sten_short
      INTEGER_T klo_sten_ht,khi_sten_ht

      REAL_T dotprod,udotn,totaludotn
      REAL_T totalwt,wt,wtnode
      REAL_T liquid_viscosity
      REAL_T nproject(SDIM)
      INTEGER_T use_DCA
      REAL_T LSTEST(num_materials)
      REAL_T VOFTEST(num_materials)
      REAL_T LSTEST_EXTEND
      REAL_T LSMAX
      REAL_T mag,mag1,mag2,mag3
      REAL_T gx

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
      REAL_T temperature_cen(num_materials)
      REAL_T grad_tension(SDIM)
      REAL_T RR
      REAL_T dnrm(SDIM)
      REAL_T dxsten(SDIM)
      INTEGER_T im_primary_sten( &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance)

      INTEGER_T crossing_status
      INTEGER_T overall_crossing_status

      INTEGER_T cell_lo(3),cell_hi(3)

      REAL_T LS_CENTER(num_materials)
      REAL_T LS_OPP(num_materials)
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
      REAL_T, INTENT(out) :: ZEYU_thet_d
      REAL_T, INTENT(out) :: ZEYU_u_cl
      REAL_T ZEYU_u_slip
      REAL_T ZEYU_mu_l
      REAL_T ZEYU_mu_g
      REAL_T ZEYU_sigma
      REAL_T angle_im
      REAL_T dist_to_CL

      REAL_T maxcurv

      INTEGER_T im_liquid,im_vapor

      INTEGER_T local_index(3)
      REAL_T local_x(SDIM)
      REAL_T local_temperature(num_materials)

      INTEGER_T nhalf_height ! in: initheightLS

      nhalf_height=2*ngrow_distance+1 

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

       ! use_DCA=-1 static angle
       ! use_DCA=0  static angle
       ! if (probtype.eq.5501).and.(sdim.eq.3) then use_DCA=[-1,0,1,2]
       ! otherwise, if ZEYU_DCA_SELECT=-1 then use_DCA=-1 (static model)
       ! else if ZEYU_DCA_SELECT in [1,...,8], then
       !  use_DCA=ZEYU_DCA_SELECT+100.
      call get_use_DCA(use_DCA)
      ZEYU_thet_d=zero
      totaludotn=zero

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
      if (ngrow_distance.ne.4) then
       print *,"expecting ngrow_distance=4 in initheightLS"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"expecting ngrow_make_distance=3 in initheightLS"
       stop
      endif
      if (im.ge.im_opp) then
       print *,"im and im_opp invalid"
       stop
      endif
      if ((im.lt.1).or.(im.gt.num_materials)) then
       print *,"im invalid31"
       stop
      endif
      if ((im_opp.lt.1).or.(im_opp.gt.num_materials)) then
       print *,"im_opp invalid init height mof im_opp=",im_opp
       stop
      endif
      if ((iten.lt.1).or.(iten.gt.num_interfaces)) then
       print *,"iten invalid"
       stop
      endif
      call get_iten(im,im_opp,iten_test)
      if (iten.ne.iten_test) then
       print *,"iten and iten_test differ"
       stop
      endif
      if (vol_sten.gt.zero) then
       ! do nothing
      else
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
       klo_sten_ht=-ngrow_distance
       khi_sten_ht=ngrow_distance
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

      do im_sort=1,num_materials
       LS_CENTER(im_sort)=lssten(0,0,0,im_sort)
       LS_OPP(im_sort)=lssten(ii,jj,kk,im_sort)
      enddo

      if (LS_CENTER(im).ge.zero) then
       ! do nothing
      else if (LS_CENTER(im_opp).ge.zero) then
       ! do nothing
      else
       print *,"LS_CENTER is corrupt"
       do im_sort=1,num_materials
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
       stop
      endif 

      call get_LS_extend(LS_CENTER,iten,LS_CENTER_EXTEND)
      call get_LS_extend(LS_OPP,iten,LS_OPP_EXTEND)

      if (LS_CENTER_EXTEND*LS_OPP_EXTEND.gt.zero) then
       print *,"level set does not change sign"
       print *,"num_materials,iten= ",num_materials,iten
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

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if ((levelrz.eq.COORDSYS_RZ).or. &
               (levelrz.eq.COORDSYS_CYLINDRICAL)) then

       if (xcenter(1).gt.zero) then
        ! do nothing
       else
        print *,"expectimg xcenter(1)>0 for RZ or RT initheightLS"
        stop
       endif

      else
       print *,"initheightLS: levelrz invalid (a)"
       stop
      endif

      do im_sort=1,num_materials
       temperature_cen(im_sort)=mgoni_temp(0,0,0,im_sort)
      enddo

      if (static_flag.eq.0) then
       call get_user_tension(xcenter,time, &
         fort_tension,user_tension,temperature_cen)
      else if (static_flag.eq.1) then
       do dir2=1,num_interfaces
        user_tension(dir2)=fort_static_tension(iten)
       enddo
      else
       print *,"static_flag invalid"
       stop
      endif

      if (unscaled_min_curvature_radius.ge.two) then
       ! do nothing
      else
       print *,"unscaled_min_curvature_radius invalid"
       stop
      endif

      do dir2=1,num_interfaces

       if (fort_tension(dir2).ge.zero) then
        ! do nothing
       else
        print *,"fort_tension(dir2) invalid"
        stop
       endif

       if (user_tension(dir2).ge.zero) then
        ! do nothing
       else
        print *,"user_tension(dir2) invalid"
        stop
       endif

       if (fort_static_tension(dir2).ge.zero) then
        ! do nothing
       else
        print *,"fort_static_tension(dir2) invalid"
        stop
       endif

      enddo ! dir2=1,num_interfaces

      do dir2=1,SDIM
       mgoni_force(dir2)=zero
      enddo

      do i=-1,1
      do j=-1,1
      do k=klo_sten_short,khi_sten_short
       do imhold=1,num_materials
        LSTEST(imhold)=lssten(i,j,k,imhold)
       enddo
       call get_LS_extend(LSTEST,iten,LS1_save(D_DECL(i,j,k)))
      enddo
      enddo
      enddo ! i,j,k (-1 ... 1)

       ! centroid in absolute coordinate system
       ! returns a volume fraction
      call getvolume( &
       bfact,dx,xsten,nhalf_height, &
       LS1_save,volpos,facearea, &
       cenpos,VOFTOL,SDIM)

      if (facearea.ge.zero) then
       delta_mgoni=facearea/vol_sten
      else
       print *,"facearea invalid"
       stop
      endif

       ! declared in GLOBALUTIL.F90
      call get_LSNRM_extend(LS_CENTER,nrmcenter,iten,nfluid)
      RR=one
      call prepare_normal(nfluid,RR,mag)
      if (mag.gt.zero) then
       ! do nothing
      else
       print *,"nfluid mag became corrupt"
       stop
      endif

       ! Marangoni force:
       ! (I-nn^T)(grad sigma) delta=
       ! (grad sigma - (grad sigma dot n)n ) delta
      if (fort_tension_slope(iten).lt.zero) then

       if (static_flag.eq.0) then

        do iofs=-1,1
        do jofs=-1,1
        do kofs=-1,1
         local_index(1)=iofs
         local_index(2)=jofs
         local_index(3)=kofs
         do dir2=1,SDIM
          if (local_index(dir2).eq.-1) then
           local_x(dir2)=xsten(-2,dir2)
          else if (local_index(dir2).eq.1) then
           local_x(dir2)=xsten(2,dir2)
          else if (local_index(dir2).eq.0) then
           local_x(dir2)=xsten(0,dir2)
          else
           print *,"local_index invalid"
           stop
          endif
         enddo !dir2=1,SDIM

         do im_sort=1,num_materials
          local_temperature(im_sort)=mgoni_temp(iofs,jofs,kofs,im_sort)
         enddo

         call get_user_tension(local_x,time, &
          fort_tension,local_tension,local_temperature)
         mgoni_tension(iofs,jofs,kofs)=local_tension(iten)
        enddo
        enddo
        enddo

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
          if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
              (levelrz.eq.COORDSYS_RZ)) then
           RR=one
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           RR=xcenter(1)
           if (RR.gt.zero) then
            ! do nothing
           else
            print *,"RR invalid"
            stop
           endif
          else
           print *,"grad_tension: levelrz invalid"
           stop
          endif
         else if ((dir2.eq.3).and.(SDIM.eq.3)) then
          kofs=1
          RR=one
         else
          print *,"dir2 invalid"
          stop
         endif 

         grad_tension(dir2)=( &
          mgoni_tension(iofs,jofs,kofs)- &
          mgoni_tension(-iofs,-jofs,-kofs))/ &
          (RR*(xsten(2,dir2)-xsten(-2,dir2)))

         dotprod=dotprod+grad_tension(dir2)*nfluid(dir2)
        enddo ! dir2
        ! for the specific case,
        ! tension=sigma_0 + slope*(T-T0),
        ! then:
        ! grad sigma=slope * grad T
        !
        ! (I-nn^T)(grad sigma) delta
        do dir2=1,SDIM
         mgoni_force(dir2)= &
          (grad_tension(dir2)- &
           nfluid(dir2)*dotprod)*delta_mgoni
        enddo

       else if (static_flag.eq.1) then
        ! do nothing
       else
        print *,"static_flag invalid"
        stop
       endif

      else if (fort_tension_slope(iten).eq.zero) then
       ! do nothing
      else
       print *,"fort_tension_slope must be non-positive"
       stop
      endif 

      if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else if (levelrz.eq.COORDSYS_RZ) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       if (xcenter(1).le.zero) then
        print *,"xcenter invalid"
        stop
       endif
      else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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
      do i=-ngrow_distance,ngrow_distance
      do j=-ngrow_distance,ngrow_distance
      do k=klo_sten_ht,khi_sten_ht

       do imhold=1,num_materials
        LSTEST(imhold)=lssten(i,j,k,imhold)
        VOFTEST(imhold)=vofsten(i,j,k,imhold)
       enddo
        ! declared in GLOBALUTIL.F90
       call get_LS_extend(LSTEST,iten,lsdata(i,j,k))
       call get_VOF_extend(VOFTEST,iten,vofdata(i,j,k))

       call get_primary_material(LSTEST,imhold)
       im_primary_sten(i,j,k)=imhold

       if ((abs(i).le.1).and.(abs(j).le.1).and.(abs(k).le.1)) then
  
        if ((imhold.ge.1).and.(imhold.le.num_materials)) then
         if ((imhold.ne.im).and. &
             (imhold.ne.im_opp)) then
          if (im3.eq.0) then
           im3=imhold
           LSMAX=LSTEST(imhold)
          else if ((im3.ge.1).and.(im3.le.num_materials)) then
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
      enddo ! i,j,k (-ngrow_distance .. ngrow_distance)
     
      cos_angle=zero
      sin_angle=zero
      iten_13=0
      iten_23=0
      if ((im3.ge.1).and.(im3.le.num_materials)) then
        ! sigma_{i,j}cos(theta_{i,k})=sigma_{j,k}-sigma_{i,k}
        ! theta_{ik}=0 => material i wets material k.
        ! im is material "i"  ("fluid" material)
        ! im_opp is material "j"
       call get_CL_iten(im,im_opp,im3,iten_13,iten_23, &
        user_tension,cos_angle,sin_angle)

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

      if (num_materials.eq.1) then
       print *,"num_materials==1 not supported"
       stop
      else if (num_materials.eq.2) then
       if (im3.ne.0) then
        print *,"im3 invalid, num_materials=",num_materials
        stop
       endif
      else if (num_materials.gt.2) then
       ! do nothing
      else
       print *,"num_materials invalid"
       stop
      endif

      ! first: standard height function technique

! normal points towards "im"
! n=grad LS/|grad LS|

      imhold=im_primary_sten(0,0,0)
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

      lmin=-ngrow_distance
      lmax=ngrow_distance

      if ((levelrz.eq.COORDSYS_RZ).or.(levelrz.eq.COORDSYS_CYLINDRICAL)) then
       if (dircrit.eq.1) then ! horizontal column
        do while (xsten(2*lmin,dircrit).lt.zero)
         lmin=lmin+1
         if (2*lmin.gt.2*ngrow_distance+1) then
          print *,"lmin too big"
          stop
         endif
        enddo
       endif
      else if (levelrz.eq.COORDSYS_CARTESIAN) then
       ! do nothing
      else
       print *,"levelrz invalid just before call to get col ht ls"
       stop
      endif 
      xbottom=xsten(2*lmin-1,dircrit)
      xtop=xsten(2*lmax+1,dircrit)

      overall_crossing_status=1

      do iwidth=itanlo,itanhi
      do jwidth=jtanlo,jtanhi

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
          if (xsten(2*iwidth,1).le.zero) then
           iwidthnew=0
          endif
         endif
        endif
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        if (xsten(-2,1).le.zero) then
         print *,"xsten cannot be negative for levelrz==COORDSYS_CYLINDRICAL"
         stop
        endif
       else
        print *,"levelrz invalid initheight ls 4"
        stop
       endif
       
       do dir2=1,SDIM
        x_col(dir2)=xsten(0,dir2)
        x_col_avg(dir2)=half*(xsten(1,dir2)+xsten(-1,dir2))
        dx_col(dir2)=xsten(1,dir2)-xsten(-1,dir2)
       enddo
       x_col(itan)=xsten(2*iwidthnew,itan)
       dx_col(itan)=xsten(2*iwidthnew+1,itan)-xsten(2*iwidthnew-1,itan)
       x_col_avg(itan)=half*(xsten(2*iwidthnew+1,itan)+ &
                             xsten(2*iwidthnew-1,itan))

       if ((SDIM.eq.3).or. &
           ((SDIM.eq.2).and.(jtan.ge.1).and.(jtan.le.SDIM))) then
        x_col(jtan)=xsten(2*jwidth,jtan)
        dx_col(jtan)=xsten(2*jwidth+1,jtan)-xsten(2*jwidth-1,jtan)
        x_col_avg(jtan)=half*(xsten(2*jwidth+1,jtan)+ &
                              xsten(2*jwidth-1,jtan))
       else if ((SDIM.eq.2).and.(jtan.eq.3)) then
        !do nothing
       else
        print *,"sdim or jtan invalid"
        stop
       endif
 
       do kheight=-ngrow_distance,ngrow_distance

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

        if (SDIM.eq.2) then
         if (kcell.eq.0) then
          ! do nothing
         else
          print *,"expecting kcell=0"
          stop
         endif
        else if (SDIM.eq.3) then
         ! do nothing
        else
         print *,"SDIM invalid"
         stop
        endif 

        columnLS(kheight)=lsdata(icell,jcell,kcell)
        columnVOF(kheight)=vofdata(icell,jcell,kcell)
          
       enddo ! kheight

        ! declared in MOF.F90
       call get_col_ht_LS( &
        vof_height_function, &
        crossing_status, &
        bfact, &
        dx, &
        xsten, &
        dx_col, &
        x_col, &
        x_col_avg, &
        columnLS, &
        columnVOF, &
        col_ht_LS, &
        col_ht_VOF, &
        dircrit, & ! 1<=dircrit<=SDIM
        n1d, & ! n1d==1 => im material on top, n1d==-1 => im on bottom.
        SDIM)

       if (crossing_status.eq.1) then
        ! do nothing
       else if (crossing_status.eq.0) then
        if (1.eq.0) then
         print *,"no crossing found iwidth,jwidth= ",iwidth,jwidth
        endif
        overall_crossing_status=0 
       else
        print *,"crossing_status invalid"
        stop
       endif

       if ((col_ht_LS.ge.xbottom).and.(col_ht_LS.le.xtop)) then
        ! do nothing
       else
        print *,"col_ht_LS out of bounds"
        stop
       endif
       if ((col_ht_VOF.ge.xbottom).and.(col_ht_VOF.le.xtop)) then
        ! do nothing
       else
        print *,"col_ht_VOF out of bounds"
        stop
       endif

       htfunc_LS(iwidth,jwidth)=col_ht_LS
       htfunc_VOF(iwidth,jwidth)=col_ht_VOF
      enddo ! jwidth=-1,1
      enddo ! iwidth=-1,1

       ! analyze_heights is declared in: GLOBALUTIL.F90
      call analyze_heights( &
        htfunc_LS, &
        htfunc_VOF, &
        xsten, &
        nhalf_height, &
        itan,jtan, &
        curvHT_LS, &
        curvHT_VOF, &
        curvHT_choice, &
        dircrit, &
        xcenter, &
        n1d, &
        overall_crossing_status, &
        vof_height_function)

      ! above: use height function 
      ! below: use finite difference 

      if ((im3.lt.0).or.(im3.gt.num_materials).or. &
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
! if num_materials=4, 12 13 14 23 24 34

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
       else if ((im3.ge.1).and.(im3.le.num_materials)) then
        nfluid_def3(dir2)=nrmcenter(SDIM*(im3-1)+dir2)
       else
        print *,"im3 invalid"
        stop
       endif
     
      enddo ! dir2=1..sdim

      RR=one
      call prepare_normal(nfluid_cen,RR,mag) ! im or im_opp fluid
      call prepare_normal(nfluid_def1,RR,mag1) ! im fluid
      call prepare_normal(nfluid_def2,RR,mag2) ! im_opp fluid
      call prepare_normal(nfluid_def3,RR,mag3) ! im3 material

      if ((mag1.gt.zero).and. &
          (mag2.gt.zero).and. &
          (mag3.gt.zero).and. &
          (mag.gt.zero)) then
       ! do nothing
      else
       print *,"err:nfluid_def1, nfluid_def2, nfluid_def3, or nfluid_cen"
       print *,"mag1,mag2,mag3,mag=",mag1,mag2,mag3,mag
       print *,"nfluid_def1 associated with material im=",im
       print *,"nfluid_def2 associated with material im_opp=",im_opp
       print *,"nfluid_def3 associated with material im3=",im3
       stop
      endif

      do i=-1,1
      do j=-1,1
      do k=klo_sten_short,khi_sten_short 

       do im_sort=1,num_materials
        LSTEST(im_sort)=lssten(i,j,k,im_sort) 
       enddo

       ! nfluid is a normal to the im,im_opp interface
       ! and points towards the im material.

       do imhold=1,num_materials
        do dir2=1,SDIM
         nrmtest(SDIM*(imhold-1)+dir2)= &
           nrmsten(i,j,k,SDIM*(imhold-1)+dir2)
        enddo
       enddo ! imhold

       call get_LSNRM_extend(LSTEST,nrmtest,iten,nfluid)

       RR=one
       call prepare_normal(nfluid,RR,mag)
       if (mag.eq.zero) then
        do dir2=1,SDIM
         nfluid(dir2)=nfluid_cen(dir2)
        enddo
       else if (mag.gt.zero) then
        ! do nothing
       else
        print *,"mag invalid LEVELSET_3D.F90 1487"
        stop
       endif  

       do dir2=1,SDIM
        nfluid_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
       enddo

       if ((im3.ge.1).and.(im3.le.num_materials)) then

        if (is_rigid(im3).eq.1) then

         if ((im3.eq.im).or.(im3.eq.im_opp)) then
          print *,"im3 invalid" 
          stop
         endif

         ! nsolid points into the solid
         do dir2=1,SDIM
          nsolid(dir2)=nrmsten(i,j,k,SDIM*(im3-1)+dir2)
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

        else if (is_rigid(im3).eq.0) then
         do dir2=1,SDIM
          nsolid_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
         enddo
        else 
         print *,"is_rigid invalid LEVELSET_3D.F90"
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
        nfluid(dir2)=nrmsten(i,j,k,SDIM*(im-1)+dir2)
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
        print *,"mag invalid LEVELSET_3D.F90 1551"
        stop
       endif 
       do dir2=1,SDIM
        nmain_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
       enddo

       do dir2=1,SDIM
        nfluid(dir2)=nrmsten(i,j,k,SDIM*(im_opp-1)+dir2)
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
        print *,"mag invalid LEVELSET_3D.F90 1570"
        stop
       endif 
       do dir2=1,SDIM
        nopp_save(D_DECL(i,j,k),dir2)=nfluid(dir2)
       enddo

      enddo
      enddo
      enddo ! i,j,k=-1,..,1

      if ((im3.ge.1).and.(im3.le.num_materials)) then

       if (is_rigid(im3).eq.1) then
        
        if ((im3.eq.im).or.(im3.eq.im_opp)) then
         print *,"im3 invalid" 
         stop
        endif
        if (user_tension(iten).eq.zero) then  
         ! do nothing
        else if (user_tension(iten).gt.zero) then

         ! implement dynamic contact angle algorithm here.
         ! first project nfluid onto the solid (im3) material

         if ((use_DCA.eq.-1).or. & !static option; probtype.ne.5501
             (use_DCA.ge.0)) then  !static option; probtype.eq.5501

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
            LSTEST_EXTEND=lssten(i2,j2,k2,im3)
  
            if (LSTEST_EXTEND.lt.zero) then
             wt=one
            else if (LSTEST_EXTEND.ge.zero) then
             wt=CURVWT
            else
             print *,"LSTEST_EXTEND is NaN"
             stop
            endif

            totalwt=totalwt+wt

            udotn=zero
            do dir2=1,SDIM
             udotn=udotn+velsten(i2,j2,k2,dir2)*nproject(dir2)
            enddo
            totaludotn=totaludotn+wt*udotn
           enddo
           enddo
           enddo  ! i2,j2,k2=-1...1

           if (totalwt.gt.zero) then
            ! do nothing
           else
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

            ! use static angle
           if (use_DCA.eq.-1) then

            ZEYU_thet_d=acos(cos_angle)

            ! modify cos_angle (initialized above as static angle)
            ! THESE CASES if probtype=5501, from Yongsheng Lian.
            ! use_DCA=-1 static angle
            ! use_DCA=0 static angle
            ! use_DCA=1 Jiang
            ! use_DCA=2 Kistler
           else if ((use_DCA.eq.0).or. & !static
                    (use_DCA.eq.1).or. & !Jiang
                    (use_DCA.eq.2)) then !Kistler

              ! DCA_select_model is declared in GLOBALUTIL.F90
            call DCA_select_model(nproject,totaludotn,cos_angle, &
             liquid_viscosity,user_tension(iten),cos_angle,use_DCA)

            ZEYU_thet_d=acos(cos_angle)

           else if ((use_DCA.ge.101).and. & ! fort_ZEYU_DCA_SELECT>=1
                    (use_DCA.le.108)) then
            if (use_DCA.eq.101) then
             ! do nothing (GNBC model) since initheightLS handles only
             ! those models in which tangential wall velocity is input,
             ! and dynamic angle is output.
             ! FOR THE OPPOSITE CASE (angle is input, and wall velocity is
             ! output), one sets ns.law_of_the_wall=2.

             ! cases 2,...,8 for Zeyu's code.
             ! case 2 Jiang 1970 ...
             ! case 8 model=Cox 1986
            else if ((use_DCA.ge.102).and. &
                     (use_DCA.le.108)) then
             ZEYU_imodel=use_DCA-100

             call dynamic_contact_angle(ZEYU_mu_l, ZEYU_mu_g, ZEYU_sigma, &
               ZEYU_thet_s, &
               ZEYU_imodel, ZEYU_ifgnbc, ZEYU_lambda, &
               ZEYU_l_macro, ZEYU_l_micro, &
               ZEYU_dgrid, ZEYU_d_closest, ZEYU_thet_d_apparent, &
               ZEYU_u_cl, ZEYU_u_slip, ZEYU_thet_d)

             if (DEBUG_DYNAMIC_CONTACT_ANGLE.eq.1) then
              print *,"ZEYU_imodel= ",ZEYU_imodel
              print *,"non GNBC dynamic model"
              print *,"ZEYU_thet_s (rad,deg) ",ZEYU_thet_s, &
                      ZEYU_thet_s*180.0d0/Pi
              print *,"ZEYU_thet_d (rad,deg) ",ZEYU_thet_d, &
                      ZEYU_thet_d*180.0d0/Pi
              print *,"ZEYU_u_cl ",ZEYU_u_cl
              print *,"im3,iten ",im3,iten
             endif
             if (im.eq.im_liquid) then
              cos_angle=cos(ZEYU_thet_d)
             else if (im.eq.im_vapor) then
              cos_angle=Pi-cos(ZEYU_thet_d)
             else
              print *,"im invalid 2240: ",im
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

       else if (is_rigid(im3).eq.0) then
        ! do nothing
       else
        print *,"is_rigid invalid LEVELSET_3D.F90"
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

       if ((im3.ge.1).and.(im3.le.num_materials)) then

        if ((im3.eq.im).or.(im3.eq.im_opp)) then
         print *,"im3 invalid" 
         stop
        endif

        if (is_rigid(im3).eq.1) then

         ! cos(theta_1)=(sigma_23-sigma_13)/sigma_12
         ! cos(theta_2)=(-sigma_23+sigma_13)/sigma_12
         if (use_DCA.eq.101) then ! GNBC
          gamma1=half*user_tension(iten)
          gamma2=half*user_tension(iten)
         else if (use_DCA.ge.-1) then  ! all other cases.
          gamma1=half*(one-cos_angle)
          gamma2=half*(one+cos_angle)
         else
          print *,"use_DCA invalid"
          stop
         endif 

        else if (is_rigid(im3).eq.0) then

         gamma1=half*(user_tension(iten)-user_tension(iten_23)+ &
           user_tension(iten_13))/user_tension(iten)
         gamma2=half*(user_tension(iten)+user_tension(iten_23)- &
           user_tension(iten_13))/user_tension(iten)

        else
         print *,"is_rigid invalid LEVELSET_3D.F90"
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

       do imhold=1,num_materials
        LSTEST(imhold)=lssten(i,j,k,imhold)
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

       else if (is_rigid(im3).eq.0) then

        LSmain=LSTEST(im)
        LSopp=LSTEST(im_opp)
        do dir2=1,SDIM
         nmain(dir2)=nmain_save(D_DECL(i,j,k),dir2)
         nopp(dir2)=nopp_save(D_DECL(i,j,k),dir2)
        enddo

       else if (is_rigid(im3).eq.1) then

        call get_LS_extend(LSTEST,iten,LSmain)
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
        print *,"is_rigid(im3) invalid"
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

       imhold=im_primary_sten(i,j,k)

       do dir2=1,SDIM 
        n1=ncurv1_save(D_DECL(i,j,k),dir2)
        n2=ncurv2_save(D_DECL(i,j,k),dir2)
        if (im3.eq.0) then
         ngrid(dir2)=gamma1*n1-gamma2*n2
        else if (is_rigid(im3).eq.0) then
         ngrid(dir2)=gamma1*n1-gamma2*n2
        else if (is_rigid(im3).eq.1) then
         if (imhold.eq.im3) then
          ngrid(dir2)=n2
         else if ((imhold.ne.im3).and. &
                  (imhold.ge.1).and. &
                  (imhold.le.num_materials)) then 
          ngrid(dir2)=n1
         else
          print *,"imhold invalid"
          stop
         endif
        else
         print *,"is_rigid(im3) invalid"
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

        if (SDIM.eq.2) then
         if (k.eq.0) then
          ! do nothing
         else
          print *,"k invalid"
          stop
         endif
        else if (SDIM.eq.3) then
         ! do nothing
        else
         print *,"dimension bust"
         stop
        endif
        imhold=im_primary_sten(i,j,k)
        if (imhold.eq.im3) then
         im3_present_node=1
        else if ((imhold.ne.im3).and. &
                 (imhold.ge.1).and. &
                 (imhold.le.num_materials)) then 
         ! do nothing
        else
         print *,"imhold invalid"
         stop
        endif

       enddo
       enddo
       enddo ! i,j,k=cell_lo ... cell_hi

       if (wtnode.gt.zero) then
        ! do nothing
       else
        print *,"wtnode invalid"
        stop
       endif

       do dir2=1,SDIM

        dxsten(dir2)=xsten(2*cell_hi(dir2),dir2)- &
                     xsten(2*cell_lo(dir2),dir2)
        if (dxsten(dir2).gt.zero) then
         ! do nothing
        else
         print *,"dxsten invalid"
         stop
        endif
        RR=one
        if (dir2.eq.1) then
         ! do nothing
        else if (dir2.eq.2) then ! theta direction in cylindrical coord.
         if (levelrz.eq.COORDSYS_CYLINDRICAL) then 
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
        print *,"mag invalid LEVELSET_3D.F90 2108"
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
        print *,"mag invalid LEVELSET_3D.F90 2119"
        stop
       endif
       
       do dir2=1,SDIM
        n_node1(dir2)=n_node1LS(dir2)
       
        if (im3.eq.0) then
         n_node2(dir2)=n_node2LS(dir2)
        else if (is_rigid(im3).eq.0) then
         n_node2(dir2)=n_node2LS(dir2)
        else if (is_rigid(im3).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(im3) invalid"
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
        else if (is_rigid(im3).eq.0) then
         ngrid(dir2)=gamma1*n1-gamma2*n2
        else if (is_rigid(im3).eq.1) then
         if (im3_present_node.eq.1) then
          ngrid(dir2)=n2
         else if (im3_present_node.eq.0) then 
          ngrid(dir2)=n1
         else
          print *,"im3_present_node invalid"
          stop
         endif
        else
         print *,"is_rigid(im3) invalid"
         stop
        endif
       enddo ! dir2=1..sdim

       do dir2=1,SDIM
        if (dir2.eq.1) then
         if (levelrz.eq.COORDSYS_CARTESIAN) then
          RR=one
         else if ((levelrz.eq.COORDSYS_RZ).or. &
                  (levelrz.eq.COORDSYS_CYLINDRICAL)) then
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

      if (totalwt.gt.zero) then
       ! do nothing
      else
       print *,"totalwt invalid in initheightLS"
       stop
      endif

       ! dxsten=(xsten(2)+xsten(0))/2-(xsten(-2)+xsten(0))/2=
       !   (xsten(2)-xsten(-2))/2
      do dir2=1,SDIM
       dxsten(dir2)=xsten_curv(1,dir2)-xsten_curv(-1,dir2)
       if (dxsten(dir2).gt.zero) then
        ! do nothing
       else
        print *,"dxsten invalid"
        stop
       endif 
      enddo ! dir2

      do dir2=1,SDIM

        if (dir2.eq.1) then
         if (levelrz.eq.COORDSYS_CARTESIAN) then
          RR=one
         else if ((levelrz.eq.COORDSYS_RZ).or. &
                  (levelrz.eq.COORDSYS_CYLINDRICAL)) then
          RR=abs(xsten_curv(0,1))
         else
          print *,"levelrz invalid initheightLS: RR 3"
          stop
         endif
        else if (dir2.eq.2) then
         if (levelrz.eq.COORDSYS_CARTESIAN) then
          RR=one
         else if (levelrz.eq.COORDSYS_RZ) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          RR=one
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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

      if (unscaled_min_curvature_radius.ge.two) then
       maxcurv=one/(unscaled_min_curvature_radius*dxmax)
       if (levelrz.eq.COORDSYS_CARTESIAN) then
        if (SDIM.eq.2) then
         ! do nothing
        else if (SDIM.eq.3) then
         maxcurv=two*maxcurv
        else
         print *,"sdim invalid"
         stop
        endif     
       else if ((levelrz.eq.COORDSYS_RZ).or. &
                (levelrz.eq.COORDSYS_CYLINDRICAL)) then
        maxcurv=two*maxcurv
       else
        print *,"initheightLS: levelrz invalid (b)"
        stop
       endif

       if (curvFD.gt.maxcurv) then
        curvFD=maxcurv
       else if (curvFD.lt.-maxcurv) then
        curvFD=-maxcurv
       else if (abs(curvFD).le.maxcurv) then
        ! do nothing
       else
        print *,"curvFD is NaN"
        stop
       endif

       if (curvHT_choice.gt.maxcurv) then
        curvHT_choice=maxcurv
       else if (curvHT_choice.lt.-maxcurv) then
        curvHT_choice=-maxcurv
       else if (abs(curvHT_choice).le.maxcurv) then
        ! do nothing
       else
        print *,"curvHT_choice is NaN"
        stop
       endif

      else
       print *,"unscaled_min_curvature_radius invalid"
       stop
      endif

      if (1.eq.0) then
       print *,"xcenter ",xcenter(1),xcenter(2),xcenter(SDIM)
       print *,"dircrit,side,signside ",dircrit,side,signside
       print *,"im3,curvFD,curvHT_choice ",im3,curvFD,curvHT_choice
      endif

      return
      end subroutine initheightLS


       ! for finding areas internal to a cell, perturb each internal 
       ! interface, find areas and volumes, then check for the difference 
       ! in volumes divided by eps times the area.
       ! 
      subroutine fort_cellfaceinit( &
         tid, &
         tessellate, &  ! = 0,1, or 3
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
         ncellfrac) &
      bind(c,name='fort_cellfaceinit')

      use global_utility_module
      use global_distance_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: tid
      INTEGER_T, INTENT(in) :: tessellate  ! =0,1, or 3
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: ncellfrac
      INTEGER_T, INTENT(in) :: ngrow
      INTEGER_T, INTENT(in) :: DIMDEC(cface)
      INTEGER_T, INTENT(in) :: DIMDEC(maskfab)
      INTEGER_T, INTENT(in) :: DIMDEC(vofrecon)

      REAL_T, INTENT(out),target :: cface(DIMV(cface),ncellfrac)
      REAL_T, pointer :: cface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: maskfab(DIMV(maskfab),2)
      REAL_T, pointer :: maskfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: vofrecon(DIMV(vofrecon),num_materials*ngeom_recon)
      REAL_T, pointer :: vofrecon_ptr(D_DECL(:,:,:),:)

      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, INTENT(in) :: rz_flag
      REAL_T, INTENT(in) :: time

      INTEGER_T i,j,k
      INTEGER_T im
      INTEGER_T im_local
      INTEGER_T im_crit
      INTEGER_T im_opp
      INTEGER_T iface

      INTEGER_T ncellfrac_test
      INTEGER_T vofcomp
      INTEGER_T vofcomp2
      REAL_T vcenter(num_materials)
      REAL_T mofdata(num_materials*ngeom_recon)
      REAL_T mofdatavalid(num_materials*ngeom_recon)
      REAL_T local_facearea(num_materials,num_materials)
      REAL_T local_facearea_dimensional(num_materials,num_materials)
      REAL_T local_dist_to_line(num_materials,num_materials)
      REAL_T local_dist(num_materials,num_materials)
      REAL_T local_normal(num_materials,num_materials,SDIM)
      REAL_T local_facefrac(num_materials)
      INTEGER_T, parameter :: nhalf=3
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T dummy_tri(SDIM+1,SDIM)
      INTEGER_T nmax,ivert
      INTEGER_T dir2
      INTEGER_T shapeflag
      REAL_T multi_volume(num_materials)
      REAL_T multi_volume_offset(num_materials)
      REAL_T multi_cen(SDIM,num_materials)
      REAL_T multi_cen_offset(SDIM,num_materials)
      REAL_T multi_area(num_materials)
      REAL_T total_facearea
      REAL_T total_facearea_mat(num_materials)
      REAL_T uncaptured_volume_fraction
      REAL_T vfrac_fluid_sum
      REAL_T vfrac_solid_sum
      INTEGER_T loop_counter
      INTEGER_T num_processed_fluid
      INTEGER_T num_processed_solid
      INTEGER_T num_materials_fluid
      INTEGER_T num_materials_rigid
      INTEGER_T testflag
      REAL_T intercept
      REAL_T volcell
      REAL_T cencell(SDIM)
      REAL_T dist_tol,dxmax
      REAL_T dpair
      REAL_T areafrac
      INTEGER_T mask1,mask2
      INTEGER_T iten
      INTEGER_T is_processed(num_interfaces)
      INTEGER_T, parameter :: nhalf_box=1
      INTEGER_T cmofsten(D_DECL(-1:1,-1:1,-1:1))
      INTEGER_T local_tessellate
 
      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if ((tessellate.ne.0).and. &
          (tessellate.ne.1).and. &
          (tessellate.ne.3)) then
       print *,"tessellate invalid1"
       stop
      endif

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

      ! (num_materials,num_materials,3+sdim)
      ! im_inside,im_outside,3+sdim --> area, dist_to_line, dist, line normal.
      ncellfrac_test=num_materials*num_materials*(3+SDIM)
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

      cface_ptr=>cface
      maskfab_ptr=>maskfab
      vofrecon_ptr=>vofrecon
      call checkbound_array(fablo,fabhi,cface_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,maskfab_ptr,ngrow,-1)
      call checkbound_array(fablo,fabhi,vofrecon_ptr,ngrow,-1)
      
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
        
        do im=1,num_materials*ngeom_recon
         mofdata(im)=vofrecon(D_DECL(i,j,k),im)
        enddo

         ! vcenter = volume fraction 
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         vcenter(im)=mofdata(vofcomp)
        enddo ! im

        call check_full_cell_vfrac(vcenter, &
          tessellate, &  !=0,1, or 3
          im_crit)

        iface=0
        do im=1,num_materials
         do im_opp=1,num_materials
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

        if (iface.ne.num_materials*num_materials) then
         print *,"iface invalid"
         stop
        endif

        if ((im_crit.ge.1).and.(im_crit.le.num_materials)) then
         ! do nothing, there are no internal faces
        else if (im_crit.eq.0) then

         ! sum F_fluid=1  sum F_solid<=1
         local_tessellate=0
         call make_vfrac_sum_ok_copy( &
           cmofsten, &
           xsten,nhalf,nhalf_box, &
           bfact,dx, &
           local_tessellate, & ! =0 
           mofdata,mofdatavalid, &
           SDIM)

         shapeflag=0

          ! base case
          ! in: FORT_CELLFACEINIT
         call multi_get_volume_grid( &
          tessellate, &  ! =0,1, or 3
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
          SDIM, &
          shapeflag)

         do iten=1,num_interfaces
          is_processed(iten)=0
         enddo

         do im=1,num_materials
          total_facearea_mat(im)=zero
         enddo

         uncaptured_volume_fraction=one

         num_materials_rigid=0
         num_materials_fluid=0
         vfrac_fluid_sum=zero
         vfrac_solid_sum=zero

         do im=1,num_materials
          if (is_rigid(im).eq.0) then
           num_materials_fluid=num_materials_fluid+1
           vfrac_fluid_sum=vfrac_fluid_sum+vcenter(im)
          else if (is_rigid(im).eq.1) then
           num_materials_rigid=num_materials_rigid+1
           vfrac_solid_sum=vfrac_solid_sum+vcenter(im)
          else
           print *,"is_rigid invalid LEVELSET_3D.F90"
           stop
          endif
         enddo ! im=1..num_materials

         if (abs(one-vfrac_fluid_sum).gt.VOFTOL) then
          print *,"vfrac_fluid_sum invalid"
          stop
         endif
         if ((vfrac_solid_sum.gt.one+VOFTOL).or. &
             (vfrac_solid_sum.lt.zero)) then
          print *,"vfrac_solid_sum invalid"
          stop
         endif

         if (num_materials_fluid+num_materials_rigid.ne.num_materials) then
          print *,"num_materials_fluid or num_materials_rigid invalid"
          stop
         endif
         num_processed_solid=0
         num_processed_fluid=0

         if ((tessellate.eq.1).and.(vfrac_solid_sum.gt.zero)) then

          loop_counter=0
          do while ((loop_counter.lt.num_materials_rigid).and. &
                    (num_processed_solid.lt.num_materials_rigid).and. &
                    (uncaptured_volume_fraction.gt. &
                     one-vfrac_solid_sum))

            ! F,CEN,ORDER,SLOPE,INTERCEPT
           do im=1,num_materials
            if (is_rigid(im).eq.1) then
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
       
                ! solid case
                ! in: FORT_CELLFACEINIT
                ! no need to compute multi_area here.
                ! also, target volume is a cube, not a tet.
               call multi_get_volume_grid_simple( &
                tessellate, &  ! =1
                bfact,dx,xsten,nhalf, &
                mofdatavalid, &
                xsten,nhalf, &
                multi_volume_offset, &
                multi_cen_offset, &
                geom_xtetlist_uncapt(1,1,1,tid+1), &
                nmax, &
                nmax, &
                SDIM)

               mofdatavalid(vofcomp+2*SDIM+2)=intercept

               if (multi_volume_offset(im).gt.multi_volume(im)) then

                if (multi_area(im).gt.zero) then

                 do im_opp=1,num_materials
                  local_facefrac(im_opp)=zero
                 enddo
                 total_facearea=zero

                 do im_opp=1,num_materials
                  if (im_opp.ne.im) then
                   call get_iten(im,im_opp,iten)
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
                 enddo !im_opp=1..num_materials

                 if (total_facearea.gt.zero) then
                  do im_opp=1,num_materials
                   if (im_opp.ne.im) then
                    call get_iten(im,im_opp,iten)
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
                       cencell, &  ! cell centroid. 
                       dist_tol, &
                       SDIM, &
                       local_dist_to_line(im,im_opp))
 
                      call dist_centroid_line( &
                       xsten,nhalf, &
                       im,im_opp, & ! distance from line(im) to point(im_opp)
                       mofdata, &
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
                  enddo ! im_opp=1..num_materials

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
                print *,"num_materials_rigid=",num_materials_rigid
                print *,"num_materials_fluid=",num_materials_fluid
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
                do im_local=1,num_materials
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
                enddo  ! im=1..num_materials
                stop
               endif
  
              endif ! uncaptured_volume_fraction>0

             else if (testflag.eq.0) then
              ! do nothing
             else
              print *,"testflag invalid"
              stop
             endif 

            else if (is_rigid(im).eq.0) then
             ! do nothing
            else
             print *,"is_rigid invalid LEVELSET_3D.F90"
             stop
            endif

           enddo ! im=1..num_materials
           loop_counter=loop_counter+1
          enddo  ! while 
                 ! loop_counter<num_materials_fluid and
                 ! num_processed_fluid<num_materials_fluid and
                 ! uncaptured_volume_fraction>0

         else if ((tessellate.eq.0).or. &
                  (tessellate.eq.3).or. &
                  (vfrac_solid_sum.eq.zero)) then
          ! do nothing
         else
          print *,"tessellate or vfrac_solid_sum invalid"
          stop
         endif

         loop_counter=0
         do while ((loop_counter.lt.num_materials_fluid).and. &
                   (num_processed_fluid.lt.num_materials_fluid).and. &
                   (uncaptured_volume_fraction.gt.zero))

           ! F,CEN,ORDER,SLOPE,INTERCEPT
          do im=1,num_materials
           if (is_rigid(im).eq.0) then
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
               tessellate, &  !=0,1, or 3
               bfact,dx,xsten,nhalf, &
               mofdatavalid, &
               xsten,nhalf, &
               multi_volume_offset, &
               multi_cen_offset, &
               geom_xtetlist_uncapt(1,1,1,tid+1), &
               nmax, &
               nmax, &
               SDIM)

              mofdatavalid(vofcomp+2*SDIM+2)=intercept

              if (multi_volume(im).gt.zero) then

               if (multi_volume_offset(im).gt.multi_volume(im)) then

                if (multi_area(im).gt.zero) then

                 do im_opp=1,num_materials
                  local_facefrac(im_opp)=zero
                 enddo
                 total_facearea=zero

                 do im_opp=1,num_materials

                  if ((is_rigid(im_opp).eq.0).or. &
                      (tessellate.eq.1)) then

                   if (im_opp.ne.im) then
                    call get_iten(im,im_opp,iten)
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
                  else if ((is_rigid(im_opp).eq.1).and. &
                           ((tessellate.eq.0).or. &
                            (tessellate.eq.3))) then
                   ! do nothing
                  else
                   print *,"is_rigid or tessellate invalid"
                   stop
                  endif

                 enddo !im_opp=1..num_materials

                 if (total_facearea.gt.zero) then
                  do im_opp=1,num_materials
                   if ((is_rigid(im_opp).eq.0).or. &
                       (tessellate.eq.1)) then

                    if (im_opp.ne.im) then
                     call get_iten(im,im_opp,iten)
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
                        cencell, &  ! cell centroid. 
                        dist_tol, &
                        SDIM, &
                        local_dist_to_line(im,im_opp))
  
                       call dist_centroid_line( &
                        xsten,nhalf, &
                        im,im_opp, & ! distance from line(im) to point(im_opp)
                        mofdata, &
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
                   else if ((is_rigid(im_opp).eq.1).and. &
                            ((tessellate.eq.0).or. &
                             (tessellate.eq.3))) then
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
                 do vofcomp=1,num_materials*ngeom_recon
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
                do vofcomp=1,num_materials*ngeom_recon
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

           else if (is_rigid(im).eq.1) then
            ! do nothing
           else
            print *,"is_rigid invalid LEVELSET_3D.F90"
            stop
           endif

          enddo ! im=1..num_materials
          loop_counter=loop_counter+1
         enddo  ! while 
                ! loop_counter<num_materials_fluid and
                ! num_processed_fluid<num_materials_fluid and
                ! uncaptured_volume_fraction>0

        else
         print *,"im_crit out of range"
         stop
        endif

        do im=1,num_materials
         do im_opp=1,num_materials
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
        do im=1,num_materials
         do im_opp=1,num_materials
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
      end subroutine fort_cellfaceinit


      subroutine fort_curvstrip( &
       static_flag, &
       caller_string, &
       caller_string_len, &
       vof_height_function, &
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
       recon,DIMS(recon), &
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
       unscaled_min_curvature_radius, &
       num_curv, & ! num_interfaces * CURVCOMP_NCOMP
       ngrow_distance_in) &
      bind(c,name='fort_curvstrip')

      use ISO_C_BINDING, ONLY: C_CHAR,C_INT

      use global_utility_module
      use global_distance_module
      use probcommon_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: static_flag

      INTEGER_T, INTENT(in) :: nhistory

      CHARACTER(KIND=C_CHAR), INTENT(in) :: caller_string(*)
      INTEGER(C_INT), INTENT(in), VALUE :: caller_string_len
      CHARACTER(:), ALLOCATABLE :: fort_caller_string
      INTEGER_T :: fort_caller_string_len
      CHARACTER(len=255) :: pattern_string
      INTEGER_T :: pattern_string_len

      INTEGER_T, INTENT(in) :: vof_height_function
      INTEGER_T :: vof_height_function_local
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: ngrow_distance_in
      INTEGER_T, INTENT(in) :: num_curv ! num_interfaces * CURVCOMP_NCOMP
      INTEGER_T icurv
      REAL_T, INTENT(in) :: visc_coef
      REAL_T, INTENT(in) :: unscaled_min_curvature_radius

      INTEGER_T, INTENT(in) :: DIMDEC(history_dat)
      INTEGER_T, INTENT(in) :: DIMDEC(maskcov)
      INTEGER_T, INTENT(in) :: DIMDEC(masknbr)
      INTEGER_T, INTENT(in) :: DIMDEC(LSPC)
      INTEGER_T, INTENT(in) :: DIMDEC(recon)
      INTEGER_T, INTENT(in) :: DIMDEC(curvfab)
      INTEGER_T, INTENT(in) :: DIMDEC(velfab)
      INTEGER_T, INTENT(in) :: DIMDEC(denfab)

      INTEGER_T, INTENT(in) :: DIMDEC(vol)
      INTEGER_T, INTENT(in) :: DIMDEC(areax)
      INTEGER_T, INTENT(in) :: DIMDEC(areay)
      INTEGER_T, INTENT(in) :: DIMDEC(areaz)

      REAL_T, INTENT(out) :: curv_min
      REAL_T, INTENT(out) :: curv_max

      REAL_T, INTENT(in),target :: maskcov(DIMV(maskcov))
      REAL_T, pointer :: maskcov_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: vol(DIMV(vol))
      REAL_T, pointer :: vol_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: areax(DIMV(areax))
      REAL_T, INTENT(in), target :: areay(DIMV(areay))
      REAL_T, INTENT(in), target :: areaz(DIMV(areaz))
      REAL_T, pointer :: areax_ptr(D_DECL(:,:,:))
      REAL_T, pointer :: areay_ptr(D_DECL(:,:,:))
      REAL_T, pointer :: areaz_ptr(D_DECL(:,:,:))

      REAL_T, INTENT(out),target :: history_dat(DIMV(history_dat),nhistory)
      REAL_T, pointer :: history_dat_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: masknbr(DIMV(masknbr),4)
      REAL_T, pointer :: masknbr_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: LSPC(DIMV(LSPC),num_materials*(1+SDIM))
      REAL_T, pointer :: LSPC_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target ::  &
           recon(DIMV(recon),num_materials*ngeom_recon)
      REAL_T, pointer :: recon_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out),target :: curvfab(DIMV(curvfab),num_curv)
      REAL_T, pointer :: curvfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: velfab(DIMV(velfab),STATE_NCOMP_VEL)
      REAL_T, pointer :: velfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target ::  &
              denfab(DIMV(denfab),num_materials*num_state_material)
      REAL_T, pointer :: denfab_ptr(D_DECL(:,:,:),:)

      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T istenlo(3),istenhi(3)
      INTEGER_T LSstenlo(3),LSstenhi(3)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: bfact_grid
      REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, INTENT(in) :: rz_flag
      REAL_T, INTENT(in) :: time

      INTEGER_T i,j,k
      INTEGER_T iside,jside,kside
      INTEGER_T iface,jface,kface
      INTEGER_T i1,j1,k1
      INTEGER_T ii,jj,kk

      REAL_T LS(num_materials)
      REAL_T LS_fixed(num_materials)
      REAL_T LS_merge(num_materials)
      REAL_T LS_merge_fixed(num_materials)
      REAL_T LSSIDE(num_materials)
      REAL_T LSSIDE_fixed(num_materials)
      REAL_T LSSIDE_merge_fixed(num_materials)

      REAL_T xcenter(SDIM)
      INTEGER_T dirloc,dircrossing,dirstar
      INTEGER_T sidestar
      INTEGER_T at_RZ_axis

      INTEGER_T im,im_opp
      INTEGER_T im_opp_merge_test
      INTEGER_T im_opp_test
      INTEGER_T im_curv
      INTEGER_T vofcomp
      INTEGER_T im_majority
      INTEGER_T im_merge_majority
      INTEGER_T im_main,im_main_opp
      INTEGER_T iten
      INTEGER_T inormal

      REAL_T nrmPROBE(SDIM*num_materials)
      REAL_T nrmPROBE_merge(SDIM*num_materials)
      REAL_T LS_PROBE(num_materials)

      REAL_T nrmFD(SDIM*num_materials)
      REAL_T nrm_local(SDIM*num_materials)
      REAL_T nrm_local_merge(SDIM*num_materials)
      REAL_T nrm_mat(SDIM)
      REAL_T nrm_test(SDIM)
      REAL_T nrm_center(SDIM)
      REAL_T curv_cellHT
      REAL_T curv_cellFD
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

      REAL_T LSCEN_hold(num_materials)
      REAL_T LSCEN_hold_merge(num_materials)
      REAL_T LSCEN_hold_fixed(num_materials)

      REAL_T vof_hold(num_materials)
      REAL_T vof_hold_merge(num_materials)

      REAL_T LCEN,LSIDE
      REAL_T LSLEFT_EXTEND,LSRIGHT_EXTEND
      REAL_T LSLEFT_fixed(num_materials)
      REAL_T LSRIGHT_fixed(num_materials)
      REAL_T XLEFT,XRIGHT,XCEN

      REAL_T LS_STAR_merge_FIXED(-1:1,SDIM,num_materials)
      REAL_T LS_STAR_FIXED(-1:1,SDIM,num_materials)
      INTEGER_T im_star_merge_majority(-1:1,SDIM)
      INTEGER_T im_star_majority(-1:1,SDIM)

      REAL_T velsten(-1:1,-1:1,-1:1,SDIM)
      REAL_T nrmsten(-1:1,-1:1,-1:1,SDIM*num_materials)
      REAL_T mgoni_temp(-1:1,-1:1,-1:1,num_materials)
      REAL_T lssten( &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        num_materials)
      REAL_T vofsten( &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        -ngrow_distance:ngrow_distance, &
        num_materials)

      REAL_T, dimension(:,:), allocatable :: xsten0
      REAL_T, dimension(:,:), allocatable :: xsten_curv

      REAL_T x1dcen,x1dside,x1dcross
      REAL_T vol_sten
      REAL_T area_sten(SDIM,2)
      INTEGER_T side_index

      INTEGER_T, parameter :: nhalf=3
      INTEGER_T nhalf_height ! in: fort_curvstrip
      INTEGER_T mask1,mask2
      INTEGER_T local_mask
      INTEGER_T ihist
      REAL_T ZEYU_thet_d,ZEYU_u_cl
      INTEGER_T test_for_post_restart
   
      allocate(CHARACTER(caller_string_len) :: fort_caller_string)
      do i=1,caller_string_len
       fort_caller_string(i:i)=caller_string(i)
      enddo
      fort_caller_string_len=caller_string_len

      if (ngrow_distance.ne.4) then
       print *,"ngrow_distance invalid in curvstrip"
       stop
      endif 
      if (ngrow_distance_in.ne.4) then
       print *,"ngrow_distance_in invalid in curvstrip"
       stop
      endif 
 
      allocate(xsten0(-nhalf:nhalf,SDIM))
      allocate(xsten_curv( &
       -(2*ngrow_distance+1):(2*ngrow_distance+1), &
       SDIM))
 
      if (bfact.lt.1) then
       print *,"bfact invalid90"
       stop
      endif

      if ((level.lt.finest_level).and.(level.ge.0)) then
       if (bfact_grid.lt.4) then
        print *,"bfact_grid invalid in curvstrip(1)"
        stop
       endif
      else if (level.eq.finest_level) then
       if (bfact_grid.lt.2) then
        print *,"bfact_grid invalid in curvstrip(2)"
        stop
       endif
      else
       print *,"level invalid"
       stop
      endif

      if (unscaled_min_curvature_radius.ge.two) then
       ! do nothing
      else
       print *,"unscaled_min_curvature radius invalid"
       stop
      endif

      do i=1,num_interfaces
       if (fort_tension(i).ge.zero) then
        ! do nothing
       else
        print *,"fort_tension(i) invalid"
        stop
       endif
       if (fort_static_tension(i).ge.zero) then
        ! do nothing
       else
        print *,"fort_static_tension(i) invalid"
        stop
       endif
      enddo ! i=1,num_interfaces

      if (nhistory.eq.num_interfaces*2) then
       ! do nothing
      else
       print *,"nhistory invalid"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in curvstrip"
       stop
      endif

      pattern_string='post_restart'
      pattern_string_len=12
      test_for_post_restart=fort_pattern_test( &
        fort_caller_string,fort_caller_string_len, &
        pattern_string,pattern_string_len)

      if (test_for_post_restart.eq.0) then
       do dirloc=1,SDIM
        if ((fablo(dirloc)/bfact_grid)*bfact_grid.ne.fablo(dirloc)) then
         print *,"fablo mod bfact_grid not 0 in fort_curvstrip"
         stop
        endif
        if (((fabhi(dirloc)+1)/bfact_grid)*bfact_grid.ne.fabhi(dirloc)+1) then
         print *,"fabhi+1 mod bfact_grid not 0 in fort_curvstrip"
         stop
        endif
       enddo ! dirloc=1..sdim
      else if (test_for_post_restart.eq.1) then
       ! do nothing
      else
       print *,"fort_caller_string invalid"
       stop
      endif

      maskcov_ptr=>maskcov
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      LSPC_ptr=>LSPC
      call checkbound_array(fablo,fabhi,LSPC_ptr,ngrow_distance,-1)
      recon_ptr=>recon
      call checkbound_array(fablo,fabhi,recon_ptr,ngrow_distance,-1)
      masknbr_ptr=>masknbr
      call checkbound_array(fablo,fabhi,masknbr_ptr,1,-1)
      curvfab_ptr=>curvfab
      call checkbound_array(fablo,fabhi,curvfab_ptr,1,-1)
      velfab_ptr=>velfab
      call checkbound_array(fablo,fabhi,velfab_ptr,2,-1)
      denfab_ptr=>denfab
      call checkbound_array(fablo,fabhi,denfab_ptr,2,-1)

      vol_ptr=>vol
      call checkbound_array1(fablo,fabhi,vol_ptr,1,-1)
      areax_ptr=>areax
      areay_ptr=>areay
      areaz_ptr=>areaz
      call checkbound_array1(fablo,fabhi,areax_ptr,1,0)
      call checkbound_array1(fablo,fabhi,areay_ptr,1,1)
      call checkbound_array1(fablo,fabhi,areaz_ptr,1,SDIM-1)

      history_dat_ptr=>history_dat
      call checkbound_array(fablo,fabhi,history_dat_ptr,1,-1)

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
       print *,"rz_flag invalid in curvstrip"
       stop
      endif

      if (num_curv.ne.num_interfaces*CURVCOMP_NCOMP) then
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

          ! num_curv=num_interfaces * CURVCOMP_NCOMP
         do icurv=1,num_curv
          curvfab(D_DECL(i,j,k),icurv)=zero
         enddo 

         do im=1,num_materials
          LS(im)=LSPC(D_DECL(i,j,k),im)
         enddo

         call merge_levelset(xcenter,time,LS,LS_merge,static_flag)

          ! FIX_LS_tessellate is declared in: MOF.F90
          ! input : fluids tessellate, solids are embedded
          ! output: fluids tessellate and one and only one fluid LS is positive
         call FIX_LS_tessellate(LS_merge,LS_merge_fixed)
         call FIX_LS_tessellate(LS,LS_fixed)

         call get_primary_material(LS_fixed,im_majority)
         call get_primary_material(LS_merge_fixed,im_merge_majority)

         if ((is_rigid(im_merge_majority).eq.1).or. &
             (is_rigid(im_majority).eq.1).or. &
             (is_ice(im_merge_majority).eq.1).or. &
             (is_ice(im_majority).eq.1).or. &
             (is_FSI_rigid(im_merge_majority).eq.1).or. &
             (is_FSI_rigid(im_majority).eq.1)) then

          ! do nothing, all interface forces are 0
  
         else if ((xcenter(1).le.VOFTOL*dx(1)).and. &
                  (levelrz.eq.COORDSYS_RZ)) then
         
          ! do nothing, all interface forces are 0
         
         else if ((is_rigid(im_merge_majority).eq.0).and. &
                  (is_rigid(im_majority).eq.0).and. &
                  (is_ice(im_merge_majority).eq.0).and. &
                  (is_ice(im_majority).eq.0).and. &
                  (is_FSI_rigid(im_merge_majority).eq.0).and. &
                  (is_FSI_rigid(im_majority).eq.0)) then

          if (vol_sten.gt.zero) then
           ! do nothing
          else
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

            do im=1,num_materials
             LSSIDE(im)=LSPC(D_DECL(iside,jside,kside),im)
            enddo
            call merge_levelset(xcenter,time,LSSIDE,LS_merge,static_flag)
            call FIX_LS_tessellate(LS_merge,LSSIDE_merge_fixed)
            call get_primary_material(LSSIDE_merge_fixed,im_opp)
            do im=1,num_materials
             LS_STAR_merge_FIXED(sidestar,dirstar,im)=LSSIDE_merge_fixed(im)
            enddo
            im_star_merge_majority(sidestar,dirstar)=im_opp

            call FIX_LS_tessellate(LSSIDE,LSSIDE_fixed)
            call get_primary_material(LSSIDE_fixed,im_opp)
            do im=1,num_materials
             LS_STAR_FIXED(sidestar,dirstar,im)=LSSIDE_fixed(im)
            enddo
            im_star_majority(sidestar,dirstar)=im_opp

           enddo ! sidestar=-1,1,2

          enddo ! dirstar=1..sdim

           ! loop through all possible interfaces involving im_merge_majority
           ! and initialize curvfab
          do im_opp=1,num_materials

           donate_flag=0
  
           if (im_opp.eq.im_merge_majority) then
            ! do nothing
           else if (is_rigid(im_opp).eq.1) then
            ! do nothing
           else if (is_ice(im_opp).eq.1) then
            ! do nothing
           else if (is_FSI_rigid(im_opp).eq.1) then
            ! do nothing
           else if ((is_rigid(im_opp).eq.0).and. &
                    (is_ice(im_opp).eq.0).and. &
                    (is_FSI_rigid(im_opp).eq.0)) then

             ! im_main < im_main_opp
            if (im_merge_majority.lt.im_opp) then
             im_main=im_merge_majority
             im_main_opp=im_opp
            else if (im_merge_majority.gt.im_opp) then
             im_main=im_opp
             im_main_opp=im_merge_majority
            else
             print *,"im_merge_majority bust"
             stop
            endif
            call get_iten(im_main,im_main_opp,iten)

            do dirloc=1,SDIM
       
             do im=1,num_materials 
              LSLEFT_fixed(im)=LS_STAR_merge_FIXED(-1,dirloc,im)
              LSRIGHT_fixed(im)=LS_STAR_merge_FIXED(1,dirloc,im)
             enddo

             call get_LS_extend(LSLEFT_fixed,iten,LSLEFT_EXTEND)
             call get_LS_extend(LSRIGHT_fixed,iten,LSRIGHT_EXTEND)

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
  
              do im=1,num_materials
               LSSIDE_merge_fixed(im)=LS_STAR_merge_FIXED(sidestar,dirstar,im)
               LSSIDE_fixed(im)=LS_STAR_FIXED(sidestar,dirstar,im)
              enddo
              im_opp_merge_test=im_star_merge_majority(sidestar,dirstar)
              im_opp_test=im_star_majority(sidestar,dirstar)

              if ((im_opp_merge_test.eq.im_opp).and. &
                  (im_opp_test.eq.im_opp)) then

               at_RZ_axis=0
               if ((levelrz.eq.COORDSYS_RZ).and. &
                   (xsten0(-1,1).le.VOFTOL*dx(1)).and. &
                   (dirstar.eq.1).and. &
                   (sidestar.eq.-1)) then
                at_RZ_axis=1
               endif

               LCEN=LS_merge_fixed(im_merge_majority)
               LSIDE=LSSIDE_merge_fixed(im_opp)

               if ((LCEN*LSIDE.ge.zero).and. &
                   (abs(LCEN)+abs(LSIDE).gt.zero).and. &
                   (at_RZ_axis.eq.0)) then

                LCEN=-LS_merge_fixed(im_opp)
                LSIDE=-LSSIDE_merge_fixed(im_merge_majority)

                if ((LCEN*LSIDE.ge.zero).and. &
                    (abs(LCEN)+abs(LSIDE).gt.zero)) then

                 call get_crossing(x1dcross,x1dcen,x1dside,LCEN,LSIDE)

                 dxside=abs(x1dcross-x1dcen)

                  ! signcrossing points to im_merge_majority material
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

              else if ((im_opp_merge_test.ne.im_opp).or. &
                       (im_opp_test.ne.im_opp)) then
               ! do nothing
              else
               print *,"im_opp_merge_test of im_opp_test invalid"
               stop
              endif 

             enddo ! sidestar=-1,1,2
            enddo ! dirstar=1..sdim

            if (donate_flag.eq.1) then

              ! sidestar points away from im_merge_majority and towards the
              ! opposite material.
              ! signcrossing points towards im_merge_majority.
             sidestar=-signcrossing

             if (im_merge_majority.lt.im_opp) then
               ! signside points to im_merge_majority=im_main
              signside=signcrossing
              if ((im_merge_majority.eq.im_main).and. &
                  (im_opp.eq.im_main_opp)) then
               ! do nothing
              else
               print *,"im_merge_majority or im_opp invalid"
               stop
              endif
             else if (im_merge_majority.gt.im_opp) then
               ! signcrossing points towards im_merge_majority.
               ! signside points away from im_merge_majority
               ! signside points towards im_opp=im_main
              signside=-signcrossing
              if ((im_merge_majority.eq.im_main_opp).and. &
                  (im_opp.eq.im_main)) then
               ! do nothing
              else
               print *,"im_merge_majority or im_opp invalid"
               stop
              endif
             else
              print *,"im_merge_majority bust"
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

             vof_height_function_local=vof_height_function

             if (mask2.eq.0) then ! mask2==0 => not interior cell 
              vof_height_function_local=0
             else if (mask2.eq.1) then ! mask2==1 => interior cell
              ! do nothing
             else
              print *,"mask2 invalid"
              stop
             endif

             LSstenlo(3)=0
             LSstenhi(3)=0
             do dirloc=1,SDIM
              LSstenlo(dirloc)=-ngrow_distance
              LSstenhi(dirloc)=ngrow_distance
             enddo
   
             nhalf_height=2*ngrow_distance+1 
             call gridsten_level(xsten_curv,i,j,k,level,nhalf_height)

             ! get normals at the cell center.
             do inormal=1,SDIM*num_materials
              nrmPROBE(inormal)=LSPC(D_DECL(i,j,k),num_materials+inormal)
             enddo ! inormal
             do im_curv=1,num_materials
              LS_PROBE(im_curv)=LSPC(D_DECL(i,j,k),im_curv)
             enddo
             call merge_normal(xcenter,time, &
               LS_PROBE, &
               nrmPROBE,nrmPROBE_merge,static_flag)

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

              do im_curv=1,num_materials
               LSRIGHT_EXTEND=LS_STAR_merge_FIXED(1,dirstar,im_curv)
               LSLEFT_EXTEND=LS_STAR_merge_FIXED(-1,dirstar,im_curv)
               LCEN=LS_merge_fixed(im_curv)

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
              enddo ! im_curv=1..num_materials

             enddo ! dirstar=1..sdim

             do im_curv=1,num_materials
              do dirloc=1,SDIM
               inormal=(im_curv-1)*SDIM+dirloc
               nrm_mat(dirloc)=nrmPROBE_merge(inormal)
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
               if (nrm_test(dircrossing)*critsign.gt.zero) then
                ! do nothing
               else
                print *,"critsign and nrm_test mismatch"
                print *,"dircrossing = ",dircrossing
                print *,"sidestar= ",sidestar
                print *,"signcrossing= ",signcrossing
                print *,"dxcrossing= ",dxcrossing
                print *,"critsign=",critsign
                print *,"im_curv=",im_curv
                print *,"im_merge_majority= ",im_merge_majority
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
                print *,"nrmPROBE_merge points towards im_curv= ",im_curv
                do dirloc=1,SDIM
                 inormal=(im_curv-1)*SDIM+dirloc
                 print *,"dirloc,nrmPROBE_merge ", &
                    dirloc,nrmPROBE_merge(inormal)
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
               if ((nrm_mat(dircrossing)*nrm_test(dircrossing).gt.zero).and. &
                   (abs(nrm_mat(dircrossing)).gt. &
                    half*abs(nrm_test(dircrossing)))) then
                ! do nothing
               else if  &
                  ((nrm_mat(dircrossing)*nrm_test(dircrossing).le.zero).or. &
                   (abs(nrm_mat(dircrossing)).le. &
                    half*abs(nrm_test(dircrossing)))) then
                do dirloc=1,SDIM
                 nrm_mat(dirloc)=nrm_test(dirloc)
                enddo
               else
                print *,"nrm_mat or nrm_test is NaN"
                stop
               endif
  
              else if ((im_curv.ne.im_main).and. &
                       (im_curv.ne.im_main_opp)) then
               ! do nothing
              else
               print *,"im_curv or im_main bust"
               stop
              endif 

               ! if R-Theta, then N(2) -> N(2)/RR + renormalize.
              RR=xcenter(1) 
               ! declared in GLOBALUTIL.F90
              call prepare_normal(nrm_mat,RR,mag)

              if (mag.eq.zero) then
               if (is_rigid(im_curv).eq.0) then
                ! do nothing
               else if (is_rigid(im_curv).eq.1) then
                 ! nrm_test obtained via finite differences.
                do dirloc=1,SDIM
                 nrm_mat(dirloc)=nrm_test(dirloc)
                enddo
                 ! if R-Theta, then N(2) -> N(2)/RR + renormalize.
                RR=xcenter(1) 
                 ! declared in GLOBALUTIL.F90
                call prepare_normal(nrm_mat,RR,mag)
                if (mag.eq.zero) then
                 ! do nothing
                else if (mag.gt.zero) then
                 ! do nothing
                else
                 print *,"mag is invalid for im_curv material: ",mag
                 stop
                endif
               else
                print *,"is_rigid invalid LEVELSET_3D.F90"
                stop
               endif

              else if (mag.gt.zero) then
               ! do nothing
              else
               print *,"mag invalid; mag=",mag
               stop
              endif

              do dirloc=1,SDIM
               inormal=(im_curv-1)*SDIM+dirloc
               nrmPROBE_merge(inormal)=nrm_mat(dirloc)
              enddo

             enddo ! im_curv=1..num_materials

             if (1.eq.0) then
              print *,"xcenter ",xcenter(1),xcenter(2),xcenter(SDIM)
              print *,"dircrossing ",dircrossing
              print *,"im_merge_majority,im_opp,im_main,im_main_opp ", &
               im_merge_majority,im_opp,im_main,im_main_opp
             endif

             ! i1,j1,k1=-ngrow_distance ... ngrow_distance
             do i1=LSstenlo(1),LSstenhi(1)
             do j1=LSstenlo(2),LSstenhi(2)
             do k1=LSstenlo(3),LSstenhi(3)

              do inormal=1,SDIM*num_materials
               call safe_data(i+i1,j+j1,k+k1,num_materials+inormal, &
                 LSPC_ptr,nrm_local(inormal))
              enddo

              do im_curv=1,num_materials
               call safe_data(i+i1,j+j1,k+k1,im_curv, &
                 LSPC_ptr,LSCEN_hold(im_curv))
               vofcomp=(im_curv-1)*ngeom_recon+1
               call safe_data(i+i1,j+j1,k+k1,vofcomp, &
                 recon_ptr,vof_hold(im_curv))
              enddo !im_curv=1..num_materials

              call merge_normal(xcenter,time, &
                 LSCEN_hold, &
                 nrm_local,nrm_local_merge,static_flag)

              call merge_levelset(xcenter,time,LSCEN_hold, &
                 LSCEN_hold_merge,static_flag)
              call merge_vof(xcenter,time,vof_hold,vof_hold_merge, &
                 static_flag)

              do im_curv=1,num_materials
               vofsten(i1,j1,k1,im_curv)=vof_hold_merge(im_curv)
              enddo

              call FIX_LS_tessellate(LSCEN_hold_merge,LSCEN_hold_fixed)
 
              do im_curv=1,num_materials

               lssten(i1,j1,k1,im_curv)=LSCEN_hold_fixed(im_curv)

               if ((abs(i1).le.1).and.(abs(j1).le.1).and.(abs(k1).le.1)) then
                do dirloc=1,SDIM
                 inormal=(im_curv-1)*SDIM+dirloc
                 nrm_mat(dirloc)=nrm_local_merge(inormal)
                enddo
                RR=one
                if (levelrz.eq.COORDSYS_CARTESIAN) then
                 ! do nothing
                else if (levelrz.eq.COORDSYS_RZ) then
                 if (SDIM.ne.2) then
                  print *,"levelrz invalid"
                  stop
                 endif
                else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
                 RR=xsten_curv(2*i1,1)
                else
                 print *,"transformed normal: levelrz invalid"
                 stop
                endif
                call prepare_normal(nrm_mat,RR,mag)
                do dirloc=1,SDIM
                 inormal=(im_curv-1)*SDIM+dirloc
                 nrmsten(i1,j1,k1,inormal)=nrm_mat(dirloc)
                enddo
               else if ((abs(i1).le.ngrow_distance).and. &
                        (abs(j1).le.ngrow_distance).and. &
                        (abs(k1).le.ngrow_distance)) then
                ! do nothing
               else
                print *,"i1,j1, or k1 invalid"
                stop
               endif

              enddo ! im_curv=1..num_materials

             enddo
             enddo
             enddo ! i1,j1,k1=LSstenlo,LSstenhi(init nrmsten,lssten,vofsten)
 
             ! i1,j1,k1=-1..1
             do i1=istenlo(1),istenhi(1)
             do j1=istenlo(2),istenhi(2)
             do k1=istenlo(3),istenhi(3)

              do dirloc=1,SDIM
               velsten(i1,j1,k1,dirloc)= &
                velfab(D_DECL(i+i1,j+j1,k+k1),dirloc)
              enddo

              do im_curv=1,num_materials
               itemperature= &
                 (im_curv-1)*num_state_material+ENUM_TEMPERATUREVAR+1
               mgoni_temp(i1,j1,k1,im_curv)= &
                denfab(D_DECL(i+i1,j+j1,k+k1),itemperature)
              enddo
  
             enddo
             enddo
             enddo  ! i1,j1,k1 (init velsten)

             ! tension used to find contact angle (scaling not 
             ! necessary)
             call initheightLS( &
              static_flag, &
              vof_height_function_local, &
              i,j,k, &
              level, &
              finest_level, &
              bfact,dx, &
              xcenter, &
              !num_materials x sdim components("nrmcenter" in initheightLS)
              nrmPROBE_merge, &
              dircrossing, &
              sidestar, &
              signside, &
              time, &
              xsten_curv, &
              velsten, &
              mgoni_temp, &
              lssten, &
              vofsten, &
              !3x3x3x num_materials x sdim components
              nrmsten, &
              vol_sten, &
              area_sten, &
              curv_cellHT, & !intent(out)
              curv_cellFD, & !intent(out)
              mgoni_force, & !intent(out) (I-nn^T)(grad sigma) delta
              ZEYU_thet_d, & !intent(out)
              ZEYU_u_cl, & !intent(out)
              im3, & !intent(out)
              visc_coef, &
              unscaled_min_curvature_radius, &
              im_main, & !intent(in)
              im_main_opp, & !intent(in) 
              iten) !intent(in)

             if (im3.eq.0) then
              ! do nothing
             else if ((im3.ge.1).and.(im3.le.num_materials)) then

               ! triple point algorithm: 3 fluids
               ! contact line algorithm: 2 fluids and "is_rigid" material.
              if ((is_ice(im3).eq.1).or. &
                  (is_FSI_rigid(im3).eq.1)) then
               curv_cellHT=zero
               curv_cellFD=zero
              else if ((is_ice(im3).eq.0).and. &
                       (is_FSI_rigid(im3).eq.0)) then
               !do nothing
              else
               print *,"is_ice or is_FSI_rigid invalid"
               stop
              endif

             else
              print *,"im3 invalid: ",im3
              stop
             endif

             if (DEBUG_CURVATURE.eq.1) then
              print *,"i,j,k,dircrossing,sidestar,nrm ", &
               i,j,k,dircrossing,sidestar, &
               nrm_center(1),nrm_center(2),nrm_center(SDIM)
              print *,"curv_cellHT,curv_cellFD ",curv_cellHT,curv_cellFD
             endif
          
             ihist=(iten-1)*2
             history_dat(D_DECL(i,j,k),ihist+1)=ZEYU_thet_d
             history_dat(D_DECL(i,j,k),ihist+2)=ZEYU_u_cl

             icurv=(iten-1)*CURVCOMP_NCOMP
             curvfab(D_DECL(i,j,k),icurv+CURVCOMP_HTFUNC_CURV+1)=curv_cellHT
             curvfab(D_DECL(i,j,k),icurv+CURVCOMP_FD_CURV+1)=curv_cellFD
             do dirloc=1,SDIM
              curvfab(D_DECL(i,j,k),icurv+CURVCOMP_MARANGONI+dirloc)= &
                 mgoni_force(dirloc)
             enddo
              ! dir=1..sdim
              ! side=-1 or 1
             curvfab(D_DECL(i,j,k),icurv+CURVCOMP_DIRSIDE_FLAG+1)= &
              dircrossing*sidestar
             curvfab(D_DECL(i,j,k),icurv+CURVCOMP_MATERIAL3_ID+1)=im3

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

          enddo ! im_opp=1...num_materials

         else
          print *,"im_majority or im_merge_majority invalid"
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
      deallocate(fort_caller_string)

      return
      end subroutine fort_curvstrip

      subroutine fort_gettypefab( &
       source_fab, &
       DIMS(source_fab), &
       typefab, &
       DIMS(typefab), &
       xlo,dx, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       type_flag, &
       ncomp_type, &
       ncomp_source, &
       zero_diag_flag) &
      bind(c,name='fort_gettypefab')

      use probf90_module
      use global_utility_module
 

      IMPLICIT NONE

      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: xlo(SDIM)
      INTEGER_T, INTENT(in) :: ncomp_type
      INTEGER_T, INTENT(in) :: ncomp_source
      INTEGER_T, INTENT(in) :: zero_diag_flag
      INTEGER_T, INTENT(out) :: type_flag(ncomp_type)

      INTEGER_T, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, INTENT(in) :: bfact

      INTEGER_T, INTENT(in) ::  DIMDEC(source_fab)
      INTEGER_T, INTENT(in) ::  DIMDEC(typefab)

      REAL_T, INTENT(in),target :: source_fab(DIMV(source_fab),ncomp_source)
      REAL_T, pointer :: source_fab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out),target :: typefab(DIMV(typefab))
      REAL_T, pointer :: typefab_ptr(D_DECL(:,:,:))

      INTEGER_T i,j,k,im,base_type


      if (bfact.lt.1) then
       print *,"bfact invalid91"
       stop
      endif
      if (zero_diag_flag.eq.0) then
       if (ncomp_source.eq.num_materials*(1+SDIM)) then
        ! do nothing
       else
        print *,"ncomp_source invalid"
        stop
       endif
       if (ncomp_type.eq.num_materials) then
        ! do nothing
       else
        print *,"ncomp_type invalid"
        stop
       endif
      else if (zero_diag_flag.eq.1) then
       if (ncomp_source.eq.1) then
        ! do nothing
       else
        print *,"ncomp_source invalid"
        stop
       endif
       if (ncomp_type.eq.2) then
        ! do nothing
       else
        print *,"ncomp_type invalid"
        stop
       endif
      else
       print *,"zero_diag_flag invalid"
       stop
      endif

      source_fab_ptr=>source_fab
      typefab_ptr=>typefab
      call checkbound_array(fablo,fabhi,source_fab_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,typefab_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       if (zero_diag_flag.eq.0) then

        base_type=1
        do im=2,num_materials
         if (source_fab(D_DECL(i,j,k),im).gt. &
             source_fab(D_DECL(i,j,k),base_type)) then
          base_type=im
         endif
        enddo

        do im=1,num_materials
         if (is_rigid(im).eq.1) then
          if (source_fab(D_DECL(i,j,k),im).ge.zero) then
           base_type=im
          endif
         else if (is_rigid(im).eq.0) then
          ! do nothing
         else
          print *,"is_rigid invalid LEVELSET_3D.F90"
          stop
         endif
        enddo ! im=1..num_materials

        typefab(D_DECL(i,j,k))=base_type
        if ((base_type.gt.num_materials).or.(base_type.lt.1)) then
         print *,"base_type invalid"
         stop
        else
         type_flag(base_type)=1
        endif

       else if (zero_diag_flag.eq.1) then

        if (source_fab(D_DECL(i,j,k),1).eq.zero) then !diagonal=0
         base_type=1
        else if (source_fab(D_DECL(i,j,k),1).eq.one) then !diagonal>0
         base_type=2
        else
         print *,"source_fab(D_DECL(i,j,k),1) invalid"
         stop
        endif
        typefab(D_DECL(i,j,k))=base_type
        type_flag(base_type)=1

       else
        print *,"zero_diag_flag invalid"
        stop
       endif

      enddo
      enddo
      enddo

      return
      end subroutine fort_gettypefab

      subroutine fort_getcolorsum( &
       tid_current, &
       operation_flag, &
       sweep_num, &
       tessellate, &
       distribute_mdot_evenly, &
       constant_volume_mdot, &
       distribute_from_target, &
       constant_density_all_time, & ! 1..num_materials
       cur_time_slab, &
       dt, &
       dx, &
       xlo, &
       nstate, &
       snew,DIMS(snew), &
       mdot, &
       DIMS(mdot), &
       mdot_complement, &
       DIMS(mdot_complement), &
       LS,DIMS(LS), &
       VEL,DIMS(VEL), &
       DEN,DIMS(DEN), &
       VOF,DIMS(VOF), &
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
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       rzflag, &
       num_colors, &
       cum_blobdata, &
       cum_mdot_data, &
       cum_mdot_complement_data, &
       level_blobdata, &
       level_blobtypedata, &
       level_mdot_data, &
       level_mdot_complement_data, &
       level_mdot_data_redistribute, &
       level_mdot_complement_data_redistribute, &
       arraysize, &
       mdot_arraysize, &
       ncomp_mdot_alloc, &
       ncomp_mdot, &
       levelbc, &
       velbc, &
       material_type_lowmach, &
       material_type_visual, &
       nface_dst, &
       ncellfrac) &
      bind(c,name='fort_getcolorsum')

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: tid_current
      INTEGER_T, INTENT(in) :: operation_flag
      INTEGER_T, INTENT(in) :: nstate
      INTEGER_T, INTENT(in) :: sweep_num
      INTEGER_T, INTENT(in) :: tessellate
      INTEGER_T, INTENT(in) :: nface_dst,ncellfrac
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      REAL_T, INTENT(in) :: cur_time_slab
      REAL_T, INTENT(in) :: dt
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: xlo(SDIM)
      INTEGER_T, INTENT(in) :: levelbc(SDIM,2)
      INTEGER_T, INTENT(in) :: velbc(SDIM,2,SDIM)
      INTEGER_T, INTENT(in) :: ncomp_mdot_alloc
      INTEGER_T, INTENT(in) :: ncomp_mdot
      INTEGER_T, INTENT(in) :: material_type_lowmach(num_materials)
      INTEGER_T, INTENT(in) :: material_type_visual(num_materials)
      INTEGER_T, INTENT(in) :: distribute_mdot_evenly(2*num_interfaces)
      INTEGER_T, INTENT(in) :: constant_volume_mdot(2*num_interfaces)
      INTEGER_T, INTENT(in) :: distribute_from_target(2*num_interfaces)
      INTEGER_T, INTENT(in) :: constant_density_all_time(num_materials)

      INTEGER_T :: i,j,k
      INTEGER_T :: ii,jj,kk
      INTEGER_T :: iface,jface,kface
      INTEGER_T :: face_index
 
      INTEGER_T, INTENT(in) :: rzflag
      INTEGER_T, INTENT(in) :: num_colors
      INTEGER_T, INTENT(in) :: arraysize
      INTEGER_T, INTENT(in) :: mdot_arraysize
      INTEGER_T, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: DIMDEC(snew)
      INTEGER_T, INTENT(in) :: DIMDEC(mdot)
      INTEGER_T, INTENT(in) :: DIMDEC(mdot_complement)
      INTEGER_T, INTENT(in) :: DIMDEC(LS)
      INTEGER_T, INTENT(in) :: DIMDEC(VEL)
      INTEGER_T, INTENT(in) :: DIMDEC(DEN)
      INTEGER_T, INTENT(in) :: DIMDEC(VOF)
      INTEGER_T, INTENT(in) :: DIMDEC(xface)
      INTEGER_T, INTENT(in) :: DIMDEC(yface)
      INTEGER_T, INTENT(in) :: DIMDEC(zface)
      INTEGER_T, INTENT(in) :: DIMDEC(areax)
      INTEGER_T, INTENT(in) :: DIMDEC(areay)
      INTEGER_T, INTENT(in) :: DIMDEC(areaz)
      INTEGER_T, INTENT(in) :: DIMDEC(cellfab)
      INTEGER_T, INTENT(in) :: DIMDEC(typefab)
      INTEGER_T, INTENT(in) :: DIMDEC(color)
      INTEGER_T, INTENT(in) :: DIMDEC(mask)
      REAL_T, INTENT(inout) ::level_blobdata(arraysize)
      REAL_T, INTENT(inout) ::level_mdot_data(mdot_arraysize)
      REAL_T, INTENT(inout) ::level_mdot_complement_data(mdot_arraysize)
      REAL_T, INTENT(inout) ::level_mdot_data_redistribute(mdot_arraysize)
      REAL_T, INTENT(inout) :: &
              level_mdot_complement_data_redistribute(mdot_arraysize)
      REAL_T, INTENT(in) :: cum_blobdata(arraysize)
      REAL_T, INTENT(in) :: cum_mdot_data(mdot_arraysize)
      REAL_T, INTENT(in) :: cum_mdot_complement_data(mdot_arraysize)
      INTEGER_T, INTENT(inout) :: level_blobtypedata(num_colors)

      REAL_T, INTENT(inout), target :: snew(DIMV(snew),nstate)
      REAL_T, pointer :: snew_ptr(D_DECL(:,:,:),:)

      REAL_T, INTENT(inout), target :: mdot(DIMV(mdot),ncomp_mdot)
      REAL_T, pointer :: mdot_ptr(D_DECL(:,:,:),:)

      REAL_T, INTENT(inout), target :: &
              mdot_complement(DIMV(mdot_complement),ncomp_mdot)
      REAL_T, pointer :: mdot_complement_ptr(D_DECL(:,:,:),:)

      REAL_T, INTENT(in), target :: typefab(DIMV(typefab))
      REAL_T, pointer :: typefab_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: LS(DIMV(LS),num_materials*(1+SDIM))
      REAL_T, pointer :: LS_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: &
           VEL(DIMV(VEL),STATE_NCOMP_VEL+STATE_NCOMP_PRES)
      REAL_T, pointer :: VEL_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: DEN(DIMV(DEN),num_materials*num_state_material)
      REAL_T, pointer :: DEN_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: VOF(DIMV(VOF),num_materials*ngeom_recon)
      REAL_T, pointer :: VOF_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: xface(DIMV(xface),nface_dst)
      REAL_T, INTENT(in), target :: yface(DIMV(yface),nface_dst)
      REAL_T, INTENT(in), target :: zface(DIMV(zface),nface_dst)
      REAL_T, pointer :: xface_ptr(D_DECL(:,:,:),:)
      REAL_T, pointer :: yface_ptr(D_DECL(:,:,:),:)
      REAL_T, pointer :: zface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: areax(DIMV(areax))
      REAL_T, INTENT(in), target :: areay(DIMV(areay))
      REAL_T, INTENT(in), target :: areaz(DIMV(areaz))
      REAL_T, pointer :: areax_ptr(D_DECL(:,:,:))
      REAL_T, pointer :: areay_ptr(D_DECL(:,:,:))
      REAL_T, pointer :: areaz_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: cellfab(DIMV(cellfab),ncellfrac)
      REAL_T, pointer :: cellfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: color(DIMV(color))
      REAL_T, pointer :: color_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))

      INTEGER_T dir,side
      INTEGER_T dir2
      INTEGER_T vofcomp
      INTEGER_T dencomp
      INTEGER_T base_type
      INTEGER_T opposite_color(num_materials)
      INTEGER_T typeside
      INTEGER_T colorside
      INTEGER_T ic
      INTEGER_T ic_center
      INTEGER_T ic_base
      INTEGER_T ic_base_mdot
      INTEGER_T icolor
      REAL_T local_facearea(num_materials,num_materials)
      REAL_T local_dist_to_line(num_materials,num_materials)
      REAL_T local_normal(num_materials,num_materials,SDIM)
      REAL_T local_dist(num_materials,num_materials)
      REAL_T vfrac
      REAL_T vol
      REAL_T mass
      REAL_T den_mat
      REAL_T TEMP_mat
      REAL_T internal_energy
      REAL_T massfrac_parm(num_species_var+1)
      REAL_T pressure_local
      REAL_T xsten_stencil(-3:3,SDIM)
      REAL_T xstencil_point(SDIM)
      INTEGER_T, parameter :: nhalf=3
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T dx_sten(SDIM)
      REAL_T cencell(SDIM)
      INTEGER_T i1,j1,k1
      INTEGER_T k1lo,k1hi
      INTEGER_T im,im_opp
      INTEGER_T im1,im2
      INTEGER_T ml,mr
      INTEGER_T local_mask
      INTEGER_T ispec
      REAL_T frac_pair(num_materials,num_materials)
      REAL_T dist_pair(num_materials,num_materials)
      REAL_T areaface
      INTEGER_T LS_change_sign(num_materials)
      INTEGER_T LS_plus(num_materials)
      INTEGER_T LS_minus(num_materials)
      REAL_T LScen(num_materials)
      REAL_T LSside(num_materials)
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
      INTEGER_T i_mdot
        !weights for BLB_MATRIX:
        !(1) weight=1 if LS>dx  weight=.001 otherwise
        !(2) weight=1 if at least one stencil cell is dominated by a solid.
        !    weight=0 otherwise.
        !(3) weight=1 if LS>0   weight=.001 otherwise
      REAL_T local_interior_wt(3) 
      REAL_T local_VEL(SDIM)
      REAL_T LS_clamped
      REAL_T VEL_clamped(SDIM)
      REAL_T temperature_clamped
      INTEGER_T prescribed_flag
      REAL_T blob_cell_count
      REAL_T blob_cellvol_count
      REAL_T blob_mass
      REAL_T blob_volume
      REAL_T blob_pressure
      REAL_T mdot_total
      REAL_T mdot_avg
      REAL_T mdot_part
      REAL_T updated_density
      REAL_T original_density
      REAL_T density_factor
      REAL_T mofdata(num_materials*ngeom_recon)
      INTEGER_T nmax
      INTEGER_T im_alt,im_mdot,im_opp_mdot
      INTEGER_T im_source,im_dest
      INTEGER_T im_evenly
      INTEGER_T iten,iten_shift
      REAL_T LL
      INTEGER_T ireverse
      INTEGER_T complement_flag
      INTEGER_T im_negate

      INTEGER_T local_material_type

      if ((tid_current.lt.0).or.(tid_current.ge.geom_nthreads)) then
       print *,"tid_current invalid"
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: fort_getcolorsum

      snew_ptr=>snew
      mdot_ptr=>mdot
      mdot_complement_ptr=>mdot_complement

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid: fort_getcolorsum"
       stop
      endif
      if (cur_time_slab.ge.zero) then
       ! do nothing
      else
       print *,"cur_time_slab invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid92"
       stop
      endif

      if (nface_dst.ne.num_materials*num_materials*2) then
       print *,"nface_dst invalid"
       stop
      endif
      if (ncellfrac.ne.num_materials*num_materials*(3+SDIM)) then
       print *,"ncellfrac invalid"
       stop
      endif
      if (ncomp_mdot.eq.0) then
       if (ncomp_mdot_alloc.eq.1) then
        ! do nothing
       else
        print *,"ncomp_mdot_alloc invalid"
        stop
       endif
      else if (ncomp_mdot.ge.1) then
       if (ncomp_mdot_alloc.eq.ncomp_mdot) then
        ! do nothing
       else
        print *,"ncomp_mdot_alloc invalid"
        stop
       endif
       if (ncomp_mdot_alloc.eq.2*num_interfaces) then
        ! do nothing
       else
        print *,"ncomp_mdot_alloc invalid"
        stop
       endif
      else
       print *,"ncomp_mdot invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid get color sum"
       stop
      endif

      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid in LEVELSET_3D.F90 "
       print *,"nstate=",nstate
       print *,"STATE_NCOMP=",STATE_NCOMP
       stop
      endif

       ! presently in: fort_getcolorsum
      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=DXMAXLS

      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,mdot_ptr,0,-1)
      call checkbound_array(fablo,fabhi,mdot_complement_ptr,0,-1)

      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      VEL_ptr=>VEL
      call checkbound_array(fablo,fabhi,VEL_ptr,1,-1)
      DEN_ptr=>DEN
      call checkbound_array(fablo,fabhi,DEN_ptr,1,-1)
      VOF_ptr=>VOF
      call checkbound_array(fablo,fabhi,VOF_ptr,1,-1)
      xface_ptr=>xface
      yface_ptr=>yface
      zface_ptr=>zface
      call checkbound_array(fablo,fabhi,xface_ptr,0,0)
      call checkbound_array(fablo,fabhi,yface_ptr,0,1)
      call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1)
      areax_ptr=>areax
      areay_ptr=>areay
      areaz_ptr=>areaz
      call checkbound_array1(fablo,fabhi,areax_ptr,0,0)
      call checkbound_array1(fablo,fabhi,areay_ptr,0,1)
      call checkbound_array1(fablo,fabhi,areaz_ptr,0,SDIM-1)
      cellfab_ptr=>cellfab
      call checkbound_array(fablo,fabhi,cellfab_ptr,0,-1)
      typefab_ptr=>typefab
      call checkbound_array1(fablo,fabhi,typefab_ptr,1,-1)
      color_ptr=>color
      call checkbound_array1(fablo,fabhi,color_ptr,1,-1)
      mask_ptr=>mask
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
  
      if (arraysize.ne.num_elements_blobclass*num_colors) then
       print *,"arraysize invalid"
       stop
      endif
       
      if (mdot_arraysize.eq.ncomp_mdot_alloc*num_colors) then
       ! do nothing
      else
       print *,"mdot_arraysize invalid"
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
         print *,"icolor invalid in fort_getcolorsum icolor=",icolor
         print *,"i,j,k ",i,j,k
         stop
        endif
        base_type=NINT(typefab(D_DECL(i,j,k)))
        if ((base_type.lt.1).or.(base_type.gt.num_materials)) then
         print *,"base_type invalid"
         stop
        endif

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do im=1,num_materials*ngeom_recon
         mofdata(im)=VOF(D_DECL(i,j,k),im)
        enddo

         ! tessellate==1:
         !  input: fluids tessellate, solids embedded 
         !  output: all materials tessellating 
         ! tessellate==3:
         !   a) if solid_vfrac>=1/2 then
         !        consider cell as F_{im_solid_max}=1
         !   b) else, only consider fluids.
        if ((tessellate.eq.1).or. &
            (tessellate.eq.3)) then
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
        else
         print *,"tessellate invalid"
         stop
        endif

        if (operation_flag.eq.OP_GATHER_MDOT) then

         if (level_blobtypedata(icolor).eq.0) then
          level_blobtypedata(icolor)=base_type
         else
          if (level_blobtypedata(icolor).ne.base_type) then
           print *,"type problems in fort_getcolorsum"
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

        else if (operation_flag.eq.OP_SCATTER_MDOT) then

         if (sweep_num.eq.0) then

          if (level_blobtypedata(icolor).ne.base_type) then
           print *,"type problems fort_getcolorsum OP_SCATTER_MDOT"
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

         else
          print *,"sweep_num invalid"
          stop
         endif
        else
         print *,"operation_flag invalid fort_getcolorsum"
         stop
        endif

        do im=1,num_materials
         opposite_color(im)=0
        enddo
        opposite_color(base_type)=icolor

        do im=1,num_materials
         LS_change_sign(im)=0
         LS_plus(im)=0
         LS_minus(im)=0
        enddo

        do im=1,num_materials
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

         call gridsten_level(xsten_stencil,i+i1,j+j1,k+k1,level,nhalf)

         do dir=1,SDIM
          local_VEL(dir)=VEL(D_DECL(i+i1,j+j1,k+k1),dir)
         enddo

         local_solid=0

         do im=1,num_materials
          LSside(im)=LS(D_DECL(i+i1,j+j1,k+k1),im)
          if (LScen(im)+LSside(im).ge.zero) then
           LS_plus(im)=1
          else if (LScen(im)+LSside(im).lt.zero) then
           LS_minus(im)=1
          else
           print *,"LScen or LSside invalid"
           stop
          endif
         enddo ! im=1..num_materials

         call get_primary_material(LSside,im_side_majority)
         if (is_rigid(im_side_majority).eq.1) then
          local_solid=1
         else if (is_rigid(im_side_majority).eq.0) then
          ! do nothing
         else
          print *,"is_rigid(im_side_majority) invalid"
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

         do dir=1,SDIM
          xstencil_point(dir)=xsten_stencil(0,dir)
         enddo

         call SUB_clamped_LS(xstencil_point,cur_time_slab,LS_clamped, &
           VEL_clamped,temperature_clamped,prescribed_flag,dx)

         if (LS_clamped.ge.zero) then
          do dir=1,SDIM
           local_VEL(dir)=VEL_clamped(dir)
          enddo
          local_solid=1
         else if (LS_clamped.le.zero) then
          ! do nothing
         else
          print *,"LS_clamped is NaN"
          stop
         endif

         if (typeside.eq.0) then
          ! do nothing
         else if (typeside.eq.base_type) then
           ! do nothing
         else if ((typeside.ge.1).and.(typeside.le.num_materials)) then
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
           solid_velocity(dir)=solid_velocity(dir)+local_VEL(dir)
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

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
         vol,cencell,SDIM)

        mass=zero
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         vfrac=mofdata(vofcomp)
         if (vfrac.ge.-VOFTOL) then
          if (vfrac.le.one+VOFTOL) then
           if (vfrac.lt.zero) then
            vfrac=zero
           endif
           if (vfrac.gt.one) then
            vfrac=one
           endif
           dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
           if (constant_density_all_time(im).eq.1) then
            den_mat=fort_denconst(im)
           else if (constant_density_all_time(im).eq.0) then
            den_mat=DEN(D_DECL(i,j,k),dencomp)
           else
            print *,"constant_density_all_time(im) invalid"
            stop
           endif
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
        enddo ! im=1..num_materials

        mass=mass*vol
        if ((mass.gt.zero).and.(mass.le.1.0D+20)) then
         ! do nothing
        else
         print *,"mass: floating point bust"
         stop
        endif

        if ((base_type.ge.1).and. &
            (base_type.le.num_materials)) then

         do im=1,num_materials
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
         enddo ! im=1..num_materials

         do dir=1,SDIM
          dx_sten(dir)=xsten(1,dir)-xsten(-1,dir)
          if (dx_sten(dir).le.zero) then
           print *,"dx_sten invalid"
           stop
          endif
         enddo
         RR=xsten(0,1)
         if (levelrz.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if (levelrz.eq.COORDSYS_RZ) then
          ! do nothing
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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
          if (levelrz.eq.COORDSYS_CARTESIAN) then
           dperim=one
          else if (levelrz.eq.COORDSYS_RZ) then
           dperim=two*Pi*RR
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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
         do im=1,num_materials
          do im_opp=1,num_materials
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

         do im=1,num_materials
          if ((opposite_color(im).ge.1).and. &
              (opposite_color(im).le.num_colors)) then

           ic_center=(opposite_color(im)-1)*num_elements_blobclass+ &
                   BLB_CEN_ACT+1

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
           ic=ic_base+BLB_MATRIX+1

           if (operation_flag.eq.OP_GATHER_MDOT) then

             ! local_interior_wt(1): avoid noisy velocity conditions near the
             ! interface.
            if (LScen(im).ge.cutoff) then ! cutoff=DXMAXLS
             local_interior_wt(1)=one
            else if (LScen(im).lt.cutoff) then
             local_interior_wt(1)=1.0E-3
            else
             print *,"LScen(im) invalid"
             stop
            endif

             ! solid_fraction=1:
             !  At least one 3x3x3 stencil cell is dominated by a solid.
             ! solid_fraction=0:
             !  Otherwise.

            local_interior_wt(2)=solid_fraction

             ! local_interior_wt(3): take into account all material im cells
             ! including those near the interface.
            if (LScen(im).ge.zero) then
             local_interior_wt(3)=one
            else if (LScen(im).lt.zero) then
             local_interior_wt(3)=1.0E-3
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
               dotprod=dotprod*local_interior_wt(veltype)
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

              dotprod=dotprod*local_interior_wt(veltype)

              level_blobdata(ic)=level_blobdata(ic)+mass*dotprod
              ic=ic+1
             enddo ! irow=1..2*sdim
            enddo ! veltype=1..3

            ic=ic_base+BLB_INT_MOM+1  

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
             dotprod=dotprod*local_interior_wt(veltype)
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
             dotprod=dotprod*local_interior_wt(veltype)
             level_blobdata(ic)=level_blobdata(ic)+mass*dotprod
             ic=ic+1
            enddo ! irow=1..2 * sdim

             ! blob_energy
            dotprod=zero
            do dir=1,SDIM
             dotprod=dotprod+fluid_velocity(dir)**2
            enddo
            dotprod=dotprod*local_interior_wt(veltype)
            level_blobdata(ic)=level_blobdata(ic)+half*mass*dotprod
            ic=ic+1

             ! blob_mass_for_velocity
            do veltype=1,3
             dotprod=local_interior_wt(veltype)
             level_blobdata(ic)=level_blobdata(ic)+mass*dotprod
             ic=ic+1
            enddo ! veltype=1..3

            if (ic.ne.ic_base+BLB_VOL+1) then
             print *,"ic invalid ic=",ic
             print *,"expecting ic=",ic_base+BLB_VOL+1
             stop
            endif

             ! blob_volume
            vofcomp=(im-1)*ngeom_recon+1
            vfrac=mofdata(vofcomp)
            level_blobdata(ic)=level_blobdata(ic)+vol*vfrac

             ! centroid integral
            ic=ic+1
            do dir=1,SDIM
             level_blobdata(ic)=level_blobdata(ic)+ &
              (cencell(dir)+mofdata(vofcomp+dir))*vfrac*vol
             ic=ic+1
            enddo

            ! centroid actual
            ic=ic+SDIM
            
             ! perimeter (internal faces)
            do im_opp=1,num_materials
             if (im_opp.ne.im) then
              level_blobdata(ic)=level_blobdata(ic)+local_facearea(im,im_opp) 
             endif
            enddo

             ! perimeter decomposed by material.
            ic=ic+1
            do im_opp=1,num_materials
             if (im_opp.ne.im) then
              level_blobdata(ic)=level_blobdata(ic)+local_facearea(im,im_opp) 
             endif
             ic=ic+1
            enddo ! im_opp

             ! contact line perimeter
            do im1=1,num_materials
            do im2=1,num_materials
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
            enddo ! im2=1..num_materials
            enddo ! im1=1..num_materials

            if (ic.eq. &
                (opposite_color(im)-1)*num_elements_blobclass+ &
                 BLB_CELL_CNT+1) then
             ! do nothing
            else
             print *,"ic invalid for blob_cell_count"
             print *,"blob_cell_count,blob_cellvol_count,"
             print *,"blob_mass,blob_pressure,"
             print *,"ic=",ic
             print *,"im=",im
             print *,"opposite_color(im)=",opposite_color(im)
             print *,"num_elements_blobclass=",num_elements_blobclass
             stop
            endif

             ! blob_cell_count  (ic)
             ! blob_cellvol_count (ic+1)
             ! blob_mass (ic+2)
             ! blob_pressure (ic+3)
            if (vfrac.ge.half) then
             level_blobdata(ic)=level_blobdata(ic)+one !blob_cell_count
             level_blobdata(ic+1)=level_blobdata(ic+1)+vol !blob_cellvol_count

             pressure_local=VEL(D_DECL(i,j,k),STATECOMP_PRES+1)
             if (is_rigid(im).eq.0) then

              if ((fort_material_type(im).ge.1).and. &
                  (fort_material_type(im).le.MAX_NUM_EOS)) then 
               local_material_type=fort_material_type(im)
              else if (fort_material_type(im).eq.0) then

               if ((material_type_lowmach(im).ge.1).and. &
                   (material_type_lowmach(im).le.MAX_NUM_EOS)) then
                local_material_type=material_type_lowmach(im)
               else if (material_type_lowmach(im).eq.0) then

                if ((material_type_visual(im).ge.1).and. &
                    (material_type_visual(im).le.MAX_NUM_EOS)) then
                 local_material_type=material_type_visual(im)
                else if (material_type_visual(im).eq.0) then
                 local_material_type=0
                else
                 print *,"material_type_visual(im) invalid"
                 stop
                endif

               else
                print *,"material_type_lowmach(im) invalid"
                stop
               endif

              else
               print *,"fort_material_type(im) invalid"
               stop
              endif


              if (local_material_type.eq.0) then
               ! do nothing
              else if ((local_material_type.ge.1).and. &
                       (local_material_type.le.MAX_NUM_EOS)) then
               dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
               
               if (constant_density_all_time(im).eq.1) then
                den_mat=fort_denconst(im)
               else if (constant_density_all_time(im).eq.0) then
                den_mat=DEN(D_DECL(i,j,k),dencomp)
               else
                print *,"constant_density_all_time(im) invalid"
                stop
               endif
               if (den_mat.gt.zero) then
                ! do nothing
               else
                print *,"den_mat has gone nonpos"
                stop
               endif
               TEMP_mat=DEN(D_DECL(i,j,k),dencomp+1)
               if (TEMP_mat.gt.zero) then
                ! do nothing
               else
                print *,"TEMP_mat has gone nonpos"
                stop
               endif
               call init_massfrac_parm(den_mat,massfrac_parm,im)
               do ispec=1,num_species_var
                massfrac_parm(ispec)=DEN(D_DECL(i,j,k),dencomp+1+ispec)
               enddo
               call INTERNAL_material(den_mat,massfrac_parm,TEMP_mat, &
                internal_energy,local_material_type,im)
               if (internal_energy.gt.zero) then
                ! do nothing
               else
                print *,"internal_energy has gone nonpos"
                stop
               endif
               call EOS_material(den_mat,massfrac_parm, &
                internal_energy, &
                pressure_local, &
                local_material_type,im)
              else
               print *,"local_material_type invalid"
               stop
              endif

             else if (is_rigid(im).eq.1) then
              ! do nothing
             else
              print *,"is_rigid(im) invalid"
              stop
             endif

             if (ic+3.eq. &
                 (opposite_color(im)-1)*num_elements_blobclass+ &
                  BLB_PRES+1) then
              level_blobdata(ic+3)=level_blobdata(ic+3)+ &
                vol*pressure_local !blob_pressure
             else
              print *,"expecting ic+3 to correspond to BLB_PRES+1"
              stop
             endif

             if (ncomp_mdot.eq.2*num_interfaces) then
              ic_base_mdot=(opposite_color(im)-1)*ncomp_mdot
              do i_mdot=1,ncomp_mdot
               level_mdot_data(ic_base_mdot+i_mdot)= &
                  level_mdot_data(ic_base_mdot+i_mdot)+ &
                  mdot(D_DECL(i,j,k),i_mdot)

                ! mdot complement
               level_mdot_complement_data(ic_base_mdot+i_mdot)= &
                  level_mdot_complement_data(ic_base_mdot+i_mdot)+ &
                  mdot_complement(D_DECL(i,j,k),i_mdot)
              enddo
             else if (ncomp_mdot.eq.0) then
              ! do nothing
             else
              print *,"ncomp_mdot invalid"
              stop
             endif

            else if (vfrac.lt.half) then
             ! do nothing
            else
             print *,"vfrac is NaN"
             stop
            endif

            ! blob_cell_count  (ic)
            ! blob_cellvol_count (ic+1)
            ! blob_mass (ic+2)
            ! blob_pressure (ic+3)
            ic=ic+3

            if (ic.eq. &
                (opposite_color(im)-1)*num_elements_blobclass+ &
                 BLB_PRES+1) then
             ! do nothing
            else
             print *,"ic invalid, blob_mass is last?"
             print *,"ic=",ic
             print *,"im=",im
             print *,"opposite_color(im)=",opposite_color(im)
             print *,"num_elements_blobclass=",num_elements_blobclass
             stop
            endif

             ! blob_mass
            dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
            if (constant_density_all_time(im).eq.1) then
             den_mat=fort_denconst(im)
            else if (constant_density_all_time(im).eq.0) then
             den_mat=DEN(D_DECL(i,j,k),dencomp)
            else
             print *,"constant_density_all_time(im) invalid"
             stop
            endif
            if (den_mat.ge.(one-VOFTOL)*fort_density_floor(im)) then
             if (den_mat.le.(one+VOFTOL)*fort_density_ceiling(im)) then
              ! blob_cell_count  (ic-3)
              ! blob_cellvol_count (ic-2)
              ! blob_mass (ic-1)
              ! blob_pressure (ic)
              level_blobdata(ic-1)=level_blobdata(ic-1)+vol*vfrac*den_mat
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

            ic=ic+1

            if (ic.eq. &
                (opposite_color(im)-1)*num_elements_blobclass+ &
                 BLB_PRES+2) then
             ! do nothing
            else
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
               do ml = 1, num_materials
               do mr = 1, num_materials
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
               do ml = 1, num_materials
               do mr = 1, num_materials
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
               do ml = 1, num_materials
               do mr = 1, num_materials
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
              ic=(opposite_color(im)-1)*num_elements_blobclass+BLB_PERIM+1

               ! ic component is blob_perim
              do im_opp=1,num_materials
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
                level_blobdata(ic)=level_blobdata(ic)+ &
                        frac_pair(ml,mr)*areaface 
               endif
              enddo !im_opp=1..num_materials

               ! perimeter decomposed (num_materials components)
              ic=ic+1
              do im_opp=1,num_materials
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
                level_blobdata(ic)=level_blobdata(ic)+ &
                        frac_pair(ml,mr)*areaface 
               endif
               ic=ic+1
              enddo ! im_opp=1..num_materials
              if (ic.eq. &
                  (opposite_color(im)-1)*num_elements_blobclass+ &
                   BLB_TRIPLE_PERIM+1) then
               ! do nothing
              else
               print *,"ic invalid for blob_cell_count"
               stop
              endif
  
             enddo ! side
            enddo ! dir

           else if (operation_flag.eq.OP_SCATTER_MDOT) then

            if (sweep_num.eq.0) then

              ! blob_volume
             vofcomp=(im-1)*ngeom_recon+1
             vfrac=mofdata(vofcomp)

             if (1.eq.0) then
              print *,"i,j,k,im,vfrac (before if) ",i,j,k,im,vfrac
             endif

             do im_alt=1,num_materials
              if (im_alt.ne.im) then
               do ireverse=0,1
                if (im_alt.lt.im) then
                 im_mdot=im_alt
                 im_opp_mdot=im
                else if (im_alt.gt.im) then
                 im_mdot=im
                 im_opp_mdot=im_alt
                else
                 print *,"im_alt bust"
                 stop
                endif
                call get_iten(im_mdot,im_opp_mdot,iten)
                iten_shift=ireverse*num_interfaces+iten
                LL=get_user_latent_heat(iten_shift,293.0d0,1)
                if (LL.eq.zero) then
                 ! do nothing
                else if (LL.ne.zero) then
                 if (ireverse.eq.0) then
                  im_source=im_mdot
                  im_dest=im_opp_mdot
                 else if (ireverse.eq.1) then
                  im_source=im_opp_mdot
                  im_dest=im_mdot
                 else
                  print *,"ireverse bust"
                  stop
                 endif

                 if (distribute_mdot_evenly(iten_shift).eq.0) then
                  ! do nothing
                 else if ((distribute_mdot_evenly(iten_shift).eq.1).or. &
                          (distribute_mdot_evenly(iten_shift).eq.2)) then

                  if (vfrac.ge.half) then

                   if (distribute_from_target(iten_shift).eq.0) then
                    im_evenly=im_dest
                   else if (distribute_from_target(iten_shift).eq.1) then
                    im_evenly=im_source
                   else
                    print *,"distribute_from_target(iten_shift) invalid"
                    stop
                   endif
              
                    ! sum alpha V_i = mdot_total
                    ! alpha=mdot_total/(sum V_i) 
                   if (im.eq.im_evenly) then 
                    ic=(opposite_color(im)-1)*num_elements_blobclass+ &
                            BLB_CELL_CNT+1
                    blob_cell_count=cum_blobdata(ic)
                    blob_cellvol_count=cum_blobdata(ic+1)
                    if ((blob_cellvol_count.gt.zero).and. &
                        (blob_cell_count.gt.zero)) then

                     if (ncomp_mdot.eq.2*num_interfaces) then
                      ic_base_mdot=(opposite_color(im)-1)*ncomp_mdot
                      mdot_total=cum_mdot_data(ic_base_mdot+iten_shift)

                      if (distribute_mdot_evenly(iten_shift).eq.1) then
                       mdot_avg=mdot_total/blob_cellvol_count
                        ! divu_cor must be a constant since 
                        ! rho_cor must be a constant.
                       mdot_part=mdot_avg*vol ! divu_cor * volume
                      else if (distribute_mdot_evenly(iten_shift).eq.2) then
                       print *,"distribute_mdot_evenly=2 not allowed"
                       print *,"only distribute_mdot_evenly=1 allowed since"
                       print *,"divu correction must be consistent with the"
                       print *,"uniform density correction requirement." 
                       stop
                       mdot_avg=mdot_total/blob_cell_count
                       mdot_part=mdot_avg
                      else
                       print *,"distribute_mdot_evenly(iten_shift) invalid"
                       stop
                      endif

                      level_mdot_data_redistribute(ic_base_mdot+iten_shift)= &
                       level_mdot_data_redistribute(ic_base_mdot+iten_shift)+ &
                       mdot_part
                 
                      if (fort_material_type(im).eq.0) then
                       mdot(D_DECL(i,j,k),iten_shift)=mdot_part
                      else if ((fort_material_type(im).gt.0).and. &
                               (fort_material_type(im).le.MAX_NUM_EOS)) then
                       print *,"phase change only for incompressible materials"
                       stop
                      else 
                       print *,"fort_material_type(im) invalid"
                       stop
                      endif
                     else
                      print *,"ncomp_mdot invalid"
                      stop
                     endif

                    else
                     print *,"blob_cell_count or blob_cellvol_count invalid"
                     stop
                    endif
                   else if (im.ne.im_evenly) then
                    ! do nothing
                   else
                    print *,"im or im_evenly bust"
                    stop
                   endif

                  else if (vfrac.lt.half) then
                   ! do nothing
                  else
                   print *,"vfrac invalid"
                   stop
                  endif

                 else
                  print *,"distribute_mdot_evenly(iten_shift) invalid"
                  stop
                 endif

                 if (constant_volume_mdot(iten_shift).eq.0) then
                  im_negate=0
                  complement_flag=0 
                 else if (constant_volume_mdot(iten_shift).eq.1) then
                  ! distribute -sum mdot to the source:
                  im_negate=im_source
                  if (distribute_from_target(iten_shift).eq.0) then
                   complement_flag=1 
                  else if (distribute_from_target(iten_shift).eq.1) then
                   complement_flag=0 
                  else
                   print *,"distribute_from_target(iten_shift) invalid"
                   stop
                  endif
                 else if (constant_volume_mdot(iten_shift).eq.-1) then
                  ! distribute -sum mdot to the dest:
                  im_negate=im_dest
                  if (distribute_from_target(iten_shift).eq.0) then
                   complement_flag=0 
                  else if (distribute_from_target(iten_shift).eq.1) then
                   complement_flag=1 
                  else
                   print *,"distribute_from_target(iten_shift) invalid"
                   stop
                  endif
                 else
                  print *,"constant_volume_mdot(iten_shift) invalid"
                  stop
                 endif

                 if (im_negate.eq.0) then
                  ! do nothing
                 else if (im_negate.eq.im) then
                  if (constant_density_all_time(im).eq.1) then
                   print *,"constant_density_all_time(im) invalid"
                   stop
                  else if (constant_density_all_time(im).eq.0) then
                   ! do nothing
                  else
                   print *,"constant_density_all_time(im) invalid"
                   stop
                  endif

                  if (1.eq.0) then
                   print *,"i,j,k,im,im_negate,vfrac ", &
                           i,j,k,im,im_negate,vfrac
                   print *,"iten_shift,im_mdot,im_opp_mdot ", &
                    iten_shift,im_mdot,im_opp_mdot
                   print *,"complement_flag ",complement_flag
                  endif
 
                  ic=(opposite_color(im)-1)*num_elements_blobclass+ &
                          BLB_CELL_CNT+1

                  blob_cell_count=cum_blobdata(ic)
                  blob_cellvol_count=cum_blobdata(ic+1)
                  blob_mass=cum_blobdata(ic+2)
                  blob_pressure=cum_blobdata(ic+3)
                  ic=ic+3

                  if (ic_base.eq. &
                      (opposite_color(im)-1)*num_elements_blobclass) then
                   ! do nothing
                  else
                   print *,"ic_base invalid"
                   stop
                  endif
                  if (ic.eq. &
                      (opposite_color(im)-1)*num_elements_blobclass+ &
                       BLB_PRES+1) then
                   ! do nothing
                  else
                   print *,"ic invalid"
                   stop
                  endif

                  ic=ic_base+BLB_VOL+1

                  blob_volume=cum_blobdata(ic)

                  if ((blob_cell_count.gt.zero).and. &
                      (blob_cellvol_count.gt.zero).and. &
                      (blob_mass.gt.zero).and. &
                      (blob_volume.gt.zero)) then

                   if (ncomp_mdot.eq.2*num_interfaces) then
                    ic_base_mdot=(opposite_color(im)-1)*ncomp_mdot

                    original_density=blob_mass/blob_volume
                    if (original_density.gt.zero) then
                     ! do nothing
                    else
                     print *,"original_density invalid"
                     stop
                    endif

                    if (complement_flag.eq.0) then
                     mdot_total=cum_mdot_data(ic_base_mdot+iten_shift)
                     density_factor=one
                    else if (complement_flag.eq.1) then
                     if (distribute_from_target(iten_shift).eq.0) then
                      density_factor=fort_denconst(im_dest)/original_density
                     else if (distribute_from_target(iten_shift).eq.1) then
                      density_factor=fort_denconst(im_source)/original_density
                     else
                      print *,"distribute_from_target(iten_shift) invalid"
                      stop
                     endif
                     mdot_total= &
                          cum_mdot_complement_data(ic_base_mdot+iten_shift)
                     if (1.eq.0) then
                      print *,"i,j,k,cell_count,mass,volume,mdot_tot ", &
                       i,j,k,blob_cell_count,blob_mass,blob_volume,mdot_total
                      print *,"i,j,k,cellvol_count ",i,j,k,blob_cellvol_count
                     endif
                    else
                     print *,"complement_flag invalid"
                     stop
                    endif
                    mdot_avg=mdot_total/blob_cellvol_count

                    if (vfrac.ge.half) then
                     level_mdot_data_redistribute(ic_base_mdot+iten_shift)= &
                      level_mdot_data_redistribute(ic_base_mdot+iten_shift)+ &
                      mdot_avg*vol

                     mdot(D_DECL(i,j,k),iten_shift)= &
                       mdot(D_DECL(i,j,k),iten_shift)-mdot_avg*vol
                    else if (vfrac.lt.half) then
                     ! do nothing
                    else
                     print *,"vfrac invalid"
                     stop
                    endif

                     ! F_t + s|grad F|=0
                     ! dF=F_t dt=-s|grad F|dt=-s dt/dx
                     ! my mdot=(den_src/den_dst-1)*dF*Vcell/dt^2=
                     !  (den_src/den_dst-1)*(-s * area)/dt
                     ! units of my mdot are "velocity" * area / seconds=
                     ! "div velocity" * Length * area / seconds =
                     ! "div velocity" * volume / seconds
                     ! distribute mdot so that rho div u=constant 
                     ! since rho=constant, make sure (div u)_correct=const.
                     ! mdot_correct=divu_cor * volume/seconds
                     ! sum mdot = sum mdot_correct
                     ! divu_cor=sum mdot/(sum volume_i)=
                     ! sum (divu_source_i vol_i)/(sum vol_i)
                     ! so, we should always have:
                     !  distribute_mdot_evenly(iten_shift).eq.1
                     ! 
                     ! mass_new-mass_old = sum_i vel_i dt * area_i * 
                     !                     (den_dst-den_src) =
                     !                     sum_i dF*Vcell*(den_dst-den_src)=
                     !                     DM
                     ! if distribute_from_target==0,
                     !  mdot=(den_src/den_dst-1)*dF*Vcell/dt^2
                     !  sum_i mdot_i=sum (den_src-den_dst)*dF*Vcell/
                     !                   (den_dst*dt^2)=-DM/(den_dst*dt^2)
                     !  rho_update-=DM/V_update
                     !  rho_update+=sum_i mdot_i*den_dst * dt^2/V_update
                     ! if distribute_from_target==1,
                     !  mdot=(1-den_dst/den_src)*dF*Vcell/dt^2
                     !  sum_i mdot_i=sum (den_src-den_dst)*dF*Vcell/
                     !                   (den_src*dt^2)=-DM/(den_src*dt^2)
                     !  rho_update-=DM/V_update
                     !  rho_update+=sum_i mdot_i*den_src * dt^2/V_update

                    if (original_density.gt.zero) then
                     updated_density=original_density* &
                         (one+dt*dt*density_factor*mdot_total/blob_volume)
                     dencomp=STATECOMP_STATES+ &
                       (im-1)*num_state_material+1+ENUM_DENVAR
                     if (updated_density.gt.zero) then
                      snew(D_DECL(i,j,k),dencomp)=updated_density
                     else
                      print *,"updated_density invalid"
                      stop
                     endif
                    else
                     print *,"original_density invalid"
                     stop
                    endif
                   else
                    print *,"ncomp_mdot invalid"
                    stop
                   endif
                  else if ((blob_cell_count.ge.zero).and. &
                           (blob_cellvol_count.ge.zero).and. &
                           (blob_mass.ge.zero).and. &
                           (blob_volume.ge.zero)) then
                   if (1.eq.0) then
                    print *,"WARNING for blob_cell_count,blob_cellvol_count"
                    print *,"blob_cell_count=",blob_cell_count
                    print *,"blob_cellvol_count=",blob_cellvol_count
                    print *,"blob_mass=",blob_mass
                    print *,"blob_volume=",blob_volume
                    print *,"im=",im
                    print *,"opposite_color(im)= (1..ncolors) ", &
                           opposite_color(im)
                   endif
                  else
                   print *,"blob_cell_count,or ..."
                   print *,"blob_cellvol_count,or ..."
                   print *,"blob_mass,or ... "
                   print *,"blob_volume bad"
                   print *,"blob_cell_count=",blob_cell_count
                   print *,"blob_cellvol_count=",blob_cellvol_count
                   print *,"blob_mass=",blob_mass
                   print *,"blob_volume=",blob_volume
                   print *,"im=",im
                   print *,"opposite_color(im)= (1..ncolors) ", &
                           opposite_color(im)
                   stop
                  endif
                 else if (im.ne.im_negate) then
                  ! do nothing
                 else
                  print *,"im or im_negate bust"
                  stop
                 endif

                else
                 print *,"latent_heat(iten_shift) invalid"
                 stop
                endif
               enddo ! ireverse=0...1
              else if (im_alt.eq.im) then
               ! do nothing
              else
               print *,"im_alt or im bust"
               stop
              endif
             enddo ! im_alt=1..num_materials

            else
             print *,"sweep_num invalid"
             stop
            endif

           else
            print *,"operation_flag invalid",operation_flag
            print *,"operation_flag options:"
            print *,"OP_GATHER_MDOT ",OP_GATHER_MDOT
            print *,"OP_SCATTER_MDOT ",OP_SCATTER_MDOT
            stop
           endif

          else if (opposite_color(im).eq.0) then
           ! do nothing
          else
           print *,"opposite_color invalid"
           stop
          endif

         enddo ! im=1..num_materials

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
      end subroutine fort_getcolorsum

      subroutine fort_get_lowmach_divu( &
       tid_current, &
       sweep_num, & !sweep_num=0: sum V_T rho DT/Dt, sweep_num=1:sum mdot=0
       constant_density_all_time, & ! 1..num_materials
       dt, &
       dx, &
       xlo, &
       mdot_local, & ! mdot_local initialized to 0
       DIMS(mdot_local), &
       mdot_global, & ! mdot_global incremented
       DIMS(mdot_global), &
       LS,DIMS(LS), &
       DEN,DIMS(DEN), &
       DTDt,DIMS(DTDt), & ! T_after_diffusion-T_after_advection
       VOF,DIMS(VOF), &
       typefab,DIMS(typefab), &
       color,DIMS(color), &
       mask,DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       rzflag, &
       num_colors, &
       level_blobdata, &
       level_mdot_data, &
       arraysize, & ! size of level_blobdata
       mdot_arraysize, & ! size of level_mdot_data
       material_type_lowmach) &
      bind(c,name='fort_get_lowmach_divu')

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: tid_current
      INTEGER_T, INTENT(in) :: sweep_num
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      REAL_T, INTENT(in) :: dt
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: xlo(SDIM)
      INTEGER_T, INTENT(in) :: material_type_lowmach(num_materials)
      INTEGER_T, INTENT(in) :: constant_density_all_time(num_materials)

      INTEGER_T, INTENT(in) :: rzflag
      INTEGER_T, INTENT(in) :: num_colors
      INTEGER_T, INTENT(in) :: arraysize
      INTEGER_T, INTENT(in) :: mdot_arraysize
      INTEGER_T, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      INTEGER_T :: growlo(3), growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: DIMDEC(mdot_local)
      INTEGER_T, INTENT(in) :: DIMDEC(mdot_global)
      INTEGER_T, INTENT(in) :: DIMDEC(LS)
      INTEGER_T, INTENT(in) :: DIMDEC(DEN)
      INTEGER_T, INTENT(in) :: DIMDEC(DTDt)
      INTEGER_T, INTENT(in) :: DIMDEC(VOF)
      INTEGER_T, INTENT(in) :: DIMDEC(typefab)
      INTEGER_T, INTENT(in) :: DIMDEC(color)
      INTEGER_T, INTENT(in) :: DIMDEC(mask)
      REAL_T, INTENT(in) ::level_blobdata(arraysize)
      REAL_T, INTENT(inout) ::level_mdot_data(mdot_arraysize)

      REAL_T, INTENT(inout), target :: mdot_local(DIMV(mdot_local))
      REAL_T, pointer :: mdot_local_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(inout), target :: mdot_global(DIMV(mdot_global))
      REAL_T, pointer :: mdot_global_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: typefab(DIMV(typefab))
      REAL_T, pointer :: typefab_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: LS(DIMV(LS),num_materials*(1+SDIM))
      REAL_T, pointer :: LS_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: DEN(DIMV(DEN),num_materials*num_state_material)
      REAL_T, pointer :: DEN_ptr(D_DECL(:,:,:),:)
       ! T_after_diffusion-T_adv
      REAL_T, INTENT(in), target :: DTDt(DIMV(DTDt),num_materials) 
      REAL_T, pointer :: DTDt_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: VOF(DIMV(VOF),num_materials*ngeom_recon)
      REAL_T, pointer :: VOF_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: color(DIMV(color))
      REAL_T, pointer :: color_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))

      INTEGER_T :: i,j,k
      INTEGER_T :: dir
      INTEGER_T vofcomp
      INTEGER_T dencomp
      INTEGER_T base_type
      INTEGER_T ic
      INTEGER_T icolor
      REAL_T vfrac
      REAL_T volcell
      REAL_T den_mat
      REAL_T TEMP_mat
      REAL_T internal_energy
      REAL_T massfrac_parm(num_species_var+1)
      REAL_T pressure_local
      REAL_T pressure_sum
      REAL_T dVdT
      INTEGER_T, parameter :: nhalf=3
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T dx_sten(SDIM)
      REAL_T cencell(SDIM)
      INTEGER_T im
      INTEGER_T local_mask
      INTEGER_T ispec
      REAL_T LScen(num_materials)
      INTEGER_T im_primary
      REAL_T blob_cellvol_count
      REAL_T mdot_local_scalar
      REAL_T mdot_sum_denom
      REAL_T mdot_sum_numerator
      REAL_T DTDt_local
      REAL_T mofdata(num_materials*ngeom_recon)

      if ((tid_current.lt.0).or.(tid_current.ge.geom_nthreads)) then
       print *,"tid_current invalid"
       stop
      endif

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid92"
       stop
      endif

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid get lowmach divu"
       stop
      endif

      mdot_local_ptr=>mdot_local
      call checkbound_array1(fablo,fabhi,mdot_local_ptr,0,-1)
      mdot_global_ptr=>mdot_global
      call checkbound_array1(fablo,fabhi,mdot_global_ptr,0,-1)
      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      DEN_ptr=>DEN
      call checkbound_array(fablo,fabhi,DEN_ptr,0,-1)
      DTDt_ptr=>DTDt
      call checkbound_array(fablo,fabhi,DTDt_ptr,0,-1)
      VOF_ptr=>VOF
      call checkbound_array(fablo,fabhi,VOF_ptr,1,-1)
      color_ptr=>color
      call checkbound_array1(fablo,fabhi,color_ptr,1,-1)
      typefab_ptr=>typefab
      call checkbound_array1(fablo,fabhi,typefab_ptr,1,-1)
      mask_ptr=>mask
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
  
      if (arraysize.ne.num_elements_blobclass*num_colors) then
       print *,"arraysize invalid"
       stop
      endif
      
        ! need sum mdot_local_{ij}=0
        ! let mdot_local_{ij}=mdot_original_{ij} - alpha Volume
        ! sum mdot_local_{ij}= sum mdot_original_{ij} -alpha sum Volume=0
        ! alpha=sum mdot_original_{ij}/sum Volume
        ! sum is over FULL cells. 
      if (mdot_arraysize.eq.2*num_colors) then
       ! do nothing
      else
       print *,"mdot_arraysize <> 2 * num_colors"
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
         print *,"icolor invalid in fort_get_lowmach_divu icolor=",icolor
         print *,"i,j,k ",i,j,k
         stop
        endif
        base_type=NINT(typefab(D_DECL(i,j,k)))
        if ((base_type.lt.1).or.(base_type.gt.num_materials)) then
         print *,"base_type invalid"
         stop
        endif

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do im=1,num_materials*ngeom_recon
         mofdata(im)=VOF(D_DECL(i,j,k),im)
        enddo

        do im=1,num_materials
         LScen(im)=LS(D_DECL(i,j,k),im)
        enddo

        call get_primary_material(LScen,im_primary)

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
         volcell,cencell,SDIM)
        if (volcell.gt.zero) then
         ! do nothing
        else
         print *,"volcell invalid"
         stop
        endif

        vofcomp=(im_primary-1)*ngeom_recon+1
        vfrac=mofdata(vofcomp)
  
        if (vfrac.ge.VOFTOL) then
         if (is_rigid(im_primary).eq.0) then

           ! anelastic approximation only appropriate for "incompressible"
           ! model.
          if (fort_material_type(im_primary).eq.0) then
          
           if ((vfrac.ge.VOFTOL).and. &
               (vfrac.le.one-VOFTOL)) then
            ! do nothing
           else if ((vfrac.ge.one-VOFTOL).and. &
                    (vfrac.le.one+VOFTOL)) then

            if (im_primary.eq.base_type) then

             if (sweep_num.eq.0) then
              level_mdot_data(2*(icolor-1)+1)= &
               level_mdot_data(2*(icolor-1)+1)+volcell
             else if (sweep_num.eq.1) then
              ! do nothing
             else
              print *,"sweep_num invalid"
              stop
             endif

             dencomp=(im_primary-1)*num_state_material+1+ENUM_DENVAR
             if (constant_density_all_time(im_primary).eq.1) then
              den_mat=fort_denconst(im_primary)
             else if (constant_density_all_time(im_primary).eq.0) then
              den_mat=DEN(D_DECL(i,j,k),dencomp)
             else
              print *,"constant_density_all_time(im_primary) invalid"
              stop
             endif
             if (den_mat.gt.zero) then
              ! do nothing
             else
              print *,"den_mat out of range"
              stop
             endif

             do dir=1,SDIM
              dx_sten(dir)=xsten(1,dir)-xsten(-1,dir)
              if (dx_sten(dir).gt.zero) then
               ! do nothing
              else
               print *,"dx_sten invalid"
               stop
              endif
             enddo ! dir=1..sdim

             ic=(icolor-1)*num_elements_blobclass+BLB_PRES+1
             pressure_sum=level_blobdata(ic) ! sum pres * vol

             if (BLB_PRES-2.eq.BLB_CELLVOL_CNT) then
              blob_cellvol_count=level_blobdata(ic-2)
             else
              print *,"BLB_PRES or BLB_CELLVOL_CNT invalid"
              stop
             endif

             if (blob_cellvol_count.gt.zero) then

              if ((pressure_sum.ge.zero).or. &
                  (pressure_sum.le.zero)) then
               pressure_sum=pressure_sum/blob_cellvol_count
 
               if (material_type_lowmach(im_primary).eq.0) then
                ! do nothing
               else if ((material_type_lowmach(im_primary).gt.0).and. &
                        (material_type_lowmach(im_primary).le.MAX_NUM_EOS)) then

                if (pressure_sum.gt.zero) then

                 TEMP_mat=DEN(D_DECL(i,j,k),dencomp+1)
                 if (TEMP_mat.gt.zero) then
                  ! do nothing
                 else
                  print *,"TEMP_mat has gone nonpos:",TEMP_mat
                  stop
                 endif
                 call init_massfrac_parm(den_mat,massfrac_parm,im_primary)
                 do ispec=1,num_species_var
                  massfrac_parm(ispec)=DEN(D_DECL(i,j,k),dencomp+1+ispec)
                 enddo
                 call INTERNAL_material(den_mat,massfrac_parm,TEMP_mat, &
                  internal_energy, &
                  material_type_lowmach(im_primary),im_primary)
                 if (internal_energy.gt.zero) then
                  ! do nothing
                 else
                  print *,"internal_energy has gone nonpos:",internal_energy
                  stop
                 endif
                 call EOS_material(den_mat,massfrac_parm, &
                  internal_energy, &
                  pressure_local, &
                  material_type_lowmach(im_primary),im_primary)

                  ! for perfect gas:dVdT=(gamma-1)Cv/p
                  ! Cv units=J/(kg Kelvin)  J=Newton meter
                  ! (kinetic energy= mass m^2/sec^2=(kg m/s^2) * (m) )
                  ! pressure=Newton/m^2
                  ! Cv/p units = (Newton Meter)/(kg Kelvin)  / (Newton/m^2)=
                  ! m^3/(kg Kelvin) = 1/((kg/m^3) * Kelvin) 
                 call dVdT_material(dVdT,massfrac_parm, &
                  pressure_sum,TEMP_mat, &
                  material_type_lowmach(im_primary),im_primary)

                 if (dVdT.gt.zero) then
                  ! do nothing
                 else
                  print *,"expecting dVdT>0: ",dVdT
                  stop
                 endif

                 !F_t + s|grad F|=0
                 ! dF=F_t dt=-s|grad F|dt=-s dt/dx
                 ! my mdot=(den_src/den_dst-1)*dF*Volume_cell/dt^2=
                 !  (den_src/den_dst-1)*(-s * area)/dt
                 ! units of my mdot are "velocity" * area / seconds=
                 ! "div velocity" * Length * area / seconds =
                 ! "div velocity" * volume / seconds = (1/s)(cm^3)/s=cm^3/s^2
                 !
                 !DTDt=T_after_diffusion-T_after_advection
                 !
                 ! units of mdot_local_scalar: (1/density)(1/Kelvin)(density)
                 !   (Kelvin)*cm^3/s^2=cm^3/s^2
                 DTDt_local=DTDt(D_DECL(i,j,k),im_primary)
                 if ((DTDt_local.ge.zero).or.(DTDt_local.le.zero)) then
                  mdot_local_scalar=dVdT*den_mat*DTDt_local*volcell/(dt*dt)
                 else
                  print *,"DTDt_local is NaN"
                  stop
                 endif

                 if (sweep_num.eq.0) then
                  mdot_local(D_DECL(i,j,k))=mdot_local_scalar
                  level_mdot_data(2*(icolor-1)+2)= &
                   level_mdot_data(2*(icolor-1)+2)+mdot_local_scalar
                 else if (sweep_num.eq.1) then
                  mdot_sum_denom=level_mdot_data(2*(icolor-1)+1)
                  mdot_sum_numerator=level_mdot_data(2*(icolor-1)+2)
                  if (mdot_sum_denom.gt.zero) then
                   mdot_sum_numerator=level_mdot_data(2*(icolor-1)+2)
                   mdot_local(D_DECL(i,j,k))= &
                    mdot_local(D_DECL(i,j,k))- &
                    mdot_sum_numerator*volcell/mdot_sum_denom
 
                   mdot_global(D_DECL(i,j,k))=mdot_global(D_DECL(i,j,k))+ &
                    mdot_local(D_DECL(i,j,k))
                  else if (mdot_sum_denom.eq.zero) then
                   print *,"since F=1 for this cell, we must have "
                   print *,"mdot_sum_denom>0"
                   stop
                  else if ((mdot_sum_denom.eq.zero).and. &
                           (mdot_sum_numerator.eq.zero)) then
                   ! do nothing
                  else
                   print *,"mdot_sum_denom or mdot_sum_numerator invalid"
                   stop
                  endif
                 else
                  print *,"sweep_num invalid"
                  stop
                 endif

                else
                 print *,"pressure_sum invalid pressure_sum=",pressure_sum
                 print *,"blob_cellvol_count=",blob_cellvol_count
                 stop
                endif

               else
                print *,"material_type_lowmach(im_primary) invalid"
                stop
               endif
              else 
               print *,"pressure_sum NaN, pressure_sum=",pressure_sum
               stop
              endif
             else if (blob_cellvol_count.eq.zero) then
              print *,"since vfrac=1>1/2 for this cell, one must have"
              print *,"blob_cellvol_count>0"
              stop
             else
              print *,"blob_cellvol_count invalid"
              print *,"blob_cellvol_count=",blob_cellvol_count
              stop
             endif
            else
             print *,"(im_primary<>base_type) even though the cell is full"
             stop
            endif
           else
            print *,"vfrac out of range"
            stop
           endif
          else if ((fort_material_type(im_primary).gt.0).and. &
                   (fort_material_type(im_primary).le.MAX_NUM_EOS)) then
           ! do nothing
          else
           print *,"fort_material_type(im_primary) invalid"
           stop
          endif
         else if (is_rigid(im_primary).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(im_primary) invalid"
          stop
         endif
        else
         print *,"vfrac out of range (2)"
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
      end subroutine fort_get_lowmach_divu


      subroutine fort_levelrecolor( &
        color, &
        DIMS(color), &
        xlo,dx, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        domaincolormap, &
        max_colors_level, &
        level,base_level,arrsize) &
      bind(c,name='fort_levelrecolor')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      REAL_T, INTENT(in) :: xlo(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact
 
      INTEGER_T, INTENT(in) ::  DIMDEC(color)
      INTEGER_T, INTENT(in) :: max_colors_level
      INTEGER_T, INTENT(in) :: level,base_level,arrsize
      INTEGER_T, INTENT(in) :: domaincolormap(arrsize)
      REAL_T, INTENT(inout),target :: color(DIMV(color))
      REAL_T, pointer :: color_ptr(D_DECL(:,:,:))

      INTEGER_T  i, j, k, m, icolor

      color_ptr=>color
      call checkbound_array1(fablo,fabhi,color_ptr,1,-1)

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
      end subroutine fort_levelrecolor

      subroutine fort_colorfill( &
       mask,DIMS(mask), &
       typefab,DIMS(typefab), &
       color,DIMS(color), &
       ijk,DIMS(ijk), &
       lo,hi, &
       ipass,number_grids,color_per_grid, &
       gridno,max_colors_grid,typedim) &
      bind(c,name='fort_colorfill')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: typedim
      INTEGER_T, INTENT(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: DIMDEC(mask)
      INTEGER_T, INTENT(in) :: DIMDEC(typefab)
      INTEGER_T, INTENT(in) :: DIMDEC(color)
      INTEGER_T, INTENT(in) :: DIMDEC(ijk)
      INTEGER_T, INTENT(in) :: number_grids,ipass
      INTEGER_T, INTENT(inout) :: color_per_grid(number_grids)
      INTEGER_T, INTENT(in) :: gridno
      INTEGER_T, INTENT(in) :: max_colors_grid

      REAL_T, INTENT(in),target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: typefab(DIMV(typefab))
      REAL_T, pointer :: typefab_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(inout),target :: color(DIMV(color))
      REAL_T, pointer :: color_ptr(D_DECL(:,:,:))
      INTEGER_T, INTENT(inout),target :: ijk(DIMV(ijk),SDIM)
      INTEGER_T, pointer :: ijk_ptr(D_DECL(:,:,:),:)

      INTEGER_T i,j,k,icolor,i1,j1,k1
      INTEGER_T istack,ii,jj,kk
      INTEGER_T iprime,jprime,kprime,base_type,test_type
      INTEGER_T ilocal,jlocal,klocal
      INTEGER_T k1lo,k1hi

      mask_ptr=>mask
      typefab_ptr=>typefab
      color_ptr=>color
      ijk_ptr=>ijk
      call checkbound_array1(lo,hi,mask_ptr,0,-1)
      call checkbound_array1(lo,hi,typefab_ptr,0,-1)
      call checkbound_array1(lo,hi,color_ptr,0,-1)
      call checkbound_array_INTEGER(lo,hi,ijk_ptr,0,-1)

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
      end subroutine fort_colorfill


      subroutine fort_gridrecolor( &
       mask, &
       DIMS(mask), &
       color, &
       DIMS(color), &
       xlo,dx, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       levelcolormap,max_colors_grid,number_grids, &
       arrsize) &
      bind(c,name='fort_gridrecolor')
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      REAL_T, INTENT(in) :: xlo(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: DIMDEC(mask)
      INTEGER_T, INTENT(in) :: DIMDEC(color)
      INTEGER_T, INTENT(in) :: max_colors_grid,number_grids,arrsize
      INTEGER_T, INTENT(in) :: levelcolormap(arrsize)
      REAL_T, INTENT(in),target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(inout),target :: color(DIMV(color))
      REAL_T, pointer :: color_ptr(D_DECL(:,:,:))

      INTEGER_T i,j,k,icolor,testsize

      mask_ptr=>mask
      color_ptr=>color

      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,color_ptr,1,-1)

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
      end subroutine fort_gridrecolor

! components are: color1,type1,color2,type2,color3,type3
      subroutine fort_levelcolorinit( &
        mask,DIMS(mask), &
        color,DIMS(color), &
        xlo,dx, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        level_color, &
        max_colors_level, &
        arrsize,check_corners) &
      bind(c,name='fort_levelcolorinit')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      REAL_T, INTENT(in) :: xlo(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: check_corners
      INTEGER_T, INTENT(in) :: DIMDEC(mask)
      INTEGER_T, INTENT(in) :: DIMDEC(color)
      INTEGER_T, INTENT(in) :: max_colors_level,arrsize
      INTEGER_T, INTENT(out) :: level_color(arrsize,arrsize)
      REAL_T, INTENT(in),target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: color(DIMV(color),6)
      REAL_T, pointer :: color_ptr(D_DECL(:,:,:),:)

      INTEGER_T i,j,k,icolor,jcolor,testsize
      INTEGER_T k1lo,k1hi
      INTEGER_T ii,jj,kk,base_type,near_type
      INTEGER_T nbase,nbase2
      REAL_T mask2

      mask_ptr=>mask
      color_ptr=>color
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      call checkbound_array(fablo,fabhi,color_ptr,1,-1)

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
      end subroutine fort_levelcolorinit

      subroutine fort_avgdowncolor( &
        problo,dxf, &
        bfact_f,bfact, &
        xlo_fine,dx, &
        crse,DIMS(crse), &
        fine,DIMS(fine), &
        typef,DIMS(typef), &
        typec,DIMS(typec), &
        clo,chi) &
      bind(c,name='fort_avgdowncolor')
      use global_utility_module
      IMPLICIT NONE

      REAL_T, INTENT(in) :: problo(SDIM)
      REAL_T, INTENT(in) :: dxf(SDIM)
      INTEGER_T, INTENT(in) :: bfact_f,bfact
      REAL_T, INTENT(in) :: xlo_fine(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      INTEGER_T, INTENT(in) ::  clo(SDIM),chi(SDIM)
      INTEGER_T  growlo(3),growhi(3)
      INTEGER_T  stenlo(3),stenhi(3)
      INTEGER_T, INTENT(in) ::  DIMDEC(crse)
      INTEGER_T, INTENT(in) ::  DIMDEC(fine)
      INTEGER_T, INTENT(in) ::  DIMDEC(typef)
      INTEGER_T, INTENT(in) ::  DIMDEC(typec)
      REAL_T, INTENT(inout),target :: crse(DIMV(crse))
      REAL_T, pointer :: crse_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: fine(DIMV(fine))
      REAL_T, pointer :: fine_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: typef(DIMV(typef))
      REAL_T, pointer :: typef_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: typec(DIMV(typec))
      REAL_T, pointer :: typec_ptr(D_DECL(:,:,:))

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


      crse_ptr=>crse
      fine_ptr=>fine 
      typef_ptr=>typef
      typec_ptr=>typec

      call checkbound_array1(clo,chi,crse_ptr,0,-1)
      call checkbound_array1(clo,chi,typec_ptr,0,-1)
 
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
      end subroutine fort_avgdowncolor


! components are: color1,type1,color2,type2,color3,type3
      subroutine fort_copyfinecoarsecolor( &
        problo,dxf,bfact_f,bfact,xlo_fine,dx, &
        crse,DIMS(crse), &
        fine,DIMS(fine), &
        typef,DIMS(typef), &
        maskf,DIMS(maskf), &
        clo,chi,flo,fhi, &
        zero_diag_flag) &
      bind(c,name='fort_copyfinecoarsecolor')

      use global_utility_module
      use probf90_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: zero_diag_flag
      REAL_T, INTENT(in) :: problo(SDIM)
      REAL_T, INTENT(in) :: dxf(SDIM)
      INTEGER_T, INTENT(in) :: bfact_f,bfact
      REAL_T, INTENT(in) :: xlo_fine(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      INTEGER_T, INTENT(in) :: flo(SDIM),fhi(SDIM)
      INTEGER_T, INTENT(in) :: clo(SDIM),chi(SDIM)
      INTEGER_T  growlo(3),growhi(3)
      INTEGER_T  stenlo(3),stenhi(3)
      INTEGER_T, INTENT(in) :: DIMDEC(crse)
      INTEGER_T, INTENT(in) :: DIMDEC(fine)
      INTEGER_T, INTENT(in) :: DIMDEC(typef)
      INTEGER_T, INTENT(in) :: DIMDEC(maskf)

      REAL_T, INTENT(out),target :: crse(DIMV(crse),6)
      REAL_T, pointer :: crse_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: fine(DIMV(fine))
      REAL_T, pointer :: fine_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: typef(DIMV(typef))
      REAL_T, pointer :: typef_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: maskf(DIMV(maskf))
      REAL_T, pointer :: maskf_ptr(D_DECL(:,:,:))

      INTEGER_T i, j, k, ic, jc, kc,n
      INTEGER_T fine_type,fine_color
      INTEGER_T icrse,jcrse,alreadyhit,crse_color,crse_type
      INTEGER_T ncomp_type
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

      ncomp_type=num_materials
      if (zero_diag_flag.eq.1) then
       ncomp_type=2
      else if (zero_diag_flag.eq.0) then
       ncomp_type=num_materials
      else
       print *,"zero_diag_flag invalid"
       stop
      endif

      crse_ptr=>crse
      fine_ptr=>fine
      typef_ptr=>typef
      maskf_ptr=>maskf

      call checkbound_array(clo,chi,crse_ptr,0,-1)
      call checkbound_array1(flo,fhi,fine_ptr,0,-1)
      call checkbound_array1(flo,fhi,typef_ptr,0,-1)
      call checkbound_array1(flo,fhi,maskf_ptr,0,-1)
 
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
              if ((fine_type.le.0).or.(fine_type.gt.ncomp_type)) then
               print *,"fine_type invalid"
               stop
              else if ((fine_type.ge.1).and.(fine_type.le.ncomp_type)) then
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
              else
               print *,"fine_type bust"
               stop
              endif 

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
      end subroutine fort_copyfinecoarsecolor



! faceden=1/rho
! facecut=A
! icefacecut=1
! mside=mass
! slope: vof,ref centroid,order,slope,intercept  x num_materials
! vof: piecewise constant interp at coarse/fine borders
! mask=1 at interior cells, physical border cells, and fine-fine border cells.
! mask=0 at coarse-fine border cells.
! mask has 3 ghost cells.

      subroutine fort_init_physics_vars( &
       caller_string, &
       caller_string_len, &
       tid, &
       FD_curv_interp, &
       curv_min, &
       curv_max, &
       isweep, &
       nrefine_vof, &
       denconst_interface_added, &
       visc_interface, &
       heatvisc_interface, &
       speciesvisc_interface, &
       freezing_model, &
       distribute_from_target, &
       solidheat_flag, &
       microlayer_size, & ! 1..num_materials
       microlayer_substrate, & ! 1..num_materials
       microlayer_temperature_substrate, & ! 1..num_materials
       spec_material_id_AMBIENT, & ! 1..num_species_var+1
       mass_fraction_id, &
       cavitation_vapor_density, &
       constant_density_all_time, &
       time, &
       dt, &
       project_option, &
       problo,probhi, &
       visc_coef, &
       xlo,dx, &
       maskcov,DIMS(maskcov), &
       masknbr,DIMS(masknbr), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       curv,DIMS(curv), &
       slope,DIMS(slope), &
       denstate, &
       DIMS(denstate), &
       mom_den, &
       DIMS(mom_den), &
       viscstate,DIMS(viscstate), & ! 3 * num_materials components
       conductstate,DIMS(conductstate), & ! num_materials components
       solxfab,DIMS(solxfab), &
       solyfab,DIMS(solyfab), &
       solzfab,DIMS(solzfab), &
        ! voltotal/DeDT_total = 1/(rho cv)
       cenDeDT, &
       DIMS(cenDeDT), & 
        ! voltotal/mass_total = (1/rho)
       cenden, &
       DIMS(cenden), &   
       cenvof,DIMS(cenvof), &   
       cenvisc,DIMS(cenvisc), &
       vol,DIMS(vol), &
       levelPC,DIMS(levelPC), &
       vofC,DIMS(vofC), &
       vofF,DIMS(vofF), & !vofF: see fort_build_semirefinevof, tessellate==3
       massF,DIMS(massF), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       presbc_arr, &
       velbc, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       num_curv, & !num_interfaces * CURVCOMP_NCOMP
       level, &
       finest_level) &
      bind(c,name='fort_init_physics_vars')

      use ISO_C_BINDING, ONLY: C_CHAR,C_INT

      use global_utility_module
      use probf90_module
      use godunov_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      CHARACTER(KIND=C_CHAR), INTENT(in) :: caller_string(*)
      INTEGER(C_INT), INTENT(in), VALUE :: caller_string_len
      CHARACTER(:), ALLOCATABLE :: fort_caller_string
      INTEGER_T :: fort_caller_string_len
      CHARACTER(len=255) :: pattern_string
      INTEGER_T :: pattern_string_len

      INTEGER_T, INTENT(in) :: tid
      INTEGER_T, INTENT(in) :: nparts
      INTEGER_T, INTENT(in) :: nparts_def
      INTEGER_T, INTENT(in) :: im_solid_map(nparts_def)
      REAL_T :: curv_min
      REAL_T :: curv_max
      INTEGER_T, INTENT(in) :: isweep
      INTEGER_T, INTENT(in) :: nrefine_vof
      INTEGER_T, INTENT(in) :: level,finest_level

      INTEGER_T, INTENT(in) :: solidheat_flag
      REAL_T, INTENT(in) :: microlayer_size(num_materials)
      INTEGER_T, INTENT(in) :: microlayer_substrate(num_materials)
      REAL_T, INTENT(in) :: microlayer_temperature_substrate(num_materials)

      INTEGER_T, INTENT(in) :: num_curv
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(in) :: dt
      INTEGER_T, INTENT(in) :: project_option
      REAL_T, INTENT(in) :: problo(SDIM),probhi(SDIM)
      REAL_T, INTENT(in) :: visc_coef
      INTEGER_T, INTENT(in) :: freezing_model(2*num_interfaces)
      INTEGER_T, INTENT(in) :: distribute_from_target(2*num_interfaces)
      INTEGER_T :: veldir
      INTEGER_T, INTENT(in) :: constant_density_all_time(num_materials)
      INTEGER_T, INTENT(in) :: spec_material_id_AMBIENT(num_species_var+1)
      INTEGER_T, INTENT(in) :: mass_fraction_id(2*num_interfaces)
      REAL_T, INTENT(in) :: cavitation_vapor_density(num_materials)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: presbc_arr(SDIM,2)
      INTEGER_T, INTENT(in) :: velbc(SDIM,2,SDIM)
      INTEGER_T, INTENT(in) :: DIMDEC(maskcov)
      INTEGER_T, INTENT(in) :: DIMDEC(masknbr)
      INTEGER_T, INTENT(in) :: DIMDEC(xface)
      INTEGER_T, INTENT(in) :: DIMDEC(yface)
      INTEGER_T, INTENT(in) :: DIMDEC(zface)
      INTEGER_T, INTENT(in) :: DIMDEC(curv)
      INTEGER_T, INTENT(in) :: DIMDEC(slope)
      INTEGER_T, INTENT(in) :: DIMDEC(denstate)
      INTEGER_T, INTENT(in) :: DIMDEC(mom_den)
      INTEGER_T, INTENT(in) :: DIMDEC(viscstate)
      INTEGER_T, INTENT(in) :: DIMDEC(conductstate)
      INTEGER_T, INTENT(in) :: DIMDEC(solxfab)
      INTEGER_T, INTENT(in) :: DIMDEC(solyfab)
      INTEGER_T, INTENT(in) :: DIMDEC(solzfab)
      INTEGER_T, INTENT(in) :: DIMDEC(cenDeDT)
      INTEGER_T, INTENT(in) :: DIMDEC(cenden)
      INTEGER_T, INTENT(in) :: DIMDEC(cenvof)
      INTEGER_T, INTENT(in) :: DIMDEC(vol)
      INTEGER_T, INTENT(in) :: DIMDEC(levelPC)
      INTEGER_T, INTENT(in) :: DIMDEC(cenvisc)
      INTEGER_T, INTENT(in) :: DIMDEC(vofC)
      INTEGER_T, INTENT(in) :: DIMDEC(vofF)
      INTEGER_T, INTENT(in) :: DIMDEC(massF)
      REAL_T, INTENT(in), target :: maskcov(DIMV(maskcov))
      REAL_T, pointer :: maskcov_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: masknbr(DIMV(masknbr),4)
      REAL_T, pointer :: masknbr_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out), target :: xface(DIMV(xface),FACECOMP_NCOMP)
      REAL_T, pointer :: xface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out), target :: yface(DIMV(yface),FACECOMP_NCOMP)
      REAL_T, pointer :: yface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out), target :: zface(DIMV(zface),FACECOMP_NCOMP)
      REAL_T, pointer :: zface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: curv(DIMV(curv),num_curv) 
      REAL_T, pointer :: curv_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: slope(DIMV(slope),num_materials*ngeom_recon) 
      REAL_T, pointer :: slope_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: &
              denstate(DIMV(denstate),num_materials*num_state_material) 
      REAL_T, pointer :: denstate_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: mom_den(DIMV(mom_den),num_materials) 
      REAL_T, pointer :: mom_den_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: viscstate(DIMV(viscstate),3*num_materials) 
      REAL_T, pointer :: viscstate_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: conductstate(DIMV(viscstate),num_materials) 
      REAL_T, pointer :: conductstate_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: solxfab(DIMV(solxfab),nparts_def*SDIM) 
      REAL_T, pointer :: solxfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: solyfab(DIMV(solyfab),nparts_def*SDIM) 
      REAL_T, pointer :: solyfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: solzfab(DIMV(solzfab),nparts_def*SDIM) 
      REAL_T, pointer :: solzfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out), target :: cenDeDT(DIMV(cenDeDT))
      REAL_T, pointer :: cenDeDT_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(out), target :: cenden(DIMV(cenden))
      REAL_T, pointer :: cenden_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(out), target :: cenvof(DIMV(cenvof),num_materials)  
      REAL_T, pointer :: cenvof_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out), target :: cenvisc(DIMV(cenvisc))
      REAL_T, pointer :: cenvisc_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: vol(DIMV(vol))
      REAL_T, pointer :: vol_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: levelPC(DIMV(levelPC),num_materials)
      REAL_T, pointer :: levelPC_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: vofC(DIMV(vofC),num_materials)
      REAL_T, pointer :: vofC_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: vofF(DIMV(vofF),nrefine_vof)
      REAL_T, pointer :: vofF_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: massF(DIMV(massF),nrefine_vof)
      REAL_T, pointer :: massF_ptr(D_DECL(:,:,:),:)

      REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)

      REAL_T, INTENT(in) :: denconst_interface_added(num_interfaces)
      REAL_T, INTENT(in) :: visc_interface(num_interfaces)
      REAL_T, INTENT(in) :: heatvisc_interface(num_interfaces)
      REAL_T, INTENT(in) :: &
          speciesvisc_interface(num_interfaces*num_species_var)

      INTEGER_T im1,jm1,km1
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
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

      REAL_T visc_total
      REAL_T DeDT,DeDT_total

      REAL_T volmat(num_materials)
      REAL_T voltotal
      REAL_T voltotal_prescribed
      REAL_T voldepart
      REAL_T mass_total
      REAL_T delta_mass
      REAL_T den
      REAL_T density_for_mass_fraction_diffusion
      REAL_T TEMPERATURE

      INTEGER_T LS_consistent
      INTEGER_T LS_consistent_tension

      REAL_T LSplus(num_materials)
      REAL_T LSminus(num_materials)
      INTEGER_T im_main,im_main_opp,im_opp
      INTEGER_T im_left_main,im_right_main
      INTEGER_T ireverse
      INTEGER_T iten_FFACE
      INTEGER_T iten_main
      INTEGER_T iten_majority
      INTEGER_T local_iten
      INTEGER_T zeroradius_flag
      INTEGER_T iten_tension
      INTEGER_T iten_micro
      INTEGER_T im_tension,im_opp_tension
      INTEGER_T im_left_tension,im_right_tension
      INTEGER_T im_fluid_micro,im_solid_micro

      REAL_T gradh,gradh_tension
      REAL_T sign_test
      REAL_T one_over_mu
      REAL_T localvisc(num_materials)
      INTEGER_T null_viscosity
       !0=>interior wall 1..num_materials=>embedded wall 
       !num_materials+1=>domain wall 
      INTEGER_T wall_flag_face 
      REAL_T wallVOF_face
      INTEGER_T imattype
      INTEGER_T dencomp,tempcomp
      REAL_T facevisc_local
      REAL_T faceheat_local
      REAL_T facespecies_local(num_species_var+1)
      REAL_T theta,visc1,visc2,heat1,heat2
      REAL_T spec1(num_species_var+1)
      REAL_T spec2(num_species_var+1)
      REAL_T spec_test
      REAL_T localvisc_plus(num_materials)
      REAL_T localvisc_minus(num_materials)
      REAL_T localheatvisc_plus(num_materials)
      REAL_T localheatvisc_minus(num_materials)
      INTEGER_T implus_majority,imminus_majority
      INTEGER_T im_secondary

      REAL_T local_face(FACECOMP_NCOMP)
      REAL_T local_volumes(2,num_materials)
      REAL_T local_masses(2,num_materials)
      INTEGER_T igridlo(3),igridhi(3)

      INTEGER_T dir_refine
      INTEGER_T noslip_wall,velbclo,velbchi
      INTEGER_T imspec,is_zero_visc
      REAL_T solid_velocity

      REAL_T LSIDE(2)
      REAL_T LSIDE_tension(2)
      REAL_T LSIDE_MAT(num_materials)
      REAL_T LSIDE_tension_MAT(num_materials)

      REAL_T DXMAXLS
      REAL_T FFACE(num_materials)
      INTEGER_T irefine
      INTEGER_T solid_present_flag
      REAL_T wtL,wtR,wtsum
      INTEGER_T, parameter :: nhalf=3
      REAL_T xstenMAC(-nhalf:nhalf,SDIM)
      REAL_T xstenMAC_center(SDIM)
      REAL_T xsten(-nhalf:nhalf,SDIM)

      REAL_T local_plus
      REAL_T local_minus
      REAL_T xclamped_minus(SDIM)
      REAL_T xclamped_plus(SDIM)
      REAL_T LS_clamped_minus
      REAL_T LS_clamped_plus
      REAL_T vel_clamped_minus(SDIM)
      REAL_T vel_clamped_plus(SDIM)
      REAL_T temperature_clamped_minus
      REAL_T temperature_clamped_plus
      INTEGER_T prescribed_flag
      INTEGER_T prescribed_exists_flag
      INTEGER_T is_clamped_face
      REAL_T vel_clamped_face(SDIM)

      INTEGER_T icurv,icurv_ofs
      INTEGER_T curv_interp_flag
      INTEGER_T dirL,sideL,im3L,orientL
      INTEGER_T dirR,sideR,im3R,orientR
      REAL_T curvL(CURVCOMP_NCOMP)
      REAL_T curvR(CURVCOMP_NCOMP)
      INTEGER_T local_maskL,local_maskR,local_masknbr
      INTEGER_T covered_face
      INTEGER_T coarse_fine_face
      INTEGER_T FD_curv_interp
      REAL_T mofdata(num_materials*ngeom_recon)
      INTEGER_T micro_table(num_materials,num_materials)
      REAL_T predict_face_afrac_solid !=0 if clamped or is_rigid
      REAL_T predict_face_afrac_prescribed !=0 if clamped or prescribed
      INTEGER_T is_solid_face
      INTEGER_T is_prescribed_face
      INTEGER_T im_solid
      INTEGER_T im_prescribed
      INTEGER_T im_solid_valid
      INTEGER_T im_prescribed_valid
      INTEGER_T partid_solid
      INTEGER_T partid_prescribed
      INTEGER_T im_prescribed_primary
      INTEGER_T ispec
      REAL_T massfrac_parm(num_species_var+1)
      INTEGER_T local_tessellate

      REAL_T local_cenvisc
      REAL_T local_cenden

      REAL_T, PARAMETER :: VISCINVTOL=1.0D-8
      INTEGER_T test_for_quasi_static

      REAL_T denconst_interface_added_max
      REAL_T rad_added_mass

! fort_init_physics_vars code starts here:

      allocate(CHARACTER(caller_string_len) :: fort_caller_string)
      do i=1,caller_string_len
       fort_caller_string(i:i)=caller_string(i)
      enddo
      fort_caller_string_len=caller_string_len

      pattern_string='static_surface_tension_advection'
      pattern_string_len=32
      test_for_quasi_static=fort_pattern_test( &
        fort_caller_string,fort_caller_string_len, &
        pattern_string,pattern_string_len)

      maskcov_ptr=>maskcov
      masknbr_ptr=>masknbr
      xface_ptr=>xface
      yface_ptr=>yface
      zface_ptr=>zface
      curv_ptr=>curv
      slope_ptr=>slope
      denstate_ptr=>denstate
      mom_den_ptr=>mom_den
      viscstate_ptr=>viscstate
      conductstate_ptr=>conductstate
      solxfab_ptr=>solxfab
      solyfab_ptr=>solyfab
      solzfab_ptr=>solzfab
      cenDeDT_ptr=>cenDeDT

      cenden_ptr=>cenden

      cenvof_ptr=>cenvof
      cenvisc_ptr=>cenvisc
      vol_ptr=>vol
      levelPC_ptr=>levelPC
      vofC_ptr=>vofC
      vofF_ptr=>vofF
      massF_ptr=>massF

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if (nrefine_vof.ne.2*num_materials*SDIM) then
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

      if (num_curv.ne.num_interfaces*CURVCOMP_NCOMP) then
       print *,"num_curv invalid"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_init_physics_vars"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_init_physics_vars"
       stop
      endif

      dir_refine=SDIM

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((solidheat_flag.lt.0).or. &
          (solidheat_flag.gt.2)) then
       print *,"solidheat_flag invalid"
       stop
      endif

      if ((project_option.eq.SOLVETYPE_PRES).or. &
          (project_option.eq.SOLVETYPE_INITPROJ)) then
       ! do nothing
      else
       print *,"project_option invalid fort_init_physics_vars"
       stop
      endif

      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid fort_init_physics_vars "
       print *,"ngeom_recon= ",ngeom_recon
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid fort_init_physics_vars "
       print *,"ngeom_raw= ",ngeom_raw
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: fort_init_physics_vars

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,masknbr_ptr,1,-1)

      call checkbound_array(fablo,fabhi,xface_ptr,0,0)
      call checkbound_array(fablo,fabhi,yface_ptr,0,1)
      call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1)

      call checkbound_array1(fablo,fabhi,cenDeDT_ptr,1,-1)

      call checkbound_array1(fablo,fabhi,cenden_ptr,1,-1)

      call checkbound_array(fablo,fabhi,cenvof_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,cenvisc_ptr,1,-1)

      call checkbound_array(fablo,fabhi,slope_ptr,1,-1)
      call checkbound_array(fablo,fabhi,curv_ptr,1,-1)
      call checkbound_array(fablo,fabhi,denstate_ptr,1,-1)
      call checkbound_array(fablo,fabhi,mom_den_ptr,1,-1)
      call checkbound_array(fablo,fabhi,viscstate_ptr,1,-1)
      call checkbound_array(fablo,fabhi,conductstate_ptr,1,-1)

      call checkbound_array(fablo,fabhi,solxfab_ptr,0,0)
      call checkbound_array(fablo,fabhi,solyfab_ptr,0,1)
      call checkbound_array(fablo,fabhi,solzfab_ptr,0,SDIM-1)

      call checkbound_array1(fablo,fabhi,vol_ptr,1,-1)
      call checkbound_array(fablo,fabhi,levelPC_ptr,2,-1) 
      call checkbound_array(fablo,fabhi,vofC_ptr,1,-1)
      call checkbound_array(fablo,fabhi,vofF_ptr,1,-1)
      call checkbound_array(fablo,fabhi,massF_ptr,1,-1)

      call get_dxmaxLS(dx,bfact,DXMAXLS)

      if (DXMAXLS.gt.zero) then
       ! do nothing
      else
       print *,"DXMAXLS invalid: fort_init_physics_vars"
       stop
      endif

      rad_added_mass=1.5d0*DXMAXLS

      do im=1,2*num_interfaces

       if (is_valid_freezing_modelF(freezing_model(im)).eq.1) then
        ! do nothing
       else
        print *,"freezing_model invalid fort_init_physics_vars"
        stop
       endif

       if ((distribute_from_target(im).lt.0).or. &
           (distribute_from_target(im).gt.1)) then
        print *,"distribute_from_target invalid"
        stop
       endif

      enddo !im=1,2*num_interfaces

      do im=1,num_interfaces

       if ((denconst_interface_added(im).ge.zero).and. &
           (visc_interface(im).ge.zero).and. &
           (heatvisc_interface(im).ge.zero)) then
        ! do nothing
       else
        print *,"den,visc, or heat interface coeff wrong"
        stop
       endif
       do imspec=1,num_species_var
        if (speciesvisc_interface(num_interfaces*(imspec-1)+im).ge.zero) then
         !do nothing
        else
         print *,"species interface coeff wrong"
         stop
        endif
       enddo

      enddo ! im=1..num_interfaces

      do im=1,num_materials

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).gt.0).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
        ! do nothing
       else 
        print *,"material_type invalid"
        stop
       endif

       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"density must be positive fort_init_physics_vars"
        print *,"im,denconst ",im,fort_denconst(im)
        stop
       endif
       if (fort_tempconst(im).gt.zero) then
        ! do nothing
       else
        print *,"fort_init_physics_vars:temperature must be positive"
        print *,"im,fort_tempconst : ",im,fort_tempconst(im)
        stop
       endif
       if (fort_energyconst(im).gt.zero) then
        ! do nothing
       else
        print *,"energy must be positive in fort_init_physics_vars"
        print *,"im= ",im
        print *,"fort_energyconst(im)= ",fort_energyconst(im)
        stop
       endif
        ! sanity check: the real viscosity coefficient(s) are derived from
        ! viscstate(D_DECL(:,:,:),3*num_materials)
       if (get_user_viscconst(im,fort_denconst(im),fort_tempconst(im)).ge. &
           zero) then
        ! do nothing
       else
        print *,"viscosity cannot be negative"
        stop
       endif

      enddo ! im=1..num_materials  (checking parameters)

       ! create look-up table for thin liquid layer model 4
      do im_fluid_micro=1,num_materials
      do im_solid_micro=1,num_materials
       micro_table(im_fluid_micro,im_solid_micro)=0
      enddo 
      enddo 

      if (solidheat_flag.eq.0) then ! diffuse in solid
       do im=1,num_materials-1
        if (is_rigid(im).eq.0) then
         do im_opp=im+1,num_materials
          if (is_rigid(im_opp).eq.0) then
           call get_iten(im,im_opp,iten_main)

           do ireverse=0,1
            local_iten=iten_main+ireverse*num_interfaces
            if (get_user_latent_heat(local_iten,293.0d0,1).ne.zero) then
             if (freezing_model(local_iten).eq.0) then
              im_solid_micro=microlayer_substrate(im)
              if ((im_solid_micro.ge.1).and. &
                  (im_solid_micro.le.num_materials)) then

               if (is_rigid(im_solid_micro).eq.1) then
                ! do nothing
               else
                print *,"expecting is_rigid(im_solid_micro)==1"
                stop
               endif 

               micro_table(im_opp,im_solid_micro)=iten_main

              else if (im_solid_micro.eq.0) then
               ! do nothing
              else
               print *,"im_solid_micro invalid"
               stop
              endif 
              im_solid_micro=microlayer_substrate(im_opp)
              if ((im_solid_micro.ge.1).and. &
                  (im_solid_micro.le.num_materials)) then

               if (is_rigid(im_solid_micro).eq.1) then
                ! do nothing
               else
                print *,"expecting is_rigid(im_solid_micro)==1"
                stop
               endif 

               micro_table(im,im_solid_micro)=iten_main
              else if (im_solid_micro.eq.0) then
               ! do nothing
              else
               print *,"im_solid_micro invalid"
               stop
              endif 
             endif
            else if (get_user_latent_heat(local_iten,293.0d0,1).eq.zero) then
             ! do nothing
            else
             print *,"get_user_latent_heat invalid"
             stop
            endif
           enddo ! ireverse=0..1

          else if (is_rigid(im_opp).eq.1) then
           ! do nothing
          else
           print *,"is_rigid invalid LEVELSET_3D.F90"
           stop
          endif
         enddo !im_opp=im+1,num_materials
        else if (is_rigid(im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid(im) invalid"
         stop
        endif
       enddo ! do im=1,num_materials-1
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

          ! first init xface,yface,zface 
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

          ! in: fort_init_physics_vars
         local_face(FACECOMP_ICEMASK+1)=one
         local_face(FACECOMP_CURV+1)=zero
         local_face(FACECOMP_FACEVEL+1)=zero

          ! veldir=0..sdim-1
         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,veldir)
         do dir2=1,SDIM
          xstenMAC_center(dir2)=xstenMAC(0,dir2)
         enddo

         wtL=(xstenMAC(0,veldir+1)-xstenMAC(-1,veldir+1))
         wtR=(xstenMAC(1,veldir+1)-xstenMAC(0,veldir+1))
         wtsum=wtL+wtR
         if ((wtL.gt.zero).and. &
             (wtR.gt.zero).and. &
             (wtsum.gt.zero)) then
          ! do nothing
         else
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
          call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)

          if (iside.eq.0) then
           do dir2=1,SDIM
            xclamped_minus(dir2)=xsten(0,dir2)
           enddo
          else if (iside.eq.1) then
           do dir2=1,SDIM
            xclamped_plus(dir2)=xsten(0,dir2)
           enddo
          else
           print *,"iside invalid"
           stop
          endif

          do im=1,num_materials
           irefine=veldir*2*num_materials+(1-iside)*num_materials+im
           local_volumes(iside+1,im)=vofF(D_DECL(icell,jcell,kcell),irefine)
           local_masses(iside+1,im)=massF(D_DECL(icell,jcell,kcell),irefine)
          enddo ! im=1..num_materials
         enddo ! iside=0..1

         do iside=0,1
         do im=1,num_materials
          local_face(FACECOMP_VOFFACE+2*(im-1)+iside+1)= &
            local_volumes(iside+1,im)
          local_face(FACECOMP_MASSFACE+2*(im-1)+iside+1)= &
            local_masses(iside+1,im)
         enddo ! im=1..num_materials
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

         if (levelrz.eq.COORDSYS_CARTESIAN) then
          im1=i-ii
         else if (levelrz.eq.COORDSYS_RZ) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          if (xstenMAC(-1,1).le.VOFTOL*dx(1)) then
           im1=i
          else
           im1=i-ii
          endif
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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

         do im=1,num_materials
          LSplus(im)=levelPC(D_DECL(i,j,k),im)
          LSminus(im)=levelPC(D_DECL(im1,jm1,km1),im)
          LSIDE_MAT(im)=half*(LSplus(im)+LSminus(im))
         enddo

          ! checks rigid and non-rigid materials.
         call get_primary_material(LSplus,implus_majority)
         call get_primary_material(LSminus,imminus_majority)

         if ((implus_majority.lt.1).or. &
             (implus_majority.gt.num_materials).or. &
             (imminus_majority.lt.1).or. &
             (imminus_majority.gt.num_materials)) then
          print *,"implus_majority or imminus_majority invalid"
          stop
         endif

         if ((time.ge.zero).and. &
             (dt.ge.zero)) then
          ! do nothing
         else
          print *,"fort_init_physics_vars:"
          print *,"time or dt invalid: time,dt ",time,dt
          stop
         endif

          ! LS>0 if clamped
         call SUB_clamped_LS(xclamped_minus,time,LS_clamped_minus, &
           vel_clamped_minus,temperature_clamped_minus,prescribed_flag,dx)
         call SUB_clamped_LS(xclamped_plus,time,LS_clamped_plus, &
           vel_clamped_plus,temperature_clamped_plus,prescribed_flag,dx)
         if ((LS_clamped_minus.ge.zero).or. &
             (LS_clamped_plus.ge.zero)) then
          is_clamped_face=1
          do dir2=1,SDIM
           if (LS_clamped_minus.lt.zero) then
            vel_clamped_face(dir2)=vel_clamped_plus(dir2)
            is_clamped_face=2
           else if (LS_clamped_plus.lt.zero) then
            vel_clamped_face(dir2)=vel_clamped_minus(dir2)
            is_clamped_face=3
           else
            vel_clamped_face(dir2)=half*(vel_clamped_plus(dir2)+ &
              vel_clamped_minus(dir2))
           endif
          enddo
         else if ((LS_clamped_minus.lt.zero).and. &
                  (LS_clamped_plus.lt.zero)) then
          is_clamped_face=0
         else
          print *,"LS_clamped plus or minus is NaN"
          stop
         endif

          ! calls fort_is_rigid_base(FSI_flag(im),im) (GLOBALUTIL.F90)
          !  FSI_PRESCRIBED_PROBF90,
          !  FSI_PRESCRIBED_NODES,
          !  FSI_SHOELE_PRESVEL,
          !  FSI_SHOELE_VELVEL.
          !
          ! Remark: 
          ! local_face(FACECOMP_ICEFACECUT+1) is initialized in
          ! in GODUNOV_3D.F90: fort_init_icemask
          !
          ! local_face(FACECOMP_FACECUT+1) is initialized in 
          ! LEVELSET_3D.F: fort_init_physics_vars (here)
          !
         if ((is_rigid(implus_majority).eq.1).or. &
             (is_rigid(imminus_majority).eq.1).or. &
             (is_clamped_face.ge.1)) then
          predict_face_afrac_solid=zero
         else if ((is_rigid(implus_majority).eq.0).and. &
                  (is_rigid(imminus_majority).eq.0).and. &
                  (is_clamped_face.eq.0)) then
          predict_face_afrac_solid=one
         else
          print *,"implus_majority, imminus_majority,or is_clamped.. invalid"
          stop
         endif

          ! calls is_rigid(im) (GLOBALUTIL.F90)
          !  is_rigid(im)==1 AND (! FSI_SHOELE_VELVEL)
         if ((is_prescribed(implus_majority).eq.1).or. &
             (is_prescribed(imminus_majority).eq.1).or. &
             (is_clamped_face.ge.1)) then

          if (predict_face_afrac_solid.eq.zero) then
           predict_face_afrac_prescribed=zero
          else
           print *,"is_prescribed==TRUE => is_rigid==TRUE"
           stop
          endif

         else if ((is_prescribed(implus_majority).eq.0).and. &
                  (is_prescribed(imminus_majority).eq.0).and. &
                  (is_clamped_face.eq.0)) then
          predict_face_afrac_prescribed=one
         else
          print *,"implus_majority or imminus_majority invalid"
          stop
         endif

         if (is_clamped_face.ge.1) then
          is_prescribed_face=1
          is_solid_face=1
          im_solid=0
          im_prescribed=0
          im_solid_valid=1
          im_prescribed_valid=1
          solid_velocity=vel_clamped_face(veldir+1)
         else if (is_clamped_face.eq.0) then
          ! fixed_face is declared in: PROB.F90
          call fixed_face( &
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
              if (veldir.eq.0) then
               solid_velocity=solxfab(D_DECL(i,j,k),dir2)
              else if (veldir.eq.1) then
               solid_velocity=solyfab(D_DECL(i,j,k),dir2)
              else if ((veldir.eq.SDIM-1).and.(SDIM.eq.3)) then
               solid_velocity=solzfab(D_DECL(i,j,k),dir2)
              else
               print *,"veldir invalid"
               stop
              endif
             else
              print *,"im_solid_map(partid_prescribed+1)+1.ne.im_prescribed"
              stop
             endif
            else 
             print *,"im_prescribed_valid bust"
             stop
            endif 
           else
            print *,"in fort_init_physics_vars"
            print *,"is_solid_face invalid(2) is_solid_face= ",is_solid_face
            print *,"tid=",tid
            print *,"isweep=",isweep
            print *,"FACECOMP_NCOMP=",FACECOMP_NCOMP
            print *,"solidheat_flag=",solidheat_flag
            print *,"time=",time
            print *,"dt=",dt
            print *,"project_option=",project_option
            print *,"num_interfaces=",num_interfaces
            print *,"num_materials=",num_materials
            print *,"nparts=",nparts
            print *,"nparts_def=",nparts_def
            print *,"num_curv=",num_curv
            print *,"level=",level
            print *,"finest_level=",finest_level
            do im=1,num_materials
             print *,"im,FSI_flag(im) ",im,FSI_flag(im)
            enddo
            print *,"is_solid_face ",is_solid_face
            print *,"is_prescribed_face ",is_prescribed_face
            print *,"partid_solid ",partid_solid
            print *,"partid_prescribed ",partid_prescribed
            print *,"im_solid_valid ",im_solid_valid
            print *,"im_prescribed_valid ",im_prescribed_valid
            print *,"im_solid ",im_solid
            print *,"im_prescribed ",im_prescribed
            print *,"predict_face_afrac_solid ",predict_face_afrac_solid
            print *,"predict_face_afrac_prescribed ", &
                   predict_face_afrac_prescribed

            stop
           endif
          else if (is_prescribed_face.eq.0) then
           solid_velocity=zero
          else
           print *,"is_prescribed_face invalid"
           stop
          endif
         else
          print *,"is_clamped_face invalid"
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

          if (presbclo.eq.INT_DIR) then
           ! do nothing
          else if (presbclo.eq.EXT_DIR) then
           ! do nothing
          else if ((presbclo.eq.REFLECT_EVEN).or. &
                   (presbclo.eq.FOEXTRAP)) then
           mask_boundary=1
          else
           print *,"presbclo invalid"
           stop
          endif

          if (velbclo.eq.REFLECT_ODD) then
           local_face(FACECOMP_FACEVEL+1)=zero
          else if (velbclo.eq.EXT_DIR) then
           iside=1
           call velbc_override(time,veldir+1,iside,veldir+1, &
            local_face(FACECOMP_FACEVEL+1), &
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

          if (presbchi.eq.INT_DIR) then
           ! do nothing
          else if (presbchi.eq.EXT_DIR) then
           ! do nothing
          else if ((presbchi.eq.REFLECT_EVEN).or. &
                   (presbchi.eq.FOEXTRAP)) then
           mask_boundary=2
          else
           print *,"presbchi invalid"
           stop
          endif

          if (velbchi.eq.REFLECT_ODD) then
           local_face(FACECOMP_FACEVEL+1)=zero
          else if (velbchi.eq.EXT_DIR) then
           iside=2
           call velbc_override(time,veldir+1,iside,veldir+1, &
            local_face(FACECOMP_FACEVEL+1), &
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
           local_face(FACECOMP_FACEVEL+1)=solid_velocity
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

          ! gradh>0 => im_main material dominates right and im_main_opp left.
          ! gradh<0 => im_main material dominates left and im_main_opp right.
          ! im_main<im_main_opp
         if (is_solid_face.eq.1) then
          gradh=zero
          gradh_tension=zero
          iten_main=0
         else if (is_solid_face.eq.0) then

          ! "merge_levelset" is NOT called inside of "fluid_interface"
          ! fluid_interface is declared in: PROB.F90
          call fluid_interface(LSminus,LSplus,gradh, &
            im_main_opp,im_main, &
            im_left_main,im_right_main)

          if (gradh.ne.zero) then
           if (im_main.ge.im_main_opp) then
            print *,"fluid_interface bust"
            stop
           endif
           call get_iten(im_main,im_main_opp,iten_main)
          else if (gradh.eq.zero) then
           iten_main=0
          else
           print *,"gradh is NaN"
           stop
          endif

          if (covered_face.eq.2) then !maskL=maskR=0
           gradh_tension=zero
          else if ((covered_face.eq.0).or. & !maskL=maskR=1
                   (covered_face.eq.1)) then !maskL=1 or maskR=1
             ! fluid_interface_tension is declared in: PROB.F90
             ! "merge_levelset" is called inside of "fluid_interface_tension"
           call fluid_interface_tension( &
            xstenMAC_center,time, &
            LSminus,LSplus,gradh_tension, &
            im_opp_tension,im_tension, & !INTENT(out)
            im_left_tension,im_right_tension, & !INTENT(out)
            test_for_quasi_static)
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
          !viscstate is initialized in fort_derviscosity
         do im=1,num_materials
          if (test_for_quasi_static.eq.1) then
           localvisc_plus(im)=zero
           localvisc_minus(im)=zero
          else if (test_for_quasi_static.eq.0) then
           localvisc_plus(im)=viscstate(D_DECL(i,j,k),im)
           localvisc_minus(im)=viscstate(D_DECL(im1,jm1,km1),im)
          else
           print *,"test_for_quasi_static invalid"
           stop
          endif
          localheatvisc_plus(im)=conductstate(D_DECL(i,j,k),im)
          localheatvisc_minus(im)=conductstate(D_DECL(im1,jm1,km1),im)
         enddo

         if (gradh.ne.zero) then
          if (im_main.ge.im_main_opp) then
           print *,"fluid_interface bust"
           stop
          endif
          if ((iten_main.ge.1).and.(iten_main.le.num_interfaces)) then
           ! do nothing
          else
           print *,"iten_main invalid"
           stop
          endif

          do im=1,num_materials
           LSIDE_MAT(im)=levelPC(D_DECL(im1,jm1,km1),im)
          enddo
          call get_LS_extend(LSIDE_MAT,iten_main,LSIDE(1))

          do im=1,num_materials
           LSIDE_MAT(im)=levelPC(D_DECL(i,j,k),im)
          enddo
          call get_LS_extend(LSIDE_MAT,iten_main,LSIDE(2))

          if (LSIDE(1)*LSIDE(2).le.zero) then
           LS_consistent=1
          else if (LSIDE(1)*LSIDE(2).gt.zero) then
           LS_consistent=0
          else
           print *,"LSIDE is NaN"
           stop
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
          call get_iten(im_tension,im_opp_tension,iten_tension)

          do im=1,num_materials
           LSIDE_tension_MAT(im)=levelPC(D_DECL(im1,jm1,km1),im)
          enddo
          call get_LS_extend(LSIDE_tension_MAT,iten_tension,LSIDE_tension(1))

          do im=1,num_materials
           LSIDE_tension_MAT(im)=levelPC(D_DECL(i,j,k),im)
          enddo
          call get_LS_extend(LSIDE_tension_MAT,iten_tension,LSIDE_tension(2))

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
         if (is_clamped_face.ge.1) then
          solid_present_flag=is_clamped_face
         else if (is_clamped_face.eq.0) then
          if ((is_rigid(implus_majority).eq.1).and. &
              (is_rigid(imminus_majority).eq.1)) then
           solid_present_flag=1
          else if ((is_rigid(implus_majority).eq.1).and. &
                   (is_rigid(imminus_majority).eq.0)) then
           solid_present_flag=2
          else if ((is_rigid(implus_majority).eq.0).and. &
                   (is_rigid(imminus_majority).eq.1)) then
           solid_present_flag=3
          else if ((is_rigid(implus_majority).eq.0).and. &
                   (is_rigid(imminus_majority).eq.0)) then
           solid_present_flag=0
          else
           print *,"implus_majority or imminus_majority invalid"
           stop
          endif
         else 
          print *,"is_clamped_face invalid"
          stop
         endif

         if (dt.eq.zero) then
          ! do nothing
         else if (dt.gt.zero) then
          ! do nothing
         else
          print *,"dt became corrupt; dt=",dt
          stop
         endif

          ! both adjoining cells are solid or clamped cells.
         if (solid_present_flag.eq.1) then

           ! dirichlet cond for velocity at the solid boundaries and 
           ! clamped boundaries.
          facevisc_local=zero  

          local_plus=localheatvisc_plus(implus_majority)
          local_minus=localheatvisc_minus(imminus_majority)
          call geom_avg(local_plus,local_minus,wtR,wtL,faceheat_local)

          do imspec=1,num_species_var
           local_plus= &
            fort_speciesviscconst((imspec-1)*num_materials+implus_majority)
           local_minus= &
            fort_speciesviscconst((imspec-1)*num_materials+imminus_majority)
           call geom_avg(local_plus,local_minus,wtR,wtL, &
                   facespecies_local(imspec))
          enddo !do imspec=1,num_species_var

          if (is_clamped_face.ge.1) then
           ! diffuse temperature and species in the clamped solid
          else if (is_clamped_face.eq.0) then
 
            ! both adjoining cells are solid cells.
           if (solidheat_flag.eq.0) then ! diffuse in solid
            ! do nothing (faceheat_local already defined)
           else if (solidheat_flag.eq.1) then ! dirichlet
            faceheat_local=zero
            do imspec=1,num_species_var
             facespecies_local(imspec)=zero
            enddo
           else if (solidheat_flag.eq.2) then ! neumann
            faceheat_local=zero
            do imspec=1,num_species_var
             facespecies_local(imspec)=zero
            enddo
           else
            print *,"solidheat_flag invalid"
            stop
           endif
          else
           print *,"is_clamped_face invalid"
           stop
          endif

           ! one adjoining cell is a solid cell.
         else if ((solid_present_flag.eq.2).or. &
                  (solid_present_flag.eq.3)) then
         
          call geom_avg(localvisc_plus(implus_majority), &
                  localvisc_minus(imminus_majority), &
                  wtR,wtL,facevisc_local)

          local_plus=localheatvisc_plus(implus_majority)
          local_minus=localheatvisc_minus(imminus_majority)
          call geom_avg(local_plus,local_minus,wtR,wtL,faceheat_local)

          do imspec=1,num_species_var
           local_plus= &
            fort_speciesviscconst((imspec-1)*num_materials+implus_majority)
           local_minus= &
            fort_speciesviscconst((imspec-1)*num_materials+imminus_majority)
           call geom_avg(local_plus,local_minus,wtR,wtL, &
                   facespecies_local(imspec))
          enddo

          if (is_clamped_face.ge.1) then
           ! do nothing special for temperature or mass fraction
          else if (is_clamped_face.eq.0) then

           iten_micro=0
           if (solidheat_flag.eq.0) then ! diffuse in solid
            im_solid_micro=0
            im_fluid_micro=0
            if (is_rigid(implus_majority).eq.1) then
             im_solid_micro=implus_majority
            else if (is_rigid(implus_majority).eq.0) then
             im_fluid_micro=implus_majority
            else
             print *,"is_rigid(implus_majority) invalid"
             stop
            endif
            if (is_rigid(imminus_majority).eq.1) then
             im_solid_micro=imminus_majority
            else if (is_rigid(imminus_majority).eq.0) then
             im_fluid_micro=imminus_majority
            else
             print *,"is_rigid(imminus_majority) invalid"
             stop
            endif
            if ((im_fluid_micro.ge.1).and. &
                (im_fluid_micro.le.num_materials).and. &
                (im_solid_micro.ge.1).and. &
                (im_solid_micro.le.num_materials)) then
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
            else if ((iten_micro.ge.1).and. &
                     (iten_micro.le.num_interfaces)) then
             faceheat_local=zero ! internal TSAT dirichlet bc 
             do imspec=1,num_species_var
              facespecies_local(imspec)=zero
             enddo
            else
             print *,"iten_micro invalid"
             stop
            endif  
           else if (solidheat_flag.eq.1) then
            ! dirichlet temperature bc
            if (is_in_probtype_list().eq.1) then
             call SUB_microcell_heat_coeff(faceheat_local,dx,veldir)
            endif
            do imspec=1,num_species_var
             facespecies_local(imspec)=zero
            enddo
           else if (solidheat_flag.eq.2) then
            faceheat_local=zero  ! neumann temperature bc
            do imspec=1,num_species_var
             facespecies_local(imspec)=zero
            enddo
           else
            print *,"solidheat_flag invalid"
            stop
           endif

            ! check if face coefficient needs to be replaced
            ! prescribed interface coefficient.
           if (implus_majority.ne.imminus_majority) then

            call get_iten(imminus_majority,implus_majority,iten_majority)

            if (heatvisc_interface(iten_majority).eq.zero) then
             ! do nothing
            else if (heatvisc_interface(iten_majority).gt.zero) then

             do ireverse=0,1
              local_iten=iten_majority+ireverse*num_interfaces
              if (get_user_latent_heat(local_iten,293.0d0,1).ne.zero) then
               if ((freezing_model(local_iten).eq.0).or. &
                   (freezing_model(local_iten).eq.5)) then
                print *,"heatvisc_interface invalid"
                stop
               endif 
              else if (get_user_latent_heat(local_iten,293.0d0,1).eq.zero) then
               !do nothing
              else
               print *,"get_user_latent_heat invalid"
               stop
              endif 
             enddo !do ireverse=0,1

             faceheat_local=heatvisc_interface(iten_majority)

            else
             print *,"heatvisc_interface invalid"
             stop
            endif
           else if (implus_majority.eq.imminus_majority) then
            ! do nothing
           else
            print *,"implus_majority or imminus_majority bust"
            stop
           endif 

          else
           print *,"is_clamped_face invalid"
           stop
          endif

           ! neither adjoining cell is a solid cell.
           ! All adjoining cells are either "fluid" cells or
           ! "ice" cells or "FSI_rigid" cells.
         else if (solid_present_flag.eq.0) then

          if (gradh.eq.zero) then
           visc1=localvisc_minus(imminus_majority)
           visc2=localvisc_plus(implus_majority)

           call geom_avg( &
             visc2, &  !intent(in)
             visc1, &  !intent(in)
             wtR,wtL, & !intent(in)
             facevisc_local) !intent(out)

           if ((is_ice_or_FSI_rigid_material(implus_majority).eq.1).and. &
               (is_ice_or_FSI_rigid_material(imminus_majority).eq.1)) then
            facevisc_local=zero
            localvisc_minus(imminus_majority)=zero
            localvisc_plus(implus_majority)=zero
           else if ((is_ice_or_FSI_rigid_material(implus_majority).eq.1).and. &
                    (is_ice_or_FSI_rigid_material(imminus_majority).eq.0)) then
            !do nothing
           else if ((is_ice_or_FSI_rigid_material(implus_majority).eq.0).and. &
                    (is_ice_or_FSI_rigid_material(imminus_majority).eq.1)) then
            !do nothing
           else if ((is_ice_or_FSI_rigid_material(implus_majority).eq.0).and. &
                    (is_ice_or_FSI_rigid_material(imminus_majority).eq.0)) then
            !do nothing
           else
            print *,"implus_majority or imminus_majority invalid"
            stop
           endif

           local_plus=localheatvisc_plus(implus_majority)
           local_minus=localheatvisc_minus(imminus_majority)
           call geom_avg(local_plus,local_minus,wtR,wtL,faceheat_local)

           do imspec=1,num_species_var
            local_plus= &
             fort_speciesviscconst((imspec-1)*num_materials+implus_majority)
            local_minus= &
             fort_speciesviscconst((imspec-1)*num_materials+imminus_majority)
            call geom_avg(local_plus,local_minus,wtR,wtL, &
                   facespecies_local(imspec))
           enddo

          else if (gradh.ne.zero) then

           if ((im_main.gt.num_materials).or. &
               (im_main_opp.gt.num_materials)) then
            print *,"im_main or im_main_opp bust 3"
            stop
           endif
           if (im_main.ge.im_main_opp) then
            print *,"im_main or im_main_opp invalid"
            stop
           endif
           if ((iten_main.ge.1).and.(iten_main.le.num_interfaces)) then
            ! do nothing
           else
            print *,"iten_main invalid"
            stop
           endif

           if (LSIDE(1).ge.LSIDE(2)) then
            visc1=localvisc_minus(im_main)
            heat1=localheatvisc_minus(im_main)
            visc2=localvisc_plus(im_main_opp)
            heat2=localheatvisc_plus(im_main_opp)
           else if (LSIDE(1).le.LSIDE(2)) then
            visc1=localvisc_plus(im_main)
            heat1=localheatvisc_plus(im_main)
            visc2=localvisc_minus(im_main_opp)
            heat2=localheatvisc_minus(im_main_opp)
           else
            print *,"LSIDE bust"
            stop
           endif

           if ((heat1.ge.zero).and.(heat2.ge.zero)) then
            ! do nothing
           else
            print *,"heat1 or heat2 invalid"
            stop
           endif

           do imspec=1,num_species_var
            spec1(imspec)= &
             fort_speciesviscconst((imspec-1)*num_materials+im_main)
            spec2(imspec)= &
             fort_speciesviscconst((imspec-1)*num_materials+im_main_opp)
           enddo
  
             ! 1/s = (wL/sL + wR/sR)/(wL+wR)=(wL sR + wR sL)/(sL sR)*1/(wL+wR)
           if ((visc1.lt.zero).or.(visc2.lt.zero)) then
            print *,"visc1 or visc2 cannot be negative"
            stop
           else if ((visc1.eq.zero).or.(visc2.eq.zero)) then
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

           if (visc_interface(iten_main).eq.zero) then
            ! do nothing
           else if (visc_interface(iten_main).gt.zero) then

            if (test_for_quasi_static.eq.1) then
             ! do nothing
            else if (test_for_quasi_static.eq.0) then
             facevisc_local=visc_interface(iten_main)
            else
             print *,"test_for_quasi_static invalid"
             stop
            endif

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

           if (heatvisc_interface(iten_main).eq.zero) then
            ! do nothing
           else if (heatvisc_interface(iten_main).gt.zero) then

            if (get_user_latent_heat(iten_main,293.0d0,1).ne.zero) then
             if ((freezing_model(iten_main).eq.0).or. &
                 (freezing_model(iten_main).eq.5)) then
              print *,"heatvisc_interface invalid"
              stop
             endif 
            endif 
            if (get_user_latent_heat( &
                 iten_main+num_interfaces,293.0d0,1).ne.zero) then
             if ((freezing_model(iten_main+num_interfaces).eq.0).or. &
                 (freezing_model(iten_main+num_interfaces).eq.5)) then
              print *,"heatvisc_interface invalid"
              stop
             endif 
            endif 

            faceheat_local=heatvisc_interface(iten_main)
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

            spec_test= &
              speciesvisc_interface((imspec-1)*num_interfaces+iten_main)
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
           print *,"gradh bust; gradh=",gradh
           stop
          endif

           ! we always use the volume fractions for the viscosity
           ! and density coefficients.
           ! mu_face=sum F_i/(sum F_i/mu_i)
          voltotal=zero
          visc_total=zero
          do im=1,num_materials
           FFACE(im)=zero
          enddo
          is_zero_visc=0
          do iside=0,1
          do im=1,num_materials
           voldepart=local_face(FACECOMP_VOFFACE+2*(im-1)+iside+1)
           voltotal=voltotal+voldepart
           FFACE(im)=FFACE(im)+voldepart
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
             print *,"visc1 invalid; visc1=",visc1
             stop
            endif
           else if (voldepart.eq.zero) then
            ! do nothing
           else 
            print *,"voldepart invalid"
            stop
           endif   
          enddo ! im=1..num_materials
          enddo ! iside=0..1

          if (voltotal.gt.zero) then
           do im=1,num_materials
            FFACE(im)=FFACE(im)/voltotal
           enddo
          else
           print *,"voltotal invalid voltotal= ",voltotal
           stop
          endif

          if (is_zero_visc.eq.0) then
           facevisc_local=voltotal/visc_total
          else if (is_zero_visc.eq.1) then
           facevisc_local=zero
          else
           print *,"is_zero_visc invalid"
           stop
          endif

          do im=1,num_materials
           do im_opp=im+1,num_materials 
 
            if ((FFACE(im).gt.VOFTOL).and. &
                (FFACE(im_opp).gt.VOFTOL)) then
             call get_iten(im,im_opp,iten_FFACE)

             if (visc_interface(iten_FFACE).eq.zero) then
              ! do nothing
             else if (visc_interface(iten_FFACE).gt.zero) then
              if (test_for_quasi_static.eq.1) then
               ! do nothing
              else if (test_for_quasi_static.eq.0) then
               facevisc_local=visc_interface(iten_FFACE)
              else
               print *,"test_for_quasi_static invalid"
               stop
              endif
             else
              print *,"visc_interface invalid"
              stop
             endif

            else if ((FFACE(im).gt.-VOFTOL).and. &
                     (FFACE(im_opp).gt.-VOFTOL)) then
             ! do nothing
            else
             print *,"FFACE invalid"
             stop
            endif

           enddo ! im_opp=1..num_materials
          enddo ! im=1..num_materials

         else
          print *,"solid_present_flag bust"
          print *,"solid_present_flag=",solid_present_flag
          stop
         endif

          ! initialize wall_flag_face

         wall_flag_face=0
         wallVOF_face=zero

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

          call gridsten_level(xsten,icell,jcell,kcell,level,nhalf)

          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if (levelrz.eq.COORDSYS_RZ) then
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
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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

          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if ((levelrz.eq.COORDSYS_RZ).or. &
                   (levelrz.eq.COORDSYS_CYLINDRICAL)) then
           if (veldir.eq.0) then
            if (iside.eq.0) then
             if (abs(xstenMAC_center(1)).le.VOFTOL*dx(1)) then
              wall_flag_face=num_materials+1 !Neumann BC, RZ axis of symmetry
             else if (xstenMAC_center(1).gt.VOFTOL*dx(1)) then
              ! do nothing
             else
              print *,"xstenMAC invalid"
              stop
             endif
            else if (iside.eq.1) then
             if (xstenMAC_center(1).gt.-VOFTOL*dx(1)) then
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

           ! volume fraction for adjoining cell 
           ! iside==0 => left cell   iside==1 => right cell
           ! this "special" volume fraction will be =1 if the 
           ! total volume fraction of solid materials exceeds 1/2
          do im=1,num_materials
           volmat(im)=vofC(D_DECL(icell,jcell,kcell),im)
          enddo
           ! voltotal_prescribed=sum F_prescribed
           ! im_prescribed_primary=argmax_im F_prescribed
           ! combine_prescribed_VOF is declared in GLOBALUTIL.F90
          call combine_prescribed_VOF( &
            volmat, & !intent(in)
            voltotal_prescribed, & ! intent(out)
            im_prescribed_primary) ! intent(out)

          if (is_clamped_face.ge.1) then !interior "wall"

           wall_flag_face=num_materials+1

          else if (is_clamped_face.eq.0) then

           if ((voltotal_prescribed.ge.one-0.01D0).and. &
               (voltotal_prescribed.le.one+VOFTOL)) then
            if ((im_prescribed_primary.ge.1).and. &
                (im_prescribed_primary.le.num_materials)) then
             if (wall_flag_face.eq.0) then
              wall_flag_face=im_prescribed_primary 
               ! wallVOF_face should be 0 or 1 since the 
               ! reconstruction on input to fort_init_physics_vars
               ! "rasterizes" the "rigid" and
               ! "prescribed" materials.
              wallVOF_face=volmat(im_prescribed_primary)
             else if ((wall_flag_face.ge.1).and. &
                      (wall_flag_face.le.num_materials)) then
              if (volmat(im_prescribed_primary).gt.wallVOF_face) then
               wall_flag_face=im_prescribed_primary
               wallVOF_face=volmat(im_prescribed_primary)
              endif
             else if (wall_flag_face.eq.num_materials+1) then
              ! do nothing
             else
              print *,"wall_flag_face invalid"
              stop
             endif
            else
             print *,"im_prescribed_primary invalid"
             stop
            endif
           else if ((voltotal_prescribed.le.one-0.01D0).and. &
                    (voltotal_prescribed.ge.-VOFTOL)) then
            ! do nothing
           else
            print *,"voltotal_prescribed invalid"
            stop
           endif

          else
           print *,"is_clamped_face invalid"
           stop
          endif

          if ((wall_flag_face.ge.0).and. &
              (wall_flag_face.le.num_materials)) then

           prescribed_exists_flag=prescribed_exists()

           if (prescribed_exists_flag.eq.1) then 

             ! at least one of the face's adjoining cells is dominated
             ! by a prescribed material.
            if (is_prescribed_face.eq.1) then 
             if ((im_prescribed.ge.1).and. &
                 (im_prescribed.le.num_materials)) then
              if (wall_flag_face.eq.0) then
               wall_flag_face=im_prescribed
               wallVOF_face=volmat(im_prescribed)
              else if ((wall_flag_face.ge.1).and. &
                       (wall_flag_face.le.num_materials)) then
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

           else if (prescribed_exists_flag.eq.0) then
            ! do nothing
           else
            print *,"prescribed_exists_flag invalid 70"
            stop
           endif

          else if (wall_flag_face.eq.num_materials+1) then
           ! do nothing
          else
           print *,"wall_flag_face invalid"
           stop
          endif

         enddo ! iside=0..1   (initializing wall_flag_side)

         voltotal=zero
         mass_total=zero
         do im=1,num_materials
          FFACE(im)=zero
         enddo

         do iside=0,1
         do im=1,num_materials

          voldepart=local_face(FACECOMP_VOFFACE+2*(im-1)+iside+1)
          voltotal=voltotal+voldepart
          FFACE(im)=FFACE(im)+voldepart
          mass_total=mass_total+  &
           local_face(FACECOMP_MASSFACE+2*(im-1)+iside+1)
          
         enddo ! im=1..num_materials
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
         else if ((mass_total.gt.zero).and. &
                  (voltotal.gt.zero)) then
          ! do nothing
         else
          print *,"mass_total or voltotal is NaN"
          stop
         endif

         if (voltotal.gt.zero) then
          do im=1,num_materials
           FFACE(im)=FFACE(im)/voltotal
          enddo
         else 
          print *,"voltotal invalid (faceden) voltotal= ",voltotal
          stop
         endif

         density_for_mass_fraction_diffusion=mass_total/voltotal

         local_face(FACECOMP_FACEDEN+1)= &
            one/density_for_mass_fraction_diffusion

         local_face(FACECOMP_FACEDEN_BASE+1)=local_face(FACECOMP_FACEDEN+1)

         do im=1,num_materials
          do im_opp=im+1,num_materials

           call get_iten(im,im_opp,iten_FFACE)

           if (((FFACE(im).gt.VOFTOL).and. &
                (FFACE(im_opp).gt.VOFTOL)).or. &
               (iten_main.eq.iten_FFACE)) then
            if (denconst_interface_added(iten_FFACE).eq.zero) then
             ! do nothing
            else if (denconst_interface_added(iten_FFACE).gt.zero) then
             if (local_face(FACECOMP_FACEDEN+1).gt.zero) then !1/rho

              local_face(FACECOMP_FACEDEN+1)=one/ &
                 denconst_interface_added(iten_FFACE)

              if (local_face(FACECOMP_FACEDEN_BASE+1).gt. &
                  local_face(FACECOMP_FACEDEN+1)) then
               ! do nothing
              else
               print *,"local_face(FACECOMP_FACEDEN_BASE+1) invalid"
               stop
              endif

             else
              print *,"local_face(FACECOMP_FACEDEN+1) invalid"
              stop
             endif
            else
             print *,"denconst_interface_added invalid"
             stop
            endif
           else if (iten_main.ne.iten_FFACE) then
            ! do nothing
           else if ((FFACE(im).gt.-VOFTOL).and. &
                    (FFACE(im_opp).gt.-VOFTOL)) then
            ! do nothing
           else
            print *,"FFACE invalid"
            stop
           endif
          enddo ! im_opp=im+1..num_materials
         enddo ! im=1..num_materials

         local_face(FACECOMP_FACEVISC+1)=facevisc_local

         local_face(FACECOMP_FACEHEAT+1)=faceheat_local
         do imspec=1,num_species_var
          local_face(FACECOMP_FACESPEC+imspec)= &
            density_for_mass_fraction_diffusion*facespecies_local(imspec)
         enddo

          ! mask_boundary=1 at left neumann boundary
          ! mask_boundary=2 at right neumann boundary
          ! mask_boundary=0 otherwise
         if (mask_boundary.eq.1) then
          wall_flag_face=num_materials+1
         else if (mask_boundary.eq.2) then
          wall_flag_face=num_materials+1
         else if (mask_boundary.eq.0) then
          ! do nothing
         else 
          print *,"mask_boundary invalid"
          stop
         endif  

         local_face(FACECOMP_ICEFACECUT+1)=one

! local_face(FACECOMP_FACECUT+1)=0.0 if presbc=REFLECT_EVEN,LO_EXTRAP
! local_face(FACECOMP_FACECUT+1)=0.0 if face has adjoining
!  prescribed solid.
         if (wall_flag_face.eq.0) then
          local_face(FACECOMP_FACECUT+1)=one
         else if ((wall_flag_face.ge.1).and. &
                  (wall_flag_face.le.num_materials)) then
          local_face(FACECOMP_FACECUT+1)=zero
         else if (wall_flag_face.eq.num_materials+1) then
          local_face(FACECOMP_FACECUT+1)=zero
         else
          print *,"wall_flag_face invalid"
          stop
         endif

         ! compute curvature at the faces
         ! also compute model based pressure forcing at the faces

         zeroradius_flag=0

         if (levelrz.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if (levelrz.eq.COORDSYS_RZ) then
          if (SDIM.eq.2) then
           if (veldir.eq.0) then
            if (abs(xstenMAC_center(1)).le.VOFTOL*dx(1)) then ! at r=0 face
             zeroradius_flag=1
            else if (xstenMAC_center(1).ge.-VOFTOL*dx(1)) then
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
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
          ! do nothing
         else
          print *,"levelrz invalid init physics vars 4"
          stop
         endif

         local_face(FACECOMP_CURV+1)=zero

         if (local_face(FACECOMP_FACECUT+1).ge.zero) then
          ! do nothing
         else
          print *,"local_face(FACECOMP_FACECUT+1) invalid"
          print *,"local_face(FACECOMP_FACECUT+1)=", &
            local_face(FACECOMP_FACECUT+1)
          stop
         endif

         if (local_face(FACECOMP_FACECUT+1).lt.zero) then
          print *,"local_face(FACECOMP_FACECUT+1).lt.zero"
          stop
         else if ((local_face(FACECOMP_FACECUT+1).ge.zero).and. &
                  (local_face(FACECOMP_FACECUT+1).le.half)) then
          ! do nothing
         else if (zeroradius_flag.eq.1) then
          ! do nothing
         else if ((zeroradius_flag.eq.0).and. &
                  (local_face(FACECOMP_FACECUT+1).ge.half).and. &
                  (local_face(FACECOMP_FACECUT+1).le.one)) then

          if (gradh_tension.ne.zero) then

           if ((im_tension.gt.num_materials).or. &
               (im_opp_tension.gt.num_materials).or. &
               (im_tension.lt.1).or. &
               (im_opp_tension.lt.1)) then
            print *,"im_tension or im_opp_tension bust 4"
            stop
           endif
           call get_iten(im_tension,im_opp_tension,iten_tension)
  
! vof,ref centroid,order,slope,intercept  x num_materials

           if (LS_consistent_tension.eq.1) then

            if (is_solid_face.eq.1) then
             local_face(FACECOMP_CURV+1)=zero
            else if (is_solid_face.eq.0) then
             icurv=(iten_tension-1)*CURVCOMP_NCOMP

             if ((level.ge.0).and.(level.le.finest_level)) then

              if ((project_option.eq.SOLVETYPE_PRES).or. &
                  (project_option.eq.SOLVETYPE_INITPROJ)) then

               do icurv_ofs=1,CURVCOMP_NCOMP
                curvL(icurv_ofs)=curv(D_DECL(im1,jm1,km1),icurv+icurv_ofs)
                curvR(icurv_ofs)=curv(D_DECL(i,j,k),icurv+icurv_ofs)
               enddo
               dirL=NINT(curvL(CURVCOMP_DIRSIDE_FLAG+1))
               sideL=1
               if (dirL.lt.0) then
                dirL=-dirL
                sideL=-sideL
               endif
               dirR=NINT(curvR(CURVCOMP_DIRSIDE_FLAG+1))
               sideR=1
               if (dirR.lt.0) then
                dirR=-dirR
                sideR=-sideR
               endif
               im3L=NINT(curvL(CURVCOMP_MATERIAL3_ID+1))
               im3R=NINT(curvR(CURVCOMP_MATERIAL3_ID+1))

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
                  do icurv_ofs=1,CURVCOMP_NCOMP
                   curvL(icurv_ofs)=zero
                   curvR(icurv_ofs)=zero
                  enddo
                 endif
                else
                 print *,"dirL or dirR invalid"
                 stop
                endif
               endif

               if ((im3L.lt.0).or.(im3L.gt.num_materials).or. &
                   (im3R.lt.0).or.(im3R.gt.num_materials).or. &
                   (dirL.lt.0).or.(dirL.gt.SDIM+1).or. &
                   (dirR.lt.0).or.(dirR.gt.SDIM+1)) then
                print *,"curvature parameters invalid"
                print *,"im3L,im3R ",im3L,im3R
                print *,"dirL,dirR ",dirL,dirR
                print *,"im_tension ",im_tension
                print *,"im_opp_tension ",im_opp_tension
                print *,"iten_tension ",iten_tension
                print *,"i,j,k,veldir ",i,j,k,veldir
                do im=1,num_materials
                 print *,"im,LS_LEFT,LS_RIGHT ",im, &
                  levelPC(D_DECL(im1,jm1,km1),im),levelPC(D_DECL(i,j,k),im)
                enddo
                stop
               endif

                ! inputs.curvature_converge, continuous_mof=1
                ! March 10, 2018: 1.99, 2.03 RZ 24x48 HT
                ! March 10, 2018: 1.00, 1.01 XY 24x48 HT
                ! March 10, 2018: 1.93, 2.07 XYZ 32x32x32 HT
               if (((im3L.ge.1).and.(im3L.le.num_materials)).or. &
                   ((im3R.ge.1).and.(im3R.le.num_materials))) then
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
                local_face(FACECOMP_CURV+1)= &
                 wtL*curvL(CURVCOMP_FD_CURV+1)+ &
                 wtR*curvR(CURVCOMP_FD_CURV+1)
               else if (curv_interp_flag.eq.1) then   
                if (wtL.gt.wtR) then
                 local_face(FACECOMP_CURV+1)=curvL(CURVCOMP_HTFUNC_CURV+1)
                else if (wtR.gt.wtL) then
                 local_face(FACECOMP_CURV+1)=curvR(CURVCOMP_HTFUNC_CURV+1)
                else if (wtR.eq.wtL) then
                 local_face(FACECOMP_CURV+1)= &
                    wtL*curvL(CURVCOMP_HTFUNC_CURV+1)+ &
                    wtR*curvR(CURVCOMP_HTFUNC_CURV+1)
                else
                 print *,"wtR or wtL is NaN"
                 stop
                endif
               else if (curv_interp_flag.eq.2) then
                local_face(FACECOMP_CURV+1)=curvL(CURVCOMP_HTFUNC_CURV+1)
               else if (curv_interp_flag.eq.3) then
                local_face(FACECOMP_CURV+1)=curvR(CURVCOMP_HTFUNC_CURV+1)
               else if (curv_interp_flag.eq.4) then
                local_face(FACECOMP_CURV+1)= &
                   wtL*curvL(CURVCOMP_HTFUNC_CURV+1)+ &
                   wtR*curvR(CURVCOMP_HTFUNC_CURV+1)
               else if (curv_interp_flag.eq.5) then
                if (wtL.gt.wtR) then
                 local_face(FACECOMP_CURV+1)=curvL(CURVCOMP_FD_CURV+1)
                else if (wtR.gt.wtL) then
                 local_face(FACECOMP_CURV+1)=curvR(CURVCOMP_FD_CURV+1)
                else if (wtR.eq.wtL) then
                 local_face(FACECOMP_CURV+1)= &
                  wtL*curvL(CURVCOMP_FD_CURV+1)+ &
                  wtR*curvR(CURVCOMP_FD_CURV+1)
                else
                 print *,"wtR or wtL is NaN"
                 stop
                endif
               else
                print *,"curv_interp_flag invalid"
                stop
               endif

               if (curv_min.gt.local_face(FACECOMP_CURV+1)) then
                curv_min=local_face(FACECOMP_CURV+1)
               endif
               if (curv_max.lt.local_face(FACECOMP_CURV+1)) then
                curv_max=local_face(FACECOMP_CURV+1)
               endif

              else
               print *,"project_option invalid fort_init_physics_vars"
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
            do im=1,num_materials
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
          print *,"zeroradius_flag or local_face(FACECOMP_FACECUT+1) invalid"
          stop
         endif 

         if (test_for_quasi_static.eq.1) then
          local_face(FACECOMP_FACEVEL+1)=zero
         else if (test_for_quasi_static.eq.0) then
          ! do nothing
         else
          print *,"test_for_quasi_static invalid"
          stop
         endif

         do im=1,FACECOMP_NCOMP
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
         enddo  ! im=1,FACECOMP_NCOMP

        enddo
        enddo
        enddo  ! i,j,k
       enddo ! veldir

      else if (isweep.eq.1) then

        ! cenden, cenvof, cenDeDT, 
        ! cenvisc, initialized in this loop.

       call growntilebox(tilelo,tilehi,fablo,fabhi,igridlo,igridhi,1) 

       do i=igridlo(1),igridhi(1)
       do j=igridlo(2),igridhi(2)
       do k=igridlo(3),igridhi(3)

        call gridsten_level(xsten,i,j,k,level,nhalf)
        do dir2=1,SDIM
         xclamped_minus(dir2)=xsten(0,dir2)
         xclamped_plus(dir2)=xsten(0,dir2)
        enddo

        do im=1,num_materials*ngeom_recon
         mofdata(im)=slope(D_DECL(i,j,k),im)
        enddo
         ! before (mofdata): fluids tessellate
         ! after  (mofdata): fluids and solids tessellate
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
         SDIM)

        voltotal=zero
        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         volmat(im)=mofdata(vofcomp)
         voltotal=voltotal+volmat(im)
        enddo ! im=1..num_materials

        if (abs(voltotal-one).gt.LSTOL) then
         print *,"voltotal invalid"
         stop
        else if (abs(voltotal-one).le.LSTOL) then
         ! do nothing
        else
         print *,"voltotal invalid"
         stop
        endif

        DeDT_total=zero

        visc_total=zero

        mass_total=zero

         ! Cell grid 1/rho and mu
        null_viscosity=0

        do im=1,num_materials
         LSIDE_MAT(im)=levelPC(D_DECL(i,j,k),im)
        enddo

         ! checks rigid and non-rigid materials.
         ! get_primary_material is declared in: GLOBALUTIL.F90
         ! get_secondary_material is declared in: MOF.F90
        call get_primary_material(LSIDE_MAT,implus_majority)
        call get_secondary_material(LSIDE_MAT,implus_majority,im_secondary)

         ! LS>0 if clamped
        call SUB_clamped_LS(xclamped_minus,time,LS_clamped_minus, &
          vel_clamped_minus,temperature_clamped_minus,prescribed_flag,dx)

        if (LS_clamped_minus.ge.zero) then
         null_viscosity=1
        else if (LS_clamped_minus.lt.zero) then

         if (is_rigid(implus_majority).eq.1) then
          null_viscosity=1
         else if (is_rigid(implus_majority).eq.0) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif

         if (is_prescribed(implus_majority).eq.1) then
          null_viscosity=1
         else if (is_prescribed(implus_majority).eq.0) then
          ! do nothing
         else
          print *,"is_prescibed invalid"
          stop
         endif

         if (is_ice_or_FSI_rigid_material(implus_majority).eq.1) then
          null_viscosity=1
         else if (is_ice_or_FSI_rigid_material(implus_majority).eq.0) then
          ! do nothing
         else
          print *,"is_ice_or_FSI_rigid_material invalid:fort_init_physics_vars"
          stop
         endif

        else
         print *,"LS_clamped_minus invalid: fort_init_physics_vars"
         stop
        endif

        denconst_interface_added_max=zero

        if ((abs(LSIDE_MAT(implus_majority)).le.rad_added_mass).and. &
            (abs(LSIDE_MAT(im_secondary)).le.rad_added_mass)) then
         call get_iten(implus_majority,im_secondary,iten_main)
         if (denconst_interface_added(iten_main).eq.zero) then
          ! do nothing
         else if (denconst_interface_added(iten_main).gt.zero) then
          if (max(fort_denconst(implus_majority),fort_denconst(im_secondary)) &
              .le.denconst_interface_added(iten_main)) then
           denconst_interface_added_max= &
             denconst_interface_added(iten_main)
          else
           print *,"denconst_interface_added(iten_main) invalid"
           stop
          endif
         else 
          print *,"denconst_interface_added(iten_main) invalid"
          stop
         endif 
        else if ((abs(LSIDE_MAT(implus_majority)).ge.rad_added_mass).or.  &
                 (abs(LSIDE_MAT(im_secondary)).ge.rad_added_mass)) then
         ! do nothing
        else
         print *,"LSIDE_MAT invalid: fort_init_physics_vars"
         stop
        endif  

        do im=1,num_materials

         dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
         tempcomp=dencomp+1

!        delta_mass=denstate(D_DECL(i,j,k),dencomp)*volmat(im)
         delta_mass=mom_den(D_DECL(i,j,k),im)*volmat(im)
         voldepart=volmat(im)
          ! if is_rigid(im), density=fort_denconst(im)
          ! if incompressible,
          !   if constant_density_all_time==1 then density=fort_denconst(im)
          !   if constant_density_all_time==0 
          !    then density=mass_depart/vol_depart
          ! if compressible,
          !   density=massdepart/voltarget
         call derive_density( &
          voldepart,voldepart,voltotal, &
          constant_density_all_time, &
          delta_mass, &
          im, &
          den) 

         if (den.gt.zero) then
          ! do nothing
         else
          print *,"density must be positive init_phyiscs_vars2"
          print *,"i,j,k,im,den ",i,j,k,im,den
          stop
         endif

         if (test_for_quasi_static.eq.1) then
          localvisc(im)=zero
         else if (test_for_quasi_static.eq.0) then
          localvisc(im)=viscstate(D_DECL(i,j,k),im)
         else
          print *,"test_for_quasi_static invalid"
          stop
         endif

         if (localvisc(im).lt.zero) then
          print *,"viscstate gone negative"
          stop
         else if (localvisc(im).eq.zero) then
          if (volmat(im).gt.zero) then
           null_viscosity=1
          else if (volmat(im).eq.zero) then
           ! do nothing
          else
           print *,"volmat(im) went negative"
           print *,"im,volmat ",im,volmat(im)
           stop
          endif
         else if (localvisc(im).gt.zero) then
          ! do nothing
         else
          print *,"localvisc NaN"
          stop
         endif

         one_over_mu=one/(localvisc(im)+VISCINVTOL)
         visc_total=visc_total+one_over_mu*volmat(im)

         imattype=fort_material_type(im)
         TEMPERATURE=denstate(D_DECL(i,j,k),tempcomp)
         if (TEMPERATURE.gt.zero) then
          ! do nothing
         else
          print *,"PHYSICS_VAR: temperature must be positive"
          print *,"num_materials ",num_materials
          print *,"num_state_material ",num_state_material
          print *,"fort_initial_temperature(im) ",fort_initial_temperature(im)
          print *,"fort_tempconst(im) ",fort_tempconst(im)
          print *,"i,j,k,im ",i,j,k,im
          print *,"den,TEMPERATURE ",den,TEMPERATURE
          stop
         endif
         call init_massfrac_parm(den,massfrac_parm,im)
         do ispec=1,num_species_var
          massfrac_parm(ispec)=denstate(D_DECL(i,j,k),tempcomp+ispec)
          if (massfrac_parm(ispec).ge.zero) then
           ! do nothing
          else
           print *,"massfrac_parm(ispec) invalid"
           stop
          endif
         enddo ! ispec=1,num_species_var

          ! DeDT = cv
          ! DeDT_material is declared in: GLOBALUTIL.F90
          ! if imattype!=18,15 then
          !  DeDT_material = (e(T+DT)-e(T))/DT
          ! e(T) is subroutine INTERNAL_material and is declared in:
          ! GLOBALUTIL.F90.
          ! if (is_in_probtype_list().eq.1) then a call to, SUB_INTERNAL,
          ! is made.
         call DeDT_material(den,massfrac_parm, &
           TEMPERATURE,DeDT,imattype,im)
         if (DeDT.gt.zero) then
          ! do nothing
         else
          print *,"DeDT must be positive"
          stop
         endif

         delta_mass=den*volmat(im)
         mass_total=mass_total+delta_mass
         DeDT_total=DeDT_total+DeDT*delta_mass
        enddo ! im=1..num_materials

        if (mass_total.gt.zero) then
         ! do nothing
        else
         print *,"mass_total invalid"
         stop
        endif

        if (DeDT_total.gt.zero) then
         ! do nothing
        else
         print *,"DeDT_total must be positive"
         stop
        endif

        do im=1,num_materials
         cenvof(D_DECL(i,j,k),im)=volmat(im)/voltotal
        enddo

        local_cenden=voltotal/mass_total
        local_cenvisc=one/((visc_total/voltotal)+VISCINVTOL)  

        if (local_cenden.gt.zero) then
         !do nothing
        else
         print *,"local_cenden invalid fort_init_physics_vars"
         stop
        endif

        if (denconst_interface_added_max.eq.zero) then
         ! do nothing
        else if (denconst_interface_added_max.gt.zero) then
         if ((local_cenden.gt.zero).and. &
             (denconst_interface_added_max.gt.one/local_cenden)) then
          !local_cenden=voltotal/mass_total
          local_cenden=one/denconst_interface_added_max
         else
          print *,"local_cenden or denconst_interface_added_max invalid"
          print *,"local_cenden=",local_cenden
          print *,"denconst_interface_added_max=",denconst_interface_added_max
          stop
         endif
        else
         print *,"denconst_interface_added_max invalid"
         stop
        endif

        cenden(D_DECL(i,j,k))=local_cenden
        cenDeDT(D_DECL(i,j,k))=voltotal/DeDT_total

        if (null_viscosity.eq.1) then
         cenvisc(D_DECL(i,j,k))=zero
        else if (null_viscosity.eq.0) then
         cenvisc(D_DECL(i,j,k))=local_cenvisc
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

      deallocate(fort_caller_string)

      return
      end subroutine fort_init_physics_vars


      subroutine fort_build_semirefinevof( &
       tid, &
       tessellate, &
       ngrow_refine, &
       nrefine_vof, &
       spec_material_id_AMBIENT, &
       mass_fraction_id, &
       cavitation_vapor_density, &
       xlo,dx, &
       slope,DIMS(slope), &
       denstate, &
       DIMS(denstate), &
       mom_den, &
       DIMS(mom_den), &
       vofF,DIMS(vofF), &
       massF,DIMS(massF), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level,finest_level) &
      bind(c,name='fort_build_semirefinevof')

      use global_utility_module
      use probf90_module
      use geometry_intersect_module
      use MOF_routines_module
      use godunov_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: tessellate
      INTEGER_T, INTENT(in) :: ngrow_refine
      INTEGER_T, INTENT(in) :: tid
      INTEGER_T, INTENT(in) :: nrefine_vof
      INTEGER_T, INTENT(in) :: spec_material_id_AMBIENT(num_species_var+1)
      INTEGER_T, INTENT(in) :: mass_fraction_id(2*num_interfaces)
      REAL_T, INTENT(in) :: cavitation_vapor_density(num_materials)
      INTEGER_T, INTENT(in) :: level,finest_level
      INTEGER_T :: veldir
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: DIMDEC(slope)
      INTEGER_T, INTENT(in) :: DIMDEC(denstate)
      INTEGER_T, INTENT(in) :: DIMDEC(mom_den)
      INTEGER_T, INTENT(in) :: DIMDEC(vofF)
      INTEGER_T, INTENT(in) :: DIMDEC(massF)
      REAL_T, INTENT(in), target :: slope(DIMV(slope),num_materials*ngeom_recon) 
      REAL_T, pointer :: slope_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: denstate(DIMV(denstate), &
              num_materials*num_state_material) 
      REAL_T, pointer :: denstate_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: mom_den(DIMV(mom_den),num_materials) 
      REAL_T, pointer :: mom_den_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out), target :: vofF(DIMV(vofF),nrefine_vof)
      REAL_T, pointer :: vofF_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out), target :: massF(DIMV(massF),nrefine_vof)
      REAL_T, pointer :: massF_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k
      INTEGER_T dir2
      INTEGER_T iside

      INTEGER_T im,nmax
      REAL_T mofdata(num_materials*ngeom_recon)

      REAL_T multi_volume(num_materials)
      REAL_T multi_cen(SDIM,num_materials)
      REAL_T voltotal_fluid
      REAL_T mass_total_fluid
      REAL_T voltotal_solid
      REAL_T mass_total_solid
      REAL_T voltotal
      REAL_T mass_total
      REAL_T den
      REAL_T mom_den_local
      REAL_T den_value

      INTEGER_T dencomp

      INTEGER_T igridlo(3),igridhi(3)

      REAL_T volrecon
      REAL_T cenrecon(SDIM)

      REAL_T voldonate
      REAL_T cendonate(SDIM)

      INTEGER_T check_donate
      INTEGER_T irefine
      REAL_T xsten_recon(-1:1,SDIM)
      REAL_T xsten_donate(-1:1,SDIM)
      REAL_T mu

      slope_ptr=>slope
      denstate_ptr=>denstate
      mom_den_ptr=>mom_den

      vofF_ptr=>vofF
      massF_ptr=>massF

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif
      if ((tessellate.ne.0).and. &
          (tessellate.ne.1).and. &
          (tessellate.ne.3)) then
       print *,"tessellate invalid2"
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
      if (nrefine_vof.ne.2*num_materials*SDIM) then
       print *,"nrefine_vof invalid"
       stop
      endif

      nmax=POLYGON_LIST_MAX ! in: BUILD_SEMIREFINEVOF

      call checkbound_array(fablo,fabhi,slope_ptr,ngrow_refine,-1)
      call checkbound_array(fablo,fabhi,denstate_ptr,ngrow_refine,-1)
      call checkbound_array(fablo,fabhi,mom_den_ptr,ngrow_refine,-1)
      call checkbound_array(fablo,fabhi,vofF_ptr,ngrow_refine,-1)
      call checkbound_array(fablo,fabhi,massF_ptr,ngrow_refine,-1)

      do im=1,num_materials

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        if (is_rigid(im).ne.1) then
         print *,"is_rigid(im).ne.1"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
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
        print *,"energy must be positive in fort_build_semirefinevof"
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

      enddo ! im=1..num_materials  (checking parameters)

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
         
          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if ((levelrz.eq.COORDSYS_RZ).or. &
                   (levelrz.eq.COORDSYS_CYLINDRICAL)) then
           if (xsten_recon(0,1).lt.VOFTOL*dx(1)) then
            check_donate=0
           endif
          else
           print *,"levelrz invalid build semi refine vof"
           stop
          endif
 
          if (check_donate.eq.0) then !RZ R=0 or RTZ R=0
           do im=1,num_materials

            if (iside.eq.-1) then
             irefine=veldir*2*num_materials+im
            else if (iside.eq.1) then
             irefine=veldir*2*num_materials+num_materials+im
            else
             print *,"iside invalid"
             stop
            endif

            vofF(D_DECL(i,j,k),irefine)=zero
            massF(D_DECL(i,j,k),irefine)=zero
           enddo ! im=1..num_materials
          else if (check_donate.eq.1) then

           call CISBOXHALF(xsten_donate,1, &
            xlo,dx,i,j,k,iside,veldir+1, &
            bfact,level, & 
            voldonate,cendonate,SDIM)

           do dir2=1,num_materials*ngeom_recon
            mofdata(dir2)=slope(D_DECL(i,j,k),dir2)
           enddo 
           ! multi_cen is "absolute" (not relative to cell centroid)
           call multi_get_volume_grid_simple( &
             tessellate, &  !=0,1, or 3
             bfact,dx,xsten_recon,1, &
             mofdata, &
             xsten_donate,1, &
             multi_volume,multi_cen, &
             geom_xtetlist(1,1,1,tid+1), &
             nmax, &
             nmax, &
             SDIM)
       
           mass_total_fluid=zero
           voltotal_fluid=zero 
           mass_total_solid=zero
           voltotal_solid=zero 
           do im=1,num_materials
            dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
            den=denstate(D_DECL(i,j,k),dencomp)
            mom_den_local=mom_den(D_DECL(i,j,k),im)

            den_value=mom_den_local

            if (den.gt.zero) then
             ! do nothing
            else
             print *,"den must be positive build_semi_refine_vof"
             print *,"im,den ",im,den
             print *,"im,fort_denconst(im) ",im,fort_denconst(im)
             print *,"level,finest_level ",level,finest_level
             stop
            endif  

            if (mom_den_local.gt.zero) then
             ! do nothing
            else
             print *,"mom_den_local must be pos build_semi_refine_vof"
             print *,"im,mom_den_local ",im,mom_den_local
             print *,"im,fort_denconst(im) ",im,fort_denconst(im)
             print *,"level,finest_level ",level,finest_level
             stop
            endif  
            if (den_value.gt.zero) then
             ! do nothing
            else
             print *,"den_value must be pos build_semi_refine_vof"
             print *,"im,den_value ",im,den_value
             print *,"im,fort_denconst(im) ",im,fort_denconst(im)
             print *,"level,finest_level ",level,finest_level
             stop
            endif  
 
            if (is_rigid(im).eq.0) then
             voltotal_fluid=voltotal_fluid+multi_volume(im)
             mass_total_fluid=mass_total_fluid+den_value*multi_volume(im)
            else if (is_rigid(im).eq.1) then
             voltotal_solid=voltotal_solid+multi_volume(im)
             mass_total_solid=mass_total_solid+den_value*multi_volume(im)
            else
             print *,"is_rigid invalid LEVELSET_3D.F90"
             stop
            endif
  
            if (iside.eq.-1) then 
             irefine=veldir*2*num_materials+im
            else if (iside.eq.1) then
             irefine=veldir*2*num_materials+num_materials+im
            else
             print *,"iside invalid"
             stop
            endif

            vofF(D_DECL(i,j,k),irefine)=multi_volume(im)
            massF(D_DECL(i,j,k),irefine)=den_value*multi_volume(im)
           enddo ! im=1,num_materials

           if (tessellate.eq.0) then
            voltotal=voltotal_fluid
            mass_total=mass_total_fluid
           else if ((tessellate.eq.1).or. &
                    (tessellate.eq.3)) then
            voltotal=voltotal_fluid+voltotal_solid
            mass_total=mass_total_fluid+mass_total_solid
           else
            print *,"tessellate invalid3"
            stop
           endif
 
           if ((voltotal.gt.zero).and.(mass_total.gt.zero)) then
            ! do nothing
           else
            print *,"voltotal or mass_total invalid"
            print *,"voltotal, mass_total, num_materials ", &
                    voltotal,mass_total,num_materials
            print *,"voldonate,volrecon ",voldonate,volrecon
            print *,"veldir ",veldir
            print *,"fablo ",fablo(1),fablo(2),fablo(SDIM)
            print *,"fabhi ",fabhi(1),fabhi(2),fabhi(SDIM)
            print *,"i,j,k ",i,j,k
            print *,"level,finest_level ",level,finest_level
            stop
           endif
           if ((voltotal_solid.ge.zero).and.(mass_total_solid.ge.zero)) then
            ! do nothing
           else
            print *,"voltotal_solid or mass_total_solid invalid"
            print *,"level,finest_level ",level,finest_level
            stop
           endif
           if ((voltotal_fluid.ge.zero).and.(mass_total_fluid.ge.zero)) then
            ! do nothing
           else
            print *,"voltotal_fluid or mass_total_fluid invalid"
            print *,"level,finest_level ",level,finest_level
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
      end subroutine fort_build_semirefinevof


!if compressible material,
! add beta * (1/cv) * (u dot u/2) to temp
      subroutine fort_inc_temp( &
       beta, &
       level, &
       finest_level, &
       ncomp_state, &
       tilelo,tilehi, &
       fablo,fabhi, &
       state,DIMS(state), &
       maskcoef,DIMS(maskcoef)) & ! 1=not covered  0=covered
       bind(c,name='fort_inc_temp')

       use probf90_module
       use global_utility_module
       IMPLICIT NONE

      REAL_T, INTENT(in) :: beta
      INTEGER_T, INTENT(in) :: ncomp_state
      INTEGER_T, INTENT(in) :: level,finest_level
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: DIMDEC(state)
      INTEGER_T, INTENT(in) :: DIMDEC(maskcoef)
      REAL_T, INTENT(inout), target :: state(DIMV(state),ncomp_state)
      REAL_T, pointer :: state_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: maskcoef(DIMV(maskcoef))
      REAL_T, pointer :: maskcoef_ptr(D_DECL(:,:,:))

      INTEGER_T i,j,k
      INTEGER_T dir
      INTEGER_T im
      INTEGER_T imattype
      INTEGER_T local_mask
      INTEGER_T vofcomp,dencomp
      REAL_T vof,KE,rho,TEMPERATURE,internal_e
      REAL_T massfrac_parm(num_species_var+1)
      INTEGER_T ispec

      if ((level.le.finest_level).and.(level.ge.0)) then
       ! do nothing
      else
       print *,"level invalid fort_inc_temp"
       stop
      endif

      do im=1,num_materials
       imattype=fort_material_type(im)
       if (imattype.eq.999) then
        ! do nothing
       else if ((imattype.ge.0).and. &
                (imattype.le.MAX_NUM_EOS)) then
        ! do nothing
       else
        print *,"imattype invalid fort_inc_temp"
        stop
       endif

      enddo ! im=1..num_materials

      if ((beta.eq.one).or.(beta.eq.-one)) then
       ! do nothing
      else
       print *,"beta invalid"
       stop
      endif 
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      do im=1,num_materials
       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"denconst invalid"
        stop
       endif
      enddo

      if (ncomp_state.ne.STATE_NCOMP) then
       print *,"ncomp_state invalid"
       stop
      endif

      state_ptr=>state
      call checkbound_array(fablo,fabhi,state_ptr,0,-1)
      maskcoef_ptr=>maskcoef
      call checkbound_array1(fablo,fabhi,maskcoef_ptr,0,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_mask=NINT(maskcoef(D_DECL(i,j,k)))
       if (local_mask.eq.1) then ! not covered

        do im=1,num_materials
         if (is_compressible_mat(im).eq.0) then
          ! do nothing
         else if (is_compressible_mat(im).eq.1) then
          vofcomp=STATECOMP_MOF+(im-1)*ngeom_raw+1
          vof=state(D_DECL(i,j,k),vofcomp)
          if ((vof.ge.-VOFTOL).and.(vof.le.one+VOFTOL)) then

           if (vof.lt.half) then
            ! do nothing
           else if (vof.ge.half) then
            KE=zero
            do dir=1,SDIM
             KE=KE+half*(state(D_DECL(i,j,k),dir)**2)
            enddo
            dencomp=STATECOMP_STATES+(im-1)*num_state_material+1+ENUM_DENVAR
            rho=state(D_DECL(i,j,k),dencomp)
            TEMPERATURE=state(D_DECL(i,j,k),dencomp+1)
            if (rho.gt.zero) then
             ! do nothing
            else
             print *,"density underflow"
             print *,"i,j,k,im,rho ",i,j,k,im,rho
             stop
            endif
             ! the first part is always to add the old kinetic energy.
            if (TEMPERATURE.gt.zero) then
             ! do nothing
            else
             print *,"temperature underflow"
             stop
            endif
            imattype=fort_material_type(im)
            if (is_compressible_mat(im).eq.1) then

             call init_massfrac_parm(rho,massfrac_parm,im)
             do ispec=1,num_species_var
              massfrac_parm(ispec)= &
                state(D_DECL(i,j,k),dencomp+1+ispec)
              if (massfrac_parm(ispec).ge.zero) then
               ! do nothing
              else
               print *,"massfrac_parm invalid"
               stop
              endif
             enddo ! ispec=1..num_species_var

             call INTERNAL_material(rho,massfrac_parm, &
               TEMPERATURE,internal_e, &
               imattype,im)
             internal_e=internal_e+beta*KE
             if (internal_e.gt.zero) then
              ! do nothing
             else
              print *,"internal_e.le.zero in fort_inc_temp"
              stop
             endif
             call TEMPERATURE_material(rho,massfrac_parm, &
               TEMPERATURE, &
               internal_e,imattype,im)
             if (TEMPERATURE.gt.zero) then
              state(D_DECL(i,j,k),dencomp+1)=TEMPERATURE
             else
              print *,"TEMPERATURE underflow in fort_inc_temp"
              stop
             endif
            else
             print *,"is_compressible_mat invalid"
             print *,"imattype invalid fort_inc_temp"
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
          print *,"is_compressible_mat invalid"
          stop
         endif
        enddo ! im=1..num_materials

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
      end subroutine fort_inc_temp

       ! OP_RHS_CELL
       ! operation_flag=100 (right hand side for solver)
       ! OP_DIV_CELL
       ! operation_flag=110 (divergence)
       ! OP_VEL_MAC_TO_CELL
       ! operation_flag=103 (mac -> cell velocity in solver or MAC_TO_CELL)
       ! OP_VEL_DIVUP_TO_CELL
       ! operation_flag=101 (copy u^{mac->cell} to u^{cell}, div(up) )
       ! OP_GRADU_MAC_TO_CELL
       ! operation_flag=105 (from face_gradients, interp grad U^T)
       ! OP_ISCHEME_CELL 
       ! operation_flag=107 (advection)
      subroutine fort_mac_to_cell( &
       ns_time_order, &
       divu_outer_sweeps, &
       num_divu_outer_sweeps, &
       operation_flag, &
       energyflag, &
       constant_density_all_time, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       level, &
       finest_level, &
       project_option, &
       enable_spectral, &
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
       xface,DIMS(xface), & !local_face_var_mf
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
       solxfab,DIMS(solxfab), & !local_fsi_ghost_mac_mf
       solyfab,DIMS(solyfab), &
       solzfab,DIMS(solzfab), &
       cterm,DIMS(cterm), &
       pold,DIMS(pold), &
       denold,DIMS(denold), &
       ustar,DIMS(ustar), &
       mdotcell,DIMS(mdotcell), & !VELADVECT_MF if OP_ISCHEME_CELL
       maskdivres,DIMS(maskdivres), & !DEN_RECON_MF if OP_ISCHEME_CELL
       maskres,DIMS(maskres), &
       SDC_outer_sweeps, &
       homflag, &
       nsolve, &
       ncomp_denold, &
       ncomp_veldest, &
       ncomp_dendest) &
       bind(c,name='fort_mac_to_cell')

       use probf90_module
       use global_utility_module
       use MOF_routines_module
       IMPLICIT NONE

      INTEGER_T, INTENT(in) :: ncomp_denold
      INTEGER_T, INTENT(in) :: ncomp_veldest
      INTEGER_T, INTENT(in) :: ncomp_dendest
      INTEGER_T, INTENT(in) :: ns_time_order
      INTEGER_T, INTENT(in) :: divu_outer_sweeps
      INTEGER_T, INTENT(in) :: num_divu_outer_sweeps
      INTEGER_T, INTENT(in) :: operation_flag
      INTEGER_T, INTENT(in) :: slab_step
      INTEGER_T, INTENT(in) :: enable_spectral
      INTEGER_T, INTENT(in) :: SDC_outer_sweeps 
      INTEGER_T, INTENT(in) :: energyflag 
      INTEGER_T, INTENT(in) :: constant_density_all_time(num_materials)
      INTEGER_T, INTENT(in) :: nparts
      INTEGER_T, INTENT(in) :: nparts_def
      INTEGER_T, INTENT(in) :: im_solid_map(nparts_def)
      INTEGER_T, INTENT(in) :: nsolve
      INTEGER_T, INTENT(in) :: homflag
      INTEGER_T, INTENT(in) :: level,finest_level
      INTEGER_T, INTENT(in) :: project_option
      INTEGER_T, INTENT(in) :: ncphys
      INTEGER_T, INTENT(in) :: velbc_in(SDIM,2,SDIM)
      INTEGER_T, INTENT(in) :: presbc_in(SDIM,2)
      REAL_T, INTENT(in) :: cur_time,dt
      REAL_T, INTENT(in) :: xlo(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: DIMDEC(xp)
      INTEGER_T, INTENT(in) :: DIMDEC(yp)
      INTEGER_T, INTENT(in) :: DIMDEC(zp)
      INTEGER_T, INTENT(in) :: DIMDEC(xvel)
      INTEGER_T, INTENT(in) :: DIMDEC(yvel)
      INTEGER_T, INTENT(in) :: DIMDEC(zvel)
      INTEGER_T, INTENT(in) :: DIMDEC(xface)
      INTEGER_T, INTENT(in) :: DIMDEC(yface)
      INTEGER_T, INTENT(in) :: DIMDEC(zface)
      INTEGER_T, INTENT(in) :: DIMDEC(ax)
      INTEGER_T, INTENT(in) :: DIMDEC(ay)
      INTEGER_T, INTENT(in) :: DIMDEC(az)
      INTEGER_T, INTENT(in) :: DIMDEC(vol)
      INTEGER_T, INTENT(in) :: DIMDEC(rhs)
      INTEGER_T, INTENT(in) :: DIMDEC(veldest)
      INTEGER_T, INTENT(in) :: DIMDEC(dendest)
      INTEGER_T, INTENT(in) :: DIMDEC(mask)
      INTEGER_T, INTENT(in) :: DIMDEC(maskcoef)
      INTEGER_T, INTENT(in) :: DIMDEC(maskSEM)
      INTEGER_T, INTENT(in) :: DIMDEC(levelPC)
      INTEGER_T, INTENT(in) :: DIMDEC(solxfab)
      INTEGER_T, INTENT(in) :: DIMDEC(solyfab)
      INTEGER_T, INTENT(in) :: DIMDEC(solzfab)
      INTEGER_T, INTENT(in) :: DIMDEC(cterm)
      INTEGER_T, INTENT(in) :: DIMDEC(pold)
      INTEGER_T, INTENT(in) :: DIMDEC(denold)
      INTEGER_T, INTENT(in) :: DIMDEC(ustar)
      INTEGER_T, INTENT(in) :: DIMDEC(mdotcell)
      INTEGER_T, INTENT(in) :: DIMDEC(maskdivres)
      INTEGER_T, INTENT(in) :: DIMDEC(maskres)

      REAL_T, INTENT(in), target :: xp(DIMV(xp),NCOMP_PEDGE)
      REAL_T, pointer :: xp_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: yp(DIMV(yp),NCOMP_PEDGE)
      REAL_T, pointer :: yp_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: zp(DIMV(zp),NCOMP_PEDGE)
      REAL_T, pointer :: zp_ptr(D_DECL(:,:,:),:)

      REAL_T, INTENT(in), target ::  xvel(DIMV(xvel),nsolve)
      REAL_T, pointer :: xvel_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target ::  yvel(DIMV(yvel),nsolve)
      REAL_T, pointer :: yvel_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target ::  zvel(DIMV(zvel),nsolve)
      REAL_T, pointer :: zvel_ptr(D_DECL(:,:,:),:)

      REAL_T, INTENT(in), target ::  xface(DIMV(xface),ncphys)
      REAL_T, pointer :: xface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target ::  yface(DIMV(yface),ncphys)
      REAL_T, pointer :: yface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target ::  zface(DIMV(zface),ncphys)
      REAL_T, pointer :: zface_ptr(D_DECL(:,:,:),:)

      REAL_T, INTENT(in), target ::  ax(DIMV(ax))
      REAL_T, INTENT(in), target ::  ay(DIMV(ay))
      REAL_T, INTENT(in), target ::  az(DIMV(az))
      REAL_T, pointer :: ax_ptr(D_DECL(:,:,:))
      REAL_T, pointer :: ay_ptr(D_DECL(:,:,:))
      REAL_T, pointer :: az_ptr(D_DECL(:,:,:))

      REAL_T, INTENT(in), target :: vol(DIMV(vol),1)
      REAL_T, pointer :: vol_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: rhs(DIMV(rhs),nsolve)
      REAL_T, pointer :: rhs_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: veldest(DIMV(veldest),ncomp_veldest)
      REAL_T, pointer :: veldest_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: dendest(DIMV(dendest),ncomp_dendest)
      REAL_T, pointer :: dendest_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: maskcoef(DIMV(maskcoef),1)
      REAL_T, pointer :: maskcoef_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: maskSEM(DIMV(maskSEM))
      REAL_T, pointer :: maskSEM_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: &
        levelPC(DIMV(levelPC),num_materials*(SDIM+1))
      REAL_T, pointer :: levelPC_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: solxfab(DIMV(solxfab),SDIM*nparts_def)
      REAL_T, INTENT(in), target :: solyfab(DIMV(solyfab),SDIM*nparts_def)
      REAL_T, INTENT(in), target :: solzfab(DIMV(solzfab),SDIM*nparts_def)
      REAL_T, pointer :: solxfab_ptr(D_DECL(:,:,:),:)
      REAL_T, pointer :: solyfab_ptr(D_DECL(:,:,:),:)
      REAL_T, pointer :: solzfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: cterm(DIMV(cterm),nsolve)
      REAL_T, pointer :: cterm_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: pold(DIMV(pold),nsolve)
      REAL_T, pointer :: pold_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: denold(DIMV(denold),ncomp_denold)
      REAL_T, pointer :: denold_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: ustar(DIMV(ustar),SDIM) 
      REAL_T, pointer :: ustar_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: &
       mdotcell(DIMV(mdotcell),nsolve) !VELADVECT_MF if OP_ISCHEME_CELL
      REAL_T, pointer :: mdotcell_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: &
       maskdivres(DIMV(maskdivres),ncomp_dendest) !DEN_RECON_MF OP_ISCHEME_CELL
      REAL_T, pointer :: maskdivres_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: maskres(DIMV(maskres))
      REAL_T, pointer :: maskres_ptr(D_DECL(:,:,:))

      REAL_T DXMAXLS,cutoff
      REAL_T Eforce_conservative
      REAL_T KE_diff

      INTEGER_T i,j,k
      INTEGER_T dir,side
      INTEGER_T veldir
      INTEGER_T im
      INTEGER_T sidecomp,ibase
      INTEGER_T ii,jj,kk
      INTEGER_T iface,jface,kface
      INTEGER_T, parameter :: nhalf=3
      REAL_T xsten(-nhalf:nhalf,SDIM)
      INTEGER_T imattype
      REAL_T AXL,AXR
      REAL_T AYL,AYR
      REAL_T AZL,AZR
      REAL_T VOLTERM,hx,RR
      REAL_T dencell
      REAL_T rho
      REAL_T TEMPERATURE,internal_e
      REAL_T NEW_TEMPERATURE
      REAL_T CC,CC_DUAL
      REAL_T MSKDV,MSKRES
      REAL_T MDOT,divu
      REAL_T local_rhs
      REAL_T local_POLD
      REAL_T local_POLD_DUAL

      REAL_T DIAG_REGULARIZE
      REAL_T uface(2)
      REAL_T ufacesolid(2)
      REAL_T aface(2)
      REAL_T pres_face(2)

       ! 0=no div(up)
       ! 1=div(up) ok
      INTEGER_T use_face_pres_cen
      INTEGER_T use_face_pres(2)  ! faces that are on either side of a cell.
      REAL_T ASIDE(2,ncphys)
      REAL_T mass_side(2)
      REAL_T masscell
      INTEGER_T local_maskSEM
      INTEGER_T stripstat
      INTEGER_T elemlo(3),elemhi(3)
      INTEGER_T ielem,jelem,kelem
      INTEGER_T scomp,scomp_bc,dcomp,ncomp ! in: mac_to_cell
      INTEGER_T ncomp_xvel
      INTEGER_T ncomp_cterm
      INTEGER_T is_rigid_near
      REAL_T LStest(num_materials)
      INTEGER_T velcomp
      INTEGER_T partid
      INTEGER_T partid_ghost
      INTEGER_T nparts_temp,im_solid
      INTEGER_T cell_velocity_override
      
      REAL_T xclamped(SDIM)
      REAL_T LS_clamped
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped
      INTEGER_T prescribed_flag

      REAL_T local_div_val

      REAL_T massfrac_parm(num_species_var+1)
      INTEGER_T ispec

      vol_ptr=>vol
      rhs_ptr=>rhs
      veldest_ptr=>veldest
      dendest_ptr=>dendest
      cterm_ptr=>cterm
      ustar_ptr=>ustar
      xface_ptr=>xface
      yface_ptr=>yface
      zface_ptr=>zface
      xp_ptr=>xp
      yp_ptr=>yp
      zp_ptr=>zp
      xvel_ptr=>xvel
      yvel_ptr=>yvel
      zvel_ptr=>zvel
      maskSEM_ptr=>maskSEM
      maskcoef_ptr=>maskcoef
      mask_ptr=>mask
      maskdivres_ptr=>maskdivres
      maskres_ptr=>maskres
      mdotcell_ptr=>mdotcell
      pold_ptr=>pold
      denold_ptr=>denold
      ax_ptr=>ax
      ay_ptr=>ay
      az_ptr=>az
      solxfab_ptr=>solxfab
      solyfab_ptr=>solyfab
      solzfab_ptr=>solzfab

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_mac_to_cell"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_mac_to_cell"
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
       print *,"num_divu_outer_sweeps invalid fort_mac_to_cell"
       stop
      endif
      if ((divu_outer_sweeps.lt.0).or. &
          (divu_outer_sweeps.ge.num_divu_outer_sweeps)) then
       print *,"divu_outer_sweeps invalid fort_mac_to_cell"
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

      if ((enable_spectral.ge.0).and. &
          (enable_spectral.le.1)) then
       ! do nothing
      else
       print *,"enable_spectral invalid mac_to_cell"
       stop
      endif

      do im=1,num_materials
       imattype=fort_material_type(im)
       if (imattype.eq.999) then
        ! do nothing
       else if ((imattype.ge.0).and. &
                (imattype.le.MAX_NUM_EOS)) then
        ! do nothing
       else
        print *,"imattype invalid fort_mac_to_cell"
        stop
       endif
      enddo ! im=1..num_materials

      if (project_option_is_validF(project_option).eq.1) then
       ! do nothing
      else
       print *,"project_option invalid"
       stop
      endif

      ! mac -> cell in solver (init_divup_cell_vel_cell) or VELMAC_TO_CELL
      if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then ! velocity
       if (ncomp_veldest.ge.SDIM) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.ge.1) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncphys.ne.FACECOMP_NCOMP) then
        print *,"ncphys invalid"
        stop
       endif
       if (nsolve.ne.1) then
        print *,"nsolve invalid 2"
        stop
       endif

       if (ncomp_denold.eq.nsolve) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

       ! div(up)
      else if (operation_flag.eq.OP_VEL_DIVUP_TO_CELL) then 
 
       if (ncomp_veldest.ge. &
           SDIM+num_state_material*num_materials) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.ge.num_state_material*num_materials) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif
       if (ncphys.ne.FACECOMP_NCOMP) then
        print *,"ncphys invalid"
        stop
       endif
       if (nsolve.ne.1) then
        print *,"nsolve invalid 2"
        stop
       endif

       if (ncomp_denold.eq.nsolve) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

       if (enable_spectral.eq.0) then
        ! do nothing
       else
        print *,"must have enable_spectral==0 if OP_VEL_DIVUP_TO_CELL"
        stop
       endif

      else if (operation_flag.eq.OP_RHS_CELL) then ! rhs of solver

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

       if (project_option_is_validF(project_option).eq.1) then

        if (ncphys.ne.FACECOMP_NCOMP) then
         print *,"ncphys invalid"
         stop
        endif

       else
        print *,"project_option invalid"
        stop
       endif

       if ((nsolve.ne.1).and. &
           (nsolve.ne.SDIM)) then
        print *,"nsolve invalid 2"
        stop
       endif

       if (ncomp_denold.eq.nsolve) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

      else if (operation_flag.eq.OP_DIV_CELL) then ! divergence

       if (ncomp_veldest.eq.1) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.eq.1) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
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

      else if (operation_flag.eq.OP_ISCHEME_CELL) then ! advection

       if (ncomp_veldest.eq.STATE_NCOMP) then
        ! do nothing
       else
        print *,"ncomp_veldest invalid"
        stop
       endif
       if (ncomp_dendest.eq.ncomp_veldest-STATECOMP_STATES) then
        ! do nothing
       else
        print *,"ncomp_dendest invalid"
        stop
       endif

       if ((nsolve.ne.NFLUXSEM).or. &
           (ncphys.ne.NFLUXSEM)) then
        print *,"nsolve or ncphys invalid"
        stop
       endif

       if (ncomp_denold.eq.num_materials*num_state_material) then
        ! do nothing
       else
        print *,"ncomp_denold invalid"
        stop
       endif

      else
       print *,"operation_flag invalid6: ",operation_flag
       stop
      endif

      if (operation_flag.eq.OP_RHS_CELL) then  ! rhs for solver

       if (energyflag.eq.SUB_OP_DEFAULT) then
        ! do nothing
       else
        print *,"energyflag invalid OP_RHS_CELL"
        stop
       endif
       if ((homflag.ge.0).and.(homflag.le.4)) then
        ! do nothing
       else
        print *,"homflag invalid in mac to cell homflag=",homflag
        stop
       endif

      else if (operation_flag.eq.OP_DIV_CELL) then ! divergence

       if (energyflag.eq.SUB_OP_DEFAULT) then
        ! do nothing
       else
        print *,"energyflag invalid OP_DIV_CELL"
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
      else if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then ! velocity

       if (homflag.eq.0) then
        ! do nothing
       else
        print *,"homflag invalid"
        stop
       endif
       if ((energyflag.ne.SUB_OP_THERMAL_DIVUP_NULL).and. &
           (energyflag.ne.SUB_OP_THERMAL_DIVUP_OK)) then
        print *,"energyflag invalid OP_VEL_MAC_TO_CELL "
        stop
       endif
       if (nsolve.ne.1) then
        print *,"nsolve invalid5"
        stop
       endif

        ! copy u^{mac->cell} to u^{cell}, div(up),
      else if (operation_flag.eq.OP_VEL_DIVUP_TO_CELL) then 
       if (homflag.ne.0) then
        print *,"homflag invalid"
        stop
       endif
       if ((energyflag.ne.SUB_OP_THERMAL_DIVUP_NULL).and. & 
           (energyflag.ne.SUB_OP_THERMAL_DIVUP_OK)) then 
        print *,"energyflag invalid OP_VEL_DIVUP_TO_CELL"
        stop
       endif
       if (nsolve.ne.1) then
        print *,"nsolve invalid6"
        stop
       endif

      else if (operation_flag.eq.OP_ISCHEME_CELL) then ! advection

        ! "source_term" 
       if ((homflag.ne.SUB_OP_SDC_ISCHEME).and. &
           (homflag.ne.SUB_OP_SDC_LOW_TIME)) then
        print *,"homflag invalid in mac to cell homflag=",homflag
        stop
       endif
        ! "advect_iter"
       if ((energyflag.ne.SUB_OP_ISCHEME_PREDICT).and. &
           (energyflag.ne.SUB_OP_ISCHEME_CORRECT)) then
        print *,"energyflag invalid OP_ISCHEME_CELL"
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

      if (project_option_is_validF(project_option).eq.1) then
       ! do nothing
      else
       print *,"project_option invalid"
       stop
      endif

      do im=1,num_materials
       if (fort_denconst(im).le.zero) then
        print *,"denconst invalid"
        stop
       endif
      enddo

      call checkbound_array(fablo,fabhi,xp_ptr,0,0)
      call checkbound_array(fablo,fabhi,yp_ptr,0,1)
      call checkbound_array(fablo,fabhi,zp_ptr,0,SDIM-1)

      call checkbound_array(fablo,fabhi,xvel_ptr,0,0)
      call checkbound_array(fablo,fabhi,yvel_ptr,0,1)
      call checkbound_array(fablo,fabhi,zvel_ptr,0,SDIM-1)

      call checkbound_array1(fablo,fabhi,ax_ptr,0,0)
      call checkbound_array1(fablo,fabhi,ay_ptr,0,1)
      call checkbound_array1(fablo,fabhi,az_ptr,0,SDIM-1)

      call checkbound_array(fablo,fabhi,vol_ptr,0,-1)
      call checkbound_array(fablo,fabhi,rhs_ptr,0,-1)
      call checkbound_array(fablo,fabhi,veldest_ptr,0,-1)
      call checkbound_array(fablo,fabhi,dendest_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      call checkbound_array(fablo,fabhi,maskcoef_ptr,1,-1)

      levelPC_ptr=>levelPC
      call checkbound_array(fablo,fabhi,levelPC_ptr,1,-1)

      call checkbound_array(fablo,fabhi,cterm_ptr,0,-1)
      call checkbound_array(fablo,fabhi,pold_ptr,0,-1)
      call checkbound_array(fablo,fabhi,denold_ptr,0,-1)
      call checkbound_array(fablo,fabhi,ustar_ptr,0,-1)

      call checkbound_array(fablo,fabhi,mdotcell_ptr,0,-1)
      call checkbound_array(fablo,fabhi,maskdivres_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,maskres_ptr,0,-1)


      if (project_option_is_validF(project_option).eq.1) then

         !local_face_var_mf
       call checkbound_array(fablo,fabhi,xface_ptr,0,0)
       call checkbound_array(fablo,fabhi,yface_ptr,0,1)
       call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1)
       call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)
         !local_fsi_ghost_mac_mf+dir
       call checkbound_array(fablo,fabhi,solxfab_ptr,0,0)
       call checkbound_array(fablo,fabhi,solyfab_ptr,0,1)
       call checkbound_array(fablo,fabhi,solzfab_ptr,0,SDIM-1)

      else
       print *,"project_option invalid"
       stop
      endif

        ! max(dx,dy,dz) if XYZ or R - Theta - Z
        ! (get_dxmax(): max(dr,r probhix,dz) if R - Theta - Z)
      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=DXMAXLS

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridsten_level(xsten,i,j,k,level,nhalf)

       do dir=1,SDIM
        xclamped(dir)=xsten(0,dir)
       enddo

       AXL=ax(D_DECL(i,j,k))
       AXR=ax(D_DECL(i+1,j,k))
       AYL=ay(D_DECL(i,j,k))
       AYR=ay(D_DECL(i,j+1,k))
       AZL=az(D_DECL(i,j,k))
       AZR=az(D_DECL(i,j,k+1))
       VOLTERM=vol(D_DECL(i,j,k),1)
       if (VOLTERM.gt.zero) then
        ! do nothing
       else
        print *,"VOLTERM invalid"
        stop
       endif
       if ((AXL.ge.zero).and. &
           (AXR.ge.zero).and. &
           (AYL.ge.zero).and. &
           (AYR.ge.zero).and. &
           (AZL.ge.zero).and. &
           (AZR.ge.zero)) then
        ! do nothing
       else
        print *,"AX,AY or AZ invalid"
        print *,"AXL,AXR,AYL,AYR,AZL,AZR ", &
                AXL,AXR,AYL,AYR,AZL,AZR
        print *,"i,j,k ",i,j,k
        print *,"xsten(0,?) : ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
        print *,"nhalf= ",nhalf
        print *,"level= ",level
        print *,"finest_level= ",finest_level
        print *,"tilelo = ",tilelo(1),tilelo(2),tilelo(SDIM)
        print *,"tilehi = ",tilehi(1),tilehi(2),tilehi(SDIM)
        print *,"fablo = ",fablo(1),fablo(2),fablo(SDIM)
        print *,"fabhi = ",fabhi(1),fabhi(2),fabhi(SDIM)
        print *,"growlo = ",growlo(1),growlo(2),growlo(SDIM)
        print *,"growhi = ",growhi(1),growhi(2),growhi(SDIM)
        print *,"homflag=",homflag
        print *,"project_option=",project_option
        print *,"operation_flag=",operation_flag
        print *,"energyflag=",energyflag
        print *,"nsolve=",nsolve
        stop
       endif

       if (operation_flag.eq.OP_DIV_CELL) then ! DIV

        divu= &
         AXR*xvel(D_DECL(i+1,j,k),1)-  &
         AXL*xvel(D_DECL(i,j,k),1)+ &
         AYR*yvel(D_DECL(i,j+1,k),1)-  &
         AYL*yvel(D_DECL(i,j,k),1)
        if (SDIM.eq.3) then
         divu=divu+ &
          AZR*zvel(D_DECL(i,j,k+1),1)-  &
          AZL*zvel(D_DECL(i,j,k),1)
        endif
        divu=divu/VOLTERM

        MSKDV=maskdivres(D_DECL(i,j,k),1)
        if (MSKDV.eq.one) then
         ! do nothing
        else if (MSKDV.eq.zero) then
         divu=zero
        else
         print *,"MSKDV invalid"
         stop
        endif

        rhs(D_DECL(i,j,k),1)=divu

       else if (operation_flag.eq.OP_RHS_CELL) then ! RHS

        ! (cterm)*p-vol grad dot grad p/rho=-vol div u/dt + mdot +
        !    cterm * p^adv 
        ! cterm=vol/(rho c^2 dt*dt)

        if (maskcoef(D_DECL(i,j,k),1).eq.one) then ! not covered

         do veldir=1,nsolve
    
          CC=cterm(D_DECL(i,j,k),veldir)
          CC_DUAL=veldest(D_DECL(i,j,k),veldir)
          if (CC_DUAL.eq.CC) then
           ! do nothing
          else
           print *,"CC_DUAL invalid"
           stop
          endif
          MSKDV=maskdivres(D_DECL(i,j,k),1)
          MSKRES=maskres(D_DECL(i,j,k))
          MDOT=mdotcell(D_DECL(i,j,k),veldir)

          DIAG_REGULARIZE=denold(D_DECL(i,j,k),veldir)

          if (DIAG_REGULARIZE.gt.zero) then
           ! check nothing
          else
           print *,"DIAG_REGULARIZE invalid"
           stop
          endif

           !SOLVETYPE_PRES, 
           !SOLVETYPE_PRESSTATIC, 
           !SOLVETYPE_PRESGRAVITY, 
           !SOLVETYPE_INITPROJ, 
           !SOLVETYPE_PRESEXTRAP
          if (project_option_singular_possibleF(project_option).eq.1) then

           if (MDOT.eq.zero) then
            ! check nothing
           else if (MDOT.ne.zero) then
            if (MSKRES.eq.zero) then
             print *,"cannot have MDOT<>0 and MSKRES==0"
             stop
            else if (MSKRES.ne.zero) then
             ! do nothing
            else
             print *,"MSKRES is NaN; fort_mac_to_cell"
             stop
            endif
           else
            print *,"MDOT is NaN; fort_mac_to_cell"
            stop
           endif
           if ((CC.ge.zero).and.(CC_DUAL.ge.zero)) then
            ! do nothing
           else
            print *,"CC or CC_DUAL invalid; fort_mac_to_cell"
            stop
           endif
          else if (project_option_singular_possibleF(project_option).eq.0) then
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

            ! AXR,AXL,AYR,AYL,AZR,AZL are face areas.
           divu= &
            AXR*xvel(D_DECL(i+1,j,k),veldir)-  &
            AXL*xvel(D_DECL(i,j,k),veldir)+ &
            AYR*yvel(D_DECL(i,j+1,k),veldir)-  &
            AYL*yvel(D_DECL(i,j,k),veldir)
           if (SDIM.eq.3) then
            divu=divu+ &
             AZR*zvel(D_DECL(i,j,k+1),veldir)-  &
             AZL*zvel(D_DECL(i,j,k),veldir)
           endif
          else
           print *,"maskdivres invalid" 
           stop
          endif

          local_div_val=divu/VOLTERM

          if ((local_div_val.ge.zero).or. &
              (local_div_val.le.zero)) then
           ! do nothing
          else
           print *,"local_div_val=NaN  ",local_div_val
           print *,"divu invalid ",divu
           print *,"VOLTERM ",VOLTERM
           print *,"operation_flag ",operation_flag
           print *,"project_option ",project_option
           print *,"AXL,AXR,AYL,AYR,AZL,AZR ", &
                AXL,AXR,AYL,AYR,AZL,AZR
           print *,"i,j,k ",i,j,k
           print *,"xsten(0,?) : ",xsten(0,1),xsten(0,2),xsten(0,SDIM)
           print *,"nhalf= ",nhalf
           print *,"level= ",level
           print *,"finest_level= ",finest_level
           print *,"tilelo = ",tilelo(1),tilelo(2),tilelo(SDIM)
           print *,"tilehi = ",tilehi(1),tilehi(2),tilehi(SDIM)
           print *,"fablo = ",fablo(1),fablo(2),fablo(SDIM)
           print *,"fabhi = ",fabhi(1),fabhi(2),fabhi(SDIM)
           print *,"growlo = ",growlo(1),growlo(2),growlo(SDIM)
           print *,"growhi = ",growhi(1),growhi(2),growhi(SDIM)
           print *,"homflag=",homflag
           print *,"energyflag=",energyflag
           print *,"nsolve=",nsolve
           print *,"veldir ",veldir
           print *,"xvel(D_DECL(i+1,j,k),veldir) ", &
                   xvel(D_DECL(i+1,j,k),veldir)
           print *,"yvel(D_DECL(i,j+1,k),veldir) ", &
                   yvel(D_DECL(i,j+1,k),veldir)
           stop
          endif

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

         enddo ! veldir=1..nsolve

        else if (maskcoef(D_DECL(i,j,k),1).eq.zero) then
         ! do nothing (covered)
        else 
         print *,"mask invalid"
         stop
        endif

       ! mac -> cell in solver (init_divup_cell_vel_cell) or VELMAC_TO_CELL
       else if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then ! velocity

         ! LS>0 if clamped
        call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
          vel_clamped,temperature_clamped,prescribed_flag,dx)
        if (project_option_is_static(project_option).eq.1) then
         do dir=1,SDIM
          vel_clamped(dir)=zero
         enddo
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

         enddo ! side=1..2

         partid=-1
         partid_ghost=0
         nparts_temp=0
         im_solid=0
         do im=1,num_materials
          LStest(im)=levelPC(D_DECL(i,j,k),im)
          if (is_lag_part(im).eq.1) then

           if (is_rigid(im).eq.1) then
            if (is_prescribed(im).eq.1) then
             if (im_solid.eq.0) then
              im_solid=im
              partid=nparts_temp
             else if ((im_solid.ge.1).and.(im_solid.le.num_materials)) then
              if (LStest(im).gt.LStest(im_solid)) then
               im_solid=im
               partid=nparts_temp
              endif
             else
              print *,"im_solid invalid 5"
              stop
             endif
            else if (is_prescribed(im).eq.0) then
             ! do nothing
            else
             print *,"is_prescribed(im) invalid"
             stop
            endif
           else if (is_rigid(im).eq.0) then
            ! do nothing
           else
            print *,"is_rigid invalid LEVELSET_3D.F90"
            stop
           endif
           nparts_temp=nparts_temp+1

          else if (is_lag_part(im).eq.0) then
           if (is_rigid(im).eq.0) then
            ! do nothing
           else
            print *,"is_rigid invalid LEVELSET_3D.F90"
            stop
           endif
          else
           print *,"is_lag_part invalid"
           stop
          endif

         enddo ! im=1..num_materials

         if (nparts_temp.ne.nparts) then
          print *,"nparts_temp invalid"
          stop
         endif
         if ((partid.ge.0).and.(partid.lt.nparts)) then
          if (im_solid_map(partid+1)+1.ne.im_solid) then
           print *,"im_solid invalid 6"
           stop
          endif
          if (LStest(im_solid).ge.zero) then
           ! do nothing
          else if (LStest(im_solid).lt.zero) then
           partid=-1
          else
           print *,"LStest(im_solid) is NaN"
           stop
          endif
         else if (partid.eq.-1) then
          ! do nothing
         else
          print *,"partid invalid"
          stop
         endif
         if (partid.eq.-1) then
          partid_ghost=0
         else if ((partid.ge.0).and.(partid.lt.nparts)) then
          partid_ghost=partid
         else
          print *,"partid invalid"
          stop
         endif
 
         velcomp=dir+1

         cell_velocity_override=0
         if ((partid.ge.0).and.(partid.lt.nparts)) then
          cell_velocity_override=1
         else if (partid.eq.-1) then
          ! do nothing
         else
          print *,"partid invalid"
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
           if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then
            uface(side)=xvel(D_DECL(iface,jface,kface),1)
            if (project_option_is_static(project_option).eq.1) then
             ufacesolid(side)=zero
            else
             ufacesolid(side)=solxfab(D_DECL(iface,jface,kface), &
                   partid_ghost*SDIM+dir+1)
            endif
           else
            print *,"operation_flag invalid"
            stop
           endif
           if (SDIM.eq.2) then
            if (levelrz.eq.COORDSYS_CARTESIAN) then
             ! do nothing
            else if (levelrz.eq.COORDSYS_RZ) then
             if (iface.eq.0) then
              if ((side.eq.1).and.(i.eq.0)) then
               if (xsten(-2,1).lt.zero) then
                uface(side)=zero
                ufacesolid(side)=zero
               else if (xsten(-2,1).gt.zero) then
                ! do nothing
               else
                print *,"xsten is NaN"
                stop
               endif
              else 
               print *,"side or i invalid"
               stop
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
            else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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
           if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then  
            uface(side)=yvel(D_DECL(iface,jface,kface),1)
            if (project_option_is_static(project_option).eq.1) then
             ufacesolid(side)=zero
            else
             ufacesolid(side)=solyfab(D_DECL(iface,jface,kface), &
                   partid_ghost*SDIM+dir+1)
            endif
           else
            print *,"operation_flag invalid"
            stop
           endif
          else if ((dir.eq.2).and.(SDIM.eq.3)) then
           if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then 
            uface(side)=zvel(D_DECL(iface,jface,kface),1)
            if (project_option_is_static(project_option).eq.1) then
             ufacesolid(side)=zero
            else
             ufacesolid(side)=solzfab(D_DECL(iface,jface,kface), &
                   partid_ghost*SDIM+dir+1)
            endif
           else
            print *,"operation_flag invalid"
            stop
           endif
          else
           print *,"dir invalid mac to cell 3"
           stop
          endif

          mass_side(side)=zero
          do im=1,num_materials 
           if (side.eq.1) then  ! left half of cell
            sidecomp=FACECOMP_MASSFACE+2*(im-1)+2
           else if (side.eq.2) then ! right half of cell
            sidecomp=FACECOMP_MASSFACE+2*(im-1)+1
           else
            print *,"side invalid"
            stop
           endif
           mass_side(side)=mass_side(side)+ASIDE(side,sidecomp) 
          enddo ! im=1..num_materials

         enddo ! side=1..2

         masscell=mass_side(1)+mass_side(2)

         if ((mass_side(1).gt.zero).and. &
             (mass_side(2).gt.zero).and. &
             (masscell.gt.zero)) then
          ! do nothing
         else
          print *,"mass invalid"
          stop
         endif

         if (LS_clamped.ge.zero) then
          veldest(D_DECL(i,j,k),velcomp)=vel_clamped(dir+1)
         else if (LS_clamped.lt.zero) then

          if (cell_velocity_override.eq.1) then
           veldest(D_DECL(i,j,k),velcomp)= &
            (mass_side(1)*ufacesolid(1)+ &
             mass_side(2)*ufacesolid(2))/masscell
          else if (cell_velocity_override.eq.0) then
           veldest(D_DECL(i,j,k),velcomp)= &
             (mass_side(1)*uface(1)+ &
              mass_side(2)*uface(2))/masscell
          else
           print *,"cell_velocity_override invalid"
           stop
          endif

         else
          print *,"LS_clamped is NaN"
          stop
         endif

         if (1.eq.0) then
          print *,"AFTER UPDATING CELL FROM FACE (velocity is scaled)"
          print *,"i,j,k,dir,velcomp ",i,j,k,dir,velcomp
          print *,"snew= ",veldest(D_DECL(i,j,k),velcomp)
         endif

        enddo  ! dir=0..sdim-1

        ! if use_face_pres.eq.1 then div(up) might be ok
        ! this routine: fort_mac_to_cell
       else if (operation_flag.eq.OP_VEL_DIVUP_TO_CELL) then ! div(up)

         ! LS>0 if clamped
        call SUB_clamped_LS(xclamped,cur_time,LS_clamped, &
          vel_clamped,temperature_clamped,prescribed_flag,dx)

        use_face_pres_cen=1

        if (LS_clamped.ge.zero) then
         use_face_pres_cen=0
        else if (LS_clamped.lt.zero) then
         ! do nothing
        else
         print *,"LS_clamped invalid"
         stop
        endif

         ! note, in fort_build_conserve, if 
         ! is_compressible_mat==0,
         ! then 
         ! (1) (rho T) is advected instead of (rho cv T + rho u dot u/2)
         ! (2) rho_t + u dot grad rho=0 instead of
         !     rho_t + div(rho u)=0
         ! 
        is_rigid_near=0
        do im=1,num_materials
         LStest(im)=levelPC(D_DECL(i,j,k),im)
         if (is_rigid(im).eq.1) then
          if (LStest(im).ge.-DXMAXLS) then
           is_rigid_near=1
          else if (LStest(im).le.-DXMAXLS) then
           ! do nothing
          else
           print *,"LStest(im) is NaN"
           stop
          endif 
         else if (is_rigid(im).eq.0) then
          ! do nothing
         else
          print *,"is_rigid(im) invalid"
          stop
         endif
        enddo ! im=1..num_materials

        if (is_rigid_near.eq.1) then
         use_face_pres_cen=0
        else if (is_rigid_near.eq.0) then
         do im=1,num_materials
          if (LStest(im).ge.-DXMAXLS) then
           if (is_compressible_mat(im).eq.0) then
            use_face_pres_cen=0
           else if (is_compressible_mat(im).eq.1) then
            ! do nothing
           else
            print *,"is_compressible_mat invalid"
            stop
           endif
          else if (LStest(im).le.-DXMAXLS) then
           ! do nothing
          else
           print *,"LStest(im) is NaN(1)fort_mac_to_cell OP_VEL_DIVUP_TO_CELL"
           stop
          endif 
         enddo ! im=1..num_materials
        else
         print *,"is_rigid_near inv(1)fort_mac_to_cell OP_VEL_DIVUP_TO_CELL"
         stop
        endif 

        Eforce_conservative=zero

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
          print *,"dir out of range in fort_mac_to_cell"
          stop
         endif

         hx=xsten(1,dir+1)-xsten(-1,dir+1)
         RR=one
         if ((levelrz.eq.COORDSYS_CARTESIAN).or.(levelrz.eq.COORDSYS_RZ)) then
          RR=one
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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
         if (hx.gt.zero) then
          ! do nothing
         else
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

          if (dir.eq.0) then
           aface(side)=ax(D_DECL(iface,jface,kface))
           do im=1,ncphys
            ASIDE(side,im)=xface(D_DECL(iface,jface,kface),im)
           enddo
           uface(side)=xvel(D_DECL(iface,jface,kface),1)
           pres_face(side)= &
                 xp(D_DECL(iface,jface,kface),PRESSURE_PEDGE+1)
           use_face_pres(side)= &
                 NINT(xp(D_DECL(iface,jface,kface),VALID_PEDGE+1))
          else if (dir.eq.1) then
           aface(side)=ay(D_DECL(iface,jface,kface))
           do im=1,ncphys
            ASIDE(side,im)=yface(D_DECL(iface,jface,kface),im)
           enddo
           uface(side)=yvel(D_DECL(iface,jface,kface),1)
           pres_face(side)= &
                 yp(D_DECL(iface,jface,kface),PRESSURE_PEDGE+1)
           use_face_pres(side)= &
                 NINT(yp(D_DECL(iface,jface,kface),VALID_PEDGE+1))
          else if ((dir.eq.2).and.(SDIM.eq.3)) then
           aface(side)=az(D_DECL(iface,jface,kface))
           do im=1,ncphys
            ASIDE(side,im)=zface(D_DECL(iface,jface,kface),im)
           enddo
           uface(side)=zvel(D_DECL(iface,jface,kface),1)
           pres_face(side)= &
                 zp(D_DECL(iface,jface,kface),PRESSURE_PEDGE+1)
           use_face_pres(side)= &
                 NINT(zp(D_DECL(iface,jface,kface),VALID_PEDGE+1))
          else
           print *,"dir invalid mac to cell 3"
           stop
          endif

          mass_side(side)=zero
          do im=1,num_materials 
           if (side.eq.1) then  ! left half of cell
            sidecomp=FACECOMP_MASSFACE+2*(im-1)+2
           else if (side.eq.2) then ! right half of cell
            sidecomp=FACECOMP_MASSFACE+2*(im-1)+1
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
         if (dencell.gt.zero) then
          ! do nothing
         else
          print *,"dencell invalid"
          stop
         endif

         if (1.eq.0) then
          if (dir.eq.SDIM-1) then
           print *,"i,j,k,dir,dencell ",i,j,k,dir,dencell
          endif
         endif

           ! use_face_pres.eq.1 => div(up) ok
         if ((use_face_pres(1).eq.1).and. &
             (use_face_pres(2).eq.1)) then
          ! do nothing
         else if ((use_face_pres(1).eq.0).or. &
                  (use_face_pres(2).eq.0)) then
          use_face_pres_cen=0
         else
          print *,"use_face_pres invalid"
          stop
         endif
        
         if (energyflag.eq.SUB_OP_THERMAL_DIVUP_NULL) then
          ! do nothing
         else if (energyflag.eq.SUB_OP_THERMAL_DIVUP_OK) then

          Eforce_conservative=Eforce_conservative- &
            dt*(aface(2)*uface(2)*pres_face(2)- &
                aface(1)*uface(1)*pres_face(1))/ &
            (dencell*VOLTERM)

         else
          print *,"energyflag invalid OP_VEL_DIVUP_TO_CELL"
          stop
         endif

        enddo ! dir=0..sdim-1 (operation_flag.eq.OP_VEL_DIVUP_TO_CELL) div(up)

        if ((use_face_pres_cen.eq.0).or. &
            (use_face_pres_cen.eq.1)) then 
         ! do nothing
        else
         print *,"use_face_pres_cen bust"
         stop
        endif

        rhs(D_DECL(i,j,k),1)=zero

         ! update the total energy in compressible cells (regular project)
        if (project_option.eq.SOLVETYPE_PRES) then

         if (use_face_pres_cen.eq.0) then

          Eforce_conservative=zero
          rhs(D_DECL(i,j,k),1)=zero

         else if (use_face_pres_cen.eq.1) then

          rhs(D_DECL(i,j,k),1)=Eforce_conservative

           ! update the temperature
          if (energyflag.eq.SUB_OP_THERMAL_DIVUP_OK) then 

           do im=1,num_materials

            KE_diff=zero
            do velcomp=1,SDIM
             KE_diff=KE_diff+ &
              half*ustar(D_DECL(i,j,k),velcomp)**2- &
              half*veldest(D_DECL(i,j,k),velcomp)**2
            enddo ! velcomp=1..sdim

            if (LStest(im).ge.-DXMAXLS) then

             if (is_compressible_mat(im).eq.1) then

              imattype=fort_material_type(im)

              ibase=(im-1)*num_state_material
           
               ! dendest is Snewfab.dataPtr(STATECOMP_STATES) 

              rho=dendest(D_DECL(i,j,k),ibase+ENUM_DENVAR+1)
              TEMPERATURE=dendest(D_DECL(i,j,k),ibase+ENUM_TEMPERATUREVAR+1)

              if (rho.gt.zero) then
               ! do nothing
              else
               print *,"density underflow"
               print *,"i,j,k,im,rho ",i,j,k,im,rho
               stop
              endif
              if (TEMPERATURE.gt.zero) then
               ! do nothing
              else
               print *,"temperature underflow"
               stop
              endif
              call init_massfrac_parm(rho,massfrac_parm,im)
              do ispec=1,num_species_var
               massfrac_parm(ispec)= &
                 dendest(D_DECL(i,j,k),ibase+ENUM_SPECIESVAR+ispec)
               if (massfrac_parm(ispec).ge.zero) then
                ! do nothing
               else
                print *,"massfrac_parm invalid"
                stop
               endif
              enddo ! ispec=1..num_species_var

              call INTERNAL_material(rho,massfrac_parm, &
               TEMPERATURE,internal_e, &
               imattype,im)

              ! sanity checks
              if (internal_e.gt.zero) then
               call TEMPERATURE_material(rho,massfrac_parm, &
                NEW_TEMPERATURE, &
                internal_e,imattype,im)
               if (abs(TEMPERATURE-NEW_TEMPERATURE).le. &
                   1.0D-3*TEMPERATURE) then
                ! do nothing 
               else
                print *,"T(rho,e) and e(rho,T) are not inverses"
                stop
               endif
              else
               print *,"internal_e must be positive"
               stop
              endif

              internal_e=internal_e+KE_diff+Eforce_conservative

              if (internal_e.le.zero) then
               NEW_TEMPERATURE=TEMPERATURE
              else if (internal_e.gt.zero) then
               call TEMPERATURE_material(rho,massfrac_parm, &
                NEW_TEMPERATURE, &
                internal_e,imattype,im)
              else
               print *,"internal_e bust"
               stop
              endif

              if (NEW_TEMPERATURE.gt.zero) then
               dendest(D_DECL(i,j,k),ibase+ENUM_TEMPERATUREVAR+1)= &
                       NEW_TEMPERATURE
              else
               print *,"NEW_TEMPERATURE must be positive"
               stop
              endif

             else if (is_compressible_mat(im).eq.0) then
              ! do nothing
             else
              print *,"is_compressible_mat(im) invalid"
              stop
             endif
            else if (LStest(im).le.-DXMAXLS) then
             ! do nothing
            else
             print *,"LStest(im) is NaN"
             stop
            endif 
           enddo ! im=1..num_materials

           !do not update temperature
          else if (energyflag.eq.SUB_OP_THERMAL_DIVUP_NULL) then
           ! do nothing
          else
           print *,"energyflag invalid OP_VEL_DIVUP_TO_CELL" 
           stop
          endif

         else 
          print *,"use_face_pres_cen invalid(2)"
          print *,"use_face_pres_cen=",use_face_pres_cen
          stop
         endif

        else if (project_option.eq.SOLVETYPE_INITPROJ) then
         if (energyflag.eq.SUB_OP_THERMAL_DIVUP_NULL) then
          ! do nothing if initial project
         else
          print *,"expecting (energyflag.eq.SUB_OP_THERMAL_DIVUP_NULL)"
          stop
         endif
        else if (project_option.eq.SOLVETYPE_PRESGRAVITY) then
         if (energyflag.eq.SUB_OP_THERMAL_DIVUP_NULL) then
          ! do nothing if solvetype_presgravity
         else
          print *,"expecting (energyflag.eq.SUB_OP_THERMAL_DIVUP_NULL)mac_cell"
          stop
         endif
        else if (project_option.eq.SOLVETYPE_PRESSTATIC) then
         if (energyflag.eq.SUB_OP_THERMAL_DIVUP_NULL) then
          ! do nothing if solvetype_presstatic
         else
          print *,"expecting (energyflag.eq.SUB_OP_THERMAL_DIVUP_NULL)mac_cell"
          stop
         endif
        else
         print *,"project_option invalid fort_mac_to_cell 5"
         print *,"operation_flag.eq.OP_VEL_DIVUP_TO_CELL"
         stop
        endif

       else if (operation_flag.eq.OP_ISCHEME_CELL) then ! advection
        ! low order approximation: CISL or sem_mac_to_cell
        ! high order approximation: sem_mac_to_cell
       else
        print *,"operation_flag invalid8"
        stop
       endif

      enddo
      enddo
      enddo

      if (ns_time_order.eq.1) then
       if (enable_spectral.eq.0) then
        ! do nothing
       else
        print *,"enable_spectral and ns_time_order are not consistent(1)"
        print *,"ns_time_order, enable_spectral ",ns_time_order, &
                enable_spectral
        stop
       endif
      else if ((ns_time_order.ge.2).and.(ns_time_order.le.32)) then

       if (enable_spectral.eq.1) then

        if (operation_flag.eq.OP_VEL_DIVUP_TO_CELL) then
         print *,"expecting enable_spectral=0"
         stop
        endif

       else if (enable_spectral.eq.0) then

        if (operation_flag.eq.OP_VEL_DIVUP_TO_CELL) then
         ! do nothing
        else if (operation_flag.eq.OP_RHS_CELL) then
         ! do nothing
        else
         print *,"expecting enable_spectral==1"
         stop
        endif

       else
        print *,"enable_spectral and ns_time_order are not consistent(2)"
        print *,"ns_time_order, enable_spectral ",ns_time_order, &
                enable_spectral
        stop
       endif
      else
       print *,"ns_time_order invalid"
       stop
      endif

       ! in: fort_mac_to_cell
      if (enable_spectral.eq.1) then

        ! The ISCHEME divergence is taken care of in SEM_MAC_TO_CELL
        ! regardless of "bfact".
       if ((bfact.ge.2).or. &
           ((bfact.eq.1).and.(operation_flag.eq.OP_ISCHEME_CELL))) then

        call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

         local_maskSEM=NINT(maskSEM(D_DECL(i,j,k)))

         if ((local_maskSEM.ge.1).and. &
             (local_maskSEM.le.num_materials)) then

          call strip_status(i,j,k,bfact,stripstat)

          if (stripstat.eq.1) then

           do dir=0,SDIM-1

            call elementbox(i,j,k,bfact,dir,elemlo,elemhi)
            do ielem=elemlo(1),elemhi(1)
            do jelem=elemlo(2),elemhi(2)
            do kelem=elemlo(3),elemhi(3)
          
             if (operation_flag.eq.OP_RHS_CELL) then ! RHS for solver
              scomp=1
              scomp_bc=dir+1
              dcomp=1
              ncomp=nsolve
              ncomp_xvel=nsolve
              ncomp_cterm=nsolve
             else if (operation_flag.eq.OP_DIV_CELL) then ! divergence
              scomp=1
              scomp_bc=dir+1
              dcomp=1
              ncomp=1
              ncomp_xvel=nsolve
              ncomp_cterm=1
             ! MAC->CELL in solver or VELMAC_to_CELL
             else if (operation_flag.eq.OP_VEL_MAC_TO_CELL) then  ! velocity
              scomp=1
              scomp_bc=dir+1
              dcomp=dir+1
              ncomp=1
              ncomp_xvel=nsolve
              ncomp_cterm=1

             else if (operation_flag.eq.OP_VEL_DIVUP_TO_CELL) then ! div(up)
              print *,"expecting enable_spectral==0"
              stop

             else if (operation_flag.eq.OP_ISCHEME_CELL) then ! advection
              scomp=1
              scomp_bc=1
              dcomp=1
              ncomp=ncphys
              ncomp_xvel=nsolve
              ncomp_cterm=NFLUXSEM
             else
              print *,"operation_flag invalid9"
              stop
             endif
             
             if (operation_flag.eq.OP_VEL_DIVUP_TO_CELL) then !div(up) 
              print *,"expecting enable_spectral==0"
              stop
             else if ((operation_flag.eq.OP_RHS_CELL).or. & ! RHS for solver
                      (operation_flag.eq.OP_DIV_CELL).or. & ! divergence
                      (operation_flag.eq.OP_VEL_MAC_TO_CELL).or. & 
                      (operation_flag.eq.OP_ISCHEME_CELL)) then ! advection

              if (dir.eq.0) then 

               call SEM_MAC_TO_CELL( &
                ncomp_denold, &
                ncomp_veldest, &
                ncomp_dendest, &
                ns_time_order, &
                divu_outer_sweeps, &
                num_divu_outer_sweeps, &
                SDC_outer_sweeps, &
                level, &
                finest_level, &
                operation_flag, & 
                project_option, &
                energyflag, &
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
                vol_ptr, &
                xface_ptr, &
                xvel_ptr, &
                maskcoef_ptr, & ! 1=not covered, 0=covered
                cterm_ptr, &
                mdotcell_ptr, & !VELADVECT_MF, OP_ISCHEME_CELL
                maskdivres_ptr, & !DEN_RECON_MF, OP_ISCHEME_CELL
                pold_ptr, &
                denold_ptr, &
                ustar_ptr, &
                veldest_ptr, &
                dendest_ptr, &
                rhs_ptr)  ! divdest

              else if (dir.eq.1) then

               call SEM_MAC_TO_CELL( &
                ncomp_denold, &
                ncomp_veldest, &
                ncomp_dendest, &
                ns_time_order, &
                divu_outer_sweeps, &
                num_divu_outer_sweeps, &
                SDC_outer_sweeps, &
                level, &
                finest_level, &
                operation_flag, & 
                project_option, &
                energyflag, &
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
                vol_ptr, &
                yface_ptr, &
                yvel_ptr, &
                maskcoef_ptr, & ! 1=not covered, 0=covered
                cterm_ptr, &
                mdotcell_ptr, & !VELADVECT_MF, OP_ISCHEME_CELL
                maskdivres_ptr, & !DEN_RECON_MF, OP_ISCHEME_CELL
                pold_ptr, &
                denold_ptr, &
                ustar_ptr, &
                veldest_ptr, &
                dendest_ptr, &
                rhs_ptr)  ! divdest

              else if ((dir.eq.2).and.(SDIM.eq.3)) then

               call SEM_MAC_TO_CELL( &
                ncomp_denold, &
                ncomp_veldest, &
                ncomp_dendest, &
                ns_time_order, &
                divu_outer_sweeps, &
                num_divu_outer_sweeps, &
                SDC_outer_sweeps, &
                level, &
                finest_level, &
                operation_flag, & 
                project_option, &
                energyflag, &
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
                vol_ptr, &
                zface_ptr, &
                zvel_ptr, &
                maskcoef_ptr, & ! 1=not covered, 0=covered
                cterm_ptr, &
                mdotcell_ptr, & !VELADVECT_MF, OP_ISCHEME_CELL
                maskdivres_ptr, & !DEN_RECON_MF, OP_ISCHEME_CELL
                pold_ptr, &
                denold_ptr, &
                ustar_ptr, &
                veldest_ptr, &
                dendest_ptr, &
                rhs_ptr) ! divdest

              else
               print *,"dir invalid mac_to_cell2 "
               stop
              endif

             else
              print *,"operation_flag invalid"
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
                (operation_flag.ne.OP_ISCHEME_CELL)) then
        ! do nothing
       else
        print *,"bfact or operation_flag invalid"
        stop
       endif

      else if (enable_spectral.eq.0) then
       ! do nothing
      else
       print *,"enable_spectral invalid"
       stop
      endif

      return
      end subroutine fort_mac_to_cell


! OP_PRESGRAD_MAC (0)
! operation_flag=0  pressure gradient on MAC grid
! OP_PRES_CELL_TO_MAC (1)
! operation_flag=1  interpolate pressure from cell to MAC grid.
! OP_POTGRAD_TO_MAC (2)
! operation_flag=2  potential gradient on MAC grid, 
!                   surface tension on MAC grid
! OP_UNEW_CELL_TO_MAC (3)
! operation_flag=3  unew^MAC=unew^CELL->MAC
! OP_UNEW_USOL_MAC_TO_MAC (4)
! operation_flag=4  unew^MAC=uSOLID^MAC or uFLUID^MAC
! OP_UMAC_PLUS_VISC_CELL_TO_MAC (5)
! operation_flag=5  unew^MAC=unew^MAC+beta diffuse_ref^CELL->MAC
! OP_UGRAD_MAC (6)
! (operation_flag=6 reserved for rate of strain tensor)
! OP_ISCHEME_MAC (7)
! operation_flag=7  advection.
! OP_UGRAD_COUPLING_MAC (8)
! operation_flag=8  reserved for coupling terms in fort_crossterm
! OP_U_COMP_CELL_MAC_TO_MAC (11)
! operation_flag=11 
!   (i) unew^{f} in incompressible non-solid regions
!   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions 
!   (iii) (unew^{c})^{c->f}   compressible regions.
!   (iv) usolid in solid regions

      subroutine fort_cell_to_mac( &
       ncomp_mgoni, &
       ncomp_xp, & !local_MF[AMRSYNC_PRES_MF]->nComp() OP_PRESGRAD_MAC
       ncomp_xgp, &
       simple_AMR_BC_flag, &
       nsolve, &
       tileloop, &
       dir, &
       operation_flag, & 
       energyflag, & 
       beta, &
       visc_coef, &
       enable_spectral, &
       ncphys, &  ! nflux for advection
       constant_density_all_time, &
       presbc_in, &  ! denbc for advection
       velbc_in, &
       slab_step, &
       dt, &
       time, &
       xlo,dx, &
       spectral_loop, &
       ncfluxreg, & !local_sem_fluxreg_ncomp
       semflux,DIMS(semflux), &
       mask,DIMS(mask), & ! 1=fine/fine  0=coarse/fine
       maskcoef,DIMS(maskcoef), & ! 1=not cov. or outside domain  0=covered
       maskSEM,DIMS(maskSEM), &
       levelPC,DIMS(levelPC), &
       solfab,DIMS(solfab), &
       xcut,DIMS(xcut), &   ! coeff*areafrac
       xface,DIMS(xface), &  ! xflux for advection
       ! Umac_old if:
       !  OP_UMAC_PLUS_VISC_CELL_TO_MAC or
       !  OP_U_COMP_CELL_MAC_TO_MAC
       xgp,DIMS(xgp), & 
       xp,DIMS(xp), & ! holds AMRSYNC_PRES if OP_PRESGRAD_MAC
       xvel,DIMS(xvel), & !Umac_new
       vel,DIMS(vel), & !primary_velfab coming from increment_face_vel
        ! hydrostatic pressure if OP_POTGRAD_TO_MAC
       pres,DIMS(pres), & ! U_old(dir) if OP_U_COMP_CELL_MAC_TO_MAC
       den,DIMS(den), & !hydrostatic density if OP_POTGRAD_TO_MAC
        ! secondary_velfab if 
        !  OP_UNEW_CELL_TO_MAC, OP_UNEW_USOL_MAC_TO_MAC,
        !  OP_UMAC_PLUS_VISC_CELL_TO_MAC, or OP_U_COMP_CELL_MAC_TO_MAC
       mgoni,DIMS(mgoni), &!DIMS(dat)=datxlo,datxhi,datylo,datyhi,datzlo,datzhi
       colorfab,DIMS(colorfab), &
       typefab,DIMS(typefab), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact,bfact_c,bfact_f, &
       level,finest_level, &
       rz_flag, &
       domlo,domhi, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       blob_array, &
       blob_array_size, &
       num_colors, &
       project_option) &
      bind(c,name='fort_cell_to_mac')

      use global_utility_module
      use MOF_routines_module
      use probf90_module
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: dir
      INTEGER_T, INTENT(in) :: ncomp_mgoni
      INTEGER_T, INTENT(in) :: ncomp_xp
      INTEGER_T, INTENT(in) :: ncomp_xgp
      INTEGER_T, INTENT(in) :: simple_AMR_BC_flag
      INTEGER_T, INTENT(in) :: nsolve
      INTEGER_T, INTENT(in) :: tileloop
      INTEGER_T, INTENT(in) :: spectral_loop
      INTEGER_T, INTENT(in) :: ncfluxreg
      INTEGER_T, INTENT(in) :: nparts
      INTEGER_T, INTENT(in) :: nparts_def
      INTEGER_T, INTENT(in) :: im_solid_map(nparts_def)

      INTEGER_T, INTENT(in) :: blob_array_size
      INTEGER_T, INTENT(in) :: num_colors
      REAL_T, INTENT(in) :: blob_array(blob_array_size)
        
      INTEGER_T, INTENT(in) :: slab_step
      INTEGER_T, INTENT(in) :: operation_flag
      INTEGER_T, INTENT(in) :: energyflag
      INTEGER_T, INTENT(in) :: enable_spectral
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: ncphys  ! nflux for advection
      INTEGER_T, INTENT(in) :: constant_density_all_time(num_materials)
      REAL_T, INTENT(in) :: dt
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(in) :: beta,visc_coef
      REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, INTENT(in) :: DIMDEC(semflux)
      INTEGER_T, INTENT(in) :: DIMDEC(mask)
      INTEGER_T, INTENT(in) :: DIMDEC(maskcoef)
      INTEGER_T, INTENT(in) :: DIMDEC(maskSEM)
      INTEGER_T, INTENT(in) :: DIMDEC(xcut)
      INTEGER_T, INTENT(in) :: DIMDEC(xface)
      INTEGER_T, INTENT(in) :: DIMDEC(xgp)
      INTEGER_T, INTENT(in) :: DIMDEC(xp)
      INTEGER_T, INTENT(in) :: DIMDEC(xvel)
      INTEGER_T, INTENT(in) :: DIMDEC(pres)
      INTEGER_T, INTENT(in) :: DIMDEC(vel)
      INTEGER_T, INTENT(in) :: DIMDEC(den)
      INTEGER_T, INTENT(in) :: DIMDEC(mgoni)
      INTEGER_T, INTENT(in) :: DIMDEC(levelPC)
      INTEGER_T, INTENT(in) :: DIMDEC(solfab)
      INTEGER_T, INTENT(in) :: DIMDEC(colorfab)
      INTEGER_T, INTENT(in) :: DIMDEC(typefab)

       ! denbc for advect
      INTEGER_T, INTENT(in) :: &
              presbc_in(SDIM,2,num_materials*num_state_material) 
      INTEGER_T, INTENT(in) :: velbc_in(SDIM,2,SDIM)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact,bfact_c,bfact_f
      INTEGER_T, INTENT(in) :: rz_flag
      INTEGER_T, INTENT(in) :: domlo(SDIM),domhi(SDIM)
      INTEGER_T, INTENT(in) :: project_option

      REAL_T, INTENT(in), target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))

      REAL_T, INTENT(in), target :: maskcoef(DIMV(maskcoef))
      REAL_T, pointer :: maskcoef_ptr(D_DECL(:,:,:))

      REAL_T, INTENT(in), target :: maskSEM(DIMV(maskSEM))
      REAL_T, pointer :: maskSEM_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: &
              levelPC(DIMV(levelPC),num_materials*(1+SDIM))
      REAL_T, pointer :: levelPC_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: solfab(DIMV(solfab),nparts_def*SDIM)
      REAL_T, pointer :: solfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: semflux(DIMV(semflux),ncfluxreg)
      REAL_T, pointer :: semflux_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: xcut(DIMV(xcut),1)
      REAL_T, pointer :: xcut_ptr(D_DECL(:,:,:),:)
       ! xflux for advection
      REAL_T, INTENT(inout), target :: xface(DIMV(xface),ncphys) 
      REAL_T, pointer :: xface_ptr(D_DECL(:,:,:),:)
       !holds Umac_old if 
       ! OP_UMAC_PLUS_VISC_CELL_TO_MAC or OP_U_COMP_CELL_MAC_TO_MAC
      REAL_T, INTENT(inout), target :: xgp(DIMV(xgp),ncomp_xgp) 
      REAL_T, pointer :: xgp_ptr(D_DECL(:,:,:),:)

        ! when xp corresponds to PEDGE_MF:
        !  xp(VALID_PEDGE+1)
        !  xp(PRESSURE_PEDGE+1)
      REAL_T, INTENT(inout), target :: xp(DIMV(xp),ncomp_xp)
      REAL_T, pointer :: xp_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: xvel(DIMV(xvel),1)
      REAL_T, pointer :: xvel_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: vel(DIMV(vel),SDIM)
      REAL_T, pointer :: vel_ptr(D_DECL(:,:,:),:)
       ! holds U_old(dir) if OP_U_COMP_CELL_MAC_TO_MAC
      REAL_T, INTENT(in), target :: pres(DIMV(pres),1)
      REAL_T, pointer :: pres_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target ::  &
              den(DIMV(den),num_materials*num_state_material)
      REAL_T, pointer :: den_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: mgoni(DIMV(mgoni),ncomp_mgoni)
      REAL_T, pointer :: mgoni_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: typefab(DIMV(typefab))
      REAL_T, pointer :: typefab_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: colorfab(DIMV(colorfab))
      REAL_T, pointer :: colorfab_ptr(D_DECL(:,:,:))
  
      INTEGER_T i,j,k,ii,jj,kk
      REAL_T pplus
      REAL_T pminus
      REAL_T pgrad
      REAL_T pgrad_gravity ! grad ppot/den_pot  ppot=dt * rho g z
      REAL_T incremental_gravity
      REAL_T pgrad_tension
      REAL_T gradh_tension
      REAL_T gradh_gravity
      REAL_T dplus,dminus
      REAL_T den_H,pres_H
      REAL_T den_im,den_im_opp
      REAL_T interp_factor
      INTEGER_T dencomp_im,dencomp_im_opp

       !OP_PRES_CELL_TO_MAC 
       !(1)use_face_pres=VALID_PEDGE+1
       !(2)face pressure=PRESSURE_PEDGE+1
      REAL_T PEDGE_local(NCOMP_PEDGE)
      INTEGER_T im1,jm1,km1
      INTEGER_T im,im_opp
      INTEGER_T im_gravity,im_opp_gravity
      INTEGER_T im_heat,tcomp
      INTEGER_T iten
      INTEGER_T im_left,im_right
      INTEGER_T im_left_tension,im_right_tension
      INTEGER_T im_left_gravity,im_right_gravity
      INTEGER_T dir2,side
      INTEGER_T velcomp,iboundary
      REAL_T cutedge,RR
      INTEGER_T, parameter :: nhalf=3
      REAL_T xstenMAC(-nhalf:nhalf,SDIM)
      REAL_T xstenMAC_center(SDIM)
      REAL_T DXMAXLS
      REAL_T local_vel_MAC
      REAL_T local_vel_old_MAC
      REAL_T uedge
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
      REAL_T velsum_primary
      REAL_T mass_sum
      REAL_T DMface
      REAL_T primary_velmaterial
      REAL_T velmaterialMAC
      REAL_T solid_velocity
      INTEGER_T mask_covered(2)
      INTEGER_T mask_coarsefine(2)
      INTEGER_T at_reflect_wall,at_wall,at_ext_wall
      INTEGER_T use_face_pres
      INTEGER_T at_coarse_fine_wallF
      INTEGER_T at_coarse_fine_wallC
      INTEGER_T static_flag
      REAL_T user_tension(num_interfaces)
      REAL_T tension_scaled
      REAL_T LSleft(num_materials)
      REAL_T LSright(num_materials)
      REAL_T LSleft_grav
      REAL_T LSright_grav
      REAL_T LSupwind(num_materials)
      REAL_T localLS(num_materials)
      REAL_T mgoni_temp(num_materials)
      INTEGER_T local_maskSEM
      INTEGER_T maskcov
      REAL_T hx
      INTEGER_T scomp,scomp_bc,dcomp
      INTEGER_T ncomp_dest,ncomp_source
      INTEGER_T nc ! in: cell_to_mac
      INTEGER_T ibase,idonate,jdonate,kdonate
      INTEGER_T stripstat
      INTEGER_T elemlo(3),elemhi(3)
      INTEGER_T ielem,jelem,kelem
      REAL_T AFACE
      REAL_T AFACE_ICE
      REAL_T denlocal,templocal
      REAL_T test_velocity_FACE
      INTEGER_T im_solid
      INTEGER_T im_prescribed
      INTEGER_T im_solid_valid
      INTEGER_T im_prescribed_valid
      INTEGER_T partid_solid
      INTEGER_T partid_prescribed
      INTEGER_T partid_check
      REAL_T cutoff
      INTEGER_T typeleft,typeright,typeface
      INTEGER_T colorleft,colorright,colorface
      INTEGER_T face_velocity_override
      INTEGER_T FSI_prescribed_flag

      REAL_T xclamped_minus_sten(-nhalf:nhalf,SDIM)
      REAL_T xclamped_plus_sten(-nhalf:nhalf,SDIM)
      REAL_T xclamped_minus(SDIM)
      REAL_T xclamped_plus(SDIM)
      REAL_T LS_clamped_plus
      REAL_T LS_clamped_minus
      REAL_T vel_clamped_plus(SDIM)
      REAL_T vel_clamped_minus(SDIM)
      REAL_T temperature_clamped_plus
      REAL_T temperature_clamped_minus
      INTEGER_T prescribed_flag
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped
      INTEGER_T is_clamped_face
      INTEGER_T local_compressible

      REAL_T test_current_icefacecut
      REAL_T test_current_icemask

      INTEGER_T :: homogeneous_rigid_velocity

      INTEGER_T, parameter :: DEBUG_PRESCRIBED=0
      REAL_T DEBUG_PRESCRIBED_VEL_TOT
      REAL_T DEBUG_PRESCRIBED_VEL_DEN

      DEBUG_PRESCRIBED_VEL_TOT=zero
      DEBUG_PRESCRIBED_VEL_DEN=zero

      homogeneous_rigid_velocity=0

      semflux_ptr=>semflux
      xcut_ptr=>xcut
      xface_ptr=>xface
      xgp_ptr=>xgp
      xp_ptr=>xp
      xvel_ptr=>xvel
      mask_ptr=>mask
      maskcoef_ptr=>maskcoef
      maskSEM_ptr=>maskSEM
      vel_ptr=>vel
      pres_ptr=>pres
      den_ptr=>den
      solfab_ptr=>solfab
      mgoni_ptr=>mgoni
      typefab_ptr=>typefab
      colorfab_ptr=>colorfab
      levelPC_ptr=>levelPC

      if (ncomp_xp.lt.1) then
       print *,"ncomp_xp invalid(1) ",ncomp_xp
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_cell_to_mac"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_cell_to_mac"
       stop
      endif
      if ((simple_AMR_BC_flag.eq.0).or.(simple_AMR_BC_flag.eq.1)) then
       ! do nothing
      else
       print *,"simple_AMR_BC_flag invalid"
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
 
      if ((enable_spectral.lt.0).or. &
          (enable_spectral.gt.1)) then
       print *,"enable_spectral invalid cell to mac"
       stop
      endif
      if ((spectral_loop.ne.0).and. &
          (spectral_loop.ne.1)) then
       print *,"spectral_loop invalid"
       stop
      endif

      if (project_option_is_validF(project_option).eq.1) then

       if ((ncfluxreg/SDIM)*SDIM.ne.ncfluxreg) then
        print *,"ncfluxreg invalid11 ",ncfluxreg
        stop
       endif
       if (ncfluxreg.lt.SDIM) then
        print *,"ncfluxreg invalid12 ",ncfluxreg
        stop
       endif

      else
       print *,"project_option invalid"
       stop
      endif

      if (blob_array_size.eq.1) then
       if (num_colors.eq.0) then
        ! do nothing
       else
        print *,"num_colors inconsistent; fort_cell_to_mac"
        stop
       endif
      else if (blob_array_size.eq.num_colors*num_elements_blobclass) then
       if (num_colors.ge.1) then
        ! do nothing
       else
        print *,"num_colors inconsistent; fort_cell_to_mac"
        stop
       endif
      else
       print *,"blob_array_size invalid fort_cell_to_mac"
       stop
      endif

      if (num_colors.ge.0) then
       ! do nothing
      else
       print *,"num_colors invalid fort_cell_to_mac"
       stop
      endif

      call fort_check_operation_flag_MAC(operation_flag)

      if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

       if (ncomp_mgoni.eq.num_materials*num_state_material) then
        ! do nothing
       else
        print *,"ncomp_mgoni invalid"
        stop
       endif

       if (ncomp_xp.ne.NFLUXSEM) then
        print *,"ncomp_xp.ne.NFLUXSEM (OP_ISCHEME_MAC) ",ncomp_xp
        stop
       endif
       if (ncomp_xgp.ne.NFLUXSEM) then
        print *,"ncomp_xp invalid(3) ",ncomp_xp
        stop
       endif
       
       if (ncphys.ne.NFLUXSEM) then
        print *,"ncphys invalid"
        stop
       endif

      else if (operation_flag.eq.OP_PRESGRAD_MAC) then

       if (ncomp_mgoni.eq.1) then
        ! do nothing
       else
        print *,"ncomp_mgoni invalid"
        stop
       endif

       if (project_option_is_validF(project_option).eq.1) then
        if (ncphys.ne.FACECOMP_NCOMP) then
         print *,"ncphys invalid"
         stop
        endif
       else
        print *,"project_option invalid"
        stop
       endif

       if (ncomp_xp.ne.nsolve) then
        print *,"ncomp_xp.ne.nsolve (OP_PRESGRAD_MAC) ",ncomp_xp
        stop
       endif

       if (nsolve.eq.1) then
        ! do nothing
       else 
        print *,"expecting nsolve==1 if OP_PRESGRAD_MAC"
        stop
       endif

      else if (operation_flag.eq.OP_PRES_CELL_TO_MAC) then ! P^{cell->mac}
       
       if (ncphys.ne.FACECOMP_NCOMP) then
        print *,"ncphys invalid"
        stop
       endif
       if (ncomp_mgoni.eq.1) then
        ! do nothing
       else
        print *,"ncomp_mgoni invalid"
        stop
       endif
       if (ncomp_xp.ne.NCOMP_PEDGE) then
        print *,"ncomp_xp.ne.NCOMP_PEDGE(OP_PRES_CELL_TO_MAC) ",ncomp_xp
        stop
       endif
       if (enable_spectral.eq.0) then
        ! do nothing
       else
        print *,"must have enable_spectral==0 if OP_PRES_CELL_TO_MAC"
        stop
       endif

       !potential gradient
      else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 

       if (ncphys.ne.FACECOMP_NCOMP) then
        print *,"ncphys invalid"
        stop
       endif
       if (ncomp_mgoni.eq.num_materials*num_state_material) then
        ! do nothing
       else
        print *,"ncomp_mgoni invalid"
        stop
       endif
       if (ncomp_xp.ne.1) then
        print *,"ncomp_xp.ne.1(OP_POTGRAD_TO_MAC) ",ncomp_xp
        stop
       endif

      else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. &
               (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC).or. &
               (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC)) then

       if (ncomp_mgoni.eq.SDIM) then
        ! do nothing
       else
        print *,"ncomp_mgoni invalid(OP_UNEW or OP_UMAC)"
        stop
       endif
       if (ncomp_xp.ne.NCOMP_AMRSYNC_VEL_MF) then
        print *,"ncomp_xp.ne.NCOMP_AMRSYNC_VEL_MF(OP_UNEW or OP_UMAC) ", &
           ncomp_xp
        stop
       endif

       if (ncphys.ne.FACECOMP_NCOMP) then
        print *,"ncphys invalid(OP_UNEW or OP_UMAC)"
        stop
       endif

      else if (operation_flag.eq.OP_UGRAD_MAC) then

       print *,"fort_face_gradients calls sem_cell_to_mac with op=6"
       print *,"(OP_UGRAD_MAC)"
       stop

      else if (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC) then

       if (ncphys.ne.FACECOMP_NCOMP) then
        print *,"ncphys invalid"
        stop
       endif

      else 
       print *,"operation_flag invalid11:",operation_flag
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (operation_flag.eq.OP_PRES_CELL_TO_MAC) then ! p^CELL->MAC
       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_PRES_CELL_TO_MAC"
        stop
       endif
       if ((project_option.eq.SOLVETYPE_PRESGRAVITY).or. &
           (project_option.eq.SOLVETYPE_PRESSTATIC)) then
        homogeneous_rigid_velocity=1
       endif
      else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 
       if ((energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
            POTGRAD_SURFTEN_INCREMENTAL_GRAV).or. &
           (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
            POTGRAD_SURFTEN_BASE_GRAV).or. &
           (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
            POTGRAD_INCREMENTAL_GRAV).or. &
           (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
            POTGRAD_BASE_GRAV).or. &
           (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
            POTGRAD_SURFTEN)) then
        ! do nothing
       else
        print *,"energyflag invalid OP_POTGRAD_TO_MAC"
        stop
       endif
       if ((ncomp_xp.ne.1).or. &
           (ncomp_xgp.ne.1)) then
        print *,"ncomp_xp or ncomp_xgp invalid"
        stop
       endif
      else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. &
               (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC).or. & 
               (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
               (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC)) then 
       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_U etc CELL/MAC to MAC"
        stop
       endif
       if (project_option_momeqnF(project_option).eq.1) then
        ! do nothing
       else
        print *,"project_option_momeqnF(project_option) invalid"
        stop
       endif

       if ((project_option.eq.SOLVETYPE_PRESGRAVITY).or. &
           (project_option.eq.SOLVETYPE_PRESSTATIC)) then
        if (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC) then
         ! do nothing
        else
         print *,"OP_UNEW_USOL_MAC_TO_MAC expected"
         stop
        endif
        homogeneous_rigid_velocity=1
       endif

      else if (operation_flag.eq.OP_PRESGRAD_MAC) then ! (grad p)_MAC
       if ((energyflag.ne.SUB_OP_FOR_MAIN).and. &
           (energyflag.ne.SUB_OP_FOR_SDC)) then
        print *,"energyflag invalid OP_PRESGRAD_MAC"
        stop
       endif
      else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

       if (energyflag.ne.SUB_OP_DEFAULT) then
        print *,"energyflag invalid OP_ISCHEME_MAC"
        stop
       endif

      else if (operation_flag.eq.OP_UGRAD_MAC) then

       print *,"fort_face_gradients calls sem_cell_to_mac with op=6"
       print *,"(OP_UGRAD_MAC)"
       stop

      else
       print *,"operation_flag invalid12"
       stop
      endif

      if (project_option_is_validF(project_option).eq.1) then
       ! do nothing
      else
       print *,"project_option invalid"
       stop
      endif

      do im=1,num_materials

       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"denconst invalid"
        stop
       endif

      enddo ! im=1..num_materials

      if ((slab_step.lt.-1).or.(slab_step.gt.bfact_time_order)) then
       print *,"slab_step invalid cell to mac"
       stop
      endif

       ! FACE_WEIGHT_MF
      call checkbound_array(fablo,fabhi,xcut_ptr,0,dir)
      call checkbound_array(fablo,fabhi,xgp_ptr,0,dir)
      call checkbound_array(fablo,fabhi,xvel_ptr,0,dir)
      call checkbound_array(fablo,fabhi,vel_ptr,1,-1)
      call checkbound_array(fablo,fabhi,pres_ptr,1,-1)
      call checkbound_array(fablo,fabhi,den_ptr,1,-1)

      call checkbound_array(fablo,fabhi,mgoni_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,typefab_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,colorfab_ptr,1,-1)
      call checkbound_array(fablo,fabhi,levelPC_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,maskcoef_ptr,1,-1)

      if (project_option_is_validF(project_option).eq.1) then
         !local_face_var_mf
       call checkbound_array(fablo,fabhi,xface_ptr,0,dir)
         !local_amrsync_pres_mf
       call checkbound_array(fablo,fabhi,xp_ptr,0,dir)
         !local_sem_fluxreg_mf
       call checkbound_array(fablo,fabhi,semflux_ptr,1,-1)

       if (project_option_is_static(project_option).eq.1) then
        ! do nothing
       else if (project_option_is_static(project_option).eq.0) then
        !local_fsi_ghost_mac_mf+dir
        call checkbound_array(fablo,fabhi,solfab_ptr,0,dir)
       else
        print *,"project_option_is_static invalid"
        stop
       endif
       call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1)
      else
       print *,"project_option invalid"
       stop
      endif

      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=DXMAXLS

      do im=1,num_materials

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).ge.1).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
        ! do nothing
       else
        print *,"fort_material_type invalid"
        stop
       endif

      enddo  ! im=1..num_materials

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
       print *,"dir out of range in fort_cell_to_mac, dir=",dir
       stop
      endif 

       ! first low order gradients.
      if (tileloop.eq.0) then

       if (spectral_loop.eq.0) then

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)
        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

          ! dir=0..sdim-1
         call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir)
         do dir2=1,SDIM
          xstenMAC_center(dir2)=xstenMAC(0,dir2)
         enddo

         is_clamped_face=-1

         local_compressible=0
         im_left=0
         im_right=0

         if (levelrz.eq.COORDSYS_CARTESIAN) then
          RR=one
         else if (levelrz.eq.COORDSYS_RZ) then
          RR=one
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
          if (dir.eq.1) then
           RR=xstenMAC_center(1)
          else
           RR=one
          endif
         else
          print *,"levelrz invalid edgegradp"
          stop
         endif 

         hx=xstenMAC(1,dir+1)-xstenMAC(-1,dir+1)

         hx=hx*RR

         if (hx.gt.zero) then
          ! do nothing
         else
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
          print *,"dir invalid fort_cell_to_mac, dir=",dir
          stop
         endif

         if (operation_flag.eq.OP_PRES_CELL_TO_MAC) then ! p^CELL->MAC

          do im=1,ncphys
           local_face(im)=xface(D_DECL(i,j,k),im)
          enddo
           ! newly projected face velocity might be overwritten
           ! with the solid velocity.
          local_vel_MAC=xvel(D_DECL(i,j,k),1)

         else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 

          do im=1,ncphys
           local_face(im)=xface(D_DECL(i,j,k),im)
          enddo

         else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. &
                  (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC).or. & 
                  (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                  (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC)) then 

          do im=1,ncphys
           local_face(im)=xface(D_DECL(i,j,k),im)
          enddo
           !Umac_new
          local_vel_MAC=xvel(D_DECL(i,j,k),1)
          if ((operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. &
              (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC)) then
           local_vel_old_MAC=xgp(D_DECL(i,j,k),1) !Umac_old
          else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. &
                   (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC)) then
           local_vel_old_MAC=zero
          else
           print *,"operation_flag invalid13"
           stop
          endif

         else if (operation_flag.eq.OP_PRESGRAD_MAC) then ! (grad p)_MAC

          ! do nothing

         else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

          local_vel_MAC=xvel(D_DECL(i,j,k),1)

         else if (operation_flag.eq.OP_UGRAD_MAC) then

          print *,"fort_face_gradients calls sem_cell_to_mac with op=6"
          print *,"(OP_UGRAD_MAC)"
          stop
 
         else
          print *,"operation_flag invalid14"
          stop
         endif

          ! set LSleft, LSright, localLS
         if ((operation_flag.eq.OP_POTGRAD_TO_MAC).or. & 
             (operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. &
             (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC).or. &
             (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. &
             (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC).or. & 
             (operation_flag.eq.OP_PRES_CELL_TO_MAC)) then  ! p^CELL->MAC

          ! levelPC() has piecewise constant BC at coarse/fine borders.
          do im=1,num_materials
           LSleft(im)=levelPC(D_DECL(im1,jm1,km1),im)
           LSright(im)=levelPC(D_DECL(i,j,k),im)
           localLS(im)=half*(LSright(im)+LSleft(im))
          enddo
          call get_primary_material(LSleft,im_left)
          call get_primary_material(LSright,im_right)

          if ((is_compressible_mat(im_left).eq.1).and. &
              (is_compressible_mat(im_right).eq.1)) then
           local_compressible=1
          else if ((is_compressible_mat(im_left).eq.0).or. &
                   (is_compressible_mat(im_right).eq.0)) then
           local_compressible=0
          else
           print *,"is_compressible invalid"
           stop
          endif

          call gridsten_level(xclamped_minus_sten,im1,jm1,km1,level,nhalf)
          call gridsten_level(xclamped_plus_sten,i,j,k,level,nhalf)
          do dir2=1,SDIM
           xclamped_minus(dir2)=xclamped_minus_sten(0,dir2)
           xclamped_plus(dir2)=xclamped_plus_sten(0,dir2)
          enddo
          call SUB_clamped_LS(xclamped_minus,time,LS_clamped_minus, &
           vel_clamped_minus,temperature_clamped_minus,prescribed_flag,dx)
          call SUB_clamped_LS(xclamped_plus,time,LS_clamped_plus, &
           vel_clamped_plus,temperature_clamped_plus,prescribed_flag,dx)
          if ((LS_clamped_plus.ge.zero).and. &
              (LS_clamped_minus.ge.zero)) then
           is_clamped_face=1
           do dir2=1,SDIM
            vel_clamped(dir2)=half*(vel_clamped_minus(dir2)+ &
             vel_clamped_plus(dir2))
           enddo
           temperature_clamped=half*(temperature_clamped_minus+ &
             temperature_clamped_plus)
          else if ((LS_clamped_plus.lt.zero).and. &
                   (LS_clamped_minus.lt.zero)) then
           is_clamped_face=0
          else if ((LS_clamped_plus.ge.zero).and. &
                   (LS_clamped_minus.lt.zero)) then
           is_clamped_face=2
           do dir2=1,SDIM
            vel_clamped(dir2)=vel_clamped_plus(dir2)
           enddo
           temperature_clamped=temperature_clamped_plus
          else if ((LS_clamped_plus.lt.zero).and. &
                   (LS_clamped_minus.ge.zero)) then
           is_clamped_face=3
           do dir2=1,SDIM
            vel_clamped(dir2)=vel_clamped_minus(dir2)
           enddo
           temperature_clamped=temperature_clamped_minus
          else
           print *,"LS_clamped_plus or LS_clamped_minus NaN"
           stop
          endif

          if (homogeneous_rigid_velocity.eq.0) then
           ! do nothing
          else if (homogeneous_rigid_velocity.eq.1) then
           do dir2=1,SDIM
            vel_clamped(dir2)=zero
           enddo
          else
           print *,"homogeneous_rigid_velocity invalid"
           stop
          endif

         else if ((operation_flag.eq.OP_PRESGRAD_MAC).or. &
                  (operation_flag.eq.OP_ISCHEME_MAC)) then ! advection

          ! do nothing

         else if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then

          print *,"fort_crosterm calls lineGRAD"
          print *,"(OP_UGRAD_COUPLING_MAC)"
          stop

         else if (operation_flag.eq.OP_UGRAD_MAC) then

          print *,"fort_face_gradients calls sem_cell_to_mac with op=6"
          print *,"(OP_UGRAD_MAC)"
          stop

         else 
          print *,"operation_flag invalid15:",operation_flag
          stop
         endif

          ! advection
          ! The low order fluxes here might be needed by the high order
          ! divergence operator if there is a spectral element next
          ! to a low order element.
         if (operation_flag.eq.OP_ISCHEME_MAC) then 

          do nc=1,ncphys
           local_face(nc)=zero
          enddo

          at_RZ_face=0
          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if (levelrz.eq.COORDSYS_RZ) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if ((dir.eq.0).and. &
               (xstenMAC_center(1).le.VOFTOL*dx(1))) then
            at_RZ_face=1
           endif
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           if ((dir.eq.0).and. &
               (xstenMAC_center(1).le.VOFTOL*dx(1))) then
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
          else if (at_RZ_face.eq.0) then
           test_velocity_FACE=local_vel_MAC
           if (test_velocity_FACE.ge.zero) then
            idonate=im1
            jdonate=jm1
            kdonate=km1
           else if (test_velocity_FACE.lt.zero) then
            idonate=i
            jdonate=j
            kdonate=k
           else
            print *,"test_velocity_FACE corrupt"
            stop
           endif
          else
           print *,"at_RZ_face invalid"
           stop
          endif

          do im=1,num_materials
           LSupwind(im)=levelPC(D_DECL(idonate,jdonate,kdonate),im)
          enddo

          call get_primary_material(LSupwind,im)
          if ((im.ge.1).and.(im.le.num_materials)) then
           ibase=num_state_material*(im-1) 
           denlocal=den(D_DECL(idonate,jdonate,kdonate),ibase+ENUM_DENVAR+1) 
           if (denlocal.gt.zero) then
            ! do nothing
           else
            print *,"denlocal invalid"
            stop
           endif
           if (constant_density_all_time(im).eq.1) then
            if (abs(denlocal-fort_denconst(im)).le. &
                VOFTOL*fort_denconst(im)) then
             ! do nothing
            else
             print *,"expecting denlocal=fort_denconst(im)"
             stop
            endif
           else if (constant_density_all_time(im).eq.0) then
            ! do nothing
           else
            print *,"constant_density_all_time invalid"
            stop
           endif
           templocal= &
             den(D_DECL(idonate,jdonate,kdonate),ibase+ENUM_TEMPERATUREVAR+1) 

           do nc=1,ncphys
            if ((nc.ge.SEM_U+1).and.(nc.le.SEM_W+1)) then
             velcomp=nc
             if (at_RZ_face.eq.1) then
              local_face(nc)=zero
             else if (at_RZ_face.eq.0) then
              local_face(nc)=vel(D_DECL(idonate,jdonate,kdonate),velcomp)
              local_face(nc)=local_face(nc)*test_velocity_FACE
             else
              print *,"at_RZ_face invalid"
              stop
             endif
            else if (nc.eq.SEM_T+1) then
             local_face(nc)=templocal 
             local_face(nc)=local_face(nc)*test_velocity_FACE
            else
             print *,"nc invalid"
             stop
            endif
           enddo !nc=1,ncphys

          else
           print *,"im invalid 15446:",im
           stop
          endif

         else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
                  (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC).or. & 
                  (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                  (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC)) then 

          face_velocity_override=0

          at_RZ_face=0
          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if (levelrz.eq.COORDSYS_RZ) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if ((dir.eq.0).and. &
               (xstenMAC_center(1).le.VOFTOL*dx(1))) then
            at_RZ_face=1
           endif
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           if ((dir.eq.0).and. &
               (xstenMAC_center(1).le.VOFTOL*dx(1))) then
            at_RZ_face=1
           endif
          else
           print *,"levelrz invalid tfrmac"
           stop
          endif 

          uedge=zero

          if (at_RZ_face.eq.1) then
           face_velocity_override=1
           uedge=zero
          else if (at_RZ_face.eq.0) then

           fluid_volface=zero
           not_prescribed_volface=zero
           volface=zero

           do side=1,2
            do im=1,num_materials
             local_volume=local_face(FACECOMP_VOFFACE+2*(im-1)+side)
             volface=volface+local_volume
             if (is_prescribed(im).eq.0) then
              not_prescribed_volface=not_prescribed_volface+local_volume
             else if (is_prescribed(im).eq.1) then
              ! do nothing
             else
              print *,"is_prescribed(im) invalid"
              stop
             endif
             if (is_rigid(im).eq.0) then
              if (is_prescribed(im).eq.0) then
               fluid_volface=fluid_volface+local_volume
              else
               print *,"is_prescribed(im) invalid"
               stop
              endif
             else if (is_rigid(im).eq.1) then
              ! do nothing
             else
              print *,"is_rigid invalid LEVELSET_3D.F90"
              stop
             endif
            enddo ! im=1..num_materials
           enddo ! side=1,2

           if (volface.gt.zero) then
            ! do nothing
           else
            print *,"not_prescribed_volface ",not_prescribed_volface
            print *,"fluid_volface ",fluid_volface
            print *,"volface ",volface
            print *,"volface bust tfrmac"
            stop
           endif

           fluid_volface=fluid_volface/volface
           not_prescribed_volface=not_prescribed_volface/volface

           if ((local_face(FACECOMP_FACECUT+1).ge.zero).and. &
               (local_face(FACECOMP_FACECUT+1).le.half)) then
            fluid_volface=zero
            not_prescribed_volface=zero
           else if ((local_face(FACECOMP_FACECUT+1).ge.half).and. &
                    (local_face(FACECOMP_FACECUT+1).le.one)) then
            fluid_volface= &
             min(fluid_volface,local_face(FACECOMP_FACECUT+1))
            not_prescribed_volface= &
             min(not_prescribed_volface,local_face(FACECOMP_FACECUT+1))
           else
            print *,"local_face(FACECOMP_FACECUT+1) invalid"
            stop
           endif
         
           ! is_solid_face==1 if:
           !   0.0<=fluid_volface<=VOFTOL_AREAFRAC  or
           !   max(LSleft(im_solid),LSright(im_solid))>=0.0
           ! is_prescribed_face==1 if:
           !   0.0<=not_prescribed_prescribed<=VOFTOL_AREAFRAC  or
           !   max(LSleft(im_prescribed),LSright(im_prescribed))>=0.0
           ! fixed_face is declared in: PROB.F90
           call fixed_face( &
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
              face_velocity_override=1
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

              if (homogeneous_rigid_velocity.eq.0) then
               uedge=solfab(D_DECL(i,j,k),velcomp)
              else if (homogeneous_rigid_velocity.eq.1) then
               uedge=zero
              else
               print *,"homogeneous_rigid_velocity invalid"
               stop
              endif

              DEBUG_PRESCRIBED_VEL_TOT=DEBUG_PRESCRIBED_VEL_TOT+uedge
              DEBUG_PRESCRIBED_VEL_DEN=DEBUG_PRESCRIBED_VEL_DEN+one
             else if (im_prescribed_valid.eq.0) then
              uedge=zero
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

            if (project_option.eq.SOLVETYPE_PRESGRAVITY) then
             test_current_icefacecut=one
            else
             test_current_icefacecut= &
                 xface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
            endif
          
            if ((test_current_icefacecut.ge.zero).and. &
                (test_current_icefacecut.le.one)) then

              ! test_current_icemask=zero for both ice materials and
              ! "is_FSI_rigid" materials.
             if (project_option.eq.SOLVETYPE_PRESGRAVITY) then
              test_current_icemask=one
             else
              test_current_icemask=xface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)
             endif

             if ((test_current_icemask.eq.zero).or. &
                 (test_current_icemask.eq.one)) then

              velsum_primary=zero
              mass_sum=zero

              do side=1,2

               partid_check=0
               do im=1,num_materials

                DMface=local_face(FACECOMP_MASSFACE+2*(im-1)+side)
                if (DMface.gt.zero) then
                 ! do nothing
                else if (DMface.eq.zero) then
                 ! do nothing
                else
                 print *,"DMface invalid"
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

                 !primary_vel_data="vel"=CURRENT_CELL_VEL_MF; 
                 !secondary_vel_data="mgoni"=CURRENT_CELL_VEL_MF; 
                if (operation_flag.eq.OP_UNEW_CELL_TO_MAC) then ! cell -> MAC
                 velcomp=dir+1
                 primary_velmaterial=vel(D_DECL(ic,jc,kc),velcomp)

                 !primary_vel_data="vel"=CURRENT_CELL_VEL_MF; 
                 !secondary_vel_data="mgoni"=CURRENT_CELL_VEL_MF; 
                else if (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC) then 
                 velcomp=1
                 primary_velmaterial=local_vel_MAC

                 !primary_vel_data="vel"=idx_velcell;  // increment
                 !secondary_vel_data="mgoni"=CURRENT_CELL_VEL_MF; 
                else if (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC) then
                 velcomp=dir+1
                  ! local_vel_MAC=xvel
                  ! local_vel_old_MAC=xgp (a copy of xvel)
                 velmaterialMAC=local_vel_old_MAC

                  !secondary_vel_data="mgoni"=CURRENT_CELL_VEL_MF; 

                 primary_velmaterial= &
                   velmaterialMAC+beta*vel(D_DECL(ic,jc,kc),velcomp)
                
                 if ((beta.eq.-one).or.(beta.eq.one)) then
                   ! do nothing
                 else
                  print *,"beta invalid"
                  stop
                 endif

                 !called after advection to update u^{advect,MAC}
                 !primary_vel_data="vel"=DELTA_CELL_VEL_MF; 
                 !secondary_vel_data="mgoni"=CURRENT_CELL_VEL_MF; 
                else if (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC) then 
                 if (local_compressible.eq.0) then
                  velcomp=1
                   !local_vel_MAC=xvel=Umac_new=UMAC^{ADVECT}
                  primary_velmaterial=local_vel_MAC
                 else if (local_compressible.eq.1) then
                  ! UMAC^{ADVECT}= 
                  !   I_{CELL}^{MAC} (U_CELL^{ADVECT})
                  velcomp=dir+1
                   !"mgoni"=secondary_vel_data=CURRENT_CELL_VEL_MF; 
                  primary_velmaterial=mgoni(D_DECL(ic,jc,kc),velcomp) 
                 else
                  print *,"local_compressible invalid"
                  stop
                 endif
                else
                 print *,"operation_flag invalid16"
                 stop
                endif

                if (is_lag_part(im).eq.1) then
                 partid_check=partid_check+1
                else if (is_lag_part(im).eq.0) then
                 ! do nothing
                else
                 print *,"is_lag_part(im) invalid"
                 stop
                endif

                mass_sum=mass_sum+DMface
                velsum_primary=velsum_primary+DMface*primary_velmaterial

               enddo ! im=1..num_materials

               if (partid_check.ne.nparts) then
                print *,"partid_check invalid"
                stop
               endif

              enddo ! side=1,2

              if (mass_sum.gt.zero) then

               uedge=velsum_primary/mass_sum

              else
               print *,"mass_sum invalid tfrmac"
               print *,"operation_flag=",operation_flag
               print *,"i,j,k,dir ",i,j,k,dir
               stop
              endif

              if (num_colors.eq.0) then
               if (blob_array_size.eq.1) then
                ! do nothing
               else
                print *,"num_colors inconsistent; fort_cell_to_mac"
                stop
               endif
              else if (num_colors.gt.0) then
               if (blob_array_size.eq.num_colors*num_elements_blobclass) then
                ! do nothing 
               else
                print *,"num_colors inconsistent; fort_cell_to_mac"
                stop
               endif
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
               if ((typeleft.ge.1).and.(typeleft.le.num_materials).and. &
                   (typeright.ge.1).and.(typeright.le.num_materials)) then

                 !is_ice_or_FSI_rigid_material=1 if "is_ice" or "is_FSI_rigid"
                if ((is_ice_or_FSI_rigid_material(typeleft).eq.1).or. &
                    (is_ice_or_FSI_rigid_material(typeright).eq.1)) then

                 if ((is_ice_or_FSI_rigid_material(typeleft).eq.1).and. &
                     (is_ice_or_FSI_rigid_material(typeright).eq.1)) then
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
                 else if (is_ice_or_FSI_rigid_material(typeleft).eq.1) then
                  typeface=typeleft
                  colorface=colorleft
                 else if (is_ice_or_FSI_rigid_material(typeright).eq.1) then
                  typeface=typeright
                  colorface=colorright
                 else
                  print *,"typeleft or typeright bust"
                  stop
                 endif

                  ! is_ice==1 or
                  ! is_FSI_rigid==1 
                 if (is_ice_or_FSI_rigid_material(typeface).eq.1) then
                  if ((colorface.ge.1).and.(colorface.le.num_colors)) then
                    ! declared in: GLOBALUTIL.F90
                   call get_rigid_velocity( &
                    FSI_prescribed_flag, &
                    colorface,dir+1,uedge_rigid, &
                    xstenMAC_center, &
                    blob_array, &
                    blob_array_size,num_colors) 

                   if (FSI_flag(typeface).eq.FSI_ICE_STATIC) then
                    uedge_rigid=zero
                   endif

                   call SUB_check_vel_rigid(xstenMAC_center, &
                     time,uedge_rigid,dir+1)

                   if (homogeneous_rigid_velocity.eq.0) then
                    ! do nothing
                   else if (homogeneous_rigid_velocity.eq.1) then
                    uedge_rigid=zero
                   else
                    print *,"homogeneous_rigid_velocity invalid"
                    stop
                   endif

                   uedge=test_current_icemask*uedge+ &
                         (one-test_current_icemask)*uedge_rigid

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

                else if ((is_ice(typeleft).eq.0).and. &
                         (is_ice(typeright).eq.0).and. &
                         (is_FSI_rigid(typeleft).eq.0).and. &
                         (is_FSI_rigid(typeright).eq.0)) then
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
           if (iboundary.eq.1) then

            if ((side.eq.1).or.(side.eq.2)) then
             if (velbc_in(dir+1,side,dir+1).eq.REFLECT_ODD) then
              face_velocity_override=1
              uedge=zero
             else if (velbc_in(dir+1,side,dir+1).eq.EXT_DIR) then
              face_velocity_override=1
              call velbc_override(time,dir+1,side,dir+1, &
               uedge, &
               xstenMAC,nhalf,dx,bfact)

              if (homogeneous_rigid_velocity.eq.0) then
               ! do nothing
              else if (homogeneous_rigid_velocity.eq.1) then
               uedge=zero
              else
               print *,"homogeneous_rigid_velocity invalid"
               stop
              endif

             endif
            else
             print *,"side invalid"
             stop
            endif

           else if (iboundary.eq.0) then

            if (side.eq.0) then
             if ((is_clamped_face.eq.1).or. &
                 (is_clamped_face.eq.2).or. &
                 (is_clamped_face.eq.3)) then
              uedge=vel_clamped(dir+1)
              face_velocity_override=1
             else if (is_clamped_face.eq.0) then
              ! do nothing
             else
              print *,"is_clamped_face invalid"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           else 
            print *,"iboundary invalid"
            stop
           endif

          else
           print *,"at_RZ_face invalid"
           stop
          endif 

          ! local_vel_MAC initialized above with the current MAC velocity
          ! contents.
         else if (operation_flag.eq.OP_PRES_CELL_TO_MAC) then ! p^CELL->MAC

          if ((project_option.eq.SOLVETYPE_PRES).or. &
              (project_option.eq.SOLVETYPE_PRESSTATIC).or. &
              (project_option.eq.SOLVETYPE_PRESGRAVITY).or. &
              (project_option.eq.SOLVETYPE_INITPROJ)) then
           ! do nothing
          else
           print *,"expecting project_option=SOLVETYPE_PRES or "
           print *,"expecting project_option=SOLVETYPE_PRESSTATIC or "
           print *,"expecting project_option=SOLVETYPE_PRESGRAVITY or "
           print *,"expecting project_option=SOLVETYPE_INITPROJ  "
           stop
          endif

          pplus=pres(D_DECL(i,j,k),1)
          pminus=pres(D_DECL(im1,jm1,km1),1)

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

          use_face_pres=1 ! use div(up)

          solid_velocity=local_face(FACECOMP_FACEVEL+1)
          if (project_option.eq.SOLVETYPE_PRESSTATIC) then
           if (solid_velocity.eq.zero) then
            ! do nothing
           else
            print *,"expecting solid_velocity==0.0 if SOLVETYPE_PRESSTATIC"
            stop
           endif
          endif

          face_velocity_override=0

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
             at_ext_wall=side
            else if (velbc_in(dir+1,side,dir+1).eq.INT_DIR) then
             if (mask_coarsefine(side).eq.0) then  
              at_coarse_fine_wallF=side
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
          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if (levelrz.eq.COORDSYS_RZ) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if ((xstenMAC_center(1).le.VOFTOL*dx(1)).and. &
               (dir.eq.0)) then
            if (at_reflect_wall.ne.1) then
             print *,"at_reflect_wall fails sanity check"
             stop
            endif
           endif
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           ! do nothing
          else
           print *,"levelrz invalid grad potential 2"
           stop
          endif

          if ((local_face(FACECOMP_FACECUT+1).ge.zero).and. &
              (local_face(FACECOMP_FACECUT+1).le.half)) then
           AFACE=zero
          else if ((local_face(FACECOMP_FACECUT+1).ge.half).and. &
                   (local_face(FACECOMP_FACECUT+1).le.one)) then
           AFACE=local_face(FACECOMP_FACECUT+1)* &
                 local_face(FACECOMP_ICEFACECUT+1)
          else
           print *,"local_face(FACECOMP_FACECUT+1) invalid"
           stop
          endif

          if ((AFACE.ge.zero).and. &
              (AFACE.le.half)) then
           if (at_reflect_wall.eq.0) then
            ! do nothing
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

          if (local_compressible.eq.0) then
           use_face_pres=0 ! do not use div(up)
          else if (local_compressible.eq.1) then
           ! do nothing
          else
           print *,"local_compressible invalid"
           stop
          endif

           ! local_face(FACECOMP_ICEMASK+1)=zero for both 
           ! ice materials and "is_FSI_rigid" materials.
          if (local_face(FACECOMP_ICEMASK+1).eq.zero) then
           use_face_pres=0 ! do not use div(up)
          else if (local_face(FACECOMP_ICEMASK+1).eq.one) then
           ! do nothing
          else
           print *,"icemask invalid in fort_cell_to_mac"
           print *,"This is the p^CELL->MAC operation"
           print *,"operation_flag (=1) = ",operation_flag
           print *,"FACECOMP_ICEMASK= ",FACECOMP_ICEMASK
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

           ! at least 1 side is covered
          if ((mask_covered(1).eq.0).or. &
              (mask_covered(2).eq.0)) then
           ! do nothing
          else if ((mask_covered(1).eq.1).and. &
                   (mask_covered(2).eq.1)) then
           ! do nothing
          else
           print *,"mask_covered invalid"
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
           do im=1,num_materials
            mass(side)=mass(side)+ &
             local_face(FACECOMP_MASSFACE+2*(im-1)+side)
            vol_local(side)=vol_local(side)+ &
             local_face(FACECOMP_VOFFACE+2*(im-1)+side)
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

           ! 1=use_face_pres=VALID_PEDGE+1 
           ! 2=face pressure=PRESSURE_PEDGE+1
          do im=1,NCOMP_PEDGE
           PEDGE_local(im)=zero
          enddo

          fluid_volface=one
          not_prescribed_volface=one
           ! fixed_face is declared in: PROB.F90
          call fixed_face( &
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
                 (im_prescribed.le.num_materials)) then
              if (im_solid_map(partid_prescribed+1)+1.eq.im_prescribed) then
               use_face_pres=0 ! do not use div(up)
               face_velocity_override=1
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

          if ((is_clamped_face.eq.1).or. &
              (is_clamped_face.eq.2).or. &
              (is_clamped_face.eq.3)) then
           use_face_pres=0 ! do not use div(up)
           face_velocity_override=1
          else if (is_clamped_face.eq.0) then
           ! do nothing
          else
           print *,"is_clamped_face invalid"
           stop
          endif 

          if (at_reflect_wall.eq.0) then
           if (denface.gt.zero) then
            ! do nothing
           else
            print *,"denface must be positive"
            stop
           endif
            ! face pressure
           PEDGE_local(PRESSURE_PEDGE+1)= &
            (den_local(1)*pplus+den_local(2)*pminus)/denface
          else if (at_reflect_wall.eq.1) then ! left wall

            ! face pressure
           PEDGE_local(PRESSURE_PEDGE+1)=pplus

          else if (at_reflect_wall.eq.2) then ! right wall

            ! face pressure
           PEDGE_local(PRESSURE_PEDGE+1)=pminus

          else
           print *,"at reflect wall invalid"
           stop
          endif
          PEDGE_local(VALID_PEDGE+1)=use_face_pres

          if (face_velocity_override.eq.1) then
           local_vel_MAC=solid_velocity
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

         else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 

           ! HYDROSTATIC_PRESDEN_MF is initialized in 
           ! NavierStokes::init_gravity_potentional()
           ! init_gravity_potential() calls fort_init_potential
           ! fort_init_potential is declared in: NAVIERSTOKES_3D.F90
           ! fort_init_potential calls general_hydrostatic_pressure_density
           ! which is declared in PROB.F90.
           ! a default expression:
           !  rhohydro=fort_denconst(1)
           !  phydro=
           !   dt*rhohydro(\vec{g}\cdot\vec{x}+||\vec{xaxis}||^2 Omega^2/2)
           !
           ! hydrostatic pressure (HYDROSTATIC_PRESDEN_MF, 1st component)

          if ((project_option.eq.SOLVETYPE_PRES).or. &
              (project_option.eq.SOLVETYPE_PRESGRAVITY)) then
           static_flag=0
          else if (project_option.eq.SOLVETYPE_PRESSTATIC) then
           static_flag=1
          else
           print *,"project_option invalid"
           stop
          endif

          pplus=pres(D_DECL(i,j,k),1)
          pminus=pres(D_DECL(im1,jm1,km1),1)

           ! hydrostatic density (HYDROSTATIC_PRESDEN_MF, 2nd component)
          dplus=den(D_DECL(i,j,k),1)
          dminus=den(D_DECL(im1,jm1,km1),1)
          if ((dplus.gt.zero).and. &
              (dminus.gt.zero)) then
           ! do nothing
          else
           print *,"hydrostatic density must be positive"
           stop
          endif

          im_gravity=0

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

          enddo ! side=1,2

           ! sanity check
          if (levelrz.eq.COORDSYS_CARTESIAN) then
           ! do nothing
          else if (levelrz.eq.COORDSYS_RZ) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           if ((xstenMAC_center(1).le.VOFTOL*dx(1)).and. &
               (dir.eq.0)) then
            if (at_reflect_wall.ne.1) then
             print *,"at_reflect_wall fails sanity check"
             stop
            endif
           endif
          else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
           ! do nothing
          else
           print *,"levelrz invalid grad potential 2"
           stop
          endif

          AFACE=local_face(FACECOMP_FACECUT+1)
          if ((AFACE.ge.zero).and.(AFACE.le.half)) then
           AFACE=zero
          else if ((AFACE.ge.half).and.(AFACE.le.one)) then
           ! do nothing
          else
           print *,"AFACE invalid"
           stop 
          endif

          AFACE_ICE=local_face(FACECOMP_ICEFACECUT+1)
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
           do im=1,num_materials
            mass(side)=mass(side)+ &
             local_face(FACECOMP_MASSFACE+2*(im-1)+side)
            vol_local(side)=vol_local(side)+ &
             local_face(FACECOMP_VOFFACE+2*(im-1)+side)
           enddo 

           if (mass(side).ge.zero) then
            ! do nothing
           else
            print *,"mass(side) invalid"
            stop
           endif

           if (vol_local(side).ge.zero) then
            ! do nothing
           else
            print *,"vol_local(side) invalid"
            stop
           endif

           if (vol_local(side).gt.zero) then
            den_local(side)=mass(side)/vol_local(side)
           else if (vol_local(side).eq.zero) then
            den_local(side)=zero
           else
            print *,"cannot define den_local(side)"
            stop
           endif

           massface=massface+mass(side)
           volface=volface+vol_local(side)
           denface=denface+den_local(side)
          enddo  ! side=1,2

          ! -(P_H/rho_H)(1/(rho_added rho_H))(rho_H grad rho-rho grad rho_H)
          !    (if gradh_gravity=0, im_left=im_right)
          ! -(P_H/rho_H)grad rho/rho_added (gradh_gravity<>0) 
          incremental_gravity=zero 
          pgrad_gravity=zero  ! grad P_H/rho_H

          pgrad_tension=zero ! -sigma kappa grad H/rho_added

          gradh_gravity=zero
          gradh_tension=zero

           ! fixed_face is declared in: PROB.F90
          call fixed_face( &
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

           ! at_wall==1 if FOEXTRAP or REFLECT_EVEN BC for pressure.
           ! gradh_tension represents (H_{m12,i}-H_{m12,i-1})
           ! m1 = "im material"
           ! m2 = "im_opp material"
          if (at_wall.eq.1) then
           ! do nothing, gradh_tension=0 on a wall
           ! do nothing, gradh_gravity=0 on a wall
           ! do nothing, incremental_gravity=0 on a wall
          else if (at_reflect_wall.eq.1) then
           ! do nothing, gradh_tension=0 at a reflecting wall 
           ! do nothing, gradh_gravity=0 at a reflecting wall 
           ! do nothing, incremental_gravity=0 at a reflecting wall 
          else if (at_reflect_wall.eq.2) then
           ! do nothing, gradh_tension=0 at a reflecting wall 
           ! do nothing, gradh_gravity=0 at a reflecting wall 
           ! do nothing, incremental_gravity=0 at a reflecting wall 
          else if ((at_reflect_wall.eq.0).and. &
                   (at_wall.eq.0)) then

           if (is_solid_face.eq.1) then
            gradh_tension=zero
            gradh_gravity=zero
            incremental_gravity=zero
           else if (is_solid_face.eq.0) then

            if ((is_clamped_face.eq.1).or. &
                (is_clamped_face.eq.2).or. &
                (is_clamped_face.eq.3)) then
             gradh_tension=zero
             gradh_gravity=zero
             incremental_gravity=zero
            else if (is_clamped_face.eq.0) then

             gradh_tension=zero
             gradh_gravity=zero
             incremental_gravity=zero

             ! fluid_interface_tension is declared in: PROB.F90
             ! "merge_levelset" is called inside of "fluid_interface_tension"
             call fluid_interface_tension( &
               xstenMAC_center,time, &
               LSleft,LSright, &
               gradh_tension, &
               im_opp,im, &
               im_left_tension,im_right_tension, &
               static_flag)

             call fluid_interface( &
               LSleft,LSright, &
               gradh_gravity, &
               im_opp_gravity,im_gravity, &
               im_left_gravity,im_right_gravity)

             if (gradh_tension.ne.zero) then

              if (im.lt.im_opp) then
               ! do nothing
              else
               print *,"im or im_opp invalid"
               stop
              endif

              do im_heat=1,num_materials
               tcomp=(im_heat-1)*num_state_material+ENUM_TEMPERATUREVAR+1
               mgoni_temp(im_heat)=half*(mgoni(D_DECL(i,j,k),tcomp)+ &
                mgoni(D_DECL(im1,jm1,km1),tcomp))
              enddo ! im_heat

              if ((project_option.eq.SOLVETYPE_PRES).or. &
                  (project_option.eq.SOLVETYPE_PRESGRAVITY)) then
               call get_user_tension(xstenMAC_center,time, &
                fort_tension,user_tension,mgoni_temp)
              else if (project_option.eq.SOLVETYPE_PRESSTATIC) then
               do iten=1,num_interfaces
                user_tension(iten)=fort_static_tension(iten)
               enddo
              else
               print *,"project_option invalid"
               stop
              endif

              call get_iten(im,im_opp,iten)
              call get_scaled_tension(user_tension(iten),tension_scaled)

               ! in: fort_cell_to_mac, OP_POTGRAD_TO_MAC,
               !     surface tension on MAC grid 
               ! pgrad_tension is combined with pgrad_gravity at the very end.
              pgrad_tension=-dt*tension_scaled*local_face(FACECOMP_CURV+1)* &
                      gradh_tension/hx

              if ((local_face(FACECOMP_FACECUT+1).ge.zero).and. &
                  (local_face(FACECOMP_FACECUT+1).le.half)) then
               pgrad_tension=zero
              else if ((local_face(FACECOMP_FACECUT+1).ge.half).and. &
                       (local_face(FACECOMP_FACECUT+1).le.one)) then
               pgrad_tension=pgrad_tension* &
                local_face(FACECOMP_FACECUT+1)* &
                local_face(FACECOMP_FACEDEN+1)
              else
               print *,"local_face(FACECOMP_FACECUT+1) invalid"
               stop
              endif

             else if (gradh_tension.eq.zero) then
               ! do nothing
             else
              print *,"gradh_tension bust"
              stop
             endif 

             if (gradh_gravity.ne.zero) then

              if (im_gravity.lt.im_opp_gravity) then
               ! do nothing
              else
               print *,"im_gravity or im_opp_gravity invalid"
               stop
              endif

              LSleft_grav=LSleft(im_gravity)-LSleft(im_opp_gravity)
              LSright_grav=LSright(im_gravity)-LSright(im_opp_gravity)
             
              dencomp_im=(im_gravity-1)*num_state_material+1+ENUM_DENVAR 
              dencomp_im_opp=(im_opp_gravity-1)*num_state_material+1+ENUM_DENVAR 
              if ((LSleft_grav.eq.zero).and. &
                  (LSright_grav.eq.zero)) then
               interp_factor=half
               den_im=mgoni(D_DECL(im1,jm1,km1),dencomp_im)   
               den_im_opp=mgoni(D_DECL(i,j,k),dencomp_im_opp)   
              else if (LSleft_grav.eq.zero) then
               interp_factor=zero
               den_im=mgoni(D_DECL(im1,jm1,km1),dencomp_im)   
               den_im_opp=mgoni(D_DECL(i,j,k),dencomp_im_opp)   
              else if (LSright_grav.eq.zero) then
               interp_factor=one
               den_im_opp=mgoni(D_DECL(im1,jm1,km1),dencomp_im_opp)
               den_im=mgoni(D_DECL(i,j,k),dencomp_im)   
              else if (LSleft_grav*LSright_grav.lt.zero) then
               interp_factor=LSleft_grav/(LSleft_grav-LSright_grav)
               if (LSleft_grav.gt.zero) then
                if (gradh_gravity.lt.zero) then
                 den_im=mgoni(D_DECL(im1,jm1,km1),dencomp_im)   
                 den_im_opp=mgoni(D_DECL(i,j,k),dencomp_im_opp)   
                else
                 print *,"gradh_gravity invalid"
                 stop
                endif 
               else if (LSright_grav.gt.zero) then
                if (gradh_gravity.gt.zero) then
                 den_im_opp=mgoni(D_DECL(im1,jm1,km1),dencomp_im_opp)
                 den_im=mgoni(D_DECL(i,j,k),dencomp_im)   
                else
                 print *,"gradh_gravity invalid"
                 stop
                endif 
               else
                print *,"LSleft_grav or LSright_grav invalid"
                stop
               endif

              else
               print *,"LSleft_grav,LSright_grav invalid"
               stop
              endif
              
              pres_H=(one-interp_factor)*pminus+interp_factor*pplus
              den_H=(one-interp_factor)*dminus+interp_factor*dplus
              if (den_H.gt.zero) then
               incremental_gravity=-(pres_H/den_H)* &
                 (den_im-den_im_opp)*gradh_gravity/hx
              else
               print *,"den_H invalid"
               stop
              endif
      
              if ((local_face(FACECOMP_FACECUT+1).ge.zero).and. &
                  (local_face(FACECOMP_FACECUT+1).le.half)) then
               incremental_gravity=zero
              else if ((local_face(FACECOMP_FACECUT+1).ge.half).and. &
                       (local_face(FACECOMP_FACECUT+1).le.one)) then
               incremental_gravity=incremental_gravity* &
                local_face(FACECOMP_FACECUT+1)* &
                local_face(FACECOMP_FACEDEN+1)
              else
               print *,"local_face(FACECOMP_FACECUT+1) invalid"
               stop
              endif

             else if (gradh_gravity.eq.zero) then

              incremental_gravity=zero

              if (im_left_gravity.eq.im_right_gravity) then

               if ((im_left_gravity.ge.1).and. &
                   (im_left_gravity.le.num_materials)) then

                if (is_rigid(im_left_gravity).eq.0) then
                 ! do nothing
                else 
                 print *,"is_rigid(im_left_gravity) invalid"
                 stop
                endif

                dencomp_im= &
                   (im_left_gravity-1)*num_state_material+1+ENUM_DENVAR 
                dencomp_im_opp=dencomp_im
                interp_factor=half
                den_im=mgoni(D_DECL(im1,jm1,km1),dencomp_im)   
                den_im_opp=mgoni(D_DECL(i,j,k),dencomp_im_opp)   
                pres_H=(one-interp_factor)*pminus+interp_factor*pplus
                den_H=(one-interp_factor)*dminus+interp_factor*dplus

                if (den_H.gt.zero) then

                 if (is_compressible_mat(im_left_gravity).eq.1) then
                  incremental_gravity=-(pres_H/den_H)* &
                   (one/den_H)* &
                   (den_H*(den_im_opp-den_im)- &
                    half*(den_im_opp+den_im)*(dplus-dminus))/hx
                 else if (is_compressible_mat(im_left_gravity).eq.0) then
                  incremental_gravity=zero
                 else
                  print *,"is_compressible_mat(im_left_gravity) invalid"
                  stop
                 endif

                 if ((local_face(FACECOMP_FACECUT+1).ge.zero).and. &
                     (local_face(FACECOMP_FACECUT+1).le.half)) then
                  incremental_gravity=zero
                 else if ((local_face(FACECOMP_FACECUT+1).ge.half).and. &
                          (local_face(FACECOMP_FACECUT+1).le.one)) then
                  incremental_gravity=incremental_gravity* &
                    local_face(FACECOMP_FACECUT+1)* &
                    local_face(FACECOMP_FACEDEN+1)
                 else
                  print *,"local_face(FACECOMP_FACECUT+1) invalid"
                  stop
                 endif
                else
                 print *,"den_H invalid"
                 stop
                endif

               else
                print *,"im_left_gravity invalid"
                stop
               endif

              else if (im_left_gravity.ne.im_right_gravity) then
               ! do nothing
              else
               print *,"im_left_gravity,im_right_gravity bust"
               stop
              endif
             else
              print *,"gradh_gravity bust"
              stop
             endif 

            else
             print *,"is_clamped_face invalid"
             stop
            endif

           else
            print *,"is_solid_face invalid 9 ",is_solid_face
            stop
           endif

          else
           print *,"at_reflect_wall or at_wall invalid"
           stop
          endif

           ! p=dt( -|g| z + (1/2)Omega^2 r^2 )
           ! force=grad p=dt( -|g| z^hat + Omega^2 r r^hat )
          cutedge=half*(dplus+dminus)
          if (cutedge.gt.zero) then
           ! do nothing
          else
           print *,"cutedge invalid"
           stop
          endif

           ! hydrostatic pressure gradient on the MAC grid
          pgrad_gravity=(pplus-pminus)/(hx*cutedge)

          ! -dt k (grad p)_MAC (energyflag=SUB_OP_FOR_MAIN)
          ! (grad p)_MAC (energyflag=SUB_OP_FOR_SDC)
         else if (operation_flag.eq.OP_PRESGRAD_MAC) then 
        
           ! xcut=(*localMF[FACE_WEIGHT_MF+dir])[mfi] 
          if (project_option_is_validF(project_option).eq.1) then
           cutedge=xcut(D_DECL(i,j,k),1)  ! e.g. A/rho
          else
           print *,"project_option invalid"
           stop
          endif

          pplus=pres(D_DECL(i,j,k),1)
          pminus=pres(D_DECL(im1,jm1,km1),1)

            ! regular solver or SDC viscosity or thermal flux.
          if (energyflag.eq.SUB_OP_FOR_MAIN) then  
           pgrad=-dt*cutedge*(pplus-pminus)/hx

           ! for SDC pressure gradient
          else if (energyflag.eq.SUB_OP_FOR_SDC) then  
           pgrad=(pplus-pminus)/hx

          else
           print *,"energyflag invalid OP_PRESGRAD_MAC"
           stop
          endif

         else
          print *,"operation_flag invalid17:",operation_flag
          stop
         endif

         if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

          if (ncphys.ne.NFLUXSEM) then
           print *,"ncphys invalid"
           stop
          endif

          do nc=1,ncphys
           xface(D_DECL(i,j,k),nc)=local_face(nc)
          enddo ! nc

         else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & !u^CELL->MAC
                  (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC).or. & 
                  (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                  (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC)) then 

          xvel(D_DECL(i,j,k),1)=uedge
  
         else if (operation_flag.eq.OP_PRESGRAD_MAC) then ! (grad p)_MAC

          xgp(D_DECL(i,j,k),1)=pgrad

         else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 
      
          if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
              POTGRAD_SURFTEN_BASE_GRAV) then 
           xgp(D_DECL(i,j,k),1)=pgrad_gravity+pgrad_tension
          else if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                   POTGRAD_SURFTEN) then 
           xgp(D_DECL(i,j,k),1)=pgrad_tension
          else if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                   POTGRAD_INCREMENTAL_GRAV) then 
           xgp(D_DECL(i,j,k),1)=incremental_gravity
          else if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                   POTGRAD_SURFTEN_INCREMENTAL_GRAV) then 
           xgp(D_DECL(i,j,k),1)=incremental_gravity+pgrad_tension
          else if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                   POTGRAD_BASE_GRAV) then 
           xgp(D_DECL(i,j,k),1)=pgrad_gravity
          else
           print *,"energyflag invalid OP_POTGRAD_TO_MAC"
           stop
          endif

         else if (operation_flag.eq.OP_PRES_CELL_TO_MAC) then !p^CELL->MAC

           ! PEDGE_local(VALID_PEDGE+1)=use_face_pres flag
           ! PEDGE_local(PRESSURE_PEDGE+1)=face pressure
          do im=1,NCOMP_PEDGE
           xp(D_DECL(i,j,k),im)=PEDGE_local(im)
          enddo
          xvel(D_DECL(i,j,k),1)=local_vel_MAC

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

        ! unew^MAC=uSOLID^MAC or uFLUID^MAC
       if (operation_flag.eq.OP_UNEW_USOL_MAC_TO_MAC) then
        ! do nothing
       else if (operation_flag.eq.OP_PRES_CELL_TO_MAC) then 
        ! do nothing ( grad(U P) term only low order )
       else if ((operation_flag.eq.OP_POTGRAD_TO_MAC).and. &
                (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                 POTGRAD_SURFTEN)) then 
        ! do nothing 
       else if ((operation_flag.eq.OP_POTGRAD_TO_MAC).and. &
                (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                 POTGRAD_INCREMENTAL_GRAV)) then 
        ! do nothing 
       else if ((operation_flag.eq.OP_POTGRAD_TO_MAC).and. &
                (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                 POTGRAD_SURFTEN_INCREMENTAL_GRAV)) then 
        ! do nothing 
       else if ((operation_flag.eq.OP_PRESGRAD_MAC).or. & ! pressure gradient
                (operation_flag.eq.OP_POTGRAD_TO_MAC).or. & 
                (operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & ! vel CELL->MAC
                (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC).or.& 
                (operation_flag.eq.OP_ISCHEME_MAC)) then ! advection

        if (operation_flag.eq.OP_POTGRAD_TO_MAC) then

         if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
             POTGRAD_SURFTEN_BASE_GRAV) then
          ! do nothing 
         else if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                  POTGRAD_BASE_GRAV) then
          ! do nothing 
         else
          print *,"spectral method not for incremental formulation"
          stop
         endif

        endif

        if (enable_spectral.eq.1) then

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

             ! local_maskSEM==0 for rigid materials or ice
            if ((local_maskSEM.ge.1).and. &
                (local_maskSEM.le.num_materials).and. &
                (maskcov.eq.1)) then

             if ((operation_flag.eq.OP_PRESGRAD_MAC).or. & ! pressure gradient
                 (operation_flag.eq.OP_POTGRAD_TO_MAC).or. & 
                 (operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & ! vel CELL->MAC
                 (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                 (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC).or.& 
                 (operation_flag.eq.OP_ISCHEME_MAC)) then ! advection
              ! do nothing
             else
              print *,"operation_flag invalid19"
              stop
             endif

             call elementbox(i,j,k,bfact,dir,elemlo,elemhi)
             do ielem=elemlo(1),elemhi(1)
             do jelem=elemlo(2),elemhi(2)
             do kelem=elemlo(3),elemhi(3)

              if (operation_flag.eq.OP_PRESGRAD_MAC) then  ! pressure gradient
               scomp=1
               dcomp=1
               ncomp_dest=1
               ncomp_source=1
               scomp_bc=1

               !grad ppot/den 
              else if (operation_flag.eq.OP_POTGRAD_TO_MAC) then 
               scomp=1
               dcomp=1
               ncomp_dest=1
               ncomp_source=1
               scomp_bc=1

               if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                   POTGRAD_SURFTEN_BASE_GRAV) then
                ! do nothing 
               else if (energyflag.eq.SUB_OP_FORCE_MASK_BASE+ &
                        POTGRAD_BASE_GRAV) then
                ! do nothing 
               else
                print *,"spectral method not for incremental formulation"
                stop
               endif

              else if (operation_flag.eq.OP_UNEW_CELL_TO_MAC) then 

               scomp=dir+1
               dcomp=1
               ncomp_dest=1
               ncomp_source=1
               scomp_bc=dir+1

              else if (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC) then 

               scomp=dir+1
               dcomp=1
               ncomp_dest=1
               ncomp_source=1
               scomp_bc=dir+1

              else if (operation_flag.eq.OP_U_COMP_CELL_MAC_TO_MAC) then 

               scomp=dir+1
               dcomp=1
               ncomp_dest=1
               ncomp_source=1
               scomp_bc=dir+1

              else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection

               scomp=1
               dcomp=1
               ncomp_dest=ncphys
               ncomp_source=SDIM
               scomp_bc=1

              else
               print *,"operation_flag invalid20: ",operation_flag
               stop
              endif
    
              call SEM_CELL_TO_MAC( &
               ncomp_xp, &
               simple_AMR_BC_flag, &
               level, &
               finest_level, &
               operation_flag, & 
               energyflag, &
               project_option, &
               beta, &
               visc_coef, &
               time, &
               dt, &
               ielem,jelem,kelem, &
               tilelo,tilehi, &
               fablo,fabhi, &
               xlo,dx, &
               dir+1, &
               bfact,bfact_c,bfact_f, &
               presbc_in, &
               velbc_in, &
               scomp, &
               scomp_bc, &
               dcomp, &
               ncomp_dest, &
               ncomp_source, &
               ncomp_xgp, &
               ncphys, &
               spectral_loop, &
               ncfluxreg, &
               semflux_ptr, &
               mask_ptr, & !mask=1.0 at interior fine bc ghost cells
               maskcoef_ptr, & ! 1=not cov. or outside domain  
               vel_ptr, &
               pres_ptr, &
               den_ptr, &
               xface_ptr, &
               ! Umac_old if: OP_UMAC_PLUS_VISC_CELL_TO_MAC, or
               !              OP_U_COMP_CELL_MAC_TO_MAC
               xgp_ptr, & 
               xcut_ptr, &   ! coeff*areafrac
               xp_ptr, &
               xvel_ptr, &
               maskSEM_ptr)

             enddo 
             enddo 
             enddo  ! ielem,jelem,kelem

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

        else if (enable_spectral.eq.0) then
         ! do nothing
        else
         print *,"enable_spectral invalid"
         stop
        endif

       else
        print *,"operation_flag invalid"
        stop
       endif

      else
       print *,"tileloop invalid"
       stop
      endif

      return
      end subroutine fort_cell_to_mac

      subroutine fort_project_to_rigid_velocity( &
       dir, &
       velbc_in, &
       slab_step, &
       time, &
       xlo,dx, &
       maskcoef,DIMS(maskcoef), & ! 1=not cov. or outside domain  0=covered
       levelPC,DIMS(levelPC), &
       velMAC,DIMS(velMAC), &
       velCELL,DIMS(velCELL), &
       colorfab,DIMS(colorfab), &
       typefab,DIMS(typefab), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level,finest_level, &
       rz_flag, &
       domlo,domhi, &
       blob_array, &
       blob_array_size, &
       num_colors) &
      bind(c,name='fort_project_to_rigid_velocity')

      use global_utility_module
      use MOF_routines_module
      use probf90_module
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: dir
      INTEGER_T, INTENT(in) :: blob_array_size
      INTEGER_T, INTENT(in) :: num_colors
      REAL_T, INTENT(in) :: blob_array(blob_array_size)
        
      INTEGER_T, INTENT(in) :: slab_step
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, INTENT(in) :: DIMDEC(maskcoef)
      INTEGER_T, INTENT(in) :: DIMDEC(velMAC)
      INTEGER_T, INTENT(in) :: DIMDEC(velCELL)
      INTEGER_T, INTENT(in) :: DIMDEC(levelPC)
      INTEGER_T, INTENT(in) :: DIMDEC(colorfab)
      INTEGER_T, INTENT(in) :: DIMDEC(typefab)

      INTEGER_T, INTENT(in) :: velbc_in(SDIM,2,SDIM)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: rz_flag
      INTEGER_T, INTENT(in) :: domlo(SDIM),domhi(SDIM)

      REAL_T, INTENT(in), target :: maskcoef(DIMV(maskcoef))
      REAL_T, pointer :: maskcoef_ptr(D_DECL(:,:,:))

      REAL_T, INTENT(in), target :: &
              levelPC(DIMV(levelPC),num_materials*(1+SDIM))
      REAL_T, pointer :: levelPC_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout), target :: velMAC(DIMV(velMAC))
      REAL_T, pointer :: velMAC_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(inout), target :: velCELL(DIMV(velCELL))
      REAL_T, pointer :: velCELL_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: typefab(DIMV(typefab))
      REAL_T, pointer :: typefab_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in), target :: colorfab(DIMV(colorfab))
      REAL_T, pointer :: colorfab_ptr(D_DECL(:,:,:))
  
      INTEGER_T i,j,k,ii,jj,kk
      INTEGER_T im1,jm1,km1
      INTEGER_T iboundary
      INTEGER_T im
      INTEGER_T im_left,im_right
      INTEGER_T dir2,side
      INTEGER_T, parameter :: nhalf=3
      REAL_T xstenMAC(-nhalf:nhalf,SDIM)
      REAL_T xstenMAC_center(SDIM)
      REAL_T uedge
      INTEGER_T idx
      INTEGER_T at_RZ_face
      REAL_T LSleft(num_materials)
      REAL_T LSright(num_materials)
      REAL_T localLS(num_materials)
      REAL_T xclamped_minus(SDIM)
      REAL_T xclamped_plus(SDIM)
      REAL_T LS_clamped_plus
      REAL_T LS_clamped_minus
      REAL_T vel_clamped_plus(SDIM)
      REAL_T vel_clamped_minus(SDIM)
      REAL_T temperature_clamped_plus
      REAL_T temperature_clamped_minus
      INTEGER_T prescribed_flag
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped
      INTEGER_T is_clamped_face
      REAL_T xclamped_minus_sten(-nhalf:nhalf,SDIM)
      REAL_T xclamped_plus_sten(-nhalf:nhalf,SDIM)

      INTEGER_T colorface,colorleft,colorright
      INTEGER_T typeface,typeleft,typeright
      INTEGER_T FSI_prescribed_flag !=1 if FSI_RIGID material touches rigid.

      velMAC_ptr=>velMAC
      velCELL_ptr=>velCELL
      maskcoef_ptr=>maskcoef
      typefab_ptr=>typefab
      colorfab_ptr=>colorfab
      levelPC_ptr=>levelPC

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      if ((level.gt.finest_level).or.(level.lt.0)) then
       print *,"level invalid fort_project_to_rigid_velocity"
       stop
      endif
 
      if (num_colors.ge.1) then
       ! do nothing
      else
       print *,"num_colors invalid; fort_project_to_rigid_velocity"
       stop
      endif

      if (blob_array_size.eq.num_colors*num_elements_blobclass) then
       ! do nothing
      else
       print *,"blob_array_size invalid fort_project_to_rigid_velocity"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      do im=1,num_materials

       if (fort_denconst(im).gt.zero) then
        ! do nothing
       else
        print *,"denconst invalid"
        stop
       endif

      enddo ! im=1..num_materials

      if ((slab_step.lt.0).or.(slab_step.ge.bfact_time_order)) then
       print *,"slab_step invalid fort_project_to_rigid_velocity"
       stop
      endif

      call checkbound_array1(fablo,fabhi,velMAC_ptr,0,dir)
      call checkbound_array1(fablo,fabhi,velCELL_ptr,1,-1)

      call checkbound_array1(fablo,fabhi,typefab_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,colorfab_ptr,1,-1)
      call checkbound_array(fablo,fabhi,levelPC_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,maskcoef_ptr,1,-1)

      do im=1,num_materials

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).ge.1).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
        ! do nothing
       else
        print *,"fort_material_type invalid"
        stop
       endif

      enddo  ! im=1..num_materials

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
       print *,"dir out of range in fort_project_to_rigid_velocity, dir=",dir
       stop
      endif 

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       ! dir=0..sdim-1
       call gridstenMAC_level(xstenMAC,i,j,k,level,nhalf,dir)
       do dir2=1,SDIM
        xstenMAC_center(dir2)=xstenMAC(0,dir2)
       enddo

       uedge=velMAC(D_DECL(i,j,k))

       is_clamped_face=-1

       im_left=0
       im_right=0

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        !do nothing 
       else if (levelrz.eq.COORDSYS_RZ) then
        !do nothing 
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        !do nothing 
       else
        print *,"levelrz invalid fort_project_to_rigid_velocity "
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
        print *,"dir invalid fort_project_to_rigid_velocity, dir=",dir
        stop
       endif

       do im=1,num_materials
        LSleft(im)=levelPC(D_DECL(im1,jm1,km1),im)
        LSright(im)=levelPC(D_DECL(i,j,k),im)
        localLS(im)=half*(LSright(im)+LSleft(im))
       enddo
       call get_primary_material(LSleft,im_left)
       call get_primary_material(LSright,im_right)

       call gridsten_level(xclamped_minus_sten,im1,jm1,km1,level,nhalf)
       call gridsten_level(xclamped_plus_sten,i,j,k,level,nhalf)
       do dir2=1,SDIM
        xclamped_minus(dir2)=xclamped_minus_sten(0,dir2)
        xclamped_plus(dir2)=xclamped_plus_sten(0,dir2)
       enddo
       call SUB_clamped_LS(xclamped_minus,time,LS_clamped_minus, &
        vel_clamped_minus,temperature_clamped_minus,prescribed_flag,dx)
       call SUB_clamped_LS(xclamped_plus,time,LS_clamped_plus, &
        vel_clamped_plus,temperature_clamped_plus,prescribed_flag,dx)
       if ((LS_clamped_plus.ge.zero).and. &
           (LS_clamped_minus.ge.zero)) then
        is_clamped_face=1
        do dir2=1,SDIM
         vel_clamped(dir2)=half*(vel_clamped_minus(dir2)+ &
          vel_clamped_plus(dir2))
        enddo
        temperature_clamped=half*(temperature_clamped_minus+ &
          temperature_clamped_plus)
       else if ((LS_clamped_plus.lt.zero).and. &
                (LS_clamped_minus.lt.zero)) then
        is_clamped_face=0
       else if ((LS_clamped_plus.ge.zero).and. &
                (LS_clamped_minus.lt.zero)) then
        is_clamped_face=2
        do dir2=1,SDIM
         vel_clamped(dir2)=vel_clamped_plus(dir2)
        enddo
        temperature_clamped=temperature_clamped_plus
       else if ((LS_clamped_plus.lt.zero).and. &
                (LS_clamped_minus.ge.zero)) then
        is_clamped_face=3
        do dir2=1,SDIM
         vel_clamped(dir2)=vel_clamped_minus(dir2)
        enddo
        temperature_clamped=temperature_clamped_minus
       else
        print *,"LS_clamped_plus or LS_clamped_minus NaN"
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
        if ((dir.eq.0).and. &
            (xstenMAC_center(1).le.VOFTOL*dx(1))) then
         at_RZ_face=1
        endif
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        if ((dir.eq.0).and. &
            (xstenMAC_center(1).le.VOFTOL*dx(1))) then
         at_RZ_face=1
        endif
       else
        print *,"levelrz invalid fort_project_to_rigid_velocity"
        stop
       endif 

       if (at_RZ_face.eq.1) then
        uedge=zero
       else if (at_RZ_face.eq.0) then

        if (num_colors.gt.0) then
         if (blob_array_size.eq.num_colors*num_elements_blobclass) then
          ! do nothing 
         else
          print *,"num_colors inconsistent; fort_project_to_rigid_velocity"
          stop
         endif
        else
         print *,"num_colors invalid"
         stop
        endif

        ! type init in FORT_GETTYPEFAB
        typeleft=NINT(typefab(D_DECL(im1,jm1,km1)))
        typeright=NINT(typefab(D_DECL(i,j,k)))
        colorleft=NINT(colorfab(D_DECL(im1,jm1,km1)))
        colorright=NINT(colorfab(D_DECL(i,j,k)))
        if ((typeleft.ge.1).and.(typeleft.le.num_materials).and. &
            (typeright.ge.1).and.(typeright.le.num_materials)) then

          !is_ice_or_FSI_rigid_material=1 if "is_ice" or "is_FSI_rigid"
         if ((is_ice_or_FSI_rigid_material(typeleft).eq.1).or. &
             (is_ice_or_FSI_rigid_material(typeright).eq.1)) then

          if ((is_ice_or_FSI_rigid_material(typeleft).eq.1).and. &
              (is_ice_or_FSI_rigid_material(typeright).eq.1)) then
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
          else if (is_ice_or_FSI_rigid_material(typeleft).eq.1) then
           typeface=typeleft
           colorface=colorleft
          else if (is_ice_or_FSI_rigid_material(typeright).eq.1) then
           typeface=typeright
           colorface=colorright
          else
           print *,"typeleft or typeright bust"
           stop
          endif

           ! is_ice==1 or
           ! is_FSI_rigid==1 
          if (is_ice_or_FSI_rigid_material(typeface).eq.1) then
           if ((colorface.ge.1).and.(colorface.le.num_colors)) then
             ! declared in: GLOBALUTIL.F90
            call get_rigid_velocity( &
             FSI_prescribed_flag, & !=1 if FSI_rigid material touches rigid.
             colorface,dir+1,uedge, &
             xstenMAC_center, &
             blob_array, &
             blob_array_size,num_colors) 

            if (FSI_flag(typeface).eq.FSI_ICE_STATIC) then
             uedge=zero
            endif

            call SUB_check_vel_rigid(xstenMAC_center, &
              time,uedge,dir+1)

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

         else if ((is_ice(typeleft).eq.0).and. &
                  (is_ice(typeright).eq.0).and. &
                  (is_FSI_rigid(typeleft).eq.0).and. &
                  (is_FSI_rigid(typeright).eq.0)) then
          ! do nothing
         else
          print *,"is_ice or is_FSI_rigid invalid"
          stop
         endif

        else
         print *,"typeleft or typeright invalid"
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

        if (iboundary.eq.1) then

         if ((side.eq.1).or.(side.eq.2)) then
          if (velbc_in(dir+1,side,dir+1).eq.REFLECT_ODD) then
           uedge=zero
          else if (velbc_in(dir+1,side,dir+1).eq.EXT_DIR) then
           call velbc_override(time,dir+1,side,dir+1, &
            uedge, &
            xstenMAC,nhalf,dx,bfact)
          endif
         else
          print *,"side invalid"
          stop
         endif

        else if (iboundary.eq.0) then

         if (side.eq.0) then
          if ((is_clamped_face.eq.1).or. &
              (is_clamped_face.eq.2).or. &
              (is_clamped_face.eq.3)) then
           uedge=vel_clamped(dir+1)
          else if (is_clamped_face.eq.0) then
           ! do nothing
          else
           print *,"is_clamped_face invalid"
           stop
          endif
         else
          print *,"side invalid"
          stop
         endif
        else 
         print *,"iboundary invalid"
         stop
        endif

       else
        print *,"at_RZ_face invalid"
        stop
       endif 

       velMAC(D_DECL(i,j,k))=uedge

      enddo
      enddo
      enddo ! i,j,k (MAC grid, zero ghost cells)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       ! dir=0..sdim-1
       call gridsten_level(xstenMAC,i,j,k,level,nhalf)
       do dir2=1,SDIM
        xstenMAC_center(dir2)=xstenMAC(0,dir2)
       enddo

       uedge=velCELL(D_DECL(i,j,k))

       is_clamped_face=-1

       im_left=0

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        !do nothing 
       else if (levelrz.eq.COORDSYS_RZ) then
        !do nothing 
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        !do nothing 
       else
        print *,"levelrz invalid fort_project_to_rigid_velocity "
        stop
       endif 

       do im=1,num_materials
        localLS(im)=levelPC(D_DECL(i,j,k),im)
       enddo
       call get_primary_material(localLS,im_left)

       call SUB_clamped_LS( &
         xstenMAC_center, & !intent(in)
         time, & !intent(in)
         LS_clamped_minus, & !intent(out)
         vel_clamped,& !intent(out)
         temperature_clamped,& !intent(out)
         prescribed_flag, & !intent(out) 
         dx) !intent(in)

       if (LS_clamped_minus.ge.zero) then
        is_clamped_face=1
       else if (LS_clamped_minus.lt.zero) then
        is_clamped_face=0
       else
        print *," LS_clamped_minus NaN"
        stop
       endif

       if (num_colors.gt.0) then
        if (blob_array_size.eq.num_colors*num_elements_blobclass) then
         ! do nothing 
        else
         print *,"num_colors inconsistent; fort_project_to_rigid_velocity"
         stop
        endif
       else
        print *,"num_colors invalid"
        stop
       endif

        ! type init in FORT_GETTYPEFAB
       typeface=NINT(typefab(D_DECL(i,j,k)))
       colorface=NINT(colorfab(D_DECL(i,j,k)))

       if ((typeface.ge.1).and.(typeface.le.num_materials)) then

         !is_ice_or_FSI_rigid_material=1 if "is_ice" or "is_FSI_rigid"
        if (is_ice_or_FSI_rigid_material(typeface).eq.1) then

         if ((colorface.ge.1).and.(colorface.le.num_colors)) then
          ! declared in: GLOBALUTIL.F90
          call get_rigid_velocity( &
           FSI_prescribed_flag, & !intent(out) (ice touches a substrate)
           colorface, & !intent(in)
           dir+1, & !intent(in)
           uedge, & !intent(out)
           xstenMAC_center, & ! intent(in)
           blob_array, & !intent(in)
           blob_array_size, & !intent(in)
           num_colors) !intent(in)

          if (FSI_flag(typeface).eq.FSI_ICE_STATIC) then
           uedge=zero
          endif

          call SUB_check_vel_rigid(xstenMAC_center, &
            time,uedge,dir+1)

         else if (colorface.eq.0) then
          ! do nothing
         else
          print *,"colorface invalid"
          stop
         endif 

        else if ((is_ice(typeface).eq.0).and. &
                 (is_FSI_rigid(typeface).eq.0)) then
         ! do nothing
        else
         print *,"is_ice or is_FSI_rigid invalid"
         stop
        endif

       else
        print *,"typeface invalid"
        stop
       endif

       if (is_clamped_face.eq.1) then
        uedge=vel_clamped(dir+1)
       else if (is_clamped_face.eq.0) then
        ! do nothing
       else
        print *,"is_clamped_face invalid"
        stop
       endif

       velCELL(D_DECL(i,j,k))=uedge

      enddo
      enddo
      enddo ! i,j,k (CELL grid, zero ghost cells)

      return
      end subroutine fort_project_to_rigid_velocity

      ! called from: NavierStokes::allocate_FACE_WEIGHT (NavierStokes3.cpp)
      !  which is called from:
      !   NavierStokes::update_SEM_forcesALL
      !   NavierStokes::multiphase_project
      !   NavierStokes::diffusion_heatingALL 
      ! mask=1 at fine-fine boundaries
      subroutine fort_buildfacewt( &
       facewt_iter, &
       level, &
       finest_level, &
       nsolve, &
       local_face_index, &
       local_face_ncomp, &
       xlo, &
       dx, &
       offdiagcheck, &
       DIMS(offdiagcheck), &
       cenden,DIMS(cenden), &
       cenvisc,DIMS(cenvisc), &
       xfwt,DIMS(xfwt), &
       yfwt,DIMS(yfwt), &
       zfwt,DIMS(zfwt), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       mask,DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       min_face_wt, & ! static variable
       max_face_wt, & ! static variable
       presbc_arr, &
       visc_coef, &
       uncoupled_viscosity, &
       project_option) &
      bind(c,name='fort_buildfacewt')

      use global_utility_module
      use probcommon_module
      use probf90_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: facewt_iter
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: nsolve
      INTEGER_T, INTENT(in) :: local_face_index
      INTEGER_T, INTENT(in) :: local_face_ncomp
      REAL_T, INTENT(in) :: visc_coef
      INTEGER_T, INTENT(in) :: uncoupled_viscosity
      INTEGER_T, INTENT(in) :: project_option
      REAL_T, INTENT(inout) :: min_face_wt(NCOMP_FACE_WT)
      REAL_T, INTENT(inout) :: max_face_wt(NCOMP_FACE_WT)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: DIMDEC(offdiagcheck)
      INTEGER_T, INTENT(in) :: DIMDEC(cenden)
      INTEGER_T, INTENT(in) :: DIMDEC(cenvisc)
      INTEGER_T, INTENT(in) :: DIMDEC(xfwt)
      INTEGER_T, INTENT(in) :: DIMDEC(yfwt)
      INTEGER_T, INTENT(in) :: DIMDEC(zfwt)
      INTEGER_T, INTENT(in) :: DIMDEC(xface)
      INTEGER_T, INTENT(in) :: DIMDEC(yface)
      INTEGER_T, INTENT(in) :: DIMDEC(zface)
      INTEGER_T, INTENT(in) :: DIMDEC(mask)
      REAL_T, INTENT(in) :: xlo(SDIM),dx(SDIM)
      REAL_T, INTENT(inout),target :: offdiagcheck(DIMV(offdiagcheck),nsolve) 
      REAL_T, pointer :: offdiagcheck_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: cenden(DIMV(cenden)) 
      REAL_T, pointer :: cenden_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(in),target :: cenvisc(DIMV(cenvisc)) 
      REAL_T, pointer :: cenvisc_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(out),target :: xfwt(DIMV(xfwt),nsolve)
      REAL_T, pointer :: xfwt_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out),target :: yfwt(DIMV(yfwt),nsolve)
      REAL_T, pointer :: yfwt_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out),target :: zfwt(DIMV(zfwt),nsolve)
      REAL_T, pointer :: zfwt_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: xface(DIMV(xface),local_face_ncomp)
      REAL_T, pointer :: xface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: yface(DIMV(yface),local_face_ncomp)
      REAL_T, pointer :: yface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: zface(DIMV(zface),local_face_ncomp)
      REAL_T, pointer :: zface_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: mask(DIMV(mask))
      REAL_T, pointer :: mask_ptr(D_DECL(:,:,:))
      INTEGER_T, INTENT(in) :: presbc_arr(SDIM,2,nsolve)
  
      INTEGER_T i,j,k
      INTEGER_T iface,jface,kface
      INTEGER_T ii,jj,kk
      INTEGER_T inorm
      REAL_T dd,dd_group
      REAL_T cc,cc_group
      REAL_T cc_ice
      INTEGER_T side
      INTEGER_T dir
      INTEGER_T veldir
      INTEGER_T velcomp
      REAL_T local_wt(nsolve)
      REAL_T local_fwt
      INTEGER_T local_presbc
      REAL_T local_mask
      INTEGER_T, parameter :: nhalf=3
      REAL_T xsten(-nhalf:nhalf,SDIM)

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level or finest_level invalid build face wt"
       stop
      endif

       ! indexes start at 0
       if ((local_face_index.ge.local_face_ncomp).or. &
           (local_face_index.lt.0)) then
       print *,"local_face_index invalid1"
       print *,"local_face_index=",local_face_index
       print *,"local_face_ncomp=",local_face_ncomp
       stop
      endif
      if (project_option_is_validF(project_option).eq.1) then

       if ((local_face_index.eq.FACECOMP_FACEDEN).or. &
           (local_face_index.eq.FACECOMP_FACEDEN_BASE).or. &
           (local_face_index.eq.FACECOMP_FACEHEAT).or. &
           (local_face_index.eq.FACECOMP_FACEVISC).or. &
           (local_face_index.eq. &
            FACECOMP_FACESPEC+project_option-SOLVETYPE_SPEC)) then
        ! do nothing
       else
        print *,"local_face_index invalid2"
        print *,"local_face_index=",local_face_index
        print *,"local_face_ncomp=",local_face_ncomp
        print *,"FACECOMP_NCOMP=",FACECOMP_NCOMP
        stop
       endif
       if (local_face_ncomp.eq.FACECOMP_NCOMP) then
        ! do nothing
       else
        print *,"local_face_ncomp invalid"
        stop
       endif

      else
       print *,"project_option invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid101"
       stop
      endif

      if ((uncoupled_viscosity.ne.0).and. &
          (uncoupled_viscosity.ne.1)) then
       print *,"uncoupled_viscosity invalid"
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

      offdiagcheck_ptr=>offdiagcheck
      call checkbound_array(fablo,fabhi,offdiagcheck_ptr,0,-1)
      cenden_ptr=>cenden
      call checkbound_array1(fablo,fabhi,cenden_ptr,1,-1)
      cenvisc_ptr=>cenvisc
      call checkbound_array1(fablo,fabhi,cenvisc_ptr,1,-1)
      xfwt_ptr=>xfwt
      yfwt_ptr=>yfwt
      zfwt_ptr=>zfwt
      call checkbound_array(fablo,fabhi,xfwt_ptr,0,0)
      call checkbound_array(fablo,fabhi,yfwt_ptr,0,1)
      call checkbound_array(fablo,fabhi,zfwt_ptr,0,SDIM-1)
      xface_ptr=>xface
      yface_ptr=>yface
      zface_ptr=>zface
      call checkbound_array(fablo,fabhi,xface_ptr,0,0)
      call checkbound_array(fablo,fabhi,yface_ptr,0,1)
      call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1)
      mask_ptr=>mask
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)

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
         print *,"dir out of range in fort_buildfacewt"
         stop
        endif 
  
        call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)

        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)
  
          ! projection: dedge is 1/rho  (FACECOMP_FACEDEN component c++) 
          ! viscosity: dedge is FACECOMP_FACEVISC component c++ ( mu )
          ! temperature: dedge is FACECOMP_FACEHEAT component c++ ( k )
          ! species: dedge is FACECOMP_FACESPEC component c++ ( rho D )

          call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir)

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

          if (project_option_projectionF(project_option).eq.1) then
           !do nothing 
           !SOLVETYPE_PRES,
           !SOLVETYPE_PRESSTATIC,
           !SOLVETYPE_PRESGRAVITY,
           !SOLVETYPE_INITPROJ
          else if (project_option.eq.SOLVETYPE_PRESEXTRAP) then 
           !do nothing
          else if (project_option.eq.SOLVETYPE_HEAT) then ! temperature
           !do nothing
          else if ((project_option.ge.SOLVETYPE_SPEC).and. &
                   (project_option.lt.SOLVETYPE_SPEC+num_species_var)) then 
           !do nothing
          else if (project_option.eq.SOLVETYPE_VISC) then ! viscosity
           !do nothing
          else if (project_option_is_static(project_option).eq.1) then
           !do nothing
          else
           print *,"project_option invalid"
           stop
          endif

          do veldir=1,nsolve

            velcomp=veldir

            if (side.eq.0) then
             ! do nothing
            else if ((side.eq.1).or.(side.eq.2)) then
             local_presbc=presbc_arr(dir+1,side,velcomp)
            else
             print *,"side invalid, side= ",side
             stop
            endif

            if (project_option_is_validF(project_option).eq.1) then

             if (dir.eq.0) then
              cc=xface(D_DECL(i,j,k),FACECOMP_FACECUT+1)
              cc_ice=xface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
             else if (dir.eq.1) then
              cc=yface(D_DECL(i,j,k),FACECOMP_FACECUT+1)
              cc_ice=yface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
             else if ((dir.eq.2).and.(SDIM.eq.3)) then
              cc=zface(D_DECL(i,j,k),FACECOMP_FACECUT+1)
              cc_ice=zface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
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

            else
             print *,"project_option invalid"
             stop
            endif

             ! eval_face_coeff is declared in: PROB.F90
             ! e.g. 1/rho for div( (1/rho) gradp ) = div( ustar )
            call eval_face_coeff( &
             xsten,nhalf, &
             level,finest_level, &
             cc,cc_ice, &
             cc_group, &  ! intent(out)
             dd, &
             dd_group, &  ! intent(out)
             visc_coef, &
             nsolve, &
             dir,veldir, &
             project_option, &
             uncoupled_viscosity, &
             side, &
             local_presbc, &
             local_wt) ! intent(out)

            if (dd_group.lt.min_face_wt(DD_COMP_FACE_WT+1)) then
             min_face_wt(DD_COMP_FACE_WT+1)=dd_group
            endif
            if (cc_group.lt.min_face_wt(CC_COMP_FACE_WT+1)) then
             min_face_wt(CC_COMP_FACE_WT+1)=cc_group
            endif
            if (dd_group.gt.max_face_wt(DD_COMP_FACE_WT+1)) then
             max_face_wt(DD_COMP_FACE_WT+1)=dd_group
            endif
            if (cc_group.gt.max_face_wt(CC_COMP_FACE_WT+1)) then
             max_face_wt(CC_COMP_FACE_WT+1)=cc_group
            endif

            if ((dd_group.ge.zero).and.(cc_group.ge.zero)) then
             ! do nothing
            else
             print *,"cannot have negative wts"
             stop
            endif

            if (local_wt(veldir).lt.min_face_wt(MERGE_COMP_FACE_WT+1)) then
             min_face_wt(MERGE_COMP_FACE_WT+1)=local_wt(veldir)
            endif
            if (local_wt(veldir).gt.max_face_wt(MERGE_COMP_FACE_WT+1)) then
             max_face_wt(MERGE_COMP_FACE_WT+1)=local_wt(veldir)
            endif

            if (DEBUG_THERMAL_COEFF.eq.1) then
             if ((j.eq.32).or.(j.eq.96)) then
              if (i.ge.25) then
               if (project_option.eq.SOLVETYPE_HEAT) then
                print *,"i,j,dir,HEATCOEFF ",i,j,dir,local_wt(veldir)
               endif
              endif
             endif
            endif

            if (dir.eq.0) then
             xfwt(D_DECL(i,j,k),velcomp)=local_wt(veldir)
             xfwt(D_DECL(i,j,k),velcomp)=local_wt(veldir)
            else if (dir.eq.1) then
             yfwt(D_DECL(i,j,k),velcomp)=local_wt(veldir)
             yfwt(D_DECL(i,j,k),velcomp)=local_wt(veldir)
            else if ((dir.eq.2).and.(SDIM.eq.3)) then
             zfwt(D_DECL(i,j,k),velcomp)=local_wt(veldir)
             zfwt(D_DECL(i,j,k),velcomp)=local_wt(veldir)
            else
             print *,"dir invalid buildfacewt 2"
             stop
            endif
          enddo ! veldir=1..nsolve

        enddo ! k
        enddo ! j
        enddo ! i
       enddo ! dir=0,..,sdim-1

      else if (facewt_iter.eq.1) then

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        do veldir=1,nsolve

          velcomp=veldir

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
            offdiagcheck(D_DECL(i,j,k),velcomp)= &
               offdiagcheck(D_DECL(i,j,k),velcomp)+local_fwt
           enddo ! side=1..2
          enddo ! dir=1..sdim
        enddo ! veldir=1,nsolve

       enddo
       enddo
       enddo
      else
       print *,"facewt_iter invalid"
       stop
      endif

      return
      end subroutine fort_buildfacewt

       ! solid: velx,vely,velz,dist  (dist<0 in solid)
       ! called from: NavierStokes::prescribe_solid_geometry
       !   (declared in NavierStokes2.cpp)
      subroutine fort_renormalize_prescribe( &
       tid, &
       level,finest_level, &
       solid_time, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       vofnew,DIMS(vofnew), &
       solxfab,DIMS(solxfab), &
       solyfab,DIMS(solyfab), &
       solzfab,DIMS(solzfab), &
       maskcov,DIMS(maskcov), &
       LS,DIMS(LS), & ! getStateDist(time)
       state_mof,DIMS(state_mof), &
       den,DIMS(den), &
       vel,DIMS(vel), &
       velnew,DIMS(velnew), &
       dennew,DIMS(dennew), &
       lsnew,DIMS(lsnew), &
       xlo,dx, &
       time, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       renormalize_only, &
       solidheat_flag, &
       num_LS_extrap, &
       num_LS_extrap_iter, &
       LS_extrap_iter, &
       ngrow_distance_in, &
       constant_density_all_time) &
      bind(c,name='fort_renormalize_prescribe')
      use global_utility_module
      use global_distance_module
      use probf90_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: tid
      INTEGER_T, INTENT(in) :: solidheat_flag
      INTEGER_T, INTENT(in) :: ngrow_distance_in

      INTEGER_T, INTENT(in) :: renormalize_only
      INTEGER_T, INTENT(inout) :: num_LS_extrap
      INTEGER_T, INTENT(in) :: num_LS_extrap_iter
      INTEGER_T, INTENT(in) :: LS_extrap_iter

      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      REAL_T, INTENT(in) :: solid_time

      REAL_T, INTENT(in) :: xlo(SDIM)
      REAL_T, INTENT(in), target :: dx(SDIM)
      REAL_T, INTENT(in) :: time
      INTEGER_T, INTENT(in) :: nparts
      INTEGER_T, INTENT(in) :: nparts_def
      INTEGER_T, INTENT(in) :: im_solid_map(nparts_def)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: constant_density_all_time(num_materials)

      INTEGER_T, INTENT(in) :: DIMDEC(vofnew)
      INTEGER_T, INTENT(in) :: DIMDEC(solxfab)
      INTEGER_T, INTENT(in) :: DIMDEC(solyfab)
      INTEGER_T, INTENT(in) :: DIMDEC(solzfab)
      INTEGER_T, INTENT(in) :: DIMDEC(maskcov)
      INTEGER_T, INTENT(in) :: DIMDEC(LS)
      INTEGER_T, INTENT(in) :: DIMDEC(state_mof)
      INTEGER_T, INTENT(in) :: DIMDEC(den)
      INTEGER_T, INTENT(in) :: DIMDEC(vel)
      INTEGER_T, INTENT(in) :: DIMDEC(velnew)
      INTEGER_T, INTENT(in) :: DIMDEC(dennew)
      INTEGER_T, INTENT(in) :: DIMDEC(lsnew)
      REAL_T, INTENT(inout),target :: vofnew(DIMV(vofnew),num_materials*ngeom_raw)
      REAL_T, pointer :: vofnew_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: &
           vel(DIMV(vel),STATE_NCOMP_VEL+STATE_NCOMP_PRES)
      REAL_T, pointer :: vel_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(out),target :: &
           velnew(DIMV(velnew),STATE_NCOMP_VEL+STATE_NCOMP_PRES)
      REAL_T, pointer :: velnew_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target ::  &
              LS(DIMV(LS),num_materials*(1+SDIM))
      REAL_T, pointer :: LS_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: state_mof(DIMV(state_mof),num_materials*ngeom_raw)
      REAL_T, pointer :: state_mof_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: den(DIMV(den),num_materials*num_state_material)
      REAL_T, pointer :: den_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout),target :: &
              dennew(DIMV(dennew),num_materials*num_state_material)
      REAL_T, pointer :: dennew_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(inout),target :: lsnew(DIMV(lsnew),num_materials*(1+SDIM))
      REAL_T, pointer :: lsnew_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: solxfab(DIMV(solxfab),SDIM*nparts_def)
      REAL_T, pointer :: solxfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: solyfab(DIMV(solyfab),SDIM*nparts_def)
      REAL_T, pointer :: solyfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: solzfab(DIMV(solzfab),SDIM*nparts_def)
      REAL_T, pointer :: solzfab_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target :: maskcov(DIMV(maskcov))
      REAL_T, pointer :: maskcov_ptr(D_DECL(:,:,:))
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)

      INTEGER_T i,j,k,dir
      INTEGER_T im,im_opp
      INTEGER_T im_primary_stencil
      INTEGER_T im_solid_max
      INTEGER_T vofcomp,vofcompraw
      INTEGER_T i1,j1,k1
      REAL_T centroid(SDIM)
      REAL_T volcell
      REAL_T cencell(SDIM)
      REAL_T censolid_new(num_materials,SDIM)

      INTEGER_T, parameter :: nhalf=9
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T xsten_debug(-nhalf:nhalf,SDIM)

      REAL_T mofnew(num_materials*ngeom_recon)
      INTEGER_T istenlo(3),istenhi(3)
      INTEGER_T LSstenlo(3),LSstenhi(3)
      REAL_T LS_solid_new(num_materials)
      INTEGER_T local_maskcov
      REAL_T vfrac_solid_new(num_materials)
      REAL_T vof_super(num_materials)
      REAL_T F_stencil
      REAL_T F_stencil_sum
      INTEGER_T statecomp
      INTEGER_T statecomp_solid
      INTEGER_T istate,ispecies
      INTEGER_T dencomp,tempcomp,speccomp
      REAL_T den_hold(num_materials*num_state_material)
      REAL_T state_stencil(num_state_material)
      REAL_T state_stencil_sum(num_state_material)
      REAL_T nslope_solid(SDIM)
      INTEGER_T nmax
      REAL_T LS_extend(D_DECL(-1:1,-1:1,-1:1),num_materials)
      REAL_T LS_temp(D_DECL(-1:1,-1:1,-1:1))
      REAL_T LSfacearea
      REAL_T LScentroid(SDIM)
      INTEGER_T nrefine_geom
      REAL_T dxmaxLS
      REAL_T ls_hold(num_materials*(1+SDIM))
      REAL_T max_solid_LS
      REAL_T sum_vfrac_solid_new
      REAL_T LS_predict(num_materials)
      REAL_T LS_virtual(num_materials)
      REAL_T LS_virtual_new(num_materials)
      REAL_T LS_virtual_max ! used for insuring tessellation property of LS
      INTEGER_T num_materials_fluid,num_materials_solid,num_materials_lag
      INTEGER_T at_center
      INTEGER_T ibase
      INTEGER_T partid
      INTEGER_T partid_max
      INTEGER_T tessellate
      INTEGER_T tessellate_transfer
      INTEGER_T LS_extrap_radius
      INTEGER_T extrap_radius
      INTEGER_T least_sqr_radius
      INTEGER_T least_sqrZ
      INTEGER_T center_stencil_im_only
      INTEGER_T center_stencil_wetting_im
      INTEGER_T im1_substencil
      INTEGER_T im2_substencil
      INTEGER_T im_fluid_critical
      INTEGER_T im_local
      INTEGER_T, PARAMETER :: continuous_mof_parm=0
      INTEGER_T cmofsten(D_DECL(-1:1,-1:1,-1:1))

      INTEGER_T :: grid_index(SDIM)
      INTEGER_T, parameter :: grid_level=-1

      REAL_T user_tension(num_interfaces)
      INTEGER_T iten
      REAL_T cos_angle,sin_angle
      REAL_T F_fluid_new
      REAL_T x_fluid_new(SDIM)
      INTEGER_T iten_13,iten_23
      INTEGER_T mof_verbose
      INTEGER_T use_ls_data
      INTEGER_T vofcomprecon
      REAL_T LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)
      REAL_T, DIMENSION(num_materials,SDIM) :: multi_centroidA
      REAL_T orderflag
      REAL_T local_temperature(num_materials)
      REAL_T local_mof(num_materials*ngeom_recon)
      type(cell_CP_parm_type) :: cell_CP_parm
      INTEGER_T cell_index(3)
      REAL_T xCP(SDIM)
      REAL_T xSOLID_BULK(SDIM)
      REAL_T local_XPOS(SDIM)
      REAL_T local_mag
      INTEGER_T, parameter :: nhalf_box=1

      if (renormalize_only.eq.1) then
       if (num_LS_extrap_iter.eq.1) then
        ! do nothing
       else
        print *,"num_LS_extrap_iter invalid"
        stop
       endif
      else if (renormalize_only.eq.0) then
       if (num_LS_extrap_iter.ge.2) then
        ! do nothing
       else
        print *,"num_LS_extrap_iter invalid"
        stop
       endif
      else
       print *,"renormalize_only invalid"
       stop
      endif

      if ((LS_extrap_iter.ge.0).and. &
          (LS_extrap_iter.lt.num_LS_extrap_iter)) then
      ! do nothing
      else
       print *,"LS_extrap_iter invalid"
       stop
      endif
      if (num_LS_extrap.ge.0) then
       ! do nothing
      else
       print *,"num_LS_extrap invalid"
       stop
      endif

      tessellate=0

      nmax=POLYGON_LIST_MAX  ! in: fort_renormalize_prescribe
      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid102"
       stop
      endif 

      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if (ngrow_distance.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance invalid"
       stop
      endif
      if (ngrow_distance_in.eq.4) then
       ! do nothing
      else
       print *,"ngrow_distance_in invalid"
       stop
      endif

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

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_renormalize_prescribe"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_renormalize_prescribe"
       stop
      endif

      extrap_radius=1
      LS_extrap_radius=1

      least_sqr_radius=1
      least_sqrZ=0
      if (SDIM.eq.2) then
       ! do nothing
      else if (SDIM.eq.3) then
       least_sqrZ=1
      else
       print *,"dimension bust"
       stop
      endif

      cell_CP_parm%least_sqrZ=least_sqrZ
      cell_CP_parm%least_sqr_radius=least_sqr_radius
      cell_CP_parm%dxmaxLS=dxmaxLS
      cell_CP_parm%bfact=bfact
      cell_CP_parm%level=level
      cell_CP_parm%finest_level=finest_level
      cell_CP_parm%fablo=>fablo
      cell_CP_parm%fabhi=>fabhi
      cell_CP_parm%dx=>dx
      cell_CP_parm%time=time
      cell_CP_parm%LS=>LS

      num_materials_fluid=0
      num_materials_solid=0
      num_materials_lag=0

      do im=1,num_materials

       if (num_state_material.ne. &
           num_state_base+num_species_var) then
        print *,"num_state_material invalid"
        stop
       endif

       if (is_lag_part(im).eq.1) then
        num_materials_lag=num_materials_lag+1
        if (is_rigid(im).eq.1) then
         num_materials_solid=num_materials_solid+1
        else if (is_rigid(im).eq.0) then
         num_materials_fluid=num_materials_fluid+1
        else
         print *,"is_rigid(im) invalid"
         stop
        endif
       else if (is_lag_part(im).eq.0) then
        if (is_rigid(im).eq.0) then
         num_materials_fluid=num_materials_fluid+1
        else
         print *,"is_rigid(im) invalid"
         stop
        endif
       else
        print *,"is_lag_part(im) invalid"
        stop
       endif

      enddo ! im=1..num_materials

      if (num_materials_lag.ne.nparts) then
       print *,"num_materials_lag invalid"
       stop
      endif

      if (num_materials_fluid+num_materials_solid.ne.num_materials) then
       print *,"num_materials_fluid and/or num_materials_solid invalid"
       stop
      endif

      vofnew_ptr=>vofnew
      call checkbound_array(fablo,fabhi,vofnew_ptr,1,-1)
      velnew_ptr=>velnew
      call checkbound_array(fablo,fabhi,velnew_ptr,1,-1)

      solxfab_ptr=>solxfab
      solyfab_ptr=>solyfab
      solzfab_ptr=>solzfab
      call checkbound_array(fablo,fabhi,solxfab_ptr,0,0)
      call checkbound_array(fablo,fabhi,solyfab_ptr,0,1)
      call checkbound_array(fablo,fabhi,solzfab_ptr,0,SDIM-1)

      maskcov_ptr=>maskcov
      call checkbound_array1(fablo,fabhi,maskcov_ptr,0,-1)
      LS_ptr=>LS
      call checkbound_array(fablo,fabhi,LS_ptr,ngrow_distance,-1)
      state_mof_ptr=>state_mof
      call checkbound_array(fablo,fabhi,state_mof_ptr,1,-1)
      vel_ptr=>vel
      call checkbound_array(fablo,fabhi,vel_ptr,1,-1)
      den_ptr=>den
      call checkbound_array(fablo,fabhi,den_ptr,1,-1)
      dennew_ptr=>dennew
      call checkbound_array(fablo,fabhi,dennew_ptr,1,-1)
      lsnew_ptr=>lsnew
      call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)

      istenlo(3)=0
      istenhi(3)=0
      do dir=1,SDIM
       istenlo(dir)=-extrap_radius
       istenhi(dir)=extrap_radius
      enddo

      LSstenlo(3)=0
      LSstenhi(3)=0
      do dir=1,SDIM
       LSstenlo(dir)=-LS_extrap_radius
       LSstenhi(dir)=LS_extrap_radius
      enddo

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_maskcov=NINT(maskcov(D_DECL(i,j,k)))

       grid_index(1)=i
       grid_index(2)=j
       if (SDIM.eq.3) then
        grid_index(SDIM)=k
       endif

       call gridsten_level(xsten,i,j,k,level,nhalf)

       call Box_volumeFAST(bfact,dx,xsten,nhalf,volcell,cencell,SDIM)

       if (local_maskcov.eq.1) then

        ! --------------------------------------------------------- 
        ! first: fluid state variable extrapolation into empty cells.
        ! ----------------------------------------------------------

        if (LS_extrap_iter.eq.0) then

         do im=1,num_materials

          vofcompraw=(im-1)*ngeom_raw+1
          F_stencil=state_mof(D_DECL(i,j,k),vofcompraw)

          if (is_rigid(im).eq.1) then
           if (constant_density_all_time(im).eq.1) then
            ! do nothing
           else
            print *,"constant_density_all_time(im) invalid"
            print *,"in: fort_renormalize_prescribe"
            stop
           endif
          else if (is_rigid(im).eq.0) then

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

             F_stencil=state_mof(D_DECL(i+i1,j+j1,k+k1),vofcompraw)

              ! in: subroutine fort_renormalize_prescribe

             istate=1
             do while (istate.le.num_state_material)

              if (istate.eq.1+ENUM_DENVAR) then
               dencomp=(im-1)*num_state_material+1+ENUM_DENVAR

               if (constant_density_all_time(im).eq.1) then 
                state_stencil(istate)=fort_denconst(im)
                if (fort_material_type(im).eq.0) then  ! incompressible
                 ! do nothing
                else
                 print *,"fort_material_type invalid"
                 stop
                endif
               else if (constant_density_all_time(im).eq.0) then 
                state_stencil(istate)=den(D_DECL(i+i1,j+j1,k+k1),dencomp)
                if ((fort_material_type(im).ge.0).and. &
                    (fort_material_type(im).le.MAX_NUM_EOS)) then 
                 ! do nothing
                else
                 print *,"fort_material_type invalid"
                 stop
                endif
               else
                print *,"constant_density_all_time(im) invalid"
                print *,"fort_renormalize_prescribe (2)"
                stop
               endif
               istate=istate+1
              else if (istate.eq.1+ENUM_TEMPERATUREVAR) then
               tempcomp=(im-1)*num_state_material+1+ENUM_TEMPERATUREVAR
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
           print *,"is_rigid invalid LEVELSET_3D.F90"
           stop
          endif

         enddo ! im=1..num_materials (extrapolation loop)

        else if (LS_extrap_iter.gt.0) then
         ! do nothing
        else
         print *,"LS_extrap_iter invalid"
         stop
        endif

        ! --------------------------------------------------------- 
        ! end: fluid state variable extrapolation into empty cells.
        ! ----------------------------------------------------------
      
        do dir=1,num_materials*ngeom_recon
         mofnew(dir)=zero
        enddo

        do im=1,num_materials

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

        enddo  ! im=1..num_materials


         ! 1. prescribe solid materials. (F,X,LS,velocity,temperature)
         ! 2. extend the fluid level set functions into the solids.
         !   (F,X,LS fluid)
        if (renormalize_only.eq.0) then

         do im=1,num_materials*(1+SDIM)
          ls_hold(im)=lsnew(D_DECL(i,j,k),im)
         enddo

         sum_vfrac_solid_new=zero
         max_solid_LS=-99999.0
         im_solid_max=0
         partid_max=0

         if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
          print *,"nparts invalid fort_renormalize_prescribe"
          stop
         endif

         do partid=1,nparts

          im=im_solid_map(partid)+1
          if ((im.lt.1).or.(im.gt.num_materials)) then
           print *,"im invalid33"
           stop
          endif

          if (is_lag_part(im).eq.0) then
           print *,"is_lag_part(im).eq.0"
           stop
          else if (is_lag_part(im).eq.1) then
           if (is_rigid(im).eq.0) then
            ! do nothing
           else if (is_rigid(im).eq.1) then

            ! positive in the rigid body
            call materialdistsolid( &
             xsten(0,1),xsten(0,2),xsten(0,SDIM), &
             LS_solid_new(im),time,im)

            if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
                (FSI_flag(im).eq.FSI_SHOELE_PRESVEL).or. & 
                (FSI_flag(im).eq.FSI_SHOELE_VELVEL)) then 
             LS_solid_new(im)=LS(D_DECL(i,j,k),im)
            else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
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

            ! find_LS_stencil_volume is in: PROB.F90, and
            ! calls find_LS_stencil_volume_coarse, which
            ! calls materialdistsolid many times.
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

            if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
                (FSI_flag(im).eq.FSI_SHOELE_PRESVEL).or. & 
                (FSI_flag(im).eq.FSI_SHOELE_VELVEL)) then 
             vofcompraw=(im-1)*ngeom_raw+1
             vfrac_solid_new(im)=vofnew(D_DECL(i,j,k),vofcompraw)
             do dir=1,SDIM
              centroid(dir)=vofnew(D_DECL(i,j,k),vofcompraw+dir)+cencell(dir)
             enddo
            else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
             ! do nothing
            else
             print *,"FSI_flag invalid"
             stop
            endif

            do dir=1,SDIM
             censolid_new(im,dir)=centroid(dir)-cencell(dir)
            enddo

            ! xCP=x-phi grad phi   grad phi=(x-xCP)/phi
            ! slope points into the solid
            ! this slope ignores 1/R term for dphi/dtheta
            call find_LS_stencil_slope( &
             bfact, &
             dx, &
             xsten,nhalf, &
             nslope_solid, &
             time,im)

            if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
                (FSI_flag(im).eq.FSI_SHOELE_PRESVEL).or. & 
                (FSI_flag(im).eq.FSI_SHOELE_VELVEL)) then 
             do dir=1,SDIM
              nslope_solid(dir)=LS(D_DECL(i,j,k),num_materials+SDIM*(im-1)+dir)
             enddo
            else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
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
             ls_hold(num_materials+SDIM*(im-1)+dir)=nslope_solid(dir)
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
              if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
                  (FSI_flag(im).eq.FSI_SHOELE_PRESVEL).or. & 
                  (FSI_flag(im).eq.FSI_SHOELE_VELVEL)) then 
               ! den_hold(statecomp) already has the solid temperature
              else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
               call tempsolid(xsten(0,1),xsten(0,2),xsten(0,SDIM), &
                den_hold(statecomp),time,im)
              else
               print *,"FSI_flag invalid"
               stop
              endif
             else if (solidheat_flag.eq.1) then ! dirichlet at solid/fluid
              if ((FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. & 
                  (FSI_flag(im).eq.FSI_SHOELE_PRESVEL).or. & 
                  (FSI_flag(im).eq.FSI_SHOELE_VELVEL)) then 
               ! den_hold(statecomp) already has the solid temperature
              else if (FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90) then 
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
            print *,"is_rigid invalid LEVELSET_3D.F90"
            stop
           endif
          else 
           print *,"is_lag_part(im) invalid"
           stop
          endif

         enddo ! partid=1..nparts

           ! velocity and temperature in rigid solid.
         if ((sum_vfrac_solid_new.ge.half).or. &
             (max_solid_LS.ge.zero)) then 

           ! number of cells in the solid region
           ! corrected with an extrapolated value from 
           ! the fluid region.
          num_LS_extrap=num_LS_extrap+1 

          if ((im_solid_max.lt.1).or. &
              (im_solid_max.gt.num_materials).or. &
              (partid_max.lt.1).or. &
              (partid_max.gt.nparts).or. &
              (is_rigid(im_solid_max).ne.1)) then
           print *,"im_solid_max or partid_max became corrupt"
           stop
          endif

           ! solid velocity
          ibase=(partid_max-1)*SDIM

           ! velocity
          if (is_prescribed(im_solid_max).eq.1) then

           dir=1
           velnew(D_DECL(i,j,k),dir)= &
              half*(solxfab(D_DECL(i,j,k),ibase+dir)+ &
                    solxfab(D_DECL(i+1,j,k),ibase+dir))
           dir=2
           velnew(D_DECL(i,j,k),dir)= &
              half*(solyfab(D_DECL(i,j,k),ibase+dir)+ &
                    solyfab(D_DECL(i,j+1,k),ibase+dir))
           if (SDIM.eq.3) then
            dir=SDIM
            velnew(D_DECL(i,j,k),dir)= &
              half*(solzfab(D_DECL(i,j,k),ibase+dir)+ &
                    solzfab(D_DECL(i,j,k+1),ibase+dir))
           endif

          else if (is_prescribed(im_solid_max).eq.0) then
           ! do nothing
          else
           print *,"is_prescribed(im_solid_max) invalid"
           stop
          endif

           ! initialize the fluid temperature with the solid temperature
           ! in the solid regions.
          do im=1,num_materials
           if (is_rigid(im).eq.1) then
            ! do nothing
           else if (is_rigid(im).eq.0) then
            istate=2
            statecomp=(im-1)*num_state_material+istate
            statecomp_solid=(im_solid_max-1)*num_state_material+istate
            den_hold(statecomp)=den_hold(statecomp_solid)
           else
            print *,"is_rigid invalid LEVELSET_3D.F90"
            stop
           endif
          enddo ! im=1..num_materials

         else if ((sum_vfrac_solid_new.le.half).and. &
                  (max_solid_LS.le.zero)) then
          ! do nothing
         else
          print *,"sum_vfrac_solid_new or max_solid_LS invalid"
          stop
         endif

         do im=1,num_materials
          do istate=1,num_state_material
           statecomp=(im-1)*num_state_material+istate
           dennew(D_DECL(i,j,k),statecomp)=den_hold(statecomp)
          enddo
         enddo

          ! solid volume fractions and centroids
         do im=1,num_materials
          vofcomp=(im-1)*ngeom_recon+1
          if (is_rigid(im).eq.0) then
           ! do nothing
          else if (is_rigid(im).eq.1) then
           mofnew(vofcomp)=vfrac_solid_new(im)
           do dir=1,SDIM
            mofnew(vofcomp+dir)=censolid_new(im,dir)
           enddo
           do istate=SDIM+2,ngeom_recon
            mofnew(vofcomp+istate-1)=zero
           enddo
          else
           print *,"is_rigid invalid LEVELSET_3D.F90"
           stop
          endif
         enddo ! im=1..num_materials

          ! extend fluid LS,F,X into the solid.
         if ((im_solid_max.ge.1).and. &
             (im_solid_max.le.num_materials)) then

           ! (i,j,k) is a "solid" cell
           ! ls_hold will be modified
          if ((LS_solid_new(im_solid_max).ge.zero).or. &
              (sum_vfrac_solid_new.ge.half)) then

           center_stencil_im_only=0
           center_stencil_wetting_im=0

           im1_substencil=0
           im2_substencil=0

           cell_CP_parm%im_solid_max=im_solid_max

            ! inner loop is needed since the volume fraction
            ! at (i,j,k) depends on the levelset function values
            ! in the (i+i1,j+j1,k+k1) node stencil.
           do i1=LSstenlo(1),LSstenhi(1)
           do j1=LSstenlo(2),LSstenhi(2)
           do k1=LSstenlo(3),LSstenhi(3)

            cell_CP_parm%i=i+i1
            cell_CP_parm%j=j+j1
            cell_CP_parm%k=k+k1
             ! xSOLID_BULK(dir)=xsten(i+i1,j+j1,k+k1,dir)
             ! xCP=xSOLID_BULK(dir)-LS_cell*nslope_cell(dir)
            call cell_xCP(cell_CP_parm,xCP,xSOLID_BULK)

            call containing_cell(bfact, &
              dx, &
              xlo, &
              fablo, &
              xCP, &
              cell_index)

             ! im1_substencil or im2_substencil is initialized
             ! for the fluid material in which 
             ! LS_FLUID_XCP_stencil>=0.0 and
             ! |XCP_stencil-xSOLID_BULK| is a minimum.
             ! (restricted to the radius=1 stencil about cell_index)

             !find primary fluid closest to xSOLID_BULK and 
             !within the radius 1 stencil about cell_index
             !LS_solid(<found stencil closest>)<0
             !LS_FLUID(<found stencil closest>)>0
            call interp_fluid_LS( &
             cell_CP_parm, &
             xCP, &
             xSOLID_BULK, &
             cell_index, & ! containing cell for xCP
             LS_virtual_new, & 
             im_solid_max, &
             im_fluid_critical) !primary fluid closest to xSOLID_BULK and 
                                !within the radius 1 stencil about cell_index
                                !LS_solid(<found stencil>)<0
                                !LS_FLUID(<found stencil>)>0


            if ((i1.eq.0).and. &
                (j1.eq.0).and. &
                (k1.eq.0)) then
             at_center=1
            else
             at_center=0
            endif

            do im=1,num_materials
             LS_predict(im)=LS(D_DECL(i+i1,j+j1,k+k1),im)
            enddo
            call get_primary_material(LS_predict,im_primary_stencil)

             !fluid stencil cell, we trust this LS value.
             !if (at_center==1) then cell is (i,j,k) cell which is solid.
            if ((is_rigid(im_primary_stencil).eq.0).and. & 
                (at_center.eq.0)) then

             do im=1,num_materials
              LS_virtual_new(im)=LS_predict(im)
              LS_extend(D_DECL(i1,j1,k1),im)=LS_virtual_new(im)
             enddo

             !solid stencil cell, extrapolate from the fluid side.
            else if ((is_rigid(im_primary_stencil).eq.1).or. & 
                     (at_center.eq.1)) then

             do im=1,num_materials
              LS_extend(D_DECL(i1,j1,k1),im)=LS_virtual_new(im)
             enddo

             if (im_fluid_critical.eq.0) then
              ! do nothing
             else if ((im_fluid_critical.ge.1).and. &
                      (im_fluid_critical.le.num_materials)) then
              if (is_rigid(im_fluid_critical).eq.0) then
               if (im1_substencil.eq.0) then
                im1_substencil=im_fluid_critical
               else if ((im1_substencil.ge.1).and. &
                        (im1_substencil.le.num_materials)) then
                if (im_fluid_critical.eq.im1_substencil) then
                 ! do nothing
                else
                 im2_substencil=im_fluid_critical
                endif
               else
                print *,"im1_substencil invalid"
                stop
               endif
              else
               print *,"is_rigid(im_fluid_critical) invalid"
               stop
              endif
             else
              print *,"im_fluid_critical invalid"
              stop
             endif

            else
             print *,"is_rigid(im_primary_stencil) or at_center invalid"
             stop
            endif
                           
           enddo
           enddo
           enddo ! i1,j1,k1=LSstenlo ... LSstenhi


           if (im1_substencil.eq.0) then

            if (abs(LS_solid_new(im_solid_max)).le.VOFTOL*dxmaxLS) then
             print *,"all materials disappeared in fort_renormalize_prescribe?"
             print *,"abs(LS_solid_new(im_solid_max)) very small, ", &
              abs(LS_solid_new(im_solid_max))
             print *,"but yet no negative values for "
             print *,"LS_solid_new(im_solid_max) were found nearby."
             print *,"level,finest_level ",level,finest_level
             cell_CP_parm%i=i
             cell_CP_parm%j=j
             cell_CP_parm%k=k
              ! xSOLID_BULK(dir)=xsten(i,j,k,dir)
              ! xCP=xSOLID_BULK(dir)-LS_cell*nslope_cell(dir)
             call cell_xCP(cell_CP_parm,xCP,xSOLID_BULK)
             local_mag=zero
             do dir=1,SDIM
              print *,"dir,xCP,xSOLID_BULK ",dir,xCP(dir),xSOLID_BULK(dir)
              local_mag=local_mag+(xCP(dir)-xSOLID_BULK(dir))**2
             enddo
             local_mag=sqrt(local_mag)
             print *,"im_solid_max,LS_solid_new(im_solid_max) ", &
                     im_solid_max,LS_solid_new(im_solid_max)
             print *,"dx,dy,dz,local_mag ",dx(1),dx(2),dx(SDIM),local_mag

             do im=1,num_materials
              print *,"im,is_rigid(im),ls_hold(im) ", &
                      im,is_rigid(im),ls_hold(im)
             enddo
             do i1=LSstenlo(1),LSstenhi(1)
             do j1=LSstenlo(2),LSstenhi(2)
             do k1=LSstenlo(3),LSstenhi(3)
              call gridsten_level(xsten_debug,i+i1,j+j1,k+k1,level,nhalf)
              print *,"i1,j1,k1,x,y,z ",i1,j1,k1, &
                 xsten_debug(0,1),xsten_debug(0,2),xsten_debug(0,SDIM)
              do im=1,num_materials
               print *,"i1,j1,k1,im,LS_extend ",i1,j1,k1,im, &
                       LS_extend(D_DECL(i1,j1,k1),im)
              enddo
             enddo
             enddo
             enddo

             stop

            else if (abs(LS_solid_new(im_solid_max)).ge.VOFTOL*dxmaxLS) then
             ! do nothing
            else
             print *,"LS_solid_new(im_solid_max) invalid"
             stop
            endif
           else if ((im1_substencil.ge.1).and. &
                    (im1_substencil.le.num_materials)) then
            if (im2_substencil.eq.0) then
             !fluid: center_stencil_im_only owns the whole cell
             center_stencil_im_only=im1_substencil 
            else if ((im2_substencil.ge.1).and. &
                     (im2_substencil.le.num_materials)) then
             if (im1_substencil.lt.im2_substencil) then
              im=im1_substencil
              im_opp=im2_substencil
             else if (im1_substencil.gt.im2_substencil) then
              im=im2_substencil
              im_opp=im1_substencil
             else
              print *,"im1_substencil or im2_substencil invalid"
              stop
             endif
             call get_iten(im,im_opp,iten)
             do im_local=1,num_materials
              dencomp=(im_local-1)*num_state_material+1+ENUM_DENVAR
              local_temperature(im_local)=den(D_DECL(i,j,k),dencomp+1)
             enddo
              ! coordinate of (i,j,k)
             do dir=1,SDIM
              local_XPOS(dir)=xsten(0,dir)
             enddo
             call get_user_tension( &
              local_XPOS, &
              time, &
              fort_tension,user_tension, &
              local_temperature)
             ! sigma_{i,j}cos(theta_{i,k})=sigma_{j,k}-sigma_{i,k}
             ! theta_{ik}=0 => material i wets material k.
             ! im is material "i"  ("fluid" material)
             ! im_opp is material "j"
             call get_CL_iten(im,im_opp,im_solid_max, &
              iten_13,iten_23, &
              user_tension,cos_angle,sin_angle)
             if ((sin_angle.eq.zero).and.(cos_angle.eq.one)) then
              center_stencil_wetting_im=im
             else if ((sin_angle.eq.zero).and.(cos_angle.eq.-one)) then 
              center_stencil_wetting_im=im_opp
             else if ((sin_angle.ge.zero).and. &
                      (sin_angle.le.one).and. &
                      (cos_angle.ge.-one).and. &
                      (cos_angle.le.one)) then
              ! do nothing
             else
              print *,"sin_angle or cos_angle invalid"
              stop
             endif
            else
             print *,"im2_substencil invalid"
             stop
            endif
           else
            print *,"im1_substencil invalid"
            stop
           endif

            ! default radius: extrap_radius=1 cell
            ! make the extension fluid level set tessellating
           do i1=istenlo(1),istenhi(1)
           do j1=istenlo(2),istenhi(2)
           do k1=istenlo(3),istenhi(3)

            do im=1,num_materials
             LS_virtual(im)=LS_extend(D_DECL(i1,j1,k1),im)
             LS_virtual_new(im)=LS_virtual(im)
            enddo
 
            do im=1,num_materials
             if (is_rigid(im).eq.1) then
              ! do nothing
             else if (is_rigid(im).eq.0) then

              LS_virtual_max=-99999.0
              im_primary_stencil=0
              do im_opp=1,num_materials
               if (is_rigid(im_opp).eq.1) then
                ! do nothing
               else if (is_rigid(im_opp).eq.0) then
                if (im.ne.im_opp) then
                 if (im_primary_stencil.eq.0) then
                  im_primary_stencil=im_opp
                  LS_virtual_max=LS_virtual(im_opp)
                 else if ((im_primary_stencil.ge.1).and. &
                          (im_primary_stencil.le.num_materials)) then
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
                print *,"is_rigid invalid LEVELSET_3D.F90"
                stop
               endif
              enddo ! im_opp=1..num_materials
 
              if (im_primary_stencil.eq.0) then
               if (num_materials_fluid.ne.1) then
                print *,"num_materials_fluid invalid"
                stop
               endif
               if (LS_virtual(im).le.zero) then
                print *,"LS_virtual(im) invalid"
                stop
               endif
              else if ((im_primary_stencil.ge.1).and. &
                       (im_primary_stencil.le.num_materials)) then
               LS_virtual_new(im)=half*(LS_virtual(im)-LS_virtual_max)
              else
               print *,"im_primary_stencil invalid"
               stop
              endif

             else
              print *,"is_rigid invalid LEVELSET_3D.F90"
              stop
             endif

            enddo ! im=1..num_materials

            if ((center_stencil_wetting_im.ge.1).and. &
                (center_stencil_wetting_im.le.num_materials)) then 
             if (LS_virtual_new(im_solid_max).ge.zero) then
              do im=1,num_materials
               if (is_rigid(im).eq.0) then
                if (im.eq.center_stencil_wetting_im) then
                 if (LS_virtual_new(im).lt.zero) then
                  LS_virtual_new(im)=zero
                 else if (LS_virtual_new(im).ge.zero) then
                  ! do nothing
                 else
                  print *,"LS_virtual_new invalid"
                  stop
                 endif
                else
                 if (LS_virtual_new(im).ge.zero) then
                  LS_virtual_new(im)=zero
                 else if (LS_virtual_new(im).lt.zero) then
                  ! do nothing
                 else
                  print *,"LS_virtual_new invalid"
                  stop
                 endif
                endif
               else if (is_rigid(im).eq.1) then
                ! do nothing
               else
                print *,"is_rigid(im) invalid"
                stop
               endif
              enddo ! im=1..num_materials
             else if (LS_virtual_new(im_solid_max).lt.zero) then
              ! do nothing
             else
              print *,"LS_virtual_new(im_solid_max) invalid"
              stop
             endif
            else if (center_stencil_wetting_im.eq.0) then
             ! do nothing
            else
             print *,"center_stencil_wetting_im invalid"
             stop
            endif

            do im=1,num_materials
             LS_extend(D_DECL(i1,j1,k1),im)=LS_virtual_new(im)
            enddo

           enddo
           enddo
           enddo ! i1,j1,k1=-extend_radius..extend_radius

           if ((center_stencil_wetting_im.ge.1).and. &
               (center_stencil_wetting_im.le.num_materials)) then 

            if (1.eq.0) then
             print *,"center_stencil_wetting_im=", &
                center_stencil_wetting_im
             print *,"i,j,k ",i,j,k
             print *,"level,finest_level ",level,finest_level
            endif

            use_ls_data=0

            do im=1,num_materials
             vofcomprecon=(im-1)*ngeom_recon+1
             vofcompraw=(im-1)*ngeom_raw+1

             if (is_rigid(im).eq.0) then
              do dir=0,SDIM
               local_mof(vofcomprecon+dir)= &
                  state_mof(D_DECL(i,j,k),vofcompraw+dir)
              enddo
             else if (is_rigid(im).eq.1) then
              local_mof(vofcomprecon)=vfrac_solid_new(im) 
              do dir=1,SDIM
               local_mof(vofcomprecon+dir)=censolid_new(im,dir)
              enddo
             else
              print *,"is_rigid invalid LEVELSET_3D.F90"
              stop
             endif

             if ((local_mof(vofcomprecon).lt.-0.1).or. &
                 (local_mof(vofcomprecon).gt.1.1)) then
              print *,"local_mof(vofcomprecon) invalid"
              print *,"local_mof(vofcomprecon)=",local_mof(vofcomprecon)
              stop
             endif

             orderflag=zero
             local_mof(vofcomprecon+SDIM+1)=orderflag

             do dir=SDIM+3,ngeom_recon
              local_mof(vofcomprecon+dir-1)=zero
             enddo

            enddo  ! im=1..num_materials

            ! sum of F_fluid=1
            ! sum of F_rigid<=1
            tessellate=0
            call make_vfrac_sum_ok_base( &
              cmofsten, &
              xsten,nhalf,nhalf_box, &
              bfact,dx, &
              tessellate, &
              local_mof, &
              SDIM)

            do im=1,num_materials
             vofcomprecon=(im-1)*ngeom_recon+1
             vof_super(im)=local_mof(vofcomprecon)
            enddo

            mof_verbose=0

            call multimaterial_MOF( &
             bfact,dx,xsten,nhalf, &
             mof_verbose, &
             use_ls_data, &
             LS_stencil, &
             geom_xtetlist(1,1,1,tid+1), &
             geom_xtetlist_old(1,1,1,tid+1), &
             nmax, &
             nmax, &
             local_mof, &
             vof_super, &
             multi_centroidA, &
             continuous_mof_parm, & !continuous_mof_parm=0
             cmofsten, &
             grid_index, &
             grid_level, &
             SDIM)
     
            tessellate_transfer=1 
            call multi_get_volume_tessellate( &
             tessellate_transfer, &
             bfact, &
             dx,xsten,nhalf, &
             local_mof, &
             geom_xtetlist(1,1,1,tid+1), &
             nmax, &
             nmax, &
             SDIM)

             ! change the solid material into the wetting fluid
             ! vfrac_solid_new(im_solid_max)
             ! censolid_new(im_solid_max,dir)
             ! vof_fluid_new=vof_fluid_old + vof_solid
             ! x_fluid_new F_fluid_new = x_fluid_old F_fluid_old + x_sol F_sol
            vofcomprecon=(center_stencil_wetting_im-1)*ngeom_recon+1

            F_fluid_new=local_mof(vofcomprecon)+vfrac_solid_new(im_solid_max)
            if (F_fluid_new.gt.zero) then
             do dir=1,SDIM
              x_fluid_new(dir)= &
               local_mof(vofcomprecon+dir)*local_mof(vofcomprecon)+ &
               censolid_new(im_solid_max,dir)*vfrac_solid_new(im_solid_max)
              x_fluid_new(dir)=x_fluid_new(dir)/F_fluid_new
              local_mof(vofcomprecon+dir)=x_fluid_new(dir)
             enddo
             local_mof(vofcomprecon)=F_fluid_new
            else
             print *,"F_fluid_new invalid"
             stop
            endif
           else if (center_stencil_wetting_im.eq.0) then
            ! do nothing
           else
            print *,"center_stencil_wetting_im invalid"
            stop
           endif

           do im=1,num_materials

            if (is_rigid(im).eq.1) then
             ! do nothing
            else if (is_rigid(im).eq.0) then

             do i1=istenlo(1),istenhi(1)
             do j1=istenlo(2),istenhi(2)
             do k1=istenlo(3),istenhi(3)
              LS_temp(D_DECL(i1,j1,k1))=LS_extend(D_DECL(i1,j1,k1),im)
             enddo 
             enddo 
             enddo 

             vofcomp=(im-1)*ngeom_recon+1

              ! if there is just one fluid material in the stencil, then
              ! "center_stencil_im_only" holds the material id
             if ((center_stencil_im_only.ge.1).and. &
                 (center_stencil_im_only.le.num_materials)) then
              if (center_stencil_im_only.eq.im) then
               mofnew(vofcomp)=one
              else 
               mofnew(vofcomp)=zero
              endif
              do dir=1,SDIM
               mofnew(vofcomp+dir)=zero
              enddo
             else if ((center_stencil_wetting_im.ge.1).and. &
                      (center_stencil_wetting_im.le.num_materials)) then 
              mofnew(vofcomp)=local_mof(vofcomp)
              do dir=1,SDIM
               mofnew(vofcomp+dir)=local_mof(vofcomp+dir)
              enddo
             else if ((center_stencil_im_only.eq.0).and. &
                      (center_stencil_wetting_im.eq.0)) then
              call getvolume(bfact,dx,xsten,nhalf, &
               LS_temp,mofnew(vofcomp),LSfacearea, &
               LScentroid,VOFTOL,SDIM)

              do dir=1,SDIM
               mofnew(vofcomp+dir)=LScentroid(dir)-cencell(dir)
              enddo
             else
              print *,"center_stencil data invalid"
              stop
             endif

              ! this is the extrapolated level set function
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
             print *,"is_rigid invalid LEVELSET_3D.F90"
             stop
            endif

           enddo ! im=1..num_materials

          else if ((LS_solid_new(im_solid_max).le.zero).and. &
                   (sum_vfrac_solid_new.le.half)) then
           ! do nothing
          else
           print *,"LS_solid_new or sum_vfrac_solid_new invalid"
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
         call make_vfrac_sum_ok_base( &
           cmofsten, &
           xsten,nhalf,nhalf_box, &
           bfact,dx, &
           tessellate,mofnew,SDIM)

         do im=1,num_materials*(1+SDIM)
          lsnew(D_DECL(i,j,k),im)=ls_hold(im)
         enddo
 
         ! above: prescribe new solid location
         ! below: just normalize the volume fractions.
        else if (renormalize_only.eq.1) then

         if (1.eq.0) then
          print *,"i,j,k ",i,j,k
         endif

         call make_vfrac_sum_ok_base( &
           cmofsten, &
           xsten,nhalf,nhalf_box, &
           bfact,dx, &
           tessellate,mofnew,SDIM)
        else
         print *,"renormalize only invalid"
         stop
        endif

        if (ngeom_raw.ne.SDIM+1) then
         print *,"ngeom_raw invalid"
         stop
        endif

        do im=1,num_materials
         vofcomp=(im-1)*ngeom_recon+1
         vofcompraw=(im-1)*ngeom_raw+1
         vofnew(D_DECL(i,j,k),vofcompraw)=mofnew(vofcomp)
         do dir=1,SDIM
          vofnew(D_DECL(i,j,k),vofcompraw+dir)=mofnew(vofcomp+dir)
         enddo
        enddo ! im=1..num_materials

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
      end subroutine fort_renormalize_prescribe


      subroutine fort_purgeflotsam( &
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
       xlo,dx) &
      bind(c,name='fort_purgeflotsam')

      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE


      INTEGER_T, INTENT(in) :: level,finest_level

      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(in) :: xlo(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(inout) :: delta_mass(num_materials)
      INTEGER_T, INTENT(in) :: DIMDEC(maskcov)
      INTEGER_T, INTENT(in) :: DIMDEC(vofnew)
      INTEGER_T, INTENT(in) :: DIMDEC(LS)
      INTEGER_T, INTENT(in) :: truncate_volume_fractions(num_materials)
      REAL_T, INTENT(in),target :: maskcov(DIMV(maskcov))
      REAL_T, pointer :: maskcov_ptr(D_DECL(:,:,:))
      REAL_T, INTENT(inout),target :: &
        vofnew(DIMV(vofnew),num_materials*ngeom_raw)
      REAL_T, pointer :: vofnew_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in),target ::  LS(DIMV(LS),num_materials)
      REAL_T, pointer :: LS_ptr(D_DECL(:,:,:),:)
      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, INTENT(in) :: bfact
      REAL_T, INTENT(in) :: truncate_thickness

      INTEGER_T i,j,k,dir
      INTEGER_T im,im2,imcrit
      INTEGER_T vofcomprecon,vofcompraw

      REAL_T mofdata(num_materials*ngeom_recon)
      REAL_T volmat(num_materials)
      REAL_T lspoint(num_materials)
      INTEGER_T sorted_list(num_materials)
      REAL_T dxmax,dxmaxLS,LSbandsize,restore_sum
      INTEGER_T, parameter :: nhalf=3
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T volgrid
      REAL_T cengrid(SDIM)
      INTEGER_T mask_test
      INTEGER_T FSI_exclude
      INTEGER_T tessellate
      INTEGER_T, parameter :: nhalf_box=1
      INTEGER_T cmofsten(D_DECL(-1:1,-1:1,-1:1))

      tessellate=0

      if (bfact.lt.1) then
       print *,"bfact invalid103"
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

      maskcov_ptr=>maskcov
      vofnew_ptr=>vofnew
      LS_ptr=>LS
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,vofnew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)

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
        do dir=1,num_materials*ngeom_recon
         mofdata(dir)=zero 
        enddo
        do im=1,num_materials
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
        call sort_volume_fraction(volmat,FSI_exclude,sorted_list)
        imcrit=sorted_list(1)
        if (is_rigid(imcrit).eq.0) then
         ! do nothing
        else
         print *,"is_rigid(imcrit) invalid"
         stop
        endif
        do im=1,num_materials
         if (is_rigid(im).eq.0) then

          if (lspoint(im).gt.LSbandsize) then
           vofcomprecon=(im-1)*ngeom_recon+1
           mofdata(vofcomprecon)=one
           do dir=1,SDIM
            mofdata(vofcomprecon+dir)=zero
           enddo

           restore_sum=zero

           do im2=1,num_materials
            if (is_rigid(im2).eq.0) then

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

             else if (im2.eq.im) then
              ! do nothing
             else
              print *,"im2 and im mismatch"
              stop
             endif 

            else if (is_rigid(im2).eq.1) then
             ! do nothing
            else
             print *,"is_rigid(im2) invalid"
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

         else if (is_rigid(im).eq.1) then
          ! do nothing
         else
          print *,"is_rigid(im) invalid"
          stop
         endif
        enddo ! im=1...num_materials

         ! sum F_fluid=1  sum F_solid <=1
        call make_vfrac_sum_ok_base( &
          cmofsten, &
          xsten,nhalf,nhalf_box, &
          bfact,dx, &
          tessellate,mofdata,SDIM)

        do im=1,num_materials
         vofcompraw=(im-1)*ngeom_raw+1
         vofcomprecon=(im-1)*ngeom_recon+1
         vofnew(D_DECL(i,j,k),vofcompraw)=mofdata(vofcomprecon)
         do dir=1,SDIM
          vofnew(D_DECL(i,j,k),vofcompraw+dir)=mofdata(vofcomprecon+dir)
         enddo

         delta_mass(im)=delta_mass(im)+ &
           volgrid*(mofdata(vofcomprecon)-volmat(im))
        enddo ! im=1..num_materials

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
      end subroutine fort_purgeflotsam

      subroutine fort_initrecalesce( &
       recalesce_material_in, &
       recalesce_state_old_in, &
       recalesce_num_state_in) &
      bind(c,name='fort_initrecalesce')

      use probcommon_module
      use probf90_module

      IMPLICIT NONE


      INTEGER_T, INTENT(in) :: recalesce_num_state_in
      INTEGER_T i
      INTEGER_T, INTENT(in) :: recalesce_material_in(num_materials)
      REAL_T, INTENT(in) :: &
           recalesce_state_old_in(recalesce_num_state*num_materials)

      if (num_materials.gt.100) then
       print *,"too many materials"
       stop
      endif
      if (recalesce_num_state.ne.recalesce_num_state_in) then
       print *,"recalesce_num_state_in invalid"
       stop
      endif

      do i=1,num_materials
       recalesce_material(i)=recalesce_material_in(i)
      enddo
      do i=1,recalesce_num_state*num_materials
       recalesce_state_old(i)=recalesce_state_old_in(i)
      enddo

      return
      end subroutine fort_initrecalesce

      end module levelset_module


      module FSI_PC_LS_module

       use local_amrex_fort_module, only : amrex_real,amrex_particle_real
       use iso_c_binding, only: c_int

       implicit none

       type, bind(C) :: particle_t
         real(amrex_particle_real) :: pos(SDIM)
         ! (insert time) is extra. 
         real(amrex_particle_real) :: extra_state(N_EXTRA_REAL)
         integer(c_int) :: id
         integer(c_int) :: cpu
         ! (material_id) is extra.
         integer(c_int) :: extra_int(N_EXTRA_INT)
       end type particle_t

       type accum_parm_type_count
        INTEGER_T :: fablo(SDIM)
        INTEGER_T :: fabhi(SDIM)
        INTEGER_T :: tilelo(SDIM)
        INTEGER_T :: tilehi(SDIM)
        INTEGER_T :: append_flag
        INTEGER_T :: bfact
        INTEGER_T :: level
        INTEGER_T :: finest_level
        REAL_T :: dx(SDIM)
        REAL_T :: xlo(SDIM)
        INTEGER_T :: Npart
!        type(particle_t), pointer, dimension(:) :: particlesptr
        INTEGER_T :: N_real_comp
!        REAL_T, pointer, dimension(:) :: real_compALLptr
        INTEGER_T :: nsubdivide
!        REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: TENSORptr
!        REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LEVELSETptr
        !cell_particle_count(i,j,k,1)=number particles in the cell.
        !cell_particle_count(i,j,k,2)=particle id of first particle in list
!        INTEGER_T, pointer, dimension(D_DECL(:,:,:),:) :: &
!           cell_particle_count
       end type accum_parm_type_count

       type grid_parm_type
        INTEGER_T :: fablo(SDIM)
        INTEGER_T :: fabhi(SDIM)
        INTEGER_T :: tilelo(SDIM)
        INTEGER_T :: tilehi(SDIM)
        INTEGER_T :: bfact
        INTEGER_T :: level
        INTEGER_T :: finest_level
        REAL_T :: dx(SDIM)
        REAL_T :: xlo(SDIM)
!        REAL_T, pointer, dimension(D_DECL(:,:,:)) :: umacptr
!        REAL_T, pointer, dimension(D_DECL(:,:,:)) :: vmacptr
!        REAL_T, pointer, dimension(D_DECL(:,:,:)) :: wmacptr
        INTEGER_T, dimension(SDIM,2,SDIM) :: velbc
        INTEGER_T, dimension(SDIM,2) :: dombc
        INTEGER_T :: domlo(SDIM)
        INTEGER_T :: domhi(SDIM)
        REAL_T :: problo(SDIM)
        REAL_T :: probhi(SDIM)
       end type grid_parm_type

      contains

      subroutine count_particles( &
       accum_PARM, &
       particlesptr, &
       cell_particle_count, &
       particle_link_data, &
       Np)

      use global_utility_module
      use mass_transfer_module

      type(accum_parm_type_count), INTENT(in) :: accum_PARM
      type(particle_t), INTENT(in), pointer, dimension(:) :: particlesptr
      INTEGER_T, value, INTENT(in) :: Np 
       ! child link 1 (particle index), parent link 1 (i,j,k index),
       ! child link 2 (particle index), parent link 2 (i,j,k index), ...
      INTEGER_T, INTENT(inout) :: particle_link_data(Np*(1+SDIM))
      INTEGER_T, INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: &
           cell_particle_count

      INTEGER_T :: interior_ID
      INTEGER_T :: dir
      REAL_T, target :: xpart(SDIM)
      INTEGER_T cell_index(SDIM)
      INTEGER_T interior_ok
      INTEGER_T i,j,k
      INTEGER_T :: local_ngrow
      INTEGER_T :: ok_to_add_link
      INTEGER_T :: previous_link
      INTEGER_T :: ibase
      INTEGER_T :: ibase_new
      INTEGER_T :: i_parent,j_parent,k_parent

      local_ngrow=1

      do interior_ID=1,accum_PARM%Npart

       do dir=1,SDIM
        xpart(dir)=particlesptr(interior_ID)%pos(dir)
       enddo
       call containing_cell(accum_PARM%bfact, &
         accum_PARM%dx, &
         accum_PARM%xlo, &
         accum_PARM%fablo, &
         xpart, &
         cell_index)

       interior_ok=1
       do dir=1,SDIM
        if ((cell_index(dir).lt.accum_PARM%tilelo(dir)).or. &
            (cell_index(dir).gt.accum_PARM%tilehi(dir))) then
         interior_ok=0
        endif
       enddo
       if (interior_ok.eq.1) then

        i=cell_index(1)
        j=cell_index(2)
        k=cell_index(SDIM)

        ok_to_add_link=1

        previous_link=cell_particle_count(D_DECL(i,j,k),2)
        if (previous_link.eq.0) then
         ! do nothing; no particles attached to cell (i,j,k)
        else if ((previous_link.ge.1).and. &
                 (previous_link.le.Np)) then
         ibase=(previous_link-1)*(SDIM+1)
         i_parent=particle_link_data(ibase+2) 
         j_parent=particle_link_data(ibase+3) 
         k_parent=particle_link_data(ibase+1+SDIM) 
         if ((i.eq.i_parent).and. &
             (j.eq.j_parent)) then
          ! do nothing
         else
          ok_to_add_link=0
         endif
         if (SDIM.eq.3) then
          if (k.eq.k_parent) then
           ! do nothing
          else
           ok_to_add_link=0
          endif
         endif
        else
         print *,"previous_link invalid"
         stop
        endif

        if (ok_to_add_link.eq.1) then

         if (previous_link.eq.interior_ID) then
          print *,"links should be unique"
          stop
         endif

         cell_particle_count(D_DECL(i,j,k),1)= &
           cell_particle_count(D_DECL(i,j,k),1)+1

         cell_particle_count(D_DECL(i,j,k),2)=interior_ID

         ibase_new=(interior_ID-1)*(SDIM+1)
         particle_link_data(ibase_new+1)=previous_link
         particle_link_data(ibase_new+2)=i 
         particle_link_data(ibase_new+3)=j
         if (SDIM.eq.3) then 
          particle_link_data(ibase_new+SDIM+1)=k
         endif

        else if (ok_to_add_link.eq.0) then
         print *,"ok_to_add_link==0"
         print *,"the parent of a particle added to (i,j,k) should"
         print *,"always be (i,j,k)"
         print *,"i,j,k ",i,j,k
         print *,"i_parent,j_parent,k_parent ", &
           i_parent,j_parent,k_parent
         print *,"previous_link,interior_ID ", &
            previous_link,interior_ID 
         stop
        else
         print *,"ok_to_add_link invalid"
         stop
        endif

       else if (interior_ok.eq.0) then
        ! do nothing
       else
        print *,"interior_ok invalid"
        stop
       endif

      enddo ! do interior_ID=1,accum_PARM%Npart

      return
      end subroutine count_particles


      subroutine containing_sub_box( &
         accum_PARM, &
         xpart, &
         i,j,k, &
         isub,jsub,ksub, &
         sub_found)
      use global_utility_module

      IMPLICIT NONE

      type(accum_parm_type_count), INTENT(in) :: accum_PARM
      REAL_T, INTENT(in) :: xpart(SDIM)
      INTEGER_T, INTENT(in) :: i,j,k
      INTEGER_T, INTENT(out) :: isub,jsub,ksub
      INTEGER_T, INTENT(out) :: sub_found
      INTEGER_T, parameter :: nhalf=3
      INTEGER_T :: dir
      REAL_T :: xsten(-nhalf:nhalf,SDIM)
      INTEGER_T :: isub_local(3)
      REAL_T :: dx_sub

      call gridsten_level(xsten,i,j,k,accum_PARM%level,nhalf)
      sub_found=1
      isub_local(3)=0
       ! x=xlo + (isub+1/2)*dx
       ! isub=(x-xlo)/dx -1/2
      do dir=1,SDIM
       dx_sub=(xsten(1,dir)-xsten(-1,dir))/accum_PARM%nsubdivide
       if (dx_sub.gt.zero) then
        if (abs(xpart(dir)-xsten(-1,dir)).le.1.0D-8*dx_sub) then
         isub_local(dir)=0
        else if (abs(xpart(dir)-xsten(1,dir)).le.1.0D-8*dx_sub) then
         isub_local(dir)=accum_PARM%nsubdivide-1
        else if ((xpart(dir).ge.xsten(-1,dir)).and. &
                 (xpart(dir).le.xsten(1,dir))) then
         isub_local(dir)=NINT((xpart(dir)-xsten(-1,dir))/dx_sub-half)
        else
         isub_local(dir)=0
         sub_found=0
        endif
       else
        print *,"dx_sub invalid"
        stop
       endif
      enddo ! dir=1..sdim

      isub=isub_local(1)
      jsub=isub_local(2)
      ksub=isub_local(3)

      end subroutine containing_sub_box

      subroutine sub_box_cell_center( &
         accum_PARM, &
         i,j,k, &
         isub,jsub,ksub, &
         xsub)
      use global_utility_module

      IMPLICIT NONE

      type(accum_parm_type_count), INTENT(in) :: accum_PARM
      REAL_T, INTENT(out) :: xsub(SDIM)
      INTEGER_T, INTENT(in) :: i,j,k
      INTEGER_T, INTENT(in) :: isub,jsub,ksub
      INTEGER_T, parameter :: nhalf=3
      INTEGER_T :: dir
      REAL_T :: xsten(-nhalf:nhalf,SDIM)
      REAL_T :: dx_sub
      INTEGER_T isub_local(3)

      isub_local(1)=isub
      isub_local(2)=jsub
      isub_local(3)=ksub

      call gridsten_level(xsten,i,j,k,accum_PARM%level,nhalf)
       ! x=xlo + (isub+1/2)*dx
      do dir=1,SDIM
       dx_sub=(xsten(1,dir)-xsten(-1,dir))/accum_PARM%nsubdivide
       if (dx_sub.gt.zero) then
        xsub(dir)=xsten(-1,dir)+(isub_local(dir)+half)*dx_sub
       else
        print *,"dx_sub invalid"
        stop
       endif
      enddo ! dir=1..sdim

      end subroutine sub_box_cell_center

       ! partition of unity interpolation:
       ! Q(xtarget)=( sum_{particles} Q_{p} w(xtarget,xp) +
       !              sum_{grid} Q_{grid} w(xtarget,xgrid) ) /
       !            ( sum_{particles} w(xtarget,xp) +
       !              sum_{grid} w(xtarget,xgrid) )
      subroutine interp_eul_lag_dist( &
         im_elastic_map, &
         accum_PARM, &
         particlesptr, &
         real_compALLptr, &
         TENSORptr, &
         LEVELSETptr, &
         cell_particle_count, &
         i,j,k, &
         xtarget, &  ! where to add the new particle
         particle_link_data, &
         Np, &
         Q_interp, &
         LS_interp, &
         X0_interp)
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

       ! 0<=im_elastic_map<num_materials
      INTEGER_T, INTENT(in) :: im_elastic_map(num_materials_viscoelastic)

      type(accum_parm_type_count), INTENT(in) :: accum_PARM
      type(particle_t), INTENT(in), pointer, dimension(:) :: particlesptr
      REAL_T, INTENT(in), pointer, dimension(:) :: real_compALLptr
      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: TENSORptr
      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: LEVELSETptr
      INTEGER_T, INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: &
         cell_particle_count
      INTEGER_T, INTENT(in) :: i,j,k
      REAL_T, target, INTENT(in) :: xtarget(SDIM)
      INTEGER_T, value, INTENT(in) :: Np
      INTEGER_T, INTENT(in) :: particle_link_data(Np*(1+SDIM))
      REAL_T, INTENT(out) :: Q_interp(NUM_CELL_ELASTIC)
      REAL_T, INTENT(out) :: LS_interp(num_materials)
      REAL_T, INTENT(out) :: X0_interp(SDIM)

      INTEGER_T, parameter :: nhalf=3
      INTEGER_T :: dir
      REAL_T :: xsten(-nhalf:nhalf,SDIM)
      REAL_T A_VEL,b_VEL(NUM_CELL_ELASTIC)
      REAL_T A_X0,b_X0(SDIM)
      INTEGER_T :: current_link
      REAL_T, target :: xpart(SDIM)
      REAL_T :: Qpart(NUM_CELL_ELASTIC)
      REAL_T :: X0part(SDIM)
      REAL_T :: LS_interp_local(num_materials)
      INTEGER_T :: ibase

      REAL_T w_p

      REAL_T :: LSlocal(num_materials)
      INTEGER_T :: im
      INTEGER_T :: im_primary
      INTEGER_T :: im_primary_part
      INTEGER_T :: im_particle_direct
      INTEGER_T :: im_primary_sub
      INTEGER_T :: ipart
      INTEGER_T :: im_map
      INTEGER_T :: dir_tensor

      type(interp_from_grid_parm_type) :: data_in 
      type(interp_from_grid_out_parm_type) :: data_out
      type(interp_from_grid_out_parm_type) :: data_out_LS

      REAL_T, target, dimension(NUM_CELL_ELASTIC) :: data_interp_local
      REAL_T, target, dimension(num_materials) :: data_interp_local_LS

      INTEGER_T :: test_count,test_cell_particle_count

      INTEGER_T :: SoA_comp

      if (1.eq.0) then
       print *,"num_materials_viscoelastic ",num_materials_viscoelastic
       print *,"i,j,k ",i,j,k
       print *,"xtarget ",xtarget(1),xtarget(2),xtarget(SDIM)
       print *,"Np= ",Np
       print *,"num_materials=",num_materials
       print *,"im_elastic_map(1) ",im_elastic_map(1)
       print *,"accum_PARM%fablo(1) ",accum_PARM%fablo(1)
       print *,"accum_PARM%fablo(2) ",accum_PARM%fablo(2)
       print *,"accum_PARM%fablo(SDIM) ",accum_PARM%fablo(SDIM)
       print *,"accum_PARM%tilelo(1) ",accum_PARM%tilelo(1)
       print *,"accum_PARM%tilelo(2) ",accum_PARM%tilelo(2)
       print *,"accum_PARM%tilelo(SDIM) ",accum_PARM%tilelo(SDIM)
       print *,"accum_PARM%Npart ",accum_PARM%Npart
       print *,"accum_PARM%N_real_comp ",accum_PARM%N_real_comp
       print *,"accum_PARM%nsubdivide ",accum_PARM%nsubdivide
       print *,"LBOUND(TENSORptr) ",LBOUND(TENSORptr)
       print *,"UBOUND(TENSORptr) ",UBOUND(TENSORptr)
       print *,"LBOUND(LEVELSETptr) ",LBOUND(LEVELSETptr)
       print *,"UBOUND(LEVELSETptr) ",UBOUND(LEVELSETptr)
      endif

      call gridsten_level(xsten,i,j,k,accum_PARM%level,nhalf)

      call checkbound_array(accum_PARM%fablo,accum_PARM%fabhi, &
         TENSORptr,1,-1)
      call checkbound_array(accum_PARM%fablo,accum_PARM%fabhi, &
         LEVELSETptr,1,-1)
      call checkbound_array_INTEGER(accum_PARM%tilelo,accum_PARM%tilehi, &
         cell_particle_count,0,-1)

      do im=1,num_materials
       LSlocal(im)=LEVELSETptr(D_DECL(i,j,k),im)
      enddo
      call get_primary_material(LSlocal,im_primary)

      data_out%data_interp=>data_interp_local
      data_out_LS%data_interp=>data_interp_local_LS

      data_in%scomp=1 
      data_in%ncomp=NUM_CELL_ELASTIC
      data_in%level=accum_PARM%level
      data_in%finest_level=accum_PARM%finest_level
      data_in%bfact=accum_PARM%bfact

      do dir=1,SDIM
       data_in%dx(dir)=accum_PARM%dx(dir)
       data_in%xlo(dir)=accum_PARM%xlo(dir)
       data_in%fablo(dir)=accum_PARM%fablo(dir)
       data_in%fabhi(dir)=accum_PARM%fabhi(dir)
      enddo

      A_VEL=zero
      A_X0=zero
      do dir=1,NUM_CELL_ELASTIC
       b_VEL(dir)=zero
      enddo
      do dir=1,SDIM
       b_X0(dir)=zero
      enddo

      test_cell_particle_count=cell_particle_count(D_DECL(i,j,k),1)

      test_count=0

      current_link=cell_particle_count(D_DECL(i,j,k),2)

      do while ((current_link.ge.1).and.(current_link.le.Np))

       do dir=1,SDIM
        xpart(dir)=particlesptr(current_link)%pos(dir)
       enddo 

       do dir=1,NUM_CELL_ELASTIC
        SoA_comp=(dir-1)*Np+current_link
        if ((SoA_comp.ge.1).and. &
            (SoA_comp.le.accum_PARM%N_real_comp)) then
         Qpart(dir)=real_compALLptr(SoA_comp)
        else
         print *,"SoA_comp invalid"
         stop
        endif
       enddo !dir=1,NUM_CELL_ELASTIC

       do dir=1,SDIM
        X0part(dir)= &
           particlesptr(current_link)%extra_state(N_EXTRA_REAL_X0+dir)
       enddo !dir=1,SDIM

       im_particle_direct= &
         particlesptr(current_link)%extra_int(N_EXTRA_INT_MATERIAL_ID+1)

       if ((im_particle_direct.ge.1).and. &
           (im_particle_direct.le.num_materials)) then
        ! do nothing
       else
        print *,"im_particle_direct invalid"
        stop
       endif

       call partition_unity_weight(xpart,xtarget,accum_PARM%dx,w_p)

       data_in%xtarget=xpart

       if (num_materials.gt.0) then
        data_in%scomp=1 
        data_in%ncomp=num_materials
        call interp_from_grid_util(data_in,LEVELSETptr,data_out_LS)
        do dir=1,num_materials
         LS_interp_local(dir)=data_out_LS%data_interp(dir)
        enddo
        call get_primary_material(LS_interp_local,im_primary_part)
       else
        print *,"num_materials invalid"
        stop
       endif

       if (accum_PARM%append_flag.eq.0) then
        print *,"there should not be any particles if append_flag==0"
        stop
       else if (accum_PARM%append_flag.eq.1) then
        ! do nothing
       else 
        print *,"accum_PARM%append_flag invalid" 
        stop
       endif

       if (w_p.gt.zero) then

        A_X0=A_X0+w_p
        do dir=1,SDIM
         b_X0(dir)=b_X0(dir)+w_p*X0part(dir)
        enddo

        if (im_primary_part.eq.im_primary) then
         if (im_primary_part.eq.im_particle_direct) then
          A_VEL=A_VEL+w_p
          do dir=1,NUM_CELL_ELASTIC
           b_VEL(dir)=b_VEL(dir)+w_p*Qpart(dir)
          enddo
         else if (im_primary_part.ne.im_particle_direct) then
          ! do nothing
         else
          print *,"im_primary_part or im_particle_direct bust"
          stop
         endif
        else if (im_primary_part.ne.im_primary) then
         ! do nothing
        else
         print *,"im_primary_part or im_primary bust"
         stop
        endif
       else
        print *,"w_p invalid"
        stop
       endif

       ibase=(current_link-1)*(1+SDIM)
       current_link=particle_link_data(ibase+1)

       test_count=test_count+1
      enddo ! while (current_link.ge.1).and.(current_link<=Np)

      if (current_link.eq.0) then
       ! do nothing
      else
       print *,"current_link invalid"
       stop
      endif
      if (test_count.eq.test_cell_particle_count) then
       ! do nothing
      else
       print *,"test_cell_particle_count invalid"
       stop
      endif

       ! xtarget might not coincide with an Eulerian grid cell.
      data_in%xtarget=xtarget

       ! bilinear interpolation
      if (NUM_CELL_ELASTIC.gt.0) then
       data_in%scomp=1 
       data_in%ncomp=NUM_CELL_ELASTIC
       call interp_from_grid_util(data_in,TENSORptr,data_out)
       do dir=1,NUM_CELL_ELASTIC
        Q_interp(dir)=data_out%data_interp(dir)
       enddo
      else if (NUM_CELL_ELASTIC.eq.0) then
       ! do nothing
      else
       print *,"NUM_CELL_ELASTIC invalid"
       stop
      endif

      if (num_materials.gt.0) then
       data_in%scomp=1 
       data_in%ncomp=num_materials
       call interp_from_grid_util(data_in,LEVELSETptr,data_out_LS)
       do dir=1,num_materials
        LS_interp(dir)=data_out_LS%data_interp(dir)
       enddo
       call get_primary_material(LS_interp,im_primary_sub)
      else
       print *,"num_materials invalid"
       stop
      endif

      if (accum_PARM%append_flag.eq.0) then
       if (A_VEL.eq.zero) then
        ! do nothing
       else
        print *,"expecting A_VEL==0 if append_flag==0"
        stop
       endif
       if (A_X0.eq.zero) then
        ! do nothing
       else
        print *,"expecting A_X0==0 if append_flag==0"
        stop
       endif
      else if (accum_PARM%append_flag.eq.1) then

       if (A_VEL.gt.zero) then
        call partition_unity_weight(xtarget,xtarget,accum_PARM%dx,w_p)

        do ipart=1,num_materials_viscoelastic
         im_map=im_elastic_map(ipart)+1
         if ((im_map.ge.1).and.(im_map.le.num_materials)) then
          if ((im_map.eq.im_primary).and. &
              (im_map.eq.im_primary_sub)) then
           do dir=1,ENUM_NUM_TENSOR_TYPE
            dir_tensor=(ipart-1)*ENUM_NUM_TENSOR_TYPE+dir
            Q_interp(dir_tensor)= &
                  (b_VEL(dir_tensor)+w_p*Q_interp(dir_tensor))/ &
                  (A_VEL+w_p)
           enddo
          else if ((im_map.ne.im_primary).or. &
                   (im_map.ne.im_primary_sub)) then
           ! do nothing
          else
           print *,"im_map, im_primary, or im_primary_sub bust"
           stop
          endif
         else
          print *,"im_map invalid"
          stop
         endif
        enddo !ipart=1,num_materials_viscoelastic

       else if (A_VEL.eq.zero) then
        ! do nothing
       else
        print *,"A_VEL invalid"
        stop
       endif

       if (A_X0.gt.zero) then

        do dir=1,SDIM
         X0_interp(dir)=b_X0(dir)/A_X0
        enddo

       else if (A_X0.eq.zero) then

        do dir=1,SDIM
         X0_interp(dir)=xtarget(dir)
        enddo

       else
        print *,"A_X0 invalid"
        stop
       endif

      else 
       print *,"accum_PARM%append_flag invalid" 
       stop
      endif

      return
      end subroutine interp_eul_lag_dist

       ! called from NavierStokes.cpp:
       !  NavierStokes::init_particle_container
      subroutine fort_init_particle_container( &
        tid, &
        single_particle_size, &
        isweep, &
        append_flag, &
        particle_nsubdivide, &
        particle_max_per_nsubdivide, &
        particle_min_per_nsubdivide, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        level, &
        finest_level, &
        cur_time_slab, &
        xlo,dx, &
        particles, & ! a list of particles in the elastic structure
        Np, & !  Np = number of particles
        real_compALL, &
        N_real_comp, & ! pass by value
        new_particles, & ! size is "new_Pdata_size"
        new_Pdata_size, &
        Np_append, & ! number of particles to add
        particle_link_data, &
        particle_delete_flag, & ! 1=> delete
        cell_particle_count, &
        DIMS(cell_particle_count), &
        tensorfab,DIMS(tensorfab), &
        lsfab,DIMS(lsfab), &
        mfiner,DIMS(mfiner), &
        viscoelastic_model, &
        im_elastic_map, &
        polymer_factor) &
      bind(c,name='fort_init_particle_container')

      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: tid
      INTEGER_T, INTENT(in) :: single_particle_size
      INTEGER_T, INTENT(in) :: isweep
      INTEGER_T, INTENT(in) :: append_flag
      INTEGER_T, INTENT(in) :: level,finest_level

      INTEGER_T, INTENT(in), target :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, INTENT(in) :: bfact
      INTEGER_T, INTENT(in) :: particle_nsubdivide
      INTEGER_T, INTENT(in) :: particle_max_per_nsubdivide
      INTEGER_T, INTENT(in) :: particle_min_per_nsubdivide
      REAL_T, INTENT(in)    :: cur_time_slab
      REAL_T, INTENT(in), target :: xlo(SDIM),dx(SDIM)
      INTEGER_T, value, INTENT(in) :: Np ! pass by value
      type(particle_t), INTENT(in), target :: particles(Np)
      type(particle_t), pointer :: particlesptr(:)
      INTEGER_T, value, INTENT(in) :: N_real_comp ! pass by value
      REAL_T, INTENT(in), target :: real_compALL(N_real_comp)
      REAL_T, pointer :: real_compALLptr(:)
      INTEGER_T, INTENT(inout) :: new_Pdata_size
      REAL_T, INTENT(out) :: new_particles(new_Pdata_size)
      INTEGER_T, INTENT(inout) :: Np_append

       ! child link 1, parent link 1,
       ! child link 2, parent link 2, ...
      INTEGER_T, INTENT(inout) :: particle_link_data(Np*(1+SDIM))
      INTEGER_T, INTENT(inout) :: particle_delete_flag(Np) ! 1=>delete

      INTEGER_T, INTENT(in) :: DIMDEC(cell_particle_count)
      INTEGER_T, INTENT(in) :: DIMDEC(tensorfab)
      INTEGER_T, INTENT(in) :: DIMDEC(lsfab)
      INTEGER_T, INTENT(in) :: DIMDEC(mfiner)
   
       ! first component: number of particles in the cell
       ! second component: link to the local particle container: 1..Np 
      INTEGER_T, INTENT(inout), target :: cell_particle_count( &
              DIMV(cell_particle_count), &
              2) 
      INTEGER_T, pointer, &
        dimension(D_DECL(:,:,:),:) :: cell_particle_count_ptr

      REAL_T, INTENT(in), target :: &
         tensorfab(DIMV(tensorfab),NUM_CELL_ELASTIC) 
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: tensorfab_ptr
      REAL_T, INTENT(in), target :: lsfab(DIMV(lsfab),num_materials) 
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: lsfab_ptr
      REAL_T, INTENT(in), target :: mfiner(DIMV(mfiner)) 
      REAL_T, pointer, dimension(D_DECL(:,:,:)) :: mfiner_ptr

      INTEGER_T, INTENT(in) :: viscoelastic_model(num_materials)
       ! 0<=im_elastic_map<num_materials
      INTEGER_T, INTENT(in) :: im_elastic_map(num_materials_viscoelastic)
      REAL_T, INTENT(in) :: polymer_factor(num_materials)

      type(accum_parm_type_count) :: accum_PARM
   
      INTEGER_T :: i,j,k
      INTEGER_T :: isub,jsub,ksub
      INTEGER_T :: dir
      INTEGER_T :: ibase
      INTEGER_T growlo(3) 
      INTEGER_T growhi(3) 
      INTEGER_T sublo(3) 
      INTEGER_T subhi(3) 
      INTEGER_T, allocatable :: sub_counter(:,:,:)
      INTEGER_T cell_count_check
      INTEGER_T cell_count_hold
      INTEGER_T current_link
      INTEGER_T local_count
      INTEGER_T sub_found
      INTEGER_T Np_append_test
      REAL_T, target :: xpart(SDIM)
      REAL_T :: xsub(SDIM)
      REAL_T :: tensor_sub(NUM_CELL_ELASTIC)
      REAL_T :: X0_sub(SDIM)
      REAL_T :: LS_sub(num_materials)
      INTEGER_T :: im_primary_sub
      INTEGER_T :: im_particle
      INTEGER_T, allocatable, dimension(:,:) :: sub_particle_data
      INTEGER_T, allocatable, dimension(:) :: sort_data_id
      REAL_T, allocatable, dimension(:) :: sort_data_time
      INTEGER_T sub_iter
      INTEGER_T cell_iter
      INTEGER_T isub_test
      INTEGER_T jsub_test
      INTEGER_T ksub_test
      INTEGER_T bubble_change
      INTEGER_T bubble_iter
      INTEGER_T ibubble
      INTEGER_T temp_id
      INTEGER_T local_mask
      REAL_T temp_time

      INTEGER_T ipart
      INTEGER_T im_map
      INTEGER_T ii,jj
      REAL_T Q(3,3)

      type(interp_from_grid_parm_type) :: data_in 
      type(interp_from_grid_out_parm_type) :: data_out_LS
      REAL_T, target, dimension(num_materials) :: data_interp_local_LS
      REAL_T time_in_future

      cell_particle_count_ptr=>cell_particle_count
      mfiner_ptr=>mfiner
      tensorfab_ptr=>tensorfab
      lsfab_ptr=>lsfab

      call checkbound_array(fablo,fabhi,tensorfab_ptr,1,-1)
      call checkbound_array(fablo,fabhi,lsfab_ptr,1,-1)

      call checkbound_array1(fablo,fabhi,mfiner_ptr,1,-1)

      call checkbound_array_INTEGER(tilelo,tilehi, &
              cell_particle_count_ptr,0,-1)

      if (cur_time_slab.ge.zero) then
       ! do nothing
      else
       print *,"cur_time_slab invalid"
       stop
      endif
      time_in_future=cur_time_slab+max(cur_time_slab,1.0D+3)

      if (NUM_CELL_ELASTIC.eq. &
          num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) then
       ! do nothing
      else
       print *,"NUM_CELL_ELASTIC invalid"
       stop
      endif

      if (Np*NUM_CELL_ELASTIC.eq.N_real_comp) then
       ! do nothing
      else
       print *,"N_real_comp invalid"
       stop
      endif

      if ((num_materials_viscoelastic.ge.0).and. &
          (num_materials_viscoelastic.le.num_materials)) then
       ! do nothing
      else
       print *,"num_materials_viscoelastic invalid"
       stop
      endif

      if (N_EXTRA_REAL.eq.1+SDIM) then
       ! do nothing
      else
       print *,"N_EXTRA_REAL invalid"
       stop
      endif
      if (N_EXTRA_INT.ge.1) then
       ! do nothing
      else
       print *,"N_EXTRA_INT invalid"
       stop
      endif

      if (single_particle_size.eq. &
          SDIM+N_EXTRA_REAL+N_EXTRA_INT+NUM_CELL_ELASTIC) then
       ! do nothing
      else
       print *,"single_particle_size invalid"
       stop
      endif

      if ((new_Pdata_size/single_particle_size)* &
          single_particle_size.eq.new_Pdata_size) then
       ! do nothing
      else
       print *,"new_Pdata_size invalid"
       stop
      endif

      accum_PARM%append_flag=append_flag

      do dir=1,SDIM
       accum_PARM%fablo(dir)=fablo(dir)
       accum_PARM%fabhi(dir)=fabhi(dir)
       accum_PARM%tilelo(dir)=tilelo(dir)
       accum_PARM%tilehi(dir)=tilehi(dir)
       accum_PARM%dx(dir)=dx(dir)
       accum_PARM%xlo(dir)=xlo(dir)
      enddo
      accum_PARM%bfact=bfact
      accum_PARM%level=level
      accum_PARM%finest_level=finest_level

      accum_PARM%nsubdivide=particle_nsubdivide

      particlesptr=>particles

      accum_PARM%Npart=Np

      real_compALLptr=>real_compALL

      accum_PARM%N_real_comp=N_real_comp

      if (isweep.eq.0) then
       if (append_flag.eq.1) then
        call count_particles( &
         accum_PARM, &
         particlesptr, &
         cell_particle_count_ptr, &
         particle_link_data, &
         Np)
       else if (append_flag.eq.0) then
        ! do nothing
       else
        print *,"append_flag invalid"
        stop
       endif
      else if (isweep.eq.1) then
       ! do nothing
      else
       print *,"isweep invalid"
       stop
      endif
       ! 1. traverse by cell
       ! 2. within each cell, initialize counts for each subdivision.
       ! 3. add particles within the subdivided cell using LS and
       !    xfoot from previous particle (Lagrangian) data in the 
       !    cell and previous Eulerian data.
      call growntilebox(tilelo,tilehi,fablo,fabhi, &
        growlo,growhi,0)

      sublo(3)=0
      subhi(3)=0
      do dir=1,SDIM
       sublo(dir)=0
       subhi(dir)=particle_nsubdivide-1
      enddo
      allocate(sub_counter(sublo(1):subhi(1), &
              sublo(2):subhi(2), &
              sublo(3):subhi(3)))

      if (isweep.eq.0) then
       Np_append=0
      else if (isweep.eq.1) then
       ! do nothing
      else
       print *,"isweep invalid"
       stop
      endif
      Np_append_test=0

      data_out_LS%data_interp=>data_interp_local_LS

      data_in%scomp=1 
      data_in%ncomp=num_materials
      data_in%level=level
      data_in%finest_level=finest_level
      data_in%bfact=bfact
      data_in%dx=dx
      data_in%xlo=xlo
      data_in%fablo=fablo
      data_in%fabhi=fabhi

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_mask=NINT(mfiner(D_DECL(i,j,k)))
       if (local_mask.eq.1) then

        do isub=sublo(1),subhi(1)
        do jsub=sublo(2),subhi(2)
        do ksub=sublo(3),subhi(3)
         sub_counter(isub,jsub,ksub)=0
        enddo
        enddo
        enddo
        cell_count_check=0
        cell_count_hold=cell_particle_count(D_DECL(i,j,k),1)

         ! isub,jsub,ksub,link
        allocate(sub_particle_data(cell_count_hold,SDIM+1))

        current_link=cell_particle_count(D_DECL(i,j,k),2)
        do while (current_link.ge.1)
         do dir=1,SDIM
          xpart(dir)=particles(current_link)%pos(dir)
         enddo 
         call containing_sub_box( &
           accum_PARM, &
           xpart, &
           i,j,k, &
           isub,jsub,ksub, &
           sub_found)
         if (sub_found.eq.1) then
          sub_counter(isub,jsub,ksub)=sub_counter(isub,jsub,ksub)+1
          cell_count_check=cell_count_check+1
          sub_particle_data(cell_count_check,1)=isub
          sub_particle_data(cell_count_check,2)=jsub
          if (SDIM.eq.3) then
           sub_particle_data(cell_count_check,SDIM)=ksub
          endif
          sub_particle_data(cell_count_check,SDIM+1)=current_link
         else
          print *,"sub_box not found"
          stop
         endif
         ibase=(current_link-1)*(1+SDIM)
         current_link=particle_link_data(ibase+1)
        enddo ! while (current_link.ge.1)

        if (cell_count_check.eq.cell_count_hold) then

         cell_count_check=0

         do isub=sublo(1),subhi(1)
         do jsub=sublo(2),subhi(2)
         do ksub=sublo(3),subhi(3)
          ! increment Np_append if isweep == 0
          ! always increment Np_append_test
          local_count=sub_counter(isub,jsub,ksub)

          cell_count_check=cell_count_check+local_count

            ! check if particles need to be deleted
          if (local_count.ge.1) then

           allocate(sort_data_time(local_count))
           allocate(sort_data_id(local_count))
           sub_iter=0
           do cell_iter=1,cell_count_hold
            isub_test=sub_particle_data(cell_iter,1)
            jsub_test=sub_particle_data(cell_iter,2)
            ksub_test=sub_particle_data(cell_iter,SDIM)
            current_link=sub_particle_data(cell_iter,SDIM+1)
            if ((isub_test.eq.isub).and. &
                (jsub_test.eq.jsub)) then
             if ((SDIM.eq.2).or. &
                 ((SDIM.eq.3).and.(ksub_test.eq.ksub))) then
              sub_iter=sub_iter+1
              sort_data_id(sub_iter)=current_link                     
              sort_data_time(sub_iter)= &
                particles(current_link)% &
                extra_state(N_EXTRA_REAL_INSERT_TIME+1) 
              do dir=1,SDIM
               xpart(dir)=particles(current_link)%pos(dir)
              enddo 

              data_in%xtarget=xpart
              call interp_from_grid_util(data_in,lsfab_ptr,data_out_LS)
              do dir=1,num_materials
               LS_sub(dir)=data_out_LS%data_interp(dir)
              enddo
              call get_primary_material(LS_sub,im_primary_sub)
              if ((im_primary_sub.ge.1).and. &
                  (im_primary_sub.le.num_materials)) then
               im_particle=particles(current_link)% &
                 extra_int(N_EXTRA_INT_MATERIAL_ID+1)
               if ((im_particle.ge.1).and.(im_particle.le.num_materials)) then
                if (im_particle.eq.im_primary_sub) then
                 ! do nothing
                else if (im_particle.ne.im_primary_sub) then
                 sort_data_time(sub_iter)=time_in_future
                else
                 print *,"im_particle or im_primary_sub invalid"
                 stop
                endif
               else
                print *,"im_particle invalid"
                stop
               endif
              else
               print *,"im_primary_sub invalid"
               stop
              endif

             else if ((SDIM.eq.3).and.(ksub_test.ne.ksub)) then
              ! do nothing
             else
              print *,"dimension or ksub bust"
              stop
             endif
            endif
           enddo ! cell_iter=1..cell_count_hold

           if (sub_iter.eq.local_count) then
            bubble_change=1
            bubble_iter=0
             ! sort from oldest particle to youngest.
             ! i.e. particle with smallest "add time" is at the top of
             ! the list.
            do while ((bubble_change.eq.1).and. &
                      (bubble_iter.lt.local_count))
             do ibubble=1,local_count-bubble_iter-1
              if (sort_data_time(ibubble).gt. &
                  sort_data_time(ibubble+1)) then
               temp_id=sort_data_id(ibubble)
               sort_data_id(ibubble)=sort_data_id(ibubble+1)
               sort_data_id(ibubble+1)=temp_id
               temp_time=sort_data_time(ibubble)
               sort_data_time(ibubble)=sort_data_time(ibubble+1)
               sort_data_time(ibubble+1)=temp_time
               bubble_change=1
              endif
             enddo ! ibubble=1..local_count-bubble_iter-1
             bubble_iter=bubble_iter+1
            enddo ! bubble_change==1 and bubble_iter<local_count

            do bubble_iter=1,local_count
               ! never delete particles that were present from the
               ! very beginning of the simulation.
             if (sort_data_time(bubble_iter).eq.zero) then
              ! do nothing
             else if (sort_data_time(bubble_iter).gt.zero) then
              if ((bubble_iter.gt.particle_max_per_nsubdivide).or. &
                  (sort_data_time(bubble_iter).ge.time_in_future-one)) then
               particle_delete_flag(sort_data_id(bubble_iter))=1
              else if ((bubble_iter.le.particle_max_per_nsubdivide).and. &
                   (sort_data_time(bubble_iter).lt.time_in_future-one)) then
               ! do nothing
              else
               print *,"bubble_iter or sort_data_time bust"
               stop
              endif
             else
              print *,"sort_data_time(bubble_iter) invalid"
              stop
             endif
            enddo ! bubble_iter
           else
            print *,"sub_iter invalid"
            stop
           endif    
           deallocate(sort_data_time)
           deallocate(sort_data_id)
          else if (local_count.eq.0) then
           ! do nothing
          else
           print *,"local_count bust"
           stop
          endif 

           ! insufficient particles in the subbox or adding the
           ! particles for the very first time.
          if ((local_count.lt.particle_min_per_nsubdivide).or. &
              (append_flag.eq.0)) then 

           call sub_box_cell_center( &
             accum_PARM, &
             i,j,k, &
             isub,jsub,ksub, &
             xsub)

             ! add bulk particles
           call interp_eul_lag_dist( &
             im_elastic_map, &
             accum_PARM, &
             particlesptr, &
             real_compALLptr, &
             tensorfab_ptr, &
             lsfab_ptr, &
             cell_particle_count_ptr, &
             i,j,k, &
             xsub, &
             particle_link_data, &
             Np, &
             tensor_sub, &
             LS_sub, &
             X0_sub)

           do ipart=1,num_materials_viscoelastic
            im_map=im_elastic_map(ipart)+1
            if ((im_map.ge.1).and.(im_map.le.num_materials)) then
             do ii=1,3
             do jj=1,3
              Q(ii,jj)=zero
             enddo
             enddo
             do dir=1,ENUM_NUM_TENSOR_TYPE
              call stress_index(dir,ii,jj)
              Q(ii,jj)=tensor_sub((ipart-1)*ENUM_NUM_TENSOR_TYPE+dir)
             enddo
             Q(2,1)=Q(1,2)
             Q(3,1)=Q(1,3)
             Q(3,2)=Q(2,3)

             if (viscoelastic_model(im_map).eq.3) then ! incremental model
              ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
              ! do nothing
             else
              do ii=1,3
               Q(ii,ii)=Q(ii,ii)+one
              enddo
             endif

             call project_A_to_positive_definite_or_traceless(Q, &
               viscoelastic_model(im_map),polymer_factor(im_map))

             if (viscoelastic_model(im_map).eq.3) then ! incremental model
              ! Maire, Abgrall, Breil, Loubere, Rebourcet JCP 2013
              ! do nothing
             else
              do ii=1,3
               Q(ii,ii)=Q(ii,ii)-one  ! Q <--  A-I
              enddo
             endif

             do dir=1,ENUM_NUM_TENSOR_TYPE
              call stress_index(dir,ii,jj)
              tensor_sub((ipart-1)*ENUM_NUM_TENSOR_TYPE+dir)=Q(ii,jj)
             enddo
            else
             print *,"im_map invalid"
             stop
            endif
           enddo ! ipart=1...num_materials_viscoelastic

           Np_append_test=Np_append_test+1

           if (isweep.eq.0) then
            ! do nothing
           else if (isweep.eq.1) then
            ibase=(Np_append_test-1)*single_particle_size
            do dir=1,SDIM
             new_particles(ibase+dir)=xsub(dir)
            enddo
            do dir=1,NUM_CELL_ELASTIC
             new_particles(ibase+SDIM+N_EXTRA_REAL+N_EXTRA_INT+dir)= &
               tensor_sub(dir)
            enddo
            new_particles(ibase+SDIM+N_EXTRA_REAL_INSERT_TIME+1)=cur_time_slab
            do dir=1,SDIM
             new_particles(ibase+SDIM+N_EXTRA_REAL_X0+dir)=X0_sub(dir)
            enddo
            call get_primary_material(LS_sub,im_primary_sub)
            if ((im_primary_sub.ge.1).and. &
                (im_primary_sub.le.num_materials)) then
             new_particles(ibase+SDIM+N_EXTRA_REAL+N_EXTRA_INT_MATERIAL_ID+1)= &
              im_primary_sub
            else
             print *,"im_primary_sub invalid"
             stop
            endif
           else
            print *,"isweep invalid"
            stop
           endif

          else if ((local_count.ge.particle_min_per_nsubdivide).and. &
                   (append_flag.eq.1)) then
           ! do nothing
          else
           print *,"local_count invalid"
           stop
          endif

         enddo ! ksub
         enddo ! jsub
         enddo ! isub

         if (cell_count_check.eq.cell_count_hold) then
          ! do nothing
         else
          print *,"cell_count_check invalid"
          stop
         endif

        else
         print *,"cell_count_check invalid"
         stop
        endif

        deallocate(sub_particle_data)

       else if (local_mask.eq.0) then
        ! do nothing
       else
        print *,"local_mask invalid"
        stop
       endif

      enddo 
      enddo 
      enddo  ! i,j,k

      if (isweep.eq.0) then
       Np_append=Np_append_test
      else if (isweep.eq.1) then
       if ((Np_append.eq.Np_append_test).and. &
           (new_Pdata_size.eq. &
            Np_append*single_particle_size)) then
        ! do nothing
       else
        print *,"Np_append or new_Pdata_size invalid"
        stop
       endif
      else
       print *,"isweep invalid"
       stop
      endif

      deallocate(sub_counter)

      return
      end subroutine fort_init_particle_container

      subroutine interp_mac_velocity( &
        grid_PARM, &
        umacptr, &
        vmacptr, &
        wmacptr, &
        xpart, &
        vel_time_slab,u)
      use global_utility_module
      use probcommon_module
      use probf90_module

      implicit none

      type(grid_parm_type), INTENT(in) :: grid_PARM
      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:)) :: umacptr
      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:)) :: vmacptr
      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:)) :: wmacptr
      REAL_T, INTENT(in) :: xpart(SDIM)
      REAL_T, INTENT(in) :: vel_time_slab
      REAL_T, INTENT(out) :: u(SDIM)

      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      INTEGER_T imac,jmac,kmac
      INTEGER_T isten,jsten,ksten
      INTEGER_T dir,dir_inner
      INTEGER_T imaclo(3)
      INTEGER_T imachi(3)
      INTEGER_T cell_index(SDIM)
      INTEGER_T, parameter :: nhalf=3
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T xstenMAC_lo(-nhalf:nhalf,SDIM)
      REAL_T xstenMAC_hi(-nhalf:nhalf,SDIM)
      REAL_T dx_inner
      REAL_T wt_dist(SDIM)
      REAL_T local_data
      REAL_T, dimension(D_DECL(2,2,2),1) :: data_stencil
      INTEGER_T ncomp_interp
      REAL_T LS_clamped
      REAL_T vel_clamped(SDIM)
      REAL_T temperature_clamped
      INTEGER_T prescribed_flag
      REAL_T, pointer, dimension(D_DECL(:,:,:)) :: local_data_fab

      if (vel_time_slab.ge.zero) then
       ! do nothing
      else
       print *,"vel_time_slab invalid"
       stop
      endif

      call SUB_clamped_LS(xpart,vel_time_slab,LS_clamped, &
       vel_clamped,temperature_clamped,prescribed_flag,grid_PARM%dx)

      if (LS_clamped.ge.zero) then

       do dir=1,SDIM
        u(dir)=vel_clamped(dir)
       enddo

      else if (LS_clamped.lt.zero) then

       call containing_cell(grid_PARM%bfact, &
          grid_PARM%dx, &
          grid_PARM%xlo, &
          grid_PARM%fablo, &
          xpart, &
          cell_index)

       i=cell_index(1)
       j=cell_index(2)
       k=cell_index(SDIM)

       call gridsten_level(xsten,i,j,k,grid_PARM%level,nhalf)

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

        imaclo(3)=0
        imachi(3)=0
        do dir_inner=1,SDIM
         if (dir_inner.eq.dir) then
          imaclo(dir_inner)=cell_index(dir_inner)
          imachi(dir_inner)=cell_index(dir_inner)+1
         else if (dir_inner.ne.dir) then
          if (xpart(dir_inner).le.xsten(0,dir_inner)) then
           imaclo(dir_inner)=cell_index(dir_inner)-1
           imachi(dir_inner)=cell_index(dir_inner)
          else if (xpart(dir_inner).ge.xsten(0,dir_inner)) then
           imaclo(dir_inner)=cell_index(dir_inner)
           imachi(dir_inner)=cell_index(dir_inner)+1
          else
           print *,"xpart or xsten invalid"
           stop
          endif
         else
          print *,"dir_inner or dir invalid"
          stop
         endif
        enddo ! dir_inner=1..sdim

         ! dir=1..sdim
        call gridstenMAC_level(xstenMAC_lo, &
         imaclo(1),imaclo(2),imaclo(3),grid_PARM%level,nhalf,dir-1)
        call gridstenMAC_level(xstenMAC_hi, &
         imachi(1),imachi(2),imachi(3),grid_PARM%level,nhalf,dir-1)

        do dir_inner=1,SDIM
         dx_inner=xstenMAC_hi(0,dir_inner)-xstenMAC_lo(0,dir_inner)
         if (dx_inner.gt.zero) then
          wt_dist(dir_inner)=(xpart(dir_inner)-xstenMAC_lo(0,dir_inner))/ &
            dx_inner
         else
          print *,"dx_inner invalid"
          stop
         endif
        enddo  ! dir_inner=1..sdim

        do imac=imaclo(1),imachi(1)
        do jmac=imaclo(2),imachi(2)
        do kmac=imaclo(3),imachi(3)

         isten=imac-imaclo(1)+1
         jsten=jmac-imaclo(2)+1
         ksten=kmac-imaclo(3)+1
         if (dir.eq.1) then
          local_data_fab=>umacptr
         else if (dir.eq.2) then
          local_data_fab=>vmacptr
         else if ((dir.eq.3).and.(SDIM.eq.3)) then
          local_data_fab=>wmacptr
         else
          print *,"dir invalid"
          stop
         endif
         call safe_data_single(imac,jmac,kmac,local_data_fab,local_data)

         data_stencil(D_DECL(isten,jsten,ksten),1)=local_data
        enddo ! kmac
        enddo ! jmac
        enddo ! imac
     
        ncomp_interp=1
        call bilinear_interp_stencil(data_stencil, &
          wt_dist,ncomp_interp,u(dir)) 

       enddo ! dir=1..sdim

      else
       print *,"LS_clamped invalid"
       stop
      endif

      end subroutine interp_mac_velocity


      subroutine check_cfl_BC(grid_PARM, xpart1, xpart2)
      use global_utility_module

      implicit none

      type(grid_parm_type), INTENT(in) :: grid_PARM
      REAL_T, INTENT(in) :: xpart1(SDIM)
      REAL_T, INTENT(inout) :: xpart2(SDIM)
      INTEGER_T bc_local
      INTEGER_T dir
      INTEGER_T dir_inner
      REAL_T factor
      REAL_T mag
      REAL_T max_travel

      max_travel=grid_PARM%dx(1)

      do dir=1,SDIM
       if (xpart2(dir).lt.grid_PARM%problo(dir)) then
        bc_local=grid_PARM%velbc(dir,1,dir)
        if (bc_local.eq.REFLECT_ODD) then
         if (xpart1(dir).ge.grid_PARM%problo(dir)) then
          factor=(xpart1(dir)-grid_PARM%problo(dir))/ &
                 (xpart1(dir)-xpart2(dir))
          if ((factor.ge.zero).and. &
              (factor.le.one)) then
           do dir_inner=1,SDIM
            xpart2(dir_inner)=xpart1(dir_inner)+ &
             factor*(xpart2(dir_inner)-xpart1(dir_inner))
           enddo
          else
           print *,"factor invalid"
           stop
          endif
         else
          print *,"xpart1(dir) invalid"
          stop
         endif
        else if ((bc_local.eq.INT_DIR).or. &
                 (bc_local.eq.EXT_DIR).or. &
                 (bc_local.eq.REFLECT_EVEN).or. &
                 (bc_local.eq.FOEXTRAP)) then
         ! do nothing
        else
         print *,"bc_local invalid"
         stop
        endif
       endif

       if (xpart2(dir).gt.grid_PARM%probhi(dir)) then
        bc_local=grid_PARM%velbc(dir,2,dir)
        if (bc_local.eq.REFLECT_ODD) then
         if (xpart1(dir).le.grid_PARM%probhi(dir)) then
          factor=(xpart1(dir)-grid_PARM%probhi(dir))/ &
                 (xpart1(dir)-xpart2(dir))
          if ((factor.ge.zero).and. &
              (factor.le.one)) then
           do dir_inner=1,SDIM
            xpart2(dir_inner)=xpart1(dir_inner)+ &
             factor*(xpart2(dir_inner)-xpart1(dir_inner))
           enddo
          else
           print *,"factor invalid"
           stop
          endif
         else
          print *,"xpart1(dir) invalid"
          stop
         endif
        else if ((bc_local.eq.INT_DIR).or. &
                 (bc_local.eq.EXT_DIR).or. &
                 (bc_local.eq.REFLECT_EVEN).or. &
                 (bc_local.eq.FOEXTRAP)) then
         ! do nothing
        else
         print *,"bc_local invalid"
         stop
        endif
       endif
      enddo !dir=1..sdim

      mag=zero
      do dir=1,SDIM
       mag=mag+(xpart1(dir)-xpart2(dir))**2
      enddo
      mag=sqrt(mag)
      if (mag.gt.max_travel) then
       factor=max_travel/mag
       if ((factor.ge.zero).and.(factor.le.one)) then
        do dir=1,SDIM
         xpart2(dir)=xpart1(dir)+ &
             factor*(xpart2(dir)-xpart1(dir))
        enddo
       else
        print *,"factor invalid"
        stop
       endif
      else if (mag.le.max_travel) then
       ! do nothing
      else
       print *,"mag invalid check_cfl_BC"
       stop
      endif

      end subroutine check_cfl_BC

       ! called from NavierStokes2.cpp
      subroutine fort_move_particle_container( &
        tid, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        level, &
        finest_level, &
        xlo,dx, &
        particles, & ! a list of particles in the elastic structure
        Np, & !  Np = number of particles
        dt, &
        vel_time_slab, &
        umac,DIMS(umac), &
        vmac,DIMS(vmac), &
        wmac,DIMS(wmac), &
        velbc_in, &
        denbc_in, &
        dombc, &
        domlo, &
        domhi) &
      bind(c,name='fort_move_particle_container')

      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: tid
      INTEGER_T, INTENT(in) :: level,finest_level

      REAL_T, INTENT(in) :: dt
      REAL_T, INTENT(in) :: vel_time_slab

      INTEGER_T, INTENT(in), target :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in), target :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, INTENT(in) :: bfact
      REAL_T, INTENT(in), target :: xlo(SDIM),dx(SDIM)
      INTEGER_T, value, INTENT(in) :: Np ! pass by value
      type(particle_t), INTENT(inout), target :: particles(Np)

      INTEGER_T, INTENT(in) :: DIMDEC(umac)
      INTEGER_T, INTENT(in) :: DIMDEC(vmac)
      INTEGER_T, INTENT(in) :: DIMDEC(wmac)

      REAL_T, INTENT(in), target :: umac(DIMV(umac)) 
      REAL_T, pointer, dimension(D_DECL(:,:,:)) :: umac_ptr
      REAL_T, INTENT(in), target :: vmac(DIMV(vmac)) 
      REAL_T, pointer, dimension(D_DECL(:,:,:)) :: vmac_ptr
      REAL_T, INTENT(in), target :: wmac(DIMV(wmac)) 
      REAL_T, pointer, dimension(D_DECL(:,:,:)) :: wmac_ptr

      INTEGER_T, INTENT(in), target :: velbc_in(SDIM,2,SDIM)
      INTEGER_T, INTENT(in) :: denbc_in(SDIM,2)
      INTEGER_T, INTENT(in), target :: dombc(SDIM,2)
      INTEGER_T, INTENT(in), target :: domlo(SDIM)
      INTEGER_T, INTENT(in), target :: domhi(SDIM)

      REAL_T, target :: problo_arr(3)
      REAL_T, target :: probhi_arr(3)

      INTEGER_T interior_ID
      INTEGER_T dir,side,veldir
      REAL_T xpart1(SDIM)
      REAL_T xpart2(SDIM)
      REAL_T xpart3(SDIM)
      REAL_T xpart4(SDIM)
      REAL_T xpart_last(SDIM)
      REAL_T u1(SDIM), u2(SDIM), u3(SDIM), u4(SDIM)
      type(grid_parm_type) grid_PARM
      INTEGER_T num_RK_stages

      REAL_T wrap_pos

      umac_ptr=>umac
      vmac_ptr=>vmac
      wmac_ptr=>wmac

      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      problo_arr(1)=problox
      problo_arr(2)=probloy
      problo_arr(3)=probloz

      probhi_arr(1)=probhix
      probhi_arr(2)=probhiy
      probhi_arr(3)=probhiz

      do dir=1,SDIM
       grid_PARM%fablo(dir)=fablo(dir)
       grid_PARM%fabhi(dir)=fabhi(dir)
       grid_PARM%tilelo(dir)=tilelo(dir)
       grid_PARM%tilehi(dir)=tilehi(dir)
       grid_PARM%dx(dir)=dx(dir)
       grid_PARM%xlo(dir)=xlo(dir)
       do side=1,2
        grid_PARM%dombc(dir,side)=dombc(dir,side)
       enddo
       grid_PARM%domlo(dir)=domlo(dir)
       grid_PARM%domhi(dir)=domhi(dir)
       grid_PARM%problo(dir)=problo_arr(dir)
       grid_PARM%probhi(dir)=probhi_arr(dir)
       do veldir=1,SDIM
        do side=1,2
         grid_PARM%velbc(dir,side,veldir)=velbc_in(dir,side,veldir)
        enddo
       enddo
      enddo ! dir=1,sdim
      grid_PARM%bfact=bfact
      grid_PARM%level=level
      grid_PARM%finest_level=finest_level

      call checkbound_array1(fablo,fabhi,umac_ptr,1,0)
      call checkbound_array1(fablo,fabhi,vmac_ptr,1,1)
      call checkbound_array1(fablo,fabhi,wmac_ptr,1,SDIM-1)

      num_RK_stages=2
      
      do interior_ID=1,Np

       !4th-RK
       do dir=1,SDIM
        xpart1(dir)=particles(interior_ID)%pos(dir)
       enddo
       call interp_mac_velocity( &
        grid_PARM, &
        umac_ptr, &
        vmac_ptr, &
        wmac_ptr, &
        xpart1, &
        vel_time_slab,u1)

       if (num_RK_stages.eq.4) then

        do dir=1,SDIM
         xpart2(dir)=xpart1(dir)+0.5d0*dt*u1(dir)
        enddo
        call check_cfl_BC(grid_PARM,xpart1,xpart2)

        call interp_mac_velocity( &
         grid_PARM, &
         umac_ptr, &
         vmac_ptr, &
         wmac_ptr, &
         xpart2, &
         vel_time_slab,u2)

        do dir=1,SDIM
         xpart3(dir)=xpart1(dir)+0.5d0*dt*u2(dir)
        enddo

        call check_cfl_BC(grid_PARM,xpart1,xpart3)

        call interp_mac_velocity( &
         grid_PARM, &
         umac_ptr, &
         vmac_ptr, &
         wmac_ptr, &
         xpart3, &
         vel_time_slab,u3)

        do dir=1,SDIM
         xpart4(dir)=xpart1(dir)+dt*u3(dir)
        enddo

        call check_cfl_BC(grid_PARM,xpart1,xpart4)

        call interp_mac_velocity( &
         grid_PARM, &
         umac_ptr, &
         vmac_ptr, &
         wmac_ptr, &
         xpart4, &
         vel_time_slab,u4)

        do dir=1,SDIM
         xpart_last(dir)=xpart1(dir)+(1.0d0/6.d0)*dt &
          *(u1(dir)+2.d0*u2(dir)+2.d0*u3(dir)+u4(dir))
        enddo
       else if (num_RK_stages.eq.2) then
        do dir=1,SDIM
         xpart2(dir)=xpart1(dir)+dt*u1(dir)
        enddo
        call check_cfl_BC(grid_PARM,xpart1,xpart2)

        call interp_mac_velocity( &
         grid_PARM, &
         umac_ptr, &
         vmac_ptr, &
         wmac_ptr, &
         xpart2, &
         vel_time_slab,u2)

        do dir=1,SDIM
         xpart_last(dir)=xpart1(dir)+0.5d0*dt &
          *(u1(dir)+u2(dir))
        enddo
       else
        print *,"num_RK_stages invalid"
        stop
       endif

       call check_cfl_BC(grid_PARM,xpart1,xpart_last)

       do dir=1,SDIM
        particles(interior_ID)%pos(dir)=xpart_last(dir)
       enddo

       do dir=1,SDIM
        wrap_pos=particles(interior_ID)%pos(dir)
        if (wrap_pos.lt.problo_arr(dir)) then
         if (dombc(dir,1).eq.INT_DIR) then
          if (dombc(dir,2).eq.INT_DIR) then
           particles(interior_ID)%pos(dir)= &
             wrap_pos+ &
             (probhi_arr(dir)-problo_arr(dir))
          else
           print *,"expecting both dombc_lo and dombc_hi to be INT_DIR"
           stop
          endif
         else if ((dombc(dir,1).eq.EXT_DIR).or. &
                  (dombc(dir,1).eq.REFLECT_EVEN).or. &
                  (dombc(dir,1).eq.FOEXTRAP)) then
          ! do nothing
         else
          print *,"dombc(dir,1) invalid in fort_move_particle_container"
          stop
         endif
        else if (wrap_pos.gt.probhi_arr(dir)) then
         if (dombc(dir,2).eq.INT_DIR) then
          if (dombc(dir,1).eq.INT_DIR) then
           particles(interior_ID)%pos(dir)= &
             wrap_pos- &
             (probhi_arr(dir)-problo_arr(dir))
          else
           print *,"expecting both dombc_lo and dombc_hi to be INT_DIR"
           stop
          endif
         else if ((dombc(dir,2).eq.EXT_DIR).or. &
                  (dombc(dir,2).eq.REFLECT_EVEN).or. &
                  (dombc(dir,2).eq.FOEXTRAP)) then
          ! do nothing
         else
          print *,"dombc(dir,2) invalid in fort_move_particle_container"
          stop
         endif

        else if ((wrap_pos.ge.problo_arr(dir)).and. &
                 (wrap_pos.le.probhi_arr(dir))) then
         ! do nothing
        else
         print *,"wrap_pos probably is NaN: ",wrap_pos
         stop
        endif
       enddo ! dir=1..sdim

      enddo!do interior_ID=1,Np

      return
      end subroutine fort_move_particle_container


      end module FSI_PC_LS_module


