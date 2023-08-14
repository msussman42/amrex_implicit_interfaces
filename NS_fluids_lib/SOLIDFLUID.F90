#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#if (AMREX_SPACEDIM == 2)
#define SDIM 2
#endif
#if (AMREX_SPACEDIM == 3)
#define SDIM 3
#endif

#include "AMReX_FORT_INTEGER.H"
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "AMReX_ArrayLim_SUSSMAN.H"

#include "EXTRAP_COMP.H"
#include "SOLIDFLUID_F.H"

      subroutine SOLIDFLUID_F90_KEYBOARD()
      use iso_c_binding
      IMPLICIT NONE

      interface 
      subroutine main_cpp_keyboard() bind(c)
      end subroutine main_cpp_keyboard
      end interface 

      call main_cpp_keyboard()

      return
      end subroutine SOLIDFLUID_F90_KEYBOARD

      module solidfluid_module

      contains

      subroutine init_3D_map(xmap3D,xslice3D,problo3D,probhi3D, &
        dx_max_level)
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, INTENT(out) :: xmap3D(3)
      REAL_T, INTENT(out) :: xslice3D(3)
      REAL_T, INTENT(out) :: problo3D(3),probhi3D(3)
      REAL_T, INTENT(in) :: dx_max_level(SDIM)
      INTEGER_T dir

      do dir=1,SDIM
       if (probhi_array(dir)-problo_array(dir).gt.zero) then
        ! do nothing
       else
        print *,"probhi_array(dir)-problo_array(dir) invalid"
        stop
       endif
       if (dx_max_level(dir).gt.zero) then
        ! do nothing
       else
        print *,"dx_max_level(dir).le.zero"
        stop
       endif
      enddo

      do dir=1,3
       xmap3D(dir)=dir
       xslice3D(dir)=zero
      enddo

      if (SDIM.eq.2) then

       call SUB_FSI_SLICE(xmap3D,xslice3D,problo3D,probhi3D, &
              dx_max_level(1))

      else if (SDIM.eq.3) then

       do dir=1,SDIM
        problo3D(dir)=problo_array(dir)
        probhi3D(dir)=probhi_array(dir)
       enddo

      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine init_3D_map

      end module solidfluid_module

      module solidfluid_cpp_module
      contains

       !FSI_operation.eq.OP_FSI_INITIALIZE_NODES  
       !  initialize node locations; post_process_nodes_elements
       !FSI_operation.eq.OP_FSI_UPDATE_NODES  update node locations
       !FSI_operation.eq.OP_FSI_MAKE_DISTANCE  make distance in narrow band
       !FSI_operation.eq.OP_FSI_MAKE_SIGN  update the sign.
       !FSI_operation.eq.OP_FSI_LAG_STRESS  copy eul fluid stress to solid
       !  FSI_sub_operation.eq.SUB_OP_FSI_CLEAR_LAG_DATA
       !  FSI_sub_operation.eq.SUB_OP_FSI_COPY_TO_LAG_DATA
       !  FSI_sub_operation.eq.SUB_OP_FSI_SYNC_LAG_DATA
      subroutine fort_headermsg( &
        tid, &
        tilenum, &
        gridno, &
        nthread_parm, &
        level, &
        finest_level, &
        max_level, &
        FSI_input_flattened, &
        FSI_output_flattened, &
        flatten_size, &
        local_caller_id, &
        FSI_operation, &
        FSI_sub_operation, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        xlo, & !problo if FSI_operation.eq.OP_FSI_INITIALIZE(UPDATE)_NODES
        dx, & !problen if FSI_operation.eq.OP_FSI_INITIALIZE(UPDATE)_NODES
        dx_max_level, & 
        velbc, &
        vofbc, &
        FSIdata, & ! drag if FSI_operation.eq.OP_FSI_LAG_STRESS
        DIMS(FSIdata), &
        drag, &
        DIMS(drag), &
        masknbr, &
        DIMS(masknbr), &
        maskfiner, &
        DIMS(maskfiner), &
        nFSI, &
        ngrow_make_distance_in, &
        nparts, &
        im_solid_map, & ! type: 0..num_materials-1
        h_small, &  ! smallest mesh size from the max_level.
        cur_time, &
        dt, &
        FSI_refine_factor, & ! 1..num_materials
        FSI_bounding_box_ngrow, & ! 1..num_materials
        touch_flag, &
        CTML_FSI_INIT, &
        iter, &
        current_step, &
        plot_interval, &
        ioproc) &
      bind(c,name='fort_headermsg')

      use CLSVOFCouplerIO, only : CLSVOF_ReadHeader,CLSVOF_ReadNodes, &
       CLSVOF_InitBox,CLSVOF_clear_lag_data,CLSVOF_sync_lag_data, &
       CLSVOF_Copy_To_LAG
      use solidfluid_module
      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: tid
      INTEGER_T, INTENT(in) :: tilenum
      INTEGER_T, INTENT(in) :: gridno
      INTEGER_T, INTENT(in) :: nthread_parm
      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: max_level

      INTEGER_T, INTENT(in) :: flatten_size
      INTEGER_T, INTENT(in) :: local_caller_id
      REAL_T, INTENT(inout) :: FSI_input_flattened(flatten_size)
      REAL_T, INTENT(inout) :: FSI_output_flattened(flatten_size)

      INTEGER_T, INTENT(in) :: FSI_operation
      INTEGER_T, INTENT(in) :: FSI_sub_operation
      INTEGER_T, INTENT(in) :: nFSI
      INTEGER_T, INTENT(in) :: ngrow_make_distance_in
      INTEGER_T, INTENT(in) :: nparts
      INTEGER_T, INTENT(in) :: im_solid_map(nparts)
      REAL_T, INTENT(in) :: h_small ! smallest mesh size from the max_level
      REAL_T, INTENT(in) :: cur_time
      REAL_T, INTENT(in) :: dt
      INTEGER_T, INTENT(in) :: FSI_refine_factor(num_materials)
      INTEGER_T, INTENT(in) :: FSI_bounding_box_ngrow(num_materials)
      INTEGER_T, INTENT(inout) :: touch_flag
      INTEGER_T, INTENT(in) :: CTML_FSI_init
      INTEGER_T, INTENT(in) :: iter
      INTEGER_T, INTENT(in) :: current_step 
      INTEGER_T, INTENT(in) :: plot_interval 
      INTEGER_T, INTENT(in) :: ioproc
      INTEGER_T isout
!drag if FSI_operation.eq.OP_FSI_LAG_STRESS
      INTEGER_T, INTENT(in) :: DIMDEC(FSIdata) 
      INTEGER_T, INTENT(in) :: DIMDEC(drag) 
      INTEGER_T, INTENT(in) :: DIMDEC(masknbr) 
      INTEGER_T, INTENT(in) :: DIMDEC(maskfiner) 
!drag if FSI_operation.eq.OP_FSI_LAG_STRESS
      REAL_T, INTENT(inout), target :: FSIdata(DIMV(FSIdata),nFSI) 
      REAL_T, pointer :: FSIdata_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: drag(DIMV(drag),N_DRAG)
      REAL_T, pointer :: drag_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: masknbr(DIMV(masknbr),2)
      REAL_T, pointer :: masknbr_ptr(D_DECL(:,:,:),:)
      REAL_T, INTENT(in), target :: maskfiner(DIMV(maskfiner),4)
      REAL_T, pointer :: maskfiner_ptr(D_DECL(:,:,:),:)

      INTEGER_T, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T tilelo3D(3),tilehi3D(3)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T growlo3D(3),growhi3D(3)
      INTEGER_T, INTENT(in) :: bfact
      REAL_T, INTENT(in) :: xlo(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: dx_max_level(SDIM)

      REAL_T problo3D(3),probhi3D(3)
      REAL_T dx3D(3)
      REAL_T vel3D(3)
      INTEGER_T, INTENT(in) :: velbc(SDIM,2)
      INTEGER_T, INTENT(in) :: vofbc(SDIM,2)
      INTEGER_T xmap3D(3)
      REAL_T xslice3D(3)
      INTEGER_T dir
      INTEGER_T im_part

      INTEGER_T i,j,k,nc
      INTEGER_T i2d,j2d,k2d

      INTEGER_T idx(3)
      INTEGER_T FSI_lo3D(3),FSI_hi3D(3)
      INTEGER_T FSI_growlo3D(3),FSI_growhi3D(3)
      INTEGER_T DIMDEC3D(FSIdata3D)
      REAL_T, allocatable,target :: FSIdata3D(:,:,:,:)
      REAL_T, pointer :: FSIdata3D_ptr(:,:,:,:)

       ! i,j,k,6*num_materials
      REAL_T, allocatable,target :: stressdata3D(:,:,:,:)
      REAL_T, pointer :: stressdata3D_ptr(:,:,:,:)

       ! i,j,k,num_materials
      REAL_T, allocatable,target :: stressflag3D(:,:,:,:)
      REAL_T, pointer :: stressflag3D_ptr(:,:,:,:)

      REAL_T, allocatable,target :: xdata3D(:,:,:,:)
      REAL_T, pointer :: xdata3D_ptr(:,:,:,:)
      REAL_T, allocatable,target :: masknbr3D(:,:,:,:)
      REAL_T, pointer :: masknbr3D_ptr(:,:,:,:)
      REAL_T, allocatable,target :: maskfiner3D(:,:,:,:)
      REAL_T, pointer :: maskfiner3D_ptr(:,:,:,:)
      INTEGER_T mask1,mask2
      INTEGER_T, parameter :: nhalf=3
      REAL_T xsten(-nhalf:nhalf,SDIM)
      INTEGER_T ibase
      INTEGER_T partid
      INTEGER_T lev77
      INTEGER_T im_local
      INTEGER_T istress,jstress
      REAL_T stress_2d(3,3)
      REAL_T stress_3d(3,3)
      REAL_T xlo3D_tile(3)
      REAL_T xhi3D_tile(3)
 
      FSIdata_ptr=>FSIdata

      if (nthread_parm.ne.geom_nthreads) then
       print *,"nthread_parm invalid"
       stop
      endif
      
      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif
      if (tilenum.lt.0) then
       print *,"tilenum invalid"
       stop
      endif
      if (gridno.lt.0) then
       print *,"gridno invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid300"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in headermsg"
       stop
      endif
      if (finest_level.le.max_level) then
       ! do nothing
      else
       print *,"finest_level or max_level invalid"
       stop
      endif
      lev77=level+1
      if ((lev77.lt.1).or.(lev77.gt.finest_level+1)) then
       print *,"lev77 invalid in headermsg"
       stop
      endif
      if ((ioproc.ne.0).and.(ioproc.ne.1)) then
       print *,"ioproc invalid"
       stop
      endif
      if (ngrow_make_distance.ne.3) then
       print *,"ngrow_make_distance invalid"
       stop
      endif
      if ((FSI_operation.eq.OP_FSI_INITIALIZE_NODES).or. &  
          (FSI_operation.eq.OP_FSI_UPDATE_NODES)) then  
       if (ngrow_make_distance_in.ne.0) then
        print *,"ngrow_make_distance_in invalid"
        print *,"fort_headermsg"
        print *,"ngrow_make_distance_in=",ngrow_make_distance_in
        stop
       endif
       if ((local_caller_id.eq.caller_initData).or. &
           (local_caller_id.eq.caller_post_restart).or. &
           (local_caller_id.eq.caller_nonlinear_advection))then
        !do nothing
       else
        print *,"local_caller_id invalid in fort_headermsg(1): ", &
          local_caller_id
        stop
       endif
      else if (FSI_operation.eq.OP_FSI_LAG_STRESS) then 

       if ((FSI_sub_operation.eq.SUB_OP_FSI_CLEAR_LAG_DATA).or. &
           (FSI_sub_operation.eq.SUB_OP_FSI_SYNC_LAG_DATA)) then
        if (ngrow_make_distance_in.ne.0) then
         print *,"expecting ngrow_make_distance_in==0"
         print *,"fort_headermsg 2a"
         print *,"ngrow_make_distance_in=",ngrow_make_distance_in
         stop
        endif
       else if (FSI_sub_operation.eq.SUB_OP_FSI_COPY_TO_LAG_DATA) then
        if (ngrow_make_distance_in.ne.3) then
         print *,"expecting ngrow_make_distance_in==3"
         print *,"fort_headermsg 2b"
         print *,"ngrow_make_distance_in=",ngrow_make_distance_in
         stop
        endif
       else
        print *,"FSI_sub_operation invalid"
        stop
       endif

       if (local_caller_id.eq.caller_nonlinear_advection) then
        !do nothing
       else
        print *,"local_caller_id invalid in fort_headermsg(2): ", &
          local_caller_id
        stop
       endif

      else if (FSI_operation.eq.OP_FSI_MAKE_DISTANCE) then 
       if (local_caller_id.eq.caller_FSI_make_distance) then
        !do nothing
       else
        print *,"local_caller_id invalid in fort_headermsg(3): ", &
          local_caller_id
        stop
       endif
      else if (FSI_operation.eq.OP_FSI_MAKE_SIGN) then 
       if (local_caller_id.eq.caller_FSI_make_distance) then
        !do nothing
       else
        print *,"local_caller_id invalid in fort_headermsg(4): ", &
          local_caller_id
        stop
       endif
      else
       print *,"FSI_operation invalid"
       stop
      endif
      if (h_small.gt.zero) then
       ! do nothing
      else
       print *,"h_small invalid: ",h_small
       stop
      endif
      do dir=1,SDIM
       if (dx_max_level(dir).gt.zero) then
        ! do nothing
       else
        print *,"dx_max_level(dir).le.zero"
        stop
       endif
      enddo
      if ((touch_flag.ne.0).and.(touch_flag.ne.1)) then
       print *,"touch_flag invalid in fort_headermsg"
       stop
      endif
      if ((CTML_FSI_init.ne.0).and. &
          (CTML_FSI_init.ne.1)) then
       print *,"CTML_FSI_init invalid"
       stop
      endif
      do im_local=1,num_materials
       if (FSI_flag(im_local).eq.FSI_SHOELE_CTML) then 
        ! do nothing
       else if (fort_FSI_flag_valid(im_local).eq.1) then
        ! do nothing
       else
        print *,"fort_FSI_flag_valid invalid"
        stop
       endif
      enddo ! im_local=1..num_materials

      if (iter.lt.0) then
       print *,"iter invalid"
       stop
      endif
      if (current_step.lt.0) then
       print *,"current_step invalid"
       stop
      endif
      if (plot_interval.lt.-1) then
       print *,"plot_interval invalid"
       stop
      endif

      if ((nparts.lt.1).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_headermsg"
       stop
      endif
 
      if (nFSI.ne.nparts*NCOMP_FSI) then 
       print *,"nFSI invalid"
       stop
      endif
  
      call checkbound_array(fablo,fabhi,FSIdata_ptr, &
        ngrow_make_distance_in,-1)
      drag_ptr=>drag
      call checkbound_array(fablo,fabhi,drag_ptr, &
        ngrow_make_distance_in,-1)
      masknbr_ptr=>masknbr
      call checkbound_array(fablo,fabhi,masknbr_ptr, &
        ngrow_make_distance_in,-1)
      maskfiner_ptr=>maskfiner
      call checkbound_array(fablo,fabhi,maskfiner_ptr, &
       ngrow_make_distance_in,-1)

      ! update ngrow_make_distance grow layers of FSIdata that do not overlap
      ! with another tile.
      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,ngrow_make_distance_in)
       ! since PCINTERP_fill_borders interpolates from coarser levels, 
       ! we only have to traverse interior values.
      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       tilelo3D,tilehi3D,0)

      do dir=1,SDIM
       if (probhi_array(dir)-problo_array(dir).gt.zero) then
        ! do nothing
       else
        print *,"probhi_array(dir)-problo_array(dir) invalid"
        stop
       endif
      enddo

      call init_3D_map(xmap3D,xslice3D,problo3D,probhi3D,dx_max_level)

      if (SDIM.eq.2) then

       do dir=1,3
        if (xmap3D(dir).eq.0) then
         dx3D(dir)=dx_max_level(1)
         FSI_lo3D(dir)=0
         FSI_hi3D(dir)=0

         xlo3D_tile(dir)=xslice3D(dir)-half*dx3D(dir)
         xhi3D_tile(dir)=xlo3D_tile(dir)+dx3D(dir)

         growlo3D(dir)=-ngrow_make_distance_in
         growhi3D(dir)=ngrow_make_distance_in
        else if ((xmap3D(dir).eq.1).or. &
                 (xmap3D(dir).eq.2)) then
         dx3D(dir)=dx(xmap3D(dir))
         problo3D(dir)=problo_array(xmap3D(dir))
         probhi3D(dir)=probhi_array(xmap3D(dir))
         FSI_lo3D(dir)=tilelo(xmap3D(dir))
         FSI_hi3D(dir)=tilehi(xmap3D(dir))

         xlo3D_tile(dir)=xlo(xmap3D(dir))+ &
           dx3D(dir)*(FSI_lo3D(dir)-fablo(xmap3D(dir)))
         xhi3D_tile(dir)=xlo3D_tile(dir)+ &
           dx3D(dir)*(FSI_hi3D(dir)-FSI_lo3D(dir)+1)

         growlo3D(dir)=growlo(xmap3D(dir))
         growhi3D(dir)=growhi(xmap3D(dir))
        else
         print *,"xmap3D(dir) invalid"
         stop
        endif
       enddo ! dir=1..3

      else if (SDIM.eq.3) then

       do dir=1,SDIM
        dx3D(dir)=dx(dir)
        FSI_lo3D(dir)=tilelo(dir)
        FSI_hi3D(dir)=tilehi(dir)

        xlo3D_tile(dir)=xlo(dir)+ &
          dx3D(dir)*(FSI_lo3D(dir)-fablo(dir))
        xhi3D_tile(dir)=xlo3D_tile(dir)+ &
          dx3D(dir)*(FSI_hi3D(dir)-FSI_lo3D(dir)+1)

        growlo3D(dir)=growlo(dir)
        growhi3D(dir)=growhi(dir)
       enddo ! dir=1..SDIM

      else
       print *,"dimension bust"
       stop
      endif

      do dir=1,3
       FSI_growlo3D(dir)=FSI_lo3D(dir)-ngrow_make_distance_in
       FSI_growhi3D(dir)=FSI_hi3D(dir)+ngrow_make_distance_in
      enddo

      ARG3D_L1(FSIdata3D)=FSI_growlo3D(1)
      ARG3D_L2(FSIdata3D)=FSI_growlo3D(2)
      ARG3D_L3(FSIdata3D)=FSI_growlo3D(3)
      ARG3D_H1(FSIdata3D)=FSI_growhi3D(1)
      ARG3D_H2(FSIdata3D)=FSI_growhi3D(2)
      ARG3D_H3(FSIdata3D)=FSI_growhi3D(3)

      if ((FSI_operation.eq.OP_FSI_INITIALIZE_NODES).or. &  
          (FSI_operation.eq.OP_FSI_UPDATE_NODES)) then  

       if ((tilenum.ne.0).or.(gridno.ne.0)) then
        print *,"tilenum or gridno invalid"
        print *,"tilenum: ",tilenum
        print *,"gridno: ",gridno
        stop
       endif

       if (FSI_sub_operation.ne.SUB_OP_FSI_DEFAULT) then
        print *,"FSI_sub_operation.ne.SUB_OP_FSI_DEFAULT"
        stop
       endif

       isout=1 ! verbose on in sci_clsvof.F90
       if (FSI_operation.eq.OP_FSI_INITIALIZE_NODES) then 

        call CLSVOF_ReadHeader( &
          FSI_input_flattened, &
          FSI_output_flattened, &
          flatten_size, &
          local_caller_id, &
          FSI_refine_factor, &
          FSI_bounding_box_ngrow, &
          nparts, &
          im_solid_map, &
          h_small, &
          dx_max_level, &
          CTML_FSI_INIT, &
          cur_time, &
          problo3D,probhi3D, &
          ioproc,isout)

       else if (FSI_operation.eq.OP_FSI_UPDATE_NODES) then 

        if (CTML_FSI_INIT.ne.1) then
         print *,"CTML_FSI_INIT.ne.1"
         stop
        endif
         ! cur_time=t^{n+1}
         ! if FSI_flag==FSI_SHOELE_CTML, then
         !  tick_fib is called (in ../StructureCodeShoele/tick.F)
        call CLSVOF_ReadNodes( &
          FSI_input_flattened, &
          FSI_output_flattened, &
          flatten_size, &
          local_caller_id, &
          FSI_refine_factor, &
          FSI_bounding_box_ngrow, &
          cur_time, &
          dt, &
          h_small, &
          problo3D,probhi3D, &
          current_step,plot_interval, &
          ioproc,isout)
       else
        print *,"FSI_operation invalid: ",FSI_operation
        stop
       endif

      else if ((FSI_operation.eq.OP_FSI_MAKE_DISTANCE).or. & 
               (FSI_operation.eq.OP_FSI_MAKE_SIGN)) then 

       isout=1 ! verbose on in sci_clsvof.F90

       if (FSI_sub_operation.ne.SUB_OP_FSI_DEFAULT) then
        print *,"FSI_sub_operation.ne.SUB_OP_FSI_DEFAULT"
        stop
       endif
       if (CTML_FSI_INIT.ne.1) then
        print *,"CTML_FSI_INIT.ne.1"
        stop
       endif
       if (nFSI.ne.nparts*NCOMP_FSI) then 
        print *,"nFSI invalid"
        stop
       endif
 
       allocate(FSIdata3D(DIMV3D(FSIdata3D),nFSI))
       FSIdata3D_ptr=>FSIdata3D
       allocate(xdata3D(DIMV3D(FSIdata3D),3))
       xdata3D_ptr=>xdata3D
       allocate(masknbr3D(DIMV3D(FSIdata3D),2))
       masknbr3D_ptr=>masknbr3D

       do i=FSI_growlo3D(1),FSI_growhi3D(1)
       do j=FSI_growlo3D(2),FSI_growhi3D(2)
       do k=FSI_growlo3D(3),FSI_growhi3D(3)

        idx(1)=i
        idx(2)=j
        idx(3)=k

        i2d=i
        j2d=j
        k2d=k

        if (SDIM.eq.3) then
         if ((k2d.lt.tilelo(SDIM)-ngrow_make_distance_in).or. &
             (k2d.gt.tilehi(SDIM)+ngrow_make_distance_in)) then
          print *,"k2d out of range"
          stop
         endif
        else if (SDIM.eq.2) then
         k2d=0
          ! dir is the coordinate on the 3D grid
          ! xmap3D(dir) is the coordinate on the 2D grid
          ! idx(1)=i  idx(2)=j   idx(3) unused
         do dir=1,3
          if (xmap3D(dir).eq.0) then
           ! do nothing
          else if (xmap3D(dir).eq.1) then
           i2d=idx(dir)
          else if (xmap3D(dir).eq.2) then
           j2d=idx(dir)
          else 
           print *,"xmap3D invalid"
           stop
          endif
         enddo ! dir=1..3
        else
         print *,"dimension bust"
         stop
        endif

        if ((i2d.lt.tilelo(1)-ngrow_make_distance_in).or. &
            (i2d.gt.tilehi(1)+ngrow_make_distance_in)) then
         print *,"i2d out of range"
         stop
        endif
        if ((j2d.lt.tilelo(2)-ngrow_make_distance_in).or. &
            (j2d.gt.tilehi(2)+ngrow_make_distance_in)) then
         print *,"j2d out of range"
         stop
        endif

        call gridsten_level(xsten,i2d,j2d,k2d,level,nhalf)

        do partid=1,nparts

         ibase=(partid-1)*NCOMP_FSI

         FSIdata3D(i,j,k,ibase+FSI_LEVELSET+1)= &
          FSIdata(D_DECL(i2d,j2d,k2d),ibase+FSI_LEVELSET+1) !LS
         FSIdata3D(i,j,k,ibase+FSI_SIGN_CONFLICT+1)= &
          FSIdata(D_DECL(i2d,j2d,k2d),ibase+FSI_SIGN_CONFLICT+1) 
         FSIdata3D(i,j,k,ibase+FSI_TEMPERATURE+1)= &
          FSIdata(D_DECL(i2d,j2d,k2d),ibase+FSI_TEMPERATURE+1) !T
         FSIdata3D(i,j,k,ibase+FSI_EXTRAP_FLAG+1)= &
          FSIdata(D_DECL(i2d,j2d,k2d),ibase+FSI_EXTRAP_FLAG+1) !flag

         ! dir is the coordinate on the 3D grid
         ! xmap3D(dir) is the coordinate on the 2D grid
         do dir=1,3
          if (SDIM.eq.3) then
           vel3D(dir)= &
            FSIdata(D_DECL(i,j,k),ibase+FSI_VELOCITY+dir)
          else if (SDIM.eq.2) then
           if (xmap3D(dir).eq.0) then
            vel3D(dir)=zero
           else if ((xmap3D(dir).eq.1).or. &
                    (xmap3D(dir).eq.2)) then
            vel3D(dir)= &
             FSIdata(D_DECL(i2d,j2d,k2d),ibase+FSI_VELOCITY+xmap3D(dir))
           else
            print *,"xmap3D(dir) invalid"
            stop
           endif
          else
           print *,"dimension bust"
           stop
          endif
           !FSI_VELOCITY is extrapolated.
          FSIdata3D(i,j,k,ibase+FSI_VELOCITY+dir)=vel3D(dir)
         enddo ! dir=1..3

        enddo ! partid=1,nparts
          
        do nc=1,2
         masknbr3D(i,j,k,nc)=masknbr(D_DECL(i2d,j2d,k2d),nc)
        enddo

        do dir=1,3
         if (SDIM.eq.3) then
          xdata3D(i,j,k,dir)=xsten(0,dir)
         else if (SDIM.eq.2) then
          if (xmap3D(dir).eq.0) then
           xdata3D(i,j,k,dir)=xslice3D(dir)+idx(dir)*dx_max_level(1)
          else if (xmap3D(dir).eq.1) then 
           xdata3D(i,j,k,dir)=xsten(0,1)
          else if (xmap3D(dir).eq.2) then 
           xdata3D(i,j,k,dir)=xsten(0,2)
          else 
           print *,"xmap3D invalid"
           stop
          endif
         else
          print *,"dimension bust"
          stop
         endif
        enddo ! dir=1..3

       enddo !k
       enddo !j
       enddo !i

       do partid=1,nparts
        im_part=im_solid_map(partid)+1
        if ((im_part.lt.1).or.(im_part.gt.num_materials)) then
         print *,"im_part invalid fort_headermsg"
         stop
        endif
        if (is_lag_part(im_part).eq.1) then

         if (fort_read_from_CAD(FSI_flag(im_part)).eq.1) then

          if (container_allocated.ne.1) then
           print *,"container_allocated.ne.1"
           stop
          endif
          if (level_container_allocated(lev77).ne.1) then
           print *,"level_container_allocated(lev77).ne.1"
           stop
          endif
          if (tilenum+1.gt.contain_elem(lev77)% &
              num_tiles_on_thread3D_proc(tid+1)) then
           print *,"tilenum+1.gt.num_tiles_on_thread3D_proc(tid+1)"
           stop
          endif
          if (contain_elem(lev77)%gridno3D(tid+1,tilenum+1).ne.gridno) then
           print *,"gridno3D(tid+1,tilenum+1).ne.gridno (2)"
           print *,"gridno= ",gridno
           print *,"lev77= ",lev77
           print *,"tid= ",tid
           print *,"tilenum= ",tilenum
           print *,"contain_elem(lev77)%gridno3D(tid+1,tilenum+1)=", &
            contain_elem(lev77)%gridno3D(tid+1,tilenum+1)
           stop
          endif
          do dir=1,3
           if (contain_elem(lev77)%tilelo3D(tid+1,tilenum+1,dir).ne. &
               FSI_lo3D(dir)) then
            print *,"tilelo3D(tid+1,tilenum+1,dir).ne.FSI_lo3D(dir)"
            stop
           endif 
           if (contain_elem(lev77)%tilehi3D(tid+1,tilenum+1,dir).ne. &
               FSI_hi3D(dir)) then
            print *,"tilehi3D(tid+1,tilenum+1,dir).ne.FSI_hi3D(dir)"
            stop
           endif 
          enddo ! dir=1..3
           ! declared in: sci_clsvof.F90
          call CLSVOF_InitBox( &
           iter, &
           SDIM, &
           lev77, &
           tid, &
           tilenum, &
           im_part, &
           nparts, &
           partid, &
           ngrow_make_distance, &
           nFSI, &
           FSI_operation, &
           touch_flag, &
           cur_time, &
           dt, &
           problo3D,probhi3D, &
           xmap3D, &
           dx3D, &
           xlo3D_tile, &
           xhi3D_tile, &
           FSI_lo3D,FSI_hi3D, &
           FSI_growlo3D,FSI_growhi3D, &
           growlo3D,growhi3D, &
           xdata3D_ptr, &
           FSIdata3D_ptr, &
           masknbr3D_ptr, &
           ioproc,isout)

         else if (FSI_flag(im_part).eq.FSI_PRESCRIBED_PROBF90) then 

          ! do nothing

         else
          print *,"FSI_flag invalid in fort_headermsg"
          print *,"im_part,FSI_flag(im_part) ",im_part,FSI_flag(im_part)
          stop
         endif

        else
         print *,"is_lag_part invalid"
         stop
        endif

       enddo ! partid=1..nparts

       if (nparts.gt.num_materials-1) then
        print *,"nparts out of range"
        stop
       endif

        ! update ngrow_make_distance grow layers of FSIdata that do not overlap
        ! with another tile.
        ! PCINTERP_fill_borders interpolates from coarser levels so
        ! that we only have to traverse interior values here.
       do i=tilelo3D(1),tilehi3D(1)
       do j=tilelo3D(2),tilehi3D(2)
       do k=tilelo3D(3),tilehi3D(3)

        mask1=NINT(masknbr(D_DECL(i,j,k),1))
        mask2=NINT(masknbr(D_DECL(i,j,k),2))

         ! mask2==1 => (i,j,k) in the interior of the tile.
         ! mask1==0 => (i,j,k) in coarse/fine or EXT_DIR ghost cell
        if ((mask1.eq.0).or.(mask2.eq.1)) then

         if (nFSI.ne.nparts*NCOMP_FSI) then 
          print *,"nFSI invalid"
          stop
         endif

         ! dir is the coordinate on the 3D grid
         ! xmap3D(dir) is the coordinate on the 2D grid
         ! idx is the 3d index.
         ! (i,j,k) is the 2D index.
         if (SDIM.eq.3) then
          idx(1)=i
          idx(2)=j
          idx(SDIM)=k
         else if (SDIM.eq.2) then
          do dir=1,3
           if (xmap3D(dir).eq.0) then
            idx(dir)=0
           else if (xmap3D(dir).eq.1) then
            idx(dir)=i
           else if (xmap3D(dir).eq.2) then
            idx(dir)=j
           else 
            print *,"xmap3D invalid"
            stop
           endif
          enddo
         else
          print *,"SDIM invalid"
          stop
         endif

         ! idx is the 3d index.
         ! (i,j,k) is the 2D index.
         do dir=1,3
          if ((idx(dir).lt.FSI_growlo3D(dir)).or. &
              (idx(dir).gt.FSI_growhi3D(dir))) then
           print *,"idx(dir) out of range"
           stop
          endif
         enddo ! dir=1..3

         do partid=1,nparts
          ibase=(partid-1)*NCOMP_FSI
          FSIdata(D_DECL(i,j,k),ibase+FSI_LEVELSET+1)= &
           FSIdata3D(idx(1),idx(2),idx(3),ibase+FSI_LEVELSET+1) ! LS
          FSIdata(D_DECL(i,j,k),ibase+FSI_SIGN_CONFLICT+1)= &
           FSIdata3D(idx(1),idx(2),idx(3),ibase+FSI_SIGN_CONFLICT+1) 
          FSIdata(D_DECL(i,j,k),ibase+FSI_TEMPERATURE+1)= &
           FSIdata3D(idx(1),idx(2),idx(3),ibase+FSI_TEMPERATURE+1) ! T
          FSIdata(D_DECL(i,j,k),ibase+FSI_EXTRAP_FLAG+1)= &
           FSIdata3D(idx(1),idx(2),idx(3),ibase+FSI_EXTRAP_FLAG+1) ! flag

          do dir=1,3
           if (SDIM.eq.3) then
            FSIdata(D_DECL(i,j,k),ibase+FSI_VELOCITY+dir)= &
             FSIdata3D(i,j,k,ibase+FSI_VELOCITY+dir) 
           else if (SDIM.eq.2) then
            if (xmap3D(dir).eq.0) then
             FSIdata(D_DECL(i,j,k),ibase+FSI_VELOCITY+3)=zero
            else if ((xmap3D(dir).eq.1).or. &
                     (xmap3D(dir).eq.2)) then
             FSIdata(D_DECL(i,j,k),ibase+FSI_VELOCITY+xmap3D(dir))= &
              FSIdata3D(idx(1),idx(2),idx(3),ibase+FSI_VELOCITY+dir) 
            else
             print *,"xmap3D(dir) invalid"
             stop
            endif
           else
            print *,"SDIM invalid"
            stop
           endif
          enddo ! dir=1..3

         enddo ! partid=1..nparts
 
        else if ((mask1.eq.1).and.(mask2.eq.0)) then
         ! do nothing
        else
         print *,"mask1 or mask2 invalid"
         stop
        endif

       enddo !k
       enddo !j
       enddo !i

       deallocate(xdata3D)
       deallocate(FSIdata3D)
       deallocate(masknbr3D)

      else if (FSI_operation.eq.OP_FSI_LAG_STRESS) then 

       isout=1 ! verbose on in sci_clsvof.F90

       if (CTML_FSI_INIT.ne.1) then
        print *,"CTML_FSI_INIT.ne.1"
        stop
       endif

        ! called on all of the processors
       if (FSI_sub_operation.eq.SUB_OP_FSI_CLEAR_LAG_DATA) then
         ! ctml_fib_frc=0.0d0
         !  for inode=1,num_nodes and dir=1..3,
         ! FSI(part_id)%NodeVel(dir,inode)=0.0 
         ! FSI(part_id)%NodeVel_old(dir,inode)=0.0 
         ! FSI(part_id)%NodeVel_new(dir,inode)=0.0 
         ! FSI(part_id)%NodeMass(inode)=1.0 
         ! FSI(part_id)%NodeDensity(inode)=1.0 
         ! FSI(part_id)%NodeForce(dir,inode)=0.0 
         ! FSI(part_id)%NodeForce_old(dir,inode)=0.0 
         ! FSI(part_id)%NodeForce_new(dir,inode)=0.0 
        call CLSVOF_clear_lag_data(ioproc,isout)

        ! called only by the processors which own a given FAB.
       else if (FSI_sub_operation.eq.SUB_OP_FSI_COPY_TO_LAG_DATA) then 

        allocate(stressdata3D(DIMV3D(FSIdata3D),6*num_materials)) 
        stressdata3D_ptr=>stressdata3D

        allocate(stressflag3D(DIMV3D(FSIdata3D),num_materials)) 
        stressflag3D_ptr=>stressflag3D

        allocate(xdata3D(DIMV3D(FSIdata3D),3))
        xdata3D_ptr=>xdata3D
        allocate(masknbr3D(DIMV3D(FSIdata3D),2))
        masknbr3D_ptr=>masknbr3D
        allocate(maskfiner3D(DIMV3D(FSIdata3D),1))
        maskfiner3D_ptr=>maskfiner3D

         ! ngrow_make_distance ghost cells
         ! in 3D:
         ! FSI_lo3D,FSI_hi3D = tilelo,tilehi
         ! FSI_growlo3D,FSI_growhi3D = 
         !   grow(FSI_lo3D,FSI_hi3D,ngrow_make_distance)
        do i=FSI_growlo3D(1),FSI_growhi3D(1)
        do j=FSI_growlo3D(2),FSI_growhi3D(2)
        do k=FSI_growlo3D(3),FSI_growhi3D(3)
         idx(1)=i
         idx(2)=j
         idx(3)=k

         if (SDIM.eq.3) then

          call gridsten_level(xsten,i,j,k,level,nhalf)
          do dir=1,6*num_materials
           if (FSI_PRESSURE_FORCE_ONLY.eq.1) then
            stressdata3D(i,j,k,dir)=drag(D_DECL(i,j,k), &
             DRAGCOMP_PSTRESS+dir)
           else
            print *,"only FSI_PRESSURE_FORCE_ONLY.eq.1 supported"
            stop
           endif
          enddo
          do dir=1,num_materials
           stressflag3D(i,j,k,dir)=drag(D_DECL(i,j,k), &
             DRAGCOMP_FLAG+dir)
          enddo
          do dir=1,SDIM
           xdata3D(i,j,k,dir)=xsten(0,dir)
          enddo
          do nc=1,2
           masknbr3D(i,j,k,nc)=masknbr(D_DECL(i,j,k),nc)
          enddo
          maskfiner3D(i,j,k,1)=maskfiner(D_DECL(i,j,k),1)

         else if (SDIM.eq.2) then

          k2d=0
          ! dir is the coordinant on the 3D grid
          ! xmap3D(dir) is the coordinant on the 2D grid
          ! idx(1)=i  idx(2)=j   idx(3) unused
          do dir=1,3
           if (xmap3D(dir).eq.0) then
            ! do nothing
           else if (xmap3D(dir).eq.1) then
            i2d=idx(dir)
           else if (xmap3D(dir).eq.2) then
            j2d=idx(dir)
           else 
            print *,"xmap3D invalid"
            stop
           endif
          enddo ! dir=1..3

          if ((i2d.lt.tilelo(1)-ngrow_make_distance).or. &
              (i2d.gt.tilehi(1)+ngrow_make_distance)) then
           print *,"i2d out of range"
           stop
          endif
          if ((j2d.lt.tilelo(2)-ngrow_make_distance).or. &
              (j2d.gt.tilehi(2)+ngrow_make_distance)) then
           print *,"j2d out of range"
           stop
          endif
          call gridsten_level(xsten,i2d,j2d,k2d,level,nhalf)

          do im_local=1,num_materials
           ibase=6*(im_local-1)
           do istress=1,3
           do jstress=1,3
            stress_2d(istress,jstress)=zero
           enddo
           enddo
           do dir=1,6
            call stress_index(dir,istress,jstress)
            if (FSI_PRESSURE_FORCE_ONLY.eq.1) then
             stress_2d(istress,jstress)= &
               drag(D_DECL(i2d,j2d,k2d),DRAGCOMP_PSTRESS+ibase+dir)
            else
             print *,"only FSI_PRESSURE_FORCE_ONLY.eq.1 supported"
             stop
            endif
           enddo
           stress_2d(2,1)=stress_2d(1,2)
           stress_2d(3,1)=stress_2d(1,3)
           stress_2d(3,2)=stress_2d(2,3)
            ! istress,jstress correspond to 3D stress
            ! xmap3D(istress),xmap3D(jstress) correspond to 2D stress
           do istress=1,3
            if (xmap3D(istress).eq.0) then
             do jstress=1,3
              stress_3d(istress,jstress)=zero
             enddo
            else if ((xmap3D(istress).eq.1).or. &
                     (xmap3D(istress).eq.2)) then
             do jstress=1,3
              if (xmap3D(jstress).eq.0) then
               stress_3d(istress,jstress)=zero
              else if ((xmap3D(jstress).eq.1).or. &
                       (xmap3D(jstress).eq.2)) then
               stress_3d(istress,jstress)= &
                   stress_2d(xmap3D(istress),xmap3D(jstress))
              else
               print *,"xmap3D(jstress) invalid"
               stop
              endif
             enddo ! jstress=1,3
            else
             print *,"xmap3D(istress) invalid"
             stop
            endif
           enddo ! istress=1,3

           do dir=1,6
            call stress_index(dir,istress,jstress)
            if (abs(stress_3d(istress,jstress)).le.1.0D+20) then
             stressdata3D(i,j,k,ibase+dir)= &
               stress_3d(istress,jstress)
            else
             print *,"stress_3d overflow SOLIDFLUID.F90 1093"
             stop
            endif
           enddo

           stressflag3D(i,j,k,im_local)= &
            drag(D_DECL(i2d,j2d,k2d),DRAGCOMP_FLAG+im_local)

          enddo ! im_local=1..num_materials

          do nc=1,2
           masknbr3D(i,j,k,nc)=masknbr(D_DECL(i2d,j2d,k2d),nc)
          enddo
          maskfiner3D(i,j,k,1)=maskfiner(D_DECL(i2d,j2d,k2d),1)
          do dir=1,3
           if (xmap3D(dir).eq.0) then
            xdata3D(i,j,k,dir)=xslice3D(dir)+idx(dir)*dx_max_level(1)
           else if (xmap3D(dir).eq.1) then 
            xdata3D(i,j,k,dir)=xsten(0,1)
           else if (xmap3D(dir).eq.2) then 
            xdata3D(i,j,k,dir)=xsten(0,2)
           else 
            print *,"xmap3D invalid"
            stop
           endif
          enddo ! dir=1..3
         else
          print *,"dimension bust"
          stop
         endif 
        enddo !k
        enddo !j
        enddo !i

        do partid=1,nparts
         im_part=im_solid_map(partid)+1
         if ((im_part.lt.1).or.(im_part.gt.num_materials)) then
          print *,"im_part invalid fort_headermsg"
          stop
         endif
         if (is_lag_part(im_part).eq.1) then

          if ((FSI_flag(im_part).eq.FSI_PRESCRIBED_NODES).or. & 
              (FSI_flag(im_part).eq.FSI_SHOELE_CTML).or. & 
              (FSI_flag(im_part).eq.FSI_ICE_NODES_INIT).or. & 
              (FSI_flag(im_part).eq.FSI_FLUID_NODES_INIT)) then 

           if (container_allocated.ne.1) then
            print *,"container_allocated.ne.1"
            stop
           endif
           if (level_container_allocated(lev77).ne.1) then
            print *,"level_container_allocated(lev77).ne.1"
            stop
           endif
           if (tilenum+1.gt.contain_elem(lev77)% &
               num_tiles_on_thread3D_proc(tid+1)) then
            print *,"tilenum+1.gt.num_tiles_on_thread3D_proc(tid+1)"
            stop
           endif
           if (contain_elem(lev77)%gridno3D(tid+1,tilenum+1).ne.gridno) then
            print *,"gridno3D(tid+1,tilenum+1).ne.gridno (1)"
            print *,"gridno= ",gridno
            print *,"lev77= ",lev77
            print *,"tid= ",tid
            print *,"contain_elem(lev77)%num_grids_on_level ", &
             contain_elem(lev77)%num_grids_on_level
            print *,"contain_elem(lev77)%num_grids_on_level_proc ", &
             contain_elem(lev77)%num_grids_on_level_proc
            print *, &
             "contain_elem(lev77)%max_num_tiles_on_thread3D_proc ", &
             contain_elem(lev77)%max_num_tiles_on_thread3D_proc
            print *, &
             "contain_elem(lev77)%num_tiles_on_thread3D_proc(tid+1) ",&
             contain_elem(lev77)%num_tiles_on_thread3D_proc(tid+1)
            print *,"contain_elem(lev77)%gridno3D(tid+1,tilenum+1)=", &
             contain_elem(lev77)%gridno3D(tid+1,tilenum+1)
            stop
           endif
           do dir=1,3
            if (contain_elem(lev77)%tilelo3D(tid+1,tilenum+1,dir).ne. &
                FSI_lo3D(dir)) then
             print *,"tilelo3D(tid+1,tilenum+1,dir).ne.FSI_lo3D(dir)"
             stop
            endif 
            if (contain_elem(lev77)%tilehi3D(tid+1,tilenum+1,dir).ne. &
                FSI_hi3D(dir)) then
             print *,"tilehi3D(tid+1,tilenum+1,dir).ne.FSI_hi3D(dir)"
             stop
            endif 
           enddo ! dir=1..3
            ! inode=1...num_nodes
            ! copies to FSI(part_id)%NodeForce(dir,inode)
           call CLSVOF_Copy_To_LAG( &
            SDIM, &
            lev77, &
            tid, &
            tilenum, &
            im_part, &
            nparts, &
            partid, &
            ngrow_make_distance, &
            nFSI, &
            FSI_operation, &
            cur_time, &
            problo3D,probhi3D, &
            xmap3D, &
            dx3D, &
            xlo3D_tile, &
            xhi3D_tile, &
            FSI_lo3D,FSI_hi3D, &
            FSI_growlo3D,FSI_growhi3D, &
            growlo3D,growhi3D, &
            xdata3D_ptr, &
            stressdata3D_ptr, &
            stressflag3D_ptr, &
            masknbr3D_ptr, &
            maskfiner3D_ptr, &
            ioproc,isout)

          else if (FSI_flag(im_part).eq.FSI_PRESCRIBED_PROBF90) then 
 
           ! do nothing

          else
           print *,"FSI_flag invalid in fort_headermsg"
           print *,"im_part,FSI_flag(im_part) ", &
             im_part,FSI_flag(im_part)
           stop
          endif

         else
          print *,"is_lag_part invalid"
          stop
         endif

        enddo ! partid=1..nparts

        if (nparts.gt.num_materials) then
         print *,"nparts out of range fort_headermsg"
         stop
        endif

        deallocate(xdata3D)
        deallocate(stressdata3D)
        deallocate(stressflag3D)
        deallocate(masknbr3D)
        deallocate(maskfiner3D)

        ! called on all of the processors
       else if (FSI_sub_operation.eq.SUB_OP_FSI_SYNC_LAG_DATA) then 
        call CLSVOF_sync_lag_data(ioproc,isout)
       else
        print *,"FSI_sub_operation invalid"
        stop
       endif

      else
       print *,"FSI_operation invalid in fort_headermsg"
       stop
      endif
 
    
      return
      end subroutine fort_headermsg

      subroutine fort_aux_tecplot(auxcomp)

      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: auxcomp
      INTEGER_T plotlo(3),plothi(3) 
      INTEGER_T lo(3),hi(3) 

      INTEGER_T ih

      character*13 newfilename !auxdata ...
      character*2 auxstr

      INTEGER_T i,j,k,dir2
      INTEGER_T nwrite

! Guibo

      character*80 Title,Varname,Zonename
      REAL*4 ZONEMARKER,EOHMARKER
      integer*4 :: iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb
      INTEGER_T strandid

      ! define zone structure
      type zone_t
         real*8, pointer :: var(:,:,:,:)
      end type zone_t
      type(zone_t), dimension(:), allocatable :: zone_gb

! Guibo

      write(auxstr,'(I2)') auxcomp
      do i=1,2
       if (auxstr(i:i).eq.' ') then
        auxstr(i:i)='0'
       endif
      enddo

      write(newfilename,'(A7,A2,A4)') 'auxdata',auxstr,'.plt'

      do dir2=1,3
       plotlo(dir2)=contain_aux(auxcomp)%lo3D(dir2)-ngrow_make_distance
       plothi(dir2)=contain_aux(auxcomp)%hi3D(dir2)+ngrow_make_distance
      enddo
      call checkbound3D_array( &
        contain_aux(auxcomp)%lo3D, &
        contain_aux(auxcomp)%hi3D, &
        contain_aux(auxcomp)%LS3D, &
        ngrow_make_distance,-1)

      nwrite=3+1

      print *,"auxdata: ",newfilename

      !--------------------------------------------------
      ! Determine nzones_gb and allocate zone_gb, lo_gb, hi_gb

      allocate(zone_gb(1))
      allocate(lo_gb(1,3))
      allocate(hi_gb(1,3))

      ! Determine lo_gb, hi_gb
      ! Allocate zone_gb%var later
      iz_gb=1
      do dir2=1,3
       lo_gb(iz_gb,dir2)=plotlo(dir2)
       hi_gb(iz_gb,dir2)=plothi(dir2)
      enddo

      !-----------------------------------------------------------
      ZONEMARKER = 299.0
      EOHMARKER  = 357.0 
       !fabdata ...
      open(unit=11,file=newfilename,form="unformatted",access="stream")

      ! +++++++ HEADER SECTION ++++++

      ! i.  Magic number, Version number
      write(11) "#!TDV112"
 
      ! ii. Integer value of 1.
      write(11) 1

      ! iii. Title and variable names.
      ! File type 0 = FULL,1 = GRID,2 = SOLUTION
      write(11) 0
      ! The TITLE
      Title = "AUX LS data"
      call dumpstring(Title)
      ! Number of variables 
      write(11) nwrite

      ! Variable names.
      Varname='X'
      call dumpstring(Varname)
      Varname='Y'
      call dumpstring(Varname)
      Varname='Z'
      call dumpstring(Varname)

      ih=1
      Varname='L'
      ih=ih+1
      Varname(ih:ih)='S'
      call dumpstring(Varname)


      ! Zones
      iz_gb=1
       ! Zone marker. Value = 299.0
      write(11) ZONEMARKER
       ! Zone name
      Zonename = "ZONE"
      call dumpstring(Zonename)

      strandid=1

      write(11) -1   ! Parent Zone
      write(11) strandid-1    ! StrandID (this does not work)
      write(11) 0.0d0 ! Solution time
      write(11) -1   ! Not used. Set to -1
      write(11) 0    ! Zone Type
      write(11) 0    ! Specify Var Location. 0 = Don't specify, 
                     ! all data is located at the nodes.
      write(11) 0    ! Are raw local 1-to-1 face neighbors supplied?
      write(11) 0    ! Number of miscellaneous user-defined  
                      !face neighbor connections

        ! ----- IMax,JMax,KMax
      do dir2=1,3
       write(11) hi_gb(iz_gb,dir2)-lo_gb(iz_gb,dir2)+1
      enddo
 
      write(11) 0

      write(11) EOHMARKER
 
      ! +++++++ DATA SECTION ++++++

      iz_gb=1

      do dir2=1,3
       lo(dir2)=lo_gb(iz_gb,dir2)
       hi(dir2)=hi_gb(iz_gb,dir2)
      enddo

      allocate(zone_gb(iz_gb)% &
        var(nwrite,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

       ! order is IMPORTANT.
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

       do dir2=1,3
        zone_gb(iz_gb)%var(dir2,i,j,k)=aux_xdata3D(i,j,k,dir2)
       enddo

       zone_gb(iz_gb)%var(3+1,i,j,k)=contain_aux(auxcomp)%LS3D(i,j,k,1)
      enddo
      enddo
      enddo


       ! Zone marker  Value = 299.0
      write(11) ZONEMARKER
       ! Data format
      do i=1,nwrite
        write(11) 2
      enddo
      write(11) 0  ! Has passive variables: 0 = no, 1 = yes.
      write(11) 0  ! Has variable sharing 0 = no, 1 = yes.
      write(11) -1 ! Share connectivity list (-1 = no sharing). 
  
      do ivar_gb=1,nwrite
       write(11) minval(zone_gb(1)%var(ivar_gb,:,:,:))
       write(11) maxval(zone_gb(1)%var(ivar_gb,:,:,:))
      enddo

       ! order is IMPORTANT
      do ivar_gb=1,nwrite
       do k=lo_gb(iz_gb,3),hi_gb(iz_gb,3)
       do j=lo_gb(iz_gb,2),hi_gb(iz_gb,2)
       do i=lo_gb(iz_gb,1),hi_gb(iz_gb,1)
        write(11) zone_gb(iz_gb)%var(ivar_gb,i,j,k)
       enddo
       enddo
       enddo
      enddo

      deallocate(zone_gb(iz_gb)%var)

      deallocate(zone_gb)
      deallocate(lo_gb)
      deallocate(hi_gb)

      close(11)
      return
      end subroutine fort_aux_tecplot


      subroutine fort_aux_tecplot_full(auxcomp,FSI_operation,iter)

      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: auxcomp
      INTEGER_T, INTENT(in) :: FSI_operation
      INTEGER_T, INTENT(in) :: iter
      INTEGER_T plotlo(3),plothi(3) 
      INTEGER_T lo(3),hi(3) 

      INTEGER_T ih

      character*11 newfilename !auxfull.plt

      INTEGER_T i,j,k,dir2
      INTEGER_T nwrite

! Guibo

      character*80 Title,Varname,Zonename
      REAL*4 ZONEMARKER,EOHMARKER
      integer*4 :: iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb
      INTEGER_T strandid

      ! define zone structure
      type zone_t
         real*8, pointer :: var(:,:,:,:)
      end type zone_t
      type(zone_t), dimension(:), allocatable :: zone_gb
      INTEGER_T iread

! Guibo

      write(newfilename,'(A11)') 'auxfull.plt'

      do dir2=1,3
       plotlo(dir2)=contain_aux(auxcomp)%lo3D(dir2)-ngrow_make_distance
       plothi(dir2)=contain_aux(auxcomp)%hi3D(dir2)+ngrow_make_distance
      enddo
      call checkbound3D_array( &
        contain_aux(auxcomp)%lo3D, &
        contain_aux(auxcomp)%hi3D, &
        contain_aux(auxcomp)%LS3D, &
        ngrow_make_distance,-1)

      call checkbound3D_array( &
        contain_aux(auxcomp)%lo3D, &
        contain_aux(auxcomp)%hi3D, &
        aux_xdata3D, &
        ngrow_make_distance,-1)

      call checkbound3D_array( &
        contain_aux(auxcomp)%lo3D, &
        contain_aux(auxcomp)%hi3D, &
        aux_FSIdata3D, &
        ngrow_make_distance,-1)

      call checkbound3D_array( &
        contain_aux(auxcomp)%lo3D, &
        contain_aux(auxcomp)%hi3D, &
        aux_masknbr3D, &
        ngrow_make_distance,-1)

      nwrite=3+3

      print *,"auxdata: ",newfilename
      print *,"auxcomp ",auxcomp
      print *,"FSI_operation ",FSI_operation
      print *,"iter ",iter

      !--------------------------------------------------
      ! Determine nzones_gb and allocate zone_gb, lo_gb, hi_gb

      allocate(zone_gb(1))
      allocate(lo_gb(1,3))
      allocate(hi_gb(1,3))

      ! Determine lo_gb, hi_gb
      ! Allocate zone_gb%var later
      iz_gb=1
      do dir2=1,3
       lo_gb(iz_gb,dir2)=plotlo(dir2)
       hi_gb(iz_gb,dir2)=plothi(dir2)
      enddo

      !-----------------------------------------------------------
      ZONEMARKER = 299.0
      EOHMARKER  = 357.0 
       !fabdata ...
      open(unit=11,file=newfilename,form="unformatted",access="stream")

      ! +++++++ HEADER SECTION ++++++

      ! i.  Magic number, Version number
      write(11) "#!TDV112"
 
      ! ii. Integer value of 1.
      write(11) 1

      ! iii. Title and variable names.
      ! File type 0 = FULL,1 = GRID,2 = SOLUTION
      write(11) 0
      ! The TITLE
      Title = "AUX SANITY data"
      call dumpstring(Title)
      ! Number of variables 
      write(11) nwrite

      ! Variable names.
      Varname='X'
      call dumpstring(Varname)
      Varname='Y'
      call dumpstring(Varname)
      Varname='Z'
      call dumpstring(Varname)

       ! FSI_LEVELSET
      ih=1
      Varname='L'
      ih=ih+1
      Varname(ih:ih)='S'
      call dumpstring(Varname)

       ! FSI_SIGN_CONFLICT
      ih=1
      Varname='S'
      ih=ih+1
      Varname(ih:ih)='C'
      call dumpstring(Varname)

       ! FSI_EXTRAP_FLAG
      ih=1
      Varname='X'
      ih=ih+1
      Varname(ih:ih)='T'
      call dumpstring(Varname)

      ! Zones
      iz_gb=1
       ! Zone marker. Value = 299.0
      write(11) ZONEMARKER
       ! Zone name
      Zonename = "ZONE"
      call dumpstring(Zonename)

      strandid=1

      write(11) -1   ! Parent Zone
      write(11) strandid-1    ! StrandID (this does not work)
      write(11) 0.0d0 ! Solution time
      write(11) -1   ! Not used. Set to -1
      write(11) 0    ! Zone Type
      write(11) 0    ! Specify Var Location. 0 = Don't specify, 
                     ! all data is located at the nodes.
      write(11) 0    ! Are raw local 1-to-1 face neighbors supplied?
      write(11) 0    ! Number of miscellaneous user-defined  
                      !face neighbor connections

        ! ----- IMax,JMax,KMax
      do dir2=1,3
       write(11) hi_gb(iz_gb,dir2)-lo_gb(iz_gb,dir2)+1
      enddo
 
      write(11) 0

      write(11) EOHMARKER
 
      ! +++++++ DATA SECTION ++++++

      iz_gb=1

      do dir2=1,3
       lo(dir2)=lo_gb(iz_gb,dir2)
       hi(dir2)=hi_gb(iz_gb,dir2)
      enddo

      allocate(zone_gb(iz_gb)% &
        var(nwrite,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

       ! order is IMPORTANT.
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

       do dir2=1,3
        zone_gb(iz_gb)%var(dir2,i,j,k)=aux_xdata3D(i,j,k,dir2)
       enddo

       zone_gb(iz_gb)%var(3+1,i,j,k)=aux_FSIdata3D(i,j,k,FSI_LEVELSET+1)
       zone_gb(iz_gb)%var(3+2,i,j,k)= &
         aux_FSIdata3D(i,j,k,FSI_SIGN_CONFLICT+1)
       zone_gb(iz_gb)%var(3+3,i,j,k)= &
         aux_FSIdata3D(i,j,k,FSI_EXTRAP_FLAG+1)
      enddo
      enddo
      enddo


       ! Zone marker  Value = 299.0
      write(11) ZONEMARKER
       ! Data format
      do i=1,nwrite
        write(11) 2
      enddo
      write(11) 0  ! Has passive variables: 0 = no, 1 = yes.
      write(11) 0  ! Has variable sharing 0 = no, 1 = yes.
      write(11) -1 ! Share connectivity list (-1 = no sharing). 
  
      do ivar_gb=1,nwrite
       write(11) minval(zone_gb(1)%var(ivar_gb,:,:,:))
       write(11) maxval(zone_gb(1)%var(ivar_gb,:,:,:))
      enddo

       ! order is IMPORTANT
      do ivar_gb=1,nwrite
       do k=lo_gb(iz_gb,3),hi_gb(iz_gb,3)
       do j=lo_gb(iz_gb,2),hi_gb(iz_gb,2)
       do i=lo_gb(iz_gb,1),hi_gb(iz_gb,1)
        write(11) zone_gb(iz_gb)%var(ivar_gb,i,j,k)
       enddo
       enddo
       enddo
      enddo

      deallocate(zone_gb(iz_gb)%var)

      deallocate(zone_gb)
      deallocate(lo_gb)
      deallocate(hi_gb)

      close(11)

      print *,"enter a digit then return"
      read(*,*) iread

      return
      end subroutine fort_aux_tecplot_full


       ! called from: NavierStokes::init_aux_data()
       ! init_aux_data() called from NavierStokes::post_restart, and
       ! init_aux_data() called from NavierStokes::initData
      subroutine fort_init_aux_data(ioproc) &
      bind(c,name='fort_init_aux_data')
      use CLSVOFCouplerIO, only : CLSVOF_Read_aux_Header, &
       CLSVOF_Init_aux_Box
      use global_utility_module
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: ioproc
      INTEGER_T auxcomp
      INTEGER_T FSI_operation
      INTEGER_T FSI_touch_flag
      INTEGER_T aux_isout
      INTEGER_T iter
      INTEGER_T dir
      INTEGER_T i,j,k
      INTEGER_T LSLO(3),LSHI(3)

      if (aux_data_allocated.eq.0) then

       do auxcomp=1,fort_num_local_aux_grids

        aux_isout=1

         ! CLSVOF_Read_aux_Header is declared in: sci_clsvof.F90
         !  aux_masknbr3D,aux_FSIdata3D, and aux_xdata3D are 
         !  allocated and initialized in this routine. 
        call CLSVOF_Read_aux_Header(auxcomp,ioproc,aux_isout)

        FSI_operation=OP_FSI_MAKE_DISTANCE
        iter=0
        FSI_touch_flag=0
        call CLSVOF_Init_aux_Box(FSI_operation,iter,auxcomp, &
          FSI_touch_flag,ioproc,aux_isout)

        if (1.eq.0) then
         call fort_aux_tecplot_full(auxcomp,FSI_operation,iter)
        endif

        do while (FSI_touch_flag.eq.1)
         FSI_operation=OP_FSI_MAKE_SIGN
         FSI_touch_flag=0
         call CLSVOF_Init_aux_Box(FSI_operation,iter,auxcomp, &
          FSI_touch_flag,ioproc,aux_isout)
         iter=iter+1

         if (1.eq.0) then
          call fort_aux_tecplot_full(auxcomp,FSI_operation,iter)
         endif

        enddo !do while (FSI_touch_flag.eq.1)

        do dir=1,3
         LSLO(dir)=contain_aux(auxcomp)%lo3D(dir)-ngrow_make_distance
         LSHI(dir)=contain_aux(auxcomp)%hi3D(dir)+ngrow_make_distance
        enddo
        do i=LSLO(1),LSHI(1)
        do j=LSLO(2),LSHI(2)
        do k=LSLO(3),LSHI(3)
         contain_aux(auxcomp)%LS3D(i,j,k,1)= &
            aux_FSIdata3D(i,j,k,FSI_LEVELSET+1)
        enddo
        enddo
        enddo
       
        if (ioproc.eq.1) then
         call fort_aux_tecplot(auxcomp)
        else if (ioproc.eq.0) then
         ! do nothing
        else
         print *,"ioproc invalid"
         stop
        endif
  
        deallocate(aux_xdata3D) 
        deallocate(aux_FSIdata3D) 
        deallocate(aux_masknbr3D) 

       enddo ! auxcomp=1,fort_num_local_aux_grids
 
       aux_data_allocated=1

      else if (aux_data_allocated.eq.1) then
       ! do nothing
      else
       print *,"aux_data_allocated invalid"
       stop
      endif

      return
      end subroutine fort_init_aux_data

      subroutine fort_fillcontainer( &
        level, &
        finest_level, &
        sci_max_level, &
        cur_time, &
        dt, &
        tilelo_array, &
        tilehi_array, &
        xlo_array, &
        dx, &
        dx_max_level, &
        num_grids_on_level, & !total number of grids on "level"
        num_grids_on_level_proc, & ! number of grids on "level" and on
                                   ! this (mpi) proc.
        gridno_array, &
        num_tiles_on_thread_proc, &
        nthread_parm, &
        max_num_tiles_on_thread_proc, &
        tile_dim, & ! nthreads x max_num_tiles_on_thread_proc
        nparts, &
        im_solid_map) &
      bind(c,name='fort_fillcontainer')

      use CLSVOFCouplerIO, only : CLSVOF_FILLCONTAINER
      use solidfluid_module
      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: level
      INTEGER_T, INTENT(in) :: finest_level
      INTEGER_T, INTENT(in) :: sci_max_level
      INTEGER_T, INTENT(in) :: nparts
      INTEGER_T, INTENT(in) :: im_solid_map(nparts)
      INTEGER_T, INTENT(in) :: nthread_parm
      INTEGER_T, INTENT(in) :: num_grids_on_level
      INTEGER_T, INTENT(in) :: num_grids_on_level_proc
      INTEGER_T, INTENT(in) :: max_num_tiles_on_thread_proc
      INTEGER_T, INTENT(in) :: tile_dim
      INTEGER_T, INTENT(in) :: tilelo_array(tile_dim*SDIM)
      INTEGER_T, INTENT(in) :: tilehi_array(tile_dim*SDIM)
      REAL_T, INTENT(in) :: cur_time,dt
      REAL_T, INTENT(in) :: xlo_array(tile_dim*SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: dx_max_level(SDIM)
      INTEGER_T, INTENT(in) :: gridno_array(tile_dim)
      INTEGER_T, INTENT(in) :: num_tiles_on_thread_proc(nthread_parm)
   
      REAL_T problo3D(3),probhi3D(3)
      REAL_T dx3D(3)
      INTEGER_T xmap3D(3)
      REAL_T xslice3D(3)
      INTEGER_T dir
      INTEGER_T ilev,tid,partid,tilenum
      INTEGER_T icomp
      INTEGER_T lo3D,hi3D
      REAL_T xlo3D
      INTEGER_T im_part
      INTEGER_T lev77
      INTEGER_T local_nelems
      INTEGER_T local_nnodes

      if (nthread_parm.ne.geom_nthreads) then
       print *,"nthread_parm invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid in fillcontainer"
       stop
      endif
      lev77=level+1
      if ((lev77.lt.1).or.(lev77.gt.finest_level+1)) then
       print *,"lev77 invalid in fillcontainer"
       stop
      endif

      if (finest_level.gt.sci_max_level) then
       print *,"sci_max_level invalid"
       stop
      endif
      if ((nparts.lt.1).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_fillcontainer"
       stop
      endif
      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid"
       stop
      endif

      do dir=1,SDIM
       if (probhi_array(dir)-problo_array(dir).gt.zero) then
        ! do nothing
       else
        print *,"probhi_array(dir)-problo_array(dir) invalid"
        stop
       endif
      enddo

      if (num_grids_on_level.ge.1) then
       ! do nothing
      else
       print *,"num_grids_on_level invalid"
       stop
      endif

      if (num_grids_on_level_proc.eq.0) then
       if (max_num_tiles_on_thread_proc.ne.0) then
        print *,"max_num_tiles_on_thread_proc.ne.0"
        stop
       endif
      else if ((num_grids_on_level_proc.gt.0).and. &
               (num_grids_on_level_proc.le. &
                num_grids_on_level)) then
       if (tile_dim.ne.nthread_parm*max_num_tiles_on_thread_proc) then
        print *,"tile_dim invalid"
        stop
       endif
       if (num_grids_on_level_proc.gt. &
           nthread_parm*max_num_tiles_on_thread_proc) then
        print *,"num_grids_on_level_proc invalid"
        stop
       endif
      else
       print *,"num_grids_on_level_proc invalid"
       stop
      endif

      do tilenum=1,tile_dim
       if ((gridno_array(tilenum).ge.0).and. &
           (gridno_array(tilenum).le.num_grids_on_level-1)) then
        ! do nothing
       else
        print *,"gridno_array(tilenum) invalid: ",gridno_array(tilenum)
        print *,"tilenum=",tilenum
        print *,"num_grids_on_level=",num_grids_on_level
        stop
       endif
      enddo ! tilenum=1,tile_dim

      call init_3D_map(xmap3D,xslice3D,problo3D,probhi3D,dx_max_level)

      if (SDIM.eq.3) then
       do dir=1,SDIM
        dx3D(dir)=dx(dir)
       enddo
      else if (SDIM.eq.2) then
       do dir=1,3
        if ((xmap3D(dir).eq.1).or. &
            (xmap3D(dir).eq.2)) then
         dx3D(dir)=dx(xmap3D(dir))
        else if (xmap3D(dir).eq.0) then
         dx3D(dir)=dx_max_level(1)
        else
         print *,"xmap3D invalid"
         stop
        endif
       enddo ! dir=1..3
      else
       print *,"dimension bust"
       stop
      endif

      if (container_allocated.eq.1) then
       if (level_container_allocated(lev77).eq.1) then

        if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.gt.0) then
          ! number_threads * max_number_tiles_per_thread * 3
         deallocate(contain_elem(lev77)%tilelo3D)
         deallocate(contain_elem(lev77)%tilehi3D)
         deallocate(contain_elem(lev77)%xlo3D)
         ! number_threads * max_number_tiles_per_thread 
         deallocate(contain_elem(lev77)%gridno3D)
         do tid=1,nthread_parm
          do partid=1,nparts
           do tilenum=1,contain_elem(lev77)% &
                       num_tiles_on_thread3D_proc(tid)

            local_nelems=contain_elem(lev77)% &
                         level_elem_data(tid,partid,tilenum)%numElems
            if (local_nelems.eq.0) then
             ! do nothing
            else if (local_nelems.gt.0) then
             deallocate(contain_elem(lev77)% &
                        level_elem_data(tid,partid,tilenum)%  &
                        ElemData)
            else
             print *,"local_nelems invalid"
             stop
            endif

            local_nnodes=contain_elem(lev77)% &
                         level_node_data(tid,partid,tilenum)%numNodes
            if (local_nnodes.eq.0) then
             ! do nothing
            else if (local_nnodes.gt.0) then
             deallocate(contain_elem(lev77)% &
                        level_node_data(tid,partid,tilenum)%  &
                        NodeData)
            else
             print *,"local_nnodes invalid"
             stop
            endif
           enddo ! tilenum
          enddo ! partid

          ! number_threads
          deallocate(contain_elem(lev77)%num_tiles_on_thread3D_proc)
          ! number_threads,nparts,max_number_tiles_per_thread_proc
          deallocate(contain_elem(lev77)%level_elem_data)
          deallocate(contain_elem(lev77)%level_node_data)
         enddo ! tid=1,nthread_parm
        else if (contain_elem(lev77)%max_num_tiles_on_thread3D_proc.eq.0) then
         ! do nothing
        else
         print *,"max_num_tiles_on_thread3D_proc invalid"
         stop
        endif
        level_container_allocated(lev77)=0
       else if (level_container_allocated(lev77).eq.0) then
        ! do nothing
       else
        print *,"level_container_allocated(lev77) invalid"
        stop
       endif
      else if (container_allocated.eq.0) then
       container_allocated=1
       allocate(level_container_allocated(sci_max_level+1))
       do ilev=1,sci_max_level+1
        level_container_allocated(ilev)=0
       enddo
       allocate(contain_elem(sci_max_level+1))
      else
       print *,"container_allocated invalid"
       stop
      endif

      if (level_container_allocated(lev77).eq.0) then

       level_container_allocated(lev77)=1

       contain_elem(lev77)%num_grids_on_level=num_grids_on_level
       contain_elem(lev77)%num_grids_on_level_proc= &
        num_grids_on_level_proc

       contain_elem(lev77)%max_num_tiles_on_thread3D_proc= &
         max_num_tiles_on_thread_proc

       if (max_num_tiles_on_thread_proc.ge.1) then
        allocate(contain_elem(lev77)% &
         tilelo3D(nthread_parm,max_num_tiles_on_thread_proc,3))
        allocate(contain_elem(lev77)% &
         tilehi3D(nthread_parm,max_num_tiles_on_thread_proc,3))
        allocate(contain_elem(lev77)% &
         xlo3D(nthread_parm,max_num_tiles_on_thread_proc,3))
        allocate(contain_elem(lev77)% &
         gridno3D(nthread_parm,max_num_tiles_on_thread_proc))
        allocate(contain_elem(lev77)% &
         num_tiles_on_thread3D_proc(nthread_parm))
        allocate(contain_elem(lev77)% &
         level_elem_data(nthread_parm,nparts,max_num_tiles_on_thread_proc))
        allocate(contain_elem(lev77)% &
         level_node_data(nthread_parm,nparts,max_num_tiles_on_thread_proc))

        do tid=1,nthread_parm
         contain_elem(lev77)%num_tiles_on_thread3D_proc(tid)= &
          num_tiles_on_thread_proc(tid)
         if (num_tiles_on_thread_proc(tid).gt. &
             max_num_tiles_on_thread_proc) then
          print *,"num_tiles_on_thread_proc(tid) invalid"
          stop
         endif
         do tilenum=1,num_tiles_on_thread_proc(tid)
          do dir=1,3
           if (SDIM.eq.3) then
            icomp=max_num_tiles_on_thread_proc*SDIM*(tid-1)+ &
             SDIM*(tilenum-1)+dir
            lo3D=tilelo_array(icomp)
            hi3D=tilehi_array(icomp)
            xlo3D=xlo_array(icomp)
           else if (SDIM.eq.2) then
            icomp=max_num_tiles_on_thread_proc*SDIM*(tid-1)+ &
             SDIM*(tilenum-1)+xmap3D(dir)
            if ((xmap3D(dir).eq.1).or. &
                (Xmap3D(dir).eq.2)) then
             lo3D=tilelo_array(icomp)
             hi3D=tilehi_array(icomp)
             xlo3D=xlo_array(icomp)
            else if (xmap3D(dir).eq.0) then
             lo3D=0
             hi3D=0
             xlo3D=xslice3D(dir)-half*dx3D(dir)
            else
             print *,"xmap3D invalid"
             stop
            endif
           else
            print *,"dimension bust"
            stop
           endif
           contain_elem(lev77)%tilelo3D(tid,tilenum,dir)=lo3D
           contain_elem(lev77)%tilehi3D(tid,tilenum,dir)=hi3D
           contain_elem(lev77)%xlo3D(tid,tilenum,dir)=xlo3D
          enddo ! dir=1..3
          icomp=max_num_tiles_on_thread_proc*(tid-1)+tilenum
          contain_elem(lev77)%gridno3D(tid,tilenum)=gridno_array(icomp)
         enddo ! tilenum=1..num_tiles_on_thread_proc(tid)
        enddo ! tid=1..nthread_parm

        do partid=1,nparts
         im_part=im_solid_map(partid)+1
         if ((im_part.lt.1).or.(im_part.gt.num_materials)) then
          print *,"im_part invalid"
          stop
         endif

         if ((FSI_flag(im_part).eq.FSI_PRESCRIBED_NODES).or. & 
             (FSI_flag(im_part).eq.FSI_SHOELE_CTML).or. & 
             (FSI_flag(im_part).eq.FSI_ICE_NODES_INIT).or. & 
             (FSI_flag(im_part).eq.FSI_FLUID_NODES_INIT)) then 
          call CLSVOF_FILLCONTAINER(lev77,sci_max_level,nthread_parm, &
           dx3D,partid,im_part,cur_time,dt)
         else if (FSI_flag(im_part).eq.FSI_PRESCRIBED_PROBF90) then 
          ! do nothing
         else
          print *,"FSI_flag(im_part) invalid in fort_fillcontainer"
          print *,"im_part,FSI_flag(im_part) ",im_part,FSI_flag(im_part)
          stop
         endif
        enddo ! partid=1..nparts
       else if (max_num_tiles_on_thread_proc.eq.0) then
        ! do nothing
       else
        print *,"max_num_tiles_on_thread_proc invalid"
        stop
       endif
      else 
       print *,"level_container_allocated(lev77) invalid"
       stop
      endif

      return
      end subroutine fort_fillcontainer

      end module solidfluid_cpp_module

