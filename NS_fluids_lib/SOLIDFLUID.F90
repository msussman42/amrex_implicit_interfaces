#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 0

#if (AMREX_SPACEDIM == 2)
#define SDIM 2
#endif
#if (AMREX_SPACEDIM == 3)
#define SDIM 3
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "SOLIDFLUID_F.H"

      module solidfluid_module

      contains

      subroutine init_3D_map(xmap3D,xslice3D,problo3D,probhi3D, &
         problo,probhi,dx_max_level,probtype_in,num_materials_in)
      use global_utility_module
#if (STANDALONE==0)
      use CAV3D_module
#endif

      IMPLICIT NONE

      INTEGER_T, intent(in) :: num_materials_in
      INTEGER_T, intent(in) :: probtype_in
      INTEGER_T, intent(out) :: xmap3D(3)
      REAL_T, intent(out) :: xslice3D(3)
      REAL_T, intent(out) :: problo3D(3),probhi3D(3)
      REAL_T, intent(in) :: problo(SDIM),probhi(SDIM)
      REAL_T, intent(in) :: dx_max_level(SDIM)
      INTEGER_T dir
      INTEGER_T nmat

      nmat=num_materials_in

      do dir=1,SDIM
       if (probhi(dir)-problo(dir).le.zero) then
        print *,"probhi(dir)-problo(dir).le.zero"
        stop
       endif
       if (dx_max_level(dir).le.zero) then
        print *,"dx_max_level(dir).le.zero"
        stop
       endif
      enddo

      do dir=1,3
       xmap3D(dir)=dir
       xslice3D(dir)=zero
      enddo

      if (SDIM.eq.2) then

       if (CTML_FSI_flagF(nmat).eq.1) then ! FSI_flag==4 or 8
        xmap3D(1)=1
        xmap3D(2)=2
        xmap3D(3)=0
        xslice3D(3)=zero
        problo3D(3)=-half*dx_max_level(1)
        probhi3D(3)=half*dx_max_level(1)
       else if (CTML_FSI_flagF(nmat).eq.0) then

         ! 537 is 6 hole injector
        if ((probtype_in.eq.538).or. &
            (probtype_in.eq.537).or. &
            (probtype_in.eq.541)) then
         xmap3D(3)=2
         xmap3D(1)=1
         xmap3D(2)=0
         xslice3D(2)=zero
         problo3D(2)=problo(1)
         probhi3D(2)=probhi(1)

          ! injector C
         if (probtype_in.eq.541) then
          if (problo(1).ne.zero) then
           print *,"problo(1).ne.zero"
           stop
          endif
          xmap3D(1)=1
          xmap3D(2)=2
          xmap3D(3)=0
          xslice3D(1)=zero
          xslice3D(2)=zero
          xslice3D(3)=zero
          problo3D(3)=-1e-3
          probhi3D(3)=1e-3
         endif

        else if (probtype_in.eq.701) then  ! flapping wing
         xmap3D(1)=1
         xmap3D(3)=2
         xmap3D(2)=0
         xslice3D(2)=0.05
         problo3D(2)=-0.1
         probhi3D(2)=0.2
        else if(probtype_in.eq.539) then ! the surface is 3D
         xmap3D(1)=1
         xmap3D(2)=2
         xmap3D(3)=0
         xslice3D(3)=0.0
         problo3D(3)=-0.014
         probhi3D(3)=0.014
        else if (probtype_in.eq.9) then ! ship wave
         xmap3D(1)=1
         xmap3D(3)=2
         xmap3D(2)=0
         xslice3D(2)=0.0
         problo3D(2)=0.0
         probhi3D(2)=0.25
        else if (probtype_in.eq.5700) then
         xmap3D(1)=1
         xmap3D(2)=2
         xmap3D(3)=0
         xslice3D(3)=0.31
         problo3D(3)=0.0
         probhi3D(3)=0.62
        else if ((probtype_in.eq.400).or. &
                 (probtype_in.eq.404)) then ! gingerbread man or Xue
         xmap3D(1)=1
         xmap3D(2)=2
         xmap3D(3)=0
         xslice3D(3)=zero
         problo3D(3)=-half*dx_max_level(1)
         probhi3D(3)=half*dx_max_level(1)
        else if (probtype_in.eq.401) then ! helix
         print *,"this geometry has no 2D analogue"
         stop
        else if (probtype_in.eq.415) then ! shock sphere
         xmap3D(1)=1
         xmap3D(2)=2
         xmap3D(3)=0
         xslice3D(3)=zero
         problo3D(3)=-half*dx_max_level(1)
         probhi3D(3)=half*dx_max_level(1)
        else if (probtype_in.eq.411) then
#if (STANDALONE==0)
         call CAV3D_SLICE(xmap3D,xslice3D,problo3D,probhi3D, &
                          dx_max_level(1))
#else
         print *,"this option not for standalone version"
         stop
#endif
        endif

       else
        print *,"CTML_FSI_flagF invalid"
        stop
       endif 

      else if (SDIM.eq.3) then

       do dir=1,SDIM
        problo3D(dir)=problo(dir)
        probhi3D(dir)=probhi(dir)
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

       !FSI_operation=0  initialize node locations; generate_new_triangles
       !FSI_operation=1  update node locations
       !FSI_operation=2  make distance in narrow band
       !FSI_operation=3  update the sign.
       !FSI_operation=4  copy eul fluid vel/pres to solid
       !  FSI_sub_operation.eq.0 (clear lagrangian data)
       !  FSI_sub_operation.eq.1 (actual copy)
       !  FSI_sub_operation.eq.2 (sync lag data)
      subroutine fort_headermsg( &
        tid, &
        tilenum, &
        gridno, &
        nthread_parm, &
        level, &
        finest_level, &
        max_level, &
        FSI_operation, &
        FSI_sub_operation, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        xlo, & ! problo if FSI_operation==0
        dx, &  ! problen if FSI_operation==1
        dx_max_level, & 
        problo, &
        probhi, &
        velbc, &
        vofbc, &
        FSIdata, & ! velfab if FSI_operation==4
        DIMS(FSIdata), &
        velfab, &
        DIMS(velfab), &
        masknbr, &
        DIMS(masknbr), &
        maskfiner, &
        DIMS(maskfiner), &
        nFSI, &
        nFSI_sub, &
        ngrowFSI, &
        nparts, &
        im_solid_map, & ! type: 0..nmat-1
        h_small, &  ! smallest mesh size from the max_level.
        cur_time, &
        dt, &
        FSI_refine_factor, &
        FSI_bounding_box_ngrow, &
        touch_flag, &
        CTML_FSI_INIT, &
        CTML_force_model, &
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

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: tilenum
      INTEGER_T, intent(in) :: gridno
      INTEGER_T, intent(in) :: nthread_parm
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: max_level
      INTEGER_T, intent(in) :: FSI_operation
      INTEGER_T, intent(in) :: FSI_sub_operation
      INTEGER_T, intent(in) :: nFSI
      INTEGER_T, intent(in) :: nFSI_sub
      INTEGER_T, intent(in) :: ngrowFSI
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: im_solid_map(nparts)
      INTEGER_T nmat
      REAL_T, intent(in) :: h_small ! smallest mesh size from the max_level
      REAL_T, intent(in) :: cur_time
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: FSI_refine_factor(num_materials)
      INTEGER_T, intent(in) :: FSI_bounding_box_ngrow(num_materials)
      INTEGER_T, intent(inout) :: touch_flag
      INTEGER_T, intent(inout) :: CTML_FSI_init
      INTEGER_T, intent(in) :: CTML_force_model(num_materials)
      INTEGER_T, intent(in) :: iter
      INTEGER_T, intent(in) :: current_step 
      INTEGER_T, intent(in) :: plot_interval 
      INTEGER_T, intent(in) :: ioproc
      INTEGER_T isout
      INTEGER_T, intent(in) :: DIMDEC(FSIdata)  ! velfab if FSI_operation==4
      INTEGER_T, intent(in) :: DIMDEC(velfab) 
      INTEGER_T, intent(in) :: DIMDEC(masknbr) 
      INTEGER_T, intent(in) :: DIMDEC(maskfiner) 
        ! velfab if FSI_operation==4
      REAL_T, intent(inout), target :: FSIdata(DIMV(FSIdata),nFSI) 
      REAL_T, pointer :: FSIdata_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: velfab(DIMV(velfab),SDIM+1)
      REAL_T, pointer :: velfab_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: masknbr(DIMV(masknbr),2)
      REAL_T, pointer :: masknbr_ptr(D_DECL(:,:,:),:)
      REAL_T, intent(in), target :: maskfiner(DIMV(maskfiner),4)
      REAL_T, pointer :: maskfiner_ptr(D_DECL(:,:,:),:)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T tilelo3D(3),tilehi3D(3)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T growlo3D(3),growhi3D(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: dx_max_level(SDIM)
      REAL_T, intent(in) :: problo(SDIM),probhi(SDIM)
      REAL_T problo3D(3),probhi3D(3)
      REAL_T dx3D(3)
      REAL_T vel3D(3)
      INTEGER_T, intent(in) :: velbc(SDIM,2)
      INTEGER_T, intent(in) :: vofbc(SDIM,2)
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
      REAL_T, allocatable :: FSIdata3D(:,:,:,:)
      REAL_T, allocatable :: veldata3D(:,:,:,:)
      REAL_T, allocatable :: xdata3D(:,:,:,:)
      REAL_T, allocatable :: masknbr3D(:,:,:,:)
      REAL_T, allocatable :: maskfiner3D(:,:,:)
      INTEGER_T mask1,mask2
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T ibase
      INTEGER_T partid
      INTEGER_T nhalf
      INTEGER_T lev77
      INTEGER_T im_local
 
      nhalf=3

      FSIdata_ptr=>FSIdata

      nmat=num_materials

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
      if ((ngrowFSI.ne.3).and.(ngrowFSI.ne.0)) then
       print *,"ngrowFSI invalid"
       stop
      endif
      if (h_small.le.zero) then
       print *,"h_small invalid"
       stop
      endif
      do dir=1,SDIM
       if (dx_max_level(dir).le.zero) then
        print *,"dx_max_level(dir).le.zero"
        stop
       endif
      enddo
      if ((touch_flag.ne.0).and.(touch_flag.ne.1)) then
       print *,"touch_flag invalid"
       stop
      endif
      if ((CTML_FSI_init.ne.0).and. &
          (CTML_FSI_init.ne.1)) then
       print *,"CTML_FSI_init invalid"
       stop
      endif
      do im_local=1,nmat
       if (FSI_flag(im_local).eq.4) then
        if ((CTML_force_model(im_local).eq.0).or. &
            (CTML_force_model(im_local).eq.1)) then
         ! do nothing
        else
         print *,"CTML_force_model invalid"
         stop
        endif
       else if (FSI_flag(im_local).eq.8) then
        if (CTML_force_model(im_local).eq.2) then
         ! do nothing
        else
         print *,"CTML_force_model invalid"
         stop
        endif
       else if (fort_FSI_flag_valid(nmat,im_local).eq.1) then
        ! do nothing
       else
        print *,"fort_FSI_flag_valid invalid"
        stop
       endif
      enddo ! im_local=1..nmat

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

      if ((nparts.lt.1).or.(nparts.gt.nmat)) then
       print *,"nparts invalid fort_headermsg"
       stop
      endif
 
       ! nparts x (velocity + LS + temperature + flag + force)
      if (nFSI.ne.nparts*nFSI_sub) then 
       print *,"nFSI invalid"
       stop
      endif
      if (nFSI_sub.ne.9) then 
       print *,"nFSI_sub invalid 1 fort_headermsg: ",nFSI_sub
       stop
      endif
  
      call checkbound_array(fablo,fabhi,FSIdata_ptr,ngrowFSI,-1,2910)
      velfab_ptr=>velfab
      call checkbound_array(fablo,fabhi,velfab_ptr,ngrowFSI,-1,2910)
      masknbr_ptr=>masknbr
      call checkbound_array(fablo,fabhi,masknbr_ptr,ngrowFSI,-1,2910)
      maskfiner_ptr=>maskfiner
      call checkbound_array(fablo,fabhi,maskfiner_ptr,ngrowFSI,-1,2910)

      ! update ngrowFSI grow layers of FSIdata that do not overlap
      ! with another tile.
      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       growlo,growhi,ngrowFSI)
       ! since PCINTERP_fill_borders interpolates from coarser levels, 
       ! we only have to traverse interior values.
      call growntilebox(tilelo,tilehi,fablo,fabhi, &
       tilelo3D,tilehi3D,0)

      do dir=1,SDIM
       if (probhi(dir)-problo(dir).le.zero) then
        print *,"probhi(dir)-problo(dir).le.zero"
        stop
       endif
      enddo

      call init_3D_map(xmap3D,xslice3D,problo3D,probhi3D, &
        problo,probhi,dx_max_level,probtype,num_materials)

      if (SDIM.eq.2) then

       do dir=1,3
        if (xmap3D(dir).eq.0) then
         dx3D(dir)=dx_max_level(1)
         FSI_lo3D(dir)=0
         FSI_hi3D(dir)=0
         growlo3D(dir)=-ngrowFSI
         growhi3D(dir)=ngrowFSI
        else if ((xmap3D(dir).ge.1).and. &
                 (xmap3D(dir).le.2)) then
         dx3D(dir)=dx(xmap3D(dir))
         problo3D(dir)=problo(xmap3D(dir))
         probhi3D(dir)=probhi(xmap3D(dir))
         FSI_lo3D(dir)=tilelo(xmap3D(dir))
         FSI_hi3D(dir)=tilehi(xmap3D(dir))
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
        growlo3D(dir)=growlo(dir)
        growhi3D(dir)=growhi(dir)
       enddo ! dir=1..SDIM

      else
       print *,"dimension bust"
       stop
      endif

      do dir=1,3
       FSI_growlo3D(dir)=FSI_lo3D(dir)-ngrowFSI
       FSI_growhi3D(dir)=FSI_hi3D(dir)+ngrowFSI
      enddo

      ARG3D_L1(FSIdata3D)=FSI_growlo3D(1)
      ARG3D_L2(FSIdata3D)=FSI_growlo3D(2)
      ARG3D_L3(FSIdata3D)=FSI_growlo3D(3)
      ARG3D_H1(FSIdata3D)=FSI_growhi3D(1)
      ARG3D_H2(FSIdata3D)=FSI_growhi3D(2)
      ARG3D_H3(FSIdata3D)=FSI_growhi3D(3)

      if ((FSI_operation.eq.0).or. &  ! initialize node locations
          (FSI_operation.eq.1)) then  ! update node locations

       if ((tilenum.ne.0).or.(gridno.ne.0)) then
        print *,"tilenum or gridno invalid"
        print *,"tilenum: ",tilenum
        print *,"gridno: ",gridno
        stop
       endif

       if (FSI_sub_operation.ne.0) then
        print *,"FSI_sub_operation.ne.0"
        stop
       endif

       isout=1 ! verbose on in sci_clsvof.F90
       if (FSI_operation.eq.0) then ! initialize node locations
        call CLSVOF_ReadHeader( &
          FSI_refine_factor, &
          FSI_bounding_box_ngrow, &
          nparts, &
          im_solid_map, &
          h_small,dx_max_level,CTML_FSI_INIT, &
          cur_time,problo3D,probhi3D,ioproc,isout)
       else if (FSI_operation.eq.1) then 
        if (CTML_FSI_INIT.ne.1) then
         print *,"CTML_FSI_INIT.ne.1"
         stop
        endif
         ! cur_time=t^{n+1}
         ! if FSI_flag==4 or 8, then
         !  a) CTML_SOLVE_SOLID is called (in CTMLFSI.F90)
         !  b) tick is called (in ../Vicar3D/distFSI/tick.F)
        call CLSVOF_ReadNodes( &
          FSI_refine_factor, &
          FSI_bounding_box_ngrow, &
          cur_time,dt,h_small,problo3D,probhi3D, &
          current_step,plot_interval,ioproc,isout)
       else
        print *,"FSI_operation invalid"
        stop
       endif

      else if ((FSI_operation.eq.2).or. & ! make distance in narrow band
               (FSI_operation.eq.3)) then ! update the sign

       isout=1 ! verbose on in sci_clsvof.F90

       if (FSI_sub_operation.ne.0) then
        print *,"FSI_sub_operation.ne.0"
        stop
       endif
       if (CTML_FSI_INIT.ne.1) then
        print *,"CTML_FSI_INIT.ne.1"
        stop
       endif
       if (nFSI_sub.ne.9) then 
        print *,"nFSI_sub invalid fort headermsg 2: ",nFSI_sub
        stop
       endif
       ! nparts x (velocity + LS + temperature + flag + force)
       if (nFSI.ne.nparts*nFSI_sub) then 
        print *,"nFSI invalid"
        stop
       endif
 
       allocate(FSIdata3D(DIMV3D(FSIdata3D),nFSI))
       allocate(xdata3D(DIMV3D(FSIdata3D),3))
       allocate(masknbr3D(DIMV3D(FSIdata3D),2))

       do i=FSI_growlo3D(1),FSI_growhi3D(1)
       do j=FSI_growlo3D(2),FSI_growhi3D(2)
       do k=FSI_growlo3D(3),FSI_growhi3D(3)
        idx(1)=i
        idx(2)=j
        idx(3)=k

        if (SDIM.eq.3) then
         call gridsten_level(xsten,i,j,k,level,nhalf)
         do nc=1,nFSI
          FSIdata3D(i,j,k,nc)=FSIdata(D_DECL(i,j,k),nc)
         enddo
         do dir=1,SDIM
          xdata3D(i,j,k,dir)=xsten(0,dir)
         enddo
         do nc=1,2
          masknbr3D(i,j,k,nc)=masknbr(D_DECL(i,j,k),nc)
         enddo
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
         if ((i2d.lt.tilelo(1)-ngrowFSI).or. &
             (i2d.gt.tilehi(1)+ngrowFSI)) then
          print *,"i2d out of range"
          stop
         endif
         if ((j2d.lt.tilelo(2)-ngrowFSI).or. &
             (j2d.gt.tilehi(2)+ngrowFSI)) then
          print *,"j2d out of range"
          stop
         endif
         call gridsten_level(xsten,i2d,j2d,k2d,level,nhalf)

         do partid=1,nparts
          ibase=(partid-1)*nFSI_sub
          FSIdata3D(i,j,k,ibase+4)=FSIdata(D_DECL(i2d,j2d,k2d),ibase+4) !LS
          FSIdata3D(i,j,k,ibase+5)=FSIdata(D_DECL(i2d,j2d,k2d),ibase+5) !T
          FSIdata3D(i,j,k,ibase+6)=FSIdata(D_DECL(i2d,j2d,k2d),ibase+6) !flag
          do dir=1,3
           if (xmap3D(dir).eq.0) then
            vel3D(dir)=zero
           else if ((xmap3D(dir).eq.1).or. &
                    (xmap3D(dir).eq.2)) then
            vel3D(dir)=FSIdata(D_DECL(i2d,j2d,k2d),ibase+xmap3D(dir))
           else
            print *,"xmap3D(dir) invalid"
            stop
           endif
           FSIdata3D(i,j,k,ibase+dir)=vel3D(dir)
          enddo ! dir=1..3
          do dir=1,3
           FSIdata3D(i,j,k,ibase+6+dir)= &
             FSIdata(D_DECL(i2d,j2d,k2d),ibase+6+dir) !force
          enddo
         enddo ! partid=1,nparts
          
         do nc=1,2
          masknbr3D(i,j,k,nc)=masknbr(D_DECL(i2d,j2d,k2d),nc)
         enddo
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
        if ((im_part.lt.1).or.(im_part.gt.nmat)) then
         print *,"im_part invalid fort_headermsg"
         stop
        endif
        if (is_lag_part(nmat,im_part).eq.1) then

         if ((FSI_flag(im_part).eq.2).or. & ! prescribed solid from CAD
             (FSI_flag(im_part).eq.4).or. & ! CTML FSI
             (FSI_flag(im_part).eq.6).or. & ! ice from CAD
             (FSI_flag(im_part).eq.7)) then ! fluid from CAD

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
           ! in: sci_clsvof.F90
          call CLSVOF_InitBox( &
           iter, &
           SDIM, &
           lev77, &
           tid, &
           tilenum, &
           im_part, &
           nparts, &
           partid, &
           ngrowFSI, &
           nmat, &
           nFSI, &
           nFSI_sub, &
           FSI_operation, &
           touch_flag, &
           h_small, &
           cur_time, &
           dt, &
           problo3D,probhi3D, &
           xmap3D, &
           xslice3D, &
           dx3D, &
           FSI_lo3D,FSI_hi3D, &
           FSI_growlo3D,FSI_growhi3D, &
           growlo3D,growhi3D, &
           xdata3D, &
           FSIdata3D, &
           masknbr3D, &
           DIMS3D(FSIdata3D), &
           CTML_force_model, &
           ioproc,isout)

         else if (FSI_flag(im_part).eq.1) then ! prescribed solid (EUL)

          ! do nothing

         else
          print *,"FSI_flag invalid"
          stop
         endif

        else
         print *,"is_lag_part invalid"
         stop
        endif

       enddo ! partid=1..nparts

       if (nparts.gt.nmat-1) then
        print *,"nparts out of range"
        stop
       endif

        ! update ngrowFSI grow layers of FSIdata that do not overlap
        ! with another tile.
        ! nparts x (velocity + LS + temperature + flag)
        ! PCINTERP_fill_borders interpolates from coarser levels so
        ! that we only have to traverse interior values here.
       do i=tilelo3D(1),tilehi3D(1)
       do j=tilelo3D(2),tilehi3D(2)
       do k=tilelo3D(3),tilehi3D(3)

        mask1=NINT(masknbr(D_DECL(i,j,k),1))
        mask2=NINT(masknbr(D_DECL(i,j,k),2))

         ! mask2==1 => (i,j,k) in the interior of the tile.
         ! mask1==0 => (i,j,k) in coarse/fine ghost cell
        if ((mask1.eq.0).or.(mask2.eq.1)) then
         ! nparts x (velocity + LS + temperature + flag + force)
         if (nFSI.ne.nparts*nFSI_sub) then 
          print *,"nFSI invalid"
          stop
         endif
         if (nFSI_sub.ne.9) then 
          print *,"nFSI_sub invalid fort headermsg 3: ",nFSI_sub
          stop
         endif

         if (SDIM.eq.3) then
          do nc=1,nFSI
           FSIdata(D_DECL(i,j,k),nc)=FSIdata3D(i,j,k,nc)
          enddo
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
           if ((idx(dir).lt.FSI_growlo3D(dir)).or. &
               (idx(dir).gt.FSI_growhi3D(dir))) then
            print *,"idx(dir) out of range"
            stop
           endif
          enddo ! dir=1..3
          do partid=1,nparts
           ibase=(partid-1)*nFSI_sub
           FSIdata(D_DECL(i,j,k),ibase+4)= &
            FSIdata3D(idx(1),idx(2),idx(3),ibase+4) ! LS
           FSIdata(D_DECL(i,j,k),ibase+5)= &
            FSIdata3D(idx(1),idx(2),idx(3),ibase+5) ! T
           FSIdata(D_DECL(i,j,k),ibase+6)= &
            FSIdata3D(idx(1),idx(2),idx(3),ibase+6) ! flag
           do dir=1,3
            if (xmap3D(dir).eq.0) then
             FSIdata(D_DECL(i,j,k),ibase+3)=zero
            else if ((xmap3D(dir).eq.1).or. &
                     (xmap3D(dir).eq.2)) then
             FSIdata(D_DECL(i,j,k),ibase+xmap3D(dir))= &
              FSIdata3D(idx(1),idx(2),idx(3),ibase+dir) 
            else
             print *,"xmap3D(dir) invalid"
             stop
            endif
           enddo ! dir=1..3
           do dir=1,3
            FSIdata(D_DECL(i,j,k),ibase+6+dir)=  &
             FSIdata3D(idx(1),idx(2),idx(3),ibase+6+dir) ! force
           enddo
          enddo ! partid=1..nparts
         else
          print *,"dimension bust"
          stop
         endif 
 
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

      else if (FSI_operation.eq.4) then ! copy Eul fluid vel/pres to solid

       isout=1 ! verbose on in sci_clsvof.F90

       if (CTML_FSI_INIT.ne.1) then
        print *,"CTML_FSI_INIT.ne.1"
        stop
       endif

       if (FSI_sub_operation.eq.0) then
        call CLSVOF_clear_lag_data(ioproc,isout)
       else if (FSI_sub_operation.eq.1) then 
        allocate(veldata3D(DIMV3D(FSIdata3D),4)) ! 3+1
        allocate(xdata3D(DIMV3D(FSIdata3D),3))
        allocate(masknbr3D(DIMV3D(FSIdata3D),2))
        allocate(maskfiner3D(DIMV3D(FSIdata3D)))

        do i=FSI_growlo3D(1),FSI_growhi3D(1)
        do j=FSI_growlo3D(2),FSI_growhi3D(2)
        do k=FSI_growlo3D(3),FSI_growhi3D(3)
         idx(1)=i
         idx(2)=j
         idx(3)=k

         if (SDIM.eq.3) then
          call gridsten_level(xsten,i,j,k,level,nhalf)
          do dir=1,SDIM+1
           veldata3D(i,j,k,dir)=velfab(D_DECL(i,j,k),dir)
          enddo
          do dir=1,SDIM
           xdata3D(i,j,k,dir)=xsten(0,dir)
          enddo
          do nc=1,2
           masknbr3D(i,j,k,nc)=masknbr(D_DECL(i,j,k),nc)
          enddo
          maskfiner3D(i,j,k)=maskfiner(D_DECL(i,j,k),1)
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
          if ((i2d.lt.tilelo(1)-ngrowFSI).or. &
              (i2d.gt.tilehi(1)+ngrowFSI)) then
           print *,"i2d out of range"
           stop
          endif
          if ((j2d.lt.tilelo(2)-ngrowFSI).or. &
              (j2d.gt.tilehi(2)+ngrowFSI)) then
           print *,"j2d out of range"
           stop
          endif
          call gridsten_level(xsten,i2d,j2d,k2d,level,nhalf)

          do dir=1,3
           if (xmap3D(dir).eq.0) then
            vel3D(dir)=zero
           else if ((xmap3D(dir).eq.1).or. &
                    (xmap3D(dir).eq.2)) then
            vel3D(dir)=velfab(D_DECL(i2d,j2d,k2d),xmap3D(dir))
           else
            print *,"xmap3D(dir) invalid"
            stop
           endif
           veldata3D(i,j,k,dir)=vel3D(dir)
          enddo ! dir=1..3
          veldata3D(i,j,k,4)=velfab(D_DECL(i2d,j2d,k2d),SDIM+1)

          do nc=1,2
           masknbr3D(i,j,k,nc)=masknbr(D_DECL(i2d,j2d,k2d),nc)
          enddo
          maskfiner3D(i,j,k)=maskfiner(D_DECL(i2d,j2d,k2d),1)
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
         if ((im_part.lt.1).or.(im_part.gt.nmat)) then
          print *,"im_part invalid fort_headermsg"
          stop
         endif
         if (is_lag_part(nmat,im_part).eq.1) then

          if ((FSI_flag(im_part).eq.2).or. & ! prescribed solid CAD
              (FSI_flag(im_part).eq.4).or. & ! CTML FSI
              (FSI_flag(im_part).eq.6).or. & ! ice from CAD
              (FSI_flag(im_part).eq.7)) then ! fluid from CAD

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
            print *,"contain_elem(lev77)%max_num_tiles_on_thread3D_proc ", &
             contain_elem(lev77)%max_num_tiles_on_thread3D_proc
            print *,"contain_elem(lev77)%num_tiles_on_thread3D_proc(tid+1) ",&
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

           call CLSVOF_Copy_To_LAG( &
            SDIM, &
            lev77, &
            tid, &
            tilenum, &
            im_part, &
            nparts, &
            partid, &
            ngrowFSI, &
            nmat, &
            nFSI, &
            nFSI_sub, &
            FSI_operation, &
            cur_time, &
            problo3D,probhi3D, &
            xmap3D, &
            xslice3D, &
            dx3D, &
            FSI_lo3D,FSI_hi3D, &
            FSI_growlo3D,FSI_growhi3D, &
            growlo3D,growhi3D, &
            xdata3D, &
            veldata3D, &
            masknbr3D, &
            maskfiner3D, &
            DIMS3D(FSIdata3D), &
            ioproc,isout)

          else if (FSI_flag(im_part).eq.1) then !prescribed solid(EUL)
 
           ! do nothing

          else
           print *,"FSI_flag invalid"
           stop
          endif

         else
          print *,"is_lag_part invalid"
          stop
         endif

        enddo ! partid=1..nparts

        if (nparts.gt.nmat) then
         print *,"nparts out of range fort_headermsg"
         stop
        endif

        deallocate(xdata3D)
        deallocate(veldata3D)
        deallocate(masknbr3D)
        deallocate(maskfiner3D)

       else if (FSI_sub_operation.eq.2) then 
        call CLSVOF_sync_lag_data(ioproc,isout)
       else
        print *,"FSI_sub_operation invalid"
        stop
       endif

      else
       print *,"FSI_operation invalid"
       stop
      endif
 
    
      return
      end subroutine fort_headermsg


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
        nmat, &
        nparts, &
        im_solid_map, &
        problo, &
        probhi) &
      bind(c,name='fort_fillcontainer')

      use CLSVOFCouplerIO, only : CLSVOF_FILLCONTAINER
      use solidfluid_module
      use global_utility_module
      use probcommon_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: sci_max_level
      INTEGER_T, intent(in) :: nparts
      INTEGER_T, intent(in) :: im_solid_map(nparts)
      INTEGER_T, intent(in) :: nthread_parm
      INTEGER_T, intent(in) :: num_grids_on_level
      INTEGER_T, intent(in) :: num_grids_on_level_proc
      INTEGER_T, intent(in) :: max_num_tiles_on_thread_proc
      INTEGER_T, intent(in) :: tile_dim
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: tilelo_array(tile_dim*SDIM)
      INTEGER_T, intent(in) :: tilehi_array(tile_dim*SDIM)
      REAL_T, intent(in) :: cur_time,dt
      REAL_T, intent(in) :: xlo_array(tile_dim*SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: dx_max_level(SDIM)
      INTEGER_T, intent(in) :: gridno_array(tile_dim)
      INTEGER_T, intent(in) :: num_tiles_on_thread_proc(nthread_parm)
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: probhi(SDIM)
   
      REAL_T problo3D(3),probhi3D(3)
      REAL_T dx3D(3)
      INTEGER_T xmap3D(3)
      REAL_T xslice3D(3)
      INTEGER_T dir
      INTEGER_T ilev,tid,partid,tilenum
      INTEGER_T icomp
      INTEGER_T lo3D,hi3D
      REAL_T xlo3D
      INTEGER_T local_flag
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
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((nparts.lt.1).or.(nparts.gt.nmat)) then
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
       if (probhi(dir)-problo(dir).le.zero) then
        print *,"probhi(dir)-problo(dir).le.zero"
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

      call init_3D_map(xmap3D,xslice3D,problo3D,probhi3D, &
       problo,probhi,dx_max_level,probtype,num_materials)

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
         if ((im_part.lt.1).or.(im_part.gt.nmat)) then
          print *,"im_part invalid"
          stop
         endif
         local_flag=FSI_flag(im_part)
         if ((local_flag.eq.2).or. & !prescribed solid from CAD
             (local_flag.eq.4).or. & !CTML FSI
             (local_flag.eq.6).or. & !ice from CAD
             (local_flag.eq.7)) then !fluid from CAD
          call CLSVOF_FILLCONTAINER(lev77,sci_max_level,nthread_parm, &
           dx3D,partid,im_part,nmat,cur_time,dt)
         else if (local_flag.eq.1) then ! prescribed solid (EUL)
          ! do nothing
         else
          print *,"local_flag invalid"
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

#undef STANDALONE

