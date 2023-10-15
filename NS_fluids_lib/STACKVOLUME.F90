#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_ArrayLim.H"
#include "EXTRAP_COMP.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

      module stackvolume_module
      use amrex_fort_module, only : amrex_real
      IMPLICIT NONE

      abstract interface
        subroutine sub_interface(xsten,nhalf,dx,bfact,dist,time)
        use amrex_fort_module, only : amrex_real
        integer, INTENT(in) :: bfact
        integer, INTENT(in) :: nhalf
        real(amrex_real), INTENT(in) :: dx(SDIM) 
        real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
        real(amrex_real), INTENT(in) :: time
        real(amrex_real), INTENT(out) :: dist(:)
        end subroutine
      end interface

      contains

        ! centroid relative to the center (or centroid?) of the cell
      subroutine extract_vof_cen_batch(fluiddata,vofdark,voflight, &
        cendark,cenlight)
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: fluiddata(num_materials,2*SDIM+2)
      real(amrex_real), INTENT(out) :: vofdark(num_materials),voflight(num_materials)
      real(amrex_real), INTENT(out) :: cendark(num_materials,SDIM)
      real(amrex_real), INTENT(out) :: cenlight(num_materials,SDIM)
      real(amrex_real) cenall(SDIM)
      real(amrex_real) volall
      integer dir,im

      do im=1,num_materials
       volall=fluiddata(im,2)
       if (volall.le.zero) then
        print *,"volume 0 bust"
        stop
       endif
       vofdark(im)=fluiddata(im,1)/volall
       if (vofdark(im).le.VOFTOL) then
        vofdark(im)=zero
       endif
       if (vofdark(im).ge.one-VOFTOL) then
        vofdark(im)=one
       endif
       voflight(im)=one-vofdark(im)
       if (voflight(im).le.VOFTOL) then
        voflight(im)=zero
       endif
       if (voflight(im).ge.one-VOFTOL) then
        voflight(im)=one
       endif
       do dir=1,SDIM
        cenall(dir)=fluiddata(im,SDIM+2+dir)/volall
        if ((vofdark(im).gt.zero).and.(vofdark(im).lt.one)) then
         cendark(im,dir)=fluiddata(im,2+dir)/fluiddata(im,1)-cenall(dir)
        else
         cendark(im,dir)=zero
        endif
         ! c v = c_1 v_1 + c_2 v_2
         ! c_2 = (c v - c_1 v_1)/v_2

        if ((voflight(im).gt.zero).and.(voflight(im).lt.one)) then
         cenlight(im,dir)=-cendark(im,dir)*vofdark(im)/voflight(im)
        else
         cenlight(im,dir)=zero
        endif
       enddo  ! dir
      enddo !im
 
      return
      end subroutine extract_vof_cen_batch

      ! voldark,volall,cendark,cenall
      ! called from: stackvolume_batch
      subroutine get_volume_data_batch(xsten,nhalf,dxin,bfact, &
        voldata,cutflag,LS_sub,time)

      use MOF_routines_module
      use geometry_intersect_module
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      procedure(sub_interface) :: LS_sub

      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: bfact,nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: dxin(SDIM)
      real(amrex_real), INTENT(out) :: voldata(num_materials,2*SDIM+2)
      integer imaterial
      integer, INTENT(out) :: cutflag
      integer i1,j1,k1,dir,minusflag,plusflag
      real(amrex_real) ltest(D_DECL(3,3,3),num_materials)
      real(amrex_real) lnode(4*(SDIM-1),num_materials)
      real(amrex_real) facearea(num_materials)
      real(amrex_real) tempvol(num_materials)
      real(amrex_real) tempcen(num_materials,SDIM)
      real(amrex_real) volall
      real(amrex_real) cenall(SDIM)
      real(amrex_real) xsten2(-1:1,SDIM)
      integer nhalf2
      integer k1lo,k1hi,isten
      real(amrex_real), allocatable, dimension(:) :: distbatch

      nhalf2=1

      if (nhalf.lt.3) then
       print *,"nhalf invalid get volume data batch"
       stop
      endif

      allocate(distbatch(num_materials))

      if (bfact.lt.1) then
       print *,"bfact invalid301"
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
       print *,"levelrz invalid get volume data batch"
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
      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi
       do isten=-1,1
        dir=1
        xsten2(isten,dir)=xsten(isten+2*i1,dir)
        dir=2
        xsten2(isten,dir)=xsten(isten+2*j1,dir)
        if (SDIM.eq.3) then
         dir=SDIM
         xsten2(isten,dir)=xsten(isten+2*k1,dir)
        endif
       enddo ! isten
      
        ! in: get_volume_data_batch
       call LS_sub(xsten2,nhalf2,dxin,bfact,distbatch,time)
       do imaterial=1,num_materials
        ltest(D_DECL(i1+2,j1+2,k1+2),imaterial)= &
         distbatch(imaterial)
       enddo
      enddo
      enddo
      enddo ! i1,j1,k1

      cutflag=0
      do imaterial=1,num_materials
       minusflag=0
       plusflag=0

       k1=0
       do i1=-1,1
       do j1=-1,1
       do k1=k1lo,k1hi
        if (ltest(D_DECL(i1+2,j1+2,k1+2),imaterial).le.zero) then
         minusflag=1
        endif
        if (ltest(D_DECL(i1+2,j1+2,k1+2),imaterial).ge.zero) then
         plusflag=1
        endif
        if (abs(ltest(D_DECL(i1+2,j1+2,k1+2),imaterial)).le.bfact*dxin(1)) then
         minusflag=1
         plusflag=1
        endif
       enddo
       enddo
       enddo
       if ((minusflag.eq.1).and.(plusflag.eq.1)) then
        cutflag=1
       endif      
      enddo ! imaterial

      call data_to_node(ltest,lnode,num_materials,xsten,nhalf,SDIM)  

      call fast_cell_intersection_grid_batch(bfact,dxin,xsten,nhalf, &
       lnode, &
       tempvol, &
       tempcen,facearea,volall, &
       cenall,SDIM)
      do imaterial=1,num_materials
       voldata(imaterial,1)=tempvol(imaterial)
       voldata(imaterial,2)=volall
       do dir=1,SDIM
        voldata(imaterial,2+dir)=tempcen(imaterial,dir)*tempvol(imaterial)
        voldata(imaterial,2+SDIM+dir)=cenall(dir)*volall
       enddo
      enddo 

      deallocate(distbatch)
      
      return
      end subroutine get_volume_data_batch


! voldark,vollight,volall,cendark,cenlight,cenall
! called from FORT_INITDATA
      recursive subroutine stackvolume_batch(xsten,nhalf,dxin,bfact, &
        voldata,level,stack_max_level,LS_sub,time)
      use probcommon_module
      IMPLICIT NONE

      procedure(sub_interface) :: LS_sub

      integer, INTENT(in) :: nhalf,bfact
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: dxin(SDIM)
      real(amrex_real), INTENT(inout) :: voldata(num_materials,2*SDIM+2)
      integer, INTENT(inout) :: stack_max_level
      integer, INTENT(in) :: level
      real(amrex_real), INTENT(in) :: time

      integer i1,j1,k1,k1lo,k1hi,dir,cutflag,im,isten
      real(amrex_real) dxsten_fine(SDIM)
      real(amrex_real) xsten_mid(SDIM)

      real(amrex_real), allocatable, dimension(:,:) :: localdata
      real(amrex_real), allocatable, dimension(:) :: dxin_fine
      real(amrex_real), allocatable, dimension(:,:) :: xsten_fine

      allocate(localdata(num_materials,2*SDIM+2))
      allocate(dxin_fine(SDIM))
      allocate(xsten_fine(-nhalf:nhalf,SDIM))
    
      if ((stack_max_level.ge.0).and.(stack_max_level.le.10)) then
       ! do nothing
      else
       print *,"stack_max_level out of range"
       stop
      endif 

      if (bfact.lt.1) then
       print *,"bfact invalid302"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid stackvolume batch"
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
       print *,"levelrz invalid stack volume batch"
       stop
      endif

      if (level.eq.0) then
       do im=1,num_materials
        do i1=1,2+2*SDIM
         voldata(im,i1)=zero
        enddo
       enddo
      endif

       ! in: stackvolume_batch
      call get_volume_data_batch(xsten,nhalf,dxin,bfact, &
        localdata,cutflag,LS_sub,time)

      if ((level.eq.stack_max_level).or.(cutflag.eq.0)) then
       do im=1,num_materials
        do i1=1,2+2*SDIM
         voldata(im,i1)=voldata(im,i1)+localdata(im,i1)
        enddo
       enddo
      else if ((level.ge.0).and. &
               (level.lt.stack_max_level).and. &
               (cutflag.eq.1)) then
       do dir=1,SDIM
        dxin_fine(dir)=half*dxin(dir)
        dxsten_fine(dir)=(xsten(1,dir)-xsten(-1,dir))/four
        xsten_mid(dir)=(xsten(1,dir)+xsten(-1,dir))/two
       enddo
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

       do i1=-1,1,2
       do j1=-1,1,2
       do k1=k1lo,k1hi,2
        do isten=-3,3
         dir=1
         xsten_fine(isten,dir)=xsten_mid(dir)+(isten+i1)*dxsten_fine(dir)
         dir=2
         xsten_fine(isten,dir)=xsten_mid(dir)+(isten+j1)*dxsten_fine(dir)
         if (SDIM.eq.3) then
          dir=SDIM
          xsten_fine(isten,dir)=xsten_mid(dir)+(isten+k1)*dxsten_fine(dir)
         endif
        enddo ! isten

        call stackvolume_batch( &
         xsten_fine,nhalf,dxin_fine,bfact, &
         voldata,level+1,stack_max_level,LS_sub,time)
       enddo 
       enddo 
       enddo
      else
       print *,"level or cutflag invalid"
       stop
      endif

      deallocate(localdata)
      deallocate(dxin_fine)
      deallocate(xsten_fine)
       
      return
      end subroutine stackvolume_batch

      end module stackvolume_module

