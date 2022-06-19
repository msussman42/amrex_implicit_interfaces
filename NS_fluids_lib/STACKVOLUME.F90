#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_ArrayLim.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

      module stackvolume_module
      IMPLICIT NONE

      abstract interface
        subroutine sub_interface(xsten,nhalf,dx,bfact,dist,nmat,time)
        INTEGER_T, intent(in) :: bfact
        INTEGER_T, intent(in) :: nhalf
        REAL_T, intent(in) :: dx(SDIM) 
        REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
        INTEGER_T, intent(in) :: nmat
        REAL_T, intent(in) :: time
        REAL_T, intent(out) :: dist(nmat)
        end subroutine
      end interface

      contains

        ! centroid relative to the center (or centroid?) of the cell
      subroutine extract_vof_cen_batch(fluiddata,vofdark,voflight, &
        cendark,cenlight,nmat)
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T nmat
      REAL_T fluiddata(nmat,2*SDIM+2)
      REAL_T vofdark(nmat),voflight(nmat)
      REAL_T cendark(nmat,SDIM)
      REAL_T cenlight(nmat,SDIM)
      REAL_T cenall(SDIM)
      REAL_T volall
      INTEGER_T dir,im


      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      do im=1,nmat
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
        voldata,nmat,cutflag,LS_sub,time)

      use MOF_routines_module
      use geometry_intersect_module
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      procedure(sub_interface) :: LS_sub

      REAL_T, intent(in) :: time
      INTEGER_T nmat,nten,bfact,nhalf
      REAL_T xsten(-nhalf:nhalf,SDIM)
      REAL_T dxin(SDIM)
      REAL_T voldata(nmat,2*SDIM+2)
      INTEGER_T imaterial,cutflag
      INTEGER_T i1,j1,k1,dir,minusflag,plusflag
      REAL_T ltest(D_DECL(3,3,3),nmat)
      REAL_T lnode(4*(SDIM-1),nmat)
      REAL_T facearea(nmat)
      REAL_T tempvol(nmat)
      REAL_T tempcen(nmat,SDIM)
      REAL_T volall
      REAL_T cenall(SDIM)
      REAL_T xsten2(-1:1,SDIM)
      INTEGER_T nhalf2
      INTEGER_T k1lo,k1hi,isten
      REAL_T, allocatable, dimension(:) :: distbatch

      nhalf2=1

      if (nhalf.lt.3) then
       print *,"nhalf invalid get volume data batch"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten=num_interfaces

      allocate(distbatch(nmat))

      if (bfact.lt.1) then
       print *,"bfact invalid301"
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
       call LS_sub(xsten2,nhalf2,dxin,bfact,distbatch,nmat,time)
       do imaterial=1,nmat
        ltest(D_DECL(i1+2,j1+2,k1+2),imaterial)= &
         distbatch(imaterial)
       enddo
      enddo
      enddo
      enddo ! i1,j1,k1

      cutflag=0
      do imaterial=1,nmat
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

      call data_to_node(ltest,lnode,nmat,xsten,nhalf,SDIM)  

      call fast_cell_intersection_grid_batch(bfact,dxin,xsten,nhalf, &
       lnode, &
       tempvol, &
       tempcen,facearea,volall, &
       cenall,nmat,SDIM)
      do imaterial=1,nmat
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
        voldata,nmat,level,stack_max_level,LS_sub,time)
      use probcommon_module
      IMPLICIT NONE

      procedure(sub_interface) :: LS_sub

      INTEGER_T, intent(in) :: nhalf,bfact
      INTEGER_T, intent(in) :: nmat 
      REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
      REAL_T, intent(in) :: dxin(SDIM)
      REAL_T, intent(inout) :: voldata(nmat,2*SDIM+2)
      INTEGER_T, intent(inout) :: stack_max_level
      INTEGER_T, intent(in) :: level
      REAL_T, intent(in) :: time

      INTEGER_T i1,j1,k1,k1lo,k1hi,dir,cutflag,im,isten
      REAL_T dxsten_fine(SDIM)
      REAL_T xsten_mid(SDIM)

      REAL_T, allocatable, dimension(:,:) :: localdata
      REAL_T, allocatable, dimension(:) :: dxin_fine
      REAL_T, allocatable, dimension(:,:) :: xsten_fine

      allocate(localdata(nmat,2*SDIM+2))
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
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
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
       print *,"levelrz invalid stack volume batch"
       stop
      endif

      if (level.eq.0) then
       do im=1,nmat
        do i1=1,2+2*SDIM
         voldata(im,i1)=zero
        enddo
       enddo
      endif

       ! in: stackvolume_batch
      call get_volume_data_batch(xsten,nhalf,dxin,bfact, &
        localdata,nmat,cutflag,LS_sub,time)

      if ((level.eq.stack_max_level).or.(cutflag.eq.0)) then
       do im=1,nmat
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
         voldata,nmat,level+1,stack_max_level,LS_sub,time)
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

