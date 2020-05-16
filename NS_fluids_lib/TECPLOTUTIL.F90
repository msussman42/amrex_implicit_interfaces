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

#include "TECPLOTUTIL_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

#if (STANDALONE==1)
      module tecplotutil_cpp_module
      contains
#endif

      subroutine FORT_TECPLOTFAB( &
       time, &
       fabdata,DIMS(fabdata), &
       growlo,growhi, &
       fablo,fabhi, &
       bfact, &
       xlo,dx, &
       dir,ncomp,interior_only,nsteps)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T dir,ncomp,interior_only,nsteps,bfact
      INTEGER_T growlo(SDIM),growhi(SDIM) 
      INTEGER_T fablo(SDIM),fabhi(SDIM) 
      INTEGER_T plotlo(SDIM),plothi(SDIM) 
      INTEGER_T lo(SDIM),hi(SDIM) 
      INTEGER_T DIMDEC(fabdata)
      REAL_T fabdata(DIMV(fabdata),ncomp)
      REAL_T xlo(SDIM)
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T dx(SDIM)
      REAL_T time
      INTEGER_T klo,khi
    

      INTEGER_T ih

      character*17 newfilename
      character*2 matstr
      character*6 stepstr

      INTEGER_T i,j,k,im,nwrite,dir2

! Guibo

      character*80 Title,Varname,Zonename
      REAL*4 ZONEMARKER,EOHMARKER
      integer*4 :: iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb

      ! define zone structure
      type zone_t
         real*8, pointer :: var(:,:,:,:)
      end type zone_t
      type(zone_t), dimension(:), allocatable :: zone_gb

! Guibo

      nhalf=1

      if (bfact.lt.1) then
       print *,"bfact invalid400"
       stop
      endif

      if (ncomp.lt.1) then
       print *,"ncomp invalid 30"
       stop
      endif
      if (nsteps.lt.0) then
       print *,"nsteps invalid in tecplotfab nsteps=",nsteps
       stop
      endif
      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif
      write(stepstr,'(I6)') nsteps
      do i=1,6
       if (stepstr(i:i).eq.' ') then
        stepstr(i:i)='0'
       endif
      enddo

      write(newfilename,'(A7,A6,A4)') 'fabdata',stepstr,'.plt'

      if (interior_only.eq.0) then
       do dir2=1,SDIM
        plotlo(dir2)=growlo(dir2)
        plothi(dir2)=growhi(dir2)
       enddo
      else if (interior_only.eq.1) then
       do dir2=1,SDIM
        plotlo(dir2)=fablo(dir2)
        plothi(dir2)=fabhi(dir2)
       enddo
      else
       print *,"interior_only invalid"
       stop
      endif

      nwrite=SDIM+ncomp

      print *,"fabdata: ",newfilename

      !--------------------------------------------------
      ! Determine nzones_gb and allocate zone_gb, lo_gb, hi_gb

      allocate(zone_gb(1))
      allocate(lo_gb(1,SDIM))
      allocate(hi_gb(1,SDIM))

      ! Determine lo_gb, hi_gb
      ! Allocate zone_gb%var later
      iz_gb=1
      do dir2=1,SDIM
       lo_gb(iz_gb,dir2)=plotlo(dir2)
       hi_gb(iz_gb,dir2)=plothi(dir2)
      enddo

      !-----------------------------------------------------------
      ZONEMARKER = 299.0
      EOHMARKER  = 357.0
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
      Title = "CLSVOF data"
      call dumpstring(Title)
      ! Number of variables 
      write(11) nwrite

      ! Variable names.
      Varname='X'
      call dumpstring(Varname)
      Varname='Y'
      call dumpstring(Varname)

      if (SDIM.eq.3) then
       Varname='Z'
       call dumpstring(Varname)
      else if (SDIM.eq.2) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif

      do im=1,ncomp
       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       ih=1
       Varname='F'
       ih=ih+1
       Varname(ih:ih)='D'
       ih=ih+1
       do i=1,2
        Varname(ih:ih)=matstr(i:i)
        ih=ih+1
       enddo
       call dumpstring(Varname)
      enddo ! im


      ! Zones
      iz_gb=1
       ! Zone marker. Value = 299.0
      write(11) ZONEMARKER
       ! Zone name
      Zonename = "ZONE"
      call dumpstring(Zonename)
 
      write(11) -1   ! Parent Zone
      write(11) 0    ! StrandID
      write(11) time ! Solution time
      write(11) -1   ! Not used. Set to -1
      write(11) 0    ! Zone Type
      write(11) 0    ! Specify Var Location. 0 = Don't specify, 
                     ! all data is located at the nodes.
      write(11) 0    ! Are raw local 1-to-1 face neighbors supplied?
      write(11) 0    ! Number of miscellaneous user-defined  
                      !face neighbor connections

        ! ----- IMax,JMax,KMax
      write(11) hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1
      write(11) hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1
      if (SDIM.eq.3) then
       write(11) hi_gb(iz_gb,SDIM)-lo_gb(iz_gb,SDIM)+1
      else if (SDIM.eq.2) then
       write(11) 1
      else
       print *,"dimension bust"
       stop
      endif
 
      write(11) 0

      write(11) EOHMARKER
 
      ! +++++++ DATA SECTION ++++++

      iz_gb=1

      do dir2=1,SDIM
       lo(dir2)=lo_gb(iz_gb,dir2)
       hi(dir2)=hi_gb(iz_gb,dir2)
      enddo

      if (SDIM.eq.3) then
       klo=lo(SDIM)
       khi=hi(SDIM)
      else if (SDIM.eq.2) then
       klo=0
       khi=0
      else
       print *,"dimension bust"
       stop
      endif

      allocate(zone_gb(iz_gb)% &
        var(nwrite,lo(1):hi(1),lo(2):hi(2),klo:khi))

       ! order is IMPORTANT.
      do k=klo,khi
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

       if (dir.eq.-1) then
        call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
       else if ((dir.ge.0).and.(dir.lt.SDIM)) then
        call gridstenMAC(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf,dir+1)
       else
        print *,"dir invalid tecplotfab 2"
        stop
       endif

       zone_gb(iz_gb)%var(1,i,j,k)=xsten(0,1)
       zone_gb(iz_gb)%var(2,i,j,k)=xsten(0,2)
       if (SDIM.eq.3) then
        zone_gb(iz_gb)%var(SDIM,i,j,k)=xsten(0,SDIM)
       else if (SDIM.eq.2) then
        ! do nothing
       else
        print *,"dimension bust"
        stop
       endif

       do im=1,ncomp
        zone_gb(iz_gb)%var(SDIM+im,i,j,k)=fabdata(D_DECL(i,j,k),im)
       enddo
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
       do k=klo,khi
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
      end subroutine FORT_TECPLOTFAB


      subroutine FORT_CELLGRID_SANITY( &
       tid, &
       data_dir, & ! data_dir=-1,0..sdim
       bfact, &
       ncomp, &
       datafab,DIMS(datafab), &
       problo, &
       probhi, &
       dx, &
       lo,hi, &
       level, &
       finest_level, &
       gridno, &
       rz_flag)
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: tid
      INTEGER_T, intent(in) :: data_dir
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: ncomp

      INTEGER_T, intent(in) :: rz_flag
      INTEGER_T, intent(in) :: lo(SDIM), hi(SDIM)
      INTEGER_T, intent(in) :: DIMDEC(datafab)
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: gridno
      REAL_T, intent(in) :: datafab(DIMV(datafab), &
        ncomp)
      REAL_T xposnd(SDIM)
      REAL_T xposndT(SDIM)
      REAL_T datand(ncomp)
      REAL_T writend(2*ncomp+2*SDIM)
      INTEGER_T scomp,iw

      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: probhi(SDIM)
      REAL_T, intent(in) :: dx(SDIM)

      character*3 levstr
      character*5 gridstr
      character*18 filename18

      INTEGER_T i,j,k
      INTEGER_T isub,jsub,ksub
      INTEGER_T isub_nrm
      INTEGER_T dir
      INTEGER_T i1,j1,k1
      INTEGER_T k1lo,k1hi
      INTEGER_T ksubhi
      INTEGER_T isubhi_local
      INTEGER_T jsubhi_local
      INTEGER_T ksubhi_local
      INTEGER_T igridlo(3),igridhi(3)
      INTEGER_T icell(3)
      INTEGER_T iproblo(3)
      INTEGER_T nhalf
      REAL_T xstenND(-3:3,SDIM)
      REAL_T dxleft,dxright
      INTEGER_T local_nd

      nhalf=3

      if (tid.ge.0) then
       ! do nothing
      else
       print *,"tid invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 39"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif

      if (data_dir.eq.-1) then 
       call checkbound(lo,hi,DIMS(datafab),0,-1,411)
      else if ((data_dir.ge.0).and.(data_dir.le.SDIM-1)) then
       call checkbound(lo,hi,DIMS(datafab),0,data_dir,411)
      else if (data_dir.eq.SDIM) then
       do dir=0,SDIM-1
        call checkbound(lo,hi,DIMS(datafab),0,dir,411)
       enddo
      else
       print *,"data_dir invalid"
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
       print *,"rz_flag invalid in cellgrid"
       stop
      endif

      write(levstr,'(I3)') level
      write(gridstr,'(I5)') gridno
      
      do i=1,3
        if (levstr(i:i).eq.' ') then
         levstr(i:i)='0'
        endif
      enddo
      do i=1,5
        if (gridstr(i:i).eq.' ') then
         gridstr(i:i)='0'
        endif
      enddo
      write(filename18,'(A10,A3,A5)') 'tempnddata',levstr,gridstr
      print *,"filename18 ",filename18

      open(unit=11,file=filename18)
      do dir=1,SDIM
       write(11,*) lo(dir),hi(dir)
      enddo

      igridlo(3)=0
      igridhi(3)=0
      iproblo(3)=0

      if (SDIM.eq.3) then
       k1lo=0
       k1hi=1
       ksubhi=2
      else if (SDIM.eq.2) then
       k1lo=0
       k1hi=0
       ksubhi=0
      else
       print *,"dimension bust"
       stop
      endif

      do dir=1,SDIM
       igridlo(dir)=lo(dir)
       igridhi(dir)=hi(dir)+1
       iproblo(dir)=0
      enddo

       ! the order k,j,i is IMPORTANT.
      do k=igridlo(3),igridhi(3) 
        if (k.eq.igridhi(3)) then
         ksubhi_local=0
        else
         ksubhi_local=ksubhi
        endif
        do ksub=0,ksubhi_local 
         do j=igridlo(2),igridhi(2) 
          if (j.eq.igridhi(2)) then
           jsubhi_local=0
          else
           jsubhi_local=2
          endif
          do jsub=0,jsubhi_local
           do i=igridlo(1),igridhi(1) 
            if (i.eq.igridhi(1)) then
             isubhi_local=0
            else
             isubhi_local=2
            endif
            do isub=0,isubhi_local

             ! iproblo=0
             call gridstenND(xstenND,problo,i,j,k,iproblo,bfact,dx,nhalf)
             do dir=1,SDIM
              dxleft=xstenND(0,dir)-xstenND(-1,dir)
              dxright=xstenND(1,dir)-xstenND(0,dir)

              if (bfact.eq.1) then
               if (abs(dxleft-dxright).le.VOFTOL*dx(dir)) then
                ! do nothing 
               else
                print *,"xstenND invalid"
                stop
               endif
              else if (bfact.gt.1) then
               if ((dxleft.gt.zero).and.(dxright.gt.zero)) then
                ! do nothing
               else
                print *,"xstenND invalid"
                stop
               endif
              else
               print *,"bfact invalid152"
               stop
              endif
               ! isub,jsub,ksub=0..2
              if (dir.eq.1) then
               isub_nrm=isub ! 0..2
              else if (dir.eq.2) then
               isub_nrm=jsub
              else if ((dir.eq.3).and.(SDIM.eq.3)) then
               isub_nrm=ksub
              else
               print *,"dir invalid"
               stop
              endif

              if (data_dir.eq.-1) then
               local_nd=0
              else if ((data_dir+1.eq.dir).or.(data_dir.eq.SDIM)) then
               local_nd=1
              else if ((data_dir+1.ne.dir).and.(data_dir.lt.SDIM)) then
               local_nd=0
              else
               print *,"data_dir invalid"
               stop
              endif
              dxright=xstenND(2,dir)-xstenND(0,dir)
              if (dxright.gt.zero) then
               if (isub_nrm.eq.0) then
                xposnd(dir)=xstenND(0,dir)
               else if (isub_nrm.eq.1) then
                if (local_nd.eq.0) then
                 xposnd(dir)=xstenND(0,dir)+dxright/100.0
                else if (local_nd.eq.1) then
                 xposnd(dir)=xstenND(0,dir)+(half-one/100.0)*dxright
                else
                 print *,"local_nd invalid"
                 stop
                endif
               else if (isub_nrm.eq.2) then
                if (local_nd.eq.0) then
                 xposnd(dir)=xstenND(2,dir)-dxright/100.0
                else if (local_nd.eq.1) then
                 xposnd(dir)=xstenND(2,dir)-(half-one/100.0)*dxright
                else
                 print *,"local_nd invalid"
                 stop
                endif
               else
                print *,"isub_nrm invalid"
                stop
               endif
              else
               print *,"dxright invalid"
               stop
              endif
             enddo ! dir = 1..sdim

             do dir=1,SDIM
              xposndT(dir)=xposnd(dir)
             enddo
             if (visual_RT_transform.eq.1) then ! defined in PROBCOMMON.F90
              call RT_transform(xposnd,xposndT)
             endif

             do dir=1,ncomp
              datand(dir)=zero
             enddo

             icell(1)=i
             icell(2)=j
             icell(3)=k

             do dir=1,SDIM

               ! isub,jsub,ksub=0..2
              if (dir.eq.1) then
               isub_nrm=isub ! 0..2
              else if (dir.eq.2) then
               isub_nrm=jsub
              else if ((dir.eq.3).and.(SDIM.eq.3)) then
               isub_nrm=ksub
              else
               print *,"dir invalid"
               stop
              endif

              if ((icell(dir).ge.lo(dir)).and. &
                  (icell(dir).le.hi(dir))) then

               if ((data_dir.eq.SDIM).or.(data_dir+1.eq.dir)) then
                if ((isub_nrm.eq.0).or.(isub_nrm.eq.1)) then
                 ! do nothing
                else if (isub_nrm.eq.2) then
                 icell(dir)=icell(dir)+1
                else
                 print *,"isub_nrm invalid"
                 stop
                endif
               else
                ! do nothing
               endif

              else if (icell(dir).eq.hi(dir)+1) then

               if ((data_dir.eq.SDIM).or.(data_dir+1.eq.dir)) then
                ! do nothing
               else
                icell(dir)=hi(dir)
               endif

              else
               print *,"icell out of range"
               stop
              endif
             enddo ! dir=1..sdim

             i1=icell(1)
             j1=icell(2)
             k1=icell(3)

             do dir=1,ncomp
              datand(dir)=datafab(D_DECL(i1,j1,k1),dir)
             enddo

             scomp=0
             do iw=1,SDIM
              writend(scomp+iw)=xposndT(iw)
             enddo

             scomp=scomp+SDIM

             do iw=1,ncomp
              writend(scomp+iw)=datand(iw)
             enddo

             scomp=scomp+ncomp 

             do iw=1,scomp
              if (iw.lt.scomp) then
               write(11,'(D25.16)',ADVANCE="NO") writend(iw)
              else if (iw.eq.scomp) then
               write(11,'(D25.16)') writend(iw)
              else
               print *,"iw invalid"
               stop
              endif
             enddo ! iw=1..scomp
            enddo ! isub=0..2
           enddo  ! i=igridlo(1),igridhi(1)
          enddo ! jsub=0..2
         enddo  ! j=igridlo(2),igridhi(2)
        enddo ! ksub=0..ksubhi
      enddo  ! k=igridlo(3),igridhi(3)  (sanity output loop to "AMR" grid)

      close(11)

      return
      end subroutine FORT_CELLGRID_SANITY


      subroutine FORT_COMBINEZONES_SANITY( &
       root_char_array, &
       n_root, &
       data_dir, &
       total_number_grids, &
       grids_per_level_array, &
       levels_array, &
       bfact_array, &
       gridno_array, &
       gridlo_array, &
       gridhi_array, &
       finest_level, &
       SDC_outer_sweeps, &
       slab_step, &
       data_id, &
       nsteps, &
       num_levels, &
       time, &
       visual_option, &
       visual_revolve, &
       ncomp)

      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: n_root
      character, dimension(n_root), intent(in) :: root_char_array
      INTEGER_T, intent(in) :: data_dir
      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: total_number_grids
      INTEGER_T, intent(in) :: num_levels
      INTEGER_T, intent(in) :: grids_per_level_array(num_levels)
      INTEGER_T, intent(in) :: levels_array(total_number_grids)
      INTEGER_T, intent(in) :: bfact_array(total_number_grids)
      INTEGER_T, intent(in) :: gridno_array(total_number_grids)
      INTEGER_T, intent(in) :: gridlo_array(total_number_grids*SDIM)
      INTEGER_T, intent(in) :: gridhi_array(total_number_grids*SDIM)
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: SDC_outer_sweeps
      INTEGER_T, intent(in) :: slab_step
      INTEGER_T, intent(in) :: data_id
      INTEGER_T, intent(in) :: nsteps
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: visual_option
      INTEGER_T, intent(in) :: visual_revolve

      INTEGER_T strandid
      INTEGER_T nwrite


      character*3 levstr
      character*5 gridstr
      character*18 filename18
      character*80 rmcommand

      character*6 stepstr
      character*3 outerstr
      character*3 slabstr
      character*3 idstr

      character(len=n_root) :: root_char_str
      character(len=n_root+36) :: newfilename40
      character(len=36) :: fname_extend
      character(len=4) :: step_chars
      character(len=2) :: dir_chars
      character(len=5) :: outer_chars
      character(len=4) :: slab_chars
      character(len=2) :: id_chars
      character(len=4) :: plt_chars

      INTEGER_T i,j,k,dir
      INTEGER_T ilev,igrid
      INTEGER_T lo(SDIM),hi(SDIM)
      INTEGER_T hi_index_shift(3)

! Guibo

      character*80 Title,Zonename
      REAL*4 ZONEMARKER,EOHMARKER
      integer*4 :: nzones_gb,iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb
      INTEGER_T bfact,testlev,testgridno
      INTEGER_T testlo(SDIM),testhi(SDIM)

      ! define zone structure
      type zone_t
         real*8, pointer :: var(:,:,:,:)
      end type zone_t
      type(zone_t), dimension(:), allocatable :: zone_gb

      INTEGER_T plot_sdim,klo_plot,khi_plot

! Guibo
      INTEGER_T sysret


      plot_sdim=SDIM

      if ((levelrz.eq.0).or.(levelrz.eq.3)) then
       if (visual_revolve.ne.0) then
        print *,"visual_revolve= ",visual_revolve
        print *,"visual_revolve invalid combine zones"
        stop
       endif
      else if (levelrz.eq.1) then
       if (SDIM.ne.2) then
        print *,"dimension bust"
        stop
       endif
       if (visual_revolve.eq.0) then
        ! do nothing
       else if ((visual_revolve.ge.1).and.(visual_revolve.le.1024)) then
        ! do nothing
       else
        print *,"visual_revolve= ",visual_revolve
        print *,"visual_revolve out of range"
        stop
       endif
      else 
       print *,"levelrz invalid combine zones"
       stop
      endif

      nwrite=plot_sdim+ncomp

      if (num_levels.ne.finest_level+1) then
       print *,"num_levels invalid"
       stop
      endif

      if ((visual_revolve.ge.1).and.(visual_revolve.le.1024)) then

       if (levelrz.ne.1) then
        print *,"levelrz invalid combine zones 2"
        stop
       endif
       if (SDIM.ne.2) then
        print *,"parameter bust"
        stop
       endif

       plot_sdim=3

       call zones_revolve_sanity( &
        root_char_array, &
        n_root, &
        data_dir, &
        plot_sdim, &
        total_number_grids, &
        grids_per_level_array, &
        levels_array, &
        bfact_array, &
        gridno_array, &
        gridlo_array, &
        gridhi_array, &
        finest_level, &
        SDC_outer_sweeps, &
        slab_step, &
        data_id, &
        nsteps, &
        num_levels, &
        time, &
        visual_option, &
        visual_revolve, &
        ncomp)
      else if (visual_revolve.eq.0) then

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
        print *,"levelrz invalid combine zones 3"
        stop
       endif

       write(stepstr,'(I6)') nsteps
       do i=1,6
        if (stepstr(i:i).eq.' ') then
         stepstr(i:i)='0'
        endif
       enddo

       write(outerstr,'(I3)') SDC_outer_sweeps
       write(slabstr,'(I3)') slab_step
       write(idstr,'(I3)') data_id
       do i=1,3
        if (outerstr(i:i).eq.' ') then
         outerstr(i:i)='0'
        endif
        if (slabstr(i:i).eq.' ') then
         slabstr(i:i)='0'
        endif
        if (idstr(i:i).eq.' ') then
         idstr(i:i)='0'
        endif
       enddo

       do i=1,n_root
        root_char_str(i:i)=root_char_array(i)
       enddo

       if (data_dir.eq.-1) then
        dir_chars='CC'
       else if (data_dir.eq.0) then
        dir_chars='XC'
       else if (data_dir.eq.1) then
        dir_chars='YC'
       else if ((data_dir.eq.SDIM-1).and.(SDIM.eq.3)) then
        dir_chars='ZC'
       else if (data_dir.eq.SDIM) then
        dir_chars='ND'
       else
        print *,"data_dir invalid"
        stop
       endif
       step_chars='step'
       outer_chars='outer'
       slab_chars='slab'
       id_chars='id'
       plt_chars='.plt'

       newfilename40(1:n_root)=root_char_str
       write(fname_extend,'(A2,A2,A3,A4,A6,A5,A3,A4,A3,A4)') &
               dir_chars,id_chars,idstr, &
               step_chars,stepstr,outer_chars,outerstr, &
               slab_chars,slabstr,plt_chars

       newfilename40(n_root+1:n_root+36)=fname_extend

       print *,"newfilename40 ",newfilename40


       !--------------------------------------------------
       ! Determine nzones_gb and allocate zone_gb, lo_gb, hi_gb
       nzones_gb=0
       do ilev=0,finest_level
       do igrid=0,grids_per_level_array(ilev+1)-1
        nzones_gb = nzones_gb+1
       enddo
       enddo

       if (nzones_gb.ne.total_number_grids) then
        print *,"nzones_gb.ne.total_number_grids: aborting"
        stop
       endif

       print *,"allocating grid structure for tecplot binary"
       print *,"finest_level= ",finest_level
       print *,"number of grids on finest level ", &
        grids_per_level_array(finest_level+1) 
       print *,"total number of grid blocks: ",nzones_gb
       print *,"num_materials ",num_materials
       print *,"num_state_material ",num_state_material

       allocate(zone_gb(nzones_gb))
       allocate(lo_gb(nzones_gb,plot_sdim))
       allocate(hi_gb(nzones_gb,plot_sdim))

        ! Determine lo_gb, hi_gb
       iz_gb=0
       do ilev=0,finest_level
       do igrid=0,grids_per_level_array(ilev+1)-1
        write(levstr,'(I3)') ilev
        write(gridstr,'(I5)') igrid
        do i=1,3
         if (levstr(i:i).eq.' ') then
          levstr(i:i)='0'
         endif
        enddo
        do i=1,5
         if (gridstr(i:i).eq.' ') then
          gridstr(i:i)='0'
         endif
        enddo
        write(filename18,'(A10,A3,A5)') 'tempnddata',levstr,gridstr
        open(unit=4,file=filename18)

        do dir=1,plot_sdim
         read(4,*) lo(dir),hi(dir)
         lo_gb(iz_gb+1,dir)=lo(dir)
         hi_gb(iz_gb+1,dir)=hi(dir)
        enddo ! dir

        bfact=bfact_array(iz_gb+1)
        testlev=levels_array(iz_gb+1)
        testgridno=gridno_array(iz_gb+1)

        if (bfact.lt.1) then
         print *,"bfact invalid156"
         stop
        endif
        if (testlev.ne.ilev) then
         print *,"testlev invalid"
         stop
        endif
        if (testgridno.ne.igrid) then
         print *,"testgridno invalid"
         stop
        endif

        do dir=1,plot_sdim
         testlo(dir)=gridlo_array(plot_sdim*iz_gb+dir)
         testhi(dir)=gridhi_array(plot_sdim*iz_gb+dir)
         if (testlo(dir).ne.lo(dir)) then
          print *,"testlo invalid"
          stop
         endif
         if (testhi(dir).ne.hi(dir)) then
          print *,"testhi invalid"
          stop
         endif
         if ((lo(dir)/bfact)*bfact.ne.lo(dir)) then
          print *,"lo not divisible by bfact"
          stop
         endif
         if (((hi(dir)+1)/bfact)*bfact.ne.hi(dir)+1) then
          print *,"hi+1 not divisible by bfact"
          stop
         endif
        enddo ! dir

        close(4)

        iz_gb=iz_gb+1
       enddo ! igrid
       enddo ! ilev

       if (iz_gb.ne.total_number_grids) then
        print *,"iz_gb or total_number_grids invalid"
        stop
       endif

       !-----------------------------------------------------------
       ZONEMARKER = 299.0
       EOHMARKER  = 357.0
       open(unit=11,file=newfilename40,form="unformatted",access="stream")

        ! +++++++ HEADER SECTION ++++++

        ! i.  Magic number, Version number
       write(11) "#!TDV112"
 
        ! ii. Integer value of 1.
       write(11) 1

        ! iii. Title and variable names.
        ! File type 0 = FULL,1 = GRID,2 = SOLUTION
       write(11) 0
        ! The TITLE
       Title = "CLSVOF data"
       call dumpstring(Title)
       ! Number of variables 
       write(11) nwrite

       ! Variable names: combinezones_sanity
       call dumpstring_headers_sanity(plot_sdim,ncomp)

        ! Zones
       do iz_gb=1,nzones_gb
        ! Zone marker. Value = 299.0
        write(11) ZONEMARKER
        ! Zone name
        Zonename = "ZONE"
        call dumpstring(Zonename)

        strandid=1
 
        write(11) -1   ! Parent Zone
        write(11) 0    ! StrandID
        write(11) time ! Solution time
        write(11) -1   ! Not used. Set to -1
        write(11) 0    ! Zone Type
        write(11) 0    ! Specify Var Location. 0 = Don't specify, 
                       ! all data is located at the nodes.
        write(11) 0    ! Are raw local 1-to-1 face neighbors supplied?
        write(11) 0    ! Number of miscellaneous user-defined  
                       !face neighbor connections

        hi_index_shift(1)=3*(hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1)+1
        hi_index_shift(2)=3*(hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1)+1
        if (plot_sdim.eq.3) then
         hi_index_shift(3)=3*(hi_gb(iz_gb,plot_sdim)-lo_gb(iz_gb,plot_sdim)+1)+1
        else if (plot_sdim.eq.2) then
         hi_index_shift(3)=1
        else
         print *,"plot_sdim invalid in combinezones_sanity"
         stop
        endif

         ! ----- IMax,JMax,KMax
        write(11) hi_index_shift(1)
        write(11) hi_index_shift(2)
        write(11) hi_index_shift(3)
 
        write(11) 0
       enddo  ! iz_gb

       write(11) EOHMARKER
 
       print *,"finished writing header section, now data section"
       print *,"nzones_gb= ",nzones_gb
       print *,"nwrite= ",nwrite

        ! +++++++ DATA SECTION ++++++

       ilev=0
       igrid=0

       do iz_gb=1,nzones_gb

        if (ilev.gt.finest_level) then
         print *,"ilev invalid"
         print *,"ilev=",ilev
         print *,"finest_level=",finest_level
         stop
        endif
      
        do while (grids_per_level_array(ilev+1).eq.0)
         ilev=ilev+1

         if (igrid.ne.0) then
          print *,"igrid should be 0"
          print *,"igrid=",igrid
          stop
         endif
         if (ilev.gt.finest_level) then
          print *,"ilev invalid in grids_per_level loop"
          print *,"ilev=",ilev
          print *,"finest_level=",finest_level
          stop
         endif
        enddo ! while grids_per_level_array(ilev+1)==0
 
        if (igrid.gt.grids_per_level_array(ilev+1)-1) then
         print *,"igrid invalid"
         print *,"igrid= ",igrid
         print *,"ilev= ",ilev
         print *,"finest_level=",finest_level
         print *,"grids_per_level_array= ", &
          grids_per_level_array(ilev+1)
         print *,"iz_gb = ",iz_gb
         print *,"nzones_gb= ",nzones_gb
         print *,"num_levels= ",num_levels
         print *,"nwrite= ",nwrite
         stop
        endif

        if (igrid.ne.gridno_array(iz_gb)) then
         print *,"igrid invalid"
         stop
        endif
        if (ilev.ne.levels_array(iz_gb)) then
         print *,"ilev invalid"
         stop
        endif

        do dir=1,plot_sdim
         lo(dir)=lo_gb(iz_gb,dir)
         hi(dir)=hi_gb(iz_gb,dir)
        enddo
      
        if (plot_sdim.eq.3) then
         klo_plot=lo(plot_sdim)
         khi_plot=hi(plot_sdim)+1
        else if (plot_sdim.eq.2) then
         klo_plot=0
         khi_plot=0
        else
         print *,"plot_sdim invalid"
         stop
        endif

        hi_index_shift(1)=3*(hi(1)-lo(1)+1)
        hi_index_shift(2)=3*(hi(2)-lo(2)+1)
        hi_index_shift(3)=3*(khi_plot-klo_plot)
 
        allocate(zone_gb(iz_gb)% &
         var(nwrite, &
            0:hi_index_shift(1), &
            0:hi_index_shift(2), &
            0:hi_index_shift(3)))
         
        write(levstr,'(I3)') ilev
        write(gridstr,'(I5)') igrid
        do i=1,3
         if (levstr(i:i).eq.' ') then
          levstr(i:i)='0'
         endif
        enddo
        do i=1,5
         if (gridstr(i:i).eq.' ') then
          gridstr(i:i)='0'
         endif
        enddo

        write(filename18,'(A10,A3,A5)') 'tempnddata',levstr,gridstr
        open(unit=4,file=filename18)
        print *,"filename18 ",filename18

        do dir=1,SDIM
         read(4,*) lo(dir),hi(dir)

         if (lo(dir).ne.lo_gb(iz_gb,dir)) then
          print *,"lo and lo_gb different"
          print *,"iz_gb,dir,lo,lo_gb ",iz_gb,dir,lo(dir), &
           lo_gb(iz_gb,dir)
          stop
         endif
         if (hi(dir).ne.hi_gb(iz_gb,dir)) then
          print *,"hi and hi_gb different"
          print *,"iz_gb,dir,hi,hi_gb ",iz_gb,dir,hi(dir), &
           hi_gb(iz_gb,dir)
          stop
         endif

        enddo  ! dir

         ! order is IMPORTANT.
        do k=0,hi_index_shift(3)
        do j=0,hi_index_shift(2)
        do i=0,hi_index_shift(1)
         read(4,*) (zone_gb(iz_gb)%var(ivar_gb,i,j,k),ivar_gb=1,nwrite)
        enddo
        enddo
        enddo

        close(4)

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
         write(11) minval(zone_gb(iz_gb)%var(ivar_gb,:,:,:))
         write(11) maxval(zone_gb(iz_gb)%var(ivar_gb,:,:,:))
        enddo

        ! order is IMPORTANT
        do ivar_gb=1,nwrite
         hi_index_shift(1)=3*(hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1)
         hi_index_shift(2)=3*(hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1)
         hi_index_shift(3)=3*(khi_plot-klo_plot)
         do k=0,hi_index_shift(3)
         do j=0,hi_index_shift(2)
         do i=0,hi_index_shift(1)
          write(11) zone_gb(iz_gb)%var(ivar_gb,i,j,k)
         enddo
         enddo
         enddo
        enddo

        deallocate(zone_gb(iz_gb)%var)

        igrid=igrid+1
        if (igrid.gt.grids_per_level_array(ilev+1)-1) then
         ilev=ilev+1
         igrid=0
        endif

       enddo  ! iz_gb

       deallocate(zone_gb)
       deallocate(lo_gb)
       deallocate(hi_gb)

       close(11)
     
       rmcommand='rm tempnddata*'

       print *,"issuing command ",rmcommand

       sysret=0

#ifdef PGIFORTRAN
       call system(rmcommand)
#else
       call execute_command_line(rmcommand,exitstat=sysret)
#endif
       if (sysret.ne.0) then
        print *,"execute_command_line has sysret=",sysret
        stop
       endif

      else
       print *,"visual_revolve invalid"
       stop
      endif

      return
      end subroutine FORT_COMBINEZONES_SANITY

      subroutine FORT_TECPLOTFAB_SANITY( &
       root_char_array, &
       n_root, &
       data_dir, &
       bfact, &
       fablo, &
       fabhi, &
       datafab,DIMS(datafab), &
       problo, &
       probhi, &
       dx, &
       SDC_outer_sweeps, &
       slab_step, &
       data_id, &
       nsteps, &
       time, &
       visual_option, &
       visual_revolve, &
       level, &
       finest_level, &
       ncomp)
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: n_root
      character, dimension(n_root), intent(in) :: root_char_array
      INTEGER_T, intent(in) :: data_dir
      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: SDC_outer_sweeps
      INTEGER_T, intent(in) :: slab_step
      INTEGER_T, intent(in) :: data_id
      INTEGER_T, intent(in) :: nsteps
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: visual_option
      INTEGER_T, intent(in) :: visual_revolve
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM) 
      INTEGER_T, intent(in) :: DIMDEC(datafab)
      REAL_T, intent(in) :: datafab(DIMV(datafab),ncomp)
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: probhi(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T :: tid_local
      INTEGER_T :: gridno_local
      INTEGER_T :: dir_local
      INTEGER_T :: total_number_grids
      INTEGER_T :: num_levels
      INTEGER_T :: grids_per_level_array(1)
      INTEGER_T :: levels_array(1)
      INTEGER_T :: bfact_array(1)
      INTEGER_T :: gridno_array(1)
      INTEGER_T :: gridlo_array(SDIM)
      INTEGER_T :: gridhi_array(SDIM)

      tid_local=0
      gridno_local=0

      call FORT_CELLGRID_SANITY( &
       tid_local, &
       data_dir, &
       bfact, &
       ncomp, &
       datafab,DIMS(datafab), &
       problo, &
       probhi, &
       dx, &
       fablo,fabhi, &
       level, &
       finest_level, &
       gridno_local, &
       levelrz)

      total_number_grids=1
      num_levels=1
      grids_per_level_array(1)=1
      levels_array(1)=0
      bfact_array(1)=bfact
      gridno_array(1)=0
      do dir_local=1,SDIM
       gridlo_array(dir_local)=fablo(dir_local)
       gridhi_array(dir_local)=fabhi(dir_local)
      enddo

      call FORT_COMBINEZONES_SANITY( &
       root_char_array, &
       n_root, &
       data_dir, &
       total_number_grids, &
       grids_per_level_array, &
       levels_array, &
       bfact_array, &
       gridno_array, &
       gridlo_array, &
       gridhi_array, &
       finest_level, &
       SDC_outer_sweeps, &
       slab_step, &
       data_id, &
       nsteps, &
       num_levels, &
       time, &
       visual_option, &
       visual_revolve, &
       ncomp)

      return
      end subroutine FORT_TECPLOTFAB_SANITY

#if (STANDALONE==1)
      end module tecplotutil_cpp_module
#endif

#undef STANDALONE
