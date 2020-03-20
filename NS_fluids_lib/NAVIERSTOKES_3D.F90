#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "NAVIERSTOKES_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

#define AREAZERO (1.0D-12)

      module navierstokesf90_module
      use probf90_module

      contains

      subroutine vinterp(valu,gridval,gridx,gridy,gridz, &
        IV0,IV1,XX,YY,ZZ,icomp)

      IMPLICIT NONE

      REAL_T valu
      REAL_T gridx(8),gridy(8),gridz(8),gridval(8)
      INTEGER_T IV0,IV1
      REAL_T XX(SDIM),YY(SDIM),ZZ(SDIM)
      INTEGER_T icomp
      REAL_T tt

      if (abs(gridval(IV0)-valu).le.1.0D-10) then
       XX(icomp)=gridx(IV0)
       YY(icomp)=gridy(IV0)
#if (AMREX_SPACEDIM==3)
       ZZ(icomp)=gridz(IV0)
#endif
      else if (abs(gridval(IV1)-valu).le.1.0D-10) then
       XX(icomp)=gridx(IV1)
       YY(icomp)=gridy(IV1)
#if (AMREX_SPACEDIM==3)
       ZZ(icomp)=gridz(IV1)
#endif
      else
       tt=(gridval(IV0)-valu)/(gridval(IV0)-gridval(IV1))
       if ((tt.lt.zero).or.(tt.gt.one)) then
        print *,"tt invalid"
        stop
       endif
 
       XX(icomp)=tt*gridx(IV1)+(one-tt)*gridx(IV0)
       YY(icomp)=tt*gridy(IV1)+(one-tt)*gridy(IV0)
#if (AMREX_SPACEDIM==3)
       ZZ(icomp)=tt*gridz(IV1)+(one-tt)*gridz(IV0)
#endif
      endif

      return
      end subroutine

! in 2d, NORMAL(3)=0
! fictitious node at x1,y1,1
      subroutine addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)

      IMPLICIT NONE

      INTEGER_T itri
      INTEGER_T imaxtri
      REAL_T geom(SDIM,imaxtri)
      REAL_T XX(SDIM),YY(SDIM),ZZ(SDIM),NORMAL(SDIM)
      REAL_T CPROD(SDIM),VEC1(3),VEC2(3),DOTPROD

#if (AMREX_SPACEDIM==3)
      VEC1(1)=XX(SDIM)-XX(1)
      VEC1(2)=YY(SDIM)-YY(1)
      VEC1(3)=ZZ(SDIM)-ZZ(1)
#elif (AMREX_SPACEDIM==2)
      VEC1(1)=zero
      VEC1(2)=zero
      VEC1(3)=one
#else
      print *,"dimension bust"
      stop
#endif

      VEC2(1)=XX(2)-XX(1)
      VEC2(2)=YY(2)-YY(1)
#if (AMREX_SPACEDIM==3)
      VEC2(3)=ZZ(2)-ZZ(1)
#elif (AMREX_SPACEDIM==2)
      VEC2(3)=zero
#else
      print *,"dimension bust"
      stop
#endif
      CPROD(1)=VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
      CPROD(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
#if (AMREX_SPACEDIM==3)
      CPROD(SDIM)=VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
#endif
      DOTPROD=CPROD(1)*NORMAL(1)+CPROD(2)*NORMAL(2)
#if (AMREX_SPACEDIM==3)
      DOTPROD=DOTPROD+CPROD(SDIM)*NORMAL(SDIM)
#endif

      if (itri+SDIM.gt.imaxtri-1) then
       print *,"itri too big"
       stop
      endif

      geom(1,itri+1)=XX(1)
      geom(2,itri+1)=YY(1)
#if (AMREX_SPACEDIM==3)
      geom(SDIM,itri+1)=ZZ(1)
#endif
#if (AMREX_SPACEDIM==3)
      if (DOTPROD.gt.0.0) then
       geom(1,itri+2)=XX(2)
       geom(2,itri+2)=YY(2)
       geom(SDIM,itri+2)=ZZ(2)
       geom(1,itri+SDIM)=XX(SDIM)
       geom(2,itri+SDIM)=YY(SDIM)
       geom(SDIM,itri+SDIM)=ZZ(SDIM)
      else
       geom(1,itri+2)=XX(SDIM)
       geom(2,itri+2)=YY(SDIM)
       geom(3,itri+2)=ZZ(SDIM)
       geom(1,itri+SDIM)=XX(2)
       geom(2,itri+SDIM)=YY(2)
       geom(SDIM,itri+SDIM)=ZZ(2)
      endif
#elif (AMREX_SPACEDIM==2)
      geom(1,itri+2)=XX(2)
      geom(2,itri+2)=YY(2)
#else
      print *,"dimension bust"
      stop
#endif


      itri=itri+SDIM

      return
      end subroutine


      subroutine polysegment(gridx,gridy,gridz,gridval,valu,geom, &
        itri,IV0,IV1,IV2,IV3,imaxtri)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T itri
      INTEGER_T IV0,IV1,IV2,IV3
      INTEGER_T imaxtri
      REAL_T geom(SDIM,imaxtri)
      REAL_T valu
      REAL_T gridx(8),gridy(8),gridz(8),gridval(8)
      REAL_T AA(4,4),sourcex(4),bb(4),NORMAL(SDIM)
      REAL_T XX(SDIM),YY(SDIM),ZZ(SDIM)
      INTEGER_T istat,idxtri

      idxtri=0
      if (gridval(IV0).lt.valu) then
       idxtri=idxtri+1
      endif
      if (gridval(IV1).lt.valu) then
       idxtri=idxtri+2
      endif
      if (gridval(IV2).lt.valu) then
       idxtri=idxtri+4
      endif

      if ((idxtri.ne.7).and.(idxtri.ne.0)) then
       AA(1,1)=gridx(IV0) 
       AA(1,2)=gridy(IV0) 
       AA(1,3)=zero
       AA(1,4)=1.0
       AA(2,1)=gridx(IV1) 
       AA(2,2)=gridy(IV1) 
       AA(2,3)=zero
       AA(2,4)=1.0
       AA(3,1)=gridx(IV2) 
       AA(3,2)=gridy(IV2) 
       AA(3,3)=zero
       AA(3,4)=1.0
       AA(4,1)=gridx(IV0) 
       AA(4,2)=gridy(IV0) 
       AA(4,3)=one
       AA(4,4)=1.0
       bb(1)=gridval(IV0)
       bb(2)=gridval(IV1)
       bb(3)=gridval(IV2)
       bb(4)=gridval(IV0)
       call matrix_solve(AA,sourcex,bb,istat,4)
       if (istat.ne.0) then
        NORMAL(1)=sourcex(1)
        NORMAL(2)=sourcex(2)
        if ((idxtri.eq.1).or.(idxtri.eq.6)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.2).or.(idxtri.eq.5)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.3).or.(idxtri.eq.4)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        endif
       endif ! istat.ne.zero
      endif ! idxtri.ne.7 or 0

      return
      end subroutine


      subroutine polytri(gridx,gridy,gridz,gridval,valu,geom, &
        itri,IV0,IV1,IV2,IV3,imaxtri)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T itri
      INTEGER_T IV0,IV1,IV2,IV3
      INTEGER_T imaxtri
      REAL_T geom(SDIM,imaxtri)
      REAL_T valu
      REAL_T gridx(8),gridy(8),gridz(8),gridval(8)
      REAL_T AA(4,4),sourcex(4),bb(4),NORMAL(SDIM)
      REAL_T XX(SDIM),YY(SDIM),ZZ(SDIM)
      INTEGER_T istat,idxtri

      idxtri=0
      if (gridval(IV0).lt.valu) then
       idxtri=idxtri+1
      endif
      if (gridval(IV1).lt.valu) then
       idxtri=idxtri+2
      endif
      if (gridval(IV2).lt.valu) then
       idxtri=idxtri+4
      endif
      if (gridval(IV3).lt.valu) then
       idxtri=idxtri+8
      endif

      if ((idxtri.ne.15).and.(idxtri.ne.0)) then
       AA(1,1)=gridx(IV0) 
       AA(1,2)=gridy(IV0) 
       AA(1,3)=gridz(IV0) 
       AA(1,4)=1.0
       AA(2,1)=gridx(IV1) 
       AA(2,2)=gridy(IV1) 
       AA(2,3)=gridz(IV1) 
       AA(2,4)=1.0
       AA(3,1)=gridx(IV2) 
       AA(3,2)=gridy(IV2) 
       AA(3,3)=gridz(IV2) 
       AA(3,4)=1.0
       AA(4,1)=gridx(IV3) 
       AA(4,2)=gridy(IV3) 
       AA(4,3)=gridz(IV3) 
       AA(4,4)=1.0
       bb(1)=gridval(IV0)
       bb(2)=gridval(IV1)
       bb(3)=gridval(IV2)
       bb(4)=gridval(IV3)
       call matrix_solve(AA,sourcex,bb,istat,4)
       if (istat.ne.0) then
        NORMAL(1)=sourcex(1)
        NORMAL(2)=sourcex(2)
        NORMAL(SDIM)=sourcex(SDIM)
        if ((idxtri.eq.1).or.(idxtri.eq.14)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.2).or.(idxtri.eq.13)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.3).or.(idxtri.eq.12)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV3,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
 
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.4).or.(idxtri.eq.11)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV1,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.5).or.(idxtri.eq.10)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
  
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.6).or.(idxtri.eq.9)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV1,IV3,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
 
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV1,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV0,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV2,IV3,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        else if ((idxtri.eq.7).or.(idxtri.eq.8)) then
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV0,XX,YY,ZZ,1)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV2,XX,YY,ZZ,2)
         call vinterp(valu,gridval,gridx,gridy,gridz, &
          IV3,IV1,XX,YY,ZZ,3)
         call addgeom(geom,itri,imaxtri,XX,YY,ZZ,NORMAL)
        endif
       endif ! istat.ne.zero
      endif ! idxtri.ne.15 or 0

      return
      end subroutine polytri


      subroutine zones_revolve( &
       plot_sdim, &
       total_number_grids, &
       grids_per_level_array, &
       levels_array, &
       bfact_array, &
       gridno_array, &
       gridlo_array, &
       gridhi_array, &
       finest_level, &
       nsteps, &
       num_levels, &
       time, &
       visual_option, &
       visual_revolve, &
       plotint, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map)
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T plot_sdim,klo_plot,khi_plot

      INTEGER_T nparts
      INTEGER_T nparts_def
      INTEGER_T im_solid_map(nparts_def)
      INTEGER_T total_number_grids
      INTEGER_T num_levels
      INTEGER_T grids_per_level_array(num_levels)
      INTEGER_T levels_array(total_number_grids)
      INTEGER_T bfact_array(total_number_grids)
      INTEGER_T gridno_array(total_number_grids)
      INTEGER_T gridlo_array(total_number_grids*SDIM)
      INTEGER_T gridhi_array(total_number_grids*SDIM)
      INTEGER_T finest_level
      INTEGER_T nsteps
      REAL_T time
      INTEGER_T visual_option
      INTEGER_T visual_revolve
      INTEGER_T plotint
      INTEGER_T nmat

      INTEGER_T strandid

      INTEGER_T nwrite3d,nwrite2d,index3d,index2d

      INTEGER_T im
      INTEGER_T partid

      character*3 levstr
      character*5 gridstr
      character*18 filename18
      character*80 rmcommand

      character*6 stepstr
      character*16 newfilename16

      INTEGER_T i,j,k,dir
      INTEGER_T ilev,igrid,ivel2d,ivel3d
      INTEGER_T lo(plot_sdim),hi(plot_sdim)
      INTEGER_T sysret


! Guibo

      character*80 Title,Zonename
      REAL*4 ZONEMARKER,EOHMARKER
      integer*4 :: nzones_gb,iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb
      INTEGER_T bfact,testlev,testgridno
      INTEGER_T testlo(SDIM),testhi(SDIM)

      ! define zone structure
      type zone3d_t
         real*8, pointer :: var(:,:,:,:)
      end type zone3d_t
      type(zone3d_t), dimension(:), allocatable :: zone3d_gb

      type zone2d_t
         real*8, pointer :: var(:,:,:)
      end type zone2d_t
      type(zone2d_t), dimension(:), allocatable :: zone2d_gb

      REAL_T theta,rr,zz,xx,yy,ur,uz,ux,uy

      if (plot_sdim.ne.3) then
       print *,"plot_sdim invalid"
       stop
      endif

      if (SDIM.ne.2) then
       print *,"must be called only in 2D"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if (levelrz.ne.1) then
       print *,"levelrz invalid zones revolve"
       stop
      endif
      if (visual_revolve.lt.1) then
       print *,"visual_revolve: ",visual_revolve
       print *,"visual_revolve invalid zones_revolve"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat bust"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid zones_revolve"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid zones_revolve"
       stop
      endif

      call get_nwrite(SDIM,nwrite2d)
      call get_nwrite(plot_sdim,nwrite3d)

      if (num_levels.ne.finest_level+1) then
       print *,"num_levels invalid"
       stop
      endif

      write(stepstr,'(I6)') nsteps
      do i=1,6
       if (stepstr(i:i).eq.' ') then
        stepstr(i:i)='0'
       endif
      enddo

      write(newfilename16,'(A6,A6,A4)') 'nddata',stepstr,'.plt'

      print *,"newfilename16 ",newfilename16


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

      allocate(zone2d_gb(nzones_gb))
      allocate(zone3d_gb(nzones_gb))
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
        if (dir.ne.3) then
         read(4,*) lo(dir),hi(dir)
        else if (dir.eq.3) then
         lo(dir)=0
         hi(dir)=visual_revolve-1
        else
         print *,"dir invalid zones revolve"
         stop
        endif
        lo_gb(iz_gb+1,dir)=lo(dir)
        hi_gb(iz_gb+1,dir)=hi(dir)
       enddo ! dir

       bfact=bfact_array(iz_gb+1)
       testlev=levels_array(iz_gb+1)
       testgridno=gridno_array(iz_gb+1)

       if (bfact.lt.1) then
        print *,"bfact invalid150"
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

       do dir=1,SDIM
        testlo(dir)=gridlo_array(SDIM*iz_gb+dir)
        testhi(dir)=gridhi_array(SDIM*iz_gb+dir)
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

      !-------------------------------------------------------------
      ZONEMARKER = 299.0
      EOHMARKER  = 357.0
      open(unit=11,file=newfilename16,form="unformatted",access="stream")

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
      write(11) nwrite3d

      ! Variable names: zones_revolve
      call dumpstring_headers(plot_sdim)

       ! Zones
      do iz_gb=1,nzones_gb
       ! Zone marker. Value = 299.0
       write(11) ZONEMARKER
       ! Zone name
       Zonename = "ZONE"
       call dumpstring(Zonename)

       if (plotint.le.0) then
        strandid=1    
       else
        strandid=(nsteps/plotint)+1
       endif

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
       write(11) hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+2
       write(11) hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+2
       if (plot_sdim.eq.3) then
        write(11) hi_gb(iz_gb,plot_sdim)-lo_gb(iz_gb,plot_sdim)+2
       else
        print *,"plot_sdim invalid in zones_revolve"
        stop
       endif
 
       write(11) 0
      enddo  ! iz_gb

      write(11) EOHMARKER

      print *,"finished writing header section, now data section"
      print *,"nzones_gb= ",nzones_gb
      print *,"nwrite3d= ",nwrite3d

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
       enddo  ! while grids_per_level_array==0
 
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
        print *,"nwrite3d= ",nwrite3d
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

       allocate(zone3d_gb(iz_gb)% &
        var(nwrite3d,lo(1):hi(1)+1, &
            lo(2):hi(2)+1, &
            klo_plot:khi_plot))

       allocate(zone2d_gb(iz_gb)% &
        var(nwrite2d,lo(1):hi(1)+1, &
            lo(SDIM):hi(SDIM)+1))
 
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

       dir=plot_sdim
       lo(dir)=0
       hi(dir)=visual_revolve-1

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

        ! order is IMPORTANT.
       do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)+1
        read(4,*) (zone2d_gb(iz_gb)%var(ivar_gb,i,j),ivar_gb=1,nwrite2d)

        do k=klo_plot,khi_plot

         index3d=0
         index2d=0

         theta=two*Pi*k/visual_revolve
         rr=zone2d_gb(iz_gb)%var(1,i,j)
         zz=zone2d_gb(iz_gb)%var(2,i,j)
         xx=rr*cos(theta)
         yy=rr*sin(theta)
         zone3d_gb(iz_gb)%var(1,i,j,k)=xx
         zone3d_gb(iz_gb)%var(2,i,j,k)=yy
         zone3d_gb(iz_gb)%var(3,i,j,k)=zz

         index3d=index3d+plot_sdim
         index2d=index2d+SDIM

         do im=1,num_materials_vel
          ivel2d=(im-1)*SDIM
          ivel3d=(im-1)*plot_sdim
          ur=zone2d_gb(iz_gb)%var(index2d+1+ivel2d,i,j)
          uz=zone2d_gb(iz_gb)%var(index2d+2+ivel2d,i,j)
          ux=ur*cos(theta)
          uy=ur*sin(theta)
          zone3d_gb(iz_gb)%var(index3d+1+ivel3d,i,j,k)=ux
          zone3d_gb(iz_gb)%var(index3d+2+ivel3d,i,j,k)=uy
          zone3d_gb(iz_gb)%var(index3d+3+ivel3d,i,j,k)=uz
         enddo ! im

         index3d=index3d+plot_sdim*num_materials_vel
         index2d=index2d+SDIM*num_materials_vel

         do partid=0,nparts_def-1
          ivel2d=partid*SDIM
          ivel3d=partid*plot_sdim
          ur=zone2d_gb(iz_gb)%var(index2d+1+ivel2d,i,j)
          uz=zone2d_gb(iz_gb)%var(index2d+2+ivel2d,i,j)
          ux=ur*cos(theta)
          uy=ur*sin(theta)
          zone3d_gb(iz_gb)%var(index3d+1+ivel3d,i,j,k)=ux
          zone3d_gb(iz_gb)%var(index3d+2+ivel3d,i,j,k)=uy
          zone3d_gb(iz_gb)%var(index3d+3+ivel3d,i,j,k)=uz
         enddo ! partid=0..nparts_def-1

         index3d=index3d+plot_sdim*nparts_def
         index2d=index2d+SDIM*nparts_def

          ! pressure,presder,div,divdat,mach,vof
         do ivar_gb=1,5*num_materials_vel+nmat
          index3d=index3d+1
          index2d=index2d+1
          zone3d_gb(iz_gb)%var(index3d,i,j,k)= &
            zone2d_gb(iz_gb)%var(index2d,i,j)
         enddo

          ! ls
         do ivar_gb=1,nmat
          index3d=index3d+1
          index2d=index2d+1
          zone3d_gb(iz_gb)%var(index3d,i,j,k)= &
            zone2d_gb(iz_gb)%var(index2d,i,j)
         enddo
          ! ls normal
         do ivar_gb=1,nmat
          ur=zone2d_gb(iz_gb)%var(index2d+1,i,j)
          uz=zone2d_gb(iz_gb)%var(index2d+2,i,j)
          ux=ur*cos(theta)
          uy=ur*sin(theta)
          zone3d_gb(iz_gb)%var(index3d+1,i,j,k)=ux
          zone3d_gb(iz_gb)%var(index3d+2,i,j,k)=uy
          zone3d_gb(iz_gb)%var(index3d+3,i,j,k)=uz
          index3d=index3d+plot_sdim
          index2d=index2d+SDIM
         enddo
          ! den,configuration tensor,visc,trace
         do ivar_gb=1,nmat*num_state_material+ &
              num_materials_viscoelastic*FORT_NUM_TENSOR_TYPE+ &
              nmat+5*nmat
          index3d=index3d+1
          index2d=index2d+1
          zone3d_gb(iz_gb)%var(index3d,i,j,k)= &
            zone2d_gb(iz_gb)%var(index2d,i,j)
         enddo
         if ((index3d.ne.nwrite3d).or. &
             (index2d.ne.nwrite2d)) then
          print *,"index mismatch in zone_revolve"
          stop
         endif

        enddo ! k

       enddo
       enddo
!      enddo

       close(4)

       ! Zone marker  Value = 299.0
       write(11) ZONEMARKER
       ! Data format
       do i=1,nwrite3d
        write(11) 2
       enddo
       write(11) 0  ! Has passive variables: 0 = no, 1 = yes.
       write(11) 0  ! Has variable sharing 0 = no, 1 = yes.
       write(11) -1 ! Share connectivity list (-1 = no sharing). 
   
       do ivar_gb=1,nwrite3d
        write(11) minval(zone3d_gb(iz_gb)%var(ivar_gb,:,:,:))
        write(11) maxval(zone3d_gb(iz_gb)%var(ivar_gb,:,:,:))
       enddo

       ! order is IMPORTANT.
       do ivar_gb=1,nwrite3d
        do k=klo_plot,khi_plot
        do j=lo_gb(iz_gb,2),hi_gb(iz_gb,2)+1
        do i=lo_gb(iz_gb,1),hi_gb(iz_gb,1)+1
         write(11) zone3d_gb(iz_gb)%var(ivar_gb,i,j,k)
        enddo
        enddo
        enddo
       enddo

       deallocate(zone3d_gb(iz_gb)%var)
       deallocate(zone2d_gb(iz_gb)%var)

       igrid=igrid+1
       if (igrid.gt.grids_per_level_array(ilev+1)-1) then
        ilev=ilev+1
        igrid=0
       endif

      enddo  ! iz_gb

      deallocate(zone3d_gb)
      deallocate(zone2d_gb)
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

      return
      end subroutine zones_revolve

      subroutine project_tangent(vel,nn)
      IMPLICIT NONE
 
      REAL_T vel(SDIM)
      REAL_T nn(SDIM)
      REAL_T velnormal
      INTEGER_T dir
    
      velnormal=zero
      do dir=1,SDIM
       velnormal=velnormal+vel(dir)*nn(dir)
      enddo
      do dir=1,SDIM
       vel(dir)=vel(dir)-velnormal*nn(dir)
      enddo

      return
      end subroutine project_tangent
 

SUBROUTINE FINDMINRAD(N,Z,F,minrad)
  IMPLICIT NONE
  INTEGER_T, INTENT(IN)  ::  N
  INTEGER_T              ::  I
  REAL_T, INTENT(OUT)    ::  minrad
  REAL_T, INTENT(IN), DIMENSION(0:N)    ::  Z
  REAL_T, INTENT(IN), DIMENSION(0:N)    ::  F

  minrad=1.0D+10
  DO I=0,N
   if (F(I).lt.minrad) then
    minrad=F(I)
   endif
  ENDDO
END SUBROUTINE FINDMINRAD

! integral of z*f(z) divided by integral of f(z).
SUBROUTINE FINDZCENTER(N,Z,F,zcenter)
  IMPLICIT NONE
  INTEGER_T, INTENT(IN)  ::  N
  INTEGER_T              ::  I
  REAL_T, INTENT(OUT)    ::  zcenter
  REAL_T, INTENT(IN), DIMENSION(0:N)    ::  Z
  REAL_T, INTENT(IN), DIMENSION(0:N)    ::  F
  REAL_T sumf,sumzf

  sumf=zero
  sumzf=zero
  DO I=0,N
   sumf=sumf+F(I)
   sumzf=sumzf+Z(I)*F(I)
  ENDDO
  zcenter=sumzf/sumf
END SUBROUTINE FINDZCENTER



! 3 points Simpson rule is used to get the integral
SUBROUTINE FINDAMPLITUDE_2D(N,Z,F,WAVELENGTH,AMPLITUDE)
  IMPLICIT NONE
  INTEGER_T, INTENT(IN)  ::  N
  INTEGER_T              ::  I
  REAL_T, INTENT(IN)     ::  WAVELENGTH
  REAL_T, INTENT(OUT)    ::  AMPLITUDE
  REAL_T                 ::  PI, H, SR, SI 
  REAL_T, INTENT(IN), DIMENSION(0:N)    ::  Z
  REAL_T, INTENT(IN), DIMENSION(0:N)    ::  F
  REAL_T, DIMENSION(0:N)   ::  FR, FI

  PI=4.0*ATAN(1.0)
  H=WAVELENGTH/N
  DO I=0,N
    FR(I)=F(I)*COS(-two*H*I*PI/WAVELENGTH)
    FI(I)=F(I)*SIN(-two*H*I*PI/WAVELENGTH)
  END DO
  CALL SIMP(N,H,FR,SR)
  CALL SIMP(N,H,FI,SI)
  AMPLITUDE=2*SQRT(SR*SR+SI*SI)/WAVELENGTH
END SUBROUTINE FINDAMPLITUDE_2D

SUBROUTINE FINDAMPLITUDE_3D(N,M,T,Z,F,WAVELENGTH,NN,AMPLITUDE)
  IMPLICIT NONE
  INTEGER_T, INTENT(IN)  ::  N,M,NN
  INTEGER_T              ::  I,J
  REAL_T, INTENT(IN)     ::  WAVELENGTH
  REAL_T, INTENT(OUT)    ::  AMPLITUDE
  REAL_T                 ::  PI, HT, HZ, SR, SI
  REAL_T, INTENT(IN), DIMENSION(0:N)    ::  T
  REAL_T, INTENT(IN), DIMENSION(0:M)    ::  Z
  REAL_T, INTENT(IN), DIMENSION(0:N,0:M)    ::  F
  REAL_T, DIMENSION(0:M)  ::  TR, TI
  REAL_T, DIMENSION(0:N)  ::  ER, EI, TER, TEI
  
  PI=4.0*ATAN(1.0)
  HT=two*PI/N
  HZ=WAVELENGTH/M
  DO I=0,N
    DO J=0,M
      TR(J)=F(I,J)*COS(-two*HZ*J*PI/WAVELENGTH)
      TI(J)=F(I,J)*SIN(-two*HZ*J*PI/WAVELENGTH)
      CALL SIMP(M,HZ,TR,TER(I))
      CALL SIMP(M,HZ,TI,TEI(I))
    END DO
    ER(I)=TER(I)*COS(-HT*I*NN)-TEI(I)*SIN(-HT*I*NN)
    EI(I)=TER(I)*SIN(-HT*I*NN)+TEI(I)*COS(-HT*I*NN)
  END DO
  CALL SIMP(N,HT,ER,SR)
  CALL SIMP(N,HT,EI,SI)
  AMPLITUDE=4*SQRT(SR*SR+SI*SI)/(two*PI*WAVELENGTH)
END SUBROUTINE FINDAMPLITUDE_3D



 
SUBROUTINE SIMP (N,H,FI,S)
!
! Subroutine for integration over f(x) with the Simpson rule.  FI:
! integrand f(x); H: interval; S: integral.  Copyright (c) Tao Pang 1997.
!
  IMPLICIT NONE
  INTEGER_T, INTENT (IN) :: N
  INTEGER_T :: I
  REAL_T, INTENT (IN) :: H
  REAL_T :: S0,S1,S2
  REAL_T, INTENT (OUT) :: S
  REAL_T, INTENT (IN), DIMENSION (N) :: FI
!
  S  = 0.0
  S0 = 0.0
  S1 = 0.0
  S2 = 0.0
  DO I = 2, N-1, 2
    S1 = S1+FI(I-1)
    S0 = S0+FI(I)
    S2 = S2+FI(I+1)
  END DO
  S = H*(S1+4.0*S0+S2)/3.0
!
! If N is even, add the last slice separately
! 
  IF (MOD(N,2).EQ.0) S = S &
     +H*(5.0*FI(N)+8.0*FI(N-1)-FI(N-2))/12.0
END SUBROUTINE SIMP



      end module navierstokesf90_module

! ZZ are the x values
! FF are the z values
! 2 files:
!   height_integral.txt
!   jet_height.txt
      subroutine FORT_COFLOW( &
       time,js,je,NN,ZZ,FF,dirx,diry,cut_flag)

      use navierstokesf90_module

      IMPLICIT NONE

      REAL_T time
      INTEGER_T dirx,diry,cut_flag
      INTEGER_T js,je,NN
      REAL_T ZZ(0:NN)
      REAL_T FF(0:NN)
      REAL_T FFMASK(0:NN)
      REAL_T major,minor,wavelength,amplitude,minrad
      REAL_T masklo,maskhi,HH,SS
      INTEGER_T j,jleft,jright
      REAL_T dzz,theta


      if (NN.ne.je-js+1) then
       print *,"NN invalid"
       stop
      endif
      if (js.ne.0) then
       print *,"js should be 0 in coflow (if 3d)"
       stop
      endif
      if (NN.lt.1) then
       print *,"NN invalid"
       stop
      endif
      if ((dirx.lt.0).or.(dirx.ge.SDIM)) then
       print *,"dirx invalid"
       stop
      endif
      if ((diry.lt.0).or.(diry.ge.SDIM).or.(diry.eq.dirx)) then
       print *,"diry invalid"
       stop
      endif
      if ((cut_flag.ne.0).and.(cut_flag.ne.1)) then
       print *,"cut_flag invalid"
       stop
      endif

      if ((cut_flag.eq.1).and.(SDIM.eq.2)) then
       major=yblob
       if (abs(major-ZZ(NN)).le.1.0D-10) then
        call FINDAMPLITUDE_2D(NN,ZZ,FF,major,minor)
        call FINDMINRAD(NN,ZZ,FF,minrad)
        wavelength=major
        amplitude=minor
        print *,"TIME,wavelength (major)",time,wavelength
        print *,"TIME,amplitude  (minor)",time,amplitude
        print *,"TIME,minimum radius ",time,minrad
        open(unit=19,file="growth.txt",access="append")
        write(19,*) time, amplitude
        close(19)
        open(unit=29,file="minrad.txt",access="append")
        write(29,*) time, minrad
        close(29)
       endif
      else if ((cut_flag.eq.1).and.(SDIM.eq.3)) then
       open(unit=19,file="height_integral.txt",access="append")
       open(unit=20,file="jet_height.txt",access="append")
       maskhi=one
       masklo=half*(xblob+maskhi)

       ZZ(NN)=probhix-problox
       dzz=(probhix-problox)/NN

         ! fill in gaps in coarse regions
       do j=js+1,je
        if (ZZ(j).eq.zero) then  ! value not initialized in summass
         ZZ(j)=j*dzz
         jleft=j-1
         jright=j+1
         do while (ZZ(jright).eq.zero)
          jright=jright+1
         enddo  
         theta=one/(jright-jleft)
         FF(j)=FF(jleft)*(one-theta)+theta*FF(jright)
        endif
       enddo

       do j=js,je+1 
        if (ZZ(j).lt.masklo) then
         FFMASK(j)=zero
        else if (ZZ(j).gt.maskhi) then
         FFMASK(j)=zero
        else 
         FFMASK(j)=FF(j)
        endif
       enddo
       HH=ZZ(NN)-ZZ(0)
       call SIMP(NN,HH,FFMASK,SS)
       if ((SS.lt.zero).or.(HH.le.zero)) then
        print *,"SS or HH invalid"
        print *,"SS= ",SS
        print *,"HH= ",HH
        print *,"NN= ",NN
        do j=js,je+1
         print *,"j,ZZ ",j,ZZ(j)
        enddo
        stop
       endif
       SS=SS/HH
       write(19,*) time,SS 
       print *,"TIME,AVERAGE HEIGHT ",time,SS
       do j=js,je+1 
        write(20,*) ZZ(j),FF(j)
       enddo
       write(20,*) ' '

       close(19) 
       close(20) 
      else if (cut_flag.eq.0) then
       ! do nothing
      else
       print *,"cut_flag invalid"
       stop
      endif


      return
      end subroutine FORT_COFLOW


      subroutine FORT_ISOGRID( &
       tid, &
       visual_tessellate_vfrac, &
       recon,DIMS(recon), &
       xlo,dx, &
       mask,DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       level,gridno,nmat)

      use probf90_module
      use navierstokesf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T tid
      INTEGER_T visual_tessellate_vfrac
      INTEGER_T nmat
      INTEGER_T tilelo(SDIM), tilehi(SDIM)
      INTEGER_T fablo(SDIM), fabhi(SDIM)
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T bfact
      INTEGER_T DIMDEC(recon)
      INTEGER_T DIMDEC(mask)
      INTEGER_T level,gridno
      INTEGER_T im
      REAL_T recon(DIMV(recon),nmat*ngeom_recon) ! vof,refcen,order,slope,int
      REAL_T mask(DIMV(mask))
      REAL_T valu
      REAL_T xlo(SDIM),dx(SDIM)

      character*2 matstr

      character*11 namestr11
      character*15 cennamestr15
      character*3 levstr
      character*5 gridstr
      character*19 filename19
      character*23 cenfilename23

      real(8), dimension(:,:), allocatable :: Node  ! dir,index
      integer(4), dimension(:,:), allocatable :: IntElem  ! 1 or 2, index
      integer(4) :: NumNodes
      integer(4) :: NumIntElems
      integer(4) :: CurNodes
      integer(4) :: CurIntElems
      INTEGER_T imaxtri,itri,curtri,ipass
      REAL_T geom(SDIM,200)
      REAL_T lnode(0:1,0:1,0:1)
      REAL_T xnode(0:1,0:1,0:1,SDIM)
      INTEGER_T i,j,k,ii,jj,kk,dir,ISUM,itemp,jtemp
      INTEGER_T N8(8)
      REAL_T gridval(8)
      REAL_T gridx(8)
      REAL_T gridy(8)
      REAL_T gridz(8)

      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      REAL_T xtarget(SDIM)
      REAL_T nn(SDIM)
      REAL_T intercept
      REAL_T xref(SDIM)
      REAL_T xrefT(SDIM)
      INTEGER_T nparticles
      REAL_T vfrac,vfrac_side
      INTEGER_T iside,jside,kside
      REAL_T xoutput(SDIM)
      REAL_T xoutputT(SDIM)
      INTEGER_T kklo,kkhi,nodehi
      INTEGER_T DIMDEC(plt)
      REAL_T, dimension(D_DECL(:,:,:),:), allocatable :: reconlocal
      INTEGER_T nmax
      INTEGER_T vofcomp
      REAL_T vfrac_sum_solid
      REAL_T mofdata(nmat*ngeom_recon)
      INTEGER_T mask_local
      INTEGER_T vcdeb,imdeb

      nhalf=3
      nmax=POLYGON_LIST_MAX ! in: ISOGRID

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid151"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat bust"
       stop
      endif 

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
      call box_to_dim( &
        DIMS(plt), &
        growlo,growhi)
      allocate(reconlocal(DIMV(plt),nmat*ngeom_recon))

      call checkbound(fablo,fabhi,DIMS(plt),1,-1,411)
      call checkbound(fablo,fabhi,DIMS(recon),1,-1,411)
      call checkbound(fablo,fabhi,DIMS(mask),0,-1,411)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)
       do dir=1,nmat*ngeom_recon
        mofdata(dir)=recon(D_DECL(i,j,k),dir)
       enddo

       if (visual_tessellate_vfrac.eq.1) then
         ! before (mofdata): fluids tessellate
         ! after  (mofdata): fluids and solids tessellate

        vfrac_sum_solid=zero
        do im=1,nmat
         if (is_rigid(nmat,im).eq.1) then
          vofcomp=(im-1)*ngeom_recon+1
          vfrac_sum_solid=vfrac_sum_solid+mofdata(vofcomp)
         else if (is_rigid(nmat,im).eq.0) then
          ! do nothing
         else
          print *,"is_rigid invalid"
          stop
         endif
        enddo ! im=1..nmat

        if ((vfrac_sum_solid.gt.VOFTOL).and. &
            (vfrac_sum_solid.le.1.1)) then
         ! before (mofdata): fluids tessellate
         ! after  (mofdata): fluids and solids tessellate
         call multi_get_volume_tessellate( &
          bfact,dx,xsten,nhalf, &
          mofdata, &
          geom_xtetlist(1,1,1,tid+1), &
          nmax, &
          nmax, &
          nmat,SDIM,7)

        else if ((vfrac_sum_solid.ge.-VOFTOL).and. &
                 (vfrac_sum_solid.le.VOFTOL)) then
         ! do nothing
        else
         print *,"vfrac_sum_solid invalid"
         stop
        endif
        
       else if (visual_tessellate_vfrac.eq.0) then
        ! do nothing
       else
        print *,"visual_tessellate_vfrac invalid"
        stop
       endif

       do dir=1,nmat*ngeom_recon
        reconlocal(D_DECL(i,j,k),dir)=mofdata(dir)
       enddo

      enddo ! k
      enddo ! j
      enddo ! i

      do im=1,nmat

       if ((im.lt.1).or. &
           (im.gt.99).or. &
           (im.gt.nmat)) then
        print *,"im out of range"
        stop
       endif

       vofcomp=(im-1)*ngeom_recon+1

       write(matstr,'(I2)') im
       do i=1,2
        if (matstr(i:i).eq.' ') then
         matstr(i:i)='0'
        endif
       enddo

       imaxtri=200

       valu=zero

       write(namestr11,'(A7,A2,A2)') 'tempmat',matstr,'ls'
       write(cennamestr15,'(A10,A2,A3)') 'temprefcen',matstr,'pos'

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
       write(filename19,'(A11,A3,A5)') namestr11,levstr,gridstr
       print *,"filename19 ",filename19
       write(cenfilename23,'(A15,A3,A5)') cennamestr15,levstr,gridstr
       print *,"cenfilename23 ",cenfilename23

       NumNodes=0
       NumIntElems=0

       open(unit=12,file=cenfilename23)
       nparticles=0

       do ipass=0,1

        if (ipass.eq.1) then
         write(12,*)nparticles
        endif

        CurNodes=0
        CurIntElems=0
        if (ipass.eq.1) then
         allocate(Node(SDIM,NumNodes))
         allocate(IntElem(SDIM,NumIntElems))
        endif

        if (SDIM.eq.3) then
         kklo=0
         kkhi=1
         nodehi=8
        else if (SDIM.eq.2) then
         kklo=0
         kkhi=0
         nodehi=4
        else
         print *,"dimension bust"
         stop
        endif

        call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
  
        do i=growlo(1),growhi(1)
        do j=growlo(2),growhi(2)
        do k=growlo(3),growhi(3)

         mask_local=NINT(mask(D_DECL(i,j,k)))

         if (mask_local.eq.1) then
          call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

          do dir=1,nmat*ngeom_recon
           mofdata(dir)=reconlocal(D_DECL(i,j,k),dir)
          enddo

           ! vfrac,x,rank,slope,int
          do dir=1,SDIM
           nn(dir)=mofdata(vofcomp+SDIM+1+dir) ! slope
          enddo
          intercept=mofdata(vofcomp+2*SDIM+2) ! intercept

          vfrac=mofdata(vofcomp) ! volume fraction

          if (1.eq.0) then
           imdeb=3
           vcdeb=(imdeb-1)*ngeom_recon+1
           if (abs(mofdata(vcdeb)-0.3942).le.0.001) then
            print *,"------sanity check isogrid"
            print *,"i,j,k,im = ",i,j,k,im
            do vcdeb=1,nmat*ngeom_recon
             print *,"vcdeb,mofdata(vcdeb) ",vcdeb,mofdata(vcdeb)
            enddo
            print *,"------end sanity check isogrid"
           endif
          endif

          itri=0
          do ii=0,1
          do jj=0,1
          do kk=kklo,kkhi

           do dir=1,SDIM
            xtarget(dir)=xsten(-1,dir)
           enddo
           if (ii.eq.1) then
            dir=1 
            xtarget(dir)=xsten(1,dir)
           endif
           if (jj.eq.1) then
            dir=2 
            xtarget(dir)=xsten(1,dir)
           endif
           if (kk.eq.1) then
            dir=SDIM 
            xtarget(dir)=xsten(1,dir)
           endif
           do dir=1,SDIM
            xnode(ii,jj,kk,dir)=xtarget(dir)
           enddo

           if (ii.eq.0) then
            iside=i-1
           else if (ii.eq.1) then
            iside=i+1
           else
            print *,"ii invalid"
            stop
           endif
           if (jj.eq.0) then
            jside=j-1
           else if (jj.eq.1) then
            jside=j+1
           else
            print *,"jj invalid"
            stop
           endif
           if (kk.eq.0) then
            kside=k-1
           else if (kk.eq.1) then
            kside=k+1
           else
            print *,"kk invalid"
            stop
           endif
           
           vfrac_side=reconlocal(D_DECL(iside,jside,kside),vofcomp)

           if ((vfrac.ge.one-VOFTOL).and. &
               (vfrac_side.le.VOFTOL)) then
            lnode(ii,jj,kk)=-VOFTOL
           else if ((vfrac.le.VOFTOL).and. &
                    (vfrac_side.ge.one-VOFTOL)) then
            lnode(ii,jj,kk)=VOFTOL
           else if (vfrac.ge.one-VOFTOL) then
            lnode(ii,jj,kk)=one
           else if (vfrac.le.VOFTOL) then
            lnode(ii,jj,kk)=-one
           else
            call distfunc(bfact,dx,xsten,nhalf, &
             intercept,nn,xtarget,lnode(ii,jj,kk),SDIM)
           endif
          enddo
          enddo
          enddo ! ii,jj,kk

          gridval(1)=lnode(0,0,0)
          gridval(2)=lnode(1,0,0)
          gridval(3)=lnode(1,1,0)
          gridval(4)=lnode(0,1,0)
          gridval(5)=lnode(0,0,kkhi)
          gridval(6)=lnode(1,0,kkhi)
          gridval(7)=lnode(1,1,kkhi)
          gridval(8)=lnode(0,1,kkhi)

          ISUM=0
          do itemp=1,nodehi
           if (gridval(itemp).ge.valu) then
            N8(itemp)=1
           else 
            N8(itemp)=0
           endif
           ISUM=ISUM+N8(itemp)
          enddo

          if ((ISUM.eq.0).or.(ISUM.eq.nodehi)) then
           goto 999
          endif

          if (ipass.eq.0) then
           nparticles=nparticles+1
          else if (ipass.eq.1) then
           do dir=1,SDIM
            xref(dir)=xsten(0,dir)+mofdata(vofcomp+dir)
            xrefT(dir)=xref(dir)
           enddo
           if (visual_RT_transform.eq.1) then
            call RT_transform(xref,xrefT)
           endif
           if (SDIM.eq.3) then
            write(12,*) xrefT(1),xrefT(2),xrefT(SDIM)
           else if (SDIM.eq.2) then
            write(12,*) xrefT(1),xrefT(2)
           else
            print *,"dimension bust"
            stop
           endif
          else
           print *,"ipass invalid"
           stop
          endif
         
          dir=1
          gridx(1)=xnode(0,0,0,dir)
          gridx(2)=xnode(1,0,0,dir)
          gridx(3)=xnode(1,1,0,dir)
          gridx(4)=xnode(0,1,0,dir)
          gridx(5)=xnode(0,0,kkhi,dir)
          gridx(6)=xnode(1,0,kkhi,dir)
          gridx(7)=xnode(1,1,kkhi,dir)
          gridx(8)=xnode(0,1,kkhi,dir)
          dir=2
          gridy(1)=xnode(0,0,0,dir)
          gridy(2)=xnode(1,0,0,dir)
          gridy(3)=xnode(1,1,0,dir)
          gridy(4)=xnode(0,1,0,dir)
          gridy(5)=xnode(0,0,kkhi,dir)
          gridy(6)=xnode(1,0,kkhi,dir)
          gridy(7)=xnode(1,1,kkhi,dir)
          gridy(8)=xnode(0,1,kkhi,dir)
          dir=SDIM
          gridz(1)=xnode(0,0,0,dir)
          gridz(2)=xnode(1,0,0,dir)
          gridz(3)=xnode(1,1,0,dir)
          gridz(4)=xnode(0,1,0,dir)
          gridz(5)=xnode(0,0,kkhi,dir)
          gridz(6)=xnode(1,0,kkhi,dir)
          gridz(7)=xnode(1,1,kkhi,dir)
          gridz(8)=xnode(0,1,kkhi,dir)

          if (SDIM.eq.2) then
           call polysegment(gridx,gridy,gridz,gridval,valu,geom, &
            itri,1,2,3,3,imaxtri)
           call polysegment(gridx,gridy,gridz,gridval,valu,geom, &
            itri,1,4,3,3,imaxtri)
          else if (SDIM.eq.3) then
           call polytri(gridx,gridy,gridz,gridval,valu,geom, &
            itri,1,3,4,8,imaxtri)
           call polytri(gridx,gridy,gridz,gridval,valu,geom, &
            itri,1,3,7,8,imaxtri)
           call polytri(gridx,gridy,gridz,gridval,valu,geom, &
            itri,1,5,7,8,imaxtri)
           call polytri(gridx,gridy,gridz,gridval,valu,geom, &
            itri,1,7,2,3,imaxtri)
           call polytri(gridx,gridy,gridz,gridval,valu,geom, &
            itri,1,7,2,5,imaxtri)
           call polytri(gridx,gridy,gridz,gridval,valu,geom, &
            itri,6,7,2,5,imaxtri)
          else
           print *,"dimension bust"
           stop
          endif

999       continue
 
          if (ipass.eq.0) then
           NumNodes=NumNodes+itri
           NumIntElems=NumIntElems+itri/SDIM
          else if (ipass.eq.1) then
           curtri=0
           do itemp=1,itri/SDIM
            CurIntElems=CurIntElems+1
            if (CurIntElems.gt.NumIntElems) then
             print *,"CurIntElems invalid"
             stop
            endif
            do jtemp=1,SDIM
             CurNodes=CurNodes+1
             if (CurNodes.gt.NumNodes) then
              print *,"CurNodes invalid"
              stop
             endif
             curtri=curtri+1
             if (curtri.gt.itri) then
              print *,"curtri invalid"
              stop
             endif
             do dir=1,SDIM
              Node(dir,CurNodes)=geom(dir,curtri)
             enddo
             IntElem(jtemp,CurIntElems)=CurNodes
            enddo
           enddo  
          else
           print *,"ipass invalid"
           stop
          endif
  
         else if (mask_local.eq.0) then
          ! do nothing
         else
          print *,"mask_local invalid"
          stop 
         endif

        enddo ! k
        enddo ! j
        enddo ! i

       enddo ! ipass

       close(12)

       open(unit=11,file=filename19)
       write(11,*) NumNodes
       write(11,*) NumIntElems
       do i=1,NumNodes 
        do dir=1,SDIM
         xoutput(dir)=Node(dir,i)
        enddo
        do dir=1,SDIM
         xoutputT(dir)=xoutput(dir)
        enddo
        if (visual_RT_transform.eq.1) then
         call RT_transform(xoutput,xoutputT)
        endif
        if (SDIM.eq.3) then
         write(11,*) xoutputT(1),xoutputT(2),xoutputT(SDIM)
        else if (SDIM.eq.2) then
         write(11,*) xoutputT(1),xoutputT(2)
        else
         print *,"dimension bust"
         stop
        endif
       enddo
       do i=1,NumIntElems
        if (SDIM.eq.3) then
         write(11,*) IntElem(1,i),IntElem(2,i),IntElem(SDIM,i)
        else if (SDIM.eq.2) then
         write(11,*) IntElem(1,i),IntElem(2,i)
        else
         print *,"dimension bust"
         stop
        endif
       enddo
       close(11)

       deallocate(Node)
       deallocate(IntElem)

      enddo ! im=1..nmat

      deallocate(reconlocal)
 
      return
      end subroutine FORT_ISOGRID


      subroutine FORT_COMBINETRIANGLES( &
       grids_per_level,finest_level,nsteps,im, &
       arrdim,time,plotint)
      use probcommon_module

      IMPLICIT NONE


      REAL_T time
      INTEGER_T    plotint,strandid
      INTEGER_T    arrdim,finest_level,nsteps,im
      INTEGER_T    grids_per_level(arrdim)

      character*11 namestr11
      character*7 newnamestr7

      character*15 cennamestr15
      character*11 newcennamestr11

      character*3 levstr
      character*5 gridstr

      character*19 filename19
      character*23 cenfilename23

      character*2 matstr
      character*6 stepstr

      character*17 newfilename17
      character*21 newcenfilename21

      character*80 rmcommand_mat
      character*80 rmcommand_refcen

      real(8), dimension(:,:), allocatable :: Node
      integer(4), dimension(:,:), allocatable :: IntElem
      integer(4) :: NumNodes
      integer(4) :: NumIntElems
      integer(4) :: PartNumNodes
      integer(4) :: PartNumIntElems
      integer(4) :: CurNumNodes
      integer(4) :: CurNumIntElems
      INTEGER_T i,dir
      INTEGER_T ilev,igrid,ipass
      REAL_T xref(SDIM)
      INTEGER_T nparticles,Part_nparticles
      INTEGER_T alloc_flag
      INTEGER_T sysret

      alloc_flag=0

      if ((im.lt.1).or.(im.gt.99).or. &
          (im.gt.num_materials)) then
       print *,"im out of range"
       stop
      endif

      write(matstr,'(I2)') im
      do i=1,2
       if (matstr(i:i).eq.' ') then
        matstr(i:i)='0'
       endif
      enddo
      
      write(namestr11,'(A7,A2,A2)') 'tempmat',matstr,'ls'
      write(newnamestr7,'(A3,A2,A2)') 'mat',matstr,'ls'
      write(cennamestr15,'(A10,A2,A3)') 'temprefcen',matstr,'pos'
      write(newcennamestr11,'(A6,A2,A3)') 'refcen',matstr,'pos'

      NumNodes=0
      NumIntElems=0
      nparticles=0

      if (arrdim.ne.finest_level+1) then
       print *,"arrdim invalid"
       stop
      endif

      do ipass=0,1
       CurNumNodes=0
       CurNumIntElems=0

       if (ipass.eq.1) then

        allocate(Node(SDIM,NumNodes))
        allocate(IntElem(SDIM,NumIntElems))
        alloc_flag=alloc_flag+1

        write(stepstr,'(I6)') nsteps
        do i=1,6
         if (stepstr(i:i).eq.' ') then
          stepstr(i:i)='0'
         endif
        enddo
        write(newfilename17,'(A7,A6,A4)') newnamestr7,stepstr,'.tec'
        print *,"newfilename17 ",newfilename17
        open(unit=11,file=newfilename17)

        if (SDIM.eq.3) then
         write(11,*) 'TITLE = "3D surface" '
         write(11,*) 'VARIABLES = "X", "Y", "Z" '
        else if (SDIM.eq.2) then
         write(11,*) 'TITLE = "2D surface" '
         write(11,*) 'VARIABLES = "X", "Y" '
        else
         print *,"dimension bust"
         stop
        endif

        write(11,'(A23,I14,A5,I14,A21)')  &
          'ZONE T="TRIANGLES", N= ',NumNodes,  &
          ', E= ',NumIntElems,', DATAPACKING=POINT, '

        if (plotint.le.0) then
         strandid=1
        else
         strandid=(nsteps/plotint)+1
        endif

        if (SDIM.eq.3) then
         write(11,'(A33,D25.16,A10,I10)')  &
          'ZONETYPE=FETRIANGLE SOLUTIONTIME=',time, &
          " STRANDID=",strandid
        else if (SDIM.eq.2) then
         write(11,'(A32,D25.16,A10,I10)')  &
          'ZONETYPE=FELINESEG SOLUTIONTIME=',time, &
          " STRANDID=",strandid
        else
         print *,"dimension bust"
         stop
        endif

        write(newcenfilename21,'(A11,A6,A4)') newcennamestr11,stepstr,'.tec'
        print *,"newcenfilename21 ",newcenfilename21
        open(unit=12,file=newcenfilename21)

        if (SDIM.eq.3) then
         write(12,*) 'TITLE = "3D moments" '
         write(12,*) 'VARIABLES = "X", "Y", "Z" '
        else if (SDIM.eq.2) then
         write(12,*) 'TITLE = "2D moments" '
         write(12,*) 'VARIABLES = "X", "Y" '
        else
         print *,"dimension bust"
         stop
        endif

        if (plotint.le.0) then
         strandid=1
        else
         strandid=(nsteps/plotint)+1
        endif

        write(12,'(A19,I14,A26,D25.16,A10,I10)') & 
          'ZONE F="POINT", I= ', nparticles,  &
          ', J=1, K=1, SOLUTIONTIME= ',time,' STRANDID=',strandid

       endif  !ipass=1

       do ilev=0,finest_level
       do igrid=0,grids_per_level(ilev+1)-1
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

          ! tempmat
         write(filename19,'(A11,A3,A5)') namestr11,levstr,gridstr
         print *,"filename19 ",filename19
         open(unit=4,file=filename19)

         read(4,*) PartNumNodes
         read(4,*) PartNumIntElems

         write(cenfilename23,'(A15,A3,A5)') cennamestr15,levstr,gridstr
         print *,"cenfilename23 ",cenfilename23
         open(unit=5,file=cenfilename23)

         read(5,*) Part_nparticles

         if (ipass.eq.0) then
          NumNodes=NumNodes+PartNumNodes
          NumIntElems=NumIntElems+PartNumIntElems
          nparticles=nparticles+Part_nparticles
         else if (ipass.eq.1) then
          do i=1,PartNumNodes 
           if (SDIM.eq.3) then
            read(4,*) Node(1,CurNumNodes+i),Node(2,CurNumNodes+i), &
             Node(SDIM,CurNumNodes+i)
           else if (SDIM.eq.2) then
            read(4,*) Node(1,CurNumNodes+i),Node(2,CurNumNodes+i)
           else
            print *,"dimension bust"
            stop
           endif
          enddo
          do i=1,PartNumIntElems
           if (SDIM.eq.3) then
            read(4,*) IntElem(1,CurNumIntElems+i), &
                     IntElem(2,CurNumIntElems+i), &
                     IntElem(SDIM,CurNumIntElems+i)
           else if (SDIM.eq.2) then
            read(4,*) IntElem(1,CurNumIntElems+i), &
                     IntElem(2,CurNumIntElems+i)
           else
            print *,"dimension bust"
            stop
           endif

           do dir=1,SDIM
            IntElem(dir,CurNumIntElems+i)=IntElem(dir,CurNumIntElems+i)+ &
             CurNumNodes
           enddo
          enddo

          CurNumIntElems=CurNumIntElems+PartNumIntElems
          if (CurNumIntElems.gt.NumIntElems) then
           print *,"CurNumIntElems invalid"
           stop
          endif
          CurNumNodes=CurNumNodes+PartNumNodes
          if (CurNumNodes.gt.NumNodes) then
           print *,"CurNumNodes invalid"
           stop
          endif

          do i=1,Part_nparticles
           if (SDIM.eq.3) then
            read(5,*) xref(1),xref(2),xref(SDIM)
            write(12,*) xref(1),xref(2),xref(SDIM)
           else if (SDIM.eq.2) then
            read(5,*) xref(1),xref(2)
            write(12,*) xref(1),xref(2)
           else
            print *,"dimension bust"
            stop
           endif
          enddo

         else
          print *,"ipass invalid"
          stop
         endif
         close(4)
         close(5)
       enddo
       enddo
      enddo
! ipass
      if (SDIM.eq.3) then
       do i=1,NumNodes
        write(11,*) Node(1,i),Node(2,i),Node(SDIM,i)
       enddo
       do i=1,NumIntElems
        write(11,*) IntElem(1,i),IntElem(2,i),IntElem(SDIM,i)
       enddo
      else if (SDIM.eq.2) then
       do i=1,NumNodes
        write(11,*) Node(1,i),Node(2,i)
       enddo
       do i=1,NumIntElems
        write(11,*) IntElem(1,i),IntElem(2,i)
       enddo
      else
       print *,"dimension bust"
       stop
      endif
 
      close(11)
      close(12)

      deallocate(Node)
      deallocate(IntElem)
      alloc_flag=alloc_flag-1

      if (alloc_flag.gt.0) then
       print *,"alloc_flag bust"
       stop
      endif

      sysret=0

      if (im.eq.num_materials) then

       rmcommand_mat='rm tempmat*'
       rmcommand_refcen='rm temprefcen*'
       print *,"issuing command ",rmcommand_mat
       print *,"and issuing command ",rmcommand_refcen

#ifdef PGIFORTRAN
       call system(rmcommand_mat)
       call system(rmcommand_refcen)
#else
       call execute_command_line(rmcommand_mat,exitstat=sysret)
       call execute_command_line(rmcommand_refcen,exitstat=sysret)
#endif
      endif

      if (sysret.ne.0) then
       print *,"execute_command_line has sysret=",sysret
       stop
      endif
 
      return
      end




      subroutine FORT_ISOGRIDSINGLE( &
       ls,DIMS(ls), &
       xlo,dx, &
       mask,DIMS(mask), &
       lo,hi, &
       level,gridno)

      use navierstokesf90_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T    lo(SDIM), hi(SDIM)
      INTEGER_T    growlo(3),growhi(3)
      INTEGER_T    DIMDEC(ls)
      INTEGER_T    DIMDEC(mask)
      INTEGER_T    level,gridno
      REAL_T ls(DIMV(ls))
      REAL_T mask(DIMV(mask))
      REAL_T valu
      REAL_T xlo(SDIM),dx(SDIM)

      character*3 levstr
      character*5 gridstr
      character*18 filename18

      real(8), dimension(:,:), allocatable :: Node  ! dir,index
      integer(4), dimension(:,:), allocatable :: IntElem  ! 1 or 2, index
      integer(4) :: NumNodes
      integer(4) :: NumIntElems
      integer(4) :: CurNodes
      integer(4) :: CurIntElems
      INTEGER_T imaxtri,itri,curtri,ipass
      REAL_T geom(SDIM,200)
      REAL_T lnode(0:1,0:1,0:1)
      REAL_T xnode(0:1,0:1,0:1,SDIM)
      INTEGER_T i,j,k,ii,jj,kk,dir,ISUM,itemp,jtemp
      INTEGER_T N8(8)
      REAL_T gridval(8)
      REAL_T gridx(8)
      REAL_T gridy(8)
      REAL_T gridz(8)

      REAL_T xtarget(SDIM)
      REAL_T xoutput(SDIM)
      REAL_T xoutputT(SDIM)
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf
      INTEGER_T kklo,kkhi,nodehi
      INTEGER_T bfact

      nhalf=3
      bfact=2
      print *,"----------WARNING: bfact override to be 2------------"

      imaxtri=200
      call checkbound(lo,hi,DIMS(ls),1,-1,411)
      call checkbound(lo,hi,DIMS(mask),0,-1,411)

      valu=zero

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
       ! isogridsingle
      write(filename18,'(A10,A3,A5)') 'templssing',levstr,gridstr
      print *,"filename18 ",filename18

      NumNodes=0
      NumIntElems=0
      do ipass=0,1
       CurNodes=0
       CurIntElems=0
       if (ipass.eq.1) then
        allocate(Node(SDIM,NumNodes))
        allocate(IntElem(SDIM,NumIntElems))
       endif

       if (SDIM.eq.3) then
        kklo=0
        kkhi=1
        nodehi=8
       else if (SDIM.eq.2) then
        kklo=0
        kkhi=0
        nodehi=4
       else
        print *,"dimension bust"
        stop
       endif

       call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        if (mask(D_DECL(i,j,k)).eq.one) then

         call gridsten(xsten,xlo,i,j,k,lo,bfact,dx,nhalf)

         itri=0
         do ii=0,1
         do jj=0,1
         do kk=kklo,kkhi

          do dir=1,SDIM
           xtarget(dir)=xsten(-1,dir)
          enddo
          if (ii.eq.1) then
           dir=1 
           xtarget(dir)=xsten(1,dir)
          endif
          if (jj.eq.1) then
           dir=2 
           xtarget(dir)=xsten(1,dir)
          endif
          if (kk.eq.1) then
           dir=SDIM 
           xtarget(dir)=xsten(1,dir)
          endif
          do dir=1,SDIM
           xnode(ii,jj,kk,dir)=xtarget(dir)
          enddo

          lnode(ii,jj,kk)=( &
            ls(D_DECL(i+ii-1,j+jj-1,k+kk-1))+ &
            ls(D_DECL(i+ii,j+jj-1,k+kk-1))+ &
            ls(D_DECL(i+ii-1,j+jj,k+kk-1))+ &
            ls(D_DECL(i+ii,j+jj,k+kk-1))+ &
            ls(D_DECL(i+ii-1,j+jj-1,k+kk))+ &
            ls(D_DECL(i+ii,j+jj-1,k+kk))+ &
            ls(D_DECL(i+ii-1,j+jj,k+kk))+ &
            ls(D_DECL(i+ii,j+jj,k+kk)))/eight

         enddo
         enddo
         enddo

         gridval(1)=lnode(0,0,0)
         gridval(2)=lnode(1,0,0)
         gridval(3)=lnode(1,1,0)
         gridval(4)=lnode(0,1,0)
         gridval(5)=lnode(0,0,kkhi)
         gridval(6)=lnode(1,0,kkhi)
         gridval(7)=lnode(1,1,kkhi)
         gridval(8)=lnode(0,1,kkhi)

         ISUM=0
         do itemp=1,nodehi
          if (gridval(itemp).ge.valu) then
           N8(itemp)=1
          else 
           N8(itemp)=0
          endif
          ISUM=ISUM+N8(itemp)
         enddo

         if ((ISUM.eq.0).or.(ISUM.eq.nodehi)) then
          goto 999
         endif

         dir=1
         gridx(1)=xnode(0,0,0,dir)
         gridx(2)=xnode(1,0,0,dir)
         gridx(3)=xnode(1,1,0,dir)
         gridx(4)=xnode(0,1,0,dir)
         gridx(5)=xnode(0,0,kkhi,dir)
         gridx(6)=xnode(1,0,kkhi,dir)
         gridx(7)=xnode(1,1,kkhi,dir)
         gridx(8)=xnode(0,1,kkhi,dir)
         dir=2
         gridy(1)=xnode(0,0,0,dir)
         gridy(2)=xnode(1,0,0,dir)
         gridy(3)=xnode(1,1,0,dir)
         gridy(4)=xnode(0,1,0,dir)
         gridy(5)=xnode(0,0,kkhi,dir)
         gridy(6)=xnode(1,0,kkhi,dir)
         gridy(7)=xnode(1,1,kkhi,dir)
         gridy(8)=xnode(0,1,kkhi,dir)
         dir=SDIM
         gridz(1)=xnode(0,0,0,dir)
         gridz(2)=xnode(1,0,0,dir)
         gridz(3)=xnode(1,1,0,dir)
         gridz(4)=xnode(0,1,0,dir)
         gridz(5)=xnode(0,0,kkhi,dir)
         gridz(6)=xnode(1,0,kkhi,dir)
         gridz(7)=xnode(1,1,kkhi,dir)
         gridz(8)=xnode(0,1,kkhi,dir)

         if (SDIM.eq.2) then
          call polysegment(gridx,gridy,gridz,gridval,valu,geom, &
           itri,1,2,3,3,imaxtri)
          call polysegment(gridx,gridy,gridz,gridval,valu,geom, &
           itri,1,4,3,3,imaxtri)
         else if (SDIM.eq.3) then
          call polytri(gridx,gridy,gridz,gridval,valu,geom, &
           itri,1,3,4,8,imaxtri)
          call polytri(gridx,gridy,gridz,gridval,valu,geom, &
           itri,1,3,7,8,imaxtri)
          call polytri(gridx,gridy,gridz,gridval,valu,geom, &
           itri,1,5,7,8,imaxtri)
          call polytri(gridx,gridy,gridz,gridval,valu,geom, &
           itri,1,7,2,3,imaxtri)
          call polytri(gridx,gridy,gridz,gridval,valu,geom, &
           itri,1,7,2,5,imaxtri)
          call polytri(gridx,gridy,gridz,gridval,valu,geom, &
           itri,6,7,2,5,imaxtri)
         else
          print *,"dimension bust"
          stop
         endif

999      continue
 
         if (ipass.eq.0) then
          NumNodes=NumNodes+itri
          NumIntElems=NumIntElems+itri/SDIM
         else if (ipass.eq.1) then
          curtri=0
          do itemp=1,itri/SDIM
           CurIntElems=CurIntElems+1
           if (CurIntElems.gt.NumIntElems) then
            print *,"CurIntElems invalid"
            stop
           endif
           do jtemp=1,SDIM
            CurNodes=CurNodes+1
            if (CurNodes.gt.NumNodes) then
             print *,"CurNodes invalid"
             stop
            endif
            curtri=curtri+1
            if (curtri.gt.itri) then
             print *,"curtri invalid"
             stop
            endif
            do dir=1,SDIM
             Node(dir,CurNodes)=geom(dir,curtri)
            enddo
            IntElem(jtemp,CurIntElems)=CurNodes
           enddo
          enddo  
         else
          print *,"ipass invalid"
          stop
         endif
  
        endif ! mask=1
       enddo
       enddo
       enddo
      enddo ! ipass

      open(unit=11,file=filename18)
      write(11,*) NumNodes
      write(11,*) NumIntElems
      do i=1,NumNodes 
       do dir=1,SDIM
        xoutput(dir)=Node(dir,i)
        xoutputT(dir)=xoutput(dir)
       enddo
       if (visual_RT_transform.eq.1) then
        call RT_transform(xoutput,xoutputT)
       endif
       if (SDIM.eq.3) then
        write(11,*) xoutputT(1),xoutputT(2),xoutputT(SDIM)
       else if (SDIM.eq.2) then
        write(11,*) xoutputT(1),xoutputT(2)
       else
        print *,"dimension bust"
        stop
       endif
      enddo
      do i=1,NumIntElems
       if (SDIM.eq.3) then
        write(11,*) IntElem(1,i),IntElem(2,i),IntElem(SDIM,i)
       else if (SDIM.eq.2) then
        write(11,*) IntElem(1,i),IntElem(2,i)
       else
        print *,"dimension bust"
        stop
       endif
      enddo
      close(11)

      deallocate(Node)
      deallocate(IntElem)
 
      return
      end subroutine FORT_ISOGRIDSINGLE


      subroutine FORT_COMBINETRIANGLESSINGLE( &
       grids_per_level,finest_level,nsteps,arrdim)
      IMPLICIT NONE

      INTEGER_T    arrdim,finest_level,nsteps
      INTEGER_T    grids_per_level(arrdim)

      character*3 levstr
      character*5 gridstr
      character*18 filename18

      character*6 stepstr
      character*16 newfilename16

      character*14 rmcommand14

      real(8), dimension(:,:), allocatable :: Node
      integer(4), dimension(:,:), allocatable :: IntElem
      integer(4) :: NumNodes
      integer(4) :: NumIntElems
      integer(4) :: PartNumNodes
      integer(4) :: PartNumIntElems
      integer(4) :: CurNumNodes
      integer(4) :: CurNumIntElems
      INTEGER_T i,dir
      INTEGER_T ilev,igrid,ipass
      INTEGER_T sysret

      NumNodes=0
      NumIntElems=0

      if (arrdim.ne.finest_level+1) then
       print *,"arrdim invalid"
       stop
      endif

      do ipass=0,1
       CurNumNodes=0
       CurNumIntElems=0

       if (ipass.eq.1) then

        allocate(Node(SDIM,NumNodes))
        allocate(IntElem(SDIM,NumIntElems))
        write(stepstr,'(I6)') nsteps
        do i=1,6
         if (stepstr(i:i).eq.' ') then
          stepstr(i:i)='0'
         endif
        enddo
        write(newfilename16,'(A6,A6,A4)') 'lssing',stepstr,'.tec'
        print *,"newfilename16 ",newfilename16
        open(unit=11,file=newfilename16)

        if (SDIM.eq.3) then
         write(11,*) 'TITLE = "3D surface" '
         write(11,*) 'VARIABLES = "X", "Y", "Z" '
        else if (SDIM.eq.2) then
         write(11,*) 'TITLE = "2D surface" '
         write(11,*) 'VARIABLES = "X", "Y" '
        else
         print *,"dimension bust"
         stop
        endif

        write(11,'(A23,I14,A5,I14,A21)')  &
          'ZONE T="TRIANGLES", N= ',NumNodes,  &
          ', E= ',NumIntElems,', DATAPACKING=POINT, '
      
        if (SDIM.eq.3) then
         write(11,'(A19)') 'ZONETYPE=FETRIANGLE'
        else if (SDIM.eq.2) then
         write(11,'(A18)') 'ZONETYPE=FELINESEG'
        else
         print *,"dimension bust"
         stop
        endif

       else if (ipass.eq.0) then
         ! do nothing
       else
        print *,"ipass invalid"
        stop 
       endif

       do ilev=0,finest_level
       do igrid=0,grids_per_level(ilev+1)-1
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
          ! combinetrianglessingle
         write(filename18,'(A10,A3,A5)') 'templssing',levstr,gridstr
         print *,"filename18 ",filename18

         open(unit=4,file=filename18)

         read(4,*) PartNumNodes
         read(4,*) PartNumIntElems

         if (ipass.eq.0) then
          NumNodes=NumNodes+PartNumNodes
          NumIntElems=NumIntElems+PartNumIntElems
         else if (ipass.eq.1) then
          do i=1,PartNumNodes 
           if (SDIM.eq.3) then
            read(4,*) Node(1,CurNumNodes+i),Node(2,CurNumNodes+i), &
             Node(SDIM,CurNumNodes+i)
           else if (SDIM.eq.2) then
            read(4,*) Node(1,CurNumNodes+i),Node(2,CurNumNodes+i)
           else
            print *,"dimension bust"
            stop
           endif
          enddo
          do i=1,PartNumIntElems
           if (SDIM.eq.3) then
            read(4,*) IntElem(1,CurNumIntElems+i), &
                     IntElem(2,CurNumIntElems+i), &
                     IntElem(SDIM,CurNumIntElems+i)
           else if (SDIM.eq.2) then
            read(4,*) IntElem(1,CurNumIntElems+i), &
                     IntElem(2,CurNumIntElems+i)
           else
            print *,"dimension bust"
            stop
           endif
           do dir=1,SDIM
            IntElem(dir,CurNumIntElems+i)=IntElem(dir,CurNumIntElems+i)+ &
             CurNumNodes
           enddo
          enddo

          CurNumIntElems=CurNumIntElems+PartNumIntElems
          if (CurNumIntElems.gt.NumIntElems) then
           print *,"CurNumIntElems invalid"
           stop
          endif
          CurNumNodes=CurNumNodes+PartNumNodes
          if (CurNumNodes.gt.NumNodes) then
           print *,"CurNumNodes invalid"
           stop
          endif
         else
          print *,"ipass invalid"
          stop
         endif
         close(4)

         if (ipass.eq.0) then
          ! do nothing
         else if (ipass.eq.1) then

          rmcommand14='rm templssing*'
          print *,"issuing command ",rmcommand14
          sysret=0

#ifdef PGIFORTRAN
          call system(rmcommand14)
#else
          call execute_command_line(rmcommand14,exitstat=sysret)
#endif
          if (sysret.ne.0) then
           print *,"execute_command_line has sysret=",sysret
           stop
          endif

         else
          print *,"ipass invalid"
          stop
         endif

       enddo ! igrid
       enddo ! ilev
      enddo ! ipass=0..1

      if (SDIM.eq.3) then
       do i=1,NumNodes
        write(11,*) Node(1,i),Node(2,i),Node(SDIM,i)
       enddo
       do i=1,NumIntElems
        write(11,*) IntElem(1,i),IntElem(2,i),IntElem(SDIM,i)
       enddo
      else if (SDIM.eq.2) then
       do i=1,NumNodes
        write(11,*) Node(1,i),Node(2,i)
       enddo
       do i=1,NumIntElems
        write(11,*) IntElem(1,i),IntElem(2,i)
       enddo
      else
       print *,"dimension bust"
       stop
      endif
 
      close(11)

      deallocate(Node)
      deallocate(IntElem)
 
      return
      end subroutine FORT_COMBINETRIANGLESSINGLE

      subroutine FORT_IO_COMPARE( &
       nmat, &
       nsteps, &
       do_input, &
       visual_compare, &
       time, &
       fabin,DIMS(fabin), &
       fabout,DIMS(fabout), &
       vislo,vishi, &
       vis_ncomp)
      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T nmat
      INTEGER_T nsteps
      INTEGER_T do_input
      INTEGER_T visual_compare
      REAL_T time
      INTEGER_T vis_ncomp
      INTEGER_T vislo(SDIM),vishi(SDIM) 
      INTEGER_T visual_ncell(SDIM)
      REAL_T visual_dx(SDIM)
      INTEGER_T gridlo(3),gridhi(3) 
      INTEGER_T DIMDEC(fabin)
      INTEGER_T DIMDEC(fabout)
      REAL_T fabout(DIMV(fabout),vis_ncomp) ! x,u,mag vort,LS
      REAL_T fabin(DIMV(fabin),vis_ncomp)
      REAL_T problo(SDIM),probhi(SDIM)
      INTEGER_T i,j,k
      INTEGER_T im
      INTEGER_T dir
      INTEGER_T n,n_data
      INTEGER_T strandid
      INTEGER_T gridsum(nmat)
      INTEGER_T LSgridsum(nmat)
      INTEGER_T gridsum_total
      REAL_T xtest(SDIM)
      REAL_T local_data(vis_ncomp+1) ! x,u,mag vort,LS,mag u
      REAL_T local_data_in(vis_ncomp+1)
      REAL_T local_data_out(vis_ncomp+1)
      REAL_T L1norm_local(vis_ncomp+1)
      REAL_T L2norm_local(vis_ncomp+1)
      INTEGER_T icrit(nmat*(vis_ncomp+1))
      INTEGER_T jcrit(nmat*(vis_ncomp+1))
      INTEGER_T kcrit(nmat*(vis_ncomp+1))
      REAL_T L1norm(nmat*(vis_ncomp+1))
      REAL_T L2norm(nmat*(vis_ncomp+1))
      REAL_T Linfnorm(nmat*(vis_ncomp+1))
      REAL_T LSnorm(nmat)
      REAL_T LS_in,LS_out
      INTEGER_T LS_base,state_ncomp,ebase
      character*2 im_str
      character*6 stepstr
      character*14 compfilename
      character*17 uniformfilename
      character*255 techeader_str1
      character*255 techeader_str2
      character*255 techeader_str3

      if (nmat.ne.num_materials) then
       print *,"nmat.ne.num_materials"
       stop
      endif

       ! x,u,mag vort,LS
      if (vis_ncomp.ne.2*SDIM+1+nmat) then
       print *,"vis_ncomp invalid" 
       stop
      endif

      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif
      if (nsteps.lt.0) then
       print *,"nsteps invalid"
       stop
      endif

      problo(1)=problox
      problo(2)=probloy

      probhi(1)=probhix
      probhi(2)=probhiy

      if (SDIM.eq.3) then
       problo(SDIM)=probloz
       probhi(SDIM)=probhiz
      endif

      gridlo(3)=0
      gridhi(3)=0

      do dir=1,SDIM
       if (vislo(dir).ne.0) then
        print *,"vislo(dir).ne.0"
        stop
       endif
       gridlo(dir)=vislo(dir)
       gridhi(dir)=vishi(dir)+1

       visual_ncell(dir)=vishi(dir)-vislo(dir)+1
       if (visual_ncell(dir).ge.1) then
        visual_dx(dir)=(probhi(dir)-problo(dir))/visual_ncell(dir)
        if (visual_dx(dir).gt.zero) then
         ! do nothing
        else
         print *,"visual_dx invalid"
         stop
        endif
       else
        print *,"visual_ncell invalid"
        stop
       endif
      enddo ! dir=1..sdim

      call checkbound(vislo,vishi,DIMS(fabout),0,0,411)
      call checkbound(vislo,vishi,DIMS(fabout),0,1,411)
      call checkbound(vislo,vishi,DIMS(fabout),0,SDIM-1,411)

      call checkbound(vislo,vishi,DIMS(fabin),0,0,411)
      call checkbound(vislo,vishi,DIMS(fabin),0,1,411)
      call checkbound(vislo,vishi,DIMS(fabin),0,SDIM-1,411)

      if (do_input.eq.1) then

       if (visual_compare.eq.1) then

        write(compfilename,'(A14)') 'COARSEDATA.tec'
        print *,"compfilename ",compfilename
        open(unit=11,file=compfilename)
        read(11,*) techeader_str1  ! title
        read(11,*) techeader_str2  ! variables
        read(11,*) techeader_str3  ! zone map

        do k=gridlo(3),gridhi(3)
        do j=gridlo(2),gridhi(2)
        do i=gridlo(1),gridhi(1)
         dir=1
         xtest(dir)=problo(dir)+i*visual_dx(dir)
         dir=2
         xtest(dir)=problo(dir)+j*visual_dx(dir)
         if (SDIM.eq.3) then
          dir=SDIM
          xtest(dir)=problo(dir)+k*visual_dx(dir)
         endif
         do n=1,vis_ncomp
          if ((n.ge.1).and.(n.lt.vis_ncomp)) then
           read(11,'(D25.16)',ADVANCE="NO") local_data(n)
          else if (n.eq.vis_ncomp) then
           read(11,'(D25.16)') local_data(n)
          else
           print *,"n invalid"
           stop
          endif
         enddo ! n=1..vis_ncomp

         do dir=1,SDIM
          if (abs(local_data(dir)-xtest(dir)).gt.VOFTOL*visual_dx(dir)) then
           print *,"local_data(dir) invalid"
           stop
          endif
         enddo
         do n=1,vis_ncomp
          if (abs(local_data(n)).le.1.0D+20) then
           fabin(D_DECL(i,j,k),n)=local_data(n)
          else
           print *,"local_data(n) overflow"
           stop
          endif
         enddo ! n=1..vis_ncomp
        enddo !i
        enddo !j
        enddo !k

        close(11)
       else
        print *,"visual_compare invalid"
        stop
       endif

      else if (do_input.eq.0) then

       write(stepstr,'(I6)') nsteps
       do i=1,6
        if (stepstr(i:i).eq.' ') then
         stepstr(i:i)='0'
        endif
       enddo
       write(uniformfilename,'(A7,A6,A4)') 'uniform',stepstr,'.tec'
       print *,"uniformfilename ",uniformfilename
       open(unit=11,file=uniformfilename)

       strandid=1

       write(11,*) '"uniform data"'

       if (SDIM.eq.2) then
        write(11,'(A34)',ADVANCE="NO")  &
         'VARIABLES="X","Y","U","V","MGVORT"'
       else if (SDIM.eq.3) then
        write(11,'(A42)',ADVANCE="NO") &
         'VARIABLES="X","Y","Z","U","V","W","MGVORT"'
       else
        print *,"dimension bust"
        stop
       endif

       do im=1,nmat
        write(im_str,'(I2)') im
        do i=1,2
         if (im_str(i:i).eq.' ') then
          im_str(i:i)='0'
         endif
        enddo
        write(11,'(A4)',ADVANCE="NO") ',"LS'
        write(11,'(A2)',ADVANCE="NO") im_str
        write(11,'(A1)',ADVANCE="NO") '"'
        if (im.eq.nmat) then
         write(11,*) ' '
        endif
       enddo ! im=1..nmat

       if (SDIM.eq.2) then
        write(11,*)'zone i=',visual_ncell(1)+1, &
                   ' j=',visual_ncell(SDIM)+1, &
                   ' f=point ', &
                   ' SOLUTIONTIME=',time, &
                   ' STRANDID=',strandid
       else if (SDIM.eq.3) then
        write(11,*)'zone i=',visual_ncell(1)+1, &
                   ' j=',visual_ncell(2)+1, &
                   ' k=',visual_ncell(SDIM)+1, &
                   ' f=point ', &
                   ' SOLUTIONTIME=',time, &
                   ' STRANDID=',strandid
       else
        print *,"dimension bust"
        stop
       endif

       do k=gridlo(3),gridhi(3)
       do j=gridlo(2),gridhi(2)
       do i=gridlo(1),gridhi(1)

        dir=1
        xtest(dir)=problo(dir)+i*visual_dx(dir)
        dir=2
        xtest(dir)=problo(dir)+j*visual_dx(dir)
        if (SDIM.eq.3) then
         dir=SDIM
         xtest(dir)=problo(dir)+k*visual_dx(dir)
        endif
        do n=1,vis_ncomp
         local_data(n)=fabout(D_DECL(i,j,k),n)
        enddo
        do dir=1,SDIM
         if (abs(local_data(dir)-xtest(dir)).gt.VOFTOL*visual_dx(dir)) then
          print *,"local_data(dir) invalid"
          stop
         endif
        enddo
        do n=1,vis_ncomp
         if (abs(local_data(n)).le.1.0D+20) then
          if ((n.ge.1).and.(n.lt.vis_ncomp)) then
           write(11,'(D25.16)',ADVANCE="NO") local_data(n)
          else if (n.eq.vis_ncomp) then
           write(11,'(D25.16)') local_data(n)
          else
           print *,"n invalid"
           stop
          endif
         else
          print *,"local_data(n) overflow"
          stop
         endif
        enddo ! n=1..vis_ncomp
 
       enddo ! i
       enddo ! j
       enddo ! k

       close(11)

       if (visual_compare.eq.1) then

        LS_base=vis_ncomp-nmat
        state_ncomp=vis_ncomp-nmat+1

        do n=1,nmat*(vis_ncomp+1)
         L1norm(n)=zero
         L2norm(n)=zero
         Linfnorm(n)=zero
         icrit(n)=-1
         jcrit(n)=-1
         kcrit(n)=-1
        enddo ! n=1..nmat*(vis_ncomp+1)

        do im=1,nmat
         gridsum(im)=0
         LSgridsum(im)=0
         LSnorm(im)=zero
        enddo

        gridsum_total=0

        do k=gridlo(3),gridhi(3)
        do j=gridlo(2),gridhi(2)
        do i=gridlo(1),gridhi(1)
           ! x,u,mag vort,LS
         do n=1,vis_ncomp
          local_data_in(n)=fabin(D_DECL(i,j,k),n)
          local_data_out(n)=fabout(D_DECL(i,j,k),n)
          if ((abs(local_data_in(n)).le.1.0D+20).and. &
              (abs(local_data_out(n)).le.1.0D+20)) then

           local_data(n)=local_data_in(n)-local_data_out(n)
           
           if (abs(local_data(n)).le.1.0D+20) then
            L1norm_local(n)=abs(local_data(n))
            L2norm_local(n)=abs(local_data(n))**2
           else
            print *,"abs(local_data(n)) overflow"
            stop
           endif
          else
           print *,"local_data_in or local_data_out corrupt"
           stop
          endif 
         enddo ! n=1..vis_ncomp

          ! x,u,mag vort,LS,mag u
         n=vis_ncomp+1
         local_data(n)=zero
         do dir=1,SDIM
          local_data(n)=local_data(n)+ &
           abs(local_data(dir+SDIM))**2
         enddo
         local_data(n)=sqrt(local_data(n))

         if (abs(local_data(n)).le.1.0D+20) then
          L1norm_local(n)=abs(local_data(n))
          L2norm_local(n)=abs(local_data(n))**2
         else
          print *,"abs(local_data(n)) overflow"
          stop
         endif

         if ((LS_base.eq.vis_ncomp-nmat).and. &
             (state_ncomp.eq.vis_ncomp-nmat+1)) then
          do im=1,nmat
           LS_in=local_data_in(LS_base+im)
           LS_out=local_data_out(LS_base+im)
           if ((LS_in.gt.zero).and.(LS_out.gt.zero)) then
            ebase=(im-1)*state_ncomp
            do n=1,state_ncomp
             if (n.lt.state_ncomp) then
              n_data=n
             else
              n_data=vis_ncomp+1
             endif
             L1norm(ebase+n)=L1norm(ebase+n)+L1norm_local(n_data)
             L2norm(ebase+n)=L2norm(ebase+n)+L2norm_local(n_data)
             if (abs(local_data(n_data)).gt.Linfnorm(ebase+n)) then
              Linfnorm(ebase+n)=abs(local_data(n_data))
              icrit(ebase+n)=i
              jcrit(ebase+n)=j
              kcrit(ebase+n)=k
             endif
            enddo ! n=1..state_ncomp
            gridsum(im)=gridsum(im)+1
            gridsum_total=gridsum_total+1
           else if ((LS_in.le.zero).or.(LS_out.le.zero)) then
            ! do nothing
           else
            print *,"LS_in or LS_out invalid"
            stop
           endif
          enddo ! im=1..nmat
          do im=1,nmat
           LS_in=local_data_in(LS_base+im)
           LS_out=local_data_out(LS_base+im)
           if ((abs(LS_in).le.visual_dx(1)).or. &
               (abs(LS_out).le.visual_dx(1))) then
            LSnorm(im)=LSnorm(im)+abs(local_data(LS_base+im))
            LSgridsum(im)=LSgridsum(im)+1
           else if ((abs(LS_in).gt.visual_dx(1)).and. &
                    (abs(LS_out).gt.visual_dx(1))) then
            ! do nothing
           else
            print *,"LS_in or LS_out invalid"
            stop
           endif
          enddo ! im=1..nmat
         else
          print *,"LS_base or state_ncomp became corrupt"
          stop
         endif
          
        enddo ! i
        enddo ! j
        enddo ! k

        if (gridsum_total.ge.1) then
         if ((LS_base.eq.vis_ncomp-nmat).and. &
             (state_ncomp.eq.vis_ncomp-nmat+1)) then
          do im=1,nmat
           if (gridsum(im).gt.0) then
            ebase=(im-1)*state_ncomp
            do n=1,state_ncomp
             L1norm(ebase+n)=L1norm(ebase+n)/gridsum(im)
             L2norm(ebase+n)=sqrt(L2norm(ebase+n)/gridsum(im))
             print *,"im,n,cells,L1,L2,Linf,ic,jc,kc ", &
              im,n,gridsum(im), &
              L1norm(ebase+n),L2norm(ebase+n),Linfnorm(ebase+n), &
              icrit(ebase+n),jcrit(ebase+n),kcrit(ebase+n)
            enddo ! n=1..state_ncomp
           else if (gridsum(im).eq.0) then
            ! do nothing
           else
            print *,"gridsum invalid"
            stop
           endif
           if (LSgridsum(im).gt.0) then
            LSnorm(im)=LSnorm(im)/LSgridsum(im)
            print *,"im,cells,LSnorm ",im,LSgridsum(im),LSnorm(im)
           else if (LSgridsum(im).eq.0) then
            ! do nothing
           else
            print *,"LSgridsum invalid"
            stop
           endif
          enddo ! im=1..nmat
         else
          print *,"LS_base or state_ncomp became corrupt"
          stop
         endif
        else
         print *,"gridsum_total invalid"
         stop
        endif
            
       else if (visual_compare.eq.0) then
        ! do nothing
       else
        print *,"visual_compare invalid"
        stop
       endif

      else
       print *,"do_input invalid"
       stop
      endif
     
      return
      end subroutine FORT_IO_COMPARE

      subroutine FORT_CELLGRID( &
       tid, &
       bfact, &
       fabout,DIMS(fabout), &
       vislo,vishi, &
       vis_ncomp, & ! x,u,mag vort,LS
       maskSEM,DIMS(maskSEM), &
       vel,DIMS(vel), &
       velsol,DIMS(velsol), &
       vof,DIMS(vof), &
       pres,DIMS(pres), &
       div,DIMS(div), &
       divdat,DIMS(divdat), &
       den,DIMS(den), &
       elastic,DIMS(elastic), &
       lsdist,DIMS(lsdist), &
       visc,DIMS(visc), &
       trace,DIMS(trace), &
       problo, &
       probhi, &
       dx, &
       lo,hi, &
       level, &
       finest_level, &
       gridno, &
       visual_tessellate_vfrac, &
       visual_option, &
       rz_flag, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       elastic_ncomp, &
       slice_data, &
       nslice, &
       nstate_slice,slice_dir, &
       xslice, &
       dxfinest, &
       do_plot,do_slice)

      use global_utility_module
      use probf90_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T tid
      INTEGER_T bfact
      INTEGER_T do_plot,do_slice

        ! nstate_slice=x,y,z,xvel,yvel,zvel,PMG,PEOS,den,Temp,KE
        ! (value of material with LS>0)
        ! nslice=domhi-domlo+3
      INTEGER_T nslice,nstate_slice,slice_dir
      REAL_T slice_data(nslice*nstate_slice)
      REAL_T xslice(SDIM)
 
      INTEGER_T rz_flag
      INTEGER_T nmat
      INTEGER_T nparts
      INTEGER_T nparts_def
      INTEGER_T im_solid_map(nparts_def) 
      INTEGER_T elastic_ncomp
      INTEGER_T visual_tessellate_vfrac
      INTEGER_T visual_option
      INTEGER_T vis_ncomp
      INTEGER_T vislo(SDIM), vishi(SDIM)
      INTEGER_T lo(SDIM), hi(SDIM)
      INTEGER_T DIMDEC(fabout)
      INTEGER_T DIMDEC(maskSEM)
      INTEGER_T DIMDEC(vel)
      INTEGER_T DIMDEC(velsol)
      INTEGER_T DIMDEC(vof)
      INTEGER_T DIMDEC(pres)
      INTEGER_T DIMDEC(div)
      INTEGER_T DIMDEC(divdat)
      INTEGER_T DIMDEC(den)
      INTEGER_T DIMDEC(elastic)
      INTEGER_T DIMDEC(lsdist)
      INTEGER_T DIMDEC(visc)
      INTEGER_T DIMDEC(trace)
      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T gridno
      REAL_T fabout(DIMV(fabout),vis_ncomp) ! x,u,mag vort,LS
      REAL_T maskSEM(DIMV(maskSEM))
      REAL_T vel(DIMV(vel), &
        num_materials_vel*(SDIM+1))
      REAL_T velsol(DIMV(velsol), &
        nparts_def*SDIM)
      REAL_T vof(DIMV(vof),nmat*ngeom_recon)
      REAL_T pres(DIMV(pres),num_materials_vel)
      REAL_T div(DIMV(div),num_materials_vel)
      REAL_T divdat(DIMV(divdat),num_materials_vel)
      REAL_T den(DIMV(den),num_state_material*nmat)
      REAL_T elastic(DIMV(elastic),elastic_ncomp)
      REAL_T visc(DIMV(visc),nmat)
      REAL_T trace(DIMV(trace),5*nmat)
      REAL_T lsdist(DIMV(lsdist),(SDIM+1)*nmat)
      REAL_T xposnd(SDIM)
      REAL_T xposndT(SDIM)
      REAL_T machnd(num_materials_vel)
      REAL_T machcell(num_materials_vel)
      REAL_T velnd(num_materials_vel*(SDIM+1))
      REAL_T velsolidnd(nparts_def*SDIM)
      REAL_T velmat(SDIM)
      REAL_T velmatT(SDIM)
      REAL_T velcell(num_materials_vel*(SDIM+1))
      REAL_T velsolidcell(nparts_def*SDIM)
      REAL_T vofnd(nmat)
      REAL_T vofcell(nmat)
      REAL_T presnd(num_materials_vel)
      REAL_T divnd(num_materials_vel)
      REAL_T divdatnd(num_materials_vel)
      REAL_T dennd(num_state_material*nmat)
      REAL_T elasticnd(elastic_ncomp)
      REAL_T dencell(num_state_material*nmat)
      REAL_T elasticcell(elastic_ncomp)
      REAL_T lsdistnd((SDIM+1)*nmat)
      REAL_T viscnd(nmat)
      REAL_T tracend(5*nmat)
      REAL_T writend(nmat*200)
      INTEGER_T scomp,iw
      INTEGER_T istate,idissolution
      INTEGER_T im
      INTEGER_T im_crit

      REAL_T problo(SDIM)
      REAL_T probhi(SDIM)
      REAL_T dx(SDIM)
      REAL_T dxfinest(SDIM)
      REAL_T dxelem

      character*3 levstr
      character*5 gridstr
      character*18 filename18

      INTEGER_T i,j,k
      INTEGER_T ii,jj
      INTEGER_T dir,dir2
      INTEGER_T n
      INTEGER_T i1,j1,k1
      INTEGER_T ic,jc,kc
      INTEGER_T iBL,jBL,kBL
      INTEGER_T iSEM,jSEM,kSEM
      INTEGER_T ilocal,jlocal,klocal
      REAL_T sumweight,localwt,int_xlo,int_xhi
      REAL_T sumweightLS,localwtLS,wt_denom
      INTEGER_T vofcomp
      INTEGER_T k1lo,k1hi
      INTEGER_T igridlo(3),igridhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      INTEGER_T stenlo_index
      INTEGER_T icrit(SDIM)
      INTEGER_T icritlo(3),icrithi(3)
      INTEGER_T icell_lo(3),icell_hi(3)
      INTEGER_T icell(3)
      INTEGER_T iproblo(3)
      REAL_T xfablo(SDIM),xfabhi(SDIM)
      REAL_T x1D
      REAL_T xc(SDIM),theta(3)
      REAL_T xstenlo(SDIM),xstenhi(SDIM)
      INTEGER_T inbox
      REAL_T psten(0:1,0:1)
      REAL_T psten1D(0:1)
      INTEGER_T DIMDEC(plt)
      REAL_T, dimension(D_DECL(:,:,:),:), allocatable :: plotfab
      REAL_T, dimension(D_DECL(:,:,:),:), allocatable :: reconfab
      INTEGER_T debug_slice
      REAL_T denslice,tempslice,eslice,KEslice
      INTEGER_T nhalf,bfact_finest
      REAL_T xstenND(-3:3,SDIM)
      REAL_T xsten_corner(-3:3,SDIM)
      REAL_T xsten(-3:3,SDIM)
      REAL_T xsten1D(-3:3)
      REAL_T xsten1DL(-3:3)
      REAL_T xsten1D_finest(-3:3)
      REAL_T xsten_fablo(-3:3,SDIM)
      REAL_T xsten_fabhi(-3:3,SDIM)
      REAL_T dxleft,dxright,dxmin
      REAL_T mofdata(nmat*ngeom_recon)
      INTEGER_T nmax
      INTEGER_T partid

      INTEGER_T visual_ncell(SDIM)
      REAL_T visual_dx(SDIM)
      REAL_T xcrit(SDIM)
      REAL_T localfab(vis_ncomp)
      REAL_T SEM_value(vis_ncomp-SDIM) ! u,mag vort,LS
      INTEGER_T ncomp_SEM  ! vis_ncomp-SDIM
      INTEGER_T VORT_comp_SEM ! SDIM+1
      INTEGER_T gridtype
      INTEGER_T SEMhi(SDIM)
      REAL_T, dimension(D_DECL(:,:,:),:),allocatable :: SEMloc
      REAL_T INTERP_TOL
      INTEGER_T local_maskSEM
      REAL_T local_data
      INTEGER_T caller_id

      caller_id=3

      nhalf=3
      nmax=POLYGON_LIST_MAX ! in: CELLGRID
      bfact_finest=2
      INTERP_TOL=1.0E-4

      debug_slice=0

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 39"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_CELLGRID"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_CELLGRID"
       stop
      endif
      if ((num_materials_viscoelastic.ge.1).and. &
          (num_materials_viscoelastic.le.nmat)) then
       if (elastic_ncomp.ne. &
           num_materials_viscoelastic*FORT_NUM_TENSOR_TYPE) then
        print *,"elastic_ncomp invalid 1"
        stop
       endif
      else if (num_materials_viscoelastic.eq.0) then
       if (elastic_ncomp.ne.num_state_material*nmat) then
        print *,"elastic_ncomp invalid 2"
        stop
       endif
      else
       print *,"num_materials_viscoelastic invalid"
       stop
      endif
 
        ! nstate_slice=x,y,z,xvel,yvel,zvel,PMG,PEOS,DIV,den,Temp,KE
        ! (value of material with LS>0)
      if (nstate_slice.ne.2*SDIM+6) then
       print *,"nstate_slice invalid in cellgrid"
       print *,"nstate_slice= ",nstate_slice
       stop
      endif 
      if ((slice_dir.lt.0).or.(slice_dir.ge.SDIM)) then
       print *,"slice_dir invalid"
       stop
      endif
       ! nslice=domhi-domlo+3
      dir=slice_dir+1
      if (nslice.lt.hi(dir)+3) then
       print *,"nslice invalid"
       stop
      endif

      if (debug_slice.eq.1) then
       print *,"debug_slice: slice_dir= ",slice_dir
       print *,"debug_slice: nslice= ",nslice
       print *,"debug_slice: nstate_slice= ",nstate_slice
       do dir=1,SDIM
        print *,"debug_slice: dir,xslice ",dir,xslice(dir)
        print *,"debug_slice: dir,problo ",dir,problo(dir)
        print *,"debug_slice: dir,dx ",dir,dx(dir)
        print *,"debug_slice: dir,dxfinest ",dir,dxfinest(dir)
       enddo
      endif

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      igridlo(3)=0
      igridhi(3)=0
      iproblo(3)=0
      do dir=1,SDIM
       igridlo(dir)=lo(dir)-1
       igridhi(dir)=hi(dir)+1
       iproblo(dir)=0
      enddo
       ! input: igridlo,igridhi
       ! output: DIMS(plt)
      call box_to_dim( &
        DIMS(plt), &
        igridlo,igridhi)
 
      allocate(plotfab(DIMV(plt),nstate_slice))
      allocate(reconfab(DIMV(plt),nmat*ngeom_recon))
 
      call checkbound(vislo,vishi,DIMS(fabout),0,0,411)
      call checkbound(vislo,vishi,DIMS(fabout),0,1,411)
      call checkbound(vislo,vishi,DIMS(fabout),0,SDIM-1,411)
       ! x,u,mag vort,LS
      if (vis_ncomp.ne.2*SDIM+1+nmat) then
       print *,"vis_ncomp invalid" 
       stop
      endif

      call checkbound(lo,hi,DIMS(maskSEM),0,-1,1264)

      call checkbound(lo,hi,DIMS(plt),1,-1,411)
      call checkbound(lo,hi,DIMS(pres),1,-1,411)
      call checkbound(lo,hi,DIMS(div),1,-1,411)
      call checkbound(lo,hi,DIMS(divdat),1,-1,411)
      call checkbound(lo,hi,DIMS(den),1,-1,411)
      call checkbound(lo,hi,DIMS(elastic),1,-1,411)
      call checkbound(lo,hi,DIMS(lsdist),1,-1,411)
      call checkbound(lo,hi,DIMS(visc),1,-1,411)
      call checkbound(lo,hi,DIMS(trace),1,-1,411)
      call checkbound(lo,hi,DIMS(vel),1,-1,411)
      call checkbound(lo,hi,DIMS(velsol),1,-1,411)
      call checkbound(lo,hi,DIMS(vof),1,-1,411)

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid cell grid"
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

      do k=igridlo(3),igridhi(3) 
      do j=igridlo(2),igridhi(2) 
      do i=igridlo(1),igridhi(1) 
         ! iproblo=0
       call gridsten(xsten,problo,i,j,k,iproblo,bfact,dx,nhalf)
       do dir=1,nmat*ngeom_recon
        mofdata(dir)=vof(D_DECL(i,j,k),dir)
       enddo

       if (visual_tessellate_vfrac.eq.1) then
         ! before (mofdata): fluids tessellate
         ! after  (mofdata): fluids and solids tessellate
        call multi_get_volume_tessellate( &
         bfact,dx,xsten,nhalf, &
         mofdata, &
         geom_xtetlist(1,1,1,tid+1), &
         nmax, &
         nmax, &
         nmat,SDIM,8)
       else if (visual_tessellate_vfrac.eq.0) then
        ! do nothing
       else
        print *,"visual_tessellate_vfrac invalid"
        stop
       endif

       do dir=1,nmat*ngeom_recon
        reconfab(D_DECL(i,j,k),dir)=mofdata(dir)
       enddo
      enddo  ! i
      enddo  ! j
      enddo  ! k

      if (do_plot.eq.1) then

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

       if (SDIM.eq.3) then
        k1lo=0
        k1hi=1
       else if (SDIM.eq.2) then
        k1lo=0
        k1hi=0
       else
        print *,"dimension bust"
        stop
       endif

       do dir=1,SDIM
        igridlo(dir)=lo(dir)
        igridhi(dir)=hi(dir)+1
       enddo

        ! the order k,j,i is IMPORTANT.
       do k=igridlo(3),igridhi(3) 
       do j=igridlo(2),igridhi(2) 
       do i=igridlo(1),igridhi(1) 

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
           ! the target interpolating box should be centered
           ! about xstenND(0,dir)
          dxmin=min(dxleft,dxright)
          xstenND(-1,dir)=xstenND(0,dir)-dxmin
          xstenND(1,dir)=xstenND(0,dir)+dxmin
         else
          print *,"bfact invalid152"
          stop
         endif
         xposnd(dir)=xstenND(0,dir)
        enddo ! dir

        do dir=1,SDIM
         xposndT(dir)=xposnd(dir)
        enddo
        if (visual_RT_transform.eq.1) then
         call RT_transform(xposnd,xposndT)
        endif

        sumweight=zero
        sumweightLS=zero

        do dir=1,num_materials_vel*(SDIM+1)
         velnd(dir)=zero
        enddo
        do dir=1,nparts_def*SDIM
         velsolidnd(dir)=zero
        enddo
        do dir=1,num_state_material*nmat
         dennd(dir)=zero
        enddo
        do dir=1,elastic_ncomp
         elasticnd(dir)=zero
        enddo
        do dir=1,num_materials_vel
         presnd(dir)=zero
         divnd(dir)=zero
         divdatnd(dir)=zero
         machnd(dir)=zero
        enddo
        do dir=1,nmat
         vofnd(dir)=zero
         viscnd(dir)=zero
        enddo
        do dir=1,nmat*(SDIM+1)
         lsdistnd(dir)=zero
        enddo
        do dir=1,5*nmat
         tracend(dir)=zero
        enddo

         ! two types of low order interpolation:
         ! (a) the weights are the area of the intersection of the node control
         !     volume with a given cell control volume.
         ! or 
         ! (b) linear interpolation (see as follows for the LS function)
         !
         ! linear interpolation for the LS functions:
         ! y=y0 (x-x1)/(x0-x1) + y1 (x-x0)/(x1-x0)
         ! yI=y0 (xI-x1)/(x0-x1)+y1 (xI-x0)/(x1-x0)=
         !   abs((xI-x1)(xI-x0)/(x1-x0))(y0/abs(xI-x0)+y1/abs(xI-x1))

        icell(1)=i
        icell(2)=j
        icell(3)=k

        do i1=0,1
        do j1=0,1
        do k1=k1lo,k1hi

         dir=1
         if ((i-i1.ge.lo(dir)).and. &
             (i-i1.le.hi(dir))) then
          icell(dir)=i-i1
         endif
         dir=2
         if ((j-j1.ge.lo(dir)).and. &
             (j-j1.le.hi(dir))) then
          icell(dir)=j-j1
         endif
         if (SDIM.eq.3) then
          dir=SDIM
          if ((k-k1.ge.lo(dir)).and. &
              (k-k1.le.hi(dir))) then
           icell(dir)=k-k1
          endif
         endif

          ! iproblo=0
         call gridsten(xsten,problo,i-i1,j-j1,k-k1,iproblo,bfact,dx,nhalf)
         
         localwt=one  ! area of intersection of node CV with given cell CV
         localwtLS=one

         do dir=1,SDIM
          dxleft=xsten(0,dir)-xsten(-1,dir)
          dxright=xsten(1,dir)-xsten(0,dir)
          
          if (bfact.eq.1) then
           if (abs(dxleft-dxright).gt.VOFTOL*dx(dir)) then
            print *,"xsten invalid"
            stop
           endif
          else if (bfact.gt.1) then
           if ((dxleft.le.zero).or.(dxright.le.zero)) then
            print *,"xsten invalid"
            stop
           endif
          else
           print *,"bfact invalid153"
           stop
          endif

           ! used for LS interpolation (bilinear interp)
          wt_denom=abs(xsten(0,dir)-xstenND(0,dir))
          if (wt_denom.le.VOFTOL*dx(dir)) then
           print *,"wt_denom invalid"
           stop
          endif

          int_xlo=max(xstenND(-1,dir),xsten(-1,dir)) 
          int_xhi=min(xstenND(1,dir),xsten(1,dir)) 
          if (int_xhi.gt.int_xlo) then
           localwt=localwt*(int_xhi-int_xlo)
          else
           localwt=zero
          endif

          localwtLS=localwtLS*(one/wt_denom)

         enddo ! dir=1..sdim

         if (localwt.le.zero) then
          print *,"localwt invalid"
          stop
         endif
         if (localwtLS.le.zero) then
          print *,"localwtLS invalid"
          stop
         endif

         do dir=1,num_materials_vel*(SDIM+1)
          velcell(dir)=vel(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         do dir=1,num_materials_vel*(SDIM+1)
          velnd(dir)=velnd(dir)+localwt*velcell(dir)
         enddo
         do dir=1,nparts_def*SDIM
          velsolidcell(dir)=velsol(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         do dir=1,nparts_def*SDIM
          velsolidnd(dir)=velsolidnd(dir)+localwt*velsolidcell(dir)
         enddo
         do dir=1,num_state_material*nmat
          dencell(dir)=den(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo

          ! e.g. transform temperature for rotating convection instability
          ! problem.
         call derive_plot_data(xsten,nhalf,dencell,nmat)

         do dir=1,num_state_material*nmat
          dennd(dir)=dennd(dir)+localwt*dencell(dir)
         enddo

         do dir=1,elastic_ncomp
          elasticcell(dir)=elastic(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo

         do dir=1,elastic_ncomp
          elasticnd(dir)=elasticnd(dir)+localwt*elasticcell(dir)
         enddo

         do dir=1,num_materials_vel
          presnd(dir)=presnd(dir)+ &
           localwt*pres(D_DECL(i-i1,j-j1,k-k1),dir)
          divnd(dir)=divnd(dir)+ &
           localwt*div(D_DECL(i-i1,j-j1,k-k1),dir)
          divdatnd(dir)=divdatnd(dir)+ &
           localwt*divdat(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         do dir=1,nmat
          vofcomp=(dir-1)*ngeom_recon+1
          vofcell(dir)=reconfab(D_DECL(i-i1,j-j1,k-k1),vofcomp)
         enddo
         do dir=1,nmat
          vofnd(dir)=vofnd(dir)+localwt*vofcell(dir)
          viscnd(dir)=viscnd(dir)+ &
            localwt*visc(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         do dir=1,nmat*(SDIM+1)
          lsdistnd(dir)=lsdistnd(dir)+ &
            localwtLS*lsdist(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo


         ! 1. \dot{gamma}
         ! 2. Tr(A) if viscoelastic
         !    \dot{gamma} o.t.
         ! 3. Tr(A) (liquid viscosity - etaS)/etaP  if FENE-CR+Carreau
         !    Tr(A) if FENE-CR
         !    \dot{gamma} o.t.
         ! 4. (3) * f(A)  if viscoelastic
         !    \dot{gamma} o.t.
         ! 5. vorticity magnitude.

         do dir=1,5*nmat
          tracend(dir)=tracend(dir)+ &
            localwt*trace(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         call get_mach_number(visual_tessellate_vfrac, &
           velcell,dencell,vofcell,machcell,nmat)

         do dir=1,num_materials_vel
          machnd(dir)=machnd(dir)+localwt*machcell(dir)
         enddo

         sumweight=sumweight+localwt
         sumweightLS=sumweightLS+localwtLS
        enddo ! k1
        enddo ! j1
        enddo ! i1

        if (sumweight.le.zero) then
         print *,"sumweight invalid"
         stop
        endif
        if (sumweightLS.le.zero) then
         print *,"sumweightLS invalid"
         stop
        endif

        do dir=1,num_materials_vel*(SDIM+1)
         velnd(dir)=velnd(dir)/sumweight
        enddo
        do dir=1,nparts_def*SDIM
         velsolidnd(dir)=velsolidnd(dir)/sumweight
        enddo
        do dir=1,5*nmat
         tracend(dir)=tracend(dir)/sumweight
        enddo

        if (bfact.eq.1) then
         ! do nothing
        else if (bfact.gt.1) then
         local_maskSEM= &
          NINT(maskSEM(D_DECL(icell(1),icell(2),icell(3))))
         if ((local_maskSEM.ge.1).and. &
             (local_maskSEM.le.nmat)) then

           ! iproblo=0
          call gridsten(xsten,problo,icell(1),icell(2),icell(3), &
            iproblo,bfact,dx,nhalf)
          do dir=1,SDIM
           xcrit(dir)=xstenND(0,dir)
          enddo
          do dir=1,SDIM
           if (xcrit(dir).le.xsten(-1,dir)) then
            xcrit(dir)=xsten(-1,dir)+1.0E-14*dx(dir)
           endif
           if (xcrit(dir).ge.xsten(1,dir)) then
            xcrit(dir)=xsten(1,dir)-1.0E-14*dx(dir)
           endif
          enddo

          ! find stencil for surrounding spectral element.
          stenlo(3)=0
          stenhi(3)=0
          do dir2=1,SDIM
           dxelem=dx(dir2)*bfact
           stenlo_index=NINT( (xsten(0,dir2)-problo(dir2))/dxelem-half )
           stenlo(dir2)=bfact*stenlo_index
           stenhi(dir2)=stenlo(dir2)+bfact-1
           if ((icell(dir2).lt.stenlo(dir2)).or. &
               (icell(dir2).gt.stenhi(dir2))) then
            print *,"SEM stencil should contain cell"
            print *,"dir2,icell(dir2),stenlo(dir2),stenhi(dir2) ", &
             dir2,icell(dir2),stenlo(dir2),stenhi(dir2)
            print *,"xcrit(dir2) ",xcrit(dir2)
            print *,"xsten(-1,dir2) ",xsten(-1,dir2)
            print *,"xsten(0,dir2) ",xsten(0,dir2)
            print *,"xsten(1,dir2) ",xsten(1,dir2)
            stop
           endif
          enddo ! dir2=1..sdim

          iSEM=stenlo(1) 
          jSEM=stenlo(2) 
          kSEM=stenlo(SDIM) 

          ! lower left hand corner of element stencil.
          ! iproblo=0
          call gridstenND(xsten_corner,problo,iSEM,jSEM,kSEM,iproblo,bfact, &
            dx,nhalf)

          do dir2=1,SDIM

           if (abs(xsten_corner(0,dir2)).le.1.0D+20) then
            xcrit(dir2)=xcrit(dir2)-xsten_corner(0,dir2)

            if ((xcrit(dir2).ge.-INTERP_TOL*dx(dir2)).and. &
                (xcrit(dir2).lt.(INTERP_TOL+bfact)*dx(dir2))) then
             ! do nothing
            else
             print *,"xcrit out of bounds CELLGRID"
             print *,"dir2,xcrit,xsten_corner ",dir2,xcrit(dir2), &
              xsten_corner(0,dir2)
             stop
            endif
           else
            print *,"xsten_corner invalid" 
            stop
           endif

          enddo ! dir2=1..sdim

          ncomp_SEM=vis_ncomp-SDIM
          VORT_comp_SEM=SDIM+1
          do n=1,ncomp_SEM
           SEM_value(n)=zero
          enddo
          gridtype=0  ! ggg  (Gauss in all directions)
          do dir2=1,SDIM
           SEMhi(dir2)=bfact-1
          enddo

          allocate(SEMloc(D_DECL(0:SEMhi(1),0:SEMhi(2),0:SEMhi(SDIM)), &
                          ncomp_SEM))

          do iSEM=stenlo(1),stenhi(1)
          do jSEM=stenlo(2),stenhi(2)
          do kSEM=stenlo(3),stenhi(3)
           ilocal=iSEM-stenlo(1)
           jlocal=jSEM-stenlo(2)
           klocal=kSEM-stenlo(3)
           do dir=1,SDIM
            local_data=vel(D_DECL(iSEM,jSEM,kSEM),dir)
            if (abs(local_data).lt.1.0D+20) then
             SEMloc(D_DECL(ilocal,jlocal,klocal),dir)=local_data
            else
             print *,"abs(local_data) overflow1"
             stop
            endif
           enddo ! dir=1..sdim
           local_data=trace(D_DECL(iSEM,jSEM,kSEM),(local_maskSEM-1)*5+5)
           if (VORT_comp_SEM.eq.SDIM+1) then
            if (abs(local_data).lt.1.0D+20) then
             SEMloc(D_DECL(ilocal,jlocal,klocal),VORT_comp_SEM)=local_data
            else
             print *,"abs(local_data) overflow2"
             stop
            endif
           else
            print *,"VORT_comp_SEM invalid"
            stop
           endif
            ! we interpolate the LS using SEM, but we discard the interpolated
            ! value.
           do im=1,nmat
            local_data=lsdist(D_DECL(iSEM,jSEM,kSEM),im)
            if (abs(local_data).lt.1.0D+20) then
             SEMloc(D_DECL(ilocal,jlocal,klocal),VORT_comp_SEM+im)=local_data
            else
             print *,"abs(local_data) overflow25"
             stop
            endif
           enddo ! im=1..nmat
            
          enddo ! kSEM
          enddo ! jSEM
          enddo ! iSEM

          call SEM_INTERP_ELEMENT( &
           ncomp_SEM,bfact,gridtype, &
           SEMhi,dx,xcrit,SEMloc,SEM_value,caller_id)

          deallocate(SEMloc)

          do dir=1,SDIM
           if (abs(SEM_value(dir)).lt.1.0D+20) then
            velnd(dir)=SEM_value(dir)
           else 
            print *,"abs(SEM_value(dir)) overflow"
            print *,"i,j,k,level,finest_level,dir ", &
             i,j,k,level,finest_level,dir 
            print *,"stenlo ",stenlo(1),stenlo(2),stenlo(3)
            print *,"stenhi ",stenhi(1),stenhi(2),stenhi(3)
            print *,"lo ",lo(1),lo(2),lo(SDIM)
            print *,"hi ",hi(1),hi(2),hi(SDIM)
            print *,"SEM_value(dir) ",SEM_value(dir)
            stop
           endif
          enddo ! dir=1..sdim
          if (vis_ncomp-SDIM.ne.SDIM+1+nmat) then
           print *,"incorrect vis_ncomp"
           stop
          endif
          if (ncomp_SEM.ne.SDIM+1+nmat) then
           print *,"incorrect ncomp_SEM"
           stop
          endif
          if (abs(SEM_value(VORT_comp_SEM)).lt.1.0D+20) then
           tracend((local_maskSEM-1)*5+5)=SEM_value(VORT_comp_SEM)
          else 
           print *,"abs(SEM_value(sdim+1)) overflow"
           print *,"i,j,k,level,finest_level ", &
            i,j,k,level,finest_level
           print *,"stenlo ",stenlo(1),stenlo(2),stenlo(3)
           print *,"stenhi ",stenhi(1),stenhi(2),stenhi(3)
           print *,"lo ",lo(1),lo(2),lo(SDIM)
           print *,"hi ",hi(1),hi(2),hi(SDIM)
           print *,"SEM_value(VORT_comp_SEM) ",SEM_value(VORT_comp_SEM)
           stop
          endif

         else if (local_maskSEM.eq.0) then
          ! do nothing
         else
          print *,"local_maskSEM invalid"
          stop
         endif

        else
         print *,"bfact invalid154"
         stop
        endif

        do dir=1,num_state_material*nmat
         dennd(dir)=dennd(dir)/sumweight
        enddo

        do dir=1,elastic_ncomp
         elasticnd(dir)=elasticnd(dir)/sumweight
        enddo

        do dir=1,num_materials_vel
         presnd(dir)=presnd(dir)/sumweight
         divnd(dir)=divnd(dir)/sumweight
         divdatnd(dir)=divdatnd(dir)/sumweight
         machnd(dir)=machnd(dir)/sumweight
        enddo

        do dir=1,nmat
         vofnd(dir)=vofnd(dir)/sumweight
         viscnd(dir)=viscnd(dir)/sumweight
        enddo
        do dir=1,nmat*(SDIM+1)
         lsdistnd(dir)=lsdistnd(dir)/sumweightLS
        enddo

        if (probtype.eq.5700) then
         if (nmat.eq.3) then
          if (lsdistnd(3).ge.zero) then
           if (lsdistnd(1).ge.-lsdistnd(3)) then
            lsdistnd(1)=-lsdistnd(3)
           endif
           if (lsdistnd(2).ge.-lsdistnd(3)) then
            lsdistnd(2)=-lsdistnd(3)
           endif
          endif
         else
          print *,"nmat invalid for 5700 option"
          stop
         endif
        endif  ! probtype=5700

        scomp=0
        do iw=1,SDIM
         writend(scomp+iw)=xposndT(iw)
        enddo

        scomp=scomp+SDIM

        do dir=1,SDIM
         velmat(dir)=velnd(dir)
         velmatT(dir)=velmat(dir)
        enddo
        if (visual_RT_transform.eq.1) then
         call RT_transformVEL(xposnd,velmat,velmatT)
        endif
        do iw=1,SDIM
         writend(scomp+iw)=velmatT(iw)
        enddo
        scomp=scomp+SDIM
  
        do partid=0,nparts_def-1
        
         do dir=1,SDIM
          velmat(dir)=velsolidnd(partid*SDIM+dir)
          velmatT(dir)=velmat(dir)
         enddo
         if (visual_RT_transform.eq.1) then
          call RT_transformVEL(xposnd,velmat,velmatT)
         endif
         do iw=1,SDIM
          writend(scomp+iw)=velmatT(iw)
         enddo
         scomp=scomp+SDIM

        enddo ! partid=0..nparts_def-1

          ! this is pressure from the projection.
        do iw=1,num_materials_vel
         writend(scomp+iw)=velnd(num_materials_vel*SDIM+iw)
        enddo
        scomp=scomp+num_materials_vel

          ! this is EOS pressure
        do iw=1,num_materials_vel
         writend(scomp+iw)=presnd(iw)
        enddo
        scomp=scomp+num_materials_vel

        do iw=1,num_materials_vel
         writend(scomp+iw)=divnd(iw)
        enddo
        scomp=scomp+num_materials_vel

        do iw=1,num_materials_vel
         writend(scomp+iw)=divdatnd(iw)
        enddo
        scomp=scomp+num_materials_vel

        do iw=1,num_materials_vel
         writend(scomp+iw)=machnd(iw)
        enddo
        scomp=scomp+num_materials_vel

        do iw=1,nmat
         writend(scomp+iw)=vofnd(iw)
        enddo
        scomp=scomp+nmat

        do iw=1,nmat*(SDIM+1)
         writend(scomp+iw)=lsdistnd(iw)
        enddo
        scomp=scomp+nmat*(SDIM+1)

        if (num_state_base.eq.2) then
         call is_dissolution(idissolution)
        else
         print *,"num_state_base invalid"
         stop
        endif

        iw=0
        do im=1,nmat
         do istate=1,num_state_material ! den,T
          iw=iw+1 
          writend(scomp+iw)=dennd(iw)
          if ((idissolution.eq.1).and.(istate.eq.2)) then
           writend(scomp+iw)=writend(scomp+iw)-one
          endif
         enddo ! istate
        enddo ! im

        scomp=scomp+num_state_material*nmat

        if ((num_materials_viscoelastic.ge.1).and. &
            (num_materials_viscoelastic.le.nmat)) then
         do iw=1,elastic_ncomp
          writend(scomp+iw)=elasticnd(iw) 
         enddo
         scomp=scomp+elastic_ncomp
        else if (num_materials_viscoelastic.eq.0) then
         ! do nothing
        else
         print *,"num_materials_viscoelastic invalid"
         stop
        endif

        do iw=1,nmat
         writend(scomp+iw)=viscnd(iw)
        enddo
        scomp=scomp+nmat

        do iw=1,5*nmat
         writend(scomp+iw)=tracend(iw)
        enddo
        scomp=scomp+5*nmat


! pgf90 will automatically break up lines if they exceed 80 chars.
! a format must be specified.  e.g. '(D25.16)'

        if (debug_slice.eq.1) then
         print *,"debug_slice: scomp= ",scomp
         do iw=1,SDIM
          print *,"debug_slice: iw,i,j,k,writend ",iw,i,j,k,writend(iw)
         enddo
        endif
  
        do iw=1,scomp
         if (iw.lt.scomp) then
          write(11,'(D25.16)',ADVANCE="NO") writend(iw)
         else if (iw.eq.scomp) then
          write(11,'(D25.16)') writend(iw)
         else
          print *,"iw invalid"
          stop
         endif
        enddo

       enddo  ! i=igridlo(1),igridhi(1)
       enddo  ! j=igridlo(2),igridhi(2)
       enddo  ! k=igridlo(3),igridhi(3)  (main output loop to "AMR" grid)

       close(11)
      else if (do_plot.eq.0) then
       ! do nothing
      else
       print *,"do_plot invalid"
       stop
      endif

      if (do_slice.eq.1) then

       if ((slice_dir.ge.0).and.(slice_dir.lt.SDIM)) then

        igridlo(3)=0
        igridhi(3)=0
        do dir=1,SDIM
         igridlo(dir)=lo(dir)-1
         igridhi(dir)=hi(dir)+1
        enddo

        do i=igridlo(1),igridhi(1) 
        do j=igridlo(2),igridhi(2) 
        do k=igridlo(3),igridhi(3) 
           ! iproblo=0
         call gridsten(xsten,problo,i,j,k,iproblo,bfact,dx,nhalf)
         do dir=1,SDIM
          xposnd(dir)=xsten(0,dir)
         enddo

         do iw=1,SDIM
          plotfab(D_DECL(i,j,k),iw)=xposnd(iw)
         enddo

         do dir=1,nmat
          lsdistnd(dir)=lsdist(D_DECL(i,j,k),dir)
         enddo

         call get_primary_material(lsdistnd,nmat,im_crit)

         do dir=1,num_materials_vel*(SDIM+1)
          velnd(dir)=vel(D_DECL(i,j,k),dir)
         enddo

         dir=num_materials_vel*SDIM+1
         plotfab(D_DECL(i,j,k),2*SDIM+1)=velnd(dir) ! pressure from solver

         presnd(1)=pres(D_DECL(i,j,k),1)
         divnd(1)=div(D_DECL(i,j,k),1)

         plotfab(D_DECL(i,j,k),2*SDIM+2)=presnd(1)
         plotfab(D_DECL(i,j,k),2*SDIM+3)=divnd(1)

         do dir=1,num_state_material*nmat
          dennd(dir)=den(D_DECL(i,j,k),dir)
         enddo

         denslice=dennd((im_crit-1)*num_state_material+1)
         tempslice=dennd((im_crit-1)*num_state_material+2)
         call INTERNAL_material(denslice,tempslice,eslice, &
           fort_material_type(im_crit),im_crit)
         KEslice=denslice*eslice

         do dir=1,SDIM
          iw=dir
          plotfab(D_DECL(i,j,k),SDIM+dir)=velnd(iw)
          KEslice=KEslice+half*denslice*(velnd(iw)**2)
         enddo ! dir

         plotfab(D_DECL(i,j,k),2*AMREX_SPACEDIM+4)=denslice
         plotfab(D_DECL(i,j,k),2*AMREX_SPACEDIM+5)=tempslice
         plotfab(D_DECL(i,j,k),2*AMREX_SPACEDIM+6)=KEslice
        enddo ! k
        enddo ! j
        enddo ! i

        igridlo(3)=0
        igridhi(3)=0
        do dir=1,SDIM
         igridlo(dir)=lo(dir)-1
         igridhi(dir)=hi(dir)+1
        enddo

          ! iproblo=0
        inbox=1
        call gridsten(xsten_fablo,problo,lo(1)-1,lo(2)-1,lo(SDIM)-1, &
         iproblo,bfact,dx,nhalf)
        call gridsten(xsten_fabhi,problo,hi(1)+1,hi(2)+1,hi(SDIM)+1, &
         iproblo,bfact,dx,nhalf)
        do dir=1,SDIM
         xfablo(dir)=xsten_fablo(0,dir)
         xfabhi(dir)=xsten_fabhi(0,dir)
        enddo

        do dir=1,SDIM
         if (slice_dir+1.ne.dir) then
          if ((xslice(dir).lt.xfablo(dir)-VOFTOL*dx(dir)).or. &
              (xslice(dir).gt.xfabhi(dir)+VOFTOL*dx(dir))) then
           inbox=0
          else
              ! iproblo=0
           do i=lo(dir)-1,hi(dir)
            call gridsten1D(xsten1D,problo,i,iproblo,bfact,dx,dir,nhalf)
            xposnd(dir)=xsten1D(0)
            if ((xslice(dir).ge.xsten1D(0)-VOFTOL*dx(dir)).and. &
                (xslice(dir).le.xsten1D(2)+VOFTOL*dx(dir))) then
             igridlo(dir)=i
             igridhi(dir)=i+1
            endif
           enddo ! i
          endif ! in the box
         endif ! tangential direction
        enddo ! dir  (finding igridlo,igridhi a thin box)

        if (inbox.eq.1) then

         do dir=1,SDIM
          if (slice_dir+1.ne.dir) then
           if (igridhi(dir)-igridlo(dir).ne.1) then
            print *,"slice: igridhi or igridlo not initialized correctly"
            stop
           endif
          endif
         enddo ! dir

         dir=slice_dir+1
          ! nslice=domhi+3
         do i=-1,nslice-2
           ! iproblo=0
          call gridsten1D(xsten1D_finest,problo,i,iproblo, &
            bfact_finest,dxfinest,dir,nhalf)
          x1D=xsten1D_finest(0)

          if (debug_slice.eq.1) then
           print *,"x1D= ",x1D
          endif

          if ((x1D.ge.xfablo(dir)-VOFTOL*dx(dir)).and. &
              (x1D.le.xfabhi(dir)+VOFTOL*dx(dir))) then
           inbox=0

            ! find bounding interval in the dir=slice_dir+1 direction.
           do j=lo(dir)-1,hi(dir)
             ! iproblo=0
            call gridsten1D(xsten1D,problo,j,iproblo,bfact,dx,dir,nhalf)
            xposnd(dir)=xsten1D(0)
            if ((x1D.ge.xsten1D(0)-VOFTOL*dx(dir)).and. &
                (x1D.le.xsten1D(2)+VOFTOL*dx(dir))) then
             igridlo(dir)=j
             igridhi(dir)=j+1
             inbox=1
            endif
           enddo ! j

           if (inbox.eq.1) then

            do n=1,nstate_slice
             i1=igridlo(1)
             j1=igridlo(2)
             k1=igridlo(3)
             if (debug_slice.eq.1) then
              print *,"i1,j1,k1 ",i1,j1,k1
             endif
             xc(dir)=x1D
             theta(3)=zero
             do dir2=1,SDIM
              if (dir2.ne.dir) then
               xc(dir2)=xslice(dir2)
              endif
               ! iproblo=0
              call gridsten1D(xsten1DL,problo,igridlo(dir2), &
                iproblo,bfact,dx,dir2,nhalf)
              xstenlo(dir2)=xsten1DL(0)
              xstenhi(dir2)=xsten1DL(2)
       
              if (xc(dir2).le.xstenlo(dir2)) then
               theta(dir2)=zero
              else if (xc(dir2).ge.xstenhi(dir2)) then
               theta(dir2)=one
              else 
               theta(dir2)=(xc(dir2)-xstenlo(dir2))/ &
                 (xstenhi(dir2)-xstenlo(dir2))
              endif
             enddo  ! dir2
             do ii=0,1
             do jj=0,1
              psten(ii,jj)= &
               (one-theta(3))*plotfab(D_DECL(i1+ii,j1+jj,k1),n)+ &
               theta(3)*plotfab(D_DECL(i1+ii,j1+jj,k1+1),n)
             enddo
             enddo
             do ii=0,1
              psten1D(ii)=(one-theta(2))*psten(ii,0)+theta(2)*psten(ii,1)
             enddo
             slice_data((i+1)*nstate_slice+n)= &
              (one-theta(1))*psten1D(0)+theta(1)*psten1D(1) 
            enddo ! n

           else if (inbox.eq.0) then
            ! do nothing
           else
            print *,"inbox invalid"
            stop
           endif 

          endif ! x1D in the box?

         enddo ! i

        endif ! inbox=1

       else
        print *,"slice_dir invalid"
        stop
       endif
      else if (do_slice.eq.0) then
       ! do nothing
      else
       print *,"do_slice invalid"
       stop
      endif

       ! interpolation to uniform grid

      igridlo(3)=0
      igridhi(3)=0
      do dir=1,SDIM
       igridlo(dir)=lo(dir)
       igridhi(dir)=hi(dir)
      enddo

      do dir=1,SDIM
       visual_ncell(dir)=vishi(dir)-vislo(dir)+1
       if (visual_ncell(dir).ge.1) then
        visual_dx(dir)=(probhi(dir)-problo(dir))/visual_ncell(dir)
        if (visual_dx(dir).gt.zero) then
         ! do nothing
        else
         print *,"visual_dx invalid"
         stop
        endif
       else
        print *,"visual_ncell invalid"
        stop
       endif
      enddo ! dir=1..sdim

      do i=igridlo(1),igridhi(1) 
      do j=igridlo(2),igridhi(2) 
      do k=igridlo(3),igridhi(3) 
        ! iproblo=0
       call gridsten(xsten,problo,i,j,k,iproblo,bfact,dx,nhalf)
       icritlo(3)=0
       icrithi(3)=0
       do dir=1,SDIM
        icrit(dir)=NINT((xsten(0,dir)-problo(dir))/visual_dx(dir))
        icritlo(dir)=icrit(dir)-1
        xcrit(dir)=icritlo(dir)*visual_dx(dir)+problo(dir)
        do while (xcrit(dir).gt.xsten(-1,dir)-visual_dx(dir)*VOFTOL) 
         icritlo(dir)=icritlo(dir)-1
         xcrit(dir)=xcrit(dir)-visual_dx(dir)
        enddo
        icrithi(dir)=icrit(dir)+1
        xcrit(dir)=icrithi(dir)*visual_dx(dir)+problo(dir)
        do while (xcrit(dir).lt.xsten(1,dir)+visual_dx(dir)*VOFTOL) 
         icrithi(dir)=icrithi(dir)+1
         xcrit(dir)=xcrit(dir)+visual_dx(dir)
        enddo
       enddo ! dir=1..sdim
       do ic=icritlo(1),icrithi(1)
       do jc=icritlo(2),icrithi(2)
       do kc=icritlo(3),icrithi(3)
        dir=1
        xcrit(dir)=problo(dir)+ic*visual_dx(dir) 
        dir=2
        xcrit(dir)=problo(dir)+jc*visual_dx(dir) 
        if (SDIM.eq.3) then 
         dir=SDIM
         xcrit(dir)=problo(dir)+kc*visual_dx(dir) 
        endif
        inbox=1
        do dir=1,SDIM
         if ((xcrit(dir).lt.xsten(-1,dir)-visual_dx(dir)*VOFTOL).or. &
             (xcrit(dir).gt.xsten(1,dir)+visual_dx(dir)*VOFTOL)) then
          inbox=0
         endif
        enddo
        if (inbox.eq.1) then
         icell(1)=i
         icell(2)=j
         icell(3)=k

         do dir=1,SDIM
          localfab(dir)=xcrit(dir)
         enddo

         icell_lo(3)=0
         icell_hi(3)=0
         do dir=1,SDIM
          if (xcrit(dir).lt.xsten(0,dir)) then
           icell_lo(dir)=icell(dir)-1
           xstenlo(dir)=xsten(-2,dir)
           xstenhi(dir)=xsten(0,dir)
          else
           icell_lo(dir)=icell(dir)
           xstenlo(dir)=xsten(0,dir)
           xstenhi(dir)=xsten(2,dir)
          endif
          icell_hi(dir)=icell_lo(dir)+1
         enddo

         do dir=SDIM+1,vis_ncomp
          localfab(dir)=zero
         enddo

         sumweight=zero
         do iBL=icell_lo(1),icell_hi(1)
         do jBL=icell_lo(2),icell_hi(2)
         do kBL=icell_lo(3),icell_hi(3)
          localwt=one
          do dir=1,SDIM
           if (xstenhi(dir)-xstenlo(dir).gt.zero) then
            theta(dir)=(xcrit(dir)-xstenlo(dir))/ &
             (xstenhi(dir)-xstenlo(dir))
            if ((theta(dir).ge.-VOFTOL).and. &
                (theta(dir).le.one+VOFTOL)) then
             if (theta(dir).lt.zero) then
              theta(dir)=zero
             endif
             if (theta(dir).gt.one) then
              theta(dir)=one
             endif
             localwt=localwt*(one-theta(dir))
            else
             print *,"theta invalid"
             stop
            endif
           else
            print *,"xstenhi or xstenlo invalid"
            stop
           endif
          enddo ! dir=1..sdim
          do dir=1,SDIM
           localfab(SDIM+dir)=localfab(SDIM+dir)+ &
              localwt*vel(D_DECL(iBL,jBL,kBL),dir)
          enddo
          ! 1. \dot{gamma}
          ! 2. Tr(A) if viscoelastic
          !    \dot{gamma} o.t.
          ! 3. Tr(A) (liquid viscosity - etaS)/etaP  if FENE-CR+Carreau
          !    Tr(A) if FENE-CR
          !    \dot{gamma} o.t.
          ! 4. (3) * f(A)  if viscoelastic
          !    \dot{gamma} o.t.
          ! 5. vorticity magnitude.
          localfab(2*SDIM+1)=localfab(2*SDIM+1)+ &
              localwt*trace(D_DECL(iBL,jBL,kBL),5)
          do im=1,nmat
           localfab(2*SDIM+1+im)=localfab(2*SDIM+1+im)+ &
              localwt*lsdist(D_DECL(iBL,jBL,kBL),im)
          enddo

          sumweight=sumweight+localwt
         enddo !kBL
         enddo !jBL
         enddo !iBL
         if (sumweight.gt.zero) then
          do dir=1,vis_ncomp-SDIM
           localfab(SDIM+dir)=localfab(SDIM+dir)/sumweight
          enddo
         else
          print *,"sumweight invalid"
          stop
         endif
         if (bfact.eq.1) then
          ! do nothing
         else if (bfact.gt.1) then
          do dir=1,SDIM
           if (xcrit(dir).le.xsten(-1,dir)) then
            xcrit(dir)=xsten(-1,dir)+1.0E-14*visual_dx(dir)
           endif
           if (xcrit(dir).ge.xsten(1,dir)) then
            xcrit(dir)=xsten(1,dir)-1.0E-14*visual_dx(dir)
           endif
          enddo

          ! find stencil for surrounding spectral element.
          stenlo(3)=0
          stenhi(3)=0
          do dir2=1,SDIM
           dxelem=dx(dir2)*bfact
           stenlo_index=NINT( (xsten(0,dir2)-problo(dir2))/dxelem-half )
           stenlo(dir2)=bfact*stenlo_index
           stenhi(dir2)=stenlo(dir2)+bfact-1
           if ((icell(dir2).lt.stenlo(dir2)).or. &
               (icell(dir2).gt.stenhi(dir2))) then
            print *,"SEM stencil should contain cell"
            print *,"dir2,icell(dir2),stenlo(dir2),stenhi(dir2) ", &
             dir2,icell(dir2),stenlo(dir2),stenhi(dir2)
            print *,"xcrit(dir2) ",xcrit(dir2)
            print *,"xsten(-1,dir2) ",xsten(-1,dir2)
            print *,"xsten(0,dir2) ",xsten(0,dir2)
            print *,"xsten(1,dir2) ",xsten(1,dir2)
            stop
           endif
          enddo ! dir2=1..sdim

          iSEM=stenlo(1) 
          jSEM=stenlo(2) 
          kSEM=stenlo(SDIM) 

          ! lower left hand corner of element stencil.
          ! iproblo=0
          call gridstenND(xstenND,problo,iSEM,jSEM,kSEM,iproblo,bfact, &
            dx,nhalf)

          do dir2=1,SDIM

           xcrit(dir2)=xcrit(dir2)-xstenND(0,dir2)

           if ((xcrit(dir2).lt.-INTERP_TOL*dx(dir2)).or. &
               (xcrit(dir2).gt.(INTERP_TOL+bfact)*dx(dir2))) then
            print *,"xcrit out of bounds CELLGRID"
            stop
           endif

          enddo ! dir2=1..sdim

          ncomp_SEM=vis_ncomp-SDIM
          VORT_comp_SEM=SDIM+1

          do n=1,ncomp_SEM
           SEM_value(n)=zero
          enddo
          gridtype=0  ! ggg  (Gauss in all directions)
          do dir2=1,SDIM
           SEMhi(dir2)=bfact-1
          enddo

          allocate(SEMloc(D_DECL(0:SEMhi(1),0:SEMhi(2),0:SEMhi(SDIM)), &
                          ncomp_SEM))

          do iSEM=stenlo(1),stenhi(1)
          do jSEM=stenlo(2),stenhi(2)
          do kSEM=stenlo(3),stenhi(3)
           ilocal=iSEM-stenlo(1)
           jlocal=jSEM-stenlo(2)
           klocal=kSEM-stenlo(3)
           do n=1,SDIM
            SEMloc(D_DECL(ilocal,jlocal,klocal),n)= &
             vel(D_DECL(iSEM,jSEM,kSEM),n)
           enddo
           SEMloc(D_DECL(ilocal,jlocal,klocal),VORT_comp_SEM)= &
            trace(D_DECL(iSEM,jSEM,kSEM),5)
           do im=1,nmat
            SEMloc(D_DECL(ilocal,jlocal,klocal),VORT_comp_SEM+im)= &
             lsdist(D_DECL(iSEM,jSEM,kSEM),im)
           enddo
          enddo
          enddo
          enddo

          caller_id=4
          call SEM_INTERP_ELEMENT( &
           ncomp_SEM,bfact,gridtype, &
           SEMhi,dx,xcrit,SEMloc,SEM_value,caller_id)

          deallocate(SEMloc)

           ! we discard the SEM interpolated values of lsdist.
          do dir=SDIM+1,vis_ncomp-nmat
           localfab(dir)=SEM_value(dir-SDIM)
          enddo

         else
          print *,"bfact invalid155"
          stop
         endif

         do dir=1,vis_ncomp
          fabout(D_DECL(ic,jc,kc),dir)=localfab(dir)
         enddo

        else if (inbox.eq.0) then
         ! do nothing
        else
         print *,"inbox invalid"
         stop
        endif
       enddo ! kc
       enddo ! jc
       enddo ! ic (visit all the uniform cells that overlap the "AMR" cells)

      enddo ! igridlo(3) .. igridhi(3)
      enddo ! igridlo(2) .. igridhi(2)
      enddo ! igridlo(1) .. igridhi(1)

      deallocate(plotfab)
      deallocate(reconfab)
 
      return
      end subroutine FORT_CELLGRID


      subroutine FORT_MEMSTATUS(procnum)
      IMPLICIT NONE

      INTEGER_T procnum,i
      character*80 pscmd 
      character*80 grepcmd 
      character*80 catcmd 
      character*6 procstr6
      character*80 rm1proc16
      character*80 rm2proc21
      character*80 pscmdproc 
      character*80 grepcmdproc 
      character*80 catcmdproc 

      character*80 cpupscmd
      character*80 cpugrepcmd
      character*80 cpucatcmd
      character*6 cpuprocstr6
      character*80 cpurm1proc16
      character*80 cpurm2proc21
      character*80 cpupscmdproc
      character*80 cpugrepcmdproc
      character*80 cpucatcmdproc

      INTEGER_T sysret

      call FLUSH(6)  ! unit=6 screen

      write(procstr6,'(I6)')procnum
      do i=1,6
       if (procstr6(i:i).eq.' ') then
        procstr6(i:i)='0'
       endif
      enddo

      write(cpuprocstr6,'(I6)')procnum
      do i=1,6
       if (cpuprocstr6(i:i).eq.' ') then
        cpuprocstr6(i:i)='0'
       endif
      enddo

      write(rm1proc16,'(A10,A6)') 'rm memtemp',procstr6
      write(rm2proc21,'(A15,A6)') 'rm memtemptrunc',procstr6
      print *,"issuing command ",rm1proc16
      print *,"and issuing command ",rm2proc21
      pscmd='ps -eo pcpu,pmem,size,comm --sort=-size > memtemp'
      write(pscmdproc,'(A49,A6)') pscmd,procstr6
      print *,"and issuing command ",pscmdproc
      grepcmd='grep "" memtemp'
      write(grepcmdproc,'(A15,A6,A21,A6)')  &
         grepcmd,procstr6,' -m 10 > memtemptrunc',procstr6
      print *,"and issuing command ",grepcmdproc
      catcmd='cat memtemptrunc'
      write(catcmdproc,'(A16,A6)') catcmd,procstr6
      print *,"and issuing command ",catcmdproc
      write(cpurm1proc16,'(A10,A6)') 'rm cputemp',cpuprocstr6
      write(cpurm2proc21,'(A15,A6)') 'rm cputemptrunc',cpuprocstr6
      print *,"and issuing command ",cpurm1proc16
      print *,"and issuing command ",cpurm2proc21
      cpupscmd='ps -eo pcpu,pmem,size,comm --sort=-pcpu > cputemp'
      write(cpupscmdproc,'(A49,A6)') cpupscmd,cpuprocstr6
      print *,"and issuing command ",cpupscmdproc
      cpugrepcmd='grep "" cputemp'
      write(cpugrepcmdproc,'(A15,A6,A21,A6)')  &
         cpugrepcmd,cpuprocstr6,' -m 10 > cputemptrunc',cpuprocstr6
      print *,"and issuing command ",cpugrepcmdproc
      cpucatcmd='cat cputemptrunc'
      write(cpucatcmdproc,'(A16,A6)') cpucatcmd,cpuprocstr6
      print *,"and issuing command ",cpucatcmdproc

      sysret=0

#ifdef PGIFORTRAN
      call system(rm1proc16)
      call system(rm2proc21)
      call system(pscmdproc)
      call system(grepcmdproc)
      call system(catcmdproc)
      call system(cpurm1proc16)
      call system(cpurm2proc21)
      call system(cpupscmdproc)
      call system(cpugrepcmdproc)
      call system(cpucatcmdproc)
#else
      call execute_command_line(rm1proc16,exitstat=sysret)
      call execute_command_line(rm2proc21,exitstat=sysret)
      call execute_command_line(pscmdproc,exitstat=sysret)
      call execute_command_line(grepcmdproc,exitstat=sysret)
      call execute_command_line(catcmdproc,exitstat=sysret)
      call execute_command_line(cpurm1proc16,exitstat=sysret)
      call execute_command_line(cpurm2proc21,exitstat=sysret)
      call execute_command_line(cpupscmdproc,exitstat=sysret)
      call execute_command_line(cpugrepcmdproc,exitstat=sysret)
      call execute_command_line(cpucatcmdproc,exitstat=sysret)
#endif

      if (sysret.ne.0) then
       print *,"execute_command_line has sysret=",sysret
       stop
      endif

      call FLUSH(6)  ! unit=6 screen

      return
      end subroutine FORT_MEMSTATUS

      subroutine FORT_OUTPUTSLICE( &
       time,nsteps,sliceint,slice_data,nslice,nstate_slice, &
       visual_option)
      IMPLICIT NONE

      INTEGER_T visual_option
      REAL_T time
      INTEGER_T nsteps,nslice,nstate_slice,n,sliceint,strandid
      REAL_T slice_data(nslice*nstate_slice)
      INTEGER_T i
      character*6 stepstr
      character*11 sfilename

        ! nstate_slice=x,y,z,xvel,yvel,zvel,PMG,PEOS,DIV,den,Temp,KE
        ! (value of material with LS>0)
      if (nstate_slice.ne.2*SDIM+6) then
       print *,"nstate_slice invalid in outputslice"
       print *,"nstate_slice= ",nstate_slice
       stop
      endif 
      if (nslice.le.1) then
       print *,"nslice invalid"
       stop
      endif

      if (sliceint.le.0) then
       strandid=1
      else
       strandid=(nsteps/sliceint)+1
      endif

      write(stepstr,'(I6)') nsteps
      do i=1,6
       if (stepstr(i:i).eq.' ') then
        stepstr(i:i)='0'
       endif
      enddo

      write(sfilename,'(A5,A6)') 'slice',stepstr

      print *,"sfilename ",sfilename
      open(unit=11,file=sfilename)

      if (visual_option.eq.-1) then

       if (SDIM.eq.2) then
        write(11,*) '# x,y,xvel,yvel,PMG,PEOS,DIV,den,temp,KE'
       else if (SDIM.eq.3) then
        write(11,*) '# x,y,z,xvel,yvel,zvel,PMG,PEOS,DIV,den,temp,KE'
       else
        print *,"dimension bust"
        stop
       endif

      else if (visual_option.eq.-2) then

       write(11,*) '"CLSVOF data"'
       if (SDIM.eq.2) then
        write(11,*) 'VARIABLES="X","Y","U","V","PMG","PEOS","DIV","D","T","KE"'
       else if (SDIM.eq.3) then
        write(11,*) 'VARIABLES="X","Y","Z","U","V","W","PMG","PEOS","DIV","D","T","KE"'
       else
        print *,"dimension bust"
        stop
       endif
       write(11,*)'zone i=',nslice-1,' SOLUTIONTIME=',time, &
        ' STRANDID=',strandid
      else
       print *,"visual_option invalid"
       stop
      endif

      do i=-1,nslice-2
       do n=1,nstate_slice-1
        write(11,'(D25.16)',ADVANCE="NO")  &
          slice_data((i+1)*nstate_slice+n)
       enddo
       n=nstate_slice
       write(11,'(D25.16)') slice_data((i+1)*nstate_slice+n)
      enddo ! i


      close(11)
      return
      end

      subroutine FORT_COMBINEZONES( &
       total_number_grids, &
       grids_per_level_array, &
       levels_array, &
       bfact_array, &
       gridno_array, &
       gridlo_array, &
       gridhi_array, &
       finest_level, &
       nsteps, &
       num_levels, &
       time, &
       visual_option, &
       visual_revolve, &
       plotint, &
       nmat, &
       nparts, &
       nparts_def, &
       im_solid_map)

      use global_utility_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T nparts
      INTEGER_T nparts_def
      INTEGER_T im_solid_map(nparts_def)
      INTEGER_T total_number_grids
      INTEGER_T num_levels
      INTEGER_T grids_per_level_array(num_levels)
      INTEGER_T levels_array(total_number_grids)
      INTEGER_T bfact_array(total_number_grids)
      INTEGER_T gridno_array(total_number_grids)
      INTEGER_T gridlo_array(total_number_grids*SDIM)
      INTEGER_T gridhi_array(total_number_grids*SDIM)
      INTEGER_T finest_level
      INTEGER_T nsteps
      REAL_T time
      INTEGER_T visual_option
      INTEGER_T visual_revolve
      INTEGER_T plotint
      INTEGER_T nmat

      INTEGER_T strandid
      INTEGER_T nwrite


      character*3 levstr
      character*5 gridstr
      character*18 filename18
      character*80 rmcommand

      character*6 stepstr
      character*16 newfilename16

      INTEGER_T i,j,k,dir
      INTEGER_T ilev,igrid
      INTEGER_T lo(SDIM),hi(SDIM)

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

      if (nmat.ne.num_materials) then
       print *,"nmat bust"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.nmat)) then
       print *,"nparts invalid FORT_COMBINEZONES"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.nmat)) then
       print *,"nparts_def invalid FORT_COMBINEZONES"
       stop
      endif

      call get_nwrite(plot_sdim,nwrite)

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

       call zones_revolve( &
        plot_sdim, &
        total_number_grids, &
        grids_per_level_array, &
        levels_array, &
        bfact_array, &
        gridno_array, &
        gridlo_array, &
        gridhi_array, &
        finest_level, &
        nsteps, &
        num_levels, &
        time, &
        visual_option, &
        visual_revolve, &
        plotint, &
        nmat, &
        nparts, &
        nparts_def, &
        im_solid_map)
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

       write(newfilename16,'(A6,A6,A4)') 'nddata',stepstr,'.plt'

       print *,"newfilename16 ",newfilename16


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
       open(unit=11,file=newfilename16,form="unformatted",access="stream")

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

       ! Variable names: combinezones
       call dumpstring_headers(plot_sdim)

        ! Zones
       do iz_gb=1,nzones_gb
        ! Zone marker. Value = 299.0
        write(11) ZONEMARKER
        ! Zone name
        Zonename = "ZONE"
        call dumpstring(Zonename)

        if (plotint.le.0) then
         strandid=1
        else
         strandid=(nsteps/plotint)+1
        endif
 
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
        write(11) hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+2
        write(11) hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+2
        if (plot_sdim.eq.3) then
         write(11) hi_gb(iz_gb,plot_sdim)-lo_gb(iz_gb,plot_sdim)+2
        else if (plot_sdim.eq.2) then
         write(11) 1
        else
         print *,"plot_sdim invalid"
         stop
        endif
 
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
 
        allocate(zone_gb(iz_gb)% &
         var(nwrite,lo(1):hi(1)+1, &
             lo(2):hi(2)+1, &
             klo_plot:khi_plot))

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
        do k=klo_plot,khi_plot
        do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)+1
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
         do k=klo_plot,khi_plot
         do j=lo_gb(iz_gb,2),hi_gb(iz_gb,2)+1
         do i=lo_gb(iz_gb,1),hi_gb(iz_gb,1)+1
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
      end subroutine FORT_COMBINEZONES


      subroutine FORT_ZALESAKNODE( &
        xlo,dx, &
        u,domlo,domhi, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        level, &
        DIMS(u), &
        time)

      use global_utility_module
      use probf90_module

      IMPLICIT NONE
      INTEGER_T bfact,level
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T domlo(SDIM), domhi(SDIM)
      INTEGER_T DIMDEC(u)
      REAL_T u(DIMV(u),SDIM)
      REAL_T time
      REAL_T xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k,dir
      REAL_T xx(SDIM)
      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf

      nhalf=1

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
      call checkbound(fablo,fabhi,DIMS(u),1,-1,411) 

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       do dir=1,SDIM
        xx(dir)=xsten(0,dir)
       enddo
 
       if (probtype.eq.31) then
        call circleuu(u(D_DECL(i,j,k),1),xx(1),xx(2),xx(SDIM))
        call circlevv(u(D_DECL(i,j,k),2),xx(1),xx(2),xx(SDIM))
        if (SDIM.eq.3) then
         call circleww(u(D_DECL(i,j,k),SDIM),xx(1),xx(2),xx(SDIM))
        endif
       else if (probtype.eq.29) then
        if (SDIM.eq.3) then
         call deform3duu(u(D_DECL(i,j,k),1),xx(1),xx(2),xx(SDIM),time,dx)
         call deform3dvv(u(D_DECL(i,j,k),2),xx(1),xx(2),xx(SDIM),time,dx)
         call deform3dww(u(D_DECL(i,j,k),SDIM),xx(1),xx(2),xx(SDIM),time,dx)
        else if (SDIM.eq.2) then
         call deformuu(u(D_DECL(i,j,k),1),xx(1),xx(2),time,dx)
         call deformvv(u(D_DECL(i,j,k),2),xx(1),xx(2),time,dx)
        else
         print *,"dimension bust"
         stop
        endif
       else if (probtype.eq.28) then
        call zalesakuu(u(D_DECL(i,j,k),1),xx(1),xx(2),xx(SDIM),time,dx)
        call zalesakvv(u(D_DECL(i,j,k),2),xx(1),xx(2),xx(SDIM),time,dx)
        if (SDIM.eq.3) then
         call zalesakww(u(D_DECL(i,j,k),SDIM),xx(1),xx(2),xx(SDIM),time,dx)
        endif
       else
        print *,"invalid probtype for ZALESAK"
       endif

       do dir=1,SDIM
        u(D_DECL(i,j,k),dir)=u(D_DECL(i,j,k),dir)/global_velocity_scale
       enddo

      enddo
      enddo
      enddo

      return 
      end subroutine FORT_ZALESAKNODE

       ! spectral_override==0 => always do low order
      subroutine FORT_AVGDOWN( &
       enable_spectral, &
       finest_level, &
       spectral_override, &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       mask,DIMS(mask), &
       lo,hi, &
       lof,hif) 

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: spectral_override
      INTEGER_T, intent(in) :: level_c
      INTEGER_T, intent(in) :: level_f
      INTEGER_T, intent(in) :: bfact_c
      INTEGER_T, intent(in) :: bfact_f
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: dxf(SDIM)
      REAL_T, intent(in) :: xlo_fine(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: DIMDEC(crse)
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      INTEGER_T, intent(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      INTEGER_T mstenlo(3), mstenhi(3)
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(out) :: crse(DIMV(crse),ncomp)
      REAL_T, intent(in) :: fine(DIMV(fine),ncomp)
      INTEGER_T flochi(SDIM)
      INTEGER_T ic,jc,kc
      INTEGER_T ifine,jfine,kfine
      INTEGER_T ilocal,jlocal,klocal
      INTEGER_T n
      INTEGER_T dir2
      INTEGER_T gridtype
      REAL_T voltotal
      REAL_T volall
      REAL_T wt(SDIM)
      REAL_T xsten(-1:1,SDIM)
      REAL_T xstenND(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T crse_value(ncomp)
      REAL_T dxelem_f
      REAL_T xcoarse(SDIM)
      INTEGER_T finelo_index
      REAL_T, dimension(D_DECL(:,:,:),:),allocatable :: ffine
      REAL_T INTERP_TOL
      INTEGER_T testmask,testmask2
      INTEGER_T nmat
      INTEGER_T local_enable_spectral
      INTEGER_T elem_test
      INTEGER_T caller_id

      caller_id=5
      nhalf=1

      INTERP_TOL=1.0E-4

      nmat=num_materials

      local_enable_spectral=enable_spectral

      if ((spectral_override.ne.0).and.(spectral_override.ne.1)) then
       print *,"spectral_override invalid"
       stop
      endif

      if (ncomp.lt.1) then
       print *,"ncomp invalid31"
       stop
      endif

      if (bfact_f.lt.1) then
       print *,"bfact_f invalid4 ",bfact_f
       stop
      endif
      if (bfact_c.lt.1) then
       print *,"bfact_c invalid4 ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level_f-1)) then
       print *,"level_c or level_f invalid"
       stop
      endif
      if (level_f.gt.finest_level) then
       print *,"level_f invalid"
       stop
      endif

      call checkbound(lo,hi,DIMS(crse),0,-1,411)
      call checkbound(lof,hif,DIMS(fine),0,-1,411)
      call checkbound(lof,hif,DIMS(mask),1,-1,411)

      gridtype=0  ! ggg  (Gauss in all directions)

      do dir2=1,SDIM
       flochi(dir2)=bfact_f-1
       if (2*lo(dir2).ne.lof(dir2)) then
        print *,"2*lo(dir2).ne.lof(dir2)"
        stop
       endif
       if (2*(hi(dir2)+1).ne.hif(dir2)+1) then
        print *,"2*(hi(dir2)+1).ne.hif(dir2)+1"
        stop
       endif
      enddo ! dir2=1..sdim

      allocate(ffine(D_DECL(0:flochi(1),0:flochi(2),0:flochi(SDIM)),ncomp))

      do dir2=1,SDIM
       elem_test=lo(dir2)/bfact_c

       if (elem_test*bfact_c.ne.lo(dir2)) then
         print *,"elem_test invalid (loc) 1"
         print *,"elem_test=",elem_test
         print *,"bfact_c=",bfact_c
         print *,"dir2=",dir2
         print *,"lo(dir2)=",lo(dir2)
         print *,"hi(dir2)=",hi(dir2)
         print *,"lof(dir2)=",lof(dir2)
         print *,"hif(dir2)=",hif(dir2)
         stop
       endif
       elem_test=(hi(dir2)+1)/bfact_c
       if (elem_test*bfact_c.ne.hi(dir2)+1) then
         print *,"elem_test invalid (hic) 1"
         print *,"elem_test=",elem_test
         print *,"bfact_c=",bfact_c
         print *,"dir2=",dir2
         print *,"lo(dir2)=",lo(dir2)
         print *,"hi(dir2)=",hi(dir2)
         print *,"lof(dir2)=",lof(dir2)
         print *,"hif(dir2)=",hif(dir2)
         stop
       endif
       elem_test=lof(dir2)/bfact_f
       if (elem_test*bfact_f.ne.lof(dir2)) then
        print *,"elem_test invalid (lof) 1"
        print *,"elem_test=",elem_test
        print *,"bfact_f=",bfact_f
        print *,"dir2=",dir2
        print *,"lo(dir2)=",lo(dir2)
        print *,"hi(dir2)=",hi(dir2)
        print *,"lof(dir2)=",lof(dir2)
        print *,"hif(dir2)=",hif(dir2)
        stop
       endif
       elem_test=(hif(dir2)+1)/bfact_f
       if (elem_test*bfact_f.ne.hif(dir2)+1) then
        print *,"elem_test invalid (hif) 1"
        print *,"elem_test=",elem_test
        print *,"bfact_f=",bfact_f
        print *,"dir2=",dir2
        print *,"lo(dir2)=",lo(dir2)
        print *,"hi(dir2)=",hi(dir2)
        print *,"lof(dir2)=",lof(dir2)
        print *,"hif(dir2)=",hif(dir2)
        stop
       endif
      enddo ! dir2=1..sdim

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)

       call gridsten_level(xsten,ic,jc,kc,level_c,nhalf)

        ! find stencil for surrounding fine level spectral element.
       stenlo(3)=0
       stenhi(3)=0
       do dir2=1,SDIM
        dxelem_f=dxf(dir2)*bfact_f
        xcoarse(dir2)=xsten(0,dir2)
        finelo_index=NINT( (xcoarse(dir2)-problo(dir2))/dxelem_f-half )
        stenlo(dir2)=bfact_f*finelo_index
        stenhi(dir2)=stenlo(dir2)+bfact_f-1
        if ((stenlo(dir2).lt.lof(dir2)).or. &
            (stenhi(dir2).gt.hif(dir2))) then
         print *,"stenlo or stenhi invalid 1"
         stop
        endif
       enddo ! dir2

       mstenlo(3)=0
       mstenhi(3)=0
       do dir2=1,SDIM
        mstenlo(dir2)=stenlo(dir2)
        mstenhi(dir2)=stenhi(dir2)
        if (stenhi(dir2).eq.stenlo(dir2)) then
         mstenlo(dir2)=stenlo(dir2)-1
         mstenhi(dir2)=stenhi(dir2)+1
        endif
       enddo ! dir2

       ifine=stenlo(1) 
       jfine=stenlo(2) 
       kfine=stenlo(SDIM) 

        ! lower left hand corner of element stencil.
       call gridstenND_level(xstenND,ifine,jfine,kfine,level_f,nhalf)

       do dir2=1,SDIM

        xcoarse(dir2)=xcoarse(dir2)-xstenND(0,dir2)

        if (stenhi(dir2).ne.stenlo(dir2)) then
         if (abs(xcoarse(dir2)).lt.INTERP_TOL*dxf(dir2)) then
          mstenlo(dir2)=stenlo(dir2)-1
         endif
         if (abs(xcoarse(dir2)-bfact_f*dxf(dir2)).lt. &
             INTERP_TOL*dxf(dir2)) then
          mstenhi(dir2)=stenhi(dir2)+1
         endif
        endif ! stenhi<>stenlo

        if ((xcoarse(dir2).lt.-INTERP_TOL*dxf(dir2)).or. &
            (xcoarse(dir2).gt.(INTERP_TOL+bfact_f)*dxf(dir2))) then
         print *,"xcoarse out of bounds avgdown"
         stop
        endif

       enddo ! dir2=1..sdim

       ifine=mstenlo(1) 
       jfine=mstenlo(2) 
       kfine=mstenlo(SDIM) 
       testmask=NINT(mask(D_DECL(ifine,jfine,kfine)))
       do ifine=mstenlo(1),mstenhi(1) 
       do jfine=mstenlo(2),mstenhi(2) 
       do kfine=mstenlo(3),mstenhi(3) 
        testmask2=NINT(mask(D_DECL(ifine,jfine,kfine)))
        if (testmask2.ne.testmask) then
         testmask=0
        endif
       enddo 
       enddo 
       enddo 

       do n=1,ncomp
        crse_value(n)=zero
       enddo
       voltotal=zero

       if ((bfact_f.ge.2).and. &
           ((local_enable_spectral.eq.1).or. &
            (local_enable_spectral.eq.2)).and. &
           (testmask.ge.1).and. &
           (testmask.le.nmat).and. &
           (spectral_override.eq.1)) then

        do ifine=stenlo(1),stenhi(1)
        do jfine=stenlo(2),stenhi(2)
        do kfine=stenlo(3),stenhi(3)
         ilocal=ifine-stenlo(1)
         jlocal=jfine-stenlo(2)
         klocal=kfine-stenlo(3)
         do n=1,ncomp
          ffine(D_DECL(ilocal,jlocal,klocal),n)= &
            fine(D_DECL(ifine,jfine,kfine),n)
         enddo
        enddo
        enddo
        enddo

        call SEM_INTERP_ELEMENT( &
         ncomp,bfact_f,gridtype, &
         flochi,dxf,xcoarse,ffine,crse_value,caller_id)

        voltotal=one

       else if ((bfact_f.eq.1).or. &
                (local_enable_spectral.eq.0).or. &
                (local_enable_spectral.eq.3).or. &
                (testmask.eq.0).or. &
                (spectral_override.eq.0)) then

        call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
         bfact_c,bfact_f)

        do ifine=stenlo(1),stenhi(1)
         call intersect_weight_avg(ic,ifine,bfact_c,bfact_f,wt(1))
         if (wt(1).gt.zero) then
          do jfine=stenlo(2),stenhi(2)
           call intersect_weight_avg(jc,jfine,bfact_c,bfact_f,wt(2))
           if (wt(2).gt.zero) then
            do kfine=stenlo(3),stenhi(3)
             if (SDIM.eq.3) then
              call intersect_weight_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
             endif
             if (wt(SDIM).gt.zero) then
              volall=wt(1)
              do dir2=2,SDIM
               volall=volall*wt(dir2)
              enddo
              do n = 1, ncomp
               crse_value(n) = crse_value(n) +  &
                volall*fine(D_DECL(ifine,jfine,kfine),n)
              enddo
              voltotal=voltotal+volall
             endif ! wt(sdim).gt.0
            enddo ! kfine
           endif
          enddo ! jfine
         endif
        enddo ! ifine

       else
        print *,"parameter bust in avgdown"
        stop
       endif

       if (voltotal.le.zero) then
        print *,"voltotal invalid"
        stop
       endif

       do n = 1,ncomp
        crse(D_DECL(ic,jc,kc),n)=crse_value(n)/voltotal
       enddo

      enddo
      enddo
      enddo ! ic,jc,kc

      deallocate(ffine)

      return
      end subroutine FORT_AVGDOWN


      subroutine FORT_AVGDOWN_COPY( &
       enable_spectral, &
       finest_level, &
       operation_flag, &
       dir, &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp_flux, &
       ncomp_den, &
       ncomp_vel, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       den_fine,DIMS(den_fine), &
       vel_fine,DIMS(vel_fine), &
       mask,DIMS(mask), &
       fine_LS,DIMS(fine_LS), &
       lo,hi, &
       lof,hif) 

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probcommon_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: operation_flag
      INTEGER_T, intent(in) :: dir
      INTEGER_T, intent(in) :: level_c
      INTEGER_T, intent(in) :: level_f
      INTEGER_T, intent(in) :: bfact_c
      INTEGER_T, intent(in) :: bfact_f
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: dxf(SDIM)
      REAL_T, intent(in) :: xlo_fine(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: ncomp_flux
      INTEGER_T, intent(in) :: ncomp_den
      INTEGER_T, intent(in) :: ncomp_vel
      INTEGER_T, intent(in) :: DIMDEC(crse)
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: DIMDEC(den_fine)
      INTEGER_T, intent(in) :: DIMDEC(vel_fine)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: DIMDEC(fine_LS)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      INTEGER_T, intent(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      INTEGER_T mstenlo(3),mstenhi(3)
      REAL_T, intent(in) :: fine_LS(DIMV(fine_LS),num_materials)
      REAL_T, intent(in) :: mask(DIMV(mask))
      REAL_T, intent(out) :: crse(DIMV(crse),ncomp_flux)
      REAL_T, intent(in) :: fine(DIMV(fine),ncomp_flux)
      REAL_T, intent(in) :: den_fine(DIMV(den_fine),ncomp_den)
      REAL_T, intent(in) :: vel_fine(DIMV(vel_fine),ncomp_vel)
      INTEGER_T flochi(SDIM)
      INTEGER_T ic,jc,kc
      INTEGER_T ifine,jfine,kfine
      INTEGER_T ilocal,jlocal,klocal
      INTEGER_T istrip,jstrip,kstrip
      INTEGER_T iface,jface,kface
      INTEGER_T ii,jj,kk
      INTEGER_T n
      INTEGER_T dir2
      INTEGER_T side
      INTEGER_T gridtype
      REAL_T voltotal
      REAL_T volall
      REAL_T wt(SDIM)
      REAL_T xsten(-1:1,SDIM)
      REAL_T xsten_fine(-1:1,SDIM)
      REAL_T xstenND(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T crse_value(ncomp_flux)
      REAL_T dxelem_f
      REAL_T xcoarse(SDIM)
      INTEGER_T finelo_index
      REAL_T, dimension(D_DECL(:,:,:),:),allocatable :: ffine
      REAL_T INTERP_TOL
      INTEGER_T testmask,testmask2
      INTEGER_T nmat
      INTEGER_T elem_test
      INTEGER_T im,imcrit
      REAL_T fine_data
      INTEGER_T dencomp
      REAL_T LS_local(num_materials)
      INTEGER_T mat_freq(num_materials)
      INTEGER_T caller_id

      caller_id=6

      nhalf=1

      INTERP_TOL=1.0E-4

      nmat=num_materials

      if (operation_flag.eq.7) then ! advection
       if (ncomp_vel.ne.SDIM) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.nmat*num_state_material) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.SDIM+num_state_base) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.1) then ! P cell to P MAC
       if (ncomp_vel.ne.1) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.1) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if ((operation_flag.eq.3).or. & !u cell to MAC
               (operation_flag.eq.5).or. & !UMAC=UMAC+beta diff_reg
               (operation_flag.eq.10).or. &
               (operation_flag.eq.11)) then 
       if (ncomp_vel.ne.AMREX_SPACEDIM) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.AMREX_SPACEDIM) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.9) then
       if (ncomp_vel.ne.nmat*num_state_material) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.nmat*num_state_material) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.0) then
       if (ncomp_vel.ne.1) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.1) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.8) then ! viscosity
       if (ncomp_vel.ne.SDIM) then
        print *,"ncomp_vel invalid visc"
        stop
       endif
       if (ncomp_den.ne.SDIM) then
        print *,"ncomp_den invalid visc"
        stop
       endif
       if (ncomp_flux.ne.SDIM) then
        print *,"ncomp_flux invalid visc"
        stop
       endif
      else
       print *,"operation_flag invalid visc"
       stop
      endif

      if (bfact_f.lt.1) then
       print *,"bfact_f invalid4 ",bfact_f
       stop
      endif
      if (bfact_c.lt.1) then
       print *,"bfact_c invalid4 ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level_f-1)) then
       print *,"level_c or level_f invalid"
       stop
      endif
      if (level_f.gt.finest_level) then
       print *,"level_f invalid"
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
       print *,"dir invalid"
       stop
      endif

      call checkbound(lo,hi,DIMS(crse),0,dir,411)
      call checkbound(lof,hif,DIMS(fine),0,dir,411)
      call checkbound(lof,hif,DIMS(mask),1,-1,411)
      call checkbound(lof,hif,DIMS(fine_LS),1,-1,411)
      call checkbound(lof,hif,DIMS(den_fine),1,-1,411)
      call checkbound(lof,hif,DIMS(vel_fine),1,-1,411)

      gridtype=0  ! ggg  (Gauss in all directions)

      do dir2=1,SDIM
       flochi(dir2)=bfact_f-1
       if (2*lo(dir2).ne.lof(dir2)) then
        print *,"2*lo(dir2).ne.lof(dir2)"
        stop
       endif
       if (2*(hi(dir2)+1).ne.hif(dir2)+1) then
        print *,"2*(hi(dir2)+1).ne.hif(dir2)+1"
        stop
       endif
      enddo ! dir2=1..sdim

      allocate(ffine(D_DECL(0:flochi(1),0:flochi(2),0:flochi(SDIM)), &
                     ncomp_flux))

      do dir2=1,SDIM
       elem_test=lo(dir2)/bfact_c
       if (elem_test*bfact_c.ne.lo(dir2)) then
         print *,"elem_test invalid (loc) 2"
         print *,"elem_test=",elem_test
         print *,"bfact_c=",bfact_c
         print *,"dir2=",dir2
         print *,"lo(dir2)=",lo(dir2)
         stop
       endif
       elem_test=(hi(dir2)+1)/bfact_c
       if (elem_test*bfact_c.ne.hi(dir2)+1) then
         print *,"elem_test invalid (hic) 2"
         print *,"elem_test=",elem_test
         print *,"bfact_c=",bfact_c
         print *,"dir2=",dir2
         print *,"hi(dir2)=",hi(dir2)
         stop
       endif
       elem_test=lof(dir2)/bfact_f
       if (elem_test*bfact_f.ne.lof(dir2)) then
        print *,"elem_test invalid (lof) 2"
        print *,"elem_test=",elem_test
        print *,"bfact_f=",bfact_f
        print *,"dir2=",dir2
        print *,"lof(dir2)=",lof(dir2)
        stop
       endif
       elem_test=(hif(dir2)+1)/bfact_f
       if (elem_test*bfact_f.ne.hif(dir2)+1) then
        print *,"elem_test invalid (hif) 2"
        print *,"elem_test=",elem_test
        print *,"bfact_f=",bfact_f
        print *,"dir2=",dir2
        print *,"hif(dir2)=",hif(dir2)
        stop
       endif
      enddo ! dir2=1..sdim

       ! coarse dimensions
      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do side=1,2

        ! strip of coarse cells below the fine grid and neighboring
        ! coarse region.
       if (side.eq.1) then
        growlo(dir+1)=lo(dir+1)
        growhi(dir+1)=lo(dir+1)
       else if (side.eq.2) then
        growlo(dir+1)=hi(dir+1)
        growhi(dir+1)=hi(dir+1)
       else
        print *,"side invalid"
        stop
       endif
 
       do ic=growlo(1),growhi(1)
       do jc=growlo(2),growhi(2)
       do kc=growlo(3),growhi(3)

        call gridsten_level(xsten,ic,jc,kc,level_c,nhalf)

         ! index for fine cell that is next to the fine grid boundary.
         ! (fine strip of cells below the fine grid, and neighboring a 
         ! coarse region)
        if (side.eq.1) then
         ifine=2*ic
         jfine=2*jc
         kfine=2*kc
         iface=ic
         jface=jc
         kface=kc
        else if (side.eq.2) then
         ifine=2*ic+1
         jfine=2*jc+1
         kfine=2*kc+1
         iface=ic+ii
         jface=jc+jj
         kface=kc+kk
        else
         print *,"side invalid"
         stop
        endif
        call gridsten_level(xsten_fine,ifine,jfine,kfine,level_f,nhalf)
        xsten(0,dir+1)=xsten_fine(0,dir+1)

         ! find stencil for surrounding fine level spectral element.
        stenlo(3)=0
        stenhi(3)=0
        do dir2=1,SDIM
         dxelem_f=dxf(dir2)*bfact_f
          ! in: AVGDOWN_COPY
          ! xsten is target location of where to interpolate to.
          ! In the normal direction, xsten corresponds to the 
          ! covering fine level coordinate (copy).  In the tangential
          ! direction, xsten corresponds to a coarse level
          ! coordinate (avgdown). 
         xcoarse(dir2)=xsten(0,dir2)
         finelo_index=NINT( (xcoarse(dir2)-problo(dir2))/dxelem_f-half )
         stenlo(dir2)=bfact_f*finelo_index
         stenhi(dir2)=stenlo(dir2)+bfact_f-1
         if ((stenlo(dir2).lt.lof(dir2)).or. &
             (stenhi(dir2).gt.hif(dir2))) then
          print *,"stenlo or stenhi invalid 2"
          stop
         endif
        enddo ! dir2=1..sdim

        mstenlo(3)=0
        mstenhi(3)=0
        do dir2=1,SDIM
         mstenlo(dir2)=stenlo(dir2)
         mstenhi(dir2)=stenhi(dir2)
         if (stenhi(dir2).eq.stenlo(dir2)) then
          if (bfact_f.eq.1) then
           mstenlo(dir2)=stenlo(dir2)-1
           mstenhi(dir2)=stenhi(dir2)+1
          else
           print *,"bfact_f invalid"
           stop
          endif
         else if (stenhi(dir2).gt.stenlo(dir2)) then
          if (bfact_f.gt.1) then
           ! do nothing
          else
           print *,"bfact_f invalid"
           stop
          endif
         else
          print *,"stenlo or stenhi bust"
          stop
         endif
        enddo ! dir2=1..sdim

        ifine=stenlo(1) 
        jfine=stenlo(2) 
        kfine=stenlo(SDIM) 

         ! lower left hand corner of element stencil.
        call gridstenND_level(xstenND,ifine,jfine,kfine,level_f,nhalf)

        do dir2=1,SDIM

         xcoarse(dir2)=xcoarse(dir2)-xstenND(0,dir2)

         if (stenhi(dir2).gt.stenlo(dir2)) then
          if (bfact_f.gt.1) then
           if (abs(xcoarse(dir2)).lt.INTERP_TOL*dxf(dir2)) then
            if (dir2.ne.dir+1) then
             mstenlo(dir2)=stenlo(dir2)-1
            else
             print *,"xcoarse invalid in normal direction to face"
             stop
            endif
           endif
           if (abs(xcoarse(dir2)-bfact_f*dxf(dir2)).lt. &
              INTERP_TOL*dxf(dir2)) then
            if (dir2.ne.dir+1) then
             mstenhi(dir2)=stenhi(dir2)+1
            else
             print *,"xcoarse invalid in normal direction to face"
             stop
            endif
           endif
          else
           print *,"bfact_f invalid"
           stop
          endif
         else if (stenhi(dir2).eq.stenlo(dir2)) then
          if (bfact_f.eq.1) then
           ! do nothing
          else
           print *,"bfact_f invalid"
           stop
          endif
         else
          print *,"stenlo or stenhi bust"
          stop
         endif 

         if ((xcoarse(dir2).lt.-INTERP_TOL*dxf(dir2)).or. &
             (xcoarse(dir2).gt.(INTERP_TOL+bfact_f)*dxf(dir2))) then
          print *,"xcoarse out of bounds avgdown_copy"
          stop
         endif

        enddo ! dir2=1..sdim

        do im=1,nmat
         mat_freq(im)=0
        enddo
        do ifine=mstenlo(1),mstenhi(1) 
        do jfine=mstenlo(2),mstenhi(2) 
        do kfine=mstenlo(3),mstenhi(3) 
         do im=1,nmat
          LS_local(im)=fine_LS(D_DECL(ifine,jfine,kfine),im)
         enddo
         call get_primary_material(LS_local,nmat,imcrit)
         if ((imcrit.ge.1).and.(imcrit.le.nmat)) then
          mat_freq(imcrit)=mat_freq(imcrit)+1
         else
          print *,"imcrit invalid"
          stop
         endif
        enddo 
        enddo 
        enddo 
        imcrit=0
        do im=1,nmat
         if (mat_freq(im).gt.0) then
          if (imcrit.eq.0) then
           imcrit=im
          else if ((imcrit.ge.1).and.(imcrit.le.nmat)) then
           if (mat_freq(im).gt.mat_freq(imcrit)) then
            imcrit=im
           endif
          else
           print *,"imcrit invalid"
           stop
          endif
         else if (mat_freq(im).eq.0) then
          ! do nothing
         else
          print *,"mat_freq(im) invalid"
          stop
         endif
        enddo ! im=1..nmat

        ifine=mstenlo(1) 
        jfine=mstenlo(2) 
        kfine=mstenlo(SDIM) 
        testmask=NINT(mask(D_DECL(ifine,jfine,kfine)))
        do ifine=mstenlo(1),mstenhi(1) 
        do jfine=mstenlo(2),mstenhi(2) 
        do kfine=mstenlo(3),mstenhi(3) 
         testmask2=NINT(mask(D_DECL(ifine,jfine,kfine)))
         if (testmask2.ne.testmask) then
          testmask=0
         endif
        enddo 
        enddo 
        enddo 

        if ((imcrit.ge.1).and.(imcrit.le.nmat)) then

         do n=1,ncomp_flux
          crse_value(n)=zero
         enddo
         voltotal=zero

         if ((bfact_f.ge.2).and. &
             ((enable_spectral.eq.1).or. &
              (enable_spectral.eq.2)).and. &
             (testmask.ge.1).and. &
             (testmask.le.nmat)) then

          do ifine=stenlo(1),stenhi(1)
          do jfine=stenlo(2),stenhi(2)
          do kfine=stenlo(3),stenhi(3)
           ilocal=ifine-stenlo(1)
           jlocal=jfine-stenlo(2)
           klocal=kfine-stenlo(3)

           istrip=ifine
           jstrip=jfine
           kstrip=kfine
           if (dir.eq.0) then
            if (side.eq.1) then
             istrip=stenlo(dir+1)
            else if (side.eq.2) then
             istrip=stenhi(dir+1)
            else
             print *,"side invalid"
             stop
            endif
           else if (dir.eq.1) then
            if (side.eq.1) then
             jstrip=stenlo(dir+1)
            else if (side.eq.2) then
             jstrip=stenhi(dir+1)
            else
             print *,"side invalid"
             stop
            endif
           else if ((dir.eq.2).and.(SDIM.eq.3)) then
            if (side.eq.1) then
             kstrip=stenlo(dir+1)
            else if (side.eq.2) then
             kstrip=stenhi(dir+1)
            else
             print *,"side invalid"
             stop
            endif
           else
            print *,"dir invalid"
            stop
           endif

           do n=1,ncomp_flux

            if (operation_flag.eq.8) then ! viscosity
             if ((n.ge.1).and.(n.le.SDIM)) then
              fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),n)
             else
              print *,"n invalid"
              stop
             endif
            else if (operation_flag.eq.7) then ! advection
             if ((n.ge.1).and.(n.le.SDIM)) then
              fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),n)
             else if ((n.gt.SDIM).and.(n.le.SDIM+num_state_base)) then
              dencomp=(imcrit-1)*num_state_material+n-SDIM
              fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),dencomp)
             else
              print *,"n invalid"
              stop
             endif
            else if (operation_flag.eq.1) then ! p cell to MAC
             if (n.eq.1) then
              fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),n)
             else
              print *,"n invalid"
              stop
             endif
            else if ((operation_flag.eq.3).or. & !u cell to MAC
                     (operation_flag.eq.5).or. & !UMAC=UMAC+beta diff_reg
                     (operation_flag.eq.10).or. &
                     (operation_flag.eq.11)) then 
             if ((dir.ge.0).and.(dir.lt.AMREX_SPACEDIM)) then
              fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),dir+1)
             else
              print *,"dir invalid"
              stop
             endif
            else if (operation_flag.eq.9) then
             if (n.eq.1) then
              dencomp=(imcrit-1)*num_state_material+1
              fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),dencomp)
             else
              print *,"n invalid"
              stop
             endif
            else if (operation_flag.eq.0) then
             if (n.eq.1) then
              fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),1)
             else
              print *,"n invalid"
              stop
             endif
            else
             print *,"operation_flag invalid"
             stop
            endif

            ffine(D_DECL(ilocal,jlocal,klocal),n)=fine_data
           enddo ! n=1,ncomp_flux
          enddo ! kfine
          enddo ! jfine
          enddo ! ifine

          call SEM_INTERP_ELEMENT( &
           ncomp_flux,bfact_f,gridtype, &
           flochi,dxf,xcoarse,ffine,crse_value,caller_id)

          voltotal=one

         else if ((bfact_f.eq.1).or. &
                  (enable_spectral.eq.0).or. &
                  (enable_spectral.eq.3).or. & ! SEM time only
                  (testmask.eq.0)) then

          call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
           bfact_c,bfact_f)
          if (side.eq.1) then
           stenlo(dir+1)=lof(dir+1)
           stenhi(dir+1)=lof(dir+1)
          else if (side.eq.2) then
           stenlo(dir+1)=hif(dir+1)
           stenhi(dir+1)=hif(dir+1)
          else
           print *,"side invalid"
           stop
          endif

          do ifine=stenlo(1),stenhi(1)
           call intersect_weight_avg_COPY(ic,ifine,bfact_c,bfact_f,wt(1), &
             dir,0)
           if (wt(1).gt.zero) then
            do jfine=stenlo(2),stenhi(2)
             call intersect_weight_avg_COPY(jc,jfine,bfact_c,bfact_f,wt(2), &
               dir,1)
             if (wt(2).gt.zero) then
              do kfine=stenlo(3),stenhi(3)
               if (SDIM.eq.3) then
                call intersect_weight_avg_COPY(kc,kfine,bfact_c,bfact_f, &
                  wt(SDIM),dir,2)
               endif
               if (wt(SDIM).gt.zero) then
                volall=wt(1)
                do dir2=2,SDIM
                 volall=volall*wt(dir2)
                enddo

                istrip=ifine
                jstrip=jfine
                kstrip=kfine
                if (dir.eq.0) then
                 if (side.eq.1) then
                  istrip=stenlo(dir+1)
                 else if (side.eq.2) then
                  istrip=stenhi(dir+1)
                 else
                  print *,"side invalid"
                  stop
                 endif
                else if (dir.eq.1) then
                 if (side.eq.1) then
                  jstrip=stenlo(dir+1)
                 else if (side.eq.2) then
                  jstrip=stenhi(dir+1)
                 else
                  print *,"side invalid"
                  stop
                 endif
                else if ((dir.eq.2).and.(SDIM.eq.3)) then
                 if (side.eq.1) then
                  kstrip=stenlo(dir+1)
                 else if (side.eq.2) then
                  kstrip=stenhi(dir+1)
                 else
                  print *,"side invalid"
                  stop
                 endif
                else
                 print *,"dir invalid"
                 stop
                endif

                do n = 1, ncomp_flux

                 if (operation_flag.eq.8) then ! viscosity
                  if ((n.ge.1).and.(n.le.SDIM)) then
                   fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),n)
                  else
                   print *,"n invalid"
                   stop
                  endif
                 else if (operation_flag.eq.7) then ! advection
                  if ((n.ge.1).and.(n.le.SDIM)) then
                   fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),n)
                  else if ((n.gt.SDIM).and.(n.le.SDIM+num_state_base)) then
                   dencomp=(imcrit-1)*num_state_material+n-SDIM
                   fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),dencomp)
                  else
                   print *,"n invalid"
                   stop
                  endif
                 else if (operation_flag.eq.1) then ! p cell to MAC
                  if (n.eq.1) then
                   fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),n)
                  else
                   print *,"n invalid"
                   stop
                  endif
                 else if ((operation_flag.eq.3).or. & !u cell to MAC
                          (operation_flag.eq.5).or. & !UMAC=UMAC+beta diff_reg
                          (operation_flag.eq.10).or. &
                          (operation_flag.eq.11)) then 
                  if ((dir.ge.0).and.(dir.lt.AMREX_SPACEDIM)) then
                   fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),dir+1)
                  else
                   print *,"dir invalid"
                   stop
                  endif
                 else if (operation_flag.eq.9) then ! den_TO_MAC
                  if (n.eq.1) then
                   dencomp=(imcrit-1)*num_state_material+1
                   fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),dencomp)
                  else
                   print *,"n invalid"
                   stop
                  endif
                 else if (operation_flag.eq.0) then
                  if (n.eq.1) then
                   fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),1)
                  else
                   print *,"n invalid"
                   stop
                  endif
                 else
                  print *,"operation_flag invalid"
                  stop
                 endif

                 crse_value(n) = crse_value(n) + volall*fine_data
                enddo ! n=1..ncomp_flux

                voltotal=voltotal+volall
               endif ! wt(sdim).gt.0
              enddo ! kfine
             endif
            enddo ! jfine
           endif
          enddo ! ifine

         else
          print *,"parameter bust in avgdown_copy"
          stop
         endif

         if (voltotal.gt.zero) then
          do n = 1,ncomp_flux
           crse(D_DECL(iface,jface,kface),n)=crse_value(n)/voltotal
          enddo
         else
          print *,"voltotal invalid"
          stop
         endif

        else
         print *,"imcrit invalid"
         stop
        endif
 
       enddo
       enddo
       enddo ! ic,jc,kc

      enddo ! side=1..2

      deallocate(ffine)

      return
      end subroutine FORT_AVGDOWN_COPY


      subroutine FORT_INTERP_COPY( &
       enable_spectral, &
       dxc,dx, &
       finest_level, &
       operation_flag, &
       tilelo,tilehi, &
       fablo,fabhi, &
       dir, &
       problo, &
       level_c,level, &
       bfact_c,bfact, &
       xlo, &
       ncomp_flux, &
       ncomp_den, &
       ncomp_vel, &
       fine,DIMS(fine), &
       den_crse,DIMS(den_crse), &
       vel_crse,DIMS(vel_crse), &
       masknbr,DIMS(masknbr), &
       masksem,DIMS(masksem), &
       cmasksem,DIMS(cmasksem), &
       coarseLS,DIMS(coarseLS), &
       velbc, &
       loc,hic) 

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probcommon_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: enable_spectral
      REAL_T, intent(in) :: dxc(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: operation_flag
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T, intent(in) :: dir
      REAL_T, intent(in) :: problo(SDIM)
      INTEGER_T, intent(in) :: level_c
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: bfact_c
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM)
      INTEGER_T, intent(in) :: ncomp_flux
      INTEGER_T, intent(in) :: ncomp_den
      INTEGER_T, intent(in) :: ncomp_vel
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: DIMDEC(den_crse)
      INTEGER_T, intent(in) :: DIMDEC(vel_crse)
      INTEGER_T, intent(in) :: DIMDEC(masknbr)
      INTEGER_T, intent(in) :: DIMDEC(masksem)
      INTEGER_T, intent(in) :: DIMDEC(cmasksem)
      INTEGER_T, intent(in) :: DIMDEC(coarseLS)
      INTEGER_T, intent(in) :: velbc(SDIM,2,SDIM)
      INTEGER_T, intent(in) :: loc(SDIM),hic(SDIM) ! coarse grid dimensions
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T :: stenlo(3),stenhi(3)
      INTEGER_T :: mstenlo(3),mstenhi(3)
      REAL_T, intent(out) :: fine(DIMV(fine),ncomp_flux)
      REAL_T, intent(in) :: den_crse(DIMV(den_crse),ncomp_den)
      REAL_T, intent(in) :: vel_crse(DIMV(vel_crse),ncomp_vel)
       ! =1 if fine-fine  =0 coarse-fine
      REAL_T, intent(in) :: masknbr(DIMV(masknbr))  
      REAL_T, intent(in) :: masksem(DIMV(masksem))
      REAL_T, intent(in) :: cmasksem(DIMV(cmasksem))
      REAL_T, intent(in) :: coarseLS(DIMV(coarseLS),num_materials)

      INTEGER_T clochi(SDIM)
      INTEGER_T ifine,jfine,kfine
      INTEGER_T icoarse,jcoarse,kcoarse
      INTEGER_T ilocal,jlocal,klocal
      INTEGER_T istrip,jstrip,kstrip
      INTEGER_T iface,jface,kface
      INTEGER_T ii,jj,kk
      INTEGER_T n
      INTEGER_T dir2
      INTEGER_T side
      INTEGER_T gridtype
      REAL_T voltotal
      REAL_T volall
      REAL_T wt(SDIM)
      REAL_T xsten(-1:1,SDIM)
      REAL_T xsten_coarse(-1:1,SDIM)
      REAL_T xstenND(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T fine_value(ncomp_flux)
      REAL_T dxelem_c
      REAL_T xfine(SDIM)
      INTEGER_T coarselo_index
      REAL_T, dimension(D_DECL(:,:,:),:),allocatable :: ccrse
      REAL_T INTERP_TOL
      INTEGER_T testmask,testmask2
      INTEGER_T nmat
      INTEGER_T elem_test
      INTEGER_T im,imcrit
      INTEGER_T mat_freq(num_materials)
      REAL_T coarseLS_local(num_materials)
      REAL_T crse_data
      INTEGER_T dencomp
      INTEGER_T local_masknbr
      INTEGER_T caller_id

      caller_id=7
      nhalf=1

      INTERP_TOL=1.0E-4

      nmat=num_materials

      if (operation_flag.eq.8) then ! viscosity
       if (ncomp_vel.ne.SDIM) then
        print *,"ncomp_vel invalid viscosity"
        stop
       endif
       if (ncomp_den.ne.SDIM) then
        print *,"ncomp_den invalid viscosity"
        stop
       endif
       if (ncomp_flux.ne.SDIM) then
        print *,"ncomp_flux invalid viscosity"
        stop
       endif
      else if (operation_flag.eq.7) then ! advection
       if (ncomp_vel.ne.SDIM) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.nmat*num_state_material) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.SDIM+num_state_base) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.1) then ! P cell to P MAC
       if (ncomp_vel.ne.1) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.1) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if ((operation_flag.eq.3).or. & !u cell to MAC
               (operation_flag.eq.5).or. & !UMAC=UMAC+beta diff_reg
               (operation_flag.eq.10).or. &
               (operation_flag.eq.11)) then 
       if (ncomp_vel.ne.AMREX_SPACEDIM) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.AMREX_SPACEDIM) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.9) then
       if (ncomp_vel.ne.nmat*num_state_material) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.nmat*num_state_material) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.0) then
       if (ncomp_vel.ne.1) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.1) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else
       print *,"operation_flag invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid4 ",bfact
       stop
      endif
      if (bfact_c.lt.1) then
       print *,"bfact_c invalid4 ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact).and. &
          (bfact_c.ne.2*bfact)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level-1)) then
       print *,"level_c or level invalid"
       stop
      endif
      if (level.gt.finest_level) then
       print *,"level invalid 40"
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
       print *,"dir invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(fine),0,dir,411)
      call checkbound(loc,hic,DIMS(den_crse),1,-1,411)
      call checkbound(loc,hic,DIMS(vel_crse),1,-1,411)
      call checkbound(fablo,fabhi,DIMS(masknbr),1,-1,411)
      call checkbound(fablo,fabhi,DIMS(masksem),1,-1,411)
      call checkbound(loc,hic,DIMS(cmasksem),1,-1,411)
      call checkbound(loc,hic,DIMS(coarseLS),1,-1,411)

      gridtype=0  ! ggg  (Gauss in all directions)

      do dir2=1,SDIM
       clochi(dir2)=bfact_c-1
       if (2*loc(dir2).ne.fablo(dir2)) then
        print *,"2*loc(dir2).ne.fablo(dir2)"
        stop
       endif
       if (2*(hic(dir2)+1).ne.fabhi(dir2)+1) then
        print *,"2*(hic(dir2)+1).ne.fabhi(dir2)+1"
        stop
       endif
      enddo ! dir2=1..sdim

      allocate(ccrse(D_DECL(0:clochi(1),0:clochi(2),0:clochi(SDIM)), &
                     ncomp_flux))

      do dir2=1,SDIM
       elem_test=loc(dir2)/bfact_c
       if (elem_test*bfact_c.ne.loc(dir2)) then
         print *,"elem_test invalid (loc) 3"
         print *,"elem_test=",elem_test
         print *,"bfact_c=",bfact_c
         print *,"dir2=",dir2
         print *,"loc(dir2)=",loc(dir2)
         stop
       endif
       elem_test=(hic(dir2)+1)/bfact_c
       if (elem_test*bfact_c.ne.hic(dir2)+1) then
         print *,"elem_test invalid (hic) 4"
         print *,"elem_test=",elem_test
         print *,"bfact_c=",bfact_c
         print *,"dir2=",dir2
         print *,"hic(dir2)=",hic(dir2)
         stop
       endif
       elem_test=fablo(dir2)/bfact
       if (elem_test*bfact.ne.fablo(dir2)) then
        print *,"elem_test invalid (fablo) 5"
        print *,"elem_test=",elem_test
        print *,"bfact=",bfact
        print *,"dir2=",dir2
        print *,"fablo(dir2)=",fablo(dir2)
        stop
       endif
       elem_test=(fabhi(dir2)+1)/bfact
       if (elem_test*bfact.ne.fabhi(dir2)+1) then
        print *,"elem_test invalid (fabhi) 6"
        print *,"elem_test=",elem_test
        print *,"bfact=",bfact
        print *,"dir2=",dir2
        print *,"fabhi(dir2)=",fabhi(dir2)
        stop
       endif
      enddo ! dir2=1..sdim

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do side=1,2

        ! strip of cells just outside the fine grid active region
       if (side.eq.1) then
        growlo(dir+1)=fablo(dir+1)-1
        growhi(dir+1)=fablo(dir+1)-1
       else if (side.eq.2) then
        growlo(dir+1)=fabhi(dir+1)+1
        growhi(dir+1)=fabhi(dir+1)+1
       else
        print *,"side invalid"
        stop
       endif

       if (velbc(dir+1,side,dir+1).eq.INT_DIR) then
 
        do ifine=growlo(1),growhi(1)
        do jfine=growlo(2),growhi(2)
        do kfine=growlo(3),growhi(3)

         local_masknbr=NINT(masknbr(D_DECL(ifine,jfine,kfine)))

          ! coarse-fine (in contrast to fine-fine)
         if (local_masknbr.eq.0) then
          call gridsten_level(xsten,ifine,jfine,kfine,level,nhalf)

          ! index for a coarse cell that is next to the fine grid boundary.
          if (side.eq.1) then
           icoarse=(ifine-1)/2
           jcoarse=(jfine-1)/2
           kcoarse=(kfine-1)/2
           iface=ifine+ii
           jface=jfine+jj
           kface=kfine+kk
          else if (side.eq.2) then
           icoarse=ifine/2
           jcoarse=jfine/2
           kcoarse=kfine/2
           iface=ifine
           jface=jfine
           kface=kfine
          else
           print *,"side invalid"
           stop
          endif
          call gridsten_level(xsten_coarse,icoarse,jcoarse,kcoarse, &
            level_c,nhalf)
          xsten(0,dir+1)=xsten_coarse(0,dir+1)

          ! find stencil for surrounding coarse level spectral element.
          stenlo(3)=0
          stenhi(3)=0
          do dir2=1,SDIM
           dxelem_c=dxc(dir2)*bfact_c
           xfine(dir2)=xsten(0,dir2)
           coarselo_index=NINT( (xfine(dir2)-problo(dir2))/dxelem_c-half )
           stenlo(dir2)=bfact_c*coarselo_index
           stenhi(dir2)=stenlo(dir2)+bfact_c-1
           if (dir2.eq.dir+1) then
            if (side.eq.1) then
             if (stenhi(dir2).ne.loc(dir2)-1) then
              print *,"stenhi invalid"
              stop
             endif
             if ((stenlo(dir2).lt.loc(dir2)-bfact_c).or. &
                 (stenhi(dir2).ge.loc(dir2))) then
              print *,"stenlo or stenhi invalid 3"
              stop
             endif
            else if (side.eq.2) then
             if (stenlo(dir2).ne.hic(dir2)+1) then
              print *,"stenlo invalid"
              stop
             endif
             if ((stenlo(dir2).le.hic(dir2)).or. &
                 (stenhi(dir2).gt.hic(dir2)+bfact_c)) then
              print *,"stenlo or stenhi invalid 4"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           else if (dir2.ne.dir+1) then
            if ((stenlo(dir2).lt.loc(dir2)).or. &
                (stenhi(dir2).gt.hic(dir2))) then
             print *,"stenlo or stenhi invalid 5"
             stop
            endif
           else
            print *,"dir2 invalid"
            stop
           endif
          enddo ! dir2=1..sdim

           ! mstenlo,mstenhi are coarse level indexes.
          mstenlo(3)=0
          mstenhi(3)=0
          do dir2=1,SDIM
           mstenlo(dir2)=stenlo(dir2)
           mstenhi(dir2)=stenhi(dir2)
           if (dir2.eq.dir+1) then
            if (side.eq.1) then
             mstenlo(dir2)=mstenhi(dir2)
             if (mstenlo(dir2).ne.loc(dir2)-1) then
              print *,"mstenlo invalid"
              stop
             endif
            else if (side.eq.2) then
             mstenhi(dir2)=mstenlo(dir2)
             if (mstenhi(dir2).ne.hic(dir2)+1) then
              print *,"mstenhi invalid"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           else if (dir2.ne.dir+1) then
            if ((mstenlo(dir2).lt.loc(dir2)).or. &
                (mstenhi(dir2).gt.hic(dir2))) then
             print *,"mstenlo or mstenhi invalid"
             stop
            endif
           else
            print *,"dir2 bust"
            stop
           endif
          enddo ! dir2=1..sdim

          icoarse=stenlo(1) 
          jcoarse=stenlo(2) 
          kcoarse=stenlo(SDIM) 

          ! lower left hand corner of element stencil.
          call gridstenND_level(xstenND,icoarse,jcoarse,kcoarse,level_c,nhalf)

          do dir2=1,SDIM

           xfine(dir2)=xfine(dir2)-xstenND(0,dir2)

           if ((xfine(dir2).lt.-INTERP_TOL*dxc(dir2)).or. &
               (xfine(dir2).gt.(INTERP_TOL+bfact_c)*dxc(dir2))) then
            print *,"fine out of bounds interp_copy"
            stop
           endif

          enddo ! dir2=1..sdim

          do im=1,nmat
           mat_freq(im)=0
          enddo
           ! mstenlo,mstenhi is a thin coarse level stencil
          do ilocal=mstenlo(1),mstenhi(1) 
          do jlocal=mstenlo(2),mstenhi(2) 
          do klocal=mstenlo(3),mstenhi(3) 
           do im=1,nmat
            coarseLS_local(im)=coarseLS(D_DECL(ilocal,jlocal,klocal),im)
           enddo
           call get_primary_material(coarseLS_local,nmat,imcrit)
           if ((imcrit.ge.1).and.(imcrit.le.nmat)) then
            mat_freq(imcrit)=mat_freq(imcrit)+1
           else
            print *,"imcrit invalid"
            stop
           endif
          enddo  !klocal
          enddo  !jlocal
          enddo  !ilocal

          imcrit=0
          do im=1,nmat
           if (mat_freq(im).gt.0) then
            if (imcrit.eq.0) then
             imcrit=im
            else if ((imcrit.ge.1).and.(imcrit.le.nmat)) then
             if (mat_freq(im).gt.mat_freq(imcrit)) then
              imcrit=im
             endif
            else
             print *,"imcrit invalid"
             stop
            endif
           else if (mat_freq(im).eq.0) then
            ! do nothing
           else
            print *,"mat_freq(im) invalid"
            stop
           endif
          enddo ! im=1..nmat

          ilocal=mstenlo(1) 
          jlocal=mstenlo(2) 
          klocal=mstenlo(SDIM) 
          testmask=NINT(cmasksem(D_DECL(ilocal,jlocal,klocal)))
          do ilocal=mstenlo(1),mstenhi(1) 
          do jlocal=mstenlo(2),mstenhi(2) 
          do klocal=mstenlo(3),mstenhi(3) 
           testmask2=NINT(cmasksem(D_DECL(ilocal,jlocal,klocal)))
           if (testmask2.ne.testmask) then
            testmask=0
           endif
          enddo 
          enddo 
          enddo 

          if ((imcrit.ge.1).and.(imcrit.le.nmat)) then

           do n=1,ncomp_flux
            fine_value(n)=zero
           enddo
           voltotal=zero

           if ((bfact_c.ge.2).and. &
               ((enable_spectral.eq.1).or. &
                (enable_spectral.eq.2)).and. &
               (testmask.ge.1).and. &
               (testmask.le.nmat)) then

            do icoarse=stenlo(1),stenhi(1)
            do jcoarse=stenlo(2),stenhi(2)
            do kcoarse=stenlo(3),stenhi(3)
             ilocal=icoarse-stenlo(1)
             jlocal=jcoarse-stenlo(2)
             klocal=kcoarse-stenlo(3)

             istrip=icoarse
             jstrip=jcoarse
             kstrip=kcoarse
             if (dir.eq.0) then
              if (side.eq.1) then
               istrip=stenhi(dir+1)
              else if (side.eq.2) then
               istrip=stenlo(dir+1)
              else
               print *,"side invalid"
               stop
              endif
             else if (dir.eq.1) then
              if (side.eq.1) then
               jstrip=stenhi(dir+1)
              else if (side.eq.2) then
               jstrip=stenlo(dir+1)
              else
               print *,"side invalid"
               stop
              endif
             else if ((dir.eq.2).and.(SDIM.eq.3)) then
              if (side.eq.1) then
               kstrip=stenhi(dir+1)
              else if (side.eq.2) then
               kstrip=stenlo(dir+1)
              else
               print *,"side invalid"
               stop
              endif
             else
              print *,"dir invalid"
              stop
             endif

             do n=1,ncomp_flux

              if (operation_flag.eq.8) then ! viscosity
               if ((n.ge.1).and.(n.le.SDIM)) then
                crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),n)
               else
                print *,"n invalid viscosity"
                stop
               endif
              else if (operation_flag.eq.7) then ! advection
               if ((n.ge.1).and.(n.le.SDIM)) then
                crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),n)
               else if ((n.gt.SDIM).and.(n.le.SDIM+num_state_base)) then
                dencomp=(imcrit-1)*num_state_material+n-SDIM
                crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),dencomp)
               else
                print *,"n invalid"
                stop
               endif
              else if (operation_flag.eq.1) then ! p cell to MAC
               if (n.eq.1) then
                crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),n)
               else
                print *,"n invalid"
                stop
               endif
              else if ((operation_flag.eq.3).or. & !u cell to MAC
                       (operation_flag.eq.5).or. & !UMAC=UMAC+beta diff_reg
                       (operation_flag.eq.10).or. &
                       (operation_flag.eq.11)) then 
               if ((dir.ge.0).and.(dir.lt.AMREX_SPACEDIM)) then
                crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),dir+1)
               else
                print *,"dir invalid"
                stop
               endif
              else if (operation_flag.eq.9) then
               if (n.eq.1) then
                dencomp=(imcrit-1)*num_state_material+1
                crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),dencomp)
               else
                print *,"n invalid"
                stop
               endif
              else if (operation_flag.eq.0) then
               if (n.eq.1) then
                crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),1)
               else
                print *,"n invalid"
                stop
               endif
              else
               print *,"operation_flag invalid"
               stop
              endif

              ccrse(D_DECL(ilocal,jlocal,klocal),n)=crse_data
             enddo ! n=1,ncomp_flux
            enddo ! kcoarse
            enddo ! jcoarse
            enddo ! icoarse

            call SEM_INTERP_ELEMENT( &
             ncomp_flux,bfact_c,gridtype, &
             clochi,dxc,xfine,ccrse,fine_value,caller_id)

            voltotal=one

           else if ((bfact_c.eq.1).or. &
                    (enable_spectral.eq.0).or. &
                    (enable_spectral.eq.3).or. &
                    (testmask.eq.0)) then

            if (side.eq.1) then
             stenlo(dir+1)=loc(dir+1)-1
             stenhi(dir+1)=loc(dir+1)-1
            else if (side.eq.2) then
             stenlo(dir+1)=hic(dir+1)+1
             stenhi(dir+1)=hic(dir+1)+1
            else
             print *,"side invalid"
             stop
            endif

            do icoarse=stenlo(1),stenhi(1)
             call intersect_weight_interp_COPY(icoarse,ifine,bfact_c,bfact, &
                    wt(1),dir,0)
             if (wt(1).gt.zero) then
              do jcoarse=stenlo(2),stenhi(2)
               call intersect_weight_interp_COPY(jcoarse,jfine,bfact_c,bfact, &
                      wt(2),dir,1)
               if (wt(2).gt.zero) then
                do kcoarse=stenlo(3),stenhi(3)
                 if (SDIM.eq.3) then
                  call intersect_weight_interp_COPY(kcoarse,kfine,bfact_c, &
                         bfact,wt(SDIM),dir,2)
                 endif
                 if (wt(SDIM).gt.zero) then
                  volall=wt(1)
                  do dir2=2,SDIM
                   volall=volall*wt(dir2)
                  enddo

                  istrip=icoarse
                  jstrip=jcoarse
                  kstrip=kcoarse
                  if (dir.eq.0) then
                   if (side.eq.1) then
                    istrip=stenhi(dir+1)
                   else if (side.eq.2) then
                    istrip=stenlo(dir+1)
                   else
                    print *,"side invalid"
                    stop
                   endif
                  else if (dir.eq.1) then
                   if (side.eq.1) then
                    jstrip=stenhi(dir+1)
                   else if (side.eq.2) then
                    jstrip=stenlo(dir+1)
                   else
                    print *,"side invalid"
                    stop
                   endif
                  else if ((dir.eq.2).and.(SDIM.eq.3)) then
                   if (side.eq.1) then
                    kstrip=stenhi(dir+1)
                   else if (side.eq.2) then
                    kstrip=stenlo(dir+1)
                   else
                    print *,"side invalid"
                    stop
                   endif
                  else
                   print *,"dir invalid"
                   stop
                  endif

                  do n = 1, ncomp_flux

                   if (operation_flag.eq.8) then ! viscosity
                    if ((n.ge.1).and.(n.le.SDIM)) then
                     crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),n)
                    else
                     print *,"n invalid"
                     stop
                    endif
                   else if (operation_flag.eq.7) then ! advection
                    if ((n.ge.1).and.(n.le.SDIM)) then
                     crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),n)
                    else if ((n.gt.SDIM).and.(n.le.SDIM+num_state_base)) then
                     dencomp=(imcrit-1)*num_state_material+n-SDIM
                     crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),dencomp)
                    else
                     print *,"n invalid"
                     stop
                    endif
                   else if (operation_flag.eq.1) then ! p cell to MAC
                    if (n.eq.1) then
                     crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),n)
                    else
                     print *,"n invalid"
                     stop
                    endif
                   else if ((operation_flag.eq.3).or. &!u cell to MAC
                            (operation_flag.eq.5).or. &!UMAC=UMAC+beta diff_reg
                            (operation_flag.eq.10).or. &
                            (operation_flag.eq.11)) then 
                    if ((dir.ge.0).and.(dir.lt.AMREX_SPACEDIM)) then
                     crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),dir+1)
                    else
                     print *,"dir invalid"
                     stop
                    endif
                   else if (operation_flag.eq.9) then
                    if (n.eq.1) then
                     dencomp=(imcrit-1)*num_state_material+1
                     crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),dencomp)
                    else
                     print *,"n invalid"
                     stop
                    endif
                   else if (operation_flag.eq.0) then
                    if (n.eq.1) then
                     crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),1)
                    else
                     print *,"n invalid"
                     stop
                    endif
                   else
                    print *,"operation_flag invalid"
                    stop
                   endif

                   fine_value(n) = fine_value(n) + volall*crse_data
                  enddo ! n=1..ncomp_flux

                  voltotal=voltotal+volall
                 endif ! wt(sdim).gt.0
                enddo ! kcoarse
               endif
              enddo ! jcoarse
             endif
            enddo ! icoarse

           else
            print *,"parameter bust in interp_copy"
            stop
           endif

           if (voltotal.gt.zero) then
            do n = 1,ncomp_flux
             fine(D_DECL(iface,jface,kface),n)=fine_value(n)/voltotal
            enddo
           else
            print *,"voltotal invalid"
            stop
           endif

          else if (imcrit.eq.0) then
           ! do nothing
          else
           print *,"imcrit invalid"
           stop
          endif
         else if (local_masknbr.eq.1) then ! fine-fine
          ! do nothing
         else
          print *,"local_masknbr invalid"
          stop
         endif
 
        enddo !kfine
        enddo !jfine
        enddo !ifine

       else if ((velbc(dir+1,side,dir+1).eq.EXT_DIR).or. &
                (velbc(dir+1,side,dir+1).eq.REFLECT_EVEN).or. &
                (velbc(dir+1,side,dir+1).eq.REFLECT_ODD).or. &
                (velbc(dir+1,side,dir+1).eq.FOEXTRAP)) then
        ! do nothing
       else
        print *,"velbc invalid"
        stop
       endif

      enddo ! side=1..2

      deallocate(ccrse)

      return
      end subroutine FORT_INTERP_COPY

      subroutine FORT_INTERP_FLUX( &
       enable_spectral, &
       dxc,dx, &
       finest_level, &
       tilelo,tilehi, &
       fablo,fabhi, &
       dir, &
       problo, &
       level_c,level, &
       bfact_c,bfact, &
       xlo, &
       ncomp_flux, &
       fine,DIMS(fine), &
       crse,DIMS(crse), &
       masknbr,DIMS(masknbr), &
       masksem,DIMS(masksem), &
       cmasksem,DIMS(cmasksem), &
       velbc, &
       loc,hic) 

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probcommon_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T enable_spectral
      REAL_T dxc(SDIM)
      REAL_T dx(SDIM)
      INTEGER_T finest_level
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T dir
      REAL_T problo(SDIM)
      INTEGER_T level_c
      INTEGER_T level
      INTEGER_T bfact_c
      INTEGER_T bfact
      REAL_T xlo(SDIM)
      INTEGER_T ncomp_flux
      INTEGER_T DIMDEC(fine)
      INTEGER_T DIMDEC(crse)
      INTEGER_T DIMDEC(masknbr)
      INTEGER_T DIMDEC(masksem)
      INTEGER_T DIMDEC(cmasksem)
      INTEGER_T velbc(SDIM,2,SDIM)
      INTEGER_T loc(SDIM),hic(SDIM) ! coarse grid dimensions
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      INTEGER_T mstenlo(3),mstenhi(3)
      REAL_T   fine(DIMV(fine),ncomp_flux)
      REAL_T   crse(DIMV(crse),ncomp_flux)
      REAL_T   masknbr(DIMV(masknbr))  ! =1 if fine-fine  =0 coarse-fine
      REAL_T   masksem(DIMV(masksem))
      REAL_T   cmasksem(DIMV(cmasksem))
      INTEGER_T clochi(SDIM)
      INTEGER_T ifine,jfine,kfine
      INTEGER_T icoarse,jcoarse,kcoarse
      INTEGER_T ilocal,jlocal,klocal
      INTEGER_T istrip,jstrip,kstrip
      INTEGER_T iface,jface,kface
      INTEGER_T ii,jj,kk
      INTEGER_T n
      INTEGER_T dir2
      INTEGER_T side
      INTEGER_T gridtype
      REAL_T voltotal
      REAL_T volall
      REAL_T wt(SDIM)
      REAL_T xsten(-1:1,SDIM)
      REAL_T xsten_coarse(-1:1,SDIM)
      REAL_T xstenND(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T fine_value(ncomp_flux)
      REAL_T dxelem_c
      REAL_T xfine(SDIM)
      INTEGER_T coarselo_index
      REAL_T, dimension(D_DECL(:,:,:),:),allocatable :: ccrse
      REAL_T INTERP_TOL
      INTEGER_T testmask,testmask2
      INTEGER_T nmat
      INTEGER_T elem_test
      REAL_T crse_data
      INTEGER_T local_masknbr
      INTEGER_T caller_id

      caller_id=17
      nhalf=1

      INTERP_TOL=1.0E-4

      nmat=num_materials

      if (ncomp_flux.ne.SDIM+num_state_base) then
       print *,"ncomp_flux invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid4 ",bfact
       stop
      endif
      if (bfact_c.lt.1) then
       print *,"bfact_c invalid4 ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact).and. &
          (bfact_c.ne.2*bfact)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level-1)) then
       print *,"level_c or level invalid"
       stop
      endif
      if (level.gt.finest_level) then
       print *,"level invalid 41"
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
       print *,"dir invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(fine),0,dir,411)
      call checkbound(loc,hic,DIMS(crse),0,dir,411)
      call checkbound(fablo,fabhi,DIMS(masknbr),1,-1,411)
      call checkbound(fablo,fabhi,DIMS(masksem),1,-1,411)
      call checkbound(loc,hic,DIMS(cmasksem),1,-1,411)

      gridtype=0  ! ggg  (Gauss in all directions)

      do dir2=1,SDIM
       clochi(dir2)=bfact_c-1
       if (2*loc(dir2).ne.fablo(dir2)) then
        print *,"2*loc(dir2).ne.fablo(dir2)"
        stop
       endif
       if (2*(hic(dir2)+1).ne.fabhi(dir2)+1) then
        print *,"2*(hic(dir2)+1).ne.fabhi(dir2)+1"
        stop
       endif
      enddo ! dir2=1..sdim

      allocate(ccrse(D_DECL(0:clochi(1),0:clochi(2),0:clochi(SDIM)), &
                     ncomp_flux))

      do dir2=1,SDIM
       elem_test=loc(dir2)/bfact_c
       if (elem_test*bfact_c.ne.loc(dir2)) then
         print *,"elem_test invalid (loc) 7"
         print *,"elem_test=",elem_test
         print *,"bfact_c=",bfact_c
         print *,"dir2=",dir2
         print *,"loc(dir2)=",loc(dir2)
         stop
       endif
       elem_test=(hic(dir2)+1)/bfact_c
       if (elem_test*bfact_c.ne.hic(dir2)+1) then
         print *,"elem_test invalid (hic) 8"
         print *,"elem_test=",elem_test
         print *,"bfact_c=",bfact_c
         print *,"dir2=",dir2
         print *,"hic(dir2)=",hic(dir2)
         stop
       endif
       elem_test=fablo(dir2)/bfact
       if (elem_test*bfact.ne.fablo(dir2)) then
        print *,"elem_test invalid (fablo) 10"
        print *,"elem_test=",elem_test
        print *,"bfact=",bfact
        print *,"dir2=",dir2
        print *,"fablo(dir2)=",fablo(dir2)
        stop
       endif
       elem_test=(fabhi(dir2)+1)/bfact
       if (elem_test*bfact.ne.fabhi(dir2)+1) then
        print *,"elem_test invalid (fabhi) 11"
        print *,"elem_test=",elem_test
        print *,"bfact=",bfact
        print *,"dir2=",dir2
        print *,"fabhi(dir2)=",fabhi(dir2)
        stop
       endif
      enddo ! dir2=1..sdim

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

      do side=1,2

        ! strip of cells just outside the fine grid active region
       if (side.eq.1) then
        growlo(dir+1)=fablo(dir+1)-1
        growhi(dir+1)=fablo(dir+1)-1
       else if (side.eq.2) then
        growlo(dir+1)=fabhi(dir+1)+1
        growhi(dir+1)=fabhi(dir+1)+1
       else
        print *,"side invalid"
        stop
       endif

       if (velbc(dir+1,side,dir+1).eq.INT_DIR) then
 
        do ifine=growlo(1),growhi(1)
        do jfine=growlo(2),growhi(2)
        do kfine=growlo(3),growhi(3)

         local_masknbr=NINT(masknbr(D_DECL(ifine,jfine,kfine)))

          ! coarse-fine (in contrast to fine-fine)
         if (local_masknbr.eq.0) then
          call gridsten_level(xsten,ifine,jfine,kfine,level,nhalf)

          ! index for a coarse cell that is next to the fine grid boundary.
          if (side.eq.1) then
           icoarse=(ifine-1)/2
           jcoarse=(jfine-1)/2
           kcoarse=(kfine-1)/2
           iface=ifine+ii
           jface=jfine+jj
           kface=kfine+kk
          else if (side.eq.2) then
           icoarse=ifine/2
           jcoarse=jfine/2
           kcoarse=kfine/2
           iface=ifine
           jface=jfine
           kface=kfine
          else
           print *,"side invalid"
           stop
          endif
          call gridsten_level(xsten_coarse,icoarse,jcoarse,kcoarse, &
            level_c,nhalf)
          xsten(0,dir+1)=xsten_coarse(0,dir+1)

          ! find stencil for surrounding coarse level spectral element.
          stenlo(3)=0
          stenhi(3)=0
          do dir2=1,SDIM
           dxelem_c=dxc(dir2)*bfact_c
           xfine(dir2)=xsten(0,dir2)
           coarselo_index=NINT( (xfine(dir2)-problo(dir2))/dxelem_c-half )
           stenlo(dir2)=bfact_c*coarselo_index
           stenhi(dir2)=stenlo(dir2)+bfact_c-1
           if (dir2.eq.dir+1) then
            if (side.eq.1) then
             if (stenhi(dir2).ne.loc(dir2)-1) then
              print *,"stenhi invalid"
              stop
             endif
             if ((stenlo(dir2).lt.loc(dir2)-bfact_c).or. &
                 (stenhi(dir2).ge.loc(dir2))) then
              print *,"stenlo or stenhi invalid 6"
              stop
             endif
            else if (side.eq.2) then
             if (stenlo(dir2).ne.hic(dir2)+1) then
              print *,"stenlo invalid"
              stop
             endif
             if ((stenlo(dir2).le.hic(dir2)).or. &
                 (stenhi(dir2).gt.hic(dir2)+bfact_c)) then
              print *,"stenlo or stenhi invalid 7"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           else if (dir2.ne.dir+1) then
            if ((stenlo(dir2).lt.loc(dir2)).or. &
                (stenhi(dir2).gt.hic(dir2))) then
             print *,"stenlo or stenhi invalid 8"
             stop
            endif
           else
            print *,"dir2 invalid"
            stop
           endif
          enddo ! dir2=1..sdim

          mstenlo(3)=0
          mstenhi(3)=0
          do dir2=1,SDIM
           mstenlo(dir2)=stenlo(dir2)
           mstenhi(dir2)=stenhi(dir2)
           if (dir2.eq.dir+1) then
            if (side.eq.1) then
             mstenlo(dir2)=mstenhi(dir2)
             if (mstenlo(dir2).ne.loc(dir2)-1) then
              print *,"mstenlo invalid"
              stop
             endif
            else if (side.eq.2) then
             mstenhi(dir2)=mstenlo(dir2)
             if (mstenhi(dir2).ne.hic(dir2)+1) then
              print *,"mstenhi invalid"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           else if (dir2.ne.dir+1) then
            if ((mstenlo(dir2).lt.loc(dir2)).or. &
                (mstenhi(dir2).gt.hic(dir2))) then
             print *,"mstenlo or mstenhi invalid"
             stop
            endif
           else
            print *,"dir2 bust"
            stop
           endif
          enddo ! dir2=1..sdim

          icoarse=stenlo(1) 
          jcoarse=stenlo(2) 
          kcoarse=stenlo(SDIM) 

          ! lower left hand corner of element stencil.
          call gridstenND_level(xstenND,icoarse,jcoarse,kcoarse,level_c,nhalf)

          do dir2=1,SDIM

           xfine(dir2)=xfine(dir2)-xstenND(0,dir2)

           if ((xfine(dir2).lt.-INTERP_TOL*dxc(dir2)).or. &
               (xfine(dir2).gt.(INTERP_TOL+bfact_c)*dxc(dir2))) then
            print *,"fine out of bounds interp_copy"
            stop
           endif

          enddo ! dir2=1..sdim

          ilocal=mstenlo(1) 
          jlocal=mstenlo(2) 
          klocal=mstenlo(SDIM) 
          testmask=NINT(cmasksem(D_DECL(ilocal,jlocal,klocal)))
          do ilocal=mstenlo(1),mstenhi(1) 
          do jlocal=mstenlo(2),mstenhi(2) 
          do klocal=mstenlo(3),mstenhi(3) 
           testmask2=NINT(cmasksem(D_DECL(ilocal,jlocal,klocal)))
           if (testmask2.ne.testmask) then
            testmask=0
           endif
          enddo 
          enddo 
          enddo 

          do n=1,ncomp_flux
           fine_value(n)=zero
          enddo
          voltotal=zero

          if ((bfact_c.ge.2).and. &
              ((enable_spectral.eq.1).or. &
               (enable_spectral.eq.2)).and. &
              (testmask.ge.1).and. &
              (testmask.le.nmat)) then

            do icoarse=stenlo(1),stenhi(1)
            do jcoarse=stenlo(2),stenhi(2)
            do kcoarse=stenlo(3),stenhi(3)
             ilocal=icoarse-stenlo(1)
             jlocal=jcoarse-stenlo(2)
             klocal=kcoarse-stenlo(3)

             istrip=icoarse
             jstrip=jcoarse
             kstrip=kcoarse
             if (dir.eq.0) then
              if (side.eq.1) then
               istrip=stenhi(dir+1)+1
              else if (side.eq.2) then
               istrip=stenlo(dir+1)
              else
               print *,"side invalid"
               stop
              endif
             else if (dir.eq.1) then
              if (side.eq.1) then
               jstrip=stenhi(dir+1)+1
              else if (side.eq.2) then
               jstrip=stenlo(dir+1)
              else
               print *,"side invalid"
               stop
              endif
             else if ((dir.eq.2).and.(SDIM.eq.3)) then
              if (side.eq.1) then
               kstrip=stenhi(dir+1)+1
              else if (side.eq.2) then
               kstrip=stenlo(dir+1)
              else
               print *,"side invalid"
               stop
              endif
             else
              print *,"dir invalid"
              stop
             endif

             do n=1,ncomp_flux
              crse_data=crse(D_DECL(istrip,jstrip,kstrip),n)
              ccrse(D_DECL(ilocal,jlocal,klocal),n)=crse_data
             enddo ! n=1,ncomp_flux
            enddo ! kcoarse
            enddo ! jcoarse
            enddo ! icoarse

            call SEM_INTERP_ELEMENT( &
             ncomp_flux,bfact_c,gridtype, &
             clochi,dxc,xfine,ccrse,fine_value,caller_id)

            voltotal=one

          else if ((bfact_c.eq.1).or. &
                   (enable_spectral.eq.0).or. &
                   (enable_spectral.eq.3).or. &
                   (testmask.eq.0)) then

            if (side.eq.1) then
             stenlo(dir+1)=loc(dir+1)-1
             stenhi(dir+1)=loc(dir+1)-1
            else if (side.eq.2) then
             stenlo(dir+1)=hic(dir+1)+1
             stenhi(dir+1)=hic(dir+1)+1
            else
             print *,"side invalid"
             stop
            endif

            do icoarse=stenlo(1),stenhi(1)
             call intersect_weight_interp_COPY(icoarse,ifine,bfact_c,bfact, &
                    wt(1),dir,0)
             if (wt(1).gt.zero) then
              do jcoarse=stenlo(2),stenhi(2)
               call intersect_weight_interp_COPY(jcoarse,jfine,bfact_c,bfact, &
                      wt(2),dir,1)
               if (wt(2).gt.zero) then
                do kcoarse=stenlo(3),stenhi(3)
                 if (SDIM.eq.3) then
                  call intersect_weight_interp_COPY(kcoarse,kfine,bfact_c, &
                         bfact,wt(SDIM),dir,2)
                 endif
                 if (wt(SDIM).gt.zero) then
                  volall=wt(1)
                  do dir2=2,SDIM
                   volall=volall*wt(dir2)
                  enddo

                  istrip=icoarse
                  jstrip=jcoarse
                  kstrip=kcoarse
                  if (dir.eq.0) then
                   if (side.eq.1) then
                    istrip=stenhi(dir+1)+1
                   else if (side.eq.2) then
                    istrip=stenlo(dir+1)
                   else
                    print *,"side invalid"
                    stop
                   endif
                  else if (dir.eq.1) then
                   if (side.eq.1) then
                    jstrip=stenhi(dir+1)+1
                   else if (side.eq.2) then
                    jstrip=stenlo(dir+1)
                   else
                    print *,"side invalid"
                    stop
                   endif
                  else if ((dir.eq.2).and.(SDIM.eq.3)) then
                   if (side.eq.1) then
                    kstrip=stenhi(dir+1)+1
                   else if (side.eq.2) then
                    kstrip=stenlo(dir+1)
                   else
                    print *,"side invalid"
                    stop
                   endif
                  else
                   print *,"dir invalid"
                   stop
                  endif

                  do n = 1, ncomp_flux
                   crse_data=crse(D_DECL(istrip,jstrip,kstrip),n)
                   fine_value(n) = fine_value(n) + volall*crse_data
                  enddo ! n=1..ncomp_flux

                  voltotal=voltotal+volall
                 endif ! wt(sdim).gt.0
                enddo ! kcoarse
               endif
              enddo ! jcoarse
             endif
            enddo ! icoarse

          else
            print *,"parameter bust in interp_flux"
            stop
          endif

          if (voltotal.gt.zero) then
           do n = 1,ncomp_flux
            fine(D_DECL(iface,jface,kface),n)=fine_value(n)/voltotal
           enddo
          else
           print *,"voltotal invalid"
           stop
          endif

         else if (local_masknbr.eq.1) then ! fine-fine
          ! do nothing
         else
          print *,"local_masknbr invalid"
          stop
         endif
 
        enddo !kfine
        enddo !jfine
        enddo !ifine

       else if ((velbc(dir+1,side,dir+1).eq.EXT_DIR).or. &
                (velbc(dir+1,side,dir+1).eq.REFLECT_EVEN).or. &
                (velbc(dir+1,side,dir+1).eq.REFLECT_ODD).or. &
                (velbc(dir+1,side,dir+1).eq.FOEXTRAP)) then
        ! do nothing
       else
        print *,"velbc invalid"
        stop
       endif

      enddo ! side=1..2

      deallocate(ccrse)

      return
      end subroutine FORT_INTERP_FLUX


      subroutine FORT_FILLBDRY_FLUX( &
       sync_iter, &
       level, &
       finest_level, &
       tilelo,tilehi, &
       fablo,fabhi, &
       dir, &
       ncomp_flux, &
       fluxtarg,DIMS(fluxtarg), &
       fluxhold,DIMS(fluxhold), &
       maskcov,DIMS(maskcov), &
       masknbr,DIMS(masknbr), &
       presbc)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T sync_iter
      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T dir
      INTEGER_T ncomp_flux
      INTEGER_T DIMDEC(fluxtarg)
      INTEGER_T DIMDEC(fluxhold)
      INTEGER_T DIMDEC(maskcov)
      INTEGER_T DIMDEC(masknbr)
      INTEGER_T presbc(SDIM,2)
       ! maskcov=tag if not covered by level+1 or outside the domain.
      REAL_T maskcov(DIMV(maskcov))
      REAL_T masknbr(DIMV(masknbr))  ! =1 if fine-fine  =0 coarse-fine
      REAL_T fluxtarg(DIMV(fluxtarg),ncomp_flux) 
      REAL_T fluxhold(DIMV(fluxhold),ncomp_flux) 
      INTEGER_T i,j,k
      INTEGER_T iface,jface,kface
      INTEGER_T i_out,j_out,k_out
      INTEGER_T ii,jj,kk
      INTEGER_T dirloc,side
      INTEGER_T gridlo(3),gridhi(3)
      INTEGER_T nc
      INTEGER_T pbc
      INTEGER_T maskcov_in
      INTEGER_T maskcov_out
      INTEGER_T masknbr_out
      REAL_T flux_orig,flux_copy


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
       print *,"dir invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(fluxtarg),0,dir,411)
      call checkbound(fablo,fabhi,DIMS(fluxhold),1,-1,411)
      call checkbound(fablo,fabhi,DIMS(masknbr),1,-1,411)
      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,411)

      if ((level.ge.0).and.(level.le.finest_level)) then
       ! do nothing
      else
       print *,"level invalid"
       stop
      endif
      if ((dir.ge.0).and.(dir.lt.SDIM)) then
       ! do nothing
      else
       print *,"dir invalid"
       stop
      endif
      if ((ncomp_flux.ge.1).and.(ncomp_flux.le.10000)) then
       ! do nothing
      else
       print *,"ncomp_flux invalid"
       stop
      endif

      do side=0,1

       gridlo(3)=0
       gridhi(3)=0
       do dirloc=1,SDIM
        gridlo(dirloc)=fablo(dirloc)
        gridhi(dirloc)=fabhi(dirloc)
        if (2*(fablo(dirloc)/2).ne.fablo(dirloc)) then
         print *,"fablo invalid"
         stop
        endif
        if (2*((fabhi(dirloc)+1)/2).ne.fabhi(dirloc)+1) then
         print *,"fabhi invalid"
         stop
        endif
        if (fabhi(dirloc).le.fablo(dirloc)) then
         print *,"fabhi invalid"
         stop
        endif
        if (dirloc.eq.dir+1) then
         if (side.eq.0) then
          gridhi(dirloc)=gridlo(dirloc)
         else if (side.eq.1) then
          gridlo(dirloc)=gridhi(dirloc)
         else
          print *,"side invalid"
          stop
         endif
        else if (dirloc.ne.dir+1) then
         ! do nothing
        else
         print *,"dirloc bust"
         stop
        endif
       enddo ! dirloc=1..sdim

       do i=gridlo(1),gridhi(1)
       do j=gridlo(2),gridhi(2)
       do k=gridlo(3),gridhi(3)

        if (side.eq.0) then
         iface=i
         jface=j
         kface=k
         i_out=i-ii
         j_out=j-jj
         k_out=k-kk
        else if (side.eq.1) then
         iface=i+ii
         jface=j+jj
         kface=k+kk
         i_out=i+ii
         j_out=j+jj
         k_out=k+kk
        else
         print *,"side invalid"
         stop
        endif

        do nc=1,ncomp_flux

         if (sync_iter.eq.0) then
          fluxhold(D_DECL(i,j,k),nc)=fluxtarg(D_DECL(iface,jface,kface),nc)
         else if (sync_iter.eq.1) then
          pbc=presbc(dir+1,side+1)
          if ((pbc.eq.EXT_DIR).or. &
              (pbc.eq.FOEXTRAP).or. &
              (pbc.eq.REFLECT_ODD).or. &
              (pbc.eq.REFLECT_EVEN)) then
           ! do nothing
          else if (pbc.eq.INT_DIR) then
           maskcov_in=NINT(maskcov(D_DECL(i,j,k)))
           maskcov_out=NINT(maskcov(D_DECL(i_out,j_out,k_out)))
           masknbr_out=NINT(masknbr(D_DECL(i_out,j_out,k_out)))
           if (masknbr_out.eq.0) then ! coarse-fine
            ! do nothing
           else if (masknbr_out.eq.1) then ! fine-fine
            if ((maskcov_in.eq.1).and.(maskcov_out.eq.1)) then
             ! do nothing (no fine face above this coarse face)
            else if ((maskcov_in.eq.0).and.(maskcov_out.eq.0)) then
             ! do nothing (not a coarse-fine interface at level+1
            else if ((maskcov_in.eq.1).and.(maskcov_out.eq.0)) then
             flux_orig=fluxtarg(D_DECL(iface,jface,kface),nc)
             flux_copy=fluxhold(D_DECL(i,j,k),nc)
             if (abs(flux_orig-flux_copy).lt.1.0E-10) then
              ! do nothing
             else
              print *,"flux_orig or flux_copy invalid"
              stop
             endif
            else if ((maskcov_in.eq.0).and.(maskcov_out.eq.1)) then
             flux_copy=fluxhold(D_DECL(i_out,j_out,k_out),nc)
             if (abs(flux_copy).lt.1.0D+20) then
              fluxtarg(D_DECL(iface,jface,kface),nc)=flux_copy
             else
              print *,"flux_copy overflow"
              stop
             endif
            else
             print *,"maskcov_in or maskcov_out invalid"
             stop
            endif
           else
            print *,"masknbr_out invalid"
            stop
           endif

          else
           print *,"pbc invalid"
           stop
          endif

         else
          print *,"sync_iter invalid"
          stop
         endif

        enddo ! nc=1..ncomp_flux

       enddo ! k
       enddo ! j
       enddo ! i

      enddo ! side=0..1

      return
      end subroutine FORT_FILLBDRY_FLUX


      subroutine FORT_AVGDOWN_LOW( &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) 

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level_c
      INTEGER_T, intent(in) :: level_f
      INTEGER_T, intent(in) :: bfact_c
      INTEGER_T, intent(in) :: bfact_f
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: dxf(SDIM)
      REAL_T, intent(in) :: xlo_fine(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: DIMDEC(crse)
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      INTEGER_T, intent(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      REAL_T, intent(out) :: crse(DIMV(crse),ncomp)
      REAL_T, intent(in) :: fine(DIMV(fine),ncomp)
      INTEGER_T ic,jc,kc
      INTEGER_T ifine,jfine,kfine
      INTEGER_T n
      INTEGER_T dir2
      REAL_T voltotal
      REAL_T volall
      REAL_T wt(SDIM)
      REAL_T crse_value(ncomp)

      if (ncomp.lt.1) then
       print *,"ncomp invalid32"
       stop
      endif

      if (bfact_f.lt.1) then
       print *,"bfact_f invalid4 ",bfact_f
       stop
      endif
      if (bfact_c.lt.1) then
       print *,"bfact_c invalid4 ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level_f-1)) then
       print *,"level_c or level_f invalid"
       stop
      endif

      call checkbound(lo,hi,DIMS(crse),0,-1,411)
      call checkbound(lof,hif,DIMS(fine),0,-1,411)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)

       do n=1,ncomp
        crse_value(n)=zero
       enddo
       voltotal=zero

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
        bfact_c,bfact_f)

       do ifine=stenlo(1),stenhi(1)
        call intersect_weight_avg(ic,ifine,bfact_c,bfact_f,wt(1))
        if (wt(1).gt.zero) then
         do jfine=stenlo(2),stenhi(2)
          call intersect_weight_avg(jc,jfine,bfact_c,bfact_f,wt(2))
          if (wt(2).gt.zero) then
           do kfine=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             volall=wt(1)
             do dir2=2,SDIM
              volall=volall*wt(dir2)
             enddo
             do n = 1, ncomp
              crse_value(n) = crse_value(n) +  &
               volall*fine(D_DECL(ifine,jfine,kfine),n)
             enddo
             voltotal=voltotal+volall
            endif ! wt(sdim).gt.0
           enddo ! kfine
          endif
         enddo ! jfine
        endif
       enddo ! ifine

       if (voltotal.le.zero) then
        print *,"voltotal invalid"
        stop
       endif

       do n = 1,ncomp
        crse(D_DECL(ic,jc,kc),n)=crse_value(n)/voltotal
       enddo

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine FORT_AVGDOWN_LOW

       ! NavierStokes::level_phase_change_redistribute
       ! NavierStokes::avgDown_tag_localMF
       ! NavierStokes::level_avgDown_tag
      subroutine FORT_AVGDOWN_TAG( &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) 

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level_c
      INTEGER_T, intent(in) :: level_f
      INTEGER_T, intent(in) :: bfact_c
      INTEGER_T, intent(in) :: bfact_f
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: dxf(SDIM)
      REAL_T, intent(in) :: xlo_fine(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: DIMDEC(crse)
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      INTEGER_T, intent(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      REAL_T, intent(out) :: crse(DIMV(crse),ncomp)
      REAL_T, intent(in) :: fine(DIMV(fine),ncomp)
      INTEGER_T ic,jc,kc
      INTEGER_T ifine,jfine,kfine
      INTEGER_T n
      INTEGER_T dir2
      REAL_T voltotal
      REAL_T volall
      REAL_T wt(SDIM)
      REAL_T crse_value(ncomp)
      INTEGER_T fine_test
      INTEGER_T coarse_test
      INTEGER_T init_flag

      if (ncomp.ne.1) then
       print *,"ncomp invalid33"
       stop
      endif

      if (bfact_f.lt.1) then
       print *,"bfact_f invalid4 ",bfact_f
       stop
      endif
      if (bfact_c.lt.1) then
       print *,"bfact_c invalid4 ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level_f-1)) then
       print *,"level_c or level_f invalid"
       stop
      endif

      call checkbound(lo,hi,DIMS(crse),0,-1,411)
      call checkbound(lof,hif,DIMS(fine),0,-1,411)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)

       do n=1,ncomp
        crse_value(n)=zero
       enddo
       init_flag=0
       voltotal=zero

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
        bfact_c,bfact_f)

       do ifine=stenlo(1),stenhi(1)
        call intersect_weight_avg(ic,ifine,bfact_c,bfact_f,wt(1))
        if (wt(1).gt.zero) then
         do jfine=stenlo(2),stenhi(2)
          call intersect_weight_avg(jc,jfine,bfact_c,bfact_f,wt(2))
          if (wt(2).gt.zero) then
           do kfine=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             volall=wt(1)
             do dir2=2,SDIM
              volall=volall*wt(dir2)
             enddo

             do n = 1, ncomp
              fine_test=NINT(fine(D_DECL(ifine,jfine,kfine),n))
              coarse_test=NINT(crse_value(n))
              if (init_flag.eq.0) then
               crse_value(n)=fine_test
              else if (init_flag.eq.1) then
               if (fine_test.eq.0) then
                crse_value(n)=zero
               else if ((fine_test.eq.1).and.(coarse_test.ne.1)) then
                crse_value(n)=zero
               else if ((fine_test.eq.2).and.(coarse_test.ne.2)) then
                crse_value(n)=zero
               else if (fine_test.eq.coarse_test) then
                ! do nothing
               else
                print *,"fine_test bad (avgdown tag):",fine_test 
                stop
               endif
              else
               print *,"init_flag invalid"
               stop
              endif

             enddo ! n

             if (init_flag.eq.0) then
              init_flag=1
             else if (init_flag.eq.1) then
              ! do nothing
             else
              print *,"init_flag invalid"
              stop
             endif

             voltotal=voltotal+volall
            endif ! wt(sdim).gt.0
           enddo ! kfine
          endif
         enddo ! jfine
        endif
       enddo ! ifine

       if (voltotal.le.zero) then
        print *,"voltotal invalid"
        stop
       endif

       do n = 1,ncomp
        crse(D_DECL(ic,jc,kc),n)=crse_value(n)
       enddo

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine FORT_AVGDOWN_TAG


       ! do_the_advance -> level_phase_change_rate
       ! avgDownBURNING_localMF 
       ! level_avgDownBURNING
       ! (note: after level_avgDownBURNING comes 
       !  level_phase_change_rate_extend)
      subroutine FORT_AVGDOWN_BURNING( &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp, &
       nmat, &
       nten, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) 

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level_c
      INTEGER_T, intent(in) :: level_f
      INTEGER_T, intent(in) :: bfact_c
      INTEGER_T, intent(in) :: bfact_f
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: dxf(SDIM)
      REAL_T, intent(in) :: xlo_fine(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nten
      INTEGER_T, intent(in) :: DIMDEC(crse)
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      INTEGER_T, intent(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      REAL_T, intent(out) :: crse(DIMV(crse),ncomp)
      REAL_T, intent(in) :: fine(DIMV(fine),ncomp)
      INTEGER_T ic,jc,kc
      INTEGER_T ifine,jfine,kfine
      INTEGER_T dir2
      INTEGER_T iten
      INTEGER_T n
      REAL_T voltotal
      REAL_T volall
      REAL_T wt(SDIM)
      REAL_T crse_value(ncomp)
      INTEGER_T fine_test
      INTEGER_T coarse_test
      REAL_T velwt(nten)
      INTEGER_T local_comp
      INTEGER_T avgdown_sweep

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif

      if (nten.eq.( (nmat-1)*(nmat-1)+nmat-1)/2) then
       ! do nothing
      else
       print *,"nten invalid"
       stop
      endif

      if (ncomp.eq.nten*(SDIM+1)) then
       ! do nothing
      else
       print *,"ncomp invalid34"
       stop
      endif

      if (bfact_f.lt.1) then
       print *,"bfact_f invalid4 ",bfact_f
       stop
      endif
      if (bfact_c.lt.1) then
       print *,"bfact_c invalid4 ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level_f-1)) then
       print *,"level_c or level_f invalid"
       stop
      endif

      call checkbound(lo,hi,DIMS(crse),0,-1,411)
      call checkbound(lof,hif,DIMS(fine),0,-1,411)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)

       do n=1,ncomp
        crse_value(n)=zero
       enddo

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
        bfact_c,bfact_f)

       do avgdown_sweep=0,1

        do iten=1,nten
         velwt(iten)=zero
        enddo
        voltotal=zero

        do ifine=stenlo(1),stenhi(1)
         call intersect_weight_avg(ic,ifine,bfact_c,bfact_f,wt(1))
         if (wt(1).gt.zero) then
          do jfine=stenlo(2),stenhi(2)
           call intersect_weight_avg(jc,jfine,bfact_c,bfact_f,wt(2))
           if (wt(2).gt.zero) then
            do kfine=stenlo(3),stenhi(3)
             if (SDIM.eq.3) then
              call intersect_weight_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
             endif
             if (wt(SDIM).gt.zero) then
              volall=wt(1)
              do dir2=2,SDIM
               volall=volall*wt(dir2)
              enddo
              if (volall.le.zero) then
               print *,"volall invalid"
               stop
              endif

              do iten=1,nten
               fine_test=NINT(fine(D_DECL(ifine,jfine,kfine),iten))
               coarse_test=NINT(crse_value(iten)) ! crse_value init to 0

               if (avgdown_sweep.eq.0) then

                if (coarse_test.eq.0) then 
                 coarse_test=fine_test
                else if (fine_test.eq.0) then
                 ! do nothing
                else if ((fine_test.eq.1).or.(fine_test.eq.-1)) then
                 coarse_test=fine_test
                else
                 print *,"fine_test invalid"
                 stop
                endif

                if ((coarse_test.eq.1).or. &
                    (coarse_test.eq.-1)) then
                 crse_value(iten)=coarse_test
                else if (coarse_test.eq.0) then
                 ! do nothing (crse_value init. to zero)
                else
                 print *,"coarse_test invalid"
                 stop
                endif

               else if (avgdown_sweep.eq.1) then

                if ((coarse_test.eq.0).or. &
                    (coarse_test.eq.1).or. &
                    (coarse_test.eq.-1)) then

                 if (fine_test.eq.0) then
                  ! do nothing
                 else if (fine_test.eq.coarse_test) then
                  velwt(iten)=velwt(iten)+volall
                  do dir2=1,SDIM
                   local_comp=nten+(iten-1)*SDIM+dir2
                   crse_value(local_comp)=crse_value(local_comp)+ &
                    volall*fine(D_DECL(ifine,jfine,kfine),local_comp)
                  enddo
                 else
                  print *,"fine_test bad (avgdown burning): ",fine_test
                  stop
                 endif
                else
                 print *,"coarse_test invalid"
                 stop
                endif
               else 
                print *,"avgdown_sweep invalid"
                stop
               endif
              enddo ! iten=1..nten

              voltotal=voltotal+volall
             endif ! wt(sdim).gt.0
            enddo ! kfine
           endif
          enddo ! jfine
         endif
        enddo ! ifine

        if (voltotal.le.zero) then
         print *,"voltotal invalid"
         stop
        endif

        if (avgdown_sweep.eq.0) then
         ! do nothing
        else if (avgdown_sweep.eq.1) then

         do iten=1,nten

          coarse_test=NINT(crse_value(iten))
          if (coarse_test.eq.0) then
           if (velwt(iten).eq.zero) then
            ! do nothing
           else
            print *,"velwt invalid"
            stop
           endif
           crse(D_DECL(ic,jc,kc),iten)=zero
           do dir2=1,SDIM
            crse(D_DECL(ic,jc,kc),nten+(iten-1)*SDIM+dir2)=zero
           enddo
          else if ((coarse_test.eq.1).or. &
                   (coarse_test.eq.-1)) then
           if (velwt(iten).gt.zero) then
            crse(D_DECL(ic,jc,kc),iten)=coarse_test
            do dir2=1,SDIM
             crse(D_DECL(ic,jc,kc),nten+(iten-1)*SDIM+dir2)= &
              crse_value(nten+(iten-1)*SDIM+dir2)/velwt(iten)
            enddo
           else
            print *,"velwt invalid"
            stop
           endif
          else
           print *,"coarse_test invalid"
           stop
          endif

         enddo ! iten=1..nten

        else 
         print *,"avgdown_sweep invalid"
         stop
        endif

       enddo ! avgdown_sweep=0..1

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine FORT_AVGDOWN_BURNING

        ! icurv=(iten-1)*(5+SDIM)
        ! dir=1..sdim
        ! side=-1 or 1
        ! curvfab(D_DECL(i,j,k),icurv+4+SDIM)=dir*side
      subroutine FORT_AVGDOWN_CURV( &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp,nmat,nten, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) 

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module
      use navierstokesf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level_c
      INTEGER_T, intent(in) :: level_f
      INTEGER_T, intent(in) :: bfact_c
      INTEGER_T, intent(in) :: bfact_f
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: dxf(SDIM)
      REAL_T, intent(in) :: xlo_fine(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: ncomp,nmat,nten
      INTEGER_T nten_test
      INTEGER_T, intent(in) :: DIMDEC(crse)
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      INTEGER_T, intent(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      REAL_T, intent(out) :: crse(DIMV(crse),ncomp)
      REAL_T, intent(in) :: fine(DIMV(fine),ncomp)
      INTEGER_T ic,jc,kc
      INTEGER_T ifine,jfine,kfine
      INTEGER_T dir2
      INTEGER_T iten
      INTEGER_T n
      REAL_T voltotal
      REAL_T volall
      REAL_T wt(SDIM)
      REAL_T crse_value(ncomp)
      INTEGER_T fine_test
      INTEGER_T coarse_test
      REAL_T velwt(nten)
      INTEGER_T icurv,istat

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten.ne.nten_test) then
       print *,"nten invalid avgdown_curv nten nten_test ",nten,nten_test
       stop
      endif

      if (ncomp.ne.nten*(SDIM+5)) then
       print *,"ncomp invalid35"
       stop
      endif

      if (bfact_f.lt.1) then
       print *,"bfact_f invalid4 ",bfact_f
       stop
      endif
      if (bfact_c.lt.1) then
       print *,"bfact_c invalid4 ",bfact_c
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level_f-1)) then
       print *,"level_c or level_f invalid"
       stop
      endif

      call checkbound(lo,hi,DIMS(crse),0,-1,411)
      call checkbound(lof,hif,DIMS(fine),0,-1,411)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)

       do n=1,ncomp
        crse_value(n)=zero
       enddo
       do iten=1,nten
        velwt(iten)=zero
       enddo
       voltotal=zero

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
        bfact_c,bfact_f)

       do ifine=stenlo(1),stenhi(1)
        call intersect_weight_avg(ic,ifine,bfact_c,bfact_f,wt(1))
        if (wt(1).gt.zero) then
         do jfine=stenlo(2),stenhi(2)
          call intersect_weight_avg(jc,jfine,bfact_c,bfact_f,wt(2))
          if (wt(2).gt.zero) then
           do kfine=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             volall=wt(1)
             do dir2=2,SDIM
              volall=volall*wt(dir2)
             enddo
             if (volall.le.zero) then
              print *,"volall invalid"
              stop
             endif

             do iten=1,nten
              icurv=(iten-1)*(5+SDIM)
              istat=icurv+4+SDIM ! dir x side  dir=1..sdim side=-1 or 1
              fine_test=NINT(fine(D_DECL(ifine,jfine,kfine),istat))
              coarse_test=NINT(crse_value(istat))
              if ((coarse_test.ne.0).and. &
                  (coarse_test.ne.SDIM+1)) then
               print *,"coarse_test invalid"
               stop
              endif
              if (fine_test.eq.0) then
               ! do nothing
              else if ((abs(fine_test).ge.1).and. &
                       (abs(fine_test).le.SDIM+1)) then
               velwt(iten)=velwt(iten)+volall
               crse_value(istat)=SDIM+1
               crse_value(istat+1)=zero ! im3
               do dir2=1,SDIM+3
                crse_value(icurv+dir2)= &
                  crse_value(icurv+dir2)+volall* &
                  fine(D_DECL(ifine,jfine,kfine),icurv+dir2)
               enddo
              else
               print *,"fine_test bad (avgdown curv): ",fine_test
               stop
              endif
             enddo ! iten

             voltotal=voltotal+volall
            endif ! wt(sdim).gt.0
           enddo ! kfine
          endif
         enddo ! jfine
        endif
       enddo ! ifine

       if (voltotal.le.zero) then
        print *,"voltotal invalid"
        stop
       endif

       do iten=1,nten
        icurv=(iten-1)*(5+SDIM)
        istat=icurv+4+SDIM ! dir x side  dir=1..sdim side=-1 or 1
        coarse_test=NINT(crse_value(istat))
        if (coarse_test.eq.0) then
         if (velwt(iten).ne.zero) then
          print *,"velwt invalid"
          stop
         endif
         do dir2=1,5+SDIM
          crse(D_DECL(ic,jc,kc),icurv+dir2)=zero
         enddo
        else if (coarse_test.eq.SDIM+1) then
         if (velwt(iten).gt.zero) then
          crse(D_DECL(ic,jc,kc),istat)=SDIM+1 ! dir x side
          crse(D_DECL(ic,jc,kc),istat+1)=zero ! im3
          do dir2=1,3+SDIM
           crse(D_DECL(ic,jc,kc),icurv+dir2)= &
            crse_value(icurv+dir2)/velwt(iten)
          enddo
         else
          print *,"velwt invalid"
          stop
         endif
        else
         print *,"coarse_test invalid"
         stop
        endif

       enddo ! iten=1..nten

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine FORT_AVGDOWN_CURV


      subroutine FORT_MOFAVGDOWN ( &
       time, &
       problo, &
       dxc, &
       dxf, &
       bfact_c,bfact_f, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi,nmat)

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probcommon_module

      IMPLICIT NONE

      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: bfact_c,bfact_f
      INTEGER_T, intent(in) :: DIMDEC(crse)
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T  growlo(3),growhi(3)
      INTEGER_T  stenlo(3),stenhi(3)
      REAL_T, intent(out) :: crse(DIMV(crse),nmat*ngeom_raw)
      REAL_T, intent(in) :: fine(DIMV(fine),nmat*ngeom_raw)
      INTEGER_T ifine,jfine,kfine
      INTEGER_T ic,jc,kc
      INTEGER_T dir
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: dxc(SDIM)
      REAL_T, intent(in) :: dxf(SDIM)
      INTEGER_T domlo(SDIM)
      INTEGER_T nhalf
      INTEGER_T nhalfgrid
      INTEGER_T n_overlap
 
      INTEGER_T im
      INTEGER_T vofcomp_raw
      INTEGER_T vofcomp_recon
      REAL_T wt(SDIM)
      REAL_T testwt

      REAL_T temp_vfrac
      REAL_T temp_cen(SDIM)

      REAL_T xstencoarse(-3:3,SDIM)
      REAL_T xstenfine(-3:3,SDIM)
      REAL_T xstengrid(-1:1,SDIM)

      REAL_T mofdatafine(nmat*ngeom_recon)
      REAL_T mofdatacoarse(nmat*ngeom_recon)
      REAL_T multi_centroidA(nmat,SDIM)

      REAL_T volcoarse
      REAL_T cencoarse(SDIM)
      REAL_T volfine
      REAL_T cenfine(SDIM)

      REAL_T multi_volume(nmat)
      REAL_T multi_cen(SDIM,nmat)
      INTEGER_T tessellate
      INTEGER_T continuous_mof
      INTEGER_T mof_verbose
      INTEGER_T use_ls_data
      INTEGER_T fine_covered
      REAL_T LS_stencil(D_DECL(-1:1,-1:1,-1:1),nmat)
      INTEGER_T nmax
      INTEGER_T tid
#ifdef _OPENMP
      INTEGER_T omp_get_thread_num
#endif

      tid=0       
#ifdef _OPENMP
      tid=omp_get_thread_num()
#endif
      if ((tid.ge.geom_nthreads).or.(tid.lt.0)) then
       print *,"tid invalid"
       stop
      endif 

      nmax=POLYGON_LIST_MAX ! in: MOFAVGDOWN

      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (bfact_f.lt.1) then
       print *,"bfact_f invalid5 ",bfact_f
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      do dir=1,SDIM
       domlo(dir)=0
      enddo
      nhalf=3
      nhalfgrid=1

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)

        ! coarse centroids and volume fractions are initialized
        ! to zero.
       do dir=1,nmat*ngeom_raw
        crse(D_DECL(ic,jc,kc),dir) = zero
       enddo
       do dir=1,nmat*ngeom_recon
        mofdatacoarse(dir) = zero
       enddo

       volcoarse=zero
       do dir=1,SDIM
        cencoarse(dir)=zero
       enddo

       n_overlap=0

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
         bfact_c,bfact_f)

       do ifine=stenlo(1),stenhi(1)
        call intersect_weight_avg(ic,ifine,bfact_c,bfact_f,wt(1))
        if (1.eq.0) then
         print *,"ic,ifine,wt(1) ",ic,ifine,wt(1)
        endif
        if (wt(1).gt.zero) then
         do jfine=stenlo(2),stenhi(2)
          call intersect_weight_avg(jc,jfine,bfact_c,bfact_f,wt(2))
          if (1.eq.0) then
           print *,"jc,jfine,wt(2) ",jc,jfine,wt(2)
          endif
          if (wt(2).gt.zero) then
           do kfine=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then
             if (SDIM.eq.2) then
              testwt=wt(1)*wt(2)
             else if (SDIM.eq.3) then
              testwt=wt(1)*wt(2)*wt(SDIM)
             else
              print *,"dimension bust"
              stop
             endif

             if (testwt.gt.zero) then

              n_overlap=n_overlap+1

              do im=1,nmat
               vofcomp_raw=(im-1)*ngeom_raw+1
               vofcomp_recon=(im-1)*ngeom_recon+1
               do dir=1,SDIM+1
                mofdatafine(vofcomp_recon+dir-1)= &
                 fine(D_DECL(ifine,jfine,kfine),vofcomp_raw+dir-1)
               enddo
               do dir=SDIM+2,ngeom_recon
                mofdatafine(vofcomp_recon+dir-1)=zero
               enddo
              enddo ! im

              call gridsten(xstencoarse,problo,ic,jc,kc, &
               domlo,bfact_c,dxc,nhalf)
              call gridsten(xstenfine,problo,ifine,jfine,kfine, &
               domlo,bfact_f,dxf,nhalf)

              call Box_volumeFAST(bfact_f,dxf,xstenfine,nhalf, &
               volfine,cenfine,SDIM)

              fine_covered=1

              do dir=1,SDIM
               xstengrid(-1,dir)=max(xstencoarse(-1,dir),xstenfine(-1,dir))
               xstengrid(1,dir)=min(xstencoarse(1,dir),xstenfine(1,dir))
               xstengrid(0,dir)=half*(xstengrid(-1,dir)+xstengrid(1,dir))

               if (xstengrid(-1,dir).gt.xstenfine(-1,dir)+VOFTOL*dxf(dir)) then
                fine_covered=0
               endif
               if (xstengrid(1,dir).lt.xstenfine(1,dir)-VOFTOL*dxf(dir)) then
                fine_covered=0
               endif

               if ((bfact_f.eq.1).and.(bfact_c.eq.1)) then
                if ((abs(xstengrid(-1,dir)-xstenfine(-1,dir)).gt. &
                     VOFTOL*dxf(dir)).or. &
                    (abs(xstengrid(1,dir)-xstenfine(1,dir)).gt. &
                     VOFTOL*dxf(dir))) then
                 print *,"fine cell should be completely covered by coarse"
                 stop
                endif
               endif
              enddo ! dir

              if ((bfact_f.eq.1).and.(bfact_c.eq.1)) then
               if (fine_covered.ne.1) then
                print *,"fine cell should be completely covered by coarse"
                stop
               endif
              endif

              if (fine_covered.eq.1) then

               do im=1,nmat
                vofcomp_raw=(im-1)*ngeom_raw+1
                multi_volume(im)= &
                 fine(D_DECL(ifine,jfine,kfine),vofcomp_raw)*volfine
                do dir=1,SDIM
                 multi_cen(dir,im)= &
                  fine(D_DECL(ifine,jfine,kfine),vofcomp_raw+dir)+cenfine(dir)
                enddo
               enddo ! im

              else if (fine_covered.eq.0) then
               use_ls_data=0
               mof_verbose=0
               continuous_mof=0

                ! sum F_fluid=1  sum F_solid<=1
               call make_vfrac_sum_ok_base(mofdatafine,nmat,SDIM,304)

               call multimaterial_MOF( &
                bfact_f,dxf,xstenfine,nhalf, &
                mof_verbose, &
                use_ls_data, &
                LS_stencil, &
                geom_xtetlist(1,1,1,tid+1), &
                geom_xtetlist(1,1,1,tid+1), &
                nmax, &
                nmax, &
                mofdatafine, &
                multi_centroidA, &
                continuous_mof, &
                nmat,SDIM,3)

               tessellate=0
               call multi_get_volume_grid_simple( &
                tessellate, &
                bfact_f,dxf,xstenfine,nhalf, &
                mofdatafine, &
                xstengrid,nhalfgrid, &
                multi_volume,multi_cen, &
                geom_xtetlist(1,1,1,tid+1), &
                nmax, &
                nmax, &
                nmat,SDIM,6)

              else
               print *,"fine_covered invalid"
               stop
              endif

              do im=1,nmat
               vofcomp_recon=(im-1)*ngeom_recon+1
               mofdatacoarse(vofcomp_recon)=mofdatacoarse(vofcomp_recon)+ &
                multi_volume(im)
               do dir=1,SDIM
                mofdatacoarse(vofcomp_recon+dir)= &
                 mofdatacoarse(vofcomp_recon+dir)+ &
                 multi_cen(dir,im)*multi_volume(im)
               enddo
               if (is_rigid(nmat,im).eq.0) then
                volcoarse=volcoarse+multi_volume(im)
                do dir=1,SDIM
                 cencoarse(dir)=cencoarse(dir)+ &
                  multi_cen(dir,im)*multi_volume(im)
                enddo
               else if (is_rigid(nmat,im).eq.1) then
                ! do nothing
               else
                print *,"is_rigid(nmat,im) invalid"
                stop
               endif
              enddo ! im=1..nmat
  
             else if (testwt.eq.zero) then
              print *,"testwt should not be 0"
              stop
             else
              print *,"testwt invalid"
              stop
             endif

            endif ! wt(sdim)>0
           enddo ! kfine
          endif ! wt(2)>0
         enddo ! jfine
        endif ! wt(1)>0
       enddo ! ifine

       if ((bfact_f.eq.1).and.(bfact_c.eq.1)) then
        if (n_overlap.ne.4*(SDIM-1)) then
         print *,"n_overlap invalid"
         stop
        endif
       else if (n_overlap.lt.1) then
        print *,"n_overlap invalid"
        stop
       endif

       if (volcoarse.le.zero) then
        print *,"volcoarse invalid:",volcoarse
        print *,"ic,jc,kc,bfact_c,bfact_f ",ic,jc,kc,bfact_c,bfact_f
        do dir=1,SDIM
         print *,"dir,stenlo,stenhi ",dir,stenlo(dir),stenhi(dir)
        enddo
        dir=3
        print *,"dir=3,stenlo,stenhi ",dir,stenlo(dir),stenhi(dir)
        stop
       endif

       do dir=1,SDIM
        cencoarse(dir)=cencoarse(dir)/volcoarse
       enddo

       do im=1,nmat
        vofcomp_recon=(im-1)*ngeom_recon+1
        vofcomp_raw=(im-1)*ngeom_raw+1

        temp_vfrac=mofdatacoarse(vofcomp_recon)/volcoarse
        if ((temp_vfrac.lt.zero).or. &
            (temp_vfrac.gt.one+VOFTOL)) then
         print *,"temp_vfrac invalid"
         stop
        endif
        if (temp_vfrac.le.VOFTOL) then
         temp_vfrac=zero
         do dir=1,SDIM
          temp_cen(dir)=zero
         enddo
        else if (temp_vfrac.ge.one-VOFTOL) then
         temp_vfrac=one
         do dir=1,SDIM
          temp_cen(dir)=zero
         enddo
        else if ((temp_vfrac.gt.zero).and.(temp_vfrac.lt.one)) then
         do dir=1,SDIM
          temp_cen(dir)= &
            mofdatacoarse(vofcomp_recon+dir)/ &
            mofdatacoarse(vofcomp_recon)-cencoarse(dir)
         enddo
        else
         print *,"temp_vfrac invalid"
         stop
        endif

        crse(D_DECL(ic,jc,kc),vofcomp_raw) = temp_vfrac
        do dir=1,SDIM
         crse(D_DECL(ic,jc,kc),vofcomp_raw+dir) = temp_cen(dir)
        enddo
       enddo ! im

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine FORT_MOFAVGDOWN


      subroutine FORT_ERRORAVGDOWN ( &
       problo, &
       dxf, &
       bfact_c,bfact_f, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi)

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: bfact_c,bfact_f
      INTEGER_T, intent(in) :: DIMDEC(crse)
      INTEGER_T, intent(in) :: DIMDEC(fine)
      INTEGER_T, intent(in) :: lo(SDIM),hi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T stenlo(3),stenhi(3)
      REAL_T, intent(out) :: crse(DIMV(crse))
      REAL_T, intent(in) :: fine(DIMV(fine))
      INTEGER_T i,j,k,ic,jc,kc
      REAL_T, intent(in) :: problo(SDIM),dxf(SDIM)
 
      REAL_T wt(SDIM)

      INTEGER_T sanity_sum
      REAL_T test_error,max_error

      if (bfact_f.lt.1) then
       print *,"bfact_f invalid5 ",bfact_f
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)

       crse(D_DECL(ic,jc,kc)) = zero

       max_error=zero
       sanity_sum=0

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi,bfact_c,bfact_f)
       do i=stenlo(1),stenhi(1)
        call intersect_weight_avg(ic,i,bfact_c,bfact_f,wt(1))
        if (1.eq.0) then
         print *,"ic,i,wt(1) ",ic,i,wt(1)
        endif
        if (wt(1).gt.zero) then
         do j=stenlo(2),stenhi(2)
          call intersect_weight_avg(jc,j,bfact_c,bfact_f,wt(2))
          if (1.eq.0) then
           print *,"jc,j,wt(2) ",jc,j,wt(2)
          endif
          if (wt(2).gt.zero) then
           do k=stenlo(3),stenhi(3)
            if (SDIM.eq.3) then
             call intersect_weight_avg(kc,k,bfact_c,bfact_f,wt(SDIM))
            endif
            if (wt(SDIM).gt.zero) then

             sanity_sum=sanity_sum+1

             test_error=fine(D_DECL(i,j,k))
             if (test_error.lt.zero) then
              print *,"error should not be <0"
              stop
             endif
             if (test_error.gt.max_error) then
              max_error=test_error
             endif
            endif ! wt(sdim)>0
           enddo ! k
          endif ! wt(2)>0
         enddo ! j
        endif ! wt(1)>0
       enddo ! i

       if (sanity_sum.le.0) then
        print *,"sanity_sum invalid"
        stop
       endif

       crse(D_DECL(ic,jc,kc)) = max_error

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine FORT_ERRORAVGDOWN


       subroutine FORT_SUMMASS( &
        tid, &
        adapt_quad_depth, &
        slice_dir,xslice, &
        problo,probhi, &
        xlo,dx, &
        cellten,DIMS(cellten), &
        lsfab,DIMS(lsfab), &
        maskSEM,DIMS(maskSEM), &
        mask,DIMS(mask), &
        drag,DIMS(drag), &
        slopes,DIMS(slopes), &
        den,DIMS(den), &
        vel,DIMS(vel), &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        time, &
        local_result, &
        resultALL, &
        sumdata_type, &
        sumdata_sweep, &
        resultsize, &
        NN,ZZ,FF, &
        dirx,diry,cut_flag, &
        nmat, &
        ntensorMM,  &
        den_ncomp, &
        isweep)

       use LegendreNodes
       use global_utility_module
       use probf90_module
       use geometry_intersect_module
       use MOF_routines_module

       IMPLICIT NONE

       INTEGER_T tid
       INTEGER_T adapt_quad_depth
       INTEGER_T max_level
       INTEGER_T slice_dir
       REAL_T xslice(SDIM)
       INTEGER_T resultsize
       INTEGER_T den_ncomp
       INTEGER_T nmat
       INTEGER_T ntensorMM
       INTEGER_T isweep
       INTEGER_T NN,dirx,diry,cut_flag
       REAL_T ZZ(0:NN)
       REAL_T FF(0:NN)

       REAL_T problo(SDIM)
       REAL_T probhi(SDIM)
       INTEGER_T bfact
       INTEGER_T tilelo(SDIM),tilehi(SDIM)
       INTEGER_T fablo(SDIM),fabhi(SDIM)
       INTEGER_T growlo(3),growhi(3)
       INTEGER_T DIMDEC(cellten)
       INTEGER_T DIMDEC(lsfab)
       INTEGER_T DIMDEC(mask)
       INTEGER_T DIMDEC(maskSEM)
       INTEGER_T DIMDEC(drag)
       INTEGER_T DIMDEC(slopes)
       INTEGER_T DIMDEC(den)
       INTEGER_T DIMDEC(vel)
       REAL_T  local_result(resultsize)
       REAL_T  resultALL(resultsize)
       INTEGER_T sumdata_type(resultsize)
       INTEGER_T sumdata_sweep(resultsize)
       REAL_T  time
       REAL_T  cellten(DIMV(cellten),ntensorMM)  
       REAL_T  lsfab(DIMV(lsfab),nmat)  
       REAL_T  maskSEM(DIMV(maskSEM))
       REAL_T  mask(DIMV(mask))
       REAL_T  drag(DIMV(drag),4*SDIM+1)
       REAL_T  slopes(DIMV(slopes),nmat*ngeom_recon)  
       REAL_T  den(DIMV(den),den_ncomp)  
       REAL_T  vel(DIMV(vel),num_materials_vel*(SDIM+1)) ! includes pressure 
       REAL_T  xlo(SDIM),dx(SDIM)

       INTEGER_T i,j,k
       INTEGER_T ii,jj,kk
       INTEGER_T dir
       INTEGER_T im
       INTEGER_T im_primary
       INTEGER_T vofcomp
       INTEGER_T dirMAC
       REAL_T xcen,xcrit
       REAL_T xsten(-3:3,SDIM)
       REAL_T xstenMAC(-3:3,SDIM)
       REAL_T mofdata(nmat*ngeom_recon)
       REAL_T mofdata_tess(nmat*ngeom_recon)
       INTEGER_T level,dir2,dir3
       REAL_T errorparm(nmat*2) ! fi,ei
       REAL_T xbottom,xtop,ls_above,ls_below,ZZgrid
       REAL_T vof_below,vof_above,vof_face
       REAL_T volgrid,distbound
       REAL_T cengrid(SDIM)
       REAL_T xboundary(SDIM)

       INTEGER_T filler_comp,FE_sum_comp,drag_sum_comp
       INTEGER_T minint_sum_comp,maxint_sum_comp
       INTEGER_T pdrag_sum_comp,minden_sum_comp,maxden_sum_comp
       INTEGER_T xnot_amp_sum_comp,cen_sum_comp
       INTEGER_T mincen_sum_comp,maxcen_sum_comp
       INTEGER_T mass_sum_comp,mom_sum_comp,energy_sum_comp
       INTEGER_T left_pressure_comp
       INTEGER_T kinetic_energy_comp
       INTEGER_T LS_F_sum_comp,LS_cen_sum_comp
       INTEGER_T torque_sum_comp,ptorque_sum_comp
       INTEGER_T step_perim_sum_comp
       INTEGER_T minint_slice
       INTEGER_T maxint_slice
       INTEGER_T vort_sum_comp
       INTEGER_T vort_error
       INTEGER_T vel_error
       INTEGER_T energy_moment
       INTEGER_T enstrophy
       INTEGER_T total_comp

       INTEGER_T isrc,idest
       INTEGER_T iside
       REAL_T dz_external
       INTEGER_T j_external

       REAL_T cen_material(SDIM)
       INTEGER_T use_vof_height

       REAL_T ldata(D_DECL(3,3,3))
       REAL_T LSvolume,LSfacearea
       REAL_T LScentroid(SDIM)
       REAL_T LScen_material(SDIM)
       REAL_T LSareacentroid(SDIM)
       INTEGER_T i1,j1,k1,k1lo,k1hi
       REAL_T KECELL,dencore,Tcore,ecore,totalE
       INTEGER_T in_slice,nhalf
       INTEGER_T ux,vx,wx,uy,vy,wy,uz,vz,wz,nbase,veldir
       REAL_T gradu(3,3)
       REAL_T vort(3)

       REAL_T local_vort
       REAL_T local_vort_error
       REAL_T vort_expect
       REAL_T local_vel_error
       REAL_T local_energy_moment
       REAL_T local_vel(SDIM)
       REAL_T vel_expect(SDIM)
       REAL_T local_xsten(SDIM)

       REAL_T local_kinetic_energy(nmat)

       REAL_T rr
       INTEGER_T mask_im
       INTEGER_T cell_index,element_index,sub_index
       INTEGER_T velcomp
       INTEGER_T prescomp
       INTEGER_T nmax

       REAL_T LS_LOCAL(nmat)

       nhalf=3
       nmax=POLYGON_LIST_MAX  ! in: SUMMASS

       if ((adapt_quad_depth.lt.1).or.(adapt_quad_depth.gt.10)) then
        print *,"adapt_quad_depth invalid"
        stop
       endif
       if ((slice_dir.lt.0).or.(slice_dir.ge.SDIM)) then
        print *,"slice_dir invalid"
        stop
       endif
       if (nmat.ne.num_materials) then
        print *,"nmat invalid"
        stop
       endif
       if (num_materials_vel.ne.1) then
        print *,"num_materials_vel invalid"
        stop
       endif
       if (ngeom_recon.ne.2*SDIM+3) then
        print *,"ngeom_recon invalid sum mass ",ngeom_recon
        stop
       endif
       if (den_ncomp.ne.nmat*num_state_material) then
        print *,"den_ncomp invalid"
        stop
       endif
       if (ntensorMM.ne.SDIM*SDIM*num_materials_vel) then
        print *,"ntensorMM invalid"
        stop
       endif
       if ((dirx.lt.0).or.(dirx.ge.SDIM)) then
        print *,"dirx invalid"
        stop
       endif
       if ((diry.lt.0).or.(diry.ge.SDIM).or.(diry.eq.dirx)) then
        print *,"diry invalid"
        stop
       endif
       if ((cut_flag.ne.0).and.(cut_flag.ne.1)) then
        print *,"cut_flag invalid"
        stop
       endif
       if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
        print *,"tid invalid"
        stop
       endif

       filler_comp=0
       FE_sum_comp=filler_comp+1
       drag_sum_comp=FE_sum_comp+2*nmat
       minint_sum_comp=drag_sum_comp+3
       maxint_sum_comp=minint_sum_comp+3*nmat
       pdrag_sum_comp=maxint_sum_comp+3*nmat
       minden_sum_comp=pdrag_sum_comp+3
       maxden_sum_comp=minden_sum_comp+2*nmat
       xnot_amp_sum_comp=maxden_sum_comp+2*nmat
       cen_sum_comp=xnot_amp_sum_comp+1
       mincen_sum_comp=cen_sum_comp+3*nmat
       maxcen_sum_comp=mincen_sum_comp+nmat
       mass_sum_comp=maxcen_sum_comp+nmat
       mom_sum_comp=mass_sum_comp+nmat
       energy_sum_comp=mom_sum_comp+3*nmat
       left_pressure_comp=energy_sum_comp+nmat
       kinetic_energy_comp=left_pressure_comp+4
       LS_F_sum_comp=kinetic_energy_comp+nmat
       LS_cen_sum_comp=LS_F_sum_comp+nmat
       torque_sum_comp=LS_cen_sum_comp+3*nmat
       ptorque_sum_comp=torque_sum_comp+3
       step_perim_sum_comp=ptorque_sum_comp+3
       minint_slice=step_perim_sum_comp+1
       maxint_slice=minint_slice+nmat
       vort_sum_comp=maxint_slice+nmat
       vort_error=vort_sum_comp+3
       vel_error=vort_error+1
       energy_moment=vel_error+1
       enstrophy=energy_moment+1 ! integral of w dot w
       total_comp=enstrophy+nmat

       if (resultsize.ne.total_comp) then
        print *,"mismatch between resultsize and total_comp"
        stop
       endif
       do idest=1,total_comp
        if ((sumdata_type(idest).ne.1).and. &
            (sumdata_type(idest).ne.2).and. &
            (sumdata_type(idest).ne.3)) then
         print *,"sumdata_type invalid"
         stop
        endif
        if ((sumdata_sweep(idest).ne.0).and. &
            (sumdata_sweep(idest).ne.1)) then
         print *,"sumdata_sweep invalid"
         stop
        endif
       enddo  ! idest

       if (isweep.eq.0) then
        max_level=adapt_quad_depth
       else if (isweep.eq.1) then
        max_level=1
       else
        print *,"isweep invalid"
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
       if (bfact.lt.1) then
        print *,"bfact too small"
        stop
       endif

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       ! compute u_x,v_x,w_x, u_y,v_y,w_y, u_z,v_z,w_z;  
       call tensorcomp_matrix(ux,uy,uz,vx,vy,vz,wx,wy,wz)

       call checkbound(fablo,fabhi,DIMS(cellten),0,-1,411) 
       call checkbound(fablo,fabhi,DIMS(lsfab),2,-1,411) 
       call checkbound(fablo,fabhi,DIMS(maskSEM),1,-1,411) 
       call checkbound(fablo,fabhi,DIMS(mask),2,-1,411) 
       call checkbound(fablo,fabhi,DIMS(drag),0,-1,413) 
       call checkbound(fablo,fabhi,DIMS(slopes),2,-1,413) 
       call checkbound(fablo,fabhi,DIMS(den),1,-1,413) 
       call checkbound(fablo,fabhi,DIMS(vel),1,-1,413) 

       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        if (mask(D_DECL(i,j,k)).gt.zero) then

         call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

         do veldir=1,3
         do dir=1,3
          gradu(veldir,dir)=zero
         enddo
         enddo

         do dir=1,SDIM
          if (dir.eq.1) then
           nbase=ux-1
          else if (dir.eq.2) then
           nbase=uy-1
          else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
           nbase=uz-1
          else
           print *,"dir invalid summass"
           stop
          endif

          do veldir=1,SDIM
           gradu(veldir,dir)=cellten(D_DECL(i,j,k),nbase+veldir) 
          enddo
         enddo ! dir

         velcomp=1

         if (levelrz.eq.0) then
          ! do nothing
         else if (levelrz.eq.1) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          rr=xsten(0,1)
          gradu(3,3)=vel(D_DECL(i,j,k),velcomp)/abs(rr)
         else if (levelrz.eq.3) then
          rr=xsten(0,1)
          gradu(2,2)=gradu(2,2)+vel(D_DECL(i,j,k),velcomp)/abs(rr)
          gradu(1,2)=gradu(1,2)-vel(D_DECL(i,j,k),velcomp+1)/abs(rr)
         else
          print *,"levelrz invalid summass"
          stop
         endif

         vort(1)=gradu(3,2)-gradu(2,3)
         vort(2)=gradu(1,3)-gradu(3,1)
         vort(3)=gradu(2,1)-gradu(1,2)

         do dir=1,nmat*ngeom_recon
          mofdata(dir)=slopes(D_DECL(i,j,k),dir)
          mofdata_tess(dir)=mofdata(dir)
         enddo

         level=0

         ! before (mofdata): fluids tessellate
         ! after  (mofdata): fluids and solids tessellate
         call multi_get_volume_tessellate( &
          bfact,dx,xsten,nhalf, &
          mofdata_tess, &
          geom_xtetlist(1,1,1,tid+1), &
          nmax, &
          nmax, &
          nmat, &
          SDIM, &
          101)

          ! tessellate==1 (internal to stackerror)
          ! in: SUMMASS
         call stackerror( &
          geom_xtetlist(1,1,1,tid+1), &
          xsten,nhalf,dx,bfact, &
          xsten,nhalf, &
          mofdata, &
          mofdata_tess, &
          errorparm,level,max_level,nmat,time)

          ! F1,E1,F2,E2,F3,E3,...
         do dir=1,2*nmat
          idest=FE_sum_comp+dir
          local_result(idest)=local_result(idest)+errorparm(dir)
         enddo

         call Box_volumeFAST(bfact,dx,xsten,nhalf, &
          volgrid,cengrid,SDIM)

         mask_im=NINT(maskSEM(D_DECL(i,j,k)))
         if ((mask_im.ge.1).and.(mask_im.le.nmat)) then
          volgrid=one
          do dir=1,SDIM
           if (dir.eq.1) then
            cell_index=i
           else if (dir.eq.2) then
            cell_index=j
           else if ((dir.eq.3).and.(SDIM.eq.3)) then
            cell_index=k
           else
            print *,"dir invalid summass 2"
            stop
           endif
           call get_element_index(cell_index,element_index,sub_index,bfact)
           volgrid=volgrid*cache_gauss_w(bfact,sub_index,SPTYPE)* &
            bfact*dx(dir)/two
          enddo ! dir=1..sdim

          if (levelrz.eq.0) then
           ! do nothing
          else if (levelrz.eq.1) then
           if (SDIM.ne.2) then
            print *,"dimension bust"
            stop
           endif
           volgrid=volgrid*two*Pi*xsten(0,1)
          else if (levelrz.eq.3) then
           volgrid=volgrid*xsten(0,1)
          else
           print *,"levelrz invalid"
           stop
          endif
          
         else if (mask_im.eq.0) then
          ! do nothing
         else
          print *,"maskSEM invalid"
          stop
         endif

         do im=1,nmat

          LS_LOCAL(im)=lsfab(D_DECL(i,j,k),im)

          do i1=-1,1
          do j1=-1,1
          do k1=k1lo,k1hi
           ldata(D_DECL(i1+2,j1+2,k1+2))= &
            lsfab(D_DECL(i+i1,j+j1,k+k1),im)
          enddo
          enddo
          enddo
           ! LSvolume is a volume fraction.
           ! LScentroid is in an absolute coordinate system.
          call getvolume(bfact,dx,xsten,nhalf, &
           ldata,LSvolume,LSfacearea, &
           LScentroid,LSareacentroid,VOFTOL,SDIM)

          vofcomp=(im-1)*ngeom_recon+1
          do dir=1,SDIM
           cen_material(dir)=mofdata_tess(vofcomp+dir)+cengrid(dir)
           LScen_material(dir)=LScentroid(dir)
          enddo

          do dir=1,SDIM
           idest=cen_sum_comp+dir+3*(im-1)
           local_result(idest)=local_result(idest)+ &
             volgrid*mofdata_tess(vofcomp)*cen_material(dir)
           idest=LS_cen_sum_comp+dir+3*(im-1)
           local_result(idest)=local_result(idest)+ &
             volgrid*LSvolume*LScen_material(dir)
          enddo
          idest=LS_F_sum_comp+im
          local_result(idest)=local_result(idest)+volgrid*LSvolume

          dencore=den(D_DECL(i,j,k),1+num_state_material*(im-1))
          Tcore=den(D_DECL(i,j,k),2+num_state_material*(im-1))

          idest=mass_sum_comp+im
          local_result(idest)=local_result(idest)+ &
            dencore*volgrid*mofdata_tess(vofcomp)

          KECELL=zero

          do dir=1,SDIM
           velcomp=dir
           idest=mom_sum_comp+3*(im-1)+dir
           local_result(idest)=local_result(idest)+ &
             vel(D_DECL(i,j,k),velcomp)*volgrid*mofdata_tess(vofcomp)*dencore

           KECELL=KECELL+vel(D_DECL(i,j,k),velcomp)**2
          enddo ! dir=1..sdim

          KECELL=half*KECELL

          idest=kinetic_energy_comp+im

          local_kinetic_energy(im)= &
            KECELL*volgrid*mofdata_tess(vofcomp)*dencore

          local_result(idest)=local_result(idest)+local_kinetic_energy(im)

          idest=energy_sum_comp+im

          call INTERNAL_material(dencore,Tcore,ecore, &
            fort_material_type(im),im)

          totalE=dencore*(KECELL+ecore)*volgrid*mofdata_tess(vofcomp)
          local_result(idest)=local_result(idest)+totalE

         enddo  ! im=1..nmat

         call get_primary_material(LS_LOCAL,nmat,im_primary)

         do dir=1,3
          idest=vort_sum_comp+dir
          local_result(idest)=local_result(idest)+volgrid*vort(dir)
          idest=enstrophy+1
          if ((im_primary.ge.1).and.(im_primary.le.nmat)) then
           local_result(idest+im_primary-1)= &
               local_result(idest+im_primary-1)+volgrid*(vort(dir)**2)
          else
           print *,"im_primary invalid"
           stop
          endif
         enddo ! dir=1..3

         local_vort=sqrt(vort(1)**2+vort(2)**2+vort(3)**2)
         do dir=1,SDIM
          local_vel(dir)=vel(D_DECL(i,j,k),dir)
          local_xsten(dir)=xsten(0,dir)
         enddo
         call get_vort_vel_error(time,local_xsten,local_vort,local_vel, &
          local_vort_error,local_vel_error,vort_expect,vel_expect, &
          local_energy_moment)
         
         local_result(energy_moment+1)=local_result(energy_moment+1)+ &
          local_energy_moment*local_kinetic_energy(1)

         if (local_result(vort_error+1).lt.local_vort_error) then
          local_result(vort_error+1)=local_vort_error
         endif
         if (local_result(vel_error+1).lt.local_vel_error) then
          local_result(vel_error+1)=local_vel_error
         endif
         if (isweep.eq.0) then
          ! do nothing
         else if (isweep.eq.1) then
          if (local_vort_error.gt.zero) then
           if (local_vort_error.gt.resultALL(vort_error+1)-VOFTOL) then
            print *,"**** POSITION OF MAX VORT ERROR ****"
            do dir=1,SDIM
             print *,"dir,xpos ",dir,local_xsten(dir)
            enddo
            print *,"time ",time
            print *,"local_vort_error ",local_vort_error
            print *,"EXPECTED mag vort ",vort_expect
            print *,"ACTUAL mag vort ",local_vort
           endif
          else if (local_vort_error.eq.zero) then
           ! do nothing
          else
           print *,"local_vort_error invalid"
           stop
          endif

          if (local_vel_error.gt.zero) then
           if (local_vel_error.gt.resultALL(vel_error+1)-VOFTOL) then
            print *,"**** POSITION OF MAX VEL ERROR ****"
            do dir=1,SDIM
             print *,"dir,xpos ",dir,local_xsten(dir)
            enddo
            print *,"time ",time
            print *,"local_vel_error ",local_vel_error
            do dir=1,SDIM
             print *,"dir, EXPECTED vel ",dir,vel_expect(dir)
             print *,"dir, ACTUAL vel ",dir,local_vel(dir)
            enddo
           endif
          else if (local_vel_error.eq.zero) then
           ! do nothing
          else
           print *,"local_vel_error invalid"
           stop
          endif

         else
          print *,"isweep invalid"
          stop
         endif

         prescomp=SDIM+1

         idest=left_pressure_comp+1
         if (xsten(-1,1).le.problox+VOFTOL*dx(1)) then
          local_result(idest)=local_result(idest)+ &
           volgrid*vel(D_DECL(i,j,k),prescomp) 
          local_result(idest+2)=local_result(idest+2)+volgrid
         endif
         idest=left_pressure_comp+2
         if (xsten(1,1).ge.probhix-VOFTOL*dx(1)) then
          local_result(idest)=local_result(idest)+ &
           volgrid*vel(D_DECL(i,j,k),prescomp) 
          local_result(idest+2)=local_result(idest+2)+volgrid
         endif

         do dir=1,SDIM
          idest=drag_sum_comp+dir
          local_result(idest)=local_result(idest)+ &
            drag(D_DECL(i,j,k),dir)
          idest=pdrag_sum_comp+dir
          local_result(idest)=local_result(idest)+ &
            drag(D_DECL(i,j,k),SDIM+dir)
          idest=torque_sum_comp+dir
          local_result(idest)=local_result(idest)+ &
            drag(D_DECL(i,j,k),2*SDIM+dir)
          idest=ptorque_sum_comp+dir
          local_result(idest)=local_result(idest)+ &
            drag(D_DECL(i,j,k),3*SDIM+dir)
         enddo ! dir
         idest=step_perim_sum_comp+1
         local_result(idest)=local_result(idest)+ &
           drag(D_DECL(i,j,k),4*SDIM+1)

         do im=1,nmat

          idest=minden_sum_comp+2*(im-1)+1
          isrc=num_state_material*(im-1)+1
          if (local_result(idest).gt.den(D_DECL(i,j,k),isrc)) then
           local_result(idest)=den(D_DECL(i,j,k),isrc)
          endif
          idest=idest+1
          isrc=isrc+1

          if (local_result(idest).gt.den(D_DECL(i,j,k),isrc)) then
           local_result(idest)=den(D_DECL(i,j,k),isrc)
          endif

          idest=maxden_sum_comp+2*(im-1)+1
          isrc=num_state_material*(im-1)+1
          if (local_result(idest).lt.den(D_DECL(i,j,k),isrc)) then
           local_result(idest)=den(D_DECL(i,j,k),isrc)
          endif
          idest=idest+1
          isrc=isrc+1

          if (local_result(idest).lt.den(D_DECL(i,j,k),isrc)) then
           local_result(idest)=den(D_DECL(i,j,k),isrc)
          endif

         enddo  ! im=1..nmat

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
           print *,"dir invalid summass 3"
           stop
          endif
          xcen=xsten(0,dir)

          do iside=0,1

           do im=1,nmat

            if (iside.eq.0) then
             xbottom=xsten(-2,dir)
             xtop=xcen
             ls_below=lsfab(D_DECL(i-ii,j-jj,k-kk),im)
             ls_above=lsfab(D_DECL(i,j,k),im)
            else if (iside.eq.1) then
             xbottom=xcen
             xtop=xsten(2,dir)
             ls_below=lsfab(D_DECL(i,j,k),im)
             ls_above=lsfab(D_DECL(i+ii,j+jj,k+kk),im)
            else
             print *,"iside invalid"
             stop
            endif

            if ((ls_below*ls_above.le.zero).and. &
                (abs(ls_below)+abs(ls_above).ge.LSTOL*dx(1))) then
 
             if (ls_below.eq.zero) then
              xcrit=xbottom
             else if (ls_above.eq.zero) then
              xcrit=xtop
             else  ! phi-phi0=slope(x-x0) x=-phi0/slope+x0
              xcrit=xbottom-ls_below*(xtop-xbottom)/ &
                 (ls_above-ls_below)
             endif

             do dir3=1,SDIM
              cen_material(dir3)=xsten(0,dir3)
             enddo
             cen_material(dir)=xcrit

             idest=minint_sum_comp+3*(im-1)+dir
             if (xcrit.lt.local_result(idest)) then
              local_result(idest)=xcrit
             endif
             idest=maxint_sum_comp+3*(im-1)+dir
             if (xcrit.gt.local_result(idest)) then
              local_result(idest)=xcrit
             endif
         
             if (slice_dir+1.eq.dir) then
              in_slice=1
              do dir3=1,SDIM
               if (dir3.ne.dir) then
                if (xslice(dir3).lt.xsten(-1,dir3)-VOFTOL*dx(dir3)) then
                 in_slice=0
                else if (xslice(dir3).gt.xsten(1,dir3)+VOFTOL*dx(dir3)) then    
                 in_slice=0
                endif
               endif ! dir3<>dir ?
              enddo ! dir3
              if (in_slice.eq.1) then
               idest=minint_slice+im
               if (xcrit.lt.local_result(idest)) then
                local_result(idest)=xcrit
               endif
               idest=maxint_slice+im
               if (xcrit.gt.local_result(idest)) then
                local_result(idest)=xcrit
               endif
              else if (in_slice.eq.0) then
               ! do nothing
              else
               print *,"in_slice invalid"
               stop
              endif
             endif ! slice_dir+1 == dir ?
             
             if ((dir.eq.SDIM).and.(im.eq.1).and.(i.eq.0)) then
              idest=xnot_amp_sum_comp+1
              if (xcrit.gt.local_result(idest)) then
               local_result(idest)=xcrit
              endif
             endif

             if (isweep.eq.0) then
              ! do nothing
             else if (isweep.eq.1) then

              do dir2=1,SDIM
               xboundary(dir2)=xsten(0,dir2)
               cengrid(dir2)=resultALL(cen_sum_comp+dir2+3*(im-1))
              enddo
              xboundary(dir)=xcrit
              distbound=zero
              do dir2=1,SDIM
               distbound=distbound+(xboundary(dir2)-cengrid(dir2))**2
              enddo
              distbound=sqrt(distbound)
              idest=mincen_sum_comp+im
              if (distbound.lt.local_result(idest)) then
               local_result(idest)=distbound
              endif
              idest=maxcen_sum_comp+im
              if (distbound.gt.local_result(idest)) then
               local_result(idest)=distbound
              endif 
 
             else
              print *,"isweep invalid"
              stop
             endif

            endif  ! change of sign
           enddo ! im
          enddo ! iside
         enddo ! dir

        endif  ! mask>0

       enddo
       enddo
       enddo  ! i,j,k

        ! cut_flag is a parameter
       if (cut_flag.eq.1) then

        if (SDIM.eq.2) then

          ! r=h(z)
         if ((dirx.eq.SDIM-1).and.(diry.eq.0)) then
        
          k=0 
          do j = growlo(2),growhi(2)+1
          do i = growlo(1),growhi(1)
           if ((mask(D_DECL(i,j,k)).gt.zero).and. &
               (mask(D_DECL(i,j-1,k)).gt.zero).and. &
               (mask(D_DECL(i+1,j,k)).gt.zero).and. &
               (mask(D_DECL(i+1,j-1,k)).gt.zero).and. &
               (mask(D_DECL(i-1,j,k)).gt.zero).and. &
               (mask(D_DECL(i-1,j-1,k)).gt.zero)) then
            dir=diry+1  !diry=0 r=h(z)  dir=1
            im=1

            do iside=0,1
             dirMAC=2
             call gridstenMAC(xstenMAC,xlo,i,j,k,fablo,bfact,dx,nhalf,dirMAC)

             xcen=xstenMAC(0,dir) ! diry=0  dir=diry+1=1

             if (iside.eq.0) then
              xbottom=xstenMAC(-2,dir)
              xtop=xcen
              ls_below=half*(lsfab(D_DECL(i-1,j,k),im)+ &
                 lsfab(D_DECL(i-1,j-1,k),im))
              ls_above=half*(lsfab(D_DECL(i,j,k),im)+ &
                 lsfab(D_DECL(i,j-1,k),im))
             else if (iside.eq.1) then
              xbottom=xcen
              xtop=xstenMAC(2,dir)
              ls_below=half*(lsfab(D_DECL(i,j,k),im)+ &
                 lsfab(D_DECL(i,j-1,k),im))
              ls_above=half*(lsfab(D_DECL(i+1,j,k),im)+ &
                 lsfab(D_DECL(i+1,j-1,k),im))
             else
              print *,"iside invalid"
              stop
             endif

              ! coordinate along streamwise direction starts at 0
              ! dirx=sdim-1
              ! dirMAC=2
              ! r=h(z)
             ZZgrid=xstenMAC(0,dirx+1)-problo(dirx+1)
             dz_external=(probhi(dirx+1)-problo(dirx+1))/NN
              ! j_external * dz_external = ZZgrid
             j_external=NINT(ZZgrid/dz_external)  ! round to nearest whole int. 
             if ((j_external.lt.0).or.(j_external.gt.NN)) then
              print *,"j_external invalid"
              stop
             endif
             ZZ(j_external)=j_external*dz_external

             if ((ls_below*ls_above.le.zero).and. &
                 (abs(ls_below)+abs(ls_above).ge.LSTOL*dx(1))) then

              if (ls_below.eq.zero) then
               xcrit=xbottom
              else if (ls_above.eq.zero) then
               xcrit=xtop
              else  ! phi-phi0=slope(x-x0) x=-phi0/slope+x0
               xcrit=xbottom-ls_below*(xtop-xbottom)/ &
                  (ls_above-ls_below)
              endif

              FF(j_external)=xcrit
             endif ! interface found
            enddo ! iside
           endif  ! mask>0
          enddo
          enddo

         else
          print *,"dirx,diry option not supported"
          stop
         endif

        endif ! 2D

        if (SDIM.eq.3) then

          ! z=h(x)
         if ((dirx.eq.0).and.(diry.eq.SDIM-1)) then
        
          j=0 
          if ((j.ge.growlo(2)).and. &
              (j.le.growhi(2))) then

           do i = growlo(1),growhi(1)+1
           do k = growlo(3),growhi(3) 
            if ((mask(D_DECL(i,j,k)).gt.zero).and. &
                (mask(D_DECL(i-1,j,k)).gt.zero).and. &
                (mask(D_DECL(i,j,k+1)).gt.zero).and. &
                (mask(D_DECL(i-1,j,k+1)).gt.zero).and. &
                (mask(D_DECL(i,j,k-1)).gt.zero).and. &
                (mask(D_DECL(i-1,j,k-1)).gt.zero)) then
             dir=diry+1  ! diry=sdim-1  dir=sdim  z=h(x)
             im=1  ! check liquid height
             vofcomp=(im-1)*ngeom_recon+1

             do iside=0,1
              dirMAC=1
              call gridstenMAC(xstenMAC,xlo,i,j,k,fablo,bfact,dx,nhalf,dirMAC)

               ! diry=sdim-1  dir=sdim  z=h(x)
              xcen=xstenMAC(0,dir)

              if (iside.eq.0) then
               xbottom=xstenMAC(-2,dir)
               xtop=xcen
               ls_below=half*(lsfab(D_DECL(i-1,j,k-1),im)+ &
                 lsfab(D_DECL(i,j,k-1),im))
               ls_above=half*(lsfab(D_DECL(i-1,j,k),im)+ &
                 lsfab(D_DECL(i,j,k),im))
               vof_below=half*(slopes(D_DECL(i-1,j,k-1),vofcomp)+ &
                 slopes(D_DECL(i,j,k-1),vofcomp))
               vof_above=half*(slopes(D_DECL(i-1,j,k),vofcomp)+ &
                 slopes(D_DECL(i,j,k),vofcomp))
              else if (iside.eq.1) then
               xbottom=xcen
               xtop=xstenMAC(2,dir)
               ls_below=half*(lsfab(D_DECL(i-1,j,k),im)+ &
                  lsfab(D_DECL(i,j,k),im))
               ls_above=half*(lsfab(D_DECL(i-1,j,k+1),im)+ &
                  lsfab(D_DECL(i,j,k+1),im))
               vof_below=half*(slopes(D_DECL(i-1,j,k),vofcomp)+ &
                 slopes(D_DECL(i,j,k),vofcomp))
               vof_above=half*(slopes(D_DECL(i-1,j,k+1),vofcomp)+ &
                 slopes(D_DECL(i,j,k+1),vofcomp))
              else
               print *,"iside invalid"
               stop
              endif

               ! dirx=0 diry=sdim-1 dir=sdim  z=h(x) 
               ! coordinate along streamwise direction starts at 0
              ZZgrid=xstenMAC(0,dirx+1)-problo(dirx+1)
              dz_external=(probhi(dirx+1)-problo(dirx+1))/NN
               ! j_external * dz_external = ZZgrid
              j_external=NINT(ZZgrid/dz_external)  ! round to nearest whole int. 
              if ((j_external.lt.0).or.(j_external.gt.NN)) then
               print *,"j_external invalid"
               stop
              endif
              ZZ(j_external)=j_external*dz_external

              use_vof_height=1

               ! diry=sdim-1  dir=sdim  z=h(x)
              if (use_vof_height.eq.0) then

               if ((ls_below*ls_above.le.zero).and. &
                   (abs(ls_below)+abs(ls_above).ge.LSTOL*dx(1))) then

                if (ls_below.eq.zero) then
                 xcrit=xbottom
                else if (ls_above.eq.zero) then
                 xcrit=xtop
                else  ! phi-phi0=slope(x-x0) x=-phi0/slope+x0
                 xcrit=xbottom-ls_below*(xtop-xbottom)/ &
                   (ls_above-ls_below)
                endif
 
                FF(j_external)=xcrit
               endif ! interface found

              else if (use_vof_height.eq.1) then
               vof_face=half*(vof_below+vof_above)
               if ((vof_face.ge.one/four).and. &
                   (vof_face.le.three/four)) then
                if (vof_below.gt.vof_above) then
                 xcrit=xstenMAC(-3,dir)+ &
                  (xstenMAC(1,dir)-xstenMAC(-3,dir))*vof_face
                else
                 xcrit=xstenMAC(3,dir)- &
                  (xstenMAC(3,dir)-xstenMAC(-1,dir))*vof_face
                endif
                FF(j_external)=xcrit
               endif
              else
               print *,"use_vof_height invalid"
               stop
              endif

             enddo ! iside
            endif  ! mask>0
           enddo
           enddo ! i,j

          endif ! j=0 in grid.

         else
          print *,"dirx,diry option not supported"
          stop
         endif

        endif ! 3D

       else if (cut_flag.eq.0) then
        ! do nothing
       else
        print *,"cut_flag invalid"
        stop
       endif

       return
       end subroutine FORT_SUMMASS



       !! FORT_FABCOM: fabz = fabx + beta * faby
       !! Added by Alan Kuhnle, 6-7-10
       subroutine FORT_FABCOM( &
        fabx,DIMS(fabx), &
        faby, DIMS(faby), &
        dotmask, &
        DIMS(dotmask), &
        mask, DIMS(mask), &
        fabz, DIMS(fabz), &
        beta, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        nsolve, &
        nsolveMM, &
        num_materials_face)

       use global_utility_module
       use probcommon_module

       IMPLICIT NONE

       INTEGER_T nsolve
       INTEGER_T nsolveMM
       INTEGER_T num_materials_face
       INTEGER_T tilelo(SDIM),tilehi(SDIM)
       INTEGER_T fablo(SDIM),fabhi(SDIM)
       INTEGER_T growlo(3),growhi(3)
       INTEGER_T bfact
       INTEGER_T nc
       INTEGER_T DIMDEC(fabx)
       INTEGER_T DIMDEC(faby)
       INTEGER_T DIMDEC(fabz)
       INTEGER_T DIMDEC(dotmask)
       INTEGER_T DIMDEC(mask)
       REAL_T  beta
       REAL_T  fabx(DIMV(fabx),nsolveMM)
       REAL_T  faby(DIMV(faby),nsolveMM)
       REAL_T  fabz(DIMV(fabz),nsolveMM)
       REAL_T  dotmask(DIMV(dotmask),num_materials_face)
       REAL_T  mask(DIMV(mask))

       INTEGER_T local_mask
       INTEGER_T local_dotmask
       INTEGER_T nc_mask

       INTEGER_T i,j,k

       if (bfact.lt.1) then
        print *,"bfact invalid157"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid"
        stop
       endif
       if (num_materials_vel.ne.1) then
        print *,"num_materials_vel invalid"
        stop
       endif
       if ((num_materials_face.ne.1).and. &
           (num_materials_face.ne.num_materials)) then
        print *,"num_materials_face invalid"
        stop
       endif
       if (nsolveMM.ne.nsolve*num_materials_face) then
        print *,"nsolveMM invalid"
        stop
       endif

       call checkbound(fablo,fabhi,DIMS(dotmask),0,-1,414) 
       call checkbound(fablo,fabhi,DIMS(mask),0,-1,414) 
       call checkbound(fablo,fabhi,DIMS(fabx),0,-1,415) 
       call checkbound(fablo,fabhi,DIMS(faby),0,-1,416) 
       call checkbound(fablo,fabhi,DIMS(fabz),0,-1,417) 
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do nc=1,nsolveMM

        if (num_materials_face.eq.1) then
         nc_mask=1
        else if (num_materials_face.eq.num_materials) then
         nc_mask=nc
        else
         print *,"num_materials_face invalid"
         stop
        endif

        if ((nc_mask.ge.1).and.(nc_mask.le.num_materials_face)) then

         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          local_mask=NINT(mask(D_DECL(i,j,k)))
          local_dotmask=NINT(dotmask(D_DECL(i,j,k),nc_mask))

          if (num_materials_face.eq.1) then
           if (local_dotmask.eq.1) then
            ! do nothing
           else
            print *,"local_dotmask invalid"
            stop
           endif
          else if (num_materials_face.eq.num_materials) then
           if ((local_dotmask.eq.1).or.(local_dotmask.eq.0)) then
            ! do nothing
           else
            print *,"local_dotmask invalid"
            stop
           endif
          else
           print *,"num_materials_face invalid"
           stop
          endif

          if (local_dotmask.eq.1) then

           if (local_mask.eq.1) then

            fabz(D_DECL(i,j,k),nc) = fabx(D_DECL(i,j,k),nc) +  &
             beta*faby(D_DECL(i,j,k),nc)

           else if (local_mask.eq.0) then
            ! do nothing
           else
            print *,"mask invalid"
            stop
           endif

          else if (local_dotmask.eq.0) then
           ! do nothing
          else
           print *,"local_dotmask invalid"
           stop
          endif

         enddo !k
         enddo !j
         enddo !i

        else
         print *,"nc_mask invalid"
         stop
        endif

       enddo !nc=1..nsolveMM

       return
       end subroutine FORT_FABCOM

       subroutine FORT_DIAGINV( &
        singular_possible, &
        offdiag_nonsing_level, &
        diag_regularization, &
        diagnonsing, &
        DIMS(diagnonsing), &
        diagsing, &
        DIMS(diagsing), &
        resid, DIMS(resid), &
        xnew, DIMS(xnew), &
        xold, DIMS(xold), &
        mask, DIMS(mask), &
        tilelo,tilehi, &
        fablo,fabhi,bfact )

       use global_utility_module

       IMPLICIT NONE

       INTEGER_T, intent(in) :: singular_possible
       REAL_T, intent(in) :: offdiag_nonsing_level
       REAL_T, intent(in) :: diag_regularization
       REAL_T :: local_diag
       INTEGER_T, intent(in) :: bfact
       INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
       INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
       INTEGER_T :: growlo(3),growhi(3)
       INTEGER_T, intent(in) :: DIMDEC(diagnonsing)
       INTEGER_T, intent(in) :: DIMDEC(diagsing)
       INTEGER_T, intent(in) :: DIMDEC(resid)
       INTEGER_T, intent(in) :: DIMDEC(xnew)
       INTEGER_T, intent(in) :: DIMDEC(xold)
       INTEGER_T, intent(in) :: DIMDEC(mask)
       REAL_T, intent(in) ::  diagnonsing(DIMV(diagnonsing))
       REAL_T, intent(in) ::  diagsing(DIMV(diagsing))
       REAL_T, intent(in) ::  resid(DIMV(resid))
       REAL_T, intent(out) :: xnew(DIMV(xnew))
       REAL_T, intent(in) ::  xold(DIMV(xold))
       REAL_T, intent(in) ::  mask(DIMV(mask))

       INTEGER_T i,j,k

       if (bfact.lt.1) then
        print *,"bfact invalid158"
        stop
       endif
       FIX ME
       if (offdiag_nonsing_level.gt.zero) then
        ! do nothing
       else
        print *,"offdiag_nonsing_level<=0"
        stop
       endif
       FIX ME
       if ((diag_regularization.gt.zero).and. &
           (diag_regularization.le.1.0D-3)) then
        ! do nothing
       else 
        print *,"diag_regularization invalid"
        stop
       endif
       if ((singular_possible.ge.0).and. &
           (singular_possible.le.1)) then
        ! do nothing
       else
        print *,"singular_possible invalid"
        stop
       endif

       call checkbound(fablo,fabhi,DIMS(mask),1,-1,414) 
       call checkbound(fablo,fabhi,DIMS(diagnonsing),0,-1,415) 
       call checkbound(fablo,fabhi,DIMS(diagsing),0,-1,415) 
       call checkbound(fablo,fabhi,DIMS(xnew),1,-1,416) 
       call checkbound(fablo,fabhi,DIMS(xold),1,-1,417) 
       call checkbound(fablo,fabhi,DIMS(resid),0,-1,417) 

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

        if (mask(D_DECL(i,j,k)).eq.one) then

                FIX ME
         if (diagsing(D_DECL(i,j,k)).eq.zero) then
           ! offdiag_nonsing_level ~ coeff_max * area / dx
          if (singular_possible.eq.1) then
           local_diag=diag_regularization*offdiag_nonsing_level
          else
           print *,"diag or singular_possible invalid 1"
           stop
          endif
         else if (diagsing(D_DECL(i,j,k)).gt.zero) then
          local_diag=diagsing(D_DECL(i,j,k))
         else
          print *,"diag invalid 2"
          stop
         endif
THIS IS OK - diagnonsing should be ok
         local_diag=diagnonsing(D_DECL(i,j,k))

         xnew(D_DECL(i,j,k))=xold(D_DECL(i,j,k))+ &
           resid(D_DECL(i,j,k))/local_diag
        else if (mask(D_DECL(i,j,k)).eq.zero) then
         xnew(D_DECL(i,j,k))=xold(D_DECL(i,j,k))
        else 
         print *,"mask invalid"
         stop
        endif

       enddo
       enddo
       enddo

       return
       end subroutine FORT_DIAGINV


       subroutine FORT_SUMDOT(mass1,  &
        rho,DIMS(rho), &
        rho2,DIMS(rho2), &
        dotmask, &
        DIMS(dotmask), &
        mask,DIMS(mask), &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        debug_dot_product, &
        levelno,gridno, &
        nsolve, &
        nsolveMM, &
        num_materials_face)
    
       use global_utility_module
       use probf90_module
 
       IMPLICIT NONE

       INTEGER_T levelno,gridno
       INTEGER_T nsolve
       INTEGER_T nsolveMM
       INTEGER_T num_materials_face
       INTEGER_T debug_dot_product
       INTEGER_T imax,jmax,kmax
       REAL_T dotmax
       INTEGER_T tilelo(SDIM),tilehi(SDIM)
       INTEGER_T fablo(SDIM),fabhi(SDIM)
       INTEGER_T growlo(3),growhi(3)
       INTEGER_T bfact
       INTEGER_T DIMDEC(rho)
       INTEGER_T DIMDEC(rho2)
       INTEGER_T DIMDEC(dotmask)
       INTEGER_T DIMDEC(mask)
       REAL_T  mass1,dm
       REAL_T  rho(DIMV(rho),nsolveMM)
       REAL_T  rho2(DIMV(rho2),nsolveMM)
       REAL_T  dotmask(DIMV(dotmask),num_materials_face)
       REAL_T  mask(DIMV(mask))
       INTEGER_T local_mask
       INTEGER_T local_dotmask
       INTEGER_T nc_mask

       INTEGER_T i,j,k,nc


       if ((levelno.lt.0).or.(gridno.lt.0)) then
        print *,"level or grid invalid"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid"
        stop
       endif
       if (num_materials_vel.ne.1) then
        print *,"num_materials_vel invalid"
        stop
       endif
       if ((num_materials_face.ne.1).and. &
           (num_materials_face.ne.num_materials)) then
        print *,"num_materials_face invalid"
        stop
       endif

       if (nsolveMM.ne.nsolve*num_materials_face) then
        print *,"nsolveMM invalid"
        stop
       endif

       call checkbound(fablo,fabhi,DIMS(dotmask),0,-1,414) 
       call checkbound(fablo,fabhi,DIMS(mask),0,-1,414) 
       call checkbound(fablo,fabhi,DIMS(rho),0,-1,415) 
       call checkbound(fablo,fabhi,DIMS(rho2),0,-1,416) 
 
       mass1=zero

       imax=0
       jmax=0
       kmax=0
       dotmax=zero
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do nc=1,nsolveMM

        if (num_materials_face.eq.1) then
         nc_mask=1
        else if (num_materials_face.eq.num_materials) then
         nc_mask=nc
        else
         print *,"num_materials_face invalid"
         stop
        endif

        if ((nc_mask.ge.1).and.(nc_mask.le.num_materials_face)) then

         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          local_mask=NINT(mask(D_DECL(i,j,k)))
          local_dotmask=NINT(dotmask(D_DECL(i,j,k),nc_mask))

          if (num_materials_face.eq.1) then
           if (local_dotmask.eq.1) then
            ! do nothing
           else
            print *,"local_dotmask invalid"
            stop
           endif
          else if (num_materials_face.eq.num_materials) then
           if ((local_dotmask.eq.1).or.(local_dotmask.eq.0)) then
            ! do nothing
           else
            print *,"local_dotmask invalid"
            stop
           endif
          else
           print *,"num_materials_face invalid"
           stop
          endif

          if (local_dotmask.eq.1) then

           if (local_mask.eq.1) then

            dm=rho(D_DECL(i,j,k),nc)*rho2(D_DECL(i,j,k),nc)

            if (debug_dot_product.eq.1) then
             if (abs(dm).gt.dotmax) then
              dotmax=abs(dm)
              imax=i
              jmax=j
              kmax=k
             endif
            else if (debug_dot_product.ne.0) then
             print *,"debug dot prod invalid"
             stop
            endif

            mass1=mass1+dm
           else if (local_mask.eq.0) then
            ! do nothing
           else 
            print *,"mask invalid"
            stop
           endif

          else if (local_dotmask.eq.0) then
           ! do nothing
          else
           print *,"local_dotmask invalid"
           stop
          endif

         enddo ! k
         enddo ! j
         enddo ! i

        else
         print *,"nc_mask invalid"
         stop
        endif

       enddo ! nc=1..nsolveMM

       if (debug_dot_product.eq.1) then
        print *,"debug dot: level,grid,i,j,k,max ",levelno,gridno, &
         imax,jmax,kmax,dotmax
       else if (debug_dot_product.ne.0) then
        print *,"debug dot prod invalid"
        stop
       endif

       return
       end subroutine FORT_SUMDOT

! coriolis force:
! R''_space=R''_earth + 2 (omega cross v)+omega cross (omega cross R)+
! w' cross R.
! coriolis term=-2 Omega cross v
! suppose Omega=(0 0 omega)  (z axis)
! coriolis force=-2 ( -omega v  omega u  0 ) 
! centrifugal force:
! F=mass r w^2   w units: radians/s   mass r w^2 = kg m/s^2
!
! gravity_normalized>0 means that gravity is directed downwards.
! if invert_gravity==1, then gravity_normalized<0 (pointing upwards)

      subroutine FORT_INITPOTENTIAL( &
       nmat, &
       ngrow, &
       override_density, &
       presden,DIMS(presden), &
       state,DIMS(state), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       presbc_arr, &
       dombcpres, &
       domlo,domhi, &
       xlo,dx,dt, &
       gravity_normalized, &
       gravity_dir_parm, &
       angular_velocity, &
       isweep)

      use global_utility_module
      use probf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: ngrow
      INTEGER_T, intent(in) :: override_density,level
      INTEGER_T, intent(in) :: gravity_dir_parm,isweep
      REAL_T, intent(in) :: dt
      REAL_T, intent(in) :: gravity_normalized  
      REAL_T, intent(in) :: angular_velocity
      INTEGER_T, intent(in) :: DIMDEC(presden)
      INTEGER_T, intent(in) :: DIMDEC(state)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T, intent(in) :: domlo(SDIM),domhi(SDIM)
       ! pressure bc at exterior walls
      INTEGER_T, intent(in) :: dombcpres(SDIM,2) 
       ! presbc at int or ext.
      INTEGER_T, intent(in) :: presbc_arr(SDIM,2,num_materials_vel) 

       !HYDROSTATIC_PRESSURE,HYDROSTATIC_DENSITY
      REAL_T, intent(inout) :: presden(DIMV(presden),2) 
      REAL_T, intent(in) :: state(DIMV(state),nmat*num_state_material) 
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      REAL_T xcell(SDIM)
      REAL_T xwall_normal
      REAL_T xsten(-1:1,SDIM)
      REAL_T den_cell,pres_cell

      INTEGER_T i,j,k,im,icomp,nhalf,dir2
      INTEGER_T ii,jj,kk,dir,side,inside
      INTEGER_T bctypepres,local_bctype,exteriorbc
      INTEGER_T for_hydro
      REAL_T liquid_temp

      nhalf=1

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid159"
       stop
      endif
      if (ngrow.ne.1) then
       print *,"ngrow invalid"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      if ((gravity_dir_parm.lt.1).or.(gravity_dir_parm.gt.SDIM)) then
       print *,"gravity dir invalid initpotential"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(presden),ngrow,-1,42)
      call checkbound(fablo,fabhi,DIMS(state),ngrow,-1,42)
     
      if (isweep.eq.0) then
 
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        call gridsten_level(xsten,i,j,k,level,nhalf)
        do dir2=1,SDIM
         xcell(dir2)=xsten(0,dir2)
        enddo 

         ! includes centrifugal force but not "coriolis force"
         ! p=dt( -|g| z + (1/2)Omega^2 r^2 )
        for_hydro=1
        im=1
        icomp=(im-1)*num_state_material+2
        liquid_temp=state(D_DECL(i,j,k),icomp) 
        call general_hydrostatic_pressure_density( &
          override_density, &
          xcell, &
          gravity_normalized, &
          gravity_dir_parm, &
          angular_velocity, &
          dt,den_cell,pres_cell, &
          for_hydro,liquid_temp)
        presden(D_DECL(i,j,k),1)=pres_cell
        presden(D_DECL(i,j,k),2)=den_cell
       enddo     
       enddo     
       enddo     

      else if (isweep.eq.1) then

       do dir=1,SDIM
       do side=1,2

        if (dir.eq.1) then
         if (side.eq.1) then
          xwall_normal=problox
         else if (side.eq.2) then
          xwall_normal=probhix
         else
          print *,"side invalid"
          stop
         endif
        else if (dir.eq.2) then
         if (side.eq.1) then
          xwall_normal=probloy
         else if (side.eq.2) then
          xwall_normal=probhiy
         else
          print *,"side invalid"
          stop
         endif
        else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
         if (side.eq.1) then
          xwall_normal=probloz
         else if (side.eq.2) then
          xwall_normal=probhiz
         else
          print *,"side invalid"
          stop
         endif
        else
         print *,"dir invalid"
         stop
        endif

        call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow)

         ! dombcpres=get_desc_lst()[State_Type].getBC(pcomp)
        bctypepres=dombcpres(dir,side)  
         ! presbc_arr=getBCArray(State_Type,gridno,pcomp,1)
        local_bctype=presbc_arr(dir,side,1) 

        exteriorbc=0

        if (side.eq.1) then
         inside=fablo(dir)
         growhi(dir)=fablo(dir)-1
         if (growlo(dir).lt.domlo(dir)) then
          exteriorbc=1
         endif
        else if (side.eq.2) then
         inside=fabhi(dir)
         growlo(dir)=fabhi(dir)+1
         if (growhi(dir).gt.domhi(dir)) then
          exteriorbc=1
         endif
        else
         print *,"side invalid"
         stop
        endif

        if (exteriorbc.eq.1) then

         do i=growlo(1),growhi(1)
         do j=growlo(2),growhi(2)
         do k=growlo(3),growhi(3)

          if (local_bctype.eq.bctypepres) then

            ! outflow wall
           if (local_bctype.eq.EXT_DIR) then

            call gridsten_level(xsten,i,j,k,level,nhalf)
            do dir2=1,SDIM
             xcell(dir2)=xsten(0,dir2)
            enddo 
            xcell(dir)=xwall_normal

            for_hydro=1
            im=1
            icomp=(im-1)*num_state_material+2
            liquid_temp=state(D_DECL(i,j,k),icomp) 
             ! p=dt( -|g| z + (1/2)Omega^2 r^2 )
            call general_hydrostatic_pressure_density( &
             override_density, &
             xcell, &
             gravity_normalized, &
             gravity_dir_parm, &
             angular_velocity, &
             dt,den_cell,pres_cell, &
             for_hydro,liquid_temp)

            ! periodic BC
           else if (local_bctype.eq.INT_DIR) then
            pres_cell=presden(D_DECL(i,j,k),1)
            den_cell=presden(D_DECL(i,j,k),2)

            ! symmetric, slip or noslip wall.
           else if ((local_bctype.eq.FOEXTRAP).or. &
                    (local_bctype.eq.REFLECT_EVEN)) then
            ii=i
            jj=j
            kk=k
            if (dir.eq.1) then
             ii=inside
            else if (dir.eq.2) then
             jj=inside
            else if ((dir.eq.3).and.(SDIM.eq.3)) then
             kk=inside
            else
             print *,"dir invalid initpotential"
             stop
            endif
            pres_cell=presden(D_DECL(ii,jj,kk),1)
            den_cell=presden(D_DECL(ii,jj,kk),2)
           else
            print *,"local_bctype invalid"
            stop
           endif

           presden(D_DECL(i,j,k),1)=pres_cell
           presden(D_DECL(i,j,k),2)=den_cell

          else 
           print *,"local_bctype should equal bctypepres outside domain"
           stop
          endif 
                  
         enddo
         enddo
         enddo

        else if (exteriorbc.eq.0) then

         if (local_bctype.eq.INT_DIR) then
          ! do nothing
         else
          print *,"local_bctype invalid"
          stop
         endif

        else
         print *,"exteriorbc invalid"
         stop
        endif

       enddo
       enddo !dir,side

      else
       print *,"isweep invalid"
       stop
      endif

      return
      end subroutine FORT_INITPOTENTIAL


      subroutine FORT_ADDGRAVITY( &
       denconst_gravity, &
       nsolveMM_FACE, &
       level, &
       finest_level, &
       facecut_index, &
       icefacecut_index, &
       vofface_index, &
       ncphys, &
       nmat, &
       nstate, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx,dir, &
       xface,DIMS(xface), &
       recon,DIMS(recon), &
       lsnew,DIMS(lsnew), &
       snew,DIMS(snew), &
       macnew,DIMS(macnew), &
       cellgrav,DIMS(cellgrav), &
       facegrav,DIMS(facegrav) )

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: nsolveMM_FACE
      INTEGER_T :: nsolveMM_FACE_test
      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: facecut_index
      INTEGER_T, intent(in) :: icefacecut_index
      INTEGER_T, intent(in) :: vofface_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nstate
      INTEGER_T, intent(in) :: dir
      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: denconst_gravity(nmat)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(recon)
      INTEGER_T, intent(in) :: DIMDEC(lsnew)
      INTEGER_T, intent(in) :: DIMDEC(snew)
      INTEGER_T, intent(in) :: DIMDEC(macnew)
      INTEGER_T, intent(in) :: DIMDEC(cellgrav)
      INTEGER_T, intent(in) :: DIMDEC(facegrav)

      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: recon(DIMV(recon),nmat*ngeom_recon) 
      REAL_T, intent(in) :: lsnew(DIMV(lsnew),nmat*(SDIM+1))
      REAL_T, intent(inout) :: snew(DIMV(snew),nstate)
      REAL_T, intent(inout) :: macnew(DIMV(macnew),nsolveMM_FACE)
      REAL_T, intent(in) :: cellgrav(DIMV(cellgrav),SDIM)
      REAL_T, intent(in) :: facegrav(DIMV(facegrav))
 
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      INTEGER_T iside

      INTEGER_T im
      INTEGER_T velcomp
      REAL_T grav_component
      REAL_T local_cut
      REAL_T local_macnewL
      REAL_T local_macnewR

      REAL_T vol_total,mass_total,volside,denface_gravity

      REAL_T xsten(-1:1,SDIM)
      INTEGER_T nhalf

      nhalf=1

      if (bfact.lt.1) then
       print *,"bfact invalid160"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid add gravity"
       stop
      endif
 
      if (nstate.ne.num_materials_vel*(SDIM+1)+ &
          nmat*(num_state_material+ngeom_raw)+1) then
       print *,"nstate invalid"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"num_materials invalid add gravity"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (facecut_index.ne.3) then
       print *,"facecut_index invalid"
       stop
      endif
      if (icefacecut_index.ne.4) then
       print *,"icefacecut_index invalid"
       stop
      endif
      if (ncphys.ne.vofface_index+2*nmat) then
       print *,"ncphys invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif
      nsolveMM_FACE_test=num_materials_vel
      if (nsolveMM_FACE_test.ne.nsolveMM_FACE) then
       print *,"nsolveMM_FACE invalid"
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
       print *,"dir out of range in addgravity"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(xface),0,dir,42)
      call checkbound(fablo,fabhi,DIMS(recon),1,-1,42)
      call checkbound(fablo,fabhi,DIMS(lsnew),1,-1,42)
      call checkbound(fablo,fabhi,DIMS(snew),1,-1,42)
      call checkbound(fablo,fabhi,DIMS(macnew),0,dir,42)
      call checkbound(fablo,fabhi,DIMS(cellgrav),0,-1,42)
      call checkbound(fablo,fabhi,DIMS(facegrav),0,dir,42)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir+1)

       vol_total=zero
       mass_total=zero

       do iside=0,1
       do im=1,nmat

        volside=xface(D_DECL(i,j,k),vofface_index+2*(im-1)+iside+1)
        vol_total=vol_total+volside
        if (denconst_gravity(im).ge.zero) then
         mass_total=mass_total+denconst_gravity(im)*volside
        else
         print *,"denconst_gravity invalid"
         stop
        endif
     
       enddo ! im=1..nmat
       enddo ! iside=0..1

       if (vol_total.eq.zero) then
        denface_gravity=one
       else if (vol_total.gt.zero) then
        denface_gravity=mass_total/vol_total
       else
        print *,"vol_total invalid"
        stop
       endif

       local_cut=xface(D_DECL(i,j,k),facecut_index+1)
       if ((local_cut.ge.zero).and.(local_cut.le.half)) then
        local_cut=zero
       else if ((local_cut.ge.half).and.(local_cut.le.one)) then
        ! do nothing
       else
        print *,"local_cut invalid"
        stop
       endif

       im=1
       local_macnewL=macnew(D_DECL(i,j,k),im)+ &
         denface_gravity*local_cut*facegrav(D_DECL(i,j,k))

       local_macnewR=local_macnewL

       if (levelrz.eq.0) then
        ! do nothing
       else if (levelrz.eq.1) then
        if (SDIM.ne.2) then
         print *,"rz in 2d only"
         stop
        endif
        if ((dir.eq.0).and. &
            (xsten(0,1).le.VOFTOL*dx(1))) then
         local_macnewL=zero
         local_macnewR=zero
        endif
       else if (levelrz.eq.3) then
        ! do nothing
       else
        print *,"levelrz invalid"
        stop 
       endif

       macnew(D_DECL(i,j,k),im)=local_macnewL

      enddo
      enddo
      enddo

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       local_cut=one
       do im=1,nmat
        if (is_prescribed(nmat,im).eq.1) then
         if (lsnew(D_DECL(i,j,k),im).ge.zero) then
          local_cut=zero
         else if (lsnew(D_DECL(i,j,k),im).lt.zero) then
          ! do nothing
         else
          print *,"lsnew bust"
          stop
         endif
        else if (is_prescribed(nmat,im).eq.0) then
         ! do nothing
        else
         print *,"is_prescribed(nmat,im) invalid"
         stop
        endif
       enddo ! im=1..nmat

       if (local_cut.eq.one) then
        grav_component=cellgrav(D_DECL(i,j,k),dir+1)

        vol_total=zero
        mass_total=zero

        do iside=0,1
        do im=1,nmat

         if (iside.eq.0) then
          volside=xface(D_DECL(i,j,k),vofface_index+2*(im-1)+2)
         else if (iside.eq.1) then
          volside=xface(D_DECL(i+ii,j+jj,k+kk),vofface_index+2*(im-1)+1)
         else
          print *,"iside invalid"
          stop
         endif
         vol_total=vol_total+volside
         if (denconst_gravity(im).ge.zero) then
          mass_total=mass_total+denconst_gravity(im)*volside
         else
          print *,"denconst_gravity invalid"
          stop
         endif

        enddo ! im=1..nmat
        enddo ! iside=0..1

        if (vol_total.gt.zero) then
         denface_gravity=mass_total/vol_total
        else
         print *,"vol_total invalid"
         stop
        endif

        velcomp=dir+1 
        snew(D_DECL(i,j,k),velcomp)= &
          snew(D_DECL(i,j,k),velcomp)+denface_gravity*grav_component

       else if (local_cut.eq.zero) then
        ! do nothing
       else
        print *,"local_cut invalid"
        stop
       endif 

      enddo
      enddo
      enddo

      return
      end subroutine FORT_ADDGRAVITY

       ! updatevelALL called from: NavierStokes::multiphase_project
       ! mac_update is called from: NavierStokes::updatevelALL
       ! correct_velocity called from: NavierStokes::multiphase_project
       ! correct_velocity called from: NavierStokes::relaxLEVEL
       ! correct_velocity called from: NavierStokes::residual_correction_form
       ! correct_velocity called from: NavierStokes::mac_update
       ! called from: NavierStokes::correct_velocity (NavierStokes3.cpp)
      subroutine FORT_FLUIDSOLIDCOR( &
       level, &
       finest_level, &
       velcomp, &
       nsolve, &
       nsolveMM, &
       nsolveMM_face, &
       facecut_index, &
       icefacecut_index, &
       ncphys, &
       project_option, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       presbc_in, &
       maskcov,DIMS(maskcov), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       alt_xface, &
       DIMS(alt_xface), &
       alt_yface, &
       DIMS(alt_yface), &
       alt_zface, &
       DIMS(alt_zface), &
       xgp,DIMS(xgp), &
       ygp,DIMS(ygp), &
       zgp,DIMS(zgp), &
       xsrc,DIMS(xsrc), &
       ysrc,DIMS(ysrc), &
       zsrc,DIMS(zsrc), &
       xdest,DIMS(xdest), &
       ydest,DIMS(ydest), &
       zdest,DIMS(zdest), &
       xlo,dx,cur_time,nmat)

      use global_utility_module
      use probcommon_module
      use probf90_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: level,finest_level
      INTEGER_T, intent(in) :: velcomp
      INTEGER_T :: veldir
      INTEGER_T, intent(in) :: nsolve
      INTEGER_T, intent(in) :: nsolveMM
      INTEGER_T, intent(in) :: nsolveMM_face
      INTEGER_T, intent(in) :: facecut_index
      INTEGER_T, intent(in) :: icefacecut_index
      INTEGER_T, intent(in) :: ncphys
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: project_option
      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(xface)
      INTEGER_T, intent(in) :: DIMDEC(yface)
      INTEGER_T, intent(in) :: DIMDEC(zface)
      INTEGER_T, intent(in) :: DIMDEC(alt_xface)
      INTEGER_T, intent(in) :: DIMDEC(alt_yface)
      INTEGER_T, intent(in) :: DIMDEC(alt_zface)
      INTEGER_T, intent(in) :: DIMDEC(xgp)
      INTEGER_T, intent(in) :: DIMDEC(ygp)
      INTEGER_T, intent(in) :: DIMDEC(zgp)
      INTEGER_T, intent(in) :: DIMDEC(xsrc)
      INTEGER_T, intent(in) :: DIMDEC(ysrc)
      INTEGER_T, intent(in) :: DIMDEC(zsrc)
      INTEGER_T, intent(in) :: DIMDEC(xdest)
      INTEGER_T, intent(in) :: DIMDEC(ydest)
      INTEGER_T, intent(in) :: DIMDEC(zdest)

      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      INTEGER_T :: dir
      INTEGER_T, intent(in) :: presbc_in(SDIM,2,nsolveMM)

      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: xface(DIMV(xface),ncphys)
      REAL_T, intent(in) :: yface(DIMV(yface),ncphys)
      REAL_T, intent(in) :: zface(DIMV(zface),ncphys)
      REAL_T, intent(in) :: alt_xface(DIMV(alt_xface))
      REAL_T, intent(in) :: alt_yface(DIMV(alt_yface))
      REAL_T, intent(in) :: alt_zface(DIMV(alt_zface))

      REAL_T, intent(in) :: xgp(DIMV(xgp))
      REAL_T, intent(in) :: ygp(DIMV(ygp))
      REAL_T, intent(in) :: zgp(DIMV(zgp))

      REAL_T, intent(in) :: xsrc(DIMV(xsrc))
      REAL_T, intent(in) :: ysrc(DIMV(ysrc))
      REAL_T, intent(in) :: zsrc(DIMV(zsrc))

      REAL_T, intent(out) :: xdest(DIMV(xdest))
      REAL_T, intent(out) :: ydest(DIMV(ydest))
      REAL_T, intent(out) :: zdest(DIMV(zdest))

      REAL_T, intent(in) :: xlo(SDIM),dx(SDIM)
      REAL_T, intent(in) :: cur_time
 
      INTEGER_T i,j,k
      INTEGER_T ii,jj,kk
      REAL_T local_dest,local_src,local_gp
      REAL_T AL
      REAL_T AL_ice
      REAL_T alt_AL_ice
      REAL_T local_AL_ice
      REAL_T local_dd,local_visc_coef,cc_group,local_dd_group
      INTEGER_T local_constant_viscosity
      INTEGER_T icrit,side,bccrit
      INTEGER_T im_vel,bc_comp
      REAL_T local_wt(nsolve)
      INTEGER_T caller_id
      REAL_T maskleft,maskright
      INTEGER_T maskface

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if ((level.ge.0).and.(level.le.finest_level)) then
       ! do nothing
      else
       print *,"level invalid"
       stop
      endif

       ! indexes start at 0
      if (facecut_index.ne.3) then
       print *,"facecut_index bust"
       stop
      endif
      if (icefacecut_index.ne.4) then
       print *,"icefacecut_index bust"
       stop
      endif
      if (ncphys.lt.9) then
       print *,"ncphys invalid"
       stop
      endif
      if ((project_option.eq.0).or. &
          (project_option.eq.1).or. &
          (project_option.eq.10).or. &
          (project_option.eq.13).or. & !FSI_material_exists 1st project
          (project_option.eq.11).or. & !FSI_material_exists 2nd project
          (project_option.eq.12).or. & !pressure extrap
          (project_option.eq.2).or. & ! thermal diffusion
          (project_option.eq.3)) then ! viscosity
       ! do nothing
      else if ((project_option.ge.100).and. &
               (project_option.lt.100+num_species_var)) then
       ! do nothing
      else
       print *,"project_option invalid fluid solid cor"
       stop
      endif 
      if (bfact.lt.1) then
       print *,"bfact invalid161"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,268)

      call checkbound(fablo,fabhi,DIMS(xface),0,0,268)
      call checkbound(fablo,fabhi,DIMS(yface),0,1,269)
      call checkbound(fablo,fabhi,DIMS(zface),0,SDIM-1,270)
      call checkbound(fablo,fabhi,DIMS(alt_xface),0,0,244)
      call checkbound(fablo,fabhi,DIMS(alt_yface),0,1,244)
      call checkbound(fablo,fabhi,DIMS(alt_zface),0,SDIM-1,244)

      call checkbound(fablo,fabhi,DIMS(xgp),0,0,2333)
      call checkbound(fablo,fabhi,DIMS(ygp),0,1,2334)
      call checkbound(fablo,fabhi,DIMS(zgp),0,SDIM-1,2335)

      call checkbound(fablo,fabhi,DIMS(xsrc),0,0,2336)
      call checkbound(fablo,fabhi,DIMS(ysrc),0,1,2337)
      call checkbound(fablo,fabhi,DIMS(zsrc),0,SDIM-1,2338)

      call checkbound(fablo,fabhi,DIMS(xdest),0,0,2339)
      call checkbound(fablo,fabhi,DIMS(ydest),0,1,2340)
      call checkbound(fablo,fabhi,DIMS(zdest),0,SDIM-1,2341)

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
        print *,"dir out of range in FLUIDSOLIDCOR"
        stop
       endif

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)

       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)

         maskleft=maskcov(D_DECL(i-ii,j-jj,k-kk))
         maskright=maskcov(D_DECL(i,j,k))
         if ((maskleft.eq.one).and.(maskright.eq.one)) then
          maskface=1
         else if ((maskleft.eq.one).and.(maskright.eq.zero)) then
          maskface=0
         else if ((maskleft.eq.zero).and.(maskright.eq.one)) then
          maskface=0
         else if ((maskleft.eq.zero).and.(maskright.eq.zero)) then
          maskface=0
         else
          print *,"maskleft or maskright invalid"
          stop
         endif

         if (dir.eq.0) then
          local_dest=xdest(D_DECL(i,j,k))
          local_src=xsrc(D_DECL(i,j,k))
          local_gp=xgp(D_DECL(i,j,k))
          AL=xface(D_DECL(i,j,k),facecut_index+1)
          AL_ice=xface(D_DECL(i,j,k),icefacecut_index+1)
          alt_AL_ice=alt_xface(D_DECL(i,j,k))
          icrit=i
         else if (dir.eq.1) then
          local_dest=ydest(D_DECL(i,j,k))
          local_src=ysrc(D_DECL(i,j,k))
          local_gp=ygp(D_DECL(i,j,k))
          AL=yface(D_DECL(i,j,k),facecut_index+1)
          AL_ice=yface(D_DECL(i,j,k),icefacecut_index+1)
          alt_AL_ice=alt_yface(D_DECL(i,j,k))
          icrit=j
         else if ((dir.eq.2).and.(SDIM.eq.3)) then
          local_dest=zdest(D_DECL(i,j,k))
          local_src=zsrc(D_DECL(i,j,k))
          local_gp=zgp(D_DECL(i,j,k))
          AL=zface(D_DECL(i,j,k),facecut_index+1)
          AL_ice=zface(D_DECL(i,j,k),icefacecut_index+1)
          alt_AL_ice=alt_zface(D_DECL(i,j,k))
          icrit=k
         else
          print *,"dir invalid fluid solid cor"
          stop
         endif

         if (project_option.eq.13) then
          local_AL_ice=alt_AL_ice
         else
          local_AL_ice=AL_ice
         endif

         if ((project_option.eq.0).or. &
             (project_option.eq.1).or. &
             (project_option.eq.10).or. &
             (project_option.eq.11).or. & !FSI_material_exists 2nd project
             (project_option.eq.13).or. & !FSI_material_exists 1st project
             (project_option.eq.12)) then !pressure extrap.
          if ((nsolve.eq.1).and.(nsolveMM.eq.1).and.(velcomp.eq.0)) then
           veldir=1
           im_vel=1
           bc_comp=1
          else
           print *,"nsolve,nsolveMM, or velcomp invalid"
           stop
          endif
         else if ((project_option.eq.2).or. & ! thermal diffusion
                  ((project_option.ge.100).and. &
                   (project_option.lt.100+num_species_var))) then
          if ((nsolve.eq.1).and. &
              ((nsolveMM.eq.1).or.(nsolveMM.eq.nmat)).and. &
              (velcomp.ge.0).and. &
              (velcomp.lt.nsolveMM_face)) then
           veldir=1
           if (nsolveMM.eq.1) then
            im_vel=1
           else if (nsolveMM.eq.nmat) then
            if (velcomp.lt.nmat) then
             im_vel=velcomp+1
            else if ((velcomp.ge.nmat).and.(velcomp.lt.nsolveMM_face)) then
             im_vel=velcomp-nmat+1
            else
             print *,"velcomp invalid"
             stop
            endif
           else
            print *,"nsolveMM invalid"
            stop
           endif
           bc_comp=im_vel
          else
           print *,"nsolve,nsolveMM, or velcomp invalid"
           stop
          endif
         else if (project_option.eq.3) then ! viscosity
          if ((nsolve.eq.SDIM).and. &
              (nsolveMM.eq.SDIM).and. &
              (velcomp.ge.0).and. &
              (velcomp.lt.SDIM)) then
           veldir=velcomp+1
           im_vel=1
           bc_comp=veldir
          else
           print *,"nsolve,nsolveMM, or velcomp invalid"
           stop
          endif
         else
          print *,"project_option invalid"
          stop
         endif

         local_dd=one
         local_visc_coef=one
         local_constant_viscosity=1

         side=0
         if (icrit.eq.fablo(dir+1)) then
          side=1
         else if (icrit.eq.fabhi(dir+1)+1) then
          side=2
         else if ((icrit.gt.fablo(dir+1)).and. &
                  (icrit.lt.fabhi(dir+1)+1)) then
          ! do nothing
         else
          print *,"icrit invalid"
          stop
         endif

         bccrit=0
         if ((side.eq.1).or.(side.eq.2)) then 
          bccrit=presbc_in(dir+1,side,bc_comp)
         else if (side.eq.0) then
          ! do nothing
         else
          print *,"side invalid"
          stop
         endif

         caller_id=2
         call eval_face_coeff( &
           caller_id, &
           level,finest_level, &
           AL,local_AL_ice,cc_group, &
           local_dd,local_dd_group, &
           local_visc_coef, &
           nsolve,nsolveMM,im_vel,dir,veldir,project_option, &
           local_constant_viscosity,side,bccrit,local_wt)

         if (local_wt(veldir).ge.zero) then

          if (local_wt(veldir).eq.zero) then
           if ((local_gp.eq.zero).or.(maskface.eq.0).or.(1.eq.0)) then
            local_dest=local_src
           else if ((local_gp.ne.zero).and.(maskface.eq.1)) then
            print *,"local_gp invalid"
            print *,"local_gp ",local_gp
            print *,"level,finest_level ",level,finest_level
            print *,"AL,AL_ice,cc_group ",AL,AL_ice,cc_group
            print *,"veldir=",veldir
            print *,"local_wt(veldir) ",local_wt(veldir)
            print *,"local_dd,local_dd_group ",local_dd,local_dd_group
            print *,"dir,project_option,side,bccrit ", &
                dir,project_option,side,bccrit
            print *,"int_dir=",INT_DIR
            print *,"ext_dir=",EXT_DIR
            print *,"reflect_odd=",REFLECT_ODD
            print *,"reflect_even=",REFLECT_EVEN
            print *,"foextrap=",FOEXTRAP
            stop
           else
            print *,"local_gp or maskface bust"
            print *,"local_gp ",local_gp
            print *,"maskface ",maskface
            stop
           endif
          else if (local_wt(veldir).gt.zero) then
           local_dest=local_src+local_gp
          else
           print *,"local_wt(veldir) invalid (fluidsolidcor):",local_wt(veldir)
           print *,"cc_group=",cc_group
           stop
          endif

          if (dir.eq.0) then
           xdest(D_DECL(i,j,k))=local_dest
          else if (dir.eq.1) then
           ydest(D_DECL(i,j,k))=local_dest
          else if ((dir.eq.2).and.(SDIM.eq.3)) then
           zdest(D_DECL(i,j,k))=local_dest
          else
           print *,"dir invalid fluid solid cor 2"
           stop
          endif

         else
          print *,"local_wt(veldir) invalid"
          print *,"local_wt(veldir)=",local_wt(veldir)
          stop
         endif

       enddo
       enddo
       enddo
      enddo ! dir

      return
      end subroutine FORT_FLUIDSOLIDCOR

      subroutine FORT_EOS_PRESSURE( &
        level, &
        finest_level, &
        xlo,dx, &
        pres,DIMS(pres), &
        recon,DIMS(recon), &
        levelpc,DIMS(levelpc), &
        den,DIMS(den), &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        nmat, &
        nden)

      use global_utility_module
      use MOF_routines_module
      use probf90_module

      IMPLICIT NONE

      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T nmat,nden
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
    
      REAL_T xlo(SDIM)
      REAL_T dx(SDIM)

      INTEGER_T DIMDEC(pres)
      INTEGER_T DIMDEC(recon)
      INTEGER_T DIMDEC(levelpc)
      INTEGER_T DIMDEC(den)
      REAL_T pres(DIMV(pres),num_materials_vel)
      REAL_T recon(DIMV(recon),nmat*ngeom_recon)
      REAL_T levelpc(DIMV(levelpc),nmat*(1+SDIM))
      REAL_T den(DIMV(den),nden) ! den,temp
      INTEGER_T i,j,k
      INTEGER_T im
      INTEGER_T im_primary
      INTEGER_T ibase
      REAL_T rho,internal_energy,TEMP
      REAL_T LS(nmat)
      REAL_T vfrac(nmat)
      INTEGER_T vofcomp

      if (bfact.lt.1) then
       print *,"bfact invalid162"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid eos pressure"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (nden.ne.nmat*num_state_material) then
       print *,"nden invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif
      call checkbound(fablo,fabhi,DIMS(recon),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(levelpc),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(pres),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(den),1,-1,44)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1)
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       do im=1,nmat
        LS(im)=levelpc(D_DECL(i,j,k),im)
        vofcomp=(im-1)*ngeom_recon+1
        vfrac(im)=recon(D_DECL(i,j,k),vofcomp)
       enddo ! im

       call get_primary_material_VFRAC(vfrac,nmat,im_primary,1)

       if (vfrac(im_primary).le.VOFTOL) then
        print *,"vfrac too small FORT_EOS_PRESSURE"
        print *,"i,j,k ",i,j,k
        print *,"level,finest_level ",level,finest_level
        print *,"im_primary= ",im_primary
        do im=1,nmat
         print *,"im,vfrac,LS ",im,vfrac(im),LS(im)
        enddo
        stop
       endif

       ibase=(im_primary-1)*num_state_material

       if (is_rigid(nmat,im_primary).eq.0) then

        ! compressible material
        if ((fort_material_type(im_primary).gt.0).and. &
            (fort_material_type(im_primary).le.fort_max_num_eos)) then

         rho=den(D_DECL(i,j,k),ibase+1)
         if (rho.le.zero) then
          print *,"density has gone nonpos"
          stop
         endif
         TEMP=den(D_DECL(i,j,k),ibase+2)
         if (TEMP.le.zero) then
          print *,"TEMP has gone nonpos"
          stop
         endif
         ! returns energy/scale
         call INTERNAL_material(rho,TEMP, &
          internal_energy,fort_material_type(im_primary),im_primary)
         if (internal_energy.le.zero) then
          print *,"internal_energy has gone nonpos"
          stop
         endif
         ! p(energy*scale)/scale
         call EOS_material(rho,internal_energy, &
          pres(D_DECL(i,j,k),1), &
          fort_material_type(im_primary),im_primary)
        else if (fort_material_type(im_primary).eq.0) then
         ! do nothing
        else
         print *,"fort material type invalid"
         stop
        endif

       else if (is_rigid(nmat,im_primary).eq.1) then
        ! do nothing
       else
        print *,"is_rigid bust"
        stop
       endif 

      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine FORT_EOS_PRESSURE

! the error is set to zero in advection.  Then
! the error is updated when the slopes are recomputed.
! this routine is called after the last slope computation of the time step.
      subroutine FORT_PRESSURE_INDICATOR( &
        pressure_error_flag, &
        vorterr, &
        pressure_error_cutoff, &
        temperature_error_cutoff, &
        xlo,dx, &
        errnew,DIMS(errnew), &
        slrecon,DIMS(slrecon), &
        den,DIMS(den), &
        vort,DIMS(vort), &
        pres,DIMS(pres), &
        maskcov,DIMS(maskcov), &
        tilelo,tilehi,  &
        fablo,fabhi, &
        bfact, &
        level,  &
        finest_level,  &
        nmat)

      use global_utility_module
      use probf90_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T level
      INTEGER_T finest_level
      INTEGER_T nmat
      INTEGER_T pressure_error_flag
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      REAL_T xlo(SDIM)
      REAL_T dx(SDIM)

      REAL_T vorterr(nmat)
      REAL_T pressure_error_cutoff(nmat)
      REAL_T temperature_error_cutoff(nmat)
      INTEGER_T DIMDEC(maskcov)
      INTEGER_T DIMDEC(errnew)
      INTEGER_T DIMDEC(slrecon)
      INTEGER_T DIMDEC(den)
      INTEGER_T DIMDEC(vort)
      INTEGER_T DIMDEC(pres)
      REAL_T maskcov(DIMV(maskcov))
      REAL_T errnew(DIMV(errnew))
      REAL_T den(DIMV(den),nmat*num_state_material)
      REAL_T slrecon(DIMV(slrecon),nmat*ngeom_recon)
      REAL_T vort(DIMV(vort),num_materials_vel)
      REAL_T pres(DIMV(pres),num_materials_vel)
      INTEGER_T i,j,k,im
      REAL_T vfrac(nmat)
      REAL_T pres_array(D_DECL(3,3,3))
      REAL_T temp_array(D_DECL(3,3,3))
      INTEGER_T i2,j2,k2,vofcomp
      INTEGER_T kstencil_lo,kstencil_hi
      INTEGER_T tcomp
      INTEGER_T pcomp
      REAL_T local_vort
      INTEGER_T im_vort
      INTEGER_T local_mask
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf

      nhalf=3

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid163"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 42"
       stop
      endif
      if ((pressure_error_flag.ne.0).and. &
          (pressure_error_flag.ne.1)) then
       print *,"pressure_error_flag invalid"
       stop
      endif
      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(slrecon),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(errnew),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(den),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(vort),0,-1,44)
      call checkbound(fablo,fabhi,DIMS(pres),1,-1,44)

      if (SDIM.eq.3) then
       kstencil_lo=-1
       kstencil_hi=1
      else if (SDIM.eq.2) then
       kstencil_lo=0
       kstencil_hi=0
      else
       print *,"dimension bust"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
 
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       ! maskcov=tag if not covered by level+1 or outside the domain.
       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         vfrac(im)=slrecon(D_DECL(i,j,k),vofcomp)
        enddo
        call get_primary_material_VFRAC(vfrac,nmat,im,2)

        if (is_rigid(nmat,im).eq.0) then

           ! only check pressure/temperature/vorticity
           ! magnitude away from interfaces
         if (vfrac(im).gt.one-VOFTOL) then

          tcomp=(im-1)*num_state_material+2

          pcomp=1
          im_vort=1

          if (vorterr(im).lt.zero) then
           print *,"vorterr cannot be negative"
           stop
          endif

          do i2=-1,1
          do j2=-1,1
          do k2=kstencil_lo,kstencil_hi
           pres_array(D_DECL(i2+2,j2+2,k2+2))= &
             pres(D_DECL(i+i2,j+j2,k+k2),pcomp)
           temp_array(D_DECL(i2+2,j2+2,k2+2))= &
             den(D_DECL(i+i2,j+j2,k+k2),tcomp)
          enddo
          enddo
          enddo

           ! vort is initialized in FORT_GETSHEAR when onlyscalar.eq.2
           ! vort is the vorticity magnitude (L2 norm)
          local_vort=vort(D_DECL(i,j,k),im_vort)

            ! error(p*scale)
            ! errnew=max(errnew,VOFTOL)
          call EOS_error_ind( &
           pressure_error_flag, &
           xsten,nhalf,bfact, &
           local_vort, &
           pres_array, &
           temp_array, &
           errnew(D_DECL(i,j,k)), &
           fort_material_type(im), &
           vorterr(im), &
           pressure_error_cutoff(im), &
           temperature_error_cutoff(im))
         endif  ! vfrac ~ 1

        else if (is_rigid(nmat,im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid"
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

      return
      end subroutine FORT_PRESSURE_INDICATOR


! 1. calculates p(rho^n+1,e_advect) and puts it in 2nd component
!    of cell_sound.
!    Equation of state to be used depends on cvof (tessellating vfracs)
! 2. calculates 1/(rho c^2 dt^2)  and puts it in 1st component of cell_sound.
!    Equation of state to be used depends on cvof (tessellating vfracs)
! 3. pressure=0.0 in incompressible regions
!
! if project_option==10, 11, mdot=0.0 (on input) (mdot corresponds
! to localMF[DIFFUSIONRHS_MF])
!
! velocity scale: V
! time scale is : 1/V
! pressure scale: V^2
! scale for "cell_sound" is 1
! called from:
! NavierStokes::init_advective_pressure
!
      subroutine FORT_ADVECTIVE_PRESSURE( &
        level, &
        finest_level, &
        xlo,dx, &
        dt, &
        make_interface_incomp, &
        maskcov,DIMS(maskcov), &
        lsnew,DIMS(lsnew), &
        csnd,DIMS(csnd), &
        cvof,DIMS(cvof), & ! tessellating
        den,DIMS(den), &
        mdot,DIMS(mdot), & ! passed from localMF[DIFFUSIONRHS_MF]
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        nmat,nden, &
        pgrad_dt_factor, &
        pressure_select_criterion, &
        project_option)

      use global_utility_module
      use probf90_module
      use MOF_routines_module

      IMPLICIT NONE


      INTEGER_T, intent(in) :: level
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: project_option
      REAL_T, intent(in) :: pgrad_dt_factor
      INTEGER_T, intent(in) :: pressure_select_criterion
      INTEGER_T, intent(in) :: nmat,nden
      INTEGER_T, intent(in) :: tilelo(SDIM),tilehi(SDIM)
      INTEGER_T, intent(in) :: fablo(SDIM),fabhi(SDIM)
      INTEGER_T :: growlo(3),growhi(3)
      INTEGER_T, intent(in) :: bfact
      REAL_T, intent(in) :: xlo(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: dt
      INTEGER_T, intent(in) :: make_interface_incomp

      INTEGER_T, intent(in) :: DIMDEC(maskcov)
      INTEGER_T, intent(in) :: DIMDEC(lsnew)
      INTEGER_T, intent(in) :: DIMDEC(csnd)
      INTEGER_T, intent(in) :: DIMDEC(cvof)
      INTEGER_T, intent(in) :: DIMDEC(den)
      INTEGER_T, intent(in) :: DIMDEC(mdot)
      REAL_T, intent(in) :: maskcov(DIMV(maskcov))
      REAL_T, intent(in) :: lsnew(DIMV(lsnew),nmat*(1+SDIM))
      REAL_T, intent(inout) :: csnd(DIMV(csnd),2) 
      REAL_T, intent(in) :: cvof(DIMV(cvof),nmat) 
      REAL_T, intent(in) :: den(DIMV(den),nden) ! den,temp
      REAL_T, intent(inout) :: mdot(DIMV(mdot),num_materials_vel) 

      INTEGER_T i,j,k
      INTEGER_T im,imcrit,im_weight,im_opp
      INTEGER_T ibase
      REAL_T temperature,internal_energy,soundsqr
      REAL_T pres(nmat)
      REAL_T rho(nmat)
      REAL_T one_over_c(nmat)
      REAL_T one_over_c2(nmat)
      REAL_T localLS(nmat)
      REAL_T vfrac(nmat)
      REAL_T vfrac_weight(nmat)
      REAL_T vfrac_solid_sum
      REAL_T vfrac_fluid_sum
      INTEGER_T infinite_weight
      INTEGER_T local_infinite_weight
      REAL_T div_hold(num_materials_vel)
      REAL_T csound_hold
      REAL_T DXMAXLS
      REAL_T cutoff
      REAL_T rmaskcov
      INTEGER_T local_mask

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid164"
       stop
      endif
      if (nmat.eq.num_materials) then
       ! do nothing
      else
       print *,"nmat invalid advective pressure"
       stop
      endif
      if (num_materials_vel.eq.1) then
       ! do nothing
      else
       print *,"num_materials_vel invalid"
       stop
      endif
      if (num_state_base.eq.2) then
       ! do nothing
      else
       print *,"num_state_base invalid"
       stop
      endif
      if (pgrad_dt_factor.ge.one) then
       ! do nothing
      else
       print *,"pgrad_dt_factor too small"
       stop
      endif
      if (nden.eq.nmat*num_state_material) then
       ! do nothing
      else
       print *,"nden invalid"
       stop
      endif
      if ((pressure_select_criterion.ge.0).and. &
          (pressure_select_criterion.le.2)) then
       ! do nothing
      else
       print *,"pressure_select_criterion invalid"
       stop
      endif
      if ((project_option.eq.0).or. &
          (project_option.eq.10).or. &
          (project_option.eq.11).or. & !FSI_material_exists 2nd project
          (project_option.eq.13)) then !FSI_material_exists 1st project
       ! do nothing
      else
       print *,"project_option invalid advective pressure"
       stop
      endif
      if ((level.ge.0).and.(level.le.finest_level)) then
       ! do nothing
      else
       print *,"level invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(maskcov),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(lsnew),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(csnd),0,-1,44)
      call checkbound(fablo,fabhi,DIMS(cvof),0,-1,44)
      call checkbound(fablo,fabhi,DIMS(den),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(mdot),0,-1,44)

      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=two*DXMAXLS

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

       ! csound (cell_sound) initialized.
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       rmaskcov=maskcov(D_DECL(i,j,k))
       local_mask=NINT(rmaskcov)
       if (((rmaskcov.eq.one).and.(local_mask.eq.1)).or. &
           ((rmaskcov.eq.zero).and.(local_mask.eq.0))) then

        imcrit=0
        vfrac_solid_sum=zero
        vfrac_fluid_sum=zero

        do im=1,nmat
         vfrac(im)=cvof(D_DECL(i,j,k),im)
         if ((vfrac(im).ge.-VOFTOL).and.(vfrac(im).le.one+VOFTOL)) then
          ! do nothing
         else
          print *,"vfrac invalid"
          stop
         endif
         if (vfrac(im).lt.VOFTOL) then
          vfrac(im)=zero
         endif
         if (is_rigid(nmat,im).eq.1) then
          vfrac_solid_sum=vfrac_solid_sum+vfrac(im)
         else if (is_rigid(nmat,im).eq.0) then
          vfrac_fluid_sum=vfrac_fluid_sum+vfrac(im)
         else
          print *,"is_rigid(nmat,im) invalid"
          stop
         endif
         if (imcrit.eq.0) then
          imcrit=im
         else if ((imcrit.ge.1).and.(imcrit.le.nmat)) then
          if (vfrac(im).gt.vfrac(imcrit)) then
           imcrit=im
          endif
         else
          print *,"imcrit invalid"
          stop
         endif
        enddo ! im=1..nmat

        if ((imcrit.ge.1).and.(imcrit.le.nmat)) then
         ! do nothing
        else
         print *,"imcrit invalid"
         stop
        endif

        if (abs(vfrac_solid_sum+vfrac_fluid_sum-one).le.1.0E-4) then
         ! do nothing
        else
         print *,"vfrac_solid_sum+vfrac_fluid_sum invalid"
         stop
        endif

        if ((vfrac_solid_sum.ge.vfrac_fluid_sum).or. &
            (is_rigid(nmat,imcrit).eq.1)) then
         imcrit=0
        endif


        if ((project_option.eq.0).or. &
            (project_option.eq.10).or. &
            (project_option.eq.11).or. & !FSI_material_exists 2nd project
            (project_option.eq.13)) then !FSI_material_exists 1st project

         if ((project_option.eq.0).or. &
             (project_option.eq.13)) then
          div_hold(1)=zero
         else if ((project_option.eq.10).or. & !sync project prior to advection
                  (project_option.eq.11)) then !FSI_material_exists 2nd project
           ! coeff_avg,p_avg
          div_hold(1)=csnd(D_DECL(i,j,k),2)   ! pavg (copied from 1st component
                                              ! of DIV_TYPE)
           ! FORT_ADVECTIVE_PRESSURE called from 
           !   NavierStokes::init_advective_pressure
           ! init_advective_pressure called from
           !   NavierStokes::multiphase_project
           ! mdot passed from localMF[DIFFUSIONRHS_MF]
           ! project_option==11 => FSI_material_exists (2nd project)
          if (mdot(D_DECL(i,j,k),1).eq.zero) then
           ! do nothing
          else
           print *,"mdot(D_DECL(i,j,k),1).ne.zero in advective_pressure"
           print *,"i,j,k,mdot ",i,j,k,mdot(D_DECL(i,j,k),1)
           print *,"level,finest_level ",level,finest_level
           print *,"div_hold(1) ",div_hold(1)
           print *,"project_option ",project_option
           print *,"dt=",dt
           print *,"make_interface_incomp ",make_interface_incomp
           print *,"nmat,nden,bfact ",nmat,nden,bfact
           print *,"pgrad_dt_factor ",pgrad_dt_factor
           print *,"pressure_select_criterion ",pressure_select_criterion
           do im=1,nmat
            localLS(im)=lsnew(D_DECL(i,j,k),im)
            print *,"im,localLS(im) ",im,localLS(im)
            print *,"im,vfrac(im) ",im,vfrac(im)
           enddo
           print *,"imcrit ",imcrit
           print *,"vfrac_solid_sum ",vfrac_solid_sum
           print *,"vfrac_fluid_sum ",vfrac_fluid_sum
           stop
          endif
         else
          print *,"project_option invalid"
          stop
         endif

         do im=1,nmat
          localLS(im)=lsnew(D_DECL(i,j,k),im)
         enddo

         do im=1,nmat

          vfrac_weight(im)=vfrac(im)

          ibase=(im-1)*num_state_material

          if (is_rigid(nmat,im).eq.0) then

           rho(im)=den(D_DECL(i,j,k),ibase+1)  ! regular density
   
           if (rho(im).le.zero) then
            print *,"cannot have non-pos density"
            stop
           endif

            ! temperature after diffusion
           temperature=den(D_DECL(i,j,k),ibase+2) 
          
            ! returns e/scale 
           call INTERNAL_material(rho(im),temperature, &
             internal_energy,fort_material_type(im),im)
          
            ! compressible material ? 
           if ((fort_material_type(im).ge.1).and. &
               (fort_material_type(im).le.fort_max_num_eos).and. &
               (vfrac(im).gt.zero)) then

              ! returns p(e*scale)/scale
            call EOS_material(rho(im),internal_energy, &
               pres(im),fort_material_type(im),im)
              ! returns c^2(e*scale)/scale
            call SOUNDSQR_material(rho(im),internal_energy, &
               soundsqr,fort_material_type(im),im)
            if (soundsqr.gt.zero) then
             ! do nothing
            else
             print *,"cannot have non-positive sound speed"
             print *,"im= ",im
             print *,"mattype= ",fort_material_type(im)
             print *,"sound speed sqr ",soundsqr
             print *,"vfrac= ",vfrac(im)
             stop
            endif

            one_over_c2(im)=one/soundsqr
            one_over_c(im)=sqrt(one_over_c2(im))

            do im_opp=1,nmat
             if (im_opp.ne.im) then
              if (make_interface_incomp.eq.0) then
               ! do nothing
              else if ((make_interface_incomp.eq.1).or. &
                       (make_interface_incomp.eq.2)) then
               if ((abs(localLS(im)).le.cutoff).and. & ! cutoff=2 * DXMAXLS
                   (abs(localLS(im_opp)).le.cutoff)) then
                one_over_c2(im)=zero
                one_over_c(im)=zero
               else if ((abs(localLS(im)).gt.cutoff).or. &
                        (abs(localLS(im_opp)).gt.cutoff)) then
                ! do nothing
               else
                print *,"localLS invalid"
                stop
               endif  
              else
               print *,"make_interface_incomp invalid"
               stop
              endif
             else if (im_opp.eq.im) then
              ! do nothing
             else
              print *,"im_opp or im invalid"
              stop
             endif
            enddo ! im_opp=1,nmat

           else if ((fort_material_type(im).eq.0).or. &
                    (vfrac(im).eq.zero)) then
            pres(im)=zero
            one_over_c2(im)=zero
            one_over_c(im)=zero
           else
            print *,"material type or vfrac invalid"
            stop
           endif
   
          else if (is_rigid(nmat,im).eq.1) then

           rho(im)=fort_denconst(im)
           if (rho(im).gt.zero) then
            ! do nothing
           else
            print *,"rho(im)=0 ",im,rho(im)
            stop
           endif
           pres(im)=zero
           one_over_c(im)=zero
           one_over_c2(im)=zero

          else
           print *,"is_rigid bust"
           stop
          endif

         enddo ! im=1,nmat

         infinite_weight=0
         im_weight=0

         if (imcrit.eq.0) then ! rigid material(s) dominate
          im_weight=imcrit
          infinite_weight=1
         else if ((imcrit.ge.1).and.(imcrit.le.nmat)) then
  
          if (is_rigid(nmat,imcrit).eq.0) then

           do im=1,nmat

            local_infinite_weight=0

            if (is_rigid(nmat,im).eq.1) then
             ! do nothing
            else if (is_rigid(nmat,im).eq.0) then
             if ((vfrac(im).gt.zero).and.(vfrac(im).le.one+VOFTOL)) then
              if (pressure_select_criterion.eq.0) then ! vol. frac.
               ! do nothing (vfrac_weight=vfrac)
              else if (pressure_select_criterion.eq.1) then ! mass frac. 
               vfrac_weight(im)=vfrac_weight(im)*rho(im)
              else if (pressure_select_criterion.eq.2) then ! impedance weight
               if (one_over_c(im).eq.zero) then
                infinite_weight=1
                local_infinite_weight=1
               else
                vfrac_weight(im)=vfrac_weight(im)*rho(im)/one_over_c(im)
               endif
              else
               print *,"pressure_select_criterion invalid"
               stop
              endif
             else if (vfrac(im).eq.zero) then
              ! do nothing (vfrac_weight is 0)
             else 
              print *,"vfrac invalid"
              stop
             endif
             if (im_weight.eq.0) then
              im_weight=im
             else if ((im_weight.ge.1).and.(im_weight.le.nmat)) then
              if ((vfrac_weight(im).gt.vfrac_weight(im_weight)).or. &
                  (local_infinite_weight.eq.1)) then
               im_weight=im
              endif
             else
              print *,"im_weight invalid"
              stop
             endif
            else
             print *,"is_rigid(nmat,im) invalid"
             stop
            endif

           enddo ! im=1,nmat

          else
           print *,"is_rigid(nmat,imcrit) invalid"
           stop
          endif

         else
          print *,"imcrit invalid"
          stop
         endif
       
         if (imcrit.eq.0) then 

          csnd(D_DECL(i,j,k),1)=zero  ! coeff
          csnd(D_DECL(i,j,k),2)=zero  ! padvect

          if (imcrit.ne.im_weight) then
           print *,"imcrit.ne.im_weight"
           stop
          endif

         else if ((imcrit.ge.1).and.(imcrit.le.nmat)) then

          if (is_rigid(nmat,imcrit).ne.0) then
           print *,"is_rigid(nmat,imcrit) invalid"
           stop
          endif

          if (is_rigid(nmat,im_weight).ne.0) then
           print *,"is_rigid(nmat,im_weight).ne.0"
           stop
          endif

          if (infinite_weight.eq.1) then ! sound speed=infinity

           csnd(D_DECL(i,j,k),1)=zero  ! coeff
           csnd(D_DECL(i,j,k),2)=zero  ! padvect
           if ((project_option.eq.10).or. &
               (project_option.eq.11)) then
            ! mdot corresponds to localMF[DIFFUSIONRHS_MF]
            mdot(D_DECL(i,j,k),1)=div_hold(1)/dt
           else if ((project_option.eq.0).or. &
                    (project_option.eq.13)) then
            ! do nothing
           else
            print *,"project_option invalid"
            stop
           endif

          else if (infinite_weight.eq.0) then

           if ((im_weight.lt.1).or.(im_weight.gt.nmat)) then
            print *,"im_weight invalid"
            stop
           endif
           if (vfrac_weight(im_weight).le.zero) then
            print *,"vfrac_weight(im_weight) invalid"
            stop
           endif

            ! im_weight is a compressible material
           if ((fort_material_type(im_weight).ge.1).and. &
               (fort_material_type(im_weight).le.fort_max_num_eos)) then

            if (rho(im_weight).le.zero) then
             print *,"rho(im_weight) invalid"
             stop
            endif

             ! padvect (zero for the incompressible materials)
            csnd(D_DECL(i,j,k),2)=pres(im_weight)

            csound_hold=one_over_c2(im_weight)/ &
              (pgrad_dt_factor*rho(im_weight)*dt*dt)

             ! coeff
            csnd(D_DECL(i,j,k),1)=csound_hold

            if ((project_option.eq.10).or. &
                (project_option.eq.11)) then

             if (csound_hold.eq.zero) then ! incomp
              csnd(D_DECL(i,j,k),2)=zero ! padvect
              mdot(D_DECL(i,j,k),1)=div_hold(1)/dt ! localMF[DIFFUSIONRHS_MF]
             else if (csound_hold.ne.zero) then
              ! (1/(rho c^2 dt^2))p=div/dt
              ! csound_hold p = div/dt
              ! p=div/(dt csound_hold)
              csnd(D_DECL(i,j,k),2)=div_hold(1)/(csound_hold*dt)
             else
              print *,"csound_hold invalid"
              stop
             endif
            else if ((project_option.eq.0).or. &
                     (project_option.eq.13)) then
             ! do nothing
            else
             print *,"project_option invalid"
             stop
            endif

            ! im_weight is an incompressible material
           else if (fort_material_type(im_weight).eq.0) then

            csnd(D_DECL(i,j,k),1)=zero ! coeff
            csnd(D_DECL(i,j,k),2)=zero ! padvect
            if ((project_option.eq.10).or. &
                (project_option.eq.11)) then
             mdot(D_DECL(i,j,k),1)=div_hold(1)/dt ! localMF[DIFFUSIONRHS_MF]
            else if ((project_option.eq.0).or. &
                     (project_option.eq.13)) then
             ! do nothing
            else
             print *,"project_option invalid"
             stop
            endif

           else
            print *,"material type invalid"
            stop
           endif
          else
           print *,"infinite_weight invalid"
           stop
          endif

         else
          print *,"imcrit invalid"
          stop
         endif

        else
         print *,"project_option invalid advective pressure 2"
         stop
        endif

       else
        print *,"rmaskcov or local_mask invalid"
        print *,"rmaskcov=",rmaskcov
        print *,"local_mask=",local_mask
        stop
       endif
         
      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine FORT_ADVECTIVE_PRESSURE


      subroutine FORT_UPDATE_DIV( &
        xlo,dx, &
        dt, &
        csound,DIMS(csound), &
        mdot,DIMS(mdot), &
        pnew,DIMS(pnew), &
        divnew,DIMS(divnew), &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        nmat)

      use global_utility_module
      use probf90_module

      IMPLICIT NONE


      INTEGER_T nmat
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      REAL_T xlo(SDIM)
      REAL_T dx(SDIM)
      REAL_T dt

      INTEGER_T DIMDEC(csound)
      INTEGER_T DIMDEC(mdot)
      INTEGER_T DIMDEC(pnew)
      INTEGER_T DIMDEC(divnew)
      REAL_T csound(DIMV(csound),2) 
      REAL_T mdot(DIMV(mdot)) 
      REAL_T pnew(DIMV(pnew),num_materials_vel) 
      REAL_T divnew(DIMV(divnew),num_materials_vel) 
      INTEGER_T i,j,k
      REAL_T coeff_hold,compress_term,mdot_term
      INTEGER_T sound_comp

      if (bfact.lt.1) then
       print *,"bfact invalid165"
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid advective pressure"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      call checkbound(fablo,fabhi,DIMS(csound),0,-1,44)
      call checkbound(fablo,fabhi,DIMS(mdot),0,-1,44)
      call checkbound(fablo,fabhi,DIMS(pnew),1,-1,44)
      call checkbound(fablo,fabhi,DIMS(divnew),1,-1,44)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif

      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)

       sound_comp=1
       coeff_hold=csound(D_DECL(i,j,k),sound_comp) ! 1/(rho c^2 dt^2)

       if (coeff_hold.gt.zero) then
        compress_term=-dt*(pnew(D_DECL(i,j,k),1)- &
          csound(D_DECL(i,j,k),sound_comp+1))*coeff_hold
       else if (coeff_hold.eq.zero) then
        compress_term=zero
       else
        print *,"coeff_hold invalid"
        stop
       endif

       mdot_term=mdot(D_DECL(i,j,k))

       divnew(D_DECL(i,j,k),1)=compress_term+mdot_term*dt
         
      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine FORT_UPDATE_DIV


       ! spectral_override==0 => always low order
      subroutine FORT_EDGEAVGDOWN( &
       enable_spectral, &
       finest_level, &
       spectral_override, &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       edge_dir, &
       crse,DIMS(c), &
       fine,DIMS(f), &
       mask,DIMS(mask), &
       loc,hic, &
       lof,hif, &
       ncomp)

      use navierstokesf90_module
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T, intent(in) :: enable_spectral
      INTEGER_T, intent(in) :: finest_level
      INTEGER_T, intent(in) :: spectral_override
      INTEGER_T, intent(in) :: level_c
      INTEGER_T, intent(in) :: level_f
      INTEGER_T, intent(in) :: bfact_c
      INTEGER_T, intent(in) :: bfact_f
      REAL_T, intent(in) :: problo(SDIM)
      REAL_T, intent(in) :: dxf(SDIM)
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: xlo_fine(SDIM)
      INTEGER_T, intent(in) :: ncomp
      INTEGER_T, intent(in) :: DIMDEC(c)
      INTEGER_T, intent(in) :: DIMDEC(f)
      INTEGER_T, intent(in) :: DIMDEC(mask)
      INTEGER_T, intent(in) :: loc(SDIM), hic(SDIM)  ! cell centered
      INTEGER_T, intent(in) :: lof(SDIM), hif(SDIM)  ! cell centered
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T stenlo(3), stenhi(3)
      INTEGER_T mstenlo(3), mstenhi(3)
      INTEGER_T edge_dir
      INTEGER_T gridtype
      INTEGER_T dir2
      REAL_T, intent(out) :: crse(DIMV(c),ncomp)
      REAL_T, intent(in) :: fine(DIMV(f),ncomp)
      REAL_T, intent(in) :: mask(DIMV(mask))
      INTEGER_T flochi(SDIM)
      INTEGER_T ic,jc,kc
      INTEGER_T ifine,jfine,kfine
      INTEGER_T ii,jj,kk
      INTEGER_T ilocal,jlocal,klocal
      INTEGER_T inorm
      INTEGER_T n
      INTEGER_T infab
      REAL_T voltotal
      REAL_T volall

      REAL_T wt(SDIM)

      REAL_T xsten(-1:1,SDIM)
      REAL_T xstenND(-1:1,SDIM)
      INTEGER_T nhalf
      REAL_T crse_value(ncomp)
      REAL_T dxelem_f
      REAL_T xcoarse(SDIM)
      INTEGER_T finelo_index
      REAL_T, dimension(D_DECL(:,:,:),:),allocatable :: ffine
      REAL_T INTERP_TOL
      INTEGER_T khi
      INTEGER_T testmask,testmask2
      INTEGER_T nmat
      INTEGER_T local_enable_spectral
      INTEGER_T caller_id

      caller_id=8
      nhalf=1

      INTERP_TOL=1.0E-4

      nmat=num_materials

      local_enable_spectral=enable_spectral

      if (level_f.gt.finest_level) then
       print *,"level_f invalid"
       stop
      endif

      if ((spectral_override.ne.0).and. &
          (spectral_override.ne.1)) then
       print *,"spectral_override invalid"
       stop
      endif

      if (ncomp.lt.1) then
       print *,"ncomp invalid35"
       stop
      endif

      if (bfact_c.lt.1) then
       print *,"bfact invalid166"
       stop
      endif
      if (bfact_f.lt.1) then
       print *,"bfact_f invalid6 ",bfact_f
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level_f-1)) then
       print *,"level_c or level_f invalid"
       stop
      endif

      if ((edge_dir.lt.0).or. &
          (edge_dir.ge.SDIM)) then
       print *,"edge_dir out of range in edge average down"
       stop
      endif

      gridtype=edge_dir+1

      ii=0
      jj=0
      kk=0
      if (gridtype.eq.1) then
       ii=1
      else if (gridtype.eq.2) then
       jj=1
      else if ((gridtype.eq.3).and.(SDIM.eq.3)) then
       kk=1
      else
       print *,"gridtype invalid"
       stop
      endif

      do dir2=1,SDIM
       flochi(dir2)=bfact_f-1
      enddo
      flochi(gridtype)=bfact_f

      if (SDIM.eq.2) then
       khi=0
      else if (SDIM.eq.3) then
       khi=flochi(SDIM)
      else
       print *,"dimension bust"
       stop
      endif

      allocate(ffine(D_DECL(0:flochi(1),0:flochi(2),0:flochi(3)),ncomp))

      call checkbound(loc,hic,DIMS(c),0,edge_dir,1301)
      call checkbound(lof,hif,DIMS(f),0,edge_dir,1301)
      call checkbound(lof,hif,DIMS(mask),1,-1,1301)

      call growntileboxMAC(loc,hic,loc,hic,growlo,growhi,0,edge_dir)

      do ic=growlo(1),growhi(1)
      do jc=growlo(2),growhi(2)
      do kc=growlo(3),growhi(3)

       call gridstenMAC_level(xsten,ic,jc,kc,level_c,nhalf,gridtype)

        ! find stencil for surrounding fine level spectral element.
       stenlo(3)=0
       stenhi(3)=0
       do dir2=1,SDIM
        dxelem_f=dxf(dir2)*bfact_f
        xcoarse(dir2)=xsten(0,dir2)
        finelo_index=NINT( (xcoarse(dir2)-problo(dir2))/dxelem_f-half )
        stenlo(dir2)=bfact_f*finelo_index
        stenhi(dir2)=stenlo(dir2)+bfact_f-1
       enddo ! dir
       stenhi(gridtype)=stenhi(gridtype)+1

        ! check if coarse face is on a coarse element boundary.
       if (gridtype.eq.1) then
        inorm=ic
       else if (gridtype.eq.2) then
        inorm=jc
       else if ((gridtype.eq.3).and.(SDIM.eq.3)) then
        inorm=kc
       else
        print *,"gridtype invalid"
        stop
       endif
       if ((inorm/bfact_c)*bfact_c.eq.inorm) then
        stenlo(gridtype)=2*inorm
        stenhi(gridtype)=stenlo(gridtype)
       endif

       mstenlo(3)=0
       mstenhi(3)=0
       do dir2=1,SDIM
        mstenlo(dir2)=stenlo(dir2)
        mstenhi(dir2)=stenlo(dir2)+bfact_f-1
        if (stenhi(dir2).eq.stenlo(dir2)) then
         mstenlo(dir2)=stenhi(dir2)-1
         mstenhi(dir2)=stenhi(dir2)
        endif
       enddo ! dir2

       ifine=stenlo(1) 
       jfine=stenlo(2) 
       kfine=stenlo(SDIM) 

        ! lower left hand corner of element stencil.
       call gridstenND_level(xstenND,ifine,jfine,kfine,level_f,nhalf)

       do dir2=1,SDIM

        xcoarse(dir2)=xcoarse(dir2)-xstenND(0,dir2)

        if (stenhi(dir2).ne.stenlo(dir2)) then
         if (abs(xcoarse(dir2)).lt.INTERP_TOL*dxf(dir2)) then
          mstenlo(dir2)=stenlo(dir2)-1
         endif
         if (abs(xcoarse(dir2)-bfact_f*dxf(dir2)).lt. &
             INTERP_TOL*dxf(dir2)) then
          mstenhi(dir2)=stenhi(dir2)+1
         endif
        endif ! stenhi<>stenlo

        if ((xcoarse(dir2).lt.-INTERP_TOL*dxf(dir2)).or. &
            (xcoarse(dir2).gt.(INTERP_TOL+bfact_f)*dxf(dir2))) then
         print *,"xcoarse out of bounds edge avgdown"
         stop
        endif

       enddo ! dir2

       ifine=mstenlo(1) 
       jfine=mstenlo(2) 
       kfine=mstenlo(SDIM) 
       testmask=NINT(mask(D_DECL(ifine,jfine,kfine)))
       do ifine=mstenlo(1),mstenhi(1) 
       do jfine=mstenlo(2),mstenhi(2) 
       do kfine=mstenlo(3),mstenhi(3) 
        testmask2=NINT(mask(D_DECL(ifine,jfine,kfine)))
        if (testmask2.ne.testmask) then
         testmask=0
        endif
       enddo 
       enddo 
       enddo 

       do ilocal=0,flochi(1)
       do jlocal=0,flochi(2)
       do klocal=0,khi
        do n=1,ncomp
         ffine(D_DECL(ilocal,jlocal,klocal),n)=zero
        enddo
       enddo
       enddo
       enddo

       do n=1,ncomp
        crse_value(n)=zero
       enddo
       voltotal=zero

       if ((bfact_f.ge.2).and. &
           ((local_enable_spectral.eq.1).or. &
            (local_enable_spectral.eq.2)).and. &
           (testmask.ge.1).and. &
           (testmask.le.nmat).and. &
           (spectral_override.eq.1)) then

        do ifine=stenlo(1),stenhi(1)
        do jfine=stenlo(2),stenhi(2)
        do kfine=stenlo(3),stenhi(3)
         ilocal=ifine-stenlo(1)
         jlocal=jfine-stenlo(2)
         klocal=kfine-stenlo(3)
         do n=1,ncomp
          ffine(D_DECL(ilocal,jlocal,klocal),n)= &
            fine(D_DECL(ifine,jfine,kfine),n)
         enddo
        enddo
        enddo
        enddo

        call SEM_INTERP_ELEMENT( &
         ncomp,bfact_f,gridtype, &
         flochi,dxf,xcoarse,ffine,crse_value,caller_id)

        voltotal=one

       else if ((bfact_f.eq.1).or. &
                (local_enable_spectral.eq.0).or. &
                (local_enable_spectral.eq.3).or. &
                (testmask.eq.0).or. &
                (spectral_override.eq.0)) then

        call fine_subelement_stencilMAC(ic,jc,kc,stenlo,stenhi, &
         bfact_c,bfact_f,edge_dir+1)

        do ifine=stenlo(1),stenhi(1)
         if (edge_dir.eq.0) then
          call intersect_weightMAC_avg(ic,ifine,bfact_c,bfact_f,wt(1))
         else
          call intersect_weight_avg(ic,ifine,bfact_c,bfact_f,wt(1))
         endif
         if (wt(1).gt.zero) then
          do jfine=stenlo(2),stenhi(2)
           if (edge_dir.eq.1) then
            call intersect_weightMAC_avg(jc,jfine,bfact_c,bfact_f,wt(2))
           else
            call intersect_weight_avg(jc,jfine,bfact_c,bfact_f,wt(2))
           endif
           if (wt(2).gt.zero) then
            do kfine=stenlo(3),stenhi(3)
             if (SDIM.eq.3) then
              if (edge_dir.eq.2) then
               call intersect_weightMAC_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
              else
               call intersect_weight_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
              endif
             endif
             if (wt(SDIM).gt.zero) then
              volall=wt(1)
              do dir2=2,SDIM
               volall=volall*wt(dir2)
              enddo
              infab=1
              if ((ifine.lt.lof(1)).or.(ifine.gt.hif(1)+1).or. &
                  (jfine.lt.lof(2)).or.(jfine.gt.hif(2)+1)) then
               infab=0
              endif
              if (SDIM.eq.3) then
               if ((kfine.lt.lof(SDIM)).or.(kfine.gt.hif(SDIM)+1)) then
                infab=0
               endif
              endif
              if (infab.eq.1) then
               do n=1,ncomp
                crse_value(n) = crse_value(n)+ &
                 volall*fine(D_DECL(ifine,jfine,kfine),n)
               enddo
               voltotal=voltotal+volall
              else if (infab.eq.0) then
               ! do nothing
              else
               print *,"infab invalid"
               stop
              endif
             endif ! wt(sdim).gt.0
            enddo ! kfine
           endif
          enddo ! jfine
         endif
        enddo ! ifine

       else
        print *,"parameter bust in edgeavgdown"
        stop
       endif

       if (voltotal.le.zero) then
        print *,"voltotal invalid edgeavgdown"
        stop
       endif

       do n=1,ncomp
        crse(D_DECL(ic,jc,kc),n)=crse_value(n)/voltotal
        if ((levelrz.eq.0).or.(levelrz.eq.3)) then
         ! do nothing
        else if (levelrz.eq.1) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
         if ((xsten(0,1).le.VOFTOL*dx(1)).and. &
             (edge_dir.eq.0)) then
          crse(D_DECL(ic,jc,kc),n) = zero
         endif
        else
         print *,"levelrz invalid edge avgdown"
         stop
        endif
       enddo ! n

      enddo 
      enddo
      enddo  ! ic,jc,kc

      deallocate(ffine)

      return
      end subroutine FORT_EDGEAVGDOWN


      subroutine FORT_VOFFLUX_CORRECT( &
       finest_level, &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       edge_dir, &
       crse,DIMS(c), &
       fine,DIMS(f), &
       loc,hic, &
       lof,hif, &
       nmat)

      use navierstokesf90_module
      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T finest_level
      INTEGER_T level_c
      INTEGER_T level_f
      INTEGER_T bfact_c
      INTEGER_T bfact_f
      REAL_T problo(SDIM)
      REAL_T dxf(SDIM)
      REAL_T dx(SDIM)
      REAL_T xlo_fine(SDIM)
      INTEGER_T nmat
      INTEGER_T DIMDEC(c)
      INTEGER_T DIMDEC(f)
      INTEGER_T loc(SDIM), hic(SDIM)  ! cell centered
      INTEGER_T lof(SDIM), hif(SDIM)  ! cell centered
      INTEGER_T growlo(3), growhi(3)
      INTEGER_T growlo_c(3), growhi_c(3)
      INTEGER_T edge_dir
      REAL_T crse(DIMV(c),nmat*num_materials_vel)
      REAL_T fine(DIMV(f),nmat*num_materials_vel)
      INTEGER_T ic,jc,kc
      INTEGER_T ifine,jfine,kfine
      INTEGER_T dir,side
      INTEGER_T im

      if (level_f.gt.finest_level) then
       print *,"level_f invalid"
       stop
      endif

      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (num_materials_vel.ne.1) then
       print *,"num_materials_vel invalid"
       stop
      endif

      if (bfact_c.ne.1) then
       print *,"bfact_c invalid ",bfact_c
       stop
      endif
      if (bfact_f.ne.1) then
       print *,"bfact_f invalid ",bfact_f
       stop
      endif
      if ((level_c.lt.0).or. &
          (level_c.ne.level_f-1)) then
       print *,"level_c or level_f invalid"
       stop
      endif

      do dir=1,SDIM
       if (abs(two*dxf(dir)-dx(dir)).gt.1.0E-12*dx(dir)) then
        print *,"dxf or dx invalid"
        stop
       endif
      enddo

      if ((edge_dir.lt.0).or. &
          (edge_dir.ge.SDIM)) then
       print *,"edge_dir out of range in vofflux correct"
       stop
      endif

      call checkbound(loc,hic,DIMS(c),0,edge_dir,1301)
      call checkbound(lof,hif,DIMS(f),0,edge_dir,1301)

      do side=0,1

       call growntileboxMAC(lof,hif,lof,hif,growlo,growhi,0,edge_dir)
       call growntileboxMAC(loc,hic,loc,hic,growlo_c,growhi_c,0,edge_dir)
       if (side.eq.0) then
        growhi(edge_dir+1)=growlo(edge_dir+1)
       else if (side.eq.1) then
        growlo(edge_dir+1)=growhi(edge_dir+1)
       else
        print *,"side invalid"
        stop
       endif
       growlo_c(edge_dir+1)=growlo(edge_dir+1)/2
       growhi_c(edge_dir+1)=growhi(edge_dir+1)/2

       if ((growlo(edge_dir+1).ne.growhi(edge_dir+1)).or. &
           (2*growlo_c(edge_dir+1).ne.growlo(edge_dir+1))) then
        print *,"growlo,growhi,growlo_c, or growhi_c invalid"
        stop
       endif

       do ic=growlo_c(1),growhi_c(1)
       do jc=growlo_c(2),growhi_c(2)
       do kc=growlo_c(3),growhi_c(3)
        do im=1,nmat*num_materials_vel
         crse(D_DECL(ic,jc,kc),im)=zero
        enddo
       enddo  
       enddo  
       enddo  

       do ifine=growlo(1),growhi(1)
       do jfine=growlo(2),growhi(2)
       do kfine=growlo(3),growhi(3)
        ic=ifine/2
        jc=jfine/2
        kc=kfine/2 

        if ((ic.lt.growlo_c(1)).or. &
            (ic.gt.growhi_c(1)).or. &
            (jc.lt.growlo_c(2)).or. &
            (jc.gt.growhi_c(2)).or. &
            (kc.lt.growlo_c(3)).or. &
            (kc.gt.growhi_c(3))) then
         print *,"ic,jc, or kc invalid"
         stop
        endif
        
        do im=1,nmat*num_materials_vel
         crse(D_DECL(ic,jc,kc),im)=crse(D_DECL(ic,jc,kc),im)+ &
           fine(D_DECL(ifine,jfine,kfine),im)
        enddo
       enddo 
       enddo 
       enddo  ! ifine,jfine,kfine

      enddo ! side=0,1

      return
      end subroutine FORT_VOFFLUX_CORRECT


      subroutine FORT_METRICS( &
       xlo,dx, &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz), &
       vol,DIMS(vol), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact,level, &
       ngrow,rzflag)

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      INTEGER_T level
      INTEGER_T tilelo(SDIM),tilehi(SDIM)
      INTEGER_T fablo(SDIM),fabhi(SDIM)
      INTEGER_T growlo(3),growhi(3)
      INTEGER_T bfact
      INTEGER_T ngrow,rzflag
      INTEGER_T DIMDEC(vol)
      INTEGER_T DIMDEC(areax)
      INTEGER_T DIMDEC(areay)
      INTEGER_T DIMDEC(areaz)
      
      REAL_T  vol(DIMV(vol))
      REAL_T  areax(DIMV(areax))
      REAL_T  areay(DIMV(areay))
      REAL_T  areaz(DIMV(areaz))
      REAL_T  xlo(SDIM),dx(SDIM)

      INTEGER_T i,j,k,dir
      REAL_T xsten(-3:3,SDIM)
      INTEGER_T nhalf,side
      REAL_T local_area
      REAL_T areacen(SDIM)

      nhalf=3
      if (bfact.lt.1) then
       print *,"bfact invalid167" 
       stop
      endif
      if (levelrz.ne.rzflag) then 
       print *,"levelrz,rzflag mismatch"
       stop
      endif
      if (SDIM.eq.3) then
       if (levelrz.eq.1) then
        print *,"cannot have levelrz=1"
        stop
       endif
      endif
      call checkbound(fablo,fabhi,DIMS(vol),ngrow,-1,1301)
      call checkbound(fablo,fabhi,DIMS(areax),ngrow,0,1301)
      call checkbound(fablo,fabhi,DIMS(areay),ngrow,1,1301)
      call checkbound(fablo,fabhi,DIMS(areaz),ngrow,SDIM-1,1301)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow)
      do i=growlo(1),growhi(1)
      do j=growlo(2),growhi(2)
      do k=growlo(3),growhi(3)
       call gridsten_level(xsten,i,j,k,level,nhalf)
       call gridvol(xsten,nhalf,rzflag,vol(D_DECL(i,j,k)))
      enddo
      enddo
      enddo
     
      do dir=0,SDIM-1
       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow,dir)
       do i=growlo(1),growhi(1)
       do j=growlo(2),growhi(2)
       do k=growlo(3),growhi(3)
        call gridsten_level(xsten,i,j,k,level,nhalf)
        side=0
        call gridarea(xsten,nhalf,rzflag,dir,side,local_area,areacen)
        if (dir.eq.0) then
         areax(D_DECL(i,j,k))=local_area
        else if (dir.eq.1) then
         areay(D_DECL(i,j,k))=local_area
        else if ((dir.eq.2).and.(SDIM.eq.3)) then
         areaz(D_DECL(i,j,k))=local_area
        else
         print *,"dir invalid metrics"
         stop
        endif 
       enddo
       enddo
       enddo ! i,j,k
      enddo ! dir  

      return
      end subroutine FORT_METRICS


