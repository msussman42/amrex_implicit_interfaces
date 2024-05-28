#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

! N_EXTRA_REAL.H is in the amrlib directory.
#include "N_EXTRA_REAL.H"
#include "INTEGRATED_QUANTITY.H"
#include "EXTRAP_COMP.H"
#include "NAVIERSTOKES_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

      module navierstokesf90_module
      use amrex_fort_module, only : amrex_real
      use probf90_module

      contains

      subroutine vinterp(valu,gridval,gridx,gridy,gridz, &
        IV0,IV1,XX,YY,ZZ,icomp)

      IMPLICIT NONE

      real(amrex_real) valu
      real(amrex_real) gridx(8),gridy(8),gridz(8),gridval(8)
      integer IV0,IV1
      real(amrex_real) XX(SDIM),YY(SDIM),ZZ(SDIM)
      integer icomp
      real(amrex_real) tt

      if (abs(gridval(IV0)-valu).le.EPS10) then
       XX(icomp)=gridx(IV0)
       YY(icomp)=gridy(IV0)
#if (AMREX_SPACEDIM==3)
       ZZ(icomp)=gridz(IV0)
#endif
      else if (abs(gridval(IV1)-valu).le.EPS10) then
       XX(icomp)=gridx(IV1)
       YY(icomp)=gridy(IV1)
#if (AMREX_SPACEDIM==3)
       ZZ(icomp)=gridz(IV1)
#endif
      else
       tt=(gridval(IV0)-valu)/(gridval(IV0)-gridval(IV1))
       if ((tt.ge.zero).and.(tt.le.one)) then
        !do nothing
       else
        print *,"tt invalid 62: ",tt
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

      integer itri
      integer imaxtri
      real(amrex_real) geom(SDIM,imaxtri)
      real(amrex_real) XX(SDIM),YY(SDIM),ZZ(SDIM),NORMAL(SDIM)
      real(amrex_real) CPROD(SDIM),VEC1(3),VEC2(3),DOTPROD

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

      integer itri
      integer IV0,IV1,IV2,IV3
      integer imaxtri
      real(amrex_real) geom(SDIM,imaxtri)
      real(amrex_real) valu
      real(amrex_real) gridx(8),gridy(8),gridz(8),gridval(8)
      real(amrex_real) AA(4,4),sourcex(4),bb(4),NORMAL(SDIM)
      real(amrex_real) XX(SDIM),YY(SDIM),ZZ(SDIM)
      integer istat,idxtri

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

      integer itri
      integer IV0,IV1,IV2,IV3
      integer imaxtri
      real(amrex_real) geom(SDIM,imaxtri)
      real(amrex_real) valu
      real(amrex_real) gridx(8),gridy(8),gridz(8),gridval(8)
      real(amrex_real) AA(4,4),sourcex(4),bb(4),NORMAL(SDIM)
      real(amrex_real) XX(SDIM),YY(SDIM),ZZ(SDIM)
      integer istat,idxtri

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
       visual_revolve, &
       plotint, &
       nparts, &
       nparts_def, &
       im_solid_map)
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: plot_sdim
      integer klo_plot,khi_plot

      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def)
      integer, INTENT(in) :: total_number_grids
      integer, INTENT(in) :: num_levels
      integer, INTENT(in) :: grids_per_level_array(num_levels)
      integer, INTENT(in) :: levels_array(total_number_grids)
      integer, INTENT(in) :: bfact_array(total_number_grids)
      integer, INTENT(in) :: gridno_array(total_number_grids)
      integer, INTENT(in) :: gridlo_array(total_number_grids*SDIM)
      integer, INTENT(in) :: gridhi_array(total_number_grids*SDIM)
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nsteps
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: visual_revolve
      integer, INTENT(in) :: plotint

      integer strandid

      integer nwrite3d,nwrite2d
      integer index3d,index2d

      character*3 levstr
      character*5 gridstr
      character*32 filename32

      character*6 stepstr
      character*16 newfilename16

      integer i,j,k,dir
      integer ilev,igrid,ivel2d,ivel3d
      integer lo(plot_sdim),hi(plot_sdim)
      integer hi_index_shift(3)

! Guibo

      character*80 Title,Zonename
      real(tecplot_real_short) :: ZONEMARKER,EOHMARKER
      integer*4 :: nzones_gb,iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb
      integer bfact,testlev,testgridno
      integer testlo(SDIM),testhi(SDIM)

      ! define zone structure
      type zone3d_t
         real(tecplot_real), pointer :: var(:,:,:,:)
      end type zone3d_t
      type(zone3d_t), dimension(:), allocatable :: zone3d_gb

      type zone2d_t
         real(tecplot_real), pointer :: var(:,:,:)
      end type zone2d_t
      type(zone2d_t), dimension(:), allocatable :: zone2d_gb

      real(amrex_real) theta,rr,zz,xx,yy,ur,uz,ux,uy
      integer plot_sdim_macro
      integer add_sub_cells

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

      if (levelrz.ne.COORDSYS_RZ) then
       print *,"levelrz invalid zones revolve"
       stop
      endif
      if (visual_revolve.lt.1) then
       print *,"visual_revolve: ",visual_revolve
       print *,"visual_revolve invalid zones_revolve"
       stop
      endif

      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid zones_revolve"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid zones_revolve"
       stop
      endif

      plot_sdim_macro=SDIM
      nwrite2d=PLOTCOMP_NCOMP

      plot_sdim_macro=plot_sdim
      nwrite3d=PLOTCOMP_NCOMP

      if (nwrite3d-nwrite2d.eq.PLOTCOMP_DIFF) then
       ! do nothing
      else
       print *,"nwrite3d-nwrite2d invalid"
       print *,"nwrite3d: ",nwrite3d
       print *,"nwrite2d: ",nwrite2d
       print *,"PLOTCOMP_DIFF: ",PLOTCOMP_DIFF
       stop
      endif

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
       write(filename32,'(A14,A10,A3,A5)') &
          './temptecplot/','tempnddata',levstr,gridstr
       open(unit=4,file=filename32)

       do dir=1,plot_sdim
        if ((dir.eq.1).or.(dir.eq.2)) then
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
       enddo ! dir=1..sdim
 
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
      ! dumpstring_headers is declared in GLOBALUTIL.F90
      add_sub_cells=0
      call dumpstring_headers(plot_sdim,add_sub_cells)

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
       write(11) strandid-1    ! StrandID (this does not work)
       write(11) round_time(time) ! Solution time
       write(11) -1   ! Not used. Set to -1
       write(11) 0    ! Zone Type
       write(11) 0    ! Specify Var Location. 0 = Don't specify, 
                      ! all data is located at the nodes.
       write(11) 0    ! Are raw local 1-to-1 face neighbors supplied?
       write(11) 0    ! Number of miscellaneous user-defined  
                      !face neighbor connections

       hi_index_shift(1)=(hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1)+1
       hi_index_shift(2)=(hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1)+1
       if (plot_sdim.eq.3) then
        hi_index_shift(3)=(hi_gb(iz_gb,plot_sdim)-lo_gb(iz_gb,plot_sdim)+1)+1
       else
        print *,"plot_sdim invalid in zones_revolve_sanity"
        stop
       endif

        ! ----- IMax,JMax,KMax
       write(11) hi_index_shift(1)
       write(11) hi_index_shift(2)
       if (plot_sdim.eq.3) then
        write(11) hi_index_shift(3)
       else
        print *,"plot_sdim invalid in zones_revolve_sanity"
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
        print *,"expecting plot_sdim==2 in zones_revolve"
        stop
       else
        print *,"plot_sdim invalid"
        stop
       endif

       hi_index_shift(1)=(hi(1)-lo(1)+1)
       hi_index_shift(2)=(hi(2)-lo(2)+1)
       hi_index_shift(3)=khi_plot-klo_plot

       allocate(zone3d_gb(iz_gb)% &
        var(nwrite3d, &
            0:hi_index_shift(1), &
            0:hi_index_shift(2), &
            0:hi_index_shift(3)))

       allocate(zone2d_gb(iz_gb)% &
        var(nwrite2d, &
            0:hi_index_shift(1), &
            0:hi_index_shift(SDIM)))
        
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

       write(filename32,'(A14,A10,A3,A5)') &
          './temptecplot/','tempnddata',levstr,gridstr
       open(unit=4,file=filename32)
       print *,"filename32 ",filename32

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

       enddo  ! dir = 1..sdim

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
       do j=0,hi_index_shift(2)
       do i=0,hi_index_shift(1)
        read(4,*) (zone2d_gb(iz_gb)%var(ivar_gb,i,j),ivar_gb=1,nwrite2d)

        do k=0,hi_index_shift(3)

         index3d=0
         index2d=0

         if (plot_sdim_macro.eq.3) then
          ! do nothing
         else
          print *,"plot_sdim_macro invalid"
          stop
         endif

         theta=two*Pi*k/visual_revolve
         rr=zone2d_gb(iz_gb)%var(PLOTCOMP_X+1,i,j)
         zz=zone2d_gb(iz_gb)%var(PLOTCOMP_Y+1,i,j)
         xx=rr*cos(theta)
         yy=rr*sin(theta)
         zone3d_gb(iz_gb)%var(PLOTCOMP_X+1,i,j,k)=xx
         zone3d_gb(iz_gb)%var(PLOTCOMP_Y+1,i,j,k)=yy
         zone3d_gb(iz_gb)%var(PLOTCOMP_Z+1,i,j,k)=zz

         index3d=index3d+plot_sdim
         index2d=index2d+SDIM

         ivel2d=0
         ivel3d=0
         ur=zone2d_gb(iz_gb)%var(index2d+1+ivel2d,i,j)
         uz=zone2d_gb(iz_gb)%var(index2d+2+ivel2d,i,j)
         ux=ur*cos(theta)
         uy=ur*sin(theta)
         zone3d_gb(iz_gb)%var(index3d+1+ivel3d,i,j,k)=ux
         zone3d_gb(iz_gb)%var(index3d+2+ivel3d,i,j,k)=uy
         zone3d_gb(iz_gb)%var(index3d+3+ivel3d,i,j,k)=uz

         index3d=index3d+plot_sdim
         index2d=index2d+SDIM

          ! pressure,presder,div,divdat,mach,vof
         do ivar_gb=1,5+num_materials
          index3d=index3d+1
          index2d=index2d+1
          zone3d_gb(iz_gb)%var(index3d,i,j,k)= &
            zone2d_gb(iz_gb)%var(index2d,i,j)
         enddo

         if (index3d.eq.PLOTCOMP_LS) then
          ! do nothing
         else
          print *,"(index3d.ne.PLOTCOMP_LS)"
          stop
         endif

          ! ls
         do ivar_gb=1,num_materials
          index3d=index3d+1
          index2d=index2d+1
          zone3d_gb(iz_gb)%var(index3d,i,j,k)= &
            zone2d_gb(iz_gb)%var(index2d,i,j)
         enddo
          ! ls normal
         do ivar_gb=1,num_materials
          ur=zone2d_gb(iz_gb)%var(index2d+1,i,j)
          uz=zone2d_gb(iz_gb)%var(index2d+2,i,j)
          ux=ur*cos(theta)
          uy=ur*sin(theta)
          zone3d_gb(iz_gb)%var(index3d+1,i,j,k)=ux
          zone3d_gb(iz_gb)%var(index3d+2,i,j,k)=uy
          zone3d_gb(iz_gb)%var(index3d+3,i,j,k)=uz
          index3d=index3d+plot_sdim
          index2d=index2d+SDIM
         enddo !ivar_gp=1..num_materials

         if (index3d.eq.PLOTCOMP_SCALARS) then
          ! do nothing
         else
          print *,"(index3d.ne.PLOTCOMP_SCALARS)"
          stop
         endif

          ! den,den_merge,mom_den,configuration tensor
         do ivar_gb=1,num_materials*num_state_material+ &
              num_state_material+ & !plotcomp_scalars_merge
              num_materials+ & !mom_den
              num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE
          index3d=index3d+1
          index2d=index2d+1
          zone3d_gb(iz_gb)%var(index3d,i,j,k)= &
            zone2d_gb(iz_gb)%var(index2d,i,j)
         enddo ! do ivar_gb=1,num_materials*num_state_material+
               !            num_state_material+ &
               !            num_materials+viscoelastic stuff

         if (index3d.eq.PLOTCOMP_VISC) then
          ! do nothing
         else
          print *,"(index3d.ne.PLOTCOMP_VISC)"
          stop
         endif

          ! visc,conduct,trace
         do ivar_gb=1,num_materials+num_materials+5*num_materials
          index3d=index3d+1
          index2d=index2d+1
          zone3d_gb(iz_gb)%var(index3d,i,j,k)= &
            zone2d_gb(iz_gb)%var(index2d,i,j)
         enddo

         if (index3d.eq.PLOTCOMP_F_ELASTIC_X) then
          ! do nothing
         else
          print *,"(index3d.ne.PLOTCOMP_F_ELASTIC_X)"
          stop
         endif

          ! elastic force
         ur=zone2d_gb(iz_gb)%var(index2d+1,i,j)
         uz=zone2d_gb(iz_gb)%var(index2d+2,i,j)
         ux=ur*cos(theta)
         uy=ur*sin(theta)
         zone3d_gb(iz_gb)%var(index3d+1,i,j,k)=ux
         zone3d_gb(iz_gb)%var(index3d+2,i,j,k)=uy
         zone3d_gb(iz_gb)%var(index3d+3,i,j,k)=uz
         index3d=index3d+plot_sdim
         index2d=index2d+SDIM

         if (index3d.eq.PLOTCOMP_GRAD_VELOCITY) then
          ! do nothing
         else
          print *,"(index3d.ne.PLOTCOMP_GRAD_VELOCITY)"
          stop
         endif

          ! GRADVEL
          ! grad velocity
         do ivar_gb=1,AMREX_SPACEDIM_SQR
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

        enddo ! do k=0,hi_index_shift(3)

       enddo ! do i=0,hi_index_shift(1)
       enddo ! do j=0,hi_index_shift(2)

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
        hi_index_shift(1)=(hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1)
        hi_index_shift(2)=(hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1)
        hi_index_shift(3)=khi_plot-klo_plot
        do k=0,hi_index_shift(3)
        do j=0,hi_index_shift(2)
        do i=0,hi_index_shift(1)
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

      enddo  ! iz_gb=1,nzones_gb

      deallocate(zone3d_gb)
      deallocate(zone2d_gb)
      deallocate(lo_gb)
      deallocate(hi_gb)

      close(11)
     
      return
      end subroutine zones_revolve

      subroutine project_tangent(vel,nn)
      IMPLICIT NONE
 
      real(amrex_real) vel(SDIM)
      real(amrex_real) nn(SDIM)
      real(amrex_real) velnormal
      integer dir
    
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
  integer, INTENT(IN)  ::  N
  integer              ::  I
  real(amrex_real), INTENT(OUT)    ::  minrad
  real(amrex_real), INTENT(IN), DIMENSION(0:N)    ::  Z
  real(amrex_real), INTENT(IN), DIMENSION(0:N)    ::  F

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
  integer, INTENT(IN)  ::  N
  integer              ::  I
  real(amrex_real), INTENT(OUT)    ::  zcenter
  real(amrex_real), INTENT(IN), DIMENSION(0:N)    ::  Z
  real(amrex_real), INTENT(IN), DIMENSION(0:N)    ::  F
  real(amrex_real) sumf,sumzf

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
  integer, INTENT(IN)  ::  N
  integer              ::  I
  real(amrex_real), INTENT(IN)     ::  WAVELENGTH
  real(amrex_real), INTENT(OUT)    ::  AMPLITUDE
  real(amrex_real)                 ::  PI, H, SR, SI 
  real(amrex_real), INTENT(IN), DIMENSION(0:N)    ::  Z
  real(amrex_real), INTENT(IN), DIMENSION(0:N)    ::  F
  real(amrex_real), DIMENSION(0:N)   ::  FR, FI

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
  integer, INTENT(IN)  ::  N,M,NN
  integer              ::  I,J
  real(amrex_real), INTENT(IN)     ::  WAVELENGTH
  real(amrex_real), INTENT(OUT)    ::  AMPLITUDE
  real(amrex_real)                 ::  PI, HT, HZ, SR, SI
  real(amrex_real), INTENT(IN), DIMENSION(0:N)    ::  T
  real(amrex_real), INTENT(IN), DIMENSION(0:M)    ::  Z
  real(amrex_real), INTENT(IN), DIMENSION(0:N,0:M)    ::  F
  real(amrex_real), DIMENSION(0:M)  ::  TR, TI
  real(amrex_real), DIMENSION(0:N)  ::  ER, EI, TER, TEI
  
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
  integer, INTENT (IN) :: N
  integer :: I
  real(amrex_real), INTENT (IN) :: H
  real(amrex_real) :: S0,S1,S2
  real(amrex_real), INTENT (OUT) :: S
  real(amrex_real), INTENT (IN), DIMENSION (N) :: FI
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

! 1. calculates p(rho^n+1,e_advect) and puts it in 2nd component
!    of cell_sound.
!    Equation of state to be used depends on cvof (tessellating vfracs)
! 2. calculates 1/(rho c^2 dt^2)  and puts it in 1st component of cell_sound.
!    Equation of state to be used depends on cvof (tessellating vfracs)
! 3. pressure=0.0 in incompressible regions
!
! velocity scale: V
! time scale is : 1/V
! pressure scale: V^2
! scale for "cell_sound" is 1
! called from:
! NavierStokes::init_advective_pressure
!
      subroutine fort_advective_pressure( &
        level, &
        finest_level, &
        xlo,dx, &
        dt, &
        maskcov,DIMS(maskcov), &
        vol,DIMS(vol), &
        lsnew,DIMS(lsnew), &
        csnd,DIMS(csnd), &
        cvof,DIMS(cvof), & ! tessellating
        den,DIMS(den), &
        mdot,DIMS(mdot), & ! passed from localMF[DIFFUSIONRHS_MF]
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        nden, &
        compressible_dt_factor, &
        pressure_select_criterion, &
        project_option) &
      bind(c,name='fort_advective_pressure')

      use global_utility_module
      use probf90_module
      use MOF_routines_module

      IMPLICIT NONE


      integer, INTENT(in) :: nden
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: project_option
      real(amrex_real), INTENT(in) :: compressible_dt_factor(num_materials)
      integer, INTENT(in) :: pressure_select_criterion
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt

      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(vol)
      integer, INTENT(in) :: DIMDEC(lsnew)
      integer, INTENT(in) :: DIMDEC(csnd)
      integer, INTENT(in) :: DIMDEC(cvof)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(mdot)
      real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: vol(DIMV(vol))
      real(amrex_real), pointer :: vol_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: &
              lsnew(DIMV(lsnew),num_materials*(1+SDIM))
      real(amrex_real), pointer :: lsnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: csnd(DIMV(csnd),2) 
      real(amrex_real), pointer :: csnd_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: cvof(DIMV(cvof),num_materials) 
      real(amrex_real), pointer :: cvof_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: den(DIMV(den),nden) ! den,temp
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: mdot(DIMV(mdot)) 
      real(amrex_real), pointer :: mdot_ptr(D_DECL(:,:,:))

      integer i,j,k
      integer im,imcrit,im_weight,im_opp
      integer iten
      integer ibase
      real(amrex_real) temperature,internal_energy,soundsqr
      real(amrex_real) pres(num_materials)
      real(amrex_real) rho(num_materials)
      real(amrex_real) one_over_c(num_materials)
      real(amrex_real) one_over_c2(num_materials)
      real(amrex_real) localLS(num_materials)
      real(amrex_real) vfrac(num_materials)
      real(amrex_real) vfrac_weight(num_materials)
      real(amrex_real) vfrac_solid_sum
      real(amrex_real) vfrac_fluid_sum
      integer infinite_weight
      integer local_infinite_weight
      real(amrex_real) csound_hold
      real(amrex_real) DXMAXLS
      real(amrex_real) cutoff
      real(amrex_real) rmaskcov
      integer local_mask
      real(amrex_real) massfrac_parm(num_species_var+1)
      integer ispec
      real(amrex_real) local_volume

      if (bfact.ge.1) then
       ! do nothing
      else
       print *,"bfact invalid164"
       stop
      endif
      if (num_state_base.eq.2) then
       ! do nothing
      else
       print *,"num_state_base invalid"
       stop
      endif

      do im=1,num_materials
       if ((compressible_dt_factor(im).ge.one).and. &
           (compressible_dt_factor(im).le.1.0D+20)) then
        ! do nothing
       else
        print *,"compressible_dt_factor invalid: ",im, &
                compressible_dt_factor(im)
        stop
       endif
      enddo

      if (nden.eq.num_materials*num_state_material) then
       ! do nothing
      else
       print *,"nden invalid: ",nden
       stop
      endif
      if ((pressure_select_criterion.ge.0).and. &
          (pressure_select_criterion.le.2)) then
       ! do nothing
      else
       print *,"pressure_select_criterion invalid: ", &
               pressure_select_criterion
       stop
      endif
      if (project_option.eq.SOLVETYPE_PRES) then
       ! do nothing
      else
       print *,"project_option invalid advective pressure"
       stop
      endif
      if ((level.ge.0).and.(level.le.finest_level)) then
       ! do nothing
      else
       print *,"level invalid: ",level
       stop
      endif

      maskcov_ptr=>maskcov
      vol_ptr=>vol
      lsnew_ptr=>lsnew
      csnd_ptr=>csnd
      cvof_ptr=>cvof
      den_ptr=>den
      mdot_ptr=>mdot

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,vol_ptr,0,-1)
      call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,csnd_ptr,0,-1)
      call checkbound_array(fablo,fabhi,cvof_ptr,0,-1)
      call checkbound_array(fablo,fabhi,den_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,mdot_ptr,0,-1)

      call get_dxmaxLS(dx,bfact,DXMAXLS)
      cutoff=two*DXMAXLS

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

       ! csound (cell_sound) initialized.
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       local_volume=vol(D_DECL(i,j,k))
       if (local_volume.gt.zero) then
        ! do nothing
       else
        print *,"local_volume must be positive: ",local_volume
        stop
       endif

       rmaskcov=maskcov(D_DECL(i,j,k))
       local_mask=NINT(rmaskcov)
       if (((rmaskcov.eq.one).and.(local_mask.eq.1)).or. &
           ((rmaskcov.eq.zero).and.(local_mask.eq.0))) then

        imcrit=0
        vfrac_solid_sum=zero
        vfrac_fluid_sum=zero

        do im=1,num_materials
         vfrac(im)=cvof(D_DECL(i,j,k),im)
         if ((vfrac(im).ge.-EPS1).and.(vfrac(im).le.one+EPS1)) then
          ! do nothing
         else
          print *,"vfrac invalid: ",im,vfrac(im)
          stop
         endif
         if (vfrac(im).lt.VOFTOL) then
          vfrac(im)=zero
         endif
         if (is_rigid(im).eq.1) then
          vfrac_solid_sum=vfrac_solid_sum+vfrac(im)
         else if (is_rigid(im).eq.0) then
          vfrac_fluid_sum=vfrac_fluid_sum+vfrac(im)
         else
          print *,"is_rigid(im) invalid"
          stop
         endif
         if (imcrit.eq.0) then
          imcrit=im
         else if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then
          if (vfrac(im).gt.vfrac(imcrit)) then
           imcrit=im
          endif
         else
          print *,"imcrit invalid"
          stop
         endif
        enddo ! im=1..num_materials

        if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then
         ! do nothing
        else
         print *,"imcrit invalid"
         stop
        endif

        if (abs(vfrac_solid_sum+vfrac_fluid_sum-one).le.EPS3) then
         ! do nothing
        else
         print *,"vfrac_solid_sum+vfrac_fluid_sum invalid", &
            vfrac_solid_sum+vfrac_fluid_sum
         stop
        endif

        if ((vfrac_solid_sum.ge.vfrac_fluid_sum).or. &
            (is_rigid(imcrit).eq.1)) then
         imcrit=0
        endif


        if (project_option.eq.SOLVETYPE_PRES) then

         do im=1,num_materials
          localLS(im)=lsnew(D_DECL(i,j,k),im)
         enddo

         do im=1,num_materials

          vfrac_weight(im)=vfrac(im)

          ibase=(im-1)*num_state_material

          if (is_rigid(im).eq.0) then

           rho(im)=den(D_DECL(i,j,k),ibase+ENUM_DENVAR+1)  ! regular density
   
           if (rho(im).gt.zero) then
            ! do nothing
           else
            print *,"cannot have non-pos density: ",im,rho(im)
            stop
           endif

            ! temperature after diffusion
           temperature=den(D_DECL(i,j,k),ibase+ENUM_TEMPERATUREVAR+1) 

           call init_massfrac_parm(rho(im),massfrac_parm,im)
           do ispec=1,num_species_var
            massfrac_parm(ispec)=den(D_DECL(i,j,k),ibase+ENUM_SPECIESVAR+ispec)
           enddo
          
            ! returns e/scale 
           call INTERNAL_material(rho(im),massfrac_parm, &
             temperature, &
             internal_energy,fort_material_type(im),im)
          
            ! compressible material ? 
           if ((fort_material_type(im).ge.1).and. &
               (fort_material_type(im).le.MAX_NUM_EOS).and. &
               (vfrac(im).gt.zero)) then

              ! returns p(e*scale)/scale
            call EOS_material(rho(im),massfrac_parm, &
               internal_energy, &
               pres(im),fort_material_type(im),im)
              ! returns c^2(e*scale)/scale
            call SOUNDSQR_material(rho(im),massfrac_parm, &
               internal_energy, &
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

            do im_opp=1,num_materials
             if (im_opp.ne.im) then
              call get_iten(im,im_opp,iten)
              if ((fort_material_type_interface(iten).eq.0).or. &
                  (fort_material_type_interface(iten).eq.999)) then
               if (localLS(im_opp).ge.-incomp_thickness*DXMAXLS) then
                pres(im)=zero
                one_over_c2(im)=zero
                one_over_c(im)=zero
               else if (localLS(im_opp).le.-incomp_thickness*DXMAXLS) then
                !do nothing
               else
                print *,"localLS(im_opp) invalid:",im,localLS(im_opp)
                stop
               endif
              else if ((fort_material_type_interface(iten).ge.1).and. &
                       (fort_material_type_interface(iten).le.MAX_NUM_EOS)) then
               !do nothing
              else
               print *,"fort_material_type_interface(iten) invalid: ", &
                 iten,fort_material_type_interface
               stop
              endif
             endif !im_opp<>im
            enddo !im_opp=1,num_materials

           else if ((fort_material_type(im).eq.0).or. &
                    (vfrac(im).eq.zero)) then
            pres(im)=zero
            one_over_c2(im)=zero
            one_over_c(im)=zero
           else
            print *,"material type or vfrac invalid: ",im, &
               fort_material_type(im),vfrac(im)
            stop
           endif
   
          else if (is_rigid(im).eq.1) then

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

         enddo ! im=1,num_materials

         infinite_weight=0
         im_weight=0

         if (imcrit.eq.0) then ! rigid material(s) dominate
          im_weight=imcrit
          infinite_weight=1

          ! fluid material dominates
         else if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then
  
          if (is_rigid(imcrit).eq.0) then

           do im=1,num_materials

            local_infinite_weight=0

            if (is_rigid(im).eq.1) then
             ! do nothing
            else if (is_rigid(im).eq.0) then
             if ((vfrac(im).gt.zero).and. &
                 (vfrac(im).le.one+EPS1)) then
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
              print *,"vfrac invalid: ",im,vfrac(im)
              stop
             endif
             if (im_weight.eq.0) then
              im_weight=im
             else if ((im_weight.ge.1).and.(im_weight.le.num_materials)) then
              if ((vfrac_weight(im).gt.vfrac_weight(im_weight)).or. &
                  (local_infinite_weight.eq.1)) then
               im_weight=im
              endif
             else
              print *,"im_weight invalid"
              stop
             endif
            else
             print *,"is_rigid(im) invalid"
             stop
            endif

           enddo ! im=1,num_materials

          else
           print *,"is_rigid(imcrit) invalid"
           stop
          endif

         else
          print *,"imcrit invalid"
          stop
         endif
       
         if (imcrit.eq.0) then ! is_rigid material(s) dominate.

          csnd(D_DECL(i,j,k),1)=zero  ! coeff
          csnd(D_DECL(i,j,k),2)=zero  ! padvect

          if (im_weight.eq.0) then
           ! do nothing
          else
           print *,"expecting im_weight==0"
           stop
          endif

         else if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then

          if (is_rigid(imcrit).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(imcrit) invalid"
           stop
          endif

          if (is_rigid(im_weight).eq.0) then
           ! do nothing
          else
           print *,"is_rigid(im_weight).ne.0"
           stop
          endif

          if (infinite_weight.eq.1) then ! sound speed=infinity

           csnd(D_DECL(i,j,k),1)=zero  ! coeff
           csnd(D_DECL(i,j,k),2)=zero  ! padvect

          else if (infinite_weight.eq.0) then

           if ((im_weight.lt.1).or.(im_weight.gt.num_materials)) then
            print *,"im_weight invalid"
            stop
           endif
           if (vfrac_weight(im_weight).gt.zero) then
            ! do nothing
           else
            print *,"vfrac_weight(im_weight) invalid"
            stop
           endif

            ! im_weight is a compressible material
           if ((fort_material_type(im_weight).ge.1).and. &
               (fort_material_type(im_weight).le.MAX_NUM_EOS)) then

            if (rho(im_weight).gt.zero) then
             ! do nothing
            else
             print *,"rho(im_weight) invalid"
             stop
            endif

             ! padvect (zero for the incompressible materials)
            csnd(D_DECL(i,j,k),2)=pres(im_weight)

            csound_hold=one_over_c2(im_weight)/ &
              (compressible_dt_factor(im_weight)*rho(im_weight)*dt*dt)

             ! coeff
            if (csound_hold.ge.zero) then
             csnd(D_DECL(i,j,k),1)=csound_hold
            else
             print *,"csound_hold invalid"
             stop
            endif

            ! im_weight is an incompressible material
           else if (fort_material_type(im_weight).eq.0) then

            csnd(D_DECL(i,j,k),1)=zero ! coeff
            csnd(D_DECL(i,j,k),2)=zero ! padvect

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
         print *,"project_option invalid advective pressure 2: ", &
           project_option
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
      end subroutine fort_advective_pressure

       ! called from NavierStokes::ADVECT_DIV which is declared in MacProj.cpp
      subroutine fort_update_div( &
        xlo,dx, &
        dt, &
        vol,DIMS(vol), &
        csound,DIMS(csound), &
        mdot,DIMS(mdot), &
        pnew,DIMS(pnew), &
        divnew,DIMS(divnew), &
        tilelo,tilehi, &
        fablo,fabhi,bfact) &
      bind(c,name='fort_update_div')

      use global_utility_module
      use probf90_module

      IMPLICIT NONE


      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dt

      integer, INTENT(in) :: DIMDEC(vol)
      integer, INTENT(in) :: DIMDEC(csound)
      integer, INTENT(in) :: DIMDEC(mdot)
      integer, INTENT(in) :: DIMDEC(pnew)
      integer, INTENT(in) :: DIMDEC(divnew)
      real(amrex_real), INTENT(in),target :: vol(DIMV(vol))
      real(amrex_real), pointer :: vol_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: csound(DIMV(csound),2) 
      real(amrex_real), pointer :: csound_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: mdot(DIMV(mdot)) 
      real(amrex_real), pointer :: mdot_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: pnew(DIMV(pnew)) 
      real(amrex_real), pointer :: pnew_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out),target :: divnew(DIMV(divnew)) 
      real(amrex_real), pointer :: divnew_ptr(D_DECL(:,:,:))

      integer i,j,k
      real(amrex_real) coeff_hold,compress_term,mdot_term
      integer sound_comp
      real(amrex_real) local_volume

      if (bfact.lt.1) then
       print *,"bfact invalid165"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      vol_ptr=>vol
      call checkbound_array1(fablo,fabhi,vol_ptr,0,-1)
      csound_ptr=>csound
      call checkbound_array(fablo,fabhi,csound_ptr,0,-1)
      mdot_ptr=>mdot
      call checkbound_array1(fablo,fabhi,mdot_ptr,0,-1)
      pnew_ptr=>pnew
      call checkbound_array1(fablo,fabhi,pnew_ptr,1,-1)
      divnew_ptr=>divnew
      call checkbound_array1(fablo,fabhi,divnew_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       local_volume=vol(D_DECL(i,j,k))
       if (local_volume.gt.zero) then
        ! do nothing
       else
        print *,"local_volume must be positive"
        stop
       endif

       sound_comp=1
       coeff_hold=csound(D_DECL(i,j,k),sound_comp) ! 1/(rho c^2 dt^2)

       if (coeff_hold.gt.zero) then ! compressible
        compress_term=-dt*(pnew(D_DECL(i,j,k))- &
          csound(D_DECL(i,j,k),sound_comp+1))*coeff_hold
       else if (coeff_hold.eq.zero) then ! incompressible
        compress_term=zero
       else
        print *,"coeff_hold invalid"
        stop
       endif

       mdot_term=mdot(D_DECL(i,j,k))

       divnew(D_DECL(i,j,k))=compress_term+mdot_term*dt/local_volume
         
      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine fort_update_div

       ! called from: NavierStokes2.cpp
      subroutine fort_cellgrid( &
       plot_grid_type, & !0=interp to nodes   1=data "as is"
       ncomp_tower, &
       tid, &
       bfact, &
       fabout,DIMS(fabout), &
       vislo,vishi, &
        ! x,u,pmg,den,T,Y1..Yn,mag vort,LS
       visual_ncomp, & 
       maskSEM,DIMS(maskSEM), &
       vel,DIMS(vel), &
       vof,DIMS(vof), &
       pres,DIMS(pres), &
       div,DIMS(div), &
       divdat,DIMS(divdat), &
       den,DIMS(den), &
       mom_den,DIMS(mom_den), &
       elastic, &
       DIMS(elastic), &
       refineden, &
       DIMS(refineden), &
       lsdist,DIMS(lsdist), &
       visc,DIMS(visc), &
       conduct,DIMS(conduct), &
       trace,DIMS(trace), &
       elasticforce, &
       DIMS(elasticforce), &
        !ux,vx,wx,uy,vy,wy,uz,vz,wz
       gradvelocity, &
       DIMS(gradvelocity), &
       towerfab, &
       DIMS(towerfab), &
       problo, &
       probhi, &
       dx, &
       tilelo,tilehi, &
       lo,hi, &
       level, &
       finest_level, &
       gridno, &
       visual_tessellate_vfrac, &
       rz_flag, &
       nparts, &
       nparts_def, &
       im_solid_map, &
       elastic_ncomp, &
       refineden_ncomp, &
       slice_data, &
       nslice, & ! number of nodes in space for the slice
       nstate_slice, & ! number of data items stored for each slice node.
       slice_dir, &
       xslice, &
       dxfinest, &
       do_plot, &
       do_slice, &
       visual_nddata_format) &
      bind(c,name='fort_cellgrid')

      use global_utility_module
      use probf90_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: plot_grid_type
      integer, INTENT(in) :: ncomp_tower
      integer, INTENT(in) :: tid
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: do_plot
      integer, INTENT(in) :: do_slice
      integer, INTENT(in) :: visual_nddata_format

        ! x,y,z,xvel,yvel,zvel,PMG,PEOS,DIV,DEN,TEMP,KE
        ! nstate_slice=SLICECOMP_NCOMP
        ! (value of material with LS>0)
        ! nslice=domhi-domlo+3
      integer, INTENT(in) :: nslice,nstate_slice,slice_dir
      real(amrex_real), INTENT(out) :: slice_data(nslice*nstate_slice)
      real(amrex_real), INTENT(in) :: xslice(SDIM)
 
      integer, INTENT(in) :: rz_flag
      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def) 
      integer, INTENT(in) :: elastic_ncomp
      integer, INTENT(in) :: refineden_ncomp
      integer, INTENT(in) :: visual_tessellate_vfrac
       ! x,u,pmg,den,temp,spec,mag vort,LS
      integer, INTENT(in) :: visual_ncomp
      integer, INTENT(in) :: vislo(SDIM), vishi(SDIM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: lo(SDIM), hi(SDIM)
      integer, INTENT(in) :: DIMDEC(fabout)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(vof)
      integer, INTENT(in) :: DIMDEC(pres)
      integer, INTENT(in) :: DIMDEC(div)
      integer, INTENT(in) :: DIMDEC(divdat)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(mom_den)
      integer, INTENT(in) :: DIMDEC(elastic)
      integer, INTENT(in) :: DIMDEC(refineden)
      integer, INTENT(in) :: DIMDEC(lsdist)
      integer, INTENT(in) :: DIMDEC(visc)
      integer, INTENT(in) :: DIMDEC(conduct)
      integer, INTENT(in) :: DIMDEC(trace)
      integer, INTENT(in) :: DIMDEC(elasticforce)
        !ux,vx,wx,uy,vy,wy,uz,vz,wz
      integer, INTENT(in) :: DIMDEC(gradvelocity)
      integer, INTENT(in) :: DIMDEC(towerfab)
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: gridno
       ! x,u,pmg,den,T,Y1..Yn,mag vort,LS
      real(amrex_real), INTENT(out), target :: fabout(DIMV(fabout),visual_ncomp) 
      real(amrex_real), pointer :: fabout_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: &
         vel(DIMV(vel),STATE_NCOMP_VEL+STATE_NCOMP_PRES)
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              vof(DIMV(vof),num_materials*ngeom_recon)
      real(amrex_real), pointer :: vof_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: pres(DIMV(pres))
      real(amrex_real), pointer :: pres_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: div(DIMV(div))
      real(amrex_real), pointer :: div_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: divdat(DIMV(divdat))
      real(amrex_real), pointer :: divdat_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: &
          den(DIMV(den),num_state_material*num_materials)
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              mom_den(DIMV(mom_den),num_materials)
      real(amrex_real), pointer :: mom_den_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: &
              elastic(DIMV(elastic),elastic_ncomp)
      real(amrex_real), pointer :: elastic_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: &
              refineden(DIMV(refineden),refineden_ncomp)
      real(amrex_real), pointer :: refineden_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: visc(DIMV(visc),num_materials)
      real(amrex_real), pointer :: visc_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              conduct(DIMV(conduct),num_materials)
      real(amrex_real), pointer :: conduct_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: trace(DIMV(trace),5*num_materials)
      real(amrex_real), pointer :: trace_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target ::  &
              elasticforce(DIMV(elasticforce),SDIM)
      real(amrex_real), pointer :: elasticforce_ptr(D_DECL(:,:,:),:)
        !ux,vx,wx,uy,vy,wy,uz,vz,wz
      real(amrex_real), INTENT(in), target :: &
              gradvelocity(DIMV(gradvelocity),AMREX_SPACEDIM_SQR)
      real(amrex_real), pointer :: gradvelocity_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout), target :: &
              towerfab(DIMV(towerfab),ncomp_tower)
      real(amrex_real), pointer :: towerfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              lsdist(DIMV(lsdist),(SDIM+1)*num_materials)
      real(amrex_real), pointer :: lsdist_ptr(D_DECL(:,:,:),:)
      real(amrex_real) xposnd(SDIM)
      real(amrex_real) xposndT(SDIM)
      real(amrex_real) machnd
      real(amrex_real) machcell
      real(amrex_real) velnd(STATE_NCOMP_VEL+STATE_NCOMP_PRES)
      real(amrex_real) velmat(SDIM)
      real(amrex_real) velmatT(SDIM)
      real(amrex_real) velcell(STATE_NCOMP_VEL+STATE_NCOMP_PRES)
      real(amrex_real) vofnd(num_materials)
      real(amrex_real) vofcell(num_materials)
      real(amrex_real) vofcell_denom
      real(amrex_real) presnd
      real(amrex_real) divnd
      real(amrex_real) divdatnd
      real(amrex_real) dennd(num_state_material*num_materials)
      real(amrex_real) dennd_merge(num_state_material)
      real(amrex_real) mom_dennd(num_materials)
      real(amrex_real) elasticnd(elastic_ncomp)
      real(amrex_real) refinedennd(refineden_ncomp)
      real(amrex_real) dencell(num_state_material*num_materials)
      real(amrex_real) dencell_merge(num_state_material)
      real(amrex_real) mom_dencell(num_materials)
      real(amrex_real) elasticcell(elastic_ncomp)
      real(amrex_real) refinedencell(refineden_ncomp)
      real(amrex_real) lsdistnd((SDIM+1)*num_materials)
      real(amrex_real) local_LS_data((SDIM+1)*num_materials)
      real(amrex_real) viscnd(num_materials)
      real(amrex_real) conductnd(num_materials)
      real(amrex_real) tracend(5*num_materials)
      real(amrex_real) elasticforcend(SDIM)
        !ux,vx,wx,uy,vy,wy,uz,vz,wz
      real(amrex_real) gradvelocitynd(AMREX_SPACEDIM_SQR)
      real(amrex_real) writend(num_materials*200)
      integer scomp,iw
      integer istate,idissolution
      integer im
      integer im_crit
      integer im_crit_SEM

      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: probhi(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: dxfinest(SDIM)
      real(amrex_real) dxelem

      character*3 levstr
      character*5 gridstr
      character*32 filename32

      integer i,j,k
      integer ii,jj
      integer dir,dir2
      integer n
      integer i1,j1,k1
      integer ic,jc,kc
      integer iBL,jBL,kBL
      integer iSEM,jSEM,kSEM
      integer ilocal,jlocal,klocal
      real(amrex_real) sumweight,localwt,int_xlo,int_xhi
      real(amrex_real) sumweightLS,localwtLS,wt_denom
      integer vofcomp
      integer i1lo,i1hi
      integer j1lo,j1hi
      integer k1lo,k1hi
      integer igridlo(3),igridhi(3)
      integer stenlo(3),stenhi(3)
      integer stenlo_index
      integer icrit(SDIM)
      integer icritlo(3),icrithi(3)
      integer icell_lo(3),icell_hi(3)
      integer icell(3)
      integer iproblo(3)
      real(amrex_real) xfablo(SDIM),xfabhi(SDIM)
      real(amrex_real) x1D
      real(amrex_real) xc(SDIM),theta(3)
      real(amrex_real) xstenlo(SDIM),xstenhi(SDIM)
      integer inbox
      real(amrex_real) psten(0:1,0:1)
      real(amrex_real) psten1D(0:1)
      integer DIMDEC(plt)

      real(amrex_real), dimension(D_DECL(:,:,:),:), allocatable, target :: plotfab
      real(amrex_real), pointer :: plotfab_ptr(D_DECL(:,:,:),:)

      real(amrex_real), dimension(D_DECL(:,:,:),:), allocatable :: reconfab
      real(amrex_real), dimension(D_DECL(:,:,:),:), allocatable :: reconfab_raster
      integer debug_slice
      real(amrex_real) denslice,tempslice,eslice,KEslice
      integer, parameter :: nhalf=3
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_corner(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten1D(-nhalf:nhalf)
      real(amrex_real) xsten1DL(-nhalf:nhalf)
      real(amrex_real) xsten1D_finest(-nhalf:nhalf)
      real(amrex_real) xsten_fablo(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_fabhi(-nhalf:nhalf,SDIM)
      real(amrex_real) dxleft,dxright,dxmin
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdata_raster(num_materials*ngeom_recon)
      integer, PARAMETER :: nmax=POLYGON_LIST_MAX ! in: fort_cellgrid
      integer, PARAMETER :: bfact_finest=2
      integer, PARAMETER :: tessellate_raster=3

      integer visual_ncell(SDIM)
      real(amrex_real) visual_dx(SDIM)
      real(amrex_real) xcrit(SDIM)
       ! x,u,pmg,den,temp,spec,mag vort,LS
      real(amrex_real) localfab(visual_ncomp)
       ! u,pmg,den,temp,spec,mag vort,LS
      real(amrex_real) vel_uniform(VISUALCOMP_NCOMP_INTERP)
       ! u,pmg,den,temp,spec,mag vort,LS
      real(amrex_real) SEM_value(VISUALCOMP_NCOMP_INTERP) 
      integer grid_type
      integer SEMhi(SDIM)
      real(amrex_real), dimension(D_DECL(:,:,:),:),allocatable :: SEMloc
      integer local_maskSEM
      real(amrex_real) local_data
      integer current_index
      real(amrex_real) massfrac_parm(num_species_var+1)
      integer ispec
      integer plot_sdim_macro
      integer nwrite3d
      integer nwrite2d

      plot_sdim_macro=2
      nwrite2d=PLOTCOMP_NCOMP
      plot_sdim_macro=3
      nwrite3d=PLOTCOMP_NCOMP
      if (nwrite3d-nwrite2d.eq.PLOTCOMP_DIFF) then
       ! do nothing
      else
       print *,"nwrite3d-nwrite2d invalid"
       print *,"nwrite3d: ",nwrite3d
       print *,"nwrite2d: ",nwrite2d
       print *,"PLOTCOMP_DIFF: ",PLOTCOMP_DIFF
       stop
      endif

      plot_sdim_macro=SDIM

      debug_slice=0

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid 39"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_cellgrid"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_cellgrid"
       stop
      endif
      if (elastic_ncomp.eq. &
          num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) then
       ! do nothing
      else
       print *,"elastic_ncomp invalid"
       stop
      endif
      if (refineden_ncomp.eq. &
          num_materials_compressible*ENUM_NUM_REFINE_DENSITY_TYPE) then
       ! do nothing
      else
       print *,"refineden_ncomp invalid"
       stop
      endif

 
        ! x,y,z,xvel,yvel,zvel,PMG,PEOS,DIV,den,Temp,KE
        ! nstate_slice
        ! (value of material with LS>0)
      if (nstate_slice.ne.SLICECOMP_NCOMP) then
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

      if (plot_grid_type.eq.0) then ! interp to nodes
       igridlo(3)=0
       igridhi(3)=0
       iproblo(3)=0
       do dir=1,SDIM
        igridlo(dir)=lo(dir)-1
        igridhi(dir)=hi(dir)+1
        iproblo(dir)=0
       enddo
      else if (plot_grid_type.eq.1) then ! data "as is"
       call growntilebox(tilelo,tilehi,lo,hi,igridlo,igridhi,1)
      else
       print *,"plot_grid_type invalid"
       stop
      endif

       ! input: igridlo,igridhi
       ! output: DIMS(plt)
      call box_to_dim( &
        DIMS(plt), &
        igridlo,igridhi)

      allocate(reconfab(DIMV(plt),num_materials*ngeom_recon))
      allocate(reconfab_raster(DIMV(plt),num_materials*ngeom_recon))

      if (plot_grid_type.eq.0) then ! interp to nodes.

         ! nstate_slice=SLICECOMP_NCOMP
         ! for the material with LS>0,
         ! x,y,z,u,v,w,pmg,peos,div,den,T,KE
       allocate(plotfab(DIMV(plt),nstate_slice))
 
       fabout_ptr=>fabout
       call checkbound_array(vislo,vishi,fabout_ptr,0,0)
       call checkbound_array(vislo,vishi,fabout_ptr,0,1)
       call checkbound_array(vislo,vishi,fabout_ptr,0,SDIM-1)
       ! x,u,pmg,den,T,Y1..Yn,mag vort,LS
       if (visual_ncomp.ne.VISUALCOMP_NCOMP) then
        print *,"visual_ncomp invalid" 
        stop
       endif

       plotfab_ptr=>plotfab
       call checkbound_array(lo,hi,plotfab_ptr,1,-1)

      else if (plot_grid_type.eq.1) then ! 1=data "as is"
       ! do nothing
      else
       print *,"plot_grid_type invalid"
       stop
      endif

      maskSEM_ptr=>maskSEM
      call checkbound_array1(lo,hi,maskSEM_ptr,0,-1)

      pres_ptr=>pres
      call checkbound_array1(lo,hi,pres_ptr,1,-1)
      div_ptr=>div
      call checkbound_array1(lo,hi,div_ptr,1,-1)
      divdat_ptr=>divdat
      call checkbound_array1(lo,hi,divdat_ptr,1,-1)
      den_ptr=>den
      call checkbound_array(lo,hi,den_ptr,1,-1)
      mom_den_ptr=>mom_den
      call checkbound_array(lo,hi,mom_den_ptr,1,-1)

      elastic_ptr=>elastic

      if ((num_materials_viscoelastic.ge.1).and. &
          (num_materials_viscoelastic.le.num_materials)) then
       call checkbound_array(lo,hi,elastic_ptr,1,-1)
      else if (num_materials_viscoelastic.eq.0) then
       ! do nothing
      else
       print *,"num_materials_viscoelastic invalid"
       stop
      endif

      refineden_ptr=>refineden

      if ((num_materials_compressible.ge.1).and. &
          (num_materials_compressible.le.num_materials)) then
       call checkbound_array(lo,hi,refineden_ptr,1,-1)
      else if (num_materials_compressible.eq.0) then
       ! do nothing
      else
       print *,"num_materials_compressible invalid"
       stop
      endif



      lsdist_ptr=>lsdist
      call checkbound_array(lo,hi,lsdist_ptr,1,-1)
      visc_ptr=>visc
      call checkbound_array(lo,hi,visc_ptr,1,-1)
      conduct_ptr=>conduct
      call checkbound_array(lo,hi,conduct_ptr,1,-1)
      trace_ptr=>trace
      call checkbound_array(lo,hi,trace_ptr,1,-1)
      elasticforce_ptr=>elasticforce
      call checkbound_array(lo,hi,elasticforce_ptr,1,-1)
        !ux,vx,wx,uy,vy,wy,uz,vz,wz
      gradvelocity_ptr=>gradvelocity
      call checkbound_array(lo,hi,gradvelocity_ptr,1,-1)
      towerfab_ptr=>towerfab
      call checkbound_array(lo,hi,towerfab_ptr,1,-1)
      vel_ptr=>vel
      call checkbound_array(lo,hi,vel_ptr,1,-1)
      vof_ptr=>vof
      call checkbound_array(lo,hi,vof_ptr,1,-1)

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
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
       print *,"rz_flag invalid in cellgrid"
       stop
      endif

      do k=igridlo(3),igridhi(3) 
      do j=igridlo(2),igridhi(2) 
      do i=igridlo(1),igridhi(1) 

       if (plot_grid_type.eq.0) then ! interp to nodes (after complement)
         ! iproblo=0
        call gridsten(xsten,problo,i,j,k,iproblo,bfact,dx,nhalf)
       else if (plot_grid_type.eq.1) then ! data "as is"
        call gridsten_level(xsten,i,j,k,level,nhalf)
       else
        print *,"plot_grid_type invalid"
        stop
       endif

       do dir=1,num_materials*ngeom_recon
        mofdata(dir)=vof(D_DECL(i,j,k),dir)
        mofdata_raster(dir)=vof(D_DECL(i,j,k),dir)
       enddo

       ! before (mofdata_raster): only fluids tessellate
       ! after  (mofdata_raster): fluids and solids tessellate (rasterized)
       ! EPS2
       call multi_get_volume_tessellate( &
        tessellate_raster, & !=3
        bfact, &
        dx,xsten,nhalf, &
        mofdata_raster, &
        geom_xtetlist(1,1,1,tid+1), &
        nmax, &
        nmax, &
        SDIM)

       if ((visual_tessellate_vfrac.eq.1).or. & 
           (visual_tessellate_vfrac.eq.3)) then ! rasterized solid.
         ! before (mofdata): only fluids tessellate
         ! after  (mofdata): fluids and solids tessellate
         ! EPS2
        call multi_get_volume_tessellate( &
         visual_tessellate_vfrac, &
         bfact, &
         dx,xsten,nhalf, &
         mofdata, &
         geom_xtetlist(1,1,1,tid+1), &
         nmax, &
         nmax, &
         SDIM)
       else if (visual_tessellate_vfrac.eq.0) then
        ! do nothing
       else
        print *,"visual_tessellate_vfrac invalid"
        stop
       endif

       do dir=1,num_materials*ngeom_recon
        reconfab(D_DECL(i,j,k),dir)=mofdata(dir)
        reconfab_raster(D_DECL(i,j,k),dir)=mofdata_raster(dir)
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

       if (plot_grid_type.eq.0) then !interp to nodes.

        igridlo(3)=0
        igridhi(3)=0

        i1lo=0
        i1hi=1

        j1lo=0
        j1hi=1

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

       else if (plot_grid_type.eq.1) then ! data "as is"
        i1lo=0
        i1hi=0

        j1lo=0
        j1hi=0

        k1lo=0
        k1hi=0

        call growntilebox(tilelo,tilehi,lo,hi,igridlo,igridhi,1)
       else
        print *,"plot_grid_type invalid"
        stop
       endif

       write(filename32,'(A14,A10,A3,A5)') &
         './temptecplot/','tempnddata',levstr,gridstr

        !visual_nddata_format.eq.0 -> tecplot nodes
       if ((visual_nddata_format.eq.0).and. &
           (plot_grid_type.eq.0)) then ! 0="interp to nodes"

        print *,"filename32 ",filename32

        open(unit=11,file=filename32)
        do dir=1,SDIM
         write(11,*) lo(dir),hi(dir)
        enddo

        !visual_nddata_format.eq.1 -> plt file cells
        !visual_nddata_format.eq.2 -> tecplot cells (piecewise const interp)
       else if ((visual_nddata_format.ne.0).or. & 
                (plot_grid_type.eq.1)) then ! data "as is"

        if (1.eq.0) then
         print *,"ignoring box for nddata*.tec output (lo123,hi123): ", &
          igridlo(1),igridlo(2),igridlo(3), &
          igridhi(1),igridhi(2),igridhi(3)
        endif

       else
        print *,"visual_nddata_format or plot_grid_type invalid"
        stop
       endif

        ! the order k,j,i is IMPORTANT.
       do k=igridlo(3),igridhi(3) 
       do j=igridlo(2),igridhi(2) 
       do i=igridlo(1),igridhi(1) 

        if (plot_grid_type.eq.0) then ! interp to nodes
          ! iproblo=0
         call gridstenND(xstenND,problo,i,j,k,iproblo,bfact,dx,nhalf)
        else if (plot_grid_type.eq.1) then ! data "as is"
         call gridsten_level(xstenND,i,j,k,level,nhalf)
        else
         print *,"plot_grid_type invalid"
         stop
        endif

        do dir=1,SDIM
         dxleft=xstenND(0,dir)-xstenND(-1,dir)
         dxright=xstenND(1,dir)-xstenND(0,dir)

         if (bfact.eq.1) then
          if (abs(dxleft-dxright).le.EPS2*dx(dir)) then
           ! do nothing 
          else
           print *,"xstenND invalid (2468) ",dxleft,dxright
           stop
          endif
         else if (bfact.gt.1) then
          if ((dxleft.gt.zero).and.(dxright.gt.zero)) then
           ! do nothing
          else
           print *,"xstenND invalid (2475) ",dxleft,dxright
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
        if (visual_RT_transform.eq.1) then ! defined in PROBCOMMON.F90
         call RT_transform(xposnd,xposndT)
        endif

        sumweight=zero
        sumweightLS=zero

        do dir=1,STATE_NCOMP_VEL+STATE_NCOMP_PRES
         velnd(dir)=zero
        enddo
        do dir=1,num_state_material*num_materials
         dennd(dir)=zero
        enddo
        do dir=1,num_state_material
         dennd_merge(dir)=zero
        enddo
        do dir=1,num_materials
         mom_dennd(dir)=zero
        enddo
        do dir=1,elastic_ncomp
         elasticnd(dir)=zero
        enddo
        do dir=1,refineden_ncomp
         refinedennd(dir)=zero
        enddo
        presnd=zero
        divnd=zero
        divdatnd=zero
        machnd=zero
        do dir=1,num_materials
         vofnd(dir)=zero
         viscnd(dir)=zero
         conductnd(dir)=zero
        enddo
        do dir=1,num_materials*(SDIM+1)
         lsdistnd(dir)=zero
        enddo
        do dir=1,5*num_materials
         tracend(dir)=zero
        enddo
        do dir=1,SDIM
         elasticforcend(dir)=zero
        enddo
        do dir=1,AMREX_SPACEDIM_SQR
         gradvelocitynd(dir)=zero
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

        do k1=k1lo,k1hi
        do j1=j1lo,j1hi
        do i1=i1lo,i1hi

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

         if (plot_grid_type.eq.0) then !interp to nodes.
          ! iproblo=0
          call gridsten(xsten,problo,i-i1,j-j1,k-k1,iproblo,bfact,dx,nhalf)
         else if (plot_grid_type.eq.1) then ! data "as is"
          call gridsten_level(xsten,i,j,k,level,nhalf)
         else
          print *,"plot_grid_type invalid"
          stop
         endif
         
         localwt=one  ! area of intersection of node CV with given cell CV
         localwtLS=one

         if (plot_grid_type.eq.0) then !interp to nodes

          do dir=1,SDIM
           dxleft=xsten(0,dir)-xsten(-1,dir)
           dxright=xsten(1,dir)-xsten(0,dir)
          
           if (bfact.eq.1) then
            if (abs(dxleft-dxright).le.EPS2*dx(dir)) then
             ! do nothing
            else
             print *,"xsten invalid (2597) : ",dxleft,dxright
             stop
            endif
           else if (bfact.gt.1) then
            if ((dxleft.gt.zero).and.(dxright.gt.zero)) then
             ! do nothing
            else
             print *,"dxleft or dxright <=0.0 (2604) :",dxleft,dxright
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

          if (localwt.gt.zero) then
           ! do nothing
          else
           print *,"localwt invalid"
           stop
          endif
          if (localwtLS.gt.zero) then
           ! do nothing
          else
           print *,"localwtLS invalid"
           stop
          endif
         else if (plot_grid_type.eq.1) then ! data "as is"
          ! do nothing
         else
          print *,"plot_grid_type invalid"
          stop
         endif

         do dir=1,STATE_NCOMP_VEL+STATE_NCOMP_PRES
          velcell(dir)=vel(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         do dir=1,STATE_NCOMP_VEL+STATE_NCOMP_PRES
          velnd(dir)=velnd(dir)+localwt*velcell(dir)
         enddo

         do dir=1,num_materials
          vofcomp=(dir-1)*ngeom_recon+1
          vofcell(dir)=reconfab_raster(D_DECL(i-i1,j-j1,k-k1),vofcomp)
         enddo

         do dir=1,num_state_material*num_materials
          dencell(dir)=den(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo

         vofcell_denom=zero
         do dir=1,num_materials
          vofcell_denom=vofcell_denom+vofcell(dir)
         enddo 
         if ((vofcell_denom.gt.half).and. &
             (vofcell_denom.lt.one+half)) then
          ! do nothing
         else
          print *,"vofcell_denom invalid"
          stop
         endif

         do dir=1,num_state_material
          dencell_merge(dir)=zero
          do dir2=1,num_materials
           current_index=(dir2-1)*num_state_material+dir
           dencell_merge(dir)=dencell_merge(dir)+ &
             vofcell(dir2)*dencell(current_index)
          enddo !dir2=1..num_materials
          dencell_merge(dir)=dencell_merge(dir)/vofcell_denom
         enddo !dir=1..num_state_material

         do dir=1,num_materials
          mom_dencell(dir)=mom_den(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo

         do dir=1,num_state_material*num_materials
          dennd(dir)=dennd(dir)+localwt*dencell(dir)
         enddo
         do dir=1,num_state_material
          dennd_merge(dir)=dennd_merge(dir)+localwt*dencell_merge(dir)
         enddo
         do dir=1,num_materials
          mom_dennd(dir)=mom_dennd(dir)+localwt*mom_dencell(dir)
         enddo

         do dir=1,elastic_ncomp
          elasticcell(dir)=elastic(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         do dir=1,refineden_ncomp
          refinedencell(dir)=refineden(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo

         do dir=1,elastic_ncomp
          elasticnd(dir)=elasticnd(dir)+localwt*elasticcell(dir)
         enddo
         do dir=1,refineden_ncomp
          refinedennd(dir)=refinedennd(dir)+localwt*refinedencell(dir)
         enddo

         presnd=presnd+localwt*pres(D_DECL(i-i1,j-j1,k-k1))
         divnd=divnd+localwt*div(D_DECL(i-i1,j-j1,k-k1))
         divdatnd=divdatnd+localwt*divdat(D_DECL(i-i1,j-j1,k-k1))

         do dir=1,num_materials
          vofcomp=(dir-1)*ngeom_recon+1
          vofcell(dir)=reconfab(D_DECL(i-i1,j-j1,k-k1),vofcomp)
         enddo

         do dir=1,num_materials
          vofnd(dir)=vofnd(dir)+localwt*vofcell(dir)
          viscnd(dir)=viscnd(dir)+ &
            localwt*visc(D_DECL(i-i1,j-j1,k-k1),dir)
          conductnd(dir)=conductnd(dir)+ &
            localwt*conduct(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         do dir=1,num_materials*(SDIM+1)
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

         do dir=1,5*num_materials
          tracend(dir)=tracend(dir)+ &
            localwt*trace(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         do dir=1,SDIM
          elasticforcend(dir)=elasticforcend(dir)+ &
            localwt*elasticforce(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         !ux,vx,wx,uy,vy,wy,uz,vz,wz
         do dir=1,AMREX_SPACEDIM_SQR
          gradvelocitynd(dir)=gradvelocitynd(dir)+ &
            localwt*gradvelocity(D_DECL(i-i1,j-j1,k-k1),dir)
         enddo
         call get_mach_number(visual_tessellate_vfrac, &
           velcell,dencell,vofcell,machcell)

         machnd=machnd+localwt*machcell

         sumweight=sumweight+localwt
         sumweightLS=sumweightLS+localwtLS
        enddo ! k1
        enddo ! j1
        enddo ! i1

        if (sumweight.gt.zero) then
         ! do nothing
        else
         print *,"sumweight invalid"
         stop
        endif
        if (sumweightLS.gt.zero) then
         ! do nothing
        else
         print *,"sumweightLS invalid"
         stop
        endif

        do dir=1,STATE_NCOMP_VEL+STATE_NCOMP_PRES
         velnd(dir)=velnd(dir)/sumweight
        enddo
        do dir=1,5*num_materials
         tracend(dir)=tracend(dir)/sumweight
        enddo
        do dir=1,SDIM
         elasticforcend(dir)=elasticforcend(dir)/sumweight
        enddo
        do dir=1,AMREX_SPACEDIM_SQR
         gradvelocitynd(dir)=gradvelocitynd(dir)/sumweight
        enddo

        do dir=1,num_state_material*num_materials
         dennd(dir)=dennd(dir)/sumweight
        enddo
        do dir=1,num_state_material
         dennd_merge(dir)=dennd_merge(dir)/sumweight
        enddo
        do dir=1,num_materials
         mom_dennd(dir)=mom_dennd(dir)/sumweight
        enddo

        do dir=1,elastic_ncomp
         elasticnd(dir)=elasticnd(dir)/sumweight
        enddo
        do dir=1,refineden_ncomp
         refinedennd(dir)=refinedennd(dir)/sumweight
        enddo

        presnd=presnd/sumweight
        divnd=divnd/sumweight
        divdatnd=divdatnd/sumweight
        machnd=machnd/sumweight

        do dir=1,num_materials
         vofnd(dir)=vofnd(dir)/sumweight
         viscnd(dir)=viscnd(dir)/sumweight
         conductnd(dir)=conductnd(dir)/sumweight
        enddo
        do dir=1,num_materials*(SDIM+1)
         lsdistnd(dir)=lsdistnd(dir)/sumweightLS
        enddo

        if (plot_grid_type.eq.0) then !interp to nodes

         if (bfact.eq.1) then
          ! do nothing
  
          ! if high order, and in the bulk, we use SEM interpolation
          ! to get values at the nodes.
         else if (bfact.gt.1) then
          local_maskSEM= &
           NINT(maskSEM(D_DECL(icell(1),icell(2),icell(3))))
          if ((local_maskSEM.ge.1).and. &
              (local_maskSEM.le.num_materials)) then

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
            stenlo_index=NINT( (xcomp_of_xphys(dir2-1,xsten(0,dir2))- &
                                problo(dir2))/dxelem-half )
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

             if ((xcrit(dir2).ge.-EPS3*dx(dir2)).and. &
                 (xcrit(dir2).le.(EPS3+bfact)*dx(dir2))) then
              ! do nothing
             else
              print *,"xcrit out of bounds fort_cellgrid"
              print *,"dir2,xcrit,xsten_corner ",dir2,xcrit(dir2), &
               xsten_corner(0,dir2)
              stop
             endif
            else
             print *,"xsten_corner invalid" 
             stop
            endif

           enddo ! dir2=1..sdim

            ! x,u,pmg,den,T,Y1..Yn,mag vort,LS
           if (visual_ncomp.ne.VISUALCOMP_NCOMP) then
            print *,"visual_ncomp invalid" 
            stop
           endif
           if (visual_ncomp.ne.VISUALCOMP_NCOMP_INTERP+AMREX_SPACEDIM) then
            print *,"visual_ncomp invalid" 
            stop
           endif
           do n=1,VISUALCOMP_NCOMP_INTERP
            SEM_value(n)=zero
           enddo
           grid_type=-1  ! ggg  (Gauss in all directions)
           do dir2=1,SDIM
            SEMhi(dir2)=bfact-1
           enddo

           allocate(SEMloc(D_DECL(0:SEMhi(1),0:SEMhi(2),0:SEMhi(SDIM)), &
                           VISUALCOMP_NCOMP_INTERP))

           do kSEM=stenlo(3),stenhi(3)
           do jSEM=stenlo(2),stenhi(2)
           do iSEM=stenlo(1),stenhi(1)
            ilocal=iSEM-stenlo(1)
            jlocal=jSEM-stenlo(2)
            klocal=kSEM-stenlo(3)
            do im=1,num_materials
             local_LS_data(im)=lsdist(D_DECL(iSEM,jSEM,kSEM),im)
            enddo
            call get_primary_material(local_LS_data,im_crit_SEM)

             ! velocity 
            do dir=1,SDIM
             local_data=vel(D_DECL(iSEM,jSEM,kSEM),dir)
             if (abs(local_data).lt.1.0D+20) then
              SEMloc(D_DECL(ilocal,jlocal,klocal),VISUALCOMP_U-SDIM+dir)= &
                  local_data
             else
              print *,"abs(local_data) overflow1"
              stop
             endif
            enddo ! dir=1..sdim
    
             ! pressure
            local_data=vel(D_DECL(iSEM,jSEM,kSEM),STATE_NCOMP_VEL+1)
            if (abs(local_data).lt.1.0D+20) then
             SEMloc(D_DECL(ilocal,jlocal,klocal),VISUALCOMP_PMG+1-SDIM)= &
                  local_data
            else
             print *,"abs(local_data) overflow1"
             stop
            endif

            do dir=1,num_state_base+num_species_var
             local_data=den(D_DECL(iSEM,jSEM,kSEM), &
               (im_crit_SEM-1)*num_state_material+dir)
             if (abs(local_data).lt.1.0D+20) then
              SEMloc(D_DECL(ilocal,jlocal,klocal), &
                 VISUALCOMP_DEN+dir-SDIM)=local_data
             else
              print *,"abs(local_data) overflow1"
              stop
             endif
            enddo ! dir=1..num_state_base+num_species_var

              ! mag vorticity component of a derived data structure.
            local_data=trace(D_DECL(iSEM,jSEM,kSEM),(local_maskSEM-1)*5+5)
            if (abs(local_data).lt.1.0D+20) then
             SEMloc(D_DECL(ilocal,jlocal,klocal),VISUALCOMP_VORTMAG+1-SDIM)= &
                 local_data
            else
             print *,"abs(local_data) overflow2"
             stop
            endif

             ! we interpolate the LS using SEM, but we discard the interpolated
             ! value.
            do im=1,num_materials
             if (abs(local_LS_data(im)).lt.1.0D+20) then
              SEMloc(D_DECL(ilocal,jlocal,klocal),VISUALCOMP_LS+im-SDIM)= &
                local_LS_data(im)
             else
              print *,"abs(local_data) overflow25"
              stop
             endif
            enddo ! im=1..num_materials
             
           enddo ! kSEM
           enddo ! jSEM
           enddo ! iSEM

           call SEM_INTERP_ELEMENT( &
            VISUALCOMP_NCOMP_INTERP, &
            bfact,grid_type, &
            SEMhi,dx,xcrit,SEMloc,SEM_value)

           deallocate(SEMloc)

            ! WE ARE IN THE BULK HERE, CAN USE SEM INTERPOLATION
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
           presnd=SEM_value(VISUALCOMP_PMG+1-SDIM)
           if (visual_ncomp.ne.VISUALCOMP_NCOMP) then
            print *,"incorrect visual_ncomp"
            stop
           endif
           if (visual_ncomp-SDIM.ne.VISUALCOMP_NCOMP_INTERP) then
            print *,"incorrect visual_ncomp"
            stop
           endif
           if (abs(SEM_value(VISUALCOMP_VORTMAG+1-SDIM)).lt.1.0D+20) then
            tracend((local_maskSEM-1)*5+5)= &
                 SEM_value(VISUALCOMP_VORTMAG+1-SDIM)
           else 
            print *,"abs(SEM_value(sdim+1)) overflow"
            print *,"i,j,k,level,finest_level ", &
             i,j,k,level,finest_level
            print *,"stenlo ",stenlo(1),stenlo(2),stenlo(3)
            print *,"stenhi ",stenhi(1),stenhi(2),stenhi(3)
            print *,"lo ",lo(1),lo(2),lo(SDIM)
            print *,"hi ",hi(1),hi(2),hi(SDIM)
            print *,"SEM_value(VISUALCOMP_VORTMAG+1-SDIM) ", &
                  SEM_value(VISUALCOMP_VORTMAG+1-SDIM)
            stop
           endif
           do dir=1,num_state_base+num_species_var
            dennd(dir)=SEM_value(VISUALCOMP_DEN+dir-SDIM)
           enddo
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

        else if (plot_grid_type.eq.1) then ! data "as is"
         ! do nothing
        else
         print *,"plot_grid_type invalid"
         stop
        endif

        if (probtype.eq.5700) then
         if (num_materials.eq.3) then
          if (lsdistnd(3).ge.zero) then
           if (lsdistnd(1).ge.-lsdistnd(3)) then
            lsdistnd(1)=-lsdistnd(3)
           endif
           if (lsdistnd(2).ge.-lsdistnd(3)) then
            lsdistnd(2)=-lsdistnd(3)
           endif
          endif
         else
          print *,"num_materials invalid for 5700 option"
          stop
         endif
        endif  ! probtype=5700

        scomp=0
        do iw=1,SDIM
         writend(scomp+iw)=xposndT(iw)
        enddo

        scomp=scomp+SDIM

        if (scomp.eq.PLOTCOMP_XVEL) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_XVEL)"
         stop
        endif

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
  
        if (scomp.eq.PLOTCOMP_PRES) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_PRES)"
         stop
        endif

          ! this is pressure from the projection.
        writend(scomp+1)=velnd(STATE_NCOMP_VEL+1)
        scomp=scomp+1

          ! this is EOS pressure
        writend(scomp+1)=presnd
        scomp=scomp+1

        writend(scomp+1)=divnd
        scomp=scomp+1

        writend(scomp+1)=divdatnd
        scomp=scomp+1

        writend(scomp+1)=machnd
        scomp=scomp+1

        if (scomp.eq.PLOTCOMP_VFRAC) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_VFRAC)"
         stop
        endif

        do iw=1,num_materials
         writend(scomp+iw)=vofnd(iw)
        enddo
        scomp=scomp+num_materials

        do iw=1,num_materials*(SDIM+1)
         writend(scomp+iw)=lsdistnd(iw)
        enddo
        scomp=scomp+num_materials*(SDIM+1)

        if (num_state_base.eq.2) then
         call is_dissolution(idissolution)
        else
         print *,"num_state_base invalid"
         stop
        endif

        iw=0

        do im=1,num_materials
         do istate=1,num_state_material ! den,T
          iw=iw+1 
          writend(scomp+iw)=dennd(iw)
          if ((idissolution.eq.1).and.(istate.eq.2)) then
           writend(scomp+iw)=writend(scomp+iw)-one
          endif
         enddo ! istate=1..num_state_material
        enddo ! im=1..num_materials

        do istate=1,num_state_material ! den_merge,T_merge,scalar1..2.. MERGE
         iw=iw+1 
         writend(scomp+iw)=dennd_merge(istate)
         if ((idissolution.eq.1).and.(istate.eq.2)) then
          writend(scomp+iw)=writend(scomp+iw)-one
         endif
        enddo ! istate=1..num_state_material


        do im=1,num_materials
         iw=iw+1 
         writend(scomp+iw)=mom_dennd(im)
        enddo

        scomp=scomp+ &
              num_state_material*num_materials+ & !den
              num_state_material+ & !dennd_merge
              num_materials !MOMDEN

        if (scomp.eq.PLOTCOMP_CONFIG_TENSOR) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_CONFIG_TENSOR)"
         stop
        endif

        do iw=1,elastic_ncomp
         writend(scomp+iw)=elasticnd(iw) 
        enddo
        scomp=scomp+elastic_ncomp


        if (scomp.eq.PLOTCOMP_REFINEDEN) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_REFINEDEN)"
         stop
        endif

        do iw=1,refineden_ncomp
         writend(scomp+iw)=refinedennd(iw) 
        enddo
        scomp=scomp+refineden_ncomp


        if (scomp.eq.PLOTCOMP_VISC) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_VISC)"
         stop
        endif

        do iw=1,num_materials
         writend(scomp+iw)=viscnd(iw)
        enddo
        scomp=scomp+num_materials

        if (scomp.eq.PLOTCOMP_THERMALCOND) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_THERMALCOND)"
         stop
        endif

        do iw=1,num_materials
         writend(scomp+iw)=conductnd(iw)
        enddo
        scomp=scomp+num_materials

        do iw=1,5*num_materials
         writend(scomp+iw)=tracend(iw)
        enddo
        scomp=scomp+5*num_materials

        if (scomp.eq.PLOTCOMP_F_ELASTIC_X) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_F_ELASTIC_X)"
         stop
        endif

        do iw=1,SDIM
         writend(scomp+iw)=elasticforcend(iw)
        enddo
        scomp=scomp+SDIM

        if (scomp.eq.PLOTCOMP_GRAD_VELOCITY) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_GRAD_VELOCITY)"
         stop
        endif

         !ux,vx,wx,uy,vy,wy,uz,vz,wz
        do iw=1,AMREX_SPACEDIM_SQR
         writend(scomp+iw)=gradvelocitynd(iw)
        enddo
        scomp=scomp+AMREX_SPACEDIM_SQR

        if (scomp.eq.PLOTCOMP_NCOMP) then
         ! do nothing
        else
         print *,"(scomp.ne.PLOTCOMP_NCOMP)"
         stop
        endif

! pgf90 will automatically break up lines if they exceed 80 chars.
! a format must be specified.  e.g. '(D25.16)' or '(E25.16)'

        if (debug_slice.eq.1) then
         print *,"debug_slice: scomp= ",scomp
         do iw=1,SDIM
          print *,"debug_slice: iw,i,j,k,writend ",iw,i,j,k,writend(iw)
         enddo
        endif

        if (ncomp_tower.eq.PLOTCOMP_NCOMP) then
         ! do nothing
        else
         print *,"(ncomp_tower.ne.PLOTCOMP_NCOMP)"
         stop
        endif

        if (scomp.ne.ncomp_tower) then
         print *,"(scomp.ne.ncomp_tower)"
         stop
        endif

        do iw=1,scomp
         towerfab(D_DECL(i,j,k),iw)=writend(iw)
        enddo        

        if ((visual_nddata_format.eq.0).and. & ! tecplot format
            (plot_grid_type.eq.0)) then ! interp to nodes
  
         do iw=1,scomp
          if (iw.lt.scomp) then
!          write(11,'(D25.16)',ADVANCE="NO") writend(iw)
           write(11,'(E25.16)',ADVANCE="NO") writend(iw)
          else if (iw.eq.scomp) then
!          write(11,'(D25.16)') writend(iw)
           write(11,'(E25.16)') writend(iw)
          else
           print *,"iw invalid"
           stop
          endif
         enddo ! iw=1..scomp

        else if ((visual_nddata_format.ne.0).or. & !=1 plt file;=2tecplot cells
                 (plot_grid_type.eq.1)) then

         ! do nothing

        else
         print *,"visual_nddata_format or plot_grid_type invalid"
         stop
        endif

       enddo  ! i=igridlo(1),igridhi(1)
       enddo  ! j=igridlo(2),igridhi(2)
       enddo  ! k=igridlo(3),igridhi(3)  (main output loop to "AMR" grid)

       if ((visual_nddata_format.eq.0).and. &
           (plot_grid_type.eq.0)) then
        close(11)
       else if ((visual_nddata_format.ne.0).or. &
                (plot_grid_type.eq.1)) then
        ! do nothing
       else
        print *,"visual_nddata_format or plot_grid_type invalid"
        stop
       endif

      else if (do_plot.eq.0) then
       ! do nothing
      else
       print *,"do_plot invalid"
       stop
      endif

      if (plot_grid_type.eq.0) then ! interp to nodes

       if (do_slice.eq.1) then

        if ((slice_dir.ge.0).and.(slice_dir.lt.SDIM)) then

         igridlo(3)=0
         igridhi(3)=0
         do dir=1,SDIM
          igridlo(dir)=lo(dir)-1
          igridhi(dir)=hi(dir)+1
         enddo

         do k=igridlo(3),igridhi(3) 
         do j=igridlo(2),igridhi(2) 
         do i=igridlo(1),igridhi(1) 
            ! iproblo=0
          call gridsten(xsten,problo,i,j,k,iproblo,bfact,dx,nhalf)
          do dir=1,SDIM
           xposnd(dir)=xsten(0,dir)
          enddo

          do iw=1,SDIM
           plotfab(D_DECL(i,j,k),SLICECOMP_X+iw)=xposnd(iw)
          enddo

          do dir=1,num_materials
           lsdistnd(dir)=lsdist(D_DECL(i,j,k),dir)
          enddo

           ! primary material w.r.t. both fluids and solids.
          call get_primary_material(lsdistnd,im_crit)

          do dir=1,STATE_NCOMP_VEL+STATE_NCOMP_PRES
           velnd(dir)=vel(D_DECL(i,j,k),dir)
          enddo

          dir=STATE_NCOMP_VEL+1
          plotfab(D_DECL(i,j,k),SLICECOMP_PMG+1)=velnd(dir)

          presnd=pres(D_DECL(i,j,k))
          divnd=div(D_DECL(i,j,k))

          plotfab(D_DECL(i,j,k),SLICECOMP_PEOS+1)=presnd
          plotfab(D_DECL(i,j,k),SLICECOMP_DIV+1)=divnd

          do dir=1,num_state_material*num_materials
           dennd(dir)=den(D_DECL(i,j,k),dir)
          enddo

          denslice=dennd((im_crit-1)*num_state_material+ENUM_DENVAR+1)
          tempslice=dennd((im_crit-1)*num_state_material+ENUM_TEMPERATUREVAR+1)
          call init_massfrac_parm(denslice,massfrac_parm,im_crit)
          do ispec=1,num_species_var
           massfrac_parm(ispec)=dennd((im_crit-1)*num_state_material+2+ispec)
          enddo
          call INTERNAL_material(denslice,massfrac_parm, &
            tempslice,eslice, &
            fort_material_type(im_crit),im_crit)

          KEslice=denslice*eslice

          do dir=1,SDIM
           plotfab(D_DECL(i,j,k),SLICECOMP_U+dir)=velnd(dir)
           KEslice=KEslice+half*denslice*(velnd(dir)**2)
          enddo ! dir

          plotfab(D_DECL(i,j,k),SLICECOMP_DEN+1)=denslice
          plotfab(D_DECL(i,j,k),SLICECOMP_TEMP+1)=tempslice
          plotfab(D_DECL(i,j,k),SLICECOMP_KE+1)=KEslice
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
             enddo ! do n=1,nstate_slice

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

       do k=igridlo(3),igridhi(3) 
       do j=igridlo(2),igridhi(2) 
       do i=igridlo(1),igridhi(1) 
         ! iproblo=0
        call gridsten(xsten,problo,i,j,k,iproblo,bfact,dx,nhalf)
        icritlo(3)=0
        icrithi(3)=0
        do dir=1,SDIM
         icrit(dir)=NINT((xcomp_of_xphys(dir-1,xsten(0,dir))- &
                          problo(dir))/visual_dx(dir))
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
        do kc=icritlo(3),icrithi(3)
        do jc=icritlo(2),icrithi(2)
        do ic=icritlo(1),icrithi(1)
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

          do dir=SDIM+1,visual_ncomp
           localfab(dir)=zero
          enddo

          sumweight=zero
          do kBL=icell_lo(3),icell_hi(3)
          do jBL=icell_lo(2),icell_hi(2)
          do iBL=icell_lo(1),icell_hi(1)
           localwt=one
           do dir=1,SDIM
            if (xstenhi(dir)-xstenlo(dir).gt.zero) then
             theta(dir)=(xcrit(dir)-xstenlo(dir))/ &
              (xstenhi(dir)-xstenlo(dir))
             if ((theta(dir).ge.-EPS1).and. &
                 (theta(dir).le.one+EPS1)) then
              if (theta(dir).lt.zero) then
               theta(dir)=zero
              endif
              if (theta(dir).gt.one) then
               theta(dir)=one
              endif
              localwt=localwt*(one-theta(dir))
             else
              print *,"theta invalid: ",dir,theta(dir)
              stop
             endif
            else
             print *,"xstenhi or xstenlo bad: ",xstenlo(dir),xstenhi(dir)
             stop
            endif
           enddo ! dir=1..sdim
           do dir=1,num_materials
            lsdistnd(dir)=lsdist(D_DECL(iBL,jBL,kBL),dir)
           enddo
            ! primary material w.r.t. both fluids and solids.
           call get_primary_material(lsdistnd,im_crit)

           do dir=1,SDIM
            vel_uniform(VISUALCOMP_U+dir-SDIM)= &
                vel(D_DECL(iBL,jBL,kBL),dir)
           enddo

           vel_uniform(VISUALCOMP_PMG+1-SDIM)= &
              vel(D_DECL(iBL,jBL,kBL),STATE_NCOMP_VEL+1)

           do dir=1,num_state_base+num_species_var
            vel_uniform(VISUALCOMP_DEN+dir-SDIM)= &
              den(D_DECL(iBL,jBL,kBL),(im_crit-1)*num_state_material+dir)
           enddo

           current_index=SDIM
           do dir=1,VISUALCOMP_VORTMAG-SDIM
            localfab(SDIM+dir)=localfab(SDIM+dir)+ &
               localwt*vel_uniform(dir)
            current_index=current_index+1
           enddo
           if (current_index.eq.VISUALCOMP_VORTMAG) then
            ! do nothing
           else
            print *,"current_index invalid"
            stop
           endif

           ! 1. \dot{gamma}
           ! 2. Tr(A) if viscoelastic
           !    \dot{gamma} o.t.
           ! 3. Tr(A) (liquid viscosity - etaS)/etaP  if FENE-CR+Carreau
           !    Tr(A) if FENE-CR
           !    \dot{gamma} o.t.
           ! 4. (3) * f(A)  if viscoelastic
           !    \dot{gamma} o.t.
           ! 5. vorticity magnitude.
           localfab(current_index+1)=localfab(current_index+1)+ &
               localwt*trace(D_DECL(iBL,jBL,kBL),5)
           do im=1,num_materials
            localfab(VISUALCOMP_LS+im)= &
               localfab(VISUALCOMP_LS+im)+ &
               localwt*lsdistnd(im)
           enddo

           sumweight=sumweight+localwt
          enddo !kBL
          enddo !jBL
          enddo !iBL
          if (sumweight.gt.zero) then
           do dir=1,VISUALCOMP_NCOMP_INTERP
            localfab(SDIM+dir)=localfab(SDIM+dir)/sumweight
           enddo
          else
           print *,"sumweight invalid: ",sumweight
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
           enddo ! dir=1..sdim

           ! find stencil for surrounding spectral element.
           stenlo(3)=0
           stenhi(3)=0
           do dir2=1,SDIM
            dxelem=dx(dir2)*bfact
            stenlo_index=NINT( (xcomp_of_xphys(dir2-1,xsten(0,dir2))- &
                                problo(dir2))/dxelem-half )
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

            if ((xcrit(dir2).ge.-EPS3*dx(dir2)).and. &
                (xcrit(dir2).le.(EPS3+bfact)*dx(dir2))) then
             !do nothing
            else
             print *,"xcrit out of bounds fort_cellgrid: ", &
                 dir2,xcrit(dir2),bfact
             stop
            endif

           enddo ! dir2=1..sdim

           do n=1,VISUALCOMP_NCOMP_INTERP
            SEM_value(n)=zero
           enddo
           grid_type=-1  ! ggg  (Gauss in all directions)
           do dir2=1,SDIM
            SEMhi(dir2)=bfact-1
           enddo

           allocate(SEMloc(D_DECL(0:SEMhi(1),0:SEMhi(2),0:SEMhi(SDIM)), &
                           VISUALCOMP_NCOMP_INTERP))

           do kSEM=stenlo(3),stenhi(3)
           do jSEM=stenlo(2),stenhi(2)
           do iSEM=stenlo(1),stenhi(1)
            ilocal=iSEM-stenlo(1)
            jlocal=jSEM-stenlo(2)
            klocal=kSEM-stenlo(3)

            do dir=1,num_materials
             lsdistnd(dir)=lsdist(D_DECL(iSEM,jSEM,kSEM),dir)
            enddo
             ! primary material w.r.t. both fluids and solids.
            call get_primary_material(lsdistnd,im_crit)

            do dir=1,SDIM
             vel_uniform(VISUALCOMP_U+dir-SDIM)= &
                   vel(D_DECL(iSEM,jSEM,kSEM),dir)
            enddo
            vel_uniform(VISUALCOMP_PMG+1-SDIM)= &
              vel(D_DECL(iSEM,jSEM,kSEM),STATE_NCOMP_VEL+1)

            do dir=1,2+num_species_var
             vel_uniform(VISUALCOMP_DEN+dir-SDIM)= &
              den(D_DECL(iSEM,jSEM,kSEM),(im_crit-1)*num_state_material+dir)
            enddo

            do n=1,VISUALCOMP_VORTMAG-SDIM
             SEMloc(D_DECL(ilocal,jlocal,klocal),n)= &
              vel_uniform(n)
            enddo
            SEMloc(D_DECL(ilocal,jlocal,klocal),VISUALCOMP_VORTMAG+1-SDIM)= &
             trace(D_DECL(iSEM,jSEM,kSEM),5)
            do im=1,num_materials
             SEMloc(D_DECL(ilocal,jlocal,klocal),VISUALCOMP_LS+im-SDIM)= &
              lsdistnd(im)
            enddo
           enddo
           enddo
           enddo

           call SEM_INTERP_ELEMENT( &
            VISUALCOMP_NCOMP_INTERP, &
            bfact,grid_type, &
            SEMhi,dx,xcrit,SEMloc,SEM_value)

           deallocate(SEMloc)

           if (VISUALCOMP_U.eq.SDIM) then
            !do nothing
           else
            print *,"expecting VISUALCOMP_U=SDIM"
            stop
           endif
           if (VISUALCOMP_PMG.eq.2*SDIM) then
            !do nothing
           else
            print *,"expecting VISUALCOMP_PMG=2 * SDIM"
            stop
           endif
           if (VISUALCOMP_TEMP.eq.2*SDIM+2) then
            !do nothing
           else
            print *,"expecting VISUALCOMP_TEMP=2 * SDIM + 2"
            stop
           endif

            ! we discard the SEM interpolated values of lsdist, density,
            ! and species.
           do dir=SDIM+1,visual_ncomp

            if ( ((dir.ge.VISUALCOMP_U+1).and. &
                  (dir.le.VISUALCOMP_U+SDIM)).or. &
                 (dir.eq.VISUALCOMP_PMG+1).or. &
                 (dir.eq.VISUALCOMP_TEMP+1).or. &
                 (dir.eq.VISUALCOMP_VORTMAG+1)) then
             localfab(dir)=SEM_value(dir-SDIM) 
            endif

           enddo !dir=SDIM+1,visual_ncomp

          else
           print *,"bfact invalid155"
           stop
          endif

          do dir=1,visual_ncomp
           fabout(D_DECL(ic,jc,kc),dir)=localfab(dir)
          enddo

         else if (inbox.eq.0) then
          ! do nothing
         else
          print *,"inbox invalid: ",inbox
          stop
         endif
        enddo ! kc
        enddo ! jc
        enddo ! ic (visit all the uniform cells that overlap the "AMR" cells)

       enddo ! igridlo(3) .. igridhi(3)
       enddo ! igridlo(2) .. igridhi(2)
       enddo ! igridlo(1) .. igridhi(1)

       deallocate(plotfab)

      else if (plot_grid_type.eq.1) then ! data "as is"
       ! do nothing
      else
       print *,"plot_grid_type invalid"
       stop
      endif

      deallocate(reconfab)
      deallocate(reconfab_raster)
 
      return
      end subroutine fort_cellgrid


      subroutine fort_memstatus(procnum) &
      bind(c,name='fort_memstatus')

      IMPLICIT NONE

      integer procnum,i
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

      integer sysret

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
      end subroutine fort_memstatus

       ! fort_combinezones is called after calls to:
       !  fort_cellgrid
      subroutine fort_combinezones( &
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
       visual_revolve, &
       plotint, &
       nparts, &
       nparts_def, &
       im_solid_map) &
      bind(c,name='fort_combinezones')

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: nparts
      integer, INTENT(in) :: nparts_def
      integer, INTENT(in) :: im_solid_map(nparts_def)
      integer, INTENT(in) :: total_number_grids
      integer, INTENT(in) :: num_levels
      integer, INTENT(in) :: grids_per_level_array(num_levels)
      integer, INTENT(in) :: levels_array(total_number_grids)
      integer, INTENT(in) :: bfact_array(total_number_grids)
      integer, INTENT(in) :: gridno_array(total_number_grids)
      integer, INTENT(in) :: gridlo_array(total_number_grids*SDIM)
      integer, INTENT(in) :: gridhi_array(total_number_grids*SDIM)
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nsteps
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: visual_revolve
      integer, INTENT(in) :: plotint

      integer strandid
      integer nwrite


      character*3 levstr
      character*5 gridstr
      character*32 filename32

      character*6 stepstr
      character*16 newfilename16

      integer i,j,k,dir
      integer ilev,igrid
      integer lo(SDIM),hi(SDIM)
      integer hi_index_shift(3)

! Guibo

      character*80 Title,Zonename
      real(tecplot_real_short) ZONEMARKER,EOHMARKER
      integer*4 :: nzones_gb,iz_gb,ivar_gb
      integer*4, dimension(:,:), allocatable :: lo_gb,hi_gb
      integer bfact,testlev,testgridno
      integer testlo(SDIM),testhi(SDIM)

      ! define zone structure
      type zone_t
         real(tecplot_real), pointer :: var(:,:,:,:)
      end type zone_t
      type(zone_t), dimension(:), allocatable :: zone_gb

      integer plot_sdim
      integer plot_sdim_macro
      integer klo_plot,khi_plot
      integer add_sub_cells
      integer nwrite3d,nwrite2d

      plot_sdim_macro=2
      nwrite2d=PLOTCOMP_NCOMP
      plot_sdim_macro=3
      nwrite3d=PLOTCOMP_NCOMP
      if (nwrite3d-nwrite2d.eq.PLOTCOMP_DIFF) then
       ! do nothing
      else
       print *,"nwrite3d-nwrite2d invalid"
       print *,"nwrite3d: ",nwrite3d
       print *,"nwrite2d: ",nwrite2d
       print *,"PLOTCOMP_DIFF: ",PLOTCOMP_DIFF
       stop
      endif

      plot_sdim=SDIM
      plot_sdim_macro=SDIM

      if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
          (levelrz.eq.COORDSYS_CYLINDRICAL)) then
       if (visual_revolve.ne.0) then
        print *,"visual_revolve= ",visual_revolve
        print *,"visual_revolve invalid combine zones"
        stop
       endif
      else if (levelrz.eq.COORDSYS_RZ) then
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

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if ((nparts.lt.0).or.(nparts.gt.num_materials)) then
       print *,"nparts invalid fort_combinezones"
       stop
      endif
      if ((nparts_def.lt.1).or.(nparts_def.gt.num_materials)) then
       print *,"nparts_def invalid fort_combinezones"
       stop
      endif

      nwrite=PLOTCOMP_NCOMP

      if (num_levels.ne.finest_level+1) then
       print *,"num_levels invalid"
       stop
      endif

      if ((visual_revolve.ge.1).and.(visual_revolve.le.1024)) then

       if (levelrz.ne.COORDSYS_RZ) then
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
        visual_revolve, &
        plotint, &
        nparts, &
        nparts_def, &
        im_solid_map)
      else if (visual_revolve.eq.0) then

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
        write(filename32,'(A14,A10,A3,A5)') &
          './temptecplot/','tempnddata',levstr,gridstr
        open(unit=4,file=filename32)

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
       add_sub_cells=0
       call dumpstring_headers(plot_sdim,add_sub_cells)

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
        write(11) strandid-1  ! StrandID (this does not work)
        write(11) round_time(time) ! Solution time
        write(11) -1   ! Not used. Set to -1
        write(11) 0    ! Zone Type
        write(11) 0    ! Specify Var Location. 0 = Don't specify, 
                       ! all data is located at the nodes.
        write(11) 0    ! Are raw local 1-to-1 face neighbors supplied?
        write(11) 0    ! Number of miscellaneous user-defined  
                       !face neighbor connections

        hi_index_shift(1)=(hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1)+1
        hi_index_shift(2)=(hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1)+1
        if (plot_sdim.eq.3) then
         hi_index_shift(3)=(hi_gb(iz_gb,plot_sdim)-lo_gb(iz_gb,plot_sdim)+1)+1
        else if (plot_sdim.eq.2) then
         hi_index_shift(3)=1
        else
         print *,"plot_sdim invalid in zones_revolve_sanity"
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
 
        hi_index_shift(1)=(hi(1)-lo(1)+1)
        hi_index_shift(2)=(hi(2)-lo(2)+1)
        hi_index_shift(3)=khi_plot-klo_plot

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

        write(filename32,'(A14,A10,A3,A5)') &
          './temptecplot/','tempnddata',levstr,gridstr
        open(unit=4,file=filename32)
        print *,"filename32 ",filename32

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
         hi_index_shift(1)=(hi_gb(iz_gb,1)-lo_gb(iz_gb,1)+1)
         hi_index_shift(2)=(hi_gb(iz_gb,2)-lo_gb(iz_gb,2)+1)
         hi_index_shift(3)=khi_plot-klo_plot
         do k=0,hi_index_shift(3)
         do j=0,hi_index_shift(2)
         do i=0,hi_index_shift(1)
          write(11) zone_gb(iz_gb)%var(ivar_gb,i,j,k)
         enddo ! i
         enddo ! j
         enddo ! k
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
     
      else
       print *,"visual_revolve invalid"
       stop
      endif

      return
      end subroutine fort_combinezones

      subroutine fort_avgdown_copy( &
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
       lof,hif) &
      bind(c,name='fort_avgdown_copy')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: operation_flag
      integer, INTENT(in) :: dir
      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level_f
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact_f
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: xlo_fine(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: ncomp_flux
      integer, INTENT(in) :: ncomp_den
      integer, INTENT(in) :: ncomp_vel
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: DIMDEC(den_fine)
      integer, INTENT(in) :: DIMDEC(vel_fine)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: DIMDEC(fine_LS)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      integer, INTENT(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      integer mstenlo(3),mstenhi(3)
      real(amrex_real), INTENT(in),target :: fine_LS(DIMV(fine_LS),num_materials)
      real(amrex_real),pointer :: fine_LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: mask(DIMV(mask))
      real(amrex_real),pointer :: mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out),target :: crse(DIMV(crse),ncomp_flux)
      real(amrex_real),pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: fine(DIMV(fine),ncomp_flux)
      real(amrex_real),pointer :: fine_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: den_fine(DIMV(den_fine),ncomp_den)
      real(amrex_real),pointer :: den_fine_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: vel_fine(DIMV(vel_fine),ncomp_vel)
      real(amrex_real),pointer :: vel_fine_ptr(D_DECL(:,:,:),:)
      integer flochi(SDIM)
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer ilocal,jlocal,klocal
      integer istrip,jstrip,kstrip
      integer iface,jface,kface
      integer ii,jj,kk
      integer n
      integer dir2
      integer side
      integer grid_type
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_fine(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) crse_value(ncomp_flux)
      real(amrex_real) dxelem_f
      real(amrex_real) xcoarse(SDIM)
      integer finelo_index
      real(amrex_real), dimension(D_DECL(:,:,:),:),allocatable :: ffine
      integer testmask,testmask2
      integer elem_test
      integer im,imcrit
      real(amrex_real) fine_data
      integer dencomp
      real(amrex_real) LS_local(num_materials)
      integer mat_freq(num_materials)


      if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection
       if (ncomp_vel.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.num_materials*num_state_material) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.NFLUXSEM) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
               (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
               (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then 
       if (ncomp_vel.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.OP_PRESGRAD_MAC) then ! grad p
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
      else if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then !viscosity
       if (ncomp_vel.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_vel invalid visc"
        stop
       endif
       if (ncomp_den.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_den invalid visc"
        stop
       endif
       if (ncomp_flux.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_flux invalid visc"
        stop
       endif
      else
       print *,"operation_flag invalid visc:",operation_flag
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

      crse_ptr=>crse
      call checkbound_array(lo,hi,crse_ptr,0,dir)
      fine_ptr=>fine
      call checkbound_array(lof,hif,fine_ptr,0,dir)
      mask_ptr=>mask
      call checkbound_array1(lof,hif,mask_ptr,1,-1)
      fine_LS_ptr=>fine_LS
      call checkbound_array(lof,hif,fine_LS_ptr,1,-1)
      den_fine_ptr=>den_fine
      call checkbound_array(lof,hif,den_fine_ptr,1,-1)
      vel_fine_ptr=>vel_fine
      call checkbound_array(lof,hif,vel_fine_ptr,1,-1)

      grid_type=-1  ! ggg  (Gauss in all directions)

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
 
       do kc=growlo(3),growhi(3)
       do jc=growlo(2),growhi(2)
       do ic=growlo(1),growhi(1)

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
         finelo_index=NINT( (xcomp_of_xphys(dir2-1,xcoarse(dir2))- &
                             problo(dir2))/dxelem_f-half )
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
           if (abs(xcoarse(dir2)).lt.EPS4*dxf(dir2)) then
            if (dir2.ne.dir+1) then
             mstenlo(dir2)=stenlo(dir2)-1
            else
             print *,"xcoarse invalid in normal direction to face"
             stop
            endif
           endif
           if (abs(xcoarse(dir2)-bfact_f*dxf(dir2)).lt. &
              EPS4*dxf(dir2)) then
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

         if ((xcoarse(dir2).ge.-EPS3*dxf(dir2)).and. &
             (xcoarse(dir2).le.(EPS3+bfact_f)*dxf(dir2))) then
          !do nothing
         else
          print *,"xcoarse out of bounds avgdown_copy: ",dir2,xcoarse(dir2)
          stop
         endif

        enddo ! dir2=1..sdim

        do im=1,num_materials
         mat_freq(im)=0
        enddo
        do kfine=mstenlo(3),mstenhi(3) 
        do jfine=mstenlo(2),mstenhi(2) 
        do ifine=mstenlo(1),mstenhi(1) 
         do im=1,num_materials
          LS_local(im)=fine_LS(D_DECL(ifine,jfine,kfine),im)
         enddo
         call get_primary_material(LS_local,imcrit)
         if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then
          mat_freq(imcrit)=mat_freq(imcrit)+1
         else
          print *,"imcrit invalid"
          stop
         endif
        enddo 
        enddo 
        enddo 
        imcrit=0
        do im=1,num_materials
         if (mat_freq(im).gt.0) then
          if (imcrit.eq.0) then
           imcrit=im
          else if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then
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
        enddo ! im=1..num_materials

        ifine=mstenlo(1) 
        jfine=mstenlo(2) 
        kfine=mstenlo(SDIM) 
        testmask=NINT(mask(D_DECL(ifine,jfine,kfine)))
        do kfine=mstenlo(3),mstenhi(3) 
        do jfine=mstenlo(2),mstenhi(2) 
        do ifine=mstenlo(1),mstenhi(1) 
         testmask2=NINT(mask(D_DECL(ifine,jfine,kfine)))
         if (testmask2.ne.testmask) then
          testmask=0
         endif
        enddo 
        enddo 
        enddo 

        if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then

         do n=1,ncomp_flux
          crse_value(n)=zero
         enddo
         voltotal=zero

         if ((bfact_f.ge.2).and. &
             (enable_spectral.eq.1).and. &
             (testmask.ge.1).and. &
             (testmask.le.num_materials)) then

          do kfine=stenlo(3),stenhi(3)
          do jfine=stenlo(2),stenhi(2)
          do ifine=stenlo(1),stenhi(1)
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

            if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then ! viscosity
             if ((n.ge.1).and.(n.le.STATE_NCOMP_VEL)) then
              fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),n)
             else
              print *,"n invalid"
              stop
             endif
            else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection
             if ((n.ge.SEM_U+1).and.(n.le.SEM_W+1)) then ! velocity
              fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),n)
             else if (n.eq.SEM_T+1) then !temperature
              if (ENUM_DENVAR+1.eq.ENUM_TEMPERATUREVAR) then
               dencomp=(imcrit-1)*num_state_material+ENUM_DENVAR+1
               fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),dencomp+1)
              else
               print *,"ENUM_TEMPERATUREVAR<>ENUM_DENVAR+1"
               stop
              endif
             else
              print *,"n invalid"
              stop
             endif
            else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
                     (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                     (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then 
             if (n.eq.1) then
              if ((dir.ge.0).and.(dir.lt.AMREX_SPACEDIM)) then
               fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),dir+1)
              else
               print *,"dir invalid"
               stop
              endif
             else
              print *,"n invalid"
              stop
             endif
            else if (operation_flag.eq.OP_PRESGRAD_MAC) then
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
           ncomp_flux,bfact_f,grid_type, &
           flochi,dxf,xcoarse,ffine,crse_value)

          voltotal=one

         else if ((bfact_f.eq.1).or. &
                  (enable_spectral.eq.0).or. &
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

                 if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then ! viscosity
                  if ((n.ge.1).and.(n.le.STATE_NCOMP_VEL)) then
                   fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),n)
                  else
                   print *,"n invalid"
                   stop
                  endif
                 else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection
                  if ((n.ge.SEM_U+1).and.(n.le.SEM_W+1)) then ! velocity
                   fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),n)
                  else if (n.eq.SEM_T+1) then ! temperature
                   if (ENUM_DENVAR+1.eq.ENUM_TEMPERATUREVAR) then
                    dencomp=(imcrit-1)*num_state_material+ENUM_DENVAR+1
                    fine_data=den_fine(D_DECL(istrip,jstrip,kstrip),dencomp+1)
                   else
                    print *,"ENUM_TEMPERATUREVAR<>ENUM_DENVAR+1"
                    stop
                   endif
                  else
                   print *,"n invalid"
                   stop
                  endif
                 else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
                      (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                      (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then 
                  if (n.eq.1) then
                   if ((dir.ge.0).and.(dir.lt.AMREX_SPACEDIM)) then
                    fine_data=vel_fine(D_DECL(istrip,jstrip,kstrip),dir+1)
                   else
                    print *,"dir invalid"
                    stop
                   endif
                  else
                   print *,"n invalid"
                   stop
                  endif
                 else if (operation_flag.eq.OP_PRESGRAD_MAC) then
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
          print *,"voltotal invalid: ",voltotal
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
      end subroutine fort_avgdown_copy


      subroutine fort_interp_copy( &
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
       loc,hic) &
      bind(c,name='fort_interp_copy')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: enable_spectral
      real(amrex_real), INTENT(in) :: dxc(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: operation_flag
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: problo(SDIM)
      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      integer, INTENT(in) :: ncomp_flux
      integer, INTENT(in) :: ncomp_den
      integer, INTENT(in) :: ncomp_vel
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: DIMDEC(den_crse)
      integer, INTENT(in) :: DIMDEC(vel_crse)
      integer, INTENT(in) :: DIMDEC(masknbr)
      integer, INTENT(in) :: DIMDEC(masksem)
      integer, INTENT(in) :: DIMDEC(cmasksem)
      integer, INTENT(in) :: DIMDEC(coarseLS)
      integer, INTENT(in) :: velbc(SDIM,2,SDIM)
      integer, INTENT(in) :: loc(SDIM),hic(SDIM) ! coarse grid dimensions
      integer :: growlo(3),growhi(3)
      integer :: stenlo(3),stenhi(3)
      integer :: mstenlo(3),mstenhi(3)
      real(amrex_real), INTENT(out),target :: fine(DIMV(fine),ncomp_flux)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: den_crse(DIMV(den_crse),ncomp_den)
      real(amrex_real), pointer :: den_crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: vel_crse(DIMV(vel_crse),ncomp_vel)
      real(amrex_real), pointer :: vel_crse_ptr(D_DECL(:,:,:),:)
       ! =1 if fine-fine  =0 coarse-fine
      real(amrex_real), INTENT(in),target :: masknbr(DIMV(masknbr))  
      real(amrex_real), pointer :: masknbr_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: masksem(DIMV(masksem))
      real(amrex_real), pointer :: masksem_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: cmasksem(DIMV(cmasksem))
      real(amrex_real), pointer :: cmasksem_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: coarseLS(DIMV(coarseLS),num_materials)
      real(amrex_real), pointer :: coarseLS_ptr(D_DECL(:,:,:),:)

      integer clochi(SDIM)
      integer ifine,jfine,kfine
      integer icoarse,jcoarse,kcoarse
      integer ilocal,jlocal,klocal
      integer istrip,jstrip,kstrip
      integer iface,jface,kface
      integer ii,jj,kk
      integer n
      integer dir2
      integer side
      integer grid_type
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_coarse(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) fine_value(ncomp_flux)
      real(amrex_real) dxelem_c
      real(amrex_real) xfine(SDIM)
      integer coarselo_index
      real(amrex_real), dimension(D_DECL(:,:,:),:),allocatable :: ccrse
      integer testmask,testmask2
      integer elem_test
      integer im,imcrit
      integer mat_freq(num_materials)
      real(amrex_real) coarseLS_local(num_materials)
      real(amrex_real) crse_data
      integer dencomp
      integer local_masknbr


      if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then ! viscosity
       if (ncomp_vel.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_vel invalid viscosity"
        stop
       endif
       if (ncomp_den.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_den invalid viscosity"
        stop
       endif
       if (ncomp_flux.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_flux invalid viscosity"
        stop
       endif
      else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection
       if (ncomp_vel.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.num_materials*num_state_material) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.NFLUXSEM) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
               (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
               (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then 
       if (ncomp_vel.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_vel invalid"
        stop
       endif
       if (ncomp_den.ne.STATE_NCOMP_VEL) then
        print *,"ncomp_den invalid"
        stop
       endif
       if (ncomp_flux.ne.1) then
        print *,"ncomp_flux invalid"
        stop
       endif
      else if (operation_flag.eq.OP_PRESGRAD_MAC) then
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

      fine_ptr=>fine
      call checkbound_array(fablo,fabhi,fine_ptr,0,dir)
      den_crse_ptr=>den_crse
      call checkbound_array(loc,hic,den_crse_ptr,1,-1)
      vel_crse_ptr=>vel_crse
      call checkbound_array(loc,hic,vel_crse_ptr,1,-1)
      masknbr_ptr=>masknbr
      call checkbound_array1(fablo,fabhi,masknbr_ptr,1,-1)
      masksem_ptr=>masksem
      call checkbound_array1(fablo,fabhi,masksem_ptr,1,-1)
      cmasksem_ptr=>cmasksem
      call checkbound_array1(loc,hic,cmasksem_ptr,1,-1)
      coarseLS_ptr=>coarseLS
      call checkbound_array(loc,hic,coarseLS_ptr,1,-1)

      grid_type=-1  ! ggg  (Gauss in all directions)

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
 
        do kfine=growlo(3),growhi(3)
        do jfine=growlo(2),growhi(2)
        do ifine=growlo(1),growhi(1)

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
           coarselo_index=NINT((xcomp_of_xphys(dir2-1,xfine(dir2))- &
                                problo(dir2))/dxelem_c-half )
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

           if ((xfine(dir2).ge.-EPS3*dxc(dir2)).and. &
               (xfine(dir2).le.(EPS3+bfact_c)*dxc(dir2))) then
            !do nothing
           else
            print *,"xfine out of bounds interp_copy: ",dir2,xfine(dir2)
            stop
           endif

          enddo ! dir2=1..sdim

          do im=1,num_materials
           mat_freq(im)=0
          enddo
           ! mstenlo,mstenhi is a thin coarse level stencil
          do klocal=mstenlo(3),mstenhi(3) 
          do jlocal=mstenlo(2),mstenhi(2) 
          do ilocal=mstenlo(1),mstenhi(1) 
           do im=1,num_materials
            coarseLS_local(im)=coarseLS(D_DECL(ilocal,jlocal,klocal),im)
           enddo
           call get_primary_material(coarseLS_local,imcrit)
           if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then
            mat_freq(imcrit)=mat_freq(imcrit)+1
           else
            print *,"imcrit invalid"
            stop
           endif
          enddo  !klocal
          enddo  !jlocal
          enddo  !ilocal

          imcrit=0
          do im=1,num_materials
           if (mat_freq(im).gt.0) then
            if (imcrit.eq.0) then
             imcrit=im
            else if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then
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
          enddo ! im=1..num_materials

          ilocal=mstenlo(1) 
          jlocal=mstenlo(2) 
          klocal=mstenlo(SDIM) 
          testmask=NINT(cmasksem(D_DECL(ilocal,jlocal,klocal)))
          do klocal=mstenlo(3),mstenhi(3) 
          do jlocal=mstenlo(2),mstenhi(2) 
          do ilocal=mstenlo(1),mstenhi(1) 
           testmask2=NINT(cmasksem(D_DECL(ilocal,jlocal,klocal)))
           if (testmask2.ne.testmask) then
            testmask=0
           endif
          enddo 
          enddo 
          enddo 

          if ((imcrit.ge.1).and.(imcrit.le.num_materials)) then

           do n=1,ncomp_flux
            fine_value(n)=zero
           enddo
           voltotal=zero

           if ((bfact_c.ge.2).and. &
               (enable_spectral.eq.1).and. &
               (testmask.ge.1).and. &
               (testmask.le.num_materials)) then

            do kcoarse=stenlo(3),stenhi(3)
            do jcoarse=stenlo(2),stenhi(2)
            do icoarse=stenlo(1),stenhi(1)
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

              if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then ! viscosity
               if ((n.ge.1).and.(n.le.STATE_NCOMP_VEL)) then
                crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),n)
               else
                print *,"n invalid viscosity"
                stop
               endif
              else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection
               if ((n.ge.SEM_U+1).and.(n.le.SEM_W+1)) then ! velocity
                crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),n)
               else if (n.eq.SEM_T+1) then ! temperature
                if (ENUM_DENVAR+1.eq.ENUM_TEMPERATUREVAR) then
                 dencomp=(imcrit-1)*num_state_material+ENUM_DENVAR+1
                 crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),dencomp+1)
                else
                 print *,"ENUM_TEMPERATUREVAR<>ENUM_DENVAR+1"
                 stop
                endif
               else
                print *,"n invalid"
                stop
               endif
              else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. & 
                       (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. & 
                       (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then 
               if (n.eq.1) then
                if ((dir.ge.0).and.(dir.lt.AMREX_SPACEDIM)) then
                 crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),dir+1)
                else
                 print *,"dir invalid"
                 stop
                endif
               else
                print *,"n invalid"
                stop
               endif
              else if (operation_flag.eq.OP_PRESGRAD_MAC) then
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
             ncomp_flux,bfact_c,grid_type, &
             clochi,dxc,xfine,ccrse,fine_value)

            voltotal=one

           else if ((bfact_c.eq.1).or. &
                    (enable_spectral.eq.0).or. &
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

                   if (operation_flag.eq.OP_UGRAD_COUPLING_MAC) then !viscosity
                    if ((n.ge.1).and.(n.le.STATE_NCOMP_VEL)) then
                     crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),n)
                    else
                     print *,"n invalid"
                     stop
                    endif
                   else if (operation_flag.eq.OP_ISCHEME_MAC) then ! advection
                    if ((n.ge.SEM_U+1).and.(n.le.SEM_W+1)) then ! velocity
                     crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),n)
                    else if (n.eq.SEM_T+1) then ! temperature
                     if (ENUM_DENVAR+1.eq.ENUM_TEMPERATUREVAR) then
                      dencomp=(imcrit-1)*num_state_material+ENUM_DENVAR+1
                      crse_data=den_crse(D_DECL(istrip,jstrip,kstrip),dencomp+1)
                     else
                      print *,"ENUM_TEMPERATUREVAR<>ENUM_DENVAR+1"
                      stop
                     endif
                    else
                     print *,"n invalid"
                     stop
                    endif
                   else if ((operation_flag.eq.OP_UNEW_CELL_TO_MAC).or. &
                      (operation_flag.eq.OP_UMAC_PLUS_VISC_CELL_TO_MAC).or. &
                      (operation_flag.eq.OP_U_SEM_CELL_MAC_TO_MAC)) then 
                    if (n.eq.1) then
                     if ((dir.ge.0).and.(dir.lt.AMREX_SPACEDIM)) then
                      crse_data=vel_crse(D_DECL(istrip,jstrip,kstrip),dir+1)
                     else
                      print *,"dir invalid"
                      stop
                     endif
                    else
                     print *,"n invalid"
                     stop
                    endif
                   else if (operation_flag.eq.OP_PRESGRAD_MAC) then
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
            print *,"voltotal invalid: ",voltotal
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
      end subroutine fort_interp_copy

      subroutine fort_interp_flux( &
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
       loc,hic) &
      bind(c,name='fort_interp_flux')

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: enable_spectral
      real(amrex_real), INTENT(in) :: dxc(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: problo(SDIM)
      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      integer, INTENT(in) :: ncomp_flux
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(masknbr)
      integer, INTENT(in) :: DIMDEC(masksem)
      integer, INTENT(in) :: DIMDEC(cmasksem)
      integer, INTENT(in) :: velbc(SDIM,2,SDIM)
      integer, INTENT(in) :: loc(SDIM),hic(SDIM) ! coarse grid dimensions
      integer :: growlo(3),growhi(3)
      integer :: stenlo(3),stenhi(3)
      integer :: mstenlo(3),mstenhi(3)
      real(amrex_real), INTENT(inout),target ::   fine(DIMV(fine),ncomp_flux)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target ::   crse(DIMV(crse),ncomp_flux)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
       ! =1 if fine-fine  =0 coarse-fine
      real(amrex_real), INTENT(in),target ::   masknbr(DIMV(masknbr)) 
      real(amrex_real), pointer :: masknbr_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target ::   masksem(DIMV(masksem))
      real(amrex_real), pointer :: masksem_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target ::   cmasksem(DIMV(cmasksem))
      real(amrex_real), pointer :: cmasksem_ptr(D_DECL(:,:,:))

      integer clochi(SDIM)
      integer ifine,jfine,kfine
      integer icoarse,jcoarse,kcoarse
      integer ilocal,jlocal,klocal
      integer istrip,jstrip,kstrip
      integer iface,jface,kface
      integer ii,jj,kk
      integer n
      integer dir2
      integer side
      integer grid_type
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_coarse(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) fine_value(ncomp_flux)
      real(amrex_real) dxelem_c
      real(amrex_real) xfine(SDIM)
      integer coarselo_index
      real(amrex_real), dimension(D_DECL(:,:,:),:),allocatable :: ccrse
      integer testmask,testmask2
      integer elem_test
      real(amrex_real) crse_data
      integer local_masknbr


      if (ncomp_flux.ne.NFLUXSEM) then
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

      fine_ptr=>fine
      call checkbound_array(fablo,fabhi,fine,0,dir)
      crse_ptr=>crse
      call checkbound_array(loc,hic,crse_ptr,0,dir)
      masknbr_ptr=>masknbr
      call checkbound_array1(fablo,fabhi,masknbr_ptr,1,-1)
      masksem_ptr=>masksem
      call checkbound_array1(fablo,fabhi,masksem_ptr,1,-1)
      cmasksem_ptr=>cmasksem
      call checkbound_array1(loc,hic,cmasksem_ptr,1,-1)

      grid_type=-1  ! ggg  (Gauss in all directions)

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
 
        do kfine=growlo(3),growhi(3)
        do jfine=growlo(2),growhi(2)
        do ifine=growlo(1),growhi(1)

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
           coarselo_index=NINT( (xcomp_of_xphys(dir2-1,xfine(dir2))- &
                                 problo(dir2))/dxelem_c-half )
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

           if ((xfine(dir2).ge.-EPS3*dxc(dir2)).and. &
               (xfine(dir2).le.(EPS3+bfact_c)*dxc(dir2))) then
            !do nothing
           else
            print *,"xfine out of bounds interp_copy: ",dir2,xfine(dir2)
            stop
           endif

          enddo ! dir2=1..sdim

          ilocal=mstenlo(1) 
          jlocal=mstenlo(2) 
          klocal=mstenlo(SDIM) 
          testmask=NINT(cmasksem(D_DECL(ilocal,jlocal,klocal)))
          do klocal=mstenlo(3),mstenhi(3) 
          do jlocal=mstenlo(2),mstenhi(2) 
          do ilocal=mstenlo(1),mstenhi(1) 
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
              (enable_spectral.eq.1).and. &
              (testmask.ge.1).and. &
              (testmask.le.num_materials)) then

            do kcoarse=stenlo(3),stenhi(3)
            do jcoarse=stenlo(2),stenhi(2)
            do icoarse=stenlo(1),stenhi(1)
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
             ncomp_flux,bfact_c,grid_type, &
             clochi,dxc,xfine,ccrse,fine_value)

            voltotal=one

          else if ((bfact_c.eq.1).or. &
                   (enable_spectral.eq.0).or. &
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

                  do n = 1,ncomp_flux
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
           print *,"voltotal invalid: ",voltotal
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
      end subroutine fort_interp_flux


      subroutine fort_fillbdry_flux( &
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
       presbc) &
       bind(c,name='fort_fillbdry_flux')

      use global_utility_module

      IMPLICIT NONE

      integer, INTENT(in) :: sync_iter
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer, INTENT(in) :: dir
      integer, INTENT(in) :: ncomp_flux
      integer, INTENT(in) :: DIMDEC(fluxtarg)
      integer, INTENT(in) :: DIMDEC(fluxhold)
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(masknbr)
      integer, INTENT(in) :: presbc(SDIM,2)
       ! maskcov=tag if not covered by level+1 or outside the domain.
      real(amrex_real), INTENT(in),target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
       ! =1 if fine-fine  =0 coarse-fine
      real(amrex_real), INTENT(in),target :: masknbr(DIMV(masknbr))  
      real(amrex_real), pointer :: masknbr_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout),target :: fluxtarg(DIMV(fluxtarg),ncomp_flux) 
      real(amrex_real), pointer :: fluxtarg_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: fluxhold(DIMV(fluxhold),ncomp_flux) 
      real(amrex_real), pointer :: fluxhold_ptr(D_DECL(:,:,:),:)
      integer i,j,k
      integer iface,jface,kface
      integer i_out,j_out,k_out
      integer ii,jj,kk
      integer dirloc,side
      integer gridlo(3),gridhi(3)
      integer nc
      integer pbc
      integer maskcov_in
      integer maskcov_out
      integer masknbr_out
      real(amrex_real) flux_orig,flux_copy

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

      fluxtarg_ptr=>fluxtarg
      call checkbound_array(fablo,fabhi,fluxtarg_ptr,0,dir)
      fluxhold_ptr=>fluxhold
      call checkbound_array(fablo,fabhi,fluxhold_ptr,1,-1)
      masknbr_ptr=>masknbr
      call checkbound_array1(fablo,fabhi,masknbr_ptr,1,-1)
      maskcov_ptr=>maskcov
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)

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

       do k=gridlo(3),gridhi(3)
       do j=gridlo(2),gridhi(2)
       do i=gridlo(1),gridhi(1)

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
      end subroutine fort_fillbdry_flux

      subroutine fort_summass( &
        tid, &
        level, &
        finest_level, &
        ncomp_sum_int_user1, &
        ncomp_sum_int_user2, &
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
        visco,DIMS(visco), &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        time, &
        local_result, & !running total
        resultALL, &    !data from previous sweep
        sumdata_type, &
        sumdata_sweep, &
        resultsize, &
        coflow_num_cells, &
        coflow_Z, &
        coflow_R_of_Z, &
        Z_dir, &
        R_dir, &
        den_ncomp, &
        isweep) &  ! isweep=0 or 1
      bind(c,name='fort_summass')

      use LegendreNodes
      use global_utility_module
      use probf90_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: ncomp_sum_int_user1
      integer, INTENT(in) :: ncomp_sum_int_user2
      integer :: ncomp_sum_int_user12
      integer, INTENT(in) :: tid
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: adapt_quad_depth
      integer :: max_level_adapt
      integer, INTENT(in) :: slice_dir
      real(amrex_real), INTENT(in) :: xslice(SDIM)
      integer, INTENT(in) :: resultsize
      integer, INTENT(in) :: den_ncomp
      integer, INTENT(in) :: isweep  ! isweep=0 or 1
      integer, INTENT(in) :: coflow_num_cells,Z_dir,R_dir
      real(amrex_real), INTENT(out) :: coflow_Z(0:coflow_num_cells)
      real(amrex_real), INTENT(out) :: coflow_R_of_Z(0:coflow_num_cells)

      real(amrex_real), INTENT(in), target :: problo(SDIM)
      real(amrex_real), INTENT(in), target :: probhi(SDIM)
      integer, INTENT(in) :: bfact
      integer, INTENT(in), target :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in), target :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: DIMDEC(cellten)
      integer, INTENT(in) :: DIMDEC(lsfab)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: DIMDEC(maskSEM)
      integer, INTENT(in) :: DIMDEC(drag)
      integer, INTENT(in) :: DIMDEC(slopes)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(vel)
      integer, INTENT(in) :: DIMDEC(visco)
      real(amrex_real), INTENT(inout) ::  local_result(resultsize)
      real(amrex_real), INTENT(in) ::  resultALL(resultsize)
      integer, INTENT(in) :: sumdata_type(resultsize)
      integer, INTENT(in) :: sumdata_sweep(resultsize)
      real(amrex_real), INTENT(in) ::  time
      real(amrex_real), INTENT(in), target :: &
              cellten(DIMV(cellten),AMREX_SPACEDIM_SQR)  
      real(amrex_real), pointer :: cellten_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: lsfab(DIMV(lsfab),num_materials)  
      real(amrex_real), pointer :: lsfab_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target ::  maskSEM(DIMV(maskSEM))
      real(amrex_real), pointer :: maskSEM_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target ::  mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target ::  drag(DIMV(drag),N_DRAG)
      real(amrex_real), pointer :: drag_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              slopes(DIMV(slopes),num_materials*ngeom_recon)  
      real(amrex_real), pointer :: slopes_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: den(DIMV(den),den_ncomp)  
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:),:)
      ! includes pressure 
      real(amrex_real), INTENT(in), target ::  &
          vel(DIMV(vel),STATE_NCOMP_VEL+STATE_NCOMP_PRES) 
      real(amrex_real), pointer :: vel_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: visco(DIMV(visco), &
         ENUM_NUM_TENSOR_TYPE*num_materials_viscoelastic) 
      real(amrex_real), pointer :: visco_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: xlo(SDIM),dx(SDIM)

      integer, parameter :: nhalf=3

      integer i,j,k
      integer ii,jj,kk
      integer dir
      integer im
      integer im_primary
      integer vofcomp
      integer dirMAC
      real(amrex_real) xcen,xcrit
      real(amrex_real), target :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenMAC(-nhalf:nhalf,SDIM)
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdata_tess(num_materials*ngeom_recon)
      integer stack_error_level
      integer dir2,dir3
      real(amrex_real) errorparm(num_materials*2) ! fi,ei
      real(amrex_real) xbottom,xtop,ls_above,ls_below,ZZgrid
      real(amrex_real) vof_below,vof_above,vof_face
      real(amrex_real) volgrid,distbound
      real(amrex_real) cengrid(SDIM)
      real(amrex_real) xboundary(SDIM)

      integer isrc,idest
      integer iside
      real(amrex_real) dz_external
      integer j_external

      real(amrex_real) cen_material(SDIM)
      integer use_vof_height

      real(amrex_real) ldata(D_DECL(3,3,3))
      real(amrex_real) LSvolume,LSfacearea
      real(amrex_real) local_flotsam
      real(amrex_real) LScentroid(SDIM)
      real(amrex_real) LScen_material(SDIM)
      integer i1,j1,k1,k1lo,k1hi
      real(amrex_real) KECELL,dencore,Tcore,ecore,totalE
      integer in_slice
      integer nbase,veldir
      real(amrex_real) gradu(3,3)
      real(amrex_real) vort(3)

      real(amrex_real) local_vort
      real(amrex_real) local_temperature
      real(amrex_real) local_vort_error
      real(amrex_real) local_temperature_error
      real(amrex_real) vort_expect
      real(amrex_real) temperature_expect
      real(amrex_real) local_vel_error
      real(amrex_real) local_energy_moment
      real(amrex_real) local_vel(SDIM)
      real(amrex_real) vel_expect(SDIM)
      real(amrex_real) local_xsten(SDIM)

      real(amrex_real) local_kinetic_energy(num_materials)

      real(amrex_real) rr
      integer mask_im
      integer cell_index,element_index,sub_index
      integer velcomp
      integer nmax

      real(amrex_real) LS_LOCAL(num_materials)

      real(amrex_real) massfrac_parm(num_species_var+1)
      integer ispec
      integer local_tessellate

      type(user_defined_sum_int_type) :: GRID_DATA_PARM
      real(amrex_real) local_user_out1(ncomp_sum_int_user1+1)
      real(amrex_real) local_user_out2(ncomp_sum_int_user2+1)

      ncomp_sum_int_user12=ncomp_sum_int_user1+ncomp_sum_int_user2

      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid fort_summass"
       stop
      endif

      nmax=POLYGON_LIST_MAX  ! in: fort_summass

      cellten_ptr=>cellten
      lsfab_ptr=>lsfab
      maskSEM_ptr=>maskSEM
      mask_ptr=>mask
      drag_ptr=>drag
      slopes_ptr=>slopes
      den_ptr=>den
      vel_ptr=>vel
      visco_ptr=>visco

      if ((adapt_quad_depth.lt.1).or.(adapt_quad_depth.gt.10)) then
       print *,"adapt_quad_depth invalid"
       stop
      endif
      if ((slice_dir.lt.0).or.(slice_dir.ge.SDIM)) then
       print *,"slice_dir invalid"
       stop
      endif
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid sum mass ",ngeom_recon
       stop
      endif
      if (den_ncomp.ne.num_materials*num_state_material) then
       print *,"den_ncomp invalid"
       stop
      endif
      if ((Z_dir.lt.0).or.(Z_dir.ge.SDIM)) then
       print *,"Z_dir invalid"
       stop
      endif
      if ((R_dir.lt.0).or.(R_dir.ge.SDIM).or.(R_dir.eq.Z_dir)) then
       print *,"R_dir invalid"
       stop
      endif
      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif

      if (resultsize.ne.IQ_TOTAL_SUM_COMP) then
       print *,"mismatch between resultsize and IQ_TOTAL_SUM_COMP"
       stop
      endif
      do idest=1,IQ_TOTAL_SUM_COMP
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
       max_level_adapt=adapt_quad_depth
      else if (isweep.eq.1) then
       max_level_adapt=1
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

      call checkbound_array(fablo,fabhi,cellten_ptr,0,-1) 
      call checkbound_array(fablo,fabhi,lsfab_ptr,2,-1) 
      call checkbound_array1(fablo,fabhi,maskSEM_ptr,1,-1) 
      call checkbound_array1(fablo,fabhi,mask_ptr,2,-1) 
       ! ngrow_distance=4
       ! ngrow_make_distance=3
      call checkbound_array(fablo,fabhi,drag_ptr,3,-1) 
      call checkbound_array(fablo,fabhi,slopes_ptr,2,-1) 
      call checkbound_array(fablo,fabhi,den_ptr,1,-1) 
      call checkbound_array(fablo,fabhi,vel_ptr,1,-1) 
      call checkbound_array(fablo,fabhi,visco_ptr,1,-1) 

      GRID_DATA_PARM%ncomp_sum_int_user1=ncomp_sum_int_user1
      GRID_DATA_PARM%ncomp_sum_int_user2=ncomp_sum_int_user2
      GRID_DATA_PARM%ncomp_sum_int_user12=ncomp_sum_int_user12
      GRID_DATA_PARM%time=time
      GRID_DATA_PARM%problo=>problo
      GRID_DATA_PARM%probhi=>probhi
      GRID_DATA_PARM%nhalf=nhalf
      GRID_DATA_PARM%bfact=bfact
      GRID_DATA_PARM%den_ncomp=den_ncomp
      GRID_DATA_PARM%level=level
      GRID_DATA_PARM%finest_level=finest_level
      GRID_DATA_PARM%tilelo=>tilelo
      GRID_DATA_PARM%tilehi=>tilehi
      GRID_DATA_PARM%fablo=>fablo
      GRID_DATA_PARM%fabhi=>fabhi
      GRID_DATA_PARM%xlo=>xlo
      GRID_DATA_PARM%dx=>dx
      GRID_DATA_PARM%cellten=>cellten
      GRID_DATA_PARM%lsfab=>lsfab
      GRID_DATA_PARM%slopes=>slopes
      GRID_DATA_PARM%den=>den
      GRID_DATA_PARM%vel=>vel
      GRID_DATA_PARM%visco=>visco

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       if (mask(D_DECL(i,j,k)).gt.zero) then

        call gridsten(xsten,xlo,i,j,k,fablo,bfact,dx,nhalf)

        do veldir=1,3
        do dir=1,3
         gradu(veldir,dir)=zero
        enddo
        enddo

        do dir=1,SDIM
         if (dir.eq.1) then
          nbase=TENSOR_TRANSPOSE_UX-1
         else if (dir.eq.2) then
          nbase=TENSOR_TRANSPOSE_UY-1
         else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
          nbase=TENSOR_TRANSPOSE_UZ-1
         else
          print *,"dir invalid summass"
          stop
         endif

         do veldir=1,SDIM
          gradu(veldir,dir)=cellten(D_DECL(i,j,k),nbase+veldir) 
         enddo
        enddo ! dir

        velcomp=1

        if (levelrz.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
         rr=xsten(0,1)
         gradu(3,3)=vel(D_DECL(i,j,k),velcomp)/abs(rr)
        else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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

        do dir=1,num_materials*ngeom_recon
         mofdata(dir)=slopes(D_DECL(i,j,k),dir)
         mofdata_tess(dir)=mofdata(dir)
        enddo

        stack_error_level=0

        ! before (mofdata): fluids tessellate
        ! after  (mofdata): fluids and solids tessellate
        ! EPS2
        local_tessellate=1
        call multi_get_volume_tessellate( &
         local_tessellate, &
         bfact, &
         dx,xsten,nhalf, &
         mofdata_tess, &
         geom_xtetlist(1,1,1,tid+1), &
         nmax, &
         nmax, &
         SDIM)

         ! tessellate==1 (internal to stackerror)
         ! in: fort_summass
        call stackerror( &
         geom_xtetlist(1,1,1,tid+1), &
         xsten,nhalf,dx,bfact, &
         xsten,nhalf, &
         mofdata, &
         mofdata_tess, &
         errorparm, &
         stack_error_level, &
         max_level_adapt, &
         time)

         ! F1,E1,F2,E2,F3,E3,...
        do dir=1,2*num_materials
         idest=IQ_FE_SUM_COMP+dir
         local_result(idest)=local_result(idest)+errorparm(dir)
        enddo

        call Box_volumeFAST(bfact,dx,xsten,nhalf, &
         volgrid,cengrid,SDIM)

        mask_im=NINT(maskSEM(D_DECL(i,j,k)))
        if ((mask_im.ge.1).and.(mask_im.le.num_materials)) then
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

         if (levelrz.eq.COORDSYS_CARTESIAN) then
          ! do nothing
         else if (levelrz.eq.COORDSYS_RZ) then
          if (SDIM.ne.2) then
           print *,"dimension bust"
           stop
          endif
          volgrid=volgrid*two*Pi*xsten(0,1)
         else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
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
      
        GRID_DATA_PARM%volgrid=volgrid
        GRID_DATA_PARM%igrid=i
        GRID_DATA_PARM%jgrid=j
        GRID_DATA_PARM%kgrid=k
        GRID_DATA_PARM%xsten=>xsten

        do im=1,num_materials
         LS_LOCAL(im)=lsfab(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(LS_LOCAL,im_primary)

        do im=1,num_materials

         vofcomp=(im-1)*ngeom_recon+1

         local_flotsam=mofdata_tess(vofcomp)

         do k1=k1lo,k1hi
         do j1=-1,1
         do i1=-1,1
          ldata(D_DECL(i1+2,j1+2,k1+2))= &
           lsfab(D_DECL(i+i1,j+j1,k+k1),im)

          if (ldata(D_DECL(i1+2,j1+2,k1+2)).ge.zero) then
           local_flotsam=zero
          endif

         enddo
         enddo
         enddo
          ! LSvolume is a volume fraction.
          ! LScentroid is in an absolute coordinate system.
         call getvolume(bfact,dx,xsten,nhalf, &
          ldata,LSvolume,LSfacearea, &
          LScentroid,VOFTOL,SDIM)

         do dir=1,SDIM
          cen_material(dir)=mofdata_tess(vofcomp+dir)+cengrid(dir)
          LScen_material(dir)=LScentroid(dir)
         enddo

         do dir=1,SDIM
          idest=IQ_CEN_SUM_COMP+dir+3*(im-1)
          local_result(idest)=local_result(idest)+ &
            volgrid*mofdata_tess(vofcomp)*cen_material(dir)
          idest=IQ_LS_CEN_SUM_COMP+dir+3*(im-1)
          local_result(idest)=local_result(idest)+ &
            volgrid*LSvolume*LScen_material(dir)
         enddo !dir=1,SDIM

         idest=IQ_F_FLOTSAM_COMP+im
         local_result(idest)=local_result(idest)+volgrid*local_flotsam

         idest=IQ_LS_F_SUM_COMP+im
         local_result(idest)=local_result(idest)+volgrid*LSvolume

         dencore=den(D_DECL(i,j,k),ENUM_DENVAR+1+num_state_material*(im-1))
         Tcore=den(D_DECL(i,j,k), &
                   ENUM_TEMPERATUREVAR+1+num_state_material*(im-1))

         idest=IQ_MASS_SUM_COMP+im
         local_result(idest)=local_result(idest)+ &
           dencore*volgrid*mofdata_tess(vofcomp)

         do ispec=1,num_species_var
          idest=IQ_SPECIES_MASS_SUM_COMP+(ispec-1)*num_materials+im
          local_result(idest)=local_result(idest)+ &
            dencore*volgrid*mofdata_tess(vofcomp)* &
            den(D_DECL(i,j,k),2+ispec+num_state_material*(im-1))
         enddo ! ispec=1..num_species_var

         KECELL=zero

         do dir=1,SDIM
          velcomp=dir
          idest=IQ_MOM_SUM_COMP+3*(im-1)+dir
          local_result(idest)=local_result(idest)+ &
            vel(D_DECL(i,j,k),velcomp)*volgrid*mofdata_tess(vofcomp)*dencore

          KECELL=KECELL+vel(D_DECL(i,j,k),velcomp)**2
         enddo ! dir=1..sdim

         KECELL=half*KECELL

         idest=IQ_KINETIC_ENERGY_SUM_COMP+im

         local_kinetic_energy(im)= &
           KECELL*volgrid*mofdata_tess(vofcomp)*dencore

         local_result(idest)=local_result(idest)+local_kinetic_energy(im)

         idest=IQ_ENERGY_SUM_COMP+im

         call init_massfrac_parm(dencore,massfrac_parm,im)
         do ispec=1,num_species_var
          massfrac_parm(ispec)= &
              den(D_DECL(i,j,k),(im-1)*num_state_material+2+ispec)
         enddo
         call INTERNAL_material(dencore,massfrac_parm, &
           Tcore,ecore, &
           fort_material_type(im),im)

         totalE=dencore*(KECELL+ecore)*volgrid*mofdata_tess(vofcomp)
         local_result(idest)=local_result(idest)+totalE

        enddo  ! im=1..num_materials

        do dir=1,3
         idest=IQ_VORT_SUM_COMP+dir
         local_result(idest)=local_result(idest)+volgrid*vort(dir)
         idest=IQ_ENSTROPHY_SUM_COMP+1
         if ((im_primary.ge.1).and.(im_primary.le.num_materials)) then
          local_result(idest+im_primary-1)= &
              local_result(idest+im_primary-1)+volgrid*(vort(dir)**2)
         else
          print *,"im_primary invalid"
          stop
         endif
        enddo ! dir=1..3
      
        if ((ncomp_sum_int_user1.ge.1).or. &
            (ncomp_sum_int_user2.ge.1)) then

         if (is_in_probtype_list().eq.1) then

          do im=1,ncomp_sum_int_user1
           idest=IQ_USER_SUM_COMP+im
           if (sumdata_sweep(idest).eq.0) then
            ! do nothing
           else
            print *,"sumdata_sweep invalid"
            stop
           endif
           if (isweep.eq.0) then
            local_user_out1(im)=zero
           else if (isweep.eq.1) then
            local_user_out1(im)=resultALL(idest)
           else
            print *,"isweep invalid"
            stop
           endif
          enddo !im=1,ncomp_sum_int_user1

          do im=1,ncomp_sum_int_user2
           idest=IQ_USER_SUM_COMP+ncomp_sum_int_user1+im
           if (sumdata_sweep(idest).eq.1) then
            ! do nothing
           else
            print *,"sumdata_sweep invalid"
            stop
           endif
           local_user_out2(im)=zero
          enddo !im=1,ncomp_sum_int_user2

          call SUB_SUMINT(GRID_DATA_PARM,local_user_out1, &
           local_user_out2,ncomp_sum_int_user1, &
           ncomp_sum_int_user2,isweep)

          if (isweep.eq.0) then
           do im=1,ncomp_sum_int_user1
            idest=IQ_USER_SUM_COMP+im
            local_result(idest)=local_result(idest)+local_user_out1(im)
           enddo
          else if (isweep.eq.1) then
           do im=1,ncomp_sum_int_user2
            idest=IQ_USER_SUM_COMP+ncomp_sum_int_user1+im
            local_result(idest)=local_result(idest)+local_user_out2(im)
           enddo
          else
           print *,"isweep invalid"
           stop
          endif

         endif ! if (is_in_probtype_list().eq.1) then

        else if ((ncomp_sum_int_user1.eq.0).and. &
                 (ncomp_sum_int_user2.eq.0)) then
         ! do nothing
        else
         print *,"ncomp_sum_int_user1 or ncomp_sum_int_user2 invalid"
         stop
        endif
  
        local_vort=sqrt(vort(1)**2+vort(2)**2+vort(3)**2)
        do dir=1,SDIM
         local_vel(dir)=vel(D_DECL(i,j,k),dir)
         local_xsten(dir)=xsten(0,dir)
        enddo
        local_temperature=den(D_DECL(i,j,k),ENUM_TEMPERATUREVAR+1)
        call get_vort_vel_error(time,local_xsten, &
         local_vort,local_vel,local_temperature, &
         local_vort_error,local_vel_error,local_temperature_error, &
         vort_expect,vel_expect,temperature_expect, &
         local_energy_moment)
        
        local_result(IQ_ENERGY_MOMENT_SUM_COMP+1)= &
             local_result(IQ_ENERGY_MOMENT_SUM_COMP+1)+ &
         local_energy_moment*local_kinetic_energy(1)

        if (local_result(IQ_VORT_ERROR_SUM_COMP+1).lt.local_vort_error) then
         local_result(IQ_VORT_ERROR_SUM_COMP+1)=local_vort_error
        endif

        if (local_result(IQ_TEMP_ERROR_SUM_COMP+1).lt. &
            local_temperature_error) then
         local_result(IQ_TEMP_ERROR_SUM_COMP+1)=local_temperature_error
        endif

        if (local_result(IQ_VEL_ERROR_SUM_COMP+1).lt.local_vel_error) then
         local_result(IQ_VEL_ERROR_SUM_COMP+1)=local_vel_error
        endif
        if (isweep.eq.0) then
         ! do nothing
        else if (isweep.eq.1) then

         if (local_vort_error.gt.zero) then
          if (local_vort_error.gt. &
              resultALL(IQ_VORT_ERROR_SUM_COMP+1)-VOFTOL) then
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

         if (local_temperature_error.gt.zero) then
          if (local_temperature_error.gt. &
              resultALL(IQ_TEMP_ERROR_SUM_COMP+1)-VOFTOL) then
           print *,"**** POSITION OF MAX TEMP ERROR ****"
           do dir=1,SDIM
            print *,"dir,xpos ",dir,local_xsten(dir)
           enddo
           print *,"time ",time
           print *,"local_temperature_error ",local_temperature_error
           print *,"EXPECTED temperature ",temperature_expect
           print *,"ACTUAL temperature ",local_temperature
          endif
         else if (local_temperature_error.eq.zero) then
          ! do nothing
         else
          print *,"local_temperature_error invalid"
          stop
         endif

         if (local_vel_error.gt.zero) then
          if (local_vel_error.gt.resultALL(IQ_VEL_ERROR_SUM_COMP+1)-VOFTOL) then
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

        idest=IQ_LEFT_PRESSURE_SUM_COMP+1
        if (xsten(-1,1).le.problox+VOFTOL*dx(1)) then
         local_result(idest)=local_result(idest)+ &
          volgrid*vel(D_DECL(i,j,k),STATECOMP_PRES+1) 
         local_result(idest+2)=local_result(idest+2)+volgrid
        endif
        idest=IQ_LEFT_PRESSURE_SUM_COMP+2
        if (xsten(1,1).ge.probhix-VOFTOL*dx(1)) then
         local_result(idest)=local_result(idest)+ &
          volgrid*vel(D_DECL(i,j,k),STATECOMP_PRES+1) 
         local_result(idest+2)=local_result(idest+2)+volgrid
        endif

        do im=1,num_materials

         idest=IQ_MINSTATE_SUM_COMP+2*(im-1)+1
         isrc=num_state_material*(im-1)+1
         if (local_result(idest).gt.den(D_DECL(i,j,k),isrc)) then
          local_result(idest)=den(D_DECL(i,j,k),isrc)
         endif
         idest=idest+1
         isrc=isrc+1

         if (local_result(idest).gt.den(D_DECL(i,j,k),isrc)) then
          local_result(idest)=den(D_DECL(i,j,k),isrc)
         endif

         idest=IQ_MAXSTATE_SUM_COMP+2*(im-1)+1
         isrc=num_state_material*(im-1)+1
         if (local_result(idest).lt.den(D_DECL(i,j,k),isrc)) then
          local_result(idest)=den(D_DECL(i,j,k),isrc)
         endif
         idest=idest+1
         isrc=isrc+1

         if (local_result(idest).lt.den(D_DECL(i,j,k),isrc)) then
          local_result(idest)=den(D_DECL(i,j,k),isrc)
         endif

        enddo  ! im=1..num_materials

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

          do im=1,num_materials

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
               (abs(ls_below)+abs(ls_above).ge.EPS2*dx(1))) then
 
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

            idest=IQ_MININT_SUM_COMP+3*(im-1)+dir
            if (xcrit.lt.local_result(idest)) then
             local_result(idest)=xcrit
            endif
            idest=IQ_MAXINT_SUM_COMP+3*(im-1)+dir
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
              idest=IQ_MININT_SLICE_SUM_COMP+im
              if (xcrit.lt.local_result(idest)) then
               local_result(idest)=xcrit
              endif
              idest=IQ_MAXINT_SLICE_SUM_COMP+im
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
             idest=IQ_XNOT_AMP_SUM_COMP+1
             if (xcrit.gt.local_result(idest)) then
              local_result(idest)=xcrit
             endif
            endif

            if (isweep.eq.0) then
             ! do nothing
            else if (isweep.eq.1) then

             do dir2=1,SDIM
              xboundary(dir2)=xsten(0,dir2)
              cengrid(dir2)=resultALL(IQ_CEN_SUM_COMP+dir2+3*(im-1))
             enddo
             xboundary(dir)=xcrit
             distbound=zero
             do dir2=1,SDIM
              distbound=distbound+(xboundary(dir2)-cengrid(dir2))**2
             enddo
             distbound=sqrt(distbound)
             idest=IQ_MINCEN_SUM_COMP+im
             if (distbound.lt.local_result(idest)) then
              local_result(idest)=distbound
             endif
             idest=IQ_MAXCEN_SUM_COMP+im
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

       ! coflow_num_cells is a parameter.
      if (coflow_num_cells.gt.0) then

       if (SDIM.eq.2) then

         ! r=h(z)
        if ((Z_dir.eq.SDIM-1).and.(R_dir.eq.0)) then
       
         k=0 
         do j = growlo(2),growhi(2)+1
         do i = growlo(1),growhi(1)
          if ((mask(D_DECL(i,j,k)).gt.zero).and. &
              (mask(D_DECL(i,j-1,k)).gt.zero).and. &
              (mask(D_DECL(i+1,j,k)).gt.zero).and. &
              (mask(D_DECL(i+1,j-1,k)).gt.zero).and. &
              (mask(D_DECL(i-1,j,k)).gt.zero).and. &
              (mask(D_DECL(i-1,j-1,k)).gt.zero)) then
           dir=R_dir+1  !R_dir=0 r=h(z)  dir=1
           im=1

           do iside=0,1
            dirMAC=1 ! y direction
            call gridstenMAC(xstenMAC,xlo,i,j,k,fablo,bfact,dx,nhalf,dirMAC)

            xcen=xstenMAC(0,dir) ! R_dir=0  dir=R_dir+1=1

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
             ! Z_dir=sdim-1
             ! dirMAC=1
             ! r=h(z)
            ZZgrid=xstenMAC(0,Z_dir+1)-problo(Z_dir+1)
            dz_external=(probhi(Z_dir+1)-problo(Z_dir+1))/coflow_num_cells
             ! j_external * dz_external = ZZgrid
            j_external=NINT(ZZgrid/dz_external)  ! round to nearest whole int. 
            if ((j_external.lt.0).or.(j_external.gt.coflow_num_cells)) then
             print *,"j_external invalid"
             stop
            endif
            coflow_Z(j_external)=j_external*dz_external

            if ((ls_below*ls_above.le.zero).and. &
                (abs(ls_below)+abs(ls_above).ge.EPS2*dx(1))) then

             if (ls_below.eq.zero) then
              xcrit=xbottom
             else if (ls_above.eq.zero) then
              xcrit=xtop
             else  ! phi-phi0=slope(x-x0) x=-phi0/slope+x0
              xcrit=xbottom-ls_below*(xtop-xbottom)/ &
                 (ls_above-ls_below)
             endif

             coflow_R_of_Z(j_external)=xcrit
            endif ! interface found
           enddo ! iside
          endif  ! mask>0
         enddo
         enddo

        else
         print *,"Z_dir,R_dir option not supported"
         stop
        endif

       else if (SDIM.eq.3) then

         ! z=h(x)
        if ((Z_dir.eq.0).and.(R_dir.eq.SDIM-1)) then
       
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
            dir=R_dir+1  ! R_dir=sdim-1  dir=sdim  z=h(x)
            im=1  ! check liquid height
            vofcomp=(im-1)*ngeom_recon+1

            do iside=0,1
             dirMAC=0 ! x direction
             call gridstenMAC(xstenMAC,xlo,i,j,k,fablo,bfact,dx,nhalf,dirMAC)

              ! R_dir=sdim-1  dir=sdim  z=h(x)
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

              ! Z_dir=0 R_dir=sdim-1 dir=sdim  z=h(x) 
              ! coordinate along streamwise direction starts at 0
             ZZgrid=xstenMAC(0,Z_dir+1)-problo(Z_dir+1)
             dz_external=(probhi(Z_dir+1)-problo(Z_dir+1))/coflow_num_cells
              ! j_external * dz_external = ZZgrid
             j_external=NINT(ZZgrid/dz_external)  ! round to nearest whole int. 
             if ((j_external.lt.0).or.(j_external.gt.coflow_num_cells)) then
              print *,"j_external invalid"
              stop
             endif
             coflow_Z(j_external)=j_external*dz_external

             use_vof_height=1

              ! R_dir=sdim-1  dir=sdim  z=h(x)
             if (use_vof_height.eq.0) then

              if ((ls_below*ls_above.le.zero).and. &
                  (abs(ls_below)+abs(ls_above).ge.EPS2*dx(1))) then

               if (ls_below.eq.zero) then
                xcrit=xbottom
               else if (ls_above.eq.zero) then
                xcrit=xtop
               else  ! phi-phi0=slope(x-x0) x=-phi0/slope+x0
                xcrit=xbottom-ls_below*(xtop-xbottom)/ &
                  (ls_above-ls_below)
               endif
 
               coflow_R_of_Z(j_external)=xcrit
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
               coflow_R_of_Z(j_external)=xcrit
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
         print *,"Z_dir,R_dir option not supported"
         stop
        endif

       else
        print *,"dimension bust"
        stop
       endif 

      else if (coflow_num_cells.eq.0) then
       ! do nothing
      else
       print *,"coflow_num_cells invalid"
       stop
      endif

      return
      end subroutine fort_summass

      subroutine fort_get_number_regions(cpp_number_regions) &
      bind(c,name='fort_get_number_regions')

      use probcommon_module
      IMPLICIT NONE
      integer, INTENT(out) :: cpp_number_regions

      cpp_number_regions=number_of_source_regions

      return
      end subroutine fort_get_number_regions

      subroutine fort_get_region_data(isweep, &
        cpp_energy_per_kelvin, &
        cpp_mass, &
        cpp_energy, &
        cpp_volume, &
        cpp_volume_raster, & ! this and above, isweep==0
        cpp_mass_after, &
        cpp_energy_after, &
        cpp_volume_after) &
      bind(c,name='fort_get_region_data')
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: isweep
      real(amrex_real), INTENT(inout) :: cpp_energy_per_kelvin(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_mass(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_energy(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_volume(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_volume_raster(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_mass_after(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_energy_after(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_volume_after(number_of_source_regions)
      integer iregions

      do iregions=1,number_of_source_regions
       if (isweep.eq.0) then
        cpp_energy_per_kelvin(iregions)= &
              regions_list(iregions,0)%region_energy_per_kelvin
        cpp_mass(iregions)= &
              regions_list(iregions,0)%region_mass
        cpp_energy(iregions)= &
              regions_list(iregions,0)%region_energy
        cpp_volume(iregions)= &
              regions_list(iregions,0)%region_volume
        cpp_volume_raster(iregions)= &
              regions_list(iregions,0)%region_volume_raster
       else if (isweep.eq.1) then
        cpp_mass_after(iregions)= &
              regions_list(iregions,0)%region_mass_after
        cpp_energy_after(iregions)= &
              regions_list(iregions,0)%region_energy_after
        cpp_volume_after(iregions)= &
              regions_list(iregions,0)%region_volume_after
       else
        print *,"isweep invalid"
        stop
       endif
      enddo ! iregions=1,number_source_regions
        
      end subroutine fort_get_region_data


      subroutine fort_put_region_data(isweep, &
        cpp_energy_per_kelvin, &
        cpp_mass, &
        cpp_energy, &
        cpp_volume, &
        cpp_volume_raster, & ! this and above, isweep==0
        cpp_mass_after, &
        cpp_energy_after, &
        cpp_volume_after) &
      bind(c,name='fort_put_region_data')
      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: isweep
      real(amrex_real), INTENT(inout) :: cpp_energy_per_kelvin(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_mass(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_energy(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_volume(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_volume_raster(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_mass_after(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_energy_after(number_of_source_regions)
      real(amrex_real), INTENT(inout) :: cpp_volume_after(number_of_source_regions)
      integer iregions

      do iregions=1,number_of_source_regions
       if (isweep.eq.0) then
        regions_list(iregions,0)%region_energy_per_kelvin= &
         cpp_energy_per_kelvin(iregions)
        regions_list(iregions,0)%region_mass= &
         cpp_mass(iregions)
        regions_list(iregions,0)%region_energy= &
         cpp_energy(iregions)
        regions_list(iregions,0)%region_volume= &
         cpp_volume(iregions)
        regions_list(iregions,0)%region_volume_raster= &
         cpp_volume_raster(iregions)
       else if (isweep.eq.1) then
        regions_list(iregions,0)%region_mass_after= &
         cpp_mass_after(iregions)
        regions_list(iregions,0)%region_energy_after= &
         cpp_energy_after(iregions)
        regions_list(iregions,0)%region_volume_after= &
         cpp_volume_after(iregions)
       else
        print *,"isweep invalid"
        stop
       endif
      enddo ! iregions=1,number_source_regions
        
      end subroutine fort_put_region_data


      subroutine fort_reduce_sum_regions(isweep) &
      bind(c,name='fort_reduce_sum_regions')

      use probcommon_module
      IMPLICIT NONE

      integer, INTENT(in) :: isweep

      integer ithreads,iregions

      do iregions=1,number_of_source_regions
       if (isweep.eq.0) then

        regions_list(iregions,0)%region_energy_per_kelvin=zero
        regions_list(iregions,0)%region_mass=zero
        regions_list(iregions,0)%region_energy=zero
        regions_list(iregions,0)%region_volume=zero
        regions_list(iregions,0)%region_volume_raster=zero

        do ithreads=1,number_of_threads_regions
         regions_list(iregions,0)%region_energy_per_kelvin= &
          regions_list(iregions,0)%region_energy_per_kelvin+ &
          regions_list(iregions,ithreads)%region_energy_per_kelvin
         regions_list(iregions,0)%region_mass= &
          regions_list(iregions,0)%region_mass+ &
          regions_list(iregions,ithreads)%region_mass
         regions_list(iregions,0)%region_energy= &
          regions_list(iregions,0)%region_energy+ &
          regions_list(iregions,ithreads)%region_energy
         regions_list(iregions,0)%region_volume= &
          regions_list(iregions,0)%region_volume+ &
          regions_list(iregions,ithreads)%region_volume
         regions_list(iregions,0)%region_volume_raster= &
          regions_list(iregions,0)%region_volume_raster+ &
          regions_list(iregions,ithreads)%region_volume_raster
        enddo ! ithreads=1,number_of_threads_regions

       else if (isweep.eq.1) then

        regions_list(iregions,0)%region_mass_after=zero
        regions_list(iregions,0)%region_energy_after=zero
        regions_list(iregions,0)%region_volume_after=zero

        do ithreads=1,number_of_threads_regions
         regions_list(iregions,0)%region_mass_after= &
          regions_list(iregions,0)%region_mass_after+ &
          regions_list(iregions,ithreads)%region_mass_after
         regions_list(iregions,0)%region_energy_after= &
          regions_list(iregions,0)%region_energy_after+ &
          regions_list(iregions,ithreads)%region_energy_after
         regions_list(iregions,0)%region_volume_after= &
          regions_list(iregions,0)%region_volume_after+ &
          regions_list(iregions,ithreads)%region_volume_after
        enddo ! ithreads=1,number_of_threads_regions

       else
        print *,"isweep invalid"
        stop
       endif

      enddo ! iregions=1,number_source_regions
        
      end subroutine fort_reduce_sum_regions

      subroutine fort_regionsum( &
       tid_current, &
       isweep, &  ! isweep=0 or 1
       constant_density_all_time, & ! 1..num_materials
       cur_time, &
       dt, &
       dx, &
       xlo, &
       nstate, &
       snew,DIMS(snew), &
       umacnew,DIMS(umacnew), &
       vmacnew,DIMS(vmacnew), &
       wmacnew,DIMS(wmacnew), &
       mdot, &
       DIMS(mdot), &
       DEN,DIMS(DEN), &
       VOF,DIMS(VOF), &
       volumefab, &
       DIMS(volumefab), &
       mask,DIMS(mask), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level) &
      bind(c,name='fort_regionsum')

      use probcommon_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid_current
      integer, INTENT(in) :: isweep ! isweep=0 or 1
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      integer, INTENT(in) :: constant_density_all_time(num_materials)

      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer :: growlo(3), growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: DIMDEC(snew)
      integer, INTENT(in) :: DIMDEC(umacnew)
      integer, INTENT(in) :: DIMDEC(vmacnew)
      integer, INTENT(in) :: DIMDEC(wmacnew)
      integer, INTENT(in) :: DIMDEC(mdot)
      integer, INTENT(in) :: DIMDEC(DEN)
      integer, INTENT(in) :: DIMDEC(VOF)
      integer, INTENT(in) :: DIMDEC(volumefab)
      integer, INTENT(in) :: DIMDEC(mask)

      real(amrex_real), INTENT(inout), target :: snew(DIMV(snew),nstate)
      real(amrex_real), pointer :: snew_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(inout), target :: umacnew(DIMV(umacnew))
      real(amrex_real), pointer :: umacnew_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: vmacnew(DIMV(vmacnew))
      real(amrex_real), pointer :: vmacnew_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout), target :: wmacnew(DIMV(wmacnew))
      real(amrex_real), pointer :: wmacnew_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(inout), target :: mdot(DIMV(mdot))
      real(amrex_real), pointer :: mdot_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in), target :: &
              DEN(DIMV(DEN),num_materials*num_state_material)
      real(amrex_real), pointer :: DEN_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: &
              VOF(DIMV(VOF),num_materials*ngeom_recon)
      real(amrex_real), pointer :: VOF_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: volumefab(DIMV(volumefab))
      real(amrex_real), pointer :: volumefab_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))

      integer, parameter :: nhalf=3

      integer :: i,j,k
      integer :: im
      integer :: dir
      integer :: local_dir
      integer :: local_mask,local_mask_L
      integer :: tessellate
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_L(-nhalf:nhalf,SDIM)
      real(amrex_real) xsten_R(-nhalf:nhalf,SDIM)
      real(amrex_real) xtarget(SDIM)
      integer iregions
      integer ispec
      integer vofcomp
      real(amrex_real) charfn
      real(amrex_real) vfrac,vfrac_L
      real(amrex_real) vfrac_raster
      real(amrex_real) region_mass_flux
      real(amrex_real) region_volume_flux
      real(amrex_real) region_energy_flux
      real(amrex_real) region_energy_per_kelvin
      real(amrex_real) region_mass
      real(amrex_real) region_volume
      real(amrex_real) region_volume_raster
      real(amrex_real) local_den
      real(amrex_real) local_temp
      real(amrex_real) DeDT
      real(amrex_real) Tflux
      real(amrex_real) density_new
      real(amrex_real) mass_new
      real(amrex_real) volume_new
      real(amrex_real) temperature_new
      real(amrex_real) temperature_prescribe
      real(amrex_real) density_flux
      real(amrex_real) divu
      integer update_density_flag
      real(amrex_real) massfrac_parm(num_species_var+1)
      integer imattype
      integer dencomp
      integer temperature_comp
      integer species_comp
      integer ii,jj,kk
      real(amrex_real) mofdata(num_materials*ngeom_recon)
      real(amrex_real) mofdata_L(num_materials*ngeom_recon)
      real(amrex_real) local_cellvol
      integer, PARAMETER :: nmax=POLYGON_LIST_MAX ! in: fort_regionsum

      if ((tid_current.lt.0).or.(tid_current.ge.geom_nthreads)) then
       print *,"tid_current invalid"
       stop
      endif

      if (number_of_source_regions.eq.0) then
       ! do nothing
      else if (number_of_source_regions.gt.0) then
       if (number_of_threads_regions.eq.geom_nthreads) then
        ! do nothing
       else
        print *,"number_of_threads_regions invalid"
        stop
       endif
      else
       print *,"number_of_source_regions invalid"
       stop
      endif

      snew_ptr=>snew
      umacnew_ptr=>umacnew
      vmacnew_ptr=>vmacnew
      wmacnew_ptr=>wmacnew
      mdot_ptr=>mdot
      DEN_ptr=>DEN
      VOF_ptr=>VOF
      mask_ptr=>mask
      volumefab_ptr=>volumefab

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
       print *,"level invalid fort_regionsum"
       stop
      endif

      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid in NAVIERSTOKES_3D.F90 "
       print *,"nstate=",nstate
       print *,"STATE_NCOMP=",STATE_NCOMP
       stop
      endif

      call checkbound_array(fablo,fabhi,snew_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,umacnew_ptr,0,0)
      call checkbound_array1(fablo,fabhi,vmacnew_ptr,0,1)
      call checkbound_array1(fablo,fabhi,wmacnew_ptr,0,SDIM-1)

      call checkbound_array1(fablo,fabhi,mdot_ptr,0,-1)
      call checkbound_array(fablo,fabhi,DEN_ptr,1,-1)
      call checkbound_array(fablo,fabhi,VOF_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,mask_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,volumefab_ptr,1,-1)
  
      do dir=1,SDIM
       if (fabhi(dir)-fablo(dir).le.0) then
        print *,"fablo,fabhi violates blocking factor"
        stop
       endif
      enddo

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       local_mask=NINT(mask(D_DECL(i,j,k)))
       if (local_mask.eq.1) then

        local_cellvol=volumefab(D_DECL(i,j,k))

        if (local_cellvol.gt.zero) then

         call gridsten_level(xsten,i,j,k,level,nhalf)

         do im=1,num_materials*ngeom_recon
          mofdata(im)=VOF(D_DECL(i,j,k),im)
         enddo

         tessellate=1 ! tessellate output (tessellate=3 => raster output)
          ! EPS2
         call multi_get_volume_tessellate( &
           tessellate, & ! tessellate=1 
           bfact, &
           dx, &
           xsten,nhalf, &
           mofdata, &
           geom_xtetlist(1,1,1,tid_current+1), &
           nmax, &
           nmax, &
           SDIM)

         do dir=1,SDIM
          xtarget(dir)=xsten(0,dir)
         enddo

         do iregions=1,number_of_source_regions
          im=regions_list(iregions,0)%region_material_id
          if ((im.ge.1).and.(im.le.num_materials)) then

           regions_list(iregions,0)%region_dt=dt

           vofcomp=(im-1)*ngeom_recon+1
           vfrac=mofdata(vofcomp)

           if ((vfrac.ge.zero).and.(vfrac.le.one)) then
            if (vfrac.le.one-VOFTOL) then
             vfrac_raster=zero
            else if (vfrac.ge.one-VOFTOL) then
             vfrac_raster=one
            else
             print *,"vfrac bust"
             stop
            endif
            call SUB_CHARFN_REGION(iregions,xtarget,cur_time,charfn)
            if ((charfn.eq.zero).or.(charfn.eq.one)) then

             imattype=fort_material_type(im)

             dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
             temperature_comp=(im-1)*num_state_material+1+ENUM_TEMPERATUREVAR
             species_comp=(im-1)*num_state_material+1+ENUM_SPECIESVAR

             local_den=DEN(D_DECL(i,j,k),dencomp)
             local_temp=DEN(D_DECL(i,j,k),temperature_comp)

             call init_massfrac_parm(local_den,massfrac_parm,im)
             do ispec=1,num_species_var
               massfrac_parm(ispec)=DEN(D_DECL(i,j,k),species_comp-1+ispec)
               if (massfrac_parm(ispec).ge.zero) then
                ! do nothing
               else
                print *,"massfrac_parm(ispec) invalid"
                stop
               endif
             enddo ! ispec=1,num_species_var

               ! DeDT = cv  ( J/(kg Kelvin)  )
             call DeDT_material(local_den,massfrac_parm, &
                local_temp,DeDT,imattype,im)

             if (DeDT.gt.zero) then
               ! do nothing
             else
               print *,"DeDT must be positive"
               stop
             endif

             if (isweep.eq.0) then
              regions_list(iregions,tid_current+1)%region_volume_raster= &
               regions_list(iregions,tid_current+1)%region_volume_raster+ &
               local_cellvol*vfrac_raster*charfn
              regions_list(iregions,tid_current+1)%region_energy_per_kelvin= &
               regions_list(iregions,tid_current+1)%region_energy_per_kelvin+ &
               local_cellvol*vfrac*charfn*local_den*DeDT
              regions_list(iregions,tid_current+1)%region_energy= &
               regions_list(iregions,tid_current+1)%region_energy+ &
               local_cellvol*vfrac*charfn*local_den*DeDT*local_temp
              regions_list(iregions,tid_current+1)%region_mass= &
               regions_list(iregions,tid_current+1)%region_mass+ &
               local_cellvol*vfrac*charfn*local_den
              regions_list(iregions,tid_current+1)%region_volume= &
               regions_list(iregions,tid_current+1)%region_volume+ &
               local_cellvol*vfrac*charfn

             else if (isweep.eq.1) then

              region_mass_flux=regions_list(iregions,0)%region_mass_flux
              region_volume_flux=regions_list(iregions,0)%region_volume_flux
              region_energy_flux=regions_list(iregions,0)%region_energy_flux
              region_volume_raster=regions_list(iregions,0)%region_volume_raster
              region_volume=regions_list(iregions,0)%region_volume
              region_mass=regions_list(iregions,0)%region_mass
              region_energy_per_kelvin= &
                   regions_list(iregions,0)%region_energy_per_kelvin
               ! mass_new = density_new * volume_new
               ! mass_old+dt * mass_flux = 
               !   density_new * (volume_old+ dt * volume_flux)
              volume_new=region_volume+dt*region_volume_flux
              if (volume_new.gt.zero) then
               mass_new=region_mass+dt*region_mass_flux
               if (mass_new.gt.zero) then
                density_new=mass_new/volume_new
                if (density_new.gt.zero) then
                 ! do nothing
                else
                 print *,"density_new invalid"
                 stop
                endif
               else
                print *,"mass_new invalid"
                stop
               endif
              else if (volume_new.eq.zero) then
               ! do nothing
              else
               print *,"volume_new invalid (1)",volume_new
               print *,"iregions=",iregions
               print *,"tid_current=",tid_current
               print *,"region_volume ",region_volume
               print *,"region_volume_flux ",region_volume_flux
               print *,"dt ",dt
               do dir=1,SDIM
                print *,"dir,dx,dx/(2 dt) ",dir,dx(dir),dx(dir)/(two*dt)
               enddo
               print *,"im ",im
               stop
              endif

               ! volume_flux=sum_p (div u)_p xi(x_p) F_raster_p vol_p
               ! assume div u is spatially uniform:
               ! div u =volume_flux/volume_raster  (units 1/seconds)
               ! units of mdot: cm^3/second^2
              if (region_volume_flux.ne.zero) then
               if (is_rigid(im).eq.1) then
                print *,"disallowed: volume_flux<>0 for is_rigid==1 material"
                stop
               else if (is_rigid(im).eq.0) then
                if ((imattype.gt.0).and.(imattype.lt.999)) then
                 print *,"disallowed: volume flux<>0 for compressible material"
                 stop
                else if (imattype.eq.0) then
                 ! do nothing
                else
                 print *,"imattype invalid"
                 stop
                endif
               else
                print *,"is_rigid invalid NAVIERSTOKES_3D.F90"
                stop
               endif

               if (region_volume_raster.gt.zero) then 
                if (charfn.eq.one) then
                 divu=region_volume_flux/region_volume_raster
                 mdot(D_DECL(i,j,k))=mdot(D_DECL(i,j,k))+ &
                  divu*local_cellvol*charfn*vfrac_raster/dt

                 temperature_prescribe= &
                   regions_list(iregions,0)%region_temperature_prescribe

                 if (region_volume_flux.lt.zero) then
                  if (temperature_prescribe.eq.zero) then
                   ! do nothing
                  else
                   print *,"cannot prescribe temperature if volume_flux<0"
                   stop
                  endif
                 else if (region_volume_flux.gt.zero) then

                  if (vfrac_raster.eq.zero) then
                   ! do nothing
                  else if (vfrac_raster.eq.one) then
                   do dir=1,SDIM
                    snew(D_DECL(i,j,k),dir)= &
                     regions_list(iregions,0)%region_velocity_prescribe(dir)
                   enddo
                  else
                   print *,"vfrac_raster invalid"
                   stop
                  endif

                  if (temperature_prescribe.eq.zero) then
                   ! do nothing
                  else if (temperature_prescribe.gt.zero) then
                   if (vfrac.eq.zero) then
                    ! do nothing
                   else if ((vfrac.gt.zero).and. &
                            (vfrac.le.one+EPS1)) then
                    snew(D_DECL(i,j,k), &
                     STATECOMP_STATES+temperature_comp)=temperature_prescribe
                   else
                    print *,"vfrac invalid: ",vfrac
                    stop
                   endif
                  else
                   print *,"temperature_prescribe invalid"
                   stop
                  endif

                 else
                  print *,"region_volume_flux bust"
                  stop
                 endif

                else if (charfn.eq.zero) then
                 ! do nothing
                else
                 print *,"charfn invalid"
                 stop
                endif
               else if (region_volume_raster.eq.zero) then
                ! do nothing
               else
                print *,"region_volume_raster invalid"
                stop
               endif
              else if (region_volume_flux.eq.zero) then
               divu=zero
              else
               print *,"region_volume_flux invalid"
               stop
              endif

              update_density_flag=0
              density_flux=zero
              if ((region_volume_flux.eq.zero).and. &
                  (region_mass_flux.eq.zero)) then
               ! do nothing
              else if ((region_volume_flux.eq.zero).and. &
                       (region_mass_flux.ne.zero)) then
               update_density_flag=1

               if (is_rigid(im).eq.1) then
                print *,"mass_flux<>0 disallowed for is_rigid==1 materials"
                stop
               else if (is_rigid(im).eq.0) then
                if ((imattype.gt.0).and.(imattype.lt.999)) then
                 ! do nothing
                else if (imattype.eq.0) then
                 if (constant_density_all_time(im).eq.0) then
                  ! do nothing
                 else
                  print *,"constant_density_all_time invalid"
                  stop
                 endif
                else
                 print *,"imattype invalid"
                 stop
                endif
               else
                print *,"is_rigid(im) invalid"
                stop
               endif
  
              else if (region_volume_flux.ne.zero) then
               density_flux=region_mass_flux/region_volume_flux
               if (is_rigid(im).eq.1) then
                print *,"volume_flux<>0 disallowed for is_rigid==1 materials"
                stop
               else if (is_rigid(im).eq.0) then
                if ((imattype.gt.0).and.(imattype.lt.999)) then
                 print *,"disallowed: volume flux<>0 for compressible material"
                 stop
                else if (imattype.eq.0) then
                 if (abs(fort_denconst(im)-density_flux).le.VOFTOL) then
                  ! do nothing
                 else if (abs(fort_denconst(im)-density_flux).gt.VOFTOL) then
                  update_density_flag=1
                  if (constant_density_all_time(im).eq.0) then
                   ! do nothing
                  else
                   print *,"constant_density_all_time invalid"
                   stop
                  endif
                 else
                  print *,"density_flux out of range"
                  stop
                 endif
                else
                 print *,"imattype invalid"
                 stop
                endif
               else
                print *,"is_rigid(im) invalid"
                stop
               endif
              else
               print *,"region_volume_flux invalid"
               stop
              endif
              if (update_density_flag.eq.1) then
               snew(D_DECL(i,j,k),STATECOMP_STATES+dencomp)=density_new
              else if (update_density_flag.eq.0) then
               ! do nothing
              else
               print *,"update_density_flag invalid"
               stop
              endif
               ! energy_flux dimensions: J/second=Watts
               ! cv dimensions: J/(kg Kelvin)
               ! rho cv T dimensions: kg/m^3   *  J/(kg Kelvin)  *  Kelvin=
               ! J/m^3
               ! energy_flux=sum_p den_p * cv_p * Tflux * charfn_p * F_p * Vol_p
               !  units of RHS: kg/m^3 * J/(kg Kelvin) * (Kelvin/sec) * m^3=
               !    J/sec
               ! energy_per_kelvin=sum_p den_p cv_p charfn_p F_p vol_p (J/K)
               !       
               ! Tflux=energy_flux/(energy_per_kelvin)=(J/s)/(J/K)=Kelvin/s

              temperature_new=local_temp
              if (region_energy_flux.ne.zero) then
               if (region_energy_per_kelvin.gt.zero) then
                if (local_den.gt.zero) then
                 region_energy_per_kelvin=(density_new/local_den)* &
                     region_energy_per_kelvin
                else
                 print *,"local_den invalid"
                 stop
                endif
                Tflux=region_energy_flux/region_energy_per_kelvin
                if (Tflux.ne.zero) then
                 temperature_new=local_temp
                 if (vfrac.gt.zero) then
                  temperature_new=temperature_new+dt*Tflux*charfn
                 else if (vfrac.eq.zero) then
                  ! do nothing
                 else
                  print *,"vfrac invalid"
                  stop
                 endif
                 
                 if (temperature_new.gt.zero) then
                  snew(D_DECL(i,j,k),STATECOMP_STATES+temperature_comp)= &
                    temperature_new
                 else
                  print *,"temperature_new invalid"
                  stop
                 endif
                else
                 print *,"Tflux invalid"
                 stop
                endif
               else
                print *,"region_energy_per_kelvin invalid"
                stop
               endif
              else if (region_energy_flux.eq.zero) then
               ! do nothing
              else
               print *,"region_energy_flux bust"
               stop
              endif
              if ((region_volume.gt.zero).and.(volume_new.ge.zero)) then

               regions_list(iregions,tid_current+1)%region_volume_after= &
                regions_list(iregions,tid_current+1)%region_volume_after+ &
                local_cellvol*vfrac*charfn* &
                (volume_new/region_volume)
               regions_list(iregions,tid_current+1)%region_mass_after= &
                regions_list(iregions,tid_current+1)%region_mass_after+ &
                local_cellvol*vfrac*charfn*density_new* &
                (volume_new/region_volume)
               regions_list(iregions,tid_current+1)%region_energy_after= &
                regions_list(iregions,tid_current+1)%region_energy_after+ &
                local_cellvol*vfrac*charfn*density_new*DeDT* &
                temperature_new* &
                (volume_new/region_volume)

              else if (region_volume.eq.zero) then
               if (volume_new.eq.zero) then
                ! do nothing
               else
                print *,"volume_new invalid (2)",volume_new
                print *,"iregions=",iregions
                print *,"tid_current=",tid_current
                print *,"region_volume ",region_volume
                print *,"region_volume_flux ",region_volume_flux
                print *,"dt ",dt
                print *,"im ",im
                stop
               endif
              else
               print *,"region_volume invalid"
               stop
              endif

             else
              print *,"isweep invalid"
              stop
             endif

            else
             print *,"charfn invalid"
             stop
            endif

           else
            print *,"vfrac out of range"
            stop
           endif
          else
           print *,"im invalid 11385:",im
           stop
          endif
         enddo !iregions=1,number_of_source_regions

        else
         print *,"local_cellvol invalid"
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

      if (isweep.eq.0) then
       ! do nothing
      else if (isweep.eq.1) then

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

        call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir-1) 

        do k=growlo(3),growhi(3)
        do j=growlo(2),growhi(2)
        do i=growlo(1),growhi(1)

         local_mask=NINT(mask(D_DECL(i,j,k)))
         local_mask_L=NINT(mask(D_DECL(i-ii,j-jj,k-kk)))
         if ((local_mask.eq.1).or.(local_mask_L.eq.1)) then

          call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir-1)
          call gridsten_level(xsten_R,i,j,k,level,nhalf)
          call gridsten_level(xsten_L,i-ii,j-jj,k-kk,level,nhalf)

          do im=1,num_materials*ngeom_recon
           mofdata(im)=VOF(D_DECL(i,j,k),im)
           mofdata_L(im)=VOF(D_DECL(i-ii,j-jj,k-kk),im)
          enddo

          tessellate=1 ! tessellate output (tessellate=3 => raster output)
           ! EPS2
          call multi_get_volume_tessellate( &
           tessellate, & ! tessellate=1 
           bfact, &
           dx, &
           xsten_R,nhalf, &
           mofdata, &
           geom_xtetlist(1,1,1,tid_current+1), &
           nmax, &
           nmax, &
           SDIM)

           ! EPS2
          call multi_get_volume_tessellate( &
           tessellate, & ! tessellate=1 
           bfact, &
           dx, &
           xsten_L,nhalf, &
           mofdata_L, &
           geom_xtetlist(1,1,1,tid_current+1), &
           nmax, &
           nmax, &
           SDIM)

          do local_dir=1,SDIM
           xtarget(local_dir)=xsten(0,local_dir)
          enddo

          do iregions=1,number_of_source_regions
           im=regions_list(iregions,0)%region_material_id
           if ((im.ge.1).and.(im.le.num_materials)) then

            vofcomp=(im-1)*ngeom_recon+1
            vfrac=mofdata(vofcomp)
            vfrac_L=mofdata_L(vofcomp)
            if ((vfrac.ge.zero).and.(vfrac.le.one).and. &
                (vfrac_L.ge.zero).and.(vfrac_L.le.one)) then
             if ((vfrac.le.one-VOFTOL).and. &
                 (vfrac_L.le.one-VOFTOL)) then
              vfrac_raster=zero
             else if ((vfrac.ge.one-VOFTOL).or. &
                      (vfrac_L.ge.one-VOFTOL)) then
              vfrac_raster=one
             else
              print *,"vfrac bust"
              stop
             endif
             call SUB_CHARFN_REGION(iregions,xtarget,cur_time,charfn)
             if ((charfn.eq.zero).or.(charfn.eq.one)) then

              imattype=fort_material_type(im)

              region_volume_flux=regions_list(iregions,0)%region_volume_flux
              region_volume_raster=regions_list(iregions,0)%region_volume_raster

              if (region_volume_flux.lt.zero) then
               ! do nothing
              else if (region_volume_flux.eq.zero) then
               ! do nothing
              else if (region_volume_flux.gt.zero) then

               if (is_rigid(im).eq.1) then
                print *,"disallowed: volume_flux>0 for is_rigid==1 material"
                stop
               else if (is_rigid(im).eq.0) then
                if ((imattype.gt.0).and.(imattype.lt.999)) then
                 print *,"disallowed: volume flux>0 for compressible material"
                else if (imattype.eq.0) then
                 ! do nothing
                else
                 print *,"imattype invalid"
                 stop
                endif
               else
                print *,"is_rigid invalid NAVIERSTOKES_3D.F90"
                stop
               endif

               if (region_volume_raster.gt.zero) then 
                if (charfn.eq.one) then
                 if (vfrac_raster.eq.zero) then
                  ! do nothing
                 else if (vfrac_raster.eq.one) then
                  if (dir.eq.1) then
                   umacnew(D_DECL(i,j,k))= &
                    regions_list(iregions,0)%region_velocity_prescribe(dir)
                  else if (dir.eq.2) then
                   vmacnew(D_DECL(i,j,k))= &
                    regions_list(iregions,0)%region_velocity_prescribe(dir)
                  else if ((dir.eq.3).and.(SDIM.eq.3)) then
                   wmacnew(D_DECL(i,j,k))= &
                    regions_list(iregions,0)%region_velocity_prescribe(dir)
                  else
                   print *,"dir invalid"
                   stop
                  endif
                 else
                  print *,"vfrac_raster invalid"
                  stop
                 endif
                else if (charfn.eq.zero) then
                 ! do nothing
                else
                 print *,"charfn invalid"
                 stop
                endif
               else if (region_volume_raster.eq.zero) then
                ! do nothing
               else
                print *,"region_volume_raster invalid"
                stop
               endif
              else
               print *,"region_volume_flux invalid"
               stop
              endif

             else
              print *,"charfn invalid"
              stop
             endif

            else
             print *,"vfrac or vfrac_L out of range"
             stop
            endif
           else
            print *,"im invalid 11582:",im
            stop
           endif
          enddo !iregions=1,number_of_source_regions

         else if ((local_mask.eq.0).and.(local_mask_L.eq.0)) then
          ! do nothing
         else
          print *,"local_mask or local_mask_L invalid"
          stop
         endif

        enddo
        enddo
        enddo

       enddo ! dir=1..sdim

      else
       print *,"isweep invalid"
       stop
      endif

      return
      end subroutine fort_regionsum

       !! fort_fabcom: fabz = fabx + beta * faby
       !! Added by Alan Kuhnle, 6-7-10
       subroutine fort_fabcom( &
        fabx,DIMS(fabx), &
        faby, DIMS(faby), &
        mask, DIMS(mask), &
        fabz, DIMS(fabz), &
        beta, &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        nsolve) &
       bind(c,name='fort_fabcom')

       use global_utility_module
       use probcommon_module

       IMPLICIT NONE

       integer, INTENT(in) :: nsolve
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer :: growlo(3),growhi(3)
       integer, INTENT(in) :: bfact
       integer, INTENT(in) :: DIMDEC(fabx)
       integer, INTENT(in) :: DIMDEC(faby)
       integer, INTENT(in) :: DIMDEC(fabz)
       integer, INTENT(in) :: DIMDEC(mask)
       real(amrex_real), INTENT(in) :: beta
       real(amrex_real), INTENT(in), target :: fabx(DIMV(fabx),nsolve)
       real(amrex_real), INTENT(in), target :: faby(DIMV(faby),nsolve)

       real(amrex_real), INTENT(out), target :: fabz(DIMV(fabz),nsolve)
       real(amrex_real), pointer :: fabz_ptr(D_DECL(:,:,:),:)

       real(amrex_real), INTENT(in), target :: mask(DIMV(mask))

       integer local_mask
       integer :: nc

       integer i,j,k

       fabz_ptr=>fabz

       if (bfact.lt.1) then
        print *,"bfact invalid157"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid"
        stop
       endif

       call checkbound_array1(fablo,fabhi,mask,0,-1) 
       call checkbound_array(fablo,fabhi,fabx,0,-1) 
       call checkbound_array(fablo,fabhi,faby,0,-1) 
       call checkbound_array(fablo,fabhi,fabz_ptr,0,-1) 

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do nc=1,nsolve

         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          local_mask=NINT(mask(D_DECL(i,j,k)))

          if (local_mask.eq.1) then

            fabz(D_DECL(i,j,k),nc) = fabx(D_DECL(i,j,k),nc) +  &
             beta*faby(D_DECL(i,j,k),nc)

          else if (local_mask.eq.0) then
           ! do nothing
          else
           print *,"mask invalid"
           stop
          endif

         enddo !k
         enddo !j
         enddo !i

       enddo !nc=1..nsolve

       return
       end subroutine fort_fabcom


       subroutine fort_sumdot(mass1,  &
        rho,DIMS(rho), &
        rho2,DIMS(rho2), &
        mask,DIMS(mask), &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        debug_dot_product, &
        levelno,gridno, &
        nsolve) &
       bind(c,name='fort_sumdot')
    
       use global_utility_module
       use probf90_module
 
       IMPLICIT NONE

       integer, INTENT(in) :: levelno,gridno
       integer, INTENT(in) :: nsolve
       integer, INTENT(in) :: debug_dot_product
       integer imax,jmax,kmax
       real(amrex_real) dotmax
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer growlo(3),growhi(3)
       integer, INTENT(in) :: bfact
       integer, INTENT(in) :: DIMDEC(rho)
       integer, INTENT(in) :: DIMDEC(rho2)
       integer, INTENT(in) :: DIMDEC(mask)
       real(amrex_real), INTENT(out) :: mass1
       real(amrex_real), INTENT(in), target :: rho(DIMV(rho),nsolve)
       real(amrex_real), INTENT(in), target :: rho2(DIMV(rho2),nsolve)
       real(amrex_real), INTENT(in), target :: mask(DIMV(mask))
       real(amrex_real) :: dm
       integer local_mask

       integer i,j,k,nc


       if ((levelno.lt.0).or.(gridno.lt.0)) then
        print *,"level or grid invalid"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid"
        stop
       endif
       call checkbound_array1(fablo,fabhi,mask,0,-1) 
       call checkbound_array(fablo,fabhi,rho,0,-1) 
       call checkbound_array(fablo,fabhi,rho2,0,-1) 
 
       mass1=zero

       imax=0
       jmax=0
       kmax=0
       dotmax=zero
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do nc=1,nsolve

         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          local_mask=NINT(mask(D_DECL(i,j,k)))

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

         enddo ! k
         enddo ! j
         enddo ! i

       enddo ! nc=1..nsolve

       if (debug_dot_product.eq.1) then
        print *,"debug dot: level,grid,i,j,k,max ",levelno,gridno, &
         imax,jmax,kmax,dotmax
       else if (debug_dot_product.ne.0) then
        print *,"debug dot prod invalid"
        stop
       endif

       return
       end subroutine fort_sumdot


       subroutine fort_sumdot_ones_size( &
        fab_sum,  &
        fab_flag,  &
        ones_fab, &
        DIMS(ones_fab), &
        type_fab, &
        DIMS(type_fab), &
        color_fab, &
        DIMS(color_fab), &
        alpha_fab, &
        DIMS(alpha_fab), &
        mask_fab, &
        DIMS(mask_fab), &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        levelno,gridno, &
        nsolve, &
        presbc, &
        type_flag, &
        color_count, &
        project_option) &
       bind(c,name='fort_sumdot_ones_size')
    
       use global_utility_module
       use probf90_module
 
       IMPLICIT NONE

       integer, INTENT(in) :: project_option
       integer, INTENT(in) :: color_count
       real(amrex_real), INTENT(inout) :: fab_sum(color_count)
       integer, INTENT(inout) :: fab_flag(color_count)
       integer, INTENT(in) :: levelno,gridno
       integer, INTENT(in) :: nsolve
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer growlo(3),growhi(3)
       integer, INTENT(in) :: bfact
       integer, INTENT(in) :: DIMDEC(ones_fab)
       integer, INTENT(in) :: DIMDEC(type_fab)
       integer, INTENT(in) :: DIMDEC(color_fab)
       integer, INTENT(in) :: DIMDEC(alpha_fab)
       integer, INTENT(in) :: DIMDEC(mask_fab)
       real(amrex_real), INTENT(in),target :: ones_fab(DIMV(ones_fab))
       real(amrex_real), INTENT(in), target :: type_fab(DIMV(type_fab))
       real(amrex_real), INTENT(in), target :: color_fab(DIMV(color_fab))
       real(amrex_real), INTENT(in), target :: alpha_fab(DIMV(alpha_fab),nsolve)
       real(amrex_real), INTENT(in), target :: mask_fab(DIMV(mask_fab))
       integer, INTENT(in) :: presbc(SDIM,2)
       integer, INTENT(in) :: type_flag(2)

       integer icolor,nc
       integer i,j,k
       integer local_mask,local_color,local_type,local_ones
       integer plus_flag,zero_flag
       real(amrex_real) local_alpha
       integer dir_local
       integer icrit
       integer side
       integer local_bc

       if ((levelno.lt.0).or.(gridno.lt.0)) then
        print *,"level or grid invalid"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid"
        stop
       endif
       if (project_option.eq.SOLVETYPE_VISC) then
        if (nsolve.eq.SDIM) then
         ! do nothing
        else
         print *,"nsolve invalid"
         stop
        endif
       else if (project_option_is_validF(project_option).eq.1) then
        if (nsolve.eq.1) then
         ! do nothing
        else
         print *,"nsolve invalid"
         stop
        endif
       else
        print *,"project_option invalid; fort_sumdot_ones_size"
        stop
       endif
       if (color_count.ge.1) then
        ! do nothing
       else
        print *,"color_count invalid"
        stop
       endif

       call checkbound_array1(fablo,fabhi,ones_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,type_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,color_fab,0,-1) 
       call checkbound_array(fablo,fabhi,alpha_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,mask_fab,0,-1) 

       do icolor=1,color_count
        fab_sum(icolor)=zero
        fab_flag(icolor)=0
       enddo 

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        local_mask=NINT(mask_fab(D_DECL(i,j,k)))

        if (local_mask.eq.1) then

         local_color=NINT(color_fab(D_DECL(i,j,k)))
         if ((local_color.ge.1).and.(local_color.le.color_count)) then
          fab_sum(local_color)=fab_sum(local_color)+1
          local_type=NINT(type_fab(D_DECL(i,j,k)))
          local_ones=NINT(ones_fab(D_DECL(i,j,k)))
           ! completely masked off cell
          if ((local_ones.eq.0).and.(local_type.eq.1)) then
           if (type_flag(local_type).eq.1) then
            ! do nothing
           else
            print *,"type_flag invalid"
            stop
           endif
          else if ((local_ones.eq.1).and.(local_type.eq.2)) then
           if (fab_flag(local_color).lt.1) then
            fab_flag(local_color)=1
           endif
           if (type_flag(local_type).eq.1) then
            ! do nothing
           else
            print *,"type_flag invalid"
            stop
           endif
           plus_flag=0
           zero_flag=0
           do nc=1,nsolve
            local_alpha=alpha_fab(D_DECL(i,j,k),nc)
            if (local_alpha.eq.zero) then
             zero_flag=1
            else if (local_alpha.gt.zero) then
             plus_flag=1
            else
             print *,"local_alpha invalid"
             stop
            endif 
           enddo !nc=1..nsolve
           if ((plus_flag.eq.1).and.(zero_flag.eq.0)) then
            fab_flag(local_color)=2
           else if ((plus_flag.eq.0).and.(zero_flag.eq.1)) then
            ! do nothing
           else
            print *,"plus_flag or zero_flag invalid"
            stop
           endif
           do dir_local=1,SDIM
            if (dir_local.eq.1) then
             icrit=i
            else if (dir_local.eq.2) then
             icrit=j
            else if ((dir_local.eq.3).and.(SDIM.eq.3)) then
             icrit=k
            else
             print *,"dir_local invalid"
             stop
            endif
            side=0
            if ((icrit.gt.fablo(dir_local)).and. &
                (icrit.lt.fabhi(dir_local))) then
             ! do nothing
            else if (icrit.eq.fablo(dir_local)) then
             side=1
            else if (icrit.eq.fabhi(dir_local)) then
             side=2
            else
             print *,"icrit invalid"
             stop
            endif
            if (side.eq.0) then
             ! do nothing
            else if ((side.eq.1).or.(side.eq.2)) then
             local_bc=presbc(dir_local,side)
             if (project_option.eq.SOLVETYPE_PRESEXTRAP) then 
              ! do nothing (all bcs are Neumann)
             else if &
               (project_option_singular_possibleF(project_option).eq.1) then
              if (local_bc.eq.INT_DIR) then
               ! do nothing
              else if (local_bc.eq.FOEXTRAP) then
               ! do nothing
              else if (local_bc.eq.REFLECT_EVEN) then
               ! do nothing
              else if (local_bc.eq.EXT_DIR) then
               fab_flag(local_color)=2
              else
               print *,"local_bc invalid"
               stop
              endif
             else if &
               (project_option_singular_possibleF(project_option).eq.0) then
              ! do nothing
             else
              print *,"project_option invalid; fort_sumdot_ones_size"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           enddo !dir_local=1..sdim
          else
           print *,"local_ones or local_type invalid"
           stop
          endif
         else
          print *,"local_color invalid"
          stop
         endif

        else if (local_mask.eq.0) then
         ! do nothing
        else 
         print *,"mask invalid"
         stop
        endif

       enddo ! k
       enddo ! j
       enddo ! i

       return
       end subroutine fort_sumdot_ones_size

       subroutine fort_sumdot_ones( &
        fab_sum,  &
        fab_flag,  &
        data_fab, &
        DIMS(data_fab), &
        ones_fab, &
        DIMS(ones_fab), &
        type_fab, &
        DIMS(type_fab), &
        color_fab, &
        DIMS(color_fab), &
        alpha_fab, &
        DIMS(alpha_fab), &
        mask_fab, &
        DIMS(mask_fab), &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        levelno, &
        gridno, &
        nsolve, &
        presbc, &
        type_flag, &
        color_count, &
        project_option) &
       bind(c,name='fort_sumdot_ones')
    
       use global_utility_module
       use probf90_module
 
       IMPLICIT NONE

       integer, INTENT(in) :: project_option
       integer, INTENT(in) :: color_count
       real(amrex_real), INTENT(inout) :: fab_sum(color_count)
       integer, INTENT(in) :: fab_flag(color_count)
       integer, INTENT(in) :: levelno,gridno
       integer, INTENT(in) :: nsolve
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer growlo(3),growhi(3)
       integer, INTENT(in) :: bfact
       integer, INTENT(in) :: DIMDEC(data_fab)
       integer, INTENT(in) :: DIMDEC(ones_fab)
       integer, INTENT(in) :: DIMDEC(type_fab)
       integer, INTENT(in) :: DIMDEC(color_fab)
       integer, INTENT(in) :: DIMDEC(alpha_fab)
       integer, INTENT(in) :: DIMDEC(mask_fab)
       real(amrex_real), INTENT(in), target :: data_fab(DIMV(data_fab))
       real(amrex_real), INTENT(in), target :: ones_fab(DIMV(ones_fab))
       real(amrex_real), INTENT(in), target :: type_fab(DIMV(type_fab))
       real(amrex_real), INTENT(in), target :: color_fab(DIMV(color_fab))
       real(amrex_real), INTENT(in), target :: alpha_fab(DIMV(alpha_fab),nsolve)
       real(amrex_real), INTENT(in), target :: mask_fab(DIMV(mask_fab))
       integer, INTENT(in) :: presbc(SDIM,2)
       integer, INTENT(in) :: type_flag(2)

       integer icolor,nc
       integer i,j,k
       integer local_mask,local_color,local_type,local_ones
       integer plus_flag,zero_flag
       real(amrex_real) local_alpha
       integer dir_local
       integer icrit
       integer side
       integer local_bc

       if ((levelno.lt.0).or.(gridno.lt.0)) then
        print *,"level or grid invalid"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid"
        stop
       endif
       if (project_option.eq.SOLVETYPE_VISC) then
        print *,"fort_sumdot_ones should not be called for viscosity"
        stop
        if (nsolve.eq.SDIM) then
         ! do nothing
        else
         print *,"nsolve invalid"
         stop
        endif
       else if (project_option_is_validF(project_option).eq.1) then
        if (nsolve.eq.1) then
         ! do nothing
        else
         print *,"nsolve invalid"
         stop
        endif
       else
        print *,"project_option invalid; fort_sumdot_ones"
        stop
       endif
       if (color_count.ge.1) then
        ! do nothing
       else
        print *,"color_count invalid"
        stop
       endif

       call checkbound_array1(fablo,fabhi,data_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,ones_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,type_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,color_fab,0,-1) 
       call checkbound_array(fablo,fabhi,alpha_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,mask_fab,0,-1) 

       do icolor=1,color_count
        fab_sum(icolor)=zero
       enddo 

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        local_mask=NINT(mask_fab(D_DECL(i,j,k)))

        if (local_mask.eq.1) then

         local_color=NINT(color_fab(D_DECL(i,j,k)))
         if ((local_color.ge.1).and.(local_color.le.color_count)) then
          fab_sum(local_color)=fab_sum(local_color)+ &
                  data_fab(D_DECL(i,j,k))
          local_type=NINT(type_fab(D_DECL(i,j,k)))
          local_ones=NINT(ones_fab(D_DECL(i,j,k)))
           ! completely masked off cell
          if ((local_ones.eq.0).and.(local_type.eq.1)) then
           if (type_flag(local_type).eq.1) then
            ! do nothing
           else
            print *,"type_flag invalid"
            stop
           endif
          else if ((local_ones.eq.1).and.(local_type.eq.2)) then
           if (fab_flag(local_color).ge.1) then
            ! do nothing
           else
            print *,"fab_flag invalid"
            stop
           endif
           if (type_flag(local_type).eq.1) then
            ! do nothing
           else
            print *,"type_flag invalid"
            stop
           endif
           plus_flag=0
           zero_flag=0
           do nc=1,nsolve
            local_alpha=alpha_fab(D_DECL(i,j,k),nc)
            if (local_alpha.eq.zero) then
             zero_flag=1
            else if (local_alpha.gt.zero) then
             plus_flag=1
            else
             print *,"local_alpha invalid"
             stop
            endif 
           enddo !nc=1..nsolve
           if ((plus_flag.eq.1).and.(zero_flag.eq.0)) then
            if (fab_flag(local_color).eq.2) then
             ! do nothing
            else
             print *,"fab_flag invalid"
             stop
            endif
           else if ((plus_flag.eq.0).and.(zero_flag.eq.1)) then
            ! do nothing
           else
            print *,"plus_flag or zero_flag invalid"
            stop
           endif
           do dir_local=1,SDIM
            if (dir_local.eq.1) then
             icrit=i
            else if (dir_local.eq.2) then
             icrit=j
            else if ((dir_local.eq.3).and.(SDIM.eq.3)) then
             icrit=k
            else
             print *,"dir_local invalid"
             stop
            endif
            side=0
            if ((icrit.gt.fablo(dir_local)).and. &
                (icrit.lt.fabhi(dir_local))) then
             ! do nothing
            else if (icrit.eq.fablo(dir_local)) then
             side=1
            else if (icrit.eq.fabhi(dir_local)) then
             side=2
            else
             print *,"icrit invalid"
             stop
            endif
            if (side.eq.0) then
             ! do nothing
            else if ((side.eq.1).or.(side.eq.2)) then
             local_bc=presbc(dir_local,side)
             if (project_option.eq.SOLVETYPE_PRESEXTRAP) then 
              ! do nothing (all bcs are Neumann)
             else if &
               (project_option_singular_possibleF(project_option).eq.1) then
              if (local_bc.eq.INT_DIR) then
               ! do nothing
              else if (local_bc.eq.FOEXTRAP) then
               ! do nothing
              else if (local_bc.eq.REFLECT_EVEN) then
               ! do nothing
              else if (local_bc.eq.EXT_DIR) then
               if (fab_flag(local_color).eq.2) then
                ! do nothing
               else
                print *,"fab_flag invalid"
                stop
               endif
              else
               print *,"local_bc invalid"
               stop
              endif
             else if &
               (project_option_singular_possibleF(project_option).eq.0) then
              ! do nothing
             else
              print *,"project_option invalid; fort_sumdot_ones"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           enddo !dir_local=1..sdim
          else
           print *,"local_ones or local_type invalid"
           stop
          endif
         else
          print *,"local_color invalid"
          stop
         endif

        else if (local_mask.eq.0) then
         ! do nothing
        else 
         print *,"mask invalid"
         stop
        endif

       enddo ! k
       enddo ! j
       enddo ! i

       return
       end subroutine fort_sumdot_ones

       subroutine fort_fabcom_ones( &
        beta,  &
        singular_patch_flag,  &
        data_fab, &
        DIMS(data_fab), &
        ones_fab, &
        DIMS(ones_fab), &
        type_fab, &
        DIMS(type_fab), &
        color_fab, &
        DIMS(color_fab), &
        alpha_fab, &
        DIMS(alpha_fab), &
        mask_fab, &
        DIMS(mask_fab), &
        tilelo,tilehi, &
        fablo,fabhi,bfact, &
        levelno, &
        gridno, &
        nsolve, &
        presbc, &
        type_flag, &
        color_count, &
        project_option) &
       bind(c,name='fort_fabcom_ones')
    
       use global_utility_module
       use probf90_module
 
       IMPLICIT NONE

       integer, INTENT(in) :: project_option
       integer, INTENT(in) :: color_count
       real(amrex_real), INTENT(in) :: beta(color_count)
       integer, INTENT(in) :: singular_patch_flag(color_count)
       integer, INTENT(in) :: levelno,gridno
       integer, INTENT(in) :: nsolve
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer growlo(3),growhi(3)
       integer, INTENT(in) :: bfact
       integer, INTENT(in) :: DIMDEC(data_fab)
       integer, INTENT(in) :: DIMDEC(ones_fab)
       integer, INTENT(in) :: DIMDEC(type_fab)
       integer, INTENT(in) :: DIMDEC(color_fab)
       integer, INTENT(in) :: DIMDEC(alpha_fab)
       integer, INTENT(in) :: DIMDEC(mask_fab)

       real(amrex_real), INTENT(inout), target :: data_fab(DIMV(data_fab))
       real(amrex_real), pointer :: data_fab_ptr(D_DECL(:,:,:))

       real(amrex_real), INTENT(in), target :: ones_fab(DIMV(ones_fab))
       real(amrex_real), INTENT(in), target :: type_fab(DIMV(type_fab))
       real(amrex_real), INTENT(in), target :: color_fab(DIMV(color_fab))
       real(amrex_real), INTENT(in), target :: alpha_fab(DIMV(alpha_fab),nsolve)
       real(amrex_real), INTENT(in), target :: mask_fab(DIMV(mask_fab))
       integer, INTENT(in) :: presbc(SDIM,2)
       integer, INTENT(in) :: type_flag(2)

       integer nc
       integer i,j,k
       integer local_mask,local_color,local_type,local_ones
       integer plus_flag,zero_flag
       real(amrex_real) local_alpha
       integer dir_local
       integer icrit
       integer side
       integer local_bc

       data_fab_ptr=>data_fab

       if ((levelno.lt.0).or.(gridno.lt.0)) then
        print *,"level or grid invalid"
        stop
       endif
       if ((nsolve.ne.1).and.(nsolve.ne.SDIM)) then
        print *,"nsolve invalid"
        stop
       endif
       if (project_option.eq.SOLVETYPE_VISC) then
        if (nsolve.eq.SDIM) then
         ! do nothing
        else
         print *,"nsolve invalid"
         stop
        endif
       else if (project_option_is_validF(project_option).eq.1) then
        if (nsolve.eq.1) then
         ! do nothing
        else
         print *,"nsolve invalid"
         stop
        endif
       else
        print *,"project_option invalid; fort_fabcom_ones"
        stop
       endif
       if (color_count.ge.1) then
        ! do nothing
       else
        print *,"color_count invalid"
        stop
       endif

       call checkbound_array1(fablo,fabhi,data_fab_ptr,0,-1) 
       call checkbound_array1(fablo,fabhi,ones_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,type_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,color_fab,0,-1) 
       call checkbound_array(fablo,fabhi,alpha_fab,0,-1) 
       call checkbound_array1(fablo,fabhi,mask_fab,0,-1) 

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        local_mask=NINT(mask_fab(D_DECL(i,j,k)))

        if (local_mask.eq.1) then

         local_color=NINT(color_fab(D_DECL(i,j,k)))

         if ((local_color.ge.1).and.(local_color.le.color_count)) then

          if (project_option_singular_possibleF(project_option).eq.1) then

            ! prescribed solid in patch.
           if (singular_patch_flag(local_color).eq.0) then
            ! do nothing (masked off region)

            !no comp or dirichlet bc in patch region.
           else if (singular_patch_flag(local_color).eq.1) then
            data_fab(D_DECL(i,j,k))= &
              data_fab(D_DECL(i,j,k))+beta(local_color)

            !comp or dirichlet bc in patch region.
           else if (singular_patch_flag(local_color).eq.2) then
            ! do nothing
           else
            print *,"singular_patch_flag invalid"
            stop
           endif

          else if (project_option_singular_possibleF(project_option).eq.0) then
           ! do nothing
          else
           print *,"project_option_singular_possible invalid; fort_fabcom_ones"
           stop
          endif

          local_type=NINT(type_fab(D_DECL(i,j,k)))
          local_ones=NINT(ones_fab(D_DECL(i,j,k)))
           ! completely masked off cell
          if ((local_ones.eq.0).and.(local_type.eq.1)) then
           if (type_flag(local_type).eq.1) then
            ! do nothing (type exists)
           else
            print *,"type_flag invalid"
            stop
           endif
          else if ((local_ones.eq.1).and.(local_type.eq.2)) then
           if ((singular_patch_flag(local_color).eq.1).or. &
               (singular_patch_flag(local_color).eq.2)) then
            ! do nothing
           else
            print *,"singular_patch_flag invalid"
            stop
           endif
           if (type_flag(local_type).eq.1) then
            ! do nothing (type exists)
           else
            print *,"type_flag invalid"
            stop
           endif
           plus_flag=0
           zero_flag=0
           do nc=1,nsolve
            local_alpha=alpha_fab(D_DECL(i,j,k),nc)
            if (local_alpha.eq.zero) then
             zero_flag=1
            else if (local_alpha.gt.zero) then
             plus_flag=1
            else
             print *,"local_alpha invalid"
             stop
            endif 
           enddo !nc=1..nsolve
           if ((plus_flag.eq.1).and.(zero_flag.eq.0)) then
            if (singular_patch_flag(local_color).eq.2) then
             ! do nothing
            else
             print *,"singular_patch_flag invalid"
             stop
            endif
           else if ((plus_flag.eq.0).and.(zero_flag.eq.1)) then
            ! do nothing
           else
            print *,"plus_flag or zero_flag invalid"
            stop
           endif
           do dir_local=1,SDIM
            if (dir_local.eq.1) then
             icrit=i
            else if (dir_local.eq.2) then
             icrit=j
            else if ((dir_local.eq.3).and.(SDIM.eq.3)) then
             icrit=k
            else
             print *,"dir_local invalid"
             stop
            endif
            side=0
            if ((icrit.gt.fablo(dir_local)).and. &
                (icrit.lt.fabhi(dir_local))) then
             ! do nothing
            else if (icrit.eq.fablo(dir_local)) then
             side=1
            else if (icrit.eq.fabhi(dir_local)) then
             side=2
            else
             print *,"icrit invalid"
             stop
            endif
            if (side.eq.0) then
             ! do nothing
            else if ((side.eq.1).or.(side.eq.2)) then
             local_bc=presbc(dir_local,side)
             if (project_option.eq.SOLVETYPE_PRESEXTRAP) then 
              ! do nothing (all bcs are Neumann)
             else if &
               (project_option_singular_possibleF(project_option).eq.1) then
              if (local_bc.eq.INT_DIR) then
               ! do nothing
              else if (local_bc.eq.FOEXTRAP) then
               ! do nothing
              else if (local_bc.eq.REFLECT_EVEN) then
               ! do nothing
              else if (local_bc.eq.EXT_DIR) then
               if (singular_patch_flag(local_color).eq.2) then
                ! do nothing
               else
                print *,"singular_patch_flag invalid"
                stop
               endif
              else
               print *,"local_bc invalid"
               stop
              endif
             else if &
               (project_option_singular_possibleF(project_option).eq.0) then
              ! do nothing
             else
              print *,"project_option invalid; fort_fabcom_ones"
              stop
             endif
            else
             print *,"side invalid"
             stop
            endif
           enddo !dir_local=1..sdim
          else
           print *,"local_ones or local_type invalid"
           stop
          endif

         else
          print *,"local_color invalid"
          stop
         endif

        else if (local_mask.eq.0) then
         ! do nothing
        else 
         print *,"mask invalid"
         stop
        endif

       enddo ! k
       enddo ! j
       enddo ! i

       return
       end subroutine fort_fabcom_ones

! coriolis force:
! R''_space=R''_earth + 2 (omega cross v)+omega cross (omega cross R)+
! w' cross R.
! coriolis term=-2 Omega cross v
! suppose Omega=(0 0 omega)  (z axis)
! coriolis force=-2 ( -omega v  omega u  0 ) 
! centrifugal force:
! F=mass r w^2   w units: radians/s   mass r w^2 = kg m/s^2
!
      subroutine fort_init_potential( &
       presden,DIMS(presden), &
       state,DIMS(state), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       presbc_arr, &
       dombcpres, &
       domlo,domhi, &
       xlo,dx, &
       dt, & ! dt_slab
       angular_velocity, & !intent(in) fort_init_potential
       centrifugal_force_factor, & !intent(in) fort_init_potential
       isweep) &
      bind(c,name='fort_init_potential')

      use global_utility_module
      use probf90_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: isweep
      real(amrex_real), INTENT(in) :: dt  ! dt_slab
      real(amrex_real), INTENT(in) :: angular_velocity
      real(amrex_real), INTENT(in) :: centrifugal_force_factor
      integer, INTENT(in) :: DIMDEC(presden)
      integer, INTENT(in) :: DIMDEC(state)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: domlo(SDIM),domhi(SDIM)
       ! pressure bc at exterior walls
      integer, INTENT(in) :: dombcpres(SDIM,2) 
       ! presbc at int or ext.
      integer, INTENT(in) :: presbc_arr(SDIM,2) 

       !HYDROSTATIC_PRESSURE,HYDROSTATIC_DENSITY
      real(amrex_real), INTENT(inout),target :: presden(DIMV(presden),2) 
      real(amrex_real), pointer :: presden_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
          state(DIMV(state),num_materials*num_state_material) 
      real(amrex_real), pointer :: state_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      real(amrex_real) den_cell,pres_cell

      integer i,j,k
      integer ii,jj,kk,dir,side,inside
      integer bctypepres,local_bctype,exteriorbc

      if (bfact.lt.1) then
       print *,"bfact invalid159"
       stop
      endif

      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif

      presden_ptr=>presden
      state_ptr=>state
      call checkbound_array(fablo,fabhi,presden_ptr,1,-1)
      call checkbound_array(fablo,fabhi,state_ptr,1,-1)
     
      if (isweep.eq.0) then
 
       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0)
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

         ! e.g. for gravity force:
         ! grad p/rho = \vec{g}
         ! grad p = rho\vec{g}
         ! p = dt rho \vec{g} \cdot \vec{x} 
         ! Du/Dt=-grad p/rho - omega cross (omega cross r)- 2 omega cross u +
         !       \vec{g} 
         ! omega z^hat cross r r^hat=omega r theta^hat
         ! omega z^hat cross (omega z^hat cross r)=
         ! omega z^hat cross omega r theta^hat = -omega^2 r r^hat
         !
         ! includes centrifugal force but not "coriolis force"
         !  den_cell=fort_denconst(1)
         !  p=dt( \vec{g} \cdot \vec{x} + (1/2)Omega^2 r^2 )
         ! general_hydrostatic_pressure_density is declared in:
         !  PROB.F90
        call general_hydrostatic_pressure_density( &
          i,j,k,level, &
          angular_velocity, &
          centrifugal_force_factor, &
          dt, &
          den_cell, &   !INTENT(out)
          pres_cell, &  !INTENT(out)
          state_ptr)
        presden(D_DECL(i,j,k),1)=pres_cell
        presden(D_DECL(i,j,k),2)=den_cell
       enddo     
       enddo     
       enddo     

      else if (isweep.eq.1) then

       do dir=1,SDIM
       do side=1,2

         ! ngrow=1
        call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1)

         ! dombcpres=get_desc_lst()[State_Type].getBC(STATECOMP_PRES)
        bctypepres=dombcpres(dir,side)  
         ! presbc_arr=getBCArray(State_Type,gridno,STATECOMP_PRES,1)
        local_bctype=presbc_arr(dir,side) 

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

         do k=growlo(3),growhi(3)
         do j=growlo(2),growhi(2)
         do i=growlo(1),growhi(1)

          if (local_bctype.eq.bctypepres) then

            ! outflow wall
           if (local_bctype.eq.EXT_DIR) then

             ! p=dt( \vec{g}\cdot\vec{x} + (1/2)Omega^2 r^2 )
            call general_hydrostatic_pressure_density( &
             i,j,k,level, &
             angular_velocity, & !intent(in)
             centrifugal_force_factor, & !intent(in)
             dt, &        ! intent(in)
             den_cell, &  ! intent(out)
             pres_cell, & ! intent(out)
             state_ptr)

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
             print *,"dir invalid fort_init_potential"
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
      end subroutine fort_init_potential

! called from:
! NavierStokes3.cpp: NavierStokes::increment_potential_force()
      subroutine fort_addgravity( &
       dt, &
       cur_time, &
       level, &
       finest_level, &
       nstate, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       xlo,dx,dir, &
       xface,DIMS(xface), &
       lsnew,DIMS(lsnew), &
       macnew,DIMS(macnew), &
       facegrav,DIMS(facegrav) ) &
      bind(c,name='fort_addgravity')

      use global_utility_module
      use probcommon_module
      use probf90_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(in) :: cur_time
      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nstate
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: DIMDEC(xface)
      integer, INTENT(in) :: DIMDEC(lsnew)
      integer, INTENT(in) :: DIMDEC(macnew)
      integer, INTENT(in) :: DIMDEC(facegrav)

      real(amrex_real), INTENT(in),target :: xface(DIMV(xface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: xface_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in),target :: &
              lsnew(DIMV(lsnew),num_materials*(SDIM+1))
      real(amrex_real), pointer :: lsnew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: macnew(DIMV(macnew))
      real(amrex_real), pointer :: macnew_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: facegrav(DIMV(facegrav))
      real(amrex_real), pointer :: facegrav_ptr(D_DECL(:,:,:))
 
      integer i,j,k
      integer ii,jj,kk

      real(amrex_real) local_cut
      real(amrex_real) local_macnew

      real(amrex_real) gravity_increment

      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      xface_ptr=>xface
      lsnew_ptr=>lsnew
      macnew_ptr=>macnew
      facegrav_ptr=>facegrav

      if (bfact.lt.1) then
       print *,"bfact invalid160"
       stop
      endif
      if ((level.lt.0).or.(level.gt.finest_level)) then
       print *,"level invalid add gravity"
       stop
      endif
 
      if (nstate.ne.STATE_NCOMP) then
       print *,"nstate invalid"
       stop
      endif
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid"
       stop
      endif
      if (dt.gt.zero) then
       ! do nothing
      else
       print *,"dt invalid in addgravity"
       stop
      endif
      if (cur_time.ge.zero) then
       ! do nothing
      else
       print *,"cur_time invalid in addgravity"
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

      call checkbound_array(fablo,fabhi,xface_ptr,0,dir)
      call checkbound_array(fablo,fabhi,lsnew_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,macnew_ptr,0,dir)
      call checkbound_array1(fablo,fabhi,facegrav_ptr,0,dir)

      call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir)

       local_cut=xface(D_DECL(i,j,k),FACECOMP_FACECUT+1)
       if ((local_cut.ge.zero).and.(local_cut.le.half)) then
        local_cut=zero
       else if ((local_cut.ge.half).and.(local_cut.le.one)) then
        ! do nothing
       else
        print *,"local_cut invalid"
        stop
       endif

        ! 1. surface tension 
        ! 2. gravity 
       gravity_increment=local_cut*facegrav(D_DECL(i,j,k))

       local_macnew=macnew(D_DECL(i,j,k))+gravity_increment

       if (levelrz.eq.COORDSYS_CARTESIAN) then
        ! do nothing
       else if (levelrz.eq.COORDSYS_RZ) then
        if (SDIM.ne.2) then
         print *,"rz in 2d only"
         stop
        endif
        if ((dir.eq.0).and. &
            (xsten(0,1).le.EPS2*dx(1))) then
         local_macnew=zero
        endif
       else if (levelrz.eq.COORDSYS_CYLINDRICAL) then
        ! do nothing
       else
        print *,"levelrz invalid"
        stop 
       endif

       macnew(D_DECL(i,j,k))=local_macnew

      enddo
      enddo
      enddo

      return
      end subroutine fort_addgravity

      subroutine fort_eos_pressure( &
        level, &
        finest_level, &
        local_material_type, &
        xlo,dx, &
        pres,DIMS(pres), &
        levelpc,DIMS(levelpc), &
        den,DIMS(den), &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        nden) &
      bind(c,name='fort_eos_pressure')

      use global_utility_module
      use probf90_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: nden
      integer, INTENT(in) :: local_material_type(num_materials)
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
    
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)

      integer, INTENT(in) :: DIMDEC(pres)
      integer, INTENT(in) :: DIMDEC(levelpc)
      integer, INTENT(in) :: DIMDEC(den)
      real(amrex_real), INTENT(inout), target :: pres(DIMV(pres))
      real(amrex_real), pointer :: pres_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: &
              levelpc(DIMV(levelpc),num_materials*(1+SDIM))
      real(amrex_real), pointer :: levelpc_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: den(DIMV(den),nden) ! den,temp,Y
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:),:)

      integer i,j,k
      integer im
      integer im_primary
      integer ibase
      real(amrex_real) rho,internal_energy,TEMP
      real(amrex_real) LS(num_materials)
      real(amrex_real) massfrac_parm(num_species_var+1)
      integer ispec

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
      if (nden.ne.num_materials*num_state_material) then
       print *,"nden invalid"
       stop
      endif

      pres_ptr=>pres
      levelpc_ptr=>levelpc
      den_ptr=>den

      call checkbound_array(fablo,fabhi,levelpc_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,pres_ptr,1,-1)
      call checkbound_array(fablo,fabhi,den_ptr,1,-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1)
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       do im=1,num_materials
        LS(im)=levelpc(D_DECL(i,j,k),im)
       enddo ! im

       call get_primary_material(LS,im_primary)

       ibase=(im_primary-1)*num_state_material

       if (is_rigid(im_primary).eq.0) then

        ! compressible material
        if ((local_material_type(im_primary).gt.0).and. &
            (local_material_type(im_primary).le.MAX_NUM_EOS)) then

         rho=den(D_DECL(i,j,k),ibase+ENUM_DENVAR+1)
         if (rho.gt.zero) then
          ! do nothing
         else
          print *,"density has gone nonpos"
          stop
         endif
         TEMP=den(D_DECL(i,j,k),ibase+ENUM_TEMPERATUREVAR+1)
         if (TEMP.gt.zero) then
          ! do nothing
         else
          print *,"TEMP has gone nonpos"
          stop
         endif
         call init_massfrac_parm(rho,massfrac_parm,im_primary)
         do ispec=1,num_species_var
          massfrac_parm(ispec)=den(D_DECL(i,j,k),ibase+ENUM_SPECIESVAR+ispec)
         enddo
         ! returns energy/scale
         call INTERNAL_material(rho,massfrac_parm,TEMP, &
          internal_energy,local_material_type(im_primary),im_primary)
         if (internal_energy.gt.zero) then
          ! do nothing
         else
          print *,"internal_energy has gone nonpos"
          stop
         endif
         ! p(energy*scale)/scale
         call EOS_material(rho,massfrac_parm, &
          internal_energy, &
          pres(D_DECL(i,j,k)), &
          local_material_type(im_primary),im_primary)
        else if (local_material_type(im_primary).eq.0) then
         ! do nothing
        else
         print *,"fort material type invalid"
         stop
        endif

       else if (is_rigid(im_primary).eq.1) then
        ! do nothing
       else
        print *,"is_rigid bust"
        stop
       endif 

      enddo
      enddo
      enddo  ! i,j,k

      return
      end subroutine fort_eos_pressure


! 2 files:
!   height_integral.txt
!   jet_height.txt
      subroutine fort_coflow( &
       time, &
       fdomlo, &
       fdomhi, &
       Z_dir, & !0..sdim-1
       R_dir, & !0..sdim-1 
       num_cells, &
       coflow_Z, &
       coflow_R_of_Z) &
      bind(c,name='fort_coflow')

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: fdomlo(SDIM)
      integer, INTENT(in) :: fdomhi(SDIM)
      integer, INTENT(inout) :: Z_dir !0..sdim-1
      integer, INTENT(inout) :: R_dir !0..sdim-1
      integer, INTENT(inout) :: num_cells
      real(amrex_real), INTENT(inout) :: coflow_Z(0:num_cells)
      real(amrex_real), INTENT(inout) :: coflow_R_of_Z(0:num_cells)
      real(amrex_real) R_of_Z_MASK(0:num_cells)
      real(amrex_real) major,minor,wavelength,amplitude,minrad
      real(amrex_real) masklo,maskhi,HH,SS
      integer j,jleft,jright
      real(amrex_real) dzz,theta
      integer local_num_cells
      integer local_Z_dir
      integer local_R_dir
      integer dir
      integer js,je

      do dir=1,SDIM
       if (fdomlo(dir).eq.0) then
        ! do nothing
       else
        print *,"fdomlo invalid"
        stop
       endif
      enddo !dir=1,sdim

      local_num_cells=0
      local_Z_dir=1
      local_R_dir=0
    
      if ((SDIM.eq.2).and. &
          (probtype.eq.41).and. &
          (axis_dir.eq.4)) then
       local_Z_dir=SDIM-1
       local_R_dir=0
       local_num_cells=fdomhi(local_Z_dir+1)-fdomlo(local_Z_dir+1)+1
      else if ((SDIM.eq.3).and. &
               (probtype.eq.53).and. &
               (axis_dir.eq.0)) then
       local_Z_dir=0
       local_R_dir=SDIM-1
       local_num_cells=fdomhi(local_Z_dir+1)-fdomlo(local_Z_dir+1)+1
      else
       local_num_cells=0
      endif
 
      if (num_cells.eq.0) then ! setup
       num_cells=local_num_cells
       Z_dir=local_Z_dir
       R_dir=local_R_dir
      else if (num_cells.gt.0) then ! diagnostics
       Z_dir=local_Z_dir
       R_dir=local_R_dir
       js=fdomlo(local_Z_dir+1)
       je=fdomhi(local_Z_dir+1)

       if (num_cells.eq.local_num_cells) then
        ! do nothing
       else 
        print *,"num_cells invalid"
        stop
       endif
       if ((Z_dir.lt.0).or.(Z_dir.ge.SDIM)) then
        print *,"Z_dir invalid"
        stop
       endif
       if ((R_dir.lt.0).or.(R_dir.ge.SDIM).or.(R_dir.eq.Z_dir)) then
        print *,"R_dir or Z_dir invalid"
        stop
       endif

       if (SDIM.eq.2) then
        major=yblob
        if (abs(major-coflow_Z(num_cells)).le.EPS10) then
         call FINDAMPLITUDE_2D(num_cells,coflow_Z,coflow_R_of_Z,major,minor)
         call FINDMINRAD(num_cells,coflow_Z,coflow_R_of_Z,minrad)
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
       else if (SDIM.eq.3) then
        open(unit=19,file="height_integral.txt",access="append")
        open(unit=20,file="jet_height.txt",access="append")
        maskhi=one
        masklo=half*(xblob+maskhi)

        coflow_Z(num_cells)=probhix-problox
        dzz=(probhix-problox)/num_cells

         ! fill in gaps in coarse regions
        do j=js+1,je
         if (coflow_Z(j).eq.zero) then  ! value not initialized in summass
          coflow_Z(j)=j*dzz
          jleft=j-1
          jright=j+1
          do while (coflow_Z(jright).eq.zero)
           jright=jright+1
          enddo  
          theta=one/(jright-jleft)
          coflow_R_of_Z(j)= &
           coflow_R_of_Z(jleft)*(one-theta)+theta*coflow_R_of_Z(jright)
         endif
        enddo

        do j=js,je+1 
         if (coflow_Z(j).lt.masklo) then
          R_of_Z_MASK(j)=zero
         else if (coflow_Z(j).gt.maskhi) then
          R_of_Z_MASK(j)=zero
         else 
          R_of_Z_MASK(j)=coflow_R_of_Z(j)
         endif
        enddo
        HH=coflow_Z(num_cells)-coflow_Z(0)
        call SIMP(num_cells,HH,R_of_Z_MASK,SS)
        if ((SS.lt.zero).or.(HH.le.zero)) then
         print *,"SS or HH invalid"
         print *,"SS= ",SS
         print *,"HH= ",HH
         print *,"num_cells= ",num_cells
         do j=js,je+1
          print *,"j,coflow_Z ",j,coflow_Z(j)
         enddo
         stop
        endif
        SS=SS/HH
        write(19,*) time,SS 
        print *,"TIME,AVERAGE HEIGHT ",time,SS
        do j=js,je+1 
         write(20,*) coflow_Z(j),coflow_R_of_Z(j)
        enddo
        write(20,*) ' '
 
        close(19) 
        close(20) 
       else 
        print *,"dimension bust"
        stop
       endif
      else
       print *,"num_cells invalid"
       stop
      endif

      return
      end subroutine fort_coflow

       ! do_the_advance -> level_phase_change_rate
       ! avgDownBURNING_localMF 
       ! level_avgDownBURNING
       ! (note: after level_avgDownBURNING comes 
       !  level_phase_change_rate_extend)
      subroutine fort_avgdown_burning( &
       velflag, &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) &
      bind(c,name='fort_avgdown_burning')

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: velflag
      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level_f
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact_f
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: xlo_fine(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      integer, INTENT(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real), target, INTENT(out) :: crse(DIMV(crse),ncomp)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: fine(DIMV(fine),ncomp)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer dir2
      integer iten
      integer n
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      real(amrex_real) crse_value(ncomp)
      integer fine_test
      integer coarse_test
      real(amrex_real) velwt(num_interfaces)
      integer local_comp
      integer avgdown_sweep
      integer ncomp_expect
      integer ncomp_per_interface

      if (velflag.eq.1) then  ! burning velocity
       ncomp_per_interface=EXTRAP_PER_BURNING
      else if (velflag.eq.0) then ! interface temperature, mass fraction
       ncomp_per_interface=EXTRAP_PER_TSAT
      else
       print *,"velflag invalid"
       stop
      endif
      ncomp_expect=num_interfaces+num_interfaces*ncomp_per_interface

      if (ncomp.eq.ncomp_expect) then
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

      crse_ptr=>crse
      call checkbound_array(lo,hi,crse_ptr,0,-1)
      fine_ptr=>fine
      call checkbound_array(lof,hif,fine_ptr,0,-1)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

       do n=1,ncomp
        crse_value(n)=zero
       enddo

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
        bfact_c,bfact_f)

       do avgdown_sweep=0,1

        do iten=1,num_interfaces
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

              do iten=1,num_interfaces
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
                  do dir2=1,ncomp_per_interface
                   local_comp=num_interfaces+(iten-1)*ncomp_per_interface+dir2
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
              enddo ! iten=1..num_interfaces

              voltotal=voltotal+volall
             endif ! wt(sdim).gt.0
            enddo ! kfine
           endif
          enddo ! jfine
         endif
        enddo ! ifine

        if (voltotal.gt.zero) then
         !do nothing
        else
         print *,"voltotal invalid: ",voltotal
         stop
        endif

        if (avgdown_sweep.eq.0) then
         ! do nothing
        else if (avgdown_sweep.eq.1) then

         do iten=1,num_interfaces

          coarse_test=NINT(crse_value(iten))
          if (coarse_test.eq.0) then
           if (velwt(iten).eq.zero) then
            ! do nothing
           else
            print *,"velwt invalid"
            stop
           endif
           crse(D_DECL(ic,jc,kc),iten)=zero
           do dir2=1,ncomp_per_interface
            crse(D_DECL(ic,jc,kc),num_interfaces+ &
                 (iten-1)*ncomp_per_interface+dir2)=zero
           enddo
          else if ((coarse_test.eq.1).or. &
                   (coarse_test.eq.-1)) then
           if (velwt(iten).gt.zero) then
            crse(D_DECL(ic,jc,kc),iten)=coarse_test
            do dir2=1,ncomp_per_interface
             crse(D_DECL(ic,jc,kc),num_interfaces+ &
                  (iten-1)*ncomp_per_interface+dir2)= &
              crse_value(num_interfaces+ &
                         (iten-1)*ncomp_per_interface+dir2)/velwt(iten)
            enddo
           else
            print *,"velwt invalid"
            stop
           endif
          else
           print *,"coarse_test invalid"
           stop
          endif

         enddo ! iten=1..num_interfaces

        else 
         print *,"avgdown_sweep invalid"
         stop
        endif

       enddo ! avgdown_sweep=0..1

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine fort_avgdown_burning

       ! 1. GetDragALL
       ! 2. GetDrag
       ! 3. avgDownDRAG_MF 
       ! 4. level_avgDownDRAG
       ! 5. level_DRAG_extend
      subroutine fort_avgdown_drag( &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) &
      bind(c,name='fort_avgdown_drag')

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level_f
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact_f
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: xlo_fine(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      integer, INTENT(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real), target, INTENT(out) :: crse(DIMV(crse),ncomp)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), target, INTENT(in) :: fine(DIMV(fine),ncomp)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer dir2
      integer n
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      real(amrex_real) crse_value(ncomp)
      integer fine_test
      integer coarse_test
      real(amrex_real) velwt(num_materials)
      integer local_comp
      integer avgdown_sweep
      integer ncomp_expect
      integer im_test
      integer drag_type,drag_im

      ncomp_expect=N_DRAG

      if (ncomp.eq.ncomp_expect) then
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

      crse_ptr=>crse
      call checkbound_array(lo,hi,crse_ptr,0,-1)
      fine_ptr=>fine
      call checkbound_array(lof,hif,fine_ptr,0,-1)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

       do n=1,ncomp
        crse_value(n)=zero
       enddo

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
        bfact_c,bfact_f)

       do avgdown_sweep=0,1

        do im_test=1,num_materials
         velwt(im_test)=zero
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

              do im_test=1,num_materials
               fine_test=NINT(fine(D_DECL(ifine,jfine,kfine), &
                       DRAGCOMP_FLAG+im_test))
                ! crse_value init to 0
               coarse_test=NINT(crse_value(DRAGCOMP_FLAG+im_test)) 

               if (avgdown_sweep.eq.0) then

                if (coarse_test.eq.0) then 
                 coarse_test=fine_test
                else if (fine_test.eq.0) then
                 ! do nothing
                else if (fine_test.eq.1) then
                 coarse_test=fine_test
                else
                 print *,"fine_test invalid"
                 stop
                endif

                if (coarse_test.eq.1) then
                 crse_value(DRAGCOMP_FLAG+im_test)=coarse_test
                else if (coarse_test.eq.0) then
                 ! do nothing (crse_value init. to zero)
                else
                 print *,"coarse_test invalid"
                 stop
                endif

               else if (avgdown_sweep.eq.1) then

                if ((coarse_test.eq.0).or. &
                    (coarse_test.eq.1)) then

                 if (fine_test.eq.0) then
                  ! do nothing
                 else if (fine_test.eq.coarse_test) then
                  velwt(im_test)=velwt(im_test)+volall
                  do local_comp=0,N_DRAG-1
                   drag_type=fort_drag_type(local_comp,drag_im)
                   if ((drag_type.ge.0).and. &
                       (drag_type.lt.DRAG_TYPE_NEXT)) then

                    if ((drag_im.ge.0).and. &
                        (drag_im.lt.num_materials)) then
                     if (drag_im.eq.im_test-1) then
                      if (drag_type.ne.DRAG_TYPE_FLAG) then
                       crse_value(local_comp+1)=crse_value(local_comp+1)+ &
                        volall*fine(D_DECL(ifine,jfine,kfine),local_comp+1)
                      else if (drag_type.eq.DRAG_TYPE_FLAG) then
                       ! do nothing
                      else
                       print *,"drag_type invalid"
                       stop
                      endif
                     else if (drag_im.ne.im_test-1) then
                      ! do nothing
                     else
                      print *,"drag_im invalid"
                      stop
                     endif
                    else
                     print *,"drag_im invalid"
                     stop
                    endif
                   else
                    print *,"drag_type invalid"
                    stop
                   endif
                  enddo !local_comp=0,N_DRAG-1

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
              enddo !do im_test=1,num_materials

              voltotal=voltotal+volall
             endif ! wt(sdim).gt.0
            enddo ! kfine
           endif
          enddo ! jfine
         endif
        enddo ! ifine

        if (voltotal.gt.zero) then
         ! do nothing
        else
         print *,"voltotal invalid: ",voltotal
         stop
        endif

        if (avgdown_sweep.eq.0) then
         ! do nothing
        else if (avgdown_sweep.eq.1) then

         do im_test=1,num_materials

          coarse_test=NINT(crse_value(DRAGCOMP_FLAG+im_test)) 
          if (coarse_test.eq.0) then
           if (velwt(im_test).eq.zero) then
            ! do nothing
           else
            print *,"velwt invalid"
            stop
           endif
           crse(D_DECL(ic,jc,kc),DRAGCOMP_FLAG+im_test)=coarse_test
           do local_comp=0,N_DRAG-1
            drag_type=fort_drag_type(local_comp,drag_im)
            if ((drag_type.ge.0).and. &
                (drag_type.lt.DRAG_TYPE_NEXT)) then
             if ((drag_im.ge.0).and. &
                 (drag_im.lt.num_materials)) then
              if (drag_im.eq.im_test-1) then
               if (drag_type.ne.DRAG_TYPE_FLAG) then
                crse(D_DECL(ic,jc,kc),local_comp+1)=zero
               else if (drag_type.eq.DRAG_TYPE_FLAG) then
                ! do nothing
               else
                print *,"drag_type invalid"
                stop
               endif
              else if (drag_im.ne.im_test-1) then
               ! do nothing
              else
               print *,"drag_im invalid"
               stop
              endif
             else
              print *,"drag_im invalid"
              stop
             endif
            else
             print *,"drag_type invalid"
             stop
            endif
           enddo !local_comp=0,N_DRAG-1
          else if (coarse_test.eq.1) then
           if (velwt(im_test).gt.zero) then
            crse(D_DECL(ic,jc,kc),DRAGCOMP_FLAG+im_test)=coarse_test

            do local_comp=0,N_DRAG-1
             drag_type=fort_drag_type(local_comp,drag_im)
             if ((drag_type.ge.0).and. &
                 (drag_type.lt.DRAG_TYPE_NEXT)) then
              if ((drag_im.ge.0).and. &
                  (drag_im.lt.num_materials)) then
               if (drag_im.eq.im_test-1) then
                if (drag_type.ne.DRAG_TYPE_FLAG) then
                 crse(D_DECL(ic,jc,kc),local_comp+1)= &
                   crse_value(local_comp+1)/velwt(im_test)
                else if (drag_type.eq.DRAG_TYPE_FLAG) then
                 ! do nothing
                else
                 print *,"drag_type invalid"
                 stop
                endif
               else if (drag_im.ne.im_test-1) then
                ! do nothing
               else
                print *,"drag_im invalid"
                stop
               endif
              else
               print *,"drag_im invalid"
               stop
              endif
             else
              print *,"drag_type invalid"
              stop
             endif
            enddo !local_comp=0,N_DRAG-1

           else
            print *,"velwt invalid"
            stop
           endif
          else
           print *,"coarse_test invalid"
           stop
          endif

         enddo !do im_test=1,num_materials

        else 
         print *,"avgdown_sweep invalid"
         stop
        endif

       enddo ! avgdown_sweep=0..1

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine fort_avgdown_drag


       ! updatevelALL called from: NavierStokes::multiphase_project
       ! mac_update is called from: NavierStokes::updatevelALL
       ! correct_velocity called from: NavierStokes::multiphase_project
       ! correct_velocity called from: NavierStokes::relaxLEVEL
       ! correct_velocity called from: NavierStokes::residual_correction_form
       ! correct_velocity called from: NavierStokes::mac_update
       ! called from: NavierStokes::correct_velocity (NavierStokes3.cpp)
      subroutine fort_fluidsolidcor( &
       level, &
       finest_level, &
       velcomp, &
       nsolve, &
       project_option, &
       tilelo,tilehi, &
       fablo,fabhi,bfact, &
       presbc_in, &
       maskcov,DIMS(maskcov), &
       xface,DIMS(xface), &
       yface,DIMS(yface), &
       zface,DIMS(zface), &
       xgp,DIMS(xgp), &
       ygp,DIMS(ygp), &
       zgp,DIMS(zgp), &
       xsrc,DIMS(xsrc), &
       ysrc,DIMS(ysrc), &
       zsrc,DIMS(zsrc), &
       xdest,DIMS(xdest), &
       ydest,DIMS(ydest), &
       zdest,DIMS(zdest), &
       xlo, &
       dx, &
       dt_slab, &
       cur_time) &
      bind(c,name='fort_fluidsolidcor')

      use global_utility_module
      use probcommon_module
      use probf90_module

      IMPLICIT NONE

      integer, INTENT(in) :: level,finest_level
      integer, INTENT(in) :: velcomp
      integer :: veldir
      integer, INTENT(in) :: nsolve
      integer, INTENT(in) :: project_option
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(xface)
      integer, INTENT(in) :: DIMDEC(yface)
      integer, INTENT(in) :: DIMDEC(zface)
      integer, INTENT(in) :: DIMDEC(xgp)
      integer, INTENT(in) :: DIMDEC(ygp)
      integer, INTENT(in) :: DIMDEC(zgp)
      integer, INTENT(in) :: DIMDEC(xsrc)
      integer, INTENT(in) :: DIMDEC(ysrc)
      integer, INTENT(in) :: DIMDEC(zsrc)
      integer, INTENT(in) :: DIMDEC(xdest)
      integer, INTENT(in) :: DIMDEC(ydest)
      integer, INTENT(in) :: DIMDEC(zdest)

      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer :: dir
      integer, INTENT(in) :: presbc_in(SDIM,2,nsolve)

      real(amrex_real), INTENT(in), target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: xface(DIMV(xface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: xface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: yface(DIMV(yface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: yface_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in), target :: zface(DIMV(zface),FACECOMP_NCOMP)
      real(amrex_real), pointer :: zface_ptr(D_DECL(:,:,:),:)

      real(amrex_real), INTENT(in), target :: xgp(DIMV(xgp))
      real(amrex_real), pointer :: xgp_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: ygp(DIMV(ygp))
      real(amrex_real), pointer :: ygp_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: zgp(DIMV(zgp))
      real(amrex_real), pointer :: zgp_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in), target :: xsrc(DIMV(xsrc))
      real(amrex_real), pointer :: xsrc_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: ysrc(DIMV(ysrc))
      real(amrex_real), pointer :: ysrc_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in), target :: zsrc(DIMV(zsrc))
      real(amrex_real), pointer :: zsrc_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(out), target :: xdest(DIMV(xdest))
      real(amrex_real), pointer :: xdest_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out), target :: ydest(DIMV(ydest))
      real(amrex_real), pointer :: ydest_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out), target :: zdest(DIMV(zdest))
      real(amrex_real), pointer :: zdest_ptr(D_DECL(:,:,:))

      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      real(amrex_real), INTENT(in) :: dt_slab
      real(amrex_real), INTENT(in) :: cur_time
 
      integer i,j,k
      integer ii,jj,kk
      real(amrex_real) local_dest,local_src,local_gp
      real(amrex_real) AL
      real(amrex_real) AL_ice
      real(amrex_real) AL_ice_mask
      real(amrex_real) local_dd,local_visc_coef
      real(amrex_real) cc_group,local_dd_group
      integer local_uncoupled_viscosity
      integer icrit,side,bccrit
      integer bc_comp
      real(amrex_real) local_wt(nsolve)
      real(amrex_real) maskleft,maskright
      integer maskface
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      if ((level.ge.0).and.(level.le.finest_level)) then
       ! do nothing
      else
       print *,"level invalid"
       stop
      endif

      if (project_option_is_validF(project_option).eq.1) then
       ! do nothing
      else
       print *,"project_option invalid fluid solid cor"
       stop
      endif 
      if (bfact.lt.1) then
       print *,"bfact invalid161"
       stop
      endif

      if (dt_slab.gt.zero) then
       ! do nothing
      else
       print *,"solve_dt_slab invalid"
       stop
      endif 
      maskcov_ptr=>maskcov
      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)

      xface_ptr=>xface
      yface_ptr=>yface
      zface_ptr=>zface
      call checkbound_array(fablo,fabhi,xface_ptr,0,0)
      call checkbound_array(fablo,fabhi,yface_ptr,0,1)
      call checkbound_array(fablo,fabhi,zface_ptr,0,SDIM-1)

      xgp_ptr=>xgp
      ygp_ptr=>ygp
      zgp_ptr=>zgp
      call checkbound_array1(fablo,fabhi,xgp_ptr,0,0)
      call checkbound_array1(fablo,fabhi,ygp_ptr,0,1)
      call checkbound_array1(fablo,fabhi,zgp_ptr,0,SDIM-1)

      xsrc_ptr=>xsrc
      ysrc_ptr=>ysrc
      zsrc_ptr=>zsrc
      call checkbound_array1(fablo,fabhi,xsrc_ptr,0,0)
      call checkbound_array1(fablo,fabhi,ysrc_ptr,0,1)
      call checkbound_array1(fablo,fabhi,zsrc_ptr,0,SDIM-1)

      xdest_ptr=>xdest
      ydest_ptr=>ydest
      zdest_ptr=>zdest
      call checkbound_array1(fablo,fabhi,xdest_ptr,0,0)
      call checkbound_array1(fablo,fabhi,ydest_ptr,0,1)
      call checkbound_array1(fablo,fabhi,zdest_ptr,0,SDIM-1)

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
        print *,"dir out of range in fort_fluidsolidcor"
        stop
       endif

       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,0,dir)

       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

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

         call gridstenMAC_level(xsten,i,j,k,level,nhalf,dir)

         if (dir.eq.0) then
          local_dest=xdest(D_DECL(i,j,k))
          local_src=xsrc(D_DECL(i,j,k))
          local_gp=xgp(D_DECL(i,j,k))
          AL=xface(D_DECL(i,j,k),FACECOMP_FACECUT+1)
          AL_ice=xface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
          AL_ice_mask=xface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)
          icrit=i
         else if (dir.eq.1) then
          local_dest=ydest(D_DECL(i,j,k))
          local_src=ysrc(D_DECL(i,j,k))
          local_gp=ygp(D_DECL(i,j,k))
          AL=yface(D_DECL(i,j,k),FACECOMP_FACECUT+1)
          AL_ice=yface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
          AL_ice_mask=yface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)
          icrit=j
         else if ((dir.eq.2).and.(SDIM.eq.3)) then
          local_dest=zdest(D_DECL(i,j,k))
          local_src=zsrc(D_DECL(i,j,k))
          local_gp=zgp(D_DECL(i,j,k))
          AL=zface(D_DECL(i,j,k),FACECOMP_FACECUT+1)
          AL_ice=zface(D_DECL(i,j,k),FACECOMP_ICEFACECUT+1)
          AL_ice_mask=zface(D_DECL(i,j,k),FACECOMP_ICEMASK+1)
          icrit=k
         else
          print *,"dir invalid fluid solid cor"
          stop
         endif

         if (project_option_singular_possibleF(project_option).eq.1) then
          if ((nsolve.eq.1).and.(velcomp.eq.0)) then
           veldir=1
           bc_comp=1
          else
           print *,"nsolve, or velcomp invalid"
           stop
          endif
         else if ((project_option.eq.SOLVETYPE_HEAT).or. & ! thermal diffusion
                  ((project_option.ge.SOLVETYPE_SPEC).and. & ! species
                   (project_option.lt.SOLVETYPE_SPEC+num_species_var))) then
          if ((nsolve.eq.1).and. &
              (velcomp.eq.0)) then
           veldir=1
           bc_comp=1
          else
           print *,"nsolve,or velcomp invalid"
           stop
          endif
         else if (project_option.eq.SOLVETYPE_VISC) then ! viscosity
          if ((nsolve.eq.SDIM).and. &
              (velcomp.ge.0).and. &
              (velcomp.lt.SDIM)) then
           veldir=velcomp+1
           bc_comp=veldir
          else
           print *,"nsolve,or velcomp invalid"
           stop
          endif

         else
          print *,"project_option invalid; fort_fluidsolidcor"
          stop
         endif

         local_dd=one
         local_visc_coef=one
         local_uncoupled_viscosity=1

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

          ! declared in: PROB.F90
         call eval_face_coeff( &
           xsten,nhalf, &
           level,finest_level, &
           AL, & !intent(in)
           AL_ice, & !intent(in)
           AL_ice_mask, & !intent(in)
           cc_group, & !intent(out)
           local_dd, & !intent(in)
           local_dd_group, & !intent(out)
           local_visc_coef, &
           nsolve,dir,veldir,project_option, &
           local_uncoupled_viscosity,side,bccrit,local_wt)

         if (local_wt(veldir).ge.zero) then

          if (local_wt(veldir).eq.zero) then
           if ((local_gp.eq.zero).or.(maskface.eq.0).or.(1.eq.0)) then
            local_dest=local_src
           else if ((local_gp.ne.zero).and.(maskface.eq.1)) then
            print *,"local_gp invalid"
            print *,"local_gp ",local_gp
            print *,"level,finest_level ",level,finest_level
            print *,"AL,AL_ice,AL_ice_mask,cc_group ", &
                  AL,AL_ice,AL_ice_mask,cc_group
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
      enddo ! dir=0..sdim-1

      return
      end subroutine fort_fluidsolidcor


       subroutine fort_diaginv( &
        diag_reg, &
        DIMS(diag_reg), &
        resid, DIMS(resid), &
        xnew, DIMS(xnew), &
        xold, DIMS(xold), &
        mask, DIMS(mask), &
        tilelo,tilehi, &
        fablo,fabhi,bfact ) &
       bind(c,name='fort_diaginv')

       use global_utility_module

       IMPLICIT NONE

       real(amrex_real) :: local_diag
       integer, INTENT(in) :: bfact
       integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
       integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
       integer :: growlo(3),growhi(3)
       integer, INTENT(in) :: DIMDEC(diag_reg)
       integer, INTENT(in) :: DIMDEC(resid)
       integer, INTENT(in) :: DIMDEC(xnew)
       integer, INTENT(in) :: DIMDEC(xold)
       integer, INTENT(in) :: DIMDEC(mask)
       real(amrex_real), INTENT(in),target :: diag_reg(DIMV(diag_reg))
       real(amrex_real), pointer :: diag_reg_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in),target :: resid(DIMV(resid))
       real(amrex_real), pointer :: resid_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(out),target :: xnew(DIMV(xnew))
       real(amrex_real), pointer :: xnew_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in),target :: xold(DIMV(xold))
       real(amrex_real), pointer :: xold_ptr(D_DECL(:,:,:))
       real(amrex_real), INTENT(in),target :: mask(DIMV(mask))
       real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))

       integer i,j,k

       if (bfact.lt.1) then
        print *,"bfact invalid158"
        stop
       endif

       mask_ptr=>mask
       diag_reg_ptr=>diag_reg
       xnew_ptr=>xnew
       xold_ptr=>xold
       resid_ptr=>resid
       call checkbound_array1(fablo,fabhi,mask_ptr,1,-1) 
       call checkbound_array1(fablo,fabhi,diag_reg_ptr,0,-1) 
       call checkbound_array1(fablo,fabhi,xnew_ptr,1,-1) 
       call checkbound_array1(fablo,fabhi,xold_ptr,1,-1) 
       call checkbound_array1(fablo,fabhi,resid_ptr,0,-1) 

       call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,0) 
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)

        ! mask=tag if not covered by level+1 or outside the domain.
        if (mask(D_DECL(i,j,k)).eq.one) then

         if (diag_reg(D_DECL(i,j,k)).eq.zero) then
          print *,"diag_reg invalid 1"
          stop
         else if (diag_reg(D_DECL(i,j,k)).gt.zero) then
          local_diag=diag_reg(D_DECL(i,j,k))
         else
          print *,"diag invalid 2"
          stop
         endif

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
       end subroutine fort_diaginv

      subroutine fort_io_compare( &
       nsteps, &
       do_input, &
       visual_compare, &
       time, &
       fabin,DIMS(fabin), & !initialized from COARSEDATA.tec
       fabout,DIMS(fabout), & !initialized by "ns_level.output_zones"
       vislo,vishi, &
       visual_ncomp) &
      bind(c,name='fort_io_compare')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      integer, INTENT(in) :: nsteps
      integer, INTENT(in) :: do_input
      integer, INTENT(in) :: visual_compare
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: visual_ncomp
      integer, INTENT(in) :: vislo(SDIM),vishi(SDIM) 

      integer, INTENT(in) :: DIMDEC(fabin)
      integer, INTENT(in) :: DIMDEC(fabout)
       ! x,u,pmg,den,T,Y1..Yn,mag vort,LS
      real(amrex_real), INTENT(inout),target :: &
              fabout(DIMV(fabout),visual_ncomp) 
      real(amrex_real), pointer :: fabout_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(inout),target :: fabin(DIMV(fabin),visual_ncomp)
      real(amrex_real), pointer :: fabin_ptr(D_DECL(:,:,:),:)

      real(amrex_real) problo(SDIM),probhi(SDIM)

      integer visual_ncell(SDIM)
      real(amrex_real) visual_dx(SDIM)
      integer gridlo(3),gridhi(3) 
      integer i,j,k
      integer im
      integer dir
      integer n
      integer strandid
      integer gridsum(num_materials)
      integer LSgridsum(num_materials)
      integer gridsum_total
      real(amrex_real) xtest(SDIM)
       ! x,u,pmg,den,T,Y1..Yn,mag vort,LS,umag
      real(amrex_real) local_data(visual_ncomp+1) 
      real(amrex_real) local_data_in(visual_ncomp+1)
      real(amrex_real) local_data_out(visual_ncomp+1)
      real(amrex_real) L1norm_local(visual_ncomp+1)
      real(amrex_real) L2norm_local(visual_ncomp+1)
      integer icrit(num_materials*(visual_ncomp+1))
      integer jcrit(num_materials*(visual_ncomp+1))
      integer kcrit(num_materials*(visual_ncomp+1))
      real(amrex_real) L1norm(num_materials*(visual_ncomp+1))
      real(amrex_real) L2norm(num_materials*(visual_ncomp+1))
      real(amrex_real) Linfnorm(num_materials*(visual_ncomp+1))
      real(amrex_real) LSnorm(num_materials)
      real(amrex_real) LS_in,LS_out
      integer ebase
      character*2 im_str
      character*6 stepstr
      character*14 compfilename
      character*17 uniformfilename
      character*255 techeader_str1
      character*255 techeader_str2
      character*255 techeader_str3

       ! x,u,pmg,den,T,Y1..Yn,mag vort,LS
      if (visual_ncomp.ne.VISUALCOMP_NCOMP) then
       print *,"visual_ncomp invalid" 
       stop
      endif
       ! x,u,pmg,den,T,Y1..Yn,mag vort,LS,umag
      if (visual_ncomp+1.ne.VISUALCOMP_NCOMP_MAGVEL) then
       print *,"visual_ncomp or VISUALCOMP_NCOMP_MAGVEL invalid" 
       stop
      endif

      if (time.ge.zero) then
       ! do nothing
      else
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

      fabout_ptr=>fabout
      fabin_ptr=>fabin

      call checkbound_array(vislo,vishi,fabout_ptr,0,0)
      call checkbound_array(vislo,vishi,fabout_ptr,0,1)
      call checkbound_array(vislo,vishi,fabout_ptr,0,SDIM-1)

      call checkbound_array(vislo,vishi,fabin_ptr,0,0)
      call checkbound_array(vislo,vishi,fabin_ptr,0,1)
      call checkbound_array(vislo,vishi,fabin_ptr,0,SDIM-1)

      if (do_input.eq.1) then

       if (visual_compare.eq.1) then

        write(compfilename,'(A14)') 'COARSEDATA.tec'
        print *,"compfilename ",compfilename
        open(unit=11,file=compfilename)
        read(11,*) techeader_str1  ! variables
        read(11,*) techeader_str2  ! zone map
        read(11,*) techeader_str3  ! SOLUTIONTIME, STRANDID

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
         do n=1,visual_ncomp
          if ((n.ge.1).and.(n.lt.visual_ncomp)) then
!          read(11,'(D25.16)',ADVANCE="NO") local_data(n)
           read(11,'(E25.16)',ADVANCE="NO") local_data(n)
          else if (n.eq.visual_ncomp) then
!          read(11,'(D25.16)') local_data(n)
           read(11,'(E25.16)') local_data(n)
          else
           print *,"n invalid"
           stop
          endif
         enddo ! n=1..visual_ncomp

         do dir=1,SDIM
          if (abs(local_data(dir)-xtest(dir)).gt.VOFTOL*visual_dx(dir)) then
           print *,"local_data(dir) invalid"
           stop
          endif
         enddo
         do n=1,visual_ncomp
          if (abs(local_data(n)).le.1.0D+20) then
           fabin(D_DECL(i,j,k),n)=local_data(n)
          else
           print *,"local_data(n) overflow"
           stop
          endif
         enddo ! n=1..visual_ncomp
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

       ! LINE 1: VARIABLES

       ! x,u,pmg,den,T,Y1..Yn,mag vort,LS
       if (SDIM.eq.2) then
         !61-11+1=51
        write(11,'(A51)',ADVANCE="NO")  &
         'VARIABLES="x_pos","y_pos","x_velocity","y_velocity"'
       else if (SDIM.eq.3) then
         !82-11+1=72
        write(11,'(A72)',ADVANCE="NO") &
         'VARIABLES="x_pos","y_pos","z_pos","x_velocity","y_velocity","z_velocity"'
       else
        print *,"dimension bust"
        stop
       endif
       write(11,'(A20)',ADVANCE="NO") ',"PRES","DEN","TEMP"'

       do im=1,num_species_var
        write(im_str,'(I2)') im
        do i=1,2
         if (im_str(i:i).eq.' ') then
          im_str(i:i)='0'
         endif
        enddo
        write(11,'(A6)',ADVANCE="NO") ',"SPEC'
        write(11,'(A2)',ADVANCE="NO") im_str
        write(11,'(A1)',ADVANCE="NO") '"'
       enddo ! im=1..num_species_var

       write(11,'(A9)',ADVANCE="NO") ',"MGVORT"'

       do im=1,num_materials
        write(im_str,'(I2)') im
        do i=1,2
         if (im_str(i:i).eq.' ') then
          im_str(i:i)='0'
         endif
        enddo
        write(11,'(A4)',ADVANCE="NO") ',"LS'
        write(11,'(A2)',ADVANCE="NO") im_str
        write(11,'(A1)',ADVANCE="NO") '"'
        if (im.eq.num_materials) then
         write(11,*) ' '
        endif
       enddo ! im=1..num_materials

       ! LINE 2: ZONE MAP

       if (SDIM.eq.2) then
        write(11,*)'zone i=',visual_ncell(1)+1, &
                   ' j=',visual_ncell(SDIM)+1, &
                   ' f=point '
       else if (SDIM.eq.3) then
        write(11,*)'zone i=',visual_ncell(1)+1, &
                   ' j=',visual_ncell(2)+1, &
                   ' k=',visual_ncell(SDIM)+1, &
                   ' f=point '
       else
        print *,"dimension bust"
        stop
       endif

       ! LINE 3: SOLUTIONTIME, STRANDID

       write(11,*) ' SOLUTIONTIME=',round_time(time),' STRANDID=',strandid

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
        do n=1,visual_ncomp
         local_data(n)=fabout(D_DECL(i,j,k),n)
        enddo
        do dir=1,SDIM
         if (abs(local_data(dir)-xtest(dir)).gt.VOFTOL*visual_dx(dir)) then
          print *,"local_data(dir) invalid"
          stop
         endif
        enddo
        do n=1,visual_ncomp
         if (abs(local_data(n)).le.1.0D+20) then
          if ((n.ge.1).and.(n.lt.visual_ncomp)) then
!          write(11,'(D25.16)',ADVANCE="NO") local_data(n)
           write(11,'(E25.16)',ADVANCE="NO") local_data(n)
          else if (n.eq.visual_ncomp) then
!          write(11,'(D25.16)') local_data(n)
           write(11,'(E25.16)') local_data(n)
          else
           print *,"n invalid"
           stop
          endif
         else
          print *,"local_data(n) overflow"
          stop
         endif
        enddo ! n=1..visual_ncomp
 
       enddo ! i
       enddo ! j
       enddo ! k

       close(11)

       if (visual_compare.eq.1) then

        do n=1,num_materials*(VISUALCOMP_NCOMP_MAGVEL)
         L1norm(n)=zero
         L2norm(n)=zero
         Linfnorm(n)=zero
         icrit(n)=-1
         jcrit(n)=-1
         kcrit(n)=-1
        enddo ! n=1..num_materials*(VISUALCOMP_NCOMP_MAGVEL)

        do im=1,num_materials
         gridsum(im)=0
         LSgridsum(im)=0
         LSnorm(im)=zero
        enddo

        gridsum_total=0

        do k=gridlo(3),gridhi(3)
        do j=gridlo(2),gridhi(2)
        do i=gridlo(1),gridhi(1)
          ! x,u,pmg,den,T,Y1..Yn,mag vort,LS
         do n=1,visual_ncomp
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
         enddo ! n=1..visual_ncomp

          ! x,u,pmg,den,T,Y1..Yn,mag vort,LS,umag
         n=VISUALCOMP_MAGVEL+1

          ! velocity magnitude.
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

         do im=1,num_materials
          LS_in=local_data_in(VISUALCOMP_LS+im)
          LS_out=local_data_out(VISUALCOMP_LS+im)
          if ((LS_in.gt.zero).and.(LS_out.gt.zero)) then
           ebase=(im-1)*(VISUALCOMP_NCOMP_MAGVEL)
           do n=1,VISUALCOMP_NCOMP_MAGVEL
            L1norm(ebase+n)=L1norm(ebase+n)+L1norm_local(n)
            L2norm(ebase+n)=L2norm(ebase+n)+L2norm_local(n)
            if (abs(local_data(n)).gt.Linfnorm(ebase+n)) then
             Linfnorm(ebase+n)=abs(local_data(n))
             icrit(ebase+n)=i
             jcrit(ebase+n)=j
             kcrit(ebase+n)=k
            endif
           enddo ! n=1..VISUALCOMP_NCOMP_MAGVEL
           gridsum(im)=gridsum(im)+1
           gridsum_total=gridsum_total+1
          else if ((LS_in.le.zero).or.(LS_out.le.zero)) then
           ! do nothing
          else
           print *,"LS_in or LS_out invalid"
           stop
          endif
         enddo ! im=1..num_materials
         do im=1,num_materials
          LS_in=local_data_in(VISUALCOMP_LS+im)
          LS_out=local_data_out(VISUALCOMP_LS+im)
          if ((abs(LS_in).le.visual_dx(1)).or. &
              (abs(LS_out).le.visual_dx(1))) then
           LSnorm(im)=LSnorm(im)+abs(local_data(VISUALCOMP_LS+im))
           LSgridsum(im)=LSgridsum(im)+1
          else if ((abs(LS_in).gt.visual_dx(1)).and. &
                   (abs(LS_out).gt.visual_dx(1))) then
           ! do nothing
          else
           print *,"LS_in or LS_out invalid"
           stop
          endif
         enddo ! im=1..num_materials
          
        enddo ! i
        enddo ! j
        enddo ! k

        if (gridsum_total.ge.1) then
         do im=1,num_materials
          if (gridsum(im).gt.0) then
           ebase=(im-1)*(VISUALCOMP_NCOMP_MAGVEL)
           do n=1,VISUALCOMP_NCOMP_MAGVEL
            L1norm(ebase+n)=L1norm(ebase+n)/gridsum(im)
            L2norm(ebase+n)=sqrt(L2norm(ebase+n)/gridsum(im))
            call print_visual_descriptor(im,n)
            print *,"im,n,cells,L1,L2,Linf,ic,jc,kc ", &
             im,n,gridsum(im), &
             L1norm(ebase+n),L2norm(ebase+n),Linfnorm(ebase+n), &
             icrit(ebase+n),jcrit(ebase+n),kcrit(ebase+n)
           enddo ! n=1..VISUALCOMP_NCOMP_MAGVEL
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
         enddo ! im=1..num_materials
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
      end subroutine fort_io_compare

      subroutine fort_outputslice( &
       time,nsteps,sliceint,slice_data,nslice,nstate_slice) &
      bind(c,name='fort_outputslice')

      use global_utility_module
      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: nsteps,nslice,nstate_slice,sliceint
      real(amrex_real), INTENT(in) :: slice_data(nslice*nstate_slice)
      integer strandid
      integer n
      integer i
      character*6 stepstr
      character*11 sfilename

        ! nstate_slice=x,y,z,xvel,yvel,zvel,PMG,PEOS,DIV,den,Temp,KE
        ! (value of material with LS>0)
      if (nstate_slice.ne.SLICECOMP_NCOMP) then
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

      ! to use gnuplot, put # at beginning of comment lines.
      write(11,*) '"CLSVOF data"'
      if (SDIM.eq.2) then
       write(11,*) 'VARIABLES="X","Y","U","V","PMG","PEOS","DIV","D","T","KE"'
      else if (SDIM.eq.3) then
       write(11,*) 'VARIABLES="X","Y","Z","U","V","W","PMG","PEOS","DIV","D","T","KE"'
      else
       print *,"dimension bust"
       stop
      endif
      write(11,*)'zone i=',nslice-1,' SOLUTIONTIME=',round_time(time), &
       ' STRANDID=',strandid

      do i=-1,nslice-2
       do n=1,nstate_slice-1
!       write(11,'(D25.16)',ADVANCE="NO")  &
        write(11,'(E25.16)',ADVANCE="NO")  &
          slice_data((i+1)*nstate_slice+n)
       enddo
       n=nstate_slice
!      write(11,'(D25.16)') slice_data((i+1)*nstate_slice+n)
       write(11,'(E25.16)') slice_data((i+1)*nstate_slice+n)
      enddo ! i


      close(11)
      return
      end subroutine fort_outputslice

      subroutine fort_zalesak_cell( &
        xlo,dx, &
        u,domlo,domhi, &
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        level, &
        DIMS(u), &
        time) &
      bind(c,name='fort_zalesak_cell')

      use global_utility_module
      use probf90_module

      IMPLICIT NONE
      integer, INTENT(in) :: bfact,level
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: domlo(SDIM), domhi(SDIM)
      integer, INTENT(in) :: DIMDEC(u)
      real(amrex_real), INTENT(inout),target :: u(DIMV(u),SDIM)
      real(amrex_real), pointer :: u_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      integer i,j,k,dir
      real(amrex_real) xx(SDIM)
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) LS
      real(amrex_real) vel(SDIM)
      real(amrex_real) temperature
      integer prescribed_flag

      if (bfact.lt.1) then
       print *,"bfact too small"
       stop
      endif
      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,1) 
      u_ptr=>u
      call checkbound_array(fablo,fabhi,u_ptr,1,-1) 

      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten_level(xsten,i,j,k,level,nhalf)

       do dir=1,SDIM
        xx(dir)=xsten(0,dir)
       enddo

       if (fort_is_passive_advect_test().eq.1) then
        call SUB_clamped_LS(xx,time,LS,vel,temperature,prescribed_flag,dx)
        if (prescribed_flag.eq.1) then
         if (LS.eq.CLAMPED_EVERYWHERE_LS) then
          do dir=1,SDIM
           u(D_DECL(i,j,k),dir)=vel(dir)
          enddo
         else
          print *,"expecting LS.eq.CLAMPED_EVERYWHERE_LS"
          stop
         endif
        else if (prescribed_flag.ne.1) then
         print *,"expecting prescribed_flag=1"
         stop
        endif
       else
        print *,"expecting fort_is_passive_advect_test().eq.1"
        stop
       endif

       do dir=1,SDIM
        u(D_DECL(i,j,k),dir)=u(D_DECL(i,j,k),dir)/global_velocity_scale
       enddo

      enddo
      enddo
      enddo

      return 
      end subroutine fort_zalesak_cell

       ! spectral_override==LOW_ORDER_AVGDOWN => always do low order
      subroutine fort_avgdown( &
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
       lof,hif) &
      bind(c,name='fort_avgdown')

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: spectral_override
      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level_f
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact_f
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: xlo_fine(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      integer, INTENT(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      integer mstenlo(3), mstenhi(3)
      real(amrex_real), INTENT(in),target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out),target :: crse(DIMV(crse),ncomp)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: fine(DIMV(fine),ncomp)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      integer flochi(SDIM)
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer ilocal,jlocal,klocal
      integer n
      integer dir2
      integer grid_type
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) crse_value(ncomp)
      real(amrex_real) dxelem_f
      real(amrex_real) xcoarse(SDIM)
      integer finelo_index
      real(amrex_real), dimension(D_DECL(:,:,:),:),allocatable :: ffine
      integer testmask,testmask2
      integer local_enable_spectral
      integer elem_test


      local_enable_spectral=enable_spectral

      if ((spectral_override.ne.LOW_ORDER_AVGDOWN).and. &
          (spectral_override.ne.SPECTRAL_ORDER_AVGDOWN)) then
       print *,"spectral_override invalid: ",spectral_override
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

      crse_ptr=>crse
      fine_ptr=>fine
      mask_ptr=>mask
      call checkbound_array(lo,hi,crse_ptr,0,-1)
      call checkbound_array(lof,hif,fine_ptr,0,-1)
      call checkbound_array1(lof,hif,mask_ptr,1,-1)

      grid_type=-1  ! ggg  (Gauss in all directions)

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

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

       call gridsten_level(xsten,ic,jc,kc,level_c,nhalf)

        ! find stencil for surrounding fine level spectral element.
       stenlo(3)=0
       stenhi(3)=0
       do dir2=1,SDIM
        dxelem_f=dxf(dir2)*bfact_f
        xcoarse(dir2)=xsten(0,dir2)
        finelo_index=NINT( (xcomp_of_xphys(dir2-1,xcoarse(dir2))- &
                            problo(dir2))/dxelem_f-half )
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
         if (abs(xcoarse(dir2)).lt.EPS4*dxf(dir2)) then
          mstenlo(dir2)=stenlo(dir2)-1
         endif
         if (abs(xcoarse(dir2)-bfact_f*dxf(dir2)).lt. &
             EPS4*dxf(dir2)) then
          mstenhi(dir2)=stenhi(dir2)+1
         endif
        endif ! stenhi<>stenlo

        if ((xcoarse(dir2).ge.-EPS3*dxf(dir2)).and. &
            (xcoarse(dir2).le.(EPS3+bfact_f)*dxf(dir2))) then
         !do nothing
        else
         print *,"xcoarse out of bounds avgdown: ",dir2,xcoarse(dir2)
         stop
        endif

       enddo ! dir2=1..sdim

       ifine=mstenlo(1) 
       jfine=mstenlo(2) 
       kfine=mstenlo(SDIM) 
       testmask=NINT(mask(D_DECL(ifine,jfine,kfine)))
       do kfine=mstenlo(3),mstenhi(3) 
       do jfine=mstenlo(2),mstenhi(2) 
       do ifine=mstenlo(1),mstenhi(1) 
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
           (local_enable_spectral.eq.1).and. &
           (testmask.ge.1).and. &
           (testmask.le.num_materials).and. &
           (spectral_override.eq.SPECTRAL_ORDER_AVGDOWN)) then

        do kfine=stenlo(3),stenhi(3)
        do jfine=stenlo(2),stenhi(2)
        do ifine=stenlo(1),stenhi(1)
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
         ncomp,bfact_f,grid_type, &
         flochi,dxf,xcoarse,ffine,crse_value)

        voltotal=one

       else if ((bfact_f.eq.1).or. &
                (local_enable_spectral.eq.0).or. &
                (testmask.eq.0).or. &
                (spectral_override.eq.LOW_ORDER_AVGDOWN)) then

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

       if (voltotal.gt.zero) then
        !do nothing
       else
        print *,"voltotal invalid: ",voltotal
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
      end subroutine fort_avgdown

      subroutine fort_avgdown_low( &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) &
      bind(c,name='fort_avgdown_low')

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level_f
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact_f
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: xlo_fine(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      integer, INTENT(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real), INTENT(out),target :: crse(DIMV(crse),ncomp)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: fine(DIMV(fine),ncomp)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer n
      integer dir2
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      real(amrex_real) crse_value(ncomp)

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

      crse_ptr=>crse
      fine_ptr=>fine
      call checkbound_array(lo,hi,crse_ptr,0,-1)
      call checkbound_array(lof,hif,fine_ptr,0,-1)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

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

       if (voltotal.gt.zero) then
        !do nothing
       else
        print *,"voltotal invalid: ",voltotal
        stop
       endif

       do n = 1,ncomp
        crse(D_DECL(ic,jc,kc),n)=crse_value(n)/voltotal
       enddo

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine fort_avgdown_low

       ! NavierStokes::level_phase_change_redistribute
       ! NavierStokes::avgDown_tag_localMF
       ! NavierStokes::level_avgDown_tag
      subroutine fort_avgdown_tag( &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) &
      bind(c,name='fort_avgdown_tag')

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level_f
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact_f
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: xlo_fine(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      integer, INTENT(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real), INTENT(out),target :: crse(DIMV(crse),ncomp)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: fine(DIMV(fine),ncomp)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer n
      integer dir2
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      real(amrex_real) crse_value(ncomp)
      integer fine_test
      integer coarse_test
      integer init_flag

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

      crse_ptr=>crse
      fine_ptr=>fine
      call checkbound_array(lo,hi,crse_ptr,0,-1)
      call checkbound_array(lof,hif,fine_ptr,0,-1)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

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

       if (voltotal.gt.zero) then
        !do nothing
       else
        print *,"voltotal invalid: ",voltotal
        stop
       endif

       do n = 1,ncomp
        crse(D_DECL(ic,jc,kc),n)=crse_value(n)
       enddo

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine fort_avgdown_tag

        ! icurv=(iten-1)*CURVCOMP_NCOMP
        ! dir=1..sdim
        ! side=-1 or 1
        ! curvfab(D_DECL(i,j,k),icurv+CURVCOMP_DIRSIDE_FLAG+1)=dir*side
      subroutine fort_avgdown_curv( &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       ncomp_curv_total, & !num_interfaces*CURVCOMP_NCOMP
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi, &
       lof,hif) &
      bind(c,name='fort_avgdown_curv')

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level_f
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact_f
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: xlo_fine(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: ncomp_curv_total !num_interfaces*CURVCOMP_NCOMP
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM) ! coarse grid dimensions
      integer, INTENT(in) :: lof(SDIM),hif(SDIM) ! fine grid dimensions
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real), INTENT(out),target :: crse(DIMV(crse),ncomp_curv_total)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: fine(DIMV(fine),ncomp_curv_total)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer dir2
      integer iten
      integer n
      real(amrex_real) voltotal
      real(amrex_real) volall
      real(amrex_real) wt(SDIM)
      real(amrex_real) crse_value(ncomp_curv_total)
      integer fine_test
      integer coarse_test
      real(amrex_real) velwt(num_interfaces)
      integer icurv
      integer idx_dirside

      if (ncomp_curv_total.ne.num_interfaces*CURVCOMP_NCOMP) then
       print *,"ncomp_curv_total invalid35"
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

      crse_ptr=>crse
      fine_ptr=>fine
      call checkbound_array(lo,hi,crse_ptr,0,-1)
      call checkbound_array(lof,hif,fine_ptr,0,-1)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

       do n=1,ncomp_curv_total
        crse_value(n)=zero
       enddo
       do iten=1,num_interfaces
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
             if (volall.gt.zero) then
              ! do nothing
             else
              print *,"volall invalid"
              stop
             endif

             do iten=1,num_interfaces
              icurv=(iten-1)*CURVCOMP_NCOMP
              ! dir x side  dir=1..sdim side=-1 or 1
              idx_dirside=icurv+CURVCOMP_DIRSIDE_FLAG+1 
              fine_test=NINT(fine(D_DECL(ifine,jfine,kfine),idx_dirside))
              coarse_test=NINT(crse_value(idx_dirside))
              if ((coarse_test.ne.0).and. &
                  (coarse_test.ne.SDIM+1)) then
               print *,"coarse_test invalid"
               stop
              endif
              if (fine_test.eq.0) then
               ! do nothing

               ! fine_test can equal SDIM+1 if the finer level (level l+1)
               ! being averaged down is not the finest level. 
               ! (i.e. level l+1 was obtained from averaging down level l+2)
              else if ((abs(fine_test).ge.1).and. &
                       (abs(fine_test).le.SDIM+1)) then
               velwt(iten)=velwt(iten)+volall
               crse_value(idx_dirside)=SDIM+1
               crse_value(icurv+CURVCOMP_MATERIAL3_ID+1)=zero 
               do dir2=1,CURVCOMP_DIRSIDE_FLAG
                crse_value(icurv+dir2)= &
                  crse_value(icurv+dir2)+volall* &
                  fine(D_DECL(ifine,jfine,kfine),icurv+dir2)
               enddo
              else
               print *,"fine_test bad (avgdown curv): ",fine_test
               stop
              endif
             enddo ! iten=1...num_interfaces

             voltotal=voltotal+volall
            endif ! wt(sdim).gt.0
           enddo ! kfine
          endif
         enddo ! jfine
        endif
       enddo ! ifine

       if (voltotal.gt.zero) then
        ! do nothing
       else
        print *,"voltotal invalid: ",voltotal
        stop
       endif

       do iten=1,num_interfaces
        icurv=(iten-1)*CURVCOMP_NCOMP
        ! dir x side  dir=1..sdim side=-1 or 1
        idx_dirside=icurv+CURVCOMP_DIRSIDE_FLAG+1 

        if ((idx_dirside.ge.1).and. &
            (idx_dirside.le.num_interfaces*CURVCOMP_NCOMP)) then
         coarse_test=NINT(crse_value(idx_dirside))
        else
         print *,"idx_dirside invalid: ",idx_dirside
         print *,"iten=",iten
         print *,"icurv=",icurv
         print *,"CURVCOMP_DIRSIDE_FLAG ",CURVCOMP_DIRSIDE_FLAG
         print *,"CURVCOMP_NCOMP=",CURVCOMP_NCOMP
         print *,"ncomp_curv_total: ",ncomp_curv_total
         print *,"num_interfaces=",num_interfaces
         stop
        endif

        if (coarse_test.eq.0) then
         if (velwt(iten).eq.zero) then
          ! do nothing
         else
          print *,"velwt invalid"
          stop
         endif
         do dir2=1,CURVCOMP_NCOMP
          crse(D_DECL(ic,jc,kc),icurv+dir2)=zero
         enddo
        else if (coarse_test.eq.SDIM+1) then
         if (velwt(iten).gt.zero) then
          crse(D_DECL(ic,jc,kc),idx_dirside)=SDIM+1 ! dir x side
          crse(D_DECL(ic,jc,kc),icurv+CURVCOMP_MATERIAL3_ID+1)=zero 
          do dir2=1,CURVCOMP_DIRSIDE_FLAG
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

       enddo ! iten=1..num_interfaces

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine fort_avgdown_curv


      subroutine fort_mofavgdown ( &
       time, &
       problo, &
       dxc, &
       dxf, &
       bfact_c,bfact_f, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi) &
      bind(c,name='fort_mofavgdown')

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probcommon_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: bfact_c,bfact_f
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
      integer  growlo(3),growhi(3)
      integer  stenlo(3),stenhi(3)
      real(amrex_real), INTENT(out),target :: &
              crse(DIMV(crse),num_materials*ngeom_raw)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
              fine(DIMV(fine),num_materials*ngeom_raw)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      integer ifine,jfine,kfine
      integer ic,jc,kc
      integer dir
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxc(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      integer domlo(SDIM)

      integer :: grid_index(SDIM)
      integer, parameter :: grid_level=-1

      integer, parameter :: nhalf=3
      integer, parameter :: nhalfgrid=1

      integer n_overlap
 
      integer im
      integer vofcomp_raw
      integer vofcomp_recon
      real(amrex_real) wt(SDIM)
      real(amrex_real) testwt

      real(amrex_real) temp_vfrac
      real(amrex_real) temp_cen(SDIM)

      real(amrex_real) xstencoarse(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenfine(-nhalf:nhalf,SDIM)
      real(amrex_real) xstengrid(-nhalfgrid:nhalfgrid,SDIM)

      real(amrex_real) mofdatafine(num_materials*ngeom_recon)
      real(amrex_real) mofdatacoarse(num_materials*ngeom_recon)
      real(amrex_real) multi_centroidA(num_materials,SDIM)

      real(amrex_real) volcoarse
      real(amrex_real) cencoarse(SDIM)
      real(amrex_real) volfine
      real(amrex_real) cenfine(SDIM)

      real(amrex_real) vof_super(num_materials)
      real(amrex_real) multi_volume(num_materials)
      real(amrex_real) multi_cen(SDIM,num_materials)
      integer tessellate
      integer, parameter :: continuous_mof=STANDARD_MOF
      integer, parameter :: mof_verbose=0
      integer, parameter :: use_ls_data=0
      integer fine_covered
      real(amrex_real) LS_stencil(D_DECL(-1:1,-1:1,-1:1),num_materials)
      integer nmax

      integer cmofsten(D_DECL(-1:1,-1:1,-1:1))

      integer tid
#ifdef _OPENMP
      integer omp_get_thread_num
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

      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid"
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

      crse_ptr=>crse
      fine_ptr=>fine
      call checkbound_array(lo,hi,crse_ptr,0,-1)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

       grid_index(1)=ic
       grid_index(2)=jc
       if (SDIM.eq.3) then
        grid_index(SDIM)=kc
       endif

        ! coarse centroids and volume fractions are initialized
        ! to zero.
       do dir=1,num_materials*ngeom_raw
        crse(D_DECL(ic,jc,kc),dir) = zero
       enddo
       do dir=1,num_materials*ngeom_recon
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

              do im=1,num_materials
               vofcomp_raw=(im-1)*ngeom_raw+1
               vofcomp_recon=(im-1)*ngeom_recon+1
               do dir=1,SDIM+1
                mofdatafine(vofcomp_recon+dir-1)= &
                 fine(D_DECL(ifine,jfine,kfine),vofcomp_raw+dir-1)
               enddo
               do dir=SDIM+2,ngeom_recon
                mofdatafine(vofcomp_recon+dir-1)=zero
               enddo
              enddo ! im=1,num_materials

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

               do im=1,num_materials
                vofcomp_raw=(im-1)*ngeom_raw+1
                multi_volume(im)= &
                 fine(D_DECL(ifine,jfine,kfine),vofcomp_raw)*volfine
                do dir=1,SDIM
                 multi_cen(dir,im)= &
                  fine(D_DECL(ifine,jfine,kfine),vofcomp_raw+dir)+cenfine(dir)
                enddo
               enddo ! im

              else if (fine_covered.eq.0) then
               tessellate=0

                ! sum F_fluid=1  sum F_solid<=1
               call make_vfrac_sum_ok_base( &
                 cmofsten, &
                 xstenfine,nhalf, &
                 continuous_mof, &
                 bfact_f,dxf, &
                 tessellate, &
                 mofdatafine, &
                 SDIM)

               do im=1,num_materials
                vofcomp_recon=(im-1)*ngeom_recon+1
                vof_super(im)=mofdatafine(vofcomp_recon)
               enddo

               call multimaterial_MOF( &
                bfact_f,dxf,xstenfine,nhalf, &
                mof_verbose, & ! =0
                use_ls_data, & ! =0
                LS_stencil, &
                geom_xtetlist(1,1,1,tid+1), &
                geom_xtetlist(1,1,1,tid+1), &
                nmax, &
                nmax, &
                mofdatafine, & !intent(inout)
                vof_super, &
                multi_centroidA, &
                continuous_mof, & !=STANDARD_MOF
                cmofsten, &
                grid_index, &
                grid_level, &
                SDIM)

                ! EPS2
               call multi_get_volume_grid_simple( &
                EPS2, &
                tessellate, &  !=0
                bfact_f,dxf,xstenfine,nhalf, &
                mofdatafine, &
                xstengrid,nhalfgrid, &
                multi_volume,multi_cen, &
                geom_xtetlist(1,1,1,tid+1), &
                nmax, &
                nmax, &
                SDIM)

              else
               print *,"fine_covered invalid"
               stop
              endif

              do im=1,num_materials
               vofcomp_recon=(im-1)*ngeom_recon+1
               mofdatacoarse(vofcomp_recon)=mofdatacoarse(vofcomp_recon)+ &
                multi_volume(im)
               do dir=1,SDIM
                mofdatacoarse(vofcomp_recon+dir)= &
                 mofdatacoarse(vofcomp_recon+dir)+ &
                 multi_cen(dir,im)*multi_volume(im)
               enddo
               if (is_rigid(im).eq.0) then
                volcoarse=volcoarse+multi_volume(im)
                do dir=1,SDIM
                 cencoarse(dir)=cencoarse(dir)+ &
                  multi_cen(dir,im)*multi_volume(im)
                enddo
               else if (is_rigid(im).eq.1) then
                ! do nothing
               else
                print *,"is_rigid(im) invalid"
                stop
               endif
              enddo ! im=1..num_materials
  
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

       do im=1,num_materials
        vofcomp_recon=(im-1)*ngeom_recon+1
        vofcomp_raw=(im-1)*ngeom_raw+1

        temp_vfrac=mofdatacoarse(vofcomp_recon)/volcoarse
        if ((temp_vfrac.ge.zero).and. &
            (temp_vfrac.le.one+EPS1)) then
         !do nothing
        else
         print *,"temp_vfrac invalid: ",temp_vfrac
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
      end subroutine fort_mofavgdown


      subroutine fort_refine_density_avgdown ( &
       time, &
       problo, &
       dxc, &
       dxf, &
       bfact_c,bfact_f, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi) &
      bind(c,name='fort_refine_density_avgdown')

      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module
      use probcommon_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: dxc(SDIM)
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: bfact_c,bfact_f
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
      integer  growlo(3),growhi(3)
      integer  stenlo(3),stenhi(3)
      real(amrex_real), INTENT(out),target :: &
              crse(DIMV(crse),ENUM_NUM_REFINE_DENSITY_TYPE)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: &
              fine(DIMV(fine),ENUM_NUM_REFINE_DENSITY_TYPE)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      integer ifine,jfine,kfine
      integer ifine2,jfine2,kfine2
      integer ic,jc,kc
      integer ic2,jc2,kc2
      integer ncrse,nfine

      real(amrex_real) wt(SDIM)
      real(amrex_real) testwt
      real(amrex_real) wttotal
      real(amrex_real) crse_value


      integer tid
#ifdef _OPENMP
      integer omp_get_thread_num
#endif

      tid=0       
#ifdef _OPENMP
      tid=omp_get_thread_num()
#endif
      if ((tid.ge.geom_nthreads).or.(tid.lt.0)) then
       print *,"tid invalid"
       stop
      endif 

      if (time.ge.zero) then
       ! do nothing
      else
       print *,"time invalid: ",time
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

      crse_ptr=>crse
      fine_ptr=>fine
      call checkbound_array(lo,hi,crse_ptr,0,-1)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

       call fine_subelement_stencil(ic,jc,kc,stenlo,stenhi, &
         bfact_c,bfact_f)

       kc2=0
#if (AMREX_SPACEDIM==3)
       do kc2=0,1
#endif
       do jc2=0,1
       do ic2=0,1

        ncrse=4*kc2+2*jc2+ic2+1
        wttotal=zero
        crse_value=zero

        do ifine=stenlo(1),stenhi(1)
        do ifine2=0,1
         call intersect_weight_avg_refine( &
           ic,ic2,ifine,ifine2,  &
           bfact_c,bfact_f,wt(1))
         if (wt(1).gt.zero) then
          do jfine=stenlo(2),stenhi(2)
          do jfine2=0,1
           call intersect_weight_avg_refine( &
             jc,jc2,jfine,jfine2,  &
             bfact_c,bfact_f,wt(2))
           if (wt(2).gt.zero) then
            do kfine=stenlo(3),stenhi(3)
            kfine2=0
#if (AMREX_SPACEDIM==3)
            do kfine2=0,1
#endif
             if (SDIM.eq.3) then
              call intersect_weight_avg_refine( &
               kc,kc2,kfine,kfine2,  &
               bfact_c,bfact_f,wt(SDIM))
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
               nfine=4*kfine2+2*jfine2+ifine2+1
               crse_value = crse_value +  &
                testwt*fine(D_DECL(ifine,jfine,kfine),nfine)
               wttotal=wttotal+testwt
              else if (testwt.eq.zero) then
               print *,"testwt should not be 0"
               stop
              else
               print *,"testwt invalid: ",testwt
               stop
              endif
             endif ! wt(sdim).gt.0
#if (AMREX_SPACEDIM==3)
            enddo ! kfine2
#endif
            enddo ! kfine
           endif ! wt(2)>0
          enddo ! jfine2
          enddo ! jfine
         endif ! wt(1)>0
        enddo ! ifin2
        enddo ! ifine

        if (wttotal.gt.zero) then
         !do nothing
        else
         print *,"wttotal invalid: ",wttotal
         stop
        endif
        crse(D_DECL(ic,jc,kc),ncrse)=crse_value/wttotal
       enddo !ic2
       enddo !jc2
#if (AMREX_SPACEDIM==3)
       enddo !kc2
#endif

      enddo
      enddo
      enddo ! ic,jc,kc

      return
      end subroutine fort_refine_density_avgdown

      subroutine fort_erroravgdown ( &
       problo, &
       dxf, &
       bfact_c,bfact_f, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       lo,hi) &
      bind(c,name='fort_erroravgdown')

      use global_utility_module
      use geometry_intersect_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: bfact_c,bfact_f
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: lo(SDIM),hi(SDIM)
      integer growlo(3),growhi(3)
      integer stenlo(3),stenhi(3)
      real(amrex_real), INTENT(out),target :: crse(DIMV(crse))
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: fine(DIMV(fine))
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:))
      integer i,j,k,ic,jc,kc
      real(amrex_real), INTENT(in) :: problo(SDIM),dxf(SDIM)
 
      real(amrex_real) wt(SDIM)

      integer sanity_sum
      real(amrex_real) test_error,max_error

      if (bfact_f.lt.1) then
       print *,"bfact_f invalid5 ",bfact_f
       stop
      endif
      if ((bfact_c.ne.bfact_f).and. &
          (bfact_c.ne.2*bfact_f)) then
       print *,"bfact_c invalid"
       stop
      endif
      crse_ptr=>crse
      fine_ptr=>fine
      call checkbound_array1(lo,hi,crse_ptr,0,-1)

      call growntilebox(lo,hi,lo,hi,growlo,growhi,0) 

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

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
      end subroutine fort_erroravgdown



! the error is set to zero in advection.  Then
! the error is updated when the slopes are recomputed.
! this routine is called after the last slope computation of the time step.
      subroutine fort_pressure_indicator( &
        pressure_error_flag, &
        vorterr, &
        pressure_error_cutoff, &
        temperature_error_cutoff, &
        xlo,dx, &
        errnew,DIMS(errnew), &
        LS,DIMS(LS), &
        den,DIMS(den), &
        vort,DIMS(vort), &
        pres,DIMS(pres), &
        maskcov,DIMS(maskcov), &
        tilelo,tilehi,  &
        fablo,fabhi, &
        bfact, &
        level,  &
        finest_level) &
      bind(c,name='fort_pressure_indicator')

      use global_utility_module
      use probf90_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: pressure_error_flag
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      real(amrex_real), INTENT(in) :: xlo(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)

      real(amrex_real), INTENT(in) :: vorterr(num_materials)
      real(amrex_real), INTENT(in) :: pressure_error_cutoff(num_materials)
      real(amrex_real), INTENT(in) :: temperature_error_cutoff(num_materials)
      integer, INTENT(in) :: DIMDEC(maskcov)
      integer, INTENT(in) :: DIMDEC(errnew)
      integer, INTENT(in) :: DIMDEC(LS)
      integer, INTENT(in) :: DIMDEC(den)
      integer, INTENT(in) :: DIMDEC(vort)
      integer, INTENT(in) :: DIMDEC(pres)
      real(amrex_real), INTENT(in),target :: maskcov(DIMV(maskcov))
      real(amrex_real), pointer :: maskcov_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(inout),target :: errnew(DIMV(errnew))
      real(amrex_real), pointer :: errnew_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: &
              den(DIMV(den),num_materials*num_state_material)
      real(amrex_real), pointer :: den_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: LS(DIMV(LS),num_materials*(1+SDIM))
      real(amrex_real), pointer :: LS_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: vort(DIMV(vort))
      real(amrex_real), pointer :: vort_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in),target :: pres(DIMV(pres))
      real(amrex_real), pointer :: pres_ptr(D_DECL(:,:,:))

      integer i,j,k,im
      real(amrex_real) LStest(num_materials)
      real(amrex_real) pres_array(D_DECL(3,3,3))
      real(amrex_real) temp_array(D_DECL(3,3,3))
      integer i2,j2,k2
      integer kstencil_lo,kstencil_hi
      integer tcomp
      real(amrex_real) local_vort
      integer local_mask
      real(amrex_real) DXMAXLS
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)

      call get_dxmaxLS(dx,bfact,DXMAXLS)

      if (bfact.lt.1) then
       print *,"bfact invalid163"
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
      maskcov_ptr=>maskcov
      errnew_ptr=>errnew
      den_ptr=>den
      LS_ptr=>LS
      vort_ptr=>vort
      pres_ptr=>pres

      call checkbound_array1(fablo,fabhi,maskcov_ptr,1,-1)
      call checkbound_array(fablo,fabhi,LS_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,errnew_ptr,1,-1)
      call checkbound_array(fablo,fabhi,den_ptr,1,-1)
      call checkbound_array1(fablo,fabhi,vort_ptr,0,-1)
      call checkbound_array1(fablo,fabhi,pres_ptr,1,-1)

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
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)

       ! maskcov=tag if not covered by level+1 or outside the domain.
       local_mask=NINT(maskcov(D_DECL(i,j,k)))

       if (local_mask.eq.1) then

        call gridsten_level(xsten,i,j,k,level,nhalf)

        do im=1,num_materials
         LStest(im)=LS(D_DECL(i,j,k),im)
        enddo
        call get_primary_material(LStest,im)

        if (is_rigid(im).eq.0) then

           ! only check pressure/temperature/vorticity
           ! magnitude away from interfaces
         if (LStest(im).gt.DXMAXLS) then

          tcomp=(im-1)*num_state_material+ENUM_TEMPERATUREVAR+1

          if (vorterr(im).lt.zero) then
           print *,"vorterr cannot be negative"
           stop
          endif

          do k2=kstencil_lo,kstencil_hi
          do j2=-1,1
          do i2=-1,1
           pres_array(D_DECL(i2+2,j2+2,k2+2))= &
             pres(D_DECL(i+i2,j+j2,k+k2))
           temp_array(D_DECL(i2+2,j2+2,k2+2))= &
             den(D_DECL(i+i2,j+j2,k+k2),tcomp)
          enddo
          enddo
          enddo

           ! vort is initialized in fort_getshear when only_scalar.eq.2
           ! vort is the vorticity magnitude (L2 norm)
          local_vort=vort(D_DECL(i,j,k))

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
         else if (LStest(im).le.DXMAXLS) then
          ! do nothing
         else
          print *,"LStest(im) is NaN"
          stop
         endif  

        else if (is_rigid(im).eq.1) then
         ! do nothing
        else
         print *,"is_rigid invalid NAVIERSTOKES_3D.F90"
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
      end subroutine fort_pressure_indicator



       ! spectral_override==LOW_ORDER_AVGDOWN => always low order
      subroutine fort_edgeavgdown( &
       enable_spectral, &
       finest_level, &
       spectral_override, &
       problo, &
       dxf, &
       level_c,level_f, &
       bfact_c,bfact_f, &
       xlo_fine,dx, &
       grid_type, &
       crse,DIMS(crse), &
       fine,DIMS(fine), &
       mask,DIMS(mask), &
       loc,hic, &
       lof,hif, &
       ncomp) &
      bind(c,name='fort_edgeavgdown')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: enable_spectral
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: spectral_override
      integer, INTENT(in) :: level_c
      integer, INTENT(in) :: level_f
      integer, INTENT(in) :: bfact_c
      integer, INTENT(in) :: bfact_f
      real(amrex_real), INTENT(in) :: problo(SDIM)
      real(amrex_real), INTENT(in) :: dxf(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: xlo_fine(SDIM)
      integer, INTENT(in) :: ncomp
      integer, INTENT(in) :: DIMDEC(crse)
      integer, INTENT(in) :: DIMDEC(fine)
      integer, INTENT(in) :: DIMDEC(mask)
      integer, INTENT(in) :: loc(SDIM), hic(SDIM)  ! cell centered
      integer, INTENT(in) :: lof(SDIM), hif(SDIM)  ! cell centered
      integer growlo(3), growhi(3)
      integer stenlo(3), stenhi(3)
      integer mstenlo(3), mstenhi(3)
      integer, INTENT(in) :: grid_type
      integer dir2
      real(amrex_real), INTENT(out),target :: crse(DIMV(crse),ncomp)
      real(amrex_real), pointer :: crse_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: fine(DIMV(fine),ncomp)
      real(amrex_real), pointer :: fine_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),target :: mask(DIMV(mask))
      real(amrex_real), pointer :: mask_ptr(D_DECL(:,:,:))

      integer flochi(SDIM)
      integer ic,jc,kc
      integer ifine,jfine,kfine
      integer ilocal,jlocal,klocal
      integer inorm
      integer n
      integer infab
      real(amrex_real) voltotal
      real(amrex_real) volall

      real(amrex_real) wt(SDIM)

      integer, parameter :: nhalf=1
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      real(amrex_real) xstenND(-nhalf:nhalf,SDIM)
      real(amrex_real) crse_value(ncomp)
      real(amrex_real) dxelem_f
      real(amrex_real) xcoarse(SDIM)
      integer finelo_index
      real(amrex_real), dimension(D_DECL(:,:,:),:),allocatable :: ffine
      real(amrex_real) icemask
      integer khi
      integer testmask,testmask2
      integer local_enable_spectral
      integer box_type(SDIM)

      local_enable_spectral=enable_spectral

      if (level_f.gt.finest_level) then
       print *,"level_f invalid"
       stop
      endif

      if ((spectral_override.ne.LOW_ORDER_AVGDOWN).and. &
          (spectral_override.ne.ICEMASK_AVGDOWN).and. &
          (spectral_override.ne.SPECTRAL_ORDER_AVGDOWN)) then
       print *,"spectral_override invalid: ",spectral_override
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
      if ((grid_type.lt.-1).or.(grid_type.gt.5)) then
       print *,"grid_type invalid edge average down"
       stop
      endif

      call grid_type_to_box_type(grid_type,box_type)

      do dir2=1,SDIM
       flochi(dir2)=bfact_f-1+box_type(dir2)
      enddo

      if (SDIM.eq.2) then
       khi=0
      else if (SDIM.eq.3) then
       khi=flochi(SDIM)
      else
       print *,"dimension bust"
       stop
      endif

      allocate(ffine(D_DECL(0:flochi(1),0:flochi(2),0:flochi(3)),ncomp))

      crse_ptr=>crse
      fine_ptr=>fine
      mask_ptr=>mask

      call checkbound_array(loc,hic,crse_ptr,0,grid_type)
      call checkbound_array(lof,hif,fine_ptr,0,grid_type)
      call checkbound_array1(lof,hif,mask_ptr,1,-1)

      call growntileboxMAC(loc,hic,loc,hic,growlo,growhi,0,grid_type)

      do kc=growlo(3),growhi(3)
      do jc=growlo(2),growhi(2)
      do ic=growlo(1),growhi(1)

       call gridstenMAC_level(xsten,ic,jc,kc,level_c,nhalf,grid_type)

        ! find stencil for surrounding fine level spectral element.
       stenlo(3)=0
       stenhi(3)=0
       do dir2=1,SDIM
        dxelem_f=dxf(dir2)*bfact_f
        xcoarse(dir2)=xsten(0,dir2)
        finelo_index=NINT( (xcomp_of_xphys(dir2-1,xcoarse(dir2))- &
                            problo(dir2))/dxelem_f-half )
        stenlo(dir2)=bfact_f*finelo_index
        stenhi(dir2)=stenlo(dir2)+bfact_f-1+box_type(dir2)
       enddo ! dir

        ! check if coarse face is on a coarse element boundary.
       do dir2=1,SDIM
        if (box_type(dir2).eq.1) then
         if (dir2.eq.1) then
          inorm=ic
         else if (dir2.eq.2) then
          inorm=jc
         else if ((dir2.eq.SDIM).and.(SDIM.eq.3)) then
          inorm=kc
         else
          print *,"dir2 invalid"
          stop
         endif
         if ((inorm/bfact_c)*bfact_c.eq.inorm) then
          stenlo(dir2)=2*inorm
          stenhi(dir2)=stenlo(dir2)
         endif
        else if (box_type(dir2).eq.0) then
         ! do nothing
        else
         print *,"box_type(dir2) invalid"
         stop
        endif
       enddo ! dir2=1..sdim

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
         if (abs(xcoarse(dir2)).lt.EPS4*dxf(dir2)) then
          mstenlo(dir2)=stenlo(dir2)-1
         endif
         if (abs(xcoarse(dir2)-bfact_f*dxf(dir2)).lt. &
             EPS4*dxf(dir2)) then
          mstenhi(dir2)=stenhi(dir2)+1
         endif
        endif ! stenhi<>stenlo

        if ((xcoarse(dir2).ge.-EPS3*dxf(dir2)).and. &
            (xcoarse(dir2).le.(EPS3+bfact_f)*dxf(dir2))) then
         !do nothing
        else
         print *,"xcoarse out of bounds edge avgdown: ",dir2,xcoarse(dir2)
         stop
        endif

       enddo ! dir2=1..sdim

       ifine=mstenlo(1) 
       jfine=mstenlo(2) 
       kfine=mstenlo(SDIM) 
       testmask=NINT(mask(D_DECL(ifine,jfine,kfine)))
       do kfine=mstenlo(3),mstenhi(3) 
       do jfine=mstenlo(2),mstenhi(2) 
       do ifine=mstenlo(1),mstenhi(1) 
        testmask2=NINT(mask(D_DECL(ifine,jfine,kfine)))
        if (testmask2.ne.testmask) then
         testmask=0
        endif
       enddo 
       enddo 
       enddo 

       do klocal=0,khi
       do jlocal=0,flochi(2)
       do ilocal=0,flochi(1)
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
           (local_enable_spectral.eq.1).and. &
           (testmask.ge.1).and. &
           (testmask.le.num_materials).and. &
           (spectral_override.eq.SPECTRAL_ORDER_AVGDOWN)) then

        if ((grid_type.ge.0).and.(grid_type.lt.SDIM)) then
         ! do nothing
        else
         print *,"grid_type must be 0,1,or sdim for SEM"
         stop
        endif

        do kfine=stenlo(3),stenhi(3)
        do jfine=stenlo(2),stenhi(2)
        do ifine=stenlo(1),stenhi(1)
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
         ncomp,bfact_f,grid_type, &
         flochi,dxf,xcoarse,ffine,crse_value)

        voltotal=one

       else if ((bfact_f.eq.1).or. &
                (local_enable_spectral.eq.0).or. &
                (testmask.eq.0).or. &
                (spectral_override.eq.LOW_ORDER_AVGDOWN).or. &
                (spectral_override.eq.ICEMASK_AVGDOWN)) then

        call fine_subelement_stencilMAC(ic,jc,kc,stenlo,stenhi, &
         bfact_c,bfact_f,grid_type)

        do ifine=stenlo(1),stenhi(1)
         if (box_type(1).eq.1) then
          call intersect_weightMAC_avg(ic,ifine,bfact_c,bfact_f,wt(1))
         else if (box_type(1).eq.0) then
          call intersect_weight_avg(ic,ifine,bfact_c,bfact_f,wt(1))
         else 
          print *,"box_type(1) invalid"
          stop
         endif
         if (wt(1).gt.zero) then
          do jfine=stenlo(2),stenhi(2)
           if (box_type(2).eq.1) then
            call intersect_weightMAC_avg(jc,jfine,bfact_c,bfact_f,wt(2))
           else if (box_type(2).eq.0) then
            call intersect_weight_avg(jc,jfine,bfact_c,bfact_f,wt(2))
           else 
            print *,"box_type(2) invalid"
            stop
           endif
           if (wt(2).gt.zero) then
            do kfine=stenlo(3),stenhi(3)
             if (SDIM.eq.3) then
              if (box_type(SDIM).eq.1) then
               call intersect_weightMAC_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
              else if (box_type(SDIM).eq.0) then
               call intersect_weight_avg(kc,kfine,bfact_c,bfact_f,wt(SDIM))
              else 
               print *,"box_type(SDIM) invalid"
               stop
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

               if (volall.gt.zero) then

                if ((spectral_override.eq.LOW_ORDER_AVGDOWN).or. &
                    (spectral_override.eq.SPECTRAL_ORDER_AVGDOWN)) then
               
                 do n=1,ncomp
                  crse_value(n) = crse_value(n)+ &
                    volall*fine(D_DECL(ifine,jfine,kfine),n)
                 enddo
                 voltotal=voltotal+volall

                else if (spectral_override.eq.ICEMASK_AVGDOWN) then

                 voltotal=one

                 do n=1,ncomp
                  icemask=fine(D_DECL(ifine,jfine,kfine),n)
                  if (icemask.eq.zero) then
                   crse_value(n)=zero
                  else if (icemask.eq.one) then
                   if (crse_value(n).eq.one) then
                    ! do nothing
                   else if (crse_value(n).eq.zero) then
                    ! do nothing
                   else
                    print *,"crse_value(n) invalid: ",crse_value(n)
                    stop
                   endif
                  else
                   print *,"icemask invalid: ",icemask
                   stop
                  endif
                 enddo  !n=1,ncomp
              
                else
                 print *,"spectral_override invalid: ",spectral_override
                 stop
                endif

               else if (volall.eq.zero) then
                ! do nothing
               else
                print *,"volall invalid: ",volall
                stop
               endif

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

       if (voltotal.gt.zero) then
        ! do nothing
       else
        print *,"voltotal invalid edgeavgdown: ",voltotal
        stop
       endif

       do n=1,ncomp
        crse(D_DECL(ic,jc,kc),n)=crse_value(n)/voltotal
        if ((levelrz.eq.COORDSYS_CARTESIAN).or. &
            (levelrz.eq.COORDSYS_CYLINDRICAL)) then
         ! do nothing
        else if (levelrz.eq.COORDSYS_RZ) then
         if (SDIM.ne.2) then
          print *,"dimension bust"
          stop
         endif
         if ((xsten(0,1).le.EPS2*dx(1)).and. &
             (box_type(1).eq.1)) then
          crse(D_DECL(ic,jc,kc),n) = zero
         endif
        else
         print *,"levelrz invalid edge avgdown: ",levelrz
         stop
        endif
       enddo ! n=1,ncomp

      enddo 
      enddo
      enddo  ! ic,jc,kc

      deallocate(ffine)

      return
      end subroutine fort_edgeavgdown


      subroutine fort_metrics( &
       xlo,dx, &
       areax,DIMS(areax), &
       areay,DIMS(areay), &
       areaz,DIMS(areaz), &
       vol,DIMS(vol), &
       tilelo,tilehi, &
       fablo,fabhi, &
       bfact, &
       level, &
       finest_level, &
       ngrow, &
       rzflag) &
      bind(c,name='fort_metrics')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      integer, INTENT(in) :: level
      integer, INTENT(in) :: finest_level
      integer, INTENT(in) :: tilelo(SDIM),tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM),fabhi(SDIM)
      integer :: growlo(3),growhi(3)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: ngrow,rzflag
      integer, INTENT(in) :: DIMDEC(vol)
      integer, INTENT(in) :: DIMDEC(areax)
      integer, INTENT(in) :: DIMDEC(areay)
      integer, INTENT(in) :: DIMDEC(areaz)
      
      real(amrex_real), INTENT(out),target :: vol(DIMV(vol))
      real(amrex_real), pointer :: vol_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(out),target :: areax(DIMV(areax))
      real(amrex_real), INTENT(out),target :: areay(DIMV(areay))
      real(amrex_real), INTENT(out),target :: areaz(DIMV(areaz))
      real(amrex_real), pointer :: areax_ptr(D_DECL(:,:,:))
      real(amrex_real), pointer :: areay_ptr(D_DECL(:,:,:))
      real(amrex_real), pointer :: areaz_ptr(D_DECL(:,:,:))
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)

      integer i,j,k,dir
      integer, parameter :: nhalf=3
      real(amrex_real) xsten(-nhalf:nhalf,SDIM)
      integer side
      real(amrex_real) local_area
      integer at_z_axis

      if (bfact.lt.1) then
       print *,"bfact invalid167" 
       stop
      endif
      if (levelrz.ne.rzflag) then 
       print *,"levelrz,rzflag mismatch"
       stop
      endif
      if (SDIM.eq.3) then
       if (levelrz.eq.COORDSYS_RZ) then
        print *,"cannot have levelrz=1"
        stop
       endif
      endif

      if (ngrow.ge.0) then
       ! do nothing
      else
       print *,"ngrow invalid in fort_metrics ",ngrow
       stop
      endif

      if ((level.ge.0).and.(level.le.finest_level)) then
       ! do nothing
      else
       print *,"level or finest_level invalid ",level,finest_level
       stop
      endif

      vol_ptr=>vol
      areax_ptr=>areax
      areay_ptr=>areay
      areaz_ptr=>areaz
      call checkbound_array1(fablo,fabhi,vol_ptr,ngrow,-1)
      call checkbound_array1(fablo,fabhi,areax_ptr,ngrow,0)
      call checkbound_array1(fablo,fabhi,areay_ptr,ngrow,1)
      call checkbound_array1(fablo,fabhi,areaz_ptr,ngrow,SDIM-1)

      call growntilebox(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow)
      do k=growlo(3),growhi(3)
      do j=growlo(2),growhi(2)
      do i=growlo(1),growhi(1)
       call gridsten_level(xsten,i,j,k,level,nhalf)
       call gridvol(xsten,nhalf,rzflag,vol(D_DECL(i,j,k)))
      enddo
      enddo
      enddo
     
      do dir=0,SDIM-1
       call growntileboxMAC(tilelo,tilehi,fablo,fabhi,growlo,growhi,ngrow,dir)
       do k=growlo(3),growhi(3)
       do j=growlo(2),growhi(2)
       do i=growlo(1),growhi(1)
        call gridsten_level(xsten,i,j,k,level,nhalf)
        side=0
        call gridarea(xsten,nhalf,rzflag,dir,side,local_area)
        at_z_axis=0
        if (rzflag.eq.COORDSYS_CARTESIAN) then
         ! do nothing
        else if (rzflag.eq.COORDSYS_RZ) then
         if (dir.eq.1) then
          ! do nothing
         else if (dir.eq.0) then
          if (abs(xsten(-1,1)).le.EPS2*dx(dir+1)) then
           at_z_axis=1
          endif
         else
          print *,"dir invalid"
          stop
         endif
        else if (rzflag.eq.COORDSYS_CYLINDRICAL) then
         if ((dir.eq.1).or.(dir.eq.SDIM-1)) then
          ! do nothing
         else if (dir.eq.0) then
          if (abs(xsten(-1,1)).le.EPS2*dx(dir+1)) then
           at_z_axis=1
          endif
         else
          print *,"dir invalid"
          stop
         endif
        else
         print *,"rzflag invalid"
         stop
        endif
        if ((local_area.gt.zero).or. &
            ((local_area.eq.zero).and.(at_z_axis.eq.1))) then
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
        else
         print *,"local_area invalid"
         stop
        endif
       enddo
       enddo
       enddo ! i,j,k
      enddo ! dir=0..sdim-1 

      return
      end subroutine fort_metrics

      end module navierstokesf90_module


      module OUTPUT_PC_module

       use iso_c_binding
       use amrex_fort_module, only : amrex_real,amrex_particle_real

       implicit none

       type, bind(C) :: particle_t
         real(amrex_particle_real) :: pos(SDIM)
         real(amrex_particle_real) :: extra_state(N_EXTRA_REAL)
         integer(c_int) :: id
         integer(c_int) :: cpu
         integer(c_int) :: extra_int(N_EXTRA_INT)
       end type particle_t

      contains

      subroutine fort_particle_grid( &
        tid, &
        xlo,dx, &
        particles, & ! a list of particles in the elastic structure
        Np, & !  Np = number of particles, pass by value
        tilelo,tilehi, &
        fablo,fabhi, &
        bfact, &
        level, &
        gridno) &
      bind(c,name='fort_particle_grid')

      use probf90_module
      use navierstokesf90_module
      use global_utility_module
      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      integer, INTENT(in) :: tid
      integer, INTENT(in) :: tilelo(SDIM), tilehi(SDIM)
      integer, INTENT(in) :: fablo(SDIM), fabhi(SDIM)
      integer, INTENT(in) :: bfact
      integer, INTENT(in) :: level,gridno
      real(amrex_real), INTENT(in) :: xlo(SDIM),dx(SDIM)
      integer, value, INTENT(in) :: Np ! pass by value
      type(particle_t), INTENT(in) :: particles(Np)

      character*28 cennamestr28
      character*3 levstr
      character*5 gridstr
      character*36 cenfilename36

      real(amrex_real) xref(SDIM)
      real(amrex_real) xrefT(SDIM)
      integer ipart_counter
      integer i,dir
      real(amrex_real) int_to_real_var

      if ((tid.lt.0).or.(tid.ge.geom_nthreads)) then
       print *,"tid invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid151"
       stop
      endif

      write(cennamestr28,'(A14,A14)') &
          './temptecplot/','tempPARCON_pos'

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
      write(cenfilename36,'(A28,A3,A5)') cennamestr28,levstr,gridstr
      print *,"cenfilename36 ",cenfilename36

      if (N_EXTRA_REAL.eq.1) then
       ! do nothing
      else
       print *,"N_EXTRA_REAL unexpected value"
       stop
      endif
      if (N_EXTRA_INT.eq.1) then
       ! do nothing
      else
       print *,"N_EXTRA_INT unexpected value"
       stop
      endif

      open(unit=12,file=cenfilename36)
      write(12,*) Np

      do ipart_counter=1,Np
       do dir=1,SDIM
        xref(dir)=particles(ipart_counter)%pos(dir)
        xrefT(dir)=xref(dir)
       enddo
       if (visual_RT_transform.eq.1) then
        call RT_transform(xref,xrefT)
       endif
       do dir=1,SDIM
        write(12,'(E25.16)',ADVANCE="NO") xrefT(dir)
       enddo

       do dir=1,N_EXTRA_REAL
        write(12,'(E25.16)',ADVANCE="NO") &
          particles(ipart_counter)%extra_state(dir)
       enddo ! dir=1..N_EXTRA_REAL

       do dir=1,N_EXTRA_INT
        int_to_real_var=particles(ipart_counter)%extra_int(dir)
        if (dir.lt.N_EXTRA_INT) then
         write(12,'(E25.16)',ADVANCE="NO") int_to_real_var
        else if (dir.eq.N_EXTRA_INT) then
         write(12,'(E25.16)') int_to_real_var
        else
         print *,"dir invalid (dir=1,N_EXTRA_INT): ",dir
         stop
        endif
       enddo ! dir=1..N_EXTRA_INT

      enddo ! ipart_counter=1,Np

      close(12)

      return
      end subroutine fort_particle_grid

      subroutine fort_combine_particles( &
       grids_per_level,finest_level,nsteps, &
       arrdim,time,plotint) &
      bind(c,name='fort_combine_particles')

      use global_utility_module
      use probcommon_module

      IMPLICIT NONE

      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: plotint
      integer :: strandid
      integer, INTENT(in) :: arrdim,finest_level,nsteps
      integer, INTENT(in) :: grids_per_level(arrdim)

      character*28 cennamestr28
      character*10 newcennamestr10

      character*3 levstr
      character*5 gridstr

      character*36 cenfilename36

      character*6 stepstr

      character*20 newcenfilename20

      integer i
      integer ilev,igrid,ipass
      real(amrex_real) xref(SDIM+N_EXTRA_REAL+N_EXTRA_INT)
      integer nparticles,Part_nparticles
      integer alloc_flag
      integer istruct

      alloc_flag=0

      write(cennamestr28,'(A14,A14)') &
          './temptecplot/','tempPARCON_pos'
      
      write(newcennamestr10,'(A10)') 'PARCON_pos'

      nparticles=0

      if (arrdim.ne.finest_level+1) then
       print *,"arrdim invalid"
       stop
      endif

      do ipass=0,1

       if (ipass.eq.1) then

        alloc_flag=alloc_flag+1

        write(stepstr,'(I6)') nsteps
        do i=1,6
         if (stepstr(i:i).eq.' ') then
          stepstr(i:i)='0'
         endif
        enddo

        if (plotint.le.0) then
         strandid=1
        else
         strandid=(nsteps/plotint)+1
        endif

        write(newcenfilename20,'(A10,A6,A4)') newcennamestr10,stepstr,'.tec'
        print *,"newcenfilename20 ",newcenfilename20
        open(unit=12,file=newcenfilename20)

        if (N_EXTRA_REAL.eq.1) then
         ! do nothing
        else
         print *,"N_EXTRA_REAL invalid"
         stop
        endif
        if (N_EXTRA_INT.eq.1) then
         ! do nothing
        else
         print *,"N_EXTRA_INT invalid"
         stop
        endif

        if (SDIM.eq.3) then
         write(12,*) 'TITLE = "3D particles" '
         write(12,'(A42)',ADVANCE="NO") &
          'VARIABLES = "X","Y","Z","LS","material id"'
        else if (SDIM.eq.2) then
         write(12,*) 'TITLE = "2D particles" '
         write(12,'(A38)',ADVANCE="NO") &
          'VARIABLES = "X","Y","LS","material id"'
        else
         print *,"dimension bust"
         stop
        endif

        write(12,*) ','

        if (plotint.le.0) then
         strandid=1
        else
         strandid=(nsteps/plotint)+1
        endif

        write(12,'(A19,I14,A26,E25.16,A10,I10)') & 
          'ZONE F="POINT", I= ', nparticles,  &
          ', J=1, K=1, SOLUTIONTIME= ',round_time(time),' STRANDID=',strandid

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

         write(cenfilename36,'(A28,A3,A5)') cennamestr28,levstr,gridstr
         print *,"cenfilename36 ",cenfilename36
         open(unit=5,file=cenfilename36)

         read(5,*) Part_nparticles

         if (ipass.eq.0) then
          nparticles=nparticles+Part_nparticles
         else if (ipass.eq.1) then

          do i=1,Part_nparticles
           read(5,*) &
             (xref(istruct),istruct=1, &
              SDIM+N_EXTRA_REAL+N_EXTRA_INT)
           write(12,*) &
             (xref(istruct),istruct=1, &
              SDIM+N_EXTRA_REAL+N_EXTRA_INT)
          enddo

         else
          print *,"ipass invalid"
          stop
         endif
         close(5)
       enddo ! igrid
       enddo ! ilev

       if (ipass.eq.1) then
        close(12)
       endif

      enddo ! ipass=0..1

      alloc_flag=alloc_flag-1

      if (alloc_flag.gt.0) then
       print *,"alloc_flag bust"
       stop
      endif

      return
      end subroutine fort_combine_particles

      end module OUTPUT_PC_module


